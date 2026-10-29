// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
// SPDX-License-Identifier: MPL-2.0

#include "alea_test.h"
#include "alea.h"
#include "alea_compute.h"
#include "core/alea_system.h"
#include "raycast/raycast.h"
#include <float.h>
#include <math.h>

static alea_system_t* sphere_scene(void) {
    alea_system_t* sys = alea_create();
    int sphere = alea_sphere_surface(sys, 1, 0, 0, 0, 2);
    int mat = alea_add_material(sys, 7);
    alea_add_cell(sys, 42, alea_halfspace(sys, sphere, -1), mat, -1, 0);
    return sys;
}

static int count_progress(void* user,const render_framebuffer_t* fb,size_t completed,size_t total) {
    (void)fb;(void)completed;(void)total;++*(int*)user;return 0;
}
TEST(progressive_requires_metal_ray_context) {
    alea_system_t* sys=sphere_scene();
    alea_compute_t* ctx=alea_compute_create(sys,ALEA_COMPUTE_CPU,ALEA_COMPUTE_FP64);ASSERT_NOT_NULL(ctx);
    render_config_t cfg;render_config_init(&cfg);cfg.width=cfg.height=2;
    render_camera_t cam;ASSERT_EQ(render_camera_setup(&cam,&cfg,sys),0);
    render_framebuffer_t* fb=render_framebuffer_create(2,2,1);ASSERT_NOT_NULL(fb);int calls=0;
    ASSERT_EQ(alea_compute_render_progressive(ctx,&cfg,&cam,fb,count_progress,&calls),-1);
    ASSERT_EQ(alea_error_code(),ALEA_ERR_UNSUPPORTED);ASSERT_EQ(calls,0);
    ASSERT_EQ(alea_compute_render_progressive(NULL,&cfg,&cam,fb,count_progress,&calls),-1);
    ASSERT_EQ(alea_error_code(),ALEA_ERR_NULL_ARG);
    render_framebuffer_free(fb);render_config_free(&cfg);alea_compute_destroy(ctx);alea_destroy(sys);
}

TEST(cpu_snapshot_and_validation) {
    alea_system_t* sys = sphere_scene();
    alea_compute_t* ctx = alea_compute_create(sys, ALEA_COMPUTE_CPU, ALEA_COMPUTE_FP64);
    ASSERT_NOT_NULL(ctx);
    alea_destroy(sys);
    double xyz[] = {0,0,0, 2,0,0, 3,0,0};
    alea_compute_result_t r[3];
    ASSERT_EQ(alea_compute_classify(ctx, xyz, 3, r), 0);
    ASSERT_EQ(r[0].cell_id, 42); ASSERT_EQ(r[0].material_id, 7);
    ASSERT_EQ(r[1].owner_count, 1); /* Surface belongs to the closed halfspace. */
    ASSERT_EQ(r[2].owner_count, 0); ASSERT_EQ(r[2].cell_id, -1);
    ASSERT_EQ(alea_compute_classify(ctx, NULL, 0, NULL), 0);
    ASSERT_EQ(alea_compute_classify(ctx, NULL, 1, r), -1);
    ASSERT_EQ(alea_compute_classify(ctx, xyz, SIZE_MAX, r), -1);
    xyz[0] = NAN;
    ASSERT_EQ(alea_compute_classify(ctx, xyz, 1, r), -1);
    xyz[0] = INFINITY;
    ASSERT_EQ(alea_compute_classify(ctx, xyz, 1, r), -1);
    xyz[0] = 0;
    alea_interrupt();
    ASSERT_EQ(alea_compute_classify(ctx, xyz, 1, r), -1);
    alea_clear_interrupt();
    alea_compute_destroy(ctx);
    alea_compute_destroy(NULL);
}

TEST(backend_precision_and_missing_root_rejections) {
    alea_system_t* sys = sphere_scene();
    ASSERT_NULL(alea_compute_create(sys, ALEA_COMPUTE_METAL, ALEA_COMPUTE_FP64));
    ASSERT_EQ(alea_error_code(), ALEA_ERR_UNSUPPORTED);
    ASSERT_NULL(alea_compute_create(sys, ALEA_COMPUTE_CPU, ALEA_COMPUTE_FP32));
    ASSERT_NULL(alea_compute_create(sys, (alea_compute_backend_t)99, ALEA_COMPUTE_FP64));
    ASSERT_NULL(alea_compute_create(NULL, ALEA_COMPUTE_CPU, ALEA_COMPUTE_FP64));
#ifndef ALEA_USE_METAL
    ASSERT_FALSE(alea_compute_available(ALEA_COMPUTE_METAL));
    ASSERT_NULL(alea_compute_create(sys, ALEA_COMPUTE_METAL, ALEA_COMPUTE_FP32));
#endif
    sys->cells.data[0].universe_id = 1;
    alea_compute_t* ctx=alea_compute_create(sys,ALEA_COMPUTE_CPU,ALEA_COMPUTE_FP64);
    ASSERT_NOT_NULL(ctx);
    double xyz[3]={0};alea_compute_result_t result;
    ASSERT_EQ(alea_compute_classify(ctx,xyz,1,&result),-1);
    ASSERT_EQ(alea_error_code(),ALEA_ERR_INVALID_STATE);
    alea_compute_destroy(ctx);
    alea_destroy(sys);
}

static void check_scene(alea_system_t* sys, const double* xyz, size_t count,
                         alea_compute_backend_t backend) {
    alea_compute_t* ctx = alea_compute_create(sys, backend,
        backend == ALEA_COMPUTE_CPU ? ALEA_COMPUTE_FP64 : ALEA_COMPUTE_FP32);
    ASSERT_MSG(ctx != NULL, alea_error());
    alea_compute_result_t* r = malloc(count * sizeof(*r));
    ASSERT_NOT_NULL(r);
    ASSERT_MSG(alea_compute_classify(ctx, xyz, 1, r) == 0, alea_error());
    ASSERT_MSG(alea_compute_classify(ctx, xyz, count, r) == 0, alea_error());
    for (size_t i = 0; i < count; ++i) {
        unsigned owners = 0;
        int cell = -1, material = -1;
        for (size_t j = 0; j < sys->cells.count; ++j) {
            const alea_cell_entry_t* c = &sys->cells.data[j];
            if (alea_point_inside(sys, c->root_node_id, xyz[3*i], xyz[3*i+1], xyz[3*i+2])) {
                ++owners; cell = c->mc_cell_id; material = c->material_id;
            }
        }
        ASSERT_EQ(r[i].owner_count, owners);
        ASSERT_EQ(r[i].cell_id, owners == 1 ? cell : -1);
        ASSERT_EQ(r[i].material_id, owners == 1 ? material : -1);
    }
    /* Reuse buffers with smaller and then original-sized batches. */
    ASSERT_EQ(alea_compute_classify(ctx, xyz, 1, r), 0);
    ASSERT_EQ(alea_compute_classify(ctx, xyz, count, r), 0);
    free(r);
    alea_compute_destroy(ctx);
}

static void check_primitives_and_csg(alea_compute_backend_t backend) {
    alea_system_t* sys = alea_create();
    int surfaces[] = {
        alea_plane_surface(sys, 1, -1, 0, 0, 0), /* canonicalization flips sense */
        alea_sphere_surface(sys, 2, 0, 0, 0, 2),
        alea_cylinder_x_surface(sys, 3, .5, -.5, 1),
        alea_cylinder_y_surface(sys, 4, -.5, .5, 1),
        alea_cylinder_z_surface(sys, 5, .5, .5, 1),
        alea_box_surface(sys, 6, -1, 1, -2, 2, -3, 3)
    };
    alea_node_id_t nodes[6];
    int mat = alea_add_material(sys, 19);
    for (int i = 0; i < 6; ++i) {
        nodes[i] = alea_halfspace(sys, surfaces[i], -1);
        alea_add_cell(sys, i + 1, nodes[i], mat, -1, 0);
        alea_add_cell(sys, i + 11, alea_halfspace(sys, surfaces[i], 1), mat, -1, 0);
    }
    alea_add_cell(sys, 21, alea_union(sys, nodes[0], nodes[1]), mat, -1, 0);
    alea_add_cell(sys, 22, alea_intersection(sys, nodes[1], nodes[2]), mat, -1, 0);
    alea_add_cell(sys, 23, alea_difference(sys, nodes[5], nodes[1]), mat, -1, 0);
    alea_add_cell(sys, 24, alea_complement(sys, nodes[3]), mat, -1, 0);
    /* Binary-exact coordinates exercise boundaries and both sides without
     * conflating the algorithm check with general FP32 rounding differences. */
    size_t count = 17 * 17 * 17;
    double* xyz = malloc(count * 3 * sizeof(double));
    ASSERT_NOT_NULL(xyz);
    size_t at = 0;
    for (int x = -8; x <= 8; ++x)
        for (int y = -8; y <= 8; ++y)
            for (int z = -8; z <= 8; ++z) {
                xyz[at++] = x * .5; xyz[at++] = y * .5; xyz[at++] = z * .5;
            }
    check_scene(sys, xyz, count, backend);
    free(xyz);
    alea_destroy(sys);
}

TEST(cpu_matches_existing_evaluator) { check_primitives_and_csg(ALEA_COMPUTE_CPU); }

TEST(empty_scene_and_overlaps) {
    alea_system_t* sys = alea_create();
    double xyz[] = {0,0,0, 1,0,0, 3,0,0};
    check_scene(sys, xyz, 3, ALEA_COMPUTE_CPU);
    int s = alea_sphere_surface(sys, 1, 0, 0, 0, 2);
    int m = alea_add_material(sys, 7);
    alea_add_cell(sys, 1, alea_halfspace(sys, s, -1), m, -1, 0);
    alea_add_cell(sys, 2, alea_halfspace(sys, s, -1), m, -1, 0);
    check_scene(sys, xyz, 3, ALEA_COMPUTE_CPU);
    alea_destroy(sys);
}

TEST(cpu_render_owned_snapshot) {
    alea_system_t* sys=sphere_scene();
    alea_compute_t* cpu=alea_compute_create(sys,ALEA_COMPUTE_CPU,ALEA_COMPUTE_FP64);
    ASSERT_NOT_NULL(cpu);
    render_config_t cfg; render_config_init(&cfg);
    cfg.width=cfg.height=1; cfg.eye_set=cfg.target_set=1; cfg.log_level=0;
    cfg.eye[0]=0; cfg.eye[1]=-10; cfg.eye[2]=0;
    cfg.target[0]=cfg.target[1]=cfg.target[2]=0;
    render_camera_t cam; ASSERT_EQ(render_camera_setup(&cam,&cfg,sys),0);
    alea_destroy(sys);
    render_framebuffer_t* fb=render_framebuffer_create(1,1,1); ASSERT_NOT_NULL(fb);
    ASSERT_MSG(alea_compute_render(cpu,&cfg,&cam,fb)==0,alea_error());
    ASSERT_EQ(fb->cell_id[0],42); ASSERT_EQ(fb->material_id[0],7);
    ASSERT_NEAR(fb->depth[0],8,1e-6);
    cfg.width=2; ASSERT_EQ(alea_compute_render(cpu,&cfg,&cam,fb),-1);
    ASSERT_EQ(alea_compute_render(NULL,&cfg,&cam,fb),-1);
    render_framebuffer_free(fb); alea_compute_destroy(cpu); render_config_free(&cfg);
}

TEST(cpu_render_nested_owned_snapshot) {
    alea_system_t* sys=sphere_scene();sys->cells.data[0].universe_id=1;
    alea_add_cell(sys,43,alea_complement(sys,sys->cells.data[0].root_node_id),ALEA_MATERIAL_VOID,0,1);
    int outer=alea_sphere_surface(sys,2,0,0,0,4);
    int container=alea_add_cell(sys,10,alea_halfspace(sys,outer,-1),ALEA_MATERIAL_VOID,0,0);
    ASSERT_EQ(alea_set_cell_fill(sys,container,1,0),0);
    alea_compute_t* ctx=alea_compute_create(sys,ALEA_COMPUTE_CPU,ALEA_COMPUTE_FP64);ASSERT_NOT_NULL(ctx);
    render_config_t cfg;render_config_init(&cfg);cfg.width=cfg.height=1;cfg.log_level=0;
    cfg.eye_set=cfg.target_set=1;cfg.eye[0]=0;cfg.eye[1]=-10;cfg.eye[2]=0;
    cfg.target[0]=cfg.target[1]=cfg.target[2]=0;
    render_camera_t cam;ASSERT_EQ(render_camera_setup(&cam,&cfg,sys),0);alea_destroy(sys);
    render_framebuffer_t* fb=render_framebuffer_create(1,1,1);ASSERT_NOT_NULL(fb);
    ASSERT_MSG(alea_compute_render(ctx,&cfg,&cam,fb)==0,alea_error());
    ASSERT_EQ(fb->cell_id[0],42);ASSERT_NEAR(fb->depth[0],8,1e-6);
    render_framebuffer_free(fb);render_config_free(&cfg);alea_compute_destroy(ctx);
}

/* Analytic nested translation/rotation: the material box is centered at
 * world (1,1,0), the terminal void at (1,2.5,0). An ancestor clips z < 0. */
static alea_system_t* point_hierarchy_scene(int overlap) {
    alea_system_t* sys=alea_create();if(!sys)abort();
    int material=alea_add_material(sys,7);
    int outer=alea_box_surface(sys,1,-4,4,-4,4,-4,4);
    int parent=alea_add_cell(sys,10,alea_halfspace(sys,outer,-1),ALEA_MATERIAL_VOID,0,0);
    double translation[]={1,0,0};
    if(alea_add_transform(sys,1,translation,3,0) || alea_set_fill(sys,parent,1,1))abort();
    if(overlap) {
        int duplicate=alea_add_cell(sys,11,alea_halfspace(sys,outer,-1),ALEA_MATERIAL_VOID,0,0);
        if(alea_set_fill(sys,duplicate,1,1))abort();
    }
    int inner=alea_box_surface(sys,2,-3,3,-3,3,0,3);
    int child=alea_add_cell(sys,20,alea_halfspace(sys,inner,-1),ALEA_MATERIAL_VOID,0,1);
    double rotation[]={0,1,0,0,1,0,-1,0,0,0,0,1};
    if(alea_add_transform(sys,2,rotation,12,0) || alea_set_fill(sys,child,2,2))abort();
    int solid=alea_box_surface(sys,3,-.5,.5,-.5,.5,-.5,.5);
    int empty=alea_box_surface(sys,4,1,2,-.5,.5,-.5,.5);
    if(alea_add_cell(sys,42,alea_halfspace(sys,solid,-1),material,-1,2)<0 ||
       alea_add_cell(sys,43,alea_halfspace(sys,empty,-1),ALEA_MATERIAL_VOID,0,2)<0 ||
       alea_add_cell(sys,99,alea_halfspace(sys,outer,-1),material,-1,99)<0)abort();
    return sys;
}

static void check_point_hierarchy(alea_compute_backend_t backend) {
    const double points[]={1,1,.25, 1,2.5,.25, 1,1,-.25, 0,0,0, 100,0,0};
    for(int overlap=0;overlap<2;++overlap) {
        alea_system_t* sys=point_hierarchy_scene(overlap);
        alea_compute_t* ctx=alea_compute_create(sys,backend,
            backend==ALEA_COMPUTE_CPU?ALEA_COMPUTE_FP64:ALEA_COMPUTE_FP32);
        ASSERT_MSG(ctx!=NULL,alea_error());
        alea_destroy(sys);
        /* Cross the native submission boundary and reuse buffers at two sizes. */
        const size_t n=4101;
        double* xyz=malloc(n*3*sizeof(*xyz));
        alea_compute_result_t* results=malloc(n*sizeof(*results));
        ASSERT_NOT_NULL(xyz);ASSERT_NOT_NULL(results);
        for(size_t i=0;i<n;++i)memcpy(xyz+3*i,points+3*(i%5),3*sizeof(double));
        for(int pass=0;pass<2;++pass) {
            size_t count=pass?n:5;
            ASSERT_MSG(alea_compute_classify(ctx,xyz,count,results)==0,alea_error());
            for(size_t i=0;i<count;++i) {
                unsigned owners=i%5<2?(unsigned)(overlap+1):0;
                ASSERT_EQ(results[i].owner_count,owners);
                ASSERT_EQ(results[i].cell_id,owners==1?(i%5==0?42:43):-1);
                ASSERT_EQ(results[i].material_id,owners==1?(i%5==0?7:0):-1);
            }
        }
        alea_interrupt();int rc=alea_compute_classify(ctx,xyz,n,results);alea_clear_interrupt();
        ASSERT_EQ(rc,-1);
        ASSERT_EQ(alea_compute_classify(ctx,xyz,n,results),0);
        free(results);free(xyz);alea_compute_destroy(ctx);
    }
    /* A root terminal plus an unused universe is not a flat scene, but the
     * unused definition must never be counted as a second owner. */
    alea_system_t* sys=sphere_scene();
    ASSERT(alea_add_cell(sys,99,sys->cells.data[0].root_node_id,0,0,99)>=0);
    alea_compute_t* ctx=alea_compute_create(sys,backend,
        backend==ALEA_COMPUTE_CPU?ALEA_COMPUTE_FP64:ALEA_COMPUTE_FP32);
    ASSERT_MSG(ctx!=NULL,alea_error());double origin[]={0,0,0};alea_compute_result_t result;
    ASSERT_EQ(alea_compute_classify(ctx,origin,1,&result),0);
    ASSERT_EQ(result.owner_count,1);ASSERT_EQ(result.cell_id,42);
    alea_compute_destroy(ctx);alea_destroy(sys);
}

static void check_point_depth_limit(alea_compute_backend_t backend) {
    for(int levels=32;levels<=33;++levels) {
        alea_system_t* sys=alea_create();int mat=alea_add_material(sys,7);
        int sphere=alea_sphere_surface(sys,1,0,0,0,2);
        for(int level=0;level<levels;++level) {
            int cell=alea_add_cell(sys,level+1,alea_halfspace(sys,sphere,-1),
                level==levels-1?mat:ALEA_MATERIAL_VOID,0,level);
            ASSERT(cell>=0);
            if(level<levels-1)ASSERT_EQ(alea_set_fill(sys,cell,level+1,0),0);
        }
        alea_compute_t* ctx=alea_compute_create(sys,backend,
            backend==ALEA_COMPUTE_CPU?ALEA_COMPUTE_FP64:ALEA_COMPUTE_FP32);
        if(levels==33 && backend==ALEA_COMPUTE_METAL) {
            ASSERT_NULL(ctx);ASSERT_EQ(alea_error_code(),ALEA_ERR_UNSUPPORTED);
        } else {
            ASSERT_MSG(ctx!=NULL,alea_error());
            double point[]={0,0,0};alea_compute_result_t result;
            int rc=alea_compute_classify(ctx,point,1,&result);
            if(levels==32) {ASSERT_EQ(rc,0);ASSERT_EQ(result.owner_count,1);ASSERT_EQ(result.cell_id,32);}
            else ASSERT_EQ(rc,-1);
            alea_compute_destroy(ctx);
        }
        alea_destroy(sys);
    }
}

TEST(cpu_point_depth_limit) {check_point_depth_limit(ALEA_COMPUTE_CPU);}

TEST(cpu_classify_transformed_terminal_occurrences) {check_point_hierarchy(ALEA_COMPUTE_CPU);}

TEST(cpu_classify_reports_hierarchy_limits) {
    alea_system_t* sys=sphere_scene();
    ASSERT_EQ(alea_set_fill(sys,0,1,0),0);
    ASSERT(alea_add_cell(sys,43,sys->cells.data[0].root_node_id,ALEA_MATERIAL_VOID,0,1)>=0);
    ASSERT_EQ(alea_set_fill(sys,1,1,0),0);
    alea_compute_t* ctx=alea_compute_create(sys,ALEA_COMPUTE_CPU,ALEA_COMPUTE_FP64);
    ASSERT_NOT_NULL(ctx);double xyz[]={0,0,0};alea_compute_result_t result;
    ASSERT_EQ(alea_compute_classify(ctx,xyz,1,&result),-1);
    alea_compute_destroy(ctx);alea_destroy(sys);
}

#ifdef ALEA_USE_METAL
#include "compute/compute_internal.h"
#include <unistd.h>
/* A requested Metal test run must actually exercise the GPU; do not skip it
 * when the device is unavailable (including restricted execution sandboxes). */
TEST(metal_point_depth_limit) {check_point_depth_limit(ALEA_COMPUTE_METAL);}

TEST(metal_classify_transformed_terminal_occurrences) {check_point_hierarchy(ALEA_COMPUTE_METAL);}

TEST(metal_matches_existing_evaluator) {
    ASSERT_TRUE(alea_compute_available(ALEA_COMPUTE_METAL));
    check_primitives_and_csg(ALEA_COMPUTE_METAL);
}

TEST(metal_empty_snapshot_and_precision) {
    ASSERT_TRUE(alea_compute_available(ALEA_COMPUTE_METAL));
    alea_system_t* sys = alea_create();
    double xyz[] = {0,0,0, 1,0,0, 3,0,0};
    check_scene(sys, xyz, 3, ALEA_COMPUTE_METAL);
    alea_destroy(sys);
    sys = sphere_scene();
    alea_compute_t* gpu = alea_compute_create(sys, ALEA_COMPUTE_METAL, ALEA_COMPUTE_FP32);
    alea_compute_t* cpu = alea_compute_create(sys, ALEA_COMPUTE_CPU, ALEA_COMPUTE_FP64);
    ASSERT_MSG(gpu != NULL, alea_error());
    ASSERT_NOT_NULL(cpu);
    alea_destroy(sys);
    /* The public precision contract deliberately exposes this difference. */
    xyz[0] = 2.0 + 1e-8;
    alea_compute_result_t a, b;
    ASSERT_EQ(alea_compute_classify(cpu, xyz, 1, &a), 0);
    ASSERT_EQ(alea_compute_classify(gpu, xyz, 1, &b), 0);
    ASSERT_EQ(a.owner_count, 0);
    ASSERT_EQ(b.owner_count, 1);
    xyz[0] = DBL_MAX;
    ASSERT_EQ(alea_compute_classify(gpu, xyz, 1, &b), -1);
    ASSERT_EQ(alea_compute_classify(gpu, NULL, 0, NULL), 0);
    alea_compute_destroy(gpu); alea_compute_destroy(cpu);
}

TEST(metal_rejects_cycles_and_accepts_deep_geometry) {
    ASSERT_TRUE(alea_compute_available(ALEA_COMPUTE_METAL));
    alea_system_t* sys = sphere_scene();
    /* Use internal representation to exercise malformed and depth-limited DAGs
     * without the construction API simplifying nested complements away. */
    alea_node_id_t root = sys->cells.data[0].root_node_id;
    alea_node_t saved = sys->nodes.data[root];
    ALEA_SET_OPERATION(&sys->nodes.data[root], ALEA_OP_COMPLEMENT);
    sys->nodes.data[root].operation.left = root;
    ASSERT_NULL(alea_compute_create(sys, ALEA_COMPUTE_METAL, ALEA_COMPUTE_FP32));
    sys->nodes.data[root] = saved;
    sys->cells.data[0].root_node_id = root;
    for (int i = 0; i < 550; ++i) {
        alea_node_t n = saved;
        ALEA_SET_OPERATION(&n, ALEA_OP_COMPLEMENT);
        n.operation.left = root;
        n.bbox=(alea_node_bbox_t){-INFINITY,INFINITY,-INFINITY,INFINITY,-INFINITY,INFINITY};
        ASSERT_EQ(alea_vec_push(&sys->nodes, n, alea_node_t), 0);
        root = (uint32_t)(sys->nodes.count - 1);
    }
    sys->cells.data[0].root_node_id = root;
    double xyz[] = {0,0,0, 2,0,0, 3,0,0};
    check_scene(sys, xyz, 3, ALEA_COMPUTE_METAL); /* Source depth 551 needs only one GPU value slot. */
    alea_node_t n = saved;
    ALEA_SET_OPERATION(&n, ALEA_OP_COMPLEMENT);
    n.operation.left = root;
    n.bbox=(alea_node_bbox_t){-INFINITY,INFINITY,-INFINITY,INFINITY,-INFINITY,INFINITY};
    ASSERT_EQ(alea_vec_push(&sys->nodes, n, alea_node_t), 0);
    sys->cells.data[0].root_node_id = (uint32_t)(sys->nodes.count - 1);
    check_scene(sys, xyz, 3, ALEA_COMPUTE_METAL);
    alea_destroy(sys);
}

TEST(metal_numeric_overflow_is_an_error) {
    alea_system_t* sys = alea_create();
    int sphere = alea_sphere_surface(sys, 1, 0, 0, 0, 1e20);
    int mat = alea_add_material(sys, 7);
    alea_add_cell(sys, 1, alea_halfspace(sys, sphere, -1), mat, -1, 0);
    alea_compute_t* ctx = alea_compute_create(sys, ALEA_COMPUTE_METAL, ALEA_COMPUTE_FP32);
    ASSERT_MSG(ctx != NULL, alea_error());
    /* Finite inputs can still overflow intermediate FP32 arithmetic. */
    double xyz[] = {0, 0, 0};
    alea_compute_result_t r;
    ASSERT_EQ(alea_compute_classify(ctx, xyz, 1, &r), -1);
    ASSERT_EQ(alea_error_code(), ALEA_ERR_OVERFLOW);
    alea_compute_destroy(ctx); alea_destroy(sys);
}
/* Use the full double-precision segment pipeline as the image oracle. The
 * CPU interactive first-visible path can miss quadric/RCC entries and return
 * zero normals; sharing that optimized path would hide GPU correctness. */
static int reference_render(alea_system_t* sys,render_config_t* cfg,
                             render_camera_t* cam,render_framebuffer_t* fb) {
    alea_raycast_result_t* result=alea_raycast_result_create();if(!result)return -1;
    for(int y=0;y<cfg->height;++y)for(int x=0;x<cfg->width;++x) {
        size_t pixel=(size_t)y*cfg->width+x;float sum[3]={0};
        for(int sy=0;sy<cfg->aa_samples;++sy)for(int sx=0;sx<cfg->aa_samples;++sx) {
            double o[3],d[3];render_camera_ray(cam,cfg->width,cfg->height,
                x+(sx+.5)/cfg->aa_samples,y+(sy+.5)/cfg->aa_samples,
                &o[0],&o[1],&o[2],&d[0],&d[1],&d[2]);
            double begin=0,end=1e5;
            if(cfg->num_clips) {
                render_clip_plane_t* c=&cfg->clips[0];double den=0,num=c->d;
                for(int a=0;a<3;++a){den+=c->normal[a]*d[a];num+=c->normal[a]*o[a];}
                if(den>0)begin=fmax(begin,-num/den);else if(den<0)end=fmin(end,-num/den);else if(num<=0)end=-1;
            }
            int cell=-1,material=0;double depth=1e30,normal[3]={0};
            if(begin<end && sys->cells.count) {
                if(alea_raycast(sys,o[0],o[1],o[2],d[0],d[1],d[2],end,result))return -1;
                for(size_t j=0;j<result->segments.count;++j) {
                    alea_ray_segment_t* seg=&result->segments.data[j];
                    if(!seg->material_id || seg->t_exit<=begin || seg->t_enter>=end)continue;
                    cell=seg->cell_id;material=seg->material_id;depth=fmax(begin,seg->t_enter);
                    if(depth==begin && begin>0) {
                        double len=0;for(int a=0;a<3;++a)len+=cfg->clips[0].normal[a]*cfg->clips[0].normal[a];
                        for(int a=0;a<3;++a)normal[a]=-cfg->clips[0].normal[a]/sqrt(len);
                    } else if(seg->enter_hit_index>=0) {
                        const alea_ray_hit_t* hit=&result->hits.data[seg->enter_hit_index];
                        normal[0]=hit->nx;normal[1]=hit->ny;normal[2]=hit->nz;
                    }
                    break;
                }
            }
            float rgb[3];memcpy(rgb,cfg->background,sizeof(rgb));
            if(cell>=0) {
                render_get_color(cfg->color_mode==RENDER_COLOR_CELL?cell:material,cfg->color_mode,cfg,&rgb[0],&rgb[1],&rgb[2]);
                if(depth==begin && begin>0){for(int a=0;a<3;++a)rgb[a]*=cfg->cross_section_tint;}
                else {
                    double n[3];memcpy(n,normal,sizeof(n));double facing=0;
                    for(int a=0;a<3;++a)facing+=n[a]*d[a];if(facing>0)for(int a=0;a<3;++a)n[a]=-n[a];
                    double ndl=0;for(int a=0;a<3;++a)ndl-=n[a]*(cfg->light_set?cfg->light_dir[a]:d[a]);ndl=fmax(ndl,0);
                    double spec=0;for(int a=0;a<3;++a)spec-=(2*ndl*n[a]+(cfg->light_set?cfg->light_dir[a]:d[a]))*d[a];
                    spec=pow(fmax(spec,0),cfg->shininess);
                    for(int a=0;a<3;++a)rgb[a]=fmin(1,fmax(0,rgb[a]*(cfg->ambient+cfg->diffuse*ndl)+cfg->specular*spec));
                }
            }
            for(int a=0;a<3;++a)sum[a]+=rgb[a];
            if(sx==cfg->aa_samples/2 && sy==cfg->aa_samples/2) {
                fb->cell_id[pixel]=cell;fb->material_id[pixel]=material;fb->depth[pixel]=depth;
                double facing=0;for(int a=0;a<3;++a)facing+=normal[a]*d[a];
                for(int a=0;a<3;++a)fb->normal[3*pixel+a]=facing>0?-normal[a]:normal[a];
            }
        }
        for(int a=0;a<3;++a)fb->color[3*pixel+a]=sum[a]/(cfg->aa_samples*cfg->aa_samples);
    }
    alea_raycast_result_destroy(result);return 0;
}

/* Compare actual images, visibility IDs and depths with the existing CPU
 * raycaster. Cameras avoid deliberately ambiguous FP32 boundary probes. */
static int compare_resume=0;
static void compare_render(alea_system_t* sys, alea_compute_t* gpu, render_config_t* cfg) {
    render_camera_t cam;
    ASSERT_EQ(render_camera_setup(&cam,cfg,sys),0);
    render_framebuffer_t* cpu=render_framebuffer_create(cfg->width,cfg->height,1);
    render_framebuffer_t* metal=render_framebuffer_create(cfg->width,cfg->height,1);
    ASSERT_NOT_NULL(cpu); ASSERT_NOT_NULL(metal);
    ASSERT_EQ(reference_render(sys,cfg,&cam,cpu),0);
    ASSERT_MSG(alea_compute_render(gpu,cfg,&cam,metal)==0,alea_error());
    if(compare_resume) {
        render_framebuffer_t* resumed=render_framebuffer_create(cfg->width,cfg->height,1);ASSERT_NOT_NULL(resumed);int calls=0;
        ASSERT_MSG(alea_compute_render_progressive(gpu,cfg,&cam,resumed,count_progress,&calls)==0,alea_error());
        size_t n=(size_t)cfg->width*cfg->height;
        ASSERT_EQ(memcmp(metal->color,resumed->color,n*3*sizeof(float)),0);
        ASSERT_EQ(memcmp(metal->depth,resumed->depth,n*sizeof(float)),0);
        ASSERT_EQ(memcmp(metal->normal,resumed->normal,n*3*sizeof(float)),0);
        ASSERT_EQ(memcmp(metal->cell_id,resumed->cell_id,n*sizeof(int)),0);
        ASSERT_EQ(memcmp(metal->material_id,resumed->material_id,n*sizeof(int)),0);
        ASSERT_TRUE(calls>2);render_framebuffer_free(resumed);
    }
    for (int i=0;i<cfg->width*cfg->height;++i) {
        if (cpu->cell_id[i]!=metal->cell_id[i])
            fprintf(stderr,"pixel %d: CPU cell %d depth %g, Metal cell %d depth %g\n",
                i,cpu->cell_id[i],cpu->depth[i],metal->cell_id[i],metal->depth[i]);
        ASSERT_EQ(cpu->cell_id[i],metal->cell_id[i]);
        ASSERT_EQ(cpu->material_id[i],metal->material_id[i]);
        if (cpu->cell_id[i]>=0) {
            ASSERT_NEAR(cpu->depth[i],metal->depth[i],.002);
            for (int a=0;a<3;++a) ASSERT_NEAR(cpu->normal[3*i+a],metal->normal[3*i+a],.002);
        }
        for (int a=0;a<3;++a) ASSERT_NEAR(cpu->color[3*i+a],metal->color[3*i+a],.002);
    }
    render_framebuffer_free(cpu); render_framebuffer_free(metal);
}

static void test_camera(render_config_t* cfg);

TEST(metal_point_and_render_share_snapshot_in_either_order) {
    for(int render_first=0;render_first<2;++render_first) {
        alea_system_t* sys=point_hierarchy_scene(0);
        alea_compute_t* gpu=alea_compute_create(sys,ALEA_COMPUTE_METAL,ALEA_COMPUTE_FP32);
        ASSERT_MSG(gpu!=NULL,alea_error());alea_destroy(sys);
        alea_system_t* reference=alea_create();int mat=alea_add_material(reference,7);
        int box=alea_box_surface(reference,1,.5,1.5,.5,1.5,0,.5);
        ASSERT(alea_add_cell(reference,42,alea_halfspace(reference,box,-1),mat,-1,0)>=0);
        ASSERT_EQ(alea_prepare_query_acceleration(reference),0);
        render_config_t cfg;test_camera(&cfg);
        if(render_first)compare_render(reference,gpu,&cfg);
        double points[]={1,1,.25,1,2.5,.25};alea_compute_result_t results[2];
        ASSERT_EQ(alea_compute_classify(gpu,points,2,results),0);
        ASSERT_EQ(results[0].cell_id,42);ASSERT_EQ(results[1].cell_id,43);
        compare_render(reference,gpu,&cfg);
        cfg.num_clips=1;cfg.clips[0]=(render_clip_plane_t){{0,1,0},-1};
        compare_render(reference,gpu,&cfg);
        ASSERT_EQ(alea_compute_classify(gpu,points,2,results),0);
        ASSERT_EQ(results[1].owner_count,1);ASSERT_EQ(results[1].material_id,0);
        render_config_free(&cfg);alea_compute_destroy(gpu);alea_destroy(reference);
    }
}

TEST(metal_render_matches_cpu_csg_and_clips) {
    alea_system_t* sys=alea_create();
    int mat=alea_add_material(sys,7);
    int sphere=alea_sphere_surface(sys,1,-3,0,0,2);
    int hole=alea_cylinder_y_surface(sys,2,-3,0,.7);
    alea_add_cell(sys,42,alea_difference(sys,alea_halfspace(sys,sphere,-1),
        alea_halfspace(sys,hole,-1)),mat,-1,0);
    int box=alea_box_surface(sys,3,0,4,-2,2,-2,2);
    int cx=alea_cylinder_x_surface(sys,4,0,0,1.7);
    int cz=alea_cylinder_z_surface(sys,5,2,0,1.2);
    int plane=alea_plane_surface(sys,6,1,.3,.7,-2.4);
    alea_node_id_t root=alea_intersection(sys,alea_halfspace(sys,box,-1),
        alea_union(sys,alea_halfspace(sys,cx,-1),alea_halfspace(sys,cz,-1)));
    root=alea_intersection(sys,root,alea_complement(sys,alea_halfspace(sys,plane,-1)));
    alea_add_cell(sys,43,root,mat,-1,0);
    /* A void region in front must not hide the solids. */
    int voidbox=alea_box_surface(sys,7,-8,8,-8,-4,-8,8);
    alea_add_cell(sys,44,alea_halfspace(sys,voidbox,-1),ALEA_MATERIAL_VOID,0,0);
    ASSERT_EQ(alea_prepare_query_acceleration(sys),0);
    ASSERT_EQ(alea_setenv("ALEA_METAL_BATCH_PIXELS","127",1),0);
    alea_compute_t* gpu=alea_compute_create(sys,ALEA_COMPUTE_METAL,ALEA_COMPUTE_FP32);
    ASSERT_EQ(alea_unsetenv("ALEA_METAL_BATCH_PIXELS"),0);
    ASSERT_MSG(gpu!=NULL,alea_error());
    render_config_t cfg; render_config_init(&cfg);
    cfg.width=63; cfg.height=47; cfg.log_level=0;
    cfg.eye_set=cfg.target_set=1; cfg.eye[0]=9; cfg.eye[1]=-16; cfg.eye[2]=8;
    cfg.target[0]=cfg.target[1]=cfg.target[2]=0;
    compare_render(sys,gpu,&cfg);
    cfg.num_clips=1; cfg.clips[0]=(render_clip_plane_t){{0,0,1},-.3};
    compare_render(sys,gpu,&cfg);
    cfg.clips[0]=(render_clip_plane_t){{1,.3,.7},.5}; cfg.color_mode=RENDER_COLOR_CELL;
    cfg.aa_samples=2; compare_render(sys,gpu,&cfg);
    cfg.fov=0; cfg.ortho_height=12; cfg.aa_samples=1;
    cfg.num_clips=0; compare_render(sys,gpu,&cfg);
    /* Reuse the same persistent device scene at a larger size, then smaller. */
    cfg.width=87; cfg.height=59; compare_render(sys,gpu,&cfg);
    cfg.width=31; cfg.height=23; cfg.light_set=1; cfg.light_dir[0]=0;
    cfg.light_dir[1]=0; cfg.light_dir[2]=-1; compare_render(sys,gpu,&cfg);
    render_camera_t cam; ASSERT_EQ(render_camera_setup(&cam,&cfg,sys),0);
    render_framebuffer_t* fb=render_framebuffer_create(cfg.width,cfg.height,1);
    ASSERT_NOT_NULL(fb);
    cfg.shadows=1; ASSERT_EQ(alea_compute_render(gpu,&cfg,&cam,fb),-1);
    cfg.shadows=0; cfg.num_clips=2; ASSERT_EQ(alea_compute_render(gpu,&cfg,&cam,fb),-1);
    cfg.num_clips=0; cfg.aa_samples=5; ASSERT_EQ(alea_compute_render(gpu,&cfg,&cam,fb),-1);
    cfg.aa_samples=1; cam.eye[0]=NAN; ASSERT_EQ(alea_compute_render(gpu,&cfg,&cam,fb),-1);
    render_framebuffer_free(fb); render_config_free(&cfg);
    alea_compute_destroy(gpu); alea_destroy(sys);
}

typedef struct {size_t last,calls;int cancel,failed;} progressive_check;
static int check_progress(void* user,const render_framebuffer_t* fb,size_t completed,size_t total) {
    progressive_check* check=user;
    check->failed|=completed<check->last || total!=(size_t)fb->width*fb->height;
    size_t written=0;for(size_t i=0;i<total;++i)written+=fb->cell_id[i]!=-777;
    check->failed|=written!=completed;
    check->last=completed;++check->calls;
    return check->failed || check->cancel==2 || (check->cancel==1 && completed>0);
}
TEST(metal_progressive_tiles_cancel_and_reuse) {
    alea_system_t* sys=sphere_scene();
    alea_compute_t* gpu=alea_compute_create(sys,ALEA_COMPUTE_METAL,ALEA_COMPUTE_FP32);ASSERT_MSG(gpu!=NULL,alea_error());
    render_config_t cfg;test_camera(&cfg);cfg.width=131;cfg.height=97;
    render_framebuffer_t* full=render_framebuffer_create(cfg.width,cfg.height,1);
    render_framebuffer_t* tiled=render_framebuffer_create(cfg.width,cfg.height,1);ASSERT_NOT_NULL(full);ASSERT_NOT_NULL(tiled);
    size_t n=(size_t)cfg.width*cfg.height;
    for(int mode=0;mode<3;++mode) {
        cfg.aa_samples=mode==1?2:1;
        cfg.num_clips=mode==1;cfg.clips[0]=(render_clip_plane_t){{0,1,0},0};
        if(mode==2){cfg.fov=0;cfg.ortho_height=10;}
        render_camera_t cam;ASSERT_EQ(render_camera_setup(&cam,&cfg,sys),0);
        ASSERT_EQ(alea_compute_render(gpu,&cfg,&cam,full),0);
        for(int cancel=2;cancel>=0;--cancel) {
            for(size_t i=0;i<n;++i)tiled->cell_id[i]=-777;
            progressive_check check={.cancel=cancel};
            ASSERT_EQ(alea_compute_render_progressive(gpu,&cfg,&cam,tiled,check_progress,&check),cancel?1:0);
            ASSERT_FALSE(check.failed);ASSERT_FALSE(alea_interrupted());
            if(cancel==2){ASSERT_EQ(check.calls,1);ASSERT_EQ(check.last,0);}
            else if(cancel==1){ASSERT_TRUE(check.last>0 && check.last<n);}
            else {
                ASSERT_EQ(check.last,n);ASSERT_TRUE(check.calls>2);
                ASSERT_EQ(memcmp(full->color,tiled->color,n*3*sizeof(float)),0);
                ASSERT_EQ(memcmp(full->depth,tiled->depth,n*sizeof(float)),0);
                ASSERT_EQ(memcmp(full->normal,tiled->normal,n*3*sizeof(float)),0);
                ASSERT_EQ(memcmp(full->cell_id,tiled->cell_id,n*sizeof(int)),0);
                ASSERT_EQ(memcmp(full->material_id,tiled->material_id,n*sizeof(int)),0);
            }
        }
        ASSERT_EQ(alea_compute_render_progressive(gpu,&cfg,&cam,tiled,NULL,NULL),-1);
    }
    render_framebuffer_free(full);render_framebuffer_free(tiled);render_config_free(&cfg);alea_compute_destroy(gpu);alea_destroy(sys);
}

TEST(metal_render_empty_and_camera_inside) {
    alea_system_t* sys=alea_create();
    alea_compute_t* gpu=alea_compute_create(sys,ALEA_COMPUTE_METAL,ALEA_COMPUTE_FP32);
    ASSERT_MSG(gpu!=NULL,alea_error());
    render_config_t cfg; render_config_init(&cfg);
    cfg.width=17; cfg.height=13; cfg.eye_set=cfg.target_set=1;
    cfg.eye[0]=cfg.eye[2]=0; cfg.eye[1]=-1;
    cfg.target[0]=cfg.target[1]=cfg.target[2]=0;
    compare_render(sys,gpu,&cfg);
    alea_compute_destroy(gpu); alea_destroy(sys);
    sys=sphere_scene();
    ASSERT_EQ(alea_prepare_query_acceleration(sys),0);
    gpu=alea_compute_create(sys,ALEA_COMPUTE_METAL,ALEA_COMPUTE_FP32);
    ASSERT_MSG(gpu!=NULL,alea_error());
    compare_render(sys,gpu,&cfg);
    alea_compute_destroy(gpu); alea_destroy(sys); render_config_free(&cfg);
}

static alea_node_id_t balanced_union(alea_system_t* sys, alea_node_id_t* nodes, int count) {
    if (count==1) return nodes[0];
    int half=count/2;
    return alea_union(sys,balanced_union(sys,nodes,half),balanced_union(sys,nodes+half,count-half));
}

TEST(metal_render_more_than_64_leaves) {
    alea_system_t* sys=alea_create();
    alea_node_id_t leaves[65];
    for (int i=0;i<65;++i) {
        int surface=alea_sphere_surface(sys,0,i*3,0,0,1);
        leaves[i]=alea_halfspace(sys,surface,-1);
    }
    int mat=alea_add_material(sys,1);
    alea_add_cell(sys,1,balanced_union(sys,leaves,65),mat,-1,0);
    alea_compute_t* gpu=alea_compute_create(sys,ALEA_COMPUTE_METAL,ALEA_COMPUTE_FP32);
    ASSERT_MSG(gpu!=NULL,alea_error());
    render_config_t cfg; render_config_init(&cfg);
    cfg.width=cfg.height=1; cfg.eye_set=cfg.target_set=1;
    cfg.eye[0]=0; cfg.eye[1]=-5; cfg.eye[2]=0;
    cfg.target[0]=cfg.target[1]=cfg.target[2]=0;
    render_camera_t cam; ASSERT_EQ(render_camera_setup(&cam,&cfg,sys),0);
    render_framebuffer_t* fb=render_framebuffer_create(1,1,1); ASSERT_NOT_NULL(fb);
    ASSERT_MSG(alea_compute_render(gpu,&cfg,&cam,fb)==0,alea_error());
    ASSERT_EQ(fb->cell_id[0],1); ASSERT_NEAR(fb->depth[0],4,1e-5);
    alea_compute_destroy(gpu); render_framebuffer_free(fb); alea_destroy(sys);
    render_config_free(&cfg);
}

TEST(metal_short_circuit_skips_nonfinite_branch) {
    for(int op=0;op<3;++op) {
        alea_system_t* sys=alea_create();int m=alea_add_material(sys,1);
        int a=alea_sphere_surface(sys,1,0,0,0,2),b=alea_sphere_surface(sys,2,0,0,0,1e20);
        alea_node_id_t left=alea_halfspace(sys,a,op==0?-1:1),right=alea_halfspace(sys,b,-1);
        alea_node_id_t root=op==0?alea_union(sys,left,right):op==1?alea_intersection(sys,left,right):alea_difference(sys,left,right);
        alea_add_cell(sys,1,root,m,-1,0);
        alea_compute_t* gpu=alea_compute_create(sys,ALEA_COMPUTE_METAL,ALEA_COMPUTE_FP32);
        ASSERT_MSG(gpu!=NULL,alea_error());
        double p[3]={0};alea_compute_result_t r;
        ASSERT_MSG(alea_compute_classify(gpu,p,1,&r)==0,alea_error());
        ASSERT_EQ(r.owner_count,op==0?1:0);
        alea_compute_destroy(gpu);alea_destroy(sys);
    }
}

static void test_camera(render_config_t* cfg) {
    render_config_init(cfg);cfg->width=39;cfg->height=31;cfg->log_level=0;
    cfg->eye_set=cfg->target_set=1;cfg->eye[0]=9;cfg->eye[1]=-16;cfg->eye[2]=8;
    cfg->target[0]=cfg->target[1]=cfg->target[2]=0;
}

/* Rounded ITER cylinder coefficients have a tiny nonzero axial eigenvalue.
 * Solving that direction moves the packing origin about 1e8 units away and
 * creates false owners at points well away from the actual boundary. */
TEST(metal_nearly_singular_quadric_classification) {
    const double tail[2][4]={{-160.6168522,-65.2955937,-471.2589331,63032.6198072},
                             {-150.6793781,-62.6328556,-443.0781518,55732.4019059}};
    const double origin[3]={863.57041810763246,152.27170774538098,207.94093005024814};
    for(int q=0;q<2;++q) {
        alea_system_t* sys=alea_create();
        int surface=alea_quadric_surface(sys,1,.1767133,.9408907,.8823960,
            -.4411980,.1667512,.6223240,tail[q][0],tail[q][1],tail[q][2],tail[q][3]);
        alea_node_id_t root=alea_halfspace(sys,surface,-1);
        alea_add_cell(sys,42,root,alea_add_material(sys,7),-1,0);
        alea_compute_t* gpu=alea_compute_create(sys,ALEA_COMPUTE_METAL,ALEA_COMPUTE_FP32);
        ASSERT_MSG(gpu!=NULL,alea_error());
        double xyz[21];alea_compute_result_t result[7];
        for(int i=0;i<7;++i) {
            xyz[3*i]=origin[0];xyz[3*i+1]=origin[1];
            xyz[3*i+2]=origin[2]+(i-3)*10;
        }
        xyz[18]=0;xyz[19]=11.1292181;xyz[20]=q==0?265.982117:250.013794;
        ASSERT_EQ(alea_compute_classify(gpu,xyz,7,result),0);
        int inside=0;
        for(int i=0;i<7;++i) {
            int expected=alea_point_inside(sys,root,xyz[3*i],xyz[3*i+1],xyz[3*i+2]);
            ASSERT_EQ(result[i].owner_count,expected);inside+=expected;
        }
        ASSERT_TRUE(inside>0 && inside<7);
        alea_compute_destroy(gpu);alea_destroy(sys);
    }
}

TEST(metal_render_quadrics_cones_and_tori) {
    for(int kind=0;kind<11;++kind) {
        alea_system_t* sys=alea_create();int m=alea_add_material(sys,7),s;
        switch(kind) {
        case 0:s=alea_quadric_surface(sys,1,1,2,3,.2,.3,.4,0,0,0,-9);break;
        case 1:s=alea_cone_x_surface_sheet(sys,1,0,0,0,.6,1);break;
        case 2:s=alea_cone_y_surface_sheet(sys,1,0,0,0,.6,-1);break;
        case 3:s=alea_cone_z_surface(sys,1,0,0,0,.6);break;
        case 4:s=alea_torus_x_surface(sys,1,0,0,0,3,.7);break;
        case 5:s=alea_torus_y_surface(sys,1,0,0,0,3,.7);break;
        case 6:s=alea_torus_z_surface(sys,1,0,0,0,3,.7);break;
        case 7:s=alea_rcc_surface(sys,1,-1,-2,-3,2,4,6,1);break;
        case 8:s=alea_quadric_surface(sys,1,1,2,0,0,0,0,0,0,-1,0);break;
        case 9:s=alea_torus_z_surface(sys,1,0,0,0,3,.015);break;
        default:s=alea_quadric_surface(sys,1,1,1,1,0,0,0,-20000,0,0,100000000-4);break;
        }
        alea_node_id_t root=alea_halfspace(sys,s,-1);
        if((kind>=1 && kind<=3) || kind==8) {
            int box=alea_box_surface(sys,2,-3,3,-3,3,-3,3);
            root=alea_intersection(sys,root,alea_halfspace(sys,box,-1));
        }
        alea_add_cell(sys,42,root,m,-1,0);
        ASSERT_EQ(alea_prepare_query_acceleration(sys),0);
        alea_compute_t* gpu=alea_compute_create(sys,ALEA_COMPUTE_METAL,ALEA_COMPUTE_FP32);
        ASSERT_MSG(gpu!=NULL,alea_error());
        render_config_t cfg;test_camera(&cfg);
        if(kind==10){cfg.eye[0]+=10000;cfg.target[0]=10000;}
        compare_render(sys,gpu,&cfg);
        cfg.num_clips=1;cfg.clips[0]=(render_clip_plane_t){{0,1,0},-.13};
        compare_render(sys,gpu,&cfg);
        render_config_free(&cfg);alea_compute_destroy(gpu);alea_destroy(sys);
    }
}

TEST(metal_render_nested_transforms_and_parent_csg) {
    alea_system_t* sys=alea_create();int m=alea_add_material(sys,7);
    int box=alea_box_surface(sys,1,-6,6,-3,3,-3,3);
    int hole=alea_cylinder_y_surface(sys,2,-3,0,.8);
    int parent=alea_add_cell(sys,10,alea_difference(sys,alea_halfspace(sys,box,-1),alea_halfspace(sys,hole,-1)),ALEA_MATERIAL_VOID,0,0);
    ASSERT_EQ(alea_set_cell_fill(sys,parent,1,0),0);
    for(int i=0;i<2;++i) {
        double x=i?3:-3;
        int bound=alea_sphere_surface(sys,3+i,x,0,0,2.6);
        int instance=alea_add_cell(sys,20+i,alea_halfspace(sys,bound,-1),ALEA_MATERIAL_VOID,0,1);
        double tr[12]={x,0,0,0,1,0,-1,0,0,0,0,1};
        ASSERT_EQ(alea_add_transform(sys,70+i,tr,12,0),0);
        ASSERT_EQ(alea_set_cell_fill(sys,instance,2,70+i),0);
    }
    int ellipsoid=alea_quadric_surface(sys,5,1,2,3,.2,.3,.4,0,0,0,-12);
    alea_add_cell(sys,42,alea_halfspace(sys,ellipsoid,-1),m,-1,2);
    ASSERT_EQ(alea_prepare_query_acceleration(sys),0);
    alea_compute_t* gpu=alea_compute_create(sys,ALEA_COMPUTE_METAL,ALEA_COMPUTE_FP32);
    ASSERT_MSG(gpu!=NULL,alea_error());
    /* Independently expand two instances into world-space quadrics for the
     * reference. This also checks fill transforms and all ancestor masks. */
    alea_system_t* reference=alea_create();int rm=alea_add_material(reference,7);
    int rb=alea_box_surface(reference,1,-6,6,-3,3,-3,3);
    int rh=alea_cylinder_y_surface(reference,2,-3,0,.8);
    alea_node_id_t mask=alea_difference(reference,alea_halfspace(reference,rb,-1),alea_halfspace(reference,rh,-1));
    alea_node_id_t parts[2];
    for(int i=0;i<2;++i) {
        double x=i?3:-3;
        int q=alea_quadric_surface(reference,3+i,2,1,3,-.2,.4,-.3,-4*x,.2*x,.3*x,2*x*x-12);
        int bound=alea_sphere_surface(reference,5+i,x,0,0,2.6);
        parts[i]=alea_intersection(reference,alea_halfspace(reference,q,-1),alea_halfspace(reference,bound,-1));
    }
    alea_add_cell(reference,42,alea_intersection(reference,mask,alea_union(reference,parts[0],parts[1])),rm,-1,0);
    ASSERT_EQ(alea_prepare_query_acceleration(reference),0);
    render_config_t cfg;test_camera(&cfg);
    compare_render(reference,gpu,&cfg);
    cfg.num_clips=1;cfg.clips[0]=(render_clip_plane_t){{0,1,0},0};compare_render(reference,gpu,&cfg);
    cfg.eye[0]=-3;cfg.eye[1]=-.5;cfg.eye[2]=0;cfg.target[0]=-3;cfg.target[1]=2;
    cfg.num_clips=0;compare_render(reference,gpu,&cfg);
    render_config_free(&cfg);alea_compute_destroy(gpu);alea_destroy(sys);alea_destroy(reference);
}

TEST(metal_rejects_missing_and_cyclic_fills) {
    alea_system_t* sys=sphere_scene();sys->cells.data[0].fill_universe=1;
    ASSERT_NULL(alea_compute_create(sys,ALEA_COMPUTE_METAL,ALEA_COMPUTE_FP32));
    alea_add_cell(sys,43,sys->cells.data[0].root_node_id,0,0,1);
    sys->cells.data[1].fill_universe=1;
    ASSERT_NULL(alea_compute_create(sys,ALEA_COMPUTE_METAL,ALEA_COMPUTE_FP32));
    alea_destroy(sys);
}

static void check_rectangular_lattices(void) {
    for(int transformed=0;transformed<2;++transformed)for(int repeating=0;repeating<2;++repeating)for(int zero=0;zero<2;++zero) {
        alea_system_t* sys=alea_create();int m=alea_add_material(sys,7);
        double dx=transformed?1:0,dy=transformed?.5:0,dz=transformed?.25:0;
        int outer=alea_box_surface(sys,1,-3+dx,3+dx,-3+dy,3+dy,-1+dz,1+dz);
        int parent=alea_add_cell(sys,10,alea_halfspace(sys,outer,-1),ALEA_MATERIAL_VOID,0,0);
        if(transformed){double tr[12]={dx,dy,dz,0,1,0,-1,0,0,0,0,1};ASSERT_EQ(alea_add_transform(sys,70,tr,12,0),0);}
        ASSERT_EQ(alea_set_cell_fill(sys,parent,1,transformed?70:0),0);
        int unit=alea_box_surface(sys,2,-1,1,-1,1,-1,1);
        int lat=alea_add_cell(sys,20,alea_halfspace(sys,unit,-1),ALEA_MATERIAL_VOID,0,1);
        int sphere=alea_sphere_surface(sys,3,0,0,0,.7);
        alea_add_cell(sys,42,alea_halfspace(sys,sphere,-1),m,-1,2);
        alea_cell_entry_t* c=&sys->cells.data[lat];c->lat_type=1;c->lat_fill_repeating=repeating;c->lat_fill_zero_element_coords=zero;
        for(int a=0;a<2;++a){c->lat_fill_dims[2*a]=repeating?0:-1;c->lat_fill_dims[2*a+1]=repeating?0:1;c->lat_pitch[a]=2;c->lat_lower_left[a]=repeating?-1:-3;}
        c->lat_pitch[2]=repeating?0:2;c->lat_lower_left[2]=-1;
        c->lat_fill_count=repeating?1:9;c->lat_fill=malloc(c->lat_fill_count*sizeof(int));
        ASSERT_NOT_NULL(c->lat_fill);for(size_t j=0;j<c->lat_fill_count;++j)c->lat_fill[j]=2;
        if(!repeating && zero)c->lat_fill[4]=0;
        alea_compute_t* gpu=alea_compute_create(sys,ALEA_COMPUTE_METAL,ALEA_COMPUTE_FP32);ASSERT_MSG(gpu!=NULL,alea_error());
        double point[]={0,0,0};alea_compute_result_t result;
        ASSERT_EQ(alea_compute_classify(gpu,point,1,&result),-1);
        ASSERT_EQ(alea_error_code(),ALEA_ERR_UNSUPPORTED);
        alea_system_t* reference=alea_create();int rm=alea_add_material(reference,7);
        alea_node_id_t parts[9];int n=0;
        for(int x=-1;x<=1;++x)for(int y=-1;y<=1;++y) {
            if(!repeating && zero && x==0 && y==0)continue;
            int s=alea_sphere_surface(reference,0,(transformed?-y:x)*2+dx,(transformed?x:y)*2+dy,dz,.7);parts[n++]=alea_halfspace(reference,s,-1);
        }
        alea_add_cell(reference,42,balanced_union(reference,parts,n),rm,-1,0);
        ASSERT_EQ(alea_prepare_query_acceleration(reference),0);
        render_config_t cfg;test_camera(&cfg);compare_render(reference,gpu,&cfg);
        cfg.num_clips=1;cfg.clips[0]=(render_clip_plane_t){{1,0,0},-.2};compare_render(reference,gpu,&cfg);
        cfg.num_clips=0;cfg.eye[0]=0;cfg.eye[1]=0;cfg.eye[2]=8;compare_render(reference,gpu,&cfg);
        render_config_free(&cfg);alea_compute_destroy(gpu);alea_destroy(sys);alea_destroy(reference);
    }
}

TEST(metal_rectangular_lattices_match_expanded_geometry) {check_rectangular_lattices();}

TEST(metal_elliptical_torus_axes_and_normals) {
    for(int axis=0;axis<3;++axis) {
        alea_system_t* sys=alea_create();int m=alea_add_material(sys,7);
        int s=axis==0?alea_torus_x_surface(sys,1,0,0,0,3,.7):
            axis==1?alea_torus_y_surface(sys,1,0,0,0,3,.7):alea_torus_z_surface(sys,1,0,0,0,3,.7);
        alea_node_id_t root=alea_halfspace(sys,s,-1);
        /* Exercise MCNP's elliptical tube payload; the public convenience
         * constructor creates only circular tubes. */
        uint32_t id=sys->nodes.data[root].primitive.primitive_id;
        alea_torus_data_t* torus=(alea_torus_data_t*)alea_primitive_payload_const(sys,id);
        torus->axial_semiwidth_B=1.4;
        sys->nodes.data[root].bbox=(alea_node_bbox_t){-4,4,-4,4,-4,4};
        alea_add_cell(sys,42,root,m,-1,0);
        alea_compute_t* gpu=alea_compute_create(sys,ALEA_COMPUTE_METAL,ALEA_COMPUTE_FP32);ASSERT_MSG(gpu!=NULL,alea_error());
        render_config_t cfg;test_camera(&cfg);cfg.width=cfg.height=1;
        memset(cfg.eye,0,sizeof(cfg.eye));memset(cfg.target,0,sizeof(cfg.target));
        int radial=(axis+1)%3;cfg.eye[axis]=10;cfg.eye[radial]=cfg.target[radial]=3;
        render_camera_t cam;ASSERT_EQ(render_camera_setup(&cam,&cfg,sys),0);
        render_framebuffer_t* fb=render_framebuffer_create(1,1,1);ASSERT_NOT_NULL(fb);
        ASSERT_MSG(alea_compute_render(gpu,&cfg,&cam,fb)==0,alea_error());
        ASSERT_EQ(fb->cell_id[0],42);ASSERT_NEAR(fb->depth[0],8.6,1e-4);
        for(int a=0;a<3;++a)ASSERT_NEAR(fb->normal[a],a==axis?1:0,1e-4);
        memset(cfg.eye,0,sizeof(cfg.eye));memset(cfg.target,0,sizeof(cfg.target));cfg.eye[radial]=10;
        ASSERT_EQ(render_camera_setup(&cam,&cfg,sys),0);
        ASSERT_MSG(alea_compute_render(gpu,&cfg,&cam,fb)==0,alea_error());
        ASSERT_EQ(fb->cell_id[0],42);ASSERT_NEAR(fb->depth[0],6.3,1e-4);
        for(int a=0;a<3;++a)ASSERT_NEAR(fb->normal[a],a==radial?1:0,1e-4);
        render_framebuffer_free(fb);render_config_free(&cfg);alea_compute_destroy(gpu);alea_destroy(sys);
    }
}

TEST(metal_right_deep_short_circuit_program) {
    alea_system_t* sys=alea_create();int m=alea_add_material(sys,1);
    alea_node_id_t root=ALEA_NODE_ID_INVALID;
    for(int i=0;i<600;++i) {
        int s=alea_sphere_surface(sys,0,i*3,0,0,1);
        alea_node_id_t leaf=alea_halfspace(sys,s,-1);
        root=i?alea_union(sys,leaf,root):leaf;
    }
    alea_add_cell(sys,1,root,m,-1,0);
    double points[]={0,0,0, 900,0,0, 1797,0,0, 1.5,0,0, -2,0,0};
    check_scene(sys,points,5,ALEA_COMPUTE_METAL);alea_destroy(sys);
}

/* The spatial experiment must preserve closed-set boundary semantics even
 * when folding a complement or either ordering of a difference. */
static float local_test_value(const compute_node* code,size_t count,
                               const compute_primitive* primitives,float x) {
    float values[COMPUTE_STACK];unsigned sp=0;
    for(size_t i=0;i<count;++i) {
        compute_node n=code[i];
        if(n.op==COMPUTE_LEAF) {
            const compute_primitive* p=&primitives[n.arg>>1];
            float v=p->p[0]*x+p->p[3];values[sp++]=(n.arg&1)?-v:v;
        } else if(n.op==COMPUTE_CONSTANT)values[sp++]=n.arg?1:-1;
        else if(n.op==COMPUTE_COMPLEMENT)values[sp-1]=-values[sp-1];
        else if(n.op==COMPUTE_SKIP_UNION || n.op==COMPUTE_SKIP_INTERSECTION || n.op==COMPUTE_SKIP_REVERSE_DIFFERENCE)continue;
        else {
            float b=values[--sp],a=values[sp-1];
            values[sp-1]=n.op==COMPUTE_UNION?fminf(a,b):n.op==COMPUTE_INTERSECTION?fmaxf(a,b):
                n.op==COMPUTE_DIFFERENCE?fmaxf(a,-b):fmaxf(b,-a);
        }
    }
    return values[0];
}
TEST(metal_local_boolean_boundaries) {
    float bounds[]={-1,1,-1,1,-1,1};
    const unsigned ops[]={COMPUTE_UNION,COMPUTE_INTERSECTION,COMPUTE_DIFFERENCE,COMPUTE_REVERSE_DIFFERENCE};
    for(int left=-1;left<=1;++left)for(int right=-1;right<=1;++right)
    for(unsigned op=0;op<4;++op)for(unsigned flip=0;flip<4;++flip)for(int complement=0;complement<2;++complement) {
        compute_primitive primitives[2]={{.type=COMPUTE_PLANE,.p={1,0,0,left*10}},
                                          {.type=COMPUTE_PLANE,.p={1,0,0,right*10}}};
        compute_node code[]={{COMPUTE_LEAF,flip&1},{COMPUTE_LEAF,2|((flip>>1)&1)},
                             {ops[op],0},{COMPUTE_COMPLEMENT,0}};
        compute_node reduced[8];size_t count=0;
        ASSERT_EQ(alea_compute_local_simplify(code,3+complement,primitives,bounds,reduced,&count),0);
        for(int j=-2;j<=2;++j) {
            float a=local_test_value(code,3+complement,primitives,j*.5f);
            float b=local_test_value(reduced,count,primitives,j*.5f);
            ASSERT_EQ((a>0)-(a<0),(b>0)-(b<0));
        }
    }
}

static alea_node_id_t local_spheres(alea_system_t* sys) {
    alea_node_id_t root=ALEA_NODE_ID_INVALID;
    for(int x=0;x<4;++x)for(int y=0;y<4;++y)for(int z=0;z<4;++z) {
        int sphere=alea_sphere_surface(sys,0,x*3-4.5,y*3-4.5,z*3-4.5,1);
        alea_node_id_t n=alea_halfspace(sys,sphere,-1);root=root==ALEA_NODE_ID_INVALID?n:alea_union(sys,root,n);
    }
    return root;
}
static alea_system_t* local_render_scene(void) {
    alea_system_t* sys=alea_create();int mat=alea_add_material(sys,7);
    alea_node_id_t root=local_spheres(sys);
    alea_add_cell(sys,42,root,mat,-1,0);return sys;
}
TEST(metal_spatial_disconnected_parent) {
    alea_system_t* sys=local_render_scene();
    sys->cells.data[0].material_id=0;
    ASSERT_EQ(alea_set_cell_fill(sys,0,1,0),0);
    int box=alea_box_surface(sys,0,-5,5,-2,2,-5,5);
    int mat=alea_add_material(sys,8);
    alea_add_cell(sys,43,alea_halfspace(sys,box,-1),mat,-1,1);
    ASSERT_EQ(alea_prepare_query_acceleration(sys),0);
    ASSERT_EQ(alea_setenv("ALEA_METAL_SPATIAL","0",1),0);
    alea_compute_t* baseline=alea_compute_create(sys,ALEA_COMPUTE_METAL,ALEA_COMPUTE_FP32);ASSERT_NOT_NULL(baseline);
    /* Flatten the parent mask and child solid into an independent reference;
     * the CPU nested segment path can miss these disconnected parent pieces. */
    alea_system_t* reference=alea_create();int rm=alea_add_material(reference,8);
    alea_node_id_t mask=local_spheres(reference);
    int rb=alea_box_surface(reference,0,-5,5,-2,2,-5,5);
    alea_add_cell(reference,43,alea_intersection(reference,mask,alea_halfspace(reference,rb,-1)),rm,-1,0);
    ASSERT_EQ(alea_prepare_query_acceleration(reference),0);
    ASSERT_EQ(alea_unsetenv("ALEA_METAL_SPATIAL"),0);
    alea_compute_t* gpu=alea_compute_create(sys,ALEA_COMPUTE_METAL,ALEA_COMPUTE_FP32);ASSERT_NOT_NULL(gpu);
    render_config_t cfg;render_config_init(&cfg);cfg.width=32;cfg.height=24;cfg.log_level=0;
    cfg.eye_set=cfg.target_set=1;cfg.target[0]=cfg.target[1]=cfg.target[2]=0;
    cfg.eye[0]=9;cfg.eye[1]=-20;cfg.eye[2]=8;compare_render(reference,baseline,&cfg);compare_render(reference,gpu,&cfg);
    cfg.num_clips=1;cfg.clips[0]=(render_clip_plane_t){{0,1,0},-1.5};compare_render(reference,baseline,&cfg);compare_render(reference,gpu,&cfg);
    render_config_free(&cfg);alea_compute_destroy(baseline);alea_compute_destroy(gpu);alea_destroy(sys);alea_destroy(reference);
}
TEST(metal_spatial_shared_content_and_original_programs) {
    alea_system_t* sys=alea_create();int mat=alea_add_material(sys,7);
    double tr[]={20,0,0};ASSERT_EQ(alea_add_transform(sys,70,tr,3,0),0);
    for(int i=0;i<3;++i) {
        double x=i==0?0:i==1?20:10;
        int box=alea_box_surface(sys,0,x-5,x+5,-5,5,-5,5);
        int cell=alea_add_cell(sys,10+i,alea_halfspace(sys,box,-1),0,0,0);
        ASSERT_EQ(alea_set_cell_fill(sys,cell,1,i==1?70:0),0);
    }
    int sphere=alea_sphere_surface(sys,0,0,0,0,2);
    alea_add_cell(sys,42,alea_halfspace(sys,sphere,-1),mat,-1,1);
    ASSERT_EQ(alea_prepare_query_acceleration(sys),0);
    compute_scene original,accelerated;
    ASSERT_EQ(alea_setenv("ALEA_METAL_SPATIAL","0",1),0);
    ASSERT_EQ(alea_compute_pack(sys,&original),0);
    ASSERT_EQ(alea_unsetenv("ALEA_METAL_SPATIAL"),0);
    ASSERT_EQ(alea_compute_pack(sys,&accelerated),0);
    ASSERT_EQ(memcmp(original.nodes,accelerated.nodes,original.node_count*sizeof(compute_node)),0);
    ASSERT_EQ(accelerated.cells[0].child,accelerated.cells[1].child);
    ASSERT_TRUE(accelerated.cells[0].spatial!=COMPUTE_NONE);
    compute_spatial_program* proof=&accelerated.spatial[accelerated.cells[0].spatial];
    ASSERT_EQ(proof->count,1);ASSERT_EQ(accelerated.nodes[proof->root].op,COMPUTE_CONSTANT);
    ASSERT_EQ(accelerated.nodes[proof->root].arg,0);
    int bounded=0,empty=0;
    for(size_t i=0;i<accelerated.render.bvh_count;++i) {
        compute_render_bvh* b=&accelerated.render.bvh[i];
        if(b->cell==0){ASSERT_TRUE(b->bounds[1]<2.1);ASSERT_TRUE(b->fallback[1]>=5);++bounded;}
        if(b->cell==2){ASSERT_TRUE(b->bounds[0]>b->bounds[1]);++empty;}
    }
    ASSERT_EQ(bounded,1);ASSERT_EQ(empty,1);
    alea_compute_scene_free(&original);alea_compute_scene_free(&accelerated);
    alea_system_t* reference=alea_create();int rm=alea_add_material(reference,7);
    int a=alea_sphere_surface(reference,0,0,0,0,2),b=alea_sphere_surface(reference,0,20,0,0,2);
    alea_add_cell(reference,42,alea_union(reference,alea_halfspace(reference,a,-1),alea_halfspace(reference,b,-1)),rm,-1,0);
    ASSERT_EQ(alea_prepare_query_acceleration(reference),0);
    alea_compute_t* gpu=alea_compute_create(sys,ALEA_COMPUTE_METAL,ALEA_COMPUTE_FP32);ASSERT_NOT_NULL(gpu);
    render_config_t cfg;render_config_init(&cfg);cfg.width=41;cfg.height=27;cfg.log_level=0;
    cfg.eye_set=cfg.target_set=1;cfg.eye[0]=10;cfg.eye[1]=-40;cfg.eye[2]=10;cfg.target[0]=10;
    int previous=compare_resume;compare_resume=1;
    compare_render(reference,gpu,&cfg);
    cfg.num_clips=1;cfg.clips[0]=(render_clip_plane_t){{0,1,0},0};compare_render(reference,gpu,&cfg);
    compare_resume=previous;render_config_free(&cfg);alea_compute_destroy(gpu);alea_destroy(reference);alea_destroy(sys);
}

TEST(metal_spatial_hint_is_not_a_clip) {
    alea_system_t* sys=alea_create();int mat=alea_add_material(sys,7);
    int near=alea_sphere_surface(sys,0,5,0,0,1),far=alea_sphere_surface(sys,0,1000,0,0,1);
    alea_node_id_t mask=alea_union(sys,alea_halfspace(sys,near,-1),alea_halfspace(sys,far,-1));
    int parent=alea_add_cell(sys,10,mask,0,0,0);ASSERT_EQ(alea_set_cell_fill(sys,parent,1,0),0);
    int box=alea_box_surface(sys,0,-2000,2000,-2000,2000,-2000,2000);
    alea_node_id_t solid=alea_halfspace(sys,box,-1);alea_add_cell(sys,42,solid,mat,-1,1);
    int hint=alea_box_surface(sys,0,-1,1,-1,1,-1,1);alea_add_cell(sys,43,alea_halfspace(sys,hint,-1),mat,-1,0);
    ASSERT_EQ(alea_prepare_query_acceleration(sys),0);
    /* Simulate conservative unbounded source caches: the finite sibling is
     * only a cache-domain hint, and must never hide either distant sphere. */
    alea_bbox_t wide={-1e10,1e10,-1e10,1e10,-1e10,1e10};
    alea_node_bbox_set(&sys->nodes.data[mask].bbox,&wide);alea_node_bbox_set(&sys->nodes.data[solid].bbox,&wide);
    compute_scene packed;ASSERT_EQ(alea_compute_pack(sys,&packed),0);
    ASSERT_TRUE(packed.cells[0].spatial!=COMPUTE_NONE);
    compute_spatial_program p=packed.spatial[packed.cells[0].spatial];
    ASSERT_EQ(p.count,1);ASSERT_EQ(packed.nodes[p.root].op,COMPUTE_CONSTANT);ASSERT_EQ(packed.nodes[p.root].arg,1);
    alea_compute_scene_free(&packed);
    alea_system_t* reference=alea_create();int rm=alea_add_material(reference,7);
    near=alea_sphere_surface(reference,0,5,0,0,1);far=alea_sphere_surface(reference,0,1000,0,0,1);
    alea_add_cell(reference,42,alea_union(reference,alea_halfspace(reference,near,-1),alea_halfspace(reference,far,-1)),rm,-1,0);
    hint=alea_box_surface(reference,0,-1,1,-1,1,-1,1);alea_add_cell(reference,43,alea_halfspace(reference,hint,-1),rm,-1,0);
    ASSERT_EQ(alea_prepare_query_acceleration(reference),0);
    alea_compute_t* gpu=alea_compute_create(sys,ALEA_COMPUTE_METAL,ALEA_COMPUTE_FP32);ASSERT_NOT_NULL(gpu);
    render_config_t cfg;render_config_init(&cfg);cfg.width=25;cfg.height=19;cfg.log_level=0;cfg.eye_set=cfg.target_set=1;
    int previous=compare_resume;compare_resume=1;
    for(int i=0;i<2;++i){cfg.eye[0]=cfg.target[0]=i?1000:5;cfg.eye[1]=-8;cfg.eye[2]=3;compare_render(reference,gpu,&cfg);}
    compare_resume=previous;render_config_free(&cfg);alea_compute_destroy(gpu);alea_destroy(reference);alea_destroy(sys);
}

TEST(metal_spatial_distant_ray_fallback) {
    alea_system_t* sys=alea_create();int mat=alea_add_material(sys,7);
    int outer=alea_box_surface(sys,0,-100,100,-100,100,-100,100);
    int parent=alea_add_cell(sys,10,alea_halfspace(sys,outer,-1),0,0,0);
    double tr[]={20,0,0,.8,.6,0,-.6,.8,0,0,0,1};ASSERT_EQ(alea_add_transform(sys,70,tr,12,0),0);
    ASSERT_EQ(alea_set_cell_fill(sys,parent,1,70),0);
    int inner=alea_box_surface(sys,0,-1,1,-1,1,-1,1);alea_add_cell(sys,42,alea_halfspace(sys,inner,-1),mat,-1,1);
    ASSERT_EQ(alea_prepare_query_acceleration(sys),0);
    ASSERT_EQ(alea_setenv("ALEA_METAL_SPATIAL","0",1),0);
    alea_compute_t* old=alea_compute_create(sys,ALEA_COMPUTE_METAL,ALEA_COMPUTE_FP32);ASSERT_NOT_NULL(old);
    ASSERT_EQ(alea_unsetenv("ALEA_METAL_SPATIAL"),0);
    alea_compute_t* gpu=alea_compute_create(sys,ALEA_COMPUTE_METAL,ALEA_COMPUTE_FP32);ASSERT_NOT_NULL(gpu);
    render_config_t cfg;render_config_init(&cfg);cfg.width=15;cfg.height=11;cfg.log_level=0;cfg.eye_set=cfg.target_set=1;
    cfg.eye[0]=cfg.target[0]=20;cfg.eye[1]=-1e6;
    render_camera_t cam;ASSERT_EQ(render_camera_setup(&cam,&cfg,sys),0);
    render_framebuffer_t* a=render_framebuffer_create(15,11,1),*b=render_framebuffer_create(15,11,1);
    ASSERT_NOT_NULL(a);ASSERT_NOT_NULL(b);ASSERT_EQ(alea_compute_render(old,&cfg,&cam,a),0);
    for(int progressive=0;progressive<2;++progressive) {
        int calls=0;ASSERT_EQ(progressive?alea_compute_render_progressive(gpu,&cfg,&cam,b,count_progress,&calls):alea_compute_render(gpu,&cfg,&cam,b),0);
        ASSERT_EQ(memcmp(a->cell_id,b->cell_id,165*sizeof(int)),0);
        ASSERT_EQ(memcmp(a->material_id,b->material_id,165*sizeof(int)),0);
        ASSERT_EQ(memcmp(a->depth,b->depth,165*sizeof(float)),0);
        ASSERT_EQ(memcmp(a->normal,b->normal,495*sizeof(float)),0);
        ASSERT_EQ(memcmp(a->color,b->color,495*sizeof(float)),0);
    }
    render_framebuffer_free(a);render_framebuffer_free(b);render_config_free(&cfg);alea_compute_destroy(old);alea_compute_destroy(gpu);alea_destroy(sys);
}
TEST(metal_spatial_thin_parent_boundary) {
    alea_system_t* sys=alea_create();int mat=alea_add_material(sys,7);
    int thin=alea_box_surface(sys,0,-.0001,.0001,-2,2,-2,2);
    int parent=alea_add_cell(sys,10,alea_halfspace(sys,thin,-1),0,0,0);ASSERT_EQ(alea_set_cell_fill(sys,parent,1,0),0);
    int sphere=alea_sphere_surface(sys,0,0,0,0,1);alea_add_cell(sys,42,alea_halfspace(sys,sphere,-1),mat,-1,1);
    ASSERT_EQ(alea_prepare_query_acceleration(sys),0);
    alea_system_t* reference=alea_create();int rm=alea_add_material(reference,7);
    thin=alea_box_surface(reference,0,-.0001,.0001,-2,2,-2,2);sphere=alea_sphere_surface(reference,0,0,0,0,1);
    alea_add_cell(reference,42,alea_intersection(reference,alea_halfspace(reference,thin,-1),alea_halfspace(reference,sphere,-1)),rm,-1,0);
    ASSERT_EQ(alea_prepare_query_acceleration(reference),0);
    alea_compute_t* gpu=alea_compute_create(sys,ALEA_COMPUTE_METAL,ALEA_COMPUTE_FP32);ASSERT_NOT_NULL(gpu);
    render_config_t cfg;render_config_init(&cfg);cfg.width=25;cfg.height=19;cfg.log_level=0;cfg.eye_set=cfg.target_set=1;
    cfg.eye[0]=4;cfg.eye[1]=.2;cfg.eye[2]=.3;
    int previous=compare_resume;compare_resume=1;compare_render(reference,gpu,&cfg);
    cfg.num_clips=1;cfg.clips[0]=(render_clip_plane_t){{0,1,0},0};compare_render(reference,gpu,&cfg);
    compare_resume=previous;render_config_free(&cfg);alea_compute_destroy(gpu);alea_destroy(reference);alea_destroy(sys);
}

TEST(metal_spatial_proved_empty_regions) {
    alea_system_t* sys=alea_create(),*reference=alea_create();
    for(int mode=0;mode<2;++mode) {
        alea_system_t* model=mode?reference:sys;int mat=alea_add_material(model,7);
        alea_node_id_t mask=ALEA_NODE_ID_INVALID;
        for(int x=0;x<4;++x)for(int y=0;y<4;++y)for(int z=0;z<4;++z) {
            int sphere=alea_sphere_surface(model,0,5.5+x*3,5.5+y*3,5.5+z*3,1);
            alea_node_id_t part=alea_halfspace(model,sphere,-1);mask=mask==ALEA_NODE_ID_INVALID?part:alea_union(model,mask,part);
        }
        if(mode)alea_add_cell(model,42,mask,mat,-1,0);
        else {
            int parent=alea_add_cell(model,10,mask,0,0,0);ASSERT_EQ(alea_set_cell_fill(model,parent,1,0),0);
            int box=alea_box_surface(model,0,-30,30,-30,30,-30,30);alea_add_cell(model,42,alea_halfspace(model,box,-1),mat,-1,1);
        }
        ASSERT_EQ(alea_prepare_query_acceleration(model),0);
        if(!mode){alea_bbox_t loose={-100,100,-100,100,-100,100};alea_node_bbox_set(&model->nodes.data[mask].bbox,&loose);}
    }
    compute_scene packed;ASSERT_EQ(alea_compute_pack(sys,&packed),0);int checked=0;
    for(size_t i=0;i<packed.render.bvh_count;++i)if(packed.render.bvh[i].cell==0) {
        ASSERT_TRUE(packed.render.bvh[i].bounds[0]>=0);ASSERT_TRUE(packed.render.bvh[i].fallback[0]<=-100);
        ASSERT_TRUE(isfinite(packed.render.bvh[i].guard));++checked;
    }
    ASSERT_EQ(checked,1);alea_compute_scene_free(&packed);
    alea_compute_t* gpu=alea_compute_create(sys,ALEA_COMPUTE_METAL,ALEA_COMPUTE_FP32);ASSERT_NOT_NULL(gpu);
    render_config_t cfg;render_config_init(&cfg);cfg.width=33;cfg.height=25;cfg.log_level=0;cfg.eye_set=cfg.target_set=1;
    cfg.eye[0]=-20;cfg.eye[1]=-25;cfg.eye[2]=20;cfg.target[0]=cfg.target[1]=cfg.target[2]=10;
    int previous=compare_resume;compare_resume=1;compare_render(reference,gpu,&cfg);
    cfg.num_clips=1;cfg.clips[0]=(render_clip_plane_t){{0,1,0},-10};compare_render(reference,gpu,&cfg);
    compare_resume=previous;render_config_free(&cfg);alea_compute_destroy(gpu);alea_destroy(reference);alea_destroy(sys);
}

TEST(metal_spatial_configuration) {
    alea_system_t* sys=sphere_scene();compute_scene packed;
    const char* bad[]={"","yes","-1","2","bounds "};
    for(unsigned i=0;i<sizeof(bad)/sizeof(*bad);++i) {
        ASSERT_EQ(alea_setenv("ALEA_METAL_SPATIAL",bad[i],1),0);
        ASSERT_EQ(alea_compute_pack(sys,&packed),-1);ASSERT_EQ(alea_error_code(),ALEA_ERR_INVALID_ARG);
    }
    ASSERT_EQ(alea_setenv("ALEA_METAL_SPATIAL","bounds",1),0);
    ASSERT_EQ(alea_compute_pack(sys,&packed),0);ASSERT_EQ(packed.spatial_count,0);alea_compute_scene_free(&packed);
    ASSERT_EQ(alea_unsetenv("ALEA_METAL_SPATIAL"),0);alea_destroy(sys);
}

TEST(metal_local_primitive_proofs) {
    /* Known interior/exterior boxes, and tangent or overflowing cases that
     * must remain unresolved. A finite corner sample cannot prove a sphere's
     * entire box outside: this box contains the sphere's center. */
    compute_node leaf={COMPUTE_LEAF,0},out[2];size_t count;
    const float inside[]={-.25,.25,-.25,.25,-.25,.25};
    const float outside[]={3,4,3,4,3,4};
    const float contains[]={-2,2,-2,2,-2,2};
    const float tangent[]={1,2,-.25,.25,-.25,.25};
    compute_primitive p={.type=COMPUTE_SPHERE,.p={0,0,0,1}};
    ASSERT_EQ(alea_compute_local_simplify(&leaf,1,&p,inside,out,&count),0);
    ASSERT_EQ(out[0].op,COMPUTE_CONSTANT);ASSERT_EQ(out[0].arg,0);
    ASSERT_EQ(alea_compute_local_simplify(&leaf,1,&p,outside,out,&count),0);
    ASSERT_EQ(out[0].op,COMPUTE_CONSTANT);ASSERT_EQ(out[0].arg,1);
    ASSERT_EQ(alea_compute_local_simplify(&leaf,1,&p,contains,out,&count),0);ASSERT_EQ(out[0].op,COMPUTE_LEAF);
    ASSERT_EQ(alea_compute_local_simplify(&leaf,1,&p,tangent,out,&count),0);ASSERT_EQ(out[0].op,COMPUTE_LEAF);
    p.p[3]=1e20f;
    ASSERT_EQ(alea_compute_local_simplify(&leaf,1,&p,inside,out,&count),0);ASSERT_EQ(out[0].op,COMPUTE_LEAF);
    p=(compute_primitive){.type=COMPUTE_PLANE,.p={1e-40f,0,0,0}};
    ASSERT_EQ(alea_compute_local_simplify(&leaf,1,&p,outside,out,&count),0);ASSERT_EQ(out[0].op,COMPUTE_LEAF);
    for(unsigned axis=COMPUTE_CYL_X;axis<=COMPUTE_CYL_Z;++axis) {
        p=(compute_primitive){.type=axis,.p={0,0,1}};
        ASSERT_EQ(alea_compute_local_simplify(&leaf,1,&p,inside,out,&count),0);ASSERT_EQ(out[0].op,COMPUTE_CONSTANT);ASSERT_EQ(out[0].arg,0);
        ASSERT_EQ(alea_compute_local_simplify(&leaf,1,&p,outside,out,&count),0);ASSERT_EQ(out[0].op,COMPUTE_CONSTANT);ASSERT_EQ(out[0].arg,1);
    }
    p=(compute_primitive){.type=COMPUTE_BOX,.p={-1,1,-1,1,-1,1}};
    ASSERT_EQ(alea_compute_local_simplify(&leaf,1,&p,inside,out,&count),0);ASSERT_EQ(out[0].op,COMPUTE_CONSTANT);ASSERT_EQ(out[0].arg,0);
    ASSERT_EQ(alea_compute_local_simplify(&leaf,1,&p,outside,out,&count),0);ASSERT_EQ(out[0].op,COMPUTE_CONSTANT);ASSERT_EQ(out[0].arg,1);
    /* An ellipsoid with a cross term; interval dependency may lose pruning,
     * but cannot assert a sign in a region crossing the boundary. */
    p=(compute_primitive){.type=COMPUTE_QUADRIC,.p={1,2,3,.5,0,0,0,0,0,-1}};
    ASSERT_EQ(alea_compute_local_simplify(&leaf,1,&p,inside,out,&count),0);ASSERT_EQ(out[0].op,COMPUTE_CONSTANT);ASSERT_EQ(out[0].arg,0);
    ASSERT_EQ(alea_compute_local_simplify(&leaf,1,&p,outside,out,&count),0);ASSERT_EQ(out[0].op,COMPUTE_CONSTANT);ASSERT_EQ(out[0].arg,1);
    ASSERT_EQ(alea_compute_local_simplify(&leaf,1,&p,contains,out,&count),0);ASSERT_EQ(out[0].op,COMPUTE_LEAF);
    p=(compute_primitive){.type=COMPUTE_TORUS_Z,.p={0,0,0,3,1,2}};
    const float tube[]={2.9,3.1,-.1,.1,-.1,.1};
    ASSERT_EQ(alea_compute_local_simplify(&leaf,1,&p,tube,out,&count),0);ASSERT_EQ(out[0].op,COMPUTE_CONSTANT);ASSERT_EQ(out[0].arg,0);
    ASSERT_EQ(alea_compute_local_simplify(&leaf,1,&p,inside,out,&count),0);ASSERT_EQ(out[0].op,COMPUTE_CONSTANT);ASSERT_EQ(out[0].arg,1);
    const float cone_inside[]={-.1,.1,-.1,.1,2,3};
    p=(compute_primitive){.type=COMPUTE_CONE_Z,.p={0,0,0,1,1}};
    ASSERT_EQ(alea_compute_local_simplify(&leaf,1,&p,cone_inside,out,&count),0);ASSERT_EQ(out[0].op,COMPUTE_CONSTANT);ASSERT_EQ(out[0].arg,0);
    p.p[4]=-1;
    ASSERT_EQ(alea_compute_local_simplify(&leaf,1,&p,cone_inside,out,&count),0);ASSERT_EQ(out[0].op,COMPUTE_CONSTANT);ASSERT_EQ(out[0].arg,1);
    p=(compute_primitive){.type=COMPUTE_RCC,.p={0,0,-1,0,0,1,2,1}};
    ASSERT_EQ(alea_compute_local_simplify(&leaf,1,&p,inside,out,&count),0);ASSERT_EQ(out[0].op,COMPUTE_CONSTANT);ASSERT_EQ(out[0].arg,0);
    ASSERT_EQ(alea_compute_local_simplify(&leaf,1,&p,outside,out,&count),0);ASSERT_EQ(out[0].op,COMPUTE_CONSTANT);ASSERT_EQ(out[0].arg,1);
}
TEST(metal_profile_counts_and_reuse) {
    char directory[]="/tmp/alea-metal-profile-XXXXXX";ASSERT_NOT_NULL(mkdtemp(directory));
    char prefix[256],path[300];snprintf(prefix,sizeof(prefix),"%s/frame",directory);
    ASSERT_EQ(alea_setenv("ALEA_METAL_PROFILE",prefix,1),0);
    alea_system_t* sys=sphere_scene();
    ASSERT_EQ(alea_prepare_query_acceleration(sys),0);
    alea_compute_t* gpu=alea_compute_create(sys,ALEA_COMPUTE_METAL,ALEA_COMPUTE_FP32);ASSERT_NOT_NULL(gpu);
    ASSERT_EQ(alea_unsetenv("ALEA_METAL_PROFILE"),0);
    render_config_t cfg;render_config_init(&cfg);cfg.width=cfg.height=1;cfg.log_level=0;
    cfg.eye_set=cfg.target_set=1;cfg.eye[0]=0;cfg.eye[1]=-10;cfg.eye[2]=0;
    cfg.target[0]=cfg.target[1]=cfg.target[2]=0;
    for(int frame=0;frame<2;++frame) {
        compare_render(sys,gpu,&cfg);
        snprintf(path,sizeof(path),"%s-rays.csv",prefix);FILE* f=fopen(path,"r");ASSERT_NOT_NULL(f);
        char header[256];ASSERT_NOT_NULL(fgets(header,sizeof(header),f));
        unsigned x,y,bvh,calls,instructions,primitives,intervals,hierarchy;
        ASSERT_EQ(fscanf(f,"%u,%u,%u,%u,%u,%u,%u,%u",&x,&y,&bvh,&calls,&instructions,&primitives,&intervals,&hierarchy),8);fclose(f);
        ASSERT_EQ(x,0);ASSERT_EQ(y,0);ASSERT_EQ(bvh,1);ASSERT_EQ(calls,1);
        ASSERT_EQ(instructions,2);ASSERT_EQ(primitives,2);ASSERT_EQ(intervals,2);
        snprintf(path,sizeof(path),"%s-cells.csv",prefix);f=fopen(path,"r");ASSERT_NOT_NULL(f);
        ASSERT_NOT_NULL(fgets(header,sizeof(header),f));unsigned index,id,mat,size;
        ASSERT_EQ(fscanf(f,"%u,%u,%u,%u,%u,%u,%u,%u",&index,&id,&mat,&size,&calls,&instructions,&primitives,&intervals),8);fclose(f);
        ASSERT_EQ(id,42);ASSERT_EQ(calls,1);ASSERT_EQ(instructions,2);ASSERT_EQ(primitives,2);ASSERT_EQ(intervals,2);
    }
    const char* suffixes[]={"-rays.csv","-cells.csv","-heatmap.png"};
    for(int i=0;i<3;++i){snprintf(path,sizeof(path),"%s%s",prefix,suffixes[i]);ASSERT_EQ(unlink(path),0);}
    ASSERT_EQ(rmdir(directory),0);alea_compute_destroy(gpu);alea_destroy(sys);render_config_free(&cfg);
}

TEST(metal_resume_feature_regressions) {
    ASSERT_EQ(alea_setenv("ALEA_METAL_RESUME_WORK","7",1),0);
    compare_resume=1;
    test_metal_render_matches_cpu_csg_and_clips();
    test_metal_render_quadrics_cones_and_tori();
    test_metal_render_nested_transforms_and_parent_csg();
    check_rectangular_lattices();
    test_metal_spatial_disconnected_parent();
    test_metal_profile_counts_and_reuse();
    compare_resume=0;
    ASSERT_EQ(alea_unsetenv("ALEA_METAL_RESUME_WORK"),0);
}

typedef struct {size_t calls,unfinished;int cancel;} resume_pause_check;
static int check_resume_pause(void* user,const render_framebuffer_t* fb,size_t done,size_t total) {
    (void)fb;(void)total;resume_pause_check* c=user;++c->calls;if(!done)++c->unfinished;
    return c->cancel && c->calls==12;
}
TEST(metal_resume_inside_expression_cancel_and_reset) {
    alea_system_t* sys=alea_create();int m=alea_add_material(sys,1);alea_node_id_t parts[128];
    for(int i=0;i<128;++i) {
        double angle=6.283185307179586*i/128;
        int sphere=alea_sphere_surface(sys,0,4*cos(angle),0,4*sin(angle),.1);
        parts[i]=alea_halfspace(sys,sphere,-1);
    }
    alea_add_cell(sys,42,balanced_union(sys,parts,128),m,-1,0);
    ASSERT_EQ(alea_setenv("ALEA_METAL_RESUME_WORK","1",1),0);
    alea_compute_t* gpu=alea_compute_create(sys,ALEA_COMPUTE_METAL,ALEA_COMPUTE_FP32);ASSERT_MSG(gpu!=NULL,alea_error());
    ASSERT_EQ(alea_unsetenv("ALEA_METAL_RESUME_WORK"),0);
    render_config_t cfg;test_camera(&cfg);cfg.width=cfg.height=1;cfg.eye[0]=cfg.eye[2]=0;cfg.eye[1]=-10;
    render_camera_t cam;ASSERT_EQ(render_camera_setup(&cam,&cfg,sys),0);
    render_framebuffer_t* fb=render_framebuffer_create(1,1,1);ASSERT_NOT_NULL(fb);fb->cell_id[0]=-777;
    resume_pause_check c={.cancel=1};
    ASSERT_EQ(alea_compute_render_progressive(gpu,&cfg,&cam,fb,check_resume_pause,&c),1);
    ASSERT_EQ(c.calls,12);ASSERT_EQ(c.unfinished,12);ASSERT_EQ(fb->cell_id[0],-777);ASSERT_FALSE(alea_interrupted());
    c=(resume_pause_check){0};
    ASSERT_MSG(alea_compute_render_progressive(gpu,&cfg,&cam,fb,check_resume_pause,&c)==0,alea_error());
    ASSERT_TRUE(c.unfinished>100);ASSERT_EQ(fb->cell_id[0],-1);
    /* A different camera must restart all continuation and AA state. */
    cfg.target[0]=cfg.eye[0]=4;cfg.aa_samples=1;ASSERT_EQ(render_camera_setup(&cam,&cfg,sys),0);
    c=(resume_pause_check){0};ASSERT_EQ(alea_compute_render_progressive(gpu,&cfg,&cam,fb,check_resume_pause,&c),0);
    ASSERT_EQ(fb->cell_id[0],42);ASSERT_NEAR(fb->depth[0],9.9,.002);
    render_framebuffer_free(fb);render_config_free(&cfg);alea_compute_destroy(gpu);alea_destroy(sys);
}

TEST(metal_resume_deep_hierarchy) {
    alea_system_t* sys=alea_create();int m=alea_add_material(sys,7);
    int box=alea_box_surface(sys,1,-4,4,-4,4,-4,4);alea_node_id_t mask=alea_halfspace(sys,box,-1);
    double tr[]={.1,0,0};ASSERT_EQ(alea_add_transform(sys,1,tr,3,0),0);
    for(int u=0;u<18;++u) {
        int cell=alea_add_cell(sys,100+u,mask,ALEA_MATERIAL_VOID,0,u);
        ASSERT_EQ(alea_set_cell_fill(sys,cell,u+1,1),0);
    }
    int sphere=alea_sphere_surface(sys,2,0,0,0,1);alea_add_cell(sys,42,alea_halfspace(sys,sphere,-1),m,-1,18);
    ASSERT_EQ(alea_prepare_query_acceleration(sys),0);
    ASSERT_EQ(alea_setenv("ALEA_METAL_RESUME_WORK","1",1),0);
    alea_compute_t* gpu=alea_compute_create(sys,ALEA_COMPUTE_METAL,ALEA_COMPUTE_FP32);ASSERT_MSG(gpu!=NULL,alea_error());
    ASSERT_EQ(alea_unsetenv("ALEA_METAL_RESUME_WORK"),0);
    render_config_t cfg;test_camera(&cfg);cfg.width=7;cfg.height=5;
    cfg.eye[0]=cfg.target[0]=1.8;cfg.eye[1]=-8;cfg.eye[2]=0;
    compare_resume=1;compare_render(sys,gpu,&cfg);compare_resume=0;
    render_config_free(&cfg);alea_compute_destroy(gpu);alea_destroy(sys);
}
TEST(metal_lazy_parent_faces_and_disjoint_intervals) {
    ASSERT_EQ(alea_setenv("ALEA_METAL_LAZY","all",1),0);
    ASSERT_EQ(alea_setenv("ALEA_METAL_RESUME_WORK","1",1),0);
    for(int mode=0;mode<3;++mode) {
        alea_system_t* sys=alea_create();int mat=alea_add_material(sys,7);
        int a=alea_box_surface(sys,0,-3,3,mode==2?3:0,mode==1?.2:10,-3,3);
        alea_node_id_t parent=alea_halfspace(sys,a,-1);
        if(mode==1) {
            int b=alea_box_surface(sys,0,-3,3,1,1.5,-3,3);
            parent=alea_union(sys,parent,alea_halfspace(sys,b,-1));
        }
        int c=alea_add_cell(sys,10,parent,0,0,0);ASSERT_EQ(alea_set_cell_fill(sys,c,1,0),0);
        int b=alea_box_surface(sys,0,-4,4,-1,8,-4,4);
        c=alea_add_cell(sys,11,alea_halfspace(sys,b,-1),0,0,1);ASSERT_EQ(alea_set_cell_fill(sys,c,2,0),0);
        int child=alea_box_surface(sys,0,-2,2,-2,2,-2,2);
        alea_add_cell(sys,42,alea_halfspace(sys,child,-1),mat,-1,2);
        ASSERT_EQ(alea_prepare_query_acceleration(sys),0);
        alea_compute_t* gpu=alea_compute_create(sys,ALEA_COMPUTE_METAL,ALEA_COMPUTE_FP32);ASSERT_MSG(gpu!=NULL,alea_error());
        render_config_t cfg;test_camera(&cfg);cfg.width=cfg.height=1;
        cfg.eye[0]=cfg.eye[2]=0;cfg.eye[1]=-10;
        render_camera_t camera;ASSERT_EQ(render_camera_setup(&camera,&cfg,sys),0);
        render_framebuffer_t* fb=render_framebuffer_create(1,1,1);ASSERT_NOT_NULL(fb);
        for(int progressive=0;progressive<2;++progressive) {
            int callbacks=0;
            ASSERT_EQ(progressive?alea_compute_render_progressive(gpu,&cfg,&camera,fb,count_progress,&callbacks):alea_compute_render(gpu,&cfg,&camera,fb),0);
            ASSERT_EQ(fb->cell_id[0],mode==2?-1:42);
            if(mode!=2){ASSERT_NEAR(fb->depth[0],10,1e-5);ASSERT_NEAR(fb->normal[1],-1,1e-5);}
        }
        render_framebuffer_free(fb);render_config_free(&cfg);alea_compute_destroy(gpu);alea_destroy(sys);
    }
    ASSERT_EQ(alea_unsetenv("ALEA_METAL_RESUME_WORK"),0);
    ASSERT_EQ(alea_unsetenv("ALEA_METAL_LAZY"),0);
}
TEST(metal_lazy_parent_feature_regressions) {
    ASSERT_EQ(alea_setenv("ALEA_METAL_LAZY","all",1),0);
    ASSERT_EQ(alea_setenv("ALEA_METAL_RESUME_WORK","7",1),0);
    int previous=compare_resume;compare_resume=1;
    test_metal_render_nested_transforms_and_parent_csg();
    test_metal_spatial_shared_content_and_original_programs();
    test_metal_spatial_thin_parent_boundary();
    test_metal_spatial_hint_is_not_a_clip();
    test_metal_spatial_distant_ray_fallback();
    test_metal_spatial_proved_empty_regions();
    test_metal_spatial_disconnected_parent();
    test_metal_resume_deep_hierarchy();
    check_rectangular_lattices();
    compare_resume=previous;
    ASSERT_EQ(alea_unsetenv("ALEA_METAL_RESUME_WORK"),0);
    ASSERT_EQ(alea_unsetenv("ALEA_METAL_LAZY"),0);
}

TEST(metal_lazy_parent_reuses_rejected_range) {
    char directory[]="/tmp/alea-guard-XXXXXX";ASSERT_NOT_NULL(mkdtemp(directory));
    char prefix[256],path[300];snprintf(prefix,sizeof(prefix),"%s/profile",directory);
    ASSERT_EQ(alea_setenv("ALEA_METAL_PROFILE",prefix,1),0);
    ASSERT_EQ(alea_setenv("ALEA_METAL_LAZY","all",1),0);
    ASSERT_EQ(alea_setenv("ALEA_METAL_RESUME_WORK","1",1),0);
    alea_system_t* sys=alea_create();int mat=alea_add_material(sys,7);
    int outer=alea_box_surface(sys,0,-5,5,-5,5,-5,5),hole=alea_sphere_surface(sys,0,0,0,0,3);
    int parent=alea_add_cell(sys,10,alea_difference(sys,alea_halfspace(sys,outer,-1),alea_halfspace(sys,hole,-1)),0,0,0);
    ASSERT_EQ(alea_set_cell_fill(sys,parent,1,0),0);
    int a=alea_box_surface(sys,0,-1,1,-2,2,-1,1),b=alea_box_surface(sys,0,-1,1,-1,4,-1,1);
    alea_add_cell(sys,42,alea_halfspace(sys,a,-1),mat,-1,1);
    alea_add_cell(sys,43,alea_halfspace(sys,b,-1),mat,-1,1);
    ASSERT_EQ(alea_prepare_query_acceleration(sys),0);
    alea_compute_t* gpu=alea_compute_create(sys,ALEA_COMPUTE_METAL,ALEA_COMPUTE_FP32);ASSERT_MSG(gpu!=NULL,alea_error());
    render_config_t cfg;test_camera(&cfg);cfg.width=cfg.height=1;cfg.eye[0]=cfg.eye[2]=0;cfg.eye[1]=-10;
    render_camera_t camera;ASSERT_EQ(render_camera_setup(&camera,&cfg,sys),0);
    render_framebuffer_t* fb=render_framebuffer_create(1,1,1);ASSERT_NOT_NULL(fb);
    for(int progressive=0;progressive<2;++progressive) {
        int callbacks=0;ASSERT_EQ(progressive?alea_compute_render_progressive(gpu,&cfg,&camera,fb,count_progress,&callbacks):alea_compute_render(gpu,&cfg,&camera,fb),0);
        ASSERT_EQ(fb->cell_id[0],43);ASSERT_NEAR(fb->depth[0],13,1e-5);ASSERT_NEAR(fb->normal[1],-1,1e-5);
        snprintf(path,sizeof(path),"%s-rays.csv",prefix);FILE* f=fopen(path,"r");ASSERT_NOT_NULL(f);
        char header[256];ASSERT_NOT_NULL(fgets(header,sizeof(header),f));unsigned values[11];
        ASSERT_EQ(fscanf(f,"%u,%u,%u,%u,%u,%u,%u,%u,%u,%u,%u",values,values+1,values+2,values+3,values+4,values+5,values+6,values+7,values+8,values+9,values+10),11);fclose(f);
        ASSERT_TRUE(values[8]>0);ASSERT_TRUE(values[9]>0);ASSERT_TRUE(values[10]>0);
    }
    const char* suffixes[]={"-rays.csv","-cells.csv","-heatmap.png"};
    for(int i=0;i<3;++i){snprintf(path,sizeof(path),"%s%s",prefix,suffixes[i]);ASSERT_EQ(unlink(path),0);}
    ASSERT_EQ(rmdir(directory),0);
    render_framebuffer_free(fb);render_config_free(&cfg);alea_compute_destroy(gpu);alea_destroy(sys);
    ASSERT_EQ(alea_unsetenv("ALEA_METAL_PROFILE"),0);ASSERT_EQ(alea_unsetenv("ALEA_METAL_LAZY"),0);ASSERT_EQ(alea_unsetenv("ALEA_METAL_RESUME_WORK"),0);
}

TEST(metal_lazy_parent_cache_backtracking) {
    ASSERT_EQ(alea_setenv("ALEA_METAL_LAZY","all",1),0);
    alea_system_t* sys=alea_create();int mat=alea_add_material(sys,7);
    int outer=alea_box_surface(sys,0,-8,8,-8,8,-8,8),hole=alea_sphere_surface(sys,0,0,0,0,3);
    int parent=alea_add_cell(sys,10,alea_difference(sys,alea_halfspace(sys,outer,-1),alea_halfspace(sys,hole,-1)),0,0,0);
    ASSERT_EQ(alea_set_cell_fill(sys,parent,1,0),0);
    int a=alea_box_surface(sys,0,-1,1,5,6,-1,1),b=alea_box_surface(sys,0,-1,1,-1,0,-1,1);
    alea_node_id_t part=alea_halfspace(sys,a,-1);
    alea_bbox_t loose={-1,1,-5,6,-1,1};alea_node_bbox_set(&sys->nodes.data[part].bbox,&loose);
    alea_add_cell(sys,42,part,mat,-1,1);alea_add_cell(sys,43,alea_halfspace(sys,b,-1),mat,-1,1);
    ASSERT_EQ(alea_prepare_query_acceleration(sys),0);
    alea_compute_t* gpu=alea_compute_create(sys,ALEA_COMPUTE_METAL,ALEA_COMPUTE_FP32);ASSERT_MSG(gpu!=NULL,alea_error());
    render_config_t cfg;test_camera(&cfg);cfg.width=cfg.height=1;cfg.eye[0]=cfg.eye[2]=0;cfg.eye[1]=-10;
    render_camera_t camera;ASSERT_EQ(render_camera_setup(&camera,&cfg,sys),0);
    render_framebuffer_t* fb=render_framebuffer_create(1,1,1);ASSERT_NOT_NULL(fb);
    for(int progressive=0;progressive<2;++progressive) {
        int calls=0;ASSERT_EQ(progressive?alea_compute_render_progressive(gpu,&cfg,&camera,fb,count_progress,&calls):alea_compute_render(gpu,&cfg,&camera,fb),0);
        ASSERT_EQ(fb->cell_id[0],42);ASSERT_NEAR(fb->depth[0],15,1e-5);
    }
    render_framebuffer_free(fb);render_config_free(&cfg);alea_compute_destroy(gpu);alea_destroy(sys);
    ASSERT_EQ(alea_unsetenv("ALEA_METAL_LAZY"),0);
}
TEST(metal_lazy_parent_preserves_deeper_coincident_face) {
    alea_system_t* sys=alea_create();int mat=alea_add_material(sys,7);
    int box=alea_box_surface(sys,0,-4,4,-4,4,-4,4),plane=alea_plane_surface(sys,0,0,1,0,0);
    alea_node_id_t mask=alea_intersection(sys,alea_halfspace(sys,box,-1),alea_halfspace(sys,plane,1));
    for(int i=0;i<32;++i) {
        int hole=alea_sphere_surface(sys,0,3,3,3+i*.01,.01);
        mask=alea_difference(sys,mask,alea_halfspace(sys,hole,-1));
    }
    int parent=alea_add_cell(sys,10,mask,0,0,0);ASSERT_EQ(alea_set_cell_fill(sys,parent,1,0),0);
    int tilted=alea_plane_surface(sys,0,1,1,0,0);
    mask=alea_intersection(sys,alea_halfspace(sys,box,-1),alea_halfspace(sys,tilted,1));
    parent=alea_add_cell(sys,11,mask,0,0,1);ASSERT_EQ(alea_set_cell_fill(sys,parent,2,0),0);
    alea_add_cell(sys,42,alea_halfspace(sys,box,-1),mat,-1,2);
    ASSERT_EQ(alea_prepare_query_acceleration(sys),0);
    render_config_t cfg;test_camera(&cfg);cfg.width=cfg.height=1;cfg.eye[0]=cfg.eye[2]=0;cfg.eye[1]=-10;
    render_camera_t camera;ASSERT_EQ(render_camera_setup(&camera,&cfg,sys),0);
    render_framebuffer_t* fb=render_framebuffer_create(1,1,1);ASSERT_NOT_NULL(fb);float normal[3]={0};
    for(int mode=0;mode<2;++mode) {
        ASSERT_EQ(alea_setenv("ALEA_METAL_LAZY",mode?"1":"0",1),0);
        alea_compute_t* gpu=alea_compute_create(sys,ALEA_COMPUTE_METAL,ALEA_COMPUTE_FP32);ASSERT_MSG(gpu!=NULL,alea_error());
        for(int progressive=0;progressive<2;++progressive) {
            int calls=0;ASSERT_EQ(progressive?alea_compute_render_progressive(gpu,&cfg,&camera,fb,count_progress,&calls):alea_compute_render(gpu,&cfg,&camera,fb),0);
            ASSERT_EQ(fb->cell_id[0],42);ASSERT_NEAR(fb->depth[0],10,1e-5);
            if(!mode && !progressive)memcpy(normal,fb->normal,sizeof(normal));
            ASSERT_EQ(memcmp(normal,fb->normal,sizeof(normal)),0);
        }
        alea_compute_destroy(gpu);
    }
    render_framebuffer_free(fb);render_config_free(&cfg);alea_destroy(sys);ASSERT_EQ(alea_unsetenv("ALEA_METAL_LAZY"),0);
}

TEST(metal_weighted_bvh_structure) {
    ASSERT_EQ(alea_setenv("ALEA_METAL_BVH","weighted",1),0);
    alea_system_t* sys=alea_create();int mat=alea_add_material(sys,7);
    for(int i=0;i<16;++i) {
        int sphere=alea_sphere_surface(sys,0,((i&1)?5:-5)+(i>>3)*2,(i&2)?5:-5,(i&4)?5:-5,.5);
        alea_add_cell(sys,100+i,alea_halfspace(sys,sphere,-1),mat,-1,0);
    }
    ASSERT_EQ(alea_prepare_query_acceleration(sys),0);
    compute_scene packed;ASSERT_EQ(alea_compute_pack(sys,&packed),0);
    for(size_t i=0;i<packed.render.bvh_count;++i) {
        compute_render_bvh* b=&packed.render.bvh[i];ASSERT_TRUE(b->escape>i);ASSERT_TRUE(b->escape<=packed.render.bvh_count);
    }
    alea_compute_scene_free(&packed);
    alea_compute_t* gpu=alea_compute_create(sys,ALEA_COMPUTE_METAL,ALEA_COMPUTE_FP32);ASSERT_MSG(gpu!=NULL,alea_error());
    render_config_t cfg;test_camera(&cfg);cfg.eye[1]=-25;cfg.width=41;cfg.height=29;
    int previous=compare_resume;compare_resume=1;compare_render(sys,gpu,&cfg);compare_resume=previous;
    render_config_free(&cfg);alea_compute_destroy(gpu);alea_destroy(sys);
    ASSERT_EQ(alea_unsetenv("ALEA_METAL_BVH"),0);
}
TEST(metal_spatial_preserves_lattice_envelope) {
    alea_system_t* sys=alea_create(),*reference=alea_create();int mat=alea_add_material(sys,7),rm=alea_add_material(reference,7);
    int unit=alea_box_surface(sys,0,-1,1,-1,1,-1,1);
    int lat=alea_add_cell(sys,10,alea_halfspace(sys,unit,-1),0,0,0);
    alea_cell_entry_t* cell=&sys->cells.data[lat];cell->lat_type=1;
    cell->lat_fill_dims[1]=7;cell->lat_fill_count=8;cell->lat_fill=malloc(8*sizeof(int));ASSERT_NOT_NULL(cell->lat_fill);
    for(int a=0;a<3;++a){cell->lat_pitch[a]=2;cell->lat_lower_left[a]=-1;}
    for(int i=0;i<8;++i)cell->lat_fill[i]=1;
    int sphere=alea_sphere_surface(sys,0,0,0,0,.4);alea_add_cell(sys,42,alea_halfspace(sys,sphere,-1),mat,-1,1);
    alea_node_id_t parts[8];
    for(int i=0;i<8;++i){sphere=alea_sphere_surface(reference,0,2*i,0,0,.4);parts[i]=alea_halfspace(reference,sphere,-1);}
    alea_add_cell(reference,42,balanced_union(reference,parts,8),rm,-1,0);
    for(int i=0;i<32;++i)for(int which=0;which<2;++which) {
        alea_system_t* model=which?reference:sys;
        sphere=alea_sphere_surface(model,0,i<16?-6:20,(i%4-1.5),(i/4%4-1.5),.08);
        alea_add_cell(model,100+i,alea_halfspace(model,sphere,-1),which?rm:mat,-1,0);
    }
    ASSERT_EQ(alea_prepare_query_acceleration(sys),0);ASSERT_EQ(alea_prepare_query_acceleration(reference),0);
    compute_scene packed;ASSERT_EQ(alea_compute_pack(sys,&packed),0);
    alea_compute_scene_free(&packed);
    alea_compute_t* gpu=alea_compute_create(sys,ALEA_COMPUTE_METAL,ALEA_COMPUTE_FP32);ASSERT_MSG(gpu!=NULL,alea_error());
    render_config_t cfg;test_camera(&cfg);cfg.width=17;cfg.height=13;
    cfg.eye[0]=cfg.target[0]=14;cfg.eye[1]=-8;cfg.eye[2]=0;
    int previous=compare_resume;compare_resume=1;compare_render(reference,gpu,&cfg);compare_resume=previous;
    render_config_free(&cfg);alea_compute_destroy(gpu);alea_destroy(sys);alea_destroy(reference);
}

TEST(metal_acceleration_feature_regressions) {
    ASSERT_EQ(alea_setenv("ALEA_METAL_LAZY","all",1),0);
    ASSERT_EQ(alea_setenv("ALEA_METAL_REFILL","1",1),0);
    ASSERT_EQ(alea_setenv("ALEA_METAL_RESUME_WORK","7",1),0);
    int previous=compare_resume;compare_resume=1;
    for(int mode=0;mode<2;++mode) {
        ASSERT_EQ(alea_setenv("ALEA_METAL_BVH",mode?"weighted":"sah",1),0);
        test_metal_render_matches_cpu_csg_and_clips();
        test_metal_render_quadrics_cones_and_tori();
        test_metal_render_nested_transforms_and_parent_csg();
        test_metal_spatial_thin_parent_boundary();
        test_metal_spatial_proved_empty_regions();
        test_metal_spatial_hint_is_not_a_clip();
        test_metal_spatial_distant_ray_fallback();
        check_rectangular_lattices();
    }
    test_metal_progressive_tiles_cancel_and_reuse();
    test_metal_resume_inside_expression_cancel_and_reset();
    compare_resume=previous;
    ASSERT_EQ(alea_unsetenv("ALEA_METAL_BVH"),0);
    ASSERT_EQ(alea_unsetenv("ALEA_METAL_LAZY"),0);
    ASSERT_EQ(alea_unsetenv("ALEA_METAL_REFILL"),0);
    ASSERT_EQ(alea_unsetenv("ALEA_METAL_RESUME_WORK"),0);
}
TEST(metal_acceleration_invalid_configuration) {
    alea_system_t* sys=sphere_scene();
    const char* flags[]={"ALEA_METAL_LAZY","ALEA_METAL_BVH","ALEA_METAL_REFILL"};
    for(int i=0;i<3;++i) {
        ASSERT_EQ(alea_setenv(flags[i],"invalid",1),0);
        ASSERT_NULL(alea_compute_create(sys,ALEA_COMPUTE_METAL,ALEA_COMPUTE_FP32));
        ASSERT_EQ(alea_error_code(),ALEA_ERR_INVALID_ARG);
        ASSERT_EQ(alea_unsetenv(flags[i]),0);
    }
    alea_destroy(sys);
}

TEST(metal_resume_errors_and_validation) {
    alea_system_t* sys=alea_create();int m=alea_add_material(sys,1);
    int sphere=alea_sphere_surface(sys,1,0,0,0,1e20);alea_add_cell(sys,42,alea_halfspace(sys,sphere,-1),m,-1,0);
    const char* bad[]={"0","-1","65537","1x",""};
    for(size_t i=0;i<sizeof(bad)/sizeof(bad[0]);++i) {
        ASSERT_EQ(alea_setenv("ALEA_METAL_RESUME_WORK",bad[i],1),0);
        ASSERT_NULL(alea_compute_create(sys,ALEA_COMPUTE_METAL,ALEA_COMPUTE_FP32));ASSERT_EQ(alea_error_code(),ALEA_ERR_INVALID_ARG);
    }
    ASSERT_EQ(alea_setenv("ALEA_METAL_RESUME_WORK","1",1),0);
    alea_compute_t* gpu=alea_compute_create(sys,ALEA_COMPUTE_METAL,ALEA_COMPUTE_FP32);ASSERT_MSG(gpu!=NULL,alea_error());
    ASSERT_EQ(alea_unsetenv("ALEA_METAL_RESUME_WORK"),0);
    render_config_t cfg;test_camera(&cfg);cfg.width=cfg.height=1;
    cfg.eye[0]=cfg.eye[2]=0;cfg.eye[1]=-10;render_camera_t cam;ASSERT_EQ(render_camera_setup(&cam,&cfg,sys),0);
    render_framebuffer_t* fb=render_framebuffer_create(1,1,1);ASSERT_NOT_NULL(fb);int calls=0;
    ASSERT_EQ(alea_compute_render_progressive(gpu,&cfg,&cam,fb,count_progress,&calls),-1);ASSERT_EQ(alea_error_code(),ALEA_ERR_OVERFLOW);
    cfg.num_clips=1;cfg.clips[0]=(render_clip_plane_t){{0,-1,0},-20};
    ASSERT_EQ(alea_compute_render_progressive(gpu,&cfg,&cam,fb,count_progress,&calls),0);ASSERT_EQ(fb->cell_id[0],-1);
    render_framebuffer_free(fb);render_config_free(&cfg);alea_compute_destroy(gpu);alea_destroy(sys);
}

#endif

TEST_MAIN()
