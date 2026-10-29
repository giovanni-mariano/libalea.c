// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
// SPDX-License-Identifier: MPL-2.0
#define _POSIX_C_SOURCE 200809L
#include "alea.h"
#include "alea_compute.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

/* Public synthetic workload: 64 translated occurrences of one sphere/void
 * universe. Point ownership and visible spheres have independent equations.
 * Timings include public API validation, dispatch, synchronization and copies;
 * oracle checks, scene creation and display are outside the timed interval. */
static double seconds(void) {
    struct timespec t; clock_gettime(CLOCK_MONOTONIC, &t);
    return t.tv_sec + t.tv_nsec * 1e-9;
}

static alea_system_t* scene(void) {
    alea_system_t* s = alea_create();
    if (!s) return NULL;
    int mat = alea_add_material(s, 7);
    int sphere = alea_sphere_surface(s, 1, 0, 0, 0, .75);
    if (mat < 0 || sphere < 0 ||
        alea_add_cell(s, 42, alea_halfspace(s, sphere, -1), mat, -1, 1) < 0 ||
        alea_add_cell(s, 43, alea_halfspace(s, sphere, 1), ALEA_MATERIAL_VOID, 0, 1) < 0)
        goto fail;
    for (int y = 0; y < 8; ++y) for (int x = 0; x < 8; ++x) {
        int id = y * 8 + x;
        double tr[] = {2*x - 7, 2*y - 7, 0};
        int box = alea_box_surface(s, 100 + id, tr[0]-1, tr[0]+1,
                                   tr[1]-1, tr[1]+1, -1, 1);
        if (box < 0) goto fail;
        int parent = alea_add_cell(s, 1000 + id, alea_halfspace(s, box, -1),
                                   ALEA_MATERIAL_VOID, 0, 0);
        if (parent < 0 || alea_add_transform(s, id+1, tr, 3, 0) ||
            alea_set_fill(s, parent, 1, id+1)) goto fail;
    }
    return s;
fail:
    alea_destroy(s); return NULL;
}

static int check_points(const double* xyz, const alea_compute_result_t* r, size_t n) {
    for (size_t i = 0; i < n; ++i) {
        double x = xyz[3*i], y = xyz[3*i+1];
        double cx = 2*floor((x+8)/2)-7, cy = 2*floor((y+8)/2)-7;
        int solid = (x-cx)*(x-cx)+(y-cy)*(y-cy) <= .75*.75;
        if (r[i].owner_count != 1 || r[i].cell_id != (solid ? 42 : 43) ||
            r[i].material_id != (solid ? 7 : 0)) return -1;
    }
    return 0;
}

static int check_render(const render_config_t* cfg, const render_camera_t* cam,
                         const render_framebuffer_t* fb) {
    for (int y = 0; y < cfg->height; ++y) for (int x = 0; x < cfg->width; ++x) {
        double o[3], d[3];
        render_camera_ray(cam, cfg->width, cfg->height, x+.5, y+.5,
                           &o[0], &o[1], &o[2], &d[0], &d[1], &d[2]);
        double begin = 0, end = INFINITY, best = INFINITY, normal[3] = {0};
        if (cfg->num_clips) {
            if (d[1] > 0) begin = fmax(begin, -o[1]/d[1]);
            else if (d[1] < 0) end = -o[1]/d[1];
            else if (o[1] < 0) end = -1;
        }
        for (int j = 0; j < 8; ++j) for (int i = 0; i < 8; ++i) {
            double p[] = {o[0]-(2*i-7), o[1]-(2*j-7), o[2]};
            double a = d[0]*d[0]+d[1]*d[1]+d[2]*d[2];
            double b = p[0]*d[0]+p[1]*d[1]+p[2]*d[2];
            double c = p[0]*p[0]+p[1]*p[1]+p[2]*p[2]-.75*.75;
            double disc = b*b-a*c;
            if (disc < 0) continue;
            double lo = (-b-sqrt(disc))/a, hi = (-b+sqrt(disc))/a;
            double hit = fmax(begin, lo);
            if (hit >= fmin(end, hi) || hit >= best) continue;
            best = hit;
            for (int k = 0; k < 3; ++k) normal[k] = (p[k]+hit*d[k])/.75;
            if (hit == begin && begin > 0) { normal[0]=normal[2]=0; normal[1]=-1; }
        }
        size_t pixel = (size_t)y*cfg->width+x;
        int material = isfinite(best) ? 7 : 0;
        if (fb->material_id[pixel] != material ||
            (material && (fb->cell_id[pixel] != 42 || fabs(fb->depth[pixel]-best) > .003)))
            return -1;
        if (material) for (int k = 0; k < 3; ++k)
            if (fabs(fb->normal[3*pixel+k]-normal[k]) > .003) return -1;
    }
    return 0;
}

enum { PASSES = 6 }; /* One initial call and five warm calls. */
static int time_order(const void* a,const void* b) {
    double x=*(const double*)a,y=*(const double*)b;return (x>y)-(x<y);
}
static void report(const char* task, const char* backend, const double time[PASSES]) {
    double warm[PASSES-1];
    for(int i=1;i<PASSES;++i)warm[i-1]=time[i];
    qsort(warm,PASSES-1,sizeof(*warm),time_order);
    double median=warm[(PASSES-1)/2];
    printf("task=%s backend=%s first_ms=%.3f warm_median_ms=%.3f warm_min_ms=%.3f warm_max_ms=%.3f\n",
           task, backend, time[0]*1000, median*1000,
           warm[0]*1000, warm[PASSES-2]*1000);
    fflush(stdout);
}

int main(void) {
    int rc = 1;
    const size_t n = 512*512;
    alea_system_t* sys = scene();
    double* xyz = malloc(n*3*sizeof(*xyz));
    alea_compute_result_t* results = malloc(n*sizeof(*results));
    alea_compute_selection_t* selected=malloc(n*sizeof(*selected));
    alea_compute_t* cpu = NULL; alea_compute_t* metal = NULL;
    render_framebuffer_t* fb = NULL;
    render_config_t cfg; render_config_init(&cfg);
    if (!sys || !xyz || !results || !selected) goto done;
    for (size_t i = 0; i < n; ++i) {
        xyz[3*i] = -8 + ((i%512)+.5)/32;
        xyz[3*i+1] = -8 + ((i/512)+.5)/32;
        xyz[3*i+2] = 0;
    }
    double start=seconds();
    cpu=alea_compute_create(sys,ALEA_COMPUTE_CPU,ALEA_COMPUTE_FP64);
    if (!cpu) goto done;
    printf("backend=cpu setup_ms=%.3f\n",(seconds()-start)*1000);
    start=seconds();
    metal=alea_compute_create(sys,ALEA_COMPUTE_METAL,ALEA_COMPUTE_FP32);
    if (!metal) goto done;
    printf("backend=metal setup_ms=%.3f\n",(seconds()-start)*1000);
    printf("Point reference enumerates each universe; it is not the core spatial lookup.\n");
    double timings[2][PASSES];
    for (int pass=0;pass<PASSES;++pass) for (int order=0;order<2;++order) {
        int backend=(pass+order)%2;
        start=seconds();
        if (alea_compute_classify(backend?metal:cpu,xyz,n,results)) goto done;
        timings[backend][pass]=seconds()-start;
        if (check_points(xyz,results,n)) {fprintf(stderr,"Point oracle mismatch\n");goto done;}
    }
    report("slice512","cpu",timings[0]);report("slice512","metal",timings[1]);
    /* Matched ordered selection: prepared production core lookup, one CPU
     * thread, versus one native Metal batch. Setup stays outside the timer. */
    if(alea_prepare_query_acceleration(sys))goto done;
    for(int pass=0;pass<PASSES;++pass)for(int order=0;order<2;++order) {
        int backend=(pass+order)%2;start=seconds();
        if(backend) {
            if(alea_compute_select(metal,xyz,n,selected))goto done;
        } else for(size_t i=0;i<n;++i) {
            if(alea_find_cell_at(sys,xyz[3*i],xyz[3*i+1],xyz[3*i+2],
                &selected[i].cell_id,&selected[i].material_id))goto done;
            selected[i].depth=1;selected[i].status=ALEA_COMPUTE_FOUND;
        }
        timings[backend][pass]=seconds()-start;
        for(size_t i=0;i<n;++i)if(selected[i].status!=ALEA_COMPUTE_FOUND || selected[i].depth!=1 ||
            selected[i].cell_id!=results[i].cell_id || selected[i].material_id!=results[i].material_id) {
            fprintf(stderr,"Selected ownership oracle mismatch\n");goto done;
        }
    }
    report("selected512","cpu-core-serial",timings[0]);report("selected512","metal",timings[1]);
    cfg.width=512;cfg.height=384;cfg.aa_samples=1;cfg.log_level=0;
    cfg.eye_set=cfg.target_set=1;
    cfg.eye[0]=16;cfg.eye[1]=-20;cfg.eye[2]=18;
    cfg.target[0]=cfg.target[1]=cfg.target[2]=0;
    render_camera_t cam;
    if (render_camera_setup(&cam,&cfg,sys)) goto done;
    fb=render_framebuffer_create(cfg.width,cfg.height,1);
    if (!fb) goto done;
    for (int view=0;view<3;++view) {
        if(view==2) {
            cfg.eye[0]=0;cfg.eye[1]=-3;cfg.eye[2]=3;
            cfg.target[0]=1;cfg.target[1]=-1;cfg.target[2]=0;
            if(render_camera_setup(&cam,&cfg,sys))goto done;
        }
        cfg.num_clips=view==1;cfg.clips[0]=(render_clip_plane_t){{0,1,0},0};
        double time[PASSES];
        for (int pass=0;pass<PASSES;++pass) {
            start=seconds();
            if (alea_compute_render(metal,&cfg,&cam,fb)) goto done;
            time[pass]=seconds()-start;
            if (check_render(&cfg,&cam,fb)) {fprintf(stderr,"Ray oracle mismatch\n");goto done;}
        }
        report(view==2?"closeup512x384":view==1?"cutaway512x384":"exterior512x384","metal",time);
    }
    puts("All sampled ownership and analytic ray checks passed.");rc=0;
done:
    if (rc) fprintf(stderr,"%s\n",alea_error());
    render_framebuffer_free(fb);render_config_free(&cfg);
    alea_compute_destroy(cpu);alea_compute_destroy(metal);alea_destroy(sys);
    free(xyz);free(results);free(selected);return rc;
}
