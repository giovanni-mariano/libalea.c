// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
// SPDX-License-Identifier: MPL-2.0
#include "alea_test.h"
#include "alea.h"
#include "compute/compute_internal.h"
/* Deterministic interruption at internal checkpoints, without sleeps or
 * production failure hooks. These implementations replace archive symbols. */
static size_t checkpoints,interrupt_at;
static int builder_interrupted(void) {++checkpoints;return interrupt_at && checkpoints>=interrupt_at;}
#define alea_interrupted builder_interrupted
#include "../../src/compute/compute_bvh.c"
#include "../../src/compute/compute_point.c"
#undef alea_interrupted

static alea_system_t* grid_scene(void) {
    alea_system_t* sys=alea_create();int mat=alea_add_material(sys,1);
    for(int i=0;i<2048;++i) {
        int sphere=alea_sphere_surface(sys,0,3*(i%16),3*(i/16%16),3*(i/256),.5);
        alea_add_cell(sys,i+1,alea_halfspace(sys,sphere,-1),mat,-1,0);
    }
    return sys;
}
TEST(bvh_interruption_preserves_scene_and_allows_retry) {
    alea_system_t* sys=grid_scene();
    ASSERT_EQ(alea_setenv("ALEA_METAL_BVH","weighted",1),0);
    for(unsigned stage=0;stage<3;++stage) {
        compute_scene s;ASSERT_EQ(alea_compute_pack_geometry(sys,&s,0),0);
        size_t bytes=s.render.bvh_count*sizeof(*s.render.bvh);
        compute_render_bvh* original=s.render.bvh,*copy=malloc(bytes);ASSERT_NOT_NULL(copy);memcpy(copy,original,bytes);
        checkpoints=0;interrupt_at=stage==0?2:stage==1?8:20;
        int result=alea_compute_bvh_rebuild(&s);interrupt_at=0;
        ASSERT_EQ(result,-1);ASSERT_EQ(alea_error_code(),ALEA_ERR_INTERRUPTED);
        ASSERT_TRUE(s.render.bvh==original);ASSERT_EQ(memcmp(copy,original,bytes),0);
        ASSERT_EQ(alea_compute_bvh_rebuild(&s),0);
        ASSERT_EQ(s.render.bvh[s.render.root].escape,s.render.bvh_count);
        free(copy);alea_compute_scene_free(&s);
    }
    ASSERT_EQ(alea_unsetenv("ALEA_METAL_BVH"),0);alea_destroy(sys);
}

TEST(point_program_certification_preserves_errors_and_retries) {
    compute_primitive primitive={.type=COMPUTE_SPHERE,.p={0,0,0,1e20}};
    compute_node node={COMPUTE_LEAF,0};
    compute_cell cell={.root=0,.count=1,.spatial=0,.lattice=COMPUTE_NONE};
    compute_spatial_program region={.bounds={-1,1,-1,1,-1,1},.child=COMPUTE_NONE};
    metal_point_cell link={.spatial=COMPUTE_NONE};
    compute_scene s={.cells=&cell,.cell_count=1,.primitives=&primitive,.nodes=&node,
        .point_cells=&link,.spatial=&region,.spatial_count=1};
    /* The original sphere overflows even in this small region. */
    ASSERT_EQ(alea_compute_point_programs_pack(&s),0);ASSERT_EQ(link.spatial,COMPUTE_NONE);
    primitive.p[3]=1;
    checkpoints=0;interrupt_at=1;
    ASSERT_EQ(alea_compute_point_programs_pack(&s),-1);ASSERT_EQ(link.spatial,COMPUTE_NONE);
    interrupt_at=0;
    ASSERT_EQ(alea_compute_point_programs_pack(&s),0);ASSERT_EQ(link.spatial,0);
    ASSERT_EQ(alea_setenv("ALEA_METAL_POINT_SPATIAL","invalid",1),0);
    ASSERT_EQ(alea_compute_point_programs_pack(&s),-1);ASSERT_EQ(alea_error_code(),ALEA_ERR_INVALID_ARG);
    ASSERT_EQ(alea_unsetenv("ALEA_METAL_POINT_SPATIAL"),0);
}
TEST(point_candidates_keep_unknown_boundaries_and_retry) {
    alea_system_t* sys=grid_scene();compute_scene s;
    ASSERT_EQ(alea_compute_pack_geometry(sys,&s,0),0);alea_destroy(sys);
    int8_t* signs=malloc(s.primitive_count);uint32_t* next=malloc(s.cell_count*sizeof(*next));
    ASSERT_NOT_NULL(signs);ASSERT_NOT_NULL(next);
    float bounds[]={-.1,.1,-.1,.1,-.1,.1};uint32_t root;compute_candidate_stats stats={0};
    ASSERT_EQ(alea_compute_point_candidates(&s,bounds,signs,next,&root,&stats),0);
    ASSERT_EQ(stats.cells,2048);ASSERT_EQ(stats.rejected,2047);ASSERT_EQ(root,0);ASSERT_EQ(next[0],COMPUTE_NONE);
    bounds[0]=bounds[1]=.5f;stats=(compute_candidate_stats){0};
    ASSERT_EQ(alea_compute_point_candidates(&s,bounds,signs,next,&root,&stats),0);ASSERT_EQ(root,0);
    bounds[0]=bounds[1]=1e30f;stats=(compute_candidate_stats){0};
    ASSERT_EQ(alea_compute_point_candidates(&s,bounds,signs,next,&root,&stats),0);ASSERT_EQ(stats.rejected,0);
    checkpoints=0;interrupt_at=2;
    ASSERT_EQ(alea_compute_point_candidates(&s,bounds,signs,next,&root,&stats),-1);
    ASSERT_EQ(alea_error_code(),ALEA_ERR_INTERRUPTED);interrupt_at=0;
    bounds[0]=-.1f;bounds[1]=.1f;stats=(compute_candidate_stats){0};
    ASSERT_EQ(alea_compute_point_candidates(&s,bounds,signs,next,&root,&stats),0);ASSERT_EQ(stats.rejected,2047);
    free(signs);free(next);alea_compute_scene_free(&s);
}
TEST(point_candidates_preserve_nonfinite_operands_and_budget_fallback) {
    compute_primitive primitives[2]={
        {.type=COMPUTE_SPHERE,.p={0,0,0,1e20}},
        {.type=COMPUTE_BOX,.p={10,11,10,11,10,11}}};
    compute_node nodes[]={{COMPUTE_LEAF,0},{COMPUTE_SKIP_INTERSECTION,2},
                          {COMPUTE_LEAF,2},{COMPUTE_INTERSECTION,0}};
    compute_cell cell={.root=0,.count=4,.lattice=COMPUTE_NONE};
    compute_scene scene={.nodes=nodes,.primitives=primitives,.primitive_count=2,.cells=&cell,.cell_count=1};
    const float bounds[]={-1,1,-1,1,-1,1};int8_t signs[2];uint32_t next,root;
    compute_candidate_stats stats={0};
    ASSERT_EQ(alea_compute_point_candidates(&scene,bounds,signs,&next,&root,&stats),0);
    ASSERT_EQ(root,0);ASSERT_EQ(stats.rejected,0);
    primitives[0].p[3]=2;stats=(compute_candidate_stats){0};
    ASSERT_EQ(alea_compute_point_candidates(&scene,bounds,signs,&next,&root,&stats),0);
    ASSERT_EQ(root,COMPUTE_NONE);ASSERT_EQ(stats.rejected,1);
    cell.count=POINT_PROOF_BUDGET+1;stats=(compute_candidate_stats){0};
    ASSERT_EQ(alea_compute_point_candidates(&scene,bounds,signs,&next,&root,&stats),0);
    ASSERT_EQ(root,0);ASSERT_EQ(stats.instructions,0);ASSERT_EQ(stats.budget_fallbacks,1);
}
TEST(compact_code_exact_in_place_and_wide_fallback) {
    compute_node original[]={{COMPUTE_LEAF,0x0fffffffu},{COMPUTE_SKIP_UNION,12345},
        {COMPUTE_CONSTANT,1},{COMPUTE_COMPLEMENT,0},{15,0}};
    enum {N=sizeof(original)/sizeof(*original)};
    compute_node nodes[N];memcpy(nodes,original,sizeof(nodes));
    uint32_t words[N];
    ASSERT_EQ(alea_compute_code_compact(nodes,N,words),1);
    compute_scene scene={.compact_nodes=words};
    for(size_t i=0;i<N;++i) {
        compute_node decoded=alea_compute_node_at(&scene,i);
        ASSERT_EQ(decoded.op,original[i].op);ASSERT_EQ(decoded.arg,original[i].arg);
    }
    ASSERT_EQ(memcmp(nodes,original,sizeof(nodes)),0);
    ASSERT_EQ(alea_compute_code_compact(nodes,N,(uint32_t*)nodes),1);
    ASSERT_EQ(memcmp(nodes,words,sizeof(words)),0);
    for(unsigned invalid=0;invalid<2;++invalid) {
        memcpy(nodes,original,sizeof(nodes));
        if(invalid)nodes[N-1].op=16;else nodes[N-1].arg=0x10000000u;
        compute_node copy[N];memcpy(copy,nodes,sizeof(copy));
        ASSERT_EQ(alea_compute_code_compact(nodes,N,(uint32_t*)nodes),0);
        ASSERT_EQ(memcmp(nodes,copy,sizeof(nodes)),0);
    }
    memcpy(nodes,original,sizeof(nodes));alea_interrupt();
    int result=alea_compute_code_compact(nodes,N,(uint32_t*)nodes);
    alea_clear_interrupt();
    ASSERT_EQ(result,-1);ASSERT_EQ(alea_error_code(),ALEA_ERR_INTERRUPTED);
    ASSERT_EQ(memcmp(nodes,original,sizeof(nodes)),0);
    ASSERT_EQ(alea_compute_code_compact(nodes,N,(uint32_t*)nodes),1);
    ASSERT_EQ(memcmp(nodes,words,sizeof(words)),0);
    ASSERT_EQ(alea_compute_code_compact(NULL,0,NULL),1);
}
TEST_MAIN()
