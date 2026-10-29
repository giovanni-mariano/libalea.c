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
TEST_MAIN()
