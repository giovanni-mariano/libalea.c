// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
// SPDX-License-Identifier: MPL-2.0
#include "alea_test.h"
#include "alea.h"
#include "core/alea_system.h"
#include "core/alea_eval.h"
#include "core/alea_spatial_hier.h"
#include <float.h>

static int add_sphere(alea_system_t* s, int id, double x, double r, int universe) {
    int surface=alea_sphere_surface(s,id,x,0,0,r);
    if(surface<0) return -1;
    return alea_add_cell(s,id,alea_halfspace(s,surface,-1),ALEA_MATERIAL_VOID,0,universe);
}
static int linear(alea_system_t* s, double x, double y, double z) {
    for(size_t i=0;i<s->cells.count;i++)
        if(s->cells.data[i].universe_id==0 &&
           alea_contains_point(s,s->cells.data[i].root_node_id,x,y,z))
            return s->cells.data[i].mc_cell_id;
    return -1;
}
static void check_point(alea_system_t* s, double x, double y, double z) {
    int expected=linear(s,x,y,z), id=-1;
    int rc=alea_find_cell_at(s,x,y,z,&id,NULL);
    ASSERT_EQ(rc,expected<0 ? -1 : 0);
    if(rc==0) ASSERT_EQ(id,expected);
}
TEST(point_query_preserves_model_order_after_cached_overlapping_hit) {
    alea_system_t* s=alea_create(); ASSERT_NOT_NULL(s);
    ASSERT(add_sphere(s,900,0,1,0)>=0);
    ASSERT(add_sphere(s,100,1.5,1,0)>=0);
    ASSERT(add_sphere(s,200,10,1,0)>=0);
    ASSERT_EQ(alea_prepare_query_acceleration(s),0);
    alea_cell_hit_t hit;
    ASSERT_EQ(alea_find_all_cells(s,1.5,0,0,&hit,1),1);
    ASSERT_EQ(hit.cell_id,100);
    int id=-1;
    ASSERT_EQ(alea_find_cell_at(s,.75,0,0,&id,NULL),0); ASSERT_EQ(id,900);
    for(int i=0;i<=160;i++) check_point(s,-2+i*.1,0,0);
    check_point(s,nextafter(1,0),0,0);
    check_point(s,1,0,0); check_point(s,nextafter(1,INFINITY),0,0);
    alea_destroy(s);
}
TEST(point_query_keeps_large_bounds_unbounded_cells_and_shared_faces) {
    alea_system_t* s=alea_create(); ASSERT_NOT_NULL(s);
    for(int i=0;i<4;i++) {
        double x=(i<2 ? 1 : -1)*1e39 + (i%2)*1e38;
        int box=alea_box_surface(s,i+1,x,x+1e37,-1,1,-1,1);
        ASSERT(box>=0);
        ASSERT(alea_add_cell(s,100-i,alea_halfspace(s,box,-1),ALEA_MATERIAL_VOID,0,0)>=0);
    }
    int plane=alea_plane_surface(s,10,0,1,0,0);
    ASSERT(alea_add_cell(s,10,alea_halfspace(s,plane,1),ALEA_MATERIAL_VOID,0,0)>=0);
    ASSERT_EQ(alea_prepare_query_acceleration(s),0);
    check_point(s,1.1e39+5e36,0,0); check_point(s,-.9e39+5e36,0,0);
    check_point(s,0,2,0); check_point(s,0,-2,0); check_point(s,0,0,0);
    alea_destroy(s);
    s=alea_create(); ASSERT_NOT_NULL(s);
    for(int i=0;i<3;i++) {
        int box=alea_box_surface(s,i+1,i,i+1,-1,1,-1,1);
        ASSERT(alea_add_cell(s,30-i,alea_halfspace(s,box,-1),ALEA_MATERIAL_VOID,0,0)>=0);
    }
    ASSERT_EQ(alea_prepare_query_acceleration(s),0);
    check_point(s,1,0,0); check_point(s,nextafter(1,0),0,0);
    check_point(s,nextafter(1,INFINITY),0,0);
    alea_destroy(s);
}
TEST(point_query_acceleration_filters_late_hits_and_keeps_linear_fallback) {
    alea_system_t* s=alea_create(); ASSERT_NOT_NULL(s);
    for(int i=0;i<1024;i++) ASSERT(add_sphere(s,i+1,i*4,1,0)>=0);
    ASSERT_EQ(alea_build_universe_index(s),0);
    check_point(s,4092,0,0); /* Universe index only: no BLAS. */
    ASSERT_EQ(alea_prepare_query_acceleration(s),0);
    alea_hier_debug_candidates_reset(); alea_perf_reset();
    int id=-1;
    ASSERT_EQ(alea_find_cell_at(s,4092,0,0,&id,NULL),0); ASSERT_EQ(id,1024);
    ASSERT(alea_hier_debug_candidates_get()<16);
    ASSERT(alea_perf_get().primitive_evaluations<16);
    check_point(s,4094,0,0); check_point(s,0,0,0);
    alea_destroy(s);
}
TEST(point_query_tracks_fill_changes_clones_and_region_edits) {
    alea_system_t* s=alea_create(); ASSERT_NOT_NULL(s);
    int root=add_sphere(s,1,0,5,0);
    ASSERT(root>=0); ASSERT(add_sphere(s,2,0,1,7)>=0);
    ASSERT_EQ(alea_prepare_query_acceleration(s),0);
    int id=-1;
    ASSERT_EQ(alea_find_cell_at(s,0,0,0,&id,NULL),0); ASSERT_EQ(id,1);
    ASSERT_EQ(alea_set_fill(s,root,7,0),0);
    ASSERT_EQ(alea_prepare_query_acceleration(s),0);
    ASSERT_EQ(alea_find_cell_at(s,0,0,0,&id,NULL),0); ASSERT_EQ(id,2);
    ASSERT_EQ(alea_find_cell_at(s,3,0,0,&id,NULL),0); ASSERT_EQ(id,1); /* Undefined fill. */
    alea_system_t* clone=alea_clone(s); ASSERT_NOT_NULL(clone);
    ASSERT_EQ(alea_set_fill(s,root,0,0),0);
    ASSERT_EQ(alea_prepare_query_acceleration(s),0);
    ASSERT_EQ(alea_prepare_query_acceleration(clone),0);
    ASSERT_EQ(alea_find_cell_at(s,0,0,0,&id,NULL),0); ASSERT_EQ(id,1);
    ASSERT_EQ(alea_find_cell_at(clone,0,0,0,&id,NULL),0); ASSERT_EQ(id,2);
    int box=alea_box_surface(s,99,20,21,-1,1,-1,1);
    ASSERT_EQ(alea_cell_set_region(s,root,alea_halfspace(s,box,-1)),0);
    ASSERT_EQ(alea_prepare_query_acceleration(s),0);
    check_point(s,0,0,0); check_point(s,20.5,0,0);
    ASSERT_EQ(alea_cell_set_universe(s,root,9),0);
    ASSERT_EQ(alea_prepare_query_acceleration(s),0);
    check_point(s,20.5,0,0);
    ASSERT_EQ(alea_cell_set_universe(s,root,0),0);
    ASSERT_EQ(alea_prepare_query_acceleration(s),0);
    check_point(s,20.5,0,0);
    alea_destroy(clone); alea_destroy(s);
}
TEST(point_query_dense_candidates_preserve_last_winner) {
    alea_system_t* s=alea_create(); ASSERT_NOT_NULL(s);
    int outer=alea_sphere_surface(s,1,0,0,0,10);
    int inner=alea_sphere_surface(s,2,0,0,0,1);
    alea_node_id_t shell=alea_intersection(s,alea_halfspace(s,outer,-1),
                                           alea_halfspace(s,inner,1));
    for(int i=0;i<128;i++)
        ASSERT(alea_add_cell(s,1000-i,shell,ALEA_MATERIAL_VOID,0,0)>=0);
    ASSERT(add_sphere(s,3,0,.5,0)>=0);
    ASSERT_EQ(alea_prepare_query_acceleration(s),0);
    check_point(s,0,0,0); check_point(s,2,0,0); check_point(s,20,0,0);
    alea_destroy(s);
}
TEST(empty_lattice_slot_does_not_reenter_root_universe) {
    alea_system_t* s=alea_create();int box=alea_box_surface(s,1,-1,1,-1,1,-1,1);
    int at=alea_add_cell(s,20,alea_halfspace(s,box,-1),ALEA_MATERIAL_VOID,0,0);
    ASSERT(at>=0);alea_cell_entry_t* c=&s->cells.data[at];
    c->lat_type=1;c->lat_fill_count=1;c->lat_fill=calloc(1,sizeof(int));ASSERT_NOT_NULL(c->lat_fill);
    for(int a=0;a<3;++a){c->lat_pitch[a]=2;c->lat_lower_left[a]=-1;c->lat_fill_dims[2*a]=c->lat_fill_dims[2*a+1]=0;}
    ASSERT(add_sphere(s,99,0,2,0)>=0);
    ASSERT_EQ(alea_prepare_query_acceleration(s),0);alea_cell_hit_t hit;
    ASSERT_EQ(alea_find_deepest_cell_hit_at_point(s,0,0,0,&hit),0);
    ASSERT_EQ(hit.cell_id,20);ASSERT_EQ(hit.depth,0);
    ASSERT_TRUE(hit.resolution_flags&ALEA_RESOLVE_UNDEFINED_FILL);
    alea_destroy(s);
}
TEST_MAIN()
