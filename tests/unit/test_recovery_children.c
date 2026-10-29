// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
// SPDX-License-Identifier: MPL-2.0
#include "alea_test.h"
#include "alea.h"
#include "alea_mcnp.h"
#include "core/alea_system.h"
#include "raycast/raycast.h"
#include "util/alea_alloc.h"

static void compare_scan(alea_system_t* sys, const alea_ray_t* ray,
                         double limit, int nearest) {
    alea_raycast_result_t prepared, fallback;
    alea_raycast_result_init(&prepared);
    alea_raycast_result_init(&fallback);
    prepared.breakpoint_nearest_only = fallback.breakpoint_nearest_only = nearest;
    int a = alea_raycast_validation_breakpoints_reuse_nocache(
        sys, ray, 0, limit, nearest ? 1 : 100000, &prepared);
    unsigned state = atomic_fetch_and(&sys->query_cache_state,
                                      ~ALEA_CACHE_RECOVERY_CHILDREN);
    int b = alea_raycast_validation_breakpoints_reuse_nocache(
        sys, ray, 0, limit, nearest ? 1 : 100000, &fallback);
    atomic_store(&sys->query_cache_state, state);
    ASSERT_EQ(a, b);
    ASSERT_EQ(prepared.breakpoint_failed, fallback.breakpoint_failed);
    ASSERT_EQ(prepared.hits.count, fallback.hits.count);
    for (size_t i = 0; i < prepared.hits.count; i++) {
        const alea_ray_hit_t* x = &prepared.hits.data[i];
        const alea_ray_hit_t* y = &fallback.hits.data[i];
        ASSERT_EQ(x->t, y->t);
        ASSERT_EQ(x->surface_id, y->surface_id);
        ASSERT_EQ(x->primitive_id, y->primitive_id);
        ASSERT_EQ(x->nx, y->nx); ASSERT_EQ(x->ny, y->ny); ASSERT_EQ(x->nz, y->nz);
    }
    alea_raycast_result_free(&prepared);
    alea_raycast_result_free(&fallback);
}

TEST(recovery_children_preserve_nested_fills_and_rectangular_hex_lattices) {
    const char* paths[] = {
        "tests/data/mcnp_deep_translated_fills.mcnp",
        "tests/data/mcnp_nested_lattice_transformed.mcnp",
        "tests/data/mcnp_lattice_repeating.mcnp",
        "tests/data/mcnp_hex_lattice.mcnp"
    };
    for (size_t p = 0; p < sizeof(paths)/sizeof(paths[0]); p++) {
        mcnp_model_t* model = mcnp_load(paths[p]);
        ASSERT_NOT_NULL(model);
        alea_system_t* sys = model->sys;
        ASSERT_EQ(alea_prepare_query_acceleration(sys), 0);
        for (size_t u = 0; u < sys->universes.count; u++) {
            const alea_universe_t* univ = &sys->universes.data[u];
            ASSERT(univ->recovery_children_ready);
            size_t next = 0;
            for (size_t c = 0; c < univ->cell_indices.count; c++) {
                size_t ci = univ->cell_indices.data[c];
                const alea_cell_entry_t* cell = &sys->cells.data[ci];
                if (cell->fill_universe > 0 || (cell->lat_type && cell->lat_fill)) {
                    ASSERT(next < univ->recovery_child_count);
                    ASSERT_EQ(univ->recovery_child_cells[next++], ci);
                }
            }
            ASSERT_EQ(next, univ->recovery_child_count);
        }
        for (int axis = 0; axis < 3; axis++) for (int j = 0; j < 9; j++) {
            double o[] = {0.1, 0.1, 0.1}, d[] = {0, 0, 0};
            o[axis] = j % 2 ? 12 : -12;
            o[(axis+1)%3] = (j-4)*0.4;
            d[axis] = j % 2 ? -1 : 1;
            alea_ray_t ray;
            ASSERT_EQ(alea_ray_init(&ray,o[0],o[1],o[2],d[0],d[1],d[2]), 0);
            compare_scan(sys, &ray, 24, 1);
            compare_scan(sys, &ray, 24, 0);
        }
        mcnp_model_destroy(model);
    }
}

TEST(recovery_children_keep_nearer_occurrences_aliases_clipping_and_failures) {
    alea_system_t* sys = alea_create(); ASSERT_NOT_NULL(sys);
    int plane = alea_plane_surface(sys, 1, 1, 0, 0, -2);
    int alias = alea_plane_surface(sys, 99, 1, 0, 0, -2);
    ASSERT(alea_add_cell(sys,1,alea_halfspace(sys,plane,-1),ALEA_MATERIAL_VOID,0,8)>=0);
    ASSERT(alea_add_cell(sys,2,alea_halfspace(sys,alias,1),ALEA_MATERIAL_VOID,0,8)>=0);
    int box = alea_box_surface(sys,3,-10,10,-2,2,-2,2);
    alea_node_id_t region = alea_halfspace(sys,box,-1);
    int later = -1;
    for (int i = 0; i < 32; i++) {
        int c = alea_add_cell(sys,100+i,region,ALEA_MATERIAL_VOID,0,0);
        ASSERT(c>=0);
        if (i == 4 || i == 28) {
            ASSERT_EQ(alea_set_fill(sys,c,8,0),0);
            later = c;
        }
    }
    /* Rotate about x (the plane remains x=2), then translate to x=1. */
    const double transform[] = {-1,0,0, 1,0,0, 0,0,-1, 0,1,0};
    ASSERT_EQ(alea_add_transform(sys,9,transform,12,0),0);
    ASSERT_EQ(alea_set_fill(sys,later,8,9),0);
    ASSERT_EQ(alea_prepare_query_acceleration(sys),0);
    ASSERT_EQ(alea_get_universe(sys,0)->recovery_child_count,2);
    ASSERT(alea_get_universe(sys,8)->recovery_children_ready);
    ASSERT_EQ(alea_get_universe(sys,8)->recovery_child_count,0);
    for (int reverse = 0; reverse < 2; reverse++) {
        alea_ray_t ray;
        ASSERT_EQ(alea_ray_init(&ray,reverse?4:.5,0,0,reverse?-1:1,0,0),0);
        compare_scan(sys,&ray,5,1);
        alea_raycast_result_t result; alea_raycast_result_init(&result);
        result.breakpoint_nearest_only = 1;
        ASSERT_EQ(alea_raycast_validation_breakpoints_reuse_nocache(sys,&ray,0,5,1,&result),0);
        ASSERT_EQ(result.hits.count,1);
        ASSERT_EQ(result.hits.data[0].t,reverse?2:.5);
        ASSERT_EQ(result.hits.data[0].surface_id,1);
        ASSERT_EQ(result.hits.data[0].nx,1);
        alea_interrupt();
        ASSERT(alea_raycast_validation_breakpoints_reuse_nocache(sys,&ray,0,5,1,&result)!=0);
        alea_clear_interrupt();
        ASSERT_EQ(alea_raycast_validation_breakpoints_reuse_nocache(sys,&ray,0,5,1,&result),0);
        ASSERT_EQ(result.hits.data[0].t,reverse?2:.5);
        alea_raycast_result_free(&result);
    }
    /* A missed parent must still reject its child; undefined fills must still fail. */
    alea_ray_t ray;
    ASSERT_EQ(alea_ray_init(&ray,.5,3,0,1,0,0),0);
    compare_scan(sys,&ray,5,1);
    alea_raycast_result_t result; alea_raycast_result_init(&result);
    result.breakpoint_nearest_only=1;
    ASSERT_EQ(alea_raycast_validation_breakpoints_reuse_nocache(sys,&ray,0,5,1,&result),0);
    /* The global diagnostic surface pass still sees the unplaced plane at
     * x=2; the clipped translated occurrence at x=1 must not replace it. */
    ASSERT_EQ(result.hits.count,1);
    ASSERT_EQ(result.hits.data[0].t,1.5);
    ASSERT_EQ(alea_set_fill(sys,later,9999,0),0);
    ASSERT_EQ(alea_prepare_query_acceleration(sys),0);
    ASSERT_EQ(alea_ray_init(&ray,.5,0,0,1,0,0),0);
    compare_scan(sys,&ray,5,1);
    ASSERT(alea_raycast_validation_breakpoints_reuse_nocache(sys,&ray,0,5,1,&result)!=0);
    alea_raycast_result_free(&result);
    alea_destroy(sys);
}

TEST(recovery_children_rebuild_clone_and_edit_ownership) {
    mcnp_model_t* model = mcnp_load("tests/data/mcnp_deep_translated_fills.mcnp");
    ASSERT_NOT_NULL(model);
    alea_system_t* sys = model->sys;
    ASSERT_EQ(alea_prepare_query_acceleration(sys),0);
    alea_system_t* clone = alea_clone(sys); ASSERT_NOT_NULL(clone);
    ASSERT_EQ(alea_prepare_query_acceleration(clone),0);
    ASSERT(alea_get_universe(sys,0)->recovery_child_cells !=
           alea_get_universe(clone,0)->recovery_child_cells);
    ASSERT_EQ(alea_build_universe_index(sys),0);
    ASSERT(!alea_system_query_cache_ready(sys,ALEA_CACHE_RECOVERY_CHILDREN));
    ASSERT(!alea_get_universe(sys,0)->recovery_children_ready);
    ASSERT_EQ(alea_prepare_query_acceleration(sys),0);
    int parent = alea_cell_find(sys,100);
    ASSERT_EQ(alea_set_fill(sys,parent,0,0),0);
    ASSERT(!alea_system_query_cache_ready(sys,ALEA_CACHE_RECOVERY_CHILDREN));
    ASSERT_EQ(alea_prepare_query_acceleration(sys),0);
    ASSERT_EQ(alea_get_universe(sys,0)->recovery_child_count,0);
    ASSERT_EQ(alea_set_fill(sys,parent,1,0),0);
    ASSERT_EQ(alea_cell_set_universe(sys,parent,99),0);
    ASSERT_EQ(alea_prepare_query_acceleration(sys),0);
    ASSERT_EQ(alea_get_universe(sys,99)->recovery_child_count,1);
    int sphere = alea_sphere_surface(sys,999,0,0,0,30);
    ASSERT_EQ(alea_cell_set_region(sys,parent,alea_halfspace(sys,sphere,-1)),0);
    ASSERT_EQ(alea_prepare_query_acceleration(sys),0);
    ASSERT_EQ(alea_get_universe(sys,99)->recovery_child_count,1);
    mcnp_model_destroy(model);
    alea_ray_t ray; ASSERT_EQ(alea_ray_init(&ray,1,0,0,1,0,0),0);
    alea_raycast_result_t result; alea_raycast_result_init(&result);
    result.breakpoint_nearest_only=1;
    ASSERT_EQ(alea_raycast_validation_breakpoints_reuse_nocache(clone,&ray,0,3,1,&result),0);
    ASSERT_EQ(result.hits.count,1); ASSERT_EQ(result.hits.data[0].t,1);
    alea_raycast_result_free(&result); alea_destroy(clone);
}

TEST(recovery_children_budget_incomplete_and_interrupted_preparation_fallback) {
    mcnp_model_t* model = mcnp_load("tests/data/mcnp_deep_translated_fills.mcnp");
    ASSERT_NOT_NULL(model);
    alea_system_t* sys = model->sys;
    ASSERT_EQ(alea_prepare_query_acceleration(sys),0);
    alea_system_invalidate_query_caches(sys,ALEA_CACHE_RECOVERY_CHILDREN);
    alea_interrupt();
    ASSERT(alea_system_prepare_query_caches(sys,ALEA_CACHE_RECOVERY_CHILDREN)!=0);
    ASSERT(!alea_system_query_cache_ready(sys,ALEA_CACHE_RECOVERY_CHILDREN));
    alea_clear_interrupt();
    alea_alloc_fail_after(0);
    int rc = alea_system_prepare_query_caches(sys,ALEA_CACHE_RECOVERY_CHILDREN);
    alea_alloc_fail_after(-1);
    ASSERT_EQ(rc,0);
    ASSERT(!alea_get_universe(sys,0)->recovery_children_ready);
    alea_ray_t oom_ray; ASSERT_EQ(alea_ray_init(&oom_ray,1,0,0,1,0,0),0);
    compare_scan(sys,&oom_ray,3,1);
    alea_system_invalidate_query_caches(sys,ALEA_CACHE_RECOVERY_CHILDREN);
    /* An oversized internal universe exercises the real cap without constructing
     * a million geometric objects. Duplicate entries retain their original order. */
    alea_universe_t* root = (alea_universe_t*)alea_get_universe(sys,0);
    size_t* original = root->cell_indices.data;
    size_t original_count = root->cell_indices.count;
    size_t count = 8*1024*1024/sizeof(size_t)+1;
    size_t* repeated = malloc(count*sizeof(size_t)); ASSERT_NOT_NULL(repeated);
    for (size_t i=0;i<count;i++) repeated[i]=original[0];
    root->cell_indices.data=repeated; root->cell_indices.count=count;
    ASSERT_EQ(alea_system_prepare_query_caches(sys,ALEA_CACHE_RECOVERY_CHILDREN),0);
    ASSERT(!root->recovery_children_ready);
    ASSERT(root->recovery_child_cells==NULL);
    size_t retained=0;
    for(size_t u=0;u<sys->universes.count;u++)
        retained += sys->universes.data[u].recovery_child_count*sizeof(size_t);
    ASSERT(retained<=8*1024*1024);
    root->cell_indices.data=original; root->cell_indices.count=original_count;
    free(repeated);
    alea_ray_t ray; ASSERT_EQ(alea_ray_init(&ray,1,0,0,1,0,0),0);
    compare_scan(sys,&ray,3,1); /* Published global cache, unavailable local list. */
    alea_system_invalidate_query_caches(sys,ALEA_CACHE_RECOVERY_CHILDREN);
    size_t saved=original[0]; original[0]=sys->cells.count;
    ASSERT_EQ(alea_system_prepare_query_caches(sys,ALEA_CACHE_RECOVERY_CHILDREN),0);
    ASSERT(!root->recovery_children_ready);
    compare_scan(sys,&ray,3,1); /* Malformed entries must not disappear. */
    original[0]=saved;
    alea_system_invalidate_query_caches(sys,ALEA_CACHE_RECOVERY_CHILDREN);
    ASSERT_EQ(alea_system_prepare_query_caches(sys,ALEA_CACHE_RECOVERY_CHILDREN),0);
    ASSERT(root->recovery_children_ready);
    compare_scan(sys,&ray,3,1);
    mcnp_model_destroy(model);
}

TEST_MAIN()
