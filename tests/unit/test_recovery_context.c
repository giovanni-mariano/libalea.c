// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
// SPDX-License-Identifier: MPL-2.0
#include "alea_test.h"
#include "alea.h"
#include "core/alea_system.h"
#include "core/alea_universe.h"
#include "raycast/recovery_context.h"
#include "raycast/recovery_planes.h"

TEST(recovery_root_records_preserve_distance_bits_and_card_order) {
    const double distances[] = {-0.0, DBL_TRUE_MIN, nextafter(1.0, 2.0), DBL_MAX};
    recovery_root_t roots[4];
    for (unsigned i = 0; i < 4; ++i) {
        roots[i] = recovery_root_make(distances[i], 100000, i);
        double actual = recovery_root_distance(&roots[i]);
        ASSERT_EQ(memcmp(&actual, &distances[i], sizeof(actual)), 0);
        ASSERT_EQ(roots[i].order, 400000u + i);
    }
    ASSERT_EQ(sizeof(roots), 4 * (sizeof(double) + sizeof(uint32_t)));
}

TEST(recovery_roots_preserve_alias_order_and_strict_current_boundary) {
    recovery_root_t roots[]={recovery_root_make(2,9,0), recovery_root_make(1,7,0),
        recovery_root_make(2,3,1), recovery_root_make(2,3,0),
        recovery_root_make(nextafter(1,INFINITY),0,0), recovery_root_make(3,0,0)};
    qsort(roots,6,sizeof(*roots),recovery_root_compare);
    recovery_context_entry_t entry={.roots=roots,.root_count=6};
    const recovery_root_t* hit=recovery_root_first(&entry,1,2,1);
    ASSERT_NOT_NULL(hit);
    ASSERT_EQ(recovery_root_distance(hit),2);
    ASSERT_EQ(hit->order / 4u,3);
    ASSERT_EQ(hit->order % 4u,0);
    ASSERT(recovery_root_first(&entry,1,nextafter(2,0),1)==NULL);
    ASSERT_EQ(recovery_root_distance(recovery_root_first(&entry,2,3,2)),3);
    ASSERT_EQ(recovery_root_first(&entry,1,1,0)->order / 4u,7);
    entry.root_count=0;
    ASSERT(recovery_root_first(&entry,0,4,0)==NULL);
}

TEST(recovery_context_separates_local_rays_and_transform_occurrences) {
    uint32_t cards[]={1,2};
    alea_ray_t ray={.ox=1,.dx=1,.inv_dx=1,.inv_dy=INFINITY,.inv_dz=INFINITY};
    recovery_context_t* table=recovery_context_create(8);
    ASSERT_NOT_NULL(table);
    recovery_context_key_t key;
    ASSERT(recovery_context_key(&key,cards,2,&ray,NULL));
    recovery_context_entry_t* original=recovery_context_find(table,&key);
    ASSERT_NOT_NULL(original);
    original->key=key; original->complete=true;
    original->lo=1; original->hi=5;
    ASSERT(recovery_context_contains(original,2,4));
    ASSERT(!recovery_context_contains(original,0,4));
    ASSERT(!recovery_context_contains(original,2,6));
    original->roots=malloc(sizeof(*original->roots));
    ASSERT_NOT_NULL(original->roots);
    original->roots[0]=recovery_root_make(3,0,0);
    original->root_count=1; original->roots_complete=true;
    ray.ox=2;
    ASSERT(recovery_context_key(&key,cards,2,&ray,NULL));
    ASSERT(recovery_context_find(table,&key)!=original);
    ray.ox=1;
    alea_matrix_t transform={0};
    ASSERT(recovery_context_key(&key,cards,2,&ray,&transform));
    ASSERT(recovery_context_find(table,&key)!=original);
    ray.ox=NAN;
    ASSERT(!recovery_context_key(&key,cards,2,&ray,NULL));
    recovery_context_destroy(table);
    recovery_context_destroy(NULL);
}

TEST(recovery_root_binary_search_matches_ordered_scan) {
    recovery_root_t original[512], sorted[512];
    for (unsigned i=0;i<512;i++) original[i]=recovery_root_make((i*73)%137,i/4,i%4);
    memcpy(sorted,original,sizeof(sorted));
    qsort(sorted,512,sizeof(*sorted),recovery_root_compare);
    recovery_context_entry_t entry={.roots=sorted,.root_count=512};
    for (unsigned k=0;k<500;k++) {
        double lo=(k*17)%137, hi=lo+(k%20), current=lo-(k%2);
        const recovery_root_t* expected=NULL;
        for (unsigned i=0;i<512;i++) {
            const recovery_root_t* r=&original[i];
            double t = recovery_root_distance(r);
            if (t >= lo && t <= hi && t > nextafter(current, INFINITY) &&
                (!expected || t < recovery_root_distance(expected))) expected = r;
        }
        const recovery_root_t* actual=recovery_root_first(&entry,lo,hi,current);
        ASSERT((expected==NULL)==(actual==NULL));
        if (expected) {
            ASSERT_EQ(recovery_root_distance(actual),recovery_root_distance(expected));
            ASSERT_EQ(actual->order,expected->order);
        }
    }
}

TEST(recovery_cached_windows_match_scan_with_planes_aliases_and_torus) {
    alea_system_t* sys=alea_create();
    ASSERT_NOT_NULL(sys);
    int enclosing=alea_sphere_surface(sys,1,0,0,0,100);
    ASSERT(enclosing>=0);
    int parent=alea_add_cell(sys,1,alea_halfspace(sys,enclosing,-1),ALEA_MATERIAL_VOID,0,0);
    ASSERT(parent>=0);
    double translation[]={20,12,0};
    ASSERT_EQ(alea_add_transform(sys,1,translation,3,0),0);
    ASSERT_EQ(alea_set_fill(sys,parent,1,1),0);
    for (int i=0;i<180;i++) {
        /* Equal roots with distinct surface cards exercise tie provenance. */
        int surface=alea_plane_surface(sys,100+i,1,0,0,-(-10+(i/2)*.2));
        ASSERT(surface>=0);
        ASSERT(alea_add_cell(sys,100+i,alea_halfspace(sys,surface,-1),ALEA_MATERIAL_VOID,0,1)>=0);
    }
    /* One full packet and a partial packet, with distinct cards sharing roots. */
    for (int i=0;i<11;i++) {
        int surf=alea_quadric_surface(sys,500+i,1,2,3,.125,0,0,14,0,0,48-(i/2)*.01);
        ASSERT(surf>=0);
        ASSERT(alea_add_cell(sys,500+i,alea_halfspace(sys,surf,-1),ALEA_MATERIAL_VOID,0,1)>=0);
    }
    int torus=alea_torus_z_surface(sys,400,0,0,0,3,.75);
    ASSERT(torus>=0);
    ASSERT(alea_add_cell(sys,400,alea_halfspace(sys,torus,-1),ALEA_MATERIAL_VOID,0,1)>=0);
    int far_torus=alea_torus_z_surface(sys,401,35,0,0,3,.75);
    ASSERT(far_torus>=0);
    ASSERT(alea_add_cell(sys,401,alea_halfspace(sys,far_torus,-1),ALEA_MATERIAL_VOID,0,1)>=0);
    int earlier_plane=alea_plane_surface(sys,402,1,0,0,42);
    ASSERT(earlier_plane>=0);
    ASSERT(alea_add_cell(sys,402,alea_halfspace(sys,earlier_plane,-1),ALEA_MATERIAL_VOID,0,1)>=0);
    ASSERT_EQ(alea_prepare_query_acceleration(sys),0);
    const alea_universe_t* filled = alea_get_universe(sys,1);
    ASSERT_NOT_NULL(filled);
    ASSERT_NOT_NULL(filled->recovery_planes);
    ASSERT_EQ(filled->recovery_planes->torus_cards,2);
    ASSERT_EQ(filled->recovery_planes->quadric_count,11);
    alea_ray_t ray;
    ASSERT_EQ(alea_ray_init(&ray,-30,12.2,0,1,0,0),0);
    recovery_context_t* table=recovery_context_create(128);
    ASSERT_NOT_NULL(table);
    for (unsigned k=0;k<65;k++) {
        /* Populate a bounded cache, then query across both ends of it.
         * The farther torus is the next crossing on the extended ray. */
        double lo=35+k*.7, hi= k < 24 ? 52 : 100;
        /* The first retry ends before any crossing. Its cache can still
         * answer later intervals on this same finite ray. */
        if (k==0) hi=lo+.01;
        if (k==63) { lo=7; hi=8; }
        if (k==64) { lo=51.99; hi=52; }
        for (size_t i=0;i<table->capacity;i++) {
            table->entries[i].lo=INFINITY; table->entries[i].hi=-INFINITY;
        }
        alea_raycast_result_t reference,cached;
        alea_raycast_result_init(&reference); alea_raycast_result_init(&cached);
        reference.breakpoint_nearest_only=cached.breakpoint_nearest_only=1;
        reference.breakpoint_t_min=cached.breakpoint_t_min=lo;
        cached.recovery_context=table;
        if (k==0) cached.recovery_root_limit=52;
        ASSERT_EQ(alea_raycast_validation_breakpoints_reuse_nocache(sys,&ray,lo,hi,1,&reference),0);
        ASSERT_EQ(alea_raycast_validation_breakpoints_reuse_nocache(sys,&ray,lo,hi,1,&cached),0);
        ASSERT_EQ(reference.hits.count,cached.hits.count);
        if (k==0) ASSERT_EQ(cached.hits.count,0);
        for (size_t i=0;i<reference.hits.count;i++) {
            const alea_ray_hit_t *a=&reference.hits.data[i], *b=&cached.hits.data[i];
            ASSERT_EQ(a->t,b->t); ASSERT_EQ(a->surface_id,b->surface_id);
            ASSERT_EQ(a->primitive_id,b->primitive_id);
            ASSERT_EQ(a->nx,b->nx); ASSERT_EQ(a->ny,b->ny); ASSERT_EQ(a->nz,b->nz);
        }
        if (k==62) {
            ASSERT_EQ(reference.hits.count,1);
            ASSERT_EQ(reference.hits.data[0].surface_id,401);
        }
        if (k==63 || k==64) {
            ASSERT_EQ(reference.hits.count,1);
            ASSERT_EQ(reference.hits.data[0].t,hi);
        }
        alea_raycast_result_free(&reference); alea_raycast_result_free(&cached);
    }
    size_t built=0, exceptional=0;
    for (size_t i=0;i<table->capacity;i++) {
        const recovery_context_entry_t* entry=&table->entries[i];
        if (!entry->roots_complete) continue;
        built++;
        exceptional+=entry->exceptional_count;
        ASSERT_EQ(entry->roots_hi,52);
        for (size_t j=0;j<entry->root_count;j++) {
            ASSERT(recovery_root_distance(&entry->roots[j])>=entry->roots_lo);
            ASSERT(recovery_root_distance(&entry->roots[j])<=entry->roots_hi);
        }
    }
    ASSERT_EQ(built,1); ASSERT_EQ(exceptional,1);
    ASSERT(table->root_bytes<=RECOVERY_ROOT_BUDGET);
    recovery_context_destroy(table);

    /* Fail the packet append after the plane scratch and final index succeed. */
    alea_system_invalidate_query_caches(sys,ALEA_CACHE_RECOVERY_SURFACES);
    alea_alloc_fail_after(2);
    int rc=alea_system_prepare_query_caches(sys,ALEA_CACHE_RECOVERY_SURFACES);
    alea_alloc_fail_after(-1);
    ASSERT_EQ(rc,0);
    ASSERT_NOT_NULL(filled->recovery_planes);
    ASSERT_EQ(filled->recovery_planes->quadric_count,0);
    ASSERT_EQ(filled->recovery_planes->torus_cards,2);
    alea_system_invalidate_query_caches(sys,ALEA_CACHE_RECOVERY_SURFACES);
    ASSERT_EQ(alea_system_prepare_query_caches(sys,ALEA_CACHE_RECOVERY_SURFACES),0);
    ASSERT_EQ(filled->recovery_planes->quadric_count,11);
    alea_destroy(sys);
}

TEST(transport_fill_filter_preserves_duplicate_hits_budget_and_unindexed_geometry) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);
    int sphere = alea_sphere_surface(sys, 1, 0, 0, 0, 100);
    int plane = alea_plane_surface(sys, 100, 1, 0, 0, 0);
    int alias = alea_plane_surface(sys, 101, 1, 0, 0, 0);
    int offset = alea_plane_surface(sys, 102, 1, 0, 0, -2);
    ASSERT(sphere >= 0 && plane >= 0 && alias >= 0 && offset >= 0);
    int parent = alea_add_cell(sys, 1, alea_halfspace(sys, sphere, -1), ALEA_MATERIAL_VOID, 0, 0);
    double translation[] = {10, 0, 0};
    ASSERT(parent >= 0);
    ASSERT_EQ(alea_add_transform(sys, 1, translation, 3, 0), 0);
    ASSERT_EQ(alea_set_fill(sys, parent, 1, 1), 0);
    for (int i = 0; i < 3; ++i) {
        alea_node_id_t leaf = alea_halfspace(sys, i == 2 ? alias : plane, -1);
        ASSERT(alea_add_cell(sys, 2+i, leaf, ALEA_MATERIAL_VOID, 0, 1) >= 0);
    }
    ASSERT_EQ(alea_prepare_query_acceleration(sys), 0);
    const alea_universe_t* universe = alea_get_universe(sys, 1);
    ASSERT_NOT_NULL(universe);
    ASSERT(universe->recovery_surfaces_exact);
    ASSERT_EQ(universe->recovery_surface_count, 2);
    alea_raycast_result_t raw;
    alea_raycast_result_init(&raw);
    raw.transport_fill_filter = true;
    alea_ray_t ray;
    ASSERT_EQ(alea_ray_init(&ray, -1, 0, 0, 1, 0, 0), 0);
    ASSERT_EQ(alea_raycast_validation_breakpoints_reuse_nocache(sys, &ray, 11, 12, 3, &raw), 0);
    ASSERT_EQ(raw.hits.count, 3);
    int aliases = 0;
    for (size_t i = 0; i < raw.hits.count; ++i) {
        ASSERT_EQ(raw.hits.data[i].t, 11);
        ASSERT(raw.hits.data[i].surface_id == 100 || raw.hits.data[i].surface_id == 101);
        aliases += raw.hits.data[i].surface_id == 101;
    }
    ASSERT_EQ(aliases, 1);
    ASSERT(alea_raycast_validation_breakpoints_reuse_nocache(sys, &ray, 11, 12, 2, &raw) < 0);

    /* An edited leaf can carry a surface ID whose registered primitive differs.
     * The cell tree then supplies a crossing absent from the indexed cards. */
    alea_node_id_t leaf = alea_add_primitive_node(sys, sys->surfaces.data[offset].primitive_id, -1, 0, 100);
    ASSERT(leaf != ALEA_NODE_ID_INVALID);
    ASSERT(alea_add_cell(sys, 5, leaf, ALEA_MATERIAL_VOID, 0, 1) >= 0);
    ASSERT_EQ(alea_prepare_query_acceleration(sys), 0);
    universe = alea_get_universe(sys, 1);
    ASSERT_NOT_NULL(universe);
    ASSERT_NOT_NULL(universe->recovery_surfaces);
    ASSERT(!universe->recovery_surfaces_exact);
    ASSERT_EQ(alea_raycast_validation_breakpoints_reuse_nocache(sys, &ray, 12, 13, 3, &raw), 0);
    ASSERT_EQ(raw.hits.count, 1);
    ASSERT_EQ(raw.hits.data[0].t, 13);
    ASSERT_EQ(raw.hits.data[0].surface_id, 100);
    alea_raycast_result_free(&raw);
    alea_destroy(sys);
}

TEST(binned_recovery_matches_sorted_roots_for_arbitrary_window_order) {
    recovery_root_t input[512], sorted[512];
    for (unsigned i=0;i<512;++i) {
        /* Crowded, sparse and empty bins, with exact ties in both scan paths. */
        double t = i%2 ? 0.5+(i%4)*0.125 : (i*37)%128;
        input[i]=recovery_root_make(t,i/4,i%4);
    }
    memcpy(sorted,input,sizeof(input));
    qsort(sorted,512,sizeof(*sorted),recovery_root_compare);
    recovery_context_entry_t reference={.roots=sorted,.root_count=512};
    size_t bytes=0;
    ASSERT(recovery_bin_roots(input,512,0,128,sizeof(input),&bytes)==NULL);
    recovery_context_entry_t binned={.roots_binned=true,.root_count=512,.roots_hi=128};
    binned.roots=recovery_bin_roots(input,512,0,128,RECOVERY_ROOT_BUDGET,&bytes);
    ASSERT_NOT_NULL(binned.roots);
    ASSERT_EQ(bytes,sizeof(input)+(RECOVERY_ROOT_BINS+1)*sizeof(uint32_t));
    ASSERT(recovery_root_first(&binned,127,128,126)==NULL);
    ASSERT_EQ(binned.sorted_bins,0);
    for (unsigned i=0;i<1200;++i) {
        double lo = i%3 ? (i*31)%128 : recovery_root_distance(&input[i%512]);
        double hi = lo+(i%5);
        double current = i%2 ? nextafter(nextafter(lo,-INFINITY),-INFINITY) : lo;
        const recovery_root_t* a = recovery_root_first(&reference,lo,hi,current);
        const recovery_root_t* b = recovery_root_first(&binned,lo,hi,current);
        ASSERT((a==NULL)==(b==NULL));
        if (a) {
            ASSERT_EQ(recovery_root_distance(a),recovery_root_distance(b));
            ASSERT_EQ(a->order,b->order);
        }
    }
    ASSERT(binned.sorted_bins & UINT64_C(1));
    ASSERT_EQ(recovery_root_bin(5,5,5),0);
    ASSERT_EQ(recovery_root_bin(DBL_MAX,0,DBL_MAX),RECOVERY_ROOT_BINS-1);
    free(binned.roots);
}

TEST(quadric_packets_preserve_scalar_roots_at_numeric_edges) {
    const alea_quadric_data_t q[RAY_QUADRIC_LANES] = {
        {.coeffs={0,0,0,0,0,0,1,0,0,-3}},                       /* linear */
        {.coeffs={0,0,0,0,0,0,0,0,0,1}},                        /* constant */
        {.coeffs={0x1p-700,0,0,0,0,0,1,0,0,-1}},                /* nearly linear */
        {.coeffs={1,0,0,0,0,0,-2,0,0,1}},                       /* tangent */
        {.coeffs={1,0,0,0,0,0,-2,0,0,1+DBL_EPSILON}},           /* tolerance */
        {.coeffs={DBL_TRUE_MIN,0,0,0,0,0,-2*DBL_TRUE_MIN,0,0,1./32}},
        {.coeffs={1,0,0,0,0,0,0,0,0,-0.0}},                     /* double zero */
        {.coeffs={2,-3,5,-7,11,-13,17,-19,23,-29}}
    };
    const alea_ray_t rays[] = {
        {.dx=1},
        {.ox=-512,.oy=34,.oz=-8,.dx=0.36,.dy=0.48,.dz=-0.8},
        {.ox=1e15,.oy=-1e15,.oz=2,.dy=-1}
    };
    ray_quadric_packet_t packet;
    for (unsigned i=0; i<RAY_QUADRIC_LANES; i++)
        for (unsigned k=0; k<10; k++)
            packet.coeffs[k*RAY_QUADRIC_LANES+i] = q[i].coeffs[k];
    for (size_t r=0; r<sizeof(rays)/sizeof(rays[0]); r++) {
        for (unsigned active=1; active<=RAY_QUADRIC_LANES; active++) {
            int hits[RAY_QUADRIC_LANES];
            double roots[RAY_QUADRIC_LANES][2];
            for (unsigned i=0; i<RAY_QUADRIC_LANES; i++) hits[i]=-1;
            ray_intersect_quadric_packet(&rays[r],&packet,active,hits,roots);
            for (unsigned i=0; i<active; i++) {
                double expected[2];
                int n=ray_intersect_quadric(&rays[r],&q[i],expected);
                ASSERT_EQ(hits[i],n);
                ASSERT_EQ(memcmp(roots[i],expected,n*sizeof(double)),0);
            }
            if (active<RAY_QUADRIC_LANES) ASSERT_EQ(hits[active],-1);
        }
    }
}

TEST_MAIN()
