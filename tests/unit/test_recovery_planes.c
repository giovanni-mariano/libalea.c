// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
// SPDX-License-Identifier: MPL-2.0
#include "alea_test.h"
#include "alea.h"
#include "alea_mcnp.h"
#include "raycast/recovery_planes.h"

TEST(recovery_plane_masks_preserve_exact_kernel_roots_and_endpoints) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);
    uint32_t cards[272];
    const double normals[][3] = {{1,0,0},{0,1,0},{0,0,1},{1,2,-3}};
    for (int i=0;i<256;++i) {
        const double* n=normals[i%4];
        ASSERT(alea_plane_surface(sys, i+1, n[0], n[1], n[2], -(i/4))>=0);
        /* Deliberately reverse source order relative to sorted offsets. */
        cards[255-i]=(uint32_t)i;
    }
    for (int i=256;i<272;++i) {
        ASSERT(alea_sphere_surface(sys,i+1,0,0,0,i)>0);
        cards[i]=(uint32_t)i;
    }
    recovery_plane_index_t* index=recovery_plane_build(sys,cards,272,1024*1024);
    ASSERT_NOT_NULL(index);
    ASSERT_EQ(index->groups,4);
    ASSERT_EQ(index->members,256);
    uint64_t mask[RECOVERY_PLANE_MASK_WORDS];
    uint64_t seed=612345;
    size_t rejected=0;
    for (unsigned k=0;k<1000;++k) {
        double v[6];
        for (unsigned j=0;j<6;++j) {
            seed=seed*UINT64_C(6364136223846793005)+1;
            v[j]=(double)(int32_t)(seed>>32)/2147483648.0;
        }
        alea_ray_t ray;
        double scale=k%5==0 ? 1e150 : k%7==0 ? 1e16 : 100;
        if (k%17==0) {v[3]=0;v[4]=1;v[5]=0;}
        alea_ray_init(&ray,scale*v[0],scale*v[1],scale*v[2],v[3],v[4],v[5]);
        /* Every sampled plane root is also tested as a closed point window,
         * where an off-by-one binary search would lose an endpoint. */
        const alea_plane_data_t* endpoint=alea_primitive_payload_const(sys,
            sys->surfaces.data[cards[k%256]].primitive_id);
        double t;
        int found=ray_intersect_plane(&ray,endpoint,&t,NULL,NULL,NULL);
        double lo=(found && k%2) ? t : 0;
        double hi=(found && k%2) ? t : 20;
        recovery_plane_candidates(index,&ray,lo,hi,mask);
        for (size_t i=0;i<272;++i) {
            bool active=(mask[i/64]>>(i%64))&1;
            if (i>=256) { ASSERT(active); continue; }
            const alea_plane_data_t* p=alea_primitive_payload_const(sys,
                sys->surfaces.data[cards[i]].primitive_id);
            double root;
            int n=ray_intersect_plane(&ray,p,&root,NULL,NULL,NULL);
            ASSERT_EQ(active,n && root>=lo && root<=hi);
            rejected+=!active;
        }
    }
    ASSERT(rejected>100000);
    /* Overflow uncertainty and nonfinite windows retain every card. */
    alea_ray_t huge={.ox=DBL_MAX,.oy=DBL_MAX,.oz=DBL_MAX,.dx=1,.dy=1,.dz=1};
    recovery_plane_candidates(index,&huge,0,10,mask);
    for (unsigned i=0;i<272;++i) ASSERT(mask[i/64]&(UINT64_C(1)<<(i%64)));
    recovery_plane_candidates(index,&huge,0,INFINITY,mask);
    for (unsigned i=0;i<272;++i) ASSERT(mask[i/64]&(UINT64_C(1)<<(i%64)));
    ASSERT(recovery_plane_build(sys,cards,272,0)==NULL);
    alea_alloc_fail_after(0);
    ASSERT(recovery_plane_build(sys,cards,272,1024*1024)==NULL);
    alea_alloc_fail_after(1);
    ASSERT(recovery_plane_build(sys,cards,272,1024*1024)==NULL);
    alea_alloc_fail_after(-1);
    ASSERT(recovery_plane_build(sys,cards,RECOVERY_PLANE_MAX_CARDS+1,1024*1024)==NULL);
    alea_interrupt();
    ASSERT(recovery_plane_build(sys,cards,272,1024*1024)==NULL);
    alea_clear_interrupt();
    free(index);
    alea_destroy(sys);
}

TEST(recovery_plane_cache_preserves_alias_ties_and_rebuild_ownership) {
    mcnp_model_t* model=mcnp_load("tests/data/mcnp_deep_translated_fills.mcnp");
    ASSERT_NOT_NULL(model);
    alea_system_t* sys=model->sys;
    /* Add distant parallel planes and aliases to the terminal universe. */
    for (int i=0;i<160;++i) {
        int surface=alea_plane_surface(sys,1000+i,1,0,0,-(double)(i%80));
        ASSERT(surface>=0);
        ASSERT(alea_add_cell(sys,1000+i,alea_halfspace(sys,surface,-1),0,0,8)>=0);
    }
    ASSERT_EQ(alea_prepare_query_acceleration(sys),0);
    alea_universe_t* u=&sys->universes.data[sys->universes.count-1];
    for (size_t i=0;i<sys->universes.count;++i)
        if (sys->universes.data[i].universe_id==8) u=&sys->universes.data[i];
    ASSERT_NOT_NULL(u->recovery_planes);
    for (int direction=-1;direction<=1;direction+=2) {
        alea_ray_t ray;
        alea_ray_init(&ray,direction>0 ? 1 : 10,0,0,direction,0,0);
        alea_raycast_result_t indexed,reference;
        alea_raycast_result_init(&indexed);alea_raycast_result_init(&reference);
        indexed.breakpoint_nearest_only=reference.breakpoint_nearest_only=1;
        ASSERT_EQ(alea_raycast_validation_breakpoints_reuse_nocache(sys,&ray,0,12,1,&indexed),0);
        recovery_plane_index_t* saved=u->recovery_planes;u->recovery_planes=NULL;
        int rc=alea_raycast_validation_breakpoints_reuse_nocache(sys,&ray,0,12,1,&reference);
        u->recovery_planes=saved;
        ASSERT_EQ(rc,0);
        ASSERT_EQ(indexed.hits.count,reference.hits.count);
        ASSERT_EQ(indexed.hits.count,1);
        ASSERT_EQ(indexed.hits.data[0].t,reference.hits.data[0].t);
        ASSERT_EQ(indexed.hits.data[0].surface_id,reference.hits.data[0].surface_id);
        ASSERT_EQ(indexed.hits.data[0].primitive_id,reference.hits.data[0].primitive_id);
        ASSERT_EQ(indexed.hits.data[0].nx,reference.hits.data[0].nx);
        ASSERT(indexed.surfaces_tested<reference.surfaces_tested);
        alea_raycast_result_free(&indexed);alea_raycast_result_free(&reference);
    }
    alea_system_t* clone=alea_clone(sys);
    ASSERT_NOT_NULL(clone);
    ASSERT_EQ(alea_prepare_query_acceleration(clone),0);
    for (size_t i=0;i<clone->universes.count;++i)
        if (clone->universes.data[i].universe_id==8) {
            ASSERT_NOT_NULL(clone->universes.data[i].recovery_planes);
            ASSERT(clone->universes.data[i].recovery_planes!=u->recovery_planes);
        }
    ASSERT_EQ(alea_build_cell_surface_index(sys),0);
    ASSERT(u->recovery_planes==NULL);
    ASSERT_EQ(alea_prepare_query_acceleration(sys),0);
    ASSERT_NOT_NULL(u->recovery_planes);
    int new_surface=alea_plane_surface(sys,9999,1,0,0,0.5);
    ASSERT(new_surface>=0);
    ASSERT(alea_add_cell(sys,9999,alea_halfspace(sys,new_surface,-1),0,0,8)>=0);
    ASSERT_EQ(alea_prepare_query_acceleration(sys),0);
    alea_ray_t edited_ray;
    alea_ray_init(&edited_ray,1,0,0,1,0,0);
    alea_raycast_result_t edited;
    alea_raycast_result_init(&edited);
    edited.breakpoint_nearest_only=1;
    ASSERT_EQ(alea_raycast_validation_breakpoints_reuse_nocache(
        sys,&edited_ray,0,12,1,&edited),0);
    ASSERT_EQ(edited.hits.count,1);
    ASSERT_EQ(edited.hits.data[0].surface_id,9999);
    ASSERT_EQ(edited.hits.data[0].t,0.5);
    alea_raycast_result_free(&edited);
    alea_destroy(clone);
    mcnp_model_destroy(model);
}
TEST(recovery_plane_blocks_keep_closed_window_hits_and_fit_preparation_budget) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);
    enum { count = 193 };
    uint32_t cards[count];
    for (unsigned i = 0; i < count; i++) {
        if (i == count - 1) {
            ASSERT(alea_sphere_surface(sys, i+1, 0, 0, 0, 3) >= 0);
        } else {
            double b = i < 16 ? 0 : i * 0.001;
            double c = i < 16 ? 0 : (int)(i % 7) * 0.05;
            ASSERT(alea_plane_surface(sys, i+1, 1, b, c, 85.0-i) >= 0);
        }
        cards[count-1-i] = i;
    }
    recovery_plane_index_t* p = recovery_plane_build(sys, cards, count, 1024*1024);
    ASSERT_NOT_NULL(p);
    ASSERT(p->blocks > 1);
    ASSERT_EQ(p->block_members, count-17);
    size_t base_bytes = (const char*)recovery_plane_block_ordinals(p) - (const char*)p;
    recovery_plane_index_t* small = recovery_plane_build(sys, cards, count, base_bytes);
    ASSERT_NOT_NULL(small);
    ASSERT_EQ(small->bytes, base_bytes);
    ASSERT_EQ(small->blocks, 0);
    free(small);

    alea_alloc_fail_after(1); /* Expanded index allocation fails after scratch. */
    small = recovery_plane_build(sys, cards, count, 1024*1024);
    alea_alloc_fail_after(-1);
    ASSERT_NOT_NULL(small);
    ASSERT_EQ(small->blocks, 0);
    ASSERT_EQ(small->bytes, base_bytes);
    free(small);

    uint64_t mask[(count+63)/64], original[(count+63)/64];
    size_t rejected = 0;
    for (unsigned k = 0; k < 600; k++) {
        alea_ray_t ray;
        double ox = k % 3 == 0 ? 1e15 : (double)(k % 151) - 75;
        ASSERT_EQ(alea_ray_init(&ray, ox, -4, 0.01*k,
                                k % 3 == 0 ? -1 : 1, (int)(k%9)-4, 0.5), 0);
        const alea_plane_data_t* target = alea_primitive_payload_const(sys,
            sys->surfaces.data[cards[1+k%(count-1)]].primitive_id);
        double t;
        double lo = 0, hi = 3;
        if (ray_intersect_plane(&ray, target, &t, NULL, NULL, NULL) && t >= 0) {
            if (k % 3 == 0) lo = hi = t;
            if (k % 3 == 1) hi = nextafter(t, INFINITY);
            if (k % 3 == 2) { lo = nextafter(t, INFINITY); hi = lo + 3; }
        }
        recovery_plane_candidates(p, &ray, lo, hi, original);
        recovery_plane_root_candidates(p, &ray, lo, hi, mask);
        ASSERT(mask[0] & 1); /* The sphere always remains active. */
        for (unsigned i = 1; i < count; i++) {
            const alea_plane_data_t* plane = alea_primitive_payload_const(sys,
                sys->surfaces.data[cards[i]].primitive_id);
            uint64_t bit = UINT64_C(1) << (i % 64);
            if (ray_intersect_plane(&ray, plane, &t, NULL, NULL, NULL) && t >= lo && t <= hi)
                ASSERT(mask[i/64] & bit);
            rejected += (original[i/64] & bit) && !(mask[i/64] & bit);
        }
    }
    ASSERT(rejected > 1000);
    alea_ray_t huge = {.ox=DBL_MAX, .oy=DBL_MAX, .dx=1};
    recovery_plane_candidates(p, &huge, 0, 1, original);
    recovery_plane_root_candidates(p, &huge, 0, 1, mask);
    ASSERT_EQ(memcmp(mask, original, sizeof(mask)), 0);
    free(p);
    alea_destroy(sys);
}

TEST_MAIN()
