// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/*
 * test_bvh.c - Unit tests for BVH acceleration structure
 */

#include "alea_test.h"
#include "alea.h"
#include "core/alea_system.h"
#include "raycast/raycast.h"
#include "raycast/bvh.h"

/* ------------------------------------------------------------------------- */
/* BVH Construction Tests                                                     */
/* ------------------------------------------------------------------------- */

TEST(bvh_empty) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);

    /* No surfaces - BVH should be NULL */
    alea_bvh_t* bvh = alea_bvh_build(sys);
    ASSERT_NULL(bvh);

    alea_destroy(sys);
}

TEST(bvh_single_surface) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);

    int surf_idx = alea_sphere_surface(sys, 1, 0, 0, 0, 5.0);
    ASSERT(surf_idx >= 0);

    alea_bvh_t* bvh = alea_bvh_build(sys);
    ASSERT_NOT_NULL(bvh);

    size_t node_count, leaf_count, max_depth;
    alea_bvh_stats(bvh, &node_count, &leaf_count, &max_depth);

    ASSERT_EQ(node_count, 1);
    ASSERT_EQ(leaf_count, 1);
    ASSERT_EQ(max_depth, 0);

    alea_bvh_free(bvh);
    alea_destroy(sys);
}

TEST(bvh_multiple_surfaces) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);

    /* Create 10 spheres along X axis */
    for (int i = 0; i < 10; i++) {
        double x = i * 10.0;
        int surf_idx = alea_sphere_surface(sys, i + 1, x, 0, 0, 2.0);
        ASSERT(surf_idx >= 0);
    }

    alea_bvh_t* bvh = alea_bvh_build(sys);
    ASSERT_NOT_NULL(bvh);

    size_t node_count, leaf_count, max_depth;
    alea_bvh_stats(bvh, &node_count, &leaf_count, &max_depth);

    /* With 10 surfaces and leaf threshold of 4, we should have a tree */
    ASSERT(node_count >= 3);
    ASSERT(leaf_count >= 2);
    ASSERT(max_depth >= 1);

    alea_bvh_free(bvh);
    alea_destroy(sys);
}

/* ------------------------------------------------------------------------- */
/* BVH Traversal Tests                                                        */
/* ------------------------------------------------------------------------- */

typedef struct {
    int count;
    uint32_t indices[100];
} traverse_result_t;

static void count_callback(uint32_t surface_idx, void* userdata) {
    traverse_result_t* res = (traverse_result_t*)userdata;
    if (res->count < 100) {
        res->indices[res->count] = surface_idx;
    }
    res->count++;
}

TEST(bvh_traversal_hit) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);

    /* Create spheres along X axis */
    for (int i = 0; i < 20; i++) {
        double x = i * 5.0;
        int surf_idx = alea_sphere_surface(sys, i + 1, x, 0, 0, 1.0);
        ASSERT(surf_idx >= 0);
    }

    alea_bvh_t* bvh = alea_bvh_build(sys);
    ASSERT_NOT_NULL(bvh);

    /* Ray along X axis should hit sphere 0 first */
    alea_ray_t ray;
    alea_ray_init(&ray, -10, 0, 0, 1, 0, 0);

    traverse_result_t result = {0};
    alea_bvh_traverse(bvh, &ray, 0, 1000, count_callback, &result);

    /* BVH should find at least some candidates */
    ASSERT(result.count >= 1);

    alea_bvh_free(bvh);
    alea_destroy(sys);
}

TEST(bvh_traversal_miss) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);

    /* Create spheres at positive Z */
    for (int i = 0; i < 10; i++) {
        alea_sphere_surface(sys, i + 1, i * 3.0, 0, 10, 1.0);
    }

    alea_bvh_t* bvh = alea_bvh_build(sys);
    ASSERT_NOT_NULL(bvh);

    /* Ray at negative Z shouldn't hit anything */
    alea_ray_t ray;
    alea_ray_init(&ray, 0, 0, -100, 1, 0, 0);

    traverse_result_t result = {0};
    alea_bvh_traverse(bvh, &ray, 0, 1000, count_callback, &result);

    ASSERT_EQ(result.count, 0);

    alea_bvh_free(bvh);
    alea_destroy(sys);
}

/* ------------------------------------------------------------------------- */
/* BVH Performance Tests                                                      */
/* ------------------------------------------------------------------------- */

TEST(bvh_vs_linear) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);

    /* Create a 5x5x5 grid of spheres */
    int surface_id = 1;
    for (int x = 0; x < 5; x++) {
        for (int y = 0; y < 5; y++) {
            for (int z = 0; z < 5; z++) {
                int surf_idx = alea_sphere_surface(sys, surface_id++,
                    x * 3.0, y * 3.0, z * 3.0, 1.0);
                ASSERT(surf_idx >= 0);
            }
        }
    }

    /* Force BVH build by doing a raycast */
    alea_raycast_result_t result;
    alea_raycast_result_init(&result);

    alea_ray_t ray;
    alea_ray_init(&ray, -5, 3, 3, 1, 0, 0);

    int rc = alea_raycast_surfaces(sys, &ray, 0, 100, &result);

    ASSERT_NOT_NULL(sys->surface_bvh);
    ASSERT_EQ(rc, 0);

    /* BVH should test fewer surfaces than total count */
    ASSERT(result.surfaces_tested < 125);

    alea_raycast_result_free(&result);
    alea_destroy(sys);
}

TEST(bvh_large_surface_count) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);

    /* Create 1000 spheres in a grid */
    int surface_id = 1;
    for (int i = 0; i < 1000; i++) {
        double x = (i % 10) * 5.0;
        double y = ((i / 10) % 10) * 5.0;
        double z = (i / 100) * 5.0;
        int surf_idx = alea_sphere_surface(sys, surface_id++, x, y, z, 1.0);
        ASSERT(surf_idx >= 0);
    }

    alea_bvh_t* bvh = alea_bvh_build(sys);
    ASSERT_NOT_NULL(bvh);

    size_t node_count, leaf_count, max_depth;
    alea_bvh_stats(bvh, &node_count, &leaf_count, &max_depth);

    /* Depth should be reasonable for 1000 surfaces (log2(1000) ~ 10) */
    ASSERT(max_depth <= 20);
    ASSERT(node_count > 0);

    alea_bvh_free(bvh);
    alea_destroy(sys);
}

TEST(bvh_rebuild_on_change) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);

    alea_sphere_surface(sys, 1, 0, 0, 0, 1.0);

    /* First raycast builds BVH */
    alea_raycast_result_t result;
    alea_raycast_result_init(&result);

    alea_ray_t ray;
    alea_ray_init(&ray, -5, 0, 0, 1, 0, 0);
    alea_raycast_surfaces(sys, &ray, 0, 100, &result);

    size_t initial_count = sys->surface_bvh ? sys->surface_bvh->surface_count : 0;

    /* Add more surfaces */
    alea_sphere_surface(sys, 2, 10, 0, 0, 1.0);
    sys->bvh_dirty = true;

    /* Next raycast should rebuild */
    alea_raycast_result_clear(&result);
    alea_raycast_surfaces(sys, &ray, 0, 100, &result);

    size_t new_count = sys->surface_bvh ? sys->surface_bvh->surface_count : 0;
    ASSERT(new_count > initial_count);

    alea_raycast_result_free(&result);
    alea_destroy(sys);
}

TEST_MAIN()
