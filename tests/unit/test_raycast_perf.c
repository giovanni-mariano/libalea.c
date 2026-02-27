// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

#define _POSIX_C_SOURCE 199309L

/*
 * test_raycast_perf.c - Raycast performance regression benchmark
 *
 * Builds representative geometries and traces many rays, reporting
 * timing for each workload.  All "tests" pass unconditionally — this
 * file is for performance tracking, not correctness.
 *
 * Run standalone:
 *   bin/tests/unit/test_raycast_perf
 *
 * Compare across commits:
 *   git stash && make test_raycast_perf && bin/tests/unit/test_raycast_perf
 *   git stash pop && make test_raycast_perf && bin/tests/unit/test_raycast_perf
 */

#include "alea_test.h"
#include "alea.h"
#include "raycast/raycast.h"
#include "raycast/ray_intersect.h"
#include "raycast/bvh.h"
#include "core/alea_system.h"

#include <math.h>
#include <time.h>
#include <stdio.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* =========================================================================
 * Timer helpers
 * ========================================================================= */

static double now_sec(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + ts.tv_nsec * 1e-9;
}

#define BENCH_START() double _t0 = now_sec()
#define BENCH_END(label, n_ops) do { \
    double _dt = now_sec() - _t0; \
    double _us = (_dt / (n_ops)) * 1e6; \
    printf("%.3f ms total, %.2f us/op (%d ops)  ", \
           _dt * 1e3, _us, (int)(n_ops)); \
} while(0)

/* =========================================================================
 * Geometry builders
 * ========================================================================= */

/**
 * Build a single sphere cell (simplest possible geometry).
 */
static alea_system_t* build_single_sphere(void) {
    alea_system_t* sys = alea_create();
    int si = alea_sphere_surface(sys, 1, 0, 0, 0, 5.0);
    alea_node_id_t interior = alea_halfspace(sys, si, -1);
    alea_add_cell(sys, 1, interior, 1, -2.7, 0);
    return sys;
}

/**
 * Build N concentric spherical shells.
 * Each shell is a cell: sphere[i] \ sphere[i-1].
 */
static alea_system_t* build_concentric_shells(int n_shells) {
    alea_system_t* sys = alea_create();

    alea_node_id_t prev_inner = 0;
    for (int i = 0; i < n_shells; i++) {
        double r = 1.0 + (double)i * 0.5;
        int si = alea_sphere_surface(sys, i + 1, 0, 0, 0, r);
        alea_node_id_t inside = alea_halfspace(sys, si, -1);

        alea_node_id_t region;
        if (i == 0) {
            region = inside;
        } else {
            alea_node_id_t outside_prev = alea_halfspace(sys,
                alea_sphere_surface(sys, 100 + i, 0, 0, 0, 1.0 + (i - 1) * 0.5), +1);
            region = alea_intersection(sys, inside, outside_prev);
        }
        alea_add_cell(sys, i + 1, region, (i % 3) + 1, -2.7, 0);
        prev_inner = inside;
    }
    (void)prev_inner;
    return sys;
}

/**
 * Build a grid of N x N x N box cells (tests BVH with many surfaces).
 * Each cell is a box [i, i+1] x [j, j+1] x [k, k+1].
 */
static alea_system_t* build_box_grid(int n) {
    alea_system_t* sys = alea_create();
    int surf_id = 1;
    int cell_id = 1;

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            for (int k = 0; k < n; k++) {
                double x0 = (double)i, x1 = (double)(i + 1);
                double y0 = (double)j, y1 = (double)(j + 1);
                double z0 = (double)k, z1 = (double)(k + 1);

                int si = alea_box_surface(sys, surf_id++,
                                          x0, x1, y0, y1, z0, z1);
                alea_node_id_t inside = alea_halfspace(sys, si, -1);
                alea_add_cell(sys, cell_id++, inside,
                              ((i + j + k) % 3) + 1, -2.7, 0);
            }
        }
    }
    return sys;
}

/**
 * Build a mixed-primitive geometry with spheres, cylinders, cones, planes.
 */
static alea_system_t* build_mixed_primitives(int n_spheres) {
    alea_system_t* sys = alea_create();
    int surf_id = 1;
    int cell_id = 1;

    /* Scatter spheres */
    for (int i = 0; i < n_spheres; i++) {
        double angle = 2.0 * M_PI * i / n_spheres;
        double cx = 10.0 * cos(angle);
        double cy = 10.0 * sin(angle);
        int si = alea_sphere_surface(sys, surf_id++, cx, cy, 0, 1.5);
        alea_node_id_t inside = alea_halfspace(sys, si, -1);
        alea_add_cell(sys, cell_id++, inside, 1, -2.7, 0);
    }

    /* Central cylinder */
    int ci = alea_cylinder_z_surface(sys, surf_id++, 0, 0, 3.0);
    alea_node_id_t cyl_in = alea_halfspace(sys, ci, -1);

    /* Cap planes */
    int pi_bot = alea_plane_surface(sys, surf_id++, 0, 0, 1, 5.0);
    int pi_top = alea_plane_surface(sys, surf_id++, 0, 0, 1, -5.0);
    alea_node_id_t below_top = alea_halfspace(sys, pi_top, -1);
    alea_node_id_t above_bot = alea_halfspace(sys, pi_bot, -1);

    alea_node_id_t capped = alea_intersection(sys,
        alea_intersection(sys, cyl_in, below_top), above_bot);
    alea_add_cell(sys, cell_id++, capped, 2, -8.0, 0);

    return sys;
}

/* =========================================================================
 * Deterministic ray generator (matches volume estimation pattern)
 * ========================================================================= */

static uint32_t bench_lcg(uint32_t* state) {
    *state = *state * 1664525u + 1013904223u;
    return *state;
}

static double bench_rand(uint32_t* state) {
    return (double)bench_lcg(state) / 4294967296.0;
}

static void bench_random_ray(uint32_t* rng, double extent,
                             double* ox, double* oy, double* oz,
                             double* dx, double* dy, double* dz) {
    /* Random origin in [-extent, extent]^3 */
    *ox = (bench_rand(rng) * 2.0 - 1.0) * extent;
    *oy = (bench_rand(rng) * 2.0 - 1.0) * extent;
    *oz = (bench_rand(rng) * 2.0 - 1.0) * extent;

    /* Random direction (uniform on sphere) */
    double phi = 2.0 * M_PI * bench_rand(rng);
    double cos_theta = 2.0 * bench_rand(rng) - 1.0;
    double sin_theta = sqrt(1.0 - cos_theta * cos_theta);
    *dx = sin_theta * cos(phi);
    *dy = sin_theta * sin(phi);
    *dz = cos_theta;
}

/* =========================================================================
 * Benchmarks
 * ========================================================================= */

TEST(perf_ray_init) {
    int N = 100000;
    alea_ray_t ray;
    uint32_t rng = 123;

    BENCH_START();
    for (int i = 0; i < N; i++) {
        double dx = bench_rand(&rng) * 2 - 1;
        double dy = bench_rand(&rng) * 2 - 1;
        double dz = bench_rand(&rng) * 2 - 1;
        alea_ray_init(&ray, 0, 0, 0, dx, dy, dz);
    }
    BENCH_END("ray_init", N);
}

TEST(perf_ray_init_normalized) {
    int N = 100000;
    alea_ray_t ray;
    uint32_t rng = 123;

    BENCH_START();
    for (int i = 0; i < N; i++) {
        double dx = bench_rand(&rng) * 2 - 1;
        double dy = bench_rand(&rng) * 2 - 1;
        double dz = bench_rand(&rng) * 2 - 1;
        double len = sqrt(dx*dx + dy*dy + dz*dz);
        if (len > 1e-10) { dx /= len; dy /= len; dz /= len; }
        alea_ray_init_normalized(&ray, 0, 0, 0, dx, dy, dz);
    }
    BENCH_END("ray_init_normalized (incl. manual normalize)", N);
}

TEST(perf_sphere_intersect) {
    int N = 1000000;
    alea_ray_t ray;
    alea_sphere_data_t sphere = {0, 0, 0, 5.0};
    double t[2];
    uint32_t rng = 456;
    int total_hits = 0;

    BENCH_START();
    for (int i = 0; i < N; i++) {
        double dx = bench_rand(&rng) * 2 - 1;
        double dy = bench_rand(&rng) * 2 - 1;
        double dz = bench_rand(&rng) * 2 - 1;
        alea_ray_init(&ray, -10, 0, 0, dx, dy, dz);
        total_hits += ray_intersect_sphere(&ray, &sphere, t);
    }
    BENCH_END("sphere intersect", N);
    printf("[%d hits]  ", total_hits);
}

TEST(perf_torus_intersect) {
    int N = 100000;
    alea_ray_t ray;
    alea_torus_data_t torus = {
        .axis = ALEA_AXIS_Z, .center_x = 0, .center_y = 0, .center_z = 0,
        .major_radius = 10.0, .minor_radius = 2.0, .axial_semiwidth_B = 0
    };
    double t[4];
    uint32_t rng = 789;
    int total_hits = 0;

    BENCH_START();
    for (int i = 0; i < N; i++) {
        double dx = bench_rand(&rng) * 2 - 1;
        double dy = bench_rand(&rng) * 2 - 1;
        double dz = bench_rand(&rng) * 2 - 1;
        alea_ray_init(&ray, -20, 0, 0, dx, dy, dz);
        total_hits += ray_intersect_torus(&ray, &torus, t);
    }
    BENCH_END("torus intersect (quartic)", N);
    printf("[%d hits]  ", total_hits);
}

TEST(perf_single_sphere_raycast) {
    int N = 10000;
    alea_system_t* sys = build_single_sphere();
    alea_raycast_ensure_caches(sys);

    alea_raycast_result_t result;
    alea_raycast_result_init(&result);
    uint32_t rng = 42;

    BENCH_START();
    for (int i = 0; i < N; i++) {
        double ox, oy, oz, dx, dy, dz;
        bench_random_ray(&rng, 10.0, &ox, &oy, &oz, &dx, &dy, &dz);
        alea_raycast_result_free(&result);
        alea_raycast(sys, ox, oy, oz, dx, dy, dz, 100, &result);
    }
    BENCH_END("single sphere full raycast", N);

    alea_raycast_result_free(&result);
    alea_destroy(sys);
}

TEST(perf_concentric_shells_20) {
    int N = 10000;
    alea_system_t* sys = build_concentric_shells(20);
    alea_raycast_ensure_caches(sys);

    alea_raycast_result_t result;
    alea_raycast_result_init(&result);
    uint32_t rng = 42;
    size_t total_segs = 0;

    BENCH_START();
    for (int i = 0; i < N; i++) {
        double ox, oy, oz, dx, dy, dz;
        bench_random_ray(&rng, 15.0, &ox, &oy, &oz, &dx, &dy, &dz);
        alea_raycast_result_free(&result);
        alea_raycast(sys, ox, oy, oz, dx, dy, dz, 100, &result);
        total_segs += result.segment_count;
    }
    BENCH_END("20 concentric shells", N);
    printf("[avg %.1f segs]  ", (double)total_segs / N);

    alea_raycast_result_free(&result);
    alea_destroy(sys);
}

TEST(perf_box_grid_5x5x5) {
    int N = 5000;
    alea_system_t* sys = build_box_grid(5);
    alea_raycast_ensure_caches(sys);

    alea_raycast_result_t result;
    alea_raycast_result_init(&result);
    uint32_t rng = 42;
    size_t total_hits = 0;

    BENCH_START();
    for (int i = 0; i < N; i++) {
        double ox, oy, oz, dx, dy, dz;
        bench_random_ray(&rng, 8.0, &ox, &oy, &oz, &dx, &dy, &dz);
        alea_raycast_result_free(&result);
        alea_raycast(sys, ox, oy, oz, dx, dy, dz, 100, &result);
        total_hits += result.hit_count;
    }
    BENCH_END("5x5x5 box grid (125 cells, BVH)", N);
    printf("[avg %.1f hits]  ", (double)total_hits / N);

    alea_raycast_result_free(&result);
    alea_destroy(sys);
}

TEST(perf_mixed_primitives_32) {
    int N = 5000;
    alea_system_t* sys = build_mixed_primitives(32);
    alea_raycast_ensure_caches(sys);

    alea_raycast_result_t result;
    alea_raycast_result_init(&result);
    uint32_t rng = 42;

    BENCH_START();
    for (int i = 0; i < N; i++) {
        double ox, oy, oz, dx, dy, dz;
        bench_random_ray(&rng, 15.0, &ox, &oy, &oz, &dx, &dy, &dz);
        alea_raycast_result_free(&result);
        alea_raycast(sys, ox, oy, oz, dx, dy, dz, 100, &result);
    }
    BENCH_END("mixed primitives (32 spheres + capped cyl)", N);

    alea_raycast_result_free(&result);
    alea_destroy(sys);
}

TEST(perf_surfaces_only_nocache) {
    int N = 10000;
    alea_system_t* sys = build_mixed_primitives(32);
    alea_raycast_ensure_caches(sys);

    alea_raycast_result_t result;
    alea_raycast_result_init(&result);
    alea_raycast_result_reserve(&result, 64, 32);
    uint32_t rng = 42;

    BENCH_START();
    for (int i = 0; i < N; i++) {
        double ox, oy, oz, dx, dy, dz;
        bench_random_ray(&rng, 15.0, &ox, &oy, &oz, &dx, &dy, &dz);
        alea_ray_t ray;
        alea_ray_init(&ray, ox, oy, oz, dx, dy, dz);
        alea_raycast_surfaces_nocache(sys, &ray, 0, 100, &result);
    }
    BENCH_END("surfaces_nocache (no segments)", N);

    alea_raycast_result_free(&result);
    alea_destroy(sys);
}

TEST(perf_bvh_build) {
    int N = 100;
    alea_system_t* sys = build_box_grid(5);

    BENCH_START();
    for (int i = 0; i < N; i++) {
        alea_bvh_t* bvh = alea_bvh_build(sys);
        alea_bvh_free(bvh);
    }
    BENCH_END("BVH build (125 cells)", N);

    alea_destroy(sys);
}

TEST(perf_occlusion_test) {
    int N = 10000;
    alea_system_t* sys = build_concentric_shells(10);
    alea_raycast_ensure_caches(sys);
    uint32_t rng = 42;
    int occluded = 0;

    BENCH_START();
    for (int i = 0; i < N; i++) {
        double ox, oy, oz, dx, dy, dz;
        bench_random_ray(&rng, 10.0, &ox, &oy, &oz, &dx, &dy, &dz);
        occluded += alea_ray_is_occluded(sys, ox, oy, oz, dx, dy, dz, 50);
    }
    BENCH_END("occlusion test (10 shells)", N);
    printf("[%d/%d occluded]  ", occluded, N);

    alea_destroy(sys);
}

TEST(perf_volume_estimation) {
    int N_RAYS = 5000;
    alea_system_t* sys = build_concentric_shells(10);

    size_t n_cells = alea_cell_count(sys);
    double* volumes = calloc(n_cells, sizeof(double));
    double* errors = calloc(n_cells, sizeof(double));

    BENCH_START();
    alea_estimate_cell_volumes(sys, 0, 0, 0, 10.0, N_RAYS, volumes, errors);
    BENCH_END("volume estimation (10 shells)", N_RAYS);

    /* Sanity: at least some cells got volume */
    int nonzero = 0;
    for (size_t i = 0; i < n_cells; i++) {
        if (volumes[i] > 0) nonzero++;
    }
    printf("[%d/%zu cells with volume]  ", nonzero, n_cells);

    free(volumes);
    free(errors);
    alea_destroy(sys);
}

TEST_MAIN()
