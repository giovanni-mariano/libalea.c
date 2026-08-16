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
#include "alea_geo_validator.h"
#include "alea_mcnp.h"
#include "alea_raycast.h"
#include "alea_slice.h"
#include "render/render3d.h"
#include "raycast/raycast.h"
#include "raycast/ray_intersect.h"
#include "raycast/bvh.h"
#include "core/alea_system.h"
#include "util/compat.h"   /* alea_monotonic_seconds */

#include <math.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* =========================================================================
 * Timer helpers
 * ========================================================================= */

static double now_sec(void) {
    return alea_monotonic_seconds();
}

#define BENCH_START() double _t0 = now_sec()
#define BENCH_END(label, n_ops) do { \
    double _dt = now_sec() - _t0; \
    double _us = (_dt / (n_ops)) * 1e6; \
    printf("%.3f ms total, %.2f us/op (%d ops)  ", \
           _dt * 1e3, _us, (int)(n_ops)); \
} while(0)

/* Audit helpers intentionally report storage owned by the query result, not
 * allocator-implementation details.  This keeps the benchmark portable while
 * making retained-capacity growth and compact-result publication visible. */
static size_t raycast_result_retained_bytes(const alea_raycast_result_t* result) {
    return result->hits.capacity * sizeof(*result->hits.data) +
           result->segments.capacity * sizeof(*result->segments.data) +
           result->paths.capacity * sizeof(*result->paths.data) +
           result->path_entries.capacity * sizeof(*result->path_entries.data);
}

static size_t batch_published_bytes(const alea_raycast_batch_result_t* result,
                                    uint32_t fields) {
    const size_t rays = alea_raycast_batch_ray_count(result);
    const size_t segments = alea_raycast_batch_segment_count(result);
    size_t bytes = (rays + 1) * sizeof(uint64_t) + 2 * segments * sizeof(double);
    if (fields & ALEA_RAY_BATCH_MATERIAL) bytes += segments * sizeof(int32_t);
    if (fields & ALEA_RAY_BATCH_DENSITY) bytes += segments * sizeof(double);
    if (fields & ALEA_RAY_BATCH_SURFACES) bytes += 2 * segments * sizeof(int32_t);
    if (fields & ALEA_RAY_BATCH_RESOLUTION_FLAGS) bytes += segments * sizeof(uint8_t);
    return bytes;
}

static int count_selected_segment(void* context,
                                  const alea_ray_segment_t* segment) {
    (void)segment;
    (*(size_t*)context)++;
    return 0;
}

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

    int m1 = alea_add_material(sys, 1);

    alea_add_cell(sys, 1, interior, m1, -2.7, 0);
    return sys;
}

/**
 * Build N concentric spherical shells.
 * Each shell is a cell: sphere[i] \ sphere[i-1].
 */
static alea_system_t* build_concentric_shells(int n_shells) {
    alea_system_t* sys = alea_create();

    int m1 = alea_add_material(sys, 1);
    int m2 = alea_add_material(sys, 2);
    int m3 = alea_add_material(sys, 3);
    int mats[3] = {m1, m2, m3};

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
        alea_add_cell(sys, i + 1, region, mats[i % 3], -2.7, 0);
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

    int m1 = alea_add_material(sys, 1);
    int m2 = alea_add_material(sys, 2);
    int m3 = alea_add_material(sys, 3);
    int mats[3] = {m1, m2, m3};

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
                              mats[(i + j + k) % 3], -2.7, 0);
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

    int m1 = alea_add_material(sys, 1);
    int m2 = alea_add_material(sys, 2);

    /* Scatter spheres */
    for (int i = 0; i < n_spheres; i++) {
        double angle = 2.0 * M_PI * i / n_spheres;
        double cx = 10.0 * cos(angle);
        double cy = 10.0 * sin(angle);
        int si = alea_sphere_surface(sys, surf_id++, cx, cy, 0, 1.5);
        alea_node_id_t inside = alea_halfspace(sys, si, -1);
        alea_add_cell(sys, cell_id++, inside, m1, -2.7, 0);
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
    alea_add_cell(sys, cell_id++, capped, m2, -8.0, 0);

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
        total_segs += result.segments.count;
    }
    BENCH_END("20 concentric shells", N);
    printf("[avg %.1f segs]  ", (double)total_segs / N);

    alea_raycast_result_free(&result);
    alea_destroy(sys);
}

TEST(perf_first_visible_vs_full_trace_20_shells) {
    const int N = 10000;
    alea_system_t* sys = build_concentric_shells(20);
    ASSERT_EQ(alea_raycast_ensure_hier_caches(sys), 0);

    alea_raycast_result_t scratch;
    alea_raycast_result_init(&scratch);
    alea_raycast_result_reserve(&scratch, 64, 32);
    alea_ray_t ray;
    ASSERT_EQ(alea_ray_init(&ray, -20, 0.125, 0, 1, 0, 0), 0);

    {
        BENCH_START();
        for (int i = 0; i < N; i++) {
            alea_raycast_result_clear(&scratch);
            ASSERT_EQ(alea_raycast_hier_with_hits_nocache(sys, &ray, 100,
                                                           &scratch), 0);
        }
        BENCH_END("20 shells full hierarchical trace", N);
    }
    const int full_steps = scratch.step_iterations;
    const size_t full_retained = raycast_result_retained_bytes(&scratch);
    printf("[%d final steps, %zu retained bytes]  ", full_steps, full_retained);

    alea_ray_first_visible_result_t visible;
    {
        BENCH_START();
        for (int i = 0; i < N; i++) {
            ASSERT_EQ(alea_raycast_hier_first_visible_nocache(
                          sys, &ray, 0, 100, -1, 0, &scratch, &visible), 0);
            ASSERT(visible.found);
        }
        BENCH_END("20 shells first-visible early stop", N);
    }
    printf("[%d final steps, %zu retained bytes]  ", scratch.step_iterations,
           raycast_result_retained_bytes(&scratch));
    ASSERT(scratch.step_iterations < full_steps);

    alea_raycast_result_free(&scratch);
    alea_destroy(sys);
}

TEST(perf_directional_event_cache_sphere) {
    const int width = 256, height = 256, iterations = 5;
    alea_system_t* sys = build_single_sphere();
    ASSERT_NOT_NULL(sys);
    ASSERT_EQ(alea_prepare_query_acceleration(sys), 0);
    alea_slice_view_t view;
    alea_slice_view_axis(&view, 2, 0.0, -8.0, 8.0, -8.0, 8.0);

    size_t total_events = 0;
    BENCH_START();
    for (int i = 0; i < iterations; i++) {
        alea_slice_directional_event_cache_t* cache =
            alea_slice_directional_event_cache_create(sys, &view, width, height);
        ASSERT_NOT_NULL(cache);
        for (int orient = ALEA_SLICE_EDGE_RIGHT;
             orient <= ALEA_SLICE_EDGE_DOWN; orient++) {
            int lines = orient == ALEA_SLICE_EDGE_RIGHT ? height : width;
            for (int reverse = 0; reverse <= 1; reverse++) {
                for (int line = 0; line < lines; line++) {
                    const alea_ray_boundary_event_t* events;
                    size_t count;
                    ASSERT_EQ(alea_slice_directional_event_cache_line_events(
                                  cache, orient, reverse, line, &events, &count), 0);
                    total_events += count;
                }
            }
        }
        alea_slice_directional_event_cache_destroy(cache);
    }
    BENCH_END("256x256 directional event cache (U/V +/-)", iterations);
    printf("[avg %.1f events/cache]  ", (double)total_events / iterations);
    alea_destroy(sys);
}

TEST(perf_directional_event_cache_lattice) {
    const int width = 128, height = 128, iterations = 3;
    mcnp_model_t* model = mcnp_load("tests/data/mcnp_lattice_eval.mcnp");
    if (!model) SKIP("lattice fixture not found");
    alea_system_t* sys = model->sys;
    ASSERT_EQ(alea_prepare_query_acceleration(sys), 0);
    alea_slice_view_t view;
    alea_slice_view_axis(&view, 2, 0.0, -2.0, 6.0, -2.0, 6.0);

    size_t total_events = 0;
    BENCH_START();
    for (int i = 0; i < iterations; i++) {
        alea_slice_directional_event_cache_t* cache =
            alea_slice_directional_event_cache_create(sys, &view, width, height);
        ASSERT_NOT_NULL(cache);
        for (int orient = ALEA_SLICE_EDGE_RIGHT;
             orient <= ALEA_SLICE_EDGE_DOWN; orient++) {
            int lines = orient == ALEA_SLICE_EDGE_RIGHT ? height : width;
            for (int reverse = 0; reverse <= 1; reverse++) {
                for (int line = 0; line < lines; line++) {
                    const alea_ray_boundary_event_t* events;
                    size_t count;
                    ASSERT_EQ(alea_slice_directional_event_cache_line_events(
                                  cache, orient, reverse, line, &events, &count), 0);
                    total_events += count;
                }
            }
        }
        alea_slice_directional_event_cache_destroy(cache);
    }
    BENCH_END("128x128 lattice directional event cache (U/V +/-)", iterations);
    printf("[avg %.1f events/cache]  ", (double)total_events / iterations);
    mcnp_model_destroy(model);
}

TEST(perf_shared_validator_event_cache_sphere) {
    const int samples = 64;
    alea_system_t* sys = build_single_sphere();
    ASSERT_NOT_NULL(sys);
    ASSERT_EQ(alea_prepare_query_acceleration(sys), 0);
    alea_slice_view_t view;
    alea_slice_view_axis(&view, 2, 0.0, -8.0, 8.0, -8.0, 8.0);
    alea_ray_slice_validation_options_t options;
    alea_ray_slice_validation_options_init(&options);
    options.checks = ALEA_RAY_SLICE_VALIDATE_FAST_BIDIRECTIONAL;
    alea_ray_slice_validation_result_t* validation =
        alea_ray_slice_validation_result_create();
    ASSERT_NOT_NULL(validation);
    BENCH_START();
    alea_slice_directional_event_cache_t* cache =
        alea_slice_directional_event_cache_create(sys, &view, samples, samples);
    ASSERT_NOT_NULL(cache);
    ASSERT_EQ(alea_validate_ray_slice_compact_with_event_cache(
                  sys, &view, samples, &options, NULL, NULL, cache,
                  validation), 0);
    BENCH_END("64x64 sphere shared validator + event cache", 1);
    printf("[%zu intervals, cache-hit traces=0x%x]  ",
           alea_ray_slice_validation_interval_count(validation),
           alea_ray_slice_validation_reused_trace_mask(validation));
    alea_slice_directional_event_cache_destroy(cache);
    alea_ray_slice_validation_result_destroy(validation);
    alea_destroy(sys);
}

TEST(perf_shared_validator_event_cache_lattice) {
    const int samples = 64;
    mcnp_model_t* model = mcnp_load("tests/data/mcnp_lattice_eval.mcnp");
    if (!model) SKIP("lattice fixture not found");
    alea_system_t* sys = model->sys;
    ASSERT_EQ(alea_prepare_query_acceleration(sys), 0);
    alea_slice_view_t view;
    alea_slice_view_axis(&view, 2, 0.0, -2.0, 6.0, -2.0, 6.0);
    alea_ray_slice_validation_options_t options;
    alea_ray_slice_validation_options_init(&options);
    options.checks = ALEA_RAY_SLICE_VALIDATE_FAST_BIDIRECTIONAL;
    alea_ray_slice_validation_result_t* validation =
        alea_ray_slice_validation_result_create();
    ASSERT_NOT_NULL(validation);
    BENCH_START();
    alea_slice_directional_event_cache_t* cache =
        alea_slice_directional_event_cache_create(sys, &view, samples, samples);
    ASSERT_NOT_NULL(cache);
    ASSERT_EQ(alea_validate_ray_slice_compact_with_event_cache(
                  sys, &view, samples, &options, NULL, NULL, cache,
                  validation), 0);
    BENCH_END("64x64 lattice shared validator + event cache", 1);
    printf("[%zu intervals, cache-hit traces=0x%x]  ",
           alea_ray_slice_validation_interval_count(validation),
           alea_ray_slice_validation_reused_trace_mask(validation));
    alea_slice_directional_event_cache_destroy(cache);
    alea_ray_slice_validation_result_destroy(validation);
    mcnp_model_destroy(model);
}

TEST(perf_nested_fill_validator_and_event_cache) {
    const int rows = 128, iterations = 3;
    mcnp_model_t* model = mcnp_load("tests/data/nested_fill.mcnp");
    if (!model) SKIP("nested-fill fixture not found");
    alea_system_t* sys = model->sys;
    ASSERT_EQ(alea_prepare_query_acceleration(sys), 0);
    alea_slice_view_t view;
    alea_slice_view_axis(&view, 2, 0.0, -10.0, 10.0, -10.0, 10.0);
    alea_ray_slice_validation_options_t options;
    alea_ray_slice_validation_options_init(&options);
    options.checks = ALEA_RAY_SLICE_VALIDATE_FAST_BIDIRECTIONAL;
    alea_ray_slice_validation_result_t* validation =
        alea_ray_slice_validation_result_create();
    ASSERT_NOT_NULL(validation);

    {
        BENCH_START();
        for (int i = 0; i < iterations; i++)
            ASSERT_EQ(alea_validate_ray_slice_compact(sys, &view, rows, &options,
                                                       NULL, NULL, validation), 0);
        BENCH_END("128-row nested-fill fast bidirectional validator", iterations);
        printf("[cache-miss fallback]  ");
    }

    {
        BENCH_START();
        for (int i = 0; i < iterations; i++) {
            alea_slice_directional_event_cache_t* cache =
                alea_slice_directional_event_cache_create(sys, &view, rows, rows);
            ASSERT_NOT_NULL(cache);
            alea_slice_directional_event_cache_destroy(cache);
        }
        BENCH_END("128x128 nested-fill directional event cache", iterations);
        printf("[canonical provenance]  ");
    }

    {
        BENCH_START();
        for (int i = 0; i < iterations; i++) {
            alea_slice_directional_event_cache_t* cache =
                alea_slice_directional_event_cache_create(sys, &view, rows, rows);
            ASSERT_NOT_NULL(cache);
            ASSERT_EQ(alea_validate_ray_slice_compact_with_event_cache(
                          sys, &view, rows, &options, NULL, NULL, cache,
                          validation), 0);
            ASSERT_NOT_NULL(
                alea_ray_slice_validation_u_enter_provenance_flags_internal(validation));
            alea_slice_directional_event_cache_destroy(cache);
        }
        BENCH_END("128-row nested-fill shared validator + event cache", iterations);
        printf("[shared consumers, cache-hit traces=0x%x]  ",
               alea_ray_slice_validation_reused_trace_mask(validation));
    }

    alea_ray_slice_validation_result_destroy(validation);
    mcnp_model_destroy(model);
}

/*
 * Native early-stop packet measurement.  Both cases use the same prepared
 * hierarchy, packed rays, acceptance policy, and requested fields.  The
 * scalar loop is the pre-packet baseline; the batch call includes output
 * allocation and direct SoA publication.
 */
TEST(perf_first_visible_native_packet_vs_scalar_20_shells) {
    const size_t n_rays = 10000;
    alea_system_t* sys = build_concentric_shells(20);
    double* origins = calloc(n_rays * 3, sizeof(*origins));
    double* directions = calloc(n_rays * 3, sizeof(*directions));
    alea_raycast_result_t scratch;
    alea_ray_first_visible_batch_result_t batch;
    const alea_ray_batch_query_t query = {
        .kind = ALEA_RAY_QUERY_FIRST_VISIBLE,
        .fields = ALEA_RAY_QUERY_FIELD_DENSITY |
                  ALEA_RAY_QUERY_FIELD_SURFACE_ID |
                  ALEA_RAY_QUERY_FIELD_SURFACE_NORMAL |
                  ALEA_RAY_QUERY_FIELD_RESOLUTION_FLAGS,
        .material_filter = -1
    };
    uint32_t rng = 42;
    size_t scalar_found = 0;

    ASSERT_TRUE(sys && origins && directions);
    for (size_t i = 0; i < n_rays; i++) {
        bench_random_ray(&rng, 15.0,
                         &origins[i * 3], &origins[i * 3 + 1], &origins[i * 3 + 2],
                         &directions[i * 3], &directions[i * 3 + 1],
                         &directions[i * 3 + 2]);
    }
    ASSERT_EQ(alea_raycast_ensure_hier_caches(sys), 0);
    alea_raycast_result_init(&scratch);
    alea_ray_first_visible_batch_result_init(&batch);

    {
        BENCH_START();
        for (size_t i = 0; i < n_rays; i++) {
            alea_ray_t ray;
            alea_ray_first_visible_result_t visible;
            ASSERT_EQ(alea_ray_init(&ray,
                                    origins[i * 3], origins[i * 3 + 1], origins[i * 3 + 2],
                                    directions[i * 3], directions[i * 3 + 1],
                                    directions[i * 3 + 2]), 0);
            ASSERT_EQ(alea_raycast_hier_first_visible_nocache(
                          sys, &ray, 0.0, 0.0, -1, 1, &scratch, &visible), 0);
            scalar_found += visible.found ? 1u : 0u;
        }
        BENCH_END("20 shells scalar first-visible", n_rays);
        printf("[%zu found]  ", scalar_found);
    }

    {
        BENCH_START();
        ASSERT_EQ(alea_raycast_hier_first_visible_batch_nocache(
                      sys, origins, directions, n_rays, &query, &batch), 0);
        BENCH_END("20 shells native first-visible packet", n_rays);
        size_t batch_found = 0;
        for (size_t i = 0; i < n_rays; i++) batch_found += batch.found[i] ? 1u : 0u;
        printf("[%zu found, %zu published bytes]  ", batch_found,
               n_rays * (sizeof(*batch.found) + sizeof(*batch.t) +
                         sizeof(*batch.cell_ids) + sizeof(*batch.material_ids) +
                         sizeof(*batch.densities) + sizeof(*batch.surface_ids) +
                         sizeof(*batch.primitive_ids) + sizeof(*batch.resolution_flags) +
                         3 * sizeof(*batch.normals_xyz)));
        ASSERT_EQ(batch_found, scalar_found);
    }

    alea_ray_first_visible_batch_result_free(&batch);
    alea_raycast_result_free(&scratch);
    free(directions);
    free(origins);
    alea_destroy(sys);
}

TEST(perf_first_visible_native_packet_vs_scalar_lattice) {
    const size_t n_rays = 10000;
    mcnp_model_t* model = mcnp_load("tests/data/mcnp_lattice_eval.mcnp");
    if (!model) SKIP("lattice fixture not found");
    alea_system_t* sys = model->sys;
    double* origins = calloc(n_rays * 3, sizeof(*origins));
    double* directions = calloc(n_rays * 3, sizeof(*directions));
    alea_raycast_result_t scratch;
    alea_ray_first_visible_batch_result_t batch;
    const alea_ray_batch_query_t query = {
        .kind = ALEA_RAY_QUERY_FIRST_VISIBLE,
        .fields = ALEA_RAY_QUERY_FIELD_SURFACE_ID |
                  ALEA_RAY_QUERY_FIELD_RESOLUTION_FLAGS,
        .material_filter = -1
    };
    size_t scalar_found = 0;

    ASSERT_TRUE(origins && directions);
    for (size_t i = 0; i < n_rays; i++) {
        origins[i * 3] = -2.0;
        origins[i * 3 + 1] = -0.95 + 0.1 * (double)(i % 50);
        origins[i * 3 + 2] = 0.0;
        directions[i * 3] = 1.0;
    }
    ASSERT_EQ(alea_prepare_query_acceleration(sys), 0);
    alea_raycast_result_init(&scratch);
    alea_ray_first_visible_batch_result_init(&batch);

    {
        BENCH_START();
        for (size_t i = 0; i < n_rays; i++) {
            alea_ray_t ray;
            alea_ray_first_visible_result_t visible;
            ASSERT_EQ(alea_ray_init(&ray,
                                    origins[i * 3], origins[i * 3 + 1], origins[i * 3 + 2],
                                    directions[i * 3], directions[i * 3 + 1],
                                    directions[i * 3 + 2]), 0);
            ASSERT_EQ(alea_raycast_hier_first_visible_nocache(
                          sys, &ray, 0.0, 0.0, -1, 0, &scratch, &visible), 0);
            scalar_found += visible.found ? 1u : 0u;
        }
        BENCH_END("lattice scalar first-visible", n_rays);
        printf("[%zu found]  ", scalar_found);
    }
    {
        BENCH_START();
        ASSERT_EQ(alea_raycast_hier_first_visible_batch_nocache(
                      sys, origins, directions, n_rays, &query, &batch), 0);
        BENCH_END("lattice native first-visible packet", n_rays);
        size_t batch_found = 0;
        for (size_t i = 0; i < n_rays; i++) batch_found += batch.found[i] ? 1u : 0u;
        printf("[%zu found]  ", batch_found);
        ASSERT_EQ(batch_found, scalar_found);
    }

    alea_ray_first_visible_batch_result_free(&batch);
    alea_raycast_result_free(&scratch);
    free(directions);
    free(origins);
    mcnp_model_destroy(model);
}

TEST(perf_any_hit_native_packet_vs_scalar_20_shells) {
    const size_t n_rays = 10000;
    alea_system_t* sys = build_concentric_shells(20);
    double* origins = calloc(n_rays * 3, sizeof(*origins));
    double* directions = calloc(n_rays * 3, sizeof(*directions));
    alea_raycast_result_t scratch;
    alea_ray_any_hit_batch_result_t batch;
    const alea_ray_batch_query_t query = {
        .kind = ALEA_RAY_QUERY_ANY_HIT,
        .material_filter = -1
    };
    uint32_t rng = 42;
    size_t scalar_hits = 0;

    ASSERT_TRUE(sys && origins && directions);
    for (size_t i = 0; i < n_rays; i++) {
        bench_random_ray(&rng, 15.0,
                         &origins[i * 3], &origins[i * 3 + 1], &origins[i * 3 + 2],
                         &directions[i * 3], &directions[i * 3 + 1],
                         &directions[i * 3 + 2]);
    }
    ASSERT_EQ(alea_raycast_ensure_hier_caches(sys), 0);
    alea_raycast_result_init(&scratch);
    alea_ray_any_hit_batch_result_init(&batch);

    {
        BENCH_START();
        for (size_t i = 0; i < n_rays; i++) {
            alea_ray_t ray;
            int hit = 0;
            ASSERT_EQ(alea_ray_init(&ray,
                                    origins[i * 3], origins[i * 3 + 1], origins[i * 3 + 2],
                                    directions[i * 3], directions[i * 3 + 1],
                                    directions[i * 3 + 2]), 0);
            ASSERT_EQ(alea_raycast_hier_any_hit_nocache(
                          sys, &ray, 0.0, 0.0, -1, &scratch, &hit), 0);
            scalar_hits += hit ? 1u : 0u;
        }
        BENCH_END("20 shells scalar any-hit", n_rays);
        printf("[%zu hits]  ", scalar_hits);
    }
    {
        BENCH_START();
        ASSERT_EQ(alea_raycast_hier_any_hit_batch_nocache(
                      sys, origins, directions, n_rays, &query, &batch), 0);
        BENCH_END("20 shells native any-hit packet", n_rays);
        size_t batch_hits = 0;
        for (size_t i = 0; i < n_rays; i++) batch_hits += batch.hits[i] ? 1u : 0u;
        printf("[%zu hits, %zu published bytes]  ", batch_hits,
               n_rays * sizeof(*batch.hits));
        ASSERT_EQ(batch_hits, scalar_hits);
    }

    alea_ray_any_hit_batch_result_free(&batch);
    alea_raycast_result_free(&scratch);
    free(directions);
    free(origins);
    alea_destroy(sys);
}

/*
 * First compact-ABI measurement gate.  Both cases use the same prepared
 * hierarchy and the identical packed rays.  The single-ray loop measures the
 * existing rich-result call boundary; the compact call includes its native
 * CSR assembly but no consumer-language conversion.
 */
TEST(perf_compact_hier_batch_20_shells) {
    const size_t n_rays = 10000;
    alea_system_t* sys = build_concentric_shells(20);
    double* origins = calloc(n_rays * 3, sizeof(*origins));
    double* directions = calloc(n_rays * 3, sizeof(*directions));
    alea_raycast_result_t* single = alea_raycast_result_create();
    alea_raycast_batch_result_t* batch = alea_raycast_batch_result_create();
    alea_raycast_batch_options_t options = {
        .struct_size = sizeof(options),
        .fields = ALEA_RAY_BATCH_MATERIAL | ALEA_RAY_BATCH_DENSITY |
                  ALEA_RAY_BATCH_SURFACES | ALEA_RAY_BATCH_RESOLUTION_FLAGS
    };
    uint32_t rng = 42;
    size_t single_segments = 0;

    ASSERT_TRUE(origins && directions && single && batch);
    for (size_t i = 0; i < n_rays; i++) {
        bench_random_ray(&rng, 15.0,
                         &origins[i * 3], &origins[i * 3 + 1], &origins[i * 3 + 2],
                         &directions[i * 3], &directions[i * 3 + 1],
                         &directions[i * 3 + 2]);
    }
    alea_raycast_ensure_caches(sys);

    {
        BENCH_START();
        for (size_t i = 0; i < n_rays; i++) {
            ASSERT_EQ(alea_raycast_hier_fast_segments(
                          sys,
                          origins[i * 3], origins[i * 3 + 1], origins[i * 3 + 2],
                          directions[i * 3], directions[i * 3 + 1], directions[i * 3 + 2],
                          100.0, single),
                      0);
            single_segments += alea_raycast_segment_count(single);
        }
        BENCH_END("20 shells hierarchical single-ray", n_rays);
        printf("[single-ray]  ");
    }
    printf("[%zu segments]  ", single_segments);

    {
        BENCH_START();
        ASSERT_EQ(alea_raycast_hier_batch(sys, origins, directions, n_rays, 100.0,
                                           &options, batch), 0);
        BENCH_END("20 shells compact batch (trace + CSR)", n_rays);
        printf("[compact]  ");
    }
    printf("[%zu segments]  ", alea_raycast_batch_segment_count(batch));
    ASSERT_EQ(alea_raycast_batch_segment_count(batch), single_segments);
    printf("[%zu published bytes, %zu known transient input bytes]  ",
           batch_published_bytes(batch, options.fields),
           n_rays * 6 * sizeof(double));
    alea_raycast_batch_work_stats_t work_stats;
    ASSERT_EQ(alea_raycast_batch_result_get_work_stats_internal(
                  batch, &work_stats), 0);
    printf("[max owner neighbor=%llu/%llu path=%llu/%llu root=%llu/%llu "
           "full=%llu/%llu boundary=%llu snapshots=%llu/%llu "
           "growths=%llu/%lluB]  ",
           (unsigned long long)work_stats.max_owner_neighbor_hits,
           (unsigned long long)work_stats.max_owner_neighbor_attempts,
           (unsigned long long)work_stats.max_owner_path_hits,
           (unsigned long long)work_stats.max_owner_path_attempts,
           (unsigned long long)work_stats.max_owner_root_hits,
           (unsigned long long)work_stats.max_owner_root_queries,
           (unsigned long long)work_stats.max_owner_full_hits,
           (unsigned long long)work_stats.max_owner_full_queries,
           (unsigned long long)work_stats.max_boundary_event_enrichments,
           (unsigned long long)work_stats.max_path_snapshot_copies,
           (unsigned long long)work_stats.max_path_snapshot_entries,
           (unsigned long long)work_stats.max_result_buffer_growths,
           (unsigned long long)work_stats.max_result_buffer_growth_bytes);

    alea_raycast_batch_result_destroy(batch);
    alea_raycast_result_destroy(single);
    free(directions);
    free(origins);
    alea_destroy(sys);
}

typedef struct {
    size_t* segment_counts;
} xray_batch_count_context_t;

static int count_xray_batch_segment(void* context, size_t ray_index,
                                    const alea_ray_segment_t* segment) {
    (void)segment;
    xray_batch_count_context_t* counts = context;
    counts->segment_counts[ray_index]++;
    return 0;
}

/* Mirrors the non-lattice, single-sample X-ray renderer's ray layout.  The
 * scalar side retains one scratch result, while the compact side streams into
 * fixed per-ray slots; timings are informational and segment totals must
 * remain identical. */
TEST(perf_xray_camera_tiles_compact_vs_reusable_scalar) {
    const int width = 96, height = 96, tile = 32;
    alea_system_t* sys = build_concentric_shells(20);
    ASSERT_NOT_NULL(sys);
    ASSERT_EQ(alea_raycast_ensure_hier_caches(sys), 0);
    render_config_t cfg;
    render_config_init(&cfg);
    cfg.width = width; cfg.height = height;
    cfg.eye[0] = 0; cfg.eye[1] = -20; cfg.eye[2] = 0;
    cfg.target[0] = 0; cfg.target[1] = 0; cfg.target[2] = 0;
    cfg.eye_set = 1; cfg.target_set = 1;
    render_camera_t cam;
    ASSERT_EQ(render_camera_setup(&cam, &cfg, sys), 0);
    const size_t tile_capacity = (size_t)tile * tile;
    double* origins = malloc(tile_capacity * 3 * sizeof(*origins));
    double* directions = malloc(tile_capacity * 3 * sizeof(*directions));
    size_t* segment_counts = calloc(tile_capacity, sizeof(*segment_counts));
    alea_raycast_result_t scratch;
    alea_raycast_result_init(&scratch);
    alea_raycast_result_reserve(&scratch, 0, 64);
    ASSERT_NOT_NULL(origins);
    ASSERT_NOT_NULL(directions);
    ASSERT_NOT_NULL(segment_counts);
    size_t scalar_segments = 0, compact_segments = 0;

    {
        BENCH_START();
        for (int y = 0; y < height; y++) for (int x = 0; x < width; x++) {
            double ox, oy, oz, dx, dy, dz;
            render_camera_ray(&cam, width, height, x + 0.5, y + 0.5,
                              &ox, &oy, &oz, &dx, &dy, &dz);
            alea_ray_t ray;
            ASSERT_EQ(alea_ray_init(&ray, ox, oy, oz, dx, dy, dz), 0);
            alea_raycast_result_clear(&scratch);
            ASSERT_EQ(alea_raycast_hier_segments_nocache(sys, &ray, 0, &scratch), 0);
            scalar_segments += scratch.segments.count;
        }
        BENCH_END("X-ray camera reusable scalar segments", width * height);
    }
    printf("[%zu segments]  ", scalar_segments);

    {
        BENCH_START();
        for (int y0 = 0; y0 < height; y0 += tile) {
            for (int x0 = 0; x0 < width; x0 += tile) {
            const int x1 = x0 + tile < width ? x0 + tile : width;
            const int y1 = y0 + tile < height ? y0 + tile : height;
            size_t ray_count = 0;
            for (int y = y0; y < y1; y++) for (int x = x0; x < x1; x++, ray_count++)
                render_camera_ray(&cam, width, height, x + 0.5, y + 0.5,
                                  &origins[ray_count * 3], &origins[ray_count * 3 + 1],
                                  &origins[ray_count * 3 + 2], &directions[ray_count * 3],
                                  &directions[ray_count * 3 + 1], &directions[ray_count * 3 + 2]);
            memset(segment_counts, 0, ray_count * sizeof(*segment_counts));
            xray_batch_count_context_t count_context = {
                .segment_counts = segment_counts
            };
            ASSERT_EQ(alea_raycast_hier_visit_segments_batch_nocache(
                          sys, origins, directions, ray_count, 0,
                          count_xray_batch_segment, &count_context), 0);
            for (size_t i = 0; i < ray_count; i++)
                compact_segments += segment_counts[i];
            }
        }
        BENCH_END("X-ray camera compact 32x32 fixed-output visitor", width * height);
    }
    printf("[%zu segments]  ", compact_segments);
    ASSERT_EQ(compact_segments, scalar_segments);
    printf("[%zu retained scalar bytes, %zu known tile-input bytes]  ",
           raycast_result_retained_bytes(&scratch), tile_capacity * 6 * sizeof(double));

    alea_raycast_result_free(&scratch);
    free(origins); free(directions); free(segment_counts);
    render_config_free(&cfg);
    alea_destroy(sys);
}

/* End-to-end scheduling measurement for the renderer's fixed-output X-ray
 * path.  This intentionally includes camera setup, tile ownership, and
 * framebuffer writes rather than extrapolating from traversal alone. */
TEST(perf_xray_render_frame_fixed_tile_scheduling) {
    const int width = 96, height = 96;
    alea_system_t* sys = build_concentric_shells(20);
    ASSERT_NOT_NULL(sys);
    render_config_t cfg;
    render_config_init(&cfg);
    cfg.width = width;
    cfg.height = height;
    cfg.tile_size = 32;
    cfg.aa_samples = 1;
    cfg.render_mode = RENDER_MODE_XRAY;
    cfg.eye[0] = 0; cfg.eye[1] = -20; cfg.eye[2] = 0;
    cfg.target[0] = 0; cfg.target[1] = 0; cfg.target[2] = 0;
    cfg.eye_set = 1; cfg.target_set = 1;
    render_camera_t cam;
    ASSERT_EQ(render_camera_setup(&cam, &cfg, sys), 0);
    render_framebuffer_t* framebuffer =
        render_framebuffer_create(width, height, 0);
    ASSERT_NOT_NULL(framebuffer);
    {
        BENCH_START();
        ASSERT_EQ(render_scene(sys, &cfg, &cam, framebuffer), 0);
        BENCH_END("X-ray render 96x96 fixed tile scheduling", width * height);
    }
    render_framebuffer_free(framebuffer);
    render_config_free(&cfg);
    alea_destroy(sys);
}

/* Two deterministic lattice-entry cases retained from the hierarchy
 * correctness suite.  They make the entry-search attribution observable in
 * the normal performance report without requiring the optional E-lite deck:
 * one transformed occurrence is pruned from trusted ancestor support before
 * DDA, while the other enters a thin active ancestor shell without stepping
 * through its large conservative AABB. */
TEST(perf_lattice_entry_attribution) {
    const char* transformed_prune_input =
        "Inactive transformed repeating lattice placement\n"
        "1 0 -1 FILL=10 (5 0 0)\n"
        "100 0 -8 9 -10 11 -12 13 LAT=1 U=10 FILL=1\n"
        "10 1 -1.0 -14 U=1\n"
        "11 2 -1.0 14 U=1\n"
        "\n"
        "1 SO 1\n"
        "2 PX -1\n"
        "3 PY -1\n"
        "4 PY 1\n"
        "5 PZ -1\n"
        "6 PZ 1\n"
        "8 PX -1\n"
        "9 PX 1\n"
        "10 PY -1\n"
        "11 PY 1\n"
        "12 PZ -1\n"
        "13 PZ 1\n"
        "14 CZ 0.3\n"
        "\n"
        "M1 1001.80c 1.0\n"
        "M2 1001.80c 1.0\n";
    const char* exact_support_input =
        "Repeating lattice exact ancestor support entry\n"
        "1 0 -1 2 FILL=10\n"
        "100 0 -3 LAT=1 U=10 FILL=1\n"
        "10 1 -1.0 -9 U=1\n"
        "\n"
        "1 SO 100\n"
        "2 SO 99\n"
        "3 RPP -0.5 0.5 -0.5 0.5 -0.5 0.5\n"
        "9 SO 0.25\n"
        "\n"
        "M1 1001.80c 1.0\n";
    const int iterations = 1000;
    alea_raycast_result_t trace;

    mcnp_model_t* transformed =
        mcnp_load_string(transformed_prune_input, strlen(transformed_prune_input));
    ASSERT_NOT_NULL(transformed);
    ASSERT_EQ(alea_prepare_query_acceleration(transformed->sys), 0);
    alea_raycast_result_init(&trace);
    {
        BENCH_START();
        for (int i = 0; i < iterations; i++) {
            alea_raycast_result_clear(&trace);
            ASSERT_EQ(alea_raycast_hier_fast_segments(
                          transformed->sys, 3.0, 0.0, 0.0,
                          1.0, 0.0, 0.0, 10000.0, &trace), 0);
        }
        BENCH_END("lattice entry transformed-support prune", iterations);
    }
    printf("[calls=%zu TLAS nodes=%zu leaves=%zu candidates=%zu DDA=%zu future=%zu]  ",
           trace.lattice_entry_calls, trace.lattice_entry_tlas_nodes_tested,
           trace.lattice_entry_tlas_leaves_visited, trace.lattice_entry_candidates,
           trace.lattice_entry_dda_steps, trace.lattice_entry_future_entry_results);
    ASSERT(trace.lattice_entry_calls > 0);
    ASSERT(trace.lattice_entry_tlas_nodes_tested > 0);
    ASSERT_EQ(trace.lattice_entry_tlas_leaves_visited, 0);
    ASSERT_EQ(trace.lattice_entry_candidates, 0);
    ASSERT_EQ(trace.lattice_entry_dda_steps, 0);
    alea_raycast_result_free(&trace);
    mcnp_model_destroy(transformed);

    mcnp_model_t* exact_support =
        mcnp_load_string(exact_support_input, strlen(exact_support_input));
    ASSERT_NOT_NULL(exact_support);
    ASSERT_EQ(alea_prepare_query_acceleration(exact_support->sys), 0);
    alea_raycast_result_init(&trace);
    {
        BENCH_START();
        for (int i = 0; i < iterations; i++) {
            alea_raycast_result_clear(&trace);
            ASSERT_EQ(alea_raycast_hier_fast_segments(
                          exact_support->sys, 0.0, 0.0, 0.0,
                          1.0, 0.0, 0.0, 110.0, &trace), 0);
        }
        BENCH_END("lattice entry exact ancestor support", iterations);
    }
    printf("[calls=%zu TLAS nodes=%zu leaves=%zu candidates=%zu DDA=%zu ancestor tests=%zu events=%zu "
           "future=%zu canonical rejects=%zu]  ",
           trace.lattice_entry_calls, trace.lattice_entry_tlas_nodes_tested,
           trace.lattice_entry_tlas_leaves_visited, trace.lattice_entry_candidates,
           trace.lattice_entry_dda_steps,
           trace.lattice_entry_ancestor_surface_tests,
           trace.lattice_entry_ancestor_events,
           trace.lattice_entry_future_entry_results,
           trace.lattice_entry_canonical_rejections);
    ASSERT(trace.lattice_entry_calls > 0);
    ASSERT(trace.lattice_entry_tlas_nodes_tested > 0);
    ASSERT(trace.lattice_entry_tlas_leaves_visited > 0);
    ASSERT(trace.lattice_entry_candidates > 0);
    ASSERT(trace.lattice_entry_ancestor_surface_tests > 0);
    ASSERT(trace.lattice_entry_ancestor_events > 0);
    ASSERT(trace.lattice_entry_dda_steps < 10);

    {
        enum { batch_rays = 8 };
        double origins[batch_rays * 3] = {0};
        double directions[batch_rays * 3] = {0};
        for (size_t i = 0; i < batch_rays; i++) directions[i * 3] = 1.0;
        alea_raycast_batch_options_t options = {
            .struct_size = sizeof(options),
            .fields = ALEA_RAY_BATCH_MATERIAL
        };
        alea_raycast_batch_result_t* batch = alea_raycast_batch_result_create();
        alea_raycast_batch_work_stats_t batch_stats;
        ASSERT_NOT_NULL(batch);
        ASSERT_EQ(alea_raycast_hier_batch(exact_support->sys, origins, directions,
                                          batch_rays, 110.0, &options, batch), 0);
        ASSERT_EQ(alea_raycast_batch_result_get_work_stats_internal(
                      batch, &batch_stats), 0);
        printf("[batch max calls=%llu TLAS nodes=%llu leaves=%llu candidates=%llu "
               "DDA=%llu ancestor tests=%llu events=%llu]  ",
               (unsigned long long)batch_stats.max_lattice_entry_calls,
               (unsigned long long)batch_stats.max_lattice_entry_tlas_nodes_tested,
               (unsigned long long)batch_stats.max_lattice_entry_tlas_leaves_visited,
               (unsigned long long)batch_stats.max_lattice_entry_candidates,
               (unsigned long long)batch_stats.max_lattice_entry_dda_steps,
               (unsigned long long)batch_stats.max_lattice_entry_ancestor_surface_tests,
               (unsigned long long)batch_stats.max_lattice_entry_ancestor_events);
        ASSERT_EQ(batch_stats.max_lattice_entry_calls, trace.lattice_entry_calls);
        ASSERT_EQ(batch_stats.max_lattice_entry_tlas_nodes_tested,
                  trace.lattice_entry_tlas_nodes_tested);
        ASSERT_EQ(batch_stats.max_lattice_entry_tlas_leaves_visited,
                  trace.lattice_entry_tlas_leaves_visited);
        ASSERT_EQ(batch_stats.max_lattice_entry_candidates,
                  trace.lattice_entry_candidates);
        ASSERT_EQ(batch_stats.max_lattice_entry_dda_steps,
                  trace.lattice_entry_dda_steps);
        ASSERT_EQ(batch_stats.max_lattice_entry_ancestor_surface_tests,
                  trace.lattice_entry_ancestor_surface_tests);
        ASSERT_EQ(batch_stats.max_lattice_entry_ancestor_events,
                  trace.lattice_entry_ancestor_events);
        alea_raycast_batch_result_destroy(batch);
    }
    alea_raycast_result_free(&trace);
    mcnp_model_destroy(exact_support);
}

/* The two paths below execute the same selected-owner walker.  Keeping the
 * rich result reserved makes their steady-state difference attributable to
 * result publication and segment-vector writes, rather than allocator growth. */
TEST(perf_selected_segment_publication_vs_streaming) {
    const int iterations = 10000;
    alea_system_t* sys = build_concentric_shells(20);
    ASSERT_NOT_NULL(sys);
    ASSERT_EQ(alea_raycast_ensure_hier_caches(sys), 0);
    alea_ray_t ray;
    ASSERT_EQ(alea_ray_init(&ray, -20.0, 0.125, 0.0, 1.0, 0.0, 0.0), 0);
    alea_raycast_result_t rich, streamed;
    alea_raycast_result_init(&rich);
    alea_raycast_result_init(&streamed);
    alea_raycast_result_reserve(&rich, 0, 64);
    size_t rich_segments = 0, streamed_segments = 0;

    {
        BENCH_START();
        for (int i = 0; i < iterations; i++) {
            alea_raycast_result_clear(&rich);
            ASSERT_EQ(alea_raycast_hier_segments_nocache(sys, &ray, 100.0,
                                                          &rich), 0);
            rich_segments += rich.segments.count;
        }
        BENCH_END("selected intervals rich segment publication", iterations);
    }
    printf("[%zu segments, %zu retained bytes]  ", rich_segments,
           raycast_result_retained_bytes(&rich));

    {
        BENCH_START();
        for (int i = 0; i < iterations; i++) {
            ASSERT_EQ(alea_raycast_hier_visit_segments_nocache(
                          sys, &ray, 100.0, &streamed,
                          count_selected_segment, &streamed_segments), 0);
        }
        BENCH_END("selected intervals streamed fixed output", iterations);
    }
    printf("[%zu segments, %zu retained bytes]  ", streamed_segments,
           raycast_result_retained_bytes(&streamed));
    ASSERT_EQ(streamed_segments, rich_segments);
    ASSERT_EQ(streamed.segments.count, 0);
    printf("[owner neighbor=%llu/%llu path=%llu/%llu root=%llu/%llu "
           "full=%llu/%llu boundary=%llu snapshots=%llu/%llu "
           "growths=%llu/%lluB]  ",
           (unsigned long long)streamed.owner_neighbor_hits,
           (unsigned long long)streamed.owner_neighbor_attempts,
           (unsigned long long)streamed.owner_path_hits,
           (unsigned long long)streamed.owner_path_attempts,
           (unsigned long long)streamed.owner_root_hits,
           (unsigned long long)streamed.owner_root_queries,
           (unsigned long long)streamed.owner_full_hits,
           (unsigned long long)streamed.owner_full_queries,
           (unsigned long long)streamed.boundary_event_enrichments,
           (unsigned long long)streamed.path_snapshot_copies,
           (unsigned long long)streamed.path_snapshot_entries,
           (unsigned long long)streamed.result_buffer_growths,
           (unsigned long long)streamed.result_buffer_growth_bytes);
    ASSERT_EQ(streamed.owner_neighbor_attempts, rich.owner_neighbor_attempts);
    ASSERT_EQ(streamed.owner_neighbor_hits, rich.owner_neighbor_hits);
    ASSERT_EQ(streamed.owner_path_attempts, rich.owner_path_attempts);
    ASSERT_EQ(streamed.owner_path_hits, rich.owner_path_hits);
    ASSERT_EQ(streamed.owner_root_queries, rich.owner_root_queries);
    ASSERT_EQ(streamed.owner_root_hits, rich.owner_root_hits);
    ASSERT_EQ(streamed.owner_full_queries, rich.owner_full_queries);
    ASSERT_EQ(streamed.owner_full_hits, rich.owner_full_hits);
    ASSERT_EQ(streamed.boundary_event_enrichments,
              rich.boundary_event_enrichments);
    ASSERT_EQ(streamed.path_snapshot_copies, 0);
    ASSERT_EQ(rich.path_snapshot_copies, 0);
    ASSERT_EQ(streamed.result_buffer_growths, 0);
    ASSERT_EQ(rich.result_buffer_growths, 0);
    ASSERT_EQ(streamed.selected_intervals_yielded,
              rich.selected_intervals_yielded);

    alea_raycast_result_free(&streamed);
    alea_raycast_result_free(&rich);
    alea_destroy(sys);
}

TEST(perf_query_lowering) {
    const int iterations = 1000000;
    const alea_ray_query_t query = {
        .kind = ALEA_RAY_QUERY_FIRST_VISIBLE,
        .backend = ALEA_RAY_QUERY_BACKEND_AUTO,
        .fields = ALEA_RAY_QUERY_FIELD_CELL_ID |
                  ALEA_RAY_QUERY_FIELD_MATERIAL_ID,
        .t_min = 0.0,
        .t_max = 100.0,
        .material_filter = -1
    };
    alea_ray_plan_t plan;
    uint64_t checksum = 0;
    int status = 0;
    {
        BENCH_START();
        for (int i = 0; i < iterations; i++) {
            status |= alea_ray_query_lower(&query, &plan);
            checksum += (uint64_t)plan.engine + plan.requirements.need_selected_owner;
        }
        BENCH_END("semantic query lowering", iterations);
    }
    ASSERT_EQ(status, 0);
    ASSERT_EQ(checksum, (uint64_t)iterations);
}

TEST(perf_query_cache_prepare_audit) {
    alea_system_t* sys = build_concentric_shells(20);
    ASSERT_NOT_NULL(sys);
    ASSERT(!alea_system_query_cache_ready(sys, ALEA_CACHE_RAYCAST));
    {
        BENCH_START();
        ASSERT_EQ(alea_raycast_ensure_caches(sys), 0);
        BENCH_END("query-cache cold preparation (20 shells)", 1);
    }
    ASSERT(alea_system_query_cache_ready(sys, ALEA_CACHE_RAYCAST));
    {
        BENCH_START();
        ASSERT_EQ(alea_raycast_ensure_caches(sys), 0);
        BENCH_END("query-cache warm preparation (cache hit)", 1);
    }
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
        total_hits += result.hits.count;
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
