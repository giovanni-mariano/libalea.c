// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
// SPDX-License-Identifier: MPL-2.0

#include "alea_test.h"
#include "alea.h"
#include "alea_raycast.h"
#include "alea_render.h"
#include "raycast/regional/interval_index.h"
#include "raycast/ray_traversal.h"
#include "alea_mcnp.h"
#include "alea_snapshot.h"
#include "util/compat.h"
#include "util/alea_parallel.h"
#include <float.h>
#include "core/alea_eval.h"
#include "core/alea_regional.h"
#include "core/alea_universe.h"
#include "util/alea_alloc.h"
#include <math.h>

static alea_node_id_t sphere(alea_system_t* sys, double cx) {
    int surface = alea_sphere_surface(sys, 0, cx, 0, 0, 1);
    return alea_halfspace(sys, surface, -1);
}

TEST(regional_retains_editable_tree_and_overlapping_owners) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);
    int mat = alea_add_material(sys, 1);
    alea_node_id_t left = sphere(sys, -3), right = sphere(sys, 3);
    alea_node_id_t root = alea_union(sys, left, right);
    for (int id = 1; id <= 2; id++)
        ASSERT(alea_add_cell(sys, id, root, mat, -1, 0) >= 0);
    FILE* original = tmpfile();
    FILE* prepared = tmpfile();
    ASSERT_NOT_NULL(original);
    ASSERT_NOT_NULL(prepared);
    ASSERT_EQ(mcnp_export_system_stream(sys, original), 0);

    alea_cell_hit_t before[4], after[4];
    int count = alea_find_all_cells(sys, -3, 0, 0, before, 4);
    ASSERT_EQ(count, 2);
    ASSERT_EQ(alea_prepare_regional_queries(sys), 0);
    ASSERT_EQ(alea_find_all_cells(sys, -3, 0, 0, after, 4), count);
    for (int i = 0; i < count; i++) ASSERT_EQ(before[i].cell_id, after[i].cell_id);
    ASSERT_EQ(alea_prepare_regional_queries(sys), 0);
    ASSERT_EQ(sys->cells.data[0].root_node_id, root);
    ASSERT_EQ(alea_node_left(sys, root), left);
    ASSERT_EQ(alea_node_right(sys, root), right);
    ASSERT_EQ(mcnp_export_system_stream(sys, prepared), 0);
    rewind(original);
    rewind(prepared);
    int ch;
    do {
        ch = fgetc(original);
        ASSERT_EQ(ch, fgetc(prepared));
    } while (ch != EOF);
    fclose(original);
    fclose(prepared);

    for (int i = -4; i <= 4; i++) {
        double xs[] = {nextafter(i, -INFINITY), i, nextafter(i, INFINITY)};
        for (int j = 0; j < 3; j++) {
            bool inside;
            ASSERT(alea_regional_contains(sys, root, xs[j], 0, 0, &inside));
            ASSERT_EQ(inside, alea_evaluate_point(sys, root, xs[j], 0, 0) <= 0);
        }
    }
    char path[4096];
    FILE* snapshot = alea_tmpfile(path);
    ASSERT_NOT_NULL(snapshot);
    fclose(snapshot);
    ASSERT_EQ(alea_snapshot_write(sys, path), 0);
    alea_system_t* copy = alea_snapshot_read(path, 0);
    remove(path);
    ASSERT_NOT_NULL(copy);
    ASSERT_NULL(copy->regional_cache);
    ASSERT_EQ(copy->cells.data[0].root_node_id, root);
    ASSERT(alea_contains_point(copy, root, -3, 0, 0));
    ASSERT_EQ(alea_prepare_regional_queries(copy), 0);
    ASSERT_NOT_NULL(copy->regional_cache);
    alea_destroy(copy);

    bool ignored;
    ASSERT(!alea_contains_point(sys, root, 5, 0, 0));
    ASSERT(!alea_contains_point(sys, root, NAN, 0, 0));
    ASSERT_EQ(alea_cell_set_region(sys, 0, right), 0);
    ASSERT_NULL(sys->regional_cache);
    ASSERT_EQ(alea_find_all_cells(sys, -3, 0, 0, after, 4), 1);
    ASSERT_EQ(after[0].cell_id, 2);
    ASSERT_EQ(alea_prepare_regional_queries(sys), 0);
    ASSERT(alea_regional_contains(sys, root, -3, 0, 0, &ignored));
    alea_destroy(sys);
}

TEST(regional_plane_folding_preserves_signed_boundaries_and_fallback) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);
    alea_node_id_t ball = sphere(sys, 0);
    int cut = alea_plane_surface(sys, 0, 1, 0, 0, 0);
    alea_node_id_t neg = alea_halfspace(sys, cut, -1);
    alea_node_id_t boundary = alea_difference(sys, ball, neg);
    /* Both diagonal planes are strictly negative throughout the ball box.
     * Keep an opaque difference in a residual with another mixed subtree. */
    int a = alea_plane_surface(sys, 0, 1, 1, 0, -5);
    int b = alea_plane_surface(sys, 0, -1, 1, 0, -5);
    alea_node_id_t clipped = alea_intersection(sys, alea_halfspace(sys, a, -1), neg);
    alea_node_id_t both = alea_union(sys, clipped, alea_halfspace(sys, b, 1));
    alea_node_id_t root = alea_intersection(sys, both, boundary);
    ASSERT(alea_add_cell(sys, 1, root, alea_add_material(sys, 1), -1, 0) >= 0);
    ASSERT_EQ(alea_prepare_regional_queries(sys), 0);

    double xs[] = {-1, -0.5, -DBL_MIN, -0.0, 0.0, DBL_MIN, 0.5, 1};
    for (size_t i = 0; i < sizeof(xs)/sizeof(*xs); i++) {
        bool inside;
        ASSERT(alea_regional_contains(sys, root, xs[i], 0, 0, &inside));
        ASSERT_EQ(inside, alea_evaluate_point(sys, root, xs[i], 0, 0) <= 0);
    }
    ASSERT(alea_contains_point(sys, root, 0, 0, 0));
    alea_perf_reset();
    (void)alea_evaluate_point(sys, root, 0, 0, 0);
    size_t full = alea_perf_get().primitive_evaluations;
    alea_perf_reset();
    ASSERT(alea_contains_point(sys, root, 0, 0, 0));
    ASSERT(alea_perf_get().primitive_evaluations < full);
    alea_destroy(sys);
}

TEST(regional_uses_fill_local_coordinates) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);
    int mat = alea_add_material(sys, 1);
    int envelope = alea_box_surface(sys, 0, 90, 110, -2, 2, -2, 2);
    int parent = alea_add_cell(sys, 1, alea_halfspace(sys, envelope, -1),
                               ALEA_MATERIAL_VOID, 0, 0);
    alea_node_id_t a = sphere(sys, -3), b = sphere(sys, 3);
    alea_node_id_t root = alea_union(sys, a, b);
    ASSERT(alea_add_cell(sys, 2, root, mat, -1, 7) >= 0);
    double shift[] = {100, 0, 0};
    ASSERT_EQ(alea_add_transform(sys, 1, shift, 3, 0), 0);
    ASSERT_EQ(alea_set_fill(sys, parent, 7, 1), 0);
    ASSERT_EQ(alea_prepare_query_acceleration(sys), 0);
    for (int mode = 0; mode < 2; mode++) {
        if (mode) ASSERT_EQ(alea_prepare_regional_queries(sys), 0);
        for (int x = 97; x <= 103; x += 6) {
            alea_cell_hit_t hits[4];
            ASSERT_EQ(alea_find_all_cells(sys, x, 0, 0, hits, 4), 2);
            ASSERT_EQ(hits[1].cell_id, 2);
            ASSERT_EQ(hits[1].local_x, x - 100);
        }
    }
    alea_destroy(sys);
}

TEST(regional_expensive_tree_falls_back) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);
    int plane = alea_plane_surface(sys, 0, 1, 1, 0, 0);
    alea_node_id_t tree = alea_halfspace(sys, plane, -1);
    /* Reused subtrees make this small DAG exceed the per-region visit budget. */
    for (int i = 0; i < 10; i++) tree = alea_union(sys, tree, tree);
    alea_node_id_t root = alea_intersection(sys, tree, sphere(sys, 0));
    ASSERT(alea_add_cell(sys, 1, root, alea_add_material(sys, 1), -1, 0) >= 0);
    ASSERT_EQ(alea_prepare_regional_queries(sys), 0);
    ASSERT_NULL(sys->regional_cache);
    ASSERT(alea_contains_point(sys, root, -0.5, 0, 0));
    ASSERT(!alea_contains_point(sys, root, 0.5, 0, 0));
    alea_destroy(sys);
}

TEST(regional_prepare_failure_can_retry) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);
    alea_node_id_t a = sphere(sys, -3), b = sphere(sys, 3);
    alea_node_id_t root = alea_union(sys, a, b);
    ASSERT(alea_add_cell(sys, 1, root, alea_add_material(sys, 1), -1, 0) >= 0);
    alea_interrupt();
    int rc = alea_prepare_regional_queries(sys);
    alea_clear_interrupt();
    ASSERT_EQ(rc, -1);
    ASSERT_NULL(sys->regional_cache);
    for (int n = 0; n < 3; n++) {
        alea_alloc_fail_after(n);
        rc = alea_prepare_regional_queries(sys);
        alea_alloc_fail_after(-1);
        ASSERT_EQ(rc, -1);
        ASSERT_NULL(sys->regional_cache);
        ASSERT(!alea_system_query_cache_ready(sys, ALEA_CACHE_REGIONAL));
        ASSERT(alea_contains_point(sys, root, 3, 0, 0));
    }
    ASSERT_EQ(alea_prepare_regional_queries(sys), 0);
    ASSERT_NOT_NULL(sys->regional_cache);
    alea_destroy(sys);
}
static bool box_excludes(const alea_system_t* sys, alea_node_id_t root,
                         const alea_bbox_t* box) {
    size_t visits = 512;
    return alea_regional_proves_outside(sys, root, box, &visits);
}

TEST(quadric_grid_proof_handles_interior_and_face_minima) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);
    /* (x-y)^2 + 2z^2 - 8 has both a mixed term and a null direction. */
    int q = alea_quadric_surface(sys, 1, 1, 1, 2, -2, 0, 0, 0, 0, 0, -8);
    ASSERT(q >= 0);
    alea_node_id_t inner = alea_halfspace(sys, q, -1);
    alea_node_id_t outer = alea_halfspace(sys, q, 1);
    alea_bbox_t far = {-1, 1, -1, 1, 5, 6};
    alea_bbox_t near = {-0.25, 0.25, -0.25, 0.25, -0.25, 0.25};
    ASSERT(box_excludes(sys, inner, &far));
    ASSERT(box_excludes(sys, outer, &near));
    ASSERT(!box_excludes(sys, outer, &far));
    ASSERT(!box_excludes(sys, inner, &near));

    /* All eight corners are outside but the middle is inside. The second
     * box has its minimum on the z=2 face, exactly on the surface. */
    alea_bbox_t enclosing = {-3, 3, -3, 3, -3, 3};
    alea_bbox_t face = {-2, 2, -1, 1, 2, 3};
    ASSERT(!box_excludes(sys, inner, &enclosing));
    ASSERT(!box_excludes(sys, inner, &face));
    ASSERT(alea_contains_point(sys, inner, 0, 0, 0));
    ASSERT(alea_contains_point(sys, inner, 0, 0, 2));

    /* A box crossing the root bbox with an extreme span cannot certify
     * finite intermediate products even though z alone suggests exclusion. */
    alea_bbox_t huge = {-DBL_MAX, DBL_MAX, 0, 0, 5, 6};
    ASSERT(!box_excludes(sys, inner, &huge));

    alea_bbox_t surface_point = {0, 0, 0, 0, 2, 2};
    ASSERT(!box_excludes(sys, inner, &surface_point));
    ASSERT(!box_excludes(sys, outer, &surface_point));
    surface_point.min_z = surface_point.max_z = nextafter(2.0, INFINITY);
    ASSERT(!box_excludes(sys, inner, &surface_point));
    alea_destroy(sys);
}

TEST(quadric_grid_proofs_agree_with_scalar_samples_in_both_senses) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);
    /* Indefinite translated quadric with all three mixed terms active. */
    int q = alea_quadric_surface(sys, 1, -2, 0.5, 3, 1, -0.75, 1.25, 4, -2, 3, -7);
    ASSERT(q >= 0);
    for (int sense = -1; sense <= 1; sense += 2) {
        alea_node_id_t root = alea_halfspace(sys, q, sense);
        unsigned exclusions = 0;
        for (int x = -4; x < 4; x++) {
            for (int y = -4; y < 4; y++) {
                for (int z = -4; z < 4; z++) {
                    alea_bbox_t box = {x, x+1, y, y+1, z, z+1};
                    if (!box_excludes(sys, root, &box)) continue;
                    exclusions++;
                    for (int i = 0; i <= 4; i++)
                        for (int j = 0; j <= 4; j++)
                            for (int k = 0; k <= 4; k++)
                                ASSERT(alea_evaluate_point(sys, root, x+i*0.25,
                                                           y+j*0.25, z+k*0.25) > 0);
                }
            }
        }
        ASSERT(exclusions > 0);
    }
    alea_destroy(sys);
}

TEST(selected_ray_batches_keep_material_thinner_than_the_forward_sample) {
    const char deck[] = "Two distinct planes\n10 0 -1\n20 1 -1 1 -2\n30 0 1 2\n\n"
                        "1 px 0\n2 p 1 1 0 0\n\nm1 1001 1\n";
    mcnp_model_t* model = mcnp_load_string(deck, sizeof(deck)-1);
    ASSERT_NOT_NULL(model);
    alea_bbox_t box = {-1024,1024,-1,1,-1,1};
    alea_prepared_queries_t* query = alea_prepared_queries_create(model->sys, &box);
    alea_raycast_result_t* prepared = alea_raycast_result_create();
    alea_raycast_batch_result_t* batch = alea_raycast_batch_result_create();
    ASSERT_NOT_NULL(query);
    ASSERT_NOT_NULL(prepared);
    ASSERT_NOT_NULL(batch);
    alea_raycast_batch_options_t options = {
        .struct_size = sizeof(options), .fields = ALEA_RAY_BATCH_MATERIAL | ALEA_RAY_BATCH_SURFACES
    };
    /* At y = -2^-39, x = 0 and x + y = 0 enclose this physical width. */
    const double width = 0x1p-39;
    for (int direction = -1; direction <= 1; direction += 2) {
        alea_segment_ray_t ray = {{-1024 * direction,-width,0}, {direction,0,0}, 2048};
        ASSERT_EQ(alea_prepared_trace_rays(query, &ray, 1, &prepared), 0);
        ASSERT_EQ(alea_raycast_hier_batch(model->sys, ray.origin, ray.direction, 1, ray.distance,
                                           &options, batch), 0);
        ASSERT_EQ(alea_raycast_segment_count(prepared), 3);
        ASSERT_EQ(alea_raycast_batch_segment_count(batch), 3);
        double lo, hi;
        int cell, material, entry_surface, exit_surface;
        ASSERT_EQ(alea_raycast_segment_get(prepared, 1, &lo, &hi, &cell, &material,
                                           NULL, &entry_surface, &exit_surface), 0);
        ASSERT_EQ(lo, direction > 0 ? 1024 : 1024-width);
        ASSERT_EQ(hi, direction > 0 ? 1024+width : 1024);
        ASSERT_EQ(cell, 20);
        ASSERT_EQ(material, 1);
        ASSERT_EQ(entry_surface, direction > 0 ? 1 : 2);
        ASSERT_EQ(exit_surface, direction > 0 ? 2 : 1);
        ASSERT_EQ(alea_raycast_batch_t_enter(batch)[1], lo);
        ASSERT_EQ(alea_raycast_batch_t_exit(batch)[1], hi);
        ASSERT_EQ(alea_raycast_batch_cell_ids(batch)[1], cell);
        ASSERT_EQ(alea_raycast_batch_material_ids(batch)[1], material);
        ASSERT_EQ(alea_raycast_batch_enter_surface_ids(batch)[1], entry_surface);
        ASSERT_EQ(alea_raycast_batch_exit_surface_ids(batch)[1], exit_surface);
    }
    alea_raycast_batch_result_destroy(batch);
    alea_raycast_result_destroy(prepared);
    alea_prepared_queries_destroy(query);
    mcnp_model_destroy(model);
}

TEST(cylindrical_void_at_a_translated_fill_boundary) {
    /* The outer radius is one binary64 step above 32. There is no world
     * x between the radii, but the ray enters the shell at t near 1 where
     * both crossings are distinct. The cylinders belong to different frames. */
    const char deck[] = "Cylindrical shell\n1 0 -2 fill=7 (4 0 0)\n2 0 2 fill=8\n"
        "11 1 -1 -1 u=7\n12 0 1 u=7\n21 1 -1 -3 u=8\n\n"
        "1 c/z -4 0 32\n2 cz 32.000000000000007\n3 so 1000\n\nm1 1001 1\n";
    mcnp_model_t* model = mcnp_load_string(deck, sizeof(deck) - 1);
    ASSERT_NOT_NULL(model);
    alea_bbox_t box = {28, 36, -2, 2, -2, 2};
    alea_prepared_queries_t* query = alea_prepared_queries_create(model->sys, &box);
    alea_raycast_result_t* prepared = alea_raycast_result_create();
    alea_raycast_batch_result_t* batch = alea_raycast_batch_result_create();
    ASSERT_NOT_NULL(query);
    ASSERT_NOT_NULL(prepared);
    ASSERT_NOT_NULL(batch);
    alea_raycast_batch_options_t options = {
        .struct_size = sizeof(options), .fields = ALEA_RAY_BATCH_MATERIAL | ALEA_RAY_BATCH_SURFACES |
                                                  ALEA_RAY_BATCH_RESOLUTION_FLAGS
    };
    double width = nextafter(32.0, INFINITY) - 32.0;
    /* The longer forward flight leaves just one interior ray distance. */
    const double offsets[] = {-1, 1, -16, 16};
    for (size_t k = 0; k < sizeof(offsets) / sizeof(offsets[0]); ++k) {
        double distance = fabs(offsets[k]);
        int sign = offsets[k] < 0 ? 1 : -1;
        alea_segment_ray_t ray = {{32 + offsets[k], 0, 0}, {sign, 0, 0}, 2 * distance};
        ASSERT_EQ(alea_prepared_trace_rays(query, &ray, 1, &prepared), 0);
        ASSERT_EQ(alea_raycast_hier_batch(model->sys, ray.origin, ray.direction, 1, ray.distance,
                                          &options, batch), 0);
        ASSERT_EQ(alea_raycast_segment_count(prepared), 3);
        ASSERT_EQ(alea_raycast_batch_segment_count(batch), 3);
        const int expected[] = {sign > 0 ? 11 : 21, 12, sign > 0 ? 21 : 11};
        for (int i = 0; i < 3; ++i) {
            int cell;
            ASSERT_EQ(alea_raycast_segment_get(prepared, i, NULL, NULL, &cell, NULL, NULL, NULL, NULL), 0);
            ASSERT_EQ(cell, expected[i]);
            ASSERT_EQ(alea_raycast_batch_cell_ids(batch)[i], expected[i]);
        }
        double lo, hi;
        int material, entry, exit;
        ASSERT_EQ(alea_raycast_segment_get(prepared, 1, &lo, &hi, NULL, &material, NULL, &entry, &exit), 0);
        ASSERT(hi > lo);
        ASSERT_NEAR(lo, sign > 0 ? distance : distance - width, distance * DBL_EPSILON);
        ASSERT_NEAR(hi, sign > 0 ? distance + width : distance, distance * DBL_EPSILON);
        ASSERT_EQ(entry, sign > 0 ? 1 : 2);
        ASSERT_EQ(exit, sign > 0 ? 2 : 1);
        ASSERT_EQ(material, 0);
        ASSERT_EQ(alea_raycast_batch_t_enter(batch)[1], lo);
        ASSERT_EQ(alea_raycast_batch_t_exit(batch)[1], hi);
        ASSERT_EQ(alea_raycast_batch_enter_surface_ids(batch)[1], entry);
        ASSERT_EQ(alea_raycast_batch_exit_surface_ids(batch)[1], exit);
        ASSERT_EQ(alea_raycast_batch_material_ids(batch)[1], 0);
        ASSERT_EQ(alea_raycast_batch_resolution_flags(batch)[1], 0);
    }
    alea_raycast_batch_result_destroy(batch);
    alea_raycast_result_destroy(prepared);
    alea_prepared_queries_destroy(query);
    mcnp_model_destroy(model);
}

TEST(adjacent_parallel_planes_at_a_fill_boundary) {
    /* The lower plane is the double immediately below 1708. The material
     * precedes the void in universe 8, so rounding a probe onto the upper
     * plane can select the material throughout the thin interval. */
    const char deck[] = "Adjacent placed planes\n1 0 -1 fill=7\n2 0 1 fill=8\n"
        "11 1 -1 -1 u=7\n21 1 -1 2 u=8\n12 0 -2 u=8\n\n"
        "1 pz 1707.999999999999800\n2 pz 1708\n\nm1 1001 1\n";
    mcnp_model_t* model = mcnp_load_string(deck, sizeof(deck) - 1);
    ASSERT_NOT_NULL(model);
    alea_bbox_t domain = {-1, 1, -1, 1, 1707, 1709};
    alea_prepared_queries_t* query = alea_prepared_queries_create(model->sys, &domain);
    ASSERT_NOT_NULL(query);
    alea_raycast_result_t* public = alea_raycast_result_create();
    alea_raycast_result_t* prepared = alea_raycast_result_create();
    ASSERT_NOT_NULL(public);
    ASSERT_NOT_NULL(prepared);
    const double origins[] = {0, 3416, -40000, 40000};
    for (size_t i = 0; i < sizeof(origins) / sizeof(origins[0]); ++i) {
        double z = origins[i], dz = z < 1708 ? 1 : -1;
        double a = fabs(nextafter(1708.0, -INFINITY) - z), b = fabs(1708.0 - z);
        alea_segment_ray_t ray = {{0, 0, z}, {0, 0, dz}, 50000};
        ASSERT_EQ(alea_raycast_hier_fast_segments(model->sys, 0, 0, z, 0, 0, dz, ray.distance, public), 0);
        ASSERT_EQ(alea_prepared_trace_rays(query, &ray, 1, &prepared), 0);
        alea_raycast_result_t* results[] = {public, prepared};
        for (int j = 0; j < 2; ++j) {
            size_t n = alea_raycast_segment_count(results[j]);
            ASSERT_EQ(n, a == b ? 2 : 3);
            for (size_t k = 0; k < n; ++k) {
                int cell, material, entry, exit;
                double lo, hi;
                ASSERT_EQ(alea_raycast_segment_get(results[j], k, &lo, &hi, &cell, &material,
                                                   NULL, &entry, &exit), 0);
                if (k == 0) {
                    ASSERT_EQ(cell, dz > 0 ? 11 : 21);
                    ASSERT(hi == fmin(a, b));
                } else if (k == n - 1) {
                    ASSERT_EQ(cell, dz > 0 ? 21 : 11);
                    ASSERT(lo == fmax(a, b));
                } else {
                    ASSERT_EQ(cell, 12);
                    ASSERT_EQ(material, 0);
                    ASSERT_EQ(entry, dz > 0 ? 1 : 2);
                    ASSERT_EQ(exit, dz > 0 ? 2 : 1);
                    ASSERT(lo == fmin(a, b) && hi == nextafter(lo, INFINITY));
                }
            }
        }
    }
    alea_raycast_result_destroy(public);
    alea_raycast_result_destroy(prepared);
    alea_prepared_queries_destroy(query);
    mcnp_model_destroy(model);
}

TEST(adjacent_internal_planes_select_the_nested_fill_boundary) {
    /* Cell 11 lies below plane 1. Its redundant plane 2 sits one step
     * above it. A rounded midpoint on plane 1 can extend cell 11 to 2. */
    const char deck[] = "Adjacent planes in a translated fill\n"
        "1 0 -3 fill=7 (1 0 0)\n11 1 -1 -1 -2 u=7\n12 0 1 u=7 fill=8\n"
        "13 0 -3 u=8 fill=9\n14 2 -1 -3 u=9\n\n"
        "1 pz 1708\n2 pz 1708.000000000000227\n3 so 10000\n\n"
        "m1 1001 1\nm2 1001 1\n";
    mcnp_model_t* model = mcnp_load_string(deck, sizeof(deck) - 1);
    ASSERT_NOT_NULL(model);
    alea_bbox_t domain = {-2, 2, -2, 2, -1, 3417};
    alea_prepared_queries_t* query = alea_prepared_queries_create(model->sys, &domain);
    ASSERT_NOT_NULL(query);
    alea_raycast_result_t* result = alea_raycast_result_create();
    ASSERT_NOT_NULL(result);
    for (int reverse = 0; reverse < 2; ++reverse) {
        alea_segment_ray_t ray = {{0, 0, reverse ? 3416 : 0}, {0, 0, reverse ? -1 : 1}, 2000};
        ASSERT_EQ(alea_prepared_trace_rays(query, &ray, 1, &result), 0);
        ASSERT_EQ(alea_raycast_segment_count(result), 2);
        int cell, surface;
        double t;
        ASSERT_EQ(alea_raycast_segment_get(result, 0, NULL, &t, &cell, NULL, NULL, NULL, &surface), 0);
        ASSERT_EQ(cell, reverse ? 14 : 11);
        ASSERT_EQ(surface, 1);
        ASSERT(t == 1708);
        ASSERT_EQ(alea_raycast_segment_get(result, 1, &t, NULL, &cell, NULL, NULL, &surface, NULL), 0);
        ASSERT_EQ(cell, reverse ? 11 : 14);
        ASSERT_EQ(surface, 1);
        ASSERT(t == 1708);
    }
    alea_raycast_result_destroy(result);
    alea_prepared_queries_destroy(query);
    mcnp_model_destroy(model);
}

TEST(lattice_plane_crossings_preserve_void_and_undefined_fill) {
    for (int explicit_void = 0; explicit_void <= 1; ++explicit_void) {
        char deck[320];
        int n = snprintf(deck, sizeof(deck), "Thin interval inside a repeating element\n"
            "1 1 -1 -1\n2 0 1 fill=7 (4 0 0)\n10 0 -3 lat=1 u=7 fill=8\n"
            "21 1 -1 2 u=8\n%s\n1 pz 32\n2 pz 32.000000000000007105\n"
            "3 rpp -100 100 -100 100 -100 100\n\nm1 1001 1\n", explicit_void ? "22 0 -2 u=8\n" : "");
        ASSERT(n > 0 && (size_t)n < sizeof(deck));
        mcnp_model_t* model = mcnp_load_string(deck, (size_t)n);
        ASSERT_NOT_NULL(model);
        alea_bbox_t domain = {-5, 5, -5, 5, -1, 65};
        alea_prepared_queries_t* query = alea_prepared_queries_create(model->sys, &domain);
        ASSERT_NOT_NULL(query);
        alea_raycast_result_t* result = alea_raycast_result_create();
        ASSERT_NOT_NULL(result);
        for (int reverse = 0; reverse < 2; ++reverse) {
            double z = reverse ? 64 : 0, dz = reverse ? -1 : 1;
            alea_segment_ray_t ray = {{0, 0, z}, {0, 0, dz}, 64};
            for (int prepared = 0; prepared <= 1; ++prepared) {
                int rc = prepared ? alea_prepared_trace_rays(query, &ray, 1, &result) :
                    alea_raycast_hier_fast_segments(model->sys, 0, 0, z, 0, 0, dz, 64, result);
                ASSERT_EQ(rc, 0);
                ASSERT_EQ(alea_raycast_segment_count(result), 3);
                ASSERT_EQ(result->segments.data[0].cell_id, reverse ? 21 : 1);
                ASSERT_EQ(result->segments.data[2].cell_id, reverse ? 1 : 21);
                const alea_ray_segment_t* middle = &result->segments.data[1];
                ASSERT_EQ(middle->cell_id, explicit_void ? 22 : 10);
                ASSERT_EQ(middle->material_id, 0);
                ASSERT_EQ(middle->resolution_flags, explicit_void ? 0 : ALEA_RESOLVE_UNDEFINED_FILL);
                ASSERT_EQ(middle->enter_surface_id, reverse ? 2 : 1);
                ASSERT_EQ(middle->exit_surface_id, reverse ? 1 : 2);
                double upper = nextafter(32.0, INFINITY);
                ASSERT(middle->t_enter == (reverse ? 64 - upper : 32));
                ASSERT(middle->t_exit == (reverse ? 32 : upper));
            }
        }
        alea_raycast_result_destroy(result);
        alea_prepared_queries_destroy(query);
        mcnp_model_destroy(model);
    }
}

TEST(one_plane_in_parent_and_child_placements) {
    /* A single plane card cuts the child at z = 32 - width and the
     * parent at z = 32. With zero displacement the cuts coincide. */
    for (int shifted = 0; shifted <= 1; ++shifted) {
        double width = shifted ? 0x1p-47 : 0;
        char deck[256];
        int n = snprintf(deck, sizeof(deck),
            "Placed plane\n1 0 -1 fill=7 (4 0 %.17g)\n2 0 1 fill=8\n"
            "11 1 -1 -1 u=7\n12 0 1 u=7\n21 1 -1 -3 u=8\n\n"
            "1 pz 32\n3 so 1000\n\nm1 1001 1\n", -width);
        ASSERT(n > 0 && (size_t)n < sizeof(deck));
        mcnp_model_t* model = mcnp_load_string(deck, (size_t)n);
        ASSERT_NOT_NULL(model);
        alea_bbox_t box = {-2, 2, -2, 2, 28, 36};
        alea_prepared_queries_t* query = alea_prepared_queries_create(model->sys, &box);
        ASSERT_NOT_NULL(query);
        alea_raycast_result_t* results[] = {alea_raycast_result_create(), alea_raycast_result_create()};
        ASSERT_NOT_NULL(results[0]);
        ASSERT_NOT_NULL(results[1]);
        for (int dz = -1; dz <= 1; dz += 2) {
            alea_segment_ray_t ray = {{0, 0, 32 - dz}, {0, 0, dz}, 2};
            ASSERT_EQ(alea_raycast_hier_fast_segments(model->sys, 0, 0, ray.origin[2],
                                                      0, 0, dz, 2, results[0]), 0);
            ASSERT_EQ(alea_prepared_trace_rays(query, &ray, 1, &results[1]), 0);
            for (int k = 0; k < 2; ++k) {
                ASSERT_EQ(alea_raycast_segment_count(results[k]), 2 + shifted);
                for (int i = 0; i < 2 + shifted; ++i) {
                    int cell, material, enter, exit;
                    double lo, hi;
                    ASSERT_EQ(alea_raycast_segment_get(results[k], i, &lo, &hi, &cell,
                                                       &material, NULL, &enter, &exit), 0);
                    int first = dz > 0 ? 11 : 21;
                    int last = dz > 0 ? 21 : 11;
                    ASSERT_EQ(cell, i == 0 ? first : i == 1 + shifted ? last : 12);
                    if (shifted && i == 1) {
                        ASSERT_EQ(material, 0);
                        ASSERT_EQ(enter, 1);
                        ASSERT_EQ(exit, 1);
                        ASSERT(lo == (dz > 0 ? 1 - width : 1));
                        ASSERT(hi == (dz > 0 ? 1 : 1 + width));
                    }
                }
            }
        }
        alea_raycast_result_destroy(results[0]);
        alea_raycast_result_destroy(results[1]);
        alea_prepared_queries_destroy(query);
        mcnp_model_destroy(model);
    }
}

TEST(regional_strict_verifies_an_adjacent_parent_and_child_crossing) {
    for (int overlap = 0; overlap <= 1; ++overlap) {
        char deck[384];
        double upper = nextafter(32.0, INFINITY);
        int n = snprintf(deck, sizeof(deck), "Adjacent parent and child boundaries\n"
            "1 1 -1 -1 -3\n2 0 1 -3 fill=7\n11 2 -1 -2 u=7\n12 3 -1 2 u=7\n"
            "%s99 0 3 imp:n=0\n\n1 pz 32\n2 pz %.17g\n3 so 1000\n\n"
            "m1 1001 1\nm2 1001 1\nm3 1001 1\n", overlap ? "13 1 -1 -2 u=7\n" : "", upper);
        ASSERT(n > 0 && (size_t)n < sizeof(deck));
        mcnp_model_t* model = mcnp_load_string(deck, (size_t)n);
        ASSERT_NOT_NULL(model);
        alea_bbox_t box = {-5, 5, -5, 5, 0, 64};
        alea_prepared_queries_t* query = alea_prepared_queries_create(model->sys, &box);
        alea_ray_navigator_t* nav = alea_ray_navigator_create(model->sys);
        ASSERT_NOT_NULL(query);
        ASSERT_NOT_NULL(nav);
        ASSERT_EQ(alea_ray_navigator_use_regional(nav, query), 0);
        double origin[] = {0, 0, 0}, direction[] = {0, 0, 1};
        alea_nav_location_t loc;
        alea_nav_event_t e;
        ASSERT_EQ(alea_ray_navigator_restart(nav, origin, direction, &loc), 0);
        ASSERT_EQ(alea_ray_navigator_advance(nav, INFINITY, 50, &e), 0);
        ASSERT_EQ(e.kind, ALEA_NAV_BOUNDARY);
        ASSERT_EQ(e.surface_id, 1);
        ASSERT(e.distance == 32.0);
        if (overlap) {
            ASSERT(e.after.kind == ALEA_NAV_UNRESOLVED || e.after.kind == ALEA_NAV_OVERLAP);
            ASSERT(alea_ray_navigator_advance(nav, INFINITY, 1, &e) < 0);
        } else {
            ASSERT_EQ(e.after.kind, ALEA_NAV_MATERIAL);
            ASSERT_EQ(e.after.cell_id, 11);
            ASSERT_EQ(alea_ray_navigator_advance(nav, INFINITY, 1, &e), 0);
            ASSERT(e.distance == upper - 32.0);
            ASSERT_EQ(e.surface_id, 2);
            ASSERT_EQ(e.after.cell_id, 12);
        }
        alea_ray_navigator_destroy(nav);
        alea_prepared_queries_destroy(query);
        mcnp_model_destroy(model);
    }
}

TEST(transport_recovery_preserves_surfaces_collisions_and_owner) {
    for (int regional = 0; regional <= 1; ++regional) {
        const char* deck = "Thin gap\n1 0 -3 fill=7\n2 0 3\n"
            "11 1 -1 -1 u=7\n12 2 -1 2 u=7\n\n1 pz 32\n2 pz 32.00000000005\n"
            "3 so 1000\n\nm1 1001 1\nm2 1001 1\n";
        mcnp_model_t* m = mcnp_load_string(deck, strlen(deck));
        ASSERT_NOT_NULL(m);
        alea_ray_navigator_t* nav = alea_ray_navigator_create(m->sys);
        ASSERT_NOT_NULL(nav);
        if (regional) {
            alea_bbox_t box = {-5,5,-5,5,0,64};
            alea_prepared_queries_t* idx = alea_prepared_queries_create(m->sys, &box);
            ASSERT_NOT_NULL(idx);
            ASSERT_EQ(alea_ray_navigator_use_regional(nav, idx), 0);
            alea_prepared_queries_destroy(idx);
        }
        const double origin[] = {0,0,0}, forward[] = {0,0,1}, backward[] = {0,0,-1};
        alea_nav_location_t loc;
        alea_nav_event_t e;
        ASSERT_EQ(alea_ray_navigator_restart(nav, origin, forward, &loc), 0);
        ASSERT_EQ(alea_ray_navigator_advance(nav, INFINITY, 50, &e), 0);
        ASSERT_EQ(e.after.kind, ALEA_NAV_UNDEFINED_FILL);
        ASSERT_EQ(alea_ray_navigator_set_validation_mode(nav, ALEA_NAV_VALIDATE_TRANSPORT), 0);
        ASSERT_EQ(alea_ray_navigator_restart(nav, origin, forward, &loc), 0);
        ASSERT_EQ(alea_ray_navigator_advance(nav, INFINITY, 50, &e), 0);
        ASSERT_EQ(e.after.cell_id, 11);
        ASSERT_EQ(e.surface_id, 1);
        ASSERT_EQ(alea_ray_navigator_advance(nav, INFINITY, 1e-11, &e), 0);
        ASSERT_EQ(e.kind, ALEA_NAV_DISTANCE_LIMIT);
        ASSERT_EQ(alea_ray_navigator_advance(nav, 1e-11, 1, &e), 0);
        ASSERT_EQ(e.kind, ALEA_NAV_COLLISION);
        ASSERT_EQ(e.before.cell_id, 11);
        alea_ray_navigator_t* child = alea_ray_navigator_clone(nav);
        ASSERT_NOT_NULL(child);
        ASSERT_EQ(alea_ray_navigator_set_direction(child, backward), 0);
        ASSERT_EQ(alea_ray_navigator_advance(child, INFINITY, 1, &e), 0);
        ASSERT_EQ(e.surface_id, 1);
        ASSERT_EQ(e.after.cell_id, 11);
        /* A banked secondary starts with a fresh ray and the parent's
         * verified recovery context, just as the fixed-source driver does. */
        alea_transport_recovery_seed_t seed;
        alea_nav_get_recovery_seed(nav, &seed);
        double birth[] = {0, 0, 32.0 + 2e-11};
        ASSERT_EQ(alea_ray_navigator_restart(child, birth, forward, &loc), 0);
        ASSERT_EQ(alea_nav_resume_recovery(child, &seed, &loc), 0);
        ASSERT_EQ(loc.cell_id, 11);
        ASSERT_EQ(alea_ray_navigator_advance(child, INFINITY, 1, &e), 0);
        ASSERT_EQ(e.after.cell_id, 12);
        alea_ray_navigator_destroy(child);

        ASSERT_EQ(alea_ray_navigator_advance(nav, INFINITY, 1, &e), 0);
        ASSERT_EQ(e.surface_id, 2);
        ASSERT_EQ(e.after.cell_id, 12);
        alea_nav_recovery_stats_t stats;
        ASSERT_EQ(alea_ray_navigator_recovery_stats(nav, &stats), 0);
        ASSERT_EQ(stats.gaps, 1u);
        ASSERT_EQ(stats.overlaps, 0u);
        ASSERT_NEAR(stats.distance_cm, 32.00000000005 - 32.0, 1e-24);
        ASSERT_EQ(stats.last_enter_surface_id, 1);
        ASSERT_EQ(stats.last_exit_surface_id, 2);
        double inside[] = {0,0,32.000000000025};
        ASSERT_EQ(alea_ray_navigator_restart(nav, inside, forward, &loc), 0);
        ASSERT_EQ(loc.kind, ALEA_NAV_UNDEFINED_FILL);
        alea_ray_navigator_destroy(nav);
        mcnp_model_destroy(m);
    }
}

TEST(transport_recovery_requires_a_unique_transmissive_destination) {
    /* The union exposes both physical crossings on the hierarchy path.
     * The incident cell remains a claimant on the overlap. */
    const char* base = "Recovery eligibility\n1 0 -3 fill=7\n2 0 3\n"
        "11 1 -1 %s u=7\n12 2 -1 2 u=7\n%s\n1 pz 32\n%s2 pz %.17g\n"
        "3 so 1000\n\nm1 1001 1\nm2 1001 1\n";
    enum { OVERLAP, THIRD, THIN_LAYER, MIRROR, TOO_WIDE };
    for (int k = OVERLAP; k <= TOO_WIDE; ++k) {
        bool overlap = k == OVERLAP || k == THIRD;
        const char* extra = k == THIRD ? "13 0 2 -1 u=7\n" :
                            k == THIN_LAYER ? "13 0 1 -2 u=7\n" : "";
        double upper = overlap ? 31.99999999995 : k == TOO_WIDE ? 32.0000000002 : 32.00000000005;
        char deck[512];
        int n = snprintf(deck, sizeof(deck), base, overlap ? "(-1 -2):(-1 2)" : "-1",
                         extra, k == MIRROR ? "*" : "", upper);
        ASSERT(n > 0 && (size_t)n < sizeof(deck));
        mcnp_model_t* m = mcnp_load_string(deck, (size_t)n);
        ASSERT_NOT_NULL(m);
        alea_ray_navigator_t* nav = alea_ray_navigator_create(m->sys);
        ASSERT_NOT_NULL(nav);
        ASSERT_EQ(alea_ray_navigator_set_validation_mode(nav, ALEA_NAV_VALIDATE_TRANSPORT), 0);
        double origin[] = {0,0,0}, direction[] = {0,0,1};
        alea_nav_location_t loc;
        alea_nav_event_t e;
        ASSERT_EQ(alea_ray_navigator_restart(nav, origin, direction, &loc), 0);
        ASSERT_EQ(alea_ray_navigator_advance(nav, INFINITY, 50, &e), 0);
        if (k == OVERLAP) {
            ASSERT_EQ(e.after.cell_id, 11);
            ASSERT_EQ(alea_ray_navigator_advance(nav, INFINITY, 1, &e), 0);
            ASSERT_EQ(e.after.cell_id, 12);
        } else if (k == THIN_LAYER) {
            ASSERT_EQ(e.after.cell_id, 13);
            ASSERT_EQ(e.after.kind, ALEA_NAV_VOID);
        } else {
            ASSERT(e.after.kind == ALEA_NAV_OVERLAP || e.after.kind == ALEA_NAV_UNDEFINED_FILL);
            ASSERT(alea_ray_navigator_advance(nav, INFINITY, 1, &e) < 0);
        }
        alea_nav_recovery_stats_t stats;
        ASSERT_EQ(alea_ray_navigator_recovery_stats(nav, &stats), 0);
        ASSERT_EQ(stats.gaps, 0u);
        ASSERT_EQ(stats.overlaps, k == OVERLAP ? 1u : 0u);
        alea_ray_navigator_destroy(nav);
        mcnp_model_destroy(m);
    }
}

TEST(regional_particle_navigator_keeps_windows_internal_and_checks_edits) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);
    int material = alea_add_material(sys, 17);
    int surf = alea_sphere_surface(sys, 1, 0, 0, 0, 128);
    ASSERT(surf >= 0);
    ASSERT(alea_add_cell(sys, 29, alea_halfspace(sys, surf, -1), material, -1, 0) >= 0);
    ASSERT_EQ(alea_surface_set_boundary(sys, 1, ALEA_BOUNDARY_VACUUM), 0);
    alea_bbox_t domain = {-200,200,-200,200,-200,200};
    alea_prepared_queries_t* idx = alea_prepared_queries_create(sys, &domain);
    alea_ray_navigator_t* nav = alea_ray_navigator_create(sys);
    ASSERT_NOT_NULL(idx);
    ASSERT_NOT_NULL(nav);
    ASSERT_EQ(alea_ray_navigator_use_regional(nav, idx), 0);
    alea_prepared_queries_destroy(idx);
    double origin[] = {0,0,0}, direction[] = {1,0,0};
    alea_nav_location_t source;
    ASSERT_EQ(alea_ray_navigator_restart(nav, origin, direction, &source), 0);
    ASSERT_EQ(source.cell_id, 29);
    alea_nav_event_t event;
    ASSERT_EQ(alea_ray_navigator_advance(nav, INFINITY, 64, &event), 0);
    ASSERT_EQ(event.kind, ALEA_NAV_DISTANCE_LIMIT);
    ASSERT_EQ(event.distance, 64);
    alea_ray_navigator_t* clone = alea_ray_navigator_clone(nav);
    ASSERT_NOT_NULL(clone);
    alea_ray_navigator_destroy(nav);
    ASSERT_EQ(alea_ray_navigator_advance(clone, 64, 128, &event), 0);
    ASSERT_EQ(event.kind, ALEA_NAV_VACUUM);
    ASSERT_EQ(event.distance, 64);
    ASSERT_EQ(event.surface_id, 1);
    ASSERT_EQ(alea_ray_navigator_restart(clone, origin, direction, &source), 0);
    ASSERT_EQ(alea_cell_set_density(sys, 0, -3), 0);
    ASSERT(alea_ray_navigator_advance(clone, 1, 2, &event) < 0);
    ASSERT(alea_ray_navigator_restart(clone, origin, direction, &source) < 0);
    alea_ray_navigator_destroy(clone);
    alea_destroy(sys);
}

TEST(prepared_point_and_ray_batches_follow_model_generations) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);
    int material = alea_add_material(sys, 17);
    alea_node_id_t ball = sphere(sys, 0);
    ASSERT(alea_add_cell(sys, 29, ball, material, -1, 0) >= 0);
    alea_bbox_t domain = {-2, 2, -2, 2, -2, 2};
    alea_prepared_queries_t* query = alea_prepared_queries_create(sys, &domain);
    ASSERT_NOT_NULL(query);
    const alea_query_point_t points[] = {{0,0,0}, {12,0,0}};
    alea_point_owner_t owners[2];
    ASSERT_EQ(alea_prepared_locate_points(query, points, 2, owners), 0);
    ASSERT_EQ(owners[0].status, 0);
    ASSERT_EQ(owners[0].cell_id, 29);
    ASSERT_EQ(owners[0].material_id, 17);
    ASSERT_EQ(owners[1].status, -1);
    ASSERT_EQ(alea_prepared_locate_points(query, NULL, 0, NULL), 0);

    alea_segment_ray_t ray = {{-3,0,0}, {2,0,0}, 6};
    alea_raycast_result_t* result = alea_raycast_result_create();
    ASSERT_NOT_NULL(result);
    ASSERT_EQ(alea_prepared_trace_rays(query, &ray, 1, &result), 0);
    ASSERT_EQ(alea_raycast_segment_count(result), 3);
    double entry, exit;
    int cell, mat;
    ASSERT_EQ(alea_raycast_segment_get(result, 1, &entry, &exit, &cell, &mat, NULL, NULL, NULL), 0);
    ASSERT_EQ(entry, 2.0);
    ASSERT_EQ(exit, 4.0);
    ASSERT_EQ(cell, 29);
    ASSERT_EQ(mat, 17);

    const double scales[] = {DBL_MIN, DBL_MAX};
    for (size_t i = 0; i < 2; ++i) {
        ray.direction[0] = scales[i];
        ASSERT_EQ(alea_prepared_trace_rays(query, &ray, 1, &result), 0);
        ASSERT_EQ(alea_raycast_segment_get(result, 1, &entry, &exit, &cell, &mat, NULL, NULL, NULL), 0);
        ASSERT_EQ(entry, 2.0);
        ASSERT_EQ(exit, 4.0);
        ASSERT_EQ(cell, 29);
    }
    ASSERT_EQ(alea_renumber_cells(sys, 59), 0);
    ASSERT_EQ(alea_prepared_locate_points(query, points, 2, owners), -1);
    alea_prepared_queries_destroy(query);
    query = alea_prepared_queries_create(sys, &domain);
    ASSERT_NOT_NULL(query);
    ASSERT_EQ(alea_prepared_locate_points(query, points, 2, owners), 0);
    ASSERT_EQ(owners[0].cell_id, 59);
    ASSERT_EQ(alea_cell_set_density(sys, 0, -3), 0);
    ASSERT_EQ(alea_prepared_trace_rays(query, &ray, 1, &result), -1);
    alea_prepared_queries_destroy(query);
    query = alea_prepared_queries_create(sys, &domain);
    ASSERT_NOT_NULL(query);
    ASSERT_EQ(alea_cell_set_region(sys, 0, sphere(sys, 8)), 0);
    ASSERT_EQ(alea_prepared_locate_points(query, points, 2, owners), -1);
    alea_prepared_queries_destroy(query);
    alea_raycast_result_destroy(result);
    alea_destroy(sys);
}

TEST(prepared_queries_report_invalid_arguments) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);
    ASSERT(alea_add_cell(sys, 1, sphere(sys, 0), ALEA_MATERIAL_VOID, 0, 0) >= 0);
    alea_bbox_t invalid = {2, -2, -2, 2, -2, 2};
    alea_error_clear();
    ASSERT_NULL(alea_prepared_queries_create(sys, &invalid));
    ASSERT_EQ(alea_error_code(), ALEA_ERR_INVALID_ARG);
    ASSERT_NULL(alea_prepared_queries_create(sys, NULL));
    ASSERT_EQ(alea_error_code(), ALEA_ERR_NULL_ARG);

    alea_bbox_t box = {-2, 2, -2, 2, -2, 2};
    alea_prepared_queries_t* query = alea_prepared_queries_create(sys, &box);
    ASSERT_NOT_NULL(query);
    alea_query_point_t point = {NAN, 0, 0};
    alea_point_owner_t owner;
    alea_error_clear();
    ASSERT_EQ(alea_prepared_locate_points(query, &point, 1, &owner), -1);
    ASSERT_EQ(alea_error_code(), ALEA_ERR_INVALID_ARG);
    ASSERT_EQ(alea_prepared_locate_points(query, NULL, 1, &owner), -1);
    ASSERT_EQ(alea_error_code(), ALEA_ERR_NULL_ARG);

    alea_segment_ray_t ray = {{0,0,0}, {0,0,0}, 1};
    alea_raycast_result_t* result = alea_raycast_result_create();
    ASSERT_NOT_NULL(result);
    alea_error_clear();
    ASSERT_EQ(alea_prepared_trace_rays(query, &ray, 1, &result), -1);
    ASSERT_EQ(alea_error_code(), ALEA_ERR_INVALID_ARG);
    ray.direction[0] = 1;
    ray.distance = INFINITY;
    alea_error_clear();
    ASSERT_EQ(alea_prepared_trace_rays(query, &ray, 1, &result), -1);
    ASSERT_EQ(alea_error_code(), ALEA_ERR_INVALID_ARG);
    alea_raycast_result_destroy(result);
    alea_prepared_queries_destroy(query);
    alea_destroy(sys);
}

TEST(point_lookup_prepares_only_the_active_point_backend) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);
    ASSERT_EQ(alea_add_cell(sys, 12, sphere(sys, 0), ALEA_MATERIAL_VOID, 0, 0), 0);
    ASSERT_EQ(alea_prepare_point_acceleration(sys), 0);
    ASSERT_EQ(alea_selected_query_stats(sys, NULL), 1);
    ASSERT_EQ(alea_prepare_point_lookup(sys), 0);
    ASSERT(sys->hier_spatial_index == NULL);
    ASSERT_EQ(alea_find_cell(sys, 0, 0, 0), 0);

    alea_clear_selected_queries(sys);
    ASSERT_EQ(alea_prepare_point_lookup(sys), 0);
    ASSERT_NOT_NULL(sys->hier_spatial_index);
    ASSERT_EQ(alea_find_cell(sys, 0, 0, 0), 0);
    ASSERT(!alea_system_query_cache_ready(sys, ALEA_CACHE_SURFACE_BVH));
    ASSERT(!alea_system_query_cache_ready(sys, ALEA_CACHE_RECOVERY_SURFACES));
    ASSERT_EQ(alea_prepare_query_acceleration(sys), 0);
    ASSERT(alea_system_query_cache_ready(sys, ALEA_CACHE_SURFACE_BVH));
    alea_destroy(sys);
}

TEST(root_point_index_shares_fill_transforms_without_building_spatial_trees) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);
    int envelope = alea_box_surface(sys, 0, -10, 10, -3, 3, -3, 3);
    int parent = alea_add_cell(sys, 11, alea_halfspace(sys, envelope, -1),
                               ALEA_MATERIAL_VOID, 0, 0);
    ASSERT(parent >= 0);
    int material = alea_add_material(sys, 2);
    int child = alea_add_cell(sys, 22, sphere(sys, 0), material, -1, 7);
    ASSERT(child >= 0);
    double shift[] = {4, 0, 0};
    ASSERT_EQ(alea_add_transform(sys, 1, shift, 3, 0), 0);
    ASSERT_EQ(alea_set_fill(sys, parent, 7, 1), 0);
    ASSERT_EQ(alea_system_prepare_query_caches(sys, ALEA_CACHE_FILL_MATRICES), 0);
    const alea_matrix_t* frame = alea_cell_fill_matrix(sys, parent);
    ASSERT_NOT_NULL(frame);
    ASSERT_EQ(alea_prepare_root_points(sys, NULL), 0);
    ASSERT_EQ(alea_selected_query_stats(sys, NULL), 1);
    ASSERT_EQ(alea_find_cell(sys, 4, 0, 0), child);
    ASSERT_EQ(alea_find_cell(sys, 0, 0, 0), parent);
    uint64_t visits;
    size_t universe = (size_t)(alea_get_universe(sys, 7) - sys->universes.data);
    ASSERT_EQ(alea_root_point_scan_tests(sys, universe, &visits, NULL), 2);
    ASSERT_EQ(visits, 2);
    ASSERT(sys->hier_spatial_index == NULL);
    ASSERT_EQ(alea_prepare_query_acceleration(sys), 0);
    ASSERT(alea_cell_fill_matrix(sys, parent) == frame);
    ASSERT_EQ(alea_set_fill(sys, parent, 7, 0), 0);
    ASSERT(alea_cell_fill_matrix(sys, parent) == NULL);
    ASSERT_EQ(alea_prepare_root_points(sys, NULL), 0);
    ASSERT_EQ(alea_find_cell(sys, 0, 0, 0), child);
    ASSERT_EQ(alea_find_cell(sys, 4, 0, 0), parent);
    ASSERT(sys->hier_spatial_index == NULL);
    alea_destroy(sys);
}

TEST(root_point_parallel_bounds_preserve_definition_order) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);
    for (int i = 0; i < 160; ++i) {
        /* Neighboring spheres overlap at half-integer x coordinates. */
        ASSERT_EQ(alea_add_cell(sys, 200 - i, sphere(sys, i),
                                ALEA_MATERIAL_VOID, 0, 0), i);
    }
    ASSERT_EQ(alea_prepare_root_points(sys, NULL), 0);
    for (int i = 0; i < 159; ++i)
        ASSERT_EQ(alea_find_cell(sys, i + 0.5, 0, 0), i);
    ASSERT(sys->hier_spatial_index == NULL);
    alea_destroy(sys);
}

TEST(child_point_indexes_preserve_local_ownership_and_uncovered_fills) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);
    const int v = ALEA_MATERIAL_VOID;
    int box = alea_box_surface(sys, 0, -4, 4, -2, 2, -2, 2);
    ASSERT_EQ(alea_add_cell(sys, 10, alea_halfspace(sys, box, -1), v, 0, 0), 0);
    ASSERT_EQ(alea_set_fill(sys, 0, 7, 0), 0);
    ASSERT_EQ(alea_add_cell(sys, 90, sphere(sys, 0), v, 0, 7), 1);
    ASSERT_EQ(alea_add_cell(sys, 20, sphere(sys, 0), v, 0, 7), 2);
    /* The winning cell is filled. Its child lies outside the indexed local box. */
    double shift[] = {50, 0, 0};
    ASSERT_EQ(alea_add_transform(sys, 1, shift, 3, 0), 0);
    ASSERT_EQ(alea_set_fill(sys, 1, 8, 1), 0);
    ASSERT_EQ(alea_add_cell(sys, 30, sphere(sys, -50), v, 0, 8), 3);
    const int ids[] = {7, 8};
    interval_point_children_t children = {ids, 2, CHILD_POINTS_BOUNDS};
    const interval_child_mode_t modes[] = {CHILD_POINTS_BOUNDS, CHILD_POINTS_REGIONS, CHILD_POINTS_ADAPTIVE};
    for (size_t i = 0; i < sizeof(modes) / sizeof(modes[0]); ++i) {
        children.mode = modes[i];
        ASSERT_EQ(alea_prepare_root_points(sys, &children), 0);
        ASSERT_EQ(alea_prepare_root_points(sys, &children), 0);
        ASSERT_EQ(alea_find_cell(sys, 0, 0, 0), 3);
        ASSERT_EQ(alea_find_cell(sys, 3, 0, 0), 0);
        ASSERT_EQ(alea_find_cell(sys, 10, 0, 0), -1);
        uint64_t visits;
        size_t u = (size_t)(alea_get_universe(sys, 7) - sys->universes.data);
        ASSERT_EQ(alea_root_point_scan_tests(sys, u, &visits, NULL), children.mode == CHILD_POINTS_BOUNDS ? 1 : 0);
        ASSERT_EQ(visits, children.mode == CHILD_POINTS_BOUNDS ? 2 : 0);
        if (children.mode == CHILD_POINTS_ADAPTIVE) {
            const interval_point_children_t different = {&ids[1], 1, CHILD_POINTS_ADAPTIVE};
            ASSERT_EQ(alea_prepare_root_points(sys, &different), -1);
        }
        ASSERT_EQ(alea_prepare_root_points(sys, NULL), -1);
        ASSERT(sys->hier_spatial_index == NULL);
        alea_clear_selected_queries(sys);
    }
    const int unknown = 42;
    children = (interval_point_children_t){&unknown, 1, CHILD_POINTS_BOUNDS};
    ASSERT_EQ(alea_prepare_root_points(sys, &children), -1);
    ASSERT_EQ(alea_selected_query_stats(sys, NULL), 0);
    children = (interval_point_children_t){ids, 2, CHILD_POINTS_REGIONS};
    ASSERT_EQ(alea_prepare_root_points(sys, &children), 0);
    ASSERT_EQ(alea_set_fill(sys, 1, 8, 0), 0);
    ASSERT_EQ(alea_selected_query_stats(sys, NULL), 0);
    ASSERT_EQ(alea_prepare_root_points(sys, &children), 0);
    ASSERT_EQ(alea_find_cell(sys, 0, 0, 0), 1);
    alea_destroy(sys);
}

static int query_adaptive_child(void* context, size_t worker, size_t begin, size_t end) {
    (void)worker;
    alea_system_t* sys = context;
    const double xs[] = {4, 0, 7, 22, 50};
    const int owners[] = {521, 2, 0, 520, -1};
    for (size_t i = begin; i < end; ++i) {
        size_t point = i % 5;
        if (alea_find_cell(sys, xs[point], 0, 0) != owners[point]) return -1;
        if (i % 64 == 0) {
            alea_selected_query_stats_t stats;
            if (alea_selected_query_stats(sys, &stats) != 1 || !stats.bytes) return -1;
        }
    }
    return 0;
}

static int promote_or_read(void* context, size_t worker, size_t begin, size_t end) {
    for (size_t i = begin; i < end; ++i) {
        if (i < 2) {
            int rc = alea_promote_point_child(context, 7);
            if (rc != 0 && rc != 1) return -1;
        }
        if (query_adaptive_child(context, worker, 0, 2000)) return -1;
    }
    return 0;
}

TEST(adaptive_child_bounds_publish_in_order_and_retry_after_allocation_failure) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);
    int world = alea_box_surface(sys, 0, -4, 8, -2, 2, -2, 2);
    int placed = alea_box_surface(sys, 0, 20, 24, -2, 2, -2, 2);
    ASSERT_EQ(alea_add_cell(sys, 1, alea_halfspace(sys, world, -1), ALEA_MATERIAL_VOID, 0, 0), 0);
    ASSERT_EQ(alea_add_cell(sys, 2, alea_halfspace(sys, placed, -1), ALEA_MATERIAL_VOID, 0, 0), 1);
    double shift[] = {120, 0, 0};
    ASSERT_EQ(alea_add_transform(sys, 1, shift, 3, 0), 0);
    ASSERT_EQ(alea_set_fill(sys, 0, 7, 0), 0);
    ASSERT_EQ(alea_set_fill(sys, 1, 7, 1), 0);
    alea_node_id_t first = sphere(sys, 0), last = sphere(sys, 4), distant = sphere(sys, -98);
    for (int i = 0; i < 521; ++i) {
        alea_node_id_t region = i < 518 ? first : i == 518 ? distant : last;
        ASSERT_EQ(alea_add_cell(sys, 1000 - i, region, ALEA_MATERIAL_VOID, 0, 7), i + 2);
    }
    ASSERT_EQ(alea_prepare_point_acceleration(sys), 0);
    size_t u = (size_t)(alea_get_universe(sys, 7) - sys->universes.data);
    ASSERT_EQ(alea_prepare_root_points(sys, NULL), -1);
    ASSERT_EQ(alea_find_cell(sys, 4, 0, 0), 521);
    ASSERT_EQ(alea_adaptive_child_progress(sys, u, NULL), 0);

    /* Fail the first optional block. Lookup must complete and a later reader
     * must be able to acquire the builder again. End the budget before asserting. */
    alea_alloc_budget_t budget;
    alea_alloc_budget_begin(&budget, 0);
    int owner = 521;
    for (int i = 0; !budget.failed && owner == 521 && i < 100000; ++i)
        owner = alea_find_cell(sys, 4, 0, 0);
    alea_alloc_budget_end(&budget);
    ASSERT(budget.failed);
    ASSERT_EQ(owner, 521);
    ASSERT_EQ(alea_adaptive_child_progress(sys, u, NULL), 0);
    size_t covered = 0;
    for (int i = 0; !covered && i < 100000; ++i) {
        ASSERT_EQ(alea_find_cell(sys, 4, 0, 0), 521);
        covered = alea_adaptive_child_progress(sys, u, NULL);
    }
    ASSERT(covered > 0 && covered < 521);
    ASSERT_EQ(alea_prepare_point_acceleration(sys), 0);
    ASSERT_EQ(alea_parallel_for(16000, 16, 8, ALEA_PARALLEL_DYNAMIC,
                                query_adaptive_child, sys, NULL), ALEA_PARALLEL_OK);
    /* The time needed to qualify the short final block depends on the host. */
    for (int i = 0; i < 100000 && alea_adaptive_child_progress(sys, u, NULL) < 521; ++i)
        ASSERT_EQ(alea_find_cell(sys, 7, 0, 0), 0);
    ASSERT_EQ(alea_adaptive_child_progress(sys, u, NULL), 521);
    ASSERT_EQ(query_adaptive_child(sys, 0, 0, 5), 0);
    alea_selected_query_stats_t before, after;
    ASSERT_EQ(alea_selected_query_stats(sys, &before), 1);
    alea_alloc_budget_begin(&budget, 0);
    int failed_build = alea_promote_point_child(sys, 7);
    alea_alloc_budget_end(&budget);
    ASSERT_EQ(failed_build, -1);
    ASSERT_EQ(alea_selected_query_stats(sys, &after), 1);
    ASSERT_EQ(after.bytes, before.bytes);
    ASSERT_EQ(query_adaptive_child(sys, 0, 0, 5), 0);
    alea_interrupt();
    int interrupted_build = alea_promote_point_child(sys, 7);
    alea_clear_interrupt();
    ASSERT_EQ(interrupted_build, -1);
    ASSERT_EQ(alea_parallel_for(8, 1, 8, ALEA_PARALLEL_STATIC_BLOCK,
                                promote_or_read, sys, NULL), ALEA_PARALLEL_OK);
    ASSERT_EQ(alea_selected_query_stats(sys, &after), 1);
    ASSERT(after.bytes > before.bytes);
    ASSERT_EQ(alea_promote_point_child(sys, 7), 0);
    ASSERT_EQ(alea_parallel_for(512, 16, 8, ALEA_PARALLEL_DYNAMIC,
                                query_adaptive_child, sys, NULL), ALEA_PARALLEL_OK);
    ASSERT_EQ(alea_selected_query_stats(sys, &before), 1);
    ASSERT_EQ(after.bytes, before.bytes);
    /* Repeat with an empty prefix so readers can also publish bounds while
     * the separate child region is being built. */
    alea_clear_selected_queries(sys);
    ASSERT_EQ(alea_prepare_point_acceleration(sys), 0);
    ASSERT_EQ(alea_parallel_for(8, 1, 8, ALEA_PARALLEL_STATIC_BLOCK,
                                promote_or_read, sys, NULL), ALEA_PARALLEL_OK);
    ASSERT(sys->hier_spatial_index == NULL);
    ASSERT_EQ(alea_cell_set_region(sys, 2, sphere(sys, 7)), 0);
    ASSERT_EQ(alea_selected_query_stats(sys, NULL), 0);
    ASSERT_EQ(alea_prepare_point_acceleration(sys), 0);
    ASSERT_EQ(alea_adaptive_child_progress(sys, u, NULL), 0);
    ASSERT_EQ(alea_find_cell(sys, 7, 0, 0), 2);
    alea_destroy(sys);
}

TEST(root_point_exterior_preserves_overlaps_and_gaps) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);
    int mat = alea_add_material(sys, 12);
    ASSERT_EQ(alea_add_cell(sys, 100, sphere(sys, 0), mat, -1, 0), 0);
    int plane = alea_plane_surface(sys, 0, 1, 0, 0, -3);
    alea_node_id_t outside = alea_halfspace(sys, plane, 1);
    ASSERT_EQ(alea_add_cell(sys, 90, outside, mat, -2, 0), 1);
    ASSERT_EQ(alea_add_cell(sys, 10, outside, ALEA_MATERIAL_VOID, 0, 0), 2);
    ASSERT_EQ(alea_prepare_root_points(sys, NULL), 0);
    ASSERT_EQ(alea_selected_query_stats(sys, NULL), 1);
    ASSERT_EQ(alea_find_cell(sys, 0, 0, 0), 0);
    ASSERT_EQ(alea_find_cell(sys, 40, 0, 0), 1);
    ASSERT_EQ(alea_find_cell(sys, -40, 0, 0), -1);
    uint64_t visits;
    size_t world = (size_t)(alea_get_universe(sys, 0) - sys->universes.data);
    ASSERT_EQ(alea_root_point_scan_tests(sys, world, &visits, NULL), 3);
    ASSERT_EQ(visits, 2);
    ASSERT(sys->hier_spatial_index == NULL);
    alea_destroy(sys);
}

typedef struct {
    const interval_index_t* index;
    size_t checks, hits, misses;
} ray_candidate_counts_t;

static int lookup_ray_root(void* context, double x, double y, double z, int surface) {
    ray_candidate_counts_t* counts = context;
    int cell = alea_regional_ray_root_candidate(counts->index, x, y, z, surface, &counts->checks);
    if (cell >= 0) counts->hits++;
    else counts->misses++;
    return cell;
}

TEST(regional_ray_root_preserves_crossed_surface_preference_and_fallback) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);
    int first = alea_sphere_surface(sys, 101, 0, 0, 0, 2);
    int second = alea_sphere_surface(sys, 202, 1, 0, 0, 2);
    int mat = alea_add_material(sys, 15);
    ASSERT_EQ(alea_add_cell(sys, 80, alea_halfspace(sys, first, -1), mat, -1, 0), 0);
    ASSERT_EQ(alea_add_cell(sys, 30, alea_halfspace(sys, second, -1), mat, -2, 0), 1);
    ASSERT_EQ(alea_prepare_query_acceleration(sys), 0);
    const alea_bbox_t domain = {-1, 2, -1, 1, -1, 1};
    ASSERT_EQ(alea_prepare_selected_queries(sys, &domain), 0);
    ray_candidate_counts_t counts = {.index = alea_system_selected_index(sys)};
    ASSERT_EQ(lookup_ray_root(&counts, 0, 0, 0, -1), 0);
    ASSERT_EQ(lookup_ray_root(&counts, 0, 0, 0, 202), 1);
    ASSERT_EQ(lookup_ray_root(&counts, 0, 0, 0, 999), -1);
    ASSERT_EQ(lookup_ray_root(&counts, -1.5, 0, 0, -1), -1);

    const alea_ray_root_lookup_t lookup = {lookup_ray_root, &counts};
    alea_raycast_result_t* ref = alea_raycast_result_create();
    alea_raycast_result_t* trial = alea_raycast_result_create();
    ASSERT_NOT_NULL(ref);
    ASSERT_NOT_NULL(trial);
    counts.hits = counts.misses = 0;
    const double origins[] = {0, -4, 4};
    for (size_t i = 0; i < 3; ++i) {
        double dx = i == 2 ? -1 : 1;
        ASSERT_EQ(alea_raycast_hier_fast_segments(sys, origins[i], 0, 0, dx, 0, 0, 9, ref), 0);
        ASSERT_EQ(alea_raycast_hier_segments_with_root_lookup(
            sys, origins[i], 0, 0, dx, 0, 0, 9, trial, &lookup), 0);
        ASSERT_EQ(alea_raycast_segment_count(ref), alea_raycast_segment_count(trial));
        for (size_t j = 0; j < alea_raycast_segment_count(ref); ++j) {
            const alea_ray_segment_t* a = &ref->segments.data[j];
            const alea_ray_segment_t* b = &trial->segments.data[j];
            ASSERT(a->t_enter == b->t_enter && a->t_exit == b->t_exit);
            ASSERT(a->cell_id == b->cell_id && a->material_id == b->material_id);
            ASSERT(a->density == b->density && a->resolution_flags == b->resolution_flags);
            ASSERT(a->enter_surface_id == b->enter_surface_id &&
                   a->exit_surface_id == b->exit_surface_id);
        }
    }
    ASSERT(counts.hits && counts.misses);
    /* Force a linear root search. Its winner in this overlap is the first
     * cell even with surface 202 as a hint. Let the hierarchy resolve it. */
    alea_system_t* linear = alea_clone(sys);
    ASSERT_NOT_NULL(linear);
    alea_setenv("ALEA_HIER_BLAS_THRESHOLD", "99", 1);
    int prepared = alea_prepare_query_acceleration(linear);
    alea_unsetenv("ALEA_HIER_BLAS_THRESHOLD");
    ASSERT_EQ(prepared, 0);
    ASSERT_EQ(alea_prepare_selected_queries(linear, &domain), 0);
    counts.index = alea_system_selected_index(linear);
    ASSERT_EQ(lookup_ray_root(&counts, 0, 0, 0, 202), -1);
    ASSERT_EQ(alea_hier_spatial_find_ordered_cell_in_universe(linear, 0, 0, 0, 0, 202), 0);
    alea_destroy(linear);
    alea_raycast_result_destroy(trial);
    alea_raycast_result_destroy(ref);
    alea_destroy(sys);
}

TEST(shared_index_preserves_ownership_and_releases_on_edit) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);
    int mat = alea_add_material(sys, 15);
    ASSERT_EQ(alea_add_cell(sys, 80, sphere(sys, 0), mat, -1, 0), 0);
    ASSERT_EQ(alea_add_cell(sys, 30, sphere(sys, 1), mat, -2, 0), 1);
    ASSERT_EQ(alea_add_cell(sys, 90, sphere(sys, 20), ALEA_MATERIAL_VOID, 0, 0), 2);
    ASSERT_EQ(alea_prepare_query_acceleration(sys), 0);
    const double xs[] = {0, .5, nextafter(1, 0), 1, nextafter(1, 2), 2.5, 20, 25};
    int ids[8], materials[8], indices[8], status[8];
    for (size_t i = 0; i < 8; ++i) {
        ids[i] = materials[i] = -91;
        indices[i] = alea_find_cell(sys, xs[i], 0, 0);
        status[i] = alea_find_cell_at(sys, xs[i], 0, 0, &ids[i], &materials[i]);
    }
    const alea_bbox_t box = {-3, 3, -2, 2, -2, 2};
    ASSERT_EQ(alea_prepare_selected_queries(sys, &box), 0);
    alea_selected_query_stats_t stats;
    ASSERT_EQ(alea_selected_query_stats(sys, &stats), 1);
    ASSERT(stats.bytes > 0);
    ASSERT(stats.domain.min_x == box.min_x);
    for (size_t i = 0; i < 8; ++i) {
        int id = -91, material = -91;
        ASSERT_EQ(alea_find_cell_at(sys, xs[i], 0, 0, &id, &material), status[i]);
        ASSERT_EQ(alea_find_cell(sys, xs[i], 0, 0), indices[i]);
        ASSERT_EQ(id, ids[i]);
        ASSERT_EQ(material, materials[i]);
        ASSERT_EQ(alea_selected_query_stats(sys, &stats), 1);
        ASSERT_EQ(stats.fallback_visits, i < 6 ? 0 : (i - 5) * 2);
    }
    ASSERT_EQ(alea_find_cell_at(sys, 0, 0, 0, NULL, NULL), 0);
    ASSERT_EQ(alea_prepare_selected_queries(sys, &box), 0);
    ASSERT_EQ(alea_prepare_selected_queries(sys, NULL), 0);
    const alea_bbox_t other = {-1, 1, -1, 1, -1, 1};
    ASSERT_EQ(alea_prepare_selected_queries(sys, &other), -1);
    ASSERT_EQ(alea_error_code(), ALEA_ERR_INVALID_STATE);

    alea_system_t* clone = alea_clone(sys);
    ASSERT_NOT_NULL(clone);
    ASSERT_EQ(alea_selected_query_stats(clone, NULL), 0);
    ASSERT_EQ(alea_prepare_selected_queries(clone, NULL), 0);
    ASSERT_EQ(alea_find_cell(clone, 20, 0, 0), 2);
    alea_destroy(clone);

    char path[4096];
    FILE* f = alea_tmpfile(path);
    ASSERT_NOT_NULL(f);
    fclose(f);
    ASSERT_EQ(alea_snapshot_write(sys, path), 0);
    clone = alea_snapshot_read(path, 0);
    remove(path);
    ASSERT_NOT_NULL(clone);
    ASSERT_EQ(alea_selected_query_stats(clone, NULL), 0);
    ASSERT_EQ(alea_prepare_selected_queries(clone, NULL), 0);
    ASSERT_EQ(alea_find_cell(clone, .5, 0, 0), 0);
    alea_destroy(clone);

    ASSERT_EQ(alea_cell_set_density(sys, 0, -3), 0);
    ASSERT_EQ(alea_selected_query_stats(sys, NULL), 0);
    ASSERT_EQ(alea_prepare_selected_queries(sys, NULL), 0);
    ASSERT_EQ(alea_renumber_cells(sys, 120), 0);
    ASSERT_EQ(alea_selected_query_stats(sys, NULL), 0);
    ASSERT_EQ(alea_prepare_selected_queries(sys, NULL), 0);
    int cell = -1;
    ASSERT_EQ(alea_find_cell_at(sys, 0, 0, 0, &cell, NULL), 0);
    ASSERT_EQ(cell, 120);
    alea_clear_selected_queries(sys);
    ASSERT_EQ(alea_selected_query_stats(sys, &stats), 0);
    ASSERT_EQ(stats.bytes, 0);
    ASSERT_EQ(alea_find_cell(sys, 20, 0, 0), 2);
    alea_destroy(sys);
}

static int prepare_and_lookup(void* context, size_t worker, size_t begin, size_t end) {
    (void)worker;
    alea_system_t* sys = context;
    if (alea_prepare_selected_queries(sys, NULL)) return -1;
    for (size_t i = begin; i < end; ++i) {
        int cell = -1;
        if (alea_find_cell_at(sys, 0, 0, 0, &cell, NULL) || cell != 29 ||
            alea_find_cell(sys, 10, 0, 0) != -1) return -1;
    }
    return 0;
}

/* Start just at the promotion threshold to exercise publication without
 * relying on machine speed or an artificial delay inside a lookup. */
static void reach_point_promotion_threshold(alea_system_t* sys) {
    atomic_store(&sys->point_lookup_count, 64);
    atomic_store(&sys->point_sample_ns, 7812500);
}

static int automatic_point_reader(void* context, size_t worker, size_t begin, size_t end) {
    (void)worker;
    alea_system_t* sys = context;
    for (size_t i = begin; i < end; ++i) {
        int id = -1;
        if (alea_find_cell_at(sys, 0, 0, 0, &id, NULL) || id != 29 ||
            alea_find_cell(sys, 0, 0, 0) != 0 || alea_material_at(sys, 0, 0, 0) != 7)
            return -1;
    }
    return 0;
}

TEST(automatic_points_publish_reuse_and_reset_on_edit) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);
    int material = alea_add_material(sys, 7);
    ASSERT_EQ(alea_add_cell(sys, 29, sphere(sys, 0), material, -1, 0), 0);
    ASSERT_EQ(alea_prepare_point_lookup(sys), 0);
    ASSERT_EQ(alea_material_at(sys, 0, 0, 0), 7);
    ASSERT_EQ(alea_selected_query_stats(sys, NULL), 0);
    reach_point_promotion_threshold(sys);
    ASSERT_EQ(alea_parallel_for(1024, 1, 8, ALEA_PARALLEL_STATIC_BLOCK,
                                automatic_point_reader, sys, NULL), ALEA_PARALLEL_OK);
    ASSERT_NOT_NULL(sys->automatic_point_index);
    ASSERT_EQ(alea_prepare_point_acceleration(sys), 0);
    ASSERT(sys->selected_index == sys->automatic_point_index);
    const alea_bbox_t box = {-2, 2, -2, 2, -2, 2};
    ASSERT_EQ(alea_prepare_selected_queries(sys, &box), -1);
    alea_clear_selected_queries(sys);
    reach_point_promotion_threshold(sys);
    ASSERT_EQ(alea_find_cell(sys, 0, 0, 0), 0);
    ASSERT_NOT_NULL(sys->automatic_point_index);
    ASSERT_EQ(alea_prepare_selected_queries(sys, &box), 0);
    ASSERT(alea_system_selected_index(sys) == sys->selected_index);
    ASSERT(alea_system_selected_index(sys) != sys->automatic_point_index);
    ASSERT_EQ(alea_renumber_cells(sys, 120), 0);
    ASSERT_EQ(alea_selected_query_stats(sys, NULL), 0);
    ASSERT(!sys->automatic_point_index && !sys->selected_index);
    ASSERT_EQ(atomic_load(&sys->point_lookup_count), 0);
    ASSERT_EQ(atomic_load(&sys->point_sample_ns), 0);
    ASSERT_EQ(atomic_load(&sys->point_promotion_attempted), 0);
    int id;
    ASSERT_EQ(alea_prepare_point_lookup(sys), 0);
    ASSERT_EQ(alea_find_cell_at(sys, 0, 0, 0, &id, NULL), 0);
    ASSERT_EQ(id, 120);
    reach_point_promotion_threshold(sys);
    alea_alloc_fail_after(0);
    int found = alea_find_cell(sys, 0, 0, 0);
    alea_alloc_fail_after(-1);
    ASSERT_EQ(found, 0);
    ASSERT_EQ(alea_selected_query_stats(sys, NULL), 0);
    ASSERT_EQ(atomic_load(&sys->point_promotion_attempted), 1);
    ASSERT_EQ(alea_prepare_point_acceleration(sys), 0);
    ASSERT_EQ(alea_selected_query_stats(sys, NULL), 1);
    alea_destroy(sys);
}

TEST(shared_index_retry_and_concurrent_readers) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);
    ASSERT_EQ(alea_add_cell(sys, 29, sphere(sys, 0), ALEA_MATERIAL_VOID, 0, 0), 0);
    ASSERT_EQ(alea_prepare_query_acceleration(sys), 0);
    alea_interrupt();
    int rc = alea_prepare_selected_queries(sys, NULL);
    alea_clear_interrupt();
    ASSERT_EQ(rc, -1);
    ASSERT_EQ(alea_error_code(), ALEA_ERR_INTERRUPTED);
    ASSERT_EQ(alea_selected_query_stats(sys, NULL), 0);
    alea_alloc_fail_after(0);
    rc = alea_prepare_selected_queries(sys, NULL);
    alea_alloc_fail_after(-1);
    ASSERT_EQ(rc, -1);
    ASSERT_EQ(alea_error_code(), ALEA_ERR_OUT_OF_MEMORY);
    alea_alloc_fail_after(1); /* Fail growing the spatial tree after allocating the index. */
    rc = alea_prepare_selected_queries(sys, NULL);
    alea_alloc_fail_after(-1);
    ASSERT_EQ(rc, -1);
    ASSERT_EQ(alea_selected_query_stats(sys, NULL), 0);
    ASSERT_EQ(alea_find_cell(sys, 0, 0, 0), 0);
    ASSERT_EQ(alea_parallel_for(256, 1, 8, ALEA_PARALLEL_STATIC_BLOCK,
                                prepare_and_lookup, sys, NULL), ALEA_PARALLEL_OK);
    alea_selected_query_stats_t stats;
    ASSERT_EQ(alea_selected_query_stats(sys, &stats), 1);
    ASSERT_EQ(stats.fallback_visits, 256);
    alea_destroy(sys);
}

typedef struct {
    alea_system_t* sys;
    alea_prepared_queries_t* render_index;
    render_config_t config;
    render_camera_t camera;
    render_framebuffer_t* frame;
} shared_reader_scene_t;

static int mixed_readers(void* context, size_t worker, size_t begin, size_t end) {
    (void)worker;
    shared_reader_scene_t* scene = context;
    for (size_t i = begin; i < end; ++i) {
        if (i == 0) {
            if (alea_render_scene_prepared(scene->sys, scene->render_index, &scene->config,
                                           &scene->camera, scene->frame)) return -1;
        } else if (i % 2) {
            double origin[] = {-2, 0, 0}, direction[] = {2, 0, 0};
            alea_raycast_batch_result_t* batch = alea_raycast_batch_result_create();
            if (!batch) return -1;
            alea_raycast_batch_options_t options = {
                .struct_size = sizeof(options), .fields = ALEA_RAY_BATCH_MATERIAL
            };
            int rc = alea_raycast_hier_batch(scene->sys, origin, direction, 1, 4, &options, batch);
            if (!rc && (alea_raycast_batch_segment_count(batch) != 3 ||
                        alea_raycast_batch_cell_ids(batch)[1] != 29)) rc = -1;
            alea_raycast_batch_result_destroy(batch);
            if (rc) return -1;
        } else if (alea_find_cell(scene->sys, 0, 0, 0) != 0) {
            return -1;
        }
    }
    return 0;
}

TEST(shared_index_coexists_with_batch_and_render_readers) {
    shared_reader_scene_t scene = {.sys = alea_create()};
    ASSERT_NOT_NULL(scene.sys);
    int material = alea_add_material(scene.sys, 7);
    ASSERT_EQ(alea_add_cell(scene.sys, 29, sphere(scene.sys, 0), material, -1, 0), 0);
    ASSERT_EQ(alea_prepare_query_acceleration(scene.sys), 0);
    ASSERT_EQ(alea_prepare_selected_queries(scene.sys, NULL), 0);
    const alea_bbox_t box = {-2, 2, -2, 2, -2, 2};
    scene.render_index = alea_prepared_queries_create(scene.sys, &box);
    ASSERT_NOT_NULL(scene.render_index);
    ASSERT_EQ(alea_set_ray_batch_accelerator(scene.sys, scene.render_index), 0);
    render_config_init(&scene.config);
    scene.config.width = scene.config.height = 9;
    scene.config.eye_set = scene.config.target_set = 1;
    scene.config.eye[1] = -3;
    ASSERT_EQ(render_camera_setup(&scene.camera, &scene.config, scene.sys), 0);
    scene.frame = render_framebuffer_create(9, 9, 1);
    ASSERT_NOT_NULL(scene.frame);
    ASSERT_EQ(alea_parallel_for(64, 1, 8, ALEA_PARALLEL_STATIC_BLOCK,
                                mixed_readers, &scene, NULL), ALEA_PARALLEL_OK);
    ASSERT_EQ(scene.frame->cell_id[40], 29);
    render_framebuffer_free(scene.frame);
    render_config_free(&scene.config);
    alea_prepared_queries_destroy(scene.render_index);
    alea_destroy(scene.sys);
}

TEST(public_ray_batch_dispatch_limits_and_invalidation) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);
    int material = alea_add_material(sys, 17);
    ASSERT_EQ(alea_add_cell(sys, 29, sphere(sys, 0), material, -1, 0), 0);
    alea_bbox_t domain = {-4, 4, -4, 4, -4, 4};
    alea_prepared_queries_t* handle = alea_prepared_queries_create(sys, &domain);
    ASSERT_NOT_NULL(handle);
    ASSERT_EQ(alea_set_ray_batch_accelerator(sys, handle), 0);
    ASSERT_EQ(alea_prepare_ray_acceleration(sys, &domain), 0);
    ASSERT(sys->ray_batch_index == handle);
    alea_system_t* other = alea_create();
    ASSERT_EQ(alea_set_ray_batch_accelerator(other, handle), -1);
    alea_destroy(other);
    alea_prepared_queries_destroy(handle); /* The model owns the remaining reference. */
    double o[] = {-3,0,0, 3,0,0}, d[] = {2,0,0, -2,0,0};
    double mins[] = {0,0}, limits[] = {6,5};
    alea_ray_batch_query_t query = {.kind = ALEA_RAY_QUERY_SEGMENTS, .t_mins = mins, .t_maxs = limits};
    alea_raycast_batch_options_t options = {.struct_size = sizeof(options), .fields = ALEA_RAY_BATCH_DENSITY};
    alea_raycast_batch_result_t* result = alea_raycast_batch_result_create();
    ASSERT_NOT_NULL(result);
    alea_raycast_batch_work_stats_t stats;
    const uint32_t fields[] = {ALEA_RAY_BATCH_DENSITY, ALEA_RAY_BATCH_DENSITY,
        ALEA_RAY_BATCH_DENSITY, ALEA_RAY_BATCH_DENSITY, ALEA_RAY_BATCH_SURFACES,
        ALEA_RAY_BATCH_PROJECTED_OWNER, ALEA_RAY_BATCH_FULL_PATHS};
    for (size_t pass = 0; pass < sizeof(fields) / sizeof(fields[0]); ++pass) {
        mins[1] = pass == 1 ? 1 : 0;
        limits[1] = pass == 2 ? 0 : pass == 3 ? DBL_MAX : 5;
        options.fields = fields[pass];
        ASSERT_EQ(alea_raycast_hier_batch_query_nocache(sys, o, d, 2, &query, &options, result), 0);
        ASSERT_EQ(alea_raycast_batch_result_get_work_stats_internal(result, &stats), 0);
        ASSERT_EQ(stats.regional_ray_count, pass == 0 ? 2 : 0);
        ASSERT_EQ(alea_raycast_batch_cell_ids(result)[1], 29);
        ASSERT_EQ(alea_raycast_batch_t_enter(result)[1], 2);
        ASSERT_EQ(alea_raycast_batch_t_exit(result)[1], 4);
        if (pass == 4) ASSERT_NOT_NULL(alea_raycast_batch_enter_surface_ids(result));
    }
    options.fields = ALEA_RAY_BATCH_DENSITY;
    options.max_segments = 1;
    ASSERT_EQ(alea_raycast_hier_batch(sys, o, d, 2, 6, &options, result), -1);
    ASSERT_EQ(alea_error_code(), ALEA_ERR_OVERFLOW);
    ASSERT_EQ(alea_raycast_batch_fields(result), ALEA_RAY_BATCH_FULL_PATHS); /* prior output survives */
    options.max_segments = 0;
    options.max_output_bytes = 1;
    ASSERT_EQ(alea_raycast_hier_batch(sys, o, d, 2, 6, &options, result), -1);
    ASSERT_EQ(alea_error_code(), ALEA_ERR_OVERFLOW);
    options.max_output_bytes = 0;
    ASSERT_EQ(alea_cell_set_density(sys, 0, -3), 0);
    ASSERT_NULL(sys->ray_batch_index);
    ASSERT_EQ(alea_raycast_hier_batch(sys, o, d, 2, 6, &options, result), 0);
    ASSERT_EQ(alea_raycast_batch_result_get_work_stats_internal(result, &stats), 0);
    ASSERT_EQ(stats.regional_ray_count, 0);
    ASSERT_EQ(alea_raycast_batch_densities(result)[1], 3);
    ASSERT_EQ(alea_prepare_ray_acceleration(sys, NULL), 0);
    ASSERT_NOT_NULL(sys->ray_batch_index);
    ASSERT_EQ(alea_cell_set_region(sys, 0, sphere(sys, 8)), 0);
    ASSERT_NULL(sys->ray_batch_index);
    ASSERT_EQ(alea_raycast_hier_batch(sys, o, d, 2, 6, &options, result), 0);
    ASSERT_EQ(alea_raycast_batch_segment_count(result), 2);
    ASSERT_EQ(alea_prepare_ray_acceleration(sys, &domain), 0);
    alea_system_reset(sys);
    ASSERT_NULL(sys->ray_batch_index);
    alea_raycast_batch_result_destroy(result);
    alea_destroy(sys);
}

TEST(shared_index_without_bounded_world) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);
    int plane = alea_plane_surface(sys, 1, 1, 0, 0, 0);
    ASSERT_EQ(alea_add_cell(sys, 11, alea_halfspace(sys, plane, 1), ALEA_MATERIAL_VOID, 0, 0), 0);
    ASSERT_EQ(alea_prepare_selected_queries(sys, NULL), 0);
    ASSERT_EQ(alea_selected_query_stats(sys, NULL), 0);
    ASSERT_EQ(alea_find_cell(sys, 2, 0, 0), 0);
    const alea_bbox_t invalid = {-INFINITY, INFINITY, -2, 2, -2, 2};
    ASSERT_EQ(alea_prepare_selected_queries(sys, &invalid), -1);
    ASSERT_EQ(alea_error_code(), ALEA_ERR_INVALID_ARG);
    const alea_bbox_t box = {-2, 2, -2, 2, -2, 2};
    ASSERT_EQ(alea_prepare_selected_queries(sys, &box), 0);
    ASSERT_EQ(alea_selected_query_stats(sys, NULL), 1);
    ASSERT_EQ(alea_find_cell(sys, 20, 0, 0), 0);
    alea_cell_hit_t expected;
    int rc = alea_find_deepest_cell_hit_at_point(sys, 1e39, 0, 0, &expected);
    ASSERT_EQ(alea_find_cell(sys, 1e39, 0, 0), rc ? -1 : expected.cell_index);
    alea_destroy(sys);
}

TEST(shared_index_follows_transformed_and_repeating_lattice_ownership) {
    const char* models[] = {"tests/data/mcnp_nested_lattice_transformed.mcnp",
                            "tests/data/mcnp_lattice_repeating.mcnp"};
    for (size_t m = 0; m < 2; ++m) {
        mcnp_model_t* model = mcnp_load(models[m]);
        ASSERT_NOT_NULL(model);
        alea_system_t* sys = model->sys;
        ASSERT_EQ(alea_prepare_selected_queries(sys, NULL), 0);
        ASSERT_EQ(alea_selected_query_stats(sys, NULL), 1);
        for (int i = -8; i <= 40; ++i) {
            double x = i * .25;
            double xs[] = {nextafter(x, -INFINITY), x, nextafter(x, INFINITY)};
            for (size_t j = 0; j < 3; ++j) {
                for (int k = -5; k <= 5; ++k) {
                    double y = k * .25;
                    alea_cell_hit_t expected;
                    int rc = alea_find_deepest_cell_hit_at_point(sys, xs[j], y, 0, &expected);
                    int id = -1, material = -1;
                    ASSERT_EQ(alea_find_cell_at(sys, xs[j], y, 0, &id, &material), rc);
                    ASSERT_EQ(alea_find_cell(sys, xs[j], y, 0), rc ? -1 : expected.cell_index);
                    ASSERT_EQ(id, rc ? -1 : expected.cell_id);
                    ASSERT_EQ(material, rc ? -1 : expected.material_id);
                }
            }
        }
        mcnp_model_destroy(model);
    }
}

TEST_MAIN()
