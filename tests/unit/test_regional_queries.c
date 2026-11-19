// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
// SPDX-License-Identifier: MPL-2.0

#include "alea_test.h"
#include "alea.h"
#include "alea_raycast.h"
#include "alea_mcnp.h"
#include "alea_snapshot.h"
#include "util/compat.h"
#include <float.h>
#include "core/alea_eval.h"
#include "core/alea_regional.h"
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

TEST_MAIN()
