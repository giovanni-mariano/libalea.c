// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

#include "alea_test.h"
#include "alea.h"
#include "alea_geo_validator.h"

static alea_system_t* build_split_box_system(void) {
    alea_system_t* sys = alea_create();
    if (!sys) return NULL;

    int plane_si = alea_plane_surface(sys, 10, 1, 0, 0, 0);
    int box_si = alea_box_surface(sys, 20, -3, 3, -3, 3, -3, 3);
    if (plane_si < 0 || box_si < 0) goto fail;

    alea_node_id_t left_half = alea_halfspace(sys, plane_si, -1);
    alea_node_id_t right_half = alea_halfspace(sys, plane_si, 1);
    alea_node_id_t box = alea_halfspace(sys, box_si, -1);
    alea_node_id_t left = alea_intersection(sys, box, left_half);
    alea_node_id_t right = alea_intersection(sys, box, right_half);
    if (left_half == ALEA_NODE_ID_INVALID ||
        right_half == ALEA_NODE_ID_INVALID ||
        box == ALEA_NODE_ID_INVALID ||
        left == ALEA_NODE_ID_INVALID ||
        right == ALEA_NODE_ID_INVALID) {
        goto fail;
    }

    int m1 = alea_add_material(sys, 1);
    int m2 = alea_add_material(sys, 2);
    if (alea_add_cell(sys, 1, left, m1, 1.0, 0) < 0 ||
        alea_add_cell(sys, 2, right, m2, 1.0, 0) < 0) {
        goto fail;
    }
    return sys;

fail:
    alea_destroy(sys);
    return NULL;
}

static int count_error_type(const alea_geom_validator_result_t* result,
                            alea_geom_error_type_t type) {
    int count = 0;
    for (size_t i = 0; i < result->error_count; i++) {
        if (result->errors[i].type == type) count++;
    }
    return count;
}

static int has_error_flag(const alea_geom_validator_result_t* result,
                          uint32_t flag) {
    for (size_t i = 0; i < result->error_count; i++) {
        if (result->errors[i].flags & flag) return 1;
    }
    return 0;
}

TEST(geo_validator_clean_adjacent_flat) {
    alea_system_t* sys = build_split_box_system();
    ASSERT_NOT_NULL(sys);

    alea_geom_validator_options_t opts;
    alea_geom_validator_options_init(&opts);
    opts.flags |= ALEA_GEOM_VALIDATE_ALLOW_EXTERIOR_VOID;
    opts.ray_count = 24;
    opts.seed = 7;

    alea_geom_validator_result_t result;
    alea_geom_validator_result_init(&result);
    ASSERT_EQ(alea_validate_geometry(sys, &opts, &result), 0);
    ASSERT_EQ(result.error_count, 0);
    ASSERT(result.crossings_checked > 0);

    alea_geom_validator_result_free(&result);
    alea_destroy(sys);
}

TEST(geo_validator_detects_nested_overlap_flat) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);

    int large_si = alea_sphere_surface(sys, 10, 0, 0, 0, 3.0);
    int small_si = alea_sphere_surface(sys, 20, 0, 0, 0, 1.0);
    ASSERT(large_si >= 0);
    ASSERT(small_si >= 0);

    alea_node_id_t large = alea_halfspace(sys, large_si, -1);
    alea_node_id_t small = alea_halfspace(sys, small_si, -1);
    ASSERT_NE(large, ALEA_NODE_ID_INVALID);
    ASSERT_NE(small, ALEA_NODE_ID_INVALID);

    int m1 = alea_add_material(sys, 1);
    int m2 = alea_add_material(sys, 2);
    ASSERT(alea_add_cell(sys, 1, large, m1, 1.0, 0) >= 0);
    ASSERT(alea_add_cell(sys, 2, small, m2, 1.0, 0) >= 0);

    alea_geom_validator_options_t opts;
    alea_geom_validator_options_init(&opts);
    opts.flags |= ALEA_GEOM_VALIDATE_ALLOW_EXTERIOR_VOID;
    opts.ray_count = 8;
    opts.seed = 11;

    alea_geom_validator_result_t result;
    alea_geom_validator_result_init(&result);
    ASSERT_EQ(alea_validate_geometry(sys, &opts, &result), 0);
    ASSERT(count_error_type(&result, ALEA_GEOM_ERR_OVERLAP_AFTER_CROSSING) > 0);

    alea_geom_validator_result_free(&result);
    alea_destroy(sys);
}

TEST(geo_validator_reports_undefined_after_crossing) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);

    int sphere_si = alea_sphere_surface(sys, 10, 0, 0, 0, 1.0);
    ASSERT(sphere_si >= 0);
    alea_node_id_t sphere = alea_halfspace(sys, sphere_si, -1);
    ASSERT_NE(sphere, ALEA_NODE_ID_INVALID);

    int m1 = alea_add_material(sys, 1);
    ASSERT(alea_add_cell(sys, 1, sphere, m1, 1.0, 0) >= 0);

    alea_geom_validator_options_t opts;
    alea_geom_validator_options_init(&opts);
    opts.sample_offset = 0.01;

    alea_geom_validator_result_t result;
    alea_geom_validator_result_init(&result);
    ASSERT_EQ(alea_validate_geometry_ray(sys, &opts,
                                         0, 0, 0, 1, 0, 0, 3,
                                         &result), 0);
    ASSERT(count_error_type(&result, ALEA_GEOM_ERR_UNDEFINED_AFTER_CROSSING) > 0);

    alea_geom_validator_result_free(&result);
    alea_destroy(sys);
}

TEST(geo_validator_allows_exterior_void) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);

    int sphere_si = alea_sphere_surface(sys, 10, 0, 0, 0, 1.0);
    ASSERT(sphere_si >= 0);
    alea_node_id_t sphere = alea_halfspace(sys, sphere_si, -1);
    ASSERT_NE(sphere, ALEA_NODE_ID_INVALID);

    int m1 = alea_add_material(sys, 1);
    ASSERT(alea_add_cell(sys, 1, sphere, m1, 1.0, 0) >= 0);

    alea_geom_validator_options_t opts;
    alea_geom_validator_options_init(&opts);
    opts.flags |= ALEA_GEOM_VALIDATE_ALLOW_EXTERIOR_VOID;
    opts.sample_offset = 0.01;

    alea_geom_validator_result_t result;
    alea_geom_validator_result_init(&result);
    ASSERT_EQ(alea_validate_geometry_ray(sys, &opts,
                                         0, 0, 0, 1, 0, 0, 3,
                                         &result), 0);
    ASSERT_EQ(result.error_count, 0);

    alea_geom_validator_result_free(&result);
    alea_destroy(sys);
}

TEST(geo_validator_reports_missing_neighbor) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);

    int inner_si = alea_sphere_surface(sys, 10, 0, 0, 0, 1.0);
    int box_si = alea_box_surface(sys, 30, 1.005, 3, -1, 1, -1, 1);
    ASSERT(inner_si >= 0);
    ASSERT(box_si >= 0);

    alea_node_id_t inner = alea_halfspace(sys, inner_si, -1);
    alea_node_id_t box = alea_halfspace(sys, box_si, -1);
    ASSERT_NE(inner, ALEA_NODE_ID_INVALID);
    ASSERT_NE(box, ALEA_NODE_ID_INVALID);

    int m1 = alea_add_material(sys, 1);
    int m2 = alea_add_material(sys, 2);
    ASSERT(alea_add_cell(sys, 1, inner, m1, 1.0, 0) >= 0);
    ASSERT(alea_add_cell(sys, 2, box, m2, 1.0, 0) >= 0);

    alea_geom_validator_options_t opts;
    alea_geom_validator_options_init(&opts);
    opts.flags |= ALEA_GEOM_VALIDATE_ALLOW_EXTERIOR_VOID;
    opts.sample_offset = 0.01;

    alea_geom_validator_result_t result;
    alea_geom_validator_result_init(&result);
    ASSERT_EQ(alea_validate_geometry_ray(sys, &opts,
                                         0, 0, 0, 1, 0, 0, 3,
                                         &result), 0);
    ASSERT(count_error_type(&result, ALEA_GEOM_ERR_MISSING_NEIGHBOR) > 0);
    ASSERT(has_error_flag(&result, ALEA_GEOM_EVENT_MISSING_ADJACENCY));

    alea_geom_validator_result_free(&result);
    alea_destroy(sys);
}

TEST(geo_validator_reports_non_adjacent_transition) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);

    int boundary_si = alea_plane_surface(sys, 10, 1, 0, 0, 0);
    int neg_si = alea_plane_surface(sys, 20, 1, 0, 0, -0.5);
    int far_si = alea_plane_surface(sys, 30, 1, 0, 0, 1.0);
    int box_si = alea_box_surface(sys, 40, -3, 100, -1, 1, -1, 1);
    int right_box_si = alea_box_surface(sys, 50, 0.5, 100, -1, 1, -1, 1);
    ASSERT(boundary_si >= 0);
    ASSERT(neg_si >= 0);
    ASSERT(far_si >= 0);
    ASSERT(box_si >= 0);
    ASSERT(right_box_si >= 0);

    alea_node_id_t left_half = alea_halfspace(sys, boundary_si, -1);
    alea_node_id_t right_half = alea_halfspace(sys, boundary_si, 1);
    alea_node_id_t impossible_left = alea_halfspace(sys, far_si, -1);
    alea_node_id_t box = alea_halfspace(sys, box_si, -1);
    alea_node_id_t right_box = alea_halfspace(sys, right_box_si, -1);
    alea_node_id_t left = alea_intersection(sys, box, left_half);
    alea_node_id_t bogus_neighbor = alea_intersection(sys, right_half,
                                                      impossible_left);
    ASSERT_NE(left, ALEA_NODE_ID_INVALID);
    ASSERT_NE(bogus_neighbor, ALEA_NODE_ID_INVALID);
    ASSERT_NE(right_box, ALEA_NODE_ID_INVALID);

    int m1 = alea_add_material(sys, 1);
    int m2 = alea_add_material(sys, 2);
    int m3 = alea_add_material(sys, 3);
    ASSERT(alea_add_cell(sys, 1, left, m1, 1.0, 0) >= 0);
    ASSERT(alea_add_cell(sys, 2, bogus_neighbor, m2, 1.0, 0) >= 0);
    ASSERT(alea_add_cell(sys, 3, right_box, m3, 1.0, 0) >= 0);

    alea_geom_validator_options_t opts;
    alea_geom_validator_options_init(&opts);
    opts.flags |= ALEA_GEOM_VALIDATE_ALLOW_EXTERIOR_VOID;
    opts.sample_offset = 0.6;

    alea_geom_validator_result_t result;
    alea_geom_validator_result_init(&result);
    ASSERT_EQ(alea_validate_geometry_ray(sys, &opts,
                                         -2, 0, 0, 1, 0, 0, 5,
                                         &result), 0);
    ASSERT(count_error_type(&result, ALEA_GEOM_ERR_NON_ADJACENT_TRANSITION) > 0);

    alea_geom_validator_result_free(&result);
    alea_destroy(sys);
}

TEST(geo_validator_reports_ambiguous_boundary) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);

    /* Genuinely coincident, geometrically DISTINCT primitives at the same
     * point: plane x=0 and a sphere centered at (1,0,0) r=1 both pass through
     * the origin.  (Two identical plane cards would instead deduplicate to one
     * primitive and validate cleanly, which is the dedup case, not ambiguity.) */
    int plane_si = alea_plane_surface(sys, 10, 1, 0, 0, 0);
    int sphere_si = alea_sphere_surface(sys, 20, 1, 0, 0, 1.0);
    int box_si = alea_box_surface(sys, 30, -3, 3, -1, 1, -1, 1);
    ASSERT(plane_si >= 0);
    ASSERT(sphere_si >= 0);
    ASSERT(box_si >= 0);

    alea_node_id_t left_half = alea_halfspace(sys, plane_si, -1);
    alea_node_id_t right_half = alea_halfspace(sys, plane_si, 1);
    alea_node_id_t sphere_in = alea_halfspace(sys, sphere_si, -1);
    alea_node_id_t sphere_out = alea_halfspace(sys, sphere_si, 1);
    alea_node_id_t box = alea_halfspace(sys, box_si, -1);
    alea_node_id_t left = alea_intersection(sys, box, left_half);
    alea_node_id_t right = alea_intersection(sys,
                                             alea_intersection(sys, box, right_half),
                                             sphere_out);
    ASSERT_NE(left, ALEA_NODE_ID_INVALID);
    ASSERT_NE(right, ALEA_NODE_ID_INVALID);
    ASSERT_NE(sphere_in, ALEA_NODE_ID_INVALID);

    int m1 = alea_add_material(sys, 1);
    int m2 = alea_add_material(sys, 2);
    int m3 = alea_add_material(sys, 3);
    ASSERT(alea_add_cell(sys, 1, left, m1, 1.0, 0) >= 0);
    ASSERT(alea_add_cell(sys, 2, right, m2, 1.0, 0) >= 0);
    ASSERT(alea_add_cell(sys, 3, sphere_in, m3, 1.0, 0) >= 0);

    alea_geom_validator_options_t opts;
    alea_geom_validator_options_init(&opts);
    opts.flags |= ALEA_GEOM_VALIDATE_ALLOW_EXTERIOR_VOID;
    opts.sample_offset = 0.01;

    alea_geom_validator_result_t result;
    alea_geom_validator_result_init(&result);
    ASSERT_EQ(alea_validate_geometry_ray(sys, &opts,
                                         -2, 0, 0, 1, 0, 0, 4,
                                         &result), 0);
    ASSERT(count_error_type(&result, ALEA_GEOM_ERR_AMBIGUOUS_BOUNDARY) > 0);

    alea_geom_validator_result_free(&result);
    alea_destroy(sys);
}

TEST(geo_validator_deduplicated_surface_cards_are_clean) {
    /* Two distinct surface cards (10 and 20) describe the SAME plane x=0 and
     * fold onto one canonical primitive, exactly as MCNP would deduplicate
     * them.  Cell 1 is written against card 10 and cell 2 against card 20.
     * Matching on canonical primitive identity (not mc_surface_id) must see a
     * clean shared boundary and report no errors. */
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);

    int p10 = alea_plane_surface(sys, 10, 1, 0, 0, 0);
    int p20 = alea_plane_surface(sys, 20, 1, 0, 0, 0);
    int box_si = alea_box_surface(sys, 30, -3, 3, -3, 3, -3, 3);
    ASSERT(p10 >= 0);
    ASSERT(p20 >= 0);
    ASSERT(box_si >= 0);

    alea_node_id_t left_half = alea_halfspace(sys, p10, -1);
    alea_node_id_t right_half = alea_halfspace(sys, p20, 1);
    alea_node_id_t box = alea_halfspace(sys, box_si, -1);
    alea_node_id_t left = alea_intersection(sys, box, left_half);
    alea_node_id_t right = alea_intersection(sys, box, right_half);
    ASSERT_NE(left, ALEA_NODE_ID_INVALID);
    ASSERT_NE(right, ALEA_NODE_ID_INVALID);

    int m1 = alea_add_material(sys, 1);
    int m2 = alea_add_material(sys, 2);
    ASSERT(alea_add_cell(sys, 1, left, m1, 1.0, 0) >= 0);
    ASSERT(alea_add_cell(sys, 2, right, m2, 1.0, 0) >= 0);

    alea_geom_validator_options_t opts;
    alea_geom_validator_options_init(&opts);
    opts.flags |= ALEA_GEOM_VALIDATE_ALLOW_EXTERIOR_VOID;
    opts.ray_count = 24;
    opts.seed = 7;

    alea_geom_validator_result_t result;
    alea_geom_validator_result_init(&result);
    ASSERT_EQ(alea_validate_geometry(sys, &opts, &result), 0);
    ASSERT_EQ(result.error_count, 0);
    ASSERT(result.crossings_checked > 0);

    alea_geom_validator_result_free(&result);
    alea_destroy(sys);
}

TEST(geo_validator_truncates_at_max_errors) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);

    int large_si = alea_sphere_surface(sys, 10, 0, 0, 0, 3.0);
    int small_si = alea_sphere_surface(sys, 20, 0, 0, 0, 1.0);
    ASSERT(large_si >= 0);
    ASSERT(small_si >= 0);

    alea_node_id_t large = alea_halfspace(sys, large_si, -1);
    alea_node_id_t small = alea_halfspace(sys, small_si, -1);
    int m1 = alea_add_material(sys, 1);
    int m2 = alea_add_material(sys, 2);
    ASSERT(alea_add_cell(sys, 1, large, m1, 1.0, 0) >= 0);
    ASSERT(alea_add_cell(sys, 2, small, m2, 1.0, 0) >= 0);

    alea_geom_validator_options_t opts;
    alea_geom_validator_options_init(&opts);
    opts.flags |= ALEA_GEOM_VALIDATE_ALLOW_EXTERIOR_VOID;
    opts.ray_count = 32;
    opts.max_errors = 1;
    opts.seed = 11;

    alea_geom_validator_result_t result;
    alea_geom_validator_result_init(&result);
    ASSERT_EQ(alea_validate_geometry(sys, &opts, &result), 0);
    ASSERT_EQ(result.error_count, 1);
    ASSERT_EQ(result.truncated, 1);

    alea_geom_validator_result_free(&result);
    alea_destroy(sys);
}

TEST(geo_validator_clean_adjacent_hier) {
    alea_system_t* sys = build_split_box_system();
    ASSERT_NOT_NULL(sys);
    alea_config_t cfg = alea_get_config(sys);
    cfg.spatial_mode = ALEA_SPATIAL_MODE_HIER;
    alea_set_config(sys, &cfg);

    alea_geom_validator_options_t opts;
    alea_geom_validator_options_init(&opts);
    opts.flags |= ALEA_GEOM_VALIDATE_ALLOW_EXTERIOR_VOID |
                  ALEA_GEOM_VALIDATE_HIERARCHICAL;
    opts.ray_count = 24;
    opts.seed = 7;

    alea_geom_validator_result_t result;
    alea_geom_validator_result_init(&result);
    ASSERT_EQ(alea_validate_geometry(sys, &opts, &result), 0);
    ASSERT_EQ(result.error_count, 0);
    ASSERT(result.crossings_checked > 0);

    alea_geom_validator_result_free(&result);
    alea_destroy(sys);
}

TEST(geo_validator_non_strict_uses_fast_adjacency_path) {
    alea_system_t* sys = build_split_box_system();
    ASSERT_NOT_NULL(sys);

    alea_geom_validator_options_t strict_opts;
    alea_geom_validator_options_init(&strict_opts);
    strict_opts.flags |= ALEA_GEOM_VALIDATE_ALLOW_EXTERIOR_VOID;
    strict_opts.sample_offset = 0.01;

    alea_geom_validator_result_t strict_result;
    alea_geom_validator_result_init(&strict_result);
    ASSERT_EQ(alea_validate_geometry_ray(sys, &strict_opts,
                                         -2, 0, 0, 1, 0, 0, 4,
                                         &strict_result), 0);

    alea_geom_validator_options_t fast_opts = strict_opts;
    fast_opts.flags &= ~ALEA_GEOM_VALIDATE_STRICT_ADJACENCY;

    alea_geom_validator_result_t fast_result;
    alea_geom_validator_result_init(&fast_result);
    ASSERT_EQ(alea_validate_geometry_ray(sys, &fast_opts,
                                         -2, 0, 0, 1, 0, 0, 4,
                                         &fast_result), 0);
    ASSERT_EQ(fast_result.error_count, 0);
    ASSERT(fast_result.exact_queries < strict_result.exact_queries);
    ASSERT(fast_result.adjacency_hits > 0);

    alea_geom_validator_result_free(&strict_result);
    alea_geom_validator_result_free(&fast_result);
    alea_destroy(sys);
}

TEST(geo_validator_rays_flag_can_disable_whole_geometry_scan) {
    alea_system_t* sys = build_split_box_system();
    ASSERT_NOT_NULL(sys);

    alea_geom_validator_options_t opts;
    alea_geom_validator_options_init(&opts);
    opts.flags &= ~ALEA_GEOM_VALIDATE_RAYS;

    alea_geom_validator_result_t result;
    alea_geom_validator_result_init(&result);
    ASSERT_EQ(alea_validate_geometry(sys, &opts, &result), 0);
    ASSERT_EQ(result.crossings_checked, 0);
    ASSERT_EQ(result.exact_queries, 0);

    alea_geom_validator_result_free(&result);
    alea_destroy(sys);
}

TEST(geo_validator_late_ladder_change_does_not_mark_thin_cell_ambiguous) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);

    int p0_si = alea_plane_surface(sys, 10, 1, 0, 0, 0);
    int pthin_si = alea_plane_surface(sys, 20, 1, 0, 0, -0.005);
    int box_si = alea_box_surface(sys, 30, -3, 3, -1, 1, -1, 1);
    ASSERT(p0_si >= 0);
    ASSERT(pthin_si >= 0);
    ASSERT(box_si >= 0);

    alea_node_id_t box = alea_halfspace(sys, box_si, -1);
    alea_node_id_t left = alea_intersection(sys, box,
                                            alea_halfspace(sys, p0_si, -1));
    alea_node_id_t thin_left = alea_halfspace(sys, p0_si, 1);
    alea_node_id_t thin_right = alea_halfspace(sys, pthin_si, -1);
    alea_node_id_t thin = alea_intersection(sys, box,
        alea_intersection(sys, thin_left, thin_right));
    alea_node_id_t right = alea_intersection(sys, box,
                                             alea_halfspace(sys, pthin_si, 1));
    ASSERT_NE(left, ALEA_NODE_ID_INVALID);
    ASSERT_NE(thin, ALEA_NODE_ID_INVALID);
    ASSERT_NE(right, ALEA_NODE_ID_INVALID);

    int m1 = alea_add_material(sys, 1);
    int m2 = alea_add_material(sys, 2);
    int m3 = alea_add_material(sys, 3);
    ASSERT(alea_add_cell(sys, 1, left, m1, 1.0, 0) >= 0);
    ASSERT(alea_add_cell(sys, 2, thin, m2, 1.0, 0) >= 0);
    ASSERT(alea_add_cell(sys, 3, right, m3, 1.0, 0) >= 0);

    alea_geom_validator_options_t opts;
    alea_geom_validator_options_init(&opts);
    opts.flags |= ALEA_GEOM_VALIDATE_ALLOW_EXTERIOR_VOID;
    opts.sample_offset = 0.001;

    alea_geom_validator_result_t result;
    alea_geom_validator_result_init(&result);
    ASSERT_EQ(alea_validate_geometry_ray(sys, &opts,
                                         -1, 0, 0, 1, 0, 0, 2,
                                         &result), 0);
    ASSERT_EQ(count_error_type(&result, ALEA_GEOM_ERR_AMBIGUOUS_BOUNDARY), 0);

    alea_geom_validator_result_free(&result);
    alea_destroy(sys);
}

TEST(geo_validator_reports_overlap_count) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);

    int s1 = alea_sphere_surface(sys, 10, 0, 0, 0, 3.0);
    int s2 = alea_sphere_surface(sys, 20, 0, 0, 0, 2.0);
    int s3 = alea_sphere_surface(sys, 30, 0, 0, 0, 1.0);
    ASSERT(s1 >= 0);
    ASSERT(s2 >= 0);
    ASSERT(s3 >= 0);

    int m1 = alea_add_material(sys, 1);
    int m2 = alea_add_material(sys, 2);
    int m3 = alea_add_material(sys, 3);
    ASSERT(alea_add_cell(sys, 1, alea_halfspace(sys, s1, -1), m1, 1.0, 0) >= 0);
    ASSERT(alea_add_cell(sys, 2, alea_halfspace(sys, s2, -1), m2, 1.0, 0) >= 0);
    ASSERT(alea_add_cell(sys, 3, alea_halfspace(sys, s3, -1), m3, 1.0, 0) >= 0);

    alea_geom_validator_options_t opts;
    alea_geom_validator_options_init(&opts);
    opts.flags |= ALEA_GEOM_VALIDATE_ALLOW_EXTERIOR_VOID;
    opts.sample_offset = 0.01;

    alea_geom_validator_result_t result;
    alea_geom_validator_result_init(&result);
    ASSERT_EQ(alea_validate_geometry_ray(sys, &opts,
                                         -4, 0, 0, 1, 0, 0, 8,
                                         &result), 0);

    int saw_three_way = 0;
    for (size_t i = 0; i < result.error_count; i++) {
        if (result.errors[i].type == ALEA_GEOM_ERR_OVERLAP_AFTER_CROSSING &&
            result.errors[i].found_cell_count >= 3) {
            saw_three_way = 1;
        }
    }
    ASSERT(saw_three_way);

    alea_geom_validator_result_free(&result);
    alea_destroy(sys);
}

TEST_MAIN()
