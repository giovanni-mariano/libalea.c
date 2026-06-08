// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/*
 * test_error_handling.c - Error handling and robustness tests
 */

#include "alea_test.h"
#include "alea.h"
#include "alea_mcnp.h"
#include "core/alea_system.h"
#include "core/alea_eval.h"
#include <string.h>
#include <stdio.h>

/* ========================================================================= */
/* NULL argument handling                                                    */
/* ========================================================================= */

TEST(null_system_create) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);
    alea_destroy(sys);
}

TEST(null_system_destroy) {
    /* Should not crash */
    alea_destroy(NULL);
}

TEST(null_find_cell) {
    int result = alea_find_cell(NULL, 0, 0, 0);
    ASSERT(result < 0);
}

TEST(null_export) {
    int result = mcnp_export(NULL, "null_test_tmp.mcnp");
    ASSERT(result != 0);
}

TEST(null_load_mcnp) {
    mcnp_model_t* model = mcnp_load(NULL);
    ASSERT_NULL(model);
}

TEST(null_load_mcnp_string) {
    mcnp_model_t* model = mcnp_load_string(NULL, 0);
    ASSERT_NULL(model);
}

/* ========================================================================= */
/* Invalid IDs                                                               */
/* ========================================================================= */

TEST(invalid_node_id) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);

    /* Evaluating invalid node should not crash */
    bool inside = alea_point_inside(sys, ALEA_NODE_ID_INVALID, 0, 0, 0);
    ASSERT_FALSE(inside);

    alea_destroy(sys);
}

TEST(invalid_cell_index) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);

    alea_cell_info_t info;
    int rc = alea_cell_get_info(sys, 99999, &info);
    ASSERT(rc != 0);

    alea_destroy(sys);
}

TEST(invalid_surface_id) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);

    int idx = alea_surface_find(sys, 99999);
    ASSERT(idx < 0);

    alea_destroy(sys);
}

/* ========================================================================= */
/* Parse error handling                                                      */
/* ========================================================================= */

TEST(parse_empty_string) {
    mcnp_model_t* model = mcnp_load_string("", 0);
    ASSERT_NULL(model);
}

TEST(parse_title_only) {
    const char* input = "Title only\n";
    mcnp_model_t* model = mcnp_load_string(input, strlen(input));
    /* Should parse but have no cells */
    if (model) {
        ASSERT_EQ(alea_cell_count(model->sys), 0);
        mcnp_model_destroy(model);
    }
}

TEST(parse_bad_surface_type) {
    const char* input =
        "Test bad surface\n"
        "1 1 -1.0 -1\n"
        "2 0 1\n"
        "\n"
        "1 XYZZY 5.0\n"
        "\n"
        "M1 92235.80c 1.0\n";
    mcnp_model_t* model = mcnp_load_string(input, strlen(input));
    /* Should either fail to parse or handle gracefully */
    if (model) mcnp_model_destroy(model);
    /* Pass if no crash */
}

/* ========================================================================= */
/* Resource limits                                                           */
/* ========================================================================= */

TEST(empty_system_export) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);

    /* Export empty system */
    FILE* f = tmpfile();
    int rc = f ? mcnp_export_system_stream(sys, f) : -1;
    if (f) fclose(f);
    /* Should succeed with empty output or return error */

    alea_destroy(sys);
    (void)rc; /* May succeed or fail, just don't crash */
}

TEST(empty_system_queries) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);

    /* Query on empty system */
    int cell = alea_find_cell(sys, 0, 0, 0);
    ASSERT(cell < 0);

    int mat = alea_material_at(sys, 0, 0, 0);
    ASSERT_EQ(mat, 0);

    alea_destroy(sys);
}

/* ========================================================================= */
/* Error message functions                                                   */
/* ========================================================================= */

TEST(error_string) {
    const char* msg = alea_error_string(ALEA_OK);
    ASSERT_NOT_NULL(msg);

    msg = alea_error_string(ALEA_ERR_NULL_ARG);
    ASSERT_NOT_NULL(msg);

    msg = alea_error_string(ALEA_ERR_PARSE_ERROR);
    ASSERT_NOT_NULL(msg);
}

TEST(error_clear) {
    alea_set_error_detail(ALEA_ERR_NULL_ARG, "test error");
    ASSERT_EQ(alea_get_last_error(), ALEA_ERR_NULL_ARG);

    alea_clear_error_detail();
    ASSERT_EQ(alea_get_last_error(), ALEA_OK);
}

/* ========================================================================= */
/* System lifecycle                                                          */
/* ========================================================================= */

TEST(system_reset) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);

    /* Add some geometry */
    alea_sphere_surface(sys, 1, 0, 0, 0, 5.0);
    int mat = alea_add_material(sys, 0);
    ASSERT(mat >= 0);
    ASSERT_EQ(alea_material_get_id(sys, mat), 1);
    int c = alea_add_cell(sys, 0, alea_halfspace(sys, 0, -1), mat, -1.0, 0);
    ASSERT(c >= 0);
    int mat_ids[] = {1};
    double fractions[] = {1.0};
    ASSERT_EQ(alea_create_mixture(sys, mat_ids, fractions, 1, 100), 100);
    ASSERT_EQ(alea_prepare_query_acceleration(sys), 0);

    ASSERT_EQ(alea_surface_count(sys), 1);
    ASSERT_EQ(alea_cell_count(sys), 1);
    ASSERT_EQ(alea_material_count(sys), 1);
    ASSERT_EQ(alea_mixture_count(sys), 1);
    ASSERT_NOT_NULL(sys->surface_bvh);
    ASSERT_NOT_NULL(sys->spatial_index);
    ASSERT_NOT_NULL(sys->prim_to_surface);

    /* Reset should clear everything */
    alea_reset(sys);
    ASSERT_EQ(alea_surface_count(sys), 0);
    ASSERT_EQ(alea_cell_count(sys), 0);
    ASSERT_EQ(alea_material_count(sys), 0);
    ASSERT_EQ(alea_mixture_count(sys), 0);
    ASSERT_NULL(sys->surface_bvh);
    ASSERT_NULL(sys->spatial_index);
    ASSERT_NULL(sys->prim_to_surface);
    ASSERT_NULL(sys->mc_id_to_surface);
    ASSERT_NULL(sys->surface_lookup);
    ASSERT_EQ(sys->surface_lookup_size, 0);
    ASSERT_EQ(atomic_load(&sys->query_cache_state), 0u);

    int mat_after_reset = alea_add_material(sys, 0);
    ASSERT(mat_after_reset >= 0);
    ASSERT_EQ(alea_material_get_id(sys, mat_after_reset), 1);
    alea_sphere_surface(sys, 0, 0, 0, 0, 1.0);
    ASSERT_EQ(sys->surfaces.data[0].mc_surface_id, 1);
    int c_after_reset = alea_add_cell(sys, 0, alea_halfspace(sys, 0, -1),
                                      mat_after_reset, -1.0, 0);
    ASSERT(c_after_reset >= 0);
    ASSERT_EQ(sys->cells.data[c_after_reset].mc_cell_id, 1);

    alea_destroy(sys);
}

TEST(system_clone) {
    const char* input =
        "Test clone\n"
        "1 1 -1.0 -1\n"
        "2 0 1\n"
        "\n"
        "1 SO 5.0\n"
        "\n"
        "M1 92235.80c 1.0\n";
    mcnp_model_t* model = mcnp_load_string(input, strlen(input));
    ASSERT_NOT_NULL(model);
    alea_build_universe_index(model->sys);

    alea_system_t* clone = alea_clone(model->sys);
    ASSERT_NOT_NULL(clone);
    alea_build_universe_index(clone);

    /* Both should give same answers */
    ASSERT_EQ(alea_material_at(model->sys, 0, 0, 0),
              alea_material_at(clone, 0, 0, 0));
    ASSERT_EQ(alea_material_at(model->sys, 10, 0, 0),
              alea_material_at(clone, 10, 0, 0));

    alea_destroy(clone);
    mcnp_model_destroy(model);
}

/* Regression: a cloned system must carry a populated primitive dedup index,
   otherwise inserting a primitive that already exists in the clone fails to
   dedup and silently creates a duplicate. */
TEST(system_clone_primitive_dedup) {
    const char* input =
        "Test clone dedup\n"
        "1 1 -1.0 -1\n"
        "2 0 1\n"
        "\n"
        "1 SO 5.0\n"
        "\n"
        "M1 92235.80c 1.0\n";
    mcnp_model_t* model = mcnp_load_string(input, strlen(input));
    ASSERT_NOT_NULL(model);

    alea_system_t* clone = alea_clone(model->sys);
    ASSERT_NOT_NULL(clone);

    size_t before = alea_vec_count(&clone->primitives);
    ASSERT(before >= 1);  /* the SO sphere must have produced a primitive */

    /* Re-submit a primitive already present in the clone. With a populated
       index this is a dedup hit (count unchanged, existing id returned);
    with the old empty-index bug it would append a duplicate. */
    alea_primitive_type_t t = clone->primitives.data[0].type;
    alea_primitive_data_t d;
    ASSERT(alea_primitive_copy_data(clone, 0, &d));
    int8_t inverted = 0;
    alea_primitive_id_t id = alea_get_or_create_primitive(clone, t, &d, &inverted);

    ASSERT_EQ(alea_vec_count(&clone->primitives), before);  /* no duplicate */
    ASSERT(id < (alea_primitive_id_t)before);               /* reused existing */

    alea_destroy(clone);
    mcnp_model_destroy(model);
}

TEST(system_merge_primitive_payload_remap_dedup) {
    alea_system_t* target = alea_create();
    alea_system_t* source = alea_create();
    ASSERT_NOT_NULL(target);
    ASSERT_NOT_NULL(source);

    int target_surface = alea_sphere_surface(target, 0, 0, 0, 0, 5.0);
    int source_surface = alea_sphere_surface(source, 0, 0, 0, 0, 5.0);
    ASSERT(target_surface >= 0);
    ASSERT(source_surface >= 0);

    alea_node_id_t target_root = alea_halfspace(target, target_surface, -1);
    alea_node_id_t source_root = alea_halfspace(source, source_surface, -1);
    ASSERT(target_root != ALEA_NODE_ID_INVALID);
    ASSERT(source_root != ALEA_NODE_ID_INVALID);

    int target_mat = alea_add_material(target, 1);
    int source_mat = alea_add_material(source, 2);
    ASSERT(target_mat >= 0);
    ASSERT(source_mat >= 0);
    ASSERT(alea_add_cell(target, 1, target_root, target_mat, -1.0, 0) >= 0);
    ASSERT(alea_add_cell(source, 2, source_root, source_mat, -1.0, 0) >= 0);

    ASSERT_EQ(alea_vec_count(&target->primitives), 1);
    ASSERT_EQ(alea_vec_count(&source->primitives), 1);

    int added = alea_merge(target, source, 100);
    ASSERT_EQ(added, 1);
    ASSERT_EQ(alea_vec_count(&target->primitives), 1);
    ASSERT_EQ(alea_vec_count(&target->surfaces), 2);
    ASSERT_EQ(target->surfaces.data[1].primitive_id, target->surfaces.data[0].primitive_id);

    alea_primitive_data_t payload;
    ASSERT(alea_primitive_copy_data(target, target->surfaces.data[1].primitive_id, &payload));
    ASSERT_NEAR(payload.sphere.radius, 5.0, 1e-12);

    const alea_cell_entry_t* merged_cell = &target->cells.data[1];
    ASSERT_EQ(merged_cell->mc_cell_id, 102);
    ASSERT(alea_point_inside(target, merged_cell->root_node_id, 0, 0, 0));
    ASSERT_FALSE(alea_point_inside(target, merged_cell->root_node_id, 6, 0, 0));

    alea_destroy(source);
    alea_destroy(target);
}

TEST_MAIN()
