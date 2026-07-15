// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/*
 * test_mcnp_model.c - Tests for mcnp_model_t and the model-aware MCNP path
 *
 * Verifies that MCNP-specific cell parameters (importances, VOL, TMP, TRCL,
 * LIKE BUT) are correctly stored in mcnp_cell_params_t and survive
 * load → export → reload roundtrips.
 */

#include "alea_test.h"
#include "alea.h"
#include "alea_openmc.h"
#include "mcnp/mcnp_model.h"
#include "core/alea_system.h"
#include <string.h>
#include <stdio.h>
#include <math.h>

/* Boltzmann constant in MeV/K */
#define KB_MEV_PER_K  8.617333262e-11

/* ========================================================================= */
/* Helper: load MCNP string via model path                                   */
/* ========================================================================= */

static mcnp_model_t* load_model(const char* input) {
    mcnp_model_t* model = mcnp_load_string(input, strlen(input));
    if (model) alea_build_universe_index(model->sys);
    return model;
}

/* ========================================================================= */
/* Basic model creation                                                      */
/* ========================================================================= */

TEST(model_load_basic) {
    const char* input =
        "Basic model test\n"
        "1 1 -7.8 -1 IMP:N=1.0 IMP:P=0.5\n"
        "2 0 1 IMP:N=0\n"
        "\n"
        "1 SO 5.0\n"
        "\n"
        "M1 26056.80c 1.0\n";
    mcnp_model_t* model = load_model(input);
    ASSERT_NOT_NULL(model);
    ASSERT_NOT_NULL(model->sys);
    ASSERT_EQ(alea_cell_count(model->sys), 2);
    ASSERT_EQ((int)model->cell_params_count, 2);

    mcnp_model_destroy(model);
}

TEST(model_load_string_uses_explicit_length) {
    const char deck[] =
        "Explicit length test\n"
        "1 1 -7.8 -1 IMP:N=1.0\n"
        "2 0 1 IMP:N=0\n"
        "\n"
        "1 SO 5.0\n"
        "\n"
        "M1 26056.80c 1.0\n";
    const char garbage[] = "this is not a valid MCNP deck";
    char input[sizeof(deck) + sizeof(garbage)];

    memcpy(input, deck, sizeof(deck) - 1);
    memcpy(input + sizeof(deck) - 1, garbage, sizeof(garbage));

    mcnp_model_t* model = mcnp_load_string(input, sizeof(deck) - 1);
    ASSERT_NOT_NULL(model);
    ASSERT_NOT_NULL(model->sys);
    ASSERT_EQ(alea_cell_count(model->sys), 2);

    mcnp_model_destroy(model);
}

/* ========================================================================= */
/* Importance parameters                                                     */
/* ========================================================================= */

TEST(model_imp_values) {
    const char* input =
        "Importance test\n"
        "1 1 -1.0 -1 IMP:N=2.5 IMP:P=0.3 IMP:E=0.1\n"
        "2 0 1 IMP:N=0 IMP:P=0\n"
        "\n"
        "1 SO 5.0\n"
        "\n"
        "M1 92235.80c 1.0\n";
    mcnp_model_t* model = load_model(input);
    ASSERT_NOT_NULL(model);

    const mcnp_cell_params_t* p0 = mcnp_cell_params_const(model, 0);
    ASSERT_NOT_NULL(p0);
    ASSERT(p0->has_imp_n);
    ASSERT(p0->has_imp_p);
    ASSERT(p0->has_imp_e);
    ASSERT_NEAR(p0->imp_n, 2.5, 1e-6);
    ASSERT_NEAR(p0->imp_p, 0.3, 1e-6);
    ASSERT_NEAR(p0->imp_e, 0.1, 1e-6);

    const mcnp_cell_params_t* p1 = mcnp_cell_params_const(model, 1);
    ASSERT_NOT_NULL(p1);
    ASSERT(p1->has_imp_n);
    ASSERT_NEAR(p1->imp_n, 0.0, 1e-6);
    ASSERT(p1->has_imp_p);
    ASSERT_NEAR(p1->imp_p, 0.0, 1e-6);

    mcnp_model_destroy(model);
}

/* ========================================================================= */
/* VOL parameter                                                             */
/* ========================================================================= */

TEST(model_vol) {
    const char* input =
        "VOL test\n"
        "1 1 -1.0 -1 VOL=123.456\n"
        "2 0 1\n"
        "\n"
        "1 SO 5.0\n"
        "\n"
        "M1 92235.80c 1.0\n";
    mcnp_model_t* model = load_model(input);
    ASSERT_NOT_NULL(model);

    const mcnp_cell_params_t* p0 = mcnp_cell_params_const(model, 0);
    ASSERT_NOT_NULL(p0);
    ASSERT(p0->has_vol);
    ASSERT_NEAR(p0->vol, 123.456, 1e-3);

    /* Cell 2 has no VOL */
    const mcnp_cell_params_t* p1 = mcnp_cell_params_const(model, 1);
    ASSERT_NOT_NULL(p1);
    ASSERT(!p1->has_vol);

    mcnp_model_destroy(model);
}

/* ========================================================================= */
/* Temperature: MeV in params, Kelvin in core cell                           */
/* ========================================================================= */

TEST(model_temperature_conversion) {
    /* TMP=2.53e-8 MeV ≈ 293.6 K (room temperature) */
    const char* input =
        "Temperature test\n"
        "1 1 -1.0 -1 TMP=2.53e-8\n"
        "2 0 1\n"
        "\n"
        "1 SO 5.0\n"
        "\n"
        "M1 92235.80c 1.0\n";
    mcnp_model_t* model = load_model(input);
    ASSERT_NOT_NULL(model);

    /* MCNP params: MeV value preserved */
    const mcnp_cell_params_t* p0 = mcnp_cell_params_const(model, 0);
    ASSERT_NOT_NULL(p0);
    ASSERT(p0->has_tmp);
    ASSERT_NEAR(p0->tmp, 2.53e-8, 1e-12);

    /* Core cell: temperature in Kelvin */
    const alea_cell_entry_t* cell = &model->sys->cells.data[0];
    ASSERT(cell->has_temperature);
    double expected_K = 2.53e-8 / KB_MEV_PER_K;
    ASSERT_NEAR(cell->temperature, expected_K, 1.0);  /* within 1 K */

    /* Verify it's roughly room temperature */
    ASSERT(cell->temperature > 280.0 && cell->temperature < 310.0);

    /* Cell 2 has no temperature */
    const alea_cell_entry_t* cell2 = &model->sys->cells.data[1];
    ASSERT(!cell2->has_temperature);

    mcnp_model_destroy(model);
}

/* ========================================================================= */
/* TRCL: inline translation via model path                                   */
/* ========================================================================= */

TEST(model_trcl_inline) {
    const char* input =
        "TRCL inline test\n"
        "1 1 -1.0 -1 TRCL=(10 0 0)\n"
        "2 0 1\n"
        "\n"
        "1 SO 5.0\n"
        "\n"
        "M1 92235.80c 1.0\n";
    mcnp_model_t* model = load_model(input);
    ASSERT_NOT_NULL(model);

    /* Params: TRCL recorded */
    const mcnp_cell_params_t* p0 = mcnp_cell_params_const(model, 0);
    ASSERT_NOT_NULL(p0);
    ASSERT(p0->has_trcl);
    ASSERT(p0->trcl > 0);  /* inline transform registered */
    ASSERT(p0->trcl_inline);
    ASSERT_EQ(p0->trcl_count, 3);
    const mcnp_inline_transform_t* trcl =
        mcnp_model_inline_transform_const(model, p0->trcl_inline_index);
    ASSERT_NOT_NULL(trcl);
    ASSERT_EQ(trcl->count, 3);
    ASSERT_NEAR(trcl->values[0], 10.0, 1e-6);
    ASSERT_NEAR(trcl->values[1], 0.0, 1e-6);
    ASSERT_NEAR(trcl->values[2], 0.0, 1e-6);

    /* Geometry: sphere translated to (10,0,0) */
    ASSERT_EQ(alea_material_at(model->sys, 0, 0, 0), 0);   /* origin: void */
    ASSERT_EQ(alea_material_at(model->sys, 10, 0, 0), 1);   /* center: mat 1 */
    ASSERT_EQ(alea_material_at(model->sys, 10, 4, 0), 1);   /* inside */
    ASSERT_EQ(alea_material_at(model->sys, 10, 6, 0), 0);   /* outside */

    mcnp_model_destroy(model);
}

/* ========================================================================= */
/* TRCL: TRn card reference via model path                                   */
/* ========================================================================= */

TEST(model_trcl_tr_card) {
    const char* input =
        "TRCL TR card test\n"
        "1 1 -1.0 -1 TRCL=1\n"
        "2 0 1\n"
        "\n"
        "1 SO 5.0\n"
        "\n"
        "TR1 20 0 0\n"
        "M1 92235.80c 1.0\n";
    mcnp_model_t* model = load_model(input);
    ASSERT_NOT_NULL(model);

    /* Params: TRCL references TR1 */
    const mcnp_cell_params_t* p0 = mcnp_cell_params_const(model, 0);
    ASSERT_NOT_NULL(p0);
    ASSERT(p0->has_trcl);
    ASSERT_EQ(p0->trcl, 1);
    ASSERT(!p0->trcl_inline);

    /* Geometry: sphere translated to (20,0,0) */
    ASSERT_EQ(alea_material_at(model->sys, 0, 0, 0), 0);
    ASSERT_EQ(alea_material_at(model->sys, 20, 0, 0), 1);

    mcnp_model_destroy(model);
}

/* ========================================================================= */
/* LIKE BUT: parameter inheritance through model                             */
/* ========================================================================= */

TEST(model_like_but_inheritance) {
    const char* input =
        "LIKE BUT test\n"
        "1 1 -1.0 -1 IMP:N=3.0 IMP:P=2.0\n"
        "2 LIKE 1 BUT TRCL=(10 0 0)\n"
        "3 0 1\n"
        "\n"
        "1 SO 5.0\n"
        "\n"
        "M1 92235.80c 1.0\n";
    mcnp_model_t* model = load_model(input);
    ASSERT_NOT_NULL(model);

    ASSERT_EQ(alea_cell_count(model->sys), 3);

    /* Cell 2 inherits IMP from cell 1 */
    const mcnp_cell_params_t* p1 = mcnp_cell_params_const(model, 1);
    ASSERT_NOT_NULL(p1);
    ASSERT(p1->has_imp_n);
    ASSERT_NEAR(p1->imp_n, 3.0, 1e-6);
    ASSERT(p1->has_imp_p);
    ASSERT_NEAR(p1->imp_p, 2.0, 1e-6);

    /* Geometry: cell 2 at (10,0,0) */
    ASSERT_EQ(alea_material_at(model->sys, 10, 0, 0), 1);

    mcnp_model_destroy(model);
}

TEST(model_like_but_override) {
    const char* input =
        "LIKE BUT override test\n"
        "1 1 -1.0 -1 IMP:N=3.0\n"
        "2 LIKE 1 BUT IMP:N=0 MAT=2 RHO=-8.0 TRCL=(10 0 0)\n"
        "3 0 1\n"
        "\n"
        "1 SO 5.0\n"
        "\n"
        "M1 92235.80c 1.0\n"
        "M2 26056.80c 1.0\n";
    mcnp_model_t* model = load_model(input);
    ASSERT_NOT_NULL(model);

    /* Cell 2: IMP:N overridden to 0 */
    const mcnp_cell_params_t* p1 = mcnp_cell_params_const(model, 1);
    ASSERT_NOT_NULL(p1);
    ASSERT(p1->has_imp_n);
    ASSERT_NEAR(p1->imp_n, 0.0, 1e-6);

    /* Cell 2: material overridden to M2 */
    const alea_cell_entry_t* cell2 = &model->sys->cells.data[1];
    ASSERT_EQ(cell2->material_id, 2);

    mcnp_model_destroy(model);
}

/* ========================================================================= */
/* FILL transform via model                                                  */
/* ========================================================================= */

TEST(model_fill_transform) {
    const char* input =
        "FILL transform test\n"
        "1 0 -1 FILL=1 (5 0 0)\n"
        "10 1 -1.0 -2 U=1\n"
        "99 0 1\n"
        "\n"
        "1 RPP -10 10 -10 10 -10 10\n"
        "2 SO 3.0\n"
        "\n"
        "M1 92235.80c 1.0\n";
    mcnp_model_t* model = load_model(input);
    ASSERT_NOT_NULL(model);

    /* Params: fill_transform_inline recorded */
    const mcnp_cell_params_t* p0 = mcnp_cell_params_const(model, 0);
    ASSERT_NOT_NULL(p0);
    ASSERT(p0->fill_transform_inline);
    ASSERT_EQ(p0->fill_transform_count, 3);
    const mcnp_inline_transform_t* fill_tr =
        mcnp_model_inline_transform_const(model, p0->fill_transform_index);
    ASSERT_NOT_NULL(fill_tr);
    ASSERT_EQ(fill_tr->count, 3);
    ASSERT_NEAR(fill_tr->values[0], 5.0, 1e-6);

    /* Core cell: fill_transform set to registered transform ID */
    const alea_cell_entry_t* cell0 = &model->sys->cells.data[0];
    ASSERT(cell0->fill_transform > 0);

    mcnp_model_destroy(model);
}

TEST(model_inline_transform_pool_boundaries) {
    const char* input =
        "Inline transform pool test\n"
        "1 1 -1.0 -1 TRCL=(10 0 0)\n"
        "2 0 -2 FILL=1 (5 0 0)\n"
        "10 1 -1.0 -3 U=1\n"
        "99 0 1 2\n"
        "\n"
        "1 SO 5.0\n"
        "2 RPP -10 10 -10 10 -10 10\n"
        "3 SO 2.0\n"
        "\n"
        "M1 92235.80c 1.0\n";
    mcnp_model_t* model = load_model(input);
    ASSERT_NOT_NULL(model);

    const mcnp_cell_params_t* p0 = mcnp_cell_params_const(model, 0);
    const mcnp_cell_params_t* p1 = mcnp_cell_params_const(model, 1);
    ASSERT_NOT_NULL(p0);
    ASSERT_NOT_NULL(p1);
    ASSERT(p0->trcl_inline);
    ASSERT(p1->fill_transform_inline);
    ASSERT(p0->trcl_inline_index != p1->fill_transform_index);

    const mcnp_inline_transform_t* trcl =
        mcnp_model_inline_transform_const(model, p0->trcl_inline_index);
    const mcnp_inline_transform_t* fill =
        mcnp_model_inline_transform_const(model, p1->fill_transform_index);
    ASSERT_NOT_NULL(trcl);
    ASSERT_NOT_NULL(fill);
    ASSERT_EQ(trcl->count, 3);
    ASSERT_EQ(fill->count, 3);
    ASSERT_NEAR(trcl->values[0], 10.0, 1e-6);
    ASSERT_NEAR(fill->values[0], 5.0, 1e-6);

    ASSERT_NULL(mcnp_model_inline_transform_const(model, MCNP_INLINE_TRANSFORM_INVALID));
    ASSERT_NULL(mcnp_model_inline_transform_const(model, p1->fill_transform_index + 1));

    mcnp_model_destroy(model);
}

/* ========================================================================= */
/* Fatal geometry errors fail closed                                         */
/* ========================================================================= */

TEST(model_fail_closed_unresolved_complement) {
    const char* input =
        "Complement failure test\n"
        "1 1 -1.0 -1 #99\n"
        "2 0 1\n"
        "\n"
        "1 SO 5.0\n"
        "\n"
        "M1 92235.80c 1.0\n";
    mcnp_model_t* model = load_model(input);
    ASSERT_NULL(model);
}

TEST(model_fail_closed_like_unknown_template) {
    const char* input =
        "LIKE failure test\n"
        "1 1 -1.0 -1\n"
        "2 LIKE 99 BUT TRCL=(10 0 0)\n"
        "3 0 1\n"
        "\n"
        "1 SO 5.0\n"
        "\n"
        "M1 92235.80c 1.0\n";
    mcnp_model_t* model = load_model(input);
    ASSERT_NULL(model);
}

TEST(model_fail_closed_unknown_trcl_transform) {
    const char* input =
        "TRCL failure test\n"
        "1 1 -1.0 -1 TRCL=99\n"
        "2 0 1\n"
        "\n"
        "1 SO 5.0\n"
        "\n"
        "M1 92235.80c 1.0\n";
    mcnp_model_t* model = load_model(input);
    ASSERT_NULL(model);
}

TEST(model_fail_closed_singular_transform_card) {
    const char* input =
        "Singular TR failure test\n"
        "1 1 -1.0 -1\n"
        "2 0 1\n"
        "\n"
        "1 1 SO 5.0\n"
        "\n"
        "TR1 0 0 0  1 0 0  1 0 0  0 0 1\n"
        "M1 92235.80c 1.0\n";
    mcnp_model_t* model = load_model(input);
    ASSERT_NULL(model);
    ASSERT_NOT_NULL(strstr(alea_error(), "orthonormal"));
}

/* ========================================================================= */
/* Roundtrip: load → export → reload via model                               */
/* ========================================================================= */

TEST(model_roundtrip_imp) {
    const char* input =
        "Roundtrip IMP test\n"
        "1 1 -1.0 -1 IMP:N=2.5 IMP:P=0.3 VOL=100.0 TMP=2.53e-8\n"
        "2 0 1 IMP:N=0 IMP:P=0\n"
        "\n"
        "1 SO 5.0\n"
        "\n"
        "M1 92235.80c 1.0\n";
    mcnp_model_t* model = load_model(input);
    ASSERT_NOT_NULL(model);

    /* Export via model path */
    const char* tmpfile = "test_model_rt_imp.mcnp";
    int rc = mcnp_export(model, tmpfile);
    ASSERT_EQ(rc, 0);
    mcnp_model_destroy(model);

    /* Reload and verify */
    mcnp_model_t* model2 = mcnp_load(tmpfile);
    ASSERT_NOT_NULL(model2);

    const mcnp_cell_params_t* p0 = mcnp_cell_params_const(model2, 0);
    ASSERT_NOT_NULL(p0);
    ASSERT(p0->has_imp_n);
    ASSERT_NEAR(p0->imp_n, 2.5, 1e-3);
    ASSERT(p0->has_imp_p);
    ASSERT_NEAR(p0->imp_p, 0.3, 1e-3);
    ASSERT(p0->has_vol);
    ASSERT_NEAR(p0->vol, 100.0, 1e-1);
    ASSERT(p0->has_tmp);
    ASSERT_NEAR(p0->tmp, 2.53e-8, 1e-12);

    /* Graveyard preserved */
    const mcnp_cell_params_t* p1 = mcnp_cell_params_const(model2, 1);
    ASSERT_NOT_NULL(p1);
    ASSERT(p1->has_imp_n);
    ASSERT_NEAR(p1->imp_n, 0.0, 1e-6);

    /* Geometry preserved */
    ASSERT_EQ(alea_prepare_query_acceleration(model2->sys), 0);
    ASSERT_EQ(alea_material_at(model2->sys, 0, 0, 0), 1);
    ASSERT_EQ(alea_material_at(model2->sys, 7, 0, 0), 0);

    mcnp_model_destroy(model2);
    remove(tmpfile);
}

TEST(model_roundtrip_trcl) {
    const char* input =
        "Roundtrip TRCL test\n"
        "1 1 -1.0 -1 TRCL=(10 0 0)\n"
        "2 0 1\n"
        "\n"
        "1 SO 5.0\n"
        "\n"
        "M1 92235.80c 1.0\n";
    mcnp_model_t* model = load_model(input);
    ASSERT_NOT_NULL(model);

    /* Export via model path */
    const char* tmpfile = "test_model_rt_trcl.mcnp";
    int rc = mcnp_export(model, tmpfile);
    ASSERT_EQ(rc, 0);
    mcnp_model_destroy(model);

    /* Reload and verify geometry */
    mcnp_model_t* model2 = mcnp_load(tmpfile);
    ASSERT_NOT_NULL(model2);
    alea_build_universe_index(model2->sys);

    /* Sphere should still be at (10,0,0) */
    ASSERT_EQ(alea_material_at(model2->sys, 0, 0, 0), 0);
    ASSERT_EQ(alea_material_at(model2->sys, 10, 0, 0), 1);
    ASSERT_EQ(alea_material_at(model2->sys, 10, 4, 0), 1);
    ASSERT_EQ(alea_material_at(model2->sys, 10, 6, 0), 0);

    /* TRCL params preserved */
    const mcnp_cell_params_t* p0 = mcnp_cell_params_const(model2, 0);
    ASSERT_NOT_NULL(p0);
    ASSERT(p0->has_trcl);

    mcnp_model_destroy(model2);
    remove(tmpfile);
}

/* ========================================================================= */
/* Hooks: parallel array grows with cell additions                           */
/* ========================================================================= */

TEST(model_hooks_sync) {
    const char* input =
        "Hook sync test\n"
        "1 1 -1.0 -1\n"
        "2 0 1\n"
        "\n"
        "1 SO 5.0\n"
        "\n"
        "M1 92235.80c 1.0\n";
    mcnp_model_t* model = load_model(input);
    ASSERT_NOT_NULL(model);
    ASSERT_EQ((int)model->cell_params_count, 2);

    /* Add a cell programmatically — hooks should grow params */
    alea_node_id_t dummy = model->sys->cells.data[0].root_node_id;
    int idx = alea_add_cell(model->sys, 999, dummy, 0, 0.0, 0);
    ASSERT(idx >= 0);
    ASSERT_EQ((int)model->cell_params_count, 3);

    /* New cell gets default importance 1.0 */
    const mcnp_cell_params_t* p = mcnp_cell_params_const(model, (size_t)idx);
    ASSERT_NOT_NULL(p);
    ASSERT_NEAR(p->imp_n, 1.0, 1e-6);

    mcnp_model_destroy(model);
}

TEST(model_wrap_hooks_sync) {
    alea_system_t* sys = alea_system_create();
    ASSERT_NOT_NULL(sys);

    int mat_idx = alea_add_material(sys, 1);
    ASSERT(mat_idx >= 0);

    int surf_idx = alea_sphere_surface(sys, 1, 0.0, 0.0, 0.0, 5.0);
    ASSERT(surf_idx >= 0);

    alea_node_id_t root = alea_halfspace(sys, surf_idx, -1);
    ASSERT(root != ALEA_NODE_ID_INVALID);

    int idx0 = alea_add_cell(sys, 1, root, mat_idx, 1.0, 0);
    ASSERT(idx0 >= 0);

    mcnp_model_t* model = mcnp_model_wrap(sys);
    ASSERT_NOT_NULL(model);
    ASSERT_EQ((int)model->cell_params_count, 1);

    int idx1 = alea_add_cell(sys, 2, root, 0, 0.0, 0);
    ASSERT(idx1 >= 0);
    ASSERT_EQ((int)model->cell_params_count, 2);

    const mcnp_cell_params_t* p1 = mcnp_cell_params_const(model, 1);
    ASSERT_NOT_NULL(p1);
    ASSERT_NEAR(p1->imp_n, 1.0, 1e-6);
    ASSERT_NEAR(p1->imp_p, 1.0, 1e-6);
    ASSERT_NEAR(p1->imp_e, 1.0, 1e-6);

    mcnp_model_destroy(model);
    alea_system_destroy(sys);
}

/* ========================================================================= */
/* Export without model: defaults used                                       */
/* ========================================================================= */

TEST(export_without_model_uses_defaults) {
    /* Load via model path, then export via core path (no model params) */
    const char* input =
        "No-model export test\n"
        "1 1 -1.0 -1 IMP:N=5.0\n"
        "2 0 1\n"
        "\n"
        "1 SO 5.0\n"
        "\n"
        "M1 92235.80c 1.0\n";
    mcnp_model_t* src_model = mcnp_load_string(input, strlen(input));
    ASSERT_NOT_NULL(src_model);
    alea_system_t* sys = src_model->sys;

    /* Export via core path (no model attached) */
    const char* tmpfile = "test_nomodel_exp.mcnp";
    FILE* f = fopen(tmpfile, "w");
    ASSERT_NOT_NULL(f);
    int rc = mcnp_export_system_stream(sys, f);
    fclose(f);
    ASSERT_EQ(rc, 0);
    mcnp_model_destroy(src_model);

    /* Reload with model: importances should be defaults (1.0) since
       the core export path doesn't have access to the MCNP params */
    mcnp_model_t* model = mcnp_load(tmpfile);
    ASSERT_NOT_NULL(model);

    const mcnp_cell_params_t* p0 = mcnp_cell_params_const(model, 0);
    ASSERT_NOT_NULL(p0);
    /* Without model, export writes IMP:N=1.0 (default) */
    ASSERT_NEAR(p0->imp_n, 1.0, 1e-6);

    mcnp_model_destroy(model);
    remove(tmpfile);
}

/* ========================================================================= */
/* OpenMC export: works without MCNP params                                  */
/* ========================================================================= */

TEST(openmc_export_after_refactoring) {
    const char* input =
        "OpenMC test\n"
        "1 1 -1.0 -1 IMP:N=2.5 TMP=2.53e-8\n"
        "2 0 1 IMP:N=0\n"
        "\n"
        "1 SO 5.0\n"
        "\n"
        "M1 92235.80c 1.0\n";
    /* Load via model path */
    mcnp_model_t* model = load_model(input);
    ASSERT_NOT_NULL(model);

    /* Export to OpenMC (uses only core alea_system_t fields) */
    const char* tmpfile = "test_openmc_exp.xml";
    int rc = openmc_export_system(model->sys, tmpfile);
    ASSERT_EQ(rc, 0);

    /* Verify file was created and has content */
    FILE* f = fopen(tmpfile, "r");
    ASSERT_NOT_NULL(f);
    fseek(f, 0, SEEK_END);
    long size = ftell(f);
    fclose(f);
    ASSERT(size > 0);

    mcnp_model_destroy(model);
    remove(tmpfile);
    remove("test_openmc_exp.dedup");
}

TEST_MAIN()
