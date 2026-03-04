// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/*
 * test_mcnp_cell_params_warn.c - Tests for cell parameter warnings
 *
 * Tests that:
 * - Duplicate parameters use the last value
 * - Invalid values fall back to 0
 * - Normal parsing is unaffected
 */

#include "alea_test.h"
#include "alea.h"
#include "alea_mcnp.h"
#include "mcnp/conversion/cell_conv.h"
#include <string.h>

/* ========================================================================= */
/* Helper                                                                    */
/* ========================================================================= */

static mcnp_model_t* parse_mcnp(const char* input) {
    mcnp_model_t* model = mcnp_load_string(input, strlen(input));
    if (model) alea_build_universe_index(model->sys);
    return model;
}

/* ========================================================================= */
/* Direct parse_cell_parameters tests                                        */
/* ========================================================================= */

TEST(params_imp_n_duplicate) {
    alea_cell_params_t p;
    parse_cell_parameters("IMP:N=2 IMP:N=5", &p, 99);
    ASSERT_NEAR(p.imp_n, 5.0, 1e-6);  /* last value wins */
    ASSERT_TRUE(p.has_imp_n);
}

TEST(params_imp_combined_duplicate) {
    alea_cell_params_t p;
    parse_cell_parameters("IMP:N,P=3 IMP:P=7", &p, 99);
    ASSERT_NEAR(p.imp_n, 3.0, 1e-6);
    ASSERT_NEAR(p.imp_p, 7.0, 1e-6);  /* last value wins for P */
}

TEST(params_universe_duplicate) {
    alea_cell_params_t p;
    parse_cell_parameters("U=1 U=2", &p, 99);
    ASSERT_EQ(p.universe_id, 2);
    ASSERT_TRUE(p.has_universe);
}

TEST(params_lat_duplicate) {
    alea_cell_params_t p;
    parse_cell_parameters("LAT=1 LAT=2", &p, 99);
    ASSERT_EQ(p.lat_type, 2);
    ASSERT_TRUE(p.has_lat);
}

TEST(params_vol_duplicate) {
    alea_cell_params_t p;
    parse_cell_parameters("VOL=10.0 VOL=20.0", &p, 99);
    ASSERT_NEAR(p.vol, 20.0, 1e-6);
    ASSERT_TRUE(p.has_vol);
}

TEST(params_tmp_duplicate) {
    alea_cell_params_t p;
    parse_cell_parameters("TMP=1e-8 TMP=2e-8", &p, 99);
    ASSERT_NEAR(p.tmp, 2e-8, 1e-15);
    ASSERT_TRUE(p.has_tmp);
}

TEST(params_mat_duplicate) {
    alea_cell_params_t p;
    parse_cell_parameters("MAT=1 MAT=2", &p, 99);
    ASSERT_EQ(p.material_id, 2);
    ASSERT_TRUE(p.has_mat);
}

TEST(params_rho_duplicate) {
    alea_cell_params_t p;
    parse_cell_parameters("RHO=-7.8 RHO=-8.5", &p, 99);
    ASSERT_NEAR(p.density, -8.5, 1e-6);
    ASSERT_TRUE(p.has_rho);
}

/* ========================================================================= */
/* Invalid value tests                                                       */
/* ========================================================================= */

TEST(params_imp_invalid_value) {
    alea_cell_params_t p;
    parse_cell_parameters("IMP:N=abc", &p, 99);
    /* Invalid parse: imp_n stays at default (1.0), has_imp_n not set */
    ASSERT_NEAR(p.imp_n, 1.0, 1e-6);
    ASSERT_FALSE(p.has_imp_n);
}

TEST(params_universe_invalid_value) {
    alea_cell_params_t p;
    parse_cell_parameters("U=xyz", &p, 99);
    ASSERT_EQ(p.universe_id, 0);
    ASSERT_TRUE(p.has_universe);  /* flag set even on bad parse */
}

TEST(params_lat_invalid_value) {
    alea_cell_params_t p;
    parse_cell_parameters("LAT=abc", &p, 99);
    ASSERT_EQ(p.lat_type, 0);
    ASSERT_TRUE(p.has_lat);
}

TEST(params_vol_invalid_value) {
    alea_cell_params_t p;
    parse_cell_parameters("VOL=xyz", &p, 99);
    ASSERT_NEAR(p.vol, 0.0, 1e-6);
    ASSERT_TRUE(p.has_vol);
}

TEST(params_mat_invalid_value) {
    alea_cell_params_t p;
    parse_cell_parameters("MAT=abc", &p, 99);
    ASSERT_EQ(p.material_id, 0);
    ASSERT_TRUE(p.has_mat);
}

TEST(params_rho_invalid_value) {
    alea_cell_params_t p;
    parse_cell_parameters("RHO=xyz", &p, 99);
    ASSERT_NEAR(p.density, 0.0, 1e-6);
    ASSERT_TRUE(p.has_rho);
}

/* ========================================================================= */
/* Normal parsing (regression)                                               */
/* ========================================================================= */

TEST(params_normal_all_simple) {
    alea_cell_params_t p;
    parse_cell_parameters("VOL=100 TMP=2.53e-8 PWT=0.5 NONU=1 PD=0.1 ELPT=0.001 UNC=1 BFLCL=2", &p, 1);
    ASSERT_NEAR(p.vol, 100.0, 1e-6);
    ASSERT_NEAR(p.tmp, 2.53e-8, 1e-15);
    ASSERT_NEAR(p.pwt, 0.5, 1e-6);
    ASSERT_EQ(p.nonu, 1);
    ASSERT_NEAR(p.pd, 0.1, 1e-6);
    ASSERT_NEAR(p.elpt, 0.001, 1e-6);
    ASSERT_EQ(p.unc, 1);
    ASSERT_EQ(p.bflcl, 2);
}

TEST(params_normal_imp_combined) {
    alea_cell_params_t p;
    parse_cell_parameters("IMP:N,P,E=0.5", &p, 1);
    ASSERT_NEAR(p.imp_n, 0.5, 1e-6);
    ASSERT_NEAR(p.imp_p, 0.5, 1e-6);
    ASSERT_NEAR(p.imp_e, 0.5, 1e-6);
    ASSERT_TRUE(p.has_imp_n);
    ASSERT_TRUE(p.has_imp_p);
    ASSERT_TRUE(p.has_imp_e);
}

TEST(params_normal_universe_fill_lat) {
    alea_cell_params_t p;
    parse_cell_parameters("U=5 FILL=3 LAT=1", &p, 1);
    ASSERT_EQ(p.universe_id, 5);
    ASSERT_EQ(p.fill_universe, 3);
    ASSERT_EQ(p.lat_type, 1);
    ASSERT_TRUE(p.has_universe);
    ASSERT_TRUE(p.has_fill);
    ASSERT_TRUE(p.has_lat);
}

TEST(params_normal_trcl_id) {
    alea_cell_params_t p;
    parse_cell_parameters("TRCL=5", &p, 1);
    ASSERT_EQ(p.trcl, 5);
    ASSERT_FALSE(p.trcl_inline);
    ASSERT_TRUE(p.has_trcl);
}

TEST(params_normal_mat_rho) {
    alea_cell_params_t p;
    parse_cell_parameters("MAT=3 RHO=-7.8", &p, 1);
    ASSERT_EQ(p.material_id, 3);
    ASSERT_NEAR(p.density, -7.8, 1e-6);
    ASSERT_TRUE(p.has_mat);
    ASSERT_TRUE(p.has_rho);
}

TEST(params_empty_string) {
    alea_cell_params_t p;
    int rc = parse_cell_parameters("", &p, 1);
    ASSERT_EQ(rc, 0);
    /* Defaults */
    ASSERT_NEAR(p.imp_n, 1.0, 1e-6);
    ASSERT_NEAR(p.imp_p, 1.0, 1e-6);
    ASSERT_NEAR(p.imp_e, 1.0, 1e-6);
    ASSERT_FALSE(p.has_imp_n);
    ASSERT_FALSE(p.has_universe);
}

TEST(params_null_string) {
    alea_cell_params_t p;
    int rc = parse_cell_parameters(NULL, &p, 1);
    ASSERT_EQ(rc, 0);
}

/* ========================================================================= */
/* Integration tests through mcnp_load_string                               */
/* ========================================================================= */

TEST(integration_duplicate_vol) {
    const char* input =
        "Test duplicate VOL\n"
        "1 1 -1.0 -1 VOL=50 VOL=100\n"
        "2 0 1\n"
        "\n"
        "1 SO 5.0\n"
        "\n"
        "M1 92235.80c 1.0\n";
    mcnp_model_t* model = parse_mcnp(input);
    ASSERT_NOT_NULL(model);
    mcnp_cell_params_t* cp = mcnp_cell_params(model, 0);
    ASSERT_NOT_NULL(cp);
    ASSERT_NEAR(cp->vol, 100.0, 1e-6);
    ASSERT_TRUE(cp->has_vol);
    mcnp_model_destroy(model);
}

TEST(integration_duplicate_imp) {
    const char* input =
        "Test duplicate IMP\n"
        "1 1 -1.0 -1 IMP:N=1 IMP:N=0\n"
        "2 0 1\n"
        "\n"
        "1 SO 5.0\n"
        "\n"
        "M1 92235.80c 1.0\n";
    mcnp_model_t* model = parse_mcnp(input);
    ASSERT_NOT_NULL(model);
    mcnp_cell_params_t* cp = mcnp_cell_params(model, 0);
    ASSERT_NOT_NULL(cp);
    ASSERT_NEAR(cp->imp_n, 0.0, 1e-6);
    mcnp_model_destroy(model);
}

TEST(integration_duplicate_universe) {
    const char* input =
        "Test duplicate U\n"
        "1 1 -1.0 -1 U=1 U=2\n"
        "2 0 1\n"
        "\n"
        "1 SO 5.0\n"
        "\n"
        "M1 92235.80c 1.0\n";
    mcnp_model_t* model = parse_mcnp(input);
    ASSERT_NOT_NULL(model);
    alea_cell_info_t info;
    ASSERT_EQ(alea_cell_get_info(model->sys, 0, &info), 0);
    ASSERT_EQ(info.universe_id, 2);
    mcnp_model_destroy(model);
}

TEST_MAIN()
