// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/**
 * @file test_nucdata.c
 * @brief Unit tests for nucdata — ZAID parsing, energy lookup, reaction classification
 */

#include "alea_nucdata.h"
#include "alea_test.h"

/* --- ZAID parsing --- */

TEST(parse_zaid_u235) {
    int Z, A, meta;
    alea_nuc_table_type_t type;
    alea_error_t err = alea_nuc_parse_zaid("92235.80c", &Z, &A, &meta, &type);
    ASSERT_EQ(err, ALEA_OK);
    ASSERT_EQ(Z, 92);
    ASSERT_EQ(A, 235);
    ASSERT_EQ(type, ALEA_NUC_TABLE_CONTINUOUS_NEUTRON);
}

TEST(parse_zaid_h1) {
    int Z, A, meta;
    alea_nuc_table_type_t type;
    alea_error_t err = alea_nuc_parse_zaid("1001.80c", &Z, &A, &meta, &type);
    ASSERT_EQ(err, ALEA_OK);
    ASSERT_EQ(Z, 1);
    ASSERT_EQ(A, 1);
}

TEST(parse_zaid_photon) {
    int Z, A, meta;
    alea_nuc_table_type_t type;
    alea_error_t err = alea_nuc_parse_zaid("82000.04p", &Z, &A, &meta, &type);
    ASSERT_EQ(err, ALEA_OK);
    ASSERT_EQ(Z, 82);
    ASSERT_EQ(A, 0);
    ASSERT_EQ(type, ALEA_NUC_TABLE_PHOTOATOMIC);
}

TEST(parse_zaid_thermal_fails) {
    int Z, A, meta;
    alea_nuc_table_type_t type;
    /* "lwtr" is not numeric — should fail */
    alea_error_t err = alea_nuc_parse_zaid("lwtr.20t", &Z, &A, &meta, &type);
    ASSERT_NE(err, ALEA_OK);
}

TEST(parse_zaid_null) {
    int Z, A, meta;
    alea_nuc_table_type_t type;
    alea_error_t err = alea_nuc_parse_zaid(NULL, &Z, &A, &meta, &type);
    ASSERT_EQ(err, ALEA_ERR_NULL_ARG);
}

/* --- Energy lookup --- */

TEST(lookup_exact_point) {
    double grid[] = {1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1.0, 10.0, 20.0};
    double frac;
    int idx = alea_nuc_energy_lookup(grid, 8, 1e-3, &frac);
    ASSERT_EQ(idx, 2);
    ASSERT_NEAR(frac, 0.0, 0.001);
}

TEST(lookup_midpoint) {
    double grid[] = {1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1.0, 10.0, 20.0};
    double frac;
    int idx = alea_nuc_energy_lookup(grid, 8, 5.0, &frac);
    ASSERT_EQ(idx, 5);
}

TEST(lookup_below_range) {
    double grid[] = {1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1.0, 10.0, 20.0};
    double frac;
    int idx = alea_nuc_energy_lookup(grid, 8, 1e-8, &frac);
    ASSERT_EQ(idx, 0);
    ASSERT_NEAR(frac, 0.0, 1e-10);
}

TEST(lookup_above_range) {
    double grid[] = {1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1.0, 10.0, 20.0};
    double frac;
    int idx = alea_nuc_energy_lookup(grid, 8, 100.0, &frac);
    ASSERT_EQ(idx, 6);
    ASSERT_NEAR(frac, 1.0, 1e-10);
}

/* --- Reaction classification --- */

TEST(classify_elastic) {
    ASSERT_EQ(alea_nuc_reaction_classify(2), ALEA_NUC_RXN_SCATTER);
}

TEST(classify_inelastic_level) {
    ASSERT_EQ(alea_nuc_reaction_classify(51), ALEA_NUC_RXN_SCATTER);
}

TEST(classify_inelastic_continuum) {
    ASSERT_EQ(alea_nuc_reaction_classify(91), ALEA_NUC_RXN_SCATTER);
}

TEST(classify_fission) {
    ASSERT_EQ(alea_nuc_reaction_classify(18), ALEA_NUC_RXN_MULTIPLY);
}

TEST(classify_n2n) {
    ASSERT_EQ(alea_nuc_reaction_classify(16), ALEA_NUC_RXN_MULTIPLY);
}

TEST(classify_n3n) {
    ASSERT_EQ(alea_nuc_reaction_classify(17), ALEA_NUC_RXN_MULTIPLY);
}

TEST(classify_n4n) {
    ASSERT_EQ(alea_nuc_reaction_classify(37), ALEA_NUC_RXN_MULTIPLY);
}

TEST(classify_capture) {
    ASSERT_EQ(alea_nuc_reaction_classify(102), ALEA_NUC_RXN_ABSORPTION);
}

TEST(classify_n_alpha) {
    ASSERT_EQ(alea_nuc_reaction_classify(107), ALEA_NUC_RXN_ABSORPTION);
}

TEST(classify_n_nalpha) {
    ASSERT_EQ(alea_nuc_reaction_classify(22), ALEA_NUC_RXN_SCATTER);
}

TEST(classify_n_np) {
    ASSERT_EQ(alea_nuc_reaction_classify(28), ALEA_NUC_RXN_SCATTER);
}

/* --- Error strings --- */

TEST(error_string_ok) {
    const char* s = alea_error_string(ALEA_OK);
    ASSERT_NOT_NULL(s);
}

TEST(error_string_not_found) {
    const char* s = alea_error_string(ALEA_ERR_NOT_FOUND);
    ASSERT_NOT_NULL(s);
}

TEST_MAIN()
