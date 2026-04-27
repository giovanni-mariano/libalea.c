// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/**
 * @file test_nucdata.c
 * @brief Unit tests for nucdata — ZAID parsing, energy lookup, reaction classification
 */

#include "alea_nucdata.h"
#include "alea_test.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

static double constant_spectrum(double E, void* ctx) {
    (void)E;
    (void)ctx;
    return 1.0;
}

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

TEST(lookup_rejects_single_point_grid) {
    double grid[] = {1.0};
    double frac = -1.0;
    int idx = alea_nuc_energy_lookup(grid, 1, 1.0, &frac);
    ASSERT_EQ(idx, -1);
}

TEST(lookup_rejects_nonascending_grid) {
    double grid[] = {1.0, 1.0, 2.0};
    double frac = -1.0;
    int idx = alea_nuc_energy_lookup(grid, 3, 1.5, &frac);
    ASSERT_EQ(idx, -1);
}

TEST(interp_loglog_rejects_single_point_grid) {
    double grid[] = {1.0};
    double values[] = {2.0};
    ASSERT_NEAR(alea_nuc_interp_loglog(grid, values, 1, 1.0), 0.0, 1e-12);
}

TEST(interp_loglog_nonpositive_energy_clamps) {
    double grid[] = {1.0, 10.0};
    double values[] = {2.0, 4.0};
    double xs = alea_nuc_interp_loglog(grid, values, 2, 0.0);
    ASSERT_TRUE(isfinite(xs));
    ASSERT_NEAR(xs, 2.0, 1e-12);
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

/* --- Multigroup validation --- */

TEST(mg_create_rejects_invalid_bounds) {
    double ascending[] = {1.0, 2.0, 3.0};
    double equal[] = {3.0, 2.0, 2.0};
    double negative[] = {3.0, -1.0};
    double not_finite[] = {3.0, NAN};

    ASSERT_NULL(alea_nuc_mg_create(2, ascending));
    ASSERT_NULL(alea_nuc_mg_create(2, equal));
    ASSERT_NULL(alea_nuc_mg_create(1, negative));
    ASSERT_NULL(alea_nuc_mg_create(1, not_finite));
}

TEST(mg_create_accepts_descending_bounds) {
    double bounds[] = {20.0, 1.0, 1e-5};
    alea_nuc_multigroup_t* mg = alea_nuc_mg_create(2, bounds);
    ASSERT_NOT_NULL(mg);
    alea_nuc_mg_destroy(mg);
}

TEST(mg_collapse_sums_all_fission_mts) {
    double bounds[] = {3.0, 1.0};
    double energy[] = {1.0, 2.0, 3.0};
    double total[] = {5.0, 5.0, 5.0};
    double absorb[] = {4.0, 4.0, 4.0};
    double f18[] = {1.0, 1.0, 1.0};
    double f19[] = {2.0, 2.0, 2.0};
    alea_nuc_reaction_t reactions[2];
    alea_nuc_nuclide_t nuc;

    memset(reactions, 0, sizeof(reactions));
    reactions[0].mt = 18;
    reactions[0].threshold_index = 1;
    reactions[0].n_energies = 3;
    reactions[0].xs = f18;
    reactions[1].mt = 19;
    reactions[1].threshold_index = 1;
    reactions[1].n_energies = 3;
    reactions[1].xs = f19;

    memset(&nuc, 0, sizeof(nuc));
    nuc.awr = 235.0;
    nuc.n_energies = 3;
    nuc.energy = energy;
    nuc.sigma_total = total;
    nuc.sigma_abs = absorb;
    nuc.n_reactions = 2;
    nuc.reactions = reactions;

    alea_nuc_multigroup_t* mg = alea_nuc_mg_create(1, bounds);
    ASSERT_NOT_NULL(mg);
    alea_nuc_mg_set_spectrum(mg, constant_spectrum, NULL);
    ASSERT_EQ(alea_nuc_mg_collapse(mg, &nuc), ALEA_OK);
    ASSERT_NEAR(mg->sigma_f[0], 3.0, 1e-12);
    alea_nuc_mg_destroy(mg);
}

TEST(mg_collapse_large_group_count_keeps_inelastic_transfer) {
    const int G = 1025;
    double* bounds = malloc((size_t)(G + 1) * sizeof(double));
    ASSERT_NOT_NULL(bounds);
    for (int i = 0; i <= G; i++)
        bounds[i] = (double)(G + 1 - i);

    alea_nuc_multigroup_t* mg = alea_nuc_mg_create(G, bounds);
    ASSERT_NOT_NULL(mg);
    alea_nuc_mg_set_spectrum(mg, constant_spectrum, NULL);

    double energy[] = {1.0, 1026.0};
    double total[] = {1.0, 1.0};
    double absorb[] = {0.0, 0.0};
    double rxn_xs[] = {1.0, 1.0};
    alea_nuc_reaction_t rxn;
    alea_nuc_nuclide_t nuc;

    memset(&rxn, 0, sizeof(rxn));
    rxn.mt = 91;
    rxn.threshold_index = 1;
    rxn.n_energies = 2;
    rxn.xs = rxn_xs;

    memset(&nuc, 0, sizeof(nuc));
    nuc.awr = 56.0;
    nuc.n_energies = 2;
    nuc.energy = energy;
    nuc.sigma_total = total;
    nuc.sigma_abs = absorb;
    nuc.n_reactions = 1;
    nuc.reactions = &rxn;

    ASSERT_EQ(alea_nuc_mg_collapse(mg, &nuc), ALEA_OK);

    double scatter_sum = 0.0;
    for (int i = 0; i < G * G; i++)
        scatter_sum += mg->scatter[i];
    ASSERT_TRUE(scatter_sum > 0.0);

    alea_nuc_mg_destroy(mg);
    free(bounds);
}

TEST(urr_factors_rejects_malformed_table) {
    alea_nuc_urr_t urr;
    alea_nuc_nuclide_t nuc;
    double factors[5] = {0.0};

    memset(&urr, 0, sizeof(urr));
    urr.n_energies = 1;
    urr.n_bands = 1;

    memset(&nuc, 0, sizeof(nuc));
    nuc.urr = &urr;

    ASSERT_EQ(alea_nuc_urr_factors(&nuc, 1.0, 0.5, factors), 0);
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
