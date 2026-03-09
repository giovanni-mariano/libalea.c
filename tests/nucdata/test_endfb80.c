// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/**
 * @file test_endfb80.c
 * @brief Tests for ENDF/B-VIII.0 (Lib80x) ACE data
 *
 * Covers: loading, XS lookup, material mixing, nuclide/reaction sampling,
 * angular & energy distribution decoding, URR probability tables,
 * Doppler broadening, multigroup collapse, reaction yield, Kalbach-Mann.
 */

#include "alea_nucdata.h"
#include "alea_test.h"

#define ENDFB80_XSDIR "Lib80x/xsdir"

/* Simple LCG for reproducible sampling tests */
static uint64_t rng_state = 12345678901234567ULL;
static double rng(void) {
    rng_state = rng_state * 6364136223846793005ULL + 1442695040888963407ULL;
    return (double)(rng_state >> 11) / (double)(1ULL << 53);
}

static alea_nuc_xsdir_t* xsdir;
static alea_nuc_nuclide_t* h1;
static alea_nuc_nuclide_t* o16;
static alea_nuc_nuclide_t* u235;
static alea_nuc_nuclide_t* u238;
static alea_nuc_nuclide_t* fe56;
static alea_nuc_nuclide_t* pu239;

static void setup(void) {
    xsdir = alea_nuc_xsdir_load(ENDFB80_XSDIR);
    if (!xsdir) return;
    h1    = alea_nuc_xsdir_get_nuclide(xsdir, "1001.00c");
    o16   = alea_nuc_xsdir_get_nuclide(xsdir, "8016.00c");
    u235  = alea_nuc_xsdir_get_nuclide(xsdir, "92235.00c");
    u238  = alea_nuc_xsdir_get_nuclide(xsdir, "92238.00c");
    fe56  = alea_nuc_xsdir_get_nuclide(xsdir, "26056.00c");
    pu239 = alea_nuc_xsdir_get_nuclide(xsdir, "94239.00c");
}

/* ================================================================
 * XSDIR
 * ================================================================ */

TEST(load_xsdir) {
    if (!xsdir) SKIP("Lib80x not downloaded (make data-endfb80)");
}

TEST(xsdir_entry_count) {
    if (!xsdir) SKIP("no data");
    size_t count = alea_nuc_xsdir_count(xsdir);
    /* Lib80x has ~550+ nuclides */
    ASSERT(count >= 500);
}

/* ================================================================
 * H-1
 * ================================================================ */

TEST(load_h1) {
    if (!xsdir) SKIP("no data");
    ASSERT_NOT_NULL(h1);
    ASSERT_EQ(h1->Z, 1);
    ASSERT_EQ(h1->A, 1);
}

TEST(h1_thermal_xs) {
    if (!h1) SKIP("no data");
    double sig_t = alea_nuc_xs_total(h1, 2.53e-8);
    /* H-1 thermal ~20 b */
    ASSERT(sig_t > 10.0 && sig_t < 50.0);
}

TEST(h1_elastic_dominates) {
    if (!h1) SKIP("no data");
    double sig_t = alea_nuc_xs_total(h1, 1.0);
    double sig_e = alea_nuc_xs_elastic(h1, 1.0);
    ASSERT(sig_e / sig_t > 0.9);
}

TEST(h1_1MeV_xs) {
    if (!h1) SKIP("no data");
    double sig_t = alea_nuc_xs_total(h1, 1.0);
    ASSERT(sig_t > 1.0 && sig_t < 10.0);
}

TEST(h1_14MeV_xs) {
    if (!h1) SKIP("no data");
    double sig_t = alea_nuc_xs_total(h1, 14.0);
    ASSERT(sig_t > 0.3 && sig_t < 3.0);
}

TEST(h1_capture_small) {
    if (!h1) SKIP("no data");
    double sig_a = alea_nuc_xs_absorption(h1, 2.53e-8);
    ASSERT(sig_a > 0.1 && sig_a < 1.0);
}

/* ================================================================
 * O-16
 * ================================================================ */

TEST(load_o16) {
    if (!xsdir) SKIP("no data");
    ASSERT_NOT_NULL(o16);
    ASSERT_EQ(o16->Z, 8);
    ASSERT_EQ(o16->A, 16);
}

TEST(o16_1MeV_xs) {
    if (!o16) SKIP("no data");
    double sig_t = alea_nuc_xs_total(o16, 1.0);
    ASSERT(sig_t > 2.0 && sig_t < 15.0);
}

/* ================================================================
 * U-235
 * ================================================================ */

TEST(load_u235) {
    if (!xsdir) SKIP("no data");
    ASSERT_NOT_NULL(u235);
    ASSERT_EQ(u235->Z, 92);
    ASSERT_EQ(u235->A, 235);
}

TEST(u235_thermal_total) {
    if (!u235) SKIP("no data");
    double sig_t = alea_nuc_xs_total(u235, 2.53e-8);
    ASSERT(sig_t > 400.0 && sig_t < 1200.0);
}

TEST(u235_thermal_fission) {
    if (!u235) SKIP("no data");
    double sig_f = alea_nuc_xs_reaction(u235, 18, 2.53e-8);
    ASSERT(sig_f > 300.0 && sig_f < 900.0);
}

TEST(u235_thermal_capture) {
    if (!u235) SKIP("no data");
    double sig_a = alea_nuc_xs_absorption(u235, 2.53e-8);
    ASSERT(sig_a > 50.0 && sig_a < 200.0);
}

TEST(u235_fission_data) {
    if (!u235) SKIP("no data");
    ASSERT_NOT_NULL(u235->fission);
    ASSERT_NOT_NULL(u235->fission->total);
}

TEST(u235_nu_bar_thermal) {
    if (!u235) SKIP("no data");
    double nu = alea_nuc_nu_bar(u235, 2.53e-8);
    ASSERT_NEAR(nu, 2.43, 0.1);
}

TEST(u235_nu_bar_1MeV) {
    if (!u235) SKIP("no data");
    double nu = alea_nuc_nu_bar(u235, 1.0);
    ASSERT(nu > 2.0 && nu < 3.5);
}

TEST(u235_nu_bar_14MeV) {
    if (!u235) SKIP("no data");
    double nu = alea_nuc_nu_bar(u235, 14.0);
    ASSERT(nu > 3.5 && nu < 6.0);
}

TEST(u235_n_reactions) {
    if (!u235) SKIP("no data");
    ASSERT(u235->n_reactions >= 10);
}

TEST(u235_1MeV_fission) {
    if (!u235) SKIP("no data");
    double sig_f = alea_nuc_xs_reaction(u235, 18, 1.0);
    ASSERT(sig_f > 0.5 && sig_f < 3.0);
}

/* ================================================================
 * U-238
 * ================================================================ */

TEST(load_u238) {
    if (!xsdir) SKIP("no data");
    ASSERT_NOT_NULL(u238);
    ASSERT_EQ(u238->Z, 92);
    ASSERT_EQ(u238->A, 238);
}

TEST(u238_thermal_total) {
    if (!u238) SKIP("no data");
    double sig_t = alea_nuc_xs_total(u238, 2.53e-8);
    ASSERT(sig_t > 5.0 && sig_t < 25.0);
}

TEST(u238_fission_threshold) {
    if (!u238) SKIP("no data");
    double sig_f_low = alea_nuc_xs_reaction(u238, 18, 0.1);
    double sig_f_high = alea_nuc_xs_reaction(u238, 18, 3.0);
    ASSERT(sig_f_low < 0.01);
    ASSERT(sig_f_high > 0.3);
}

TEST(u238_capture_resonance_6eV) {
    if (!u238) SKIP("no data");
    double sig_a = alea_nuc_xs_absorption(u238, 6.67e-6);
    ASSERT(sig_a > 1000.0);
}

/* ================================================================
 * Fe-56
 * ================================================================ */

TEST(load_fe56) {
    if (!xsdir) SKIP("no data");
    ASSERT_NOT_NULL(fe56);
    ASSERT_EQ(fe56->Z, 26);
    ASSERT_EQ(fe56->A, 56);
}

TEST(fe56_1MeV_xs) {
    if (!fe56) SKIP("no data");
    double sig_t = alea_nuc_xs_total(fe56, 1.0);
    ASSERT(sig_t > 1.0 && sig_t < 10.0);
}

TEST(fe56_has_inelastic) {
    if (!fe56) SKIP("no data");
    ASSERT(fe56->n_reactions >= 20);
}

/* ================================================================
 * Pu-239
 * ================================================================ */

TEST(load_pu239) {
    if (!xsdir) SKIP("no data");
    ASSERT_NOT_NULL(pu239);
    ASSERT_EQ(pu239->Z, 94);
    ASSERT_EQ(pu239->A, 239);
}

TEST(pu239_thermal_fission) {
    if (!pu239) SKIP("no data");
    double sig_f = alea_nuc_xs_reaction(pu239, 18, 2.53e-8);
    ASSERT(sig_f > 400.0 && sig_f < 1200.0);
}

TEST(pu239_nu_bar) {
    if (!pu239) SKIP("no data");
    double nu = alea_nuc_nu_bar(pu239, 2.53e-8);
    ASSERT_NEAR(nu, 2.88, 0.15);
}

/* ================================================================
 * MATERIAL
 * ================================================================ */

TEST(water_material) {
    if (!h1 || !o16) SKIP("no data");
    alea_nuc_material_t* water = alea_nuc_material_create();
    ASSERT_NOT_NULL(water);
    alea_nuc_material_add(water, h1, 6.676e-2);
    alea_nuc_material_add(water, o16, 3.338e-2);

    double sig_t = alea_nuc_mat_xs_total(water, 1.0);
    ASSERT(sig_t > 0.3 && sig_t < 1.0);

    double mfp = alea_nuc_mean_free_path(water, 1.0);
    ASSERT(mfp > 1.0 && mfp < 3.0);

    alea_nuc_material_destroy(water);
}

TEST(water_macroscopic_xs_thermal) {
    if (!h1 || !o16) SKIP("no data");
    alea_nuc_material_t* water = alea_nuc_material_create();
    alea_nuc_material_add(water, h1, 0.06694);
    alea_nuc_material_add(water, o16, 0.03347);
    double sigma = alea_nuc_mat_xs_total(water, 2.53e-8);
    /* Expected ~2 cm^-1 */
    ASSERT(sigma > 1.0 && sigma < 5.0);
    alea_nuc_material_destroy(water);
}

TEST(water_mfp_thermal) {
    if (!h1 || !o16) SKIP("no data");
    alea_nuc_material_t* water = alea_nuc_material_create();
    alea_nuc_material_add(water, h1, 0.06694);
    alea_nuc_material_add(water, o16, 0.03347);
    double mfp = alea_nuc_mean_free_path(water, 2.53e-8);
    /* lambda ~ 0.5 cm */
    ASSERT(mfp > 0.1 && mfp < 2.0);
    alea_nuc_material_destroy(water);
}

TEST(sample_distance_positive) {
    if (!h1 || !o16) SKIP("no data");
    alea_nuc_material_t* water = alea_nuc_material_create();
    alea_nuc_material_add(water, h1, 0.06694);
    alea_nuc_material_add(water, o16, 0.03347);
    double d = alea_nuc_sample_distance(water, 1.0, 0.5);
    ASSERT(d > 0.0 && isfinite(d));
    alea_nuc_material_destroy(water);
}

/* ================================================================
 * NUCLIDE SAMPLING
 * ================================================================ */

TEST(sample_nuclide_water) {
    if (!h1 || !o16) SKIP("no data");
    alea_nuc_material_t* water = alea_nuc_material_create();
    alea_nuc_material_add(water, h1, 0.06694);
    alea_nuc_material_add(water, o16, 0.03347);

    int count_h = 0;
    int N = 10000;
    for (int i = 0; i < N; i++) {
        alea_nuc_nuclide_t* target = NULL;
        alea_nuc_sample_nuclide(water, 2.53e-8, rng(), &target);
        if (target == h1) count_h++;
    }
    double frac_h = (double)count_h / N;
    /* H should dominate at thermal */
    ASSERT(frac_h > 0.7 && frac_h < 1.0);
    alea_nuc_material_destroy(water);
}

/* ================================================================
 * REACTION SAMPLING
 * ================================================================ */

TEST(sample_reaction_h1_thermal) {
    if (!h1) SKIP("no data");
    int count_elastic = 0;
    int N = 10000;
    for (int i = 0; i < N; i++) {
        int mt;
        alea_nuc_sample_reaction(h1, 2.53e-8, rng(), &mt);
        if (mt == 2) count_elastic++;
    }
    double frac = (double)count_elastic / N;
    ASSERT(frac > 0.95);
}

TEST(sample_reaction_u235_thermal) {
    if (!u235) SKIP("no data");
    int count_fission = 0;
    int N = 10000;
    for (int i = 0; i < N; i++) {
        int mt;
        alea_nuc_sample_reaction(u235, 2.53e-8, rng(), &mt);
        if (mt == 18 || mt == 19 || mt == 20 || mt == 21) count_fission++;
    }
    double frac_f = (double)count_fission / N;
    /* Fission should be dominant (~80-85%) */
    ASSERT(frac_f > 0.5);
}

TEST(u235_reaction_sampling) {
    if (!u235) SKIP("no data");
    int mt;
    int rc = alea_nuc_sample_reaction(u235, 2.53e-8, 0.5, &mt);
    ASSERT(rc >= 0);
    ASSERT(mt > 0);
}

/* ================================================================
 * ANGULAR DISTRIBUTIONS
 * ================================================================ */

TEST(u235_elastic_angular) {
    if (!u235) SKIP("no data");
    ASSERT_NOT_NULL(u235->elastic_angular);
    ASSERT(u235->elastic_angular->n_energies > 0);
}

TEST(h1_elastic_angular) {
    if (!h1) SKIP("no data");
    ASSERT_NOT_NULL(h1->elastic_angular);
    ASSERT(h1->elastic_angular->n_energies > 0);
}

TEST(u235_has_angular_data) {
    if (!u235) SKIP("no data");
    int count_with_ang = 0;
    for (int i = 0; i < u235->n_reactions; i++) {
        if (u235->reactions[i].angular)
            count_with_ang++;
    }
    ASSERT(count_with_ang >= 0); /* just ensure no crash */
}

/* ================================================================
 * ENERGY DISTRIBUTIONS
 * ================================================================ */

TEST(u235_has_energy_data) {
    if (!u235) SKIP("no data");
    int count_with_edist = 0;
    for (int i = 0; i < u235->n_reactions; i++) {
        if (u235->reactions[i].energy)
            count_with_edist++;
    }
    ASSERT(count_with_edist >= 1);
}

TEST(u235_law4_has_tabular_data) {
    if (!u235) SKIP("no data");
    int found = 0;
    for (int i = 0; i < u235->n_reactions && !found; i++) {
        alea_nuc_energy_dist_t* ed = u235->reactions[i].energy;
        while (ed && !found) {
            if ((ed->law == ALEA_NUC_ELAW_CONT_TABULAR || ed->law == ALEA_NUC_ELAW_KALBACH)
                && ed->tab.n_ein > 0) {
                found = 1;
            }
            ed = ed->next;
        }
    }
    ASSERT(found);
}

TEST(u235_law61_decode) {
    if (!u235) SKIP("no data");
    int count_61 = 0;
    for (int i = 0; i < u235->n_reactions; i++) {
        alea_nuc_energy_dist_t* ed = u235->reactions[i].energy;
        while (ed) {
            if (ed->law == ALEA_NUC_ELAW_CORRELATED && ed->tab.n_ein > 0)
                count_61++;
            ed = ed->next;
        }
    }
    ASSERT(count_61 > 0);
}

/* ================================================================
 * URR PROBABILITY TABLES
 * ================================================================ */

TEST(u238_has_urr) {
    if (!u238) SKIP("no data");
    ASSERT_NOT_NULL(u238->urr);
}

TEST(u235_has_urr) {
    if (!u235) SKIP("no data");
    ASSERT_NOT_NULL(u235->urr);
    ASSERT(u235->urr->n_energies > 0);
    ASSERT(u235->urr->n_bands > 0);
}

TEST(urr_factors_in_range) {
    if (!u238 || !u238->urr) SKIP("no data");
    int ok = 1;
    double factors[5];
    for (int i = 0; i < 100 && ok; i++) {
        double xi = rng();
        int in_urr = alea_nuc_urr_factors(u238, 0.05, xi, factors);
        if (!in_urr) { ok = 0; break; }
        for (int q = 0; q < 5; q++) {
            if (factors[q] < 0.0) { ok = 0; break; }
        }
    }
    ASSERT(ok);
}

TEST(urr_factors_outside_range) {
    if (!u238 || !u238->urr) SKIP("no data");
    double factors[5];
    int in_urr = alea_nuc_urr_factors(u238, 1.0, 0.5, factors);
    ASSERT(!in_urr);
}

TEST(urr_factor_variation) {
    if (!u238 || !u238->urr) SKIP("no data");
    double f1[5], f2[5];
    alea_nuc_urr_factors(u238, 0.05, 0.01, f1);
    alea_nuc_urr_factors(u238, 0.05, 0.99, f2);
    /* Extreme bands should differ */
    ASSERT(f1[0] != f2[0]);
}

/* ================================================================
 * DOPPLER BROADENING
 * ================================================================ */

TEST(doppler_rejects_lower_temp) {
    if (!xsdir) SKIP("no data");
    alea_nuc_nuclide_t* h1_copy = alea_nuc_load_nuclide(xsdir, "1001.00c");
    if (!h1_copy) SKIP("could not load");
    alea_error_t e = alea_nuc_doppler_broaden(h1_copy, 0.0);
    ASSERT(e != ALEA_OK);
    alea_nuc_nuclide_free(h1_copy);
}

TEST(doppler_broadens_u238_resonances) {
    if (!xsdir) SKIP("no data");
    alea_nuc_nuclide_t* u = alea_nuc_load_nuclide(xsdir, "92238.00c");
    if (!u) SKIP("could not load");
    double E_res = 6.67e-6; /* 6.67 eV */
    double sig_cold = alea_nuc_xs_total(u, E_res);
    double kT_1200K = 1.034e-7;
    alea_error_t e = alea_nuc_doppler_broaden(u, kT_1200K);
    double sig_hot = alea_nuc_xs_total(u, E_res);
    ASSERT(e == ALEA_OK);
    ASSERT(sig_hot < sig_cold);
    ASSERT(sig_hot > 0.0);
    alea_nuc_nuclide_free(u);
}

TEST(doppler_total_xs_positive) {
    if (!xsdir) SKIP("no data");
    alea_nuc_nuclide_t* u = alea_nuc_load_nuclide(xsdir, "92238.00c");
    if (!u) SKIP("could not load");
    alea_nuc_doppler_broaden(u, 1.034e-7);
    int ok = 1;
    for (int i = 0; i < u->n_energies && ok; i++) {
        if (u->sigma_total[i] < 0.0) ok = 0;
    }
    ASSERT(ok);
    alea_nuc_nuclide_free(u);
}

/* ================================================================
 * MULTIGROUP
 * ================================================================ */

TEST(mg_create_destroy) {
    double bounds[] = {20.0, 6.25e-7, 1e-11};
    alea_nuc_multigroup_t* mg = alea_nuc_mg_create(2, bounds);
    ASSERT_NOT_NULL(mg);
    ASSERT_EQ(mg->n_groups, 2);
    alea_nuc_mg_destroy(mg);
}

TEST(mg_collapse_h1) {
    if (!h1) SKIP("no data");
    double bounds[] = {20.0, 6.25e-7, 1e-11};
    alea_nuc_multigroup_t* mg = alea_nuc_mg_create(2, bounds);
    alea_error_t e = alea_nuc_mg_collapse(mg, h1);
    ASSERT(e == ALEA_OK);
    ASSERT(mg->sigma_t[0] > 0.1);
    ASSERT(mg->sigma_t[1] > 1.0);
    alea_nuc_mg_destroy(mg);
}

TEST(mg_collapse_u235_fission) {
    if (!u235) SKIP("no data");
    double bounds[] = {20.0, 6.25e-7, 1e-11};
    alea_nuc_multigroup_t* mg = alea_nuc_mg_create(2, bounds);
    alea_nuc_mg_collapse(mg, u235);
    ASSERT(mg->nu_sigma_f[0] > 0.0);
    ASSERT(mg->nu_sigma_f[1] > 0.0);
    ASSERT(mg->chi[0] > 0.5); /* most fission neutrons born fast */
    alea_nuc_mg_destroy(mg);
}

TEST(mg_scatter_matrix_h1) {
    if (!h1) SKIP("no data");
    double bounds[] = {20.0, 6.25e-7, 1e-11};
    alea_nuc_multigroup_t* mg = alea_nuc_mg_create(2, bounds);
    alea_nuc_mg_collapse(mg, h1);
    double s00 = alea_nuc_mg_scatter(mg, 0, 0);
    double s01 = alea_nuc_mg_scatter(mg, 0, 1);
    ASSERT(s00 >= 0.0);
    ASSERT(s01 >= 0.0);
    ASSERT((s00 + s01) > 0.0);
    alea_nuc_mg_destroy(mg);
}

TEST(mg_adjoint_transpose) {
    if (!h1) SKIP("no data");
    double bounds[] = {20.0, 1.0, 6.25e-7, 1e-11};
    alea_nuc_multigroup_t* mg = alea_nuc_mg_create(3, bounds);
    alea_nuc_mg_collapse(mg, h1);
    int ok = 1;
    for (int i = 0; i < 3 && ok; i++)
        for (int j = 0; j < 3 && ok; j++)
            if (alea_nuc_mg_scatter_adjoint(mg, i, j) != alea_nuc_mg_scatter(mg, j, i))
                ok = 0;
    ASSERT(ok);
    alea_nuc_mg_destroy(mg);
}

TEST(mg_sample_scatter_forward) {
    if (!h1) SKIP("no data");
    double bounds4[] = {20.0, 1.0, 0.1, 0.01, 1e-5};
    alea_nuc_multigroup_t* mg = alea_nuc_mg_create(4, bounds4);
    alea_nuc_mg_collapse(mg, h1);
    int counts[4] = {0};
    for (int i = 0; i < 10000; i++) {
        int g = alea_nuc_mg_sample_scatter(mg, 0, rng(), 0);
        if (g >= 0 && g < 4) counts[g]++;
    }
    /* H-1 (A=1) can scatter to any energy: downscatter should exist */
    ASSERT(counts[0] > 0);
    ASSERT(counts[1] > 0);
    alea_nuc_mg_destroy(mg);
}

/* ================================================================
 * REACTION YIELD (TYR)
 * ================================================================ */

TEST(yield_elastic_is_1) {
    if (!h1) SKIP("no data");
    double y = alea_nuc_reaction_yield(h1, 2, 1.0);
    ASSERT_NEAR(y, 1.0, 0.01);
}

TEST(yield_capture_is_0) {
    if (!u235) SKIP("no data");
    double y = alea_nuc_reaction_yield(u235, 102, 1.0);
    ASSERT(y < 0.01);
}

TEST(yield_fission_is_nubar) {
    if (!u235) SKIP("no data");
    double y = alea_nuc_reaction_yield(u235, 18, 2.53e-8);
    double nu = alea_nuc_nu_bar(u235, 2.53e-8);
    ASSERT_NEAR(y, nu, 0.01);
}

TEST(yield_n2n_positive) {
    if (!u238) SKIP("no data");
    double y = alea_nuc_reaction_yield(u238, 16, 14.0);
    ASSERT(y >= 1.5 && y <= 2.5);
}

/* ================================================================
 * KALBACH-MANN ANGULAR
 * ================================================================ */

TEST(kalbach_reaction_exists) {
    if (!u235 && !u238 && !h1 && !o16) SKIP("no data");
    int count_44 = 0;
    alea_nuc_nuclide_t* nucs[] = {h1, o16, u235, u238};
    for (int n = 0; n < 4; n++) {
        if (!nucs[n]) continue;
        for (int i = 0; i < nucs[n]->n_reactions; i++) {
            alea_nuc_energy_dist_t* ed = nucs[n]->reactions[i].energy;
            while (ed) {
                if (ed->law == ALEA_NUC_ELAW_KALBACH && ed->tab.n_ein > 0)
                    count_44++;
                ed = ed->next;
            }
        }
    }
    ASSERT(count_44 > 0);
}

/* ================================================================
 * CACHING
 * ================================================================ */

TEST(nuclide_caching) {
    if (!xsdir || !h1) SKIP("no data");
    alea_nuc_nuclide_t* h1_again = alea_nuc_xsdir_get_nuclide(xsdir, "1001.00c");
    ASSERT(h1_again == h1);
}

/* ================================================================
 * RUNNER
 * ================================================================ */

alea_test_entry_t *alea_test_list = NULL;
alea_test_entry_t **alea_test_tail = &alea_test_list;
int alea_test_passed = 0;
int alea_test_failed = 0;
int alea_test_current_failed = 0;
const char *alea_test_current_name = NULL;

int main(int argc, char **argv) {
    const char *filter = argc > 1 ? argv[1] : NULL;

    printf("ENDF/B-VIII.0 (Lib80x) integration tests\n");
    printf("=========================================\n");

    setup();

    for (alea_test_entry_t *t = alea_test_list; t; t = t->next) {
        if (filter && strstr(t->name, filter) == NULL) continue;
        alea_test_current_failed = 0;
        alea_test_current_name = t->name;
        printf("  %-50s ", t->name);
        fflush(stdout);
        t->fn();
        if (alea_test_current_failed) {
            alea_test_failed++;
        } else {
            printf("OK\n");
            alea_test_passed++;
        }
    }

    printf("\n----------------------------------------\n");
    printf("Results: %d passed, %d failed\n", alea_test_passed, alea_test_failed);
    printf("----------------------------------------\n\n");

    alea_nuc_xsdir_free(xsdir);
    return alea_test_failed > 0 ? 1 : 0;
}
