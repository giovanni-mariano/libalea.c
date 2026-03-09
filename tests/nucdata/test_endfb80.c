// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/**
 * @file test_endfb80.c
 * @brief Tests for ENDF/B-VIII.0 (Lib80x) ACE data
 */

#include "alea_nucdata.h"
#include "alea_test.h"

#define ENDFB80_XSDIR "Lib80x/Lib80x/xsdir"

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
    h1    = alea_nuc_xsdir_get_nuclide(xsdir, "1001.80c");
    o16   = alea_nuc_xsdir_get_nuclide(xsdir, "8016.80c");
    u235  = alea_nuc_xsdir_get_nuclide(xsdir, "92235.80c");
    u238  = alea_nuc_xsdir_get_nuclide(xsdir, "92238.80c");
    fe56  = alea_nuc_xsdir_get_nuclide(xsdir, "26056.80c");
    pu239 = alea_nuc_xsdir_get_nuclide(xsdir, "94239.80c");
}

/* --- xsdir --- */

TEST(load_xsdir) {
    if (!xsdir) SKIP("Lib80x not downloaded (make data-endfb80)");
}

TEST(xsdir_entry_count) {
    if (!xsdir) SKIP("no data");
    size_t count = alea_nuc_xsdir_count(xsdir);
    /* Lib80x has ~550+ nuclides */
    ASSERT(count >= 500);
}

/* --- H-1 --- */

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
    /* For H-1, elastic is almost all of total */
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
    /* H-1 capture (n,gamma) is ~0.33 b at thermal */
    double sig_a = alea_nuc_xs_absorption(h1, 2.53e-8);
    ASSERT(sig_a > 0.1 && sig_a < 1.0);
}

/* --- O-16 --- */

TEST(load_o16) {
    if (!xsdir) SKIP("no data");
    ASSERT_NOT_NULL(o16);
    ASSERT_EQ(o16->Z, 8);
    ASSERT_EQ(o16->A, 16);
}

TEST(o16_1MeV_xs) {
    if (!o16) SKIP("no data");
    double sig_t = alea_nuc_xs_total(o16, 1.0);
    /* O-16 at 1 MeV ~8 b */
    ASSERT(sig_t > 2.0 && sig_t < 15.0);
}

/* --- U-235 --- */

TEST(load_u235) {
    if (!xsdir) SKIP("no data");
    ASSERT_NOT_NULL(u235);
    ASSERT_EQ(u235->Z, 92);
    ASSERT_EQ(u235->A, 235);
}

TEST(u235_thermal_total) {
    if (!u235) SKIP("no data");
    double sig_t = alea_nuc_xs_total(u235, 2.53e-8);
    /* U-235 thermal total ~700 b */
    ASSERT(sig_t > 400.0 && sig_t < 1200.0);
}

TEST(u235_thermal_fission) {
    if (!u235) SKIP("no data");
    double sig_f = alea_nuc_xs_reaction(u235, 18, 2.53e-8);
    /* U-235 thermal fission ~585 b */
    ASSERT(sig_f > 300.0 && sig_f < 900.0);
}

TEST(u235_thermal_capture) {
    if (!u235) SKIP("no data");
    double sig_a = alea_nuc_xs_absorption(u235, 2.53e-8);
    /* U-235 capture ~99 b */
    ASSERT(sig_a > 50.0 && sig_a < 200.0);
}

TEST(u235_fission_data) {
    if (!u235) SKIP("no data");
    ASSERT_NOT_NULL(u235->fission);
    ASSERT_NOT_NULL(u235->fission->total);
}

TEST(u235_nu_bar) {
    if (!u235) SKIP("no data");
    double nu = alea_nuc_nu_bar(u235, 2.53e-8);
    /* nu-bar ~2.43 at thermal */
    ASSERT_NEAR(nu, 2.43, 0.1);
}

TEST(u235_n_reactions) {
    if (!u235) SKIP("no data");
    ASSERT(u235->n_reactions >= 10);
}

TEST(u235_1MeV_fission) {
    if (!u235) SKIP("no data");
    double sig_f = alea_nuc_xs_reaction(u235, 18, 1.0);
    /* U-235 fission at 1 MeV ~1.2 b */
    ASSERT(sig_f > 0.5 && sig_f < 3.0);
}

/* --- U-238 --- */

TEST(load_u238) {
    if (!xsdir) SKIP("no data");
    ASSERT_NOT_NULL(u238);
    ASSERT_EQ(u238->Z, 92);
    ASSERT_EQ(u238->A, 238);
}

TEST(u238_thermal_total) {
    if (!u238) SKIP("no data");
    double sig_t = alea_nuc_xs_total(u238, 2.53e-8);
    /* U-238 thermal total ~12 b */
    ASSERT(sig_t > 5.0 && sig_t < 25.0);
}

TEST(u238_fission_threshold) {
    if (!u238) SKIP("no data");
    /* U-238 fission threshold ~1 MeV; below should be negligible */
    double sig_f_low = alea_nuc_xs_reaction(u238, 18, 0.1);
    double sig_f_high = alea_nuc_xs_reaction(u238, 18, 3.0);
    ASSERT(sig_f_low < 0.01);
    ASSERT(sig_f_high > 0.3);
}

TEST(u238_capture_resonance_6eV) {
    if (!u238) SKIP("no data");
    /* U-238 has a famous resonance at 6.67 eV; capture should be large */
    double sig_a = alea_nuc_xs_absorption(u238, 6.67e-6);
    ASSERT(sig_a > 1000.0);
}

TEST(u238_has_urr) {
    if (!u238) SKIP("no data");
    /* U-238 should have unresolved resonance data */
    ASSERT_NOT_NULL(u238->urr);
}

/* --- Fe-56 --- */

TEST(load_fe56) {
    if (!xsdir) SKIP("no data");
    ASSERT_NOT_NULL(fe56);
    ASSERT_EQ(fe56->Z, 26);
    ASSERT_EQ(fe56->A, 56);
}

TEST(fe56_1MeV_xs) {
    if (!fe56) SKIP("no data");
    double sig_t = alea_nuc_xs_total(fe56, 1.0);
    /* Fe-56 at 1 MeV ~3-4 b */
    ASSERT(sig_t > 1.0 && sig_t < 10.0);
}

TEST(fe56_has_inelastic) {
    if (!fe56) SKIP("no data");
    /* Fe-56 should have many inelastic levels */
    ASSERT(fe56->n_reactions >= 20);
}

/* --- Pu-239 --- */

TEST(load_pu239) {
    if (!xsdir) SKIP("no data");
    ASSERT_NOT_NULL(pu239);
    ASSERT_EQ(pu239->Z, 94);
    ASSERT_EQ(pu239->A, 239);
}

TEST(pu239_thermal_fission) {
    if (!pu239) SKIP("no data");
    double sig_f = alea_nuc_xs_reaction(pu239, 18, 2.53e-8);
    /* Pu-239 thermal fission ~748 b */
    ASSERT(sig_f > 400.0 && sig_f < 1200.0);
}

TEST(pu239_nu_bar) {
    if (!pu239) SKIP("no data");
    double nu = alea_nuc_nu_bar(pu239, 2.53e-8);
    /* Pu-239 nu-bar ~2.88 at thermal */
    ASSERT_NEAR(nu, 2.88, 0.15);
}

/* --- Water material --- */

TEST(water_material) {
    if (!h1 || !o16) SKIP("no data");
    alea_nuc_material_t* water = alea_nuc_material_create();
    ASSERT_NOT_NULL(water);
    alea_nuc_material_add(water, h1, 6.676e-2);
    alea_nuc_material_add(water, o16, 3.338e-2);

    double sig_t = alea_nuc_mat_xs_total(water, 1.0);
    /* Water macroscopic total at 1 MeV ~0.5-0.6 cm^-1 */
    ASSERT(sig_t > 0.3 && sig_t < 1.0);

    double mfp = alea_nuc_mean_free_path(water, 1.0);
    ASSERT(mfp > 1.0 && mfp < 3.0);

    alea_nuc_material_destroy(water);
}

/* --- Reaction sampling --- */

TEST(u235_reaction_sampling) {
    if (!u235) SKIP("no data");
    int mt;
    int rc = alea_nuc_sample_reaction(u235, 2.53e-8, 0.5, &mt);
    ASSERT(rc >= 0);
    ASSERT(mt > 0);
}

/* --- Nuclide caching --- */

TEST(nuclide_caching) {
    if (!xsdir || !h1) SKIP("no data");
    alea_nuc_nuclide_t* h1_again = alea_nuc_xsdir_get_nuclide(xsdir, "1001.80c");
    ASSERT(h1_again == h1);
}

/* --- Runner --- */

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
