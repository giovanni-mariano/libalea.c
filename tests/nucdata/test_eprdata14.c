// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/**
 * @file test_eprdata14.c
 * @brief Tests for photoatomic data loading from eprdata14
 */

#include "alea_nucdata.h"
#include "alea_test.h"

#define EPRDATA14_XSDIR "eprdata14/eprdata14/xsdir"

static alea_nuc_xsdir_t* xsdir;
static alea_nuc_nuclide_t* h;
static alea_nuc_nuclide_t* pb;

static void setup(void) {
    xsdir = alea_nuc_xsdir_load(EPRDATA14_XSDIR);
    if (!xsdir) return;
    h  = alea_nuc_xsdir_get_nuclide(xsdir, "1000.14p");
    pb = alea_nuc_xsdir_get_nuclide(xsdir, "82000.14p");
}

/* --- xsdir --- */

TEST(load_xsdir) {
    if (!xsdir) SKIP("eprdata14 not downloaded (make data-eprdata14)");
}

TEST(xsdir_has_100_entries) {
    if (!xsdir) SKIP("no data");
    ASSERT_EQ(alea_nuc_xsdir_count(xsdir), 100);
}

/* --- Hydrogen (Z=1) --- */

TEST(load_hydrogen) {
    if (!xsdir) SKIP("no data");
    ASSERT_NOT_NULL(h);
}

TEST(hydrogen_is_photon) {
    if (!h) SKIP("no data");
    ASSERT_EQ(h->particle, ALEA_NUC_PARTICLE_PHOTON);
}

TEST(hydrogen_energy_grid_larger_than_mcplib84) {
    if (!h) SKIP("no data");
    ASSERT(h->n_energies > 278);
}

TEST(hydrogen_compton_1MeV) {
    if (!h) SKIP("no data");
    double sig_c = alea_nuc_photon_xs_incoherent(h, 1.0);
    ASSERT(sig_c > 0.10 && sig_c < 0.30);
}

TEST(hydrogen_total_1MeV) {
    if (!h) SKIP("no data");
    double sig_t = alea_nuc_xs_total(h, 1.0);
    ASSERT(sig_t > 0.10 && sig_t < 0.35);
}

TEST(hydrogen_coherent_1MeV_small) {
    if (!h) SKIP("no data");
    double sig_coh = alea_nuc_photon_xs_coherent(h, 1.0);
    ASSERT(sig_coh >= 0.0 && sig_coh < 0.01);
}

TEST(hydrogen_pe_10keV) {
    if (!h) SKIP("no data");
    double sig_pe = alea_nuc_photon_xs_photoelectric(h, 0.01);
    ASSERT(sig_pe >= 0.0);
}

TEST(hydrogen_pair_below_threshold) {
    if (!h) SKIP("no data");
    double sig_pp = alea_nuc_photon_xs_pair(h, 0.5);
    ASSERT(sig_pp < 0.001);
}

TEST(hydrogen_heating_1MeV) {
    if (!h) SKIP("no data");
    double heat = alea_nuc_xs_heating(h, 1.0);
    ASSERT(heat > 0.0);
}

/* --- Lead (Z=82) --- */

TEST(load_lead) {
    if (!xsdir) SKIP("no data");
    ASSERT_NOT_NULL(pb);
}

TEST(lead_energy_grid) {
    if (!pb) SKIP("no data");
    ASSERT(pb->n_energies > 0);
}

TEST(lead_total_1MeV) {
    if (!pb) SKIP("no data");
    double sig_t = alea_nuc_xs_total(pb, 1.0);
    ASSERT(sig_t > 5.0 && sig_t < 50.0);
}

TEST(lead_pe_dominates_50keV) {
    if (!pb) SKIP("no data");
    double sig_pe = alea_nuc_photon_xs_photoelectric(pb, 0.05);
    double sig_tot = alea_nuc_xs_total(pb, 0.05);
    ASSERT(sig_tot > 0.0);
    ASSERT(sig_pe / sig_tot > 0.5);
}

TEST(lead_pair_5MeV) {
    if (!pb) SKIP("no data");
    double sig_pp = alea_nuc_photon_xs_pair(pb, 5.0);
    ASSERT(sig_pp > 0.0);
}

TEST(lead_incoherent_ff) {
    if (!pb) SKIP("no data");
    ASSERT_NOT_NULL(pb->photon);
    ASSERT(pb->photon->n_incoherent_ff > 0);
}

TEST(lead_coherent_ff) {
    if (!pb) SKIP("no data");
    ASSERT_NOT_NULL(pb->photon);
    ASSERT(pb->photon->n_coherent_ff > 0);
}

TEST(lead_heating_per_collision_1MeV) {
    if (!pb) SKIP("no data");
    double hpc = alea_nuc_heating_per_collision(pb, 1.0);
    ASSERT(hpc > 0.0 && hpc < 1.0);
}

/* --- Bulk load all 100 elements --- */

TEST(load_all_100_elements) {
    if (!xsdir) SKIP("no data");
    int loaded = 0;
    for (int z = 1; z <= 100; z++) {
        char zaid[24];
        snprintf(zaid, sizeof(zaid), "%d000.14p", z);
        alea_nuc_nuclide_t* nuc = alea_nuc_xsdir_get_nuclide(xsdir, zaid);
        if (!nuc || !nuc->photon || nuc->n_energies <= 0) continue;
        double sig = alea_nuc_xs_total(nuc, 1.0);
        if (sig > 0.0) loaded++;
    }
    ASSERT_EQ(loaded, 100);
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

    printf("nucdata eprdata14 photon tests\n");
    printf("===============================\n");

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
