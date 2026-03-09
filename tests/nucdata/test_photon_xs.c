// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/**
 * @file test_photon_xs.c
 * @brief Tests for photoatomic cross-section lookup from mcplib84
 */

#include "alea_nucdata.h"
#include "alea_test.h"

#define MCPLIB84_XSDIR "mcplib84/xsdir"

static alea_nuc_xsdir_t* xsdir;
static alea_nuc_nuclide_t* h;
static alea_nuc_nuclide_t* pb;

static void setup(void) {
    xsdir = alea_nuc_xsdir_load(MCPLIB84_XSDIR);
    if (!xsdir) return;
    h  = alea_nuc_xsdir_get_nuclide(xsdir, "1000.84p");
    pb = alea_nuc_xsdir_get_nuclide(xsdir, "82000.84p");
}

/* --- xsdir --- */

TEST(load_xsdir) {
    if (!xsdir) SKIP("mcplib84 not downloaded (make data-mcplib84)");
}

TEST(xsdir_has_entries) {
    if (!xsdir) SKIP("no data");
    ASSERT(alea_nuc_xsdir_count(xsdir) > 0);
}

TEST(xsdir_find_hydrogen) {
    if (!xsdir) SKIP("no data");
    ASSERT_NOT_NULL(alea_nuc_xsdir_find(xsdir, "1000.84p"));
}

TEST(xsdir_find_lead) {
    if (!xsdir) SKIP("no data");
    ASSERT_NOT_NULL(alea_nuc_xsdir_find(xsdir, "82000.84p"));
}

/* --- Hydrogen --- */

TEST(load_hydrogen) {
    if (!xsdir) SKIP("no data");
    ASSERT_NOT_NULL(h);
}

TEST(hydrogen_is_photon) {
    if (!h) SKIP("no data");
    ASSERT_EQ(h->particle, ALEA_NUC_PARTICLE_PHOTON);
}

TEST(hydrogen_has_energy_grid) {
    if (!h) SKIP("no data");
    ASSERT(h->n_energies > 0);
    ASSERT_NOT_NULL(h->energy);
}

TEST(hydrogen_has_photon_data) {
    if (!h) SKIP("no data");
    ASSERT_NOT_NULL(h->photon);
}

TEST(hydrogen_has_total_xs) {
    if (!h) SKIP("no data");
    ASSERT_NOT_NULL(h->sigma_total);
}

TEST(hydrogen_xs_1MeV) {
    if (!h) SKIP("no data");
    double sig_t = alea_nuc_xs_total(h, 1.0);
    ASSERT(sig_t > 0.05 && sig_t < 2.0);
}

TEST(hydrogen_compton_xs_1MeV) {
    if (!h || !h->photon) SKIP("no data");
    double sig_inc = alea_nuc_photon_xs_incoherent(h, 1.0);
    ASSERT(sig_inc > 0.05 && sig_inc < 1.0);
}

TEST(hydrogen_pe_xs_low_energy) {
    if (!h || !h->photon) SKIP("no data");
    double sig_pe = alea_nuc_photon_xs_photoelectric(h, 0.01);
    ASSERT(sig_pe >= 0.0);
}

/* --- Lead --- */

TEST(load_lead) {
    if (!xsdir) SKIP("no data");
    ASSERT_NOT_NULL(pb);
}

TEST(lead_has_photon_data) {
    if (!pb) SKIP("no data");
    ASSERT_NOT_NULL(pb->photon);
}

TEST(lead_xs_1MeV) {
    if (!pb) SKIP("no data");
    double sig_t = alea_nuc_xs_total(pb, 1.0);
    ASSERT(sig_t > 1.0 && sig_t < 100.0);
}

TEST(lead_pe_dominates_low_energy) {
    if (!pb || !pb->photon) SKIP("no data");
    double sig_pe = alea_nuc_photon_xs_photoelectric(pb, 0.05);
    double sig_tot = alea_nuc_xs_total(pb, 0.05);
    ASSERT(sig_tot > 0.0);
    ASSERT(sig_pe / sig_tot > 0.5);
}

TEST(lead_incoherent_form_factors) {
    if (!pb || !pb->photon) SKIP("no data");
    ASSERT(pb->photon->n_incoherent_ff > 0);
    ASSERT_NOT_NULL(pb->photon->incoherent_momentum);
    ASSERT_NOT_NULL(pb->photon->incoherent_ff);
}

TEST(lead_coherent_form_factors) {
    if (!pb || !pb->photon) SKIP("no data");
    ASSERT(pb->photon->n_coherent_ff > 0);
    ASSERT_NOT_NULL(pb->photon->coherent_momentum);
    ASSERT_NOT_NULL(pb->photon->coherent_ff);
}

TEST(lead_pair_production_5MeV) {
    if (!pb || !pb->photon) SKIP("no data");
    double sig_pp = alea_nuc_photon_xs_pair(pb, 5.0);
    ASSERT(sig_pp > 0.0);
}

/* --- Extra elements --- */

TEST(load_iron) {
    if (!xsdir) SKIP("no data");
    alea_nuc_nuclide_t* fe = alea_nuc_xsdir_get_nuclide(xsdir, "26000.84p");
    ASSERT_NOT_NULL(fe);
    ASSERT_NOT_NULL(fe->photon);
    double sig = alea_nuc_xs_total(fe, 1.0);
    ASSERT(sig > 0.0);
}

TEST(load_copper) {
    if (!xsdir) SKIP("no data");
    alea_nuc_nuclide_t* cu = alea_nuc_xsdir_get_nuclide(xsdir, "29000.84p");
    ASSERT_NOT_NULL(cu);
    ASSERT_NOT_NULL(cu->photon);
    double sig = alea_nuc_xs_total(cu, 1.0);
    ASSERT(sig > 0.0);
}

TEST(load_uranium) {
    if (!xsdir) SKIP("no data");
    alea_nuc_nuclide_t* u = alea_nuc_xsdir_get_nuclide(xsdir, "92000.84p");
    ASSERT_NOT_NULL(u);
    ASSERT_NOT_NULL(u->photon);
    double sig = alea_nuc_xs_total(u, 1.0);
    ASSERT(sig > 0.0);
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

    printf("nucdata photon tests (mcplib84)\n");
    printf("================================\n");

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
