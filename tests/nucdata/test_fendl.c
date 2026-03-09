// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/**
 * @file test_fendl.c
 * @brief Test loading real FENDL-3.2c ACE data
 */

#include "alea_nucdata.h"
#include "alea_test.h"

#define FENDL_DIR "fendl-FENDL-3.2c-neutron-ace/neutron/ace"

static alea_nuc_xsdir_t* xsdir;
static alea_nuc_nuclide_t* h1;
static alea_nuc_nuclide_t* u235;
static alea_nuc_nuclide_t* fe56;

static void setup(void) {
    xsdir = alea_nuc_xsdir_load_dir(FENDL_DIR);
    if (!xsdir) return;
    h1   = alea_nuc_xsdir_get_nuclide(xsdir, "1001.32c");
    u235 = alea_nuc_xsdir_get_nuclide(xsdir, "92235.32c");
    fe56 = alea_nuc_xsdir_get_nuclide(xsdir, "26056.32c");
}

/* --- xsdir --- */

TEST(load_xsdir_dir) {
    if (!xsdir) SKIP("FENDL data not downloaded (make data-fendl)");
}

TEST(xsdir_entry_count) {
    if (!xsdir) SKIP("no data");
    size_t count = alea_nuc_xsdir_count(xsdir);
    ASSERT(count >= 100);
}

TEST(find_xsdir_u235) {
    if (!xsdir) SKIP("no data");
    const alea_nuc_xsdir_entry_t* e = alea_nuc_xsdir_find(xsdir, "92235.32c");
    ASSERT_NOT_NULL(e);
    ASSERT_NEAR(e->awr, 233.0, 1.0);
}

TEST(find_xsdir_h1) {
    if (!xsdir) SKIP("no data");
    ASSERT_NOT_NULL(alea_nuc_xsdir_find(xsdir, "1001.32c"));
}

TEST(find_xsdir_fe56) {
    if (!xsdir) SKIP("no data");
    ASSERT_NOT_NULL(alea_nuc_xsdir_find(xsdir, "26056.32c"));
}

/* --- H-1 --- */

TEST(load_h1) {
    if (!xsdir) SKIP("no data");
    ASSERT_NOT_NULL(h1);
}

TEST(h1_xs_thermal) {
    if (!h1) SKIP("no data");
    double sig_t = alea_nuc_xs_total(h1, 2.53e-8);
    ASSERT(sig_t > 10.0 && sig_t < 100.0);
}

TEST(h1_xs_1MeV) {
    if (!h1) SKIP("no data");
    double sig_t = alea_nuc_xs_total(h1, 1.0);
    ASSERT(sig_t > 1.0 && sig_t < 10.0);
}

TEST(h1_xs_14MeV) {
    if (!h1) SKIP("no data");
    double sig_t = alea_nuc_xs_total(h1, 14.0);
    ASSERT(sig_t > 0.3 && sig_t < 3.0);
}

/* --- U-235 --- */

TEST(load_u235) {
    if (!xsdir) SKIP("no data");
    ASSERT_NOT_NULL(u235);
}

TEST(u235_xs_thermal) {
    if (!u235) SKIP("no data");
    double sig_t = alea_nuc_xs_total(u235, 2.53e-8);
    ASSERT(sig_t > 400.0 && sig_t < 1200.0);
}

TEST(u235_fission_thermal) {
    if (!u235) SKIP("no data");
    double sig_f = alea_nuc_xs_reaction(u235, 18, 2.53e-8);
    ASSERT(sig_f > 300.0 && sig_f < 900.0);
}

TEST(u235_absorption_thermal) {
    if (!u235) SKIP("no data");
    double sig_a = alea_nuc_xs_absorption(u235, 2.53e-8);
    /* ACE absorption = disappearance = capture only (not fission) */
    ASSERT(sig_a > 50.0 && sig_a < 200.0);
}

TEST(u235_has_fission_data) {
    if (!u235) SKIP("no data");
    ASSERT_NOT_NULL(u235->fission);
    ASSERT_NOT_NULL(u235->fission->total);
}

TEST(u235_n_reactions) {
    if (!u235) SKIP("no data");
    ASSERT(u235->n_reactions >= 10);
}

/* --- Fe-56 --- */

TEST(load_fe56) {
    if (!xsdir) SKIP("no data");
    ASSERT_NOT_NULL(fe56);
}

/* --- Caching --- */

TEST(nuclide_caching) {
    if (!xsdir || !h1) SKIP("no data");
    alea_nuc_nuclide_t* h1_again = alea_nuc_xsdir_get_nuclide(xsdir, "1001.32c");
    ASSERT(h1_again == h1);
}

/* --- Runner --- */

/* Custom main to do setup/teardown around the auto-registered tests */
alea_test_entry_t *alea_test_list = NULL;
alea_test_entry_t **alea_test_tail = &alea_test_list;
int alea_test_passed = 0;
int alea_test_failed = 0;
int alea_test_current_failed = 0;
const char *alea_test_current_name = NULL;

int main(int argc, char **argv) {
    const char *filter = argc > 1 ? argv[1] : NULL;

    printf("FENDL-3.2c integration tests\n");
    printf("============================\n");

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
