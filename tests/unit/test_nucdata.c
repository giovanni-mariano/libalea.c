// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/**
 * @file test_nucdata.c
 * @brief Unit tests for nucdata — ZAID parsing, energy lookup, reaction classification
 */

#include "alea_nucdata.h"
#include "alea.h"
#include "nucdata/nuclear_internal.h"
#include "alea_test.h"
#include <math.h>
#include <limits.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static double constant_spectrum(double E, void* ctx) {
    (void)E;
    (void)ctx;
    return 1.0;
}

static int write_minimal_ascii_ace(const char* path, int threshold) {
    FILE* fp = fopen(path, "w");
    if (!fp) return 0;
    fprintf(fp, "1001.80c 1.0 2.53e-8 01/01/26\nsynthetic test fixture\n");
    for (int row = 0; row < 4; row++)
        fprintf(fp, "0 0.0 0 0.0 0 0.0 0 0.0\n");
    fprintf(fp, "18 1001 2 1 0 0 0 0\n0 0 0 0 0 0 0 0\n");
    fprintf(fp, "1 0 11 12 13 14 15 0\n");
    for (int row = 1; row < 4; row++)
        fprintf(fp, "0 0 0 0 0 0 0 0\n");
    fprintf(fp, "1 3 10 10\n2 2 8 8\n3 3 102 0\n0 1 %d 2\n2 2\n", threshold);
    return fclose(fp) == 0;
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

TEST(parse_zaid_mcnp_metastable_conventions) {
    int Z, A, meta;
    ASSERT_EQ(alea_nuc_parse_zaid("95242.80c", &Z, &A, &meta, NULL), ALEA_OK);
    ASSERT_EQ(Z, 95); ASSERT_EQ(A, 242); ASSERT_EQ(meta, 1);
    ASSERT_EQ(alea_nuc_parse_zaid("95642.80c", &Z, &A, &meta, NULL), ALEA_OK);
    ASSERT_EQ(Z, 95); ASSERT_EQ(A, 242); ASSERT_EQ(meta, 0);
    ASSERT_EQ(alea_nuc_parse_zaid("1095242.80c", &Z, &A, &meta, NULL), ALEA_OK);
    ASSERT_EQ(Z, 95); ASSERT_EQ(A, 242); ASSERT_EQ(meta, 0);
}

TEST(parse_zaid_rejects_nonnumeric_tail) {
    ASSERT_NE(alea_nuc_parse_zaid("1001x.80c", NULL, NULL, NULL, NULL), ALEA_OK);
    ASSERT_NE(alea_nuc_parse_zaid("1001.80x", NULL, NULL, NULL, NULL), ALEA_OK);
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

TEST(interp_loglog_clamps_above_grid) {
    double grid[] = {1.0, 100.0};
    double values[] = {100.0, 1.0};
    ASSERT_NEAR(alea_nuc_interp_loglog(grid, values, 2, 1000.0), 1.0, 1e-12);
}

TEST(xss_copy_rejects_overflowing_range) {
    double value = 1.0;
    alea_nuc_ace_table_t table;
    memset(&table, 0, sizeof(table));
    table.xss = &value;
    table.xss_length = 1;
    ASSERT_NULL(xss_copy(&table, INT_MAX - 7, 16));
    ASSERT_TRUE(table.decode_error);

    table.decode_error = false;
    ASSERT_NEAR(xss(&table, INT_MIN), 0.0, 0.0);
    ASSERT_TRUE(table.decode_error);
}

TEST(ace_type2_reads_direct_access_records) {
    const char* path = "nucdata_binary_fixture.tmp";
    FILE* fp = fopen(path, "wb");
    ASSERT_NOT_NULL(fp);
    char zaid[10] = {'1','0','0','1','.','8','0','c',' ',' '};
    char date[10] = "01/01/26";
    char comment[70] = {0};
    char mat[10] = {0};
    double awr = 1.0, temperature = 2.53e-8;
    ASSERT_EQ(fwrite(zaid, 1, 10, fp), 10);
    ASSERT_EQ(fwrite(&awr, sizeof(awr), 1, fp), 1);
    ASSERT_EQ(fwrite(&temperature, sizeof(temperature), 1, fp), 1);
    ASSERT_EQ(fwrite(date, 1, 10, fp), 10);
    ASSERT_EQ(fwrite(comment, 1, 70, fp), 70);
    ASSERT_EQ(fwrite(mat, 1, 10, fp), 10);
    for (int i = 0; i < 16; i++) {
        int32_t iz = i;
        double mass = (double)i;
        ASSERT_EQ(fwrite(&iz, sizeof(iz), 1, fp), 1);
        ASSERT_EQ(fwrite(&mass, sizeof(mass), 1, fp), 1);
    }
    int32_t nxs[16] = {2};
    int32_t jxs[32] = {1};
    ASSERT_EQ(fwrite(nxs, sizeof(nxs), 1, fp), 1);
    ASSERT_EQ(fwrite(jxs, sizeof(jxs), 1, fp), 1);
    for (long pos = ftell(fp); pos < 4096; pos++) fputc(0, fp);
    double xss[] = {3.0, 4.0};
    ASSERT_EQ(fwrite(xss, sizeof(xss), 1, fp), 1);
    ASSERT_EQ(fclose(fp), 0);

    alea_nuc_ace_table_t table;
    ASSERT_EQ(alea_nuc_ace_read(path, 1, 2, &table), ALEA_OK);
    ASSERT_STR_EQ(table.zaid, "1001.80c");
    ASSERT_NEAR(table.awr, 1.0, 1e-12);
    ASSERT_EQ(table.xss_length, 2);
    ASSERT_NEAR(table.xss[1], 4.0, 1e-12);
    alea_nuc_ace_free(&table);
    remove(path);
}

TEST(ace_versioned_text_header_is_explicitly_unsupported) {
    const char* path = "nucdata_v2_header.tmp";
    FILE* fp = fopen(path, "w");
    ASSERT_NOT_NULL(fp);
    ASSERT_TRUE(fputs("2.0.0\n", fp) >= 0);
    ASSERT_EQ(fclose(fp), 0);

    alea_nuc_ace_table_t table;
    ASSERT_EQ(alea_nuc_ace_read(path, 1, 1, &table), ALEA_ERR_UNSUPPORTED);
    remove(path);
}

TEST(ace_decode_converts_heating_number_to_cross_section) {
    const char* path = "nucdata_ascii_fixture.tmp";
    ASSERT_TRUE(write_minimal_ascii_ace(path, 1));
    alea_nuc_xsdir_entry_t entry;
    alea_nuc_xsdir_t xsdir;
    memset(&entry, 0, sizeof(entry));
    snprintf(entry.zaid, sizeof(entry.zaid), "1001.80c");
    snprintf(entry.filename, sizeof(entry.filename), "%s", path);
    entry.file_type = 1; entry.address = 1;
    memset(&xsdir, 0, sizeof(xsdir)); xsdir.count = 1; xsdir.entries = &entry;
    alea_nuc_nuclide_t* nuc = alea_nuc_load_nuclide(&xsdir, "1001.80c");
    ASSERT_NOT_NULL(nuc);
    ASSERT_NEAR(alea_nuc_xs_heating(nuc, 1.0), 30.0, 1e-12);
    ASSERT_NEAR(alea_nuc_heating_per_collision(nuc, 1.0), 3.0, 1e-12);
    alea_nuc_nuclide_free(nuc);
    remove(path);
}

TEST(ace_decode_rejects_invalid_reaction_threshold) {
    const char* path = "nucdata_bad_threshold.tmp";
    ASSERT_TRUE(write_minimal_ascii_ace(path, 0));
    alea_nuc_xsdir_entry_t entry;
    alea_nuc_xsdir_t xsdir;
    memset(&entry, 0, sizeof(entry));
    snprintf(entry.zaid, sizeof(entry.zaid), "1001.80c");
    snprintf(entry.filename, sizeof(entry.filename), "%s", path);
    entry.file_type = 1; entry.address = 1;
    memset(&xsdir, 0, sizeof(xsdir)); xsdir.count = 1; xsdir.entries = &entry;
    ASSERT_NULL(alea_nuc_load_nuclide(&xsdir, "1001.80c"));
    remove(path);
}

TEST(material_bridge_normalizes_fraction_and_library_suffix) {
    const char* path = "nucdata_material_fixture.tmp";
    ASSERT_TRUE(write_minimal_ascii_ace(path, 1));
    alea_nuc_xsdir_t* xsdir = calloc(1, sizeof(*xsdir));
    ASSERT_NOT_NULL(xsdir);
    xsdir->entries = calloc(1, sizeof(*xsdir->entries));
    ASSERT_NOT_NULL(xsdir->entries);
    xsdir->count = xsdir->capacity = 1;
    snprintf(xsdir->entries[0].zaid, sizeof(xsdir->entries[0].zaid), "1001.80c");
    snprintf(xsdir->entries[0].filename, sizeof(xsdir->entries[0].filename), "%s", path);
    xsdir->entries[0].file_type = 1; xsdir->entries[0].address = 1;
    xsdir->entries[0].type = ALEA_NUC_TABLE_CONTINUOUS_NEUTRON;

    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);
    int material = alea_add_material(sys, 1);
    ASSERT_TRUE(material >= 0);
    ASSERT_EQ(alea_material_add_nuclide(sys, material, 1001, "80c", 2.0), 0);
    int surface = alea_sphere_surface(sys, 1, 0.0, 0.0, 0.0, 1.0);
    int cell = alea_add_cell(sys, 1, alea_halfspace(sys, surface, -1),
                             material, 0.1, 0);
    alea_nuc_material_t* nmat = alea_nuc_material_from_cell(sys, cell, xsdir);
    ASSERT_NOT_NULL(nmat);
    ASSERT_EQ(nmat->n_components, 1);
    ASSERT_NEAR(nmat->components[0].number_density, 0.1, 1e-12);
    alea_nuc_material_destroy(nmat);
    alea_destroy(sys);
    alea_nuc_xsdir_free(xsdir);
    remove(path);
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

TEST(mg_collapse_does_not_double_count_partial_fission) {
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
    ASSERT_NEAR(mg->sigma_f[0], 1.0, 1e-12);
    alea_nuc_mg_destroy(mg);
}

TEST(mg_collapse_uses_full_group_flux_for_threshold_reaction) {
    double bounds[] = {4.0, 0.0};
    double energy[] = {0.0, 2.0, 3.0, 4.0};
    double total[] = {0.0, 0.0, 1.0, 1.0};
    double absorb[] = {0.0, 0.0, 0.0, 0.0};
    double xs[] = {0.0, 1.0, 1.0};
    alea_nuc_reaction_t rxn;
    alea_nuc_nuclide_t nuc;
    memset(&rxn, 0, sizeof(rxn));
    rxn.mt = 51; rxn.ty = 1; rxn.threshold_index = 2;
    rxn.n_energies = 3; rxn.xs = xs;
    memset(&nuc, 0, sizeof(nuc));
    nuc.awr = 56.0; nuc.n_energies = 4; nuc.energy = energy;
    nuc.sigma_total = total; nuc.sigma_abs = absorb;
    nuc.n_reactions = 1; nuc.reactions = &rxn;
    alea_nuc_multigroup_t* mg = alea_nuc_mg_create(1, bounds);
    ASSERT_NOT_NULL(mg);
    alea_nuc_mg_set_spectrum(mg, constant_spectrum, NULL);
    ASSERT_EQ(alea_nuc_mg_collapse(mg, &nuc), ALEA_OK);
    ASSERT_NEAR(mg->sigma_s[0], mg->sigma_t[0], 1e-12);
    alea_nuc_mg_destroy(mg);
}

TEST(photon_total_matches_component_sum) {
    double energy[] = {1.0, 100.0};
    double total[] = {100.0, 1.0};
    double ln_energy[] = {0.0, log(100.0)};
    double ln_values[] = {log(100.0), 0.0};
    double ln_zero[] = {-HUGE_VAL, -HUGE_VAL};
    alea_nuc_photon_data_t photon;
    alea_nuc_nuclide_t nuc;
    memset(&photon, 0, sizeof(photon));
    photon.n_energies = 2; photon.ln_energy = ln_energy;
    photon.ln_sigma_incoherent = ln_values;
    photon.ln_sigma_coherent = ln_zero;
    photon.ln_sigma_photoelectric = ln_zero;
    photon.ln_sigma_pair = ln_zero;
    memset(&nuc, 0, sizeof(nuc));
    nuc.particle = ALEA_NUC_PARTICLE_PHOTON;
    nuc.n_energies = 2; nuc.energy = energy; nuc.sigma_total = total;
    nuc.photon = &photon;
    ASSERT_NEAR(alea_nuc_xs_total(&nuc, 10.0), 10.0, 1e-12);
}

TEST(urr_loglog_uses_log_energy_fraction) {
    double energy[] = {1.0, 100.0};
    double table[] = {1, 1, 1, 1, 1, 1, 1, 100, 100, 100, 100, 100};
    double factors[5];
    alea_nuc_urr_t urr;
    alea_nuc_nuclide_t nuc;
    memset(&urr, 0, sizeof(urr));
    urr.n_energies = 2; urr.n_bands = 1; urr.interp = 5;
    urr.multiply_smooth = true;
    urr.energy = energy; urr.table = table;
    memset(&nuc, 0, sizeof(nuc)); nuc.urr = &urr;
    ASSERT_EQ(alea_nuc_urr_factors(&nuc, 10.0, 0.5, factors), 1);
    ASSERT_NEAR(factors[0], 10.0, 1e-12);
}

TEST(urr_absolute_values_are_normalized_to_factors) {
    double energy[] = {1.0, 3.0};
    double total[] = {10.0, 10.0};
    double elastic[] = {3.0, 3.0};
    double heating[] = {5.0, 5.0};
    double fission[] = {2.0, 2.0};
    double capture[] = {4.0, 4.0};
    double urr_energy[] = {1.0, 3.0};
    double table[] = {1, 20, 6, 4, 8, 10, 1, 20, 6, 4, 8, 10};
    double factors[5];
    alea_nuc_reaction_t reactions[2];
    alea_nuc_urr_t urr;
    alea_nuc_nuclide_t nuc;

    memset(reactions, 0, sizeof(reactions));
    reactions[0] = (alea_nuc_reaction_t){.mt=18, .threshold_index=1,
        .n_energies=2, .xs=fission};
    reactions[1] = (alea_nuc_reaction_t){.mt=102, .threshold_index=1,
        .n_energies=2, .xs=capture};
    memset(&urr, 0, sizeof(urr));
    urr.n_energies = 2; urr.n_bands = 1; urr.interp = 2;
    urr.multiply_smooth = false; urr.energy = urr_energy; urr.table = table;
    memset(&nuc, 0, sizeof(nuc));
    nuc.n_energies = 2; nuc.energy = energy; nuc.sigma_total = total;
    nuc.sigma_elastic = elastic; nuc.heating = heating;
    nuc.n_reactions = 2; nuc.reactions = reactions; nuc.urr = &urr;

    ASSERT_EQ(alea_nuc_urr_factors(&nuc, 2.0, 0.5, factors), 1);
    for (int i = 0; i < 5; i++) ASSERT_NEAR(factors[i], 2.0, 1e-12);
}

TEST(tabulated_tyr_uses_101_based_locator) {
    double raw[] = {0, 2, 1, 3, 2, 4};
    alea_nuc_reaction_t rxn;
    alea_nuc_nuclide_t nuc;
    memset(&rxn, 0, sizeof(rxn)); rxn.mt = 16; rxn.ty = 101;
    memset(&nuc, 0, sizeof(nuc));
    nuc.n_reactions = 1; nuc.reactions = &rxn;
    nuc.raw.xss = raw; nuc.raw.xss_length = 6; nuc.raw.jxs[10] = 1;
    ASSERT_NEAR(alea_nuc_reaction_yield(&nuc, 16, 2.0), 3.0, 1e-12);
}

TEST(angular_decode_stops_at_secondary_neutron_count) {
    double raw[] = {0, 2, 1, 3, 0, 0};
    alea_nuc_reaction_t rxn;
    alea_nuc_nuclide_t nuc;
    memset(&rxn, 0, sizeof(rxn)); rxn.mt = 102;
    memset(&nuc, 0, sizeof(nuc));
    nuc.n_reactions = 1; nuc.reactions = &rxn;
    nuc.raw.xss = raw; nuc.raw.xss_length = 6;
    nuc.raw.jxs[7] = 1; nuc.raw.jxs[8] = 2; nuc.raw.nxs[4] = 0;
    alea_nuc_decode_all_angular(&nuc);
    ASSERT_NULL(rxn.angular);
    ASSERT_FALSE(nuc.raw.decode_error);
}

TEST(watt_decode_keeps_independent_a_and_b_grids) {
    double raw[] = {
        0, 11, 10, 0, 2, 1, 3, 1, 1,
        0, 2, 1, 3, 1, 1,
        0, 1, 1, 2, 0
    };
    alea_nuc_ace_table_t table;
    memset(&table, 0, sizeof(table));
    table.xss = raw; table.xss_length = (int)(sizeof(raw) / sizeof(raw[0]));
    table.jxs[10] = 1;

    alea_nuc_energy_dist_t* ed = alea_nuc_decode_energy_dist(&table, 1);
    ASSERT_NOT_NULL(ed);
    ASSERT_EQ(ed->law, ALEA_NUC_ELAW_WATT);
    ASSERT_EQ(ed->n_temp, 2);
    ASSERT_EQ(ed->n_watt_b, 1);
    ASSERT_NEAR(ed->temp_energy[1], 3.0, 1e-12);
    ASSERT_NEAR(ed->watt_b_energy[0], 1.0, 1e-12);
    ASSERT_NEAR(ed->temp_C[0], 2.0, 1e-12);

    free(ed->nbt); free(ed->interp); free(ed->energy); free(ed->probability);
    free(ed->temp_energy); free(ed->temp_T);
    free(ed->watt_b_energy); free(ed->temp_C); free(ed);
}

TEST(doppler_preserves_high_energy_constant_and_absorption) {
    double energy[] = {1.0, 2.0, 3.0};
    double total[] = {3.0, 3.0, 3.0};
    double elastic[] = {1.0, 1.0, 1.0};
    double absorption[] = {0.0, 0.0, 0.0};
    alea_nuc_nuclide_t nuc;
    memset(&nuc, 0, sizeof(nuc));
    nuc.awr = 56.0; nuc.temperature = 2.53e-8;
    nuc.n_energies = 3; nuc.energy = energy;
    nuc.sigma_total = total; nuc.sigma_elastic = elastic; nuc.sigma_abs = absorption;
    ASSERT_EQ(alea_nuc_doppler_broaden(&nuc, 5.06e-8), ALEA_OK);
    ASSERT_NEAR(total[1], 3.0, 1e-6);
    ASSERT_NEAR(absorption[1], 0.0, 1e-12);
}

TEST(sample_reaction_ignores_partial_fission_when_total_exists) {
    double energy[] = {1.0, 3.0};
    double total[] = {4.0, 4.0};
    double elastic[] = {1.0, 1.0};
    double fission[] = {2.0, 2.0};
    double capture[] = {1.0, 1.0};
    alea_nuc_reaction_t rxn[3];
    alea_nuc_nuclide_t nuc;
    memset(rxn, 0, sizeof(rxn));
    rxn[0] = (alea_nuc_reaction_t){.mt=18,.threshold_index=1,.n_energies=2,.xs=fission};
    rxn[1] = (alea_nuc_reaction_t){.mt=19,.threshold_index=1,.n_energies=2,.xs=fission};
    rxn[2] = (alea_nuc_reaction_t){.mt=102,.threshold_index=1,.n_energies=2,.xs=capture};
    memset(&nuc, 0, sizeof(nuc));
    nuc.n_energies=2; nuc.energy=energy; nuc.sigma_total=total;
    nuc.sigma_elastic=elastic; nuc.n_reactions=3; nuc.reactions=rxn;
    int mt = 0;
    ASSERT_EQ(alea_nuc_sample_reaction(&nuc, 2.0, 0.9, &mt), 3);
    ASSERT_EQ(mt, 102);
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
