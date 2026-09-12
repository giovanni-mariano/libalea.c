// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/**
 * @file test_nucdata.c
 * @brief Unit tests for nucdata — ZAID parsing, energy lookup, reaction classification
 */

#include "alea_nucdata.h"
#include "alea.h"
#include "rng/alea_rng_distribution.h"
#include "nucdata/nuclear_internal.h"
#include "util/alea_parallel.h"
#include "alea_test.h"
#include <math.h>
#include <limits.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

TEST(nuclear_rng_matches_core_event_and_isolates_histories) {
    alea_nuc_rng_t rng, other, replay;
    ASSERT_EQ(alea_nuc_rng_init(&rng, 42, 7, 3, 9,
                               ALEA_NUC_RNG_COLLISION), ALEA_OK);
    replay = rng;
    alea_rng_event_t core;
    ASSERT_EQ(alea_rng_event_init(&core, ALEA_RNG_PHILOX4X32_10, 42,
        ALEA_RNG_DOMAIN_TRANSPORT_REACTION,
        alea_rng_transport_entity_id(7, 3), 9), 0);
    ASSERT_EQ(alea_nuc_rng_init(&other, 42, 8, 0, 9,
                               ALEA_NUC_RNG_FLIGHT), ALEA_OK);
    for (int i = 0; i < 20; i++) {
        double expected;
        ASSERT_EQ(alea_rng_event_next_uniform(&core, &expected), 0);
        ASSERT_EQ(alea_nuc_rng_uniform(&rng), expected);
        (void)alea_nuc_rng_uniform(&other);
        ASSERT_EQ(alea_nuc_rng_uniform(&replay), expected);
    }
    rng.local_draw = UINT64_C(0x100000000) - 2;
    ASSERT_TRUE(isfinite(alea_nuc_rng_uniform(&rng)));
    ASSERT_TRUE(isnan(alea_nuc_rng_uniform(&rng)));
    ASSERT_EQ(rng.local_draw, UINT64_C(0x100000000));
    ASSERT_TRUE(isnan(alea_nuc_rng_uniform(NULL)));
}

static double constant_spectrum(double E, void* ctx) {
    (void)E;
    (void)ctx;
    return 1.0;
}

typedef struct {
    const double* values;
    int count;
    int position;
} sequence_rng_t;

static double sequence_rng(void* context) {
    sequence_rng_t* rng = context;
    if (!rng || rng->position >= rng->count) return NAN;
    return rng->values[rng->position++];
}

typedef struct {
    alea_error_t evaluation_status;
    alea_error_t collision_status;
    alea_nuc_collision_result_t result;
} scheduled_collision_t;

typedef struct {
    const alea_nuc_prepared_material_t* prepared;
    scheduled_collision_t* collisions;
} scheduled_collision_context_t;

static int run_scheduled_collisions(void* opaque, size_t worker_index,
                                    size_t begin, size_t end) {
    scheduled_collision_context_t* context = opaque;
    (void)worker_index;
    for (size_t history = begin; history < end; history++) {
        scheduled_collision_t* output = &context->collisions[history];
        memset(output, 0, sizeof(*output));
        alea_nuc_particle_state_t incident = {
            ALEA_NUC_PARTICLE_NEUTRON, 2.0, {0.0, 0.0, 1.0}, 1.0, 0.0
        };
        alea_nuc_evaluation_t evaluation;
        output->evaluation_status =
            alea_nuc_evaluate(context->prepared, &incident, &evaluation);
        if (output->evaluation_status != ALEA_OK) continue;
        alea_nuc_rng_t rng;
        output->collision_status = alea_nuc_rng_init(
            &rng, UINT64_C(0x123456789abcdef0), history, 0, 0,
            ALEA_NUC_RNG_COLLISION);
        if (output->collision_status != ALEA_OK) continue;
        output->collision_status = alea_nuc_collide(
            &evaluation, alea_nuc_rng_uniform, &rng, &output->result);
    }
    return 0;
}

static int scheduled_collisions_equal(const scheduled_collision_t* a,
                                      const scheduled_collision_t* b) {
    const alea_nuc_collision_result_t* x = &a->result;
    const alea_nuc_collision_result_t* y = &b->result;
    if (a->evaluation_status != b->evaluation_status ||
        a->collision_status != b->collision_status ||
        x->outcome != y->outcome ||
        x->component_index != y->component_index || x->mt != y->mt ||
        x->deposition_available != y->deposition_available ||
        x->outgoing.type != y->outgoing.type || x->n_emitted != y->n_emitted)
        return 0;
    const double* xd[] = {
        &x->mu_cm, &x->mu_lab, &x->local_energy_deposition,
        &x->outgoing.energy, &x->outgoing.direction[0],
        &x->outgoing.direction[1], &x->outgoing.direction[2],
        &x->outgoing.weight, &x->outgoing.time
    };
    const double* yd[] = {
        &y->mu_cm, &y->mu_lab, &y->local_energy_deposition,
        &y->outgoing.energy, &y->outgoing.direction[0],
        &y->outgoing.direction[1], &y->outgoing.direction[2],
        &y->outgoing.weight, &y->outgoing.time
    };
    for (size_t i = 0; i < sizeof(xd) / sizeof(xd[0]); i++)
        if (memcmp(xd[i], yd[i], sizeof(double)) != 0) return 0;
    return 1;
}

typedef struct { uint64_t state; } nucdata_rng_t;

static double nucdata_rng(void* context) {
    nucdata_rng_t* rng = context;
    uint64_t x = rng->state;
    x ^= x >> 12;
    x ^= x << 25;
    x ^= x >> 27;
    rng->state = x;
    return (double)((x * UINT64_C(2685821657736338717)) >> 11) * 0x1p-53;
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

static int write_minimal_photoatomic_ace(const char* path) {
    FILE* fp = fopen(path, "w");
    if (!fp) return 0;
    fprintf(fp, "92000.31p 235.984 0.0 01/01/26\nsynthetic photon fixture\n");
    for (int row = 0; row < 4; row++)
        fprintf(fp, "0 0.0 0 0.0 0 0.0 0 0.0\n");
    fprintf(fp, "167 92 2 6 0 0 0 0\n0 0 0 0 0 0 0 0\n");
    fprintf(fp, "1 11 32 142 166 0 0 0\n");
    for (int row = 1; row < 4; row++)
        fprintf(fp, "0 0 0 0 0 0 0 0\n");

    double xss[167] = {0};
    xss[0] = log(0.001); xss[1] = log(20.0);
    for (int block = 1; block < 4; block++) {
        xss[2 * block] = log(1.0);
        xss[2 * block + 1] = log(1.0);
    }
    for (int i = 0; i < 21; i++) xss[10 + i] = (double)i;
    for (int i = 0; i < 55; i++) {
        xss[31 + i] = (double)i;
        xss[31 + 55 + i] = (double)i;
    }
    const double fluorescence[24] = {
        0.02, 0.02, 0.12, 0.12, 0.12, 0.12,
        0.2, 1.0, 2.0, 3.0, 4.0, 5.0,
        0.0, 0.3, 1.0, 2.0, 3.0, 4.0,
        0.0, 0.015, 0.09, 0.10, 0.11, 0.115
    };
    memcpy(&xss[141], fluorescence, sizeof(fluorescence));
    xss[165] = 0.5; xss[166] = 0.5;
    for (int i = 0; i < 167; i++)
        fprintf(fp, "%.17g%c", xss[i], (i % 4 == 3) ? '\n' : ' ');
    if (167 % 4) fputc('\n', fp);
    return fclose(fp) == 0;
}

static int write_minimal_epr_ace(const char* path) {
    FILE* fp = fopen(path, "w");
    if (!fp) return 0;
    fprintf(fp, "82000.14p 205.42 0.0 01/01/26\nsynthetic EPR fixture\n");
    for (int row = 0; row < 4; row++)
        fprintf(fp, "0 0.0 0 0.0 0 0.0 0 0.0\n");
    fprintf(fp, "70 82 2 0 1 3 3 0\n0 0 0 0 2 2 0 0\n");
    fprintf(fp, "1 11 15 21 21 23 24 25\n26 27 35 38 41 44 47 50\n");
    fprintf(fp, "56 59 0 0 0 0 0 0\n0 0 0 0 0 0 0 0\n");

    double xss[70] = {0};
    xss[0] = log(0.1); xss[1] = log(2.0);
    xss[2] = xss[3] = -HUGE_VAL;
    xss[4] = xss[5] = -HUGE_VAL;
    xss[6] = xss[7] = log(2.0);
    xss[10] = 0.0; xss[11] = 1.0;
    xss[12] = 0.0; xss[13] = 1.0;
    xss[14] = 0.0; xss[15] = 1.0;
    xss[16] = 0.0; xss[17] = 1.0;
    xss[18] = 1.0; xss[19] = 0.0;
    xss[20] = xss[21] = 0.5;
    xss[22] = 82; xss[23] = 0.005; xss[24] = 1.0; xss[25] = 1;
    xss[26] = 2; xss[27] = 2;
    xss[28] = 0.0; xss[29] = 1.0;
    xss[30] = xss[31] = 1.0;
    xss[32] = 0.0; xss[33] = 1.0;
    xss[34] = 1; xss[35] = 2; xss[36] = 3;
    xss[37] = 2; xss[38] = 2; xss[39] = 4;
    xss[40] = 0.1; xss[41] = 0.02; xss[42] = 0.005;
    xss[43] = 0.2; xss[44] = 0.7; xss[45] = 1.0;
    xss[46] = 2; xss[47] = 1; xss[48] = 0;
    xss[49] = xss[50] = log(2.0);
    xss[55] = 0; xss[56] = 8; xss[57] = 12;
    const double transitions[] = {
        2, 0, 0.080, 0.5, 2, 3, 0.070, 1.0,
        3, 0, 0.015, 1.0
    };
    memcpy(&xss[58], transitions, sizeof(transitions));
    for (int i = 0; i < 70; i++)
        fprintf(fp, "%.17g%c", xss[i], (i % 4 == 3) ? '\n' : ' ');
    if (70 % 4) fputc('\n', fp);
    return fclose(fp) == 0;
}

typedef struct {
    size_t fail_at;
    size_t calls;
} nuc_alloc_failure_t;

static bool fail_nuc_allocation(void* context) {
    nuc_alloc_failure_t* failure = context;
    return failure->calls++ == failure->fail_at;
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

TEST(xsdir_temperature_selection_uses_family_type_and_closest_table) {
    alea_nuc_xsdir_entry_t entries[5];
    memset(entries, 0, sizeof(entries));
    snprintf(entries[0].zaid, sizeof(entries[0].zaid), "1001.80c");
    entries[0].type = ALEA_NUC_TABLE_CONTINUOUS_NEUTRON;
    entries[0].temperature = 2.5300e-8;
    snprintf(entries[1].zaid, sizeof(entries[1].zaid), "1001.81c");
    entries[1].type = ALEA_NUC_TABLE_CONTINUOUS_NEUTRON;
    entries[1].temperature = 5.1700e-8;
    snprintf(entries[2].zaid, sizeof(entries[2].zaid), "1002.80c");
    entries[2].type = ALEA_NUC_TABLE_CONTINUOUS_NEUTRON;
    entries[2].temperature = 5.1600e-8;
    snprintf(entries[3].zaid, sizeof(entries[3].zaid), "1001.80p");
    entries[3].type = ALEA_NUC_TABLE_PHOTOATOMIC;
    entries[3].temperature = 5.1650e-8;
    snprintf(entries[4].zaid, sizeof(entries[4].zaid), "1001.bad");
    entries[4].temperature = 5.1650e-8;
    alea_nuc_xsdir_t xsdir;
    memset(&xsdir, 0, sizeof(xsdir));
    xsdir.entries = entries; xsdir.count = 5;

    const alea_nuc_xsdir_entry_t* selected = (const void*)1;
    ASSERT_EQ(alea_nuc_xsdir_find_temperature(
                  &xsdir, "1001.99c", 5.1650e-8, 1.0e-10, &selected),
              ALEA_OK);
    ASSERT_TRUE(selected == &entries[1]);
}

TEST(xsdir_temperature_selection_supports_symbolic_thermal_names) {
    alea_nuc_xsdir_entry_t entries[2];
    memset(entries, 0, sizeof(entries));
    snprintf(entries[0].zaid, sizeof(entries[0].zaid), "lwtr.20t");
    entries[0].type = ALEA_NUC_TABLE_THERMAL_SAB;
    entries[0].temperature = 2.5300e-8;
    snprintf(entries[1].zaid, sizeof(entries[1].zaid), "lwtr.21t");
    entries[1].type = ALEA_NUC_TABLE_THERMAL_SAB;
    entries[1].temperature = 5.1700e-8;
    alea_nuc_xsdir_t xsdir;
    memset(&xsdir, 0, sizeof(xsdir));
    xsdir.entries = entries; xsdir.count = 2;

    const alea_nuc_xsdir_entry_t* selected = NULL;
    ASSERT_EQ(alea_nuc_xsdir_find_temperature(
                  &xsdir, "lwtr.99t", 2.53e-8, 0.0, &selected), ALEA_OK);
    ASSERT_TRUE(selected == &entries[0]);
}

TEST(xsdir_temperature_selection_fails_closed) {
    alea_nuc_xsdir_entry_t entries[2];
    memset(entries, 0, sizeof(entries));
    snprintf(entries[0].zaid, sizeof(entries[0].zaid), "1001.80c");
    entries[0].temperature = 2.0e-8;
    snprintf(entries[1].zaid, sizeof(entries[1].zaid), "1001.81c");
    entries[1].temperature = 4.0e-8;
    alea_nuc_xsdir_t xsdir;
    memset(&xsdir, 0, sizeof(xsdir));
    xsdir.entries = entries; xsdir.count = 2;

    const alea_nuc_xsdir_entry_t* selected = (const void*)1;
    ASSERT_EQ(alea_nuc_xsdir_find_temperature(
                  &xsdir, "1001.99c", 3.0e-8, 2.0e-8, &selected),
              ALEA_ERR_INVALID_STATE);
    ASSERT_NULL(selected);
    selected = (const void*)1;
    ASSERT_EQ(alea_nuc_xsdir_find_temperature(
                  &xsdir, "1001.99c", 7.0e-8, 1.0e-9, &selected),
              ALEA_ERR_NOT_FOUND);
    ASSERT_NULL(selected);
    ASSERT_EQ(alea_nuc_xsdir_find_temperature(
                  &xsdir, "1001", 2.0e-8, 0.0, &selected),
              ALEA_ERR_INVALID_ARG);
    ASSERT_EQ(alea_nuc_xsdir_find_temperature(
                  &xsdir, "1001.80c", NAN, 0.0, &selected),
              ALEA_ERR_INVALID_ARG);
}

TEST(xsdir_temperature_bracket_is_bounded_and_unambiguous) {
    alea_nuc_xsdir_entry_t entries[4];
    memset(entries, 0, sizeof(entries));
    for (int i = 0; i < 3; i++) {
        snprintf(entries[i].zaid, sizeof(entries[i].zaid),
                 "1001.%02dc", 80 + i);
        entries[i].type = ALEA_NUC_TABLE_CONTINUOUS_NEUTRON;
    }
    entries[0].temperature = 2.0e-8;
    entries[1].temperature = 4.0e-8;
    entries[2].temperature = 8.0e-8;
    snprintf(entries[3].zaid, sizeof(entries[3].zaid), "1002.80c");
    entries[3].temperature = 4.0e-8;
    alea_nuc_xsdir_t xsdir = {0};
    xsdir.entries = entries;
    xsdir.count = 4;

    const alea_nuc_xsdir_entry_t *lower = NULL, *upper = NULL;
    double fraction = -1.0;
    ASSERT_EQ(alea_nuc_xsdir_find_temperature_bracket(
                  &xsdir, "1001.99c", 3.0e-8,
                  &lower, &upper, &fraction), ALEA_OK);
    ASSERT_TRUE(lower == &entries[0]);
    ASSERT_TRUE(upper == &entries[1]);
    ASSERT_NEAR(fraction, 0.5, 1e-15);

    ASSERT_EQ(alea_nuc_xsdir_find_temperature_bracket(
                  &xsdir, "1001.99c", 4.0e-8,
                  &lower, &upper, &fraction), ALEA_OK);
    ASSERT_TRUE(lower == &entries[1]);
    ASSERT_TRUE(upper == &entries[1]);
    ASSERT_EQ(fraction, 0.0);
    ASSERT_EQ(alea_nuc_xsdir_find_temperature_bracket(
                  &xsdir, "1001.99c", 1.0e-8,
                  &lower, &upper, &fraction), ALEA_ERR_NOT_FOUND);
    ASSERT_NULL(lower);
    ASSERT_NULL(upper);

    snprintf(entries[3].zaid, sizeof(entries[3].zaid), "1001.83c");
    entries[3].temperature = 4.0e-8;
    ASSERT_EQ(alea_nuc_xsdir_find_temperature_bracket(
                  &xsdir, "1001.99c", 4.0e-8,
                  &lower, &upper, &fraction), ALEA_ERR_INVALID_STATE);
    ASSERT_NULL(lower);
    ASSERT_NULL(upper);
}

TEST(material_temperature_mix_interpolates_and_selects_tables) {
    double energy[] = {1.0, 2.0};
    double lower_xs[] = {2.0, 2.0};
    double upper_xs[] = {6.0, 6.0};
    alea_nuc_nuclide_t lower = {0}, upper = {0};
    lower.Z = upper.Z = 1;
    lower.A = upper.A = 1;
    lower.particle = upper.particle = ALEA_NUC_PARTICLE_NEUTRON;
    lower.temperature = 2.0e-8;
    upper.temperature = 6.0e-8;
    lower.n_energies = upper.n_energies = 2;
    lower.energy = upper.energy = energy;
    lower.sigma_total = lower_xs;
    upper.sigma_total = upper_xs;

    alea_nuc_material_t* material = alea_nuc_material_create();
    ASSERT_NOT_NULL(material);
    ASSERT_EQ(alea_nuc_material_add_temperature_mix(
                  material, &lower, &upper, 0.25, 0.4), ALEA_OK);
    ASSERT_EQ(material->n_components, 2);
    ASSERT_NEAR(material->components[0].number_density, 0.3, 1e-15);
    ASSERT_NEAR(material->components[1].number_density, 0.1, 1e-15);
    ASSERT_NEAR(alea_nuc_mat_xs_total(material, 1.5), 1.2, 1e-14);

    nucdata_rng_t rng = {UINT64_C(0x12fb791a33d96e41)};
    int upper_selected = 0;
    const int samples = 50000;
    for (int i = 0; i < samples; i++) {
        alea_nuc_nuclide_t* selected = NULL;
        int component = alea_nuc_sample_nuclide(
            material, 1.5, nucdata_rng(&rng), &selected);
        ASSERT_TRUE(component == 0 || component == 1);
        ASSERT_TRUE(selected == (component == 0 ? &lower : &upper));
        if (component == 1) upper_selected++;
    }
    ASSERT_NEAR((double)upper_selected / samples, 0.5, 0.01);
    alea_nuc_material_destroy(material);

    material = alea_nuc_material_create();
    ASSERT_NOT_NULL(material);
    ASSERT_EQ(alea_nuc_material_add_temperature_mix(
                  material, &lower, &lower, 0.0, 0.4), ALEA_OK);
    ASSERT_EQ(material->n_components, 1);
    ASSERT_TRUE(material->components[0].nuclide == &lower);
    ASSERT_EQ(alea_nuc_material_add_temperature_mix(
                  material, &lower, &upper, NAN, 0.4), ALEA_ERR_INVALID_ARG);
    ASSERT_EQ(material->n_components, 1);
    alea_nuc_material_destroy(material);
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

TEST(ace_decode_cleans_every_injected_allocation_failure) {
    const char* path = "nucdata_allocation_fixture.tmp";
    ASSERT_TRUE(write_minimal_ascii_ace(path, 1));
    alea_nuc_xsdir_entry_t entry;
    alea_nuc_xsdir_t xsdir;
    memset(&entry, 0, sizeof(entry));
    snprintf(entry.zaid, sizeof(entry.zaid), "1001.80c");
    snprintf(entry.filename, sizeof(entry.filename), "%s", path);
    entry.file_type = 1;
    entry.address = 1;
    memset(&xsdir, 0, sizeof(xsdir));
    xsdir.count = 1;
    xsdir.entries = &entry;

    size_t successful_at = SIZE_MAX;
    for (size_t fail_at = 0; fail_at < 128; fail_at++) {
        nuc_alloc_failure_t failure = {fail_at, 0};
        alea_nuc_set_alloc_failure(fail_nuc_allocation, &failure);
        alea_nuc_nuclide_t* nuc =
            alea_nuc_load_nuclide(&xsdir, "1001.80c");
        alea_nuc_set_alloc_failure(NULL, NULL);
        if (nuc) {
            successful_at = fail_at;
            alea_nuc_nuclide_free(nuc);
            break;
        }
        ASSERT_TRUE(failure.calls > fail_at);
    }
    ASSERT_NE(successful_at, SIZE_MAX);
    ASSERT_TRUE(successful_at >= 8);

    alea_nuc_nuclide_t* clean =
        alea_nuc_load_nuclide(&xsdir, "1001.80c");
    ASSERT_NOT_NULL(clean);
    alea_nuc_nuclide_free(clean);
    remove(path);
}

TEST(photoatomic_decode_cleans_every_injected_allocation_failure) {
    const char* path = "nucdata_photon_allocation_fixture.tmp";
    ASSERT_TRUE(write_minimal_photoatomic_ace(path));
    alea_nuc_xsdir_entry_t entry;
    alea_nuc_xsdir_t xsdir;
    memset(&entry, 0, sizeof(entry));
    snprintf(entry.zaid, sizeof(entry.zaid), "92000.31p");
    snprintf(entry.filename, sizeof(entry.filename), "%s", path);
    entry.file_type = 1;
    entry.address = 1;
    memset(&xsdir, 0, sizeof(xsdir));
    xsdir.count = 1;
    xsdir.entries = &entry;

    size_t successful_at = SIZE_MAX;
    for (size_t fail_at = 0; fail_at < 128; fail_at++) {
        nuc_alloc_failure_t failure = {fail_at, 0};
        alea_nuc_set_alloc_failure(fail_nuc_allocation, &failure);
        alea_nuc_nuclide_t* photon =
            alea_nuc_load_nuclide(&xsdir, "92000.31p");
        alea_nuc_set_alloc_failure(NULL, NULL);
        if (photon) {
            successful_at = fail_at;
            ASSERT_EQ(photon->photon->n_fluorescence, 6);
            alea_nuc_nuclide_free(photon);
            break;
        }
        ASSERT_TRUE(failure.calls > fail_at);
    }
    ASSERT_NE(successful_at, SIZE_MAX);
    ASSERT_TRUE(successful_at >= 20);
    remove(path);
}

TEST(epr_decode_cleans_every_nested_allocation_failure) {
    const char* path = "nucdata_epr_allocation_fixture.tmp";
    ASSERT_TRUE(write_minimal_epr_ace(path));
    alea_nuc_xsdir_entry_t entry;
    alea_nuc_xsdir_t xsdir;
    memset(&entry, 0, sizeof(entry));
    snprintf(entry.zaid, sizeof(entry.zaid), "82000.14p");
    snprintf(entry.filename, sizeof(entry.filename), "%s", path);
    entry.file_type = 1;
    entry.address = 1;
    memset(&xsdir, 0, sizeof(xsdir));
    xsdir.count = 1;
    xsdir.entries = &entry;

    size_t successful_at = SIZE_MAX;
    for (size_t fail_at = 0; fail_at < 128; fail_at++) {
        nuc_alloc_failure_t failure = {fail_at, 0};
        alea_nuc_set_alloc_failure(fail_nuc_allocation, &failure);
        alea_nuc_nuclide_t* photon =
            alea_nuc_load_nuclide(&xsdir, "82000.14p");
        alea_nuc_set_alloc_failure(NULL, NULL);
        if (photon) {
            successful_at = fail_at;
            ASSERT_EQ(photon->photon->n_subshells, 3);
            ASSERT_EQ(photon->photon->n_compton_profiles, 1);
            ASSERT_EQ(photon->photon->compton_profiles[0].n_momenta, 2);
            ASSERT_EQ(photon->photon->subshells[0].max_relaxation_photons, 2);
            alea_nuc_nuclide_free(photon);
            break;
        }
        ASSERT_TRUE(failure.calls > fail_at);
    }
    ASSERT_NE(successful_at, SIZE_MAX);
    ASSERT_TRUE(successful_at >= 25);
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

static void synthetic_photon(alea_nuc_nuclide_t* nuc,
                             alea_nuc_photon_data_t* photon,
                             double* ln_energy, double* coherent,
                             double* incoherent, double* photoelectric,
                             double* pair) {
    memset(photon, 0, sizeof(*photon));
    photon->n_energies = 2; photon->ln_energy = ln_energy;
    photon->ln_sigma_coherent = coherent;
    photon->ln_sigma_incoherent = incoherent;
    photon->ln_sigma_photoelectric = photoelectric;
    photon->ln_sigma_pair = pair;
    memset(nuc, 0, sizeof(*nuc));
    nuc->particle = ALEA_NUC_PARTICLE_PHOTON; nuc->photon = photon;
}

TEST(photon_collision_samples_coherent_angle_and_preserves_energy) {
    double ln_energy[] = {log(1.0e-4), log(10.0)};
    double active[] = {0.0, 0.0}, absent[] = {-HUGE_VAL, -HUGE_VAL};
    double momentum[] = {0.0, 1.0}, integral[] = {0.0, 1.0};
    alea_nuc_photon_data_t photon;
    alea_nuc_nuclide_t nuc;
    synthetic_photon(&nuc, &photon, ln_energy, active, absent, absent, absent);
    photon.n_coherent_ff = 2; photon.coherent_momentum = momentum;
    photon.coherent_ff_cumulative = integral;
    alea_nuc_particle_state_t incident = {
        ALEA_NUC_PARTICLE_PHOTON, 1.0 / 80.65543896,
        {0.0, 0.0, 1.0}, 1.0, 2.0
    };
    double draws[] = {0.5, 0.25, 0.0, 0.5};
    sequence_rng_t rng = {draws, 4, 0};
    alea_nuc_collision_result_t result;
    ASSERT_EQ(alea_nuc_sample_photon_collision(
                  &nuc, &incident, sequence_rng, &rng, &result), ALEA_OK);
    ASSERT_EQ(result.mt, 502);
    ASSERT_EQ(result.outcome, ALEA_NUC_OUTCOME_SCATTERED);
    ASSERT_NEAR(result.mu_lab, 0.5, 1e-10);
    ASSERT_NEAR(result.outgoing.energy, incident.energy, 1e-15);
    ASSERT_NEAR(result.local_energy_deposition, 0.0, 1e-15);
}

TEST(photon_collision_samples_form_factor_corrected_compton_event) {
    double ln_energy[] = {log(1.0e-4), log(10.0)};
    double active[] = {0.0, 0.0}, absent[] = {-HUGE_VAL, -HUGE_VAL};
    double momentum[] = {0.0, 1.0}, factor[] = {0.0, 1.0};
    alea_nuc_photon_data_t photon;
    alea_nuc_nuclide_t nuc;
    synthetic_photon(&nuc, &photon, ln_energy, absent, active, absent, absent);
    photon.n_incoherent_ff = 2; photon.incoherent_momentum = momentum;
    photon.incoherent_ff = factor;
    alea_nuc_particle_state_t incident = {
        ALEA_NUC_PARTICLE_PHOTON, 0.01, {0.0, 0.0, 1.0}, 1.0, 0.0
    };
    double draws[] = {0.5, 0.5, 0.0, 0.0, 0.5};
    sequence_rng_t rng = {draws, 5, 0};
    alea_nuc_collision_result_t result;
    ASSERT_EQ(alea_nuc_sample_photon_collision(
                  &nuc, &incident, sequence_rng, &rng, &result), ALEA_OK);
    ASSERT_EQ(result.mt, 504);
    ASSERT_EQ(result.outcome, ALEA_NUC_OUTCOME_SCATTERED);
    ASSERT_NEAR(result.mu_lab, 0.0, 1e-12);
    ASSERT_NEAR(result.outgoing.energy,
                0.01 / (1.0 + 0.01 / 0.51099895069), 1e-14);
    ASSERT_NEAR(result.local_energy_deposition,
                incident.energy - result.outgoing.energy, 1e-14);
}

TEST(bound_compton_broadens_energy_and_emits_relaxation_photon) {
    double ln_energy[] = {log(1.0e-4), log(10.0)};
    double active[] = {0.0, 0.0}, absent[] = {-HUGE_VAL, -HUGE_VAL};
    double factor_momentum[] = {0.0, 100.0}, factor[] = {1.0, 1.0};
    double profile_momentum[] = {0.0, 1.0};
    double profile_pdf[] = {1.0, 1.0}, profile_cdf[] = {0.0, 1.0};
    alea_nuc_compton_shell_t compton_shell = {2.0, 0.005, 1.0, 0};
    alea_nuc_compton_profile_t profile = {
        2, 2, profile_momentum, profile_pdf, profile_cdf
    };
    alea_nuc_atomic_transition_t transition = {2, 0, 0.003, 1.0};
    alea_nuc_atomic_subshell_t shells[] = {
        {1, 1.0, 0.005, 1.0, 1, &transition, absent, 1},
        {2, 1.0, 0.001, 1.0, 0, NULL, absent, 0}
    };
    alea_nuc_photon_data_t photon;
    alea_nuc_nuclide_t nuc;
    synthetic_photon(&nuc, &photon, ln_energy,
                     absent, active, absent, absent);
    photon.n_incoherent_ff = 2;
    photon.incoherent_momentum = factor_momentum;
    photon.incoherent_ff = factor;
    photon.epr_format = 3;
    photon.n_compton_shells = 1;
    photon.compton_shells = &compton_shell;
    photon.n_compton_profiles = 1;
    photon.compton_profiles = &profile;
    photon.n_subshells = 2;
    photon.subshells = shells;
    alea_nuc_particle_state_t incident = {
        ALEA_NUC_PARTICLE_PHOTON, 0.1, {0.0, 0.0, 1.0}, 0.75, 3.0
    };
    alea_nuc_particle_state_t particle;
    alea_nuc_secondary_buffer_t buffer = {&particle, 1, 0};
    double draws[] = {
        0.5, 0.5, 0.0, 0.0, 0.25, 0.5,
        0.25, 0.0, 0.5, 0.0, 0.5, 0.25
    };
    sequence_rng_t rng = {draws, 12, 0};
    alea_nuc_collision_result_t result;
    ASSERT_EQ(alea_nuc_photon_secondary_capacity(&nuc, incident.energy), 1);
    ASSERT_EQ(alea_nuc_sample_photon_collision_with_secondaries(
                  &nuc, &incident, sequence_rng, &rng, &buffer, &result),
              ALEA_OK);
    ASSERT_EQ(result.mt, 504);
    ASSERT_EQ(result.outcome, ALEA_NUC_OUTCOME_SCATTERED);
    ASSERT_NEAR(result.mu_lab, 0.0, 1e-14);
    ASSERT_NEAR(result.outgoing.energy, 0.08323633170254346, 2e-14);
    ASSERT_TRUE(fabs(result.outgoing.energy -
                     0.1 / (1.0 + 0.1 / 0.51099895069)) > 1e-4);
    ASSERT_EQ(result.n_emitted, 1);
    ASSERT_EQ(buffer.count, 1);
    ASSERT_NEAR(particle.energy, 0.003, 1e-15);
    ASSERT_NEAR(result.outgoing.energy + particle.energy +
                result.local_energy_deposition, incident.energy, 2e-14);

    buffer.capacity = 0;
    buffer.count = 0;
    rng.position = 0;
    memset(&result, 0x5a, sizeof(result));
    alea_nuc_collision_result_t before = result;
    ASSERT_EQ(alea_nuc_sample_photon_collision_with_secondaries(
                  &nuc, &incident, sequence_rng, &rng, &buffer, &result),
              ALEA_ERR_OUT_OF_MEMORY);
    ASSERT_EQ(rng.position, 0);
    ASSERT_EQ(buffer.count, 0);
    ASSERT_EQ(memcmp(&result, &before, sizeof(result)), 0);
}

TEST(photon_absorption_uses_named_local_deposition_model) {
    double ln_energy[] = {log(1.0e-4), log(10.0)};
    double active[] = {0.0, 0.0}, absent[] = {-HUGE_VAL, -HUGE_VAL};
    alea_nuc_photon_data_t photon;
    alea_nuc_nuclide_t nuc;
    alea_nuc_particle_state_t incident = {
        ALEA_NUC_PARTICLE_PHOTON, 2.0, {0.0, 0.0, 1.0}, 1.0, 0.0
    };
    double draw = 0.5;
    sequence_rng_t rng = {&draw, 1, 0};
    alea_nuc_collision_result_t result;

    synthetic_photon(&nuc, &photon, ln_energy, absent, absent, active, absent);
    ASSERT_EQ(alea_nuc_sample_photon_collision(
                  &nuc, &incident, sequence_rng, &rng, &result), ALEA_OK);
    ASSERT_EQ(result.mt, 522);
    ASSERT_EQ(result.outcome, ALEA_NUC_OUTCOME_ABSORBED);
    ASSERT_NEAR(result.local_energy_deposition, 2.0, 1e-15);

    rng.position = 0;
    synthetic_photon(&nuc, &photon, ln_energy, absent, absent, absent, active);
    ASSERT_EQ(alea_nuc_sample_photon_collision(
                  &nuc, &incident, sequence_rng, &rng, &result), ALEA_OK);
    ASSERT_EQ(result.mt, 517);
    ASSERT_EQ(result.outcome, ALEA_NUC_OUTCOME_ABSORBED);
    ASSERT_NEAR(result.local_energy_deposition, 2.0, 1e-15);
}

TEST(photon_pair_production_emits_back_to_back_annihilation_photons) {
    double ln_energy[] = {log(1.1), log(10.0)};
    double active[] = {0.0, 0.0}, absent[] = {-HUGE_VAL, -HUGE_VAL};
    alea_nuc_photon_data_t photon;
    alea_nuc_nuclide_t nuc;
    synthetic_photon(&nuc, &photon, ln_energy,
                     absent, absent, absent, active);
    alea_nuc_particle_state_t incident = {
        ALEA_NUC_PARTICLE_PHOTON, 2.0, {0.0, 0.0, 1.0}, 0.75, 3.0
    };
    alea_nuc_particle_state_t particles[2];
    alea_nuc_secondary_buffer_t buffer = {particles, 2, 0};
    double draws[] = {0.5, 0.25, 0.125};
    sequence_rng_t rng = {draws, 3, 0};
    alea_nuc_collision_result_t result;
    ASSERT_EQ(alea_nuc_photon_secondary_capacity(&nuc, incident.energy), 2);
    ASSERT_EQ(alea_nuc_sample_photon_collision_with_secondaries(
                  &nuc, &incident, sequence_rng, &rng, &buffer, &result),
              ALEA_OK);
    ASSERT_EQ(result.mt, 517);
    ASSERT_EQ(result.outcome, ALEA_NUC_OUTCOME_REPLACED);
    ASSERT_EQ(result.n_emitted, 2);
    ASSERT_EQ(buffer.count, 2);
    for (int i = 0; i < 2; i++) {
        ASSERT_EQ(particles[i].type, ALEA_NUC_PARTICLE_PHOTON);
        ASSERT_NEAR(particles[i].energy, 0.51099895069, 1e-15);
        ASSERT_EQ(particles[i].weight, incident.weight);
        ASSERT_EQ(particles[i].time, incident.time);
    }
    for (int i = 0; i < 3; i++)
        ASSERT_NEAR(particles[0].direction[i],
                    -particles[1].direction[i], 1e-15);
    ASSERT_NEAR(result.local_energy_deposition,
                2.0 - 2.0 * 0.51099895069, 1e-15);

    buffer.capacity = 1;
    buffer.count = 0;
    rng.position = 0;
    memset(&result, 0x5a, sizeof(result));
    alea_nuc_collision_result_t before = result;
    ASSERT_EQ(alea_nuc_sample_photon_collision_with_secondaries(
                  &nuc, &incident, sequence_rng, &rng, &buffer, &result),
              ALEA_ERR_OUT_OF_MEMORY);
    ASSERT_EQ(rng.position, 0);
    ASSERT_EQ(buffer.count, 0);
    ASSERT_EQ(memcmp(&result, &before, sizeof(result)), 0);
}

TEST(photoelectric_relaxation_emits_epr_cascade_and_conserves_energy) {
    double ln_energy[] = {log(0.1), log(2.0)};
    double active[] = {0.0, 0.0}, absent[] = {-HUGE_VAL, -HUGE_VAL};
    alea_nuc_atomic_transition_t inner_transitions[] = {
        {2, 0, 0.080, 0.5}, {2, 3, 0.070, 1.0}
    };
    alea_nuc_atomic_transition_t middle_transitions[] = {
        {3, 0, 0.015, 1.0}
    };
    alea_nuc_atomic_subshell_t shells[] = {
        {1, 2.0, 0.100, 0.2, 2, inner_transitions, active, 2},
        {2, 2.0, 0.020, 0.7, 1, middle_transitions, absent, 1},
        {3, 4.0, 0.005, 1.0, 0, NULL, absent, 0}
    };
    alea_nuc_photon_data_t photon;
    alea_nuc_nuclide_t nuc;
    synthetic_photon(&nuc, &photon, ln_energy,
                     absent, absent, active, absent);
    photon.epr_format = 3;
    photon.n_subshells = 3;
    photon.subshells = shells;
    alea_nuc_particle_state_t incident = {
        ALEA_NUC_PARTICLE_PHOTON, 1.0, {0.0, 0.0, 1.0}, 0.75, 3.0
    };
    alea_nuc_particle_state_t particles[2];
    alea_nuc_secondary_buffer_t buffer = {particles, 2, 0};
    double draws[] = {0.5, 0.25, 0.25, 0.25, 0.0, 0.5, 0.75, 0.25};
    sequence_rng_t rng = {draws, 8, 0};
    alea_nuc_collision_result_t result;
    ASSERT_EQ(alea_nuc_photon_secondary_capacity(&nuc, incident.energy), 2);
    ASSERT_EQ(alea_nuc_photon_xs_photoelectric_subshell(&nuc, 1, 1.0), 1.0);
    ASSERT_EQ(alea_nuc_photon_xs_photoelectric_subshell(&nuc, 2, 1.0), 0.0);
    ASSERT_EQ(alea_nuc_sample_photon_collision_with_secondaries(
                  &nuc, &incident, sequence_rng, &rng, &buffer, &result),
              ALEA_OK);
    ASSERT_EQ(result.mt, 522);
    ASSERT_EQ(result.outcome, ALEA_NUC_OUTCOME_REPLACED);
    ASSERT_EQ(result.n_emitted, 2);
    ASSERT_EQ(buffer.count, 2);
    ASSERT_NEAR(particles[0].energy, 0.080, 1e-15);
    ASSERT_NEAR(particles[1].energy, 0.015, 1e-15);
    ASSERT_EQ(particles[0].weight, incident.weight);
    ASSERT_EQ(particles[1].time, incident.time);
    ASSERT_NEAR(result.local_energy_deposition, 0.905, 1e-15);
    for (int i = 0; i < 2; i++) {
        double norm2 = 0.0;
        for (int axis = 0; axis < 3; axis++)
            norm2 += particles[i].direction[axis] * particles[i].direction[axis];
        ASSERT_NEAR(norm2, 1.0, 1e-14);
    }

    buffer.capacity = 1;
    buffer.count = 0;
    rng.position = 0;
    memset(&result, 0x5a, sizeof(result));
    alea_nuc_collision_result_t before = result;
    ASSERT_EQ(alea_nuc_sample_photon_collision_with_secondaries(
                  &nuc, &incident, sequence_rng, &rng, &buffer, &result),
              ALEA_ERR_OUT_OF_MEMORY);
    ASSERT_EQ(rng.position, 0);
    ASSERT_EQ(buffer.count, 0);
    ASSERT_EQ(memcmp(&result, &before, sizeof(result)), 0);
}

TEST(prepared_photon_collision_banks_annihilation_photons) {
    double energy[] = {1.1, 10.0};
    double ln_energy[] = {log(1.1), log(10.0)};
    double active[] = {0.0, 0.0}, absent[] = {-HUGE_VAL, -HUGE_VAL};
    alea_nuc_photon_data_t photon;
    alea_nuc_nuclide_t element;
    synthetic_photon(&element, &photon, ln_energy,
                     absent, absent, absent, active);
    photon.energy = energy;
    element.n_energies = 2;
    element.energy = energy;
    alea_nuc_mat_component_t component = {&element, 0.5};
    alea_nuc_material_t material = {&component, 1, 1};
    alea_nuc_prepare_requirements_t requirements = {
        .required_capabilities = ALEA_NUC_CAP_PHOTON
    };
    alea_nuc_capability_report_t report;
    alea_nuc_prepared_material_t* prepared = NULL;
    ASSERT_EQ(alea_nuc_prepare_material(&material, &requirements, &report,
                                        &prepared), ALEA_OK);
    alea_nuc_particle_state_t incident = {
        ALEA_NUC_PARTICLE_PHOTON, 2.0, {0.0, 0.0, 1.0}, 1.0, 0.0
    };
    alea_nuc_evaluation_t evaluation;
    ASSERT_EQ(alea_nuc_evaluate(prepared, &incident, &evaluation), ALEA_OK);
    alea_nuc_particle_state_t particles[2];
    alea_nuc_secondary_buffer_t buffer = {particles, 2, 0};
    double draws[] = {0.5, 0.5, 0.25, 0.125};
    sequence_rng_t rng = {draws, 4, 0};
    alea_nuc_collision_result_t result;
    ASSERT_EQ(alea_nuc_collide_with_secondaries(
                  &evaluation, sequence_rng, &rng, &buffer, &result), ALEA_OK);
    ASSERT_EQ(result.component_index, 0);
    ASSERT_EQ(result.mt, 517);
    ASSERT_EQ(result.n_emitted, 2);
    ASSERT_EQ(buffer.count, 2);
    alea_nuc_prepared_material_free(prepared);
}

TEST(prepared_photon_mixture_evaluates_and_selects_element) {
    double energy[] = {0.1, 10.0};
    double ln_energy[] = {log(0.1), log(10.0)};
    double one[] = {0.0, 0.0}, three[] = {log(3.0), log(3.0)};
    double absent[] = {-HUGE_VAL, -HUGE_VAL};
    alea_nuc_photon_data_t photons[2];
    alea_nuc_nuclide_t elements[2];
    synthetic_photon(&elements[0], &photons[0], ln_energy,
                     absent, absent, one, absent);
    synthetic_photon(&elements[1], &photons[1], ln_energy,
                     absent, absent, three, absent);
    for (int i = 0; i < 2; i++) {
        photons[i].energy = energy;
        elements[i].n_energies = 2; elements[i].energy = energy;
    }
    alea_nuc_mat_component_t components[] = {
        {&elements[0], 1.0}, {&elements[1], 1.0}
    };
    alea_nuc_material_t material = {components, 2, 2};
    alea_nuc_prepare_requirements_t requirements = {
        .required_capabilities = ALEA_NUC_CAP_PHOTON
    };
    alea_nuc_capability_report_t report;
    alea_nuc_prepared_material_t* prepared = NULL;
    ASSERT_EQ(alea_nuc_prepare_material(&material, &requirements, &report,
                                        &prepared), ALEA_OK);
    alea_nuc_particle_state_t incident = {
        ALEA_NUC_PARTICLE_PHOTON, 1.0, {0.0, 0.0, 1.0}, 1.0, 0.0
    };
    alea_nuc_evaluation_t evaluation;
    ASSERT_EQ(alea_nuc_evaluate(prepared, &incident, &evaluation), ALEA_OK);
    ASSERT_NEAR(evaluation.macro_total, 4.0, 1e-12);
    ASSERT_NEAR(evaluation.macro_elastic, 0.0, 1e-12);
    ASSERT_NEAR(evaluation.macro_absorption, 4.0, 1e-12);
    double draws[] = {0.9, 0.5};
    sequence_rng_t rng = {draws, 2, 0};
    alea_nuc_collision_result_t result;
    ASSERT_EQ(alea_nuc_collide(&evaluation, sequence_rng, &rng, &result),
              ALEA_OK);
    ASSERT_EQ(result.component_index, 1);
    ASSERT_EQ(result.mt, 522);
    alea_nuc_prepared_material_free(prepared);
}

TEST(preparation_rejects_mixed_neutron_and_photon_tables) {
    double energy[] = {0.1, 10.0}, total[] = {1.0, 1.0};
    double ln_energy[] = {log(0.1), log(10.0)};
    double active[] = {0.0, 0.0}, absent[] = {-HUGE_VAL, -HUGE_VAL};
    alea_nuc_photon_data_t photon;
    alea_nuc_nuclide_t tables[2];
    memset(&tables[0], 0, sizeof(tables[0]));
    tables[0].particle = ALEA_NUC_PARTICLE_NEUTRON; tables[0].awr = 1.0;
    tables[0].n_energies = 2; tables[0].energy = energy;
    tables[0].sigma_total = total; tables[0].sigma_elastic = total;
    tables[0].sigma_abs = (double[2]){0.0, 0.0};
    synthetic_photon(&tables[1], &photon, ln_energy,
                     absent, absent, active, absent);
    photon.energy = energy; tables[1].n_energies = 2; tables[1].energy = energy;
    alea_nuc_mat_component_t components[] = {
        {&tables[0], 1.0}, {&tables[1], 1.0}
    };
    alea_nuc_material_t material = {components, 2, 2};
    alea_nuc_prepare_requirements_t requirements = {
        .required_capabilities = ALEA_NUC_CAP_RESTRICTED_NEUTRON
    };
    alea_nuc_capability_report_t report;
    alea_nuc_prepared_material_t* prepared = NULL;
    ASSERT_EQ(alea_nuc_prepare_material(&material, &requirements, &report,
                                        &prepared), ALEA_ERR_UNSUPPORTED);
    ASSERT_NULL(prepared);
    ASSERT_EQ(report.issue, ALEA_NUC_PREP_UNSUPPORTED_PARTICLE);
    ASSERT_EQ(report.component_index, 1);
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

TEST(urr_clamps_negative_collision_cross_sections_but_keeps_signed_heating) {
    double energy[] = {1.0, 3.0};
    double table[] = {
        1, 10, 4, 0, -0.25, -2,
        1, 10, 4, 0, -0.50, -4
    };
    double factors[5];
    alea_nuc_urr_t urr = {
        .n_energies = 2, .n_bands = 1, .interp = 2,
        .multiply_smooth = true, .energy = energy, .table = table
    };
    alea_nuc_nuclide_t nuc = {0};
    nuc.urr = &urr;

    ASSERT_EQ(alea_nuc_urr_factors(&nuc, 2.0, 0.5, factors), 1);
    ASSERT_NEAR(factors[0], 10.0, 1e-12);
    ASSERT_NEAR(factors[1], 4.0, 1e-12);
    ASSERT_NEAR(factors[2], 0.0, 1e-12);
    ASSERT_NEAR(factors[3], 0.0, 1e-12);
    ASSERT_NEAR(factors[4], -3.0, 1e-12);
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

TEST(energy_decode_accepts_empty_block_at_end_of_xss) {
    double raw[] = {0.0};
    alea_nuc_reaction_t reaction = {.mt = 2, .ty = 1};
    alea_nuc_nuclide_t nuc = {0};
    nuc.n_reactions = 1;
    nuc.reactions = &reaction;
    nuc.raw.xss = raw;
    nuc.raw.xss_length = 1;
    nuc.raw.jxs[9] = 2;
    nuc.raw.jxs[10] = 2;

    alea_nuc_decode_all_energy(&nuc);
    ASSERT_NULL(reaction.energy);
    ASSERT_FALSE(nuc.raw.decode_error);
}

TEST(watt_decode_keeps_independent_a_and_b_grids) {
    double raw[] = {
        0, 11, 10, 0, 2, 1, 3, 1, 1,
        1, 2, 5, 2, 1, 3, 1, 1,
        1, 1, 2, 1, 1, 2, 0
    };
    alea_nuc_ace_table_t table;
    memset(&table, 0, sizeof(table));
    table.xss = raw; table.xss_length = (int)(sizeof(raw) / sizeof(raw[0]));
    table.jxs[10] = 1;

    alea_nuc_energy_dist_t* ed = alea_nuc_decode_energy_dist(&table, 1);
    ASSERT_NOT_NULL(ed);
    ASSERT_EQ(ed->law, ALEA_NUC_ELAW_WATT);
    ASSERT_EQ(ed->n_temp, 2);
    ASSERT_EQ(ed->n_temp_regions, 1);
    ASSERT_EQ(ed->temp_nbt[0], 2);
    ASSERT_EQ(ed->temp_interp[0], 5);
    ASSERT_EQ(ed->n_watt_b, 1);
    ASSERT_EQ(ed->n_watt_b_regions, 1);
    ASSERT_EQ(ed->watt_b_nbt[0], 1);
    ASSERT_EQ(ed->watt_b_interp[0], 2);
    ASSERT_NEAR(ed->temp_energy[1], 3.0, 1e-12);
    ASSERT_NEAR(ed->watt_b_energy[0], 1.0, 1e-12);
    ASSERT_NEAR(ed->temp_C[0], 2.0, 1e-12);

    free(ed->nbt); free(ed->interp); free(ed->energy); free(ed->probability);
    free(ed->temp_energy); free(ed->temp_T);
    free(ed->temp_nbt); free(ed->temp_interp);
    free(ed->watt_b_energy); free(ed->temp_C);
    free(ed->watt_b_nbt); free(ed->watt_b_interp); free(ed);
}

TEST(general_evaporation_decode_preserves_temperature_and_bins) {
    double raw[] = {
        0, 5, 10, 0, 2, 1, 3, 1, 1,
        1, 2, 2, 2, 1, 3, 2, 4, 3, 0, 1, 3
    };
    alea_nuc_ace_table_t table;
    memset(&table, 0, sizeof(table));
    table.xss = raw; table.xss_length = (int)(sizeof(raw) / sizeof(raw[0]));
    table.jxs[10] = 1;
    alea_nuc_energy_dist_t* ed = alea_nuc_decode_energy_dist(&table, 1);
    ASSERT_NOT_NULL(ed);
    ASSERT_EQ(ed->law, ALEA_NUC_ELAW_GENERAL_EVAP);
    ASSERT_EQ(ed->n_temp_regions, 1);
    ASSERT_EQ(ed->temp_nbt[0], 2);
    ASSERT_EQ(ed->temp_interp[0], 2);
    ASSERT_EQ(ed->n_temp, 2);
    ASSERT_NEAR(ed->temp_energy[1], 3.0, 1e-12);
    ASSERT_NEAR(ed->temp_T[1], 4.0, 1e-12);
    ASSERT_EQ(ed->n_general_evap, 3);
    ASSERT_NEAR(ed->general_evap_x[0], 0.0, 1e-12);
    ASSERT_NEAR(ed->general_evap_x[2], 3.0, 1e-12);
    ASSERT_EQ(alea_nuc_energy_dist_validate(ed, NULL), ALEA_OK);
    alea_nuc_energy_dist_free(ed);
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

TEST(doppler_preserves_one_over_v_cross_section) {
    enum { N = 401 };
    double energy[N], total[N];
    for (int i = 0; i < N; i++) {
        energy[i] = pow(10.0, -10.0 + 8.0 * i / (N - 1));
        total[i] = 1.0 / sqrt(energy[i]);
    }
    alea_nuc_nuclide_t nuc = {0};
    nuc.awr = 238.0;
    nuc.temperature = 0.0;
    nuc.n_energies = N;
    nuc.energy = energy;
    nuc.sigma_total = total;
    ASSERT_EQ(alea_nuc_doppler_broaden(&nuc, 2.53e-8), ALEA_OK);
    for (int i = 50; i < N - 50; i++)
        ASSERT_NEAR(total[i] * sqrt(energy[i]), 1.0, 2.0e-4);
}

TEST(doppler_narrow_resonance_matches_independent_quadrature) {
    enum { N = 403 };
    double energy[N], total[N];
    energy[0] = 1.0e-8;
    for (int i = 1; i < N - 1; i++)
        energy[i] = 9.8e-6 + (i - 1) * 1.0e-9;
    energy[N - 1] = 1.0e-3;
    for (int i = 0; i < N; i++) {
        double offset = (energy[i] - 1.0e-5) / 2.0e-8;
        total[i] = 1.0 + 1000.0 / (1.0 + offset * offset);
    }
    alea_nuc_nuclide_t nuc = {0};
    nuc.awr = 238.0;
    nuc.temperature = 2.53e-8;
    nuc.n_energies = N;
    nuc.energy = energy;
    nuc.sigma_total = total;
    ASSERT_EQ(alea_nuc_doppler_broaden(&nuc, 5.06e-8), ALEA_OK);

    const int indices[] = {101, 151, 201, 251, 301};
    const double reference[] = {
        91.08437061341653, 265.3889497482132, 397.8544755861514,
        263.7873032117916, 90.88081703210231
    };
    for (int i = 0; i < 5; i++)
        ASSERT_NEAR(total[indices[i]], reference[i], 0.25);
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

TEST(nu_bar_honors_loglog_interpolation_metadata) {
    int nbt[] = {2};
    int interp[] = {5};
    double energy[] = {1.0, 100.0};
    double values[] = {1.0, 100.0};
    alea_nuc_nu_bar_t total;
    alea_nuc_fission_t fission;
    alea_nuc_nuclide_t nuc;
    memset(&total, 0, sizeof(total));
    total.type = ALEA_NUC_NU_TABULAR;
    total.n_regions = 1; total.nbt = nbt; total.interp = interp;
    total.n_energies = 2; total.energy = energy; total.nu = values;
    memset(&fission, 0, sizeof(fission)); fission.total = &total;
    memset(&nuc, 0, sizeof(nuc)); nuc.fission = &fission;
    ASSERT_NEAR(alea_nuc_nu_bar(&nuc, 10.0), 10.0, 1e-12);
}

TEST(fission_collision_yield_uses_prompt_nu_bar) {
    double total_coeff[] = {2.5};
    double prompt_coeff[] = {2.0};
    alea_nuc_nu_bar_t total, prompt;
    alea_nuc_fission_t fission;
    alea_nuc_reaction_t reaction;
    alea_nuc_nuclide_t nuc;
    memset(&total, 0, sizeof(total)); memset(&prompt, 0, sizeof(prompt));
    total.type=ALEA_NUC_NU_POLYNOMIAL; total.n_coeffs=1;
    total.coeffs=total_coeff;
    prompt.type=ALEA_NUC_NU_POLYNOMIAL; prompt.n_coeffs=1;
    prompt.coeffs=prompt_coeff;
    memset(&fission, 0, sizeof(fission));
    fission.total=&total; fission.prompt=&prompt;
    memset(&reaction, 0, sizeof(reaction)); reaction.mt=18; reaction.ty=19;
    memset(&nuc, 0, sizeof(nuc)); nuc.fission=&fission;
    nuc.n_reactions=1; nuc.reactions=&reaction;
    ASSERT_NEAR(alea_nuc_nu_bar(&nuc, 1.0), 2.5, 1e-12);
    ASSERT_NEAR(alea_nuc_prompt_nu_bar(&nuc, 1.0), 2.0, 1e-12);
    ASSERT_NEAR(alea_nuc_reaction_yield(&nuc, 18, 1.0), 2.0, 1e-12);
}

TEST(shared_interpolation_honors_all_endf_schemes) {
    double x[] = {1.0, 100.0};
    double linear_y[] = {1.0, 3.0};
    double log_y[] = {1.0, 100.0};
    int nbt[] = {2};
    int interp[] = {1};
    double value = 0.0;

    ASSERT_EQ(alea_nuc_interp_eval(x, linear_y, 2, nbt, interp, 1,
                                   10.0, &value), ALEA_OK);
    ASSERT_NEAR(value, 1.0, 1e-12);
    interp[0] = 2;
    ASSERT_EQ(alea_nuc_interp_eval(x, linear_y, 2, nbt, interp, 1,
                                   50.5, &value), ALEA_OK);
    ASSERT_NEAR(value, 2.0, 1e-12);
    interp[0] = 3;
    ASSERT_EQ(alea_nuc_interp_eval(x, linear_y, 2, nbt, interp, 1,
                                   10.0, &value), ALEA_OK);
    ASSERT_NEAR(value, 2.0, 1e-12);
    interp[0] = 4;
    ASSERT_EQ(alea_nuc_interp_eval(x, log_y, 2, nbt, interp, 1,
                                   50.5, &value), ALEA_OK);
    ASSERT_NEAR(value, 10.0, 1e-12);
    interp[0] = 5;
    ASSERT_EQ(alea_nuc_interp_eval(x, log_y, 2, nbt, interp, 1,
                                   10.0, &value), ALEA_OK);
    ASSERT_NEAR(value, 10.0, 1e-12);
}

TEST(shared_interpolation_accepts_endf_discontinuities) {
    double x[] = {0.8, 1.0, 1.0, 2.0};
    double y[] = {0.1, 0.2, 0.7, 0.9};
    double value = 0.0;

    ASSERT_EQ(alea_nuc_interp_eval(x, y, 4, NULL, NULL, 0,
                                   0.9, &value), ALEA_OK);
    ASSERT_NEAR(value, 0.15, 1e-14);
    ASSERT_EQ(alea_nuc_interp_eval(x, y, 4, NULL, NULL, 0,
                                   1.0, &value), ALEA_OK);
    ASSERT_NEAR(value, 0.7, 1e-14);
    ASSERT_EQ(alea_nuc_interp_eval(x, y, 4, NULL, NULL, 0,
                                   1.5, &value), ALEA_OK);
    ASSERT_NEAR(value, 0.8, 1e-14);
}

TEST(shared_tabular_sampler_handles_discrete_continuous_mixture) {
    double x[] = {2.0, 0.0, 2.0};
    double pdf[] = {0.25, 0.375, 0.375};
    double cdf[] = {0.25, 0.25, 1.0};
    double value = -1.0;
    int sampled = -1;
    ASSERT_TRUE(alea_nuc_tabular_pdf_valid(x, pdf, cdf, 3, 1, 1));
    ASSERT_EQ(alea_nuc_tabular_pdf_sample(x, pdf, cdf, 3, 1, 1, 0.1,
                                          &value, &sampled), ALEA_OK);
    ASSERT_NEAR(value, 2.0, 1e-12);
    ASSERT_EQ(sampled, 0);
    ASSERT_EQ(alea_nuc_tabular_pdf_sample(x, pdf, cdf, 3, 1, 1, 0.625,
                                          &value, &sampled), ALEA_OK);
    ASSERT_NEAR(value, 1.0, 1e-12);
    ASSERT_EQ(sampled, 1);
}

TEST(shared_tabular_sampler_accepts_pure_discrete_interpolation_zero) {
    double x[] = {1.25};
    double pdf[] = {1.0};
    double cdf[] = {1.0};
    double value = -1.0;
    int sampled = -1;
    ASSERT_TRUE(alea_nuc_tabular_pdf_valid(x, pdf, cdf, 1, 0, 1));
    ASSERT_EQ(alea_nuc_tabular_pdf_sample(x, pdf, cdf, 1, 0, 1, 0.5,
                                          &value, &sampled), ALEA_OK);
    ASSERT_NEAR(value, 1.25, 1e-14);
    ASSERT_EQ(sampled, 0);
}

TEST(shared_tabular_sampler_accepts_zero_mass_duplicate_endpoint) {
    double x[] = {0.0, 0.2, 0.2};
    double pdf[] = {0.0, 10.0, 0.0};
    double cdf[] = {0.0, 1.0, 1.0};
    double value = -1.0;
    ASSERT_TRUE(alea_nuc_tabular_pdf_valid(x, pdf, cdf, 3, 2, 0));
    ASSERT_EQ(alea_nuc_tabular_pdf_sample(x, pdf, cdf, 3, 2, 0, 0.5,
                                          &value, NULL), ALEA_OK);
    ASSERT_TRUE(value >= 0.0 && value <= 0.2);

    cdf[1] = 0.9;
    ASSERT_FALSE(alea_nuc_tabular_pdf_valid(x, pdf, cdf, 3, 2, 0));
}

TEST(shared_tabular_sampler_accepts_rounded_and_unreachable_tail) {
    double x[] = {0.0, 1.0, 0.999};
    double pdf[] = {1.0, 0.0, 0.0};
    double cdf[] = {0.0, 1.0, 1.00000001};
    double value = -1.0;
    ASSERT_TRUE(alea_nuc_tabular_pdf_valid(x, pdf, cdf, 3, 1, 0));
    ASSERT_EQ(alea_nuc_tabular_pdf_sample(x, pdf, cdf, 3, 1, 0, 0.75,
                                          &value, NULL), ALEA_OK);
    ASSERT_NEAR(value, 0.75, 1e-12);
}

TEST(energy_law_applicability_accepts_right_continuous_boundary) {
    double applicability_energy[] = {4.0, 4.0, 20.0};
    double probability[] = {0.0, 0.5, 0.5};
    double temp_energy[] = {4.0, 20.0};
    double temperature[] = {1.0, 1.0};
    alea_nuc_energy_dist_t law = {0};
    law.law = ALEA_NUC_ELAW_EVAPORATION;
    law.n_energies = 3;
    law.energy = applicability_energy;
    law.probability = probability;
    law.n_temp = 2;
    law.temp_energy = temp_energy;
    law.temp_T = temperature;
    ASSERT_EQ(alea_nuc_energy_dist_validate(&law, NULL), ALEA_OK);
}

TEST(energy_sampler_level_uses_ace_threshold_and_mass_ratio) {
    alea_nuc_energy_dist_t law;
    memset(&law, 0, sizeof(law));
    law.law = ALEA_NUC_ELAW_LEVEL;
    law.level_A = 1.0;
    law.level_Q = 0.8;
    double unused_draw = 0.5;
    sequence_rng_t rng = {&unused_draw, 1, 0};
    double energy = -1.0;
    ASSERT_EQ(alea_nuc_sample_energy_distribution(
                  &law, 6.0, sequence_rng, &rng, &energy), ALEA_OK);
    ASSERT_NEAR(energy, 4.0, 1e-12);
    ASSERT_EQ(rng.position, 0);
}

TEST(energy_sampler_selects_applicable_law) {
    int nbt[] = {2};
    int interp[] = {2};
    double incident[] = {1.0, 3.0};
    double probability[] = {0.25, 0.25};
    alea_nuc_energy_dist_t first, second;
    memset(&first, 0, sizeof(first));
    memset(&second, 0, sizeof(second));
    first.law = ALEA_NUC_ELAW_LEVEL;
    first.n_regions = 1; first.nbt = nbt; first.interp = interp;
    first.n_energies = 2; first.energy = incident;
    first.probability = probability; first.level_Q = 0.5;
    first.next = &second;
    second.law = ALEA_NUC_ELAW_LEVEL; second.level_Q = 0.25;
    double draws[] = {0.75};
    sequence_rng_t rng = {draws, 1, 0};
    double energy = -1.0;
    ASSERT_EQ(alea_nuc_sample_energy_distribution(
                  &first, 2.0, sequence_rng, &rng, &energy), ALEA_OK);
    ASSERT_NEAR(energy, 0.5, 1e-12);
}

TEST(energy_sampler_continuous_tabular_applies_unit_base_scaling) {
    int ein_nbt[] = {2};
    int ein_interp[] = {12};
    double ein[] = {1.0, 3.0};
    int outgoing_interp[] = {1, 1};
    int n_discrete[] = {0, 0};
    int n_eout[] = {2, 2};
    double eout0[] = {0.0, 2.0}, eout1[] = {2.0, 6.0};
    double pdf0[] = {0.5, 0.5}, pdf1[] = {0.25, 0.25};
    double cdf0[] = {0.0, 1.0}, cdf1[] = {0.0, 1.0};
    double* eout[] = {eout0, eout1};
    double* pdf[] = {pdf0, pdf1};
    double* cdf[] = {cdf0, cdf1};
    alea_nuc_energy_dist_t law;
    memset(&law, 0, sizeof(law)); law.law = ALEA_NUC_ELAW_CONT_TABULAR;
    law.tab.n_ein = 2; law.tab.n_regions = 1; law.tab.nbt = ein_nbt;
    law.tab.interp = ein_interp; law.tab.ein = ein;
    law.tab.interpolation = outgoing_interp; law.tab.n_discrete = n_discrete;
    law.tab.n_eout = n_eout; law.tab.eout = eout;
    law.tab.pdf = pdf; law.tab.cdf = cdf;
    /* Select the lower table, then its midpoint. At Ein=2 the interpolated
       outgoing range is [1,4], so the unit-base midpoint is 2.5 MeV. */
    double draws[] = {0.75, 0.5};
    sequence_rng_t rng = {draws, 2, 0};
    double energy = -1.0;
    ASSERT_EQ(alea_nuc_sample_energy_distribution(
                  &law, 2.0, sequence_rng, &rng, &energy), ALEA_OK);
    ASSERT_NEAR(energy, 2.5, 1e-12);
}

TEST(energy_sampler_code22_interpolates_discrete_line_energy) {
    int ein_nbt[] = {2};
    int ein_interp[] = {22};
    double ein[] = {1.0, 3.0};
    int outgoing_interp[] = {1, 1};
    int n_discrete[] = {1, 1};
    int n_eout[] = {3, 3};
    double eout0[] = {5.0, 0.0, 2.0};
    double eout1[] = {7.0, 2.0, 6.0};
    double pdf0[] = {0.25, 0.375, 0.375};
    double pdf1[] = {0.25, 0.1875, 0.1875};
    double cdf0[] = {0.25, 0.25, 1.0};
    double cdf1[] = {0.25, 0.25, 1.0};
    double* eout[] = {eout0, eout1};
    double* pdf[] = {pdf0, pdf1};
    double* cdf[] = {cdf0, cdf1};
    alea_nuc_energy_dist_t law = {0};
    law.law = ALEA_NUC_ELAW_CONT_TABULAR;
    law.tab.n_ein = 2;
    law.tab.n_regions = 1;
    law.tab.nbt = ein_nbt;
    law.tab.interp = ein_interp;
    law.tab.ein = ein;
    law.tab.interpolation = outgoing_interp;
    law.tab.n_discrete = n_discrete;
    law.tab.n_eout = n_eout;
    law.tab.eout = eout;
    law.tab.pdf = pdf;
    law.tab.cdf = cdf;
    double draws[] = {0.75, 0.1};
    sequence_rng_t rng = {draws, 2, 0};
    double sampled = -1.0;
    ASSERT_EQ(alea_nuc_energy_dist_validate(&law, NULL), ALEA_OK);
    ASSERT_EQ(alea_nuc_sample_energy_distribution(
                  &law, 2.0, sequence_rng, &rng, &sampled), ALEA_OK);
    ASSERT_NEAR(sampled, 6.0, 1e-14);
}

TEST(general_evaporation_sampler_interpolates_temperature_and_bin) {
    int nbt[] = {2}, interp[] = {2};
    double incident[] = {1.0, 3.0}, temperature[] = {2.0, 4.0};
    double boundaries[] = {0.0, 1.0, 3.0};
    alea_nuc_energy_dist_t law;
    memset(&law, 0, sizeof(law));
    law.law = ALEA_NUC_ELAW_GENERAL_EVAP;
    law.n_temp_regions = 1; law.temp_nbt = nbt; law.temp_interp = interp;
    law.n_temp = 2; law.temp_energy = incident; law.temp_T = temperature;
    law.n_general_evap = 3; law.general_evap_x = boundaries;
    double draw = 0.75;
    sequence_rng_t rng = {&draw, 1, 0};
    double energy = -1.0;
    ASSERT_EQ(alea_nuc_energy_dist_validate(&law, NULL), ALEA_OK);
    ASSERT_EQ(alea_nuc_sample_energy_distribution(
                  &law, 2.0, sequence_rng, &rng, &energy), ALEA_OK);
    ASSERT_NEAR(energy, 6.0, 1e-12);
}

TEST(energy_sampler_maxwell_evaporation_and_watt_obey_cutoff) {
    double grid[] = {0.0, 10.0};
    double one[] = {1.0, 1.0};
    alea_nuc_energy_dist_t law;
    memset(&law, 0, sizeof(law));
    law.n_temp = 2; law.temp_energy = grid; law.temp_T = one;
    law.restriction_energy = 1.0;
    double draws[] = {
        0.5, 0.5, 0.5,
        0.5, 0.5,
        0.5, 0.5, 0.5, 0.5
    };
    sequence_rng_t rng = {draws, 9, 0};
    double energy = -1.0;

    law.law = ALEA_NUC_ELAW_MAXWELL;
    ASSERT_EQ(alea_nuc_sample_energy_distribution(
                  &law, 5.0, sequence_rng, &rng, &energy), ALEA_OK);
    ASSERT_TRUE(energy >= 0.0 && energy <= 4.0);

    law.law = ALEA_NUC_ELAW_EVAPORATION;
    ASSERT_EQ(alea_nuc_sample_energy_distribution(
                  &law, 5.0, sequence_rng, &rng, &energy), ALEA_OK);
    ASSERT_TRUE(energy >= 0.0 && energy <= 4.0);

    law.law = ALEA_NUC_ELAW_WATT;
    law.n_watt_b = 2; law.watt_b_energy = grid; law.temp_C = one;
    ASSERT_EQ(alea_nuc_sample_energy_distribution(
                  &law, 5.0, sequence_rng, &rng, &energy), ALEA_OK);
    ASSERT_TRUE(energy >= 0.0 && energy <= 4.0);
}

TEST(energy_sampler_accepts_negative_fission_restriction) {
    double grid[] = {0.0, 20.0};
    double a[] = {1.0, 1.0};
    double b[] = {2.0, 2.0};
    alea_nuc_energy_dist_t law;
    memset(&law, 0, sizeof(law));
    law.law = ALEA_NUC_ELAW_WATT;
    law.n_temp = 2;
    law.temp_energy = grid;
    law.temp_T = a;
    law.n_watt_b = 2;
    law.watt_b_energy = grid;
    law.temp_C = b;
    law.restriction_energy = -30.0;
    double draws[] = {0.5, 0.5, 0.5, 0.5};
    sequence_rng_t rng = {draws, 4, 0};
    double energy = -1.0;

    ASSERT_EQ(alea_nuc_energy_dist_validate(&law, NULL), ALEA_OK);
    ASSERT_EQ(alea_nuc_sample_energy_distribution(
                  &law, 0.025, sequence_rng, &rng, &energy), ALEA_OK);
    ASSERT_TRUE(energy >= 0.0 && energy <= 30.025);
}

TEST(kalbach_sampler_preserves_energy_angle_correlation) {
    double ein[] = {1.0};
    int interpolation[] = {1}, n_discrete[] = {0}, n_eout[] = {2};
    double eout0[] = {0.0, 2.0}, pdf0[] = {0.5, 0.5};
    double cdf0[] = {0.0, 1.0}, r0[] = {1.0, 1.0}, a0[] = {1.0, 1.0};
    double* eout[] = {eout0}; double* pdf[] = {pdf0};
    double* cdf[] = {cdf0}; double* r[] = {r0}; double* a[] = {a0};
    alea_nuc_energy_dist_t law;
    memset(&law, 0, sizeof(law)); law.law = ALEA_NUC_ELAW_KALBACH;
    law.tab.n_ein=1; law.tab.ein=ein; law.tab.interpolation=interpolation;
    law.tab.n_discrete=n_discrete; law.tab.n_eout=n_eout;
    law.tab.eout=eout; law.tab.pdf=pdf; law.tab.cdf=cdf;
    law.tab.precompound_r=r; law.tab.precompound_a=a;
    ASSERT_EQ(alea_nuc_energy_dist_validate(&law, NULL), ALEA_OK);
    double draws[] = {0.5, 0.1, 0.5};
    sequence_rng_t rng = {draws, 3, 0};
    double energy = -1.0, mu = -2.0;
    bool correlated = false;
    ASSERT_EQ(alea_nuc_sample_energy_angle_distribution(
                  &law, 1.0, sequence_rng, &rng, &energy, &mu,
                  &correlated), ALEA_OK);
    ASSERT_NEAR(energy, 1.0, 1e-12);
    ASSERT_TRUE(correlated);
    ASSERT_NEAR(mu, 1.0 + log(0.5 + 0.5 * exp(-2.0)), 1e-12);

    r0[0] = r0[1] = -5e-4;
    ASSERT_EQ(alea_nuc_energy_dist_validate(&law, NULL), ALEA_OK);
    rng.position = 0;
    ASSERT_EQ(alea_nuc_sample_energy_angle_distribution(
                  &law, 1.0, sequence_rng, &rng, &energy, &mu,
                  &correlated), ALEA_OK);
    r0[0] = -0.01;
    ASSERT_EQ(alea_nuc_energy_dist_validate(&law, NULL),
              ALEA_ERR_INVALID_ARG);
}

TEST(law61_sampler_selects_angle_conditioned_on_outgoing_energy) {
    double ein[] = {1.0};
    int interpolation[] = {2}, n_discrete[] = {0}, n_eout[] = {2};
    double eout0[] = {0.0, 2.0}, pdf0[] = {0.5, 0.5};
    double cdf0[] = {0.0, 1.0};
    double* eout[] = {eout0}; double* pdf[] = {pdf0};
    double* cdf[] = {cdf0};
    double mu_grid[] = {0.0, 1.0}, mu_pdf[] = {1.0, 1.0};
    double mu_cdf[] = {0.0, 1.0};
    alea_nuc_angular_point_t angles[2];
    memset(angles, 0, sizeof(angles));
    angles[0].type = ALEA_NUC_ANG_ISOTROPIC;
    angles[1].type = ALEA_NUC_ANG_TABULAR;
    angles[1].interpolation = 1; angles[1].n_cosines = 2;
    angles[1].cosine = mu_grid; angles[1].pdf = mu_pdf;
    angles[1].cdf = mu_cdf;
    alea_nuc_angular_point_t* correlated[] = {angles};
    alea_nuc_energy_dist_t law;
    memset(&law, 0, sizeof(law)); law.law = ALEA_NUC_ELAW_CORRELATED;
    law.tab.n_ein=1; law.tab.ein=ein; law.tab.interpolation=interpolation;
    law.tab.n_discrete=n_discrete; law.tab.n_eout=n_eout;
    law.tab.eout=eout; law.tab.pdf=pdf; law.tab.cdf=cdf;
    law.tab.correlated_mu=correlated;
    ASSERT_EQ(alea_nuc_energy_dist_validate(&law, NULL), ALEA_OK);
    double draws[] = {0.75, 0.5, 0.5};
    sequence_rng_t rng = {draws, 3, 0};
    double sampled_energy = -1.0, mu = -2.0;
    bool correlated_angle = false;
    ASSERT_EQ(alea_nuc_sample_energy_angle_distribution(
                  &law, 1.0, sequence_rng, &rng, &sampled_energy, &mu,
                  &correlated_angle), ALEA_OK);
    ASSERT_NEAR(sampled_energy, 1.5, 1e-12);
    ASSERT_TRUE(correlated_angle);
    ASSERT_NEAR(mu, 0.5, 1e-12);
}

TEST(law67_sampler_selects_energy_conditioned_on_angle) {
    double incident_grid[] = {1.0};
    double cosine[] = {-1.0, 1.0};
    double low_energy[] = {0.0, 2.0}, low_pdf[] = {0.5, 0.5};
    double high_energy[] = {10.0, 14.0}, high_pdf[] = {0.25, 0.25};
    double cdf[] = {0.0, 1.0};
    alea_nuc_law67_energy_t spectra[2] = {
        {1, 2, low_energy, low_pdf, cdf},
        {1, 2, high_energy, high_pdf, cdf}
    };
    alea_nuc_law67_incident_t incident = {
        2, 2, cosine, spectra
    };
    alea_nuc_energy_dist_t law;
    memset(&law, 0, sizeof(law));
    law.law = ALEA_NUC_ELAW_LAB_ANGLE_ENERGY;
    law.law67.n_ein = 1;
    law.law67.ein = incident_grid;
    law.law67.incident = &incident;
    ASSERT_EQ(alea_nuc_energy_dist_validate(&law, NULL), ALEA_OK);

    double draws[] = {0.5, 0.75, 0.1, 0.2, 0.5};
    sequence_rng_t rng = {draws, 5, 0};
    double energy = -1.0, mu = -2.0;
    bool correlated = false;
    ASSERT_EQ(alea_nuc_sample_energy_angle_distribution(
                  &law, 1.0, sequence_rng, &rng, &energy, &mu,
                  &correlated), ALEA_OK);
    ASSERT_TRUE(correlated);
    ASSERT_NEAR(mu, 0.5, 1e-12);
    ASSERT_NEAR(energy, 12.0, 1e-12);
}

TEST(law67_decode_cleans_every_injected_allocation_failure) {
    double encoded[] = {
        0, 67, 8, 0, 1, 1.0, 1.0,
        0, 1, 1.0, 12,
        2, 2, -1.0, 1.0, 18, 26,
        1, 2, 0.0, 2.0, 0.5, 0.5, 0.0, 1.0,
        1, 2, 10.0, 14.0, 0.25, 0.25, 0.0, 1.0
    };
    alea_nuc_ace_table_t table;
    memset(&table, 0, sizeof(table));
    table.xss = encoded;
    table.xss_length = (int)(sizeof(encoded) / sizeof(encoded[0]));

    size_t successful_at = SIZE_MAX;
    for (size_t fail_at = 0; fail_at < 64; fail_at++) {
        table.allocation_error = false;
        nuc_alloc_failure_t failure = {fail_at, 0};
        alea_nuc_set_alloc_failure(fail_nuc_allocation, &failure);
        alea_nuc_energy_dist_t* distribution =
            alea_nuc_decode_energy_dist_base(&table, 1, 1);
        alea_nuc_set_alloc_failure(NULL, NULL);
        if (!table.allocation_error) {
            successful_at = fail_at;
            ASSERT_NOT_NULL(distribution);
            ASSERT_EQ(alea_nuc_energy_dist_validate(distribution, NULL),
                      ALEA_OK);
            alea_nuc_energy_dist_free(distribution);
            break;
        }
        ASSERT_TRUE(failure.calls > fail_at);
        alea_nuc_energy_dist_free(distribution);
    }
    ASSERT_NE(successful_at, SIZE_MAX);
    ASSERT_TRUE(successful_at >= 14);
}

TEST(nbody_sampler_uses_available_center_of_mass_energy) {
    alea_nuc_energy_dist_t law;
    memset(&law, 0, sizeof(law)); law.law = ALEA_NUC_ELAW_NBODY;
    law.nbody_particles = 3; law.nbody_total_mass = 3.0;
    law.nbody_target_awr = 2.0; law.nbody_q_value = 0.0;
    ASSERT_EQ(alea_nuc_energy_dist_validate(&law, NULL), ALEA_OK);
    double draws[] = {0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5};
    sequence_rng_t rng = {draws, 7, 0};
    double energy = -1.0, mu = -2.0;
    bool correlated = false;
    ASSERT_EQ(alea_nuc_sample_energy_angle_distribution(
                  &law, 3.0, sequence_rng, &rng, &energy, &mu,
                  &correlated), ALEA_OK);
    ASSERT_TRUE(correlated);
    ASSERT_NEAR(energy, 2.0 / 3.0, 1e-12);
    ASSERT_NEAR(mu, 0.0, 1e-12);
}

TEST(low_level_elastic_collision_samples_tabular_angle) {
    double energy[] = {0.5, 2.0};
    double incident_energy[] = {1.0};
    double cosine[] = {-1.0, 1.0};
    double pdf[] = {0.5, 0.5};
    double cdf[] = {0.0, 1.0};
    alea_nuc_angular_point_t point;
    alea_nuc_angular_dist_t angular;
    alea_nuc_nuclide_t nuc;
    double xi[] = {0.0, 0.25, 0.0};
    alea_nuc_interaction_t result;

    memset(&point, 0, sizeof(point));
    point.type = ALEA_NUC_ANG_TABULAR; point.interpolation = 1;
    point.n_cosines = 2; point.cosine = cosine; point.pdf = pdf; point.cdf = cdf;
    memset(&angular, 0, sizeof(angular));
    angular.n_energies = 1; angular.energy = incident_energy; angular.data = &point;
    memset(&nuc, 0, sizeof(nuc));
    nuc.awr = 1.0; nuc.n_energies = 2; nuc.energy = energy;
    nuc.elastic_angular = &angular;

    ASSERT_EQ(alea_nuc_sample_collision(&nuc, 2, 1.0, xi, &result), ALEA_OK);
    ASSERT_NEAR(result.mu, -0.5, 1e-12);
    ASSERT_NEAR(result.energy_out, 0.25, 1e-12);
    ASSERT_EQ(result.n_secondary, 1);
}

TEST(hydrogen_isotropic_elastic_samples_expected_energy_mean) {
    double energy[] = {0.5, 2.0};
    alea_nuc_nuclide_t nuc;
    memset(&nuc, 0, sizeof(nuc));
    nuc.awr = 1.0; nuc.n_energies = 2; nuc.energy = energy;
    double sum = 0.0;
    for (int i = 0; i < 1000; i++) {
        double xi[] = {0.0, (i + 0.5) / 1000.0, 0.0};
        alea_nuc_interaction_t result;
        ASSERT_EQ(alea_nuc_sample_collision(&nuc, 2, 1.0, xi, &result),
                  ALEA_OK);
        ASSERT_TRUE(result.energy_out >= 0.0 && result.energy_out <= 1.0);
        sum += result.energy_out;
    }
    ASSERT_NEAR(sum / 1000.0, 0.5, 1e-12);
}

TEST(low_level_elastic_inverts_linear_tabular_pdf) {
    double energy[] = {0.5, 2.0};
    double incident_energy[] = {1.0};
    double cosine[] = {-1.0, 1.0};
    double pdf[] = {0.0, 1.0};
    double cdf[] = {0.0, 1.0};
    alea_nuc_angular_point_t point = {
        ALEA_NUC_ANG_TABULAR, 2, 2, cosine, pdf, cdf
    };
    alea_nuc_angular_dist_t angular = {1, incident_energy, &point};
    alea_nuc_nuclide_t nuc;
    memset(&nuc, 0, sizeof(nuc));
    nuc.awr = 1.0; nuc.n_energies = 2; nuc.energy = energy;
    nuc.elastic_angular = &angular;
    double xi[] = {0.0, 0.25, 0.0};
    alea_nuc_interaction_t result;
    ASSERT_EQ(alea_nuc_sample_collision(&nuc, 2, 1.0, xi, &result), ALEA_OK);
    ASSERT_NEAR(result.mu, 0.0, 1e-12);
    ASSERT_NEAR(result.energy_out, 0.5, 1e-12);
}

TEST(prepared_collision_evaluates_flight_and_elastic_scatter) {
    double energy[] = {1.0, 3.0};
    double total[] = {2.0, 2.0};
    double elastic[] = {1.5, 1.5};
    double absorption[] = {0.5, 0.5};
    double capture[] = {0.5, 0.5};
    alea_nuc_reaction_t reaction;
    alea_nuc_nuclide_t nuc;
    alea_nuc_mat_component_t component;
    alea_nuc_material_t material;
    alea_nuc_prepare_requirements_t requirements = {
        .required_capabilities = ALEA_NUC_CAP_RESTRICTED_NEUTRON
    };
    alea_nuc_capability_report_t report;
    alea_nuc_prepared_material_t* prepared = NULL;

    memset(&reaction, 0, sizeof(reaction));
    reaction.mt = 102; reaction.threshold_index = 1;
    reaction.n_energies = 2; reaction.xs = capture;
    memset(&nuc, 0, sizeof(nuc));
    nuc.particle = ALEA_NUC_PARTICLE_NEUTRON; nuc.awr = 1.0;
    nuc.n_energies = 2; nuc.energy = energy; nuc.sigma_total = total;
    nuc.sigma_elastic = elastic; nuc.sigma_abs = absorption;
    nuc.n_reactions = 1; nuc.reactions = &reaction;
    component.nuclide = &nuc; component.number_density = 0.1;
    memset(&material, 0, sizeof(material));
    material.components = &component; material.n_components = 1;

    ASSERT_EQ(alea_nuc_prepare_material(&material, &requirements, &report,
                                        &prepared), ALEA_OK);
    ASSERT_NOT_NULL(prepared);
    ASSERT_EQ(report.available_capabilities, ALEA_NUC_CAP_RESTRICTED_NEUTRON);

    alea_nuc_particle_state_t incident = {
        .type = ALEA_NUC_PARTICLE_NEUTRON,
        .energy = 2.0,
        .direction = {0.0, 0.0, 1.0},
        .weight = 1.0,
        .time = 0.0
    };
    alea_nuc_evaluation_t evaluation;
    ASSERT_EQ(alea_nuc_evaluate(prepared, &incident, &evaluation), ALEA_OK);
    ASSERT_NEAR(evaluation.macro_total, 0.2, 1e-12);
    ASSERT_NEAR(evaluation.macro_elastic, 0.15, 1e-12);
    ASSERT_NEAR(evaluation.macro_absorption, 0.05, 1e-12);

    double flight_draws[] = {1.0 - exp(-2.0)};
    sequence_rng_t flight_rng = {flight_draws, 1, 0};
    double distance = 0.0;
    ASSERT_EQ(alea_nuc_sample_flight(&evaluation, sequence_rng, &flight_rng,
                                     &distance), ALEA_OK);
    ASSERT_NEAR(distance, 10.0, 1e-12);

    double collision_draws[] = {0.25, 0.25, 0.0, 0.75, 0.0};
    sequence_rng_t collision_rng = {collision_draws, 5, 0};
    alea_nuc_collision_result_t result;
    ASSERT_EQ(alea_nuc_collide(&evaluation, sequence_rng, &collision_rng,
                               &result), ALEA_OK);
    ASSERT_EQ(result.outcome, ALEA_NUC_OUTCOME_SCATTERED);
    ASSERT_EQ(result.component_index, 0);
    ASSERT_EQ(result.mt, 2);
    ASSERT_NEAR(result.mu_cm, 0.5, 1e-12);
    ASSERT_NEAR(result.outgoing.energy, 1.5, 1e-12);
    ASSERT_NEAR(result.outgoing.direction[0] * result.outgoing.direction[0] +
                result.outgoing.direction[1] * result.outgoing.direction[1] +
                result.outgoing.direction[2] * result.outgoing.direction[2],
                1.0, 1e-12);
    ASSERT_TRUE(result.deposition_available);
    ASSERT_NEAR(result.local_energy_deposition, 0.5, 1e-12);

    evaluation.macro_total = 0.25;
    result.mt = -1;
    ASSERT_EQ(alea_nuc_collide(&evaluation, sequence_rng, &collision_rng,
                               &result), ALEA_ERR_INVALID_STATE);
    ASSERT_EQ(result.mt, -1);

    alea_nuc_prepared_material_free(prepared);
}

TEST(prepared_collision_selects_absorption_channel) {
    double energy[] = {1.0, 3.0};
    double total[] = {2.0, 2.0};
    double elastic[] = {1.0, 1.0};
    double absorption[] = {1.0, 1.0};
    double capture[] = {1.0, 1.0};
    alea_nuc_reaction_t reaction = {
        .mt = 102, .threshold_index = 1, .n_energies = 2, .xs = capture
    };
    alea_nuc_nuclide_t nuc;
    memset(&nuc, 0, sizeof(nuc));
    nuc.particle = ALEA_NUC_PARTICLE_NEUTRON; nuc.awr = 12.0;
    nuc.n_energies = 2; nuc.energy = energy; nuc.sigma_total = total;
    nuc.sigma_elastic = elastic; nuc.sigma_abs = absorption;
    nuc.n_reactions = 1; nuc.reactions = &reaction;
    alea_nuc_mat_component_t component = {&nuc, 0.1};
    alea_nuc_material_t material = {&component, 1, 1};
    alea_nuc_prepare_requirements_t requirements = {
        .required_capabilities = ALEA_NUC_CAP_RESTRICTED_NEUTRON
    };
    alea_nuc_capability_report_t report;
    alea_nuc_prepared_material_t* prepared = NULL;
    ASSERT_EQ(alea_nuc_prepare_material(&material, &requirements, &report,
                                        &prepared), ALEA_OK);
    alea_nuc_particle_state_t incident = {
        ALEA_NUC_PARTICLE_NEUTRON, 2.0, {1.0, 0.0, 0.0}, 1.0, 0.0
    };
    alea_nuc_evaluation_t evaluation;
    ASSERT_EQ(alea_nuc_evaluate(prepared, &incident, &evaluation), ALEA_OK);
    double draws[] = {0.2, 0.9};
    sequence_rng_t rng = {draws, 2, 0};
    alea_nuc_collision_result_t result;
    ASSERT_EQ(alea_nuc_collide(&evaluation, sequence_rng, &rng, &result), ALEA_OK);
    ASSERT_EQ(result.outcome, ALEA_NUC_OUTCOME_ABSORBED);
    ASSERT_EQ(result.mt, 102);
    ASSERT_FALSE(result.deposition_available);
    ASSERT_TRUE(isnan(result.local_energy_deposition));
    alea_nuc_prepared_material_free(prepared);
}

TEST(prepared_collision_selects_target_by_macroscopic_total) {
    double energy[] = {1.0, 3.0};
    double total[] = {1.0, 1.0};
    double elastic[] = {1.0, 1.0};
    double absorption[] = {0.0, 0.0};
    alea_nuc_nuclide_t nuclides[2];
    memset(nuclides, 0, sizeof(nuclides));
    for (int i = 0; i < 2; i++) {
        nuclides[i].particle = ALEA_NUC_PARTICLE_NEUTRON;
        nuclides[i].awr = 1.0 + i;
        nuclides[i].n_energies = 2; nuclides[i].energy = energy;
        nuclides[i].sigma_total = total; nuclides[i].sigma_elastic = elastic;
        nuclides[i].sigma_abs = absorption;
    }
    alea_nuc_mat_component_t components[] = {
        {&nuclides[0], 0.25}, {&nuclides[1], 0.75}
    };
    alea_nuc_material_t material = {components, 2, 2};
    alea_nuc_prepare_requirements_t requirements = {
        .required_capabilities = ALEA_NUC_CAP_RESTRICTED_NEUTRON
    };
    alea_nuc_capability_report_t report;
    alea_nuc_prepared_material_t* prepared = NULL;
    ASSERT_EQ(alea_nuc_prepare_material(&material, &requirements, &report,
                                        &prepared), ALEA_OK);
    alea_nuc_particle_state_t incident = {
        ALEA_NUC_PARTICLE_NEUTRON, 2.0, {0.0, 0.0, 1.0}, 1.0, 0.0
    };
    alea_nuc_evaluation_t evaluation;
    ASSERT_EQ(alea_nuc_evaluate(prepared, &incident, &evaluation), ALEA_OK);
    double draws[] = {0.5, 0.5, 0.0, 0.5, 0.0};
    sequence_rng_t rng = {draws, 5, 0};
    alea_nuc_collision_result_t result;
    ASSERT_EQ(alea_nuc_collide(&evaluation, sequence_rng, &rng, &result),
              ALEA_OK);
    ASSERT_EQ(result.component_index, 1);
    alea_nuc_prepared_material_free(prepared);
}

TEST(preparation_rejects_inelastic_without_distribution_and_urr) {
    double energy[] = {1.0, 3.0};
    double total[] = {1.0, 1.0};
    double elastic[] = {0.0, 0.0};
    double absorption[] = {0.0, 0.0};
    double scatter[] = {1.0, 1.0};
    alea_nuc_reaction_t reaction = {
        .mt = 51, .ty = 1, .threshold_index = 1,
        .n_energies = 2, .xs = scatter
    };
    alea_nuc_nuclide_t nuc;
    memset(&nuc, 0, sizeof(nuc));
    nuc.particle = ALEA_NUC_PARTICLE_NEUTRON; nuc.awr = 56.0;
    nuc.n_energies = 2; nuc.energy = energy; nuc.sigma_total = total;
    nuc.sigma_elastic = elastic; nuc.sigma_abs = absorption;
    nuc.n_reactions = 1; nuc.reactions = &reaction;
    alea_nuc_mat_component_t component = {&nuc, 0.1};
    alea_nuc_material_t material = {&component, 1, 1};
    alea_nuc_prepare_requirements_t requirements = {
        .required_capabilities = ALEA_NUC_CAP_RESTRICTED_NEUTRON
    };
    alea_nuc_capability_report_t report;
    alea_nuc_prepared_material_t* prepared = NULL;
    ASSERT_EQ(alea_nuc_prepare_material(&material, &requirements, &report,
                                        &prepared), ALEA_ERR_UNSUPPORTED);
    ASSERT_NULL(prepared);
    ASSERT_EQ(report.issue, ALEA_NUC_PREP_INVALID_ENERGY_DISTRIBUTION);
    ASSERT_EQ(report.mt, 51);

    nuc.n_reactions = 0; nuc.reactions = NULL;
    nuc.sigma_total = elastic;
    alea_nuc_urr_t urr;
    memset(&urr, 0, sizeof(urr)); nuc.urr = &urr;
    ASSERT_EQ(alea_nuc_capabilities(&nuc, &report), ALEA_ERR_UNSUPPORTED);
    ASSERT_EQ(report.issue, ALEA_NUC_PREP_UNSUPPORTED_URR);
}

TEST(prepared_collision_emits_multiple_neutrons_into_caller_buffer) {
    double energy[] = {1.0, 3.0};
    double total[] = {1.0, 1.0};
    double elastic[] = {0.0, 0.0};
    double absorption[] = {0.0, 0.0};
    double production[] = {1.0, 1.0};
    alea_nuc_energy_dist_t law;
    memset(&law, 0, sizeof(law));
    law.law = ALEA_NUC_ELAW_LEVEL; law.level_Q = 0.5;
    alea_nuc_reaction_t reaction = {
        .mt = 16, .ty = 2, .threshold_index = 1,
        .n_energies = 2, .xs = production, .energy = &law
    };
    alea_nuc_nuclide_t nuc;
    memset(&nuc, 0, sizeof(nuc));
    nuc.particle = ALEA_NUC_PARTICLE_NEUTRON; nuc.awr = 12.0;
    nuc.n_energies = 2; nuc.energy = energy; nuc.sigma_total = total;
    nuc.sigma_elastic = elastic; nuc.sigma_abs = absorption;
    nuc.n_reactions = 1; nuc.reactions = &reaction;
    alea_nuc_mat_component_t component = {&nuc, 0.1};
    alea_nuc_material_t material = {&component, 1, 1};
    alea_nuc_prepare_requirements_t requirements = {
        .required_capabilities = ALEA_NUC_CAP_CONTINUOUS_NEUTRON
    };
    alea_nuc_capability_report_t report;
    alea_nuc_prepared_material_t* prepared = NULL;
    ASSERT_EQ(alea_nuc_prepare_material(&material, &requirements, &report,
                                        &prepared), ALEA_OK);
    ASSERT_TRUE(report.available_capabilities & ALEA_NUC_CAP_NEUTRON_EMISSION);
    alea_nuc_particle_state_t incident = {
        ALEA_NUC_PARTICLE_NEUTRON, 2.0, {0.0, 0.0, 1.0}, 1.0, 0.0
    };
    alea_nuc_evaluation_t evaluation;
    nuc_alloc_failure_t allocation = {0, 0};
    alea_nuc_set_alloc_failure(fail_nuc_allocation, &allocation);
    alea_error_t evaluation_status =
        alea_nuc_evaluate(prepared, &incident, &evaluation);
    alea_nuc_set_alloc_failure(NULL, NULL);
    ASSERT_EQ(evaluation_status, ALEA_OK);
    ASSERT_EQ(allocation.calls, 0);
    ASSERT_NEAR(evaluation.macro_neutron_emission, 0.1, 1e-12);

    alea_nuc_particle_state_t particles[2];
    alea_nuc_secondary_buffer_t buffer = {particles, 2, 0};
    double draws[] = {0.5, 0.5, 0.1, 0.25, 0.0, 0.9, 0.75, 0.5};
    sequence_rng_t rng = {draws, 8, 0};
    alea_nuc_collision_result_t result;
    allocation.calls = 0;
    alea_nuc_set_alloc_failure(fail_nuc_allocation, &allocation);
    alea_error_t collision_status = alea_nuc_collide_with_secondaries(
        &evaluation, sequence_rng, &rng, &buffer, &result);
    alea_nuc_set_alloc_failure(NULL, NULL);
    ASSERT_EQ(collision_status, ALEA_OK);
    ASSERT_EQ(allocation.calls, 0);
    ASSERT_EQ(result.outcome, ALEA_NUC_OUTCOME_REPLACED);
    ASSERT_EQ(result.mt, 16);
    ASSERT_EQ(result.n_emitted, 2);
    ASSERT_EQ(buffer.count, 2);
    ASSERT_NEAR(particles[0].energy, 1.0, 1e-12);
    ASSERT_NEAR(particles[1].energy, 1.0, 1e-12);
    alea_nuc_prepared_material_free(prepared);
}

TEST(prepared_collisions_are_worker_schedule_deterministic) {
    enum { HISTORIES = 257 };
    double energy[] = {1.0, 3.0};
    double elastic[] = {1.0, 1.0};
    double zero[] = {0.0, 0.0};
    alea_nuc_nuclide_t nuc;
    memset(&nuc, 0, sizeof(nuc));
    nuc.particle = ALEA_NUC_PARTICLE_NEUTRON;
    nuc.awr = 12.0;
    nuc.n_energies = 2;
    nuc.energy = energy;
    nuc.sigma_total = elastic;
    nuc.sigma_elastic = elastic;
    nuc.sigma_abs = zero;
    alea_nuc_mat_component_t component = {&nuc, 0.1};
    alea_nuc_material_t material = {&component, 1, 1};
    alea_nuc_prepare_requirements_t requirements = {
        .required_capabilities = ALEA_NUC_CAP_RESTRICTED_NEUTRON
    };
    alea_nuc_capability_report_t report;
    alea_nuc_prepared_material_t* prepared = NULL;
    ASSERT_EQ(alea_nuc_prepare_material(
                  &material, &requirements, &report, &prepared), ALEA_OK);

    scheduled_collision_t serial[HISTORIES];
    scheduled_collision_t parallel[HISTORIES];
    scheduled_collision_context_t serial_context = {prepared, serial};
    scheduled_collision_context_t parallel_context = {prepared, parallel};
    ASSERT_EQ(run_scheduled_collisions(
                  &serial_context, 0, 0, HISTORIES), 0);
    size_t actual_workers = 0;
    ASSERT_EQ(alea_parallel_for(
                  HISTORIES, 1, 4, ALEA_PARALLEL_DYNAMIC,
                  run_scheduled_collisions, &parallel_context,
                  &actual_workers), ALEA_PARALLEL_OK);
    ASSERT_TRUE(actual_workers >= 1 && actual_workers <= 4);
    for (size_t history = 0; history < HISTORIES; history++) {
        ASSERT_TRUE(scheduled_collisions_equal(
            &serial[history], &parallel[history]));
        ASSERT_EQ(serial[history].evaluation_status, ALEA_OK);
        ASSERT_EQ(serial[history].collision_status, ALEA_OK);
        ASSERT_EQ(serial[history].result.outcome, ALEA_NUC_OUTCOME_SCATTERED);
        ASSERT_EQ(serial[history].result.mt, 2);
    }
    alea_nuc_prepared_material_free(prepared);
}

TEST(secondary_capacity_failure_precedes_rng_and_publication) {
    double energy[] = {1.0, 3.0}, total[] = {1.0, 1.0};
    double zero[] = {0.0, 0.0}, production[] = {1.0, 1.0};
    alea_nuc_energy_dist_t law;
    memset(&law, 0, sizeof(law));
    law.law = ALEA_NUC_ELAW_LEVEL; law.level_Q = 0.5;
    alea_nuc_reaction_t reaction = {
        .mt=16, .ty=2, .threshold_index=1, .n_energies=2,
        .xs=production, .energy=&law
    };
    alea_nuc_nuclide_t nuc;
    memset(&nuc, 0, sizeof(nuc)); nuc.particle=ALEA_NUC_PARTICLE_NEUTRON;
    nuc.awr=12.0; nuc.n_energies=2; nuc.energy=energy;
    nuc.sigma_total=total; nuc.sigma_elastic=zero; nuc.sigma_abs=zero;
    nuc.n_reactions=1; nuc.reactions=&reaction;
    alea_nuc_mat_component_t component = {&nuc, 0.1};
    alea_nuc_material_t material = {&component, 1, 1};
    alea_nuc_prepare_requirements_t requirements = {
        .required_capabilities = ALEA_NUC_CAP_NEUTRON_EMISSION
    };
    alea_nuc_capability_report_t report;
    alea_nuc_prepared_material_t* prepared = NULL;
    ASSERT_EQ(alea_nuc_prepare_material(&material, &requirements, &report,
                                        &prepared), ALEA_OK);
    alea_nuc_particle_state_t incident = {
        ALEA_NUC_PARTICLE_NEUTRON, 2.0, {0,0,1}, 1.0, 0.0
    };
    alea_nuc_evaluation_t evaluation;
    ASSERT_EQ(alea_nuc_evaluate(prepared, &incident, &evaluation), ALEA_OK);
    alea_nuc_particle_state_t particle;
    alea_nuc_secondary_buffer_t buffer = {&particle, 1, 0};
    double draw = 0.5;
    sequence_rng_t rng = {&draw, 1, 0};
    alea_nuc_collision_result_t result;
    memset(&result, 0x5a, sizeof(result));
    alea_nuc_collision_result_t before = result;
    ASSERT_EQ(alea_nuc_collide_with_secondaries(
                  &evaluation, sequence_rng, &rng, &buffer, &result),
              ALEA_ERR_OUT_OF_MEMORY);
    ASSERT_EQ(rng.position, 0);
    ASSERT_EQ(buffer.count, 0);
    ASSERT_EQ(memcmp(&result, &before, sizeof(result)), 0);
    alea_nuc_prepared_material_free(prepared);
}

TEST(coordinated_urr_evaluation_drives_flight_and_reaction_probabilities) {
    double energy[] = {1.0, 3.0};
    double total[] = {2.0, 2.0}, elastic[] = {1.0, 1.0};
    double absorption[] = {1.0, 1.0}, capture[] = {1.0, 1.0};
    double urr_energy[] = {1.0, 3.0};
    double urr_table[] = {
        1.0, 0.7505, 0.5, 0.0, 1.0, 1.0,
        1.0, 0.7505, 0.5, 0.0, 1.0, 1.0
    };
    alea_nuc_urr_t urr;
    memset(&urr, 0, sizeof(urr)); urr.n_energies=2; urr.n_bands=1;
    urr.interp=2; urr.multiply_smooth=true;
    urr.energy=urr_energy; urr.table=urr_table;
    alea_nuc_reaction_t reaction = {
        .mt=102, .threshold_index=1, .n_energies=2, .xs=capture
    };
    alea_nuc_nuclide_t nuc;
    memset(&nuc, 0, sizeof(nuc)); nuc.particle=ALEA_NUC_PARTICLE_NEUTRON;
    nuc.awr=12.0; nuc.n_energies=2; nuc.energy=energy;
    nuc.sigma_total=total; nuc.sigma_elastic=elastic;
    nuc.sigma_abs=absorption; nuc.n_reactions=1; nuc.reactions=&reaction;
    nuc.urr=&urr;
    alea_nuc_mat_component_t component = {&nuc, 0.1};
    alea_nuc_material_t material = {&component, 1, 1};
    alea_nuc_prepare_requirements_t requirements = {
        .required_capabilities = ALEA_NUC_CAP_RESTRICTED_NEUTRON |
                                 ALEA_NUC_CAP_URR
    };
    alea_nuc_capability_report_t report;
    alea_nuc_prepared_material_t* prepared = NULL;
    ASSERT_EQ(alea_nuc_prepare_material(&material, &requirements, &report,
                                        &prepared), ALEA_OK);
    alea_nuc_particle_state_t incident = {
        ALEA_NUC_PARTICLE_NEUTRON, 2.0, {0,0,1}, 1.0, 0.0
    };
    alea_nuc_evaluation_t evaluation;
    ASSERT_EQ(alea_nuc_evaluate(prepared, &incident, &evaluation),
              ALEA_ERR_UNSUPPORTED);
    alea_nuc_urr_sample_t sample;
    alea_nuc_evaluation_workspace_t workspace = {&sample, 1};
    double urr_draw = 0.5;
    sequence_rng_t urr_rng = {&urr_draw, 1, 0};
    ASSERT_EQ(alea_nuc_evaluate_urr(prepared, &incident, sequence_rng,
                                    &urr_rng, &workspace, &evaluation), ALEA_OK);
    ASSERT_TRUE(sample.active);
    ASSERT_NEAR(evaluation.macro_total, 0.15, 1e-12);
    ASSERT_NEAR(evaluation.macro_elastic, 0.05, 1e-12);
    ASSERT_NEAR(evaluation.macro_absorption, 0.10, 1e-12);
    double collision_draws[] = {0.5, 0.9};
    sequence_rng_t collision_rng = {collision_draws, 2, 0};
    alea_nuc_collision_result_t result;
    ASSERT_EQ(alea_nuc_collide(&evaluation, sequence_rng, &collision_rng,
                               &result), ALEA_OK);
    ASSERT_EQ(result.outcome, ALEA_NUC_OUTCOME_ABSORBED);
    ASSERT_EQ(result.mt, 102);
    sample.factors[1] = 1.0;
    double distance_draw = 0.5, distance;
    sequence_rng_t distance_rng = {&distance_draw, 1, 0};
    ASSERT_EQ(alea_nuc_sample_flight(&evaluation, sequence_rng, &distance_rng,
                                     &distance), ALEA_ERR_INVALID_STATE);
    ASSERT_EQ(distance_rng.position, 0);
    alea_nuc_prepared_material_free(prepared);
}

TEST(delayed_fission_emission_selects_group_spectrum_and_time) {
    double energy[] = {1.0, 3.0}, total_xs[] = {1.0, 1.0};
    double zero[] = {0.0, 0.0}, fission_xs[] = {1.0, 1.0};
    double prompt_coeff[] = {1.0}, delayed_coeff[] = {1.0};
    alea_nuc_nu_bar_t prompt_nu, delayed_nu;
    memset(&prompt_nu, 0, sizeof(prompt_nu));
    memset(&delayed_nu, 0, sizeof(delayed_nu));
    prompt_nu.type=ALEA_NUC_NU_POLYNOMIAL; prompt_nu.n_coeffs=1;
    prompt_nu.coeffs=prompt_coeff;
    delayed_nu.type=ALEA_NUC_NU_POLYNOMIAL; delayed_nu.n_coeffs=1;
    delayed_nu.coeffs=delayed_coeff;
    alea_nuc_energy_dist_t prompt_spectrum, delayed_spectrum;
    memset(&prompt_spectrum, 0, sizeof(prompt_spectrum));
    memset(&delayed_spectrum, 0, sizeof(delayed_spectrum));
    prompt_spectrum.law=ALEA_NUC_ELAW_LEVEL; prompt_spectrum.level_Q=0.5;
    delayed_spectrum.law=ALEA_NUC_ELAW_LEVEL; delayed_spectrum.level_Q=0.25;
    double group_energy[] = {1.0}, group_probability[] = {1.0};
    alea_nuc_delayed_group_t group;
    memset(&group, 0, sizeof(group)); group.decay_rate=2.0;
    group.n_energies=1; group.energy=group_energy;
    group.probability=group_probability; group.spectrum=&delayed_spectrum;
    alea_nuc_fission_t fission;
    memset(&fission, 0, sizeof(fission)); fission.prompt=&prompt_nu;
    fission.delayed=&delayed_nu; fission.n_delayed_groups=1;
    fission.delayed_groups=&group;
    alea_nuc_reaction_t reaction = {
        .mt=18, .ty=19, .threshold_index=1, .n_energies=2,
        .xs=fission_xs, .energy=&prompt_spectrum
    };
    alea_nuc_nuclide_t nuc;
    memset(&nuc, 0, sizeof(nuc)); nuc.particle=ALEA_NUC_PARTICLE_NEUTRON;
    nuc.awr=235.0; nuc.n_energies=2; nuc.energy=energy;
    nuc.sigma_total=total_xs; nuc.sigma_elastic=zero; nuc.sigma_abs=zero;
    nuc.n_reactions=1; nuc.reactions=&reaction; nuc.fission=&fission;
    alea_nuc_mat_component_t component = {&nuc, 0.1};
    alea_nuc_material_t material = {&component, 1, 1};
    alea_nuc_prepare_requirements_t requirements = {
        .required_capabilities = ALEA_NUC_CAP_CONTINUOUS_NEUTRON |
                                 ALEA_NUC_CAP_FISSION |
                                 ALEA_NUC_CAP_DELAYED_NEUTRON
    };
    alea_nuc_capability_report_t report;
    alea_nuc_prepared_material_t* prepared = NULL;
    ASSERT_EQ(alea_nuc_prepare_material(&material, &requirements, &report,
                                        &prepared), ALEA_OK);
    alea_nuc_particle_state_t incident = {
        ALEA_NUC_PARTICLE_NEUTRON, 2.0, {0,0,1}, 1.0, 3.0
    };
    alea_nuc_evaluation_t evaluation;
    ASSERT_EQ(alea_nuc_evaluate(prepared, &incident, &evaluation), ALEA_OK);
    alea_nuc_particle_state_t particles[2];
    alea_nuc_secondary_buffer_t buffer = {particles, 2, 0};
    double draws[] = {0.5, 0.5, 0.1, 0.5, 0.0,
                      0.2, 0.5, 0.25, 0.5};
    sequence_rng_t rng = {draws, 9, 0};
    alea_nuc_collision_result_t result;
    ASSERT_EQ(alea_nuc_collide_with_secondaries(
                  &evaluation, sequence_rng, &rng, &buffer, &result), ALEA_OK);
    ASSERT_EQ(result.outcome, ALEA_NUC_OUTCOME_REPLACED);
    ASSERT_EQ(result.n_emitted, 2);
    ASSERT_NEAR(particles[0].energy, 1.0, 1e-12);
    ASSERT_NEAR(particles[0].time, 3.0, 1e-12);
    ASSERT_NEAR(particles[1].energy, 0.5, 1e-12);
    ASSERT_NEAR(particles[1].time, 3.0 + log(2.0) / 2.0, 1e-12);
    alea_nuc_prepared_material_free(prepared);
}

TEST(delayed_ace_blocks_decode_group_probability_and_spectrum) {
    double encoded[] = {
        2, 0, 1, 1.0, 0.1,
        2.0e-8, 0, 1, 1.0, 1.0,
        1,
        0, 3, 8, 0, 0, 0, 0, 0.0, 0.5
    };
    alea_nuc_nuclide_t* nuc = calloc(1, sizeof(*nuc));
    ASSERT_NOT_NULL(nuc);
    nuc->awr = 235.0;
    nuc->raw.awr = 235.0;
    nuc->raw.xss_length = (int)(sizeof(encoded) / sizeof(encoded[0]));
    nuc->raw.xss = malloc(sizeof(encoded));
    ASSERT_NOT_NULL(nuc->raw.xss);
    memcpy(nuc->raw.xss, encoded, sizeof(encoded));
    nuc->raw.nxs[7] = 1;
    nuc->raw.jxs[23] = 1;
    nuc->raw.jxs[24] = 6;
    nuc->raw.jxs[25] = 11;
    nuc->raw.jxs[26] = 12;
    ASSERT_EQ(alea_nuc_decode_delayed_neutrons(nuc, &nuc->raw), ALEA_OK);
    ASSERT_NOT_NULL(nuc->fission);
    ASSERT_EQ(nuc->fission->n_delayed_groups, 1);
    ASSERT_NEAR(nuc->fission->delayed_groups[0].decay_rate, 2.0, 1e-12);
    ASSERT_NEAR(nuc->fission->delayed_groups[0].probability[0], 1.0, 1e-12);
    ASSERT_NOT_NULL(nuc->fission->delayed_groups[0].spectrum);
    ASSERT_EQ(nuc->fission->delayed_groups[0].spectrum->law,
              ALEA_NUC_ELAW_LEVEL);
    ASSERT_NEAR(nuc->fission->delayed_groups[0].spectrum->level_Q,
                0.5, 1e-12);
    alea_nuc_nuclide_free(nuc);
}

TEST(delayed_ace_decode_cleans_every_injected_allocation_failure) {
    const double encoded[] = {
        2, 0, 1, 1.0, 0.1,
        2.0e-8, 1, 2, 2, 2, 1.0, 2.0, 0.25, 1.0,
        1,
        0, 3, 8, 0, 0, 0, 0, 0.0, 0.5
    };
    size_t successful_at = SIZE_MAX;

    for (size_t fail_at = 0; fail_at < 64; fail_at++) {
        alea_nuc_nuclide_t* nuc = calloc(1, sizeof(*nuc));
        ASSERT_NOT_NULL(nuc);
        nuc->awr = 235.0;
        nuc->raw.awr = 235.0;
        nuc->raw.xss_length = (int)(sizeof(encoded) / sizeof(encoded[0]));
        nuc->raw.xss = malloc(sizeof(encoded));
        ASSERT_NOT_NULL(nuc->raw.xss);
        memcpy(nuc->raw.xss, encoded, sizeof(encoded));
        nuc->raw.nxs[7] = 1;
        nuc->raw.jxs[23] = 1;
        nuc->raw.jxs[24] = 6;
        nuc->raw.jxs[25] = 15;
        nuc->raw.jxs[26] = 16;

        nuc_alloc_failure_t failure = {fail_at, 0};
        alea_nuc_set_alloc_failure(fail_nuc_allocation, &failure);
        alea_error_t err =
            alea_nuc_decode_delayed_neutrons(nuc, &nuc->raw);
        alea_nuc_set_alloc_failure(NULL, NULL);
        if (err == ALEA_OK) {
            successful_at = fail_at;
            ASSERT_EQ(nuc->fission->n_delayed_groups, 1);
            ASSERT_NOT_NULL(nuc->fission->delayed_groups[0].nbt);
            ASSERT_NOT_NULL(nuc->fission->delayed_groups[0].interp);
            ASSERT_NOT_NULL(nuc->fission->delayed_groups[0].energy);
            ASSERT_NOT_NULL(nuc->fission->delayed_groups[0].probability);
            ASSERT_NOT_NULL(nuc->fission->delayed_groups[0].spectrum);
            alea_nuc_nuclide_free(nuc);
            break;
        }
        ASSERT_EQ(err, ALEA_ERR_OUT_OF_MEMORY);
        ASSERT_TRUE(failure.calls > fail_at);
        alea_nuc_nuclide_free(nuc);
    }
    ASSERT_NE(successful_at, SIZE_MAX);
    ASSERT_TRUE(successful_at >= 9);
}

TEST(preparation_rejects_nonfinite_reaction_cross_section) {
    double energy[] = {1.0, 3.0};
    double total[] = {1.0, 1.0};
    double elastic[] = {0.0, 0.0};
    double absorption[] = {1.0, 1.0};
    double capture[] = {1.0, NAN};
    alea_nuc_reaction_t reaction = {
        .mt = 102, .threshold_index = 1, .n_energies = 2, .xs = capture
    };
    alea_nuc_nuclide_t nuc;
    memset(&nuc, 0, sizeof(nuc));
    nuc.particle = ALEA_NUC_PARTICLE_NEUTRON; nuc.awr = 12.0;
    nuc.n_energies = 2; nuc.energy = energy; nuc.sigma_total = total;
    nuc.sigma_elastic = elastic; nuc.sigma_abs = absorption;
    nuc.n_reactions = 1; nuc.reactions = &reaction;
    alea_nuc_capability_report_t report;
    ASSERT_EQ(alea_nuc_capabilities(&nuc, &report), ALEA_ERR_UNSUPPORTED);
    ASSERT_EQ(report.issue, ALEA_NUC_PREP_INVALID_CROSS_SECTIONS);
    ASSERT_EQ(report.mt, 102);
}

TEST(preparation_rejects_unknown_zero_yield_reaction) {
    double energy[] = {1.0, 3.0};
    double total[] = {1.0, 1.0};
    double elastic[] = {0.0, 0.0};
    double absorption[] = {1.0, 1.0};
    double unknown[] = {1.0, 1.0};
    alea_nuc_reaction_t reaction = {
        .mt = 999, .threshold_index = 1, .n_energies = 2, .xs = unknown
    };
    alea_nuc_nuclide_t nuc;
    memset(&nuc, 0, sizeof(nuc));
    nuc.particle = ALEA_NUC_PARTICLE_NEUTRON; nuc.awr = 12.0;
    nuc.n_energies = 2; nuc.energy = energy; nuc.sigma_total = total;
    nuc.sigma_elastic = elastic; nuc.sigma_abs = absorption;
    nuc.n_reactions = 1; nuc.reactions = &reaction;
    alea_nuc_capability_report_t report;
    ASSERT_EQ(alea_nuc_capabilities(&nuc, &report), ALEA_ERR_UNSUPPORTED);
    ASSERT_EQ(report.issue, ALEA_NUC_PREP_UNSUPPORTED_REACTION);
    ASSERT_EQ(report.mt, 999);
}

TEST(preparation_excludes_derived_responses_from_event_channels) {
    double energy[] = {1.0, 3.0};
    double total[] = {1.0, 1.0};
    double elastic[] = {1.0, 1.0};
    double zero[] = {0.0, 0.0};
    double proton_production[] = {2.0, 3.0};
    double damage[] = {50.0, 60.0};
    alea_nuc_reaction_t reactions[] = {
        {.mt = 3, .ty = 0, .threshold_index = 1,
         .n_energies = 2, .xs = proton_production},
        {.mt = 203, .ty = 0, .threshold_index = 1,
         .n_energies = 2, .xs = proton_production},
        {.mt = 444, .ty = 0, .threshold_index = 1,
         .n_energies = 2, .xs = damage}
    };
    alea_nuc_nuclide_t nuc;
    memset(&nuc, 0, sizeof(nuc));
    nuc.particle = ALEA_NUC_PARTICLE_NEUTRON;
    nuc.awr = 12.0;
    nuc.n_energies = 2;
    nuc.energy = energy;
    nuc.sigma_total = total;
    nuc.sigma_elastic = elastic;
    nuc.sigma_abs = zero;
    nuc.n_reactions = 3;
    nuc.reactions = reactions;
    alea_nuc_mat_component_t component = {&nuc, 0.1};
    alea_nuc_material_t material = {&component, 1, 1};
    alea_nuc_prepare_requirements_t requirements = {
        .required_capabilities = ALEA_NUC_CAP_RESTRICTED_NEUTRON
    };
    alea_nuc_capability_report_t report;
    alea_nuc_prepared_material_t* prepared = NULL;
    ASSERT_EQ(alea_nuc_prepare_material(
                  &material, &requirements, &report, &prepared), ALEA_OK);
    alea_nuc_particle_state_t incident = {
        ALEA_NUC_PARTICLE_NEUTRON, 2.0, {0.0, 0.0, 1.0}, 1.0, 0.0
    };
    alea_nuc_evaluation_t evaluation;
    ASSERT_EQ(alea_nuc_evaluate(prepared, &incident, &evaluation), ALEA_OK);
    ASSERT_NEAR(evaluation.macro_total, 0.1, 1e-14);
    double draws[] = {0.5, 0.5, 0.5, 0.5, 0.5};
    sequence_rng_t rng = {draws, 5, 0};
    alea_nuc_collision_result_t result;
    ASSERT_EQ(alea_nuc_collide(
                  &evaluation, sequence_rng, &rng, &result), ALEA_OK);
    ASSERT_EQ(result.mt, 2);
    alea_nuc_prepared_material_free(prepared);
}

TEST(preparation_accepts_right_continuous_main_grid) {
    double energy[] = {1.0, 2.0, 2.0, 3.0};
    double total[] = {1.0, 1.0, 2.0, 2.0};
    double zero[] = {0.0, 0.0, 0.0, 0.0};
    alea_nuc_nuclide_t nuc = {0};
    nuc.particle = ALEA_NUC_PARTICLE_NEUTRON;
    nuc.awr = 12.0;
    nuc.n_energies = 4;
    nuc.energy = energy;
    nuc.sigma_total = total;
    nuc.sigma_elastic = total;
    nuc.sigma_abs = zero;
    alea_nuc_mat_component_t component = {&nuc, 0.1};
    alea_nuc_material_t material = {&component, 1, 1};
    alea_nuc_prepare_requirements_t requirements = {
        .required_capabilities = ALEA_NUC_CAP_RESTRICTED_NEUTRON
    };
    alea_nuc_capability_report_t report;
    alea_nuc_prepared_material_t* prepared = NULL;
    ASSERT_EQ(alea_nuc_prepare_material(
                  &material, &requirements, &report, &prepared), ALEA_OK);
    alea_nuc_particle_state_t incident = {
        ALEA_NUC_PARTICLE_NEUTRON, 2.0, {0.0, 0.0, 1.0}, 1.0, 0.0
    };
    alea_nuc_evaluation_t evaluation;
    ASSERT_EQ(alea_nuc_evaluate(prepared, &incident, &evaluation), ALEA_OK);
    ASSERT_NEAR(evaluation.macro_total, 0.2, 1e-14);
    alea_nuc_prepared_material_free(prepared);
}

TEST(free_gas_elastic_converges_to_stationary_target_limit) {
    double energy[] = {1.0e-6, 2.0};
    alea_nuc_nuclide_t nuc;
    memset(&nuc, 0, sizeof(nuc));
    nuc.awr = 1.0; nuc.n_energies = 2; nuc.energy = energy;
    alea_nuc_particle_state_t incident = {
        ALEA_NUC_PARTICLE_NEUTRON, 1.0, {0.0, 0.0, 1.0}, 1.0, 4.0
    };
    double draws[] = {0.5, 0.5, 0.5, 0.5, 0.5,
                      0.5, 0.5, 0.5, 0.5, 0.5};
    sequence_rng_t rng = {draws, 10, 0};
    alea_nuc_free_gas_result_t result;
    ASSERT_EQ(alea_nuc_sample_free_gas_elastic(
                  &nuc, &incident, 1.0e-20, sequence_rng, &rng, &result),
              ALEA_OK);
    ASSERT_NEAR(result.mu_cm, 0.0, 1e-12);
    ASSERT_NEAR(result.outgoing.energy, 0.5, 1e-9);
    ASSERT_NEAR(result.mu_lab, sqrt(0.5), 1e-9);
    ASSERT_NEAR(result.outgoing.time, 4.0, 1e-12);
    ASSERT_NEAR(result.outgoing.direction[0] * result.outgoing.direction[0] +
                result.outgoing.direction[1] * result.outgoing.direction[1] +
                result.outgoing.direction[2] * result.outgoing.direction[2],
                1.0, 1e-12);
}

TEST(free_gas_elastic_rejects_invalid_angular_data_without_publishing) {
    double energy[] = {1.0e-6, 2.0};
    alea_nuc_angular_dist_t angular = {0};
    alea_nuc_nuclide_t nuc;
    memset(&nuc, 0, sizeof(nuc));
    nuc.awr = 1.0; nuc.n_energies = 2; nuc.energy = energy;
    nuc.elastic_angular = &angular;
    alea_nuc_particle_state_t incident = {
        ALEA_NUC_PARTICLE_NEUTRON, 1.0, {0.0, 0.0, 1.0}, 1.0, 0.0
    };
    double draw = 0.5;
    sequence_rng_t rng = {&draw, 1, 0};
    alea_nuc_free_gas_result_t result;
    memset(&result, 0, sizeof(result));
    result.mu_cm = 7.0;
    ASSERT_EQ(alea_nuc_sample_free_gas_elastic(
                  &nuc, &incident, 1.0e-8, sequence_rng, &rng, &result),
              ALEA_ERR_INVALID_STATE);
    ASSERT_EQ(rng.position, 0);
    ASSERT_NEAR(result.mu_cm, 7.0, 0.0);
}

TEST(free_gas_elastic_upscatters_subthermal_neutrons) {
    double energy[] = {1.0e-12, 1.0};
    alea_nuc_nuclide_t nuc;
    memset(&nuc, 0, sizeof(nuc));
    nuc.awr = 1.0; nuc.n_energies = 2; nuc.energy = energy;
    alea_nuc_particle_state_t incident = {
        ALEA_NUC_PARTICLE_NEUTRON, 1.0e-10, {0.0, 0.0, 1.0}, 1.0, 0.0
    };
    nucdata_rng_t rng = {UINT64_C(0x68e31da4f27c5b19)};
    int upscattered = 0;
    double energy_sum = 0.0;
    for (int i = 0; i < 2000; i++) {
        alea_nuc_free_gas_result_t result;
        ASSERT_EQ(alea_nuc_sample_free_gas_elastic(
                      &nuc, &incident, 2.5e-8, nucdata_rng, &rng, &result),
                  ALEA_OK);
        ASSERT_TRUE(isfinite(result.outgoing.energy));
        ASSERT_TRUE(result.outgoing.energy > 0.0);
        if (result.outgoing.energy > incident.energy) upscattered++;
        energy_sum += result.outgoing.energy;
    }
    ASSERT_TRUE(upscattered > 1800);
    ASSERT_TRUE(energy_sum / 2000.0 > 1.0e-8);
}

static double free_gas_mean_relative_speed_factor(double energy, double kT,
                                                   double awr) {
    double y = sqrt(awr * energy / kT);
    return (y + 0.5 / y) * erf(y) +
           exp(-y * y) / sqrt(3.14159265358979323846);
}

TEST(free_gas_elastic_preserves_thermal_equilibrium) {
    double energy_grid[] = {1.0e-12, 1.0};
    alea_nuc_nuclide_t nuc;
    memset(&nuc, 0, sizeof(nuc));
    nuc.awr = 1.0;
    nuc.n_energies = 2;
    nuc.energy = energy_grid;
    const double kT = 2.5e-8;
    double collision_energy = 0.0;
    double residence_weight = 0.0;
    double residence_energy = 0.0;
    int recorded = 0;

    for (uint32_t history = 0; history < 64; history++) {
        alea_nuc_rng_t rng;
        ASSERT_EQ(alea_nuc_rng_init(&rng, 9173, history, 0, 0,
                                   ALEA_NUC_RNG_COLLISION), ALEA_OK);
        alea_nuc_particle_state_t particle = {
            ALEA_NUC_PARTICLE_NEUTRON,
            kT * (0.1 + (history % 8)),
            {0.0, 0.0, 1.0}, 1.0, 0.0
        };
        for (int collision = 0; collision < 1200; collision++) {
            alea_nuc_free_gas_result_t result;
            ASSERT_EQ(alea_nuc_sample_free_gas_elastic(
                          &nuc, &particle, kT, alea_nuc_rng_uniform, &rng,
                          &result), ALEA_OK);
            particle = result.outgoing;
            if (collision < 200) continue;
            double weight = 1.0 / free_gas_mean_relative_speed_factor(
                particle.energy, kT, nuc.awr);
            collision_energy += particle.energy;
            residence_weight += weight;
            residence_energy += weight * particle.energy;
            recorded++;
        }
    }

    /* Collision states carry the mean-relative-speed bias. Dividing out that
     * analytic rate recovers the Maxwell equilibrium mean, 3 kT / 2. */
    ASSERT_NEAR(collision_energy / recorded / kT, 1.75, 0.04);
    ASSERT_NEAR(residence_energy / residence_weight / kT, 1.5, 0.03);
}

TEST(prepared_free_gas_elastic_uses_table_temperature_below_cutoff) {
    double energy[] = {1.0e-9, 1.0e-3};
    double total[] = {1.0, 1.0}, zero[] = {0.0, 0.0};
    alea_nuc_nuclide_t nuc;
    memset(&nuc, 0, sizeof(nuc));
    nuc.particle = ALEA_NUC_PARTICLE_NEUTRON;
    nuc.awr = 1.0; nuc.temperature = 2.5e-8;
    nuc.n_energies = 2; nuc.energy = energy;
    nuc.sigma_total = total; nuc.sigma_elastic = total; nuc.sigma_abs = zero;
    alea_nuc_mat_component_t component = {&nuc, 0.1};
    alea_nuc_material_t material = {&component, 1, 1};
    alea_nuc_prepare_requirements_t requirements = {
        .required_capabilities = ALEA_NUC_CAP_RESTRICTED_NEUTRON |
                                 ALEA_NUC_CAP_FREE_GAS
    };
    alea_nuc_capability_report_t report;
    alea_nuc_prepared_material_t* prepared = NULL;
    ASSERT_EQ(alea_nuc_prepare_material(&material, &requirements, &report,
                                        &prepared), ALEA_OK);
    alea_nuc_particle_state_t incident = {
        ALEA_NUC_PARTICLE_NEUTRON, 1.0e-6, {0.0, 0.0, 1.0}, 1.0, 0.0
    };
    alea_nuc_evaluation_t evaluation;
    ASSERT_EQ(alea_nuc_evaluate(prepared, &incident, &evaluation), ALEA_OK);
    double collision_draws[] = {0.5, 0.5, 0.5, 0.5, 0.5, 0.5,
                                0.5, 0.5, 0.5, 0.5, 0.5, 0.5};
    sequence_rng_t collision_rng = {collision_draws, 12, 0};
    alea_nuc_collision_result_t collision;
    ASSERT_EQ(alea_nuc_collide(&evaluation, sequence_rng, &collision_rng,
                               &collision), ALEA_OK);

    sequence_rng_t direct_rng = {collision_draws + 2, 10, 0};
    alea_nuc_free_gas_result_t direct;
    ASSERT_EQ(alea_nuc_sample_free_gas_elastic(
                  &nuc, &incident, nuc.temperature, sequence_rng,
                  &direct_rng, &direct), ALEA_OK);
    ASSERT_EQ(collision_rng.position, 12);
    ASSERT_NEAR(collision.outgoing.energy, direct.outgoing.energy, 1e-15);
    ASSERT_NEAR(collision.mu_cm, direct.mu_cm, 1e-15);
    ASSERT_NEAR(collision.mu_lab, direct.mu_lab, 1e-15);
    alea_nuc_prepared_material_free(prepared);
}

TEST(preparation_requires_temperature_for_every_free_gas_component) {
    double energy[] = {1.0e-9, 1.0e-3};
    double total[] = {1.0, 1.0}, zero[] = {0.0, 0.0};
    alea_nuc_nuclide_t nuclides[2];
    memset(nuclides, 0, sizeof(nuclides));
    for (int i = 0; i < 2; i++) {
        nuclides[i].particle = ALEA_NUC_PARTICLE_NEUTRON;
        nuclides[i].awr = i + 1.0; nuclides[i].n_energies = 2;
        nuclides[i].energy = energy; nuclides[i].sigma_total = total;
        nuclides[i].sigma_elastic = total; nuclides[i].sigma_abs = zero;
    }
    nuclides[0].temperature = 2.5e-8;
    alea_nuc_mat_component_t components[] = {
        {&nuclides[0], 0.1}, {&nuclides[1], 0.2}
    };
    alea_nuc_material_t material = {components, 2, 2};
    alea_nuc_prepare_requirements_t requirements = {
        .required_capabilities = ALEA_NUC_CAP_RESTRICTED_NEUTRON |
                                 ALEA_NUC_CAP_FREE_GAS
    };
    alea_nuc_capability_report_t report;
    alea_nuc_prepared_material_t* prepared = NULL;
    ASSERT_EQ(alea_nuc_prepare_material(&material, &requirements, &report,
                                        &prepared), ALEA_ERR_UNSUPPORTED);
    ASSERT_NULL(prepared);
    ASSERT_EQ(report.issue, ALEA_NUC_PREP_UNSUPPORTED_CAPABILITY);
    ASSERT_EQ(report.missing_capabilities, ALEA_NUC_CAP_FREE_GAS);
    ASSERT_EQ(report.component_index, 1);
    ASSERT_EQ(report.mt, 2);
}

static void synthetic_discrete_thermal_ace(alea_nuc_ace_table_t* table,
                                           double xss_values[43]) {
    memset(table, 0, sizeof(*table));
    memset(xss_values, 0, 43 * sizeof(*xss_values));
    snprintf(table->zaid, sizeof(table->zaid), "lwtr.20t");
    table->type = ALEA_NUC_TABLE_THERMAL_SAB;
    table->awr = 0.999167; table->temperature = 2.53e-8;
    table->iz[0] = 1001;
    table->xss = xss_values; table->xss_length = 43;
    table->nxs[0] = 43;
    table->nxs[2] = 1; /* two inelastic cosines */
    table->nxs[3] = 4; /* four outgoing energies */
    table->nxs[4] = 5; /* mixed elastic */
    table->nxs[6] = 0; /* equiprobable outgoing energies */
    table->nxs[7] = 1; /* two mixed-incoherent elastic cosines */
    table->jxs[0] = 1;  /* ITIE */
    table->jxs[2] = 6;  /* ITXE */
    table->jxs[3] = 30; /* ITCE: coherent */
    table->jxs[6] = 35; /* ITCE2: incoherent */
    table->jxs[8] = 40; /* ITCA2 */

    double header[] = {2.0, 1.0e-8, 1.0e-6, 4.0, 2.0};
    memcpy(xss_values, header, sizeof(header));
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 4; j++) {
            int base = 5 + (i * 4 + j) * 3;
            xss_values[base] = (double)(1 + j + 4 * i) * 1.0e-8;
            xss_values[base + 1] = 0.4 + 0.4 * i - 0.1 * j;
            xss_values[base + 2] = -0.6 + 0.2 * i + 0.1 * j;
        }
    }
    double coherent[] = {2.0, 1.0e-8, 5.0e-7, 1.0e-7, 3.0e-7};
    memcpy(&xss_values[29], coherent, sizeof(coherent));
    double incoherent[] = {2.0, 1.0e-8, 1.0e-6, 1.0, 3.0};
    memcpy(&xss_values[34], incoherent, sizeof(incoherent));
    double elastic_mu[] = {0.2, -0.8, 0.6, -0.4};
    memcpy(&xss_values[39], elastic_mu, sizeof(elastic_mu));
}

static void synthetic_continuous_thermal_ace(alea_nuc_ace_table_t* table,
                                             double xss_values[39]) {
    memset(table, 0, sizeof(*table));
    memset(xss_values, 0, 39 * sizeof(*xss_values));
    snprintf(table->zaid, sizeof(table->zaid), "lwtr.20t");
    table->type = ALEA_NUC_TABLE_THERMAL_SAB;
    table->awr = 0.999167; table->temperature = 2.53e-8;
    table->iz[0] = 1001;
    table->xss = xss_values; table->xss_length = 39;
    table->nxs[0] = 39;
    table->nxs[2] = 3; /* two continuous-inelastic cosines */
    table->nxs[3] = 2;
    table->nxs[4] = 0;
    table->nxs[6] = 2;
    table->jxs[0] = 1; /* ITIE */
    table->jxs[2] = 6; /* ITXE */

    double header[] = {2.0, 1.0e-8, 1.0e-6, 4.0, 2.0,
                       9.0, 24.0, 3.0, 3.0};
    memcpy(xss_values, header, sizeof(header));
    const double first[] = {
        1.0e-8, 1.25e7, 0.0, -0.5, 0.5,
        5.0e-8, 1.25e7, 0.5, -0.4, 0.6,
        9.0e-8, 1.25e7, 1.0, -0.3, 0.7
    };
    const double second[] = {
        2.0e-8, 1.25e7, 0.0, 0.0, 0.8,
        6.0e-8, 1.25e7, 0.5, 0.1, 0.9,
        1.0e-7, 1.25e7, 1.0, 0.2, 1.0
    };
    memcpy(&xss_values[9], first, sizeof(first));
    memcpy(&xss_values[24], second, sizeof(second));
}

static int write_discrete_thermal_ace(const char* path) {
    alea_nuc_ace_table_t table;
    double xss_values[43];
    synthetic_discrete_thermal_ace(&table, xss_values);
    FILE* fp = fopen(path, "w");
    if (!fp) return 0;
    fprintf(fp, "lwtr.20t %.12g %.12g 01/01/26\n", table.awr,
            table.temperature);
    fprintf(fp, "synthetic discrete thermal fixture\n");
    for (int row = 0; row < 4; row++) {
        for (int col = 0; col < 4; col++) {
            int i = row * 4 + col;
            fprintf(fp, "%d %.12g%c", table.iz[i],
                    i == 0 ? table.awr : 0.0, col == 3 ? '\n' : ' ');
        }
    }
    for (int row = 0; row < 2; row++)
        for (int col = 0; col < 8; col++)
            fprintf(fp, "%d%c", table.nxs[row * 8 + col],
                    col == 7 ? '\n' : ' ');
    for (int row = 0; row < 4; row++)
        for (int col = 0; col < 8; col++)
            fprintf(fp, "%d%c", table.jxs[row * 8 + col],
                    col == 7 ? '\n' : ' ');
    for (int i = 0; i < 43; i++)
        fprintf(fp, "%.17g%c", xss_values[i],
                i % 4 == 3 || i == 42 ? '\n' : ' ');
    return fclose(fp) == 0;
}

TEST(prepared_thermal_replaces_free_atom_elastic_below_table_cutoff) {
    alea_nuc_ace_table_t raw;
    double thermal_xss[43];
    synthetic_discrete_thermal_ace(&raw, thermal_xss);
    alea_nuc_thermal_t* thermal = NULL;
    ASSERT_EQ(alea_nuc_decode_thermal_internal(&raw, &thermal), ALEA_OK);

    double energy[] = {1.0e-8, 2.0e-6};
    double total[] = {10.0, 10.0};
    double elastic[] = {8.0, 8.0};
    double absorption[] = {2.0, 2.0};
    double capture[] = {2.0, 2.0};
    alea_nuc_reaction_t reaction = {
        .mt = 102, .threshold_index = 1, .n_energies = 2, .xs = capture
    };
    alea_nuc_nuclide_t nuc;
    memset(&nuc, 0, sizeof(nuc));
    nuc.Z = 1; nuc.A = 1; nuc.particle = ALEA_NUC_PARTICLE_NEUTRON;
    nuc.awr = 1.0; nuc.temperature = thermal->temperature;
    nuc.n_energies = 2; nuc.energy = energy; nuc.sigma_total = total;
    nuc.sigma_elastic = elastic; nuc.sigma_abs = absorption;
    nuc.n_reactions = 1; nuc.reactions = &reaction;
    alea_nuc_mat_component_t component = {&nuc, 0.1};
    alea_nuc_material_t material = {&component, 1, 1};
    alea_nuc_thermal_association_t association = {0, thermal};
    alea_nuc_prepare_requirements_t requirements = {
        .required_capabilities = ALEA_NUC_CAP_RESTRICTED_NEUTRON |
                                 ALEA_NUC_CAP_THERMAL_SAB,
        .thermal_associations = &association,
        .n_thermal_associations = 1,
        .thermal_temperature_tolerance = 0.0
    };
    alea_nuc_capability_report_t report;
    alea_nuc_prepared_material_t* prepared = NULL;
    ASSERT_EQ(alea_nuc_prepare_material(&material, &requirements, &report,
                                        &prepared), ALEA_OK);
    ASSERT_TRUE(report.available_capabilities & ALEA_NUC_CAP_THERMAL_SAB);

    alea_nuc_particle_state_t incident = {
        ALEA_NUC_PARTICLE_NEUTRON, 1.0e-6,
        {0.0, 0.0, 1.0}, 1.0, 0.0
    };
    alea_nuc_evaluation_t evaluation;
    ASSERT_EQ(alea_nuc_evaluate(prepared, &incident, &evaluation), ALEA_OK);
    ASSERT_NEAR(evaluation.macro_total, 0.73, 1e-12);
    ASSERT_NEAR(evaluation.macro_elastic, 0.0, 1e-12);
    ASSERT_NEAR(evaluation.macro_thermal, 0.53, 1e-12);
    ASSERT_NEAR(evaluation.macro_absorption, 0.2, 1e-12);

    double draws[] = {0.0, 0.1, 0.0, 0.6, 0.75, 0.0};
    sequence_rng_t rng = {draws, 6, 0};
    alea_nuc_collision_result_t result;
    ASSERT_EQ(alea_nuc_collide(&evaluation, sequence_rng, &rng, &result),
              ALEA_OK);
    ASSERT_EQ(result.component_index, 0);
    ASSERT_EQ(result.outcome, ALEA_NUC_OUTCOME_SCATTERED);
    ASSERT_EQ(result.mt, 4);
    ASSERT_EQ(rng.position, 6);

    sequence_rng_t secondary_rng = {draws, 6, 0};
    alea_nuc_secondary_buffer_t secondaries = {NULL, 0, 0};
    ASSERT_EQ(alea_nuc_collide_with_secondaries(
                  &evaluation, sequence_rng, &secondary_rng,
                  &secondaries, &result), ALEA_OK);
    ASSERT_EQ(result.component_index, 0);
    ASSERT_EQ(result.mt, 4);
    ASSERT_EQ(secondaries.count, 0);
    ASSERT_EQ(secondary_rng.position, 6);

    incident.energy = 1.5e-6;
    ASSERT_EQ(alea_nuc_evaluate(prepared, &incident, &evaluation), ALEA_OK);
    ASSERT_NEAR(evaluation.macro_total, 1.0, 1e-12);
    ASSERT_NEAR(evaluation.macro_elastic, 0.8, 1e-12);
    ASSERT_NEAR(evaluation.macro_thermal, 0.0, 1e-12);

    alea_nuc_prepared_material_free(prepared);
    alea_nuc_thermal_free(thermal);
}

TEST(prepared_thermal_rejects_wrong_nuclide_temperature_and_duplicates) {
    alea_nuc_ace_table_t raw;
    double thermal_xss[43];
    synthetic_discrete_thermal_ace(&raw, thermal_xss);
    alea_nuc_thermal_t* thermal = NULL;
    ASSERT_EQ(alea_nuc_decode_thermal_internal(&raw, &thermal), ALEA_OK);

    double energy[] = {1.0e-8, 2.0e-6};
    double total[] = {1.0, 1.0}, zero[] = {0.0, 0.0};
    alea_nuc_nuclide_t nuc;
    memset(&nuc, 0, sizeof(nuc));
    nuc.Z = 8; nuc.A = 16; nuc.particle = ALEA_NUC_PARTICLE_NEUTRON;
    nuc.awr = 16.0; nuc.temperature = thermal->temperature;
    nuc.n_energies = 2; nuc.energy = energy; nuc.sigma_total = total;
    nuc.sigma_elastic = total; nuc.sigma_abs = zero;
    alea_nuc_mat_component_t component = {&nuc, 0.1};
    alea_nuc_material_t material = {&component, 1, 1};
    alea_nuc_thermal_association_t associations[] = {
        {0, thermal}, {0, thermal}
    };
    alea_nuc_prepare_requirements_t requirements = {
        .required_capabilities = ALEA_NUC_CAP_RESTRICTED_NEUTRON |
                                 ALEA_NUC_CAP_THERMAL_SAB,
        .thermal_associations = associations,
        .n_thermal_associations = 1
    };
    alea_nuc_capability_report_t report;
    alea_nuc_prepared_material_t* prepared = NULL;
    requirements.n_thermal_associations = 0;
    ASSERT_EQ(alea_nuc_prepare_material(&material, &requirements, &report,
                                        &prepared), ALEA_ERR_UNSUPPORTED);
    ASSERT_EQ(report.missing_capabilities, ALEA_NUC_CAP_THERMAL_SAB);
    ASSERT_NULL(prepared);

    requirements.n_thermal_associations = 1;
    ASSERT_EQ(alea_nuc_prepare_material(&material, &requirements, &report,
                                        &prepared), ALEA_ERR_UNSUPPORTED);
    ASSERT_EQ(report.issue, ALEA_NUC_PREP_INVALID_THERMAL_ASSOCIATION);
    ASSERT_NULL(prepared);

    nuc.Z = 1; nuc.A = 1; nuc.temperature += 1.0e-10;
    ASSERT_EQ(alea_nuc_prepare_material(&material, &requirements, &report,
                                        &prepared), ALEA_ERR_UNSUPPORTED);
    ASSERT_EQ(report.issue, ALEA_NUC_PREP_INVALID_THERMAL_ASSOCIATION);
    requirements.thermal_temperature_tolerance = 1.0e-9;
    requirements.n_thermal_associations = 2;
    ASSERT_EQ(alea_nuc_prepare_material(&material, &requirements, &report,
                                        &prepared), ALEA_ERR_UNSUPPORTED);
    ASSERT_EQ(report.issue, ALEA_NUC_PREP_INVALID_THERMAL_ASSOCIATION);
    ASSERT_NULL(prepared);

    thermal->applicable_zaids[0] = 6012;
    nuc.Z = 6; nuc.A = 13; nuc.temperature = thermal->temperature;
    requirements.thermal_temperature_tolerance = 0.0;
    requirements.n_thermal_associations = 1;
    ASSERT_EQ(alea_nuc_prepare_material(&material, &requirements, &report,
                                        &prepared), ALEA_OK);
    alea_nuc_prepared_material_free(prepared);

    nuc.A = 14;
    ASSERT_EQ(alea_nuc_prepare_material(&material, &requirements, &report,
                                        &prepared), ALEA_ERR_UNSUPPORTED);
    ASSERT_EQ(report.issue, ALEA_NUC_PREP_INVALID_THERMAL_ASSOCIATION);
    ASSERT_NULL(prepared);

    alea_nuc_thermal_free(thermal);
}

TEST(load_symbolic_discrete_thermal_ace_table) {
    const char* path = "nucdata_thermal_fixture.tmp";
    ASSERT_TRUE(write_discrete_thermal_ace(path));
    alea_nuc_xsdir_entry_t entry;
    memset(&entry, 0, sizeof(entry));
    snprintf(entry.zaid, sizeof(entry.zaid), "lwtr.20t");
    snprintf(entry.filename, sizeof(entry.filename), "%s", path);
    entry.type = ALEA_NUC_TABLE_THERMAL_SAB;
    entry.file_type = 1; entry.address = 1;
    alea_nuc_xsdir_t xsdir;
    memset(&xsdir, 0, sizeof(xsdir));
    xsdir.entries = &entry; xsdir.count = 1;
    alea_nuc_thermal_t* thermal = alea_nuc_load_thermal(&xsdir, "lwtr.20t");
    ASSERT_NOT_NULL(thermal);
    ASSERT_EQ(thermal->elastic_mode, ALEA_NUC_THERMAL_ELASTIC_MIXED);
    ASSERT_NEAR(thermal->temperature, 2.53e-8, 1e-20);
    alea_nuc_thermal_free(thermal);
    remove(path);
}

TEST(thermal_discrete_ace_decode_preserves_channels_and_applicability) {
    alea_nuc_ace_table_t table;
    double xss_values[43];
    synthetic_discrete_thermal_ace(&table, xss_values);
    alea_nuc_thermal_t* thermal = NULL;
    ASSERT_EQ(alea_nuc_decode_thermal_internal(&table, &thermal), ALEA_OK);
    ASSERT_NOT_NULL(thermal);
    ASSERT_STR_EQ(thermal->zaid, "lwtr.20t");
    ASSERT_EQ(thermal->n_applicable_zaids, 1);
    ASSERT_EQ(thermal->applicable_zaids[0], 1001);
    ASSERT_EQ(thermal->elastic_mode, ALEA_NUC_THERMAL_ELASTIC_MIXED);
    ASSERT_EQ(thermal->n_inelastic_energies, 2);
    ASSERT_EQ(thermal->n_inelastic_outgoing, 4);
    ASSERT_EQ(thermal->n_inelastic_cosines, 2);
    ASSERT_NEAR(thermal->inelastic_mu[0], -0.6, 1e-15);
    ASSERT_NEAR(thermal->inelastic_mu[1], 0.4, 1e-15);
    ASSERT_NEAR(thermal->incoherent_mu[0], -0.8, 1e-15);
    ASSERT_NEAR(thermal->incoherent_mu[1], 0.2, 1e-15);
    ASSERT_NEAR(alea_nuc_thermal_xs_inelastic(thermal, 1.0e-6), 2.0,
                1e-14);
    ASSERT_NEAR(alea_nuc_thermal_xs_elastic(thermal, 1.0e-6), 3.3,
                1e-14);
    ASSERT_NEAR(alea_nuc_thermal_xs_elastic(thermal, 1.0e-8), 11.0,
                1e-14);
    ASSERT_NEAR(alea_nuc_thermal_xs_elastic(thermal, 5.0e-7), 2.5898989899,
                1e-10);
    ASSERT_NEAR(alea_nuc_thermal_xs_total(thermal, 2.0e-6), 0.15, 1e-14);
    alea_nuc_thermal_free(thermal);
}

TEST(thermal_discrete_collision_interpolates_correlated_quantiles) {
    alea_nuc_ace_table_t table;
    double xss_values[43];
    synthetic_discrete_thermal_ace(&table, xss_values);
    alea_nuc_thermal_t* thermal = NULL;
    ASSERT_EQ(alea_nuc_decode_thermal_internal(&table, &thermal), ALEA_OK);
    alea_nuc_particle_state_t incident = {
        ALEA_NUC_PARTICLE_NEUTRON, 5.05e-7, {0.0, 0.0, 1.0}, 2.0, 7.0
    };
    double draws[] = {0.0, 0.60, 0.75, 0.0};
    sequence_rng_t rng = {draws, 4, 0};
    alea_nuc_collision_result_t result;
    memset(&result, 0xA5, sizeof(result));
    ASSERT_EQ(alea_nuc_sample_thermal_collision(
                  thermal, &incident, sequence_rng, &rng, &result), ALEA_OK);
    ASSERT_EQ(rng.position, 4);
    ASSERT_EQ(result.outcome, ALEA_NUC_OUTCOME_SCATTERED);
    ASSERT_EQ(result.mt, 4);
    ASSERT_NEAR(result.outgoing.energy, 5.0e-8, 1e-18);
    ASSERT_NEAR(result.mu_lab, 0.4, 1e-14);
    ASSERT_NEAR(result.outgoing.direction[2], 0.4, 1e-14);
    ASSERT_NEAR(result.outgoing.weight, 2.0, 1e-15);
    ASSERT_NEAR(result.outgoing.time, 7.0, 1e-15);
    ASSERT_TRUE(result.deposition_available);
    ASSERT_NEAR(result.local_energy_deposition, 4.55e-7, 1e-18);
    alea_nuc_thermal_free(thermal);
}

TEST(thermal_coherent_collision_samples_bragg_edge) {
    alea_nuc_ace_table_t table;
    double xss_values[43];
    synthetic_discrete_thermal_ace(&table, xss_values);
    alea_nuc_thermal_t* thermal = NULL;
    ASSERT_EQ(alea_nuc_decode_thermal_internal(&table, &thermal), ALEA_OK);
    alea_nuc_particle_state_t incident = {
        ALEA_NUC_PARTICLE_NEUTRON, 1.0e-6, {1.0, 0.0, 0.0}, 1.0, 0.0
    };
    double draws[] = {2.1 / 5.3, 0.75, 0.25};
    sequence_rng_t rng = {draws, 3, 0};
    alea_nuc_collision_result_t result;
    ASSERT_EQ(alea_nuc_sample_thermal_collision(
                  thermal, &incident, sequence_rng, &rng, &result), ALEA_OK);
    ASSERT_EQ(result.mt, 2);
    ASSERT_NEAR(result.outgoing.energy, incident.energy, 1e-18);
    ASSERT_NEAR(result.mu_lab, 0.0, 1e-14);
    ASSERT_NEAR(result.outgoing.direction[0], 0.0, 1e-14);
    alea_nuc_thermal_free(thermal);
}

TEST(thermal_mixed_collision_samples_incoherent_elastic_quantile) {
    alea_nuc_ace_table_t table;
    double xss_values[43];
    synthetic_discrete_thermal_ace(&table, xss_values);
    alea_nuc_thermal_t* thermal = NULL;
    ASSERT_EQ(alea_nuc_decode_thermal_internal(&table, &thermal), ALEA_OK);
    alea_nuc_particle_state_t incident = {
        ALEA_NUC_PARTICLE_NEUTRON, 1.0e-6, {0.0, 1.0, 0.0}, 1.0, 0.0
    };
    double draws[] = {0.9, 0.75, 0.0};
    sequence_rng_t rng = {draws, 3, 0};
    alea_nuc_collision_result_t result;
    ASSERT_EQ(alea_nuc_sample_thermal_collision(
                  thermal, &incident, sequence_rng, &rng, &result), ALEA_OK);
    ASSERT_EQ(result.mt, 2);
    ASSERT_NEAR(result.mu_lab, 0.6, 1e-14);
    ASSERT_NEAR(result.outgoing.energy, incident.energy, 1e-18);
    alea_nuc_thermal_free(thermal);
}

TEST(thermal_skewed_inelastic_uses_ace_edge_probabilities) {
    alea_nuc_ace_table_t table;
    double xss_values[43];
    synthetic_discrete_thermal_ace(&table, xss_values);
    table.nxs[6] = 1;
    table.nxs[4] = 0;
    alea_nuc_thermal_t* thermal = NULL;
    ASSERT_EQ(alea_nuc_decode_thermal_internal(&table, &thermal), ALEA_OK);
    alea_nuc_particle_state_t incident = {
        ALEA_NUC_PARTICLE_NEUTRON, 1.0e-8, {0.0, 0.0, 1.0}, 1.0, 0.0
    };
    const double quantiles[] = {0.05, 0.30, 0.55, 0.80};
    const double expected[] = {1.0e-8, 2.0e-8, 4.0e-8, 3.0e-8};
    for (int i = 0; i < 4; i++) {
        double draws[] = {0.0, quantiles[i], 0.0, 0.0};
        sequence_rng_t rng = {draws, 4, 0};
        alea_nuc_collision_result_t result;
        ASSERT_EQ(alea_nuc_sample_thermal_collision(
                      thermal, &incident, sequence_rng, &rng, &result),
                  ALEA_OK);
        ASSERT_EQ(result.mt, 4);
        ASSERT_NEAR(result.outgoing.energy, expected[i], 1e-18);
    }
    alea_nuc_thermal_free(thermal);
}

TEST(thermal_continuous_decode_and_correlated_sampling) {
    alea_nuc_ace_table_t table;
    double xss_values[39];
    synthetic_continuous_thermal_ace(&table, xss_values);
    alea_nuc_thermal_t* thermal = NULL;
    ASSERT_EQ(alea_nuc_decode_thermal_internal(&table, &thermal), ALEA_OK);
    ASSERT_TRUE(thermal->inelastic_continuous);
    ASSERT_EQ(thermal->n_inelastic_outgoing, 0);
    ASSERT_EQ(thermal->n_inelastic_outgoing_total, 6);
    ASSERT_EQ(thermal->inelastic_outgoing_offset[0], 0);
    ASSERT_EQ(thermal->inelastic_outgoing_offset[1], 3);
    ASSERT_EQ(thermal->inelastic_outgoing_offset[2], 6);
    ASSERT_NEAR(thermal->inelastic_pdf[1], 1.25e7, 1e-8);

    alea_nuc_particle_state_t incident = {
        ALEA_NUC_PARTICLE_NEUTRON, 5.05e-7,
        {0.0, 0.0, 1.0}, 1.0, 0.0
    };
    double draws[] = {0.0, 0.25, 0.1, 0.5, 0.0};
    sequence_rng_t rng = {draws, 5, 0};
    alea_nuc_collision_result_t result;
    ASSERT_EQ(alea_nuc_sample_thermal_collision(
                  thermal, &incident, sequence_rng, &rng, &result), ALEA_OK);
    ASSERT_EQ(rng.position, 5);
    ASSERT_EQ(result.mt, 4);
    ASSERT_NEAR(result.outgoing.energy, 5.25e-7, 1e-18);
    ASSERT_NEAR(result.mu_lab, -0.45, 1e-14);
    ASSERT_NEAR(result.outgoing.direction[2], -0.45, 1e-14);

    double neutron_energy[] = {1.0e-8, 2.0e-6};
    double one[] = {1.0, 1.0}, zero[] = {0.0, 0.0};
    alea_nuc_nuclide_t nuc;
    memset(&nuc, 0, sizeof(nuc));
    nuc.Z = 1; nuc.A = 1; nuc.particle = ALEA_NUC_PARTICLE_NEUTRON;
    nuc.awr = 1.0; nuc.temperature = thermal->temperature;
    nuc.n_energies = 2; nuc.energy = neutron_energy;
    nuc.sigma_total = one; nuc.sigma_elastic = one; nuc.sigma_abs = zero;
    alea_nuc_mat_component_t component = {&nuc, 0.1};
    alea_nuc_material_t material = {&component, 1, 1};
    alea_nuc_thermal_association_t association = {0, thermal};
    alea_nuc_prepare_requirements_t requirements = {
        .required_capabilities = ALEA_NUC_CAP_RESTRICTED_NEUTRON |
                                 ALEA_NUC_CAP_THERMAL_SAB,
        .thermal_associations = &association,
        .n_thermal_associations = 1
    };
    alea_nuc_capability_report_t report;
    alea_nuc_prepared_material_t* prepared = NULL;
    ASSERT_EQ(alea_nuc_prepare_material(&material, &requirements, &report,
                                        &prepared), ALEA_OK);
    alea_nuc_evaluation_t evaluation;
    ASSERT_EQ(alea_nuc_evaluate(prepared, &incident, &evaluation), ALEA_OK);
    ASSERT_NEAR(evaluation.macro_total, 0.3, 1e-14);
    ASSERT_NEAR(evaluation.macro_thermal, 0.3, 1e-14);
    ASSERT_NEAR(evaluation.macro_elastic, 0.0, 1e-14);
    double prepared_draws[] = {0.0, 0.0, 0.0, 0.25, 0.1, 0.5, 0.0};
    sequence_rng_t prepared_rng = {prepared_draws, 7, 0};
    ASSERT_EQ(alea_nuc_collide(&evaluation, sequence_rng, &prepared_rng,
                               &result), ALEA_OK);
    ASSERT_EQ(result.component_index, 0);
    ASSERT_NEAR(result.outgoing.energy, 5.25e-7, 1e-18);
    ASSERT_EQ(prepared_rng.position, 7);
    alea_nuc_prepared_material_free(prepared);
    alea_nuc_thermal_free(thermal);
}

TEST(thermal_continuous_inserts_zero_cdf_endpoint_and_rejects_bad_locator) {
    alea_nuc_ace_table_t table;
    double xss_values[39];
    synthetic_continuous_thermal_ace(&table, xss_values);
    xss_values[11] = 0.1;
    alea_nuc_thermal_t* thermal = NULL;
    ASSERT_EQ(alea_nuc_decode_thermal_internal(&table, &thermal), ALEA_OK);
    ASSERT_EQ(thermal->n_inelastic_outgoing_total, 7);
    ASSERT_EQ(thermal->inelastic_outgoing_offset[1], 4);
    ASSERT_EQ(thermal->inelastic_energy_out[0], 0.0);
    ASSERT_EQ(thermal->inelastic_cdf[0], 0.0);
    ASSERT_NEAR(thermal->inelastic_mu[0], -0.5, 1e-15);
    ASSERT_NEAR(thermal->inelastic_mu[1], 0.5, 1e-15);
    alea_nuc_thermal_free(thermal);

    synthetic_continuous_thermal_ace(&table, xss_values);
    xss_values[5] = 39.0;
    thermal = (void*)1;
    ASSERT_EQ(alea_nuc_decode_thermal_internal(&table, &thermal),
              ALEA_ERR_PARSE_ERROR);
    ASSERT_NULL(thermal);

    synthetic_continuous_thermal_ace(&table, xss_values);
    xss_values[36] = 0.8;
    thermal = (void*)1;
    ASSERT_EQ(alea_nuc_decode_thermal_internal(&table, &thermal),
              ALEA_ERR_PARSE_ERROR);
    ASSERT_NULL(thermal);
}

TEST(thermal_continuous_clamps_only_roundoff_negative_pdf) {
    alea_nuc_ace_table_t table;
    double xss_values[39];
    synthetic_continuous_thermal_ace(&table, xss_values);
    xss_values[10] = -1.0e-27;
    alea_nuc_thermal_t* thermal = NULL;
    ASSERT_EQ(alea_nuc_decode_thermal_internal(&table, &thermal), ALEA_OK);
    ASSERT_EQ(thermal->inelastic_pdf[0], 0.0);
    alea_nuc_thermal_free(thermal);

    synthetic_continuous_thermal_ace(&table, xss_values);
    xss_values[10] = -1.0e-6;
    thermal = (void*)1;
    ASSERT_EQ(alea_nuc_decode_thermal_internal(&table, &thermal),
              ALEA_ERR_PARSE_ERROR);
    ASSERT_NULL(thermal);
}

static alea_nuc_nuclide_t* photon_production_allocation_fixture(void) {
    static const double values[] = {
        102001,                 /* MTRP */
        1,                      /* LSIGP */
        12, 102,                /* MF=12, parent MT */
        1, 2, 2, 2,             /* one lin-lin yield region */
        1.0, 2.0, 0.5, 1.5,    /* yield energy and values */
        1,                      /* LANDP */
        1,                      /* LDLWP */
        0, 2, 10, 0, 2,        /* law-2 header and applicability */
        1.0, 2.0, 1.0, 1.0,
        2, 0.9,                 /* primary discrete photon */
        1, 102,                 /* YP */
        1, 1.0, 4,              /* one angular incident energy */
        2, 2, -1.0, 1.0,       /* lin-lin angular table */
        0.5, 0.5, 0.0, 1.0,
        0.25, 0.5               /* GPD */
    };
    alea_nuc_nuclide_t* nuc = calloc(1, sizeof(*nuc));
    if (!nuc) return NULL;
    nuc->raw.xss = malloc(sizeof(values));
    if (!nuc->raw.xss) {
        free(nuc);
        return NULL;
    }
    memcpy(nuc->raw.xss, values, sizeof(values));
    nuc->particle = ALEA_NUC_PARTICLE_NEUTRON;
    nuc->awr = 10.0;
    nuc->n_energies = 2;
    nuc->raw.awr = 10.0;
    nuc->raw.xss_length = (int)(sizeof(values) / sizeof(values[0]));
    nuc->raw.nxs[0] = nuc->raw.xss_length;
    nuc->raw.nxs[5] = 1;
    nuc->raw.jxs[11] = 39;
    nuc->raw.jxs[12] = 1;
    nuc->raw.jxs[13] = 2;
    nuc->raw.jxs[14] = 3;
    nuc->raw.jxs[15] = 13;
    nuc->raw.jxs[16] = 28;
    nuc->raw.jxs[17] = 14;
    nuc->raw.jxs[18] = 15;
    nuc->raw.jxs[19] = 26;
    return nuc;
}

TEST(photon_production_decode_cleans_every_injected_allocation_failure) {
    size_t successful_at = SIZE_MAX;
    for (size_t fail_at = 0; fail_at < 128; fail_at++) {
        alea_nuc_nuclide_t* nuc = photon_production_allocation_fixture();
        ASSERT_NOT_NULL(nuc);
        nuc_alloc_failure_t failure = {fail_at, 0};
        alea_nuc_set_alloc_failure(fail_nuc_allocation, &failure);
        alea_error_t status =
            alea_nuc_decode_photon_production(nuc, &nuc->raw);
        alea_nuc_set_alloc_failure(NULL, NULL);
        if (status == ALEA_OK) {
            successful_at = fail_at;
            ASSERT_EQ(nuc->n_photon_productions, 1);
            ASSERT_NOT_NULL(nuc->photon_productions[0].angular);
            ASSERT_NOT_NULL(nuc->photon_productions[0].spectrum);
            alea_nuc_nuclide_free(nuc);
            break;
        }
        ASSERT_EQ(status, ALEA_ERR_OUT_OF_MEMORY);
        ASSERT_TRUE(failure.calls > fail_at);
        alea_nuc_nuclide_free(nuc);
    }
    ASSERT_NE(successful_at, SIZE_MAX);
    ASSERT_TRUE(successful_at >= 16);
}

TEST(photon_production_decodes_yield_and_cross_section_forms) {
    alea_nuc_nuclide_t* nuc = calloc(1, sizeof(*nuc));
    ASSERT_NOT_NULL(nuc);
    nuc->particle = ALEA_NUC_PARTICLE_NEUTRON;
    nuc->awr = 10.0;
    nuc->n_energies = 2;
    nuc->energy = malloc(2 * sizeof(double));
    nuc->sigma_total = malloc(2 * sizeof(double));
    nuc->sigma_elastic = calloc(2, sizeof(double));
    nuc->sigma_abs = calloc(2, sizeof(double));
    nuc->heating = calloc(2, sizeof(double));
    nuc->reactions = calloc(1, sizeof(*nuc->reactions));
    ASSERT_NOT_NULL(nuc->energy);
    ASSERT_NOT_NULL(nuc->sigma_total);
    ASSERT_NOT_NULL(nuc->sigma_elastic);
    ASSERT_NOT_NULL(nuc->sigma_abs);
    ASSERT_NOT_NULL(nuc->heating);
    ASSERT_NOT_NULL(nuc->reactions);
    nuc->energy[0] = 1.0; nuc->energy[1] = 2.0;
    nuc->sigma_total[0] = 2.0; nuc->sigma_total[1] = 4.0;
    nuc->sigma_abs[0] = 2.0; nuc->sigma_abs[1] = 4.0;
    nuc->n_reactions = 1;
    nuc->reactions[0].mt = 102;
    nuc->reactions[0].threshold_index = 1;
    nuc->reactions[0].n_energies = 2;
    nuc->reactions[0].xs = malloc(2 * sizeof(double));
    ASSERT_NOT_NULL(nuc->reactions[0].xs);
    nuc->reactions[0].xs[0] = 2.0;
    nuc->reactions[0].xs[1] = 4.0;

    double values[] = {
        102001, 102002, /* MTRP */
        1, 9,           /* LSIGP */
        12, 102, 0, 2, 1.0, 2.0, 0.5, 1.5,
        13, 1, 2, 1.0, 2.0,
        0, 0,           /* LANDP: isotropic */
        1, 1,           /* LDLWP */
        0, 2, 10, 0, 2, 1.0, 2.0, 1.0, 1.0,
        2, 0.9,         /* primary discrete-photon law */
        1, 102          /* YP */
    };
    nuc->raw.xss_length = (int)(sizeof(values) / sizeof(values[0]));
    nuc->raw.xss = malloc(sizeof(values));
    ASSERT_NOT_NULL(nuc->raw.xss);
    memcpy(nuc->raw.xss, values, sizeof(values));
    nuc->raw.nxs[0] = nuc->raw.xss_length;
    nuc->raw.nxs[5] = 2;
    nuc->raw.awr = 10.0;
    nuc->raw.jxs[12] = 1;
    nuc->raw.jxs[13] = 3;
    nuc->raw.jxs[14] = 5;
    nuc->raw.jxs[15] = 18;
    nuc->raw.jxs[17] = 20;
    nuc->raw.jxs[18] = 22;
    nuc->raw.jxs[19] = (int)(sizeof(values) / sizeof(values[0])) - 1;

    ASSERT_EQ(alea_nuc_decode_photon_production(nuc, &nuc->raw), ALEA_OK);
    ASSERT_EQ(nuc->n_photon_productions, 2);
    ASSERT_EQ(nuc->n_photon_yield_multipliers, 1);
    ASSERT_EQ(nuc->photon_yield_multipliers[0], 102);
    ASSERT_EQ(nuc->photon_productions[0].parent_mt, 102);
    ASSERT_EQ(nuc->photon_productions[0].mf, 12);
    ASSERT_NEAR(alea_nuc_photon_production_yield(
                    nuc, &nuc->photon_productions[0], 1.5), 1.0, 1e-14);
    ASSERT_EQ(nuc->photon_productions[1].mf, 13);
    ASSERT_TRUE(nuc->photon_productions[1].production_xs);
    ASSERT_NEAR(alea_nuc_photon_production_yield(
                    nuc, &nuc->photon_productions[1], 1.5), 0.5, 1e-14);
    ASSERT_NEAR(alea_nuc_xs_photon_production_total(nuc, 1.5),
                4.5, 1e-14);
    ASSERT_NOT_NULL(nuc->photon_productions[0].spectrum);
    ASSERT_EQ(nuc->photon_productions[0].spectrum->law,
              ALEA_NUC_ELAW_DISCRETE_PHOTON);
    double draws_data[] = {0.5};
    sequence_rng_t draws = {draws_data, 1, 0};
    double photon_energy = 0.0;
    ASSERT_EQ(alea_nuc_sample_energy_distribution(
                  nuc->photon_productions[0].spectrum, 1.1,
                  sequence_rng, &draws, &photon_energy), ALEA_OK);
    ASSERT_NEAR(photon_energy, 1.9, 1e-14);
    alea_nuc_particle_state_t incident = {
        ALEA_NUC_PARTICLE_NEUTRON, 1.1, {0.0, 0.0, 1.0}, 0.75, 2.0
    };
    double photon_draw_data[] = {0.1, 0.75, 0.0};
    sequence_rng_t photon_draws = {photon_draw_data, 3, 0};
    alea_nuc_particle_state_t photon = {0};
    ASSERT_EQ(alea_nuc_sample_photon_production(
                  nuc, &nuc->photon_productions[0], &incident,
                  sequence_rng, &photon_draws, &photon), ALEA_OK);
    ASSERT_EQ(photon.type, ALEA_NUC_PARTICLE_PHOTON);
    ASSERT_NEAR(photon.energy, 1.9, 1e-14);
    ASSERT_NEAR(photon.direction[2], 0.5, 1e-14);
    ASSERT_NEAR(photon.weight, 0.75, 1e-14);
    ASSERT_NEAR(photon.time, 2.0, 1e-14);

    incident.energy = 1.5;
    alea_nuc_mat_component_t component = {nuc, 0.1};
    alea_nuc_material_t material = {&component, 1, 1};
    alea_nuc_prepare_requirements_t requirements = {
        .required_capabilities = ALEA_NUC_CAP_RESTRICTED_NEUTRON |
                                 ALEA_NUC_CAP_PHOTON_PRODUCTION
    };
    alea_nuc_capability_report_t report;
    alea_nuc_prepared_material_t* prepared = NULL;
    ASSERT_EQ(alea_nuc_prepare_material(
                  &material, &requirements, &report, &prepared), ALEA_OK);
    alea_nuc_evaluation_t evaluation;
    ASSERT_EQ(alea_nuc_evaluate(prepared, &incident, &evaluation), ALEA_OK);
    alea_nuc_particle_state_t particles[2];
    alea_nuc_secondary_buffer_t buffer = {particles, 1, 0};
    double collision_data[] = {0.0, 0.0, 0.1, 0.75, 0.0, 0.9};
    sequence_rng_t collision_draws = {collision_data, 6, 0};
    alea_nuc_collision_result_t result;
    ASSERT_EQ(alea_nuc_collide_with_secondaries(
                  &evaluation, sequence_rng, &collision_draws, &buffer,
                  &result), ALEA_ERR_OUT_OF_MEMORY);
    ASSERT_EQ(collision_draws.position, 0);
    ASSERT_EQ(buffer.count, 0);
    buffer.capacity = 2;
    ASSERT_EQ(alea_nuc_collide_with_secondaries(
                  &evaluation, sequence_rng, &collision_draws, &buffer,
                  &result), ALEA_OK);
    ASSERT_EQ(result.mt, 102);
    ASSERT_EQ(result.outcome, ALEA_NUC_OUTCOME_REPLACED);
    ASSERT_EQ(result.n_emitted, 1);
    ASSERT_EQ(buffer.count, 1);
    ASSERT_EQ(buffer.particles[0].type, ALEA_NUC_PARTICLE_PHOTON);
    ASSERT_NEAR(buffer.particles[0].energy, 0.9 + 15.0 / 11.0, 1e-14);
    alea_nuc_prepared_material_free(prepared);
    alea_nuc_nuclide_free(nuc);
}

TEST(prepared_photon_production_accepts_mt3_aggregate_parent) {
    double energy[] = {1.0, 2.0};
    double total[] = {2.0, 4.0};
    double zero[] = {0.0, 0.0};
    double reaction_xs[] = {2.0, 4.0};
    double yield[] = {0.25, 0.25};
    alea_nuc_energy_dist_t spectrum = {0};
    spectrum.law = ALEA_NUC_ELAW_DISCRETE_PHOTON;
    spectrum.discrete_photon_primary = 0;
    spectrum.discrete_photon_energy = 0.75;
    spectrum.discrete_photon_awr = 10.0;
    alea_nuc_reaction_t reaction = {0};
    reaction.mt = 102;
    reaction.threshold_index = 1;
    reaction.n_energies = 2;
    reaction.xs = reaction_xs;
    alea_nuc_photon_production_t production = {0};
    production.mt = 3001;
    production.parent_mt = 3;
    production.mf = 12;
    production.n_energies = 2;
    production.energy = energy;
    production.values = yield;
    production.spectrum = &spectrum;
    alea_nuc_nuclide_t nuc = {0};
    nuc.particle = ALEA_NUC_PARTICLE_NEUTRON;
    nuc.awr = 10.0;
    nuc.n_energies = 2;
    nuc.energy = energy;
    nuc.sigma_total = total;
    nuc.sigma_elastic = zero;
    nuc.sigma_abs = total;
    nuc.heating = zero;
    nuc.n_reactions = 1;
    nuc.reactions = &reaction;
    nuc.n_photon_productions = 1;
    nuc.photon_productions = &production;

    alea_nuc_mat_component_t component = {&nuc, 0.1};
    alea_nuc_material_t material = {&component, 1, 1};
    alea_nuc_prepare_requirements_t requirements = {
        .required_capabilities = ALEA_NUC_CAP_RESTRICTED_NEUTRON |
                                 ALEA_NUC_CAP_PHOTON_PRODUCTION
    };
    alea_nuc_capability_report_t report;
    alea_nuc_prepared_material_t* prepared = NULL;
    ASSERT_EQ(alea_nuc_prepare_material(
                  &material, &requirements, &report, &prepared), ALEA_OK);
    alea_nuc_particle_state_t incident = {
        ALEA_NUC_PARTICLE_NEUTRON, 1.5, {0.0, 0.0, 1.0}, 1.0, 0.0
    };
    alea_nuc_evaluation_t evaluation;
    ASSERT_EQ(alea_nuc_evaluate(prepared, &incident, &evaluation), ALEA_OK);
    alea_nuc_particle_state_t particles[1];
    alea_nuc_secondary_buffer_t buffer = {particles, 1, 0};
    double draws_data[] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    sequence_rng_t draws = {draws_data, 6, 0};
    alea_nuc_collision_result_t result;
    ASSERT_EQ(alea_nuc_collide_with_secondaries(
                  &evaluation, sequence_rng, &draws, &buffer, &result),
              ALEA_OK);
    ASSERT_EQ(result.mt, 102);
    ASSERT_EQ(result.n_emitted, 1);
    ASSERT_EQ(buffer.count, 1);
    ASSERT_EQ(buffer.particles[0].type, ALEA_NUC_PARTICLE_PHOTON);
    ASSERT_NEAR(buffer.particles[0].energy, 0.75, 1e-14);

    nucdata_rng_t statistical_rng = {UINT64_C(0x489ad512a55d71f3)};
    int emitted = 0;
    const int samples = 50000;
    for (int i = 0; i < samples; i++) {
        buffer.count = 0;
        ASSERT_EQ(alea_nuc_collide_with_secondaries(
                      &evaluation, nucdata_rng, &statistical_rng, &buffer,
                      &result), ALEA_OK);
        ASSERT_EQ(result.mt, 102);
        ASSERT_TRUE(result.n_emitted == 0 || result.n_emitted == 1);
        ASSERT_EQ(buffer.count, result.n_emitted);
        emitted += result.n_emitted;
    }
    ASSERT_NEAR((double)emitted / samples, 0.25, 0.01);
    alea_nuc_prepared_material_free(prepared);
}

TEST(prepared_photon_production_conditions_partial_parent_on_aggregate) {
    double energy[] = {1.0, 2.0};
    double total[] = {4.0, 4.0};
    double zero[] = {0.0, 0.0};
    double aggregate_xs[] = {4.0, 4.0};
    double partial_xs[] = {1.0, 1.0};
    double yield[] = {1.0, 1.0};
    alea_nuc_energy_dist_t spectrum = {0};
    spectrum.law = ALEA_NUC_ELAW_DISCRETE_PHOTON;
    spectrum.discrete_photon_energy = 0.75;
    spectrum.discrete_photon_awr = 10.0;
    alea_nuc_reaction_t reactions[2] = {0};
    reactions[0].mt = 103;
    reactions[0].threshold_index = 1;
    reactions[0].n_energies = 2;
    reactions[0].xs = aggregate_xs;
    reactions[1].mt = 601;
    reactions[1].threshold_index = 1;
    reactions[1].n_energies = 2;
    reactions[1].xs = partial_xs;
    alea_nuc_photon_production_t production = {0};
    production.mt = 601001;
    production.parent_mt = 601;
    production.mf = 12;
    production.n_energies = 2;
    production.energy = energy;
    production.values = yield;
    production.spectrum = &spectrum;
    alea_nuc_nuclide_t nuc = {0};
    nuc.particle = ALEA_NUC_PARTICLE_NEUTRON;
    nuc.awr = 10.0;
    nuc.n_energies = 2;
    nuc.energy = energy;
    nuc.sigma_total = total;
    nuc.sigma_elastic = zero;
    nuc.sigma_abs = total;
    nuc.heating = zero;
    nuc.n_reactions = 2;
    nuc.reactions = reactions;
    nuc.n_photon_productions = 1;
    nuc.photon_productions = &production;

    alea_nuc_mat_component_t component = {&nuc, 0.1};
    alea_nuc_material_t material = {&component, 1, 1};
    alea_nuc_prepare_requirements_t requirements = {
        .required_capabilities = ALEA_NUC_CAP_RESTRICTED_NEUTRON |
                                 ALEA_NUC_CAP_PHOTON_PRODUCTION
    };
    alea_nuc_capability_report_t report;
    alea_nuc_prepared_material_t* prepared = NULL;
    ASSERT_EQ(alea_nuc_prepare_material(
                  &material, &requirements, &report, &prepared), ALEA_OK);
    alea_nuc_particle_state_t incident = {
        ALEA_NUC_PARTICLE_NEUTRON, 1.5, {0.0, 0.0, 1.0}, 1.0, 0.0
    };
    alea_nuc_evaluation_t evaluation;
    ASSERT_EQ(alea_nuc_evaluate(prepared, &incident, &evaluation), ALEA_OK);
    alea_nuc_particle_state_t particle;
    alea_nuc_secondary_buffer_t buffer = {&particle, 1, 0};
    alea_nuc_collision_result_t result;
    nucdata_rng_t rng = {UINT64_C(0x68c327aab9417f25)};
    int emitted = 0;
    const int samples = 50000;
    for (int i = 0; i < samples; i++) {
        buffer.count = 0;
        ASSERT_EQ(alea_nuc_collide_with_secondaries(
                      &evaluation, nucdata_rng, &rng, &buffer, &result),
                  ALEA_OK);
        ASSERT_EQ(result.mt, 103);
        ASSERT_TRUE(result.n_emitted == 0 || result.n_emitted == 1);
        ASSERT_EQ(buffer.count, result.n_emitted);
        emitted += result.n_emitted;
    }
    ASSERT_NEAR((double)emitted / samples, 0.25, 0.01);
    alea_nuc_prepared_material_free(prepared);

    /* Production cross sections for omitted composite parents are normalized
     * by the inclusive MT=5 event that carries the neutron channel. */
    reactions[0].mt = 5;
    reactions[0].ty = 1;
    reactions[0].energy = &spectrum;
    partial_xs[0] = 0.0;
    partial_xs[1] = 0.0;
    nuc.sigma_abs = zero;
    production.mt = 28001;
    production.parent_mt = 28;
    production.mf = 13;
    production.production_xs = true;
    production.threshold_index = 1;
    production.energy = NULL;
    ASSERT_NEAR(alea_nuc_photon_production_yield(
                    &nuc, &production, 1.5), 0.0, 1e-14);
    prepared = NULL;
    ASSERT_EQ(alea_nuc_prepare_material(
                  &material, &requirements, &report, &prepared), ALEA_OK);
    ASSERT_EQ(alea_nuc_evaluate(prepared, &incident, &evaluation), ALEA_OK);
    alea_nuc_particle_state_t emitted_particles[2];
    buffer.particles = emitted_particles;
    buffer.capacity = 2;
    rng.state = UINT64_C(0x7e81c559fb2a416d);
    emitted = 0;
    for (int i = 0; i < samples; i++) {
        buffer.count = 0;
        ASSERT_EQ(alea_nuc_collide_with_secondaries(
                      &evaluation, nucdata_rng, &rng, &buffer, &result),
                  ALEA_OK);
        ASSERT_EQ(result.mt, 5);
        ASSERT_TRUE(result.n_emitted == 1 || result.n_emitted == 2);
        emitted += result.n_emitted == 2;
    }
    ASSERT_NEAR((double)emitted / samples, 0.25, 0.01);
    alea_nuc_prepared_material_free(prepared);
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
