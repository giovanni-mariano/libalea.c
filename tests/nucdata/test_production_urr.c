// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/** @file test_production_urr.c
 *  @brief Production U-238 unresolved-resonance transport benchmark
 */

#include "alea_nucdata.h"
#include "alea_test.h"

#include <math.h>
#include <stdlib.h>

#define URANIUM_ZAID "92238.00c"

static alea_nuc_xsdir_t* xsdir;
static alea_nuc_nuclide_t* uranium;
static int teardown_registered;

typedef struct {
    double value;
    size_t calls;
} fixed_rng_t;

static double fixed_uniform(void* context) {
    fixed_rng_t* rng = context;
    rng->calls++;
    return rng->value;
}

static void teardown(void) {
    alea_nuc_xsdir_free(xsdir);
}

static void setup(void) {
    if (xsdir) return;
    const char* path = getenv("ALEA_ENDFB80_XSDIR");
    if (!path || !path[0]) return;
    xsdir = alea_nuc_xsdir_load(path);
    if (xsdir) uranium = alea_nuc_xsdir_get_nuclide(xsdir, URANIUM_ZAID);
    if (!teardown_registered) {
        atexit(teardown);
        teardown_registered = 1;
    }
}

static int production_data_requested(void) {
    const char* path = getenv("ALEA_ENDFB80_XSDIR");
    return path && path[0];
}

static alea_nuc_prepared_material_t* prepare(double number_density,
                                              alea_nuc_material_t** material) {
    *material = alea_nuc_material_create();
    if (!*material ||
        alea_nuc_material_add(*material, uranium, number_density) != ALEA_OK)
        return NULL;
    alea_nuc_prepare_requirements_t requirements = {
        .required_capabilities = ALEA_NUC_CAP_CONTINUOUS_NEUTRON |
                                 ALEA_NUC_CAP_URR
    };
    alea_nuc_capability_report_t report;
    alea_nuc_prepared_material_t* prepared = NULL;
    if (alea_nuc_prepare_material(*material, &requirements, &report,
                                  &prepared) != ALEA_OK)
        return NULL;
    return prepared;
}

TEST(production_u238_urr_metadata_and_sampled_rates) {
    setup();
    if (!production_data_requested()) SKIP("ALEA_ENDFB80_XSDIR is not set");
    ASSERT_NOT_NULL(xsdir);
    ASSERT_NOT_NULL(uranium);
    ASSERT_NOT_NULL(uranium->urr);
    ASSERT_EQ(uranium->urr->n_energies, 83);
    ASSERT_EQ(uranium->urr->n_bands, 16);
    ASSERT_EQ(uranium->urr->interp, 2);
    ASSERT_EQ(uranium->urr->inelastic_flag, 51);
    ASSERT_EQ(uranium->urr->absorption_flag, 0);
    ASSERT_TRUE(uranium->urr->multiply_smooth);
    ASSERT_NEAR(uranium->urr->energy[0], 0.02000001, 1e-14);
    ASSERT_NEAR(uranium->urr->energy[82], 0.1490086, 1e-13);

    const double energy = 0.05;
    const double xi = 0.5;
    double factors[5];
    ASSERT_TRUE(alea_nuc_urr_factors(uranium, energy, xi, factors));
    const double expected[] = {
        0.9666967, 0.9650968, 1.0, 1.022523, 1.0
    };
    for (int i = 0; i < 5; i++) ASSERT_NEAR(factors[i], expected[i], 1e-13);

    alea_nuc_material_t* material = NULL;
    alea_nuc_prepared_material_t* prepared = prepare(0.01, &material);
    ASSERT_NOT_NULL(prepared);
    alea_nuc_particle_state_t incident = {
        ALEA_NUC_PARTICLE_NEUTRON, energy, {0.0, 0.0, 1.0}, 1.0, 0.0
    };
    alea_nuc_urr_sample_t sample;
    alea_nuc_evaluation_workspace_t workspace = {&sample, 1};
    alea_nuc_evaluation_t evaluation;
    fixed_rng_t urr_rng = {xi, 0};
    ASSERT_EQ(alea_nuc_evaluate_urr(prepared, &incident, fixed_uniform,
                                    &urr_rng, &workspace, &evaluation), ALEA_OK);
    ASSERT_EQ(urr_rng.calls, 1);
    ASSERT_TRUE(sample.active);
    ASSERT_NEAR(evaluation.macro_total, 0.125996254677189, 2e-14);
    ASSERT_NEAR(evaluation.macro_elastic, 0.121930329712, 2e-14);
    ASSERT_NEAR(evaluation.macro_absorption, 0.0033072995180265974, 2e-14);
    ASSERT_NEAR(evaluation.macro_neutron_emission,
                0.00075862544716241004, 2e-14);
    ASSERT_NEAR(evaluation.macro_total,
                evaluation.macro_elastic + evaluation.macro_absorption +
                evaluation.macro_neutron_emission, 2e-14);

    alea_nuc_prepared_material_free(prepared);
    alea_nuc_material_destroy(material);
}

TEST(production_u238_urr_transmission_across_boundaries) {
    setup();
    if (!production_data_requested()) SKIP("ALEA_ENDFB80_XSDIR is not set");
    ASSERT_NOT_NULL(uranium);
    alea_nuc_material_t* material = NULL;
    alea_nuc_prepared_material_t* prepared = prepare(0.01, &material);
    ASSERT_NOT_NULL(prepared);
    alea_nuc_particle_state_t incident = {
        ALEA_NUC_PARTICLE_NEUTRON, 0.05, {0.0, 0.0, 1.0}, 1.0, 0.0
    };
    alea_nuc_urr_sample_t sample;
    alea_nuc_evaluation_workspace_t workspace = {&sample, 1};
    alea_nuc_evaluation_t evaluation;
    fixed_rng_t urr_rng = {0.5, 0};
    ASSERT_EQ(alea_nuc_evaluate_urr(prepared, &incident, fixed_uniform,
                                    &urr_rng, &workspace, &evaluation), ALEA_OK);

    const double first_length = 3.0;
    const double second_length = 4.0;
    const int histories = 131072;
    int transmitted = 0;
    for (int i = 0; i < histories; i++) {
        alea_nuc_rng_t flight_rng;
        ASSERT_EQ(alea_nuc_rng_init(&flight_rng, 0x238u, (uint64_t)i, 0, 0,
                                    ALEA_NUC_RNG_FLIGHT), ALEA_OK);
        double distance;
        ASSERT_EQ(alea_nuc_sample_flight(&evaluation, alea_nuc_rng_uniform,
                                         &flight_rng, &distance), ALEA_OK);
        if (distance < first_length) continue;
        ASSERT_EQ(alea_nuc_sample_flight(&evaluation, alea_nuc_rng_uniform,
                                         &flight_rng, &distance), ALEA_OK);
        if (distance >= second_length) transmitted++;
    }
    double expected = exp(-evaluation.macro_total *
                          (first_length + second_length));
    ASSERT_NEAR((double)transmitted / histories, expected, 0.006);
    ASSERT_EQ(urr_rng.calls, 1);
    const double expected_factors[] = {
        0.9666967, 0.9650968, 1.0, 1.022523, 1.0
    };
    for (int i = 0; i < 5; i++)
        ASSERT_NEAR(sample.factors[i], expected_factors[i], 1e-13);

    alea_nuc_material_t* dense_material = NULL;
    alea_nuc_prepared_material_t* dense_prepared =
        prepare(0.025, &dense_material);
    ASSERT_NOT_NULL(dense_prepared);
    alea_nuc_urr_sample_t dense_sample;
    alea_nuc_evaluation_workspace_t dense_workspace = {&dense_sample, 1};
    alea_nuc_evaluation_t dense_evaluation;
    fixed_rng_t dense_rng = {0.5, 0};
    ASSERT_EQ(alea_nuc_evaluate_urr(dense_prepared, &incident, fixed_uniform,
                                    &dense_rng, &dense_workspace,
                                    &dense_evaluation), ALEA_OK);
    ASSERT_EQ(dense_rng.calls, 1);
    ASSERT_NEAR(dense_evaluation.macro_total,
                2.5 * evaluation.macro_total, 3e-14);
    for (int i = 0; i < 5; i++)
        ASSERT_NEAR(dense_sample.factors[i], sample.factors[i], 1e-14);

    alea_nuc_prepared_material_free(dense_prepared);
    alea_nuc_material_destroy(dense_material);
    alea_nuc_prepared_material_free(prepared);
    alea_nuc_material_destroy(material);
}

TEST_MAIN()
