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
#define URANIUM_HOT_ZAID "92238.01c"

static alea_nuc_xsdir_t* xsdir;
static alea_nuc_nuclide_t* uranium;
static alea_nuc_nuclide_t* uranium_hot;
static int teardown_registered;

typedef struct {
    double value;
    size_t calls;
} fixed_rng_t;

typedef struct {
    double target;
    size_t calls;
} elastic_collision_rng_t;

static double fixed_uniform(void* context) {
    fixed_rng_t* rng = context;
    rng->calls++;
    return rng->value;
}

static double elastic_collision_uniform(void* context) {
    elastic_collision_rng_t* rng = context;
    size_t call = rng->calls++;
    if (call == 0) return rng->target;
    if (call == 1) return 0.0;
    return 0.5;
}

static void teardown(void) {
    alea_nuc_xsdir_free(xsdir);
}

static void setup(void) {
    if (xsdir) return;
    const char* path = getenv("ALEA_ENDFB80_XSDIR");
    if (!path || !path[0]) return;
    xsdir = alea_nuc_xsdir_load(path);
    if (xsdir) {
        uranium = alea_nuc_xsdir_get_nuclide(xsdir, URANIUM_ZAID);
        uranium_hot = alea_nuc_xsdir_get_nuclide(xsdir, URANIUM_HOT_ZAID);
    }
    if (!teardown_registered) {
        atexit(teardown);
        teardown_registered = 1;
    }
}

static int production_data_requested(void) {
    const char* path = getenv("ALEA_ENDFB80_XSDIR");
    return path && path[0];
}

static alea_nuc_prepared_material_t* prepare_nuclide(
    alea_nuc_nuclide_t* nuclide, double number_density,
    alea_nuc_material_t** material) {
    *material = alea_nuc_material_create();
    if (!*material ||
        alea_nuc_material_add(*material, nuclide, number_density) != ALEA_OK)
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

static alea_nuc_prepared_material_t* prepare(double number_density,
                                              alea_nuc_material_t** material) {
    return prepare_nuclide(uranium, number_density, material);
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

TEST(production_u238_temperature_mix_correlates_urr_tables) {
    setup();
    if (!production_data_requested()) SKIP("ALEA_ENDFB80_XSDIR is not set");
    ASSERT_NOT_NULL(uranium);
    ASSERT_NOT_NULL(uranium_hot);
    ASSERT_NOT_NULL(uranium->urr);
    ASSERT_NOT_NULL(uranium_hot->urr);
    ASSERT_NEAR(uranium->temperature, 2.5301e-8, 1e-16);
    ASSERT_NEAR(uranium_hot->temperature, 5.1704e-8, 1e-16);

    const double fraction = 0.5;
    const double density = 0.01;
    alea_nuc_material_t* material = alea_nuc_material_create();
    ASSERT_NOT_NULL(material);
    ASSERT_EQ(alea_nuc_material_add_temperature_mix(
                  material, uranium, uranium_hot, fraction, density), ALEA_OK);
    alea_nuc_prepare_requirements_t requirements = {
        .required_capabilities = ALEA_NUC_CAP_CONTINUOUS_NEUTRON |
                                 ALEA_NUC_CAP_URR
    };
    alea_nuc_capability_report_t report;
    alea_nuc_prepared_material_t* prepared = NULL;
    ASSERT_EQ(alea_nuc_prepare_material(material, &requirements, &report,
                                        &prepared), ALEA_OK);
    alea_nuc_particle_state_t incident = {
        ALEA_NUC_PARTICLE_NEUTRON, 0.05, {0.0, 0.0, 1.0}, 1.0, 0.0
    };
    alea_nuc_urr_sample_t samples[2];
    alea_nuc_evaluation_workspace_t workspace = {samples, 2};
    alea_nuc_evaluation_t evaluation;
    fixed_rng_t rng = {0.5, 0};
    ASSERT_EQ(alea_nuc_evaluate_urr(prepared, &incident, fixed_uniform, &rng,
                                    &workspace, &evaluation), ALEA_OK);
    ASSERT_EQ(rng.calls, 1);
    double lower_factors[5], upper_factors[5];
    ASSERT_TRUE(alea_nuc_urr_factors(
        uranium, incident.energy, 0.5, lower_factors));
    ASSERT_TRUE(alea_nuc_urr_factors(
        uranium_hot, incident.energy, 0.5, upper_factors));
    for (int i = 0; i < 5; i++) {
        ASSERT_NEAR(samples[0].factors[i], lower_factors[i], 1e-14);
        ASSERT_NEAR(samples[1].factors[i], upper_factors[i], 1e-14);
    }

    alea_nuc_material_t *lower_material = NULL, *upper_material = NULL;
    alea_nuc_prepared_material_t* lower_prepared =
        prepare_nuclide(uranium, density, &lower_material);
    alea_nuc_prepared_material_t* upper_prepared =
        prepare_nuclide(uranium_hot, density, &upper_material);
    ASSERT_NOT_NULL(lower_prepared);
    ASSERT_NOT_NULL(upper_prepared);
    alea_nuc_urr_sample_t lower_sample, upper_sample;
    alea_nuc_evaluation_workspace_t lower_workspace = {&lower_sample, 1};
    alea_nuc_evaluation_workspace_t upper_workspace = {&upper_sample, 1};
    alea_nuc_evaluation_t lower_evaluation, upper_evaluation;
    fixed_rng_t lower_rng = {0.5, 0}, upper_rng = {0.5, 0};
    ASSERT_EQ(alea_nuc_evaluate_urr(
                  lower_prepared, &incident, fixed_uniform, &lower_rng,
                  &lower_workspace, &lower_evaluation), ALEA_OK);
    ASSERT_EQ(alea_nuc_evaluate_urr(
                  upper_prepared, &incident, fixed_uniform, &upper_rng,
                  &upper_workspace, &upper_evaluation), ALEA_OK);
    ASSERT_NEAR(evaluation.macro_total,
                (1.0 - fraction) * lower_evaluation.macro_total +
                fraction * upper_evaluation.macro_total, 3e-14);
    ASSERT_NEAR(evaluation.macro_elastic,
                (1.0 - fraction) * lower_evaluation.macro_elastic +
                fraction * upper_evaluation.macro_elastic, 3e-14);
    ASSERT_NEAR(evaluation.macro_absorption,
                (1.0 - fraction) * lower_evaluation.macro_absorption +
                fraction * upper_evaluation.macro_absorption, 3e-14);
    ASSERT_NEAR(evaluation.macro_neutron_emission,
                (1.0 - fraction) * lower_evaluation.macro_neutron_emission +
                fraction * upper_evaluation.macro_neutron_emission, 3e-14);

    double expected_upper = fraction * upper_evaluation.macro_total /
                            evaluation.macro_total;
    alea_nuc_rng_t selection_rng;
    ASSERT_EQ(alea_nuc_rng_init(&selection_rng, 0x23801u, 238, 0, 17,
                                ALEA_NUC_RNG_COLLISION), ALEA_OK);
    const int histories = 32768;
    int selected_upper = 0;
    for (int i = 0; i < histories; i++) {
        elastic_collision_rng_t collision_rng = {
            alea_nuc_rng_uniform(&selection_rng), 0
        };
        alea_nuc_collision_result_t result;
        ASSERT_EQ(alea_nuc_collide(&evaluation, elastic_collision_uniform,
                                   &collision_rng, &result), ALEA_OK);
        ASSERT_TRUE(result.component_index == 0 || result.component_index == 1);
        ASSERT_EQ(result.mt, 2);
        if (result.component_index == 1) selected_upper++;
    }
    ASSERT_NEAR((double)selected_upper / histories, expected_upper, 0.01);

    alea_nuc_prepared_material_free(lower_prepared);
    alea_nuc_prepared_material_free(upper_prepared);
    alea_nuc_material_destroy(lower_material);
    alea_nuc_material_destroy(upper_material);
    alea_nuc_prepared_material_free(prepared);
    alea_nuc_material_destroy(material);
}

TEST_MAIN()
