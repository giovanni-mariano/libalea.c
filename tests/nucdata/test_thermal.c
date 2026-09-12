// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/**
 * @file test_thermal.c
 * @brief ENDF/B-VII.1 thermal ACE checks against OpenMC 0.15.3 references
 */

#include "alea_nucdata.h"
#include "alea_test.h"

#define TSL_XSDIR "endfb71-tsl/xsdir"
#define NJOY74_XSDIR "njoy-test74/xsdir"
#define NJOY25_XSDIR "njoy-test25/xsdir"

static alea_nuc_thermal_t* load_table(const char* name,
                                      alea_nuc_xsdir_t** xsdir) {
    *xsdir = alea_nuc_xsdir_load(TSL_XSDIR);
    if (!*xsdir) return NULL;
    return alea_nuc_load_thermal(*xsdir, name);
}

static double sequence_random(void* context) {
    double** values = context;
    return *(*values)++;
}

typedef struct {
    alea_nuc_rng_t core;
    int force_inelastic;
} thermal_rng_t;

static double thermal_random(void* context) {
    thermal_rng_t* rng = context;
    if (rng->force_inelastic) {
        rng->force_inelastic = 0;
        return 0.0;
    }
    return alea_nuc_rng_uniform(&rng->core);
}

static double skewed_outgoing_probability(int index, int count);

TEST(endfb71_all_thermal_tables_decode) {
    static const char* names[] = {
        "Al.71t", "BeBeO.71t", "Be.71t", "Benz.71t", "DD2O.71t",
        "Fe.71t", "Graph.71t", "HCH2.71t", "HH2O.71t", "HZrH.71t",
        "lCH4.71t", "OBeO.71t", "orthoD.71t", "orthoH.71t",
        "OUO2.71t", "paraD.71t", "paraH.71t", "sCH4.71t",
        "UUO2.71t", "ZrZrH.71t"
    };
    alea_nuc_xsdir_t* xsdir = alea_nuc_xsdir_load(TSL_XSDIR);
    if (!xsdir) SKIP("ENDF/B-VII.1 thermal data not downloaded");
    for (size_t i = 0; i < sizeof(names) / sizeof(names[0]); i++) {
        alea_nuc_thermal_t* thermal = alea_nuc_load_thermal(xsdir, names[i]);
        ASSERT_NOT_NULL(thermal);
        ASSERT_TRUE(thermal->n_inelastic_energies > 1);
        alea_nuc_thermal_free(thermal);
    }
    alea_nuc_xsdir_free(xsdir);
}

TEST(endfb71_light_water_matches_openmc_reference) {
    alea_nuc_xsdir_t* xsdir = NULL;
    alea_nuc_thermal_t* thermal = load_table("HH2O.71t", &xsdir);
    if (!xsdir) SKIP("ENDF/B-VII.1 thermal data not downloaded");
    ASSERT_NOT_NULL(thermal);
    ASSERT_EQ(thermal->n_inelastic_energies, 108);
    ASSERT_EQ(thermal->n_inelastic_outgoing, 64);
    ASSERT_EQ(thermal->n_inelastic_cosines, 16);
    ASSERT_TRUE(thermal->inelastic_skewed);
    ASSERT_EQ(thermal->elastic_mode, ALEA_NUC_THERMAL_ELASTIC_NONE);
    ASSERT_NEAR(thermal->inelastic_energy[107], 4.46e-6, 1e-18);
    ASSERT_NEAR(alea_nuc_thermal_xs_inelastic(thermal, 2.53e-8),
                52.1483, 1e-10);
    size_t record = (size_t)50 * 64 + 20;
    ASSERT_NEAR(thermal->inelastic_energy_out[record],
                2.81524646514e-8, 1e-20);
    ASSERT_NEAR(thermal->inelastic_mu[record * 16 + 7],
                0.478054098149, 1e-12);
    alea_nuc_thermal_free(thermal);
    alea_nuc_xsdir_free(xsdir);
}

TEST(endfb71_graphite_matches_openmc_reference) {
    alea_nuc_xsdir_t* xsdir = NULL;
    alea_nuc_thermal_t* thermal = load_table("Graph.71t", &xsdir);
    if (!xsdir) SKIP("ENDF/B-VII.1 thermal data not downloaded");
    ASSERT_NOT_NULL(thermal);
    ASSERT_EQ(thermal->elastic_mode, ALEA_NUC_THERMAL_ELASTIC_COHERENT);
    ASSERT_EQ(thermal->n_coherent_edges, 166);
    ASSERT_NEAR(thermal->coherent_edge[0], 1.822197e-9, 1e-21);
    ASSERT_NEAR(thermal->coherent_factor[20], 1.29051403659e-7, 1e-19);
    ASSERT_NEAR(alea_nuc_thermal_xs_inelastic(thermal, 2.53e-8),
                0.3308013, 1e-10);
    ASSERT_NEAR(alea_nuc_thermal_xs_elastic(thermal, 2.53e-8),
                4.710462146561265, 1e-11);
    size_t record = (size_t)50 * 64 + 20;
    ASSERT_NEAR(thermal->inelastic_energy_out[record],
                2.87251692868e-8, 1e-20);
    ASSERT_NEAR(thermal->inelastic_mu[record * 16 + 7],
                -0.431445115511, 1e-12);
    alea_nuc_thermal_free(thermal);
    alea_nuc_xsdir_free(xsdir);
}

TEST(endfb71_graphite_applies_to_natural_carbon_isotopes) {
    alea_nuc_xsdir_t* xsdir = NULL;
    alea_nuc_thermal_t* thermal = load_table("Graph.71t", &xsdir);
    if (!xsdir) SKIP("ENDF/B-VII.1 thermal data not downloaded");
    ASSERT_NOT_NULL(thermal);

    double energy[] = {1.0e-11, 1.0e-5};
    double one[] = {1.0, 1.0}, zero[] = {0.0, 0.0};
    alea_nuc_nuclide_t carbon13 = {0};
    carbon13.Z = 6;
    carbon13.A = 13;
    carbon13.particle = ALEA_NUC_PARTICLE_NEUTRON;
    carbon13.awr = 13.0;
    carbon13.temperature = thermal->temperature;
    carbon13.n_energies = 2;
    carbon13.energy = energy;
    carbon13.sigma_total = one;
    carbon13.sigma_elastic = one;
    carbon13.sigma_abs = zero;
    alea_nuc_mat_component_t component = {&carbon13, 0.1};
    alea_nuc_material_t material = {&component, 1, 1, NULL};
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
    alea_nuc_prepared_material_free(prepared);
    alea_nuc_thermal_free(thermal);
    alea_nuc_xsdir_free(xsdir);
}

TEST(endfb71_graphite_collision_statistics_match_direct_integration) {
    alea_nuc_xsdir_t* xsdir = NULL;
    alea_nuc_thermal_t* thermal = load_table("Graph.71t", &xsdir);
    if (!xsdir) SKIP("ENDF/B-VII.1 thermal data not downloaded");
    ASSERT_NOT_NULL(thermal);
    ASSERT_TRUE(thermal->inelastic_skewed);
    ASSERT_EQ(thermal->elastic_mode, ALEA_NUC_THERMAL_ELASTIC_COHERENT);

    const int incident_index = 50;
    const double incident_energy = thermal->inelastic_energy[incident_index];
    double expected_energy = 0.0, expected_inelastic_mu = 0.0;
    for (int j = 0; j < thermal->n_inelastic_outgoing; j++) {
        double probability = skewed_outgoing_probability(
            j, thermal->n_inelastic_outgoing);
        size_t record = (size_t)incident_index *
            thermal->n_inelastic_outgoing + (size_t)j;
        expected_energy += probability * thermal->inelastic_energy_out[record];
        double mu = 0.0;
        for (int k = 0; k < thermal->n_inelastic_cosines; k++)
            mu += thermal->inelastic_mu[
                record * thermal->n_inelastic_cosines + (size_t)k];
        expected_inelastic_mu += probability * mu /
            thermal->n_inelastic_cosines;
    }

    int last_edge = -1;
    for (int i = 0; i < thermal->n_coherent_edges; i++) {
        if (thermal->coherent_edge[i] > incident_energy) break;
        last_edge = i;
    }
    ASSERT_TRUE(last_edge >= 0);
    double expected_coherent_mu = 0.0;
    double previous_factor = 0.0;
    for (int i = 0; i <= last_edge; i++) {
        double increment = thermal->coherent_factor[i] - previous_factor;
        expected_coherent_mu += increment *
            (1.0 - 2.0 * thermal->coherent_edge[i] / incident_energy);
        previous_factor = thermal->coherent_factor[i];
    }
    expected_coherent_mu /= thermal->coherent_factor[last_edge];

    double inelastic_xs = alea_nuc_thermal_xs_inelastic(
        thermal, incident_energy);
    double coherent_xs = alea_nuc_thermal_xs_elastic(
        thermal, incident_energy);
    double expected_inelastic_fraction =
        inelastic_xs / (inelastic_xs + coherent_xs);
    alea_nuc_particle_state_t incident = {
        ALEA_NUC_PARTICLE_NEUTRON, incident_energy,
        {0.0, 0.0, 1.0}, 1.0, 0.0
    };
    alea_nuc_rng_t rng;
    ASSERT_EQ(alea_nuc_rng_init(&rng, 0x671u, 6, 0, 71,
                               ALEA_NUC_RNG_COLLISION), ALEA_OK);
    const int samples = 131072;
    int inelastic = 0, coherent = 0;
    double energy_sum = 0.0, inelastic_mu_sum = 0.0, coherent_mu_sum = 0.0;
    for (int i = 0; i < samples; i++) {
        alea_nuc_collision_result_t result;
        ASSERT_EQ(alea_nuc_sample_thermal_collision(
                      thermal, &incident, alea_nuc_rng_uniform, &rng,
                      &result), ALEA_OK);
        if (result.mt == 4) {
            inelastic++;
            energy_sum += result.outgoing.energy;
            inelastic_mu_sum += result.mu_lab;
        } else {
            ASSERT_EQ(result.mt, 2);
            coherent++;
            ASSERT_NEAR(result.outgoing.energy, incident_energy, 1e-20);
            coherent_mu_sum += result.mu_lab;
        }
    }
    ASSERT_TRUE(inelastic > 0 && coherent > 0);
    ASSERT_NEAR((double)inelastic / samples, expected_inelastic_fraction,
                0.004);
    ASSERT_NEAR(energy_sum / inelastic, expected_energy, 5e-10);
    ASSERT_NEAR(inelastic_mu_sum / inelastic, expected_inelastic_mu, 0.012);
    ASSERT_NEAR(coherent_mu_sum / coherent, expected_coherent_mu, 0.008);

    alea_nuc_thermal_free(thermal);
    alea_nuc_xsdir_free(xsdir);
}

TEST(njoy74_continuous_thermal_matches_openmc_reference) {
    alea_nuc_xsdir_t* xsdir = alea_nuc_xsdir_load(NJOY74_XSDIR);
    if (!xsdir) {
        if (getenv("ALEA_NUCDATA_REQUIRED")) ASSERT_NOT_NULL(xsdir);
        SKIP("NJOY2016 test 74 data not downloaded");
    }
    alea_nuc_thermal_t* thermal = alea_nuc_load_thermal(xsdir, "hzrh.10t");
    ASSERT_NOT_NULL(thermal);
    ASSERT_TRUE(thermal->inelastic_continuous);
    ASSERT_EQ(thermal->n_inelastic_energies, 91);
    ASSERT_EQ(thermal->n_inelastic_cosines, 20);
    ASSERT_EQ(thermal->n_inelastic_outgoing_total, 99806);
    ASSERT_EQ(thermal->elastic_mode, ALEA_NUC_THERMAL_ELASTIC_INCOHERENT);
    ASSERT_EQ(thermal->n_incoherent_energies, 91);
    ASSERT_NEAR(alea_nuc_thermal_xs_inelastic(thermal, 2.5507e-8),
                2.8609020673826664, 1e-12);
    ASSERT_NEAR(alea_nuc_thermal_xs_elastic(thermal, 2.5507e-8),
                54.85262486393934, 1e-11);

    int begin = thermal->inelastic_outgoing_offset[45];
    ASSERT_EQ(thermal->inelastic_outgoing_offset[46] - begin, 903);
    ASSERT_NEAR(thermal->inelastic_energy_out[begin + 3],
                1.33669e-11, 1e-23);
    ASSERT_NEAR(thermal->inelastic_pdf[begin + 3],
                409609.995797, 1e-6);
    ASSERT_NEAR(thermal->inelastic_cdf[begin + 3],
                3.59428830443e-6, 1e-18);
    ASSERT_NEAR(thermal->inelastic_mu[(size_t)(begin + 3) * 20 + 5],
                -0.467951116, 1e-12);

    alea_nuc_particle_state_t incident = {
        ALEA_NUC_PARTICLE_NEUTRON, thermal->inelastic_energy[45],
        {0.0, 0.0, 1.0}, 1.0, 0.0
    };
    double draw_values[] = {0.0, 0.25, 0.3, 0.7, 0.0};
    double* draws = draw_values;
    alea_nuc_collision_result_t result;
    ASSERT_EQ(alea_nuc_sample_thermal_collision(
                  thermal, &incident, sequence_random, &draws, &result),
              ALEA_OK);
    ASSERT_NEAR(result.outgoing.energy, 1.8083420334178604e-8, 1e-20);
    ASSERT_NEAR(result.mu_lab, -0.5662829328095581, 1e-12);
    alea_nuc_thermal_free(thermal);
    alea_nuc_xsdir_free(xsdir);
}

TEST(njoy74_continuous_thermal_spectrum_statistics) {
    alea_nuc_xsdir_t* xsdir = alea_nuc_xsdir_load(NJOY74_XSDIR);
    if (!xsdir) {
        if (getenv("ALEA_NUCDATA_REQUIRED")) ASSERT_NOT_NULL(xsdir);
        SKIP("NJOY2016 test 74 data not downloaded");
    }
    alea_nuc_thermal_t* thermal = alea_nuc_load_thermal(xsdir, "hzrh.10t");
    ASSERT_NOT_NULL(thermal);

    const int incident_index = 45;
    alea_nuc_particle_state_t incident = {
        ALEA_NUC_PARTICLE_NEUTRON,
        thermal->inelastic_energy[incident_index],
        {0.0, 0.0, 1.0}, 1.0, 0.0
    };
    thermal_rng_t rng;
    ASSERT_EQ(alea_nuc_rng_init(&rng.core, 0x7416u, 74, 0, 2,
                               ALEA_NUC_RNG_COLLISION), ALEA_OK);
    const int samples = 131072;
    int upscattered = 0;
    double energy_sum = 0.0;
    double mu_sum = 0.0;
    for (int i = 0; i < samples; i++) {
        alea_nuc_collision_result_t result;
        rng.force_inelastic = 1;
        ASSERT_EQ(alea_nuc_sample_thermal_collision(
                      thermal, &incident, thermal_random, &rng, &result),
                  ALEA_OK);
        ASSERT_EQ(result.mt, 4);
        energy_sum += result.outgoing.energy;
        mu_sum += result.mu_lab;
        if (result.outgoing.energy > incident.energy) upscattered++;
    }

    /* Independently integrated from the pinned outgoing PDF and cosine sets. */
    ASSERT_NEAR(energy_sum / samples, 4.415954963033656e-8, 8e-10);
    ASSERT_NEAR(mu_sum / samples, -0.1969222486289254, 0.009);
    ASSERT_NEAR((double)upscattered / samples,
                0.7039972202844251, 0.008);
    alea_nuc_thermal_free(thermal);
    alea_nuc_xsdir_free(xsdir);
}

TEST(njoy74_prepared_moderator_channel_statistics) {
    alea_nuc_xsdir_t* xsdir = alea_nuc_xsdir_load(NJOY74_XSDIR);
    if (!xsdir) {
        if (getenv("ALEA_NUCDATA_REQUIRED")) ASSERT_NOT_NULL(xsdir);
        SKIP("NJOY2016 test 74 data not downloaded");
    }
    alea_nuc_thermal_t* thermal = alea_nuc_load_thermal(xsdir, "hzrh.10t");
    ASSERT_NOT_NULL(thermal);

    double energy[] = {1.0e-11, 1.0e-5};
    double free_elastic[] = {20.0, 20.0};
    double zero[] = {0.0, 0.0};
    alea_nuc_nuclide_t hydrogen = {0};
    hydrogen.Z = 1;
    hydrogen.A = 1;
    hydrogen.particle = ALEA_NUC_PARTICLE_NEUTRON;
    hydrogen.awr = 0.999167;
    hydrogen.temperature = thermal->temperature;
    hydrogen.n_energies = 2;
    hydrogen.energy = energy;
    hydrogen.sigma_total = free_elastic;
    hydrogen.sigma_elastic = free_elastic;
    hydrogen.sigma_abs = zero;
    hydrogen.heating = zero;
    alea_nuc_mat_component_t component = {&hydrogen, 0.1};
    alea_nuc_material_t material = {&component, 1, 1, NULL};
    alea_nuc_thermal_association_t association = {0, thermal};
    alea_nuc_prepare_requirements_t requirements = {
        .required_capabilities = ALEA_NUC_CAP_RESTRICTED_NEUTRON |
                                 ALEA_NUC_CAP_THERMAL_SAB,
        .thermal_associations = &association,
        .n_thermal_associations = 1
    };
    alea_nuc_capability_report_t report;
    alea_nuc_prepared_material_t* prepared = NULL;
    ASSERT_EQ(alea_nuc_prepare_material(
                  &material, &requirements, &report, &prepared), ALEA_OK);

    alea_nuc_particle_state_t incident = {
        ALEA_NUC_PARTICLE_NEUTRON, thermal->inelastic_energy[45],
        {0.0, 0.0, 1.0}, 1.0, 0.0
    };
    alea_nuc_evaluation_t evaluation;
    ASSERT_EQ(alea_nuc_evaluate(prepared, &incident, &evaluation), ALEA_OK);
    ASSERT_NEAR(evaluation.macro_total, 6.159646592073, 1e-11);
    ASSERT_NEAR(evaluation.macro_elastic, 0.0, 1e-14);
    ASSERT_NEAR(evaluation.macro_thermal, 6.159646592073, 1e-11);

    alea_nuc_rng_t rng;
    ASSERT_EQ(alea_nuc_rng_init(&rng, 0x74a11u, 74, 1, 0,
                               ALEA_NUC_RNG_COLLISION), ALEA_OK);
    const int samples = 131072;
    int inelastic = 0;
    for (int i = 0; i < samples; i++) {
        alea_nuc_collision_result_t result;
        ASSERT_EQ(alea_nuc_collide(
                      &evaluation, alea_nuc_rng_uniform, &rng, &result),
                  ALEA_OK);
        ASSERT_TRUE(result.mt == 2 || result.mt == 4);
        if (result.mt == 4) inelastic++;
    }
    ASSERT_NEAR((double)inelastic / samples, 0.0410123, 0.004);

    alea_nuc_prepared_material_free(prepared);
    alea_nuc_thermal_free(thermal);
    alea_nuc_xsdir_free(xsdir);
}

static double skewed_outgoing_probability(int index, int count) {
    if (index == 0) return 0.1 / (count - 3);
    if (index == 1) return 0.4 / (count - 3);
    if (index == count - 2) return 0.4 / (count - 3);
    if (index == count - 1) return 0.1 / (count - 3);
    return 1.0 / (count - 3);
}

TEST(njoy25_light_water_moderator_matches_discrete_expectation) {
    alea_nuc_xsdir_t* xsdir = alea_nuc_xsdir_load(NJOY25_XSDIR);
    if (!xsdir) {
        if (getenv("ALEA_NUCDATA_REQUIRED")) ASSERT_NOT_NULL(xsdir);
        SKIP("NJOY2016 test 25 data not downloaded");
    }
    alea_nuc_thermal_t* thermal = alea_nuc_load_thermal(xsdir, "lwtr.10t");
    ASSERT_NOT_NULL(thermal);
    ASSERT_EQ(thermal->n_inelastic_energies, 117);
    ASSERT_EQ(thermal->n_inelastic_outgoing, 20);
    ASSERT_EQ(thermal->n_inelastic_cosines, 20);
    ASSERT_TRUE(thermal->inelastic_skewed);
    ASSERT_EQ(thermal->elastic_mode, ALEA_NUC_THERMAL_ELASTIC_NONE);

    const int incident_index = 48;
    const double incident_energy = thermal->inelastic_energy[incident_index];
    double expected_energy = 0.0, expected_mu = 0.0;
    double expected_upscatter = 0.0;
    for (int j = 0; j < thermal->n_inelastic_outgoing; j++) {
        double probability = skewed_outgoing_probability(
            j, thermal->n_inelastic_outgoing);
        size_t record = (size_t)incident_index *
            thermal->n_inelastic_outgoing + (size_t)j;
        double outgoing = thermal->inelastic_energy_out[record];
        expected_energy += probability * outgoing;
        if (outgoing > incident_energy) expected_upscatter += probability;
        double mu = 0.0;
        for (int k = 0; k < thermal->n_inelastic_cosines; k++)
            mu += thermal->inelastic_mu[
                record * thermal->n_inelastic_cosines + (size_t)k];
        expected_mu += probability * mu / thermal->n_inelastic_cosines;
    }

    alea_nuc_particle_state_t incident = {
        ALEA_NUC_PARTICLE_NEUTRON, incident_energy,
        {0.0, 0.0, 1.0}, 1.0, 0.0
    };
    alea_nuc_rng_t rng;
    ASSERT_EQ(alea_nuc_rng_init(&rng, 0x2510u, 25, 0, 0,
                               ALEA_NUC_RNG_COLLISION), ALEA_OK);
    const int samples = 131072;
    int upscattered = 0;
    double energy_sum = 0.0, mu_sum = 0.0;
    for (int i = 0; i < samples; i++) {
        alea_nuc_collision_result_t result;
        ASSERT_EQ(alea_nuc_sample_thermal_collision(
                      thermal, &incident, alea_nuc_rng_uniform, &rng,
                      &result), ALEA_OK);
        ASSERT_EQ(result.mt, 4);
        energy_sum += result.outgoing.energy;
        mu_sum += result.mu_lab;
        if (result.outgoing.energy > incident_energy) upscattered++;
    }
    ASSERT_NEAR(energy_sum / samples, expected_energy, 4e-10);
    ASSERT_NEAR(mu_sum / samples, expected_mu, 0.009);
    ASSERT_NEAR((double)upscattered / samples, expected_upscatter, 0.008);

    double grid[] = {1.0e-11, 1.0e-5};
    double free_elastic[] = {20.0, 20.0}, zero[] = {0.0, 0.0};
    alea_nuc_nuclide_t hydrogen = {0};
    hydrogen.Z = 1; hydrogen.A = 1;
    hydrogen.particle = ALEA_NUC_PARTICLE_NEUTRON;
    hydrogen.awr = 0.999167; hydrogen.temperature = thermal->temperature;
    hydrogen.n_energies = 2; hydrogen.energy = grid;
    hydrogen.sigma_total = free_elastic;
    hydrogen.sigma_elastic = free_elastic; hydrogen.sigma_abs = zero;
    alea_nuc_mat_component_t component = {&hydrogen, 0.1};
    alea_nuc_material_t material = {&component, 1, 1, NULL};
    alea_nuc_thermal_association_t association = {0, thermal};
    alea_nuc_prepare_requirements_t requirements = {
        .required_capabilities = ALEA_NUC_CAP_RESTRICTED_NEUTRON |
                                 ALEA_NUC_CAP_THERMAL_SAB,
        .thermal_associations = &association,
        .n_thermal_associations = 1
    };
    alea_nuc_capability_report_t report;
    alea_nuc_prepared_material_t* prepared = NULL;
    ASSERT_EQ(alea_nuc_prepare_material(
                  &material, &requirements, &report, &prepared), ALEA_OK);
    alea_nuc_evaluation_t evaluation;
    ASSERT_EQ(alea_nuc_evaluate(prepared, &incident, &evaluation), ALEA_OK);
    double expected_macro = 0.1 *
        alea_nuc_thermal_xs_inelastic(thermal, incident_energy);
    ASSERT_NEAR(evaluation.macro_total, expected_macro, 1e-12);
    ASSERT_EQ(evaluation.macro_elastic, 0.0);
    ASSERT_NEAR(evaluation.macro_thermal, expected_macro, 1e-12);

    alea_nuc_prepared_material_free(prepared);
    alea_nuc_thermal_free(thermal);
    alea_nuc_xsdir_free(xsdir);
}

TEST_MAIN()
