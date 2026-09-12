// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/** @file test_photon_production.c
 *  @brief Independent photon-production checks with pinned NJOY2016 fixtures
 */

#include "alea_nucdata.h"
#include "alea_test.h"

#define NJOY07_XSDIR "njoy-test07/xsdir"
#define NJOY08_XSDIR "njoy-test08/xsdir"

static alea_nuc_xsdir_t* xsdir07;
static alea_nuc_xsdir_t* xsdir08;
static alea_nuc_nuclide_t* u235;
static alea_nuc_nuclide_t* ni61;
static int teardown_registered;
static void teardown(void);

static void register_teardown(void) {
    if (!teardown_registered) {
        atexit(teardown);
        teardown_registered = 1;
    }
}

static void setup07(void) {
    register_teardown();
    if (xsdir07) return;
    xsdir07 = alea_nuc_xsdir_load(NJOY07_XSDIR);
    if (xsdir07)
        u235 = alea_nuc_xsdir_get_nuclide(xsdir07, "92235.00c");
}

static void setup08(void) {
    register_teardown();
    if (xsdir08) return;
    xsdir08 = alea_nuc_xsdir_load(NJOY08_XSDIR);
    if (xsdir08)
        ni61 = alea_nuc_xsdir_get_nuclide(xsdir08, "28061.00c");
}

static const alea_nuc_photon_production_t* production(
    const alea_nuc_nuclide_t* nuc, int mt) {
    if (!nuc) return NULL;
    for (int i = 0; i < nuc->n_photon_productions; i++)
        if (nuc->photon_productions[i].mt == mt)
            return &nuc->photon_productions[i];
    return NULL;
}

static double fixed_random(void* context) {
    (void)context;
    return 0.5;
}

static double tabular_mean(const alea_nuc_law67_energy_t* spectrum) {
    double integral = 0.0;
    for (int i = 0; i + 1 < spectrum->n_points; i++) {
        double x0 = spectrum->energy[i];
        double dx = spectrum->energy[i + 1] - x0;
        double p0 = spectrum->pdf[i];
        if (spectrum->interpolation == 1) {
            integral += p0 * (x0 * dx + 0.5 * dx * dx);
        } else {
            double slope = (spectrum->pdf[i + 1] - p0) / dx;
            integral += x0 * p0 * dx +
                0.5 * (x0 * slope + p0) * dx * dx +
                slope * dx * dx * dx / 3.0;
        }
    }
    return integral;
}

TEST(njoy07_photon_blocks_decode) {
    setup07();
    if (!xsdir07) SKIP("NJOY2016 test 07 data not downloaded");
    ASSERT_NOT_NULL(u235);
    ASSERT_EQ(u235->n_photon_productions, 33);
    ASSERT_EQ(u235->n_photon_yield_multipliers, 3);
    ASSERT_EQ(u235->photon_yield_multipliers[0], 4);
    ASSERT_EQ(u235->photon_yield_multipliers[1], 18);
    ASSERT_EQ(u235->photon_yield_multipliers[2], 102);
    ASSERT_NOT_NULL(u235->total_photon_production_xs);
}

TEST(njoy07_gpd_matches_independent_reference) {
    setup07();
    if (!u235) SKIP("no data");
    ASSERT_NEAR(alea_nuc_xs_photon_production_total(u235, 1.0),
                12.11682, 1e-10);
    ASSERT_NEAR(alea_nuc_xs_photon_production_total(u235, 2.0),
                14.83, 1e-10);
    ASSERT_NEAR(alea_nuc_xs_photon_production_total(u235, 14.0),
                26.59129, 1e-10);
}

TEST(njoy07_repeated_yield_knots_match_independent_reference) {
    setup07();
    if (!u235) SKIP("no data");
    const alea_nuc_photon_production_t* inelastic_line =
        production(u235, 4001);
    const alea_nuc_photon_production_t* fission = production(u235, 18001);
    ASSERT_NOT_NULL(inelastic_line);
    ASSERT_NOT_NULL(fission);
    ASSERT_NEAR(alea_nuc_photon_production_yield(
                    u235, inelastic_line, 1.0), 0.119, 1e-14);
    ASSERT_NEAR(alea_nuc_photon_production_yield(
                    u235, inelastic_line, 1.09), 0.0, 1e-14);
    ASSERT_NEAR(alea_nuc_photon_production_yield(
                    u235, fission, 1.0), 7.17, 1e-14);
    ASSERT_NEAR(alea_nuc_photon_production_yield(
                    u235, fission, 1.09), 0.0, 1e-14);
}

TEST(njoy07_aggregate_parent_cross_sections_match_reference) {
    setup07();
    if (!u235) SKIP("no data");
    ASSERT_NEAR(alea_nuc_xs_reaction(u235, 3, 2.0), 3.1021611, 1e-10);
    ASSERT_NEAR(alea_nuc_xs_reaction(u235, 4, 2.0), 1.75, 1e-12);
    ASSERT_NEAR(alea_nuc_xs_reaction(u235, 18, 2.0), 1.298, 1e-12);
    const alea_nuc_photon_production_t* continuum = production(u235, 3001);
    ASSERT_NOT_NULL(continuum);
    ASSERT_TRUE(continuum->production_xs);
    ASSERT_NEAR(alea_nuc_photon_production_yield(
                    u235, continuum, 2.0), 4.780538315692245, 1e-12);
}

TEST(njoy07_discrete_photon_energy_matches_independent_reference) {
    setup07();
    if (!u235) SKIP("no data");
    const alea_nuc_photon_production_t* line = production(u235, 4001);
    ASSERT_NOT_NULL(line);
    ASSERT_NOT_NULL(line->spectrum);
    ASSERT_EQ(line->spectrum->law, ALEA_NUC_ELAW_DISCRETE_PHOTON);
    double energy = 0.0;
    ASSERT_EQ(alea_nuc_sample_energy_distribution(
                  line->spectrum, 1.0, fixed_random, NULL, &energy), ALEA_OK);
    ASSERT_NEAR(energy, 0.771, 1e-14);
}

TEST(njoy07_negative_watt_restriction_is_collision_ready) {
    setup07();
    if (!u235) SKIP("NJOY2016 test 07 data not downloaded");
    const alea_nuc_reaction_t* chance_fission = NULL;
    for (int i = 0; i < u235->n_reactions; i++)
        if (u235->reactions[i].mt == 19)
            chance_fission = &u235->reactions[i];
    ASSERT_NOT_NULL(chance_fission);
    ASSERT_NOT_NULL(chance_fission->energy);
    ASSERT_EQ(chance_fission->energy->law, ALEA_NUC_ELAW_WATT);
    ASSERT_NEAR(chance_fission->energy->restriction_energy, -30.0, 1e-14);

    alea_nuc_capability_report_t report;
    ASSERT_EQ(alea_nuc_capabilities(u235, &report), ALEA_OK);
    ASSERT_TRUE(report.available_capabilities & ALEA_NUC_CAP_FISSION);
    ASSERT_TRUE(report.available_capabilities & ALEA_NUC_CAP_DELAYED_NEUTRON);

    alea_nuc_rng_t rng;
    ASSERT_EQ(alea_nuc_rng_init(&rng, 0x235u, 19, 0, 11,
                               ALEA_NUC_RNG_COLLISION), ALEA_OK);
    const int samples = 65536;
    double sum = 0.0;
    for (int i = 0; i < samples; i++) {
        double energy = -1.0;
        ASSERT_EQ(alea_nuc_sample_energy_distribution(
                      chance_fission->energy, 1.0, alea_nuc_rng_uniform,
                      &rng, &energy), ALEA_OK);
        ASSERT_TRUE(energy >= 0.0 && energy <= 31.0);
        sum += energy;
    }
    /* For an effectively untruncated Watt law, E[E'] = 3a/2 + a^2 b/4. */
    ASSERT_NEAR(sum / samples, 2.03084, 0.02);
}

TEST(njoy08_mf16_channels_and_yields_match_reference) {
    setup08();
    if (!xsdir08) SKIP("NJOY2016 test 08 data not downloaded");
    ASSERT_NOT_NULL(ni61);
    ASSERT_EQ(ni61->n_photon_productions, 36);

    const int mts[] = {16001, 28001, 91001, 103001, 107001};
    const double energies[] = {10.0, 14.0, 5.0, 10.0, 20.0};
    const double expected[] = {
        0.23651, 0.7764133333333333, 2.1605, 2.0155, 4.0495
    };
    int mf16_count = 0;
    for (int i = 0; i < ni61->n_photon_productions; i++)
        if (ni61->photon_productions[i].mf == 16) mf16_count++;
    ASSERT_EQ(mf16_count, 5);
    for (int i = 0; i < 5; i++) {
        const alea_nuc_photon_production_t* channel =
            production(ni61, mts[i]);
        ASSERT_NOT_NULL(channel);
        ASSERT_EQ(channel->mf, 16);
        ASSERT_NEAR(alea_nuc_photon_production_yield(
                        ni61, channel, energies[i]), expected[i], 1e-12);
    }
}

TEST(njoy08_law67_joint_moments_match_direct_integration) {
    setup08();
    if (!ni61) SKIP("NJOY2016 test 08 data not downloaded");
    const alea_nuc_reaction_t* reaction = NULL;
    for (int i = 0; i < ni61->n_reactions; i++)
        if (ni61->reactions[i].mt == 16) reaction = &ni61->reactions[i];
    ASSERT_NOT_NULL(reaction);
    ASSERT_NOT_NULL(reaction->energy);
    const alea_nuc_energy_dist_t* law = reaction->energy;
    ASSERT_EQ(law->law, ALEA_NUC_ELAW_LAB_ANGLE_ENERGY);
    ASSERT_EQ(law->law67.n_ein, 13);

    const int incident_index = 8;
    const alea_nuc_law67_incident_t* distribution =
        &law->law67.incident[incident_index];
    ASSERT_EQ(distribution->n_cosines, 2);
    double expected_mu = 0.0, expected_energy = 0.0;
    for (int i = 0; i + 1 < distribution->n_cosines; i++) {
        expected_mu += 0.5 *
            (distribution->cosine[i] + distribution->cosine[i + 1]);
        expected_energy += 0.5 *
            (tabular_mean(&distribution->spectrum[i]) +
             tabular_mean(&distribution->spectrum[i + 1]));
    }
    expected_mu /= distribution->n_cosines - 1;
    expected_energy /= distribution->n_cosines - 1;

    alea_nuc_capability_report_t report;
    ASSERT_EQ(alea_nuc_capabilities(ni61, &report), ALEA_OK);
    alea_nuc_rng_t rng;
    ASSERT_EQ(alea_nuc_rng_init(&rng, 0x6716u, 61, 16, 67,
                               ALEA_NUC_RNG_COLLISION), ALEA_OK);
    const int samples = 131072;
    double energy_sum = 0.0, mu_sum = 0.0;
    for (int i = 0; i < samples; i++) {
        double energy = -1.0, mu = -2.0;
        bool correlated = false;
        ASSERT_EQ(alea_nuc_sample_energy_angle_distribution(
                      law, law->law67.ein[incident_index],
                      alea_nuc_rng_uniform, &rng, &energy, &mu,
                      &correlated), ALEA_OK);
        ASSERT_TRUE(correlated);
        ASSERT_TRUE(energy >= 0.0);
        ASSERT_TRUE(mu >= -1.0 && mu <= 1.0);
        energy_sum += energy;
        mu_sum += mu;
    }
    ASSERT_NEAR(energy_sum / samples, expected_energy, 0.015);
    ASSERT_NEAR(mu_sum / samples, expected_mu, 0.008);
}

TEST(njoy08_continuous_spectrum_and_isotropic_angle_statistics) {
    setup08();
    if (!ni61) SKIP("NJOY2016 test 08 data not downloaded");
    const alea_nuc_photon_production_t* channel = production(ni61, 16001);
    ASSERT_NOT_NULL(channel);
    ASSERT_NOT_NULL(channel->spectrum);
    ASSERT_EQ(channel->spectrum->law, ALEA_NUC_ELAW_CONT_TABULAR);

    alea_nuc_particle_state_t incident = {
        ALEA_NUC_PARTICLE_NEUTRON, 13.0, {0.0, 0.0, 1.0}, 1.0, 0.0
    };
    alea_nuc_rng_t rng;
    ASSERT_EQ(alea_nuc_rng_init(&rng, 0x8e71u, 61, 0, 16,
                               ALEA_NUC_RNG_COLLISION), ALEA_OK);
    const int samples = 65536;
    double energy_sum = 0.0;
    double mu_sum = 0.0;
    for (int i = 0; i < samples; i++) {
        alea_nuc_particle_state_t photon;
        ASSERT_EQ(alea_nuc_sample_photon_production(
                      ni61, channel, &incident, alea_nuc_rng_uniform, &rng,
                      &photon), ALEA_OK);
        energy_sum += photon.energy;
        mu_sum += photon.direction[2];
    }
    ASSERT_NEAR(energy_sum / samples, 1.2823561086446267, 0.015);
    ASSERT_NEAR(mu_sum / samples, 0.0, 0.008);
}

TEST(njoy08_core_rng_replays_photon_samples) {
    setup08();
    if (!ni61) SKIP("NJOY2016 test 08 data not downloaded");
    const alea_nuc_photon_production_t* channel = production(ni61, 16001);
    ASSERT_NOT_NULL(channel);
    alea_nuc_particle_state_t incident = {
        ALEA_NUC_PARTICLE_NEUTRON, 13.0, {0.0, 0.0, 1.0}, 1.0, 0.0
    };
    alea_nuc_rng_t first, replay;
    ASSERT_EQ(alea_nuc_rng_init(&first, 1234, 9, 2, 7,
                               ALEA_NUC_RNG_COLLISION), ALEA_OK);
    replay = first;
    for (int i = 0; i < 256; i++) {
        alea_nuc_particle_state_t a, b;
        ASSERT_EQ(alea_nuc_sample_photon_production(
                      ni61, channel, &incident, alea_nuc_rng_uniform, &first,
                      &a), ALEA_OK);
        ASSERT_EQ(alea_nuc_sample_photon_production(
                      ni61, channel, &incident, alea_nuc_rng_uniform, &replay,
                      &b), ALEA_OK);
        ASSERT_TRUE(a.energy == b.energy);
        ASSERT_TRUE(a.direction[0] == b.direction[0]);
        ASSERT_TRUE(a.direction[1] == b.direction[1]);
        ASSERT_TRUE(a.direction[2] == b.direction[2]);
    }
}

static void teardown(void) {
    alea_nuc_xsdir_free(xsdir07);
    alea_nuc_xsdir_free(xsdir08);
}

TEST_MAIN()
