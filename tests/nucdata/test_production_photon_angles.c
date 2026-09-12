// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/** @file test_production_photon_angles.c
 *  @brief Production-library anisotropic photon angular benchmark
 */

#include "alea_nucdata.h"
#include "alea_test.h"

#include <math.h>
#include <stdlib.h>

#define CARBON_ZAID "6000.71c"
#define CARBON_PHOTON_MT 51001

static alea_nuc_xsdir_t* xsdir;
static alea_nuc_nuclide_t* carbon;
static int teardown_registered;

static void teardown(void) {
    alea_nuc_xsdir_free(xsdir);
}

static void setup(void) {
    if (xsdir) return;
    const char* path = getenv("ALEA_PRODUCTION_NEUTRON_XSDIR");
    if (!path || !path[0]) return;
    xsdir = alea_nuc_xsdir_load(path);
    if (xsdir) carbon = alea_nuc_xsdir_get_nuclide(xsdir, CARBON_ZAID);
    if (!teardown_registered) {
        atexit(teardown);
        teardown_registered = 1;
    }
}

static int production_data_requested(void) {
    const char* path = getenv("ALEA_PRODUCTION_NEUTRON_XSDIR");
    return path && path[0];
}

static const alea_nuc_photon_production_t* carbon_channel(void) {
    if (!carbon) return NULL;
    for (int i = 0; i < carbon->n_photon_productions; i++)
        if (carbon->photon_productions[i].mt == CARBON_PHOTON_MT)
            return &carbon->photon_productions[i];
    return NULL;
}

static int point_moments(const alea_nuc_angular_point_t* point,
                         double* first, double* second) {
    if (point->type == ALEA_NUC_ANG_ISOTROPIC) {
        *first = 0.0;
        *second = 1.0 / 3.0;
        return 1;
    }
    if (point->type != ALEA_NUC_ANG_EQUIPROBABLE ||
        point->n_cosines != 33 || !point->cosine)
        return 0;
    double m1 = 0.0, m2 = 0.0;
    for (int i = 0; i < 32; i++) {
        double lo = point->cosine[i];
        double hi = point->cosine[i + 1];
        m1 += 0.5 * (lo + hi) / 32.0;
        m2 += (lo * lo + lo * hi + hi * hi) / (3.0 * 32.0);
    }
    *first = m1;
    *second = m2;
    return 1;
}

static int distribution_moments(const alea_nuc_angular_dist_t* angular,
                                double energy, double* first,
                                double* second) {
    if (!angular || angular->n_energies < 2 || !angular->energy ||
        !angular->data || energy < angular->energy[0] ||
        energy > angular->energy[angular->n_energies - 1])
        return 0;
    int lo = 0;
    while (lo + 2 < angular->n_energies &&
           angular->energy[lo + 1] <= energy)
        lo++;
    int hi = lo + 1;
    double width = angular->energy[hi] - angular->energy[lo];
    double fraction = width > 0.0
        ? (energy - angular->energy[lo]) / width : 1.0;
    double lo_first, lo_second, hi_first, hi_second;
    if (!point_moments(&angular->data[lo], &lo_first, &lo_second) ||
        !point_moments(&angular->data[hi], &hi_first, &hi_second))
        return 0;
    *first = lo_first + fraction * (hi_first - lo_first);
    *second = lo_second + fraction * (hi_second - lo_second);
    return 1;
}

TEST(production_carbon_photon_angular_metadata) {
    setup();
    if (!production_data_requested())
        SKIP("ALEA_PRODUCTION_NEUTRON_XSDIR is not set");
    ASSERT_NOT_NULL(xsdir);
    ASSERT_NOT_NULL(carbon);
    const alea_nuc_photon_production_t* channel = carbon_channel();
    ASSERT_NOT_NULL(channel);
    ASSERT_EQ(channel->parent_mt, 51);
    ASSERT_EQ(channel->mf, 12);
    ASSERT_NOT_NULL(channel->angular);
    ASSERT_EQ(channel->angular->n_energies, 72);
    int isotropic = 0, equiprobable = 0;
    for (int i = 0; i < channel->angular->n_energies; i++) {
        isotropic += channel->angular->data[i].type ==
                     ALEA_NUC_ANG_ISOTROPIC;
        equiprobable += channel->angular->data[i].type ==
                        ALEA_NUC_ANG_EQUIPROBABLE;
    }
    ASSERT_EQ(isotropic, 2);
    ASSERT_EQ(equiprobable, 70);
}

TEST(production_carbon_photon_angular_moments) {
    setup();
    if (!production_data_requested())
        SKIP("ALEA_PRODUCTION_NEUTRON_XSDIR is not set");
    ASSERT_NOT_NULL(xsdir);
    const alea_nuc_photon_production_t* channel = carbon_channel();
    ASSERT_NOT_NULL(channel);
    const double incident_energy = 14.0;
    double expected_first, expected_second;
    ASSERT_TRUE(distribution_moments(
        channel->angular, incident_energy, &expected_first, &expected_second));
    ASSERT_NEAR(expected_first, 3.9700889284399743e-8, 1e-14);
    ASSERT_NEAR(expected_second, 0.42500765627646159, 1e-14);

    alea_nuc_particle_state_t incident = {
        ALEA_NUC_PARTICLE_NEUTRON, incident_energy,
        {0.0, 0.0, 1.0}, 1.0, 0.0
    };
    alea_nuc_rng_t rng;
    ASSERT_EQ(alea_nuc_rng_init(&rng, 0xc12u, 51, 0, 23,
                               ALEA_NUC_RNG_COLLISION), ALEA_OK);
    const int samples = 131072;
    double first = 0.0, second = 0.0;
    for (int i = 0; i < samples; i++) {
        alea_nuc_particle_state_t photon;
        ASSERT_EQ(alea_nuc_sample_photon_production(
                      carbon, channel, &incident, alea_nuc_rng_uniform, &rng,
                      &photon), ALEA_OK);
        double mu = photon.direction[2];
        ASSERT_TRUE(mu >= -1.0 && mu <= 1.0);
        first += mu;
        second += mu * mu;
    }
    ASSERT_NEAR(first / samples, expected_first, 0.006);
    ASSERT_NEAR(second / samples, expected_second, 0.006);
}

TEST_MAIN()
