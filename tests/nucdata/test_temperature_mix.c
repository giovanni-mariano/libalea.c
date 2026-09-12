// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/** @file test_temperature_mix.c
 *  @brief Temperature-table mixing over resonance-bearing evaluated data
 */

#include "alea_nucdata.h"
#include "alea_test.h"

#include <math.h>

#define XSDIR "njoy-test07/xsdir"
#define ZAID "92235.00c"

TEST(real_resonance_table_mix_preserves_expectation_and_event_weights) {
    alea_nuc_xsdir_t* xsdir = alea_nuc_xsdir_load(XSDIR);
    if (!xsdir) SKIP("NJOY2016 test 07 data not downloaded");
    alea_nuc_nuclide_t* lower = alea_nuc_load_nuclide(xsdir, ZAID);
    alea_nuc_nuclide_t* upper = alea_nuc_load_nuclide(xsdir, ZAID);
    ASSERT_NOT_NULL(lower);
    ASSERT_NOT_NULL(upper);

    const double upper_temperature = 3.0 * lower->temperature;
    ASSERT_EQ(alea_nuc_doppler_broaden(upper, upper_temperature), ALEA_OK);
    ASSERT_EQ(upper->temperature, upper_temperature);
    const double fraction = 0.5;
    const double density = 0.02;
    alea_nuc_material_t* material = alea_nuc_material_create();
    ASSERT_NOT_NULL(material);
    ASSERT_EQ(alea_nuc_material_add_temperature_mix(
                  material, lower, upper, fraction, density), ALEA_OK);

    const double energies[] = {1.0e-8, 1.0e-6, 1.0e-5, 1.0e-3, 1.0};
    for (size_t i = 0; i < sizeof(energies) / sizeof(energies[0]); i++) {
        double lower_xs = alea_nuc_xs_total(lower, energies[i]);
        double upper_xs = alea_nuc_xs_total(upper, energies[i]);
        double expected = density *
            ((1.0 - fraction) * lower_xs + fraction * upper_xs);
        ASSERT_NEAR(alea_nuc_mat_xs_total(material, energies[i]), expected,
                    fmax(1.0e-13, fabs(expected) * 2.0e-14));
    }

    const double sample_energy = 1.0e-5;
    double lower_xs = alea_nuc_xs_total(lower, sample_energy);
    double upper_xs = alea_nuc_xs_total(upper, sample_energy);
    double expected_upper = upper_xs / (lower_xs + upper_xs);
    alea_nuc_rng_t rng;
    ASSERT_EQ(alea_nuc_rng_init(&rng, 0x71e5u, 235, 0, 4,
                                ALEA_NUC_RNG_COLLISION), ALEA_OK);
    const int samples = 65536;
    int selected_upper = 0;
    for (int i = 0; i < samples; i++) {
        alea_nuc_nuclide_t* selected = NULL;
        int component = alea_nuc_sample_nuclide(
            material, sample_energy, alea_nuc_rng_uniform(&rng), &selected);
        ASSERT_TRUE(component == 0 || component == 1);
        ASSERT_TRUE(selected == (component == 0 ? lower : upper));
        if (component == 1) selected_upper++;
    }
    ASSERT_NEAR((double)selected_upper / samples, expected_upper, 0.006);

    alea_nuc_material_destroy(material);
    alea_nuc_nuclide_free(lower);
    alea_nuc_nuclide_free(upper);
    alea_nuc_xsdir_free(xsdir);
}

TEST_MAIN()
