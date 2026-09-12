// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/** @file test_photoatomic_relaxation.c
 *  @brief Pinned validation of decoded photoatomic fluorescence metadata
 */

#include "alea_nucdata.h"
#include "alea_test.h"

#define XSDIR "njoy-test59/xsdir"

TEST(njoy59_decodes_averaged_fluorescence_block) {
    alea_nuc_xsdir_t* xsdir = alea_nuc_xsdir_load(XSDIR);
    if (!xsdir) SKIP("NJOY2016 test 59 data not downloaded");
    alea_nuc_nuclide_t* uranium =
        alea_nuc_xsdir_get_nuclide(xsdir, "92000.31p");
    ASSERT_NOT_NULL(uranium);
    ASSERT_NOT_NULL(uranium->photon);
    const alea_nuc_photon_data_t* photon = uranium->photon;
    ASSERT_EQ(photon->n_fluorescence, 6);

    const double edge[] = {
        0.019962, 0.019962, 0.11561, 0.11561, 0.11561, 0.11561
    };
    const double phi[] = {
        0.25821941822, 1.0, 2.57079117311,
        3.55254965665, 4.10152778392, 4.2424694148
    };
    const double yield[] = {
        0.0, 0.329284393225, 1.82993490556,
        2.76785485095, 3.29231941509, 3.42696756834
    };
    const double energy[] = {
        0.0, 0.014836536586, 0.098928,
        0.095066, 0.111541736024, 0.115012586119
    };
    for (int i = 0; i < 6; i++) {
        ASSERT_NEAR(photon->fluorescence_edge[i], edge[i], 1e-14);
        ASSERT_NEAR(photon->fluorescence_phi[i], phi[i], 1e-13);
        ASSERT_NEAR(photon->fluorescence_yield[i], yield[i], 1e-13);
        ASSERT_NEAR(photon->fluorescence_energy[i], energy[i], 1e-13);
    }
    alea_nuc_xsdir_free(xsdir);
}

TEST_MAIN()
