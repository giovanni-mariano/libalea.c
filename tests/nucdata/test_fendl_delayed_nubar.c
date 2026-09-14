// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/**
 * @file test_fendl_delayed_nubar.c
 * @brief Regression test for repeated energies in FENDL delayed nubar data
 */

#include "alea_nucdata.h"

#include <math.h>
#include <stdio.h>

#define FENDL_XSDIR "fendl-FENDL-3.2c-neutron-ace/neutron/ace/xsdir"

static int check_repeated_knot(const alea_nuc_xsdir_t* xsdir,
                               const char* zaid, double expected_energy) {
    alea_nuc_nuclide_t* nuc = alea_nuc_load_nuclide(xsdir, zaid);
    if (!nuc) {
        fprintf(stderr, "FAIL: could not load %s\n", zaid);
        return 0;
    }

    const alea_nuc_nu_bar_t* delayed =
        nuc->fission ? nuc->fission->delayed : NULL;
    int knot = -1;
    if (delayed && delayed->type == ALEA_NUC_NU_TABULAR) {
        for (int i = 1; i < delayed->n_energies; i++) {
            if (delayed->energy[i - 1] == expected_energy &&
                delayed->energy[i] == expected_energy) {
                knot = i;
                break;
            }
        }
    }

    int ok = knot > 0 &&
        fabs(alea_nuc_delayed_nu_bar(nuc, expected_energy) -
             delayed->nu[knot]) <= 1e-14;
    if (!ok)
        fprintf(stderr, "FAIL: %s delayed nubar knot at %.0f MeV\n",
                zaid, expected_energy);
    alea_nuc_nuclide_free(nuc);
    return ok;
}

int main(void) {
    alea_nuc_xsdir_t* xsdir = alea_nuc_xsdir_load(FENDL_XSDIR);
    if (!xsdir) {
        printf("SKIP: FENDL data not downloaded (%s)\n", FENDL_XSDIR);
        return 0;
    }

    int ok = check_repeated_knot(xsdir, "92235.32c", 20.0) &&
             check_repeated_knot(xsdir, "92238.32c", 30.0);
    alea_nuc_xsdir_free(xsdir);
    if (ok) printf("PASS: FENDL delayed nubar repeated knots\n");
    return ok ? 0 : 1;
}
