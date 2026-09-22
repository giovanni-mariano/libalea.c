// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/** @file alea_multigroup.h
 * Material data shared by forward and adjoint multigroup transport.
 */
#ifndef ALEA_MULTIGROUP_H
#define ALEA_MULTIGROUP_H

#include "alea_nucdata.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    const double* total;    /* G macroscopic totals, 1/cm */
    const double* transfer; /* G*G forward expected production, 1/cm;
                             * transfer[incoming*G + outgoing] */
} alea_mg_material_t;

/** Combine already collapsed neutron nuclides into macroscopic cell data.
 * All inputs must use identical descending group bounds. Number densities
 * are atoms/(barn cm); nuclide group cross sections are barns. Fission data
 * are rejected because the current collapse does not include fission in its
 * production matrix. Thermal and angular approximations of the supplied
 * collapsed data remain in effect.
 *
 * The caller owns total[G] and transfer[G*G], which are overwritten only on
 * success. Set an alea_mg_material_t to these arrays for transport. */
alea_error_t alea_mg_neutron_material_build(
    const alea_nuc_multigroup_t* const* nuclides,
    const double* number_densities, size_t nuclide_count,
    size_t groups, double* total, double* transfer);

#ifdef __cplusplus
}
#endif

#endif
