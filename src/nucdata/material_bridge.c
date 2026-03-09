// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/**
 * @file material_bridge.c
 * @brief Bridge between core materials and nucdata materials
 *
 * Converts a core alea_material_t (with fractions and density) into a
 * nucdata alea_nuc_material_t (with number densities and loaded ACE nuclides).
 *
 * Nuclides are loaded via the xsdir cache (alea_nuc_xsdir_get_nuclide),
 * so multiple materials sharing the same nuclide reuse the same data.
 * The nuclides are owned by the xsdir, not by the returned material.
 */

#include "alea_nucdata.h"
#include "core/alea_materials.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* Avogadro's number in barn⁻¹·cm⁻¹ units: 0.602214076 × 10²⁴ */
#define N_AVOGADRO 0.602214076

/* Neutron mass in amu */
#define NEUTRON_MASS_AMU 1.00866491595

alea_nuc_material_t* alea_nuc_material_from_core(
    const alea_material_t* mat,
    double cell_density,
    bool is_mass_density,
    alea_nuc_xsdir_t* xsdir)
{
    if (!mat || !xsdir) return NULL;

    int n = (int)mat->nuclides.count;
    if (n == 0) return NULL;

    /* Load all nuclides (cached) and collect atomic masses */
    alea_nuc_nuclide_t** nuclides = malloc((size_t)n * sizeof(*nuclides));
    double* A = malloc((size_t)n * sizeof(double));
    if (!nuclides || !A) { free(nuclides); free(A); return NULL; }

    for (int i = 0; i < n; i++) {
        const alea_nuclide_t* cn = &mat->nuclides.data[i];

        /* Form ZAID string: "92235.80c" */
        char zaid_str[32];
        if (cn->library && cn->library[0])
            snprintf(zaid_str, sizeof(zaid_str), "%d%s", cn->zaid, cn->library);
        else
            snprintf(zaid_str, sizeof(zaid_str), "%d", cn->zaid);

        nuclides[i] = alea_nuc_xsdir_get_nuclide(xsdir, zaid_str);
        if (!nuclides[i]) {
            /* No cleanup needed — cached nuclides are owned by xsdir */
            free(nuclides);
            free(A);
            return NULL;
        }

        A[i] = nuclides[i]->awr * NEUTRON_MASS_AMU;
    }

    alea_nuc_material_t* nmat = alea_nuc_material_create();
    if (!nmat) { free(nuclides); free(A); return NULL; }

    bool weight_frac = mat->is_weight_fraction;
    double rho = fabs(cell_density);

    if (weight_frac && is_mass_density) {
        /* N_i = ρ · N_A · w_i / A_i */
        for (int i = 0; i < n; i++) {
            double Ni = rho * N_AVOGADRO * mat->nuclides.data[i].fraction / A[i];
            alea_nuc_material_add(nmat, nuclides[i], Ni);
        }
    } else if (!weight_frac && is_mass_density) {
        /* Ā = Σ f_i · A_i, then N_i = (ρ · N_A / Ā) · f_i */
        double A_avg = 0.0;
        for (int i = 0; i < n; i++)
            A_avg += mat->nuclides.data[i].fraction * A[i];
        if (A_avg <= 0.0) A_avg = 1.0;
        for (int i = 0; i < n; i++) {
            double Ni = rho * N_AVOGADRO / A_avg * mat->nuclides.data[i].fraction;
            alea_nuc_material_add(nmat, nuclides[i], Ni);
        }
    } else if (!weight_frac && !is_mass_density) {
        /* N_i = N_total · f_i */
        for (int i = 0; i < n; i++) {
            double Ni = rho * mat->nuclides.data[i].fraction;
            alea_nuc_material_add(nmat, nuclides[i], Ni);
        }
    } else {
        /* Weight fractions + atom density: N_i = N_total · (w_i/A_i) / Σ(w_j/A_j) */
        double sum_wA = 0.0;
        for (int i = 0; i < n; i++)
            sum_wA += mat->nuclides.data[i].fraction / A[i];
        if (sum_wA <= 0.0) sum_wA = 1.0;
        for (int i = 0; i < n; i++) {
            double Ni = rho * (mat->nuclides.data[i].fraction / A[i]) / sum_wA;
            alea_nuc_material_add(nmat, nuclides[i], Ni);
        }
    }

    free(nuclides);
    free(A);
    return nmat;
}
