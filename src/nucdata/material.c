// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/**
 * @file material.c
 * @brief Material composition and macroscopic cross-section lookup
 *
 * A material is a set of nuclides with number densities.
 * Macroscopic XS: Σ = Σᵢ Nᵢ · σᵢ(E)
 */

#include "alea_nucdata.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

alea_nuc_material_t* alea_nuc_material_create(void) {
    alea_nuc_material_t* mat = calloc(1, sizeof(*mat));
    return mat;
}

void alea_nuc_material_destroy(alea_nuc_material_t* mat) {
    if (!mat) return;
    free(mat->components);
    free(mat);
}

alea_error_t alea_nuc_material_add(alea_nuc_material_t* mat, alea_nuc_nuclide_t* nuclide,
                              double number_density) {
    if (!mat || !nuclide) return ALEA_ERR_NULL_ARG;

    /* Grow array (doubling, matches alea_vec pattern) */
    if (mat->n_components >= mat->capacity) {
        int new_cap = mat->capacity ? mat->capacity * 2 : 4;
        alea_nuc_mat_component_t* p = realloc(mat->components,
                                          (size_t)new_cap * sizeof(*p));
        if (!p) return ALEA_ERR_OUT_OF_MEMORY;
        mat->components = p;
        mat->capacity = new_cap;
    }

    mat->components[mat->n_components].nuclide = nuclide;
    mat->components[mat->n_components].number_density = number_density;
    mat->n_components++;

    return ALEA_OK;
}

double alea_nuc_mat_xs_total(const alea_nuc_material_t* mat, double energy) {
    if (!mat) return 0.0;
    double sigma = 0.0;
    for (int i = 0; i < mat->n_components; i++) {
        sigma += mat->components[i].number_density *
                 alea_nuc_xs_total(mat->components[i].nuclide, energy);
    }
    return sigma;
}

double alea_nuc_mat_xs_absorption(const alea_nuc_material_t* mat, double energy) {
    if (!mat) return 0.0;
    double sigma = 0.0;
    for (int i = 0; i < mat->n_components; i++) {
        sigma += mat->components[i].number_density *
                 alea_nuc_xs_absorption(mat->components[i].nuclide, energy);
    }
    return sigma;
}

double alea_nuc_mat_xs_elastic(const alea_nuc_material_t* mat, double energy) {
    if (!mat) return 0.0;
    double sigma = 0.0;
    for (int i = 0; i < mat->n_components; i++) {
        sigma += mat->components[i].number_density *
                 alea_nuc_xs_elastic(mat->components[i].nuclide, energy);
    }
    return sigma;
}

int alea_nuc_sample_nuclide(const alea_nuc_material_t* mat, double energy, double xi,
                        alea_nuc_nuclide_t** out_nuclide) {
    if (!mat || mat->n_components <= 0 || !out_nuclide) return -1;

    /* Single pass: accumulate cumulative macroscopic XS */
    int nc = mat->n_components;
    double cumul = 0.0;
    double partials[64];
    double* p = (nc <= 64) ? partials : malloc((size_t)nc * sizeof(double));
    if (!p) return -1;

    for (int i = 0; i < nc; i++) {
        cumul += mat->components[i].number_density *
                 alea_nuc_xs_total(mat->components[i].nuclide, energy);
        p[i] = cumul;
    }
    if (cumul <= 0.0) { if (p != partials) free(p); return -1; }

    double threshold = xi * cumul;
    for (int i = 0; i < nc; i++) {
        if (p[i] >= threshold) {
            if (p != partials) free(p);
            *out_nuclide = mat->components[i].nuclide;
            return i;
        }
    }

    /* Rounding: return last */
    if (p != partials) free(p);
    *out_nuclide = mat->components[nc - 1].nuclide;
    return nc - 1;
}

int alea_nuc_sample_reaction(const alea_nuc_nuclide_t* nuc, double energy, double xi,
                         int* out_mt) {
    if (!nuc || !out_mt) return -1;

    double sigma_t = alea_nuc_xs_total(nuc, energy);
    if (sigma_t <= 0.0) return -1;

    /* Single binary search on main energy grid — reuse for all reactions */
    double f;
    int ie = alea_nuc_energy_lookup(nuc->energy, nuc->n_energies, energy, &f);
    if (ie < 0) { *out_mt = 2; return 0; }

    double threshold = xi * sigma_t;
    double cumul = 0.0;

    /* Elastic (direct interpolation, no binary search) */
    if (nuc->sigma_elastic) {
        cumul += nuc->sigma_elastic[ie] +
                 f * (nuc->sigma_elastic[ie + 1] - nuc->sigma_elastic[ie]);
        if (cumul >= threshold) { *out_mt = 2; return 0; }
    }

    /* Non-elastic reactions (direct sub-grid access, no per-reaction search) */
    for (int i = 0; i < nuc->n_reactions; i++) {
        const alea_nuc_reaction_t* r = &nuc->reactions[i];
        if (r->n_energies <= 0 || !r->xs) continue;

        int ie_start = r->threshold_index - 1;
        int ri = ie - ie_start;
        double sig_r = 0.0;
        if (ri >= 0 && ri < r->n_energies - 1)
            sig_r = r->xs[ri] + f * (r->xs[ri + 1] - r->xs[ri]);
        else if (ri == r->n_energies - 1)
            sig_r = r->xs[ri];

        cumul += sig_r;
        if (cumul >= threshold) {
            *out_mt = r->mt;
            return i + 1;
        }
    }

    /* Rounding fallback */
    if (nuc->n_reactions > 0) {
        *out_mt = nuc->reactions[nuc->n_reactions - 1].mt;
        return nuc->n_reactions;
    }
    *out_mt = 2;
    return 0;
}

double alea_nuc_mean_free_path(const alea_nuc_material_t* mat, double energy) {
    double sigma_t = alea_nuc_mat_xs_total(mat, energy);
    if (sigma_t <= 0.0) return INFINITY;
    return 1.0 / sigma_t;
}

double alea_nuc_sample_distance(const alea_nuc_material_t* mat, double energy, double xi) {
    double sigma_t = alea_nuc_mat_xs_total(mat, energy);
    if (sigma_t <= 0.0) return INFINITY;
    return -log(1.0 - xi) / sigma_t;
}
