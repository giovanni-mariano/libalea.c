// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

#include "alea_multigroup.h"
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

alea_error_t alea_mg_neutron_material_build(
    const alea_nuc_multigroup_t* const* nuclides,
    const double* number_densities, size_t nuclide_count,
    size_t groups, double* total, double* transfer) {
    if (!nuclides || !number_densities || !total || !transfer)
        return ALEA_ERR_NULL_ARG;
    if (!nuclide_count || !groups || groups > INT32_MAX ||
        groups > SIZE_MAX/groups)
        return ALEA_ERR_INVALID_ARG;
    const alea_nuc_multigroup_t* reference = nuclides[0];
    if (!reference || reference->n_groups != (int)groups ||
        !reference->bounds) return ALEA_ERR_INVALID_ARG;
    for (size_t g = 0; g <= groups; ++g) {
        if (!isfinite(reference->bounds[g]) || reference->bounds[g] < 0.0 ||
            (g && reference->bounds[g - 1] <= reference->bounds[g]))
            return ALEA_ERR_INVALID_ARG;
    }
    for (size_t i = 0; i < nuclide_count; ++i) {
        const alea_nuc_multigroup_t* mg = nuclides[i];
        if (!mg || mg->n_groups != (int)groups || !mg->bounds ||
            !mg->sigma_t || !mg->sigma_f || !mg->nu_sigma_f || !mg->scatter ||
            !isfinite(number_densities[i]) || number_densities[i] < 0.0)
            return ALEA_ERR_INVALID_ARG;
        for (size_t g = 0; g <= groups; ++g)
            if (mg->bounds[g] != reference->bounds[g])
                return ALEA_ERR_INVALID_ARG;
        for (size_t g = 0; g < groups; ++g) {
            if (!isfinite(mg->sigma_t[g]) || mg->sigma_t[g] < 0.0 ||
                !isfinite(mg->sigma_f[g]) || mg->sigma_f[g] < 0.0 ||
                !isfinite(mg->nu_sigma_f[g]) || mg->nu_sigma_f[g] < 0.0)
                return ALEA_ERR_INVALID_ARG;
            if (mg->sigma_f[g] > 0.0 || mg->nu_sigma_f[g] > 0.0)
                return ALEA_ERR_UNSUPPORTED;
            for (size_t h = 0; h < groups; ++h) {
                double xs = mg->scatter[g*groups + h];
                if (!isfinite(xs) || xs < 0.0 ||
                    (mg->sigma_t[g] == 0.0 && xs > 0.0))
                    return ALEA_ERR_INVALID_ARG;
            }
        }
    }
    size_t matrix_size = groups*groups;
    if (groups > SIZE_MAX - matrix_size ||
        groups + matrix_size > SIZE_MAX/sizeof(double))
        return ALEA_ERR_OVERFLOW;
    double* scratch = calloc(groups + matrix_size, sizeof(double));
    if (!scratch) return ALEA_ERR_OUT_OF_MEMORY;
    double* sum_total = scratch;
    double* sum_transfer = scratch + groups;
    alea_error_t err = ALEA_OK;
    for (size_t i = 0; i < nuclide_count; ++i) {
        const alea_nuc_multigroup_t* mg = nuclides[i];
        double density = number_densities[i];
        for (size_t g = 0; g < groups; ++g) {
            sum_total[g] += density * mg->sigma_t[g];
            if (!isfinite(sum_total[g])) { err = ALEA_ERR_OVERFLOW; break; }
        }
        if (err != ALEA_OK) break;
        for (size_t j = 0; j < matrix_size; ++j) {
            sum_transfer[j] += density * mg->scatter[j];
            if (!isfinite(sum_transfer[j])) { err = ALEA_ERR_OVERFLOW; break; }
        }
        if (err != ALEA_OK) break;
    }
    if (err == ALEA_OK) {
        memcpy(total, sum_total, groups*sizeof(double));
        memcpy(transfer, sum_transfer, matrix_size*sizeof(double));
    }
    free(scratch);
    return err;
}
