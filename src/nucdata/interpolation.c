// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/**
 * @file interpolation.c
 * @brief Shared validated ENDF interpolation and ACE tabular-PDF primitives
 */

#include "nuclear_internal.h"

#include <math.h>

bool alea_nuc_interp_regions_valid(const int* nbt, const int* interp,
                                   int n_regions, int n_points) {
    if (n_points <= 0 || n_regions < 0) return false;
    if (n_regions == 0) return true;
    if (!nbt || !interp || n_points < 2) return false;
    for (int i = 0; i < n_regions; i++) {
        if (nbt[i] < 2 || nbt[i] > n_points ||
            (i > 0 && nbt[i] <= nbt[i - 1]) ||
            interp[i] < 1 || interp[i] > 5) return false;
    }
    return nbt[n_regions - 1] == n_points;
}

int alea_nuc_interp_code_for_interval(const int* nbt, const int* interp,
                                      int n_regions, int interval) {
    if (!nbt || !interp || n_regions <= 0 || interval < 0) return 2;
    int upper_point = interval + 2;
    for (int i = 0; i < n_regions; i++)
        if (upper_point <= nbt[i]) return interp[i];
    return interp[n_regions - 1];
}

alea_error_t alea_nuc_interp_pair(double x0, double x1,
                                  double y0, double y1, double query,
                                  int interpolation, double* value) {
    if (!isfinite(x0) || !isfinite(x1) || !isfinite(y0) || !isfinite(y1) ||
        !isfinite(query) || x1 <= x0 || !value)
        return ALEA_ERR_INVALID_ARG;
    if (query <= x0) { *value = y0; return ALEA_OK; }
    if (query >= x1) { *value = y1; return ALEA_OK; }

    double f;
    switch (interpolation) {
    case 1:
        *value = y0;
        return ALEA_OK;
    case 2:
        f = (query - x0) / (x1 - x0);
        *value = y0 + f * (y1 - y0);
        return ALEA_OK;
    case 3:
        if (x0 <= 0.0 || query <= 0.0) return ALEA_ERR_INVALID_ARG;
        f = log(query / x0) / log(x1 / x0);
        *value = y0 + f * (y1 - y0);
        return ALEA_OK;
    case 4:
        if (y0 <= 0.0 || y1 <= 0.0) return ALEA_ERR_INVALID_ARG;
        f = (query - x0) / (x1 - x0);
        *value = y0 * exp(f * log(y1 / y0));
        return ALEA_OK;
    case 5:
        if (x0 <= 0.0 || query <= 0.0 || y0 <= 0.0 || y1 <= 0.0)
            return ALEA_ERR_INVALID_ARG;
        f = log(query / x0) / log(x1 / x0);
        *value = y0 * exp(f * log(y1 / y0));
        return ALEA_OK;
    default:
        return ALEA_ERR_UNSUPPORTED;
    }
}

alea_error_t alea_nuc_interp_eval(const double* x, const double* y, int n,
                                  const int* nbt, const int* interp,
                                  int n_regions, double query, double* value) {
    if (!x || !y || !value || n <= 0 || !isfinite(query) ||
        !alea_nuc_interp_regions_valid(nbt, interp, n_regions, n))
        return ALEA_ERR_INVALID_ARG;
    for (int i = 0; i < n; i++) {
        if (!isfinite(x[i]) || !isfinite(y[i]) ||
            (i > 0 && x[i] < x[i - 1])) return ALEA_ERR_INVALID_ARG;
    }
    if (n == 1 || query <= x[0]) { *value = y[0]; return ALEA_OK; }
    if (query >= x[n - 1]) { *value = y[n - 1]; return ALEA_OK; }
    double ignored;
    int interval = alea_nuc_energy_lookup_trusted(x, n, query, &ignored);
    int code = alea_nuc_interp_code_for_interval(nbt, interp, n_regions,
                                                  interval);
    return alea_nuc_interp_pair(x[interval], x[interval + 1], y[interval],
                                y[interval + 1], query, code, value);
}

bool alea_nuc_tabular_pdf_valid(const double* x, const double* pdf,
                                const double* cdf, int n, int interpolation,
                                int n_discrete) {
    if (!x || !pdf || !cdf || n < 1 ||
        (interpolation != 1 && interpolation != 2 &&
         !(interpolation == 0 && n_discrete == n)) ||
        n_discrete < 0 || n_discrete > n ||
        (n_discrete < n && n - n_discrete < 2)) return false;
    for (int i = 0; i < n; i++) {
        if (!isfinite(x[i]) || !isfinite(pdf[i]) || !isfinite(cdf[i]) ||
            pdf[i] < 0.0 || cdf[i] < -1e-8 || cdf[i] > 1.0 + 1e-6 ||
            (i > 0 && cdf[i] < cdf[i - 1])) return false;
        if (i > n_discrete && x[i] < x[i - 1] &&
            (cdf[i - 1] < 1.0 - 1e-8 || cdf[i] < 1.0 - 1e-8))
            return false;
        if (i > n_discrete && x[i] == x[i - 1] &&
            fabs(cdf[i] - cdf[i - 1]) > 1e-8) return false;
    }
    if (n_discrete == 0 && fabs(cdf[0]) > 1e-8) return false;
    return fabs(cdf[n - 1] - 1.0) <= 1e-6;
}

alea_error_t alea_nuc_tabular_pdf_sample(const double* x, const double* pdf,
                                         const double* cdf, int n,
                                         int interpolation, int n_discrete,
                                         double xi, double* value,
                                         int* sampled_index) {
    if (!value || !isfinite(xi) || xi < 0.0 || xi >= 1.0 ||
        !alea_nuc_tabular_pdf_valid(x, pdf, cdf, n, interpolation, n_discrete))
        return ALEA_ERR_INVALID_ARG;

    for (int i = 0; i < n_discrete; i++) {
        if (xi < cdf[i]) {
            *value = x[i];
            if (sampled_index) *sampled_index = i;
            return ALEA_OK;
        }
    }
    if (n_discrete == n) {
        *value = x[n - 1];
        if (sampled_index) *sampled_index = n - 1;
        return ALEA_OK;
    }

    int lo = n_discrete;
    while (lo + 1 < n && cdf[lo + 1] <= xi) lo++;
    if (lo + 1 >= n) lo = n - 2;
    double probability = xi - cdf[lo];
    double width = x[lo + 1] - x[lo];
    double offset = 0.0;
    if (probability > 0.0) {
        if (interpolation == 1 || fabs(pdf[lo + 1] - pdf[lo]) <= 1e-14) {
            if (pdf[lo] <= 0.0) return ALEA_ERR_INVALID_STATE;
            offset = probability / pdf[lo];
        } else {
            double slope = (pdf[lo + 1] - pdf[lo]) / width;
            double discriminant = pdf[lo] * pdf[lo] +
                                  2.0 * slope * probability;
            if (discriminant < 0.0) return ALEA_ERR_INVALID_STATE;
            double denominator = pdf[lo] + sqrt(discriminant);
            if (denominator <= 0.0) return ALEA_ERR_INVALID_STATE;
            offset = 2.0 * probability / denominator;
        }
    }
    *value = fmin(x[lo + 1], fmax(x[lo], x[lo] + offset));
    if (sampled_index) *sampled_index = lo;
    return ALEA_OK;
}
