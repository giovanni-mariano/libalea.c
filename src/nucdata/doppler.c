// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/**
 * @file doppler.c
 * @brief On-the-fly Doppler broadening of pointwise cross sections
 *
 * Broadens cross sections from temperature T₀ to T using the exact
 * kernel method:
 *
 *   σ_D(E) = 1/(y²√π) ∫₀^∞ y'²·σ(E')·[exp(-(y'-y)²) - exp(-(y'+y)²)] dy'
 *
 * where y = √(α·E), α = AWR / ΔkT, ΔkT = kT_new - kT_old.
 *
 * The second exponential (y'+y)² is negligible for y > ~4 (above ~1 eV
 * for most nuclides at room temperature).
 *
 * Uses composite Simpson quadrature in transformed y-space. Sampling the
 * existing ACE knots directly is unstable when the Doppler kernel is much
 * narrower than a grid interval, so sigma is interpolated between knots on a
 * quadrature mesh that resolves the Gaussian.
 *
 * Strategy for consistency:
 * - Broaden σ_total, σ_absorption, and σ_elastic independently.
 * - Broaden per-reaction XS on the full grid (padded with zeros below
 *   threshold) to properly capture the kernel tail near threshold.
 */

#include "nuclear_internal.h"
#include "constants.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "util/math.h"

/**
 * Doppler-broaden a single cross-section array.
 *
 * @param y_grid   Precomputed y = √(α·E) for each grid point
 * @param sigma_in Input cross sections (barns)
 * @param sigma_out Output broadened cross sections (barns), may alias sigma_in
 * @param n        Number of grid points
 * @param buf      Work buffer of size n
 */
static void broaden_array(const double* y_grid,
                          const double* sigma_in, double* sigma_out, int n,
                          double* buf) {

    for (int j = 0; j < n; j++) {
        double yj = y_grid[j];
        if (yj < 1e-6) {
            buf[j] = sigma_in[j];
            continue;
        }

        double x_lo = yj > 8.0 ? yj - 8.0 : 0.0;
        double x_hi = yj + 8.0;
        int panels = (int)ceil((x_hi - x_lo) / 0.125);
        if (panels < 16) panels = 16;
        if (panels & 1) panels++;
        double h = (x_hi - x_lo) / panels;
        double sum = 0.0;

        for (int k = 0; k <= panels; k++) {
            double x = x_lo + k * h;
            double sigma;
            if (x <= y_grid[0]) {
                sigma = sigma_in[0];
            } else if (x >= y_grid[n - 1]) {
                sigma = sigma_in[n - 1];
            } else {
                int lo = 0, hi = n - 2;
                while (lo < hi) {
                    int mid = lo + (hi - lo + 1) / 2;
                    if (y_grid[mid] <= x) lo = mid;
                    else hi = mid - 1;
                }
                /* ACE cross sections are linear in energy, and E is
                 * proportional to y squared. */
                double f = (x * x - y_grid[lo] * y_grid[lo]) /
                           (y_grid[lo + 1] * y_grid[lo + 1] -
                            y_grid[lo] * y_grid[lo]);
                sigma = sigma_in[lo] + f * (sigma_in[lo + 1] - sigma_in[lo]);
            }

            double d = x - yj;
            double kernel_diff = exp(-d * d) * -expm1(-4.0 * x * yj);
            double integrand = x * x * sigma * kernel_diff;
            int weight = (k == 0 || k == panels) ? 1 : (k & 1 ? 4 : 2);
            sum += weight * integrand;
        }

        buf[j] = (h / 3.0) * sum / (yj * yj * sqrt(M_PI));

        /* Ensure non-negative */
        if (buf[j] < 0.0) buf[j] = 0.0;
    }

    /* Copy result */
    for (int j = 0; j < n; j++)
        sigma_out[j] = buf[j];
}

alea_error_t alea_nuc_doppler_broaden(alea_nuc_nuclide_t* nuc, double kT_target) {
    if (!nuc) return ALEA_ERR_NULL_ARG;

    double kT_old = nuc->temperature;
    double dkT = kT_target - kT_old;

    /* Can only broaden to higher temperature */
    if (!isfinite(kT_target) || dkT <= 0.0) return ALEA_ERR_UNSUPPORTED;

    /* Probability tables are processed for the source table temperature.
     * Broadening only the smooth pointwise data would create an internally
     * inconsistent transport table, so reject before allocation or mutation. */
    if (nuc->urr) return ALEA_ERR_UNSUPPORTED;

    int n = nuc->n_energies;
    if (n <= 1 || !nuc->energy || !isfinite(nuc->awr) || nuc->awr <= 0.0)
        return ALEA_ERR_INVALID_ARG;
    for (int i = 0; i < n; i++) {
        if (!isfinite(nuc->energy[i]) || nuc->energy[i] < 0.0 ||
            (i > 0 && nuc->energy[i] <= nuc->energy[i - 1]))
            return ALEA_ERR_INVALID_ARG;
    }
    for (int r = 0; r < nuc->n_reactions; r++) {
        const alea_nuc_reaction_t* rxn = &nuc->reactions[r];
        if (!rxn->xs) continue;
        int ie_start = rxn->threshold_index - 1;
        if (ie_start < 0 || ie_start >= n || rxn->n_energies <= 0 ||
            rxn->n_energies > n - ie_start)
            return ALEA_ERR_INVALID_ARG;
    }

    /* α = AWR / ΔkT — controls the kernel width */
    double alpha = nuc->awr / dkT;
    if (!isfinite(alpha) || alpha <= 0.0) return ALEA_ERR_INVALID_ARG;

    /* Precompute y = √(α·E) for each grid point */
    double* y_grid = alea_nuc_malloc((size_t)n * sizeof(double));
    double* work = alea_nuc_malloc((size_t)n * sizeof(double));
    double* full_xs = nuc->n_reactions > 0 ? alea_nuc_calloc((size_t)n, sizeof(double)) : NULL;
    if (!y_grid || !work || (nuc->n_reactions > 0 && !full_xs)) {
        free(y_grid); free(work); free(full_xs);
        return ALEA_ERR_OUT_OF_MEMORY;
    }

    for (int i = 0; i < n; i++)
        y_grid[i] = sqrt(alpha * nuc->energy[i]);

    /* Broaden total and elastic on the full grid */
    if (nuc->sigma_total)
        broaden_array(y_grid, nuc->sigma_total, nuc->sigma_total, n, work);
    if (nuc->sigma_elastic)
        broaden_array(y_grid, nuc->sigma_elastic, nuc->sigma_elastic, n, work);
    if (nuc->sigma_abs)
        broaden_array(y_grid, nuc->sigma_abs, nuc->sigma_abs, n, work);
    if (nuc->heating)
        broaden_array(y_grid, nuc->heating, nuc->heating, n, work);

    /* Broaden per-reaction cross sections on the FULL energy grid.
     * Embed the reaction XS (which starts at threshold_index) into a
     * full-grid array padded with zeros, broaden, then extract back.
     * This properly captures the kernel tail reaching into the zero
     * region below threshold. */
    if (full_xs) {
        for (int r = 0; r < nuc->n_reactions; r++) {
            alea_nuc_reaction_t* rxn = &nuc->reactions[r];
            if (!rxn->xs || rxn->n_energies <= 1) continue;

            int ie_start = rxn->threshold_index - 1; /* 0-based */
            int nr = rxn->n_energies;
            /* Embed into full grid (zeros below threshold) */
            memset(full_xs, 0, (size_t)n * sizeof(double));
            for (int i = 0; i < nr; i++)
                full_xs[ie_start + i] = rxn->xs[i];

            /* Broaden on full grid */
            broaden_array(y_grid, full_xs, full_xs, n, work);

            /* Extract back into reaction sub-grid */
            for (int i = 0; i < nr; i++)
                rxn->xs[i] = full_xs[ie_start + i];
        }
        free(full_xs);
    }

    free(y_grid);
    free(work);

    /* Update nuclide temperature */
    nuc->temperature = kT_target;

    return ALEA_OK;
}
