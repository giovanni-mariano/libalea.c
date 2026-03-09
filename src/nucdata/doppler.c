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
 *   σ_D(E) = 1/(y√π) ∫₀^∞ y'·σ(E')·[exp(-(y'-y)²) - exp(-(y'+y)²)] dy'
 *
 * where y = √(α·E), α = AWR / ΔkT, ΔkT = kT_new - kT_old.
 *
 * The second exponential (y'+y)² is negligible for y > ~4 (above ~1 eV
 * for most nuclides at room temperature).
 *
 * Uses trapezoidal quadrature on the existing ACE energy grid
 * transformed to y-space. Integration window is ±4σ of the Gaussian
 * kernel, found via binary search for O(N·W) complexity where W is
 * the local window width.
 *
 * Strategy for consistency:
 * - Broaden σ_total and σ_elastic independently on the full grid.
 * - Re-derive σ_abs = σ_total - σ_elastic after broadening to maintain
 *   the identity σ_total = σ_elastic + σ_abs.
 * - Broaden per-reaction XS on the full grid (padded with zeros below
 *   threshold) to properly capture the kernel tail near threshold.
 */

#include "alea_nucdata.h"
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

        /* At very low y (< 1e-6), the Gaussian kernel is so wide that
         * all grid points contribute equally — σ_D ≈ σ_0. This threshold
         * corresponds to E < ~1e-12/α MeV, well below any resonance
         * structure. */
        if (yj < 1e-6) {
            buf[j] = sigma_in[j];
            continue;
        }

        /* Find integration bounds in index space: y_grid[i] ∈ [yj-4, yj+4] */
        double y_lo = yj - ALEA_NUC_KERNEL_HALFWIDTH;
        if (y_lo < 0.0) y_lo = 0.0;
        double y_hi = yj + ALEA_NUC_KERNEL_HALFWIDTH;

        /* Binary search for lower bound */
        int lo = 0, hi = n - 1;
        while (lo < hi) {
            int mid = (lo + hi) / 2;
            if (y_grid[mid] < y_lo)
                lo = mid + 1;
            else
                hi = mid;
        }
        int i_lo = lo;

        /* Binary search for upper bound */
        lo = i_lo;
        hi = n - 1;
        while (lo < hi) {
            int mid = (lo + hi + 1) / 2;
            if (y_grid[mid] > y_hi)
                hi = mid - 1;
            else
                lo = mid;
        }
        int i_hi = lo;

        /* Trapezoidal integration in y-space */
        double sum = 0.0;
        for (int i = i_lo; i <= i_hi; i++) {
            double yi = y_grid[i];
            double d1 = yi - yj;
            double d2 = yi + yj;
            double kernel = yi * (exp(-d1 * d1) - exp(-d2 * d2));

            /* Trapezoidal weight */
            double dy;
            if (i == 0)
                dy = (y_grid[1] - y_grid[0]);
            else if (i == n - 1)
                dy = (y_grid[n - 1] - y_grid[n - 2]);
            else
                dy = (y_grid[i + 1] - y_grid[i - 1]) / 2.0;

            sum += sigma_in[i] * kernel * dy;
        }

        buf[j] = sum / (yj * sqrt(M_PI));

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
    if (dkT <= 0.0) return ALEA_ERR_UNSUPPORTED;

    int n = nuc->n_energies;
    if (n <= 1 || !nuc->energy) return ALEA_ERR_INVALID_ARG;

    /* α = AWR / ΔkT — controls the kernel width */
    double alpha = nuc->awr / dkT;

    /* Precompute y = √(α·E) for each grid point */
    double* y_grid = malloc((size_t)n * sizeof(double));
    double* work = malloc((size_t)n * sizeof(double));
    if (!y_grid || !work) { free(y_grid); free(work); return ALEA_ERR_OUT_OF_MEMORY; }

    for (int i = 0; i < n; i++)
        y_grid[i] = sqrt(alpha * nuc->energy[i]);

    /* Broaden total and elastic on the full grid */
    if (nuc->sigma_total)
        broaden_array(y_grid, nuc->sigma_total, nuc->sigma_total, n, work);
    if (nuc->sigma_elastic)
        broaden_array(y_grid, nuc->sigma_elastic, nuc->sigma_elastic, n, work);
    if (nuc->heating)
        broaden_array(y_grid, nuc->heating, nuc->heating, n, work);

    /* Re-derive absorption = total - elastic to maintain consistency.
     * This ensures σ_total = σ_elastic + σ_abs exactly after broadening,
     * rather than broadening σ_abs independently which breaks the identity. */
    if (nuc->sigma_abs && nuc->sigma_total && nuc->sigma_elastic) {
        for (int i = 0; i < n; i++) {
            nuc->sigma_abs[i] = nuc->sigma_total[i] - nuc->sigma_elastic[i];
            if (nuc->sigma_abs[i] < 0.0) nuc->sigma_abs[i] = 0.0;
        }
    }

    /* Broaden per-reaction cross sections on the FULL energy grid.
     * Embed the reaction XS (which starts at threshold_index) into a
     * full-grid array padded with zeros, broaden, then extract back.
     * This properly captures the kernel tail reaching into the zero
     * region below threshold. */
    double* full_xs = calloc((size_t)n, sizeof(double));
    if (full_xs) {
        for (int r = 0; r < nuc->n_reactions; r++) {
            alea_nuc_reaction_t* rxn = &nuc->reactions[r];
            if (!rxn->xs || rxn->n_energies <= 1) continue;

            int ie_start = rxn->threshold_index - 1; /* 0-based */
            int nr = rxn->n_energies;

            /* Embed into full grid (zeros below threshold) */
            memset(full_xs, 0, (size_t)n * sizeof(double));
            for (int i = 0; i < nr && ie_start + i < n; i++)
                full_xs[ie_start + i] = rxn->xs[i];

            /* Broaden on full grid */
            broaden_array(y_grid, full_xs, full_xs, n, work);

            /* Extract back into reaction sub-grid */
            for (int i = 0; i < nr && ie_start + i < n; i++)
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
