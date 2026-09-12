// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/**
 * @file lookup.c
 * @brief Energy grid binary search and cross-section interpolation
 */

#include "nuclear_internal.h"
#include <math.h>

/* Binary search for grids already validated during decoding. */
int alea_nuc_energy_lookup_trusted(const double* energy, int n, double E,
                                   double* frac) {
    if (n < 2 || !energy || !isfinite(E)) return -1;
    if (E <= energy[0]) {
        if (frac) *frac = 0.0;
        return 0;
    }
    if (E >= energy[n - 1]) {
        if (frac) *frac = 1.0;
        return n - 2;
    }
    int lo = 0, hi = n - 2;
    while (lo < hi) {
        int mid = lo + (hi - lo + 1) / 2;
        if (energy[mid] <= E) lo = mid;
        else hi = mid - 1;
    }
    if (frac)
        *frac = (E - energy[lo]) / (energy[lo + 1] - energy[lo]);
    return lo;
}

int alea_nuc_energy_lookup(const double* energy, int n, double E, double* frac) {
    if (n < 2 || !energy || !isfinite(E)) return -1;

    for (int i = 0; i < n; i++) {
        if (!isfinite(energy[i])) return -1;
        if (i > 0 && !(energy[i] > energy[i - 1])) return -1;
    }

    return alea_nuc_energy_lookup_trusted(energy, n, E, frac);
}

/** Linear interpolation on a grid */
static double interp(const double* grid, const double* values, int n,
                     double E) {
    if (!grid || !values || n < 2) return 0.0;
    double f;
    int i = alea_nuc_energy_lookup_trusted(grid, n, E, &f);
    if (i < 0) return 0.0;
    return values[i] + f * (values[i + 1] - values[i]);
}

/* Photoatomic grids may repeat an energy to encode a right-continuous shell
 * threshold. The trusted lookup selects the last copy at the exact knot. */
static double interp_photon_values_ll(const alea_nuc_photon_data_t* photon,
                                      const double* values, double energy) {
    if (!photon || !photon->ln_energy || !values ||
        photon->n_energies < 2 || energy <= 0.0 || !isfinite(energy)) return 0.0;
    double fraction;
    int i = alea_nuc_energy_lookup_trusted(
        photon->ln_energy, photon->n_energies, log(energy), &fraction);
    if (i < 0) return 0.0;
    double v0 = values[i], v1 = values[i + 1];
    if (fraction <= 0.0) return v0;
    if (fraction >= 1.0) return v1;
    if (v0 <= 0.0 || v1 <= 0.0)
        return v0 + fraction * (v1 - v0);
    return exp(log(v0) + fraction * (log(v1) - log(v0)));
}

/** Log-log interpolation (common for cross sections) */
double alea_nuc_interp_loglog(const double* grid, const double* values, int n,
                         double E) {
    if (!grid || !values || n < 2) return 0.0;
    double f;
    int i = alea_nuc_energy_lookup(grid, n, E, &f);
    if (i < 0) return 0.0;

    double v0 = values[i];
    double v1 = values[i + 1];

    if (f <= 0.0) return v0;
    if (f >= 1.0) return v1;

    if (E <= 0.0)
        return v0 + f * (v1 - v0);

    /* Fall back to lin-lin if values are zero or negative */
    if (v0 <= 0.0 || v1 <= 0.0)
        return v0 + f * (v1 - v0);

    double e0 = grid[i];
    double e1 = grid[i + 1];
    if (e0 <= 0.0 || e1 <= 0.0)
        return v0 + f * (v1 - v0);

    /* log-log: log(σ) = log(σ0) + [log(E/E0)/log(E1/E0)] * log(σ1/σ0) */
    double log_ratio = log(E / e0) / log(e1 / e0);
    return v0 * pow(v1 / v0, log_ratio);
}

double alea_nuc_xs_total(const alea_nuc_nuclide_t* nuc, double energy) {
    if (!nuc) return 0.0;
    if (nuc->particle == ALEA_NUC_PARTICLE_PHOTON && nuc->photon) {
        return alea_nuc_photon_xs_incoherent(nuc, energy) +
               alea_nuc_photon_xs_coherent(nuc, energy) +
               alea_nuc_photon_xs_photoelectric(nuc, energy) +
               alea_nuc_photon_xs_pair(nuc, energy);
    }
    if (!nuc->sigma_total) return 0.0;
    return interp(nuc->energy, nuc->sigma_total, nuc->n_energies, energy);
}

double alea_nuc_xs_absorption(const alea_nuc_nuclide_t* nuc, double energy) {
    if (!nuc || !nuc->sigma_abs) return 0.0;
    return interp(nuc->energy, nuc->sigma_abs, nuc->n_energies, energy);
}

double alea_nuc_xs_elastic(const alea_nuc_nuclide_t* nuc, double energy) {
    if (!nuc || !nuc->sigma_elastic) return 0.0;
    return interp(nuc->energy, nuc->sigma_elastic, nuc->n_energies, energy);
}

double alea_nuc_xs_heating(const alea_nuc_nuclide_t* nuc, double energy) {
    if (!nuc) return 0.0;
    if (nuc->particle == ALEA_NUC_PARTICLE_PHOTON && nuc->photon) {
        if (!nuc->photon->heating) return 0.0;
        return interp_photon_values_ll(nuc->photon, nuc->photon->heating,
                                       energy);
    }
    if (!nuc->heating) return 0.0;
    return interp(nuc->energy, nuc->heating, nuc->n_energies, energy);
}

double alea_nuc_heating_per_collision(const alea_nuc_nuclide_t* nuc, double energy) {
    double h = alea_nuc_xs_heating(nuc, energy);
    double sig_t = alea_nuc_xs_total(nuc, energy);
    return (sig_t > 0.0) ? h / sig_t : 0.0;
}

/** Fast log-log interpolation using pre-stored ln values.
 *  Binary search on ln_energy grid, linear interp in log space, single exp(). */
static double interp_photon_ll(const double* ln_grid, const double* ln_values,
                                int n, double ln_E) {
    if (!ln_grid || !ln_values || n < 2 || !isfinite(ln_E)) return 0.0;
    double f;
    int i = alea_nuc_energy_lookup_trusted(ln_grid, n, ln_E, &f);
    if (i < 0) return 0.0;
    double ln_v0 = ln_values[i];
    double ln_v1 = ln_values[i + 1];
    if (f <= 0.0) return ln_v0 <= -1e30 ? 0.0 : exp(ln_v0);
    if (f >= 1.0) return ln_v1 <= -1e30 ? 0.0 : exp(ln_v1);
    if (ln_v0 <= -1e30 || ln_v1 <= -1e30) return 0.0;
    return exp(ln_v0 + f * (ln_v1 - ln_v0));
}

double alea_nuc_photon_xs_incoherent(const alea_nuc_nuclide_t* nuc, double energy) {
    if (!nuc || !nuc->photon) return 0.0;
    if (energy <= 0.0 || !isfinite(energy)) return 0.0;
    const alea_nuc_photon_data_t* ph = nuc->photon;
    if (!ph->ln_energy || !ph->ln_sigma_incoherent) return 0.0;
    return interp_photon_ll(ph->ln_energy, ph->ln_sigma_incoherent,
                             ph->n_energies, log(energy));
}

double alea_nuc_photon_xs_coherent(const alea_nuc_nuclide_t* nuc, double energy) {
    if (!nuc || !nuc->photon) return 0.0;
    if (energy <= 0.0 || !isfinite(energy)) return 0.0;
    const alea_nuc_photon_data_t* ph = nuc->photon;
    if (!ph->ln_energy || !ph->ln_sigma_coherent) return 0.0;
    return interp_photon_ll(ph->ln_energy, ph->ln_sigma_coherent,
                             ph->n_energies, log(energy));
}

double alea_nuc_photon_xs_photoelectric(const alea_nuc_nuclide_t* nuc, double energy) {
    if (!nuc || !nuc->photon) return 0.0;
    if (energy <= 0.0 || !isfinite(energy)) return 0.0;
    const alea_nuc_photon_data_t* ph = nuc->photon;
    if (!ph->ln_energy || !ph->ln_sigma_photoelectric) return 0.0;
    return interp_photon_ll(ph->ln_energy, ph->ln_sigma_photoelectric,
                             ph->n_energies, log(energy));
}

double alea_nuc_photon_xs_photoelectric_subshell(
    const alea_nuc_nuclide_t* nuc, int designator, double energy) {
    if (!nuc || !nuc->photon || designator <= 0 ||
        energy <= 0.0 || !isfinite(energy)) return 0.0;
    const alea_nuc_photon_data_t* ph = nuc->photon;
    if (!ph->ln_energy || ph->n_subshells < 0 ||
        (ph->n_subshells > 0 && !ph->subshells)) return 0.0;
    for (int i = 0; i < ph->n_subshells; i++) {
        const alea_nuc_atomic_subshell_t* shell = &ph->subshells[i];
        if (shell->designator != designator || !shell->ln_photoelectric_xs)
            continue;
        return interp_photon_ll(ph->ln_energy, shell->ln_photoelectric_xs,
                                ph->n_energies, log(energy));
    }
    return 0.0;
}

double alea_nuc_photon_xs_pair(const alea_nuc_nuclide_t* nuc, double energy) {
    if (!nuc || !nuc->photon) return 0.0;
    if (energy <= 0.0 || !isfinite(energy)) return 0.0;
    const alea_nuc_photon_data_t* ph = nuc->photon;
    if (!ph->ln_energy || !ph->ln_sigma_pair) return 0.0;
    return interp_photon_ll(ph->ln_energy, ph->ln_sigma_pair,
                             ph->n_energies, log(energy));
}

/** Interpolate reaction XS at a pre-computed energy grid index */
static inline double reaction_xs_at(const alea_nuc_reaction_t* r, int ie, double f) {
    if (r->n_energies <= 0 || !r->xs) return 0.0;
    int ie_start = r->threshold_index - 1;
    int ri = ie - ie_start;
    if (ri < 0) return 0.0;
    if (ri >= r->n_energies - 1) return r->xs[r->n_energies - 1];
    return r->xs[ri] + f * (r->xs[ri + 1] - r->xs[ri]);
}

double alea_nuc_xs_reaction(const alea_nuc_nuclide_t* nuc, int mt, double energy) {
    if (!nuc) return 0.0;
    if (mt == 1) return alea_nuc_xs_total(nuc, energy);
    if (mt == 2) return alea_nuc_xs_elastic(nuc, energy);
    if (mt == 3) {
        double value = alea_nuc_xs_total(nuc, energy) -
                       alea_nuc_xs_elastic(nuc, energy);
        return value > 0.0 ? value : 0.0;
    }
    if (mt == 27 || mt == 101) return alea_nuc_xs_absorption(nuc, energy);

    /* O(1) lookup via MT table */
    const alea_nuc_reaction_t* r = NULL;
    if (nuc->mt_to_rxn && mt >= 0 && mt < ALEA_NUC_MT_TABLE_SIZE) {
        int idx = nuc->mt_to_rxn[mt];
        if (idx >= 0) r = &nuc->reactions[idx];
    } else {
        for (int i = 0; i < nuc->n_reactions; i++) {
            if (nuc->reactions[i].mt == mt) { r = &nuc->reactions[i]; break; }
        }
    }

    /* Some ACE photon-production channels reference redundant aggregate
     * reactions that are omitted from MTR. Reconstruct the two aggregates
     * used by conventional neutron tables from their sampled children. */
    if (!r && (mt == 4 || mt == 18)) {
        double sum = 0.0;
        for (int i = 0; i < nuc->n_reactions; i++) {
            int child = nuc->reactions[i].mt;
            if ((mt == 4 && child >= 51 && child <= 91) ||
                (mt == 18 && (child == 19 || child == 20 || child == 21 ||
                              child == 38)))
                sum += alea_nuc_xs_reaction(nuc, child, energy);
        }
        return sum;
    }
    if (!r) return 0.0;

    double f;
    int ie = alea_nuc_energy_lookup_trusted(nuc->energy, nuc->n_energies, energy, &f);
    if (ie < 0) return 0.0;
    return reaction_xs_at(r, ie, f);
}

int alea_nuc_urr_factors(const alea_nuc_nuclide_t* nuc, double energy, double xi,
                     double factors[5]) {
    if (!nuc || !nuc->urr || !factors || !isfinite(energy) ||
        !isfinite(xi) || xi < 0.0 || xi >= 1.0) return 0;

    const alea_nuc_urr_t* urr = nuc->urr;
    int N = urr->n_energies;
    int M = urr->n_bands;
    if (N < 2 || M <= 0 || !urr->energy || !urr->table) return 0;

    /* Check if energy is in URR range */
    if (energy < urr->energy[0] || energy > urr->energy[N - 1])
        return 0;

    /* Find energy bracket */
    double f;
    int ie = alea_nuc_energy_lookup(urr->energy, N, energy, &f);
    if (ie < 0) return 0;

    /* For each energy, table layout: [M cum_prob] [M σ_t] [M σ_el] [M σ_f] [M σ_c] [M heat] */
    /* Find probability band using xi on lower energy */
    const double* tab_lo = &urr->table[ie * 6 * M];
    int band_lo = M - 1;
    for (int j = 0; j < M; j++) {
        if (tab_lo[j] >= xi) { band_lo = j; break; }
    }

    /* Extract factors at lower energy */
    double f_lo[5];
    for (int q = 0; q < 5; q++)
        f_lo[q] = tab_lo[(q + 1) * M + band_lo];

    if (f <= 0.0 || ie >= N - 1) {
        /* At or below first energy, use lower bracket directly */
        for (int q = 0; q < 5; q++)
            factors[q] = f_lo[q];
        goto normalize;
    }

    /* Find probability band at upper energy */
    const double* tab_hi = &urr->table[(ie + 1) * 6 * M];
    int band_hi = M - 1;
    for (int j = 0; j < M; j++) {
        if (tab_hi[j] >= xi) { band_hi = j; break; }
    }

    double f_hi[5];
    for (int q = 0; q < 5; q++)
        f_hi[q] = tab_hi[(q + 1) * M + band_hi];

    /* Interpolate between energies */
    if (urr->interp == 5) {
        /* Log-log interpolation */
        if (energy > 0.0 && urr->energy[ie] > 0.0 && urr->energy[ie + 1] > 0.0)
            f = log(energy / urr->energy[ie]) /
                log(urr->energy[ie + 1] / urr->energy[ie]);
        for (int q = 0; q < 5; q++) {
            if (f_lo[q] > 0.0 && f_hi[q] > 0.0)
                factors[q] = f_lo[q] * pow(f_hi[q] / f_lo[q], f);
            else
                factors[q] = f_lo[q] + f * (f_hi[q] - f_lo[q]);
        }
    } else {
        /* Lin-lin interpolation */
        for (int q = 0; q < 5; q++)
            factors[q] = f_lo[q] + f * (f_hi[q] - f_lo[q]);
    }

normalize:
    /* Some processed probability tables contain small negative capture
     * values in their highest probability bands. A negative reaction rate
     * cannot enter transport selection, so clamp collision cross sections;
     * keep heating signed because negative kerma is meaningful. */
    for (int q = 0; q < 4; q++)
        if (factors[q] < 0.0) factors[q] = 0.0;
    if (!urr->multiply_smooth) {
        double smooth[5];
        smooth[0] = alea_nuc_xs_total(nuc, energy);
        smooth[1] = alea_nuc_xs_elastic(nuc, energy);
        smooth[2] = alea_nuc_xs_reaction(nuc, 18, energy);
        if (smooth[2] <= 0.0) {
            for (int r = 0; r < nuc->n_reactions; r++) {
                int mt = nuc->reactions[r].mt;
                if (mt == 19 || mt == 20 || mt == 21 || mt == 38)
                    smooth[2] += alea_nuc_xs_reaction(nuc, mt, energy);
            }
        }
        smooth[3] = alea_nuc_xs_reaction(nuc, 102, energy);
        smooth[4] = alea_nuc_xs_heating(nuc, energy);
        for (int q = 0; q < 5; q++)
            factors[q] = smooth[q] > 0.0 ? factors[q] / smooth[q] : 0.0;
    }

    return 1;
}
