// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/**
 * @file multigroup.c
 * @brief Multigroup cross sections and adjoint scattering matrix
 *
 * Collapses continuous-energy cross sections into group constants
 * using a weighting spectrum (default: 1/E in slowing-down range,
 * Maxwellian below 0.625 eV, fission spectrum above 0.1 MeV).
 *
 * Group convention: group 0 = highest energy, group G-1 = lowest.
 * Group boundaries are in descending order: bounds[0] > bounds[1] > ... > bounds[G].
 *
 * The scattering matrix scatter[g_from * n_groups + g_to] gives the
 * group transfer cross section from group g_from to g_to,
 * including elastic, inelastic, and (n,Xn) contributions.
 * The adjoint matrix is the transpose.
 */

#include "alea_nucdata.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

alea_nuc_multigroup_t* alea_nuc_mg_create(int n_groups, const double* bounds) {
    if (n_groups <= 0 || !bounds) return NULL;

    for (int i = 0; i <= n_groups; i++) {
        if (!isfinite(bounds[i]) || bounds[i] < 0.0) return NULL;
        if (i > 0 && !(bounds[i - 1] > bounds[i])) return NULL;
    }

    size_t groups = (size_t)n_groups;
    if (groups > SIZE_MAX / groups) return NULL;

    alea_nuc_multigroup_t* mg = calloc(1, sizeof(*mg));
    if (!mg) return NULL;

    mg->n_groups = n_groups;
    mg->bounds = malloc((groups + 1) * sizeof(double));
    mg->sigma_t = calloc(groups, sizeof(double));
    mg->sigma_a = calloc(groups, sizeof(double));
    mg->sigma_s = calloc(groups, sizeof(double));
    mg->sigma_f = calloc(groups, sizeof(double));
    mg->nu_sigma_f = calloc(groups, sizeof(double));
    mg->chi = calloc(groups, sizeof(double));
    mg->scatter = calloc(groups * groups, sizeof(double));

    if (!mg->bounds || !mg->sigma_t || !mg->sigma_a ||
        !mg->sigma_s || !mg->sigma_f || !mg->nu_sigma_f ||
        !mg->chi || !mg->scatter) {
        alea_nuc_mg_destroy(mg);
        return NULL;
    }

    for (int i = 0; i <= n_groups; i++)
        mg->bounds[i] = bounds[i];

    return mg;
}

void alea_nuc_mg_destroy(alea_nuc_multigroup_t* mg) {
    if (!mg) return;
    free(mg->bounds);
    free(mg->sigma_t);
    free(mg->sigma_a);
    free(mg->sigma_s);
    free(mg->sigma_f);
    free(mg->nu_sigma_f);
    free(mg->chi);
    free(mg->scatter);
    free(mg);
}

/* ============================================================================
 * Default weighting spectrum: Maxwellian + 1/E + fission
 *
 *   φ(E) = C₁ · E · exp(-E/kT)     for E < E_thermal   (0.625 eV)
 *   φ(E) = C₂ / E                   for E_thermal ≤ E ≤ E_fission (0.1 MeV)
 *   φ(E) = C₃ · exp(-E/a) · sinh(√(bE))  for E > E_fission
 *
 * The constants C₁,C₂,C₃ ensure continuity. For group averaging we only
 * need the shape within each group, so absolute normalization cancels.
 * ============================================================================ */

#define E_THERMAL  0.625e-6   /* 0.625 eV in MeV */
#define E_FISSION  0.1        /* 100 keV in MeV */
#define KT_DEFAULT 2.53e-8    /* 293.6 K in MeV */
#define WATT_A_DEFAULT 0.988  /* U-235 thermal Watt parameter */
#define WATT_B_DEFAULT 2.249

static double default_weight(double E, void* ctx) {
    (void)ctx;
    if (E <= 0.0) return 0.0;
    if (E < E_THERMAL) {
        /* Maxwellian: E·exp(-E/kT) */
        return E * exp(-E / KT_DEFAULT);
    } else if (E < E_FISSION) {
        /* 1/E spectrum, continuous with Maxwellian at E_THERMAL:
         * Maxwellian at E_THERMAL = E_T · exp(-E_T/kT)
         * 1/E at E_THERMAL = C/E_T => C = E_T² · exp(-E_T/kT) */
        double C = E_THERMAL * E_THERMAL * exp(-E_THERMAL / KT_DEFAULT);
        return C / E;
    } else {
        /* Fission (Watt) spectrum, continuous with 1/E at E_FISSION:
         * 1/E at E_FISSION = C/E_F */
        double C = E_THERMAL * E_THERMAL * exp(-E_THERMAL / KT_DEFAULT);
        double val_1overE = C / E_FISSION;
        double watt = exp(-E / WATT_A_DEFAULT) * sinh(sqrt(WATT_B_DEFAULT * E));
        double watt_ref = exp(-E_FISSION / WATT_A_DEFAULT) *
                          sinh(sqrt(WATT_B_DEFAULT * E_FISSION));
        return (watt_ref > 0.0) ? val_1overE * watt / watt_ref : 0.0;
    }
}

/** Evaluate the weighting spectrum, using default if none set. */
static inline double weight(const alea_nuc_multigroup_t* mg, double E) {
    if (mg->spectrum)
        return mg->spectrum(E, mg->spectrum_ctx);
    return default_weight(E, NULL);
}

void alea_nuc_mg_set_spectrum(alea_nuc_multigroup_t* mg, alea_nuc_spectrum_fn fn,
                              void* ctx) {
    if (!mg) return;
    mg->spectrum = fn;
    mg->spectrum_ctx = ctx;
}

/**
 * Integrate σ(E) · φ(E) dE over [E_lo, E_hi] on the nuclide grid.
 *
 * Uses trapezoidal rule on the pointwise grid segments that overlap [E_lo, E_hi].
 * Both numerator ∫ σ·φ dE and denominator ∫ φ dE use the same quadrature
 * for consistency.
 *
 * Returns the flux-weighted group average: ∫ σ·φ dE / ∫ φ dE.
 */
static double group_average_xs(const alea_nuc_multigroup_t* mg,
                                const double* egrid, const double* xs,
                                int n, double E_lo, double E_hi) {
    if (E_hi <= E_lo || n <= 0) return 0.0;

    /* Find starting index */
    int i_lo = 0;
    while (i_lo < n - 1 && egrid[i_lo + 1] < E_lo) i_lo++;

    double sum_xs = 0.0;   /* ∫ σ(E)·φ(E) dE */
    double sum_flux = 0.0; /* ∫ φ(E) dE */

    for (int i = i_lo; i < n - 1; i++) {
        if (egrid[i] >= E_hi) break;

        double e0 = egrid[i];
        double e1 = egrid[i + 1];

        /* Clamp to group bounds */
        if (e0 < E_lo) e0 = E_lo;
        if (e1 > E_hi) e1 = E_hi;
        if (e1 <= e0) continue;

        /* Interpolate XS at segment endpoints */
        double dE_grid = egrid[i + 1] - egrid[i];
        double f0 = (dE_grid > 0.0) ? (e0 - egrid[i]) / dE_grid : 0.0;
        double f1 = (dE_grid > 0.0) ? (e1 - egrid[i]) / dE_grid : 0.0;
        double s0 = xs[i] + f0 * (xs[i + 1] - xs[i]);
        double s1 = xs[i] + f1 * (xs[i + 1] - xs[i]);

        double phi0 = weight(mg, e0);
        double phi1 = weight(mg, e1);

        /* Trapezoidal: ∫ σ·φ dE ≈ (σ₀·φ₀ + σ₁·φ₁)/2 · ΔE */
        double de = e1 - e0;
        sum_xs += (s0 * phi0 + s1 * phi1) / 2.0 * de;
        sum_flux += (phi0 + phi1) / 2.0 * de;
    }

    return (sum_flux > 0.0) ? sum_xs / sum_flux : 0.0;
}

/**
 * Compute ∫ ν(E) · σ_f(E) · φ(E) dE / ∫ φ(E) dE over a group.
 *
 * Unlike the single-midpoint approach, this properly integrates the
 * energy-dependent ν̄(E) together with σ_f(E).
 */
static double group_average_nu_sigma_f(const alea_nuc_multigroup_t* mg,
                                        const alea_nuc_nuclide_t* nuc,
                                        const double* fission_xs,
                                        double E_lo, double E_hi) {
    if (E_hi <= E_lo || !fission_xs) return 0.0;

    const double* egrid = nuc->energy;
    int n = nuc->n_energies;

    int i_lo = 0;
    while (i_lo < n - 1 && egrid[i_lo + 1] < E_lo) i_lo++;

    double sum_nusf = 0.0;
    double sum_flux = 0.0;

    for (int i = i_lo; i < n - 1; i++) {
        if (egrid[i] >= E_hi) break;

        double e0 = egrid[i];
        double e1 = egrid[i + 1];
        if (e0 < E_lo) e0 = E_lo;
        if (e1 > E_hi) e1 = E_hi;
        if (e1 <= e0) continue;

        double dE_grid = egrid[i + 1] - egrid[i];
        double f0 = (dE_grid > 0.0) ? (e0 - egrid[i]) / dE_grid : 0.0;
        double f1 = (dE_grid > 0.0) ? (e1 - egrid[i]) / dE_grid : 0.0;
        double sf0 = fission_xs[i] + f0 * (fission_xs[i + 1] - fission_xs[i]);
        double sf1 = fission_xs[i] + f1 * (fission_xs[i + 1] - fission_xs[i]);

        double nu0 = alea_nuc_nu_bar(nuc, e0);
        double nu1 = alea_nuc_nu_bar(nuc, e1);

        double phi0 = weight(mg, e0);
        double phi1 = weight(mg, e1);

        double de = e1 - e0;
        sum_nusf += (nu0 * sf0 * phi0 + nu1 * sf1 * phi1) / 2.0 * de;
        sum_flux += (phi0 + phi1) / 2.0 * de;
    }

    return (sum_flux > 0.0) ? sum_nusf / sum_flux : 0.0;
}

/**
 * Find which group an energy belongs to.
 * Returns -1 if outside the group structure.
 */
static int find_group(const alea_nuc_multigroup_t* mg, double E) {
    if (E >= mg->bounds[0] || E <= mg->bounds[mg->n_groups]) return -1;
    for (int g = 0; g < mg->n_groups; g++) {
        if (E <= mg->bounds[g] && E > mg->bounds[g + 1])
            return g;
    }
    return -1;
}

static int is_fission_mt(int mt) {
    return mt == 18 || mt == 19 || mt == 20 || mt == 21 || mt == 38;
}

/**
 * Build elastic scattering transfer matrix by integrating over
 * incident energies within each group.
 *
 * For elastic scattering off a nucleus with mass ratio A, the outgoing
 * energy range is [αE, E] where α = ((A-1)/(A+1))².
 * Assuming isotropic CM, the transfer probability is uniform in
 * the interval [αE, E].
 *
 * We integrate over the pointwise energy grid within each source group:
 *   σ_s(g→g') = ∫_g σ_el(E) · P(E→g') · φ(E) dE / ∫_g φ(E) dE
 *
 * where P(E→g') = overlap([αE,E], [E'_lo,E'_hi]) / ((1-α)E).
 */
static void build_elastic_scatter(alea_nuc_multigroup_t* mg,
                                   const alea_nuc_nuclide_t* nuc) {
    int G = mg->n_groups;
    double A = nuc->awr;
    double alpha = (A - 1.0) * (A - 1.0) / ((A + 1.0) * (A + 1.0));

    const double* egrid = nuc->energy;
    const double* sigma_el = nuc->sigma_elastic;
    int n = nuc->n_energies;
    if (!sigma_el || n <= 1) return;

    for (int gf = 0; gf < G; gf++) {
        double Eg_hi = mg->bounds[gf];
        double Eg_lo = mg->bounds[gf + 1];

        /* Integrate over pointwise grid within source group */
        int i_lo = 0;
        while (i_lo < n - 1 && egrid[i_lo + 1] < Eg_lo) i_lo++;

        double sum_flux = 0.0;

        for (int i = i_lo; i < n - 1; i++) {
            if (egrid[i] >= Eg_hi) break;

            double e0 = egrid[i];
            double e1 = egrid[i + 1];
            if (e0 < Eg_lo) e0 = Eg_lo;
            if (e1 > Eg_hi) e1 = Eg_hi;
            if (e1 <= e0) continue;

            /* Use midpoint of this sub-interval for the transfer calc */
            double E_mid = 0.5 * (e0 + e1);
            double dE_grid = egrid[i + 1] - egrid[i];
            double f_mid = (dE_grid > 0.0) ? (E_mid - egrid[i]) / dE_grid : 0.0;
            double sig = sigma_el[i] + f_mid * (sigma_el[i + 1] - sigma_el[i]);
            double phi = weight(mg, E_mid);
            double de = e1 - e0;
            double w = sig * phi * de;

            sum_flux += phi * de;

            if (w <= 0.0) continue;

            /* Outgoing energy range for elastic at E_mid */
            double Eout_lo = alpha * E_mid;
            double Eout_hi = E_mid;
            double Eout_range = Eout_hi - Eout_lo;

            if (Eout_range <= 0.0) {
                /* Heavy target: no energy loss */
                mg->scatter[gf * G + gf] += w;
                continue;
            }

            /* Distribute into target groups */
            for (int gt = 0; gt < G; gt++) {
                double Et_hi = mg->bounds[gt];
                double Et_lo = mg->bounds[gt + 1];

                double ov_lo = (Eout_lo > Et_lo) ? Eout_lo : Et_lo;
                double ov_hi = (Eout_hi < Et_hi) ? Eout_hi : Et_hi;

                if (ov_hi > ov_lo) {
                    double frac = (ov_hi - ov_lo) / Eout_range;
                    mg->scatter[gf * G + gt] += w * frac;
                }
            }
        }

        /* Normalize: divide by ∫ φ dE to get group-averaged transfer XS */
        if (sum_flux > 0.0) {
            for (int gt = 0; gt < G; gt++)
                mg->scatter[gf * G + gt] /= sum_flux;
        }
    }
}

/**
 * Add inelastic and (n,Xn) scattering contributions to the transfer matrix.
 *
 * For each non-elastic reaction that produces neutrons:
 * - Inelastic levels (MT 51-90): discrete Q-value determines E_out
 * - (n,2n) MT 16, (n,3n) MT 17, etc.: use Q-value for energy transfer
 *
 * For level scattering (MT 51-90):
 *   E_out = (A/(A+1))² · E_in + Q · A/(A+1)
 *
 * The contribution is weighted by the reaction cross section and
 * neutron yield.
 */
static alea_error_t build_inelastic_scatter(alea_nuc_multigroup_t* mg,
                                             const alea_nuc_nuclide_t* nuc) {
    int G = mg->n_groups;
    double A = nuc->awr;
    const double* egrid = nuc->energy;
    int n = nuc->n_energies;
    double* transfer = calloc((size_t)G, sizeof(double));
    if (!transfer) return ALEA_ERR_OUT_OF_MEMORY;

    for (int r = 0; r < nuc->n_reactions; r++) {
        const alea_nuc_reaction_t* rxn = &nuc->reactions[r];
        int mt = rxn->mt;

        /* Skip non-scattering reactions */
        alea_nuc_reaction_class_t cls = alea_nuc_reaction_classify(mt);
        if (cls == ALEA_NUC_RXN_ABSORPTION) continue;

        /* Skip fission — handled separately */
        if (is_fission_mt(mt)) continue;

        if (!rxn->xs || rxn->n_energies <= 0) continue;

        int ie_start = rxn->threshold_index - 1;
        double Q = rxn->q_value;
        double yield = 1.0;

        /* Determine neutron yield */
        if (mt == 16) yield = 2.0;       /* (n,2n) */
        else if (mt == 17) yield = 3.0;  /* (n,3n) */
        else if (mt == 24) yield = 2.0;  /* (n,2nα) */
        else if (mt == 25) yield = 3.0;  /* (n,3nα) */
        else if (mt == 37) yield = 4.0;  /* (n,4n) */
        else if (mt == 30) yield = 2.0;  /* (n,2n2α) */

        for (int gf = 0; gf < G; gf++) {
            double Eg_hi = mg->bounds[gf];
            double Eg_lo = mg->bounds[gf + 1];

            int i_start = ie_start;
            if (i_start < 0) i_start = 0;
            while (i_start < n - 1 && egrid[i_start + 1] < Eg_lo) i_start++;

            double sum_flux = 0.0;
            memset(transfer, 0, (size_t)G * sizeof(double));

            for (int i = i_start; i < n - 1; i++) {
                if (egrid[i] >= Eg_hi) break;

                double e0 = egrid[i];
                double e1 = egrid[i + 1];
                if (e0 < Eg_lo) e0 = Eg_lo;
                if (e1 > Eg_hi) e1 = Eg_hi;
                if (e1 <= e0) continue;

                double E_mid = 0.5 * (e0 + e1);
                double de = e1 - e0;
                double phi = weight(mg, E_mid);
                sum_flux += phi * de;

                /* Get reaction XS at E_mid */
                double dE_grid = egrid[i + 1] - egrid[i];
                double f_mid = (dE_grid > 0.0) ? (E_mid - egrid[i]) / dE_grid : 0.0;
                int ri = i - ie_start;
                double sig = 0.0;
                if (ri >= 0 && ri < rxn->n_energies - 1)
                    sig = rxn->xs[ri] + f_mid * (rxn->xs[ri + 1] - rxn->xs[ri]);
                else if (ri == rxn->n_energies - 1 && ri >= 0)
                    sig = rxn->xs[ri];
                if (sig <= 0.0) continue;

                double w = sig * phi * de * yield;

                /* Determine outgoing energy */
                double E_out;

                /* Level scattering: discrete Q-value */
                if (mt >= 51 && mt <= 90) {
                    /* E_out = ((A/(A+1))^2) * E_in + Q * A/(A+1) */
                    double mass_frac = A / (A + 1.0);
                    E_out = mass_frac * mass_frac * E_mid + Q * mass_frac;
                    if (E_out <= 0.0) continue;

                    int gt = find_group(mg, E_out);
                    if (gt >= 0)
                        transfer[gt] += w;
                } else if (mt == 91) {
                    /* Continuum inelastic — spread over lower groups
                     * Use a simple model: E_out ∈ [0, E_in + Q] uniformly */
                    double E_max = E_mid + Q;
                    if (E_max <= 0.0) continue;
                    double E_min = 0.0;

                    for (int gt = 0; gt < G; gt++) {
                        double Et_hi = mg->bounds[gt];
                        double Et_lo = mg->bounds[gt + 1];
                        double ov_lo = (E_min > Et_lo) ? E_min : Et_lo;
                        double ov_hi = (E_max < Et_hi) ? E_max : Et_hi;
                        if (ov_hi > ov_lo) {
                            double frac = (ov_hi - ov_lo) / (E_max - E_min);
                            transfer[gt] += w * frac;
                        }
                    }
                } else {
                    /* (n,Xn) reactions: approximate E_out per emitted neutron.
                     * Available energy shared equally among emitted neutrons. */
                    E_out = (E_mid + Q) / yield;
                    if (E_out <= 0.0) continue;

                    int gt = find_group(mg, E_out);
                    if (gt >= 0)
                        transfer[gt] += w;
                }
            }

            /* Normalize and add to scatter matrix */
            if (sum_flux > 0.0) {
                for (int gt = 0; gt < G; gt++)
                    mg->scatter[gf * G + gt] += transfer[gt] / sum_flux;
            }
        }
    }

    free(transfer);
    return ALEA_OK;
}

alea_error_t alea_nuc_mg_collapse(alea_nuc_multigroup_t* mg, const alea_nuc_nuclide_t* nuc) {
    if (!mg || !nuc) return ALEA_ERR_NULL_ARG;
    if (nuc->n_energies < 2 || !nuc->energy ||
        !nuc->sigma_total || !nuc->sigma_abs) {
        return ALEA_ERR_INVALID_ARG;
    }
    for (int i = 0; i < nuc->n_energies; i++) {
        if (!isfinite(nuc->energy[i])) return ALEA_ERR_INVALID_ARG;
        if (i > 0 && !(nuc->energy[i] > nuc->energy[i - 1]))
            return ALEA_ERR_INVALID_ARG;
    }

    int G = mg->n_groups;

    /* Build full-grid fission XS once (outside group loop) */
    double* fission_xs = NULL;
    for (int r = 0; r < nuc->n_reactions; r++) {
        if (is_fission_mt(nuc->reactions[r].mt)) {
            const alea_nuc_reaction_t* rxn = &nuc->reactions[r];
            if (!rxn->xs || rxn->n_energies <= 0) continue;
            if (!fission_xs) {
                fission_xs = calloc((size_t)nuc->n_energies, sizeof(double));
                if (!fission_xs) return ALEA_ERR_OUT_OF_MEMORY;
            }
            int ie_start = rxn->threshold_index - 1;
            for (int i = 0; i < rxn->n_energies; i++)
                if (ie_start + i >= 0 && ie_start + i < nuc->n_energies)
                    fission_xs[ie_start + i] += rxn->xs[i];
        }
    }

    for (int g = 0; g < G; g++) {
        double E_hi = mg->bounds[g];
        double E_lo = mg->bounds[g + 1];

        mg->sigma_t[g] = group_average_xs(mg, nuc->energy, nuc->sigma_total,
                                           nuc->n_energies, E_lo, E_hi);
        mg->sigma_a[g] = group_average_xs(mg, nuc->energy, nuc->sigma_abs,
                                           nuc->n_energies, E_lo, E_hi);

        mg->sigma_f[g] = 0.0;
        mg->nu_sigma_f[g] = 0.0;
        if (fission_xs) {
            mg->sigma_f[g] = group_average_xs(mg, nuc->energy, fission_xs,
                                               nuc->n_energies, E_lo, E_hi);
            mg->nu_sigma_f[g] = group_average_nu_sigma_f(mg, nuc, fission_xs,
                                                          E_lo, E_hi);
        }
    }
    free(fission_xs);

    /* Fission spectrum χ: fraction of fission neutrons born in each group.
     * Use Watt spectrum: χ(E) ∝ exp(-E/a) * sinh(√(bE)).
     * Extract a,b from nuclide's fission energy distribution (Law 11)
     * if available; otherwise fall back to U-235 thermal defaults.
     *
     * For energy-dependent chi, we average Watt parameters over the
     * fission source spectrum (groups weighted by nu-sigma-f). */
    {
        double a = WATT_A_DEFAULT, b = WATT_B_DEFAULT;

        /* Try to extract Watt parameters from nuclide's fission reaction */
        for (int r = 0; r < nuc->n_reactions; r++) {
            if (is_fission_mt(nuc->reactions[r].mt)) {
                const alea_nuc_energy_dist_t* ed = nuc->reactions[r].energy;
                while (ed) {
                    if (ed->law == ALEA_NUC_ELAW_WATT && ed->n_temp > 0 &&
                        ed->temp_T && ed->temp_C) {
                        /* Average Watt parameters over fission groups */
                        double sum_w = 0.0, sum_a = 0.0, sum_b = 0.0;
                        for (int g = 0; g < G; g++) {
                            double nsf = mg->nu_sigma_f[g];
                            if (nsf <= 0.0) continue;
                            double E_mid = sqrt(mg->bounds[g] * mg->bounds[g + 1]);
                            double f;
                            int ie = alea_nuc_energy_lookup(ed->temp_energy,
                                                       ed->n_temp, E_mid, &f);
                            if (ie >= 0) {
                                double a_g = ed->temp_T[ie] +
                                    f * (ed->temp_T[ie + 1] - ed->temp_T[ie]);
                                double b_g = ed->temp_C[ie] +
                                    f * (ed->temp_C[ie + 1] - ed->temp_C[ie]);
                                sum_a += a_g * nsf;
                                sum_b += b_g * nsf;
                                sum_w += nsf;
                            }
                        }
                        if (sum_w > 0.0) {
                            a = sum_a / sum_w;
                            b = sum_b / sum_w;
                        }
                        goto watt_found;
                    }
                    ed = ed->next;
                }
                break;
            }
        }
        watt_found:

        double chi_sum = 0.0;
        for (int g = 0; g < G; g++) {
            double E_hi = mg->bounds[g];
            double E_lo = mg->bounds[g + 1];
            /* Integrate Watt spectrum over group with midpoint quadrature */
            int npts = 100;
            double sum = 0.0;
            for (int k = 0; k < npts; k++) {
                double E = E_lo + (E_hi - E_lo) * (k + 0.5) / npts;
                if (E > 0.0) {
                    double w = exp(-E / a) * sinh(sqrt(b * E));
                    sum += w;
                }
            }
            mg->chi[g] = sum * (E_hi - E_lo) / npts;
            chi_sum += mg->chi[g];
        }
        if (chi_sum > 0.0) {
            for (int g = 0; g < G; g++)
                mg->chi[g] /= chi_sum;
        }
    }

    /* Scattering matrix: elastic + inelastic + (n,Xn) */
    memset(mg->scatter, 0, (size_t)G * (size_t)G * sizeof(double));
    build_elastic_scatter(mg, nuc);
    alea_error_t scatter_err = build_inelastic_scatter(mg, nuc);
    if (scatter_err != ALEA_OK) return scatter_err;

    /* sigma_s = sum of scatter matrix row (total scattering XS per group) */
    for (int g = 0; g < G; g++) {
        double sum = 0.0;
        for (int gp = 0; gp < G; gp++)
            sum += mg->scatter[g * G + gp];
        mg->sigma_s[g] = sum;
    }

    return ALEA_OK;
}

double alea_nuc_mg_scatter(const alea_nuc_multigroup_t* mg, int g_from, int g_to) {
    if (!mg || g_from < 0 || g_from >= mg->n_groups ||
        g_to < 0 || g_to >= mg->n_groups)
        return 0.0;
    return mg->scatter[g_from * mg->n_groups + g_to];
}

double alea_nuc_mg_scatter_adjoint(const alea_nuc_multigroup_t* mg, int g_from, int g_to) {
    /* Adjoint scattering matrix = transpose of forward */
    return alea_nuc_mg_scatter(mg, g_to, g_from);
}

int alea_nuc_mg_sample_scatter(const alea_nuc_multigroup_t* mg, int g_from,
                           double xi, int adjoint) {
    if (!mg || g_from < 0 || g_from >= mg->n_groups) return -1;

    int G = mg->n_groups;

    /* Build CDF from the row (forward) or column (adjoint) */
    double cum = 0.0;
    for (int g = 0; g < G; g++) {
        double val = adjoint ? mg->scatter[g * G + g_from]
                             : mg->scatter[g_from * G + g];
        cum += val;
    }

    if (cum <= 0.0) return g_from; /* no scattering data */

    /* Sample from CDF */
    double target = xi * cum;
    double running = 0.0;
    for (int g = 0; g < G; g++) {
        double val = adjoint ? mg->scatter[g * G + g_from]
                             : mg->scatter[g_from * G + g];
        running += val;
        if (running >= target)
            return g;
    }

    return G - 1;
}
