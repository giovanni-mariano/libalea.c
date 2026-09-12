// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

#include "nuclear_internal.h"
#include <stdio.h>

static int photon_shell_index(const alea_nuc_photon_data_t* ph,
                              int designator) {
    for (int i = 0; i < ph->n_subshells; i++)
        if (ph->subshells[i].designator == designator) return i;
    return -1;
}

static bool relaxation_bound(alea_nuc_photon_data_t* ph, int index,
                             unsigned char* state, size_t* bound) {
    if (state[index] == 2) {
        *bound = ph->subshells[index].max_relaxation_photons;
        return true;
    }
    if (state[index] == 1) return false;
    state[index] = 1;
    alea_nuc_atomic_subshell_t* shell = &ph->subshells[index];
    size_t maximum = 0;
    for (int i = 0; i < shell->n_transitions; i++) {
        const alea_nuc_atomic_transition_t* transition =
            &shell->transitions[i];
        int primary = photon_shell_index(ph, transition->primary_designator);
        int secondary = transition->secondary_designator == 0 ? -1 :
            photon_shell_index(ph, transition->secondary_designator);
        size_t first = 0, second = 0;
        if (primary < 0 ||
            !relaxation_bound(ph, primary, state, &first) ||
            (secondary >= 0 &&
             !relaxation_bound(ph, secondary, state, &second))) return false;
        size_t radiative = transition->secondary_designator == 0 ? 1u : 0u;
        if (first > SIZE_MAX - radiative ||
            second > SIZE_MAX - radiative - first) return false;
        size_t candidate = radiative + first + second;
        if (candidate > maximum) maximum = candidate;
    }
    shell->max_relaxation_photons = maximum;
    state[index] = 2;
    *bound = maximum;
    return true;
}

static alea_error_t decode_compton_profiles(alea_nuc_photon_data_t* ph,
                                            const alea_nuc_ace_table_t* t) {
    int nd = t->nxs[4];
    if (nd == 0) return ALEA_OK;
    if (nd < 0 || nd > 128) return ALEA_ERR_PARSE_ERROR;
    int lneps = t->jxs[5], lswd = t->jxs[8], swd = t->jxs[9];
    if (lneps <= 0 || lswd <= lneps || swd <= 0 ||
        (lswd - lneps) % 3 != 0) return ALEA_ERR_PARSE_ERROR;
    int ns = (lswd - lneps) / 3;
    if (ns < nd || ns > 128 || !xss_range_valid(t, lneps, 3 * ns) ||
        !xss_range_valid(t, lswd, nd)) return ALEA_ERR_PARSE_ERROR;

    ph->compton_shells = alea_nuc_calloc(
        (size_t)ns, sizeof(*ph->compton_shells));
    ph->compton_profiles = alea_nuc_calloc(
        (size_t)nd, sizeof(*ph->compton_profiles));
    if (!ph->compton_shells || !ph->compton_profiles)
        return ALEA_ERR_OUT_OF_MEMORY;
    ph->n_compton_shells = ns;
    ph->n_compton_profiles = nd;

    double probability_sum = 0.0;
    int cumulative = t->nxs[5] == 3;
    for (int i = 0; i < ns; i++) {
        alea_nuc_compton_shell_t* shell = &ph->compton_shells[i];
        shell->electron_count = xss(t, lneps + i);
        shell->binding_energy = xss(t, lneps + ns + i);
        double probability = xss(t, lneps + 2 * ns + i);
        if (!isfinite(shell->electron_count) || shell->electron_count <= 0.0 ||
            !isfinite(shell->binding_energy) || shell->binding_energy < 0.0 ||
            !isfinite(probability) || probability < 0.0)
            return ALEA_ERR_PARSE_ERROR;
        if (cumulative) {
            if (probability < probability_sum || probability > 1.0)
                return ALEA_ERR_PARSE_ERROR;
            probability_sum = probability;
        } else {
            probability_sum += probability;
            if (!isfinite(probability_sum)) return ALEA_ERR_PARSE_ERROR;
        }
        shell->cumulative_probability = probability_sum;
        /* Some later EPR tables retain more relativistic selection bins than
         * profile records. Their trailing bins use the final profile. */
        shell->profile_index = i < nd ? i : nd - 1;
    }
    if (!(probability_sum > 0.0) || fabs(probability_sum - 1.0) > 1e-8)
        return ALEA_ERR_PARSE_ERROR;
    for (int i = 0; i < ns; i++)
        ph->compton_shells[i].cumulative_probability /= probability_sum;

    for (int i = 0; i < nd; i++) {
        int offset = xss_int(t, lswd + i);
        int pos = xss_relative_loc(t, swd, offset);
        if (pos == 0 || !xss_range_valid(t, pos, 2))
            return ALEA_ERR_PARSE_ERROR;
        alea_nuc_compton_profile_t* profile = &ph->compton_profiles[i];
        profile->interpolation = xss_int(t, pos);
        profile->n_momenta = xss_int(t, pos + 1);
        int n = profile->n_momenta;
        if (profile->interpolation != 2 || n < 2 || n > 100000 ||
            n > t->xss_length / 3 || !xss_range_valid(t, pos + 2, 3 * n))
            return ALEA_ERR_PARSE_ERROR;
        profile->momentum = xss_copy(t, pos + 2, n);
        profile->pdf = xss_copy(t, pos + 2 + n, n);
        profile->cdf = xss_copy(t, pos + 2 + 2 * n, n);
        if (!profile->momentum || !profile->pdf || !profile->cdf)
            return ALEA_ERR_OUT_OF_MEMORY;
        for (int j = 0; j < n; j++) {
            if (!isfinite(profile->momentum[j]) || profile->momentum[j] < 0.0 ||
                !isfinite(profile->pdf[j]) || profile->pdf[j] < 0.0 ||
                !isfinite(profile->cdf[j]) || profile->cdf[j] < 0.0 ||
                profile->cdf[j] > 1.0 ||
                (j > 0 && (profile->momentum[j] <= profile->momentum[j - 1] ||
                           profile->cdf[j] < profile->cdf[j - 1])))
                return ALEA_ERR_PARSE_ERROR;
        }
        if (profile->cdf[0] > 1e-10 || profile->cdf[n - 1] < 1.0 - 1e-10)
            return ALEA_ERR_PARSE_ERROR;
    }
    return t->decode_error ? ALEA_ERR_PARSE_ERROR : ALEA_OK;
}

static alea_error_t decode_epr_shells(alea_nuc_photon_data_t* ph,
                                      const alea_nuc_ace_table_t* t) {
    int format = t->nxs[5];
    if (format != 1 && format != 3) return ALEA_OK;
    alea_error_t compton_error = decode_compton_profiles(ph, t);
    if (compton_error != ALEA_OK) return compton_error;
    int ns = t->nxs[6];
    if (ns <= 0 || ns > 128) return ALEA_ERR_PARSE_ERROR;
    int subsh = t->jxs[10], occup = t->jxs[11], bind = t->jxs[12];
    int cprob = t->jxs[13], ntr = t->jxs[14], sphel = t->jxs[15];
    int relo = t->jxs[16], xprob = t->jxs[17];
    if (!xss_range_valid(t, subsh, ns) || !xss_range_valid(t, occup, ns) ||
        !xss_range_valid(t, bind, ns) || !xss_range_valid(t, cprob, ns) ||
        !xss_range_valid(t, ntr, ns) || !xss_range_valid(t, relo, ns) ||
        ph->n_energies > t->xss_length / ns ||
        !xss_range_valid(t, sphel, ns * ph->n_energies) || xprob <= 0)
        return ALEA_ERR_PARSE_ERROR;

    ph->subshells = alea_nuc_calloc((size_t)ns, sizeof(*ph->subshells));
    if (!ph->subshells) return ALEA_ERR_OUT_OF_MEMORY;
    ph->epr_format = format;
    ph->n_subshells = ns;
    double previous_vacancy_cdf = 0.0;
    for (int i = 0; i < ns; i++) {
        alea_nuc_atomic_subshell_t* shell = &ph->subshells[i];
        shell->designator = xss_int(t, subsh + i);
        shell->occupancy = xss(t, occup + i);
        shell->binding_energy = xss(t, bind + i);
        shell->compton_vacancy_probability = xss(t, cprob + i);
        shell->n_transitions = xss_int(t, ntr + i);
        int offset = xss_int(t, relo + i);
        if (shell->designator <= 0 || !isfinite(shell->occupancy) ||
            shell->occupancy < 0.0 || !isfinite(shell->binding_energy) ||
            shell->binding_energy <= 0.0 ||
            !isfinite(shell->compton_vacancy_probability) ||
            shell->compton_vacancy_probability < previous_vacancy_cdf ||
            shell->compton_vacancy_probability > 1.0 ||
            shell->n_transitions < 0 || shell->n_transitions > 100000 ||
            offset < 0 ||
            photon_shell_index(ph, shell->designator) != i)
            return ALEA_ERR_PARSE_ERROR;
        previous_vacancy_cdf = shell->compton_vacancy_probability;
        shell->ln_photoelectric_xs = xss_copy(
            t, sphel + i * ph->n_energies, ph->n_energies);
        if (!shell->ln_photoelectric_xs) return ALEA_ERR_OUT_OF_MEMORY;
        for (int j = 0; j < ph->n_energies; j++) {
            if (!isfinite(shell->ln_photoelectric_xs[j]))
                return ALEA_ERR_PARSE_ERROR;
            if (shell->ln_photoelectric_xs[j] == 0.0)
                shell->ln_photoelectric_xs[j] = -HUGE_VAL;
        }
        if (shell->n_transitions == 0) continue;
        if (offset > t->xss_length - xprob ||
            shell->n_transitions > t->xss_length / 4 ||
            !xss_range_valid(t, xprob + offset, 4 * shell->n_transitions))
            return ALEA_ERR_PARSE_ERROR;
        shell->transitions = alea_nuc_malloc(
            (size_t)shell->n_transitions * sizeof(*shell->transitions));
        if (!shell->transitions) return ALEA_ERR_OUT_OF_MEMORY;
        double previous = 0.0;
        for (int j = 0; j < shell->n_transitions; j++) {
            int pos = xprob + offset + 4 * j;
            alea_nuc_atomic_transition_t* transition = &shell->transitions[j];
            transition->primary_designator = xss_int(t, pos);
            transition->secondary_designator = xss_int(t, pos + 1);
            transition->energy = xss(t, pos + 2);
            transition->cumulative_probability = xss(t, pos + 3);
            if (transition->primary_designator <= 0 ||
                transition->secondary_designator < 0 ||
                !isfinite(transition->energy) ||
                (transition->secondary_designator == 0 &&
                 transition->energy < 0.0) ||
                !isfinite(transition->cumulative_probability) ||
                transition->cumulative_probability < previous ||
                transition->cumulative_probability > 1.0)
                return ALEA_ERR_PARSE_ERROR;
            previous = transition->cumulative_probability;
        }
        if (previous < 1.0 - 1e-10) return ALEA_ERR_PARSE_ERROR;
    }

    for (int i = 0; i < ns; i++) {
        alea_nuc_atomic_subshell_t* shell = &ph->subshells[i];
        for (int j = 0; j < shell->n_transitions; j++) {
            alea_nuc_atomic_transition_t* transition = &shell->transitions[j];
            int first = photon_shell_index(ph, transition->primary_designator);
            int second = transition->secondary_designator == 0 ? -1 :
                photon_shell_index(ph, transition->secondary_designator);
            if (first <= i ||
                (transition->secondary_designator != 0 && second <= i))
                return ALEA_ERR_PARSE_ERROR;
        }
    }
    unsigned char state[128] = {0};
    for (int i = 0; i < ns; i++) {
        size_t ignored;
        if (!relaxation_bound(ph, i, state, &ignored))
            return ALEA_ERR_PARSE_ERROR;
    }
    return t->decode_error ? ALEA_ERR_PARSE_ERROR : ALEA_OK;
}

/**
 * @file xs_decode.c
 * @brief Decode cross-section data from raw ACE XSS array
 *
 * Decodes the ESZ (principal cross sections) and SIG (reaction cross sections)
 * blocks from a loaded ACE table into a alea_nuc_nuclide_t.
 *
 * ACE block layout for continuous neutron (.c) tables:
 *   ESZ block at JXS[1]: 5 arrays of NXS[3] doubles each:
 *     energy grid, σ_total, σ_absorption, σ_elastic, heating
 *
 *   MTR block at JXS[3]: NXS[4] MT numbers
 *   LQR block at JXS[4]: NXS[4] Q-values
 *   TYR block at JXS[5]: NXS[4] reaction types
 *   LSIG block at JXS[6]: NXS[4] locators into SIG block
 *   SIG block at JXS[7]: cross-section sub-arrays
 *
 * All JXS values are 1-based (Fortran convention).
 */

/**
 * Decode ESZ block — principal cross sections on the main energy grid.
 */
static alea_error_t decode_esz(alea_nuc_nuclide_t* nuc, const alea_nuc_ace_table_t* t) {
    int ne = t->nxs[2];    /* NXS[3]: number of energies */
    int esz = t->jxs[0];   /* JXS[1]: start of ESZ block */

    if (ne <= 0 || esz <= 0 || ne > t->xss_length / 5 ||
        !xss_range_valid(t, esz, 5 * ne))
        return ALEA_ERR_INVALID_ARG;

    nuc->n_energies = ne;
    nuc->energy       = xss_copy(t, esz,          ne);
    nuc->sigma_total  = xss_copy(t, esz + ne,     ne);
    nuc->sigma_abs    = xss_copy(t, esz + 2 * ne, ne);
    nuc->sigma_elastic= xss_copy(t, esz + 3 * ne, ne);
    nuc->heating      = xss_copy(t, esz + 4 * ne, ne);

    if (!nuc->energy || !nuc->sigma_total || !nuc->sigma_abs ||
        !nuc->sigma_elastic || !nuc->heating)
        return ALEA_ERR_OUT_OF_MEMORY;

    for (int i = 0; i < ne; i++) {
        if (!isfinite(nuc->energy[i]) ||
            (i > 0 && nuc->energy[i] < nuc->energy[i - 1]))
            return ALEA_ERR_PARSE_ERROR;
        /* ACE stores heating numbers (MeV/collision). Internally/publicly the
         * heating array is a heating cross section (MeV-barn). */
        nuc->heating[i] *= nuc->sigma_total[i];
    }

    return ALEA_OK;
}

/**
 * Decode SIG block — per-reaction cross sections.
 */
static alea_error_t decode_reactions(alea_nuc_nuclide_t* nuc, const alea_nuc_ace_table_t* t) {
    int nr = t->nxs[3];    /* NXS[4]: number of reactions (excluding elastic) */
    if (nr <= 0) {
        nuc->n_reactions = 0;
        nuc->reactions = NULL;
        return ALEA_OK;
    }

    int mtr  = t->jxs[2];  /* JXS[3]: MT numbers */
    int lqr  = t->jxs[3];  /* JXS[4]: Q-values */
    int tyr  = t->jxs[4];  /* JXS[5]: reaction types */
    int lsig = t->jxs[5];  /* JXS[6]: XS locators (relative to JXS[7]) */
    int sig  = t->jxs[6];  /* JXS[7]: SIG block start */

    if (mtr <= 0 || lqr <= 0 || tyr <= 0 || lsig <= 0 || sig <= 0 ||
        !xss_range_valid(t, mtr, nr) || !xss_range_valid(t, lqr, nr) ||
        !xss_range_valid(t, tyr, nr) || !xss_range_valid(t, lsig, nr))
        return ALEA_ERR_INVALID_ARG;

    nuc->reactions = alea_nuc_calloc((size_t)nr, sizeof(alea_nuc_reaction_t));
    if (!nuc->reactions) return ALEA_ERR_OUT_OF_MEMORY;
    nuc->n_reactions = nr;

    for (int i = 0; i < nr; i++) {
        alea_nuc_reaction_t* r = &nuc->reactions[i];

        r->mt      = xss_int(t, mtr + i);
        r->q_value = xss(t, lqr + i);
        r->ty      = xss_int(t, tyr + i);
        r->center_of_mass = r->ty < 0;

        /* Cross-section data location */
        int loc = xss_int(t, lsig + i);  /* relative to SIG block */
        int abs_loc = xss_relative_loc(t, sig, loc);
        if (abs_loc == 0 || !xss_range_valid(t, abs_loc, 2))
            return ALEA_ERR_PARSE_ERROR;

        /* SIG sub-array format: threshold_index, n_energies, xs[n_energies] */
        r->threshold_index = xss_int(t, abs_loc);
        r->n_energies      = xss_int(t, abs_loc + 1);

        if (r->threshold_index < 1 || r->threshold_index > nuc->n_energies ||
            r->n_energies <= 0 ||
            r->n_energies > nuc->n_energies - r->threshold_index + 1)
            return ALEA_ERR_PARSE_ERROR;

        r->xs = xss_copy(t, abs_loc + 2, r->n_energies);
        if (!r->xs) return t->decode_error ? ALEA_ERR_PARSE_ERROR : ALEA_ERR_OUT_OF_MEMORY;

        r->angular = NULL;
        r->energy = NULL;
    }

    return ALEA_OK;
}

/**
 * Decode a single ν̄ representation at the given XSS location.
 * Returns allocated nu_bar, or NULL on failure.
 */
static alea_nuc_nu_bar_t* decode_nu_block(const alea_nuc_ace_table_t* t, int loc) {
    int lnu = xss_int(t, loc);  /* type flag: 1=polynomial, 2=tabular */

    alea_nuc_nu_bar_t* nu = alea_nuc_calloc(1, sizeof(*nu));
    if (!nu) {
        xss_mark_allocation_error(t);
        return NULL;
    }

    if (lnu == 1) {
        nu->type = ALEA_NUC_NU_POLYNOMIAL;
        nu->n_coeffs = xss_int(t, loc + 1);
        if (nu->n_coeffs <= 0 || nu->n_coeffs > 100000) {
            xss_mark_corrupt(t);
            free(nu);
            return NULL;
        }
        nu->coeffs = xss_copy(t, loc + 2, nu->n_coeffs);
        if (!nu->coeffs) { free(nu); return NULL; }
        for (int i = 0; i < nu->n_coeffs; i++) {
            if (!isfinite(nu->coeffs[i])) {
                xss_mark_corrupt(t);
                free(nu->coeffs);
                free(nu);
                return NULL;
            }
        }
    } else if (lnu == 2) {
        nu->type = ALEA_NUC_NU_TABULAR;
        int nr_interp = xss_int(t, loc + 1);
        if (nr_interp < 0 || nr_interp > 100000 ||
            !xss_range_valid(t, loc + 2, 2 * nr_interp + 1)) {
            free(nu);
            return NULL;
        }
        nu->n_regions = nr_interp;
        if (nr_interp > 0) {
            nu->nbt = alea_nuc_malloc((size_t)nr_interp * sizeof(int));
            nu->interp = alea_nuc_malloc((size_t)nr_interp * sizeof(int));
            if (!nu->nbt || !nu->interp) {
                xss_mark_allocation_error(t);
                free(nu->nbt);
                free(nu->interp);
                free(nu);
                return NULL;
            }
            for (int i = 0; i < nr_interp; i++) {
                nu->nbt[i] = xss_int(t, loc + 2 + i);
                nu->interp[i] = xss_int(t, loc + 2 + nr_interp + i);
            }
        }
        int base = loc + 2 + 2 * nr_interp; /* skip NBT/INT pairs */
        nu->n_energies = xss_int(t, base);
        if (nu->n_energies <= 0 || nu->n_energies > 100000 ||
            !xss_range_valid(t, base + 1, 2 * nu->n_energies)) {
            free(nu->nbt); free(nu->interp); free(nu);
            return NULL;
        }
        nu->energy = xss_copy(t, base + 1, nu->n_energies);
        nu->nu     = xss_copy(t, base + 1 + nu->n_energies, nu->n_energies);
        if (!nu->energy || !nu->nu) {
            free(nu->nbt); free(nu->interp);
            free(nu->energy); free(nu->nu); free(nu);
            return NULL;
        }
        for (int i = 0; i < nr_interp; i++) {
            if (nu->nbt[i] < 2 || nu->nbt[i] > nu->n_energies ||
                (i > 0 && nu->nbt[i] <= nu->nbt[i - 1]) ||
                nu->interp[i] < 1 || nu->interp[i] > 5) {
                xss_mark_corrupt(t);
                free(nu->nbt); free(nu->interp);
                free(nu->energy); free(nu->nu); free(nu);
                return NULL;
            }
        }
        if (nr_interp > 0 && nu->nbt[nr_interp - 1] != nu->n_energies) {
            xss_mark_corrupt(t);
            free(nu->nbt); free(nu->interp);
            free(nu->energy); free(nu->nu); free(nu);
            return NULL;
        }
        for (int i = 0; i < nu->n_energies; i++) {
            if (!isfinite(nu->energy[i]) || !isfinite(nu->nu[i]) ||
                nu->nu[i] < 0.0 ||
                (i > 0 && nu->energy[i] <= nu->energy[i - 1])) {
                xss_mark_corrupt(t);
                free(nu->nbt); free(nu->interp);
                free(nu->energy); free(nu->nu); free(nu);
                return NULL;
            }
        }
    } else {
        free(nu);
        return NULL;
    }

    return nu;
}

/**
 * Decode NU block — ν̄ (average neutrons per fission).
 */
static alea_error_t decode_nu(alea_nuc_nuclide_t* nuc, const alea_nuc_ace_table_t* t) {
    int nu_loc = t->jxs[1];  /* JXS[2]: NU block */
    if (nu_loc <= 0) return ALEA_OK;  /* non-fissile nuclide */

    nuc->fission = alea_nuc_calloc(1, sizeof(alea_nuc_fission_t));
    if (!nuc->fission) return ALEA_ERR_OUT_OF_MEMORY;

    /*
     * NU block format:
     *   If first value > 0: single representation (prompt or total)
     *   If first value < 0: both prompt and total present
     *     At nu_loc: KNU (negative flag)
     *     At nu_loc+1: prompt ν̄ block (|KNU| words)
     *     After prompt: total ν̄ block
     */
    int knu = xss_int(t, nu_loc);

    if (knu == INT_MIN) return ALEA_ERR_PARSE_ERROR;
    if (knu < 0) {
        /* Both prompt and total ν̄ present */
        if (nu_loc >= t->xss_length) return ALEA_ERR_PARSE_ERROR;
        int prompt_loc = nu_loc + 1;
        alea_nuc_nu_bar_t* prompt = decode_nu_block(t, prompt_loc);
        if (!prompt)
            return t->allocation_error ? ALEA_ERR_OUT_OF_MEMORY
                                       : ALEA_ERR_PARSE_ERROR;
        nuc->fission->prompt = prompt;

        if (abs(knu) > t->xss_length - prompt_loc)
            return ALEA_ERR_PARSE_ERROR;
        int total_loc = prompt_loc + abs(knu);
        alea_nuc_nu_bar_t* total = decode_nu_block(t, total_loc);
        if (!total)
            return t->allocation_error ? ALEA_ERR_OUT_OF_MEMORY
                                       : ALEA_ERR_PARSE_ERROR;
        nuc->fission->total = total;
    } else {
        /* Single ν̄ representation (total) */
        alea_nuc_nu_bar_t* total = decode_nu_block(t, nu_loc);
        if (!total)
            return t->allocation_error ? ALEA_ERR_OUT_OF_MEMORY
                                       : ALEA_ERR_PARSE_ERROR;
        nuc->fission->total = total;
    }

    return ALEA_OK;
}

alea_error_t alea_nuc_decode_delayed_neutrons(
    alea_nuc_nuclide_t* nuc, const alea_nuc_ace_table_t* t) {
    int groups = t->nxs[7];       /* NXS(8): precursor families */
    int dnu = t->jxs[23];         /* JXS(24): delayed nubar */
    int bdd = t->jxs[24];         /* JXS(25): precursor data */
    int dnedl = t->jxs[25];       /* JXS(26): spectrum locators */
    int dned = t->jxs[26];        /* JXS(27): spectra */
    if (groups == 0 && dnu == 0 && bdd == 0 && dnedl == 0 && dned == 0)
        return ALEA_OK;
    if (groups <= 0 || groups > 64 || dnu <= 0 || bdd <= 0 ||
        dnedl <= 0 || dned <= 0 || !xss_range_valid(t, dnedl, groups))
        return ALEA_ERR_PARSE_ERROR;
    if (!nuc->fission) {
        nuc->fission = alea_nuc_calloc(1, sizeof(*nuc->fission));
        if (!nuc->fission) return ALEA_ERR_OUT_OF_MEMORY;
    }
    nuc->fission->delayed = decode_nu_block(t, dnu);
    if (!nuc->fission->delayed)
        return t->allocation_error ? ALEA_ERR_OUT_OF_MEMORY
                                   : ALEA_ERR_PARSE_ERROR;
    nuc->fission->delayed_groups = alea_nuc_calloc(
        (size_t)groups, sizeof(*nuc->fission->delayed_groups));
    if (!nuc->fission->delayed_groups) return ALEA_ERR_OUT_OF_MEMORY;
    nuc->fission->n_delayed_groups = groups;

    int pos = bdd;
    for (int g = 0; g < groups; g++) {
        alea_nuc_delayed_group_t* group = &nuc->fission->delayed_groups[g];
        if (!xss_range_valid(t, pos, 3)) return ALEA_ERR_PARSE_ERROR;
        double inverse_shakes = xss(t, pos);
        int nr = xss_int(t, pos + 1);
        if (!isfinite(inverse_shakes) || inverse_shakes <= 0.0 || nr < 0 ||
            nr > (t->xss_length - pos - 2) / 2)
            return ALEA_ERR_PARSE_ERROR;
        int ne_pos = pos + 2 + 2 * nr;
        int ne = xss_int(t, ne_pos);
        if (ne <= 0 || ne > 100000 ||
            !xss_range_valid(t, ne_pos + 1, 2 * ne))
            return ALEA_ERR_PARSE_ERROR;
        group->decay_rate = inverse_shakes * 1.0e8;
        group->n_regions = nr;
        group->n_energies = ne;
        if (nr > 0) {
            group->nbt = alea_nuc_malloc((size_t)nr * sizeof(int));
            group->interp = alea_nuc_malloc((size_t)nr * sizeof(int));
            if (!group->nbt || !group->interp) return ALEA_ERR_OUT_OF_MEMORY;
            for (int i = 0; i < nr; i++) {
                group->nbt[i] = xss_int(t, pos + 2 + i);
                group->interp[i] = xss_int(t, pos + 2 + nr + i);
            }
        }
        group->energy = xss_copy(t, ne_pos + 1, ne);
        group->probability = xss_copy(t, ne_pos + 1 + ne, ne);
        if (!group->energy || !group->probability)
            return ALEA_ERR_OUT_OF_MEMORY;
        if (!alea_nuc_interp_regions_valid(group->nbt, group->interp,
                                           nr, ne))
            return ALEA_ERR_PARSE_ERROR;
        for (int i = 0; i < ne; i++) {
            if (!isfinite(group->energy[i]) ||
                !isfinite(group->probability[i]) ||
                group->probability[i] < 0.0 ||
                (i > 0 && group->energy[i] <= group->energy[i - 1]))
                return ALEA_ERR_PARSE_ERROR;
        }
        int spectrum_locator = xss_int(t, dnedl + g);
        group->spectrum = alea_nuc_decode_energy_dist_base(
            t, spectrum_locator, dned);
        if (!group->spectrum)
            return t->allocation_error ? ALEA_ERR_OUT_OF_MEMORY
                                       : ALEA_ERR_PARSE_ERROR;
        pos += 3 + 2 * nr + 2 * ne;
    }
    return ALEA_OK;
}

/**
 * Decode photoatomic (.p) table cross sections.
 *
 * mcplib ESZG block layout at JXS[1], 5 arrays of NXS[3] doubles each,
 * all stored as natural logarithms:
 *   ln(E), ln(σ_incoherent), ln(σ_coherent), ln(σ_PE), ln(σ_pair)
 *
 * Below pair production threshold, ln(σ_pair) = 0.0 is a sentinel
 * meaning zero cross section (exp(0) = 1 would be wrong).
 *
 * Additional blocks:
 *   JXS[2] = JINC: 21-point incoherent scattering function S(q,Z)
 *   JXS[3] = JCOH: 55-point coherent form factor (momentum + integrated FF)
 *   JXS[4] = JFLO: fluorescence edge data
 *   JXS[5] = LHNM: heating numbers (NE values)
 */
static alea_error_t decode_photon(alea_nuc_nuclide_t* nuc, const alea_nuc_ace_table_t* t) {
    int ne = t->nxs[2];    /* NXS[3]: number of energies */
    if (ne <= 0 || ne > t->xss_length / 5) return ALEA_ERR_INVALID_ARG;

    nuc->photon = alea_nuc_calloc(1, sizeof(alea_nuc_photon_data_t));
    if (!nuc->photon) return ALEA_ERR_OUT_OF_MEMORY;

    alea_nuc_photon_data_t* ph = nuc->photon;
    ph->n_energies = ne;

    /* ESZG block: 5 arrays of NE starting at JXS[1], stored as ln values */
    int esz = t->jxs[0]; /* JXS[1] */
    if (!xss_range_valid(t, esz, 5 * ne)) return ALEA_ERR_INVALID_ARG;

    ph->energy              = xss_copy(t, esz,          ne);
    ph->sigma_incoherent    = xss_copy(t, esz + ne,     ne);
    ph->sigma_coherent      = xss_copy(t, esz + 2 * ne, ne);
    ph->sigma_photoelectric = xss_copy(t, esz + 3 * ne, ne);
    ph->sigma_pair          = xss_copy(t, esz + 4 * ne, ne);

    if (!ph->energy || !ph->sigma_incoherent || !ph->sigma_coherent ||
        !ph->sigma_photoelectric || !ph->sigma_pair)
        return ALEA_ERR_OUT_OF_MEMORY;

    for (int i = 0; i < ne; i++) {
        if (!isfinite(ph->energy[i]) ||
            (i > 0 && ph->energy[i] < ph->energy[i - 1]))
            return ALEA_ERR_PARSE_ERROR;
    }

    /* Save log-space copies before converting to linear (for fast log-log interp) */
    ph->ln_energy              = alea_nuc_malloc((size_t)ne * sizeof(double));
    ph->ln_sigma_incoherent    = alea_nuc_malloc((size_t)ne * sizeof(double));
    ph->ln_sigma_coherent      = alea_nuc_malloc((size_t)ne * sizeof(double));
    ph->ln_sigma_photoelectric = alea_nuc_malloc((size_t)ne * sizeof(double));
    ph->ln_sigma_pair          = alea_nuc_malloc((size_t)ne * sizeof(double));
    if (!ph->ln_energy || !ph->ln_sigma_incoherent ||
        !ph->ln_sigma_coherent || !ph->ln_sigma_photoelectric ||
        !ph->ln_sigma_pair) return ALEA_ERR_OUT_OF_MEMORY;
    memcpy(ph->ln_energy, ph->energy, (size_t)ne * sizeof(double));
    memcpy(ph->ln_sigma_incoherent, ph->sigma_incoherent,
           (size_t)ne * sizeof(double));
    memcpy(ph->ln_sigma_coherent, ph->sigma_coherent,
           (size_t)ne * sizeof(double));
    memcpy(ph->ln_sigma_photoelectric, ph->sigma_photoelectric,
           (size_t)ne * sizeof(double));
    for (int i = 0; i < ne; i++)
        ph->ln_sigma_pair[i] = (ph->sigma_pair[i] == 0.0)
                            ? -HUGE_VAL : ph->sigma_pair[i];

    /* Convert from natural log to linear */
    for (int i = 0; i < ne; i++) {
        ph->energy[i] = exp(ph->energy[i]);
        ph->sigma_incoherent[i] = exp(ph->sigma_incoherent[i]);
        ph->sigma_coherent[i] = exp(ph->sigma_coherent[i]);
        ph->sigma_photoelectric[i] = exp(ph->sigma_photoelectric[i]);

        /* Pair production: ln(σ) = 0.0 is sentinel for "zero below threshold" */
        if (ph->sigma_pair[i] != 0.0)
            ph->sigma_pair[i] = exp(ph->sigma_pair[i]);
        if (!isfinite(ph->energy[i])) return ALEA_ERR_PARSE_ERROR;
    }

    /* Heating numbers at JXS[5] if present */
    int lhnm = t->jxs[4]; /* JXS[5]: heating numbers */
    if (lhnm > 0 && xss_range_valid(t, lhnm, ne)) {
        ph->heating = xss_copy(t, lhnm, ne);
    } else {
        ph->heating = alea_nuc_calloc((size_t)ne, sizeof(double));
    }
    if (!ph->heating) return ALEA_ERR_OUT_OF_MEMORY;

    /* Build total XS and nuclide-level arrays */
    nuc->n_energies = ne;
    nuc->energy = alea_nuc_malloc((size_t)ne * sizeof(double));
    nuc->sigma_total = alea_nuc_malloc((size_t)ne * sizeof(double));
    if (!nuc->energy || !nuc->sigma_total) return ALEA_ERR_OUT_OF_MEMORY;

    for (int i = 0; i < ne; i++) {
        nuc->energy[i] = ph->energy[i];
        nuc->sigma_total[i] = ph->sigma_incoherent[i] +
                               ph->sigma_coherent[i] +
                               ph->sigma_photoelectric[i] +
                               ph->sigma_pair[i];
    }

    int epr_format = t->nxs[5];

    /* Incoherent scattering function S(q,Z) at JXS[2]. */
    int jinc = t->jxs[1]; /* JXS[2]: JINC */
    if (jinc > 0) {
        if (epr_format == 1 || epr_format == 3) {
            int n = epr_format == 3 ? t->nxs[12] :
                (t->jxs[2] - jinc) / 2;
            if (n < 2 || n > t->xss_length / 2 ||
                (epr_format == 1 && t->jxs[2] - jinc != 2 * n) ||
                !xss_range_valid(t, jinc, 2 * n)) return ALEA_ERR_PARSE_ERROR;
            ph->n_incoherent_ff = n;
            ph->incoherent_momentum = xss_copy(t, jinc, n);
            ph->incoherent_ff = xss_copy(t, jinc + n, n);
        } else {
            ph->n_incoherent_ff = 21;
            ph->incoherent_ff = xss_copy(t, jinc, 21);
            ph->incoherent_momentum = alea_nuc_malloc(21 * sizeof(double));
            if (ph->incoherent_momentum) {
            static const double jinc_grid[21] = {
                0.0, 0.005, 0.01, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5,
                0.6, 0.7, 0.8, 0.9, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 8.0
            };
            memcpy(ph->incoherent_momentum, jinc_grid, 21 * sizeof(double));
            }
        }
        if (!ph->incoherent_ff || !ph->incoherent_momentum)
            return t->decode_error ? ALEA_ERR_PARSE_ERROR
                                   : ALEA_ERR_OUT_OF_MEMORY;
    }

    /* Coherent form-factor data at JXS[3]. */
    int jcoh = t->jxs[2]; /* JXS[3]: JCOH */
    if (jcoh > 0) {
        if (epr_format == 1 || epr_format == 3) {
            int n = epr_format == 3 ? t->nxs[13] :
                (t->jxs[3] - jcoh) / 3;
            if (n < 2 || n > t->xss_length / 3 ||
                (epr_format == 1 && t->jxs[3] - jcoh != 3 * n) ||
                !xss_range_valid(t, jcoh, 3 * n)) return ALEA_ERR_PARSE_ERROR;
            ph->n_coherent_ff = n;
            ph->coherent_momentum = xss_copy(t, jcoh, n);
            ph->coherent_ff_cumulative = xss_copy(t, jcoh + n, n);
            ph->coherent_ff = xss_copy(t, jcoh + 2 * n, n);
        } else {
            ph->n_coherent_ff = 55;
            ph->coherent_momentum = alea_nuc_malloc(55 * sizeof(double));
            if (ph->coherent_momentum) {
            /* Standard MCNP momentum transfer grid for coherent (inverse Angstroms) */
            static const double jcoh_grid[55] = {
                0.0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.08, 0.10, 0.12,
                0.15, 0.18, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.55,
                0.60, 0.70, 0.80, 0.90, 1.00, 1.10, 1.20, 1.30, 1.40, 1.50,
                1.60, 1.70, 1.80, 1.90, 2.00, 2.20, 2.40, 2.60, 2.80, 3.00,
                3.20, 3.40, 3.60, 3.80, 4.00, 4.20, 4.40, 4.60, 4.80, 5.00,
                5.20, 5.40, 5.60, 5.80, 6.00
            };
            memcpy(ph->coherent_momentum, jcoh_grid, 55 * sizeof(double));
            }
            ph->coherent_ff_cumulative = xss_copy(t, jcoh, 55);
            ph->coherent_ff = xss_copy(t, jcoh + 55, 55);
        }
        if (!ph->coherent_momentum || !ph->coherent_ff ||
            !ph->coherent_ff_cumulative)
            return t->decode_error ? ALEA_ERR_PARSE_ERROR
                                   : ALEA_ERR_OUT_OF_MEMORY;
    }

    /* JFLO stores the Cashwell-Everett averaged fluorescence representation
     * as four arrays of NXS(4) values: edge, phi, cumulative yield, and
     * representative photon energy. Preserve it for inspection and future
     * relaxation sampling; it is not a shell-resolved transition table. */
    alea_error_t epr_error = decode_epr_shells(ph, t);
    if (epr_error != ALEA_OK) return epr_error;

    int nflo = (epr_format == 1 || epr_format == 3) ? 0 : t->nxs[3];
    int jflo = t->jxs[3];
    if (nflo < 0 || nflo > 1000) return ALEA_ERR_PARSE_ERROR;
    if (nflo > 0) {
        if (jflo <= 0 || !xss_range_valid(t, jflo, 4 * nflo))
            return ALEA_ERR_PARSE_ERROR;
        ph->n_fluorescence = nflo;
        ph->fluorescence_edge = xss_copy(t, jflo, nflo);
        ph->fluorescence_phi = xss_copy(t, jflo + nflo, nflo);
        ph->fluorescence_yield = xss_copy(t, jflo + 2 * nflo, nflo);
        ph->fluorescence_energy = xss_copy(t, jflo + 3 * nflo, nflo);
        if (!ph->fluorescence_edge || !ph->fluorescence_phi ||
            !ph->fluorescence_yield || !ph->fluorescence_energy)
            return ALEA_ERR_OUT_OF_MEMORY;
        for (int i = 0; i < nflo; i++) {
            if (!isfinite(ph->fluorescence_edge[i]) ||
                !isfinite(ph->fluorescence_phi[i]) ||
                !isfinite(ph->fluorescence_yield[i]) ||
                !isfinite(ph->fluorescence_energy[i]) ||
                ph->fluorescence_edge[i] < 0.0 ||
                ph->fluorescence_phi[i] < 0.0 ||
                ph->fluorescence_yield[i] < 0.0 ||
                ph->fluorescence_energy[i] < 0.0 ||
                ph->fluorescence_energy[i] > ph->fluorescence_edge[i] ||
                (i > 0 &&
                 (ph->fluorescence_edge[i] < ph->fluorescence_edge[i - 1] ||
                  ph->fluorescence_phi[i] < ph->fluorescence_phi[i - 1] ||
                  ph->fluorescence_yield[i] <
                      ph->fluorescence_yield[i - 1])))
                return ALEA_ERR_PARSE_ERROR;
        }
    } else if (epr_format != 1 && epr_format != 3 &&
               jflo > 0 && jflo != lhnm) {
        return ALEA_ERR_PARSE_ERROR;
    }

    return ALEA_OK;
}

/**
 * Decode URR probability tables from JXS[23].
 *
 * Table layout (per energy): 6 blocks of M values each:
 *   [M cumulative probs] [M σ_total] [M σ_elastic] [M σ_fission] [M σ_capture] [M heating]
 */
static alea_error_t decode_urr(alea_nuc_nuclide_t* nuc, const alea_nuc_ace_table_t* t) {
    int loc = t->jxs[22]; /* JXS[23]: URR block */
    if (loc <= 0) return ALEA_OK; /* no URR data */

    int N   = xss_int(t, loc);      /* number of energies */
    int M   = xss_int(t, loc + 1);  /* number of probability bands */
    int interp = xss_int(t, loc + 2);
    int ilf = xss_int(t, loc + 3);  /* inelastic flag */
    int ioa = xss_int(t, loc + 4);  /* absorption flag */
    int iff = xss_int(t, loc + 5);  /* multiply smooth cross sections */

    if (N <= 0 || M <= 0 || N > 10000 || M > 1000) return ALEA_OK;

    nuc->urr = alea_nuc_calloc(1, sizeof(alea_nuc_urr_t));
    if (!nuc->urr) return ALEA_ERR_OUT_OF_MEMORY;

    nuc->urr->n_energies = N;
    nuc->urr->n_bands = M;
    nuc->urr->interp = interp;
    nuc->urr->inelastic_flag = ilf;
    nuc->urr->absorption_flag = ioa;
    nuc->urr->multiply_smooth = (iff == 1);

    nuc->urr->energy = xss_copy(t, loc + 6, N);
    if (M > INT_MAX / 6 || N > INT_MAX / (6 * M))
        return ALEA_ERR_PARSE_ERROR;
    int table_size = N * 6 * M;
    nuc->urr->table = xss_copy(t, loc + 6 + N, table_size);

    if (!nuc->urr->energy || !nuc->urr->table) return ALEA_ERR_OUT_OF_MEMORY;

    return ALEA_OK;
}

/**
 * Full nuclide decode: load ACE table + decode all blocks.
 * Caller owns the returned nuclide.
 */
alea_nuc_nuclide_t* alea_nuc_load_nuclide(const alea_nuc_xsdir_t* xsdir, const char* zaid) {
    if (!xsdir || !zaid) return NULL;

    /* Find in xsdir */
    const alea_nuc_xsdir_entry_t* entry = alea_nuc_xsdir_find(xsdir, zaid);
    if (!entry) {
        ALEA_LOG_WARN("ZAID '%s' not found in xsdir", zaid);
        return NULL;
    }

    /* Resolve file path */
    char filepath_buf[1024];
    char* filepath = filepath_buf;
    char* filepath_alloc = NULL;
    if (entry->filename[0] == '/') {
        size_t len = strlen(entry->filename);
        if (len >= sizeof(filepath_buf)) {
            filepath_alloc = alea_nuc_malloc(len + 1);
            if (!filepath_alloc) return NULL;
            filepath = filepath_alloc;
        }
        memcpy(filepath, entry->filename, len + 1);
    } else if (xsdir->datapath[0]) {
        size_t len = strlen(xsdir->datapath) + 1 + strlen(entry->filename);
        if (len >= sizeof(filepath_buf)) {
            filepath_alloc = alea_nuc_malloc(len + 1);
            if (!filepath_alloc) return NULL;
            filepath = filepath_alloc;
        }
        snprintf(filepath, len + 1, "%s/%s",
                 xsdir->datapath, entry->filename);
    } else {
        size_t len = strlen(entry->filename);
        if (len >= sizeof(filepath_buf)) {
            filepath_alloc = alea_nuc_malloc(len + 1);
            if (!filepath_alloc) return NULL;
            filepath = filepath_alloc;
        }
        memcpy(filepath, entry->filename, len + 1);
    }

    /* Read raw ACE table */
    alea_nuc_ace_table_t raw;
    alea_error_t err = alea_nuc_ace_read(filepath, entry->address, entry->file_type, &raw);
    if (err != ALEA_OK) {
        ALEA_LOG_ERROR("failed to read ACE file '%s': %s",
                      filepath, alea_error_string(err));
        free(filepath_alloc);
        return NULL;
    }
    free(filepath_alloc);

    /* Allocate nuclide */
    alea_nuc_nuclide_t* nuc = alea_nuc_calloc(1, sizeof(*nuc));
    if (!nuc) { alea_nuc_ace_free(&raw); return NULL; }

    strncpy(nuc->zaid, zaid, sizeof(nuc->zaid) - 1);
    nuc->awr = raw.awr;
    nuc->temperature = raw.temperature;
    nuc->raw = raw; /* transfer ownership of XSS */

    /* Parse Z, A from ZAID */
    alea_nuc_parse_zaid(zaid, &nuc->Z, &nuc->A, &nuc->metastable, NULL);

    /* Decode based on table type */
    if (raw.type == ALEA_NUC_TABLE_CONTINUOUS_NEUTRON) {
        nuc->particle = ALEA_NUC_PARTICLE_NEUTRON;
        err = decode_esz(nuc, &raw);
        if (err == ALEA_OK) err = decode_reactions(nuc, &raw);
        if (err == ALEA_OK) err = decode_nu(nuc, &raw);
        if (err == ALEA_OK) err = alea_nuc_decode_delayed_neutrons(nuc, &raw);
        if (err == ALEA_OK) err = alea_nuc_decode_photon_production(nuc, &raw);
        if (err == ALEA_OK) err = decode_urr(nuc, &raw);
        if (err == ALEA_OK) {
            /* Decode angular and energy distributions */
            alea_nuc_decode_all_angular(nuc);
            alea_nuc_decode_all_energy(nuc);

            /* Build MT → reaction index lookup table */
            nuc->mt_to_rxn = alea_nuc_malloc(ALEA_NUC_MT_TABLE_SIZE * sizeof(int));
            if (nuc->mt_to_rxn) {
                memset(nuc->mt_to_rxn, 0xFF,
                       ALEA_NUC_MT_TABLE_SIZE * sizeof(int)); /* -1 */
                for (int i = 0; i < nuc->n_reactions; i++) {
                    int mt = nuc->reactions[i].mt;
                    if (mt >= 0 && mt < ALEA_NUC_MT_TABLE_SIZE)
                        nuc->mt_to_rxn[mt] = i;
                }
            } else {
                err = ALEA_ERR_OUT_OF_MEMORY;
            }
        }
    } else if (raw.type == ALEA_NUC_TABLE_PHOTOATOMIC) {
        nuc->particle = ALEA_NUC_PARTICLE_PHOTON;
        err = decode_photon(nuc, &raw);
    } else {
        ALEA_LOG_WARN("unsupported table type for '%s'", zaid);
        alea_nuc_nuclide_free(nuc);
        return NULL;
    }

    if (err != ALEA_OK) {
        ALEA_LOG_ERROR("failed to decode '%s': %s", zaid, alea_error_string(err));
        alea_nuc_nuclide_free(nuc);
        return NULL;
    }

    /* Fail closed: if any decode step touched the XSS array out of bounds, the
     * data is corrupt and the decoded values cannot be trusted. Reject the
     * nuclide rather than returning fabricated zero-valued physics. The phase-1
     * decoders operate on `raw`; the angular/energy decoders on `nuc->raw`
     * (a copy sharing the same XSS) — inspect both. */
    if (raw.decode_error || nuc->raw.decode_error) {
        ALEA_LOG_ERROR("rejecting '%s': out-of-bounds XSS access during decode "
                       "(corrupt or malformed ACE data)", zaid);
        alea_nuc_nuclide_free(nuc);
        return NULL;
    }
    if (raw.allocation_error || nuc->raw.allocation_error) {
        ALEA_LOG_ERROR("rejecting '%s': allocation failed during ACE decode", zaid);
        alea_nuc_nuclide_free(nuc);
        return NULL;
    }

    ALEA_LOG_INFO("loaded %s (Z=%d A=%d, %d energies, %d reactions)",
                 zaid, nuc->Z, nuc->A, nuc->n_energies, nuc->n_reactions);

    return nuc;
}
