// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

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

#include "nuclear_internal.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

/**
 * Decode ESZ block — principal cross sections on the main energy grid.
 */
static alea_error_t decode_esz(alea_nuc_nuclide_t* nuc, const alea_nuc_ace_table_t* t) {
    int ne = t->nxs[2];    /* NXS[3]: number of energies */
    int esz = t->jxs[0];   /* JXS[1]: start of ESZ block */

    if (ne <= 0 || esz <= 0) return ALEA_ERR_INVALID_ARG;

    nuc->n_energies = ne;
    nuc->energy       = xss_copy(t, esz,          ne);
    nuc->sigma_total  = xss_copy(t, esz + ne,     ne);
    nuc->sigma_abs    = xss_copy(t, esz + 2 * ne, ne);
    nuc->sigma_elastic= xss_copy(t, esz + 3 * ne, ne);
    nuc->heating      = xss_copy(t, esz + 4 * ne, ne);

    if (!nuc->energy || !nuc->sigma_total || !nuc->sigma_abs ||
        !nuc->sigma_elastic || !nuc->heating)
        return ALEA_ERR_OUT_OF_MEMORY;

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

    if (mtr <= 0 || lqr <= 0 || tyr <= 0 || lsig <= 0 || sig <= 0)
        return ALEA_ERR_INVALID_ARG;

    nuc->reactions = calloc((size_t)nr, sizeof(alea_nuc_reaction_t));
    if (!nuc->reactions) return ALEA_ERR_OUT_OF_MEMORY;
    nuc->n_reactions = nr;

    for (int i = 0; i < nr; i++) {
        alea_nuc_reaction_t* r = &nuc->reactions[i];

        r->mt      = xss_int(t, mtr + i);
        r->q_value = xss(t, lqr + i);
        r->ty      = xss_int(t, tyr + i);

        /* Cross-section data location */
        int loc = xss_int(t, lsig + i);  /* relative to SIG block */
        int abs_loc = sig + loc - 1;     /* 1-based absolute position */

        /* SIG sub-array format: threshold_index, n_energies, xs[n_energies] */
        r->threshold_index = xss_int(t, abs_loc);
        r->n_energies      = xss_int(t, abs_loc + 1);

        if (r->n_energies > 0) {
            r->xs = xss_copy(t, abs_loc + 2, r->n_energies);
            if (!r->xs) return ALEA_ERR_OUT_OF_MEMORY;
        }

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

    alea_nuc_nu_bar_t* nu = calloc(1, sizeof(*nu));
    if (!nu) return NULL;

    if (lnu == 1) {
        nu->type = ALEA_NUC_NU_POLYNOMIAL;
        nu->n_coeffs = xss_int(t, loc + 1);
        nu->coeffs = xss_copy(t, loc + 2, nu->n_coeffs);
        if (!nu->coeffs) { free(nu); return NULL; }
    } else if (lnu == 2) {
        nu->type = ALEA_NUC_NU_TABULAR;
        int nr_interp = xss_int(t, loc + 1);
        int base = loc + 2 + 2 * nr_interp; /* skip NBT/INT pairs */
        nu->n_energies = xss_int(t, base);
        nu->energy = xss_copy(t, base + 1, nu->n_energies);
        nu->nu     = xss_copy(t, base + 1 + nu->n_energies, nu->n_energies);
        if (!nu->energy || !nu->nu) {
            free(nu->energy); free(nu->nu); free(nu);
            return NULL;
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

    nuc->fission = calloc(1, sizeof(alea_nuc_fission_t));
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

    if (knu < 0) {
        /* Both prompt and total ν̄ present */
        int prompt_loc = nu_loc + 1;
        alea_nuc_nu_bar_t* prompt = decode_nu_block(t, prompt_loc);
        if (prompt) nuc->fission->prompt = prompt;

        int total_loc = prompt_loc + abs(knu);
        alea_nuc_nu_bar_t* total = decode_nu_block(t, total_loc);
        if (total) nuc->fission->total = total;
    } else {
        /* Single ν̄ representation (total) */
        alea_nuc_nu_bar_t* total = decode_nu_block(t, nu_loc);
        if (total) nuc->fission->total = total;
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
    if (ne <= 0) return ALEA_ERR_INVALID_ARG;

    nuc->photon = calloc(1, sizeof(alea_nuc_photon_data_t));
    if (!nuc->photon) return ALEA_ERR_OUT_OF_MEMORY;

    alea_nuc_photon_data_t* ph = nuc->photon;
    ph->n_energies = ne;

    /* ESZG block: 5 arrays of NE starting at JXS[1], stored as ln values */
    int esz = t->jxs[0]; /* JXS[1] */

    ph->energy              = xss_copy(t, esz,          ne);
    ph->sigma_incoherent    = xss_copy(t, esz + ne,     ne);
    ph->sigma_coherent      = xss_copy(t, esz + 2 * ne, ne);
    ph->sigma_photoelectric = xss_copy(t, esz + 3 * ne, ne);
    ph->sigma_pair          = xss_copy(t, esz + 4 * ne, ne);

    if (!ph->energy || !ph->sigma_incoherent || !ph->sigma_coherent ||
        !ph->sigma_photoelectric || !ph->sigma_pair)
        return ALEA_ERR_OUT_OF_MEMORY;

    /* Save log-space copies before converting to linear (for fast log-log interp) */
    ph->ln_energy              = malloc((size_t)ne * sizeof(double));
    ph->ln_sigma_incoherent    = malloc((size_t)ne * sizeof(double));
    ph->ln_sigma_coherent      = malloc((size_t)ne * sizeof(double));
    ph->ln_sigma_photoelectric = malloc((size_t)ne * sizeof(double));
    ph->ln_sigma_pair          = malloc((size_t)ne * sizeof(double));
    if (ph->ln_energy) memcpy(ph->ln_energy, ph->energy, (size_t)ne * sizeof(double));
    if (ph->ln_sigma_incoherent) memcpy(ph->ln_sigma_incoherent, ph->sigma_incoherent, (size_t)ne * sizeof(double));
    if (ph->ln_sigma_coherent) memcpy(ph->ln_sigma_coherent, ph->sigma_coherent, (size_t)ne * sizeof(double));
    if (ph->ln_sigma_photoelectric) memcpy(ph->ln_sigma_photoelectric, ph->sigma_photoelectric, (size_t)ne * sizeof(double));
    if (ph->ln_sigma_pair) {
        for (int i = 0; i < ne; i++)
            ph->ln_sigma_pair[i] = (ph->sigma_pair[i] == 0.0) ? -HUGE_VAL : ph->sigma_pair[i];
    }

    /* Convert from natural log to linear */
    for (int i = 0; i < ne; i++) {
        ph->energy[i] = exp(ph->energy[i]);
        ph->sigma_incoherent[i] = exp(ph->sigma_incoherent[i]);
        ph->sigma_coherent[i] = exp(ph->sigma_coherent[i]);
        ph->sigma_photoelectric[i] = exp(ph->sigma_photoelectric[i]);

        /* Pair production: ln(σ) = 0.0 is sentinel for "zero below threshold" */
        if (ph->sigma_pair[i] != 0.0)
            ph->sigma_pair[i] = exp(ph->sigma_pair[i]);
    }

    /* Heating numbers at JXS[5] if present */
    int lhnm = t->jxs[4]; /* JXS[5]: heating numbers */
    if (lhnm > 0 && lhnm + ne - 1 <= t->xss_length) {
        ph->heating = xss_copy(t, lhnm, ne);
    } else {
        ph->heating = calloc((size_t)ne, sizeof(double));
    }
    if (!ph->heating) return ALEA_ERR_OUT_OF_MEMORY;

    /* Build total XS and nuclide-level arrays */
    nuc->n_energies = ne;
    nuc->energy = malloc((size_t)ne * sizeof(double));
    nuc->sigma_total = malloc((size_t)ne * sizeof(double));
    if (!nuc->energy || !nuc->sigma_total) return ALEA_ERR_OUT_OF_MEMORY;

    for (int i = 0; i < ne; i++) {
        nuc->energy[i] = ph->energy[i];
        nuc->sigma_total[i] = ph->sigma_incoherent[i] +
                               ph->sigma_coherent[i] +
                               ph->sigma_photoelectric[i] +
                               ph->sigma_pair[i];
    }

    /* Incoherent scattering function S(q,Z) at JXS[2] — 21 fixed points */
    int jinc = t->jxs[1]; /* JXS[2]: JINC */
    if (jinc > 0) {
        ph->n_incoherent_ff = 21;
        ph->incoherent_ff = xss_copy(t, jinc, 21);
        /* Standard MCNP momentum transfer grid for incoherent (inverse Angstroms) */
        ph->incoherent_momentum = malloc(21 * sizeof(double));
        if (ph->incoherent_momentum) {
            static const double jinc_grid[21] = {
                0.0, 0.005, 0.01, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5,
                0.6, 0.7, 0.8, 0.9, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 8.0
            };
            memcpy(ph->incoherent_momentum, jinc_grid, 21 * sizeof(double));
        }
    }

    /* Coherent form factor at JXS[3] — 55 momentum values + 55 integrated FF */
    int jcoh = t->jxs[2]; /* JXS[3]: JCOH */
    if (jcoh > 0) {
        ph->n_coherent_ff = 55;
        ph->coherent_momentum = malloc(55 * sizeof(double));
        if (ph->coherent_momentum) {
            /* Standard MCNP momentum transfer grid for coherent (inverse Angstroms) */
            static const double jcoh_grid[55] = {
                0.0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.08, 0.10, 0.12,
                0.15, 0.18, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.55,
                0.60, 0.70, 0.80, 0.90, 1.00, 1.10, 1.20, 1.30, 1.40, 1.50,
                1.60, 1.70, 1.80, 1.90, 2.00, 2.20, 2.40, 2.60, 2.80, 3.00,
                3.20, 3.40, 3.60, 3.80, 4.00, 4.20, 4.40, 4.60, 4.80, 5.00,
                5.50, 6.00, 7.00, 8.00, 10.0
            };
            memcpy(ph->coherent_momentum, jcoh_grid, 55 * sizeof(double));
        }
        ph->coherent_ff = xss_copy(t, jcoh, 55);
        ph->coherent_ff_cumulative = xss_copy(t, jcoh + 55, 55);
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

    if (N <= 0 || M <= 0 || N > 10000 || M > 1000) return ALEA_OK;

    nuc->urr = calloc(1, sizeof(alea_nuc_urr_t));
    if (!nuc->urr) return ALEA_ERR_OUT_OF_MEMORY;

    nuc->urr->n_energies = N;
    nuc->urr->n_bands = M;
    nuc->urr->interp = interp;
    nuc->urr->inelastic_flag = ilf;
    nuc->urr->absorption_flag = ioa;

    nuc->urr->energy = xss_copy(t, loc + 6, N);
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
            filepath_alloc = malloc(len + 1);
            if (!filepath_alloc) return NULL;
            filepath = filepath_alloc;
        }
        memcpy(filepath, entry->filename, len + 1);
    } else if (xsdir->datapath[0]) {
        size_t len = strlen(xsdir->datapath) + 1 + strlen(entry->filename);
        if (len >= sizeof(filepath_buf)) {
            filepath_alloc = malloc(len + 1);
            if (!filepath_alloc) return NULL;
            filepath = filepath_alloc;
        }
        snprintf(filepath, len + 1, "%s/%s",
                 xsdir->datapath, entry->filename);
    } else {
        size_t len = strlen(entry->filename);
        if (len >= sizeof(filepath_buf)) {
            filepath_alloc = malloc(len + 1);
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
    alea_nuc_nuclide_t* nuc = calloc(1, sizeof(*nuc));
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
        if (err == ALEA_OK) err = decode_urr(nuc, &raw);
        if (err == ALEA_OK) {
            /* Decode angular and energy distributions */
            alea_nuc_decode_all_angular(nuc);
            alea_nuc_decode_all_energy(nuc);

            /* Propagate Q-values to Law 66 (N-body phase space)
             * energy distributions, which need Q for the available
             * energy calculation but don't store it in the ACE data. */
            for (int i = 0; i < nuc->n_reactions; i++) {
                alea_nuc_energy_dist_t* ed = nuc->reactions[i].energy;
                while (ed) {
                    if (ed->law == ALEA_NUC_ELAW_NBODY)
                        ed->level_Q = nuc->reactions[i].q_value;
                    ed = ed->next;
                }
            }

            /* Build MT → reaction index lookup table */
            nuc->mt_to_rxn = malloc(ALEA_NUC_MT_TABLE_SIZE * sizeof(int));
            if (nuc->mt_to_rxn) {
                memset(nuc->mt_to_rxn, 0xFF,
                       ALEA_NUC_MT_TABLE_SIZE * sizeof(int)); /* -1 */
                for (int i = 0; i < nuc->n_reactions; i++) {
                    int mt = nuc->reactions[i].mt;
                    if (mt >= 0 && mt < ALEA_NUC_MT_TABLE_SIZE)
                        nuc->mt_to_rxn[mt] = i;
                }
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

    ALEA_LOG_INFO("loaded %s (Z=%d A=%d, %d energies, %d reactions)",
                 zaid, nuc->Z, nuc->A, nuc->n_energies, nuc->n_reactions);

    return nuc;
}
