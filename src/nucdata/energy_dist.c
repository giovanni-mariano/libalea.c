// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/**
 * @file energy_dist.c
 * @brief Energy distribution decoding from ACE DLW block
 *
 * ACE DLW block (JXS[9]):
 *   LDLW[NR] = locators for each reaction's energy distribution
 *
 *   At each locator (relative to JXS[10]):
 *     LNW = next law locator (0 if none)
 *     LAW = law number
 *     IDAT = locator for law data (relative to JXS[10])
 *     NR = number of interpolation regions for applicability
 *     NBT[NR], INT[NR] = interpolation breakpoints/types
 *     NE = number of energy points
 *     E[NE] = incident energies
 *     P[NE] = applicability probability
 *
 *   Supported laws:
 *     3  - Level scattering (Q-value shift)
 *     4  - Continuous tabular
 *     7  - Maxwell fission spectrum
 *     9  - Evaporation spectrum
 *     11 - Watt fission spectrum
 *     44 - Kalbach-Mann
 */

#include "nuclear_internal.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

/** Free a single energy distribution node (not its ->next chain). */
static void ed_free_node(alea_nuc_energy_dist_t* ed) {
    if (!ed) return;
    free(ed->nbt);
    free(ed->interp);
    free(ed->energy);
    free(ed->probability);
    free(ed->data);
    free(ed->temp_energy);
    free(ed->temp_T);
    free(ed->temp_C);
    if (ed->tab.n_ein > 0) {
        for (int j = 0; j < ed->tab.n_ein; j++) {
            if (ed->tab.eout) free(ed->tab.eout[j]);
            if (ed->tab.pdf) free(ed->tab.pdf[j]);
            if (ed->tab.cdf) free(ed->tab.cdf[j]);
            if (ed->tab.precompound_r) free(ed->tab.precompound_r[j]);
            if (ed->tab.precompound_a) free(ed->tab.precompound_a[j]);
            if (ed->tab.ang_lc) free(ed->tab.ang_lc[j]);
        }
        free(ed->tab.ein);
        free(ed->tab.n_eout);
        free(ed->tab.eout);
        free(ed->tab.pdf);
        free(ed->tab.cdf);
        free(ed->tab.precompound_r);
        free(ed->tab.precompound_a);
        free(ed->tab.ang_lc);
    }
    free(ed);
}

/** Free an entire energy distribution linked list. */
static void ed_free_chain(alea_nuc_energy_dist_t* ed) {
    while (ed) {
        alea_nuc_energy_dist_t* next = ed->next;
        ed_free_node(ed);
        ed = next;
    }
}

/**
 * Decode energy distribution for one reaction.
 * Returns linked list of energy distribution laws.
 */
alea_nuc_energy_dist_t* alea_nuc_decode_energy_dist(const alea_nuc_ace_table_t* t, int ldlw_loc) {
    if (ldlw_loc <= 0) return NULL;

    int dlw_base = t->jxs[10]; /* JXS[11]: DLW data block (0-indexed: jxs[10]) */
    if (dlw_base <= 0) return NULL;

    int abs_loc = dlw_base + ldlw_loc - 1;

    alea_nuc_energy_dist_t* head = NULL;
    alea_nuc_energy_dist_t* prev = NULL;
    int max_laws = 10; /* safety limit on chained laws */

    while (abs_loc > 0 && abs_loc <= t->xss_length && max_laws-- > 0) {
        alea_nuc_energy_dist_t* ed = calloc(1, sizeof(*ed));
        if (!ed) { ed_free_chain(head); return NULL; }

        int lnw = xss_int(t, abs_loc);      /* next law locator (0 = none) */
        int law = xss_int(t, abs_loc + 1);   /* law number */
        int idat = xss_int(t, abs_loc + 2);  /* data locator rel to DLW */

        ed->law = (alea_nuc_energy_law_t)law;

        /* Applicability region */
        int pos = abs_loc + 3;
        int nr_interp = xss_int(t, pos);
        pos++;

        if (nr_interp > 0) {
            ed->n_regions = nr_interp;
            ed->nbt = malloc((size_t)nr_interp * sizeof(int));
            ed->interp = malloc((size_t)nr_interp * sizeof(int));
            if (ed->nbt && ed->interp) {
                for (int i = 0; i < nr_interp; i++)
                    ed->nbt[i] = xss_int(t, pos + i);
                for (int i = 0; i < nr_interp; i++)
                    ed->interp[i] = xss_int(t, pos + nr_interp + i);
            }
            pos += 2 * nr_interp;
        }

        int ne = xss_int(t, pos);
        pos++;
        ed->n_energies = ne;

        if (ne > 0) {
            ed->energy = malloc((size_t)ne * sizeof(double));
            ed->probability = malloc((size_t)ne * sizeof(double));
            if (ed->energy && ed->probability) {
                for (int i = 0; i < ne; i++)
                    ed->energy[i] = xss(t, pos + i);
                for (int i = 0; i < ne; i++)
                    ed->probability[i] = xss(t, pos + ne + i);
            }
        }

        /* Decode law-specific data */
        int data_loc = dlw_base + idat - 1;

        switch (law) {
        case 3: /* Level scattering */
            ed->level_A = xss(t, data_loc);     /* (A+1)/A factor */
            ed->level_Q = xss(t, data_loc + 1); /* Q-value */
            break;

        case 7: /* Maxwell fission spectrum */
        case 9: /* Evaporation spectrum */
        {
            int nr2 = xss_int(t, data_loc);
            int base = data_loc + 1 + 2 * nr2;
            int nt = xss_int(t, base);
            ed->n_temp = nt;
            if (nt > 0) {
                ed->temp_energy = malloc((size_t)nt * sizeof(double));
                ed->temp_T = malloc((size_t)nt * sizeof(double));
                if (ed->temp_energy && ed->temp_T) {
                    for (int i = 0; i < nt; i++)
                        ed->temp_energy[i] = xss(t, base + 1 + i);
                    for (int i = 0; i < nt; i++)
                        ed->temp_T[i] = xss(t, base + 1 + nt + i);
                }
                /* Restriction energy U */
                ed->level_Q = xss(t, base + 1 + 2 * nt); /* using level_Q for U */
            }
            break;
        }

        case 11: /* Watt fission spectrum */
        {
            int nr2 = xss_int(t, data_loc);
            int base = data_loc + 1 + 2 * nr2;
            int na = xss_int(t, base);
            /* Read 'a' parameter table */
            ed->n_temp = na;
            ed->temp_energy = malloc((size_t)na * sizeof(double));
            ed->temp_T = malloc((size_t)na * sizeof(double));
            if (ed->temp_energy && ed->temp_T) {
                for (int i = 0; i < na; i++)
                    ed->temp_energy[i] = xss(t, base + 1 + i);
                for (int i = 0; i < na; i++)
                    ed->temp_T[i] = xss(t, base + 1 + na + i);
            }
            /* 'b' parameter table follows */
            int bbase = base + 1 + 2 * na;
            int nr3 = xss_int(t, bbase);
            int bbase2 = bbase + 1 + 2 * nr3;
            int nb = xss_int(t, bbase2);
            /* Read b parameter (store as temp_C for first entry) */
            ed->temp_C = malloc((size_t)nb * sizeof(double));
            if (ed->temp_C) {
                for (int i = 0; i < nb; i++)
                    ed->temp_C[i] = xss(t, bbase2 + 1 + nb + i);
            }
            /* Restriction energy */
            ed->level_Q = xss(t, bbase2 + 1 + 2 * nb);
            break;
        }

        case 4:  /* Continuous tabular */
        case 44: /* Kalbach-Mann */
        case 61: /* Correlated energy-angle */
        {
            int pos2 = data_loc;
            int nr2 = xss_int(t, pos2);
            pos2 += 1 + 2 * nr2; /* skip interpolation data */

            int n_ein = xss_int(t, pos2);
            pos2++;

            if (n_ein <= 0 || n_ein > 100000) break;

            ed->tab.n_ein = n_ein;
            ed->tab.dlw_base = dlw_base;
            ed->tab.ein = malloc((size_t)n_ein * sizeof(double));
            ed->tab.n_eout = malloc((size_t)n_ein * sizeof(int));
            ed->tab.eout = calloc((size_t)n_ein, sizeof(double*));
            ed->tab.pdf = calloc((size_t)n_ein, sizeof(double*));
            ed->tab.cdf = calloc((size_t)n_ein, sizeof(double*));
            if (law == 44) {
                ed->tab.precompound_r = calloc((size_t)n_ein, sizeof(double*));
                ed->tab.precompound_a = calloc((size_t)n_ein, sizeof(double*));
            }
            if (law == 61) {
                ed->tab.ang_lc = calloc((size_t)n_ein, sizeof(int*));
            }

            if (!ed->tab.ein || !ed->tab.n_eout || !ed->tab.eout ||
                !ed->tab.pdf || !ed->tab.cdf) break;

            /* Read incident energies */
            for (int j = 0; j < n_ein; j++)
                ed->tab.ein[j] = xss(t, pos2 + j);

            /* Read locators */
            int* locs = malloc((size_t)n_ein * sizeof(int));
            if (!locs) break;
            for (int j = 0; j < n_ein; j++)
                locs[j] = xss_int(t, pos2 + n_ein + j);

            /* Decode each outgoing distribution */
            for (int j = 0; j < n_ein; j++) {
                int dloc = dlw_base + locs[j] - 1;
                /* int intt = xss_int(t, dloc); */ /* interpolation type */
                int np = xss_int(t, dloc + 1);

                if (np <= 0 || np > 100000) { ed->tab.n_eout[j] = 0; continue; }
                ed->tab.n_eout[j] = np;

                ed->tab.eout[j] = malloc((size_t)np * sizeof(double));
                ed->tab.pdf[j] = malloc((size_t)np * sizeof(double));
                ed->tab.cdf[j] = malloc((size_t)np * sizeof(double));

                if (!ed->tab.eout[j] || !ed->tab.pdf[j] || !ed->tab.cdf[j]) continue;

                for (int k = 0; k < np; k++)
                    ed->tab.eout[j][k] = xss(t, dloc + 2 + k);
                for (int k = 0; k < np; k++)
                    ed->tab.pdf[j][k] = xss(t, dloc + 2 + np + k);
                for (int k = 0; k < np; k++)
                    ed->tab.cdf[j][k] = xss(t, dloc + 2 + 2 * np + k);

                /* Kalbach-Mann: R and A arrays follow CDF */
                if (law == 44 && ed->tab.precompound_r && ed->tab.precompound_a) {
                    ed->tab.precompound_r[j] = malloc((size_t)np * sizeof(double));
                    ed->tab.precompound_a[j] = malloc((size_t)np * sizeof(double));
                    if (ed->tab.precompound_r[j] && ed->tab.precompound_a[j]) {
                        for (int k = 0; k < np; k++)
                            ed->tab.precompound_r[j][k] = xss(t, dloc + 2 + 3 * np + k);
                        for (int k = 0; k < np; k++)
                            ed->tab.precompound_a[j][k] = xss(t, dloc + 2 + 4 * np + k);
                    }
                }

                /* Law 61: angular locators follow CDF */
                if (law == 61 && ed->tab.ang_lc) {
                    ed->tab.ang_lc[j] = malloc((size_t)np * sizeof(int));
                    if (ed->tab.ang_lc[j]) {
                        for (int k = 0; k < np; k++)
                            ed->tab.ang_lc[j][k] = xss_int(t, dloc + 2 + 3 * np + k);
                    }
                }
            }
            free(locs);
            break;
        }

        case 66: /* N-body phase space */
        {
            ed->tab.n_ein = xss_int(t, data_loc);   /* NPSX */
            ed->level_A = xss(t, data_loc + 1);      /* Ap (total mass) */
            break;
        }

        default:
            /* Unsupported law: store data location for debugging */
            ed->data_length = data_loc;
            break;
        }

        /* Link list */
        if (prev)
            prev->next = ed;
        else
            head = ed;
        prev = ed;

        /* Next law? */
        if (lnw > 0)
            abs_loc = dlw_base + lnw - 1;
        else
            abs_loc = 0;
    }

    return head;
}

/**
 * Decode energy distributions for all reactions.
 */
void alea_nuc_decode_all_energy(alea_nuc_nuclide_t* nuc) {
    const alea_nuc_ace_table_t* t = &nuc->raw;

    int ldlw_base = t->jxs[9]; /* JXS[10]: LDLW locator array (0-indexed: jxs[9]) */
    int dlw_base = t->jxs[10]; /* JXS[11]: DLW data block */
    if (ldlw_base <= 0 || dlw_base <= 0) return;
    /* If LDLW == DLW, there is no separate locator array */
    if (ldlw_base == dlw_base && nuc->n_reactions > 0) {
        /* Check if first value looks like a valid locator (small integer).
         * Use the bounds-checked accessor: ldlw_base may exceed xss_length. */
        double first = xss(t, ldlw_base);
        if (first != (int)first || first < 0 || first > t->xss_length)
            return; /* not valid locators */
    }

    for (int i = 0; i < nuc->n_reactions; i++) {
        /* Only decode if reaction produces secondary particles (TYR != 0) */
        if (nuc->reactions[i].ty == 0) continue;

        int ldlw = xss_int(t, ldlw_base + i);
        if (ldlw <= 0) continue;
        nuc->reactions[i].energy = alea_nuc_decode_energy_dist(t, ldlw);
    }
}
