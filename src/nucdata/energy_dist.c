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
 *     5  - General evaporation
 *     7  - Maxwell fission spectrum
 *     9  - Evaporation spectrum
 *     11 - Watt fission spectrum
 *     44 - Kalbach-Mann
 *     61 - correlated energy-angle
 *     66 - N-body phase space
 *     67 - laboratory angle-energy
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
    free(ed->temp_nbt);
    free(ed->temp_interp);
    free(ed->temp_C);
    free(ed->general_evap_x);
    free(ed->watt_b_energy);
    free(ed->watt_b_nbt);
    free(ed->watt_b_interp);
    free(ed->tab.nbt);
    free(ed->tab.interp);
    if (ed->tab.n_ein > 0) {
        for (int j = 0; j < ed->tab.n_ein; j++) {
            if (ed->tab.eout) free(ed->tab.eout[j]);
            if (ed->tab.pdf) free(ed->tab.pdf[j]);
            if (ed->tab.cdf) free(ed->tab.cdf[j]);
            if (ed->tab.precompound_r) free(ed->tab.precompound_r[j]);
            if (ed->tab.precompound_a) free(ed->tab.precompound_a[j]);
            if (ed->tab.ang_lc) free(ed->tab.ang_lc[j]);
            if (ed->tab.correlated_mu) {
                for (int k = 0; k < ed->tab.n_eout[j]; k++) {
                    free(ed->tab.correlated_mu[j][k].cosine);
                    free(ed->tab.correlated_mu[j][k].pdf);
                    free(ed->tab.correlated_mu[j][k].cdf);
                }
                free(ed->tab.correlated_mu[j]);
            }
        }
        free(ed->tab.ein);
        free(ed->tab.interpolation);
        free(ed->tab.n_discrete);
        free(ed->tab.n_eout);
        free(ed->tab.eout);
        free(ed->tab.pdf);
        free(ed->tab.cdf);
        free(ed->tab.precompound_r);
        free(ed->tab.precompound_a);
        free(ed->tab.ang_lc);
        free(ed->tab.correlated_mu);
    }
    if (ed->law67.incident) {
        for (int i = 0; i < ed->law67.n_ein; i++) {
            alea_nuc_law67_incident_t* incident = &ed->law67.incident[i];
            if (incident->spectrum)
                for (int j = 0; j < incident->n_cosines; j++) {
                    free(incident->spectrum[j].energy);
                    free(incident->spectrum[j].pdf);
                    free(incident->spectrum[j].cdf);
                }
            free(incident->cosine);
            free(incident->spectrum);
        }
    }
    free(ed->law67.nbt);
    free(ed->law67.interp);
    free(ed->law67.ein);
    free(ed->law67.incident);
    free(ed);
}

/** Free an entire energy distribution linked list. */
void alea_nuc_energy_dist_free(alea_nuc_energy_dist_t* ed) {
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
alea_nuc_energy_dist_t* alea_nuc_decode_energy_dist_base(
    const alea_nuc_ace_table_t* t, int ldlw_loc, int dlw_base) {
    if (ldlw_loc <= 0) return NULL;
    if (dlw_base <= 0) return NULL;

    int abs_loc = xss_relative_loc(t, dlw_base, ldlw_loc);
    if (abs_loc == 0) return NULL;

    alea_nuc_energy_dist_t* head = NULL;
    alea_nuc_energy_dist_t* prev = NULL;
    int max_laws = 10; /* safety limit on chained laws */

    while (abs_loc > 0 && abs_loc <= t->xss_length && max_laws-- > 0) {
        alea_nuc_energy_dist_t* ed = alea_nuc_calloc(1, sizeof(*ed));
        if (!ed) {
            xss_mark_allocation_error(t);
            alea_nuc_energy_dist_free(head);
            return NULL;
        }

        if (!xss_range_valid(t, abs_loc, 4)) {
            ed_free_node(ed);
            alea_nuc_energy_dist_free(head);
            return NULL;
        }

        int lnw = xss_int(t, abs_loc);      /* next law locator (0 = none) */
        int law = xss_int(t, abs_loc + 1);   /* law number */
        int idat = xss_int(t, abs_loc + 2);  /* data locator rel to DLW */

        ed->law = (alea_nuc_energy_law_t)law;

        /* Applicability region */
        int pos = abs_loc + 3;
        int nr_interp = xss_int(t, pos);
        pos++;

        if (nr_interp < 0 || nr_interp > 100000 ||
            !xss_range_valid(t, pos, 2 * nr_interp + 1)) {
            xss_mark_corrupt(t);
            ed_free_node(ed);
            alea_nuc_energy_dist_free(head);
            return NULL;
        }

        if (nr_interp > 0) {
            ed->n_regions = nr_interp;
            ed->nbt = alea_nuc_malloc((size_t)nr_interp * sizeof(int));
            ed->interp = alea_nuc_malloc((size_t)nr_interp * sizeof(int));
            if (!ed->nbt || !ed->interp) xss_mark_allocation_error(t);
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
        if (ne < 0 || ne > 100000 || !xss_range_valid(t, pos, 2 * ne)) {
            xss_mark_corrupt(t);
            ed_free_node(ed);
            alea_nuc_energy_dist_free(head);
            return NULL;
        }
        ed->n_energies = ne;

        if (ne > 0) {
            ed->energy = alea_nuc_malloc((size_t)ne * sizeof(double));
            ed->probability = alea_nuc_malloc((size_t)ne * sizeof(double));
            if (!ed->energy || !ed->probability) xss_mark_allocation_error(t);
            if (ed->energy && ed->probability) {
                for (int i = 0; i < ne; i++)
                    ed->energy[i] = xss(t, pos + i);
                for (int i = 0; i < ne; i++)
                    ed->probability[i] = xss(t, pos + ne + i);
            }
        }

        /* Decode law-specific data */
        int data_loc = xss_relative_loc(t, dlw_base, idat);
        if (data_loc == 0) {
            ed_free_node(ed);
            alea_nuc_energy_dist_free(head);
            return NULL;
        }
        switch (law) {
        case 2: /* Discrete photon energy */
            if (!xss_range_valid(t, data_loc, 2)) break;
            ed->discrete_photon_primary = xss_int(t, data_loc);
            ed->discrete_photon_energy = xss(t, data_loc + 1);
            ed->discrete_photon_awr = t->awr;
            break;

        case 3: /* Level scattering */
            if (!xss_range_valid(t, data_loc, 2)) break;
            ed->level_A = xss(t, data_loc);     /* threshold energy */
            ed->level_Q = xss(t, data_loc + 1); /* mass-ratio factor */
            break;

        case 5: /* General evaporation spectrum */
        {
            int nr2 = xss_int(t, data_loc);
            if (nr2 < 0 || data_loc >= t->xss_length ||
                nr2 > (t->xss_length - data_loc - 1) / 2) {
                xss_mark_corrupt(t);
                break;
            }
            ed->n_temp_regions = nr2;
            if (nr2 > 0) {
                ed->temp_nbt = alea_nuc_malloc((size_t)nr2 * sizeof(int));
                ed->temp_interp = alea_nuc_malloc((size_t)nr2 * sizeof(int));
                if (!ed->temp_nbt || !ed->temp_interp) {
                    xss_mark_allocation_error(t);
                    break;
                }
                for (int i = 0; i < nr2; i++) {
                    ed->temp_nbt[i] = xss_int(t, data_loc + 1 + i);
                    ed->temp_interp[i] =
                        xss_int(t, data_loc + 1 + nr2 + i);
                }
            }
            int base = data_loc + 1 + 2 * nr2;
            int nt = xss_int(t, base);
            if (nt <= 0 || nt > 100000 ||
                !xss_range_valid(t, base + 1, 2 * nt + 1)) {
                xss_mark_corrupt(t);
                break;
            }
            ed->n_temp = nt;
            ed->temp_energy = alea_nuc_malloc((size_t)nt * sizeof(double));
            ed->temp_T = alea_nuc_malloc((size_t)nt * sizeof(double));
            if (!ed->temp_energy || !ed->temp_T) {
                xss_mark_allocation_error(t);
                break;
            }
            for (int i = 0; i < nt; i++) {
                ed->temp_energy[i] = xss(t, base + 1 + i);
                ed->temp_T[i] = xss(t, base + 1 + nt + i);
            }
            int net_loc = base + 1 + 2 * nt;
            int net = xss_int(t, net_loc);
            if (net < 2 || net > 100000 ||
                !xss_range_valid(t, net_loc + 1, net)) {
                xss_mark_corrupt(t);
                break;
            }
            ed->n_general_evap = net;
            ed->general_evap_x = xss_copy(t, net_loc + 1, net);
            break;
        }

        case 7: /* Maxwell fission spectrum */
        case 9: /* Evaporation spectrum */
        {
            int nr2 = xss_int(t, data_loc);
            if (nr2 < 0 || data_loc >= t->xss_length ||
                nr2 > (t->xss_length - data_loc - 1) / 2) {
                xss_mark_corrupt(t);
                break;
            }
            ed->n_temp_regions = nr2;
            if (nr2 > 0) {
                ed->temp_nbt = alea_nuc_malloc((size_t)nr2 * sizeof(int));
                ed->temp_interp = alea_nuc_malloc((size_t)nr2 * sizeof(int));
                if (!ed->temp_nbt || !ed->temp_interp) {
                    xss_mark_allocation_error(t);
                    break;
                }
                for (int i = 0; i < nr2; i++) {
                    ed->temp_nbt[i] = xss_int(t, data_loc + 1 + i);
                    ed->temp_interp[i] = xss_int(t, data_loc + 1 + nr2 + i);
                }
            }
            int base = data_loc + 1 + 2 * nr2;
            int nt = xss_int(t, base);
            if (nt < 0 || nt > 100000 ||
                !xss_range_valid(t, base + 1, 2 * nt + 1)) {
                xss_mark_corrupt(t);
                break;
            }
            ed->n_temp = nt;
            if (nt > 0) {
                ed->temp_energy = alea_nuc_malloc((size_t)nt * sizeof(double));
                ed->temp_T = alea_nuc_malloc((size_t)nt * sizeof(double));
                if (!ed->temp_energy || !ed->temp_T) xss_mark_allocation_error(t);
                if (ed->temp_energy && ed->temp_T) {
                    for (int i = 0; i < nt; i++)
                        ed->temp_energy[i] = xss(t, base + 1 + i);
                    for (int i = 0; i < nt; i++)
                        ed->temp_T[i] = xss(t, base + 1 + nt + i);
                }
                /* Restriction energy U */
                ed->restriction_energy = xss(t, base + 1 + 2 * nt);
            }
            break;
        }

        case 11: /* Watt fission spectrum */
        {
            int nr2 = xss_int(t, data_loc);
            if (nr2 < 0 || data_loc >= t->xss_length ||
                nr2 > (t->xss_length - data_loc - 1) / 2) {
                xss_mark_corrupt(t);
                break;
            }
            ed->n_temp_regions = nr2;
            if (nr2 > 0) {
                ed->temp_nbt = alea_nuc_malloc((size_t)nr2 * sizeof(int));
                ed->temp_interp = alea_nuc_malloc((size_t)nr2 * sizeof(int));
                if (!ed->temp_nbt || !ed->temp_interp) {
                    xss_mark_allocation_error(t);
                    break;
                }
                for (int i = 0; i < nr2; i++) {
                    ed->temp_nbt[i] = xss_int(t, data_loc + 1 + i);
                    ed->temp_interp[i] = xss_int(t, data_loc + 1 + nr2 + i);
                }
            }
            int base = data_loc + 1 + 2 * nr2;
            int na = xss_int(t, base);
            if (na <= 0 || na > 100000 ||
                !xss_range_valid(t, base + 1, 2 * na + 1)) {
                xss_mark_corrupt(t);
                break;
            }
            /* Read 'a' parameter table */
            ed->n_temp = na;
            ed->temp_energy = alea_nuc_malloc((size_t)na * sizeof(double));
            ed->temp_T = alea_nuc_malloc((size_t)na * sizeof(double));
            if (!ed->temp_energy || !ed->temp_T) xss_mark_allocation_error(t);
            if (ed->temp_energy && ed->temp_T) {
                for (int i = 0; i < na; i++)
                    ed->temp_energy[i] = xss(t, base + 1 + i);
                for (int i = 0; i < na; i++)
                    ed->temp_T[i] = xss(t, base + 1 + na + i);
            }
            /* 'b' parameter table follows */
            int bbase = base + 1 + 2 * na;
            int nr3 = xss_int(t, bbase);
            if (nr3 < 0 || bbase >= t->xss_length ||
                nr3 > (t->xss_length - bbase - 1) / 2) {
                xss_mark_corrupt(t);
                break;
            }
            ed->n_watt_b_regions = nr3;
            if (nr3 > 0) {
                ed->watt_b_nbt = alea_nuc_malloc((size_t)nr3 * sizeof(int));
                ed->watt_b_interp = alea_nuc_malloc((size_t)nr3 * sizeof(int));
                if (!ed->watt_b_nbt || !ed->watt_b_interp) {
                    xss_mark_allocation_error(t);
                    break;
                }
                for (int i = 0; i < nr3; i++) {
                    ed->watt_b_nbt[i] = xss_int(t, bbase + 1 + i);
                    ed->watt_b_interp[i] = xss_int(t, bbase + 1 + nr3 + i);
                }
            }
            int bbase2 = bbase + 1 + 2 * nr3;
            int nb = xss_int(t, bbase2);
            if (nb <= 0 || nb > 100000 ||
                !xss_range_valid(t, bbase2 + 1, 2 * nb + 1)) {
                xss_mark_corrupt(t);
                break;
            }
            ed->n_watt_b = nb;
            ed->watt_b_energy = alea_nuc_malloc((size_t)nb * sizeof(double));
            ed->temp_C = alea_nuc_malloc((size_t)nb * sizeof(double));
            if (!ed->watt_b_energy || !ed->temp_C) xss_mark_allocation_error(t);
            if (ed->watt_b_energy && ed->temp_C) {
                for (int i = 0; i < nb; i++) {
                    ed->watt_b_energy[i] = xss(t, bbase2 + 1 + i);
                    ed->temp_C[i] = xss(t, bbase2 + 1 + nb + i);
                }
            }
            /* Restriction energy */
            ed->restriction_energy = xss(t, bbase2 + 1 + 2 * nb);
            break;
        }

        case 4:  /* Continuous tabular */
        case 44: /* Kalbach-Mann */
        case 61: /* Correlated energy-angle */
        {
            int pos2 = data_loc;
            int nr2 = xss_int(t, pos2);
            if (nr2 < 0 || pos2 >= t->xss_length ||
                nr2 > (t->xss_length - pos2 - 1) / 2) {
                xss_mark_corrupt(t);
                break;
            }
            ed->tab.n_regions = nr2;
            if (nr2 > 0) {
                ed->tab.nbt = alea_nuc_malloc((size_t)nr2 * sizeof(int));
                ed->tab.interp = alea_nuc_malloc((size_t)nr2 * sizeof(int));
                if (!ed->tab.nbt || !ed->tab.interp) {
                    xss_mark_allocation_error(t);
                    break;
                }
                for (int i = 0; i < nr2; i++) {
                    ed->tab.nbt[i] = xss_int(t, pos2 + 1 + i);
                    ed->tab.interp[i] = xss_int(t, pos2 + 1 + nr2 + i);
                }
            }
            pos2 += 1 + 2 * nr2; /* skip interpolation data */

            int n_ein = xss_int(t, pos2);
            pos2++;

            if (n_ein <= 0 || n_ein > 100000) break;
            if (!xss_range_valid(t, pos2, 2 * n_ein)) break;

            ed->tab.n_ein = n_ein;
            ed->tab.dlw_base = dlw_base;
            ed->tab.ein = alea_nuc_malloc((size_t)n_ein * sizeof(double));
            ed->tab.interpolation = alea_nuc_calloc((size_t)n_ein, sizeof(int));
            ed->tab.n_discrete = alea_nuc_calloc((size_t)n_ein, sizeof(int));
            ed->tab.n_eout = alea_nuc_malloc((size_t)n_ein * sizeof(int));
            ed->tab.eout = alea_nuc_calloc((size_t)n_ein, sizeof(double*));
            ed->tab.pdf = alea_nuc_calloc((size_t)n_ein, sizeof(double*));
            ed->tab.cdf = alea_nuc_calloc((size_t)n_ein, sizeof(double*));
            if (law == 44) {
                ed->tab.precompound_r = alea_nuc_calloc((size_t)n_ein, sizeof(double*));
                ed->tab.precompound_a = alea_nuc_calloc((size_t)n_ein, sizeof(double*));
            }
            if (law == 61) {
                ed->tab.ang_lc = alea_nuc_calloc((size_t)n_ein, sizeof(int*));
                ed->tab.correlated_mu = alea_nuc_calloc(
                    (size_t)n_ein, sizeof(*ed->tab.correlated_mu));
            }

            if (!ed->tab.ein || !ed->tab.interpolation ||
                !ed->tab.n_discrete || !ed->tab.n_eout || !ed->tab.eout ||
                !ed->tab.pdf || !ed->tab.cdf) {
                xss_mark_allocation_error(t);
                break;
            }
            if ((law == 44 && (!ed->tab.precompound_r || !ed->tab.precompound_a)) ||
                (law == 61 && (!ed->tab.ang_lc ||
                               !ed->tab.correlated_mu))) {
                xss_mark_allocation_error(t);
                break;
            }

            /* Read incident energies */
            for (int j = 0; j < n_ein; j++)
                ed->tab.ein[j] = xss(t, pos2 + j);

            /* Read locators */
            int* locs = alea_nuc_malloc((size_t)n_ein * sizeof(int));
            if (!locs) {
                xss_mark_allocation_error(t);
                break;
            }
            for (int j = 0; j < n_ein; j++)
                locs[j] = xss_int(t, pos2 + n_ein + j);

            /* Decode each outgoing distribution */
            for (int j = 0; j < n_ein; j++) {
                int dloc = xss_relative_loc(t, dlw_base, locs[j]);
                if (dloc == 0) continue;
                int intt = xss_int(t, dloc);
                ed->tab.interpolation[j] = intt % 10;
                ed->tab.n_discrete[j] = intt / 10;
                int np = xss_int(t, dloc + 1);

                if (np <= 0 || np > 100000) { ed->tab.n_eout[j] = 0; continue; }
                int arrays = law == 44 ? 5 : (law == 61 ? 4 : 3);
                if (!xss_range_valid(t, dloc + 2, arrays * np)) {
                    ed->tab.n_eout[j] = 0;
                    continue;
                }
                ed->tab.n_eout[j] = np;

                ed->tab.eout[j] = alea_nuc_malloc((size_t)np * sizeof(double));
                ed->tab.pdf[j] = alea_nuc_malloc((size_t)np * sizeof(double));
                ed->tab.cdf[j] = alea_nuc_malloc((size_t)np * sizeof(double));

                if (!ed->tab.eout[j] || !ed->tab.pdf[j] || !ed->tab.cdf[j]) {
                    xss_mark_allocation_error(t);
                    continue;
                }

                for (int k = 0; k < np; k++)
                    ed->tab.eout[j][k] = xss(t, dloc + 2 + k);
                for (int k = 0; k < np; k++)
                    ed->tab.pdf[j][k] = xss(t, dloc + 2 + np + k);
                for (int k = 0; k < np; k++)
                    ed->tab.cdf[j][k] = xss(t, dloc + 2 + 2 * np + k);

                /* Kalbach-Mann: R and A arrays follow CDF */
                if (law == 44 && ed->tab.precompound_r && ed->tab.precompound_a) {
                    ed->tab.precompound_r[j] = alea_nuc_malloc((size_t)np * sizeof(double));
                    ed->tab.precompound_a[j] = alea_nuc_malloc((size_t)np * sizeof(double));
                    if (!ed->tab.precompound_r[j] || !ed->tab.precompound_a[j])
                        xss_mark_allocation_error(t);
                    if (ed->tab.precompound_r[j] && ed->tab.precompound_a[j]) {
                        for (int k = 0; k < np; k++)
                            ed->tab.precompound_r[j][k] = xss(t, dloc + 2 + 3 * np + k);
                        for (int k = 0; k < np; k++)
                            ed->tab.precompound_a[j][k] = xss(t, dloc + 2 + 4 * np + k);
                    }
                }

                /* Law 61: angular locators follow CDF */
                if (law == 61 && ed->tab.ang_lc) {
                    ed->tab.ang_lc[j] = alea_nuc_malloc((size_t)np * sizeof(int));
                    ed->tab.correlated_mu[j] = alea_nuc_calloc(
                        (size_t)np, sizeof(*ed->tab.correlated_mu[j]));
                    if (!ed->tab.ang_lc[j] || !ed->tab.correlated_mu[j])
                        xss_mark_allocation_error(t);
                    if (ed->tab.ang_lc[j]) {
                        for (int k = 0; k < np; k++) {
                            ed->tab.ang_lc[j][k] = xss_int(t, dloc + 2 + 3 * np + k);
                            alea_nuc_angular_point_t* point =
                                ed->tab.correlated_mu[j]
                                    ? &ed->tab.correlated_mu[j][k] : NULL;
                            if (!point) continue;
                            int locator = ed->tab.ang_lc[j][k];
                            if (locator <= 0) {
                                point->type = ALEA_NUC_ANG_ISOTROPIC;
                                continue;
                            }
                            int aloc = xss_relative_loc(t, dlw_base, locator);
                            if (aloc == 0 || !xss_range_valid(t, aloc, 2))
                                continue;
                            int aint = xss_int(t, aloc);
                            int ncos = xss_int(t, aloc + 1);
                            if ((aint != 1 && aint != 2) || ncos < 2 ||
                                ncos > 100000 ||
                                !xss_range_valid(t, aloc + 2, 3 * ncos)) {
                                xss_mark_corrupt(t);
                                continue;
                            }
                            point->type = ALEA_NUC_ANG_TABULAR;
                            point->interpolation = aint;
                            point->n_cosines = ncos;
                            point->cosine = xss_copy(t, aloc + 2, ncos);
                            point->pdf = xss_copy(t, aloc + 2 + ncos, ncos);
                            point->cdf = xss_copy(t, aloc + 2 + 2 * ncos, ncos);
                        }
                    }
                }
            }
            free(locs);
            break;
        }

        case 66: /* N-body phase space */
        {
            if (!xss_range_valid(t, data_loc, 2)) break;
            ed->nbody_particles = xss_int(t, data_loc); /* NPSX */
            ed->nbody_total_mass = xss(t, data_loc + 1); /* Ap */
            ed->nbody_target_awr = t->awr;
            break;
        }

        case 67: /* Laboratory angle-energy law */
        {
            int cursor = data_loc;
            int nr = xss_int(t, cursor);
            if (nr < 0 || nr > 100000 ||
                !xss_range_valid(t, cursor + 1, 2 * nr + 1)) {
                xss_mark_corrupt(t);
                break;
            }
            ed->law67.n_regions = nr;
            if (nr > 0) {
                ed->law67.nbt = alea_nuc_malloc((size_t)nr * sizeof(int));
                ed->law67.interp = alea_nuc_malloc((size_t)nr * sizeof(int));
                if (!ed->law67.nbt || !ed->law67.interp) {
                    xss_mark_allocation_error(t);
                    break;
                }
                for (int i = 0; i < nr; i++) {
                    ed->law67.nbt[i] = xss_int(t, cursor + 1 + i);
                    ed->law67.interp[i] = xss_int(t, cursor + 1 + nr + i);
                }
            }
            cursor += 1 + 2 * nr;
            int ne = xss_int(t, cursor++);
            if (ne <= 0 || ne > 100000 ||
                !xss_range_valid(t, cursor, 2 * ne)) {
                xss_mark_corrupt(t);
                break;
            }
            ed->law67.n_ein = ne;
            ed->law67.ein = alea_nuc_malloc((size_t)ne * sizeof(double));
            ed->law67.incident = alea_nuc_calloc(
                (size_t)ne, sizeof(*ed->law67.incident));
            int* locators = alea_nuc_malloc((size_t)ne * sizeof(int));
            if (!ed->law67.ein || !ed->law67.incident || !locators) {
                xss_mark_allocation_error(t);
                free(locators);
                break;
            }
            for (int i = 0; i < ne; i++) {
                ed->law67.ein[i] = xss(t, cursor + i);
                locators[i] = xss_int(t, cursor + ne + i);
            }
            for (int i = 0; i < ne; i++) {
                int aloc = xss_relative_loc(t, dlw_base, locators[i]);
                if (aloc == 0 || !xss_range_valid(t, aloc, 2)) {
                    xss_mark_corrupt(t);
                    continue;
                }
                alea_nuc_law67_incident_t* incident =
                    &ed->law67.incident[i];
                incident->interpolation = xss_int(t, aloc);
                int nmu = xss_int(t, aloc + 1);
                if ((incident->interpolation != 1 &&
                     incident->interpolation != 2) ||
                    nmu < 2 || nmu > 100000 ||
                    !xss_range_valid(t, aloc + 2, 2 * nmu)) {
                    xss_mark_corrupt(t);
                    continue;
                }
                incident->n_cosines = nmu;
                incident->cosine = xss_copy(t, aloc + 2, nmu);
                incident->spectrum = alea_nuc_calloc(
                    (size_t)nmu, sizeof(*incident->spectrum));
                if (!incident->cosine || !incident->spectrum) {
                    xss_mark_allocation_error(t);
                    continue;
                }
                for (int j = 0; j < nmu; j++) {
                    int locator = xss_int(t, aloc + 2 + nmu + j);
                    int eloc = xss_relative_loc(t, dlw_base, locator);
                    if (eloc == 0 || !xss_range_valid(t, eloc, 2)) {
                        xss_mark_corrupt(t);
                        continue;
                    }
                    alea_nuc_law67_energy_t* spectrum =
                        &incident->spectrum[j];
                    spectrum->interpolation = xss_int(t, eloc);
                    int np = xss_int(t, eloc + 1);
                    if ((spectrum->interpolation != 1 &&
                         spectrum->interpolation != 2) ||
                        np < 2 || np > 100000 ||
                        !xss_range_valid(t, eloc + 2, 3 * np)) {
                        xss_mark_corrupt(t);
                        continue;
                    }
                    spectrum->n_points = np;
                    spectrum->energy = xss_copy(t, eloc + 2, np);
                    spectrum->pdf = xss_copy(t, eloc + 2 + np, np);
                    spectrum->cdf = xss_copy(t, eloc + 2 + 2 * np, np);
                }
            }
            free(locators);
            break;
        }

        default:
            /* Unsupported law: store data location for debugging */
            ed->data_length = data_loc;
            break;
        }

        if (t->allocation_error || t->decode_error) {
            ed_free_node(ed);
            alea_nuc_energy_dist_free(head);
            return NULL;
        }

        /* Link list */
        if (prev)
            prev->next = ed;
        else
            head = ed;
        prev = ed;

        /* Next law? */
        if (lnw > 0) {
            abs_loc = xss_relative_loc(t, dlw_base, lnw);
            if (abs_loc == 0) {
                alea_nuc_energy_dist_free(head);
                return NULL;
            }
        } else
            abs_loc = 0;
    }

    if (abs_loc > 0) {
        xss_mark_corrupt(t);
        alea_nuc_energy_dist_free(head);
        return NULL;
    }

    return head;
}

alea_nuc_energy_dist_t* alea_nuc_decode_energy_dist(
    const alea_nuc_ace_table_t* t, int ldlw_loc) {
    return alea_nuc_decode_energy_dist_base(t, ldlw_loc, t->jxs[10]);
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
        if (ldlw_base > t->xss_length) return;
        /* Check if the first value looks like a valid locator. */
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
        for (alea_nuc_energy_dist_t* ed = nuc->reactions[i].energy;
             ed; ed = ed->next) {
            if (ed->law == ALEA_NUC_ELAW_NBODY) {
                ed->nbody_target_awr = nuc->awr;
                ed->nbody_q_value = nuc->reactions[i].q_value;
            }
        }
    }
}
