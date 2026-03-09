// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

#include "nuclear_internal.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

void alea_nuc_nuclide_free(alea_nuc_nuclide_t* nuc) {
    if (!nuc) return;

    free(nuc->energy);
    free(nuc->sigma_total);
    free(nuc->sigma_abs);
    free(nuc->sigma_elastic);
    free(nuc->heating);

    /* Free elastic angular distribution */
    if (nuc->elastic_angular) {
        for (int j = 0; j < nuc->elastic_angular->n_energies; j++) {
            free(nuc->elastic_angular->data[j].cosine);
            free(nuc->elastic_angular->data[j].pdf);
            free(nuc->elastic_angular->data[j].cdf);
        }
        free(nuc->elastic_angular->energy);
        free(nuc->elastic_angular->data);
        free(nuc->elastic_angular);
    }

    for (int i = 0; i < nuc->n_reactions; i++) {
        alea_nuc_reaction_t* r = &nuc->reactions[i];
        free(r->xs);

        if (r->angular) {
            for (int j = 0; j < r->angular->n_energies; j++) {
                free(r->angular->data[j].cosine);
                free(r->angular->data[j].pdf);
                free(r->angular->data[j].cdf);
            }
            free(r->angular->energy);
            free(r->angular->data);
            free(r->angular);
        }

        alea_nuc_energy_dist_t* ed = r->energy;
        while (ed) {
            alea_nuc_energy_dist_t* next = ed->next;
            free(ed->nbt);
            free(ed->interp);
            free(ed->energy);
            free(ed->probability);
            free(ed->data);
            free(ed->temp_energy);
            free(ed->temp_T);
            free(ed->temp_C);
            /* Free tabular data (law 4/44) */
            if (ed->tab.n_ein > 0) {
                for (int j = 0; j < ed->tab.n_ein; j++) {
                    free(ed->tab.eout ? ed->tab.eout[j] : NULL);
                    free(ed->tab.pdf ? ed->tab.pdf[j] : NULL);
                    free(ed->tab.cdf ? ed->tab.cdf[j] : NULL);
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
            ed = next;
        }
    }
    free(nuc->reactions);

    if (nuc->fission) {
        alea_nuc_nu_bar_t* bars[] = {
            nuc->fission->total,
            nuc->fission->prompt,
            nuc->fission->delayed
        };
        for (int i = 0; i < 3; i++) {
            if (bars[i]) {
                free(bars[i]->coeffs);
                free(bars[i]->energy);
                free(bars[i]->nu);
                free(bars[i]);
            }
        }
        free(nuc->fission);
    }

    if (nuc->urr) {
        free(nuc->urr->energy);
        free(nuc->urr->table);
        free(nuc->urr);
    }

    free(nuc->mt_to_rxn);

    if (nuc->photon) {
        free(nuc->photon->energy);
        free(nuc->photon->sigma_incoherent);
        free(nuc->photon->sigma_coherent);
        free(nuc->photon->sigma_photoelectric);
        free(nuc->photon->sigma_pair);
        free(nuc->photon->heating);
        free(nuc->photon->ln_energy);
        free(nuc->photon->ln_sigma_incoherent);
        free(nuc->photon->ln_sigma_coherent);
        free(nuc->photon->ln_sigma_photoelectric);
        free(nuc->photon->ln_sigma_pair);
        free(nuc->photon->incoherent_momentum);
        free(nuc->photon->incoherent_ff);
        free(nuc->photon->coherent_momentum);
        free(nuc->photon->coherent_ff);
        free(nuc->photon->coherent_ff_cumulative);
        free(nuc->photon);
    }

    alea_nuc_ace_free(&nuc->raw);
    free(nuc);
}

alea_error_t alea_nuc_parse_zaid(const char* zaid, int* Z, int* A, int* meta,
                            alea_nuc_table_type_t* type) {
    if (!zaid) return ALEA_ERR_NULL_ARG;

    /* Parse ZZAAA or ZZZAAA before the dot */
    const char* dot = strchr(zaid, '.');
    if (!dot) return ALEA_ERR_PARSE_ERROR;

    /* Extract numeric part: ZZZAAA (up to 6 digits) */
    char buf[16];
    size_t len = (size_t)(dot - zaid);
    if (len == 0 || len >= sizeof(buf)) return ALEA_ERR_PARSE_ERROR;
    memcpy(buf, zaid, len);
    buf[len] = '\0';

    int za = atoi(buf);
    if (za <= 0) return ALEA_ERR_PARSE_ERROR;

    if (Z) *Z = za / 1000;
    if (A) *A = za % 1000;
    if (meta) *meta = 0; /* TODO: metastable parsing from suffix */

    /* Parse table type from suffix character */
    const char* suffix = dot + 1;
    /* Skip digits to find the letter */
    while (*suffix >= '0' && *suffix <= '9') suffix++;

    if (type) {
        switch (*suffix) {
        case 'c': *type = ALEA_NUC_TABLE_CONTINUOUS_NEUTRON; break;
        case 'p': *type = ALEA_NUC_TABLE_PHOTOATOMIC; break;
        case 'u': *type = ALEA_NUC_TABLE_PHOTONUCLEAR; break;
        case 't': *type = ALEA_NUC_TABLE_THERMAL_SAB; break;
        case 'e': *type = ALEA_NUC_TABLE_ELECTRON; break;
        default:  return ALEA_ERR_PARSE_ERROR;
        }
    }

    return ALEA_OK;
}

