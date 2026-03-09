// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/**
 * @file angular.c
 * @brief Angular distribution decoding from ACE AND block
 *
 * ACE AND block (JXS[8]):
 *   For elastic (first entry) and each reaction in MTR:
 *     LAND[i] = locator (0 = isotropic, >0 = offset into AND block)
 *
 *   At each locator:
 *     NE = number of incident energies
 *     E[NE] = incident energy grid
 *     LC[NE] = locators for angular data at each energy
 *       LC > 0: tabular (32 equiprobable bins or full tabular)
 *       LC = 0: isotropic at this energy
 *       LC < 0: |LC| points to 32 equiprobable cosine bins
 *
 *   Tabular format (LC > 0, at JXS[8] + LC - 1):
 *     JJ = interpolation flag (histogram=1, lin-lin=2)
 *     NP = number of cosine bins
 *     CSOUT[NP] = cosine values
 *     PDF[NP] = probability density
 *     CDF[NP] = cumulative distribution
 *
 *   Equiprobable format (LC < 0, at JXS[8] + |LC| - 1):
 *     33 cosine values defining 32 equiprobable bins
 */

#include "nuclear_internal.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

/**
 * Decode angular distribution for one reaction from the AND block.
 * Returns NULL if isotropic (locator = 0).
 */
alea_nuc_angular_dist_t* alea_nuc_decode_angular(const alea_nuc_ace_table_t* t, int land_loc) {
    if (land_loc <= 0) return NULL; /* isotropic */

    int and_base = t->jxs[8]; /* JXS[9]: AND data block (0-indexed: jxs[8]) */
    int abs_loc = and_base + land_loc - 1;

    int ne = xss_int(t, abs_loc);
    if (ne <= 0 || ne > 100000) return NULL; /* sanity check */

    alea_nuc_angular_dist_t* ang = calloc(1, sizeof(*ang));
    if (!ang) return NULL;

    ang->n_energies = ne;
    ang->energy = malloc((size_t)ne * sizeof(double));
    ang->data = calloc((size_t)ne, sizeof(alea_nuc_angular_point_t));
    if (!ang->energy || !ang->data) goto fail;

    /* Read incident energy grid */
    for (int i = 0; i < ne; i++)
        ang->energy[i] = xss(t, abs_loc + 1 + i);

    /* Read LC locators */
    for (int i = 0; i < ne; i++) {
        int lc = xss_int(t, abs_loc + 1 + ne + i);
        alea_nuc_angular_point_t* pt = &ang->data[i];

        if (lc == 0) {
            /* Isotropic at this energy */
            pt->type = ALEA_NUC_ANG_ISOTROPIC;
            pt->n_cosines = 0;
            pt->cosine = NULL;
            pt->pdf = NULL;
            pt->cdf = NULL;
        } else {
            /* LC != 0: angular data at |LC| offset in AND block.
             *
             * ACE format (LA-UR-03-1987) says LC < 0 means 32 equiprobable
             * cosine bins (33 raw values, no header). In practice, many ACE
             * files (FENDL, JEFF) use the same tabular format (JJ, NP,
             * cosines, PDF, CDF) regardless of sign. The negative sign
             * indicates shared/reused data from another energy point.
             *
             * Detection: read the first value at the locator. If it is
             * a valid interpolation flag (1 or 2), treat as tabular.
             * Otherwise treat as 32 equiprobable bins. */
            int dloc = and_base + abs(lc) - 1;
            int jj = xss_int(t, dloc);

            if (jj == 1 || jj == 2) {
                /* Tabular format: JJ, NP, cosines[NP], PDF[NP], CDF[NP] */
                int np = xss_int(t, dloc + 1);
                pt->type = ALEA_NUC_ANG_TABULAR;
                pt->n_cosines = np;
                pt->cosine = malloc((size_t)np * sizeof(double));
                pt->pdf = malloc((size_t)np * sizeof(double));
                pt->cdf = malloc((size_t)np * sizeof(double));
                if (!pt->cosine || !pt->pdf || !pt->cdf) goto fail;

                for (int j = 0; j < np; j++)
                    pt->cosine[j] = xss(t, dloc + 2 + j);
                for (int j = 0; j < np; j++)
                    pt->pdf[j] = xss(t, dloc + 2 + np + j);
                for (int j = 0; j < np; j++)
                    pt->cdf[j] = xss(t, dloc + 2 + 2 * np + j);
            } else {
                /* 32 equiprobable cosine bins (33 boundary values, no header) */
                pt->type = ALEA_NUC_ANG_EQUIPROBABLE;
                pt->n_cosines = 33;
                pt->cosine = malloc(33 * sizeof(double));
                if (!pt->cosine) goto fail;
                for (int j = 0; j < 33; j++)
                    pt->cosine[j] = xss(t, dloc + j);
                pt->pdf = NULL;
                pt->cdf = NULL;
            }
        }
    }

    return ang;

fail:
    if (ang) {
        if (ang->data) {
            for (int i = 0; i < ne; i++) {
                free(ang->data[i].cosine);
                free(ang->data[i].pdf);
                free(ang->data[i].cdf);
            }
            free(ang->data);
        }
        free(ang->energy);
        free(ang);
    }
    return NULL;
}

/**
 * Decode angular distributions for elastic and all reactions.
 * Called during nuclide loading.
 */
void alea_nuc_decode_all_angular(alea_nuc_nuclide_t* nuc) {
    const alea_nuc_ace_table_t* t = &nuc->raw;

    int land_base = t->jxs[7]; /* JXS[8]: LAND locator array */
    if (land_base <= 0) return;
    if (t->jxs[8] <= 0) return; /* JXS[9]: AND data block required */

    int nr = nuc->n_reactions;

    /* LAND array at JXS[8]: first entry is elastic, then NR reactions */

    /* Elastic angular distribution (index 0 in LAND) */
    int land_elastic = xss_int(t, land_base);
    nuc->elastic_angular = alea_nuc_decode_angular(t, land_elastic);

    /* Non-elastic reactions (indices 1..NR in LAND) */
    for (int i = 0; i < nr; i++) {
        int land_loc = xss_int(t, land_base + 1 + i);
        nuc->reactions[i].angular = alea_nuc_decode_angular(t, land_loc);
    }
}
