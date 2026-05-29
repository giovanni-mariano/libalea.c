// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/**
 * @file reaction.c
 * @brief Reaction classification, neutron yield, and nu-bar evaluation
 *
 * Data-level functions that classify reactions by MT number,
 * evaluate neutron yields from the TYR field, and look up
 * average fission neutron multiplicity.
 */

#include "nuclear_internal.h"
#include <math.h>
#include <stdlib.h>

/* ============================================================================
 * REACTION CLASSIFICATION
 * ============================================================================ */

alea_nuc_reaction_class_t alea_nuc_reaction_classify(int mt) {
    /* Elastic */
    if (mt == 2) return ALEA_NUC_RXN_SCATTER;

    /* Inelastic levels and continuum */
    if (mt >= 51 && mt <= 91) return ALEA_NUC_RXN_SCATTER;

    /* Fission */
    if (mt == 18 || mt == 19 || mt == 20 || mt == 21 || mt == 38)
        return ALEA_NUC_RXN_MULTIPLY;

    /* (n,Xn) reactions — produce multiple neutrons */
    if (mt == 16) return ALEA_NUC_RXN_MULTIPLY;  /* (n,2n) */
    if (mt == 17) return ALEA_NUC_RXN_MULTIPLY;  /* (n,3n) */
    if (mt == 24) return ALEA_NUC_RXN_MULTIPLY;  /* (n,2nα) */
    if (mt == 25) return ALEA_NUC_RXN_MULTIPLY;  /* (n,3nα) */
    if (mt == 37) return ALEA_NUC_RXN_MULTIPLY;  /* (n,4n) */

    /* (n,nX) reactions — one neutron out plus charged particle */
    if (mt == 22) return ALEA_NUC_RXN_SCATTER;   /* (n,nα) */
    if (mt == 23) return ALEA_NUC_RXN_SCATTER;   /* (n,n3α) */
    if (mt == 28) return ALEA_NUC_RXN_SCATTER;   /* (n,np) */
    if (mt == 29) return ALEA_NUC_RXN_SCATTER;   /* (n,n2α) */
    if (mt == 30) return ALEA_NUC_RXN_MULTIPLY;   /* (n,2n2α) — 2 neutrons */
    if (mt == 32) return ALEA_NUC_RXN_SCATTER;   /* (n,nd) */
    if (mt == 33) return ALEA_NUC_RXN_SCATTER;   /* (n,nt) */
    if (mt == 34) return ALEA_NUC_RXN_SCATTER;   /* (n,nHe3) */
    if (mt == 35) return ALEA_NUC_RXN_SCATTER;   /* (n,nd2α) */
    if (mt == 36) return ALEA_NUC_RXN_SCATTER;   /* (n,nt2α) */

    /* Capture and charged-particle production — absorption */
    if (mt == 102) return ALEA_NUC_RXN_ABSORPTION; /* (n,γ) */
    if (mt >= 103 && mt <= 107) return ALEA_NUC_RXN_ABSORPTION; /* (n,p)/(n,d)/(n,t)/(n,He3)/(n,α) */
    if (mt == 108) return ALEA_NUC_RXN_ABSORPTION; /* (n,2α) */
    if (mt == 111) return ALEA_NUC_RXN_ABSORPTION; /* (n,2p) */
    if (mt == 112) return ALEA_NUC_RXN_ABSORPTION; /* (n,pα) */

    /* Default: absorption (safe for neutron balance) */
    return ALEA_NUC_RXN_ABSORPTION;
}

/**
 * Evaluate tabulated yield from XSS at the given locator.
 * Format: NR, NBT[NR], INT[NR], NE, E[NE], Y[NE]
 */
static double eval_tabulated_yield(const alea_nuc_nuclide_t* nuc, int loc, double energy) {
    const alea_nuc_ace_table_t* t = &nuc->raw;
    if (loc <= 0 || loc > t->xss_length) return 1.0;

    int nr = (int)t->xss[loc - 1];
    /* Bound nr against the remaining space *before* forming loc + 2*nr, so the
     * multiply cannot overflow int and wrap past the check. */
    if (nr < 0 || nr > (t->xss_length - loc) / 2) return 1.0;
    int base = loc + 2 * nr; /* skip NBT/INT pairs */
    /* base is the index of NE, read below — it must be a valid index, not one
     * past the end (the old `> xss_length` allowed base == xss_length). */
    if (base >= t->xss_length) return 1.0;

    int ne = (int)t->xss[base];
    /* Likewise bound ne before forming base + 1 + 2*ne. */
    if (ne <= 0 || ne > (t->xss_length - base - 1) / 2) return 1.0;

    const double* egrid = &t->xss[base];     /* NE, E[0..NE-1] */
    const double* evals = &t->xss[base + ne]; /* Y[0..NE-1] */

    double f;
    int ie = alea_nuc_energy_lookup(egrid + 1, ne, energy, &f);
    if (ie < 0) return evals[1]; /* below grid, use first value */
    return evals[ie + 1] + f * (evals[ie + 2] - evals[ie + 1]);
}

double alea_nuc_reaction_yield(const alea_nuc_nuclide_t* nuc, int mt, double energy) {
    if (!nuc) return 0.0;

    /* Find the reaction via O(1) MT table */
    const alea_nuc_reaction_t* r = NULL;
    if (nuc->mt_to_rxn && mt >= 0 && mt < ALEA_NUC_MT_TABLE_SIZE) {
        int idx = nuc->mt_to_rxn[mt];
        if (idx >= 0) r = &nuc->reactions[idx];
    } else {
        for (int i = 0; i < nuc->n_reactions; i++) {
            if (nuc->reactions[i].mt == mt) { r = &nuc->reactions[i]; break; }
        }
    }

    if (r) {

        int ty = r->ty;
        if (ty == 0) return 0.0; /* absorption */

        int abs_ty = abs(ty);

        /* TYR=19: fission, use ν̄ */
        if (abs_ty == 19) return alea_nuc_nu_bar(nuc, energy);

        /* |TYR| = 1-4: fixed integer yield */
        if (abs_ty >= 1 && abs_ty <= 4) return (double)abs_ty;

        /* |TYR| > 4: locator to tabulated yield in DLW block */
        if (abs_ty > 4) {
            int dlw_base = nuc->raw.jxs[10]; /* JXS[11] */
            if (dlw_base > 0)
                return eval_tabulated_yield(nuc, dlw_base + abs_ty - 1, energy);
        }

        return 1.0; /* fallback */
    }

    /* Reaction not in table — use MT classification */
    if (mt == 2) return 1.0;
    if (mt >= 51 && mt <= 91) return 1.0;

    alea_nuc_reaction_class_t cls = alea_nuc_reaction_classify(mt);
    if (cls == ALEA_NUC_RXN_ABSORPTION) return 0.0;
    if (cls == ALEA_NUC_RXN_SCATTER) return 1.0;

    /* Fission without nubar data */
    if (mt == 18 || mt == 19 || mt == 20 || mt == 21 || mt == 38)
        return alea_nuc_nu_bar(nuc, energy);

    /* (n,2n), (n,3n), (n,4n) */
    if (mt == 16) return 2.0;
    if (mt == 17) return 3.0;
    if (mt == 37) return 4.0;

    return 1.0;
}

/**
 * Evaluate ν̄ at given energy.
 */
double alea_nuc_nu_bar(const alea_nuc_nuclide_t* nuc, double energy) {
    if (!nuc || !nuc->fission) return 0.0;

    const alea_nuc_nu_bar_t* nu = nuc->fission->total;
    if (!nu) nu = nuc->fission->prompt;
    if (!nu) return 0.0;

    if (nu->type == ALEA_NUC_NU_POLYNOMIAL) {
        double result = 0.0;
        double En = 1.0;
        for (int i = 0; i < nu->n_coeffs; i++) {
            result += nu->coeffs[i] * En;
            En *= energy;
        }
        return result;
    }

    if (nu->type == ALEA_NUC_NU_TABULAR && nu->n_energies > 0) {
        double f;
        int ie = alea_nuc_energy_lookup(nu->energy, nu->n_energies, energy, &f);
        if (ie < 0) return nu->nu[0];
        return nu->nu[ie] + f * (nu->nu[ie + 1] - nu->nu[ie]);
    }

    return 0.0;
}
