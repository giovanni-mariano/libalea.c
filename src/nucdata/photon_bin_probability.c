// SPDX-FileCopyrightText: 2026 F4E
// SPDX-License-Identifier: EUPL-1.2

/** @file photon_bin_probability.c
 *  @brief Exact photon-energy bin probabilities for supported ACE laws.
 */
#include "nuclear_internal.h"
#include <math.h>

typedef struct {
    const alea_nuc_energy_dist_t* law;
    int lower;
    int upper;
    double fraction;
    double mapped_min;
    double mapped_max;
    int map_continuum;
} photon_tabular_context_t;

/* Reproduce the tabular sampler's incident table selection and unit-base
 * mapping. The production channel was validated during preparation. */
static int photon_tabular_prepare(const alea_nuc_energy_dist_t* law,
                                  double energy, photon_tabular_context_t* ctx) {
    if (!law || law->next ||
        (law->law != ALEA_NUC_ELAW_CONT_TABULAR &&
         law->law != ALEA_NUC_ELAW_KALBACH &&
         law->law != ALEA_NUC_ELAW_CORRELATED)) return 0;
    int n = law->tab.n_ein;
    int lower = 0;
    double f = 0.0;
    if (n > 1) {
        lower = alea_nuc_energy_lookup(law->tab.ein, n, energy, &f);
        if (lower < 0) return 0;
        int code = 2;
        for (int i = 0; i < law->tab.n_regions; i++) {
            if (lower + 2 <= law->tab.nbt[i]) {
                code = law->tab.interp[i];
                break;
            }
        }
        if (code >= 10) code %= 10;
        if (code == 1) f = 0.0;
        else if (code != 2) return 0;
    }
    int upper = n == 1 ? lower : lower + 1;
    int dl = law->tab.n_discrete[lower];
    int du = law->tab.n_discrete[upper];
    ctx->law = law;
    ctx->lower = lower;
    ctx->upper = upper;
    ctx->fraction = f;
    ctx->mapped_min = 0.0;
    ctx->mapped_max = 0.0;
    ctx->map_continuum = 0;
    if (lower != upper && dl < law->tab.n_eout[lower] &&
        du < law->tab.n_eout[upper]) {
        ctx->mapped_min = law->tab.eout[lower][dl] + f *
            (law->tab.eout[upper][du] - law->tab.eout[lower][dl]);
        ctx->mapped_max =
            law->tab.eout[lower][law->tab.n_eout[lower] - 1] + f *
            (law->tab.eout[upper][law->tab.n_eout[upper] - 1] -
             law->tab.eout[lower][law->tab.n_eout[lower] - 1]);
        if (!(ctx->mapped_max > ctx->mapped_min)) return 0;
        ctx->map_continuum = 1;
    }
    if (lower != upper && dl != du) return 0;
    return 1;
}

static double photon_tabular_continuous_cdf(
        const alea_nuc_energy_dist_t* law, int table, double x) {
    int d = law->tab.n_discrete[table];
    int n = law->tab.n_eout[table];
    const double* energy = law->tab.eout[table];
    const double* pdf = law->tab.pdf[table];
    const double* cdf = law->tab.cdf[table];
    if (x <= energy[d]) return cdf[d];
    if (x >= energy[n - 1]) return cdf[n - 1];
    int lo = d, hi = n - 1;
    while (lo + 1 < hi) {
        int mid = lo + (hi - lo) / 2;
        if (energy[mid] <= x) lo = mid;
        else hi = mid;
    }
    double dx = x - energy[lo];
    double area = pdf[lo] * dx;
    if (law->tab.interpolation[table] == 2) {
        double width = energy[lo + 1] - energy[lo];
        if (width > 0.0)
            area += 0.5 * (pdf[lo + 1] - pdf[lo]) * dx * dx / width;
    }
    return fmin(cdf[lo + 1], fmax(cdf[lo], cdf[lo] + area));
}

static double photon_tabular_bin_probability(
        const photon_tabular_context_t* ctx, double low, double high,
        int include_high) {
    const alea_nuc_energy_dist_t* law = ctx->law;
    double probability = 0.0;
    for (int choice = 0; choice < (ctx->lower == ctx->upper ? 1 : 2);
         choice++) {
        int table = choice == 0 ? ctx->lower : ctx->upper;
        double weight = ctx->lower == ctx->upper ? 1.0 :
            (choice == 0 ? 1.0 - ctx->fraction : ctx->fraction);
        if (weight <= 0.0) continue;
        int d = law->tab.n_discrete[table];
        int n = law->tab.n_eout[table];
        const double* out = law->tab.eout[table];
        const double* cdf = law->tab.cdf[table];
        double mass = 0.0;
        for (int i = 0; i < d; i++) {
            double line = out[i];
            if (ctx->lower != ctx->upper)
                line = law->tab.eout[ctx->lower][i] + ctx->fraction *
                    (law->tab.eout[ctx->upper][i] -
                     law->tab.eout[ctx->lower][i]);
            if (line >= low && (line < high ||
                                (include_high && line == high)))
                mass += cdf[i] - (i ? cdf[i - 1] : 0.0);
        }
        if (d < n) {
            double raw_low = low;
            double raw_high = high;
            if (ctx->map_continuum) {
                double raw_min = out[d];
                double scale = (out[n - 1] - raw_min) /
                    (ctx->mapped_max - ctx->mapped_min);
                raw_low = raw_min + (low - ctx->mapped_min) * scale;
                raw_high = raw_min + (high - ctx->mapped_min) * scale;
            }
            mass += photon_tabular_continuous_cdf(law, table, raw_high) -
                photon_tabular_continuous_cdf(law, table, raw_low);
        }
        probability += weight * fmax(0.0, mass);
    }
    return probability;
}

alea_error_t alea_nuc_prepared_photon_bin_probabilities(
        const alea_nuc_prepared_photon_production_t* prepared,
        double incident_energy, const double* edges, size_t n_edges,
        double* probabilities) {
    if (!prepared || !prepared->nuclide || !prepared->production ||
        !prepared->production->spectrum || !edges || !probabilities ||
        !isfinite(incident_energy) || incident_energy < 0.0 || n_edges < 2)
        return ALEA_ERR_INVALID_ARG;
    for (size_t i = 0; i < n_edges; i++)
        if (!isfinite(edges[i]) || (i && !(edges[i] > edges[i - 1])))
            return ALEA_ERR_INVALID_ARG;

    const alea_nuc_energy_dist_t* law = prepared->production->spectrum;
    if (law->next) return ALEA_ERR_UNSUPPORTED;
    if (law->law == ALEA_NUC_ELAW_DISCRETE_PHOTON ||
        law->law == ALEA_NUC_ELAW_LEVEL) {
        double line;
        if (law->law == ALEA_NUC_ELAW_DISCRETE_PHOTON) {
            line = law->discrete_photon_energy;
            if (law->discrete_photon_primary == 2)
                line += law->discrete_photon_awr /
                    (law->discrete_photon_awr + 1.0) * incident_energy;
        } else {
            line = law->level_Q * (incident_energy - law->level_A);
        }
        if (!isfinite(line) || line < 0.0) return ALEA_ERR_INVALID_STATE;
        for (size_t i = 0; i + 1 < n_edges; i++)
            probabilities[i] = line >= edges[i] &&
                (line < edges[i + 1] ||
                 (i + 2 == n_edges && line == edges[i + 1])) ? 1.0 : 0.0;
        return ALEA_OK;
    }

    photon_tabular_context_t context;
    if (!photon_tabular_prepare(law, incident_energy, &context))
        return ALEA_ERR_UNSUPPORTED;
    for (size_t i = 0; i + 1 < n_edges; i++)
        probabilities[i] = photon_tabular_bin_probability(
            &context, edges[i], edges[i + 1], i + 2 == n_edges);
    return ALEA_OK;
}
