// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/**
 * @file energy_sample.c
 * @brief Sampling of ACE secondary-energy distributions
 */

#include "nuclear_internal.h"

#include <float.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define ALEA_NUC_MAX_REJECTION_ATTEMPTS 100000
#define ALEA_NUC_KALBACH_R_TOL 1e-3

static alea_error_t draw_uniform(alea_nuc_random_fn random, void* context,
                                 double* value) {
    if (!random || !value) return ALEA_ERR_NULL_ARG;
    double u = random(context);
    if (!isfinite(u) || u < 0.0 || u >= 1.0) return ALEA_ERR_INVALID_STATE;
    *value = u;
    return ALEA_OK;
}

static int grid_valid(const double* x, int n) {
    if (!x || n <= 0) return 0;
    for (int i = 0; i < n; i++)
        if (!isfinite(x[i]) || (i > 0 && x[i] < x[i - 1])) return 0;
    return 1;
}

static int table_valid(const double* x, const double* y, int n,
                       const int* nbt, const int* interp, int n_regions) {
    if (!grid_valid(x, n) || !y ||
        !alea_nuc_interp_regions_valid(nbt, interp, n_regions, n)) return 0;
    for (int i = 0; i < n; i++)
        if (!isfinite(y[i])) return 0;
    return 1;
}

/* ACE permits 11/12 and 21/22 for histogram/linear incident-energy
 * interpolation with unit-base mapping of the outgoing distribution. */
static int incident_regions_valid(const int* nbt, const int* interp,
                                  int n_regions, int n_points) {
    if (n_points <= 0 || n_regions < 0) return 0;
    if (n_regions == 0) return 1;
    if (!nbt || !interp || n_points < 2) return 0;
    for (int i = 0; i < n_regions; i++) {
        int code = interp[i];
        if (nbt[i] < 2 || nbt[i] > n_points ||
            (i > 0 && nbt[i] <= nbt[i - 1]) ||
            (code != 1 && code != 2 && code != 11 && code != 12 &&
             code != 21 && code != 22))
            return 0;
    }
    return nbt[n_regions - 1] == n_points;
}

static int angular_point_valid(const alea_nuc_angular_point_t* point) {
    if (!point) return 0;
    if (point->type == ALEA_NUC_ANG_ISOTROPIC) return 1;
    if (point->type != ALEA_NUC_ANG_TABULAR || point->n_cosines < 2)
        return 0;
    if (!alea_nuc_tabular_pdf_valid(point->cosine, point->pdf, point->cdf,
                                    point->n_cosines, point->interpolation, 0))
        return 0;
    for (int i = 0; i < point->n_cosines; i++)
        if (point->cosine[i] < -1.0 || point->cosine[i] > 1.0) return 0;
    return 1;
}

static alea_error_t table_value(const double* x, const double* y, int n,
                                const int* nbt, const int* interp,
                                int n_regions, double query, double* value) {
    if (!table_valid(x, y, n, nbt, interp, n_regions))
        return ALEA_ERR_INVALID_ARG;
    return alea_nuc_interp_eval(x, y, n, nbt, interp, n_regions,
                                query, value);
}

alea_error_t alea_nuc_energy_dist_validate(
    const alea_nuc_energy_dist_t* distribution, int* unsupported_law) {
    if (unsupported_law) *unsupported_law = 0;
    if (!distribution) return ALEA_ERR_INVALID_ARG;
    int count = 0;
    for (const alea_nuc_energy_dist_t* law = distribution; law;
         law = law->next) {
        if (++count > 10) return ALEA_ERR_INVALID_ARG;
        if (law->next &&
            !table_valid(law->energy, law->probability, law->n_energies,
                         law->nbt, law->interp, law->n_regions))
            return ALEA_ERR_INVALID_ARG;
        if (law->next) {
            for (int i = 0; i < law->n_energies; i++)
                if (law->probability[i] < 0.0 || law->probability[i] > 1.0)
                    return ALEA_ERR_INVALID_ARG;
        }
        if (law->law == ALEA_NUC_ELAW_DISCRETE_PHOTON) {
            if ((law->discrete_photon_primary != 0 &&
                 law->discrete_photon_primary != 1 &&
                 law->discrete_photon_primary != 2) ||
                !isfinite(law->discrete_photon_energy) ||
                law->discrete_photon_energy < 0.0 ||
                !isfinite(law->discrete_photon_awr) ||
                law->discrete_photon_awr <= 0.0)
                return ALEA_ERR_INVALID_ARG;
            continue;
        }
        if (law->law == ALEA_NUC_ELAW_LEVEL) {
            if (!isfinite(law->level_A) || !isfinite(law->level_Q) ||
                law->level_Q < 0.0) return ALEA_ERR_INVALID_ARG;
            continue;
        }
        if (law->law == ALEA_NUC_ELAW_CONT_TABULAR ||
            law->law == ALEA_NUC_ELAW_KALBACH ||
            law->law == ALEA_NUC_ELAW_CORRELATED) {
            if (!grid_valid(law->tab.ein, law->tab.n_ein) ||
                !incident_regions_valid(
                    law->tab.nbt, law->tab.interp, law->tab.n_regions,
                    law->tab.n_ein) || !law->tab.interpolation ||
                !law->tab.n_discrete || !law->tab.n_eout ||
                !law->tab.eout || !law->tab.pdf || !law->tab.cdf)
                return ALEA_ERR_INVALID_ARG;
            for (int i = 0; i + 1 < law->tab.n_ein; i++) {
                int code = alea_nuc_interp_code_for_interval(
                    law->tab.nbt, law->tab.interp, law->tab.n_regions, i);
                if (code >= 10) code %= 10;
                if (code != 1 && code != 2) return ALEA_ERR_UNSUPPORTED;
            }
            for (int i = 0; i < law->tab.n_ein; i++)
                if (!alea_nuc_tabular_pdf_valid(
                        law->tab.eout[i], law->tab.pdf[i], law->tab.cdf[i],
                        law->tab.n_eout[i], law->tab.interpolation[i],
                        law->tab.n_discrete[i])) return ALEA_ERR_INVALID_ARG;
            if (law->law == ALEA_NUC_ELAW_KALBACH) {
                if (!law->tab.precompound_r || !law->tab.precompound_a)
                    return ALEA_ERR_INVALID_ARG;
                for (int i = 0; i < law->tab.n_ein; i++) {
                    if (!law->tab.precompound_r[i] ||
                        !law->tab.precompound_a[i]) return ALEA_ERR_INVALID_ARG;
                    for (int j = 0; j < law->tab.n_eout[i]; j++) {
                        double r = law->tab.precompound_r[i][j];
                        double a = law->tab.precompound_a[i][j];
                        if (!isfinite(r) || r < -ALEA_NUC_KALBACH_R_TOL ||
                            r > 1.0 + ALEA_NUC_KALBACH_R_TOL ||
                            !isfinite(a) || a < 0.0) return ALEA_ERR_INVALID_ARG;
                    }
                }
            }
            if (law->law == ALEA_NUC_ELAW_CORRELATED) {
                if (!law->tab.correlated_mu) return ALEA_ERR_INVALID_ARG;
                for (int i = 0; i < law->tab.n_ein; i++) {
                    if (!law->tab.correlated_mu[i]) return ALEA_ERR_INVALID_ARG;
                    for (int j = 0; j < law->tab.n_eout[i]; j++)
                        if (!angular_point_valid(
                                &law->tab.correlated_mu[i][j]))
                            return ALEA_ERR_INVALID_ARG;
                }
            }
            continue;
        }
        if (law->law == ALEA_NUC_ELAW_GENERAL_EVAP) {
            if (!table_valid(law->temp_energy, law->temp_T, law->n_temp,
                             law->temp_nbt, law->temp_interp,
                             law->n_temp_regions) ||
                law->n_general_evap < 2 || !law->general_evap_x)
                return ALEA_ERR_INVALID_ARG;
            for (int i = 0; i < law->n_temp; i++)
                if (law->temp_T[i] <= 0.0) return ALEA_ERR_INVALID_ARG;
            for (int i = 0; i < law->n_general_evap; i++)
                if (!isfinite(law->general_evap_x[i]) ||
                    law->general_evap_x[i] < 0.0 ||
                    (i > 0 && law->general_evap_x[i] <
                                      law->general_evap_x[i - 1]))
                    return ALEA_ERR_INVALID_ARG;
            continue;
        }
        if (law->law == ALEA_NUC_ELAW_MAXWELL ||
            law->law == ALEA_NUC_ELAW_EVAPORATION ||
            law->law == ALEA_NUC_ELAW_WATT) {
            if (!table_valid(law->temp_energy, law->temp_T, law->n_temp,
                             law->temp_nbt, law->temp_interp,
                             law->n_temp_regions) ||
                !isfinite(law->restriction_energy))
                return ALEA_ERR_INVALID_ARG;
            for (int i = 0; i < law->n_temp; i++)
                if (law->temp_T[i] <= 0.0) return ALEA_ERR_INVALID_ARG;
            if (law->law == ALEA_NUC_ELAW_WATT) {
                if (!table_valid(law->watt_b_energy, law->temp_C,
                                 law->n_watt_b, law->watt_b_nbt,
                                 law->watt_b_interp,
                                 law->n_watt_b_regions))
                    return ALEA_ERR_INVALID_ARG;
                for (int i = 0; i < law->n_watt_b; i++)
                    if (law->temp_C[i] < 0.0) return ALEA_ERR_INVALID_ARG;
            }
            continue;
        }
        if (law->law == ALEA_NUC_ELAW_NBODY) {
            if (law->nbody_particles < 3 || law->nbody_particles > 5 ||
                !isfinite(law->nbody_total_mass) ||
                law->nbody_total_mass <= 1.0 ||
                !isfinite(law->nbody_target_awr) ||
                law->nbody_target_awr <= 0.0 ||
                !isfinite(law->nbody_q_value)) return ALEA_ERR_INVALID_ARG;
            continue;
        }
        if (law->law == ALEA_NUC_ELAW_LAB_ANGLE_ENERGY) {
            if (!grid_valid(law->law67.ein, law->law67.n_ein) ||
                !incident_regions_valid(
                    law->law67.nbt, law->law67.interp,
                    law->law67.n_regions, law->law67.n_ein) ||
                !law->law67.incident)
                return ALEA_ERR_INVALID_ARG;
            for (int i = 0; i < law->law67.n_ein; i++) {
                const alea_nuc_law67_incident_t* incident =
                    &law->law67.incident[i];
                if ((incident->interpolation != 1 &&
                     incident->interpolation != 2) ||
                    !grid_valid(incident->cosine, incident->n_cosines) ||
                    !incident->spectrum || incident->cosine[0] < -1.0 ||
                    incident->cosine[incident->n_cosines - 1] > 1.0)
                    return ALEA_ERR_INVALID_ARG;
                for (int j = 0; j < incident->n_cosines; j++) {
                    const alea_nuc_law67_energy_t* spectrum =
                        &incident->spectrum[j];
                    if (!alea_nuc_tabular_pdf_valid(
                            spectrum->energy, spectrum->pdf, spectrum->cdf,
                            spectrum->n_points, spectrum->interpolation, 0))
                        return ALEA_ERR_INVALID_ARG;
                }
            }
            continue;
        }
        if (unsupported_law) *unsupported_law = (int)law->law;
        return ALEA_ERR_UNSUPPORTED;
    }
    return ALEA_OK;
}

static alea_error_t choose_law(const alea_nuc_energy_dist_t* head,
                               double incident_energy,
                               alea_nuc_random_fn random, void* context,
                               const alea_nuc_energy_dist_t** selected) {
    const alea_nuc_energy_dist_t* law = head;
    if (!law || !selected) return ALEA_ERR_NULL_ARG;

    while (law->next) {
        double probability;
        alea_error_t err = table_value(law->energy, law->probability,
                                       law->n_energies, law->nbt, law->interp,
                                       law->n_regions, incident_energy,
                                       &probability);
        if (err != ALEA_OK || probability < 0.0 || probability > 1.0)
            return ALEA_ERR_INVALID_ARG;
        double u;
        err = draw_uniform(random, context, &u);
        if (err != ALEA_OK) return err;
        if (u < probability) break;
        law = law->next;
    }
    *selected = law;
    return ALEA_OK;
}

static alea_error_t incident_interval(const alea_nuc_energy_dist_t* law,
                                      double incident_energy, int* lower,
                                      double* fraction) {
    const int n = law->tab.n_ein;
    if (!grid_valid(law->tab.ein, n) || !lower || !fraction ||
        !incident_regions_valid(law->tab.nbt, law->tab.interp,
                                law->tab.n_regions, n))
        return ALEA_ERR_INVALID_ARG;
    if (n == 1) {
        *lower = 0;
        *fraction = 0.0;
        return ALEA_OK;
    }
    int i = alea_nuc_energy_lookup_trusted(law->tab.ein, n,
                                           incident_energy, fraction);
    if (i < 0) return ALEA_ERR_INVALID_ARG;
    int code = alea_nuc_interp_code_for_interval(
        law->tab.nbt, law->tab.interp, law->tab.n_regions, i);
    if (code >= 10) code %= 10;
    if (code == 1) *fraction = 0.0;
    else if (code != 2) return ALEA_ERR_UNSUPPORTED;
    *lower = i;
    return ALEA_OK;
}

static alea_error_t sample_continuous_tabular_details(
    const alea_nuc_energy_dist_t* law, double incident_energy,
    alea_nuc_random_fn random, void* context, double* energy,
    int* selected_table, int* selected_index, double* raw_energy) {
    if (!law->tab.interpolation || !law->tab.n_discrete ||
        !law->tab.n_eout || !law->tab.eout || !law->tab.pdf || !law->tab.cdf)
        return ALEA_ERR_INVALID_ARG;

    int lower;
    double f;
    alea_error_t err = incident_interval(law, incident_energy, &lower, &f);
    if (err != ALEA_OK) return err;
    int upper = law->tab.n_ein == 1 ? lower : lower + 1;
    int selected = lower;
    if (upper != lower) {
        double u;
        err = draw_uniform(random, context, &u);
        if (err != ALEA_OK) return err;
        if (u < f) selected = upper;
    }

    int n = law->tab.n_eout[selected];
    if (!alea_nuc_tabular_pdf_valid(law->tab.eout[selected],
                                    law->tab.pdf[selected],
                                    law->tab.cdf[selected], n,
                                    law->tab.interpolation[selected],
                                    law->tab.n_discrete[selected]))
        return ALEA_ERR_INVALID_ARG;
    double u, sampled;
    int sampled_index = -1;
    err = draw_uniform(random, context, &u);
    if (err != ALEA_OK) return err;
    err = alea_nuc_tabular_pdf_sample(
        law->tab.eout[selected], law->tab.pdf[selected],
        law->tab.cdf[selected], n, law->tab.interpolation[selected],
        law->tab.n_discrete[selected], u, &sampled, &sampled_index);
    if (err != ALEA_OK) return err;
    if (selected_table) *selected_table = selected;
    if (selected_index) *selected_index = sampled_index;
    if (raw_energy) *raw_energy = sampled;
    if (upper == lower) {
        *energy = sampled;
        return ALEA_OK;
    }

    int selected_discrete = law->tab.n_discrete[selected];
    if (sampled_index >= 0 && sampled_index < selected_discrete) {
        int lower_discrete = law->tab.n_discrete[lower];
        int upper_discrete = law->tab.n_discrete[upper];
        if (sampled_index >= lower_discrete ||
            sampled_index >= upper_discrete)
            return ALEA_ERR_UNSUPPORTED;
        *energy = law->tab.eout[lower][sampled_index] + f *
            (law->tab.eout[upper][sampled_index] -
             law->tab.eout[lower][sampled_index]);
        return ALEA_OK;
    }

    if (law->tab.n_discrete[lower] >= law->tab.n_eout[lower] ||
        law->tab.n_discrete[upper] >= law->tab.n_eout[upper])
        return ALEA_ERR_UNSUPPORTED;
    const double selected_min =
        law->tab.eout[selected][selected_discrete];
    const double selected_max = law->tab.eout[selected][n - 1];
    const double interpolated_min =
        law->tab.eout[lower][law->tab.n_discrete[lower]] + f *
        (law->tab.eout[upper][law->tab.n_discrete[upper]] -
         law->tab.eout[lower][law->tab.n_discrete[lower]]);
    const double interpolated_max =
        law->tab.eout[lower][law->tab.n_eout[lower] - 1] +
        f * (law->tab.eout[upper][law->tab.n_eout[upper] - 1] -
             law->tab.eout[lower][law->tab.n_eout[lower] - 1]);
    if (selected_max > selected_min) {
        *energy = interpolated_min + (sampled - selected_min) *
                  (interpolated_max - interpolated_min) /
                  (selected_max - selected_min);
        return ALEA_OK;
    }

    /* A pure, single discrete line has no unit-base width. */
    if (n == 1 && law->tab.n_eout[lower] == 1 &&
        law->tab.n_eout[upper] == 1 && sampled_index == 0) {
        *energy = law->tab.eout[lower][0] +
                  f * (law->tab.eout[upper][0] - law->tab.eout[lower][0]);
        return ALEA_OK;
    }
    return ALEA_ERR_UNSUPPORTED;
}

static alea_error_t sample_continuous_tabular(
    const alea_nuc_energy_dist_t* law, double incident_energy,
    alea_nuc_random_fn random, void* context, double* energy) {
    return sample_continuous_tabular_details(law, incident_energy, random,
                                             context, energy, NULL, NULL, NULL);
}

static alea_error_t sample_maxwell(double temperature, double limit,
                                   alea_nuc_random_fn random, void* context,
                                   double* energy) {
    if (!(temperature > 0.0) || !(limit >= 0.0) || isnan(limit))
        return ALEA_ERR_INVALID_ARG;
    for (int attempt = 0; attempt < ALEA_NUC_MAX_REJECTION_ATTEMPTS; attempt++) {
        double u1, u2, u3;
        alea_error_t err = draw_uniform(random, context, &u1);
        if (err != ALEA_OK) return err;
        if ((err = draw_uniform(random, context, &u2)) != ALEA_OK) return err;
        if ((err = draw_uniform(random, context, &u3)) != ALEA_OK) return err;
        if (u1 == 0.0) u1 = DBL_MIN;
        if (u2 == 0.0) u2 = DBL_MIN;
        double c = cos(0.5 * M_PI * u3);
        double sampled = -temperature * (log(u1) + log(u2) * c * c);
        if (sampled <= limit) {
            *energy = sampled;
            return ALEA_OK;
        }
    }
    return ALEA_ERR_INVALID_STATE;
}

static alea_error_t sample_evaporation(double temperature, double limit,
                                       alea_nuc_random_fn random,
                                       void* context, double* energy) {
    if (!(temperature > 0.0) || !(limit >= 0.0) || !isfinite(limit))
        return ALEA_ERR_INVALID_ARG;
    double g = -expm1(-limit / temperature);
    for (int attempt = 0; attempt < ALEA_NUC_MAX_REJECTION_ATTEMPTS; attempt++) {
        double u1, u2;
        alea_error_t err = draw_uniform(random, context, &u1);
        if (err != ALEA_OK) return err;
        if ((err = draw_uniform(random, context, &u2)) != ALEA_OK) return err;
        double sampled = -temperature * log((1.0 - g * u1) *
                                            (1.0 - g * u2));
        if (sampled <= limit) {
            *energy = sampled;
            return ALEA_OK;
        }
    }
    return ALEA_ERR_INVALID_STATE;
}

static alea_error_t negative_log_product(int count,
                                         alea_nuc_random_fn random,
                                         void* context, double* value) {
    double sum = 0.0;
    for (int i = 0; i < count; i++) {
        double u;
        alea_error_t err = draw_uniform(random, context, &u);
        if (err != ALEA_OK) return err;
        if (u == 0.0) u = DBL_MIN;
        sum -= log(u);
    }
    *value = sum;
    return ALEA_OK;
}

static alea_error_t sample_nbody(const alea_nuc_energy_dist_t* law,
                                 double incident_energy,
                                 alea_nuc_random_fn random, void* context,
                                 double* energy) {
    double available = law->nbody_target_awr /
        (law->nbody_target_awr + 1.0) * incident_energy +
        law->nbody_q_value;
    double maximum = (law->nbody_total_mass - 1.0) /
                     law->nbody_total_mass * available;
    if (!(maximum >= 0.0) || !isfinite(maximum)) return ALEA_ERR_INVALID_ARG;
    double x, y;
    alea_error_t err = sample_maxwell(1.0, INFINITY, random, context, &x);
    if (err != ALEA_OK) return err;
    if (law->nbody_particles == 3) {
        err = sample_maxwell(1.0, INFINITY, random, context, &y);
    } else if (law->nbody_particles == 4) {
        err = negative_log_product(3, random, context, &y);
    } else if (law->nbody_particles == 5) {
        err = negative_log_product(4, random, context, &y);
        if (err == ALEA_OK) {
            double u1, u2;
            if ((err = draw_uniform(random, context, &u1)) != ALEA_OK)
                return err;
            if ((err = draw_uniform(random, context, &u2)) != ALEA_OK)
                return err;
            if (u1 == 0.0) u1 = DBL_MIN;
            double c = cos(0.5 * M_PI * u2);
            y -= log(u1) * c * c;
        }
    } else {
        return ALEA_ERR_UNSUPPORTED;
    }
    if (err != ALEA_OK) return err;
    if (x + y <= 0.0) return ALEA_ERR_INVALID_STATE;
    *energy = maximum * x / (x + y);
    return ALEA_OK;
}

static int law67_cosine_interval(const alea_nuc_law67_incident_t* incident,
                                 double mu, double* fraction) {
    int n = incident->n_cosines;
    if (mu <= incident->cosine[0]) {
        *fraction = 0.0;
        return 0;
    }
    if (mu >= incident->cosine[n - 1]) {
        *fraction = 1.0;
        return n - 2;
    }
    int lo = 0, hi = n - 1;
    while (hi - lo > 1) {
        int mid = lo + (hi - lo) / 2;
        if (incident->cosine[mid] <= mu) lo = mid;
        else hi = mid;
    }
    if (incident->interpolation == 1) *fraction = 0.0;
    else *fraction = (mu - incident->cosine[lo]) /
        (incident->cosine[lo + 1] - incident->cosine[lo]);
    return lo;
}

static alea_error_t law67_choose_spectrum(
    const alea_nuc_law67_incident_t* incident, double mu,
    alea_nuc_random_fn random, void* context, int* selected) {
    double fraction, u;
    int lower = law67_cosine_interval(incident, mu, &fraction);
    alea_error_t err = draw_uniform(random, context, &u);
    if (err != ALEA_OK) return err;
    *selected = u < fraction ? lower + 1 : lower;
    return ALEA_OK;
}

static alea_error_t sample_law67(const alea_nuc_energy_dist_t* law,
                                 double incident_energy,
                                 alea_nuc_random_fn random, void* context,
                                 double* energy, double* mu) {
    int lower;
    double f;
    const int n = law->law67.n_ein;
    if (!grid_valid(law->law67.ein, n) || n < 1)
        return ALEA_ERR_INVALID_ARG;
    if (n == 1) {
        lower = 0;
        f = 0.0;
    } else {
        lower = alea_nuc_energy_lookup_trusted(
            law->law67.ein, n, incident_energy, &f);
        if (lower < 0) return ALEA_ERR_INVALID_ARG;
        int code = alea_nuc_interp_code_for_interval(
            law->law67.nbt, law->law67.interp,
            law->law67.n_regions, lower);
        if (code >= 10) code -= 10;
        if (code == 1) f = 0.0;
        else if (code != 2) return ALEA_ERR_UNSUPPORTED;
    }
    int upper = n == 1 ? lower : lower + 1;

    double u;
    alea_error_t err = draw_uniform(random, context, &u);
    if (err != ALEA_OK) return err;
    int angular_table = (upper != lower && u < f) ? upper : lower;
    const alea_nuc_law67_incident_t* angular =
        &law->law67.incident[angular_table];
    if (angular->n_cosines < 2) return ALEA_ERR_INVALID_ARG;
    if ((err = draw_uniform(random, context, &u)) != ALEA_OK) return err;
    double scaled = u * (angular->n_cosines - 1);
    int bin = (int)scaled;
    if (bin >= angular->n_cosines - 1) bin = angular->n_cosines - 2;
    double bin_fraction = scaled - bin;
    *mu = angular->cosine[bin] + bin_fraction *
        (angular->cosine[bin + 1] - angular->cosine[bin]);

    int lower_spectrum;
    err = law67_choose_spectrum(&law->law67.incident[lower], *mu,
                                random, context, &lower_spectrum);
    if (err != ALEA_OK) return err;
    int upper_spectrum = lower_spectrum;
    if (upper != lower) {
        err = law67_choose_spectrum(&law->law67.incident[upper], *mu,
                                    random, context, &upper_spectrum);
        if (err != ALEA_OK) return err;
    }
    if ((err = draw_uniform(random, context, &u)) != ALEA_OK) return err;
    int incident_table = upper != lower && u < f ? upper : lower;
    int spectrum_index = incident_table == lower
        ? lower_spectrum : upper_spectrum;
    const alea_nuc_law67_energy_t* spectrum =
        &law->law67.incident[incident_table].spectrum[spectrum_index];
    double raw;
    if ((err = draw_uniform(random, context, &u)) != ALEA_OK) return err;
    err = alea_nuc_tabular_pdf_sample(
        spectrum->energy, spectrum->pdf, spectrum->cdf,
        spectrum->n_points, spectrum->interpolation, 0, u, &raw, NULL);
    if (err != ALEA_OK) return err;

    double low_min = law->law67.incident[lower]
        .spectrum[lower_spectrum].energy[0];
    int low_count = law->law67.incident[lower]
        .spectrum[lower_spectrum].n_points;
    double low_max = law->law67.incident[lower]
        .spectrum[lower_spectrum].energy[low_count - 1];
    double high_min = low_min, high_max = low_max;
    if (upper != lower) {
        const alea_nuc_law67_energy_t* high =
            &law->law67.incident[upper].spectrum[upper_spectrum];
        high_min = high->energy[0];
        high_max = high->energy[high->n_points - 1];
    }
    double selected_min = spectrum->energy[0];
    double selected_max = spectrum->energy[spectrum->n_points - 1];
    double interpolated_min = low_min + f * (high_min - low_min);
    double interpolated_max = low_max + f * (high_max - low_max);
    if (!(selected_max > selected_min)) return ALEA_ERR_INVALID_ARG;
    *energy = interpolated_min + (raw - selected_min) *
        (interpolated_max - interpolated_min) /
        (selected_max - selected_min);
    return isfinite(*energy) && *energy >= 0.0
        ? ALEA_OK : ALEA_ERR_INVALID_STATE;
}

static alea_error_t sample_one(const alea_nuc_energy_dist_t* law,
                               double incident_energy,
                               alea_nuc_random_fn random, void* context,
                               double* energy) {
    if (!isfinite(incident_energy) || incident_energy < 0.0)
        return ALEA_ERR_INVALID_ARG;
    if (law->law == ALEA_NUC_ELAW_DISCRETE_PHOTON) {
        double sampled = law->discrete_photon_energy;
        if (law->discrete_photon_primary == 2)
            sampled += law->discrete_photon_awr /
                       (law->discrete_photon_awr + 1.0) * incident_energy;
        if (!isfinite(sampled) || sampled < 0.0)
            return ALEA_ERR_INVALID_ARG;
        *energy = sampled;
        return ALEA_OK;
    }
    if (law->law == ALEA_NUC_ELAW_LEVEL) {
        double sampled = law->level_Q * (incident_energy - law->level_A);
        if (!isfinite(law->level_A) || !isfinite(law->level_Q) ||
            law->level_Q < 0.0 || sampled < 0.0 || !isfinite(sampled))
            return ALEA_ERR_INVALID_ARG;
        *energy = sampled;
        return ALEA_OK;
    }
    if (law->law == ALEA_NUC_ELAW_CONT_TABULAR ||
        law->law == ALEA_NUC_ELAW_KALBACH ||
        law->law == ALEA_NUC_ELAW_CORRELATED)
        return sample_continuous_tabular(law, incident_energy, random,
                                         context, energy);
    if (law->law == ALEA_NUC_ELAW_NBODY)
        return sample_nbody(law, incident_energy, random, context, energy);
    if (law->law == ALEA_NUC_ELAW_LAB_ANGLE_ENERGY) {
        double mu;
        return sample_law67(law, incident_energy, random, context, energy,
                            &mu);
    }

    double a;
    alea_error_t err = table_value(law->temp_energy, law->temp_T,
                                   law->n_temp, law->temp_nbt,
                                   law->temp_interp, law->n_temp_regions,
                                   incident_energy, &a);
    if (err != ALEA_OK || !isfinite(law->restriction_energy))
        return ALEA_ERR_INVALID_ARG;
    if (law->law == ALEA_NUC_ELAW_GENERAL_EVAP) {
        double u;
        if ((err = draw_uniform(random, context, &u)) != ALEA_OK) return err;
        double scaled = u * (law->n_general_evap - 1);
        int bin = (int)scaled;
        if (bin >= law->n_general_evap - 1)
            bin = law->n_general_evap - 2;
        double fraction = scaled - bin;
        *energy = a * (law->general_evap_x[bin] +
                       fraction * (law->general_evap_x[bin + 1] -
                                   law->general_evap_x[bin]));
        return isfinite(*energy) ? ALEA_OK : ALEA_ERR_INVALID_STATE;
    }
    double limit = incident_energy - law->restriction_energy;
    if (!isfinite(limit) || limit < 0.0) return ALEA_ERR_INVALID_ARG;
    if (law->law == ALEA_NUC_ELAW_MAXWELL)
        return sample_maxwell(a, limit, random, context, energy);
    if (law->law == ALEA_NUC_ELAW_EVAPORATION)
        return sample_evaporation(a, limit, random, context, energy);
    if (law->law == ALEA_NUC_ELAW_WATT) {
        double b, w;
        err = table_value(law->watt_b_energy, law->temp_C, law->n_watt_b,
                          law->watt_b_nbt, law->watt_b_interp,
                          law->n_watt_b_regions, incident_energy, &b);
        if (err != ALEA_OK || b < 0.0) return ALEA_ERR_INVALID_ARG;
        for (int attempt = 0; attempt < ALEA_NUC_MAX_REJECTION_ATTEMPTS;
             attempt++) {
            err = sample_maxwell(a, INFINITY, random, context, &w);
            if (err != ALEA_OK) return err;
            double u;
            if ((err = draw_uniform(random, context, &u)) != ALEA_OK) return err;
            double sampled = w + 0.25 * a * a * b +
                             (2.0 * u - 1.0) * sqrt(a * a * b * w);
            if (sampled >= 0.0 && sampled <= limit) {
                *energy = sampled;
                return ALEA_OK;
            }
        }
        return ALEA_ERR_INVALID_STATE;
    }
    return ALEA_ERR_UNSUPPORTED;
}

static double interpolate_outgoing_parameter(
    const alea_nuc_energy_dist_t* law, int table, int index,
    double raw_energy, double* const* parameter) {
    int n = law->tab.n_eout[table];
    if (index < 0) index = 0;
    if (index >= n - 1 || index < law->tab.n_discrete[table] ||
        law->tab.interpolation[table] == 1)
        return parameter[table][index];
    double x0 = law->tab.eout[table][index];
    double x1 = law->tab.eout[table][index + 1];
    double f = x1 > x0 ? (raw_energy - x0) / (x1 - x0) : 0.0;
    return parameter[table][index] +
           f * (parameter[table][index + 1] - parameter[table][index]);
}

static alea_error_t sample_kalbach_mu(double r, double a,
                                      alea_nuc_random_fn random,
                                      void* context, double* mu) {
    double choose, u;
    r = fmin(1.0, fmax(0.0, r));
    alea_error_t err = draw_uniform(random, context, &choose);
    if (err != ALEA_OK) return err;
    if ((err = draw_uniform(random, context, &u)) != ALEA_OK) return err;
    if (a < 1e-10) {
        *mu = 2.0 * u - 1.0;
        return ALEA_OK;
    }
    double tail = exp(-2.0 * a);
    double magnitude = 1.0 + log(u + (1.0 - u) * tail) / a;
    *mu = choose < 0.5 * (1.0 + r) ? magnitude : -magnitude;
    *mu = fmin(1.0, fmax(-1.0, *mu));
    return ALEA_OK;
}

static alea_error_t sample_correlated_mu(
    const alea_nuc_energy_dist_t* law, int table, int index,
    double raw_energy, alea_nuc_random_fn random, void* context, double* mu) {
    int selected = index;
    int n = law->tab.n_eout[table];
    if (index >= law->tab.n_discrete[table] && index + 1 < n &&
        law->tab.interpolation[table] == 2) {
        double x0 = law->tab.eout[table][index];
        double x1 = law->tab.eout[table][index + 1];
        double f = x1 > x0 ? (raw_energy - x0) / (x1 - x0) : 0.0;
        double choose;
        alea_error_t err = draw_uniform(random, context, &choose);
        if (err != ALEA_OK) return err;
        if (choose < f) selected++;
    }
    const alea_nuc_angular_point_t* point =
        &law->tab.correlated_mu[table][selected];
    double u;
    alea_error_t err = draw_uniform(random, context, &u);
    if (err != ALEA_OK) return err;
    if (point->type == ALEA_NUC_ANG_ISOTROPIC) {
        *mu = 2.0 * u - 1.0;
        return ALEA_OK;
    }
    return alea_nuc_tabular_pdf_sample(
        point->cosine, point->pdf, point->cdf, point->n_cosines,
        point->interpolation, 0, u, mu, NULL);
}

alea_error_t alea_nuc_sample_energy_angle_distribution(
    const alea_nuc_energy_dist_t* distribution, double incident_energy,
    alea_nuc_random_fn random, void* random_context, double* energy_out,
    double* mu_out, bool* angle_is_correlated) {
    if (!distribution || !random || !energy_out || !mu_out ||
        !angle_is_correlated) return ALEA_ERR_NULL_ARG;
    const alea_nuc_energy_dist_t* selected;
    alea_error_t err = choose_law(distribution, incident_energy, random,
                                  random_context, &selected);
    if (err != ALEA_OK) return err;
    double energy;
    if (selected->law != ALEA_NUC_ELAW_KALBACH &&
        selected->law != ALEA_NUC_ELAW_CORRELATED &&
        selected->law != ALEA_NUC_ELAW_NBODY &&
        selected->law != ALEA_NUC_ELAW_LAB_ANGLE_ENERGY) {
        err = sample_one(selected, incident_energy, random, random_context,
                         &energy);
        if (err != ALEA_OK) return err;
        *energy_out = energy;
        *angle_is_correlated = false;
        return ALEA_OK;
    }
    if (selected->law == ALEA_NUC_ELAW_NBODY) {
        err = sample_nbody(selected, incident_energy, random, random_context,
                           &energy);
        if (err != ALEA_OK) return err;
        double u;
        if ((err = draw_uniform(random, random_context, &u)) != ALEA_OK)
            return err;
        *energy_out = energy;
        *mu_out = 2.0 * u - 1.0;
        *angle_is_correlated = true;
        return ALEA_OK;
    }
    if (selected->law == ALEA_NUC_ELAW_LAB_ANGLE_ENERGY) {
        double mu;
        err = sample_law67(selected, incident_energy, random, random_context,
                           &energy, &mu);
        if (err != ALEA_OK) return err;
        *energy_out = energy;
        *mu_out = mu;
        *angle_is_correlated = true;
        return ALEA_OK;
    }

    int table, index;
    double raw_energy;
    err = sample_continuous_tabular_details(
        selected, incident_energy, random, random_context, &energy,
        &table, &index, &raw_energy);
    if (err != ALEA_OK) return err;
    double mu;
    if (selected->law == ALEA_NUC_ELAW_KALBACH) {
        double r = interpolate_outgoing_parameter(
            selected, table, index, raw_energy, selected->tab.precompound_r);
        double a = interpolate_outgoing_parameter(
            selected, table, index, raw_energy, selected->tab.precompound_a);
        err = sample_kalbach_mu(r, a, random, random_context, &mu);
    } else {
        err = sample_correlated_mu(selected, table, index, raw_energy,
                                   random, random_context, &mu);
    }
    if (err != ALEA_OK) return err;
    *energy_out = energy;
    *mu_out = mu;
    *angle_is_correlated = true;
    return ALEA_OK;
}

alea_error_t alea_nuc_sample_energy_distribution(
    const alea_nuc_energy_dist_t* distribution, double incident_energy,
    alea_nuc_random_fn random, void* random_context, double* energy_out) {
    if (!distribution || !random || !energy_out) return ALEA_ERR_NULL_ARG;
    const alea_nuc_energy_dist_t* selected;
    alea_error_t err = choose_law(distribution, incident_energy, random,
                                  random_context, &selected);
    if (err != ALEA_OK) return err;
    double sampled;
    err = sample_one(selected, incident_energy, random, random_context,
                     &sampled);
    if (err == ALEA_OK) *energy_out = sampled;
    return err;
}
