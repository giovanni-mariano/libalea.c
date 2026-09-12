// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/** @file photon_sample.c @brief Photoatomic collision sampling. */

#include "nuclear_internal.h"
#include <float.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define ELECTRON_REST_MEV 0.51099895069
#define INV_PLANCK_C_PER_ANGSTROM_MEV 80.65543896
#define PHOTON_MAX_ATTEMPTS 100000

static alea_error_t uniform(alea_nuc_random_fn random, void* context,
                            double* value) {
    double u = random(context);
    if (!isfinite(u) || u < 0.0 || u >= 1.0)
        return ALEA_ERR_INVALID_STATE;
    *value = u;
    return ALEA_OK;
}

static int particle_valid(const alea_nuc_particle_state_t* particle) {
    if (!particle || particle->type != ALEA_NUC_PARTICLE_PHOTON ||
        !isfinite(particle->energy) || particle->energy <= 0.0 ||
        !isfinite(particle->weight) || particle->weight <= 0.0 ||
        !isfinite(particle->time)) return 0;
    double norm2 = 0.0;
    for (int i = 0; i < 3; i++) {
        if (!isfinite(particle->direction[i])) return 0;
        norm2 += particle->direction[i] * particle->direction[i];
    }
    return fabs(norm2 - 1.0) <= 1e-8;
}

static int factor_valid(const double* x, const double* y, int n,
                        int cumulative) {
    if (!x || !y || n < 2) return 0;
    for (int i = 0; i < n; i++)
        if (!isfinite(x[i]) || !isfinite(y[i]) || x[i] < 0.0 || y[i] < 0.0 ||
            (i > 0 && (x[i] <= x[i - 1] ||
                       (cumulative && y[i] < y[i - 1])))) return 0;
    return 1;
}

static double interpolate(const double* x, const double* y, int n,
                          double query) {
    if (query <= x[0]) return y[0];
    if (query >= x[n - 1]) return y[n - 1];
    double f;
    int i = alea_nuc_energy_lookup_trusted(x, n, query, &f);
    return y[i] + f * (y[i + 1] - y[i]);
}

static double interpolate_square(const double* x, const double* y, int n,
                                 double query) {
    if (query <= x[0]) return y[0];
    if (query >= x[n - 1]) return y[n - 1];
    double ignored;
    int i = alea_nuc_energy_lookup_trusted(x, n, query, &ignored);
    double q2 = query * query;
    double x02 = x[i] * x[i], x12 = x[i + 1] * x[i + 1];
    double f = x12 > x02 ? (q2 - x02) / (x12 - x02) : 0.0;
    return y[i] + f * (y[i + 1] - y[i]);
}

static double inverse_cumulative_square(const double* x, const double* cdf,
                                        int n, double value) {
    if (value <= cdf[0]) return x[0];
    int lo = 0, hi = n - 1;
    while (hi - lo > 1) {
        int mid = lo + (hi - lo) / 2;
        if (cdf[mid] <= value) lo = mid;
        else hi = mid;
    }
    if (cdf[hi] <= cdf[lo]) return x[lo];
    double f = (value - cdf[lo]) / (cdf[hi] - cdf[lo]);
    double x2 = x[lo] * x[lo] +
                f * (x[hi] * x[hi] - x[lo] * x[lo]);
    return sqrt(fmax(0.0, x2));
}

static alea_error_t sample_coherent(const alea_nuc_photon_data_t* photon,
                                    double energy, alea_nuc_random_fn random,
                                    void* context, double* mu) {
    int n = photon->n_coherent_ff;
    if (!factor_valid(photon->coherent_momentum,
                      photon->coherent_ff_cumulative, n, 1))
        return ALEA_ERR_INVALID_STATE;
    double xmax = energy * INV_PLANCK_C_PER_ANGSTROM_MEV;
    if (!(xmax > 0.0) || !isfinite(xmax)) return ALEA_ERR_INVALID_ARG;
    double capped = fmin(xmax, photon->coherent_momentum[n - 1]);
    double integral = interpolate_square(photon->coherent_momentum,
                                         photon->coherent_ff_cumulative,
                                         n, capped);
    if (!(integral > 0.0)) return ALEA_ERR_INVALID_STATE;
    for (int attempt = 0; attempt < PHOTON_MAX_ATTEMPTS; attempt++) {
        double u, accept;
        alea_error_t err = uniform(random, context, &u);
        if (err != ALEA_OK) return err;
        double x = inverse_cumulative_square(
            photon->coherent_momentum, photon->coherent_ff_cumulative,
            n, u * integral);
        double candidate = 1.0 - 2.0 * (x / xmax) * (x / xmax);
        candidate = fmin(1.0, fmax(-1.0, candidate));
        if ((err = uniform(random, context, &accept)) != ALEA_OK) return err;
        if (accept < 0.5 * (1.0 + candidate * candidate)) {
            *mu = candidate;
            return ALEA_OK;
        }
    }
    return ALEA_ERR_INVALID_STATE;
}

static alea_error_t sample_incoherent(const alea_nuc_photon_data_t* photon,
                                      double energy, alea_nuc_random_fn random,
                                      void* context, double* mu,
                                      double* energy_out, int* vacancy) {
    *vacancy = -1;
    int n = photon->n_incoherent_ff;
    if (!factor_valid(photon->incoherent_momentum,
                      photon->incoherent_ff, n, 0))
        return ALEA_ERR_INVALID_STATE;
    double alpha = energy / ELECTRON_REST_MEV;
    double xmax = energy * INV_PLANCK_C_PER_ANGSTROM_MEV;
    double smax = interpolate(photon->incoherent_momentum,
                              photon->incoherent_ff, n, xmax);
    if (!(smax > 0.0)) return ALEA_ERR_INVALID_STATE;
    int angle_accepted = 0;
    for (int attempt = 0; attempt < PHOTON_MAX_ATTEMPTS; attempt++) {
        double angle, kn_accept, factor_accept;
        alea_error_t err = uniform(random, context, &angle);
        if (err != ALEA_OK) return err;
        double candidate = 2.0 * angle - 1.0;
        double ratio = 1.0 / (1.0 + alpha * (1.0 - candidate));
        double kernel = ratio * ratio *
            (ratio + 1.0 / ratio - (1.0 - candidate * candidate));
        if ((err = uniform(random, context, &kn_accept)) != ALEA_OK) return err;
        if (kn_accept >= 0.5 * kernel) continue;
        double x = xmax * sqrt(0.5 * (1.0 - candidate));
        double factor = interpolate(photon->incoherent_momentum,
                                    photon->incoherent_ff, n, x);
        if ((err = uniform(random, context, &factor_accept)) != ALEA_OK)
            return err;
        if (factor_accept < factor / smax) {
            *mu = candidate;
            angle_accepted = 1;
            break;
        }
    }
    if (!angle_accepted) return ALEA_ERR_INVALID_STATE;
    if (photon->n_compton_profiles == 0) {
        *energy_out = energy / (1.0 + alpha * (1.0 - *mu));
        return ALEA_OK;
    }
    if (photon->n_compton_shells <= 0 || !photon->compton_shells ||
        !photon->compton_profiles) return ALEA_ERR_INVALID_STATE;

    /* Namito impulse-approximation sampling as used with the ACE SWD block.
     * Momentum is tabulated in atomic units, hence the fine-structure factor. */
    const double inv_rest = 1.0 / ELECTRON_REST_MEV;
    const double inv_alpha = 137.03605;
    for (int attempt = 0; attempt < PHOTON_MAX_ATTEMPTS; attempt++) {
        double u;
        alea_error_t err = uniform(random, context, &u);
        if (err != ALEA_OK) return err;
        int selected = photon->n_compton_shells - 1;
        for (int i = 0; i < photon->n_compton_shells; i++)
            if (u < photon->compton_shells[i].cumulative_probability) {
                selected = i;
                break;
            }
        const alea_nuc_compton_shell_t* shell =
            &photon->compton_shells[selected];
        int selected_vacancy = -1;
        for (int i = 0; i < photon->n_subshells; i++)
            if (u < photon->subshells[i].compton_vacancy_probability) {
                selected_vacancy = i;
                break;
            }
        double maximum = energy - shell->binding_energy;
        if (!(maximum > 0.0)) continue;
        if (shell->profile_index < 0 ||
            shell->profile_index >= photon->n_compton_profiles)
            return ALEA_ERR_INVALID_STATE;
        const alea_nuc_compton_profile_t* profile =
            &photon->compton_profiles[shell->profile_index];
        if (profile->n_momenta < 2 || !profile->momentum ||
            !profile->pdf || !profile->cdf)
            return ALEA_ERR_INVALID_STATE;
        if ((err = uniform(random, context, &u)) != ALEA_OK) return err;
        int bin = profile->n_momenta - 2;
        for (int i = 0; i < profile->n_momenta - 1; i++)
            if (u < profile->cdf[i + 1]) { bin = i; break; }
        if (!(profile->pdf[bin] > 0.0)) continue;
        double pz = profile->momentum[bin] +
                    (u - profile->cdf[bin]) / profile->pdf[bin];
        if (!isfinite(pz) || pz < profile->momentum[bin] - 1e-12 ||
            pz > profile->momentum[bin + 1] * (1.0 + 1e-10)) continue;

        double p2 = (pz / inv_alpha) * (pz / inv_alpha);
        double transfer = energy * (1.0 - *mu) * inv_rest;
        double a = p2 - 1.0 - transfer * transfer - 2.0 * transfer;
        double b = 2.0 * energy +
            2.0 * energy * energy * (1.0 - *mu) * inv_rest -
            2.0 * p2 * energy * *mu;
        double c = p2 * energy * energy - energy * energy;
        double discriminant = b * b - 4.0 * a * c;
        if (!(discriminant >= 0.0) || !isfinite(discriminant) ||
            fabs(a) <= DBL_EPSILON * (fabs(b) + fabs(c) + 1.0)) continue;
        if ((err = uniform(random, context, &u)) != ALEA_OK) return err;
        double scattered = (-b + (u < 0.5 ? 1.0 : -1.0) *
                                      sqrt(discriminant)) / (2.0 * a);
        if (!isfinite(scattered) || scattered < 0.0 || scattered > maximum)
            continue;
        if ((err = uniform(random, context, &u)) != ALEA_OK) return err;
        if (scattered >= maximum * u) {
            *energy_out = scattered;
            *vacancy = selected_vacancy;
            return ALEA_OK;
        }
    }
    return ALEA_ERR_INVALID_STATE;
}

static int shell_index(const alea_nuc_photon_data_t* photon, int designator) {
    for (int i = 0; i < photon->n_subshells; i++)
        if (photon->subshells[i].designator == designator) return i;
    return -1;
}

size_t alea_nuc_photon_secondary_capacity(
    const alea_nuc_nuclide_t* element, double energy) {
    if (!element || !element->photon || energy <= 0.0 || !isfinite(energy))
        return 0;
    size_t required = alea_nuc_photon_xs_pair(element, energy) > 0.0 ? 2u : 0u;
    const alea_nuc_photon_data_t* photon = element->photon;
    if (photon->n_subshells < 0 || photon->n_subshells > 128 ||
        (photon->n_subshells > 0 && !photon->subshells)) return 0;
    for (int i = 0; i < photon->n_subshells; i++) {
        const alea_nuc_atomic_subshell_t* shell = &photon->subshells[i];
        if (alea_nuc_photon_xs_photoelectric_subshell(
                element, shell->designator, energy) > 0.0 &&
            shell->max_relaxation_photons > required)
            required = shell->max_relaxation_photons;
    }
    if (alea_nuc_photon_xs_incoherent(element, energy) > 0.0) {
        double previous = 0.0;
        for (int i = 0; i < photon->n_subshells; i++) {
            const alea_nuc_atomic_subshell_t* shell = &photon->subshells[i];
            if (shell->compton_vacancy_probability > previous &&
                shell->max_relaxation_photons > required)
                required = shell->max_relaxation_photons;
            previous = shell->compton_vacancy_probability;
        }
    }
    return required;
}

static alea_error_t sample_relaxation_from_shell(
    const alea_nuc_photon_data_t* photon, int initial,
    const alea_nuc_particle_state_t* incident,
    alea_nuc_random_fn random, void* random_context,
    alea_nuc_secondary_buffer_t* secondaries,
    size_t* emitted, double* emitted_energy) {
    *emitted = 0;
    *emitted_energy = 0.0;
    if (initial < 0 || initial >= photon->n_subshells)
        return ALEA_ERR_INVALID_STATE;
    alea_error_t err;
    int vacancies[129];
    size_t n_vacancies = 1;
    vacancies[0] = initial;
    size_t start = secondaries->count;
    size_t count = 0;
    double energy_sum = 0.0;
    while (n_vacancies > 0) {
        const alea_nuc_atomic_subshell_t* shell =
            &photon->subshells[vacancies[--n_vacancies]];
        if (shell->n_transitions == 0) continue;
        if (shell->n_transitions < 0 || !shell->transitions)
            return ALEA_ERR_INVALID_STATE;
        double u;
        if ((err = uniform(random, random_context, &u)) != ALEA_OK) return err;
        int chosen = shell->n_transitions - 1;
        for (int i = 0; i < shell->n_transitions; i++)
            if (u < shell->transitions[i].cumulative_probability) {
                chosen = i;
                break;
            }
        const alea_nuc_atomic_transition_t* transition =
            &shell->transitions[chosen];
        int primary = shell_index(photon, transition->primary_designator);
        int secondary = transition->secondary_designator == 0 ? -1 :
            shell_index(photon, transition->secondary_designator);
        if (primary < 0 || n_vacancies >= 128) return ALEA_ERR_INVALID_STATE;
        vacancies[n_vacancies++] = primary;
        if (secondary >= 0) {
            if (n_vacancies >= 129) return ALEA_ERR_INVALID_STATE;
            vacancies[n_vacancies++] = secondary;
        } else if (transition->energy > 0.0) {
            if (count >= secondaries->capacity - start)
                return ALEA_ERR_INVALID_STATE;
            double polar, azimuth;
            if ((err = uniform(random, random_context, &polar)) != ALEA_OK ||
                (err = uniform(random, random_context, &azimuth)) != ALEA_OK)
                return err;
            alea_nuc_particle_state_t particle = *incident;
            particle.energy = transition->energy;
            double mu = 2.0 * polar - 1.0;
            double radial = sqrt(fmax(0.0, 1.0 - mu * mu));
            particle.direction[0] = radial * cos(2.0 * M_PI * azimuth);
            particle.direction[1] = radial * sin(2.0 * M_PI * azimuth);
            particle.direction[2] = mu;
            secondaries->particles[start + count] = particle;
            count++;
            energy_sum += transition->energy;
        }
    }
    if (!isfinite(energy_sum) ||
        energy_sum > incident->energy * (1.0 + 1e-12))
        return ALEA_ERR_INVALID_STATE;
    secondaries->count = start + count;
    *emitted = count;
    *emitted_energy = energy_sum;
    return ALEA_OK;
}

static alea_error_t sample_photoelectric_relaxation(
    const alea_nuc_nuclide_t* element,
    const alea_nuc_particle_state_t* incident,
    alea_nuc_random_fn random, void* random_context,
    alea_nuc_secondary_buffer_t* secondaries,
    size_t* emitted, double* emitted_energy) {
    const alea_nuc_photon_data_t* photon = element->photon;
    double total = 0.0;
    for (int i = 0; i < photon->n_subshells; i++)
        total += alea_nuc_photon_xs_photoelectric_subshell(
            element, photon->subshells[i].designator, incident->energy);
    *emitted = 0;
    *emitted_energy = 0.0;
    if (!(total > 0.0)) return ALEA_OK;
    double selection;
    alea_error_t err = uniform(random, random_context, &selection);
    if (err != ALEA_OK) return err;
    double threshold = selection * total, sum = 0.0;
    int initial = photon->n_subshells - 1;
    for (int i = 0; i < photon->n_subshells; i++) {
        sum += alea_nuc_photon_xs_photoelectric_subshell(
            element, photon->subshells[i].designator, incident->energy);
        if (threshold < sum) { initial = i; break; }
    }
    return sample_relaxation_from_shell(
        photon, initial, incident, random, random_context, secondaries,
        emitted, emitted_energy);
}

static alea_error_t sample_compton_relaxation(
    const alea_nuc_photon_data_t* photon, int vacancy,
    const alea_nuc_particle_state_t* incident,
    alea_nuc_random_fn random, void* random_context,
    alea_nuc_secondary_buffer_t* secondaries,
    size_t* emitted, double* emitted_energy) {
    *emitted = 0;
    *emitted_energy = 0.0;
    if (vacancy < 0) return ALEA_OK;
    return sample_relaxation_from_shell(
        photon, vacancy, incident, random, random_context, secondaries,
        emitted, emitted_energy);
}

static alea_error_t sample_photon_collision(
    const alea_nuc_nuclide_t* element,
    const alea_nuc_particle_state_t* incident,
    alea_nuc_random_fn random, void* random_context,
    alea_nuc_secondary_buffer_t* secondaries,
    alea_nuc_collision_result_t* result) {
    if (!element || !incident || !random || !result) return ALEA_ERR_NULL_ARG;
    if (!element->photon || element->particle != ALEA_NUC_PARTICLE_PHOTON ||
        !particle_valid(incident)) return ALEA_ERR_INVALID_ARG;
    double xs[] = {
        alea_nuc_photon_xs_coherent(element, incident->energy),
        alea_nuc_photon_xs_incoherent(element, incident->energy),
        alea_nuc_photon_xs_photoelectric(element, incident->energy),
        alea_nuc_photon_xs_pair(element, incident->energy)
    };
    double total = 0.0;
    for (int i = 0; i < 4; i++) {
        if (!isfinite(xs[i]) || xs[i] < 0.0) return ALEA_ERR_INVALID_STATE;
        total += xs[i];
    }
    if (!(total > 0.0) || !isfinite(total)) return ALEA_ERR_INVALID_STATE;
    if (xs[3] > 0.0 && incident->energy < 2.0 * ELECTRON_REST_MEV)
        return ALEA_ERR_INVALID_STATE;
    if (secondaries) {
        if (secondaries->count > secondaries->capacity ||
            (secondaries->capacity > 0 && !secondaries->particles))
            return ALEA_ERR_INVALID_STATE;
        size_t required = alea_nuc_photon_secondary_capacity(
            element, incident->energy);
        if (secondaries->capacity - secondaries->count < required)
            return ALEA_ERR_OUT_OF_MEMORY;
    }
    double select;
    alea_error_t err = uniform(random, random_context, &select);
    if (err != ALEA_OK) return err;
    double threshold = select * total, sum = 0.0;
    int channel = 3;
    for (int i = 0; i < 4; i++) {
        sum += xs[i];
        if (threshold < sum) { channel = i; break; }
    }

    static const int mts[] = {502, 504, 522, 517};
    alea_nuc_collision_result_t candidate;
    memset(&candidate, 0, sizeof(candidate));
    candidate.component_index = -1;
    candidate.mt = mts[channel];
    candidate.deposition_available = true;
    if (channel == 2 && secondaries && element->photon->n_subshells > 0) {
        size_t emitted;
        double emitted_energy;
        err = sample_photoelectric_relaxation(
            element, incident, random, random_context, secondaries,
            &emitted, &emitted_energy);
        if (err != ALEA_OK) return err;
        candidate.outcome = emitted > 0 ? ALEA_NUC_OUTCOME_REPLACED
                                        : ALEA_NUC_OUTCOME_ABSORBED;
        candidate.n_emitted = emitted;
        candidate.local_energy_deposition =
            fmax(0.0, incident->energy - emitted_energy);
        *result = candidate;
        return ALEA_OK;
    }
    if (channel == 2 || (channel == 3 && !secondaries)) {
        candidate.outcome = ALEA_NUC_OUTCOME_ABSORBED;
        candidate.local_energy_deposition = incident->energy;
        *result = candidate;
        return ALEA_OK;
    }

    if (channel == 3) {
        double polar, azimuth;
        if ((err = uniform(random, random_context, &polar)) != ALEA_OK ||
            (err = uniform(random, random_context, &azimuth)) != ALEA_OK)
            return err;
        alea_nuc_particle_state_t photons[2] = {*incident, *incident};
        photons[0].energy = ELECTRON_REST_MEV;
        photons[1].energy = ELECTRON_REST_MEV;
        alea_nuc_rotate_direction_internal(
            incident->direction, 2.0 * polar - 1.0, 2.0 * M_PI * azimuth,
            photons[0].direction);
        for (int i = 0; i < 3; i++)
            photons[1].direction[i] = -photons[0].direction[i];
        size_t start = secondaries->count;
        secondaries->particles[start] = photons[0];
        secondaries->particles[start + 1] = photons[1];
        secondaries->count = start + 2;
        candidate.outcome = ALEA_NUC_OUTCOME_REPLACED;
        candidate.n_emitted = 2;
        candidate.local_energy_deposition =
            fmax(0.0, incident->energy - 2.0 * ELECTRON_REST_MEV);
        *result = candidate;
        return ALEA_OK;
    }

    double mu, energy_out = incident->energy;
    int compton_vacancy = -1;
    if (channel == 0)
        err = sample_coherent(element->photon, incident->energy, random,
                              random_context, &mu);
    else
        err = sample_incoherent(element->photon, incident->energy, random,
                                random_context, &mu, &energy_out,
                                &compton_vacancy);
    if (err != ALEA_OK) return err;
    double azimuth;
    if ((err = uniform(random, random_context, &azimuth)) != ALEA_OK) return err;
    candidate.outcome = ALEA_NUC_OUTCOME_SCATTERED;
    candidate.mu_lab = mu;
    candidate.mu_cm = NAN;
    candidate.outgoing = *incident;
    candidate.outgoing.energy = energy_out;
    alea_nuc_rotate_direction_internal(
        incident->direction, mu, 2.0 * M_PI * azimuth,
        candidate.outgoing.direction);
    double emitted_energy = 0.0;
    if (channel == 1 && secondaries && element->photon->n_subshells > 0) {
        size_t emitted = 0;
        err = sample_compton_relaxation(
            element->photon, compton_vacancy, incident, random,
            random_context, secondaries, &emitted, &emitted_energy);
        if (err != ALEA_OK) return err;
        candidate.n_emitted = emitted;
    }
    if (emitted_energy > incident->energy - energy_out + 1e-10)
        return ALEA_ERR_INVALID_STATE;
    candidate.local_energy_deposition =
        fmax(0.0, incident->energy - energy_out - emitted_energy);
    *result = candidate;
    return ALEA_OK;
}

alea_error_t alea_nuc_sample_photon_collision(
    const alea_nuc_nuclide_t* element,
    const alea_nuc_particle_state_t* incident,
    alea_nuc_random_fn random, void* random_context,
    alea_nuc_collision_result_t* result) {
    return sample_photon_collision(element, incident, random, random_context,
                                   NULL, result);
}

alea_error_t alea_nuc_sample_photon_collision_with_secondaries(
    const alea_nuc_nuclide_t* element,
    const alea_nuc_particle_state_t* incident,
    alea_nuc_random_fn random, void* random_context,
    alea_nuc_secondary_buffer_t* secondaries,
    alea_nuc_collision_result_t* result) {
    if (!secondaries) return ALEA_ERR_NULL_ARG;
    return sample_photon_collision(element, incident, random, random_context,
                                   secondaries, result);
}
