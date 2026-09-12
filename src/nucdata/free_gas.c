// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/**
 * @file free_gas.c
 * @brief Collision-conditioned free-gas elastic scattering
 */

#include "nuclear_internal.h"

#include <float.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define ALEA_NUC_FREE_GAS_MAX_ATTEMPTS 100000

static alea_error_t uniform(alea_nuc_random_fn random, void* context,
                            double* value) {
    if (!random || !value) return ALEA_ERR_NULL_ARG;
    double u = random(context);
    if (!isfinite(u) || u < 0.0 || u >= 1.0)
        return ALEA_ERR_INVALID_STATE;
    *value = u;
    return ALEA_OK;
}

static double positive(double value) {
    return value == 0.0 ? DBL_MIN : value;
}

static alea_error_t sample_maxwell_x(alea_nuc_random_fn random, void* context,
                                     double* x) {
    double u1, u2, u3;
    alea_error_t err = uniform(random, context, &u1);
    if (err != ALEA_OK) return err;
    if ((err = uniform(random, context, &u2)) != ALEA_OK) return err;
    if ((err = uniform(random, context, &u3)) != ALEA_OK) return err;
    double c = cos(0.5 * M_PI * u3);
    *x = sqrt(-log(positive(u1)) - log(positive(u2)) * c * c);
    return ALEA_OK;
}

static alea_error_t sample_target_velocity(
    const alea_nuc_particle_state_t* incident, double awr, double kT,
    alea_nuc_random_fn random, void* context, double target[3]) {
    double neutron_speed = sqrt(incident->energy);
    double y = sqrt(awr * incident->energy / kT);
    double branch_probability = 2.0 / (sqrt(M_PI) * y + 2.0);

    for (int attempt = 0; attempt < ALEA_NUC_FREE_GAS_MAX_ATTEMPTS;
         attempt++) {
        double choose, x;
        alea_error_t err = uniform(random, context, &choose);
        if (err != ALEA_OK) return err;
        if (choose < branch_probability) {
            double u1, u2;
            if ((err = uniform(random, context, &u1)) != ALEA_OK) return err;
            if ((err = uniform(random, context, &u2)) != ALEA_OK) return err;
            x = sqrt(-log(positive(u1)) - log(positive(u2)));
        } else {
            if ((err = sample_maxwell_x(random, context, &x)) != ALEA_OK)
                return err;
        }
        double target_speed = x * sqrt(kT / awr);
        double mu, phi_draw, accept;
        if ((err = uniform(random, context, &mu)) != ALEA_OK) return err;
        if ((err = uniform(random, context, &phi_draw)) != ALEA_OK) return err;
        if ((err = uniform(random, context, &accept)) != ALEA_OK) return err;
        mu = 2.0 * mu - 1.0;
        double relative = sqrt(
            neutron_speed * neutron_speed + target_speed * target_speed -
            2.0 * neutron_speed * target_speed * mu);
        if (accept >= relative / (neutron_speed + target_speed)) continue;
        double direction[3];
        alea_nuc_rotate_direction_internal(
            incident->direction, mu, 2.0 * M_PI * phi_draw, direction);
        for (int i = 0; i < 3; i++) target[i] = target_speed * direction[i];
        return ALEA_OK;
    }
    return ALEA_ERR_INVALID_STATE;
}

alea_error_t alea_nuc_sample_free_gas_elastic(
    const alea_nuc_nuclide_t* nuc,
    const alea_nuc_particle_state_t* incident, double kT,
    alea_nuc_random_fn random, void* random_context,
    alea_nuc_free_gas_result_t* result) {
    if (!nuc || !incident || !random || !result) return ALEA_ERR_NULL_ARG;
    if (incident->type != ALEA_NUC_PARTICLE_NEUTRON ||
        !isfinite(incident->energy) || incident->energy <= 0.0 ||
        !isfinite(incident->weight) || incident->weight <= 0.0 ||
        !isfinite(incident->time) ||
        !isfinite(kT) || kT <= 0.0 || !isfinite(nuc->awr) || nuc->awr <= 0.0)
        return ALEA_ERR_INVALID_ARG;
    double norm2 = 0.0;
    for (int i = 0; i < 3; i++) {
        if (!isfinite(incident->direction[i])) return ALEA_ERR_INVALID_ARG;
        norm2 += incident->direction[i] * incident->direction[i];
    }
    if (fabs(norm2 - 1.0) > 1e-8) return ALEA_ERR_INVALID_ARG;
    if (!alea_nuc_validate_angular_internal(nuc->elastic_angular))
        return ALEA_ERR_INVALID_STATE;

    double target[3];
    alea_error_t err = sample_target_velocity(
        incident, nuc->awr, kT, random, random_context, target);
    if (err != ALEA_OK) return err;
    double neutron[3], center[3], neutron_cm[3];
    double cm_speed2 = 0.0, relative_energy = 0.0;
    for (int i = 0; i < 3; i++) {
        neutron[i] = sqrt(incident->energy) * incident->direction[i];
        center[i] = (neutron[i] + nuc->awr * target[i]) / (nuc->awr + 1.0);
        neutron_cm[i] = neutron[i] - center[i];
        cm_speed2 += neutron_cm[i] * neutron_cm[i];
        double relative = neutron[i] - target[i];
        relative_energy += relative * relative;
    }
    if (!(cm_speed2 > 0.0) || !isfinite(cm_speed2))
        return ALEA_ERR_INVALID_STATE;
    double cm_speed = sqrt(cm_speed2);
    double cm_direction[3];
    for (int i = 0; i < 3; i++) cm_direction[i] = neutron_cm[i] / cm_speed;

    double select, sample, azimuth;
    if ((err = uniform(random, random_context, &select)) != ALEA_OK ||
        (err = uniform(random, random_context, &sample)) != ALEA_OK ||
        (err = uniform(random, random_context, &azimuth)) != ALEA_OK)
        return err;
    double mu_cm = alea_nuc_sample_angular_mu_internal(
        nuc->elastic_angular, relative_energy, select, sample);
    double outgoing_cm_direction[3];
    alea_nuc_rotate_direction_internal(
        cm_direction, mu_cm, 2.0 * M_PI * azimuth, outgoing_cm_direction);
    double outgoing_velocity[3], outgoing_energy = 0.0;
    for (int i = 0; i < 3; i++) {
        outgoing_velocity[i] = center[i] + cm_speed * outgoing_cm_direction[i];
        outgoing_energy += outgoing_velocity[i] * outgoing_velocity[i];
    }
    if (!(outgoing_energy > 0.0) || !isfinite(outgoing_energy))
        return ALEA_ERR_INVALID_STATE;

    alea_nuc_free_gas_result_t candidate;
    memset(&candidate, 0, sizeof(candidate));
    candidate.outgoing = *incident;
    candidate.outgoing.energy = outgoing_energy;
    candidate.mu_cm = mu_cm;
    candidate.mu_lab = 0.0;
    double outgoing_speed = sqrt(outgoing_energy);
    for (int i = 0; i < 3; i++) {
        candidate.outgoing.direction[i] = outgoing_velocity[i] / outgoing_speed;
        candidate.mu_lab += incident->direction[i] *
                            candidate.outgoing.direction[i];
    }
    candidate.mu_lab = fmin(1.0, fmax(-1.0, candidate.mu_lab));
    *result = candidate;
    return ALEA_OK;
}
