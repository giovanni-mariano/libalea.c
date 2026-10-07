// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

#ifndef ALEA_CELL_PARAMETER_H
#define ALEA_CELL_PARAMETER_H

#include "alea_model.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

static inline int alea_parameter_nonnegative(double value) {
    return isfinite(value) && value >= 0;
}
static inline int alea_parameter_probability(double value) {
    return alea_parameter_nonnegative(value) && value <= 1;
}

static inline int alea_secondary_state_valid(alea_secondary_collision_state_t state) {
    return state == ALEA_SECONDARY_COLLIDED || state == ALEA_SECONDARY_UNCOLLIDED;
}

static inline int alea_cell_parameters_valid(const alea_model_cell_metadata_t* m) {
    if (!m) return 1;
    if ((m->parameter_flags & ~((1u << 7) - 1)) ||
        (m->energy_cutoff_particles & ~((1u << ALEA_PARTICLE_COUNT) - 1)) ||
        (m->secondary_state_particles & ~((1u << ALEA_PARTICLE_COUNT) - 1))) return 0;
    if ((m->has_importance_neutron && !alea_parameter_nonnegative(m->importance_neutron)) ||
        (m->has_importance_photon && !alea_parameter_nonnegative(m->importance_photon)) ||
        (m->has_importance_electron && !alea_parameter_nonnegative(m->importance_electron))) return 0;
    if ((m->parameter_flags & ALEA_CELL_PARAM_VOLUME) && !alea_parameter_nonnegative(m->user_volume)) return 0;
    if (m->parameter_flags & ALEA_CELL_PARAM_PHOTON_PRODUCTION) {
        const alea_photon_production_t* p = &m->photon_production;
        if (p->mode == ALEA_PHOTON_PRODUCTION_THRESHOLD) {
            if (!alea_parameter_nonnegative(p->weight_threshold) || p->weight_threshold == 0 ||
                (p->weight_basis != ALEA_PHOTON_WEIGHT_ABSOLUTE && p->weight_basis != ALEA_PHOTON_WEIGHT_SOURCE_RELATIVE)) return 0;
        } else if (p->mode != ALEA_PHOTON_PRODUCTION_OFF && p->mode != ALEA_PHOTON_PRODUCTION_ONE_PER_COLLISION) return 0;
    }
    if ((m->parameter_flags & ALEA_CELL_PARAM_FISSION_MODE) &&
        m->fission_mode != ALEA_FISSION_NORMAL &&
        m->fission_mode != ALEA_FISSION_CAPTURE_WITH_PHOTONS &&
        m->fission_mode != ALEA_FISSION_CAPTURE_WITHOUT_PHOTONS) return 0;
    if ((m->parameter_flags & ALEA_CELL_PARAM_DETECTOR_PROBABILITY) && !alea_parameter_probability(m->detector_contribution)) return 0;
    if ((m->parameter_flags & ALEA_CELL_PARAM_ENERGY_CUTOFF) && !alea_parameter_nonnegative(m->energy_cutoff)) return 0;
    if ((m->parameter_flags & ALEA_CELL_PARAM_SECONDARY_STATE) && !alea_secondary_state_valid(m->secondary_collision_state)) return 0;
    if ((m->parameter_flags & ALEA_CELL_PARAM_MAGNETIC_FIELD) && m->magnetic_field < 0) return 0;
    for (int i = 0; i < ALEA_PARTICLE_COUNT; i++) {
        if ((m->energy_cutoff_particles & (1u << i)) && !alea_parameter_nonnegative(m->particle_energy_cutoff[i])) return 0;
        if ((m->secondary_state_particles & (1u << i)) && !alea_secondary_state_valid(m->particle_secondary_state[i])) return 0;
    }
    if (m->detector_probability_count && !m->detector_probabilities) return 0;
    for (size_t i = 0; i < m->detector_probability_count; i++) {
        const alea_detector_probability_t* p = &m->detector_probabilities[i];
        if (p->tally <= 0 || !alea_parameter_probability(p->probability)) return 0;
        for (size_t j = 0; j < i; j++)
            if (m->detector_probabilities[j].tally == p->tally) return 0;
    }
    return 1;
}

static inline int alea_photon_production_from_mcnp(
    double value, alea_photon_production_t* production) {
    if (!isfinite(value)) return -1;
    memset(production, 0, sizeof(*production));
    if (value == -1.0e6) production->mode = ALEA_PHOTON_PRODUCTION_OFF;
    else if (value == 0) production->mode = ALEA_PHOTON_PRODUCTION_ONE_PER_COLLISION;
    else {
        production->mode = ALEA_PHOTON_PRODUCTION_THRESHOLD;
        production->weight_basis = value < 0 ? ALEA_PHOTON_WEIGHT_SOURCE_RELATIVE
                                           : ALEA_PHOTON_WEIGHT_ABSOLUTE;
        production->weight_threshold = fabs(value);
    }
    return 0;
}

static inline int alea_photon_production_to_mcnp(
    const alea_photon_production_t* production, double* value) {
    if (production->mode == ALEA_PHOTON_PRODUCTION_OFF) *value = -1.0e6;
    else if (production->mode == ALEA_PHOTON_PRODUCTION_ONE_PER_COLLISION) *value = 0;
    else if (production->mode == ALEA_PHOTON_PRODUCTION_THRESHOLD &&
             isfinite(production->weight_threshold) && production->weight_threshold > 0 &&
             (production->weight_basis == ALEA_PHOTON_WEIGHT_ABSOLUTE ||
              production->weight_basis == ALEA_PHOTON_WEIGHT_SOURCE_RELATIVE)) {
        *value = production->weight_basis == ALEA_PHOTON_WEIGHT_SOURCE_RELATIVE
            ? -production->weight_threshold : production->weight_threshold;
        /* This magnitude collides with MCNP's production-off sentinel. */
        if (*value == -1.0e6) return -1;
    } else return -1;
    return 0;
}

static inline int alea_fission_mode_from_mcnp(int value, alea_fission_mode_t* mode) {
    if (value == 0) *mode = ALEA_FISSION_CAPTURE_WITH_PHOTONS;
    else if (value == 1) *mode = ALEA_FISSION_NORMAL;
    else if (value == 2) *mode = ALEA_FISSION_CAPTURE_WITHOUT_PHOTONS;
    else return -1;
    return 0;
}

static inline int alea_fission_mode_to_mcnp(alea_fission_mode_t mode, int* value) {
    if (mode == ALEA_FISSION_NORMAL) *value = 1;
    else if (mode == ALEA_FISSION_CAPTURE_WITH_PHOTONS) *value = 0;
    else if (mode == ALEA_FISSION_CAPTURE_WITHOUT_PHOTONS) *value = 2;
    else return -1;
    return 0;
}

/* Sorted insertion keeps exports deterministic and supports per-tally overrides. */
static inline int alea_detector_probability_store(
    alea_detector_probability_t** entries, size_t* count,
    int tally, double probability) {
    size_t index = 0;
    while (index < *count && (*entries)[index].tally < tally) index++;
    if (index < *count && (*entries)[index].tally == tally) {
        (*entries)[index].probability = probability;
        return 0;
    }
    if (*count >= SIZE_MAX / sizeof(**entries)) return -1;
    alea_detector_probability_t* grown = realloc(
        *entries, (*count + 1) * sizeof(**entries));
    if (!grown) return -1;
    *entries = grown;
    memmove(grown + index + 1, grown + index,
            (*count - index) * sizeof(*grown));
    grown[index] = (alea_detector_probability_t){tally, probability};
    (*count)++;
    return 0;
}

static inline alea_detector_probability_t* alea_detector_probability_copy(
    const alea_detector_probability_t* entries, size_t count) {
    if (!count) return NULL;
    if (!entries || count > SIZE_MAX / sizeof(*entries)) return NULL;
    alea_detector_probability_t* copy = malloc(count * sizeof(*copy));
    if (copy) memcpy(copy, entries, count * sizeof(*copy));
    return copy;
}

#endif
