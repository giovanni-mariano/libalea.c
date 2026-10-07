// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

#ifndef MCNP_MODEL_H
#define MCNP_MODEL_H

/**
 * @file mcnp_model.h
 * @brief Internal MCNP model header
 *
 * Re-exports the public API from alea_mcnp.h and adds internal declarations.
 */

#include "alea_mcnp.h"
#include "util/cell_parameter.h"

static inline void mcnp_scoped_params_free(mcnp_scoped_cell_params_t* params) {
    free(params->detector_probabilities);
    memset(params, 0, sizeof(*params));
}

static inline int mcnp_scoped_params_copy(mcnp_scoped_cell_params_t* dst,
                                         const mcnp_scoped_cell_params_t* src) {
    if (dst == src) return 0;
    alea_detector_probability_t* copy = alea_detector_probability_copy(
        src->detector_probabilities, src->detector_probability_count);
    if (src->detector_probability_count && !copy) return -1;
    mcnp_scoped_params_free(dst);
    *dst = *src;
    dst->detector_probabilities = copy;
    return 0;
}

/* Merge only explicit entries, leaving defaults and other scopes untouched. */
static inline int mcnp_scoped_params_merge(mcnp_scoped_cell_params_t* dst,
                                          const mcnp_scoped_cell_params_t* src,
                                          int override) {
    for (int i = 0; i < ALEA_PARTICLE_COUNT; i++) {
        uint32_t bit = 1u << i;
        if ((src->energy_cutoff_particles & bit) &&
            (override || !(dst->energy_cutoff_particles & bit))) {
            dst->energy_cutoff[i] = src->energy_cutoff[i];
            dst->energy_cutoff_particles |= bit;
        }
        if ((src->secondary_state_particles & bit) &&
            (override || !(dst->secondary_state_particles & bit))) {
            dst->secondary_state[i] = src->secondary_state[i];
            dst->secondary_state_particles |= bit;
        }
    }
    for (size_t i = 0; i < src->detector_probability_count; i++) {
        const alea_detector_probability_t* entry = &src->detector_probabilities[i];
        int present = 0;
        for (size_t j = 0; j < dst->detector_probability_count; j++)
            if (dst->detector_probabilities[j].tally == entry->tally) present = 1;
        if ((override || !present) && alea_detector_probability_store(
            &dst->detector_probabilities, &dst->detector_probability_count,
            entry->tally, entry->probability)) return -1;
    }
    return 0;
}

/* Internal storage and lifecycle helpers. */
int mcnp_model_reserve_params(mcnp_model_t* model, size_t cap);
int mcnp_model_add_params(mcnp_model_t* model);
void mcnp_model_register_hooks(mcnp_model_t* model);
uint32_t mcnp_model_add_inline_transform(mcnp_model_t* model,
                                         const double* values,
                                         int count,
                                         int degrees);

/**
 * @brief Internal: Full MCNP file conversion returning model
 *
 * Called by mcnp_load(). Parses the file, creates the system and model,
 * converts all geometry/materials/transforms, and populates cell params.
 */
mcnp_model_t* mcnp_convert_to_model(const char* filename);

/**
 * @brief Internal: Full MCNP buffer conversion returning model
 */
mcnp_model_t* mcnp_convert_buffer_to_model(const char* input, size_t len,
                                           const char* source_name);


#endif /* MCNP_MODEL_H */
