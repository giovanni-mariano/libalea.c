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
