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

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief Internal: Full MCNP file conversion returning model
 *
 * Called by mcnp_load(). Parses the file, creates the system and model,
 * converts all geometry/materials/transforms, and populates cell params.
 */
mcnp_model_t* mcnp_convert_to_model(const char* filename);

#ifdef __cplusplus
}
#endif

#endif /* MCNP_MODEL_H */
