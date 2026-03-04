// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

#ifndef OPENMC_MODEL_H
#define OPENMC_MODEL_H

/**
 * @file openmc_model.h
 * @brief Internal OpenMC model header
 *
 * Re-exports the public API from alea_openmc.h and adds internal declarations.
 */

#include "alea_openmc.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief Internal: Full OpenMC file conversion returning model
 *
 * Called by openmc_load(). Parses the file, creates the system and model.
 */
openmc_model_t* openmc_convert_to_model(const char* filename);

/**
 * @brief Internal: OpenMC string conversion returning model
 *
 * Called by openmc_load_string().
 */
openmc_model_t* openmc_convert_string_to_model(const char* input, size_t len);

#ifdef __cplusplus
}
#endif

#endif /* OPENMC_MODEL_H */
