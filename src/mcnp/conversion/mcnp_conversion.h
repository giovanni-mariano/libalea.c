// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

#ifndef MCNP_CONVERSION_H
#define MCNP_CONVERSION_H

#include "core/alea_system.h"
#include "mcnp/mcnp_model.h"
#include <stddef.h>

// ============================================================================
// CONVERSION FUNCTIONS
// ============================================================================

// Convert an MCNP file to mcnp_model_t (owns system + cell params)
mcnp_model_t* mcnp_convert_to_model(const char* filename);

// Convert an MCNP input buffer to mcnp_model_t (owns system + cell params)
mcnp_model_t* mcnp_convert_buffer_to_model(const char* input, size_t len,
                                           const char* source_name);

#endif // MCNP_CONVERSION_H
