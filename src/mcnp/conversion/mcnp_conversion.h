// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

#ifndef MCNP_CONVERSION_H
#define MCNP_CONVERSION_H

#include "core/alea_system.h"

// ============================================================================
// CONVERSION FUNCTIONS
// ============================================================================

// Convert an MCNP file to CSG system
alea_system_t* mcnp_convert_file(const char* filename);

#endif // MCNP_CONVERSION_H
