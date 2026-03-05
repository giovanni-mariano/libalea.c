// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/**
 * @file mcnp_export.h
 * @brief MCNP format export functions
 *
 * MCNP-specific export functionality including CSG tree to expression
 * conversion, cell card writing, and surface card writing.
 */

#ifndef MCNP_EXPORT_H
#define MCNP_EXPORT_H

#include "alea_types.h"
#include "mcnp_str.h"
#include <stdio.h>


/* Forward declarations */
typedef struct export_context export_context_t;

/**
 * @brief Core MCNP export function
 *
 * Exports cells, surfaces, materials, and transforms in MCNP format.
 *
 * @param sys CSG system
 * @param ctx Export context
 * @return 0 on success, -1 on error
 */
int export_mcnp(const alea_system_t* sys, export_context_t* ctx);


#endif /* MCNP_EXPORT_H */
