// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/**
 * @file openmc_export.h
 * @brief Export CSG geometry to OpenMC XML format
 */

#ifndef OPENMC_EXPORT_H
#define OPENMC_EXPORT_H

#include "core/alea_export.h"
#include "core/alea_system.h"


/**
 * @brief Export CSG system to OpenMC model.xml format
 *
 * Generates a complete OpenMC model.xml file containing:
 * - <materials> section: all materials with normalized fractions
 * - <geometry> section: surfaces and cells with region expressions
 * - <settings> section: basic simulation settings
 *
 * Material fractions are automatically normalized to sum to 1.0.
 * ZAID numbers are converted to OpenMC nuclide names (e.g., 92235 -> "U235").
 *
 * Note: Macrobodies (BOX, RCC) must be decomposed before export,
 * as OpenMC does not support macrobody surfaces.
 *
 * @param sys CSG system to export
 * @param ctx Export context with output stream
 * @return 0 on success, -1 on error
 */
int export_openmc(const alea_system_t* sys, export_context_t* ctx);


#endif /* OPENMC_EXPORT_H */
