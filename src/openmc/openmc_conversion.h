// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/**
 * @file openmc_conversion.h
 * @brief Convert OpenMC XML geometry to CSG system
 */

#ifndef OPENMC_CONVERSION_H
#define OPENMC_CONVERSION_H

#include "core/alea_system.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ============================================================================
 * CONVERSION FUNCTIONS
 * ============================================================================ */

/**
 * @brief Convert OpenMC XML file to CSG system
 *
 * Parses an OpenMC geometry XML file (typically geometry.xml or model.xml)
 * and creates a CSG system with equivalent geometry.
 *
 * @param filename Path to XML file
 * @return New CSG system or NULL on error
 */
alea_system_t* openmc_convert_file(const char* filename);

/**
 * @brief Convert OpenMC XML string to CSG system
 */
alea_system_t* openmc_convert_string(const char* xml_content, size_t length);

#ifdef __cplusplus
}
#endif

#endif /* OPENMC_CONVERSION_H */
