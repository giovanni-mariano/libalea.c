// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/**
 * @file openmc_conversion.h
 * @brief Convert OpenMC XML geometry to CSG system
 */

#ifndef OPENMC_CONVERSION_H
#define OPENMC_CONVERSION_H

#include "openmc_model.h"


/* ============================================================================
 * CONVERSION FUNCTIONS
 * ============================================================================ */

/**
 * @brief Convert OpenMC XML file to model
 *
 * Parses an OpenMC geometry XML file (typically geometry.xml or model.xml)
 * and creates a model with equivalent geometry.
 *
 * @param filename Path to XML file
 * @return New model or NULL on error
 */
openmc_model_t* openmc_convert_to_model(const char* filename);

/**
 * @brief Convert OpenMC XML string to model
 */
openmc_model_t* openmc_convert_string_to_model(const char* xml_content, size_t length);


#endif /* OPENMC_CONVERSION_H */
