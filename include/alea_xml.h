// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/** @file alea_xml.h @brief Native, lossless ALEA model XML I/O. */

#ifndef ALEA_XML_H
#define ALEA_XML_H

#include "alea_model.h"
#include <stddef.h>
#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

/** Load an ALEA XML model from a file. The returned model owns its system. */
alea_model_t* alea_xml_load(const char* filename);

/** Load an ALEA XML model from memory (length 0 means NUL-terminated). */
alea_model_t* alea_xml_load_string(const char* xml, size_t length);

/** Write a complete ALEA model to a file. */
int alea_xml_export(const alea_model_t* model, const char* filename);

/** Write a complete ALEA model to an already-open stream. */
int alea_xml_export_stream(const alea_model_t* model, FILE* stream);

/** Convenience export for a geometry system without model metadata. */
int alea_xml_export_system(const alea_system_t* sys, const char* filename);

/** Convenience stream export for a geometry system without model metadata. */
int alea_xml_export_system_stream(const alea_system_t* sys, FILE* stream);

#ifdef __cplusplus
}
#endif

#endif
