// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/**
 * @file alea_serpent.h
 * @brief Public API for Serpent export
 *
 * Users who want to export Serpent input files should include this header
 * and link against libalea_serpent.a.
 */

#ifndef ALEA_SERPENT_H
#define ALEA_SERPENT_H

#include "alea.h"
#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief Export a system to a Serpent input file
 * @return 0 on success, -1 on error
 */
int serpent_export_system(alea_system_t* sys, const char* filename);

/**
 * @brief Export a system to a Serpent input stream
 * @return 0 on success, -1 on error
 */
int serpent_export_system_stream(alea_system_t* sys, FILE* out);

#ifdef __cplusplus
}
#endif

#endif /* ALEA_SERPENT_H */
