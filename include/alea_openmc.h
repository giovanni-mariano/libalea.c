// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/**
 * @file alea_openmc.h
 * @brief Public API for OpenMC module
 *
 * Users who want to load/export OpenMC files should include this header
 * and link against libalea_openmc.a.
 *
 * Basic usage:
 *   #include "alea.h"
 *   #include "alea_openmc.h"
 *
 *   openmc_model_t* model = openmc_load("model.xml");
 *   alea_system_t* sys = model->sys;
 *   // ... work with sys (geometry, queries, void generation) ...
 *   openmc_export(model, "output.xml");
 *   openmc_model_destroy(model);
 */

#ifndef ALEA_OPENMC_H
#define ALEA_OPENMC_H

#include "alea.h"
#include <stddef.h>
#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

/* ============================================================================
 * OPENMC MODEL
 * ============================================================================ */

/**
 * @brief OpenMC model: geometry system wrapper
 */
typedef struct {
    alea_system_t* sys;
    int owns_sys;               /* 1 if destroy should free sys */
} openmc_model_t;

/* ============================================================================
 * LOADING
 * ============================================================================ */

/**
 * @brief Load an OpenMC XML file into an openmc_model_t
 * @param filename Path to OpenMC XML file (geometry.xml or model.xml)
 * @return Model on success, NULL on error
 */
openmc_model_t* openmc_load(const char* filename);

/**
 * @brief Load OpenMC from a string buffer
 * @param input OpenMC XML text
 * @param len Length of input (0 = use strlen)
 * @return Model on success, NULL on error
 */
openmc_model_t* openmc_load_string(const char* input, size_t len);

/* ============================================================================
 * EXPORT
 * ============================================================================ */

/**
 * @brief Export OpenMC model to file
 * @return 0 on success, -1 on error
 */
int openmc_export(const openmc_model_t* model, const char* filename);

/**
 * @brief Export OpenMC model to stream
 * @return 0 on success, -1 on error
 */
int openmc_export_stream(const openmc_model_t* model, FILE* out);

/**
 * @brief Export a bare system to OpenMC file (convenience, uses default config)
 * @return 0 on success, -1 on error
 */
int openmc_export_system(alea_system_t* sys, const char* filename);

/**
 * @brief Export a bare system to OpenMC stream (convenience, uses default config)
 * @return 0 on success, -1 on error
 */
int openmc_export_system_stream(alea_system_t* sys, FILE* out);

/* ============================================================================
 * MODEL MANAGEMENT
 * ============================================================================ */

/**
 * @brief Destroy an OpenMC model
 *
 * If owns_sys is set, also destroys the system.
 */
void openmc_model_destroy(openmc_model_t* model);

/**
 * @brief Create a non-owning model wrapper around an existing system
 *
 * The returned model does NOT own the system — openmc_model_destroy()
 * will free the model but not the system.
 *
 * @param sys System to wrap (caller retains ownership)
 * @return Model on success, NULL on failure
 */
openmc_model_t* openmc_model_wrap(alea_system_t* sys);

#ifdef __cplusplus
}
#endif

#endif /* ALEA_OPENMC_H */
