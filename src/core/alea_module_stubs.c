// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/**
 * @file alea_module_stubs.c
 * @brief Weak stubs for optional module functions
 *
 * These stubs allow the core library to compile and run without
 * the optional modules (raycast, slice, MCNP, OpenMC). When a
 * module is linked, its implementation overrides these weak symbols.
 */

#include "alea_system.h"
#include <stdio.h>
#include "util/alea_log.h"

/* Forward declarations for types used in stubs */
struct alea_bvh;
struct export_context;

/* Raycast module stubs */

__attribute__((weak))
void alea_bvh_free(struct alea_bvh* bvh) {
    (void)bvh;  /* No-op when raycast module not linked */
}

/* MCNP module stubs */

__attribute__((weak))
int export_mcnp(const alea_system_t* sys, struct export_context* ctx) {
    (void)sys; (void)ctx;
    ALEA_LOG_ERROR("MCNP module not linked");
    return -1;
}

__attribute__((weak))
alea_system_t* mcnp_convert_file(const char* filename) {
    (void)filename;
    ALEA_LOG_ERROR("MCNP module not linked");
    return NULL;
}

/* Raycast volume estimation stubs */

__attribute__((weak))
int alea_estimate_cell_volumes(const alea_system_t* sys,
                               double ox, double oy, double oz,
                               double radius, int n_rays,
                               double* volumes, double* rel_errors) {
    (void)sys; (void)ox; (void)oy; (void)oz;
    (void)radius; (void)n_rays; (void)volumes; (void)rel_errors;
    ALEA_LOG_ERROR("Raycast module not linked");
    return -1;
}

__attribute__((weak))
int alea_estimate_instance_volumes(const alea_system_t* sys,
                                   int n_rays,
                                   double* volumes, double* rel_errors) {
    (void)sys; (void)n_rays; (void)volumes; (void)rel_errors;
    ALEA_LOG_ERROR("Raycast module not linked");
    return -1;
}

__attribute__((weak))
int alea_remove_cells_by_volume(alea_system_t* sys,
                                const double* volumes,
                                double threshold) {
    (void)sys; (void)volumes; (void)threshold;
    ALEA_LOG_ERROR("Raycast module not linked");
    return -1;
}

/* OpenMC module stubs */

__attribute__((weak))
int export_openmc(const alea_system_t* sys, struct export_context* ctx) {
    (void)sys; (void)ctx;
    ALEA_LOG_ERROR("OpenMC module not linked");
    return -1;
}

__attribute__((weak))
alea_system_t* openmc_convert_file(const char* filename) {
    (void)filename;
    ALEA_LOG_ERROR("OpenMC module not linked");
    return NULL;
}
