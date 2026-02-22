// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/**
 * @file alea_mesh.h
 * @brief Alea Mesh Export Module API
 *
 * Samples CSG geometry onto a structured hexahedral grid and exports
 * to Gmsh (.msh v2.2) or VTK (.vtk legacy) format for use with
 * mesh-based transport codes.
 *
 * Example workflow:
 * @code
 *   alea_system_t* sys = alea_load_mcnp("model.i");
 *   alea_build_universe_index(sys);
 *
 *   alea_mesh_config_t cfg;
 *   alea_mesh_config_init(&cfg);
 *   cfg.nx = cfg.ny = cfg.nz = 100;
 *
 *   alea_mesh_export_system(sys, &cfg, "output.msh");
 *   alea_destroy(sys);
 * @endcode
 */

#ifndef ALEA_MESH_H
#define ALEA_MESH_H

#include "alea.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ============================================================================
 * TYPES
 * ============================================================================ */

/** Output mesh format */
typedef enum {
    ALEA_MESH_GMSH = 0,   /**< Gmsh .msh v2.2 ASCII */
    ALEA_MESH_VTK          /**< VTK legacy .vtk ASCII */
} alea_mesh_format_t;

/** Mesh sampling configuration */
typedef struct {
    double x_min, x_max;        /**< X bounds (0,0 = auto-detect) */
    double y_min, y_max;        /**< Y bounds (0,0 = auto-detect) */
    double z_min, z_max;        /**< Z bounds (0,0 = auto-detect) */
    int nx, ny, nz;             /**< Elements per axis */
    const double *x_nodes;      /**< Custom X node positions (nx+1), NULL=uniform */
    const double *y_nodes;      /**< Custom Y node positions (ny+1), NULL=uniform */
    const double *z_nodes;      /**< Custom Z node positions (nz+1), NULL=uniform */
    alea_mesh_format_t format;
    int void_material_id;       /**< Material ID for void regions (default 0) */
    double auto_pad;            /**< Fractional padding for auto-bounds (default 0.01) */
} alea_mesh_config_t;

/** Mesh sampling result */
typedef struct {
    int nx, ny, nz;
    double *x_nodes, *y_nodes, *z_nodes;   /**< (n+1) positions each */
    int *material_ids;                      /**< nx*ny*nz, Z-major [k][j][i] */
    int *cell_ids;                          /**< nx*ny*nz, Z-major [k][j][i] */
    int num_materials;
    int *unique_materials;                  /**< Sorted array of unique material IDs */
} alea_mesh_result_t;

/* ============================================================================
 * API FUNCTIONS
 * ============================================================================ */

/**
 * @brief Initialize mesh config with defaults
 *
 * Sets nx=ny=nz=10, format=GMSH, void_material_id=0, auto_pad=0.01,
 * bounds and custom nodes to zero/NULL.
 */
void alea_mesh_config_init(alea_mesh_config_t *cfg);

/**
 * @brief Sample CSG geometry onto a structured grid
 *
 * Queries material at each cell center. Auto-detects bounds if all zero.
 * Builds spatial and universe indices if needed.
 *
 * @param sys CSG system (must have cells loaded)
 * @param cfg Sampling configuration
 * @return Allocated result (caller must free with alea_mesh_result_free), or NULL on error
 */
alea_mesh_result_t *alea_mesh_sample(const alea_system_t *sys,
                                             const alea_mesh_config_t *cfg);

/**
 * @brief Export mesh result to a file
 * @return 0 on success, -1 on error
 */
int alea_mesh_export(const alea_mesh_result_t *mesh,
                         alea_mesh_format_t fmt, const char *filename);

/**
 * @brief Export mesh result to an open stream
 * @return 0 on success, -1 on error
 */
int alea_mesh_export_stream(const alea_mesh_result_t *mesh,
                                alea_mesh_format_t fmt, FILE *out);

/**
 * @brief One-shot: sample CSG geometry and export to file
 * @return 0 on success, -1 on error
 */
int alea_mesh_export_system(const alea_system_t *sys,
                             const alea_mesh_config_t *cfg,
                             const char *filename);

/**
 * @brief Free mesh result allocated by alea_mesh_sample
 */
void alea_mesh_result_free(alea_mesh_result_t *mesh);

#ifdef __cplusplus
}
#endif

#endif /* ALEA_MESH_H */
