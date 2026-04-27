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
 *   mcnp_model_t* model = mcnp_load("model.i");
 *   alea_system_t* sys = model->sys;
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

/** Voxel composition sampling mode */
typedef enum {
    ALEA_MESH_SAMPLE_CENTER = 0,   /**< Sample only the voxel center */
    ALEA_MESH_SAMPLE_CORNERS,      /**< Probe near voxel corners */
    ALEA_MESH_SAMPLE_SUBCELL       /**< Probe an NxNxN subcell lattice */
} alea_mesh_sampling_mode_t;

/** One material fraction entry for a sampled voxel */
typedef struct {
    int material_id;
    double fraction;
} alea_mesh_material_fraction_t;

/** Span into alea_mesh_result_t::fractions for one voxel */
typedef struct {
    uint32_t offset;
    uint32_t count;
} alea_mesh_fraction_span_t;

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
    alea_mesh_sampling_mode_t sampling_mode;
                                  /**< Composition sampling mode */
    int subsamples_per_axis;     /**< Subsample lattice size for SAMPLE_SUBCELL */
    double mixed_threshold;      /**< Tolerance for declaring a voxel mixed */
} alea_mesh_config_t;

/** Mesh sampling result */
typedef struct {
    int nx, ny, nz;
    double *x_nodes, *y_nodes, *z_nodes;   /**< (n+1) positions each */
    int *material_ids;                      /**< Dominant sampled material per voxel */
    int *cell_ids;                          /**< Center-sampled cell ID per voxel */
    int num_materials;
    int *unique_materials;                  /**< Sorted array of unique material IDs */
    unsigned char *mixed_flags;             /**< nx*ny*nz flags, 1=mixed */
    double *dominant_fractions;             /**< nx*ny*nz dominant sampled fraction */
    int mixed_count;                        /**< Number of voxels flagged as mixed */
    alea_mesh_fraction_span_t *fraction_spans;
                                            /**< nx*ny*nz spans into fractions */
    alea_mesh_material_fraction_t *fractions;
                                            /**< Packed material fractions */
    size_t fraction_count;                  /**< Number of packed fraction entries */
} alea_mesh_result_t;

/* ============================================================================
 * API FUNCTIONS
 * ============================================================================ */

/**
 * @brief Initialize mesh config with defaults
 *
 * Sets nx=ny=nz=10, format=GMSH, void_material_id=0, auto_pad=0.01,
 * sampling_mode=SUBCELL, subsamples_per_axis=2, mixed_threshold=0,
 * bounds and custom nodes to zero/NULL.
 */
void alea_mesh_config_init(alea_mesh_config_t *cfg);

/**
 * @brief Sample CSG geometry onto a structured grid
 *
 * Samples voxel composition, assigns the dominant material to material_ids,
 * and records center-sampled cell IDs. Auto-detects bounds if all zero.
 * Builds spatial and universe indices if needed.
 *
 * @param sys CSG system (must have cells loaded)
 * @param cfg Sampling configuration
 * @return Allocated result (caller must free with alea_mesh_result_free), or NULL on error
 */
alea_mesh_result_t *alea_mesh_sample(alea_system_t *sys,
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
int alea_mesh_export_system(alea_system_t *sys,
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
