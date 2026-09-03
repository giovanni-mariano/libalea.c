// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/**
 * @file alea_mesh.h
 * @brief Alea Mesh Export Module API
 *
 * Samples CSG geometry onto a structured rectilinear grid and exports
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

/** Optional cell data written by mesh exporters. */
typedef enum {
    ALEA_MESH_EXPORT_MIXED_FLAG = 1u << 0,
    ALEA_MESH_EXPORT_DOMINANT_FRACTION = 1u << 1,
    ALEA_MESH_EXPORT_TIE_FLAG = 1u << 2,
    ALEA_MESH_EXPORT_SAMPLE_COUNT = 1u << 3,
    ALEA_MESH_EXPORT_MATERIAL_FRACTIONS = 1u << 4,
    ALEA_MESH_EXPORT_ESTIMATED_ERROR = 1u << 5,
    ALEA_MESH_EXPORT_REFINEMENT_FLAG = 1u << 6
} alea_mesh_export_field_t;

typedef struct {
    uint32_t fields;              /**< ALEA_MESH_EXPORT_* mask */
    int max_fraction_materials;   /**< Safety limit for dense exported arrays */
} alea_mesh_export_options_t;

/** Voxel composition sampling mode */
typedef enum {
    ALEA_MESH_SAMPLE_CENTER = 0,   /**< Sample only the voxel center */
    ALEA_MESH_SAMPLE_CORNERS,      /**< Probe near voxel corners */
    ALEA_MESH_SAMPLE_SUBCELL,      /**< Probe an NxNxN subcell lattice */
    ALEA_MESH_SAMPLE_STRATIFIED,   /**< Deterministic jitter within subcells */
    ALEA_MESH_SAMPLE_ADAPTIVE      /**< Refine stratified estimates to tolerance */
} alea_mesh_sampling_mode_t;

typedef enum {
    ALEA_MESH_BOUNDS_LEGACY = 0, /**< All-zero bounds mean auto; otherwise explicit */
    ALEA_MESH_BOUNDS_AUTO,
    ALEA_MESH_BOUNDS_EXPLICIT
} alea_mesh_bounds_mode_t;

typedef enum {
    ALEA_MESH_BOUNDS_SOURCE_EXPLICIT = 0,
    ALEA_MESH_BOUNDS_SOURCE_CUSTOM_NODES,
    ALEA_MESH_BOUNDS_SOURCE_INFERRED_ROOT_AABB
} alea_mesh_bounds_source_t;

/** Deterministic dominant-selection diagnostics for one voxel. */
typedef enum {
    ALEA_MESH_TIE_NONE = 0,
    ALEA_MESH_TIE_MATERIAL = 1u << 0, /**< Multiple materials share max count */
    ALEA_MESH_TIE_CELL = 1u << 1      /**< Multiple cells share max count */
} alea_mesh_tie_flag_t;

/** Optional arrays materialized by alea_mesh_sample(). */
typedef enum {
    ALEA_MESH_FIELD_MATERIAL_ID = 1u << 0,
    ALEA_MESH_FIELD_CELL_ID = 1u << 1,
    ALEA_MESH_FIELD_MIXED_FLAG = 1u << 2,
    ALEA_MESH_FIELD_DOMINANT_FRACTION = 1u << 3,
    ALEA_MESH_FIELD_SAMPLED_FRACTIONS = 1u << 4,
    ALEA_MESH_FIELD_SAMPLE_COUNT = 1u << 5,
    ALEA_MESH_FIELD_TIE_FLAG = 1u << 6,
    ALEA_MESH_FIELD_ESTIMATED_ERROR = 1u << 7,
    ALEA_MESH_FIELD_REFINEMENT_FLAG = 1u << 8,
    ALEA_MESH_FIELD_CELL_FRACTIONS = 1u << 9
} alea_mesh_result_field_t;

#define ALEA_MESH_REFINEMENT_LIMIT_REACHED 0x01u

typedef int (*alea_mesh_progress_fn)(size_t completed_voxels,
                                     size_t total_voxels,
                                     void *user_data);

typedef struct alea_mesh_material_fraction alea_mesh_material_fraction_t;
typedef struct alea_mesh_cell_fraction alea_mesh_cell_fraction_t;

typedef struct {
    int i, j, k;
    double x_min, x_max, y_min, y_max, z_min, z_max;
    int material_id;
    int cell_id;
    unsigned char mixed;
    uint8_t tie_flags;
    double dominant_fraction;
    double estimated_error;
    uint32_t sample_count;
    uint8_t refinement_flags;
    const alea_mesh_material_fraction_t *fractions; /**< Valid during callback */
    uint32_t fraction_count;
    const alea_mesh_cell_fraction_t *cell_fractions; /**< Valid during callback */
    uint32_t cell_fraction_count;
} alea_mesh_voxel_sample_t;

typedef int (*alea_mesh_voxel_visit_fn)(const alea_mesh_voxel_sample_t *sample,
                                        void *user_data);

/** One sampled material-fraction estimate for a voxel (not exact volume). */
struct alea_mesh_material_fraction {
    int material_id;
    double fraction;              /**< Point-count fraction in [0,1] */
};

/** One sampled concrete-cell fraction estimate for a voxel. */
struct alea_mesh_cell_fraction {
    int cell_id;                  /**< MCNP/OpenMC cell ID; -1 for void */
    int material_id;              /**< Material associated with this owner */
    double fraction;              /**< Point-count fraction in [0,1] */
};

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
    double mixed_threshold;      /**< Max tolerated non-dominant fraction */
    double target_error;         /**< Adaptive empirical L1 tolerance */
    int max_refine_depth;
    uint32_t max_samples_per_voxel;
    uint64_t max_total_samples;  /**< 0 means unlimited */
    uint64_t sampling_seed;
    int workers;                 /**< 0=runtime default, 1=serial */
    alea_mesh_bounds_mode_t bounds_mode;
    uint32_t fields;             /**< ALEA_MESH_FIELD_* arrays to retain */
    alea_mesh_progress_fn progress; /**< Optional; nonzero return cancels */
    void *progress_user_data;
    alea_mesh_voxel_visit_fn visit; /**< Optional streaming voxel callback */
    void *visit_user_data;
} alea_mesh_config_t;

/** Mesh sampling result */
typedef struct {
    int nx, ny, nz;
    uint32_t fields;                         /**< Materialized ALEA_MESH_FIELD_* */
    alea_mesh_bounds_source_t bounds_source;
    double bounds_padding;                   /**< Fractional padding actually used */
    alea_mesh_sampling_mode_t sampling_mode;
    uint64_t sampling_seed;
    double target_error;                     /**< Requested empirical tolerance */
    double *x_nodes, *y_nodes, *z_nodes;   /**< (n+1) positions each */
    int *material_ids;                      /**< Dominant sampled material; X fastest */
    int *cell_ids;                          /**< Dominant cell in dominant material */
    int num_materials;
    int *unique_materials;                  /**< Sorted array of unique material IDs */
    unsigned char *mixed_flags;             /**< nx*ny*nz flags, 1=mixed */
    double *dominant_fractions;             /**< nx*ny*nz dominant sampled fraction */
    double *estimated_errors;                /**< Empirical refinement differences */
    uint32_t *sample_counts;                 /**< nx*ny*nz point sample counts */
    uint8_t *tie_flags;                      /**< nx*ny*nz ALEA_MESH_TIE_* bits */
    uint8_t *refinement_flags;               /**< ALEA_MESH_REFINEMENT_* bits */
    int mixed_count;                        /**< Number of voxels flagged as mixed */
    alea_mesh_fraction_span_t *fraction_spans;
                                            /**< nx*ny*nz spans into fractions */
    alea_mesh_material_fraction_t *fractions;
                                            /**< Packed material fractions */
    size_t fraction_count;                  /**< Number of packed fraction entries */
    alea_mesh_fraction_span_t *cell_fraction_spans;
                                            /**< nx*ny*nz spans into cell_fractions */
    alea_mesh_cell_fraction_t *cell_fractions;
                                            /**< Packed concrete-cell fractions */
    size_t cell_fraction_count;             /**< Number of packed cell entries */
} alea_mesh_result_t;

/** One node in a nonconforming adaptive hexahedral grid. IDs are 1-based. */
typedef struct {
    uint64_t id;
    uint64_t parent_id;                    /**< 0 for a root voxel */
    uint64_t child_ids[8];                 /**< 0 for a leaf voxel */
    uint32_t level;
    unsigned char is_leaf;
    uint8_t flags;                         /**< ALEA_ADAPTIVE_GRID_* bits */
    double x_min, x_max, y_min, y_max, z_min, z_max;
    int material_id;
    int cell_id;
    unsigned char mixed;
    uint8_t tie_flags;
    uint8_t refinement_flags;
    double dominant_fraction;
    double estimated_error;
    uint32_t sample_count;
} alea_adaptive_grid_cell_t;

#define ALEA_ADAPTIVE_GRID_DEPTH_LIMIT_REACHED 0x01u
#define ALEA_ADAPTIVE_GRID_CELL_LIMIT_REACHED  0x02u

typedef struct {
    alea_mesh_config_t sampling;            /**< Initial grid and estimator */
    uint32_t max_grid_depth;                /**< Octree levels below roots */
    size_t max_cells;                       /**< Includes internal nodes */
    int refine_mixed;
    int refine_high_error;
} alea_adaptive_grid_config_t;

typedef struct {
    alea_adaptive_grid_cell_t *cells;       /**< Stable append-order IDs */
    size_t cell_count;                      /**< Internal plus leaf cells */
    size_t leaf_count;
    size_t root_count;
    uint32_t max_level;
    int balanced;                           /**< Currently 0: nonconforming octree */
} alea_adaptive_grid_result_t;

/* ============================================================================
 * API FUNCTIONS
 * ============================================================================ */

/**
 * @brief Initialize mesh config with defaults
 *
 * Sets nx=ny=nz=10, format=GMSH, void_material_id=0, auto_pad=0.01,
 * sampling_mode=SUBCELL, subsamples_per_axis=2, mixed_threshold=0,
 * target_error=0.05, max_refine_depth=3, max_samples_per_voxel=32768,
 * workers=1, all current result fields enabled, bounds_mode=LEGACY, and bounds/custom
 * nodes/callbacks to zero/NULL.
 */
void alea_mesh_config_init(alea_mesh_config_t *cfg);

/** Initialize exporter options with diagnostic fields and no dense fractions. */
void alea_mesh_export_options_init(alea_mesh_export_options_t *options);

/**
 * @brief Sample CSG geometry onto a structured grid
 *
 * Samples voxel composition, assigns the dominant material to material_ids,
 * and records a cell observed with that material. Ties use the lowest material
 * ID and then the lowest cell ID, and are recorded in tie_flags. Reported
 * fractions are point-count estimates, not exact material volumes. Auto-detects
 * bounds if all zero.
 * Builds spatial and universe indices if needed.
 *
 * @param sys CSG system (must have cells loaded)
 * @param cfg Sampling configuration
 * @return Allocated result (caller must free with alea_mesh_result_free), or NULL on error
 */
alea_mesh_result_t *alea_mesh_sample(alea_system_t *sys,
                                             const alea_mesh_config_t *cfg);

/** Visit voxels without retaining per-voxel result arrays. */
int alea_mesh_visit(alea_system_t *sys, const alea_mesh_config_t *cfg,
                    alea_mesh_voxel_visit_fn visit, void *user_data);

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

/** Export to an open stream with explicit cell-data options. */
int alea_mesh_export_stream_ex(const alea_mesh_result_t *mesh,
                               alea_mesh_format_t fmt, FILE *out,
                               const alea_mesh_export_options_t *options);

/** Export to a file with explicit cell-data options. */
int alea_mesh_export_ex(const alea_mesh_result_t *mesh,
                        alea_mesh_format_t fmt, const char *filename,
                        const alea_mesh_export_options_t *options);

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

/** Initialize adaptive-grid defaults (adaptive composition, depth 4). */
void alea_adaptive_grid_config_init(alea_adaptive_grid_config_t *cfg);

/** Build a separate, nonconforming octree grid by refining mixed/high-error cells. */
alea_adaptive_grid_result_t *alea_adaptive_grid_sample(
    alea_system_t *sys, const alea_adaptive_grid_config_t *cfg);

void alea_adaptive_grid_result_free(alea_adaptive_grid_result_t *grid);

/** Export leaf cells as an unstructured Gmsh 2.2 or VTK legacy hex grid. */
int alea_adaptive_grid_export_stream(const alea_adaptive_grid_result_t *grid,
                                     alea_mesh_format_t fmt, FILE *out);
int alea_adaptive_grid_export(const alea_adaptive_grid_result_t *grid,
                              alea_mesh_format_t fmt, const char *filename);

#ifdef __cplusplus
}
#endif

#endif /* ALEA_MESH_H */
