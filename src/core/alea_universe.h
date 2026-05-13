// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

#ifndef ALEA_UNIVERSE_H
#define ALEA_UNIVERSE_H

#include "alea_types.h"
#include "util/alea_vec.h"
#include <stdbool.h>


/**
 * @file alea_universe.h
 * @brief Universe types, point queries, and flattening
 *
 * Key concepts:
 * - Universe: a collection of cells (defined by U= parameter)
 * - Point queries: recursive descent through FILL hierarchy
 * - Flattening: expand all FILL references into a single universe
 */

/* Forward declarations */
typedef struct alea_system alea_system_t;
typedef struct alea_material alea_material_t;

// ============================================================================
// UNIVERSE TYPE
// ============================================================================

/**
 * @brief Universe information
 */
typedef struct {
    alea_bbox_t bbox;
    uint32_t left_or_first;
    uint32_t right_child;
    uint16_t count;
    uint8_t axis;
    uint8_t pad;
} alea_universe_bvh_node_t;

ALEA_VEC_DEFINE(alea_universe_bvh_node_vec, alea_universe_bvh_node_t);

typedef struct {
    int universe_id;
    alea_size_vec_t cell_indices;
    alea_bbox_t bbox;         // Bounding box of all cells
    alea_universe_bvh_node_vec_t point_bvh_nodes;
    uint32_t* point_bvh_indices;
    bool point_bvh_built;
    bool point_bvh_disabled;
} alea_universe_t;

ALEA_VEC_DEFINE(alea_universe_vec, alea_universe_t);

typedef struct {
    size_t bvh_builds;
    size_t bvh_queries;
    size_t bvh_node_visits;
    size_t bvh_bbox_tests;
    size_t bvh_candidates;
    size_t linear_universe_scans;
    size_t linear_cell_tests;
    size_t exact_cell_tests;
} alea_universe_point_bvh_stats_t;

void alea_universe_point_bvh_stats_reset(void);
void alea_universe_point_bvh_stats_get(alea_universe_point_bvh_stats_t* out);

/**
 * @brief Hash table for tracking primitive ID remapping during flatten
 *
 * Maps: src_primitive_id -> dst_primitive_id
 * Needed to update surface references after primitives are copied.
 */
typedef struct {
    uint32_t* map;          /* Array: map[src_id] = dst_id */
    int8_t* inverted;       /* Array: cached inverted flags per src_id */
    size_t capacity;        /* Size of map/inverted arrays */
} primitive_remap_t;

/* ============================================================================
 * TRANSFORM OPERATIONS
 * ============================================================================ */

/**
 * @brief 4x4 homogeneous transform matrix (row-major)
 *
 * [ R00 R01 R02 Tx ]
 * [ R10 R11 R12 Ty ]
 * [ R20 R21 R22 Tz ]
 * [  0   0   0   1 ]
 *
 * Stored as 12 values: R00,R01,R02,Tx, R10,R11,R12,Ty, R20,R21,R22,Tz
 */
typedef struct {
    double m[12];       /* 3x4 matrix (bottom row implicit: 0,0,0,1) */
    double inv[12];     /* Inverse matrix */
    bool has_inverse;   /* True if inverse is computed */
} alea_matrix_t;

/**
 * @brief Initialize matrix from MCNP transform data
 */
bool alea_matrix_from_mcnp(alea_matrix_t* mat, const double* data,
                          int count, bool degrees);

/**
 * @brief Set matrix to identity
 */
void alea_matrix_identity(alea_matrix_t* mat);

/**
 * @brief Compose two matrices: result = a * b (b applied first)
 */
void alea_matrix_multiply(alea_matrix_t* result,
                         const alea_matrix_t* a,
                         const alea_matrix_t* b);

/**
 * @brief Compute inverse of matrix
 * @return true on success, false if singular
 */
bool alea_matrix_invert(alea_matrix_t* mat);

/**
 * @brief Transform a point: p' = M * p
 */
void alea_matrix_transform_point(const alea_matrix_t* mat,
                                double* x, double* y, double* z);

/**
 * @brief Inverse transform a point: p = M^-1 * p'
 */
void alea_matrix_transform_point_inverse(const alea_matrix_t* mat,
                                        double* x, double* y, double* z);

/* ============================================================================
 * POINT QUERY WITH LAZY EVALUATION
 * ============================================================================ */

/**
 * @brief Find cell containing point (lazy, handles FILL)
 *
 * Recursively descends into FILL cells by inverse-transforming the point.
 *
 * @param sys CSG system
 * @param x,y,z Point coordinates
 * @param out_cell_id Output: MCNP cell ID containing point
 * @param out_material Output: material ID (0 = void)
 * @param out_depth Output: nesting depth (0 = base universe)
 * @return 0 if found, -1 if void/error
 */
int alea_find_cell_lazy(const alea_system_t* sys,
                       double x, double y, double z,
                       int* out_cell_id,
                       int* out_material,
                       int* out_depth);

/* alea_cell_hit_t is defined in alea_types.h */

/**
 * @brief Find ALL cells containing a point at all nesting levels
 *
 * Traverses the universe hierarchy and returns every cell the point
 * passes through, from outermost (base universe) to innermost (terminal).
 *
 * @param sys CSG system
 * @param x,y,z Point in global coordinates
 * @param out_hits Array to fill with cell information
 * @param max_hits Maximum entries in out_hits
 * @return Number of cells found (hierarchy depth), or -1 on error
 */
int alea_find_all_cells_at_point(alea_system_t* sys,
                                double x, double y, double z,
                                alea_cell_hit_t* out_hits,
                                size_t max_hits);

/**
 * @brief Pre-build per-universe point BVHs for all universes.
 *
 * Called from query-cache preparation so subsequent concurrent point queries
 * (OpenMP grid render etc.) never enter the lazy BVH build path. Returns 0
 * on success; per-universe errors are recorded as disabled flags and do not
 * abort the prebuild loop.
 */
int alea_prebuild_universe_point_bvhs(alea_system_t* sys);

/**
 * @brief Find the deepest hierarchy-aware cell hit for a point
 *
 * Uses the same fill/lattice-aware semantics as `alea_find_all_cells_at_point()`
 * and returns the innermost hit in the ordered stack.
 *
 * @param sys CSG system
 * @param x,y,z Point in global coordinates
 * @param out_hit Receives the deepest hit on success
 * @return 0 on success, -1 if the point is in void or on error
 */
int alea_find_deepest_cell_hit_at_point(alea_system_t* sys,
                                       double x, double y, double z,
                                       alea_cell_hit_t* out_hit);

/* Debug version that always uses recursive path (for comparison) */
int alea_find_all_cells_at_point_recursive(const alea_system_t* sys,
                                          double x, double y, double z,
                                          alea_cell_hit_t* out_hits,
                                          size_t max_hits);

/**
 * @brief Enable/disable debug point trace output
 * @param enable Non-zero to enable, 0 to disable
 */
void alea_set_debug_point_trace(int enable);

/* ============================================================================
 * EXPLICIT FLATTENING
 * ============================================================================ */

/**
 * @brief Flatten configuration
 */
typedef struct alea_flatten_config {
    int starting_universe_id;   /* Universe to start flattening from (default: 0) */
    int max_depth;              /* Max recursion depth (0 = unlimited) */
    int starting_cell_id;       /* First cell ID in output (default: 1) */
    bool copy_materials;        /* Copy material definitions to new system */
    bool copy_transforms;       /* Copy transform definitions (usually not needed) */
} alea_flatten_config_t;

/**
 * @brief Default flatten configuration
 */
extern const alea_flatten_config_t ALEA_FLATTEN_DEFAULT;

/**
 * @brief Flatten in-place (destructive)
 *
 * Expands all FILL references in the system, then removes:
 * - All cells with FILL (replaced by their expanded content)
 * - All non-zero universes (content moved to universe 0)
 *
 * After this call, the system contains only terminal cells in universe 0.
 *
 * @param sys System to flatten (modified in place)
 * @param config Configuration (NULL = defaults)
 * @return Number of terminal cells, or -1 on error
 */
int alea_flatten_in_place(alea_system_t* sys,
                         const alea_flatten_config_t* config);

/**
 * @brief Clone a CSG tree with transformed primitives (same system)
 *
 * @param sys CSG system
 * @param root_id Root of tree to clone
 * @param mat Transform matrix
 * @param prim_remap Optional: primitive remapping cache (can be NULL)
 * @param prim_inverted Optional: inverted flags from canonicalization (can be NULL)
 * @param remap_size Size of remap arrays
 * @return New root node ID, or ALEA_NODE_ID_INVALID on error
 */
alea_node_id_t alea_clone_tree_transformed(alea_system_t* sys,
                                         alea_node_id_t root_id,
                                         const alea_matrix_t* mat,
                                         uint32_t* prim_remap,
                                         int8_t* prim_inverted,
                                         size_t remap_size);

/* ============================================================================
 * REMAP TABLE AND COPY HELPERS
 * ============================================================================ */

/**
 * @brief Create a primitive remap table
 */
primitive_remap_t* alea_create_remap_table(size_t src_prim_count);

/**
 * @brief Destroy a primitive remap table
 */
void alea_destroy_remap_table(primitive_remap_t* remap);

/**
 * @brief Copy surfaces from src to dst, remapping primitive IDs
 */
int alea_copy_surfaces_with_remap(alea_system_t* dst,
                                 const alea_system_t* src,
                                 primitive_remap_t* remap);

/**
 * @brief Copy only materials referenced by cells in dst
 */
int alea_copy_referenced_materials(alea_system_t* dst, const alea_system_t* src);

/**
 * @brief Copy only mixtures whose mc_material_id matches a material in dst cells
 */
int alea_copy_referenced_mixtures(alea_system_t* dst, const alea_system_t* src);

/**
 * @brief Copy only transforms referenced by cells and surfaces in dst
 */
int alea_copy_referenced_transforms(alea_system_t* dst, const alea_system_t* src);

/**
 * @brief Copy cell_refs where both referencing and referenced cells exist in dst
 */
int alea_copy_referenced_cell_refs(alea_system_t* dst, const alea_system_t* src);

/**
 * @brief Clone a CSG tree from one system to another
 */
alea_node_id_t alea_clone_tree_to_system(alea_system_t* dst,
                                        const alea_system_t* src,
                                        alea_node_id_t src_root,
                                        primitive_remap_t* remap);


/* Lattice index lookup (used by raycast DDA) */
struct alea_cell_entry;
int lattice_hex_lookup(const struct alea_cell_entry* cell,
                       double px, double py, double pz,
                       double* ox, double* oy, double* oz);

#endif /* ALEA_UNIVERSE_H */
