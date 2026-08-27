// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/**
 * @file alea.h
 * @brief Alea Public API
 *
 * Consolidated API for Constructive Solid Geometry operations.
 * Prefix: alea_ for all public functions.
 *
 * Basic usage:
 *   // Programmatic geometry
 *   alea_system_t* sys = alea_create();
 *   alea_build_universe_index(sys);
 *   int cell_id, material_id;
 *   alea_find_cell_at(sys, 0.0, 0.0, 0.0, &cell_id, &material_id);
 *   alea_destroy(sys);
 *
 * For MCNP I/O, include "alea_mcnp.h" and link libalea_mcnp.a:
 *   mcnp_model_t* model = mcnp_load("input.inp");
 *   mcnp_export(model, "output.inp");
 *   mcnp_model_destroy(model);
 */

#ifndef ALEA_H
#define ALEA_H

#include <stddef.h>
#include <stdint.h>
#include <stdbool.h>
#include <stdio.h>
#include "alea_types.h"
#include "alea_log.h"

/* alea_cell_hit_t is defined in core/alea_universe.h (included via alea_system.h) */
/* void_result_t is forward-declared in alea_types.h */


#define ALEA_MATERIAL_VOID (-1)

/* Mark a public API as deprecated; callers get a compile-time warning with the
 * given message (GCC/Clang). No effect on other compilers. */
#if defined(__GNUC__) || defined(__clang__)
#  define ALEA_DEPRECATED(msg) __attribute__((deprecated(msg)))
#else
#  define ALEA_DEPRECATED(msg)
#endif

#ifdef __cplusplus
extern "C" {
#endif

/* ============================================================================
 * VERSION
 * ============================================================================ */

#define ALEA_VERSION_MAJOR 0
#define ALEA_VERSION_MINOR 1
#define ALEA_VERSION_PATCH 0

const char* alea_version(void);
int alea_openmp_enabled(void);
int alea_openmp_max_threads(void);

/* ============================================================================
 * ERROR HANDLING
 * ============================================================================ */

const char* alea_error(void);
int alea_error_code(void);
void alea_error_clear(void);

/* ============================================================================
 * INTERRUPT SUPPORT
 *
 * Cooperative interrupt mechanism for cancelling long-running operations.
 * The Python binding's SIGINT handler calls alea_interrupt(), then the
 * next iteration of any long-running loop checks the flag and returns
 * ALEA_ERR_INTERRUPTED.
 *
 * Usage from Python C extension:
 *   static void sigint_handler(int sig) { alea_interrupt(); }
 *   // Before long operation:
 *   old_handler = signal(SIGINT, sigint_handler);
 *   Py_BEGIN_ALLOW_THREADS
 *   result = mcnp_load(filename);  // from alea_mcnp.h
 *   Py_END_ALLOW_THREADS
 *   signal(SIGINT, old_handler);
 *   alea_clear_interrupt();
 * ============================================================================ */

/**
 * @brief Request interruption of the current operation
 *
 * Signal-safe: can be called from a signal handler.
 */
void alea_interrupt(void);

/**
 * @brief Clear the interrupt flag
 *
 * Call after handling an interrupted operation, before starting a new one.
 */
void alea_clear_interrupt(void);

/**
 * @brief Check if an interrupt has been requested
 * @return true if interrupted
 */
bool alea_interrupted(void);

/* ============================================================================
 * SYSTEM LIFECYCLE
 * ============================================================================ */

alea_system_t* alea_create(void);
void alea_destroy(alea_system_t* sys);
alea_system_t* alea_clone(const alea_system_t* sys);
void alea_reset(alea_system_t* sys);

/* ============================================================================
 * CONFIGURATION
 * ============================================================================ */

/**
 * Get the current configuration of a system.
 *
 * Usage: start from the returned config or ALEA_CONFIG_DEFAULT, tweak fields,
 * then call alea_set_config().
 */
alea_config_t alea_get_config(const alea_system_t* sys);

/**
 * Set the configuration of a system.
 *
 * Example:
 *   alea_config_t cfg = alea_get_config(sys);
 *   cfg.dedup = false;
 *   alea_set_config(sys, &cfg);
 */
void alea_set_config(alea_system_t* sys, const alea_config_t* config);

/* ============================================================================
 * LOADING
 *
 * For MCNP: use mcnp_load() from alea_mcnp.h (link libalea_mcnp.a)
 * For OpenMC: use openmc_load() from alea_openmc.h (link libalea_openmc.a)
 * ============================================================================ */

/* ============================================================================
 * INDEXING
 * ============================================================================ */

int alea_build_universe_index(alea_system_t* sys);

/**
 * @brief Prepare query caches (hierarchical spatial index, surface BVH, etc.).
 *
 * This is the query setup API for raycast, slice/grid, mesh, render, and
 * geometry queries.
 */
int alea_prepare_query_acceleration(alea_system_t* sys);

typedef struct {
    bool built;

    /* Hierarchical spatial index shape. */
    size_t hier_universe_count;
    size_t hier_blas_count;
    size_t hier_linear_universe_count;
    size_t hier_blas_cell_count;
    size_t hier_blas_node_count;
    size_t hier_fill_cell_count;
    size_t hier_lattice_cell_count;
    size_t hier_transform_count;
    size_t hier_placement_count;
    size_t hier_root_placement_count;
    size_t hier_fill_placement_count;
    size_t hier_lattice_placement_count;
    int hier_max_placement_depth;
    int hier_max_universe_cells;
    int hier_largest_universe_id;
    size_t memory_bytes;
} alea_query_acceleration_stats_t;

/**
 * @brief Inspect the currently built query-acceleration structure.
 *
 * This does not build caches. Call alea_prepare_query_acceleration() first when
 * a built index is required. Callers can inspect placement, BLAS, transform,
 * and memory counts of the hierarchical spatial index.
 */
int alea_query_acceleration_stats(const alea_system_t* sys,
                                  alea_query_acceleration_stats_t* out_stats);

/** One exact MCNP surface-card reference in the retained reverse index. */
typedef struct {
    uint32_t cell_index;
    int cell_id;
    int universe_id;
    int8_t sense;
} alea_surface_cell_reference_t;

typedef struct {
    bool built;
    size_t surface_count;
    size_t reference_count;
    size_t max_references_per_surface;
    size_t memory_bytes;
} alea_surface_reference_stats_t;

/** Inspect the exact surface-to-cell reverse index without building it. */
int alea_surface_reference_stats(const alea_system_t* sys,
                                 alea_surface_reference_stats_t* out_stats);

/** Copy references to an exact MCNP surface card.
 *
 * ``out_count`` receives the full number of references even when ``capacity``
 * is smaller. ``out_refs`` may be NULL when capacity is zero, allowing a
 * count-only query. The adjacency query cache must already be prepared.
 */
int alea_surface_cell_references(const alea_system_t* sys,
                                 int surface_id,
                                 alea_surface_cell_reference_t* out_refs,
                                 size_t capacity,
                                 size_t* out_count);

/* ============================================================================
 * GEOMETRY QUERIES
 * ============================================================================ */

/* Legacy convenience wrapper.
 * Prefer `alea_find_cell_at()` for single-point queries, or
 * `alea_find_all_cells()` when hierarchy detail is needed. */
int alea_find_cell(alea_system_t* sys, double x, double y, double z);
int alea_find_all_cells(alea_system_t* sys, double x, double y, double z,
                            alea_cell_hit_t* hits, size_t max_hits);
/** Complete diagnostic ownership query with concrete occurrence ancestry.
 * Returns the retained hit count, which equals max_hits when the caller must
 * treat the owner set as potentially truncated. */
int alea_find_all_cells_coverage_chain(alea_system_t* sys,
                                       double x, double y, double z,
                                       alea_cell_hit_t* hits,
                                       uint64_t* occurrence_keys,
                                       uint64_t* parent_occurrence_keys,
                                       size_t max_hits);

/** Complete diagnostic ownership query starting in a universe-local frame. */
int alea_find_all_cells_in_universe_coverage_chain(
    alea_system_t* sys, int universe_id,
    double x, double y, double z,
    alea_cell_hit_t* hits,
    uint64_t* occurrence_keys,
    uint64_t* parent_occurrence_keys,
    size_t max_hits);

typedef enum {
    ALEA_POINT_COVERAGE_UNIQUE = 0,
    ALEA_POINT_COVERAGE_GAP,
    ALEA_POINT_COVERAGE_OVERLAP,
    ALEA_POINT_COVERAGE_UNDEFINED_FILL,
    ALEA_POINT_COVERAGE_UNRESOLVED
} alea_point_coverage_kind_t;

typedef struct {
    alea_point_coverage_kind_t kind;
    int target_depth; /**< -1 for a resolved leaf set spanning depths. */
    size_t owner_count;
} alea_point_coverage_classification_t;

/** Classify an occurrence tree returned by a complete point query.
 *
 * ``universe_depth < 0`` selects every terminal leaf, so competing branches
 * remain overlaps even when their leaves occur at different depths. At an
 * explicit depth, only claimants projected to that depth are selected.
 * ``out_owner_mask`` is parallel to ``hits`` and may be NULL.
 */
int alea_classify_point_coverage_chain(
    const alea_cell_hit_t* hits,
    const uint64_t* occurrence_keys,
    const uint64_t* parent_occurrence_keys,
    size_t hit_count,
    int universe_depth,
    uint8_t* out_owner_mask,
    alea_point_coverage_classification_t* out_classification);

/** Allocation-free form of alea_classify_point_coverage_chain().
 *
 * ``child_count_scratch`` must contain at least ``hit_count`` elements. It is
 * overwritten by the classifier. This form is intended for hot native loops
 * that already own bounded per-worker scratch storage.
 */
int alea_classify_point_coverage_chain_with_scratch(
    const alea_cell_hit_t* hits,
    const uint64_t* occurrence_keys,
    const uint64_t* parent_occurrence_keys,
    size_t hit_count,
    int universe_depth,
    uint8_t* out_owner_mask,
    size_t* child_count_scratch,
    alea_point_coverage_classification_t* out_classification);
bool alea_point_inside(const alea_system_t* sys, alea_node_id_t node,
                           double x, double y, double z);
/* Legacy convenience wrapper.
 * Prefer `alea_find_cell_at()` for single-point queries. */
int alea_material_at(alea_system_t* sys, double x, double y, double z);
int alea_find_overlaps(alea_system_t* sys, int* pairs, size_t max_pairs);

/* ============================================================================
 * CSG CONSTRUCTION - SURFACES
 *
 * These functions create a complete surface entry with both positive and
 * negative sense nodes registered in the system. Returns the index into
 * sys->surfaces array. Use alea_halfspace() to get +/- sense node.
 *
 * Pass surface_id=0 for automatic ID assignment.
 *
 * Example:
 *   int idx = alea_sphere_surface(sys, 0, 0, 0, 0, 5.0);
 *   alea_node_id_t interior = alea_halfspace(sys, idx, -1);
 * ============================================================================ */

/**
 * @brief Create a plane surface with automatic registration
 *
 * The plane is defined by the equation  a·x + b·y + c·z + d = 0.
 * The negative half-space (sense=-1) is where a·x + b·y + c·z + d < 0.
 *
 * To place a plane at x = K with outward normal (+x): use a=1, b=0, c=0, d=-K.
 * Example: plane at x = 4  →  alea_plane_surface(sys, id, 1, 0, 0, -4).
 *
 * @param surface_id MCNP-style surface ID (0 = auto-assign)
 * @param a,b,c     Normal direction (need not be unit length; will be normalised)
 * @param d         Constant term in the plane equation (d = −K for plane at position K)
 * @return Index into surfaces array, or -1 on error
 */
int alea_plane_surface(alea_system_t* sys, int surface_id,
                           double a, double b, double c, double d);

/**
 * @brief Create a sphere surface with automatic registration
 */
int alea_sphere_surface(alea_system_t* sys, int surface_id,
                            double cx, double cy, double cz, double r);

/**
 * @brief Create a Z-axis cylinder surface with automatic registration
 */
int alea_cylinder_z_surface(alea_system_t* sys, int surface_id,
                                double cx, double cy, double r);

/**
 * @brief Create an X-axis cylinder surface with automatic registration
 */
int alea_cylinder_x_surface(alea_system_t* sys, int surface_id,
                                double cy, double cz, double r);

/**
 * @brief Create a Y-axis cylinder surface with automatic registration
 */
int alea_cylinder_y_surface(alea_system_t* sys, int surface_id,
                                double cx, double cz, double r);

/**
 * @brief Create a box surface with automatic registration
 */
int alea_box_surface(alea_system_t* sys, int surface_id,
                         double xmin, double xmax,
                         double ymin, double ymax,
                         double zmin, double zmax);

/**
 * @brief Create a Z-axis cone surface with automatic registration
 */
int alea_cone_z_surface(alea_system_t* sys, int surface_id,
                            double cx, double cy, double cz,
                            double t_squared);

/**
 * @brief Create an X-axis cone surface with automatic registration
 */
int alea_cone_x_surface(alea_system_t* sys, int surface_id,
                            double cx, double cy, double cz,
                            double t_squared);

/**
 * @brief Create a Y-axis cone surface with automatic registration
 */
int alea_cone_y_surface(alea_system_t* sys, int surface_id,
                            double cx, double cy, double cz,
                            double t_squared);

/**
 * @brief Create a Z-axis torus surface with automatic registration
 * @param major_radius Distance from axis to tube center
 * @param minor_radius Tube radius
 */
int alea_torus_z_surface(alea_system_t* sys, int surface_id,
                             double cx, double cy, double cz,
                             double major_radius, double minor_radius);

/**
 * @brief Create an X-axis torus surface with automatic registration
 */
int alea_torus_x_surface(alea_system_t* sys, int surface_id,
                             double cx, double cy, double cz,
                             double major_radius, double minor_radius);

/**
 * @brief Create a Y-axis torus surface with automatic registration
 */
int alea_torus_y_surface(alea_system_t* sys, int surface_id,
                             double cx, double cy, double cz,
                             double major_radius, double minor_radius);

/**
 * @brief Create a general quadric surface with automatic registration
 * Ax² + By² + Cz² + Dxy + Eyz + Fxz + Gx + Hy + Iz + J = 0
 */
int alea_quadric_surface(alea_system_t* sys, int surface_id,
                             double A, double B, double C,
                             double D, double E, double F,
                             double G, double H, double I, double J);

/* --- Macrobody surfaces --- */

/**
 * @brief Create RCC (Right Circular Cylinder) macrobody surface
 * @param base_x, base_y, base_z Base center point
 * @param height_x, height_y, height_z Height vector (axis direction * length)
 * @param radius Cylinder radius
 */
int alea_rcc_surface(alea_system_t* sys, int surface_id,
                         double base_x, double base_y, double base_z,
                         double height_x, double height_y, double height_z,
                         double radius);

/**
 * @brief Create BOX (general oriented box) macrobody surface
 * @param corner_x, corner_y, corner_z Corner vertex
 * @param v1_x, v1_y, v1_z First edge vector
 * @param v2_x, v2_y, v2_z Second edge vector
 * @param v3_x, v3_y, v3_z Third edge vector
 */
int alea_box_general_surface(alea_system_t* sys, int surface_id,
                                 double corner_x, double corner_y, double corner_z,
                                 double v1_x, double v1_y, double v1_z,
                                 double v2_x, double v2_y, double v2_z,
                                 double v3_x, double v3_y, double v3_z);

/**
 * @brief Create SPH (sphere) macrobody surface
 */
int alea_sph_surface(alea_system_t* sys, int surface_id,
                         double cx, double cy, double cz, double r);

/**
 * @brief Create TRC (Truncated Right Cone) macrobody surface
 * @param base_x, base_y, base_z Base center point
 * @param height_x, height_y, height_z Height vector
 * @param base_radius Radius at base
 * @param top_radius Radius at top
 */
int alea_trc_surface(alea_system_t* sys, int surface_id,
                         double base_x, double base_y, double base_z,
                         double height_x, double height_y, double height_z,
                         double base_radius, double top_radius);

/**
 * @brief Create ELL (Ellipsoid) macrobody surface
 * @param v1_x, v1_y, v1_z First focus
 * @param v2_x, v2_y, v2_z Second focus
 * @param major_axis_len Length of major axis (2a)
 */
int alea_ell_surface(alea_system_t* sys, int surface_id,
                         double v1_x, double v1_y, double v1_z,
                         double v2_x, double v2_y, double v2_z,
                         double major_axis_len);

/**
 * @brief Create REC (Right Elliptical Cylinder) macrobody surface
 * @param base_x, base_y, base_z Base center point
 * @param height_x, height_y, height_z Height vector
 * @param axis1_x, axis1_y, axis1_z First semi-axis vector
 * @param axis2_x, axis2_y, axis2_z Second semi-axis vector (perpendicular)
 */
int alea_rec_surface(alea_system_t* sys, int surface_id,
                         double base_x, double base_y, double base_z,
                         double height_x, double height_y, double height_z,
                         double axis1_x, double axis1_y, double axis1_z,
                         double axis2_x, double axis2_y, double axis2_z);

/**
 * @brief Create WED (Wedge) macrobody surface
 * @param vertex_x, vertex_y, vertex_z Vertex of right angle
 * @param v1_x, v1_y, v1_z First leg of right triangle
 * @param v2_x, v2_y, v2_z Second leg (perpendicular to v1)
 * @param v3_x, v3_y, v3_z Height/extrusion vector
 */
int alea_wed_surface(alea_system_t* sys, int surface_id,
                         double vertex_x, double vertex_y, double vertex_z,
                         double v1_x, double v1_y, double v1_z,
                         double v2_x, double v2_y, double v2_z,
                         double v3_x, double v3_y, double v3_z);

/**
 * @brief Create RHP (Right Hexagonal Prism) macrobody surface
 * @param base_x, base_y, base_z Base center
 * @param height_x, height_y, height_z Height vector
 * @param r1_x, r1_y, r1_z Vector to first vertex pair
 * @param r2_x, r2_y, r2_z Vector to second vertex pair
 * @param r3_x, r3_y, r3_z Vector to third vertex pair
 */
int alea_rhp_surface(alea_system_t* sys, int surface_id,
                         double base_x, double base_y, double base_z,
                         double height_x, double height_y, double height_z,
                                           double r1_x, double r1_y, double r1_z,
                                           double r2_x, double r2_y, double r2_z,
                                           double r3_x, double r3_y, double r3_z);

/**
 * @brief Get sense node from a surface index
 * @param surface_index Index returned by alea_*_surface()
 * @param sense Requested sense: -1 for negative (inside), +1 for positive (outside)
 * @return Node ID or ALEA_NODE_ID_INVALID if out of range
 */
alea_node_id_t alea_halfspace(const alea_system_t* sys, int surface_index, int sense);

/* ============================================================================
 * CSG CONSTRUCTION - OPERATIONS
 * ============================================================================ */

alea_node_id_t alea_union(alea_system_t* sys, alea_node_id_t a, alea_node_id_t b);
alea_node_id_t alea_intersection(alea_system_t* sys, alea_node_id_t a, alea_node_id_t b);
alea_node_id_t alea_difference(alea_system_t* sys, alea_node_id_t a, alea_node_id_t b);
alea_node_id_t alea_complement(alea_system_t* sys, alea_node_id_t a);

alea_node_id_t alea_union_n(alea_system_t* sys, const alea_node_id_t* nodes, size_t count);
alea_node_id_t alea_intersection_n(alea_system_t* sys, const alea_node_id_t* nodes, size_t count);

/* ============================================================================
 * CSG CONSTRUCTION - CELLS
 * ============================================================================ */

int alea_add_cell(alea_system_t* sys, int cell_id, alea_node_id_t root,
                      int material_index, double density, int universe);
int alea_set_fill(alea_system_t* sys, int cell_index, int fill_universe, int transform);
int alea_cell_set_comment(alea_system_t* sys, int cell_index, const char* comment);
int alea_cell_set_inline_comment(alea_system_t* sys, int cell_index, const char* comment);

/**
 * @brief Set cell material
 * @param cell_index Cell index (0 to alea_cell_count-1)
 * @param material_index Material index from alea_add_material(), or ALEA_MATERIAL_VOID (-1) for void
 * @return 0 on success, -1 on error
 */
int alea_cell_set_material(alea_system_t* sys, int cell_index, int material_index);

/**
 * @brief Set cell material to a mixture
 * @param cell_index Cell index (0 to alea_cell_count-1)
 * @param mixture_id Mixture ID from alea_create_mixture()
 * @return 0 on success, -1 on error
 */
int alea_cell_set_mixture(alea_system_t* sys, int cell_index, int mixture_id);

/**
 * @brief Set cell density
 * @param cell_index Cell index (0 to alea_cell_count-1)
 * @param density Signed density: negative = g/cm3, positive = atoms/b-cm, 0 = void
 * @return 0 on success, -1 on error
 */
int alea_cell_set_density(alea_system_t* sys, int cell_index, double density);

/**
 * @brief Set cell universe membership
 * @param cell_index Cell index (0 to alea_cell_count-1)
 * @param universe_id Universe ID (0 = base universe)
 * @return 0 on success, -1 on error
 */
int alea_cell_set_universe(alea_system_t* sys, int cell_index, int universe_id);

/**
 * @brief Remove a cell by index
 *
 * Frees per-cell allocations, compacts the array, and rebuilds the
 * cell ID hashmap. Triggers on_cell_removed hook before removal.
 * Invalidates the universe index.
 *
 * @param cell_index Cell index (0 to alea_cell_count-1)
 * @return 0 on success, -1 on error
 */
int alea_cell_remove(alea_system_t* sys, int cell_index);

/* ============================================================================
 * UNIVERSE OPERATIONS
 * ============================================================================ */

int alea_flatten(alea_system_t* sys, int universe_id);

/**
 * @brief Statistics from simplification pass
 */
#ifndef ALEA_SIMPLIFY_STATS_DEFINED
#define ALEA_SIMPLIFY_STATS_DEFINED
typedef struct {
    size_t nodes_before;
    size_t nodes_after;
    size_t complements_eliminated;
    size_t double_negations;
    size_t idempotent_reductions;
    size_t absorption_reductions;
    size_t subtrees_deduplicated;
    size_t cell_complements_expanded;
    size_t contradictions_found;
    size_t tautologies_found;
    size_t empty_cells_removed;
    size_t union_branches_absorbed;
    size_t union_common_factors;
    size_t union_branches_subsumed;
} alea_simplify_stats_t;
#endif

/**
 * @brief Simplify all cells and remove cells proven empty
 *
 * Applies NNF conversion, associative Boolean normalization, balancing,
 * contradiction detection, and empty-cell removal.
 *
 * @param sys CSG system
 * @param stats Optional accumulated stats output (can be NULL)
 */
void alea_simplify_and_prune_cells(alea_system_t* sys,
                                   alea_simplify_stats_t* stats);

/**
 * @brief Carve a region out of every ordinary cell in one universe.
 *
 * Each affected cell is replaced by ``cell_region - carve_root``.  The
 * carve tree must belong to @p sys.  When @p simplify is non-zero, each new
 * tree is simplified before it is installed and cells proved empty are
 * removed.  Lattice cells are rejected: their containment is defined by
 * lattice-specific lookup rules, so attaching a generic CSG mask to their
 * root would not be reliable.
 *
 * This is intended for same-universe component insertion: construct the
 * union of the incoming component regions and carve the pre-existing cells.
 * Callers that register incoming cells first should pass their original cell
 * count as @p cell_limit so those cells are excluded from the rewrite.
 *
 * @param sys System to modify
 * @param universe_id Target universe
 * @param carve_root Root of the region to remove
 * @param simplify Whether to simplify each rewritten cell
 * @param cell_limit Only cells with an index below this value are rewritten;
 *        pass -1 to rewrite every cell in the universe.
 * @param out_modified Optional number of rewritten cells
 * @param out_removed Optional number of cells removed as empty
 * @return 0 on success, -1 on error (the target cells are unchanged on
 *         validation or construction failure)
 */
int alea_carve_universe(alea_system_t* sys, int universe_id,
                        alea_node_id_t carve_root, int simplify, int cell_limit,
                        int* out_modified, int* out_removed);

/**
 * @brief Split cells with top-level unions into multiple simpler cells.
 *
 * For each cell whose root is a union (T1 ∪ T2 ∪ ... ∪ Tk),
 * creates k new cells (one per branch), inheriting material, density,
 * and universe from the original. The original cell is removed.
 *
 * New cells receive auto-assigned cell IDs.
 * Call alea_build_universe_index() after splitting.
 *
 * @return Number of new cells created, or -1 on error.
 */
int alea_split_union_cells(alea_system_t* sys);

alea_system_t* alea_extract_universe(const alea_system_t* sys, int universe_id);
int alea_merge(alea_system_t* target, const alea_system_t* source, int id_offset);

/* ============================================================================
 * INFORMATION
 * ============================================================================ */

size_t alea_cell_count(const alea_system_t* sys);
size_t alea_surface_count(const alea_system_t* sys);
size_t alea_universe_count(const alea_system_t* sys);

void alea_stats(const alea_system_t* sys, alea_stats_t* stats);
void alea_print_summary(const alea_system_t* sys);
void alea_tree_print(const alea_system_t* sys, alea_node_id_t node_id);

/* ============================================================================
 * VOLUME ESTIMATION
 * ============================================================================ */

/**
 * @brief Compute a tight bounding sphere for the entire model
 *
 * For each cell, binary-searches each axis to tighten the bounding box
 * using interval arithmetic.  The union of all tight bboxes gives the
 * model extent.  The bounding sphere is centered on this extent.
 *
 * @param sys      CSG system (read-only)
 * @param tol      Spatial tolerance for bbox tightening (e.g. 1.0)
 * @param cx,cy,cz Output: center of bounding sphere
 * @param radius   Output: radius of bounding sphere
 * @return 0 on success, -1 on error (no bounded cells)
 */
int alea_compute_bounding_sphere(alea_system_t* sys,
                                      double tol,
                                      double* cx, double* cy, double* cz,
                                      double* radius);

/**
 * @brief Estimate cell volumes using random ray tracing (Cauchy-Crofton method)
 *
 * Generates @p n_rays uniform random lines through a sphere of the given
 * @p radius centered at (@p ox, @p oy, @p oz). Accumulates track lengths
 * per cell and converts to volume: V_cell = pi * R^2 * (sum_L / N).
 *
 * If @p rel_errors is non-NULL, it receives the 1-sigma relative statistical
 * error for each cell: sigma_V / V.  Set to -1 for zero-volume cells.
 *
 * Requires the raycast module to be linked; returns -1 otherwise.
 *
 * @param sys       CSG system (read-only)
 * @param ox,oy,oz  Center of the bounding sphere
 * @param radius    Bounding sphere radius
 * @param n_rays    Number of random rays (e.g. 10000)
 * @param volumes    Output array of size alea_cell_count(sys)
 * @param rel_errors Output array of size alea_cell_count(sys), or NULL
 * @return 0 on success, -1 on error
 */
int alea_estimate_cell_volumes(alea_system_t* sys,
                                   double ox, double oy, double oz,
                                   double radius, int n_rays,
                                   double* volumes,
                                   double* rel_errors);

#define ALEA_VOLUME_PATH_MAX_DEPTH 16

typedef struct {
    int lattice_cell_index;
    int fill_universe;
    int i, j, k;
    int linear_index;
} alea_volume_lattice_step_t;

typedef struct {
    uint64_t path_id;

    int terminal_cell_index;
    int terminal_cell_id;
    int material_id;
    int universe_id;
    int depth;

    uint8_t ancestor_count;
    uint8_t lattice_step_count;
    uint8_t truncated;
    uint8_t reserved;

    int ancestor_cell_indices[ALEA_VOLUME_PATH_MAX_DEPTH];
    int ancestor_universe_ids[ALEA_VOLUME_PATH_MAX_DEPTH];
    alea_volume_lattice_step_t lattice_steps[ALEA_VOLUME_PATH_MAX_DEPTH];
    double world_to_local[12];
} alea_volume_path_t;

typedef enum {
    ALEA_VOLUME_PATH_NONE = 0,
    ALEA_VOLUME_PATH_UNIQUE = 1,
    ALEA_VOLUME_PATH_MULTIPLE = 2
} alea_volume_path_status_t;

/** A contiguous point range sharing one target cell occurrence. */
typedef struct {
    size_t point_offset;
    size_t point_count;
    int target_cell_id;
    int target_universe_id;
} alea_volume_path_point_set_t;

/** Aggregate occurrence evidence across one point set. */
typedef struct {
    alea_volume_path_status_t status;
    size_t tested_point_count;
    size_t matching_occurrence_count;
    alea_volume_path_t unique_path;
} alea_volume_path_point_set_result_t;

typedef struct {
    size_t point_set_count;
    size_t completed_point_set_count;
    size_t requested_workers;
    size_t actual_workers;
    uint64_t reserved_scratch_bytes_per_worker;
    uint64_t reserved_parallel_scratch_bytes;
} alea_volume_path_point_set_batch_stats_t;

/**
 * @brief Return the number of concrete hierarchical volume paths.
 *
 * Hierarchical mode only in this implementation. The count is stable until
 * geometry/config mutation invalidates query caches. Use this to size arrays
 * for alea_volume_paths_get() and alea_estimate_path_volumes().
 */
size_t alea_volume_path_count(alea_system_t* sys);

/**
 * @brief Enumerate concrete hierarchical volume paths.
 *
 * Writes up to @p max_paths entries into @p out_paths and returns the total
 * number of available paths. If the return value is greater than @p max_paths,
 * the output was truncated.
 */
size_t alea_volume_paths_get(alea_system_t* sys,
                             alea_volume_path_t* out_paths,
                             size_t max_paths);

/**
 * @brief Resolve a world-space point to a concrete volume path.
 *
 * Hierarchical mode only in this implementation. The returned path has the same
 * structural identity as paths from alea_volume_paths_get(), including the
 * dense path_id when the path exists in the current enumeration.
 */
int alea_volume_path_at_point(alea_system_t* sys,
                              double x, double y, double z,
                              alea_volume_path_t* out_path);

/** Resolve one world-space point without enumerating all concrete paths.
 *
 * This walks only the hierarchy branch containing the point and therefore
 * does not build the potentially large global volume-path index. The
 * structural path, ancestors, lattice steps, and world-to-local transform are
 * populated normally; ``path_id`` is set to ``UINT64_MAX`` because no dense
 * global identifier is assigned.
 */
int alea_volume_path_resolve_at_point(alea_system_t* sys,
                                      double x, double y, double z,
                                      alea_volume_path_t* out_path);

/** Resolve a specified cell occurrence at one world-space point.
 *
 * Unlike alea_volume_path_resolve_at_point(), this explores every hierarchy
 * branch that contains the point and can therefore distinguish a unique
 * occurrence from an overlap.  The target may be a terminal cell or a
 * fill/lattice container at any hierarchy depth.  @p target_universe_id is
 * the universe in which the MCNP cell card is defined.
 *
 * The traversal stops after finding a second distinct occurrence: callers
 * only need the zero/unique/multiple classification, so this is a semantic
 * early exit rather than a model-size cap.  On a unique result, @p out_path
 * describes the target occurrence itself and its world-to-target-universe
 * transform; ``path_id`` is ``UINT64_MAX`` because no global path index is
 * built.
 *
 * @return 0 when absent, 1 when unique, 2 when multiple, or -1 on error.
 */
int alea_volume_path_resolve_cell_at_point(alea_system_t* sys,
                                           double x, double y, double z,
                                           int target_cell_id,
                                           int target_universe_id,
                                           alea_volume_path_t* out_path);

/** Resolve target-cell occurrences across independent point sets.
 *
 * ``points_xyz`` contains ``3 * point_count`` doubles. Each point-set
 * descriptor selects a contiguous point range and one target cell/universe.
 * A result is unique only when every matching point identifies the same
 * structural occurrence and no individual point is ambiguous. The hierarchy
 * query cache is prepared once before any parallel work. A zero worker request
 * uses the OpenMP runtime maximum; a zero parallel scratch budget selects the
 * serial path. Output order is always point-set input order.
 */
int alea_volume_path_resolve_cell_point_sets(
    alea_system_t* sys,
    const double* points_xyz,
    size_t point_count,
    const alea_volume_path_point_set_t* point_sets,
    size_t point_set_count,
    size_t requested_workers,
    uint64_t max_parallel_scratch_bytes,
    alea_volume_path_point_set_result_t* results,
    alea_volume_path_point_set_batch_stats_t* out_stats);

/** Resolve a cell occurrence from corresponding world/local coordinates.
 *
 * This identifies the placement of @p target_universe_id by comparing the
 * supplied world-space point/direction with the same event printed in that
 * universe's local frame.  It deliberately does not require the target cell
 * to contain the point, which makes it suitable for diagnosing a reported
 * boundary loss.  Only ordinary FILL placements and the root universe are
 * supported; concrete lattice elements require containment-based resolution.
 *
 * Each precision is one last-printed-place for the corresponding coordinate.
 * The comparison allows half a printed unit from each frame plus a small
 * floating-point allowance.  No global volume-path index is built.
 *
 * @return 0 when absent, 1 when unique, 2 when multiple, or -1 on error.
 */
int alea_volume_path_resolve_cell_from_transform_evidence(
    alea_system_t* sys,
    int target_cell_id,
    int target_universe_id,
    const double world_point[3],
    const double world_direction[3],
    const double local_point[3],
    const double local_direction[3],
    const double world_point_precision[3],
    const double world_direction_precision[3],
    const double local_point_precision[3],
    const double local_direction_precision[3],
    alea_volume_path_t* out_path);

/**
 * @brief Estimate volumes for concrete hierarchy paths.
 *
 * Output arrays must be sized to alea_volume_path_count(sys).
 */
int alea_estimate_path_volumes(alea_system_t* sys,
                               int n_rays,
                               double* volumes,
                               double* rel_errors);

/**
 * @brief Remove cells whose estimated volume is below a threshold
 *
 * Requires the raycast module to be linked; returns -1 otherwise.
 *
 * @param sys        CSG system (modified in place)
 * @param volumes    Volume array from alea_estimate_cell_volumes
 * @param threshold  Cells with volume <= threshold are removed
 * @return Number of cells removed, or -1 on error
 */
int alea_remove_cells_by_volume(alea_system_t* sys,
                                    const double* volumes,
                                    double threshold);

/* ============================================================================
 * BOUNDING BOX TIGHTENING
 * ============================================================================ */

/**
 * @brief Tighten a single cell's bounding box via interval arithmetic
 *
 * @param sys   CSG system
 * @param cell_index  Cell index (0 to alea_cell_count-1)
 * @param tol   Spatial tolerance
 * @param out   Output tightened bbox
 * @return 0 on success, -1 on error
 */
int alea_tighten_cell_bbox(const alea_system_t* sys,
                                size_t cell_index,
                                double tol,
                                alea_bbox_t* out);

/**
 * @brief Tighten all cell bboxes in the system
 *
 * @param sys  CSG system (modified: node bboxes updated)
 * @param tol  Spatial tolerance
 * @return Number of cells tightened, or -1 on error
 */
int alea_tighten_all_bboxes(alea_system_t* sys, double tol);

/**
 * @brief Numerically tighten a cell's bounding box
 *
 * Uses LP vertex enumeration (for plane-only cells) and/or octree
 * recursive bisection (for general cells) to find a finite bounding
 * box for cells that have infinite analytical bboxes.
 *
 * @param sys        CSG system (modified: node bbox updated on success)
 * @param cell_index Cell index (0 to alea_cell_count-1)
 * @return 0 on success, -1 on error (no finite bbox found)
 */
int alea_tighten_cell_bbox_numerical(alea_system_t* sys, int cell_index);

/* ============================================================================
 * EXPORT
 *
 * For MCNP: use mcnp_export() from alea_mcnp.h (link libalea_mcnp.a)
 * For OpenMC: use openmc_export() from alea_openmc.h (link libalea_openmc.a)
 * For Serpent: use serpent_export_system() from alea_serpent.h
 * (link libalea_serpent.a)
 *
 * Generic export options are read from sys->config (dedup, universe_depth,
 * fill_depth, export_materials, export_transforms).
 * ============================================================================ */

/* ============================================================================
 * VOID GENERATION
 *
 * Void options are read from sys->config (void_max_depth, void_min_size,
 * void_probes_per_axis, merge_cell_weight, merge_surface_weight, etc.).
 * ============================================================================ */

/* Generate voids inside an existing finite CSG region. The supplied bounds
 * region is preserved and reused during consolidation when possible. */
void_result_t* alea_void_generate_in_region(alea_system_t* sys, alea_node_id_t bounds_region);

/* Generate voids inside a bbox. If bounds is NULL, universe 0 bbox plus a
 * margin is used. The bbox path creates an internal RPP bounds region. */
void_result_t* alea_void_generate_in_bbox(alea_system_t* sys, const alea_bbox_t* bounds);
int alea_void_add_cells(alea_system_t* sys, void_result_t* result);
int alea_void_add_graveyard(alea_system_t* sys, void_result_t* result);
size_t alea_void_count(const void_result_t* result);
int alea_void_get(const void_result_t* result, size_t index, alea_bbox_t* box);
void alea_void_free(void_result_t* result);

/**
 * @brief Convert void result to CSG node (union of boxes)
 */
alea_node_id_t alea_void_to_node(alea_system_t* sys, const void_result_t* result);

/**
 * @brief Merge void cells to reduce count while balancing complexity
 */
int alea_void_merge(alea_system_t* sys, void_result_t* result);

/* ============================================================================
 * LOGGING
 * ============================================================================ */

/* ============================================================================
 * MATERIAL OPERATIONS
 * ============================================================================ */

/**
 * @brief Register a material in the system
 *
 * Creates an empty material entry and returns its index. Use
 * alea_material_add_nuclide() etc. to populate it.
 *
 * @param sys System
 * @param material_id MCNP-style material ID (0 or negative for auto-assign)
 * @return Material index (>= 0), or -1 on error
 */
int alea_add_material(alea_system_t* sys, int material_id);

/**
 * @brief Find a material by its MCNP ID
 * @return Material index (>= 0), or -1 if not found
 */
int alea_find_material_by_id(const alea_system_t* sys, int material_id);

/**
 * @brief Get number of materials in the system
 */
size_t alea_material_count(const alea_system_t* sys);

/**
 * @brief Get material MCNP ID by index
 * @return Material ID, or -1 on error
 */
int alea_material_get_id(const alea_system_t* sys, int mat_index);

/**
 * @brief Add a nuclide to a material
 * @param mat_index Material index from alea_add_material()
 * @param zaid ZAID (ZZAAA or ZZAAAI format, e.g. 92235 for U-235)
 * @param library Cross-section library suffix (e.g. ".80c"), or NULL
 * @param fraction Atom or weight fraction (always positive)
 * @return 0 on success, -1 on error
 */
int alea_material_add_nuclide(alea_system_t* sys, int mat_index,
                              int zaid, const char* library, double fraction);

/**
 * @brief Add an element to a material (natural isotopic composition)
 * @param mat_index Material index from alea_add_material()
 * @param Z Atomic number (1-118)
 * @param library Library suffix for generated nuclides, or NULL
 * @param fraction Atom or weight fraction
 * @return 0 on success, -1 on error
 */
int alea_material_add_element(alea_system_t* sys, int mat_index,
                              int Z, const char* library, double fraction);

/**
 * @brief Set material standard density
 * @param mat_index Material index
 * @param density Density (g/cm3 if positive, atoms/b-cm if negative)
 * @return 0 on success, -1 on error
 */
int alea_material_set_density(alea_system_t* sys, int mat_index, double density);

/**
 * @brief Set whether fractions are weight or atom fractions
 * @param mat_index Material index
 * @param is_weight true = weight fractions, false = atom fractions
 * @return 0 on success, -1 on error
 */
int alea_material_set_weight_fraction(alea_system_t* sys, int mat_index, bool is_weight);

/**
 * @brief Expand element components to explicit nuclides
 *
 * Replaces element entries with individual nuclide entries based on
 * natural isotopic abundances.
 *
 * @return 0 on success, -1 on error
 */
int alea_material_expand_elements(alea_system_t* sys, int mat_index);

/**
 * @brief Get nuclide count for a material
 * @return Number of nuclides, or 0 on error
 */
size_t alea_material_nuclide_count(const alea_system_t* sys, int mat_index);

/**
 * @brief Get nuclide data from a material
 * @param nuc_index Nuclide index (0 to nuclide_count-1)
 * @param zaid Output: ZAID (can be NULL)
 * @param library Output: library suffix pointer (can be NULL, do not free)
 * @param fraction Output: fraction (can be NULL)
 * @return 0 on success, -1 on error
 */
int alea_material_nuclide_get(const alea_system_t* sys, int mat_index,
                              size_t nuc_index, int* zaid,
                              const char** library, double* fraction);

/**
 * @brief Get element count for a material
 * @return Number of elements, or 0 on error
 */
size_t alea_material_element_count(const alea_system_t* sys, int mat_index);

/**
 * @brief Get element data from a material
 * @param elem_index Element index (0 to element_count-1)
 * @param Z Output: atomic number (can be NULL)
 * @param library Output: library suffix pointer (can be NULL, do not free)
 * @param fraction Output: fraction (can be NULL)
 * @return 0 on success, -1 on error
 */
int alea_material_element_get(const alea_system_t* sys, int mat_index,
                              size_t elem_index, int* Z,
                              const char** library, double* fraction);

/**
 * @brief Get material density
 * @param density Output: density value
 * @param has_density Output: whether density was set (can be NULL)
 * @return 0 on success, -1 on error
 */
int alea_material_get_density(const alea_system_t* sys, int mat_index,
                              double* density, bool* has_density);

/**
 * @brief Check if material uses weight fractions
 * @return true if weight fractions, false if atom fractions or on error
 */
bool alea_material_is_weight_fraction(const alea_system_t* sys, int mat_index);

/* ============================================================================
 * MIXTURE OPERATIONS
 * ============================================================================ */

/**
 * @brief Create a mixture of materials
 *
 * Creates a new mixture as a weighted combination of existing materials.
 *
 * @param mat_ids Array of material IDs to mix
 * @param fractions Weight fractions
 * @param count Number of materials to mix
 * @param new_mat_id ID for the new mixture (0 for auto-assign)
 * @return Assigned mixture ID on success, -1 on error
 */
int alea_create_mixture(alea_system_t* sys, const int* mat_ids,
                            const double* fractions, size_t count, int new_mat_id);

/**
 * @brief Find mixture index by ID
 * @return Mixture index, or -1 if not found
 */
int alea_find_mixture_by_id(const alea_system_t* sys, int mixture_id);

/**
 * @brief Get number of mixtures in the system
 */
size_t alea_mixture_count(const alea_system_t* sys);

/**
 * @brief Get mixture ID by index
 * @return Mixture ID, or -1 on error
 */
int alea_mixture_get_id(const alea_system_t* sys, int mix_index);

/**
 * @brief Get number of components in a mixture
 * @return Number of components, or 0 on error
 */
size_t alea_mixture_component_count(const alea_system_t* sys, int mix_index);

/**
 * @brief Get mixture component data
 * @param comp_index Component index (0 to component_count-1)
 * @param material_id Output: material ID of this component (can be NULL)
 * @param fraction Output: fraction of this component (can be NULL)
 * @return 0 on success, -1 on error
 */
int alea_mixture_component_get(const alea_system_t* sys, int mix_index,
                               size_t comp_index, int* material_id,
                               double* fraction);

/* ============================================================================
 * EXTRACT / FILTER
 * ============================================================================ */

/**
 * @brief Extract cells in a bounding box region
 *
 * Creates a new system containing cells whose bounding boxes intersect
 * the specified region.
 */
alea_system_t* alea_extract_region(const alea_system_t* sys, const alea_bbox_t* bbox);

/**
 * @brief Get cell indices whose bounding box intersects a region
 */
size_t alea_get_cells_in_bbox(const alea_system_t* sys, const alea_bbox_t* bbox,
                                  int* out_indices, size_t max_count);

/* ============================================================================
 * CELL INFO
 * ============================================================================ */

/**
 * @brief Get cell info by index
 * @param sys System
 * @param index Cell index (0 to alea_cell_count-1)
 * @param cell_id Output: MCNP cell ID (can be NULL)
 * @param material_id Output: Material ID (can be NULL)
 * @param density Output: Density (can be NULL)
 * @param universe_id Output: Universe ID (can be NULL)
 * @param fill_universe Output: Fill universe ID (can be NULL)
 * @param root Output: Root node ID (can be NULL)
 * @return 0 on success, -1 on error
 */
int alea_cell_get(const alea_system_t* sys, size_t index,
                      int* cell_id, int* material_id, double* density,
                      int* universe_id, int* fill_universe, alea_node_id_t* root);

/**
 * @brief Find cell index by MCNP cell ID
 * @return Cell index or -1 if not found
 */
int alea_cell_find(const alea_system_t* sys, int cell_id);

/**
 * @brief Get cell info by index (struct-based)
 * @param info Output cell info structure
 * @return 0 on success, -1 on error
 */
int alea_cell_get_info(const alea_system_t* sys, size_t index, alea_cell_info_t* info);

/**
 * @brief Find cell by ID and get info (struct-based)
 * @param info Output cell info structure
 * @return 0 on success, -1 if not found
 */
int alea_cell_find_info(const alea_system_t* sys, int cell_id, alea_cell_info_t* info);

/**
 * @brief Get CSG expression string for a cell
 *
 * Walks the CSG tree and builds a human-readable expression using the
 * provided operator symbols. The caller must free() the returned string.
 *
 * Example (MCNP-style): alea_cell_expr(sys, i, ":", " ", "#")
 * Example (OpenMC-style): alea_cell_expr(sys, i, " | ", " ", "~")
 *
 * @param sys System
 * @param cell_index Cell index (0 to alea_cell_count-1)
 * @param union_op String to emit between union operands (e.g. ":" or " | ")
 * @param inter_op String to emit between intersection operands (e.g. " ")
 * @param compl_op String to emit before complement operands (e.g. "#" or "~")
 * @return malloc'd string (caller frees), or NULL on error
 */
char* alea_cell_expr(const alea_system_t* sys, size_t cell_index,
                     const char* union_op, const char* inter_op,
                     const char* compl_op);

/**
 * @brief Get cell indices in a universe
 * @return Number of cells found
 */
int alea_cells_in_universe(const alea_system_t* sys, int universe_id,
                                int* out_indices, size_t max_count);

/* ============================================================================
 * SURFACE INFO
 * ============================================================================ */

int alea_surface_get(const alea_system_t* sys, size_t index,
                          int* surface_id, alea_primitive_type_t* type,
                          alea_node_id_t* pos_node, alea_node_id_t* neg_node,
                          alea_boundary_type_t* boundary_type);

int alea_surface_find(const alea_system_t* sys, int surface_id);

/**
 * @brief Get the surface ID at index `idx` in [0, alea_surface_count(sys)).
 * @return The surface ID, or -1 if `sys` is NULL or `idx` is out of range.
 */
int alea_surface_id_at(const alea_system_t* sys, size_t idx);

/**
 * @brief Fill `out_ids` with all surface IDs in the system (storage order).
 *
 * Caller provides a buffer of at least `alea_surface_count(sys)` ints.
 * @return The number of IDs written, or 0 if `sys` or `out_ids` is NULL.
 */
size_t alea_get_surface_ids(const alea_system_t* sys, int* out_ids);

/* ============================================================================
 * UNIVERSE INFO
 * ============================================================================ */

int alea_universe_get(const alea_system_t* sys, size_t index,
                           int* universe_id, size_t* cell_count, alea_bbox_t* bbox);

int alea_universe_find(const alea_system_t* sys, int universe_id);

/* ============================================================================
 * EXTENDED GEOMETRY QUERIES
 * ============================================================================ */

/**
 * @brief Find cell and get both cell ID and material
 */
int alea_find_cell_at(alea_system_t* sys, double x, double y, double z,
                           int* out_cell_id, int* out_material);

/**
 * @brief Enable/disable debug tracing for point queries
 */
void alea_set_debug_trace(int enable);

/* ============================================================================
 * CSG TREE CONSTRUCTION - EXTENDED
 * ============================================================================ */

int alea_get_cell_id(const alea_system_t* sys, int cell_index);

/* ============================================================================
 * RENUMBERING
 * ============================================================================ */

int alea_renumber_cells(alea_system_t* sys, int start_id);
int alea_renumber_surfaces(alea_system_t* sys, int start_id);
int alea_offset_cell_ids(alea_system_t* sys, int offset);
int alea_offset_surface_ids(alea_system_t* sys, int offset);
int alea_offset_material_ids(alea_system_t* sys, int offset);

/* ============================================================================
 * FILTER OPERATIONS
 * ============================================================================ */

size_t alea_get_cells_by_material(const alea_system_t* sys, int material_id,
                                       int* out_indices, size_t max_count);
size_t alea_get_cells_by_universe(const alea_system_t* sys, int universe_id,
                                       int* out_indices, size_t max_count);
size_t alea_get_cells_filling_universe(const alea_system_t* sys, int universe_id,
                                            int* out_indices, size_t max_count);
/* ============================================================================
 * VALIDATION
 * ============================================================================ */

/**
 * @brief Validate system integrity
 * @return 0 if valid, number of issues found otherwise
 */
int alea_validate(const alea_system_t* sys);

/* ============================================================================
 * MACROBODY EXPANSION
 * ============================================================================ */

bool alea_is_macrobody(alea_primitive_type_t type);
alea_node_id_t alea_expand_macrobody(alea_system_t* sys, alea_node_id_t node_id);
alea_node_id_t alea_expand_all_macrobodies(alea_system_t* sys, alea_node_id_t root_id);
int alea_expand_macrobodies_in_system(alea_system_t* sys);

/* ============================================================================
 * CSG NODE INSPECTION
 *
 * Walk the CSG tree starting from a cell's root node. Enables reconstruction
 * of region descriptions like (-S1 & +S2) | -S3 from loaded models.
 *
 * Usage pattern:
 *   alea_node_id_t root;
 *   alea_cell_get(sys, idx, NULL, NULL, NULL, NULL, NULL, &root);
 *   alea_operation_t op = alea_node_operation(sys, root);
 *   if (op == ALEA_OP_PRIMITIVE) {
 *       int sid = alea_node_surface_id(sys, root);
 *       int sense = alea_node_sense(sys, root);
 *   } else {
 *       alea_node_id_t left = alea_node_left(sys, root);
 *       alea_node_id_t right = alea_node_right(sys, root);
 *   }
 * ============================================================================ */

/**
 * @brief Get operation type of a node
 * @return Operation type, or ALEA_OP_PRIMITIVE if node is invalid
 */
alea_operation_t alea_node_operation(const alea_system_t* sys, alea_node_id_t node);

/**
 * @brief Get left child of a boolean node
 * @return Left child node ID, or ALEA_NODE_ID_INVALID if node is a primitive or invalid
 */
alea_node_id_t alea_node_left(const alea_system_t* sys, alea_node_id_t node);

/**
 * @brief Get right child of a boolean node
 * @return Right child node ID, or ALEA_NODE_ID_INVALID if node is primitive/complement/invalid
 */
alea_node_id_t alea_node_right(const alea_system_t* sys, alea_node_id_t node);

/**
 * @brief Get primitive type for a leaf node
 * @return Primitive type, or 0 if not a primitive node
 */
alea_primitive_type_t alea_node_primitive_type(const alea_system_t* sys, alea_node_id_t node);

/**
 * @brief Get primitive ID for a leaf node
 *
 * Two nodes with the same primitive_id reference the same geometric surface
 * (possibly with different sense). Useful for detecting +S / -S pairs.
 *
 * @return Primitive ID, or ALEA_PRIMITIVE_ID_INVALID if not a primitive node
 */
alea_primitive_id_t alea_node_primitive_id(const alea_system_t* sys, alea_node_id_t node);

/**
 * @brief Get primitive data for a leaf node
 * @param out Output primitive data (caller-provided)
 * @return 0 on success, -1 if not a primitive node or invalid
 */
int alea_node_primitive_data(const alea_system_t* sys, alea_node_id_t node,
                                 alea_primitive_data_t* out);

/**
 * @brief Get sense for a primitive node
 * @return +1 (positive/outside) or -1 (negative/inside), or 0 if not a primitive node
 */
int alea_node_sense(const alea_system_t* sys, alea_node_id_t node);

/**
 * @brief Get MCNP surface ID associated with a primitive node
 * @return Surface ID (>0), or 0 if not registered or not a primitive node
 */
int alea_node_surface_id(const alea_system_t* sys, alea_node_id_t node);

#ifdef __cplusplus
}
#endif

#endif /* ALEA_H */
