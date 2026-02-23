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
 *   alea_system_t* sys = alea_load_mcnp("input.inp");
 *   if (!sys) {
 *       fprintf(stderr, "Error: %s\n", alea_error());
 *       return 1;
 *   }
 *
 *   alea_build_universe_index(sys);
 *   int cell = alea_find_cell(sys, 0.0, 0.0, 0.0);
 *
 *   alea_export_mcnp(sys, "output.inp");
 *   alea_destroy(sys);
 */

#ifndef ALEA_H
#define ALEA_H

#include <stddef.h>
#include <stdint.h>
#include <stdbool.h>
#include <stdio.h>
#include "alea_types.h"

/* alea_cell_hit_t is defined in core/alea_universe.h (included via alea_system.h) */
/* void_result_t is forward-declared in alea_types.h */


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
 *   result = alea_load_mcnp(filename);
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
 * ============================================================================ */

alea_system_t* alea_load_mcnp(const char* filename);
alea_system_t* alea_load_mcnp_string(const char* input, size_t length);
alea_system_t* alea_load_openmc(const char* filename);

/* ============================================================================
 * INDEXING
 * ============================================================================ */

int alea_build_universe_index(alea_system_t* sys);
int alea_build_spatial_index(alea_system_t* sys);

/* ============================================================================
 * GEOMETRY QUERIES
 * ============================================================================ */

int alea_find_cell(const alea_system_t* sys, double x, double y, double z);
int alea_find_all_cells(const alea_system_t* sys, double x, double y, double z,
                            alea_cell_hit_t* hits, size_t max_hits);
bool alea_point_inside(const alea_system_t* sys, alea_node_id_t node,
                           double x, double y, double z);
int alea_material_at(const alea_system_t* sys, double x, double y, double z);
int alea_find_overlaps(const alea_system_t* sys, int* pairs, size_t max_pairs);

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
 * @brief Access surface entry by index (stable after vector growth)
 */
#define alea_surface_at(sys, idx) (&(sys)->surfaces.data[idx])

/**
 * @brief Create a plane surface with automatic registration
 * @param surface_id MCNP-style surface ID (will be registered)
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
                      int material, double density, int universe);
int alea_set_fill(alea_system_t* sys, int cell_index, int fill_universe, int transform);

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
 * @brief Flatten and simplify all cells (full optimization)
 *
 * Applies NNF conversion, flattening, balancing, contradiction detection,
 * and empty cell removal. Cells that simplify to empty are removed.
 *
 * @param sys CSG system
 * @param stats Optional accumulated stats output (can be NULL)
 */
void alea_flatten_all_cells(alea_system_t* sys, alea_simplify_stats_t* stats);

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
int alea_compute_bounding_sphere(const alea_system_t* sys,
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
int alea_estimate_cell_volumes(const alea_system_t* sys,
                                   double ox, double oy, double oz,
                                   double radius, int n_rays,
                                   double* volumes,
                                   double* rel_errors);

/**
 * @brief Estimate volumes per cell instance (spatial-index aware)
 *
 * Like alea_estimate_cell_volumes but resolves each ray segment to
 * a specific cell instance via the spatial index. Distinguishes the same
 * cell appearing in different fill contexts.
 *
 * If @p rel_errors is non-NULL, it receives the 1-sigma relative statistical
 * error for each instance: sigma_V / V.  Set to -1 for zero-volume instances.
 *
 * Requires the raycast module to be linked; returns -1 otherwise.
 *
 * @param sys        CSG system (must have spatial index built)
 * @param n_rays     Number of random rays
 * @param volumes    Output array of size alea_spatial_index_instance_count(sys)
 * @param rel_errors Output array of same size, or NULL
 * @return 0 on success, -1 on error
 */
int alea_estimate_instance_volumes(const alea_system_t* sys,
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
 * All export options are read from sys->config (surface_policy, dedup,
 * export_materials, export_transforms, universe_depth, fill_depth).
 * ============================================================================ */

int alea_export_mcnp(const alea_system_t* sys, const char* filename);
int alea_export_mcnp_stream(const alea_system_t* sys, FILE* out);
int alea_export_openmc(const alea_system_t* sys, const char* filename);
int alea_export_stream(const alea_system_t* sys, alea_export_format_t format, FILE* out);

/* ============================================================================
 * VOID GENERATION
 *
 * Void options are read from sys->config (void_max_depth, void_min_size,
 * void_samples, merge_cell_weight, merge_surface_weight, etc.).
 * ============================================================================ */

void_result_t* alea_void_generate(alea_system_t* sys, const alea_bbox_t* bounds);
int alea_void_add_cells(alea_system_t* sys, void_result_t* result);
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

/* Log level: 0=none, 1=error, 2=warn, 3=info, 4=debug, 5=trace */
#ifndef ALEA_LOG_LEVEL_DEFINED
#define ALEA_LOG_LEVEL_DEFINED
typedef enum {
    ALEA_LOG_LEVEL_NONE  = 0,
    ALEA_LOG_LEVEL_ERROR = 1,
    ALEA_LOG_LEVEL_WARN  = 2,
    ALEA_LOG_LEVEL_INFO  = 3,
    ALEA_LOG_LEVEL_DEBUG = 4,
    ALEA_LOG_LEVEL_TRACE = 5
} alea_log_level_t;
#endif
void alea_log_set_level(alea_log_level_t level);
alea_log_level_t alea_log_get_level(void);

/* ============================================================================
 * MATERIAL OPERATIONS
 * ============================================================================ */

/**
 * @brief Create a mixture of materials
 *
 * Creates a new material as a weighted mixture of existing materials.
 *
 * @param sys System
 * @param mat_ids Array of material IDs to mix
 * @param fractions Weight fractions
 * @param count Number of materials to mix
 * @param new_mat_id ID for the new mixture (0 for auto-assign)
 * @return Assigned material ID on success, -1 on error
 */
int alea_create_mixture(alea_system_t* sys, const int* mat_ids,
                            const double* fractions, size_t count, int new_mat_id);

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
int alea_find_cell_at(const alea_system_t* sys, double x, double y, double z,
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
size_t alea_spatial_index_instance_count(const alea_system_t* sys);

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
