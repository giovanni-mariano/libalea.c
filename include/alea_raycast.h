// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/**
 * @file alea_raycast.h
 * @brief Alea Raycast Module API
 *
 * Optional raycast functionality. Requires linking with libalea_raycast.a
 */

#ifndef ALEA_RAYCAST_H
#define ALEA_RAYCAST_H

#include "alea.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ============================================================================
 * RAYCAST TYPES
 * ============================================================================ */

/** Opaque raycast result type */
typedef struct alea_raycast_result alea_raycast_result_t;

/** One resolved hierarchy entry attached to an opt-in ray segment path.
 *
 * ``occurrence_key`` identifies this concrete occurrence, including its
 * ancestor placement/lattice context.  It is stable for equivalent paths
 * traced in opposite directions through an unchanged system.
 */
typedef struct {
    uint32_t cell_index;
    int cell_id;
    int material_id;
    int universe_id;
    int fill_universe;
    int depth;
    uint8_t is_lattice;
    double lattice_origin[3];
    uint64_t occurrence_key;
} alea_raycast_path_entry_t;

/* ============================================================================
 * RAYCAST FUNCTIONS
 * ============================================================================ */

/**
 * @brief Create a new raycast result object
 * @return New result object (caller must destroy with alea_raycast_result_destroy)
 */
alea_raycast_result_t* alea_raycast_result_create(void);

/**
 * @brief Cast a ray through geometry and find all cell intersections
 *
 * @param sys System
 * @param ox, oy, oz Ray origin
 * @param dx, dy, dz Ray direction (will be normalized)
 * @param t_max Maximum ray distance (0 = infinite)
 * @param result Output result
 * @return 0 on success, -1 on error
 */
int alea_raycast(alea_system_t* sys,
                     double ox, double oy, double oz,
                     double dx, double dy, double dz,
                     double t_max,
                     alea_raycast_result_t* result);

/**
 * @brief Cell-aware raycast using per-cell surface index
 *
 * Semantic equivalent of alea_raycast(). For non-lattice models this can
 * track through cells one at a time using the per-cell surface index. For
 * lattice models it uses the canonical lattice-aware pipeline so synthetic
 * DDA element-boundary hits are included.
 */
int alea_raycast_cell_aware(alea_system_t* sys,
                                double ox, double oy, double oz,
                                double dx, double dy, double dz,
                                double t_max,
                                alea_raycast_result_t* result);

/**
 * @brief Fast hierarchical raycast for material/path segments.
 *
 * Uses the hierarchical spatial index and per-ray path state to step through
 * cells without building a flat spatial index or collecting global surface
 * hits. The returned segments contain boundary surface IDs where available,
 * but the result's full hit list is not part of this function's contract.
 */
int alea_raycast_hier_fast_segments(alea_system_t* sys,
                                    double ox, double oy, double oz,
                                    double dx, double dy, double dz,
                                    double t_max,
                                    alea_raycast_result_t* result);

/**
 * @brief Find first cell along ray
 *
 * @param sys System
 * @param ox, oy, oz Ray origin
 * @param dx, dy, dz Ray direction
 * @param t_max Maximum distance (0 = infinite)
 * @param out_t Output: distance to first hit (can be NULL)
 * @return Cell ID of first cell hit, or -1 if none
 */
int alea_ray_first_cell(alea_system_t* sys,
                            double ox, double oy, double oz,
                            double dx, double dy, double dz,
                            double t_max, double* out_t);

/**
 * @brief Get number of segments in raycast result
 */
size_t alea_raycast_segment_count(const alea_raycast_result_t* result);

/**
 * @brief Get segment data
 * @param result Raycast result
 * @param index Segment index
 * @param t_enter Output: entry distance
 * @param t_exit Output: exit distance
 * @param cell_id Output: cell ID (-1 for void)
 * @param material_id Output: material ID
 * @param density Output: material density
 * @param enter_surface_id Output: boundary crossed at t_enter (-1 none, 0 synthetic, >0 physical)
 * @param exit_surface_id Output: boundary crossed at t_exit (-1 none, 0 synthetic, >0 physical)
 * @return 0 on success, -1 on invalid index
 */
int alea_raycast_segment_get(const alea_raycast_result_t* result, size_t index,
                                 double* t_enter, double* t_exit,
                                 int* cell_id, int* material_id, double* density,
                                 int* enter_surface_id, int* exit_surface_id);

/** Enable or disable hierarchy-path capture for subsequent hierarchical traces.
 * Disabled by default.  The setting survives result-buffer reuse. */
void alea_raycast_result_set_path_capture(alea_raycast_result_t* result,
                                          int enabled);

/** Return the number of hierarchy entries attached to a segment. */
size_t alea_raycast_segment_path_count(const alea_raycast_result_t* result,
                                       size_t segment_index);

/** Copy one hierarchy path entry attached to a segment. */
int alea_raycast_segment_path_get(const alea_raycast_result_t* result,
                                  size_t segment_index,
                                  size_t path_entry_index,
                                  alea_raycast_path_entry_t* out_entry);

/**
 * @brief Calculate total path length through material
 *
 * @param result Raycast result
 * @param material_id Material to sum (-1 = all materials)
 * @return Total path length
 */
double alea_raycast_path_length(const alea_raycast_result_t* result, int material_id);

/**
 * @brief Free raycast result internal buffers (for stack-allocated results)
 */
void alea_raycast_result_free(alea_raycast_result_t* result);

/**
 * @brief Destroy a heap-allocated raycast result (from alea_raycast_result_create)
 */
void alea_raycast_result_destroy(alea_raycast_result_t* result);

#ifdef __cplusplus
}
#endif

#endif /* ALEA_RAYCAST_H */
