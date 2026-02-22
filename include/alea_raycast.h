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

/* ============================================================================
 * RAYCAST FUNCTIONS
 * ============================================================================ */

/**
 * @brief Create a new raycast result object
 * @return New result object (caller must free with alea_raycast_result_free)
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
int alea_raycast(const alea_system_t* sys,
                     double ox, double oy, double oz,
                     double dx, double dy, double dz,
                     double t_max,
                     alea_raycast_result_t* result);

/**
 * @brief Cell-aware raycast using per-cell surface index
 *
 * More efficient than global surface testing. Tracks through cells
 * one at a time, testing only surfaces belonging to each cell.
 */
int alea_raycast_cell_aware(const alea_system_t* sys,
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
int alea_ray_first_cell(const alea_system_t* sys,
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
 * @return 0 on success, -1 on invalid index
 */
int alea_raycast_segment_get(const alea_raycast_result_t* result, size_t index,
                                 double* t_enter, double* t_exit,
                                 int* cell_id, int* material_id, double* density);

/**
 * @brief Calculate total path length through material
 *
 * @param result Raycast result
 * @param material_id Material to sum (-1 = all materials)
 * @return Total path length
 */
double alea_raycast_path_length(const alea_raycast_result_t* result, int material_id);

/**
 * @brief Free raycast result memory
 */
void alea_raycast_result_free(alea_raycast_result_t* result);

#ifdef __cplusplus
}
#endif

#endif /* ALEA_RAYCAST_H */
