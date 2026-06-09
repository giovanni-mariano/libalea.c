// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

#ifndef ALEA_SPATIAL_HIER_H
#define ALEA_SPATIAL_HIER_H

#include <stddef.h>
#include <stdint.h>
#include "alea_types.h"
#include "core/alea_spatial.h"

typedef struct alea_system alea_system_t;
typedef struct alea_hier_spatial_index alea_hier_spatial_index_t;

typedef struct {
    alea_cell_hit_t hit;
    alea_matrix_t transform;
    int lattice_cell_index;
    alea_matrix_t lattice_transform;
} alea_hier_cell_hit_t;

#define ALEA_HIER_SPATIAL_HIT_CHAIN_MAX 16
#define ALEA_HIER_RAY_PATH_MAX 64

typedef struct {
    alea_spatial_hit_t hit;
    uint8_t ancestor_count;
    uint8_t chain_truncated;
    uint32_t ancestor_cell_indices[ALEA_HIER_SPATIAL_HIT_CHAIN_MAX];
    alea_matrix_t ancestor_transforms[ALEA_HIER_SPATIAL_HIT_CHAIN_MAX];
} alea_hier_spatial_chain_hit_t;

typedef struct {
    uint32_t cell_index;
    int cell_id;
    int material_id;
    int universe_id;
    int fill_universe;
    int depth;
    uint8_t is_lattice;
    int lat_fill_universe;
    double lat_ox;
    double lat_oy;
    double lat_oz;
    alea_matrix_t transform;
} alea_hier_ray_path_entry_t;

typedef struct {
    int count;
    alea_hier_ray_path_entry_t entries[ALEA_HIER_RAY_PATH_MAX];
} alea_hier_ray_path_t;

typedef struct {
    uint32_t cell_index;
    int cell_id;
    double t_enter;
    double t_exit;
} alea_hier_ray_candidate_t;

typedef struct {
    uint32_t placement_index;
    uint32_t parent_cell_index;
    int universe_id;
    int depth;
    uint32_t flags;
    alea_matrix_t transform;
    double t_enter;
    double t_exit;
} alea_hier_placement_ray_candidate_t;

typedef struct {
    size_t universe_count;
    size_t blas_count;
    size_t linear_universe_count;
    size_t blas_cell_count;
    size_t blas_node_count;
    size_t fill_cell_count;
    size_t lattice_cell_count;
    size_t transform_count;
    size_t placement_count;
    size_t root_placement_count;
    size_t fill_placement_count;
    size_t lattice_placement_count;
    int max_placement_depth;
    int max_universe_cells;
    int largest_universe_id;
    size_t memory_bytes;
    size_t point_queries;
} alea_hier_spatial_stats_t;

int alea_hier_spatial_index_build(alea_system_t* sys);
void alea_hier_spatial_index_free(alea_hier_spatial_index_t* idx);
const alea_hier_spatial_stats_t*
alea_hier_spatial_index_stats(const alea_hier_spatial_index_t* idx);
int alea_hier_spatial_find_cells_at_point(alea_system_t* sys,
                                          double x,
                                          double y,
                                          double z,
                                          alea_cell_hit_t* out_hits,
                                          size_t max_hits);
int alea_hier_spatial_find_cells_at_point_uncached(alea_system_t* sys,
                                                   double x,
                                                   double y,
                                                   double z,
                                                   alea_cell_hit_t* out_hits,
                                                   size_t max_hits);
int alea_hier_spatial_find_deepest_cell_at_point(alea_system_t* sys,
                                                 double x,
                                                 double y,
                                                 double z,
                                                 alea_hier_cell_hit_t* out_hit);
int alea_hier_spatial_find_path_at_point(alea_system_t* sys,
                                         double x,
                                         double y,
                                         double z,
                                         alea_hier_cell_hit_t* out_hit,
                                         alea_hier_ray_path_t* out_path);
int alea_hier_spatial_check_path_containment(alea_system_t* sys,
                                             const alea_hier_ray_path_t* path,
                                             uint32_t cell_index,
                                             double x,
                                             double y,
                                             double z,
                                             alea_matrix_t* out_transform,
                                             int* out_lattice_cell_index,
                                             alea_matrix_t* out_lattice_transform);
int alea_hier_spatial_find_path_from_parent(alea_system_t* sys,
                                            const alea_hier_ray_path_t* path,
                                            int parent_entry,
                                            double x,
                                            double y,
                                            double z,
                                            alea_hier_cell_hit_t* out_hit,
                                            alea_hier_ray_path_t* out_path);

/**
 * @brief Reset the hier point-query coherence cache.
 *
 * Mirror of alea_spatial_reset_cache for the hierarchical query path.
 * Called automatically by alea_spatial_reset_cache so existing callers
 * that reset between renders do not need to learn a second API.
 */
void alea_hier_spatial_reset_cache(void);

/**
 * @brief Test whether a world-space point lies inside a cell, using the
 * transform cached by the most recent alea_hier_spatial_find_cells_at_point.
 *
 * Lets raycast's Tier 2 coherence check work for cells with universe_id
 * != 0, where the caller has no other access to the world→local transform.
 *
 * @return 1 if the point is inside the cell,
 *         0 if it is outside,
 *        -1 if the cell is not present in the current cache or the cache
 *           has been invalidated. Caller must fall back to a full lookup.
 */
int alea_hier_spatial_check_cached_containment(alea_system_t* sys,
                                               uint32_t cell_index,
                                               double x, double y, double z);

int alea_hier_spatial_get_cached_cell_state(alea_system_t* sys,
                                            uint32_t cell_index,
                                            double x, double y, double z,
                                            alea_matrix_t* out_transform,
                                            int* out_lattice_cell_index,
                                            alea_matrix_t* out_lattice_transform);

/**
 * @brief Find the cell of `universe_id` containing (lx, ly, lz) in that
 *        universe's local frame, using the per-universe BLAS for pruning.
 *
 * Does NOT descend into fills; returns the cell directly owning the point.
 * Returns the cell index on success, -1 if none, -2 on error (no index).
 */
int alea_hier_spatial_find_cell_in_universe(alea_system_t* sys,
                                            int universe_id,
                                            double lx,
                                            double ly,
                                            double lz);
int alea_hier_spatial_find_ordered_cell_in_universe(alea_system_t* sys,
                                                    int universe_id,
                                                    double lx,
                                                    double ly,
                                                    double lz);
int alea_hier_spatial_query_universe_region(alea_system_t* sys,
                                            int universe_id,
                                            const alea_bbox_t* local_bbox,
                                            alea_spatial_hit_t* out_hits,
                                            size_t max_hits);
int alea_hier_spatial_query_universe_ray(alea_system_t* sys,
                                         int universe_id,
                                         double ox,
                                         double oy,
                                         double oz,
                                         double dx,
                                         double dy,
                                         double dz,
                                         double inv_dx,
                                         double inv_dy,
                                         double inv_dz,
                                         double t_min,
                                         double t_max,
                                         alea_hier_ray_candidate_t* out_hits,
                                         size_t max_hits);
int alea_hier_spatial_query_placements_ray(alea_system_t* sys,
                                           double ox,
                                           double oy,
                                           double oz,
                                           double dx,
                                           double dy,
                                           double dz,
                                           double inv_dx,
                                           double inv_dy,
                                           double inv_dz,
                                           double t_min,
                                           double t_max,
                                           alea_hier_placement_ray_candidate_t* out_hits,
                                           size_t max_hits);
int alea_hier_spatial_check_placement_chain(alea_system_t* sys,
                                            uint32_t placement_index,
                                            double x,
                                            double y,
                                            double z);

/**
 * @brief Return the precomputed per-cell fill transform (forward + inverse)
 *        for `cell_index`, or NULL if the hier index is not built or this
 *        cell has no fill.
 *
 * Cached at hier-index build time so callers in hot loops can skip the
 * per-pixel `alea_matrix_from_mcnp` + `alea_matrix_invert` work. The returned
 * pointer is owned by the hier index and valid until it is destroyed.
 */
const alea_matrix_t*
alea_hier_spatial_get_cell_fill_matrix(const alea_system_t* sys,
                                       uint32_t cell_index);
bool alea_hier_spatial_index_needs_rebuild(const alea_system_t* sys);

/* Debug: total candidates pushed across all universes in the current query. */
void alea_hier_debug_candidates_reset(void);
size_t alea_hier_debug_candidates_get(void);
int alea_hier_spatial_query_region(alea_system_t* sys,
                                   const alea_bbox_t* query_bbox,
                                   alea_spatial_hit_t* out_hits,
                                   size_t max_hits);
int alea_hier_spatial_query_region_direct(alea_system_t* sys,
                                          const alea_bbox_t* query_bbox,
                                          alea_spatial_hit_t* out_hits,
                                          size_t max_hits);
int alea_hier_spatial_query_region_chain(alea_system_t* sys,
                                         const alea_bbox_t* query_bbox,
                                         alea_hier_spatial_chain_hit_t* out_hits,
                                         size_t max_hits);
int alea_hier_spatial_query_slice_z(alea_system_t* sys,
                                    double z,
                                    double x_min,
                                    double x_max,
                                    double y_min,
                                    double y_max,
                                    alea_spatial_hit_t* out_hits,
                                    size_t max_hits);

#endif /* ALEA_SPATIAL_HIER_H */
