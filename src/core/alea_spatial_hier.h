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
    size_t point_blas_queries;
    size_t point_blas_node_visits;
    size_t point_bbox_tests;
    size_t point_candidates;
    size_t point_exact_tests;
    size_t point_linear_scans;
    size_t point_linear_cell_tests;
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
int alea_hier_spatial_find_deepest_cell_at_point(alea_system_t* sys,
                                                 double x,
                                                 double y,
                                                 double z,
                                                 alea_hier_cell_hit_t* out_hit);

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
int alea_hier_spatial_query_region(alea_system_t* sys,
                                   const alea_bbox_t* query_bbox,
                                   alea_spatial_hit_t* out_hits,
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
