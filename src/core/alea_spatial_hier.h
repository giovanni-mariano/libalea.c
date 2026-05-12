// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

#ifndef ALEA_SPATIAL_HIER_H
#define ALEA_SPATIAL_HIER_H

#include <stddef.h>
#include <stdint.h>
#include "alea_types.h"

typedef struct alea_system alea_system_t;
typedef struct alea_hier_spatial_index alea_hier_spatial_index_t;

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

#endif /* ALEA_SPATIAL_HIER_H */
