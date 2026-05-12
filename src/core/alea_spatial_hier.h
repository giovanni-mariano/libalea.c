// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

#ifndef ALEA_SPATIAL_HIER_H
#define ALEA_SPATIAL_HIER_H

#include <stddef.h>
#include <stdint.h>

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
    int max_universe_cells;
    int largest_universe_id;
    size_t memory_bytes;
} alea_hier_spatial_stats_t;

int alea_hier_spatial_index_build(alea_system_t* sys);
void alea_hier_spatial_index_free(alea_hier_spatial_index_t* idx);
const alea_hier_spatial_stats_t*
alea_hier_spatial_index_stats(const alea_hier_spatial_index_t* idx);

#endif /* ALEA_SPATIAL_HIER_H */
