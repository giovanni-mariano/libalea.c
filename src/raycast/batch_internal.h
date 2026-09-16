// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

#ifndef ALEA_RAYCAST_BATCH_INTERNAL_H
#define ALEA_RAYCAST_BATCH_INTERNAL_H

#include "raycast.h"

/* Shared only by the compact batch producer and the optional cluster assembler.
 * Applications continue to use the opaque public result and its accessors. */
struct alea_raycast_batch_result {
    size_t ray_count;
    size_t segment_count;
    uint32_t fields;
    uint64_t* ray_offsets;
    double* t_enter;
    double* t_exit;
    int32_t* cell_ids;
    int32_t* material_ids;
    double* densities;
    int32_t* enter_surface_ids;
    int32_t* exit_surface_ids;
    uint8_t* resolution_flags;
    int32_t* projected_cell_ids;
    int32_t* projected_material_ids;
    int32_t* projected_universe_ids;
    int32_t* projected_fill_universes;
    int32_t* projected_depths;
    uint8_t* projected_is_lattice;
    uint64_t* projected_occurrence_keys;
    size_t path_entry_count;
    uint64_t* segment_path_offsets;
    int32_t* path_cell_ids;
    int32_t* path_material_ids;
    int32_t* path_universe_ids;
    int32_t* path_fill_universes;
    int32_t* path_depths;
    uint8_t* path_is_lattice;
    double* path_lattice_origins_xyz;
    uint64_t* path_occurrence_keys;
    struct {
        uint8_t valid;
        uint64_t system_id;
        uint64_t geometry_generation;
        double origin[3];
        double u_axis[3];
        double v_axis[3];
        double u_min, u_max, v_min, v_max;
        size_t row_count;
        int projected_depth;
    } fast_slice_cache;
    alea_raycast_batch_work_stats_t work_stats;
};

/* Steal allocated arrays from source after freeing destination's old arrays. */
void alea_raycast_batch_result_replace_internal(
    alea_raycast_batch_result_t* destination,
    alea_raycast_batch_result_t* source);

#endif
