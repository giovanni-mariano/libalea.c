// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

#ifndef ALEA_MESH_CLUSTER_INTERNAL_H
#define ALEA_MESH_CLUSTER_INTERNAL_H

#include "alea_mesh.h"

/* Resolve the same inferred root AABB and padding used by mesh_sample. */
int alea_mesh_cluster_auto_bounds(const alea_system_t* sys,
    const alea_mesh_config_t* config, double bounds[6]);

/* Keep seeded voxel samples tied to global Z indices in a cluster slab. */
alea_mesh_result_t* alea_mesh_sample_with_z_offset(alea_system_t* sys,
    const alea_mesh_config_t* config, int z_index_offset,
    double ray_global_z_min, double ray_global_z_max,
    uint64_t initial_sample_work, uint64_t* final_sample_work);

#endif
