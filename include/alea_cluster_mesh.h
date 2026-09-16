// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

#ifndef ALEA_CLUSTER_MESH_H
#define ALEA_CLUSTER_MESH_H

#include "alea_cluster_base.h"
#include "alea_mesh.h"

#ifdef __cplusplus
extern "C" {
#endif

/** Collectively sample a structured mesh by contiguous Z slabs. Supports
 * CENTER, CORNERS, SUBCELL, STRATIFIED, and ADAPTIVE modes with explicit or
 * inferred bounds and custom nodes. A nonzero global adaptive sample budget
 * executes slabs in rank order to preserve serial refinement decisions. RAY
 * mode requires max_total_samples=0. Z-directed rays retrace full columns on
 * each rank before clipping contributions to their local slabs.
 * Every rank supplies equivalent configuration except worker counts. Root
 * receives an allocated result and frees it with alea_mesh_result_free();
 * other ranks may pass NULL for root_result. On failure *root_result is NULL. */
alea_cluster_status_t alea_cluster_mesh_sample(
    alea_cluster_t* cluster, alea_system_t* sys,
    const alea_mesh_config_t* config,
    alea_mesh_result_t** root_result);

/** Consume a rank-owned mesh slab. first_z is its global first Z voxel;
 * the borrowed result has local nz and global X/Y/Z coordinates. Return
 * nonzero to cancel collectively. Ranks without Z voxels are not called. */
typedef int (*alea_cluster_mesh_slab_callback_t)(
    size_t first_z, const alea_mesh_result_t* slab, void* user_data);

/** Sample fixed-grid mesh slabs without gathering a complete root mesh.
 * Every rank supplies an equivalent configuration and a callback. Supports
 * the same modes and bounds as alea_cluster_mesh_sample(). Callback side
 * effects on other ranks may occur before a cancellation is reported. */
alea_cluster_status_t alea_cluster_mesh_sample_shards(
    alea_cluster_t* cluster, alea_system_t* sys,
    const alea_mesh_config_t* config,
    alea_cluster_mesh_slab_callback_t callback, void* user_data);

/** Visit voxels on root in global Z/Y/X order while retaining at most one
 * transferred slab in addition to root's own slab. Every rank supplies an
 * equivalent configuration. Root supplies the callback; workers pass NULL.
 * All voxel fields are sampled regardless of config->fields. A callback
 * cancellation returns ALEA_CLUSTER_INTERRUPTED collectively. Earlier
 * callback side effects are not rolled back on failure. */
alea_cluster_status_t alea_cluster_mesh_stream_root(
    alea_cluster_t* cluster, alea_system_t* sys,
    const alea_mesh_config_t* config,
    alea_mesh_voxel_visit_fn callback, void* user_data);

/** Visit sampled voxels on rank zero in global Z/Y/X order. The callback
 * receives all mesh result fields and may cancel by returning nonzero. This
 * path assembles the complete mesh on root before visiting; use rank-owned
 * slabs for meshes that exceed root memory. */
alea_cluster_status_t alea_cluster_mesh_visit_root(
    alea_cluster_t* cluster, alea_system_t* sys,
    const alea_mesh_config_t* config,
    alea_mesh_voxel_visit_fn callback, void* user_data);

#ifdef __cplusplus
}
#endif

#endif /* ALEA_CLUSTER_MESH_H */
