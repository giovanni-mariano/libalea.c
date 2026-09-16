// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

#ifndef ALEA_CLUSTER_RAYCAST_H
#define ALEA_CLUSTER_RAYCAST_H

#include "alea_cluster_base.h"
#include "alea_raycast.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Collectively trace rays and return the first concrete-cell segment for each
 * input ray, skipping outside/void segments with negative cell IDs.
 * Only rank zero supplies packed XYZ origin/direction arrays, ray_count, t_max,
 * and output arrays sized to ray_count. Other ranks may pass NULL arrays and
 * zero scalar values. A miss sets hit[i] to zero; other output fields for that
 * ray are then unspecified. Results are published on rank zero only. Output
 * may be partially written if an error occurs. Geometry must be equivalent on
 * all ranks. MPI calls stay on the initializing thread.
 */
alea_cluster_status_t alea_cluster_raycast_first_segments(
    alea_cluster_t* cluster, alea_system_t* sys,
    const double* origins_xyz, const double* directions_xyz,
    size_t ray_count, double t_max,
    unsigned char* hit, int32_t* cell_ids,
    double* t_enter, double* t_exit);

/**
 * Collectively trace root-owned packed XYZ rays into the compact batch result.
 * Root supplies ray_count, t_max, options, input arrays, and an allocated result
 * from alea_raycast_batch_result_create(). Other ranks may pass zero/NULL for
 * these arguments, but every rank supplies an equivalent system. The result is
 * replaced on root only after complete success; it is unchanged on failure.
 * All requested local batch fields, including full hierarchy paths, are
 * assembled in input order. The global option limits apply to the final result.
 * The operation uses bounded ray batches and root-owned output storage.
 */
alea_cluster_status_t alea_cluster_raycast_batch(
    alea_cluster_t* cluster, alea_system_t* sys,
    const double* origins_xyz, const double* directions_xyz,
    size_t ray_count, double t_max,
    const alea_raycast_batch_options_t* options,
    alea_raycast_batch_result_t* result);

/** Consume one completed compact ray batch on rank zero. The result covers
 * consecutive input rays beginning at first_ray and is borrowed until this
 * callback returns. Its CSR offsets start at zero. Return nonzero to stop the
 * collective operation with ALEA_CLUSTER_INTERRUPTED. Do not call another
 * cluster collective from inside the callback. */
typedef int (*alea_cluster_ray_batch_callback_t)(
    size_t first_ray, const alea_raycast_batch_result_t* batch,
    void* user_data);

/** Collectively trace root-owned rays and stream completed batches to rank
 * zero. Only root supplies ray inputs, options, callback, and user_data; other
 * ranks may pass zero/NULL. A callback is invoked once per nonempty batch, in
 * input order, with at most 1024 rays. Empty input invokes no callback.
 * Segment, path-entry, and output-byte limits apply to each assembled batch,
 * rather than to the whole stream. Side effects of earlier callbacks are not
 * rolled back if a later batch fails. */
alea_cluster_status_t alea_cluster_raycast_batch_stream(
    alea_cluster_t* cluster, alea_system_t* sys,
    const double* origins_xyz, const double* directions_xyz,
    size_t ray_count, double t_max,
    const alea_raycast_batch_options_t* options,
    alea_cluster_ray_batch_callback_t callback, void* user_data);

/** Consume a rank-owned contiguous shard of one ray batch. Called on every
 * rank with work, with a borrowed local result and global first-ray index.
 * Return nonzero to stop the collective operation. Callback side effects on
 * other ranks may already have occurred when cancellation is reported. */
typedef int (*alea_cluster_ray_shard_callback_t)(
    size_t first_ray, const alea_raycast_batch_result_t* shard,
    void* user_data);

/** Collectively trace root-owned rays in bounded batches, retaining each
 * shard only on its computing rank. Root supplies the rays, count, t_max, and
 * options; workers may pass NULL/zero. Every rank supplies its callback and
 * user_data. Limits in options apply to each rank's shard in each batch.
 * Callbacks may write rank-specific output files. No full result is gathered. */
alea_cluster_status_t alea_cluster_raycast_batch_shards(
    alea_cluster_t* cluster, alea_system_t* sys,
    const double* origins_xyz, const double* directions_xyz,
    size_t ray_count, double t_max,
    const alea_raycast_batch_options_t* options,
    alea_cluster_ray_shard_callback_t callback, void* user_data);

/** Consume rank-owned complete-coverage rows. first_row is the global input
 * row index; the borrowed result remains valid until the callback returns.
 * Return nonzero to cancel the operation collectively. */
typedef int (*alea_cluster_coverage_shard_callback_t)(
    size_t first_row, const alea_ray_coverage_slice_result_t* shard,
    void* user_data);

/** Collectively classify root-owned rays into complete-coverage rows. Fixed
 * input rows are assigned in contiguous shards, with packed interval/owner
 * arrays retained on their computing ranks. Root supplies rays, row count,
 * options and optional provenance arrays; every rank supplies a callback.
 * Adaptive refinement is not supported. Resource limits apply per shard.
 * A nonzero callback return cancels collectively. */
alea_cluster_status_t alea_cluster_coverage_shards(
    alea_cluster_t* cluster, alea_system_t* sys,
    const double* origins_xyz, const double* directions_xyz,
    size_t row_count, const uint8_t* direction_tags,
    const double* transverse_coordinates,
    const alea_ray_coverage_slice_options_t* options,
    alea_cluster_coverage_shard_callback_t callback, void* user_data);

/** Gather one bounded coverage batch to rank zero and invoke a callback with
 * globally ordered rows. CSR offsets are rebased to zero for each batch.
 * Root supplies the callback; workers may pass NULL. Resource limits apply
 * to each assembled batch, not the entire stream. */
alea_cluster_status_t alea_cluster_coverage_stream(
    alea_cluster_t* cluster, alea_system_t* sys,
    const double* origins_xyz, const double* directions_xyz,
    size_t row_count, const uint8_t* direction_tags,
    const double* transverse_coordinates,
    const alea_ray_coverage_slice_options_t* options,
    alea_cluster_coverage_shard_callback_t callback, void* user_data);

/** Assemble all fixed complete-coverage rows into one root-owned result.
 * Root supplies inputs, options, and a result created with
 * alea_ray_coverage_slice_result_create(); workers may pass NULL/zero for
 * these arguments. The result is replaced only on success. Resource limits
 * apply to the entire result, and root memory grows with total output size.
 * Adaptive refinement is not supported. */
alea_cluster_status_t alea_cluster_coverage(
    alea_cluster_t* cluster, alea_system_t* sys,
    const double* origins_xyz, const double* directions_xyz,
    size_t row_count, const uint8_t* direction_tags,
    const double* transverse_coordinates,
    const alea_ray_coverage_slice_options_t* options,
    alea_ray_coverage_slice_result_t* root_result);

#ifdef __cplusplus
}
#endif

#endif /* ALEA_CLUSTER_RAYCAST_H */
