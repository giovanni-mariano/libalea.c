// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

#ifndef ALEA_CLUSTER_H
#define ALEA_CLUSTER_H

#include "alea.h"
#include "alea_raycast.h"
#include "alea_render.h"
#include "alea_slice.h"
#include "alea_mesh.h"
#include "alea_geo_validator.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct alea_cluster alea_cluster_t;

typedef enum {
    ALEA_CLUSTER_OK = 0,
    ALEA_CLUSTER_INVALID_ARGUMENT,
    ALEA_CLUSTER_INVALID_STATE,
    ALEA_CLUSTER_OUT_OF_MEMORY,
    ALEA_CLUSTER_BACKEND_ERROR,
    ALEA_CLUSTER_MODEL_MISMATCH,
    ALEA_CLUSTER_COMPUTE_ERROR,
    ALEA_CLUSTER_INTERRUPTED,
    ALEA_CLUSTER_IO_ERROR,
    ALEA_CLUSTER_OUTPUT_LIMIT
} alea_cluster_status_t;

typedef struct {
    alea_volume_estimate_stats_t volume;
    int rank_count;
    size_t local_rays_completed;
    size_t local_workers;
} alea_cluster_volume_stats_t;

/**
 * Initialize the process-wide cluster runtime.
 *
 * In an MPI build this collectively attaches to MPI_COMM_WORLD, initializing
 * MPI with MPI_THREAD_FUNNELED when needed. MPI calls must be made by the
 * thread that calls this function. A local build creates a one-rank runtime.
 */
alea_cluster_status_t alea_cluster_initialize(int* argc, char*** argv);

/**
 * Finalize a runtime initialized by Alea. All cluster contexts must first be
 * destroyed. MPI is finalized only when Alea initialized it.
 */
alea_cluster_status_t alea_cluster_finalize(void);

/** Create a cluster context. V1 permits one live context per process. */
alea_cluster_t* alea_cluster_create(void);
void alea_cluster_destroy(alea_cluster_t* cluster);

int alea_cluster_rank(const alea_cluster_t* cluster);
int alea_cluster_size(const alea_cluster_t* cluster);
int alea_cluster_is_root(const alea_cluster_t* cluster);
const char* alea_cluster_backend(const alea_cluster_t* cluster);

/**
 * Collectively read a file on rank zero and copy its bytes to every rank.
 * Only rank zero uses path. Each rank receives a malloc-allocated, NUL-terminated
 * buffer in *data (free it with free()), and its byte count in *length. Binary
 * input is supported; the terminator is not counted. Outputs are reset on error.
 * All ranks must call this function in the same order.
 */
alea_cluster_status_t alea_cluster_read_file(
    alea_cluster_t* cluster, const char* path, char** data, size_t* length);

/** Collectively read an MCNP input, expanding READ FILE= cards on rank zero.
 * Paths in nested cards are resolved relative to the file containing the card.
 * The returned malloc-allocated, NUL-terminated text is identical on every
 * rank and can be passed to mcnp_load_string(); free it with free(). Only root
 * uses path. Missing files return ALEA_CLUSTER_IO_ERROR; malformed cards or
 * nesting deeper than 32 files return ALEA_CLUSTER_COMPUTE_ERROR. */
alea_cluster_status_t alea_cluster_read_mcnp_input(
    alea_cluster_t* cluster, const char* path, char** data, size_t* length);

/** Collectively agree on a local status before the next cluster operation. */
alea_cluster_status_t alea_cluster_agree(
    alea_cluster_t* cluster, alea_cluster_status_t local_status);

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

/** Collectively render statically assigned image tiles. Every rank supplies
 * an equivalent system, render configuration, and prepared camera. Only rank
 * zero supplies a framebuffer with cfg width/height; workers may pass NULL.
 * The framebuffer receives color and cell IDs, plus material IDs, depth and
 * normals when all three auxiliary arrays are present. On success, edges are
 * darkened on root when cfg->edges is set. On failure the framebuffer may be
 * partially written. MPI calls stay on the initializing thread. */
alea_cluster_status_t alea_cluster_render_scene(
    alea_cluster_t* cluster, alea_system_t* sys,
    const render_config_t* cfg, const render_camera_t* cam,
    render_framebuffer_t* root_framebuffer);

/** Collectively trace and rasterize centered rows of one slice view. Every
 * rank supplies an equivalent system and view. Root supplies output and
 * options; other ranks may pass NULL. Requested raster arrays are filled on
 * root in global row order. Trace segment and byte limits apply to all rows
 * combined. Output arrays may be partially written if a gather fails. */
alea_cluster_status_t alea_cluster_slice_raster(
    alea_cluster_t* cluster, alea_system_t* sys,
    const alea_slice_view_t* view,
    const alea_slice_raster_options_t* options,
    alea_slice_raster_t* root_output);

/** Consume one completed raster plane on rank zero. The raster uses the
 * caller's reusable output arrays; copy any data that must outlive the next
 * plane. Return nonzero to cancel the stack collectively. */
typedef int (*alea_cluster_slice_plane_callback_t)(
    size_t plane_index, const alea_slice_raster_t* raster, void* user_data);

/** Process an offline stack of slice views in order. Every rank supplies the
 * same view_count and equivalent views. Root supplies one output descriptor,
 * options, and callback; the descriptor's dimensions and fields are reused
 * for every plane. Limits apply to each plane. Earlier callbacks are not
 * rolled back if a later plane fails or the callback cancels. */
alea_cluster_status_t alea_cluster_slice_stack_stream(
    alea_cluster_t* cluster, alea_system_t* sys,
    const alea_slice_view_t* views, size_t view_count,
    const alea_slice_raster_options_t* options,
    alea_slice_raster_t* root_output,
    alea_cluster_slice_plane_callback_t callback, void* user_data);

/** Collectively sample a structured mesh by contiguous Z slabs. Supports
 * CENTER, CORNERS, SUBCELL, STRATIFIED, and ADAPTIVE modes with explicit or
 * inferred bounds and custom nodes. ADAPTIVE and horizontal X/Y RAY modes
 * require max_total_samples=0. Z-directed rays are not supported.
 * Every rank supplies equivalent configuration except workers; callbacks and
 * Z-directed ray sampling is not accepted. Root receives an allocated
 * result and frees it with alea_mesh_result_free(); other ranks may pass NULL
 * for root_result. On failure *root_result is NULL. */
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

/** Visit sampled voxels on rank zero in global Z/Y/X order. The callback
 * receives all mesh result fields and may cancel by returning nonzero. This
 * path assembles the complete mesh on root before visiting; use rank-owned
 * slabs for meshes that exceed root memory. */
alea_cluster_status_t alea_cluster_mesh_visit_root(
    alea_cluster_t* cluster, alea_system_t* sys,
    const alea_mesh_config_t* config,
    alea_mesh_voxel_visit_fn callback, void* user_data);

/** Collectively validate the existing seeded random-ray sequence. Every rank
 * supplies equivalent options and geometry; root supplies an initialized
 * result. Root merges per-ray findings in serial ray order, including global
 * signature sampling and truncation. Intermediate findings are capped at 4096
 * per ray; a ray exceeding that bound returns ALEA_CLUSTER_OUTPUT_LIMIT rather
 * than silently dropping evidence. The root result may contain a committed
 * prefix on failure. */
alea_cluster_status_t alea_cluster_validate_geometry(
    alea_cluster_t* cluster, alea_system_t* sys,
    const alea_geom_validator_options_t* options,
    alea_geom_validator_result_t* root_result);

/**
 * Collectively estimate physical volumes on every rank.
 *
 * Every rank must call this function in the same order with equivalent
 * systems and identical scalar options except requested_workers, which may
 * differ by rank. Output arrays are sized to
 * alea_volume_path_count(sys), and receive the final result on every rank.
 * An explicit sampling sphere in options is used identically on every rank;
 * the caller is responsible for its coverage. The progress callback, when
 * supplied, is invoked on rank zero only.
 */
alea_cluster_status_t alea_cluster_estimate_volumes(
    alea_cluster_t* cluster,
    alea_system_t* sys,
    const alea_volume_estimate_options_t* options,
    double* volumes,
    double* rel_errors,
    alea_cluster_volume_stats_t* out_stats);

const char* alea_cluster_status_string(alea_cluster_status_t status);

#ifdef __cplusplus
}
#endif

#endif /* ALEA_CLUSTER_H */
