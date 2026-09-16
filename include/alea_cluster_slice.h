// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

#ifndef ALEA_CLUSTER_SLICE_H
#define ALEA_CLUSTER_SLICE_H

#include "alea_cluster_base.h"
#include "alea_slice.h"

#ifdef __cplusplus
extern "C" {
#endif

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

#ifdef __cplusplus
}
#endif

#endif /* ALEA_CLUSTER_SLICE_H */
