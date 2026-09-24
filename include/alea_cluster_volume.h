// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

#ifndef ALEA_CLUSTER_VOLUME_H
#define ALEA_CLUSTER_VOLUME_H

#include "alea_cluster_base.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    alea_volume_estimate_stats_t volume;
    int rank_count;
    size_t local_rays_completed;
    size_t local_workers;
} alea_cluster_volume_stats_t;

/** Collectively estimate physical volumes on every rank. Systems and scalar
 * options must match, except requested_workers and max_parallel_scratch_bytes
 * may differ by rank. Arrays
 * have alea_volume_path_count(sys) entries. Root alone receives progress
 * callbacks; final volumes and uncertainties are replicated. */
alea_cluster_status_t alea_cluster_estimate_volumes(
    alea_cluster_t* cluster, alea_system_t* sys,
    const alea_volume_estimate_options_t* options,
    double* volumes, double* rel_errors,
    alea_cluster_volume_stats_t* out_stats);

#ifdef __cplusplus
}
#endif

#endif /* ALEA_CLUSTER_VOLUME_H */
