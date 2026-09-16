// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

#ifndef ALEA_CLUSTER_VALIDATOR_H
#define ALEA_CLUSTER_VALIDATOR_H

#include "alea_cluster_base.h"
#include "alea_geo_validator.h"
#include "alea_slice.h"

#ifdef __cplusplus
extern "C" {
#endif

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

/** Collectively validate caller-provided slice curves in global curve order.
 * Every rank supplies equivalent geometry, view, curves, and options; root
 * supplies an initialized result. Findings retain original curve indices and
 * event order. Root applies signature sampling and global crossing/error
 * limits, rerunning a boundary curve when necessary. A curve producing more
 * than 4096 intermediate findings returns ALEA_CLUSTER_OUTPUT_LIMIT. */
alea_cluster_status_t alea_cluster_validate_slice_curves(
    alea_cluster_t* cluster, alea_system_t* sys,
    const alea_slice_view_t* view,
    const alea_slice_curves_t* curves,
    const alea_geom_validator_options_t* options,
    alea_geom_validator_result_t* root_result);

#ifdef __cplusplus
}
#endif

#endif /* ALEA_CLUSTER_VALIDATOR_H */
