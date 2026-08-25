// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
// SPDX-License-Identifier: MPL-2.0

#ifndef ALEA_TRANSITION_SLICE_CRITICAL_H
#define ALEA_TRANSITION_SLICE_CRITICAL_H

#include "alea_geo_validator.h"

typedef int (*alea_transition_slice_critical_finding_sink_t)(
    const alea_transition_slice_critical_finding_t* finding, void* userdata);

int alea_transition_slice_enumerate_critical_tiles(
    alea_system_t* sys,
    const alea_slice_view_t* view,
    const alea_transition_slice_options_t* options,
    const alea_transition_slice_critical_tile_t* tiles,
    size_t tile_count,
    alea_transition_slice_critical_finding_sink_t finding_sink,
    void* finding_sink_userdata,
    alea_transition_slice_stats_t* stats);

#endif
