// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
// SPDX-License-Identifier: MPL-2.0

#ifndef ALEA_TRANSITION_VALIDATION_H
#define ALEA_TRANSITION_VALIDATION_H

#include "alea_geo_validator.h"

/* Internal reusable storage for the exact-coverage fallback.  One logical
 * worker owns one workspace; it is not safe to use concurrently. */
typedef struct alea_transition_workspace {
    alea_cell_hit_t* hits;
    uint64_t* keys;
    uint64_t* parents;
    uint8_t* owners;
    size_t* child_counts;
    size_t capacity;
} alea_transition_workspace_t;

void alea_transition_workspace_init(alea_transition_workspace_t* workspace);
void alea_transition_workspace_free(alea_transition_workspace_t* workspace);

int alea_check_transition_local_reuse(
    alea_system_t* sys, int universe_id, int current_cell_id,
    int primary_surface_id, const int* tied_surface_ids,
    size_t tied_surface_count, const double point[3],
    const double direction[3], const alea_transition_options_t* options,
    alea_transition_result_t* result,
    alea_transition_workspace_t* workspace);

#endif
