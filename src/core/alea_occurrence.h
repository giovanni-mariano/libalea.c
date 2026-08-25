// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
// SPDX-License-Identifier: MPL-2.0

#ifndef ALEA_OCCURRENCE_H
#define ALEA_OCCURRENCE_H

#include <stdint.h>

#define ALEA_OCCURRENCE_ROOT UINT64_C(1469598103934665603)

static inline uint64_t alea_occurrence_mix(uint64_t key, uint64_t value) {
    key ^= value;
    return key * UINT64_C(1099511628211);
}

static inline uint64_t alea_occurrence_cell(uint64_t universe_occurrence,
                                            uint32_t cell_index) {
    return alea_occurrence_mix(
        universe_occurrence, ((uint64_t)cell_index << 1) | UINT64_C(1));
}

static inline uint64_t alea_occurrence_lattice(uint64_t cell_occurrence,
                                               int i, int j, int k) {
    return alea_occurrence_mix(
        alea_occurrence_mix(cell_occurrence, (uint32_t)i),
        ((uint64_t)(uint32_t)j << 32) | (uint32_t)k);
}

typedef struct {
    uint64_t universe_key;
    uint64_t last_cell_key;
} alea_occurrence_state_t;

static inline void alea_occurrence_state_init(alea_occurrence_state_t* state) {
    state->universe_key = ALEA_OCCURRENCE_ROOT;
    state->last_cell_key = 0;
}

static inline uint64_t alea_occurrence_state_step(
    alea_occurrence_state_t* state, uint32_t cell_index,
    int is_lattice, int i, int j, int k) {
    const uint64_t cell_key = alea_occurrence_cell(
        state->universe_key, cell_index);
    state->last_cell_key = cell_key;
    state->universe_key = is_lattice
        ? alea_occurrence_lattice(cell_key, i, j, k)
        : cell_key;
    return cell_key;
}

#endif
