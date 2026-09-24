// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

#ifndef ALEA_VOLUME_INTERNAL_H
#define ALEA_VOLUME_INTERNAL_H

#include "alea.h"

typedef struct {
    size_t path_count;
    double cx, cy, cz, radius;
} alea_volume_problem_t;

#define ALEA_VOLUME_DEFAULT_PARALLEL_SCRATCH_BYTES \
    ((size_t)256 * 1024 * 1024)

int alea_volume_parallel_scratch_layout(
    size_t path_count, int collect_second_moment,
    size_t scratch_limit, size_t requested_workers, size_t ray_count,
    size_t* out_worker_bytes, size_t* out_worker_limit,
    size_t* out_parallel_bytes);

int alea_volume_problem_prepare(alea_system_t* sys,
                                const alea_volume_estimate_options_t* options,
                                alea_volume_problem_t* problem);

int alea_volume_accumulate_ray_range(
    alea_system_t* sys, const alea_volume_problem_t* problem,
    size_t ray_begin, size_t ray_end, alea_rng_algorithm_t rng_algorithm,
    uint64_t seed, size_t requested_workers, size_t max_parallel_scratch_bytes,
    double* sum_l, double* sum_l2,
    size_t* out_actual_workers);

void alea_volume_compute_errors(const double* sum_l, const double* sum_l2,
                                double* rel_errors, size_t count,
                                size_t ray_count);

#endif
