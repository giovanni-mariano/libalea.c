// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/** @file ray_coverage_api.c
 * @brief Installed compact complete-coverage API adapters. */

#include "alea.h"
#include "alea_raycast.h"
#include "raycast.h"
#include "core/alea_system.h"
#include <limits.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#ifdef _OPENMP
#include <omp.h>
#endif

static int coverage_validate_ray_inputs(const double* origins_xyz,
                                        const double* directions_xyz,
                                        size_t row_count) {
    if (row_count == 0) return 0;
    if (!origins_xyz || !directions_xyz) {
        alea_set_error_detail(ALEA_ERR_NULL_ARG,
                              "coverage ray origins and directions are required");
        return -1;
    }
    if (row_count > SIZE_MAX / 3) {
        alea_set_error_detail(ALEA_ERR_OVERFLOW, "coverage row count is too large");
        return -1;
    }
    for (size_t row = 0; row < row_count; row++) {
        const double* o = &origins_xyz[row * 3];
        const double* d = &directions_xyz[row * 3];
        if (!isfinite(o[0]) || !isfinite(o[1]) || !isfinite(o[2]) ||
            !isfinite(d[0]) || !isfinite(d[1]) || !isfinite(d[2]) ||
            (d[0] == 0.0 && d[1] == 0.0 && d[2] == 0.0)) {
            alea_set_error_detail(ALEA_ERR_INVALID_ARG,
                                  "coverage row %zu has non-finite data or zero direction",
                                  row);
            return -1;
        }
    }
    return 0;
}

void alea_ray_coverage_slice_options_init(
    alea_ray_coverage_slice_options_t* options) {
    if (options) memset(options, 0, sizeof(*options));
    if (options) options->struct_size = sizeof(*options);
}

alea_ray_coverage_slice_result_t* alea_ray_coverage_slice_result_create(void) {
    alea_ray_coverage_slice_result_t* result = calloc(1, sizeof(*result));
    if (result) alea_ray_coverage_slice_result_init(result);
    return result;
}

void alea_ray_coverage_slice_result_destroy(
    alea_ray_coverage_slice_result_t* result) {
    if (!result) return;
    alea_ray_coverage_slice_result_free(result);
    free(result);
}

static int coverage_public_limit(size_t* out, uint64_t value,
                                 const char* name) {
    if (value > SIZE_MAX) {
        alea_set_error_detail(ALEA_ERR_OVERFLOW,
                              "coverage %s limit overflows size_t", name);
        return -1;
    }
    *out = (size_t)value;
    return 0;
}

int alea_ray_coverage_query(
    alea_system_t* sys, double ox, double oy, double oz,
    double dx, double dy, double dz,
    const alea_ray_coverage_slice_options_t* options,
    alea_ray_coverage_slice_result_t* result) {
    const double origin[3] = {ox, oy, oz};
    const double direction[3] = {dx, dy, dz};
    return alea_ray_coverage_slice_query(sys, origin, direction, 1, NULL,
                                         NULL, options, result);
}

int alea_ray_coverage_slice_query(
    alea_system_t* sys, const double* origins_xyz, const double* directions_xyz,
    size_t row_count, const uint8_t* direction_tags,
    const double* transverse_coordinates,
    const alea_ray_coverage_slice_options_t* options,
    alea_ray_coverage_slice_result_t* result) {
    const uint32_t known_flags = ALEA_RAY_COVERAGE_DOMAIN |
                                 ALEA_RAY_COVERAGE_REPORT_EXTERIOR;
    alea_ray_coverage_slice_limits_t limits = {0};
    alea_ray_coverage_row_t* rows = NULL;
    alea_ray_coverage_executor_t executor;
    size_t worker_count = 1;
    int rc = -1;
    alea_ray_coverage_executor_init(&executor);
    if (!sys || !result || !options) {
        alea_set_error_detail(ALEA_ERR_NULL_ARG,
                              "coverage query requires system, options, and result");
        return -1;
    }
    if (options->struct_size < sizeof(*options) || options->flags & ~known_flags ||
        !isfinite(options->t_max) || options->t_max <= 0.0 ||
        ((options->flags & ALEA_RAY_COVERAGE_REPORT_EXTERIOR) &&
         !(options->flags & ALEA_RAY_COVERAGE_DOMAIN)) ||
        ((options->flags & ALEA_RAY_COVERAGE_DOMAIN) &&
         (!isfinite(options->domain_t_min) || !isfinite(options->domain_t_max) ||
          options->domain_t_min < 0.0 ||
          options->domain_t_min > options->domain_t_max))) {
        alea_set_error_detail(ALEA_ERR_INVALID_ARG, "invalid coverage query options");
        return -1;
    }
    if (coverage_validate_ray_inputs(origins_xyz, directions_xyz, row_count) != 0)
        return -1;
    if (coverage_public_limit(&limits.max_rows, options->max_rows, "row") ||
        coverage_public_limit(&limits.max_intervals, options->max_intervals,
                              "interval") ||
        coverage_public_limit(&limits.max_owners, options->max_owners, "owner") ||
        coverage_public_limit(&limits.max_bytes, options->max_output_bytes,
                              "output-byte"))
        return -1;
    if (row_count > SIZE_MAX / sizeof(*rows)) {
        alea_set_error_detail(ALEA_ERR_OVERFLOW, "coverage row input overflows");
        return -1;
    }
    rows = calloc(row_count ? row_count : 1, sizeof(*rows));
    if (!rows) {
        alea_set_error_detail(ALEA_ERR_OUT_OF_MEMORY, "failed to allocate coverage rows");
        return -1;
    }
    for (size_t row = 0; row < row_count; row++) {
        const double* o = &origins_xyz[row * 3];
        const double* d = &directions_xyz[row * 3];
        if (alea_ray_init(&rows[row].ray, o[0], o[1], o[2], d[0], d[1], d[2]) != 0)
            goto cleanup;
        rows[row].t_max = options->t_max;
        rows[row].use_domain = (options->flags & ALEA_RAY_COVERAGE_DOMAIN) != 0;
        rows[row].domain = (alea_ray_coverage_domain_t){
            .t_min = options->domain_t_min, .t_max = options->domain_t_max,
            .has_domain = rows[row].use_domain,
            .report_allowed_exterior =
                (options->flags & ALEA_RAY_COVERAGE_REPORT_EXTERIOR) != 0
        };
        rows[row].direction_tag = direction_tags ? direction_tags[row] : 0;
        rows[row].transverse_coordinate =
            transverse_coordinates ? transverse_coordinates[row] : 0.0;
    }
    if (alea_system_prepare_query_caches(
            sys, ALEA_CACHE_HIER_SPATIAL | ALEA_CACHE_CELL_SURFACES) != 0)
        goto cleanup;
#ifdef _OPENMP
    worker_count = (size_t)omp_get_max_threads();
#endif
    if (worker_count == 0) worker_count = 1;
    if (row_count != 0 && worker_count > row_count) worker_count = row_count;
    if (alea_ray_coverage_executor_prepare(&executor, worker_count) != 0) goto cleanup;
    if (options->max_refinement_depth != 0) {
        alea_ray_coverage_refinement_policy_t policy;
        alea_ray_coverage_refinement_policy_init(&policy);
        if (options->refinement_signals != 0)
            policy.signals = options->refinement_signals;
        policy.min_transverse_spacing = options->min_transverse_spacing;
        policy.endpoint_displacement = options->endpoint_displacement;
        policy.crossing_density = options->crossing_density;
        rc = alea_ray_coverage_slice_build_adaptive_policy_executor_nocache(
            sys, rows, row_count, options->max_refinement_depth, &policy,
            &limits, &executor, result);
    } else {
        rc = alea_ray_coverage_slice_build_executor_nocache(
            sys, rows, row_count, &limits, &executor, result);
    }
cleanup:
    alea_ray_coverage_executor_free(&executor);
    free(rows);
    return rc;
}

size_t alea_ray_coverage_slice_row_count(const alea_ray_coverage_slice_result_t* r) { return r ? r->row_count : 0; }
size_t alea_ray_coverage_slice_interval_count(const alea_ray_coverage_slice_result_t* r) { return r ? r->interval_count : 0; }
size_t alea_ray_coverage_slice_owner_count(const alea_ray_coverage_slice_result_t* r) { return r ? r->owner_count : 0; }
int alea_ray_coverage_slice_refinement_status(const alea_ray_coverage_slice_result_t* r) { return r ? (int)r->refinement_status : -1; }
const size_t* alea_ray_coverage_slice_row_offsets(const alea_ray_coverage_slice_result_t* r) { return r ? r->row_offsets : NULL; }
const uint8_t* alea_ray_coverage_slice_row_direction_tags(const alea_ray_coverage_slice_result_t* r) { return r ? r->row_direction_tags : NULL; }
const double* alea_ray_coverage_slice_row_transverse_coordinates(const alea_ray_coverage_slice_result_t* r) { return r ? r->row_transverse_coordinates : NULL; }
const double* alea_ray_coverage_slice_t_enter(const alea_ray_coverage_slice_result_t* r) { return r ? r->t_enter : NULL; }
const double* alea_ray_coverage_slice_t_exit(const alea_ray_coverage_slice_result_t* r) { return r ? r->t_exit : NULL; }
const uint8_t* alea_ray_coverage_slice_kinds(const alea_ray_coverage_slice_result_t* r) { return r ? r->kinds : NULL; }
const size_t* alea_ray_coverage_slice_owner_offsets(const alea_ray_coverage_slice_result_t* r) { return r ? r->owner_offsets : NULL; }
const size_t* alea_ray_coverage_slice_owner_count_lower_bounds(const alea_ray_coverage_slice_result_t* r) { return r ? r->owner_count_lower_bounds : NULL; }
const int* alea_ray_coverage_slice_owner_cell_ids(const alea_ray_coverage_slice_result_t* r) { return r ? r->owner_cell_ids : NULL; }
const int* alea_ray_coverage_slice_owner_material_ids(const alea_ray_coverage_slice_result_t* r) { return r ? r->owner_material_ids : NULL; }
const int* alea_ray_coverage_slice_owner_universe_ids(const alea_ray_coverage_slice_result_t* r) { return r ? r->owner_universe_ids : NULL; }
const int* alea_ray_coverage_slice_owner_fill_universes(const alea_ray_coverage_slice_result_t* r) { return r ? r->owner_fill_universes : NULL; }
const int* alea_ray_coverage_slice_owner_depths(const alea_ray_coverage_slice_result_t* r) { return r ? r->owner_depths : NULL; }
const uint64_t* alea_ray_coverage_slice_owner_occurrence_keys(const alea_ray_coverage_slice_result_t* r) { return r ? r->owner_occurrence_keys : NULL; }
const uint64_t* alea_ray_coverage_slice_owner_parent_occurrence_keys(const alea_ray_coverage_slice_result_t* r) { return r ? r->owner_parent_occurrence_keys : NULL; }
const uint8_t* alea_ray_coverage_slice_owner_resolution_flags(const alea_ray_coverage_slice_result_t* r) { return r ? r->owner_resolution_flags : NULL; }
