// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/**
 * @file alea_geo_validator.h
 * @brief Transport-style geometry validation API
 */

#ifndef ALEA_GEO_VALIDATOR_H
#define ALEA_GEO_VALIDATOR_H

#include "alea.h"
#include "alea_raycast.h"
#include "alea_slice.h"
#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef enum {
    ALEA_GEOM_ERR_UNDEFINED_AFTER_CROSSING = 1,
    ALEA_GEOM_ERR_OVERLAP_AFTER_CROSSING,
    ALEA_GEOM_ERR_NON_ADJACENT_TRANSITION,
    ALEA_GEOM_ERR_MISSING_NEIGHBOR,
    ALEA_GEOM_ERR_AMBIGUOUS_BOUNDARY
} alea_geom_error_type_t;

typedef enum {
    ALEA_GEOM_VALIDATE_RAYS                = 1u << 0,
    ALEA_GEOM_VALIDATE_STRICT_ADJACENCY    = 1u << 1,
    ALEA_GEOM_VALIDATE_ALLOW_EXTERIOR_VOID = 1u << 2,
    ALEA_GEOM_VALIDATE_HIERARCHICAL        = 1u << 3
} alea_geom_validate_flags_t;

typedef enum {
    ALEA_GEOM_EVENT_INITIAL_POINT          = 1u << 0,
    ALEA_GEOM_EVENT_PREVIOUS_NO_SURFACE    = 1u << 1,
    ALEA_GEOM_EVENT_MISSING_ADJACENCY      = 1u << 2,
    ALEA_GEOM_EVENT_EXTERIOR_ALLOWED       = 1u << 3,
    ALEA_GEOM_EVENT_COINCIDENT_SURFACES    = 1u << 4,
    ALEA_GEOM_EVENT_TRUNCATED_COVERAGE     = 1u << 5,
    ALEA_GEOM_EVENT_FOUND_WITHOUT_ADJACENCY = 1u << 6
} alea_geom_event_flags_t;

typedef enum {
    ALEA_GEOM_EVENT_SOURCE_UNKNOWN = 0,
    ALEA_GEOM_EVENT_SOURCE_RAY,
    ALEA_GEOM_EVENT_SOURCE_SLICE_CURVE,
    ALEA_GEOM_EVENT_SOURCE_INITIAL_POINT
} alea_geom_event_source_t;

typedef struct {
    unsigned flags;
    int universe_depth;
    size_t max_errors;
    size_t max_crossings;
    double sample_offset;
    double t_max;
    uint64_t seed;
    int ray_count;
} alea_geom_validator_options_t;

typedef struct {
    alea_geom_error_type_t type;
    alea_geom_event_source_t source;
    int previous_cell_id;
    int found_cell_id;
    int expected_neighbor_cell_id;
    int secondary_cell_id;
    int found_cell_count;
    int surface_id;
    uint32_t primitive_id;
    int universe_id;
    int universe_depth;
    double crossing_point[3];
    double sample_point[3];
    double direction[3];
    double t;                       /**< Ray distance for ray events; curve parameter for slice events */
    double offset;
    size_t curve_index;             /**< Index into caller-provided curves for slice events, SIZE_MAX otherwise */
    uint32_t component_index;        /**< Component/edge/branch index for slice events */
    double uv[2];                    /**< Slice-plane coordinate for slice events */
    uint32_t flags;
} alea_geom_error_t;

typedef struct {
    alea_geom_error_t* errors;
    size_t error_count;
    size_t error_capacity;
    size_t crossings_checked;
    size_t adjacency_hits;
    size_t exact_queries;
    size_t ambiguous_crossings;
    int truncated;
} alea_geom_validator_result_t;

void alea_geom_validator_options_init(alea_geom_validator_options_t* options);
void alea_geom_validator_result_init(alea_geom_validator_result_t* result);
void alea_geom_validator_result_free(alea_geom_validator_result_t* result);

const char* alea_geom_error_type_name(alea_geom_error_type_t type);
size_t alea_geom_validator_error_count(const alea_geom_validator_result_t* result);
int alea_geom_validator_error_get(const alea_geom_validator_result_t* result,
                                  size_t index,
                                  alea_geom_error_t* out_error);

int alea_validate_geometry(alea_system_t* sys,
                           const alea_geom_validator_options_t* options,
                           alea_geom_validator_result_t* result);

int alea_validate_geometry_ray(alea_system_t* sys,
                               const alea_geom_validator_options_t* options,
                               double ox, double oy, double oz,
                               double dx, double dy, double dz,
                               double t_max,
                               alea_geom_validator_result_t* result);

/* Surface/slice-driven validation: samples caller-provided analytical boundary
 * curves on `view`'s plane and emits the same structured events as the
 * ray-driven path, augmented with curve_index/t/uv provenance. Catches hidden
 * nested overlaps that random rays can miss. */
int alea_validate_geometry_slice(alea_system_t* sys,
                                 const alea_slice_view_t* view,
                                 const alea_slice_curves_t* curves,
                                 const alea_geom_validator_options_t* options,
                                 alea_geom_validator_result_t* result);

/* ==========================================================================
 * COMPACT RAY-SLICE DIRECTIONAL VALIDATION
 * ==========================================================================
 *
 * This is a viewport-local diagnostic. A forward/reverse disagreement is a
 * ray-tracing consistency signal; it is not by itself a geometry error.
 */

typedef struct alea_ray_slice_validation_result
    alea_ray_slice_validation_result_t;

#define ALEA_RAY_SLICE_VALIDATE_FAST_BIDIRECTIONAL (1u << 0)

#define ALEA_RAY_SLICE_DIAG_FAST_DIRECTION_MISMATCH (1u << 0)

#define ALEA_RAY_SLICE_VALIDATION_INCLUDE_AGREEMENTS (1u << 0)

#define ALEA_RAY_SLICE_TRACE_FAST_FORWARD (1u << 0)
#define ALEA_RAY_SLICE_TRACE_FAST_REVERSE (1u << 1)

#define ALEA_RAY_SLICE_VALIDATION_FIELD_FAST_FORWARD (1u << 0)
#define ALEA_RAY_SLICE_VALIDATION_FIELD_FAST_REVERSE (1u << 1)
#define ALEA_RAY_SLICE_VALIDATION_FIELD_OCCURRENCE_KEYS (1u << 2)

typedef struct {
    size_t struct_size;
    uint32_t checks;
    uint32_t flags;
    int projected_depth;       /**< -1 = resolved leaf */
    double absolute_tolerance; /**< 0 = native default */
    double relative_tolerance; /**< 0 = native default */
    uint64_t max_trace_intervals;
    uint64_t max_path_entries;
    uint64_t max_output_intervals;
    uint64_t max_output_bytes;
} alea_ray_slice_validation_options_t;

void alea_ray_slice_validation_options_init(
    alea_ray_slice_validation_options_t* options);

alea_ray_slice_validation_result_t*
alea_ray_slice_validation_result_create(void);

void alea_ray_slice_validation_result_destroy(
    alea_ray_slice_validation_result_t* result);

size_t alea_ray_slice_validation_row_count(
    const alea_ray_slice_validation_result_t* result);
size_t alea_ray_slice_validation_interval_count(
    const alea_ray_slice_validation_result_t* result);
uint32_t alea_ray_slice_validation_fields(
    const alea_ray_slice_validation_result_t* result);
uint32_t alea_ray_slice_validation_executed_trace_mask(
    const alea_ray_slice_validation_result_t* result);
uint32_t alea_ray_slice_validation_reused_trace_mask(
    const alea_ray_slice_validation_result_t* result);

const uint64_t* alea_ray_slice_validation_row_offsets(
    const alea_ray_slice_validation_result_t* result);
const double* alea_ray_slice_validation_u_enter(
    const alea_ray_slice_validation_result_t* result);
const double* alea_ray_slice_validation_u_exit(
    const alea_ray_slice_validation_result_t* result);
const uint32_t* alea_ray_slice_validation_diagnostic_flags(
    const alea_ray_slice_validation_result_t* result);
const int32_t* alea_ray_slice_validation_fast_forward_cell_ids(
    const alea_ray_slice_validation_result_t* result);
const int32_t* alea_ray_slice_validation_fast_reverse_cell_ids(
    const alea_ray_slice_validation_result_t* result);
const uint64_t* alea_ray_slice_validation_fast_forward_occurrence_keys(
    const alea_ray_slice_validation_result_t* result);
const uint64_t* alea_ray_slice_validation_fast_reverse_occurrence_keys(
    const alea_ray_slice_validation_result_t* result);

/* Optional provenance attached by validation with a directional trace cache.
 * Each array has interval_count entries; a surface ID of -1 means no physical
 * surface is reportable at that endpoint. */
const int32_t* alea_ray_slice_validation_u_enter_forward_surface_ids(
    const alea_ray_slice_validation_result_t* result);
const int32_t* alea_ray_slice_validation_u_enter_reverse_surface_ids(
    const alea_ray_slice_validation_result_t* result);
const int32_t* alea_ray_slice_validation_u_exit_forward_surface_ids(
    const alea_ray_slice_validation_result_t* result);
const int32_t* alea_ray_slice_validation_u_exit_reverse_surface_ids(
    const alea_ray_slice_validation_result_t* result);
const uint32_t* alea_ray_slice_validation_u_enter_provenance_flags(
    const alea_ray_slice_validation_result_t* result);
const uint32_t* alea_ray_slice_validation_u_exit_provenance_flags(
    const alea_ray_slice_validation_result_t* result);

#define ALEA_RAY_SLICE_BOUNDARY_PROVENANCE_COINCIDENT (1u << 0)
#define ALEA_RAY_SLICE_BOUNDARY_PROVENANCE_SYNTHETIC  (1u << 1)
#define ALEA_RAY_SLICE_BOUNDARY_PROVENANCE_UNRESOLVED (1u << 2)

/**
 * Validate centered U-directed rows. `inout_fast_forward` is optional. When
 * it holds a compatible result from alea_trace_ray_slice_compact(), the
 * forward trace is reused and only the reverse trace is executed. Otherwise a
 * replacement forward result is published there only after full success.
 */
int alea_validate_ray_slice_compact(
    alea_system_t* sys,
    const alea_slice_view_t* view,
    size_t row_count,
    const alea_ray_slice_validation_options_t* validation_options,
    const alea_raycast_batch_options_t* render_options,
    alea_raycast_batch_result_t* inout_fast_forward,
    alea_ray_slice_validation_result_t* out_validation);

/** As above, but reuses a matching public directional trace cache and emits
 * canonical boundary provenance.  A mismatch is rejected transactionally. */
int alea_validate_ray_slice_compact_with_directional_cache(
    alea_system_t* sys, const alea_slice_view_t* view, size_t row_count,
    const alea_ray_slice_validation_options_t* validation_options,
    const alea_raycast_batch_options_t* render_options,
    alea_raycast_batch_result_t* inout_fast_forward,
    const alea_slice_directional_trace_cache_t* cache,
    alea_ray_slice_validation_result_t* out_validation);

#ifdef __cplusplus
}
#endif

#endif /* ALEA_GEO_VALIDATOR_H */
