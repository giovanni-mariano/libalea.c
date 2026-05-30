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

#ifdef __cplusplus
}
#endif

#endif /* ALEA_GEO_VALIDATOR_H */
