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
    ALEA_GEOM_ERR_AMBIGUOUS_BOUNDARY,
    ALEA_GEOM_ERR_INTERIOR_GAP
} alea_geom_error_type_t;

typedef enum {
    ALEA_GEOM_VALIDATE_RAYS                = 1u << 0,
    ALEA_GEOM_VALIDATE_STRICT_ADJACENCY    = 1u << 1,
    ALEA_GEOM_VALIDATE_ALLOW_EXTERIOR_VOID = 1u << 2,
    /* validation_bounds defines a closed world-space AABB inside which
     * unowned coverage is an interior-gap finding. */
    ALEA_GEOM_VALIDATE_DOMAIN_BOUNDS       = 1u << 4
} alea_geom_validate_flags_t;

typedef enum {
    ALEA_GEOM_EVENT_INITIAL_POINT          = 1u << 0,
    ALEA_GEOM_EVENT_PREVIOUS_NO_SURFACE    = 1u << 1,
    ALEA_GEOM_EVENT_MISSING_ADJACENCY      = 1u << 2,
    ALEA_GEOM_EVENT_EXTERIOR_ALLOWED       = 1u << 3,
    ALEA_GEOM_EVENT_COINCIDENT_SURFACES    = 1u << 4,
    ALEA_GEOM_EVENT_TRUNCATED_COVERAGE     = 1u << 5,
    ALEA_GEOM_EVENT_FOUND_WITHOUT_ADJACENCY = 1u << 6,
    /* A slice sample lies on the observation-window boundary.  It has no
     * inward/outward support, so transition evidence is inconclusive. */
    ALEA_GEOM_EVENT_VIEWPORT_EDGE          = 1u << 7
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
    /* Retain at most this many raw samples with the same causal signature.
     * Zero disables this local cap. */
    size_t max_samples_per_signature;
    size_t max_crossings;
    double sample_offset;
    double t_max;
    uint64_t seed;
    int ray_count;
    /* [min_x, max_x, min_y, max_y, min_z, max_z]; consulted only with
     * ALEA_GEOM_VALIDATE_DOMAIN_BOUNDS. */
    double validation_bounds[6];
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
    size_t suppressed_samples;
    /* Private validator bookkeeping for bounded signature sampling. */
    void* signature_table;
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
/** Classify each row with the complete-coverage sweep.  This is the only
 *  check that can report gaps and overlaps: a bidirectional comparison can
 *  agree on the same incomplete ownership interpretation in both directions. */
#define ALEA_RAY_SLICE_VALIDATE_COVERAGE (1u << 1)

/** Trace-consistency evidence.  Not a geometry verdict. */
#define ALEA_RAY_SLICE_DIAG_FAST_DIRECTION_MISMATCH (1u << 0)

/** Confirmed complete-coverage classifications. */
#define ALEA_RAY_SLICE_DIAG_COVERAGE_GAP             (1u << 1)
#define ALEA_RAY_SLICE_DIAG_COVERAGE_OVERLAP         (1u << 2)
#define ALEA_RAY_SLICE_DIAG_COVERAGE_UNDEFINED_FILL  (1u << 3)
/** Coverage could not be established; neither a finding nor a clean result. */
#define ALEA_RAY_SLICE_DIAG_COVERAGE_UNRESOLVED      (1u << 4)
#define ALEA_RAY_SLICE_DIAG_COVERAGE_TRUNCATED       (1u << 5)
/** Unowned space outside the configured validation domain. */
#define ALEA_RAY_SLICE_DIAG_COVERAGE_ALLOWED_EXTERIOR (1u << 6)

#define ALEA_RAY_SLICE_VALIDATION_INCLUDE_AGREEMENTS (1u << 0)

/** Coverage domain policy.  Exactly one of these is required when
 *  ALEA_RAY_SLICE_VALIDATE_COVERAGE is requested: the viewport extent never
 *  defines by itself where ownership is required, so unowned space cannot be
 *  classified without an explicit domain or an explicit uniform policy. */
#define ALEA_RAY_SLICE_COVERAGE_HAS_DOMAIN          (1u << 0)
/** Treat all observed unowned space as allowed exterior; never a gap. */
#define ALEA_RAY_SLICE_COVERAGE_UNOWNED_IS_EXTERIOR (1u << 1)
/** Emit allowed-exterior intervals rather than omitting them. */
#define ALEA_RAY_SLICE_COVERAGE_REPORT_EXTERIOR     (1u << 2)

#define ALEA_RAY_SLICE_TRACE_FAST_FORWARD (1u << 0)
#define ALEA_RAY_SLICE_TRACE_FAST_REVERSE (1u << 1)

#define ALEA_RAY_SLICE_VALIDATION_FIELD_FAST_FORWARD (1u << 0)
#define ALEA_RAY_SLICE_VALIDATION_FIELD_FAST_REVERSE (1u << 1)
#define ALEA_RAY_SLICE_VALIDATION_FIELD_OCCURRENCE_KEYS (1u << 2)
#define ALEA_RAY_SLICE_VALIDATION_FIELD_COVERAGE (1u << 3)
/** Published rows carry transverse coordinates, direction tags, and base-row
 *  provenance.  Set whenever adaptive refinement was requested. */
#define ALEA_RAY_SLICE_VALIDATION_FIELD_ADAPTIVE_ROWS (1u << 4)

/** Adaptive refinement probe signals.  Signature difference alone cannot see a
 *  boundary that moves without changing owner identity, a region already
 *  carrying findings, or a dense row whose neighbours happen to agree. */
#define ALEA_RAY_SLICE_REFINE_SIGNATURE    (1u << 0)
#define ALEA_RAY_SLICE_REFINE_DISPLACEMENT (1u << 1)
#define ALEA_RAY_SLICE_REFINE_DENSITY      (1u << 2)
#define ALEA_RAY_SLICE_REFINE_FINDING      (1u << 3)

/** Why adaptive refinement stopped.  Every value but CONVERGED means the
 *  adaptive criteria did not converge: completed rows are published, but the
 *  report must not claim the viewport was sampled to the requested criteria. */
typedef enum {
    ALEA_RAY_SLICE_REFINEMENT_NOT_REQUESTED = 0,
    ALEA_RAY_SLICE_REFINEMENT_CONVERGED,
    ALEA_RAY_SLICE_REFINEMENT_MAX_DEPTH,
    ALEA_RAY_SLICE_REFINEMENT_MAX_ROWS,
    ALEA_RAY_SLICE_REFINEMENT_MIN_SPACING
} alea_ray_slice_refinement_status_t;

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
    uint32_t coverage_flags;   /**< ALEA_RAY_SLICE_COVERAGE_* */
    /** Explicit validation domain in view U coordinates, used only with
     *  ALEA_RAY_SLICE_COVERAGE_HAS_DOMAIN.  Unowned space inside it is an
     *  interior gap; unowned space outside it is allowed exterior.  Viewport
     *  clipping narrows observation only and never reclassifies either. */
    double coverage_domain_u_min;
    double coverage_domain_u_max;
    /** Publication budget for retained concrete owner records across the whole
     *  slice.  Exhausting it is an operation failure that publishes nothing,
     *  which is distinct from per-interval owner-set saturation: that is a
     *  successful result carrying ALEA_RAY_SLICE_DIAG_COVERAGE_TRUNCATED, and
     *  it prevents a successful validation verdict.  0 = unlimited. */
    uint64_t max_coverage_owners;

    /** Maximum adaptive refinement waves.  0 (the default) samples exactly the
     *  requested rows and leaves the published row count equal to row_count.
     *  Above 0, refined rows are published between the requested rows in
     *  transverse order, so the published row count grows and the caller must
     *  read row provenance rather than assume a one-to-one row mapping.
     *  Refinement drives coverage probes only: a refined row carries no
     *  bidirectional or event-cache evidence, because both are defined on the
     *  requested viewport rows. */
    uint32_t coverage_max_refinement_depth;
    /** ALEA_RAY_SLICE_REFINE_* signals.  0 selects signature difference. */
    uint32_t coverage_refine_signals;
    /** Required, and only used, with ALEA_RAY_SLICE_REFINE_DISPLACEMENT. */
    double coverage_endpoint_displacement;
    /** Required, and only used, with ALEA_RAY_SLICE_REFINE_DENSITY. */
    uint64_t coverage_crossing_density;
    /** Transverse spacing below which a row pair is never split.  0 = no
     *  spacing limit. */
    double coverage_min_transverse_spacing;
    /** Publication budget for total sampled rows.  Reaching it publishes the
     *  completed rows with MAX_ROWS status rather than failing.  0 = unlimited. */
    uint64_t max_coverage_rows;
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

/* Adaptive row provenance.  Present only with
 * ALEA_RAY_SLICE_VALIDATION_FIELD_ADAPTIVE_ROWS.  Each array has row_count
 * entries and is ordered first by direction ordinal, then by transverse view
 * coordinate.  A base-row index of SIZE_MAX marks a refined row, which carries
 * coverage evidence only. */
const double* alea_ray_slice_validation_row_transverse_coordinates(
    const alea_ray_slice_validation_result_t* result);
const uint8_t* alea_ray_slice_validation_row_direction_tags(
    const alea_ray_slice_validation_result_t* result);
const size_t* alea_ray_slice_validation_row_base_indices(
    const alea_ray_slice_validation_result_t* result);

/* Why refinement stopped.  Anything but CONVERGED (or NOT_REQUESTED) means the
 * published rows are an explicitly refinement-limited sample. */
alea_ray_slice_refinement_status_t alea_ray_slice_validation_refinement_status(
    const alea_ray_slice_validation_result_t* result);

/* Retained concrete owner records backing each interval's coverage
 * classification.  Present only when ALEA_RAY_SLICE_VALIDATION_FIELD_COVERAGE
 * is set.  For a TRUNCATED interval this is the retained prefix, not the
 * complete owner count. */
const uint32_t* alea_ray_slice_validation_coverage_owner_counts(
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
