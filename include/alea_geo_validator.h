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
    /* Bound analytical sampling work per curve. Zero uses the hard safety
     * ceiling (2000); interactive callers use the default 512. */
    size_t max_samples_per_curve;
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
    size_t sample_limited_curves;
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
 * SUPPLIED UNIVERSE-LOCAL TRANSITION DIAGNOSTIC
 * ========================================================================== */

typedef enum {
    ALEA_TRANSITION_VALID = 0,
    ALEA_TRANSITION_GAP,
    ALEA_TRANSITION_OVERLAP,
    ALEA_TRANSITION_UNDEFINED_FILL,
    ALEA_TRANSITION_MISSING_NEIGHBOR,
    ALEA_TRANSITION_NON_ADJACENT,
    ALEA_TRANSITION_AMBIGUOUS_NEIGHBOR,
    ALEA_TRANSITION_SURFACE_CHAIN_CORNER,
    ALEA_TRANSITION_AMBIGUOUS_BOUNDARY,
    ALEA_TRANSITION_UNRESOLVED,
    ALEA_TRANSITION_TRUNCATED
} alea_transition_kind_t;

#define ALEA_TRANSITION_EVIDENCE_CAPACITY 16

#define ALEA_TRANSITION_FLAG_CURRENT_BEFORE_CONTAINS (1u << 0)
#define ALEA_TRANSITION_FLAG_CURRENT_AFTER_CONTAINS  (1u << 1)
#define ALEA_TRANSITION_FLAG_COVERAGE_FALLBACK       (1u << 2)
#define ALEA_TRANSITION_FLAG_PRIMARY_MISSING         (1u << 3)
#define ALEA_TRANSITION_FLAG_TIED_SURFACE_CONNECTS   (1u << 4)
#define ALEA_TRANSITION_FLAG_OFFSET_STABLE           (1u << 5)
#define ALEA_TRANSITION_FLAG_CANDIDATES_TRUNCATED    (1u << 6)
#define ALEA_TRANSITION_FLAG_OWNERS_TRUNCATED        (1u << 7)

typedef struct {
    size_t struct_size;
    /** Central offset in model units. Zero uses the validator default. */
    double probe_distance;
    /** Upper bound for the offset ladder. Zero uses 8*probe_distance. */
    double max_probe_distance;
    /** Complete-coverage hit budget, including ancestry. Zero uses 256. */
    size_t max_coverage_hits;
    /** Maximum complete-coverage calls in this transition. Zero is unlimited. */
    size_t max_coverage_fallbacks;
} alea_transition_options_t;

typedef struct {
    alea_transition_kind_t kind;
    alea_point_coverage_kind_t after_coverage_kind;
    uint32_t flags;
    int universe_id;
    int current_cell_id;
    int primary_surface_id;
    int connecting_surface_id;
    int after_cell_id;
    int current_sense;
    int occurrence_depth;
    uint64_t current_occurrence_key;
    uint64_t current_parent_occurrence_key;
    uint64_t before_occurrence_key;
    uint64_t before_parent_occurrence_key;
    uint64_t selected_after_occurrence_key;
    uint64_t selected_after_parent_occurrence_key;
    double crossing_point[3];
    double direction[3];
    double before_point[3];
    double after_point[3];
    double probe_distance;
    size_t offset_attempts;
    size_t coverage_fallbacks;
    size_t primary_candidate_count;
    size_t primary_containing_count;
    size_t candidate_cell_count;
    int candidate_cell_ids[ALEA_TRANSITION_EVIDENCE_CAPACITY];
    size_t after_owner_count;
    size_t owner_cell_count;
    int owner_cell_ids[ALEA_TRANSITION_EVIDENCE_CAPACITY];
} alea_transition_result_t;

void alea_transition_options_init(alea_transition_options_t* options);
const char* alea_transition_kind_name(alea_transition_kind_t kind);

/** Check a supplied crossing witness in an active universe's local frame.
 * ``tied_surface_ids`` excludes or may repeat the primary; duplicates are
 * ignored. Query acceleration must already be prepared.
 */
int alea_check_transition_local(
    alea_system_t* sys,
    int universe_id,
    int current_cell_id,
    int primary_surface_id,
    const int* tied_surface_ids,
    size_t tied_surface_count,
    const double point[3],
    const double direction[3],
    const alea_transition_options_t* options,
    alea_transition_result_t* result);

/* ==========================================================================
 * MEMORY-BOUNDED 2D TRANSITION SCREEN
 * ========================================================================== */

typedef enum {
    ALEA_TRANSITION_SLICE_HORIZONTAL = 0,
    ALEA_TRANSITION_SLICE_VERTICAL = 1
} alea_transition_slice_orientation_t;

typedef enum {
    ALEA_TRANSITION_SLICE_STOP_NONE = 0,
    ALEA_TRANSITION_SLICE_STOP_MAX_RAYS,
    ALEA_TRANSITION_SLICE_STOP_MAX_EVENTS,
    ALEA_TRANSITION_SLICE_STOP_MAX_FINDINGS,
    ALEA_TRANSITION_SLICE_STOP_MAX_COMPONENTS,
    ALEA_TRANSITION_SLICE_STOP_MAX_OUTPUT_BYTES,
    ALEA_TRANSITION_SLICE_STOP_MAX_COVERAGE_FALLBACKS,
    ALEA_TRANSITION_SLICE_STOP_MAX_COVERAGE_PROBES,
    ALEA_TRANSITION_SLICE_STOP_MAX_SCRATCH_BYTES,
    ALEA_TRANSITION_SLICE_STOP_INTERRUPTED,
    ALEA_TRANSITION_SLICE_STOP_MAX_COMPONENT_LINKS
} alea_transition_slice_stop_reason_t;

#define ALEA_TRANSITION_SLICE_REFINE_SIGNATURE (1u << 0)
#define ALEA_TRANSITION_SLICE_REFINE_FINDING   (1u << 1)

typedef enum {
    ALEA_TRANSITION_SLICE_REFINEMENT_NOT_REQUESTED = 0,
    ALEA_TRANSITION_SLICE_REFINEMENT_CONVERGED,
    ALEA_TRANSITION_SLICE_REFINEMENT_MAX_DEPTH,
    ALEA_TRANSITION_SLICE_REFINEMENT_MIN_SPACING,
    ALEA_TRANSITION_SLICE_REFINEMENT_STOPPED
} alea_transition_slice_refinement_status_t;

typedef struct {
    size_t struct_size;
    size_t horizontal_rays;
    size_t vertical_rays;
    uint64_t max_rays;
    uint64_t max_events;
    uint64_t max_events_per_ray;
    uint64_t max_findings;
    uint64_t max_components;
    uint64_t max_output_bytes;
    uint64_t max_scratch_bytes;
    uint64_t max_coverage_fallbacks;
    size_t max_coverage_hits;
    double probe_distance;
    double max_probe_distance;
    int include_void_transitions;
    uint32_t max_refinement_depth;
    uint32_t refine_signals;
    double min_transverse_spacing;
    uint64_t max_row_scratch_bytes;
    size_t coverage_uniform_probes_per_ray;
    int coverage_probe_selected_intervals;
    int report_unowned_coverage;
    uint64_t max_coverage_probes;
    uint64_t max_coverage_findings;
    uint64_t max_coverage_components;
    uint64_t max_component_links;
    int enable_critical_refinement;
    int critical_full_view;
    size_t max_refinement_frontiers;
    size_t max_critical_tiles;
    size_t max_critical_tile_sources;
    uint64_t max_critical_scratch_bytes;
    double critical_tile_padding;
    size_t max_curves_per_tile;
    size_t max_critical_points;
    uint64_t max_active_boundary_tests;
    uint64_t max_critical_probes;
    uint64_t max_critical_findings;
    uint64_t max_curve_pairs;
    uint64_t max_critical_sector_witnesses;
    double critical_probe_radius;
    uint64_t max_critical_boundary_evidence;
} alea_transition_slice_options_t;

#define ALEA_TRANSITION_SLICE_BOUNDARY_PIECE_CAPACITY 3
#define ALEA_TRANSITION_SLICE_BOUNDARY_POINT_CAPACITY 17

#define ALEA_TRANSITION_SLICE_BOUNDARY_ROLE_SOURCE     (1u << 0)
#define ALEA_TRANSITION_SLICE_BOUNDARY_ROLE_PRIMARY    (1u << 1)
#define ALEA_TRANSITION_SLICE_BOUNDARY_ROLE_CONNECTING (1u << 2)

typedef struct {
    int surface_id;
    uint32_t role_flags;
    size_t point_count;
    double uv[ALEA_TRANSITION_SLICE_BOUNDARY_POINT_CAPACITY][2];
} alea_transition_slice_boundary_piece_t;

typedef struct {
    alea_transition_result_t transition;
    alea_transition_slice_orientation_t orientation;
    size_t ray_index;
    size_t event_index;
    size_t base_ray_index;
    uint32_t refinement_depth;
    double transverse_coordinate;
    double ray_t;
    double uv[2];
    double world_point[3];
} alea_transition_slice_finding_t;

typedef struct {
    alea_transition_kind_t kind;
    alea_transition_slice_orientation_t orientation;
    int universe_id;
    int current_cell_id;
    int after_cell_id;
    int primary_surface_id;
    int connecting_surface_id;
    uint64_t current_occurrence_key;
    size_t first_finding_index;
    size_t finding_count;
    uint32_t max_refinement_depth;
    double uv_min[2];
    double uv_max[2];
    double world_min[3];
    double world_max[3];
} alea_transition_slice_component_t;

#define ALEA_TRANSITION_SLICE_COVERAGE_OWNER_CAPACITY 16

typedef struct {
    alea_point_coverage_kind_t kind;
    int truncated;
    int target_depth;
    size_t owner_count;
    size_t owner_count_lower_bound;
    int owner_cell_ids[ALEA_TRANSITION_SLICE_COVERAGE_OWNER_CAPACITY];
    int owner_universe_ids[ALEA_TRANSITION_SLICE_COVERAGE_OWNER_CAPACITY];
    int owner_depths[ALEA_TRANSITION_SLICE_COVERAGE_OWNER_CAPACITY];
    uint64_t owner_occurrence_keys[
        ALEA_TRANSITION_SLICE_COVERAGE_OWNER_CAPACITY];
    uint64_t owner_parent_occurrence_keys[
        ALEA_TRANSITION_SLICE_COVERAGE_OWNER_CAPACITY];
    alea_transition_slice_orientation_t orientation;
    size_t ray_index;
    size_t base_ray_index;
    uint32_t refinement_depth;
    double transverse_coordinate;
    double ray_t;
    double bracket_t_enter;
    double bracket_t_exit;
    double uv[2];
    double world_point[3];
} alea_transition_slice_coverage_finding_t;

typedef struct {
    alea_point_coverage_kind_t kind;
    int truncated;
    alea_transition_slice_orientation_t orientation;
    size_t first_finding_index;
    size_t finding_count;
    size_t owner_count_lower_bound;
    uint32_t max_refinement_depth;
    double uv_min[2];
    double uv_max[2];
    double world_min[3];
    double world_max[3];
} alea_transition_slice_coverage_component_t;

#define ALEA_TRANSITION_SLICE_LINK_ENTER (1u << 0)
#define ALEA_TRANSITION_SLICE_LINK_EXIT  (1u << 1)

typedef struct {
    size_t transition_component_index;
    size_t coverage_component_index;
    uint32_t boundary_sides;
    size_t witness_pair_count;
} alea_transition_slice_component_link_t;

typedef enum {
    ALEA_TRANSITION_SLICE_CRITICAL_DISABLED = 0,
    ALEA_TRANSITION_SLICE_CRITICAL_NONE,
    ALEA_TRANSITION_SLICE_CRITICAL_MAX_FRONTIERS,
    ALEA_TRANSITION_SLICE_CRITICAL_MAX_TILE_SOURCES,
    ALEA_TRANSITION_SLICE_CRITICAL_MAX_TILES,
    ALEA_TRANSITION_SLICE_CRITICAL_MAX_SCRATCH_BYTES,
    ALEA_TRANSITION_SLICE_CRITICAL_MAX_OUTPUT_BYTES,
    ALEA_TRANSITION_SLICE_CRITICAL_MAX_CURVES,
    ALEA_TRANSITION_SLICE_CRITICAL_CHAIN_TRUNCATED,
    ALEA_TRANSITION_SLICE_CRITICAL_MAX_POINTS,
    ALEA_TRANSITION_SLICE_CRITICAL_MAX_PROBES,
    ALEA_TRANSITION_SLICE_CRITICAL_MAX_FINDINGS,
    ALEA_TRANSITION_SLICE_CRITICAL_UNSUPPORTED_CURVE,
    ALEA_TRANSITION_SLICE_CRITICAL_MAX_CURVE_PAIRS,
    ALEA_TRANSITION_SLICE_CRITICAL_MAX_SECTOR_WITNESSES
} alea_transition_slice_critical_stop_reason_t;

typedef enum {
    ALEA_TRANSITION_SLICE_TILE_SOURCE_TRANSITION_COMPONENT = 0,
    ALEA_TRANSITION_SLICE_TILE_SOURCE_COVERAGE_COMPONENT,
    ALEA_TRANSITION_SLICE_TILE_SOURCE_REFINEMENT_FRONTIER,
    ALEA_TRANSITION_SLICE_TILE_SOURCE_FULL_VIEW
} alea_transition_slice_tile_source_kind_t;

typedef struct {
    alea_transition_slice_orientation_t orientation;
    uint32_t refinement_depth;
    double transverse_min;
    double transverse_max;
    uint64_t signature_a[2];
    uint64_t signature_b[2];
    size_t max_event_count;
    double uv_min[2];
    double uv_max[2];
} alea_transition_slice_refinement_frontier_t;

typedef struct {
    double uv_min[2];
    double uv_max[2];
    uint32_t source_flags;
    size_t first_source_index;
    size_t source_count;
} alea_transition_slice_critical_tile_t;

typedef struct {
    alea_transition_slice_tile_source_kind_t kind;
    size_t source_index;
} alea_transition_slice_critical_tile_source_t;

typedef struct {
    alea_transition_result_t transition;
    size_t tile_index;
    size_t point_index;
    int source_cell_id;
    int source_surface_id;
    uint64_t source_occurrence_key;
    uint64_t source_universe_occurrence_key;
    double uv[2];
    double world_point[3];
    double direction[3];
    double radius;
    size_t boundary_piece_count;
    int boundary_evidence_truncated;
    alea_transition_slice_boundary_piece_t boundary_pieces[
        ALEA_TRANSITION_SLICE_BOUNDARY_PIECE_CAPACITY];
} alea_transition_slice_critical_finding_t;

typedef struct {
    size_t requested_rays;
    size_t executed_rays;
    size_t horizontal_rays_executed;
    size_t vertical_rays_executed;
    size_t events_checked;
    size_t physical_events_seen;
    size_t valid_transitions;
    size_t findings;
    size_t components;
    size_t coverage_probes;
    size_t unique_coverage_probes;
    size_t coverage_findings;
    size_t coverage_components;
    size_t truncated_coverage_probes;
    size_t skipped_unowned_coverage_probes;
    size_t coverage_fallbacks;
    size_t skipped_void_transitions;
    size_t peak_live_events;
    size_t peak_live_event_bytes;
    size_t retained_output_bytes;
    size_t refined_rays_executed;
    size_t peak_row_scratch_bytes;
    size_t peak_scratch_bytes;
    uint32_t max_refinement_depth_reached;
    int complete;
    int converged;
    alea_transition_slice_stop_reason_t stop_reason;
    alea_transition_slice_refinement_status_t refinement_status;
    size_t component_links;
    size_t refinement_frontiers;
    size_t omitted_refinement_frontiers;
    size_t critical_tile_seeds;
    size_t critical_tiles;
    size_t critical_tile_sources;
    size_t omitted_critical_tile_sources;
    size_t peak_critical_scratch_bytes;
    size_t critical_tiles_processed;
    size_t critical_tiles_saturated;
    size_t critical_region_hits;
    size_t critical_region_candidates_scanned;
    size_t critical_occurrence_seed_points;
    size_t critical_occurrence_paths;
    size_t critical_occurrence_universe_queries;
    size_t critical_root_region_fallbacks;
    size_t critical_chain_truncated_hits;
    size_t critical_surface_references;
    size_t critical_duplicate_surface_occurrences;
    size_t critical_curves;
    size_t critical_curves_culled;
    size_t critical_ranked_curves_omitted;
    size_t critical_active_boundary_tests;
    size_t critical_active_boundary_fallbacks;
    size_t critical_active_capacity_fallbacks;
    size_t critical_active_test_budget_fallbacks;
    size_t critical_active_unsupported_parabola_fallbacks;
    size_t critical_active_unsupported_hyperbola_fallbacks;
    size_t critical_active_unsupported_quartic_fallbacks;
    size_t critical_active_unsupported_polygon_fallbacks;
    size_t critical_active_unsupported_point_fallbacks;
    size_t critical_active_unsupported_other_fallbacks;
    size_t critical_active_evaluation_line_fallbacks;
    size_t critical_active_evaluation_closed_conic_fallbacks;
    size_t critical_active_evaluation_general_conic_fallbacks;
    size_t critical_active_evaluation_quartic_fallbacks;
    size_t critical_active_evaluation_polygon_fallbacks;
    size_t critical_active_evaluation_other_fallbacks;
    size_t critical_active_open_conic_canonical_fallbacks;
    size_t critical_active_open_conic_breakpoint_fallbacks;
    size_t critical_active_segments;
    size_t critical_whole_curve_fallbacks;
    size_t peak_critical_curves;
    size_t critical_point_candidates;
    size_t critical_points;
    size_t critical_duplicate_points;
    size_t critical_unsupported_curves;
    size_t critical_unsupported_parabola_curves;
    size_t critical_unsupported_hyperbola_curves;
    size_t critical_unsupported_quartic_curves;
    size_t critical_unsupported_other_curves;
    size_t critical_probes;
    size_t critical_probe_events;
    size_t critical_probe_findings;
    size_t critical_findings;
    size_t omitted_critical_findings;
    size_t critical_curve_pair_candidates;
    size_t critical_curve_pairs_tested;
    size_t critical_unsupported_curve_pairs;
    size_t critical_unsupported_general_conic_pairs;
    size_t critical_unsupported_quartic_pairs;
    size_t critical_unsupported_quartic_line_pairs;
    size_t critical_unsupported_quartic_closed_conic_pairs;
    size_t critical_unsupported_quartic_general_conic_pairs;
    size_t critical_unsupported_quartic_quartic_pairs;
    size_t critical_unsupported_quartic_other_pairs;
    size_t critical_unsupported_polygon_pairs;
    size_t critical_unsupported_other_pairs;
    size_t critical_pair_algebraic_points;
    size_t critical_pair_domain_rejections;
    /* Pair points newly inserted into the unique critical-point set.  A
     * valid pair point may already have been emitted by single-curve logic. */
    size_t critical_pair_intersection_points;
    size_t critical_sector_witnesses;
    size_t critical_sector_gap_witnesses;
    size_t critical_sector_overlap_witnesses;
    size_t critical_sector_unresolved_witnesses;
    int critical_enabled;
    int critical_complete;
    alea_transition_slice_critical_stop_reason_t critical_stop_reason;
    size_t critical_boundary_evidence;
    size_t omitted_critical_boundary_evidence;
} alea_transition_slice_stats_t;

typedef struct alea_transition_slice_result alea_transition_slice_result_t;

typedef struct {
    size_t page_count;
    size_t completed_page_count;
    size_t requested_workers;
    size_t actual_workers;
    uint64_t reserved_scratch_bytes_per_worker;
    uint64_t reserved_parallel_scratch_bytes;
} alea_transition_slice_batch_stats_t;

void alea_transition_slice_options_init(
    alea_transition_slice_options_t* options);
const char* alea_transition_slice_stop_reason_name(
    alea_transition_slice_stop_reason_t reason);
const char* alea_transition_slice_refinement_status_name(
    alea_transition_slice_refinement_status_t status);
const char* alea_transition_slice_critical_stop_reason_name(
    alea_transition_slice_critical_stop_reason_t reason);
alea_transition_slice_result_t* alea_transition_slice_result_create(void);
void alea_transition_slice_result_destroy(
    alea_transition_slice_result_t* result);
int alea_transition_slice_screen(
    alea_system_t* sys, const alea_slice_view_t* view,
    const alea_transition_slice_options_t* options,
    alea_transition_slice_result_t* result);
/**
 * Screen independent slice pages with bounded OpenMP workers.
 *
 * Results are written by page ordinal and are therefore deterministic with
 * respect to worker scheduling. Every entry in `results` must be a distinct,
 * initialized result object. A zero worker request uses the OpenMP runtime
 * maximum; a zero parallel scratch budget selects the serial path. The actual
 * worker count never exceeds page count or
 * `max_parallel_scratch_bytes / reserved_scratch_bytes_per_worker`. Nested
 * OpenMP execution is disabled by falling back to one worker.
 */
int alea_transition_slice_screen_batch(
    alea_system_t* sys, const alea_slice_view_t* views, size_t page_count,
    const alea_transition_slice_options_t* options,
    size_t requested_workers, uint64_t max_parallel_scratch_bytes,
    alea_transition_slice_result_t* const* results,
    alea_transition_slice_batch_stats_t* out_stats);
size_t alea_transition_slice_finding_count(
    const alea_transition_slice_result_t* result);
int alea_transition_slice_finding_get(
    const alea_transition_slice_result_t* result, size_t index,
    alea_transition_slice_finding_t* out_finding);
size_t alea_transition_slice_component_count(
    const alea_transition_slice_result_t* result);
int alea_transition_slice_component_get(
    const alea_transition_slice_result_t* result, size_t index,
    alea_transition_slice_component_t* out_component);
size_t alea_transition_slice_coverage_finding_count(
    const alea_transition_slice_result_t* result);
int alea_transition_slice_coverage_finding_get(
    const alea_transition_slice_result_t* result, size_t index,
    alea_transition_slice_coverage_finding_t* out_finding);
size_t alea_transition_slice_coverage_component_count(
    const alea_transition_slice_result_t* result);
int alea_transition_slice_coverage_component_get(
    const alea_transition_slice_result_t* result, size_t index,
    alea_transition_slice_coverage_component_t* out_component);
size_t alea_transition_slice_component_link_count(
    const alea_transition_slice_result_t* result);
int alea_transition_slice_component_link_get(
    const alea_transition_slice_result_t* result, size_t index,
    alea_transition_slice_component_link_t* out_link);
size_t alea_transition_slice_refinement_frontier_count(
    const alea_transition_slice_result_t* result);
int alea_transition_slice_refinement_frontier_get(
    const alea_transition_slice_result_t* result, size_t index,
    alea_transition_slice_refinement_frontier_t* out_frontier);
size_t alea_transition_slice_critical_tile_count(
    const alea_transition_slice_result_t* result);
int alea_transition_slice_critical_tile_get(
    const alea_transition_slice_result_t* result, size_t index,
    alea_transition_slice_critical_tile_t* out_tile);
size_t alea_transition_slice_critical_tile_source_count(
    const alea_transition_slice_result_t* result);
int alea_transition_slice_critical_tile_source_get(
    const alea_transition_slice_result_t* result, size_t index,
    alea_transition_slice_critical_tile_source_t* out_source);
size_t alea_transition_slice_critical_finding_count(
    const alea_transition_slice_result_t* result);
int alea_transition_slice_critical_finding_get(
    const alea_transition_slice_result_t* result, size_t index,
    alea_transition_slice_critical_finding_t* out_finding);
int alea_transition_slice_stats(
    const alea_transition_slice_result_t* result,
    alea_transition_slice_stats_t* out_stats);

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
