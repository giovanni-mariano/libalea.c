// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

#ifndef ALEA_RAYCAST_INTERNAL_H
#define ALEA_RAYCAST_INTERNAL_H

#include "alea_types.h"
#include "alea_raycast.h"
#include "alea_slice.h"
#include "alea_geo_validator.h"
#include "util/alea_vec.h"
#include "util/alea_atomic.h"
#include <stddef.h>


/**
 * @file raycast.h
 * @brief Ray casting for CSG geometry
 */

/* ============================================================================
 * DATA STRUCTURES
 * ============================================================================ */

typedef struct {
    double ox, oy, oz;
    double dx, dy, dz;
    double inv_dx, inv_dy, inv_dz;  /* 1.0/d* — IEEE 754 inf for zero components */
} alea_ray_t;

typedef struct {
    double t;
    int surface_id;        /* MCNP/user surface ID for reporting (0 synthetic, -1 none) */
    uint32_t primitive_id; /* Canonical deduplicated primitive (ALEA_PRIMITIVE_ID_INVALID if none) */
    double nx, ny, nz;
} alea_ray_hit_t;

typedef struct {
    double t_enter;
    double t_exit;
    int cell_id;
    int material_id;
    double density;
    int enter_surface_id; /* Surface crossed at t_enter: -1 none, 0 synthetic, >0 physical */
    int exit_surface_id;  /* Surface crossed at t_exit: -1 none, 0 synthetic, >0 physical */
    int enter_hit_index;  /* Index into hits[] for the surface at t_enter, or -1 */
    uint8_t resolution_flags; /* ALEA_RESOLVE_* bits (alea_types.h) */
    uint32_t path_index;  /* Result-local hierarchy path, UINT32_MAX if absent */
} alea_ray_segment_t;

/* Internal, ordered ownership boundary.  This is intentionally distinct from
 * alea_ray_hit_t: hits are mathematical intersections, whereas events are the
 * subset that changes resolved ownership (plus synthetic lattice crossings). */
typedef enum {
    ALEA_RAY_BOUNDARY_EVENT_PHYSICAL,
    ALEA_RAY_BOUNDARY_EVENT_SYNTHETIC_LATTICE,
    ALEA_RAY_BOUNDARY_EVENT_UNRESOLVED
} alea_ray_boundary_event_kind_t;

typedef struct {
    double t;
    alea_ray_boundary_event_kind_t kind;
    int surface_id;          /* -1 unresolved, 0 synthetic, >0 physical */
    uint32_t primitive_id;   /* UINT32_MAX when no physical primitive applies */
    int cell_before;
    int cell_after;
    int material_before;
    int material_after;
    uint8_t resolution_flags;
    double nx, ny, nz;       /* zero for synthetic/unresolved events */
} alea_ray_boundary_event_t;

/* Internal query contract.  This deliberately lives below the public API
 * while consumers converge on its result and budget semantics. */
typedef enum {
    ALEA_RAY_QUERY_ANY_HIT,
    ALEA_RAY_QUERY_FIRST_CELL,
    ALEA_RAY_QUERY_FIRST_VISIBLE,
    ALEA_RAY_QUERY_SEGMENTS,
    ALEA_RAY_QUERY_BOUNDARY_EVENTS
} alea_ray_query_kind_t;

/* Traversal selection is deliberately independent of query kind.  AUTO
 * preserves the legacy policy; the FAST modes use the hierarchical stepper.
 * FAST_FORWARD_REVERSE returns the forward trace and records whether it
 * disagrees with a normalized reverse trace.  The disagreement is diagnostic
 * evidence, not a traversal failure or a geometry-validity verdict. */
typedef enum {
    ALEA_RAY_QUERY_BACKEND_AUTO = 0,
    ALEA_RAY_QUERY_BACKEND_GLOBAL,
    ALEA_RAY_QUERY_BACKEND_FAST_FORWARD,
    ALEA_RAY_QUERY_BACKEND_FAST_REVERSE,
    ALEA_RAY_QUERY_BACKEND_FAST_FORWARD_REVERSE
} alea_ray_query_backend_t;

enum {
    ALEA_RAY_QUERY_FIELD_CELL_ID          = 1u << 0,
    ALEA_RAY_QUERY_FIELD_MATERIAL_ID      = 1u << 1,
    ALEA_RAY_QUERY_FIELD_DENSITY          = 1u << 2,
    ALEA_RAY_QUERY_FIELD_SURFACE_ID       = 1u << 3,
    ALEA_RAY_QUERY_FIELD_SURFACE_NORMAL   = 1u << 4,
    ALEA_RAY_QUERY_FIELD_RESOLUTION_FLAGS = 1u << 5,
    ALEA_RAY_QUERY_FIELD_PRIMITIVE_ID     = 1u << 6
};

typedef struct {
    alea_ray_query_kind_t kind;
    alea_ray_query_backend_t backend;
    uint32_t fields;
    double t_min;
    double t_max;             /* <= 0 means unbounded */
    int material_filter;      /* -1 accepts every non-void material */
    uint64_t max_events;      /* 0 means unbounded */
    uint64_t max_output_bytes;/* 0 means unbounded */
} alea_ray_query_t;

/* Ownership semantics are intentionally independent of a consumer product.
 * TRACK_COHERENT is the selected-owner policy used by the fast hierarchical
 * walker. SELECT_CANONICAL is a deterministic point/compatibility selection.
 * COMPLETE_COVERAGE is reserved for the independent diagnostic sweep and must
 * never be silently substituted by either selected-owner policy. */
typedef enum {
    ALEA_RAY_OWNERSHIP_TRACK_COHERENT,
    ALEA_RAY_OWNERSHIP_SELECT_CANONICAL,
    ALEA_RAY_OWNERSHIP_COMPLETE_COVERAGE
} alea_ray_ownership_policy_t;

typedef enum {
    ALEA_RAY_ENGINE_SELECTED_WALKER,
    ALEA_RAY_ENGINE_GLOBAL_BREAKPOINTS,
    ALEA_RAY_ENGINE_GLOBAL_COVERAGE
} alea_ray_engine_t;

/* Engine work requirements, deliberately not a public field mask. */
typedef struct {
    uint8_t need_selected_owner;
    uint8_t need_complete_coverage;
    uint8_t need_density;
    uint8_t need_surface_identity;
    uint8_t need_all_coincident_primitives;
    uint8_t need_normal;
    uint8_t need_occurrence_key;
    uint8_t need_projected_owner;
    uint8_t need_full_path;
} alea_ray_requirements_t;

/* Immutable, allocation-free lowering of an internal query descriptor.
 * Consumers use this rather than interpreting descriptor field bits in their
 * hot traversal path. */
typedef struct {
    alea_ray_query_kind_t product;
    alea_ray_query_backend_t backend;
    alea_ray_engine_t engine;
    alea_ray_ownership_policy_t ownership;
    alea_ray_requirements_t requirements;
    double t_min;
    double t_max;
    int material_filter;
    uint64_t max_events;
    uint64_t max_output_bytes;
} alea_ray_plan_t;

int alea_ray_query_lower(const alea_ray_query_t* query,
                         alea_ray_plan_t* out_plan);

typedef struct {
    bool found;
    double t;
    int cell_id;
    int material_id;
    double density;
    int surface_id;
    uint32_t primitive_id;
    uint8_t resolution_flags;
    double nx, ny, nz;
} alea_ray_first_visible_result_t;

/* Query policies publish only the scalar answer they own. Segment output is
 * returned through the supplied reusable trace; boundary output through the
 * supplied reusable event vector. */
typedef struct {
    bool any_hit;
    bool directional_mismatch;
    int first_cell_id;
    double first_cell_t;
    alea_ray_first_visible_result_t first_visible;
} alea_ray_query_output_t;

typedef struct {
    uint32_t offset;
    uint16_t count;
} alea_ray_path_t;

typedef struct {
    uint32_t cell_index;
    int cell_id;
    int material_id;
    int universe_id;
    int fill_universe;
    int depth;
    uint8_t is_lattice;
    double lattice_origin[3];
    uint64_t occurrence_key;
} alea_ray_path_entry_t;

ALEA_VEC_DEFINE(alea_ray_hit_vec, alea_ray_hit_t);
ALEA_VEC_DEFINE(alea_ray_segment_vec, alea_ray_segment_t);
ALEA_VEC_DEFINE(alea_ray_boundary_event_vec, alea_ray_boundary_event_t);
ALEA_VEC_DEFINE(alea_ray_path_vec, alea_ray_path_t);
ALEA_VEC_DEFINE(alea_ray_path_entry_vec, alea_ray_path_entry_t);

/* Number of histogram bins for primitive-type instrumentation. Sized to
 * comfortably exceed the alea_primitive_type_t range (currently <= 22). */
#define ALEA_RAYCAST_PRIM_TYPE_BINS 32

/* Bounded, diagnostic-only attribution of work performed to find enclosing
 * hierarchy boundaries.  Sixteen cells is sufficient to expose the dominant
 * owners in a single ray without allocating or maintaining a per-system map
 * on the ray-walk hot path. */
#define ALEA_RAYCAST_ANCESTOR_HOT_CELLS 16

typedef struct {
    uint32_t cell_index;
    int cell_id;
    uint64_t surface_tests;
    uint64_t queries;
    uint64_t winning_events;
} alea_raycast_ancestor_cell_stat_t;

/* Use struct tag matching the public API forward declaration */
struct alea_raycast_result {
    alea_ray_t ray;
    alea_ray_hit_vec_t hits;
    alea_ray_segment_vec_t segments;
    alea_ray_path_vec_t paths;
    alea_ray_path_entry_vec_t path_entries;
    uint8_t capture_paths;
    /* Optional batch-wide live segment budget. The public single-ray APIs
     * leave these unset; the compact batch API reserves before every append. */
    atomic_uint_fast64_t* segment_counter;
    uint64_t segment_limit;
    uint8_t segment_limit_exceeded;
    /* Optional batch-wide budget for the flattened full-path CSR output. */
    atomic_uint_fast64_t* path_entry_counter;
    uint64_t path_entry_limit;
    uint8_t path_entry_limit_exceeded;
    int surfaces_tested;
    int bbox_culled;
    int point_lookups;
    int step_iterations;
    size_t blas_placement_candidates;
    size_t blas_placements_pruned;
    size_t blas_universe_queries;
    size_t blas_cell_candidates;
    size_t blas_cells_tested;
    size_t blas_hits_before_dedup;

    /* Selected-walker ownership and optional-work attribution.  These count
     * coarse decisions only; no timers or allocator hooks enter the hot loop. */
    uint64_t owner_neighbor_attempts;
    uint64_t owner_neighbor_hits;
    uint64_t owner_path_attempts;
    uint64_t owner_path_hits;
    uint64_t owner_root_queries;
    uint64_t owner_root_hits;
    uint64_t owner_full_queries;
    uint64_t owner_full_hits;
    uint64_t boundary_event_enrichments;
    uint64_t path_snapshot_copies;
    uint64_t path_snapshot_entries;
    uint64_t selected_intervals_yielded;
    uint64_t result_buffer_growths;
    uint64_t result_buffer_growth_bytes;

    /* Phase 2 per-ray surface-test instrumentation (cheap, always-on).
     * Attributes the surfaces_tested total across the cell-aware stepper's
     * sources and characterises the crossed cells, so an acceleration
     * strategy (per-cell BVH vs type-specialised loops) can be chosen from
     * measured data rather than guessed. */
    int terminal_surfaces_tested;  /* tests against the terminal cell's surfaces */
    int lattice_surfaces_tested;   /* tests against lattice wrapper surfaces */
    int ancestor_surfaces_tested;  /* tests against active-path ancestor surfaces */
    size_t crossed_cell_count;       /* stepper iterations landing in a surfaced cell */
    uint32_t max_cell_surface_count; /* largest crossed surface_index_count */
    uint64_t sum_cell_surface_count; /* sum of crossed surface_index_count (for avg) */
    uint32_t prim_type_tests[ALEA_RAYCAST_PRIM_TYPE_BINS]; /* tests by primitive type */

    /* Ancestor-boundary diagnostics.  These are deliberately bounded: they
     * identify the first distinct active ancestors seen by a ray and account
     * any remaining ones in the unattributed totals.  They are evidence for
     * choosing an acceleration strategy, not a correctness cache key. */
    uint64_t ancestor_surface_queries;
    uint64_t ancestor_path_entries_examined;
    uint32_t ancestor_max_path_depth;
    uint64_t ancestor_unattributed_surface_tests;
    uint64_t ancestor_unattributed_queries;
    uint64_t ancestor_unattributed_winning_events;
    alea_raycast_ancestor_cell_stat_t
        ancestor_hot_cells[ALEA_RAYCAST_ANCESTOR_HOT_CELLS];

    /* Private lattice-entry work counters. These identify work performed while
     * the hierarchical walker is unresolved (void/container), rather than the
     * DDA/surface work of an already-resolved lattice segment. */
    uint64_t lattice_entry_calls;
    uint64_t lattice_entry_tlas_nodes_tested;
    uint64_t lattice_entry_tlas_leaves_visited;
    uint64_t lattice_entry_candidates;
    uint64_t lattice_entry_dda_steps;
    uint64_t lattice_entry_no_entry_results;
    uint64_t lattice_entry_future_entry_results;
    uint64_t lattice_entry_already_inside_results;
    uint64_t lattice_entry_ancestor_surface_tests;
    uint64_t lattice_entry_ancestor_events;
    uint64_t lattice_entry_canonical_rejections;
};

/* Typedef for internal use */
typedef struct alea_raycast_result alea_raycast_result_t;

/* Internal canonical event cache for a slice view.  It is owned by the slice
 * module and intentionally opaque; validator/provenance users share it only
 * through identity-checked accessors. */
typedef struct alea_slice_directional_event_cache alea_slice_directional_event_cache_t;

alea_slice_directional_event_cache_t* alea_slice_directional_event_cache_create(
    alea_system_t* sys, const alea_slice_view_t* view, int width, int height);
void alea_slice_directional_event_cache_destroy(
    alea_slice_directional_event_cache_t* cache);
int alea_slice_directional_event_cache_matches(
    const alea_slice_directional_event_cache_t* cache,
    const alea_system_t* sys, const alea_slice_view_t* view,
    int width, int height);
int alea_slice_directional_event_cache_line_events(
    const alea_slice_directional_event_cache_t* cache,
    alea_slice_edge_orientation_t orientation, int reverse, size_t line,
    const alea_ray_boundary_event_t** out_events, size_t* out_count);
int alea_slice_directional_event_cache_dimensions(
    const alea_slice_directional_event_cache_t* cache,
    int* out_width, int* out_height);
/* Cache-contract inspection and ownership-trace reuse are private to slice
 * consumers.  The ownership contract retains occurrence keys at the declared
 * projected depth.  A nonzero complete result means the requested semantics
 * were fully materialized (rather than merely having an event stream available). */
int alea_slice_directional_event_cache_contract(
    const alea_slice_directional_event_cache_t* cache,
    alea_slice_edge_orientation_t orientation, int reverse,
    uint32_t required_event_fields, int projected_depth,
    uint64_t max_events, uint64_t max_output_bytes, int* out_complete);
const alea_raycast_batch_result_t*
alea_slice_directional_event_cache_ownership_trace(
    const alea_slice_directional_event_cache_t* cache,
    alea_slice_edge_orientation_t orientation, int reverse,
    int projected_depth, uint32_t required_fields);
int alea_slice_surface_boundary_map_create_with_event_cache(
    alea_system_t* sys, const alea_slice_view_t* view, int width, int height,
    const int* grid_ids, alea_slice_classify_point_fn classify,
    void* classify_userdata, const alea_slice_directional_event_cache_t* cache,
    alea_slice_surface_boundary_map_t** out_map);

/* Private Phase 6 validator extension.  A matching canonical event cache
 * augments mismatch intervals with boundary evidence; it never changes the
 * ownership comparison itself. */
int alea_validate_ray_slice_compact_with_event_cache(
    alea_system_t* sys, const alea_slice_view_t* view, size_t row_count,
    const alea_ray_slice_validation_options_t* validation_options,
    const alea_raycast_batch_options_t* render_options,
    alea_raycast_batch_result_t* inout_fast_forward,
    const alea_slice_directional_event_cache_t* event_cache,
    alea_ray_slice_validation_result_t* out_validation);

const int32_t* alea_ray_slice_validation_u_enter_forward_surface_ids_internal(
    const alea_ray_slice_validation_result_t* result);
const int32_t* alea_ray_slice_validation_u_enter_reverse_surface_ids_internal(
    const alea_ray_slice_validation_result_t* result);
const int32_t* alea_ray_slice_validation_u_exit_forward_surface_ids_internal(
    const alea_ray_slice_validation_result_t* result);
const int32_t* alea_ray_slice_validation_u_exit_reverse_surface_ids_internal(
    const alea_ray_slice_validation_result_t* result);
const uint32_t* alea_ray_slice_validation_u_enter_provenance_flags_internal(
    const alea_ray_slice_validation_result_t* result);
const uint32_t* alea_ray_slice_validation_u_exit_provenance_flags_internal(
    const alea_ray_slice_validation_result_t* result);

/* ============================================================================
 * MAIN API
 * ============================================================================ */

/**
 * @brief Initialize a ray
 *
 * Direction is automatically normalized. Precomputes inverse direction
 * for efficient slab tests (IEEE 754 inf for zero components).
 *
 * @return 0 on success, -1 if direction is zero-length (fallback (0,0,1) used)
 */
int alea_ray_init(alea_ray_t* ray,
                  double ox, double oy, double oz,
                  double dx, double dy, double dz);

/**
 * @brief Initialize a ray with pre-normalized direction (skip sqrt + divides)
 *
 * Use when the direction is already unit-length (e.g., camera rays).
 * Caller must ensure |dx,dy,dz| = 1 and direction is non-zero.
 */
void alea_ray_init_normalized(alea_ray_t* ray,
                              double ox, double oy, double oz,
                              double dx, double dy, double dz);

/**
 * @brief Build all raycast caches (BVH, spatial index, cell adjacency)
 *
 * Call once before tracing rays. Avoids lazy-build warnings.
 */
int alea_raycast_ensure_caches(alea_system_t* sys);

/**
 * @brief Build raycast caches for the hierarchical spatial path.
 *
 * This prepares surface/cell caches plus the hierarchical spatial index, but
 * intentionally does not build the flat spatial instance index.
 */
int alea_raycast_ensure_hier_caches(alea_system_t* sys);

/**
 * @brief Initialize raycast result (call before first use)
 */
void alea_raycast_result_init(alea_raycast_result_t* result);

/**
 * @brief Free raycast result memory
 */
void alea_raycast_result_free(alea_raycast_result_t* result);

/**
 * @brief Clear result for reuse (keeps allocated memory)
 */
void alea_raycast_result_clear(alea_raycast_result_t* result);

/**
 * @brief Pre-allocate result buffers to avoid per-ray realloc
 */
void alea_raycast_result_reserve(alea_raycast_result_t* result,
                                size_t hit_cap, size_t seg_cap);

/* Internal compact slice-cache contract.  The public batch result remains
 * opaque; these helpers are consumed by the geo-validator implementation. */
int alea_raycast_batch_result_matches_fast_slice_cache(
    const alea_raycast_batch_result_t* result,
    const alea_system_t* sys,
    const alea_slice_view_t* view,
    size_t row_count,
    const alea_raycast_batch_options_t* render_options,
    int projected_depth);
int alea_raycast_batch_result_get_compact_slice_provenance(
    const alea_raycast_batch_result_t* result,
    const alea_slice_view_t* view,
    size_t row_count,
    int* out_projected_depth);
void alea_raycast_batch_result_swap_internal(
    alea_raycast_batch_result_t* a,
    alea_raycast_batch_result_t* b);

/* Private batch-level maxima of per-ray lattice-entry work.  This preserves
 * the public batch ABI while allowing benchmark and diagnostic code to expose
 * outlier rays instead of hiding them behind aggregate throughput. */
typedef struct {
    uint64_t max_owner_neighbor_attempts;
    uint64_t max_owner_neighbor_hits;
    uint64_t max_owner_path_attempts;
    uint64_t max_owner_path_hits;
    uint64_t max_owner_root_queries;
    uint64_t max_owner_root_hits;
    uint64_t max_owner_full_queries;
    uint64_t max_owner_full_hits;
    uint64_t max_boundary_event_enrichments;
    uint64_t max_path_snapshot_copies;
    uint64_t max_path_snapshot_entries;
    uint64_t max_selected_intervals_yielded;
    uint64_t max_result_buffer_growths;
    uint64_t max_result_buffer_growth_bytes;
    uint64_t max_lattice_entry_calls;
    uint64_t max_lattice_entry_tlas_nodes_tested;
    uint64_t max_lattice_entry_tlas_leaves_visited;
    uint64_t max_lattice_entry_candidates;
    uint64_t max_lattice_entry_dda_steps;
    uint64_t max_lattice_entry_no_entry_results;
    uint64_t max_lattice_entry_future_entry_results;
    uint64_t max_lattice_entry_already_inside_results;
    uint64_t max_lattice_entry_ancestor_surface_tests;
    uint64_t max_lattice_entry_ancestor_events;
    uint64_t max_lattice_entry_canonical_rejections;
} alea_raycast_batch_work_stats_t;

int alea_raycast_batch_result_get_work_stats_internal(
    const alea_raycast_batch_result_t* result,
    alea_raycast_batch_work_stats_t* out_stats);

/* Phase 5 internal compact executor.  Public alea_raycast_hier_batch() is a
 * scalar-t_max SEGMENTS adapter.  Ranges are per input ray; t_max <= 0 means
 * unbounded, matching the scalar API.  Only SEGMENTS is materialized in this
 * first increment; other query kinds remain deliberately unsupported rather
 * than silently producing a segment CSR with different semantics. */
typedef struct {
    alea_ray_query_kind_t kind;
    uint32_t fields;
    int material_filter;
    const double* t_mins;
    const double* t_maxs;
    uint64_t max_events;
    uint64_t max_output_bytes;
} alea_ray_batch_query_t;

int alea_raycast_hier_batch_query_nocache(
    alea_system_t* sys, const double* origins_xyz, const double* directions_xyz,
    size_t ray_count, const alea_ray_batch_query_t* query,
    const alea_raycast_batch_options_t* options,
    alea_raycast_batch_result_t* result);

/* Private SoA writer for FIRST_VISIBLE batch queries.  Cell/material/distance
 * are always present; surface, normal, and resolution arrays are allocated
 * only when requested by alea_ray_batch_query_t.fields. */
typedef struct {
    size_t ray_count;
    uint8_t* found;
    double* t;
    int32_t* cell_ids;
    int32_t* material_ids;
    double* densities;
    int32_t* surface_ids;
    uint32_t* primitive_ids;
    uint8_t* resolution_flags;
    double* normals_xyz;
} alea_ray_first_visible_batch_result_t;

void alea_ray_first_visible_batch_result_init(
    alea_ray_first_visible_batch_result_t* result);
void alea_ray_first_visible_batch_result_free(
    alea_ray_first_visible_batch_result_t* result);
int alea_raycast_hier_first_visible_batch_nocache(
    alea_system_t* sys, const double* origins_xyz, const double* directions_xyz,
    size_t ray_count, const alea_ray_batch_query_t* query,
    alea_ray_first_visible_batch_result_t* result);

/* Internal packet executor. The API layer owns validation/allocation and
 * checks statuses after the parallel region; this function owns only bounded
 * lane traversal and direct SoA writes. */
int alea_raycast_hier_first_visible_batch_execute_nocache(
    alea_system_t* sys, const double* origins_xyz, const double* directions_xyz,
    size_t ray_count, const alea_ray_batch_query_t* query,
    alea_ray_first_visible_batch_result_t* result, int* statuses);

/** Private compact writer for ANY_HIT batch queries. */
typedef struct {
    size_t ray_count;
    uint8_t* hits;
} alea_ray_any_hit_batch_result_t;

void alea_ray_any_hit_batch_result_init(alea_ray_any_hit_batch_result_t* result);
void alea_ray_any_hit_batch_result_free(alea_ray_any_hit_batch_result_t* result);
int alea_raycast_hier_any_hit_batch_nocache(
    alea_system_t* sys, const double* origins_xyz, const double* directions_xyz,
    size_t ray_count, const alea_ray_batch_query_t* query,
    alea_ray_any_hit_batch_result_t* result);
int alea_raycast_hier_any_hit_batch_execute_nocache(
    alea_system_t* sys, const double* origins_xyz, const double* directions_xyz,
    size_t ray_count, const alea_ray_batch_query_t* query,
    alea_ray_any_hit_batch_result_t* result, int* statuses);

/* Private CSR writer for the canonical BOUNDARY_EVENTS contract. */
typedef struct {
    size_t ray_count;
    size_t event_count;
    uint64_t* ray_offsets;
    double* t;
    uint8_t* kinds;
    int32_t* surface_ids;
    int32_t* cell_before;
    int32_t* cell_after;
    int32_t* material_before;
    int32_t* material_after;
    uint8_t* resolution_flags;
    uint32_t* primitive_ids;
    double* normals_xyz;
} alea_ray_boundary_event_batch_result_t;

void alea_ray_boundary_event_batch_result_init(
    alea_ray_boundary_event_batch_result_t* result);
void alea_ray_boundary_event_batch_result_free(
    alea_ray_boundary_event_batch_result_t* result);
int alea_raycast_boundary_events_batch_nocache(
    alea_system_t* sys, const double* origins_xyz, const double* directions_xyz,
    size_t ray_count, const alea_ray_batch_query_t* query,
    alea_ray_boundary_event_batch_result_t* result);

/**
 * @brief Cast ray and find all surface intersections
 *
 * Finds where ray intersects each surface in the geometry.
 * Results are sorted by distance.
 *
 * @param sys CSG system
 * @param ray Ray to cast
 * @param t_min Minimum distance (usually 0 or small epsilon)
 * @param t_max Maximum distance
 * @param result Output: surface hits
 * @return 0 on success, -1 on error
 */
int alea_raycast_surfaces(alea_system_t* sys,
                         const alea_ray_t* ray,
                         double t_min, double t_max,
                         alea_raycast_result_t* result);

/**
 * @brief Cast ray and find surface intersections (skip cache check)
 *
 * Like alea_raycast_surfaces but assumes caches are already built.
 * Use from tight loops (render, volume estimation) after calling
 * alea_raycast_ensure_caches() once up front.
 */
int alea_raycast_surfaces_nocache(alea_system_t* sys,
                                  const alea_ray_t* ray,
                                  double t_min, double t_max,
                                  alea_raycast_result_t* result);

/**
 * @brief Convert surface hits to cell segments
 *
 * Takes raw surface intersections and determines which cell
 * the ray passes through between each pair of hits.
 *
 * @param sys CSG system
 * @param result Raycast result (hits must be populated)
 * @return 0 on success, -1 on error
 */
int alea_raycast_to_segments(alea_system_t* sys,
                            double t_max,
                            alea_raycast_result_t* result);

/**
 * @brief Full raycast: surfaces + segments in one call
 *
 * Convenience function that does alea_raycast_surfaces + alea_raycast_to_segments.
 */
int alea_raycast(alea_system_t* sys,
                double ox, double oy, double oz,
                double dx, double dy, double dz,
                double t_max,
                alea_raycast_result_t* result);

/**
 * @brief Fast hierarchical material/path segment raycast.
 *
 * Segment output is the primary contract. The complete ordered surface-hit list
 * is intentionally not reconstructed on this fast path; use segment
 * enter/exit surface IDs for boundary reporting.
 */
int alea_raycast_hier(alea_system_t* sys,
                      double ox, double oy, double oz,
                      double dx, double dy, double dz,
                      double t_max,
                      alea_raycast_result_t* result);

/**
 * @brief Deprecated alias for alea_raycast_hier_fast_segments().
 *
 * Kept for source compatibility with earlier experiments.
 */
int alea_raycast_hier_cell_aware(alea_system_t* sys,
                                 double ox, double oy, double oz,
                                 double dx, double dy, double dz,
                                 double t_max,
                                 alea_raycast_result_t* result);

/**
 * @brief Fast hierarchical material/path segment raycast.
 *
 * Segment output is the primary contract. The complete ordered surface-hit list
 * is intentionally not reconstructed on this fast path.
 */
int alea_raycast_hier_fast_segments(alea_system_t* sys,
                                    double ox, double oy, double oz,
                                    double dx, double dy, double dz,
                                    double t_max,
                                    alea_raycast_result_t* result);

int alea_raycast_hier_blas_experimental(alea_system_t* sys,
                                        double ox, double oy, double oz,
                                        double dx, double dy, double dz,
                                        double t_max,
                                        alea_raycast_result_t* result);

/**
 * @brief Hierarchical raycast that also emits boundary surface hits.
 *
 * Returns the same material/path segments as alea_raycast_hier_fast_segments()
 * (segment parity), and additionally populates result->hits with one
 * alea_ray_hit_t per physical surface crossing along the stepped path. Each
 * segment's enter_hit_index links back to the hit at its entry boundary
 * (-1 when that boundary has no physical surface, e.g. synthetic lattice DDA
 * boundaries or the initial void entry).
 *
 * Normals are reported in world space. For non-lattice/non-fill geometry the
 * frame transform is identity, so hits match the flat alea_raycast() path.
 * This path may be slower than fast_segments but remains step-based rather
 * than candidate-list based.
 */
int alea_raycast_hier_with_hits(alea_system_t* sys,
                                double ox, double oy, double oz,
                                double dx, double dy, double dz,
                                double t_max,
                                alea_raycast_result_t* result);

/**
 * @brief Buffer-reuse hierarchical raycast with boundary hits.
 *
 * Same contract as alea_raycast_hier_with_hits() (segment + hit parity), but
 * takes a pre-normalized ray, assumes query caches are already built
 * (ALEA_CACHE_RAYCAST), and does NOT free the result buffers. The caller must
 * call alea_raycast_result_clear() between rays. Intended for per-pixel render
 * loops that reuse a thread-local result to avoid malloc/free churn.
 */
int alea_raycast_hier_with_hits_nocache(alea_system_t* sys,
                                        const alea_ray_t* ray,
                                        double t_max,
                                        alea_raycast_result_t* result);

/**
 * Hierarchical first-visible query. Stops after resolving the first qualifying
 * material interval (including lattice DDA/fill traversal) and does not append
 * hit or segment vectors. Caches must be prepared; result is scratch only.
 */
int alea_raycast_hier_first_visible_nocache(
    alea_system_t* sys, const alea_ray_t* ray, double t_min, double t_max,
    int material_filter, int include_normal, alea_raycast_result_t* scratch,
    alea_ray_first_visible_result_t* out_visible);

/* Resolve the first owned cell interval without materializing hits or
 * segments. Supports hierarchy, fills, and lattice DDA. */
int alea_raycast_hier_first_cell_nocache(
    alea_system_t* sys, const alea_ray_t* ray, double t_min, double t_max,
    int material_filter,
    alea_raycast_result_t* scratch, int* out_cell_id, double* out_t);

/** Hierarchical occlusion policy. Stops at the first non-void material
 * interval and leaves scratch hit/segment/path vectors empty. */
int alea_raycast_hier_any_hit_nocache(
    alea_system_t* sys, const alea_ray_t* ray, double t_min, double t_max,
    int material_filter, alea_raycast_result_t* scratch, int* out_hit);

/**
 * @brief Buffer-reuse hierarchical raycast, segments only (no hit list).
 *
 * Like alea_raycast_hier_with_hits_nocache() but skips boundary-hit
 * reconstruction; produces material/path segments only. Use when surface
 * normals are not needed (e.g. x-ray accumulation).
 */
int alea_raycast_hier_segments_nocache(alea_system_t* sys,
                                       const alea_ray_t* ray,
                                       double t_max,
                                       alea_raycast_result_t* result);

/* Internal scalar streaming adapter for consumers that need selected material
 * intervals but do not need to publish a segment vector.  Return a positive
 * value from the callback to stop successfully; return a negative value to
 * abort the walk as an error. */
typedef int (*alea_raycast_selected_segment_callback_t)(
    void* context, const alea_ray_segment_t* segment);

int alea_raycast_hier_visit_segments_nocache(
    alea_system_t* sys, const alea_ray_t* ray, double t_max,
    alea_raycast_result_t* scratch,
    alea_raycast_selected_segment_callback_t callback, void* context);

/* Internal fixed-output batch companion to the scalar visitor.  The callback
 * may run concurrently for different ray_index values and must therefore use
 * per-ray output slots or otherwise provide its own synchronization. */
typedef int (*alea_raycast_batch_selected_segment_callback_t)(
    void* context, size_t ray_index, const alea_ray_segment_t* segment);

int alea_raycast_hier_visit_segments_batch_nocache(
    alea_system_t* sys, const double* origins_xyz, const double* directions_xyz,
    size_t ray_count, double t_max,
    alea_raycast_batch_selected_segment_callback_t callback, void* context);

/**
 * Buffer-reuse canonical global raycast.
 *
 * Matches alea_raycast()'s global surface/fill/lattice pipeline but clears
 * logical result contents rather than freeing vector capacity. The caller
 * supplies a normalized alea_ray_t and must have prepared raycast caches.
 * Intended for high-frequency internal loops that require canonical lattice
 * semantics without per-ray allocation churn.
 */
int alea_raycast_global_reuse_nocache(alea_system_t* sys,
                                      const alea_ray_t* ray,
                                      double t_max,
                                      alea_raycast_result_t* result);

/**
 * Enumerate the global diagnostic breakpoint set without constructing
 * selected-owner segments.  The returned hit vector contains physical,
 * fill-transformed, and synthetic lattice crossings in deterministic order.
 * This is deliberately separate from the selected-owner walker: callers that
 * need complete coverage must not accidentally inherit deck-precedence
 * segments from a production trace.
 *
 * Caches must already be prepared and result is caller-owned reusable
 * scratch.  Only result->hits is materialized.
 */
int alea_raycast_global_breakpoints_reuse_nocache(
    alea_system_t* sys, const alea_ray_t* ray, double t_max,
    alea_raycast_result_t* result);

/**
 * Complete-coverage diagnostic sweep behind the legacy public interval
 * classifier.  It uses global breakpoints and recursive all-owner point
 * resolution for each open interval; it never consumes selected segments.
 * Caches and the universe index must already be prepared.
 */
int alea_ray_coverage_classify_reuse_nocache(
    alea_system_t* sys, const alea_ray_t* ray, double t_max,
    alea_raycast_result_t* breakpoint_scratch,
    alea_ray_interval_finding_t* out, size_t max_out);

typedef enum {
    ALEA_RAY_COVERAGE_UNIQUE,
    ALEA_RAY_COVERAGE_GAP,
    ALEA_RAY_COVERAGE_ALLOWED_EXTERIOR,
    ALEA_RAY_COVERAGE_OVERLAP,
    ALEA_RAY_COVERAGE_UNDEFINED_FILL,
    /* The ownership records cannot be normalized into complete chains. */
    ALEA_RAY_COVERAGE_UNRESOLVED,
    /* The owner set exceeded the diagnostic budget.  owners contains the
     * retained prefix only and must not be treated as a complete set. */
    ALEA_RAY_COVERAGE_TRUNCATED
} alea_ray_coverage_kind_t;

typedef struct {
    int cell_id;
    int cell_index;
    int material_id;
    int universe_id;
    int fill_universe;
    int depth;
    uint64_t occurrence_key;
    uint64_t parent_occurrence_key; /* zero = root-universe claimant */
    uint8_t resolution_flags;
} alea_ray_coverage_owner_t;

/* The owners span is valid only for the duration of the callback. */
typedef struct {
    double t_enter;
    double t_exit;
    alea_ray_coverage_kind_t kind;
    const alea_ray_coverage_owner_t* owners;
    size_t owner_count;             /* retained owner records */
    size_t owner_count_lower_bound; /* exact unless kind is TRUNCATED */
} alea_ray_coverage_interval_t;

typedef int (*alea_ray_coverage_interval_callback_t)(
    void* context, const alea_ray_coverage_interval_t* interval);

/* Internal diagnostic-domain contract.  The observation interval is always
 * [0, t_max] for the current scalar sweep.  When has_domain is set, only the
 * closed t-domain requires ownership; unowned portions outside it are either
 * omitted or emitted as ALLOWED_EXTERIOR. */
typedef struct {
    double t_min;
    double t_max;
    uint8_t has_domain;
    uint8_t report_allowed_exterior;
} alea_ray_coverage_domain_t;

/* Stream merged complete-coverage intervals without allocating public output.
 * The callback may stop the sweep by returning nonzero. */
int alea_ray_coverage_sweep_reuse_nocache(
    alea_system_t* sys, const alea_ray_t* ray, double t_max,
    alea_raycast_result_t* breakpoint_scratch,
    alea_ray_coverage_interval_callback_t callback, void* context);

/* Domain-aware internal coverage sweep.  The supplied domain is intersected
 * with [0, t_max]; its entry/exit become interval boundaries even when they
 * are not geometry crossings. */
int alea_ray_coverage_sweep_domain_reuse_nocache(
    alea_system_t* sys, const alea_ray_t* ray, double t_max,
    const alea_ray_coverage_domain_t* domain,
    alea_raycast_result_t* breakpoint_scratch,
    alea_ray_coverage_interval_callback_t callback, void* context);

/* One serial coverage-row specification.  direction_tag and
 * transverse_coordinate are caller-owned provenance used by slice consumers;
 * this scalar adapter preserves input row order and does not interpret them. */
typedef struct {
    alea_ray_t ray;
    double t_max;
    alea_ray_coverage_domain_t domain;
    uint8_t use_domain;
    uint8_t direction_tag;
    double transverse_coordinate;
} alea_ray_coverage_row_t;

typedef int (*alea_ray_coverage_row_interval_callback_t)(
    void* context, size_t row_index,
    const alea_ray_coverage_interval_t* interval);

/* Run complete coverage for ordered slice rows with one caller-owned reusable
 * breakpoint scratch result.  The callback span is valid only until the next
 * interval.  Any nonzero callback return aborts the operation with -1. */
int alea_ray_coverage_rows_serial_reuse_nocache(
    alea_system_t* sys, const alea_ray_coverage_row_t* rows, size_t row_count,
    alea_raycast_result_t* breakpoint_scratch,
    alea_ray_coverage_row_interval_callback_t callback, void* context);

/* Compact scalar coverage-slice publication.  These are internal result
 * types until the Phase 11 installed API stabilizes.  All arrays are owned by
 * the result and are published transactionally by the serial builder. */
typedef struct {
    size_t max_rows;
    size_t max_intervals;
    size_t max_owners;
    size_t max_bytes;
} alea_ray_coverage_slice_limits_t;

typedef enum {
    ALEA_RAY_COVERAGE_REFINEMENT_COMPLETE = 0,
    ALEA_RAY_COVERAGE_REFINEMENT_MAX_DEPTH,
    ALEA_RAY_COVERAGE_REFINEMENT_MAX_ROWS,
    /* Candidate pairs remained, but splitting them would have produced rows
     * closer than the configured minimum transverse spacing. */
    ALEA_RAY_COVERAGE_REFINEMENT_MIN_SPACING
} alea_ray_coverage_refinement_status_t;

typedef struct {
    size_t row_count;
    size_t interval_count;
    size_t owner_count;
    alea_ray_coverage_refinement_status_t refinement_status;
    size_t* row_offsets;              /* row_count + 1 */
    uint8_t* row_direction_tags;      /* row_count */
    double* row_transverse_coordinates; /* row_count */
    double* t_enter;                  /* interval_count */
    double* t_exit;                   /* interval_count */
    uint8_t* kinds;                   /* alea_ray_coverage_kind_t */
    size_t* owner_offsets;            /* interval_count + 1 */
    size_t* owner_count_lower_bounds; /* interval_count */
    int* owner_cell_ids;
    int* owner_cell_indices;
    int* owner_material_ids;
    int* owner_universe_ids;
    int* owner_fill_universes;
    int* owner_depths;
    uint64_t* owner_occurrence_keys;
    uint64_t* owner_parent_occurrence_keys;
    uint8_t* owner_resolution_flags;
} alea_ray_coverage_slice_result_t;

void alea_ray_coverage_slice_limits_init(
    alea_ray_coverage_slice_limits_t* limits);
void alea_ray_coverage_slice_result_init(
    alea_ray_coverage_slice_result_t* result);
void alea_ray_coverage_slice_result_free(
    alea_ray_coverage_slice_result_t* result);

/* Build scalar row/interval/owner CSR output.  On failure, result retains its
 * previous successful publication. */
int alea_ray_coverage_slice_build_serial_nocache(
    alea_system_t* sys, const alea_ray_coverage_row_t* rows, size_t row_count,
    const alea_ray_coverage_slice_limits_t* limits,
    alea_raycast_result_t* breakpoint_scratch,
    alea_ray_coverage_slice_result_t* result);

/* Compare the adaptive-refinement signature of two published rows.  Endpoints
 * and row coordinates are intentionally excluded; kind, owner-chain identity,
 * retained owner count, truncation lower bound, and resolution flags are not.
 * Returns 1 when equal, 0 when different, and -1 for malformed input. */
int alea_ray_coverage_slice_rows_same_signature(
    const alea_ray_coverage_slice_result_t* result,
    size_t first_row, size_t second_row);

/* Adaptive refinement probe signals.  Signature difference alone cannot see a
 * defect that displaces boundaries without changing owner identity, nor a
 * dense or already-suspect region whose neighbours happen to agree. */
#define ALEA_RAY_COVERAGE_REFINE_SIGNATURE    (1u << 0)
#define ALEA_RAY_COVERAGE_REFINE_DISPLACEMENT (1u << 1)
#define ALEA_RAY_COVERAGE_REFINE_DENSITY      (1u << 2)
#define ALEA_RAY_COVERAGE_REFINE_FINDING      (1u << 3)

/* Probe-selection policy.  This governs which rows are sampled; it never
 * changes complete-coverage classification on any sampled ray. */
typedef struct {
    uint32_t signals;
    /* Absolute transverse distance below which a pair is never split.  A pair
     * spanning g is split into gaps of g/2, so a pair is suppressed when
     * g < 2 * min_transverse_spacing.  Zero disables the limit. */
    double min_transverse_spacing;
    /* Absolute scan-coordinate tolerance for ALEA_RAY_COVERAGE_REFINE_
     * DISPLACEMENT.  Paired endpoints further apart than this refine even when
     * owner identity matches.  Must be > 0 when that signal is selected. */
    double endpoint_displacement;
    /* Interval count at or above which a row is considered dense for
     * ALEA_RAY_COVERAGE_REFINE_DENSITY.  Must be > 0 when selected. */
    size_t crossing_density;
} alea_ray_coverage_refinement_policy_t;

/* Defaults to signature-only refinement with no spacing limit, matching the
 * behaviour of alea_ray_coverage_slice_mark_refinement_boundaries(). */
void alea_ray_coverage_refinement_policy_init(
    alea_ray_coverage_refinement_policy_t* policy);

/* Mark adjacent same-direction row pairs selected by the policy.  The output
 * has row_count - 1 entries and is ordered with the published rows.  A marked
 * pair is a deterministic candidate for a later midpoint refinement wave;
 * this helper never generates rays or changes classification.  When
 * out_spacing_limited is non-NULL it receives the number of pairs a signal
 * selected but the minimum transverse spacing suppressed.  Returns the number
 * of marked pairs, or -1 for malformed/non-monotonic row provenance or an
 * inconsistent policy. */
int alea_ray_coverage_slice_mark_refinement_boundaries_policy(
    const alea_ray_coverage_slice_result_t* result,
    const alea_ray_coverage_refinement_policy_t* policy,
    uint8_t* out_refine_between, size_t* out_spacing_limited);

/* Signature-only adapter over the policy marker. */
int alea_ray_coverage_slice_mark_refinement_boundaries(
    const alea_ray_coverage_slice_result_t* result,
    uint8_t* out_refine_between);

/* Materialize one deterministic midpoint-refinement wave selected by
 * out_refine_between.  Input rows and markers are in existing published order;
 * marked pairs must be adjacent rows of one direction with increasing
 * transverse coordinates.  max_rows == 0 is unlimited.  The caller supplies
 * all output storage, and failure leaves it untouched. */
int alea_ray_coverage_rows_refine_midpoints(
    const alea_ray_coverage_row_t* rows, size_t row_count,
    const uint8_t* refine_between, size_t max_rows,
    alea_ray_coverage_row_t* out_rows, size_t out_capacity,
    size_t* out_row_count);

/* Run bounded deterministic refinement waves serially under an explicit probe
 * policy.  Reaching the depth, row-selection, or minimum-spacing limit
 * publishes the completed sampled rows with the matching refinement_status; an
 * output/materialization failure leaves result intact.  A NULL policy selects
 * signature-only refinement. */
int alea_ray_coverage_slice_build_adaptive_policy_serial_nocache(
    alea_system_t* sys, const alea_ray_coverage_row_t* initial_rows,
    size_t initial_row_count, size_t max_refinement_depth,
    const alea_ray_coverage_refinement_policy_t* policy,
    const alea_ray_coverage_slice_limits_t* limits,
    alea_raycast_result_t* breakpoint_scratch,
    alea_ray_coverage_slice_result_t* result);

/* Signature-only adapter over the policy builder. */
int alea_ray_coverage_slice_build_adaptive_serial_nocache(
    alea_system_t* sys, const alea_ray_coverage_row_t* initial_rows,
    size_t initial_row_count, size_t max_refinement_depth,
    const alea_ray_coverage_slice_limits_t* limits,
    alea_raycast_result_t* breakpoint_scratch,
    alea_ray_coverage_slice_result_t* result);

/* Phase 10 executor storage.  One worker owns reusable scalar breakpoint
 * scratch and a variable-output coverage arena.  The arena is internal staging
 * only: callers receive CSR output after deterministic transactional
 * compaction, never a worker-local layout. */
typedef struct {
    alea_raycast_result_t breakpoint_scratch;
    alea_ray_coverage_slice_result_t arena;
} alea_ray_coverage_worker_scratch_t;

typedef struct {
    alea_ray_coverage_worker_scratch_t* workers;
    size_t worker_count;
} alea_ray_coverage_executor_t;

void alea_ray_coverage_executor_init(alea_ray_coverage_executor_t* executor);
void alea_ray_coverage_executor_free(alea_ray_coverage_executor_t* executor);
int alea_ray_coverage_executor_prepare(alea_ray_coverage_executor_t* executor,
                                        size_t worker_count);
/* Stable row ownership for serial and future OpenMP schedulers. */
alea_ray_coverage_worker_scratch_t*
alea_ray_coverage_executor_worker_for_row(
    alea_ray_coverage_executor_t* executor, size_t row_index);

/* Execute independent coverage rows into worker-local reusable arenas, then
 * compact them in input row order.  Each row is assigned to
 * row_index % worker_count in both serial and OpenMP builds.  Failure leaves
 * result's previous publication intact. */
int alea_ray_coverage_slice_build_executor_nocache(
    alea_system_t* sys, const alea_ray_coverage_row_t* rows, size_t row_count,
    const alea_ray_coverage_slice_limits_t* limits,
    alea_ray_coverage_executor_t* executor,
    alea_ray_coverage_slice_result_t* result);

/** Reusable ordered boundary-event storage for internal query consumers. */
typedef struct {
    alea_ray_boundary_event_vec_t events;
} alea_ray_boundary_event_result_t;

typedef struct {
    /* Preserve every reportable physical surface in a coincident crossing.
     * The default reports the lowest positive surface ID as the deterministic
     * canonical representative. Synthetic lattice events are always emitted. */
    bool include_all_coincident_physical;
    /* Zero means unlimited. These are enforced while materializing events. */
    uint64_t max_events;
    uint64_t max_output_bytes;
} alea_ray_boundary_event_options_internal_t;

void alea_ray_boundary_event_result_init(alea_ray_boundary_event_result_t* result);
void alea_ray_boundary_event_result_clear(alea_ray_boundary_event_result_t* result);
void alea_ray_boundary_event_result_free(alea_ray_boundary_event_result_t* result);

/**
 * Trace canonical global hits/segments and derive ownership boundary events.
 * Physical events are emitted only when adjacent resolved ownership changes;
 * synthetic lattice crossings are retained even when ownership is unchanged.
 */
int alea_raycast_boundary_events_reuse_nocache(
    alea_system_t* sys, const alea_ray_t* ray, double t_max,
    alea_raycast_result_t* trace,
    alea_ray_boundary_event_result_t* events);

int alea_raycast_boundary_events_with_options(
    alea_system_t* sys, const alea_ray_t* ray, double t_max,
    const alea_ray_boundary_event_options_internal_t* options,
    alea_raycast_result_t* trace,
    alea_ray_boundary_event_result_t* events);

/**
 * Execute an internal semantic ray query using reusable caller-owned storage.
 * Caches must already be prepared.  The current implementation deliberately
 * shares the canonical global trace; query kind controls materialization and
 * answer selection, so consumers can migrate before traversal early-stop
 * specializations land.  SEGMENTS clips its reusable segment vector to
 * [t_min, t_max]; a leading clipped interval has no reportable enter surface.
 * On failure every supplied output is cleared.
 */
int alea_raycast_query_reuse_nocache(
    alea_system_t* sys, const alea_ray_t* ray,
    const alea_ray_query_t* query, alea_raycast_result_t* trace,
    alea_ray_boundary_event_result_t* events,
    alea_ray_query_output_t* output);

/* ============================================================================
 * QUERY HELPERS
 * ============================================================================ */

/**
 * @brief Get point along ray at distance t
 */
static inline void alea_ray_point_at(const alea_ray_t* ray, double t,
                                    double* x, double* y, double* z) {
    *x = ray->ox + t * ray->dx;
    *y = ray->oy + t * ray->dy;
    *z = ray->oz + t * ray->dz;
}

/**
 * @brief Find first cell hit by ray
 *
 * @param sys CSG system
 * @param ox, oy, oz Ray origin
 * @param dx, dy, dz Ray direction
 * @param t_max Maximum distance
 * @param out_t Output: distance to first hit (if found)
 * @return Cell ID of first cell hit, or -1 if none
 */
int alea_ray_first_cell(alea_system_t* sys,
                       double ox, double oy, double oz,
                       double dx, double dy, double dz,
                       double t_max,
                       double* out_t);

/**
 * @brief Quick occlusion test — is any non-void cell along the ray?
 *
 * Like alea_ray_first_cell but returns 1/0 and early-exits without
 * building full segment arrays. Uses thread-local buffers.
 *
 * @return 1 if occluded, 0 if clear
 */
int alea_ray_is_occluded(alea_system_t* sys,
                        double ox, double oy, double oz,
                        double dx, double dy, double dz,
                        double t_max);

/**
 * @brief Compute total path length through cells with given material
 *
 * @param result Raycast result with segments
 * @param material_id Material to sum (0 for void, -1 for all)
 * @return Total path length
 */
double alea_raycast_path_length(const alea_raycast_result_t* result,
                               int material_id);

/**
 * @brief Cell-aware raycast using per-cell surface index
 *
 * Semantic equivalent of alea_raycast(). For non-lattice models this can
 * track through cells one at a time using only surfaces belonging to each
 * cell. For lattice models it delegates to the canonical DDA-aware raycast
 * pipeline so lattice element-boundary hits are not skipped.
 *
 * Requires alea_build_cell_surface_index() to have been called first.
 *
 * @param sys CSG system (must have cell surface index built)
 * @param ox, oy, oz Ray origin
 * @param dx, dy, dz Ray direction
 * @param t_max Maximum distance (0 for infinite)
 * @param result Output: cell segments
 * @return 0 on success, -1 on error
 */
int alea_raycast_cell_aware(alea_system_t* sys,
                           double ox, double oy, double oz,
                           double dx, double dy, double dz,
                           double t_max,
                           alea_raycast_result_t* result);


#endif /* ALEA_RAYCAST_INTERNAL_H */
