// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/**
 * @file alea_raycast.h
 * @brief Alea Raycast Module API
 *
 * Optional raycast functionality. Requires linking with libalea_raycast.a
 */

#ifndef ALEA_RAYCAST_H
#define ALEA_RAYCAST_H

#include "alea.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ============================================================================
 * RAYCAST TYPES
 * ============================================================================ */

/** Opaque raycast result type */
typedef struct alea_raycast_result alea_raycast_result_t;

/** Opaque compact result for a batch of hierarchical ray traces. */
typedef struct alea_raycast_batch_result alea_raycast_batch_result_t;

/** Opaque reusable answer for a first-visible ray query. */
typedef struct alea_ray_first_visible_query_result
    alea_ray_first_visible_query_result_t;
typedef struct alea_ray_boundary_event_query_result
    alea_ray_boundary_event_query_result_t;
typedef struct alea_ray_coverage_slice_result
    alea_ray_coverage_slice_result_t;

#define ALEA_RAY_FIRST_VISIBLE_SURFACE_ID      (1u << 0)
#define ALEA_RAY_FIRST_VISIBLE_SURFACE_NORMAL  (1u << 1)

/** Public first-visible query options.  Set struct_size to sizeof(*options).
 * A zero t_max is unbounded; material_filter < 0 accepts every non-void
 * material.  Older struct prefixes are accepted safely. */
typedef struct {
    size_t struct_size;
    uint32_t fields;
    double t_min;
    double t_max;
    int material_filter;
} alea_ray_first_visible_options_t;

#define ALEA_RAY_EVENT_PHYSICAL          0
#define ALEA_RAY_EVENT_SYNTHETIC_LATTICE 1
#define ALEA_RAY_EVENT_UNRESOLVED        2
#define ALEA_RAY_BOUNDARY_EVENT_PRIMITIVE_ID (1u << 0)
#define ALEA_RAY_BOUNDARY_EVENT_NORMAL       (1u << 1)

typedef struct {
    size_t struct_size;
    uint32_t fields;
    double t_min;
    double t_max;
    uint64_t max_events;
    uint64_t max_output_bytes;
    int include_all_coincident_physical;
} alea_ray_boundary_event_options_t;

/* Optional fields in alea_raycast_batch_result_t. Distances and cell IDs are
 * always present; an accessor for an unrequested optional field returns NULL. */
#define ALEA_RAY_BATCH_MATERIAL          (1u << 0)
#define ALEA_RAY_BATCH_DENSITY           (1u << 1)
#define ALEA_RAY_BATCH_SURFACES          (1u << 2)
#define ALEA_RAY_BATCH_RESOLUTION_FLAGS  (1u << 3)
#define ALEA_RAY_BATCH_PROJECTED_OWNER   (1u << 4)
#define ALEA_RAY_BATCH_FULL_PATHS         (1u << 5)

/** Options for alea_raycast_hier_batch(). Set struct_size to sizeof(*options).
 * max_segments == 0 permits an unbounded result; nonzero rejects batches whose
 * total segment count exceeds the limit. */
typedef struct {
    size_t struct_size;
    uint32_t fields;
    int projected_depth;      /**< -1 = leaf; >= 0 = requested path depth */
    uint64_t max_segments;
    uint64_t max_path_entries; /**< Applies to flattened full-path CSR output */
    uint64_t max_output_bytes; /**< 0 = unbounded compact-array allocation */
} alea_raycast_batch_options_t;

/* Complete-ownership diagnostic coverage is published as input-order CSR:
 * rows -> intervals -> concrete owner occurrences.  Returned arrays are
 * borrowed and remain valid until the next successful query on the result or
 * until its destruction. */
typedef enum {
    ALEA_RAY_COVERAGE_UNIQUE,
    ALEA_RAY_COVERAGE_GAP,
    ALEA_RAY_COVERAGE_ALLOWED_EXTERIOR,
    ALEA_RAY_COVERAGE_OVERLAP,
    ALEA_RAY_COVERAGE_UNDEFINED_FILL,
    ALEA_RAY_COVERAGE_UNRESOLVED,
    ALEA_RAY_COVERAGE_TRUNCATED
} alea_ray_coverage_kind_t;

typedef enum {
    ALEA_RAY_COVERAGE_REFINEMENT_COMPLETE = 0,
    ALEA_RAY_COVERAGE_REFINEMENT_MAX_DEPTH,
    ALEA_RAY_COVERAGE_REFINEMENT_MAX_ROWS,
    ALEA_RAY_COVERAGE_REFINEMENT_MIN_SPACING
} alea_ray_coverage_refinement_status_t;

#define ALEA_RAY_COVERAGE_DOMAIN             (1u << 0)
#define ALEA_RAY_COVERAGE_REPORT_EXTERIOR    (1u << 1)
#define ALEA_RAY_COVERAGE_REFINE_SIGNATURE    (1u << 0)
#define ALEA_RAY_COVERAGE_REFINE_DISPLACEMENT (1u << 1)
#define ALEA_RAY_COVERAGE_REFINE_DENSITY      (1u << 2)
#define ALEA_RAY_COVERAGE_REFINE_FINDING      (1u << 3)

/** Options for compact complete-coverage rows. Set struct_size to
 * sizeof(*options). t_max must be finite and positive. When DOMAIN is set,
 * [domain_t_min, domain_t_max] is the ownership-validation domain; otherwise
 * unowned intervals are reported as gaps. REPORT_EXTERIOR requires DOMAIN.
 * Zero resource limits are unlimited. A nonzero max_refinement_depth enables
 * deterministic midpoint refinement; a zero refinement_signals field uses
 * signature-only selection.
 */
typedef struct {
    size_t struct_size;
    uint32_t flags;
    double t_max;
    double domain_t_min;
    double domain_t_max;
    uint64_t max_rows;
    uint64_t max_intervals;
    uint64_t max_owners;
    uint64_t max_output_bytes;
    size_t max_refinement_depth;
    uint32_t refinement_signals;
    double min_transverse_spacing;
    double endpoint_displacement;
    size_t crossing_density;
} alea_ray_coverage_slice_options_t;

/** One resolved hierarchy entry attached to an opt-in ray segment path.
 *
 * ``occurrence_key`` identifies this concrete occurrence, including its
 * ancestor placement/lattice context.  It is stable for equivalent paths
 * traced in opposite directions through an unchanged system.
 */
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
} alea_raycast_path_entry_t;

/* ============================================================================
 * RAYCAST FUNCTIONS
 * ============================================================================ */

/**
 * @brief Create a new raycast result object
 * @return New result object (caller must destroy with alea_raycast_result_destroy)
 */
alea_raycast_result_t* alea_raycast_result_create(void);

/**
 * @brief Cast a ray through geometry and find all cell intersections
 *
 * @param sys System
 * @param ox, oy, oz Ray origin
 * @param dx, dy, dz Ray direction (will be normalized)
 * @param t_max Maximum ray distance (0 = infinite)
 * @param result Output result
 * @return 0 on success, -1 on error
 */
int alea_raycast(alea_system_t* sys,
                     double ox, double oy, double oz,
                     double dx, double dy, double dz,
                     double t_max,
                     alea_raycast_result_t* result);

/**
 * @brief Cell-aware raycast using per-cell surface index
 *
 * Semantic equivalent of alea_raycast(). For non-lattice models this can
 * track through cells one at a time using the per-cell surface index. For
 * lattice models it uses the canonical lattice-aware pipeline so synthetic
 * DDA element-boundary hits are included.
 */
int alea_raycast_cell_aware(alea_system_t* sys,
                                double ox, double oy, double oz,
                                double dx, double dy, double dz,
                                double t_max,
                                alea_raycast_result_t* result);

/**
 * @brief Fast hierarchical raycast for material/path segments.
 *
 * Uses the hierarchical spatial index and per-ray path state to step through
 * cells without building a flat spatial index or collecting global surface
 * hits. The returned segments contain boundary surface IDs where available,
 * but the result's full hit list is not part of this function's contract.
 */
int alea_raycast_hier_fast_segments(alea_system_t* sys,
                                    double ox, double oy, double oz,
                                    double dx, double dy, double dz,
                                    double t_max,
                                    alea_raycast_result_t* result);

/* ============================================================================
 * COMPACT HIERARCHICAL BATCH RAYCAST
 * ============================================================================ */

/** Create/destroy a reusable compact batch result. */
alea_raycast_batch_result_t* alea_raycast_batch_result_create(void);
void alea_raycast_batch_result_destroy(alea_raycast_batch_result_t* result);

/**
 * Generate packed Cauchy-Crofton rays through a bounding sphere.
 *
 * Each ray has a uniform direction and a uniform disk offset perpendicular
 * to that direction.  Its origin is placed 2 * radius behind the disk, so a
 * t_max of 4 * radius spans the enclosing sphere.  The output buffers hold
 * ray_count XYZ triples and can be passed directly to
 * alea_raycast_hier_batch().  rng_state is advanced deterministically.
 *
 * @param cx,cy,cz Sphere center
 * @param radius Positive finite sphere radius
 * @param rng_state In/out LCG state (must not be NULL)
 * @param ray_count Number of rays to generate
 * @param origins_xyz Output packed XYZ origins; required when ray_count > 0
 * @param directions_xyz Output packed XYZ directions; required when ray_count > 0
 * @return 0 on success, -1 on invalid input
 */
int alea_generate_cauchy_crofton_rays(double cx, double cy, double cz,
                                      double radius, uint32_t* rng_state,
                                      size_t ray_count, double* origins_xyz,
                                      double* directions_xyz);

/**
 * Trace packed XYZ origin/direction pairs into compact CSR arrays.
 *
 * The output ray order matches input order. ray_offsets has ray_count + 1
 * entries; a ray's segments occupy [offsets[i], offsets[i + 1]). Directions
 * are normalized using the same rules as the single-ray API. The function
 * prepares required hierarchical caches before parallel tracing.
 */
int alea_raycast_hier_batch(alea_system_t* sys,
                            const double* origins_xyz,
                            const double* directions_xyz,
                            size_t ray_count,
                            double t_max,
                            const alea_raycast_batch_options_t* options,
                            alea_raycast_batch_result_t* result);

size_t alea_raycast_batch_ray_count(const alea_raycast_batch_result_t* result);
size_t alea_raycast_batch_segment_count(const alea_raycast_batch_result_t* result);
uint32_t alea_raycast_batch_fields(const alea_raycast_batch_result_t* result);

const uint64_t* alea_raycast_batch_ray_offsets(
    const alea_raycast_batch_result_t* result);
const double* alea_raycast_batch_t_enter(
    const alea_raycast_batch_result_t* result);
const double* alea_raycast_batch_t_exit(
    const alea_raycast_batch_result_t* result);
const int32_t* alea_raycast_batch_cell_ids(
    const alea_raycast_batch_result_t* result);
const int32_t* alea_raycast_batch_material_ids(
    const alea_raycast_batch_result_t* result);
const double* alea_raycast_batch_densities(
    const alea_raycast_batch_result_t* result);
const int32_t* alea_raycast_batch_enter_surface_ids(
    const alea_raycast_batch_result_t* result);
const int32_t* alea_raycast_batch_exit_surface_ids(
    const alea_raycast_batch_result_t* result);
const uint8_t* alea_raycast_batch_resolution_flags(
    const alea_raycast_batch_result_t* result);
const int32_t* alea_raycast_batch_projected_cell_ids(
    const alea_raycast_batch_result_t* result);
const int32_t* alea_raycast_batch_projected_material_ids(
    const alea_raycast_batch_result_t* result);
const int32_t* alea_raycast_batch_projected_universe_ids(
    const alea_raycast_batch_result_t* result);
const int32_t* alea_raycast_batch_projected_fill_universes(
    const alea_raycast_batch_result_t* result);
const int32_t* alea_raycast_batch_projected_depths(
    const alea_raycast_batch_result_t* result);
const uint8_t* alea_raycast_batch_projected_is_lattice(
    const alea_raycast_batch_result_t* result);
const uint64_t* alea_raycast_batch_projected_occurrence_keys(
    const alea_raycast_batch_result_t* result);

/** Full hierarchy paths use a second CSR layer when ALEA_RAY_BATCH_FULL_PATHS
 * is requested. segment_path_offsets has segment_count + 1 entries. */
size_t alea_raycast_batch_path_entry_count(
    const alea_raycast_batch_result_t* result);
const uint64_t* alea_raycast_batch_segment_path_offsets(
    const alea_raycast_batch_result_t* result);
const int32_t* alea_raycast_batch_path_cell_ids(
    const alea_raycast_batch_result_t* result);
const int32_t* alea_raycast_batch_path_material_ids(
    const alea_raycast_batch_result_t* result);
const int32_t* alea_raycast_batch_path_universe_ids(
    const alea_raycast_batch_result_t* result);
const int32_t* alea_raycast_batch_path_fill_universes(
    const alea_raycast_batch_result_t* result);
const int32_t* alea_raycast_batch_path_depths(
    const alea_raycast_batch_result_t* result);
const uint8_t* alea_raycast_batch_path_is_lattice(
    const alea_raycast_batch_result_t* result);
const double* alea_raycast_batch_path_lattice_origins_xyz(
    const alea_raycast_batch_result_t* result);
const uint64_t* alea_raycast_batch_path_occurrence_keys(
    const alea_raycast_batch_result_t* result);

/* ==========================================================================
 * COMPACT COMPLETE-COVERAGE DIAGNOSTICS
 * ========================================================================== */

void alea_ray_coverage_slice_options_init(
    alea_ray_coverage_slice_options_t* options);
alea_ray_coverage_slice_result_t* alea_ray_coverage_slice_result_create(void);
void alea_ray_coverage_slice_result_destroy(
    alea_ray_coverage_slice_result_t* result);

/** Scalar adapter over alea_ray_coverage_slice_query(). It publishes one row
 * of compact CSR data using the same ownership and resource-limit contract. */
int alea_ray_coverage_query(alea_system_t* sys,
    double ox, double oy, double oz, double dx, double dy, double dz,
    const alea_ray_coverage_slice_options_t* options,
    alea_ray_coverage_slice_result_t* result);

/** Build complete ownership coverage for packed XYZ rays. direction_tags and
 * transverse_coordinates are optional row provenance arrays; NULL supplies
 * zeros. Output row order always matches the input. On failure, result keeps
 * its prior successful publication. */
int alea_ray_coverage_slice_query(alea_system_t* sys,
    const double* origins_xyz, const double* directions_xyz, size_t row_count,
    const uint8_t* direction_tags, const double* transverse_coordinates,
    const alea_ray_coverage_slice_options_t* options,
    alea_ray_coverage_slice_result_t* result);

size_t alea_ray_coverage_slice_row_count(
    const alea_ray_coverage_slice_result_t* result);
size_t alea_ray_coverage_slice_interval_count(
    const alea_ray_coverage_slice_result_t* result);
size_t alea_ray_coverage_slice_owner_count(
    const alea_ray_coverage_slice_result_t* result);
int alea_ray_coverage_slice_refinement_status(
    const alea_ray_coverage_slice_result_t* result);
const size_t* alea_ray_coverage_slice_row_offsets(
    const alea_ray_coverage_slice_result_t* result);
const uint8_t* alea_ray_coverage_slice_row_direction_tags(
    const alea_ray_coverage_slice_result_t* result);
const double* alea_ray_coverage_slice_row_transverse_coordinates(
    const alea_ray_coverage_slice_result_t* result);
const double* alea_ray_coverage_slice_t_enter(
    const alea_ray_coverage_slice_result_t* result);
const double* alea_ray_coverage_slice_t_exit(
    const alea_ray_coverage_slice_result_t* result);
const uint8_t* alea_ray_coverage_slice_kinds(
    const alea_ray_coverage_slice_result_t* result);
const size_t* alea_ray_coverage_slice_owner_offsets(
    const alea_ray_coverage_slice_result_t* result);
const size_t* alea_ray_coverage_slice_owner_count_lower_bounds(
    const alea_ray_coverage_slice_result_t* result);
const int* alea_ray_coverage_slice_owner_cell_ids(
    const alea_ray_coverage_slice_result_t* result);
const int* alea_ray_coverage_slice_owner_material_ids(
    const alea_ray_coverage_slice_result_t* result);
const int* alea_ray_coverage_slice_owner_universe_ids(
    const alea_ray_coverage_slice_result_t* result);
const int* alea_ray_coverage_slice_owner_fill_universes(
    const alea_ray_coverage_slice_result_t* result);
const int* alea_ray_coverage_slice_owner_depths(
    const alea_ray_coverage_slice_result_t* result);
const uint64_t* alea_ray_coverage_slice_owner_occurrence_keys(
    const alea_ray_coverage_slice_result_t* result);
const uint64_t* alea_ray_coverage_slice_owner_parent_occurrence_keys(
    const alea_ray_coverage_slice_result_t* result);
const uint8_t* alea_ray_coverage_slice_owner_resolution_flags(
    const alea_ray_coverage_slice_result_t* result);

/**
 * @brief Find first cell along ray
 *
 * @param sys System
 * @param ox, oy, oz Ray origin
 * @param dx, dy, dz Ray direction
 * @param t_max Maximum distance (0 = infinite)
 * @param out_t Output: distance to first hit (can be NULL)
 * @return Cell ID of first cell hit, or -1 if none
 */
int alea_ray_first_cell(alea_system_t* sys,
                            double ox, double oy, double oz,
                            double dx, double dy, double dz,
                            double t_max, double* out_t);

/* ==========================================================================
 * FIRST-VISIBLE QUERY
 * ========================================================================== */

void alea_ray_first_visible_options_init(
    alea_ray_first_visible_options_t* options);
alea_ray_first_visible_query_result_t* alea_ray_first_visible_query_result_create(void);
void alea_ray_first_visible_query_result_destroy(
    alea_ray_first_visible_query_result_t* result);

/** Execute a reusable first-visible query.  On failure, `result` is cleared.
 * Returned pointers remain valid until the next query on `result` or its
 * destruction. */
int alea_ray_first_visible_query(alea_system_t* sys,
    double ox, double oy, double oz, double dx, double dy, double dz,
    const alea_ray_first_visible_options_t* options,
    alea_ray_first_visible_query_result_t* result);

int alea_ray_first_visible_found(
    const alea_ray_first_visible_query_result_t* result);
double alea_ray_first_visible_t(
    const alea_ray_first_visible_query_result_t* result);
int alea_ray_first_visible_cell_id(
    const alea_ray_first_visible_query_result_t* result);
int alea_ray_first_visible_material_id(
    const alea_ray_first_visible_query_result_t* result);
double alea_ray_first_visible_density(
    const alea_ray_first_visible_query_result_t* result);
int alea_ray_first_visible_surface_id(
    const alea_ray_first_visible_query_result_t* result);
int alea_ray_first_visible_normal(
    const alea_ray_first_visible_query_result_t* result,
    double* nx, double* ny, double* nz);

/* ==========================================================================
 * BOUNDARY-EVENT QUERY
 * ========================================================================== */

void alea_ray_boundary_event_options_init(
    alea_ray_boundary_event_options_t* options);
alea_ray_boundary_event_query_result_t* alea_ray_boundary_event_query_result_create(void);
void alea_ray_boundary_event_query_result_destroy(
    alea_ray_boundary_event_query_result_t* result);
int alea_ray_boundary_event_query(alea_system_t* sys,
    double ox, double oy, double oz, double dx, double dy, double dz,
    const alea_ray_boundary_event_options_t* options,
    alea_ray_boundary_event_query_result_t* result);
size_t alea_ray_boundary_event_count(
    const alea_ray_boundary_event_query_result_t* result);
int alea_ray_boundary_event_get(
    const alea_ray_boundary_event_query_result_t* result, size_t index,
    double* t, int* kind, int* surface_id, int* cell_before, int* cell_after,
    int* material_before, int* material_after, uint32_t* resolution_flags,
    uint32_t* primitive_id, double* nx, double* ny, double* nz);

/**
 * @brief Get number of segments in raycast result
 */
size_t alea_raycast_segment_count(const alea_raycast_result_t* result);

/**
 * @brief Get segment data
 * @param result Raycast result
 * @param index Segment index
 * @param t_enter Output: entry distance
 * @param t_exit Output: exit distance
 * @param cell_id Output: cell ID (-1 for void)
 * @param material_id Output: material ID
 * @param density Output: material density
 * @param enter_surface_id Output: boundary crossed at t_enter (-1 none, 0 synthetic, >0 physical)
 * @param exit_surface_id Output: boundary crossed at t_exit (-1 none, 0 synthetic, >0 physical)
 * @return 0 on success, -1 on invalid index
 */
int alea_raycast_segment_get(const alea_raycast_result_t* result, size_t index,
                                 double* t_enter, double* t_exit,
                                 int* cell_id, int* material_id, double* density,
                                 int* enter_surface_id, int* exit_surface_id);

/** Get per-segment resolution flags. Returns 0 on success, -1 on error. */
int alea_raycast_segment_resolution_flags(const alea_raycast_result_t* result,
                                          size_t index,
                                          uint8_t* out_flags);

/** Return the number of physical boundary hits retained by a raycast result. */
size_t alea_raycast_hit_count(const alea_raycast_result_t* result);

/**
 * @brief Get one retained physical boundary hit.
 *
 * @return 0 on success, -1 for a null result, invalid index, or null output.
 */
int alea_raycast_hit_get(const alea_raycast_result_t* result,
                         size_t index,
                         double* out_t,
                         int* out_surface_id);

/** Enable or disable hierarchy-path capture for subsequent hierarchical traces.
 * Disabled by default.  The setting survives result-buffer reuse. */
void alea_raycast_result_set_path_capture(alea_raycast_result_t* result,
                                          int enabled);

/** Return the number of hierarchy entries attached to a segment. */
size_t alea_raycast_segment_path_count(const alea_raycast_result_t* result,
                                       size_t segment_index);

/** Copy one hierarchy path entry attached to a segment. */
int alea_raycast_segment_path_get(const alea_raycast_result_t* result,
                                  size_t segment_index,
                                  size_t path_entry_index,
                                  alea_raycast_path_entry_t* out_entry);

/* ============================================================================
 * INTERVAL DEFECT CLASSIFICATION (geometry error detection, oracle-grade)
 * ============================================================================ */

/** Kinds for alea_ray_interval_finding_t. */
#define ALEA_INTERVAL_OK             0  /**< exactly one containing chain */
#define ALEA_INTERVAL_GAP            1  /**< no cell contains the interval */
#define ALEA_INTERVAL_OVERLAP        2  /**< claimants cannot form one ownership chain */
#define ALEA_INTERVAL_UNDEFINED_FILL 3  /**< container with no fill content */
#define ALEA_INTERVAL_UNRESOLVED     4  /**< coverage ancestry/numerics indeterminate */
#define ALEA_INTERVAL_TRUNCATED      5  /**< owner budget prevented complete coverage */

typedef struct {
    double t_enter;
    double t_exit;
    int kind;            /**< ALEA_INTERVAL_* */
    int cell_id;         /**< display owner (first containing chain); -1 gap */
    int overlap_cell_id; /**< second claimant at the overlap depth, else -1 */
    int depth;           /**< depth of the anomaly (overlap) or of the owner */
} alea_ray_interval_finding_t;

/**
 * @brief Classify elementary ray intervals by cell ownership.
 *
 * Collects ALL surface crossings along the ray (global pipeline, including
 * fill-transformed and synthetic lattice boundaries), then determines the
 * complete containing-cell set of each interval with an uncached recursive
 * query at an interior point. Within such an interval the owner set is
 * invariant, so overlaps and gaps are detected regardless of their size —
 * including overlaps in isolation (a cell fully inside another) that
 * produce no ownership transition and are invisible to trace segments.
 *
 * Consecutive intervals with identical owner sets are merged. This is the
 * completeness-first ("oracle") detector: cost is one recursive query per
 * elementary interval; use trace()/segments for rendering speed.
 *
 * @param sys System
 * @param ox,oy,oz Ray origin
 * @param dx,dy,dz Ray direction (normalized internally)
 * @param t_max Maximum distance (must be > 0)
 * @param out Output findings array (may be NULL if max_out is 0)
 * @param max_out Capacity of out
 * @return Total number of merged intervals (may exceed max_out, in which
 *         case only the first max_out are written), or -1 on error
 */
int alea_ray_classify_intervals(alea_system_t* sys,
                                double ox, double oy, double oz,
                                double dx, double dy, double dz,
                                double t_max,
                                alea_ray_interval_finding_t* out,
                                size_t max_out);

/**
 * @brief Calculate total path length through material
 *
 * @param result Raycast result
 * @param material_id Material to sum (-1 = all materials)
 * @return Total path length
 */
double alea_raycast_path_length(const alea_raycast_result_t* result, int material_id);

/**
 * @brief Free raycast result internal buffers (for stack-allocated results)
 */
void alea_raycast_result_free(alea_raycast_result_t* result);

/**
 * @brief Destroy a heap-allocated raycast result (from alea_raycast_result_create)
 */
void alea_raycast_result_destroy(alea_raycast_result_t* result);

#ifdef __cplusplus
}
#endif

#endif /* ALEA_RAYCAST_H */
