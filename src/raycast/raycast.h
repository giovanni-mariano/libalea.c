// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

#ifndef ALEA_RAYCAST_INTERNAL_H
#define ALEA_RAYCAST_INTERNAL_H

#include "alea_types.h"
#include "alea_raycast.h"
#include "alea_slice.h"
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
void alea_raycast_batch_result_swap_internal(
    alea_raycast_batch_result_t* a,
    alea_raycast_batch_result_t* b);

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

/** Reusable ordered boundary-event storage for internal query consumers. */
typedef struct {
    alea_ray_boundary_event_vec_t events;
} alea_ray_boundary_event_result_t;

typedef struct {
    /* Preserve every reportable physical surface in a coincident crossing.
     * The default reports the lowest positive surface ID as the deterministic
     * canonical representative. Synthetic lattice events are always emitted. */
    bool include_all_coincident_physical;
} alea_ray_boundary_event_options_t;

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
    const alea_ray_boundary_event_options_t* options,
    alea_raycast_result_t* trace,
    alea_ray_boundary_event_result_t* events);

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
