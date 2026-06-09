// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

#ifndef ALEA_RAYCAST_INTERNAL_H
#define ALEA_RAYCAST_INTERNAL_H

#include "alea_types.h"
#include "util/alea_vec.h"
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
} alea_ray_segment_t;

ALEA_VEC_DEFINE(alea_ray_hit_vec, alea_ray_hit_t);
ALEA_VEC_DEFINE(alea_ray_segment_vec, alea_ray_segment_t);

/* Use struct tag matching the public API forward declaration */
struct alea_raycast_result {
    alea_ray_t ray;
    alea_ray_hit_vec_t hits;
    alea_ray_segment_vec_t segments;
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
};

/* Typedef for internal use */
typedef struct alea_raycast_result alea_raycast_result_t;

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
