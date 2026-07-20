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
#define ALEA_INTERVAL_OVERLAP        2  /**< >=2 cells at the same depth */
#define ALEA_INTERVAL_UNDEFINED_FILL 3  /**< container with no fill content */

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
