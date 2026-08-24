// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

#ifndef ALEA_SPATIAL_HIER_H
#define ALEA_SPATIAL_HIER_H

#include <stddef.h>
#include <stdint.h>
#include "alea_types.h"
#include "alea_universe.h"  /* alea_matrix_t */

typedef struct alea_system alea_system_t;
typedef struct alea_hier_spatial_index alea_hier_spatial_index_t;

/**
 * @brief Query result for region/slice/point queries.
 *
 * `instance_index` is a synthetic per-query index in hierarchical mode.
 */
typedef struct {
    uint32_t instance_index;
    uint32_t cell_index;       /* Index into sys->cells */
    int cell_id;               /* MCNP cell ID */
    int material_id;           /* Material ID */
    int universe_id;           /* Universe this instance belongs to */
    int depth;                 /* Nesting depth (0 = base universe) */
    bool is_terminal;          /* True if no FILL (actual geometry) */
    alea_matrix_t transform;    /* Transform to apply */
} alea_spatial_hit_t;

typedef struct {
    alea_cell_hit_t hit;
    alea_matrix_t transform;
    int lattice_cell_index;
    alea_matrix_t lattice_transform;
} alea_hier_cell_hit_t;

#define ALEA_HIER_SPATIAL_HIT_CHAIN_MAX 16
#define ALEA_HIER_RAY_PATH_MAX 64

typedef struct {
    alea_spatial_hit_t hit;
    uint8_t ancestor_count;
    uint8_t chain_truncated;
    uint32_t ancestor_cell_indices[ALEA_HIER_SPATIAL_HIT_CHAIN_MAX];
    alea_matrix_t ancestor_transforms[ALEA_HIER_SPATIAL_HIT_CHAIN_MAX];
    /* A lattice ancestor constrains a terminal occurrence to one concrete
     * element.  Its transform alone is insufficient because identical child
     * universes may be instantiated in adjacent elements. */
    uint8_t ancestor_is_lattice[ALEA_HIER_SPATIAL_HIT_CHAIN_MAX];
    int ancestor_lattice_fill_universes[ALEA_HIER_SPATIAL_HIT_CHAIN_MAX];
    int ancestor_lattice_i[ALEA_HIER_SPATIAL_HIT_CHAIN_MAX];
    int ancestor_lattice_j[ALEA_HIER_SPATIAL_HIT_CHAIN_MAX];
    int ancestor_lattice_k[ALEA_HIER_SPATIAL_HIT_CHAIN_MAX];
    double ancestor_lattice_ox[ALEA_HIER_SPATIAL_HIT_CHAIN_MAX];
    double ancestor_lattice_oy[ALEA_HIER_SPATIAL_HIT_CHAIN_MAX];
    double ancestor_lattice_oz[ALEA_HIER_SPATIAL_HIT_CHAIN_MAX];
} alea_hier_spatial_chain_hit_t;

typedef struct {
    uint32_t cell_index;
    int cell_id;
    int material_id;
    int universe_id;
    int fill_universe;
    int depth;
    uint8_t is_lattice;
    int lat_fill_universe;
    int lat_i;
    int lat_j;
    int lat_k;
    double lat_ox;
    double lat_oy;
    double lat_oz;
    alea_matrix_t transform;
} alea_hier_ray_path_entry_t;

typedef struct {
    int count;
    alea_hier_ray_path_entry_t entries[ALEA_HIER_RAY_PATH_MAX];
} alea_hier_ray_path_t;

/* Worker-local state for coherent point/path queries.  It is deliberately
 * bounded by one path per active worker; callers swap two workspaces rather
 * than retaining a path for every grid column. */
typedef struct {
    alea_hier_ray_path_t path;
    alea_hier_cell_hit_t deepest;
    uint64_t system_id;
    uint64_t geometry_generation;
    uint8_t complete;
    uint8_t overflowed;
} alea_hier_coherence_state_t;

typedef enum {
    ALEA_HIER_COH_ROOT_QUERY,
    ALEA_HIER_COH_PATH_REUSED,
    ALEA_HIER_COH_LATTICE_TRANSITION,
    ALEA_HIER_COH_PREFIX_RESTART,
    ALEA_HIER_COH_FULL_FALLBACK
} alea_hier_coherence_kind_t;

typedef struct {
    uint32_t cell_index;
    int cell_id;
    double t_enter;
    double t_exit;
} alea_hier_ray_candidate_t;

typedef struct {
    uint32_t placement_index;
    uint32_t parent_cell_index;
    int universe_id;
    int depth;
    uint32_t flags;
    alea_matrix_t transform;
    double t_enter;
    double t_exit;
} alea_hier_placement_ray_candidate_t;

typedef struct {
    size_t universe_count;
    size_t blas_count;
    size_t linear_universe_count;
    size_t blas_cell_count;
    size_t blas_node_count;
    size_t fill_cell_count;
    size_t lattice_cell_count;
    size_t transform_count;
    size_t placement_count;
    size_t root_placement_count;
    size_t fill_placement_count;
    size_t lattice_placement_count;
    int max_placement_depth;
    int max_universe_cells;
    int largest_universe_id;
    size_t memory_bytes;
} alea_hier_spatial_stats_t;

int alea_hier_spatial_index_build(alea_system_t* sys);
void alea_hier_spatial_index_free(alea_hier_spatial_index_t* idx);
const alea_hier_spatial_stats_t*
alea_hier_spatial_index_stats(const alea_hier_spatial_index_t* idx);
int alea_hier_spatial_find_cells_at_point(alea_system_t* sys,
                                          double x,
                                          double y,
                                          double z,
                                          alea_cell_hit_t* out_hits,
                                          size_t max_hits);
int alea_hier_spatial_find_cells_at_point_uncached(alea_system_t* sys,
                                                   double x,
                                                   double y,
                                                   double z,
                                                   alea_cell_hit_t* out_hits,
                                                   size_t max_hits);
int alea_hier_spatial_find_deepest_cell_at_point(alea_system_t* sys,
                                                 double x,
                                                 double y,
                                                 double z,
                                                 alea_hier_cell_hit_t* out_hit);
int alea_hier_spatial_find_path_at_point(alea_system_t* sys,
                                         double x,
                                         double y,
                                         double z,
                                         alea_hier_cell_hit_t* out_hit,
                                         alea_hier_ray_path_t* out_path);
void alea_hier_coherence_state_clear(alea_hier_coherence_state_t* state);

/* How a reused path resolves ownership.  The two modes differ only where more
 * than one cell of a universe contains the point, which is illegal geometry:
 * elsewhere the owner is unique and both agree. */
typedef enum {
    /* Accept the cached path once it still contains the point.  Ownership in
     * an overlap then follows the scan that produced the path, the way a
     * tracked particle keeps the cell it entered.  Cheap: no owner lookup. */
    ALEA_HIER_COH_OWNERSHIP_COHERENT = 0,
    /* Re-derive the deck-first owner at every level, so the result matches a
     * from-scratch query point for point regardless of how the sweep was
     * split.  Costs a BVH descent per level per point. */
    ALEA_HIER_COH_OWNERSHIP_CANONICAL = 1
} alea_hier_coherence_ownership_t;

/* Resolve a point from a worker-local prior path when its validated prefix
 * still applies.  Callers that report overlaps separately (the grid's boundary
 * pass, the geometry validator) want ALEA_HIER_COH_OWNERSHIP_COHERENT; ask for
 * CANONICAL only when the result must be independent of the sweep order. */
int alea_hier_spatial_resolve_coherent(
    alea_system_t* sys,
    double x, double y, double z,
    const alea_hier_coherence_state_t* previous,
    alea_hier_coherence_ownership_t ownership,
    alea_hier_coherence_state_t* current,
    alea_hier_cell_hit_t* out_hit,
    alea_hier_coherence_kind_t* out_kind);
/* Validate the ancestor prefix through target_entry.  The entry index is
 * required because one cell definition may occur multiple times in a path
 * through distinct placements. */
int alea_hier_spatial_check_path_containment(alea_system_t* sys,
                                             const alea_hier_ray_path_t* path,
                                             int target_entry,
                                             double x,
                                             double y,
                                             double z,
                                             alea_matrix_t* out_transform,
                                             int* out_lattice_cell_index,
                                             alea_matrix_t* out_lattice_transform);
int alea_hier_spatial_find_path_from_parent(alea_system_t* sys,
                                            const alea_hier_ray_path_t* path,
                                            int parent_entry,
                                            double x,
                                            double y,
                                            double z,
                                            alea_hier_cell_hit_t* out_hit,
                                            alea_hier_ray_path_t* out_path);
/* Restart below an already validated ordinary fill without re-deriving its
 * deck-order owner.  This is the particle-tracking form of the operation:
 * in illegal overlaps, it retains the path's current owner while that owner
 * still contains the point.  Keep the canonical function above for callers
 * that deliberately need point-query/deck-order semantics. */
int alea_hier_spatial_find_path_from_parent_coherent(
    alea_system_t* sys,
    const alea_hier_ray_path_t* path,
    int parent_entry,
    double x,
    double y,
    double z,
    alea_hier_cell_hit_t* out_hit,
    alea_hier_ray_path_t* out_path);
/* Replace a lattice entry's selected element, retain its validated ancestors,
 * and descend canonically from that concrete child placement. `location` must
 * be the canonical element containing (x,y,z) at `lattice_entry`. */
int alea_hier_path_enter_lattice_location(
    alea_system_t* sys,
    const alea_hier_ray_path_t* path,
    int lattice_entry,
    double x, double y, double z,
    const alea_lattice_location_t* location,
    alea_hier_cell_hit_t* out_hit,
    alea_hier_ray_path_t* out_path);

/**
 * @brief Reset the hier point-query coherence cache.
 *
 * Mirror of alea_spatial_reset_cache for the hierarchical query path.
 * Called automatically by alea_spatial_reset_cache so existing callers
 * that reset between renders do not need to learn a second API.
 */
void alea_hier_spatial_reset_cache(void);

/**
 * @brief Test whether a world-space point lies inside a cell, using the
 * transform cached by the most recent alea_hier_spatial_find_cells_at_point.
 *
 * Lets raycast's Tier 2 coherence check work for cells with universe_id
 * != 0, where the caller has no other access to the world→local transform.
 *
 * @return 1 if the point is inside the cell,
 *         0 if it is outside,
 *        -1 if the cell is not present in the current cache or the cache
 *           has been invalidated. Caller must fall back to a full lookup.
 */
int alea_hier_spatial_check_cached_containment(alea_system_t* sys,
                                               uint32_t cell_index,
                                               double x, double y, double z);

int alea_hier_spatial_get_cached_cell_state(alea_system_t* sys,
                                            uint32_t cell_index,
                                            double x, double y, double z,
                                            alea_matrix_t* out_transform,
                                            int* out_lattice_cell_index,
                                            alea_matrix_t* out_lattice_transform);

/**
 * @brief Find the cell of `universe_id` containing (lx, ly, lz) in that
 *        universe's local frame, using the per-universe BLAS for pruning.
 *
 * Does NOT descend into fills; returns the cell directly owning the point.
 * Returns the cell index on success, -1 if none, -2 on error (no index).
 */
int alea_hier_spatial_find_cell_in_universe(alea_system_t* sys,
                                            int universe_id,
                                            double lx,
                                            double ly,
                                            double lz);
/* crossed_mc_surface_id: if > 0, the expensive containment test runs first only
 * on candidates that reference that surface (the cell just across a crossed
 * boundary shares it), falling back to the full scan if none contains the point.
 * Pass -1 to disable the filter (exhaustive scan, e.g. for geometry validation). */
int alea_hier_spatial_find_ordered_cell_in_universe(alea_system_t* sys,
                                                    int universe_id,
                                                    double lx,
                                                    double ly,
                                                    double lz,
                                                    int crossed_mc_surface_id);
int alea_hier_spatial_query_universe_region(alea_system_t* sys,
                                            int universe_id,
                                            const alea_bbox_t* local_bbox,
                                            alea_spatial_hit_t* out_hits,
                                            size_t max_hits);
int alea_hier_spatial_query_universe_ray(alea_system_t* sys,
                                         int universe_id,
                                         double ox,
                                         double oy,
                                         double oz,
                                         double dx,
                                         double dy,
                                         double dz,
                                         double inv_dx,
                                         double inv_dy,
                                         double inv_dz,
                                         double t_min,
                                         double t_max,
                                         alea_hier_ray_candidate_t* out_hits,
                                         size_t max_hits);
int alea_hier_spatial_query_placements_ray(alea_system_t* sys,
                                           double ox,
                                           double oy,
                                           double oz,
                                           double dx,
                                           double dy,
                                           double dz,
                                           double inv_dx,
                                           double inv_dy,
                                           double inv_dz,
                                           double t_min,
                                           double t_max,
                                           alea_hier_placement_ray_candidate_t* out_hits,
                                           size_t max_hits);
typedef int (*alea_hier_lattice_placement_ray_visitor_t)(
    uint32_t placement_index,
    uint32_t lattice_cell_index,
    const alea_matrix_t* transform,
    double t_enter,
    double t_exit,
    void* userdata);
/* Optional traversal attribution for a single lattice-placement ray query.
 * Nodes include every TLAS node whose bounds are tested; leaves include only
 * leaves whose bounds overlap the ray interval. */
typedef struct {
    uint64_t nodes_tested;
    uint64_t leaves_visited;
} alea_hier_lattice_placement_ray_stats_t;
/* Visit the exact enclosing fill cells that constrain a lattice occurrence.
 * `transform` maps the visited cell's local frame to world space.  The
 * lattice fundamental cell is intentionally excluded: it describes one
 * repeated element rather than the occurrence support. */
typedef int (*alea_hier_lattice_ancestor_visitor_t)(
    uint32_t cell_index,
    const alea_matrix_t* transform,
    void* userdata);
int alea_hier_spatial_visit_lattice_placements_ray(
    alea_system_t* sys,
    double ox, double oy, double oz,
    double inv_dx, double inv_dy, double inv_dz,
    double t_min, double t_max,
    alea_hier_lattice_placement_ray_visitor_t visitor,
    alea_hier_lattice_placement_ray_stats_t* stats,
    void* userdata);
int alea_hier_spatial_visit_lattice_placement_ancestors(
    alea_system_t* sys,
    uint32_t placement_index,
    alea_hier_lattice_ancestor_visitor_t visitor,
    void* userdata);
int alea_hier_spatial_check_placement_chain(alea_system_t* sys,
                                            uint32_t placement_index,
                                            double x,
                                            double y,
                                            double z);
/* Validate the enclosing fill-placement chain for a lattice placement.
 * The lattice element itself is intentionally omitted: its finite extent and
 * occupancy are established by lattice DDA, whereas the lattice cell's CSG
 * describes one repeated template element rather than the full array. */
int alea_hier_spatial_check_lattice_placement_ancestors(
    alea_system_t* sys,
    uint32_t placement_index,
    double x,
    double y,
    double z);
/* As above, but also requires every enclosing cell to be the canonical
 * deck-order owner in its universe.  Use this before publishing a synthetic
 * lattice-entry event in potentially overlapping input. */
int alea_hier_spatial_check_lattice_placement_canonical_ancestors(
    alea_system_t* sys,
    uint32_t placement_index,
    double x,
    double y,
    double z);

/**
 * @brief Return the precomputed per-cell fill transform (forward + inverse)
 *        for `cell_index`, or NULL if the hier index is not built or this
 *        cell has no fill.
 *
 * Cached at hier-index build time so callers in hot loops can skip the
 * per-pixel `alea_matrix_from_mcnp` + `alea_matrix_invert` work. The returned
 * pointer is owned by the hier index and valid until it is destroyed.
 */
const alea_matrix_t*
alea_hier_spatial_get_cell_fill_matrix(const alea_system_t* sys,
                                       uint32_t cell_index);
bool alea_hier_spatial_index_needs_rebuild(const alea_system_t* sys);

/* Debug: total candidates pushed across all universes in the current query. */
void alea_hier_debug_candidates_reset(void);
size_t alea_hier_debug_candidates_get(void);
int alea_hier_spatial_query_region(alea_system_t* sys,
                                   const alea_bbox_t* query_bbox,
                                   alea_spatial_hit_t* out_hits,
                                   size_t max_hits);
int alea_hier_spatial_query_region_direct(alea_system_t* sys,
                                          const alea_bbox_t* query_bbox,
                                          alea_spatial_hit_t* out_hits,
                                          size_t max_hits);
int alea_hier_spatial_query_region_chain(alea_system_t* sys,
                                         const alea_bbox_t* query_bbox,
                                         alea_hier_spatial_chain_hit_t* out_hits,
                                         size_t max_hits);
int alea_hier_spatial_query_slice_z(alea_system_t* sys,
                                    double z,
                                    double x_min,
                                    double x_max,
                                    double y_min,
                                    double y_max,
                                    alea_spatial_hit_t* out_hits,
                                    size_t max_hits);

#endif /* ALEA_SPATIAL_HIER_H */
