// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/**
 * @file raycast.c
 * @brief Ray casting through CSG geometry
 *
 * Main raycast logic: find surface intersections, determine cell segments.
 */

#include "util/ray_profile.h"
#include "raycast.h"
#include "ray_traversal.h"
#include "regional/interval_index.h"
#include "ray_substep.h"
#include "core/alea_occurrence.h"
#include "alea.h"
#include "ray_intersect.h"
#include "recovery_context.h"
#include "recovery_planes.h"
#include "ray_epsilon.h"
#include "ray_bbox.h"
#include "bvh.h"
#include "core/alea_system.h"
#include "core/alea_universe.h"
#include "util/compat.h"
#include "core/alea_spatial_hier.h"
#include "geo_validator/transition_validation.h"
#include "core/alea_eval.h"
#include "util/alea_log.h"
#include "util/alea_parallel.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <limits.h>
#include "util/math.h"
#include <stdio.h>

#define INITIAL_CAPACITY 32
#define MAX_RAY_WALK_ITERATIONS 10000000

/* ============================================================================
 * CACHE PRE-BUILD (thread safety)
 *
 * All shared caches (surface BVH, spatial index, cell adjacency) must be
 * built before any concurrent access.  Call once from each public entry
 * point so that worker threads never race on lazy initialisation.
 * ============================================================================ */

int alea_raycast_ensure_caches(alea_system_t* sys) {
    return alea_system_prepare_query_caches(sys, ALEA_CACHE_RAYCAST);
}

int alea_raycast_ensure_hier_caches(alea_system_t* sys) {
    return alea_system_prepare_query_caches(sys, ALEA_CACHE_RAYCAST_HIER);
}

static bool raycast_prefers_hier_mode(const alea_system_t* sys) {
    (void)sys;
    return true;  /* hierarchical spatial index is the only backend */
}

static int raycast_flat_cell_aware_impl(alea_system_t* sys,
                                        const alea_ray_t* ray,
                                        double effective_t_max,
                                        alea_raycast_result_t* result);
static int ray_selected_interval_trace(alea_system_t* sys,
                                       const alea_ray_t* ray, double t_max,
                                       bool emit_hits,
                                       alea_raycast_result_t* result);
static void boundary_event_world_normal(const alea_matrix_t* mat,
                                        double nlx, double nly, double nlz,
                                        double* nx, double* ny, double* nz);

/* ============================================================================
 * RAY UTILITIES
 * ============================================================================ */

int alea_ray_init(alea_ray_t* ray,
                  double ox, double oy, double oz,
                  double dx, double dy, double dz) {
    ray->ox = ox;
    ray->oy = oy;
    ray->oz = oz;

    /* Normalize direction */
    double len = sqrt(dx * dx + dy * dy + dz * dz);
    if (len > RAY_EPSILON) {
        ray->dx = dx / len;
        ray->dy = dy / len;
        ray->dz = dz / len;
    } else {
        /* Zero-length direction: set fallback so struct is valid */
        ray->dx = 0;
        ray->dy = 0;
        ray->dz = 1;
    }

    /* IEEE 754: 1.0/0.0 = +inf, 1.0/(-0.0) = -inf — correct for slab tests */
    ray->inv_dx = 1.0 / ray->dx;
    ray->inv_dy = 1.0 / ray->dy;
    ray->inv_dz = 1.0 / ray->dz;

    return (len > RAY_EPSILON) ? 0 : -1;
}

static void ray_init_parametric(alea_ray_t* ray,
                                double ox, double oy, double oz,
                                double dx, double dy, double dz) {
    ray->ox = ox; ray->oy = oy; ray->oz = oz;
    ray->dx = dx; ray->dy = dy; ray->dz = dz;
    ray->inv_dx = 1.0 / dx;
    ray->inv_dy = 1.0 / dy;
    ray->inv_dz = 1.0 / dz;
}

void alea_ray_init_normalized(alea_ray_t* ray,
                              double ox, double oy, double oz,
                              double dx, double dy, double dz) {
    ray_init_parametric(ray, ox, oy, oz, dx, dy, dz);
}

/* Sample beyond the rounding uncertainty of O + t D. Near a small cell,
 * cancellation can make the resulting coordinate much smaller than either
 * operand: scaling by that coordinate alone can leave the sample on the
 * wrong side of the boundary, even after advancing t by one ULP. */
static double raycast_forward_sample_t(const alea_ray_t* ray, double t) {
    double x, y, z;
    alea_ray_point_at(ray, t, &x, &y, &z);
    double coordinate_scale = fmax(1.0,
        fmax(fabs(x), fmax(fabs(y), fabs(z))));
    const double direction_scale = fmax(fabs(ray->dx),
        fmax(fabs(ray->dy), fabs(ray->dz)));
    coordinate_scale = fmax(coordinate_scale,
        fmax(fabs(ray->ox), fmax(fabs(ray->oy), fabs(ray->oz))));
    coordinate_scale = fmax(coordinate_scale, fabs(t) * direction_scale);
    double dt = 16.0 * DBL_EPSILON * coordinate_scale / direction_scale;
    double sample = t + dt;
    if (!(sample > t)) sample = nextafter(t, INFINITY);
    return sample;
}

/* ============================================================================
 * RESULT MANAGEMENT
 * ============================================================================ */

void alea_raycast_result_init(alea_raycast_result_t* result) {
    memset(result, 0, sizeof(*result));
}

/* Zero every per-ray counter (everything in the result except the ray and the
 * hit/segment vectors). Shared by free and clear so new counters only need to
 * be added in one place. */
static void raycast_result_reset_counters(alea_raycast_result_t* result) {
    result->surfaces_tested = 0;
    result->bbox_culled = 0;
    result->point_lookups = 0;
    result->step_iterations = 0;
    result->blas_placement_candidates = 0;
    result->blas_placements_pruned = 0;
    result->blas_universe_queries = 0;
    result->blas_cell_candidates = 0;
    result->blas_cells_tested = 0;
    result->blas_hits_before_dedup = 0;
    result->owner_neighbor_attempts = 0;
    result->owner_neighbor_hits = 0;
    result->owner_path_attempts = 0;
    result->owner_path_hits = 0;
    result->owner_root_queries = 0;
    result->owner_root_hits = 0;
    result->owner_full_queries = 0;
    result->owner_full_hits = 0;
    result->boundary_event_enrichments = 0;
    result->path_snapshot_copies = 0;
    result->path_snapshot_entries = 0;
    result->selected_intervals_yielded = 0;
    result->result_buffer_growths = 0;
    result->result_buffer_growth_bytes = 0;
    result->terminal_surfaces_tested = 0;
    result->lattice_surfaces_tested = 0;
    result->ancestor_surfaces_tested = 0;
    result->crossed_cell_count = 0;
    result->max_cell_surface_count = 0;
    result->sum_cell_surface_count = 0;
    memset(result->prim_type_tests, 0, sizeof(result->prim_type_tests));
    result->ancestor_surface_queries = 0;
    result->ancestor_path_entries_examined = 0;
    result->ancestor_max_path_depth = 0;
    result->ancestor_unattributed_surface_tests = 0;
    result->ancestor_unattributed_queries = 0;
    result->ancestor_unattributed_winning_events = 0;
    memset(result->ancestor_hot_cells, 0, sizeof(result->ancestor_hot_cells));
    result->lattice_entry_calls = 0;
    result->lattice_entry_tlas_nodes_tested = 0;
    result->lattice_entry_tlas_leaves_visited = 0;
    result->lattice_entry_candidates = 0;
    result->lattice_entry_dda_steps = 0;
    result->lattice_entry_no_entry_results = 0;
    result->lattice_entry_future_entry_results = 0;
    result->lattice_entry_already_inside_results = 0;
    result->lattice_entry_ancestor_surface_tests = 0;
    result->lattice_entry_ancestor_events = 0;
    result->lattice_entry_canonical_rejections = 0;
}

void alea_raycast_result_free(alea_raycast_result_t* result) {
    recovery_context_destroy(result->recovery_roots);
    result->recovery_roots = NULL;
    alea_vec_free(&result->hits);
    alea_vec_free(&result->segments);
    alea_vec_free(&result->paths);
    alea_vec_free(&result->path_entries);
    raycast_result_reset_counters(result);
}

void alea_raycast_result_clear(alea_raycast_result_t* result) {
    recovery_context_destroy(result->recovery_roots);
    result->recovery_roots = NULL;
    alea_vec_clear(&result->hits);
    alea_vec_clear(&result->segments);
    alea_vec_clear(&result->paths);
    alea_vec_clear(&result->path_entries);
    raycast_result_reset_counters(result);
}

void alea_raycast_result_reserve(alea_raycast_result_t* result,
                                size_t hit_cap, size_t seg_cap) {
    alea_vec_reserve(&result->hits, hit_cap, alea_ray_hit_t);
    alea_vec_reserve(&result->segments, seg_cap, alea_ray_segment_t);
}

static void record_result_buffer_growth(alea_raycast_result_t* result,
                                        size_t old_capacity,
                                        size_t new_capacity,
                                        size_t element_size) {
    if (new_capacity <= old_capacity) return;
    result->result_buffer_growths++;
    result->result_buffer_growth_bytes +=
        (uint64_t)(new_capacity - old_capacity) * element_size;
}

static int add_hit(alea_raycast_result_t* result, const alea_ray_hit_t* hit) {
    /* Preserve the first bounded-enumeration error while unwinding. */
    if (result->breakpoint_hit_limit && result->breakpoint_failed) return -1;
    if (result->breakpoint_hit_limit &&
        (!isfinite(hit->t) || hit->t < 0.0)) {
        result->breakpoint_failed = 1;
        alea_set_error_detail(ALEA_ERR_INVALID_STATE,
                              "particle interval breakpoint is not finite");
        return -1;
    }
    if (result->breakpoint_nearest_only) {
        if (hit->t <= nextafter(result->breakpoint_t_min, INFINITY)) return 0;
        if (result->hits.count) {
            if (hit->t < result->hits.data[0].t) result->hits.data[0] = *hit;
            return 0;
        }
    }
    if (result->breakpoint_hit_limit &&
        result->hits.count >= result->breakpoint_hit_limit) {
        result->breakpoint_failed = 1;
        alea_set_error_detail(ALEA_ERR_OVERFLOW,
                              "particle interval breakpoint budget exceeded");
        return -1;
    }
    size_t old_capacity = result->hits.capacity;
    int res = alea_vec_push(&result->hits, *hit, alea_ray_hit_t);
    if (res != 0) {
        result->breakpoint_failed = 1;
        if (result->breakpoint_hit_limit)
            alea_set_error_detail(ALEA_ERR_OUT_OF_MEMORY,
                                  "failed to store particle interval breakpoint");
    }
    record_result_buffer_growth(result, old_capacity, result->hits.capacity,
                                sizeof(alea_ray_hit_t));
    return res != 0 ? -1 : 0;
}

static bool raycast_primitive_copy_payload(const alea_system_t* sys,
                                           uint32_t primitive_id,
                                           alea_primitive_type_t type,
                                           alea_primitive_data_t* out) {
    if (!out) return false;
    const void* payload = alea_primitive_payload_const(sys, primitive_id);
    if (!payload) return false;

    switch (type) {
        case ALEA_PRIMITIVE_PLANE:       out->plane = *(const alea_plane_data_t*)payload; return true;
        case ALEA_PRIMITIVE_SPHERE:      out->sphere = *(const alea_sphere_data_t*)payload; return true;
        case ALEA_PRIMITIVE_CYLINDER_X:  out->cyl_x = *(const alea_cylinder_x_data_t*)payload; return true;
        case ALEA_PRIMITIVE_CYLINDER_Y:  out->cyl_y = *(const alea_cylinder_y_data_t*)payload; return true;
        case ALEA_PRIMITIVE_CYLINDER_Z:  out->cyl_z = *(const alea_cylinder_z_data_t*)payload; return true;
        case ALEA_PRIMITIVE_CONE_X:      out->cone_x = *(const alea_cone_x_data_t*)payload; return true;
        case ALEA_PRIMITIVE_CONE_Y:      out->cone_y = *(const alea_cone_y_data_t*)payload; return true;
        case ALEA_PRIMITIVE_CONE_Z:      out->cone_z = *(const alea_cone_z_data_t*)payload; return true;
        case ALEA_PRIMITIVE_RPP:         out->box = *(const alea_box_data_t*)payload; return true;
        case ALEA_PRIMITIVE_QUADRIC:     out->quadric = *(const alea_quadric_data_t*)payload; return true;
        case ALEA_PRIMITIVE_TORUS_X:
        case ALEA_PRIMITIVE_TORUS_Y:
        case ALEA_PRIMITIVE_TORUS_Z:     out->torus = *(const alea_torus_data_t*)payload; return true;
        case ALEA_PRIMITIVE_RCC:         out->rcc = *(const alea_rcc_data_t*)payload; return true;
        case ALEA_PRIMITIVE_BOX:         out->box_general = *(const alea_box_general_data_t*)payload; return true;
        case ALEA_PRIMITIVE_SPH:         out->sph = *(const alea_sph_data_t*)payload; return true;
        case ALEA_PRIMITIVE_TRC:         out->trc = *(const alea_trc_data_t*)payload; return true;
        case ALEA_PRIMITIVE_ELL:         out->ell = *(const alea_ell_data_t*)payload; return true;
        case ALEA_PRIMITIVE_REC:         out->rec = *(const alea_rec_data_t*)payload; return true;
        case ALEA_PRIMITIVE_WED:         out->wed = *(const alea_wed_data_t*)payload; return true;
        case ALEA_PRIMITIVE_RHP:         out->rhp = *(const alea_rhp_data_t*)payload; return true;
        case ALEA_PRIMITIVE_ARB:         out->arb = *(const alea_arb_data_t*)payload; return true;
        default: return false;
    }
}

static int reserve_batch_budget(atomic_uint_fast64_t* counter, uint64_t limit,
                                uint64_t amount) {
    if (!counter || limit == 0 || amount == 0) return 0;
    uint_fast64_t used = atomic_load_explicit(counter, memory_order_relaxed);
    for (;;) {
        if (used > limit || amount > limit - used) return -1;
        if (atomic_compare_exchange_strong(counter, &used, used + amount)) {
            return 0;
        }
    }
}

static int add_segment(alea_raycast_result_t* result, const alea_ray_segment_t* seg) {
    if (reserve_batch_budget(result->segment_counter, result->segment_limit, 1) != 0) {
        result->segment_limit_exceeded = 1;
        return -1;
    }
    if (result->path_entry_counter && result->path_entry_limit != 0 &&
        seg->path_index != UINT32_MAX && seg->path_index < result->paths.count) {
        uint64_t path_count = result->paths.data[seg->path_index].count;
        if (reserve_batch_budget(result->path_entry_counter,
                                 result->path_entry_limit, path_count) != 0) {
            result->path_entry_limit_exceeded = 1;
            return -1;
        }
    }
    size_t old_capacity = result->segments.capacity;
    int res = alea_vec_push(&result->segments, *seg, alea_ray_segment_t);
    record_result_buffer_growth(result, old_capacity, result->segments.capacity,
                                sizeof(alea_ray_segment_t));
    return res != 0 ? -1 : 0;
}

static uint64_t path_entry_occurrence_key(
    alea_occurrence_state_t* state,
    const alea_hier_ray_path_entry_t* entry) {
    return alea_occurrence_state_step(
        state, entry->cell_index, entry->is_lattice,
        entry->lat_i, entry->lat_j, entry->lat_k);
}

static int ray_path_entry_equal(const alea_ray_path_entry_t* a,
                                const alea_ray_path_entry_t* b) {
    return a->cell_index == b->cell_index &&
           a->cell_id == b->cell_id &&
           a->material_id == b->material_id &&
           a->universe_id == b->universe_id &&
           a->fill_universe == b->fill_universe &&
           a->depth == b->depth &&
           a->is_lattice == b->is_lattice &&
           a->lattice_origin[0] == b->lattice_origin[0] &&
           a->lattice_origin[1] == b->lattice_origin[1] &&
           a->lattice_origin[2] == b->lattice_origin[2] &&
           a->occurrence_key == b->occurrence_key;
}

static int capture_hier_path(alea_raycast_result_t* result,
                             const alea_hier_ray_path_t* path,
                             uint32_t* out_path_index) {
    *out_path_index = UINT32_MAX;
    if (!result->capture_paths || !path || path->count <= 0) return 0;

    size_t count = (size_t)path->count;
    alea_ray_path_entry_t entries[ALEA_HIER_RAY_PATH_MAX];
    alea_occurrence_state_t occurrence_state;
    alea_occurrence_state_init(&occurrence_state);
    for (size_t i = 0; i < count; ++i) {
        const alea_hier_ray_path_entry_t* src = &path->entries[i];
        alea_ray_path_entry_t* dst = &entries[i];
        dst->cell_index = src->cell_index;
        dst->cell_id = src->cell_id;
        dst->material_id = src->material_id;
        dst->universe_id = src->universe_id;
        dst->fill_universe = src->fill_universe;
        dst->depth = src->depth;
        dst->is_lattice = src->is_lattice;
        dst->lattice_origin[0] = src->lat_ox;
        dst->lattice_origin[1] = src->lat_oy;
        dst->lattice_origin[2] = src->lat_oz;
        dst->occurrence_key = path_entry_occurrence_key(
            &occurrence_state, src);
    }

    for (size_t i = 0; i < result->paths.count; ++i) {
        const alea_ray_path_t* candidate = &result->paths.data[i];
        if (candidate->count != count) continue;
        int same = 1;
        for (size_t j = 0; j < count; ++j) {
            if (!ray_path_entry_equal(&result->path_entries.data[candidate->offset + j],
                                      &entries[j])) {
                same = 0;
                break;
            }
        }
        if (same) {
            *out_path_index = (uint32_t)i;
            return 0;
        }
    }

    alea_ray_path_t record = {
        .offset = (uint32_t)result->path_entries.count,
        .count = (uint16_t)count,
    };
    for (size_t i = 0; i < count; ++i) {
        size_t old_capacity = result->path_entries.capacity;
        int push_result = alea_vec_push(&result->path_entries, entries[i],
                                        alea_ray_path_entry_t);
        record_result_buffer_growth(result, old_capacity,
                                    result->path_entries.capacity,
                                    sizeof(alea_ray_path_entry_t));
        if (push_result != 0)
            return -1;
    }
    size_t old_capacity = result->paths.capacity;
    int push_result = alea_vec_push(&result->paths, record, alea_ray_path_t);
    record_result_buffer_growth(result, old_capacity, result->paths.capacity,
                                sizeof(alea_ray_path_t));
    if (push_result != 0)
        return -1;
    *out_path_index = (uint32_t)(result->paths.count - 1);
    return 0;
}

/* ============================================================================
 * COMPARISON FOR SORTING
 * ============================================================================ */

static int compare_hits(const void* a, const void* b) {
    const alea_ray_hit_t* ha = (const alea_ray_hit_t*)a;
    const alea_ray_hit_t* hb = (const alea_ray_hit_t*)b;
    if (ha->t < hb->t) return -1;
    if (ha->t > hb->t) return 1;
    /* When t values are equal, put DDA boundary hits (surface_id=0) LAST.
     * This ensures the segment builder correctly identifies intervals after
     * coincident surface+DDA hits as crossing a DDA boundary, not a surface.
     * Critical for correct lattice element transitions on all platforms. */
    if (ha->surface_id == 0 && hb->surface_id != 0) return 1;
    if (ha->surface_id != 0 && hb->surface_id == 0) return -1;
    return 0;
}

/* Inline comparison for insertion sort (avoids function-call overhead) */
static inline int hit_less(const alea_ray_hit_t* a, const alea_ray_hit_t* b) {
    if (a->t != b->t) return a->t < b->t;
    /* DDA boundary hits (surface_id=0) sort after real surface hits */
    return (a->surface_id != 0 && b->surface_id == 0);
}

/**
 * Sort hits by distance. Uses insertion sort for small arrays (<=64),
 * qsort for larger arrays. Insertion sort has lower overhead and is
 * cache-friendly for the small arrays typical of raycast results.
 */
static void sort_hits(alea_ray_hit_t* hits, size_t count) {
    if (count <= 1) return;

    if (count <= 64) {
        /* Insertion sort — O(n^2) but faster than qsort for small n */
        for (size_t i = 1; i < count; i++) {
            alea_ray_hit_t key = hits[i];
            size_t j = i;
            while (j > 0 && hit_less(&key, &hits[j - 1])) {
                hits[j] = hits[j - 1];
                j--;
            }
            hits[j] = key;
        }
    } else {
        qsort(hits, count, sizeof(alea_ray_hit_t), compare_hits);
    }
}

/* ray_bbox_slab and ray_bbox_slab_enter_exit are in ray_bbox.h */

/* ============================================================================
 * MAIN RAYCAST - SURFACE INTERSECTIONS
 * ============================================================================ */

/* Context for BVH traversal callback */
typedef struct {
    alea_system_t* sys;
    const alea_ray_t* ray;
    double t_min;
    double t_max;
    alea_raycast_result_t* result;
    bool failed;
} bvh_raycast_ctx_t;

/* BVH batch callback: test surface intersections for a leaf node */
static void bvh_surface_batch_callback(const uint32_t* surface_indices,
                                        uint32_t count, void* userdata) {
    bvh_raycast_ctx_t* ctx = (bvh_raycast_ctx_t*)userdata;
    if (ctx->failed) return;

    for (uint32_t si = 0; si < count; si++) {
        uint32_t surface_idx = surface_indices[si];
        const alea_surface_entry_t* surf = &ctx->sys->surfaces.data[surface_idx];
        const alea_primitive_entry_t* prim = &ctx->sys->primitives.data[surf->primitive_id];
        alea_primitive_data_t prim_data;
        if (!raycast_primitive_copy_payload(ctx->sys, surf->primitive_id,
                                            prim->type, &prim_data)) continue;

        ctx->result->surfaces_tested++;

        /* Find intersections */
        double t[4];
        int n = ray_intersect_primitive(ctx->ray, prim->type, &prim_data, t);

        /* Add valid hits */
        for (int j = 0; j < n; j++) {
            if (t[j] >= ctx->t_min && t[j] <= ctx->t_max) {
                alea_ray_hit_t hit;
                hit.t = t[j];
                hit.surface_id = surf->mc_surface_id;
                hit.primitive_id = surf->primitive_id;

                /* Compute normal at hit point */
                double px, py, pz;
                alea_ray_point_at(ctx->ray, t[j], &px, &py, &pz);
                primitive_normal_at(prim->type, &prim_data, px, py, pz,
                                   &hit.nx, &hit.ny, &hit.nz);

                if (add_hit(ctx->result, &hit) != 0) {
                    if (!ctx->result->breakpoint_hit_limit)
                        alea_set_error_detail(ALEA_ERR_OUT_OF_MEMORY,
                                              "failed to store BVH surface hit");
                    ctx->failed = true;
                    return;
                }
            }
        }
    }
}

/* Linear scan fallback (used when BVH not available) */
static int raycast_surfaces_linear(alea_system_t* sys,
                                   const alea_ray_t* ray,
                                   double t_min, double t_max,
                                   alea_raycast_result_t* result) {
    for (size_t i = 0; i < alea_vec_count(&sys->surfaces); i++) {
        
        const alea_surface_entry_t* surf = &sys->surfaces.data[i];
        const alea_primitive_entry_t* prim = &sys->primitives.data[surf->primitive_id];
        alea_primitive_data_t prim_data;
        if (!raycast_primitive_copy_payload(sys, surf->primitive_id,
                                            prim->type, &prim_data)) continue;

        result->surfaces_tested++;

        double t[4];
        int count = ray_intersect_primitive(ray, prim->type, &prim_data, t);

        for (int j = 0; j < count; j++) {
            if (t[j] >= t_min && t[j] <= t_max) {
                alea_ray_hit_t hit;
                hit.t = t[j];
                hit.surface_id = surf->mc_surface_id;
                hit.primitive_id = surf->primitive_id;

                double px, py, pz;
                alea_ray_point_at(ray, t[j], &px, &py, &pz);
                primitive_normal_at(prim->type, &prim_data, px, py, pz,
                                   &hit.nx, &hit.ny, &hit.nz);

                if (add_hit(result, &hit) != 0) {
                    if (!result->breakpoint_hit_limit)
                        alea_set_error_detail(ALEA_ERR_OUT_OF_MEMORY,
                                              "failed to store linear surface hit");
                    return -1;
                }
            }
        }
    }
    return 0;
}

/* Shared implementation: intersect, sort, dedup */
static int raycast_surfaces_impl(alea_system_t* sys,
                                  const alea_ray_t* ray,
                                  double t_min, double t_max,
                                  alea_raycast_result_t* result) {
    alea_raycast_result_clear(result);
    result->ray = *ray;

    /* Use BVH if available, otherwise fall back to linear scan */
    if (!BVH_DISABLED && sys->surface_bvh) {
        bvh_raycast_ctx_t ctx = {
            .sys = sys,
            .ray = ray,
            .t_min = t_min,
            .t_max = t_max,
            .result = result
        };
        int traversed = alea_bvh_traverse_batch(
            sys->surface_bvh, ray, t_min, t_max,
            bvh_surface_batch_callback, &ctx);
        if (traversed < 0 || ctx.failed) {
            result->hits.count = 0;
            return -1;
        }
    } else {
        if (raycast_surfaces_linear(sys, ray, t_min, t_max, result) != 0) {
            result->hits.count = 0;
            return -1;
        }
    }

    /* Sort hits by distance */
    sort_hits(result->hits.data, result->hits.count);

    /* Remove duplicate hits (same t AND same surface_id within epsilon) */
    if (!result->breakpoint_hit_limit && result->hits.count > 1) {
        size_t write = 1;
        for (size_t read = 1; read < result->hits.count; read++) {
            int same_t = fabs(result->hits.data[read].t - result->hits.data[write - 1].t) <= DEDUP_EPSILON;
            int same_surf = result->hits.data[read].surface_id == result->hits.data[write - 1].surface_id;
            if (!(same_t && same_surf)) {
                result->hits.data[write++] = result->hits.data[read];
            }
        }
        result->hits.count = write;
    }

    return 0;
}

int alea_raycast_surfaces(alea_system_t* sys,
                         const alea_ray_t* ray,
                         double t_min, double t_max,
                         alea_raycast_result_t* result) {
    if (!sys || !ray || !result) return -1;

    if (!alea_system_query_cache_ready(sys, ALEA_CACHE_RAYCAST)) {
        alea_set_error_detail(ALEA_ERR_INVALID_STATE,
                              "raycast caches are not prepared; call alea_prepare_query_acceleration()");
        return -1;
    }

    return raycast_surfaces_impl(sys, ray, t_min, t_max, result);
}

int alea_raycast_surfaces_nocache(alea_system_t* sys,
                                  const alea_ray_t* ray,
                                  double t_min, double t_max,
                                  alea_raycast_result_t* result) {
    if (!sys || !ray || !result) return -1;
    return raycast_surfaces_impl(sys, ray, t_min, t_max, result);
}

/* ============================================================================
 * CONVERT HITS TO CELL SEGMENTS (with cell adjacency optimization)
 * ============================================================================ */

/**
 * @brief Find cell via neighbor lookup (O(1) when adjacency is built)
 * @return 1 if found, 0 if full lookup needed
 */
static int raycast_find_neighbor(alea_system_t* sys,
                                 int current_cell_idx,
                                 int surface_id,
                                 double px, double py, double pz,
                                 int* out_cell_id, int* out_cell_idx,
                                 int* out_material_id, double* out_density) {
    if (!sys->cell_adjacency_built || current_cell_idx < 0) {
        return 0;
    }

    /* A surface may bound several cells, and crossing a surface referenced
     * by a non-convex cell does not always leave that cell, so the first
     * adjacency entry for surface_id is only a candidate: accept a neighbor
     * only if it contains the sample point past the crossing. Callers gate
     * this on universe-0 cells, whose local frame is the world frame. */
    const alea_cell_entry_t* cell = &sys->cells.data[current_cell_idx];
    for (size_t i = 0; i < cell->neighbor_count; i++) {
        if (cell->neighbors[i].surface_id != surface_id) continue;
        uint32_t nb_idx = cell->neighbors[i].neighbor_index;
        if (nb_idx >= alea_vec_count(&sys->cells)) continue;
        const alea_cell_entry_t* neighbor = &sys->cells.data[nb_idx];
        if (neighbor->root_node_id == ALEA_NODE_ID_INVALID) continue;
        if (!alea_contains_point(sys, neighbor->root_node_id, px, py, pz))
            continue;
        *out_cell_id = neighbor->mc_cell_id;
        *out_cell_idx = (int)nb_idx;
        *out_material_id = neighbor->material_id;
        *out_density = neighbor->density;
        return 1;
    }
    return 0;  /* No containing neighbor (re-entrant crossing/exterior) */
}

static int cell_references_surface_id(const alea_system_t* sys,
                                      const alea_cell_entry_t* cell,
                                      int surface_id) {
    if (!sys || !cell || surface_id <= 0 ||
        !cell->surface_indices || cell->surface_index_count == 0) {
        return 0;
    }

    for (size_t i = 0; i < cell->surface_index_count; i++) {
        uint32_t surf_idx = cell->surface_indices[i];
        if (surf_idx >= alea_vec_count(&sys->surfaces))
            continue;
        if (sys->surfaces.data[surf_idx].mc_surface_id == surface_id)
            return 1;
    }
    return 0;
}

/* Hits from the all-cells queries are in DFS preorder. Follow the first
 * strictly-deepening chain and take its deepest hit, so overlapping
 * same-depth cells resolve to the first containing cell in deck order —
 * matching the canonical resolver and the hier descent. */
static int deepest_first_chain_index(const alea_cell_hit_t* hits, int n) {
    int idx = 0;
    while (idx + 1 < n && hits[idx + 1].depth == hits[idx].depth + 1) idx++;
    return idx;
}

/**
 * @brief Full cell lookup via point-in-cell search
 */
static void raycast_find_cell_full(alea_system_t* sys,
                                   double px, double py, double pz,
                                   int* out_cell_id, int* out_cell_idx,
                                   int* out_material_id, double* out_density) {
    *out_cell_id = -1;
    *out_cell_idx = -1;
    *out_material_id = 0;
    *out_density = 0;

    alea_cell_hit_t hits[32];
    int num_hits = alea_find_all_cells_at_point(sys, px, py, pz, hits, 32);

    if (num_hits > 0) {
        const alea_cell_hit_t* hit =
            &hits[deepest_first_chain_index(hits, num_hits)];
        *out_cell_id = hit->cell_id;
        *out_cell_idx = hit->cell_index;
        *out_material_id = hit->material_id;

        if (hit->cell_index >= 0 && (size_t)hit->cell_index < alea_vec_count(&sys->cells)) {
            *out_density = sys->cells.data[hit->cell_index].density;
        }
    }
}

static void raycast_find_cell_full_hier(alea_system_t* sys,
                                        double px, double py, double pz,
                                        int* out_cell_id, int* out_cell_idx,
                                        int* out_material_id, double* out_density) {
    *out_cell_id = -1;
    *out_cell_idx = -1;
    *out_material_id = 0;
    *out_density = 0;

    alea_cell_hit_t hits[32];
    int num_hits = alea_hier_spatial_find_cells_at_point(sys, px, py, pz,
                                                         hits, 32);

    if (num_hits > 0) {
        const alea_cell_hit_t* hit =
            &hits[deepest_first_chain_index(hits, num_hits)];
        *out_cell_id = hit->cell_id;
        *out_cell_idx = hit->cell_index;
        *out_material_id = hit->material_id;

        if (hit->cell_index >= 0 && (size_t)hit->cell_index < alea_vec_count(&sys->cells)) {
            *out_density = sys->cells.data[hit->cell_index].density;
        }
    }
}

/**
 * 3-tier cell lookup after crossing a surface:
 * 1) Neighbor lookup (O(1) via adjacency)
 * 2) Coherence check (O(tree_depth) — is point still in previous cell?)
 * 3) Full point-in-cell search (O(n_candidates) via spatial index)
 */
static void find_cell_after_crossing(alea_system_t* sys,
                                     const alea_ray_t* ray,
                                     double t_prev, double t_curr,
                                     int prev_cell_idx,
                                     int crossed_surface_id,
                                     bool use_hier_lookup,
                                     int* out_cell_id, int* out_cell_idx,
                                     int* out_material_id, double* out_density) {
    *out_cell_id = -1;
    *out_cell_idx = -1;
    *out_material_id = 0;
    *out_density = 0;
    int found = 0;

    /* Tier 1: neighbor lookup */
    if (prev_cell_idx >= 0 && crossed_surface_id > 0) {
        const alea_cell_entry_t* prev_cell = &sys->cells.data[prev_cell_idx];
        if (prev_cell->universe_id == 0 &&
            cell_references_surface_id(sys, prev_cell, crossed_surface_id)) {
            double t_sample = t_prev +
                fmin(0.5 * (t_curr - t_prev), SURFACE_SAMPLE_OFFSET);
            double px, py, pz;
            alea_ray_point_at(ray, t_sample, &px, &py, &pz);
            found = raycast_find_neighbor(sys, prev_cell_idx, crossed_surface_id,
                                          px, py, pz,
                                          out_cell_id, out_cell_idx,
                                          out_material_id, out_density);
            /* In hier mode a filled container is not a terminal answer:
             * reject it so the full lookup descends into the fill. */
            if (found && use_hier_lookup) {
                const alea_cell_entry_t* nb = &sys->cells.data[*out_cell_idx];
                if (nb->fill_universe > 0 ||
                    (nb->lat_type != 0 && nb->lat_fill)) {
                    found = 0;
                    *out_cell_id = -1;
                    *out_cell_idx = -1;
                    *out_material_id = 0;
                    *out_density = 0;
                }
            }
        }
    }

    int crossed_dda = (crossed_surface_id == 0);

    /* Tier 2: coherence check (skip for DDA boundaries) */
    if (!found && prev_cell_idx >= 0 && !crossed_dda) {
        const alea_cell_entry_t* prev_cell = &sys->cells.data[prev_cell_idx];

        /* Nested-cell case: in flat mode Tier 3 is already fast via the
         * flat point-query coherence cache, so skip Tier 2 there. In hier
         * mode we consult the hier cache for the cell's world→local
         * transform. */
        if (prev_cell->universe_id != 0 && !use_hier_lookup)
            goto full_lookup;

        double t_sample = t_prev + fmin(0.5 * (t_curr - t_prev), SURFACE_SAMPLE_OFFSET);
        double px, py, pz;
        alea_ray_point_at(ray, t_sample, &px, &py, &pz);

        int in_cell;
        if (prev_cell->universe_id == 0) {
            in_cell = alea_contains_point(sys, prev_cell->root_node_id,
                                          px, py, pz) ? 1 : 0;
        } else {
            in_cell = alea_hier_spatial_check_cached_containment(
                sys, (uint32_t)prev_cell_idx, px, py, pz);
            if (in_cell < 0) goto full_lookup;
        }

        if (in_cell) {
            *out_cell_id = prev_cell->mc_cell_id;
            *out_cell_idx = prev_cell_idx;
            *out_material_id = prev_cell->material_id;
            *out_density = prev_cell->density;
            found = 1;
        }
    }

    /* Tier 3: full lookup */
full_lookup:
    if (!found) {
        double max_offset = crossed_dda ? DDA_SAMPLE_OFFSET : SURFACE_SAMPLE_OFFSET;
        double t_sample = t_prev + fmin(0.5 * (t_curr - t_prev), max_offset);
        double px, py, pz;
        alea_ray_point_at(ray, t_sample, &px, &py, &pz);
        if (use_hier_lookup) {
            raycast_find_cell_full_hier(sys, px, py, pz,
                                        out_cell_id, out_cell_idx,
                                        out_material_id, out_density);
        } else {
            raycast_find_cell_full(sys, px, py, pz,
                                   out_cell_id, out_cell_idx,
                                   out_material_id, out_density);
        }
    }
}

static int raycast_to_segments_impl(alea_system_t* sys,
                                    double t_max,
                                    alea_raycast_result_t* result,
                                    bool use_hier_lookup) {
    if (!sys || !result) return -1;

    result->segments.count = 0;

    const alea_ray_t* ray = &result->ray;
    double t_prev = 0;
    int prev_cell_id = -2;  /* Use -2 as "no previous" since -1 is valid (void) */
    int prev_cell_idx = -1; /* Track cell index for neighbor lookup */
    double t_terminus = (t_max > 0) ? t_max : DBL_MAX;

    /* Process intervals between hits */
    for (size_t i = 0; i <= result->hits.count; i++) {
        double t_curr = (i < result->hits.count) ? result->hits.data[i].t : t_terminus;
        if (t_curr > t_terminus) t_curr = t_terminus;

        /* Only process if there's a real interval */
        if (t_curr > t_prev + RAY_EPSILON) {
            int cell_id, cell_idx, material_id;
            double density;
            int crossed_surface = (i > 0) ? result->hits.data[i - 1].surface_id : -1;
            find_cell_after_crossing(sys, ray, t_prev, t_curr,
                                     prev_cell_idx, crossed_surface,
                                     use_hier_lookup,
                                     &cell_id, &cell_idx, &material_id, &density);

            /* Update tracking for next iteration */
            prev_cell_idx = cell_idx;

            /* Extend previous segment or start new one */
            if (cell_id == prev_cell_id && result->segments.count > 0) {
                /* Extend */
                alea_ray_segment_t* prev_seg =
                    &result->segments.data[result->segments.count - 1];
                prev_seg->t_exit = t_curr;
                prev_seg->exit_surface_id =
                    (i < result->hits.count &&
                     result->hits.data[i].t <= t_terminus + RAY_EPSILON)
                        ? result->hits.data[i].surface_id
                        : -1;
            } else {
                /* New segment */
                alea_ray_segment_t seg;
                seg.t_enter = t_prev;
                seg.t_exit = t_curr;
                seg.cell_id = cell_id;
                seg.material_id = material_id;
                seg.density = density;
                seg.enter_surface_id = (i > 0) ? result->hits.data[i - 1].surface_id : -1;
                seg.exit_surface_id =
                    (i < result->hits.count &&
                     result->hits.data[i].t <= t_terminus + RAY_EPSILON)
                        ? result->hits.data[i].surface_id
                        : -1;
                seg.enter_hit_index = (i > 0) ? (int)(i - 1) : -1;
                seg.path_index = UINT32_MAX;
                seg.resolution_flags =
                    (cell_idx >= 0 &&
                     (size_t)cell_idx < alea_vec_count(&sys->cells) &&
                     alea_cell_entry_is_container(&sys->cells.data[cell_idx]))
                        ? ALEA_RESOLVE_UNDEFINED_FILL : 0;

                if (add_segment(result, &seg) != 0) return -1;
                prev_cell_id = cell_id;
            }
        }

        t_prev = t_curr;
    }

    return 0;
}

int alea_raycast_to_segments(alea_system_t* sys,
                            double t_max,
                            alea_raycast_result_t* result) {
    return raycast_to_segments_impl(sys, t_max, result, false);
}

/* ============================================================================
 * LATTICE DDA RAYCAST
 * ============================================================================ */

/**
 * Walk a CSG tree and test ray intersection with every primitive leaf.
 * Adds hits to result. Duplicates from shared surfaces are handled by
 * the global dedup pass after all hits are collected.
 */
void alea_raycast_tree_primitives(alea_system_t* sys,
                                    const alea_ray_t* ray,
                                    alea_node_id_t node_id,
                                    double t_min, double t_max,
                                    alea_raycast_result_t* result,
                                    const alea_matrix_t* local_to_world) {
    if (alea_ray_interrupted(result)) { result->breakpoint_failed = 1; return; }
    if (result->breakpoint_hit_limit && result->breakpoint_failed) return;
    if (result->breakpoint_nearest_only && result->hits.count)
        t_max = fmin(t_max, result->hits.data[0].t);
    if (node_id >= alea_vec_count(&sys->nodes)) {
        if (result->breakpoint_hit_limit) result->breakpoint_failed = 1;
        return;
    }
    const alea_node_t* node = &sys->nodes.data[node_id];
    alea_operation_t op = ALEA_GET_OPERATION(node);

    if (op == ALEA_OP_PRIMITIVE) {
        if (node->primitive.primitive_id >= alea_vec_count(&sys->primitives)) {
            if (result->breakpoint_hit_limit) result->breakpoint_failed = 1;
            return;
        }
        const alea_primitive_entry_t* prim =
            &sys->primitives.data[node->primitive.primitive_id];
        alea_primitive_data_t prim_data;
        if (!raycast_primitive_copy_payload(sys, node->primitive.primitive_id,
                                            prim->type, &prim_data)) {
            if (result->breakpoint_hit_limit) result->breakpoint_failed = 1;
            return;
        }
        double t[4];
        result->surfaces_tested++;
        int count = ray_intersect_primitive(ray, prim->type, &prim_data, t);
        for (int j = 0; j < count; j++) {
            if (t[j] >= t_min && t[j] <= t_max) {
                alea_ray_hit_t hit;
                hit.t = t[j];
                hit.surface_id = node->primitive.mc_surface_id;
                hit.primitive_id = node->primitive.primitive_id;
                double px, py, pz;
                alea_ray_point_at(ray, t[j], &px, &py, &pz);
                primitive_normal_at(prim->type, &prim_data, px, py, pz,
                                   &hit.nx, &hit.ny, &hit.nz);
                if (local_to_world) {
                    boundary_event_world_normal(
                        local_to_world, hit.nx, hit.ny, hit.nz,
                        &hit.nx, &hit.ny, &hit.nz);
                }
                if (add_hit(result, &hit) != 0) {
                    if (!result->breakpoint_hit_limit)
                        ALEA_LOG_WARN("add_hit failed (out of memory) - raycast results may be incomplete");
                    return;
                }
            }
        }
        return;
    }

    if (op == ALEA_OP_COMPLEMENT) {
        alea_raycast_tree_primitives(sys, ray, node->operation.left,
                                t_min, t_max, result, local_to_world);
    } else {
        alea_raycast_tree_primitives(sys, ray, node->operation.left,
                                t_min, t_max, result, local_to_world);
        alea_raycast_tree_primitives(sys, ray, node->operation.right,
                                t_min, t_max, result, local_to_world);
    }
}

#ifdef ALEA_RAY_PROFILE
static void raycast_indexed_surface_scan_profile_impl(alea_system_t* sys,
                                              const alea_ray_t* ray,
                                              const uint32_t* indices, size_t count,
                                              double t_min,
                                              double t_max,
                                              alea_raycast_result_t* result,
                                              const alea_matrix_t* local_to_world,
                                              uint8_t* seen_surfaces,
                                              const recovery_plane_index_t* plane_index);
static void raycast_indexed_surface_scan(alea_system_t* sys,
                                              const alea_ray_t* ray,
                                              const uint32_t* indices, size_t count,
                                              double t_min,
                                              double t_max,
                                              alea_raycast_result_t* result,
                                              const alea_matrix_t* local_to_world,
                                              uint8_t* seen_surfaces,
                                              const recovery_plane_index_t* plane_index) {
    int previous=alea_ray_profile_enter(RAY_RECOVERY_SCAN);
    uint64_t tests_before = result->surfaces_tested;
    raycast_indexed_surface_scan_profile_impl(sys, ray, indices, count, t_min, t_max, result, local_to_world, seen_surfaces, plane_index);
    alea_ray_profile_list_record(indices, count, result->surfaces_tested - tests_before);
    alea_ray_profile_leave(previous);
}
static void raycast_indexed_surface_scan_profile_impl(alea_system_t* sys,
                                              const alea_ray_t* ray,
                                              const uint32_t* indices, size_t count,
                                              double t_min,
                                              double t_max,
                                              alea_raycast_result_t* result,
                                              const alea_matrix_t* local_to_world,
                                              uint8_t* seen_surfaces,
                                              const recovery_plane_index_t* plane_index) {
#else
static void raycast_indexed_surface_scan(alea_system_t* sys,
                                              const alea_ray_t* ray,
                                              const uint32_t* indices, size_t count,
                                              double t_min,
                                              double t_max,
                                              alea_raycast_result_t* result,
                                              const alea_matrix_t* local_to_world,
                                              uint8_t* seen_surfaces,
                                              const recovery_plane_index_t* plane_index) {
#endif
    if (result->breakpoint_hit_limit && result->breakpoint_failed) return;
    if (result->breakpoint_nearest_only && result->hits.count)
        t_max = fmin(t_max, result->hits.data[0].t);
    uint64_t plane_mask[RECOVERY_PLANE_MASK_WORDS];
    bool use_planes = result->breakpoint_nearest_only && plane_index &&
                      plane_index->cards == count;
    if (use_planes) recovery_plane_candidates(plane_index, ray, t_min, t_max, plane_mask);
    for (size_t i = 0; i < count; i++) {
        if ((i & 1023) == 0 && alea_ray_interrupted(result)) {
            result->breakpoint_failed = 1;
            return;
        }
        if (result->breakpoint_hit_limit && result->breakpoint_failed) return;
        if (use_planes && !(plane_mask[i/64] & (UINT64_C(1) << (i%64)))) continue;
        uint32_t surf_idx = indices[i];
        if (surf_idx >= alea_vec_count(&sys->surfaces)) {
            if (result->breakpoint_hit_limit) result->breakpoint_failed = 1;
            continue;
        }

        /* The same occurrence uses one local ray and a shrinking hit bound.
         * Repeated cards cannot add an earlier hit; distinct aliases remain. */
        if (seen_surfaces) {
            const uint8_t mask = (uint8_t)(1u << (surf_idx & 7));
            if (seen_surfaces[surf_idx >> 3] & mask) continue;
            seen_surfaces[surf_idx >> 3] |= mask;
        }
        const alea_surface_entry_t* surf = &sys->surfaces.data[surf_idx];
        if (surf->primitive_id >= alea_vec_count(&sys->primitives)) {
            if (result->breakpoint_hit_limit) result->breakpoint_failed = 1;
            continue;
        }

        const alea_primitive_entry_t* prim =
            &sys->primitives.data[surf->primitive_id];
        alea_primitive_data_t prim_data;
        const alea_plane_data_t* plane = NULL;
        double t[4];
        int count;
        if (result->breakpoint_nearest_only && prim->type == ALEA_PRIMITIVE_PLANE) {
            /* Most recovery planes miss the window. Use the existing kernel
             * directly and defer the normal payload copy until a root is kept. */
            plane = alea_primitive_payload_const(sys, surf->primitive_id);
            if (!plane) {
                if (result->breakpoint_hit_limit) result->breakpoint_failed = 1;
                continue;
            }
            result->surfaces_tested++;
            count = ray_intersect_plane(ray, plane, t, NULL, NULL, NULL);
        } else {
            if (!raycast_primitive_copy_payload(sys, surf->primitive_id,
                                                prim->type, &prim_data)) {
                if (result->breakpoint_hit_limit) result->breakpoint_failed = 1;
                continue;
            }

            if (result->breakpoint_nearest_only &&
                prim->type >= ALEA_PRIMITIVE_TORUS_X &&
                prim->type <= ALEA_PRIMITIVE_TORUS_Z &&
                !ray_torus_window_may_hit(ray, &prim_data.torus, t_min, t_max))
                continue;

            result->surfaces_tested++;
            count = ray_intersect_primitive(ray, prim->type, &prim_data, t);
        }
        RAY_PROFILE_COUNT(scan_tests[(unsigned)prim->type < RAY_PROFILE_TYPES ? prim->type : 0]);
        for (int j = 0; j < count; j++) {
            if (t[j] >= t_min && t[j] <= t_max) {
                RAY_PROFILE_COUNT(scan_window_roots[(unsigned)prim->type < RAY_PROFILE_TYPES ? prim->type : 0]);
                alea_ray_hit_t hit;
                hit.t = t[j];
                hit.surface_id = surf->mc_surface_id;
                hit.primitive_id = surf->primitive_id;
                double px, py, pz;
                alea_ray_point_at(ray, t[j], &px, &py, &pz);
                if (plane) prim_data.plane = *plane;
                primitive_normal_at(prim->type, &prim_data, px, py, pz,
                                    &hit.nx, &hit.ny, &hit.nz);
                if (local_to_world) {
                    boundary_event_world_normal(
                        local_to_world, hit.nx, hit.ny, hit.nz,
                        &hit.nx, &hit.ny, &hit.nz);
                }
                if (add_hit(result, &hit) != 0) {
                    if (!result->breakpoint_hit_limit)
                        ALEA_LOG_WARN("add_hit failed (out of memory) - raycast results may be incomplete");
                    return;
                }
            }
        }
    }
}

/* Build for the remaining ray, so later windows can reuse the roots even if
 * this visit has already clipped its search to a nearer hit. */
static void recovery_roots_build(alea_system_t* sys, const alea_ray_t* ray,
    const uint32_t* indices, size_t count, recovery_context_t* table,
    recovery_context_entry_t* entry, alea_raycast_result_t* result,
    const recovery_plane_index_t* plane_index) {
    entry->roots_attempted = true;
    if (count < 128 || count > RECOVERY_ROOT_BUDGET / (4*sizeof(recovery_root_t))) return;
    const double lo = result->breakpoint_t_min;
    const double hi = fmax(result->breakpoint_t_max, result->recovery_root_limit);
    if (!isfinite(lo) || !isfinite(hi) || lo < 0 || hi < lo) return;
    uint64_t plane_mask[RECOVERY_PLANE_MASK_WORDS];
    bool use_planes = plane_index && plane_index->cards == count;
    if (use_planes) recovery_plane_root_candidates(plane_index, ray, lo, hi, plane_mask);
    /* The prepared index validated this exact card list and counted its tori. */
    size_t exceptional_count = use_planes ? plane_index->torus_cards : 0;
    if (use_planes) { RAY_PROFILE_ADD(prepared_card_checks_reused,count); }
    if (!use_planes) {
        for (size_t i=0;i<count;i++) {
            if (indices[i]>=sys->surfaces.count) return;
            uint32_t id=sys->surfaces.data[indices[i]].primitive_id;
            if (id>=sys->primitives.count) return;
            alea_primitive_type_t type=sys->primitives.data[id].type;
            if (type>=ALEA_PRIMITIVE_TORUS_X && type<=ALEA_PRIMITIVE_TORUS_Z) exceptional_count++;
        }
    }
    size_t available=RECOVERY_ROOT_BUDGET-table->root_bytes;
    size_t exceptional_bytes=exceptional_count*sizeof(uint32_t);
    if (exceptional_bytes>=available) return;
    size_t capacity=count<512 ? count : 512;
    if (capacity>(available-exceptional_bytes)/sizeof(recovery_root_t)) return;
    size_t root_bytes=capacity*sizeof(recovery_root_t);
    recovery_root_t* roots=malloc(root_bytes);
    if (!roots) return; /* Optional scratch allocation: ordinary scan remains. */
    uint32_t* exceptional=exceptional_count ? malloc(exceptional_count*sizeof(*exceptional)) : NULL;
    if (exceptional_count && !exceptional) { free(roots); return; }
    uint32_t quads_left = use_planes ? plane_index->quadric_count : 0;
    const recovery_quadric_group_t* group = quads_left
        ? (const recovery_quadric_group_t*)((const char*)plane_index + plane_index->quadric_offset) : NULL;
    unsigned lane = 0;
    int quad_hits[RAY_QUADRIC_LANES];
    double quad_roots[RAY_QUADRIC_LANES][2];
    size_t n=0, ne=0;
    for (size_t i=0;i<count;i++) {
        if ((i&1023)==0 && alea_ray_interrupted(result)) { free(roots); free(exceptional); return; }
        if (use_planes && !(plane_mask[i/64] & (UINT64_C(1) << (i%64)))) continue;
        double single[4];
        const double* t = single;
        alea_primitive_type_t type;
        int nr;
        if (quads_left && group->ordinals[lane] == i) {
            if (!lane) {
                unsigned active = quads_left < RAY_QUADRIC_LANES ? quads_left : RAY_QUADRIC_LANES;
                ray_intersect_quadric_packet(ray, &group->packet, active, quad_hits, quad_roots);
            }
            type = ALEA_PRIMITIVE_QUADRIC;
            t = quad_roots[lane];
            nr = quad_hits[lane++];
            quads_left--;
            if (lane == RAY_QUADRIC_LANES) { group++; lane = 0; }
        } else {
            uint32_t id=sys->surfaces.data[indices[i]].primitive_id;
            type=sys->primitives.data[id].type;
            if (type>=ALEA_PRIMITIVE_TORUS_X && type<=ALEA_PRIMITIVE_TORUS_Z) {
                const alea_torus_data_t* torus = alea_primitive_payload_const(sys, id);
                if (ray_torus_window_may_hit(ray, torus, lo, hi))
                    exceptional[ne++]=(uint32_t)i;
                continue;
            }
            alea_primitive_data_t data;
            if (!raycast_primitive_copy_payload(sys,id,type,&data)) { free(roots); free(exceptional); return; }
            nr=type==ALEA_PRIMITIVE_PLANE ?
                ray_intersect_plane(ray,&data.plane,single,NULL,NULL,NULL) :
                ray_intersect_primitive(ray,type,&data,single);
        }
        result->surfaces_tested++;
        RAY_PROFILE_COUNT(recovery_root_tests);
        RAY_PROFILE_COUNT(root_tests_by_type[type]);
        for (int j=0;j<nr;j++) if (t[j]>=lo && t[j]<=hi) {
            if (n==capacity) {
                size_t next=capacity*2;
                size_t maximum=(available-exceptional_bytes)/sizeof(*roots);
                if (next>maximum) next=maximum;
                if (next<=capacity) { free(roots); free(exceptional); return; }
                recovery_root_t* grown=realloc(roots,next*sizeof(*roots));
                if (!grown) { free(roots); free(exceptional); return; }
                roots=grown; capacity=next; root_bytes=next*sizeof(*roots);
            }
            roots[n++]=recovery_root_make(t[j], (uint32_t)i, j);
            RAY_PROFILE_COUNT(roots_retained_by_type[type]);
        }
    }
    size_t indexed_bytes=0;
    recovery_root_t* indexed = recovery_bin_roots(roots,n,lo,hi,
        available-exceptional_bytes-root_bytes,&indexed_bytes);
    if (indexed) {
        free(roots);
        roots=indexed;
        root_bytes=indexed_bytes;
        entry->roots_binned=true;
        RAY_PROFILE_COUNT(indexed_root_lists);
    } else if (n) {
        qsort(roots,n,sizeof(*roots),recovery_root_compare);
        RAY_PROFILE_ADD(roots_sorted,n);
        recovery_root_t* compact=realloc(roots,n*sizeof(*roots));
        if (compact) { roots=compact; root_bytes=n*sizeof(*roots); }
    } else {
        free(roots);
        roots=NULL;
        root_bytes=0;
    }
    entry->roots=roots;
    entry->root_count=n;
    entry->exceptional=exceptional;
    entry->exceptional_count=ne;
    entry->roots_lo=lo;
    entry->roots_hi=hi;
    entry->roots_complete=true;
    table->root_bytes+=root_bytes+exceptional_bytes;
#ifdef ALEA_RAY_PROFILE
    if (alea_ray_profile.enabled && table->root_bytes > alea_ray_profile.peak_root_bytes)
        alea_ray_profile.peak_root_bytes = table->root_bytes;
#endif
    RAY_PROFILE_COUNT(recovery_root_builds);
}

/* Keep the larger key and table work off the ordinary indexed scan's stack. */
#if defined(__GNUC__) || defined(__clang__)
__attribute__((noinline))
#elif defined(_MSC_VER)
__declspec(noinline)
#endif
static void raycast_indexed_surface_hits_cached(alea_system_t* sys,
    const alea_ray_t* ray, const uint32_t* indices, size_t count,
    double t_min, double t_max, alea_raycast_result_t* result,
    const alea_matrix_t* local_to_world, const recovery_plane_index_t* plane_index) {
    if (alea_ray_interrupted(result)) { result->breakpoint_failed = 1; return; }
    if (result->hits.count) t_max = fmin(t_max, result->hits.data[0].t);
    recovery_context_key_t key;
    recovery_context_entry_t* entry = NULL;
    if (isfinite(t_min) && isfinite(t_max) && t_min >= 0 && t_max >= t_min &&
        recovery_context_key(&key, indices, count, ray, local_to_world)) {
#ifdef ALEA_RAY_PROFILE
        int previous=alea_ray_profile_enter(RAY_RECOVERY_CACHE);
#endif
        entry = recovery_context_find(result->recovery_context, &key);
#ifdef ALEA_RAY_PROFILE
        alea_ray_profile_leave(previous);
#endif
        if (recovery_context_contains(entry, t_min, t_max)) {
            RAY_PROFILE_COUNT(recovery_cache_hits);
            return;
        }
        /* An indexed universe can filter the remaining ray on its first visit.
         * Other lists keep the second-miss threshold. */
        if (entry && !entry->roots_attempted && (plane_index || entry->complete)) {
            recovery_roots_build(sys, ray, indices, count, result->recovery_context,
                                 entry, result, plane_index);
            if (entry->roots_complete && !entry->complete) {
                entry->key = key;
                entry->complete = true;
                entry->lo = INFINITY;
                entry->hi = -INFINITY;
            }
        }
        if (entry && entry->roots_complete &&
            t_min >= entry->roots_lo && t_max <= entry->roots_hi) {
            RAY_PROFILE_COUNT(recovery_root_queries);
            const recovery_root_t* root=recovery_root_first(entry,t_min,t_max,result->breakpoint_t_min);
            if (root) { RAY_PROFILE_COUNT(cached_root_candidates); }
            else if (!entry->exceptional_count) { RAY_PROFILE_COUNT(cached_empty_windows); }
            recovery_root_t selected;
            if (root) { selected=*root; root=&selected; }
            /* Merge exceptional roots by original card/root order, retaining
             * the same window test even if a cached root is nearer. */
            for (size_t k=0;k<entry->exceptional_count;k++) {
                uint32_t ordinal=entry->exceptional[k];
                uint32_t id=sys->surfaces.data[indices[ordinal]].primitive_id;
                alea_primitive_data_t data;
                alea_primitive_type_t type=sys->primitives.data[id].type;
                if (!raycast_primitive_copy_payload(sys,id,type,&data)) { result->breakpoint_failed=1; return; }
                if (!ray_torus_window_may_hit(ray,&data.torus,t_min,t_max)) continue;
                double t[4];
                int nr=ray_intersect_primitive(ray,type,&data,t);
                result->surfaces_tested++;
                for (int j=0;j<nr;j++) if (t[j]>=t_min && t[j]<=t_max &&
                        t[j]>nextafter(result->breakpoint_t_min,INFINITY)) {
                    recovery_root_t candidate = recovery_root_make(t[j], ordinal, j);
                    if (!root || recovery_root_compare(&candidate,root)<0) {
                        selected=candidate;
                        root=&selected;
                    }
                }
            }
            if (root) {
                const alea_surface_entry_t* surf=&sys->surfaces.data[indices[root->order / 4u]];
                alea_primitive_type_t type=sys->primitives.data[surf->primitive_id].type;
                alea_primitive_data_t data;
                if (!raycast_primitive_copy_payload(sys,surf->primitive_id,type,&data)) {
                    result->breakpoint_failed=1;
                    return;
                }
                alea_ray_hit_t hit={.t=recovery_root_distance(root),.surface_id=surf->mc_surface_id,
                                    .primitive_id=surf->primitive_id};
                double x,y,z;
                alea_ray_point_at(ray,hit.t,&x,&y,&z);
                primitive_normal_at(type,&data,x,y,z,&hit.nx,&hit.ny,&hit.nz);
                if (local_to_world) boundary_event_world_normal(local_to_world,
                    hit.nx,hit.ny,hit.nz,&hit.nx,&hit.ny,&hit.nz);
                if (add_hit(result,&hit)) return;
            }
            entry->lo=t_min;
            entry->hi=t_max;
            return;
        }
    }
    raycast_indexed_surface_scan(sys, ray, indices, count, t_min, t_max,
                                result, local_to_world, NULL, plane_index);
    if (entry && !result->breakpoint_failed && !alea_ray_interrupted(result)) {
        entry->key = key;
        entry->lo = t_min;
        entry->hi = t_max;
        entry->complete = true;
    }
}

static void raycast_indexed_surface_hits(alea_system_t* sys,
    const alea_ray_t* ray, const uint32_t* indices, size_t count,
    double t_min, double t_max, alea_raycast_result_t* result,
    const alea_matrix_t* local_to_world, uint8_t* seen_surfaces,
    const recovery_plane_index_t* plane_index) {
    /* Short recoveries should not pay for a large optional table. Attempt
     * allocation once, only after enough scans to amortize its initialization. */
    if (result->breakpoint_nearest_only && result->recovery_context_warmup &&
        --result->recovery_context_warmup == 0)
        result->recovery_context = recovery_context_create(2048);
    if (result->breakpoint_nearest_only && result->recovery_context && !seen_surfaces) {
        raycast_indexed_surface_hits_cached(sys, ray, indices, count, t_min, t_max,
                                           result, local_to_world, plane_index);
    } else {
        raycast_indexed_surface_scan(sys, ray, indices, count, t_min, t_max,
                                    result, local_to_world, seen_surfaces, plane_index);
    }
}

static void raycast_cell_indexed_surface_hits(alea_system_t* sys,
                                              const alea_ray_t* ray,
                                              const alea_cell_entry_t* cell,
                                              double t_min,
                                              double t_max,
                                              alea_raycast_result_t* result,
                                              const alea_matrix_t* local_to_world,
                                              uint8_t* seen_surfaces) {
    if (!cell->surface_indices || !cell->surface_index_count) {
        alea_raycast_tree_primitives(sys, ray, cell->root_node_id,
                                t_min, t_max, result, local_to_world);
        return;
    }
    raycast_indexed_surface_hits(sys, ray, cell->surface_indices,
        cell->surface_index_count, t_min, t_max, result, local_to_world, seen_surfaces, NULL);
}

/**
 * Raycast all surfaces in a universe using a translated ray.
 * The ray is in element-local coordinates (origin shifted by -element_center).
 * Since translation preserves t, hits can be added directly to the result.
 */
static void raycast_universe_surfaces(alea_system_t* sys,
                                      const alea_ray_t* local_ray,
                                      int universe_id,
                                      double t_min, double t_max,
                                      alea_raycast_result_t* result,
                                      const alea_matrix_t* local_to_world) {
    if (result->breakpoint_hit_limit && result->breakpoint_failed) return;
    if (result->breakpoint_nearest_only && result->hits.count)
        t_max = fmin(t_max, result->hits.data[0].t);
    const alea_universe_t* univ = alea_get_universe(sys, universe_id);
    if (!univ) {
        if (result->breakpoint_hit_limit) result->breakpoint_failed = 1;
        return;
    }

    /* Readers that prepare only a subset of caches must not observe a list
     * while another thread is still publishing the optional recovery cache. */
    if ((result->breakpoint_nearest_only || result->transport_fill_filter) &&
        alea_system_query_cache_ready(sys, ALEA_CACHE_RECOVERY_SURFACES) &&
        univ->recovery_surfaces &&
        (result->breakpoint_nearest_only || univ->recovery_surfaces_exact)) {
        size_t first = result->hits.count;
        raycast_indexed_surface_hits(sys, local_ray, univ->recovery_surfaces,
            univ->recovery_surface_count, t_min, t_max, result, local_to_world, NULL,
            univ->recovery_planes);
        if (result->breakpoint_nearest_only || result->breakpoint_failed ||
            result->hits.count == first) return;
        /* A card hit may represent several tree leaves. Replay this occurrence
         * to preserve all raw hits, including their order and budget cost. */
        result->hits.count = first;
    }

    /* Optional scratch is bounded and scoped to this occurrence, so distinct
     * transforms and lattice elements never share solved-surface state.
     * Allocation failure retains the original scan. */
    const size_t seen_bytes = sys->surfaces.count / 8 +
                              (sys->surfaces.count % 8 != 0);
    uint8_t* seen_surfaces = result->breakpoint_nearest_only &&
        univ->cell_indices.count > 1 && seen_bytes <= 1024 * 1024
        ? calloc(seen_bytes, 1) : NULL;
    for (size_t c = 0; c < univ->cell_indices.count; c++) {
        if (alea_ray_interrupted(result)) { result->breakpoint_failed = 1; break; }
        if (result->breakpoint_hit_limit && result->breakpoint_failed) break;
        size_t cell_idx = univ->cell_indices.data[c];
        if (cell_idx >= alea_vec_count(&sys->cells)) {
            if (result->breakpoint_hit_limit) result->breakpoint_failed = 1;
            continue;
        }
        const alea_cell_entry_t* cell = &sys->cells.data[cell_idx];
        if (result->breakpoint_nearest_only) {
            /* Prepared surface lists retain first tree-occurrence order,
             * including distinct cards sharing a primitive. Nearest-only
             * recovery can skip repeated references without changing ties. */
            raycast_cell_indexed_surface_hits(sys, local_ray, cell,
                t_min, t_max, result, local_to_world, seen_surfaces);
        } else {
            alea_raycast_tree_primitives(sys, local_ray, cell->root_node_id,
                                    t_min, t_max, result, local_to_world);
        }
    }
    free(seen_surfaces);
}

static bool system_has_lattice_cells(const alea_system_t* sys);
static bool system_has_fill_cells(const alea_system_t* sys);
static int dedup_sorted_hits(alea_raycast_result_t* result);
static void raycast_lattice_rect(alea_system_t* sys,
                                 const alea_ray_t* ray,
                                 const alea_cell_entry_t* lat_cell,
                                 double t_min, double t_max, int depth,
                                 alea_raycast_result_t* result,
                                 const alea_matrix_t* frame_to_world);
static void raycast_lattice_hex(alea_system_t* sys,
                                const alea_ray_t* ray,
                                const alea_cell_entry_t* lat_cell,
                                double t_min, double t_max, int depth,
                                alea_raycast_result_t* result,
                                const alea_matrix_t* frame_to_world);
static void raycast_lattice_universe_hits_recursive(
    alea_system_t* sys, const alea_ray_t* ray, int universe_id,
    double t_min, double t_max, int depth, alea_raycast_result_t* result,
    const alea_matrix_t* frame_to_world);

static int raycast_root_blas_surfaces(alea_system_t* sys,
                                      const alea_ray_t* ray,
                                      double t_min,
                                      double t_max,
                                      alea_raycast_result_t* result) {
    alea_raycast_result_clear(result);
    result->ray = *ray;

    const alea_hier_spatial_stats_t* stats =
        alea_hier_spatial_index_stats(sys->hier_spatial_index);
    if (!stats || stats->placement_count == 0 ||
        stats->max_universe_cells <= 0) {
        return 0;
    }

    alea_hier_placement_ray_candidate_t* placements =
        malloc(stats->placement_count * sizeof(*placements));
    alea_hier_ray_candidate_t* candidates =
        malloc((size_t)stats->max_universe_cells * sizeof(*candidates));
    if (!placements || !candidates) {
        free(placements);
        free(candidates);
        return -1;
    }

    int placement_count = alea_hier_spatial_query_placements_ray(
        sys,
        ray,
        t_min, t_max,
        placements, stats->placement_count);
    if (placement_count < 0) {
        free(placements);
        free(candidates);
        return -1;
    }
    result->blas_placement_candidates = (size_t)placement_count;

    for (int p = 0; p < placement_count; p++) {
        alea_hier_placement_ray_candidate_t* placement = &placements[p];
        const alea_universe_t* univ =
            alea_get_universe(sys, placement->universe_id);
        if (!univ || univ->cell_indices.count == 0) continue;

        alea_ray_t local_ray;
        alea_ray_transform_inverse(&placement->transform, ray, &local_ray);

        size_t max_candidates = univ->cell_indices.count;
        if (max_candidates > (size_t)stats->max_universe_cells) {
            max_candidates = (size_t)stats->max_universe_cells;
        }

        double placement_min = fmax(t_min, placement->t_enter - RAY_EPSILON);
        double placement_max = fmin(t_max, placement->t_exit + RAY_EPSILON);
        bool lattice_placement = false;
        if ((size_t)placement->parent_cell_index < alea_vec_count(&sys->cells)) {
            const alea_cell_entry_t* parent_cell =
                &sys->cells.data[placement->parent_cell_index];
            lattice_placement = parent_cell->lat_type != 0 && parent_cell->lat_fill;
        }
        if (placement->placement_index != 0 && !lattice_placement) {
            double samples[3] = {
                placement_min + RAY_EPSILON,
                0.5 * (placement_min + placement_max),
                placement_max - RAY_EPSILON
            };
            int chain_valid = 0;
            for (int s = 0; s < 3; s++) {
                if (samples[s] < placement_min) samples[s] = placement_min;
                if (samples[s] > placement_max) samples[s] = placement_max;
                double sx, sy, sz;
                alea_ray_point_at(ray, samples[s], &sx, &sy, &sz);
                int valid = alea_hier_spatial_check_placement_chain(
                    sys, placement->placement_index, sx, sy, sz);
                if (valid < 0) {
                    free(placements);
                    free(candidates);
                    return -1;
                }
                if (valid) {
                    chain_valid = 1;
                    break;
                }
            }
            if (!chain_valid) {
                result->blas_placements_pruned++;
                continue;
            }
        }

        result->blas_universe_queries++;
        int count = alea_hier_spatial_query_universe_ray(
            sys, placement->universe_id,
            &local_ray,
            placement_min, placement_max,
            candidates, max_candidates);
        if (count < 0) {
            free(placements);
            free(candidates);
            return -1;
        }
        result->blas_cell_candidates += (size_t)count;

        for (int i = 0; i < count; i++) {
            uint32_t cell_index = candidates[i].cell_index;
            if ((size_t)cell_index >= alea_vec_count(&sys->cells)) continue;
            const alea_cell_entry_t* cell = &sys->cells.data[cell_index];
            if (cell->root_node_id == ALEA_NODE_ID_INVALID) continue;
            result->blas_cells_tested++;

            double local_min = fmax(placement_min,
                                    candidates[i].t_enter - RAY_EPSILON);
            double local_max = fmin(placement_max,
                                    candidates[i].t_exit + RAY_EPSILON);
            raycast_cell_indexed_surface_hits(
                sys, &local_ray, cell, local_min, local_max, result,
                &placement->transform, NULL);
            if (cell->lat_type == 1 && cell->lat_fill) {
                raycast_lattice_rect(sys, &local_ray, cell,
                                     local_min, local_max, 0, result,
                                     &placement->transform);
            } else if (cell->lat_type == 2 && cell->lat_fill) {
                raycast_lattice_hex(sys, &local_ray, cell,
                                    local_min, local_max, 0, result,
                                    &placement->transform);
            }
        }
    }

    free(placements);
    free(candidates);
    result->blas_hits_before_dedup = result->hits.count;
    dedup_sorted_hits(result);
    return 0;
}

void alea_ray_transform_inverse(const alea_matrix_t* mat,
                                const alea_ray_t* ray,
                                alea_ray_t* local_ray) {
    double ox = ray->ox;
    double oy = ray->oy;
    double oz = ray->oz;
    alea_matrix_transform_point_inverse(mat, &ox, &oy, &oz);

    double dx = mat->inv[0] * ray->dx + mat->inv[1] * ray->dy + mat->inv[2] * ray->dz;
    double dy = mat->inv[4] * ray->dx + mat->inv[5] * ray->dy + mat->inv[6] * ray->dz;
    double dz = mat->inv[8] * ray->dx + mat->inv[9] * ray->dy + mat->inv[10] * ray->dz;

    /* Keep O_local + t D_local at the inverse image of the world point.
     * Rounded rotation cards can change |D_local|. All walker limits, roots
     * and lattice steps use the original world-distance parameter t. */
    ray_init_parametric(local_ray, ox, oy, oz, dx, dy, dz);
}

static void raycast_lattice_rect(alea_system_t* sys,
                                 const alea_ray_t* ray,
                                 const alea_cell_entry_t* lat_cell,
                                 double t_min, double t_max, int depth,
                                 alea_raycast_result_t* result,
                                 const alea_matrix_t* frame_to_world);
static void raycast_lattice_hex(alea_system_t* sys,
                                const alea_ray_t* ray,
                                const alea_cell_entry_t* lat_cell,
                                double t_min, double t_max, int depth,
                                alea_raycast_result_t* result,
                                const alea_matrix_t* frame_to_world);

static void raycast_fill_universe_hits_recursive(alea_system_t* sys,
                                                 const alea_ray_t* global_ray,
                                                 int universe_id,
                                                 const alea_matrix_t* accumulated,
                                                 double t_min, double t_max,
                                                 int depth,
                                                 alea_raycast_result_t* result) {
    if (!sys || !global_ray || !accumulated || !result)
        return;
    if (result->breakpoint_hit_limit && result->breakpoint_failed) return;
    if (result->breakpoint_nearest_only && result->hits.count)
        t_max = fmin(t_max, result->hits.data[0].t);
    if (depth >= MAX_FILL_RAYCAST_DEPTH) {
        if (result->breakpoint_hit_limit) result->breakpoint_failed = 1;
        return;
    }

    const alea_universe_t* univ = alea_get_universe(sys, universe_id);
    if (!univ) {
        if (result->breakpoint_hit_limit) result->breakpoint_failed = 1;
        return;
    }

    bool use_children = result->breakpoint_nearest_only &&
        alea_system_query_cache_ready(sys, ALEA_CACHE_RECOVERY_CHILDREN) &&
        univ->recovery_children_ready;
    size_t count = use_children ? univ->recovery_child_count : univ->cell_indices.count;
    const size_t* indices = use_children ? univ->recovery_child_cells : univ->cell_indices.data;
    if (count == 0) return;

    alea_ray_t parent_local_ray;
    alea_ray_transform_inverse(accumulated, global_ray, &parent_local_ray);

    /* Siblings can place the same universe in exactly the same frame. Once a
     * traversal leaves an existing nearest hit untouched, repeating it with
     * that same limit cannot improve the hit or change its tie provenance. */
    int previous_universe = 0, previous_transform = 0;
    double previous_limit = -1.0;
    const bool nearest = result->breakpoint_nearest_only && result->breakpoint_hit_limit;
    const bool prepared_matrices = nearest &&
        alea_system_query_cache_ready(sys, ALEA_CACHE_HIER_SPATIAL);

    for (size_t i = 0; i < count; i++) {
        if ((i & 1023) == 0 && alea_ray_interrupted(result)) {
            result->breakpoint_failed = 1;
            return;
        }
        if (result->breakpoint_hit_limit && result->breakpoint_failed) return;
        size_t cell_idx = indices[i];
        if (cell_idx >= alea_vec_count(&sys->cells)) {
            if (result->breakpoint_hit_limit) result->breakpoint_failed = 1;
            continue;
        }

        const alea_cell_entry_t* cell = &sys->cells.data[cell_idx];

        /* Global diagnostic breakpoints must use the same occurrence frame as
         * selected traversal.  A lattice nested below a transformed fill is
         * expressed in this universe's local coordinates, not world space. */
        if (cell->lat_type == 1 && cell->lat_fill) {
            raycast_lattice_rect(sys, &parent_local_ray, cell,
                                 t_min, t_max, depth, result, accumulated);
        } else if (cell->lat_type == 2 && cell->lat_fill) {
            raycast_lattice_hex(sys, &parent_local_ray, cell,
                                t_min, t_max, depth, result, accumulated);
        }

        if (cell->fill_universe <= 0)
            continue;

        RAY_PROFILE_COUNT(fill_children);
        const double limit = nearest && result->hits.count ? result->hits.data[0].t : -1.0;
        if (limit >= 0 && limit == previous_limit &&
            cell->fill_universe == previous_universe && cell->fill_transform == previous_transform) {
            RAY_PROFILE_COUNT(fill_duplicates);
            continue;
        }

        if (cell->root_node_id != ALEA_NODE_ID_INVALID &&
            cell->root_node_id < alea_vec_count(&sys->nodes)) {
            const alea_bbox_t bbox_v = alea_node_bbox_get(&sys->nodes.data[cell->root_node_id].bbox);
            if (!ray_bbox_slab(&parent_local_ray, &bbox_v, t_min, t_max)) {
                RAY_PROFILE_COUNT(fill_bbox_rejected);
                continue;
            }
        }

        alea_matrix_t fill_transform;
        const alea_matrix_t* matrix = prepared_matrices && cell_idx <= UINT32_MAX
            ? alea_cell_fill_matrix(sys, (uint32_t)cell_idx) : NULL;
        if (!matrix && cell->fill_transform > 0) {
            const alea_transform_t* tr = alea_get_transform(sys, cell->fill_transform);
            if (!tr ||
                !alea_matrix_from_mcnp(&fill_transform, tr->cosines,
                                       tr->value_count, false)) {
                if (result->breakpoint_hit_limit) result->breakpoint_failed = 1;
                continue;
            }
        } else if (!matrix) {
            alea_matrix_identity(&fill_transform);
        }

        alea_matrix_t fill_accumulated;
        alea_matrix_multiply(&fill_accumulated, accumulated, matrix ? matrix : &fill_transform);
        if (!fill_accumulated.has_inverse && !alea_matrix_invert(&fill_accumulated)) {
            if (result->breakpoint_hit_limit) result->breakpoint_failed = 1;
            continue;
        }

        alea_ray_t fill_local_ray;
        alea_ray_transform_inverse(&fill_accumulated, global_ray, &fill_local_ray);

        RAY_PROFILE_COUNT(fill_descents);
        raycast_universe_surfaces(
            sys, &fill_local_ray, cell->fill_universe,
            t_min, t_max, result, &fill_accumulated);
        raycast_fill_universe_hits_recursive(sys, global_ray, cell->fill_universe,
                                             &fill_accumulated, t_min, t_max,
                                             depth + 1, result);
        if (limit >= 0 && !result->breakpoint_failed && result->hits.data[0].t == limit) {
            previous_universe = cell->fill_universe;
            previous_transform = cell->fill_transform;
            previous_limit = limit;
        }
    }
}

static void raycast_add_fill_hits(alea_system_t* sys,
                                  const alea_ray_t* ray,
                                  double t_min, double t_max,
                                  alea_raycast_result_t* result) {
    alea_matrix_t identity;
    alea_matrix_identity(&identity);
    raycast_fill_universe_hits_recursive(sys, ray, 0, &identity,
                                         t_min, t_max, 0, result);
}

/* ray_bbox_slab_enter_exit in ray_bbox.h replaces ray_bbox_slab_enter_exit */

static int lattice_rect_bounds(const alea_cell_entry_t* cell,
                               alea_bbox_t* out_bounds) {
    if (!cell || !out_bounds || cell->lat_pitch[0] <= 0.0 ||
        cell->lat_pitch[1] <= 0.0) {
        return -1;
    }
    const int ni = cell->lat_fill_dims[1] - cell->lat_fill_dims[0] + 1;
    const int nj = cell->lat_fill_dims[3] - cell->lat_fill_dims[2] + 1;
    const int nk = cell->lat_fill_dims[5] - cell->lat_fill_dims[4] + 1;
    if (ni <= 0 || nj <= 0 || nk <= 0) return -1;
    *out_bounds = (alea_bbox_t){
        .min_x = cell->lat_lower_left[0],
        .max_x = cell->lat_lower_left[0] + ni * cell->lat_pitch[0],
        .min_y = cell->lat_lower_left[1],
        .max_y = cell->lat_lower_left[1] + nj * cell->lat_pitch[1],
        .min_z = nk == 1 && cell->lat_pitch[2] <= 0.0 ? -1e30
                                                      : cell->lat_lower_left[2],
        .max_z = nk == 1 && cell->lat_pitch[2] <= 0.0 ? 1e30
            : cell->lat_lower_left[2] + nk * cell->lat_pitch[2]
    };
    return 0;
}

static int lattice_hex_bounds(const alea_cell_entry_t* cell,
                              alea_bbox_t* out_bounds) {
    if (!cell || !out_bounds || cell->lat_pitch[0] <= 0.0) return -1;
    const int imin = cell->lat_fill_dims[0], imax = cell->lat_fill_dims[1];
    const int jmin = cell->lat_fill_dims[2], jmax = cell->lat_fill_dims[3];
    const int nk = cell->lat_fill_dims[5] - cell->lat_fill_dims[4] + 1;
    if (imax < imin || jmax < jmin || nk <= 0) return -1;

    double min_x = DBL_MAX, max_x = -DBL_MAX;
    double min_y = DBL_MAX, max_y = -DBL_MAX;
    const double pitch = cell->lat_pitch[0];
    const double radius = pitch / M_SQRT3;
    for (int i = imin; i <= imax; i++) {
        for (int j = jmin; j <= jmax; j++) {
            const double x = i * pitch + j * pitch * 0.5;
            const double y = j * pitch * M_SQRT3 * 0.5;
            if (x - radius < min_x) min_x = x - radius;
            if (x + radius > max_x) max_x = x + radius;
            if (y - pitch * 0.5 < min_y) min_y = y - pitch * 0.5;
            if (y + pitch * 0.5 > max_y) max_y = y + pitch * 0.5;
        }
    }
    *out_bounds = (alea_bbox_t){
        .min_x = min_x, .max_x = max_x,
        .min_y = min_y, .max_y = max_y,
        .min_z = nk > 1 ? cell->lat_lower_left[2] : -1e30,
        .max_z = nk > 1 ? cell->lat_lower_left[2] + nk * cell->lat_pitch[2]
                         : 1e30
    };
    return 0;
}

/* The geometry-specific DDA walkers only decide the next interval and its
 * canonical lattice location.  Publishing the synthetic transition and
 * tracing the selected element are shared so full-trace and future early-stop
 * steppers cannot acquire different element-local ray conventions. */
static int raycast_lattice_element_step(
    alea_system_t* sys,
    const alea_ray_t* ray,
    const alea_lattice_location_t* location,
    double t_enter,
    double t_exit,
    bool emit_synthetic_transition,
    int depth,
    alea_raycast_result_t* result,
    const alea_matrix_t* frame_to_world) {
    if (!sys || !ray || !result) return -1;
    if (emit_synthetic_transition) {
        alea_ray_hit_t boundary = {
            .t = t_enter,
            .surface_id = 0,
            .primitive_id = ALEA_PRIMITIVE_ID_INVALID
        };
        if (add_hit(result, &boundary) != 0) {
            if (!result->breakpoint_hit_limit)
                ALEA_LOG_WARN("add_hit failed (out of memory) - lattice raycast incomplete");
            return -1;
        }
    }
    if (!location || t_exit <= t_enter) return 0;

    alea_ray_t local_ray = *ray;
    local_ray.ox -= location->ox;
    local_ray.oy -= location->oy;
    local_ray.oz -= location->oz;
    raycast_universe_surfaces(sys, &local_ray, location->fill_universe,
                              t_enter, t_exit, result, frame_to_world);
    raycast_lattice_universe_hits_recursive(
        sys, &local_ray, location->fill_universe, t_enter, t_exit,
        depth + 1, result, frame_to_world);
    return result->breakpoint_hit_limit && result->breakpoint_failed ? -1 : 0;
}

static bool raycast_has_hit_at(const alea_raycast_result_t* result, double t) {
    if (!result) return false;
    for (size_t i = 0; i < result->hits.count; i++) {
        /* Parsed lattice pitch/bounds can differ slightly from the source CSG
         * plane used to infer them.  Treat a nearby physical hit as the same
         * ownership boundary so the DDA does not publish a micro-segment. */
        if (result->hits.data[i].surface_id > 0 &&
            fabs(result->hits.data[i].t - t) <= SURFACE_SAMPLE_OFFSET)
            return true;
    }
    return false;
}

int alea_lattice_raycast_interval(const alea_ray_t* ray,
                                    const alea_cell_entry_t* cell,
                                    double t_min, double t_max,
                                    double* out_enter, double* out_exit) {
    if (!ray || !cell || !out_enter || !out_exit) return -1;
    *out_enter = t_min;
    *out_exit = t_max;

    alea_bbox_t bounds;
    if (cell->lat_type == 1) {
        if (lattice_rect_bounds(cell, &bounds) != 0) return -1;
        if (cell->lat_fill_repeating || cell->lat_outer_universe > 0) return 0;
    } else if (cell->lat_type == 2) {
        if (lattice_hex_bounds(cell, &bounds) != 0) return -1;
        if (cell->lat_outer_universe > 0) return 0;
    } else {
        return -1;
    }
    return ray_bbox_slab_enter_exit(ray, &bounds, t_min, t_max,
                                    out_enter, out_exit) ? 0 : 1;
}

int alea_lattice_raycast_step_limit(const alea_ray_t* ray,
                                      const alea_cell_entry_t* cell,
                                      double t_enter, double t_exit) {
    const int ni = cell->lat_fill_dims[1] - cell->lat_fill_dims[0] + 1;
    const int nj = cell->lat_fill_dims[3] - cell->lat_fill_dims[2] + 1;
    const int nk = cell->lat_fill_dims[5] - cell->lat_fill_dims[4] + 1;
    if (ni <= 0 || nj <= 0 || nk <= 0) return 0;
    const int unbounded = cell->lat_fill_repeating ||
                          cell->lat_outer_universe > 0;
    if (cell->lat_type == 1 && unbounded) {
        const double span = t_exit - t_enter;
        const double crossings = fabs(ray->dx) * span / cell->lat_pitch[0]
            + fabs(ray->dy) * span / cell->lat_pitch[1]
            + ((cell->lat_pitch[2] > 0.0 &&
                (nk > 1 || cell->lat_fill_repeating))
                ? fabs(ray->dz) * span / cell->lat_pitch[2] : 0.0);
        return crossings < (double)INT_MAX - 8.0
            ? (int)ceil(crossings) + 8 : INT_MAX;
    }
    if (cell->lat_type == 2 && unbounded) {
        const double span = t_exit - t_enter;
        const double p = cell->lat_pitch[0];
        if (p <= 0.0) return 0;
        const double dq = fabs(ray->dx / p - ray->dy / (p * M_SQRT3));
        const double dr = fabs(2.0 * ray->dy / (p * M_SQRT3));
        const double ds = fabs(-ray->dx / p - ray->dy / (p * M_SQRT3));
        const double dz = (nk > 1 && cell->lat_pitch[2] > 0.0)
            ? fabs(ray->dz) / cell->lat_pitch[2] : 0.0;
        const double crossings = (dq + dr + ds + dz) * span;
        return crossings < (double)INT_MAX - 12.0
            ? (int)ceil(crossings) + 12 : INT_MAX;
    }
    return cell->lat_type == 1 ? 2 * (ni + nj + nk) + 4
                               : 3 * (ni + nj + nk) + 10;
}

/* Full global tracing uses the same boundary query as the cell-aware
 * first-visible walker.  Geometry-specific rectangular/hex math is confined
 * to alea_lattice_next_boundary(); this loop owns interval construction, canonical
 * location lookup, synthetic event publication, and element-local tracing. */
static void raycast_lattice_walk(alea_system_t* sys,
                                 const alea_ray_t* ray,
                                 const alea_cell_entry_t* cell,
                                 double t_min, double t_max, int depth,
                                 alea_raycast_result_t* result,
                                 const alea_matrix_t* frame_to_world) {
    if (result->breakpoint_hit_limit && result->breakpoint_failed) return;
    if (result->breakpoint_nearest_only && result->hits.count)
        t_max = fmin(t_max, result->hits.data[0].t);
    double t_enter, t_exit;
    const int interval = alea_lattice_raycast_interval(ray, cell, t_min, t_max,
                                                  &t_enter, &t_exit);
    if (interval != 0 || t_exit <= t_enter) {
        if (interval < 0 && result->breakpoint_hit_limit)
            result->breakpoint_failed = 1;
        return;
    }

    const int max_steps = alea_lattice_raycast_step_limit(ray, cell, t_enter, t_exit);
    if (max_steps <= 0) {
        if (result->breakpoint_hit_limit) result->breakpoint_failed = 1;
        return;
    }
    alea_lattice_location_t previous = {0};
    int have_previous = 0;
    double t_current = t_enter;
    if (t_enter > (result->breakpoint_hit_limit ? t_min :
                   t_min + RAY_EPSILON) &&
        (result->breakpoint_hit_limit ||
         !raycast_has_hit_at(result, t_enter))) {
        if (raycast_lattice_element_step(sys, ray, NULL, t_enter, t_enter,
                                         true, depth, result,
                                         frame_to_world) != 0) {
            if (result->breakpoint_hit_limit) result->breakpoint_failed = 1;
            return;
        }
    }
    for (int step = 0; step < max_steps && t_current < t_exit; step++) {
        if (result->breakpoint_hit_limit && result->breakpoint_failed) return;
        if (result->breakpoint_nearest_only && result->hits.count) {
            t_exit = fmin(t_exit, result->hits.data[0].t);
            if (t_current >= t_exit) return;
        }
        double t_next = alea_lattice_next_boundary(ray, cell, t_current, t_exit);
        if (result->breakpoint_hit_limit &&
            t_next > t_current && t_next <= t_current + RAY_EPSILON &&
            t_next < t_exit) {
            result->breakpoint_failed = 1;
            return;
        }
        if (t_next <= t_current + RAY_EPSILON || t_next > t_exit)
            t_next = t_exit;
        const double t_sample = t_current +
            0.5 * (t_next - t_current);
        double x, y, z;
        alea_ray_point_at(ray, t_sample, &x, &y, &z);
        alea_lattice_location_t location;
        const int locate_status = alea_lattice_locate_point(
            sys, cell, x, y, z, &location);
        if (locate_status < 0 && result->breakpoint_hit_limit) {
            result->breakpoint_failed = 1;
            return;
        }
        const int has_location = locate_status == 1;

        const bool changed = !have_previous || !has_location ||
            location.i != previous.i || location.j != previous.j ||
            location.k != previous.k ||
            location.fill_universe != previous.fill_universe;
        const bool emit_synthetic = step > 0 &&
            t_current > (result->breakpoint_hit_limit ? t_enter :
                         t_enter + RAY_EPSILON) && changed;
        if (raycast_lattice_element_step(
                sys, ray, has_location ? &location : NULL,
                t_current, t_next, emit_synthetic, depth, result,
                frame_to_world) != 0) {
            if (result->breakpoint_hit_limit) result->breakpoint_failed = 1;
            return;
        }
        if (has_location) {
            previous = location;
            have_previous = 1;
        } else {
            have_previous = 0;
        }
        t_current = t_next;
    }
    if (t_current < t_exit && result->breakpoint_hit_limit)
        result->breakpoint_failed = 1;
    if (t_exit < (result->breakpoint_hit_limit ? t_max :
                  t_max - RAY_EPSILON) &&
        (result->breakpoint_hit_limit ||
         !raycast_has_hit_at(result, t_exit))) {
        if (raycast_lattice_element_step(sys, ray, NULL, t_exit, t_exit,
                                         true, depth, result,
                                         frame_to_world) != 0 &&
            result->breakpoint_hit_limit)
            result->breakpoint_failed = 1;
    }
}

/**
 * DDA raycast through a rectangular lattice.
 * Steps through elements along the ray, raycasting base universe surfaces
 * in each element's local coordinate system.
 */
static void raycast_lattice_rect(alea_system_t* sys,
                                 const alea_ray_t* ray,
                                 const alea_cell_entry_t* lat_cell,
                                 double t_min, double t_max, int depth,
                                 alea_raycast_result_t* result,
                                 const alea_matrix_t* frame_to_world) {
    raycast_lattice_walk(sys, ray, lat_cell, t_min, t_max, depth, result,
                         frame_to_world);
}

static void raycast_lattice_hex(alea_system_t* sys,
                                const alea_ray_t* ray,
                                const alea_cell_entry_t* lat_cell,
                                double t_min, double t_max, int depth,
                                alea_raycast_result_t* result,
                                const alea_matrix_t* frame_to_world) {
    raycast_lattice_walk(sys, ray, lat_cell, t_min, t_max, depth, result,
                         frame_to_world);
}

/* A lattice element can itself contain another lattice.  Descend with the
 * element-local ray so global diagnostic breakpoints retain the same nested
 * lattice coordinate convention as the selected walker. */
static void raycast_lattice_universe_hits_recursive(
    alea_system_t* sys, const alea_ray_t* ray, int universe_id,
    double t_min, double t_max, int depth, alea_raycast_result_t* result,
    const alea_matrix_t* frame_to_world) {
    if (!sys || !ray || !result) return;
    if (result->breakpoint_hit_limit && result->breakpoint_failed) return;
    if (result->breakpoint_nearest_only && result->hits.count)
        t_max = fmin(t_max, result->hits.data[0].t);
    if (depth >= MAX_FILL_RAYCAST_DEPTH) {
        if (result->breakpoint_hit_limit) result->breakpoint_failed = 1;
        return;
    }
    const alea_universe_t* universe = alea_get_universe(sys, universe_id);
    if (!universe) {
        if (result->breakpoint_hit_limit) result->breakpoint_failed = 1;
        return;
    }
    for (size_t i = 0; i < universe->cell_indices.count; i++) {
        if (result->breakpoint_hit_limit && result->breakpoint_failed) return;
        const size_t cell_index = universe->cell_indices.data[i];
        if (cell_index >= alea_vec_count(&sys->cells)) {
            if (result->breakpoint_hit_limit) result->breakpoint_failed = 1;
            continue;
        }
        const alea_cell_entry_t* cell = &sys->cells.data[cell_index];
        if (cell->lat_type == 1 && cell->lat_fill) {
            raycast_lattice_rect(sys, ray, cell, t_min, t_max, depth, result,
                                 frame_to_world);
        } else if (cell->lat_type == 2 && cell->lat_fill) {
            raycast_lattice_hex(sys, ray, cell, t_min, t_max, depth, result,
                                frame_to_world);
        }
    }
}
/**
 * Add lattice surface hits for all lattice cells the ray may cross.
 * Called after alea_raycast_surfaces() and before sorting/dedup.
 */
static void raycast_add_lattice_hits(alea_system_t* sys,
                                     const alea_ray_t* ray,
                                     double t_min, double t_max,
                                     alea_raycast_result_t* result) {
    if (result->breakpoint_hit_limit && result->breakpoint_failed) return;
    if (result->breakpoint_nearest_only && result->hits.count)
        t_max = fmin(t_max, result->hits.data[0].t);
    alea_matrix_t identity;
    alea_matrix_identity(&identity);
    const size_t* indices = alea_system_query_cache_ready(sys, ALEA_CACHE_UNIVERSE)
                          ? sys->lattice_cells : NULL;
    const size_t count = indices ? sys->lattice_cell_count : sys->cells.count;
    for (size_t i = 0; i < count; i++) {
        if (result->breakpoint_hit_limit && result->breakpoint_failed) return;
        const alea_cell_entry_t* cell = &sys->cells.data[indices ? indices[i] : i];
        if (cell->lat_type == 0 || !cell->lat_fill) continue;

        if (cell->lat_type == 1) {
            raycast_lattice_rect(sys, ray, cell, t_min, t_max, 0, result,
                                 &identity);
        } else if (cell->lat_type == 2) {
            raycast_lattice_hex(sys, ray, cell, t_min, t_max, 0, result,
                                &identity);
        }
    }
}

static bool system_has_lattice_cells(const alea_system_t* sys) {
    if (!sys) return false;
    if (alea_system_query_cache_ready(sys, ALEA_CACHE_UNIVERSE)) return sys->has_populated_lattice;
    for (size_t i = 0; i < alea_vec_count(&sys->cells); i++) {
        const alea_cell_entry_t* cell = &sys->cells.data[i];
        if (cell->lat_type != 0 && cell->lat_fill && cell->lat_fill_count > 0)
            return true;
    }
    return false;
}

static bool system_has_fill_cells(const alea_system_t* sys) {
    if (!sys) return false;
    if (alea_system_query_cache_ready(sys, ALEA_CACHE_UNIVERSE)) return sys->has_fill;
    for (size_t i = 0; i < alea_vec_count(&sys->cells); i++) {
        if (sys->cells.data[i].fill_universe > 0)
            return true;
    }
    return false;
}

static int dedup_sorted_hits(alea_raycast_result_t* result) {
    if (!result) return -1;
    sort_hits(result->hits.data, result->hits.count);

    if (result->hits.count > 1) {
        size_t write = 1;
        for (size_t read = 1; read < result->hits.count; read++) {
            int same_t = fabs(result->hits.data[read].t -
                              result->hits.data[write - 1].t) <= DEDUP_EPSILON;
            int same_surf = result->hits.data[read].surface_id ==
                            result->hits.data[write - 1].surface_id;
            if (!(same_t && same_surf)) {
                result->hits.data[write++] = result->hits.data[read];
            }
        }
        result->hits.count = write;
    }

    return 0;
}

static int raycast_global_pipeline(alea_system_t* sys,
                                   const alea_ray_t* ray,
                                   double t_min, double t_max,
                                   bool include_lattice_hits,
                                   bool use_hier_lookup,
                                   alea_raycast_result_t* result) {
    if (use_hier_lookup) {
        raycast_surfaces_impl(sys, ray, t_min, t_max, result);
    } else {
        int rc = alea_raycast_surfaces(sys, ray, t_min, t_max, result);
        if (rc != 0) return rc;
    }

    if (system_has_fill_cells(sys))
        raycast_add_fill_hits(sys, ray, t_min, t_max, result);

    if (include_lattice_hits)
        raycast_add_lattice_hits(sys, ray, t_min, t_max, result);

    dedup_sorted_hits(result);
    return raycast_to_segments_impl(sys, t_max, result, use_hier_lookup);
}

int alea_raycast_global_reuse_nocache(alea_system_t* sys,
                                      const alea_ray_t* ray,
                                      double t_max,
                                      alea_raycast_result_t* result) {
    if (!sys || !ray || !result) return -1;

    /* This is the reusable counterpart to alea_raycast(). It intentionally
     * assumes the caller prepared the raycast caches, as render loops do. */
    alea_raycast_result_clear(result);
    double effective_t_max = (t_max <= 0) ? DBL_MAX : t_max;
    return raycast_global_pipeline(sys, ray, 0, effective_t_max,
                                   system_has_lattice_cells(sys), false,
                                   result);
}

int alea_raycast_global_breakpoints_reuse_nocache(
    alea_system_t* sys, const alea_ray_t* ray, double t_max,
    alea_raycast_result_t* result) {
    if (!sys || !ray || !result) return -1;

    /* Keep this deliberately below raycast_global_pipeline(): that helper
     * materializes precedence-selected segments as a compatibility service.
     * Coverage diagnostics need only the complete geometric partition. */
    alea_raycast_result_clear(result);
    const double effective_t_max = t_max <= 0 ? DBL_MAX : t_max;
    if (raycast_surfaces_impl(sys, ray, 0, effective_t_max, result) != 0)
        return -1;
    if (system_has_fill_cells(sys))
        raycast_add_fill_hits(sys, ray, 0, effective_t_max, result);
    if (system_has_lattice_cells(sys))
        raycast_add_lattice_hits(sys, ray, 0, effective_t_max, result);
    dedup_sorted_hits(result);
    return 0;
}

int alea_raycast_validation_breakpoints_reuse_nocache(
    alea_system_t* sys, const alea_ray_t* ray, double t_min, double t_max,
    size_t hit_limit, alea_raycast_result_t* result) {
    if (!sys || !ray || !result || !(t_max > t_min) || t_min < 0.0 ||
        !isfinite(t_max) || hit_limit == 0) return -1;
    alea_clear_error_detail();
    result->breakpoint_hit_limit = hit_limit;
    result->breakpoint_t_max = t_max;
    result->breakpoint_failed = 0;
    int rc = raycast_surfaces_impl(sys, ray, t_min, t_max, result);
    if (rc == 0 && !result->breakpoint_failed && system_has_fill_cells(sys)) {
#ifdef ALEA_RAY_PROFILE
        const double world_t = result->hits.count ? result->hits.data[0].t : INFINITY;
#endif
        raycast_add_fill_hits(sys, ray, t_min, t_max, result);
#ifdef ALEA_RAY_PROFILE
        if (result->breakpoint_nearest_only && !result->breakpoint_failed) {
            RAY_PROFILE_COUNT(fill_passes);
            if (result->hits.count && result->hits.data[0].t < world_t)
                RAY_PROFILE_COUNT(fill_improved);
        }
#endif
    }
    if (rc == 0 && !result->breakpoint_failed && system_has_lattice_cells(sys))
        raycast_add_lattice_hits(sys, ray, t_min, t_max, result);
    if (result->breakpoint_failed) rc = -1;
    if (rc == 0) sort_hits(result->hits.data, result->hits.count);
    result->breakpoint_hit_limit = 0;
    if (rc != 0 && alea_get_last_error() == ALEA_OK)
        alea_set_error_detail(ALEA_ERR_INVALID_STATE,
                              "particle interval breakpoint enumeration incomplete");
    return rc;
}

void alea_ray_boundary_event_result_init(alea_ray_boundary_event_result_t* result) {
    if (!result) return;
    memset(result, 0, sizeof(*result));
}

void alea_ray_boundary_event_result_clear(alea_ray_boundary_event_result_t* result) {
    if (!result) return;
    alea_vec_clear(&result->events);
    result->breakpoint_hits = 0;
    result->selected_segments = 0;
}

void alea_ray_boundary_event_result_free(alea_ray_boundary_event_result_t* result) {
    if (!result) return;
    alea_vec_free(&result->events);
}

static const alea_ray_hit_t* boundary_event_find_hit(
    const alea_raycast_result_t* crossings, double t, int surface_id) {
    for (size_t i = 0; i < crossings->hits.count; i++) {
        const alea_ray_hit_t* hit = &crossings->hits.data[i];
        if (hit->surface_id == surface_id &&
            fabs(hit->t - t) <= RAY_EPSILON)
            return hit;
    }
    return NULL;
}

static int boundary_event_append(alea_ray_boundary_event_result_t* events,
                                 const alea_ray_segment_t* before,
                                 const alea_ray_segment_t* after,
                                 double t,
                                 alea_ray_boundary_event_kind_t kind,
                                 int surface_id,
                                 const alea_ray_hit_t* hit,
                                 const alea_ray_boundary_event_options_internal_t* options) {
    if (options &&
        ((options->max_events && events->events.count >= options->max_events) ||
         (options->max_output_bytes &&
          events->events.count >=
              options->max_output_bytes / sizeof(*events->events.data)))) {
        alea_set_error_detail(ALEA_ERR_OVERFLOW,
                              "boundary-event query output limit exceeded");
        return -1;
    }
    alea_ray_boundary_event_t event = {
        .t = t,
        .kind = kind,
        .surface_id = surface_id,
        .primitive_id = hit ? hit->primitive_id : UINT32_MAX,
        .cell_before = before->cell_id,
        .cell_after = after->cell_id,
        .material_before = before->material_id,
        .material_after = after->material_id,
        .resolution_flags = before->resolution_flags | after->resolution_flags,
        .nx = hit ? hit->nx : 0,
        .ny = hit ? hit->ny : 0,
        .nz = hit ? hit->nz : 0
    };
    return alea_vec_push(&events->events, event, alea_ray_boundary_event_t);
}

static int boundary_events_append_group(
    const alea_raycast_result_t* crossings,
    const alea_ray_segment_t* before, const alea_ray_segment_t* after,
    bool include_all_coincident_physical,
    const alea_ray_boundary_event_options_internal_t* options,
    alea_ray_boundary_event_result_t* events) {
    const double t = before->t_exit;
    int ownership_changed = before->cell_id != after->cell_id ||
                            before->material_id != after->material_id;
    int lowest_physical_id = -1;
    int* physical_ids = NULL;
    size_t physical_count = 0;
    size_t physical_capacity = 0;
    int saw_synthetic = 0;

    for (size_t i = 0; i < crossings->hits.count; i++) {
        const alea_ray_hit_t* hit = &crossings->hits.data[i];
        if (fabs(hit->t - t) > RAY_EPSILON) continue;
        if (hit->surface_id == 0) {
            saw_synthetic = 1;
            continue;
        }
        if (hit->surface_id <= 0) continue;
        if (!include_all_coincident_physical) {
            if (lowest_physical_id < 0 || hit->surface_id < lowest_physical_id)
                lowest_physical_id = hit->surface_id;
            continue;
        }
        size_t j;
        for (j = 0; j < physical_count; j++)
            if (physical_ids[j] == hit->surface_id) break;
        if (j == physical_count) {
            if (physical_count == physical_capacity) {
                size_t next_capacity = physical_capacity ? physical_capacity * 2 : 16;
                if (next_capacity < physical_capacity ||
                    next_capacity > (size_t)-1 / sizeof(*physical_ids)) {
                    alea_set_error_detail(ALEA_ERR_OUT_OF_MEMORY,
                                          "too many coincident boundary surfaces");
                    free(physical_ids);
                    return -1;
                }
                int* next_ids = realloc(physical_ids,
                                        next_capacity * sizeof(*physical_ids));
                if (!next_ids) {
                    alea_set_error_detail(ALEA_ERR_OUT_OF_MEMORY,
                                          "failed to grow coincident boundary surface list");
                    free(physical_ids);
                    return -1;
                }
                physical_ids = next_ids;
                physical_capacity = next_capacity;
            }
            physical_ids[physical_count++] = hit->surface_id;
        }
    }

    /* A segment can carry a synthetic boundary generated by traversal even
     * when its t differs microscopically from the surface hit group. */
    if (before->exit_surface_id == 0 || after->enter_surface_id == 0)
        saw_synthetic = 1;

    if (ownership_changed &&
        (include_all_coincident_physical ? physical_count != 0 : lowest_physical_id >= 0)) {
        /* Surface ID is the deterministic reporting tie-breaker: the lowest
         * positive MC surface ID wins in the normal path.  Diagnostic mode
         * preserves every physical participant at this crossing. */
        if (include_all_coincident_physical) {
            for (size_t i = 0; i < physical_count; i++) {
                for (size_t j = i + 1; j < physical_count; j++) {
                    if (physical_ids[j] < physical_ids[i]) {
                        int tmp = physical_ids[i];
                        physical_ids[i] = physical_ids[j];
                        physical_ids[j] = tmp;
                    }
                }
            }
        }
        size_t emit_count = include_all_coincident_physical ? physical_count : 1;
        for (size_t i = 0; i < emit_count; i++) {
            int surface_id = include_all_coincident_physical ? physical_ids[i] :
                                                            lowest_physical_id;
            const alea_ray_hit_t* hit = boundary_event_find_hit(crossings, t,
                                                                 surface_id);
            if (boundary_event_append(events, before, after, t,
                                      ALEA_RAY_BOUNDARY_EVENT_PHYSICAL,
                                      surface_id, hit, options) != 0) {
                free(physical_ids);
                return -1;
            }
        }
    } else if (ownership_changed && !saw_synthetic) {
        if (boundary_event_append(events, before, after, t,
                                  ALEA_RAY_BOUNDARY_EVENT_UNRESOLVED,
                                  -1, NULL, options) != 0) {
            free(physical_ids);
            return -1;
        }
    }

    /* Physical events are emitted first; the synthetic event follows at the
     * same distance. This matches hit ordering and remains deterministic. */
    if (saw_synthetic &&
        boundary_event_append(events, before, after, t,
                              ALEA_RAY_BOUNDARY_EVENT_SYNTHETIC_LATTICE,
                              0, NULL, options) != 0) {
        free(physical_ids);
        return -1;
    }
    free(physical_ids);
    return 0;
}

/* Ownership transitions come from adjacent selected segments; physical
 * provenance comes from the independently enumerated crossing trace. */
static int boundary_events_from_segments_and_crossings(
    const alea_raycast_result_t* selected_trace,
    const alea_raycast_result_t* crossings,
    const alea_ray_boundary_event_options_internal_t* options,
    alea_ray_boundary_event_result_t* events);

int alea_raycast_boundary_events_with_options(
    alea_system_t* sys, const alea_ray_t* ray, double t_max,
    const alea_ray_boundary_event_options_internal_t* options,
    alea_raycast_result_t* trace, alea_ray_boundary_event_result_t* events) {
    if (!trace || !events)
        return -1;
    alea_ray_boundary_event_result_clear(events);
    /* Boundary events need the same ownership resolver as public first-hit
     * queries.  The legacy global path eventually falls back to a recursive
     * point lookup, which is both slower and cannot resolve filled/lattice
     * models reliably once the hierarchical cache is the active backend. */
    alea_raycast_result_clear(trace);
    double effective_t_max = (t_max <= 0) ? DBL_MAX : t_max;
    alea_raycast_result_t crossings;
    alea_raycast_result_init(&crossings);
    const bool use_hier_blas = options && options->use_hier_blas;
    if ((use_hier_blas
            ? raycast_root_blas_surfaces(
                sys, ray, 0.0, effective_t_max, &crossings)
            : alea_raycast_global_breakpoints_reuse_nocache(
                sys, ray, effective_t_max, &crossings)) != 0) {
        alea_raycast_result_free(&crossings);
        alea_set_error_detail(ALEA_ERR_INVALID_STATE,
                              "boundary-event breakpoint enumeration failed");
        return -1;
    }
    int trace_rc;
    if (use_hier_blas) {
        trace_rc = raycast_root_blas_surfaces(
            sys, ray, 0.0, effective_t_max, trace);
        if (trace_rc == 0)
            trace_rc = raycast_to_segments_impl(
                sys, effective_t_max, trace, true);
    } else {
        trace_rc = raycast_global_pipeline(
            sys, ray, 0, effective_t_max,
            system_has_lattice_cells(sys), true, trace);
    }
    if (trace_rc != 0) {
        alea_raycast_result_free(&crossings);
        alea_set_error_detail(ALEA_ERR_INVALID_STATE,
                              "boundary-event trace failed while resolving ray intervals");
        return -1;
    }
    const int rc = boundary_events_from_segments_and_crossings(
        trace, &crossings, options, events);
    events->breakpoint_hits = crossings.hits.count;
    events->selected_segments = trace->segments.count;
    alea_raycast_result_free(&crossings);
    return rc;
}

/* This bridge accepts separately-produced selected ownership and crossing
 * provenance traces. Fast compatibility backends may provide one trace for
 * both; global diagnostic backends deliberately do not. */
static int boundary_events_from_segments_and_crossings(
    const alea_raycast_result_t* selected_trace,
    const alea_raycast_result_t* crossings,
    const alea_ray_boundary_event_options_internal_t* options,
    alea_ray_boundary_event_result_t* events) {
    if (!selected_trace || !crossings || !events)
        return -1;
    alea_ray_boundary_event_result_clear(events);

    for (size_t i = 0; i + 1 < selected_trace->segments.count; i++) {
        const alea_ray_segment_t* before = &selected_trace->segments.data[i];
        const alea_ray_segment_t* after = &selected_trace->segments.data[i + 1];
        if (fabs(before->t_exit - after->t_enter) > RAY_EPSILON)
            continue;  /* Incomplete trace: do not invent an event. */

        if (boundary_events_append_group(
                crossings, before, after,
                options && options->include_all_coincident_physical,
                options,
                events) != 0) {
            alea_ray_boundary_event_result_clear(events);
            return -1;
        }
    }
    return 0;
}

int alea_raycast_boundary_events_reuse_nocache(
    alea_system_t* sys, const alea_ray_t* ray, double t_max,
    alea_raycast_result_t* trace, alea_ray_boundary_event_result_t* events) {
    return alea_raycast_boundary_events_with_options(sys, ray, t_max, NULL,
                                                      trace, events);
}

static void ray_query_output_clear(alea_ray_query_output_t* output) {
    if (!output) return;
    memset(output, 0, sizeof(*output));
    output->first_cell_id = -1;
    output->first_visible.cell_id = -1;
    output->first_visible.material_id = 0;
    output->first_visible.surface_id = -1;
    output->first_visible.primitive_id = UINT32_MAX;
}

static bool ray_plan_accepts_material(const alea_ray_plan_t* plan,
                                       const alea_ray_segment_t* seg,
                                       bool require_nonvoid) {
    if (seg->cell_id < 0) return false;
    if (require_nonvoid && seg->material_id == 0) return false;
    return plan->material_filter < 0 ||
           seg->material_id == plan->material_filter;
}

static const alea_ray_hit_t* ray_query_enter_hit(
    const alea_raycast_result_t* trace, const alea_ray_segment_t* seg) {
    if (seg->enter_hit_index >= 0 &&
        (size_t)seg->enter_hit_index < trace->hits.count)
        return &trace->hits.data[seg->enter_hit_index];
    return boundary_event_find_hit(trace, seg->t_enter, seg->enter_surface_id);
}

/* Keep scalar SEGMENTS queries consistent with the compact executor: an
 * interval that begins before the requested range becomes a cross-section,
 * not a surface entry.  This mutates only the reusable logical contents and
 * therefore retains all caller-owned capacity. */
static void ray_query_clip_segments(alea_raycast_result_t* trace,
                                    double t_min) {
    if (t_min <= 0.0) return;
    size_t write = 0;
    for (size_t i = 0; i < trace->segments.count; i++) {
        alea_ray_segment_t segment = trace->segments.data[i];
        if (segment.t_exit <= t_min + RAY_EPSILON) continue;
        if (segment.t_enter < t_min) {
            segment.t_enter = t_min;
            segment.enter_surface_id = -1;
            segment.enter_hit_index = -1;
        }
        trace->segments.data[write++] = segment;
    }
    trace->segments.count = write;
}

static int ray_query_fast_forward_trace(alea_system_t* sys,
                                        const alea_ray_t* ray,
                                        double t_max, bool emit_hits,
                                        alea_raycast_result_t* trace) {
    alea_raycast_result_clear(trace);
    trace->ray = *ray;
    return ray_selected_interval_trace(sys, ray, t_max, emit_hits, trace);
}

static int ray_query_fast_reverse_trace(alea_system_t* sys,
                                        const alea_ray_t* ray,
                                        double t_max, bool emit_hits,
                                        alea_raycast_result_t* trace) {
    if (!isfinite(t_max) || t_max <= 0 || t_max == DBL_MAX) {
        alea_set_error_detail(ALEA_ERR_INVALID_ARG,
                              "FAST_REVERSE requires a finite positive t_max");
        return -1;
    }
    double ex, ey, ez;
    alea_ray_point_at(ray, t_max, &ex, &ey, &ez);
    alea_ray_t reverse_ray;
    alea_ray_init_normalized(&reverse_ray, ex, ey, ez,
                             -ray->dx, -ray->dy, -ray->dz);

    alea_raycast_result_t reverse;
    alea_raycast_result_init(&reverse);
    int rc = ray_query_fast_forward_trace(sys, &reverse_ray, t_max,
                                           emit_hits, &reverse);
    if (rc != 0) {
        alea_raycast_result_free(&reverse);
        return -1;
    }
    alea_raycast_result_clear(trace);
    trace->ray = *ray;
    if (alea_vec_reserve(&trace->segments, reverse.segments.count,
                         alea_ray_segment_t) != 0) {
        alea_raycast_result_free(&reverse);
        return -1;
    }
    for (size_t i = 0; i < reverse.segments.count; i++) {
        const alea_ray_segment_t* src =
            &reverse.segments.data[reverse.segments.count - 1 - i];
        alea_ray_segment_t dst = *src;
        dst.t_enter = t_max - src->t_exit;
        dst.t_exit = t_max - src->t_enter;
        dst.enter_surface_id = src->exit_surface_id;
        dst.exit_surface_id = src->enter_surface_id;
        dst.enter_hit_index = -1;
        trace->segments.data[trace->segments.count++] = dst;
    }
    /* Reverse-hit normalization is intentionally deferred.  Fast boundary
     * provenance requests retain the verified forward trace in bidirectional
     * mode, and FAST_REVERSE exposes segment ownership only for now. */
    alea_raycast_result_free(&reverse);
    return 0;
}

static int ray_query_segments_match(const alea_raycast_result_t* forward,
                                    const alea_raycast_result_t* reverse) {
    if (forward->segments.count != reverse->segments.count)
        return 0;
    for (size_t i = 0; i < forward->segments.count; i++) {
        const alea_ray_segment_t* a = &forward->segments.data[i];
        const alea_ray_segment_t* b = &reverse->segments.data[i];
        if (fabs(a->t_enter - b->t_enter) > RAY_EPSILON ||
            fabs(a->t_exit - b->t_exit) > RAY_EPSILON ||
            a->cell_id != b->cell_id || a->material_id != b->material_id)
            return 0;
    }
    return 1;
}

int alea_raycast_query_reuse_nocache(
    alea_system_t* sys, const alea_ray_t* ray,
    const alea_ray_query_t* query, alea_raycast_result_t* trace,
    alea_ray_boundary_event_result_t* events,
    alea_ray_query_output_t* output) {
    /* Clear every valid caller-owned output before validation too: malformed
     * descriptors must not leave a previous ray's answer observable. */
    ray_query_output_clear(output);
    if (events) alea_ray_boundary_event_result_clear(events);
    if (trace) alea_raycast_result_clear(trace);
    if (!sys || !ray || !query || !trace) return -1;
    alea_ray_plan_t plan;
    if (alea_ray_query_lower(query, &plan) != 0)
        return -1;
    if (plan.product == ALEA_RAY_QUERY_BOUNDARY_EVENTS && !events)
        return -1;

    const double t_max = plan.t_max;
    const bool selected_scalar_backend =
        plan.backend == ALEA_RAY_QUERY_BACKEND_AUTO ||
        plan.backend == ALEA_RAY_QUERY_BACKEND_FAST_FORWARD;
    if (selected_scalar_backend &&
        plan.product == ALEA_RAY_QUERY_FIRST_CELL && output) {
        if (alea_raycast_hier_first_cell_nocache(
                sys, ray, plan.t_min, t_max, plan.material_filter,
                trace, &output->first_cell_id,
                &output->first_cell_t) != 0)
            goto fail;
        return 0;
    }
    if (selected_scalar_backend &&
        plan.product == ALEA_RAY_QUERY_FIRST_VISIBLE) {
        /* The hierarchical stepper validates each interval before returning,
         * then stops before allocating its segment/hit representation. */
        if (!output || alea_raycast_hier_first_visible_nocache(
                            sys, ray, plan.t_min, t_max,
                            plan.material_filter,
                            plan.requirements.need_normal != 0,
                            trace,
                            &output->first_visible) != 0)
            goto fail;
        return 0;
    }
    if (selected_scalar_backend &&
        plan.product == ALEA_RAY_QUERY_ANY_HIT) {
        int any_hit = 0;
        if (alea_raycast_hier_any_hit_nocache(
                sys, ray, plan.t_min, t_max, plan.material_filter,
                trace, &any_hit) != 0)
            goto fail;
        if (output) output->any_hit = any_hit != 0;
        return 0;
    }
    if (plan.product == ALEA_RAY_QUERY_BOUNDARY_EVENTS &&
        plan.engine == ALEA_RAY_ENGINE_GLOBAL_BREAKPOINTS) {
        const alea_ray_boundary_event_options_internal_t event_options = {
            .include_all_coincident_physical = true
        };
        if (alea_raycast_boundary_events_with_options(
                sys, ray, t_max, &event_options, trace, events) != 0)
            goto fail;
        size_t write = 0;
        for (size_t i = 0; i < events->events.count; i++) {
            const alea_ray_boundary_event_t event = events->events.data[i];
            if (event.t + RAY_EPSILON < plan.t_min ||
                event.t > t_max + RAY_EPSILON)
                continue;
            events->events.data[write++] = event;
        }
        events->events.count = write;
        if ((plan.max_events && write > plan.max_events) ||
            (plan.max_output_bytes &&
             write > plan.max_output_bytes / sizeof(*events->events.data)))
            goto fail;
        return 0;
    }
    if (plan.backend == ALEA_RAY_QUERY_BACKEND_AUTO &&
        plan.product == ALEA_RAY_QUERY_SEGMENTS) {
        /* Selected intervals always use the coherent hierarchical walker.
         * Lattices used to fall through to the global tracer here, which
         * made AUTO depend on model topology rather than query semantics. */
        if (ray_query_fast_forward_trace(sys, ray, t_max, false, trace) != 0)
            goto fail;
    } else if (plan.backend == ALEA_RAY_QUERY_BACKEND_FAST_FORWARD ||
        plan.backend == ALEA_RAY_QUERY_BACKEND_FAST_REVERSE ||
        plan.backend == ALEA_RAY_QUERY_BACKEND_FAST_FORWARD_REVERSE) {
        const bool emit_hits = plan.product == ALEA_RAY_QUERY_BOUNDARY_EVENTS ||
            (plan.product == ALEA_RAY_QUERY_FIRST_VISIBLE &&
             plan.requirements.need_surface_identity);
        if (plan.backend == ALEA_RAY_QUERY_BACKEND_FAST_FORWARD) {
            if (ray_query_fast_forward_trace(sys, ray, t_max, emit_hits, trace) != 0)
                goto fail;
        } else if (plan.backend == ALEA_RAY_QUERY_BACKEND_FAST_REVERSE) {
            if (ray_query_fast_reverse_trace(sys, ray, t_max, emit_hits, trace) != 0)
                goto fail;
        } else {
            alea_raycast_result_t reverse;
            alea_raycast_result_init(&reverse);
            int rc = ray_query_fast_forward_trace(sys, ray, t_max, emit_hits, trace);
            if (rc == 0)
                rc = ray_query_fast_reverse_trace(sys, ray, t_max, false, &reverse);
            if (rc == 0 && !ray_query_segments_match(trace, &reverse) && output)
                output->directional_mismatch = true;
            alea_raycast_result_free(&reverse);
            if (rc != 0) goto fail;
        }
    } else if (alea_raycast_global_reuse_nocache(sys, ray, t_max, trace) != 0) {
        goto fail;
    }

    if (plan.product == ALEA_RAY_QUERY_BOUNDARY_EVENTS) {
        /* Apply output limits only after the requested ray range has clipped
         * the event stream, matching the global boundary adapter above. */
        const alea_ray_boundary_event_options_internal_t event_options = {0};
        if (boundary_events_from_segments_and_crossings(
                trace, trace, &event_options, events) != 0)
            goto fail;
        size_t write = 0;
        for (size_t i = 0; i < events->events.count; i++) {
            const alea_ray_boundary_event_t event = events->events.data[i];
            if (event.t + RAY_EPSILON < plan.t_min ||
                event.t > t_max + RAY_EPSILON)
                continue;
            events->events.data[write++] = event;
        }
        events->events.count = write;
        if ((plan.max_events && write > plan.max_events) ||
            (plan.max_output_bytes &&
             write > plan.max_output_bytes / sizeof(*events->events.data)))
            goto fail;
        return 0;
    }

    if (plan.product == ALEA_RAY_QUERY_SEGMENTS) {
        ray_query_clip_segments(trace, plan.t_min);
        const size_t count = trace->segments.count;
        if ((plan.max_events && count > plan.max_events) ||
            (plan.max_output_bytes &&
             count > plan.max_output_bytes / sizeof(*trace->segments.data)))
            goto fail;
        return 0;
    }

    for (size_t i = 0; i < trace->segments.count; i++) {
        const alea_ray_segment_t* seg = &trace->segments.data[i];
        if (seg->t_exit <= plan.t_min + RAY_EPSILON)
            continue;

        if (plan.product == ALEA_RAY_QUERY_FIRST_CELL) {
            if (!ray_plan_accepts_material(&plan, seg, false)) continue;
            if (output) {
                output->first_cell_id = seg->cell_id;
                output->first_cell_t = fmax(seg->t_enter, plan.t_min);
            }
            return 0;
        }

        if (!ray_plan_accepts_material(&plan, seg, true)) continue;
        if (plan.product == ALEA_RAY_QUERY_ANY_HIT) {
            if (output) output->any_hit = true;
            return 0;
        }

        /* FIRST_VISIBLE: an interval clipped by t_min has no reportable
         * geometry boundary, matching the renderer's cross-section rule. */
        if (output) {
            alea_ray_first_visible_result_t* hit = &output->first_visible;
            const bool crossed_boundary =
                seg->t_enter >= plan.t_min - RAY_EPSILON;
            const alea_ray_hit_t* enter = crossed_boundary
                ? ray_query_enter_hit(trace, seg) : NULL;
            hit->found = true;
            hit->t = fmax(seg->t_enter, plan.t_min);
            hit->cell_id = seg->cell_id;
            hit->material_id = seg->material_id;
            hit->density = seg->density;
            hit->resolution_flags = seg->resolution_flags;
            if (enter) {
                hit->surface_id = enter->surface_id;
                hit->primitive_id = enter->primitive_id;
                if (plan.requirements.need_normal) {
                    hit->nx = enter->nx;
                    hit->ny = enter->ny;
                    hit->nz = enter->nz;
                }
            }
        }
        return 0;
    }
    return 0;

fail:
    alea_raycast_result_clear(trace);
    if (events) alea_ray_boundary_event_result_clear(events);
    ray_query_output_clear(output);
    return -1;
}

/* ============================================================================
 * CONVENIENCE FUNCTIONS
 * ============================================================================ */

int alea_raycast(alea_system_t* sys,
                double ox, double oy, double oz,
                double dx, double dy, double dz,
                double t_max,
                alea_raycast_result_t* result) {
    if (!result) return -1;

    /* alea_raycast() is the hit-producing surface+segment pipeline (uses the
     * surface BVH). For segment-only fast tracing use alea_raycast_hier(); for
     * hit-producing hierarchical tracing use alea_raycast_hier_with_hits(). */

    /* Free any prior allocations then reinitialize (safe for reuse) */
    alea_raycast_result_free(result);

    alea_ray_t ray;
    if (alea_ray_init(&ray, ox, oy, oz, dx, dy, dz) != 0) {
        return -1;  /* Zero-length direction */
    }

    /* t_max=0 means infinite */
    double effective_t_max = (t_max <= 0) ? DBL_MAX : t_max;
    return raycast_global_pipeline(sys, &ray, 0, effective_t_max,
                                   system_has_lattice_cells(sys), false,
                                   result);
}

int alea_raycast_hier(alea_system_t* sys,
                      double ox, double oy, double oz,
                      double dx, double dy, double dz,
                      double t_max,
                      alea_raycast_result_t* result) {
    if (!sys || !result) return -1;

    if (alea_system_prepare_query_caches(sys,
            ALEA_CACHE_HIER_SPATIAL | ALEA_CACHE_CELL_SURFACES) != 0)
        return -1;

    alea_raycast_result_free(result);

    alea_ray_t ray;
    if (alea_ray_init(&ray, ox, oy, oz, dx, dy, dz) != 0) {
        return -1;
    }

    double effective_t_max = (t_max <= 0) ? DBL_MAX : t_max;
    result->ray = ray;
    return ray_selected_interval_trace(sys, &ray, effective_t_max, false, result);
}

int alea_raycast_hier_blas_experimental(alea_system_t* sys,
                                        double ox, double oy, double oz,
                                        double dx, double dy, double dz,
                                        double t_max,
                                        alea_raycast_result_t* result) {
    if (!sys || !result) return -1;
    if (alea_raycast_ensure_hier_caches(sys) != 0) return -1;

    alea_raycast_result_free(result);

    alea_ray_t ray;
    if (alea_ray_init(&ray, ox, oy, oz, dx, dy, dz) != 0) {
        return -1;
    }
    result->ray = ray;

    double effective_t_max = (t_max <= 0) ? DBL_MAX : t_max;
    int rc = raycast_root_blas_surfaces(sys, &ray, 0.0, effective_t_max,
                                        result);
    if (rc == 1) {
        return alea_raycast_hier(sys, ox, oy, oz, dx, dy, dz, t_max, result);
    }
    if (rc != 0) return rc;

    return raycast_to_segments_impl(sys, effective_t_max, result, true);
}

int alea_ray_first_cell(alea_system_t* sys,
                        double ox, double oy, double oz,
                       double dx, double dy, double dz,
                       double t_max,
                       double* out_t) {
    if (!sys) return -1;
    alea_ray_t ray;
    if (alea_ray_init(&ray, ox, oy, oz, dx, dy, dz) != 0) return -1;
    if (alea_raycast_ensure_hier_caches(sys) != 0) return -1;
    static ALEA_THREAD_LOCAL alea_raycast_result_t hier_tls_result;
    static ALEA_THREAD_LOCAL int hier_tls_init = 0;
    if (!hier_tls_init) {
        alea_raycast_result_init(&hier_tls_result);
        hier_tls_init = 1;
    }
    int first_cell = -1;
    if (alea_raycast_hier_first_cell_nocache(
            sys, &ray, 0, t_max, -1, &hier_tls_result, &first_cell, out_t) != 0)
        return -1;
    return first_cell;
}

int alea_ray_is_occluded(alea_system_t* sys,
                        double ox, double oy, double oz,
                        double dx, double dy, double dz,
                        double t_max) {
    if (!sys) return 0;

    if (raycast_prefers_hier_mode(sys)) {
        if (alea_raycast_ensure_hier_caches(sys) != 0) return 0;
        static ALEA_THREAD_LOCAL alea_raycast_result_t hier_tls_result;
        static ALEA_THREAD_LOCAL int hier_tls_init = 0;
        if (!hier_tls_init) {
            alea_raycast_result_init(&hier_tls_result);
            hier_tls_init = 1;
        }
        alea_ray_t ray;
        if (alea_ray_init(&ray, ox, oy, oz, dx, dy, dz) != 0) return 0;
        int occluded = 0;
        if (alea_raycast_hier_any_hit_nocache(
                sys, &ray, 0, t_max, -1, &hier_tls_result, &occluded) != 0)
            return 0;
        return occluded;
    }

    /* Reuse a thread-local result to avoid malloc per shadow ray */
    static ALEA_THREAD_LOCAL alea_raycast_result_t tls_result;
    static ALEA_THREAD_LOCAL int tls_init = 0;
    if (!tls_init) {
        alea_raycast_result_init(&tls_result);
        alea_raycast_result_reserve(&tls_result, 64, 32);
        tls_init = 1;
    }
    alea_raycast_result_clear(&tls_result);

    alea_ray_t ray;
    alea_ray_init(&ray, ox, oy, oz, dx, dy, dz);

    double effective_t_max = (t_max <= 0) ? DBL_MAX : t_max;

    /* Get surface intersections */
    if (alea_raycast_surfaces(sys, &ray, 0, effective_t_max, &tls_result) != 0)
        return 0;

    /* Sort hits by distance */
    sort_hits(tls_result.hits.data, tls_result.hits.count);

    /* Walk intervals and early-out on first non-void cell */
    double t_prev = 0;
    int prev_cell_idx = -1;

    for (size_t i = 0; i <= tls_result.hits.count; i++) {
        double t_curr = (i < tls_result.hits.count) ? tls_result.hits.data[i].t : effective_t_max;

        if (t_curr > t_prev + RAY_EPSILON) {
            int cell_id, cell_idx, material_id;
            double density;
            int crossed_surface = (i > 0) ? tls_result.hits.data[i - 1].surface_id : -1;
            find_cell_after_crossing(sys, &ray, t_prev, t_curr,
                                     prev_cell_idx, crossed_surface,
                                     false,
                                     &cell_id, &cell_idx, &material_id, &density);

            prev_cell_idx = cell_idx;

            /* Early out: found a non-void cell with non-void material */
            if (cell_id >= 0 && material_id != 0) {
                return 1;
            }
        }

        t_prev = t_curr;
    }

    return 0;
}

double alea_raycast_path_length(const alea_raycast_result_t* result,
                               int material_id) {
    if (!result) return 0;

    double total = 0;
    for (size_t i = 0; i < result->segments.count; i++) {
        const alea_ray_segment_t* seg = &result->segments.data[i];

        /* Skip infinite segments */
        if (seg->t_exit >= DBL_MAX - 1) continue;

        int match = (material_id < 0) ||  /* -1 = all materials */
                    (seg->material_id == material_id);

        if (match) {
            total += seg->t_exit - seg->t_enter;
        }
    }

    return total;
}

/* ============================================================================
 * CELL-AWARE RAYCAST
 * ============================================================================ */

/**
 * One boundary crossing produced by a single stepper iteration.
 *
 * The segment-only path consumes only `t` / `surface_id`; the hit-producing
 * path (alea_raycast_hier_with_hits) also needs to materialize a world-space
 * surface normal. Rather than carry the primitive payload and local ray by
 * value through the hot loop, we keep the winning primitive id plus the frame
 * transform (primitive-local -> world). At emit time the payload is fetched
 * from `primitive_id` and the local hit point is recovered from the world hit
 * point via the transform, so the segment-only path pays nothing for this.
 */
typedef struct {
    double t;                /* distance of the boundary along the world ray */
    int surface_id;          /* MCNP id: -1 none, 0 synthetic lattice, >0 physical */
    uint32_t primitive_id;   /* winning primitive, ALEA_PRIMITIVE_ID_INVALID if none */
    alea_matrix_t transform; /* primitive-local -> world frame of the winning surface */
    bool has_physical_surface;
    bool is_synthetic_lattice_boundary;
    uint8_t local_surface_count;
    uint8_t local_surface_complete;
    int local_surface_ids[16];
    uint8_t provenance_flags;
    int active_cell_id;
    int active_universe_id;
    int active_depth;
    uint64_t active_occurrence_key;
    uint64_t active_parent_occurrence_key;
    double local_point[3];
    double local_direction[3];
} alea_raycast_boundary_event_t;

/* Persistent portion of one hierarchical ray walk.  Keeping this explicit
 * lets the scalar adapter run to completion while a future packet executor
 * advances the same walk one verified interval at a time. */
typedef struct {
    double t_current;
    int prev_cell_idx;
    int prev_surface_id;
    int pending_enter_hit_index;
    alea_raycast_boundary_event_t enter_event;
    alea_hier_ray_path_t current_path;
    double pending_lattice_entry_sample;
    int iterations_remaining;
    /* Particle transport must reject a containment retry that never verifies
     * its selected owner; general ray queries retain legacy compatibility. */
    bool fail_unverified_ownership;
    size_t recovery_proof_budget;
    bool require_interior_sample;
    /* Navigation-only policy for competing hierarchy depths. Legacy ray
     * queries leave this zero and retain their ordinary numerical ordering. */
    double boundary_distance_tolerance;
    const alea_ray_root_lookup_t* root_lookup;
} alea_ray_walk_t;

#define ALEA_NAV_DEFAULT_BOUNDARY_DISTANCE_TOLERANCE 1.0e-6

static bool ray_crossing_distances_near(double reference_t,
                                        double first_t,
                                        double second_t,
                                        double relative_tolerance) {
    if (!(relative_tolerance > 0.0) || !isfinite(reference_t) ||
        !isfinite(first_t) || !isfinite(second_t)) return false;
    const double first_distance = first_t - reference_t;
    const double second_distance = second_t - reference_t;
    if (!(first_distance > 0.0) || !(second_distance > 0.0)) return false;
    const double scale = fmax(first_distance, second_distance);
    return fabs(first_distance - second_distance) <=
           relative_tolerance * scale;
}

static bool ray_boundaries_same_world_surface(
    uint32_t first_primitive, const alea_matrix_t* first_transform,
    uint32_t second_primitive, const alea_matrix_t* second_transform) {
    if (first_primitive == ALEA_PRIMITIVE_ID_INVALID ||
        first_primitive != second_primitive ||
        !first_transform || !second_transform) return false;
    for (int i = 0; i < 12; ++i) {
        if (first_transform->m[i] != second_transform->m[i]) return false;
    }
    return true;
}

/* One verified open interval produced by the selected-owner walk.  This is
 * deliberately private while compatibility APIs still publish their legacy
 * segment/hit vectors. Keeping it separate from those vectors makes the
 * geometric step independent of path capture and hit materialization. */
typedef struct {
    double t_enter;
    double t_exit;
    int cell_index;
    int cell_id;
    int material_id;
    double density;
    int enter_surface_id;
    int exit_surface_id;
    alea_raycast_boundary_event_t enter_event;
    alea_raycast_boundary_event_t exit_event;
    const alea_hier_ray_path_t* path;
    uint8_t owner_provenance_complete;
    uint64_t owner_occurrence_key;
    uint64_t owner_parent_occurrence_key;
    uint8_t resolution_flags;
    bool containment_verified;
} alea_ray_selected_interval_t;

static void alea_ray_walk_init(alea_ray_walk_t* state) {
    memset(state, 0, sizeof(*state));
    state->prev_cell_idx = -2;
    state->prev_surface_id = -1;
    state->pending_enter_hit_index = -1;
    state->enter_event.surface_id = -1;
    state->enter_event.primitive_id = ALEA_PRIMITIVE_ID_INVALID;
    state->enter_event.active_cell_id = -1;
    state->enter_event.active_universe_id = -1;
    state->enter_event.active_depth = -1;
    alea_matrix_identity(&state->enter_event.transform);
    state->iterations_remaining = MAX_RAY_WALK_ITERATIONS;
}

static void ray_path_copy_live(alea_hier_ray_path_t* dst,
                               const alea_hier_ray_path_t* src) {
    int count = src->count;
    if (count < 0) count = 0;
    if (count > ALEA_HIER_RAY_PATH_MAX) count = ALEA_HIER_RAY_PATH_MAX;
    dst->count = count;
    if (count > 0) {
        memcpy(dst->entries, src->entries,
               (size_t)count * sizeof(*dst->entries));
    }
}

/* Map a primitive-local normal back to world space.
 *
 * Normals transform by the inverse-transpose of the local->world linear part.
 * `mat->inv` is the world->local linear map, so its transpose is exactly that
 * inverse-transpose. For an identity transform this is a no-op, which keeps
 * non-lattice/non-fill normals bit-identical to the flat path. */
static void boundary_event_world_normal(const alea_matrix_t* mat,
                                        double nlx, double nly, double nlz,
                                        double* nx, double* ny, double* nz) {
    double wx = mat->inv[0] * nlx + mat->inv[4] * nly + mat->inv[8] * nlz;
    double wy = mat->inv[1] * nlx + mat->inv[5] * nly + mat->inv[9] * nlz;
    double wz = mat->inv[2] * nlx + mat->inv[6] * nly + mat->inv[10] * nlz;
    double len = sqrt(wx * wx + wy * wy + wz * wz);
    if (len > RAY_EPSILON) {
        *nx = wx / len;
        *ny = wy / len;
        *nz = wz / len;
    } else {
        *nx = nlx;
        *ny = nly;
        *nz = nlz;
    }
}

/* Emit one hit for the winning boundary event of a stepper iteration.
 * Returns the index of the pushed hit in result->hits, or -1 if no physical
 * surface hit was produced. */
static int boundary_event_emit_hit(alea_system_t* sys,
                                    const alea_ray_t* ray,
                                    const alea_raycast_boundary_event_t* ev,
                                    alea_raycast_result_t* result) {
    if (!ev->has_physical_surface ||
        ev->primitive_id == ALEA_PRIMITIVE_ID_INVALID ||
        ev->primitive_id >= alea_vec_count(&sys->primitives)) {
        return -1;
    }

    const alea_primitive_entry_t* prim = &sys->primitives.data[ev->primitive_id];
    alea_primitive_data_t prim_data;
    if (!raycast_primitive_copy_payload(sys, ev->primitive_id,
                                        prim->type, &prim_data)) {
        return -1;
    }

    double wx, wy, wz;
    alea_ray_point_at(ray, ev->t, &wx, &wy, &wz);

    /* Recover the primitive-local hit point and evaluate the local normal. */
    double lx = wx, ly = wy, lz = wz;
    alea_matrix_transform_point_inverse(&ev->transform, &lx, &ly, &lz);
    double nlx, nly, nlz;
    primitive_normal_at(prim->type, &prim_data, lx, ly, lz, &nlx, &nly, &nlz);

    alea_ray_hit_t hit;
    hit.t = ev->t;
    hit.surface_id = ev->surface_id;
    hit.primitive_id = ev->primitive_id;
    boundary_event_world_normal(&ev->transform, nlx, nly, nlz,
                                &hit.nx, &hit.ny, &hit.nz);

    int hit_index = (int)alea_vec_count(&result->hits);
    if (add_hit(result, &hit) != 0) {
        ALEA_LOG_WARN("add_hit failed (out of memory) - hier hits may be incomplete");
        return -1;
    }
    return hit_index;
}

/* Populate the one normal needed by FIRST_VISIBLE without allocating or
 * appending a hit record. */
static void boundary_event_first_visible_normal(
    alea_system_t* sys, const alea_ray_t* ray,
    const alea_raycast_boundary_event_t* ev,
    bool include_normal, alea_ray_first_visible_result_t* visible) {
    if (!ev->has_physical_surface ||
        ev->primitive_id == ALEA_PRIMITIVE_ID_INVALID ||
        ev->primitive_id >= alea_vec_count(&sys->primitives))
        return;
    visible->surface_id = ev->surface_id;
    visible->primitive_id = ev->primitive_id;
    if (!include_normal) return;
    alea_raycast_visible_normal(sys, ray, ev->t, &ev->transform, visible);
}

void alea_raycast_visible_normal(alea_system_t* sys, const alea_ray_t* ray, double t,
                                 const alea_matrix_t* frame, alea_ray_first_visible_result_t* visible) {
    const alea_primitive_entry_t* prim = &sys->primitives.data[visible->primitive_id];
    alea_primitive_data_t prim_data;
    if (!raycast_primitive_copy_payload(sys, visible->primitive_id,
                                        prim->type, &prim_data))
        return;
    double wx, wy, wz;
    alea_ray_point_at(ray, t, &wx, &wy, &wz);
    alea_matrix_transform_point_inverse(frame, &wx, &wy, &wz);
    double nlx, nly, nlz;
    primitive_normal_at(prim->type, &prim_data, wx, wy, wz,
                        &nlx, &nly, &nlz);
    boundary_event_world_normal(frame, nlx, nly, nlz,
                                &visible->nx, &visible->ny, &visible->nz);
}

int alea_ray_append_selected_segment(alea_raycast_result_t* result, int* previous_cell,
                                     int cell, const alea_ray_segment_t* segment) {
    if (result->segments.count && *previous_cell == cell) {
        alea_ray_segment_t* previous = &result->segments.data[result->segments.count - 1];
        if (previous->path_index == segment->path_index) {
            previous->t_exit = segment->t_exit;
            previous->exit_surface_id = segment->exit_surface_id;
            return 0;
        }
    }
    if (add_segment(result, segment)) return -1;
    *previous_cell = cell;
    return 0;
}

/* Compatibility publication adapter for a verified selected interval. The
 * walker does not need public paths, hit indices, or segment storage to
 * establish ownership and endpoints. */
static int ray_selected_interval_publish(
    alea_system_t* sys, const alea_ray_t* ray,
    const alea_ray_selected_interval_t* interval, bool use_hier_lookup,
    bool emit_hits, alea_raycast_result_t* result, int* io_prev_cell_idx,
    int* io_pending_enter_hit_index) {
    int exit_hit_index = -1;
    if (emit_hits) {
        exit_hit_index = boundary_event_emit_hit(
            sys, ray, &interval->exit_event, result);
    }

    uint32_t path_index = UINT32_MAX;
    if (use_hier_lookup && capture_hier_path(result, interval->path,
                                              &path_index) != 0) {
        return -1;
    }

    const alea_ray_segment_t segment = {
        .t_enter = interval->t_enter,
        .t_exit = interval->t_exit,
        .cell_id = interval->cell_id,
        .material_id = interval->material_id,
        .density = interval->density,
        .enter_surface_id = interval->enter_surface_id,
        .exit_surface_id = interval->exit_surface_id,
        .enter_hit_index = *io_pending_enter_hit_index,
        .resolution_flags = interval->resolution_flags,
        .path_index = path_index
    };
    if (alea_ray_append_selected_segment(result, io_prev_cell_idx, interval->cell_index, &segment))
        return -1;
    if (emit_hits) *io_pending_enter_hit_index = exit_hit_index;
    return 0;
}

static bool ray_selected_interval_first_visible(
    alea_system_t* sys, const alea_ray_t* ray,
    const alea_ray_selected_interval_t* interval, double t_min,
    int material_filter, bool include_normal,
    alea_ray_first_visible_result_t* out_visible) {
    if (interval->cell_id < 0 || interval->material_id == 0 ||
        (material_filter >= 0 && interval->material_id != material_filter) ||
        interval->t_exit <= t_min + RAY_EPSILON) {
        return false;
    }
    out_visible->found = true;
    out_visible->t = fmax(interval->t_enter, t_min);
    out_visible->cell_id = interval->cell_id;
    out_visible->material_id = interval->material_id;
    out_visible->density = interval->density;
    out_visible->resolution_flags = interval->resolution_flags;
    if (interval->t_enter >= t_min - RAY_EPSILON) {
        out_visible->surface_id = interval->enter_event.surface_id;
        boundary_event_first_visible_normal(
            sys, ray, &interval->enter_event, include_normal, out_visible);
    }
    return true;
}

/**
 * Test ray against a specific cell's surfaces only.
 * Returns closest hit distance, or DBL_MAX if no hit.
 *
 * When `out_primitive_id` is non-NULL it receives the canonical primitive id
 * of the closest surface (ALEA_PRIMITIVE_ID_INVALID if none), so the caller
 * can build a boundary event without re-scanning the cell.
 */
static double raycast_tree_closest_surface(alea_system_t* sys,
                                           const alea_ray_t* ray,
                                           alea_node_id_t node_id,
                                           double t_min,
                                           double t_max,
                                           alea_raycast_result_t* result,
                                           int* out_surface_id,
                                           uint32_t* out_primitive_id) {
    if (node_id >= alea_vec_count(&sys->nodes)) return DBL_MAX;
    const alea_node_t* node = &sys->nodes.data[node_id];
    alea_operation_t op = ALEA_GET_OPERATION(node);
    if (op == ALEA_OP_PRIMITIVE) {
        if (node->primitive.primitive_id >= alea_vec_count(&sys->primitives))
            return DBL_MAX;
        const alea_primitive_entry_t* prim =
            &sys->primitives.data[node->primitive.primitive_id];
        alea_primitive_data_t prim_data;
        if (!raycast_primitive_copy_payload(sys, node->primitive.primitive_id,
                                            prim->type, &prim_data)) {
            return DBL_MAX;
        }
        if (result) {
            result->surfaces_tested++;
            if ((unsigned)prim->type < ALEA_RAYCAST_PRIM_TYPE_BINS)
                result->prim_type_tests[prim->type]++;
        }
        double closest = DBL_MAX;
        double intersections[4];
        int count = ray_intersect_primitive(ray, prim->type, &prim_data,
                                            intersections);
        for (int i = 0; i < count; i++) {
            if (intersections[i] > t_min && intersections[i] < t_max &&
                intersections[i] < closest) {
                closest = intersections[i];
            }
        }
        if (closest < DBL_MAX) {
            *out_surface_id = node->primitive.mc_surface_id;
            if (out_primitive_id) *out_primitive_id = node->primitive.primitive_id;
        }
        return closest;
    }

    double left = raycast_tree_closest_surface(
        sys, ray, node->operation.left, t_min, t_max, result,
        out_surface_id, out_primitive_id);
    if (op == ALEA_OP_COMPLEMENT) return left;
    int right_surface_id = -1;
    uint32_t right_primitive_id = ALEA_PRIMITIVE_ID_INVALID;
    double right = raycast_tree_closest_surface(
        sys, ray, node->operation.right, t_min, t_max, result,
        &right_surface_id, &right_primitive_id);
    if (right < left) {
        *out_surface_id = right_surface_id;
        if (out_primitive_id) *out_primitive_id = right_primitive_id;
        return right;
    }
    return left;
}

static double raycast_cell_surfaces(alea_system_t* sys,
                                    const alea_ray_t* ray,
                                    const alea_cell_entry_t* cell,
                                    double t_min, double t_max,
                                    alea_raycast_result_t* result,
                                    int* out_surface_id,
                                    uint32_t* out_primitive_id) {
    double closest_t = DBL_MAX;
    *out_surface_id = -1;
    if (out_primitive_id) *out_primitive_id = ALEA_PRIMITIVE_ID_INVALID;

    /* A few imported/container cells deliberately have no usable surface
     * index.  Their CSG is still exact occurrence support, so retain a
     * bounded tree fallback for ancestor-event queries and traversal. */
    if (!cell->surface_indices || cell->surface_index_count == 0) {
        return raycast_tree_closest_surface(
            sys, ray, cell->root_node_id, t_min, t_max, result,
            out_surface_id, out_primitive_id);
    }

    /* Test each surface belonging to this cell */
    for (size_t i = 0; i < cell->surface_index_count; i++) {
        uint32_t surf_idx = cell->surface_indices[i];
        if (surf_idx >= alea_vec_count(&sys->surfaces)) continue;

        const alea_surface_entry_t* surf = &sys->surfaces.data[surf_idx];
        const alea_primitive_entry_t* prim = &sys->primitives.data[surf->primitive_id];
        alea_primitive_data_t prim_data;
        if (!raycast_primitive_copy_payload(sys, surf->primitive_id,
                                            prim->type, &prim_data)) continue;
        if (result) {
            result->surfaces_tested++;
            if ((unsigned)prim->type < ALEA_RAYCAST_PRIM_TYPE_BINS)
                result->prim_type_tests[prim->type]++;
        }

        double t[4];
        int count = ray_intersect_primitive(ray, prim->type, &prim_data, t);

        for (int j = 0; j < count; j++) {
            if (t[j] > t_min && t[j] < t_max && t[j] < closest_t) {
                closest_t = t[j];
                *out_surface_id = surf->mc_surface_id;
                if (out_primitive_id) *out_primitive_id = surf->primitive_id;
            }
        }
    }

    return closest_t;
}

int alea_ray_cell_surface_at(alea_system_t* sys, const alea_ray_t* local,
                             uint32_t cell, double t, alea_raycast_result_t* result,
                             uint32_t* primitive) {
    int surface;
    const double hit = raycast_cell_surfaces(sys, local, &sys->cells.data[cell],
        nextafter(t, -INFINITY), nextafter(t, INFINITY), result, &surface, primitive);
    return hit == t ? surface : -1;
}

static void boundary_event_add_local_surface(alea_raycast_boundary_event_t* event,
                                             int surface_id) {
    if (!event || surface_id <= 0) return;
    for (size_t i = 0; i < event->local_surface_count; i++)
        if (event->local_surface_ids[i] == surface_id) return;
    if (event->local_surface_count >=
        sizeof(event->local_surface_ids) / sizeof(event->local_surface_ids[0])) {
        event->local_surface_complete = 0;
        return;
    }
    event->local_surface_ids[event->local_surface_count++] = surface_id;
}

static int boundary_event_cell_uses_exact_surface(
    const alea_system_t* sys, uint32_t cell_index, int surface_id) {
    if (!sys || surface_id <= 0 || !sys->mc_id_to_surface ||
        (size_t)surface_id >= sys->mc_id_to_surface_size ||
        !sys->surface_cell_offsets || !sys->surface_cell_refs)
        return 0;
    const uint32_t surface_index = sys->mc_id_to_surface[surface_id];
    if (surface_index == UINT32_MAX ||
        surface_index >= alea_vec_count(&sys->surfaces)) return 0;
    const size_t begin = sys->surface_cell_offsets[surface_index];
    const size_t end = sys->surface_cell_offsets[surface_index + 1];
    for (size_t i = begin; i < end; i++)
        if (sys->surface_cell_refs[i].cell_index == cell_index) return 1;
    return 0;
}

static void boundary_event_path_keys_at(
    const alea_hier_ray_path_t* path, int index,
    uint64_t* out_key, uint64_t* out_parent_key) {
    uint64_t key = 0, parent = 0;
    alea_occurrence_state_t occurrence_state;
    alea_occurrence_state_init(&occurrence_state);
    if (!path || index < 0 || index >= path->count) {
        if (out_key) *out_key = 0;
        if (out_parent_key) *out_parent_key = 0;
        return;
    }
    for (int i = 0; i <= index; i++) {
        parent = key;
        key = path_entry_occurrence_key(
            &occurrence_state, &path->entries[i]);
    }
    if (out_key) *out_key = key;
    if (out_parent_key) *out_parent_key = parent;
}

static void boundary_interval_capture_owner(
    const alea_hier_ray_path_t* path, alea_ray_selected_interval_t* interval) {
    if (!path || path->count <= 0 || !interval) return;
    boundary_event_path_keys_at(
        path, path->count - 1, &interval->owner_occurrence_key,
        &interval->owner_parent_occurrence_key);
    interval->owner_provenance_complete =
        interval->owner_occurrence_key != 0;
}

static bool boundary_event_collect_cell_surfaces_at(
    alea_system_t* sys, const alea_ray_t* local_ray,
    const alea_cell_entry_t* cell, double t, alea_raycast_boundary_event_t* event) {
    if (!cell || !cell->surface_indices || !cell->surface_index_count) return false;
    bool scanned = true;
    for (size_t i = 0; i < cell->surface_index_count; i++) {
        const uint32_t surface_index = cell->surface_indices[i];
        if (surface_index >= alea_vec_count(&sys->surfaces)) {
            scanned = false;
            continue;
        }
        const alea_surface_entry_t* surface = &sys->surfaces.data[surface_index];
        if (surface->primitive_id >= alea_vec_count(&sys->primitives)) {
            scanned = false;
            continue;
        }
        const alea_primitive_entry_t* primitive =
            &sys->primitives.data[surface->primitive_id];
        alea_primitive_data_t data;
        if (!raycast_primitive_copy_payload(sys, surface->primitive_id,
                                            primitive->type, &data)) {
            event->local_surface_complete = 0;
            continue;
        }
        double hits[4];
        const int count = ray_intersect_primitive(local_ray, primitive->type,
                                                  &data, hits);
        for (int hit = 0; hit < count; hit++)
            if (fabs(hits[hit] - t) <= RAY_EPSILON)
                boundary_event_add_local_surface(event, surface->mc_surface_id);
    }
    return scanned && event->local_surface_complete;
}

static void boundary_event_collect_local_group(
    alea_system_t* sys, const alea_ray_t* ray,
    const alea_hier_ray_path_t* path, alea_raycast_boundary_event_t* event) {
    if (!event || !event->has_physical_surface ||
        event->is_synthetic_lattice_boundary || !path) return;
    event->local_surface_count = 0;
    event->local_surface_complete = 1;
    for (int i = 0; i < path->count; i++) {
        const alea_hier_ray_path_entry_t* entry = &path->entries[i];
        if ((size_t)entry->cell_index >= alea_vec_count(&sys->cells)) {
            event->local_surface_complete = 0;
            continue;
        }
        alea_ray_t local_ray;
        alea_ray_transform_inverse(&entry->transform, ray, &local_ray);
        boundary_event_collect_cell_surfaces_at(
            sys, &local_ray, &sys->cells.data[entry->cell_index], event->t, event);
        if (event->surface_id > 0 &&
            boundary_event_cell_uses_exact_surface(
                sys, entry->cell_index, event->surface_id)) {
            /* Prefer the deepest matching entry if an exact card is present
             * at more than one hierarchy level. */
            event->active_cell_id = entry->cell_id;
            event->active_universe_id = entry->universe_id;
            event->active_depth = entry->depth;
            boundary_event_path_keys_at(
                path, i, &event->active_occurrence_key,
                &event->active_parent_occurrence_key);
            alea_ray_point_at(
                &local_ray, event->t, &event->local_point[0],
                &event->local_point[1], &event->local_point[2]);
            /* The event reports a unit direction in its active frame. */
            double length = sqrt(local_ray.dx * local_ray.dx + local_ray.dy * local_ray.dy +
                                 local_ray.dz * local_ray.dz);
            event->local_direction[0] = local_ray.dx / length;
            event->local_direction[1] = local_ray.dy / length;
            event->local_direction[2] = local_ray.dz / length;
            event->provenance_flags |= ALEA_BOUNDARY_PROVENANCE_ACTIVE_FRAME;
        }
    }
    if (event->local_surface_count == 0) event->local_surface_complete = 0;
}

/* The incident path is empty for a ray entering from exterior space. The
 * global walker searches world-space surfaces there, so use the same frame to
 * retain every physical surface tied at its chosen crossing distance. */
static void boundary_event_collect_global_group(
    alea_system_t* sys, const alea_ray_t* ray,
    alea_raycast_boundary_event_t* event) {
    event->local_surface_count = 0;
    event->local_surface_complete = 1;
    for (size_t i = 0; i < alea_vec_count(&sys->surfaces); i++) {
        const alea_surface_entry_t* surface = &sys->surfaces.data[i];
        if (surface->primitive_id >= alea_vec_count(&sys->primitives)) {
            event->local_surface_complete = 0;
            continue;
        }
        const alea_primitive_entry_t* primitive =
            &sys->primitives.data[surface->primitive_id];
        alea_primitive_data_t data;
        if (!raycast_primitive_copy_payload(sys, surface->primitive_id,
                                            primitive->type, &data)) {
            event->local_surface_complete = 0;
            continue;
        }
        double hits[4];
        const int count = ray_intersect_primitive(ray, primitive->type,
                                                  &data, hits);
        for (int hit = 0; hit < count; hit++)
            if (fabs(hits[hit] - event->t) <= RAY_EPSILON)
                boundary_event_add_local_surface(event, surface->mc_surface_id);
    }
    if (event->local_surface_count == 0) event->local_surface_complete = 0;
}

/* A selected ownership transition is only a useful shortcut when both open
 * sides have unambiguous coverage.  The selected walker deliberately picks
 * one deck-precedence owner; it therefore cannot, by itself, expose a second
 * coincident or overlapping claimant hidden behind that choice.  Use the
 * diagnostic recursive coverage query at one interior point on each side.
 *
 * The hit list contains ancestry as well as terminal occurrences.  Counting
 * only occurrence keys that are not a parent of another returned occurrence
 * distinguishes a normal fill chain from genuinely competing terminals.
 */
#define SELECTED_BOUNDARY_COVERAGE_HIT_CAP 64
static int boundary_event_side_has_unique_terminal_owner(
    const alea_system_t* sys, const alea_ray_t* ray, double t_sample) {
    alea_cell_hit_t hits[SELECTED_BOUNDARY_COVERAGE_HIT_CAP];
    uint64_t occurrence_keys[SELECTED_BOUNDARY_COVERAGE_HIT_CAP];
    uint64_t parent_keys[SELECTED_BOUNDARY_COVERAGE_HIT_CAP];
    double x, y, z;
    alea_ray_point_at(ray, t_sample, &x, &y, &z);
    const int count = alea_find_all_cells_at_point_coverage_chain_recursive(
        sys, x, y, z, hits, occurrence_keys, parent_keys,
        SELECTED_BOUNDARY_COVERAGE_HIT_CAP);
    if (count < 0 || count >= SELECTED_BOUNDARY_COVERAGE_HIT_CAP) return 0;
    int terminal_count = 0;
    for (int i = 0; i < count; i++) {
        int is_parent = 0;
        for (int j = 0; j < count; j++) {
            if (parent_keys[j] == occurrence_keys[i] &&
                occurrence_keys[i] != 0) {
                is_parent = 1;
                break;
            }
        }
        if (!is_parent && ++terminal_count > 1) return 0;
    }
    return 1;
}

static void boundary_event_require_unambiguous_open_sides(
    const alea_system_t* sys, const alea_ray_t* ray,
    const alea_ray_selected_interval_t* before,
    const alea_ray_selected_interval_t* after,
    alea_ray_boundary_event_t* event) {
    if (!event || !event->local_surface_complete) return;
    const double left_span = before->t_exit - before->t_enter;
    const double right_span = after->t_exit - after->t_enter;
    /* The walker can identify all locally tied primitive intersections, but
     * not yet prove that every one is causally active in the cell's CSG
     * expression.  A multi-surface group therefore remains canonical-only;
     * publishing it would turn incidental tangencies into labels. */
    if (event->local_surface_count != 1 ||
        left_span <= RAY_EPSILON || right_span <= RAY_EPSILON ||
        !boundary_event_side_has_unique_terminal_owner(
            sys, ray, before->t_enter + 0.5 * left_span) ||
        !boundary_event_side_has_unique_terminal_owner(
            sys, ray, after->t_enter + 0.5 * right_span))
        event->local_surface_complete = 0;
}

int alea_raycast_shared_terminal_surface_nocache(
    alea_system_t* sys, int cell_a_id, int cell_b_id,
    const alea_ray_t* ray, double t_max, int* out_surface_id,
    double* out_fraction) {
    if (!sys || !ray || !out_surface_id || !out_fraction || t_max <= 0.0)
        return 0;
    const int a_index = alea_find_cell_by_id(sys, cell_a_id);
    const int b_index = alea_find_cell_by_id(sys, cell_b_id);
    if (a_index < 0 || b_index < 0 || a_index == b_index) return 0;
    const alea_cell_entry_t* a = &sys->cells.data[a_index];
    const alea_cell_entry_t* b = &sys->cells.data[b_index];
    /* Occurrence-specific placement is deliberately left to canonical rays. */
    if (a->universe_id != b->universe_id ||
        alea_cell_entry_is_container(a) || alea_cell_entry_is_container(b) ||
        !a->surface_indices || !b->surface_indices)
        return 0;
    int found_id = -1;
    double found_t = 0.0;
    for (size_t i = 0; i < a->surface_index_count; i++) {
        const uint32_t index = a->surface_indices[i];
        if (index >= alea_vec_count(&sys->surfaces)) continue;
        const alea_surface_entry_t* surface = &sys->surfaces.data[index];
        bool shared = false;
        int shared_surface_id = surface->mc_surface_id;
        for (size_t j = 0; j < b->surface_index_count; j++) {
            const uint32_t other_index = b->surface_indices[j];
            if (other_index >= alea_vec_count(&sys->surfaces)) continue;
            const alea_surface_entry_t* other = &sys->surfaces.data[other_index];
            /* MCNP conversion deduplicates equal geometry to one primitive
             * while retaining the original surface-card entries.  Treat the
             * primitive as the shared boundary and publish its lowest local
             * card ID; the sparse caller canonicalizes across all aliases. */
            if (other->primitive_id == surface->primitive_id) {
                shared = true;
                if (other->mc_surface_id < shared_surface_id)
                    shared_surface_id = other->mc_surface_id;
                break;
            }
        }
        if (!shared) continue;
        /* MCNP surface transforms are baked into the registered primitive at
         * conversion time.  transform_applied is retained only for round-trip
         * export and does not require another ray-coordinate transform. */
        if (surface->primitive_id >= alea_vec_count(&sys->primitives))
            return 0;
        const alea_primitive_entry_t* primitive =
            &sys->primitives.data[surface->primitive_id];
        alea_primitive_data_t data;
        if (!raycast_primitive_copy_payload(sys, surface->primitive_id,
                                            primitive->type, &data))
            return 0;
        double ts[4];
        const int count = ray_intersect_primitive(ray, primitive->type, &data, ts);
        for (int hit = 0; hit < count; hit++) {
            if (ts[hit] <= RAY_EPSILON || ts[hit] >= t_max - RAY_EPSILON)
                continue;
            if (found_id >= 0 && (found_id != shared_surface_id ||
                                  fabs(found_t - ts[hit]) > RAY_EPSILON))
                return 0;
            found_id = shared_surface_id;
            found_t = ts[hit];
        }
    }
    if (found_id <= 0) return 0;
    *out_surface_id = found_id;
    *out_fraction = found_t / t_max;
    return 1;
}

int alea_raycast_terminal_surface_nocache(
    alea_system_t* sys, int cell_id, const alea_ray_t* ray, double t_max,
    int* out_surface_id, double* out_fraction) {
    if (!sys || !ray || !out_surface_id || !out_fraction || t_max <= 0.0)
        return 0;
    const int index = alea_find_cell_by_id(sys, cell_id);
    if (index < 0) return 0;
    const alea_cell_entry_t* cell = &sys->cells.data[index];
    /* Container-region surfaces are exactly the causal boundary for many
     * transitions between hierarchy paths, so they are valid local evidence. */
    if (!cell->surface_indices) return 0;
    int found = -1; double found_t = 0.0;
    for (size_t i = 0; i < cell->surface_index_count; i++) {
        uint32_t si = cell->surface_indices[i];
        if (si >= alea_vec_count(&sys->surfaces)) return 0;
        const alea_surface_entry_t* s = &sys->surfaces.data[si];
        if (s->primitive_id >= alea_vec_count(&sys->primitives)) return 0;
        const alea_primitive_entry_t* p = &sys->primitives.data[s->primitive_id];
        alea_primitive_data_t data;
        if (!raycast_primitive_copy_payload(sys, s->primitive_id, p->type, &data)) return 0;
        double ts[4]; int n = ray_intersect_primitive(ray, p->type, &data, ts);
        for (int h = 0; h < n; h++) if (ts[h] > RAY_EPSILON && ts[h] < t_max - RAY_EPSILON) {
            if (found >= 0 && (found != s->mc_surface_id || fabs(found_t - ts[h]) > RAY_EPSILON)) return 0;
            found = s->mc_surface_id; found_t = ts[h];
        }
    }
    if (found < 0) return 0;
    *out_surface_id = found; *out_fraction = found_t / t_max; return 1;
}

static alea_raycast_ancestor_cell_stat_t*
raycast_ancestor_cell_stat(alea_raycast_result_t* result,
                           uint32_t cell_index, int cell_id) {
    for (size_t i = 0; i < ALEA_RAYCAST_ANCESTOR_HOT_CELLS; ++i) {
        alea_raycast_ancestor_cell_stat_t* stat = &result->ancestor_hot_cells[i];
        if (stat->queries != 0 && stat->cell_index == cell_index)
            return stat;
    }
    for (size_t i = 0; i < ALEA_RAYCAST_ANCESTOR_HOT_CELLS; ++i) {
        alea_raycast_ancestor_cell_stat_t* stat = &result->ancestor_hot_cells[i];
        if (stat->queries == 0) {
            stat->cell_index = cell_index;
            stat->cell_id = cell_id;
            return stat;
        }
    }
    return NULL;
}

static double raycast_hier_path_ancestor_surfaces(alea_system_t* sys,
                                                  const alea_ray_t* ray,
                                                  const alea_hier_ray_path_t* path,
                                                  uint32_t terminal_cell_index,
                                                  int already_tested_lattice_cell,
                                                  double t_min,
                                                  double t_max,
                                                  double reference_t,
                                                  double relative_tolerance,
                                                  alea_raycast_result_t* result,
                                                  int* out_surface_id,
                                                  uint32_t* out_primitive_id,
                                                  alea_matrix_t* out_transform) {
    double closest_t = DBL_MAX;
    *out_surface_id = -1;
    if (out_primitive_id) *out_primitive_id = ALEA_PRIMITIVE_ID_INVALID;

    if (!path || path->count <= 1) return closest_t;
    if (result) {
        uint32_t depth = (uint32_t)(path->count - 1);
        result->ancestor_surface_queries++;
        result->ancestor_path_entries_examined += depth;
        if (depth > result->ancestor_max_path_depth)
            result->ancestor_max_path_depth = depth;
    }

    for (int i = 0; i < path->count - 1; i++) {
        const alea_hier_ray_path_entry_t* entry = &path->entries[i];
        if (entry->cell_index == terminal_cell_index) continue;
        if ((int)entry->cell_index == already_tested_lattice_cell) continue;
        if ((size_t)entry->cell_index >= alea_vec_count(&sys->cells)) continue;

        const alea_cell_entry_t* cell = &sys->cells.data[entry->cell_index];
        if (!cell->surface_indices || cell->surface_index_count == 0) continue;

        alea_ray_t local_ray;
        alea_ray_transform_inverse(&entry->transform, ray, &local_ray);

        int surface_id = -1;
        uint32_t prim_id = ALEA_PRIMITIVE_ID_INVALID;
        int tested_before = result ? result->surfaces_tested : 0;
        double t = raycast_cell_surfaces(sys, &local_ray, cell, t_min, t_max,
                                         result, &surface_id, &prim_id);
        int tested = result ? result->surfaces_tested - tested_before : 0;
        alea_raycast_ancestor_cell_stat_t* stat = result
            ? raycast_ancestor_cell_stat(result, entry->cell_index, cell->mc_cell_id)
            : NULL;
        if (stat) {
            stat->queries++;
            stat->surface_tests += (uint64_t)tested;
        } else if (result) {
            result->ancestor_unattributed_queries++;
            result->ancestor_unattributed_surface_tests += (uint64_t)tested;
        }
        /* Entries are visited from root to leaf. Retain the shallower
         * crossing when distances are equal under the navigation policy. */
        const bool same_near_boundary = ray_crossing_distances_near(
                reference_t, t, closest_t, relative_tolerance) &&
            ray_boundaries_same_world_surface(
                prim_id, &entry->transform,
                out_primitive_id ? *out_primitive_id
                                 : ALEA_PRIMITIVE_ID_INVALID,
                out_transform);
        if (t < closest_t && !same_near_boundary) {
            closest_t = t;
            *out_surface_id = surface_id;
            if (out_primitive_id) *out_primitive_id = prim_id;
            if (out_transform) *out_transform = entry->transform;
            if (stat) stat->winning_events++;
            else if (result) result->ancestor_unattributed_winning_events++;
        }
    }

    return closest_t;
}

/* An unresolved ordinary fill owns the current point, but none of its
 * immediate children does. Its own surfaces cannot bound that undefined span:
 * a child may start before the container exits. Search only that universe in
 * its placed frame, stopping at the already known ancestor/container boundary.
 * Further descent happens after ownership is resolved at the next crossing. */
typedef struct {
    alea_system_t* sys;
    const alea_ray_t* ray;
    alea_raycast_result_t* result;
    double t_min;
    double closest_t;
    int surface_id;
    uint32_t primitive_id;
} undefined_fill_crossing_t;

static int visit_undefined_fill_child(uint32_t cell_index,
                                      double t_enter, double t_exit,
                                      void* context) {
    (void)t_exit;
    undefined_fill_crossing_t* query = context;
    if (t_enter > query->closest_t) return 0;
    if ((size_t)cell_index >= query->sys->cells.count) return -1;
    const alea_cell_entry_t* child = &query->sys->cells.data[cell_index];
    int surface_id = -1;
    uint32_t primitive_id = ALEA_PRIMITIVE_ID_INVALID;
    double t = raycast_cell_surfaces(query->sys, query->ray, child,
        query->t_min, query->closest_t, query->result,
        &surface_id, &primitive_id);
    if (t < query->closest_t) {
        query->closest_t = t;
        query->surface_id = surface_id;
        query->primitive_id = primitive_id;
    }
    return 0;
}

static int raycast_undefined_fill_crossing(
    alea_system_t* sys, const alea_ray_t* ray,
    const alea_hier_ray_path_t* path, double t_min, double t_max,
    alea_raycast_result_t* result, alea_raycast_boundary_event_t* event) {
    if (!path || path->count == 0) return 0;
    const alea_hier_ray_path_entry_t* parent = &path->entries[path->count - 1];
    if ((size_t)parent->cell_index >= sys->cells.count) return -1;
    const alea_cell_entry_t* cell = &sys->cells.data[parent->cell_index];
    /* Lattice element transitions retain their dedicated DDA handling. */
    if (parent->is_lattice || cell->fill_universe <= 0) return 0;
    const alea_matrix_t* fill =
        alea_cell_fill_matrix(sys, parent->cell_index);
    if (!fill) return -1;
    alea_matrix_t transform;
    alea_matrix_multiply(&transform, &parent->transform, fill);
    if (!alea_matrix_invert(&transform)) return -1;
    alea_ray_t local_ray;
    alea_ray_transform_inverse(&transform, ray, &local_ray);
    undefined_fill_crossing_t query = {
        .sys = sys, .ray = &local_ray, .result = result,
        .t_min = t_min, .closest_t = t_max, .surface_id = -1,
        .primitive_id = ALEA_PRIMITIVE_ID_INVALID
    };
    if (alea_hier_spatial_visit_universe_ray(sys, cell->fill_universe,
            &local_ray,
            t_min, t_max, visit_undefined_fill_child, &query) < 0) return -1;
    if (!(query.closest_t < t_max)) return 0;
    event->t = query.closest_t;
    event->surface_id = query.surface_id;
    event->primitive_id = query.primitive_id;
    event->transform = transform;
    event->has_physical_surface = query.surface_id > 0;
    event->is_synthetic_lattice_boundary = false;
    return 1;
}

static double lattice_rect_next_boundary(const alea_ray_t* ray,
                                         const alea_cell_entry_t* lat_cell,
                                         double t_min,
                                         double t_max) {
    int ni = lat_cell->lat_fill_dims[1] - lat_cell->lat_fill_dims[0] + 1;
    int nj = lat_cell->lat_fill_dims[3] - lat_cell->lat_fill_dims[2] + 1;
    int nk = lat_cell->lat_fill_dims[5] - lat_cell->lat_fill_dims[4] + 1;
    double px = lat_cell->lat_pitch[0];
    double py = lat_cell->lat_pitch[1];
    double pz = lat_cell->lat_pitch[2];
    const double* ll = lat_cell->lat_lower_left;

    if (px <= 0.0 || py <= 0.0 || (nk > 1 && pz <= 0.0)) return t_max;

    alea_bbox_t lat_bbox;
    if (lattice_rect_bounds(lat_cell, &lat_bbox) != 0) return t_max;

    double t_enter = 0.0, t_exit = t_max;
    const int unbounded = lat_cell->lat_fill_repeating ||
                          lat_cell->lat_outer_universe > 0;
    if (!unbounded &&
        !ray_bbox_slab_enter_exit(ray, &lat_bbox, 0.0, t_max,
                                  &t_enter, &t_exit)) return t_max;

    double sample_t = t_min + RAY_EPSILON;
    double sx = ray->ox + sample_t * ray->dx;
    double sy = ray->oy + sample_t * ray->dy;
    double sz = ray->oz + sample_t * ray->dz;

    int i = (int)floor((sx - ll[0]) / px);
    int j = (lat_cell->lat_fill_repeating || nj > 1)
        ? (int)floor((sy - ll[1]) / py) : 0;
    const int axial_tiling = nk > 1 ||
        (lat_cell->lat_fill_repeating && pz > 0.0);
    int k = axial_tiling
        ? (int)floor((sz - ll[2]) / pz) : 0;

    if (!unbounded) {
        if (i < 0) i = 0;
        if (i >= ni) i = ni - 1;
        if (j < 0) j = 0;
        if (j >= nj) j = nj - 1;
        if (k < 0) k = 0;
        if (k >= nk) k = nk - 1;
    }

    (void)t_enter;
    (void)t_exit;
    double t_next = t_max;
    if (fabs(ray->dx) > RAY_EPSILON) {
        double boundary = ll[0] + ((ray->dx > 0.0) ? (i + 1) : i) * px;
        double t = (boundary - ray->ox) / ray->dx;
        if (t > t_min + RAY_EPSILON && t < t_next) t_next = t;
    }
    if ((lat_cell->lat_fill_repeating || nj > 1) &&
        fabs(ray->dy) > RAY_EPSILON) {
        double boundary = ll[1] + ((ray->dy > 0.0) ? (j + 1) : j) * py;
        double t = (boundary - ray->oy) / ray->dy;
        if (t > t_min + RAY_EPSILON && t < t_next) t_next = t;
    }
    if (axial_tiling &&
        fabs(ray->dz) > RAY_EPSILON) {
        double boundary = ll[2] + ((ray->dz > 0.0) ? (k + 1) : k) * pz;
        double t = (boundary - ray->oz) / ray->dz;
        if (t > t_min + RAY_EPSILON && t < t_next) t_next = t;
    }

    return (t_next < t_max) ? t_next : t_max;
}

static double lattice_hex_next_boundary(const alea_ray_t* ray,
                                        const alea_cell_entry_t* lat_cell,
                                        double t_min,
                                        double t_max) {
    double p = lat_cell->lat_pitch[0];
    if (p <= 0.0) return t_max;

    alea_bbox_t lat_bbox;
    if (lattice_hex_bounds(lat_cell, &lat_bbox) != 0) return t_max;

    double t_enter, t_exit;
    if (lat_cell->lat_outer_universe <= 0 &&
        !ray_bbox_slab_enter_exit(ray, &lat_bbox, 0.0, t_max,
                                  &t_enter, &t_exit)) {
        return t_max;
    }

    double inv_p = 1.0 / p;
    double inv_ps = inv_p / M_SQRT3;
    double dq = ray->dx * inv_p - ray->dy * inv_ps;
    double dr = 2.0 * ray->dy * inv_ps;
    double ds = -dq - dr;

    double sample_t = t_min + RAY_EPSILON;
    double sx = ray->ox + sample_t * ray->dx;
    double sy = ray->oy + sample_t * ray->dy;
    double q = sx * inv_p - sy * inv_ps;
    double rr = 2.0 * sy * inv_ps;
    double s = -q - rr;

    (void)t_enter;
    (void)t_exit;
    double t_next = t_max;
#define CONSIDER_NEXT_HALF(val, rate) do {                         \
        if (fabs(rate) > RAY_EPSILON) {                            \
            double bnd;                                            \
            if ((rate) > 0.0) {                                    \
                bnd = floor((val) + 0.5) + 0.5;                   \
                if (bnd <= (val) + RAY_EPSILON) bnd += 1.0;        \
            } else {                                               \
                bnd = ceil((val) - 0.5) - 0.5;                    \
                if (bnd >= (val) - RAY_EPSILON) bnd -= 1.0;        \
            }                                                      \
            double t = sample_t + (bnd - (val)) / (rate);          \
            if (t > t_min + RAY_EPSILON && t < t_next) t_next = t; \
        }                                                          \
    } while (0)

    CONSIDER_NEXT_HALF(q, dq);
    CONSIDER_NEXT_HALF(rr, dr);
    CONSIDER_NEXT_HALF(s, ds);
#undef CONSIDER_NEXT_HALF

    return (t_next < t_max) ? t_next : t_max;
}

double alea_lattice_next_boundary(const alea_ray_t* ray,
                                    const alea_cell_entry_t* lat_cell,
                                    double t_min,
                                    double t_max) {
    if (lat_cell->lat_type == 1)
        return lattice_rect_next_boundary(ray, lat_cell, t_min, t_max);
    if (lat_cell->lat_type == 2)
        return lattice_hex_next_boundary(ray, lat_cell, t_min, t_max);
    return t_max;
}

typedef enum {
    LATTICE_ENTRY_NO_ENTRY,
    LATTICE_ENTRY_FUTURE_ENTRY,
    LATTICE_ENTRY_ALREADY_INSIDE,
    LATTICE_ENTRY_ERROR
} lattice_entry_kind_t;

typedef struct {
    lattice_entry_kind_t kind;
    double t;
    double sample;
} lattice_entry_query_t;

typedef struct {
    const alea_system_t* sys;
    const alea_ray_t* ray;
    double t_min;
    double t_max;
    double closest;
    alea_raycast_result_t* result;
} lattice_ancestor_event_ctx_t;

static int visit_lattice_ancestor_surface(uint32_t cell_index,
                                          const alea_matrix_t* transform,
                                          void* userdata) {
    lattice_ancestor_event_ctx_t* ctx = userdata;
    if (!ctx || !transform || cell_index >= alea_vec_count(&ctx->sys->cells)) {
        return -1;
    }
    const alea_cell_entry_t* cell = &ctx->sys->cells.data[cell_index];
    alea_ray_t local_ray;
    alea_ray_transform_inverse(transform, ctx->ray, &local_ray);
    int surface_id = -1;
    uint32_t primitive_id = ALEA_PRIMITIVE_ID_INVALID;
    int tested_before = ctx->result ? ctx->result->surfaces_tested : 0;
    double t = raycast_cell_surfaces((alea_system_t*)ctx->sys, &local_ray,
                                     cell, ctx->t_min, ctx->t_max,
                                     ctx->result, &surface_id, &primitive_id);
    if (ctx->result) {
        ctx->result->lattice_entry_ancestor_surface_tests +=
            ctx->result->surfaces_tested - tested_before;
    }
    if (t < ctx->closest) ctx->closest = t;
    return 0;
}

/* Return the next exact enclosing-cell surface for this occurrence.  The
 * support AABB only found the candidate; these surfaces determine when a
 * repeating lattice can actually become active. */
static int lattice_next_ancestor_event(const alea_system_t* sys,
                                       uint32_t placement_index,
                                       const alea_ray_t* ray,
                                       double t_min,
                                       double t_max,
                                       alea_raycast_result_t* result,
                                       double* out_t) {
    if (!out_t) return -1;
    *out_t = DBL_MAX;
    lattice_ancestor_event_ctx_t ctx = {
        .sys = sys,
        .ray = ray,
        .t_min = t_min,
        .t_max = t_max,
        .closest = DBL_MAX,
        .result = result
    };
    int rc = alea_hier_spatial_visit_lattice_placement_ancestors(
        (alea_system_t*)sys, placement_index,
        visit_lattice_ancestor_surface, &ctx);
    if (rc != 0) return -1;
    *out_t = ctx.closest;
    if (result && ctx.closest < DBL_MAX) {
        result->lattice_entry_ancestor_events++;
    }
    return 0;
}

static int lattice_placement_ancestors_active(const alea_system_t* sys,
                                               uint32_t placement_index,
                                               const alea_ray_t* ray,
                                               double t,
                                               double t_max) {
    double sample = t + RAY_EPSILON;
    if (sample >= t_max) return 0;
    double x, y, z;
    alea_ray_point_at(ray, sample, &x, &y, &z);
    return alea_hier_spatial_check_lattice_placement_ancestors(
        (alea_system_t*)sys, placement_index, x, y, z);
}

/* Find the first interval after t_min where both the lattice element and its
 * enclosing occurrence are active.  A repeating lattice is locally valid in
 * every pitch, so while the exact ancestor CSG is false we jump directly to
 * its next surface rather than walking pitches that cannot contribute. */
static lattice_entry_query_t lattice_first_valid_entry(
    const alea_system_t* sys,
    const alea_ray_t* local_ray,
    const alea_ray_t* world_ray,
    const alea_cell_entry_t* cell,
    uint32_t placement_index,
    double t_min,
    double t_max,
    alea_raycast_result_t* result) {
    lattice_entry_query_t out = {
        .kind = LATTICE_ENTRY_NO_ENTRY,
        .t = DBL_MAX,
        .sample = -1.0
    };
    double t_enter, t_exit;
    if (alea_lattice_raycast_interval(local_ray, cell, t_min, t_max,
                                 &t_enter, &t_exit) != 0 ||
        t_exit <= t_enter + RAY_EPSILON) {
        if (result) result->lattice_entry_no_entry_results++;
        return out;
    }

    const int max_steps = alea_lattice_raycast_step_limit(local_ray, cell,
                                                      t_enter, t_exit);
    double t_current = t_enter;
    int ancestor_active = lattice_placement_ancestors_active(
        sys, placement_index, world_ray, t_current, t_exit);
    if (ancestor_active < 0) {
        out.kind = LATTICE_ENTRY_ERROR;
        return out;
    }
    double t_ancestor = DBL_MAX;
    if (lattice_next_ancestor_event(sys, placement_index, world_ray,
                                    t_current + RAY_EPSILON, t_exit,
                                    result, &t_ancestor) != 0) {
        out.kind = LATTICE_ENTRY_ERROR;
        return out;
    }

    for (int step = 0; step < max_steps && t_current < t_exit; step++) {
        /* Outside the occurrence support, lattice periodicity is irrelevant.
         * Move to the next exact ancestor event in one step. */
        if (!ancestor_active) {
            if (t_ancestor >= t_exit - RAY_EPSILON) break;
            t_current = t_ancestor;
            ancestor_active = lattice_placement_ancestors_active(
                sys, placement_index, world_ray, t_current, t_exit);
            if (ancestor_active < 0) {
                out.kind = LATTICE_ENTRY_ERROR;
                return out;
            }
            if (lattice_next_ancestor_event(
                    sys, placement_index, world_ray,
                    t_current + RAY_EPSILON, t_exit, result,
                    &t_ancestor) != 0) {
                out.kind = LATTICE_ENTRY_ERROR;
                return out;
            }
            continue;
        }

        if (result) result->lattice_entry_dda_steps++;
        double t_lattice = alea_lattice_next_boundary(local_ray, cell,
                                                 t_current, t_exit);
        double t_next = fmin(t_lattice, t_ancestor);
        if (t_next <= t_current + RAY_EPSILON || t_next > t_exit)
            t_next = t_exit;
        double x, y, z;
        alea_ray_point_at(local_ray, t_current + 0.5 * (t_next - t_current),
                          &x, &y, &z);
        alea_lattice_location_t location;
        const int valid = alea_lattice_locate_point(sys, cell, x, y, z,
                                                     &location) == 1;
        if (valid) {
            /* CSG containment alone is insufficient in overlapping decks:
             * publishing a synthetic transition for a shadowed fill would
             * make the next resolver sample ambiguous.  Leave the ordinary
             * global surface walker in charge until ownership changes. */
            int canonical =
                alea_hier_spatial_check_lattice_placement_canonical_ancestors(
                    (alea_system_t*)sys, placement_index,
                    world_ray->ox + (t_current + RAY_EPSILON) * world_ray->dx,
                    world_ray->oy + (t_current + RAY_EPSILON) * world_ray->dy,
                    world_ray->oz + (t_current + RAY_EPSILON) * world_ray->dz);
            if (canonical < 0) {
                out.kind = LATTICE_ENTRY_ERROR;
                return out;
            }
            if (!canonical) {
                if (result) result->lattice_entry_canonical_rejections++;
                return out;
            }
            out.t = t_current;
            out.sample = t_current +
                fmin(DDA_SAMPLE_OFFSET, 0.25 * (t_next - t_current));
            if (t_current <= t_min + RAY_EPSILON) {
                out.kind = LATTICE_ENTRY_ALREADY_INSIDE;
                if (result) result->lattice_entry_already_inside_results++;
            } else {
                out.kind = LATTICE_ENTRY_FUTURE_ENTRY;
                if (result) result->lattice_entry_future_entry_results++;
            }
            return out;
        }

        if (t_ancestor <= t_lattice + RAY_EPSILON) {
            ancestor_active = lattice_placement_ancestors_active(
                sys, placement_index, world_ray, t_next, t_exit);
            if (ancestor_active < 0) {
                out.kind = LATTICE_ENTRY_ERROR;
                return out;
            }
            if (lattice_next_ancestor_event(
                    sys, placement_index, world_ray,
                    t_next + RAY_EPSILON, t_exit, result,
                    &t_ancestor) != 0) {
                out.kind = LATTICE_ENTRY_ERROR;
                return out;
            }
        }
        t_current = t_next;
    }
    if (result) result->lattice_entry_no_entry_results++;
    return out;
}

typedef struct {
    const alea_system_t* sys;
    const alea_ray_t* ray;
    double t_min;
    double t_max;
    lattice_entry_query_t closest;
    alea_raycast_result_t* result;
} lattice_entry_visit_ctx_t;

static int visit_lattice_entry(uint32_t placement_index,
                               uint32_t lattice_cell_index,
                               const alea_matrix_t* transform,
                               double placement_enter,
                               double placement_exit,
                               void* userdata) {
    (void)placement_index;
    lattice_entry_visit_ctx_t* ctx = userdata;
    if (!ctx || !transform ||
        lattice_cell_index >= alea_vec_count(&ctx->sys->cells)) {
        return -1;
    }
    const alea_cell_entry_t* cell =
        &ctx->sys->cells.data[lattice_cell_index];
    if (ctx->result) ctx->result->lattice_entry_candidates++;
    alea_ray_t local_ray;
    alea_ray_transform_inverse(transform, ctx->ray, &local_ray);
    lattice_entry_query_t entry = lattice_first_valid_entry(
        ctx->sys, &local_ray, ctx->ray, cell, placement_index,
        fmax(ctx->t_min, placement_enter),
        fmin(ctx->t_max, placement_exit), ctx->result);
    if (entry.kind == LATTICE_ENTRY_ERROR) return -1;
    if (entry.kind != LATTICE_ENTRY_FUTURE_ENTRY ||
        entry.t >= ctx->closest.t) {
        return 0;
    }
    {
        double x, y, z;
        alea_ray_point_at(ctx->ray, entry.sample, &x, &y, &z);
        int valid = alea_hier_spatial_check_lattice_placement_canonical_ancestors(
            (alea_system_t*)ctx->sys, placement_index, x, y, z);
        if (valid < 0) return -1;
        if (!valid) return 0;
        ctx->closest = entry;
    }
    return 0;
}

static int system_first_lattice_entry(alea_system_t* sys,
                                      const alea_ray_t* ray,
                                      double t_min,
                                      double t_max,
                                      alea_raycast_result_t* result,
                                      lattice_entry_query_t* out_entry) {
    if (!out_entry) return -1;
    out_entry->kind = LATTICE_ENTRY_NO_ENTRY;
    out_entry->t = DBL_MAX;
    out_entry->sample = -1.0;
    if (result) result->lattice_entry_calls++;
    alea_hier_lattice_placement_ray_stats_t placement_stats = {0};
    lattice_entry_visit_ctx_t ctx = {
        .sys = sys,
        .ray = ray,
        .t_min = t_min,
        .t_max = t_max,
        .closest = {
            .kind = LATTICE_ENTRY_NO_ENTRY,
            .t = DBL_MAX,
            .sample = -1.0
        },
        .result = result
    };
    if (alea_hier_spatial_visit_lattice_placements_ray(
            sys, ray,
            t_min, t_max, visit_lattice_entry, &placement_stats, &ctx) != 0) {
        return -1;
    }
    if (result) {
        result->lattice_entry_tlas_nodes_tested += placement_stats.nodes_tested;
        result->lattice_entry_tlas_leaves_visited += placement_stats.leaves_visited;
    }
    *out_entry = ctx.closest;
    return 0;
}

/**
 * Find the cell at a point, returning cell index or -1 for void.
 */
static int find_cell_at_point(alea_system_t* sys,
                              double px, double py, double pz,
                              int* out_material_id,
                              double* out_density) {
    alea_cell_hit_t hits[32];
    int num_hits = alea_find_all_cells_at_point(sys, px, py, pz, hits, 32);

    if (num_hits > 0) {
        /* Deepest cell along the first containing chain */
        const alea_cell_hit_t* hit =
            &hits[deepest_first_chain_index(hits, num_hits)];
        if (out_material_id) *out_material_id = hit->material_id;
        if (out_density && hit->cell_index >= 0 &&
            (size_t)hit->cell_index < alea_vec_count(&sys->cells)) {
            *out_density = sys->cells.data[hit->cell_index].density;
        }
        return hit->cell_index;
    }

    if (out_material_id) *out_material_id = 0;
    if (out_density) *out_density = 0;
    return -1;
}

static int find_cell_at_point_hier_path(alea_system_t* sys,
                                        double px, double py, double pz,
                                        int* out_material_id,
                                        double* out_density,
                                        alea_matrix_t* out_transform,
                                        int* out_lattice_cell_index,
                                        alea_matrix_t* out_lattice_transform,
                                        alea_hier_ray_path_t* out_path);

static int find_cell_at_point_hier_path(alea_system_t* sys,
                                        double px, double py, double pz,
                                        int* out_material_id,
                                        double* out_density,
                                        alea_matrix_t* out_transform,
                                        int* out_lattice_cell_index,
                                        alea_matrix_t* out_lattice_transform,
                                        alea_hier_ray_path_t* out_path) {
    alea_hier_cell_hit_t hit_with_transform;
    int found = alea_hier_spatial_find_path_at_point(sys, px, py, pz,
                                                     &hit_with_transform,
                                                     out_path);

    if (found > 0) {
        const alea_cell_hit_t* hit = &hit_with_transform.hit;
        if (out_material_id) *out_material_id = hit->material_id;
        if (out_density && hit->cell_index >= 0 &&
            (size_t)hit->cell_index < alea_vec_count(&sys->cells)) {
            *out_density = sys->cells.data[hit->cell_index].density;
        }
        if (out_transform) *out_transform = hit_with_transform.transform;
        if (out_lattice_cell_index)
            *out_lattice_cell_index = hit_with_transform.lattice_cell_index;
        if (out_lattice_transform)
            *out_lattice_transform = hit_with_transform.lattice_transform;
        return hit->cell_index;
    }

    if (out_material_id) *out_material_id = 0;
    if (out_density) *out_density = 0;
    if (out_transform) alea_matrix_identity(out_transform);
    if (out_lattice_cell_index) *out_lattice_cell_index = -1;
    if (out_lattice_transform) alea_matrix_identity(out_lattice_transform);
    return -1;
}

static int find_cell_from_existing_hier_path(alea_system_t* sys,
                                             const alea_hier_ray_path_t* current_path,
                                             double px, double py, double pz,
                                             int* out_material_id,
                                             double* out_density,
                                             alea_matrix_t* out_transform,
                                             int* out_lattice_cell_index,
                                             alea_matrix_t* out_lattice_transform,
                                             alea_hier_ray_path_t* out_path) {
    if (!sys || !current_path || current_path->count <= 1) return -1;

    for (int parent = current_path->count - 2; parent >= 0; parent--) {
        alea_hier_cell_hit_t hit_with_transform;
        alea_hier_ray_path_t candidate_path;
        const alea_hier_ray_path_entry_t* parent_entry =
            &current_path->entries[parent];
        int found;
        if (parent_entry->is_lattice &&
            (size_t)parent_entry->cell_index < alea_vec_count(&sys->cells)) {
            const alea_cell_entry_t* lattice =
                &sys->cells.data[parent_entry->cell_index];
            double lx = px, ly = py, lz = pz;
            alea_matrix_transform_point_inverse(&parent_entry->transform,
                                                &lx, &ly, &lz);
            alea_lattice_location_t location;
            if (alea_lattice_locate_point(sys, lattice, lx, ly, lz,
                                          &location) != 1) {
                continue;
            }
            found = alea_hier_path_enter_lattice_location(
                sys, current_path, parent, px, py, pz, &location,
                &hit_with_transform, &candidate_path);
        } else {
            /* A terminal-child boundary does not change a still-containing
             * ordinary fill owner.  Preserve that tracked owner rather than
             * issuing a deck-order universe query for every ancestor restart.
             * Lattice transitions above remain canonical because their
             * selected element can change at the boundary. */
            found = alea_hier_spatial_find_path_from_parent_coherent(
                sys, current_path, parent, px, py, pz,
                &hit_with_transform, &candidate_path);
        }
        if (found < 0) return -2;
        if (found == 0) continue;

        const alea_cell_hit_t* hit = &hit_with_transform.hit;
        if (out_material_id) *out_material_id = hit->material_id;
        if (out_density && hit->cell_index >= 0 &&
            (size_t)hit->cell_index < alea_vec_count(&sys->cells)) {
            *out_density = sys->cells.data[hit->cell_index].density;
        }
        if (out_transform) *out_transform = hit_with_transform.transform;
        if (out_lattice_cell_index) {
            *out_lattice_cell_index = hit_with_transform.lattice_cell_index;
        }
        if (out_lattice_transform) {
            *out_lattice_transform = hit_with_transform.lattice_transform;
        }
        if (out_path) *out_path = candidate_path;
        return hit->cell_index;
    }

    return -1;
}

static int find_cell_from_root_universe(alea_system_t* sys,
                                        const alea_ray_root_lookup_t* lookup,
                                        double px, double py, double pz,
                                        int crossed_surface_id,
                                        int* out_material_id,
                                        double* out_density,
                                        alea_matrix_t* out_transform,
                                        int* out_lattice_cell_index,
                                        alea_matrix_t* out_lattice_transform,
                                        alea_hier_ray_path_t* out_path) {
    if (!sys) return -2;

    int root_cell = lookup ? lookup->find(lookup->context, px, py, pz, crossed_surface_id) : -1;
    if (root_cell < 0)
        root_cell = alea_hier_spatial_find_ordered_cell_in_universe(
            sys, 0, px, py, pz, crossed_surface_id);
    if (root_cell < 0) return root_cell == -2 ? -2 : -1;
    if ((size_t)root_cell >= alea_vec_count(&sys->cells)) return -2;

    const alea_cell_entry_t* cell = &sys->cells.data[root_cell];
    alea_hier_ray_path_t root_path;
    root_path.count = 1;
    alea_hier_ray_path_entry_t* ent = &root_path.entries[0];
    ent->cell_index = (uint32_t)root_cell;
    ent->cell_id = cell->mc_cell_id;
    ent->material_id = cell->material_id;
    ent->universe_id = cell->universe_id;
    ent->fill_universe = cell->fill_universe;
    ent->depth = 0;
    ent->is_lattice = 0;
    ent->lat_fill_universe = 0;
    ent->lat_i = 0;
    ent->lat_j = 0;
    ent->lat_k = 0;
    ent->lat_ox = 0.0;
    ent->lat_oy = 0.0;
    ent->lat_oz = 0.0;
    alea_matrix_identity(&ent->transform);

    if (cell->lat_type != 0 && cell->lat_fill) {
        alea_lattice_location_t location;
        if (alea_lattice_locate_point(sys, cell, px, py, pz, &location) != 1)
            return -1;
        ent->fill_universe = location.fill_universe;
        ent->is_lattice = 1;
        ent->lat_fill_universe = location.fill_universe;
        ent->lat_i = location.i;
        ent->lat_j = location.j;
        ent->lat_k = location.k;
        ent->lat_ox = location.ox;
        ent->lat_oy = location.oy;
        ent->lat_oz = location.oz;
    }

    if (cell->fill_universe > 0 || ent->is_lattice) {
        alea_hier_cell_hit_t hit_with_transform;
        alea_hier_ray_path_t candidate_path;
        int found = alea_hier_spatial_find_path_from_parent(
            sys, &root_path, 0, px, py, pz,
            &hit_with_transform, &candidate_path);
        if (found < 0) return -2;
        if (found > 0) {
            const alea_cell_hit_t* hit = &hit_with_transform.hit;
            if (out_material_id) *out_material_id = hit->material_id;
            if (out_density && hit->cell_index >= 0 &&
                (size_t)hit->cell_index < alea_vec_count(&sys->cells)) {
                *out_density = sys->cells.data[hit->cell_index].density;
            }
            if (out_transform) *out_transform = hit_with_transform.transform;
            if (out_lattice_cell_index) {
                *out_lattice_cell_index = hit_with_transform.lattice_cell_index;
            }
            if (out_lattice_transform) {
                *out_lattice_transform = hit_with_transform.lattice_transform;
            }
            if (out_path) *out_path = candidate_path;
            return hit->cell_index;
        }
    }

    if (out_material_id) *out_material_id = cell->material_id;
    if (out_density) *out_density = cell->density;
    if (out_transform) alea_matrix_identity(out_transform);
    if (out_lattice_cell_index) *out_lattice_cell_index = -1;
    if (out_lattice_transform) alea_matrix_identity(out_lattice_transform);
    if (out_path) *out_path = root_path;
    return root_cell;
}

/**
 * Context for finding closest intersection during BVH traversal.
 */
typedef struct {
    alea_system_t* sys;
    const alea_ray_t* ray;
    double t_min;
    double closest_t;
    int closest_surface_id;
    uint32_t closest_primitive_id;
} closest_hit_ctx_t;

/**
 * BVH callback that finds the closest intersection.
 */
static void find_closest_callback(uint32_t surface_idx, void* userdata) {
    closest_hit_ctx_t* ctx = (closest_hit_ctx_t*)userdata;
    const alea_surface_entry_t* surf = &ctx->sys->surfaces.data[surface_idx];
    const alea_primitive_entry_t* prim = &ctx->sys->primitives.data[surf->primitive_id];
    alea_primitive_data_t prim_data;
    if (!raycast_primitive_copy_payload(ctx->sys, surf->primitive_id,
                                        prim->type, &prim_data)) return;

    double t[4];
    int count = ray_intersect_primitive(ctx->ray, prim->type, &prim_data, t);

    for (int j = 0; j < count; j++) {
        if (t[j] > ctx->t_min && t[j] < ctx->closest_t) {
            ctx->closest_t = t[j];
            ctx->closest_surface_id = surf->mc_surface_id;
            ctx->closest_primitive_id = surf->primitive_id;
        }
    }
}

/**
 * Find the closest surface intersection using BVH or linear scan.
 * Return -1 when traversal fails, leaving out_t unpublished.
 */
static int find_closest_intersection(alea_system_t* sys,
                                        const alea_ray_t* ray,
                                        double t_min, double t_max,
                                        int* out_surface_id,
                                        uint32_t* out_primitive_id,
                                        double* out_t) {
    *out_surface_id = -1;
    if (out_primitive_id) *out_primitive_id = ALEA_PRIMITIVE_ID_INVALID;

    /* Use BVH if available */
#if !BVH_DISABLED
    if (sys->surface_bvh) {
        closest_hit_ctx_t ctx = {
            .sys = sys,
            .ray = ray,
            .t_min = t_min,
            .closest_t = t_max,
            .closest_surface_id = -1,
            .closest_primitive_id = ALEA_PRIMITIVE_ID_INVALID
        };
        if (alea_bvh_traverse(sys->surface_bvh, ray, t_min, t_max,
                              find_closest_callback, &ctx) < 0) return -1;
        *out_surface_id = ctx.closest_surface_id;
        if (out_primitive_id) *out_primitive_id = ctx.closest_primitive_id;
        *out_t = ctx.closest_t;
        return 0;
    }
#endif

    /* Fallback: linear scan */
    double closest_t = t_max;
    for (size_t i = 0; i < alea_vec_count(&sys->surfaces); i++) {
        const alea_surface_entry_t* surf = &sys->surfaces.data[i];
        const alea_primitive_entry_t* prim = &sys->primitives.data[surf->primitive_id];
        alea_primitive_data_t prim_data;
        if (!raycast_primitive_copy_payload(sys, surf->primitive_id,
                                            prim->type, &prim_data)) continue;

        double t[4];
        int count = ray_intersect_primitive(ray, prim->type, &prim_data, t);

        for (int j = 0; j < count; j++) {
            if (t[j] > t_min && t[j] < closest_t) {
                closest_t = t[j];
                *out_surface_id = surf->mc_surface_id;
                if (out_primitive_id) *out_primitive_id = surf->primitive_id;
            }
        }
    }
    *out_t = closest_t;
    return 0;
}

static double ray_interval_probe(double enter, double exit) {
    /* DBL_MAX denotes an open ray. A fraction of that sentinel can exceed
     * every model bound even when the next physical step is short. */
    if (exit == DBL_MAX)
        return fmax(enter + 0.381966011250105, nextafter(enter, INFINITY));
    return enter + 0.381966011250105 * (exit - enter);
}

/* A failed boundary-side lookup must be retried inside the very next
 * geometric interval. Probing a fraction of a stale cell's proposed span can
 * jump across other cells (or the whole model). Find only the nearest crossing
 * on this exceptional path, including transformed fills and lattice boundaries.
 * Retaining all later hits would make first-visible rendering depend on the
 * complexity of geometry behind the visible surface. */
#ifdef ALEA_RAY_PROFILE
static int raycast_recovery_sample_profile_impl(alea_system_t* sys, const alea_ray_t* ray,
                                   double t_current, double t_max,
                                   double* out_sample, alea_ray_hit_t* out_hit,
                                   alea_raycast_result_t* owner);
static int raycast_recovery_sample(alea_system_t* sys, const alea_ray_t* ray,
                                   double t_current, double t_max,
                                   double* out_sample, alea_ray_hit_t* out_hit,
                                   alea_raycast_result_t* owner) {
    int previous=alea_ray_profile_enter(RAY_RECOVERY);
    int result=raycast_recovery_sample_profile_impl(sys, ray, t_current, t_max, out_sample, out_hit, owner);
    alea_ray_profile_leave(previous);
    return result;
}
static int raycast_recovery_sample_profile_impl(alea_system_t* sys, const alea_ray_t* ray,
                                   double t_current, double t_max,
                                   double* out_sample, alea_ray_hit_t* out_hit,
                                   alea_raycast_result_t* owner) {
#else
static int raycast_recovery_sample(alea_system_t* sys, const alea_ray_t* ray,
                                   double t_current, double t_max,
                                   double* out_sample, alea_ray_hit_t* out_hit,
                                   alea_raycast_result_t* owner) {
#endif
    alea_raycast_result_t crossings;
    alea_raycast_result_init(&crossings);
    crossings.breakpoint_nearest_only = 1;
    crossings.breakpoint_t_min = t_current;
    recovery_context_t* table=owner->recovery_roots;
    /* Empty-window certificates depend on the accumulating nearest result.
     * Roots do not. Reset certificates before starting a new recovery. */
    if (table) for (size_t i=0;i<table->capacity;i++) {
        table->entries[i].lo=INFINITY;
        table->entries[i].hi=-INFINITY;
    }
    crossings.cancelled = owner->cancelled;
    crossings.cancel_user = owner->cancel_user;
    crossings.recovery_context=table;
    crossings.recovery_context_warmup = table ? 0 : 16;
    crossings.recovery_root_limit = owner->recovery_root_limit;
    const int rc = alea_raycast_validation_breakpoints_reuse_nocache(
        sys, ray, t_current, t_max, 1, &crossings);
    double next = t_max;
    *out_hit = (alea_ray_hit_t){.t = t_max, .surface_id = -1,
                               .primitive_id = ALEA_PRIMITIVE_ID_INVALID};
    if (rc == 0) {
        for (size_t i = 0; i < crossings.hits.count; i++) {
            const double t = crossings.hits.data[i].t;
            if (t > nextafter(t_current, INFINITY) && t < next) {
                next = t;
                *out_hit = crossings.hits.data[i];
            }
        }
    }
    owner->recovery_roots=crossings.recovery_context;
    crossings.recovery_context = NULL;
    alea_raycast_result_free(&crossings);
    if (rc != 0) return -1;
    *out_sample = ray_interval_probe(t_current, next);
    if (!(*out_sample > t_current && *out_sample < next)) {
        alea_set_error_detail(ALEA_ERR_INVALID_STATE,
            "ray boundary recovery has no representable interior sample");
        return -1;
    }
    return 0;
}

#ifdef ALEA_RAY_PROFILE
static void record_owner_probe_recovery(const alea_hier_ray_path_t* path,
        const alea_cell_entry_t* owner, ray_profile_owner_t source, int attempt,
        double current_t, double proposed_t, int proposed_surface, const alea_ray_hit_t* hit) {
    if (!alea_ray_profile.enabled) return;
    if (alea_ray_profile.recovery_count == RAY_PROFILE_RECOVERIES) {
        alea_ray_profile.dropped_recoveries++;
        return;
    }
    ray_profile_recovery_event* e = &alea_ray_profile.recoveries[alea_ray_profile.recovery_count++];
    *e = (ray_profile_recovery_event){
        .ray = alea_ray_profile.ray_ordinal, .attempt = attempt, .source = source,
        .path_length = path ? path->count : 0,
        .root_cell = path && path->count ? path->entries[0].cell_id : -1,
        .owner_cell = owner->mc_cell_id, .owner_universe = owner->universe_id,
        .current_t = current_t, .proposed_t = proposed_t, .proposed_surface = proposed_surface,
        .recovered_t = hit->t, .recovered_surface = hit->surface_id
    };
    if (path) boundary_event_path_keys_at(path, path->count - 1, &e->occurrence, NULL);
}
#endif

/* A midpoint rounded back onto an endpoint in world space cannot establish
 * that the claimed occurrence excludes the open interval. */
static bool interior_excludes_owner(alea_system_t* sys, const alea_ray_t* ray,
                                    const alea_hier_ray_path_t* path, int cell_idx,
                                    double lo, double sample, double hi) {
    if (cell_idx < 0 || !path->count ||
        path->entries[path->count - 1].cell_index != (uint32_t)cell_idx) return false;
    double x, y, z, ax, ay, az, bx, by, bz;
    alea_ray_point_at(ray, sample, &x, &y, &z);
    alea_ray_point_at(ray, lo, &ax, &ay, &az);
    alea_ray_point_at(ray, hi, &bx, &by, &bz);
    if ((x == ax && y == ay && z == az) || (x == bx && y == by && z == bz)) return false;
    return alea_hier_spatial_check_path_containment(
        sys, path, path->count - 1, x, y, z, NULL, NULL, NULL) == 0;
}

typedef struct {
    int cell_index;
    double covered;
} root_owner_span_t;

static int check_root_owner_span(void* context, const alea_ray_coverage_interval_t* interval) {
    root_owner_span_t* span = context;
    if (interval->t_enter != span->covered || interval->kind != ALEA_RAY_COVERAGE_UNIQUE ||
        interval->owner_count != 1 || interval->owners[0].cell_index != span->cell_index)
        return 1;
    span->covered = interval->t_exit;
    return 0;
}

/* Recovery bounds can be roots of a different cell's infinite surfaces.
 * In a terminal world cell, the proposed exit is already the nearest local
 * crossing. Accept that exit only if independent all-owner coverage confirms
 * the same cell throughout the entire extension. */
static int verify_recovered_root_span(alea_system_t* sys, const alea_ray_t* ray,
        int cell_index, int guard_surface, double lo, double hi, size_t budget) {
    const alea_cell_entry_t* cell = &sys->cells.data[cell_index];
    if (cell->universe_id != 0 || alea_cell_entry_is_container(cell) ||
        !cell->surface_indices || !cell->surface_index_count || hi == DBL_MAX ||
        cell_references_surface_id(sys, cell, guard_surface))
        return 0;
    alea_raycast_result_t scratch;
    alea_raycast_result_init(&scratch);
    root_owner_span_t span = {.cell_index = cell_index, .covered = lo};
    int rc = alea_ray_coverage_sweep_strict_reuse_nocache(sys, ray, lo, hi, budget,
            NULL, &scratch, check_root_owner_span, &span);
    alea_raycast_result_free(&scratch);
    if (rc >= 0) return span.covered == hi;
    /* The optional proof can extend beyond the actual requested flight.
     * An inconclusive sweep keeps the guard, allowing shorter advances. */
    alea_error_t error = alea_get_last_error();
    if (error != ALEA_ERR_INVALID_STATE && error != ALEA_ERR_OVERFLOW) return -1;
    alea_clear_error_detail();
    return 0;
}

static int raycast_cell_aware_resume(alea_system_t* sys,
                                     const alea_ray_t* ray,
                                     double effective_t_max,
                                     bool use_hier_lookup,
                                     bool emit_hits,
                                     bool enrich_local_group,
                                     alea_raycast_result_t* result,
                                     alea_ray_walk_t* state,
                                     int step_budget,
                                     alea_ray_selected_interval_t* out_selected) {
    if (!state || step_budget <= 0) return -1;
    /* Current position along ray */
    double t_current = state->t_current;
    int prev_cell_idx = state->prev_cell_idx;
    int prev_surface_id = state->prev_surface_id;
    /* Hit index of the surface crossed at t_current (the next segment's
     * enter boundary), or -1. Only maintained when emit_hits is set. */
    int pending_enter_hit_index = state->pending_enter_hit_index;
    /* The physical or synthetic event that entered the current interval.
     * FIRST_VISIBLE reports an entering surface, so it must not reuse the
     * current iteration's exit event. */
    alea_raycast_boundary_event_t enter_event = state->enter_event;
    /* The lane owns the authoritative path.  Operating on it in place avoids
     * copying the fixed 64-entry path into and out of a local value at every
     * one-interval yield.  The retry path below still takes an explicit
     * snapshot only when midpoint verification actually re-resolves. */
    alea_hier_ray_path_t* current_path = &state->current_path;
    double pending_lattice_entry_sample = state->pending_lattice_entry_sample;
    const bool need_boundary_event = emit_hits || out_selected;

    while (t_current < effective_t_max && step_budget-- > 0 &&
           state->iterations_remaining-- > 0) {
        if (alea_ray_interrupted(result)) {
            alea_set_error_detail(ALEA_ERR_INTERRUPTED, "ray traversal interrupted");
            return -1;
        }
        result->step_iterations++;
        /* Cell resolution samples just past the crossing at t_current, where
         * quadric sign evaluation is numerically noisy. Midpoint verification
         * after the crossing search jumps back here once with a sample point
         * in the segment interior when the resolved cell fails to contain it. */
        double t_sample = pending_lattice_entry_sample > t_current
            ? pending_lattice_entry_sample
            : raycast_forward_sample_t(ray, t_current);
        pending_lattice_entry_sample = -1.0;
        bool incident_lattice = current_path->count >= 2 &&
            current_path->entries[current_path->count - 2].is_lattice;
        int resolve_attempt = 0;
        alea_ray_hit_t recovery_hit = {.t = effective_t_max};
        /* Preserve the initial path in case recovery finds only a container. */
        int saved_valid = 0;
        int saved_cell_idx = -1;
        int saved_cell_id = -1;
        int saved_material_id = 0;
        double saved_density = 0;
        double saved_t_next = 0;
        int saved_hit_surface_id = -1;
        int saved_next_enter_surface_id = -1;
        alea_raycast_boundary_event_t saved_bevent;
        alea_hier_ray_path_t saved_path;
        memset(&saved_bevent, 0, sizeof(saved_bevent));
        saved_path.count = 0;
resolve_cell:;
#ifdef ALEA_RAY_PROFILE
        uint64_t path_hits_before = result->owner_path_hits;
        uint64_t root_hits_before = result->owner_root_hits;
        uint64_t full_hits_before = result->owner_full_hits;
#endif
        int cell_idx = -1;
        int cell_id = -1;
        int material_id = 0;
        double density = 0;
        int found_via_neighbor = 0;
        alea_matrix_t cell_transform;
        alea_matrix_t lattice_transform;
        int lattice_cell_index = -1;
        alea_matrix_identity(&cell_transform);
        alea_matrix_identity(&lattice_transform);

        /* Try neighbor lookup first if we have previous cell and surface info */
        if (prev_cell_idx >= 0 && prev_surface_id > 0) {
            result->owner_neighbor_attempts++;
            const alea_cell_entry_t* prev_cell = &sys->cells.data[prev_cell_idx];
            if (!use_hier_lookup || prev_cell->universe_id == 0) {
                double px, py, pz;
                alea_ray_point_at(ray, t_sample, &px, &py, &pz);
                found_via_neighbor = raycast_find_neighbor(sys, prev_cell_idx,
                                                           prev_surface_id,
                                                           px, py, pz,
                                                           &cell_id, &cell_idx,
                                                           &material_id,
                                                           &density);
                /* In hier mode a filled container is not a terminal answer:
                 * reject it so the full lookup descends into the fill. */
                if (found_via_neighbor && use_hier_lookup) {
                    const alea_cell_entry_t* nb = &sys->cells.data[cell_idx];
                    if (nb->fill_universe > 0 ||
                        (nb->lat_type != 0 && nb->lat_fill)) {
                        found_via_neighbor = 0;
                        cell_idx = -1;
                        cell_id = -1;
                        material_id = 0;
                        density = 0;
                    } else if (nb->universe_id == 0) {
                        /* Neighbor resolution at the root bypasses the full
                         * hierarchy lookup.  Refresh its path explicitly:
                         * retaining the previous root entry would attach the
                         * wrong owner occurrence to this segment. */
                        current_path->count = 1;
                        alea_hier_ray_path_entry_t* entry =
                            &current_path->entries[0];
                        entry->cell_index = (uint32_t)cell_idx;
                        entry->cell_id = nb->mc_cell_id;
                        entry->material_id = nb->material_id;
                        entry->universe_id = nb->universe_id;
                        entry->fill_universe = nb->fill_universe;
                        entry->depth = 0;
                        entry->is_lattice = 0;
                        entry->lat_fill_universe = 0;
                        entry->lat_i = 0;
                        entry->lat_j = 0;
                        entry->lat_k = 0;
                        entry->lat_ox = 0.0;
                        entry->lat_oy = 0.0;
                        entry->lat_oz = 0.0;
                        alea_matrix_identity(&entry->transform);
                    }
                }
            }
            if (found_via_neighbor) result->owner_neighbor_hits++;
        }

        /* In hierarchical mode, avoid a full root-to-deepest lookup when the
         * sampled point remains inside the cached placement path. */
        if (!found_via_neighbor && use_hier_lookup && prev_cell_idx >= 0 &&
            current_path->count > 0 &&
            current_path->entries[current_path->count - 1].cell_index ==
                (uint32_t)prev_cell_idx) {
            result->owner_path_attempts++;
            const alea_cell_entry_t* prev_cell = &sys->cells.data[prev_cell_idx];
            double px, py, pz;
            alea_ray_point_at(ray, t_sample, &px, &py, &pz);
            int in_cell = 0;
            if (prev_cell->universe_id == 0) {
                in_cell = prev_cell->root_node_id != ALEA_NODE_ID_INVALID &&
                    alea_contains_point(sys, prev_cell->root_node_id,
                                        px, py, pz);
                alea_matrix_identity(&cell_transform);
                lattice_cell_index = -1;
                alea_matrix_identity(&lattice_transform);
            } else {
                in_cell = alea_hier_spatial_check_path_containment(
                    sys, current_path, current_path->count - 1, px, py, pz,
                    &cell_transform, &lattice_cell_index, &lattice_transform);
                if (in_cell < 0) in_cell = 0;
            }
            if (in_cell && !alea_cell_entry_is_container(prev_cell)) {
                cell_idx = prev_cell_idx;
                cell_id = prev_cell->mc_cell_id;
                material_id = prev_cell->material_id;
                density = prev_cell->density;
                found_via_neighbor = 1;
                result->owner_path_hits++;
            }
        }

        /* Fall back to full lookup if neighbor lookup failed or not available */
        if (!found_via_neighbor) {
            double px, py, pz;
            alea_ray_point_at(ray, t_sample, &px, &py, &pz);
            if (use_hier_lookup) {
                cell_idx = find_cell_from_existing_hier_path(
                    sys, current_path, px, py, pz,
                    &material_id, &density,
                    &cell_transform, &lattice_cell_index,
                    &lattice_transform, current_path);
                if (cell_idx == -2) return -1;
            }
            if (use_hier_lookup && cell_idx < 0) {
                result->owner_root_queries++;
                cell_idx = find_cell_from_root_universe(
                    sys, state->root_lookup, px, py, pz, prev_surface_id,
                    &material_id, &density,
                    &cell_transform, &lattice_cell_index,
                    &lattice_transform, current_path);
                if (cell_idx == -2) return -1;
                if (cell_idx >= 0) result->owner_root_hits++;
            }
            if (cell_idx < 0) {
                result->point_lookups++;
                if (use_hier_lookup) result->owner_full_queries++;
                cell_idx = use_hier_lookup
                    ? find_cell_at_point_hier_path(
                        sys, px, py, pz, &material_id, &density,
                        &cell_transform, &lattice_cell_index,
                        &lattice_transform, current_path)
                    : find_cell_at_point(sys, px, py, pz,
                                         &material_id, &density);
                if (use_hier_lookup && cell_idx >= 0)
                    result->owner_full_hits++;
                if (use_hier_lookup && cell_idx < 0) current_path->count = 0;
            }
            if (cell_idx >= 0 && (size_t)cell_idx < alea_vec_count(&sys->cells)) {
                cell_id = sys->cells.data[cell_idx].mc_cell_id;
            }
        }

        /* Find next surface crossing:
         * Use per-cell surface index when inside a cell (the whole point of cell-aware),
         * fall back to global search only for void regions. */
        int hit_surface_id = -1;
        int next_enter_surface_id = -1;
        double t_next;
        double t_lattice_next = DBL_MAX;
        /* Winning boundary event for this step. Only populated/consumed when
         * emit_hits is set; the segment-only path never touches it, so it adds
         * no work to the hot loop. */
        alea_raycast_boundary_event_t bevent;
        if (need_boundary_event) {
            memset(&bevent, 0, sizeof(bevent));
            bevent.t = 0;
            bevent.surface_id = -1;
            bevent.primitive_id = ALEA_PRIMITIVE_ID_INVALID;
            bevent.active_cell_id = -1;
            bevent.active_universe_id = -1;
            bevent.active_depth = -1;
            bevent.has_physical_surface = false;
            bevent.is_synthetic_lattice_boundary = false;
            alea_matrix_identity(&bevent.transform);
        }
        bool can_prove_interval = !resolve_attempt && use_hier_lookup && need_boundary_event &&
            enter_event.has_physical_surface && !state->fail_unverified_ownership;
        const double crossing_t_min = nextafter(t_current, INFINITY);
        int adjacent_surface = -1;
        uint32_t adjacent_primitive = ALEA_PRIMITIVE_ID_INVALID;
        bool in_lattice_element = current_path->count >= 2 &&
            current_path->entries[current_path->count - 2].is_lattice;
        bool left_lattice = incident_lattice && !in_lattice_element;
        bool containment_unverified = false;
        if (cell_idx >= 0 && (size_t)cell_idx < alea_vec_count(&sys->cells) &&
            sys->cells.data[cell_idx].surface_indices) {
            alea_ray_t local_ray;
            const alea_ray_t* surface_ray = ray;
            if (use_hier_lookup) {
                alea_ray_transform_inverse(&cell_transform, ray, &local_ray);
                surface_ray = &local_ray;
            }
            /* Instrumentation: attribute this cell's surface tests and record
             * its size. Snapshotting surfaces_tested keeps raycast_cell_surfaces
             * free of a category parameter. */
            uint32_t terminal_surf_count =
                sys->cells.data[cell_idx].surface_index_count;
            result->crossed_cell_count++;
            result->sum_cell_surface_count += terminal_surf_count;
            if (terminal_surf_count > result->max_cell_surface_count)
                result->max_cell_surface_count = terminal_surf_count;
            int surfaces_tested_before = result->surfaces_tested;

            uint32_t terminal_prim_id = ALEA_PRIMITIVE_ID_INVALID;
            bool check_adjacent = can_prove_interval && t_sample >= crossing_t_min &&
                ((lattice_cell_index < 0 && current_path->count == 2) || in_lattice_element || incident_lattice);
            t_next = raycast_cell_surfaces(sys, surface_ray,
                                           &sys->cells.data[cell_idx],
                                           check_adjacent ? t_current : crossing_t_min, effective_t_max,
                                           result,
                                           &hit_surface_id,
                                           &terminal_prim_id);
            if (check_adjacent && t_next == crossing_t_min) {
                /* Keep the ordinary next hit available if the two crossings
                 * cannot be distinguished in their placement. */
                adjacent_surface = hit_surface_id;
                adjacent_primitive = terminal_prim_id;
                t_next = raycast_cell_surfaces(sys, surface_ray, &sys->cells.data[cell_idx],
                    crossing_t_min, effective_t_max, result, &hit_surface_id, &terminal_prim_id);
            }
            result->terminal_surfaces_tested +=
                result->surfaces_tested - surfaces_tested_before;
            if (need_boundary_event) {
                bevent.t = t_next;
                bevent.surface_id = hit_surface_id;
                bevent.primitive_id = terminal_prim_id;
                bevent.has_physical_surface = (hit_surface_id > 0);
                bevent.is_synthetic_lattice_boundary = false;
                bevent.transform = cell_transform;
            }
            if (use_hier_lookup &&
                lattice_cell_index >= 0 &&
                (size_t)lattice_cell_index < alea_vec_count(&sys->cells)) {
                const alea_cell_entry_t* lattice_cell =
                    &sys->cells.data[lattice_cell_index];
                alea_ray_t lattice_ray;
                alea_ray_transform_inverse(&lattice_transform, ray, &lattice_ray);
                int lattice_surface_id = -1;
                uint32_t lattice_prim_id = ALEA_PRIMITIVE_ID_INVALID;
                int lattice_tested_before = result->surfaces_tested;
                double t_lattice_surface =
                    raycast_cell_surfaces(sys, &lattice_ray, lattice_cell,
                                          crossing_t_min,
                                          effective_t_max,
                                          result,
                                          &lattice_surface_id,
                                          &lattice_prim_id);
                result->lattice_surfaces_tested +=
                    result->surfaces_tested - lattice_tested_before;
                t_lattice_next = alea_lattice_next_boundary(&lattice_ray, lattice_cell,
                                                       t_current + RAY_EPSILON,
                                                       effective_t_max);
                if (t_lattice_surface < t_next - RAY_EPSILON) {
                    t_next = t_lattice_surface;
                    hit_surface_id = lattice_surface_id;
                    if (fabs(t_lattice_surface - t_lattice_next) <=
                        SURFACE_SAMPLE_OFFSET) {
                        const double post_t =
                            fmax(t_lattice_surface, t_lattice_next) +
                            10.0 * RAY_EPSILON;
                        double post_x, post_y, post_z;
                        alea_ray_point_at(&lattice_ray, post_t,
                                          &post_x, &post_y, &post_z);
                        alea_lattice_location_t post_location;
                        if (alea_lattice_locate_point(
                                sys, lattice_cell, post_x, post_y, post_z,
                                &post_location) != 1) {
                            pending_lattice_entry_sample = post_t;
                        }
                    }
                    if (need_boundary_event) {
                        bevent.t = t_lattice_surface;
                        bevent.surface_id = lattice_surface_id;
                        bevent.primitive_id = lattice_prim_id;
                        bevent.has_physical_surface = (lattice_surface_id > 0);
                        bevent.is_synthetic_lattice_boundary = false;
                        bevent.transform = lattice_transform;
                    }
                }
                if (t_lattice_next < effective_t_max - RAY_EPSILON &&
                    t_lattice_next < t_next - RAY_EPSILON) {
                    t_next = t_lattice_next;
                    if (fabs(t_lattice_surface - t_lattice_next) <= RAY_EPSILON) {
                        hit_surface_id = lattice_surface_id;
                        next_enter_surface_id = 0;
                        if (need_boundary_event) {
                            bevent.t = t_lattice_next;
                            bevent.surface_id = lattice_surface_id;
                            bevent.primitive_id = lattice_prim_id;
                            bevent.has_physical_surface = (lattice_surface_id > 0);
                            bevent.is_synthetic_lattice_boundary = true;
                            bevent.transform = lattice_transform;
                        }
                    } else {
                        hit_surface_id = 0;
                        if (need_boundary_event) {
                            /* Pure synthetic DDA boundary: no physical surface. */
                            bevent.t = t_lattice_next;
                            bevent.surface_id = 0;
                            bevent.primitive_id = ALEA_PRIMITIVE_ID_INVALID;
                            bevent.has_physical_surface = false;
                            bevent.is_synthetic_lattice_boundary = true;
                        }
                    }
                } else if (fabs(t_lattice_surface - t_lattice_next) <= RAY_EPSILON &&
                           fabs(t_lattice_next - t_next) <= RAY_EPSILON &&
                           lattice_surface_id >= 0) {
                    hit_surface_id = lattice_surface_id;
                    next_enter_surface_id = 0;
                    if (need_boundary_event) {
                        bevent.surface_id = lattice_surface_id;
                        bevent.primitive_id = lattice_prim_id;
                        bevent.has_physical_surface = (lattice_surface_id > 0);
                        bevent.is_synthetic_lattice_boundary = true;
                        bevent.transform = lattice_transform;
                    }
                }
            }
            if (use_hier_lookup && current_path->count > 1) {
                int ancestor_surface_id = -1;
                uint32_t ancestor_prim_id = ALEA_PRIMITIVE_ID_INVALID;
                alea_matrix_t ancestor_transform;
                alea_matrix_identity(&ancestor_transform);
                int ancestor_tested_before = result->surfaces_tested;
                double t_ancestor = raycast_hier_path_ancestor_surfaces(
                    sys, ray, current_path, (uint32_t)cell_idx,
                    lattice_cell_index, crossing_t_min,
                    effective_t_max, t_current,
                    state->boundary_distance_tolerance,
                    result, &ancestor_surface_id,
                    &ancestor_prim_id, &ancestor_transform);
                result->ancestor_surfaces_tested +=
                    result->surfaces_tested - ancestor_tested_before;
                const bool same_near_boundary = need_boundary_event &&
                    ray_crossing_distances_near(
                        t_current, t_ancestor, t_next,
                        state->boundary_distance_tolerance) &&
                    ray_boundaries_same_world_surface(
                        ancestor_prim_id, &ancestor_transform,
                        bevent.primitive_id, &bevent.transform);
                if (t_ancestor < t_next - RAY_EPSILON ||
                    same_near_boundary) {
                    t_next = t_ancestor;
                    hit_surface_id = ancestor_surface_id;
                    next_enter_surface_id = ancestor_surface_id;
                    if (need_boundary_event) {
                        bevent.t = t_ancestor;
                        bevent.surface_id = ancestor_surface_id;
                        bevent.primitive_id = ancestor_prim_id;
                        bevent.has_physical_surface = (ancestor_surface_id > 0);
                        bevent.is_synthetic_lattice_boundary = false;
                        bevent.transform = ancestor_transform;
                    }
                }
            }
        } else {
            uint32_t void_prim_id = ALEA_PRIMITIVE_ID_INVALID;
            if (find_closest_intersection(sys, ray,
                    nextafter(t_current, INFINITY), effective_t_max,
                    &hit_surface_id, &void_prim_id, &t_next) != 0) return -1;
            if (need_boundary_event) {
                /* Void-region global search runs in the world frame. */
                bevent.t = t_next;
                bevent.surface_id = hit_surface_id;
                bevent.primitive_id = void_prim_id;
                bevent.has_physical_surface = (hit_surface_id > 0);
                bevent.is_synthetic_lattice_boundary = false;
                alea_matrix_identity(&bevent.transform);
            }
        }

        if (use_hier_lookup && cell_idx >= 0 &&
            alea_cell_entry_is_container(&sys->cells.data[cell_idx])) {
            alea_raycast_boundary_event_t child_event = {0};
            if (need_boundary_event) child_event = bevent;
            int child_crossing = raycast_undefined_fill_crossing(
                sys, ray, current_path, nextafter(t_current, INFINITY),
                fmin(t_next, effective_t_max), result, &child_event);
            if (child_crossing < 0) return -1;
            if (child_crossing > 0) {
                t_next = child_event.t;
                hit_surface_id = child_event.surface_id;
                next_enter_surface_id = child_event.surface_id;
                if (need_boundary_event) bevent = child_event;
            }
        }

        if (use_hier_lookup && lattice_cell_index < 0 &&
            (cell_idx < 0 ||
             ((size_t)cell_idx < alea_vec_count(&sys->cells) &&
              alea_cell_entry_is_container(&sys->cells.data[cell_idx])))) {
            lattice_entry_query_t lattice_entry;
            if (system_first_lattice_entry(
                sys, ray, t_current + RAY_EPSILON, effective_t_max,
                result, &lattice_entry) != 0) {
                return -1;
            }
            if (lattice_entry.kind == LATTICE_ENTRY_FUTURE_ENTRY &&
                lattice_entry.t < t_next - SURFACE_SAMPLE_OFFSET) {
                t_next = lattice_entry.t;
                hit_surface_id = 0;
                next_enter_surface_id = 0;
                pending_lattice_entry_sample = lattice_entry.sample;
                if (need_boundary_event) {
                    bevent.t = lattice_entry.t;
                    bevent.surface_id = 0;
                    bevent.primitive_id = ALEA_PRIMITIVE_ID_INVALID;
                    bevent.has_physical_surface = false;
                    bevent.is_synthetic_lattice_boundary = true;
                }
            }
        }

#ifdef ALEA_RAY_PROFILE
        /* Capture this attempt's local proposal before applying a prior
         * recovery bound. The open endpoint is represented by surface -1. */
        double proposed_t = fmin(t_next, effective_t_max);
        int proposed_surface = t_next < effective_t_max ? hit_surface_id : -1;
#endif
        /* A recovery probe is valid only up to its independently enumerated
         * next crossing. In particular a lattice DDA step must not overrun a
         * nearby physical ancestor boundary. */
        bool clip_to_recovery = resolve_attempt && recovery_hit.t < t_next;
        if (clip_to_recovery && state->recovery_proof_budget && recovery_hit.surface_id > 0 &&
            cell_idx >= 0 && current_path->count == 1 && bevent.has_physical_surface) {
            int verified = verify_recovered_root_span(sys, ray, cell_idx, recovery_hit.surface_id,
                t_current, t_next, state->recovery_proof_budget);
            if (verified < 0) return -1;
            if (verified) clip_to_recovery = false;
        }
        if (clip_to_recovery) {
            t_next = recovery_hit.t;
            hit_surface_id = recovery_hit.surface_id;
            next_enter_surface_id = 0;
            if (need_boundary_event) {
                bevent.t = t_next;
                bevent.surface_id = hit_surface_id;
                bevent.primitive_id = recovery_hit.primitive_id;
                bevent.has_physical_surface = false;
                bevent.is_synthetic_lattice_boundary = hit_surface_id == 0;
                for (int p = current_path->count - 1; p >= 0; p--) {
                    if (boundary_event_cell_uses_exact_surface(sys,
                            current_path->entries[p].cell_index, hit_surface_id)) {
                        bevent.transform = current_path->entries[p].transform;
                        bevent.has_physical_surface = hit_surface_id > 0;
                        break;
                    }
                }
            }
        }

        /* Never skip geometry to manufacture progress.  Every accepted
         * intersection is strictly beyond t_current; failure here indicates
         * a numerical or traversal defect that callers must see. */
        if (!(t_next > t_current)) {
            alea_set_error_detail(
                ALEA_ERR_INVALID_STATE,
                "ray traversal failed to advance at t=%.17g", t_current);
            return -1;
        }

        /* raycast_cell_surfaces returns DBL_MAX when no hit lies before t_max;
         * clamp so the terminal segment stops at the user's max_distance. */
        if (t_next > effective_t_max) t_next = effective_t_max;

        /* A forward ownership sample can overrun a very thin interval.
         * Try to prove ownership between the known crossings before retrying
         * at an interior sample. Keep the original entry distance. */
        bool slab_verified = false;
        if (adjacent_surface > 0 && crossing_t_min < t_next &&
            alea_ray_resolve_short_interval(sys, ray, enter_event.primitive_id, &enter_event.transform,
                adjacent_primitive, &cell_transform, t_current, crossing_t_min, true, left_lattice, current_path)) {
            slab_verified = true;
            t_next = crossing_t_min;
            hit_surface_id = next_enter_surface_id = adjacent_surface;
            bevent = (alea_raycast_boundary_event_t){.t = t_next, .surface_id = hit_surface_id,
                .primitive_id = adjacent_primitive, .transform = cell_transform, .has_physical_surface = true,
                .active_cell_id = -1, .active_depth = -1, .active_universe_id = -1};
        } else if (can_prove_interval && t_sample >= t_next && bevent.has_physical_surface) {
            slab_verified = alea_ray_resolve_short_interval(sys, ray, enter_event.primitive_id,
                &enter_event.transform, bevent.primitive_id, &bevent.transform, t_current, t_next,
                false, left_lattice, current_path);
        }
        if (slab_verified) {
            cell_idx = (int)current_path->entries[current_path->count - 1].cell_index;
            const alea_cell_entry_t* owner = &sys->cells.data[cell_idx];
            cell_id = owner->mc_cell_id;
            material_id = owner->material_id;
            density = owner->density;
        } else if (t_sample >= t_next) {
            const double sample = t_current + 0.5 * (t_next - t_current);
            const bool strict = state->require_interior_sample || state->fail_unverified_ownership;
            if (resolve_attempt < 2 && sample > t_current && sample < t_next &&
                (strict || (use_hier_lookup && !resolve_attempt &&
                            interior_excludes_owner(sys, ray, current_path, cell_idx,
                                                    t_current, sample, t_next)))) {
                recovery_hit.t = t_next;
                recovery_hit.surface_id = hit_surface_id;
                recovery_hit.primitive_id = need_boundary_event ? bevent.primitive_id : ALEA_PRIMITIVE_ID_INVALID;
                t_sample = sample;
                resolve_attempt++;
                goto resolve_cell;
            }
            if (strict) {
                alea_set_error_detail(ALEA_ERR_INVALID_STATE,
                    "ray crossing has no stable interior sample");
                return -1;
            }
        }

        /* Verify the resolved cell actually owns this segment. No surface of
         * the chosen cell lies inside the open interval, so containment is
         * constant across it: probing an interior point — well away from the
         * numerically noisy boundary sample — catches both bogus adjacency
         * transitions (non-convex cells re-entered across their own surface)
         * and on-boundary sign noise picking a coincident sibling cell. On
         * failure, redo the resolution once, sampling at the probe point. */
        if (!slab_verified &&
            cell_idx >= 0 && (size_t)cell_idx < alea_vec_count(&sys->cells) &&
            (t_next - t_current > 1e-6 || resolve_attempt > 0)) {
            /* Probe at an irrational fraction of the interval, not the exact
             * midpoint: a segment spanning lattice elements has periodic
             * internal boundaries, and its midpoint can land exactly on one,
             * where containment is numerically undefined. */
            double t_probe = ray_interval_probe(t_current, t_next);
            double mx, my, mz;
            alea_ray_point_at(ray, t_probe, &mx, &my, &mz);
            const alea_cell_entry_t* seg_cell = &sys->cells.data[cell_idx];
            int probe_ok;
            if (!use_hier_lookup && seg_cell->universe_id != 0) {
                /* Flat mode resolves nested cells through the flat spatial
                 * query; their CSG lives in a local frame we don't track
                 * here, so containment against world coordinates would be
                 * meaningless. Trust the resolution. */
                probe_ok = 1;
            } else if (!use_hier_lookup || seg_cell->universe_id == 0) {
                probe_ok = seg_cell->root_node_id != ALEA_NODE_ID_INVALID &&
                    alea_contains_point(sys, seg_cell->root_node_id,
                                        mx, my, mz);
            } else {
                probe_ok = alea_hier_spatial_check_path_containment(
                    sys, current_path, current_path->count - 1, mx, my, mz,
                    NULL, NULL, NULL) > 0;
            }
            /* Leaving a fill can lead directly to a shallower terminal cell.
             * Accept it only when the new path verifies its ownership. */
            bool unverified_ascent = saved_valid && use_hier_lookup &&
                current_path->count < saved_path.count &&
                (!probe_ok || alea_cell_entry_is_container(seg_cell));
            if (!probe_ok || unverified_ascent) {
                if (resolve_attempt < 2 && !probe_ok && !unverified_ascent) {
                    /* Re-resolve inside the next geometric interval, whose
                     * crossings are independent of the possibly stale owner. */
                    if (!saved_valid) {
                        saved_valid = 1;
                        saved_cell_idx = cell_idx;
                        saved_cell_id = cell_id;
                        saved_material_id = material_id;
                        saved_density = density;
                        saved_t_next = t_next;
                        saved_hit_surface_id = hit_surface_id;
                        saved_next_enter_surface_id = next_enter_surface_id;
                        if (need_boundary_event) saved_bevent = bevent;
                        ray_path_copy_live(&saved_path, current_path);
                        result->path_snapshot_copies++;
                        result->path_snapshot_entries +=
                            (uint64_t)saved_path.count;
                    }
                    resolve_attempt++;
                    /* The active path already supplies an upper bound. A
                     * broader search may miss this placement's surface and
                     * select a point beyond the interval being repaired. */
#ifdef ALEA_RAY_PROFILE
                    double recovery_start = alea_ray_profile.enabled ? alea_monotonic_seconds() : 0;
#endif
                    int recovery_rc = raycast_recovery_sample(sys, ray, t_current,
                            t_next, &t_sample, &recovery_hit, result);
#ifdef ALEA_RAY_PROFILE
                    if (alea_ray_profile.enabled) {
                        double elapsed = alea_monotonic_seconds() - recovery_start;
                        const alea_universe_t* u = alea_get_universe(sys, seg_cell->universe_id);
                        if (u) {
                            ray_profile_universe* row = &alea_ray_profile.universes[u - sys->universes.data];
                            row->recoveries++;
                            row->recovery_seconds += elapsed;
                        }
                    }
#endif
                    if (recovery_rc != 0) return -1;
                    if (recovery_hit.t == t_next) {
                        recovery_hit.surface_id = hit_surface_id;
                        if (need_boundary_event)
                            recovery_hit.primitive_id = bevent.primitive_id;
                    }
#ifdef ALEA_RAY_PROFILE
                    ray_profile_owner_t source = RAY_OWNER_PATH_DESCENT;
                    if (result->owner_path_hits != path_hits_before) source = RAY_OWNER_CACHED_PATH;
                    else if (found_via_neighbor) source = RAY_OWNER_NEIGHBOR;
                    else if (!use_hier_lookup || result->owner_full_hits != full_hits_before) source = RAY_OWNER_FULL;
                    else if (result->owner_root_hits != root_hits_before) source = RAY_OWNER_ROOT;
                    record_owner_probe_recovery(use_hier_lookup ? current_path : NULL, seg_cell,
                        source, resolve_attempt, t_current, proposed_t, proposed_surface, &recovery_hit);
#endif
                    goto resolve_cell;
                }
                if (unverified_ascent) {
                    cell_idx = saved_cell_idx;
                    cell_id = saved_cell_id;
                    material_id = saved_material_id;
                    density = saved_density;
                    t_next = saved_t_next;
                    hit_surface_id = saved_hit_surface_id;
                    next_enter_surface_id = saved_next_enter_surface_id;
                    bevent = saved_bevent;
                    ray_path_copy_live(current_path, &saved_path);
                }
                containment_unverified = true;
            }
        } else if (resolve_attempt >= 1 && saved_valid && cell_idx >= 0) {
            /* A degenerate owned interval still cannot establish ownership.
             * Exterior void, however, is a legitimate result of the bounded
             * recovery probe and must not restore the cell just exited. */
            cell_idx = saved_cell_idx;
            cell_id = saved_cell_id;
            material_id = saved_material_id;
            density = saved_density;
            t_next = saved_t_next;
            hit_surface_id = saved_hit_surface_id;
            next_enter_surface_id = saved_next_enter_surface_id;
            bevent = saved_bevent;
            ray_path_copy_live(current_path, &saved_path);
            containment_unverified = true;
        }

        if (containment_unverified && state->fail_unverified_ownership) {
            alea_set_error_detail(ALEA_ERR_INVALID_STATE,
                "ray traversal could not verify cell ownership at t=%.17g",
                t_current);
            return -1;
        }

        bevent.t = t_next;
        if (need_boundary_event && enrich_local_group)
            boundary_event_collect_local_group(sys, ray, current_path, &bevent);
        alea_ray_selected_interval_t selected = {
            .t_enter = t_current,
            .containment_verified = !containment_unverified,
            .t_exit = t_next,
            .cell_index = cell_idx,
            .cell_id = cell_id,
            .material_id = material_id,
            .density = density,
            .enter_surface_id = prev_surface_id,
            .exit_surface_id = hit_surface_id,
            .enter_event = enter_event,
            .exit_event = bevent,
            .path = current_path,
            .resolution_flags =
                (cell_idx >= 0 &&
                 (size_t)cell_idx < alea_vec_count(&sys->cells) &&
                 alea_cell_entry_is_container(&sys->cells.data[cell_idx]))
                    ? ALEA_RESOLVE_UNDEFINED_FILL : 0
        };
        boundary_interval_capture_owner(current_path, &selected);

        /* The selected-interval facade consumes the geometric result before
         * any compatibility publication or query-specific stop policy. */
        if (out_selected) {
            result->selected_intervals_yielded++;
            if (need_boundary_event && enrich_local_group)
                result->boundary_event_enrichments++;
            prev_cell_idx = cell_idx;
            if (next_enter_surface_id < 0)
                next_enter_surface_id = hit_surface_id;
            prev_surface_id = next_enter_surface_id;
            enter_event = bevent;
            t_current = t_next;
            state->t_current = t_current;
            state->prev_cell_idx = prev_cell_idx;
            state->prev_surface_id = prev_surface_id;
            state->pending_enter_hit_index = pending_enter_hit_index;
            state->enter_event = enter_event;
            state->pending_lattice_entry_sample = pending_lattice_entry_sample;
            *out_selected = selected;
            out_selected->path = &state->current_path;
            return 0;
        }

        if (ray_selected_interval_publish(
                sys, ray, &selected, use_hier_lookup, emit_hits, result,
                &prev_cell_idx, &pending_enter_hit_index) != 0)
            return -1;

        if (next_enter_surface_id < 0)
            next_enter_surface_id = hit_surface_id;
        prev_surface_id = next_enter_surface_id;
        enter_event = bevent;

        /* Move past the intersection */
        t_current = t_next;

        /* Persist only the state needed to resume at the next verified
         * interval. Per-step candidate and retry variables remain local. */
        state->t_current = t_current;
        state->prev_cell_idx = prev_cell_idx;
        state->prev_surface_id = prev_surface_id;
        state->pending_enter_hit_index = pending_enter_hit_index;
        state->enter_event = enter_event;
        state->pending_lattice_entry_sample = pending_lattice_entry_sample;

        /* If we hit nothing, we're done */
        if (t_next >= effective_t_max) {
            break;
        }
    }

    if (t_current < effective_t_max && state->iterations_remaining <= 0) {
        alea_set_error_detail(
            ALEA_ERR_OVERFLOW,
            "ray traversal exceeded %d intervals before t_max",
            MAX_RAY_WALK_ITERATIONS);
        return -1;
    }
    if (t_current < effective_t_max && step_budget <= 0) {
        return 1;  /* yielded after the requested number of intervals */
    }
    return 0;
}

/* Legacy flat-cell compatibility adapter. Hierarchical consumers stream
 * selected intervals through their dedicated adapters instead. */
static int raycast_flat_cell_aware_impl(alea_system_t* sys,
                                        const alea_ray_t* ray,
                                        double effective_t_max,
                                        alea_raycast_result_t* result) {
    alea_ray_walk_t state;
    alea_ray_walk_init(&state);
    return raycast_cell_aware_resume(
        sys, ray, effective_t_max, false, false, false,
        result, &state, INT_MAX, NULL);
}

/* Advance the coherent hierarchical walk by one verified interval. The
 * returned interval is scratch-backed through `walk` and must be consumed
 * before the next call. Returns 1 when another interval may follow, 0 for
 * the terminal interval, and 2 when the walk was already exhausted. */
static int alea_ray_walk_next_selected_mode(
    alea_system_t* sys, const alea_ray_t* ray, double t_max,
    alea_raycast_result_t* scratch, alea_ray_walk_t* walk,
    alea_ray_selected_interval_t* out_interval,
    bool enrich_local_group) {
    if (!isfinite(t_max) || !(t_max > 0.0)) {
        alea_set_error_detail(ALEA_ERR_INVALID_ARG,
                              "ray traversal requires a finite positive t_max");
        return -1;
    }
    if (!out_interval || walk->t_current >= t_max)
        return 2;
    const int rc = raycast_cell_aware_resume(
        sys, ray, t_max, true, false, enrich_local_group,
        scratch, walk, 1, out_interval);
    if (rc != 0) return rc;
    return out_interval->t_exit >= t_max ? 0 : 1;
}

static int alea_ray_walk_next_selected(
    alea_system_t* sys, const alea_ray_t* ray, double t_max,
    alea_raycast_result_t* scratch, alea_ray_walk_t* walk,
    alea_ray_selected_interval_t* out_interval) {
    return alea_ray_walk_next_selected_mode(
        sys, ray, t_max, scratch, walk, out_interval, true);
}

/* Persistent particle navigator. The walker keeps the exact selected path;
 * the pending interval is retained across collision and distance-limit stops.
 * Neither step nor clone publishes a segment vector. */
struct alea_ray_navigator {
    interval_index_t* regional_index;
    regional_worker_t* regional_worker;
    size_t regional_spans;
    alea_nav_recovery_stats_t recovery_stats;
    alea_transport_recovery_seed_t recovery_seed;
    double recovery_remaining;
    uint64_t recovered_outgoing; /* Nonzero while using the incident material. */
    alea_uint32_vec_t action_surfaces;
    bool window_end;
    alea_system_t* sys;
    alea_nav_validation_mode_t validation_mode;
    size_t max_interval_breakpoints;
    uint32_t event_fields;
    double boundary_distance_tolerance;
    alea_ray_t ray;
    alea_ray_walk_t walk;
    alea_ray_selected_interval_t pending;
    alea_raycast_result_t scratch;
    alea_raycast_result_t validation_scratch;
    double current_t;
    alea_nav_location_t location;
    bool initialized;
    bool have_pending;
    alea_ray_substep_fill_t substep_fill;
    bool awaiting_action;
    bool leaked;
    alea_boundary_type_t action_type;
    int action_surface_id;
    double action_normal[3];
};

static bool transport_span_verified(alea_ray_navigator_t* nav,
        const alea_nav_location_t* incident, int enter_surface, uint64_t* outgoing);

int alea_ray_navigator_recovery_stats(const alea_ray_navigator_t* nav,
                                      alea_nav_recovery_stats_t* out) {
    if (!nav || !out) return -1;
    *out = nav->recovery_stats;
    return 0;
}

static alea_nav_location_t navigator_interval_location(
    const alea_ray_selected_interval_t* interval) {
    alea_nav_location_t location = {
        .kind = interval->cell_id < 0 ? ALEA_NAV_GAP
            : interval->material_id == 0 ? ALEA_NAV_VOID : ALEA_NAV_MATERIAL,
        .cell_id = interval->cell_id,
        .cell_index = interval->cell_id < 0 ? -1 : interval->cell_index,
        .material_id = interval->material_id,
        .density = interval->density,
        .occurrence_key = interval->owner_provenance_complete
            ? interval->owner_occurrence_key : 0
    };
    if (interval->resolution_flags & ALEA_RESOLVE_UNDEFINED_FILL)
        location.kind = ALEA_NAV_UNDEFINED_FILL;
    else if (interval->cell_id >= 0 &&
             !interval->owner_provenance_complete)
        location.kind = ALEA_NAV_UNRESOLVED;
    return location;
}

#define ALEA_NAV_OWNER_CAPACITY 64
static alea_nav_location_kind_t navigator_verify_point(
    const alea_ray_navigator_t* navigator, double sample) {
    const alea_ray_selected_interval_t* interval = &navigator->pending;
    double x, y, z;
    alea_ray_point_at(&navigator->ray, sample, &x, &y, &z);
    if (!isfinite(x) || !isfinite(y) || !isfinite(z))
        return ALEA_NAV_UNRESOLVED;
    alea_cell_hit_t hits[ALEA_NAV_OWNER_CAPACITY];
    uint64_t keys[ALEA_NAV_OWNER_CAPACITY];
    uint64_t parents[ALEA_NAV_OWNER_CAPACITY];
    uint8_t owners[ALEA_NAV_OWNER_CAPACITY];
    size_t child_counts[ALEA_NAV_OWNER_CAPACITY];
    int count = alea_find_all_cells_at_point_coverage_chain_recursive(
        navigator->sys, x, y, z, hits, keys, parents,
        ALEA_NAV_OWNER_CAPACITY);
    if (count < 0 || count >= ALEA_NAV_OWNER_CAPACITY)
        return ALEA_NAV_UNRESOLVED;
    alea_point_coverage_classification_t classification;
    if (alea_classify_point_coverage_chain_with_scratch(
        hits, keys, parents, (size_t)count, -1, owners, child_counts,
        &classification) != 0) return ALEA_NAV_UNRESOLVED;
    switch (classification.kind) {
        case ALEA_POINT_COVERAGE_GAP:
            return interval->cell_id < 0
                ? ALEA_NAV_GAP : ALEA_NAV_UNRESOLVED;
        case ALEA_POINT_COVERAGE_OVERLAP:
            return ALEA_NAV_OVERLAP;
        case ALEA_POINT_COVERAGE_UNDEFINED_FILL:
            return ALEA_NAV_UNDEFINED_FILL;
        case ALEA_POINT_COVERAGE_UNRESOLVED:
            return ALEA_NAV_UNRESOLVED;
        case ALEA_POINT_COVERAGE_UNIQUE:
            for (int i = 0; i < count; i++) {
                if (!owners[i]) continue;
                if (!interval->owner_provenance_complete ||
                    interval->owner_occurrence_key != keys[i] ||
                    interval->cell_id != hits[i].cell_id)
                    return ALEA_NAV_UNRESOLVED;
                return interval->material_id == 0
                    ? ALEA_NAV_VOID : ALEA_NAV_MATERIAL;
            }
            return ALEA_NAV_UNRESOLVED;
    }
    return ALEA_NAV_UNRESOLVED;
}

static alea_nav_location_kind_t navigator_verify_interval(
    alea_ray_navigator_t* navigator) {
    const alea_ray_selected_interval_t* interval = &navigator->pending;
    if (navigator->substep_fill.enter_surface_id || !isfinite(interval->t_enter))
        return ALEA_NAV_UNRESOLVED;
    alea_nav_location_kind_t kind = ALEA_NAV_UNRESOLVED;
    const double span = interval->t_exit - interval->t_enter;
    double sample = interval->t_enter + 0.381966011250105 * fmin(span, 1.0);
    if (nextafter(interval->t_enter, INFINITY) < interval->t_exit) {
        kind = navigator_verify_point(navigator, sample);
    } else if (interval->resolution_flags & ALEA_RESOLVE_UNDEFINED_FILL) {
        alea_ray_substep_fill_t gap;
        if (alea_substep_fill_at(navigator->sys, &navigator->ray,
                interval->t_enter, &gap) == 1 && gap.parent_cell_id == interval->cell_id)
            navigator->substep_fill = gap;
    }
    if (navigator->regional_worker && (navigator->validation_mode == ALEA_NAV_VALIDATE_STRICT ||
        navigator->validation_mode == ALEA_NAV_VALIDATE_TRANSPORT) &&
        (kind == ALEA_NAV_OVERLAP || kind == ALEA_NAV_UNRESOLVED) &&
        interval->owner_provenance_complete && !interval->resolution_flags &&
        alea_ray_verify_short_owner(navigator->sys, &navigator->ray,
            interval->t_enter, interval->t_exit, sample, interval->path))
        return interval->material_id ? ALEA_NAV_MATERIAL : ALEA_NAV_VOID;
    return kind;
}

#define ALEA_NAV_DEFAULT_FLIGHT_BREAKPOINTS 8192
typedef struct {
    alea_ray_navigator_t* navigator;
    int failed;
} alea_nav_interval_check_t;

static int navigator_check_coverage_interval(
    void* context, const alea_ray_coverage_interval_t* interval) {
    alea_nav_interval_check_t* check = context;
    const alea_ray_navigator_t* navigator = check->navigator;
    alea_nav_location_kind_t observed = ALEA_NAV_UNRESOLVED;
    const alea_ray_coverage_owner_t* owner = NULL;
    if (interval->kind == ALEA_RAY_COVERAGE_GAP) {
        observed = ALEA_NAV_GAP;
    } else if (interval->kind == ALEA_RAY_COVERAGE_OVERLAP) {
        observed = ALEA_NAV_OVERLAP;
    } else if (interval->kind == ALEA_RAY_COVERAGE_UNDEFINED_FILL) {
        observed = ALEA_NAV_UNDEFINED_FILL;
    } else if (interval->kind == ALEA_RAY_COVERAGE_UNIQUE) {
        for (size_t i = 0; i < interval->owner_count; i++) {
            if (!owner || interval->owners[i].depth > owner->depth)
                owner = &interval->owners[i];
        }
        if (owner) observed = owner->material_id == 0
            ? ALEA_NAV_VOID : ALEA_NAV_MATERIAL;
    }
    if (observed != navigator->location.kind ||
        (owner && (owner->cell_id != navigator->pending.cell_id ||
                   owner->material_id != navigator->pending.material_id ||
                   !navigator->pending.owner_provenance_complete ||
                   owner->occurrence_key !=
                       navigator->pending.owner_occurrence_key))) {
        alea_set_error_detail(ALEA_ERR_INVALID_STATE,
            "particle interval ownership differs at t=%.17g: selected "
            "cell=%d occurrence=%llu kind=%d, observed cell=%d kind=%d",
            interval->t_enter, navigator->pending.cell_id,
            (unsigned long long)navigator->pending.owner_occurrence_key,
            (int)navigator->location.kind, owner ? owner->cell_id : -1,
            (int)observed);
        check->failed = 1;
        return 1;
    }
    return 0;
}

static int navigator_validate_flight(alea_ray_navigator_t* navigator,
                                     double t_end) {
    alea_raycast_result_t* scratch = &navigator->validation_scratch;
    if (!(t_end > navigator->current_t) || !isfinite(t_end)) {
        alea_set_error_detail(ALEA_ERR_INVALID_STATE,
                              "particle interval endpoint cannot advance");
        return -1;
    }
    alea_nav_interval_check_t check = { .navigator = navigator };
    const int count = alea_ray_coverage_sweep_strict_reuse_nocache(
        navigator->sys, &navigator->ray, navigator->current_t, t_end,
        navigator->max_interval_breakpoints, NULL, scratch,
        navigator_check_coverage_interval, &check);
    return count < 0 || check.failed ? -1 : 0;
}

/* A Boolean union can keep the same owner across a reflective, white,
 * periodic or vacuum card. Preserve these actions inside a selected span.
 * In a gap the native navigator considers world-frame surfaces; inside
 * a cell it tests cards along the concrete incident placement path. */
static int navigator_check_action_surfaces(const alea_ray_navigator_t* nav,
        double start, regional_nav_span_t* span) {
    if (!nav->action_surfaces.count) return 0;
    alea_matrix_t identity;
    alea_matrix_identity(&identity);
    const alea_hier_ray_path_t* path = span->path;
    int frames = path->count ? path->count : 1;
    for (int i = 0; i < frames; ++i) {
        const alea_cell_entry_t* cell = path->count
            ? &nav->sys->cells.data[path->entries[i].cell_index] : NULL;
        const alea_matrix_t* frame = cell ? &path->entries[i].transform : &identity;
        const uint32_t* indices = cell ? cell->surface_indices : nav->action_surfaces.data;
        size_t count = cell ? cell->surface_index_count : nav->action_surfaces.count;
        alea_ray_t local;
        alea_ray_transform_inverse(frame, &nav->ray, &local);
        for (size_t j = 0; j < count; ++j) {
            const alea_surface_entry_t* surface = &nav->sys->surfaces.data[indices[j]];
            if (surface->boundary_type == ALEA_BOUNDARY_TRANSMISSIVE) continue;
            uint32_t id = surface->primitive_id;
            if (id >= nav->sys->primitives.count) return -1;
            alea_primitive_type_t type = nav->sys->primitives.data[id].type;
            alea_primitive_data_t data;
            if (!raycast_primitive_copy_payload(nav->sys, id, type, &data)) return -1;
            double roots[4];
            int n = ray_intersect_primitive(&local, type, &data, roots);
            for (int k = 0; k < n; ++k) {
                if (roots[k] > start && roots[k] < span->end) {
                    span->end = roots[k];
                    span->surface = surface->mc_surface_id;
                    span->primitive = id;
                    span->exit_frame = *frame;
                }
            }
        }
    }
    return 0;
}

/* A search-window end certifies a span, but carries no physical action.
 * When a request reaches it, repeat the search with a longer lookahead. */
static int navigator_regional_span(alea_ray_navigator_t* nav, double start, double limit) {
    regional_nav_span_t span = {.path=&nav->walk.current_path};
    if (alea_regional_first_span(nav->regional_worker, &nav->ray, start, limit,
                                 &nav->scratch, &span) || navigator_check_action_surfaces(nav, start, &span)) {
        alea_set_error_detail(ALEA_ERR_INVALID_STATE, "Regional navigator could not establish the next span");
        return -1;
    }
    ++nav->regional_spans;
    const alea_cell_entry_t* cell = span.cell >= 0 ? &nav->sys->cells.data[span.cell] : NULL;
    nav->window_end = span.end == limit;
    nav->pending = (alea_ray_selected_interval_t){
        .t_enter=start, .t_exit=span.end, .cell_index=span.cell,
        .cell_id=cell ? cell->mc_cell_id : -1, .material_id=cell ? cell->material_id : 0,
        .density=cell ? cell->density : 0, .resolution_flags=span.resolution_flags,
        .path=span.path, .exit_surface_id=span.surface,
        .exit_event={.t=span.end, .surface_id=span.surface, .primitive_id=span.primitive,
            .has_physical_surface=span.surface > 0, .is_synthetic_lattice_boundary=span.surface == 0,
            .transform=span.exit_frame}
    };
    boundary_interval_capture_owner(span.path, &nav->pending);
    if (!nav->window_end && span.surface < 0) {
        alea_set_error_detail(ALEA_ERR_INVALID_STATE, "Regional navigator cannot attribute a crossing");
        return -1;
    }
    return 0;
}

int alea_ray_navigator_use_regional(alea_ray_navigator_t* nav, interval_index_t* idx) {
    if (!nav || nav->initialized || nav->regional_worker ||
        !alea_regional_nav_ready(idx, nav->sys)) return -1;
    nav->regional_worker = alea_regional_worker_create(idx);
    if (!nav->regional_worker) return -1;
    for (size_t i = 0; i < nav->sys->surfaces.count; ++i) {
        if (nav->sys->surfaces.data[i].boundary_type != ALEA_BOUNDARY_TRANSMISSIVE &&
            alea_vec_push(&nav->action_surfaces, (uint32_t)i, uint32_t)) {
            alea_regional_worker_destroy(nav->regional_worker);
            nav->regional_worker = NULL;
            alea_vec_free(&nav->action_surfaces);
            return -1;
        }
    }
    alea_regional_index_retain(idx);
    nav->regional_index = idx;
    return 0;
}

void alea_ray_navigator_regional_counts(const alea_ray_navigator_t* nav,
                                       size_t* spans, size_t* visits) {
    *spans = nav->regional_spans;
    *visits = nav->regional_worker ? alea_regional_worker_candidate_tree_visits(nav->regional_worker) : 0;
}

static bool navigator_index_current(const alea_ray_navigator_t* nav) {
    return !nav->regional_index || alea_regional_nav_ready(nav->regional_index, nav->sys);
}

static int navigator_load_interval(alea_ray_navigator_t* navigator) {
    navigator->walk.recovery_proof_budget = navigator->max_interval_breakpoints;
    int rc;
    if (navigator->regional_worker) {
        double start = navigator->current_t;
        double limit = fmax(start + 64.0, nextafter(start, INFINITY));
        rc = navigator_regional_span(navigator, start, limit);
    } else {
        rc = alea_ray_walk_next_selected_mode(
            navigator->sys, &navigator->ray, DBL_MAX, &navigator->scratch,
            &navigator->walk, &navigator->pending,
            navigator->validation_mode != ALEA_NAV_VALIDATE_FAST);
    }
    if (rc < 0) return -1;
    if (rc == 2) {
        alea_set_error_detail(ALEA_ERR_INVALID_STATE,
                              "particle navigator has no next interval");
        return -1;
    }
    navigator->have_pending = true;
    navigator->location = navigator_interval_location(&navigator->pending);
    if (navigator->validation_mode != ALEA_NAV_VALIDATE_FAST)
        navigator->location.kind = navigator_verify_interval(navigator);
    return 0;
}

alea_ray_navigator_t* alea_ray_navigator_create(alea_system_t* sys) {
    if (!sys) return NULL;
    alea_ray_navigator_t* navigator = calloc(1, sizeof(*navigator));
    if (!navigator) return NULL;
    navigator->sys = sys;
    navigator->validation_mode = ALEA_NAV_VALIDATE_STRICT;
    navigator->max_interval_breakpoints = ALEA_NAV_DEFAULT_FLIGHT_BREAKPOINTS;
    navigator->event_fields = ALEA_NAV_EVENT_NORMAL;
    navigator->boundary_distance_tolerance =
        ALEA_NAV_DEFAULT_BOUNDARY_DISTANCE_TOLERANCE;
    alea_raycast_result_init(&navigator->scratch);
    alea_raycast_result_init(&navigator->validation_scratch);
    return navigator;
}

void alea_ray_navigator_destroy(alea_ray_navigator_t* navigator) {
    if (!navigator) return;
    alea_raycast_result_free(&navigator->scratch);
    alea_raycast_result_free(&navigator->validation_scratch);
    alea_vec_free(&navigator->action_surfaces);
    alea_regional_worker_destroy(navigator->regional_worker);
    alea_prepared_queries_destroy(navigator->regional_index);
    free(navigator);
}

int alea_ray_navigator_set_validation_mode(
    alea_ray_navigator_t* navigator, alea_nav_validation_mode_t mode) {
    if (!navigator || !navigator_index_current(navigator) || (mode != ALEA_NAV_VALIDATE_FAST &&
                       mode != ALEA_NAV_VALIDATE_STRICT &&
                       mode != ALEA_NAV_VALIDATE_INTERVAL && mode != ALEA_NAV_VALIDATE_TRANSPORT)) return -1;
    if (mode == ALEA_NAV_VALIDATE_TRANSPORT && navigator->recovered_outgoing) return 0;
    navigator->recovered_outgoing = 0;
    navigator->validation_mode = mode;
    if (navigator->initialized && navigator->have_pending) {
        navigator->location = navigator_interval_location(&navigator->pending);
        if (mode != ALEA_NAV_VALIDATE_FAST)
            navigator->location.kind = navigator_verify_interval(navigator);
    }
    return 0;
}

int alea_ray_navigator_substep_fill(const alea_ray_navigator_t* navigator,
                                    alea_ray_substep_fill_t* out) {
    if (!navigator || !out) return -1;
    if (!navigator->initialized || !navigator->substep_fill.enter_surface_id ||
        navigator->location.kind != ALEA_NAV_UNRESOLVED) return 0;
    *out = navigator->substep_fill;
    return 1;
}

int alea_ray_navigator_set_interval_budget(
    alea_ray_navigator_t* navigator, size_t max_breakpoints) {
    if (!navigator || max_breakpoints == 0) return -1;
    navigator->max_interval_breakpoints = max_breakpoints;
    return 0;
}

int alea_ray_navigator_set_boundary_distance_tolerance(
    alea_ray_navigator_t* navigator, double relative_tolerance) {
    if (!navigator || !isfinite(relative_tolerance) ||
        relative_tolerance < 0.0 || navigator->initialized) return -1;
    navigator->boundary_distance_tolerance = relative_tolerance;
    navigator->walk.boundary_distance_tolerance = relative_tolerance;
    return 0;
}

int alea_ray_navigator_set_event_fields(
    alea_ray_navigator_t* navigator, uint32_t fields) {
    if (!navigator || (fields & ~ALEA_NAV_EVENT_NORMAL)) return -1;
    navigator->event_fields = fields;
    return 0;
}

int alea_ray_navigator_restart(alea_ray_navigator_t* navigator,
                               const double position[3],
                               const double direction[3],
                               alea_nav_location_t* location) {
    if (!navigator || !navigator_index_current(navigator) || !position || !direction || !location ||
        !isfinite(position[0]) || !isfinite(position[1]) ||
        !isfinite(position[2])) return -1;
    alea_ray_t ray;
    if (alea_ray_init(&ray, position[0], position[1], position[2],
                      direction[0], direction[1], direction[2]) != 0 ||
        alea_raycast_ensure_hier_caches(navigator->sys) != 0) return -1;
    navigator->recovered_outgoing = 0;
    navigator->ray = ray;
    alea_ray_walk_init(&navigator->walk);
    navigator->walk.fail_unverified_ownership = true;
    navigator->walk.boundary_distance_tolerance =
        navigator->boundary_distance_tolerance;
    alea_raycast_result_clear(&navigator->scratch);
    navigator->scratch.ray = ray;
    navigator->substep_fill.enter_surface_id = 0;
    navigator->current_t = 0.0;
    navigator->have_pending = false;
    navigator->awaiting_action = false;
    navigator->leaked = false;
    navigator->action_type = ALEA_BOUNDARY_TRANSMISSIVE;
    navigator->action_surface_id = -1;
    navigator->initialized = false;
    if (navigator_load_interval(navigator) != 0) return -1;
    navigator->initialized = true;
    *location = navigator->location;
    return 0;
}

alea_ray_navigator_t* alea_ray_navigator_clone(
    const alea_ray_navigator_t* navigator) {
    if (!navigator) return NULL;
    alea_ray_navigator_t* copy = malloc(sizeof(*copy));
    if (!copy) return NULL;
    *copy = *navigator;
    copy->recovery_stats = (alea_nav_recovery_stats_t){0};
    alea_raycast_result_init(&copy->scratch);
    alea_raycast_result_init(&copy->validation_scratch);
    copy->scratch.ray = copy->ray;
    if (copy->regional_index) {
        copy->regional_worker = NULL;
        memset(&copy->action_surfaces, 0, sizeof(copy->action_surfaces));
        alea_regional_index_retain(copy->regional_index);
        if (!navigator_index_current(copy) ||
            !(copy->regional_worker = alea_regional_worker_create(copy->regional_index)) ||
            alea_vec_reserve(&copy->action_surfaces, navigator->action_surfaces.count, uint32_t)) {
            alea_ray_navigator_destroy(copy);
            return NULL;
        }
        if (navigator->action_surfaces.count) {
            copy->action_surfaces.count = navigator->action_surfaces.count;
            memcpy(copy->action_surfaces.data, navigator->action_surfaces.data,
                   copy->action_surfaces.count * sizeof(uint32_t));
        }
        copy->regional_spans = 0;
    }
    if (navigator->pending.path) copy->pending.path = &copy->walk.current_path;
    return copy;
}

int alea_ray_navigator_set_direction(alea_ray_navigator_t* navigator,
                                     const double direction[3]) {
    if (!navigator || !navigator_index_current(navigator) || !direction || !navigator->initialized ||
        navigator->leaked || (navigator->substep_fill.enter_surface_id &&
            navigator->validation_mode != ALEA_NAV_VALIDATE_FAST)) return -1;
    double position[3];
    alea_ray_point_at(&navigator->ray, navigator->current_t,
                      &position[0], &position[1], &position[2]);
    alea_ray_t ray;
    if (alea_ray_init(&ray, position[0], position[1], position[2],
                      direction[0], direction[1], direction[2]) != 0)
        return -1;
    /* A collision keeps the owner. At a reflective/white boundary the old
     * path still identifies the incident side. The walker checks containment
     * of this hint in the new direction before it trusts it. */
    const uint64_t previous_outgoing = navigator->recovered_outgoing;
    navigator->recovered_outgoing = 0;
    const int hinted_cell = navigator->have_pending
        ? navigator->pending.cell_index : navigator->walk.prev_cell_idx;
    alea_hier_ray_path_t path;
    ray_path_copy_live(&path, &navigator->walk.current_path);
    navigator->ray = ray;
    alea_ray_walk_init(&navigator->walk);
    navigator->walk.fail_unverified_ownership = true;
    navigator->walk.boundary_distance_tolerance =
        navigator->boundary_distance_tolerance;
    ray_path_copy_live(&navigator->walk.current_path, &path);
    navigator->walk.prev_cell_idx = hinted_cell;
    navigator->scratch.ray = ray;
    navigator->substep_fill.enter_surface_id = 0;
    navigator->current_t = 0.0;
    navigator->have_pending = false;
    navigator->awaiting_action = false;
    navigator->action_type = ALEA_BOUNDARY_TRANSMISSIVE;
    navigator->action_surface_id = -1;
    navigator->initialized = false;
    if (navigator_load_interval(navigator) != 0) return -1;
    if (previous_outgoing && alea_nav_resume_recovery(navigator, &navigator->recovery_seed,
                                                        &navigator->location)) return -1;
    navigator->initialized = true;
    return 0;
}

int alea_ray_navigator_reflect_specular(alea_ray_navigator_t* navigator) {
    if (!navigator || !navigator->initialized ||
        !navigator->awaiting_action ||
        navigator->action_type != ALEA_BOUNDARY_REFLECTIVE) return -1;
    const double* normal = navigator->action_normal;
    double norm2 = normal[0] * normal[0] + normal[1] * normal[1] +
                   normal[2] * normal[2];
    if (!isfinite(norm2) || norm2 < 0.5 || norm2 > 1.5) return -1;
    const double dot = navigator->ray.dx * normal[0] +
                       navigator->ray.dy * normal[1] +
                       navigator->ray.dz * normal[2];
    const double direction[3] = {
        navigator->ray.dx - 2.0 * dot * normal[0] / norm2,
        navigator->ray.dy - 2.0 * dot * normal[1] / norm2,
        navigator->ray.dz - 2.0 * dot * normal[2] / norm2
    };
    return alea_ray_navigator_set_direction(navigator, direction);
}

static const alea_surface_entry_t* navigator_periodic_partner(
    const alea_system_t* sys, int source_id) {
    const int source_index = alea_surface_find(sys, source_id);
    if (source_index < 0) return NULL;
    const alea_surface_entry_t* source = &sys->surfaces.data[source_index];
    int partner_id = source->periodic_surface_id;
    if (partner_id == 0) {
        for (size_t i = 0; i < alea_vec_count(&sys->surfaces); i++) {
            const alea_surface_entry_t* candidate = &sys->surfaces.data[i];
            if (candidate->periodic_surface_id != source_id) continue;
            if (partner_id != 0) return NULL;
            partner_id = candidate->mc_surface_id;
        }
    }
    if (partner_id <= 0 || partner_id == source_id) return NULL;
    const int partner_index = alea_surface_find(sys, partner_id);
    if (partner_index < 0) return NULL;
    const alea_surface_entry_t* partner = &sys->surfaces.data[partner_index];
    if (partner->boundary_type != ALEA_BOUNDARY_PERIODIC ||
        (partner->periodic_surface_id != 0 &&
         partner->periodic_surface_id != source_id)) return NULL;
    if (partner->periodic_surface_id == 0) {
        size_t references = 0;
        for (size_t i = 0; i < alea_vec_count(&sys->surfaces); i++)
            if (sys->surfaces.data[i].periodic_surface_id == partner_id)
                references++;
        if (references != 1) return NULL;
    }
    return partner;
}

int alea_ray_navigator_apply_periodic(alea_ray_navigator_t* navigator,
                                      double out_position[3],
                                      alea_nav_location_t* out_location) {
    if (!navigator || !navigator_index_current(navigator) || !out_position || !out_location ||
        !navigator->initialized || !navigator->awaiting_action ||
        navigator->action_type != ALEA_BOUNDARY_PERIODIC ||
        navigator->action_surface_id <= 0) return -1;
    if (navigator->pending.path && navigator->pending.path->count > 1) {
        alea_set_error_detail(ALEA_ERR_UNSUPPORTED,
            "periodic mapping inside a nested occurrence is unsupported");
        return -1;
    }
    const int source_index = alea_surface_find(navigator->sys,
                                               navigator->action_surface_id);
    const alea_surface_entry_t* partner = navigator_periodic_partner(
        navigator->sys, navigator->action_surface_id);
    if (source_index < 0 || !partner) {
        alea_set_error_detail(ALEA_ERR_INVALID_STATE,
            "periodic boundary has no unique paired surface");
        return -1;
    }
    const alea_surface_entry_t* source =
        &navigator->sys->surfaces.data[source_index];
    if (source->primitive_id >= alea_vec_count(&navigator->sys->primitives) ||
        partner->primitive_id >= alea_vec_count(&navigator->sys->primitives) ||
        navigator->sys->primitives.data[source->primitive_id].type !=
            ALEA_PRIMITIVE_PLANE ||
        navigator->sys->primitives.data[partner->primitive_id].type !=
            ALEA_PRIMITIVE_PLANE) {
        alea_set_error_detail(ALEA_ERR_UNSUPPORTED,
            "periodic transport requires a pair of planes");
        return -1;
    }
    alea_primitive_data_t first, second;
    if (!raycast_primitive_copy_payload(navigator->sys, source->primitive_id,
            ALEA_PRIMITIVE_PLANE, &first) ||
        !raycast_primitive_copy_payload(navigator->sys, partner->primitive_id,
            ALEA_PRIMITIVE_PLANE, &second)) {
        alea_set_error_detail(ALEA_ERR_INVALID_STATE,
            "periodic plane data is unavailable");
        return -1;
    }
    const double a[3] = {first.plane.a, first.plane.b, first.plane.c};
    const double b[3] = {second.plane.a, second.plane.b, second.plane.c};
    const double an = sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2]);
    const double bn = sqrt(b[0]*b[0] + b[1]*b[1] + b[2]*b[2]);
    const double dot = a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
    if (!(an > 0.0) || !(bn > 0.0) ||
        !isfinite(an) || !isfinite(bn) ||
        !isfinite(dot) || !isfinite(an * bn) ||
        fabs(fabs(dot / (an * bn)) - 1.0) > 1e-10) {
        alea_set_error_detail(ALEA_ERR_UNSUPPORTED,
            "periodic transport requires parallel planes");
        return -1;
    }
    double position[3];
    alea_ray_point_at(&navigator->ray, navigator->current_t,
                      &position[0], &position[1], &position[2]);
    const double distance = -(
        b[0]*position[0] + b[1]*position[1] + b[2]*position[2] +
        second.plane.d) * an / dot;
    if (!isfinite(distance) || distance == 0.0) {
        alea_set_error_detail(ALEA_ERR_INVALID_STATE,
            "periodic planes have no distinct finite translation");
        return -1;
    }
    double mapped[3];
    for (int i = 0; i < 3; i++)
        mapped[i] = position[i] + distance * a[i] / an;
    if (!isfinite(mapped[0]) || !isfinite(mapped[1]) || !isfinite(mapped[2])) {
        alea_set_error_detail(ALEA_ERR_OVERFLOW,
            "periodic mapping produced a nonfinite position");
        return -1;
    }
    if (mapped[0] == position[0] && mapped[1] == position[1] &&
        mapped[2] == position[2]) {
        alea_set_error_detail(ALEA_ERR_INVALID_STATE,
            "periodic translation cannot advance the position");
        return -1;
    }
    alea_nav_location_t location;
    const double direction[3] = {
        navigator->ray.dx, navigator->ray.dy, navigator->ray.dz
    };
    if (alea_ray_navigator_restart(navigator, mapped, direction,
                                   &location) != 0) return -1;
    if (location.kind != ALEA_NAV_MATERIAL &&
        location.kind != ALEA_NAV_VOID) {
        alea_set_error_detail(ALEA_ERR_INVALID_STATE,
            "periodic destination has no unique material or void owner");
        return -1;
    }
    memcpy(out_position, mapped, sizeof(mapped));
    *out_location = location;
    return 0;
}

static bool boundary_group_contains(const alea_raycast_boundary_event_t* group, int id) {
    for (size_t i = 0; i < group->local_surface_count; ++i)
        if (group->local_surface_ids[i] == id) return true;
    return false;
}

/* An entering cell can supply the regional root. On an incident lookup miss,
 * collect cards from both concrete placements before resolving the action.
 * Keep all speculative state separate from the particle's pending interval. */
static bool navigator_two_sided_group(alea_ray_navigator_t* nav,
        alea_raycast_boundary_event_t* group) {
    const double t = group->t;
    alea_hier_ray_path_t outgoing;
    if (nav->regional_worker) {
        regional_nav_span_t span = {.path = &outgoing};
        double limit = fmax(t + 64.0, nextafter(t, INFINITY));
        ++nav->regional_spans;
        if (alea_regional_first_span(nav->regional_worker, &nav->ray, t, limit,
                                     &nav->scratch, &span) || !outgoing.count)
            return false;
    } else {
        alea_ray_walk_t walk = nav->walk;
        alea_raycast_result_t scratch;
        alea_raycast_result_init(&scratch);
        scratch.ray = nav->ray;
        alea_ray_selected_interval_t next;
        int rc = alea_ray_walk_next_selected_mode(nav->sys, &nav->ray, DBL_MAX,
                                                   &scratch, &walk, &next, true);
        alea_raycast_result_free(&scratch);
        if (rc < 0 || rc == 2 || next.t_enter != t || !next.path || !next.path->count)
            return false;
        ray_path_copy_live(&outgoing, next.path);
    }
    const int index = alea_surface_find(nav->sys, group->surface_id);
    if (index < 0) return false;
    const alea_surface_entry_t* selected = &nav->sys->surfaces.data[index];
    if (selected->primitive_id != group->primitive_id) return false;

    bool found_outgoing = false;
    const alea_hier_ray_path_t* paths[] = {nav->pending.path, &outgoing};
    group->local_surface_complete = 1;
    group->local_surface_count = 0;
    for (int side = 0; side < 2; ++side) {
        for (int i = 0; i < paths[side]->count; ++i) {
            const alea_hier_ray_path_entry_t* entry = &paths[side]->entries[i];
            if ((size_t)entry->cell_index >= nav->sys->cells.count) return false;
            alea_ray_t local;
            alea_ray_transform_inverse(&entry->transform, &nav->ray, &local);
            alea_raycast_boundary_event_t cards = {.local_surface_complete = 1};
            if (!boundary_event_collect_cell_surfaces_at(nav->sys, &local,
                    &nav->sys->cells.data[entry->cell_index], t, &cards))
                return false;
            if (side == 1 && boundary_group_contains(&cards, selected->mc_surface_id)) {
                bool same_frame = memcmp(entry->transform.m, group->transform.m,
                                         sizeof(group->transform.m)) == 0;
                if (selected->boundary_type != ALEA_BOUNDARY_TRANSMISSIVE &&
                    (!same_frame || (i > 0 && selected->boundary_type == ALEA_BOUNDARY_PERIODIC)))
                    return false;
                found_outgoing |= same_frame;
            }
            for (size_t j = 0; j < cards.local_surface_count; ++j)
                boundary_event_add_local_surface(group, cards.local_surface_ids[j]);
        }
    }
    return found_outgoing && group->local_surface_complete;
}

/* Resolve conditions from the incident path, including physical surfaces
 * coincident with a synthetic lattice crossing. Multiple physical actions
 * leave the action normal or destination ambiguous, so stop before moving. */
static int navigator_boundary_condition(
    alea_ray_navigator_t* navigator,
    const alea_raycast_boundary_event_t* crossing,
    bool verify_outgoing, alea_boundary_type_t* out_type) {
    alea_raycast_boundary_event_t group = *crossing;
    if (!navigator->pending.path || navigator->pending.path->count == 0) {
        boundary_event_collect_global_group(navigator->sys, &navigator->ray,
                                            &group);
    } else if (navigator->regional_worker || navigator->validation_mode == ALEA_NAV_VALIDATE_FAST ||
        crossing->is_synthetic_lattice_boundary ||
        !crossing->has_physical_surface) {
        group.has_physical_surface = true;
        group.is_synthetic_lattice_boundary = false;
        boundary_event_collect_local_group(navigator->sys, &navigator->ray,
            navigator->pending.path, &group);
    }
    if ((navigator->regional_worker || verify_outgoing) &&
        navigator->pending.path && navigator->pending.path->count &&
        crossing->surface_id > 0 && !boundary_group_contains(&group, crossing->surface_id)) {
        if (!navigator_two_sided_group(navigator, &group)) {
            alea_set_error_detail(ALEA_ERR_INVALID_STATE,
                "cannot verify particle boundary %d on both paths at t=%.17g",
                crossing->surface_id, crossing->t);
            return -1;
        }
    }
    if (crossing->surface_id > 0 &&
        (!group.local_surface_complete || group.local_surface_count == 0)) {
        alea_set_error_detail(ALEA_ERR_INVALID_STATE,
            "particle boundary surface group is incomplete at t=%.17g",
            crossing->t);
        return -1;
    }
    if (crossing->surface_id <= 0 && group.local_surface_count > 0 &&
        !group.local_surface_complete) {
        alea_set_error_detail(ALEA_ERR_INVALID_STATE,
            "particle lattice crossing has an incomplete physical surface group");
        return -1;
    }
    *out_type = ALEA_BOUNDARY_TRANSMISSIVE;
    int contains_selected = crossing->surface_id <= 0;
    for (size_t i = 0; i < group.local_surface_count; i++) {
        const int id = group.local_surface_ids[i];
        if (id == crossing->surface_id) contains_selected = 1;
        const int index = alea_surface_find(navigator->sys, id);
        alea_boundary_type_t type;
        if (index < 0 || alea_surface_get(navigator->sys, (size_t)index,
                NULL, NULL, NULL, NULL, &type) != 0) {
            alea_set_error_detail(ALEA_ERR_INVALID_STATE,
                "particle boundary surface %d is unavailable", id);
            return -1;
        }
        if (type != ALEA_BOUNDARY_TRANSMISSIVE) {
            if (crossing->surface_id <= 0 ||
                group.local_surface_count != 1 ||
                !group.local_surface_complete) {
                alea_set_error_detail(ALEA_ERR_INVALID_STATE,
                    "particle boundary action is ambiguous at t=%.17g",
                    crossing->t);
                return -1;
            }
            *out_type = type;
        }
    }
    if (!contains_selected) {
        alea_set_error_detail(ALEA_ERR_INVALID_STATE,
            "selected particle boundary is absent from its surface group");
        return -1;
    }
    return 0;
}

/* An RPP has one surface ID for all six faces. Its normal routine chooses the
 * first nearest face, which cannot define a reflection or periodic map at an
 * edge or corner. Count incident faces before accepting a directional action. */
static int navigator_rpp_action_face_unique(
    const alea_ray_navigator_t* navigator,
    const alea_raycast_boundary_event_t* crossing) {
    if (crossing->primitive_id >=
        alea_vec_count(&navigator->sys->primitives)) return 0;
    const alea_primitive_entry_t* primitive =
        &navigator->sys->primitives.data[crossing->primitive_id];
    if (primitive->type != ALEA_PRIMITIVE_RPP) return 1;
    alea_primitive_data_t data;
    if (!raycast_primitive_copy_payload(navigator->sys, crossing->primitive_id,
                                        primitive->type, &data)) return 0;
    alea_ray_t local_ray;
    alea_ray_transform_inverse(&crossing->transform, &navigator->ray, &local_ray);
    double p[3];
    alea_ray_point_at(&local_ray, crossing->t, &p[0], &p[1], &p[2]);
    const double mins[3] = {
        data.box.min_x, data.box.min_y, data.box.min_z
    };
    const double maxs[3] = {
        data.box.max_x, data.box.max_y, data.box.max_z
    };
    int faces = 0;
    for (int axis = 0; axis < 3; axis++) {
        const double scale = fmax(1.0,
            fmax(fabs(p[axis]), fmax(fabs(mins[axis]), fabs(maxs[axis]))));
        const double tolerance = RAY_EPSILON + 64.0 * DBL_EPSILON * scale;
        if (fabs(p[axis] - mins[axis]) <= tolerance) faces++;
        if (fabs(p[axis] - maxs[axis]) <= tolerance) faces++;
    }
    return faces == 1;
}

/* A recovery must be bounded by physical crossings on the original ray.
 * Sweep a short window on both sides so the entry and exit are observed
 * breakpoints, not merely the clipped ends of a query. */
typedef struct {
    double lo, hi;
    alea_ray_coverage_kind_t defect;
    uint64_t incident, outgoing, overlap_other;
    bool before, inside, after;
} transport_coverage_t;

static uint64_t unique_terminal_key(const alea_ray_coverage_interval_t* row) {
    if (row->kind != ALEA_RAY_COVERAGE_UNIQUE || !row->owner_count) return 0;
    const alea_ray_coverage_owner_t* leaf = row->owners;
    for (size_t i = 1; i < row->owner_count; ++i)
        if (row->owners[i].depth > leaf->depth) leaf = row->owners + i;
    return leaf->occurrence_key;
}

static int verify_transport_row(void* context, const alea_ray_coverage_interval_t* row) {
    transport_coverage_t* p = context;
    if (row->t_exit == p->lo) p->before = unique_terminal_key(row) == p->incident;
    if (row->t_enter == p->hi) {
        p->outgoing = unique_terminal_key(row);
        p->after = p->outgoing != 0;
    }
    if (row->t_enter != p->lo || row->t_exit != p->hi || row->kind != p->defect)
        return 0;
    if (row->kind != ALEA_RAY_COVERAGE_OVERLAP) {
        p->inside = true;
        return 0;
    }
    size_t leaves = 0;
    bool has_incident = false;
    for (size_t i = 0; i < row->owner_count; ++i) {
        const alea_ray_coverage_owner_t* owner = &row->owners[i];
        bool parent = false;
        for (size_t j = 0; j < row->owner_count; ++j)
            if (row->owners[j].parent_occurrence_key == owner->occurrence_key) parent = true;
        if (parent) continue;
        if (owner->resolution_flags) return 0;
        ++leaves;
        if (owner->occurrence_key == p->incident) has_incident = true;
        else p->overlap_other = owner->occurrence_key;
    }
    p->inside = leaves == 2 && has_incident;
    return 0;
}

static bool transport_span_verified(alea_ray_navigator_t* nav,
        const alea_nav_location_t* incident, int enter_surface, uint64_t* outgoing) {
    const double lo = nav->pending.t_enter, hi = nav->pending.t_exit;
    const double pad = ALEA_NAV_TRANSPORT_RECOVERY_CM;
    const int exit_surface = nav->pending.exit_surface_id;
    if (!(hi > lo && hi - lo <= pad) || nav->window_end || exit_surface <= 0 ||
        !incident->occurrence_key || nav->substep_fill.enter_surface_id) return false;
    transport_coverage_t proof = {.lo = lo, .hi = hi, .incident = incident->occurrence_key};
    switch (nav->location.kind) {
    case ALEA_NAV_GAP: proof.defect = ALEA_RAY_COVERAGE_GAP; break;
    case ALEA_NAV_UNDEFINED_FILL: proof.defect = ALEA_RAY_COVERAGE_UNDEFINED_FILL; break;
    case ALEA_NAV_OVERLAP: proof.defect = ALEA_RAY_COVERAGE_OVERLAP; break;
    default: return false;
    }
    /* A zero entry ID denotes continuing a proven interval after scattering. */
    const double start = fmax(0, lo - pad);
    alea_raycast_result_t* raw = &nav->validation_scratch;
    raw->transport_fill_filter = true;
    int rows = alea_ray_coverage_sweep_strict_reuse_nocache(nav->sys, &nav->ray,
            start, hi + pad, nav->max_interval_breakpoints, NULL, raw,
            verify_transport_row, &proof);
    raw->transport_fill_filter = false;
    if (rows < 0) {
        alea_clear_error_detail();
        return false;
    }
    if (proof.defect == ALEA_RAY_COVERAGE_UNDEFINED_FILL &&
        !(nav->pending.resolution_flags & ALEA_RESOLVE_UNDEFINED_FILL)) return false;
    if (proof.defect == ALEA_RAY_COVERAGE_OVERLAP &&
        nav->pending.owner_occurrence_key != incident->occurrence_key &&
        nav->pending.owner_occurrence_key != proof.overlap_other) return false;
    if (!enter_surface && proof.overlap_other && proof.overlap_other != nav->recovery_seed.other)
        return false;
    if (!proof.inside || !proof.after || (enter_surface && !proof.before) ||
        (proof.overlap_other && proof.outgoing != proof.overlap_other &&
         (enter_surface || proof.outgoing != proof.incident))) return false;
    bool saw_entry = !enter_surface, saw_exit = false;
    for (size_t i = 0; i < raw->hits.count; ++i) {
        const alea_ray_hit_t* hit = &raw->hits.data[i];
        if (hit->t < lo || hit->t > hi) continue;
        if (hit->t > lo && hit->t < hi) return false;
        int index = alea_surface_find(nav->sys, hit->surface_id);
        if (index < 0 || nav->sys->surfaces.data[index].boundary_type != ALEA_BOUNDARY_TRANSMISSIVE)
            return false;
        if (hit->t == lo && hit->surface_id == enter_surface) saw_entry = true;
        if (hit->t == hi && hit->surface_id == exit_surface) saw_exit = true;
    }
    if (!saw_entry || !saw_exit) return false;
    alea_boundary_type_t condition;
    if (navigator_boundary_condition(nav, &nav->pending.exit_event, true, &condition) ||
        condition != ALEA_BOUNDARY_TRANSMISSIVE) {
        alea_clear_error_detail();
        return false;
    }
    *outgoing = proof.outgoing;
    return true;
}

void alea_nav_get_recovery_seed(const alea_ray_navigator_t* nav, alea_transport_recovery_seed_t* seed) {
    *seed = nav->recovered_outgoing ? nav->recovery_seed : (alea_transport_recovery_seed_t){0};
}

int alea_nav_resume_recovery(alea_ray_navigator_t* nav, const alea_transport_recovery_seed_t* seed,
                             alea_nav_location_t* location) {
    uint64_t outgoing;
    nav->recovery_seed = *seed;
    if (nav->location.kind != seed->defect ||
        !transport_span_verified(nav, &seed->incident, 0, &outgoing) ||
        (outgoing != seed->other && outgoing != seed->incident.occurrence_key)) {
        nav->initialized = false;
        alea_set_error_detail(ALEA_ERR_INVALID_STATE,
            "cannot verify the new flight from a recovered interval");
        return -1;
    }
    nav->recovered_outgoing = outgoing;
    nav->recovery_remaining = nav->pending.t_exit - nav->current_t;
    *location = nav->location = seed->incident;
    return 0;
}

int alea_ray_navigator_advance(alea_ray_navigator_t* navigator,
                               double collision_distance, double max_distance,
                               alea_nav_event_t* event) {
    if (!navigator || !navigator_index_current(navigator) || !event || !navigator->initialized || navigator->leaked ||
        navigator->awaiting_action ||
        !(collision_distance > 0.0) || isnan(collision_distance) ||
        !isfinite(max_distance) || !(max_distance > 0.0)) return -1;
    if (navigator->location.kind != ALEA_NAV_MATERIAL &&
        navigator->location.kind != ALEA_NAV_VOID &&
        navigator->location.kind != ALEA_NAV_GAP) {
        alea_set_error_detail(ALEA_ERR_INVALID_STATE,
                              "particle navigator has ambiguous or undefined ownership");
        return -1;
    }
    if (navigator->location.kind != ALEA_NAV_MATERIAL &&
        isfinite(collision_distance)) {
        alea_set_error_detail(ALEA_ERR_INVALID_ARG,
                              "collision distance requires a material owner");
        return -1;
    }
    if (!navigator->have_pending && navigator_load_interval(navigator) != 0)
        return -1;
    while (navigator->window_end && navigator->pending.t_exit - navigator->current_t <=
            fmin(max_distance, collision_distance)) {
        const double start = navigator->pending.t_enter, end = navigator->pending.t_exit;
        const double limit = end + (end - start);
        const uint64_t owner = navigator->pending.owner_occurrence_key;
        const int cell = navigator->pending.cell_index;
        if (!(limit > end) || !isfinite(limit) ||
            navigator_regional_span(navigator, start, limit) ||
            navigator->pending.t_exit < end || navigator->pending.cell_index != cell ||
            navigator->pending.owner_occurrence_key != owner) {
            navigator->initialized = false;
            alea_set_error_detail(ALEA_ERR_INVALID_STATE, "Regional navigator cannot extend the certified span");
            return -1;
        }
    }
    const double boundary_distance = navigator->recovered_outgoing ? navigator->recovery_remaining :
        navigator->pending.t_exit - navigator->current_t;
    if (!(boundary_distance > 0.0)) {
        alea_set_error_detail(ALEA_ERR_INVALID_STATE,
                              "particle navigator boundary did not advance");
        return -1;
    }
    const bool reaches_boundary = boundary_distance <= max_distance &&
                                  boundary_distance <= collision_distance;
    if (navigator->validation_mode == ALEA_NAV_VALIDATE_INTERVAL) {
        /* Reconstructing the absolute root from its rounded flight distance
         * can move the validation endpoint one ULP beyond the boundary. */
        double end = reaches_boundary ? navigator->pending.t_exit
            : navigator->current_t + fmin(collision_distance, max_distance);
        if (navigator_validate_flight(navigator, end) != 0)
            return -1;
    }
    const uint64_t recovered_outgoing = navigator->recovered_outgoing;
    alea_nav_event_t next = {0};
    next.before = navigator->location;
    next.after = navigator->location;
    next.surface_id = -1;
    next.boundary_type = ALEA_BOUNDARY_TRANSMISSIVE;
    if (reaches_boundary) {
        next.kind = ALEA_NAV_BOUNDARY;
        next.distance = boundary_distance;
        const alea_raycast_boundary_event_t* crossing =
            &navigator->pending.exit_event;
        next.surface_id = crossing->surface_id;
        if (navigator_boundary_condition(navigator, crossing,
                recovered_outgoing != 0, &next.boundary_type) != 0) return -1;
        if (next.surface_id > 0 &&
            next.boundary_type != ALEA_BOUNDARY_TRANSMISSIVE &&
            next.boundary_type != ALEA_BOUNDARY_VACUUM &&
            !navigator_rpp_action_face_unique(navigator, crossing)) {
            alea_set_error_detail(ALEA_ERR_INVALID_STATE,
                "particle boundary action has no unique face normal at t=%.17g",
                crossing->t);
            return -1;
        }
        if (next.surface_id > 0) {
            if (next.boundary_type == ALEA_BOUNDARY_VACUUM)
                next.kind = ALEA_NAV_VACUUM;
            else if (next.boundary_type != ALEA_BOUNDARY_TRANSMISSIVE)
                next.kind = ALEA_NAV_BOUNDARY_ACTION;
            if ((navigator->event_fields & ALEA_NAV_EVENT_NORMAL) ||
                next.boundary_type == ALEA_BOUNDARY_REFLECTIVE ||
                next.boundary_type == ALEA_BOUNDARY_WHITE) {
                alea_ray_first_visible_result_t visible = {0};
                boundary_event_first_visible_normal(navigator->sys,
                    &navigator->ray, crossing, true, &visible);
                next.normal[0] = visible.nx;
                next.normal[1] = visible.ny;
                next.normal[2] = visible.nz;
            }
        }
        navigator->current_t = navigator->pending.t_exit;
        navigator->have_pending = false;
        if (next.kind == ALEA_NAV_VACUUM) {
            navigator->leaked = true;
            next.after.kind = ALEA_NAV_LEAKED;
            navigator->location = next.after;
        } else if (next.kind == ALEA_NAV_BOUNDARY_ACTION) {
            navigator->awaiting_action = true;
            navigator->action_type = next.boundary_type;
            navigator->action_surface_id = next.surface_id;
            memcpy(navigator->action_normal, next.normal,
                   sizeof(navigator->action_normal));
        } else {
            const int incident = navigator->pending.cell_index;
            const uint64_t parent = navigator->pending.owner_parent_occurrence_key;
            navigator->substep_fill.enter_surface_id = 0;
            if (navigator_load_interval(navigator) != 0) {
                navigator->initialized = false;
                return -1;
            }
            if (navigator->validation_mode != ALEA_NAV_VALIDATE_FAST && parent &&
                parent == navigator->pending.owner_parent_occurrence_key &&
                next.surface_id > 0) {
                if (alea_substep_fill_transition(navigator->sys, &navigator->ray,
                        navigator->current_t, incident, next.surface_id,
                        navigator->pending.path, &navigator->substep_fill)) {
                    navigator->location.kind = ALEA_NAV_UNRESOLVED;
                }
            }
            if (recovered_outgoing) {
                navigator->recovered_outgoing = 0;
                if ((navigator->location.kind != ALEA_NAV_MATERIAL && navigator->location.kind != ALEA_NAV_VOID) ||
                    navigator->location.occurrence_key != recovered_outgoing) {
                    navigator->initialized = false;
                    alea_set_error_detail(ALEA_ERR_INVALID_STATE, "recovered crossing has inconsistent outgoing owner");
                    return -1;
                }
            } else if (navigator->validation_mode == ALEA_NAV_VALIDATE_TRANSPORT &&
                       (next.before.kind == ALEA_NAV_MATERIAL || next.before.kind == ALEA_NAV_VOID) &&
                       next.surface_id > 0 && transport_span_verified(navigator, &next.before,
                           next.surface_id, &navigator->recovered_outgoing)) {
                if (navigator->location.kind == ALEA_NAV_OVERLAP) navigator->recovery_stats.overlaps++;
                else navigator->recovery_stats.gaps++;
                navigator->recovery_seed = (alea_transport_recovery_seed_t){
                    .incident = next.before, .other = navigator->recovered_outgoing,
                    .defect = navigator->location.kind
                };
                navigator->recovery_remaining = navigator->pending.t_exit - navigator->current_t;
                navigator->recovery_stats.last_enter_surface_id = next.surface_id;
                navigator->recovery_stats.last_exit_surface_id = navigator->pending.exit_surface_id;
                navigator->location = next.before;
            }
            next.after = navigator->location;
        }
    } else {
        next.kind = collision_distance <= max_distance
            ? ALEA_NAV_COLLISION : ALEA_NAV_DISTANCE_LIMIT;
        next.distance = fmin(collision_distance, max_distance);
        if (recovered_outgoing) {
            navigator->recovery_remaining -= next.distance;
            navigator->current_t = navigator->pending.t_exit - navigator->recovery_remaining;
        } else navigator->current_t += next.distance;
        navigator->substep_fill.enter_surface_id = 0;
    }
    alea_ray_point_at(&navigator->ray, navigator->current_t,
                      &next.position[0], &next.position[1], &next.position[2]);
    if (recovered_outgoing) navigator->recovery_stats.distance_cm += next.distance;
    *event = next;
    return 0;
}

/* Stream selected intervals into the legacy segment vector without routing
 * normal no-hit production traces through the combined compatibility loop. */
static int ray_selected_interval_trace_mode(
    alea_system_t* sys, const alea_ray_t* ray, double t_max,
    bool emit_hits, bool collect_boundary_groups, alea_raycast_result_t* result,
    const alea_ray_root_lookup_t* lookup) {
    alea_ray_walk_t walk;
    alea_ray_walk_init(&walk);
    /* This consumer walks the entire segment. Prepare reusable roots through
     * its finite end even when the first retry searches only a short span. */
    result->recovery_root_limit = t_max < DBL_MAX ? t_max : 0;
    walk.root_lookup = lookup;
    int published_prev_cell_idx = -2;
    int pending_enter_hit_index = -1;
    int status = 0;
    for (;;) {
        alea_ray_selected_interval_t interval;
        const int rc = alea_ray_walk_next_selected_mode(
            sys, ray, t_max, result, &walk, &interval,
            collect_boundary_groups);
        if (rc < 0) { status = -1; break; }
        if (rc == 2) break;
        if (ray_selected_interval_publish(
                sys, ray, &interval, true, emit_hits, result,
                &published_prev_cell_idx, &pending_enter_hit_index) != 0) {
            status = -1;
            break;
        }
        walk.pending_enter_hit_index = pending_enter_hit_index;
        if (rc == 0) break;
    }
    recovery_context_destroy(result->recovery_roots);
    result->recovery_roots = NULL;
    result->recovery_root_limit = 0;
    return status;
}

static int ray_selected_interval_trace(
    alea_system_t* sys, const alea_ray_t* ray, double t_max,
    bool emit_hits, alea_raycast_result_t* result) {
    return ray_selected_interval_trace_mode(
        sys, ray, t_max, emit_hits, true, result, NULL);
}

int alea_raycast_selected_boundary_events_with_options_nocache(
    alea_system_t* sys, const alea_ray_t* ray, double t_max,
    const alea_ray_boundary_event_options_internal_t* options,
    alea_raycast_result_t* scratch, alea_ray_boundary_event_result_t* events) {
    if (!sys || !ray || !scratch || !events || !(t_max > 0.0)) return -1;
    alea_ray_boundary_event_result_clear(events);
    alea_raycast_result_clear(scratch);
    alea_ray_walk_t walk;
    alea_ray_walk_init(&walk);
    alea_ray_selected_interval_t previous = {0};
    int have_previous = 0;
    for (;;) {
        alea_ray_selected_interval_t current;
        const int rc = alea_ray_walk_next_selected(
            sys, ray, t_max, scratch, &walk, &current);
        if (rc < 0) return -1;
        if (rc == 2) return 0;
        const int occurrence_changed =
            previous.owner_provenance_complete &&
            current.owner_provenance_complete &&
            previous.owner_occurrence_key != current.owner_occurrence_key;
        if (have_previous &&
            fabs(previous.t_exit - current.t_enter) <= RAY_EPSILON &&
            (previous.cell_id != current.cell_id ||
             previous.material_id != current.material_id ||
             occurrence_changed)) {
            const alea_raycast_boundary_event_t* source = &previous.exit_event;
            alea_ray_boundary_event_t event = {
                .t = previous.t_exit,
                .kind = source->has_physical_surface
                    ? ALEA_RAY_BOUNDARY_EVENT_PHYSICAL
                    : (source->is_synthetic_lattice_boundary
                        ? ALEA_RAY_BOUNDARY_EVENT_SYNTHETIC_LATTICE
                        : ALEA_RAY_BOUNDARY_EVENT_UNRESOLVED),
                .surface_id = source->surface_id,
                .primitive_id = source->primitive_id,
                .cell_before = previous.cell_id,
                .cell_after = current.cell_id,
                .material_before = previous.material_id,
                .material_after = current.material_id,
                .resolution_flags = previous.resolution_flags |
                                    current.resolution_flags,
                .local_surface_count = source->local_surface_count,
                .local_surface_complete = source->local_surface_complete,
                .provenance_flags = source->provenance_flags,
                .active_cell_id = source->active_cell_id,
                .active_universe_id = source->active_universe_id,
                .active_depth = source->active_depth,
                .active_occurrence_key = source->active_occurrence_key,
                .active_parent_occurrence_key =
                    source->active_parent_occurrence_key
            };
            memcpy(event.local_surface_ids, source->local_surface_ids,
                   sizeof(event.local_surface_ids));
            memcpy(event.local_point, source->local_point,
                   sizeof(event.local_point));
            memcpy(event.local_direction, source->local_direction,
                   sizeof(event.local_direction));
            alea_ray_first_visible_result_t normal = {0};
            boundary_event_first_visible_normal(
                sys, ray, source, true, &normal);
            event.nx = normal.nx;
            event.ny = normal.ny;
            event.nz = normal.nz;
            if (previous.owner_provenance_complete) {
                event.before_occurrence_key = previous.owner_occurrence_key;
                event.before_parent_occurrence_key =
                    previous.owner_parent_occurrence_key;
                event.provenance_flags |= ALEA_BOUNDARY_PROVENANCE_BEFORE_OWNER;
            }
            if (current.owner_provenance_complete) {
                event.after_occurrence_key = current.owner_occurrence_key;
                event.after_parent_occurrence_key =
                    current.owner_parent_occurrence_key;
                event.provenance_flags |= ALEA_BOUNDARY_PROVENANCE_AFTER_OWNER;
            }
            if (!(options && options->skip_open_side_coverage))
                boundary_event_require_unambiguous_open_sides(
                    sys, ray, &previous, &current, &event);
            if (options &&
                ((options->max_events &&
                  events->events.count >= options->max_events) ||
                 (options->max_output_bytes &&
                  events->events.count >=
                    options->max_output_bytes / sizeof(*events->events.data)))) {
                alea_set_error_detail(
                    ALEA_ERR_OVERFLOW,
                    "selected boundary-event query output limit exceeded");
                return -1;
            }
            if (alea_vec_push(&events->events, event,
                              alea_ray_boundary_event_t) != 0)
                return -1;
        }
        previous = current;
        have_previous = 1;
        if (rc == 0) return 0;
    }
}

int alea_raycast_selected_boundary_events_nocache(
    alea_system_t* sys, const alea_ray_t* ray, double t_max,
    alea_raycast_result_t* scratch, alea_ray_boundary_event_result_t* events) {
    return alea_raycast_selected_boundary_events_with_options_nocache(
        sys, ray, t_max, NULL, scratch, events);
}

static int selected_boundary_groups_equal(
    const alea_ray_boundary_event_t* a,
    const alea_ray_boundary_event_t* b) {
    if (a->local_surface_count != b->local_surface_count) return 0;
    for (size_t i = 0; i < a->local_surface_count; i++) {
        int found = 0;
        for (size_t j = 0; j < b->local_surface_count; j++)
            if (a->local_surface_ids[i] == b->local_surface_ids[j]) {
                found = 1;
                break;
            }
        if (!found) return 0;
    }
    return 1;
}

int alea_raycast_selected_boundary_events_bidirectional_nocache(
    alea_system_t* sys, const alea_ray_t* ray, double t_max,
    alea_raycast_result_t* forward_scratch,
    alea_raycast_result_t* reverse_scratch,
    alea_ray_boundary_event_result_t* events) {
    if (!sys || !ray || !forward_scratch || !reverse_scratch || !events ||
        !(t_max > 0.0)) return -1;
    alea_ray_boundary_event_result_t forward, reverse;
    alea_ray_boundary_event_result_init(&forward);
    alea_ray_boundary_event_result_init(&reverse);
    double ex, ey, ez;
    alea_ray_point_at(ray, t_max, &ex, &ey, &ez);
    alea_ray_t reverse_ray;
    if (alea_ray_init(&reverse_ray, ex, ey, ez,
                      -ray->dx, -ray->dy, -ray->dz) != 0 ||
        alea_raycast_selected_boundary_events_nocache(
            sys, ray, t_max, forward_scratch, &forward) != 0 ||
        alea_raycast_selected_boundary_events_nocache(
            sys, &reverse_ray, t_max, reverse_scratch, &reverse) != 0)
        goto fail;
    if (forward.events.count != reverse.events.count) goto incomplete;
    alea_ray_boundary_event_result_clear(events);
    for (size_t i = 0; i < forward.events.count; i++) {
        const alea_ray_boundary_event_t* a = &forward.events.data[i];
        const alea_ray_boundary_event_t* b =
            &reverse.events.data[reverse.events.count - 1 - i];
        if (a->kind != ALEA_RAY_BOUNDARY_EVENT_PHYSICAL ||
            b->kind != ALEA_RAY_BOUNDARY_EVENT_PHYSICAL ||
            fabs(a->t - (t_max - b->t)) > RAY_EPSILON ||
            a->cell_before != b->cell_after || a->cell_after != b->cell_before ||
            a->material_before != b->material_after ||
            a->material_after != b->material_before)
            goto incomplete;
        alea_ray_boundary_event_t merged = *a;
        const int terminal_to_void = merged.cell_before < 0 || merged.cell_after < 0;
        if (a->local_surface_complete && b->local_surface_complete &&
            !selected_boundary_groups_equal(a, b))
            goto incomplete;
        if (!a->local_surface_complete && b->local_surface_complete) {
            merged.local_surface_count = b->local_surface_count;
            memcpy(merged.local_surface_ids, b->local_surface_ids,
                   sizeof(merged.local_surface_ids));
        }
        merged.local_surface_complete = terminal_to_void
            ? (a->local_surface_complete || b->local_surface_complete)
            : (a->local_surface_complete && b->local_surface_complete);
        if (!merged.local_surface_complete || merged.local_surface_count == 0)
            goto incomplete;
        if (alea_vec_push(&events->events, merged,
                          alea_ray_boundary_event_t) != 0)
            goto fail;
    }
    alea_ray_boundary_event_result_free(&forward);
    alea_ray_boundary_event_result_free(&reverse);
    return 0;
incomplete:
    alea_ray_boundary_event_result_clear(events);
    alea_ray_boundary_event_result_free(&forward);
    alea_ray_boundary_event_result_free(&reverse);
    return 1;
fail:
    alea_ray_boundary_event_result_clear(events);
    alea_ray_boundary_event_result_free(&forward);
    alea_ray_boundary_event_result_free(&reverse);
    return -1;
}

static void selected_event_transition_receipt(
    const alea_ray_boundary_event_t* event,
    alea_transition_result_t* result) {
    result->occurrence_depth = event->active_depth;
    result->current_occurrence_key = event->active_occurrence_key;
    result->current_parent_occurrence_key =
        event->active_parent_occurrence_key;
    result->before_occurrence_key = event->before_occurrence_key;
    result->before_parent_occurrence_key =
        event->before_parent_occurrence_key;
    result->selected_after_occurrence_key = event->after_occurrence_key;
    result->selected_after_parent_occurrence_key =
        event->after_parent_occurrence_key;
}

int alea_check_selected_boundary_event_transition_reuse_nocache(
    alea_system_t* sys, const alea_ray_boundary_event_t* event,
    const alea_transition_options_t* options,
    alea_transition_result_t* result,
    alea_transition_workspace_t* workspace) {
    if (!sys || !event || !result || !workspace) return -1;
    if (event->kind != ALEA_RAY_BOUNDARY_EVENT_PHYSICAL ||
        event->surface_id <= 0) {
        alea_set_error_detail(
            ALEA_ERR_INVALID_ARG,
            "selected boundary transition requires a physical event");
        return -1;
    }
    if (!(event->provenance_flags & ALEA_BOUNDARY_PROVENANCE_ACTIVE_FRAME)) {
        memset(result, 0, sizeof(*result));
        result->kind = ALEA_TRANSITION_UNRESOLVED;
        result->after_coverage_kind = ALEA_POINT_COVERAGE_UNRESOLVED;
        result->current_cell_id = event->cell_before;
        result->primary_surface_id = event->surface_id;
        result->connecting_surface_id = -1;
        result->after_cell_id = event->cell_after;
        result->occurrence_depth = -1;
        return 0;
    }
    if (!event->local_surface_complete || event->local_surface_count == 0) {
        memset(result, 0, sizeof(*result));
        result->kind = ALEA_TRANSITION_AMBIGUOUS_BOUNDARY;
        result->after_coverage_kind = ALEA_POINT_COVERAGE_UNRESOLVED;
        result->universe_id = event->active_universe_id;
        result->current_cell_id = event->active_cell_id;
        result->primary_surface_id = event->surface_id;
        result->connecting_surface_id = -1;
        result->after_cell_id = event->cell_after;
        memcpy(result->crossing_point, event->local_point,
               sizeof(result->crossing_point));
        memcpy(result->direction, event->local_direction,
               sizeof(result->direction));
        selected_event_transition_receipt(event, result);
        return 0;
    }
    int rc = alea_check_transition_local_reuse(
        sys, event->active_universe_id, event->active_cell_id,
        event->surface_id, event->local_surface_ids,
        event->local_surface_count, event->local_point,
        event->local_direction, options, result, workspace);
    if (rc == 0) selected_event_transition_receipt(event, result);
    return rc;
}

int alea_check_selected_boundary_event_transition_nocache(
    alea_system_t* sys, const alea_ray_boundary_event_t* event,
    const alea_transition_options_t* options,
    alea_transition_result_t* result) {
    alea_transition_workspace_t workspace;
    alea_transition_workspace_init(&workspace);
    const int rc = alea_check_selected_boundary_event_transition_reuse_nocache(
        sys, event, options, result, &workspace);
    alea_transition_workspace_free(&workspace);
    return rc;
}

static int alea_ray_walk_next_first_visible(
    alea_system_t* sys, const alea_ray_t* ray, double t_max,
    double t_min, int material_filter, bool include_normal,
    alea_raycast_result_t* scratch, alea_ray_walk_t* walk,
    alea_ray_first_visible_result_t* visible) {
    alea_ray_selected_interval_t interval;
    const int rc = alea_ray_walk_next_selected(
        sys, ray, t_max, scratch, walk, &interval);
    if (rc < 0) return rc;
    if (rc == 2) return 0;
    if (ray_selected_interval_first_visible(
            sys, ray, &interval, t_min, material_filter, include_normal,
            visible))
        return 0;
    return rc == 0 ? 0 : 1;
}

#ifndef RAYCAST_FIRST_VISIBLE_PACKET_WIDTH
#define RAYCAST_FIRST_VISIBLE_PACKET_WIDTH 4
#endif

typedef struct {
    alea_ray_t ray;
    alea_ray_walk_t state;
    alea_raycast_result_t scratch;
    alea_ray_first_visible_result_t visible;
    size_t ray_index;
    uint8_t active;
} raycast_first_visible_lane_t;

static void raycast_first_visible_lane_init(
    raycast_first_visible_lane_t* lane,
    size_t ray_index,
    const double* origins_xyz,
    const double* directions_xyz) {
    memset(lane, 0, sizeof(*lane));
    lane->ray_index = ray_index;
    const double* o = &origins_xyz[ray_index * 3];
    const double* d = &directions_xyz[ray_index * 3];
    if (alea_ray_init(&lane->ray, o[0], o[1], o[2],
                      d[0], d[1], d[2]) != 0) {
        return;
    }
    alea_ray_walk_init(&lane->state);
    alea_raycast_result_init(&lane->scratch);
    lane->scratch.ray = lane->ray;
    lane->visible.cell_id = -1;
    lane->visible.surface_id = -1;
    lane->visible.primitive_id = UINT32_MAX;
    lane->active = 1;
}

static void raycast_first_visible_lane_store(
    const raycast_first_visible_lane_t* lane,
    alea_ray_first_visible_batch_result_t* result) {
    const size_t i = lane->ray_index;
    const alea_ray_first_visible_result_t* visible = &lane->visible;
    result->found[i] = visible->found ? 1 : 0;
    result->t[i] = visible->t;
    result->cell_ids[i] = visible->cell_id;
    result->material_ids[i] = visible->material_id;
    if (result->densities) result->densities[i] = visible->density;
    if (result->surface_ids) result->surface_ids[i] = visible->surface_id;
    if (result->primitive_ids) result->primitive_ids[i] = visible->primitive_id;
    if (result->resolution_flags)
        result->resolution_flags[i] = visible->resolution_flags;
    if (result->normals_xyz) {
        result->normals_xyz[i * 3] = visible->nx;
        result->normals_xyz[i * 3 + 1] = visible->ny;
        result->normals_xyz[i * 3 + 2] = visible->nz;
    }
}

int alea_raycast_hier_first_visible_batch_execute_nocache(
    alea_system_t* sys, const double* origins_xyz, const double* directions_xyz,
    size_t ray_count, const alea_ray_batch_query_t* query,
    alea_ray_first_visible_batch_result_t* result, int* statuses,
    size_t begin, size_t end) {
    if (!sys || !origins_xyz || !directions_xyz || !query || !result ||
        !statuses) {
        return -1;
    }

    if (begin > end || end > ray_count) return -1;
    for (size_t packet_start = begin; packet_start < end;
         packet_start += RAYCAST_FIRST_VISIBLE_PACKET_WIDTH) {
        raycast_first_visible_lane_t lanes[RAYCAST_FIRST_VISIBLE_PACKET_WIDTH];
        const size_t remaining = end - packet_start;
        const size_t packet_end = packet_start +
            (remaining < RAYCAST_FIRST_VISIBLE_PACKET_WIDTH
                 ? remaining : RAYCAST_FIRST_VISIBLE_PACKET_WIDTH);
        size_t active_count = 0;
        for (size_t i = packet_start; i < packet_end; i++) {
            raycast_first_visible_lane_init(&lanes[i - packet_start], i,
                                            origins_xyz, directions_xyz);
            if (!lanes[i - packet_start].active) {
                statuses[i] = -1;
                continue;
            }
            active_count++;
        }

        while (active_count > 0 && !alea_interrupted()) {
            for (size_t slot = 0; slot < packet_end - packet_start; slot++) {
                raycast_first_visible_lane_t* lane = &lanes[slot];
                if (!lane->active) continue;
                const size_t i = lane->ray_index;
                const double t_min = query->t_mins ? query->t_mins[i] : 0.0;
                const double t_max = query->t_maxs ? query->t_maxs[i] : 0.0;
                const double effective_t_max = t_max <= 0.0 ? DBL_MAX : t_max;
                const int rc = alea_ray_walk_next_first_visible(
                    sys, &lane->ray, effective_t_max, t_min,
                    query->material_filter,
                    (query->fields & ALEA_RAY_QUERY_FIELD_SURFACE_NORMAL) != 0,
                    &lane->scratch, &lane->state, &lane->visible);
                if (rc < 0) {
                    statuses[i] = -1;
                    lane->active = 0;
                    active_count--;
                    continue;
                }
                if (rc == 0) {
                    raycast_first_visible_lane_store(lane, result);
                    lane->active = 0;
                    active_count--;
                }
            }
        }
        for (size_t slot = 0; slot < packet_end - packet_start; slot++) {
            raycast_first_visible_lane_t* lane = &lanes[slot];
            if (lane->active) {
                statuses[lane->ray_index] = -1;
                lane->active = 0;
            }
            alea_raycast_result_free(&lane->scratch);
        }
    }
    return alea_interrupted() ? -1 : 0;
}

int alea_raycast_hier_any_hit_batch_execute_nocache(
    alea_system_t* sys, const double* origins_xyz, const double* directions_xyz,
    size_t ray_count, const alea_ray_batch_query_t* query,
    alea_ray_any_hit_batch_result_t* result, int* statuses,
    size_t begin, size_t end) {
    if (!sys || !origins_xyz || !directions_xyz || !query || !result ||
        !statuses) {
        return -1;
    }

    if (begin > end || end > ray_count) return -1;
    for (size_t packet_start = begin; packet_start < end;
         packet_start += RAYCAST_FIRST_VISIBLE_PACKET_WIDTH) {
        raycast_first_visible_lane_t lanes[RAYCAST_FIRST_VISIBLE_PACKET_WIDTH];
        const size_t remaining = end - packet_start;
        const size_t packet_end = packet_start +
            (remaining < RAYCAST_FIRST_VISIBLE_PACKET_WIDTH
                 ? remaining : RAYCAST_FIRST_VISIBLE_PACKET_WIDTH);
        size_t active_count = 0;
        for (size_t i = packet_start; i < packet_end; i++) {
            raycast_first_visible_lane_init(&lanes[i - packet_start], i,
                                            origins_xyz, directions_xyz);
            if (!lanes[i - packet_start].active) {
                statuses[i] = -1;
                continue;
            }
            active_count++;
        }

        while (active_count > 0 && !alea_interrupted()) {
            for (size_t slot = 0; slot < packet_end - packet_start; slot++) {
                raycast_first_visible_lane_t* lane = &lanes[slot];
                if (!lane->active) continue;
                const size_t i = lane->ray_index;
                const double t_min = query->t_mins ? query->t_mins[i] : 0.0;
                const double t_max = query->t_maxs ? query->t_maxs[i] : 0.0;
                const double effective_t_max = t_max <= 0.0 ? DBL_MAX : t_max;
                const int rc = alea_ray_walk_next_first_visible(
                    sys, &lane->ray, effective_t_max, t_min,
                    query->material_filter, false,
                    &lane->scratch, &lane->state, &lane->visible);
                if (rc < 0) {
                    statuses[i] = -1;
                    lane->active = 0;
                    active_count--;
                    continue;
                }
                if (rc == 0) {
                    result->hits[i] = lane->visible.found ? 1 : 0;
                    lane->active = 0;
                    active_count--;
                }
            }
        }
        for (size_t slot = 0; slot < packet_end - packet_start; slot++) {
            raycast_first_visible_lane_t* lane = &lanes[slot];
            if (lane->active) {
                statuses[lane->ray_index] = -1;
                lane->active = 0;
            }
            alea_raycast_result_free(&lane->scratch);
        }
    }
    return alea_interrupted() ? -1 : 0;
}

int alea_raycast_cell_aware(alea_system_t* sys,
                           double ox, double oy, double oz,
                           double dx, double dy, double dz,
                           double t_max,
                           alea_raycast_result_t* result) {
    if (!sys || !result) return -1;

    if (system_has_lattice_cells(sys)) {
        /* Lattice transport requires synthetic DDA boundary hits and
         * element-local universe raycasts. Use the canonical lattice-aware
         * pipeline so this public entry point cannot diverge semantically. */
        return alea_raycast(sys, ox, oy, oz, dx, dy, dz, t_max, result);
    }

    alea_raycast_result_free(result);

    alea_ray_t ray;
    if (alea_ray_init(&ray, ox, oy, oz, dx, dy, dz) != 0) {
        return -1;
    }
    result->ray = ray;

    double effective_t_max = (t_max <= 0) ? DBL_MAX : t_max;

    if (!alea_system_query_cache_ready(sys, ALEA_CACHE_RAYCAST)) {
        alea_set_error_detail(ALEA_ERR_INVALID_STATE,
                              "raycast caches are not prepared; call alea_prepare_query_acceleration()");
        return -1;
    }

    return raycast_flat_cell_aware_impl(sys, &ray, effective_t_max, result);
}

int alea_raycast_hier_cell_aware(alea_system_t* sys,
                                 double ox, double oy, double oz,
                                 double dx, double dy, double dz,
                                 double t_max,
                                 alea_raycast_result_t* result) {
    return alea_raycast_hier_fast_segments(sys, ox, oy, oz,
                                           dx, dy, dz, t_max,
                                           result);
}

int alea_raycast_hier_fast_segments(alea_system_t* sys,
                                    double ox, double oy, double oz,
                                    double dx, double dy, double dz,
                                    double t_max,
                                    alea_raycast_result_t* result) {
    return alea_raycast_hier_segments_with_root_lookup(
        sys, ox, oy, oz, dx, dy, dz, t_max, result, NULL);
}

int alea_raycast_hier_segments_with_root_lookup(
    alea_system_t* sys, double ox, double oy, double oz,
    double dx, double dy, double dz, double t_max, alea_raycast_result_t* result,
    const alea_ray_root_lookup_t* lookup) {
    if (!sys || !result) return -1;

    if (alea_system_prepare_query_caches(sys,
            ALEA_CACHE_HIER_SPATIAL | ALEA_CACHE_CELL_SURFACES) != 0)
        return -1;

    alea_raycast_result_free(result);

    alea_ray_t ray;
    if (alea_ray_init(&ray, ox, oy, oz, dx, dy, dz) != 0) {
        return -1;
    }
    result->ray = ray;

    double effective_t_max = (t_max <= 0) ? DBL_MAX : t_max;
    return ray_selected_interval_trace_mode(sys, &ray, effective_t_max, false, true, result, lookup);
}

int alea_raycast_hier_with_hits(alea_system_t* sys,
                                double ox, double oy, double oz,
                                double dx, double dy, double dz,
                                double t_max,
                                alea_raycast_result_t* result) {
    if (!sys || !result) return -1;

    if (alea_system_prepare_query_caches(sys,
            ALEA_CACHE_HIER_SPATIAL | ALEA_CACHE_CELL_SURFACES) != 0)
        return -1;

    alea_raycast_result_free(result);

    alea_ray_t ray;
    if (alea_ray_init(&ray, ox, oy, oz, dx, dy, dz) != 0) {
        return -1;
    }
    result->ray = ray;

    double effective_t_max = (t_max <= 0) ? DBL_MAX : t_max;
    return ray_selected_interval_trace(sys, &ray, effective_t_max, true, result);
}

/* Buffer-reuse hierarchical variants: take a pre-normalized ray, assume query
 * caches are already built (ALEA_CACHE_RAYCAST), and do NOT free the result
 * buffers. The caller must clear the result (count reset, capacity retained)
 * between rays. These are the hot-loop entry points for renderers that trace
 * one ray per pixel and must avoid per-pixel malloc/free churn. */
int alea_raycast_hier_with_hits_nocache(alea_system_t* sys,
                                        const alea_ray_t* ray,
                                        double t_max,
                                        alea_raycast_result_t* result) {
    if (!sys || !ray || !result) return -1;
    double effective_t_max = (t_max <= 0) ? DBL_MAX : t_max;
    result->ray = *ray;
    return ray_selected_interval_trace(sys, ray, effective_t_max, true, result);
}

int alea_raycast_hier_first_visible_nocache(
    alea_system_t* sys, const alea_ray_t* ray, double t_min, double t_max,
    int material_filter, int include_normal, alea_raycast_result_t* scratch,
    alea_ray_first_visible_result_t* out_visible) {
    if (!sys || !ray || !scratch || !out_visible || t_min < 0 ||
        (t_max > 0 && t_min > t_max))
        return -1;
    alea_raycast_result_clear(scratch);
    memset(out_visible, 0, sizeof(*out_visible));
    out_visible->cell_id = -1;
    out_visible->surface_id = -1;
    out_visible->primitive_id = UINT32_MAX;
    const double effective_t_max = t_max <= 0 ? DBL_MAX : t_max;
    scratch->ray = *ray;
    alea_ray_walk_t walk;
    alea_ray_walk_init(&walk);
    for (;;) {
        alea_ray_selected_interval_t interval;
        const int rc = alea_ray_walk_next_selected(
            sys, ray, effective_t_max, scratch, &walk, &interval);
        if (rc < 0) return -1;
        if (rc == 2) return 0;
        if (ray_selected_interval_first_visible(
                sys, ray, &interval, t_min, material_filter,
                include_normal != 0, out_visible))
            return 0;
        if (rc == 0) return 0;
    }
}

int alea_raycast_render_first_visible_nocache(
    alea_system_t* sys, const alea_ray_t* ray, double t_max,
    int include_normal, alea_raycast_render_accept_t accept, void* context,
    alea_raycast_result_t* scratch,
    alea_ray_first_visible_result_t* out_visible) {
    if (!sys || !ray || !scratch || !out_visible) return -1;
    alea_raycast_result_clear(scratch);
    memset(out_visible, 0, sizeof(*out_visible));
    out_visible->cell_id = -1;
    out_visible->surface_id = -1;
    out_visible->primitive_id = UINT32_MAX;
    scratch->ray = *ray;
    const double effective_t_max = t_max <= 0 ? DBL_MAX : t_max;
    alea_ray_walk_t walk;
    alea_ray_walk_init(&walk);
    /* A provisional void owner can hide a nearer material surface. Verify
     * ownership as for selected intervals, then stop at the first accepted
     * hit. Rendering still needs no coincident-surface group expansion. */
    walk.require_interior_sample = true;
    for (;;) {
        alea_ray_selected_interval_t interval;
        const int rc = alea_ray_walk_next_selected_mode(
            sys, ray, effective_t_max, scratch, &walk, &interval, false);
        if (rc < 0) return -1;
        if (rc == 2) return 0;
        if (interval.cell_id >= 0 && interval.material_id != 0 &&
            (!accept || accept(context, interval.cell_id,
                                interval.material_id, interval.density)) &&
            ray_selected_interval_first_visible(
                sys, ray, &interval, 0.0, -1, include_normal != 0,
                out_visible))
            return 0;
        if (rc == 0) return 0;
    }
}

int alea_raycast_hier_first_cell_nocache(
    alea_system_t* sys, const alea_ray_t* ray, double t_min, double t_max,
    int material_filter, alea_raycast_result_t* scratch,
    int* out_cell_id, double* out_t) {
    if (!sys || !ray || !scratch || !out_cell_id || t_min < 0 ||
        (t_max > 0 && t_min > t_max))
        return -1;
    *out_cell_id = -1;
    alea_raycast_result_clear(scratch);
    scratch->ray = *ray;
    const double effective_t_max = t_max <= 0 ? DBL_MAX : t_max;
    alea_ray_walk_t walk;
    alea_ray_walk_init(&walk);
    for (;;) {
        alea_ray_selected_interval_t interval;
        const int rc = alea_ray_walk_next_selected(
            sys, ray, effective_t_max, scratch, &walk, &interval);
        if (rc < 0) return -1;
        if (rc == 2) return 0;
        if (interval.cell_id >= 0 &&
            (material_filter < 0 || interval.material_id == material_filter) &&
            interval.t_exit > t_min + RAY_EPSILON) {
            *out_cell_id = interval.cell_id;
            if (out_t) *out_t = fmax(interval.t_enter, t_min);
            return 0;
        }
        if (rc == 0) return 0;
    }
}

int alea_raycast_hier_any_hit_nocache(
    alea_system_t* sys, const alea_ray_t* ray, double t_min, double t_max,
    int material_filter, alea_raycast_result_t* scratch, int* out_hit) {
    if (!sys || !ray || !scratch || !out_hit || t_min < 0 ||
        (t_max > 0 && t_min > t_max))
        return -1;
    *out_hit = 0;
    alea_raycast_result_clear(scratch);
    scratch->ray = *ray;
    const double effective_t_max = t_max <= 0 ? DBL_MAX : t_max;
    alea_ray_walk_t walk;
    alea_ray_walk_init(&walk);
    for (;;) {
        alea_ray_selected_interval_t interval;
        const int rc = alea_ray_walk_next_selected(
            sys, ray, effective_t_max, scratch, &walk, &interval);
        if (rc < 0) return -1;
        if (rc == 2) return 0;
        if (interval.cell_id >= 0 && interval.material_id != 0 &&
            (material_filter < 0 || interval.material_id == material_filter) &&
            interval.t_exit > t_min + RAY_EPSILON) {
            *out_hit = 1;
            return 0;
        }
        if (rc == 0) return 0;
    }
}

int alea_raycast_hier_segments_nocache(alea_system_t* sys,
                                       const alea_ray_t* ray,
                                       double t_max, bool collect_boundary_groups,
                                       alea_raycast_result_t* result) {
    if (!sys || !ray || !result) return -1;
    double effective_t_max = (t_max <= 0) ? DBL_MAX : t_max;
    result->ray = *ray;
    return ray_selected_interval_trace_mode(
        sys, ray, effective_t_max, false, collect_boundary_groups, result, NULL);
}

int alea_raycast_hier_visit_intervals_nocache(
    alea_system_t* sys, const alea_ray_t* ray, double t_max,
    alea_raycast_result_t* scratch,
    alea_raycast_selected_interval_callback_t callback, void* context) {
    if (!sys || !ray || !scratch || !callback) return -1;
    alea_raycast_result_clear(scratch);
    scratch->ray = *ray;
    const double effective_t_max = t_max <= 0 ? DBL_MAX : t_max;
    alea_ray_walk_t walk;
    alea_ray_walk_init(&walk);
    for (;;) {
        alea_ray_selected_interval_t interval;
        const int rc = alea_ray_walk_next_selected(
            sys, ray, effective_t_max, scratch, &walk, &interval);
        if (rc < 0) return -1;
        if (rc == 2) return 0;
        const alea_raycast_selected_interval_view_t view = {
            .t_enter = interval.t_enter,
            .containment_verified = interval.containment_verified,
            .t_exit = interval.t_exit,
            .cell_index = interval.cell_index,
            .cell_id = interval.cell_id,
            .material_id = interval.material_id,
            .density = interval.density,
            .resolution_flags = interval.resolution_flags,
            .owner_provenance_complete = interval.owner_provenance_complete,
            .owner_occurrence_key = interval.owner_occurrence_key,
            .owner_parent_occurrence_key =
                interval.owner_parent_occurrence_key,
            .path = interval.path
        };
        const int callback_rc = callback(context, &view);
        if (callback_rc < 0) return -1;
        if (callback_rc > 0 || rc == 0) return 0;
    }
}

typedef struct {
    alea_raycast_selected_segment_callback_t callback;
    void* context;
    alea_ray_segment_t pending;
    int pending_cell_index;
    bool have_pending;
    bool stopped;
} raycast_segment_adapter_t;

static int raycast_segment_adapter_visit(
        void* context,
        const alea_raycast_selected_interval_view_t* interval) {
    raycast_segment_adapter_t* adapter = context;
    const alea_ray_segment_t segment = {
        .t_enter = interval->t_enter,
        .t_exit = interval->t_exit,
        .cell_id = interval->cell_id,
        .material_id = interval->material_id,
        .density = interval->density,
        .enter_surface_id = -1,
        .exit_surface_id = -1,
        .enter_hit_index = -1,
        .resolution_flags = interval->resolution_flags,
        .path_index = UINT32_MAX
    };
    /* Match ray_selected_interval_publish(): adjacent intervals owned by the
     * same cell are one public segment when hierarchy paths are not captured.
     * Delay publication by one interval so the streaming and vector APIs keep
     * identical segment boundaries even when the walker observes a harmless
     * same-cell split at a tangent/coincident surface. */
    if (!adapter->have_pending) {
        adapter->pending = segment;
        adapter->pending_cell_index = interval->cell_index;
        adapter->have_pending = true;
        return 0;
    }
    if (adapter->pending_cell_index == interval->cell_index) {
        adapter->pending.t_exit = interval->t_exit;
        return 0;
    }
    const int callback_rc =
        adapter->callback(adapter->context, &adapter->pending);
    if (callback_rc != 0) {
        adapter->stopped = callback_rc > 0;
        return callback_rc;
    }
    adapter->pending = segment;
    adapter->pending_cell_index = interval->cell_index;
    return 0;
}

int alea_raycast_hier_visit_segments_nocache(
    alea_system_t* sys, const alea_ray_t* ray, double t_max,
    alea_raycast_result_t* scratch,
    alea_raycast_selected_segment_callback_t callback, void* context) {
    if (!callback) return -1;
    raycast_segment_adapter_t adapter = {
        .callback = callback,
        .context = context,
        .pending_cell_index = -2
    };
    const int rc = alea_raycast_hier_visit_intervals_nocache(
        sys, ray, t_max, scratch, raycast_segment_adapter_visit, &adapter);
    if (rc != 0) return rc;
    if (!adapter.have_pending || adapter.stopped) return 0;
    return adapter.callback(adapter.context, &adapter.pending) < 0 ? -1 : 0;
}

typedef struct {
    alea_raycast_batch_selected_segment_callback_t callback;
    void* context;
    size_t ray_index;
} raycast_batch_segment_bridge_t;

static int raycast_batch_segment_bridge(void* context,
                                        const alea_ray_segment_t* segment) {
    raycast_batch_segment_bridge_t* bridge = context;
    return bridge->callback(bridge->context, bridge->ray_index, segment);
}

typedef struct {
    alea_system_t* sys;
    const double* origins_xyz;
    const double* directions_xyz;
    double t_max;
    alea_raycast_batch_selected_segment_callback_t callback;
    void* callback_context;
    atomic_int* failed;
} raycast_visit_segments_batch_context_t;

static int raycast_visit_segments_batch_range(void* opaque,
                                              size_t worker_index,
                                              size_t begin, size_t end) {
    (void)worker_index;
    raycast_visit_segments_batch_context_t* context = opaque;
    for (size_t i = begin; i < end; i++) {
        if (atomic_load(context->failed) || alea_interrupted()) continue;
        const double* origin = &context->origins_xyz[i * 3];
        const double* direction = &context->directions_xyz[i * 3];
        alea_ray_t ray;
        alea_raycast_result_t scratch;
        alea_raycast_result_init(&scratch);
        raycast_batch_segment_bridge_t bridge = {
            .callback = context->callback,
            .context = context->callback_context,
            .ray_index = i
        };
        int rc = alea_ray_init(&ray, origin[0], origin[1], origin[2],
                               direction[0], direction[1], direction[2]);
        if (rc == 0) {
            rc = alea_raycast_hier_visit_segments_nocache(
                context->sys, &ray, context->t_max, &scratch,
                raycast_batch_segment_bridge, &bridge);
        }
        alea_raycast_result_free(&scratch);
        if (rc != 0) atomic_store(context->failed, 1);
    }
    return 0;
}

int alea_raycast_hier_visit_segments_batch_nocache(
    alea_system_t* sys, const double* origins_xyz, const double* directions_xyz,
    size_t ray_count, double t_max,
    alea_raycast_batch_selected_segment_callback_t callback, void* context) {
    if (!sys || !origins_xyz || !directions_xyz || !callback ||
        !isfinite(t_max) || t_max < 0.0)
        return -1;
    if (alea_system_prepare_query_caches(sys,
            ALEA_CACHE_HIER_SPATIAL | ALEA_CACHE_CELL_SURFACES) != 0)
        return -1;
    atomic_int failed;
    atomic_init(&failed, 0);
    raycast_visit_segments_batch_context_t parallel_context = {
        sys, origins_xyz, directions_xyz, t_max, callback, context, &failed
    };
    alea_parallel_status_t parallel_status = alea_parallel_for(
        ray_count, 1, 0, ALEA_PARALLEL_STATIC_BLOCK,
        raycast_visit_segments_batch_range, &parallel_context, NULL);
    if (parallel_status != ALEA_PARALLEL_OK) return -1;
    return atomic_load(&failed) || alea_interrupted() ? -1 : 0;
}
