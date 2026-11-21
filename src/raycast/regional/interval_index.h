// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
// SPDX-License-Identifier: MPL-2.0
#ifndef ALEA_INTERVAL_INDEX_INTERNAL_H
#define ALEA_INTERVAL_INDEX_INTERNAL_H
#include "raycast/raycast.h"
#include "core/alea_spatial_hier.h"

typedef alea_prepared_queries_t interval_index_t;
interval_index_t* interval_index_build(alea_system_t* sys, const alea_bbox_t* domain, bool selective);
typedef enum {
    CHILD_POINTS_BOUNDS,
    CHILD_POINTS_REGIONS,
    CHILD_POINTS_ADAPTIVE
} interval_child_mode_t;

/* IDs are in increasing order. In adaptive mode the listed children get full
 * regions at preparation and other children accumulate bounds on demand. */
typedef struct {
    const int* ids;
    size_t count;
    interval_child_mode_t mode;
} interval_point_children_t;
int alea_prepare_automatic_points(alea_system_t* sys);
int alea_prepare_root_points(alea_system_t* sys, const interval_point_children_t* children);
/* Benchmark only. May run alongside point readers on a prepared model.
 * Returns 0 when ready, 1 when another caller is building the child, and -1
 * on error. Join builders and readers before clearing, editing or destroying. */
int alea_promote_point_child(alea_system_t* sys, int universe_id);
uint64_t alea_root_point_scan_tests(const alea_system_t* sys, size_t universe,
                                    uint64_t* visits, double* local_seconds);
size_t alea_adaptive_child_progress(const alea_system_t* sys, size_t universe, double* seconds);
/* Try the preferred root owner in the covered region. A miss leaves the
 * usual surface-hinted and then unhinted hierarchy search to the walker. */
int alea_regional_ray_root_candidate(const interval_index_t* idx, double x, double y,
                                     double z, int crossed_surface, size_t* candidates);
bool interval_index_is_regional(const interval_index_t* index);
typedef struct { size_t visits,candidates,intersections; } interval_fallback_stats_t;
interval_fallback_stats_t interval_index_fallback_stats(const interval_index_t* index);
size_t interval_index_bytes(const interval_index_t* index);
typedef alea_query_point_t interval_point_t;
typedef alea_point_owner_t interval_owner_t;
typedef alea_segment_ray_t interval_line_t;
int alea_hierarchy_point_owner(alea_system_t* sys, double x, double y, double z,
                               int* cell, int* material, int* index);
int alea_regional_point_owner(const interval_index_t* idx, double x, double y, double z,
                              int* cell, int* material, int* index);
int alea_regional_prepare_render_bounds(interval_index_t* idx);
bool alea_regional_render_ready(const interval_index_t* idx, const alea_system_t* sys);
typedef struct regional_worker regional_worker_t;
regional_worker_t* alea_regional_worker_create(const interval_index_t* idx);
void alea_regional_worker_destroy(regional_worker_t* worker);
int alea_regional_trace_segments(regional_worker_t* worker, const alea_ray_t* ray, double limit,
                                 alea_raycast_result_t* result);
/* Experimental CPU transport interval. The caller supplies path storage. */
typedef struct {
    alea_hier_ray_path_t* path;
    int cell, surface;
    uint32_t primitive;
    double end;
    alea_matrix_t exit_frame;
    uint8_t resolution_flags;
} regional_nav_span_t;
int alea_regional_first_span(regional_worker_t* worker, const alea_ray_t* ray,
                            double start, double limit, alea_raycast_result_t* scratch,
                            regional_nav_span_t* span);
bool alea_regional_nav_ready(const interval_index_t* idx, const alea_system_t* sys);
void alea_regional_index_retain(interval_index_t* idx);
size_t alea_regional_worker_candidate_tree_visits(const regional_worker_t* worker);
/* Private replay selection. The navigator retains the prepared index. */
int alea_ray_navigator_use_regional(alea_ray_navigator_t* nav, interval_index_t* idx);
void alea_ray_navigator_regional_counts(const alea_ray_navigator_t* nav,
                                       size_t* spans, size_t* candidate_tree_visits);
/* Reserve once before dispatch. Each TinyPar worker exclusively uses its slot. */
int alea_regional_reserve_workers(interval_index_t* idx, size_t count);
regional_worker_t* alea_regional_worker(interval_index_t* idx, size_t slot);
int alea_regional_render_visible(regional_worker_t* worker, const alea_ray_t* ray, double limit,
                                int normal, alea_raycast_render_accept_t accept, void* context,
                                alea_raycast_result_t* scratch, alea_ray_first_visible_result_t* hit);
#endif
