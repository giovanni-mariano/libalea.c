// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
// SPDX-License-Identifier: MPL-2.0
#ifndef ALEA_INTERVAL_INDEX_INTERNAL_H
#define ALEA_INTERVAL_INDEX_INTERNAL_H
#include "raycast/raycast.h"

typedef alea_prepared_queries_t interval_index_t;
interval_index_t* interval_index_build(alea_system_t* sys, const alea_bbox_t* domain, bool selective);
bool interval_index_is_regional(const interval_index_t* index);
typedef struct { size_t visits,candidates,intersections; } interval_fallback_stats_t;
interval_fallback_stats_t interval_index_fallback_stats(const interval_index_t* index);
size_t interval_index_bytes(const interval_index_t* index);
typedef alea_query_point_t interval_point_t;
typedef alea_point_owner_t interval_owner_t;
typedef alea_segment_ray_t interval_line_t;
bool alea_regional_render_ready(const interval_index_t* idx, const alea_system_t* sys);
typedef struct regional_worker regional_worker_t;
/* Reserve once before dispatch. Each TinyPar worker exclusively uses its slot. */
int alea_regional_reserve_workers(interval_index_t* idx, size_t count);
regional_worker_t* alea_regional_worker(interval_index_t* idx, size_t slot);
int alea_regional_render_visible(regional_worker_t* worker, const alea_ray_t* ray, double limit,
                                int normal, alea_raycast_render_accept_t accept, void* context,
                                alea_raycast_result_t* scratch, alea_ray_first_visible_result_t* hit);
#endif
