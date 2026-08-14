// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/**
 * @file ray_request.c
 * @brief Validation and lowering of internal semantic ray requests.
 */

#include "raycast.h"

#include <float.h>
#include <math.h>
#include <string.h>

int alea_ray_query_lower(const alea_ray_query_t* query,
                         alea_ray_plan_t* out_plan) {
    if (!query || !out_plan) return -1;
    if (query->kind < ALEA_RAY_QUERY_ANY_HIT ||
        query->kind > ALEA_RAY_QUERY_BOUNDARY_EVENTS ||
        query->backend < ALEA_RAY_QUERY_BACKEND_AUTO ||
        query->backend > ALEA_RAY_QUERY_BACKEND_FAST_FORWARD_REVERSE ||
        !isfinite(query->t_min) || !isfinite(query->t_max) ||
        query->t_min < 0 ||
        (query->t_max > 0 && query->t_min > query->t_max))
        return -1;

    alea_ray_plan_t plan = {
        .product = query->kind,
        .backend = query->backend,
        .engine = ALEA_RAY_ENGINE_GLOBAL_BREAKPOINTS,
        .ownership = ALEA_RAY_OWNERSHIP_SELECT_CANONICAL,
        .t_min = query->t_min,
        .t_max = query->t_max <= 0.0 ? DBL_MAX : query->t_max,
        .material_filter = query->material_filter,
        .max_events = query->max_events,
        .max_output_bytes = query->max_output_bytes
    };

    plan.requirements.need_density =
        (query->fields & ALEA_RAY_QUERY_FIELD_DENSITY) != 0;
    plan.requirements.need_surface_identity =
        (query->fields & (ALEA_RAY_QUERY_FIELD_SURFACE_ID |
                          ALEA_RAY_QUERY_FIELD_PRIMITIVE_ID)) != 0;
    plan.requirements.need_normal =
        (query->fields & ALEA_RAY_QUERY_FIELD_SURFACE_NORMAL) != 0;
    if (plan.requirements.need_normal)
        plan.requirements.need_surface_identity = 1;

    if (query->kind == ALEA_RAY_QUERY_BOUNDARY_EVENTS) {
        plan.engine = ALEA_RAY_ENGINE_GLOBAL_BREAKPOINTS;
        plan.ownership = ALEA_RAY_OWNERSHIP_SELECT_CANONICAL;
        plan.requirements.need_selected_owner = 1;
        plan.requirements.need_all_coincident_primitives = 1;
        plan.requirements.need_surface_identity = 1;
    } else {
        plan.requirements.need_selected_owner = 1;
    }

    if (query->backend == ALEA_RAY_QUERY_BACKEND_FAST_FORWARD ||
        query->backend == ALEA_RAY_QUERY_BACKEND_FAST_REVERSE ||
        query->backend == ALEA_RAY_QUERY_BACKEND_FAST_FORWARD_REVERSE ||
        (query->backend == ALEA_RAY_QUERY_BACKEND_AUTO &&
         (query->kind == ALEA_RAY_QUERY_FIRST_CELL ||
          query->kind == ALEA_RAY_QUERY_FIRST_VISIBLE ||
          query->kind == ALEA_RAY_QUERY_ANY_HIT ||
          query->kind == ALEA_RAY_QUERY_SEGMENTS))) {
        plan.engine = ALEA_RAY_ENGINE_SELECTED_WALKER;
        plan.ownership = ALEA_RAY_OWNERSHIP_TRACK_COHERENT;
    }

    *out_plan = plan;
    return 0;
}
