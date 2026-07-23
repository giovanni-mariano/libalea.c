// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/**
 * @file slice_directional_trace.c
 * @brief Private directional boundary-event caches for slice consumers
 */

#include "alea.h"
#include "alea_slice.h"
#include "raycast/raycast.h"
#include "core/alea_system.h"
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

typedef struct {
    size_t line_count;
    size_t* offsets;
    alea_ray_boundary_event_t* events;
    size_t event_count;
} slice_directional_event_stream_t;

typedef struct {
    uint32_t event_fields;
    int projected_depth;
    uint8_t retain_occurrence_keys;
    uint64_t max_events;
    uint64_t max_output_bytes;
    uint8_t complete;
} slice_directional_contract_t;

struct alea_slice_directional_event_cache {
    const alea_system_t* sys;
    uint64_t geometry_generation;
    alea_slice_view_t view;
    int width, height;
    slice_directional_event_stream_t streams[2][2];
    slice_directional_contract_t contracts[2][2];
    alea_raycast_batch_result_t* ownership[2][2];
};

static void slice_world_point(const alea_slice_view_t* view, double u, double v,
                              double out[3]) {
    for (int i = 0; i < 3; i++)
        out[i] = view->plane.origin[i] + u * view->plane.u_axis[i] +
                 v * view->plane.v_axis[i];
}

static void slice_directional_event_stream_free(
    slice_directional_event_stream_t* cache) {
    if (!cache) return;
    free(cache->offsets);
    free(cache->events);
    memset(cache, 0, sizeof(*cache));
}

static int slice_directional_event_stream_build(
    alea_system_t* sys, const alea_slice_view_t* view, int width, int height,
    int orient, int reverse, slice_directional_event_stream_t* cache) {
    const int line_count = orient == ALEA_SLICE_EDGE_RIGHT ? height : width;
    const int samples = orient == ALEA_SLICE_EDGE_RIGHT ? width : height;
    const double step = orient == ALEA_SLICE_EDGE_RIGHT
        ? (view->u_max - view->u_min) / width
        : (view->v_max - view->v_min) / height;
    memset(cache, 0, sizeof(*cache));
    cache->line_count = (size_t)line_count;
    cache->offsets = calloc((size_t)line_count + 1, sizeof(*cache->offsets));
    if (!cache->offsets) return -1;
    if (samples < 2) return 0;

    const alea_ray_boundary_event_options_internal_t options = {
        .include_all_coincident_physical = true
    };
    const double length = (samples - 1) * step;
    alea_ray_boundary_event_t** line_events = calloc(
        (size_t)line_count, sizeof(*line_events));
    size_t* line_counts = calloc((size_t)line_count, sizeof(*line_counts));
    int failed = !line_events || !line_counts;
    if (failed) goto cleanup;

    #pragma omp parallel for schedule(static)
    for (int line = 0; line < line_count; line++) {
        alea_raycast_result_t trace;
        alea_ray_boundary_event_result_t events;
        alea_raycast_result_init(&trace);
        alea_ray_boundary_event_result_init(&events);
        double u0, v0, u1, v1;
        if (orient == ALEA_SLICE_EDGE_RIGHT) {
            u0 = reverse ? view->u_max - 0.5 * step : view->u_min + 0.5 * step;
            u1 = reverse ? view->u_min + 0.5 * step : view->u_max - 0.5 * step;
            v0 = v1 = view->v_min + (line + 0.5) *
                (view->v_max - view->v_min) / height;
        } else {
            v0 = reverse ? view->v_min + 0.5 * step : view->v_max - 0.5 * step;
            v1 = reverse ? view->v_max - 0.5 * step : view->v_min + 0.5 * step;
            u0 = u1 = view->u_min + (line + 0.5) *
                (view->u_max - view->u_min) / width;
        }
        double start[3], end[3];
        slice_world_point(view, u0, v0, start);
        slice_world_point(view, u1, v1, end);
        alea_ray_t ray;
        alea_ray_init_normalized(&ray, start[0], start[1], start[2],
                                 (end[0] - start[0]) / length,
                                 (end[1] - start[1]) / length,
                                 (end[2] - start[2]) / length);
        if (alea_raycast_boundary_events_with_options(
                sys, &ray, length, &options, &trace, &events) != 0) {
            #pragma omp atomic write
            failed = 1;
        } else if (events.events.count != 0) {
            line_events[line] = malloc(events.events.count * sizeof(*line_events[line]));
            if (!line_events[line]) {
                #pragma omp atomic write
                failed = 1;
            } else {
                memcpy(line_events[line], events.events.data,
                       events.events.count * sizeof(*line_events[line]));
                line_counts[line] = events.events.count;
            }
        }
        alea_ray_boundary_event_result_free(&events);
        alea_raycast_result_free(&trace);
    }
    if (failed) goto cleanup;
    for (int line = 0; line < line_count; line++) {
        if (line_counts[line] > SIZE_MAX - cache->event_count) {
            failed = 1;
            goto cleanup;
        }
        cache->offsets[line] = cache->event_count;
        cache->event_count += line_counts[line];
    }
    cache->offsets[line_count] = cache->event_count;
    if (cache->event_count != 0) {
        cache->events = malloc(cache->event_count * sizeof(*cache->events));
        if (!cache->events) { failed = 1; goto cleanup; }
        #pragma omp parallel for schedule(static)
        for (int line = 0; line < line_count; line++) {
            if (line_counts[line] != 0)
                memcpy(cache->events + cache->offsets[line], line_events[line],
                       line_counts[line] * sizeof(*cache->events));
        }
    }

cleanup:
    for (int line = 0; line < line_count; line++) free(line_events ? line_events[line] : NULL);
    free(line_events);
    free(line_counts);
    if (failed) slice_directional_event_stream_free(cache);
    return failed ? -1 : 0;
}

static int slice_directional_ownership_build(
    alea_system_t* sys, const alea_slice_view_t* view, int width, int height,
    int orient, int reverse, alea_raycast_batch_result_t** out) {
    const size_t count = (size_t)(orient == ALEA_SLICE_EDGE_RIGHT ? height : width);
    const double span = orient == ALEA_SLICE_EDGE_RIGHT ?
        view->u_max - view->u_min : view->v_max - view->v_min;
    double *origins = NULL, *directions = NULL;
    alea_raycast_batch_options_t options = {
        .struct_size = sizeof(options),
        .fields = ALEA_RAY_BATCH_PROJECTED_OWNER,
        .projected_depth = -1
    };
    alea_raycast_batch_result_t* result = NULL;
    if (count > SIZE_MAX / (3 * sizeof(*origins))) return -1;
    origins = malloc(count * 3 * sizeof(*origins));
    directions = malloc(count * 3 * sizeof(*directions));
    result = alea_raycast_batch_result_create();
    if (!origins || !directions || !result) goto fail;
    for (size_t line = 0; line < count; line++) {
        double u, v;
        if (orient == ALEA_SLICE_EDGE_RIGHT) {
            u = reverse ? view->u_max : view->u_min;
            v = view->v_min + ((double)line + 0.5) *
                (view->v_max - view->v_min) / (double)height;
        } else {
            u = view->u_min + ((double)line + 0.5) *
                (view->u_max - view->u_min) / (double)width;
            v = reverse ? view->v_min : view->v_max;
        }
        double p[3];
        slice_world_point(view, u, v, p);
        for (int axis = 0; axis < 3; axis++) {
            origins[line * 3 + axis] = p[axis];
            directions[line * 3 + axis] = orient == ALEA_SLICE_EDGE_RIGHT ?
                (reverse ? -view->plane.u_axis[axis] : view->plane.u_axis[axis]) :
                (reverse ? view->plane.v_axis[axis] : -view->plane.v_axis[axis]);
        }
    }
    if (alea_raycast_hier_batch(sys, origins, directions, count, span,
                                &options, result) != 0) goto fail;
    /* U-forward traces use slice coordinates, matching alea_trace_ray_slice_compact.
     * Reverse traces deliberately retain distance from U-max for validation. */
    if (orient == ALEA_SLICE_EDGE_RIGHT && !reverse) {
        double* enter = (double*)alea_raycast_batch_t_enter(result);
        double* exit = (double*)alea_raycast_batch_t_exit(result);
        for (size_t i = 0; i < alea_raycast_batch_segment_count(result); i++) {
            enter[i] += view->u_min;
            exit[i] += view->u_min;
        }
    }
    free(origins); free(directions);
    *out = result;
    return 0;
fail:
    free(origins); free(directions);
    alea_raycast_batch_result_destroy(result);
    return -1;
}

alea_slice_directional_event_cache_t* alea_slice_directional_event_cache_create(
    alea_system_t* sys, const alea_slice_view_t* view, int width, int height) {
    if (!sys || !view || width <= 0 || height <= 0) return NULL;
    alea_slice_directional_event_cache_t* cache = calloc(1, sizeof(*cache));
    if (!cache) return NULL;
    cache->sys = sys;
    cache->geometry_generation = alea_system_geometry_generation(sys);
    cache->view = *view;
    cache->width = width;
    cache->height = height;
    for (int orient = ALEA_SLICE_EDGE_RIGHT; orient <= ALEA_SLICE_EDGE_DOWN; orient++)
        for (int reverse = 0; reverse <= 1; reverse++)
            if (slice_directional_event_stream_build(sys, view, width, height,
                                                     orient, reverse,
                                                     &cache->streams[orient][reverse]) != 0) {
                alea_slice_directional_event_cache_destroy(cache);
                return NULL;
            } else if (slice_directional_ownership_build(sys, view, width, height,
                                                          orient, reverse,
                                                          &cache->ownership[orient][reverse]) != 0) {
                alea_slice_directional_event_cache_destroy(cache);
                return NULL;
            } else {
                cache->contracts[orient][reverse] = (slice_directional_contract_t){
                    .event_fields = ALEA_RAY_BOUNDARY_EVENT_PRIMITIVE_ID |
                                    ALEA_RAY_BOUNDARY_EVENT_NORMAL,
                    .projected_depth = -1,
                    .retain_occurrence_keys = 1,
                    .complete = 1
                };
            }
    return cache;
}

alea_slice_directional_trace_cache_t*
alea_slice_directional_trace_cache_create(
    alea_system_t* sys, const alea_slice_view_t* view, int width, int height) {
    return alea_slice_directional_event_cache_create(sys, view, width, height);
}

void alea_slice_directional_event_cache_destroy(
    alea_slice_directional_event_cache_t* cache) {
    if (!cache) return;
    for (int orient = ALEA_SLICE_EDGE_RIGHT; orient <= ALEA_SLICE_EDGE_DOWN; orient++)
        for (int reverse = 0; reverse <= 1; reverse++) {
            slice_directional_event_stream_free(&cache->streams[orient][reverse]);
            alea_raycast_batch_result_destroy(cache->ownership[orient][reverse]);
        }
    free(cache);
}

void alea_slice_directional_trace_cache_destroy(
    alea_slice_directional_trace_cache_t* cache) {
    alea_slice_directional_event_cache_destroy(cache);
}

int alea_slice_directional_event_cache_contract(
    const alea_slice_directional_event_cache_t* cache,
    alea_slice_edge_orientation_t orientation, int reverse,
    uint32_t required_event_fields, int projected_depth,
    uint64_t max_events, uint64_t max_output_bytes, int* out_complete) {
    if (!cache || !out_complete || (orientation != ALEA_SLICE_EDGE_RIGHT &&
        orientation != ALEA_SLICE_EDGE_DOWN) || (reverse != 0 && reverse != 1)) return -1;
    const slice_directional_contract_t* contract = &cache->contracts[orientation][reverse];
    const slice_directional_event_stream_t* stream =
        &cache->streams[orientation][reverse];
    *out_complete = contract->complete &&
        (contract->event_fields & required_event_fields) == required_event_fields &&
        contract->projected_depth == projected_depth &&
        (!max_events || stream->event_count <= max_events) &&
        (!max_output_bytes ||
         stream->event_count <= max_output_bytes / sizeof(*stream->events));
    return 0;
}

const alea_raycast_batch_result_t*
alea_slice_directional_event_cache_ownership_trace(
    const alea_slice_directional_event_cache_t* cache,
    alea_slice_edge_orientation_t orientation, int reverse,
    int projected_depth, uint32_t required_fields) {
    int complete = 0;
    if (alea_slice_directional_event_cache_contract(
            cache, orientation, reverse, 0, projected_depth, 0, 0, &complete) != 0 ||
        !complete) return NULL;
    /* Ownership traces intentionally materialize only the requested owner
     * projection; richer surface/path requests must use another producer. */
    if (required_fields & ~ALEA_RAY_BATCH_PROJECTED_OWNER) return NULL;
    return cache->ownership[orientation][reverse];
}

int alea_slice_directional_event_cache_matches(
    const alea_slice_directional_event_cache_t* cache,
    const alea_system_t* sys, const alea_slice_view_t* view,
    int width, int height) {
    if (!cache || !sys || !view || cache->sys != sys ||
        cache->geometry_generation != alea_system_geometry_generation(sys) ||
        cache->width != width || cache->height != height)
        return 0;
    for (int axis = 0; axis < 3; axis++) {
        if (cache->view.plane.origin[axis] != view->plane.origin[axis] ||
            cache->view.plane.normal[axis] != view->plane.normal[axis] ||
            cache->view.plane.u_axis[axis] != view->plane.u_axis[axis] ||
            cache->view.plane.v_axis[axis] != view->plane.v_axis[axis])
            return 0;
    }
    return cache->view.u_min == view->u_min && cache->view.u_max == view->u_max &&
           cache->view.v_min == view->v_min && cache->view.v_max == view->v_max;
}

int alea_slice_directional_trace_cache_matches(
    const alea_slice_directional_trace_cache_t* cache,
    const alea_system_t* sys, const alea_slice_view_t* view,
    int width, int height) {
    return alea_slice_directional_event_cache_matches(cache, sys, view, width, height);
}

int alea_slice_directional_event_cache_line_events(
    const alea_slice_directional_event_cache_t* cache,
    alea_slice_edge_orientation_t orientation, int reverse, size_t line,
    const alea_ray_boundary_event_t** out_events, size_t* out_count) {
    if (!cache || !out_events || !out_count ||
        (orientation != ALEA_SLICE_EDGE_RIGHT &&
         orientation != ALEA_SLICE_EDGE_DOWN) ||
        (reverse != 0 && reverse != 1))
        return -1;
    const slice_directional_event_stream_t* stream =
        &cache->streams[orientation][reverse];
    if (line >= stream->line_count) return -1;
    *out_count = stream->offsets[line + 1] - stream->offsets[line];
    *out_events = *out_count ? stream->events + stream->offsets[line] : NULL;
    return 0;
}

int alea_slice_directional_event_cache_dimensions(
    const alea_slice_directional_event_cache_t* cache,
    int* out_width, int* out_height) {
    if (!cache || !out_width || !out_height) return -1;
    *out_width = cache->width;
    *out_height = cache->height;
    return 0;
}
