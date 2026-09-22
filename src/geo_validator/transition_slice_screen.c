// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

#include "alea_geo_validator.h"
#include "transition_slice_critical.h"
#include "transition_validation.h"

#include "raycast/raycast.h"
#include "raycast/ray_epsilon.h"
#include "core/alea_system.h"
#include "primitives/bbox.h"
#include "util/alea_parallel.h"

#include <float.h>
#include <math.h>
#include <stdatomic.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

typedef struct {
    size_t cell_index;
    alea_bbox_t box;
} slice_error_index_cell_t;

struct alea_slice_error_query {
    alea_system_t* sys;
    alea_slice_error_query_options_t options;
    uint64_t query_id;
    uint64_t geometry_generation;
    size_t page_count;
    slice_error_index_cell_t* indexed_cells;
    size_t indexed_count;
    size_t* uncertain_cells;
    size_t uncertain_count;
    size_t uncached_cell_begin;
    size_t index_bytes;
};

struct alea_slice_error_page {
    alea_slice_error_page_receipt_t receipt;
    alea_transition_slice_critical_finding_t* findings;
    size_t finding_count;
    size_t finding_capacity;
    alea_slice_error_interval_t* intervals;
    size_t interval_count;
    alea_slice_error_region_t* regions;
    size_t region_count;
    int populated;
};

static atomic_uint_fast64_t slice_error_next_query_id = 1;
static int slice_error_build_index(alea_slice_error_query_t* query);
static int slice_error_slice_axis_for_world(const alea_slice_view_t* view,
                                           int world_axis);
static double slice_error_slice_axis_sign(const alea_slice_view_t* view,
                                          int world_axis);

void alea_slice_error_query_options_init(
    alea_slice_error_query_options_t* options) {
    if (!options) return;
    memset(options, 0, sizeof(*options));
    options->struct_size = sizeof(*options);
    options->tile_columns = 1;
    options->tile_rows = 1;
    options->max_index_bytes = 32u * 1024u * 1024u;
    alea_transition_slice_options_init(&options->scan_options);
    options->scan_options.occurrence_discovery =
        ALEA_TRANSITION_SLICE_OCCURRENCE_EXHAUSTIVE;
}

static int slice_error_finite_view(const alea_slice_view_t* view) {
    for (size_t i = 0; i < 3; ++i) {
        if (!isfinite(view->plane.origin[i]) ||
            !isfinite(view->plane.normal[i]) ||
            !isfinite(view->plane.u_axis[i]) ||
            !isfinite(view->plane.v_axis[i])) return 0;
    }
    return isfinite(view->u_min) && isfinite(view->u_max) &&
        isfinite(view->v_min) && isfinite(view->v_max) &&
        view->u_min < view->u_max && view->v_min < view->v_max;
}

alea_slice_error_query_t* alea_slice_error_query_create(
    alea_system_t* sys, const alea_slice_error_query_options_t* input) {
    if (!sys || !input || input->struct_size <
            offsetof(alea_slice_error_query_options_t, scan_options))
        return NULL;
    alea_slice_error_query_options_t options;
    alea_slice_error_query_options_init(&options);
    const size_t copy_size = input->struct_size < sizeof(options)
        ? input->struct_size : sizeof(options);
    memcpy(&options, input, copy_size);
    if (!options.max_index_bytes)
        options.max_index_bytes = 32u * 1024u * 1024u;
    const alea_slice_view_t* view = &options.view;
    if (!slice_error_finite_view(view) ||
        !isfinite(options.required_uv_min[0]) ||
        !isfinite(options.required_uv_min[1]) ||
        !isfinite(options.required_uv_max[0]) ||
        !isfinite(options.required_uv_max[1]) ||
        !(options.required_uv_max[0] > options.required_uv_min[0]) ||
        !(options.required_uv_max[1] > options.required_uv_min[1]) ||
        !isfinite(options.required_uv_max[0] -
                  options.required_uv_min[0]) ||
        !isfinite(options.required_uv_max[1] -
                  options.required_uv_min[1]) ||
        !options.tile_columns || !options.tile_rows ||
        options.tile_columns > SIZE_MAX / options.tile_rows)
        return NULL;
    /* A candidate scan may use the existing critical machinery, but its
     * sampled discovery mode cannot establish exhaustive occurrence coverage. */
    options.scan_options.occurrence_discovery =
        ALEA_TRANSITION_SLICE_OCCURRENCE_EXHAUSTIVE;
    if (!options.scan_options.max_curves_per_tile ||
        !options.scan_options.max_critical_points ||
        !options.scan_options.max_coverage_hits ||
        !options.scan_options.max_exhaustive_occurrence_hits ||
        !options.scan_options.max_critical_scratch_bytes ||
        !options.scan_options.max_critical_findings ||
        !options.scan_options.max_output_bytes)
        return NULL;
    alea_slice_error_query_t* query = calloc(1, sizeof(*query));
    if (!query) return NULL;
    query->sys = sys;
    query->options = options;
    query->query_id = atomic_fetch_add(&slice_error_next_query_id, 1);
    query->geometry_generation = alea_system_geometry_generation(sys);
    query->page_count = options.tile_columns * options.tile_rows;
    if (slice_error_build_index(query) != 0) {
        alea_slice_error_query_destroy(query);
        return NULL;
    }
    return query;
}

void alea_slice_error_query_destroy(alea_slice_error_query_t* query) {
    if (!query) return;
    free(query->indexed_cells);
    free(query->uncertain_cells);
    free(query);
}

size_t alea_slice_error_query_page_count(const alea_slice_error_query_t* query) {
    return query ? query->page_count : 0;
}

uint64_t alea_slice_error_query_id(const alea_slice_error_query_t* query) {
    return query ? query->query_id : 0;
}

alea_slice_error_page_t* alea_slice_error_page_create(void) {
    return calloc(1, sizeof(alea_slice_error_page_t));
}

void alea_slice_error_page_destroy(alea_slice_error_page_t* page) {
    if (!page) return;
    free(page->findings);
    free(page->intervals);
    free(page->regions);
    free(page);
}

int alea_slice_error_page_receipt(const alea_slice_error_page_t* page,
                                  alea_slice_error_page_receipt_t* out_receipt) {
    if (!page || !out_receipt || !page->populated) return -1;
    *out_receipt = page->receipt;
    return 0;
}

size_t alea_slice_error_page_context_finding_count(
    const alea_slice_error_page_t* page) {
    return page && page->populated ? page->finding_count : 0;
}

int alea_slice_error_page_context_finding_get(
    const alea_slice_error_page_t* page, size_t index,
    alea_transition_slice_critical_finding_t* out_finding) {
    if (!page || !page->populated || !out_finding ||
        index >= page->finding_count) return -1;
    *out_finding = page->findings[index];
    return 0;
}

size_t alea_slice_error_page_interval_count(const alea_slice_error_page_t* page) {
    return page && page->populated ? page->interval_count : 0;
}

int alea_slice_error_page_interval_get(const alea_slice_error_page_t* page,
                                       size_t index,
                                       alea_slice_error_interval_t* out_interval) {
    if (!page || !page->populated || !out_interval ||
        index >= page->interval_count) return -1;
    *out_interval = page->intervals[index];
    return 0;
}

size_t alea_slice_error_page_region_count(const alea_slice_error_page_t* page) {
    return page && page->populated ? page->region_count : 0;
}

int alea_slice_error_page_region_get(const alea_slice_error_page_t* page,
                                     size_t index,
                                     alea_slice_error_region_t* out_region) {
    if (!page || !page->populated || !out_region ||
        index >= page->region_count) return -1;
    *out_region = page->regions[index];
    return 0;
}

typedef struct {
    alea_slice_error_page_t* page;
    size_t max_findings;
    size_t max_output_bytes;
    size_t omitted_findings;
} slice_error_finding_sink_t;

static int slice_error_retain_finding(
    const alea_transition_slice_critical_finding_t* finding, void* userdata) {
    slice_error_finding_sink_t* sink = userdata;
    alea_slice_error_page_t* page = sink->page;
    if (page->finding_count >= sink->max_findings ||
        page->finding_count >=
            sink->max_output_bytes / sizeof(*page->findings)) {
        if (sink->omitted_findings < SIZE_MAX) sink->omitted_findings++;
        return 1;
    }
    if (page->finding_count == page->finding_capacity) {
        size_t capacity = page->finding_capacity
            ? (page->finding_capacity > SIZE_MAX / 2u
                ? SIZE_MAX : page->finding_capacity * 2u) : 4u;
        if (capacity > sink->max_findings) capacity = sink->max_findings;
        const size_t byte_capacity =
            sink->max_output_bytes / sizeof(*page->findings);
        if (capacity > byte_capacity) capacity = byte_capacity;
        void* next = realloc(page->findings,
                             capacity * sizeof(*page->findings));
        if (!next) return -1;
        page->findings = next;
        page->finding_capacity = capacity;
    }
    page->findings[page->finding_count++] = *finding;
    return 0;
}

typedef struct {
    int axis;
    int coefficient_sign;
    double coordinate;
    double uncertainty;
    int surface_id;
    uint32_t primitive_id;
} slice_error_axis_plane_t;

/* The analytic CSG bbox is conservative for these operations. Reject
 * complement and malformed trees before using it to exclude a cell. */
static int slice_error_bbox_tree_safe(const alea_system_t* sys,
                                     alea_node_id_t id, size_t depth,
                                     size_t* work) {
    if (!*work || depth > 128 || id >= sys->nodes.count) return 0;
    --*work;
    const alea_node_t* node = &sys->nodes.data[id];
    const alea_operation_t op = ALEA_GET_OPERATION(node);
    if (op == ALEA_OP_PRIMITIVE) {
        const uint32_t id = node->primitive.primitive_id;
        if (id >= sys->primitives.count) return 0;
        const alea_primitive_entry_t* primitive = &sys->primitives.data[id];
        if (node->primitive.prim_type != primitive->type) return 0;
        if (primitive->type == ALEA_PRIMITIVE_SPHERE ||
            primitive->type == ALEA_PRIMITIVE_SPH) return 1;
        if (primitive->type != ALEA_PRIMITIVE_PLANE ||
            primitive->payload_index >= sys->primitive_planes.count)
            return 0;
        const alea_plane_data_t* p =
            &sys->primitive_planes.data[primitive->payload_index];
        if (!isfinite(p->a) || !isfinite(p->b) ||
            !isfinite(p->c) || !isfinite(p->d)) return 0;
        return (p->a != 0.0 && p->b == 0.0 && p->c == 0.0) ||
               (p->b != 0.0 && p->a == 0.0 && p->c == 0.0) ||
               (p->c != 0.0 && p->a == 0.0 && p->b == 0.0);
    }
    if (op != ALEA_OP_UNION && op != ALEA_OP_INTERSECTION &&
        op != ALEA_OP_DIFFERENCE) return 0;
    return slice_error_bbox_tree_safe(sys, node->operation.left,
                                     depth + 1, work) &&
           slice_error_bbox_tree_safe(sys, node->operation.right,
                                     depth + 1, work);
}

static int slice_error_index_cell_compare(const void* lhs, const void* rhs) {
    const slice_error_index_cell_t* a = lhs;
    const slice_error_index_cell_t* b = rhs;
    if (a->box.min_x < b->box.min_x) return -1;
    if (a->box.min_x > b->box.min_x) return 1;
    return (a->cell_index > b->cell_index) -
           (a->cell_index < b->cell_index);
}

static int slice_error_size_compare(const void* lhs, const void* rhs) {
    const size_t a = *(const size_t*)lhs;
    const size_t b = *(const size_t*)rhs;
    return (a > b) - (a < b);
}

static double slice_error_trusted_lower(double value) {
    return isfinite(value) && value > -1e10 && value < 1e10
        ? value : -INFINITY;
}

static double slice_error_trusted_upper(double value) {
    return isfinite(value) && value > -1e10 && value < 1e10
        ? value : INFINITY;
}

/* Build once per model revision. A cell with uncertain analytic support is
 * always visited; it is never excluded by the sorted bounded-cell index. */
static int slice_error_build_index(alea_slice_error_query_t* query) {
    const alea_system_t* sys = query->sys;
    size_t roots = 0;
    for (size_t i = 0; i < sys->cells.count; ++i)
        roots += sys->cells.data[i].universe_id == 0;
    const size_t per_root = sizeof(*query->indexed_cells) +
        sizeof(*query->uncertain_cells);
    size_t capacity = query->options.max_index_bytes / per_root;
    if (capacity > roots) capacity = roots;
    query->indexed_cells = capacity
        ? malloc(capacity * sizeof(*query->indexed_cells)) : NULL;
    query->uncertain_cells = capacity
        ? malloc(capacity * sizeof(*query->uncertain_cells)) : NULL;
    if (capacity && (!query->indexed_cells || !query->uncertain_cells))
        return -1;
    query->index_bytes = capacity * per_root;
    query->uncached_cell_begin = sys->cells.count;
    for (size_t i = 0; i < sys->cells.count; ++i) {
        const alea_cell_entry_t* cell = &sys->cells.data[i];
        if (cell->universe_id != 0) continue;
        if (query->indexed_count + query->uncertain_count == capacity) {
            query->uncached_cell_begin = i;
            break;
        }
        size_t work = query->options.scan_options.max_active_boundary_tests;
        if (work > 4096) work = 4096;
        if (cell->original_root_node_id != ALEA_NODE_ID_INVALID ||
            !slice_error_bbox_tree_safe(sys, cell->root_node_id, 0,
                                        &work)) {
            query->uncertain_cells[query->uncertain_count++] = i;
            continue;
        }
        alea_bbox_t box = alea_get_bbox(sys, cell->root_node_id);
        if (!(box.min_x <= box.max_x && box.min_y <= box.max_y &&
              box.min_z <= box.max_z)) {
            query->uncertain_cells[query->uncertain_count++] = i;
            continue;
        }
        box.min_x = slice_error_trusted_lower(box.min_x);
        box.max_x = slice_error_trusted_upper(box.max_x);
        box.min_y = slice_error_trusted_lower(box.min_y);
        box.max_y = slice_error_trusted_upper(box.max_y);
        box.min_z = slice_error_trusted_lower(box.min_z);
        box.max_z = slice_error_trusted_upper(box.max_z);
        query->indexed_cells[query->indexed_count++] =
            (slice_error_index_cell_t){i, box};
    }
    if (query->indexed_count > 1)
        qsort(query->indexed_cells, query->indexed_count,
              sizeof(*query->indexed_cells), slice_error_index_cell_compare);
    return 0;
}

/* Finite analytic bounds exclude a root cell only when separated from the
 * whole tile slab by more than coordinate rounding. */
static int slice_error_box_outside_tile(
    const alea_bbox_t* box, const alea_slice_error_query_options_t* o,
    const alea_transition_slice_critical_tile_t* tile) {
    double lower[3], upper[3];
    for (int world_axis = 0; world_axis < 2; ++world_axis) {
        const int slice_axis = slice_error_slice_axis_for_world(
            &o->view, world_axis);
        const double sign = slice_error_slice_axis_sign(
            &o->view, world_axis);
        lower[world_axis] = o->view.plane.origin[world_axis] +
            (sign > 0.0 ? tile->uv_min[slice_axis]
                        : -tile->uv_max[slice_axis]);
        upper[world_axis] = o->view.plane.origin[world_axis] +
            (sign > 0.0 ? tile->uv_max[slice_axis]
                        : -tile->uv_min[slice_axis]);
    }
    lower[2] = upper[2] = o->view.plane.origin[2];
    const double lo[3] = {box->min_x, box->min_y, box->min_z};
    const double hi[3] = {box->max_x, box->max_y, box->max_z};
    for (int axis = 0; axis < 3; ++axis) {
        if (!isfinite(lower[axis]) || !isfinite(upper[axis])) return 0;
        const double margin = 64.0 * DBL_EPSILON *
            fmax(1.0, fmax(fabs(lower[axis]), fabs(upper[axis])));
        if (hi[axis] < lower[axis] - margin ||
            lo[axis] > upper[axis] + margin) return 1;
    }
    return 0;
}

static int slice_error_plane_index(const slice_error_axis_plane_t* planes,
                                   size_t count, uint32_t primitive_id) {
    for (size_t i = 0; i < count; ++i)
        if (planes[i].primitive_id == primitive_id) return (int)i;
    return -1;
}

static int slice_error_collect_primitives(const alea_system_t* sys,
                                          alea_node_id_t id, size_t depth,
                                          slice_error_axis_plane_t* planes,
                                          size_t* count, size_t max_count,
                                          size_t* work) {
    if (!*work || depth > 128 || id >= sys->nodes.count) return 0;
    --*work;
    const alea_node_t* node = &sys->nodes.data[id];
    const alea_operation_t op = ALEA_GET_OPERATION(node);
    if (op == ALEA_OP_PRIMITIVE) {
        const uint32_t primitive_id = node->primitive.primitive_id;
        if (primitive_id >= sys->primitives.count) return 0;
        if (slice_error_plane_index(planes, *count, primitive_id) >= 0)
            return 1;
        if (*count == max_count) return 0;
        planes[*count].primitive_id = primitive_id;
        ++*count;
        return 1;
    }
    if (op != ALEA_OP_UNION && op != ALEA_OP_INTERSECTION &&
        op != ALEA_OP_DIFFERENCE && op != ALEA_OP_COMPLEMENT) return 0;
    if (!slice_error_collect_primitives(sys, node->operation.left, depth + 1,
                                        planes, count, max_count, work)) return 0;
    return op == ALEA_OP_COMPLEMENT ||
        slice_error_collect_primitives(sys, node->operation.right, depth + 1,
                                       planes, count, max_count, work);
}

typedef struct {
    double coordinate;
    double uncertainty;
    int surface_id;
    uint32_t primitive_id;
} slice_error_grid_line_t;

typedef struct {
    alea_point_coverage_kind_t kind;
    size_t owner_count;
    int owner_cell_ids[ALEA_SLICE_ERROR_OWNER_CAPACITY];
} slice_error_face_t;

static int slice_error_grid_line_compare(const void* aa, const void* bb) {
    const slice_error_grid_line_t* a = aa;
    const slice_error_grid_line_t* b = bb;
    return (a->coordinate > b->coordinate) -
           (a->coordinate < b->coordinate);
}

static int slice_error_axis_view_supported(const alea_slice_view_t* view) {
    const alea_slice_plane_t* p = &view->plane;
    return p->normal[0] == 0.0 && p->normal[1] == 0.0 &&
        fabs(p->normal[2]) == 1.0 &&
        p->u_axis[2] == 0.0 && p->v_axis[2] == 0.0 &&
        ((fabs(p->u_axis[0]) == 1.0 && p->u_axis[1] == 0.0 &&
          p->v_axis[0] == 0.0 && fabs(p->v_axis[1]) == 1.0) ||
         (p->u_axis[0] == 0.0 && fabs(p->u_axis[1]) == 1.0 &&
          fabs(p->v_axis[0]) == 1.0 && p->v_axis[1] == 0.0));
}

static int slice_error_slice_axis_for_world(const alea_slice_view_t* view,
                                           int world_axis) {
    return view->plane.u_axis[world_axis] != 0.0 ? 0 : 1;
}

static double slice_error_slice_axis_sign(const alea_slice_view_t* view,
                                          int world_axis) {
    return view->plane.u_axis[world_axis] != 0.0
        ? view->plane.u_axis[world_axis]
        : view->plane.v_axis[world_axis];
}

/* This first proof path intentionally accepts only root cells formed from
 * vertical X/Y planes. In each open grid face every primitive sign is fixed,
 * so the Boolean cell ownership is invariant throughout that face. */
static int slice_error_axis_node_inside(
    const alea_system_t* sys, alea_node_id_t id,
    const slice_error_axis_plane_t* planes, size_t plane_count,
    const slice_error_grid_line_t* x, size_t xi,
    const slice_error_grid_line_t* y, size_t yi,
    size_t depth, size_t* work_remaining, int* out_inside) {
    if (*work_remaining == 0) return 0;
    --*work_remaining;
    if (id >= sys->nodes.count || depth > 128) return 0;
    const alea_node_t* node = &sys->nodes.data[id];
    const alea_operation_t op = ALEA_GET_OPERATION(node);
    if (op == ALEA_OP_PRIMITIVE) {
        const uint32_t primitive_id = node->primitive.primitive_id;
        if (primitive_id >= sys->primitives.count) return 0;
        const int plane_index = slice_error_plane_index(
            planes, plane_count, primitive_id);
        if (plane_index < 0) return 0;
        const slice_error_axis_plane_t* plane = &planes[plane_index];
        const slice_error_grid_line_t* lines = plane->axis == 0 ? x : y;
        const size_t index = plane->axis == 0 ? xi : yi;
        int side;
        if (lines[index + 1].coordinate <= plane->coordinate)
            side = -1;
        else if (lines[index].coordinate >= plane->coordinate)
            side = 1;
        else
            return 0;
        const int raw_negative = side * plane->coefficient_sign < 0;
        const int flip = (node->primitive.sense > 0) !=
                         (node->primitive.inverted != 0);
        *out_inside = raw_negative != flip;
        return 1;
    }
    if (op != ALEA_OP_UNION && op != ALEA_OP_INTERSECTION &&
        op != ALEA_OP_DIFFERENCE && op != ALEA_OP_COMPLEMENT) return 0;
    int left = 0, right = 0;
    if (!slice_error_axis_node_inside(sys, node->operation.left, planes,
                                      plane_count,
                                      x, xi, y, yi, depth + 1,
                                      work_remaining, &left)) return 0;
    if (op == ALEA_OP_COMPLEMENT) {
        *out_inside = !left;
        return 1;
    }
    if (!slice_error_axis_node_inside(sys, node->operation.right, planes,
                                      plane_count,
                                      x, xi, y, yi, depth + 1,
                                      work_remaining, &right)) return 0;
    *out_inside = op == ALEA_OP_UNION ? left || right
        : op == ALEA_OP_INTERSECTION ? left && right
        : left && !right;
    return 1;
}

static int slice_error_axis_grid_valid(
    const slice_error_grid_line_t* lines, size_t count) {
    for (size_t i = 1; i < count; ++i) {
        const double a = lines[i - 1].coordinate;
        const double b = lines[i].coordinate;
        const double margin = lines[i - 1].uncertainty +
            lines[i].uncertainty +
            4.0 * DBL_EPSILON * fmax(1.0, fmax(fabs(a), fabs(b)));
        if (!(b > a) || !(b - a > margin)) return 0;
    }
    return 1;
}

static int slice_error_axis_is_defect(alea_point_coverage_kind_t kind) {
    return kind == ALEA_POINT_COVERAGE_GAP ||
           kind == ALEA_POINT_COVERAGE_OVERLAP;
}

typedef struct {
    size_t count;
    double uv[ALEA_SLICE_ERROR_POLYGON_CAPACITY][2];
    double uncertainty[ALEA_SLICE_ERROR_POLYGON_CAPACITY][2];
} slice_error_polygon_t;

static int slice_error_oblique_node_inside(
    const alea_system_t* sys, alea_node_id_t id,
    const uint32_t* primitive_ids, const int* raw_negative,
    size_t plane_count, size_t depth, size_t* work, int* inside) {
    if (!*work || depth > 128 || id >= sys->nodes.count) return 0;
    --*work;
    const alea_node_t* node = &sys->nodes.data[id];
    const alea_operation_t op = ALEA_GET_OPERATION(node);
    if (op == ALEA_OP_PRIMITIVE) {
        size_t plane = 0;
        while (plane < plane_count &&
               primitive_ids[plane] != node->primitive.primitive_id)
            ++plane;
        if (plane == plane_count) return 0;
        const int flip = (node->primitive.sense > 0) !=
                         (node->primitive.inverted != 0);
        *inside = raw_negative[plane] != flip;
        return 1;
    }
    if (op != ALEA_OP_UNION && op != ALEA_OP_INTERSECTION &&
        op != ALEA_OP_DIFFERENCE && op != ALEA_OP_COMPLEMENT) return 0;
    int left = 0, right = 0;
    if (!slice_error_oblique_node_inside(
            sys, node->operation.left, primitive_ids, raw_negative,
            plane_count, depth + 1, work, &left)) return 0;
    if (op == ALEA_OP_COMPLEMENT) {
        *inside = !left;
        return 1;
    }
    if (!slice_error_oblique_node_inside(
            sys, node->operation.right, primitive_ids, raw_negative,
            plane_count, depth + 1, work, &right)) return 0;
    *inside = op == ALEA_OP_UNION ? left || right
        : op == ALEA_OP_INTERSECTION ? left && right
        : left && !right;
    return 1;
}

static int slice_error_polygon_add(slice_error_polygon_t* polygon,
                                   double u, double v,
                                   double uncertainty) {
    if (polygon->count == ALEA_SLICE_ERROR_POLYGON_CAPACITY ||
        !isfinite(u) || !isfinite(v) || !isfinite(uncertainty)) return 0;
    const size_t index = polygon->count++;
    polygon->uv[index][0] = u;
    polygon->uv[index][1] = v;
    polygon->uncertainty[index][0] = uncertainty;
    polygon->uncertainty[index][1] = uncertainty;
    return 1;
}

/* A single oblique line partitions a rectangular core into two convex faces.
 * This deliberately excludes line corners and multiple-line arrangements. */
static int slice_error_classify_single_oblique_tile(
    const alea_slice_error_query_t* query,
    const alea_transition_slice_critical_tile_t* tile,
    const size_t* cells, size_t cell_count,
    const slice_error_axis_plane_t* selected,
    size_t contextual_bytes, size_t* work,
    alea_slice_error_page_t* out,
    alea_slice_error_unresolved_reason_t* reason,
    int* output_omitted) {
    const alea_slice_error_query_options_t* o = &query->options;
    const alea_system_t* sys = query->sys;
    if (selected->surface_id <= 0) return 0;
    const alea_primitive_entry_t* primitive =
        &sys->primitives.data[selected->primitive_id];
    if (primitive->type != ALEA_PRIMITIVE_PLANE ||
        primitive->payload_index >= sys->primitive_planes.count)
        return 0;
    const alea_plane_data_t* p =
        &sys->primitive_planes.data[primitive->payload_index];
    if (!isfinite(p->a) || !isfinite(p->b) ||
        !isfinite(p->c) || !isfinite(p->d) ||
        p->a == 0.0 || p->b == 0.0 || p->c != 0.0) return 0;
    const alea_slice_plane_t* frame = &o->view.plane;
    long double a = (long double)p->a * frame->u_axis[0] +
                    (long double)p->b * frame->u_axis[1];
    long double b = (long double)p->a * frame->v_axis[0] +
                    (long double)p->b * frame->v_axis[1];
    long double c = (long double)p->a * frame->origin[0] +
                    (long double)p->b * frame->origin[1] +
                    (long double)p->d;
    const long double norm = fmaxl(fabsl(a), fabsl(b));
    if (!(norm > 0.0L) || !isfinite(norm)) {
        *reason = ALEA_SLICE_ERROR_UNRESOLVED_NUMERICAL;
        return 0;
    }
    a /= norm; b /= norm; c /= norm;
    const double corner[4][2] = {
        {tile->uv_min[0], tile->uv_min[1]},
        {tile->uv_max[0], tile->uv_min[1]},
        {tile->uv_max[0], tile->uv_max[1]},
        {tile->uv_min[0], tile->uv_max[1]}
    };
    long double value[4];
    const long double scale = 1.0L + fabsl(c) +
        fabsl(a) * fmaxl(fabsl(tile->uv_min[0]),
                         fabsl(tile->uv_max[0])) +
        fabsl(b) * fmaxl(fabsl(tile->uv_min[1]),
                         fabsl(tile->uv_max[1]));
    const long double tolerance = scale *
        (128.0L * LDBL_EPSILON + 16.0L * DBL_EPSILON);
    if (!isfinite(tolerance)) {
        *reason = ALEA_SLICE_ERROR_UNRESOLVED_NUMERICAL;
        return 0;
    }
    for (int i = 0; i < 4; ++i) {
        value[i] = a * corner[i][0] + b * corner[i][1] + c;
        if (!isfinite(value[i]) || fabsl(value[i]) <= tolerance) {
            *reason = ALEA_SLICE_ERROR_UNRESOLVED_NUMERICAL;
            return 0;
        }
    }
    slice_error_polygon_t polygon[2] = {{0}, {0}};
    slice_error_polygon_t crossings = {0};
    for (int i = 0; i < 4; ++i) {
        const int next = (i + 1) & 3;
        const int side = value[i] > 0.0L;
        if (!slice_error_polygon_add(&polygon[side],
                corner[i][0], corner[i][1], 0.0)) {
            *reason = ALEA_SLICE_ERROR_UNRESOLVED_CANDIDATE_LIMIT;
            return 0;
        }
        if ((value[i] > 0.0L) == (value[next] > 0.0L)) continue;
        const long double t = value[i] / (value[i] - value[next]);
        const long double u = (long double)corner[i][0] + t *
            ((long double)corner[next][0] - corner[i][0]);
        const long double v = (long double)corner[i][1] + t *
            ((long double)corner[next][1] - corner[i][1]);
        const double du = (double)u, dv = (double)v;
        const double uncertainty = (double)(
            4.0L * tolerance / fabsl(value[i] - value[next]) *
            fmaxl(fabsl((long double)corner[next][0] - corner[i][0]),
                  fabsl((long double)corner[next][1] - corner[i][1])) +
            8.0L * DBL_EPSILON *
                fmaxl(1.0L, fmaxl(fabsl(u), fabsl(v))));
        const long double edge_length = fmaxl(
            fabsl((long double)corner[next][0] - corner[i][0]),
            fabsl((long double)corner[next][1] - corner[i][1]));
        const long double corner_distance =
            fminl(t, 1.0L - t) * edge_length;
        if (!(t > 0.0L && t < 1.0L) || !isfinite(du) ||
            !isfinite(dv) || !isfinite(uncertainty) ||
            !((long double)uncertainty * 4.0L < corner_distance) ||
            !slice_error_polygon_add(&crossings, du, dv, uncertainty) ||
            !slice_error_polygon_add(&polygon[0], du, dv, uncertainty) ||
            !slice_error_polygon_add(&polygon[1], du, dv, uncertainty)) {
            *reason = ALEA_SLICE_ERROR_UNRESOLVED_NUMERICAL;
            return 0;
        }
    }
    if (crossings.count != 0 && crossings.count != 2) {
        *reason = ALEA_SLICE_ERROR_UNRESOLVED_NUMERICAL;
        return 0;
    }
    slice_error_face_t face[2] = {{0}, {0}};
    for (int side = 0; side < 2; ++side) {
        if (polygon[side].count == 0) continue;
        if (polygon[side].count < 3) {
            *reason = ALEA_SLICE_ERROR_UNRESOLVED_NUMERICAL;
            return 0;
        }
        double u = 0.0, v = 0.0;
        for (size_t j = 0; j < polygon[side].count; ++j) {
            u += polygon[side].uv[j][0] / polygon[side].count;
            v += polygon[side].uv[j][1] / polygon[side].count;
        }
        const long double residual = a * (long double)u +
                                     b * (long double)v + c;
        if (!isfinite(u) || !isfinite(v) ||
            (side ? residual <= tolerance : residual >= -tolerance)) {
            *reason = ALEA_SLICE_ERROR_UNRESOLVED_NUMERICAL;
            return 0;
        }
        const double wx = frame->origin[0] + frame->u_axis[0] * u +
                          frame->v_axis[0] * v;
        const double wy = frame->origin[1] + frame->u_axis[1] * u +
                          frame->v_axis[1] * v;
        if (!isfinite(wx) || !isfinite(wy)) {
            *reason = ALEA_SLICE_ERROR_UNRESOLVED_NUMERICAL;
            return 0;
        }
        for (size_t ci = 0; ci < cell_count; ++ci) {
            const alea_cell_entry_t* cell = &sys->cells.data[cells[ci]];
            int inside = 0;
            const uint32_t primitive_id = selected->primitive_id;
            const int raw_negative = side == 0;
            if (!slice_error_oblique_node_inside(sys, cell->root_node_id,
                    &primitive_id, &raw_negative, 1, 0, work, &inside)) {
                *reason = *work == 0
                    ? ALEA_SLICE_ERROR_UNRESOLVED_CANDIDATE_LIMIT
                    : ALEA_SLICE_ERROR_UNRESOLVED_UNSUPPORTED_GEOMETRY;
                return 0;
            }
            if (inside != alea_point_inside(sys, cell->root_node_id,
                                            wx, wy, frame->origin[2])) {
                *reason = ALEA_SLICE_ERROR_UNRESOLVED_NUMERICAL;
                return 0;
            }
            if (inside) {
                if (face[side].owner_count == ALEA_SLICE_ERROR_OWNER_CAPACITY) {
                    *reason = ALEA_SLICE_ERROR_UNRESOLVED_CANDIDATE_LIMIT;
                    return 0;
                }
                face[side].owner_cell_ids[face[side].owner_count++] =
                    cell->mc_cell_id;
            }
        }
        face[side].kind = face[side].owner_count == 0
            ? ALEA_POINT_COVERAGE_GAP
            : face[side].owner_count == 1
                ? ALEA_POINT_COVERAGE_UNIQUE
                : ALEA_POINT_COVERAGE_OVERLAP;
    }
    size_t region_count = 0;
    for (int side = 0; side < 2; ++side)
        if (polygon[side].count && slice_error_axis_is_defect(face[side].kind))
            ++region_count;
    const size_t interval_count = crossings.count == 2 &&
        slice_error_axis_is_defect(face[0].kind) !=
        slice_error_axis_is_defect(face[1].kind) ? 1u : 0u;
    const size_t output_limit = o->scan_options.max_output_bytes;
    if (contextual_bytes > output_limit ||
        region_count > (output_limit - contextual_bytes) /
            sizeof(alea_slice_error_region_t) ||
        interval_count > (output_limit - contextual_bytes -
            region_count * sizeof(alea_slice_error_region_t)) /
            sizeof(alea_slice_error_interval_t)) {
        *reason = ALEA_SLICE_ERROR_UNRESOLVED_CANDIDATE_LIMIT;
        *output_omitted = 1;
        return 0;
    }
    out->regions = region_count
        ? calloc(region_count, sizeof(*out->regions)) : NULL;
    out->intervals = interval_count
        ? calloc(interval_count, sizeof(*out->intervals)) : NULL;
    if ((region_count && !out->regions) ||
        (interval_count && !out->intervals)) return -1;
    for (int side = 0; side < 2; ++side) {
        if (!polygon[side].count ||
            !slice_error_axis_is_defect(face[side].kind)) continue;
        alea_slice_error_region_t* region =
            &out->regions[out->region_count++];
        region->kind = face[side].kind;
        region->owner_count = face[side].owner_count;
        memcpy(region->owner_cell_ids, face[side].owner_cell_ids,
               region->owner_count * sizeof(int));
        region->polygon_vertex_count = polygon[side].count;
        region->uv_min[0] = region->uv_min[1] = INFINITY;
        region->uv_max[0] = region->uv_max[1] = -INFINITY;
        for (size_t j = 0; j < polygon[side].count; ++j) {
            for (int axis = 0; axis < 2; ++axis) {
                const double coordinate = polygon[side].uv[j][axis];
                region->polygon_uv[j][axis] = coordinate;
                region->polygon_uv_uncertainty[j][axis] =
                    polygon[side].uncertainty[j][axis];
                region->uv_min[axis] = fmin(region->uv_min[axis], coordinate);
                region->uv_max[axis] = fmax(region->uv_max[axis], coordinate);
                if (region->uv_min[axis] == coordinate)
                    region->uv_min_uncertainty[axis] = fmax(
                        region->uv_min_uncertainty[axis],
                        polygon[side].uncertainty[j][axis]);
                if (region->uv_max[axis] == coordinate)
                    region->uv_max_uncertainty[axis] = fmax(
                        region->uv_max_uncertainty[axis],
                        polygon[side].uncertainty[j][axis]);
            }
        }
    }
    if (interval_count) {
        alea_slice_error_interval_t* interval = &out->intervals[0];
        interval->evidence_scope =
            ALEA_SLICE_BOUNDARY_EVIDENCE_VERIFIED_INTERVAL;
        interval->surface_id = selected->surface_id;
        interval->primitive_id = selected->primitive_id;
        interval->axis = -1;
        const int reverse = crossings.uv[0][0] > crossings.uv[1][0] ||
            (crossings.uv[0][0] == crossings.uv[1][0] &&
             crossings.uv[0][1] > crossings.uv[1][1]);
        const int first = reverse ? 1 : 0, last = reverse ? 0 : 1;
        memcpy(interval->uv_start, crossings.uv[first],
               sizeof(interval->uv_start));
        memcpy(interval->uv_end, crossings.uv[last],
               sizeof(interval->uv_end));
        interval->endpoint_uncertainty[0] =
            crossings.uncertainty[first][0];
        interval->endpoint_uncertainty[1] =
            crossings.uncertainty[last][0];
        interval->negative_side_kind = face[0].kind;
        interval->positive_side_kind = face[1].kind;
        interval->negative_owner_count = face[0].owner_count;
        interval->positive_owner_count = face[1].owner_count;
        memcpy(interval->negative_owner_cell_ids, face[0].owner_cell_ids,
               face[0].owner_count * sizeof(int));
        memcpy(interval->positive_owner_cell_ids, face[1].owner_cell_ids,
               face[1].owner_count * sizeof(int));
        out->interval_count = 1;
    }
    *reason = ALEA_SLICE_ERROR_RESOLVED;
    return 1;
}

typedef struct {
    int plane;
    int negative_face, positive_face;
    double uv[2][2];
    double uncertainty[2];
} slice_error_oblique_segment_t;

static int slice_error_publish_oblique_faces(
    const alea_slice_error_query_options_t* options,
    const slice_error_axis_plane_t* selected,
    const slice_error_polygon_t* polygon,
    const slice_error_face_t* face, size_t face_count,
    const slice_error_oblique_segment_t* segments, size_t segment_count,
    size_t contextual_bytes, alea_slice_error_page_t* out,
    alea_slice_error_unresolved_reason_t* reason, int* output_omitted) {
    size_t region_count = 0, interval_count = 0;
    for (size_t i = 0; i < segment_count; ++i)
        if (segments[i].negative_face < 0 ||
            segments[i].positive_face < 0 ||
            (size_t)segments[i].negative_face >= face_count ||
            (size_t)segments[i].positive_face >= face_count ||
            !polygon[segments[i].negative_face].count ||
            !polygon[segments[i].positive_face].count) {
            *reason = ALEA_SLICE_ERROR_UNRESOLVED_NUMERICAL;
            return 0;
        }
    for (size_t i = 0; i < face_count; ++i)
        if (polygon[i].count && slice_error_axis_is_defect(face[i].kind))
            ++region_count;
    for (size_t i = 0; i < segment_count; ++i)
        if (slice_error_axis_is_defect(
                face[segments[i].negative_face].kind) !=
            slice_error_axis_is_defect(
                face[segments[i].positive_face].kind)) ++interval_count;
    const size_t output_limit = options->scan_options.max_output_bytes;
    if (contextual_bytes > output_limit ||
        region_count > (output_limit - contextual_bytes) /
            sizeof(alea_slice_error_region_t) ||
        interval_count > (output_limit - contextual_bytes -
            region_count * sizeof(alea_slice_error_region_t)) /
            sizeof(alea_slice_error_interval_t)) {
        *reason = ALEA_SLICE_ERROR_UNRESOLVED_CANDIDATE_LIMIT;
        *output_omitted = 1;
        return 0;
    }
    out->regions = region_count
        ? calloc(region_count, sizeof(*out->regions)) : NULL;
    out->intervals = interval_count
        ? calloc(interval_count, sizeof(*out->intervals)) : NULL;
    if ((region_count && !out->regions) ||
        (interval_count && !out->intervals)) return -1;
    for (size_t i = 0; i < face_count; ++i) {
        if (!polygon[i].count ||
            !slice_error_axis_is_defect(face[i].kind)) continue;
        alea_slice_error_region_t* region =
            &out->regions[out->region_count++];
        region->kind = face[i].kind;
        region->owner_count = face[i].owner_count;
        memcpy(region->owner_cell_ids, face[i].owner_cell_ids,
               region->owner_count * sizeof(int));
        region->polygon_vertex_count = polygon[i].count;
        region->uv_min[0] = region->uv_min[1] = INFINITY;
        region->uv_max[0] = region->uv_max[1] = -INFINITY;
        for (size_t j = 0; j < polygon[i].count; ++j)
            for (int axis = 0; axis < 2; ++axis) {
                const double coordinate = polygon[i].uv[j][axis];
                const double uncertainty = polygon[i].uncertainty[j][axis];
                region->polygon_uv[j][axis] = coordinate;
                region->polygon_uv_uncertainty[j][axis] = uncertainty;
                region->uv_min[axis] =
                    fmin(region->uv_min[axis], coordinate);
                region->uv_max[axis] =
                    fmax(region->uv_max[axis], coordinate);
                if (region->uv_min[axis] == coordinate)
                    region->uv_min_uncertainty[axis] = fmax(
                        region->uv_min_uncertainty[axis], uncertainty);
                if (region->uv_max[axis] == coordinate)
                    region->uv_max_uncertainty[axis] = fmax(
                        region->uv_max_uncertainty[axis], uncertainty);
            }
    }
    for (size_t i = 0; i < segment_count; ++i) {
        const slice_error_oblique_segment_t* segment = &segments[i];
        const slice_error_face_t* negative =
            &face[segment->negative_face];
        const slice_error_face_t* positive =
            &face[segment->positive_face];
        if (slice_error_axis_is_defect(negative->kind) ==
            slice_error_axis_is_defect(positive->kind)) continue;
        alea_slice_error_interval_t* interval =
            &out->intervals[out->interval_count++];
        interval->evidence_scope =
            ALEA_SLICE_BOUNDARY_EVIDENCE_VERIFIED_INTERVAL;
        interval->surface_id = selected[segment->plane].surface_id;
        interval->primitive_id = selected[segment->plane].primitive_id;
        interval->axis = -1;
        const int reverse = segment->uv[0][0] > segment->uv[1][0] ||
            (segment->uv[0][0] == segment->uv[1][0] &&
             segment->uv[0][1] > segment->uv[1][1]);
        const int first = reverse ? 1 : 0, last = reverse ? 0 : 1;
        memcpy(interval->uv_start, segment->uv[first],
               sizeof(interval->uv_start));
        memcpy(interval->uv_end, segment->uv[last],
               sizeof(interval->uv_end));
        interval->endpoint_uncertainty[0] =
            segment->uncertainty[first];
        interval->endpoint_uncertainty[1] =
            segment->uncertainty[last];
        interval->negative_side_kind = negative->kind;
        interval->positive_side_kind = positive->kind;
        interval->negative_owner_count = negative->owner_count;
        interval->positive_owner_count = positive->owner_count;
        memcpy(interval->negative_owner_cell_ids, negative->owner_cell_ids,
               negative->owner_count * sizeof(int));
        memcpy(interval->positive_owner_cell_ids, positive->owner_cell_ids,
               positive->owner_count * sizeof(int));
    }
    *reason = ALEA_SLICE_ERROR_RESOLVED;
    return 1;
}

/* Two parallel oblique planes divide the core into at most three convex
 * bands. The primitive signs are fixed in each band, including a thin band
 * between the lines. Nonparallel or nearly coincident lines remain
 * unresolved. */
static int slice_error_classify_parallel_oblique_tile(
    const alea_slice_error_query_t* query,
    const alea_transition_slice_critical_tile_t* tile,
    const size_t* cells, size_t cell_count,
    const slice_error_axis_plane_t selected[2],
    size_t contextual_bytes, size_t* work,
    alea_slice_error_page_t* out,
    alea_slice_error_unresolved_reason_t* reason,
    int* output_omitted) {
    const alea_slice_error_query_options_t* o = &query->options;
    const alea_system_t* sys = query->sys;
    const alea_slice_plane_t* frame = &o->view.plane;
    long double coefficient[2][2], threshold[2];
    double normal[2][2];
    int orientation[2] = {1, 0};
    uint32_t primitive_ids[2];
    for (int i = 0; i < 2; ++i) {
        if (selected[i].surface_id <= 0) return 0;
        primitive_ids[i] = selected[i].primitive_id;
        const alea_primitive_entry_t* primitive =
            &sys->primitives.data[primitive_ids[i]];
        if (primitive->type != ALEA_PRIMITIVE_PLANE ||
            primitive->payload_index >= sys->primitive_planes.count)
            return 0;
        const alea_plane_data_t* p =
            &sys->primitive_planes.data[primitive->payload_index];
        if (!isfinite(p->a) || !isfinite(p->b) || !isfinite(p->c) ||
            !isfinite(p->d) || p->a == 0.0 || p->b == 0.0 ||
            p->c != 0.0) return 0;
        normal[i][0] = p->a;
        normal[i][1] = p->b;
        const long double a = (long double)p->a * frame->u_axis[0] +
                              (long double)p->b * frame->u_axis[1];
        const long double b = (long double)p->a * frame->v_axis[0] +
                              (long double)p->b * frame->v_axis[1];
        const long double c = (long double)p->a * frame->origin[0] +
                              (long double)p->b * frame->origin[1] +
                              (long double)p->d;
        const long double norm = fmaxl(fabsl(a), fabsl(b));
        if (!(norm > 0.0L) || !isfinite(norm)) {
            *reason = ALEA_SLICE_ERROR_UNRESOLVED_NUMERICAL;
            return 0;
        }
        coefficient[i][0] = a / norm;
        coefficient[i][1] = b / norm;
        threshold[i] = -c / norm;
        if (!isfinite(threshold[i])) {
            *reason = ALEA_SLICE_ERROR_UNRESOLVED_NUMERICAL;
            return 0;
        }
    }
    /* Exact matching normals avoid treating two nearly parallel lines as
     * disjoint bands when their intersection is hidden inside the tile. */
    if (normal[1][0] == normal[0][0] &&
        normal[1][1] == normal[0][1])
        orientation[1] = 1;
    else if (normal[1][0] == -normal[0][0] &&
             normal[1][1] == -normal[0][1]) {
        orientation[1] = -1;
        threshold[1] = -threshold[1];
    } else return 0;
    const int line[2] = {
        threshold[0] < threshold[1] ? 0 : 1,
        threshold[0] < threshold[1] ? 1 : 0
    };
    const long double bound[2] = {
        threshold[line[0]], threshold[line[1]]
    };
    const double corner[4][2] = {
        {tile->uv_min[0], tile->uv_min[1]},
        {tile->uv_max[0], tile->uv_min[1]},
        {tile->uv_max[0], tile->uv_max[1]},
        {tile->uv_min[0], tile->uv_max[1]}
    };
    long double value[4];
    const long double a = coefficient[0][0], b = coefficient[0][1];
    const long double scale = 1.0L +
        fmaxl(fabsl(bound[0]), fabsl(bound[1])) +
        fabsl(a) * fmaxl(fabsl(tile->uv_min[0]),
                         fabsl(tile->uv_max[0])) +
        fabsl(b) * fmaxl(fabsl(tile->uv_min[1]),
                         fabsl(tile->uv_max[1]));
    const long double tolerance = scale *
        (128.0L * LDBL_EPSILON + 16.0L * DBL_EPSILON);
    if (!isfinite(tolerance) || !(bound[1] - bound[0] >
                                  8.0L * tolerance)) {
        *reason = ALEA_SLICE_ERROR_UNRESOLVED_NUMERICAL;
        return 0;
    }
    for (int i = 0; i < 4; ++i) {
        value[i] = a * corner[i][0] + b * corner[i][1];
        if (!isfinite(value[i]) ||
            fabsl(value[i] - bound[0]) <= tolerance ||
            fabsl(value[i] - bound[1]) <= tolerance) {
            *reason = ALEA_SLICE_ERROR_UNRESOLVED_NUMERICAL;
            return 0;
        }
    }
    slice_error_polygon_t polygon[3] = {{0}, {0}, {0}};
    slice_error_polygon_t crossings[2] = {{0}, {0}};
    for (int i = 0; i < 4; ++i) {
        const int next = (i + 1) & 3;
        const int band = value[i] < bound[0] ? 0
            : value[i] < bound[1] ? 1 : 2;
        if (!slice_error_polygon_add(&polygon[band],
                corner[i][0], corner[i][1], 0.0)) {
            *reason = ALEA_SLICE_ERROR_UNRESOLVED_CANDIDATE_LIMIT;
            return 0;
        }
        const int increasing = value[next] > value[i];
        for (int step = 0; step < 2; ++step) {
            const int k = increasing ? step : 1 - step;
            if ((value[i] < bound[k]) == (value[next] < bound[k]))
                continue;
            const long double t =
                (bound[k] - value[i]) / (value[next] - value[i]);
            const long double du_edge =
                (long double)corner[next][0] - corner[i][0];
            const long double dv_edge =
                (long double)corner[next][1] - corner[i][1];
            const long double u = (long double)corner[i][0] + t * du_edge;
            const long double v = (long double)corner[i][1] + t * dv_edge;
            const double du = (double)u, dv = (double)v;
            const long double edge_length =
                fmaxl(fabsl(du_edge), fabsl(dv_edge));
            const double uncertainty = (double)(
                4.0L * tolerance / fabsl(value[next] - value[i]) *
                    edge_length +
                8.0L * DBL_EPSILON *
                    fmaxl(1.0L, fmaxl(fabsl(u), fabsl(v))));
            const long double corner_distance =
                fminl(t, 1.0L - t) * edge_length;
            if (!(t > 0.0L && t < 1.0L) || !isfinite(du) ||
                !isfinite(dv) || !isfinite(uncertainty) ||
                !((long double)uncertainty * 4.0L < corner_distance) ||
                !slice_error_polygon_add(&crossings[k], du, dv,
                                         uncertainty) ||
                !slice_error_polygon_add(&polygon[k], du, dv,
                                         uncertainty) ||
                !slice_error_polygon_add(&polygon[k + 1], du, dv,
                                         uncertainty)) {
                *reason = ALEA_SLICE_ERROR_UNRESOLVED_NUMERICAL;
                return 0;
            }
        }
    }
    for (int k = 0; k < 2; ++k)
        if (crossings[k].count != 0 && crossings[k].count != 2) {
            *reason = ALEA_SLICE_ERROR_UNRESOLVED_NUMERICAL;
            return 0;
        }
    slice_error_face_t face[3] = {{0}, {0}, {0}};
    for (int band = 0; band < 3; ++band) {
        if (!polygon[band].count) continue;
        if (polygon[band].count < 3) {
            *reason = ALEA_SLICE_ERROR_UNRESOLVED_NUMERICAL;
            return 0;
        }
        double u = 0.0, v = 0.0;
        for (size_t j = 0; j < polygon[band].count; ++j) {
            u += polygon[band].uv[j][0] / polygon[band].count;
            v += polygon[band].uv[j][1] / polygon[band].count;
        }
        const long double s = a * (long double)u + b * (long double)v;
        if (!isfinite(u) || !isfinite(v) ||
            (band > 0 && !(s - bound[band - 1] > tolerance)) ||
            (band < 2 && !(bound[band] - s > tolerance))) {
            *reason = ALEA_SLICE_ERROR_UNRESOLVED_NUMERICAL;
            return 0;
        }
        const double wx = frame->origin[0] + frame->u_axis[0] * u +
                          frame->v_axis[0] * v;
        const double wy = frame->origin[1] + frame->u_axis[1] * u +
                          frame->v_axis[1] * v;
        if (!isfinite(wx) || !isfinite(wy)) {
            *reason = ALEA_SLICE_ERROR_UNRESOLVED_NUMERICAL;
            return 0;
        }
        int raw_negative[2];
        for (int i = 0; i < 2; ++i)
            raw_negative[i] = orientation[i] > 0
                ? s < threshold[i] : s > threshold[i];
        for (size_t ci = 0; ci < cell_count; ++ci) {
            const alea_cell_entry_t* cell = &sys->cells.data[cells[ci]];
            int inside = 0;
            if (!slice_error_oblique_node_inside(sys,
                    cell->root_node_id, primitive_ids, raw_negative,
                    2, 0, work, &inside)) {
                *reason = *work == 0
                    ? ALEA_SLICE_ERROR_UNRESOLVED_CANDIDATE_LIMIT
                    : ALEA_SLICE_ERROR_UNRESOLVED_UNSUPPORTED_GEOMETRY;
                return 0;
            }
            if (inside != alea_point_inside(sys, cell->root_node_id,
                                            wx, wy, frame->origin[2])) {
                *reason = ALEA_SLICE_ERROR_UNRESOLVED_NUMERICAL;
                return 0;
            }
            if (inside) {
                if (face[band].owner_count ==
                    ALEA_SLICE_ERROR_OWNER_CAPACITY) {
                    *reason = ALEA_SLICE_ERROR_UNRESOLVED_CANDIDATE_LIMIT;
                    return 0;
                }
                face[band].owner_cell_ids[face[band].owner_count++] =
                    cell->mc_cell_id;
            }
        }
        face[band].kind = face[band].owner_count == 0
            ? ALEA_POINT_COVERAGE_GAP
            : face[band].owner_count == 1
                ? ALEA_POINT_COVERAGE_UNIQUE
                : ALEA_POINT_COVERAGE_OVERLAP;
    }
    slice_error_oblique_segment_t segments[2] = {{0}, {0}};
    size_t segment_count = 0;
    for (int k = 0; k < 2; ++k) {
        if (crossings[k].count != 2) continue;
        const int plane = line[k];
        slice_error_oblique_segment_t* segment =
            &segments[segment_count++];
        segment->plane = plane;
        segment->negative_face = orientation[plane] > 0 ? k : k + 1;
        segment->positive_face = orientation[plane] > 0 ? k + 1 : k;
        memcpy(segment->uv, crossings[k].uv, sizeof(segment->uv));
        segment->uncertainty[0] = crossings[k].uncertainty[0][0];
        segment->uncertainty[1] = crossings[k].uncertainty[1][0];
    }
    return slice_error_publish_oblique_faces(
        o, selected, polygon, face, 3, segments, segment_count,
        contextual_bytes, out, reason, output_omitted);
}

/* Two nonparallel vertical lines partition a rectangular core into at most four convex
 * faces. An interior crossing is a vertex of all four faces. A crossing well
 * outside the core has no face vertex. Boundary or poorly conditioned
 * crossings are withheld. */
static int slice_error_classify_crossing_oblique_tile(
    const alea_slice_error_query_t* query,
    const alea_transition_slice_critical_tile_t* tile,
    const size_t* cells, size_t cell_count,
    const slice_error_axis_plane_t selected[2],
    size_t contextual_bytes, size_t* work,
    alea_slice_error_page_t* out,
    alea_slice_error_unresolved_reason_t* reason,
    int* output_omitted) {
    const alea_slice_error_query_options_t* o = &query->options;
    const alea_system_t* sys = query->sys;
    const alea_slice_plane_t* frame = &o->view.plane;
    uint32_t primitive_ids[2];
    long double a[2], b[2], c[2], tolerance[2];
    const double corner[4][2] = {
        {tile->uv_min[0], tile->uv_min[1]},
        {tile->uv_max[0], tile->uv_min[1]},
        {tile->uv_max[0], tile->uv_max[1]},
        {tile->uv_min[0], tile->uv_max[1]}
    };
    long double value[2][4];
    for (int i = 0; i < 2; ++i) {
        if (selected[i].surface_id <= 0) return 0;
        primitive_ids[i] = selected[i].primitive_id;
        const alea_primitive_entry_t* primitive =
            &sys->primitives.data[primitive_ids[i]];
        if (primitive->type != ALEA_PRIMITIVE_PLANE ||
            primitive->payload_index >= sys->primitive_planes.count)
            return 0;
        const alea_plane_data_t* p =
            &sys->primitive_planes.data[primitive->payload_index];
        if (!isfinite(p->a) || !isfinite(p->b) || !isfinite(p->c) ||
            !isfinite(p->d) || (p->a == 0.0 && p->b == 0.0) ||
            p->c != 0.0) return 0;
        a[i] = (long double)p->a * frame->u_axis[0] +
               (long double)p->b * frame->u_axis[1];
        b[i] = (long double)p->a * frame->v_axis[0] +
               (long double)p->b * frame->v_axis[1];
        c[i] = (long double)p->a * frame->origin[0] +
               (long double)p->b * frame->origin[1] +
               (long double)p->d;
        const long double norm = fmaxl(fabsl(a[i]), fabsl(b[i]));
        if (!(norm > 0.0L) || !isfinite(norm)) {
            *reason = ALEA_SLICE_ERROR_UNRESOLVED_NUMERICAL;
            return 0;
        }
        a[i] /= norm; b[i] /= norm; c[i] /= norm;
        const long double scale = 1.0L + fabsl(c[i]) +
            fabsl(a[i]) * fmaxl(fabsl(tile->uv_min[0]),
                                fabsl(tile->uv_max[0])) +
            fabsl(b[i]) * fmaxl(fabsl(tile->uv_min[1]),
                                fabsl(tile->uv_max[1]));
        tolerance[i] = scale *
            (128.0L * LDBL_EPSILON + 16.0L * DBL_EPSILON);
        if (!isfinite(tolerance[i])) {
            *reason = ALEA_SLICE_ERROR_UNRESOLVED_NUMERICAL;
            return 0;
        }
        for (int j = 0; j < 4; ++j) {
            value[i][j] = a[i] * corner[j][0] +
                          b[i] * corner[j][1] + c[i];
            if (!isfinite(value[i][j]) ||
                fabsl(value[i][j]) <= tolerance[i]) {
                *reason = ALEA_SLICE_ERROR_UNRESOLVED_NUMERICAL;
                return 0;
            }
        }
    }
    const long double det = a[0] * b[1] - a[1] * b[0];
    if (!isfinite(det) || fabsl(det) <=
        128.0L * LDBL_EPSILON) return 0;
    const long double center_u = (b[0] * c[1] - b[1] * c[0]) / det;
    const long double center_v = (a[1] * c[0] - a[0] * c[1]) / det;
    const double crossing[2] = {(double)center_u, (double)center_v};
    const double center_uncertainty = (double)(
        32.0L * (tolerance[0] + tolerance[1]) / fabsl(det) +
        8.0L * DBL_EPSILON *
            fmaxl(1.0L, fmaxl(fabsl(center_u), fabsl(center_v))));
    if (!isfinite(center_u) || !isfinite(center_v) ||
        !isfinite(crossing[0]) || !isfinite(crossing[1]) ||
        !isfinite(center_uncertainty)) {
        *reason = ALEA_SLICE_ERROR_UNRESOLVED_NUMERICAL;
        return 0;
    }
    int center_inside = 1, center_outside = 0;
    for (int axis = 0; axis < 2; ++axis) {
        if (!(crossing[axis] - tile->uv_min[axis] >
              4.0 * center_uncertainty &&
              tile->uv_max[axis] - crossing[axis] >
              4.0 * center_uncertainty)) center_inside = 0;
        if (crossing[axis] < tile->uv_min[axis] -
                4.0 * center_uncertainty ||
            crossing[axis] > tile->uv_max[axis] +
                4.0 * center_uncertainty) center_outside = 1;
    }
    if (!center_inside && !center_outside) {
        *reason = ALEA_SLICE_ERROR_UNRESOLVED_NUMERICAL;
        return 0;
    }
    slice_error_polygon_t polygon[4] = {{0}, {0}, {0}, {0}};
    slice_error_polygon_t edge_hit[2] = {{0}, {0}};
    for (int j = 0; j < 4; ++j) {
        const int mask = (value[0][j] > 0.0L ? 1 : 0) |
                         (value[1][j] > 0.0L ? 2 : 0);
        if (!slice_error_polygon_add(&polygon[mask],
                corner[j][0], corner[j][1], 0.0)) {
            *reason = ALEA_SLICE_ERROR_UNRESOLVED_CANDIDATE_LIMIT;
            return 0;
        }
        const int next = (j + 1) & 3;
        for (int i = 0; i < 2; ++i) {
            if ((value[i][j] > 0.0L) == (value[i][next] > 0.0L))
                continue;
            const long double t = value[i][j] /
                (value[i][j] - value[i][next]);
            const long double du_edge =
                (long double)corner[next][0] - corner[j][0];
            const long double dv_edge =
                (long double)corner[next][1] - corner[j][1];
            const long double u = corner[j][0] + t * du_edge;
            const long double v = corner[j][1] + t * dv_edge;
            const double du = (double)u, dv = (double)v;
            const long double edge_length =
                fmaxl(fabsl(du_edge), fabsl(dv_edge));
            const double uncertainty = (double)(
                4.0L * tolerance[i] /
                    fabsl(value[i][j] - value[i][next]) * edge_length +
                8.0L * DBL_EPSILON *
                    fmaxl(1.0L, fmaxl(fabsl(u), fabsl(v))));
            const int other = 1 - i;
            const long double other_value =
                a[other] * u + b[other] * v + c[other];
            if (!(t > 0.0L && t < 1.0L) || !isfinite(du) ||
                !isfinite(dv) || !isfinite(uncertainty) ||
                !((long double)uncertainty * 4.0L <
                  fminl(t, 1.0L - t) * edge_length) ||
                !(fabsl(other_value) > 4.0L *
                    (tolerance[other] +
                     (fabsl(a[other]) + fabsl(b[other])) *
                        uncertainty))) {
                *reason = ALEA_SLICE_ERROR_UNRESOLVED_NUMERICAL;
                return 0;
            }
            const int other_bit = other_value > 0.0L ? 1 : 0;
            const int negative_face = i == 0 ? other_bit << 1 : other_bit;
            const int positive_face = negative_face | (1 << i);
            if (!slice_error_polygon_add(&edge_hit[i], du, dv,
                                         uncertainty) ||
                !slice_error_polygon_add(&polygon[negative_face], du, dv,
                                         uncertainty) ||
                !slice_error_polygon_add(&polygon[positive_face], du, dv,
                                         uncertainty)) {
                *reason = ALEA_SLICE_ERROR_UNRESOLVED_CANDIDATE_LIMIT;
                return 0;
            }
        }
    }
    for (int i = 0; i < 2; ++i) {
        if (edge_hit[i].count != 0 && edge_hit[i].count != 2) {
            *reason = ALEA_SLICE_ERROR_UNRESOLVED_NUMERICAL;
            return 0;
        }
        if (center_inside && edge_hit[i].count != 2) return 0;
        if (center_inside)
            for (int j = 0; j < 2; ++j) {
                const double separation = fmax(
                    fabs(crossing[0] - edge_hit[i].uv[j][0]),
                    fabs(crossing[1] - edge_hit[i].uv[j][1]));
                if (!(separation > 4.0 * (center_uncertainty +
                            edge_hit[i].uncertainty[j][0]))) {
                    *reason = ALEA_SLICE_ERROR_UNRESOLVED_NUMERICAL;
                    return 0;
                }
            }
    }
    slice_error_face_t face[4] = {{0}, {0}, {0}, {0}};
    for (int mask = 0; mask < 4; ++mask) {
        if (center_inside && !slice_error_polygon_add(&polygon[mask],
                crossing[0], crossing[1], center_uncertainty)) {
            *reason = ALEA_SLICE_ERROR_UNRESOLVED_NUMERICAL;
            return 0;
        }
        if (!polygon[mask].count) continue;
        if (polygon[mask].count < 3) {
            *reason = ALEA_SLICE_ERROR_UNRESOLVED_NUMERICAL;
            return 0;
        }
        double centroid_u = 0.0, centroid_v = 0.0;
        for (size_t j = 0; j < polygon[mask].count; ++j) {
            centroid_u += polygon[mask].uv[j][0] / polygon[mask].count;
            centroid_v += polygon[mask].uv[j][1] / polygon[mask].count;
        }
        /* Every point belongs to this convex face. Angular order gives its
         * boundary, including a face with no rectangle corner. */
        for (size_t j = 1; j < polygon[mask].count; ++j) {
            const double u = polygon[mask].uv[j][0];
            const double v = polygon[mask].uv[j][1];
            const double uncertainty = polygon[mask].uncertainty[j][0];
            const double angle = atan2(v - centroid_v, u - centroid_u);
            size_t k = j;
            while (k > 0 && atan2(
                    polygon[mask].uv[k - 1][1] - centroid_v,
                    polygon[mask].uv[k - 1][0] - centroid_u) > angle) {
                memcpy(polygon[mask].uv[k], polygon[mask].uv[k - 1],
                       sizeof(polygon[mask].uv[k]));
                memcpy(polygon[mask].uncertainty[k],
                       polygon[mask].uncertainty[k - 1],
                       sizeof(polygon[mask].uncertainty[k]));
                --k;
            }
            polygon[mask].uv[k][0] = u;
            polygon[mask].uv[k][1] = v;
            polygon[mask].uncertainty[k][0] = uncertainty;
            polygon[mask].uncertainty[k][1] = uncertainty;
        }
        for (size_t j = 0; j < polygon[mask].count; ++j)
            for (size_t k = j + 1; k < polygon[mask].count; ++k) {
                const double separation = fmax(
                    fabs(polygon[mask].uv[j][0] - polygon[mask].uv[k][0]),
                    fabs(polygon[mask].uv[j][1] - polygon[mask].uv[k][1]));
                if (!(separation > 4.0 *
                    (polygon[mask].uncertainty[j][0] +
                     polygon[mask].uncertainty[k][0]))) {
                    *reason = ALEA_SLICE_ERROR_UNRESOLVED_NUMERICAL;
                    return 0;
                }
            }
        const long double residual[2] = {
            a[0] * centroid_u + b[0] * centroid_v + c[0],
            a[1] * centroid_u + b[1] * centroid_v + c[1]
        };
        for (int i = 0; i < 2; ++i)
            if (((mask & (1 << i)) != 0)
                    ? !(residual[i] > tolerance[i])
                    : !(residual[i] < -tolerance[i])) {
                *reason = ALEA_SLICE_ERROR_UNRESOLVED_NUMERICAL;
                return 0;
            }
        const double wx = frame->origin[0] +
            frame->u_axis[0] * centroid_u +
            frame->v_axis[0] * centroid_v;
        const double wy = frame->origin[1] +
            frame->u_axis[1] * centroid_u +
            frame->v_axis[1] * centroid_v;
        if (!isfinite(wx) || !isfinite(wy)) {
            *reason = ALEA_SLICE_ERROR_UNRESOLVED_NUMERICAL;
            return 0;
        }
        const int raw_negative[2] = {
            (mask & 1) == 0, (mask & 2) == 0
        };
        for (size_t ci = 0; ci < cell_count; ++ci) {
            const alea_cell_entry_t* cell = &sys->cells.data[cells[ci]];
            int inside = 0;
            if (!slice_error_oblique_node_inside(sys,
                    cell->root_node_id, primitive_ids, raw_negative,
                    2, 0, work, &inside)) {
                *reason = *work == 0
                    ? ALEA_SLICE_ERROR_UNRESOLVED_CANDIDATE_LIMIT
                    : ALEA_SLICE_ERROR_UNRESOLVED_UNSUPPORTED_GEOMETRY;
                return 0;
            }
            if (inside != alea_point_inside(sys, cell->root_node_id,
                                            wx, wy, frame->origin[2])) {
                *reason = ALEA_SLICE_ERROR_UNRESOLVED_NUMERICAL;
                return 0;
            }
            if (inside) {
                if (face[mask].owner_count ==
                    ALEA_SLICE_ERROR_OWNER_CAPACITY) {
                    *reason = ALEA_SLICE_ERROR_UNRESOLVED_CANDIDATE_LIMIT;
                    return 0;
                }
                face[mask].owner_cell_ids[face[mask].owner_count++] =
                    cell->mc_cell_id;
            }
        }
        face[mask].kind = face[mask].owner_count == 0
            ? ALEA_POINT_COVERAGE_GAP
            : face[mask].owner_count == 1
                ? ALEA_POINT_COVERAGE_UNIQUE
                : ALEA_POINT_COVERAGE_OVERLAP;
    }
    slice_error_oblique_segment_t segments[4] = {{0}, {0}, {0}, {0}};
    size_t segment_count = 0;
    for (int i = 0; i < 2; ++i) {
        const int other = 1 - i;
        const int parts = center_inside ? 2
            : edge_hit[i].count == 2 ? 1 : 0;
        for (int j = 0; j < parts; ++j) {
            const double* start = edge_hit[i].uv[j];
            const double* end = center_inside
                ? crossing : edge_hit[i].uv[1];
            const long double mid_u =
                ((long double)start[0] + end[0]) * 0.5L;
            const long double mid_v =
                ((long double)start[1] + end[1]) * 0.5L;
            const long double other_value =
                a[other] * mid_u + b[other] * mid_v + c[other];
            if (!(fabsl(other_value) > 4.0L * tolerance[other])) {
                *reason = ALEA_SLICE_ERROR_UNRESOLVED_NUMERICAL;
                return 0;
            }
            const int other_bit = other_value > 0.0L ? 1 : 0;
            slice_error_oblique_segment_t* segment =
                &segments[segment_count++];
            segment->plane = i;
            segment->negative_face = i == 0
                ? other_bit << 1 : other_bit;
            segment->positive_face =
                segment->negative_face | (1 << i);
            memcpy(segment->uv[0], start,
                   sizeof(segment->uv[0]));
            memcpy(segment->uv[1], end,
                   sizeof(segment->uv[1]));
            segment->uncertainty[0] =
                edge_hit[i].uncertainty[j][0];
            segment->uncertainty[1] = center_inside
                ? center_uncertainty : edge_hit[i].uncertainty[1][0];
        }
    }
    return slice_error_publish_oblique_faces(
        o, selected, polygon, face, 4, segments, segment_count,
        contextual_bytes, out, reason, output_omitted);
}

/* Returns 1 for a complete proof, 0 with a whole-tile unresolved reason for
 * unsupported/numerically uncertain work, or -1 on allocation failure. */
static int slice_error_classify_axis_tile(
    const alea_slice_error_query_t* query,
    const alea_transition_slice_critical_tile_t* tile,
    size_t contextual_bytes, alea_slice_error_page_t* out,
    alea_slice_error_unresolved_reason_t* reason,
    int* output_omitted, size_t* peak_scratch) {
    *output_omitted = 0;
    *peak_scratch = 0;
    alea_system_t* sys = query->sys;
    const alea_slice_error_query_options_t* options = &query->options;
    size_t primitive_count = 0;
    size_t cell_count = 0;
    if (!slice_error_axis_view_supported(&options->view)) {
        *reason = ALEA_SLICE_ERROR_UNRESOLVED_UNSUPPORTED_GEOMETRY;
        return 0;
    }
    const size_t candidate_limit = options->scan_options.max_curves_per_tile;
    const size_t scratch_limit =
        options->scan_options.max_critical_scratch_bytes;
    size_t cell_capacity = sys->cells.count;
    if (cell_capacity > candidate_limit) cell_capacity = candidate_limit;
    if (cell_capacity > scratch_limit / sizeof(size_t))
        cell_capacity = scratch_limit / sizeof(size_t);
    size_t* cells = cell_capacity
        ? calloc(cell_capacity, sizeof(*cells)) : NULL;
    if (cell_capacity && !cells) return -1;
    *peak_scratch = cell_capacity * sizeof(*cells);
    int status = 0;
    *reason = ALEA_SLICE_ERROR_UNRESOLVED_UNSUPPORTED_GEOMETRY;
    size_t discovery_work = options->scan_options.max_active_boundary_tests;
    const int x_slice_axis = slice_error_slice_axis_for_world(
        &options->view, 0);
    const double x_slice_sign = slice_error_slice_axis_sign(
        &options->view, 0);
    const double world_max_x = options->view.plane.origin[0] +
        (x_slice_sign > 0.0 ? tile->uv_max[x_slice_axis]
                            : -tile->uv_min[x_slice_axis]);
    const double x_margin = 64.0 * DBL_EPSILON *
        fmax(1.0, fabs(world_max_x));
    const double x_cutoff = world_max_x + x_margin;
    for (size_t i = 0; i < query->indexed_count; ++i) {
        const slice_error_index_cell_t* indexed =
            &query->indexed_cells[i];
        if (isfinite(x_cutoff) && indexed->box.min_x > x_cutoff)
            break;
        if (slice_error_box_outside_tile(&indexed->box, options, tile))
            continue;
        if (cell_count == cell_capacity) {
            *reason = ALEA_SLICE_ERROR_UNRESOLVED_CANDIDATE_LIMIT;
            goto selection_done;
        }
        cells[cell_count++] = indexed->cell_index;
    }
    for (size_t i = 0; i < query->uncertain_count; ++i) {
        if (cell_count == cell_capacity) {
            *reason = ALEA_SLICE_ERROR_UNRESOLVED_CANDIDATE_LIMIT;
            goto selection_done;
        }
        cells[cell_count++] = query->uncertain_cells[i];
    }
    for (size_t ci = query->uncached_cell_begin;
         ci < sys->cells.count; ++ci) {
        const alea_cell_entry_t* cell = &sys->cells.data[ci];
        if (cell->universe_id != 0) continue;
        if (cell->original_root_node_id == ALEA_NODE_ID_INVALID &&
            slice_error_bbox_tree_safe(sys, cell->root_node_id, 0,
                                       &discovery_work)) {
            alea_bbox_t box = alea_get_bbox(sys, cell->root_node_id);
            if (box.min_x <= box.max_x && box.min_y <= box.max_y &&
                box.min_z <= box.max_z) {
                box.min_x = slice_error_trusted_lower(box.min_x);
                box.max_x = slice_error_trusted_upper(box.max_x);
                box.min_y = slice_error_trusted_lower(box.min_y);
                box.max_y = slice_error_trusted_upper(box.max_y);
                box.min_z = slice_error_trusted_lower(box.min_z);
                box.max_z = slice_error_trusted_upper(box.max_z);
                if (slice_error_box_outside_tile(&box, options, tile))
                    continue;
            }
        }
        if (cell_count == cell_capacity) {
            *reason = ALEA_SLICE_ERROR_UNRESOLVED_CANDIDATE_LIMIT;
            goto selection_done;
        }
        cells[cell_count++] = ci;
    }
    if (cell_count > 1)
        qsort(cells, cell_count, sizeof(*cells), slice_error_size_compare);
    if (cell_count == 0) {
        const size_t output_limit = options->scan_options.max_output_bytes;
        if (contextual_bytes > output_limit ||
            sizeof(alea_slice_error_region_t) >
                output_limit - contextual_bytes) {
            *reason = ALEA_SLICE_ERROR_UNRESOLVED_CANDIDATE_LIMIT;
            *output_omitted = 1;
            goto selection_done;
        }
        out->regions = calloc(1, sizeof(*out->regions));
        if (!out->regions) { status = -1; goto selection_done; }
        out->region_count = 1;
        memcpy(out->regions[0].uv_min, tile->uv_min,
               sizeof(tile->uv_min));
        memcpy(out->regions[0].uv_max, tile->uv_max,
               sizeof(tile->uv_max));
        out->regions[0].kind = ALEA_POINT_COVERAGE_GAP;
        *reason = ALEA_SLICE_ERROR_RESOLVED;
        status = 1;
        goto selection_done;
    }
    size_t plane_capacity = sys->primitives.count;
    if (plane_capacity > candidate_limit) plane_capacity = candidate_limit;
    const size_t remaining_scratch = scratch_limit -
        cell_capacity * sizeof(*cells);
    const size_t per_plane = sizeof(slice_error_axis_plane_t) +
        2 * sizeof(slice_error_grid_line_t);
    if (remaining_scratch < 4 * sizeof(slice_error_grid_line_t) ||
        plane_capacity > (remaining_scratch -
            4 * sizeof(slice_error_grid_line_t)) / per_plane)
        plane_capacity = remaining_scratch <
            4 * sizeof(slice_error_grid_line_t) ? 0 :
            (remaining_scratch - 4 * sizeof(slice_error_grid_line_t)) /
                per_plane;
    if (!plane_capacity) {
        *reason = ALEA_SLICE_ERROR_UNRESOLVED_CANDIDATE_LIMIT;
        goto selection_done;
    }
    slice_error_axis_plane_t* planes =
        calloc(plane_capacity, sizeof(*planes));
    slice_error_grid_line_t* x =
        calloc(plane_capacity + 2, sizeof(*x));
    slice_error_grid_line_t* y =
        calloc(plane_capacity + 2, sizeof(*y));
    if (!planes || !x || !y) {
        free(planes); free(x); free(y);
        status = -1;
        goto selection_done;
    }
    *peak_scratch = cell_capacity * sizeof(*cells) +
        plane_capacity * sizeof(*planes) +
        2 * (plane_capacity + 2) * sizeof(*x);
    size_t nx = 1, ny = 1;
    for (size_t i = 0; i < cell_count; ++i) {
        const alea_cell_entry_t* cell = &sys->cells.data[cells[i]];
        if (cell->fill_universe > 0 || cell->lat_type != 0 ||
            cell->original_root_node_id != ALEA_NODE_ID_INVALID)
            goto done;
        if (!slice_error_collect_primitives(
                sys, cell->root_node_id, 0, planes, &primitive_count,
                plane_capacity, &discovery_work)) {
            *reason = !discovery_work || primitive_count == plane_capacity
                ? ALEA_SLICE_ERROR_UNRESOLVED_CANDIDATE_LIMIT
                : ALEA_SLICE_ERROR_UNRESOLVED_UNSUPPORTED_GEOMETRY;
            goto done;
        }
    }
    for (size_t si = 0; si < sys->surfaces.count; ++si) {
        const alea_surface_entry_t* surface = &sys->surfaces.data[si];
        const int index = slice_error_plane_index(
            planes, primitive_count, surface->primitive_id);
        if (index < 0) continue;
        slice_error_axis_plane_t* plane = &planes[index];
        if (plane->surface_id != 0) goto done; /* Coincident cards. */
        plane->surface_id = surface->mc_surface_id;
    }
    if (primitive_count == 1) {
        const alea_primitive_entry_t* primitive =
            &sys->primitives.data[planes[0].primitive_id];
        if (primitive->type == ALEA_PRIMITIVE_PLANE &&
            primitive->payload_index < sys->primitive_planes.count) {
            const alea_plane_data_t* p =
                &sys->primitive_planes.data[primitive->payload_index];
            if (p->a != 0.0 && p->b != 0.0 && p->c == 0.0) {
                status = slice_error_classify_single_oblique_tile(
                    query, tile, cells, cell_count, &planes[0],
                    contextual_bytes, &discovery_work, out, reason,
                    output_omitted);
                goto done;
            }
        }
    }
    if (primitive_count == 2) {
        int vertical = 1, oblique = 0;
        for (size_t pi = 0; pi < 2; ++pi) {
            const alea_primitive_entry_t* primitive =
                &sys->primitives.data[planes[pi].primitive_id];
            if (primitive->type != ALEA_PRIMITIVE_PLANE ||
                primitive->payload_index >= sys->primitive_planes.count) {
                vertical = 0;
                break;
            }
            const alea_plane_data_t* p =
                &sys->primitive_planes.data[primitive->payload_index];
            if ((p->a == 0.0 && p->b == 0.0) || p->c != 0.0) {
                vertical = 0;
                break;
            }
            oblique += p->a != 0.0 && p->b != 0.0;
        }
        if (vertical && oblique) {
            if (oblique == 2) {
                status = slice_error_classify_parallel_oblique_tile(
                    query, tile, cells, cell_count, planes,
                    contextual_bytes, &discovery_work, out, reason,
                    output_omitted);
            }
            if (oblique == 1 ||
                (status == 0 && *reason ==
                 ALEA_SLICE_ERROR_UNRESOLVED_UNSUPPORTED_GEOMETRY))
                status = slice_error_classify_crossing_oblique_tile(
                    query, tile, cells, cell_count, planes,
                    contextual_bytes, &discovery_work, out, reason,
                    output_omitted);
            goto done;
        }
    }
    for (size_t pi = 0; pi < primitive_count; ++pi) {
        const uint32_t primitive_id = planes[pi].primitive_id;
        const alea_primitive_entry_t* primitive =
            &sys->primitives.data[primitive_id];
        if (primitive->type != ALEA_PRIMITIVE_PLANE ||
            primitive->payload_index >= sys->primitive_planes.count ||
            planes[pi].surface_id <= 0) goto done;
        const alea_plane_data_t* p =
            &sys->primitive_planes.data[primitive->payload_index];
        if (!isfinite(p->a) || !isfinite(p->b) ||
            !isfinite(p->c) || !isfinite(p->d) || p->c != 0.0)
            goto done;
        const int world_axis = p->a != 0.0 && p->b == 0.0 ? 0
            : p->b != 0.0 && p->a == 0.0 ? 1 : -1;
        if (world_axis < 0) goto done;
        const int axis = slice_error_slice_axis_for_world(
            &options->view, world_axis);
        const double slice_sign = slice_error_slice_axis_sign(
            &options->view, world_axis);
        const double coefficient = world_axis == 0 ? p->a : p->b;
        const long double quotient =
            -(long double)p->d / (long double)coefficient;
        const long double origin =
            (long double)options->view.plane.origin[world_axis];
        const long double coordinate =
            (quotient - origin) / (long double)slice_sign;
        if (!isfinite(coordinate) || fabsl(coordinate) > DBL_MAX) {
            *reason = ALEA_SLICE_ERROR_UNRESOLVED_NUMERICAL;
            goto done;
        }
        const double coord = (double)coordinate;
        const long double uncertainty_ld =
            64.0L * LDBL_EPSILON *
                fmaxl(1.0L, fmaxl(fabsl(quotient), fabsl(origin))) +
            fabsl(coordinate - (long double)coord) +
            8.0L * DBL_EPSILON * fmaxl(1.0L, fabsl(coordinate));
        if (!isfinite(uncertainty_ld) || uncertainty_ld > DBL_MAX) {
            *reason = ALEA_SLICE_ERROR_UNRESOLVED_NUMERICAL;
            goto done;
        }
        const double uncertainty = (double)uncertainty_ld;
        planes[pi] = (slice_error_axis_plane_t){
            .axis = axis,
            .coefficient_sign = coefficient * slice_sign > 0.0 ? 1 : -1,
            .coordinate = coord,
            .uncertainty = uncertainty,
            .surface_id = planes[pi].surface_id,
            .primitive_id = primitive_id
        };
        const double lo = tile->uv_min[axis], hi = tile->uv_max[axis];
        if ((fabsl(coordinate - (long double)lo) <= uncertainty &&
             coordinate != (long double)lo) ||
            (fabsl(coordinate - (long double)hi) <= uncertainty &&
             coordinate != (long double)hi)) {
            *reason = ALEA_SLICE_ERROR_UNRESOLVED_NUMERICAL;
            goto done;
        }
    }
    double support_min[2] = {tile->uv_min[0], tile->uv_min[1]};
    for (int axis = 0; axis < 2; ++axis) {
        const double lo = tile->uv_min[axis];
        if (lo == options->required_uv_min[axis]) continue;
        int seam_plane = 0;
        double previous = options->required_uv_min[axis];
        double previous_uncertainty = 0.0;
        for (size_t pi = 0; pi < primitive_count; ++pi) {
            const slice_error_axis_plane_t* plane = &planes[pi];
            if (plane->axis != axis) continue;
            if (plane->coordinate == lo) seam_plane = 1;
            if (plane->coordinate < lo &&
                plane->coordinate > previous) {
                previous = plane->coordinate;
                previous_uncertainty = plane->uncertainty;
            }
        }
        if (!seam_plane) continue;
        support_min[axis] = previous + (lo - previous) * 0.5;
        if (!(support_min[axis] > previous && support_min[axis] < lo) ||
            support_min[axis] - previous <= previous_uncertainty +
                4.0 * DBL_EPSILON * fmax(1.0, fabs(previous))) {
            *reason = ALEA_SLICE_ERROR_UNRESOLVED_NUMERICAL;
            goto done;
        }
    }
    x[0].coordinate = support_min[0];
    y[0].coordinate = support_min[1];
    for (size_t pi = 0; pi < primitive_count; ++pi) {
        const slice_error_axis_plane_t* plane = &planes[pi];
        const int axis = plane->axis;
        const double coord = plane->coordinate;
        if (coord > support_min[axis] &&
            coord < tile->uv_max[axis]) {
            slice_error_grid_line_t line = {
                .coordinate = coord,
                .uncertainty = plane->uncertainty,
                .surface_id = plane->surface_id,
                .primitive_id = plane->primitive_id
            };
            if (axis == 0) x[nx++] = line;
            else y[ny++] = line;
        }
    }
    x[nx++].coordinate = tile->uv_max[0];
    y[ny++].coordinate = tile->uv_max[1];
    qsort(x, nx, sizeof(*x), slice_error_grid_line_compare);
    qsort(y, ny, sizeof(*y), slice_error_grid_line_compare);
    if (!slice_error_axis_grid_valid(x, nx) ||
        !slice_error_axis_grid_valid(y, ny)) {
        *reason = ALEA_SLICE_ERROR_UNRESOLVED_NUMERICAL;
        goto done;
    }
    const size_t x_faces = nx - 1, y_faces = ny - 1;
    if (x_faces > SIZE_MAX / y_faces ||
        x_faces * y_faces > options->scan_options.max_active_boundary_tests ||
        (nx - 2) > options->scan_options.max_curve_pairs /
            (ny - 2 ? ny - 2 : 1)) {
        *reason = ALEA_SLICE_ERROR_UNRESOLVED_CANDIDATE_LIMIT;
        goto done;
    }
    const size_t face_count = x_faces * y_faces;
    if (face_count > discovery_work / cell_count) {
        *reason = ALEA_SLICE_ERROR_UNRESOLVED_CANDIDATE_LIMIT;
        goto done;
    }
    const size_t scratch_bytes = cell_capacity * sizeof(*cells) +
        plane_capacity * sizeof(*planes) +
        2 * (plane_capacity + 2) * sizeof(*x);
    if (face_count > SIZE_MAX / sizeof(slice_error_face_t) ||
        scratch_bytes > options->scan_options.max_critical_scratch_bytes ||
        face_count * sizeof(slice_error_face_t) >
            options->scan_options.max_critical_scratch_bytes - scratch_bytes) {
        *reason = ALEA_SLICE_ERROR_UNRESOLVED_CANDIDATE_LIMIT;
        goto done;
    }
    *peak_scratch = scratch_bytes +
        face_count * sizeof(slice_error_face_t);
    slice_error_face_t* faces = calloc(face_count, sizeof(*faces));
    if (!faces) { status = -1; goto done; }
    size_t work_remaining = discovery_work;
    for (size_t yi = 0; yi < y_faces; ++yi) {
        for (size_t xi = 0; xi < x_faces; ++xi) {
            slice_error_face_t* face = &faces[yi * x_faces + xi];
            const double u = x[xi].coordinate +
                0.5 * (x[xi + 1].coordinate - x[xi].coordinate);
            const double v = y[yi].coordinate +
                0.5 * (y[yi + 1].coordinate - y[yi].coordinate);
            const double wx = options->view.plane.origin[0] +
                options->view.plane.u_axis[0] * u +
                options->view.plane.v_axis[0] * v;
            const double wy = options->view.plane.origin[1] +
                options->view.plane.u_axis[1] * u +
                options->view.plane.v_axis[1] * v;
            if (!(u > x[xi].coordinate && u < x[xi + 1].coordinate) ||
                !(v > y[yi].coordinate && v < y[yi + 1].coordinate) ||
                !isfinite(wx) || !isfinite(wy)) {
                *reason = ALEA_SLICE_ERROR_UNRESOLVED_NUMERICAL;
                free(faces); goto done;
            }
            for (size_t ci = 0; ci < cell_count; ++ci) {
                const alea_cell_entry_t* cell = &sys->cells.data[cells[ci]];
                int inside = 0;
                if (!slice_error_axis_node_inside(
                        sys, cell->root_node_id, planes, primitive_count,
                        x, xi, y, yi, 0, &work_remaining, &inside)) {
                    *reason = work_remaining == 0
                        ? ALEA_SLICE_ERROR_UNRESOLVED_CANDIDATE_LIMIT
                        : ALEA_SLICE_ERROR_UNRESOLVED_NUMERICAL;
                    free(faces); goto done;
                }
                if (inside != alea_point_inside(
                        sys, cell->root_node_id, wx, wy,
                        options->view.plane.origin[2])) {
                    *reason = ALEA_SLICE_ERROR_UNRESOLVED_NUMERICAL;
                    free(faces); goto done;
                }
                if (inside) {
                    if (face->owner_count == ALEA_SLICE_ERROR_OWNER_CAPACITY) {
                        *reason = ALEA_SLICE_ERROR_UNRESOLVED_CANDIDATE_LIMIT;
                        free(faces); goto done;
                    }
                    face->owner_cell_ids[face->owner_count++] =
                        cell->mc_cell_id;
                }
            }
            face->kind = face->owner_count == 0
                ? ALEA_POINT_COVERAGE_GAP
                : face->owner_count == 1
                    ? ALEA_POINT_COVERAGE_UNIQUE
                    : ALEA_POINT_COVERAGE_OVERLAP;
        }
    }
    size_t region_count = 0, interval_count = 0;
    for (size_t yi = 0; yi < y_faces; ++yi)
        for (size_t xi = 0; xi < x_faces; ++xi)
            if (x[xi + 1].coordinate > tile->uv_min[0] &&
                x[xi].coordinate < tile->uv_max[0] &&
                y[yi + 1].coordinate > tile->uv_min[1] &&
                y[yi].coordinate < tile->uv_max[1])
                region_count += slice_error_axis_is_defect(
                    faces[yi * x_faces + xi].kind);
    for (size_t yi = 0; yi < y_faces; ++yi)
        for (size_t xi = 1; xi < x_faces; ++xi) {
            if (x[xi].coordinate < tile->uv_min[0] ||
                x[xi].coordinate >= tile->uv_max[0] ||
                y[yi + 1].coordinate <= tile->uv_min[1] ||
                y[yi].coordinate >= tile->uv_max[1]) continue;
            const slice_error_face_t* a = &faces[yi * x_faces + xi - 1];
            const slice_error_face_t* b = &faces[yi * x_faces + xi];
            interval_count += slice_error_axis_is_defect(a->kind) !=
                              slice_error_axis_is_defect(b->kind);
        }
    for (size_t yi = 1; yi < y_faces; ++yi)
        for (size_t xi = 0; xi < x_faces; ++xi) {
            if (y[yi].coordinate < tile->uv_min[1] ||
                y[yi].coordinate >= tile->uv_max[1] ||
                x[xi + 1].coordinate <= tile->uv_min[0] ||
                x[xi].coordinate >= tile->uv_max[0]) continue;
            const slice_error_face_t* a = &faces[(yi - 1) * x_faces + xi];
            const slice_error_face_t* b = &faces[yi * x_faces + xi];
            interval_count += slice_error_axis_is_defect(a->kind) !=
                              slice_error_axis_is_defect(b->kind);
        }
    const size_t output_limit = options->scan_options.max_output_bytes;
    if (contextual_bytes > output_limit ||
        region_count > (output_limit - contextual_bytes) /
            sizeof(alea_slice_error_region_t) ||
        interval_count > (output_limit - contextual_bytes -
            region_count * sizeof(alea_slice_error_region_t)) /
            sizeof(alea_slice_error_interval_t)) {
        *reason = ALEA_SLICE_ERROR_UNRESOLVED_CANDIDATE_LIMIT;
        *output_omitted = 1;
        free(faces); goto done;
    }
    out->regions = region_count
        ? calloc(region_count, sizeof(*out->regions)) : NULL;
    out->intervals = interval_count
        ? calloc(interval_count, sizeof(*out->intervals)) : NULL;
    if ((region_count && !out->regions) ||
        (interval_count && !out->intervals)) {
        free(faces); status = -1; goto done;
    }
    for (size_t yi = 0; yi < y_faces; ++yi)
        for (size_t xi = 0; xi < x_faces; ++xi) {
            const slice_error_face_t* face = &faces[yi * x_faces + xi];
            if (!slice_error_axis_is_defect(face->kind) ||
                x[xi + 1].coordinate <= tile->uv_min[0] ||
                x[xi].coordinate >= tile->uv_max[0] ||
                y[yi + 1].coordinate <= tile->uv_min[1] ||
                y[yi].coordinate >= tile->uv_max[1]) continue;
            alea_slice_error_region_t* region =
                &out->regions[out->region_count++];
            region->uv_min[0] = fmax(x[xi].coordinate, tile->uv_min[0]);
            region->uv_max[0] = fmin(x[xi + 1].coordinate, tile->uv_max[0]);
            region->uv_min[1] = fmax(y[yi].coordinate, tile->uv_min[1]);
            region->uv_max[1] = fmin(y[yi + 1].coordinate, tile->uv_max[1]);
            region->uv_min_uncertainty[0] =
                region->uv_min[0] == x[xi].coordinate
                    ? x[xi].uncertainty : 0.0;
            region->uv_max_uncertainty[0] =
                region->uv_max[0] == x[xi + 1].coordinate
                    ? x[xi + 1].uncertainty : 0.0;
            region->uv_min_uncertainty[1] =
                region->uv_min[1] == y[yi].coordinate
                    ? y[yi].uncertainty : 0.0;
            region->uv_max_uncertainty[1] =
                region->uv_max[1] == y[yi + 1].coordinate
                    ? y[yi + 1].uncertainty : 0.0;
            region->kind = face->kind;
            region->owner_count = face->owner_count;
            memcpy(region->owner_cell_ids, face->owner_cell_ids,
                   face->owner_count * sizeof(int));
        }
    for (int axis = 0; axis < 2; ++axis) {
        const size_t major_faces = axis == 0 ? x_faces : y_faces;
        const size_t minor_faces = axis == 0 ? y_faces : x_faces;
        const slice_error_grid_line_t* major = axis == 0 ? x : y;
        const slice_error_grid_line_t* minor = axis == 0 ? y : x;
        for (size_t mi = 1; mi < major_faces; ++mi)
            for (size_t ni = 0; ni < minor_faces; ++ni) {
                if (major[mi].coordinate < tile->uv_min[axis] ||
                    major[mi].coordinate >= tile->uv_max[axis] ||
                    minor[ni + 1].coordinate <= tile->uv_min[1 - axis] ||
                    minor[ni].coordinate >= tile->uv_max[1 - axis]) continue;
                const size_t negative_index = axis == 0
                    ? ni * x_faces + mi - 1
                    : (mi - 1) * x_faces + ni;
                const size_t positive_index = axis == 0
                    ? ni * x_faces + mi
                    : mi * x_faces + ni;
                const slice_error_face_t* negative =
                    &faces[negative_index];
                const slice_error_face_t* positive =
                    &faces[positive_index];
                if (slice_error_axis_is_defect(negative->kind) ==
                    slice_error_axis_is_defect(positive->kind)) continue;
                alea_slice_error_interval_t* interval =
                    &out->intervals[out->interval_count++];
                interval->evidence_scope =
                    ALEA_SLICE_BOUNDARY_EVIDENCE_VERIFIED_INTERVAL;
                interval->surface_id = major[mi].surface_id;
                interval->primitive_id = major[mi].primitive_id;
                interval->axis = axis;
                interval->uv_start[axis] = major[mi].coordinate;
                interval->uv_end[axis] = major[mi].coordinate;
                interval->uv_start[1 - axis] = fmax(
                    minor[ni].coordinate, tile->uv_min[1 - axis]);
                interval->uv_end[1 - axis] = fmin(
                    minor[ni + 1].coordinate, tile->uv_max[1 - axis]);
                interval->endpoint_uncertainty[0] = fmax(
                    major[mi].uncertainty, minor[ni].uncertainty);
                interval->endpoint_uncertainty[1] = fmax(
                    major[mi].uncertainty, minor[ni + 1].uncertainty);
                interval->negative_side_kind = negative->kind;
                interval->positive_side_kind = positive->kind;
                interval->negative_owner_count = negative->owner_count;
                interval->positive_owner_count = positive->owner_count;
                memcpy(interval->negative_owner_cell_ids,
                       negative->owner_cell_ids,
                       negative->owner_count * sizeof(int));
                memcpy(interval->positive_owner_cell_ids,
                       positive->owner_cell_ids,
                       positive->owner_count * sizeof(int));
            }
    }
    free(faces);
    *reason = ALEA_SLICE_ERROR_RESOLVED;
    status = 1;
done:
    free(planes); free(x); free(y);
    if (status != 1) {
        free(out->intervals); out->intervals = NULL;
        free(out->regions); out->regions = NULL;
        out->interval_count = out->region_count = 0;
    }
selection_done:
    free(cells);
    return status;
}

int alea_slice_error_query_run_page(alea_slice_error_query_t* query,
                                    size_t page_index,
                                    alea_slice_error_page_t* page) {
    if (!query || !page || page_index >= query->page_count ||
        alea_system_geometry_generation(query->sys) !=
            query->geometry_generation) return -1;
    const alea_slice_error_query_options_t* options = &query->options;
    const size_t col = page_index % options->tile_columns;
    const size_t row = page_index / options->tile_columns;
    alea_transition_slice_critical_tile_t tile = {0};
    for (size_t axis = 0; axis < 2; ++axis) {
        const size_t index = axis == 0 ? col : row;
        const size_t count = axis == 0
            ? options->tile_columns : options->tile_rows;
        const double lo = options->required_uv_min[axis];
        const double hi = options->required_uv_max[axis];
        const double width = hi - lo;
        if (!isfinite(width)) return -1;
        tile.uv_min[axis] = index == 0 ? lo : lo + width *
            ((double)index / (double)count);
        tile.uv_max[axis] = index + 1 == count ? hi : lo + width *
            ((double)(index + 1) / (double)count);
        if (!(tile.uv_max[axis] > tile.uv_min[axis])) return -1;
    }
    if (alea_raycast_ensure_hier_caches(query->sys) != 0) return -1;
    alea_transition_slice_stats_t stats = {0};
    stats.critical_stop_reason = ALEA_TRANSITION_SLICE_CRITICAL_NONE;
    alea_slice_error_page_t candidate = {0};
    slice_error_finding_sink_t sink = {
        .page = &candidate,
        .max_findings = options->scan_options.max_critical_findings,
        .max_output_bytes = options->scan_options.max_output_bytes
    };
    alea_slice_view_t analysis_view = options->view;
    analysis_view.u_min = options->required_uv_min[0];
    analysis_view.u_max = options->required_uv_max[0];
    analysis_view.v_min = options->required_uv_min[1];
    analysis_view.v_max = options->required_uv_max[1];
    if (alea_transition_slice_enumerate_critical_tiles(
            query->sys, &analysis_view, &options->scan_options,
            &tile, 1, slice_error_retain_finding, &sink, &stats) != 0 ||
        alea_system_geometry_generation(query->sys) !=
            query->geometry_generation) {
        free(candidate.findings);
        return -1;
    }
    alea_slice_error_page_receipt_t receipt = {0};
    receipt.query_id = query->query_id;
    receipt.geometry_generation = query->geometry_generation;
    receipt.page_index = page_index;
    memcpy(receipt.core_uv_min, tile.uv_min, sizeof(tile.uv_min));
    memcpy(receipt.core_uv_max, tile.uv_max, sizeof(tile.uv_max));
    receipt.occurrence_paths = stats.critical_occurrence_paths;
    receipt.candidate_curves = stats.critical_curves;
    receipt.candidate_pairs_tested = stats.critical_curve_pairs_tested;
    receipt.peak_scratch_bytes = stats.peak_critical_scratch_bytes;
    receipt.query_index_bytes = query->index_bytes;
    receipt.contextual_finding_count = candidate.finding_count;
    receipt.omitted_contextual_findings = sink.omitted_findings;
    receipt.omitted_contextual_boundary_evidence =
        stats.omitted_critical_boundary_evidence;
    receipt.scan_stop_reason = stats.critical_stop_reason;
    receipt.unresolved_reason = stats.critical_stop_reason ==
            ALEA_TRANSITION_SLICE_CRITICAL_NONE
        ? ALEA_SLICE_ERROR_UNRESOLVED_CLASSIFIER_PENDING
        : stats.critical_stop_reason ==
            ALEA_TRANSITION_SLICE_CRITICAL_UNSUPPORTED_CURVE ||
          stats.critical_stop_reason ==
            ALEA_TRANSITION_SLICE_CRITICAL_UNSUPPORTED_OCCURRENCE_TRAVERSAL
        ? ALEA_SLICE_ERROR_UNRESOLVED_UNSUPPORTED_GEOMETRY
        : ALEA_SLICE_ERROR_UNRESOLVED_CANDIDATE_LIMIT;
    receipt.requested_work_complete = stats.critical_stop_reason ==
        ALEA_TRANSITION_SLICE_CRITICAL_NONE && !sink.omitted_findings;
    receipt.scope_classified = 0;
    receipt.output_complete = sink.omitted_findings == 0 &&
        stats.omitted_critical_boundary_evidence == 0;
    if (receipt.requested_work_complete) {
        int verified_output_omitted = 0;
        size_t classifier_peak_scratch = 0;
        const int verified = slice_error_classify_axis_tile(
            query, &tile,
            candidate.finding_count * sizeof(*candidate.findings),
            &candidate, &receipt.unresolved_reason,
            &verified_output_omitted, &classifier_peak_scratch);
        if (receipt.peak_scratch_bytes < classifier_peak_scratch)
            receipt.peak_scratch_bytes = classifier_peak_scratch;
        if (verified < 0) {
            free(candidate.findings);
            free(candidate.intervals);
            free(candidate.regions);
            return -1;
        }
        if (verified > 0) receipt.scope_classified = 1;
        if (verified_output_omitted) receipt.output_complete = 0;
    }
    receipt.verified_interval_count = candidate.interval_count;
    receipt.region_count = candidate.region_count;
    candidate.receipt = receipt;
    candidate.populated = 1;
    free(page->findings);
    free(page->intervals);
    free(page->regions);
    *page = candidate;
    return 0;
}

struct alea_transition_slice_result {
    alea_transition_slice_finding_t* findings;
    size_t finding_count;
    size_t finding_capacity;
    alea_transition_slice_component_t* components;
    size_t component_count;
    size_t component_capacity;
    alea_transition_slice_coverage_finding_t* coverage_findings;
    size_t coverage_finding_count;
    size_t coverage_finding_capacity;
    alea_transition_slice_coverage_component_t* coverage_components;
    size_t coverage_component_count;
    size_t coverage_component_capacity;
    alea_transition_slice_component_link_t* component_links;
    size_t component_link_count;
    size_t component_link_capacity;
    alea_transition_slice_refinement_frontier_t* refinement_frontiers;
    size_t refinement_frontier_count;
    size_t refinement_frontier_capacity;
    alea_transition_slice_critical_tile_t* critical_tiles;
    size_t critical_tile_count;
    alea_transition_slice_critical_tile_source_t* critical_tile_sources;
    size_t critical_tile_source_count;
    alea_transition_slice_critical_finding_t* critical_findings;
    size_t critical_finding_count;
    size_t critical_finding_capacity;
    alea_transition_slice_stats_t stats;
};

typedef struct {
    double transverse_coordinate;
    size_t base_ray_index;
    uint32_t refinement_depth;
    uint64_t signature_a;
    uint64_t signature_b;
    size_t event_count;
    size_t finding_count;
} transition_slice_row_t;

typedef struct {
    alea_cell_hit_t* hits;
    uint64_t* occurrence_keys;
    uint64_t* parent_occurrence_keys;
    uint8_t* owner_mask;
    size_t capacity;
} transition_slice_coverage_scratch_t;

void alea_transition_slice_options_init(
    alea_transition_slice_options_t* options) {
    if (!options) return;
    memset(options, 0, sizeof(*options));
    options->struct_size = sizeof(*options);
    options->horizontal_rays = 64;
    options->vertical_rays = 64;
    options->max_rays = 4096;
    options->max_events = 1000000;
    options->max_events_per_ray = 8192;
    options->max_findings = 10000;
    options->max_components = 10000;
    options->max_output_bytes = 64u * 1024u * 1024u;
    options->max_scratch_bytes = 16u * 1024u * 1024u;
    options->max_coverage_fallbacks = 10000;
    options->max_coverage_hits = 256;
    options->refine_signals = ALEA_TRANSITION_SLICE_REFINE_SIGNATURE;
    options->max_row_scratch_bytes = 1024u * 1024u;
    options->max_coverage_probes = 10000;
    options->max_coverage_findings = 10000;
    options->max_coverage_components = 10000;
    options->max_component_links = 10000;
    options->max_refinement_frontiers = 1024;
    options->max_critical_tiles = 256;
    options->max_critical_tile_sources = 1024;
    options->max_critical_scratch_bytes = 4u * 1024u * 1024u;
    options->max_curves_per_tile = 1024;
    options->max_critical_points = 2048;
    options->max_active_boundary_tests = 100000;
    options->max_critical_probes = 4096;
    options->max_critical_findings = 1024;
    options->max_critical_boundary_evidence = 1024;
    options->max_curve_pairs = 100000;
    options->max_critical_sector_witnesses = 8192;
    options->occurrence_discovery = ALEA_TRANSITION_SLICE_OCCURRENCE_SAMPLED;
    options->max_exhaustive_occurrence_hits = 256;
}

const char* alea_transition_slice_critical_stop_reason_name(
    alea_transition_slice_critical_stop_reason_t reason) {
    switch (reason) {
    case ALEA_TRANSITION_SLICE_CRITICAL_DISABLED: return "disabled";
    case ALEA_TRANSITION_SLICE_CRITICAL_NONE: return "none";
    case ALEA_TRANSITION_SLICE_CRITICAL_MAX_FRONTIERS:
        return "max_refinement_frontiers";
    case ALEA_TRANSITION_SLICE_CRITICAL_MAX_TILE_SOURCES:
        return "max_critical_tile_sources";
    case ALEA_TRANSITION_SLICE_CRITICAL_MAX_TILES:
        return "max_critical_tiles";
    case ALEA_TRANSITION_SLICE_CRITICAL_MAX_SCRATCH_BYTES:
        return "max_critical_scratch_bytes";
    case ALEA_TRANSITION_SLICE_CRITICAL_MAX_OUTPUT_BYTES:
        return "max_output_bytes";
    case ALEA_TRANSITION_SLICE_CRITICAL_MAX_CURVES:
        return "max_curves_per_tile";
    case ALEA_TRANSITION_SLICE_CRITICAL_CHAIN_TRUNCATED:
        return "occurrence_chain_truncated";
    case ALEA_TRANSITION_SLICE_CRITICAL_MAX_POINTS:
        return "max_critical_points";
    case ALEA_TRANSITION_SLICE_CRITICAL_MAX_PROBES:
        return "max_critical_probes";
    case ALEA_TRANSITION_SLICE_CRITICAL_MAX_FINDINGS:
        return "max_critical_findings";
    case ALEA_TRANSITION_SLICE_CRITICAL_UNSUPPORTED_CURVE:
        return "unsupported_curve";
    case ALEA_TRANSITION_SLICE_CRITICAL_MAX_CURVE_PAIRS:
        return "max_curve_pairs";
    case ALEA_TRANSITION_SLICE_CRITICAL_MAX_SECTOR_WITNESSES:
        return "max_critical_sector_witnesses";
    case ALEA_TRANSITION_SLICE_CRITICAL_MAX_OCCURRENCE_HITS:
        return "max_exhaustive_occurrence_hits";
    case ALEA_TRANSITION_SLICE_CRITICAL_UNSUPPORTED_OCCURRENCE_TRAVERSAL:
        return "unsupported_occurrence_traversal";
    case ALEA_TRANSITION_SLICE_CRITICAL_NUMERICAL_UNRESOLVED:
        return "numerical_unresolved";
    }
    return "unknown";
}

const char* alea_transition_slice_refinement_status_name(
    alea_transition_slice_refinement_status_t status) {
    switch (status) {
    case ALEA_TRANSITION_SLICE_REFINEMENT_NOT_REQUESTED:
        return "not_requested";
    case ALEA_TRANSITION_SLICE_REFINEMENT_CONVERGED: return "converged";
    case ALEA_TRANSITION_SLICE_REFINEMENT_MAX_DEPTH: return "max_depth";
    case ALEA_TRANSITION_SLICE_REFINEMENT_MIN_SPACING: return "min_spacing";
    case ALEA_TRANSITION_SLICE_REFINEMENT_STOPPED: return "stopped";
    }
    return "unknown";
}

const char* alea_transition_slice_stop_reason_name(
    alea_transition_slice_stop_reason_t reason) {
    switch (reason) {
    case ALEA_TRANSITION_SLICE_STOP_NONE: return "none";
    case ALEA_TRANSITION_SLICE_STOP_MAX_RAYS: return "max_rays";
    case ALEA_TRANSITION_SLICE_STOP_MAX_EVENTS: return "max_events";
    case ALEA_TRANSITION_SLICE_STOP_MAX_FINDINGS: return "max_findings";
    case ALEA_TRANSITION_SLICE_STOP_MAX_COMPONENTS: return "max_components";
    case ALEA_TRANSITION_SLICE_STOP_MAX_OUTPUT_BYTES: return "max_output_bytes";
    case ALEA_TRANSITION_SLICE_STOP_MAX_COVERAGE_FALLBACKS:
        return "max_coverage_fallbacks";
    case ALEA_TRANSITION_SLICE_STOP_MAX_COVERAGE_PROBES:
        return "max_coverage_probes";
    case ALEA_TRANSITION_SLICE_STOP_MAX_SCRATCH_BYTES:
        return "max_scratch_bytes";
    case ALEA_TRANSITION_SLICE_STOP_INTERRUPTED: return "interrupted";
    case ALEA_TRANSITION_SLICE_STOP_MAX_COMPONENT_LINKS:
        return "max_component_links";
    }
    return "unknown";
}

alea_transition_slice_result_t* alea_transition_slice_result_create(void) {
    return calloc(1, sizeof(alea_transition_slice_result_t));
}

void alea_transition_slice_result_destroy(
    alea_transition_slice_result_t* result) {
    if (!result) return;
    free(result->findings);
    free(result->components);
    free(result->coverage_findings);
    free(result->coverage_components);
    free(result->component_links);
    free(result->refinement_frontiers);
    free(result->critical_tiles);
    free(result->critical_tile_sources);
    free(result->critical_findings);
    free(result);
}

size_t alea_transition_slice_component_count(
    const alea_transition_slice_result_t* result) {
    return result ? result->component_count : 0;
}

int alea_transition_slice_component_get(
    const alea_transition_slice_result_t* result, size_t index,
    alea_transition_slice_component_t* out_component) {
    if (!result || !out_component || index >= result->component_count)
        return -1;
    *out_component = result->components[index];
    return 0;
}

size_t alea_transition_slice_coverage_finding_count(
    const alea_transition_slice_result_t* result) {
    return result ? result->coverage_finding_count : 0;
}

int alea_transition_slice_coverage_finding_get(
    const alea_transition_slice_result_t* result, size_t index,
    alea_transition_slice_coverage_finding_t* out_finding) {
    if (!result || !out_finding || index >= result->coverage_finding_count)
        return -1;
    *out_finding = result->coverage_findings[index];
    return 0;
}

size_t alea_transition_slice_coverage_component_count(
    const alea_transition_slice_result_t* result) {
    return result ? result->coverage_component_count : 0;
}

int alea_transition_slice_coverage_component_get(
    const alea_transition_slice_result_t* result, size_t index,
    alea_transition_slice_coverage_component_t* out_component) {
    if (!result || !out_component || index >= result->coverage_component_count)
        return -1;
    *out_component = result->coverage_components[index];
    return 0;
}

size_t alea_transition_slice_component_link_count(
    const alea_transition_slice_result_t* result) {
    return result ? result->component_link_count : 0;
}

int alea_transition_slice_component_link_get(
    const alea_transition_slice_result_t* result, size_t index,
    alea_transition_slice_component_link_t* out_link) {
    if (!result || !out_link || index >= result->component_link_count)
        return -1;
    *out_link = result->component_links[index];
    return 0;
}

size_t alea_transition_slice_refinement_frontier_count(
    const alea_transition_slice_result_t* result) {
    return result ? result->refinement_frontier_count : 0;
}

int alea_transition_slice_refinement_frontier_get(
    const alea_transition_slice_result_t* result, size_t index,
    alea_transition_slice_refinement_frontier_t* out_frontier) {
    if (!result || !out_frontier ||
        index >= result->refinement_frontier_count) return -1;
    *out_frontier = result->refinement_frontiers[index];
    return 0;
}

size_t alea_transition_slice_critical_tile_count(
    const alea_transition_slice_result_t* result) {
    return result ? result->critical_tile_count : 0;
}

int alea_transition_slice_critical_tile_get(
    const alea_transition_slice_result_t* result, size_t index,
    alea_transition_slice_critical_tile_t* out_tile) {
    if (!result || !out_tile || index >= result->critical_tile_count) return -1;
    *out_tile = result->critical_tiles[index];
    return 0;
}

size_t alea_transition_slice_critical_tile_source_count(
    const alea_transition_slice_result_t* result) {
    return result ? result->critical_tile_source_count : 0;
}

int alea_transition_slice_critical_tile_source_get(
    const alea_transition_slice_result_t* result, size_t index,
    alea_transition_slice_critical_tile_source_t* out_source) {
    if (!result || !out_source ||
        index >= result->critical_tile_source_count) return -1;
    *out_source = result->critical_tile_sources[index];
    return 0;
}

size_t alea_transition_slice_critical_finding_count(
    const alea_transition_slice_result_t* result) {
    return result ? result->critical_finding_count : 0;
}

int alea_transition_slice_critical_finding_get(
    const alea_transition_slice_result_t* result, size_t index,
    alea_transition_slice_critical_finding_t* out_finding) {
    if (!result || !out_finding || index >= result->critical_finding_count)
        return -1;
    *out_finding = result->critical_findings[index];
    return 0;
}

size_t alea_transition_slice_finding_count(
    const alea_transition_slice_result_t* result) {
    return result ? result->finding_count : 0;
}

int alea_transition_slice_finding_get(
    const alea_transition_slice_result_t* result, size_t index,
    alea_transition_slice_finding_t* out_finding) {
    if (!result || !out_finding || index >= result->finding_count) return -1;
    *out_finding = result->findings[index];
    return 0;
}

int alea_transition_slice_stats(
    const alea_transition_slice_result_t* result,
    alea_transition_slice_stats_t* out_stats) {
    if (!result || !out_stats) return -1;
    *out_stats = result->stats;
    return 0;
}

static int transition_slice_add_bytes(size_t* total, size_t count,
                                      size_t item_size) {
    if (count > SIZE_MAX / item_size ||
        *total > SIZE_MAX - count * item_size)
        return -1;
    *total += count * item_size;
    return 0;
}

static size_t transition_slice_retained_bytes(
    const alea_transition_slice_result_t* result) {
    size_t bytes = 0;
    if (transition_slice_add_bytes(&bytes, result->finding_count,
                                   sizeof(*result->findings)) ||
        transition_slice_add_bytes(&bytes, result->component_count,
                                   sizeof(*result->components)) ||
        transition_slice_add_bytes(&bytes, result->coverage_finding_count,
                                   sizeof(*result->coverage_findings)) ||
        transition_slice_add_bytes(&bytes, result->coverage_component_count,
                                   sizeof(*result->coverage_components)) ||
        transition_slice_add_bytes(&bytes, result->component_link_count,
                                   sizeof(*result->component_links)) ||
        transition_slice_add_bytes(&bytes, result->refinement_frontier_count,
                                   sizeof(*result->refinement_frontiers)) ||
        transition_slice_add_bytes(&bytes, result->critical_tile_count,
                                   sizeof(*result->critical_tiles)) ||
        transition_slice_add_bytes(&bytes, result->critical_tile_source_count,
                                   sizeof(*result->critical_tile_sources)) ||
        transition_slice_add_bytes(&bytes, result->critical_finding_count,
                                   sizeof(*result->critical_findings)))
        return SIZE_MAX;
    return bytes;
}

typedef struct {
    alea_transition_slice_result_t* result;
    const alea_transition_slice_options_t* options;
} transition_slice_critical_sink_context_t;

static int transition_slice_append_critical_finding(
    const alea_transition_slice_critical_finding_t* finding,
    void* userdata) {
    transition_slice_critical_sink_context_t* context = userdata;
    alea_transition_slice_result_t* result = context->result;
    const alea_transition_slice_options_t* options = context->options;
    if (options->max_critical_findings &&
        result->critical_finding_count >= options->max_critical_findings) {
        result->stats.omitted_critical_findings++;
        return 1;
    }
    const size_t retained = transition_slice_retained_bytes(result);
    if (options->max_output_bytes &&
        (retained > options->max_output_bytes ||
         sizeof(*result->critical_findings) >
            options->max_output_bytes - retained)) {
        result->stats.omitted_critical_findings++;
        return 1;
    }
    if (result->critical_finding_count == result->critical_finding_capacity) {
        size_t capacity = result->critical_finding_capacity
            ? result->critical_finding_capacity * 2u : 16u;
        if (options->max_critical_findings &&
            capacity > options->max_critical_findings)
            capacity = (size_t)options->max_critical_findings;
        void* memory = realloc(
            result->critical_findings,
            capacity * sizeof(*result->critical_findings));
        if (!memory) return -1;
        result->critical_findings = memory;
        result->critical_finding_capacity = capacity;
    }
    result->critical_findings[result->critical_finding_count++] = *finding;
    result->stats.critical_findings = result->critical_finding_count;
    result->stats.retained_output_bytes = transition_slice_retained_bytes(result);
    return 0;
}

static size_t transition_slice_coverage_scratch_bytes(size_t capacity) {
    size_t bytes = 0;
    if (transition_slice_add_bytes(&bytes, capacity, sizeof(alea_cell_hit_t)) ||
        transition_slice_add_bytes(&bytes, capacity, sizeof(uint64_t)) ||
        transition_slice_add_bytes(&bytes, capacity, sizeof(uint64_t)) ||
        transition_slice_add_bytes(&bytes, capacity, sizeof(uint8_t)))
        return SIZE_MAX;
    return bytes;
}

static int transition_slice_coverage_scratch_init(
    transition_slice_coverage_scratch_t* scratch, size_t capacity) {
    memset(scratch, 0, sizeof(*scratch));
    if (capacity == 0) return 0;
    scratch->hits = calloc(capacity, sizeof(*scratch->hits));
    scratch->occurrence_keys = calloc(capacity, sizeof(*scratch->occurrence_keys));
    scratch->parent_occurrence_keys =
        calloc(capacity, sizeof(*scratch->parent_occurrence_keys));
    scratch->owner_mask = calloc(capacity, sizeof(*scratch->owner_mask));
    if (!scratch->hits || !scratch->occurrence_keys ||
        !scratch->parent_occurrence_keys || !scratch->owner_mask) {
        free(scratch->hits); free(scratch->occurrence_keys);
        free(scratch->parent_occurrence_keys); free(scratch->owner_mask);
        memset(scratch, 0, sizeof(*scratch));
        return -1;
    }
    scratch->capacity = capacity;
    return 0;
}

static void transition_slice_coverage_scratch_free(
    transition_slice_coverage_scratch_t* scratch) {
    free(scratch->hits); free(scratch->occurrence_keys);
    free(scratch->parent_occurrence_keys); free(scratch->owner_mask);
    memset(scratch, 0, sizeof(*scratch));
}

static int transition_slice_append(
    alea_transition_slice_result_t* result,
    const alea_transition_slice_options_t* options,
    const alea_transition_slice_finding_t* finding) {
    if (options->max_findings &&
        result->finding_count >= options->max_findings) {
        result->stats.stop_reason = ALEA_TRANSITION_SLICE_STOP_MAX_FINDINGS;
        return 1;
    }
    size_t retained = transition_slice_retained_bytes(result);
    if (options->max_output_bytes &&
        (retained > options->max_output_bytes ||
         sizeof(*result->findings) > options->max_output_bytes - retained)) {
        result->stats.stop_reason =
            ALEA_TRANSITION_SLICE_STOP_MAX_OUTPUT_BYTES;
        return 1;
    }
    if (result->finding_count == result->finding_capacity) {
        size_t capacity = result->finding_capacity
            ? result->finding_capacity * 2 : 16;
        size_t hard_capacity = SIZE_MAX / sizeof(*result->findings);
        if (options->max_findings && hard_capacity > options->max_findings)
            hard_capacity = (size_t)options->max_findings;
        if (options->max_output_bytes) {
            size_t byte_capacity = (size_t)(options->max_output_bytes /
                                             sizeof(*result->findings));
            if (hard_capacity > byte_capacity) hard_capacity = byte_capacity;
        }
        if (capacity > hard_capacity) capacity = hard_capacity;
        if (capacity <= result->finding_capacity) {
            result->stats.stop_reason = options->max_findings &&
                result->finding_capacity >= options->max_findings
                ? ALEA_TRANSITION_SLICE_STOP_MAX_FINDINGS
                : ALEA_TRANSITION_SLICE_STOP_MAX_OUTPUT_BYTES;
            return 1;
        }
        void* memory = realloc(
            result->findings, capacity * sizeof(*result->findings));
        if (!memory) return -1;
        result->findings = memory;
        result->finding_capacity = capacity;
    }
    result->findings[result->finding_count++] = *finding;
    result->stats.findings = result->finding_count;
    result->stats.retained_output_bytes = transition_slice_retained_bytes(result);
    return 0;
}

static int transition_slice_coverage_append(
    alea_transition_slice_result_t* result,
    const alea_transition_slice_options_t* options,
    const alea_transition_slice_coverage_finding_t* finding) {
    if (options->max_coverage_findings &&
        result->coverage_finding_count >= options->max_coverage_findings) {
        result->stats.stop_reason = ALEA_TRANSITION_SLICE_STOP_MAX_FINDINGS;
        return 1;
    }
    size_t retained = transition_slice_retained_bytes(result);
    if (options->max_output_bytes &&
        (retained > options->max_output_bytes ||
         sizeof(*result->coverage_findings) >
            options->max_output_bytes - retained)) {
        result->stats.stop_reason = ALEA_TRANSITION_SLICE_STOP_MAX_OUTPUT_BYTES;
        return 1;
    }
    if (result->coverage_finding_count == result->coverage_finding_capacity) {
        size_t capacity = result->coverage_finding_capacity
            ? result->coverage_finding_capacity * 2 : 16;
        size_t hard_capacity = SIZE_MAX / sizeof(*result->coverage_findings);
        if (options->max_coverage_findings &&
            hard_capacity > options->max_coverage_findings)
            hard_capacity = (size_t)options->max_coverage_findings;
        if (capacity > hard_capacity) capacity = hard_capacity;
        if (capacity <= result->coverage_finding_capacity) {
            result->stats.stop_reason = ALEA_TRANSITION_SLICE_STOP_MAX_FINDINGS;
            return 1;
        }
        void* memory = realloc(
            result->coverage_findings,
            capacity * sizeof(*result->coverage_findings));
        if (!memory) return -1;
        result->coverage_findings = memory;
        result->coverage_finding_capacity = capacity;
    }
    result->coverage_findings[result->coverage_finding_count++] = *finding;
    result->stats.coverage_findings = result->coverage_finding_count;
    result->stats.retained_output_bytes = transition_slice_retained_bytes(result);
    return 0;
}

static void transition_slice_world_point(
    const alea_slice_view_t* view, double u, double v, double point[3]) {
    for (int axis = 0; axis < 3; axis++)
        point[axis] = view->plane.origin[axis] +
            u * view->plane.u_axis[axis] +
            v * view->plane.v_axis[axis];
}

static int transition_slice_probe_coverage(
    alea_system_t* sys, const alea_slice_view_t* view,
    const alea_transition_slice_options_t* options,
    alea_transition_slice_orientation_t orientation,
    size_t ray_index, size_t base_ray_index, uint32_t refinement_depth,
    double transverse_coordinate, const alea_ray_t* ray,
    double ray_t, double bracket_t_enter, double bracket_t_exit,
    transition_slice_coverage_scratch_t* scratch,
    alea_transition_slice_result_t* result) {
    if (options->max_coverage_probes &&
        result->stats.coverage_probes >= options->max_coverage_probes) {
        result->stats.stop_reason =
            ALEA_TRANSITION_SLICE_STOP_MAX_COVERAGE_PROBES;
        return 1;
    }
    double point[3];
    alea_ray_point_at(ray, ray_t, &point[0], &point[1], &point[2]);
    int hit_count = alea_find_all_cells_coverage_chain(
        sys, point[0], point[1], point[2], scratch->hits,
        scratch->occurrence_keys, scratch->parent_occurrence_keys,
        scratch->capacity);
    if (hit_count < 0) return -1;
    result->stats.coverage_probes++;

    alea_transition_slice_coverage_finding_t finding;
    memset(&finding, 0, sizeof(finding));
    finding.orientation = orientation;
    finding.ray_index = ray_index;
    finding.base_ray_index = base_ray_index;
    finding.refinement_depth = refinement_depth;
    finding.transverse_coordinate = transverse_coordinate;
    finding.ray_t = ray_t;
    finding.bracket_t_enter = bracket_t_enter;
    finding.bracket_t_exit = bracket_t_exit;
    memcpy(finding.world_point, point, sizeof(point));
    if (orientation == ALEA_TRANSITION_SLICE_HORIZONTAL) {
        finding.uv[0] = view->u_min + ray_t;
        finding.uv[1] = transverse_coordinate;
    } else {
        finding.uv[0] = transverse_coordinate;
        finding.uv[1] = view->v_min + ray_t;
    }

    if ((size_t)hit_count >= scratch->capacity) {
        finding.kind = ALEA_POINT_COVERAGE_UNRESOLVED;
        finding.truncated = 1;
        finding.target_depth = -1;
        finding.owner_count_lower_bound = scratch->capacity;
        result->stats.truncated_coverage_probes++;
    } else {
        alea_point_coverage_classification_t classification;
        memset(scratch->owner_mask, 0,
               scratch->capacity * sizeof(*scratch->owner_mask));
        if (alea_classify_point_coverage_chain(
                scratch->hits, scratch->occurrence_keys,
                scratch->parent_occurrence_keys, (size_t)hit_count, -1,
                scratch->owner_mask, &classification) != 0)
            return -1;
        finding.kind = classification.kind;
        finding.target_depth = classification.target_depth;
        finding.owner_count_lower_bound = classification.owner_count;
        if (classification.kind == ALEA_POINT_COVERAGE_UNIQUE) {
            result->stats.unique_coverage_probes++;
            return 0;
        }
        if (classification.kind == ALEA_POINT_COVERAGE_GAP &&
            !options->report_unowned_coverage) {
            result->stats.skipped_unowned_coverage_probes++;
            return 0;
        }
        for (int hit = 0; hit < hit_count &&
             finding.owner_count <
                ALEA_TRANSITION_SLICE_COVERAGE_OWNER_CAPACITY; hit++) {
            if (!scratch->owner_mask[hit]) continue;
            size_t owner = finding.owner_count++;
            finding.owner_cell_ids[owner] = scratch->hits[hit].cell_id;
            finding.owner_universe_ids[owner] = scratch->hits[hit].universe_id;
            finding.owner_depths[owner] = scratch->hits[hit].depth;
            finding.owner_occurrence_keys[owner] =
                scratch->occurrence_keys[hit];
            finding.owner_parent_occurrence_keys[owner] =
                scratch->parent_occurrence_keys[hit];
        }
    }
    return transition_slice_coverage_append(result, options, &finding);
}

static size_t transition_slice_event_cap(
    const alea_transition_slice_options_t* options,
    const alea_transition_slice_stats_t* stats,
    alea_transition_slice_stop_reason_t* limiting_reason) {
    uint64_t cap = options->max_events_per_ray;
    *limiting_reason = ALEA_TRANSITION_SLICE_STOP_MAX_EVENTS;
    if (options->max_events) {
        uint64_t remaining = stats->events_checked < options->max_events
            ? options->max_events - stats->events_checked : 0;
        if (!cap || remaining < cap) cap = remaining;
    }
    if (options->max_scratch_bytes) {
        /* Vector growth is geometric. Half the byte-derived count guarantees
         * its capacity remains under the caller's byte ceiling. */
        uint64_t scratch_cap = options->max_scratch_bytes /
            (2u * sizeof(alea_ray_boundary_event_t));
        if (scratch_cap == 0 && options->max_scratch_bytes >=
                sizeof(alea_ray_boundary_event_t)) scratch_cap = 1;
        if (!cap || scratch_cap < cap) {
            cap = scratch_cap;
            *limiting_reason = ALEA_TRANSITION_SLICE_STOP_MAX_SCRATCH_BYTES;
        }
    }
    return cap > SIZE_MAX ? SIZE_MAX : (size_t)cap;
}

static uint64_t transition_slice_hash(uint64_t hash, uint64_t value) {
    value += UINT64_C(0x9e3779b97f4a7c15);
    value = (value ^ (value >> 30)) * UINT64_C(0xbf58476d1ce4e5b9);
    value = (value ^ (value >> 27)) * UINT64_C(0x94d049bb133111eb);
    value ^= value >> 31;
    return hash ^ (value + (hash << 6) + (hash >> 2));
}

static void transition_slice_row_signature(
    const alea_ray_boundary_event_result_t* events,
    transition_slice_row_t* row) {
    uint64_t a = UINT64_C(0xcbf29ce484222325);
    uint64_t b = UINT64_C(0x84222325cbf29ce4);
    for (size_t i = 0; i < events->events.count; i++) {
        const alea_ray_boundary_event_t* event = &events->events.data[i];
#define HASH_EVENT(VALUE) \
        do { \
            uint64_t value_ = (uint64_t)(VALUE); \
            a = transition_slice_hash(a, value_); \
            b = transition_slice_hash(b, value_ ^ (uint64_t)i); \
        } while (0)
        HASH_EVENT(event->kind);
        HASH_EVENT((uint32_t)event->surface_id);
        HASH_EVENT((uint32_t)event->cell_before);
        HASH_EVENT((uint32_t)event->cell_after);
        HASH_EVENT((uint32_t)event->active_cell_id);
        HASH_EVENT((uint32_t)event->active_universe_id);
        HASH_EVENT((uint32_t)event->active_depth);
        HASH_EVENT(event->active_occurrence_key);
        HASH_EVENT(event->active_parent_occurrence_key);
        HASH_EVENT(event->before_occurrence_key);
        HASH_EVENT(event->after_occurrence_key);
        HASH_EVENT(event->local_surface_complete);
        HASH_EVENT(event->local_surface_count);
        for (size_t surface = 0; surface < event->local_surface_count;
             surface++)
            HASH_EVENT((uint32_t)event->local_surface_ids[surface]);
#undef HASH_EVENT
    }
    row->signature_a = transition_slice_hash(a, events->events.count);
    row->signature_b = transition_slice_hash(b, events->events.count);
    row->event_count = events->events.count;
}

static int transition_slice_scan_ray_reuse(
    alea_system_t* sys, const alea_slice_view_t* view,
    const alea_transition_slice_options_t* options,
    alea_transition_slice_orientation_t orientation, size_t ray_index,
    size_t base_ray_index, uint32_t refinement_depth,
    double transverse_coordinate, alea_raycast_result_t* scratch,
    alea_ray_boundary_event_result_t* events,
    transition_slice_coverage_scratch_t* coverage_scratch,
    alea_transition_workspace_t* transition_workspace,
    alea_transition_slice_result_t* result, transition_slice_row_t* row) {
    const double u_span = view->u_max - view->u_min;
    const double v_span = view->v_max - view->v_min;
    double u, v, origin[3], direction[3], length;
    if (orientation == ALEA_TRANSITION_SLICE_HORIZONTAL) {
        u = view->u_min;
        v = transverse_coordinate;
        direction[0] = view->plane.u_axis[0];
        direction[1] = view->plane.u_axis[1];
        direction[2] = view->plane.u_axis[2];
        length = u_span;
    } else {
        u = transverse_coordinate;
        v = view->v_min;
        direction[0] = view->plane.v_axis[0];
        direction[1] = view->plane.v_axis[1];
        direction[2] = view->plane.v_axis[2];
        length = v_span;
    }
    transition_slice_world_point(view, u, v, origin);
    alea_ray_t ray;
    if (alea_ray_init(&ray, origin[0], origin[1], origin[2],
                      direction[0], direction[1], direction[2]) != 0)
        return -1;

    alea_transition_slice_stop_reason_t limiting_reason;
    size_t event_cap = transition_slice_event_cap(
        options, &result->stats, &limiting_reason);
    if (event_cap == 0) {
        result->stats.stop_reason = limiting_reason;
        return 1;
    }
    const alea_ray_boundary_event_options_internal_t event_options = {
        .max_events = event_cap,
        .max_output_bytes = options->max_scratch_bytes,
        .skip_open_side_coverage = true
    };
    if (alea_raycast_selected_boundary_events_with_options_nocache(
            sys, &ray, length, &event_options, scratch, events) != 0) {
        if (alea_error_code() != ALEA_ERR_OVERFLOW) return -1;
        result->stats.stop_reason = limiting_reason;
        return 1;
    }
    result->stats.executed_rays++;
    if (orientation == ALEA_TRANSITION_SLICE_HORIZONTAL)
        result->stats.horizontal_rays_executed++;
    else
        result->stats.vertical_rays_executed++;
    if (events->events.count > result->stats.peak_live_events)
        result->stats.peak_live_events = events->events.count;
    const size_t live_bytes = events->events.capacity * sizeof(*events->events.data);
    if (live_bytes > result->stats.peak_live_event_bytes)
        result->stats.peak_live_event_bytes = live_bytes;
    row->transverse_coordinate = transverse_coordinate;
    row->base_ray_index = base_ray_index;
    row->refinement_depth = refinement_depth;
    row->finding_count = 0;
    transition_slice_row_signature(events, row);
    const size_t findings_before = result->finding_count;
    const size_t coverage_findings_before = result->coverage_finding_count;

    if (options->coverage_uniform_probes_per_ray) {
        const size_t count = options->coverage_uniform_probes_per_ray;
        size_t selected_boundary = 0;
        for (size_t probe = 0; probe < count; probe++) {
            double bin_enter = length * (double)probe / (double)count;
            double bin_exit = length * (double)(probe + 1) / (double)count;
            double ray_t = bin_enter + 0.381966011250105 *
                (bin_exit - bin_enter);
            while (selected_boundary < events->events.count &&
                   events->events.data[selected_boundary].t < ray_t)
                selected_boundary++;
            double t_enter = selected_boundary
                ? events->events.data[selected_boundary - 1].t : 0.0;
            double t_exit = selected_boundary < events->events.count
                ? events->events.data[selected_boundary].t : length;
            int rc = transition_slice_probe_coverage(
                sys, view, options, orientation, ray_index, base_ray_index,
                refinement_depth, transverse_coordinate, &ray,
                ray_t, t_enter, t_exit, coverage_scratch, result);
            if (rc != 0) return rc;
        }
    }
    if (options->coverage_probe_selected_intervals) {
        double t_enter = 0.0;
        for (size_t boundary = 0; boundary <= events->events.count; boundary++) {
            double t_exit = boundary < events->events.count
                ? events->events.data[boundary].t : length;
            if (t_exit - t_enter > 64.0 * RAY_EPSILON) {
                double ray_t = t_enter + 0.381966011250105 *
                    (t_exit - t_enter);
                int rc = transition_slice_probe_coverage(
                    sys, view, options, orientation, ray_index, base_ray_index,
                    refinement_depth, transverse_coordinate, &ray,
                    ray_t, t_enter, t_exit, coverage_scratch, result);
                if (rc != 0) return rc;
            }
            if (t_exit > t_enter) t_enter = t_exit;
        }
    }

    for (size_t event_index = 0; event_index < events->events.count;
         event_index++) {
        const alea_ray_boundary_event_t* event = &events->events.data[event_index];
        result->stats.events_checked++;
        if (event->kind != ALEA_RAY_BOUNDARY_EVENT_PHYSICAL) continue;
        result->stats.physical_events_seen++;
        if (!options->include_void_transitions &&
            (event->cell_before < 0 || event->cell_after < 0)) {
            result->stats.skipped_void_transitions++;
            continue;
        }
        if (options->max_coverage_fallbacks &&
            result->stats.coverage_fallbacks >=
                options->max_coverage_fallbacks) {
            result->stats.stop_reason =
                ALEA_TRANSITION_SLICE_STOP_MAX_COVERAGE_FALLBACKS;
            return 1;
        }

        double previous_t = event_index ? events->events.data[event_index - 1].t : 0.0;
        double next_t = event_index + 1 < events->events.count
            ? events->events.data[event_index + 1].t : length;
        double safe_max = 0.25 * fmin(event->t - previous_t, next_t - event->t);
        alea_transition_result_t transition;
        memset(&transition, 0, sizeof(transition));
        if (!(safe_max > 64.0 * RAY_EPSILON)) {
            transition.kind = ALEA_TRANSITION_AMBIGUOUS_BOUNDARY;
            transition.after_coverage_kind = ALEA_POINT_COVERAGE_UNRESOLVED;
            transition.universe_id = event->active_universe_id;
            transition.current_cell_id = event->active_cell_id;
            transition.primary_surface_id = event->surface_id;
            transition.connecting_surface_id = -1;
            transition.after_cell_id = event->cell_after;
            transition.occurrence_depth = event->active_depth;
            memcpy(transition.crossing_point, event->local_point,
                   sizeof(transition.crossing_point));
            memcpy(transition.direction, event->local_direction,
                   sizeof(transition.direction));
        } else {
            alea_transition_options_t transition_options;
            alea_transition_options_init(&transition_options);
            if (options->probe_distance > 0.0)
                transition_options.probe_distance = options->probe_distance;
            if (transition_options.probe_distance > safe_max)
                transition_options.probe_distance = safe_max;
            transition_options.max_probe_distance = safe_max;
            if (options->max_probe_distance > 0.0 &&
                transition_options.max_probe_distance >
                    options->max_probe_distance)
                transition_options.max_probe_distance =
                    options->max_probe_distance;
            transition_options.max_coverage_hits = options->max_coverage_hits;
            if (options->max_coverage_fallbacks)
                transition_options.max_coverage_fallbacks =
                    options->max_coverage_fallbacks -
                    result->stats.coverage_fallbacks;
            if (alea_check_selected_boundary_event_transition_reuse_nocache(
                    sys, event, &transition_options, &transition,
                    transition_workspace) != 0)
                return -1;
        }
        result->stats.coverage_fallbacks += transition.coverage_fallbacks;
        if (transition.kind == ALEA_TRANSITION_VALID) {
            result->stats.valid_transitions++;
            continue;
        }
        alea_transition_slice_finding_t finding;
        memset(&finding, 0, sizeof(finding));
        finding.transition = transition;
        finding.orientation = orientation;
        finding.ray_index = ray_index;
        finding.event_index = event_index;
        finding.base_ray_index = base_ray_index;
        finding.refinement_depth = refinement_depth;
        finding.transverse_coordinate = transverse_coordinate;
        finding.ray_t = event->t;
        if (orientation == ALEA_TRANSITION_SLICE_HORIZONTAL) {
            finding.uv[0] = view->u_min + event->t;
            finding.uv[1] = v;
        } else {
            finding.uv[0] = u;
            finding.uv[1] = view->v_min + event->t;
        }
        alea_ray_point_at(&ray, event->t,
                          &finding.world_point[0], &finding.world_point[1],
                          &finding.world_point[2]);
        int append_rc = transition_slice_append(result, options, &finding);
        if (append_rc != 0) return append_rc;
        if (transition.kind == ALEA_TRANSITION_TRUNCATED &&
            options->max_coverage_fallbacks &&
            result->stats.coverage_fallbacks >=
                options->max_coverage_fallbacks) {
            result->stats.stop_reason =
                ALEA_TRANSITION_SLICE_STOP_MAX_COVERAGE_FALLBACKS;
            return 1;
        }
    }
    row->finding_count = result->finding_count - findings_before +
        result->coverage_finding_count - coverage_findings_before;
    return 0;
}

static int transition_slice_rows_differ(
    const transition_slice_row_t* first,
    const transition_slice_row_t* second, uint32_t signals) {
    if ((signals & ALEA_TRANSITION_SLICE_REFINE_SIGNATURE) &&
        (first->event_count != second->event_count ||
         first->signature_a != second->signature_a ||
         first->signature_b != second->signature_b))
        return 1;
    return (signals & ALEA_TRANSITION_SLICE_REFINE_FINDING) &&
        (first->finding_count != 0 || second->finding_count != 0);
}

static int transition_slice_append_frontier(
    const alea_slice_view_t* view,
    const alea_transition_slice_options_t* options,
    alea_transition_slice_orientation_t orientation,
    const transition_slice_row_t* first,
    const transition_slice_row_t* second,
    alea_transition_slice_result_t* result) {
    if (!options->enable_critical_refinement) return 0;
    if (options->max_refinement_frontiers &&
        result->refinement_frontier_count >=
            options->max_refinement_frontiers) {
        result->stats.omitted_refinement_frontiers++;
        if (result->stats.critical_stop_reason ==
            ALEA_TRANSITION_SLICE_CRITICAL_NONE)
            result->stats.critical_stop_reason =
                ALEA_TRANSITION_SLICE_CRITICAL_MAX_FRONTIERS;
        return 0;
    }
    size_t retained = transition_slice_retained_bytes(result);
    if (options->max_output_bytes &&
        (retained > options->max_output_bytes ||
         sizeof(*result->refinement_frontiers) >
             options->max_output_bytes - retained)) {
        result->stats.omitted_refinement_frontiers++;
        if (result->stats.critical_stop_reason ==
            ALEA_TRANSITION_SLICE_CRITICAL_NONE)
            result->stats.critical_stop_reason =
                ALEA_TRANSITION_SLICE_CRITICAL_MAX_OUTPUT_BYTES;
        return 0;
    }
    if (result->refinement_frontier_count ==
        result->refinement_frontier_capacity) {
        size_t capacity = result->refinement_frontier_capacity
            ? result->refinement_frontier_capacity * 2 : 16;
        if (options->max_refinement_frontiers &&
            capacity > options->max_refinement_frontiers)
            capacity = options->max_refinement_frontiers;
        void* memory = realloc(
            result->refinement_frontiers,
            capacity * sizeof(*result->refinement_frontiers));
        if (!memory) return -1;
        result->refinement_frontiers = memory;
        result->refinement_frontier_capacity = capacity;
    }
    alea_transition_slice_refinement_frontier_t* frontier =
        &result->refinement_frontiers[result->refinement_frontier_count++];
    memset(frontier, 0, sizeof(*frontier));
    frontier->orientation = orientation;
    frontier->refinement_depth = first->refinement_depth >
            second->refinement_depth
        ? first->refinement_depth : second->refinement_depth;
    frontier->transverse_min = first->transverse_coordinate;
    frontier->transverse_max = second->transverse_coordinate;
    frontier->signature_a[0] = first->signature_a;
    frontier->signature_a[1] = second->signature_a;
    frontier->signature_b[0] = first->signature_b;
    frontier->signature_b[1] = second->signature_b;
    frontier->max_event_count = first->event_count > second->event_count
        ? first->event_count : second->event_count;
    if (orientation == ALEA_TRANSITION_SLICE_HORIZONTAL) {
        frontier->uv_min[0] = view->u_min;
        frontier->uv_max[0] = view->u_max;
        frontier->uv_min[1] = frontier->transverse_min;
        frontier->uv_max[1] = frontier->transverse_max;
    } else {
        frontier->uv_min[0] = frontier->transverse_min;
        frontier->uv_max[0] = frontier->transverse_max;
        frontier->uv_min[1] = view->v_min;
        frontier->uv_max[1] = view->v_max;
    }
    result->stats.refinement_frontiers = result->refinement_frontier_count;
    result->stats.retained_output_bytes = transition_slice_retained_bytes(result);
    return 0;
}

static int transition_slice_retain_differing_frontiers(
    const alea_slice_view_t* view,
    const alea_transition_slice_options_t* options,
    alea_transition_slice_orientation_t orientation,
    const transition_slice_row_t* rows, size_t row_count,
    alea_transition_slice_result_t* result) {
    for (size_t row = 0; row + 1 < row_count; row++) {
        if (!transition_slice_rows_differ(
                &rows[row], &rows[row + 1], options->refine_signals))
            continue;
        if (transition_slice_append_frontier(
                view, options, orientation, &rows[row], &rows[row + 1],
                result) != 0)
            return -1;
    }
    return 0;
}

static void transition_slice_record_row_scratch(
    alea_transition_slice_result_t* result, size_t row_bytes) {
    if (row_bytes > result->stats.peak_row_scratch_bytes)
        result->stats.peak_row_scratch_bytes = row_bytes;
    size_t total = row_bytes;
    if (SIZE_MAX - total < result->stats.peak_live_event_bytes)
        total = SIZE_MAX;
    else
        total += result->stats.peak_live_event_bytes;
    if (total > result->stats.peak_scratch_bytes)
        result->stats.peak_scratch_bytes = total;
}

static int transition_slice_status_rank(
    alea_transition_slice_refinement_status_t status) {
    switch (status) {
    case ALEA_TRANSITION_SLICE_REFINEMENT_STOPPED: return 4;
    case ALEA_TRANSITION_SLICE_REFINEMENT_MAX_DEPTH: return 3;
    case ALEA_TRANSITION_SLICE_REFINEMENT_MIN_SPACING: return 2;
    case ALEA_TRANSITION_SLICE_REFINEMENT_CONVERGED: return 1;
    case ALEA_TRANSITION_SLICE_REFINEMENT_NOT_REQUESTED: return 0;
    }
    return 4;
}

static void transition_slice_merge_status(
    alea_transition_slice_result_t* result,
    alea_transition_slice_refinement_status_t status) {
    if (transition_slice_status_rank(status) >
        transition_slice_status_rank(result->stats.refinement_status))
        result->stats.refinement_status = status;
}

static int transition_slice_compare_int(int first, int second) {
    return first < second ? -1 : first > second ? 1 : 0;
}

static int transition_slice_compare_size(size_t first, size_t second) {
    return first < second ? -1 : first > second ? 1 : 0;
}

static int transition_slice_finding_compare(const void* first_pointer,
                                            const void* second_pointer) {
    const alea_transition_slice_finding_t* first = first_pointer;
    const alea_transition_slice_finding_t* second = second_pointer;
#define COMPARE_INT(FIELD) \
    do { \
        int comparison = transition_slice_compare_int( \
            (int)first->FIELD, (int)second->FIELD); \
        if (comparison) return comparison; \
    } while (0)
    COMPARE_INT(transition.kind);
    COMPARE_INT(orientation);
    COMPARE_INT(transition.universe_id);
    COMPARE_INT(transition.current_cell_id);
    COMPARE_INT(transition.after_cell_id);
    COMPARE_INT(transition.primary_surface_id);
    COMPARE_INT(transition.connecting_surface_id);
#undef COMPARE_INT
    if (first->transition.current_occurrence_key !=
        second->transition.current_occurrence_key)
        return first->transition.current_occurrence_key <
                second->transition.current_occurrence_key ? -1 : 1;
    for (int axis = 0; axis < 2; axis++) {
        if (first->uv[axis] != second->uv[axis])
            return first->uv[axis] < second->uv[axis] ? -1 : 1;
    }
    int comparison = transition_slice_compare_size(
        first->ray_index, second->ray_index);
    if (comparison) return comparison;
    return transition_slice_compare_size(
        first->event_index, second->event_index);
}

static int transition_slice_findings_same_component(
    const alea_transition_slice_finding_t* first,
    const alea_transition_slice_finding_t* second) {
    return first->transition.kind == second->transition.kind &&
        first->orientation == second->orientation &&
        first->transition.universe_id == second->transition.universe_id &&
        first->transition.current_cell_id == second->transition.current_cell_id &&
        first->transition.after_cell_id == second->transition.after_cell_id &&
        first->transition.primary_surface_id ==
            second->transition.primary_surface_id &&
        first->transition.connecting_surface_id ==
            second->transition.connecting_surface_id &&
        first->transition.current_occurrence_key ==
            second->transition.current_occurrence_key;
}

static int transition_slice_build_components(
    const alea_transition_slice_options_t* options,
    alea_transition_slice_result_t* result) {
    if (result->finding_count == 0) return 0;
    qsort(result->findings, result->finding_count,
          sizeof(*result->findings), transition_slice_finding_compare);
    size_t required = 1;
    for (size_t i = 1; i < result->finding_count; i++)
        if (!transition_slice_findings_same_component(
                &result->findings[i - 1], &result->findings[i]))
            required++;
    size_t capacity = required;
    if (options->max_components && capacity > options->max_components)
        capacity = (size_t)options->max_components;
    if (options->max_output_bytes) {
        size_t retained = transition_slice_retained_bytes(result);
        size_t available = retained < options->max_output_bytes
            ? (size_t)(options->max_output_bytes - retained) : 0;
        size_t byte_capacity = available / sizeof(*result->components);
        if (capacity > byte_capacity) capacity = byte_capacity;
    }
    if (capacity > 0) {
        result->components = calloc(capacity, sizeof(*result->components));
        if (!result->components) return -1;
        result->component_capacity = capacity;
    }
    size_t finding = 0;
    while (finding < result->finding_count &&
           result->component_count < capacity) {
        size_t end = finding + 1;
        while (end < result->finding_count &&
               transition_slice_findings_same_component(
                   &result->findings[finding], &result->findings[end]))
            end++;
        const alea_transition_slice_finding_t* representative =
            &result->findings[finding];
        alea_transition_slice_component_t* component =
            &result->components[result->component_count++];
        component->kind = representative->transition.kind;
        component->orientation = representative->orientation;
        component->universe_id = representative->transition.universe_id;
        component->current_cell_id = representative->transition.current_cell_id;
        component->after_cell_id = representative->transition.after_cell_id;
        component->primary_surface_id =
            representative->transition.primary_surface_id;
        component->connecting_surface_id =
            representative->transition.connecting_surface_id;
        component->current_occurrence_key =
            representative->transition.current_occurrence_key;
        component->first_finding_index = finding;
        component->finding_count = end - finding;
        for (int axis = 0; axis < 2; axis++)
            component->uv_min[axis] = component->uv_max[axis] =
                representative->uv[axis];
        for (int axis = 0; axis < 3; axis++)
            component->world_min[axis] = component->world_max[axis] =
                representative->world_point[axis];
        for (size_t i = finding; i < end; i++) {
            const alea_transition_slice_finding_t* item = &result->findings[i];
            if (component->max_refinement_depth < item->refinement_depth)
                component->max_refinement_depth = item->refinement_depth;
            for (int axis = 0; axis < 2; axis++) {
                if (component->uv_min[axis] > item->uv[axis])
                    component->uv_min[axis] = item->uv[axis];
                if (component->uv_max[axis] < item->uv[axis])
                    component->uv_max[axis] = item->uv[axis];
            }
            for (int axis = 0; axis < 3; axis++) {
                if (component->world_min[axis] > item->world_point[axis])
                    component->world_min[axis] = item->world_point[axis];
                if (component->world_max[axis] < item->world_point[axis])
                    component->world_max[axis] = item->world_point[axis];
            }
        }
        finding = end;
    }
    result->stats.components = result->component_count;
    result->stats.retained_output_bytes = transition_slice_retained_bytes(result);
    if (result->component_count < required &&
        result->stats.stop_reason == ALEA_TRANSITION_SLICE_STOP_NONE) {
        result->stats.stop_reason = options->max_components &&
                result->component_count >= options->max_components
            ? ALEA_TRANSITION_SLICE_STOP_MAX_COMPONENTS
            : ALEA_TRANSITION_SLICE_STOP_MAX_OUTPUT_BYTES;
        return 1;
    }
    return 0;
}

static int transition_slice_coverage_key_compare(
    const alea_transition_slice_coverage_finding_t* first,
    const alea_transition_slice_coverage_finding_t* second) {
    int comparison = transition_slice_compare_int(
        (int)first->kind, (int)second->kind);
    if (comparison) return comparison;
    comparison = transition_slice_compare_int(first->truncated, second->truncated);
    if (comparison) return comparison;
    comparison = transition_slice_compare_int(
        (int)first->orientation, (int)second->orientation);
    if (comparison) return comparison;
    comparison = transition_slice_compare_size(
        first->owner_count_lower_bound, second->owner_count_lower_bound);
    if (comparison) return comparison;
    comparison = transition_slice_compare_size(
        first->owner_count, second->owner_count);
    if (comparison) return comparison;
    for (size_t owner = 0; owner < first->owner_count; owner++) {
        comparison = transition_slice_compare_int(
            first->owner_cell_ids[owner], second->owner_cell_ids[owner]);
        if (comparison) return comparison;
        comparison = transition_slice_compare_int(
            first->owner_universe_ids[owner],
            second->owner_universe_ids[owner]);
        if (comparison) return comparison;
        if (first->owner_occurrence_keys[owner] !=
            second->owner_occurrence_keys[owner])
            return first->owner_occurrence_keys[owner] <
                second->owner_occurrence_keys[owner] ? -1 : 1;
    }
    return 0;
}

static int transition_slice_coverage_finding_compare(
    const void* first_pointer, const void* second_pointer) {
    const alea_transition_slice_coverage_finding_t* first = first_pointer;
    const alea_transition_slice_coverage_finding_t* second = second_pointer;
    int comparison = transition_slice_coverage_key_compare(first, second);
    if (comparison) return comparison;
    for (int axis = 0; axis < 2; axis++) {
        if (first->uv[axis] != second->uv[axis])
            return first->uv[axis] < second->uv[axis] ? -1 : 1;
    }
    comparison = transition_slice_compare_size(
        first->ray_index, second->ray_index);
    if (comparison) return comparison;
    return first->ray_t < second->ray_t ? -1 : first->ray_t > second->ray_t;
}

static int transition_slice_build_coverage_components(
    const alea_transition_slice_options_t* options,
    alea_transition_slice_result_t* result) {
    if (result->coverage_finding_count == 0) return 0;
    qsort(result->coverage_findings, result->coverage_finding_count,
          sizeof(*result->coverage_findings),
          transition_slice_coverage_finding_compare);
    size_t required = 1;
    for (size_t i = 1; i < result->coverage_finding_count; i++)
        if (transition_slice_coverage_key_compare(
                &result->coverage_findings[i - 1],
                &result->coverage_findings[i]) != 0)
            required++;
    size_t capacity = required;
    if (options->max_coverage_components &&
        capacity > options->max_coverage_components)
        capacity = (size_t)options->max_coverage_components;
    if (options->max_output_bytes) {
        size_t retained = transition_slice_retained_bytes(result);
        size_t available = retained < options->max_output_bytes
            ? (size_t)(options->max_output_bytes - retained) : 0;
        size_t byte_capacity = available / sizeof(*result->coverage_components);
        if (capacity > byte_capacity) capacity = byte_capacity;
    }
    if (capacity) {
        result->coverage_components =
            calloc(capacity, sizeof(*result->coverage_components));
        if (!result->coverage_components) return -1;
        result->coverage_component_capacity = capacity;
    }
    size_t finding = 0;
    while (finding < result->coverage_finding_count &&
           result->coverage_component_count < capacity) {
        size_t end = finding + 1;
        while (end < result->coverage_finding_count &&
               transition_slice_coverage_key_compare(
                   &result->coverage_findings[finding],
                   &result->coverage_findings[end]) == 0)
            end++;
        const alea_transition_slice_coverage_finding_t* representative =
            &result->coverage_findings[finding];
        alea_transition_slice_coverage_component_t* component =
            &result->coverage_components[result->coverage_component_count++];
        component->kind = representative->kind;
        component->truncated = representative->truncated;
        component->orientation = representative->orientation;
        component->first_finding_index = finding;
        component->finding_count = end - finding;
        component->owner_count_lower_bound =
            representative->owner_count_lower_bound;
        for (int axis = 0; axis < 2; axis++)
            component->uv_min[axis] = component->uv_max[axis] =
                representative->uv[axis];
        for (int axis = 0; axis < 3; axis++)
            component->world_min[axis] = component->world_max[axis] =
                representative->world_point[axis];
        for (size_t i = finding; i < end; i++) {
            const alea_transition_slice_coverage_finding_t* item =
                &result->coverage_findings[i];
            if (component->max_refinement_depth < item->refinement_depth)
                component->max_refinement_depth = item->refinement_depth;
            for (int axis = 0; axis < 2; axis++) {
                if (component->uv_min[axis] > item->uv[axis])
                    component->uv_min[axis] = item->uv[axis];
                if (component->uv_max[axis] < item->uv[axis])
                    component->uv_max[axis] = item->uv[axis];
            }
            for (int axis = 0; axis < 3; axis++) {
                if (component->world_min[axis] > item->world_point[axis])
                    component->world_min[axis] = item->world_point[axis];
                if (component->world_max[axis] < item->world_point[axis])
                    component->world_max[axis] = item->world_point[axis];
            }
        }
        finding = end;
    }
    result->stats.coverage_components = result->coverage_component_count;
    result->stats.retained_output_bytes = transition_slice_retained_bytes(result);
    if (result->coverage_component_count < required &&
        result->stats.stop_reason == ALEA_TRANSITION_SLICE_STOP_NONE) {
        result->stats.stop_reason = options->max_coverage_components &&
                result->coverage_component_count >=
                    options->max_coverage_components
            ? ALEA_TRANSITION_SLICE_STOP_MAX_COMPONENTS
            : ALEA_TRANSITION_SLICE_STOP_MAX_OUTPUT_BYTES;
        return 1;
    }
    return 0;
}

typedef struct {
    alea_transition_slice_orientation_t orientation;
    size_t ray_index;
    double ray_t;
    size_t component_index;
} transition_slice_component_boundary_ref_t;

typedef struct {
    alea_transition_slice_orientation_t orientation;
    size_t ray_index;
    double ray_t;
    size_t component_index;
    uint32_t side;
} transition_slice_coverage_endpoint_ref_t;

static int transition_slice_boundary_ref_compare(const void* first_pointer,
                                                 const void* second_pointer) {
    const transition_slice_component_boundary_ref_t* first = first_pointer;
    const transition_slice_component_boundary_ref_t* second = second_pointer;
    int comparison = transition_slice_compare_int(
        (int)first->orientation, (int)second->orientation);
    if (comparison) return comparison;
    comparison = transition_slice_compare_size(first->ray_index,
                                                second->ray_index);
    if (comparison) return comparison;
    if (first->ray_t != second->ray_t)
        return first->ray_t < second->ray_t ? -1 : 1;
    return transition_slice_compare_size(first->component_index,
                                         second->component_index);
}

static int transition_slice_endpoint_ref_compare(const void* first_pointer,
                                                 const void* second_pointer) {
    const transition_slice_coverage_endpoint_ref_t* first = first_pointer;
    const transition_slice_coverage_endpoint_ref_t* second = second_pointer;
    int comparison = transition_slice_compare_int(
        (int)first->orientation, (int)second->orientation);
    if (comparison) return comparison;
    comparison = transition_slice_compare_size(first->ray_index,
                                                second->ray_index);
    if (comparison) return comparison;
    if (first->ray_t != second->ray_t)
        return first->ray_t < second->ray_t ? -1 : 1;
    comparison = transition_slice_compare_size(first->component_index,
                                                second->component_index);
    if (comparison) return comparison;
    return transition_slice_compare_int((int)first->side, (int)second->side);
}

static int transition_slice_link_compare(const void* first_pointer,
                                         const void* second_pointer) {
    const alea_transition_slice_component_link_t* first = first_pointer;
    const alea_transition_slice_component_link_t* second = second_pointer;
    int comparison = transition_slice_compare_size(
        first->transition_component_index,
        second->transition_component_index);
    if (comparison) return comparison;
    comparison = transition_slice_compare_size(
        first->coverage_component_index, second->coverage_component_index);
    if (comparison) return comparison;
    return transition_slice_compare_int((int)first->boundary_sides,
                                        (int)second->boundary_sides);
}

static int transition_slice_same_ray(
    alea_transition_slice_orientation_t first_orientation,
    size_t first_ray, alea_transition_slice_orientation_t second_orientation,
    size_t second_ray) {
    return first_orientation == second_orientation && first_ray == second_ray;
}

static int transition_slice_ray_key_compare(
    alea_transition_slice_orientation_t first_orientation,
    size_t first_ray, alea_transition_slice_orientation_t second_orientation,
    size_t second_ray) {
    int comparison = transition_slice_compare_int(
        (int)first_orientation, (int)second_orientation);
    return comparison ? comparison :
        transition_slice_compare_size(first_ray, second_ray);
}

static int transition_slice_build_component_links(
    const alea_transition_slice_options_t* options,
    alea_transition_slice_result_t* result) {
    if (result->component_count == 0 ||
        result->coverage_component_count == 0)
        return 0;
    if (result->coverage_finding_count > SIZE_MAX / 2) return -1;
    const size_t boundary_count = result->finding_count;
    const size_t endpoint_count = result->coverage_finding_count * 2;
    size_t scratch_bytes = 0;
    if (transition_slice_add_bytes(
            &scratch_bytes, boundary_count,
            sizeof(transition_slice_component_boundary_ref_t)) ||
        transition_slice_add_bytes(
            &scratch_bytes, endpoint_count,
            sizeof(transition_slice_coverage_endpoint_ref_t)) ||
        transition_slice_add_bytes(
            &scratch_bytes, endpoint_count,
            sizeof(alea_transition_slice_component_link_t)))
        return -1;
    if (options->max_scratch_bytes &&
        scratch_bytes > options->max_scratch_bytes) {
        if (result->stats.stop_reason == ALEA_TRANSITION_SLICE_STOP_NONE)
            result->stats.stop_reason =
                ALEA_TRANSITION_SLICE_STOP_MAX_SCRATCH_BYTES;
        return 1;
    }
    if (scratch_bytes > result->stats.peak_scratch_bytes)
        result->stats.peak_scratch_bytes = scratch_bytes;

    transition_slice_component_boundary_ref_t* boundaries =
        malloc(boundary_count * sizeof(*boundaries));
    transition_slice_coverage_endpoint_ref_t* endpoints =
        malloc(endpoint_count * sizeof(*endpoints));
    alea_transition_slice_component_link_t* raw_links =
        malloc(endpoint_count * sizeof(*raw_links));
    if (!boundaries || !endpoints || !raw_links) {
        free(boundaries); free(endpoints); free(raw_links);
        return -1;
    }
    size_t boundary = 0;
    for (size_t component_index = 0;
         component_index < result->component_count; component_index++) {
        const alea_transition_slice_component_t* component =
            &result->components[component_index];
        const size_t end = component->first_finding_index +
            component->finding_count;
        for (size_t finding_index = component->first_finding_index;
             finding_index < end; finding_index++) {
            const alea_transition_slice_finding_t* finding =
                &result->findings[finding_index];
            boundaries[boundary++] =
                (transition_slice_component_boundary_ref_t){
                    .orientation = finding->orientation,
                    .ray_index = finding->ray_index,
                    .ray_t = finding->ray_t,
                    .component_index = component_index
                };
        }
    }
    size_t endpoint = 0;
    for (size_t component_index = 0;
         component_index < result->coverage_component_count;
         component_index++) {
        const alea_transition_slice_coverage_component_t* component =
            &result->coverage_components[component_index];
        const size_t end = component->first_finding_index +
            component->finding_count;
        for (size_t finding_index = component->first_finding_index;
             finding_index < end; finding_index++) {
            const alea_transition_slice_coverage_finding_t* finding =
                &result->coverage_findings[finding_index];
            endpoints[endpoint++] =
                (transition_slice_coverage_endpoint_ref_t){
                    .orientation = finding->orientation,
                    .ray_index = finding->ray_index,
                    .ray_t = finding->bracket_t_enter,
                    .component_index = component_index,
                    .side = ALEA_TRANSITION_SLICE_LINK_ENTER
                };
            endpoints[endpoint++] =
                (transition_slice_coverage_endpoint_ref_t){
                    .orientation = finding->orientation,
                    .ray_index = finding->ray_index,
                    .ray_t = finding->bracket_t_exit,
                    .component_index = component_index,
                    .side = ALEA_TRANSITION_SLICE_LINK_EXIT
                };
        }
    }
    qsort(boundaries, boundary, sizeof(*boundaries),
          transition_slice_boundary_ref_compare);
    qsort(endpoints, endpoint, sizeof(*endpoints),
          transition_slice_endpoint_ref_compare);

    size_t boundary_cursor = 0, raw_count = 0;
    for (size_t endpoint_index = 0; endpoint_index < endpoint;
         endpoint_index++) {
        const transition_slice_coverage_endpoint_ref_t* current =
            &endpoints[endpoint_index];
        while (boundary_cursor < boundary &&
               transition_slice_ray_key_compare(
                   boundaries[boundary_cursor].orientation,
                   boundaries[boundary_cursor].ray_index,
                   current->orientation, current->ray_index) < 0)
            boundary_cursor++;
        size_t candidate = boundary_cursor;
        while (candidate < boundary &&
               transition_slice_same_ray(
                   boundaries[candidate].orientation,
                   boundaries[candidate].ray_index,
                   current->orientation, current->ray_index)) {
            const double scale = 1.0 + fmax(fabs(current->ray_t),
                                            fabs(boundaries[candidate].ray_t));
            const double tolerance = 64.0 * RAY_EPSILON * scale;
            if (boundaries[candidate].ray_t < current->ray_t - tolerance) {
                candidate++;
                continue;
            }
            if (boundaries[candidate].ray_t > current->ray_t + tolerance)
                break;
            raw_links[raw_count++] =
                (alea_transition_slice_component_link_t){
                    .transition_component_index =
                        boundaries[candidate].component_index,
                    .coverage_component_index = current->component_index,
                    .boundary_sides = current->side,
                    .witness_pair_count = 1
                };
            break;
        }
        boundary_cursor = candidate;
    }
    free(boundaries);
    free(endpoints);
    if (raw_count == 0) {
        free(raw_links);
        return 0;
    }
    qsort(raw_links, raw_count, sizeof(*raw_links),
          transition_slice_link_compare);
    size_t required = 1;
    for (size_t i = 1; i < raw_count; i++)
        if (raw_links[i - 1].transition_component_index !=
                raw_links[i].transition_component_index ||
            raw_links[i - 1].coverage_component_index !=
                raw_links[i].coverage_component_index)
            required++;
    size_t capacity = required;
    if (options->max_component_links &&
        capacity > options->max_component_links)
        capacity = (size_t)options->max_component_links;
    if (options->max_output_bytes) {
        const size_t retained = transition_slice_retained_bytes(result);
        const size_t available = retained < options->max_output_bytes
            ? (size_t)(options->max_output_bytes - retained) : 0;
        const size_t byte_capacity =
            available / sizeof(*result->component_links);
        if (capacity > byte_capacity) capacity = byte_capacity;
    }
    if (capacity) {
        result->component_links =
            calloc(capacity, sizeof(*result->component_links));
        if (!result->component_links) {
            free(raw_links);
            return -1;
        }
        result->component_link_capacity = capacity;
    }
    size_t raw = 0;
    while (raw < raw_count && result->component_link_count < capacity) {
        size_t end = raw + 1;
        alea_transition_slice_component_link_t link = raw_links[raw];
        while (end < raw_count &&
               raw_links[end].transition_component_index ==
                   link.transition_component_index &&
               raw_links[end].coverage_component_index ==
                   link.coverage_component_index) {
            link.boundary_sides |= raw_links[end].boundary_sides;
            link.witness_pair_count += raw_links[end].witness_pair_count;
            end++;
        }
        result->component_links[result->component_link_count++] = link;
        raw = end;
    }
    free(raw_links);
    result->stats.component_links = result->component_link_count;
    result->stats.retained_output_bytes = transition_slice_retained_bytes(result);
    if (result->component_link_count < required &&
        result->stats.stop_reason == ALEA_TRANSITION_SLICE_STOP_NONE) {
        result->stats.stop_reason = options->max_component_links &&
                result->component_link_count >= options->max_component_links
            ? ALEA_TRANSITION_SLICE_STOP_MAX_COMPONENT_LINKS
            : ALEA_TRANSITION_SLICE_STOP_MAX_OUTPUT_BYTES;
        return 1;
    }
    return 0;
}

typedef struct {
    double uv_min[2];
    double uv_max[2];
    alea_transition_slice_tile_source_kind_t kind;
    size_t source_index;
    int priority;
} transition_slice_tile_seed_t;

static int transition_slice_tile_seed_compare(const void* first_pointer,
                                              const void* second_pointer) {
    const transition_slice_tile_seed_t* first = first_pointer;
    const transition_slice_tile_seed_t* second = second_pointer;
    int comparison = transition_slice_compare_int(first->priority,
                                                  second->priority);
    if (comparison) return comparison;
    for (int axis = 0; axis < 2; axis++) {
        if (first->uv_min[axis] != second->uv_min[axis])
            return first->uv_min[axis] < second->uv_min[axis] ? -1 : 1;
    }
    for (int axis = 0; axis < 2; axis++) {
        if (first->uv_max[axis] != second->uv_max[axis])
            return first->uv_max[axis] < second->uv_max[axis] ? -1 : 1;
    }
    comparison = transition_slice_compare_int((int)first->kind,
                                               (int)second->kind);
    return comparison ? comparison :
        transition_slice_compare_size(first->source_index,
                                      second->source_index);
}

static int transition_slice_rectangles_touch(
    const double first_min[2], const double first_max[2],
    const double second_min[2], const double second_max[2]) {
    return first_min[0] <= second_max[0] && second_min[0] <= first_max[0] &&
        first_min[1] <= second_max[1] && second_min[1] <= first_max[1];
}

static void transition_slice_critical_stop(
    alea_transition_slice_result_t* result,
    alea_transition_slice_critical_stop_reason_t reason) {
    if (result->stats.critical_stop_reason ==
        ALEA_TRANSITION_SLICE_CRITICAL_NONE)
        result->stats.critical_stop_reason = reason;
}

static int transition_slice_build_critical_tiles(
    const alea_slice_view_t* view,
    const alea_transition_slice_options_t* options,
    alea_transition_slice_result_t* result) {
    if (!options->enable_critical_refinement) return 0;
    size_t seed_count = result->component_count;
    if (seed_count > SIZE_MAX - result->coverage_component_count ||
        seed_count + result->coverage_component_count >
            SIZE_MAX - result->refinement_frontier_count)
        return -1;
    seed_count += result->coverage_component_count +
        result->refinement_frontier_count;
    if (options->critical_full_view) {
        if (seed_count == SIZE_MAX) return -1;
        seed_count++;
    }
    result->stats.critical_tile_seeds = seed_count;
    if (seed_count == 0) return 0;
    size_t scratch_bytes = 0;
    if (transition_slice_add_bytes(&scratch_bytes, seed_count,
                                   sizeof(transition_slice_tile_seed_t)) ||
        transition_slice_add_bytes(&scratch_bytes, seed_count,
                                   sizeof(alea_transition_slice_critical_tile_t)))
        return -1;
    if (options->max_critical_scratch_bytes &&
        scratch_bytes > options->max_critical_scratch_bytes) {
        transition_slice_critical_stop(
            result, ALEA_TRANSITION_SLICE_CRITICAL_MAX_SCRATCH_BYTES);
        result->stats.omitted_critical_tile_sources = seed_count;
        return 0;
    }
    result->stats.peak_critical_scratch_bytes = scratch_bytes;
    if (result->stats.peak_scratch_bytes < scratch_bytes)
        result->stats.peak_scratch_bytes = scratch_bytes;
    transition_slice_tile_seed_t* seeds =
        calloc(seed_count, sizeof(*seeds));
    alea_transition_slice_critical_tile_t* merged =
        calloc(seed_count, sizeof(*merged));
    if (!seeds || !merged) { free(seeds); free(merged); return -1; }
    size_t seed = 0;
    if (options->critical_full_view) {
        transition_slice_tile_seed_t* item = &seeds[seed++];
        item->uv_min[0] = view->u_min;
        item->uv_min[1] = view->v_min;
        item->uv_max[0] = view->u_max;
        item->uv_max[1] = view->v_max;
        item->kind = ALEA_TRANSITION_SLICE_TILE_SOURCE_FULL_VIEW;
        item->priority = -1;
    }
    for (size_t i = 0; i < result->component_count; i++) {
        const alea_transition_slice_component_t* source =
            &result->components[i];
        transition_slice_tile_seed_t* item = &seeds[seed++];
        memcpy(item->uv_min, source->uv_min, sizeof(item->uv_min));
        memcpy(item->uv_max, source->uv_max, sizeof(item->uv_max));
        item->kind = ALEA_TRANSITION_SLICE_TILE_SOURCE_TRANSITION_COMPONENT;
        item->source_index = i;
        item->priority = source->kind == ALEA_TRANSITION_TRUNCATED ||
                source->kind == ALEA_TRANSITION_UNRESOLVED ||
                source->kind == ALEA_TRANSITION_AMBIGUOUS_BOUNDARY
            ? 0 : 1;
    }
    for (size_t i = 0; i < result->coverage_component_count; i++) {
        const alea_transition_slice_coverage_component_t* source =
            &result->coverage_components[i];
        transition_slice_tile_seed_t* item = &seeds[seed++];
        memcpy(item->uv_min, source->uv_min, sizeof(item->uv_min));
        memcpy(item->uv_max, source->uv_max, sizeof(item->uv_max));
        const size_t end = source->first_finding_index + source->finding_count;
        for (size_t j = source->first_finding_index; j < end; j++) {
            const alea_transition_slice_coverage_finding_t* finding =
                &result->coverage_findings[j];
            const int axis = finding->orientation ==
                ALEA_TRANSITION_SLICE_HORIZONTAL ? 0 : 1;
            const double offset = axis == 0 ? view->u_min : view->v_min;
            const double low = offset + finding->bracket_t_enter;
            const double high = offset + finding->bracket_t_exit;
            if (item->uv_min[axis] > low) item->uv_min[axis] = low;
            if (item->uv_max[axis] < high) item->uv_max[axis] = high;
        }
        item->kind = ALEA_TRANSITION_SLICE_TILE_SOURCE_COVERAGE_COMPONENT;
        item->source_index = i;
        item->priority = source->truncated ||
            source->kind == ALEA_POINT_COVERAGE_UNRESOLVED ? 0 : 2;
    }
    for (size_t i = 0; i < result->refinement_frontier_count; i++) {
        const alea_transition_slice_refinement_frontier_t* source =
            &result->refinement_frontiers[i];
        transition_slice_tile_seed_t* item = &seeds[seed++];
        memcpy(item->uv_min, source->uv_min, sizeof(item->uv_min));
        memcpy(item->uv_max, source->uv_max, sizeof(item->uv_max));
        item->kind = ALEA_TRANSITION_SLICE_TILE_SOURCE_REFINEMENT_FRONTIER;
        item->source_index = i;
        item->priority = 0;
    }
    const double padding = options->critical_tile_padding > 0.0
        ? options->critical_tile_padding : 0.0;
    for (size_t i = 0; i < seed_count; i++) {
        seeds[i].uv_min[0] = fmax(view->u_min, seeds[i].uv_min[0] - padding);
        seeds[i].uv_max[0] = fmin(view->u_max, seeds[i].uv_max[0] + padding);
        seeds[i].uv_min[1] = fmax(view->v_min, seeds[i].uv_min[1] - padding);
        seeds[i].uv_max[1] = fmin(view->v_max, seeds[i].uv_max[1] + padding);
    }
    qsort(seeds, seed_count, sizeof(*seeds),
          transition_slice_tile_seed_compare);
    size_t selected = seed_count;
    if (options->max_critical_tile_sources &&
        selected > options->max_critical_tile_sources) {
        selected = options->max_critical_tile_sources;
        result->stats.omitted_critical_tile_sources += seed_count - selected;
        transition_slice_critical_stop(
            result, ALEA_TRANSITION_SLICE_CRITICAL_MAX_TILE_SOURCES);
    }
    size_t merged_count = 0;
    for (size_t i = 0; i < selected; i++) {
        size_t target = SIZE_MAX;
        for (size_t tile = 0; tile < merged_count; tile++) {
            if (transition_slice_rectangles_touch(
                    merged[tile].uv_min, merged[tile].uv_max,
                    seeds[i].uv_min, seeds[i].uv_max)) {
                target = tile;
                break;
            }
        }
        if (target == SIZE_MAX) {
            if (options->max_critical_tiles &&
                merged_count >= options->max_critical_tiles) {
                transition_slice_critical_stop(
                    result, ALEA_TRANSITION_SLICE_CRITICAL_MAX_TILES);
                continue;
            }
            target = merged_count++;
            memcpy(merged[target].uv_min, seeds[i].uv_min,
                   sizeof(merged[target].uv_min));
            memcpy(merged[target].uv_max, seeds[i].uv_max,
                   sizeof(merged[target].uv_max));
        } else {
            for (int axis = 0; axis < 2; axis++) {
                if (merged[target].uv_min[axis] > seeds[i].uv_min[axis])
                    merged[target].uv_min[axis] = seeds[i].uv_min[axis];
                if (merged[target].uv_max[axis] < seeds[i].uv_max[axis])
                    merged[target].uv_max[axis] = seeds[i].uv_max[axis];
            }
        }
        merged[target].source_flags |= 1u << (unsigned)seeds[i].kind;
        for (size_t tile = 0; tile < merged_count;) {
            if (tile == target || !transition_slice_rectangles_touch(
                    merged[target].uv_min, merged[target].uv_max,
                    merged[tile].uv_min, merged[tile].uv_max)) {
                tile++;
                continue;
            }
            for (int axis = 0; axis < 2; axis++) {
                if (merged[target].uv_min[axis] > merged[tile].uv_min[axis])
                    merged[target].uv_min[axis] = merged[tile].uv_min[axis];
                if (merged[target].uv_max[axis] < merged[tile].uv_max[axis])
                    merged[target].uv_max[axis] = merged[tile].uv_max[axis];
            }
            merged[target].source_flags |= merged[tile].source_flags;
            if (tile < target) target--;
            memmove(&merged[tile], &merged[tile + 1],
                    (merged_count - tile - 1) * sizeof(*merged));
            merged_count--;
            tile = 0;
        }
    }
    size_t retained = transition_slice_retained_bytes(result);
    size_t tile_capacity = merged_count;
    if (options->max_output_bytes) {
        size_t available = retained < options->max_output_bytes
            ? options->max_output_bytes - retained : 0;
        size_t cap = available / sizeof(*result->critical_tiles);
        if (tile_capacity > cap) tile_capacity = cap;
    }
    if (tile_capacity < merged_count)
        transition_slice_critical_stop(
            result, ALEA_TRANSITION_SLICE_CRITICAL_MAX_OUTPUT_BYTES);
    if (tile_capacity) {
        result->critical_tiles = calloc(
            tile_capacity, sizeof(*result->critical_tiles));
        if (!result->critical_tiles) { free(seeds); free(merged); return -1; }
        memcpy(result->critical_tiles, merged,
               tile_capacity * sizeof(*result->critical_tiles));
        result->critical_tile_count = tile_capacity;
    }
    retained = transition_slice_retained_bytes(result);
    size_t source_capacity = selected;
    if (options->max_output_bytes) {
        size_t available = retained < options->max_output_bytes
            ? options->max_output_bytes - retained : 0;
        size_t cap = available / sizeof(*result->critical_tile_sources);
        if (source_capacity > cap) source_capacity = cap;
    }
    if (source_capacity) {
        result->critical_tile_sources = calloc(
            source_capacity, sizeof(*result->critical_tile_sources));
        if (!result->critical_tile_sources) {
            free(seeds); free(merged); return -1;
        }
    }
    for (size_t tile = 0; tile < result->critical_tile_count; tile++) {
        result->critical_tiles[tile].first_source_index =
            result->critical_tile_source_count;
        for (size_t i = 0; i < selected; i++) {
            if (result->critical_tile_source_count >= source_capacity) break;
            if (seeds[i].uv_min[0] < merged[tile].uv_min[0] ||
                seeds[i].uv_max[0] > merged[tile].uv_max[0] ||
                seeds[i].uv_min[1] < merged[tile].uv_min[1] ||
                seeds[i].uv_max[1] > merged[tile].uv_max[1]) continue;
            result->critical_tile_sources[result->critical_tile_source_count++] =
                (alea_transition_slice_critical_tile_source_t){
                    .kind = seeds[i].kind,
                    .source_index = seeds[i].source_index
                };
            result->critical_tiles[tile].source_count++;
        }
    }
    if (result->critical_tile_source_count < selected) {
        result->stats.omitted_critical_tile_sources +=
            selected - result->critical_tile_source_count;
        transition_slice_critical_stop(
            result, ALEA_TRANSITION_SLICE_CRITICAL_MAX_OUTPUT_BYTES);
    }
    result->stats.critical_tiles = result->critical_tile_count;
    result->stats.critical_tile_sources = result->critical_tile_source_count;
    result->stats.retained_output_bytes = transition_slice_retained_bytes(result);
    free(seeds);
    free(merged);
    return 0;
}

static int transition_slice_scan_orientation(
    alea_system_t* sys, const alea_slice_view_t* view,
    const alea_transition_slice_options_t* options,
    alea_transition_slice_orientation_t orientation,
    size_t reserved_base_rays,
    alea_raycast_result_t* scratch,
    alea_ray_boundary_event_result_t* events,
    transition_slice_coverage_scratch_t* coverage_scratch,
    alea_transition_workspace_t* transition_workspace,
    size_t coverage_scratch_bytes,
    alea_transition_slice_result_t* result) {
    const size_t base_count = orientation == ALEA_TRANSITION_SLICE_HORIZONTAL
        ? options->horizontal_rays : options->vertical_rays;
    if (base_count == 0) return 0;
    if (base_count > SIZE_MAX / sizeof(transition_slice_row_t)) return -1;
    size_t row_bytes = base_count * sizeof(transition_slice_row_t);
    size_t base_live_bytes = row_bytes > SIZE_MAX - coverage_scratch_bytes
        ? SIZE_MAX : row_bytes + coverage_scratch_bytes;
    if (options->max_row_scratch_bytes &&
        base_live_bytes > options->max_row_scratch_bytes) {
        result->stats.stop_reason =
            ALEA_TRANSITION_SLICE_STOP_MAX_SCRATCH_BYTES;
        return 1;
    }
    transition_slice_row_t* rows = calloc(base_count, sizeof(*rows));
    if (!rows) return -1;
    transition_slice_record_row_scratch(result, base_live_bytes);
    const double low = orientation == ALEA_TRANSITION_SLICE_HORIZONTAL
        ? view->v_min : view->u_min;
    const double span = orientation == ALEA_TRANSITION_SLICE_HORIZONTAL
        ? view->v_max - view->v_min : view->u_max - view->u_min;
    for (size_t base = 0; base < base_count; base++) {
        if (alea_interrupted()) {
            result->stats.stop_reason = ALEA_TRANSITION_SLICE_STOP_INTERRUPTED;
            free(rows);
            return 1;
        }
        const double coordinate = low + ((double)base + 0.5) * span /
            (double)base_count;
        const size_t ray_index = orientation == ALEA_TRANSITION_SLICE_HORIZONTAL
            ? result->stats.horizontal_rays_executed
            : result->stats.vertical_rays_executed;
        int rc = transition_slice_scan_ray_reuse(
            sys, view, options, orientation, ray_index, base, 0, coordinate,
            scratch, events, coverage_scratch, transition_workspace, result,
            &rows[base]);
        transition_slice_record_row_scratch(result, base_live_bytes);
        if (rc != 0) { free(rows); return rc; }
    }
    if (options->max_refinement_depth == 0) {
        if (transition_slice_retain_differing_frontiers(
                view, options, orientation, rows, base_count, result) != 0) {
            free(rows);
            return -1;
        }
        free(rows);
        return 0;
    }

    size_t row_count = base_count;
    for (uint32_t depth = 0;; depth++) {
        size_t marked = 0, spacing_limited = 0;
        for (size_t row = 0; row + 1 < row_count; row++) {
            if (!transition_slice_rows_differ(
                    &rows[row], &rows[row + 1], options->refine_signals))
                continue;
            const double spacing = rows[row + 1].transverse_coordinate -
                rows[row].transverse_coordinate;
            if (options->min_transverse_spacing > 0.0 &&
                spacing < 2.0 * options->min_transverse_spacing) {
                spacing_limited++;
                continue;
            }
            marked++;
        }
        if (marked == 0) {
            if (spacing_limited &&
                transition_slice_retain_differing_frontiers(
                    view, options, orientation, rows, row_count,
                    result) != 0) {
                free(rows);
                return -1;
            }
            transition_slice_merge_status(
                result, spacing_limited
                    ? ALEA_TRANSITION_SLICE_REFINEMENT_MIN_SPACING
                    : ALEA_TRANSITION_SLICE_REFINEMENT_CONVERGED);
            free(rows);
            return 0;
        }
        if (depth >= options->max_refinement_depth) {
            if (transition_slice_retain_differing_frontiers(
                    view, options, orientation, rows, row_count,
                    result) != 0) {
                free(rows);
                return -1;
            }
            transition_slice_merge_status(
                result, ALEA_TRANSITION_SLICE_REFINEMENT_MAX_DEPTH);
            free(rows);
            return 0;
        }
        uint64_t available_rays = options->max_rays &&
                result->stats.executed_rays + reserved_base_rays <
                    options->max_rays
            ? options->max_rays - result->stats.executed_rays -
                reserved_base_rays
            : 0;
        if (options->max_rays && marked > available_rays) {
            if (transition_slice_retain_differing_frontiers(
                    view, options, orientation, rows, row_count,
                    result) != 0) {
                free(rows);
                return -1;
            }
            result->stats.stop_reason = ALEA_TRANSITION_SLICE_STOP_MAX_RAYS;
            free(rows);
            return 1;
        }
        if (marked > SIZE_MAX - row_count) { free(rows); return -1; }
        size_t next_count = row_count + marked;
        if (next_count > SIZE_MAX / sizeof(*rows)) { free(rows); return -1; }
        size_t next_bytes = next_count * sizeof(*rows);
        size_t live_bytes = row_bytes > SIZE_MAX - next_bytes
            ? SIZE_MAX : row_bytes + next_bytes;
        live_bytes = live_bytes > SIZE_MAX - coverage_scratch_bytes
            ? SIZE_MAX : live_bytes + coverage_scratch_bytes;
        if (options->max_row_scratch_bytes &&
            live_bytes > options->max_row_scratch_bytes) {
            if (transition_slice_retain_differing_frontiers(
                    view, options, orientation, rows, row_count,
                    result) != 0) {
                free(rows);
                return -1;
            }
            result->stats.stop_reason =
                ALEA_TRANSITION_SLICE_STOP_MAX_SCRATCH_BYTES;
            free(rows);
            return 1;
        }
        transition_slice_row_t* next = calloc(next_count, sizeof(*next));
        if (!next) { free(rows); return -1; }
        transition_slice_record_row_scratch(result, live_bytes);
        size_t output = 0;
        for (size_t row = 0; row < row_count; row++) {
            next[output++] = rows[row];
            if (row + 1 == row_count ||
                !transition_slice_rows_differ(
                    &rows[row], &rows[row + 1], options->refine_signals))
                continue;
            const double spacing = rows[row + 1].transverse_coordinate -
                rows[row].transverse_coordinate;
            if (options->min_transverse_spacing > 0.0 &&
                spacing < 2.0 * options->min_transverse_spacing)
                continue;
            const double coordinate = 0.5 *
                (rows[row].transverse_coordinate +
                 rows[row + 1].transverse_coordinate);
            const size_t ray_index =
                orientation == ALEA_TRANSITION_SLICE_HORIZONTAL
                    ? result->stats.horizontal_rays_executed
                    : result->stats.vertical_rays_executed;
            int rc = transition_slice_scan_ray_reuse(
                sys, view, options, orientation, ray_index, SIZE_MAX,
                depth + 1, coordinate, scratch, events, coverage_scratch,
                transition_workspace, result, &next[output]);
            transition_slice_record_row_scratch(result, live_bytes);
            if (rc != 0) {
                free(next); free(rows); return rc;
            }
            result->stats.refined_rays_executed++;
            if (result->stats.max_refinement_depth_reached < depth + 1)
                result->stats.max_refinement_depth_reached = depth + 1;
            output++;
        }
        free(rows);
        rows = next;
        row_count = next_count;
        row_bytes = next_bytes;
    }
}

int alea_transition_slice_screen(
    alea_system_t* sys, const alea_slice_view_t* view,
    const alea_transition_slice_options_t* input,
    alea_transition_slice_result_t* result) {
    if (!sys || !view || !result) return -1;
    alea_transition_slice_options_t defaults, options;
    alea_transition_slice_options_init(&defaults);
    options = defaults;
    if (input) {
        if (input->struct_size < sizeof(input->struct_size)) return -1;
        size_t bytes = input->struct_size < sizeof(options)
            ? input->struct_size : sizeof(options);
        memcpy(&options, input, bytes);
    }
    if (!(view->u_max > view->u_min) || !(view->v_max > view->v_min) ||
        options.max_coverage_hits == 0 ||
        options.max_coverage_hits > 16384 ||
        options.horizontal_rays > SIZE_MAX - options.vertical_rays ||
        !isfinite(options.min_transverse_spacing) ||
        options.min_transverse_spacing < 0.0 ||
        !isfinite(options.critical_tile_padding) ||
        options.critical_tile_padding < 0.0 ||
        options.occurrence_discovery <
            ALEA_TRANSITION_SLICE_OCCURRENCE_SAMPLED ||
        options.occurrence_discovery >
            ALEA_TRANSITION_SLICE_OCCURRENCE_EXHAUSTIVE ||
        (options.occurrence_discovery ==
             ALEA_TRANSITION_SLICE_OCCURRENCE_EXHAUSTIVE &&
         options.max_exhaustive_occurrence_hits == 0) ||
        (options.refine_signals &
         ~(ALEA_TRANSITION_SLICE_REFINE_SIGNATURE |
           ALEA_TRANSITION_SLICE_REFINE_FINDING)))
        return -1;

    alea_transition_slice_result_t candidate = {0};
    candidate.stats.occurrence_discovery = options.occurrence_discovery;
    candidate.stats.critical_enabled = options.enable_critical_refinement != 0;
    candidate.stats.critical_stop_reason = options.enable_critical_refinement
        ? ALEA_TRANSITION_SLICE_CRITICAL_NONE
        : ALEA_TRANSITION_SLICE_CRITICAL_DISABLED;
    candidate.stats.requested_rays =
        options.horizontal_rays + options.vertical_rays;
    candidate.stats.refinement_status = options.max_refinement_depth
        ? ALEA_TRANSITION_SLICE_REFINEMENT_CONVERGED
        : ALEA_TRANSITION_SLICE_REFINEMENT_NOT_REQUESTED;
    int scan_incomplete = 0;
    if (options.max_rays && candidate.stats.requested_rays > options.max_rays) {
        candidate.stats.stop_reason = ALEA_TRANSITION_SLICE_STOP_MAX_RAYS;
        scan_incomplete = 1;
    } else if (alea_raycast_ensure_hier_caches(sys) != 0) {
        return -1;
    } else {
        alea_raycast_result_t scratch;
        alea_ray_boundary_event_result_t events;
        transition_slice_coverage_scratch_t coverage_scratch;
        alea_transition_workspace_t transition_workspace;
        memset(&coverage_scratch, 0, sizeof(coverage_scratch));
        const int coverage_enabled =
            options.coverage_uniform_probes_per_ray != 0 ||
            options.coverage_probe_selected_intervals;
        size_t coverage_scratch_bytes = coverage_enabled
            ? transition_slice_coverage_scratch_bytes(
                options.max_coverage_hits) : 0;
        alea_raycast_result_init(&scratch);
        alea_ray_boundary_event_result_init(&events);
        alea_transition_workspace_init(&transition_workspace);
        int stop = 0, failed = 0;
        if (coverage_scratch_bytes == SIZE_MAX ||
            (options.max_row_scratch_bytes && coverage_enabled &&
             coverage_scratch_bytes > options.max_row_scratch_bytes)) {
            candidate.stats.stop_reason =
                ALEA_TRANSITION_SLICE_STOP_MAX_SCRATCH_BYTES;
            stop = 1;
            scan_incomplete = 1;
        } else if (coverage_enabled &&
                   transition_slice_coverage_scratch_init(
                       &coverage_scratch, options.max_coverage_hits) != 0) {
            failed = 1;
        }
        for (int orientation = ALEA_TRANSITION_SLICE_HORIZONTAL;
             orientation <= ALEA_TRANSITION_SLICE_VERTICAL && !stop && !failed;
             orientation++) {
            int rc = transition_slice_scan_orientation(
                sys, view, &options,
                (alea_transition_slice_orientation_t)orientation,
                orientation == ALEA_TRANSITION_SLICE_HORIZONTAL
                    ? options.vertical_rays : 0,
                &scratch, &events, &coverage_scratch, &transition_workspace,
                coverage_scratch_bytes, &candidate);
            if (rc < 0) failed = 1;
            if (rc > 0) { stop = 1; scan_incomplete = 1; }
            if (failed) break;
        }
        alea_ray_boundary_event_result_free(&events);
        alea_raycast_result_free(&scratch);
        transition_slice_coverage_scratch_free(&coverage_scratch);
        alea_transition_workspace_free(&transition_workspace);
        if (!failed) {
            int component_rc = transition_slice_build_components(
                &options, &candidate);
            if (component_rc < 0) failed = 1;
            if (component_rc > 0) stop = 1;
        }
        if (!failed) {
            int component_rc = transition_slice_build_coverage_components(
                &options, &candidate);
            if (component_rc < 0) failed = 1;
            if (component_rc > 0) stop = 1;
        }
        if (!failed) {
            int component_rc = transition_slice_build_component_links(
                &options, &candidate);
            if (component_rc < 0) failed = 1;
            if (component_rc > 0) stop = 1;
        }
        if (!failed && transition_slice_build_critical_tiles(
                view, &options, &candidate) != 0)
            failed = 1;
        if (!failed && options.enable_critical_refinement) {
            transition_slice_critical_sink_context_t critical_sink = {
                .result = &candidate, .options = &options
            };
            if (alea_transition_slice_enumerate_critical_tiles(
                    sys, view, &options, candidate.critical_tiles,
                    candidate.critical_tile_count,
                    transition_slice_append_critical_finding,
                    &critical_sink, &candidate.stats) != 0)
                failed = 1;
        }
        if (failed) {
            free(candidate.findings);
            free(candidate.components);
            free(candidate.coverage_findings);
            free(candidate.coverage_components);
            free(candidate.component_links);
            free(candidate.refinement_frontiers);
            free(candidate.critical_tiles);
            free(candidate.critical_tile_sources);
            free(candidate.critical_findings);
            return -1;
        }
    }
    candidate.stats.critical_complete = options.enable_critical_refinement &&
        candidate.stats.critical_stop_reason ==
            ALEA_TRANSITION_SLICE_CRITICAL_NONE;
    candidate.stats.occurrence_enumeration_complete =
        options.occurrence_discovery ==
            ALEA_TRANSITION_SLICE_OCCURRENCE_EXHAUSTIVE &&
        candidate.stats.critical_complete;
    candidate.stats.complete =
        candidate.stats.stop_reason == ALEA_TRANSITION_SLICE_STOP_NONE &&
        (!options.enable_critical_refinement ||
         candidate.stats.critical_complete);
    if (scan_incomplete && options.max_refinement_depth) {
        candidate.stats.refinement_status =
            ALEA_TRANSITION_SLICE_REFINEMENT_STOPPED;
    }
    candidate.stats.converged = candidate.stats.refinement_status ==
            ALEA_TRANSITION_SLICE_REFINEMENT_CONVERGED;
    free(result->findings);
    free(result->components);
    free(result->coverage_findings);
    free(result->coverage_components);
    free(result->component_links);
    free(result->refinement_frontiers);
    free(result->critical_tiles);
    free(result->critical_tile_sources);
    free(result->critical_findings);
    *result = candidate;
    return 0;
}

typedef struct {
    alea_system_t* sys;
    const alea_slice_view_t* views;
    const alea_transition_slice_options_t* options;
    alea_transition_slice_result_t* const* results;
    int* statuses;
} alea_transition_slice_batch_context_t;

static int alea_transition_slice_batch_run(
    void* opaque, size_t worker_index, size_t begin, size_t end) {
    (void)worker_index;
    alea_transition_slice_batch_context_t* context = opaque;
    for (size_t page = begin; page < end; page++) {
        context->statuses[page] = alea_transition_slice_screen(
            context->sys, &context->views[page], context->options,
            context->results[page]);
    }
    return 0;
}

int alea_transition_slice_screen_batch(
    alea_system_t* sys, const alea_slice_view_t* views, size_t page_count,
    const alea_transition_slice_options_t* input,
    size_t requested_workers, uint64_t max_parallel_scratch_bytes,
    alea_transition_slice_result_t* const* results,
    alea_transition_slice_batch_stats_t* out_stats) {
    if (!sys || (!views && page_count) || (!results && page_count) ||
            !out_stats) {
        alea_set_error_detail(ALEA_ERR_NULL_ARG,
                              "transition slice batch argument is null");
        return -1;
    }
    memset(out_stats, 0, sizeof(*out_stats));
    out_stats->page_count = page_count;
    out_stats->requested_workers = requested_workers;
    if (page_count == 0) return 0;
    for (size_t page = 0; page < page_count; page++) {
        if (!results[page]) {
            alea_set_error_detail(ALEA_ERR_NULL_ARG,
                                  "transition slice batch result is null");
            return -1;
        }
    }

    alea_transition_slice_options_t defaults, options;
    alea_transition_slice_options_init(&defaults);
    options = defaults;
    if (input) {
        if (input->struct_size < sizeof(input->struct_size)) {
            alea_set_error_detail(ALEA_ERR_INVALID_ARG,
                                  "transition slice options are too small");
            return -1;
        }
        size_t bytes = input->struct_size < sizeof(options)
            ? input->struct_size : sizeof(options);
        memcpy(&options, input, bytes);
    }
    uint64_t worker_scratch = options.max_scratch_bytes;
    if (worker_scratch > UINT64_MAX - options.max_critical_scratch_bytes) {
        alea_set_error_detail(ALEA_ERR_OVERFLOW,
                              "transition slice worker scratch overflows");
        return -1;
    }
    worker_scratch += options.max_critical_scratch_bytes;
    if (worker_scratch == 0) worker_scratch = 1;
    out_stats->reserved_scratch_bytes_per_worker = worker_scratch;

    size_t workers = alea_parallel_effective_workers(page_count, 1,
                                                     requested_workers);
    if (max_parallel_scratch_bytes == 0) {
        workers = 1;
    } else {
        uint64_t budget_workers = max_parallel_scratch_bytes / worker_scratch;
        if (budget_workers == 0) budget_workers = 1;
        if ((uint64_t)workers > budget_workers)
            workers = (size_t)budget_workers;
    }
    if (workers > UINT64_MAX / worker_scratch) {
        alea_set_error_detail(ALEA_ERR_OVERFLOW,
                              "transition slice parallel scratch overflows");
        return -1;
    }
    out_stats->reserved_parallel_scratch_bytes =
        (uint64_t)workers * worker_scratch;

    if (alea_raycast_ensure_hier_caches(sys) != 0) return -1;
    int* statuses = calloc(page_count, sizeof(*statuses));
    if (!statuses) {
        alea_set_error_detail(ALEA_ERR_OUT_OF_MEMORY,
                              "transition slice batch status allocation failed");
        return -1;
    }
    size_t actual_workers = 1;
    alea_transition_slice_batch_context_t context = {
        sys, views, &options, results, statuses
    };
    alea_parallel_status_t parallel_status = alea_parallel_for(
        page_count, 1, workers, ALEA_PARALLEL_STATIC_BLOCK,
        alea_transition_slice_batch_run, &context, &actual_workers);
    if (parallel_status != ALEA_PARALLEL_OK) {
        free(statuses);
        alea_set_error_detail(ALEA_ERR_INVALID_STATE,
                              alea_parallel_status_string(parallel_status));
        return -1;
    }
    out_stats->actual_workers = actual_workers;
    int failed = 0;
    for (size_t page = 0; page < page_count; page++) {
        if (statuses[page] == 0) out_stats->completed_page_count++;
        else failed = 1;
    }
    free(statuses);
    return failed ? -1 : 0;
}
