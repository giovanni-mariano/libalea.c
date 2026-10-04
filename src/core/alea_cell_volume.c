// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/* Read-only, universe-local cell volume estimation. */

#include "alea.h"
#include "core/alea_system.h"
#include "util/alea_parallel.h"
#include "core/alea_eval.h"
#include "primitives/bbox.h"

#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#define CELL_VOLUME_UNBOUNDED_EXTENT 9e5
#define CELL_VOLUME_DISCOVERY_DEPTH 6
#define CELL_VOLUME_DISCOVERY_EXPANSIONS 24
#define CELL_VOLUME_BATCH_PARENTS 64
#define CELL_VOLUME_BATCH_CHILDREN (8 * CELL_VOLUME_BATCH_PARENTS)

typedef struct {
    alea_bbox_t bbox;
    int level;
    int axis_depth[3];
    unsigned active_axes;
    double interval_width;
} cell_volume_task_t;

typedef struct {
    unsigned active_axes;
    uint8_t relation; /* 0 outside, 1 inside, 2 mixed */
    double interval_width;
} cell_volume_classification_t;

typedef struct {
    alea_bbox_t envelope;
    bool has_content;
    bool touches_boundary;
} bbox_discovery_scan_t;

static bool cell_volume_bbox_valid(const alea_bbox_t* b) {
    return b && isfinite(b->min_x) && isfinite(b->max_x) &&
        isfinite(b->min_y) && isfinite(b->max_y) &&
        isfinite(b->min_z) && isfinite(b->max_z) &&
        b->min_x < b->max_x && b->min_y < b->max_y &&
        b->min_z < b->max_z;
}

static bool cell_volume_bbox_finite_stored(const alea_bbox_t* b) {
    if (!cell_volume_bbox_valid(b)) return false;
    return b->max_x - b->min_x <= CELL_VOLUME_UNBOUNDED_EXTENT &&
           b->max_y - b->min_y <= CELL_VOLUME_UNBOUNDED_EXTENT &&
           b->max_z - b->min_z <= CELL_VOLUME_UNBOUNDED_EXTENT;
}

static bool cell_volume_bbox_contains(const alea_bbox_t* outer,
                                      const alea_bbox_t* inner) {
    return outer->min_x <= inner->min_x && outer->max_x >= inner->max_x &&
           outer->min_y <= inner->min_y && outer->max_y >= inner->max_y &&
           outer->min_z <= inner->min_z && outer->max_z >= inner->max_z;
}

static double cell_volume_box_volume(const alea_bbox_t* b) {
    return (b->max_x - b->min_x) * (b->max_y - b->min_y) *
           (b->max_z - b->min_z);
}

static double cell_volume_min_extent(const alea_bbox_t* b) {
    double x = b->max_x - b->min_x;
    double y = b->max_y - b->min_y;
    double z = b->max_z - b->min_z;
    return fmin(x, fmin(y, z));
}

static alea_bbox_t cell_volume_child_bbox(const alea_bbox_t* b, int child) {
    double mx = (b->min_x + b->max_x) * 0.5;
    double my = (b->min_y + b->max_y) * 0.5;
    double mz = (b->min_z + b->max_z) * 0.5;
    return (alea_bbox_t){
        (child & 1) ? mx : b->min_x, (child & 1) ? b->max_x : mx,
        (child & 2) ? my : b->min_y, (child & 2) ? b->max_y : my,
        (child & 4) ? mz : b->min_z, (child & 4) ? b->max_z : mz,
    };
}

/* For volume, boxes whose interval merely touches zero have zero-measure
 * ambiguity and can be accepted. This makes an axis-aligned box exact at its
 * own bbox while preserving conservative volume bounds. */
static uint8_t cell_volume_relation(const alea_system_t* sys,
                                    alea_node_id_t root,
                                    const alea_bbox_t* bbox) {
    alea_interval_t iv = alea_evaluate_interval(sys, root, bbox);
    if (iv.max <= 0.0) return 1;
    if (iv.min >= 0.0) return 0;
    return 2;
}

static double cell_volume_sample_fraction(const alea_system_t* sys,
                                          alea_node_id_t root,
                                          const alea_bbox_t* b, int n) {
    size_t inside = 0;
    size_t total = (size_t)n * (size_t)n * (size_t)n;
    double dx = (b->max_x - b->min_x) / (double)n;
    double dy = (b->max_y - b->min_y) / (double)n;
    double dz = (b->max_z - b->min_z) / (double)n;
    for (int k = 0; k < n; k++) {
        if (alea_interrupted()) return 0.0;
        double z = b->min_z + ((double)k + 0.5) * dz;
        for (int j = 0; j < n; j++) {
            double y = b->min_y + ((double)j + 0.5) * dy;
            for (int i = 0; i < n; i++) {
                double x = b->min_x + ((double)i + 0.5) * dx;
                if (alea_contains_point(sys, root, x, y, z)) inside++;
            }
        }
    }
    return total ? (double)inside / (double)total : 0.0;
}

static bool discovery_touches(const alea_bbox_t* leaf,
                              const alea_bbox_t* domain) {
    double eps = 1e-12 * fmax(1.0, cell_volume_min_extent(domain));
    return leaf->min_x <= domain->min_x + eps ||
           leaf->max_x >= domain->max_x - eps ||
           leaf->min_y <= domain->min_y + eps ||
           leaf->max_y >= domain->max_y - eps ||
           leaf->min_z <= domain->min_z + eps ||
           leaf->max_z >= domain->max_z - eps;
}

static void discovery_add_leaf(bbox_discovery_scan_t* scan,
                               const alea_bbox_t* leaf,
                               const alea_bbox_t* domain) {
    if (!scan->has_content) {
        scan->envelope = *leaf;
        scan->has_content = true;
    } else {
        scan->envelope = alea_bbox_union(&scan->envelope, leaf);
    }
    if (discovery_touches(leaf, domain)) scan->touches_boundary = true;
}

static void discovery_scan_box(const alea_system_t* sys, alea_node_id_t root,
                               const alea_bbox_t* box,
                               const alea_bbox_t* domain, int depth,
                               bbox_discovery_scan_t* scan) {
    uint8_t relation = cell_volume_relation(sys, root, box);
    if (relation == 0) return;
    if (relation == 1 || depth == 0) {
        discovery_add_leaf(scan, box, domain);
        return;
    }
    for (int child = 0; child < 8; child++) {
        alea_bbox_t next = cell_volume_child_bbox(box, child);
        discovery_scan_box(sys, root, &next, domain, depth - 1, scan);
    }
}

static bool discovery_envelopes_stable(const alea_bbox_t* a,
                                       const alea_bbox_t* b,
                                       double tolerance) {
    return fabs(a->min_x - b->min_x) <= tolerance &&
           fabs(a->max_x - b->max_x) <= tolerance &&
           fabs(a->min_y - b->min_y) <= tolerance &&
           fabs(a->max_y - b->max_y) <= tolerance &&
           fabs(a->min_z - b->min_z) <= tolerance &&
           fabs(a->max_z - b->max_z) <= tolerance;
}

static int cell_volume_discover_bbox(const alea_system_t* sys,
                                     alea_node_id_t root,
                                     const alea_bbox_t* stored,
                                     alea_bbox_t* out,
                                     alea_cell_volume_bounds_source_t* source,
                                     size_t* expansions) {
    double tol = 1e-6;
    if (alea_tighten_bbox_plane_constraints(sys, root, tol, out) == 0 &&
        cell_volume_bbox_valid(out)) {
        *source = ALEA_CELL_VOLUME_BOUNDS_PLANE_CONSTRAINTS;
        *expansions = 0;
        return 0;
    }

    double scale = 1.0;
    const double values[6] = {stored->min_x, stored->max_x, stored->min_y,
                              stored->max_y, stored->min_z, stored->max_z};
    for (int i = 0; i < 6; i++) {
        if (isfinite(values[i]) && fabs(values[i]) < CELL_VOLUME_UNBOUNDED_EXTENT)
            scale = fmax(scale, fabs(values[i]));
    }

    bool have_isolated = false;
    alea_bbox_t isolated = {0};
    for (size_t expansion = 0; expansion < CELL_VOLUME_DISCOVERY_EXPANSIONS;
         expansion++) {
        if (!isfinite(scale) || scale > 1e15) break;
        alea_bbox_t domain = {-scale, scale, -scale, scale, -scale, scale};
        bbox_discovery_scan_t scan = {0};
        discovery_scan_box(sys, root, &domain, &domain,
                           CELL_VOLUME_DISCOVERY_DEPTH, &scan);
        if (scan.has_content && !scan.touches_boundary) {
            double leaf = (2.0 * scale) /
                (double)(1u << CELL_VOLUME_DISCOVERY_DEPTH);
            if (have_isolated &&
                discovery_envelopes_stable(&isolated, &scan.envelope,
                                           2.0 * leaf)) {
                alea_bbox_t candidate = alea_bbox_union(&isolated, &scan.envelope);
                alea_tighten_tree_bbox(sys, root, &candidate,
                                       fmax(1e-9, leaf * 0.01), out);
                if (!cell_volume_bbox_valid(out)) return -1;
                *source = ALEA_CELL_VOLUME_BOUNDS_ADAPTIVE_SEARCH;
                *expansions = expansion;
                return 0;
            }
            isolated = scan.envelope;
            have_isolated = true;
        } else {
            have_isolated = false;
        }
        scale *= 2.0;
    }
    return -1;
}

static size_t cell_volume_select_workers(size_t requested, size_t tasks,
                                         uint64_t budget,
                                         uint64_t scratch_per_worker) {
    size_t workers = requested;
    if (alea_parallel_in_region()) return 1;
    if (workers == 0) workers = alea_parallel_max_workers();
    if (budget == 0 || tasks < 2) return 1;
    if (workers < 1) workers = 1;
    if (workers > tasks) workers = tasks;
    /* Small fixed batches should not wake hundreds of backend threads. */
    size_t useful_workers = tasks / 64;
    if (useful_workers < 1) useful_workers = 1;
    if (workers > useful_workers) workers = useful_workers;
    if (scratch_per_worker != 0) {
        uint64_t by_budget = budget / scratch_per_worker;
        if (by_budget == 0) return 1;
        if ((uint64_t)workers > by_budget) workers = (size_t)by_budget;
    }
    return workers ? workers : 1;
}

/* The heap prefers the largest unresolved volume, then uses box coordinates
 * for a stable order independent of worker scheduling. */
static bool task_precedes(const cell_volume_task_t* a, const cell_volume_task_t* b) {
    double va = cell_volume_box_volume(&a->bbox), vb = cell_volume_box_volume(&b->bbox);
    if (va != vb) return va > vb;
    const double ca[] = {a->bbox.min_x, a->bbox.min_y, a->bbox.min_z,
                         a->bbox.max_x, a->bbox.max_y, a->bbox.max_z};
    const double cb[] = {b->bbox.min_x, b->bbox.min_y, b->bbox.min_z,
                         b->bbox.max_x, b->bbox.max_y, b->bbox.max_z};
    for (int i = 0; i < 6; i++) if (ca[i] != cb[i]) return ca[i] < cb[i];
    return false;
}

static void heap_push(cell_volume_task_t* heap, size_t* count, cell_volume_task_t t) {
    size_t i = (*count)++;
    while (i && task_precedes(&t, &heap[(i-1)/2])) {
        heap[i] = heap[(i-1)/2]; i = (i-1)/2;
    }
    heap[i] = t;
}

static cell_volume_task_t heap_pop(cell_volume_task_t* heap, size_t* count) {
    cell_volume_task_t result = heap[0], tail = heap[--(*count)];
    if (*count == 0) return result;
    size_t i = 0;
    while (2*i+1 < *count) {
        size_t j = 2*i+1;
        if (j+1 < *count && task_precedes(&heap[j+1], &heap[j])) j++;
        if (!task_precedes(&heap[j], &tail)) break;
        heap[i] = heap[j]; i = j;
    }
    heap[i] = tail;
    return result;
}

static double axis_min(const alea_bbox_t* b, int a) {
    return a == 0 ? b->min_x : a == 1 ? b->min_y : b->min_z;
}
static double axis_max(const alea_bbox_t* b, int a) {
    return a == 0 ? b->max_x : a == 1 ? b->max_y : b->max_z;
}
static double axis_mid(const alea_bbox_t* b, int a) {
    /* Avoid overflow when finite endpoints have the same large sign. */
    return axis_min(b, a) * 0.5 + axis_max(b, a) * 0.5;
}
static alea_bbox_t split_bbox(alea_bbox_t b, int mask, int child) {
    for (int a = 0; a < 3; a++) if (mask & (1 << a)) {
        double mid = axis_mid(&b, a);
        if (a == 0) { if (child & 1) b.min_x = mid; else b.max_x = mid; }
        if (a == 1) { if (child & 2) b.min_y = mid; else b.max_y = mid; }
        if (a == 2) { if (child & 4) b.min_z = mid; else b.max_z = mid; }
    }
    return b;
}
static uint8_t interval_relation(alea_interval_t iv) {
    return iv.max <= 0.0 ? 1 : iv.min >= 0.0 ? 0 : 2;
}
static int eligible_axes(const cell_volume_task_t* t,
                         const alea_cell_volume_options_t* o) {
    int mask = 0;
    for (int a = 0; a < 3; a++) {
        double lo = axis_min(&t->bbox, a), hi = axis_max(&t->bbox, a);
        double mid = axis_mid(&t->bbox, a);
        if ((o->split_strategy != ALEA_CELL_VOLUME_SPLIT_ADAPTIVE ||
             (t->active_axes & (1u << a))) && t->axis_depth[a] < o->max_depth &&
            (o->min_size == 0.0 || hi-lo > o->min_size) && mid > lo && mid < hi)
            mask |= 1 << a;
    }
    /* Octree retains the existing rule: all axes must be eligible. */
    if (o->split_strategy == ALEA_CELL_VOLUME_SPLIT_OCTREE && mask != 7) return 0;
    return mask;
}

static int select_axes(const alea_system_t* sys, alea_node_id_t root,
                        const cell_volume_task_t* t, int eligible,
                        const alea_cell_volume_options_t* o,
                        alea_cell_volume_result_t* out) {
    if (o->split_strategy == ALEA_CELL_VOLUME_SPLIT_OCTREE) return 7;
    int longest = -1;
    for (int a = 0; a < 3; a++) if (eligible & (1 << a)) {
        if (longest < 0 || axis_max(&t->bbox,a)-axis_min(&t->bbox,a) >
                           axis_max(&t->bbox,longest)-axis_min(&t->bbox,longest))
            longest = a;
    }
    if (o->split_strategy == ALEA_CELL_VOLUME_SPLIT_LONGEST_AXIS) return 1 << longest;
    double scores[3] = {0}, best = 0.0;
    for (int a = 0; a < 3; a++) if (eligible & (1 << a)) {
        double width = 0.0, resolved = 0.0;
        for (int side = 0; side < 2; side++) {
            alea_bbox_t child = split_bbox(t->bbox, 1 << a, side << a);
            alea_interval_t iv = alea_evaluate_interval(sys, root, &child);
            out->interval_evaluations++;
            if (interval_relation(iv) != 2) resolved += 0.5;
            width += 0.5 * (iv.max - iv.min);
        }
        /* A split through a symmetric extremum (e.g. [-1,1]^2) can leave
         * both child widths unchanged. Probe the central half as well, so
         * such an axis is not starved by a greedy one-step score. */
        alea_bbox_t central = t->bbox;
        double lo = axis_min(&central, a), hi = axis_max(&central, a);
        double qlo = lo * 0.75 + hi * 0.25, qhi = lo * 0.25 + hi * 0.75;
        if (a == 0) { central.min_x = qlo; central.max_x = qhi; }
        if (a == 1) { central.min_y = qlo; central.max_y = qhi; }
        if (a == 2) { central.min_z = qlo; central.max_z = qhi; }
        alea_interval_t civ = alea_evaluate_interval(sys, root, &central);
        out->interval_evaluations++;
        width = fmin(width, civ.max - civ.min);
        double reduction = 0.0;
        if (isfinite(t->interval_width) && t->interval_width > 0.0 && isfinite(width))
            reduction = fmax(0.0, 1.0 - width / t->interval_width);
        scores[a] = resolved + reduction;
        if (scores[a] > best) best = scores[a];
    }
    if (best > 1e-12) {
        int mask = 0;
        /* Refine equally useful axes together (XY for an aligned cylinder,
         * XYZ for a symmetric sphere). Probes are included in the work budget. */
        for (int a = 0; a < 3; a++)
            if ((eligible & (1 << a)) && scores[a] >= best * 0.75) mask |= 1 << a;
        return mask;
    }
    /* Zero-gain ties are common before a boundary is resolved. Fair axis
     * counts prevent repeatedly bisecting one unproductive direction. */
    int fair = longest;
    for (int a = 0; a < 3; a++) if ((eligible & (1 << a)) &&
            t->axis_depth[a] < t->axis_depth[fair]) fair = a;
    return 1 << fair;
}

typedef struct {
    const alea_system_t* sys;
    alea_node_id_t root;
    const cell_volume_task_t* tasks;
    cell_volume_classification_t* classes;
} cell_volume_parallel_context_t;

static int cell_volume_parallel_range(void* opaque, size_t worker,
                                      size_t begin, size_t end) {
    cell_volume_parallel_context_t* c = opaque;
    (void)worker;
    for (size_t i = begin; i < end; i++) {
        alea_interval_t iv = alea_evaluate_interval_axes(c->sys, c->root,
            &c->tasks[i].bbox, &c->classes[i].active_axes);
        c->classes[i].relation = interval_relation(iv);
        c->classes[i].interval_width = iv.max - iv.min;
    }
    return 0;
}

static long double sample_pending(const alea_system_t* sys, alea_node_id_t root,
                                  const cell_volume_task_t* heap, size_t count, int n) {
    long double estimate = 0.0L;
    for (size_t i = 0; i < count; i++) {
        if ((i & 1023) == 0 && alea_interrupted()) break;
        estimate += (long double)cell_volume_sample_fraction(sys, root, &heap[i].bbox, n) *
                    cell_volume_box_volume(&heap[i].bbox);
    }
    return estimate;
}

void alea_cell_volume_options_init(alea_cell_volume_options_t* options) {
    if (!options) return;
    *options = (alea_cell_volume_options_t){
        .has_bounds = false,
        .bounds = {0},
        .relative_tolerance = 1e-3,
        .absolute_tolerance = 0.0,
        .max_depth = 14,
        .min_size = 0.0,
        .samples_per_axis = 2,
        .requested_workers = 0,
        .max_parallel_scratch_bytes = 64u * 1024u * 1024u,
        .split_strategy = ALEA_CELL_VOLUME_SPLIT_ADAPTIVE,
        .max_evaluations = 16000000,
        .max_memory_bytes = 64u * 1024u * 1024u,
    };
}

static int cell_volume_options_valid(const alea_cell_volume_options_t* o) {
    return o && isfinite(o->relative_tolerance) &&
        o->relative_tolerance >= 0.0 && isfinite(o->absolute_tolerance) &&
        o->absolute_tolerance >= 0.0 && o->max_depth >= 0 && o->max_depth <= INT_MAX / 3 &&
        isfinite(o->min_size) && o->min_size >= 0.0 &&
        o->samples_per_axis >= 1 && o->samples_per_axis <= 32 &&
        o->split_strategy >= ALEA_CELL_VOLUME_SPLIT_OCTREE &&
        o->split_strategy <= ALEA_CELL_VOLUME_SPLIT_ADAPTIVE &&
        o->max_evaluations >= 1 && o->max_memory_bytes >= 1024 &&
        (!o->has_bounds || cell_volume_bbox_valid(&o->bounds));
}

int alea_cell_estimate_volume(
        const alea_system_t* sys, size_t cell_index,
        const alea_cell_volume_options_t* options,
        alea_cell_volume_result_t* out) {
    if (!sys || !out || !cell_volume_options_valid(options) ||
        cell_index >= alea_vec_count(&sys->cells)) {
        alea_set_error_detail(ALEA_ERR_INVALID_ARG,
                              "invalid cell volume estimate arguments");
        return -1;
    }
    const alea_cell_entry_t* cell = &sys->cells.data[cell_index];
    if (cell->root_node_id == ALEA_NODE_ID_INVALID ||
        cell->root_node_id >= alea_vec_count(&sys->nodes)) {
        alea_set_error_detail(ALEA_ERR_INVALID_STATE,
                              "cell %zu has no valid CSG root", cell_index);
        return -1;
    }

    memset(out, 0, sizeof(*out));
    out->requested_workers = options->requested_workers;
    out->actual_workers = 1;
    out->scratch_bytes_per_worker = sizeof(cell_volume_classification_t);

    alea_bbox_t stored = alea_node_bbox_get(
        &sys->nodes.data[cell->root_node_id].bbox);
    alea_bbox_t geometric = alea_get_bbox(sys, cell->root_node_id);
    if (options->has_bounds) {
        out->bounds = options->bounds;
        out->bounds_source = ALEA_CELL_VOLUME_BOUNDS_EXPLICIT;
        out->complete_cell_domain = cell_volume_bbox_valid(&geometric) &&
            cell_volume_bbox_finite_stored(&geometric) &&
            cell_volume_bbox_contains(&options->bounds, &geometric);
    } else if (cell_volume_bbox_finite_stored(&stored)) {
        /* Node boxes use compact, conservatively padded storage. Recompute the
         * same CSG bound in double precision for the integration domain. */
        out->bounds = cell_volume_bbox_valid(&geometric) ? geometric : stored;
        out->bounds_source = ALEA_CELL_VOLUME_BOUNDS_STORED;
        out->complete_cell_domain = true;
    } else {
        if (cell_volume_discover_bbox(sys, cell->root_node_id, &stored,
                &out->bounds, &out->bounds_source,
                &out->bounds_search_expansions) != 0) {
            alea_set_error_detail(ALEA_ERR_INVALID_STATE,
                "cell %zu is unbounded or a finite bbox could not be isolated",
                cell_index);
            return -1;
        }
        out->complete_cell_domain =
            out->bounds_source != ALEA_CELL_VOLUME_BOUNDS_ADAPTIVE_SEARCH;
    }

    double root_volume = cell_volume_box_volume(&out->bounds);
    if (!isfinite(root_volume) || root_volume <= 0.0) {
        alea_set_error_detail(ALEA_ERR_OVERFLOW,
                              "cell volume integration bbox is not finite");
        return -1;
    }

    /* All explicit integration storage is charged to the memory budget.
     * Use fixed batch storage and a fixed-capacity heap: no hidden queue growth. */
    size_t batch_capacity = CELL_VOLUME_BATCH_CHILDREN;
    uint64_t per_batch = sizeof(cell_volume_task_t) + sizeof(cell_volume_classification_t);
    while (batch_capacity > 1 &&
           batch_capacity * per_batch + sizeof(cell_volume_task_t) > options->max_memory_bytes)
        batch_capacity /= 2;
    uint64_t batch_bytes = batch_capacity * per_batch;
    uint64_t capacity64 = (options->max_memory_bytes - batch_bytes) / sizeof(cell_volume_task_t);
    /* There can never be more queued leaves than classified boxes. */
    if (capacity64 > options->max_evaluations) capacity64 = options->max_evaluations;
    if (capacity64 > SIZE_MAX / sizeof(cell_volume_task_t))
        capacity64 = SIZE_MAX / sizeof(cell_volume_task_t);
    size_t capacity = (size_t)capacity64;
    cell_volume_task_t* heap = malloc(capacity * sizeof(*heap));
    cell_volume_task_t* children = malloc(batch_capacity * sizeof(*children));
    cell_volume_classification_t* classes = calloc(batch_capacity, sizeof(*classes));
    if (!heap || !children || !classes) {
        free(heap); free(children); free(classes);
        alea_set_error_detail(ALEA_ERR_OUT_OF_MEMORY, "failed to allocate cell volume storage");
        return -1;
    }
    out->peak_memory_bytes = capacity * sizeof(*heap) + batch_bytes;
    size_t count = 0;
    long double lower = 0.0L, terminal_gap = 0.0L, terminal_estimate = 0.0L;
    long double pending_gap = 0.0L;
    children[0] = (cell_volume_task_t){.bbox = out->bounds};
    size_t child_count = 1;
    bool stop = false;
    while (true) {
        if (alea_interrupted()) goto interrupted;
        size_t workers = cell_volume_select_workers(options->requested_workers,
            child_count, options->max_parallel_scratch_bytes, sizeof(*classes));
        size_t actual = 1;
        cell_volume_parallel_context_t context = {sys, cell->root_node_id, children, classes};
        if (alea_parallel_for(child_count, 1, workers, ALEA_PARALLEL_STATIC_BLOCK,
                              cell_volume_parallel_range, &context, &actual) != ALEA_PARALLEL_OK) {
            free(heap); free(children); free(classes);
            alea_set_error_detail(ALEA_ERR_INVALID_STATE, "cell volume parallel execution failed");
            return -1;
        }
        if (actual > out->actual_workers) out->actual_workers = actual;
        if (actual > 1) out->parallel_batch_count++;
        uint64_t scratch = workers * sizeof(*classes);
        if (scratch > out->reserved_parallel_scratch_bytes)
            out->reserved_parallel_scratch_bytes = scratch;
        out->total_nodes += child_count;
        out->interval_evaluations += child_count;
        for (size_t i = 0; i < child_count; i++) {
            double volume = cell_volume_box_volume(&children[i].bbox);
            size_t level = (size_t)children[i].level;
            if (level > out->deepest_level) out->deepest_level = level;
            for (int a = 0; a < 3; a++)
                if (children[i].axis_depth[a] > out->max_axis_depth[a])
                    out->max_axis_depth[a] = children[i].axis_depth[a];
            if (classes[i].relation == 1) { lower += volume; out->inside_nodes++; }
            else if (classes[i].relation == 0) out->outside_nodes++;
            else {
                children[i].active_axes = classes[i].active_axes;
                if (eligible_axes(&children[i], options)) {
                    children[i].interval_width = classes[i].interval_width;
                    heap_push(heap, &count, children[i]);
                    pending_gap += volume;
                } else {
                    terminal_gap += volume;
                    terminal_estimate += (long double)volume * cell_volume_sample_fraction(
                        sys, cell->root_node_id, &children[i].bbox, options->samples_per_axis);
                    out->unresolved_leaf_nodes++;
                    bool depth = false, size = false;
                    for (int a = 0; a < 3; a++) {
                        if (options->split_strategy == ALEA_CELL_VOLUME_SPLIT_ADAPTIVE &&
                            !(children[i].active_axes & (1u << a))) continue;
                        depth |= children[i].axis_depth[a] >= options->max_depth;
                        size |= options->min_size > 0.0 && axis_max(&children[i].bbox,a)-
                                axis_min(&children[i].bbox,a) <= options->min_size;
                    }
                    if (depth) out->max_depth_reached++;
                    if (size) out->min_size_reached++;
                    if (!depth && !size) out->precision_limit_reached++;
                }
            }
        }
        if (count == 0) pending_gap = 0.0L;
        if (count > out->frontier_task_count) out->frontier_task_count = count;
        long double gap = terminal_gap + pending_gap;
        /* Delay point samples until they can affect stopping or final output.
         * The midpoint estimate is only a trigger; actual samples confirm the
         * existing tolerance formula before convergence is declared. */
        double provisional = (double)(lower + terminal_estimate + pending_gap * 0.5L);
        double target = fmax(options->absolute_tolerance,
                             options->relative_tolerance * fabs(provisional));
        if (stop || count == 0 || gap <= target) {
            long double pending_estimate = sample_pending(sys, cell->root_node_id,
                heap, count, options->samples_per_axis);
            if (alea_interrupted()) goto interrupted;
            double estimate = (double)(lower + terminal_estimate + pending_estimate);
            target = fmax(options->absolute_tolerance, options->relative_tolerance * fabs(estimate));
            if (stop || count == 0 || gap <= target) {
                out->volume = estimate;
                out->lower_bound = (double)lower;
                out->unresolved_volume = (double)gap;
                out->upper_bound = (double)(lower + gap);
                out->unresolved_leaf_nodes += count;
                break;
            }
        }
        child_count = 0;
        /* Reserve both worst-case child storage and evaluation work before
         * removing a parent. An incomplete replacement never loses volume. */
        size_t reserve_children = options->split_strategy == ALEA_CELL_VOLUME_SPLIT_LONGEST_AXIS ? 2 : 8;
        while (count && child_count + reserve_children <= batch_capacity) {
            int eligible = eligible_axes(&heap[0], options);
            uint64_t probes = 0;
            if (options->split_strategy == ALEA_CELL_VOLUME_SPLIT_ADAPTIVE)
                for (int a = 0; a < 3; a++) if (eligible & (1 << a)) probes += 3;
            if (child_count + probes + reserve_children > options->max_evaluations - out->interval_evaluations) {
                out->evaluation_limit_reached = true; stop = true; break;
            }
            int mask = select_axes(sys, cell->root_node_id, &heap[0], eligible, options, out);
            size_t replacements = 1;
            for (int a = 0; a < 3; a++) if (mask & (1 << a)) replacements *= 2;
            if (count + child_count + replacements - 1 > capacity) {
                out->memory_limit_reached = true; stop = true; break;
            }
            cell_volume_task_t parent = heap_pop(heap, &count);
            for (int a = 0; a < 3; a++) if (mask & (1 << a)) out->axis_splits[a]++;
            pending_gap -= cell_volume_box_volume(&parent.bbox);
            for (int child = 0; child < 8; child++) {
                if (child & ~mask) continue;
                cell_volume_task_t next = parent;
                next.bbox = split_bbox(parent.bbox, mask, child);
                next.level++;
                for (int a = 0; a < 3; a++) if (mask & (1 << a)) next.axis_depth[a]++;
                children[child_count++] = next;
            }
        }
        if (child_count == 0 && !stop) {
            out->memory_limit_reached = true; stop = true;
        }
        /* If stop was set after building a partial batch, classify that batch
         * first. Its parents have already been removed from pending volume. */
    }
    free(heap); free(children); free(classes);
    if (alea_interrupted()) {
        alea_set_error_detail(ALEA_ERR_INTERRUPTED, "cell volume estimation interrupted");
        return -1;
    }
    if (!isfinite(out->volume) || !isfinite(out->lower_bound) || !isfinite(out->upper_bound)) {
        alea_set_error_detail(ALEA_ERR_OVERFLOW, "cell volume accumulation is not finite");
        return -1;
    }
    if (out->volume < out->lower_bound) out->volume = out->lower_bound;
    if (out->volume > out->upper_bound) out->volume = out->upper_bound;
    out->relative_uncertainty = out->unresolved_volume == 0.0 ? 0.0 :
        out->unresolved_volume / fmax(fabs(out->volume), DBL_MIN);
    out->converged = out->unresolved_volume <= fmax(options->absolute_tolerance,
                        options->relative_tolerance * fabs(out->volume));
    out->resource_limit_reached = out->evaluation_limit_reached || out->memory_limit_reached;
    return 0;
interrupted:
    free(heap); free(children); free(classes);
    alea_set_error_detail(ALEA_ERR_INTERRUPTED, "cell volume estimation interrupted");
    return -1;
}
