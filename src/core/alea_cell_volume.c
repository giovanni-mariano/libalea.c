// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/* Read-only, universe-local cell volume estimation. */

#include "alea.h"
#include "core/alea_system.h"
#include "core/alea_eval.h"
#include "primitives/bbox.h"

#include <float.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#define CELL_VOLUME_UNBOUNDED_EXTENT 9e5
#define CELL_VOLUME_DISCOVERY_DEPTH 6
#define CELL_VOLUME_DISCOVERY_EXPANSIONS 24
#define CELL_VOLUME_MAX_FRONTIER ((size_t)1024 * 1024)

typedef struct {
    alea_bbox_t bbox;
    int depth;
} cell_volume_task_t;

typedef struct {
    uint8_t relation; /* 0 outside, 1 inside, 2 mixed */
    double sample_fraction;
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
#ifdef _OPENMP
    if (omp_in_parallel()) return 1;
    if (workers == 0) workers = (size_t)omp_get_max_threads();
#else
    (void)requested;
    workers = 1;
#endif
    if (budget == 0 || tasks < 2) return 1;
    if (workers < 1) workers = 1;
    if (workers > tasks) workers = tasks;
    if (scratch_per_worker != 0) {
        uint64_t by_budget = budget / scratch_per_worker;
        if (by_budget == 0) return 1;
        if ((uint64_t)workers > by_budget) workers = (size_t)by_budget;
    }
    return workers ? workers : 1;
}

void alea_cell_volume_options_init(alea_cell_volume_options_t* options) {
    if (!options) return;
    *options = (alea_cell_volume_options_t){
        .has_bounds = false,
        .bounds = {0},
        .relative_tolerance = 1e-3,
        .absolute_tolerance = 0.0,
        .max_depth = 10,
        .min_size = 0.0,
        .samples_per_axis = 2,
        .requested_workers = 0,
        .max_parallel_scratch_bytes = 64u * 1024u * 1024u,
    };
}

static int cell_volume_options_valid(const alea_cell_volume_options_t* o) {
    return o && isfinite(o->relative_tolerance) &&
        o->relative_tolerance >= 0.0 && isfinite(o->absolute_tolerance) &&
        o->absolute_tolerance >= 0.0 && o->max_depth >= 0 &&
        isfinite(o->min_size) && o->min_size >= 0.0 &&
        o->samples_per_axis >= 1 && o->samples_per_axis <= 32 &&
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

    cell_volume_task_t* frontier = malloc(sizeof(*frontier));
    if (!frontier) {
        alea_set_error_detail(ALEA_ERR_OUT_OF_MEMORY,
                              "failed to allocate cell volume frontier");
        return -1;
    }
    frontier[0] = (cell_volume_task_t){out->bounds, 0};
    size_t frontier_count = 1;
    double lower = 0.0, terminal_gap = 0.0, terminal_estimate = 0.0;

    while (frontier_count != 0) {
        if (alea_interrupted()) {
            free(frontier);
            alea_set_error_detail(ALEA_ERR_INTERRUPTED,
                                  "cell volume estimation interrupted");
            return -1;
        }
        if (frontier_count > out->frontier_task_count)
            out->frontier_task_count = frontier_count;
        if (frontier_count > SIZE_MAX / sizeof(cell_volume_classification_t)) {
            free(frontier);
            alea_set_error_detail(ALEA_ERR_OVERFLOW,
                                  "cell volume classification allocation overflows");
            return -1;
        }
        cell_volume_classification_t* classes =
            calloc(frontier_count, sizeof(*classes));
        if (!classes) {
            free(frontier);
            alea_set_error_detail(ALEA_ERR_OUT_OF_MEMORY,
                                  "failed to allocate cell volume classifications");
            return -1;
        }
        size_t workers = cell_volume_select_workers(
            options->requested_workers, frontier_count,
            options->max_parallel_scratch_bytes,
            sizeof(cell_volume_classification_t));
        size_t actual_workers = 1;
#ifdef _OPENMP
        #pragma omp parallel num_threads(workers) if(workers > 1)
        {
            #pragma omp single
            actual_workers = (size_t)omp_get_num_threads();
            #pragma omp for schedule(static)
            for (size_t i = 0; i < frontier_count; i++) {
                classes[i].relation = cell_volume_relation(
                    sys, cell->root_node_id, &frontier[i].bbox);
                if (classes[i].relation == 2) {
                    classes[i].sample_fraction = cell_volume_sample_fraction(
                        sys, cell->root_node_id, &frontier[i].bbox,
                        options->samples_per_axis);
                }
            }
        }
#else
        for (size_t i = 0; i < frontier_count; i++) {
            classes[i].relation = cell_volume_relation(
                sys, cell->root_node_id, &frontier[i].bbox);
            if (classes[i].relation == 2) {
                classes[i].sample_fraction = cell_volume_sample_fraction(
                    sys, cell->root_node_id, &frontier[i].bbox,
                    options->samples_per_axis);
            }
        }
#endif
        if (actual_workers > out->actual_workers) out->actual_workers = actual_workers;
        if (actual_workers > 1) out->parallel_batch_count++;
        uint64_t reserved = (uint64_t)workers * sizeof(*classes);
        if (reserved > out->reserved_parallel_scratch_bytes)
            out->reserved_parallel_scratch_bytes = reserved;

        size_t expandable = 0;
        double pending_gap = 0.0, pending_estimate = 0.0;
        out->total_nodes += frontier_count;
        for (size_t i = 0; i < frontier_count; i++) {
            double box_volume = cell_volume_box_volume(&frontier[i].bbox);
            if (classes[i].relation == 1) {
                lower += box_volume;
                out->inside_nodes++;
            } else if (classes[i].relation == 0) {
                out->outside_nodes++;
            } else {
                bool depth_stop = frontier[i].depth >= options->max_depth;
                bool size_stop = options->min_size > 0.0 &&
                    cell_volume_min_extent(&frontier[i].bbox) <= options->min_size;
                if (depth_stop || size_stop) {
                    terminal_gap += box_volume;
                    terminal_estimate += classes[i].sample_fraction * box_volume;
                    out->unresolved_leaf_nodes++;
                    if (depth_stop) out->max_depth_reached++;
                } else {
                    expandable++;
                    pending_gap += box_volume;
                    pending_estimate += classes[i].sample_fraction * box_volume;
                }
            }
        }

        double estimate = lower + terminal_estimate + pending_estimate;
        double gap = terminal_gap + pending_gap;
        double target = fmax(options->absolute_tolerance,
                             options->relative_tolerance * fabs(estimate));
        bool converged = gap <= target;
        if (converged || expandable == 0) {
            if (converged) out->unresolved_leaf_nodes += expandable;
            out->converged = converged;
            out->volume = estimate;
            out->lower_bound = lower;
            out->unresolved_volume = gap;
            out->upper_bound = lower + gap;
            free(classes);
            free(frontier);
            frontier = NULL;
            frontier_count = 0;
            break;
        }

        if (expandable > CELL_VOLUME_MAX_FRONTIER / 8 ||
            expandable > SIZE_MAX / (8 * sizeof(cell_volume_task_t))) {
            out->resource_limit_reached = true;
            out->unresolved_leaf_nodes += expandable;
            out->converged = false;
            out->volume = estimate;
            out->lower_bound = lower;
            out->unresolved_volume = gap;
            out->upper_bound = lower + gap;
            free(classes);
            free(frontier);
            frontier = NULL;
            frontier_count = 0;
            break;
        }
        size_t next_count = expandable * 8;
        cell_volume_task_t* next = malloc(next_count * sizeof(*next));
        if (!next) {
            free(classes);
            free(frontier);
            alea_set_error_detail(ALEA_ERR_OUT_OF_MEMORY,
                                  "failed to grow cell volume frontier");
            return -1;
        }
        size_t write = 0;
        for (size_t i = 0; i < frontier_count; i++) {
            if (classes[i].relation != 2 ||
                frontier[i].depth >= options->max_depth ||
                (options->min_size > 0.0 &&
                 cell_volume_min_extent(&frontier[i].bbox) <= options->min_size))
                continue;
            for (int child = 0; child < 8; child++) {
                next[write++] = (cell_volume_task_t){
                    cell_volume_child_bbox(&frontier[i].bbox, child),
                    frontier[i].depth + 1};
            }
        }
        free(classes);
        free(frontier);
        frontier = next;
        frontier_count = write;
    }

    if (!isfinite(out->volume) || !isfinite(out->lower_bound) ||
        !isfinite(out->upper_bound)) {
        alea_set_error_detail(ALEA_ERR_OVERFLOW,
                              "cell volume accumulation is not finite");
        return -1;
    }
    if (out->volume < out->lower_bound) out->volume = out->lower_bound;
    if (out->volume > out->upper_bound) out->volume = out->upper_bound;
    out->relative_uncertainty = out->unresolved_volume == 0.0 ? 0.0 :
        out->unresolved_volume / fmax(fabs(out->volume), DBL_MIN);
    return 0;
}
