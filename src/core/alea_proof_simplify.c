// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/* Exact, proof-assisted simplification of one universe-local cell. */

#include "alea.h"
#include "core/alea_eval.h"
#include "core/alea_ops.h"
#include "core/alea_system.h"
#include "primitives/bbox.h"

#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#define PROOF_UNBOUNDED_EXTENT 9e5

typedef struct {
    alea_bbox_t bbox;
    int depth;
} proof_task_t;

typedef struct {
    uint8_t relation; /* 0 outside, 1 counterexample, 2 mixed */
    double witness[3];
} proof_classification_t;

typedef struct {
    size_t nodes;
    size_t surfaces;
    size_t depth;
} tree_size_t;

typedef struct {
    const alea_system_t* sys;
    const alea_node_id_t* terms;
    const uint8_t* keep;
    size_t count;
    size_t removed;
    alea_operation_t operation;
} proof_expression_t;

static bool proof_bbox_valid(const alea_bbox_t* b) {
    return b && isfinite(b->min_x) && isfinite(b->max_x) &&
        isfinite(b->min_y) && isfinite(b->max_y) &&
        isfinite(b->min_z) && isfinite(b->max_z) &&
        b->min_x < b->max_x && b->min_y < b->max_y &&
        b->min_z < b->max_z;
}

static bool proof_bbox_finite(const alea_bbox_t* b) {
    return proof_bbox_valid(b) &&
        b->max_x - b->min_x <= PROOF_UNBOUNDED_EXTENT &&
        b->max_y - b->min_y <= PROOF_UNBOUNDED_EXTENT &&
        b->max_z - b->min_z <= PROOF_UNBOUNDED_EXTENT;
}

static bool proof_bbox_empty(const alea_bbox_t* b) {
    return b && (b->min_x > b->max_x || b->min_y > b->max_y ||
                 b->min_z > b->max_z);
}

static bool proof_bbox_contains(const alea_bbox_t* outer,
                                const alea_bbox_t* inner) {
    return outer->min_x <= inner->min_x && outer->max_x >= inner->max_x &&
           outer->min_y <= inner->min_y && outer->max_y >= inner->max_y &&
           outer->min_z <= inner->min_z && outer->max_z >= inner->max_z;
}

static alea_bbox_t proof_child_bbox(const alea_bbox_t* b, int child) {
    double mx = b->min_x + (b->max_x - b->min_x) * 0.5;
    double my = b->min_y + (b->max_y - b->min_y) * 0.5;
    double mz = b->min_z + (b->max_z - b->min_z) * 0.5;
    return (alea_bbox_t){
        (child & 1) ? mx : b->min_x, (child & 1) ? b->max_x : mx,
        (child & 2) ? my : b->min_y, (child & 2) ? b->max_y : my,
        (child & 4) ? mz : b->min_z, (child & 4) ? b->max_z : mz,
    };
}

static size_t flatten_count(const alea_system_t* sys, alea_node_id_t id,
                            alea_operation_t operation) {
    const alea_node_t* node = &sys->nodes.data[id];
    if (ALEA_GET_OPERATION(node) != operation) return 1;
    return flatten_count(sys, node->operation.left, operation) +
           flatten_count(sys, node->operation.right, operation);
}

static void flatten_fill(const alea_system_t* sys, alea_node_id_t id,
                         alea_operation_t operation, alea_node_id_t* terms,
                         size_t* cursor) {
    const alea_node_t* node = &sys->nodes.data[id];
    if (ALEA_GET_OPERATION(node) != operation) {
        terms[(*cursor)++] = id;
        return;
    }
    flatten_fill(sys, node->operation.left, operation, terms, cursor);
    flatten_fill(sys, node->operation.right, operation, terms, cursor);
}

static tree_size_t tree_size(const alea_system_t* sys, alea_node_id_t id) {
    tree_size_t out = {1, 0, 1};
    const alea_node_t* node = &sys->nodes.data[id];
    alea_operation_t op = ALEA_GET_OPERATION(node);
    if (op == ALEA_OP_PRIMITIVE) {
        out.surfaces = 1;
    } else if (op == ALEA_OP_COMPLEMENT) {
        tree_size_t child = tree_size(sys, node->operation.left);
        out.nodes += child.nodes;
        out.surfaces = child.surfaces;
        out.depth += child.depth;
    } else {
        tree_size_t left = tree_size(sys, node->operation.left);
        tree_size_t right = tree_size(sys, node->operation.right);
        out.nodes += left.nodes + right.nodes;
        out.surfaces = left.surfaces + right.surfaces;
        out.depth += left.depth > right.depth ? left.depth : right.depth;
    }
    return out;
}

static tree_size_t selected_size_range(const alea_system_t* sys,
                                       const alea_node_id_t* selected,
                                       size_t begin, size_t end) {
    if (end - begin == 1) return tree_size(sys, selected[begin]);
    size_t middle = begin + (end - begin) / 2;
    tree_size_t left = selected_size_range(sys, selected, begin, middle);
    tree_size_t right = selected_size_range(sys, selected, middle, end);
    return (tree_size_t){left.nodes + right.nodes + 1,
                         left.surfaces + right.surfaces,
                         (left.depth > right.depth ? left.depth : right.depth) + 1};
}

static void collect_unique_primitives(const alea_system_t* sys,
                                      alea_node_id_t id,
                                      alea_primitive_id_t* ids,
                                      size_t* count) {
    const alea_node_t* node = &sys->nodes.data[id];
    alea_operation_t op = ALEA_GET_OPERATION(node);
    if (op == ALEA_OP_PRIMITIVE) {
        alea_primitive_id_t primitive = node->primitive.primitive_id;
        for (size_t i = 0; i < *count; i++) if (ids[i] == primitive) return;
        ids[(*count)++] = primitive;
        return;
    }
    collect_unique_primitives(sys, node->operation.left, ids, count);
    if (op != ALEA_OP_COMPLEMENT)
        collect_unique_primitives(sys, node->operation.right, ids, count);
}

static size_t unique_surface_count(const alea_system_t* sys,
                                   const alea_node_id_t* roots,
                                   size_t root_count, size_t leaf_upper_bound) {
    if (leaf_upper_bound == 0) return 0;
    alea_primitive_id_t* ids = malloc(leaf_upper_bound * sizeof(*ids));
    if (!ids) return leaf_upper_bound;
    size_t count = 0;
    for (size_t i = 0; i < root_count; i++)
        collect_unique_primitives(sys, roots[i], ids, &count);
    free(ids);
    return count;
}

static alea_interval_t proof_interval_node(const alea_system_t* sys,
                                           alea_node_id_t id,
                                           const alea_bbox_t* box) {
    alea_bbox_t expanded = {
        nextafter(box->min_x, -INFINITY), nextafter(box->max_x, INFINITY),
        nextafter(box->min_y, -INFINITY), nextafter(box->max_y, INFINITY),
        nextafter(box->min_z, -INFINITY), nextafter(box->max_z, INFINITY),
    };
    alea_interval_t iv = alea_evaluate_interval(sys, id, &expanded);
    if (!isfinite(iv.min) || !isfinite(iv.max) || iv.min > iv.max)
        return (alea_interval_t){-INFINITY, INFINITY};
    iv.min = nextafter(iv.min, -INFINITY);
    iv.max = nextafter(iv.max, INFINITY);
    return iv;
}

static alea_interval_t retained_interval(const proof_expression_t* expr,
                                         const alea_bbox_t* box) {
    alea_interval_t result = expr->operation == ALEA_OP_INTERSECTION
        ? (alea_interval_t){-INFINITY, -INFINITY}
        : (alea_interval_t){INFINITY, INFINITY};
    for (size_t i = 0; i < expr->count; i++) {
        if (!expr->keep[i] || i == expr->removed) continue;
        alea_interval_t iv = proof_interval_node(expr->sys, expr->terms[i], box);
        if (expr->operation == ALEA_OP_INTERSECTION) {
            result.min = fmax(result.min, iv.min);
            result.max = fmax(result.max, iv.max);
        } else {
            result.min = fmin(result.min, iv.min);
            result.max = fmin(result.max, iv.max);
        }
    }
    return result;
}

/* One-term removal gives a subset relation in one direction, so the full
 * symmetric difference reduces to a single directed set difference. */
static alea_interval_t difference_interval(const proof_expression_t* expr,
                                           const alea_bbox_t* box) {
    alea_interval_t retained = retained_interval(expr, box);
    alea_interval_t removed = proof_interval_node(
        expr->sys, expr->terms[expr->removed], box);
    alea_interval_t left = expr->operation == ALEA_OP_INTERSECTION
        ? retained : removed;
    alea_interval_t right = expr->operation == ALEA_OP_INTERSECTION
        ? removed : retained;
    return (alea_interval_t){fmax(left.min, -right.max),
                             fmax(left.max, -right.min)};
}

static bool retained_contains(const proof_expression_t* expr,
                              double x, double y, double z) {
    bool value = expr->operation == ALEA_OP_INTERSECTION;
    for (size_t i = 0; i < expr->count; i++) {
        if (!expr->keep[i] || i == expr->removed) continue;
        bool inside = alea_contains_point(expr->sys, expr->terms[i], x, y, z);
        if (expr->operation == ALEA_OP_INTERSECTION) value = value && inside;
        else value = value || inside;
    }
    return value;
}

static bool difference_contains(const proof_expression_t* expr,
                                double x, double y, double z) {
    bool retained = retained_contains(expr, x, y, z);
    bool removed = alea_contains_point(expr->sys, expr->terms[expr->removed],
                                       x, y, z);
    return expr->operation == ALEA_OP_INTERSECTION
        ? retained && !removed : removed && !retained;
}

static bool find_box_witness(const proof_expression_t* expr,
                             const alea_bbox_t* b, double out[3]) {
    const double xs[3] = {b->min_x, b->min_x + (b->max_x-b->min_x)*0.5,
                          b->max_x};
    const double ys[3] = {b->min_y, b->min_y + (b->max_y-b->min_y)*0.5,
                          b->max_y};
    const double zs[3] = {b->min_z, b->min_z + (b->max_z-b->min_z)*0.5,
                          b->max_z};
    static const uint8_t order[9][3] = {
        {1,1,1}, {0,0,0}, {2,0,0}, {0,2,0}, {2,2,0},
        {0,0,2}, {2,0,2}, {0,2,2}, {2,2,2}
    };
    for (size_t i = 0; i < 9; i++) {
        double x = xs[order[i][0]], y = ys[order[i][1]], z = zs[order[i][2]];
        if (difference_contains(expr, x, y, z)) {
            out[0] = x; out[1] = y; out[2] = z;
            return true;
        }
    }
    return false;
}

static proof_classification_t classify_difference_box(
        const proof_expression_t* expr, const alea_bbox_t* box) {
    proof_classification_t out = {0};
    alea_interval_t iv = difference_interval(expr, box);
    if (iv.min > 0.0) return out;
    if (find_box_witness(expr, box, out.witness)) {
        out.relation = 1;
        return out;
    }
    out.relation = 2;
    return out;
}

static size_t proof_select_workers(size_t requested, size_t tasks,
                                   uint64_t budget) {
    size_t workers = requested;
#ifdef _OPENMP
    if (omp_in_parallel()) return 1;
    if (workers == 0) workers = (size_t)omp_get_max_threads();
#else
    (void)requested;
    workers = 1;
#endif
    if (workers == 0) workers = 1;
    if (workers > tasks) workers = tasks;
    size_t budget_workers = budget / sizeof(proof_classification_t);
    if (workers > budget_workers) workers = budget_workers;
    return workers ? workers : 1;
}

static int automatic_domain(const proof_expression_t* expr, alea_bbox_t* out,
                            alea_proof_bounds_source_t* source) {
    if (expr->operation == ALEA_OP_UNION) {
        *out = alea_get_bbox(expr->sys, expr->terms[expr->removed]);
        if (proof_bbox_finite(out)) {
            *source = ALEA_PROOF_BOUNDS_STORED;
            return 0;
        }
        if (alea_tighten_bbox_plane_constraints(expr->sys,
                expr->terms[expr->removed], 1e-8, out) == 0 &&
            proof_bbox_finite(out)) {
            *source = ALEA_PROOF_BOUNDS_PLANE_CONSTRAINTS;
            return 0;
        }
        return -1;
    }

    bool have = false;
    size_t retained_count = 0;
    alea_node_id_t sole = ALEA_NODE_ID_INVALID;
    for (size_t i = 0; i < expr->count; i++) {
        if (!expr->keep[i] || i == expr->removed) continue;
        retained_count++;
        sole = expr->terms[i];
        alea_bbox_t term = alea_get_bbox(expr->sys, expr->terms[i]);
        if (!have) { *out = term; have = true; }
        else *out = alea_bbox_intersection(out, &term);
    }
    if (have && proof_bbox_finite(out)) {
        *source = ALEA_PROOF_BOUNDS_STORED;
        return 0;
    }
    if (retained_count == 1 &&
        alea_tighten_bbox_plane_constraints(expr->sys, sole, 1e-8, out) == 0 &&
        proof_bbox_finite(out)) {
        *source = ALEA_PROOF_BOUNDS_PLANE_CONSTRAINTS;
        return 0;
    }
    return -1;
}

static alea_proof_status_t prove_difference_empty(
    const proof_expression_t* expr, const alea_bbox_t* bounds,
    const alea_cell_simplify_proof_options_t* options,
    alea_cell_simplify_proof_result_t* result) {
    if (result->proof_nodes >= options->max_nodes) {
        result->last_limit = ALEA_PROOF_LIMIT_NODES;
        result->mixed_leaf_nodes++;
        return ALEA_PROOF_INCONCLUSIVE;
    }
    proof_classification_t root_class = classify_difference_box(expr, bounds);
    result->proof_nodes++;
    if (root_class.relation == 0) return ALEA_PROOF_PROVEN;
    if (root_class.relation == 1) {
        result->has_witness = true;
        memcpy(result->witness, root_class.witness, sizeof(result->witness));
        return ALEA_PROOF_DISPROVEN;
    }
    if (options->max_depth == 0) {
        result->last_limit = ALEA_PROOF_LIMIT_DEPTH;
        result->mixed_leaf_nodes++;
        return ALEA_PROOF_INCONCLUSIVE;
    }

    size_t capacity = 64;
    if (capacity > options->max_nodes) capacity = options->max_nodes;
    proof_task_t* stack = capacity ? malloc(capacity * sizeof(*stack)) : NULL;
    if (!stack) {
        result->last_limit = ALEA_PROOF_LIMIT_MEMORY;
        return ALEA_PROOF_INCONCLUSIVE;
    }
    size_t count = 1;
    stack[0] = (proof_task_t){*bounds, 0};
    bool incomplete = false;
    while (count) {
        if (alea_interrupted()) {
            free(stack);
            return (alea_proof_status_t)-1;
        }
        proof_task_t task = stack[--count];
        if (task.depth >= options->max_depth) {
            result->mixed_leaf_nodes++;
            result->last_limit = ALEA_PROOF_LIMIT_DEPTH;
            incomplete = true;
            continue;
        }
        if (result->proof_nodes + 8 > options->max_nodes) {
            result->mixed_leaf_nodes++;
            result->last_limit = ALEA_PROOF_LIMIT_NODES;
            incomplete = true;
            continue;
        }
        if (count + 8 > capacity) {
            size_t next = capacity < options->max_nodes / 2
                ? capacity * 2 : options->max_nodes;
            if (next < count + 8) next = count + 8;
            proof_task_t* grown = realloc(stack, next * sizeof(*grown));
            if (!grown) {
                result->last_limit = ALEA_PROOF_LIMIT_MEMORY;
                incomplete = true;
                continue;
            }
            stack = grown;
            capacity = next;
        }
        proof_task_t children[8];
        proof_classification_t classes[8];
        for (int child = 0; child < 8; child++)
            children[child] = (proof_task_t){
                proof_child_bbox(&task.bbox, child), task.depth + 1};
        size_t workers = proof_select_workers(options->requested_workers, 8,
            options->max_parallel_scratch_bytes);
        if (workers > result->actual_workers) result->actual_workers = workers;
        uint64_t reserved = workers > 1
            ? (uint64_t)workers * sizeof(*classes) : 0;
        if (reserved > result->reserved_parallel_scratch_bytes)
            result->reserved_parallel_scratch_bytes = reserved;
#ifdef _OPENMP
        if (workers > 1) {
            result->parallel_batch_count++;
            #pragma omp parallel for num_threads((int)workers) schedule(static)
            for (int child = 0; child < 8; child++)
                classes[child] = classify_difference_box(
                    expr, &children[child].bbox);
        } else
#endif
        {
            for (int child = 0; child < 8; child++)
                classes[child] = classify_difference_box(
                    expr, &children[child].bbox);
        }
        result->proof_nodes += 8;
        /* Stable serial reduction chooses the lowest child ordinal witness. */
        for (int child = 0; child < 8; child++) {
            if (classes[child].relation == 1) {
                result->has_witness = true;
                memcpy(result->witness, classes[child].witness,
                       sizeof(result->witness));
                free(stack);
                return ALEA_PROOF_DISPROVEN;
            }
        }
        /* Reverse push makes child zero the next deterministic mixed task. */
        for (int child = 7; child >= 0; child--)
            if (classes[child].relation == 2) stack[count++] = children[child];
        if (count > result->frontier_task_count)
            result->frontier_task_count = count;
    }
    free(stack);
    return incomplete ? ALEA_PROOF_INCONCLUSIVE : ALEA_PROOF_PROVEN;
}

alea_proof_status_t alea_prove_union_branch_redundant(
    const alea_system_t* sys, alea_node_id_t branch,
    const alea_node_id_t* remaining, size_t remaining_count,
    const alea_bbox_t* complete_branch_bounds, int max_depth,
    size_t max_nodes) {
    if (!sys || branch == ALEA_NODE_ID_INVALID || !remaining ||
        remaining_count == 0 || !proof_bbox_valid(complete_branch_bounds) ||
        max_depth < 0 || max_nodes == 0) return ALEA_PROOF_INCONCLUSIVE;
    size_t count = remaining_count + 1;
    if (count < remaining_count) return ALEA_PROOF_INCONCLUSIVE;
    alea_node_id_t* terms = malloc(count * sizeof(*terms));
    uint8_t* keep = malloc(count);
    if (!terms || !keep) {
        free(terms); free(keep);
        return ALEA_PROOF_INCONCLUSIVE;
    }
    terms[0] = branch;
    memcpy(terms + 1, remaining, remaining_count * sizeof(*remaining));
    memset(keep, 1, count);
    proof_expression_t expr = {
        sys, terms, keep, count, 0, ALEA_OP_UNION
    };
    alea_cell_simplify_proof_options_t options;
    alea_cell_simplify_proof_options_init(&options);
    options.max_depth = max_depth;
    options.max_nodes = max_nodes;
    options.requested_workers = 1;
    options.max_parallel_scratch_bytes = 0;
    alea_cell_simplify_proof_result_t receipt;
    memset(&receipt, 0, sizeof(receipt));
    receipt.actual_workers = 1;
    alea_proof_status_t status = prove_difference_empty(
        &expr, complete_branch_bounds, &options, &receipt);
    free(terms); free(keep);
    return (int)status < 0 ? ALEA_PROOF_INCONCLUSIVE : status;
}

static void analyze_patterns(const alea_system_t* sys,
                             const alea_node_id_t* terms, size_t count,
                             alea_operation_t operation, const alea_bbox_t* box,
                             size_t max_patterns, size_t* scores,
                             size_t* patterns_collected) {
    if (!proof_bbox_finite(box) || max_patterns == 0) return;
    size_t side = 1;
    while ((side + 1) * (side + 1) * (side + 1) <= max_patterns && side < 8)
        side++;
    uint64_t* hashes = malloc(max_patterns * sizeof(*hashes));
    if (!hashes) return;
    for (size_t k = 0; k < side; k++) for (size_t j = 0; j < side; j++)
    for (size_t i = 0; i < side; i++) {
        double x = box->min_x + ((double)i + 0.5) *
            (box->max_x-box->min_x) / (double)side;
        double y = box->min_y + ((double)j + 0.5) *
            (box->max_y-box->min_y) / (double)side;
        double z = box->min_z + ((double)k + 0.5) *
            (box->max_z-box->min_z) / (double)side;
        uint64_t hash = UINT64_C(1469598103934665603);
        size_t inside_count = 0, decisive = SIZE_MAX;
        for (size_t term = 0; term < count; term++) {
            bool inside = alea_contains_point(sys, terms[term], x, y, z);
            hash ^= inside ? UINT64_C(1) : UINT64_C(0);
            hash *= UINT64_C(1099511628211);
            if (inside) { inside_count++; decisive = term; }
            else if (operation == ALEA_OP_INTERSECTION) decisive = term;
        }
        bool unique = true;
        for (size_t h = 0; h < *patterns_collected; h++)
            if (hashes[h] == hash) { unique = false; break; }
        if (unique && *patterns_collected < max_patterns)
            hashes[(*patterns_collected)++] = hash;
        if (operation == ALEA_OP_UNION && inside_count == 1)
            scores[decisive]++;
        if (operation == ALEA_OP_INTERSECTION && inside_count + 1 == count)
            scores[decisive]++;
    }
    free(hashes);
}

void alea_cell_simplify_proof_options_init(
    alea_cell_simplify_proof_options_t* options) {
    if (!options) return;
    memset(options, 0, sizeof(*options));
    options->max_depth = 12;
    options->max_nodes = 250000;
    options->max_patterns = 64;
    options->max_candidates = 256;
    options->max_parallel_scratch_bytes = 64u * 1024u * 1024u;
}

typedef struct {
    alea_operation_t operation;
    alea_node_id_t* selected;
    size_t selected_count;
} cell_proof_analysis_t;

enum {
    CELL_ANALYSIS_OK = 0,
    CELL_ANALYSIS_MEMORY = 1,
    CELL_ANALYSIS_INTERRUPTED = 2
};

static void cell_proof_analysis_free(cell_proof_analysis_t* analysis) {
    if (!analysis) return;
    free(analysis->selected);
    memset(analysis, 0, sizeof(*analysis));
}

static int analyze_cell_proven(
    const alea_system_t* sys, size_t cell_index,
    const alea_cell_simplify_proof_options_t* supplied,
    alea_cell_simplify_proof_result_t* result,
    cell_proof_analysis_t* analysis) {
    const alea_cell_entry_t* cell = &sys->cells.data[cell_index];
    memset(analysis, 0, sizeof(*analysis));
    memset(result, 0, sizeof(*result));
    result->complete = true;
    result->last_limit = ALEA_PROOF_LIMIT_NONE;
    result->root_node_id = cell->root_node_id;
    result->requested_workers = supplied->requested_workers;
    result->actual_workers = 1;
    tree_size_t before = tree_size(sys, cell->root_node_id);
    result->nodes_before = before.nodes;
    result->nodes_after = before.nodes;
    result->surfaces_before = unique_surface_count(
        sys, &cell->root_node_id, 1, before.surfaces);
    result->surfaces_after = result->surfaces_before;
    result->depth_before = before.depth;
    result->depth_after = before.depth;

    alea_bbox_t root_support = alea_get_bbox(sys, cell->root_node_id);
    if (proof_bbox_empty(&root_support)) {
        /* Boolean bbox intersection is conservative. A strictly empty
         * enclosure therefore certifies that the represented set is empty;
         * equality is intentionally not treated as empty because boundaries
         * belong to libalea half-spaces. */
        result->proven_empty = true;
        return 0;
    }
    if (supplied->has_bounds) {
        result->bounds = supplied->bounds;
        result->bounds_source = ALEA_PROOF_BOUNDS_EXPLICIT;
        result->bounds_verified = proof_bbox_finite(&root_support) &&
            proof_bbox_contains(&supplied->bounds, &root_support);
    } else if (proof_bbox_finite(&root_support)) {
        result->bounds = root_support;
        result->bounds_source = ALEA_PROOF_BOUNDS_STORED;
        result->bounds_verified = true;
    }

    alea_operation_t operation = ALEA_GET_OPERATION(
        &sys->nodes.data[cell->root_node_id]);
    if (operation != ALEA_OP_INTERSECTION && operation != ALEA_OP_UNION)
        return 0;

    size_t count = flatten_count(sys, cell->root_node_id, operation);
    if (count < 2) return 0;
    alea_node_id_t* terms = malloc(count * sizeof(*terms));
    alea_node_id_t* selected = malloc(count * sizeof(*selected));
    uint8_t* keep = malloc(count);
    size_t* order = malloc(count * sizeof(*order));
    size_t* scores = calloc(count, sizeof(*scores));
    if (!terms || !selected || !keep || !order || !scores) {
        free(terms); free(selected); free(keep); free(order); free(scores);
        return CELL_ANALYSIS_MEMORY;
    }
    size_t cursor = 0;
    flatten_fill(sys, cell->root_node_id, operation, terms, &cursor);
    memset(keep, 1, count);
    for (size_t i = 0; i < count; i++) order[i] = i;

    alea_bbox_t discovery = supplied->has_bounds ? supplied->bounds : root_support;
    analyze_patterns(sys, terms, count, operation, &discovery,
                     supplied->max_patterns, scores,
                     &result->patterns_collected);
    /* Low observed necessity is the promising removal order; ties retain the
     * source-tree order, making proposals and witnesses reproducible. */
    for (size_t i = 1; i < count; i++) {
        size_t item = order[i], j = i;
        while (j && scores[order[j-1]] > scores[item]) {
            order[j] = order[j-1]; j--;
        }
        order[j] = item;
    }

    size_t retained_count = count;
    for (size_t candidate = 0;
         candidate < count && result->candidates_proposed < supplied->max_candidates;
         candidate++) {
        if (retained_count <= 1) break;
        size_t removed = order[candidate];
        if (!keep[removed]) continue;
        proof_expression_t expr = {sys, terms, keep, count, removed, operation};
        alea_bbox_t domain, automatic;
        alea_proof_bounds_source_t source = ALEA_PROOF_BOUNDS_STORED;
        bool have_automatic = automatic_domain(&expr, &automatic, &source) == 0;
        if (supplied->has_bounds) {
            domain = supplied->bounds;
            result->bounds_source = ALEA_PROOF_BOUNDS_EXPLICIT;
            result->bounds_verified = result->bounds_verified && have_automatic &&
                proof_bbox_contains(&domain, &automatic);
        } else if (have_automatic) {
            domain = automatic;
            result->bounds_source = source;
            result->bounds_verified = true;
        } else {
            result->candidates_proposed++;
            result->candidates_inconclusive++;
            result->complete = false;
            result->last_limit = ALEA_PROOF_LIMIT_DOMAIN;
            continue;
        }
        result->bounds = domain;
        result->candidates_proposed++;
        alea_proof_status_t status = prove_difference_empty(
            &expr, &domain, supplied, result);
        if ((int)status < 0) {
            free(terms); free(selected); free(keep); free(order); free(scores);
            return CELL_ANALYSIS_INTERRUPTED;
        }
        if (status == ALEA_PROOF_PROVEN) {
            keep[removed] = 0;
            retained_count--;
            result->candidates_proven++;
            result->changed = true;
        } else if (status == ALEA_PROOF_DISPROVEN) {
            result->candidates_disproven++;
        } else {
            result->candidates_inconclusive++;
            result->complete = false;
        }
    }

    cursor = 0;
    for (size_t i = 0; i < count; i++) if (keep[i]) selected[cursor++] = terms[i];
    tree_size_t after = selected_size_range(sys, selected, 0, cursor);
    result->nodes_after = after.nodes;
    result->surfaces_after = unique_surface_count(
        sys, selected, cursor, after.surfaces);
    result->depth_after = after.depth;
    if (result->changed) {
        analysis->operation = operation;
        analysis->selected = selected;
        analysis->selected_count = cursor;
        selected = NULL;
    }
    free(terms); free(selected); free(keep); free(order); free(scores);
    return CELL_ANALYSIS_OK;
}

static alea_node_id_t materialize_cell_analysis(
        alea_system_t* sys, const cell_proof_analysis_t* analysis) {
    if (!analysis || !analysis->selected || analysis->selected_count == 0)
        return ALEA_NODE_ID_INVALID;
    if (analysis->selected_count == 1) return analysis->selected[0];
    return analysis->operation == ALEA_OP_INTERSECTION
        ? alea_create_intersection_many(sys, analysis->selected,
                                        analysis->selected_count)
        : alea_create_union_many(sys, analysis->selected,
                                 analysis->selected_count);
}

static int validate_cell_options(
        const alea_cell_simplify_proof_options_t* options) {
    return options && options->max_depth >= 0 && options->max_nodes != 0 &&
        options->max_candidates != 0 &&
        (!options->has_bounds || proof_bbox_valid(&options->bounds));
}

int alea_cell_simplify_proven(
    alea_system_t* sys, size_t cell_index,
    const alea_cell_simplify_proof_options_t* supplied,
    alea_cell_simplify_proof_result_t* result) {
    alea_cell_simplify_proof_options_t defaults;
    if (!supplied) {
        alea_cell_simplify_proof_options_init(&defaults);
        supplied = &defaults;
    }
    if (!sys || !result || cell_index >= alea_vec_count(&sys->cells) ||
        !validate_cell_options(supplied)) {
        alea_set_error_detail(ALEA_ERR_INVALID_ARG,
                              "invalid cell proof-simplification arguments");
        return -1;
    }
    alea_cell_entry_t* cell = &sys->cells.data[cell_index];
    if (cell->root_node_id == ALEA_NODE_ID_INVALID ||
        cell->root_node_id >= alea_vec_count(&sys->nodes)) {
        alea_set_error_detail(ALEA_ERR_INVALID_STATE,
                              "cell %zu has no valid CSG root", cell_index);
        return -1;
    }
    cell_proof_analysis_t analysis;
    int rc = analyze_cell_proven(sys, cell_index, supplied, result, &analysis);
    if (rc != CELL_ANALYSIS_OK) {
        cell_proof_analysis_free(&analysis);
        if (rc == CELL_ANALYSIS_INTERRUPTED)
            alea_set_error_detail(ALEA_ERR_INTERRUPTED,
                                  "cell proof simplification interrupted");
        else
            alea_set_error_detail(ALEA_ERR_OUT_OF_MEMORY,
                                  "failed to allocate proof simplifier workspace");
        return -1;
    }
    if (result->changed && supplied->apply) {
        alea_node_id_t root = materialize_cell_analysis(sys, &analysis);
        if (root == ALEA_NODE_ID_INVALID) {
            cell_proof_analysis_free(&analysis);
            return -1;
        }
        cell->root_node_id = root;
        cell->original_root_node_id = ALEA_NODE_ID_INVALID;
        result->root_node_id = root;
        result->applied = true;
        alea_system_invalidate_query_caches(sys, ALEA_CACHE_ALL);
    }
    cell_proof_analysis_free(&analysis);
    return 0;
}

void alea_cells_simplify_proof_options_init(
        alea_cells_simplify_proof_options_t* options) {
    if (!options) return;
    memset(options, 0, sizeof(*options));
    options->max_depth = 12;
    options->max_nodes_per_cell = 250000;
    options->max_patterns_per_cell = 64;
    options->max_candidates_per_cell = 256;
    options->max_parallel_scratch_bytes = 64u * 1024u * 1024u;
}

static uint64_t cell_analysis_scratch_estimate(
        const alea_system_t* sys, size_t cell_index) {
    alea_node_id_t root = sys->cells.data[cell_index].root_node_id;
    alea_operation_t operation = ALEA_GET_OPERATION(&sys->nodes.data[root]);
    size_t terms = (operation == ALEA_OP_INTERSECTION || operation == ALEA_OP_UNION)
        ? flatten_count(sys, root, operation) : 1;
    const uint64_t per_term = 2u * sizeof(alea_node_id_t) +
        sizeof(uint8_t) + 2u * sizeof(size_t);
    if (terms > (UINT64_MAX - 8u * sizeof(proof_classification_t)) / per_term)
        return UINT64_MAX;
    return (uint64_t)terms * per_term +
        8u * sizeof(proof_classification_t);
}

static size_t batch_select_workers(size_t requested, size_t tasks,
                                   uint64_t budget, uint64_t per_worker) {
    size_t workers = requested;
#ifdef _OPENMP
    if (omp_in_parallel()) return 1;
    if (workers == 0) workers = (size_t)omp_get_max_threads();
#else
    (void)requested;
    workers = 1;
#endif
    if (workers == 0) workers = 1;
    if (workers > tasks) workers = tasks;
    if (budget == 0 || per_worker == UINT64_MAX) return 1;
    uint64_t by_budget = per_worker ? budget / per_worker : workers;
    if ((uint64_t)workers > by_budget) workers = (size_t)by_budget;
    return workers ? workers : 1;
}

static void batch_free_analyses(cell_proof_analysis_t* analyses, size_t count) {
    if (!analyses) return;
    for (size_t i = 0; i < count; i++) cell_proof_analysis_free(&analyses[i]);
    free(analyses);
}

int alea_cells_simplify_proven(
    alea_system_t* sys, const alea_cell_simplify_request_t* requests,
    size_t request_count, const alea_cells_simplify_proof_options_t* supplied,
    alea_cell_simplify_proof_result_t* results,
    alea_cells_simplify_proof_summary_t* summary) {
    alea_cells_simplify_proof_options_t defaults;
    if (!supplied) {
        alea_cells_simplify_proof_options_init(&defaults);
        supplied = &defaults;
    }
    if (!sys || (!requests && request_count) || (!results && request_count) ||
        !summary || supplied->max_depth < 0 ||
        supplied->max_nodes_per_cell == 0 ||
        supplied->max_candidates_per_cell == 0) {
        alea_set_error_detail(ALEA_ERR_INVALID_ARG,
                              "invalid batch proof-simplification arguments");
        return -1;
    }
    memset(summary, 0, sizeof(*summary));
    summary->selected_cells = request_count;
    summary->requested_workers = supplied->requested_workers;
    summary->actual_workers = 1;
    if (request_count == 0) return 0;

    size_t cell_count = alea_vec_count(&sys->cells);
    uint64_t max_workspace = 1;
    for (size_t i = 0; i < request_count; i++) {
        const alea_cell_simplify_request_t* request = &requests[i];
        if (request->cell_index >= cell_count ||
            (request->has_bounds && !proof_bbox_valid(&request->bounds))) {
            alea_set_error_detail(ALEA_ERR_INVALID_ARG,
                                  "invalid proof request at ordinal %zu", i);
            return -1;
        }
        for (size_t j = 0; j < i; j++) {
            if (requests[j].cell_index == request->cell_index) {
                alea_set_error_detail(ALEA_ERR_INVALID_ARG,
                                      "duplicate proof request for cell %zu",
                                      request->cell_index);
                return -1;
            }
        }
        alea_node_id_t root = sys->cells.data[request->cell_index].root_node_id;
        if (root == ALEA_NODE_ID_INVALID || root >= alea_vec_count(&sys->nodes)) {
            alea_set_error_detail(ALEA_ERR_INVALID_STATE,
                                  "cell %zu has no valid CSG root",
                                  request->cell_index);
            return -1;
        }
        uint64_t workspace = cell_analysis_scratch_estimate(
            sys, request->cell_index);
        if (workspace > max_workspace) max_workspace = workspace;
    }

    cell_proof_analysis_t* analyses = calloc(request_count, sizeof(*analyses));
    int* errors = calloc(request_count, sizeof(*errors));
    if (!analyses || !errors) {
        free(analyses); free(errors);
        alea_set_error_detail(ALEA_ERR_OUT_OF_MEMORY,
                              "failed to allocate batch proof workspace");
        return -1;
    }
    uint64_t generation = alea_system_geometry_generation(sys);
    size_t outer_workers = request_count == 1 ? 1 : batch_select_workers(
        supplied->requested_workers, request_count,
        supplied->max_parallel_scratch_bytes, max_workspace);
    size_t observed_workers = 1;

#ifdef _OPENMP
    #pragma omp parallel if(outer_workers > 1) num_threads((int)outer_workers)
    {
        #pragma omp single
        observed_workers = (size_t)omp_get_num_threads();
        #pragma omp for schedule(static)
#endif
        for (size_t i = 0; i < request_count; i++) {
            alea_cell_simplify_proof_options_t local;
            alea_cell_simplify_proof_options_init(&local);
            local.has_bounds = requests[i].has_bounds;
            local.bounds = requests[i].bounds;
            local.apply = false;
            local.max_depth = supplied->max_depth;
            local.max_nodes = supplied->max_nodes_per_cell;
            local.max_patterns = supplied->max_patterns_per_cell;
            local.max_candidates = supplied->max_candidates_per_cell;
            local.requested_workers = request_count == 1
                ? supplied->requested_workers : 1;
            local.max_parallel_scratch_bytes = request_count == 1
                ? supplied->max_parallel_scratch_bytes : 0;
            errors[i] = analyze_cell_proven(sys, requests[i].cell_index, &local,
                                            &results[i], &analyses[i]);
            results[i].requested_workers = supplied->requested_workers;
        }
#ifdef _OPENMP
    }
#endif
    summary->actual_workers = observed_workers;
    summary->reserved_parallel_scratch_bytes = observed_workers > 1
        ? (uint64_t)observed_workers * max_workspace : 0;
    if (observed_workers > 1) summary->parallel_batch_count = 1;

    for (size_t i = 0; i < request_count; i++) {
        if (errors[i] != CELL_ANALYSIS_OK) {
            int error = errors[i];
            batch_free_analyses(analyses, request_count);
            free(errors);
            if (error == CELL_ANALYSIS_INTERRUPTED)
                alea_set_error_detail(ALEA_ERR_INTERRUPTED,
                    "batch proof simplification interrupted at request %zu", i);
            else
                alea_set_error_detail(ALEA_ERR_OUT_OF_MEMORY,
                    "batch proof workspace allocation failed at request %zu", i);
            return -1;
        }
        if (results[i].changed) summary->changed_cells++;
        if (results[i].proven_empty) summary->proven_empty_cells++;
        if (results[i].complete) summary->complete_cells++;
        else summary->inconclusive_cells++;
    }
    free(errors);

    if (request_count == 1) {
        summary->actual_workers = results[0].actual_workers;
        summary->parallel_batch_count = results[0].parallel_batch_count;
        summary->reserved_parallel_scratch_bytes =
            results[0].reserved_parallel_scratch_bytes;
    }
    if (!supplied->apply || summary->changed_cells == 0) {
        batch_free_analyses(analyses, request_count);
        return 0;
    }
    if (alea_system_geometry_generation(sys) != generation) {
        batch_free_analyses(analyses, request_count);
        alea_set_error_detail(ALEA_ERR_INVALID_STATE,
                              "geometry changed during batch proof analysis");
        return -1;
    }

    alea_node_id_t* roots = malloc(request_count * sizeof(*roots));
    if (!roots) {
        batch_free_analyses(analyses, request_count);
        alea_set_error_detail(ALEA_ERR_OUT_OF_MEMORY,
                              "failed to allocate batch root transaction");
        return -1;
    }
    size_t node_count_before = alea_vec_count(&sys->nodes);
    for (size_t i = 0; i < request_count; i++) {
        roots[i] = results[i].root_node_id;
        if (!results[i].changed) continue;
        roots[i] = materialize_cell_analysis(sys, &analyses[i]);
        if (roots[i] == ALEA_NODE_ID_INVALID) {
            sys->nodes.count = node_count_before;
            free(roots);
            batch_free_analyses(analyses, request_count);
            alea_set_error_detail(ALEA_ERR_OUT_OF_MEMORY,
                                  "failed to materialize batch request %zu", i);
            return -1;
        }
    }
    for (size_t i = 0; i < request_count; i++) {
        if (!results[i].changed) continue;
        alea_cell_entry_t* cell = &sys->cells.data[requests[i].cell_index];
        cell->root_node_id = roots[i];
        cell->original_root_node_id = ALEA_NODE_ID_INVALID;
        results[i].root_node_id = roots[i];
        results[i].applied = true;
        summary->applied_cells++;
    }
    alea_system_invalidate_query_caches(sys, ALEA_CACHE_ALL);
    free(roots);
    batch_free_analyses(analyses, request_count);
    return 0;
}
