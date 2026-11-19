// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
// SPDX-License-Identifier: MPL-2.0

#include "alea_regional.h"
#include "alea.h"
#include "core/alea_eval.h"
#include "primitives/bbox.h"
#include "primitives/primitive_desc.h"
#include "util/alea_alloc.h"
#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

/* Two closed half-boxes per root. Unprofitable regions have no program.
 * Prefix instructions carry their length so the VM can short-circuit in
 * original operand order. A REF always names an original editable node. */
enum { REF, INSIDE, OUTSIDE, BOTH, EITHER };
typedef struct {
    uint32_t node;
    uint16_t length, op;
} instruction_t;

enum { MAX_DEPTH = 48, MAX_VISITS = 511, MAX_ROOTS = 4096,
       MAX_CODE = 512 * 1024, BUILD_WORK = 2 * 1024 * 1024 };

typedef struct {
    uint32_t key, program[2]; /* Root and program indices plus one. */
    unsigned axis;
    double middle;
} region_entry_t;

struct alea_regional_cache {
    region_entry_t* table;
    size_t slots, roots, count, capacity;
    instruction_t* code;
};

typedef struct {
    instruction_t code[MAX_VISITS];
    size_t count, visits;
} compiler_t;

static region_entry_t* lookup(const struct alea_regional_cache* c, uint32_t root) {
    size_t i = (root * 2654435761u) & (c->slots - 1);
    while (c->table[i].key && c->table[i].key != root + 1)
        i = (i + 1) & (c->slots - 1);
    return &c->table[i];
}

/* Choose extrema of the affine expression over the closed domain. Use the
 * scalar primitive kernel and leave a conservative margin for rounding
 * (three products and three additions, including possible contraction). */
static int plane_relation(const alea_system_t* sys, const alea_node_t* node,
                          const alea_bbox_t* b) {
    uint32_t id = node->primitive.primitive_id;
    if (id >= sys->primitives.count || sys->primitives.data[id].type != ALEA_PRIMITIVE_PLANE)
        return REF;
    const alea_plane_data_t* p = alea_primitive_payload_const(sys, id);
    if (!p) return REF;
    double scale = fabs(p->d) + fabs(p->a) * fmax(fabs(b->min_x), fabs(b->max_x))
                              + fabs(p->b) * fmax(fabs(b->min_y), fabs(b->max_y))
                              + fabs(p->c) * fmax(fabs(b->min_z), fabs(b->max_z));
    if (!isfinite(scale)) return REF;
    double margin = 32 * (scale * DBL_EPSILON + DBL_MIN);
    bool flip = (node->primitive.sense > 0) != (node->primitive.inverted != 0);
    double lo = alea_primitive_eval_payload(ALEA_PRIMITIVE_PLANE, p,
        p->a < 0 ? b->max_x : b->min_x, p->b < 0 ? b->max_y : b->min_y,
        p->c < 0 ? b->max_z : b->min_z);
    if (lo > margin) return flip ? INSIDE : OUTSIDE;
    double hi = alea_primitive_eval_payload(ALEA_PRIMITIVE_PLANE, p,
        p->a < 0 ? b->min_x : b->max_x, p->b < 0 ? b->min_y : b->max_y,
        p->c < 0 ? b->min_z : b->max_z);
    if (hi < -margin) return flip ? OUTSIDE : INSIDE;
    return REF;
}

/* The trilinear interpolant of a quadric matches its constant, linear and
 * mixed terms exactly. For a diagonal term A*x*x the interpolant exceeds
 * the polynomial by A*(x-lo)*(hi-x), at most max(A,0)*half_width^2.
 * Thus all eight signed corners above the correction imply strict exclusion,
 * including when the minimum is on a face or in the interior. */
static bool quadric_proves_outside(const alea_quadric_data_t* q, double sign,
                                   const alea_bbox_t* box) {
    if (!q) return false;
    const double lo[] = {box->min_x, box->min_y, box->min_z};
    const double hi[] = {box->max_x, box->max_y, box->max_z};
    double coord = 1.0, coefficient = 1.0;
    for (unsigned a = 0; a < 3; a++)
        coord = fmax(coord, fmax(fabs(lo[a]), fabs(hi[a])));
    for (unsigned i = 0; i < 10; i++) {
        if (!isfinite(q->coeffs[i])) return false;
        coefficient = fmax(coefficient, fabs(q->coeffs[i]));
    }

    /* Each of the ten monomials has magnitude <= scale. At most two
     * products per term and nine additions cost <128 eps*scale per scalar
     * evaluation. Cover both corner and query rounding, the three curvature
     * corrections and contraction. Flooring each scale factor at one also
     * covers underflow. Decline overflow-prone sums before evaluating them. */
    double scale = coefficient * coord * coord;
    if (!(scale < DBL_MAX / 64)) return false;
    double threshold = 512 * DBL_EPSILON * scale;
    for (unsigned a = 0; a < 3; a++) {
        double radius = (hi[a] - lo[a]) * 0.5;
        threshold += fmax(0.0, sign * q->coeffs[a]) * radius * radius;
    }
    for (unsigned corner = 0; corner < 8; corner++) {
        double value = alea_primitive_eval_payload(ALEA_PRIMITIVE_QUADRIC, q,
            (corner & 1) ? hi[0] : lo[0], (corner & 2) ? hi[1] : lo[1],
            (corner & 4) ? hi[2] : lo[2]);
        if (!(sign * value > threshold)) return false;
    }
    return true;
}

static bool proves_outside(const alea_system_t* sys, alea_node_id_t id,
                           const alea_bbox_t* box, unsigned depth, size_t* remaining) {
    if (!*remaining || depth == MAX_DEPTH || id >= sys->nodes.count) return false;
    --*remaining;
    const alea_node_t* n = &sys->nodes.data[id];
    const alea_node_bbox_t* b = &n->bbox;
    if (box->max_x < b->min_x || box->min_x > b->max_x ||
        box->max_y < b->min_y || box->min_y > b->max_y ||
        box->max_z < b->min_z || box->min_z > b->max_z) return true;

    /* A partial bbox overlap can only turn more points outside. Positive
     * proofs remain valid there. Difference can safely use its left operand. */
    switch (ALEA_GET_OPERATION(n)) {
    case ALEA_OP_PRIMITIVE: {
        uint32_t id = n->primitive.primitive_id;
        if (id < sys->primitives.count && sys->primitives.data[id].type == ALEA_PRIMITIVE_QUADRIC) {
            bool flip = (n->primitive.sense > 0) != (n->primitive.inverted != 0);
            return quadric_proves_outside(alea_primitive_payload_const(sys, id),
                                          flip ? -1.0 : 1.0, box);
        }
        return plane_relation(sys, n, box) == OUTSIDE;
    }
    case ALEA_OP_UNION:
        return proves_outside(sys, n->operation.left, box, depth + 1, remaining) &&
               proves_outside(sys, n->operation.right, box, depth + 1, remaining);
    case ALEA_OP_INTERSECTION:
        if (proves_outside(sys, n->operation.right, box, depth + 1, remaining)) return true;
        /* fall through */
    case ALEA_OP_DIFFERENCE:
        return proves_outside(sys, n->operation.left, box, depth + 1, remaining);
    default:
        return false;
    }
}

bool alea_regional_proves_outside(const alea_system_t* sys, alea_node_id_t root,
                                  const alea_bbox_t* box, size_t* remaining) {
    return proves_outside(sys, root, box, 0, remaining);
}

static bool compile(compiler_t* c, const alea_system_t* sys, uint32_t id,
                    const alea_bbox_t* box, unsigned depth) {
    if (++c->visits > MAX_VISITS) return false;
    size_t start = c->count++;
    instruction_t* out = &c->code[start];
    *out = (instruction_t){id, 1, REF};
    if (id >= sys->nodes.count || depth == MAX_DEPTH) return true;
    const alea_node_t* n = &sys->nodes.data[id];
    const alea_node_bbox_t* b = &n->bbox;
    if (box->max_x < b->min_x || box->min_x > b->max_x ||
        box->max_y < b->min_y || box->min_y > b->max_y ||
        box->max_z < b->min_z || box->min_z > b->max_z) {
        out->op = OUTSIDE;
        return true;
    }
    /* The full evaluator checks each node's bbox before its operation. */
    if (!alea_node_bbox_contains_point(b, box->min_x, box->min_y, box->min_z) ||
        !alea_node_bbox_contains_point(b, box->max_x, box->max_y, box->max_z))
        return true;
    alea_operation_t op = ALEA_GET_OPERATION(n);
    if (op == ALEA_OP_PRIMITIVE) {
        out->op = plane_relation(sys, n, box);
        return true;
    }
    if (op != ALEA_OP_INTERSECTION && op != ALEA_OP_UNION) return true;

    size_t left = c->count;
    if (!compile(c, sys, n->operation.left, box, depth + 1)) return false;
    int absorb = op == ALEA_OP_INTERSECTION ? OUTSIDE : INSIDE;
    int neutral = op == ALEA_OP_INTERSECTION ? INSIDE : OUTSIDE;
    int lop = c->code[left].op;
    if (lop == absorb) {
        out->op = absorb;
        c->count = start + 1;
        return true;
    }

    size_t right = c->count;
    if (!compile(c, sys, n->operation.right, box, depth + 1)) return false;
    int rop = c->code[right].op;
    if (rop == absorb) {
        out->op = absorb;
        c->count = start + 1;
    } else if (lop == neutral || rop == neutral) {
        size_t keep = lop == neutral ? right : left;
        size_t length = c->code[keep].length;
        memmove(out, &c->code[keep], length * sizeof(*out));
        c->count = start + length;
    } else if (lop == REF && c->code[left].node == n->operation.left &&
               rop == REF && c->code[right].node == n->operation.right) {
        c->count = start + 1;
    } else {
        out->op = op == ALEA_OP_INTERSECTION ? BOTH : EITHER;
        out->length = (uint16_t)(c->count - start);
    }
    return true;
}

void alea_regional_free(struct alea_regional_cache* c) {
    if (c) {
        free(c->code);
        free(c->table);
        free(c);
    }
}

void alea_regional_stats(const alea_system_t* sys, size_t* roots, size_t* bytes) {
    const struct alea_regional_cache* c = sys->regional_cache;
    *roots = c ? c->roots : 0;
    *bytes = c ? sizeof(*c) + c->slots * sizeof(*c->table) +
                 c->capacity * sizeof(*c->code) : 0;
}

static int save_program(struct alea_regional_cache* c, const compiler_t* p,
                        region_entry_t* entry, unsigned half) {
    if (half && entry->program[0]) {
        const instruction_t* prev = c->code + entry->program[0] - 1;
        if (prev->length == p->count && !memcmp(prev, p->code, p->count * sizeof(*prev))) {
            entry->program[1] = entry->program[0];
            return 0;
        }
    }
    size_t need = c->count + p->count;
    if (need > MAX_CODE) return 0;
    if (need > c->capacity) {
        size_t cap = c->capacity ? 2 * c->capacity : MAX_VISITS;
        if (cap > MAX_CODE) cap = MAX_CODE;
        instruction_t* data = alea_realloc(c->code, cap * sizeof(*data));
        if (!data) return -1;
        c->code = data;
        c->capacity = cap;
    }
    entry->program[half] = (uint32_t)c->count + 1;
    memcpy(c->code + c->count, p->code, p->count * sizeof(*p->code));
    c->count = need;
    return 0;
}

int alea_regional_build(alea_system_t* sys) {
    struct alea_regional_cache* c = alea_calloc(1, sizeof(*c));
    if (!c) goto oom;
    size_t limit = sys->cells.count < MAX_ROOTS ? sys->cells.count : MAX_ROOTS;
    c->slots = 2;
    while (c->slots < 2 * limit) c->slots *= 2;
    c->table = alea_calloc(c->slots, sizeof(*c->table));
    if (!c->table) goto oom;

    alea_bbox_t unbounded = alea_bbox_infinite();
    size_t work = 0;
    for (size_t i = 0; i < sys->cells.count && c->roots < limit && work < BUILD_WORK; i++) {
        work++;
        if ((i & 255) == 0 && alea_interrupted()) {
            alea_regional_free(c);
            alea_set_error_detail(ALEA_ERR_INTERRUPTED, "Preparing regional expressions interrupted");
            return -1;
        }
        uint32_t root = sys->cells.data[i].root_node_id;
        if (root >= sys->nodes.count) continue;
        const alea_node_t* n = &sys->nodes.data[root];
        alea_operation_t op = ALEA_GET_OPERATION(n);
        if (op != ALEA_OP_INTERSECTION && op != ALEA_OP_UNION) continue;
        alea_bbox_t box = alea_node_bbox_get(&n->bbox);
        double lo[] = {box.min_x, box.min_y, box.min_z};
        double hi[] = {box.max_x, box.max_y, box.max_z};
        unsigned axis = 0;
        bool finite = true;
        for (unsigned k = 0; k < 3; k++) {
            if (!(lo[k] > unbounded.min_x && hi[k] < unbounded.max_x && hi[k] > lo[k]))
                finite = false;
            if (hi[k] - lo[k] > hi[axis] - lo[axis]) axis = k;
        }
        if (!finite) continue;
        region_entry_t* entry = lookup(c, root);
        if (entry->key) continue;
        double middle = (lo[axis] + hi[axis]) * 0.5;
        for (unsigned half = 0; half < 2; half++) {
            alea_bbox_t region = box;
            double* low[] = {&region.min_x, &region.min_y, &region.min_z};
            double* high[] = {&region.max_x, &region.max_y, &region.max_z};
            *(half ? low[axis] : high[axis]) = middle;
            compiler_t p;
            p.visits = p.count = 0;
            bool ok = compile(&p, sys, root, &region, 0);
            work += p.visits;
            if (ok && (p.code[0].op != REF || p.code[0].node != root) &&
                save_program(c, &p, entry, half) != 0) goto oom;
        }
        if (entry->program[0] || entry->program[1]) {
            entry->key = root + 1;
            entry->axis = axis;
            entry->middle = middle;
            c->roots++;
        }
    }
    if (!c->roots) {
        alea_regional_free(c);
        c = NULL;
    }
    sys->regional_cache = c;
    return 0;
oom:
    alea_regional_free(c);
    alea_set_error_detail(ALEA_ERR_OUT_OF_MEMORY, "Preparing regional expressions");
    return -1;
}

/* An opaque signed subtree can return NaN. The original evaluator's ordered
 * min/max then matters, so propagate uncertainty to the full root. */
static int evaluate(const instruction_t* p, const alea_system_t* sys,
                    double x, double y, double z) {
    switch (p->op) {
    case INSIDE: return 1;
    case OUTSIDE: return 0;
    case REF: {
        double v = alea_evaluate_point(sys, p->node, x, y, z);
        return isnan(v) ? -1 : v <= 0;
    }
    default: {
        int left = evaluate(p + 1, sys, x, y, z);
        if (left < 0 || left == (p->op == EITHER)) return left;
        return evaluate(p + 1 + p[1].length, sys, x, y, z);
    }
    }
}

bool alea_regional_contains(const alea_system_t* sys, alea_node_id_t root,
                            double x, double y, double z, bool* inside) {
    const struct alea_regional_cache* c = sys->regional_cache;
    if (!c || root >= sys->nodes.count) return false;
    /* Use the full evaluator's first rejection before paying for a lookup. */
    if (!alea_node_bbox_contains_point(&sys->nodes.data[root].bbox, x, y, z)) {
        *inside = false;
        return true;
    }
    const region_entry_t* e = lookup(c, root);
    if (!e->key) return false;
    double coord = e->axis == 0 ? x : e->axis == 1 ? y : z;
    uint32_t program = e->program[coord > e->middle];
    if (!program) return false;
    int result = evaluate(c->code + program - 1, sys, x, y, z);
    if (result < 0) return false;
    *inside = result != 0;
    return true;
}
