// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
// SPDX-License-Identifier: MPL-2.0
/* Conservative emptiness tests for the untraced part of a camera or shadow
 * ray. Treat either side of a boundary as possible until strictly excluded. */
#include "raycast/ray_bbox.h"
#include "primitives/primitive_codec.h"

static bool render_box_interval(const alea_ray_t* ray, alea_bbox_t box,
                                double lo, double hi, double* enter, double* exit) {
    if (box.min_x > box.max_x) return false;
    const double pad = 128 * DBL_EPSILON * (1 + fabs(ray->ox) + fabs(ray->oy) +
        fabs(ray->oz) + fabs(box.min_x) + fabs(box.max_x) + fabs(box.min_y) +
        fabs(box.max_y) + fabs(box.min_z) + fabs(box.max_z));
    box.min_x -= pad; box.max_x += pad;
    box.min_y -= pad; box.max_y += pad;
    box.min_z -= pad; box.max_z += pad;
    return ray_bbox_slab_enter_exit(ray, &box, lo, hi, enter, exit);
}

/* A one-sheet cone evaluates to +1 on the rejected sheet, independently
 * of the double cone's polynomial there. */
static int render_cone_sign(const alea_system_t* sys, const alea_node_t* n, interval_box_t box) {
    alea_primitive_data_t data;
    if (!alea_primitive_copy_data(sys, n->primitive.primitive_id, &data)) return 0;
    unsigned axis = n->primitive.prim_type - ALEA_PRIMITIVE_CONE_X;
    /* The cone codec exports apex x/y/z, squared slope and sheet. */
    double values[5];
    if (alea_primitive_to_values(n->primitive.prim_type, &data, values) != 5) return 0;
    const double* apex = values;
    double slope = values[3], sheet = values[4];
    if (!(slope >= 0) || !isfinite(slope)) return 0;
    if ((sheet > 0 && box.hi[axis] < apex[axis]) ||
        (sheet < 0 && box.lo[axis] > apex[axis])) return 1;
    double near2[3], far2[3];
    for (unsigned k = 0; k < 3; ++k) {
        double a = box.lo[k]-apex[k], b = box.hi[k]-apex[k];
        double near = fmax(0, fmax(a, -b)), far = fmax(fabs(a), fabs(b));
        near2[k] = near*near;
        far2[k] = far*far;
    }
    unsigned j=(axis+1)%3, k=(axis+2)%3;
    double slack=256*(DBL_EPSILON*(far2[j]+far2[k]+slope*far2[axis])+DBL_MIN);
    if (!isfinite(slack)) return 0;
    if (near2[j]+near2[k]-slope*far2[axis] > slack) return 1;
    if ((sheet > 0 && box.lo[axis] < apex[axis]) ||
        (sheet < 0 && box.hi[axis] > apex[axis])) return 0;
    return far2[j]+far2[k]-slope*near2[axis] < -slack ? -1 : 0;
}

/* Only prove that the nonpositive side is absent. A difference is contained
 * in its left operand. Complements and unsupported primitives stay possible.
 * Iterate left-deep unions without consuming stack for every term. */
static bool render_node_excluded(const alea_system_t* sys, uint32_t id,
                                 interval_box_t box, unsigned depth, unsigned* budget) {
    while (depth < 128 && *budget && id < sys->nodes.count) {
        --*budget;
        const alea_node_t* n = sys->nodes.data + id;
        const alea_bbox_t bound = alea_node_bbox_get(&n->bbox);
        if (box.hi[0] < bound.min_x || box.lo[0] > bound.max_x ||
            box.hi[1] < bound.min_y || box.lo[1] > bound.max_y ||
            box.hi[2] < bound.min_z || box.lo[2] > bound.max_z) return true;
        switch (ALEA_GET_OPERATION(n)) {
        case ALEA_OP_PRIMITIVE: {
            guard_surface_t surface = {0};
            guard_surface_init(&surface, sys, n->primitive.primitive_id, true);
            int sign = surface.kind ? guard_box_sign(&surface, box) : 0;
            if (n->primitive.prim_type >= ALEA_PRIMITIVE_CONE_X &&
                n->primitive.prim_type <= ALEA_PRIMITIVE_CONE_Z) sign = render_cone_sign(sys, n, box);
            bool flipped = (n->primitive.sense > 0) != (n->primitive.inverted != 0);
            return flipped ? sign < 0 : sign > 0;
        }
        case ALEA_OP_UNION:
            if (!render_node_excluded(sys, n->operation.right, box, depth + 1, budget)) return false;
            break;
        case ALEA_OP_INTERSECTION:
            if (render_node_excluded(sys, n->operation.right, box, depth + 1, budget)) return true;
            break;
        case ALEA_OP_DIFFERENCE: break;
        default: return false;
        }
        id = n->operation.left;
    }
    return false;
}

static interval_box_t render_ray_box(const alea_ray_t* ray, double lo, double hi) {
    double a[3], b[3];
    alea_ray_point_at(ray, lo, a, a+1, a+2);
    alea_ray_point_at(ray, hi, b, b+1, b+2);
    const double pad = 64 * DBL_EPSILON *
        (1 + fabs(ray->ox) + fabs(ray->oy) + fabs(ray->oz) + hi);
    interval_box_t box;
    for (unsigned k = 0; k < 3; ++k) {
        box.lo[k] = fmin(a[k], b[k]) - pad;
        box.hi[k] = fmax(a[k], b[k]) + pad;
    }
    return box;
}

/* Subdivide long tails geometrically so distant endpoints do not obscure
 * short intervals near the model. */
static double render_split_distance(double lo, double hi) {
    double near = fmax(1.0, lo);
    return hi > 4 * near ? sqrt(near) * sqrt(hi) : lo + (hi-lo)*0.5;
}

static bool render_cell_tail_empty(const alea_system_t* sys, const render_bound_t* bound,
                                   const alea_ray_t* ray, double lo, double hi,
                                   unsigned depth, unsigned* budget) {
    interval_box_t box = render_ray_box(ray, lo, hi);
    if (render_node_excluded(sys, bound->branch, box, 0, budget) ||
        (bound->cell_root != bound->branch && render_node_excluded(sys, bound->cell_root, box, 0, budget))) return true;
    if (!*budget || depth >= 32) return false;
    double mid = render_split_distance(lo, hi);
    return mid > lo && mid < hi &&
        render_cell_tail_empty(sys, bound, ray, lo, mid, depth+1, budget) &&
        render_cell_tail_empty(sys, bound, ray, mid, hi, depth+1, budget);
}

static bool render_node_encloses_box(const alea_node_t* node, interval_box_t box) {
    return alea_node_bbox_contains_point(&node->bbox, box.lo[0], box.lo[1], box.lo[2]) &&
           alea_node_bbox_contains_point(&node->bbox, box.hi[0], box.hi[1], box.hi[2]);
}

/* At points not owned by the leading exterior void, the opposite side of
 * its sphere is known. Other literals must be proven over the entire box.
 * Check every node's own bounds just as the scalar evaluator does. */
static bool render_void_contains(const alea_system_t* sys, uint32_t id,
                                 const alea_node_t* exterior, interval_box_t box,
                                 unsigned depth, unsigned* work) {
    if (!*work || depth >= 128 || id >= sys->nodes.count) return false;
    --*work;
    const alea_node_t* n = sys->nodes.data + id;
    if (!render_node_encloses_box(n, box)) return false;
    switch (ALEA_GET_OPERATION(n)) {
    case ALEA_OP_PRIMITIVE: {
        bool flip = (n->primitive.sense > 0) != (n->primitive.inverted != 0);
        if (!flip && n->primitive.primitive_id == exterior->primitive.primitive_id) return true;
        guard_surface_t surface = {0};
        guard_surface_init(&surface, sys, n->primitive.primitive_id, true);
        int sign = surface.kind ? guard_box_sign(&surface, box) : 0;
        return flip ? sign > 0 : sign < 0;
    }
    case ALEA_OP_INTERSECTION:
        return render_void_contains(sys, n->operation.left, exterior, box, depth+1, work) &&
               render_void_contains(sys, n->operation.right, exterior, box, depth+1, work);
    case ALEA_OP_UNION:
        return render_void_contains(sys, n->operation.left, exterior, box, depth+1, work) ||
               render_void_contains(sys, n->operation.right, exterior, box, depth+1, work);
    default:
        return false;
    }
}

static bool render_void_covers_ray(const interval_index_t* idx, const alea_universe_t* world,
                                   const alea_ray_t* ray, double lo, double hi,
                                   unsigned depth, unsigned* work) {
    interval_box_t box = render_ray_box(ray, lo, hi);
    const alea_node_t* sphere = idx->sys->nodes.data + idx->render_exterior_void_root;
    if (!render_node_encloses_box(sphere, box)) return false;
    guard_surface_t surface = {0};
    guard_surface_init(&surface, idx->sys, sphere->primitive.primitive_id, true);
    if (guard_box_sign(&surface, box) > 0) return true;
    for (size_t i = 0; i < idx->render_void_prefix_count && *work; ++i) {
        uint32_t root = idx->sys->cells.data[world->cell_indices.data[i]].root_node_id;
        if (root == idx->render_exterior_void_root) {
            --*work;
            continue;
        }
        if (render_void_contains(idx->sys, root, sphere, box, 0, work)) return true;
    }
    if (!*work || depth >= 32) return false;
    --*work;
    double mid = render_split_distance(lo, hi);
    return mid > lo && mid < hi &&
        render_void_covers_ray(idx, world, ray, lo, mid, depth+1, work) &&
        render_void_covers_ray(idx, world, ray, mid, hi, depth+1, work);
}

static bool render_tail_empty(const interval_index_t* idx, const alea_ray_t* ray,
                              double lo, double hi) {
    if (!idx->render_bounds_ready || !isfinite(hi) ||
        render_box_interval(ray, idx->render_core, lo, hi, NULL, NULL)) return false;
    if (idx->render_void_prefix_count && idx->render_exterior_void_root != UINT32_MAX) {
        const alea_node_t* sphere = idx->sys->nodes.data + idx->render_exterior_void_root;
        alea_bbox_t bounds = alea_node_bbox_get(&sphere->bbox);
        double exit;
        if (ray_bbox_slab_enter_exit(ray, &bounds, lo, hi, NULL, &exit)) {
            /* Leave room for the padded ray box at the edge of the node.
             * The containment proof checks this choice before we advance. */
            double end = exit == hi ? hi : exit - 1e-10 * (1 + fabs(exit));
            unsigned work = 2048;
            const alea_universe_t* world = alea_get_universe(idx->sys, 0);
            if (world && end > lo &&
                render_void_covers_ray(idx, world, ray, lo, end, 0, &work)) {
                if (end == hi) return true;
                lo = end;
            }
        }
    }
    unsigned budget = 100000;
    for (size_t i = 0; i < idx->render_boxes.count; ++i) {
        double enter, exit;
        alea_bbox_t box = alea_node_bbox_get(&idx->render_boxes.data[i].box);
        if (render_box_interval(ray, box, lo, hi, &enter, &exit) &&
            !render_cell_tail_empty(idx->sys, idx->render_boxes.data+i, ray,
                                    enter, exit, 0, &budget)) return false;
    }
    return true;
}
