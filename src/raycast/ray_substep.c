// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
// SPDX-License-Identifier: MPL-2.0
#include "ray_substep.h"
#include "ray_traversal.h"
#include "ray_wide.h"
#include "ray_box_sign.h"
#include "primitives/bbox.h"
#include "core/alea_universe.h"
#include <float.h>

/* Enclose O+tD for a short neighborhood of the hint in each local frame.
 * fma rounds each endpoint once. Expand outwards to include that rounding. */
static alea_bbox_t neighborhood(const alea_ray_t* r, double t, double radius) {
    const double o[] = {r->ox, r->oy, r->oz}, d[] = {r->dx, r->dy, r->dz};
    double lo[3], hi[3];
    for (int i = 0; i < 3; ++i) {
        double a = fma(t - radius, d[i], o[i]), b = fma(t + radius, d[i], o[i]);
        lo[i] = nextafter(fmin(a, b), -INFINITY);
        hi[i] = nextafter(fmax(a, b), INFINITY);
    }
    return (alea_bbox_t){lo[0], hi[0], lo[1], hi[1], lo[2], hi[2]};
}

/* Keep the low parts of O + tD when a STRICT interior sample rounds onto a face. */
typedef struct {
    ray_wide_t xyz[3];
    double scale[3];
} ray_sample_t;

static int sample_plane_sign(const ray_sample_t* point, const alea_plane_data_t* p) {
    if (!p) return 0;
    const double normal[] = {p->a, p->b, p->c};
    ray_wide_t value = {p->d, 0};
    double scale = fabs(p->d);
    for (int i = 0; i < 3; ++i) {
        value = ray_wide_term(normal[i], point->xyz[i], value);
        scale += fabs(normal[i]) * point->scale[i];
    }
    double error = 256 * (DBL_EPSILON * DBL_EPSILON * scale + DBL_MIN);
    if (!isfinite(error) || !isfinite(value.hi) || !isfinite(value.lo)) return 0;
    return value.hi > error ? 1 : value.hi < -error ? -1 : 0;
}

typedef struct {
    alea_bbox_t box;
    uint32_t first, second;
    int first_sign, second_sign;
    const ray_sample_t* point;
} slab_t;

/* Cone interval evaluation must include the sheet test. The discarded
 * half of a one-sheet cone has positive sign regardless of radius. */
static int cone_sign(const alea_system_t* sys, uint32_t id, interval_box_t b) {
    alea_primitive_data_t data;
    if (!alea_primitive_copy_data(sys, id, &data)) return 0;
    /* The three union members share their entire initial sequence. */
    const alea_cone_x_data_t* c = &data.cone_x;
    const double apex[] = {c->apex_x, c->apex_y, c->apex_z};
    unsigned axis = sys->primitives.data[id].type - ALEA_PRIMITIVE_CONE_X;
    int sheet = c->sheet_selection;
    if ((sheet > 0 && b.hi[axis] < apex[axis]) ||
        (sheet < 0 && b.lo[axis] > apex[axis])) return 1;
    double lo = 0, hi = 0, magnitude = 0;
    for (unsigned i = 0; i < 3; ++i) {
        double pad = 64 * DBL_EPSILON * (1 + fabs(apex[i]) + fabs(b.lo[i]) + fabs(b.hi[i]));
        double near = fmax(0, fmax(b.lo[i] - apex[i], apex[i] - b.hi[i])-pad);
        double far = fmax(fabs(b.lo[i] - apex[i]), fabs(b.hi[i] - apex[i])) + pad;
        if (i == axis) {
            lo -= c->tan_angle_sq * far * far;
            hi -= c->tan_angle_sq * near * near;
            magnitude += c->tan_angle_sq * far * far;
        } else {
            lo += near * near;
            hi += far * far;
            magnitude += far * far;
        }
    }
    double error = 64 * (DBL_EPSILON * magnitude + DBL_MIN);
    if (lo > error) return 1;
    if ((sheet > 0 && b.lo[axis] < apex[axis]) ||
        (sheet < 0 && b.hi[axis] > apex[axis])) return 0;
    return hi < -error ? -1 : 0;
}

/* CSG signs on an open slab. Its bounding surfaces have known signs even
 * when a Cartesian binary64 enclosure straddles both of them. Other leaves
 * must resolve over the whole enclosure or remain unknown (zero). */
static int cell_sign(const alea_system_t* sys, alea_node_id_t id,
                     const slab_t* slab, unsigned depth) {
    if (id >= sys->nodes.count || depth > 256) return 0;
    const alea_node_t* n = &sys->nodes.data[id];
    alea_operation_t op = ALEA_GET_OPERATION(n);
    if (op == ALEA_OP_PRIMITIVE) {
        uint32_t p = n->primitive.primitive_id;
        if (p >= sys->primitives.count) return 0;
        int sign;
        if (p == slab->first || p == slab->second) {
            sign = p == slab->first ? slab->first_sign : slab->second_sign;
        } else {
            alea_primitive_type_t type = sys->primitives.data[p].type;
            interval_box_t box = interval_box_from(slab->box);
            if (slab->point && type == ALEA_PRIMITIVE_PLANE) {
                sign = sample_plane_sign(slab->point, alea_primitive_payload_const(sys, p));
            } else if (type >= ALEA_PRIMITIVE_CONE_X && type <= ALEA_PRIMITIVE_CONE_Z) {
                sign = cone_sign(sys, p, box);
            } else {
                guard_surface_t guard = {0};
                guard_surface_init(&guard, sys, p, false);
                sign = guard.kind ? guard_box_sign(&guard, box) : 0;
            }
        }
        return (n->primitive.sense > 0) != !!n->primitive.inverted ? -sign : sign;
    }
    int a = cell_sign(sys, n->operation.left, slab, depth + 1);
    if (op == ALEA_OP_COMPLEMENT) return -a;
    if (op == ALEA_OP_UNION && a < 0) return -1;
    if ((op == ALEA_OP_INTERSECTION || op == ALEA_OP_DIFFERENCE) && a > 0) return 1;
    int b = cell_sign(sys, n->operation.right, slab, depth + 1);
    if (op == ALEA_OP_UNION) return a == 0 ? 0 : b;
    if (op == ALEA_OP_DIFFERENCE) b = -b;
    if (op != ALEA_OP_DIFFERENCE && op != ALEA_OP_INTERSECTION) return 0;
    return a > b ? a : b;
}

static const alea_plane_data_t* as_plane(const alea_system_t* s, uint32_t id) {
    if (id >= s->primitives.count || s->primitives.data[id].type != ALEA_PRIMITIVE_PLANE)
        return NULL;
    return alea_primitive_payload_const(s, id);
}

/* Cancel O+tD and the plane offset with two-component arithmetic. Reject
 * ill-conditioned divisions and roots too close to a rounding boundary. */
static bool rounded_cut(const alea_ray_t* r, const alea_plane_data_t* p,
                        double t, double slope, double* cut) {
    double normal[] = {p->a, p->b, p->c}, o[] = {r->ox, r->oy, r->oz}, d[] = {r->dx, r->dy, r->dz};
    ray_wide_t f = {p->d, 0};
    double scale = fabs(p->d), slope_scale = 0;
    for (int i = 0; i < 3; ++i) {
        ray_wide_t x = ray_wide_term(d[i], (ray_wide_t){t, 0}, (ray_wide_t){o[i], 0});
        f = ray_wide_term(normal[i], x, f);
        scale += fabs(normal[i]) * (fabs(o[i]) + fabs(t * d[i]));
        slope_scale += fabs(normal[i] * d[i]);
    }
    double slope_error = 8 * DBL_EPSILON * slope_scale + DBL_MIN;
    if (!(fabs(slope) > 2 * slope_error)) return false;
    double shift = -(f.hi + f.lo) / slope;
    double error = (128 * DBL_EPSILON * DBL_EPSILON * scale + DBL_MIN +
                    2 * fabs(shift) * slope_error) / fabs(slope) + 4 * DBL_EPSILON * fabs(shift);
    double low = t + (shift - error), high = t + (shift + error);
    if (!isfinite(low) || low != high) return false;
    *cut = low;
    return true;
}

bool alea_substep_plane_cuts(const alea_ray_t* r, const alea_plane_data_t* p,
    const alea_plane_data_t* q, double t, double* enter, double* exit, double* width) {
    int orientation = 1;
    if (p->a != q->a || p->b != q->b || p->c != q->c) {
        if (p->a != -q->a || p->b != -q->b || p->c != -q->c) return false;
        orientation = -1;
    }
    double slope = p->a * r->dx + p->b * r->dy + p->c * r->dz;
    if (!isfinite(slope) || slope == 0) return false;
    double step = nextafter(t, INFINITY) - t;
    *width = (p->d - orientation * q->d) / slope;
    if (!isfinite(step) || !(*width > 0 && *width <= 2 * step)) return false;
    return rounded_cut(r, p, t, slope, enter) && fabs(*enter - t) <= 4 * step &&
           rounded_cut(r, q, t, orientation * slope, exit) &&
           *exit >= *enter && *exit <= nextafter(*enter, INFINITY);
}

/* Return the only cell containing the open slab, -1 for an empty
 * universe, or -2 when any claimant or exclusion is indeterminate. */
static int unique_slab_cell(const alea_system_t* sys, const alea_universe_t* u,
                            const slab_t* slab) {
    if (!u) return -2;
    int owner = -1;
    for (size_t i = 0; i < u->cell_indices.count; ++i) {
        uint32_t index = u->cell_indices.data[i];
        const alea_cell_entry_t* c = &sys->cells.data[index];
        /* A repeated lattice can own space outside its zero-element CSG. */
        if (c->lat_type || c->root_node_id >= sys->nodes.count) return -2;
        /* Only the original outward-rounded cell bound can exclude a
         * claimant here. Placement bounds may be clipped by ancestors. */
        alea_bbox_t box = alea_node_bbox_get(&sys->nodes.data[c->root_node_id].bbox);
        if (alea_bbox_is_valid(&box) && !alea_bbox_intersects(&box, &slab->box)) continue;
        int sign = cell_sign(sys, c->root_node_id, slab, 0);
        if (sign == 0 || (sign < 0 && owner >= 0)) return -2;
        if (sign < 0) owner = (int)index;
    }
    return owner;
}

bool alea_ray_verify_short_owner(const alea_system_t* sys, const alea_ray_t* ray,
        double lo, double hi, double sample_t, const alea_hier_ray_path_t* path) {
    const double width = hi - lo, half = 0.5 * width;
    /* Bound the cost of the exhaustive search. Width does not decide ownership. */
    if (!path || !path->count || !isfinite(hi) || !(half > 0 && half < width) ||
        width > 4 * (hi - nextafter(hi, lo))) return false;

    const double origin[] = {ray->ox, ray->oy, ray->oz};
    const double direction[] = {ray->dx, ray->dy, ray->dz};
    ray_sample_t point;
    /* Use a normalized two-component parameter even when its low part
     * cannot be represented by any single distance between lo and hi. */
    ray_wide_t t = {sample_t, 0};
    if (!(sample_t > lo && sample_t < hi))
        t = ray_wide_term(1, (ray_wide_t){half, 0}, (ray_wide_t){lo, 0});
    for (int k = 0; k < 3; ++k) {
        point.xyz[k] = ray_wide_term(direction[k], t, (ray_wide_t){origin[k], 0});
        point.scale[k] = fabs(origin[k]) + fabs(hi * direction[k]);
    }
    slab_t sample = {.box = neighborhood(ray, lo, width), .first = UINT32_MAX,
                     .second = UINT32_MAX, .point = &point};
    int universe = 0;
    for (int k = 0; k < path->count; ++k) {
        int owner = unique_slab_cell(sys, alea_get_universe(sys, universe), &sample);
        if (owner < 0 || (uint32_t)owner != path->entries[k].cell_index) return false;
        const alea_cell_entry_t* cell = &sys->cells.data[owner];
        if (k + 1 == path->count) return !alea_cell_entry_is_container(cell);
        /* Only ordinary fills sharing the world coordinate system can
         * use this witness without propagating transform roundoff. */
        if (cell->fill_universe <= 0 || cell->fill_transform > 0) return false;
        universe = cell->fill_universe;
    }
    return false;
}

static bool load_cylinder(const alea_system_t* sys, uint32_t id, guard_surface_t* cylinder) {
    if (id >= sys->primitives.count) return false;
    alea_primitive_type_t type = sys->primitives.data[id].type;
    if (type < ALEA_PRIMITIVE_CYLINDER_X || type > ALEA_PRIMITIVE_CYLINDER_Z) return false;
    guard_surface_init(cylinder, sys, id, false);
    return cylinder->kind == GUARD_RADIAL;
}

/* Keep the low part of each projected coordinate. The sign of the radial
 * residual can be lost even when the two ray distances are many ULPs apart. */
static int cylinder_side(const guard_surface_t* cylinder, const alea_ray_t* ray, ray_wide_t t) {
    const double o[] = {ray->ox, ray->oy, ray->oz}, d[] = {ray->dx, ray->dy, ray->dz};
    double radius = cylinder->k[3];
    ray_wide_t residual = ray_wide_term(-radius, (ray_wide_t){radius, 0}, (ray_wide_t){0, 0});
    double scale = radius * radius;
    for (unsigned i = 0; i < 3; ++i) {
        if (i == cylinder->axis) continue;
        ray_wide_t x = ray_wide_term(d[i], t, (ray_wide_t){o[i], 0});
        x = ray_wide_term(-1, (ray_wide_t){cylinder->k[i], 0}, x);
        residual = ray_wide_madd(x, x, residual);
        double magnitude = fabs(o[i]) + (fabs(t.hi) + fabs(t.lo)) * fabs(d[i]) + fabs(cylinder->k[i]);
        scale += magnitude * magnitude;
    }
    /* Bound the discarded second-order product and sum terms using the
     * magnitudes before cancellation. Decline overflow or an uncertain sign. */
    double error = 256 * (DBL_EPSILON * DBL_EPSILON * scale + DBL_MIN);
    double value = residual.hi + residual.lo;
    if (!isfinite(value) || !isfinite(error)) return 0;
    return value > error ? 1 : value < -error ? -1 : 0;
}

/* Regional runs have already been partitioned at the cell's candidate
 * crossings. Evaluate the bounding predicates at their midpoint in
 * this cell's frame, including when a parent supplied one of the endpoints. */
int alea_ray_interval_midpoint_sign(const alea_system_t* sys, alea_node_id_t root,
    const alea_ray_t* local, uint32_t first, uint32_t second, double t) {
    if (first == second) {
        /* A parent and child can reference the same plane. Only the
         * occurrence in this cell's frame determines its leaf sign. */
        const alea_plane_data_t* p = as_plane(sys, first);
        if (!p) return 0;
        double slope = p->a * local->dx + p->b * local->dy + p->c * local->dz;
        double cut;
        if (!rounded_cut(local, p, t, slope, &cut) || cut == t) return 0;
        int sign = (t > cut) == (slope > 0) ? 1 : -1;
        slab_t sample = {neighborhood(local, t, 0), first, UINT32_MAX, sign, 0, NULL};
        return cell_sign(sys, root, &sample, 0);
    }
    uint32_t ids[] = {first, second};
    int signs[2];
    for (int i = 0; i < 2; ++i) {
        guard_surface_t cylinder = {0};
        if (!load_cylinder(sys, ids[i], &cylinder)) return 0;
        signs[i] = cylinder_side(&cylinder, local, (ray_wide_t){t, 0});
        if (!signs[i]) return 0;
    }
    slab_t sample = {neighborhood(local, t, 0), first, second, signs[0], signs[1], NULL};
    return cell_sign(sys, root, &sample, 0);
}

static bool same_inverse(const alea_matrix_t* a, const alea_matrix_t* b) {
    for (int i = 0; i < 12; ++i)
        if (a->inv[i] != b->inv[i]) return false;
    return true;
}

static alea_hier_ray_path_entry_t ordinary_entry(const alea_system_t* sys, int index,
                                                  int depth, const alea_matrix_t* frame) {
    const alea_cell_entry_t* c = &sys->cells.data[index];
    return (alea_hier_ray_path_entry_t){
        .cell_index = (uint32_t)index, .cell_id = c->mc_cell_id, .material_id = c->material_id,
        .universe_id = c->universe_id, .fill_universe = c->fill_universe,
        .depth = depth, .transform = *frame
    };
}

static int placed_slab_owner(const alea_system_t* sys, int universe, const alea_ray_t* world,
                             const alea_matrix_t* frame, double mid, double radius,
                             slab_t slab, const alea_matrix_t* first, const alea_matrix_t* second) {
    alea_ray_t local;
    alea_ray_transform_inverse(frame, world, &local);
    slab.box = neighborhood(&local, mid, radius);
    /* A reused primitive can cross elsewhere in a different placement. */
    if (!same_inverse(frame, first)) slab.first = UINT32_MAX;
    if (!same_inverse(frame, second)) slab.second = UINT32_MAX;
    return unique_slab_cell(sys, alea_get_universe(sys, universe), &slab);
}

/* Limit the cost of exhaustive ownership checks. The width alone never
 * decides whether a physical interval exists or which cell contains it. */
static bool within_lattice_proof_window(double lo, double hi) {
    double step = nextafter(hi, INFINITY) - hi;
    return hi > lo && isfinite(step) && hi - lo <= 32 * step;
}

/* The forward sample nominates an element. Verify its entire open
 * interval before replacing the leaf, including every enclosing fill. */
static bool resolve_lattice_planes(const alea_system_t* sys, const alea_ray_t* world,
    uint32_t first, const alea_matrix_t* first_frame, uint32_t second, const alea_matrix_t* second_frame,
    double lo, double hi, alea_hier_ray_path_t* path) {
    int leaf = path->count - 1;
    if (leaf < 1 || !path->entries[leaf - 1].is_lattice || path->entries[leaf].is_lattice ||
        !within_lattice_proof_window(lo, hi)) return false;
    const alea_hier_ray_path_entry_t* element = &path->entries[leaf - 1];
    const alea_cell_entry_t* lattice = &sys->cells.data[element->cell_index];
    const alea_universe_t* u = alea_get_universe(sys, element->universe_id);
    if (!first_frame->has_inverse || !second_frame->has_inverse ||
        lattice->lat_type != 1 || !u || u->cell_indices.count != 1 ||
        (!same_inverse(first_frame, &path->entries[leaf].transform) &&
         !same_inverse(second_frame, &path->entries[leaf].transform))) return false;

    uint32_t ids[] = {first, second};
    const alea_matrix_t* frames[] = {first_frame, second_frame};
    double distances[] = {lo, hi}, cuts[2];
    int signs[2];
    for (int i = 0; i < 2; ++i) {
        const alea_plane_data_t* p = as_plane(sys, ids[i]);
        if (!p) return false;
        alea_ray_t local;
        alea_ray_transform_inverse(frames[i], world, &local);
        double slope = p->a * local.dx + p->b * local.dy + p->c * local.dz;
        if (!rounded_cut(&local, p, distances[i], slope, &cuts[i]) ||
            cuts[i] < nextafter(distances[i], -INFINITY) ||
            cuts[i] > nextafter(distances[i], INFINITY)) return false;
        int after_crossing = slope > 0 ? 1 : -1;
        signs[i] = i == 0 ? after_crossing : -after_crossing;
    }
    if (!(cuts[0] < cuts[1])) return false;
    double radius = hi - lo + 4 * (nextafter(hi, INFINITY) - hi);
    alea_ray_t local;
    alea_ray_transform_inverse(&element->transform, world, &local);
    alea_bbox_t box = neighborhood(&local, lo, radius);
    /* Rectangular indices are monotone on each axis. These two corners
     * bound the element index everywhere in the enclosure. */
    alea_lattice_location_t corners[2];
    if (alea_lattice_locate_point(sys, lattice, box.min_x, box.min_y, box.min_z, &corners[0]) != 1 ||
        alea_lattice_locate_point(sys, lattice, box.max_x, box.max_y, box.max_z, &corners[1]) != 1)
        return false;
    for (int i = 0; i < 2; ++i) {
        const alea_lattice_location_t* c = &corners[i];
        if (c->i != element->lat_i || c->j != element->lat_j || c->k != element->lat_k ||
            c->fill_universe != path->entries[leaf].universe_id) return false;
    }
    slab_t slab = {.first=first, .second=second, .first_sign=signs[0], .second_sign=signs[1]};
    for (int i = 0; i < leaf - 1; ++i) {
        const alea_hier_ray_path_entry_t* entry = &path->entries[i];
        if (placed_slab_owner(sys, entry->universe_id, world, &entry->transform, lo, radius,
                slab, first_frame, second_frame) != (int)entry->cell_index) return false;
    }
    const alea_hier_ray_path_entry_t* entry = &path->entries[leaf];
    int owner = placed_slab_owner(sys, entry->universe_id, world, &entry->transform,
                                  lo, radius, slab, first_frame, second_frame);
    if (owner < -1 || (owner >= 0 && alea_cell_entry_is_container(&sys->cells.data[owner]))) return false;
    if (owner < 0) path->count = leaf;
    else path->entries[leaf] = ordinary_entry(sys, owner, leaf, &entry->transform);
    return true;
}

static bool resolve_placed_slab(const alea_system_t* sys, const alea_ray_t* world,
    double mid, double radius, slab_t slab, const alea_matrix_t* first_frame,
    const alea_matrix_t* second_frame, int expected_parent, alea_hier_ray_path_t* path) {
    slab_t root_slab = slab;
    if (expected_parent >= 0) {
        /* Both plane cuts came from the current child. Verify the parent
         * independently, even if its inverse happens to match the child's. */
        root_slab.first = root_slab.second = UINT32_MAX;
    }
    alea_matrix_t identity;
    alea_matrix_identity(&identity);
    int root = placed_slab_owner(sys, 0, world, &identity, mid, radius, root_slab, first_frame, second_frame);
    if (root < 0 || (expected_parent >= 0 && root != expected_parent)) return false;
    const alea_cell_entry_t* parent = &sys->cells.data[root];
    if (parent->lat_type || parent->fill_universe <= 0) return false;
    const alea_matrix_t* fill = alea_cell_fill_matrix(sys, (uint32_t)root);
    if (!fill || !fill->has_inverse) return false;
    int child = placed_slab_owner(sys, parent->fill_universe, world, fill, mid, radius, slab,
                                  first_frame, second_frame);
    if (child < -1 || (child >= 0 && alea_cell_entry_is_container(&sys->cells.data[child]))) return false;
    path->entries[0] = ordinary_entry(sys, root, 0, &identity);
    path->count = 1;
    if (child >= 0) {
        path->entries[1] = ordinary_entry(sys, child, 1, fill);
        path->count = 2;
    }
    return true;
}

/* The ordinary plane solver may shift both reported endpoints by one
 * distance step. Establish the physical open interval before using either
 * plane's sign in a cell expression. Both planes must use the local ray. */
static bool adjacent_plane_slab(const alea_system_t* sys, const alea_ray_t* local,
    uint32_t first, uint32_t second, double lo, double hi, slab_t* slab) {
    const alea_plane_data_t* p = as_plane(sys, first);
    const alea_plane_data_t* q = as_plane(sys, second);
    if (!p || !q || first == second) return false;
    double enter, exit, width;
    if (!alea_substep_plane_cuts(local, p, q, lo, &enter, &exit, &width) ||
        !(enter < exit) || enter < nextafter(lo, -INFINITY) || enter > hi ||
        exit > nextafter(hi, INFINITY)) return false;
    double slope = p->a * local->dx + p->b * local->dy + p->c * local->dz;
    double other = q->a * local->dx + q->b * local->dy + q->c * local->dz;
    *slab = (slab_t){.first = first, .second = second, .first_sign = slope > 0 ? 1 : -1,
                     .second_sign = other > 0 ? -1 : 1};
    return true;
}

int alea_ray_adjacent_plane_cell_sign(const alea_system_t* sys, alea_node_id_t node,
    const alea_ray_t* local, uint32_t first, uint32_t second, double lo, double hi) {
    slab_t slab;
    if (!adjacent_plane_slab(sys, local, first, second, lo, hi, &slab)) return 0;
    slab.box = neighborhood(local, lo, 4 * (hi - lo));
    return cell_sign(sys, node, &slab, 0);
}

static bool resolve_adjacent_planes(const alea_system_t* sys, const alea_ray_t* world,
    uint32_t first, const alea_matrix_t* first_frame, uint32_t second, const alea_matrix_t* second_frame,
    double lo, double hi, alea_hier_ray_path_t* path) {
    if (!first_frame->has_inverse || !second_frame->has_inverse ||
        !same_inverse(first_frame, second_frame)) return false;
    alea_ray_t local;
    alea_ray_transform_inverse(first_frame, world, &local);
    slab_t slab;
    if (!adjacent_plane_slab(sys, &local, first, second, lo, hi, &slab)) return false;
    return resolve_placed_slab(sys, world, lo, 4 * (hi - lo), slab, first_frame, second_frame, -1, path);
}

static bool resolve_cylinder_interval(const alea_system_t* sys, const alea_ray_t* world,
    uint32_t first, const alea_matrix_t* first_frame, uint32_t second, const alea_matrix_t* second_frame,
    double lo, double hi, alea_hier_ray_path_t* path) {
    uint32_t ids[] = {first, second};
    const alea_matrix_t* frames[] = {first_frame, second_frame};
    double span = hi - lo, mid = lo + 0.5 * span, radius = 2 * span;
    /* Test on the exact rounding boundary between two neighboring ray
     * distances. Cuts strictly on opposite sides must round distinctly,
     * even when the interval has only one representable interior point. */
    double lower = mid < hi ? mid : nextafter(hi, lo);
    ray_wide_t witness = {lower, 0.5 * (nextafter(lower, hi) - lower)};
    if (!(witness.lo > 0) || !isfinite(mid + radius)) return false;

    int signs[2];
    for (int i = 0; i < 2; ++i) {
        guard_surface_t cylinder = {0};
        if (!load_cylinder(sys, ids[i], &cylinder)) return false;
        alea_ray_t local;
        alea_ray_transform_inverse(frames[i], world, &local);
        int before = cylinder_side(&cylinder, &local, (ray_wide_t){mid - radius, 0});
        int after = cylinder_side(&cylinder, &local, (ray_wide_t){mid + radius, 0});
        /* A quadratic with opposite signs at the ends has exactly one cut
         * inside. Verify that the witness lies between the two cuts. */
        signs[i] = i == 0 ? after : before;
        if (!before || after != -before ||
            cylinder_side(&cylinder, &local, witness) != signs[i]) return false;
    }

    slab_t slab = {.first = first, .second = second, .first_sign = signs[0], .second_sign = signs[1]};
    return resolve_placed_slab(sys, world, mid, radius, slab, first_frame, second_frame, -1, path);
}

static bool resolve_root_fill_interval(const alea_system_t* sys, const alea_ray_t* world,
    uint32_t first, const alea_matrix_t* first_frame, uint32_t second, const alea_matrix_t* second_frame,
    double lo, double hi, alea_hier_ray_path_t* path) {
    if (first >= sys->primitives.count || second >= sys->primitives.count ||
        !first_frame->has_inverse || !second_frame->has_inverse ||
        (first == second && same_inverse(first_frame, second_frame))) return false;
    alea_primitive_type_t type = sys->primitives.data[first].type;
    if (type >= ALEA_PRIMITIVE_CYLINDER_X && type <= ALEA_PRIMITIVE_CYLINDER_Z)
        return first != second && resolve_cylinder_interval(sys, world, first, first_frame, second, second_frame, lo, hi, path);

    const alea_plane_data_t* planes[] = {as_plane(sys, first), as_plane(sys, second)};
    if (!planes[0] || !planes[1]) return false;
    /* For distinct planes the existing proof covers a pair in the same
     * ordinary child placement. For a reused plane the two occurrences
     * themselves supply its constraints in the parent and child frames. */
    if (first != second && (path->count < 2 ||
        !same_inverse(first_frame, &path->entries[1].transform) ||
        !same_inverse(second_frame, &path->entries[1].transform))) return false;
    const alea_matrix_t* frames[] = {first_frame, second_frame};
    double distances[] = {lo, hi};
    int signs[2];
    for (int i = 0; i < 2; ++i) {
        alea_ray_t local;
        alea_ray_transform_inverse(frames[i], world, &local);
        const alea_plane_data_t* p = planes[i];
        double slope = p->a * local.dx + p->b * local.dy + p->c * local.dz;
        double cut;
        if (!rounded_cut(&local, p, distances[i], slope, &cut) || cut != distances[i]) return false;
        int after_crossing = slope > 0 ? 1 : -1;
        signs[i] = i == 0 ? after_crossing : -after_crossing;
    }
    /* Distinct rounded cuts establish their physical order. Enclose both
     * rounding intervals and use each plane sign only in its own placement. */
    double radius = (hi - lo) + 4 * (nextafter(hi, INFINITY) - hi);
    if (!isfinite(radius)) return false;
    slab_t slab = {.first = first, .second = second, .first_sign = signs[0], .second_sign = signs[1]};
    int expected_parent = first == second ? -1 : (int)path->entries[0].cell_index;
    return resolve_placed_slab(sys, world, lo, radius, slab, first_frame, second_frame, expected_parent, path);
}

/* A forward probe may have left the entire lattice placement. In that
 * case a point on the incident side can nominate its element. Ownership
 * still requires the plane ordering and all enclosing cells to verify. */
bool alea_ray_resolve_short_interval(alea_system_t* sys, const alea_ray_t* world,
    uint32_t first, const alea_matrix_t* first_frame, uint32_t second, const alea_matrix_t* second_frame,
    double lo, double hi, bool adjacent, bool left_lattice, alea_hier_ray_path_t* path) {
    if (resolve_lattice_planes(sys, world, first, first_frame, second, second_frame, lo, hi, path))
        return true;
    if (left_lattice && within_lattice_proof_window(lo, hi) &&
        as_plane(sys, first) && as_plane(sys, second)) {
        double radius = 64 * DBL_EPSILON * (1 + lo + fabs(world->ox) + fabs(world->oy) + fabs(world->oz));
        double x, y, z;
        alea_ray_point_at(world, fmax(0, lo - radius), &x, &y, &z);
        alea_hier_ray_path_t incident;
        alea_hier_cell_hit_t hit;
        if (alea_hier_spatial_find_path_at_point(sys, x, y, z, &hit, &incident) > 0 &&
            resolve_lattice_planes(sys, world, first, first_frame, second, second_frame, lo, hi, &incident)) {
            path->count = incident.count;
            memcpy(path->entries, incident.entries, (size_t)path->count * sizeof(*path->entries));
            return true;
        }
    }
    if (adjacent)
        return resolve_adjacent_planes(sys, world, first, first_frame, second, second_frame, lo, hi, path);
    return resolve_root_fill_interval(sys, world, first, first_frame, second, second_frame, lo, hi, path);
}

static bool empty_slab(const alea_system_t* sys, const alea_universe_t* u,
                       uint32_t incident, uint32_t outgoing, slab_t slab) {
    for (size_t i = 0; i < u->cell_indices.count; ++i) {
        const alea_cell_entry_t* cell = &sys->cells.data[u->cell_indices.data[i]];
        if (cell_sign(sys, cell->root_node_id, &slab, 0) != 1) return false;
    }
    slab.first_sign = -slab.first_sign;
    if (cell_sign(sys, sys->cells.data[incident].root_node_id, &slab, 0) != -1) return false;
    slab.first_sign = -slab.first_sign;
    slab.second_sign = -slab.second_sign;
    return cell_sign(sys, sys->cells.data[outgoing].root_node_id, &slab, 0) == -1;
}

int alea_substep_fill_transition(alea_system_t* sys, const alea_ray_t* world,
    double t, int incident, int surface, const alea_hier_ray_path_t* path,
    alea_ray_substep_fill_t* out) {
    if (incident < 0 || !path || path->count < 2) return 0;
    const alea_hier_ray_path_entry_t* leaf = &path->entries[path->count - 1];
    const alea_hier_ray_path_entry_t* parent = leaf - 1;
    const alea_cell_entry_t* before = &sys->cells.data[incident];
    const alea_cell_entry_t* after = &sys->cells.data[leaf->cell_index];
    if ((uint32_t)incident == leaf->cell_index || before->universe_id != after->universe_id ||
        parent->fill_universe != after->universe_id ||
        alea_cell_entry_is_container(before) || alea_cell_entry_is_container(after)) return 0;
    for (int i = 0; i < path->count; ++i) {
        if (path->entries[i].is_lattice) return 0;
    }

    alea_ray_t local;
    alea_ray_transform_inverse(&leaf->transform, world, &local);
    double step = nextafter(t, INFINITY) - t;
    for (size_t i = 0; i < before->surface_index_count; ++i) {
        const alea_surface_entry_t* a = &sys->surfaces.data[before->surface_indices[i]];
        if (surface > 0 && a->mc_surface_id != surface) continue;
        const alea_plane_data_t* p = as_plane(sys, a->primitive_id);
        if (!p) continue;
        double slope = p->a*local.dx + p->b*local.dy + p->c*local.dz;
        if (!isfinite(slope) || slope == 0) continue;
        for (size_t j = 0; j < after->surface_index_count; ++j) {
            const alea_surface_entry_t* b = &sys->surfaces.data[after->surface_indices[j]];
            const alea_plane_data_t* q = as_plane(sys, b->primitive_id);
            if (!q || a->primitive_id == b->primitive_id) continue;
            double enter, exit, width;
            if (!alea_substep_plane_cuts(&local, p, q, t, &enter, &exit, &width)) continue;
            double exit_slope = q->a*local.dx + q->b*local.dy + q->c*local.dz;
            slab_t slab = {neighborhood(&local, t, 8 * step), a->primitive_id, b->primitive_id,
                           slope > 0 ? 1 : -1, exit_slope > 0 ? -1 : 1, NULL};
            const alea_universe_t* u = alea_get_universe(sys, after->universe_id);
            if (!u || !empty_slab(sys, u, (uint32_t)incident, leaf->cell_index, slab)) continue;
            for (int k = 0; k < path->count - 1; ++k) {
                alea_ray_t frame;
                alea_ray_transform_inverse(&path->entries[k].transform, world, &frame);
                alea_bbox_t box = neighborhood(&frame, t, 8 * step);
                const alea_cell_entry_t* c = &sys->cells.data[path->entries[k].cell_index];
                slab_t ancestor = {.box=box, .first=UINT32_MAX, .second=UINT32_MAX};
                if (cell_sign(sys, c->root_node_id, &ancestor, 0) != -1) return 0;
            }
            *out = (alea_ray_substep_fill_t){parent->cell_id, after->universe_id,
                a->mc_surface_id, b->mc_surface_id, enter, exit, width};
            return 1;
        }
    }
    return 0;
}

int alea_substep_fill_at(alea_system_t* sys, const alea_ray_t* ray, double t,
                         alea_ray_substep_fill_t* out) {
    /* These probes only nominate siblings. The interval proof uses the
     * original ray and the original cell expressions in the placement. */
    double radius = 64 * DBL_EPSILON * (1 + t + fabs(ray->ox) + fabs(ray->oy) + fabs(ray->oz));
    alea_hier_ray_path_t paths[2];
    for (int k = 0; k < 2; ++k) {
        double x, y, z;
        alea_hier_cell_hit_t hit;
        alea_ray_point_at(ray, k ? t + radius : fmax(0, t - radius), &x, &y, &z);
        if (!isfinite(x) || !isfinite(y) || !isfinite(z)) return -1;
        int n = alea_hier_spatial_find_path_at_point(sys, x, y, z, &hit, &paths[k]);
        if (n <= 0) return n;
    }
    int n = paths[0].count;
    if (n != paths[1].count || n < 2) return 0;
    for (int i = 0; i < n - 1; ++i)
        if (paths[0].entries[i].cell_index != paths[1].entries[i].cell_index) return 0;
    return alea_substep_fill_transition(sys, ray, t, paths[0].entries[n - 1].cell_index,
                                        -1, &paths[1], out);
}

int alea_ray_diagnose_substep_fill(alea_system_t* sys, const double origin[3],
    const double direction[3], double t, alea_ray_substep_fill_t* out) {
    alea_ray_t ray;
    alea_clear_error_detail();
    if (!sys || !origin || !direction || !out || !(t > 0) || !isfinite(t)) goto invalid;
    for (int i = 0; i < 3; ++i) {
        if (!isfinite(origin[i]) || !isfinite(direction[i])) goto invalid;
    }
    if (alea_ray_init(&ray, origin[0], origin[1], origin[2], direction[0], direction[1], direction[2]) ||
        !(ray.dx || ray.dy || ray.dz)) goto invalid;
    if (alea_raycast_ensure_hier_caches(sys)) return -1;
    return alea_substep_fill_at(sys, &ray, t, out);
invalid:
    alea_set_error_detail(ALEA_ERR_INVALID_ARG,
                          "substep fill query requires a finite ray and a positive finite distance");
    return -1;
}
