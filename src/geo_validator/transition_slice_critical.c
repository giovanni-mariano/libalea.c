// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
// SPDX-License-Identifier: MPL-2.0

#include "transition_slice_critical.h"

#include "core/alea_occurrence.h"
#include "core/alea_eval.h"
#include "core/alea_spatial_hier.h"
#include "core/alea_system.h"
#include "core/alea_universe.h"
#include "primitives/bbox.h"
#include "primitives/primitive_desc.h"
#include "slice/curve_intersect.h"
#include "raycast/raycast.h"
#include "raycast/ray_epsilon.h"
#include "util/poly_solve.h"

#include <float.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

typedef struct {
    alea_curve_2d_t curve;
    uint64_t occurrence_key;
    uint64_t universe_occurrence_key;
    uint32_t surface_index;
    int cell_id;
    double bbox[4];
    double priority;
    uint8_t priority_class;
    uint8_t has_active_point;
    uint8_t has_parameter_domain;
    uint8_t has_scanline_domain;
    uint8_t scanline_root;
    uint8_t component_index;
    int8_t conic_branch;
    double active_uv[2];
    double parameter_min, parameter_max;
    double parameter_origin[2];
    double parameter_axis_x[2];
    double parameter_axis_y[2];
    double parameter_scale[2];
    double scanline_v_min, scanline_v_max;
    double scanline_endpoint_u[2];
} critical_curve_t;

typedef struct {
    double uv[2];
    size_t curve_index;
} critical_point_t;

typedef struct {
    int64_t qu, qv;
    uint8_t occupied;
} critical_point_slot_t;

typedef struct {
    size_t curve_index;
    double u_min, u_max, v_min, v_max;
} critical_curve_order_t;

typedef struct {
    int universe_id;
    uint64_t universe_occurrence_key;
} critical_universe_occurrence_t;

static void record_unsupported_curve(
    alea_curve_type_t type, alea_transition_slice_stats_t* stats);

static double conic_value(const alea_conic_2d_t* q, double u, double v);
static int curve_is_general_conic(alea_curve_type_t type);

static int canonicalize_open_conic(const alea_curve_2d_t* curve,
    critical_curve_t* item) {
    const alea_conic_2d_t* q = &curve->data.conic;
    /* Conic coefficients are homogeneous and their constant/linear terms can
     * grow quadratically under a large slice-coordinate translation.  Rank
     * tests must therefore use coefficients of the same degree; comparing a
     * quadratic eigenvalue with |F| incorrectly rejects valid distant
     * parabolas and hyperbolas. */
    const double quadratic_scale = fmax(
        fabs(q->A)+fabs(q->B)+fabs(q->C), DBL_MIN);
    const double linear_scale = fmax(fabs(q->D)+fabs(q->E), DBL_MIN);
    const double quadratic_tolerance = 1e-12*quadratic_scale;
    const double linear_tolerance = 1e-12*linear_scale;
    const double angle = 0.5*atan2(q->B, q->A-q->C);
    const double c = cos(angle), s = sin(angle);
    double eigenvalues[2] = {
        q->A*c*c + q->B*c*s + q->C*s*s,
        q->A*s*s - q->B*c*s + q->C*c*c};
    double axes[2][2] = {{c, s}, {-s, c}};
    const double linear[2] = {q->D*c+q->E*s, -q->D*s+q->E*c};

    memset(&item->parameter_origin, 0, sizeof(item->parameter_origin));
    item->has_parameter_domain = 0;
    item->conic_branch = 0;
    if (curve->type == ALEA_CURVE_PARABOLA) {
        const int quadratic = fabs(eigenvalues[1]) > fabs(eigenvalues[0]);
        const int axial = 1-quadratic;
        const double lambda = eigenvalues[quadratic];
        const double axial_linear = linear[axial];
        if (fabs(lambda) <= quadratic_tolerance ||
            fabs(axial_linear) <= linear_tolerance)
            return 0;
        const double x0 = -linear[quadratic]/(2.0*lambda);
        const long double completed = (long double)q->F -
            (long double)linear[quadratic]*linear[quadratic]/
                (4.0L*lambda);
        const double y0 = (double)(-completed/axial_linear);
        item->parameter_origin[0] =
            x0*axes[quadratic][0] + y0*axes[axial][0];
        item->parameter_origin[1] =
            x0*axes[quadratic][1] + y0*axes[axial][1];
        memcpy(item->parameter_axis_x, axes[quadratic], 2*sizeof(double));
        memcpy(item->parameter_axis_y, axes[axial], 2*sizeof(double));
        item->parameter_scale[0] = -lambda/axial_linear;
        item->parameter_scale[1] = 0.0;
    } else if (curve->type == ALEA_CURVE_HYPERBOLA) {
        const long double denominator =
            4.0L*q->A*q->C-(long double)q->B*q->B;
        if (fabsl(denominator) <=
            1e-12*quadratic_scale*quadratic_scale) return 0;
        const long double center_ld[2] = {
            ((long double)q->B*q->E-2.0L*q->C*q->D)/denominator,
            ((long double)q->B*q->D-2.0L*q->A*q->E)/denominator};
        const double center[2] = {
            (double)center_ld[0], (double)center_ld[1]};
        const long double constant_ld =
            (long double)q->A*center_ld[0]*center_ld[0] +
            (long double)q->B*center_ld[0]*center_ld[1] +
            (long double)q->C*center_ld[1]*center_ld[1] +
            (long double)q->D*center_ld[0] +
            (long double)q->E*center_ld[1] + q->F;
        const double constant = (double)constant_ld;
        int positive_axis = -1;
        for (int axis = 0; axis < 2; axis++) {
            if (fabs(eigenvalues[axis]) <= quadratic_tolerance) continue;
            const double squared_scale = -constant/eigenvalues[axis];
            if (squared_scale > 0.0 && isfinite(squared_scale))
                positive_axis = axis;
        }
        if (positive_axis < 0) return 0;
        const int negative_axis = 1-positive_axis;
        if (fabs(eigenvalues[negative_axis]) <= quadratic_tolerance) return 0;
        const double a2 = -constant/eigenvalues[positive_axis];
        const double b2 = constant/eigenvalues[negative_axis];
        if (!(a2 > 0.0) || !(b2 > 0.0) ||
            !isfinite(a2) || !isfinite(b2)) return 0;
        memcpy(item->parameter_origin, center, 2*sizeof(double));
        memcpy(item->parameter_axis_x, axes[positive_axis], 2*sizeof(double));
        memcpy(item->parameter_axis_y, axes[negative_axis], 2*sizeof(double));
        item->parameter_scale[0] = sqrt(a2);
        item->parameter_scale[1] = sqrt(b2);
    } else {
        return 0;
    }
    item->has_parameter_domain = 1;
    return 1;
}

static int open_conic_parameter(const critical_curve_t* item,
                                double u, double v,
                                double* parameter, int* branch) {
    if (!item->has_parameter_domain) return 0;
    const double du = u-item->parameter_origin[0];
    const double dv = v-item->parameter_origin[1];
    const double x = du*item->parameter_axis_x[0] +
                     dv*item->parameter_axis_x[1];
    const double y = du*item->parameter_axis_y[0] +
                     dv*item->parameter_axis_y[1];
    if (item->curve.type == ALEA_CURVE_PARABOLA) {
        *parameter = x;
        *branch = 0;
        return 1;
    }
    if (item->curve.type == ALEA_CURVE_HYPERBOLA) {
        if (!(item->parameter_scale[0] > 0.0) ||
            !(item->parameter_scale[1] > 0.0)) return 0;
        *parameter = asinh(y/item->parameter_scale[1]);
        *branch = x < 0.0 ? -1 : 1;
        return 1;
    }
    return 0;
}

static int open_conic_eval(const critical_curve_t* item, double parameter,
                           int branch, double* u, double* v) {
    if (!item->has_parameter_domain) return 0;
    double x, y;
    if (item->curve.type == ALEA_CURVE_PARABOLA) {
        x = parameter;
        y = item->parameter_scale[0]*parameter*parameter;
    } else if (item->curve.type == ALEA_CURVE_HYPERBOLA) {
        x = (branch < 0 ? -1.0 : 1.0)*item->parameter_scale[0]*cosh(parameter);
        y = item->parameter_scale[1]*sinh(parameter);
    } else {
        return 0;
    }
    *u = item->parameter_origin[0] + x*item->parameter_axis_x[0] +
         y*item->parameter_axis_y[0];
    *v = item->parameter_origin[1] + x*item->parameter_axis_x[1] +
         y*item->parameter_axis_y[1];
    return isfinite(*u) && isfinite(*v);
}

static void sort_small_doubles(double* values, int count) {
    for (int i = 1; i < count; i++) {
        const double value = values[i];
        int j = i;
        while (j > 0 && values[j-1] > value) {
            values[j] = values[j-1];
            j--;
        }
        values[j] = value;
    }
}

static int critical_curve_contains_point(const critical_curve_t* item,
                                         double u, double v,
                                         double tolerance) {
    if (item->has_scanline_domain) {
        if (v < item->scanline_v_min-tolerance ||
            v > item->scanline_v_max+tolerance) return 0;
        if (fabs(v-item->scanline_v_min) <= tolerance ||
            fabs(v-item->scanline_v_max) <= tolerance) {
            const int endpoint = fabs(v-item->scanline_v_max) <
                                 fabs(v-item->scanline_v_min);
            const double endpoint_tolerance = fmax(
                64.0*tolerance, 1e-6*fmax(1.0, fabs(u)));
            return fabs(u-item->scanline_endpoint_u[endpoint]) <=
                   endpoint_tolerance;
        }
        double roots[8];
        const int count = alea_curve_scanline_intersect(
            &item->curve, v, roots, (int)(sizeof(roots)/sizeof(roots[0])));
        if (count <= 0) return 0;
        sort_small_doubles(roots, count);
        int closest = 0;
        for (int i = 1; i < count; i++)
            if (fabs(roots[i]-u) < fabs(roots[closest]-u)) closest = i;
        return closest == item->scanline_root;
    }
    if (!item->has_parameter_domain) return 1;
    double parameter;
    int branch;
    if (!open_conic_parameter(item, u, v, &parameter, &branch)) return 0;
    if (item->curve.type == ALEA_CURVE_HYPERBOLA &&
        branch != item->conic_branch) return 0;
    return parameter >= item->parameter_min-tolerance &&
           parameter <= item->parameter_max+tolerance;
}

static int checked_add(size_t* total, size_t count, size_t item_size) {
    if (item_size && count > (SIZE_MAX - *total) / item_size) return -1;
    *total += count * item_size;
    return 0;
}

static size_t hash_capacity(size_t max_curves) {
    size_t wanted = max_curves > SIZE_MAX / 2 ? SIZE_MAX : max_curves * 2;
    size_t capacity = 8;
    while (capacity < wanted && capacity <= SIZE_MAX / 2) capacity <<= 1;
    return capacity < wanted ? 0 : capacity;
}

static int point_insert(critical_point_t* points, size_t* count,
                        size_t capacity, critical_point_slot_t* slots,
                        size_t slot_capacity, double tolerance,
                        double u, double v, size_t curve_index,
                        alea_transition_slice_stats_t* stats) {
    if (!isfinite(u) || !isfinite(v)) return 0;
    stats->critical_point_candidates++;
    const int64_t qu = (int64_t)llround(u / tolerance);
    const int64_t qv = (int64_t)llround(v / tolerance);
    size_t slot = (size_t)alea_occurrence_mix(
        (uint64_t)qu, (uint64_t)qv) & (slot_capacity - 1);
    for (size_t probe = 0; probe < slot_capacity; probe++) {
        critical_point_slot_t* entry = &slots[slot];
        if (!entry->occupied) {
            if (*count == capacity) return -1;
            entry->occupied = 1; entry->qu = qu; entry->qv = qv;
            points[*count].uv[0] = u; points[*count].uv[1] = v;
            points[*count].curve_index = curve_index;
            (*count)++;
            return 1;
        }
        if (entry->qu == qu && entry->qv == qv) {
            stats->critical_duplicate_points++;
            return 0;
        }
        slot = (slot + 1) & (slot_capacity - 1);
    }
    return -1;
}

static int point_in_tile(const alea_transition_slice_critical_tile_t* tile,
                         double u, double v, double tolerance) {
    return u >= tile->uv_min[0] - tolerance &&
           u <= tile->uv_max[0] + tolerance &&
           v >= tile->uv_min[1] - tolerance &&
           v <= tile->uv_max[1] + tolerance;
}

static int append_line_box_points(
    const alea_line_2d_t* line,
    const alea_transition_slice_critical_tile_t* tile,
    critical_point_t* points, size_t* point_count, size_t point_capacity,
    critical_point_slot_t* slots, size_t slot_capacity, double tolerance,
    size_t curve_index, alea_transition_slice_stats_t* stats) {
    const double bounds[4] = {tile->uv_min[0], tile->uv_max[0],
                              tile->uv_min[1], tile->uv_max[1]};
    for (int side = 0; side < 4; side++) {
        double u, v;
        if (side < 2) {
            if (fabs(line->b) <= tolerance) continue;
            u = bounds[side]; v = -(line->a * u + line->c) / line->b;
        } else {
            if (fabs(line->a) <= tolerance) continue;
            v = bounds[side]; u = -(line->b * v + line->c) / line->a;
        }
        if (!point_in_tile(tile, u, v, tolerance)) continue;
        if (point_insert(points, point_count, point_capacity, slots,
                         slot_capacity, tolerance, u, v, curve_index,
                         stats) < 0) return -1;
    }
    return 0;
}

static int generate_single_curve_points(
    const critical_curve_t* item, size_t curve_index,
    const alea_transition_slice_critical_tile_t* tile,
    critical_point_t* points, size_t* point_count, size_t point_capacity,
    critical_point_slot_t* slots, size_t slot_capacity, double tolerance,
    alea_transition_slice_stats_t* stats) {
    const alea_curve_2d_t* curve = &item->curve;
#define ADD_POINT(U, V) do { \
    const double pu_ = (U), pv_ = (V); \
    if (point_in_tile(tile, pu_, pv_, tolerance) && \
        point_insert(points, point_count, point_capacity, slots, \
                     slot_capacity, tolerance, pu_, pv_, curve_index, \
                     stats) < 0) return -1; \
} while (0)
    if (item->has_active_point)
        ADD_POINT(item->active_uv[0], item->active_uv[1]);
    if (item->has_scanline_domain) {
        ADD_POINT(item->scanline_endpoint_u[0], item->scanline_v_min);
        ADD_POINT(item->scanline_endpoint_u[1], item->scanline_v_max);
        return 0;
    }
    if (item->has_parameter_domain &&
        curve_is_general_conic(curve->type)) {
        double u, v;
        if (open_conic_eval(item, item->parameter_min,
                            item->conic_branch, &u, &v))
            ADD_POINT(u, v);
        if (open_conic_eval(item, item->parameter_max,
                            item->conic_branch, &u, &v))
            ADD_POINT(u, v);
        return 0;
    }
    switch (curve->type) {
    case ALEA_CURVE_POINT:
        ADD_POINT(curve->data.point[0], curve->data.point[1]);
        break;
    case ALEA_CURVE_LINE:
    case ALEA_CURVE_LINE_SEGMENT:
    case ALEA_CURVE_RAY:
        if (curve->type == ALEA_CURVE_LINE &&
            append_line_box_points(&curve->data.line, tile, points,
                point_count, point_capacity, slots, slot_capacity, tolerance,
                curve_index, stats) != 0) return -1;
        if (curve->type != ALEA_CURVE_LINE) {
            double u, v;
            if (alea_curve_eval(curve, curve->bounds.t_min, &u, &v))
                ADD_POINT(u, v);
            if (alea_curve_eval(curve, curve->bounds.t_max, &u, &v))
                ADD_POINT(u, v);
        }
        break;
    case ALEA_CURVE_CIRCLE: {
        const alea_circle_2d_t* c = &curve->data.circle;
        ADD_POINT(c->center[0] + c->radius, c->center[1]);
        ADD_POINT(c->center[0] - c->radius, c->center[1]);
        ADD_POINT(c->center[0], c->center[1] + c->radius);
        ADD_POINT(c->center[0], c->center[1] - c->radius);
        break;
    }
    case ALEA_CURVE_ARC: {
        double u, v;
        if (alea_curve_eval(curve, curve->bounds.theta_start, &u, &v))
            ADD_POINT(u, v);
        if (alea_curve_eval(curve, curve->bounds.theta_end, &u, &v))
            ADD_POINT(u, v);
        break;
    }
    case ALEA_CURVE_ELLIPSE: {
        const alea_ellipse_2d_t* e = &curve->data.ellipse;
        const double ct = cos(e->angle), st = sin(e->angle);
        const double du = hypot(e->semi_a * ct, e->semi_b * st);
        const double dv = hypot(e->semi_a * st, e->semi_b * ct);
        ADD_POINT(e->center[0] + du, e->center[1]);
        ADD_POINT(e->center[0] - du, e->center[1]);
        ADD_POINT(e->center[0], e->center[1] + dv);
        ADD_POINT(e->center[0], e->center[1] - dv);
        break;
    }
    case ALEA_CURVE_ELLIPSE_ARC: {
        double u, v;
        if (alea_curve_eval(curve, curve->bounds.theta_start, &u, &v))
            ADD_POINT(u, v);
        if (alea_curve_eval(curve, curve->bounds.theta_end, &u, &v))
            ADD_POINT(u, v);
        break;
    }
    case ALEA_CURVE_PARABOLA:
    case ALEA_CURVE_HYPERBOLA: {
        const alea_conic_2d_t* q = &curve->data.conic;
        const double fixed_v[2] = {tile->uv_min[1], tile->uv_max[1]};
        const double fixed_u[2] = {tile->uv_min[0], tile->uv_max[0]};
        double roots[2];
        for (int side = 0; side < 2; side++) {
            const double v = fixed_v[side];
            const int n = alea_solve_quadratic(
                q->A, q->B*v + q->D,
                q->C*v*v + q->E*v + q->F, roots);
            for (int i = 0; i < n; i++) ADD_POINT(roots[i], v);
        }
        for (int side = 0; side < 2; side++) {
            const double u = fixed_u[side];
            const int n = alea_solve_quadratic(
                q->C, q->B*u + q->E,
                q->A*u*u + q->D*u + q->F, roots);
            for (int i = 0; i < n; i++) ADD_POINT(u, roots[i]);
        }
        break;
    }
    case ALEA_CURVE_QUARTIC: {
        double roots[8];
        for (int side = 0; side < 2; side++) {
            const double v = side ? tile->uv_max[1] : tile->uv_min[1];
            const int n = alea_curve_scanline_intersect(
                curve, v, roots, (int)(sizeof(roots)/sizeof(roots[0])));
            for (int i = 0; i < n; i++) ADD_POINT(roots[i], v);
        }
        alea_curve_2d_t swapped = *curve;
        alea_torus_2d_t* torus = &swapped.data.torus;
        double temporary;
#define SWAP_VALUE(A, B) do { temporary = (A); (A) = (B); (B) = temporary; } while (0)
        SWAP_VALUE(torus->Ux, torus->Vx);
        SWAP_VALUE(torus->Uy, torus->Vy);
        SWAP_VALUE(torus->Uz, torus->Vz);
        SWAP_VALUE(torus->c1u, torus->c1v);
        SWAP_VALUE(torus->c2u, torus->c2v);
#undef SWAP_VALUE
        for (int side = 0; side < 2; side++) {
            const double u = side ? tile->uv_max[0] : tile->uv_min[0];
            const int n = alea_curve_scanline_intersect(
                &swapped, u, roots,
                (int)(sizeof(roots)/sizeof(roots[0])));
            for (int i = 0; i < n; i++) ADD_POINT(u, roots[i]);
        }
        break;
    }
    case ALEA_CURVE_POLYGON:
        for (int i = 0; i < curve->data.polygon.vertex_count; i++)
            ADD_POINT(curve->data.polygon.vertices[i][0],
                      curve->data.polygon.vertices[i][1]);
        break;
    default:
        record_unsupported_curve(curve->type, stats);
        return 1;
    }
#undef ADD_POINT
    return 0;
}

static int curve_normal_at(const alea_curve_2d_t* curve,
                           const double uv[2], double normal[2]) {
    switch (curve->type) {
    case ALEA_CURVE_LINE:
    case ALEA_CURVE_LINE_SEGMENT:
    case ALEA_CURVE_RAY:
        normal[0] = curve->data.line.a; normal[1] = curve->data.line.b;
        break;
    case ALEA_CURVE_CIRCLE:
    case ALEA_CURVE_ARC:
        normal[0] = uv[0] - curve->data.circle.center[0];
        normal[1] = uv[1] - curve->data.circle.center[1];
        break;
    case ALEA_CURVE_ELLIPSE:
    case ALEA_CURVE_ELLIPSE_ARC: {
        const alea_ellipse_2d_t* e = &curve->data.ellipse;
        const double ct = cos(e->angle), st = sin(e->angle);
        const double x = uv[0] - e->center[0], y = uv[1] - e->center[1];
        const double a = e->semi_a, b = e->semi_b;
        normal[0] = 2.0 * ((ct*ct/(a*a) + st*st/(b*b))*x +
                           ct*st*(1.0/(a*a)-1.0/(b*b))*y);
        normal[1] = 2.0 * (ct*st*(1.0/(a*a)-1.0/(b*b))*x +
                           (st*st/(a*a) + ct*ct/(b*b))*y);
        break;
    }
    case ALEA_CURVE_PARABOLA:
    case ALEA_CURVE_HYPERBOLA: {
        const alea_conic_2d_t* q = &curve->data.conic;
        normal[0] = 2.0*q->A*uv[0] + q->B*uv[1] + q->D;
        normal[1] = q->B*uv[0] + 2.0*q->C*uv[1] + q->E;
        break;
    }
    case ALEA_CURVE_QUARTIC: {
        const alea_torus_2d_t* t = &curve->data.torus;
        if (t->mode == 1 || t->mode == 2) {
            const double distances[2] = {
                fabs(hypot(uv[0]-t->c1u, uv[1]-t->c1v)-t->r1),
                fabs(hypot(uv[0]-t->c2u, uv[1]-t->c2v)-t->r2)};
            const int circle = distances[1] < distances[0] ? 1 : 0;
            normal[0] = uv[0] - (circle ? t->c2u : t->c1u);
            normal[1] = uv[1] - (circle ? t->c2v : t->c1v);
        } else {
            const double x = t->Ox + uv[0]*t->Ux + uv[1]*t->Vx;
            const double y = t->Oy + uv[0]*t->Uy + uv[1]*t->Vy;
            const double z = t->Oz + uv[0]*t->Uz + uv[1]*t->Vz;
            const double sigma = x*x + y*y + z*z +
                                 t->R*t->R - t->r*t->r;
            const double gx = 4.0*x*(sigma - 2.0*t->R*t->R);
            const double gy = 4.0*y*(sigma - 2.0*t->R*t->R);
            const double gz = 4.0*z*sigma;
            normal[0] = gx*t->Ux + gy*t->Uy + gz*t->Uz;
            normal[1] = gx*t->Vx + gy*t->Vy + gz*t->Vz;
        }
        break;
    }
    default: return 0;
    }
    const double length = hypot(normal[0], normal[1]);
    if (!(length > 0.0) || !isfinite(length)) return 0;
    normal[0] /= length; normal[1] /= length;
    return 1;
}

static int curve_order_compare(const void* lhs, const void* rhs) {
    const critical_curve_order_t* a = lhs;
    const critical_curve_order_t* b = rhs;
    if (a->u_min < b->u_min) return -1;
    if (a->u_min > b->u_min) return 1;
    return (a->curve_index > b->curve_index) -
           (a->curve_index < b->curve_index);
}

static int curve_priority_compare(const void* lhs, const void* rhs) {
    const critical_curve_t* a = lhs;
    const critical_curve_t* b = rhs;
    if (a->priority_class != b->priority_class)
        return a->priority_class > b->priority_class ? -1 : 1;
    if (a->priority < b->priority) return -1;
    if (a->priority > b->priority) return 1;
    return (a->surface_index > b->surface_index) -
           (a->surface_index < b->surface_index);
}

static int critical_curves_same_source(const critical_curve_t* first,
                                       const critical_curve_t* second) {
    return first->surface_index == second->surface_index &&
           first->universe_occurrence_key ==
               second->universe_occurrence_key;
}

static int curve_is_line(alea_curve_type_t type) {
    return type == ALEA_CURVE_LINE || type == ALEA_CURVE_LINE_SEGMENT ||
           type == ALEA_CURVE_RAY;
}

static int curve_is_circle(alea_curve_type_t type) {
    return type == ALEA_CURVE_CIRCLE || type == ALEA_CURVE_ARC;
}

static int curve_is_ellipse(alea_curve_type_t type) {
    return type == ALEA_CURVE_ELLIPSE || type == ALEA_CURVE_ELLIPSE_ARC;
}

static int curve_is_general_conic(alea_curve_type_t type) {
    return type == ALEA_CURVE_PARABOLA || type == ALEA_CURVE_HYPERBOLA;
}

static int curve_is_closed_conic(alea_curve_type_t type) {
    return curve_is_circle(type) || curve_is_ellipse(type);
}

static void record_unsupported_curve(
    alea_curve_type_t type, alea_transition_slice_stats_t* stats) {
    stats->critical_unsupported_curves++;
    if (type == ALEA_CURVE_PARABOLA)
        stats->critical_unsupported_parabola_curves++;
    else if (type == ALEA_CURVE_HYPERBOLA)
        stats->critical_unsupported_hyperbola_curves++;
    else if (type == ALEA_CURVE_QUARTIC)
        stats->critical_unsupported_quartic_curves++;
    else
        stats->critical_unsupported_other_curves++;
}

static void record_unsupported_pair(
    alea_curve_type_t first, alea_curve_type_t second,
    alea_transition_slice_stats_t* stats) {
    stats->critical_unsupported_curve_pairs++;
    if (first == ALEA_CURVE_QUARTIC || second == ALEA_CURVE_QUARTIC) {
        stats->critical_unsupported_quartic_pairs++;
        const alea_curve_type_t other = first == ALEA_CURVE_QUARTIC
            ? second : first;
        if (curve_is_line(other))
            stats->critical_unsupported_quartic_line_pairs++;
        else if (curve_is_closed_conic(other))
            stats->critical_unsupported_quartic_closed_conic_pairs++;
        else if (curve_is_general_conic(other))
            stats->critical_unsupported_quartic_general_conic_pairs++;
        else if (other == ALEA_CURVE_QUARTIC)
            stats->critical_unsupported_quartic_quartic_pairs++;
        else
            stats->critical_unsupported_quartic_other_pairs++;
    } else if (first == ALEA_CURVE_PARABOLA || second == ALEA_CURVE_PARABOLA ||
             first == ALEA_CURVE_HYPERBOLA || second == ALEA_CURVE_HYPERBOLA)
        stats->critical_unsupported_general_conic_pairs++;
    else if (first == ALEA_CURVE_POLYGON || second == ALEA_CURVE_POLYGON)
        stats->critical_unsupported_polygon_pairs++;
    else
        stats->critical_unsupported_other_pairs++;
}

static alea_ellipse_2d_t closed_conic_parameterization(
    const alea_curve_2d_t* curve) {
    if (curve_is_ellipse(curve->type)) return curve->data.ellipse;
    alea_ellipse_2d_t ellipse;
    ellipse.center[0] = curve->data.circle.center[0];
    ellipse.center[1] = curve->data.circle.center[1];
    ellipse.semi_a = curve->data.circle.radius;
    ellipse.semi_b = curve->data.circle.radius;
    ellipse.angle = 0.0;
    return ellipse;
}

static alea_conic_2d_t ellipse_implicit(const alea_ellipse_2d_t* ellipse) {
    const double ct = cos(ellipse->angle), st = sin(ellipse->angle);
    const double ia = 1.0 / (ellipse->semi_a * ellipse->semi_a);
    const double ib = 1.0 / (ellipse->semi_b * ellipse->semi_b);
    const double cx = ellipse->center[0], cy = ellipse->center[1];
    alea_conic_2d_t conic;
    conic.A = ct*ct*ia + st*st*ib;
    conic.B = 2.0*ct*st*(ia - ib);
    conic.C = st*st*ia + ct*ct*ib;
    conic.D = -2.0*conic.A*cx - conic.B*cy;
    conic.E = -conic.B*cx - 2.0*conic.C*cy;
    conic.F = conic.A*cx*cx + conic.B*cx*cy +
              conic.C*cy*cy - 1.0;
    return conic;
}

static double conic_value(const alea_conic_2d_t* conic,
                          double u, double v) {
    return conic->A*u*u + conic->B*u*v + conic->C*v*v +
           conic->D*u + conic->E*v + conic->F;
}

static double conic_value_scale(const alea_conic_2d_t* conic,
                                double u, double v) {
    return fabs(conic->A*u*u) + fabs(conic->B*u*v) +
           fabs(conic->C*v*v) + fabs(conic->D*u) +
           fabs(conic->E*v) + fabs(conic->F) + 1.0;
}

static int parametric_closed_conic_intersections(
    const alea_curve_2d_t* first, const alea_conic_2d_t* q,
    double uv[4][2]) {
    const alea_ellipse_2d_t parameter =
        closed_conic_parameterization(first);
    if (!(parameter.semi_a > 0.0) || !(parameter.semi_b > 0.0)) return -1;
    const double ct = cos(parameter.angle), st = sin(parameter.angle);
    const double cx = parameter.center[0], cy = parameter.center[1];
    const double p = parameter.semi_a*ct;
    const double r = parameter.semi_a*st;
    const double s = parameter.semi_b*ct;
    const double qsin = -parameter.semi_b*st;

    const double cc = q->A*p*p + q->B*p*r + q->C*r*r;
    const double cs = 2.0*q->A*p*qsin + q->B*(p*s + qsin*r) +
                      2.0*q->C*r*s;
    const double ss = q->A*qsin*qsin + q->B*qsin*s + q->C*s*s;
    const double lc = 2.0*q->A*cx*p + q->B*(cx*r + cy*p) +
                      2.0*q->C*cy*r + q->D*p + q->E*r;
    const double ls = 2.0*q->A*cx*qsin + q->B*(cx*s + cy*qsin) +
                      2.0*q->C*cy*s + q->D*qsin + q->E*s;
    const double constant = conic_value(q, cx, cy);
    double coefficients[5] = {
        cc + lc + constant,
        2.0*(cs + ls),
        -2.0*cc + 4.0*ss + 2.0*constant,
        2.0*(-cs + ls),
        cc - lc + constant
    };
    double scale = 0.0;
    for (int i = 0; i < 5; i++) scale = fmax(scale, fabs(coefficients[i]));
    if (!(scale > 1e-14) || !isfinite(scale)) return 0;
    for (int i = 0; i < 5; i++) coefficients[i] /= scale;

    double roots[4];
    const int root_count = alea_solve_quartic(
        coefficients[4], coefficients[3], coefficients[2],
        coefficients[1], coefficients[0], roots);
    int count = 0;
    for (int i = 0; i < root_count && count < 4; i++) {
        const double theta = 2.0*atan(roots[i]);
        if (!alea_curve_eval(first, theta, &uv[count][0], &uv[count][1]))
            return -1;
        int duplicate = 0;
        for (int j = 0; j < count; j++)
            if (hypot(uv[count][0] - uv[j][0],
                      uv[count][1] - uv[j][1]) <= 1e-10)
                duplicate = 1;
        if (!duplicate) count++;
    }

    if (count < 4 && fabs(coefficients[4]) <= 1e-11) {
        double u, v;
        if (!alea_curve_eval(first, 3.14159265358979323846, &u, &v))
            return -1;
        int duplicate = 0;
        for (int i = 0; i < count; i++)
            if (hypot(u - uv[i][0], v - uv[i][1]) <= 1e-10)
                duplicate = 1;
        if (!duplicate && fabs(conic_value(q, u, v)) <=
                1e-9*conic_value_scale(q, u, v)) {
            uv[count][0] = u; uv[count][1] = v; count++;
        }
    }
    return count;
}

static int closed_conic_intersections(const alea_curve_2d_t* first,
                                      const alea_curve_2d_t* second,
                                      double uv[4][2]) {
    const alea_ellipse_2d_t target =
        closed_conic_parameterization(second);
    if (!(target.semi_a > 0.0) || !(target.semi_b > 0.0)) return -1;
    const alea_conic_2d_t q = ellipse_implicit(&target);
    return parametric_closed_conic_intersections(first, &q, uv);
}

static double positive_angle(double angle) {
    const double turn = 6.28318530717958647693;
    angle = fmod(angle, turn);
    return angle < 0.0 ? angle + turn : angle;
}

static int angle_in_arc(const alea_curve_2d_t* curve, double angle,
                        double tolerance) {
    if (curve->type != ALEA_CURVE_ARC &&
        curve->type != ALEA_CURVE_ELLIPSE_ARC) return 1;
    const double start = curve->bounds.theta_start;
    const double end = curve->bounds.theta_end;
    const double turn = 6.28318530717958647693;
    while (angle < start - tolerance) angle += turn;
    while (angle > end + tolerance && angle - turn >= start - tolerance)
        angle -= turn;
    return angle >= start - tolerance && angle <= end + tolerance;
}

static double ellipse_parameter(const alea_ellipse_2d_t* ellipse,
                                double u, double v) {
    const double ct = cos(ellipse->angle), st = sin(ellipse->angle);
    const double x = u - ellipse->center[0], y = v - ellipse->center[1];
    return positive_angle(atan2(
        (-st*x + ct*y) / ellipse->semi_b,
        (ct*x + st*y) / ellipse->semi_a));
}

static int add_pair_point(
    const critical_curve_t* first_item,
    const critical_curve_t* second_item,
    double u, double v, size_t first_curve,
    const alea_transition_slice_critical_tile_t* tile,
    critical_point_t* points, size_t* point_count, size_t point_capacity,
    critical_point_slot_t* slots, size_t slot_capacity, double tolerance,
    alea_transition_slice_stats_t* stats) {
    stats->critical_pair_algebraic_points++;
    if (!point_in_tile(tile, u, v, tolerance)) return 0;
    if (!critical_curve_contains_point(first_item, u, v, tolerance) ||
        !critical_curve_contains_point(second_item, u, v, tolerance)) {
        stats->critical_pair_domain_rejections++;
        return 0;
    }
    const alea_curve_2d_t* first = &first_item->curve;
    const alea_curve_2d_t* second = &second_item->curve;
    const alea_curve_2d_t* domains[2] = {first, second};
    for (int i = 0; i < 2; i++) {
        const alea_curve_2d_t* curve = domains[i];
        if (curve->type != ALEA_CURVE_LINE_SEGMENT &&
            curve->type != ALEA_CURVE_RAY) {
            if (curve->type == ALEA_CURVE_ARC) {
                const double angle = atan2(
                    v - curve->data.circle.center[1],
                    u - curve->data.circle.center[0]);
                if (!angle_in_arc(curve, angle, tolerance)) {
                    stats->critical_pair_domain_rejections++;
                    return 0;
                }
            } else if (curve->type == ALEA_CURVE_ELLIPSE_ARC) {
                const double angle = ellipse_parameter(
                    &curve->data.ellipse, u, v);
                if (!angle_in_arc(curve, angle, tolerance)) {
                    stats->critical_pair_domain_rejections++;
                    return 0;
                }
            }
            continue;
        }
        const double dx = u - curve->data.line.point[0];
        const double dy = v - curve->data.line.point[1];
        const double t = dx*curve->data.line.direction[0] +
                         dy*curve->data.line.direction[1];
        if (t < curve->bounds.t_min - tolerance ||
            t > curve->bounds.t_max + tolerance) {
            stats->critical_pair_domain_rejections++;
            return 0;
        }
    }
    const int rc = point_insert(points, point_count, point_capacity, slots,
                                slot_capacity, tolerance, u, v, first_curve,
                                stats);
    if (rc > 0) stats->critical_pair_intersection_points++;
    return rc < 0 ? -1 : 0;
}

/* Intersect two circles/ellipses by parameterizing the first curve and
 * substituting it into the second curve's implicit quadratic.  The
 * tan(theta/2) substitution produces a real quartic and preserves tangent
 * roots, unlike angular sign-change sampling. */
static int solve_closed_conic_pair(
    const critical_curve_t* first_item,
    const critical_curve_t* second_item,
    size_t first_curve, const alea_transition_slice_critical_tile_t* tile,
    critical_point_t* points, size_t* point_count, size_t point_capacity,
    critical_point_slot_t* slots, size_t slot_capacity, double tolerance,
    alea_transition_slice_stats_t* stats) {
    const alea_curve_2d_t* first = &first_item->curve;
    const alea_curve_2d_t* second = &second_item->curve;
    double uv[4][2];
    const int count = closed_conic_intersections(first, second, uv);
    if (count < 0) return 1;
    for (int i = 0; i < count; i++) {
        if (add_pair_point(first_item, second_item,
                uv[i][0], uv[i][1], first_curve, tile,
                points, point_count, point_capacity, slots, slot_capacity,
                tolerance, stats) != 0) return -1;
    }
    return 0;
}

static void polynomial_add_product(double out[5], const double* first,
                                   int first_degree, const double* second,
                                   int second_degree, double factor) {
    for (int i = 0; i <= first_degree; i++)
        for (int j = 0; j <= second_degree && i + j <= 4; j++)
            out[i + j] += factor*first[i]*second[j];
}

static int solve_real_polynomial(double coefficients[5], double roots[4]) {
    double scale = 0.0;
    for (int i = 0; i < 5; i++) scale = fmax(scale, fabs(coefficients[i]));
    if (!(scale > 0.0) || !isfinite(scale)) return 0;
    for (int i = 0; i < 5; i++) coefficients[i] /= scale;
    int degree = 4;
    while (degree > 0 && fabs(coefficients[degree]) <= 1e-12) degree--;
    if (degree == 4)
        return alea_solve_quartic(
            coefficients[4], coefficients[3], coefficients[2],
            coefficients[1], coefficients[0], roots);
    if (degree == 3)
        return alea_solve_cubic(
            coefficients[3], coefficients[2], coefficients[1],
            coefficients[0], roots);
    if (degree == 2)
        return alea_solve_quadratic(
            coefficients[2], coefficients[1], coefficients[0], roots);
    if (degree == 1) {
        roots[0] = -coefficients[0]/coefficients[1];
        return 1;
    }
    return 0;
}

/* Resultant of two implicit quadratics, eliminating v.  Coefficients are
 * polynomials in u, stored low degree first. */
static int general_conic_intersections(const alea_conic_2d_t* first,
                                       const alea_conic_2d_t* second,
                                       double uv[4][2]) {
    const double a = first->C, d = second->C;
    const double b[2] = {first->E, first->B};
    const double c[3] = {first->F, first->D, first->A};
    const double e[2] = {second->E, second->B};
    const double f[3] = {second->F, second->D, second->A};
    const double coefficient_scale = fmax(1.0, fmax(
        fabs(first->A) + fabs(first->B) + fabs(first->C),
        fabs(second->A) + fabs(second->B) + fabs(second->C)));
    const double epsilon = 1e-14*coefficient_scale;
    double resultant[5] = {0};
    if (fabs(a) <= epsilon && fabs(d) <= epsilon) {
        /* (b*v+c, e*v+f): b*f-c*e. */
        polynomial_add_product(resultant, b, 1, f, 2, 1.0);
        polynomial_add_product(resultant, c, 2, e, 1, -1.0);
    } else if (fabs(a) <= epsilon || fabs(d) <= epsilon) {
        /* (b*v+c, d*v^2+e*v+f): d*c^2-e*b*c+f*b^2. */
        const double* linear_b = fabs(a) <= epsilon ? b : e;
        const double* linear_c = fabs(a) <= epsilon ? c : f;
        const double quadratic_a = fabs(a) <= epsilon ? d : a;
        const double* quadratic_b = fabs(a) <= epsilon ? e : b;
        const double* quadratic_c = fabs(a) <= epsilon ? f : c;
        double square[5] = {0}, product[5] = {0};
        polynomial_add_product(square, linear_c, 2, linear_c, 2, 1.0);
        for (int i = 0; i < 5; i++) resultant[i] += quadratic_a*square[i];
        polynomial_add_product(product, quadratic_b, 1, linear_b, 1, 1.0);
        double triple[5] = {0};
        polynomial_add_product(triple, product, 2, linear_c, 2, 1.0);
        for (int i = 0; i < 5; i++) resultant[i] -= triple[i];
        memset(square, 0, sizeof(square));
        polynomial_add_product(square, linear_b, 1, linear_b, 1, 1.0);
        polynomial_add_product(resultant, quadratic_c, 2, square, 2, 1.0);
    } else {
        /* (a*f-c*d)^2 - (a*e-b*d)*(b*f-c*e). */
        double p[3] = {0}, q[2] = {0}, r[5] = {0};
        for (int i = 0; i < 3; i++) p[i] = a*f[i] - d*c[i];
        for (int i = 0; i < 2; i++) q[i] = a*e[i] - d*b[i];
        polynomial_add_product(r, b, 1, f, 2, 1.0);
        polynomial_add_product(r, c, 2, e, 1, -1.0);
        polynomial_add_product(resultant, p, 2, p, 2, 1.0);
        polynomial_add_product(resultant, q, 1, r, 3, -1.0);
    }

    double u_roots[4];
    const int u_count = solve_real_polynomial(resultant, u_roots);
    int count = 0;
    for (int ui = 0; ui < u_count && count < 4; ui++) {
        const double u = u_roots[ui];
        const double first_a = first->C;
        const double first_b = first->B*u + first->E;
        const double first_c = first->A*u*u + first->D*u + first->F;
        const double second_a = second->C;
        const double second_b = second->B*u + second->E;
        const double second_c = second->A*u*u + second->D*u + second->F;
        const double first_strength = fabs(first_a) + fabs(first_b);
        const double second_strength = fabs(second_a) + fabs(second_b);
        double v_roots[2];
        const int v_count = first_strength >= second_strength
            ? alea_solve_quadratic(first_a, first_b, first_c, v_roots)
            : alea_solve_quadratic(second_a, second_b, second_c, v_roots);
        for (int vi = 0; vi < v_count && count < 4; vi++) {
            const double v = v_roots[vi];
            if (fabs(conic_value(first, u, v)) >
                    1e-8*conic_value_scale(first, u, v) ||
                fabs(conic_value(second, u, v)) >
                    1e-8*conic_value_scale(second, u, v))
                continue;
            int duplicate = 0;
            for (int i = 0; i < count; i++)
                if (hypot(u - uv[i][0], v - uv[i][1]) <= 1e-9)
                    duplicate = 1;
            if (!duplicate) { uv[count][0] = u; uv[count][1] = v; count++; }
        }
    }
    return count;
}

static int torus_line_intersections(const alea_torus_2d_t* torus,
                                    const alea_line_2d_t* line,
                                    double uv[4][2]) {
    if (torus->mode == 1 || torus->mode == 2) {
        const double centers[2][2] = {
            {torus->c1u, torus->c1v}, {torus->c2u, torus->c2v}};
        const double radii[2] = {torus->r1, torus->r2};
        int count = 0;
        for (int circle = 0; circle < 2 && count < 4; circle++) {
            if (!(radii[circle] > 0.0)) continue;
            const double x = line->point[0] - centers[circle][0];
            const double y = line->point[1] - centers[circle][1];
            const double b = x*line->direction[0] +
                             y*line->direction[1];
            double discriminant = b*b -
                (x*x + y*y - radii[circle]*radii[circle]);
            if (discriminant < -1e-12) continue;
            if (discriminant < 0.0) discriminant = 0.0;
            const double root = sqrt(discriminant);
            const double parameters[2] = {-b-root, -b+root};
            const int n = root > 1e-12 ? 2 : 1;
            for (int i = 0; i < n && count < 4; i++) {
                uv[count][0] = line->point[0] +
                    parameters[i]*line->direction[0];
                uv[count][1] = line->point[1] +
                    parameters[i]*line->direction[1];
                count++;
            }
        }
        return count;
    }

    const double u = line->point[0], v = line->point[1];
    const double du = line->direction[0], dv = line->direction[1];
    double px = torus->Ox + u*torus->Ux + v*torus->Vx;
    double py = torus->Oy + u*torus->Uy + v*torus->Vy;
    double pz = torus->Oz + u*torus->Uz + v*torus->Vz;
    const double dx = du*torus->Ux + dv*torus->Vx;
    const double dy = du*torus->Uy + dv*torus->Vy;
    const double dz = du*torus->Uz + dv*torus->Vz;
    const double direction2 = dx*dx + dy*dy + dz*dz;
    if (!(direction2 > 0.0)) return 0;
    const double shift = -(px*dx + py*dy + pz*dz)/direction2;
    px += shift*dx; py += shift*dy; pz += shift*dz;
    const double rho0 = px*px + py*py;
    const double rho1 = 2.0*(px*dx + py*dy);
    const double rho2 = dx*dx + dy*dy;
    const double sigma0 = rho0 + pz*pz;
    const double sigma1 = rho1 + 2.0*pz*dz;
    const double sigma2 = rho2 + dz*dz;
    const double s0 = sigma0 + torus->R*torus->R - torus->r*torus->r;
    const double s1 = sigma1, s2 = sigma2;
    const double four_r2 = 4.0*torus->R*torus->R;
    double roots[4];
    const int n = alea_solve_quartic(
        s2*s2, 2.0*s1*s2,
        s1*s1 + 2.0*s0*s2 - four_r2*rho2,
        2.0*s0*s1 - four_r2*rho1,
        s0*s0 - four_r2*rho0, roots);
    for (int i = 0; i < n; i++) {
        const double parameter = roots[i] + shift;
        uv[i][0] = u + parameter*du;
        uv[i][1] = v + parameter*dv;
    }
    return n;
}

static int torus_circle_components(const alea_torus_2d_t* torus,
                                   alea_curve_2d_t circles[2]) {
    if (torus->mode != 1 && torus->mode != 2) return 0;
    const double centers[2][2] = {
        {torus->c1u, torus->c1v}, {torus->c2u, torus->c2v}};
    const double radii[2] = {torus->r1, torus->r2};
    int count = 0;
    for (int i = 0; i < 2; i++) {
        if (!(radii[i] > 0.0)) continue;
        memset(&circles[count], 0, sizeof(circles[count]));
        circles[count].type = ALEA_CURVE_CIRCLE;
        circles[count].data.circle.center[0] = centers[i][0];
        circles[count].data.circle.center[1] = centers[i][1];
        circles[count].data.circle.radius = radii[i];
        count++;
    }
    return count;
}

#define CRITICAL_POLY_MAX_DEGREE 16

static double polynomial_eval_degree(const double* coefficients, int degree,
                                     double x) {
    double value = coefficients[degree];
    for (int i = degree-1; i >= 0; i--)
        value = value*x+coefficients[i];
    return value;
}

static double polynomial_eval_scale_degree(
    const double* coefficients, int degree, double x) {
    const double ax = fabs(x);
    double value = fabs(coefficients[degree]);
    for (int i = degree-1; i >= 0; i--)
        value = value*ax+fabs(coefficients[i]);
    return fmax(1.0, value);
}

static void insert_sorted_unique_root(double* roots, int* count,
                                      int capacity, double root) {
    if (!isfinite(root) || *count >= capacity) return;
    int position = 0;
    while (position < *count && roots[position] < root) position++;
    const double tolerance = 1e-9*fmax(1.0, fabs(root));
    if ((position > 0 && fabs(roots[position-1]-root) <= tolerance) ||
        (position < *count && fabs(roots[position]-root) <= tolerance))
        return;
    memmove(&roots[position+1], &roots[position],
            (size_t)(*count-position)*sizeof(double));
    roots[position] = root;
    (*count)++;
}

/* Isolate all distinct real roots of a degree <= 16 polynomial.  Derivative
 * roots partition the real line into monotonic intervals, so sign-changing
 * roots are bracketed and even-multiplicity roots are detected at critical
 * points.  Return -1 for an identically-zero polynomial. */
static int solve_real_polynomial_bounded(const double* input, int degree,
                                         double* roots, int capacity) {
    double coefficients[CRITICAL_POLY_MAX_DEGREE+1] = {0};
    double scale = 0.0;
    for (int i = 0; i <= degree; i++) scale = fmax(scale, fabs(input[i]));
    if (!(scale > 0.0) || !isfinite(scale)) return -1;
    for (int i = 0; i <= degree; i++) coefficients[i] = input[i]/scale;
    while (degree > 0 && fabs(coefficients[degree]) <= 1e-13) degree--;
    if (degree == 0) return 0;
    /* Do not delegate the recursive partition to the analytic cubic and
     * quartic solvers.  Their discriminant branches are intentionally
     * tolerance-based and can classify a nearly multiple derivative root
     * differently across libm implementations.  Recursing to the linear
     * base case keeps every partition bracketed by the same algorithm. */
    if (degree == 1) {
        roots[0] = -coefficients[0]/coefficients[1];
        return 1;
    }

    double derivative[CRITICAL_POLY_MAX_DEGREE] = {0};
    for (int i = 1; i <= degree; i++) derivative[i-1] = i*coefficients[i];
    double critical[CRITICAL_POLY_MAX_DEGREE] = {0};
    const int critical_count = solve_real_polynomial_bounded(
        derivative, degree-1, critical, CRITICAL_POLY_MAX_DEGREE);
    if (critical_count < 0) return -1;
    double bound = 1.0;
    for (int i = 0; i < degree; i++)
        bound = fmax(bound, 1.0+fabs(coefficients[i]/coefficients[degree]));
    if (!isfinite(bound)) return -1;

    double points[CRITICAL_POLY_MAX_DEGREE+1];
    int point_count = 0;
    points[point_count++] = -bound;
    for (int i = 0; i < critical_count; i++)
        if (critical[i] > -bound && critical[i] < bound)
            points[point_count++] = critical[i];
    points[point_count++] = bound;
    int count = 0;
    double previous = polynomial_eval_degree(coefficients, degree, points[0]);
    for (int i = 1; i < point_count; i++) {
        const double endpoint = points[i];
        double value = polynomial_eval_degree(coefficients, degree, endpoint);
        if ((previous < 0.0 && value > 0.0) ||
            (previous > 0.0 && value < 0.0)) {
            double lo = points[i-1], hi = endpoint, flo = previous;
            for (int iteration = 0; iteration < 100; iteration++) {
                const double mid = 0.5*(lo+hi);
                const double fm = polynomial_eval_degree(
                    coefficients, degree, mid);
                if (fabs(fm) <= 4.0e-15*
                    polynomial_eval_scale_degree(coefficients, degree, mid)) {
                    lo = mid; hi = mid;
                    break;
                }
                if ((flo < 0.0 && fm > 0.0) ||
                    (flo > 0.0 && fm < 0.0)) {
                    hi = mid;
                } else {
                    lo = mid; flo = fm;
                }
                if (hi-lo <= 8.0e-15*fmax(1.0, fabs(mid))) break;
            }
            insert_sorted_unique_root(
                roots, &count, capacity, 0.5*(lo+hi));
        }
        if (i < point_count-1 &&
            fabs(value) <= 1e-9*polynomial_eval_scale_degree(
                coefficients, degree, endpoint))
            insert_sorted_unique_root(roots, &count, capacity, endpoint);
        previous = value;
    }
    return count;
}

static void polynomial_multiply_degree8(
    const double* first, int first_degree,
    const double* second, int second_degree, double out[9]) {
    memset(out, 0, 9*sizeof(double));
    for (int i = 0; i <= first_degree; i++)
        for (int j = 0; j <= second_degree && i+j <= 8; j++)
            out[i+j] += first[i]*second[j];
}

/* Substitute rational u=U/D, v=V/D into a torus implicit equation and
 * multiply by D^4.  U, V, and D have degree <= 2, producing degree <= 8. */
static void torus_rational_polynomial(
    const alea_torus_2d_t* torus, const double U[3], const double V[3],
    const double D[3], double polynomial[9]) {
    double X[3], Y[3], Z[3];
    for (int i = 0; i < 3; i++) {
        X[i] = torus->Ox*D[i]+torus->Ux*U[i]+torus->Vx*V[i];
        Y[i] = torus->Oy*D[i]+torus->Uy*U[i]+torus->Vy*V[i];
        Z[i] = torus->Oz*D[i]+torus->Uz*U[i]+torus->Vz*V[i];
    }
    double x2[9], y2[9], z2[9], d2[9];
    polynomial_multiply_degree8(X, 2, X, 2, x2);
    polynomial_multiply_degree8(Y, 2, Y, 2, y2);
    polynomial_multiply_degree8(Z, 2, Z, 2, z2);
    polynomial_multiply_degree8(D, 2, D, 2, d2);
    double rho[9], sigma[9];
    for (int i = 0; i < 9; i++) {
        rho[i] = x2[i]+y2[i];
        sigma[i] = rho[i]+z2[i]+
            (torus->R*torus->R-torus->r*torus->r)*d2[i];
    }
    double sigma2[9], rho_d2[9];
    polynomial_multiply_degree8(sigma, 4, sigma, 4, sigma2);
    polynomial_multiply_degree8(rho, 4, d2, 4, rho_d2);
    for (int i = 0; i < 9; i++)
        polynomial[i] = sigma2[i]-4.0*torus->R*torus->R*rho_d2[i];
}

static double torus_implicit_value(const alea_torus_2d_t* torus,
                                   double u, double v) {
    const double x = torus->Ox+u*torus->Ux+v*torus->Vx;
    const double y = torus->Oy+u*torus->Uy+v*torus->Vy;
    const double z = torus->Oz+u*torus->Uz+v*torus->Vz;
    const double sigma = x*x+y*y+z*z+
                         torus->R*torus->R-torus->r*torus->r;
    return sigma*sigma-4.0*torus->R*torus->R*(x*x+y*y);
}

static double torus_implicit_scale(const alea_torus_2d_t* torus,
                                   double u, double v) {
    const double x = torus->Ox+u*torus->Ux+v*torus->Vx;
    const double y = torus->Oy+u*torus->Uy+v*torus->Vy;
    const double z = torus->Oz+u*torus->Uz+v*torus->Vz;
    const double sigma = fabs(x*x+y*y+z*z)+
                         torus->R*torus->R+torus->r*torus->r;
    return fmax(1.0, sigma*sigma+
        4.0*torus->R*torus->R*(x*x+y*y));
}

static void torus_implicit_gradient(const alea_torus_2d_t* torus,
                                    double u, double v,
                                    double* derivative_u,
                                    double* derivative_v) {
    const double x = torus->Ox+u*torus->Ux+v*torus->Vx;
    const double y = torus->Oy+u*torus->Uy+v*torus->Vy;
    const double z = torus->Oz+u*torus->Uz+v*torus->Vz;
    const double sigma = x*x+y*y+z*z+
                         torus->R*torus->R-torus->r*torus->r;
    const double gx = 4.0*x*(sigma-2.0*torus->R*torus->R);
    const double gy = 4.0*y*(sigma-2.0*torus->R*torus->R);
    const double gz = 4.0*z*sigma;
    *derivative_u = gx*torus->Ux+gy*torus->Uy+gz*torus->Uz;
    *derivative_v = gx*torus->Vx+gy*torus->Vy+gz*torus->Vz;
}

static void torus_bivariate_polynomial(const alea_torus_2d_t* torus,
                                       double polynomial[5][5]) {
    memset(polynomial, 0, 25*sizeof(double));
    const double affine[3][3] = {
        {torus->Ox, torus->Ux, torus->Vx},
        {torus->Oy, torus->Uy, torus->Vy},
        {torus->Oz, torus->Uz, torus->Vz}};
    const int u_degree[3] = {0, 1, 0};
    const int v_degree[3] = {0, 0, 1};
    double rho[3][3] = {{0}}, sigma[3][3] = {{0}};
    for (int coordinate = 0; coordinate < 3; coordinate++) {
        for (int first = 0; first < 3; first++) {
            for (int second = 0; second < 3; second++) {
                const int ui = u_degree[first]+u_degree[second];
                const int vi = v_degree[first]+v_degree[second];
                const double product = affine[coordinate][first]*
                                       affine[coordinate][second];
                sigma[ui][vi] += product;
                if (coordinate < 2) rho[ui][vi] += product;
            }
        }
    }
    sigma[0][0] += torus->R*torus->R-torus->r*torus->r;
    for (int ui = 0; ui <= 2; ui++) {
        for (int vi = 0; vi <= 2; vi++) {
            for (int uj = 0; uj <= 2; uj++) {
                for (int vj = 0; vj <= 2; vj++) {
                    if (ui+uj <= 4 && vi+vj <= 4)
                        polynomial[ui+uj][vi+vj] +=
                            sigma[ui][vi]*sigma[uj][vj];
                }
            }
            polynomial[ui][vi] -=
                4.0*torus->R*torus->R*rho[ui][vi];
        }
    }
}

/* Eliminate u between a quartic torus section and an implicit quadratic
 * conic.  The 6x6 Sylvester determinant has polynomial entries in v.  A
 * subset dynamic program evaluates that determinant without materializing
 * all 720 permutations. */
static void torus_conic_resultant(
    const double torus[5][5], const alea_conic_2d_t* conic,
    double resultant[17]) {
    double conic_polynomial[3][5] = {{0}};
    conic_polynomial[0][0] = conic->F;
    conic_polynomial[0][1] = conic->E;
    conic_polynomial[0][2] = conic->C;
    conic_polynomial[1][0] = conic->D;
    conic_polynomial[1][1] = conic->B;
    conic_polynomial[2][0] = conic->A;
    double matrix[6][6][5] = {{{0}}};
    for (int row = 0; row < 2; row++)
        for (int u_degree = 0; u_degree <= 4; u_degree++)
            memcpy(matrix[row][row+u_degree], torus[u_degree],
                   5*sizeof(double));
    for (int shift = 0; shift < 4; shift++)
        for (int u_degree = 0; u_degree <= 2; u_degree++)
            memcpy(matrix[2+shift][shift+u_degree],
                   conic_polynomial[u_degree], 5*sizeof(double));

    double determinants[64][17] = {{0}};
    determinants[0][0] = 1.0;
    for (int row = 0; row < 6; row++) {
        for (unsigned mask = 0; mask < 64; mask++) {
            int selected = 0;
            for (int bit = 0; bit < 6; bit++)
                if (mask & (1u << bit)) selected++;
            if (selected != row) continue;
            for (int column = 0; column < 6; column++) {
                if (mask & (1u << column)) continue;
                int inversions = 0;
                for (int prior = column+1; prior < 6; prior++)
                    if (mask & (1u << prior)) inversions++;
                const double sign = inversions & 1 ? -1.0 : 1.0;
                double* destination =
                    determinants[mask | (1u << column)];
                for (int first_degree = 0; first_degree <= 16;
                     first_degree++) {
                    if (determinants[mask][first_degree] == 0.0) continue;
                    for (int second_degree = 0; second_degree <= 4 &&
                         first_degree+second_degree <= 16; second_degree++)
                        destination[first_degree+second_degree] += sign*
                            determinants[mask][first_degree]*
                            matrix[row][column][second_degree];
                }
            }
        }
    }
    memcpy(resultant, determinants[63], 17*sizeof(double));
}

static int refine_torus_conic_point(const alea_torus_2d_t* torus,
                                    const alea_conic_2d_t* conic,
                                    double* u, double* v) {
    for (int iteration = 0; iteration < 12; iteration++) {
        const double f = torus_implicit_value(torus, *u, *v);
        const double g = conic_value(conic, *u, *v);
        double fu, fv;
        torus_implicit_gradient(torus, *u, *v, &fu, &fv);
        const double gu = 2.0*conic->A*(*u)+conic->B*(*v)+conic->D;
        const double gv = conic->B*(*u)+2.0*conic->C*(*v)+conic->E;
        const double determinant = fu*gv-fv*gu;
        const double jacobian_scale =
            fmax(1.0, (fabs(fu)+fabs(fv))*(fabs(gu)+fabs(gv)));
        if (fabs(determinant) <= 1e-14*jacobian_scale) break;
        const double du = (-f*gv+fv*g)/determinant;
        const double dv = (-fu*g+f*gu)/determinant;
        if (!isfinite(du) || !isfinite(dv)) break;
        *u += du; *v += dv;
        if (hypot(du, dv) <= 1e-13*fmax(1.0, hypot(*u, *v))) break;
    }
    return fabs(torus_implicit_value(torus, *u, *v)) <=
               1e-7*torus_implicit_scale(torus, *u, *v) &&
           fabs(conic_value(conic, *u, *v)) <=
               1e-8*conic_value_scale(conic, *u, *v);
}

static int torus_implicit_conic_intersections(
    const alea_torus_2d_t* torus, const alea_conic_2d_t* conic,
    double uv[8][2]) {
    double torus_polynomial[5][5], resultant[17], v_roots[16];
    torus_bivariate_polynomial(torus, torus_polynomial);
    torus_conic_resultant(torus_polynomial, conic, resultant);
    const int v_count = solve_real_polynomial_bounded(
        resultant, 16, v_roots, 16);
    if (v_count < 0) return -1;
    int count = 0;
    for (int vi = 0; vi < v_count && count < 8; vi++) {
        const double v = v_roots[vi];
        const double u_polynomial[3] = {
            conic->C*v*v+conic->E*v+conic->F,
            conic->B*v+conic->D, conic->A};
        double u_roots[2];
        const int u_count = solve_real_polynomial_bounded(
            u_polynomial, 2, u_roots, 2);
        if (u_count < 0) return -1;
        for (int ui = 0; ui < u_count && count < 8; ui++) {
            double u = u_roots[ui], refined_v = v;
            if (!refine_torus_conic_point(
                    torus, conic, &u, &refined_v)) continue;
            int duplicate = 0;
            for (int i = 0; i < count; i++)
                if (hypot(u-uv[i][0], refined_v-uv[i][1]) <=
                    1e-8*fmax(1.0, hypot(u, refined_v))) duplicate = 1;
            if (!duplicate) {
                uv[count][0] = u; uv[count][1] = refined_v; count++;
            }
        }
    }
    return count;
}

static void polynomial_add_product_degree16(
    double* out, const double* first, int first_degree,
    const double* second, int second_degree, double factor) {
    for (int i = 0; i <= first_degree; i++)
        for (int j = 0; j <= second_degree && i+j <= 16; j++)
            out[i+j] += factor*first[i]*second[j];
}

static void quartic_bezout_resultant(const double first[5][5],
                                     const double second[5][5],
                                     double resultant[17]) {
    double bezout[4][4][9] = {{{0}}};
    for (int high = 1; high <= 4; high++) {
        for (int low = 0; low < high; low++) {
            double cross[9] = {0};
            polynomial_add_product_degree16(
                cross, first[high], 4, second[low], 4, 1.0);
            polynomial_add_product_degree16(
                cross, first[low], 4, second[high], 4, -1.0);
            for (int offset = 0; offset < high-low; offset++) {
                const int row = high-1-offset;
                const int column = low+offset;
                for (int degree = 0; degree <= 8; degree++)
                    bezout[row][column][degree] += cross[degree];
            }
        }
    }
    static const unsigned char permutations[24][4] = {
        {0,1,2,3},{0,1,3,2},{0,2,1,3},{0,2,3,1},
        {0,3,1,2},{0,3,2,1},{1,0,2,3},{1,0,3,2},
        {1,2,0,3},{1,2,3,0},{1,3,0,2},{1,3,2,0},
        {2,0,1,3},{2,0,3,1},{2,1,0,3},{2,1,3,0},
        {2,3,0,1},{2,3,1,0},{3,0,1,2},{3,0,2,1},
        {3,1,0,2},{3,1,2,0},{3,2,0,1},{3,2,1,0}};
    memset(resultant, 0, 17*sizeof(double));
    for (int permutation = 0; permutation < 24; permutation++) {
        int inversions = 0;
        for (int i = 0; i < 4; i++)
            for (int j = i+1; j < 4; j++)
                if (permutations[permutation][i] >
                    permutations[permutation][j]) inversions++;
        double term[17] = {1.0};
        int term_degree = 0;
        for (int row = 0; row < 4; row++) {
            double product[17] = {0};
            polynomial_add_product_degree16(
                product, term, term_degree,
                bezout[row][permutations[permutation][row]], 8, 1.0);
            memcpy(term, product, sizeof(term));
            term_degree = term_degree+8 < 16 ? term_degree+8 : 16;
        }
        const double sign = inversions & 1 ? -1.0 : 1.0;
        for (int degree = 0; degree <= 16; degree++)
            resultant[degree] += sign*term[degree];
    }
}

static int refine_torus_pair_point(const alea_torus_2d_t* first,
                                   const alea_torus_2d_t* second,
                                   double* u, double* v) {
    for (int iteration = 0; iteration < 12; iteration++) {
        const double f = torus_implicit_value(first, *u, *v);
        const double g = torus_implicit_value(second, *u, *v);
        double fu, fv, gu, gv;
        torus_implicit_gradient(first, *u, *v, &fu, &fv);
        torus_implicit_gradient(second, *u, *v, &gu, &gv);
        const double determinant = fu*gv-fv*gu;
        const double jacobian_scale =
            fmax(1.0, (fabs(fu)+fabs(fv))*(fabs(gu)+fabs(gv)));
        if (fabs(determinant) <= 1e-14*jacobian_scale) break;
        const double du = (-f*gv+fv*g)/determinant;
        const double dv = (-fu*g+f*gu)/determinant;
        if (!isfinite(du) || !isfinite(dv)) break;
        *u += du; *v += dv;
        if (hypot(du, dv) <= 1e-13*fmax(1.0, hypot(*u, *v))) break;
    }
    return fabs(torus_implicit_value(first, *u, *v)) <=
               1e-7*torus_implicit_scale(first, *u, *v) &&
           fabs(torus_implicit_value(second, *u, *v)) <=
               1e-7*torus_implicit_scale(second, *u, *v);
}

static int torus_torus_intersections(const alea_torus_2d_t* first,
                                     const alea_torus_2d_t* second,
                                     double uv[16][2]) {
    /* Use a fixed orthogonal change of slice coordinates before elimination.
     * It prevents ordinary axis symmetries from projecting several distinct
     * intersections onto one high-multiplicity resultant root. */
    alea_torus_2d_t rotated_first = *first, rotated_second = *second;
    alea_torus_2d_t* rotated[2] = {&rotated_first, &rotated_second};
    const alea_torus_2d_t* original[2] = {first, second};
    const double rotation_c = 0.8, rotation_s = 0.6;
    for (int torus_index = 0; torus_index < 2; torus_index++) {
        rotated[torus_index]->Ux = rotation_c*original[torus_index]->Ux+
                                   rotation_s*original[torus_index]->Vx;
        rotated[torus_index]->Uy = rotation_c*original[torus_index]->Uy+
                                   rotation_s*original[torus_index]->Vy;
        rotated[torus_index]->Uz = rotation_c*original[torus_index]->Uz+
                                   rotation_s*original[torus_index]->Vz;
        rotated[torus_index]->Vx = -rotation_s*original[torus_index]->Ux+
                                    rotation_c*original[torus_index]->Vx;
        rotated[torus_index]->Vy = -rotation_s*original[torus_index]->Uy+
                                    rotation_c*original[torus_index]->Vy;
        rotated[torus_index]->Vz = -rotation_s*original[torus_index]->Uz+
                                    rotation_c*original[torus_index]->Vz;
    }
    double first_poly[5][5], second_poly[5][5], resultant[17];
    torus_bivariate_polynomial(&rotated_first, first_poly);
    torus_bivariate_polynomial(&rotated_second, second_poly);
    quartic_bezout_resultant(first_poly, second_poly, resultant);
    double v_roots[16];
    const int v_count = solve_real_polynomial_bounded(
        resultant, 16, v_roots, 16);
    if (v_count < 0) return -1;
    int count = 0;
    for (int vi = 0; vi < v_count && count < 16; vi++) {
        double first_u[5], second_u[5];
        double first_strength = 0.0, second_strength = 0.0;
        for (int ui = 0; ui <= 4; ui++) {
            first_u[ui] = polynomial_eval_degree(first_poly[ui], 4,
                                                  v_roots[vi]);
            second_u[ui] = polynomial_eval_degree(second_poly[ui], 4,
                                                   v_roots[vi]);
            first_strength += fabs(first_u[ui]);
            second_strength += fabs(second_u[ui]);
        }
        const double* source = first_strength >= second_strength
            ? first_u : second_u;
        double u_roots[4];
        const int u_count = solve_real_polynomial_bounded(
            source, 4, u_roots, 4);
        if (u_count < 0) continue;
        for (int ui = 0; ui < u_count && count < 16; ui++) {
            double rotated_u = u_roots[ui], rotated_v = v_roots[vi];
            if (!refine_torus_pair_point(&rotated_first, &rotated_second,
                                         &rotated_u, &rotated_v)) continue;
            const double u = rotation_c*rotated_u-rotation_s*rotated_v;
            const double v = rotation_s*rotated_u+rotation_c*rotated_v;
            if (fabs(torus_implicit_value(first, u, v)) >
                    1e-7*torus_implicit_scale(first, u, v) ||
                fabs(torus_implicit_value(second, u, v)) >
                    1e-7*torus_implicit_scale(second, u, v)) continue;
            int duplicate = 0;
            for (int i = 0; i < count; i++)
                if (hypot(u-uv[i][0], v-uv[i][1]) <=
                    1e-8*fmax(1.0, hypot(u, v))) duplicate = 1;
            if (!duplicate) { uv[count][0] = u; uv[count][1] = v; count++; }
        }
    }
    return count;
}

static void polynomial_add_discriminant_term(
    double out[17], double factor, const double coefficient[5][5],
    int power_a, int power_b, int power_c, int power_d, int power_e) {
    const int powers[5] = {power_e, power_d, power_c, power_b, power_a};
    double term[17] = {1.0};
    int degree = 0;
    for (int coefficient_index = 0; coefficient_index < 5;
         coefficient_index++) {
        for (int repeat = 0; repeat < powers[coefficient_index]; repeat++) {
            double product[17] = {0};
            polynomial_add_product_degree16(
                product, term, degree, coefficient[coefficient_index],
                4, 1.0);
            memcpy(term, product, sizeof(term));
            degree = degree+4 < 16 ? degree+4 : 16;
        }
    }
    for (int i = 0; i <= 16; i++) out[i] += factor*term[i];
}

static int torus_vertical_tangencies(const alea_torus_2d_t* torus,
                                     double values[16]) {
    double coefficient[5][5], discriminant[17] = {0};
    torus_bivariate_polynomial(torus, coefficient);
#define ADD_DISC(F, A, B, C, D, E) \
    polynomial_add_discriminant_term( \
        discriminant, (F), coefficient, (A), (B), (C), (D), (E))
    ADD_DISC( 256.0, 3,0,0,0,3);
    ADD_DISC(-192.0, 2,1,0,1,2);
    ADD_DISC(-128.0, 2,0,2,0,2);
    ADD_DISC( 144.0, 2,0,1,2,1);
    ADD_DISC( -27.0, 2,0,0,4,0);
    ADD_DISC( 144.0, 1,2,1,0,2);
    ADD_DISC(  -6.0, 1,2,0,2,1);
    ADD_DISC( -80.0, 1,1,2,1,1);
    ADD_DISC(  18.0, 1,1,1,3,0);
    ADD_DISC(  16.0, 1,0,4,0,1);
    ADD_DISC(  -4.0, 1,0,3,2,0);
    ADD_DISC( -27.0, 0,4,0,0,2);
    ADD_DISC(  18.0, 0,3,1,1,1);
    ADD_DISC(  -4.0, 0,3,0,3,0);
    ADD_DISC(  -4.0, 0,2,3,0,1);
    ADD_DISC(   1.0, 0,2,2,2,0);
#undef ADD_DISC
    return solve_real_polynomial_bounded(discriminant, 16, values, 16);
}

static int torus_closed_conic_intersections(
    const alea_torus_2d_t* torus, const alea_curve_2d_t* closed,
    double uv[8][2]) {
    const alea_ellipse_2d_t ellipse =
        closed_conic_parameterization(closed);
    if (!(ellipse.semi_a > 0.0) || !(ellipse.semi_b > 0.0)) return -1;
    const double c = cos(ellipse.angle), s = sin(ellipse.angle);
    const double pu = ellipse.semi_a*c;
    const double qu = -ellipse.semi_b*s;
    const double pv = ellipse.semi_a*s;
    const double qv = ellipse.semi_b*c;
    const double U[3] = {
        ellipse.center[0]+pu, 2.0*qu, ellipse.center[0]-pu};
    const double V[3] = {
        ellipse.center[1]+pv, 2.0*qv, ellipse.center[1]-pv};
    const double D[3] = {1.0, 0.0, 1.0};
    double polynomial[9], roots[8];
    torus_rational_polynomial(torus, U, V, D, polynomial);
    const int root_count = solve_real_polynomial_bounded(
        polynomial, 8, roots, 8);
    if (root_count < 0) return -1;
    int count = 0;
    for (int i = 0; i < root_count && count < 8; i++) {
        const double theta = 2.0*atan(roots[i]);
        double u, v;
        if (!alea_curve_eval(closed, theta, &u, &v)) return -1;
        if (fabs(torus_implicit_value(torus, u, v)) >
            1e-7*torus_implicit_scale(torus, u, v)) continue;
        int duplicate = 0;
        for (int j = 0; j < count; j++)
            if (hypot(u-uv[j][0], v-uv[j][1]) <= 1e-8)
                duplicate = 1;
        if (!duplicate) { uv[count][0] = u; uv[count][1] = v; count++; }
    }
    if (count < 8) {
        double u, v;
        if (!alea_curve_eval(closed, 3.14159265358979323846, &u, &v))
            return -1;
        if (fabs(torus_implicit_value(torus, u, v)) <=
            1e-7*torus_implicit_scale(torus, u, v)) {
            int duplicate = 0;
            for (int i = 0; i < count; i++)
                if (hypot(u-uv[i][0], v-uv[i][1]) <= 1e-8)
                    duplicate = 1;
            if (!duplicate) { uv[count][0] = u; uv[count][1] = v; count++; }
        }
    }
    return count;
}

static int torus_open_conic_intersections(
    const alea_torus_2d_t* torus, const critical_curve_t* open_item,
    double uv[8][2]) {
    critical_curve_t canonical = *open_item;
    if (!canonical.has_parameter_domain &&
        !canonicalize_open_conic(&canonical.curve, &canonical))
        return -1;
    const int branches[2] = {
        canonical.curve.type == ALEA_CURVE_HYPERBOLA
            ? (open_item->has_parameter_domain ? canonical.conic_branch : -1)
            : 0,
        1};
    const int branch_count = canonical.curve.type == ALEA_CURVE_HYPERBOLA &&
                             !open_item->has_parameter_domain ? 2 : 1;
    int count = 0;
    for (int bi = 0; bi < branch_count; bi++) {
        const int branch = branches[bi];
        double U[3] = {0}, V[3] = {0}, D[3] = {0};
        if (canonical.curve.type == ALEA_CURVE_PARABOLA) {
            D[0] = 1.0;
            U[0] = canonical.parameter_origin[0];
            U[1] = canonical.parameter_axis_x[0];
            U[2] = canonical.parameter_axis_y[0]*canonical.parameter_scale[0];
            V[0] = canonical.parameter_origin[1];
            V[1] = canonical.parameter_axis_x[1];
            V[2] = canonical.parameter_axis_y[1]*canonical.parameter_scale[0];
        } else if (canonical.curve.type == ALEA_CURVE_HYPERBOLA) {
            D[1] = 1.0;
            const double sign = branch < 0 ? -1.0 : 1.0;
            for (int coordinate = 0; coordinate < 2; coordinate++) {
                const double positive = 0.5*sign*canonical.parameter_scale[0]*
                    canonical.parameter_axis_x[coordinate];
                const double negative = 0.5*canonical.parameter_scale[1]*
                    canonical.parameter_axis_y[coordinate];
                double* numerator = coordinate ? V : U;
                numerator[0] = positive-negative;
                numerator[1] = canonical.parameter_origin[coordinate];
                numerator[2] = positive+negative;
            }
        } else {
            return -1;
        }
        double polynomial[9], roots[8];
        torus_rational_polynomial(torus, U, V, D, polynomial);
        const int root_count = solve_real_polynomial_bounded(
            polynomial, 8, roots, 8);
        if (root_count < 0) return -1;
        for (int i = 0; i < root_count && count < 8; i++) {
            double parameter = roots[i];
            if (canonical.curve.type == ALEA_CURVE_HYPERBOLA) {
                if (!(parameter > 0.0)) continue;
                parameter = log(parameter);
            }
            double u, v;
            if (!open_conic_eval(&canonical, parameter, branch, &u, &v))
                continue;
            if (fabs(torus_implicit_value(torus, u, v)) >
                    1e-7*torus_implicit_scale(torus, u, v) ||
                fabs(conic_value(&canonical.curve.data.conic, u, v)) >
                    1e-8*conic_value_scale(
                        &canonical.curve.data.conic, u, v))
                continue;
            int duplicate = 0;
            for (int j = 0; j < count; j++)
                if (hypot(u-uv[j][0], v-uv[j][1]) <= 1e-8)
                    duplicate = 1;
            if (!duplicate) { uv[count][0] = u; uv[count][1] = v; count++; }
        }
    }
    return count;
}

static int solve_curve_pair(
    const critical_curve_t* first_item,
    const critical_curve_t* second_item,
    size_t first_curve, const alea_transition_slice_critical_tile_t* tile,
    critical_point_t* points, size_t* point_count, size_t point_capacity,
    critical_point_slot_t* slots, size_t slot_capacity, double tolerance,
    alea_transition_slice_stats_t* stats) {
    const alea_curve_2d_t* first = &first_item->curve;
    const alea_curve_2d_t* second = &second_item->curve;
    const alea_curve_2d_t* line_curve = NULL;
    const alea_curve_2d_t* other = NULL;
    if (curve_is_line(first->type)) { line_curve = first; other = second; }
    else if (curve_is_line(second->type)) { line_curve = second; other = first; }

    if (curve_is_line(first->type) && curve_is_line(second->type)) {
        const alea_line_2d_t* a = &first->data.line;
        const alea_line_2d_t* b = &second->data.line;
        const double det = a->a * b->b - b->a * a->b;
        if (fabs(det) <= tolerance) return 0;
        const double u = (a->b * b->c - b->b * a->c) / det;
        const double v = (b->a * a->c - a->a * b->c) / det;
        return add_pair_point(first_item, second_item, u, v, first_curve, tile,
                              points, point_count,
                              point_capacity, slots, slot_capacity, tolerance,
                              stats);
    }
    if (line_curve && curve_is_circle(other->type)) {
        const alea_line_2d_t* l = &line_curve->data.line;
        const alea_circle_2d_t* c = &other->data.circle;
        const double px = l->point[0] - c->center[0];
        const double py = l->point[1] - c->center[1];
        const double dx = l->direction[0], dy = l->direction[1];
        const double b = 2.0 * (px*dx + py*dy);
        const double cc = px*px + py*py - c->radius*c->radius;
        double disc = b*b - 4.0*cc;
        if (disc < -tolerance) return 0;
        if (disc < 0.0) disc = 0.0;
        const double root = sqrt(disc);
        const double ts[2] = {(-b-root)/2.0, (-b+root)/2.0};
        const int count = root <= tolerance ? 1 : 2;
        for (int i = 0; i < count; i++)
            if (add_pair_point(first_item, second_item,
                    l->point[0] + ts[i]*dx,
                    l->point[1] + ts[i]*dy, first_curve, tile, points,
                    point_count, point_capacity, slots, slot_capacity,
                    tolerance, stats) != 0) return -1;
        return 0;
    }
    if (curve_is_circle(first->type) && curve_is_circle(second->type)) {
        const alea_circle_2d_t* a = &first->data.circle;
        const alea_circle_2d_t* b = &second->data.circle;
        const double dx = b->center[0]-a->center[0];
        const double dy = b->center[1]-a->center[1];
        const double d = hypot(dx, dy);
        if (d <= tolerance || d > a->radius+b->radius+tolerance ||
            d < fabs(a->radius-b->radius)-tolerance) return 0;
        const double x = (a->radius*a->radius - b->radius*b->radius + d*d) /
                         (2.0*d);
        double h2 = a->radius*a->radius - x*x;
        if (h2 < -tolerance) return 0;
        if (h2 < 0.0) h2 = 0.0;
        const double h = sqrt(h2), ux = dx/d, uy = dy/d;
        const double pu = a->center[0] + x*ux;
        const double pv = a->center[1] + x*uy;
        if (add_pair_point(first_item, second_item, pu-h*uy, pv+h*ux,
                first_curve, tile, points,
                point_count, point_capacity, slots, slot_capacity, tolerance,
                stats) != 0) return -1;
        if (h > tolerance && add_pair_point(first_item, second_item,
                pu+h*uy, pv-h*ux, first_curve, tile, points, point_count,
                point_capacity, slots, slot_capacity, tolerance, stats) != 0)
            return -1;
        return 0;
    }
    if (line_curve && curve_is_ellipse(other->type)) {
        const alea_line_2d_t* l = &line_curve->data.line;
        const alea_ellipse_2d_t* e = &other->data.ellipse;
        const double ct = cos(e->angle), st = sin(e->angle);
        const double x = l->point[0]-e->center[0];
        const double y = l->point[1]-e->center[1];
        const double px = ct*x + st*y, py = -st*x + ct*y;
        const double dx = ct*l->direction[0] + st*l->direction[1];
        const double dy = -st*l->direction[0] + ct*l->direction[1];
        const double aa = dx*dx/(e->semi_a*e->semi_a) +
                          dy*dy/(e->semi_b*e->semi_b);
        const double bb = 2.0*(px*dx/(e->semi_a*e->semi_a) +
                               py*dy/(e->semi_b*e->semi_b));
        const double cc = px*px/(e->semi_a*e->semi_a) +
                          py*py/(e->semi_b*e->semi_b) - 1.0;
        double disc = bb*bb - 4.0*aa*cc;
        if (!(aa > tolerance) || disc < -tolerance) return 0;
        if (disc < 0.0) disc = 0.0;
        const double root = sqrt(disc);
        const double ts[2] = {(-bb-root)/(2.0*aa), (-bb+root)/(2.0*aa)};
        const int count = root <= tolerance ? 1 : 2;
        for (int i = 0; i < count; i++)
            if (add_pair_point(first_item, second_item,
                    l->point[0] + ts[i]*l->direction[0],
                    l->point[1] + ts[i]*l->direction[1], first_curve, tile,
                    points, point_count, point_capacity, slots, slot_capacity,
                    tolerance, stats) != 0) return -1;
        return 0;
    }
    if (line_curve && curve_is_general_conic(other->type)) {
        const alea_line_2d_t* l = &line_curve->data.line;
        const alea_conic_2d_t* q = &other->data.conic;
        const double u = l->point[0], v = l->point[1];
        const double du = l->direction[0], dv = l->direction[1];
        const double aa = q->A*du*du + q->B*du*dv + q->C*dv*dv;
        const double bb = 2.0*q->A*u*du + q->B*(u*dv + v*du) +
                          2.0*q->C*v*dv + q->D*du + q->E*dv;
        const double cc = conic_value(q, u, v);
        double roots[2];
        const int n = alea_solve_quadratic(aa, bb, cc, roots);
        for (int i = 0; i < n; i++)
            if (add_pair_point(first_item, second_item,
                    u + roots[i]*du, v + roots[i]*dv,
                    first_curve, tile, points, point_count, point_capacity,
                    slots, slot_capacity, tolerance, stats) != 0)
                return -1;
        return 0;
    }
    if (line_curve && other->type == ALEA_CURVE_QUARTIC) {
        double uv[4][2];
        const int n = torus_line_intersections(
            &other->data.torus, &line_curve->data.line, uv);
        for (int i = 0; i < n; i++)
            if (add_pair_point(first_item, second_item, uv[i][0], uv[i][1],
                    first_curve, tile, points, point_count, point_capacity,
                    slots, slot_capacity, tolerance, stats) != 0)
                return -1;
        return 0;
    }
    if (first->type == ALEA_CURVE_QUARTIC ||
        second->type == ALEA_CURVE_QUARTIC) {
        if (first->type == ALEA_CURVE_QUARTIC &&
            second->type == ALEA_CURVE_QUARTIC) {
            alea_curve_2d_t first_circles[2], second_circles[2];
            const int first_count = torus_circle_components(
                &first->data.torus, first_circles);
            const int second_count = torus_circle_components(
                &second->data.torus, second_circles);
            if (first_count == 0 && second_count == 0) {
                double uv[16][2];
                const int n = torus_torus_intersections(
                    &first->data.torus, &second->data.torus, uv);
                if (n < 0) return 1;
                for (int i = 0; i < n; i++)
                    if (add_pair_point(first_item, second_item,
                            uv[i][0], uv[i][1], first_curve, tile,
                            points, point_count, point_capacity, slots,
                            slot_capacity, tolerance, stats) != 0)
                        return -1;
                return 0;
            }
            if ((first_count == 0) != (second_count == 0)) {
                const alea_torus_2d_t* general = first_count == 0
                    ? &first->data.torus : &second->data.torus;
                const alea_curve_2d_t* circles = first_count == 0
                    ? second_circles : first_circles;
                const int circle_count = first_count == 0
                    ? second_count : first_count;
                for (int circle = 0; circle < circle_count; circle++) {
                    double uv[8][2];
                    const int n = torus_closed_conic_intersections(
                        general, &circles[circle], uv);
                    if (n < 0) return 1;
                    for (int i = 0; i < n; i++)
                        if (add_pair_point(first_item, second_item,
                                uv[i][0], uv[i][1], first_curve, tile,
                                points, point_count, point_capacity, slots,
                                slot_capacity, tolerance, stats) != 0)
                            return -1;
                }
                return 0;
            }
        }
        const alea_curve_2d_t* torus_curve =
            first->type == ALEA_CURVE_QUARTIC ? first : second;
        const alea_curve_2d_t* companion =
            torus_curve == first ? second : first;
        alea_curve_2d_t torus_circles[2];
        const int torus_count = torus_circle_components(
            &torus_curve->data.torus, torus_circles);
        if (torus_count == 0 && curve_is_closed_conic(companion->type)) {
            double uv[8][2];
            const int n = torus_closed_conic_intersections(
                &torus_curve->data.torus, companion, uv);
            if (n < 0) return 1;
            for (int i = 0; i < n; i++)
                if (add_pair_point(first_item, second_item,
                        uv[i][0], uv[i][1], first_curve, tile,
                        points, point_count, point_capacity, slots,
                        slot_capacity, tolerance, stats) != 0)
                    return -1;
            return 0;
        }
        if (torus_count == 0 && curve_is_general_conic(companion->type)) {
            const critical_curve_t* open_item =
                torus_curve == first ? second_item : first_item;
            double uv[8][2];
            const int n = torus_open_conic_intersections(
                &torus_curve->data.torus, open_item, uv);
            if (n < 0) return 1;
            for (int i = 0; i < n; i++)
                if (add_pair_point(first_item, second_item,
                        uv[i][0], uv[i][1], first_curve, tile,
                        points, point_count, point_capacity, slots,
                        slot_capacity, tolerance, stats) != 0)
                    return -1;
            return 0;
        }
        if (torus_count > 0 && curve_is_closed_conic(companion->type)) {
            for (int circle = 0; circle < torus_count; circle++) {
                double uv[4][2];
                const int n = closed_conic_intersections(
                    &torus_circles[circle], companion, uv);
                if (n < 0) return 1;
                for (int i = 0; i < n; i++)
                    if (add_pair_point(first_item, second_item,
                            uv[i][0], uv[i][1],
                            first_curve, tile, points, point_count,
                            point_capacity, slots, slot_capacity, tolerance,
                            stats) != 0)
                        return -1;
            }
            return 0;
        }
        if (torus_count > 0 &&
            curve_is_general_conic(companion->type)) {
            for (int circle = 0; circle < torus_count; circle++) {
                double uv[4][2];
                const int n = parametric_closed_conic_intersections(
                    &torus_circles[circle], &companion->data.conic, uv);
                if (n < 0) return 1;
                for (int i = 0; i < n; i++)
                    if (add_pair_point(first_item, second_item,
                            uv[i][0], uv[i][1],
                            first_curve, tile, points, point_count,
                            point_capacity, slots, slot_capacity, tolerance,
                            stats) != 0)
                        return -1;
            }
            return 0;
        }
        if (torus_count > 0 && companion->type == ALEA_CURVE_QUARTIC) {
            alea_curve_2d_t companion_circles[2];
            const int companion_count = torus_circle_components(
                &companion->data.torus, companion_circles);
            if (companion_count == 0) return 1;
            for (int a = 0; a < torus_count; a++) {
                for (int b = 0; b < companion_count; b++) {
                    double uv[4][2];
                    const int n = closed_conic_intersections(
                        &torus_circles[a], &companion_circles[b], uv);
                    if (n < 0) return 1;
                    for (int i = 0; i < n; i++)
                        if (add_pair_point(first_item, second_item,
                                uv[i][0], uv[i][1], first_curve, tile,
                                points, point_count, point_capacity, slots,
                                slot_capacity, tolerance, stats) != 0)
                            return -1;
                }
            }
            return 0;
        }
    }
    if ((curve_is_closed_conic(first->type) &&
         curve_is_general_conic(second->type)) ||
        (curve_is_closed_conic(second->type) &&
         curve_is_general_conic(first->type))) {
        const alea_curve_2d_t* closed = curve_is_closed_conic(first->type)
            ? first : second;
        const alea_curve_2d_t* conic = closed == first ? second : first;
        double uv[4][2];
        const int n = parametric_closed_conic_intersections(
            closed, &conic->data.conic, uv);
        if (n < 0) return 1;
        for (int i = 0; i < n; i++)
            if (add_pair_point(first_item, second_item, uv[i][0], uv[i][1],
                    first_curve, tile, points, point_count, point_capacity,
                    slots, slot_capacity, tolerance, stats) != 0)
                return -1;
        return 0;
    }
    if (curve_is_general_conic(first->type) &&
        curve_is_general_conic(second->type)) {
        double uv[4][2];
        const int n = general_conic_intersections(
            &first->data.conic, &second->data.conic, uv);
        for (int i = 0; i < n; i++)
            if (add_pair_point(first_item, second_item, uv[i][0], uv[i][1],
                    first_curve, tile, points, point_count, point_capacity,
                    slots, slot_capacity, tolerance, stats) != 0)
                return -1;
        return 0;
    }
    if (curve_is_closed_conic(first->type) &&
        curve_is_closed_conic(second->type))
        return solve_closed_conic_pair(
            first_item, second_item, first_curve, tile, points, point_count,
            point_capacity, slots, slot_capacity, tolerance, stats);
    return 1;
}

static void chain_occurrences(const alea_hier_spatial_chain_hit_t* hit,
                              uint64_t* terminal_key,
                              uint64_t* universe_key) {
    alea_occurrence_state_t state;
    alea_occurrence_state_init(&state);
    for (uint8_t i = 0; i < hit->ancestor_count; i++) {
        alea_occurrence_state_step(
            &state, hit->ancestor_cell_indices[i],
            hit->ancestor_is_lattice[i], hit->ancestor_lattice_i[i],
            hit->ancestor_lattice_j[i], hit->ancestor_lattice_k[i]);
    }
    *universe_key = state.universe_key;
    *terminal_key = alea_occurrence_state_step(
        &state, hit->hit.cell_index, 0, 0, 0, 0);
}

static void tile_world_point(const alea_slice_view_t* view,
                             double u, double v, double world[3]) {
    for (int axis = 0; axis < 3; axis++)
        world[axis] = view->plane.origin[axis] +
            u * view->plane.u_axis[axis] + v * view->plane.v_axis[axis];
}

static alea_bbox_t tile_local_bbox(
    const alea_slice_view_t* view,
    const alea_transition_slice_critical_tile_t* tile,
    const alea_matrix_t* local_to_world) {
    alea_bbox_t bbox = alea_bbox_empty();
    for (int iu = 0; iu < 2; iu++) {
        const double u = iu ? tile->uv_max[0] : tile->uv_min[0];
        for (int iv = 0; iv < 2; iv++) {
            const double v = iv ? tile->uv_max[1] : tile->uv_min[1];
            double local[3];
            tile_world_point(view, u, v, local);
            alea_matrix_transform_point_inverse(
                local_to_world, &local[0], &local[1], &local[2]);
            if (local[0] < bbox.min_x) bbox.min_x = local[0];
            if (local[0] > bbox.max_x) bbox.max_x = local[0];
            if (local[1] < bbox.min_y) bbox.min_y = local[1];
            if (local[1] > bbox.max_y) bbox.max_y = local[1];
            if (local[2] < bbox.min_z) bbox.min_z = local[2];
            if (local[2] > bbox.max_z) bbox.max_z = local[2];
        }
    }
    const double scale = fmax(1.0, fmax(
        fabs(tile->uv_max[0] - tile->uv_min[0]),
        fabs(tile->uv_max[1] - tile->uv_min[1])));
    const double epsilon = 1e-9 * scale;
    bbox.min_x -= epsilon; bbox.max_x += epsilon;
    bbox.min_y -= epsilon; bbox.max_y += epsilon;
    bbox.min_z -= epsilon; bbox.max_z += epsilon;
    return bbox;
}

static uint64_t path_universe_occurrence_key(
    const alea_hier_ray_path_t* path, int level) {
    alea_occurrence_state_t state;
    alea_occurrence_state_init(&state);
    for (int i = 0; i < level; i++) {
        const alea_hier_ray_path_entry_t* entry = &path->entries[i];
        alea_occurrence_state_step(
            &state, entry->cell_index, entry->is_lattice,
            entry->lat_i, entry->lat_j, entry->lat_k);
    }
    return state.universe_key;
}

static void occurrence_chain_hit(
    const alea_hier_ray_path_t* path, int level,
    const alea_spatial_hit_t* candidate,
    alea_hier_spatial_chain_hit_t* out) {
    memset(out, 0, sizeof(*out));
    out->hit = *candidate;
    out->hit.depth = path->entries[level].depth;
    out->hit.universe_id = path->entries[level].universe_id;
    out->hit.transform = path->entries[level].transform;
    int ancestor_count = level;
    if (ancestor_count > ALEA_HIER_SPATIAL_HIT_CHAIN_MAX) {
        ancestor_count = ALEA_HIER_SPATIAL_HIT_CHAIN_MAX;
        out->chain_truncated = 1;
    }
    out->ancestor_count = (uint8_t)ancestor_count;
    for (int i = 0; i < ancestor_count; i++) {
        const alea_hier_ray_path_entry_t* ancestor = &path->entries[i];
        out->ancestor_cell_indices[i] = ancestor->cell_index;
        out->ancestor_transforms[i] = ancestor->transform;
        out->ancestor_is_lattice[i] = ancestor->is_lattice;
        out->ancestor_lattice_fill_universes[i] =
            ancestor->lat_fill_universe;
        out->ancestor_lattice_i[i] = ancestor->lat_i;
        out->ancestor_lattice_j[i] = ancestor->lat_j;
        out->ancestor_lattice_k[i] = ancestor->lat_k;
        out->ancestor_lattice_ox[i] = ancestor->lat_ox;
        out->ancestor_lattice_oy[i] = ancestor->lat_oy;
        out->ancestor_lattice_oz[i] = ancestor->lat_oz;
    }
}

typedef struct {
    alea_system_t* sys;
    const alea_slice_view_t* view;
    const alea_transition_slice_critical_tile_t* tile;
    const alea_hier_ray_path_t* path;
    int level;
    uint64_t universe_key;
    critical_curve_t* curves;
    size_t* curve_count;
    size_t curve_capacity;
    critical_curve_t* cell_curves;
    size_t cell_curve_capacity;
    double* breakpoints;
    size_t breakpoint_capacity;
    uint64_t max_active_boundary_tests;
    alea_transition_slice_stats_t* stats;
    int stop_code;
} critical_region_visit_t;

static int emit_cell_degenerate_conic_lines(
    critical_region_visit_t* ctx, const alea_cell_entry_t* cell,
    size_t curve_index, size_t cell_curve_count,
    uint64_t occurrence_key, uint64_t universe_occurrence_key);
static int emit_cell_implicit_conic_pieces(
    critical_region_visit_t* ctx, const alea_cell_entry_t* cell,
    size_t curve_index, size_t cell_curve_count,
    uint64_t occurrence_key, uint64_t universe_occurrence_key);

enum {
    CRITICAL_COLLECT_OK = 0,
    CRITICAL_COLLECT_MAX_CURVES = 1,
    CRITICAL_COLLECT_CHAIN_TRUNCATED = 2
};

static int double_compare(const void* lhs, const void* rhs) {
    const double a = *(const double*)lhs;
    const double b = *(const double*)rhs;
    return (a > b) - (a < b);
}

static size_t compact_sorted_breakpoints(double* values, size_t count) {
    if (!count) return 0;
    size_t kept = 1;
    for (size_t i = 1; i < count; i++) {
        if (fabs(values[i] - values[kept - 1]) > 1e-10)
            values[kept++] = values[i];
    }
    return kept;
}

static void rotate_closed_breakpoints(double* values, size_t count) {
    if (count <= 2) return;
    const double first_physical_break = values[1];
    for (size_t i = 0; i + 1 < count; i++) values[i] = values[i + 1];
    values[count - 1] = first_physical_break + 6.28318530717958647693;
}

static int clip_line_to_tile(
    const alea_curve_2d_t* curve,
    const alea_transition_slice_critical_tile_t* tile,
    double* out_min, double* out_max);

static int make_cell_curve(critical_region_visit_t* ctx,
                           const alea_cell_entry_t* cell,
                           uint32_t surface_index,
                           critical_curve_t* out) {
    if (surface_index >= alea_vec_count(&ctx->sys->surfaces)) return 0;
    const alea_surface_entry_t* surface =
        &ctx->sys->surfaces.data[surface_index];
    if (surface->primitive_id >= alea_vec_count(&ctx->sys->primitives))
        return 0;
    const alea_primitive_entry_t* primitive =
        &ctx->sys->primitives.data[surface->primitive_id];
    alea_primitive_data_t data, transformed;
    alea_primitive_type_t transformed_type;
    const alea_matrix_t* transform =
        &ctx->path->entries[ctx->level].transform;
    if (!alea_primitive_copy_data(ctx->sys, surface->primitive_id, &data) ||
        !alea_primitive_transform(primitive->type, &data, transform->m,
                                  &transformed_type, &transformed))
        return 0;
    alea_curve_2d_t curve;
    if (!alea_intersect_primitive_plane(
            transformed_type, &transformed, &ctx->view->plane, &curve))
        return 0;
    double u0, u1, v0, v1;
    alea_curve_bbox(&curve, &u0, &u1, &v0, &v1);
    if (u1 < ctx->tile->uv_min[0] || u0 > ctx->tile->uv_max[0] ||
        v1 < ctx->tile->uv_min[1] || v0 > ctx->tile->uv_max[1])
        return 0;
    if (curve_is_line(curve.type)) {
        double t_min, t_max;
        if (!clip_line_to_tile(&curve, ctx->tile, &t_min, &t_max)) return 0;
    }
    curve.surface_id = surface->mc_surface_id;
    curve.primitive_id = surface->primitive_id;
    curve.universe_id = cell->universe_id;
    memset(out, 0, sizeof(*out));
    out->curve = curve;
    out->surface_index = surface_index;
    out->cell_id = cell->mc_cell_id;
    out->bbox[0] = u0; out->bbox[1] = u1;
    out->bbox[2] = v0; out->bbox[3] = v1;
    const double center_u = 0.5 * (ctx->tile->uv_min[0] +
                                   ctx->tile->uv_max[0]);
    const double center_v = 0.5 * (ctx->tile->uv_min[1] +
                                   ctx->tile->uv_max[1]);
    if (curve_is_line(curve.type)) {
        const double norm = hypot(curve.data.line.a, curve.data.line.b);
        out->priority = norm > 0.0 ? fabs(
            curve.data.line.a*center_u + curve.data.line.b*center_v +
            curve.data.line.c) / norm : INFINITY;
    } else if (curve.type == ALEA_CURVE_CIRCLE) {
        out->priority = fabs(hypot(
            center_u - curve.data.circle.center[0],
            center_v - curve.data.circle.center[1]) -
            curve.data.circle.radius);
    } else {
        const double du = center_u < u0 ? u0 - center_u :
                          center_u > u1 ? center_u - u1 : 0.0;
        const double dv = center_v < v0 ? v0 - center_v :
                          center_v > v1 ? center_v - v1 : 0.0;
        out->priority = hypot(du, dv);
    }
    return 1;
}

static int clip_line_to_tile(
    const alea_curve_2d_t* curve,
    const alea_transition_slice_critical_tile_t* tile,
    double* out_min, double* out_max) {
    const alea_line_2d_t* line = &curve->data.line;
    double lo = -INFINITY, hi = INFINITY;
    const double mins[2] = {tile->uv_min[0], tile->uv_min[1]};
    const double maxs[2] = {tile->uv_max[0], tile->uv_max[1]};
    for (int axis = 0; axis < 2; axis++) {
        const double p = line->point[axis];
        const double d = line->direction[axis];
        if (fabs(d) <= 1e-15) {
            if (p < mins[axis] || p > maxs[axis]) return 0;
            continue;
        }
        double a = (mins[axis] - p) / d;
        double b = (maxs[axis] - p) / d;
        if (a > b) { const double tmp = a; a = b; b = tmp; }
        if (a > lo) lo = a;
        if (b < hi) hi = b;
    }
    if (curve->type != ALEA_CURVE_LINE) {
        if (isfinite(curve->bounds.t_min) && curve->bounds.t_min > lo)
            lo = curve->bounds.t_min;
        if (isfinite(curve->bounds.t_max) && curve->bounds.t_max < hi)
            hi = curve->bounds.t_max;
    }
    if (!(hi > lo) || !isfinite(lo) || !isfinite(hi)) return 0;
    *out_min = lo; *out_max = hi;
    return 1;
}

static int split_parallel_lines(
    critical_region_visit_t* ctx, const critical_curve_t* source,
    critical_curve_t out[2]) {
    const alea_parallel_lines_2d_t* pair =
        &source->curve.data.parallel_lines;
    const double* points[2] = {pair->point1, pair->point2};
    int count = 0;
    for (int i = 0; i < 2; i++) {
        critical_curve_t item = *source;
        item.curve.type = ALEA_CURVE_LINE;
        item.curve.data.line.point[0] = points[i][0];
        item.curve.data.line.point[1] = points[i][1];
        item.curve.data.line.direction[0] = pair->direction[0];
        item.curve.data.line.direction[1] = pair->direction[1];
        item.curve.data.line.a = -pair->direction[1];
        item.curve.data.line.b = pair->direction[0];
        item.curve.data.line.c = -(
            item.curve.data.line.a*points[i][0] +
            item.curve.data.line.b*points[i][1]);
        item.curve.bounds.t_min = -INFINITY;
        item.curve.bounds.t_max = INFINITY;
        double lo, hi;
        if (!clip_line_to_tile(&item.curve, ctx->tile, &lo, &hi)) continue;
        const double u0 = points[i][0] + lo*pair->direction[0];
        const double v0 = points[i][1] + lo*pair->direction[1];
        const double u1 = points[i][0] + hi*pair->direction[0];
        const double v1 = points[i][1] + hi*pair->direction[1];
        item.bbox[0] = fmin(u0, u1); item.bbox[1] = fmax(u0, u1);
        item.bbox[2] = fmin(v0, v1); item.bbox[3] = fmax(v0, v1);
        const double center_u = 0.5 * (ctx->tile->uv_min[0] +
                                       ctx->tile->uv_max[0]);
        const double center_v = 0.5 * (ctx->tile->uv_min[1] +
                                       ctx->tile->uv_max[1]);
        item.priority = fabs(item.curve.data.line.a*center_u +
                             item.curve.data.line.b*center_v +
                             item.curve.data.line.c);
        out[count++] = item;
    }
    return count;
}

static int split_special_torus(
    critical_region_visit_t* ctx, const critical_curve_t* source,
    critical_curve_t out[2]) {
    alea_curve_2d_t circles[2];
    const int count = torus_circle_components(
        &source->curve.data.torus, circles);
    const double center_u = 0.5*(ctx->tile->uv_min[0]+ctx->tile->uv_max[0]);
    const double center_v = 0.5*(ctx->tile->uv_min[1]+ctx->tile->uv_max[1]);
    for (int i = 0; i < count; i++) {
        out[i] = *source;
        out[i].component_index = (uint8_t)(i+1);
        out[i].curve.type = ALEA_CURVE_CIRCLE;
        out[i].curve.data.circle = circles[i].data.circle;
        const alea_circle_2d_t* circle = &out[i].curve.data.circle;
        out[i].bbox[0] = circle->center[0]-circle->radius;
        out[i].bbox[1] = circle->center[0]+circle->radius;
        out[i].bbox[2] = circle->center[1]-circle->radius;
        out[i].bbox[3] = circle->center[1]+circle->radius;
        out[i].priority = fabs(hypot(
            center_u-circle->center[0], center_v-circle->center[1])-
            circle->radius);
    }
    return count;
}

/* Add the parameters where a line meets another supported slice curve.  An
 * unsupported family returns -1, forcing conservative retention rather than
 * allowing a false inactive classification. */
static int append_line_curve_breaks(const alea_curve_2d_t* line_curve,
                                    const alea_curve_2d_t* other,
                                    double lo, double hi,
                                    double* values, size_t* count,
                                    size_t capacity) {
    const alea_line_2d_t* line = &line_curve->data.line;
#define ADD_T(T) do { \
    const double t_ = (T); \
    if (isfinite(t_) && t_ > lo && t_ < hi) { \
        if (*count == capacity) return -1; \
        values[(*count)++] = t_; \
    } \
} while (0)
    if (curve_is_line(other->type)) {
        const alea_line_2d_t* b = &other->data.line;
        const double det = line->direction[0] * b->direction[1] -
                           line->direction[1] * b->direction[0];
        if (fabs(det) > 1e-14) {
            const double dx = b->point[0] - line->point[0];
            const double dy = b->point[1] - line->point[1];
            ADD_T((dx * b->direction[1] - dy * b->direction[0]) / det);
        }
        return 0;
    }
    if (other->type == ALEA_CURVE_CIRCLE) {
        const double px = line->point[0] - other->data.circle.center[0];
        const double py = line->point[1] - other->data.circle.center[1];
        const double b = px * line->direction[0] +
                         py * line->direction[1];
        const double c = px*px + py*py -
                         other->data.circle.radius * other->data.circle.radius;
        double disc = b*b - c;
        if (disc >= -1e-12) {
            if (disc < 0.0) disc = 0.0;
            const double root = sqrt(disc);
            ADD_T(-b - root);
            if (root > 1e-12) ADD_T(-b + root);
        }
        return 0;
    }
    if (other->type == ALEA_CURVE_ELLIPSE) {
        const alea_ellipse_2d_t* e = &other->data.ellipse;
        const double ct = cos(e->angle), st = sin(e->angle);
        const double x = line->point[0] - e->center[0];
        const double y = line->point[1] - e->center[1];
        const double px = ct*x + st*y, py = -st*x + ct*y;
        const double dx = ct*line->direction[0] + st*line->direction[1];
        const double dy = -st*line->direction[0] + ct*line->direction[1];
        const double aa = dx*dx/(e->semi_a*e->semi_a) +
                          dy*dy/(e->semi_b*e->semi_b);
        const double bb = 2.0*(px*dx/(e->semi_a*e->semi_a) +
                               py*dy/(e->semi_b*e->semi_b));
        const double cc = px*px/(e->semi_a*e->semi_a) +
                          py*py/(e->semi_b*e->semi_b) - 1.0;
        double disc = bb*bb - 4.0*aa*cc;
        if (aa > 0.0 && disc >= -1e-12) {
            if (disc < 0.0) disc = 0.0;
            const double root = sqrt(disc);
            ADD_T((-bb - root) / (2.0*aa));
            if (root > 1e-12) ADD_T((-bb + root) / (2.0*aa));
        }
        return 0;
    }
    if (curve_is_general_conic(other->type)) {
        const alea_conic_2d_t* q = &other->data.conic;
        const double u = line->point[0], v = line->point[1];
        const double du = line->direction[0], dv = line->direction[1];
        const double aa = q->A*du*du + q->B*du*dv + q->C*dv*dv;
        const double bb = 2.0*q->A*u*du + q->B*(u*dv + v*du) +
                          2.0*q->C*v*dv + q->D*du + q->E*dv;
        const double cc = conic_value(q, u, v);
        double roots[2];
        const int n = alea_solve_quadratic(aa, bb, cc, roots);
        for (int i = 0; i < n; i++) ADD_T(roots[i]);
        return 0;
    }
    if (other->type == ALEA_CURVE_QUARTIC) {
        double uv[4][2];
        const int n = torus_line_intersections(
            &other->data.torus, line, uv);
        for (int i = 0; i < n; i++) {
            const double t =
                (uv[i][0]-line->point[0])*line->direction[0] +
                (uv[i][1]-line->point[1])*line->direction[1];
            ADD_T(t);
        }
        return 0;
    }
    if (other->type == ALEA_CURVE_POINT) return 0;
    if (other->type == ALEA_CURVE_POLYGON) {
        const alea_polygon_2d_t* polygon = &other->data.polygon;
        const int edge_count = polygon->closed
            ? polygon->vertex_count : polygon->vertex_count - 1;
        for (int i = 0; i < edge_count; i++) {
            const int j = (i + 1) % polygon->vertex_count;
            const double ex = polygon->vertices[j][0] -
                              polygon->vertices[i][0];
            const double ey = polygon->vertices[j][1] -
                              polygon->vertices[i][1];
            const double det = line->direction[0]*ey -
                               line->direction[1]*ex;
            if (fabs(det) <= 1e-14) continue;
            const double qx = polygon->vertices[i][0] - line->point[0];
            const double qy = polygon->vertices[i][1] - line->point[1];
            const double t = (qx*ey - qy*ex) / det;
            const double s = (qx*line->direction[1] -
                              qy*line->direction[0]) / det;
            if (s >= 0.0 && s <= 1.0) ADD_T(t);
        }
        return 0;
    }
#undef ADD_T
    return -1;
}

static int append_circle_angle(const alea_circle_2d_t* circle,
                               double u, double v,
                               double* values, size_t* count,
                               size_t capacity) {
    if (*count == capacity) return -1;
    values[(*count)++] = positive_angle(atan2(
        v - circle->center[1], u - circle->center[0]));
    return 0;
}

/* Add angular parameters where a circle meets another supported curve. */
static int append_circle_curve_breaks(const alea_curve_2d_t* circle_curve,
                                      const alea_curve_2d_t* other,
                                      double* values, size_t* count,
                                      size_t capacity) {
    const alea_circle_2d_t* circle = &circle_curve->data.circle;
    if (curve_is_line(other->type)) {
        const alea_line_2d_t* line = &other->data.line;
        const double px = line->point[0] - circle->center[0];
        const double py = line->point[1] - circle->center[1];
        const double b = px*line->direction[0] +
                         py*line->direction[1];
        double disc = b*b - (px*px + py*py - circle->radius*circle->radius);
        if (disc < -1e-12) return 0;
        if (disc < 0.0) disc = 0.0;
        const double root = sqrt(disc);
        const double parameters[2] = {-b-root, -b+root};
        const int n = root > 1e-12 ? 2 : 1;
        for (int i = 0; i < n; i++) {
            const double t = parameters[i];
            if (other->type != ALEA_CURVE_LINE &&
                (t < other->bounds.t_min || t > other->bounds.t_max))
                continue;
            if (append_circle_angle(
                    circle, line->point[0] + t*line->direction[0],
                    line->point[1] + t*line->direction[1],
                    values, count, capacity) != 0) return -1;
        }
        return 0;
    }
    if (curve_is_circle(other->type)) {
        const alea_circle_2d_t* b = &other->data.circle;
        const double dx = b->center[0]-circle->center[0];
        const double dy = b->center[1]-circle->center[1];
        const double d = hypot(dx, dy);
        if (d <= 1e-14 || d > circle->radius+b->radius+1e-12 ||
            d < fabs(circle->radius-b->radius)-1e-12) return 0;
        const double x = (circle->radius*circle->radius -
                          b->radius*b->radius + d*d) / (2.0*d);
        double h2 = circle->radius*circle->radius - x*x;
        if (h2 < -1e-12) return 0;
        if (h2 < 0.0) h2 = 0.0;
        const double h = sqrt(h2), ux = dx/d, uy = dy/d;
        const double pu = circle->center[0] + x*ux;
        const double pv = circle->center[1] + x*uy;
        if (append_circle_angle(circle, pu-h*uy, pv+h*ux,
                                values, count, capacity) != 0) return -1;
        if (h > 1e-12 && append_circle_angle(
                circle, pu+h*uy, pv-h*ux,
                values, count, capacity) != 0) return -1;
        return 0;
    }
    if (curve_is_ellipse(other->type)) {
        double uv[4][2];
        const int n = closed_conic_intersections(circle_curve, other, uv);
        if (n < 0) return -1;
        for (int i = 0; i < n; i++)
            if (append_circle_angle(circle, uv[i][0], uv[i][1],
                                    values, count, capacity) != 0)
                return -1;
        return 0;
    }
    if (curve_is_general_conic(other->type)) {
        double uv[4][2];
        const int n = parametric_closed_conic_intersections(
            circle_curve, &other->data.conic, uv);
        if (n < 0) return -1;
        for (int i = 0; i < n; i++)
            if (append_circle_angle(circle, uv[i][0], uv[i][1],
                                    values, count, capacity) != 0)
                return -1;
        return 0;
    }
    if (other->type == ALEA_CURVE_QUARTIC) {
        double uv[8][2];
        alea_curve_2d_t components[2];
        const int component_count = torus_circle_components(
            &other->data.torus, components);
        if (component_count == 0) {
            const int n = torus_closed_conic_intersections(
                &other->data.torus, circle_curve, uv);
            if (n < 0) return -1;
            for (int i = 0; i < n; i++)
                if (append_circle_angle(
                        circle, uv[i][0], uv[i][1],
                        values, count, capacity) != 0) return -1;
        } else {
            for (int component = 0; component < component_count; component++) {
                double component_uv[4][2];
                const int n = closed_conic_intersections(
                    circle_curve, &components[component], component_uv);
                if (n < 0) return -1;
                for (int i = 0; i < n; i++)
                    if (append_circle_angle(
                            circle, component_uv[i][0], component_uv[i][1],
                            values, count, capacity) != 0) return -1;
            }
        }
        return 0;
    }
    if (other->type == ALEA_CURVE_POINT) return 0;
    if (other->type == ALEA_CURVE_POLYGON) {
        const alea_polygon_2d_t* polygon = &other->data.polygon;
        const int edge_count = polygon->closed
            ? polygon->vertex_count : polygon->vertex_count - 1;
        for (int i = 0; i < edge_count; i++) {
            alea_curve_2d_t edge;
            memset(&edge, 0, sizeof(edge));
            edge.type = ALEA_CURVE_LINE_SEGMENT;
            edge.data.line.point[0] = polygon->vertices[i][0];
            edge.data.line.point[1] = polygon->vertices[i][1];
            const int j = (i + 1) % polygon->vertex_count;
            double dx = polygon->vertices[j][0] - edge.data.line.point[0];
            double dy = polygon->vertices[j][1] - edge.data.line.point[1];
            const double length = hypot(dx, dy);
            if (!(length > 0.0)) continue;
            edge.data.line.direction[0] = dx / length;
            edge.data.line.direction[1] = dy / length;
            edge.bounds.t_min = 0.0; edge.bounds.t_max = length;
            if (append_circle_curve_breaks(
                    circle_curve, &edge, values, count, capacity) != 0)
                return -1;
        }
        return 0;
    }
    return -1;
}

static int append_ellipse_curve_breaks(const alea_curve_2d_t* ellipse_curve,
                                       const alea_curve_2d_t* other,
                                       double* values, size_t* count,
                                       size_t capacity) {
    const alea_ellipse_2d_t* ellipse = &ellipse_curve->data.ellipse;
    if (curve_is_line(other->type)) {
        const alea_line_2d_t* line = &other->data.line;
        const double ct = cos(ellipse->angle), st = sin(ellipse->angle);
        const double x = line->point[0]-ellipse->center[0];
        const double y = line->point[1]-ellipse->center[1];
        const double px = ct*x + st*y, py = -st*x + ct*y;
        const double dx = ct*line->direction[0] + st*line->direction[1];
        const double dy = -st*line->direction[0] + ct*line->direction[1];
        const double aa = dx*dx/(ellipse->semi_a*ellipse->semi_a) +
                          dy*dy/(ellipse->semi_b*ellipse->semi_b);
        const double bb = 2.0*(px*dx/(ellipse->semi_a*ellipse->semi_a) +
                               py*dy/(ellipse->semi_b*ellipse->semi_b));
        const double cc = px*px/(ellipse->semi_a*ellipse->semi_a) +
                          py*py/(ellipse->semi_b*ellipse->semi_b) - 1.0;
        double disc = bb*bb - 4.0*aa*cc;
        if (!(aa > 0.0) || disc < -1e-12) return 0;
        if (disc < 0.0) disc = 0.0;
        const double root = sqrt(disc);
        const double parameters[2] = {
            (-bb-root)/(2.0*aa), (-bb+root)/(2.0*aa)};
        const int n = root > 1e-12 ? 2 : 1;
        for (int i = 0; i < n; i++) {
            const double t = parameters[i];
            if (other->type != ALEA_CURVE_LINE &&
                (t < other->bounds.t_min || t > other->bounds.t_max))
                continue;
            if (*count == capacity) return -1;
            values[(*count)++] = ellipse_parameter(
                ellipse, line->point[0] + t*line->direction[0],
                line->point[1] + t*line->direction[1]);
        }
        return 0;
    }
    if (curve_is_closed_conic(other->type)) {
        double uv[4][2];
        const int n = closed_conic_intersections(ellipse_curve, other, uv);
        if (n < 0) return -1;
        for (int i = 0; i < n; i++) {
            if (*count == capacity) return -1;
            values[(*count)++] = ellipse_parameter(
                ellipse, uv[i][0], uv[i][1]);
        }
        return 0;
    }
    if (curve_is_general_conic(other->type)) {
        double uv[4][2];
        const int n = parametric_closed_conic_intersections(
            ellipse_curve, &other->data.conic, uv);
        if (n < 0) return -1;
        for (int i = 0; i < n; i++) {
            if (*count == capacity) return -1;
            values[(*count)++] = ellipse_parameter(
                ellipse, uv[i][0], uv[i][1]);
        }
        return 0;
    }
    if (other->type == ALEA_CURVE_QUARTIC) {
        double uv[8][2];
        alea_curve_2d_t components[2];
        const int component_count = torus_circle_components(
            &other->data.torus, components);
        if (component_count == 0) {
            const int n = torus_closed_conic_intersections(
                &other->data.torus, ellipse_curve, uv);
            if (n < 0) return -1;
            for (int i = 0; i < n; i++) {
                if (*count == capacity) return -1;
                values[(*count)++] = ellipse_parameter(
                    ellipse, uv[i][0], uv[i][1]);
            }
        } else {
            for (int component = 0; component < component_count; component++) {
                double component_uv[4][2];
                const int n = closed_conic_intersections(
                    ellipse_curve, &components[component], component_uv);
                if (n < 0) return -1;
                for (int i = 0; i < n; i++) {
                    if (*count == capacity) return -1;
                    values[(*count)++] = ellipse_parameter(
                        ellipse, component_uv[i][0], component_uv[i][1]);
                }
            }
        }
        return 0;
    }
    if (other->type == ALEA_CURVE_POINT) return 0;
    if (other->type == ALEA_CURVE_POLYGON) {
        const alea_polygon_2d_t* polygon = &other->data.polygon;
        const int edge_count = polygon->closed
            ? polygon->vertex_count : polygon->vertex_count - 1;
        for (int i = 0; i < edge_count; i++) {
            const int j = (i + 1) % polygon->vertex_count;
            alea_curve_2d_t edge;
            memset(&edge, 0, sizeof(edge));
            edge.type = ALEA_CURVE_LINE_SEGMENT;
            edge.data.line.point[0] = polygon->vertices[i][0];
            edge.data.line.point[1] = polygon->vertices[i][1];
            const double dx = polygon->vertices[j][0] - edge.data.line.point[0];
            const double dy = polygon->vertices[j][1] - edge.data.line.point[1];
            const double length = hypot(dx, dy);
            if (!(length > 0.0)) continue;
            edge.data.line.direction[0] = dx/length;
            edge.data.line.direction[1] = dy/length;
            edge.bounds.t_min = 0.0; edge.bounds.t_max = length;
            if (append_ellipse_curve_breaks(
                    ellipse_curve, &edge, values, count, capacity) != 0)
                return -1;
        }
        return 0;
    }
    return -1;
}

static int active_partition_pair_supported(alea_curve_type_t source,
                                           alea_curve_type_t other) {
    if (curve_is_line(source))
        return curve_is_line(other) || curve_is_closed_conic(other) ||
               curve_is_general_conic(other) ||
               other == ALEA_CURVE_QUARTIC ||
               other == ALEA_CURVE_POINT || other == ALEA_CURVE_POLYGON;
    if (curve_is_closed_conic(source))
        return curve_is_line(other) || curve_is_closed_conic(other) ||
               curve_is_general_conic(other) ||
               other == ALEA_CURVE_QUARTIC ||
               other == ALEA_CURVE_POINT || other == ALEA_CURVE_POLYGON;
    if (curve_is_general_conic(source))
        return curve_is_line(other) || curve_is_closed_conic(other) ||
               curve_is_general_conic(other) ||
               other == ALEA_CURVE_QUARTIC ||
               other == ALEA_CURVE_POINT || other == ALEA_CURVE_POLYGON;
    if (source == ALEA_CURVE_QUARTIC)
        return curve_is_line(other) || curve_is_closed_conic(other) ||
               curve_is_general_conic(other) ||
               other == ALEA_CURVE_QUARTIC ||
               other == ALEA_CURVE_POINT || other == ALEA_CURVE_POLYGON;
    if (source == ALEA_CURVE_POLYGON)
        return curve_is_line(other) || curve_is_closed_conic(other) ||
               curve_is_general_conic(other) ||
               other == ALEA_CURVE_QUARTIC ||
               other == ALEA_CURVE_POINT || other == ALEA_CURVE_POLYGON;
    if (source == ALEA_CURVE_POINT) return 1;
    return 0;
}

static void record_active_unsupported(
    alea_curve_type_t type, alea_transition_slice_stats_t* stats) {
    if (type == ALEA_CURVE_PARABOLA)
        stats->critical_active_unsupported_parabola_fallbacks++;
    else if (type == ALEA_CURVE_HYPERBOLA)
        stats->critical_active_unsupported_hyperbola_fallbacks++;
    else if (type == ALEA_CURVE_QUARTIC)
        stats->critical_active_unsupported_quartic_fallbacks++;
    else if (type == ALEA_CURVE_POLYGON)
        stats->critical_active_unsupported_polygon_fallbacks++;
    else if (type == ALEA_CURVE_POINT)
        stats->critical_active_unsupported_point_fallbacks++;
    else
        stats->critical_active_unsupported_other_fallbacks++;
}

static void record_active_evaluation_failure(
    alea_curve_type_t type, alea_transition_slice_stats_t* stats) {
    if (curve_is_line(type))
        stats->critical_active_evaluation_line_fallbacks++;
    else if (curve_is_closed_conic(type))
        stats->critical_active_evaluation_closed_conic_fallbacks++;
    else if (curve_is_general_conic(type))
        stats->critical_active_evaluation_general_conic_fallbacks++;
    else if (type == ALEA_CURVE_QUARTIC)
        stats->critical_active_evaluation_quartic_fallbacks++;
    else if (type == ALEA_CURVE_POLYGON)
        stats->critical_active_evaluation_polygon_fallbacks++;
    else
        stats->critical_active_evaluation_other_fallbacks++;
}

static int retain_ranked_curve(critical_region_visit_t* ctx,
                               const critical_curve_t* item,
                               uint64_t occurrence_key,
                               uint64_t universe_occurrence_key);

static void emit_active_line_segment(
    critical_region_visit_t* ctx, const critical_curve_t* item,
    const alea_line_2d_t* line, double t_min, double t_max,
    uint64_t occurrence_key, uint64_t universe_occurrence_key) {
    critical_curve_t segment = *item;
    segment.curve.type = ALEA_CURVE_LINE_SEGMENT;
    segment.curve.bounds.t_min = t_min;
    segment.curve.bounds.t_max = t_max;
    const double u0 = line->point[0] + t_min*line->direction[0];
    const double v0 = line->point[1] + t_min*line->direction[1];
    const double u1 = line->point[0] + t_max*line->direction[0];
    const double v1 = line->point[1] + t_max*line->direction[1];
    segment.bbox[0] = fmin(u0, u1);
    segment.bbox[1] = fmax(u0, u1);
    segment.bbox[2] = fmin(v0, v1);
    segment.bbox[3] = fmax(v0, v1);
    const double mid = 0.5 * (t_min + t_max);
    segment.has_active_point = 1;
    segment.active_uv[0] = line->point[0] + mid*line->direction[0];
    segment.active_uv[1] = line->point[1] + mid*line->direction[1];
    const double center_u = 0.5 * (ctx->tile->uv_min[0] +
                                   ctx->tile->uv_max[0]);
    const double center_v = 0.5 * (ctx->tile->uv_min[1] +
                                   ctx->tile->uv_max[1]);
    double center_t =
        (center_u-line->point[0])*line->direction[0] +
        (center_v-line->point[1])*line->direction[1];
    if (center_t < t_min) center_t = t_min;
    if (center_t > t_max) center_t = t_max;
    segment.priority = hypot(
        center_u - (line->point[0] + center_t*line->direction[0]),
        center_v - (line->point[1] + center_t*line->direction[1]));
    retain_ranked_curve(ctx, &segment, occurrence_key,
                        universe_occurrence_key);
    ctx->stats->critical_active_segments++;
}

static void emit_active_circle_arc(
    critical_region_visit_t* ctx, const critical_curve_t* item,
    double theta_min, double theta_max,
    uint64_t occurrence_key, uint64_t universe_occurrence_key) {
    critical_curve_t arc = *item;
    arc.curve.type = ALEA_CURVE_ARC;
    arc.curve.bounds.theta_start = theta_min;
    arc.curve.bounds.theta_end = theta_max;
    arc.curve.bounds.t_min = theta_min;
    arc.curve.bounds.t_max = theta_max;
    arc.bbox[0] = fmax(item->bbox[0], ctx->tile->uv_min[0]);
    arc.bbox[1] = fmin(item->bbox[1], ctx->tile->uv_max[0]);
    arc.bbox[2] = fmax(item->bbox[2], ctx->tile->uv_min[1]);
    arc.bbox[3] = fmin(item->bbox[3], ctx->tile->uv_max[1]);
    const double mid = 0.5 * (theta_min + theta_max);
    const alea_circle_2d_t* circle = &arc.curve.data.circle;
    arc.has_active_point = 1;
    arc.active_uv[0] = circle->center[0] + circle->radius*cos(mid);
    arc.active_uv[1] = circle->center[1] + circle->radius*sin(mid);
    const double center_u = 0.5 * (ctx->tile->uv_min[0] +
                                   ctx->tile->uv_max[0]);
    const double center_v = 0.5 * (ctx->tile->uv_min[1] +
                                   ctx->tile->uv_max[1]);
    arc.priority = hypot(center_u - arc.active_uv[0],
                         center_v - arc.active_uv[1]);
    retain_ranked_curve(ctx, &arc, occurrence_key, universe_occurrence_key);
    ctx->stats->critical_active_segments++;
}

static int emit_cell_circle_arcs(
    critical_region_visit_t* ctx, const alea_cell_entry_t* cell,
    size_t curve_index, size_t cell_curve_count,
    uint64_t occurrence_key, uint64_t universe_occurrence_key) {
    const critical_curve_t* item = &ctx->cell_curves[curve_index];
    if (item->curve.type != ALEA_CURVE_CIRCLE) return 2;
    const alea_circle_2d_t* circle = &item->curve.data.circle;
    const double turn = 6.28318530717958647693;
    size_t break_count = 2;
    ctx->breakpoints[0] = 0.0; ctx->breakpoints[1] = turn;
    alea_curve_2d_t edge;
    memset(&edge, 0, sizeof(edge));
    edge.type = ALEA_CURVE_LINE;
    for (int axis = 0; axis < 2; axis++) {
        for (int side = 0; side < 2; side++) {
            const double coordinate = side
                ? ctx->tile->uv_max[axis] : ctx->tile->uv_min[axis];
            edge.data.line.point[0] = axis ? 0.0 : coordinate;
            edge.data.line.point[1] = axis ? coordinate : 0.0;
            edge.data.line.direction[0] = axis ? 1.0 : 0.0;
            edge.data.line.direction[1] = axis ? 0.0 : 1.0;
            if (append_circle_curve_breaks(
                    &item->curve, &edge, ctx->breakpoints, &break_count,
                    ctx->breakpoint_capacity) != 0) return 2;
        }
    }
    for (size_t i = 0; i < cell_curve_count; i++) {
        if (i == curve_index) continue;
        if (append_circle_curve_breaks(
                &item->curve, &ctx->cell_curves[i].curve,
                ctx->breakpoints, &break_count,
                ctx->breakpoint_capacity) != 0) return 2;
    }
    qsort(ctx->breakpoints, break_count, sizeof(double), double_compare);
    break_count = compact_sorted_breakpoints(ctx->breakpoints, break_count);
    rotate_closed_breakpoints(ctx->breakpoints, break_count);
    const double tile_scale = fmax(1.0, fmax(
        ctx->tile->uv_max[0] - ctx->tile->uv_min[0],
        ctx->tile->uv_max[1] - ctx->tile->uv_min[1]));
    int emitted = 0, active_run = 0;
    double active_min = 0.0, active_max = 0.0;
    for (size_t i = 1; i < break_count; i++) {
        const double span = ctx->breakpoints[i] - ctx->breakpoints[i - 1];
        if (!(span > 1e-10)) continue;
        const double mid = 0.5 * (ctx->breakpoints[i - 1] +
                                  ctx->breakpoints[i]);
        const double radial_u = cos(mid), radial_v = sin(mid);
        const double u0 = circle->center[0] + circle->radius*radial_u;
        const double v0 = circle->center[1] + circle->radius*radial_v;
        int active = 0;
        if (point_in_tile(ctx->tile, u0, v0, 1e-12 * tile_scale)) {
            if (ctx->stats->critical_active_boundary_tests >=
                ctx->max_active_boundary_tests) return 2;
            double radius = 1e-6 * tile_scale;
            if (!(radius > 64.0 * RAY_EPSILON))
                radius = 64.0 * RAY_EPSILON;
            const double arc_span = circle->radius * span;
            if (radius > 0.125 * arc_span) radius = 0.125 * arc_span;
            int contains[2];
            for (int side = 0; side < 2; side++) {
                const double sign = side ? 1.0 : -1.0;
                double local[3];
                tile_world_point(ctx->view,
                    u0 + sign*radius*radial_u,
                    v0 + sign*radius*radial_v, local);
                alea_matrix_transform_point_inverse(
                    &ctx->path->entries[ctx->level].transform,
                    &local[0], &local[1], &local[2]);
                contains[side] = alea_contains_point(
                    ctx->sys, cell->root_node_id,
                    local[0], local[1], local[2]);
            }
            ctx->stats->critical_active_boundary_tests++;
            active = contains[0] != contains[1];
        }
        if (active) {
            if (!active_run) {
                active_min = ctx->breakpoints[i - 1];
                active_run = 1;
            }
            active_max = ctx->breakpoints[i];
        } else if (active_run) {
            emit_active_circle_arc(ctx, item, active_min, active_max,
                                   occurrence_key, universe_occurrence_key);
            emitted = 1; active_run = 0;
        }
    }
    if (active_run) {
        emit_active_circle_arc(ctx, item, active_min, active_max,
                               occurrence_key, universe_occurrence_key);
        emitted = 1;
    }
    return emitted;
}

static void emit_active_ellipse_arc(
    critical_region_visit_t* ctx, const critical_curve_t* item,
    double theta_min, double theta_max,
    uint64_t occurrence_key, uint64_t universe_occurrence_key) {
    critical_curve_t arc = *item;
    arc.curve.type = ALEA_CURVE_ELLIPSE_ARC;
    arc.curve.bounds.theta_start = theta_min;
    arc.curve.bounds.theta_end = theta_max;
    arc.curve.bounds.t_min = theta_min;
    arc.curve.bounds.t_max = theta_max;
    arc.bbox[0] = fmax(item->bbox[0], ctx->tile->uv_min[0]);
    arc.bbox[1] = fmin(item->bbox[1], ctx->tile->uv_max[0]);
    arc.bbox[2] = fmax(item->bbox[2], ctx->tile->uv_min[1]);
    arc.bbox[3] = fmin(item->bbox[3], ctx->tile->uv_max[1]);
    const double mid = 0.5 * (theta_min + theta_max);
    arc.has_active_point = alea_curve_eval(
        &arc.curve, mid, &arc.active_uv[0], &arc.active_uv[1]);
    const double center_u = 0.5 * (ctx->tile->uv_min[0] +
                                   ctx->tile->uv_max[0]);
    const double center_v = 0.5 * (ctx->tile->uv_min[1] +
                                   ctx->tile->uv_max[1]);
    arc.priority = hypot(center_u - arc.active_uv[0],
                         center_v - arc.active_uv[1]);
    retain_ranked_curve(ctx, &arc, occurrence_key, universe_occurrence_key);
    ctx->stats->critical_active_segments++;
}

static int emit_cell_ellipse_arcs(
    critical_region_visit_t* ctx, const alea_cell_entry_t* cell,
    size_t curve_index, size_t cell_curve_count,
    uint64_t occurrence_key, uint64_t universe_occurrence_key) {
    const critical_curve_t* item = &ctx->cell_curves[curve_index];
    if (item->curve.type != ALEA_CURVE_ELLIPSE) return 2;
    const double turn = 6.28318530717958647693;
    size_t break_count = 2;
    ctx->breakpoints[0] = 0.0; ctx->breakpoints[1] = turn;
    alea_curve_2d_t edge;
    memset(&edge, 0, sizeof(edge));
    edge.type = ALEA_CURVE_LINE;
    for (int axis = 0; axis < 2; axis++) {
        for (int side = 0; side < 2; side++) {
            const double coordinate = side
                ? ctx->tile->uv_max[axis] : ctx->tile->uv_min[axis];
            edge.data.line.point[0] = axis ? 0.0 : coordinate;
            edge.data.line.point[1] = axis ? coordinate : 0.0;
            edge.data.line.direction[0] = axis ? 1.0 : 0.0;
            edge.data.line.direction[1] = axis ? 0.0 : 1.0;
            if (append_ellipse_curve_breaks(
                    &item->curve, &edge, ctx->breakpoints, &break_count,
                    ctx->breakpoint_capacity) != 0) return 2;
        }
    }
    for (size_t i = 0; i < cell_curve_count; i++) {
        if (i == curve_index) continue;
        if (append_ellipse_curve_breaks(
                &item->curve, &ctx->cell_curves[i].curve,
                ctx->breakpoints, &break_count,
                ctx->breakpoint_capacity) != 0) return 2;
    }
    qsort(ctx->breakpoints, break_count, sizeof(double), double_compare);
    break_count = compact_sorted_breakpoints(ctx->breakpoints, break_count);
    rotate_closed_breakpoints(ctx->breakpoints, break_count);
    const double tile_scale = fmax(1.0, fmax(
        ctx->tile->uv_max[0] - ctx->tile->uv_min[0],
        ctx->tile->uv_max[1] - ctx->tile->uv_min[1]));
    int emitted = 0, active_run = 0;
    double active_min = 0.0, active_max = 0.0;
    for (size_t i = 1; i < break_count; i++) {
        const double span = ctx->breakpoints[i] - ctx->breakpoints[i - 1];
        if (!(span > 1e-10)) continue;
        const double mid = 0.5 * (ctx->breakpoints[i - 1] +
                                  ctx->breakpoints[i]);
        double u, v, normal[2];
        int active = 0;
        if (alea_curve_eval(&item->curve, mid, &u, &v) &&
            point_in_tile(ctx->tile, u, v, 1e-12*tile_scale) &&
            curve_normal_at(&item->curve, (double[2]){u, v}, normal)) {
            if (ctx->stats->critical_active_boundary_tests >=
                ctx->max_active_boundary_tests) return 2;
            double radius = fmax(1e-6*tile_scale, 64.0*RAY_EPSILON);
            int contains[2];
            for (int side = 0; side < 2; side++) {
                const double sign = side ? 1.0 : -1.0;
                double local[3];
                tile_world_point(ctx->view, u + sign*radius*normal[0],
                                 v + sign*radius*normal[1], local);
                alea_matrix_transform_point_inverse(
                    &ctx->path->entries[ctx->level].transform,
                    &local[0], &local[1], &local[2]);
                contains[side] = alea_contains_point(
                    ctx->sys, cell->root_node_id,
                    local[0], local[1], local[2]);
            }
            ctx->stats->critical_active_boundary_tests++;
            active = contains[0] != contains[1];
        }
        if (active) {
            if (!active_run) { active_min = ctx->breakpoints[i-1]; active_run = 1; }
            active_max = ctx->breakpoints[i];
        } else if (active_run) {
            emit_active_ellipse_arc(ctx, item, active_min, active_max,
                                    occurrence_key, universe_occurrence_key);
            emitted = 1; active_run = 0;
        }
    }
    if (active_run) {
        emit_active_ellipse_arc(ctx, item, active_min, active_max,
                                occurrence_key, universe_occurrence_key);
        emitted = 1;
    }
    return emitted;
}

static int append_open_conic_parameter(
    const critical_curve_t* item, int wanted_branch, double u, double v,
    double* values, size_t* count, size_t capacity) {
    double parameter;
    int branch;
    if (!open_conic_parameter(item, u, v, &parameter, &branch)) return -1;
    if (branch != wanted_branch) return 0;
    if (*count == capacity) return -1;
    values[(*count)++] = parameter;
    return 0;
}

static int append_open_conic_curve_breaks(
    const critical_curve_t* item, const alea_curve_2d_t* other,
    int branch, double* values, size_t* count, size_t capacity) {
    double uv[16][2];
    int n = 0;
    if (curve_is_line(other->type)) {
        const alea_line_2d_t* line = &other->data.line;
        const alea_conic_2d_t* q = &item->curve.data.conic;
        const double u = line->point[0], v = line->point[1];
        const double du = line->direction[0], dv = line->direction[1];
        double roots[2];
        n = alea_solve_quadratic(
            q->A*du*du + q->B*du*dv + q->C*dv*dv,
            2.0*q->A*u*du + q->B*(u*dv+v*du) +
                2.0*q->C*v*dv + q->D*du + q->E*dv,
            conic_value(q, u, v), roots);
        for (int i = 0; i < n; i++) {
            const double t = roots[i];
            if (other->type != ALEA_CURVE_LINE &&
                (t < other->bounds.t_min || t > other->bounds.t_max))
                continue;
            if (append_open_conic_parameter(
                    item, branch, u+t*du, v+t*dv,
                    values, count, capacity) != 0)
                return -1;
        }
        return 0;
    }
    if (curve_is_closed_conic(other->type)) {
        n = parametric_closed_conic_intersections(
            other, &item->curve.data.conic, uv);
    } else if (curve_is_general_conic(other->type)) {
        n = general_conic_intersections(
            &item->curve.data.conic, &other->data.conic, uv);
    } else if (other->type == ALEA_CURVE_QUARTIC) {
        alea_curve_2d_t components[2];
        const int component_count = torus_circle_components(
            &other->data.torus, components);
        if (component_count == 0) {
            n = torus_open_conic_intersections(
                &other->data.torus, item, uv);
        } else {
            for (int component = 0; component < component_count; component++) {
                double component_uv[4][2];
                const int component_intersections =
                    parametric_closed_conic_intersections(
                        &components[component], &item->curve.data.conic,
                        component_uv);
                if (component_intersections < 0) return -1;
                for (int i = 0; i < component_intersections; i++)
                    if (append_open_conic_parameter(
                            item, branch, component_uv[i][0],
                            component_uv[i][1], values, count,
                            capacity) != 0) return -1;
            }
            return 0;
        }
    } else if (other->type == ALEA_CURVE_POINT) {
        return 0;
    } else if (other->type == ALEA_CURVE_POLYGON) {
        const alea_polygon_2d_t* polygon = &other->data.polygon;
        const int edge_count = polygon->closed
            ? polygon->vertex_count : polygon->vertex_count-1;
        for (int i = 0; i < edge_count; i++) {
            const int j = (i+1)%polygon->vertex_count;
            alea_curve_2d_t edge;
            memset(&edge, 0, sizeof(edge));
            edge.type = ALEA_CURVE_LINE_SEGMENT;
            edge.data.line.point[0] = polygon->vertices[i][0];
            edge.data.line.point[1] = polygon->vertices[i][1];
            const double du = polygon->vertices[j][0]-polygon->vertices[i][0];
            const double dv = polygon->vertices[j][1]-polygon->vertices[i][1];
            const double length = hypot(du, dv);
            if (!(length > 0.0)) continue;
            edge.data.line.direction[0] = du/length;
            edge.data.line.direction[1] = dv/length;
            edge.bounds.t_min = 0.0; edge.bounds.t_max = length;
            if (append_open_conic_curve_breaks(
                    item, &edge, branch, values, count, capacity) != 0)
                return -1;
        }
        return 0;
    } else {
        return -1;
    }
    if (n < 0) return -1;
    for (int i = 0; i < n; i++)
        if (append_open_conic_parameter(
                item, branch, uv[i][0], uv[i][1],
                values, count, capacity) != 0)
            return -1;
    return 0;
}

static int open_conic_piece_bbox(
    const critical_curve_t* item, double parameter_min,
    double parameter_max, int branch, double bbox[4]) {
    double parameters[4] = {parameter_min, parameter_max, 0.0, 0.0};
    int parameter_count = 2;
    for (int coordinate = 0; coordinate < 2; coordinate++) {
        double extremum = 0.0;
        int has_extremum = 0;
        if (item->curve.type == ALEA_CURVE_PARABOLA) {
            const double denominator = 2.0*item->parameter_axis_y[coordinate]*
                                       item->parameter_scale[0];
            if (fabs(denominator) > 1e-15) {
                extremum = -item->parameter_axis_x[coordinate]/denominator;
                has_extremum = 1;
            }
        } else if (item->curve.type == ALEA_CURVE_HYPERBOLA) {
            const double denominator = (branch < 0 ? -1.0 : 1.0)*
                item->parameter_scale[0]*item->parameter_axis_x[coordinate];
            if (fabs(denominator) > 1e-15) {
                const double ratio = -item->parameter_scale[1]*
                    item->parameter_axis_y[coordinate]/denominator;
                if (fabs(ratio) < 1.0) {
                    extremum = atanh(ratio);
                    has_extremum = 1;
                }
            }
        }
        if (has_extremum && extremum > parameter_min &&
            extremum < parameter_max && parameter_count < 4)
            parameters[parameter_count++] = extremum;
    }
    bbox[0] = bbox[2] = INFINITY;
    bbox[1] = bbox[3] = -INFINITY;
    for (int i = 0; i < parameter_count; i++) {
        double u, v;
        if (!open_conic_eval(item, parameters[i], branch, &u, &v)) return 0;
        bbox[0] = fmin(bbox[0], u); bbox[1] = fmax(bbox[1], u);
        bbox[2] = fmin(bbox[2], v); bbox[3] = fmax(bbox[3], v);
    }
    return isfinite(bbox[0]) && isfinite(bbox[1]) &&
           isfinite(bbox[2]) && isfinite(bbox[3]);
}

static void emit_active_open_conic_piece(
    critical_region_visit_t* ctx, const critical_curve_t* canonical,
    double parameter_min, double parameter_max, int branch,
    uint64_t occurrence_key, uint64_t universe_occurrence_key) {
    critical_curve_t piece = *canonical;
    piece.has_parameter_domain = 1;
    piece.conic_branch = (int8_t)branch;
    piece.parameter_min = parameter_min;
    piece.parameter_max = parameter_max;
    const double mid = 0.5*(parameter_min+parameter_max);
    piece.has_active_point = open_conic_eval(
        &piece, mid, branch, &piece.active_uv[0], &piece.active_uv[1]);
    if (!open_conic_piece_bbox(
            &piece, parameter_min, parameter_max, branch, piece.bbox))
        memcpy(piece.bbox, canonical->bbox, sizeof(piece.bbox));
    const double padding = 1e-12*fmax(1.0, fmax(
        ctx->tile->uv_max[0]-ctx->tile->uv_min[0],
        ctx->tile->uv_max[1]-ctx->tile->uv_min[1]));
    piece.bbox[0] = fmax(piece.bbox[0]-padding, ctx->tile->uv_min[0]);
    piece.bbox[1] = fmin(piece.bbox[1]+padding, ctx->tile->uv_max[0]);
    piece.bbox[2] = fmax(piece.bbox[2]-padding, ctx->tile->uv_min[1]);
    piece.bbox[3] = fmin(piece.bbox[3]+padding, ctx->tile->uv_max[1]);
    const double center_u = 0.5*(ctx->tile->uv_min[0]+ctx->tile->uv_max[0]);
    const double center_v = 0.5*(ctx->tile->uv_min[1]+ctx->tile->uv_max[1]);
    piece.priority = hypot(center_u-piece.active_uv[0],
                           center_v-piece.active_uv[1]);
    retain_ranked_curve(ctx, &piece, occurrence_key,
                        universe_occurrence_key);
    ctx->stats->critical_active_segments++;
}

static int emit_cell_open_conic_pieces(
    critical_region_visit_t* ctx, const alea_cell_entry_t* cell,
    size_t curve_index, size_t cell_curve_count,
    uint64_t occurrence_key, uint64_t universe_occurrence_key) {
    critical_curve_t canonical = ctx->cell_curves[curve_index];
    if (!curve_is_general_conic(canonical.curve.type) ||
        !canonicalize_open_conic(&canonical.curve, &canonical)) {
        const int degenerate = emit_cell_degenerate_conic_lines(
            ctx, cell, curve_index, cell_curve_count,
            occurrence_key, universe_occurrence_key);
        if (degenerate != 2) return degenerate;
        const int implicit = emit_cell_implicit_conic_pieces(
            ctx, cell, curve_index, cell_curve_count,
            occurrence_key, universe_occurrence_key);
        if (implicit != 2) return implicit;
        ctx->stats->critical_active_open_conic_canonical_fallbacks++;
        return 2;
    }
    const int branches[2] = {canonical.curve.type == ALEA_CURVE_HYPERBOLA
        ? -1 : 0, 1};
    const int branch_count = canonical.curve.type == ALEA_CURVE_HYPERBOLA
        ? 2 : 1;
    const double tile_scale = fmax(1.0, fmax(
        ctx->tile->uv_max[0]-ctx->tile->uv_min[0],
        ctx->tile->uv_max[1]-ctx->tile->uv_min[1]));
    int emitted = 0;
    for (int bi = 0; bi < branch_count; bi++) {
        const int branch = branches[bi];
        size_t break_count = 0;
        alea_curve_2d_t edge;
        memset(&edge, 0, sizeof(edge));
        edge.type = ALEA_CURVE_LINE;
        for (int axis = 0; axis < 2; axis++) {
            for (int side = 0; side < 2; side++) {
                const double coordinate = side
                    ? ctx->tile->uv_max[axis] : ctx->tile->uv_min[axis];
                edge.data.line.point[0] = axis ? 0.0 : coordinate;
                edge.data.line.point[1] = axis ? coordinate : 0.0;
                edge.data.line.direction[0] = axis ? 1.0 : 0.0;
                edge.data.line.direction[1] = axis ? 0.0 : 1.0;
                if (append_open_conic_curve_breaks(
                        &canonical, &edge, branch, ctx->breakpoints,
                        &break_count, ctx->breakpoint_capacity) != 0) {
                    ctx->stats->
                        critical_active_open_conic_breakpoint_fallbacks++;
                    return 2;
                }
            }
        }
        for (size_t i = 0; i < cell_curve_count; i++) {
            if (i == curve_index) continue;
            if (append_open_conic_curve_breaks(
                    &canonical, &ctx->cell_curves[i].curve, branch,
                    ctx->breakpoints, &break_count,
                    ctx->breakpoint_capacity) != 0) {
                ctx->stats->
                    critical_active_open_conic_breakpoint_fallbacks++;
                return 2;
            }
        }
        qsort(ctx->breakpoints, break_count, sizeof(double), double_compare);
        break_count = compact_sorted_breakpoints(ctx->breakpoints, break_count);
        int active_run = 0;
        double active_min = 0.0, active_max = 0.0;
        for (size_t i = 1; i < break_count; i++) {
            const double lo = ctx->breakpoints[i-1];
            const double hi = ctx->breakpoints[i];
            if (!(hi-lo > 1e-12)) continue;
            const double mid = 0.5*(lo+hi);
            double u, v, normal[2];
            int active = 0;
            if (open_conic_eval(&canonical, mid, branch, &u, &v) &&
                point_in_tile(ctx->tile, u, v, 1e-12*tile_scale) &&
                curve_normal_at(&canonical.curve,
                                (double[2]){u, v}, normal)) {
                if (ctx->stats->critical_active_boundary_tests >=
                    ctx->max_active_boundary_tests) return 2;
                double radius = fmax(
                    1e-6*tile_scale, 64.0*RAY_EPSILON);
                double ulo, vlo, uhi, vhi;
                if (open_conic_eval(
                        &canonical, lo, branch, &ulo, &vlo) &&
                    open_conic_eval(
                        &canonical, hi, branch, &uhi, &vhi)) {
                    const double piece_chord = hypot(uhi-ulo, vhi-vlo);
                    if (radius > 0.125*piece_chord)
                        radius = 0.125*piece_chord;
                }
                if (!(radius > 0.0)) continue;
                int contains[2];
                for (int side = 0; side < 2; side++) {
                    const double sign = side ? 1.0 : -1.0;
                    double local[3];
                    tile_world_point(ctx->view,
                        u+sign*radius*normal[0],
                        v+sign*radius*normal[1], local);
                    alea_matrix_transform_point_inverse(
                        &ctx->path->entries[ctx->level].transform,
                        &local[0], &local[1], &local[2]);
                    contains[side] = alea_contains_point(
                        ctx->sys, cell->root_node_id,
                        local[0], local[1], local[2]);
                }
                ctx->stats->critical_active_boundary_tests++;
                active = contains[0] != contains[1];
            }
            if (active) {
                if (!active_run) { active_min = lo; active_run = 1; }
                active_max = hi;
            } else if (active_run) {
                emit_active_open_conic_piece(
                    ctx, &canonical, active_min, active_max, branch,
                    occurrence_key, universe_occurrence_key);
                emitted = 1; active_run = 0;
            }
        }
        if (active_run) {
            emit_active_open_conic_piece(
                ctx, &canonical, active_min, active_max, branch,
                occurrence_key, universe_occurrence_key);
            emitted = 1;
        }
    }
    return emitted;
}

static int append_torus_break_value(double value, double* values,
                                    size_t* count, size_t capacity) {
    if (!isfinite(value)) return 0;
    if (*count == capacity) return -1;
    values[(*count)++] = value;
    return 0;
}

static int append_torus_curve_breaks(
    const alea_torus_2d_t* torus, const alea_curve_2d_t* other,
    double* values, size_t* count, size_t capacity) {
    double uv[16][2];
    int n = 0;
    if (curve_is_line(other->type)) {
        n = torus_line_intersections(torus, &other->data.line, uv);
        for (int i = 0; i < n; i++) {
            if (other->type != ALEA_CURVE_LINE) {
                const double du = uv[i][0]-other->data.line.point[0];
                const double dv = uv[i][1]-other->data.line.point[1];
                const double parameter =
                    du*other->data.line.direction[0] +
                    dv*other->data.line.direction[1];
                if (parameter < other->bounds.t_min-1e-12 ||
                    parameter > other->bounds.t_max+1e-12) continue;
            }
            if (append_torus_break_value(
                    uv[i][1], values, count, capacity) != 0) return -1;
        }
        return 0;
    }
    if (curve_is_closed_conic(other->type)) {
        n = torus_closed_conic_intersections(torus, other, uv);
    } else if (curve_is_general_conic(other->type)) {
        critical_curve_t open_item;
        memset(&open_item, 0, sizeof(open_item));
        open_item.curve = *other;
        n = torus_open_conic_intersections(torus, &open_item, uv);
    } else if (other->type == ALEA_CURVE_QUARTIC) {
        alea_curve_2d_t circles[2];
        const int circle_count = torus_circle_components(
            &other->data.torus, circles);
        if (circle_count == 0) {
            n = torus_torus_intersections(torus, &other->data.torus, uv);
        } else {
            for (int circle = 0; circle < circle_count; circle++) {
                double component_uv[8][2];
                const int component_count = torus_closed_conic_intersections(
                    torus, &circles[circle], component_uv);
                if (component_count < 0) return -1;
                for (int i = 0; i < component_count; i++)
                    if (append_torus_break_value(
                            component_uv[i][1], values, count,
                            capacity) != 0) return -1;
            }
            return 0;
        }
    } else if (other->type == ALEA_CURVE_POINT) {
        return 0;
    } else if (other->type == ALEA_CURVE_POLYGON) {
        const alea_polygon_2d_t* polygon = &other->data.polygon;
        const int edge_count = polygon->closed
            ? polygon->vertex_count : polygon->vertex_count-1;
        for (int i = 0; i < edge_count; i++) {
            const int j = (i+1)%polygon->vertex_count;
            alea_curve_2d_t edge;
            memset(&edge, 0, sizeof(edge));
            edge.type = ALEA_CURVE_LINE_SEGMENT;
            edge.data.line.point[0] = polygon->vertices[i][0];
            edge.data.line.point[1] = polygon->vertices[i][1];
            const double du = polygon->vertices[j][0]-polygon->vertices[i][0];
            const double dv = polygon->vertices[j][1]-polygon->vertices[i][1];
            const double length = hypot(du, dv);
            if (!(length > 0.0)) continue;
            edge.data.line.direction[0] = du/length;
            edge.data.line.direction[1] = dv/length;
            edge.bounds.t_min = 0.0; edge.bounds.t_max = length;
            if (append_torus_curve_breaks(
                    torus, &edge, values, count, capacity) != 0) return -1;
        }
        return 0;
    } else {
        return -1;
    }
    if (n < 0) return -1;
    for (int i = 0; i < n; i++)
        if (append_torus_break_value(
                uv[i][1], values, count, capacity) != 0) return -1;
    return 0;
}

static int track_scanline_root(const alea_curve_2d_t* curve,
                               double start_u, double start_v,
                               double end_v, double* end_u) {
    double current_u = start_u;
    for (int step = 1; step < 16; step++) {
        const double fraction = (double)step/16.0;
        const double v = start_v+fraction*(end_v-start_v);
        double roots[8];
        const int count = alea_curve_scanline_intersect(
            curve, v, roots, (int)(sizeof(roots)/sizeof(roots[0])));
        if (count <= 0) return 0;
        int closest = 0;
        for (int i = 1; i < count; i++)
            if (fabs(roots[i]-current_u) < fabs(roots[closest]-current_u))
                closest = i;
        current_u = roots[closest];
    }
    /* A repeated root at a vertical tangency can be lost at a different
     * distance from the exact endpoint on different libm implementations.
     * Approach geometrically from the last proven 15/16 sample and retain
     * the closest successfully solved inside limit.  Do not require one
     * arbitrarily chosen near-endpoint sample to succeed. */
    double roots[8];
    for (int exponent = 5; exponent <= 40; exponent++) {
        const double remaining = ldexp(1.0, -exponent);
        const double v = end_v+(start_v-end_v)*remaining;
        const int count = alea_curve_scanline_intersect(
            curve, v, roots, (int)(sizeof(roots)/sizeof(roots[0])));
        if (count <= 0) break;
        int closest = 0;
        for (int i = 1; i < count; i++)
            if (fabs(roots[i]-current_u) < fabs(roots[closest]-current_u))
                closest = i;
        current_u = roots[closest];
    }
    int count = alea_curve_scanline_intersect(
        curve, end_v, roots, (int)(sizeof(roots)/sizeof(roots[0])));
    if (count > 0) {
        int closest = 0;
        for (int i = 1; i < count; i++)
            if (fabs(roots[i]-current_u) < fabs(roots[closest]-current_u))
                closest = i;
        current_u = roots[closest];
    }
    *end_u = current_u;
    return 1;
}

static void emit_active_scanline_piece(
    critical_region_visit_t* ctx, const critical_curve_t* item,
    double v_min, double v_max, int root_index,
    double midpoint_u, double endpoint_min_u, double endpoint_max_u,
    uint64_t occurrence_key, uint64_t universe_occurrence_key) {
    critical_curve_t piece = *item;
    piece.has_scanline_domain = 1;
    piece.scanline_root = (uint8_t)root_index;
    piece.scanline_v_min = v_min;
    piece.scanline_v_max = v_max;
    piece.scanline_endpoint_u[0] = endpoint_min_u;
    piece.scanline_endpoint_u[1] = endpoint_max_u;
    piece.has_active_point = 1;
    piece.active_uv[0] = midpoint_u;
    piece.active_uv[1] = 0.5*(v_min+v_max);
    piece.bbox[0] = fmax(item->bbox[0], ctx->tile->uv_min[0]);
    piece.bbox[1] = fmin(item->bbox[1], ctx->tile->uv_max[0]);
    piece.bbox[2] = v_min;
    piece.bbox[3] = v_max;
    const double center_u = 0.5*(ctx->tile->uv_min[0]+ctx->tile->uv_max[0]);
    const double center_v = 0.5*(ctx->tile->uv_min[1]+ctx->tile->uv_max[1]);
    piece.priority = hypot(center_u-midpoint_u,
                           center_v-piece.active_uv[1]);
    retain_ranked_curve(ctx, &piece, occurrence_key,
                        universe_occurrence_key);
    ctx->stats->critical_active_segments++;
}

static int emit_cell_torus_pieces(
    critical_region_visit_t* ctx, const alea_cell_entry_t* cell,
    size_t curve_index, size_t cell_curve_count,
    uint64_t occurrence_key, uint64_t universe_occurrence_key) {
    const critical_curve_t* item = &ctx->cell_curves[curve_index];
    if (item->curve.type != ALEA_CURVE_QUARTIC) return 2;
    if (item->curve.data.torus.mode != 0) {
        const critical_curve_t original = *item;
        alea_curve_2d_t components[2];
        const int component_count = torus_circle_components(
            &original.curve.data.torus, components);
        if (component_count <= 0) return 2;
        int emitted = 0;
        for (int component = 0; component < component_count; component++) {
            critical_curve_t circle = original;
            const int surface_id = circle.curve.surface_id;
            const int primitive_id = circle.curve.primitive_id;
            const int universe_id = circle.curve.universe_id;
            circle.curve = components[component];
            circle.curve.surface_id = surface_id;
            circle.curve.primitive_id = primitive_id;
            circle.curve.universe_id = universe_id;
            circle.component_index = (uint8_t)(component+1);
            alea_curve_bbox(&circle.curve, &circle.bbox[0], &circle.bbox[1],
                            &circle.bbox[2], &circle.bbox[3]);
            ctx->cell_curves[curve_index] = circle;
            const int activity = emit_cell_circle_arcs(
                ctx, cell, curve_index, cell_curve_count,
                occurrence_key, universe_occurrence_key);
            if (activity == 2) {
                ctx->cell_curves[curve_index] = original;
                return 2;
            }
            if (activity == 1) emitted = 1;
        }
        ctx->cell_curves[curve_index] = original;
        return emitted;
    }
    size_t break_count = 2;
    ctx->breakpoints[0] = ctx->tile->uv_min[1];
    ctx->breakpoints[1] = ctx->tile->uv_max[1];
    double tangencies[16];
    const int tangent_count = torus_vertical_tangencies(
        &item->curve.data.torus, tangencies);
    if (tangent_count < 0) return 2;
    for (int i = 0; i < tangent_count; i++) {
        if (tangencies[i] <= ctx->tile->uv_min[1] ||
            tangencies[i] >= ctx->tile->uv_max[1]) continue;
        if (append_torus_break_value(
                tangencies[i], ctx->breakpoints, &break_count,
                ctx->breakpoint_capacity) != 0) return 2;
    }
    for (size_t i = 0; i < cell_curve_count; i++) {
        if (i == curve_index) continue;
        if (append_torus_curve_breaks(
                &item->curve.data.torus, &ctx->cell_curves[i].curve,
                ctx->breakpoints, &break_count,
                ctx->breakpoint_capacity) != 0) return 2;
    }
    qsort(ctx->breakpoints, break_count, sizeof(double), double_compare);
    break_count = compact_sorted_breakpoints(ctx->breakpoints, break_count);
    const double tile_scale = fmax(1.0, fmax(
        ctx->tile->uv_max[0]-ctx->tile->uv_min[0],
        ctx->tile->uv_max[1]-ctx->tile->uv_min[1]));
    int emitted = 0;
    for (size_t interval = 1; interval < break_count; interval++) {
        const double lo = ctx->breakpoints[interval-1];
        const double hi = ctx->breakpoints[interval];
        if (!(hi-lo > 1e-12*tile_scale)) continue;
        const double mid = 0.5*(lo+hi);
        double roots[8];
        const int root_count = alea_curve_scanline_intersect(
            &item->curve, mid, roots,
            (int)(sizeof(roots)/sizeof(roots[0])));
        sort_small_doubles(roots, root_count);
        for (int root = 0; root < root_count; root++) {
            const double u = roots[root];
            if (!point_in_tile(ctx->tile, u, mid, 1e-12*tile_scale)) continue;
            double normal[2];
            if (!curve_normal_at(
                    &item->curve, (double[2]){u, mid}, normal)) return 2;
            if (ctx->stats->critical_active_boundary_tests >=
                ctx->max_active_boundary_tests) return 2;
            double radius = fmax(1e-6*tile_scale, 64.0*RAY_EPSILON);
            if (radius > 0.125*(hi-lo)) radius = 0.125*(hi-lo);
            if (!(radius > 0.0)) continue;
            int contains[2];
            for (int side = 0; side < 2; side++) {
                const double sign = side ? 1.0 : -1.0;
                double local[3];
                tile_world_point(ctx->view, u+sign*radius*normal[0],
                                 mid+sign*radius*normal[1], local);
                alea_matrix_transform_point_inverse(
                    &ctx->path->entries[ctx->level].transform,
                    &local[0], &local[1], &local[2]);
                contains[side] = alea_contains_point(
                    ctx->sys, cell->root_node_id,
                    local[0], local[1], local[2]);
            }
            ctx->stats->critical_active_boundary_tests++;
            if (contains[0] == contains[1]) continue;
            double endpoint_u[2];
            if (!track_scanline_root(
                    &item->curve, u, mid, lo, &endpoint_u[0]) ||
                !track_scanline_root(
                    &item->curve, u, mid, hi, &endpoint_u[1])) return 2;
            emit_active_scanline_piece(
                ctx, item, lo, hi, root, u,
                endpoint_u[0], endpoint_u[1], occurrence_key,
                universe_occurrence_key);
            emitted = 1;
        }
    }
    return emitted;
}

static int append_implicit_conic_breaks(
    const alea_curve_2d_t* source, const alea_curve_2d_t* other,
    double* values, size_t* count, size_t capacity) {
    const alea_conic_2d_t* q = &source->data.conic;
    double uv[16][2];
    int n = 0;
    if (curve_is_line(other->type)) {
        const alea_line_2d_t* line = &other->data.line;
        const double u = line->point[0], v = line->point[1];
        const double du = line->direction[0], dv = line->direction[1];
        double roots[2];
        n = alea_solve_quadratic(
            q->A*du*du + q->B*du*dv + q->C*dv*dv,
            2.0*q->A*u*du + q->B*(u*dv+v*du) +
                2.0*q->C*v*dv + q->D*du + q->E*dv,
            conic_value(q, u, v), roots);
        for (int i = 0; i < n; i++) {
            const double parameter = roots[i];
            if (other->type != ALEA_CURVE_LINE &&
                (parameter < other->bounds.t_min-1e-12 ||
                 parameter > other->bounds.t_max+1e-12)) continue;
            if (append_torus_break_value(
                    v+parameter*dv, values, count, capacity) != 0) return -1;
        }
        return 0;
    }
    if (curve_is_closed_conic(other->type)) {
        n = parametric_closed_conic_intersections(other, q, uv);
    } else if (curve_is_general_conic(other->type)) {
        n = general_conic_intersections(q, &other->data.conic, uv);
    } else if (other->type == ALEA_CURVE_POINT) {
        const double u = other->data.point[0], v = other->data.point[1];
        const double scale = fabs(q->A*u*u)+fabs(q->B*u*v)+
            fabs(q->C*v*v)+fabs(q->D*u)+fabs(q->E*v)+fabs(q->F)+1.0;
        if (fabs(conic_value(q, u, v)) <= 1e-10*scale)
            return append_torus_break_value(v, values, count, capacity);
        return 0;
    } else if (other->type == ALEA_CURVE_POLYGON) {
        const alea_polygon_2d_t* polygon = &other->data.polygon;
        const int edge_count = polygon->closed
            ? polygon->vertex_count : polygon->vertex_count-1;
        for (int edge_index = 0; edge_index < edge_count; edge_index++) {
            const int next = (edge_index+1)%polygon->vertex_count;
            const double du = polygon->vertices[next][0]-
                              polygon->vertices[edge_index][0];
            const double dv = polygon->vertices[next][1]-
                              polygon->vertices[edge_index][1];
            const double length = hypot(du, dv);
            if (!(length > 0.0)) continue;
            alea_curve_2d_t edge;
            memset(&edge, 0, sizeof(edge));
            edge.type = ALEA_CURVE_LINE_SEGMENT;
            edge.data.line.point[0] = polygon->vertices[edge_index][0];
            edge.data.line.point[1] = polygon->vertices[edge_index][1];
            edge.data.line.direction[0] = du/length;
            edge.data.line.direction[1] = dv/length;
            edge.bounds.t_min = 0.0; edge.bounds.t_max = length;
            if (append_implicit_conic_breaks(
                    source, &edge, values, count, capacity) != 0) return -1;
        }
        return 0;
    } else if (other->type == ALEA_CURVE_QUARTIC) {
        alea_curve_2d_t circles[2];
        const int circle_count = torus_circle_components(
            &other->data.torus, circles);
        if (circle_count > 0) {
            for (int circle = 0; circle < circle_count; circle++) {
                double circle_uv[4][2];
                const int circle_points =
                    parametric_closed_conic_intersections(
                        &circles[circle], q, circle_uv);
                if (circle_points < 0) return -1;
                for (int i = 0; i < circle_points; i++)
                    if (append_torus_break_value(
                            circle_uv[i][1], values, count,
                            capacity) != 0) return -1;
            }
            return 0;
        }
        n = torus_implicit_conic_intersections(
            &other->data.torus, q, uv);
    } else {
        return -1;
    }
    if (n < 0) return -1;
    for (int i = 0; i < n; i++)
        if (append_torus_break_value(
                uv[i][1], values, count, capacity) != 0) return -1;
    return 0;
}

static int emit_cell_implicit_conic_pieces(
    critical_region_visit_t* ctx, const alea_cell_entry_t* cell,
    size_t curve_index, size_t cell_curve_count,
    uint64_t occurrence_key, uint64_t universe_occurrence_key) {
    const critical_curve_t* item = &ctx->cell_curves[curve_index];
    if (!curve_is_general_conic(item->curve.type)) return 2;
    const alea_conic_2d_t* q = &item->curve.data.conic;
    size_t break_count = 2;
    ctx->breakpoints[0] = ctx->tile->uv_min[1];
    ctx->breakpoints[1] = ctx->tile->uv_max[1];
    double roots[2];
    int root_count = alea_solve_quadratic(
        q->B*q->B-4.0*q->A*q->C,
        2.0*q->B*q->D-4.0*q->A*q->E,
        q->D*q->D-4.0*q->A*q->F, roots);
    for (int i = 0; i < root_count; i++)
        if (roots[i] > ctx->tile->uv_min[1] &&
            roots[i] < ctx->tile->uv_max[1] &&
            append_torus_break_value(
                roots[i], ctx->breakpoints, &break_count,
                ctx->breakpoint_capacity) != 0) return 2;
    const double quadratic_scale = fmax(
        fabs(q->A)+fabs(q->B)+fabs(q->C), DBL_MIN);
    if (fabs(q->A) <= 1e-12*quadratic_scale && fabs(q->B) > DBL_MIN) {
        const double degree_loss = -q->D/q->B;
        if (degree_loss > ctx->tile->uv_min[1] &&
            degree_loss < ctx->tile->uv_max[1] &&
            append_torus_break_value(
                degree_loss, ctx->breakpoints, &break_count,
                ctx->breakpoint_capacity) != 0) return 2;
    }
    alea_curve_2d_t edge;
    memset(&edge, 0, sizeof(edge));
    edge.type = ALEA_CURVE_LINE;
    for (int side = 0; side < 2; side++) {
        const double u = side
            ? ctx->tile->uv_max[0] : ctx->tile->uv_min[0];
        edge.data.line.point[0] = u;
        edge.data.line.point[1] = 0.0;
        edge.data.line.direction[0] = 0.0;
        edge.data.line.direction[1] = 1.0;
        if (append_implicit_conic_breaks(
                &item->curve, &edge, ctx->breakpoints, &break_count,
                ctx->breakpoint_capacity) != 0) return 2;
    }
    for (size_t i = 0; i < cell_curve_count; i++) {
        if (i == curve_index) continue;
        if (append_implicit_conic_breaks(
                &item->curve, &ctx->cell_curves[i].curve,
                ctx->breakpoints, &break_count,
                ctx->breakpoint_capacity) != 0) return 2;
    }
    qsort(ctx->breakpoints, break_count, sizeof(double), double_compare);
    break_count = compact_sorted_breakpoints(ctx->breakpoints, break_count);
    const double tile_scale = fmax(1.0, fmax(
        ctx->tile->uv_max[0]-ctx->tile->uv_min[0],
        ctx->tile->uv_max[1]-ctx->tile->uv_min[1]));
    int emitted = 0;
    for (size_t interval = 1; interval < break_count; interval++) {
        const double lo = ctx->breakpoints[interval-1];
        const double hi = ctx->breakpoints[interval];
        if (!(hi-lo > 1e-12*tile_scale)) continue;
        const double mid = 0.5*(lo+hi);
        double scanline_roots[8];
        const int scanline_count = alea_curve_scanline_intersect(
            &item->curve, mid, scanline_roots,
            (int)(sizeof(scanline_roots)/sizeof(scanline_roots[0])));
        sort_small_doubles(scanline_roots, scanline_count);
        for (int root = 0; root < scanline_count; root++) {
            const double u = scanline_roots[root];
            if (!point_in_tile(ctx->tile, u, mid, 1e-12*tile_scale)) continue;
            double normal[2];
            if (!curve_normal_at(
                    &item->curve, (double[2]){u, mid}, normal)) return 2;
            if (ctx->stats->critical_active_boundary_tests >=
                ctx->max_active_boundary_tests) return 2;
            double radius = fmax(1e-6*tile_scale, 64.0*RAY_EPSILON);
            if (radius > 0.125*(hi-lo)) radius = 0.125*(hi-lo);
            if (!(radius > 0.0)) continue;
            int contains[2];
            for (int side = 0; side < 2; side++) {
                const double sign = side ? 1.0 : -1.0;
                double local[3];
                tile_world_point(ctx->view, u+sign*radius*normal[0],
                                 mid+sign*radius*normal[1], local);
                alea_matrix_transform_point_inverse(
                    &ctx->path->entries[ctx->level].transform,
                    &local[0], &local[1], &local[2]);
                contains[side] = alea_contains_point(
                    ctx->sys, cell->root_node_id,
                    local[0], local[1], local[2]);
            }
            ctx->stats->critical_active_boundary_tests++;
            if (contains[0] == contains[1]) continue;
            double endpoint_u[2];
            if (!track_scanline_root(
                    &item->curve, u, mid, lo, &endpoint_u[0]) ||
                !track_scanline_root(
                    &item->curve, u, mid, hi, &endpoint_u[1])) return 2;
            emit_active_scanline_piece(
                ctx, item, lo, hi, root, u,
                endpoint_u[0], endpoint_u[1], occurrence_key,
                universe_occurrence_key);
            emitted = 1;
        }
    }
    return emitted;
}

/* Return 0 when inactive, 1 after publishing all proven active segments, and
 * 2 when the caller must conservatively retain the whole analytical curve. */
static int emit_cell_line_segments(
    critical_region_visit_t* ctx, const alea_cell_entry_t* cell,
    size_t curve_index, size_t cell_curve_count,
    uint64_t occurrence_key, uint64_t universe_occurrence_key) {
    const critical_curve_t* item = &ctx->cell_curves[curve_index];
    const alea_curve_2d_t* curve = &item->curve;
    int emitted = 0;
    if (!curve_is_line(curve->type)) return 2;
    double lo, hi;
    if (!clip_line_to_tile(curve, ctx->tile, &lo, &hi)) return 2;
    size_t break_count = 2;
    ctx->breakpoints[0] = lo;
    ctx->breakpoints[1] = hi;
    for (size_t i = 0; i < cell_curve_count; i++) {
        if (i == curve_index) continue;
        if (append_line_curve_breaks(
                curve, &ctx->cell_curves[i].curve, lo, hi,
                ctx->breakpoints, &break_count,
                ctx->breakpoint_capacity) != 0)
            return 2;
    }
    qsort(ctx->breakpoints, break_count, sizeof(double), double_compare);
    break_count = compact_sorted_breakpoints(ctx->breakpoints, break_count);
    const alea_line_2d_t* line = &curve->data.line;
    const double normal[2] = {-line->direction[1], line->direction[0]};
    const double tile_scale = fmax(1.0, fmax(
        ctx->tile->uv_max[0] - ctx->tile->uv_min[0],
        ctx->tile->uv_max[1] - ctx->tile->uv_min[1]));
    int active_run = 0;
    double active_min = 0.0, active_max = 0.0;
    for (size_t i = 1; i < break_count; i++) {
        const double span = ctx->breakpoints[i] - ctx->breakpoints[i - 1];
        if (!(span > 1e-12 * tile_scale)) continue;
        const double mid = 0.5 * (ctx->breakpoints[i - 1] +
                                  ctx->breakpoints[i]);
        double radius = 1e-6 * tile_scale;
        if (ctx->stats->critical_active_boundary_tests >=
            ctx->max_active_boundary_tests) {
            return 2;
        }
        if (!(radius > 64.0 * RAY_EPSILON)) radius = 64.0 * RAY_EPSILON;
        if (radius > 0.125 * span) radius = 0.125 * span;
        int contains[2];
        for (int side = 0; side < 2; side++) {
            const double sign = side ? 1.0 : -1.0;
            const double u = line->point[0] + mid*line->direction[0] +
                             sign*radius*normal[0];
            const double v = line->point[1] + mid*line->direction[1] +
                             sign*radius*normal[1];
            double local[3];
            tile_world_point(ctx->view, u, v, local);
            alea_matrix_transform_point_inverse(
                &ctx->path->entries[ctx->level].transform,
                &local[0], &local[1], &local[2]);
            contains[side] = alea_contains_point(
                ctx->sys, cell->root_node_id, local[0], local[1], local[2]);
        }
        ctx->stats->critical_active_boundary_tests++;
        if (contains[0] != contains[1]) {
            if (!active_run) {
                active_min = ctx->breakpoints[i - 1];
                active_run = 1;
            }
            active_max = ctx->breakpoints[i];
        } else if (active_run) {
            emit_active_line_segment(
                ctx, item, line, active_min, active_max,
                occurrence_key, universe_occurrence_key);
            emitted = 1;
            active_run = 0;
        }
    }
    if (active_run) {
        emit_active_line_segment(
            ctx, item, line, active_min, active_max,
            occurrence_key, universe_occurrence_key);
        emitted = 1;
    }
    return emitted;
}

/* Factor a real degenerate conic into one or two infinite lines.  Return the
 * line count, zero for a proven empty real locus, or -1 when the conic is not
 * a supported degenerate line family. */
static int factor_degenerate_conic_lines(
    const alea_curve_2d_t* curve, alea_curve_2d_t lines[2]) {
    const alea_conic_2d_t* q = &curve->data.conic;
    const double angle = 0.5*atan2(q->B, q->A-q->C);
    const double c = cos(angle), s = sin(angle);
    const double eigenvalues[2] = {
        q->A*c*c + q->B*c*s + q->C*s*s,
        q->A*s*s - q->B*c*s + q->C*c*c};
    const double axes[2][2] = {{c, s}, {-s, c}};
    const double linear[2] = {q->D*c+q->E*s, -q->D*s+q->E*c};
    const double quadratic_scale = fmax(
        fabs(q->A)+fabs(q->B)+fabs(q->C), DBL_MIN);
    const double linear_scale = fmax(fabs(q->D)+fabs(q->E), DBL_MIN);
    const double quadratic_tolerance = 1e-12*quadratic_scale;
    const double axial_linear_tolerance =
        64.0*DBL_EPSILON*linear_scale;
    const long double determinant =
        4.0L*q->A*q->C-(long double)q->B*q->B;
    memset(lines, 0, 2*sizeof(*lines));

    if (fabsl(determinant) >
        1e-12*quadratic_scale*quadratic_scale) {
        const long double center[2] = {
            ((long double)q->B*q->E-2.0L*q->C*q->D)/determinant,
            ((long double)q->B*q->D-2.0L*q->A*q->E)/determinant};
        const long double terms[6] = {
            (long double)q->A*center[0]*center[0],
            (long double)q->B*center[0]*center[1],
            (long double)q->C*center[1]*center[1],
            (long double)q->D*center[0],
            (long double)q->E*center[1], q->F};
        long double constant = 0.0L, scale = LDBL_MIN;
        for (int i = 0; i < 6; i++) {
            constant += terms[i];
            scale += fabsl(terms[i]);
        }
        if (fabsl(constant) > 1e-10L*scale ||
            eigenvalues[0]*eigenvalues[1] >= 0.0) return -1;
        const double slope = sqrt(-eigenvalues[0]/eigenvalues[1]);
        for (int i = 0; i < 2; i++) {
            const double sign = i ? 1.0 : -1.0;
            double du = axes[0][0]+sign*slope*axes[1][0];
            double dv = axes[0][1]+sign*slope*axes[1][1];
            const double length = hypot(du, dv);
            if (!(length > 0.0)) return -1;
            du /= length; dv /= length;
            lines[i].type = ALEA_CURVE_LINE;
            lines[i].data.line.point[0] = (double)center[0];
            lines[i].data.line.point[1] = (double)center[1];
            lines[i].data.line.direction[0] = du;
            lines[i].data.line.direction[1] = dv;
            lines[i].data.line.a = -dv;
            lines[i].data.line.b = du;
            lines[i].data.line.c = dv*(double)center[0]-
                                    du*(double)center[1];
        }
        return 2;
    }

    const int quadratic = fabs(eigenvalues[1]) > fabs(eigenvalues[0]);
    const int axial = 1-quadratic;
    const double lambda = eigenvalues[quadratic];
    if (fabs(lambda) <= quadratic_tolerance ||
        fabs(linear[axial]) > axial_linear_tolerance) return -1;
    const double x0 = -linear[quadratic]/(2.0*lambda);
    const long double completed = (long double)q->F -
        (long double)linear[quadratic]*linear[quadratic]/
            (4.0L*lambda);
    const long double completed_scale = fabsl(q->F) +
        fabsl((long double)linear[quadratic]*linear[quadratic]/
              (4.0L*lambda)) + LDBL_MIN;
    long double squared_offset = -completed/lambda;
    if (squared_offset < -1e-10L*completed_scale/fabs(lambda)) return 0;
    int count = 2;
    if (squared_offset <= 1e-10L*completed_scale/fabs(lambda)) {
        squared_offset = 0.0L;
        count = 1;
    }
    const double offset = sqrt((double)squared_offset);
    for (int i = 0; i < count; i++) {
        const double coordinate = x0+(i ? offset : -offset);
        lines[i].type = ALEA_CURVE_LINE;
        lines[i].data.line.point[0] = coordinate*axes[quadratic][0];
        lines[i].data.line.point[1] = coordinate*axes[quadratic][1];
        lines[i].data.line.direction[0] = axes[axial][0];
        lines[i].data.line.direction[1] = axes[axial][1];
        lines[i].data.line.a = -axes[axial][1];
        lines[i].data.line.b = axes[axial][0];
        lines[i].data.line.c =
            -(lines[i].data.line.a*lines[i].data.line.point[0] +
              lines[i].data.line.b*lines[i].data.line.point[1]);
    }
    return count;
}

static int emit_cell_degenerate_conic_lines(
    critical_region_visit_t* ctx, const alea_cell_entry_t* cell,
    size_t curve_index, size_t cell_curve_count,
    uint64_t occurrence_key, uint64_t universe_occurrence_key) {
    const critical_curve_t original = ctx->cell_curves[curve_index];
    alea_curve_2d_t lines[2];
    const int line_count = factor_degenerate_conic_lines(
        &original.curve, lines);
    if (line_count < 0) return 2;
    int emitted = 0;
    for (int line_index = 0; line_index < line_count; line_index++) {
        critical_curve_t line = original;
        const int surface_id = line.curve.surface_id;
        const int primitive_id = line.curve.primitive_id;
        const int universe_id = line.curve.universe_id;
        line.curve = lines[line_index];
        line.curve.surface_id = surface_id;
        line.curve.primitive_id = primitive_id;
        line.curve.universe_id = universe_id;
        line.component_index = (uint8_t)(line_index+1);
        alea_curve_bbox(&line.curve, &line.bbox[0], &line.bbox[1],
                        &line.bbox[2], &line.bbox[3]);
        double clipped_min, clipped_max;
        if (!clip_line_to_tile(
                &line.curve, ctx->tile, &clipped_min, &clipped_max))
            continue;
        ctx->cell_curves[curve_index] = line;
        const int activity = emit_cell_line_segments(
            ctx, cell, curve_index, cell_curve_count,
            occurrence_key, universe_occurrence_key);
        if (activity == 2) {
            ctx->cell_curves[curve_index] = original;
            return 2;
        }
        if (activity == 1) emitted = 1;
    }
    ctx->cell_curves[curve_index] = original;
    return emitted;
}

static int emit_cell_polygon_edges(
    critical_region_visit_t* ctx, const alea_cell_entry_t* cell,
    size_t curve_index, size_t cell_curve_count,
    uint64_t occurrence_key, uint64_t universe_occurrence_key) {
    const critical_curve_t original = ctx->cell_curves[curve_index];
    if (original.curve.type != ALEA_CURVE_POLYGON) return 2;
    const alea_polygon_2d_t* polygon = &original.curve.data.polygon;
    const int edge_count = polygon->closed
        ? polygon->vertex_count : polygon->vertex_count-1;
    if (edge_count <= 0 || edge_count > ALEA_MAX_POLYGON_VERTICES) return 2;
    int emitted = 0;
    for (int edge_index = 0; edge_index < edge_count; edge_index++) {
        const int next = (edge_index+1)%polygon->vertex_count;
        const double du = polygon->vertices[next][0] -
                          polygon->vertices[edge_index][0];
        const double dv = polygon->vertices[next][1] -
                          polygon->vertices[edge_index][1];
        const double length = hypot(du, dv);
        if (!(length > 0.0)) continue;
        critical_curve_t edge = original;
        edge.curve.type = ALEA_CURVE_LINE_SEGMENT;
        edge.curve.data.line.point[0] = polygon->vertices[edge_index][0];
        edge.curve.data.line.point[1] = polygon->vertices[edge_index][1];
        edge.curve.data.line.direction[0] = du/length;
        edge.curve.data.line.direction[1] = dv/length;
        edge.curve.bounds.t_min = 0.0;
        edge.curve.bounds.t_max = length;
        edge.component_index = (uint8_t)(edge_index+1);
        alea_curve_bbox(&edge.curve, &edge.bbox[0], &edge.bbox[1],
                        &edge.bbox[2], &edge.bbox[3]);
        double clipped_min, clipped_max;
        if (!clip_line_to_tile(
                &edge.curve, ctx->tile, &clipped_min, &clipped_max))
            continue;
        ctx->cell_curves[curve_index] = edge;
        const int activity = emit_cell_line_segments(
            ctx, cell, curve_index, cell_curve_count,
            occurrence_key, universe_occurrence_key);
        if (activity == 2) {
            ctx->cell_curves[curve_index] = original;
            return 2;
        }
        if (activity == 1) emitted = 1;
    }
    ctx->cell_curves[curve_index] = original;
    return emitted;
}

static int retain_ranked_curve(critical_region_visit_t* ctx,
                               const critical_curve_t* item,
                               uint64_t occurrence_key,
                               uint64_t universe_occurrence_key) {
    for (size_t i = 0; i < *ctx->curve_count; i++) {
        const critical_curve_t* saved = &ctx->curves[i];
        if (saved->surface_index == item->surface_index &&
            saved->universe_occurrence_key == universe_occurrence_key &&
            saved->curve.type == item->curve.type &&
            saved->component_index == item->component_index &&
            saved->has_parameter_domain == item->has_parameter_domain &&
            saved->has_scanline_domain == item->has_scanline_domain &&
            (!item->has_parameter_domain ||
             (saved->conic_branch == item->conic_branch &&
              fabs(saved->parameter_min-item->parameter_min) <= 1e-12 &&
              fabs(saved->parameter_max-item->parameter_max) <= 1e-12)) &&
            (!item->has_scanline_domain ||
             (saved->scanline_root == item->scanline_root &&
              fabs(saved->scanline_v_min-item->scanline_v_min) <= 1e-12 &&
              fabs(saved->scanline_v_max-item->scanline_v_max) <= 1e-12)) &&
            ((item->curve.type != ALEA_CURVE_LINE_SEGMENT &&
              item->curve.type != ALEA_CURVE_ARC &&
              item->curve.type != ALEA_CURVE_ELLIPSE_ARC) ||
             (fabs(saved->curve.bounds.t_min - item->curve.bounds.t_min) <=
                  1e-12 &&
              fabs(saved->curve.bounds.t_max - item->curve.bounds.t_max) <=
                  1e-12 &&
              fabs(saved->active_uv[0] - item->active_uv[0]) <= 1e-12 &&
              fabs(saved->active_uv[1] - item->active_uv[1]) <= 1e-12))) {
            ctx->stats->critical_duplicate_surface_occurrences++;
            return 0;
        }
    }
    critical_curve_t ranked = *item;
    ranked.occurrence_key = occurrence_key;
    ranked.universe_occurrence_key = universe_occurrence_key;
    if (*ctx->curve_count < ctx->curve_capacity) {
        ctx->curves[(*ctx->curve_count)++] = ranked;
        return 0;
    }
    size_t worst = 0;
    for (size_t i = 1; i < *ctx->curve_count; i++) {
        const critical_curve_t* a = &ctx->curves[i];
        const critical_curve_t* b = &ctx->curves[worst];
        if (a->priority_class < b->priority_class ||
            (a->priority_class == b->priority_class &&
             (a->priority > b->priority ||
              (a->priority == b->priority &&
               a->surface_index > b->surface_index))))
            worst = i;
    }
    const critical_curve_t* old = &ctx->curves[worst];
    if (ranked.priority_class > old->priority_class ||
        (ranked.priority_class == old->priority_class &&
         (ranked.priority < old->priority ||
          (ranked.priority == old->priority &&
           ranked.surface_index < old->surface_index))))
        ctx->curves[worst] = ranked;
    ctx->stats->critical_ranked_curves_omitted++;
    return 0;
}

static int emit_cell_point(
    critical_region_visit_t* ctx, const critical_curve_t* item,
    uint64_t occurrence_key, uint64_t universe_occurrence_key) {
    if (item->curve.type != ALEA_CURVE_POINT) return 2;
    if (!point_in_tile(ctx->tile, item->curve.data.point[0],
                       item->curve.data.point[1], 1e-12)) return 0;
    retain_ranked_curve(ctx, item, occurrence_key, universe_occurrence_key);
    return 1;
}

static int cell_has_both_sense_card(const alea_system_t* sys,
                                    const alea_cell_entry_t* cell,
                                    uint32_t cell_index) {
    if (!sys->surface_cell_offsets || !sys->surface_cell_refs) return 0;
    for (size_t si = 0; si < cell->surface_index_count; si++) {
        const uint32_t surface_index = cell->surface_indices[si];
        if (surface_index >= alea_vec_count(&sys->surfaces)) continue;
        int sense = 0;
        const size_t begin = sys->surface_cell_offsets[surface_index];
        const size_t end = sys->surface_cell_offsets[surface_index + 1];
        for (size_t ri = begin; ri < end; ri++) {
            const alea_surface_cell_ref_t* ref = &sys->surface_cell_refs[ri];
            if (ref->cell_index != cell_index) continue;
            if (sense && sense != ref->sense) return 1;
            sense = ref->sense;
        }
    }
    return 0;
}

static int collect_cell_curves(critical_region_visit_t* ctx,
                               const alea_spatial_hit_t* candidate,
                               int* out_relevant) {
    *out_relevant = 0;
    if (candidate->cell_index >= alea_vec_count(&ctx->sys->cells)) return 0;
    alea_hier_spatial_chain_hit_t chain;
    occurrence_chain_hit(ctx->path, ctx->level, candidate, &chain);
    if (chain.chain_truncated) return CRITICAL_COLLECT_CHAIN_TRUNCATED;
    uint64_t occurrence_key, universe_occurrence_key;
    chain_occurrences(&chain, &occurrence_key, &universe_occurrence_key);
    const alea_cell_entry_t* cell =
        &ctx->sys->cells.data[candidate->cell_index];
    const int both_sense_cell = cell_has_both_sense_card(
        ctx->sys, cell, candidate->cell_index);
    size_t cell_curve_count = 0;
    int conservative = 0;
    const size_t inspect_count = cell->surface_index_count;
    for (size_t si = 0; si < inspect_count; si++) {
        ctx->stats->critical_surface_references++;
        critical_curve_t item;
        if (!make_cell_curve(
                ctx, cell, cell->surface_indices[si], &item))
            continue;
        critical_curve_t items[2];
        int item_count = 1;
        items[0] = item;
        if (item.curve.type == ALEA_CURVE_PARALLEL_LINES)
            item_count = split_parallel_lines(ctx, &item, items);
        else if (!both_sense_cell && item.curve.type == ALEA_CURVE_QUARTIC) {
            const int special_count = split_special_torus(ctx, &item, items);
            if (special_count > 0) item_count = special_count;
        }
        for (int ii = 0; ii < item_count; ii++) {
            *out_relevant = 1;
            items[ii].priority_class = both_sense_cell ? 1 : 0;
            if (!conservative &&
                cell_curve_count == ctx->cell_curve_capacity) {
                conservative = 1;
                ctx->stats->critical_active_boundary_fallbacks++;
                ctx->stats->critical_active_capacity_fallbacks++;
                for (size_t saved = 0; saved < cell_curve_count; saved++) {
                    ctx->stats->critical_whole_curve_fallbacks++;
                    retain_ranked_curve(
                        ctx, &ctx->cell_curves[saved], occurrence_key,
                        universe_occurrence_key);
                }
            }
            if (!conservative) {
                ctx->cell_curves[cell_curve_count++] = items[ii];
                continue;
            }
            ctx->stats->critical_whole_curve_fallbacks++;
            retain_ranked_curve(ctx, &items[ii], occurrence_key,
                                universe_occurrence_key);
        }
    }
    if (conservative) return CRITICAL_COLLECT_OK;

    for (size_t ci = 0; ci < cell_curve_count; ci++) {
        critical_curve_t* item = &ctx->cell_curves[ci];
        alea_curve_type_t unsupported_type = ALEA_CURVE_NONE;
        if (!active_partition_pair_supported(
                item->curve.type, ALEA_CURVE_POINT)) {
            unsupported_type = item->curve.type;
        } else {
            for (size_t oi = 0; oi < cell_curve_count; oi++) {
                if (oi == ci) continue;
                if (!active_partition_pair_supported(
                        item->curve.type, ctx->cell_curves[oi].curve.type)) {
                    unsupported_type = ctx->cell_curves[oi].curve.type;
                    break;
                }
            }
        }
        const int activity = unsupported_type != ALEA_CURVE_NONE ? 2 :
            item->curve.type == ALEA_CURVE_CIRCLE
            ? emit_cell_circle_arcs(
                ctx, cell, ci, cell_curve_count,
                occurrence_key, universe_occurrence_key)
            : item->curve.type == ALEA_CURVE_ELLIPSE
            ? emit_cell_ellipse_arcs(
                ctx, cell, ci, cell_curve_count,
                occurrence_key, universe_occurrence_key)
            : curve_is_general_conic(item->curve.type)
            ? emit_cell_open_conic_pieces(
                ctx, cell, ci, cell_curve_count,
                occurrence_key, universe_occurrence_key)
            : item->curve.type == ALEA_CURVE_QUARTIC
            ? emit_cell_torus_pieces(
                ctx, cell, ci, cell_curve_count,
                occurrence_key, universe_occurrence_key)
            : item->curve.type == ALEA_CURVE_POLYGON
            ? emit_cell_polygon_edges(
                ctx, cell, ci, cell_curve_count,
                occurrence_key, universe_occurrence_key)
            : item->curve.type == ALEA_CURVE_POINT
            ? emit_cell_point(
                ctx, item, occurrence_key, universe_occurrence_key)
            : emit_cell_line_segments(
                ctx, cell, ci, cell_curve_count,
                occurrence_key, universe_occurrence_key);
        if (activity == 1) continue;
        if (activity == 0) {
            ctx->stats->critical_curves_culled++;
            continue;
        }
        ctx->stats->critical_active_boundary_fallbacks++;
        if (unsupported_type != ALEA_CURVE_NONE)
            record_active_unsupported(unsupported_type, ctx->stats);
        else if (ctx->stats->critical_active_boundary_tests >=
                 ctx->max_active_boundary_tests)
            ctx->stats->critical_active_test_budget_fallbacks++;
        else {
            ctx->stats->critical_active_unsupported_other_fallbacks++;
            record_active_evaluation_failure(
                item->curve.type, ctx->stats);
        }
        ctx->stats->critical_whole_curve_fallbacks++;
        retain_ranked_curve(ctx, item, occurrence_key,
                            universe_occurrence_key);
    }
    return CRITICAL_COLLECT_OK;
}

static int visit_region_cell_curves(const alea_spatial_hit_t* candidate,
                                    void* userdata) {
    critical_region_visit_t* ctx = userdata;
    int relevant = 0;
    ctx->stop_code = collect_cell_curves(ctx, candidate, &relevant);
    if (relevant) ctx->stats->critical_region_hits++;
    if (ctx->stop_code != CRITICAL_COLLECT_OK)
        return 1;
    return 0;
}

static int collect_path_cell(critical_region_visit_t* ctx,
                             uint32_t cell_index) {
    if (cell_index >= alea_vec_count(&ctx->sys->cells)) return 0;
    const alea_cell_entry_t* cell = &ctx->sys->cells.data[cell_index];
    if (cell->universe_id != ctx->path->entries[ctx->level].universe_id)
        return 0;
    alea_spatial_hit_t candidate;
    memset(&candidate, 0, sizeof(candidate));
    candidate.cell_index = cell_index;
    candidate.cell_id = cell->mc_cell_id;
    candidate.material_id = cell->material_id;
    candidate.universe_id = cell->universe_id;
    candidate.is_terminal = !alea_cell_entry_is_container(cell);
    int relevant = 0;
    const int rc = collect_cell_curves(ctx, &candidate, &relevant);
    if (relevant) ctx->stats->critical_region_hits++;
    return rc;
}

static int collect_ancestor_neighborhood(critical_region_visit_t* ctx) {
    const uint32_t active_index =
        ctx->path->entries[ctx->level].cell_index;
    int rc = collect_path_cell(ctx, active_index);
    if (rc != 0 || active_index >= alea_vec_count(&ctx->sys->cells) ||
        !ctx->sys->surface_cell_offsets || !ctx->sys->surface_cell_refs)
        return rc;
    const alea_cell_entry_t* active = &ctx->sys->cells.data[active_index];
    for (size_t si = 0; si < active->surface_index_count; si++) {
        const uint32_t surface_index = active->surface_indices[si];
        if (surface_index >= alea_vec_count(&ctx->sys->surfaces)) continue;
        const size_t begin = ctx->sys->surface_cell_offsets[surface_index];
        const size_t end = ctx->sys->surface_cell_offsets[surface_index + 1];
        for (size_t ri = begin; ri < end; ri++) {
            rc = collect_path_cell(
                ctx, ctx->sys->surface_cell_refs[ri].cell_index);
            if (rc != 0) return rc;
        }
    }
    return 0;
}

/* Enumerate only universe occurrences reached by a fixed 3x3 set of points
 * in the suspicious tile.  The deepest occurrence gets a local BLAS query so
 * card-distinct overlap/gap partners are retained.  Ancestors use the active
 * path-cell neighborhood so broad root bounds cannot starve the leaf query. */
static int collect_tile_occurrence_curves(
    alea_system_t* sys, const alea_slice_view_t* view,
    const alea_transition_slice_critical_tile_t* tile,
    critical_universe_occurrence_t* occurrences, size_t occurrence_capacity,
    critical_curve_t* curves, size_t* curve_count, size_t curve_capacity,
    critical_curve_t* cell_curves, double* breakpoints,
    size_t breakpoint_capacity,
    uint64_t max_active_boundary_tests,
    int* out_found_path, alea_transition_slice_stats_t* stats) {
    static const double fractions[3] = {0.0, 0.5, 1.0};
    size_t occurrence_count = 0;
    *out_found_path = 0;
    for (int iu = 0; iu < 3; iu++) {
        const double u = tile->uv_min[0] + fractions[iu] *
            (tile->uv_max[0] - tile->uv_min[0]);
        for (int iv = 0; iv < 3; iv++) {
            const double v = tile->uv_min[1] + fractions[iv] *
                (tile->uv_max[1] - tile->uv_min[1]);
            double world[3];
            tile_world_point(view, u, v, world);
            stats->critical_occurrence_seed_points++;
            alea_hier_cell_hit_t deepest;
            alea_hier_ray_path_t path;
            const int found = alea_hier_spatial_find_path_at_point(
                sys, world[0], world[1], world[2], &deepest, &path);
            if (found < 0) return -1;
            if (!found || path.count <= 0) continue;
            *out_found_path = 1;
            stats->critical_occurrence_paths++;
            /* Query the leaf before considering another seed.  This is the
             * only broad same-level query required for this path. */
            {
                const int level = path.count - 1;
                const alea_hier_ray_path_entry_t* entry = &path.entries[level];
                const uint64_t universe_key =
                    path_universe_occurrence_key(&path, level);
                int duplicate_occurrence = 0;
                for (size_t oi = 0; oi < occurrence_count; oi++) {
                    if (occurrences[oi].universe_id == entry->universe_id &&
                        occurrences[oi].universe_occurrence_key == universe_key) {
                        duplicate_occurrence = 1;
                        break;
                    }
                }
                if (duplicate_occurrence) continue;
                if (occurrence_count == occurrence_capacity)
                    return CRITICAL_COLLECT_CHAIN_TRUNCATED;
                occurrences[occurrence_count++] =
                    (critical_universe_occurrence_t){
                        .universe_id = entry->universe_id,
                        .universe_occurrence_key = universe_key
                    };

                const alea_bbox_t local_bbox =
                    tile_local_bbox(view, tile, &entry->transform);
                critical_region_visit_t visit = {
                    .sys = sys, .view = view, .tile = tile,
                    .path = &path, .level = level,
                    .universe_key = universe_key,
                    .curves = curves, .curve_count = curve_count,
                    .curve_capacity = curve_capacity,
                    .cell_curves = cell_curves,
                    .cell_curve_capacity = curve_capacity,
                    .breakpoints = breakpoints,
                    .breakpoint_capacity = breakpoint_capacity,
                    .max_active_boundary_tests = max_active_boundary_tests,
                    .stats = stats
                };
                size_t visited = 0;
                if (alea_hier_spatial_visit_universe_region(
                        sys, entry->universe_id, &local_bbox,
                        visit_region_cell_curves, &visit, &visited) != 0)
                    return -1;
                stats->critical_occurrence_universe_queries++;
                stats->critical_region_candidates_scanned += visited;
                if (visit.stop_code != CRITICAL_COLLECT_OK)
                    return visit.stop_code;
            }
            for (int level = path.count - 2; level >= 0; level--) {
                const uint64_t universe_key =
                    path_universe_occurrence_key(&path, level);
                critical_region_visit_t visit = {
                    .sys = sys, .view = view, .tile = tile,
                    .path = &path, .level = level,
                    .universe_key = universe_key,
                    .curves = curves, .curve_count = curve_count,
                    .curve_capacity = curve_capacity,
                    .cell_curves = cell_curves,
                    .cell_curve_capacity = curve_capacity,
                    .breakpoints = breakpoints,
                    .breakpoint_capacity = breakpoint_capacity,
                    .max_active_boundary_tests = max_active_boundary_tests,
                    .stats = stats
                };
                const int rc = collect_ancestor_neighborhood(&visit);
                if (rc != CRITICAL_COLLECT_OK) return rc;
            }
        }
    }
    if (!*out_found_path) {
        alea_hier_ray_path_t root_path;
        memset(&root_path, 0, sizeof(root_path));
        root_path.count = 1;
        root_path.entries[0].universe_id = 0;
        alea_matrix_identity(&root_path.entries[0].transform);
        const alea_bbox_t root_bbox =
            tile_local_bbox(view, tile, &root_path.entries[0].transform);
        critical_region_visit_t visit = {
            .sys = sys, .view = view, .tile = tile,
            .path = &root_path, .level = 0,
            .universe_key = path_universe_occurrence_key(&root_path, 0),
            .curves = curves, .curve_count = curve_count,
            .curve_capacity = curve_capacity,
            .cell_curves = cell_curves,
            .cell_curve_capacity = curve_capacity,
            .breakpoints = breakpoints,
            .breakpoint_capacity = breakpoint_capacity,
            .max_active_boundary_tests = max_active_boundary_tests,
            .stats = stats
        };
        size_t visited = 0;
        if (alea_hier_spatial_visit_universe_region(
                sys, 0, &root_bbox, visit_region_cell_curves,
                &visit, &visited) != 0)
            return -1;
        stats->critical_root_region_fallbacks++;
        stats->critical_occurrence_universe_queries++;
        stats->critical_region_candidates_scanned += visited;
        if (visit.stop_code != CRITICAL_COLLECT_OK) return visit.stop_code;
    }
    return CRITICAL_COLLECT_OK;
}

static int sample_boundary_piece(
    const critical_curve_t* curve,
    alea_transition_slice_boundary_piece_t* piece) {
    const size_t count = ALEA_TRANSITION_SLICE_BOUNDARY_POINT_CAPACITY;
    if (curve->has_scanline_domain) {
        double previous_u = curve->scanline_endpoint_u[0];
        for (size_t i = 0; i < count; i++) {
            const double fraction = (double)i / (double)(count - 1u);
            const double v = curve->scanline_v_min + fraction *
                (curve->scanline_v_max - curve->scanline_v_min);
            double roots[8];
            const int root_count = alea_curve_scanline_intersect(
                &curve->curve, v, roots,
                (int)(sizeof(roots) / sizeof(roots[0])));
            if (root_count <= 0) return 0;
            sort_small_doubles(roots, root_count);
            int selected = 0;
            for (int root = 1; root < root_count; root++)
                if (fabs(roots[root] - previous_u) <
                    fabs(roots[selected] - previous_u)) selected = root;
            previous_u = roots[selected];
            piece->uv[i][0] = previous_u;
            piece->uv[i][1] = v;
        }
    } else {
        double parameter_min, parameter_max;
        if (curve->has_parameter_domain) {
            parameter_min = curve->parameter_min;
            parameter_max = curve->parameter_max;
        } else if (curve_is_line(curve->curve.type) ||
                   curve_is_circle(curve->curve.type) ||
                   curve_is_ellipse(curve->curve.type)) {
            parameter_min = curve->curve.bounds.t_min;
            parameter_max = curve->curve.bounds.t_max;
        } else {
            return 0;
        }
        if (!isfinite(parameter_min) || !isfinite(parameter_max) ||
            parameter_max < parameter_min) return 0;
        for (size_t i = 0; i < count; i++) {
            const double fraction = (double)i / (double)(count - 1u);
            const double parameter = parameter_min + fraction *
                (parameter_max - parameter_min);
            int evaluated;
            if (curve->has_parameter_domain)
                evaluated = open_conic_eval(
                    curve, parameter, curve->conic_branch,
                    &piece->uv[i][0], &piece->uv[i][1]);
            else
                evaluated = alea_curve_eval(
                    &curve->curve, parameter,
                    &piece->uv[i][0], &piece->uv[i][1]);
            if (!evaluated) return 0;
        }
    }
    piece->point_count = count;
    return 1;
}

static const critical_curve_t* find_boundary_curve(
    const critical_curve_t* curves, size_t curve_count,
    int surface_id, uint64_t universe_occurrence_key,
    const double uv[2], double tolerance) {
    for (size_t i = 0; i < curve_count; i++)
        if (curves[i].curve.surface_id == surface_id &&
            curves[i].universe_occurrence_key == universe_occurrence_key &&
            critical_curve_contains_point(
                &curves[i], uv[0], uv[1], tolerance))
            return &curves[i];
    return NULL;
}

static void retain_boundary_evidence(
    alea_transition_slice_critical_finding_t* finding,
    const critical_curve_t* curves, size_t curve_count,
    const critical_curve_t* source, const double uv[2], double tolerance,
    const alea_transition_slice_options_t* options,
    alea_transition_slice_stats_t* stats) {
    if (options->max_critical_boundary_evidence &&
        stats->critical_boundary_evidence >=
            options->max_critical_boundary_evidence) {
        finding->boundary_evidence_truncated = 1;
        stats->omitted_critical_boundary_evidence++;
        return;
    }
    const int surface_ids[3] = {
        source->curve.surface_id,
        finding->transition.primary_surface_id,
        finding->transition.connecting_surface_id};
    const uint32_t roles[3] = {
        ALEA_TRANSITION_SLICE_BOUNDARY_ROLE_SOURCE,
        ALEA_TRANSITION_SLICE_BOUNDARY_ROLE_PRIMARY,
        ALEA_TRANSITION_SLICE_BOUNDARY_ROLE_CONNECTING};
    for (size_t role = 0; role < 3; role++) {
        const int surface_id = surface_ids[role];
        if (surface_id <= 0) continue;
        size_t existing = finding->boundary_piece_count;
        for (size_t i = 0; i < finding->boundary_piece_count; i++)
            if (finding->boundary_pieces[i].surface_id == surface_id) {
                existing = i;
                break;
            }
        if (existing < finding->boundary_piece_count) {
            finding->boundary_pieces[existing].role_flags |= roles[role];
            continue;
        }
        if (finding->boundary_piece_count >=
            ALEA_TRANSITION_SLICE_BOUNDARY_PIECE_CAPACITY) {
            finding->boundary_evidence_truncated = 1;
            continue;
        }
        const critical_curve_t* curve = role == 0 ? source :
            find_boundary_curve(
                curves, curve_count, surface_id,
                source->universe_occurrence_key, uv, tolerance);
        if (!curve) {
            finding->boundary_evidence_truncated = 1;
            continue;
        }
        alea_transition_slice_boundary_piece_t* piece =
            &finding->boundary_pieces[finding->boundary_piece_count];
        piece->surface_id = surface_id;
        piece->role_flags = roles[role];
        if (!sample_boundary_piece(curve, piece)) {
            memset(piece, 0, sizeof(*piece));
            finding->boundary_evidence_truncated = 1;
            continue;
        }
        finding->boundary_piece_count++;
    }
    if (finding->boundary_piece_count > 0)
        stats->critical_boundary_evidence++;
    if (finding->boundary_evidence_truncated)
        stats->omitted_critical_boundary_evidence++;
}

static void set_saturated(alea_transition_slice_stats_t* stats,
                          alea_transition_slice_critical_stop_reason_t reason) {
    stats->critical_tiles_saturated++;
    if (stats->critical_stop_reason == ALEA_TRANSITION_SLICE_CRITICAL_NONE)
        stats->critical_stop_reason = reason;
}

int alea_transition_slice_enumerate_critical_tiles(
    alea_system_t* sys, const alea_slice_view_t* view,
    const alea_transition_slice_options_t* options,
    const alea_transition_slice_critical_tile_t* tiles, size_t tile_count,
    alea_transition_slice_critical_finding_sink_t finding_sink,
    void* finding_sink_userdata,
    alea_transition_slice_stats_t* stats) {
    if (!sys || !view || !options || !stats || (tile_count && !tiles)) return -1;
    if (!tile_count) return 0;

    const size_t max_curves = options->max_curves_per_tile;
    const size_t max_points = options->max_critical_points;
    const size_t coverage_capacity = options->max_coverage_hits;
    const size_t point_slot_capacity = hash_capacity(max_points);
    const size_t occurrence_capacity = 9u * ALEA_HIER_RAY_PATH_MAX;
    if (!max_curves || max_curves > (SIZE_MAX - 10) / 2 ||
        !max_points || !coverage_capacity ||
        !point_slot_capacity) {
        set_saturated(stats, ALEA_TRANSITION_SLICE_CRITICAL_MAX_SCRATCH_BYTES);
        return 0;
    }
    size_t scratch_bytes = 0;
    if (checked_add(&scratch_bytes, occurrence_capacity,
                    sizeof(critical_universe_occurrence_t)) != 0 ||
        checked_add(&scratch_bytes, max_curves, sizeof(critical_curve_t)) != 0 ||
        checked_add(&scratch_bytes, max_curves, sizeof(critical_curve_t)) != 0 ||
        checked_add(&scratch_bytes, 2 * max_curves + 10, sizeof(double)) != 0 ||
        checked_add(&scratch_bytes, max_points,
                    sizeof(critical_point_t)) != 0 ||
        checked_add(&scratch_bytes, point_slot_capacity,
                    sizeof(critical_point_slot_t)) != 0 ||
        checked_add(&scratch_bytes, max_curves,
                    sizeof(critical_curve_order_t)) != 0 ||
        checked_add(&scratch_bytes, coverage_capacity,
                    sizeof(alea_cell_hit_t) + 2*sizeof(uint64_t) +
                    sizeof(uint8_t)) != 0 ||
        scratch_bytes > options->max_critical_scratch_bytes) {
        set_saturated(stats, ALEA_TRANSITION_SLICE_CRITICAL_MAX_SCRATCH_BYTES);
        return 0;
    }

    critical_universe_occurrence_t* occurrences =
        calloc(occurrence_capacity, sizeof(*occurrences));
    critical_curve_t* curves = calloc(max_curves, sizeof(*curves));
    critical_curve_t* cell_curves = calloc(max_curves, sizeof(*cell_curves));
    const size_t breakpoint_capacity = 2 * max_curves + 10;
    double* breakpoints = calloc(breakpoint_capacity, sizeof(*breakpoints));
    critical_point_t* points = calloc(max_points, sizeof(*points));
    critical_point_slot_t* point_slots =
        calloc(point_slot_capacity, sizeof(*point_slots));
    critical_curve_order_t* order = calloc(max_curves, sizeof(*order));
    alea_cell_hit_t* coverage_hits = calloc(coverage_capacity,
                                             sizeof(*coverage_hits));
    uint64_t* coverage_keys = calloc(coverage_capacity, sizeof(*coverage_keys));
    uint64_t* coverage_parent_keys = calloc(
        coverage_capacity, sizeof(*coverage_parent_keys));
    uint8_t* coverage_mask = calloc(coverage_capacity, sizeof(*coverage_mask));
    if (!occurrences || !curves || !cell_curves || !breakpoints ||
        !points ||
        !point_slots || !order ||
        !coverage_hits || !coverage_keys || !coverage_parent_keys ||
        !coverage_mask) {
        free(occurrences); free(curves); free(cell_curves); free(breakpoints);
        free(points); free(point_slots);
        free(order);
        free(coverage_hits); free(coverage_keys); free(coverage_parent_keys);
        free(coverage_mask);
        return -1;
    }
    alea_raycast_result_t ray_scratch;
    alea_ray_boundary_event_result_t ray_events;
    alea_raycast_result_init(&ray_scratch);
    alea_ray_boundary_event_result_init(&ray_events);
    if (scratch_bytes > stats->peak_critical_scratch_bytes)
        stats->peak_critical_scratch_bytes = scratch_bytes;

    for (size_t ti = 0; ti < tile_count; ti++) {
        memset(point_slots, 0, point_slot_capacity * sizeof(*point_slots));
        size_t curve_count = 0;
        size_t point_count = 0;
        int found_path = 0;
        const int collect_rc = collect_tile_occurrence_curves(
            sys, view, &tiles[ti], occurrences, occurrence_capacity,
            curves, &curve_count, max_curves, cell_curves, breakpoints,
            breakpoint_capacity,
            options->max_active_boundary_tests,
            &found_path, stats);
        if (collect_rc < 0) {
            free(occurrences); free(curves); free(cell_curves);
            free(breakpoints); free(points);
            free(point_slots); free(order);
            free(coverage_hits); free(coverage_keys);
            free(coverage_parent_keys); free(coverage_mask);
            alea_raycast_result_free(&ray_scratch);
            alea_ray_boundary_event_result_free(&ray_events);
            return -1;
        }
        stats->critical_tiles_processed++;
        if (collect_rc == CRITICAL_COLLECT_MAX_CURVES) {
            set_saturated(stats, ALEA_TRANSITION_SLICE_CRITICAL_MAX_CURVES);
        }
        if (stats->critical_ranked_curves_omitted > 0)
            set_saturated(stats, ALEA_TRANSITION_SLICE_CRITICAL_MAX_CURVES);
        if (collect_rc == CRITICAL_COLLECT_CHAIN_TRUNCATED) {
            stats->critical_chain_truncated_hits++;
            set_saturated(stats,
                          ALEA_TRANSITION_SLICE_CRITICAL_CHAIN_TRUNCATED);
            continue;
        }
        stats->critical_curves += curve_count;
        if (curve_count > stats->peak_critical_curves)
            stats->peak_critical_curves = curve_count;
        qsort(curves, curve_count, sizeof(*curves), curve_priority_compare);
        const double tile_scale = fmax(1.0, fmax(
            tiles[ti].uv_max[0] - tiles[ti].uv_min[0],
            tiles[ti].uv_max[1] - tiles[ti].uv_min[1]));
        const double point_tolerance = 1e-10 * tile_scale;
        const size_t single_point_capacity = max_points / 4u > 0
            ? max_points / 4u : max_points;
        size_t single_curve_index = 0;
        for (; single_curve_index < curve_count; single_curve_index++) {
            if (point_count >= single_point_capacity) {
                set_saturated(
                    stats, ALEA_TRANSITION_SLICE_CRITICAL_MAX_POINTS);
                break;
            }
            const int point_rc = generate_single_curve_points(
                &curves[single_curve_index], single_curve_index, &tiles[ti],
                points, &point_count,
                single_point_capacity,
                point_slots, point_slot_capacity, point_tolerance, stats);
            if (point_rc < 0) {
                set_saturated(
                    stats, ALEA_TRANSITION_SLICE_CRITICAL_MAX_POINTS);
                break;
            }
            if (point_rc > 0 && stats->critical_stop_reason ==
                    ALEA_TRANSITION_SLICE_CRITICAL_NONE)
                stats->critical_stop_reason =
                    ALEA_TRANSITION_SLICE_CRITICAL_UNSUPPORTED_CURVE;
        }
        for (size_t ci = 0; ci < curve_count; ci++) {
            order[ci].curve_index = ci;
            order[ci].u_min = curves[ci].bbox[0];
            order[ci].u_max = curves[ci].bbox[1];
            order[ci].v_min = curves[ci].bbox[2];
            order[ci].v_max = curves[ci].bbox[3];
        }
        qsort(order, curve_count, sizeof(*order), curve_order_compare);
        int pair_stopped = 0;
        const uint64_t local_pair_budget = options->max_curve_pairs
            ? options->max_curve_pairs / 2u : UINT64_MAX;
        uint64_t local_pairs = 0;
        int local_pair_stopped = 0;
        for (size_t first = 0; first < curve_count && !pair_stopped &&
             !local_pair_stopped; first++) {
            for (size_t second = first + 1; second < curve_count; second++) {
                if (curves[first].cell_id != curves[second].cell_id) continue;
                if (critical_curves_same_source(
                        &curves[first], &curves[second])) continue;
                if (curves[first].bbox[1] < curves[second].bbox[0] ||
                    curves[second].bbox[1] < curves[first].bbox[0] ||
                    curves[first].bbox[3] < curves[second].bbox[2] ||
                    curves[second].bbox[3] < curves[first].bbox[2]) continue;
                stats->critical_curve_pair_candidates++;
                if (local_pairs >= local_pair_budget) {
                    set_saturated(
                        stats, ALEA_TRANSITION_SLICE_CRITICAL_MAX_CURVE_PAIRS);
                    local_pair_stopped = 1;
                    break;
                }
                local_pairs++;
                stats->critical_curve_pairs_tested++;
                const int pair_rc = solve_curve_pair(
                    &curves[first], &curves[second], first,
                    &tiles[ti], points, &point_count, max_points, point_slots,
                    point_slot_capacity, point_tolerance, stats);
                if (pair_rc < 0) {
                    set_saturated(stats,
                                  ALEA_TRANSITION_SLICE_CRITICAL_MAX_POINTS);
                    pair_stopped = 1;
                    break;
                }
                if (pair_rc > 0) {
                    record_unsupported_pair(
                        curves[first].curve.type, curves[second].curve.type,
                        stats);
                    if (stats->critical_stop_reason ==
                            ALEA_TRANSITION_SLICE_CRITICAL_NONE)
                        stats->critical_stop_reason =
                            ALEA_TRANSITION_SLICE_CRITICAL_UNSUPPORTED_CURVE;
                }
            }
        }
        for (size_t oi = 0; oi < curve_count && !pair_stopped; oi++) {
            for (size_t oj = oi + 1; oj < curve_count; oj++) {
                if (order[oj].u_min > order[oi].u_max + point_tolerance)
                    break;
                if (order[oj].v_min > order[oi].v_max + point_tolerance ||
                    order[oi].v_min > order[oj].v_max + point_tolerance)
                    continue;
                stats->critical_curve_pair_candidates++;
                const size_t first = order[oi].curve_index;
                const size_t second = order[oj].curve_index;
                if (curves[first].cell_id == curves[second].cell_id)
                    continue;
                if (critical_curves_same_source(
                        &curves[first], &curves[second])) continue;
                if (options->max_curve_pairs &&
                    stats->critical_curve_pairs_tested >=
                        options->max_curve_pairs) {
                    set_saturated(stats,
                        ALEA_TRANSITION_SLICE_CRITICAL_MAX_CURVE_PAIRS);
                    pair_stopped = 1;
                    break;
                }
                stats->critical_curve_pairs_tested++;
                const int pair_rc = solve_curve_pair(
                    &curves[first], &curves[second], first,
                    &tiles[ti], points, &point_count, max_points, point_slots,
                    point_slot_capacity, point_tolerance, stats);
                if (pair_rc < 0) {
                    set_saturated(stats,
                                  ALEA_TRANSITION_SLICE_CRITICAL_MAX_POINTS);
                    pair_stopped = 1;
                    break;
                }
                if (pair_rc > 0) {
                    record_unsupported_pair(
                        curves[first].curve.type, curves[second].curve.type,
                        stats);
                    if (stats->critical_stop_reason ==
                            ALEA_TRANSITION_SLICE_CRITICAL_NONE)
                        stats->critical_stop_reason =
                            ALEA_TRANSITION_SLICE_CRITICAL_UNSUPPORTED_CURVE;
                }
            }
        }
        stats->critical_points += point_count;

        for (size_t pi = 0; pi < point_count; pi++) {
            if (options->max_critical_probes && stats->critical_probes >=
                    options->max_critical_probes) {
                set_saturated(stats,
                              ALEA_TRANSITION_SLICE_CRITICAL_MAX_PROBES);
                break;
            }
            double normal[2];
            if (!curve_normal_at(&curves[points[pi].curve_index].curve,
                                 points[pi].uv, normal))
                continue;
            double radius = options->critical_probe_radius > 0.0
                ? options->critical_probe_radius : 1e-6 * tile_scale;
            if (!(radius > 64.0 * RAY_EPSILON))
                radius = 64.0 * RAY_EPSILON;
            const double start_u = points[pi].uv[0] - radius * normal[0];
            const double start_v = points[pi].uv[1] - radius * normal[1];
            double origin[3], direction[3];
            for (int axis = 0; axis < 3; axis++) {
                origin[axis] = view->plane.origin[axis] +
                    start_u * view->plane.u_axis[axis] +
                    start_v * view->plane.v_axis[axis];
                direction[axis] = normal[0] * view->plane.u_axis[axis] +
                    normal[1] * view->plane.v_axis[axis];
            }
            alea_ray_t ray;
            if (alea_ray_init(&ray, origin[0], origin[1], origin[2],
                              direction[0], direction[1], direction[2]) != 0)
                continue;
            const alea_ray_boundary_event_options_internal_t event_options = {
                .max_events = 64,
                .max_output_bytes = options->max_critical_scratch_bytes,
                .skip_open_side_coverage = true
            };
            if (alea_raycast_selected_boundary_events_with_options_nocache(
                    sys, &ray, 2.0 * radius, &event_options, &ray_scratch,
                    &ray_events) != 0) {
                if (alea_error_code() == ALEA_ERR_OVERFLOW) {
                    set_saturated(stats,
                                  ALEA_TRANSITION_SLICE_CRITICAL_MAX_PROBES);
                    break;
                }
                free(occurrences); free(curves); free(cell_curves);
                free(breakpoints); free(points);
                free(point_slots); free(order);
                free(coverage_hits); free(coverage_keys);
                free(coverage_parent_keys); free(coverage_mask);
                alea_raycast_result_free(&ray_scratch);
                alea_ray_boundary_event_result_free(&ray_events);
                return -1;
            }
            stats->critical_probes++;
            const double center_tolerance = fmax(32.0 * RAY_EPSILON,
                                                 radius * 1e-3);
            for (size_t ei = 0; ei < ray_events.events.count; ei++) {
                const alea_ray_boundary_event_t* event =
                    &ray_events.events.data[ei];
                if (event->kind != ALEA_RAY_BOUNDARY_EVENT_PHYSICAL ||
                    fabs(event->t - radius) > center_tolerance)
                    continue;
                stats->critical_probe_events++;
                alea_transition_options_t transition_options;
                alea_transition_options_init(&transition_options);
                transition_options.probe_distance = fmin(radius * 0.25,
                    options->probe_distance > 0.0
                        ? options->probe_distance : radius * 0.25);
                transition_options.max_probe_distance = radius * 0.25;
                transition_options.max_coverage_hits = options->max_coverage_hits;
                alea_transition_result_t transition;
                if (alea_check_selected_boundary_event_transition_nocache(
                        sys, event, &transition_options, &transition) != 0) {
                    free(occurrences); free(curves); free(cell_curves);
                    free(breakpoints); free(points);
                    free(point_slots); free(order);
                    free(coverage_hits); free(coverage_keys);
                    free(coverage_parent_keys); free(coverage_mask);
                    alea_raycast_result_free(&ray_scratch);
                    alea_ray_boundary_event_result_free(&ray_events);
                    return -1;
                }
                const critical_curve_t* source =
                    &curves[points[pi].curve_index];
                if (source->curve.surface_id > 0 &&
                    source->curve.surface_id != transition.primary_surface_id &&
                    (event->provenance_flags &
                        ALEA_BOUNDARY_PROVENANCE_ACTIVE_FRAME) &&
                    transition.current_cell_id > 0) {
                    alea_transition_result_t source_transition;
                    if (alea_check_transition_local(
                            sys, event->active_universe_id,
                            transition.current_cell_id,
                            source->curve.surface_id, NULL, 0,
                            event->local_point, event->local_direction,
                            &transition_options, &source_transition) == 0 &&
                        source_transition.kind != ALEA_TRANSITION_VALID) {
                        source_transition.occurrence_depth =
                            transition.occurrence_depth;
                        source_transition.current_occurrence_key =
                            transition.current_occurrence_key;
                        source_transition.current_parent_occurrence_key =
                            transition.current_parent_occurrence_key;
                        source_transition.before_occurrence_key =
                            transition.before_occurrence_key;
                        source_transition.before_parent_occurrence_key =
                            transition.before_parent_occurrence_key;
                        source_transition.selected_after_occurrence_key =
                            transition.selected_after_occurrence_key;
                        source_transition.selected_after_parent_occurrence_key =
                            transition.selected_after_parent_occurrence_key;
                        transition = source_transition;
                    }
                }
                stats->coverage_fallbacks += transition.coverage_fallbacks;
                if (transition.kind != ALEA_TRANSITION_VALID) {
                    stats->critical_probe_findings++;
                    if (finding_sink) {
                        alea_transition_slice_critical_finding_t finding;
                        memset(&finding, 0, sizeof(finding));
                        finding.transition = transition;
                        finding.tile_index = ti;
                        finding.point_index = pi;
                        finding.source_cell_id = source->cell_id;
                        finding.source_surface_id = source->curve.surface_id;
                        finding.source_occurrence_key = source->occurrence_key;
                        finding.source_universe_occurrence_key =
                            source->universe_occurrence_key;
                        finding.uv[0] = points[pi].uv[0];
                        finding.uv[1] = points[pi].uv[1];
                        tile_world_point(
                            view, finding.uv[0], finding.uv[1],
                            finding.world_point);
                        memcpy(finding.direction, direction,
                               sizeof(finding.direction));
                        finding.radius = radius;
                        retain_boundary_evidence(
                            &finding, curves, curve_count, source,
                            finding.uv, point_tolerance, options, stats);
                        const int sink_rc = finding_sink(
                            &finding, finding_sink_userdata);
                        if (sink_rc < 0) {
                            free(occurrences); free(curves); free(cell_curves);
                            free(breakpoints); free(points);
                            free(point_slots); free(order);
                            free(coverage_hits); free(coverage_keys);
                            free(coverage_parent_keys); free(coverage_mask);
                            alea_raycast_result_free(&ray_scratch);
                            alea_ray_boundary_event_result_free(&ray_events);
                            return -1;
                        }
                        if (sink_rc > 0) {
                            if (finding.boundary_piece_count > 0)
                                stats->critical_boundary_evidence--;
                            if (finding.boundary_evidence_truncated)
                                stats->omitted_critical_boundary_evidence--;
                            set_saturated(
                                stats,
                                ALEA_TRANSITION_SLICE_CRITICAL_MAX_FINDINGS);
                            pi = point_count;
                            break;
                        }
                    }
                }
            }
        }
        for (size_t pi = 0; pi < point_count; pi++) {
            double normal[2];
            if (!curve_normal_at(&curves[points[pi].curve_index].curve,
                                 points[pi].uv, normal)) continue;
            const double tangent[2] = {-normal[1], normal[0]};
            const double radius = options->critical_probe_radius > 0.0
                ? options->critical_probe_radius * 0.5 : 5e-7 * tile_scale;
            for (int sector = 0; sector < 4; sector++) {
                if (options->max_critical_sector_witnesses &&
                    stats->critical_sector_witnesses >=
                        options->max_critical_sector_witnesses) {
                    set_saturated(stats,
                        ALEA_TRANSITION_SLICE_CRITICAL_MAX_SECTOR_WITNESSES);
                    pi = point_count;
                    break;
                }
                const double sn = (sector & 1) ? 1.0 : -1.0;
                const double st = (sector & 2) ? 1.0 : -1.0;
                const double u = points[pi].uv[0] +
                    radius * (sn*normal[0] + st*tangent[0]) *
                    0.70710678118654752440;
                const double v = points[pi].uv[1] +
                    radius * (sn*normal[1] + st*tangent[1]) *
                    0.70710678118654752440;
                double world[3];
                for (int axis = 0; axis < 3; axis++)
                    world[axis] = view->plane.origin[axis] +
                        u*view->plane.u_axis[axis] +
                        v*view->plane.v_axis[axis];
                const int count =
                    alea_find_all_cells_at_point_coverage_chain_recursive(
                        sys, world[0], world[1], world[2], coverage_hits,
                        coverage_keys, coverage_parent_keys, coverage_capacity);
                if (count < 0) {
                    free(occurrences); free(curves); free(cell_curves);
                    free(breakpoints); free(points);
                    free(point_slots); free(order); free(coverage_hits);
                    free(coverage_keys); free(coverage_parent_keys);
                    free(coverage_mask);
                    alea_raycast_result_free(&ray_scratch);
                    alea_ray_boundary_event_result_free(&ray_events);
                    return -1;
                }
                stats->critical_sector_witnesses++;
                if ((size_t)count >= coverage_capacity) {
                    stats->critical_sector_unresolved_witnesses++;
                    continue;
                }
                alea_point_coverage_classification_t classification;
                memset(coverage_mask, 0, coverage_capacity);
                if (alea_classify_point_coverage_chain(
                        coverage_hits, coverage_keys, coverage_parent_keys,
                        (size_t)count, -1, coverage_mask, &classification) != 0)
                    continue;
                if (classification.kind == ALEA_POINT_COVERAGE_GAP)
                    stats->critical_sector_gap_witnesses++;
                else if (classification.kind == ALEA_POINT_COVERAGE_OVERLAP)
                    stats->critical_sector_overlap_witnesses++;
                else if (classification.kind != ALEA_POINT_COVERAGE_UNIQUE)
                    stats->critical_sector_unresolved_witnesses++;
            }
        }
    }
    free(occurrences); free(curves); free(cell_curves); free(breakpoints);
    free(points); free(point_slots);
    free(order);
    free(coverage_hits); free(coverage_keys); free(coverage_parent_keys);
    free(coverage_mask);
    alea_raycast_result_free(&ray_scratch);
    alea_ray_boundary_event_result_free(&ray_events);
    return 0;
}
