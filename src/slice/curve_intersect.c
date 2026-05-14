// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/**
 * @file curve_intersect.c
 * @brief Analytical surface-plane intersection for 2D slice rendering
 *
 * The math here is straightforward algebraic substitution:
 * Given a slice plane P(u,v) = O + u*U + v*V and a surface f(x,y,z) = 0,
 * substitute the plane equation into the surface equation to get f(u,v) = 0,
 * which is the 2D intersection curve.
 */

#include "curve_intersect.h"
#include "core/alea_system.h"
#include "core/alea_spatial.h"
#include "core/alea_spatial_hier.h"
#include "core/alea_universe.h"
#include "primitives/primitive_desc.h"
#include "util/poly_solve.h"
#include "util/math.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <stdio.h>
#include "util/alea_log.h"

#include "slice/slice_epsilon.h"
#define EPSILON SLICE_DISCRIMINANT_TOL  /* legacy alias — use named constants for new code */

/* Solve quadratic ax² + bx + c = 0, return number of real roots.
 * Uses Vieta's form (same as raycast module) to avoid catastrophic
 * cancellation when b is large: q = -0.5*(b + copysign(sqrt(disc), b)),
 * then x1 = q/a, x2 = c/q. */
static int solve_quadratic(double a, double b, double c, double* x1, double* x2) {
    if (fabs(a) < EPSILON) {
        /* Linear case */
        if (fabs(b) < EPSILON) return 0;
        *x1 = -c / b;
        return 1;
    }

    double disc = b*b - 4*a*c;
    if (disc < -EPSILON) return 0;
    if (disc < EPSILON) {
        *x1 = -b / (2*a);
        return 1;
    }

    double sqrt_disc = sqrt(disc);
    double q = -0.5 * (b + copysign(sqrt_disc, b));
    *x1 = q / a;
    *x2 = c / q;
    if (*x1 > *x2) { double tmp = *x1; *x1 = *x2; *x2 = tmp; }
    return 2;
}

/**
 * @brief Convert general conic to canonical ellipse form
 *
 * Given: A*u² + B*u*v + C*v² + D*u + E*v + F = 0
 * Compute: center, semi-axes, rotation angle
 *
 * Returns true if conic is an ellipse, false otherwise.
 */
bool alea_conic_to_ellipse(const alea_conic_2d_t* conic, alea_ellipse_2d_t* ellipse) {
    double A = conic->A;
    double B = conic->B;
    double C = conic->C;
    double D = conic->D;
    double E = conic->E;
    double F = conic->F;

    /* Discriminant: B² - 4AC < 0 for ellipse */
    double disc = B*B - 4*A*C;
    if (disc >= -EPSILON) {
        return false;  /* Not an ellipse */
    }

    /* Find center by solving: ∂f/∂u = 0, ∂f/∂v = 0 */
    /* 2Au + Bv + D = 0 */
    /* Bu + 2Cv + E = 0 */
    double denom = 4*A*C - B*B;
    if (fabs(denom) < EPSILON) return false;

    double cx = (B*E - 2*C*D) / denom;
    double cy = (B*D - 2*A*E) / denom;

    /* Translate to center: substitute u' = u - cx, v' = v - cy */
    /* New constant term: F' = A*cx² + B*cx*cy + C*cy² + D*cx + E*cy + F */
    double F_prime = A*cx*cx + B*cx*cy + C*cy*cy + D*cx + E*cy + F;

    if (fabs(F_prime) < EPSILON) return false;

    /* Normalize: divide by -F' so RHS = 1 */
    double A_n = -A / F_prime;
    double B_n = -B / F_prime;
    double C_n = -C / F_prime;

    /* Find rotation angle: tan(2θ) = B / (A - C) */
    double theta;
    if (fabs(A_n - C_n) < EPSILON) {
        theta = (B_n > 0) ? M_PI/4 : -M_PI/4;
    } else {
        theta = 0.5 * atan2(B_n, A_n - C_n);
    }

    /* Rotate to eliminate cross term */
    double cos_t = cos(theta);
    double sin_t = sin(theta);
    double cos2 = cos_t * cos_t;
    double sin2 = sin_t * sin_t;
    double sincos = sin_t * cos_t;

    /* A' = A*cos² + B*sin*cos + C*sin² */
    /* C' = A*sin² - B*sin*cos + C*cos² */
    double A_rot = A_n*cos2 + B_n*sincos + C_n*sin2;
    double C_rot = A_n*sin2 - B_n*sincos + C_n*cos2;

    /* Semi-axes: a² = 1/A', b² = 1/C' */
    if (A_rot <= EPSILON || C_rot <= EPSILON) return false;

    double a = sqrt(1.0 / A_rot);
    double b = sqrt(1.0 / C_rot);

    ellipse->center[0] = cx;
    ellipse->center[1] = cy;
    if (a >= b) {
        ellipse->semi_a = a;
        ellipse->semi_b = b;
        ellipse->angle = theta;
    } else {
        ellipse->semi_a = b;
        ellipse->semi_b = a;
        ellipse->angle = theta + M_PI/2;
    }

    /* Normalize angle to [-π, π] */
    while (ellipse->angle > M_PI) ellipse->angle -= 2*M_PI;
    while (ellipse->angle < -M_PI) ellipse->angle += 2*M_PI;

    return true;
}

/**
 * @brief Try to convert a conic-form ellipse curve to canonical ellipse form.
 *
 * If the curve type is ALEA_CURVE_ELLIPSE and the conic can be converted
 * to canonical form (center, semi-axes, angle), overwrite data.ellipse.
 * On failure (degenerate conic), fall back to HYPERBOLA type so scanline
 * rendering uses the conic path.
 */
static void finalize_conic_as_ellipse(alea_curve_2d_t* curve) {
    if (curve->type != ALEA_CURVE_ELLIPSE) return;
    alea_conic_2d_t conic = curve->data.conic;  /* save before union overwrite */
    if (alea_conic_to_ellipse(&conic, &curve->data.ellipse)) {
        /* Successfully converted — type stays ALEA_CURVE_ELLIPSE,
         * data.ellipse is now valid for eval/scanline/bbox. */
    } else {
        /* Degenerate — restore conic and reclassify */
        curve->data.conic = conic;
        curve->type = ALEA_CURVE_HYPERBOLA;
    }
}

/**
 * @brief Compute bounding box of a conic section
 */
static void conic_bbox(const alea_conic_2d_t* conic,
                       double* u_min, double* u_max,
                       double* v_min, double* v_max) {
    /* Parabola/hyperbola: unbounded, use large extent */
    (void)conic;
    *u_min = *v_min = -1e10;
    *u_max = *v_max = 1e10;
}

/* ============================================================================
 * SLICE PLANE INITIALIZATION
 * ============================================================================ */

void alea_slice_plane_init(alea_slice_plane_t* plane,
                          double ox, double oy, double oz,
                          double nx, double ny, double nz,
                          double ux, double uy, double uz) {
    /* Set origin */
    plane->origin[0] = ox;
    plane->origin[1] = oy;
    plane->origin[2] = oz;

    /* Normalize the normal */
    plane->normal[0] = nx;
    plane->normal[1] = ny;
    plane->normal[2] = nz;
    v3arr_normalize(plane->normal);

    /* Orthogonalize up vector to get V (vertical in image) */
    double up[3] = {ux, uy, uz};
    double dot = v3arr_dot(up, plane->normal);

    plane->v_axis[0] = up[0] - dot * plane->normal[0];
    plane->v_axis[1] = up[1] - dot * plane->normal[1];
    plane->v_axis[2] = up[2] - dot * plane->normal[2];
    v3arr_normalize(plane->v_axis);

    /* U = V × N (horizontal in image, right-handed) */
    v3arr_cross(plane->v_axis, plane->normal, plane->u_axis);
    v3arr_normalize(plane->u_axis);
}

void alea_slice_plane_init_axis(alea_slice_plane_t* plane, int axis, double value) {
    memset(plane, 0, sizeof(*plane));

    switch (axis) {
        case 0: /* X-axis: YZ plane at x=value */
            plane->origin[0] = value;
            plane->normal[0] = 1.0;
            plane->u_axis[1] = 1.0;  /* Y is horizontal */
            plane->v_axis[2] = 1.0;  /* Z is vertical */
            break;
        case 1: /* Y-axis: XZ plane at y=value */
            plane->origin[1] = value;
            plane->normal[1] = 1.0;
            plane->u_axis[0] = 1.0;  /* X is horizontal */
            plane->v_axis[2] = 1.0;  /* Z is vertical */
            break;
        case 2: /* Z-axis: XY plane at z=value */
        default:
            plane->origin[2] = value;
            plane->normal[2] = 1.0;
            plane->u_axis[0] = 1.0;  /* X is horizontal */
            plane->v_axis[1] = 1.0;  /* Y is vertical */
            break;
    }
}

void alea_plane_to_2d(const alea_slice_plane_t* plane,
                     double x, double y, double z,
                     double* u, double* v) {
    /* Vector from origin to point */
    double p[3] = {x - plane->origin[0],
                   y - plane->origin[1],
                   z - plane->origin[2]};

    /* Project onto basis vectors */
    *u = v3arr_dot(p, plane->u_axis);
    *v = v3arr_dot(p, plane->v_axis);
}

void alea_plane_to_3d(const alea_slice_plane_t* plane,
                     double u, double v,
                     double* x, double* y, double* z) {
    *x = plane->origin[0] + u * plane->u_axis[0] + v * plane->v_axis[0];
    *y = plane->origin[1] + u * plane->u_axis[1] + v * plane->v_axis[1];
    *z = plane->origin[2] + u * plane->u_axis[2] + v * plane->v_axis[2];
}

/* ============================================================================
 * PRIMITIVE-SPECIFIC INTERSECTION FUNCTIONS
 * ============================================================================ */

/**
 * @brief Intersect plane (surface) with slice plane
 *
 * Surface plane: a*x + b*y + c*z + d = 0
 * Slice plane: P(u,v) = O + u*U + v*V
 *
 * Substituting: a*(Ox + u*Ux + v*Vx) + b*(Oy + u*Uy + v*Vy) + c*(Oz + u*Uz + v*Vz) + d = 0
 * Rearranging: (a*Ux + b*Uy + c*Uz)*u + (a*Vx + b*Vy + c*Vz)*v + (a*Ox + b*Oy + c*Oz + d) = 0
 *
 * This is a 2D line: A*u + B*v + C = 0
 */
static bool intersect_plane(const alea_plane_data_t* surf,
                            const alea_slice_plane_t* slice,
                            alea_curve_2d_t* curve) {
    double A = surf->a * slice->u_axis[0] + surf->b * slice->u_axis[1] + surf->c * slice->u_axis[2];
    double B = surf->a * slice->v_axis[0] + surf->b * slice->v_axis[1] + surf->c * slice->v_axis[2];
    double C = surf->a * slice->origin[0] + surf->b * slice->origin[1] +
               surf->c * slice->origin[2] + surf->d;

    double len = sqrt(A*A + B*B);
    if (len < EPSILON) {
        /* Surface plane is parallel to slice plane */
        if (fabs(C) < EPSILON) {
            /* Coincident planes - the entire slice is on the surface */
            /* This is a degenerate case, we could return a special marker */
            return false;
        }
        /* No intersection */
        return false;
    }

    curve->type = ALEA_CURVE_LINE;
    curve->data.line.a = A;
    curve->data.line.b = B;
    curve->data.line.c = C;

    /* Direction vector perpendicular to (A, B) */
    curve->data.line.direction[0] = -B / len;
    curve->data.line.direction[1] = A / len;

    /* Find a point on the line */
    if (fabs(A) > fabs(B)) {
        curve->data.line.point[0] = -C / A;
        curve->data.line.point[1] = 0;
    } else {
        curve->data.line.point[0] = 0;
        curve->data.line.point[1] = -C / B;
    }

    return true;
}

/**
 * @brief Intersect sphere with slice plane
 *
 * Sphere: (x-cx)² + (y-cy)² + (z-cz)² = r²
 * Slice plane: P(u,v) = O + u*U + v*V
 *
 * The intersection is a circle if the plane passes through the sphere.
 * Substituting and simplifying gives a 2D circle.
 */
static bool intersect_sphere(const alea_sphere_data_t* sphere,
                             const alea_slice_plane_t* slice,
                             alea_curve_2d_t* curve) {
    /* Vector from sphere center to slice origin */
    double d[3] = {
        slice->origin[0] - sphere->center_x,
        slice->origin[1] - sphere->center_y,
        slice->origin[2] - sphere->center_z
    };

    /* Distance from sphere center to slice plane (signed) */
    double dist_to_plane = v3arr_dot(d, slice->normal);
    double abs_dist = fabs(dist_to_plane);

    if (abs_dist > sphere->radius - EPSILON) {
        /* Plane doesn't intersect sphere (or tangent) */
        if (fabs(abs_dist - sphere->radius) < EPSILON) {
            /* Tangent - single point */
            curve->type = ALEA_CURVE_POINT;
            /* Project center onto plane */
            alea_plane_to_2d(slice, sphere->center_x, sphere->center_y, sphere->center_z,
                           &curve->data.point[0], &curve->data.point[1]);
            return true;
        }
        return false;
    }

    /* Intersection is a circle */
    /* Radius of intersection circle: sqrt(R² - d²) */
    double r_intersect = sqrt(sphere->radius * sphere->radius - dist_to_plane * dist_to_plane);

    /* Center of intersection circle: project sphere center onto plane */
    double center_3d[3] = {
        sphere->center_x - dist_to_plane * slice->normal[0],
        sphere->center_y - dist_to_plane * slice->normal[1],
        sphere->center_z - dist_to_plane * slice->normal[2]
    };

    curve->type = ALEA_CURVE_CIRCLE;
    alea_plane_to_2d(slice, center_3d[0], center_3d[1], center_3d[2],
                   &curve->data.circle.center[0], &curve->data.circle.center[1]);
    curve->data.circle.radius = r_intersect;

    return true;
}

/**
 * @brief Intersect axis-aligned cylinder with slice plane (unified)
 *
 * Cylinder equation for axis a: p² + q² = r²  (p, q = perpendicular coords)
 *
 * @param center_p  Cylinder center along first perpendicular axis
 * @param center_q  Cylinder center along second perpendicular axis
 * @param radius    Cylinder radius
 * @param axis      Cylinder axis (0=X, 1=Y, 2=Z)
 */
static bool intersect_cylinder(double center_p, double center_q, double radius,
                               int axis,
                               const alea_slice_plane_t* slice,
                               alea_curve_2d_t* curve) {
    /* Perpendicular coordinate indices: axis=0→{1,2}, axis=1→{0,2}, axis=2→{0,1} */
    static const int perp[3][2] = {{1,2}, {0,2}, {0,1}};
    int pi = perp[axis][0];
    int qi = perp[axis][1];

    /* Check perpendicular case: slice normal parallel to cylinder axis */
    double cos_angle = fabs(slice->normal[axis]);

    if (cos_angle > 1.0 - EPSILON) {
        /* Perpendicular cut → circle */
        double center_3d[3];
        center_3d[axis] = slice->origin[axis];
        center_3d[pi] = center_p;
        center_3d[qi] = center_q;

        curve->type = ALEA_CURVE_CIRCLE;
        alea_plane_to_2d(slice, center_3d[0], center_3d[1], center_3d[2],
                       &curve->data.circle.center[0], &curve->data.circle.center[1]);
        curve->data.circle.radius = radius;
        return true;
    }

    /* General case: conic from substituting P(u,v) = O + u*U + v*V
     * into cylinder eq: (P_p - center_p)² + (P_q - center_q)² = r² */
    double Op = slice->origin[pi] - center_p;
    double Oq = slice->origin[qi] - center_q;
    double Up = slice->u_axis[pi];
    double Uq = slice->u_axis[qi];
    double Vp = slice->v_axis[pi];
    double Vq = slice->v_axis[qi];
    double r2 = radius * radius;

    double A = Up*Up + Uq*Uq;
    double B = 2*(Up*Vp + Uq*Vq);
    double C = Vp*Vp + Vq*Vq;
    double D = 2*(Op*Up + Oq*Uq);
    double E = 2*(Op*Vp + Oq*Vq);
    double F = Op*Op + Oq*Oq - r2;

    if (fabs(A) < EPSILON && fabs(C) < EPSILON) {
        /* Degenerate: plane contains cylinder axis */
        double dist_sq = Op*Op + Oq*Oq;
        if (dist_sq > r2 + EPSILON) return false;
        curve->type = ALEA_CURVE_HYPERBOLA;
        curve->data.conic.A = A;
        curve->data.conic.B = B;
        curve->data.conic.C = C;
        curve->data.conic.D = D;
        curve->data.conic.E = E;
        curve->data.conic.F = F;
        return true;
    }

    /* Check for parallel lines (plane parallel to cylinder axis) */
    if (fabs(B) < EPSILON) {
        if (fabs(C) < EPSILON && fabs(E) < EPSILON && fabs(A) > EPSILON) {
            /* Vertical parallel lines: A*u² + D*u + F = 0 */
            double disc = D*D - 4*A*F;
            if (disc < 0) return false;
            double sqrt_disc = sqrt(disc);
            double u1 = (-D - sqrt_disc) / (2*A);
            double u2 = (-D + sqrt_disc) / (2*A);
            curve->type = ALEA_CURVE_PARALLEL_LINES;
            curve->data.parallel_lines.point1[0] = u1;
            curve->data.parallel_lines.point1[1] = 0;
            curve->data.parallel_lines.point2[0] = u2;
            curve->data.parallel_lines.point2[1] = 0;
            curve->data.parallel_lines.direction[0] = 0;
            curve->data.parallel_lines.direction[1] = 1;
            return true;
        }
        if (fabs(A) < EPSILON && fabs(D) < EPSILON && fabs(C) > EPSILON) {
            /* Horizontal parallel lines: C*v² + E*v + F = 0 */
            double disc = E*E - 4*C*F;
            if (disc < 0) return false;
            double sqrt_disc = sqrt(disc);
            double v1 = (-E - sqrt_disc) / (2*C);
            double v2 = (-E + sqrt_disc) / (2*C);
            curve->type = ALEA_CURVE_PARALLEL_LINES;
            curve->data.parallel_lines.point1[0] = 0;
            curve->data.parallel_lines.point1[1] = v1;
            curve->data.parallel_lines.point2[0] = 0;
            curve->data.parallel_lines.point2[1] = v2;
            curve->data.parallel_lines.direction[0] = 1;
            curve->data.parallel_lines.direction[1] = 0;
            return true;
        }
    }

    curve->type = ALEA_CURVE_ELLIPSE;
    curve->data.conic.A = A;
    curve->data.conic.B = B;
    curve->data.conic.C = C;
    curve->data.conic.D = D;
    curve->data.conic.E = E;
    curve->data.conic.F = F;
    finalize_conic_as_ellipse(curve);

    return true;
}

/**
 * @brief Intersect axis-aligned cone with slice plane (unified)
 *
 * Cone equation for axis a: p² + q² = k²·a²  (p, q = perpendicular coords)
 *
 * @param apex_x, apex_y, apex_z  Cone apex
 * @param k2       tan²(half-angle)
 * @param axis     Cone axis (0=X, 1=Y, 2=Z)
 */
static bool intersect_cone(double apex_x, double apex_y, double apex_z,
                           double k2, int axis,
                           const alea_slice_plane_t* slice,
                           alea_curve_2d_t* curve) {
    /* Perpendicular coordinate indices */
    static const int perp[3][2] = {{1,2}, {0,2}, {0,1}};
    int pi = perp[axis][0];
    int qi = perp[axis][1];

    double apex[3] = {apex_x, apex_y, apex_z};

    /* Perpendicular case: slice normal parallel to cone axis → circle */
    double cos_angle = fabs(slice->normal[axis]);
    if (cos_angle > 1.0 - EPSILON) {
        /* Distance from apex to slice plane along axis */
        double d_axis = slice->origin[axis] - apex[axis];
        if (fabs(d_axis) < EPSILON) {
            /* Plane through apex — single point */
            curve->type = ALEA_CURVE_POINT;
            alea_plane_to_2d(slice, apex_x, apex_y, apex_z,
                           &curve->data.point[0], &curve->data.point[1]);
            return true;
        }
        double radius = fabs(d_axis) * sqrt(k2);
        double center_3d[3];
        center_3d[axis] = slice->origin[axis];
        center_3d[pi] = apex[pi];
        center_3d[qi] = apex[qi];

        curve->type = ALEA_CURVE_CIRCLE;
        alea_plane_to_2d(slice, center_3d[0], center_3d[1], center_3d[2],
                       &curve->data.circle.center[0], &curve->data.circle.center[1]);
        curve->data.circle.radius = radius;
        return true;
    }

    double O[3] = {slice->origin[0] - apex[0],
                   slice->origin[1] - apex[1],
                   slice->origin[2] - apex[2]};

    /* Cone: (O_p + u*U_p + v*V_p)² + (O_q + u*U_q + v*V_q)²
     *     - k²*(O_a + u*U_a + v*V_a)² = 0 */
    double Up = slice->u_axis[pi], Uq = slice->u_axis[qi], Ua = slice->u_axis[axis];
    double Vp = slice->v_axis[pi], Vq = slice->v_axis[qi], Va = slice->v_axis[axis];
    double Op = O[pi], Oq = O[qi], Oa = O[axis];

    double A = Up*Up + Uq*Uq - k2*Ua*Ua;
    double B = 2*(Up*Vp + Uq*Vq - k2*Ua*Va);
    double C = Vp*Vp + Vq*Vq - k2*Va*Va;
    double D = 2*(Op*Up + Oq*Uq - k2*Oa*Ua);
    double E = 2*(Op*Vp + Oq*Vq - k2*Oa*Va);
    double F = Op*Op + Oq*Oq - k2*Oa*Oa;

    double disc = B*B - 4*A*C;
    curve->type = (disc < -EPSILON) ? ALEA_CURVE_ELLIPSE :
                  (disc > EPSILON)  ? ALEA_CURVE_HYPERBOLA : ALEA_CURVE_PARABOLA;

    curve->data.conic.A = A;
    curve->data.conic.B = B;
    curve->data.conic.C = C;
    curve->data.conic.D = D;
    curve->data.conic.E = E;
    curve->data.conic.F = F;
    finalize_conic_as_ellipse(curve);

    return true;
}

/**
 * @brief Intersect axis-aligned box with slice plane
 *
 * Returns a polygon (up to 6 vertices) formed by the intersection.
 */
static bool intersect_box(const alea_box_data_t* box,
                          const alea_slice_plane_t* slice,
                          alea_curve_2d_t* curve) {
    /* 12 edges of the box */
    double edges[12][2][3] = {
        /* Bottom face (z = min_z) */
        {{box->min_x, box->min_y, box->min_z}, {box->max_x, box->min_y, box->min_z}},
        {{box->max_x, box->min_y, box->min_z}, {box->max_x, box->max_y, box->min_z}},
        {{box->max_x, box->max_y, box->min_z}, {box->min_x, box->max_y, box->min_z}},
        {{box->min_x, box->max_y, box->min_z}, {box->min_x, box->min_y, box->min_z}},
        /* Top face (z = max_z) */
        {{box->min_x, box->min_y, box->max_z}, {box->max_x, box->min_y, box->max_z}},
        {{box->max_x, box->min_y, box->max_z}, {box->max_x, box->max_y, box->max_z}},
        {{box->max_x, box->max_y, box->max_z}, {box->min_x, box->max_y, box->max_z}},
        {{box->min_x, box->max_y, box->max_z}, {box->min_x, box->min_y, box->max_z}},
        /* Vertical edges */
        {{box->min_x, box->min_y, box->min_z}, {box->min_x, box->min_y, box->max_z}},
        {{box->max_x, box->min_y, box->min_z}, {box->max_x, box->min_y, box->max_z}},
        {{box->max_x, box->max_y, box->min_z}, {box->max_x, box->max_y, box->max_z}},
        {{box->min_x, box->max_y, box->min_z}, {box->min_x, box->max_y, box->max_z}},
    };

    double intersections[12][2];  /* (u, v) for each valid intersection */
    int count = 0;

    /* Test each edge against the slice plane */
    for (int i = 0; i < 12 && count < 12; i++) {
        double* p0 = edges[i][0];
        double* p1 = edges[i][1];

        /* Edge direction */
        double d[3] = {p1[0] - p0[0], p1[1] - p0[1], p1[2] - p0[2]};

        /* Vector from plane origin to edge start */
        double w[3] = {p0[0] - slice->origin[0],
                       p0[1] - slice->origin[1],
                       p0[2] - slice->origin[2]};

        /* Signed distance of endpoints from plane */
        double d0 = v3arr_dot(w, slice->normal);
        double dd = v3arr_dot(d, slice->normal);

        if (fabs(dd) < EPSILON) {
            /* Edge parallel to plane */
            continue;
        }

        double t = -d0 / dd;
        if (t < -EPSILON || t > 1.0 + EPSILON) {
            /* Intersection outside edge */
            continue;
        }

        /* Intersection point in 3D */
        double px = p0[0] + t * d[0];
        double py = p0[1] + t * d[1];
        double pz = p0[2] + t * d[2];

        /* Convert to 2D */
        alea_plane_to_2d(slice, px, py, pz,
                       &intersections[count][0], &intersections[count][1]);
        count++;
    }

    if (count < 3) {
        return false;  /* Not enough points for a polygon */
    }

    /* Sort intersection points by angle around centroid */
    double cx = 0, cy = 0;
    for (int i = 0; i < count; i++) {
        cx += intersections[i][0];
        cy += intersections[i][1];
    }
    cx /= count;
    cy /= count;

    /* Pre-compute angles, then insertion sort (max 12 points) */
    double angles[12];
    for (int i = 0; i < count; i++)
        angles[i] = atan2(intersections[i][1] - cy, intersections[i][0] - cx);

    for (int i = 1; i < count; i++) {
        double ka = angles[i];
        double ku = intersections[i][0], kv = intersections[i][1];
        int j = i - 1;
        while (j >= 0 && angles[j] > ka) {
            angles[j + 1] = angles[j];
            intersections[j + 1][0] = intersections[j][0];
            intersections[j + 1][1] = intersections[j][1];
            j--;
        }
        angles[j + 1] = ka;
        intersections[j + 1][0] = ku;
        intersections[j + 1][1] = kv;
    }

    curve->type = ALEA_CURVE_POLYGON;
    curve->data.polygon.vertex_count = count;
    curve->data.polygon.closed = true;
    for (int i = 0; i < count && i < ALEA_MAX_POLYGON_VERTICES; i++) {
        curve->data.polygon.vertices[i][0] = intersections[i][0];
        curve->data.polygon.vertices[i][1] = intersections[i][1];
    }

    return true;
}

/**
 * @brief Intersect general quadric with slice plane
 *
 * Quadric: A*x² + B*y² + C*z² + D*xy + E*yz + F*xz + G*x + H*y + I*z + J = 0
 * Result is a general conic section.
 */
static bool intersect_quadric(const alea_quadric_data_t* quad,
                              const alea_slice_plane_t* slice,
                              alea_curve_2d_t* curve) {
    double A = quad->coeffs[0];
    double B = quad->coeffs[1];
    double C = quad->coeffs[2];
    double D = quad->coeffs[3];
    double E = quad->coeffs[4];
    double F = quad->coeffs[5];
    double G = quad->coeffs[6];
    double H = quad->coeffs[7];
    double I = quad->coeffs[8];
    double J = quad->coeffs[9];

    double Ox = slice->origin[0];
    double Oy = slice->origin[1];
    double Oz = slice->origin[2];
    double Ux = slice->u_axis[0];
    double Uy = slice->u_axis[1];
    double Uz = slice->u_axis[2];
    double Vx = slice->v_axis[0];
    double Vy = slice->v_axis[1];
    double Vz = slice->v_axis[2];

    /* Substitute x = Ox + u*Ux + v*Vx, etc. into quadric equation */
    /* After much algebra, we get: A'*u² + B'*uv + C'*v² + D'*u + E'*v + F' = 0 */

    double Ap = A*Ux*Ux + B*Uy*Uy + C*Uz*Uz + D*Ux*Uy + E*Uy*Uz + F*Ux*Uz;
    double Bp = 2*A*Ux*Vx + 2*B*Uy*Vy + 2*C*Uz*Vz + D*(Ux*Vy + Uy*Vx) +
                E*(Uy*Vz + Uz*Vy) + F*(Ux*Vz + Uz*Vx);
    double Cp = A*Vx*Vx + B*Vy*Vy + C*Vz*Vz + D*Vx*Vy + E*Vy*Vz + F*Vx*Vz;
    double Dp = 2*A*Ox*Ux + 2*B*Oy*Uy + 2*C*Oz*Uz + D*(Ox*Uy + Oy*Ux) +
                E*(Oy*Uz + Oz*Uy) + F*(Ox*Uz + Oz*Ux) + G*Ux + H*Uy + I*Uz;
    double Ep = 2*A*Ox*Vx + 2*B*Oy*Vy + 2*C*Oz*Vz + D*(Ox*Vy + Oy*Vx) +
                E*(Oy*Vz + Oz*Vy) + F*(Ox*Vz + Oz*Vx) + G*Vx + H*Vy + I*Vz;
    double Fp = A*Ox*Ox + B*Oy*Oy + C*Oz*Oz + D*Ox*Oy + E*Oy*Oz + F*Ox*Oz +
                G*Ox + H*Oy + I*Oz + J;

    /* Classify and check for degenerate cases */
    double disc = Bp*Bp - 4*Ap*Cp;

    /* Check for parallel lines: disc ≈ 0 means quadratic form is a perfect square.
     * The conic factors as (a*u + b*v + c1)(a*u + b*v + c2) = 0
     * where a² = Ap, b² = Cp, 2ab = Bp.
     */
    if (fabs(disc) < EPSILON) {
        /* Degenerate case - could be parallel lines, single line, or parabola */
        double a = 0, b = 0;

        if (fabs(Ap) > EPSILON) {
            a = sqrt(fabs(Ap));
            if (fabs(Bp) > EPSILON) {
                /* Choose sign of b so that 2ab = Bp */
                b = Bp / (2 * a);
            } else if (fabs(Cp) > EPSILON) {
                b = sqrt(fabs(Cp));
            }
        } else if (fabs(Cp) > EPSILON) {
            b = sqrt(fabs(Cp));
            if (fabs(Bp) > EPSILON) {
                a = Bp / (2 * b);
            }
        }

        double norm_sq = a*a + b*b;
        if (norm_sq > EPSILON) {
            /* Normal direction: (a, b), line direction: (-b, a) */
            double norm = sqrt(norm_sq);
            double nx = a / norm;
            double ny = b / norm;

            /* Project linear terms onto normal and tangent directions */
            /* Linear part: Dp*u + Ep*v = (Dp*nx + Ep*ny)*w + (Dp*(-ny) + Ep*nx)*t
             * where w = nx*u + ny*v (normal dir) and t = -ny*u + nx*v (tangent dir) */
            double lin_tangent = -Dp * ny + Ep * nx; /* Coefficient in tangent direction */

            /* For parallel lines, tangent coefficient must be ~0 */
            if (fabs(lin_tangent) < EPSILON * (fabs(Dp) + fabs(Ep) + 1)) {
                /* Equation in normal direction: w² + lin_normal*w + Fp = 0
                 * where w = a*u + b*v (unnormalized) = norm * (nx*u + ny*v)
                 * So: (norm*w')² + lin_normal*(norm*w'/norm)*norm + Fp = 0...
                 * Actually, let's work with the original: (a*u + b*v)² + Dp*u + Ep*v + Fp = 0
                 * Let w = a*u + b*v. Then u = (w*a - t*b)/norm_sq, v = (w*b + t*a)/norm_sq
                 * where t = -b*u + a*v (perpendicular component).
                 * Dp*u + Ep*v = (Dp*a + Ep*b)*w/norm_sq + (-Dp*b + Ep*a)*t/norm_sq
                 * For parallel lines: -Dp*b + Ep*a ≈ 0 (which is lin_tangent*norm)
                 * Equation becomes: w² + (Dp*a + Ep*b)/norm_sq * w + Fp = 0 */
                double coef = (Dp * a + Ep * b) / norm_sq;
                double w_disc = coef * coef - 4 * Fp;

                if (w_disc >= -EPSILON) {
                    /* Two parallel lines (or one line if w_disc ≈ 0) */
                    double sqrt_w_disc = (w_disc > 0) ? sqrt(w_disc) : 0;
                    double w1 = (-coef - sqrt_w_disc) / 2;
                    double w2 = (-coef + sqrt_w_disc) / 2;

                    /* Lines: a*u + b*v = w1 and a*u + b*v = w2 */
                    /* Direction of lines: perpendicular to (a,b), i.e., (-b, a) normalized */
                    double dir_len = norm;
                    double dir_u = -b / dir_len;
                    double dir_v = a / dir_len;

                    /* Points on lines: solve a*u + b*v = w for a point */
                    double p1_u, p1_v, p2_u, p2_v;
                    if (fabs(a) > fabs(b)) {
                        p1_u = w1 / a; p1_v = 0;
                        p2_u = w2 / a; p2_v = 0;
                    } else {
                        p1_u = 0; p1_v = w1 / b;
                        p2_u = 0; p2_v = w2 / b;
                    }

                    curve->type = ALEA_CURVE_PARALLEL_LINES;
                    curve->data.parallel_lines.point1[0] = p1_u;
                    curve->data.parallel_lines.point1[1] = p1_v;
                    curve->data.parallel_lines.point2[0] = p2_u;
                    curve->data.parallel_lines.point2[1] = p2_v;
                    curve->data.parallel_lines.direction[0] = dir_u;
                    curve->data.parallel_lines.direction[1] = dir_v;
                    return true;
                }
            }
        }

        /* Fall through to parabola if not parallel lines */
        curve->type = ALEA_CURVE_PARABOLA;
    } else if (disc < -EPSILON) {
        curve->type = ALEA_CURVE_ELLIPSE;
    } else {
        curve->type = ALEA_CURVE_HYPERBOLA;
    }

    curve->data.conic.A = Ap;
    curve->data.conic.B = Bp;
    curve->data.conic.C = Cp;
    curve->data.conic.D = Dp;
    curve->data.conic.E = Ep;
    curve->data.conic.F = Fp;
    finalize_conic_as_ellipse(curve);

    return true;
}

/**
 * @brief Intersect RCC (Right Circular Cylinder macrobody) with slice plane
 *
 * RCC is a finite cylinder with arbitrary orientation.
 * The intersection can be an ellipse, circle, or bounded ellipse arc.
 *
 * Method: Transform the problem so cylinder axis is along Z, compute intersection,
 * then transform the result back to slice plane coordinates.
 */
static bool intersect_rcc(const alea_rcc_data_t* rcc,
                          const alea_slice_plane_t* slice,
                          alea_curve_2d_t* curve) {
    /* Cylinder axis direction */
    double axis[3] = {rcc->height_x, rcc->height_y, rcc->height_z};
    double h_len = v3arr_length(axis);
    if (h_len < EPSILON) return false;

    double axis_n[3] = {axis[0]/h_len, axis[1]/h_len, axis[2]/h_len};

    /* Check angle between slice normal and cylinder axis */
    double cos_angle = v3arr_dot(slice->normal, axis_n);
    double abs_cos = fabs(cos_angle);

    /* Vector from cylinder base to slice origin */
    double d_base[3] = {
        slice->origin[0] - rcc->base_x,
        slice->origin[1] - rcc->base_y,
        slice->origin[2] - rcc->base_z
    };

    if (abs_cos > 1.0 - EPSILON) {
        /* Perpendicular case: slice normal parallel to axis → circle */
        double da = v3arr_dot(axis_n, slice->normal);
        if (fabs(da) < EPSILON) return false;

        /* Find where axis intersects slice plane */
        double dp = v3arr_dot(d_base, slice->normal);
        double t = dp / da;

        /* Check if intersection is within cylinder bounds */
        if (t < -EPSILON || t > h_len + EPSILON) return false;

        double center_3d[3] = {
            rcc->base_x + t * axis_n[0],
            rcc->base_y + t * axis_n[1],
            rcc->base_z + t * axis_n[2]
        };

        curve->type = ALEA_CURVE_CIRCLE;
        alea_plane_to_2d(slice, center_3d[0], center_3d[1], center_3d[2],
                       &curve->data.circle.center[0], &curve->data.circle.center[1]);
        curve->data.circle.radius = rcc->radius;
        return true;
    }

    if (abs_cos < EPSILON) {
        /* Parallel case: slice normal perpendicular to axis */
        /* Check if plane intersects cylinder (distance from axis to plane < r) */
        double dist_to_plane = v3arr_dot(d_base, slice->normal);

        if (fabs(dist_to_plane) > rcc->radius + EPSILON) {
            return false;  /* Plane doesn't intersect cylinder */
        }

        /* Intersection is two parallel lines (or one line if tangent) */
        /* Project axis onto slice plane to get line direction */
        double dir_3d[3] = {
            axis_n[0] - cos_angle * slice->normal[0],
            axis_n[1] - cos_angle * slice->normal[1],
            axis_n[2] - cos_angle * slice->normal[2]
        };
        v3arr_normalize(dir_3d);

        /* Project base onto plane */
        double base_proj[3] = {
            rcc->base_x - dist_to_plane * slice->normal[0],
            rcc->base_y - dist_to_plane * slice->normal[1],
            rcc->base_z - dist_to_plane * slice->normal[2]
        };

        double base_2d[2], dir_2d[2];
        alea_plane_to_2d(slice, base_proj[0], base_proj[1], base_proj[2],
                       &base_2d[0], &base_2d[1]);

        /* Direction in 2D */
        double end_3d[3] = {base_proj[0] + dir_3d[0], base_proj[1] + dir_3d[1], base_proj[2] + dir_3d[2]};
        double end_2d[2];
        alea_plane_to_2d(slice, end_3d[0], end_3d[1], end_3d[2], &end_2d[0], &end_2d[1]);
        dir_2d[0] = end_2d[0] - base_2d[0];
        dir_2d[1] = end_2d[1] - base_2d[1];
        double dir_len = sqrt(dir_2d[0]*dir_2d[0] + dir_2d[1]*dir_2d[1]);
        if (dir_len > EPSILON) {
            dir_2d[0] /= dir_len;
            dir_2d[1] /= dir_len;
        }

        /* Perpendicular direction for offset */
        double perp_2d[2] = {-dir_2d[1], dir_2d[0]};

        /* Offset = sqrt(r² - d²) where d is distance from plane to axis */
        double offset = sqrt(rcc->radius * rcc->radius - dist_to_plane * dist_to_plane);

        if (offset < EPSILON) {
            /* Tangent — single line */
            curve->type = ALEA_CURVE_LINE;
            curve->data.line.point[0] = base_2d[0];
            curve->data.line.point[1] = base_2d[1];
            curve->data.line.direction[0] = dir_2d[0];
            curve->data.line.direction[1] = dir_2d[1];
            curve->data.line.a = -dir_2d[1];
            curve->data.line.b = dir_2d[0];
            curve->data.line.c = -(curve->data.line.a * base_2d[0] +
                                   curve->data.line.b * base_2d[1]);
            return true;
        }

        curve->type = ALEA_CURVE_PARALLEL_LINES;
        curve->data.parallel_lines.point1[0] = base_2d[0] + offset * perp_2d[0];
        curve->data.parallel_lines.point1[1] = base_2d[1] + offset * perp_2d[1];
        curve->data.parallel_lines.point2[0] = base_2d[0] - offset * perp_2d[0];
        curve->data.parallel_lines.point2[1] = base_2d[1] - offset * perp_2d[1];
        curve->data.parallel_lines.direction[0] = dir_2d[0];
        curve->data.parallel_lines.direction[1] = dir_2d[1];

        return true;
    }

    /* General oblique case: intersection is an ellipse */
    /* Build local coordinate system with cylinder axis as Z */
    double local_z[3] = {axis_n[0], axis_n[1], axis_n[2]};

    /* Pick an arbitrary perpendicular for local X */
    double local_x[3];
    if (fabs(local_z[0]) < 0.9) {
        local_x[0] = 1; local_x[1] = 0; local_x[2] = 0;
    } else {
        local_x[0] = 0; local_x[1] = 1; local_x[2] = 0;
    }
    /* Orthogonalize: local_x = local_x - (local_x·local_z)*local_z */
    double dot_xz = v3arr_dot(local_x, local_z);
    local_x[0] -= dot_xz * local_z[0];
    local_x[1] -= dot_xz * local_z[1];
    local_x[2] -= dot_xz * local_z[2];
    v3arr_normalize(local_x);

    /* local_y = local_z × local_x */
    double local_y[3];
    v3arr_cross(local_z, local_x, local_y);

    /* In local coords, cylinder is x² + y² = r² for z in [0, h] */
    /* Slice plane in local coords: transform origin and normal */

    /* Transform slice origin to local coords (relative to cylinder base) */
    double local_origin[3] = {
        v3arr_dot(d_base, local_x),
        v3arr_dot(d_base, local_y),
        v3arr_dot(d_base, local_z)
    };

    /* Transform slice basis vectors to local coords */
    double local_u[3] = {
        v3arr_dot(slice->u_axis, local_x),
        v3arr_dot(slice->u_axis, local_y),
        v3arr_dot(slice->u_axis, local_z)
    };
    double local_v[3] = {
        v3arr_dot(slice->v_axis, local_x),
        v3arr_dot(slice->v_axis, local_y),
        v3arr_dot(slice->v_axis, local_z)
    };

    /* Cylinder: x² + y² = r² */
    /* Slice plane: P(u,v) = O + u*U + v*V (in local coords) */
    /* Substitute: (Ox + u*Ux + v*Vx)² + (Oy + u*Uy + v*Vy)² = r² */

    double Ox = local_origin[0];
    double Oy = local_origin[1];
    double Ux = local_u[0];
    double Uy = local_u[1];
    double Vx = local_v[0];
    double Vy = local_v[1];
    double r2 = rcc->radius * rcc->radius;

    /* Expand to get conic coefficients */
    double A = Ux*Ux + Uy*Uy;
    double B = 2*(Ux*Vx + Uy*Vy);
    double C = Vx*Vx + Vy*Vy;
    double D = 2*(Ox*Ux + Oy*Uy);
    double E = 2*(Ox*Vx + Oy*Vy);
    double F = Ox*Ox + Oy*Oy - r2;

    /* Check discriminant */
    double disc = B*B - 4*A*C;

    if (disc >= -EPSILON) {
        /* Degenerate or hyperbolic - shouldn't happen for cylinder */
        curve->type = ALEA_CURVE_HYPERBOLA;
    } else {
        curve->type = ALEA_CURVE_ELLIPSE;
    }

    curve->data.conic.A = A;
    curve->data.conic.B = B;
    curve->data.conic.C = C;
    curve->data.conic.D = D;
    curve->data.conic.E = E;
    curve->data.conic.F = F;
    finalize_conic_as_ellipse(curve);

    return true;
}

/**
 * @brief Intersect TRC (Truncated Right Cone) with slice plane
 *
 * TRC: frustum from base (radius r1) to top (radius r2), arbitrary axis.
 * Lateral surface equation in local coords (axis=Z, base at origin):
 *   x² + y² = (r1 + z·(r2-r1)/h)²
 *
 * Three cases:
 * 1. Perpendicular (slice normal ∥ axis) → circle with radius r(z)
 * 2. Parallel (slice normal ⊥ axis) → general conic
 * 3. Oblique → general conic (ellipse/parabola/hyperbola)
 */
static bool intersect_trc(const alea_trc_data_t* trc,
                          const alea_slice_plane_t* slice,
                          alea_curve_2d_t* curve) {
    /* Axis direction and length */
    double axis[3] = {trc->height_x, trc->height_y, trc->height_z};
    double h_len = v3arr_length(axis);
    if (h_len < EPSILON) return false;

    double axis_n[3] = {axis[0]/h_len, axis[1]/h_len, axis[2]/h_len};
    double r1 = trc->base_radius;
    double r2 = trc->top_radius;
    double dr = (r2 - r1) / h_len;  /* radius slope per unit length */

    /* Angle between slice normal and TRC axis */
    double cos_angle = v3arr_dot(slice->normal, axis_n);
    double abs_cos = fabs(cos_angle);

    /* Vector from TRC base to slice origin */
    double d_base[3] = {
        slice->origin[0] - trc->base_x,
        slice->origin[1] - trc->base_y,
        slice->origin[2] - trc->base_z
    };

    if (abs_cos > 1.0 - EPSILON) {
        /* Perpendicular case: slice normal parallel to axis → circle */
        double da = v3arr_dot(axis_n, slice->normal);
        if (fabs(da) < EPSILON) return false;

        /* Find where axis intersects slice plane */
        double dp = v3arr_dot(d_base, slice->normal);
        double t = dp / da;

        /* Check if intersection is within TRC bounds */
        if (t < -EPSILON || t > h_len + EPSILON) return false;

        /* Radius at this height */
        double t_clamped = (t < 0) ? 0 : (t > h_len ? h_len : t);
        double radius = r1 + t_clamped * dr;
        if (radius < EPSILON) return false;

        double center_3d[3] = {
            trc->base_x + t * axis_n[0],
            trc->base_y + t * axis_n[1],
            trc->base_z + t * axis_n[2]
        };

        curve->type = ALEA_CURVE_CIRCLE;
        alea_plane_to_2d(slice, center_3d[0], center_3d[1], center_3d[2],
                       &curve->data.circle.center[0], &curve->data.circle.center[1]);
        curve->data.circle.radius = radius;
        return true;
    }

    /* General case (including parallel): build local coordinate system */
    double local_z[3] = {axis_n[0], axis_n[1], axis_n[2]};

    double local_x[3];
    if (fabs(local_z[0]) < 0.9) {
        local_x[0] = 1; local_x[1] = 0; local_x[2] = 0;
    } else {
        local_x[0] = 0; local_x[1] = 1; local_x[2] = 0;
    }
    double dot_xz = v3arr_dot(local_x, local_z);
    local_x[0] -= dot_xz * local_z[0];
    local_x[1] -= dot_xz * local_z[1];
    local_x[2] -= dot_xz * local_z[2];
    v3arr_normalize(local_x);

    double local_y[3];
    v3arr_cross(local_z, local_x, local_y);

    /* Transform to local coords (relative to TRC base) */
    double local_origin[3] = {
        v3arr_dot(d_base, local_x),
        v3arr_dot(d_base, local_y),
        v3arr_dot(d_base, local_z)
    };
    double local_u[3] = {
        v3arr_dot(slice->u_axis, local_x),
        v3arr_dot(slice->u_axis, local_y),
        v3arr_dot(slice->u_axis, local_z)
    };
    double local_v[3] = {
        v3arr_dot(slice->v_axis, local_x),
        v3arr_dot(slice->v_axis, local_y),
        v3arr_dot(slice->v_axis, local_z)
    };

    /* TRC lateral surface: x² + y² = (r1 + z·dr)²
     * Slice plane: P(u,v) = O + u*U + v*V (in local coords)
     *
     * Substitute: (Ox+u·Ux+v·Vx)² + (Oy+u·Uy+v·Vy)² = (R0 + u·Ru + v·Rv)²
     * where R0 = r1 + Oz·dr, Ru = Uz·dr, Rv = Vz·dr
     *
     * Expanding Left - Right = 0 gives conic: Au² + Buv + Cv² + Du + Ev + F = 0
     */
    double Ox = local_origin[0], Oy = local_origin[1], Oz = local_origin[2];
    double Ux = local_u[0], Uy = local_u[1], Uz = local_u[2];
    double Vx = local_v[0], Vy = local_v[1], Vz = local_v[2];

    double R0 = r1 + Oz * dr;
    double Ru = Uz * dr;
    double Rv = Vz * dr;

    double A = Ux*Ux + Uy*Uy - Ru*Ru;
    double B = 2*(Ux*Vx + Uy*Vy - Ru*Rv);
    double C = Vx*Vx + Vy*Vy - Rv*Rv;
    double D = 2*(Ox*Ux + Oy*Uy - R0*Ru);
    double E = 2*(Ox*Vx + Oy*Vy - R0*Rv);
    double F = Ox*Ox + Oy*Oy - R0*R0;

    /* Classify conic */
    double disc = B*B - 4*A*C;

    if (disc < -EPSILON) {
        curve->type = ALEA_CURVE_ELLIPSE;
    } else if (disc > EPSILON) {
        curve->type = ALEA_CURVE_HYPERBOLA;
    } else {
        curve->type = ALEA_CURVE_PARABOLA;
    }

    curve->data.conic.A = A;
    curve->data.conic.B = B;
    curve->data.conic.C = C;
    curve->data.conic.D = D;
    curve->data.conic.E = E;
    curve->data.conic.F = F;
    finalize_conic_as_ellipse(curve);

    return true;
}

/**
 * @brief Intersect torus with slice plane
 *
 * Torus equation (axis along Z):
 *   (sqrt(x² + y²) - R)² + z² = r²
 * where R = major radius, r = minor radius
 *
 * Intersection with arbitrary plane is a quartic curve.
 * For axis-aligned planes through the center, we get circles.
 *
 * We handle special cases analytically and fall back to quartic for general case.
 */
static bool intersect_torus(const alea_torus_data_t* torus,
                            alea_primitive_type_t prim_type,
                            const alea_slice_plane_t* slice,
                            alea_curve_2d_t* curve) {
    double R = torus->major_radius;
    double r = torus->minor_radius;
    double cx = torus->center_x;
    double cy = torus->center_y;
    double cz = torus->center_z;

    /* Determine torus axis based on primitive type */
    double torus_axis[3] = {0, 0, 0};
    switch (prim_type) {
        case ALEA_PRIMITIVE_TORUS_X: torus_axis[0] = 1; break;
        case ALEA_PRIMITIVE_TORUS_Y: torus_axis[1] = 1; break;
        case ALEA_PRIMITIVE_TORUS_Z: torus_axis[2] = 1; break;
        default: return false;
    }

    /* Angle between slice normal and torus axis */
    double cos_angle = v3arr_dot(slice->normal, torus_axis);
    double abs_cos = fabs(cos_angle);

    /* Distance from torus center to slice plane */
    double d_center[3] = {
        slice->origin[0] - cx,
        slice->origin[1] - cy,
        slice->origin[2] - cz
    };
    double dist_to_plane = v3arr_dot(d_center, slice->normal);

    /* Special case 1: Slice perpendicular to torus axis (through center region) */
    if (abs_cos > 1.0 - EPSILON) {
        /* Slice plane perpendicular to torus axis */
        double d = fabs(dist_to_plane);  /* Distance from torus center to plane */

        if (d > r + EPSILON) {
            return false;  /* Plane doesn't intersect torus */
        }

        /* Intersection is two concentric circles (or one if tangent) */
        /* Radii: R ± sqrt(r² - d²) */
        double r_offset = sqrt(r*r - d*d);
        double r_inner = R - r_offset;
        double r_outer = R + r_offset;

        /* Project center onto plane */
        double center_3d[3] = {
            cx - dist_to_plane * slice->normal[0],
            cy - dist_to_plane * slice->normal[1],
            cz - dist_to_plane * slice->normal[2]
        };
        double center_2d[2];
        alea_plane_to_2d(slice, center_3d[0], center_3d[1], center_3d[2],
                       &center_2d[0], &center_2d[1]);

        if (r_inner > EPSILON) {
            /* Two concentric circles */
            curve->type = ALEA_CURVE_QUARTIC;
            curve->data.torus.mode = 1;  /* concentric circles */
            curve->data.torus.c1u = center_2d[0];
            curve->data.torus.c1v = center_2d[1];
            curve->data.torus.r1 = r_inner;
            curve->data.torus.c2u = center_2d[0];
            curve->data.torus.c2v = center_2d[1];
            curve->data.torus.r2 = r_outer;
        } else {
            /* Single circle (tangent case) */
            curve->type = ALEA_CURVE_CIRCLE;
            curve->data.circle.center[0] = center_2d[0];
            curve->data.circle.center[1] = center_2d[1];
            curve->data.circle.radius = r_outer;
        }
        return true;
    }

    /* Special case 2: Slice parallel to torus axis, passing through center */
    if (abs_cos < EPSILON && fabs(dist_to_plane) < EPSILON) {
        /* Slice through the center, parallel to axis */
        /* This produces two circles of radius r, centered at distance R from center */

        /* Find the two circle centers: project ±R along the radial direction in the plane */
        double radial[3];
        v3arr_cross(torus_axis, slice->normal, radial);
        double rad_len = v3arr_length(radial);
        if (rad_len < EPSILON) {
            return false;
        }
        radial[0] /= rad_len;
        radial[1] /= rad_len;
        radial[2] /= rad_len;

        /* Two circle centers at cx ± R*radial */
        double c1[3] = {cx + R*radial[0], cy + R*radial[1], cz + R*radial[2]};
        double c2[3] = {cx - R*radial[0], cy - R*radial[1], cz - R*radial[2]};

        double c1_2d[2], c2_2d[2];
        alea_plane_to_2d(slice, c1[0], c1[1], c1[2], &c1_2d[0], &c1_2d[1]);
        alea_plane_to_2d(slice, c2[0], c2[1], c2[2], &c2_2d[0], &c2_2d[1]);

        curve->type = ALEA_CURVE_QUARTIC;
        curve->data.torus.mode = 2;  /* two separate circles */
        curve->data.torus.c1u = c1_2d[0];
        curve->data.torus.c1v = c1_2d[1];
        curve->data.torus.r1 = r;
        curve->data.torus.c2u = c2_2d[0];
        curve->data.torus.c2v = c2_2d[1];
        curve->data.torus.r2 = r;
        return true;
    }

    /* General case: store parameters for quartic evaluation during scanline */
    /* Transform to local coordinates where torus axis is Z */
    double local_z[3] = {torus_axis[0], torus_axis[1], torus_axis[2]};

    /* Build local coordinate system */
    double local_x[3];
    if (fabs(local_z[0]) < 0.9) {
        local_x[0] = 1; local_x[1] = 0; local_x[2] = 0;
    } else {
        local_x[0] = 0; local_x[1] = 1; local_x[2] = 0;
    }
    double dot_xz = v3arr_dot(local_x, local_z);
    local_x[0] -= dot_xz * local_z[0];
    local_x[1] -= dot_xz * local_z[1];
    local_x[2] -= dot_xz * local_z[2];
    v3arr_normalize(local_x);

    double local_y[3];
    v3arr_cross(local_z, local_x, local_y);

    /* Transform slice plane to local coords (relative to torus center) */
    curve->type = ALEA_CURVE_QUARTIC;
    curve->data.torus.mode = 0;  /* general quartic */
    curve->data.torus.R = R;
    curve->data.torus.r = r;
    curve->data.torus.Ox = v3arr_dot(d_center, local_x);
    curve->data.torus.Oy = v3arr_dot(d_center, local_y);
    curve->data.torus.Oz = v3arr_dot(d_center, local_z);
    curve->data.torus.Ux = v3arr_dot(slice->u_axis, local_x);
    curve->data.torus.Uy = v3arr_dot(slice->u_axis, local_y);
    curve->data.torus.Uz = v3arr_dot(slice->u_axis, local_z);
    curve->data.torus.Vx = v3arr_dot(slice->v_axis, local_x);
    curve->data.torus.Vy = v3arr_dot(slice->v_axis, local_y);
    curve->data.torus.Vz = v3arr_dot(slice->v_axis, local_z);

    return true;
}

/* ============================================================================
 * MAIN DISPATCH FUNCTION
 * ============================================================================ */

bool alea_intersect_primitive_plane(alea_primitive_type_t type,
                                   const alea_primitive_data_t* data,
                                   const alea_slice_plane_t* plane,
                                   alea_curve_2d_t* curve) {
    memset(curve, 0, sizeof(*curve));
    curve->type = ALEA_CURVE_NONE;

    switch (type) {
        case ALEA_PRIMITIVE_PLANE:
            return intersect_plane(&data->plane, plane, curve);

        case ALEA_PRIMITIVE_SPHERE:
        case ALEA_PRIMITIVE_SPH:
            return intersect_sphere(&data->sphere, plane, curve);

        case ALEA_PRIMITIVE_CYLINDER_X:
            return intersect_cylinder(data->cyl_x.center_y, data->cyl_x.center_z,
                                      data->cyl_x.radius, 0, plane, curve);

        case ALEA_PRIMITIVE_CYLINDER_Y:
            return intersect_cylinder(data->cyl_y.center_x, data->cyl_y.center_z,
                                      data->cyl_y.radius, 1, plane, curve);

        case ALEA_PRIMITIVE_CYLINDER_Z:
            return intersect_cylinder(data->cyl_z.center_x, data->cyl_z.center_y,
                                      data->cyl_z.radius, 2, plane, curve);

        case ALEA_PRIMITIVE_CONE_X:
            return intersect_cone(data->cone_x.apex_x, data->cone_x.apex_y,
                                  data->cone_x.apex_z, data->cone_x.tan_angle_sq,
                                  0, plane, curve);

        case ALEA_PRIMITIVE_CONE_Y:
            return intersect_cone(data->cone_y.apex_x, data->cone_y.apex_y,
                                  data->cone_y.apex_z, data->cone_y.tan_angle_sq,
                                  1, plane, curve);

        case ALEA_PRIMITIVE_CONE_Z:
            return intersect_cone(data->cone_z.apex_x, data->cone_z.apex_y,
                                  data->cone_z.apex_z, data->cone_z.tan_angle_sq,
                                  2, plane, curve);

        case ALEA_PRIMITIVE_RPP:
            return intersect_box(&data->box, plane, curve);

        case ALEA_PRIMITIVE_QUADRIC:
            return intersect_quadric(&data->quadric, plane, curve);

        case ALEA_PRIMITIVE_RCC:
            return intersect_rcc(&data->rcc, plane, curve);

        case ALEA_PRIMITIVE_TRC:
            return intersect_trc(&data->trc, plane, curve);

        case ALEA_PRIMITIVE_TORUS_X:
        case ALEA_PRIMITIVE_TORUS_Y:
        case ALEA_PRIMITIVE_TORUS_Z:
            return intersect_torus(&data->torus, type, plane, curve);

        default:
            return false;
    }
}

/* ============================================================================
 * SYSTEM-WIDE CURVE COMPUTATION
 * ============================================================================ */

/** Ensure room for one more curve. Returns 0 on success, -1 on OOM. */
static int ensure_curve_capacity(alea_curve_collection_t* r) {
    if (r->curves.count < r->curves.capacity) return 0;
    int res = alea_vec_reserve(&r->curves, r->curves.capacity * 2, alea_curve_2d_t);
    return res != 0 ? -1 : 0;
}

int alea_compute_slice_curves(alea_system_t* sys,
                             const alea_slice_plane_t* plane,
                             alea_curve_collection_t* result) {
    if (!sys || !plane || !result) return -1;

    memset(result, 0, sizeof(*result));

    /* Initial allocation */
    alea_vec_init(&result->curves);
    int vres = alea_vec_reserve(&result->curves, 64, alea_curve_2d_t);
    if (vres != 0) return -1;

    result->u_min = result->v_min = DBL_MAX;
    result->u_max = result->v_max = -DBL_MAX;

    /* Iterate over all surfaces */
    for (size_t i = 0; i < alea_vec_count(&sys->surfaces); i++) {
        const alea_surface_entry_t* surf = &sys->surfaces.data[i];
        const alea_primitive_entry_t* prim = &sys->primitives.data[surf->primitive_id];

        alea_curve_2d_t curve;
        if (alea_intersect_primitive_plane(prim->type, &prim->data, plane, &curve)) {
            /* Store surface ID */
            curve.surface_id = surf->mc_surface_id;

            if (ensure_curve_capacity(result) != 0) {
                alea_curve_collection_free(result);
                return -1;
            }
            result->curves.data[result->curves.count++] = curve;

            /* Update bounding box */
            double cu_min, cu_max, cv_min, cv_max;
            alea_curve_bbox(&curve, &cu_min, &cu_max, &cv_min, &cv_max);
            if (cu_min < result->u_min) result->u_min = cu_min;
            if (cu_max > result->u_max) result->u_max = cu_max;
            if (cv_min < result->v_min) result->v_min = cv_min;
            if (cv_max > result->v_max) result->v_max = cv_max;
        }
    }

    return 0;
}

void alea_curve_collection_free(alea_curve_collection_t* result) {
    if (result) {
        alea_vec_free(&result->curves);
        memset(result, 0, sizeof(*result));
    }
}

/* ============================================================================
 * CURVE EVALUATION AND RASTERIZATION
 * ============================================================================ */

bool alea_curve_eval(const alea_curve_2d_t* curve, double t,
                    double* u, double* v) {
    switch (curve->type) {
        case ALEA_CURVE_POINT:
            *u = curve->data.point[0];
            *v = curve->data.point[1];
            return true;

        case ALEA_CURVE_LINE:
        case ALEA_CURVE_LINE_SEGMENT:
        case ALEA_CURVE_RAY:
            *u = curve->data.line.point[0] + t * curve->data.line.direction[0];
            *v = curve->data.line.point[1] + t * curve->data.line.direction[1];
            return true;

        case ALEA_CURVE_CIRCLE:
        case ALEA_CURVE_ARC:
            *u = curve->data.circle.center[0] + curve->data.circle.radius * cos(t);
            *v = curve->data.circle.center[1] + curve->data.circle.radius * sin(t);
            return true;

        case ALEA_CURVE_ELLIPSE:
        case ALEA_CURVE_ELLIPSE_ARC: {
            /* Canonical form pre-computed by finalize_conic_as_ellipse */
            const alea_ellipse_2d_t* ell = &curve->data.ellipse;
            double cos_a = cos(ell->angle);
            double sin_a = sin(ell->angle);
            double cos_t = cos(t);
            double sin_t = sin(t);
            *u = ell->center[0] + ell->semi_a * cos_t * cos_a - ell->semi_b * sin_t * sin_a;
            *v = ell->center[1] + ell->semi_a * cos_t * sin_a + ell->semi_b * sin_t * cos_a;
            return true;
        }

        case ALEA_CURVE_QUARTIC: {
            /* For quartic (torus), we can't easily parameterize */
            /* Return false - use scanline intersection instead */
            return false;
        }

        case ALEA_CURVE_POLYGON: {
            int n = curve->data.polygon.vertex_count;
            if (n < 2) return false;
            int idx = (int)t;
            if (idx < 0) idx = 0;
            if (idx >= n) idx = n - 1;
            *u = curve->data.polygon.vertices[idx][0];
            *v = curve->data.polygon.vertices[idx][1];
            return true;
        }

        default:
            return false;
    }
}

int alea_curve_scanline_intersect(const alea_curve_2d_t* curve,
                                 double v_scanline,
                                 double* u_out,
                                 int max_intersections) {
    int count = 0;

    switch (curve->type) {
        case ALEA_CURVE_LINE:
        case ALEA_CURVE_LINE_SEGMENT:
        case ALEA_CURVE_RAY: {
            /* Line: a*u + b*v + c = 0, solve for u when v = v_scanline */
            double a = curve->data.line.a;
            double b = curve->data.line.b;
            double c = curve->data.line.c;

            if (fabs(a) < EPSILON) {
                /* Line is horizontal - no intersection (or infinite) */
                return 0;
            }

            double u = -(b * v_scanline + c) / a;
            u_out[0] = u;
            return 1;
        }

        case ALEA_CURVE_CIRCLE:
        case ALEA_CURVE_ARC: {
            /* Circle: (u - cx)² + (v - cy)² = r² */
            /* With v = v_scanline: (u - cx)² = r² - (v_scanline - cy)² */
            double cx = curve->data.circle.center[0];
            double cy = curve->data.circle.center[1];
            double r = curve->data.circle.radius;

            double dv = v_scanline - cy;
            double d2 = r*r - dv*dv;

            if (d2 < -EPSILON) return 0;
            if (d2 < EPSILON) {
                /* Tangent */
                u_out[0] = cx;
                return 1;
            }

            double du = sqrt(d2);
            u_out[0] = cx - du;
            if (count < max_intersections - 1) {
                u_out[1] = cx + du;
                return 2;
            }
            return 1;
        }

        case ALEA_CURVE_ELLIPSE:
        case ALEA_CURVE_ELLIPSE_ARC: {
            /* Canonical ellipse form: solve scanline intersection directly.
             * Ellipse: ((u-cx)·cosθ + dv·sinθ)²/a² + (-(u-cx)·sinθ + dv·cosθ)²/b² = 1
             * where dv = v_scanline - cy. This is a quadratic in du = u - cx. */
            const alea_ellipse_2d_t* ell = &curve->data.ellipse;
            double dv = v_scanline - ell->center[1];
            double ct = cos(ell->angle), st = sin(ell->angle);
            double inv_a2 = 1.0 / (ell->semi_a * ell->semi_a);
            double inv_b2 = 1.0 / (ell->semi_b * ell->semi_b);

            double qa = ct*ct*inv_a2 + st*st*inv_b2;
            double qb = 2*dv*ct*st*(inv_a2 - inv_b2);
            double qc = dv*dv*(st*st*inv_a2 + ct*ct*inv_b2) - 1.0;

            double x1, x2;
            int n = solve_quadratic(qa, qb, qc, &x1, &x2);

            if (n >= 1) u_out[0] = ell->center[0] + x1;
            if (n >= 2 && max_intersections > 1) u_out[1] = ell->center[0] + x2;
            return (n > max_intersections) ? max_intersections : n;
        }

        case ALEA_CURVE_PARABOLA:
        case ALEA_CURVE_HYPERBOLA: {
            /* General conic: A*u² + B*u*v + C*v² + D*u + E*v + F = 0 */
            /* With v = v_scanline: A*u² + (B*v_scanline + D)*u + (C*v_scanline² + E*v_scanline + F) = 0 */
            double A = curve->data.conic.A;
            double B = curve->data.conic.B;
            double C = curve->data.conic.C;
            double D = curve->data.conic.D;
            double E = curve->data.conic.E;
            double F = curve->data.conic.F;

            double a = A;
            double b = B * v_scanline + D;
            double c = C * v_scanline * v_scanline + E * v_scanline + F;

            double x1, x2;
            int n = solve_quadratic(a, b, c, &x1, &x2);

            if (n >= 1) u_out[0] = x1;
            if (n >= 2 && max_intersections > 1) u_out[1] = x2;
            return (n > max_intersections) ? max_intersections : n;
        }

        case ALEA_CURVE_POLYGON: {
            /* Check each edge */
            int n = curve->data.polygon.vertex_count;
            for (int i = 0; i < n && count < max_intersections; i++) {
                int j = (i + 1) % n;
                double v0 = curve->data.polygon.vertices[i][1];
                double v1 = curve->data.polygon.vertices[j][1];
                double u0 = curve->data.polygon.vertices[i][0];
                double u1 = curve->data.polygon.vertices[j][0];

                /* Check if scanline crosses this edge */
                if ((v0 <= v_scanline && v1 > v_scanline) ||
                    (v1 <= v_scanline && v0 > v_scanline)) {
                    /* Linear interpolation */
                    double t = (v_scanline - v0) / (v1 - v0);
                    u_out[count++] = u0 + t * (u1 - u0);
                }
            }
            return count;
        }

        case ALEA_CURVE_QUARTIC: {
            /* Quartic (torus) scanline intersection */
            const alea_torus_2d_t* tor = &curve->data.torus;

            if (tor->mode == 2) {
                /* Two separate circles */
                /* Circle 1 */
                double dv1 = v_scanline - tor->c1v;
                double d2_1 = tor->r1*tor->r1 - dv1*dv1;
                if (d2_1 > EPSILON && count < max_intersections) {
                    double du1 = sqrt(d2_1);
                    u_out[count++] = tor->c1u - du1;
                    if (count < max_intersections) u_out[count++] = tor->c1u + du1;
                }

                /* Circle 2 */
                double dv2 = v_scanline - tor->c2v;
                double d2_2 = tor->r2*tor->r2 - dv2*dv2;
                if (d2_2 > EPSILON && count < max_intersections) {
                    double du2 = sqrt(d2_2);
                    u_out[count++] = tor->c2u - du2;
                    if (count < max_intersections) u_out[count++] = tor->c2u + du2;
                }
                return count;

            } else if (tor->mode == 1) {
                /* Two concentric circles */
                double dv = v_scanline - tor->c1v;

                /* Outer circle (r2) */
                double d2_outer = tor->r2*tor->r2 - dv*dv;
                if (d2_outer > EPSILON) {
                    double du_outer = sqrt(d2_outer);
                    u_out[count++] = tor->c1u - du_outer;
                    if (count < max_intersections) u_out[count++] = tor->c1u + du_outer;
                }

                /* Inner circle (r1) */
                if (tor->r1 > EPSILON) {
                    double d2_inner = tor->r1*tor->r1 - dv*dv;
                    if (d2_inner > EPSILON && count < max_intersections) {
                        double du_inner = sqrt(d2_inner);
                        u_out[count++] = tor->c1u - du_inner;
                        if (count < max_intersections) u_out[count++] = tor->c1u + du_inner;
                    }
                }
                return count;

            } else {
                /* General quartic torus: solve quartic equation */
                /*
                 * Torus in local coords: (sqrt(x² + y²) - R)² + z² = r²
                 * Squared form: (x² + y² + z² + R² - r²)² = 4R²(x² + y²)
                 *
                 * Substitute: x = Ox + u*Ux + v*Vx (and similar for y, z)
                 * With v = v_scanline, we get a quartic in u.
                 */
                double R = tor->R;
                double rr = tor->r;
                double R2 = R * R;
                double r2 = rr * rr;
                double K = R2 - r2;

                /* Point on slice: P = O + u*U + v_scanline*V */
                double v = v_scanline;

                /* Precompute: Pv = O + v*V (the v-dependent part) */
                double Pvx = tor->Ox + v * tor->Vx;
                double Pvy = tor->Oy + v * tor->Vy;
                double Pvz = tor->Oz + v * tor->Vz;

                /* Full point: P = Pv + u*U */
                /* x = Pvx + u*Ux, y = Pvy + u*Uy, z = Pvz + u*Uz */

                double Ux = tor->Ux, Uy = tor->Uy, Uz = tor->Uz;

                /* Let's define:
                 * ρ² = x² + y² = (Pvx + u*Ux)² + (Pvy + u*Uy)²
                 * σ² = x² + y² + z² = ρ² + (Pvz + u*Uz)²
                 *
                 * Torus equation: (σ² + K)² = 4R²ρ²
                 */

                /* ρ² = (Pvx + u*Ux)² + (Pvy + u*Uy)² */
                /* = Pvx² + Pvy² + 2u(Pvx*Ux + Pvy*Uy) + u²(Ux² + Uy²) */
                /* = rho0 + rho1*u + rho2*u² */
                double rho0 = Pvx*Pvx + Pvy*Pvy;
                double rho1 = 2*(Pvx*Ux + Pvy*Uy);
                double rho2 = Ux*Ux + Uy*Uy;

                /* σ² = ρ² + (Pvz + u*Uz)² */
                /* = rho0 + rho1*u + rho2*u² + Pvz² + 2u*Pvz*Uz + u²*Uz² */
                /* = sig0 + sig1*u + sig2*u² */
                double sig0 = rho0 + Pvz*Pvz;
                double sig1 = rho1 + 2*Pvz*Uz;
                double sig2 = rho2 + Uz*Uz;

                /* (σ² + K)² = 4R²ρ² */
                /* (sig0 + sig1*u + sig2*u² + K)² = 4R²(rho0 + rho1*u + rho2*u²) */

                /* Let S = σ² + K = (sig0 + K) + sig1*u + sig2*u² = s0 + s1*u + s2*u² */
                double s0 = sig0 + K;
                double s1 = sig1;
                double s2 = sig2;

                /* S² = s0² + s1²u² + s2²u⁴ + 2s0*s1*u + 2s0*s2*u² + 2s1*s2*u³ */
                /* = s0² + 2s0*s1*u + (s1² + 2s0*s2)*u² + 2s1*s2*u³ + s2²*u⁴ */

                /* RHS = 4R²ρ² = 4R²*rho0 + 4R²*rho1*u + 4R²*rho2*u² */

                /* Quartic: a4*u⁴ + a3*u³ + a2*u² + a1*u + a0 = 0 */
                double R2_4 = 4 * R2;
                double a4 = s2 * s2;
                double a3 = 2 * s1 * s2;
                double a2 = s1*s1 + 2*s0*s2 - R2_4*rho2;
                double a1 = 2*s0*s1 - R2_4*rho1;
                double a0 = s0*s0 - R2_4*rho0;

                /* Solve quartic */
                double roots[4];
                int n_roots = alea_solve_quartic(a4, a3, a2, a1, a0, roots);

                /* Copy roots to output (up to max_intersections) */
                for (int i = 0; i < n_roots && count < max_intersections; i++) {
                    u_out[count++] = roots[i];
                }
                return count;
            }
        }

        case ALEA_CURVE_PARALLEL_LINES: {
            /* Two parallel lines: each line is point + t*direction */
            const alea_parallel_lines_2d_t* pl = &curve->data.parallel_lines;
            double dv = pl->direction[1];

            /* For each line, find intersection with scanline v = v_scanline */
            /* Line 1: (p1x + t*dx, p1y + t*dy) where p1y + t*dy = v_scanline */
            if (fabs(dv) > EPSILON) {
                double t1 = (v_scanline - pl->point1[1]) / dv;
                double u1 = pl->point1[0] + t1 * pl->direction[0];
                if (count < max_intersections) u_out[count++] = u1;

                double t2 = (v_scanline - pl->point2[1]) / dv;
                double u2 = pl->point2[0] + t2 * pl->direction[0];
                if (count < max_intersections) u_out[count++] = u2;
            }
            /* If dv ≈ 0, lines are horizontal - infinite intersections if on scanline, none otherwise */
            return count;
        }

        default:
            return 0;
    }
}

void alea_curve_bbox(const alea_curve_2d_t* curve,
                    double* u_min, double* u_max,
                    double* v_min, double* v_max) {
    *u_min = *v_min = DBL_MAX;
    *u_max = *v_max = -DBL_MAX;

    switch (curve->type) {
        case ALEA_CURVE_POINT:
            *u_min = *u_max = curve->data.point[0];
            *v_min = *v_max = curve->data.point[1];
            break;

        case ALEA_CURVE_LINE:
            /* Infinite line - use large bounds */
            *u_min = *v_min = -1e10;
            *u_max = *v_max = 1e10;
            break;

        case ALEA_CURVE_PARALLEL_LINES: {
            /* Two parallel infinite lines - compute bbox based on points and direction */
            const alea_parallel_lines_2d_t* pl = &curve->data.parallel_lines;
            double du = pl->direction[0];
            double dv = pl->direction[1];
            /* Span between the two points */
            *u_min = fmin(pl->point1[0], pl->point2[0]);
            *u_max = fmax(pl->point1[0], pl->point2[0]);
            *v_min = fmin(pl->point1[1], pl->point2[1]);
            *v_max = fmax(pl->point1[1], pl->point2[1]);
            /* Extend infinitely along direction */
            if (fabs(du) > EPSILON) {
                *u_min = -1e10;
                *u_max = 1e10;
            }
            if (fabs(dv) > EPSILON) {
                *v_min = -1e10;
                *v_max = 1e10;
            }
            break;
        }

        case ALEA_CURVE_CIRCLE:
        case ALEA_CURVE_ARC:
            *u_min = curve->data.circle.center[0] - curve->data.circle.radius;
            *u_max = curve->data.circle.center[0] + curve->data.circle.radius;
            *v_min = curve->data.circle.center[1] - curve->data.circle.radius;
            *v_max = curve->data.circle.center[1] + curve->data.circle.radius;
            break;

        case ALEA_CURVE_ELLIPSE:
        case ALEA_CURVE_ELLIPSE_ARC: {
            /* Canonical form already stored in data.ellipse */
            const alea_ellipse_2d_t* ell = &curve->data.ellipse;
            double cos_t = cos(ell->angle);
            double sin_t = sin(ell->angle);
            double du = fabs(ell->semi_a * cos_t) + fabs(ell->semi_b * sin_t);
            double dv = fabs(ell->semi_a * sin_t) + fabs(ell->semi_b * cos_t);
            *u_min = ell->center[0] - du;
            *u_max = ell->center[0] + du;
            *v_min = ell->center[1] - dv;
            *v_max = ell->center[1] + dv;
            break;
        }

        case ALEA_CURVE_PARABOLA:
        case ALEA_CURVE_HYPERBOLA:
            conic_bbox(&curve->data.conic, u_min, u_max, v_min, v_max);
            break;

        case ALEA_CURVE_QUARTIC: {
            /* Quartic (torus) bounding box */
            const alea_torus_2d_t* tor = &curve->data.torus;

            if (tor->mode == 2) {
                /* Two separate circles */
                *u_min = fmin(tor->c1u - tor->r1, tor->c2u - tor->r2);
                *u_max = fmax(tor->c1u + tor->r1, tor->c2u + tor->r2);
                *v_min = fmin(tor->c1v - tor->r1, tor->c2v - tor->r2);
                *v_max = fmax(tor->c1v + tor->r1, tor->c2v + tor->r2);
            } else if (tor->mode == 1) {
                /* Two concentric circles - use outer radius */
                *u_min = tor->c1u - tor->r2;
                *u_max = tor->c1u + tor->r2;
                *v_min = tor->c1v - tor->r2;
                *v_max = tor->c1v + tor->r2;
            } else {
                /* General quartic: compute conservative bound */
                /* The torus fits in a box of size 2(R+r) around origin in local coords */
                /* Transform the bounding box corners to slice coords */
                double R = tor->R;
                double rr = tor->r;
                double bound = R + rr;

                /* Sample the extreme points of the torus bounding box */
                /* In local coords, torus is bounded by |x|,|y| <= R+r and |z| <= r */
                /* We need to find the extent in (u,v) space */
                /* Conservative: sample corners of the local bounding box */
                double corners[8][3] = {
                    {-bound, -bound, -rr}, {-bound, -bound, rr},
                    {-bound, bound, -rr}, {-bound, bound, rr},
                    {bound, -bound, -rr}, {bound, -bound, rr},
                    {bound, bound, -rr}, {bound, bound, rr}
                };

                *u_min = *v_min = 1e10;
                *u_max = *v_max = -1e10;

                for (int i = 0; i < 8; i++) {
                    /* Local point = Ox + corners[i], so in slice coords: */
                    /* (u,v) projection of local point onto slice plane */
                    /* u = corners[i] · U_local, v = corners[i] · V_local */
                    /* But we stored the inverse: U_local = (Ux, Uy, Uz) */
                    /* Point in slice space: u = Ox*?... this is complex */

                    /* Simpler: the point in slice (u,v) that maps to local corner: */
                    /* We have: local = O_local + u*U_local + v*V_local */
                    /* So: corner - O_local = u*U_local + v*V_local */
                    /* Solve for u,v (least squares if not on plane) */

                    double dx = corners[i][0] - tor->Ox;
                    double dy = corners[i][1] - tor->Oy;
                    double dz = corners[i][2] - tor->Oz;

                    /* Project onto U and V (they're unit vectors in the plane) */
                    double u_coord = dx * tor->Ux + dy * tor->Uy + dz * tor->Uz;
                    double v_coord = dx * tor->Vx + dy * tor->Vy + dz * tor->Vz;

                    if (u_coord < *u_min) *u_min = u_coord;
                    if (u_coord > *u_max) *u_max = u_coord;
                    if (v_coord < *v_min) *v_min = v_coord;
                    if (v_coord > *v_max) *v_max = v_coord;
                }

                /* Add some margin */
                double margin = (R + rr) * 0.1;
                *u_min -= margin;
                *u_max += margin;
                *v_min -= margin;
                *v_max += margin;
            }
            break;
        }

        case ALEA_CURVE_POLYGON:
            for (int i = 0; i < curve->data.polygon.vertex_count; i++) {
                double u = curve->data.polygon.vertices[i][0];
                double v = curve->data.polygon.vertices[i][1];
                if (u < *u_min) *u_min = u;
                if (u > *u_max) *u_max = u;
                if (v < *v_min) *v_min = v;
                if (v > *v_max) *v_max = v;
            }
            break;

        default:
            break;
    }
}

/* ============================================================================
 * SLICE CURVES WITH SPATIAL INDEX
 * ============================================================================ */

#define SLICE_SPATIAL_TOL 1e-10
#define MAX_SPATIAL_HITS 262144  /* 256K - increased for large models */

/* Debug flag for slice curve generation - set to 1 to enable verbose output */
static int g_slice_curve_debug = 0;

void alea_slice_curve_set_debug(int enable) {
    g_slice_curve_debug = enable;
}

/* Forward declaration for point trace debug function */
extern void alea_set_debug_point_trace(int enable);

void alea_slice_point_trace_set_debug(int enable) {
    alea_set_debug_point_trace(enable);
}

/* Compare two 3x4 matrices element-wise within tolerance */
static bool matrices_equal(const double a[12], const double b[12], double tol) {
    for (int i = 0; i < 12; i++) {
        if (fabs(a[i] - b[i]) > tol) return false;
    }
    return true;
}

/* ---- Open-addressing hash set for seen surfaces ---- */

#define SEEN_HASH_EMPTY  0xFFFFFFFF  /* sentinel: surf_idx can never be this */
#define SEEN_HASH_LOAD_NUM   3       /* resize when count * 4 > capacity * 3 */
#define SEEN_HASH_LOAD_DEN   4

typedef struct {
    uint32_t surf_idx;               /* SEEN_HASH_EMPTY = slot unused */
    double m[12];                    /* Full 3x4 transform */
    int first_cell_id;
    int first_depth;
} seen_entry_t;

typedef struct {
    seen_entry_t* slots;
    size_t capacity;                 /* Always a power of two */
    size_t count;
} seen_set_t;

static void seen_set_init(seen_set_t* s, size_t initial_cap) {
    /* Round up to power of two */
    size_t cap = 64;
    while (cap < initial_cap) cap *= 2;
    s->slots = malloc(cap * sizeof(seen_entry_t));
    s->capacity = cap;
    s->count = 0;
    for (size_t i = 0; i < cap; i++)
        s->slots[i].surf_idx = SEEN_HASH_EMPTY;
}

static void seen_set_free(seen_set_t* s) {
    free(s->slots);
    s->slots = NULL;
    s->capacity = s->count = 0;
}

/* FNV-1a hash of (surf_idx, quantized matrix) */
static uint64_t seen_hash(uint32_t surf_idx, const double m[12]) {
    uint64_t h = 14695981039346656037ULL;
    /* Hash surf_idx */
    const uint8_t* p = (const uint8_t*)&surf_idx;
    for (size_t i = 0; i < sizeof(surf_idx); i++) {
        h ^= p[i]; h *= 1099511628211ULL;
    }
    /* Quantize doubles to 1e-8 grid and hash the int64 representation.
     * This ensures values within 1e-9 of each other may hash to the same
     * or adjacent buckets (we still do exact comparison in the probe). */
    for (int i = 0; i < 12; i++) {
        int64_t q = (int64_t)(m[i] * 1e8);
        p = (const uint8_t*)&q;
        for (size_t j = 0; j < sizeof(q); j++) {
            h ^= p[j]; h *= 1099511628211ULL;
        }
    }
    return h;
}

static void seen_set_grow(seen_set_t* s);

/* Returns true if already present. If not, inserts and returns false. */
static bool seen_set_insert(seen_set_t* s,
                             uint32_t surf_idx, const alea_matrix_t* transform,
                             int cell_id, int depth,
                             int* out_first_cell_id, int* out_first_depth) {
    /* Grow if needed */
    if (s->count * SEEN_HASH_LOAD_DEN >= s->capacity * SEEN_HASH_LOAD_NUM)
        seen_set_grow(s);

    uint64_t h = seen_hash(surf_idx, transform->m);
    size_t mask = s->capacity - 1;
    size_t idx = (size_t)(h & mask);

    for (;;) {
        seen_entry_t* e = &s->slots[idx];
        if (e->surf_idx == SEEN_HASH_EMPTY) {
            /* Empty slot — insert here */
            e->surf_idx = surf_idx;
            memcpy(e->m, transform->m, sizeof(e->m));
            e->first_cell_id = cell_id;
            e->first_depth = depth;
            s->count++;
            return false;
        }
        if (e->surf_idx == surf_idx &&
            matrices_equal(e->m, transform->m, SLICE_TRANSFORM_TOL)) {
            /* Already seen */
            if (out_first_cell_id) *out_first_cell_id = e->first_cell_id;
            if (out_first_depth) *out_first_depth = e->first_depth;
            return true;
        }
        idx = (idx + 1) & mask;
    }
}

static void seen_set_grow(seen_set_t* s) {
    size_t old_cap = s->capacity;
    seen_entry_t* old_slots = s->slots;

    size_t new_cap = old_cap * 2;
    s->slots = malloc(new_cap * sizeof(seen_entry_t));
    s->capacity = new_cap;
    s->count = 0;
    for (size_t i = 0; i < new_cap; i++)
        s->slots[i].surf_idx = SEEN_HASH_EMPTY;

    for (size_t i = 0; i < old_cap; i++) {
        if (old_slots[i].surf_idx == SEEN_HASH_EMPTY) continue;
        alea_matrix_t tmp;
        memcpy(tmp.m, old_slots[i].m, sizeof(tmp.m));
        seen_set_insert(s, old_slots[i].surf_idx, &tmp,
                        old_slots[i].first_cell_id, old_slots[i].first_depth,
                        NULL, NULL);
    }
    free(old_slots);
}

/**
 * Intersect 3D plane (a*x + b*y + c*z + d = 0) with slice plane,
 * viewport-cull, and append a ALEA_CURVE_LINE with surface_id=0.
 * Returns 0 on success (or cull/skip), -1 on OOM.
 */
static int emit_grid_line(const alea_slice_plane_t* plane,
                           alea_curve_collection_t* result,
                           double u_min, double u_max,
                           double v_min, double v_max,
                           double pa, double pb, double pc, double pd) {
    double A = pa * plane->u_axis[0] + pb * plane->u_axis[1]
             + pc * plane->u_axis[2];
    double B = pa * plane->v_axis[0] + pb * plane->v_axis[1]
             + pc * plane->v_axis[2];
    double C = pa * plane->origin[0] + pb * plane->origin[1]
             + pc * plane->origin[2] + pd;
    double len = sqrt(A * A + B * B);
    if (len < EPSILON) return 0;  /* parallel to slice */

    /* Viewport cull: evaluate A*u+B*v+C at four corners */
    double v0 = A * u_min + B * v_min + C;
    double v1 = A * u_max + B * v_min + C;
    double v2 = A * u_min + B * v_max + C;
    double v3 = A * u_max + B * v_max + C;
    if ((v0 > 0 && v1 > 0 && v2 > 0 && v3 > 0) ||
        (v0 < 0 && v1 < 0 && v2 < 0 && v3 < 0))
        return 0;  /* outside viewport */

    alea_curve_2d_t gc;
    memset(&gc, 0, sizeof(gc));
    gc.type = ALEA_CURVE_LINE;
    gc.surface_id = 0;
    gc.data.line.a = A;
    gc.data.line.b = B;
    gc.data.line.c = C;
    gc.data.line.direction[0] = -B / len;
    gc.data.line.direction[1] =  A / len;
    if (fabs(A) > fabs(B)) {
        gc.data.line.point[0] = -C / A;
        gc.data.line.point[1] = 0;
    } else {
        gc.data.line.point[0] = 0;
        gc.data.line.point[1] = -C / B;
    }
    if (ensure_curve_capacity(result) != 0) return -1;
    result->curves.data[result->curves.count++] = gc;
    return 0;
}

int alea_compute_slice_curves_spatial(alea_system_t* sys,
                                     const alea_slice_plane_t* plane,
                                     double u_min, double u_max,
                                     double v_min, double v_max,
                                     alea_curve_collection_t* result) {
    if (!sys || !plane || !result) return -1;

    bool use_hier = alea_spatial_mode_is_hierarchical(sys);

    /* Use the prepared spatial index if available, fall back to brute-force. */
    if (use_hier) {
        if (alea_hier_spatial_index_needs_rebuild(sys)) {
            return alea_compute_slice_curves(sys, plane, result);
        }
    } else {
        if (alea_spatial_index_needs_rebuild(sys)) {
            return alea_compute_slice_curves(sys, plane, result);
        }
    }

    memset(result, 0, sizeof(*result));

    /* Initial allocation */
    alea_vec_init(&result->curves);
    int vres = alea_vec_reserve(&result->curves, 64, alea_curve_2d_t);
    if (vres != 0) return -1;

    result->u_min = result->v_min = DBL_MAX;
    result->u_max = result->v_max = -DBL_MAX;

    /* Build query bbox from plane and viewport */
    alea_bbox_t query_bbox;
    double nx = plane->normal[0];
    double ny = plane->normal[1];
    double nz = plane->normal[2];

    if (fabs(nz) > 0.99) {
        /* Z-plane */
        double z = plane->origin[2];
        query_bbox.min_x = u_min; query_bbox.max_x = u_max;
        query_bbox.min_y = v_min; query_bbox.max_y = v_max;
        query_bbox.min_z = z - SLICE_SPATIAL_TOL;
        query_bbox.max_z = z + SLICE_SPATIAL_TOL;
    } else if (fabs(ny) > 0.99) {
        /* Y-plane */
        double y = plane->origin[1];
        query_bbox.min_x = u_min; query_bbox.max_x = u_max;
        query_bbox.min_y = y - SLICE_SPATIAL_TOL;
        query_bbox.max_y = y + SLICE_SPATIAL_TOL;
        query_bbox.min_z = v_min; query_bbox.max_z = v_max;
    } else if (fabs(nx) > 0.99) {
        /* X-plane */
        double x = plane->origin[0];
        query_bbox.min_x = x - SLICE_SPATIAL_TOL;
        query_bbox.max_x = x + SLICE_SPATIAL_TOL;
        query_bbox.min_y = u_min; query_bbox.max_y = u_max;
        query_bbox.min_z = v_min; query_bbox.max_z = v_max;
    } else {
        /* Arbitrary plane: compute 3D bbox of the viewport quad.
         * Four corners: origin + u*U + v*V at (u_min,v_min), (u_max,v_min),
         * (u_min,v_max), (u_max,v_max). */
        double uvals[2] = {u_min, u_max};
        double vvals[2] = {v_min, v_max};
        query_bbox.min_x = query_bbox.min_y = query_bbox.min_z = DBL_MAX;
        query_bbox.max_x = query_bbox.max_y = query_bbox.max_z = -DBL_MAX;
        for (int ui = 0; ui < 2; ui++) {
            for (int vi = 0; vi < 2; vi++) {
                double px = plane->origin[0] + uvals[ui] * plane->u_axis[0]
                                             + vvals[vi] * plane->v_axis[0];
                double py = plane->origin[1] + uvals[ui] * plane->u_axis[1]
                                             + vvals[vi] * plane->v_axis[1];
                double pz = plane->origin[2] + uvals[ui] * plane->u_axis[2]
                                             + vvals[vi] * plane->v_axis[2];
                if (px < query_bbox.min_x) query_bbox.min_x = px;
                if (px > query_bbox.max_x) query_bbox.max_x = px;
                if (py < query_bbox.min_y) query_bbox.min_y = py;
                if (py > query_bbox.max_y) query_bbox.max_y = py;
                if (pz < query_bbox.min_z) query_bbox.min_z = pz;
                if (pz > query_bbox.max_z) query_bbox.max_z = pz;
            }
        }
        /* Small margin for surfaces tangent to the bbox edges */
        query_bbox.min_x -= SLICE_SPATIAL_TOL;
        query_bbox.max_x += SLICE_SPATIAL_TOL;
        query_bbox.min_y -= SLICE_SPATIAL_TOL;
        query_bbox.max_y += SLICE_SPATIAL_TOL;
        query_bbox.min_z -= SLICE_SPATIAL_TOL;
        query_bbox.max_z += SLICE_SPATIAL_TOL;
    }

    /* Query instances that intersect the slice */
    size_t max_hits = use_hier ? MAX_SPATIAL_HITS : alea_spatial_get_instance_count(sys);
    if (max_hits == 0) max_hits = MAX_SPATIAL_HITS;
    if (max_hits > MAX_SPATIAL_HITS) max_hits = MAX_SPATIAL_HITS;

    alea_spatial_hit_t* hits = malloc(max_hits * sizeof(alea_spatial_hit_t));
    if (!hits) {
        alea_vec_free(&result->curves);
        return -1;
    }

    int hit_count = use_hier
        ? alea_hier_spatial_query_region(sys, &query_bbox, hits, max_hits)
        : alea_spatial_query_region(sys, &query_bbox, hits, max_hits);
    if (hit_count < 0) {
        free(hits);
        alea_vec_free(&result->curves);
        return -1;
    }

    /* Deduplicate hits - remove duplicate (cell_index, transform) pairs */
    /* This handles the case where the spatial index contains multiple instances
     * of the same cell with the same transform (e.g., from multiple parent fills) */
    int unique_count = 0;
    for (int h = 0; h < hit_count; h++) {
        bool is_duplicate = false;
        for (int k = 0; k < unique_count; k++) {
            if (hits[k].cell_index == hits[h].cell_index &&
                matrices_equal(hits[k].transform.m, hits[h].transform.m, SLICE_TRANSFORM_TOL)) {
                is_duplicate = true;
                break;
            }
        }
        if (!is_duplicate) {
            if (unique_count != h) {
                hits[unique_count] = hits[h];
            }
            unique_count++;
        }
    }

    if (g_slice_curve_debug && unique_count < hit_count) {
        ALEA_LOG_DEBUG("Deduplicated %d -> %d cell instances", hit_count, unique_count);
    }
    hit_count = unique_count;

    /* Track seen surfaces to avoid duplicates (hash set, O(1) amortized) */
    seen_set_t seen;
    seen_set_init(&seen, 256);

    if (g_slice_curve_debug) {
        ALEA_LOG_DEBUG("Processing %d cell instances for slice curves", hit_count);
    }

    /* For each hit instance, compute curves */
    for (int h = 0; h < hit_count; h++) {
        alea_spatial_hit_t* hit = &hits[h];
        const alea_cell_entry_t* cell = &sys->cells.data[hit->cell_index];

        if (g_slice_curve_debug) {
            ALEA_LOG_DEBUG("Hit %d: cell_id=%d cell_idx=%d universe=%d depth=%d terminal=%d transform=(%.3f,%.3f,%.3f)",
                   h, hit->cell_id, hit->cell_index, hit->universe_id, hit->depth, hit->is_terminal,
                   hit->transform.m[3], hit->transform.m[7], hit->transform.m[11]);
        }

        /* Get surfaces for this cell */
        if (!cell->surface_indices || cell->surface_index_count == 0) {
            if (g_slice_curve_debug) {
                ALEA_LOG_DEBUG("  -> No surfaces for this cell");
            }
            continue;
        }

        for (size_t s = 0; s < cell->surface_index_count; s++) {
            uint32_t surf_idx = cell->surface_indices[s];
            if (surf_idx >= alea_vec_count(&sys->surfaces)) continue;

            const alea_surface_entry_t* surf = &sys->surfaces.data[surf_idx];

            /* Skip if we've already processed this surface with this transform */
            int first_cell_id = -1, first_depth = -1;
            if (seen_set_insert(&seen, surf_idx, &hit->transform,
                                hit->cell_id, hit->depth,
                                &first_cell_id, &first_depth)) {
                if (g_slice_curve_debug) {
                    ALEA_LOG_DEBUG("  -> SKIPPED surface %d (mcnp=%d): already seen from cell_id=%d depth=%d",
                           surf_idx, surf->mc_surface_id, first_cell_id, first_depth);
                }
                continue;
            }

            const alea_primitive_entry_t* prim = &sys->primitives.data[surf->primitive_id];

            /* Transform primitive */
            alea_primitive_data_t transformed_data;
            alea_primitive_type_t transformed_type;

            if (!alea_primitive_transform(prim->type, &prim->data,
                                        hit->transform.m,
                                        &transformed_type, &transformed_data)) {
                if (g_slice_curve_debug) {
                    ALEA_LOG_DEBUG("  -> surface %d (mcnp=%d): transform failed",
                           surf_idx, surf->mc_surface_id);
                }
                continue;
            }

            /* Compute curve intersection */
            alea_curve_2d_t curve;
            if (!alea_intersect_primitive_plane(transformed_type, &transformed_data,
                                               plane, &curve)) {
                if (g_slice_curve_debug) {
                    ALEA_LOG_DEBUG("  -> surface %d (mcnp=%d): no intersection with plane",
                           surf_idx, surf->mc_surface_id);
                }
                continue;
            }

            /* Cull curves outside viewport */
            double cu_min, cu_max, cv_min, cv_max;
            alea_curve_bbox(&curve, &cu_min, &cu_max, &cv_min, &cv_max);

            if (cu_max < u_min || cu_min > u_max ||
                cv_max < v_min || cv_min > v_max) {
                if (g_slice_curve_debug) {
                    ALEA_LOG_DEBUG("  -> surface %d (mcnp=%d): culled (outside viewport)",
                           surf_idx, surf->mc_surface_id);
                }
                continue;
            }

            curve.surface_id = surf->mc_surface_id;
            curve.universe_id = cell->universe_id;

            /* Add curve to result */
            if (ensure_curve_capacity(result) != 0) {
                seen_set_free(&seen);
                free(hits);
                alea_curve_collection_free(result);
                return -1;
            }
            result->curves.data[result->curves.count++] = curve;

            if (g_slice_curve_debug) {
                const char* type_names[] = {
                    "NONE", "POINT", "LINE", "LINE_SEGMENT", "RAY", "CIRCLE", "ARC",
                    "ELLIPSE", "ELLIPSE_ARC", "PARABOLA", "HYPERBOLA", "POLYGON",
                    "QUARTIC", "PARALLEL_LINES"
                };
                const char* type_name = (curve.type < 14) ? type_names[curve.type] : "UNKNOWN";
                ALEA_LOG_DEBUG("  -> ADDED surface %d (mcnp=%d) as curve type=%s bbox=(%.2f,%.2f)-(%.2f,%.2f)",
                       surf_idx, surf->mc_surface_id, type_name, cu_min, cv_min, cu_max, cv_max);
                if (curve.type == ALEA_CURVE_QUARTIC) {
                    ALEA_LOG_DEBUG("     QUARTIC mode=%d: c1=(%.2f,%.2f) r1=%.2f, c2=(%.2f,%.2f) r2=%.2f",
                           curve.data.torus.mode,
                           curve.data.torus.c1u, curve.data.torus.c1v, curve.data.torus.r1,
                           curve.data.torus.c2u, curve.data.torus.c2v, curve.data.torus.r2);
                }
            }

            /* (surface+transform already inserted by seen_set_insert above) */

            /* Update bounds */
            if (cu_min < result->u_min) result->u_min = cu_min;
            if (cu_max > result->u_max) result->u_max = cu_max;
            if (cv_min < result->v_min) result->v_min = cv_min;
            if (cv_max > result->v_max) result->v_max = cv_max;
        }
    }

    /* ---- Lattice element curve stamping ----
     * The spatial index doesn't expand lattice fills.  For each lattice cell
     * visible in the viewport, compute which elements are in range and stamp
     * the fill-universe surfaces at each element's translation offset. */
    for (size_t ci = 0; ci < alea_vec_count(&sys->cells); ci++) {
        const alea_cell_entry_t* cell = &sys->cells.data[ci];
        if (cell->lat_type == 0 || !cell->lat_fill) continue;

        int imin = cell->lat_fill_dims[0], imax = cell->lat_fill_dims[1];
        int jmin = cell->lat_fill_dims[2], jmax = cell->lat_fill_dims[3];
        int kmin = cell->lat_fill_dims[4], kmax = cell->lat_fill_dims[5];
        int nj = jmax - jmin + 1;
        int nk = kmax - kmin + 1;
        double p = cell->lat_pitch[0];


        for (int ei = imin; ei <= imax; ei++) {
            for (int ej = jmin; ej <= jmax; ej++) {
                for (int ek = kmin; ek <= kmax; ek++) {
                    /* Element center */
                    double cx, cy, cz;
                    if (cell->lat_type == 2) {
                        /* Hex: axial (ei, ej) */
                        cx = ei * p + ej * p * 0.5;
                        cy = ej * p * M_SQRT3 * 0.5;
                    } else {
                        /* Rect */
                        cx = cell->lat_lower_left[0]
                           + (ei - imin + 0.5) * cell->lat_pitch[0];
                        cy = cell->lat_lower_left[1]
                           + (ej - jmin + 0.5) * cell->lat_pitch[1];
                    }
                    cz = (nk == 1) ? 0.0
                       : cell->lat_lower_left[2]
                       + (ek - kmin + 0.5) * cell->lat_pitch[2];

                    /* Quick viewport cull: element center ± pitch */
                    double margin = (cell->lat_type == 2)
                                  ? p / M_SQRT3 : cell->lat_pitch[0];
                    int in_viewport = 0;
                    if (fabs(nz) > 0.99) {
                        if (cx + margin >= u_min && cx - margin <= u_max &&
                            cy + margin >= v_min && cy - margin <= v_max)
                            in_viewport = 1;
                    } else if (fabs(ny) > 0.99) {
                        if (cx + margin >= u_min && cx - margin <= u_max &&
                            cz + margin >= v_min && cz - margin <= v_max)
                            in_viewport = 1;
                    } else if (fabs(nx) > 0.99) {
                        if (cy + margin >= u_min && cy - margin <= u_max &&
                            cz + margin >= v_min && cz - margin <= v_max)
                            in_viewport = 1;
                    } else {
                        in_viewport = 1;  /* arbitrary plane — can't cull */
                    }
                    if (!in_viewport) continue;

                    /* Fill universe */
                    int oi = ei - imin, oj = ej - jmin, ok = ek - kmin;
                    size_t idx = (size_t)(oi * nj * nk + oj * nk + ok);
                    if (idx >= cell->lat_fill_count) continue;
                    int univ_id = cell->lat_fill[idx];

                    const alea_universe_t* fill_univ =
                        alea_get_universe(sys, univ_id);
                    if (!fill_univ) continue;

                    /* Translation transform for this element */
                    alea_matrix_t elem_tr;
                    alea_matrix_identity(&elem_tr);
                    elem_tr.m[3]  = cx;
                    elem_tr.m[7]  = cy;
                    elem_tr.m[11] = cz;

                    /* Stamp each surface of every cell in the fill universe */
                    for (size_t fi = 0; fi < fill_univ->cell_indices.count; fi++) {
                        size_t fc = fill_univ->cell_indices.data[fi];
                        const alea_cell_entry_t* fc_cell =
                            &sys->cells.data[fc];
                        if (!fc_cell->surface_indices) continue;

                        for (size_t s = 0; s < fc_cell->surface_index_count;
                             s++) {
                            uint32_t surf_idx = fc_cell->surface_indices[s];
                            if (surf_idx >= alea_vec_count(&sys->surfaces))
                                continue;

                            if (seen_set_insert(&seen, surf_idx, &elem_tr,
                                               fc_cell->mc_cell_id, 0,
                                               NULL, NULL))
                                continue;

                            const alea_surface_entry_t* surf =
                                &sys->surfaces.data[surf_idx];
                            const alea_primitive_entry_t* prim =
                                &sys->primitives.data[surf->primitive_id];

                            alea_primitive_data_t td;
                            alea_primitive_type_t tt;
                            if (!alea_primitive_transform(prim->type,
                                    &prim->data, elem_tr.m, &tt, &td))
                                continue;

                            alea_curve_2d_t curve;
                            if (!alea_intersect_primitive_plane(
                                    tt, &td, plane, &curve))
                                continue;

                            /* Viewport cull */
                            double cu0, cu1, cv0, cv1;
                            alea_curve_bbox(&curve, &cu0, &cu1, &cv0, &cv1);
                            if (cu1 < u_min || cu0 > u_max ||
                                cv1 < v_min || cv0 > v_max)
                                continue;

                            curve.surface_id = surf->mc_surface_id;
                            curve.universe_id = fc_cell->universe_id;

                            if (ensure_curve_capacity(result) != 0) goto lat_done;
                            result->curves.data[result->curves.count++] = curve;

                            if (cu0 < result->u_min) result->u_min = cu0;
                            if (cu1 > result->u_max) result->u_max = cu1;
                            if (cv0 < result->v_min) result->v_min = cv0;
                            if (cv1 > result->v_max) result->v_max = cv1;
                        }
                    }
                }
            }
        }
    }

    /* ---- Lattice grid lines ----
     * Emit analytical line curves at every lattice element boundary so the
     * grid is visible in the rendered slice. */
    for (size_t ci = 0; ci < alea_vec_count(&sys->cells); ci++) {
        const alea_cell_entry_t* cell = &sys->cells.data[ci];
        if (cell->lat_type == 0 || !cell->lat_fill) continue;

        int imin = cell->lat_fill_dims[0], imax = cell->lat_fill_dims[1];
        int jmin = cell->lat_fill_dims[2], jmax = cell->lat_fill_dims[3];
        int kmin = cell->lat_fill_dims[4], kmax = cell->lat_fill_dims[5];
        int ni = imax - imin + 1;
        int nj_grid = jmax - jmin + 1;
        int nk_grid = kmax - kmin + 1;
        double px = cell->lat_pitch[0];
        double py = cell->lat_pitch[1];
        double pz = cell->lat_pitch[2];


        #define GRID_LINE(a, b, c, d)                                         \
            if (emit_grid_line(plane, result, u_min, u_max, v_min, v_max,     \
                               (a), (b), (c), (d)) < 0)                       \
                goto lat_done

        if (cell->lat_type == 1) {
            /* Rectangular lattice — axis-aligned boundary planes */
            for (int n = 0; n <= ni; n++) {
                double xval = cell->lat_lower_left[0] + n * px;
                GRID_LINE(1.0, 0.0, 0.0, -xval);
            }
            for (int n = 0; n <= nj_grid; n++) {
                double yval = cell->lat_lower_left[1] + n * py;
                GRID_LINE(0.0, 1.0, 0.0, -yval);
            }
            if (nk_grid > 1) {
                for (int n = 0; n <= nk_grid; n++) {
                    double zval = cell->lat_lower_left[2] + n * pz;
                    GRID_LINE(0.0, 0.0, 1.0, -zval);
                }
            }
        } else {
            /* Hexagonal lattice — three families of oblique planes.
             *
             * Hex element centers: cx = i*p + j*p/2, cy = j*p*sqrt3/2
             * Oblique coordinates: fi = x/p - y/(p*sqrt3), fj = 2y/(p*sqrt3)
             * Boundaries at half-integer values of fi, fj, fi+fj.
             *
             * Family 1: 1*x + (-1/sqrt3)*y - (n+0.5)*p = 0
             * Family 2: 0*x + 1*y - (n+0.5)*p*sqrt3/2 = 0
             * Family 3: 1*x + (1/sqrt3)*y - (n+0.5)*p = 0
             */
            double inv_sqrt3 = 1.0 / M_SQRT3;

            for (int n = imin - 1; n <= imax; n++) {
                double d = (n + 0.5) * px;
                GRID_LINE(1.0, -inv_sqrt3, 0.0, -d);
            }
            for (int n = jmin - 1; n <= jmax; n++) {
                double d = (n + 0.5) * px * M_SQRT3 * 0.5;
                GRID_LINE(0.0, 1.0, 0.0, -d);
            }
            for (int n = imin + jmin - 1; n <= imax + jmax; n++) {
                double d = (n + 0.5) * px;
                GRID_LINE(1.0, inv_sqrt3, 0.0, -d);
            }
        }

        #undef GRID_LINE
    }

lat_done:
    seen_set_free(&seen);
    free(hits);
    return 0;
}
