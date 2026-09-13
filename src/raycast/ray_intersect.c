// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/**
 * @file ray_intersect.c
 * @brief Ray-primitive intersection math
 *
 * Standard analytical ray-surface intersection formulas.
 * All functions assume ray direction is normalized.
 */

#include "ray_intersect.h"
#include "ray_epsilon.h"
#include "primitives/primitive_eval.h"
#include "primitives/primitive_desc.h"
#include "util/poly_solve.h"
#include <math.h>
#include <float.h>

/**
 * Numerically stable quadratic solver using Vieta's formula.
 * Avoids catastrophic cancellation in (-b ± sqrt(disc)) / (2a).
 * Returns number of real roots (0, 1, or 2) in t_out[].
 */
static inline int solve_quadratic_stable(double a, double b, double c,
                                         double* t_out) {
    if (fabs(a) < RAY_EPSILON) {
        if (fabs(b) < RAY_EPSILON) return 0;
        t_out[0] = -c / b;
        return 1;
    }

    double discriminant = b * b - 4.0 * a * c;
    if (discriminant < -DISCRIMINANT_TOL) return 0;
    if (discriminant < 0) discriminant = 0;

    double sqrt_disc = sqrt(discriminant);
    double q = -0.5 * (b + copysign(sqrt_disc, b));

    if (fabs(q) < DISCRIMINANT_TOL) {
        /* Double root fallback */
        t_out[0] = -b / (2.0 * a);
        t_out[1] = t_out[0];
        return 2;
    }

    t_out[0] = q / a;
    t_out[1] = c / q;

    /* Ensure t_out[0] <= t_out[1] */
    if (t_out[0] > t_out[1]) {
        double tmp = t_out[0];
        t_out[0] = t_out[1];
        t_out[1] = tmp;
    }
    return 2;
}

/* ============================================================================
 * PLANE
 * ============================================================================ */

int ray_intersect_plane(const alea_ray_t* ray,
                        const alea_plane_data_t* plane,
                        double* restrict t_out,
                        double* restrict nx, double* restrict ny, double* restrict nz) {
    /*
     * Plane: ax + by + cz + d = 0
     * Ray: P = O + t*D
     * Substitute: a(ox + t*dx) + b(oy + t*dy) + c(oz + t*dz) + d = 0
     * Solve: t = -(a*ox + b*oy + c*oz + d) / (a*dx + b*dy + c*dz)
     */
    double denom = plane->a * ray->dx + plane->b * ray->dy + plane->c * ray->dz;

    if (fabs(denom) < RAY_EPSILON) {
        return 0;  /* Ray parallel to plane */
    }

    double numer = plane->a * ray->ox + plane->b * ray->oy +
                   plane->c * ray->oz + plane->d;

    t_out[0] = -numer / denom;

    /* Normal: planes are pre-normalized by alea_canonicalize_primitive() */
    if (nx) {
        *nx = plane->a;
        *ny = plane->b;
        *nz = plane->c;
    }

    return 1;
}

/* ============================================================================
 * SPHERE
 * ============================================================================ */

int ray_intersect_sphere(const alea_ray_t* ray,
                         const alea_sphere_data_t* sphere,
                         double* restrict t_out) {
    /*
     * Sphere: |P - C|^2 = r^2
     * Ray: P = O + t*D
     * Let L = O - C
     * |L + t*D|^2 = r^2
     * t^2|D|^2 + 2t(L.D) + |L|^2 - r^2 = 0
     *
     * Since D is normalized, |D|^2 = 1
     * t^2 + 2t(L.D) + |L|^2 - r^2 = 0
     */
    double lx = ray->ox - sphere->center_x;
    double ly = ray->oy - sphere->center_y;
    double lz = ray->oz - sphere->center_z;

    /* a = |D|^2 = 1 (normalized), b_half = L.D, c = |L|^2 - r^2 */
    double b_half = lx * ray->dx + ly * ray->dy + lz * ray->dz;
    double c = lx * lx + ly * ly + lz * lz - sphere->radius * sphere->radius;

    double discriminant = b_half * b_half - c;

    if (discriminant < -DISCRIMINANT_TOL) return 0;
    if (discriminant < 0) discriminant = 0;

    /* Stable form: q = -(b_half + sign(b_half)*sqrt(disc)) */
    double sqrt_disc = sqrt(discriminant);
    double q = -(b_half + copysign(sqrt_disc, b_half));

    if (fabs(q) < DISCRIMINANT_TOL) {
        t_out[0] = -b_half;
        t_out[1] = -b_half;
        return 2;
    }

    t_out[0] = q;       /* q / a where a = 1 */
    t_out[1] = c / q;   /* Vieta: t0*t1 = c/a = c */

    if (t_out[0] > t_out[1]) {
        double tmp = t_out[0];
        t_out[0] = t_out[1];
        t_out[1] = tmp;
    }
    return 2;
}

/* ============================================================================
 * AXIS-ALIGNED CYLINDERS (INFINITE) — consolidated
 * ============================================================================ */

/**
 * Unified cylinder intersection: axis selects which 2D plane to project onto.
 * The center is a full 3D point so coordinate and center are extracted with
 * the same index — the perpendicular-axis order can never get out of sync
 * with the call site (a swapped (x,z) pair for Y-cylinders once misplaced
 * every C/Y intersection whose center had unequal components).
 */
static int ray_intersect_cylinder(const alea_ray_t* ray,
                                  int axis, const double center[3],
                                  double radius, double* restrict t_out) {
    double o[3] = { ray->ox, ray->oy, ray->oz };
    double d[3] = { ray->dx, ray->dy, ray->dz };
    int a0 = (axis + 1) % 3;
    int a1 = (axis + 2) % 3;

    double p0 = o[a0] - center[a0];
    double p1 = o[a1] - center[a1];

    double a = d[a0] * d[a0] + d[a1] * d[a1];
    double b = 2.0 * (p0 * d[a0] + p1 * d[a1]);
    double c = p0 * p0 + p1 * p1 - radius * radius;

    return solve_quadratic_stable(a, b, c, t_out);
}

int ray_intersect_cylinder_x(const alea_ray_t* ray,
                             const alea_cylinder_x_data_t* cyl,
                             double* restrict t_out) {
    double center[3] = { 0.0, cyl->center_y, cyl->center_z };
    return ray_intersect_cylinder(ray, 0, center, cyl->radius, t_out);
}

int ray_intersect_cylinder_y(const alea_ray_t* ray,
                             const alea_cylinder_y_data_t* cyl,
                             double* restrict t_out) {
    double center[3] = { cyl->center_x, 0.0, cyl->center_z };
    return ray_intersect_cylinder(ray, 1, center, cyl->radius, t_out);
}

int ray_intersect_cylinder_z(const alea_ray_t* ray,
                             const alea_cylinder_z_data_t* cyl,
                             double* restrict t_out) {
    double center[3] = { cyl->center_x, cyl->center_y, 0.0 };
    return ray_intersect_cylinder(ray, 2, center, cyl->radius, t_out);
}

/* ============================================================================
 * AXIS-ALIGNED CONES
 * ============================================================================ */

/**
 * Unified cone intersection: axis selects which component is the cone axis.
 * axis=0 (X): perp=(Y,Z), axis=1 (Y): perp=(X,Z), axis=2 (Z): perp=(X,Y).
 * All three cone types share identical field layout (apex_x/y/z, tan_angle_sq, sheet_selection).
 */
static int ray_intersect_cone(const alea_ray_t* ray, int axis,
                              double apex_x, double apex_y, double apex_z,
                              double k2, int sheet_selection,
                              double* restrict t_out) {
    double p[3] = { ray->ox - apex_x, ray->oy - apex_y, ray->oz - apex_z };
    double d[3] = { ray->dx, ray->dy, ray->dz };
    int a0 = (axis + 1) % 3;
    int a1 = (axis + 2) % 3;

    double a = d[a0] * d[a0] + d[a1] * d[a1] - k2 * d[axis] * d[axis];
    double b = 2.0 * (p[a0] * d[a0] + p[a1] * d[a1] - k2 * p[axis] * d[axis]);
    double c = p[a0] * p[a0] + p[a1] * p[a1] - k2 * p[axis] * p[axis];

    double t_raw[2];
    int raw_count = solve_quadratic_stable(a, b, c, t_raw);

    /* Filter by sheet_selection */
    int count = 0;
    for (int i = 0; i < raw_count; i++) {
        double axis_val = p[axis] + t_raw[i] * d[axis];
        if (sheet_selection > 0 && axis_val < -RAY_EPSILON) continue;
        if (sheet_selection < 0 && axis_val > RAY_EPSILON) continue;
        t_out[count++] = t_raw[i];
    }
    return count;
}

int ray_intersect_cone_x(const alea_ray_t* ray,
                         const alea_cone_x_data_t* cone,
                         double* restrict t_out) {
    return ray_intersect_cone(ray, 0, cone->apex_x, cone->apex_y, cone->apex_z,
                              cone->tan_angle_sq, cone->sheet_selection, t_out);
}

int ray_intersect_cone_y(const alea_ray_t* ray,
                         const alea_cone_y_data_t* cone,
                         double* restrict t_out) {
    return ray_intersect_cone(ray, 1, cone->apex_x, cone->apex_y, cone->apex_z,
                              cone->tan_angle_sq, cone->sheet_selection, t_out);
}

int ray_intersect_cone_z(const alea_ray_t* ray,
                         const alea_cone_z_data_t* cone,
                         double* restrict t_out) {
    return ray_intersect_cone(ray, 2, cone->apex_x, cone->apex_y, cone->apex_z,
                              cone->tan_angle_sq, cone->sheet_selection, t_out);
}

/* ============================================================================
 * AXIS-ALIGNED BOX
 * ============================================================================ */

int ray_intersect_box(const alea_ray_t* ray,
                      const alea_box_data_t* box,
                      double* restrict t_out) {
    /*
     * Slab method: intersect ray with each pair of parallel planes,
     * find overlapping interval.
     */
    double t_min = -DBL_MAX;
    double t_max = DBL_MAX;

    /* X slabs */
    if (fabs(ray->dx) < RAY_EPSILON) {
        if (ray->ox < box->min_x || ray->ox > box->max_x) {
            return 0;
        }
    } else {
        double inv_d = 1.0 / ray->dx;
        double t1 = (box->min_x - ray->ox) * inv_d;
        double t2 = (box->max_x - ray->ox) * inv_d;
        if (t1 > t2) { double tmp = t1; t1 = t2; t2 = tmp; }
        if (t1 > t_min) t_min = t1;
        if (t2 < t_max) t_max = t2;
        if (t_min > t_max) return 0;
    }

    /* Y slabs */
    if (fabs(ray->dy) < RAY_EPSILON) {
        if (ray->oy < box->min_y || ray->oy > box->max_y) {
            return 0;
        }
    } else {
        double inv_d = 1.0 / ray->dy;
        double t1 = (box->min_y - ray->oy) * inv_d;
        double t2 = (box->max_y - ray->oy) * inv_d;
        if (t1 > t2) { double tmp = t1; t1 = t2; t2 = tmp; }
        if (t1 > t_min) t_min = t1;
        if (t2 < t_max) t_max = t2;
        if (t_min > t_max) return 0;
    }

    /* Z slabs */
    if (fabs(ray->dz) < RAY_EPSILON) {
        if (ray->oz < box->min_z || ray->oz > box->max_z) {
            return 0;
        }
    } else {
        double inv_d = 1.0 / ray->dz;
        double t1 = (box->min_z - ray->oz) * inv_d;
        double t2 = (box->max_z - ray->oz) * inv_d;
        if (t1 > t2) { double tmp = t1; t1 = t2; t2 = tmp; }
        if (t1 > t_min) t_min = t1;
        if (t2 < t_max) t_max = t2;
        if (t_min > t_max) return 0;
    }

    t_out[0] = t_min;
    t_out[1] = t_max;

    return 2;
}

/* ============================================================================
 * GENERAL QUADRIC
 * ============================================================================ */

int ray_intersect_quadric(const alea_ray_t* ray,
                          const alea_quadric_data_t* q,
                          double* restrict t_out) {
    /*
     * Quadric: Ax^2 + By^2 + Cz^2 + Dxy + Eyz + Fxz + Gx + Hy + Iz + J = 0
     * Coefficients: q->coeffs[0..9] = A, B, C, D, E, F, G, H, I, J
     *
     * Substitute ray P = O + t*D, get quadratic in t.
     */
    double A = q->coeffs[0], B = q->coeffs[1], C = q->coeffs[2];
    double D = q->coeffs[3], E = q->coeffs[4], F = q->coeffs[5];
    double G = q->coeffs[6], H = q->coeffs[7], I = q->coeffs[8];
    double J = q->coeffs[9];

    double ox = ray->ox, oy = ray->oy, oz = ray->oz;
    double dx = ray->dx, dy = ray->dy, dz = ray->dz;

    /* Quadratic coefficients: a*t^2 + b*t + c = 0 */
    double a = A * dx * dx + B * dy * dy + C * dz * dz +
               D * dx * dy + E * dy * dz + F * dx * dz;

    double b = 2.0 * A * ox * dx + 2.0 * B * oy * dy + 2.0 * C * oz * dz +
               D * (ox * dy + oy * dx) + E * (oy * dz + oz * dy) +
               F * (ox * dz + oz * dx) + G * dx + H * dy + I * dz;

    double c = A * ox * ox + B * oy * oy + C * oz * oz +
               D * ox * oy + E * oy * oz + F * ox * oz +
               G * ox + H * oy + I * oz + J;

    return solve_quadratic_stable(a, b, c, t_out);
}

/* ============================================================================
 * RCC MACROBODY (FINITE CYLINDER, ARBITRARY ORIENTATION)
 * ============================================================================ */

int ray_intersect_rcc(const alea_ray_t* ray,
                      const alea_rcc_data_t* rcc,
                      double* restrict t_out) {
    /*
     * RCC: cylinder from base to base+height with given radius.
     * Need to intersect with infinite cylinder along axis, then clip to caps.
     *
     * Transform to cylinder-local coords where axis is along Z, then use
     * standard cylinder intersection.
     */

    /* Axis direction and length */
    double ax = rcc->height_x;
    double ay = rcc->height_y;
    double az = rcc->height_z;
    double height = sqrt(ax * ax + ay * ay + az * az);

    if (height < RAY_EPSILON) {
        return 0;
    }

    /* Normalized axis */
    double inv_h = 1.0 / height;
    ax *= inv_h;
    ay *= inv_h;
    az *= inv_h;

    /* Vector from base to ray origin */
    double dx = ray->ox - rcc->base_x;
    double dy = ray->oy - rcc->base_y;
    double dz = ray->oz - rcc->base_z;

    /* Project onto axis */
    double d_dot_a = dx * ax + dy * ay + dz * az;
    double rd_dot_a = ray->dx * ax + ray->dy * ay + ray->dz * az;

    /* Perpendicular components */
    double perp_dx = dx - d_dot_a * ax;
    double perp_dy = dy - d_dot_a * ay;
    double perp_dz = dz - d_dot_a * az;

    double perp_rdx = ray->dx - rd_dot_a * ax;
    double perp_rdy = ray->dy - rd_dot_a * ay;
    double perp_rdz = ray->dz - rd_dot_a * az;

    /* Quadratic for infinite cylinder */
    double A = perp_rdx * perp_rdx + perp_rdy * perp_rdy + perp_rdz * perp_rdz;
    double B = 2.0 * (perp_dx * perp_rdx + perp_dy * perp_rdy + perp_dz * perp_rdz);
    double C = perp_dx * perp_dx + perp_dy * perp_dy + perp_dz * perp_dz -
               rcc->radius * rcc->radius;

    int count = 0;
    int cyl_hit_count = 0;
    double t_cyl[2];

    cyl_hit_count = solve_quadratic_stable(A, B, C, t_cyl);

    /* Intersect with caps (planes at s=0 and s=height) */
    double t_cap[2];
    int cap_count = 0;

    if (fabs(rd_dot_a) > RAY_EPSILON) {
        /* Bottom cap: s = 0 */
        double t_bot = -d_dot_a / rd_dot_a;
        /* Check if hit is inside radius using cross product (avoids cancellation near axis) */
        double bx = ray->ox + t_bot * ray->dx - rcc->base_x;
        double by = ray->oy + t_bot * ray->dy - rcc->base_y;
        double bz = ray->oz + t_bot * ray->dz - rcc->base_z;
        double cx0 = by * az - bz * ay;
        double cy0 = bz * ax - bx * az;
        double cz0 = bx * ay - by * ax;
        double perp_sq = cx0 * cx0 + cy0 * cy0 + cz0 * cz0;
        if (perp_sq <= rcc->radius * rcc->radius + RAY_EPSILON) {
            t_cap[cap_count++] = t_bot;
        }

        /* Top cap: s = height */
        double t_top = (height - d_dot_a) / rd_dot_a;
        double tx = ray->ox + t_top * ray->dx - rcc->base_x - rcc->height_x;
        double ty = ray->oy + t_top * ray->dy - rcc->base_y - rcc->height_y;
        double tz = ray->oz + t_top * ray->dz - rcc->base_z - rcc->height_z;
        double cx1 = ty * az - tz * ay;
        double cy1 = tz * ax - tx * az;
        double cz1 = tx * ay - ty * ax;
        perp_sq = cx1 * cx1 + cy1 * cy1 + cz1 * cz1;
        if (perp_sq <= rcc->radius * rcc->radius + RAY_EPSILON) {
            t_cap[cap_count++] = t_top;
        }
    }

    /* Collect valid hits: cylinder side hits that are within height bounds */
    for (int i = 0; i < cyl_hit_count; i++) {
        double s = d_dot_a + t_cyl[i] * rd_dot_a;
        if (s >= -RAY_EPSILON && s <= height + RAY_EPSILON) {
            t_out[count++] = t_cyl[i];
        }
    }

    /* Add cap hits */
    for (int i = 0; i < cap_count && count < 2; i++) {
        t_out[count++] = t_cap[i];
    }

    /* Sort */
    if (count == 2 && t_out[0] > t_out[1]) {
        double tmp = t_out[0];
        t_out[0] = t_out[1];
        t_out[1] = tmp;
    }

    return count;
}

/* ============================================================================
 * TRC MACROBODY (TRUNCATED CONE)
 * ============================================================================ */

int ray_intersect_trc(const alea_ray_t* ray,
                      const alea_trc_data_t* trc,
                      double* restrict t_out) {
    /*
     * TRC: frustum from base (radius r1) to top (radius r2).
     * If r1 == r2, it's a cylinder.
     *
     * Similar to RCC but with varying radius along axis.
     */

    /* Axis direction and length */
    double ax = trc->height_x;
    double ay = trc->height_y;
    double az = trc->height_z;
    double height = sqrt(ax * ax + ay * ay + az * az);

    if (height < RAY_EPSILON) {
        return 0;
    }

    /* Normalized axis */
    double inv_h = 1.0 / height;
    ax *= inv_h;
    ay *= inv_h;
    az *= inv_h;

    /* Radius at position s along axis: r(s) = r1 + s*(r2-r1)/height */
    double r1 = trc->base_radius;
    double r2 = trc->top_radius;
    double dr = (r2 - r1) / height;

    /* Vector from base to ray origin */
    double dx = ray->ox - trc->base_x;
    double dy = ray->oy - trc->base_y;
    double dz = ray->oz - trc->base_z;

    /* Project onto axis */
    double d_dot_a = dx * ax + dy * ay + dz * az;
    double rd_dot_a = ray->dx * ax + ray->dy * ay + ray->dz * az;

    /* Perpendicular components */
    double perp_dx = dx - d_dot_a * ax;
    double perp_dy = dy - d_dot_a * ay;
    double perp_dz = dz - d_dot_a * az;

    double perp_rdx = ray->dx - rd_dot_a * ax;
    double perp_rdy = ray->dy - rd_dot_a * ay;
    double perp_rdz = ray->dz - rd_dot_a * az;

    /* Quadratic: |perp(O + t*D)|^2 = r(s)^2 where s = d_dot_a + t*rd_dot_a */
    /* |perp_d + t*perp_rd|^2 = (r1 + s*dr)^2 */

    double A = perp_rdx * perp_rdx + perp_rdy * perp_rdy + perp_rdz * perp_rdz -
               dr * dr * rd_dot_a * rd_dot_a;
    double B = 2.0 * (perp_dx * perp_rdx + perp_dy * perp_rdy + perp_dz * perp_rdz) -
               2.0 * dr * rd_dot_a * (r1 + dr * d_dot_a);
    double C = perp_dx * perp_dx + perp_dy * perp_dy + perp_dz * perp_dz -
               (r1 + dr * d_dot_a) * (r1 + dr * d_dot_a);

    int count = 0;
    double t_cone[2];
    int cone_count = 0;

    cone_count = solve_quadratic_stable(A, B, C, t_cone);

    /* Check cone surface hits are within height bounds */
    for (int i = 0; i < cone_count; i++) {
        double s = d_dot_a + t_cone[i] * rd_dot_a;
        if (s >= -RAY_EPSILON && s <= height + RAY_EPSILON) {
            t_out[count++] = t_cone[i];
        }
    }

    /* Intersect with caps */
    if (fabs(rd_dot_a) > RAY_EPSILON) {
        /* Bottom cap */
        double t_bot = -d_dot_a / rd_dot_a;
        double px = perp_dx + t_bot * perp_rdx;
        double py = perp_dy + t_bot * perp_rdy;
        double pz = perp_dz + t_bot * perp_rdz;
        if (px * px + py * py + pz * pz <= r1 * r1 + RAY_EPSILON) {
            if (count < 2) t_out[count++] = t_bot;
        }

        /* Top cap */
        double t_top = (height - d_dot_a) / rd_dot_a;
        px = perp_dx + t_top * perp_rdx;
        py = perp_dy + t_top * perp_rdy;
        pz = perp_dz + t_top * perp_rdz;
        if (px * px + py * py + pz * pz <= r2 * r2 + RAY_EPSILON) {
            if (count < 2) t_out[count++] = t_top;
        }
    }

    /* Sort */
    if (count == 2 && t_out[0] > t_out[1]) {
        double tmp = t_out[0];
        t_out[0] = t_out[1];
        t_out[1] = tmp;
    }

    return count;
}

static void ray_quadric_from_matrix(const double q[3][3],
                                    double cx, double cy, double cz,
                                    alea_quadric_data_t* out) {
    const double c[3] = {cx, cy, cz};
    double qc[3] = {0.0, 0.0, 0.0};
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            qc[i] += q[i][j] * c[j];
    out->coeffs[0] = q[0][0]; out->coeffs[1] = q[1][1]; out->coeffs[2] = q[2][2];
    out->coeffs[3] = 2.0*q[0][1]; out->coeffs[4] = 2.0*q[1][2];
    out->coeffs[5] = 2.0*q[0][2];
    out->coeffs[6] = -2.0*qc[0]; out->coeffs[7] = -2.0*qc[1];
    out->coeffs[8] = -2.0*qc[2];
    out->coeffs[9] = cx*qc[0] + cy*qc[1] + cz*qc[2] - 1.0;
}

static int ell_quadric(const alea_ell_data_t* ell, alea_quadric_data_t* out) {
    double center[3], u[3], a, b_sq;
    if (ell->major_axis_len < 0.0) {
        center[0]=ell->v1_x; center[1]=ell->v1_y; center[2]=ell->v1_z;
        u[0]=ell->v2_x; u[1]=ell->v2_y; u[2]=ell->v2_z;
        a=sqrt(u[0]*u[0]+u[1]*u[1]+u[2]*u[2]);
        b_sq=ell->major_axis_len*ell->major_axis_len;
    } else {
        center[0]=0.5*(ell->v1_x+ell->v2_x);
        center[1]=0.5*(ell->v1_y+ell->v2_y);
        center[2]=0.5*(ell->v1_z+ell->v2_z);
        u[0]=0.5*(ell->v2_x-ell->v1_x);
        u[1]=0.5*(ell->v2_y-ell->v1_y);
        u[2]=0.5*(ell->v2_z-ell->v1_z);
        a=0.5*ell->major_axis_len;
        const double focal_sq=u[0]*u[0]+u[1]*u[1]+u[2]*u[2];
        b_sq=a*a-focal_sq;
    }
    const double axial_sq=u[0]*u[0]+u[1]*u[1]+u[2]*u[2];
    const double a_sq = a*a;
    if (!(a > 0.0) || !(b_sq > 1e-20) || !isfinite(b_sq)) return 0;
    double q[3][3] = {{0}};
    if (axial_sq <= 1e-20) {
        q[0][0] = q[1][1] = q[2][2] = 1.0/a_sq;
    } else {
        const double inv_c = 1.0/sqrt(axial_sq);
        for (int i = 0; i < 3; i++) u[i] *= inv_c;
        const double inv_b2 = 1.0/b_sq;
        const double delta = 1.0/a_sq - inv_b2;
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
                q[i][j] = (i == j ? inv_b2 : 0.0) + delta*u[i]*u[j];
    }
    ray_quadric_from_matrix(q, center[0], center[1], center[2], out);
    return 1;
}

static int rec_quadric(const alea_rec_data_t* rec, alea_quadric_data_t* out) {
    const double a[3] = {rec->axis1_x, rec->axis1_y, rec->axis1_z};
    const double b[3] = {rec->axis2_x, rec->axis2_y, rec->axis2_z};
    const double a2 = a[0]*a[0] + a[1]*a[1] + a[2]*a[2];
    const double b2 = b[0]*b[0] + b[1]*b[1] + b[2]*b[2];
    if (a2 <= 1e-20 || b2 <= 1e-20) return 0;
    double q[3][3];
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            q[i][j] = a[i]*a[j]/(a2*a2) + b[i]*b[j]/(b2*b2);
    ray_quadric_from_matrix(q, rec->base_x, rec->base_y, rec->base_z, out);
    return 1;
}

int ray_intersect_ell(const alea_ray_t* ray, const alea_ell_data_t* ell,
                      double* restrict t_out) {
    alea_quadric_data_t q;
    if (!ell_quadric(ell, &q)) return 0;
    return ray_intersect_quadric(ray, &q, t_out);
}

static void sort_unique_hits(double* hits, int* count) {
    for (int i = 1; i < *count; i++) {
        double v = hits[i];
        int j = i;
        while (j > 0 && hits[j-1] > v) { hits[j] = hits[j-1]; j--; }
        hits[j] = v;
    }
    int n = 0;
    for (int i = 0; i < *count; i++) {
        if (n == 0 || fabs(hits[i] - hits[n-1]) > 1e-9*(1.0 + fabs(hits[i])))
            hits[n++] = hits[i];
    }
    *count = n;
}

int ray_intersect_rec(const alea_ray_t* ray, const alea_rec_data_t* rec,
                      double* restrict t_out) {
    const double h[3] = {rec->height_x, rec->height_y, rec->height_z};
    const double h2 = h[0]*h[0] + h[1]*h[1] + h[2]*h[2];
    const double a2 = rec->axis1_x*rec->axis1_x + rec->axis1_y*rec->axis1_y +
                      rec->axis1_z*rec->axis1_z;
    const double b2 = rec->axis2_x*rec->axis2_x + rec->axis2_y*rec->axis2_y +
                      rec->axis2_z*rec->axis2_z;
    if (h2 <= 1e-20 || a2 <= 1e-20 || b2 <= 1e-20) return 0;
    alea_quadric_data_t q;
    if (!rec_quadric(rec, &q)) return 0;

    const double op[3] = {ray->ox-rec->base_x, ray->oy-rec->base_y,
                          ray->oz-rec->base_z};
    const double d[3] = {ray->dx, ray->dy, ray->dz};
    const double oh = (op[0]*h[0] + op[1]*h[1] + op[2]*h[2])/h2;
    const double dh = (d[0]*h[0] + d[1]*h[1] + d[2]*h[2])/h2;
    double hits[4];
    int count = 0;
    double side[2];
    int side_count = ray_intersect_quadric(ray, &q, side);
    for (int i = 0; i < side_count; i++) {
        const double s = oh + side[i]*dh;
        if (s >= -RAY_EPSILON && s <= 1.0 + RAY_EPSILON) hits[count++] = side[i];
    }
    if (fabs(dh) > RAY_EPSILON) {
        for (int cap = 0; cap <= 1; cap++) {
            const double t = ((double)cap - oh)/dh;
            const double px = op[0] + t*d[0] - cap*h[0];
            const double py = op[1] + t*d[1] - cap*h[1];
            const double pz = op[2] + t*d[2] - cap*h[2];
            const double ea = (px*rec->axis1_x + py*rec->axis1_y +
                               pz*rec->axis1_z)/a2;
            const double eb = (px*rec->axis2_x + py*rec->axis2_y +
                               pz*rec->axis2_z)/b2;
            if (ea*ea + eb*eb <= 1.0 + RAY_EPSILON) hits[count++] = t;
        }
    }
    sort_unique_hits(hits, &count);
    if (count > 2) count = 2;
    for (int i = 0; i < count; i++) t_out[i] = hits[i];
    return count;
}

/* ============================================================================
 * TORUS (QUARTIC - EXPENSIVE)
 * ============================================================================ */

int ray_intersect_torus(const alea_ray_t* ray,
                        const alea_torus_data_t* torus,
                        double* restrict t_out) {
    /*
     * Torus intersection by quartic equation.
     *
     * Torus implicit equation (Z-axis): (√(x² + y²) - R)² + z² = r²
     * This can be rewritten as: (x² + y² + z² + R² - r²)² = 4R²(x² + y²)
     *
     * Squaring and rearranging gives a quartic in the ray parameter t.
     */

    double R = torus->major_radius;
    double r = torus->minor_radius;
    double R2 = R * R;
    double r2 = r * r;

    /* Transform ray to torus local coordinates (centered, axis-aligned) */
    double ox, oy, oz, dx, dy, dz;

    if (torus->axis == ALEA_AXIS_Z) {
        ox = ray->ox - torus->center_x;
        oy = ray->oy - torus->center_y;
        oz = ray->oz - torus->center_z;
        dx = ray->dx;
        dy = ray->dy;
        dz = ray->dz;
    } else if (torus->axis == ALEA_AXIS_X) {
        /* Rotate so X-axis becomes Z-axis */
        ox = ray->oz - torus->center_z;
        oy = ray->oy - torus->center_y;
        oz = ray->ox - torus->center_x;
        dx = ray->dz;
        dy = ray->dy;
        dz = ray->dx;
    } else { /* ALEA_AXIS_Y */
        /* Rotate so Y-axis becomes Z-axis */
        ox = ray->ox - torus->center_x;
        oy = ray->oz - torus->center_z;
        oz = ray->oy - torus->center_y;
        dx = ray->dx;
        dy = ray->dz;
        dz = ray->dy;
    }

    /* Shift the origin to the closest approach to the torus center so that
     * sum_od = 0. This zeroes the cubic coefficient, so Ferrari's depressed-
     * quartic reduction (whose constant term otherwise cancels catastrophically
     * as ~B^4 terms) stays well-conditioned for rays starting far from the
     * torus. Roots are mapped back by t_shift after solving. */
    double t_shift = -(ox*dx + oy*dy + oz*dz) / (dx*dx + dy*dy + dz*dz);
    ox += t_shift * dx;
    oy += t_shift * dy;
    oz += t_shift * dz;

    /*
     * For ray P + t*D, define:
     * sum_d2 = dx² + dy² + dz²  (should be 1 if normalized)
     * sum_od = ox*dx + oy*dy + oz*dz
     * sum_o2 = ox² + oy² + oz²
     * k = sum_o2 + R² - r²
     *
     * The quartic coefficients come from expanding:
     * (|P + tD|² + R² - r²)² = 4R²((px + t*dx)² + (py + t*dy)²)
     */

    double sum_d2 = dx*dx + dy*dy + dz*dz;
    double sum_od = ox*dx + oy*dy + oz*dz;
    double sum_o2 = ox*ox + oy*oy + oz*oz;
    double k = sum_o2 + R2 - r2;

    /* Quartic coefficients: c4*t⁴ + c3*t³ + c2*t² + c1*t + c0 = 0 */
    double c4 = sum_d2 * sum_d2;
    double c3 = 4.0 * sum_d2 * sum_od;
    double c2 = 2.0 * sum_d2 * k + 4.0 * sum_od * sum_od - 4.0 * R2 * (dx*dx + dy*dy);
    double c1 = 4.0 * sum_od * k - 8.0 * R2 * (ox*dx + oy*dy);
    double c0 = k * k - 4.0 * R2 * (ox*ox + oy*oy);

    /* Solve quartic */
    double roots[4];
    int n = alea_solve_quartic(c4, c3, c2, c1, c0, roots);

    /* Map roots back to the original (unshifted) ray parameter */
    for (int i = 0; i < n; i++)
        roots[i] += t_shift;

    /* Filter positive roots */
    n = alea_filter_positive_roots(roots, n, RAY_EPSILON);

    if (n == 0) {
        return 0;
    }

    /* Return all positive roots (sorted ascending by alea_filter_positive_roots) */
    for (int i = 0; i < n; i++)
        t_out[i] = roots[i];
    return n;
}

typedef struct {
    double a, b, c, d;
} ray_halfspace_t;

static void orient_ray_halfspace(ray_halfspace_t* p,
                                 double x, double y, double z) {
    if (p->a*x + p->b*y + p->c*z + p->d > 0.0) {
        p->a = -p->a; p->b = -p->b; p->c = -p->c; p->d = -p->d;
    }
}

static int ray_clip_halfspaces(const alea_ray_t* ray,
                               const ray_halfspace_t* planes, int plane_count,
                               double* restrict t_out) {
    double enter = -INFINITY, leave = INFINITY;
    for (int i = 0; i < plane_count; i++) {
        const ray_halfspace_t* p = &planes[i];
        const double value = p->a*ray->ox + p->b*ray->oy + p->c*ray->oz + p->d;
        const double slope = p->a*ray->dx + p->b*ray->dy + p->c*ray->dz;
        if (fabs(slope) <= RAY_EPSILON) {
            if (value > RAY_EPSILON) return 0;
            continue;
        }
        const double t = -value/slope;
        if (slope < 0.0) enter = fmax(enter, t);
        else leave = fmin(leave, t);
        if (enter > leave + RAY_EPSILON) return 0;
    }
    if (!isfinite(enter) || !isfinite(leave)) return 0;
    t_out[0] = enter;
    t_out[1] = leave;
    return 2;
}

static int ray_intersect_box_general(const alea_ray_t* ray,
                                     const alea_box_general_data_t* box,
                                     double* restrict t_out) {
    const double v1[3] = {box->v1_x, box->v1_y, box->v1_z};
    const double v2[3] = {box->v2_x, box->v2_y, box->v2_z};
    const double v3[3] = {box->v3_x, box->v3_y, box->v3_z};
    double n[3][3] = {
        {v2[1]*v3[2]-v2[2]*v3[1], v2[2]*v3[0]-v2[0]*v3[2], v2[0]*v3[1]-v2[1]*v3[0]},
        {v3[1]*v1[2]-v3[2]*v1[1], v3[2]*v1[0]-v3[0]*v1[2], v3[0]*v1[1]-v3[1]*v1[0]},
        {v1[1]*v2[2]-v1[2]*v2[1], v1[2]*v2[0]-v1[0]*v2[2], v1[0]*v2[1]-v1[1]*v2[0]}
    };
    const double c[3] = {box->corner_x, box->corner_y, box->corner_z};
    const double inside[3] = {c[0]+0.5*(v1[0]+v2[0]+v3[0]),
                              c[1]+0.5*(v1[1]+v2[1]+v3[1]),
                              c[2]+0.5*(v1[2]+v2[2]+v3[2])};
    ray_halfspace_t p[6];
    for (int i = 0; i < 3; i++) {
        const double len2 = n[i][0]*n[i][0] + n[i][1]*n[i][1] + n[i][2]*n[i][2];
        if (len2 <= 1e-20) return 0;
        const double* edge = i == 0 ? v1 : (i == 1 ? v2 : v3);
        p[2*i] = (ray_halfspace_t){n[i][0], n[i][1], n[i][2],
                                  -(n[i][0]*c[0]+n[i][1]*c[1]+n[i][2]*c[2])};
        const double top[3] = {c[0]+edge[0], c[1]+edge[1], c[2]+edge[2]};
        p[2*i+1] = (ray_halfspace_t){n[i][0], n[i][1], n[i][2],
                                    -(n[i][0]*top[0]+n[i][1]*top[1]+n[i][2]*top[2])};
        orient_ray_halfspace(&p[2*i], inside[0], inside[1], inside[2]);
        orient_ray_halfspace(&p[2*i+1], inside[0], inside[1], inside[2]);
    }
    return ray_clip_halfspaces(ray, p, 6, t_out);
}

static int ray_intersect_wed(const alea_ray_t* ray, const alea_wed_data_t* w,
                             double* restrict t_out) {
    const double c[3] = {w->vertex_x, w->vertex_y, w->vertex_z};
    const double v1[3] = {w->v1_x, w->v1_y, w->v1_z};
    const double v2[3] = {w->v2_x, w->v2_y, w->v2_z};
    const double v3[3] = {w->v3_x, w->v3_y, w->v3_z};
    const double inside[3] = {c[0]+(v1[0]+v2[0])/3.0+0.5*v3[0],
                              c[1]+(v1[1]+v2[1])/3.0+0.5*v3[1],
                              c[2]+(v1[2]+v2[2])/3.0+0.5*v3[2]};
    const double points[5][3] = {
        {c[0], c[1], c[2]}, {c[0]+v3[0], c[1]+v3[1], c[2]+v3[2]},
        {c[0], c[1], c[2]}, {c[0], c[1], c[2]},
        {c[0]+v1[0], c[1]+v1[1], c[2]+v1[2]}
    };
    double normals[5][3] = {
        {v1[1]*v2[2]-v1[2]*v2[1], v1[2]*v2[0]-v1[0]*v2[2], v1[0]*v2[1]-v1[1]*v2[0]},
        {v1[1]*v2[2]-v1[2]*v2[1], v1[2]*v2[0]-v1[0]*v2[2], v1[0]*v2[1]-v1[1]*v2[0]},
        {v3[1]*v1[2]-v3[2]*v1[1], v3[2]*v1[0]-v3[0]*v1[2], v3[0]*v1[1]-v3[1]*v1[0]},
        {v2[1]*v3[2]-v2[2]*v3[1], v2[2]*v3[0]-v2[0]*v3[2], v2[0]*v3[1]-v2[1]*v3[0]},
        {(v2[1]-v1[1])*v3[2]-(v2[2]-v1[2])*v3[1],
         (v2[2]-v1[2])*v3[0]-(v2[0]-v1[0])*v3[2],
         (v2[0]-v1[0])*v3[1]-(v2[1]-v1[1])*v3[0]}
    };
    ray_halfspace_t p[5];
    for (int i = 0; i < 5; i++) {
        const double len2 = normals[i][0]*normals[i][0] + normals[i][1]*normals[i][1] + normals[i][2]*normals[i][2];
        if (len2 <= 1e-20) return 0;
        p[i] = (ray_halfspace_t){normals[i][0], normals[i][1], normals[i][2],
                                -(normals[i][0]*points[i][0]+normals[i][1]*points[i][1]+normals[i][2]*points[i][2])};
        orient_ray_halfspace(&p[i], inside[0], inside[1], inside[2]);
    }
    return ray_clip_halfspaces(ray, p, 5, t_out);
}

static int ray_intersect_rhp(const alea_ray_t* ray, const alea_rhp_data_t* r,
                             double* restrict t_out) {
    const double hlen = sqrt(r->height_x*r->height_x+r->height_y*r->height_y+r->height_z*r->height_z);
    if (hlen <= 1e-20) return 0;
    const double axis[3] = {r->height_x/hlen, r->height_y/hlen, r->height_z/hlen};
    const double c[3] = {r->base_x, r->base_y, r->base_z};
    ray_halfspace_t p[8] = {
        {-axis[0],-axis[1],-axis[2], axis[0]*c[0]+axis[1]*c[1]+axis[2]*c[2]},
        { axis[0], axis[1], axis[2],-(axis[0]*(c[0]+r->height_x)+axis[1]*(c[1]+r->height_y)+axis[2]*(c[2]+r->height_z))}
    };
    double radial[3][3];
    if (!alea_rhp_resolve_radials(r, radial)) return 0;
    for (int i = 0; i < 3; i++) {
        const double len = sqrt(radial[i][0]*radial[i][0]+radial[i][1]*radial[i][1]+radial[i][2]*radial[i][2]);
        if (len <= 1e-20) return 0;
        const double nx=radial[i][0]/len, ny=radial[i][1]/len, nz=radial[i][2]/len;
        p[2+2*i] = (ray_halfspace_t){nx,ny,nz,-(nx*c[0]+ny*c[1]+nz*c[2])-len};
        p[3+2*i] = (ray_halfspace_t){-nx,-ny,-nz,nx*c[0]+ny*c[1]+nz*c[2]-len};
    }
    return ray_clip_halfspaces(ray, p, 8, t_out);
}

static int ray_intersect_arb(const alea_ray_t* ray, const alea_arb_data_t* arb,
                             double* restrict t_out) {
    if (arb->num_corners < 4 || arb->num_corners > 8 ||
        arb->num_faces < 4 || arb->num_faces > 6) return 0;
    double inside[3] = {0};
    for (int i=0;i<arb->num_corners;i++) for (int j=0;j<3;j++) inside[j]+=arb->corners[i][j];
    for (int j=0;j<3;j++) inside[j]/=arb->num_corners;
    ray_halfspace_t p[6];
    for (int f=0;f<arb->num_faces;f++) {
        int i0=arb->faces[f][0]-1, i1=arb->faces[f][1]-1, i2=arb->faces[f][2]-1;
        if (i0<0||i1<0||i2<0||i0>=arb->num_corners||i1>=arb->num_corners||i2>=arb->num_corners) return 0;
        const double* a=arb->corners[i0]; const double* b=arb->corners[i1]; const double* c=arb->corners[i2];
        const double e1[3]={b[0]-a[0],b[1]-a[1],b[2]-a[2]};
        const double e2[3]={c[0]-a[0],c[1]-a[1],c[2]-a[2]};
        const double nx=e1[1]*e2[2]-e1[2]*e2[1], ny=e1[2]*e2[0]-e1[0]*e2[2], nz=e1[0]*e2[1]-e1[1]*e2[0];
        if (nx*nx+ny*ny+nz*nz <= 1e-20) return 0;
        p[f]=(ray_halfspace_t){nx,ny,nz,-(nx*a[0]+ny*a[1]+nz*a[2])};
        orient_ray_halfspace(&p[f],inside[0],inside[1],inside[2]);
    }
    return ray_clip_halfspaces(ray,p,arb->num_faces,t_out);
}

/* ============================================================================
 * DISPATCH
 * ============================================================================ */

int ray_intersect_primitive(const alea_ray_t* ray,
                            alea_primitive_type_t type,
                            const alea_primitive_data_t* data,
                            double* restrict t_out) {
    switch (type) {
        case ALEA_PRIMITIVE_PLANE:
            return ray_intersect_plane(ray, &data->plane, t_out, NULL, NULL, NULL);

        case ALEA_PRIMITIVE_SPHERE:
        case ALEA_PRIMITIVE_SPH:
            return ray_intersect_sphere(ray, &data->sphere, t_out);

        case ALEA_PRIMITIVE_CYLINDER_X:
            return ray_intersect_cylinder_x(ray, &data->cyl_x, t_out);

        case ALEA_PRIMITIVE_CYLINDER_Y:
            return ray_intersect_cylinder_y(ray, &data->cyl_y, t_out);

        case ALEA_PRIMITIVE_CYLINDER_Z:
            return ray_intersect_cylinder_z(ray, &data->cyl_z, t_out);

        case ALEA_PRIMITIVE_CONE_X:
            return ray_intersect_cone_x(ray, &data->cone_x, t_out);

        case ALEA_PRIMITIVE_CONE_Y:
            return ray_intersect_cone_y(ray, &data->cone_y, t_out);

        case ALEA_PRIMITIVE_CONE_Z:
            return ray_intersect_cone_z(ray, &data->cone_z, t_out);

        case ALEA_PRIMITIVE_RPP:
            return ray_intersect_box(ray, &data->box, t_out);

        case ALEA_PRIMITIVE_QUADRIC:
            return ray_intersect_quadric(ray, &data->quadric, t_out);

        case ALEA_PRIMITIVE_TORUS_X:
        case ALEA_PRIMITIVE_TORUS_Y:
        case ALEA_PRIMITIVE_TORUS_Z:
            return ray_intersect_torus(ray, &data->torus, t_out);

        case ALEA_PRIMITIVE_RCC:
            return ray_intersect_rcc(ray, &data->rcc, t_out);

        case ALEA_PRIMITIVE_BOX:
            return ray_intersect_box_general(ray, &data->box_general, t_out);

        case ALEA_PRIMITIVE_TRC:
            return ray_intersect_trc(ray, &data->trc, t_out);

        case ALEA_PRIMITIVE_ELL:
            return ray_intersect_ell(ray, &data->ell, t_out);

        case ALEA_PRIMITIVE_REC:
            return ray_intersect_rec(ray, &data->rec, t_out);

        case ALEA_PRIMITIVE_WED:
            return ray_intersect_wed(ray, &data->wed, t_out);

        case ALEA_PRIMITIVE_RHP:
            return ray_intersect_rhp(ray, &data->rhp, t_out);

        case ALEA_PRIMITIVE_ARB:
            return ray_intersect_arb(ray, &data->arb, t_out);

        default:
            return 0;  /* Unknown primitive */
    }
}

/* ============================================================================
 * SURFACE NORMALS
 * ============================================================================ */

void primitive_normal_at(alea_primitive_type_t type,
                         const alea_primitive_data_t* data,
                         double px, double py, double pz,
                         double* nx, double* ny, double* nz) {
    switch (type) {
        case ALEA_PRIMITIVE_PLANE: {
            /* Planes are pre-normalized by alea_canonicalize_primitive() */
            const alea_plane_data_t* p = &data->plane;
            *nx = p->a;
            *ny = p->b;
            *nz = p->c;
            break;
        }

        case ALEA_PRIMITIVE_SPHERE:
        case ALEA_PRIMITIVE_SPH: {
            const alea_sphere_data_t* s = &data->sphere;
            double dx = px - s->center_x;
            double dy = py - s->center_y;
            double dz = pz - s->center_z;
            double len = sqrt(dx * dx + dy * dy + dz * dz);
            if (len > RAY_EPSILON) {
                *nx = dx / len;
                *ny = dy / len;
                *nz = dz / len;
            } else {
                *nx = 0; *ny = 0; *nz = 1;
            }
            break;
        }

        case ALEA_PRIMITIVE_CYLINDER_Z: {
            const alea_cylinder_z_data_t* c = &data->cyl_z;
            double dx = px - c->center_x;
            double dy = py - c->center_y;
            double len = sqrt(dx * dx + dy * dy);
            if (len > RAY_EPSILON) {
                *nx = dx / len;
                *ny = dy / len;
                *nz = 0;
            } else {
                *nx = 1; *ny = 0; *nz = 0;
            }
            break;
        }

        case ALEA_PRIMITIVE_CYLINDER_X: {
            const alea_cylinder_x_data_t* c = &data->cyl_x;
            double dy = py - c->center_y;
            double dz = pz - c->center_z;
            double len = sqrt(dy * dy + dz * dz);
            if (len > RAY_EPSILON) {
                *nx = 0;
                *ny = dy / len;
                *nz = dz / len;
            } else {
                *nx = 0; *ny = 1; *nz = 0;
            }
            break;
        }

        case ALEA_PRIMITIVE_CYLINDER_Y: {
            const alea_cylinder_y_data_t* c = &data->cyl_y;
            double dx = px - c->center_x;
            double dz = pz - c->center_z;
            double len = sqrt(dx * dx + dz * dz);
            if (len > RAY_EPSILON) {
                *nx = dx / len;
                *ny = 0;
                *nz = dz / len;
            } else {
                *nx = 1; *ny = 0; *nz = 0;
            }
            break;
        }

        case ALEA_PRIMITIVE_CONE_Z: {
            const alea_cone_z_data_t* c = &data->cone_z;
            double dx = px - c->apex_x;
            double dy = py - c->apex_y;
            double dz = pz - c->apex_z;
            /* Gradient of (x-ax)^2 + (y-ay)^2 - k^2*(z-az)^2 */
            double gx = 2.0 * dx;
            double gy = 2.0 * dy;
            double gz = -2.0 * c->tan_angle_sq * dz;
            double len = sqrt(gx * gx + gy * gy + gz * gz);
            if (len > RAY_EPSILON) {
                *nx = gx / len; *ny = gy / len; *nz = gz / len;
            } else {
                *nx = 0; *ny = 0; *nz = 1;
            }
            break;
        }

        case ALEA_PRIMITIVE_CONE_X: {
            const alea_cone_x_data_t* c = &data->cone_x;
            double dx = px - c->apex_x;
            double dy = py - c->apex_y;
            double dz = pz - c->apex_z;
            /* Gradient of (y-ay)^2 + (z-az)^2 - k^2*(x-ax)^2 */
            double gx = -2.0 * c->tan_angle_sq * dx;
            double gy = 2.0 * dy;
            double gz = 2.0 * dz;
            double len = sqrt(gx * gx + gy * gy + gz * gz);
            if (len > RAY_EPSILON) {
                *nx = gx / len; *ny = gy / len; *nz = gz / len;
            } else {
                *nx = 1; *ny = 0; *nz = 0;
            }
            break;
        }

        case ALEA_PRIMITIVE_CONE_Y: {
            const alea_cone_y_data_t* c = &data->cone_y;
            double dx = px - c->apex_x;
            double dy = py - c->apex_y;
            double dz = pz - c->apex_z;
            /* Gradient of (x-ax)^2 + (z-az)^2 - k^2*(y-ay)^2 */
            double gx = 2.0 * dx;
            double gy = -2.0 * c->tan_angle_sq * dy;
            double gz = 2.0 * dz;
            double len = sqrt(gx * gx + gy * gy + gz * gz);
            if (len > RAY_EPSILON) {
                *nx = gx / len; *ny = gy / len; *nz = gz / len;
            } else {
                *nx = 0; *ny = 1; *nz = 0;
            }
            break;
        }

        case ALEA_PRIMITIVE_RPP: {
            const alea_box_data_t* b = &data->box;
            /* Determine closest face and return its outward normal */
            double dx_min = fabs(px - b->min_x);
            double dx_max = fabs(px - b->max_x);
            double dy_min = fabs(py - b->min_y);
            double dy_max = fabs(py - b->max_y);
            double dz_min = fabs(pz - b->min_z);
            double dz_max = fabs(pz - b->max_z);
            double dmin = dx_min;
            *nx = -1; *ny = 0; *nz = 0;
            if (dx_max < dmin) { dmin = dx_max; *nx = 1; *ny = 0; *nz = 0; }
            if (dy_min < dmin) { dmin = dy_min; *nx = 0; *ny = -1; *nz = 0; }
            if (dy_max < dmin) { dmin = dy_max; *nx = 0; *ny = 1; *nz = 0; }
            if (dz_min < dmin) { dmin = dz_min; *nx = 0; *ny = 0; *nz = -1; }
            if (dz_max < dmin) { *nx = 0; *ny = 0; *nz = 1; }
            break;
        }

        case ALEA_PRIMITIVE_QUADRIC: {
            const alea_quadric_data_t* q = &data->quadric;
            double A = q->coeffs[0], B = q->coeffs[1], C = q->coeffs[2];
            double D = q->coeffs[3], E = q->coeffs[4], F = q->coeffs[5];
            double G = q->coeffs[6], H = q->coeffs[7], I = q->coeffs[8];
            /* Gradient: (2Ax + Dy + Fz + G, 2By + Dx + Ez + H, 2Cz + Ey + Fx + I) */
            double gx = 2.0*A*px + D*py + F*pz + G;
            double gy = 2.0*B*py + D*px + E*pz + H;
            double gz = 2.0*C*pz + E*py + F*px + I;
            double len = sqrt(gx * gx + gy * gy + gz * gz);
            if (len > RAY_EPSILON) {
                *nx = gx / len; *ny = gy / len; *nz = gz / len;
            } else {
                *nx = 0; *ny = 0; *nz = 1;
            }
            break;
        }

        case ALEA_PRIMITIVE_TORUS_X:
        case ALEA_PRIMITIVE_TORUS_Y:
        case ALEA_PRIMITIVE_TORUS_Z: {
            /* Analytical gradient of torus implicit equation */
            const alea_torus_data_t* t = &data->torus;
            double R = t->major_radius;
            double r = t->minor_radius;
            double R2 = R * R, r2 = r * r;
            /* Transform to local coords (torus centered at origin, axis along Z) */
            double lx, ly, lz;
            if (t->axis == ALEA_AXIS_Z) {
                lx = px - t->center_x; ly = py - t->center_y; lz = pz - t->center_z;
            } else if (t->axis == ALEA_AXIS_X) {
                lx = pz - t->center_z; ly = py - t->center_y; lz = px - t->center_x;
            } else {
                lx = px - t->center_x; ly = pz - t->center_z; lz = py - t->center_y;
            }
            /* f(x,y,z) = (x^2+y^2+z^2+R^2-r^2)^2 - 4R^2(x^2+y^2) */
            double sum2 = lx*lx + ly*ly + lz*lz;
            double k = sum2 + R2 - r2;
            double gx_l = 4.0 * k * lx - 8.0 * R2 * lx;
            double gy_l = 4.0 * k * ly - 8.0 * R2 * ly;
            double gz_l = 4.0 * k * lz;
            /* Transform gradient back to world coords */
            double gx, gy, gz;
            if (t->axis == ALEA_AXIS_Z) {
                gx = gx_l; gy = gy_l; gz = gz_l;
            } else if (t->axis == ALEA_AXIS_X) {
                gx = gz_l; gy = gy_l; gz = gx_l;
            } else {
                gx = gx_l; gy = gz_l; gz = gy_l;
            }
            double len = sqrt(gx * gx + gy * gy + gz * gz);
            if (len > RAY_EPSILON) {
                *nx = gx / len; *ny = gy / len; *nz = gz / len;
            } else {
                *nx = 0; *ny = 0; *nz = 1;
            }
            break;
        }

        case ALEA_PRIMITIVE_RCC:
        case ALEA_PRIMITIVE_BOX:
        case ALEA_PRIMITIVE_TRC:
        case ALEA_PRIMITIVE_ELL:
        case ALEA_PRIMITIVE_REC:
        case ALEA_PRIMITIVE_WED:
        case ALEA_PRIMITIVE_RHP:
        case ALEA_PRIMITIVE_ARB: {
            /* The macrobody evaluators are signed classifiers whose local
             * gradient gives the outward normal, including caps and faces. */
            const double scale = 1.0 + fmax(fabs(px), fmax(fabs(py), fabs(pz)));
            const double h = 1e-6 * scale;
            const double gx = alea_primitive_eval(type, data, px+h, py, pz) -
                              alea_primitive_eval(type, data, px-h, py, pz);
            const double gy = alea_primitive_eval(type, data, px, py+h, pz) -
                              alea_primitive_eval(type, data, px, py-h, pz);
            const double gz = alea_primitive_eval(type, data, px, py, pz+h) -
                              alea_primitive_eval(type, data, px, py, pz-h);
            const double len = sqrt(gx*gx + gy*gy + gz*gz);
            if (len > RAY_EPSILON) {
                *nx = gx/len; *ny = gy/len; *nz = gz/len;
            } else {
                *nx = 0; *ny = 0; *nz = 1;
            }
            break;
        }

        default:
            /* Fallback for types without dedicated normal computation */
            *nx = 0;
            *ny = 0;
            *nz = 1;
            break;
    }
}
