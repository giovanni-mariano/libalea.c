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
#include "util/poly_solve.h"
#include <math.h>
#include <float.h>

#define EPSILON 1e-10

/* ============================================================================
 * PLANE
 * ============================================================================ */

int ray_intersect_plane(const alea_ray_t* ray,
                        const alea_plane_data_t* plane,
                        double* t_out,
                        double* nx, double* ny, double* nz) {
    /*
     * Plane: ax + by + cz + d = 0
     * Ray: P = O + t*D
     * Substitute: a(ox + t*dx) + b(oy + t*dy) + c(oz + t*dz) + d = 0
     * Solve: t = -(a*ox + b*oy + c*oz + d) / (a*dx + b*dy + c*dz)
     */
    double denom = plane->a * ray->dx + plane->b * ray->dy + plane->c * ray->dz;

    if (fabs(denom) < EPSILON) {
        return 0;  /* Ray parallel to plane */
    }

    double numer = plane->a * ray->ox + plane->b * ray->oy +
                   plane->c * ray->oz + plane->d;

    t_out[0] = -numer / denom;

    /* Normal points in direction of (a,b,c) */
    if (nx) {
        double len = sqrt(plane->a * plane->a + plane->b * plane->b +
                         plane->c * plane->c);
        *nx = plane->a / len;
        *ny = plane->b / len;
        *nz = plane->c / len;
    }

    return 1;
}

/* ============================================================================
 * SPHERE
 * ============================================================================ */

int ray_intersect_sphere(const alea_ray_t* ray,
                         const alea_sphere_data_t* sphere,
                         double* t_out) {
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

    double b = lx * ray->dx + ly * ray->dy + lz * ray->dz;  /* L.D */
    double c = lx * lx + ly * ly + lz * lz - sphere->radius * sphere->radius;

    double discriminant = b * b - c;

    /* Allow small negative discriminant due to floating-point rounding errors */
    if (discriminant < -1e-12) {
        return 0;
    }
    if (discriminant < 0) discriminant = 0;

    double sqrt_disc = sqrt(discriminant);
    t_out[0] = -b - sqrt_disc;
    t_out[1] = -b + sqrt_disc;

    return 2;
}

/* ============================================================================
 * AXIS-ALIGNED CYLINDERS (INFINITE)
 * ============================================================================ */

int ray_intersect_cylinder_z(const alea_ray_t* ray,
                             const alea_cylinder_z_data_t* cyl,
                             double* t_out) {
    /*
     * Cylinder along Z: (x - cx)^2 + (y - cy)^2 = r^2
     * Substitute ray, solve quadratic in XY plane only
     */
    double px = ray->ox - cyl->center_x;
    double py = ray->oy - cyl->center_y;

    double a = ray->dx * ray->dx + ray->dy * ray->dy;
    double b = 2.0 * (px * ray->dx + py * ray->dy);
    double c = px * px + py * py - cyl->radius * cyl->radius;

    if (fabs(a) < EPSILON) {
        return 0;  /* Ray parallel to cylinder axis */
    }

    double discriminant = b * b - 4.0 * a * c;

    /* Allow small negative discriminant due to floating-point rounding errors */
    if (discriminant < -1e-12) {
        return 0;
    }
    if (discriminant < 0) discriminant = 0;

    double sqrt_disc = sqrt(discriminant);
    double inv_2a = 0.5 / a;

    double t_raw[2];
    t_raw[0] = (-b - sqrt_disc) * inv_2a;
    t_raw[1] = (-b + sqrt_disc) * inv_2a;

    t_out[0] = t_raw[0];
    t_out[1] = t_raw[1];
    return 2;
}

int ray_intersect_cylinder_x(const alea_ray_t* ray,
                             const alea_cylinder_x_data_t* cyl,
                             double* t_out) {
    /* Cylinder along X: (y - cy)^2 + (z - cz)^2 = r^2 */
    double py = ray->oy - cyl->center_y;
    double pz = ray->oz - cyl->center_z;

    double a = ray->dy * ray->dy + ray->dz * ray->dz;
    double b = 2.0 * (py * ray->dy + pz * ray->dz);
    double c = py * py + pz * pz - cyl->radius * cyl->radius;

    if (fabs(a) < EPSILON) {
        return 0;
    }

    double discriminant = b * b - 4.0 * a * c;

    /* Allow small negative discriminant due to floating-point rounding errors */
    if (discriminant < -1e-12) {
        return 0;
    }
    if (discriminant < 0) discriminant = 0;

    double sqrt_disc = sqrt(discriminant);
    double inv_2a = 0.5 / a;

    double t_raw[2];
    t_raw[0] = (-b - sqrt_disc) * inv_2a;
    t_raw[1] = (-b + sqrt_disc) * inv_2a;

    t_out[0] = t_raw[0];
    t_out[1] = t_raw[1];
    return 2;
}

int ray_intersect_cylinder_y(const alea_ray_t* ray,
                             const alea_cylinder_y_data_t* cyl,
                             double* t_out) {
    /* Cylinder along Y: (x - cx)^2 + (z - cz)^2 = r^2 */
    double px = ray->ox - cyl->center_x;
    double pz = ray->oz - cyl->center_z;

    double a = ray->dx * ray->dx + ray->dz * ray->dz;
    double b = 2.0 * (px * ray->dx + pz * ray->dz);
    double c = px * px + pz * pz - cyl->radius * cyl->radius;

    if (fabs(a) < EPSILON) {
        return 0;
    }

    double discriminant = b * b - 4.0 * a * c;

    /* Allow small negative discriminant due to floating-point rounding errors */
    if (discriminant < -1e-12) {
        return 0;
    }
    if (discriminant < 0) discriminant = 0;

    double sqrt_disc = sqrt(discriminant);
    double inv_2a = 0.5 / a;

    double t_raw[2];
    t_raw[0] = (-b - sqrt_disc) * inv_2a;
    t_raw[1] = (-b + sqrt_disc) * inv_2a;

    t_out[0] = t_raw[0];
    t_out[1] = t_raw[1];
    return 2;
}

/* ============================================================================
 * AXIS-ALIGNED CONES
 * ============================================================================ */

int ray_intersect_cone_z(const alea_ray_t* ray,
                         const alea_cone_z_data_t* cone,
                         double* t_out) {
    /*
     * Cone along Z: (x - ax)^2 + (y - ay)^2 = k^2 * (z - az)^2
     * where k^2 = tan^2(half-angle)
     *
     * Substituting ray and rearranging gives quadratic in t.
     */
    double px = ray->ox - cone->apex_x;
    double py = ray->oy - cone->apex_y;
    double pz = ray->oz - cone->apex_z;

    double k2 = cone->tan_angle_sq;

    double a = ray->dx * ray->dx + ray->dy * ray->dy - k2 * ray->dz * ray->dz;
    double b = 2.0 * (px * ray->dx + py * ray->dy - k2 * pz * ray->dz);
    double c = px * px + py * py - k2 * pz * pz;

    double t_raw[2];
    int raw_count = 0;

    if (fabs(a) < EPSILON) {
        /* Degenerate case: ray along cone surface */
        if (fabs(b) < EPSILON) {
            return 0;
        }
        t_raw[0] = -c / b;
        raw_count = 1;
    } else {
        double discriminant = b * b - 4.0 * a * c;
        /* Allow small negative discriminant due to floating-point rounding errors */
        if (discriminant < -1e-12) {
            return 0;
        }
        if (discriminant < 0) discriminant = 0;
        double sqrt_disc = sqrt(discriminant);
        double inv_2a = 0.5 / a;
        t_raw[0] = (-b - sqrt_disc) * inv_2a;
        t_raw[1] = (-b + sqrt_disc) * inv_2a;
        raw_count = 2;
    }

    /* Filter by sheet_selection */
    int count = 0;
    for (int i = 0; i < raw_count; i++) {
        double axis_val = pz + t_raw[i] * ray->dz;  /* z relative to apex */
        /* Sheet selection: +1 = positive nappe, -1 = negative nappe, 0 = both */
        if (cone->sheet_selection > 0 && axis_val < -EPSILON) continue;
        if (cone->sheet_selection < 0 && axis_val > EPSILON) continue;
        t_out[count++] = t_raw[i];
    }

    return count;
}

int ray_intersect_cone_x(const alea_ray_t* ray,
                         const alea_cone_x_data_t* cone,
                         double* t_out) {
    /* Cone along X: (y - ay)^2 + (z - az)^2 = k^2 * (x - ax)^2 */
    double px = ray->ox - cone->apex_x;
    double py = ray->oy - cone->apex_y;
    double pz = ray->oz - cone->apex_z;

    double k2 = cone->tan_angle_sq;

    double a = ray->dy * ray->dy + ray->dz * ray->dz - k2 * ray->dx * ray->dx;
    double b = 2.0 * (py * ray->dy + pz * ray->dz - k2 * px * ray->dx);
    double c = py * py + pz * pz - k2 * px * px;

    double t_raw[2];
    int raw_count = 0;

    if (fabs(a) < EPSILON) {
        if (fabs(b) < EPSILON) {
            return 0;
        }
        t_raw[0] = -c / b;
        raw_count = 1;
    } else {
        double discriminant = b * b - 4.0 * a * c;
        /* Allow small negative discriminant due to floating-point rounding errors */
        if (discriminant < -1e-12) {
            return 0;
        }
        if (discriminant < 0) discriminant = 0;
        double sqrt_disc = sqrt(discriminant);
        double inv_2a = 0.5 / a;
        t_raw[0] = (-b - sqrt_disc) * inv_2a;
        t_raw[1] = (-b + sqrt_disc) * inv_2a;
        raw_count = 2;
    }

    /* Filter by sheet_selection */
    int count = 0;
    for (int i = 0; i < raw_count; i++) {
        double axis_val = px + t_raw[i] * ray->dx;  /* x relative to apex */
        if (cone->sheet_selection > 0 && axis_val < -EPSILON) continue;
        if (cone->sheet_selection < 0 && axis_val > EPSILON) continue;
        t_out[count++] = t_raw[i];
    }

    return count;
}

int ray_intersect_cone_y(const alea_ray_t* ray,
                         const alea_cone_y_data_t* cone,
                         double* t_out) {
    /* Cone along Y: (x - ax)^2 + (z - az)^2 = k^2 * (y - ay)^2 */
    double px = ray->ox - cone->apex_x;
    double py = ray->oy - cone->apex_y;
    double pz = ray->oz - cone->apex_z;

    double k2 = cone->tan_angle_sq;

    double a = ray->dx * ray->dx + ray->dz * ray->dz - k2 * ray->dy * ray->dy;
    double b = 2.0 * (px * ray->dx + pz * ray->dz - k2 * py * ray->dy);
    double c = px * px + pz * pz - k2 * py * py;

    double t_raw[2];
    int raw_count = 0;

    if (fabs(a) < EPSILON) {
        if (fabs(b) < EPSILON) {
            return 0;
        }
        t_raw[0] = -c / b;
        raw_count = 1;
    } else {
        double discriminant = b * b - 4.0 * a * c;
        /* Allow small negative discriminant due to floating-point rounding errors */
        if (discriminant < -1e-12) {
            return 0;
        }
        if (discriminant < 0) discriminant = 0;
        double sqrt_disc = sqrt(discriminant);
        double inv_2a = 0.5 / a;
        t_raw[0] = (-b - sqrt_disc) * inv_2a;
        t_raw[1] = (-b + sqrt_disc) * inv_2a;
        raw_count = 2;
    }

    /* Filter by sheet_selection */
    int count = 0;
    for (int i = 0; i < raw_count; i++) {
        double axis_val = py + t_raw[i] * ray->dy;  /* y relative to apex */
        if (cone->sheet_selection > 0 && axis_val < -EPSILON) continue;
        if (cone->sheet_selection < 0 && axis_val > EPSILON) continue;
        t_out[count++] = t_raw[i];
    }

    return count;
}

/* ============================================================================
 * AXIS-ALIGNED BOX
 * ============================================================================ */

int ray_intersect_box(const alea_ray_t* ray,
                      const alea_box_data_t* box,
                      double* t_out) {
    /*
     * Slab method: intersect ray with each pair of parallel planes,
     * find overlapping interval.
     */
    double t_min = -DBL_MAX;
    double t_max = DBL_MAX;

    /* X slabs */
    if (fabs(ray->dx) < EPSILON) {
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
    if (fabs(ray->dy) < EPSILON) {
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
    if (fabs(ray->dz) < EPSILON) {
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
                          double* t_out) {
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

    if (fabs(a) < EPSILON) {
        /* Linear case */
        if (fabs(b) < EPSILON) {
            return 0;
        }
        t_out[0] = -c / b;
        return 1;
    }

    double discriminant = b * b - 4.0 * a * c;

    /* Allow small negative discriminant due to floating-point rounding errors */
    if (discriminant < -1e-12) {
        return 0;
    }
    if (discriminant < 0) discriminant = 0;

    double sqrt_disc = sqrt(discriminant);
    double inv_2a = 0.5 / a;

    t_out[0] = (-b - sqrt_disc) * inv_2a;
    t_out[1] = (-b + sqrt_disc) * inv_2a;

    return 2;
}

/* ============================================================================
 * RCC MACROBODY (FINITE CYLINDER, ARBITRARY ORIENTATION)
 * ============================================================================ */

int ray_intersect_rcc(const alea_ray_t* ray,
                      const alea_rcc_data_t* rcc,
                      double* t_out) {
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

    if (height < EPSILON) {
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

    if (fabs(A) > EPSILON) {
        double discriminant = B * B - 4.0 * A * C;
        if (discriminant >= 0) {
            double sqrt_disc = sqrt(discriminant);
            double inv_2A = 0.5 / A;
            t_cyl[0] = (-B - sqrt_disc) * inv_2A;
            t_cyl[1] = (-B + sqrt_disc) * inv_2A;
            cyl_hit_count = 2;
        }
    }

    /* Intersect with caps (planes at s=0 and s=height) */
    double t_cap[2];
    int cap_count = 0;

    if (fabs(rd_dot_a) > EPSILON) {
        /* Bottom cap: s = 0 */
        double t_bot = -d_dot_a / rd_dot_a;
        /* Check if hit is inside radius */
        double px = ray->ox + t_bot * ray->dx - rcc->base_x;
        double py = ray->oy + t_bot * ray->dy - rcc->base_y;
        double pz = ray->oz + t_bot * ray->dz - rcc->base_z;
        double perp_sq = px * px + py * py + pz * pz -
                        (px * ax + py * ay + pz * az) *
                        (px * ax + py * ay + pz * az);
        if (perp_sq <= rcc->radius * rcc->radius + EPSILON) {
            t_cap[cap_count++] = t_bot;
        }

        /* Top cap: s = height */
        double t_top = (height - d_dot_a) / rd_dot_a;
        px = ray->ox + t_top * ray->dx - rcc->base_x - rcc->height_x;
        py = ray->oy + t_top * ray->dy - rcc->base_y - rcc->height_y;
        pz = ray->oz + t_top * ray->dz - rcc->base_z - rcc->height_z;
        perp_sq = px * px + py * py + pz * pz -
                 (px * ax + py * ay + pz * az) *
                 (px * ax + py * ay + pz * az);
        if (perp_sq <= rcc->radius * rcc->radius + EPSILON) {
            t_cap[cap_count++] = t_top;
        }
    }

    /* Collect valid hits: cylinder side hits that are within height bounds */
    for (int i = 0; i < cyl_hit_count; i++) {
        double s = d_dot_a + t_cyl[i] * rd_dot_a;
        if (s >= -EPSILON && s <= height + EPSILON) {
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
                      double* t_out) {
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

    if (height < EPSILON) {
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

    if (fabs(A) > EPSILON) {
        double discriminant = B * B - 4.0 * A * C;
        if (discriminant >= 0) {
            double sqrt_disc = sqrt(discriminant);
            double inv_2A = 0.5 / A;
            t_cone[0] = (-B - sqrt_disc) * inv_2A;
            t_cone[1] = (-B + sqrt_disc) * inv_2A;
            cone_count = 2;
        }
    } else if (fabs(B) > EPSILON) {
        t_cone[0] = -C / B;
        cone_count = 1;
    }

    /* Check cone surface hits are within height bounds */
    for (int i = 0; i < cone_count; i++) {
        double s = d_dot_a + t_cone[i] * rd_dot_a;
        if (s >= -EPSILON && s <= height + EPSILON) {
            t_out[count++] = t_cone[i];
        }
    }

    /* Intersect with caps */
    if (fabs(rd_dot_a) > EPSILON) {
        /* Bottom cap */
        double t_bot = -d_dot_a / rd_dot_a;
        double px = perp_dx + t_bot * perp_rdx;
        double py = perp_dy + t_bot * perp_rdy;
        double pz = perp_dz + t_bot * perp_rdz;
        if (px * px + py * py + pz * pz <= r1 * r1 + EPSILON) {
            if (count < 2) t_out[count++] = t_bot;
        }

        /* Top cap */
        double t_top = (height - d_dot_a) / rd_dot_a;
        px = perp_dx + t_top * perp_rdx;
        py = perp_dy + t_top * perp_rdy;
        pz = perp_dz + t_top * perp_rdz;
        if (px * px + py * py + pz * pz <= r2 * r2 + EPSILON) {
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

/* ============================================================================
 * TORUS (QUARTIC - EXPENSIVE)
 * ============================================================================ */

int ray_intersect_torus(const alea_ray_t* ray,
                        const alea_torus_data_t* torus,
                        double* t_out) {
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

    /* Filter positive roots */
    n = alea_filter_positive_roots(roots, n, EPSILON);

    if (n == 0) {
        return 0;
    }

    /* Return all positive roots (sorted ascending by alea_filter_positive_roots) */
    for (int i = 0; i < n; i++)
        t_out[i] = roots[i];
    return n;
}

/* ============================================================================
 * DISPATCH
 * ============================================================================ */

int ray_intersect_primitive(const alea_ray_t* ray,
                            alea_primitive_type_t type,
                            const alea_primitive_data_t* data,
                            double* t_out) {
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

        case ALEA_PRIMITIVE_TRC:
            return ray_intersect_trc(ray, &data->trc, t_out);

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
            const alea_plane_data_t* p = &data->plane;
            double len = sqrt(p->a * p->a + p->b * p->b + p->c * p->c);
            *nx = p->a / len;
            *ny = p->b / len;
            *nz = p->c / len;
            break;
        }

        case ALEA_PRIMITIVE_SPHERE:
        case ALEA_PRIMITIVE_SPH: {
            const alea_sphere_data_t* s = &data->sphere;
            double dx = px - s->center_x;
            double dy = py - s->center_y;
            double dz = pz - s->center_z;
            double len = sqrt(dx * dx + dy * dy + dz * dz);
            if (len > EPSILON) {
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
            if (len > EPSILON) {
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
            if (len > EPSILON) {
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
            if (len > EPSILON) {
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
            if (len > EPSILON) {
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
            if (len > EPSILON) {
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
            if (len > EPSILON) {
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
            if (len > EPSILON) {
                *nx = gx / len; *ny = gy / len; *nz = gz / len;
            } else {
                *nx = 0; *ny = 0; *nz = 1;
            }
            break;
        }

        case ALEA_PRIMITIVE_TORUS_X:
        case ALEA_PRIMITIVE_TORUS_Y:
        case ALEA_PRIMITIVE_TORUS_Z: {
            /* Numerical gradient for torus (analytical is complex) */
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
            if (len > EPSILON) {
                *nx = gx / len; *ny = gy / len; *nz = gz / len;
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
