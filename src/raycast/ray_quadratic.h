// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

#ifndef ALEA_RAY_QUADRATIC_H
#define ALEA_RAY_QUADRATIC_H

#include "raycast.h"
#include <float.h>
#include <math.h>

/**
 * Numerically stable quadratic solver using Vieta's formula.
 * Avoids catastrophic cancellation in (-b ± sqrt(disc)) / (2a).
 * Returns number of real roots (0, 1, or 2) in t_out[].
 */
static inline int solve_quadratic_stable(double a, double b, double c,
                                         double* t_out) {
    /* A small quadratic coefficient is not necessarily negligible: for a ray
     * almost parallel to a cylinder/cone axis it represents a real, distant
     * crossing.  Only reduce an actually linear equation. */
    if (a == 0.0) {
        if (b == 0.0) return 0;
        t_out[0] = -c / b;
        return 1;
    }

    double ac4 = 4.0 * a * c;
    double discriminant = fma(b, b, -ac4);
    double discriminant_scale = fabs(b * b) + fabs(ac4);
    double discriminant_tol = 16.0 * DBL_EPSILON * discriminant_scale;
    if (discriminant < -discriminant_tol) return 0;
    if (discriminant < 0) discriminant = 0;

    double sqrt_disc = sqrt(discriminant);
    double q = -0.5 * (b + copysign(sqrt_disc, b));

    if (q == 0.0) {
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

typedef struct { double a, b, c; } ray_quadratic_t;

/* Keep the operation order identical for one surface and for independent
 * lanes. The stride selects the scalar or coefficient-major layout. */
static inline ray_quadratic_t quadric_polynomial(const alea_ray_t* ray,
                                                  const double* coeffs, size_t stride) {
    /*
     * Quadric: Ax^2 + By^2 + Cz^2 + Dxy + Eyz + Fxz + Gx + Hy + Iz + J = 0
     * Coefficients: A, B, C, D, E, F, G, H, I, J
     *
     * Substitute ray P = O + t*D, get quadratic in t.
     */
    double A = coeffs[0], B = coeffs[stride], C = coeffs[2 * stride];
    double D = coeffs[3 * stride], E = coeffs[4 * stride], F = coeffs[5 * stride];
    double G = coeffs[6 * stride], H = coeffs[7 * stride], I = coeffs[8 * stride];
    double J = coeffs[9 * stride];

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

    return (ray_quadratic_t){a, b, c};
}

#endif
