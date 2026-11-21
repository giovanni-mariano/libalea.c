// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
// SPDX-License-Identifier: MPL-2.0

/* Keep the product and addition residuals when evaluating F(O + tD).
 * Near coincident quadrics even rounding the point to three doubles can
 * lose the sign. This cold path uses ordinary binary64 and C99 fma. */
#include "raycast/ray_wide.h"

static double refine_quadric_crossing(const alea_ray_t* ray, const double k[10], double t) {
    ray_quadratic_t poly = quadric_polynomial(ray, k, 1);
    double slope = fma(2 * poly.a, t, poly.b);
    if (!isfinite(slope) || slope == 0) return t;
    ray_wide_t parameter = {t, 0};
    ray_wide_t x = ray_wide_term(ray->dx, parameter, (ray_wide_t){ray->ox, 0});
    ray_wide_t y = ray_wide_term(ray->dy, parameter, (ray_wide_t){ray->oy, 0});
    ray_wide_t z = ray_wide_term(ray->dz, parameter, (ray_wide_t){ray->oz, 0});
    /* x(Ax + Dy + Fz + G) + y(By + Ez + H) + z(Cz + I) + J. */
    ray_wide_t fx = ray_wide_term(k[5], z, (ray_wide_t){k[6], 0});
    fx = ray_wide_term(k[0], x, ray_wide_term(k[3], y, fx));
    ray_wide_t fy = ray_wide_term(k[1], y, ray_wide_term(k[4], z, (ray_wide_t){k[7], 0}));
    ray_wide_t fz = ray_wide_term(k[2], z, (ray_wide_t){k[8], 0});
    ray_wide_t f = ray_wide_madd(z, fz, (ray_wide_t){k[9], 0});
    f = ray_wide_madd(x, fx, ray_wide_madd(y, fy, f));
    double step = (f.hi + f.lo) / slope;
    double refined = t - step;
    /* Stay on the same simple root. Near a double root the Newton step is
     * unreliable. Its two solutions are separated by |F'(t)/a|. */
    if (isfinite(refined) && fabs(4 * poly.a * step) < fabs(slope)) return refined;
    return t;
}
