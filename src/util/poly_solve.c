// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/**
 * @file poly_solve.c
 * @brief Polynomial root solvers implementation
 *
 * Implements quadratic, cubic, and quartic equation solvers.
 * Based on standard analytical formulas with numerical stability improvements.
 */

#include "poly_solve.h"
#include <math.h>
#include <float.h>
#include <string.h>
#include "util/math.h"

#define EPSILON 1e-12
#define CBRT(x) ((x) >= 0 ? pow((x), 1.0/3.0) : -pow(-(x), 1.0/3.0))

/* Sort roots in ascending order */
static void sort_roots(double* roots, int n) {
    for (int i = 0; i < n - 1; i++) {
        for (int j = i + 1; j < n; j++) {
            if (roots[j] < roots[i]) {
                double tmp = roots[i];
                roots[i] = roots[j];
                roots[j] = tmp;
            }
        }
    }
}

int alea_solve_quadratic(double a, double b, double c, double roots[2]) {
    if (fabs(a) < EPSILON) {
        /* Linear equation: bx + c = 0 */
        if (fabs(b) < EPSILON) {
            return 0;
        }
        roots[0] = -c / b;
        return 1;
    }

    double discriminant = b * b - 4 * a * c;

    if (discriminant < -EPSILON) {
        return 0;  /* No real roots */
    }

    if (discriminant < EPSILON) {
        /* One repeated root */
        roots[0] = -b / (2 * a);
        return 1;
    }

    /* Two distinct roots - use numerically stable formula */
    double sqrt_disc = sqrt(discriminant);
    double q = (b >= 0) ? -0.5 * (b + sqrt_disc) : -0.5 * (b - sqrt_disc);

    roots[0] = q / a;
    roots[1] = c / q;

    sort_roots(roots, 2);
    return 2;
}

int alea_solve_cubic(double a, double b, double c, double d, double roots[3]) {
    if (fabs(a) < EPSILON) {
        /* Reduce to quadratic */
        return alea_solve_quadratic(b, c, d, roots);
    }

    /* Normalize: x³ + px² + qx + r = 0 */
    double p = b / a;
    double q = c / a;
    double r = d / a;

    /* Substitution x = t - p/3 to get depressed cubic: t³ + At + B = 0 */
    double p2 = p * p;
    double A = q - p2 / 3.0;
    double B = (2.0 * p2 * p - 9.0 * p * q + 27.0 * r) / 27.0;

    double A3 = A * A * A;
    double discriminant = B * B / 4.0 + A3 / 27.0;

    int num_roots;
    double offset = -p / 3.0;

    if (discriminant > EPSILON) {
        /* One real root */
        double sqrt_disc = sqrt(discriminant);
        double u = CBRT(-B / 2.0 + sqrt_disc);
        double v = CBRT(-B / 2.0 - sqrt_disc);
        roots[0] = u + v + offset;
        num_roots = 1;
    } else if (discriminant < -EPSILON) {
        /* Three distinct real roots (casus irreducibilis) */
        double m = 2.0 * sqrt(-A / 3.0);
        double cos_arg = 3.0 * B / (A * m);
        if (cos_arg > 1.0) cos_arg = 1.0;
        else if (cos_arg < -1.0) cos_arg = -1.0;
        double theta = acos(cos_arg) / 3.0;
        roots[0] = m * cos(theta) + offset;
        roots[1] = m * cos(theta - 2.0 * M_PI / 3.0) + offset;
        roots[2] = m * cos(theta - 4.0 * M_PI / 3.0) + offset;
        num_roots = 3;
    } else {
        /* One single and one double root */
        double u = CBRT(-B / 2.0);
        roots[0] = 2.0 * u + offset;
        roots[1] = -u + offset;
        num_roots = 2;
    }

    sort_roots(roots, num_roots);
    return num_roots;
}

/* Evaluate monic quartic x⁴ + Bx³ + Cx² + Dx + E (Horner). */
static inline double quartic_eval(double B, double C, double D, double E,
                                  double x) {
    return (((x + B) * x + C) * x + D) * x + E;
}

/* Evaluate derivative 4x³ + 3Bx² + 2Cx + D (Horner). */
static inline double quartic_deriv(double B, double C, double D, double x) {
    return ((4.0 * x + 3.0 * B) * x + 2.0 * C) * x + D;
}

/* Absolute-value evaluation |x|⁴ + |B||x|³ + ... — the natural error scale of
 * quartic_eval at x; used to decide whether a residual is "numerically zero". */
static inline double quartic_eval_scale(double B, double C, double D, double E,
                                        double x) {
    double ax = fabs(x);
    return (((ax + fabs(B)) * ax + fabs(C)) * ax + fabs(D)) * ax + fabs(E);
}

/**
 * Refine a root of the monic quartic inside a bracket [lo, hi] where
 * f(lo) and f(hi) have opposite signs. Safeguarded Newton: take the Newton
 * step when it stays inside the bracket, otherwise bisect. Converges
 * quadratically near the root and never escapes the bracket.
 */
static double quartic_refine(double B, double C, double D, double E,
                             double lo, double hi, double flo) {
    double x = 0.5 * (lo + hi);
    for (int iter = 0; iter < 100; iter++) {
        double f = quartic_eval(B, C, D, E, x);
        if (f == 0.0) return x;
        /* Shrink bracket around the sign change */
        if ((f > 0) == (flo > 0)) {
            lo = x;
        } else {
            hi = x;
        }
        double fp = quartic_deriv(B, C, D, x);
        double x_next;
        if (fp != 0.0) {
            x_next = x - f / fp;
            if (!(x_next > lo && x_next < hi)) {
                x_next = 0.5 * (lo + hi);  /* Newton left bracket: bisect */
            }
        } else {
            x_next = 0.5 * (lo + hi);
        }
        if (fabs(x_next - x) <= 4.0 * DBL_EPSILON * (fabs(x_next) + 1.0)) {
            return x_next;
        }
        x = x_next;
    }
    return x;
}

int alea_solve_quartic(double a, double b, double c, double d, double e, double roots[4]) {
    if (fabs(a) < EPSILON) {
        /* Reduce to cubic */
        return alea_solve_cubic(b, c, d, e, roots);
    }

    /* Normalize: x⁴ + Bx³ + Cx² + Dx + E = 0 */
    double B = b / a;
    double C = c / a;
    double D = d / a;
    double E = e / a;

    /* Robust real-root isolation. Closed-form (Ferrari) solutions are
     * numerically fragile: the depressed-quartic reduction cancels
     * catastrophically for large |B|, and near-biquadratic inputs make the
     * resolvent-cubic factorization pick wrong branches (observed as lost
     * torus intersections). Instead, isolate roots between the critical
     * points of f (roots of the derivative cubic, where f is monotonic on
     * each interval) and refine each bracketed root with safeguarded Newton.
     */

    /* Critical points: roots of f'(x) = 4x³ + 3Bx² + 2Cx + D */
    double crit[3];
    int n_crit = alea_solve_cubic(4.0, 3.0 * B, 2.0 * C, D, crit);

    /* Polish critical points (the cubic solver itself is closed-form) */
    for (int i = 0; i < n_crit; i++) {
        double x = crit[i];
        for (int it = 0; it < 3; it++) {
            double g = quartic_deriv(B, C, D, x);
            double gp = (12.0 * x + 6.0 * B) * x + 2.0 * C;
            if (fabs(gp) > 1e-30) x -= g / gp;
        }
        crit[i] = x;
    }
    sort_roots(crit, n_crit);

    /* Outer bound beyond all real roots (Fujiwara), also covering the
     * critical points so every interval endpoint has a defined sign. */
    double M = fabs(B);
    double t2 = sqrt(fabs(C));       if (t2 > M) M = t2;
    double t3 = cbrt(fabs(D));       if (t3 > M) M = t3;
    double t4 = sqrt(sqrt(fabs(E))); if (t4 > M) M = t4;
    M = 2.0 * M + 1.0;
    if (n_crit > 0) {
        double cmax = fabs(crit[0]);
        if (fabs(crit[n_crit - 1]) > cmax) cmax = fabs(crit[n_crit - 1]);
        if (cmax + 1.0 > M) M = cmax + 1.0;
    }

    /* Interval endpoints: -M, critical points, +M (f is monotonic between
     * consecutive endpoints, so each sign change brackets exactly one root). */
    double pts[5];
    int n_pts = 0;
    pts[n_pts++] = -M;
    for (int i = 0; i < n_crit; i++) {
        if (crit[i] > -M && crit[i] < M) pts[n_pts++] = crit[i];
    }
    pts[n_pts++] = M;

    int num_roots = 0;
    double f_prev = quartic_eval(B, C, D, E, pts[0]);
    for (int i = 1; i < n_pts; i++) {
        double f_here = quartic_eval(B, C, D, E, pts[i]);

        if ((f_prev > 0) != (f_here > 0)) {
            /* Sign change: exactly one root in (pts[i-1], pts[i]) */
            roots[num_roots++] = quartic_refine(B, C, D, E,
                                                pts[i - 1], pts[i], f_prev);
        } else if (i < n_pts - 1 &&
                   fabs(f_here) <= 4.0 * DBL_EPSILON *
                       quartic_eval_scale(B, C, D, E, pts[i])) {
            /* Interior critical point with f ≈ 0: double (tangent) root.
             * Report it twice so crossing parity is preserved. */
            if (num_roots <= 2) {
                roots[num_roots++] = pts[i];
                roots[num_roots++] = pts[i];
            }
            /* Skip the sign-change test on the adjacent interval: treat the
             * residual sign at this point as unreliable. */
            f_here = f_prev;
        }
        f_prev = f_here;
    }

    sort_roots(roots, num_roots);
    return num_roots;
}

int alea_filter_positive_roots(double* roots, int count, double min_t) {
    int pos_count = 0;
    for (int i = 0; i < count; i++) {
        if (roots[i] > min_t) {
            roots[pos_count++] = roots[i];
        }
    }
    sort_roots(roots, pos_count);
    return pos_count;
}
