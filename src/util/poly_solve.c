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

/**
 * Polish roots via Newton-Raphson iteration (Horner evaluation).
 * 2 iterations per root; guards against near-degenerate derivative.
 */
static void newton_polish_quartic(double a, double b, double c, double d, double e,
                                  double* roots, int count) {
    for (int i = 0; i < count; i++) {
        double x = roots[i];
        for (int iter = 0; iter < 2; iter++) {
            /* f(x) = ax^4 + bx^3 + cx^2 + dx + e  (Horner) */
            double f = (((a * x + b) * x + c) * x + d) * x + e;
            /* f'(x) = 4ax^3 + 3bx^2 + 2cx + d  (Horner) */
            double fp = ((4.0 * a * x + 3.0 * b) * x + 2.0 * c) * x + d;
            if (fabs(fp) > 1e-30) {
                x -= f / fp;
            }
        }
        roots[i] = x;
    }
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

    /* Check for biquadratic: x⁴ + Cx² + E = 0 (B = D = 0) */
    if (fabs(B) < EPSILON && fabs(D) < EPSILON) {
        double quad_roots[2];
        int n = alea_solve_quadratic(1.0, C, E, quad_roots);
        int num_roots = 0;
        for (int i = 0; i < n; i++) {
            if (quad_roots[i] >= 0) {
                double sq = sqrt(quad_roots[i]);
                roots[num_roots++] = -sq;
                roots[num_roots++] = sq;
            } else if (quad_roots[i] > -EPSILON) {
                roots[num_roots++] = 0;
            }
        }
        sort_roots(roots, num_roots);
        return num_roots;
    }

    /* Ferrari's method: reduce to depressed quartic then resolve cubic */

    /* Substitution x = t - B/4 to get: t⁴ + pt² + qt + r = 0 */
    double B2 = B * B;
    double B3 = B2 * B;
    double B4 = B2 * B2;

    double p = C - 3.0 * B2 / 8.0;
    double q = B3 / 8.0 - B * C / 2.0 + D;
    double r = -3.0 * B4 / 256.0 + B2 * C / 16.0 - B * D / 4.0 + E;

    double offset = -B / 4.0;

    /* Check for depressed biquadratic (q = 0) */
    if (fabs(q) < EPSILON) {
        double quad_roots[2];
        int n = alea_solve_quadratic(1.0, p, r, quad_roots);
        int num_roots = 0;
        for (int i = 0; i < n; i++) {
            if (quad_roots[i] >= 0) {
                double sq = sqrt(quad_roots[i]);
                roots[num_roots++] = -sq + offset;
                roots[num_roots++] = sq + offset;
            } else if (quad_roots[i] > -EPSILON) {
                roots[num_roots++] = offset;
            }
        }
        sort_roots(roots, num_roots);
        return num_roots;
    }

    /* Resolve cubic: 8m³ + 8pm² + (2p² - 8r)m - q² = 0 */
    double cubic_roots[3];
    int n_cubic = alea_solve_cubic(
        8.0,
        8.0 * p,
        2.0 * p * p - 8.0 * r,
        -q * q,
        cubic_roots
    );

    /* Find a suitable root m (we need 2m + p > 0 for real factorization) */
    double m = 0;
    for (int i = n_cubic - 1; i >= 0; i--) {
        if (2.0 * cubic_roots[i] + p > EPSILON) {
            m = cubic_roots[i];
            break;
        }
    }

    if (2.0 * m + p <= EPSILON) {
        /* Fallback: use largest root and hope for the best */
        m = cubic_roots[n_cubic - 1];
        if (2.0 * m + p < 0) {
            /* No real solution possible with this method */
            return 0;
        }
    }

    /* Factor into two quadratics:
     * t² + sqrt(2m+p)*t + (m + sqrt(m² - r)) = 0
     * t² - sqrt(2m+p)*t + (m - sqrt(m² - r)) = 0
     */
    double sqrt_2mp = sqrt(2.0 * m + p);
    double m2_r = m * m - r;
    double sqrt_m2r = (m2_r >= 0) ? sqrt(m2_r) : 0;

    /* Determine signs based on q */
    double sign = (q >= 0) ? 1.0 : -1.0;

    double q1_roots[2], q2_roots[2];
    int n1 = alea_solve_quadratic(1.0, sqrt_2mp, m - sign * sqrt_m2r, q1_roots);
    int n2 = alea_solve_quadratic(1.0, -sqrt_2mp, m + sign * sqrt_m2r, q2_roots);

    int num_roots = 0;
    for (int i = 0; i < n1; i++) {
        roots[num_roots++] = q1_roots[i] + offset;
    }
    for (int i = 0; i < n2; i++) {
        roots[num_roots++] = q2_roots[i] + offset;
    }

    /* Polish roots via Newton-Raphson for improved accuracy */
    newton_polish_quartic(a, b, c, d, e, roots, num_roots);

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
