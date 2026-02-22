// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/*
 * test_poly.c - Unit tests for polynomial solvers
 */

#include "alea_test.h"
#include "util/poly_solve.h"

#define EPS 1e-8

/* ------------------------------------------------------------------------- */
/* Quadratic Tests                                                            */
/* ------------------------------------------------------------------------- */

TEST(quadratic_two_roots) {
    double roots[2];

    /* x^2 - 5x + 6 = 0  =>  x = 2, 3 */
    int n = alea_solve_quadratic(1, -5, 6, roots);
    ASSERT_EQ(n, 2);
    ASSERT_NEAR(roots[0], 2.0, EPS);
    ASSERT_NEAR(roots[1], 3.0, EPS);

    /* x^2 - 1 = 0  =>  x = -1, 1 */
    n = alea_solve_quadratic(1, 0, -1, roots);
    ASSERT_EQ(n, 2);
    ASSERT_NEAR(roots[0], -1.0, EPS);
    ASSERT_NEAR(roots[1], 1.0, EPS);
}

TEST(quadratic_one_root) {
    double roots[2];

    /* x^2 - 2x + 1 = 0  =>  x = 1 (double root) */
    int n = alea_solve_quadratic(1, -2, 1, roots);
    ASSERT_EQ(n, 1);
    ASSERT_NEAR(roots[0], 1.0, EPS);
}

TEST(quadratic_no_roots) {
    double roots[2];

    /* x^2 + 1 = 0  =>  no real roots */
    int n = alea_solve_quadratic(1, 0, 1, roots);
    ASSERT_EQ(n, 0);
}

/* ------------------------------------------------------------------------- */
/* Cubic Tests                                                                */
/* ------------------------------------------------------------------------- */

TEST(cubic_three_roots) {
    double roots[3];

    /* x^3 - 6x^2 + 11x - 6 = 0  =>  x = 1, 2, 3 */
    int n = alea_solve_cubic(1, -6, 11, -6, roots);
    ASSERT_EQ(n, 3);
    ASSERT_NEAR(roots[0], 1.0, EPS);
    ASSERT_NEAR(roots[1], 2.0, EPS);
    ASSERT_NEAR(roots[2], 3.0, EPS);

    /* x^3 - 3x = 0  =>  x = -sqrt(3), 0, sqrt(3) */
    n = alea_solve_cubic(1, 0, -3, 0, roots);
    ASSERT_EQ(n, 3);
    ASSERT_NEAR(roots[0], -sqrt(3.0), EPS);
    ASSERT_NEAR(roots[1], 0.0, EPS);
    ASSERT_NEAR(roots[2], sqrt(3.0), EPS);
}

TEST(cubic_one_root) {
    double roots[3];

    /* x^3 + 1 = 0  =>  x = -1 */
    int n = alea_solve_cubic(1, 0, 0, 1, roots);
    ASSERT_EQ(n, 1);
    ASSERT_NEAR(roots[0], -1.0, EPS);
}

/* ------------------------------------------------------------------------- */
/* Quartic Tests                                                              */
/* ------------------------------------------------------------------------- */

TEST(quartic_four_roots) {
    double roots[4];

    /* (x-1)(x-2)(x-3)(x-4) = x^4 - 10x^3 + 35x^2 - 50x + 24 = 0 */
    int n = alea_solve_quartic(1, -10, 35, -50, 24, roots);
    ASSERT_EQ(n, 4);
    ASSERT_NEAR(roots[0], 1.0, EPS);
    ASSERT_NEAR(roots[1], 2.0, EPS);
    ASSERT_NEAR(roots[2], 3.0, EPS);
    ASSERT_NEAR(roots[3], 4.0, EPS);
}

TEST(quartic_biquadratic) {
    double roots[4];

    /* x^4 - 5x^2 + 4 = 0  =>  (x^2-1)(x^2-4) = 0  =>  x = -2, -1, 1, 2 */
    int n = alea_solve_quartic(1, 0, -5, 0, 4, roots);
    ASSERT_EQ(n, 4);
    ASSERT_NEAR(roots[0], -2.0, EPS);
    ASSERT_NEAR(roots[1], -1.0, EPS);
    ASSERT_NEAR(roots[2], 1.0, EPS);
    ASSERT_NEAR(roots[3], 2.0, EPS);
}

TEST(quartic_two_roots) {
    double roots[4];

    /* x^4 - 1 = 0  =>  x = -1, 1 (and +/-i) */
    int n = alea_solve_quartic(1, 0, 0, 0, -1, roots);
    ASSERT_EQ(n, 2);
    ASSERT_NEAR(roots[0], -1.0, EPS);
    ASSERT_NEAR(roots[1], 1.0, EPS);
}

TEST(quartic_no_roots) {
    double roots[4];

    /* x^4 + 1 = 0  =>  no real roots */
    int n = alea_solve_quartic(1, 0, 0, 0, 1, roots);
    ASSERT_EQ(n, 0);
}

TEST(quartic_torus_like) {
    /* Quartic arising from torus intersection (realistic coefficients) */
    double roots[4];

    /* x^4 - 2x^2 + 0.5 = 0 */
    /* Roots: +/- sqrt(1 +/- sqrt(0.5)) */
    int n = alea_solve_quartic(1, 0, -2, 0, 0.5, roots);
    ASSERT_EQ(n, 4);

    double expected[] = {
        -sqrt(1 + sqrt(0.5)),
        -sqrt(1 - sqrt(0.5)),
        sqrt(1 - sqrt(0.5)),
        sqrt(1 + sqrt(0.5))
    };

    for (int i = 0; i < 4; i++) {
        ASSERT_NEAR(roots[i], expected[i], EPS);
    }
}

/* ------------------------------------------------------------------------- */
/* Filter Tests                                                               */
/* ------------------------------------------------------------------------- */

TEST(filter_positive) {
    double roots[] = {-2.0, -0.5, 0.0001, 1.0, 3.0};
    int n = alea_filter_positive_roots(roots, 5, 1e-6);

    ASSERT_EQ(n, 3);
    ASSERT_NEAR(roots[0], 0.0001, EPS);
    ASSERT_NEAR(roots[1], 1.0, EPS);
    ASSERT_NEAR(roots[2], 3.0, EPS);
}

TEST_MAIN()
