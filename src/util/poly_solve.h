// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

#ifndef ALEA_POLY_SOLVE_H
#define ALEA_POLY_SOLVE_H


/**
 * @file poly_solve.h
 * @brief Polynomial root solvers for ray intersection
 *
 * Provides solvers for quadratic, cubic, and quartic equations.
 * Quartic solver is needed for torus ray intersection.
 */

/**
 * @brief Solve quadratic equation: ax² + bx + c = 0
 *
 * @param a, b, c Coefficients (a != 0)
 * @param roots Output array for up to 2 real roots (sorted ascending)
 * @return Number of real roots found (0, 1, or 2)
 */
int alea_solve_quadratic(double a, double b, double c, double roots[2]);

/**
 * @brief Solve cubic equation: ax³ + bx² + cx + d = 0
 *
 * Uses Cardano's formula with careful handling of edge cases.
 *
 * @param a, b, c, d Coefficients (a != 0)
 * @param roots Output array for up to 3 real roots (sorted ascending)
 * @return Number of real roots found (1, 2, or 3)
 */
int alea_solve_cubic(double a, double b, double c, double d, double roots[3]);

/**
 * @brief Solve quartic equation: ax⁴ + bx³ + cx² + dx + e = 0
 *
 * Uses Ferrari's method to reduce to cubic, then quadratics.
 *
 * @param a, b, c, d, e Coefficients (a != 0)
 * @param roots Output array for up to 4 real roots (sorted ascending)
 * @return Number of real roots found (0, 1, 2, 3, or 4)
 */
int alea_solve_quartic(double a, double b, double c, double d, double e, double roots[4]);

/**
 * @brief Filter and sort positive roots
 *
 * Given roots array, filter to only positive values and sort ascending.
 *
 * @param roots Input/output array
 * @param count Number of roots in input
 * @param min_t Minimum value (typically small epsilon)
 * @return Number of positive roots
 */
int alea_filter_positive_roots(double* roots, int count, double min_t);


#endif /* ALEA_POLY_SOLVE_H */
