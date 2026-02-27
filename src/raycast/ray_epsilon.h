// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

#ifndef RAY_EPSILON_H
#define RAY_EPSILON_H

/**
 * @file ray_epsilon.h
 * @brief Named epsilon constants for raycast module
 *
 * Centralizes all floating-point tolerance constants used across
 * ray intersection, BVH traversal, and segment building.
 */

/* Parallel-ray checks, hit dedup, axis-aligned direction checks */
#define RAY_EPSILON 1e-10

/* Quadratic/quartic near-tangent discriminant clamping */
#define DISCRIMINANT_TOL 1e-12

/* Hit deduplication after sorting (same t AND same surface_id) */
#define DEDUP_EPSILON 1e-10

/* Offset for cell-finding sample point after surface crossing */
#define SURFACE_SAMPLE_OFFSET 0.001

/* Offset for cell-finding sample point after lattice DDA boundary crossing */
#define DDA_SAMPLE_OFFSET 0.01

#endif /* RAY_EPSILON_H */
