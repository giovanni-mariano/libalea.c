// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

#ifndef SLICE_EPSILON_H
#define SLICE_EPSILON_H

/**
 * @file slice_epsilon.h
 * @brief Named epsilon constants for slice module
 *
 * Centralizes all floating-point tolerance constants used across
 * surface-plane intersection, curve evaluation, and spatial queries.
 */

/* Parallel/perpendicular direction checks: cos_angle > 1-TOL, len < TOL */
#define SLICE_PARALLEL_TOL 1e-10

/* Quadratic discriminant classification and near-zero coefficient guards */
#define SLICE_DISCRIMINANT_TOL 1e-10

/* Tangent detection: |distance - radius| < TOL */
#define SLICE_TANGENT_TOL 1e-10

/* Edge parameter bounds: t < -TOL or t > 1+TOL */
#define SLICE_EDGE_TOL 1e-10

/* Spatial index query bbox margin */
#define SLICE_SPATIAL_TOL 1e-10

/* Transform comparison tolerance for surface deduplication */
#define SLICE_TRANSFORM_TOL 1e-9

#endif /* SLICE_EPSILON_H */
