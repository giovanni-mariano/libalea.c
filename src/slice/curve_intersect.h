// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/**
 * @file curve_intersect.h
 * @brief Analytical surface-plane intersection for 2D slice rendering
 *
 * Computes the exact 2D curves where geometric surfaces intersect a slice plane.
 * These curves can then be rasterized for pixel-perfect boundary rendering.
 */

#ifndef CURVE_INTERSECT_H
#define CURVE_INTERSECT_H

#include "alea_types.h"
#include "alea_slice.h"
#include "util/alea_vec.h"
#include <stddef.h>
#include <stdbool.h>


/* alea_slice_plane_t is defined in alea_slice.h (included above).
 * All internal code uses alea_slice_plane_t directly — no separate
 * internal type needed since the structs were byte-identical. */

/* ============================================================================
 * 2D CURVE TYPES
 * ============================================================================ */

/**
 * @brief Types of 2D curves that can result from plane-surface intersection
 *
 * Extended version of the enum in alea_slice.h (public header).
 * If already defined via alea_slice.h, we just add the extra values.
 */
#ifndef ALEA_CURVE_TYPE_DEFINED
#define ALEA_CURVE_TYPE_DEFINED
typedef enum {
    ALEA_CURVE_NONE = 0,     /* No intersection (surface doesn't cross plane) */
    ALEA_CURVE_POINT,        /* Single point (tangent) */
    ALEA_CURVE_LINE,         /* Straight line */
    ALEA_CURVE_LINE_SEGMENT, /* Bounded line segment */
    ALEA_CURVE_RAY,          /* Half-line (one endpoint) */
    ALEA_CURVE_CIRCLE,       /* Full circle */
    ALEA_CURVE_ARC,          /* Circular arc */
    ALEA_CURVE_ELLIPSE,      /* Full ellipse */
    ALEA_CURVE_ELLIPSE_ARC,  /* Elliptical arc */
    ALEA_CURVE_PARABOLA,     /* Parabola (unbounded) */
    ALEA_CURVE_HYPERBOLA,    /* Hyperbola (two branches, unbounded) */
    ALEA_CURVE_POLYGON,      /* Closed polygon (for boxes, wedges) */
    ALEA_CURVE_QUARTIC,      /* Quartic curve (torus intersection - complex) */
    ALEA_CURVE_PARALLEL_LINES, /* Two parallel lines (degenerate conic) */
} alea_curve_type_t;
#endif

/**
 * @brief 2D line: a*u + b*v + c = 0 (implicit form)
 *
 * Parametric form: P(t) = point + t * direction
 */
typedef struct {
    double a, b, c;         /* Implicit: a*u + b*v + c = 0 */
    double point[2];        /* A point on the line */
    double direction[2];    /* Direction vector (unit) */
} alea_line_2d_t;

/**
 * @brief 2D circle: (u - cx)² + (v - cy)² = r²
 */
typedef struct {
    double center[2];       /* Center (cu, cv) */
    double radius;
} alea_circle_2d_t;

/**
 * @brief 2D ellipse in standard form (can be rotated)
 *
 * Parametric: P(θ) = center + cos(θ)*a*axis_u + sin(θ)*b*axis_v
 * where axis_u and axis_v are unit vectors along semi-axes
 */
typedef struct {
    double center[2];       /* Center point */
    double semi_a;          /* Semi-major axis length */
    double semi_b;          /* Semi-minor axis length */
    double angle;           /* Rotation angle (radians, from u-axis) */
} alea_ellipse_2d_t;

/**
 * @brief 2D conic section in general form
 *
 * A*u² + B*u*v + C*v² + D*u + E*v + F = 0
 *
 * Discriminant = B² - 4AC determines type:
 * - < 0: ellipse (or circle if A=C, B=0)
 * - = 0: parabola
 * - > 0: hyperbola
 */
typedef struct {
    double A, B, C, D, E, F;
} alea_conic_2d_t;

/**
 * @brief 2D polygon (for box/wedge intersections)
 */
#define ALEA_MAX_POLYGON_VERTICES 16

typedef struct {
    double vertices[ALEA_MAX_POLYGON_VERTICES][2];
    int vertex_count;
    bool closed;            /* True if polygon is closed */
} alea_polygon_2d_t;

/**
 * @brief Torus quartic curve data
 *
 * Stores parameters needed to evaluate torus-plane intersection.
 * The torus equation in local coords (axis along Z):
 *   (sqrt(x² + y²) - R)² + z² = r²
 *
 * For scanline intersection, we substitute P(u,v) = O + u*U + v*V
 * and solve the resulting quartic in u for fixed v.
 */
typedef struct {
    double R;               /* Major radius */
    double r;               /* Minor radius */
    double Ox, Oy, Oz;      /* Slice origin in torus-local coords */
    double Ux, Uy, Uz;      /* Slice U-axis in torus-local coords */
    double Vx, Vy, Vz;      /* Slice V-axis in torus-local coords */
    int mode;               /* 0=general, 1=concentric circles, 2=separate circles */
    double c1u, c1v, r1;    /* First circle (for mode 1,2) */
    double c2u, c2v, r2;    /* Second circle (for mode 1,2) */
} alea_torus_2d_t;

/**
 * @brief Two parallel lines (degenerate conic)
 *
 * Occurs when slicing a cylinder parallel to its axis.
 * The lines share a common direction and are offset perpendicular to it.
 *
 * For axis-aligned cases (direction along u or v axis):
 * - vertical=true, direction=(0,1): lines at u=offset1 and u=offset2
 * - vertical=false, direction=(1,0): lines at v=offset1 and v=offset2
 *
 * For arbitrary orientation:
 * - Lines pass through point1 and point2 with common direction
 */
typedef struct {
    double point1[2];       /* A point on line 1 */
    double point2[2];       /* A point on line 2 */
    double direction[2];    /* Common direction vector (unit) */
} alea_parallel_lines_2d_t;

/**
 * @brief Union type for 2D curve data
 */
typedef union {
    double point[2];        /* For ALEA_CURVE_POINT */
    alea_line_2d_t line;     /* For ALEA_CURVE_LINE, LINE_SEGMENT, RAY */
    alea_circle_2d_t circle; /* For ALEA_CURVE_CIRCLE, ARC */
    alea_ellipse_2d_t ellipse; /* For ALEA_CURVE_ELLIPSE, ELLIPSE_ARC */
    alea_conic_2d_t conic;   /* For general conics */
    alea_polygon_2d_t polygon; /* For ALEA_CURVE_POLYGON */
    alea_torus_2d_t torus;   /* For ALEA_CURVE_QUARTIC (torus) */
    alea_parallel_lines_2d_t parallel_lines; /* For ALEA_CURVE_PARALLEL_LINES */
} alea_curve_data_t;

/**
 * @brief Bounds for parametric curves (arcs, segments)
 */
typedef struct {
    double t_min, t_max;    /* Parameter bounds */
    double theta_start;     /* For arcs: start angle */
    double theta_end;       /* For arcs: end angle */
} alea_curve_bounds_t;

/**
 * @brief Complete 2D curve description
 */
typedef struct {
    alea_curve_type_t type;
    alea_curve_data_t data;
    alea_curve_bounds_t bounds;
    int surface_id;         /* MCNP surface ID this curve came from */
    int universe_id;        /* Universe ID this curve belongs to */
    int sense;              /* +1 or -1: which side is "inside" */
} alea_curve_2d_t;

/* ============================================================================
 * CURVE INTERSECTION RESULT
 * ============================================================================ */

/**
 * @brief Result of intersecting all surfaces with a slice plane
 */
ALEA_VEC_DEFINE(alea_curve_vec, alea_curve_2d_t);

typedef struct {
    alea_curve_vec_t curves;

    /* Bounding box in plane coordinates */
    double u_min, u_max;
    double v_min, v_max;
} alea_curve_collection_t;

/* ============================================================================
 * CONIC UTILITIES
 * ============================================================================ */

/**
 * @brief Convert general conic to canonical ellipse form
 *
 * Given: A*u² + B*u*v + C*v² + D*u + E*v + F = 0
 * Compute: center, semi-axes, rotation angle
 *
 * @param conic Input conic coefficients
 * @param ellipse Output ellipse parameters
 * @return true if conic is an ellipse, false otherwise
 */
bool alea_conic_to_ellipse(const alea_conic_2d_t* conic, alea_ellipse_2d_t* ellipse);

/* ============================================================================
 * API FUNCTIONS
 * ============================================================================ */

/**
 * @brief Initialize a slice plane from origin, normal, and up vector
 *
 * Computes orthonormal basis vectors U and V from the given normal and up hint.
 *
 * @param plane Output plane structure
 * @param ox, oy, oz Origin point
 * @param nx, ny, nz Normal vector (will be normalized)
 * @param ux, uy, uz Up vector hint (will be orthogonalized)
 */
void alea_slice_plane_init(alea_slice_plane_t* plane,
                          double ox, double oy, double oz,
                          double nx, double ny, double nz,
                          double ux, double uy, double uz);

/**
 * @brief Initialize axis-aligned slice plane
 *
 * @param plane Output plane structure
 * @param axis 0=X (YZ plane), 1=Y (XZ plane), 2=Z (XY plane)
 * @param value Coordinate value along the axis
 */
void alea_slice_plane_init_axis(alea_slice_plane_t* plane, int axis, double value);

/**
 * @brief Transform 3D point to 2D plane coordinates
 *
 * Projects point onto plane and returns (u, v) coordinates.
 * If point is not on plane, projects orthogonally.
 *
 * @param plane Slice plane
 * @param x, y, z 3D point
 * @param u, v Output 2D coordinates
 */
void alea_plane_to_2d(const alea_slice_plane_t* plane,
                     double x, double y, double z,
                     double* u, double* v);

/**
 * @brief Transform 2D plane coordinates to 3D point
 *
 * @param plane Slice plane
 * @param u, v 2D coordinates
 * @param x, y, z Output 3D point
 */
void alea_plane_to_3d(const alea_slice_plane_t* plane,
                     double u, double v,
                     double* x, double* y, double* z);

/**
 * @brief Compute intersection of a single primitive with slice plane
 *
 * @param type Primitive type
 * @param data Primitive data
 * @param plane Slice plane
 * @param curve Output curve (caller provides storage)
 * @return true if intersection exists, false if none
 */
bool alea_intersect_primitive_plane(alea_primitive_type_t type,
                                   const alea_primitive_data_t* data,
                                   const alea_slice_plane_t* plane,
                                   alea_curve_2d_t* curve);

/**
 * @brief Compute all surface-plane intersections for a CSG system
 *
 * Iterates over all surfaces in the system and computes their intersection
 * with the given slice plane.
 *
 * @param sys CSG system
 * @param plane Slice plane
 * @param result Output structure (caller provides, curves will be allocated)
 * @return 0 on success, -1 on error
 */
int alea_compute_slice_curves(const alea_system_t* sys,
                             const alea_slice_plane_t* plane,
                             alea_curve_collection_t* result);

/**
 * @brief Free curves allocated by alea_compute_slice_curves
 */
void alea_curve_collection_free(alea_curve_collection_t* result);

/**
 * @brief Compute slice curves using spatial index (no flattening needed)
 *
 * Uses the spatial index for efficient curve computation:
 * 1. Query spatial index for instances intersecting slice plane
 * 2. For each instance, get cell's surfaces
 * 3. Transform surfaces and compute curve intersections
 * 4. Cull curves outside viewport
 *
 * @param sys CSG system (spatial index will be built if needed)
 * @param plane Slice plane definition
 * @param u_min, u_max, v_min, v_max Viewport in plane coordinates
 * @param result Output curves (caller provides, will be populated)
 * @return 0 on success, -1 on error
 */
int alea_compute_slice_curves_spatial(const alea_system_t* sys,
                                     const alea_slice_plane_t* plane,
                                     double u_min, double u_max,
                                     double v_min, double v_max,
                                     alea_curve_collection_t* result);

/**
 * @brief Enable/disable debug output for slice curve generation
 *
 * When enabled, prints detailed information about:
 * - Which cell instances are being processed
 * - Which surfaces are skipped (deduplication) and why
 * - Which surfaces produce curves and their types
 *
 * @param enable 1 to enable debug output, 0 to disable
 */
void alea_slice_curve_set_debug(int enable);

/**
 * @brief Enable/disable point trace debugging for cell lookup
 * @param enable 1 to enable debug output, 0 to disable
 */
void alea_slice_point_trace_set_debug(int enable);

/* ============================================================================
 * CURVE EVALUATION (for rasterization)
 * ============================================================================ */

/**
 * @brief Evaluate parametric curve at parameter t
 *
 * For lines: P(t) = point + t * direction
 * For circles: P(t) = center + r*(cos(t), sin(t))
 * For ellipses: P(t) = center + a*cos(t)*u + b*sin(t)*v
 *
 * @param curve Curve to evaluate
 * @param t Parameter value
 * @param u, v Output 2D point
 * @return true if valid, false if t is out of bounds
 */
bool alea_curve_eval(const alea_curve_2d_t* curve, double t,
                    double* u, double* v);

/**
 * @brief Find curve intersection with horizontal scanline
 *
 * Returns all u-values where the curve crosses v = v_scanline.
 *
 * @param curve Curve to intersect
 * @param v_scanline V coordinate of scanline
 * @param u_out Output array for u values (caller provides, size >= 4)
 * @param max_intersections Maximum intersections to return
 * @return Number of intersections found
 */
int alea_curve_scanline_intersect(const alea_curve_2d_t* curve,
                                 double v_scanline,
                                 double* u_out,
                                 int max_intersections);

/**
 * @brief Get bounding box of curve in 2D
 *
 * @param curve Curve to bound
 * @param u_min, u_max Output horizontal bounds
 * @param v_min, v_max Output vertical bounds
 */
void alea_curve_bbox(const alea_curve_2d_t* curve,
                    double* u_min, double* u_max,
                    double* v_min, double* v_max);


#endif /* CURVE_INTERSECT_H */
