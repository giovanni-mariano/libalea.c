// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/**
 * @file alea_slice.h
 * @brief Alea Slice/Render Module API
 *
 * 2D slice functionality for CSG geometry visualization.
 * Requires linking with libalea_slice.a
 *
 * ## Recommended Usage
 *
 * The recommended approach uses two main operations:
 *
 * 1. **Grid-based cell queries** (`alea_find_cells_grid()`)
 *    - Returns cell/material IDs for each pixel in a grid
 *    - Supports universe hierarchy traversal and error detection
 *    - Very fast due to spatial coherence optimization
 *
 * 2. **Analytical curve extraction** (`alea_get_slice_curves()`)
 *    - Returns exact surface boundaries as curves (lines, circles, ellipses)
 *    - For rendering surface contours with subpixel accuracy
 *
 * Example workflow (see examples/plotter.c):
 * @code
 *   // 1. Set up a slice view
 *   alea_slice_view_t view;
 *   alea_slice_view_axis(&view, 2, z, x_min, x_max, y_min, y_max);
 *
 *   // 2. Get cell/material IDs for coloring
 *   int* cell_ids = malloc(width * height * sizeof(int));
 *   alea_find_cells_grid(sys, &view, width, height, -1,
 *                            cell_ids, NULL, NULL);
 *
 *   // 3. Get curves for drawing boundaries
 *   alea_slice_curves_t* curves = alea_get_slice_curves(sys, &view);
 *
 *   // 4. Optionally get label positions
 *   alea_label_position_t* labels;
 *   int label_count;
 *   alea_find_label_positions(cell_ids, width, height, 100, &labels, &label_count);
 * @endcode
 */

#ifndef ALEA_SLICE_H
#define ALEA_SLICE_H

#include "alea.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ============================================================================
 * CORE TYPES
 * ============================================================================ */

/** Slice plane definition for arbitrary orientation */
typedef struct {
    double origin[3];   /**< Point on the plane */
    double normal[3];   /**< Unit normal vector */
    double u_axis[3];   /**< Unit vector for horizontal (U) direction in plane */
    double v_axis[3];   /**< Unit vector for vertical (V) direction in plane */
} alea_slice_plane_t;

/** Slice view: plane + viewport bounds */
typedef struct {
    alea_slice_plane_t plane;   /**< Slice plane (orientation + position) */
    double u_min, u_max;            /**< Horizontal viewport bounds */
    double v_min, v_max;            /**< Vertical viewport bounds */
} alea_slice_view_t;

/** Opaque slice curves result */
typedef struct alea_slice_curves alea_slice_curves_t;

/** Curve types from surface intersections */
#ifndef ALEA_CURVE_TYPE_DEFINED
#define ALEA_CURVE_TYPE_DEFINED
typedef enum {
    ALEA_CURVE_NONE = 0,
    ALEA_CURVE_POINT,           /**< Single point (tangent) */
    ALEA_CURVE_LINE,            /**< Infinite line */
    ALEA_CURVE_LINE_SEGMENT,    /**< Finite line segment */
    ALEA_CURVE_RAY,             /**< Half-line (one endpoint) */
    ALEA_CURVE_CIRCLE,          /**< Full circle */
    ALEA_CURVE_ARC,             /**< Circular arc */
    ALEA_CURVE_ELLIPSE,         /**< Full ellipse */
    ALEA_CURVE_ELLIPSE_ARC,     /**< Elliptical arc */
    ALEA_CURVE_PARABOLA,        /**< Parabola (unbounded) */
    ALEA_CURVE_HYPERBOLA,       /**< Hyperbola (two branches, unbounded) */
    ALEA_CURVE_POLYGON,         /**< Polygon (from polyhedron faces) */
    ALEA_CURVE_QUARTIC,         /**< Quartic curve (torus intersection) */
    ALEA_CURVE_PARALLEL_LINES   /**< Two parallel lines (degenerate conic) */
} alea_curve_type_t;
#endif

/** 2D curve from slice intersection */
typedef struct {
    alea_curve_type_t type;
    int surface_id;                 /**< Surface ID that generated this curve */
    union {
        struct { double point[2]; double direction[2]; } line;
        struct { double center[2]; double radius; } circle;
        struct { double center[2]; double semi_a, semi_b, angle; } ellipse;
        struct { double vertices[16][2]; int count; int closed; } polygon;
        struct { double point1[2]; double point2[2]; double direction[2]; } parallel_lines;
    } data;
    double t_min, t_max;            /**< Parameter range (0 to 2*PI for full circles) */
} alea_curve_t;

/** Error codes for grid cell queries */
typedef enum {
    ALEA_GRID_OK        = 0,  /**< No error - valid cell */
    ALEA_GRID_VOID      = 0,  /**< Void region (cell_id = -1) */
    ALEA_GRID_OVERLAP   = 1,  /**< Multiple cells claim this point */
    ALEA_GRID_UNDEFINED = 2   /**< No cell claims this point (geometry error) */
} alea_grid_error_t;

/** Label position information */
typedef struct {
    int id;             /**< Cell/material/surface ID */
    int px, py;         /**< Pixel coordinates for label placement */
    int pixel_count;    /**< Region size (for cells/materials) or 0 (for surfaces) */
} alea_label_position_t;

/* ============================================================================
 * SLICE VIEW SETUP
 * ============================================================================ */

/**
 * @brief Initialize a slice view for an axis-aligned slice
 * @param view Output view structure
 * @param axis 0=X (YZ plane), 1=Y (XZ plane), 2=Z (XY plane)
 * @param value Coordinate value along the axis
 * @param u_min, u_max Horizontal viewport bounds
 * @param v_min, v_max Vertical viewport bounds
 */
void alea_slice_view_axis(alea_slice_view_t* view,
                               int axis, double value,
                               double u_min, double u_max,
                               double v_min, double v_max);

/**
 * @brief Initialize a slice view with arbitrary orientation
 * @param view Output view structure
 * @param ox, oy, oz Origin point
 * @param nx, ny, nz Normal vector (will be normalized)
 * @param ux, uy, uz Up vector hint (will be orthogonalized)
 * @param u_min, u_max Horizontal viewport bounds
 * @param v_min, v_max Vertical viewport bounds
 */
void alea_slice_view_init(alea_slice_view_t* view,
                               double ox, double oy, double oz,
                               double nx, double ny, double nz,
                               double ux, double uy, double uz,
                               double u_min, double u_max,
                               double v_min, double v_max);

/* ============================================================================
 * GRID-BASED CELL QUERIES
 * ============================================================================ */

/**
 * @brief Find cells on a 2D grid defined by a slice view
 *
 * Supports:
 * - Universe hierarchy traversal (universe_depth parameter)
 * - Geometry error detection (overlaps, undefined regions)
 *
 * For universe_depth:
 * - -1 = innermost cell (follows all fills)
 * -  0 = root level only (fastest, uses adjacency optimization)
 * -  N = cell at depth N
 *
 * @param sys CSG system
 * @param view Slice view (plane + viewport bounds)
 * @param nu Number of pixels in horizontal direction
 * @param nv Number of pixels in vertical direction
 * @param universe_depth Which level of universe hierarchy to return (-1 for innermost)
 * @param out_cell_ids Output cell IDs (size nu*nv, -1 for void/error)
 * @param out_material_ids Output material IDs (size nu*nv, can be NULL)
 * @param out_errors Output error codes (size nu*nv, can be NULL)
 *                   0=ok, 1=overlap, 2=undefined
 */
int alea_find_cells_grid(alea_system_t* sys,
                              const alea_slice_view_t* view,
                              int nu, int nv,
                              int universe_depth,
                              int* out_cell_ids,
                              int* out_material_ids,
                              uint8_t* out_errors);

/**
 * @brief Check grid for overlapping cells (comprehensive)
 *
 * Re-queries every non-void pixel with the full hierarchy search to detect
 * all overlaps, including fully-nested ones (e.g., concentric spheres).
 * Call after alea_find_cells_grid() for comprehensive overlap detection.
 *
 * Note: alea_find_cells_grid() already detects overlaps at cell boundaries
 * when out_errors is provided. This function is for cases where you need to
 * detect all overlaps including fully-nested geometry errors.
 *
 * @param sys CSG system
 * @param view Slice view (plane + viewport bounds)
 * @param nu Number of pixels in horizontal direction
 * @param nv Number of pixels in vertical direction
 * @param universe_depth Which level of universe hierarchy (-1 for innermost)
 * @param cell_ids Cell ID grid from alea_find_cells_grid()
 * @param errors Error grid to update (size nu*nv, must be pre-initialized)
 * @return 0 on success, -1 on error
 */
int alea_check_grid_overlaps(alea_system_t* sys,
                                  const alea_slice_view_t* view,
                                  int nu, int nv,
                                  int universe_depth,
                                  const int* cell_ids,
                                  uint8_t* errors);

/**
 * @brief Detect nested overlaps by probing pixels along surface curves
 *
 * Efficient alternative to alea_check_grid_overlaps(): instead of
 * re-querying every pixel, only probes pixels where a surface curve
 * crosses the grid without a cell-ID transition (i.e., same cell on
 * both sides). These are the only places where a fully-nested overlap
 * (e.g., sphere inside sphere) can hide.
 *
 * Call after alea_find_cells_grid() and alea_get_slice_curves().
 *
 * @param sys CSG system
 * @param view Slice view (plane + viewport bounds)
 * @param curves Previously computed curves for this view
 * @param nu Grid width
 * @param nv Grid height
 * @param universe_depth Which level of universe hierarchy (-1 for innermost)
 * @param cell_ids Cell ID grid from alea_find_cells_grid()
 * @param errors Error grid to update (size nu*nv)
 * @return Number of new overlap pixels found, or -1 on error
 */
int alea_check_grid_overlaps_curves(alea_system_t* sys,
                                    const alea_slice_view_t* view,
                                    const alea_slice_curves_t* curves,
                                    int nu, int nv,
                                    int universe_depth,
                                    const int* cell_ids,
                                    uint8_t* errors);

/* ============================================================================
 * LABEL POSITION COMPUTATION
 * ============================================================================ */

/**
 * @brief Find optimal label positions for regions in a cell/material grid
 *
 * For each unique region, finds a point that is:
 * 1. Guaranteed to be inside the region (not just centroid)
 * 2. Close to the region's geometric center
 *
 * Handles non-convex regions correctly (C-shapes, donuts, etc.)
 * Uses reservoir sampling + closest-to-centroid selection.
 *
 * @param ids          ID grid (cell_ids or material_ids from alea_find_cells_grid)
 * @param width        Grid width in pixels
 * @param height       Grid height in pixels
 * @param min_pixels   Minimum region size to include (e.g., 100 to skip tiny regions)
 * @param out_labels   Output array of label positions (caller frees with free())
 * @param out_count    Number of labels returned
 * @return 0 on success, -1 on error
 */
int alea_find_label_positions(
    const int* ids,
    int width, int height,
    int min_pixels,
    alea_label_position_t** out_labels,
    int* out_count
);

/**
 * @brief Find optimal label positions for surfaces along curve boundaries
 *
 * For each unique surface in the curves, finds a point along the curve
 * suitable for label placement.
 *
 * @param curves       Slice curves from alea_get_slice_curves()
 * @param x_min, x_max World coordinate bounds (for pixel conversion)
 * @param y_min, y_max World coordinate bounds (for pixel conversion)
 * @param width        Image width in pixels
 * @param height       Image height in pixels
 * @param margin       Margin from image edge (labels within margin are skipped)
 * @param out_labels   Output array of label positions (caller frees with free())
 * @param out_count    Number of labels returned
 * @return 0 on success, -1 on error
 */
int alea_find_surface_label_positions(
    const alea_slice_curves_t* curves,
    double x_min, double x_max,
    double y_min, double y_max,
    int width, int height,
    int margin,
    alea_label_position_t** out_labels,
    int* out_count
);

/**
 * @brief Find surface label positions only where contours are drawn
 *
 * Like alea_find_surface_label_positions(), but filters candidate positions
 * against the same cell/material ID grid used for contour drawing. Pass
 * cell_ids for cell contours or material_ids for material contours.
 *
 * @param curves       Slice curves from alea_get_slice_curves()
 * @param boundary_ids Cell/material ID grid in math order (width*height)
 * @param x_min, x_max World coordinate bounds (for pixel conversion)
 * @param y_min, y_max World coordinate bounds (for pixel conversion)
 * @param width        Image width in pixels
 * @param height       Image height in pixels
 * @param margin       Margin from image edge (labels within margin are skipped)
 * @param out_labels   Output array of label positions (caller frees with free())
 * @param out_count    Number of labels returned
 * @return 0 on success, -1 on error
 */
int alea_find_surface_label_positions_on_boundaries(
    const alea_slice_curves_t* curves,
    const int* boundary_ids,
    double x_min, double x_max,
    double y_min, double y_max,
    int width, int height,
    int margin,
    alea_label_position_t** out_labels,
    int* out_count
);

/* ============================================================================
 * ANALYTICAL CURVE API
 * ============================================================================ */

/**
 * @brief Get analytical curves from a slice view
 *
 * Returns exact surface boundaries as parametric curves for accurate rendering.
 *
 * @param sys CSG system
 * @param view Slice view (plane + viewport bounds)
 * @return Opaque curves object (must be freed with alea_slice_curves_free)
 */
alea_slice_curves_t* alea_get_slice_curves(alea_system_t* sys,
                                                    const alea_slice_view_t* view);

/**
 * @brief Get number of curves in result
 */
size_t alea_slice_curves_count(const alea_slice_curves_t* curves);

/**
 * @brief Get curve data at index
 * @param curves Curves result
 * @param index Curve index (0 to count-1)
 * @param out Output curve structure
 * @return 0 on success, -1 on error
 */
int alea_slice_curves_get(const alea_slice_curves_t* curves, size_t index, alea_curve_t* out);

/**
 * @brief Get bounding box of all curves
 */
void alea_slice_curves_bounds(const alea_slice_curves_t* curves,
                                  double* u_min, double* u_max,
                                  double* v_min, double* v_max);

/**
 * @brief Free curves result
 */
void alea_slice_curves_free(alea_slice_curves_t* curves);

/* ============================================================================
 * ANALYTICAL ERROR LINE CHECKING
 * ============================================================================ */

/** Error type for analytical error lines */
typedef enum {
    ALEA_SLICE_ERR_OVERLAP = 1,  /**< Multiple cells on one/both sides */
    ALEA_SLICE_ERR_GAP     = 2   /**< No cell on one/both sides */
} alea_slice_error_type_t;

/** A segment of a curve identified as an error */
typedef struct {
    size_t curve_index;            /**< Index into curves collection */
    int surface_id;                /**< Surface ID for reference */
    alea_slice_error_type_t type;  /**< OVERLAP or GAP */
    double t_start, t_end;         /**< Parameter range of error segment */
} alea_slice_error_t;

/** Result of error line checking */
typedef struct {
    alea_slice_error_t* errors;
    size_t error_count;
} alea_slice_error_result_t;

/**
 * @brief Check surface curves for geometry errors (overlaps/gaps)
 *
 * Samples points along each curve, offsets to both sides, and checks
 * if exactly one cell exists on each side. Returns segments where
 * errors are detected.
 *
 * Uses CSG point-containment queries. Prefer alea_check_slice_errors_grid()
 * when a cell grid is already available.
 *
 * @param sys CSG system
 * @param view Slice view
 * @param curves Previously computed curves from alea_get_slice_curves()
 * @param universe_depth Which level of universe hierarchy (-1 for innermost)
 * @return Error result (must be freed with alea_slice_errors_free), or NULL
 */
alea_slice_error_result_t* alea_check_slice_errors(
    alea_system_t* sys,
    const alea_slice_view_t* view,
    const alea_slice_curves_t* curves,
    int universe_depth);

/**
 * @brief Check surface curves for geometry errors using a pre-computed grid
 *
 * Fast version that uses O(1) pixel lookups in the cell grid for the common
 * case. Falls back to a single CSG query only when both sides of a curve
 * show the same cell (possible nested overlap, e.g., sphere inside sphere).
 *
 * The grid must come from alea_find_cells_grid() for the same view.
 *
 * Classification per sample:
 * - Both sides have different cells → OK (real boundary)
 * - One/both sides void (cell_id < 0) → GAP
 * - Grid pixel flagged as overlap → OVERLAP
 * - Same cell both sides → CSG fallback to check for nested overlap
 *
 * @param sys CSG system (for nested overlap fallback, can be NULL to skip)
 * @param view Slice view (must match grid)
 * @param curves Previously computed curves
 * @param cell_ids Cell ID grid from alea_find_cells_grid() (size nu*nv)
 * @param grid_errors Error grid from alea_find_cells_grid() (size nu*nv, can be NULL)
 * @param nu Grid width
 * @param nv Grid height
 * @return Error result (must be freed with alea_slice_errors_free), or NULL
 */
alea_slice_error_result_t* alea_check_slice_errors_grid(
    alea_system_t* sys,
    const alea_slice_view_t* view,
    const alea_slice_curves_t* curves,
    const int* cell_ids,
    const uint8_t* grid_errors,
    int nu, int nv);

/** Free error result */
void alea_slice_errors_free(alea_slice_error_result_t* result);

/* ============================================================================
 * DEBUG UTILITIES
 * ============================================================================ */

/**
 * @brief Enable/disable debug output for slice curve generation
 *
 * When enabled, prints detailed information about:
 * - Which cell instances are being processed
 * - Which surfaces are skipped (deduplication) and why
 * - Which surfaces produce curves and their types
 * - For torus curves, the mode and circle parameters
 *
 * Useful for debugging issues where surfaces or curve components
 * are unexpectedly missing from the slice output.
 *
 * @param enable 1 to enable debug output to stdout, 0 to disable
 */
void alea_slice_curve_set_debug(int enable);

/**
 * @brief Enable/disable point trace debugging for cell lookup
 *
 * When enabled, prints detailed trace information for each point lookup,
 * showing which cells are tested and why they pass/fail containment tests.
 * Useful for debugging issues where cells are not found at expected locations.
 *
 * @param enable 1 to enable debug output to stdout, 0 to disable
 */
void alea_slice_point_trace_set_debug(int enable);

#ifdef __cplusplus
}
#endif

#endif /* ALEA_SLICE_H */
