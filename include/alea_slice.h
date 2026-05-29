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
#include <stddef.h>
#include <stdint.h>

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
    int surface_id;                 /**< MCNP/user surface ID that generated this curve (display) */
    uint32_t primitive_id;          /**< Canonical deduplicated primitive (internal match-only) */
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

/** Coverage class for plot-time geometry diagnostics */
typedef enum {
    ALEA_COVERAGE_NONE  = 0,  /**< No cell claims this point */
    ALEA_COVERAGE_ONE   = 1,  /**< Exactly one known cell claims this point */
    ALEA_COVERAGE_MULTI = 2   /**< Multiple cells claim this point */
} alea_coverage_class_t;

/** Flags for alea_find_cells_grid_coverage() */
#define ALEA_GRID_COVERAGE_NONE_FLAG       0u
#define ALEA_GRID_COVERAGE_FAST            1u
#define ALEA_GRID_COVERAGE_EXACT           2u
#define ALEA_GRID_SECONDARY_CELL_IDS       4u
#define ALEA_GRID_PATH_IDS                 8u

/** Interned concrete universe path used by path-conditioned slice refinement */
typedef struct {
    int universe_id;           /**< Universe at the target/deepest depth */
    int depth;                 /**< Nesting depth */
    uint64_t chain_hash;       /**< Stable hash for diagnostics/grouping */
    double world_to_local[12]; /**< Concrete 3x4 transform values */
} alea_slice_path_record_t;

/** Table of interned slice paths. Caller owns records and frees with helper. */
typedef struct {
    alea_slice_path_record_t* records;
    size_t count;
    size_t capacity;
} alea_slice_path_table_t;

/** Diagnostic counters from the most recent tile coverage refinement pass */
typedef struct {
    size_t tiles;                    /**< Number of tiles visited */
    size_t fallback_tiles;           /**< Tiles refined with exact per-pixel fallback */
    size_t query_errors;             /**< Spatial candidate query failures */
    size_t pixels;                   /**< Pixels evaluated by tile candidates */
    size_t exact_fallback_pixels;    /**< Pixels evaluated by exact fallback */
    size_t candidate_total;          /**< Raw spatial candidates across tiles */
    size_t candidate_max;            /**< Maximum raw candidates in one tile */
    size_t dedup_candidate_total;    /**< Deduplicated candidates across tiles */
    size_t dedup_candidate_max;      /**< Maximum deduplicated candidates in one tile */
    size_t candidate_pixel_tests;    /**< Deduplicated candidate/pixel containment tests */
    size_t refined_pixels;           /**< Pixels newly classified as overlap */
    size_t path_groups;              /**< Path-local groups evaluated */
    size_t path_group_pixels_max;    /**< Maximum pixels in one path-local group */
    size_t path_group_candidates_max; /**< Maximum candidates in one path group */
    size_t bbox_pixel_tests;         /**< Candidate bbox-to-pixel tests */
    size_t bbox_pixel_rejects;       /**< Candidate/pixel pairs rejected by bbox */
    size_t early_multi_skips;        /**< Candidate/pixel pairs skipped after MULTI */
    size_t primary_cell_skips;       /**< Candidate/pixel pairs skipped as primary */
    size_t path_2d_verify_queries;   /**< 2D buckets checked against 3D query */
    size_t path_2d_missing_candidates; /**< 3D candidates absent from 2D buckets */
    size_t path_2d_missing_tiles;     /**< Path/tile groups with missing candidates */
} alea_tile_coverage_stats_t;

/** Diagnostic counters from boundary-ambiguity filtering */
typedef struct {
    size_t checked;      /**< Provisional MULTI pixels tested */
    size_t suppressed;   /**< Boundary-only overlaps downgraded */
    size_t retained;     /**< Overlaps confirmed by off-boundary samples */
    size_t inconclusive; /**< Conservative keeps due to ambiguous samples */
} alea_boundary_filter_stats_t;

/** Diagnostic counters from exact point-coverage refinement */
typedef struct {
    size_t queries;                  /**< Total exact coverage point queries */
    size_t spatial_queries;          /**< Queries handled by spatial candidates */
    size_t spatial_multi_early_exit; /**< Spatial queries stopped at second claimant */
    size_t recursive_fallbacks;      /**< Queries handled by recursive all-hit fallback */
    size_t lattice_fallbacks;        /**< Fallbacks caused by lattice semantics */
    size_t truncated_fallbacks;      /**< Fallbacks caused by candidate truncation */
    size_t query_errors;             /**< Point coverage query errors */
    size_t candidate_total;          /**< Deduplicated spatial candidates scanned */
    size_t candidate_max;            /**< Max deduplicated spatial candidates */
    size_t contains_tests;           /**< Exact containment tests in spatial path */
} alea_point_coverage_stats_t;

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
 * @brief Find cells and coverage classes on a slice grid
 *
 * Compatibility extension of alea_find_cells_grid(). The cell/material/error
 * arrays have the same meaning as alea_find_cells_grid(). out_coverage, when
 * provided, receives ALEA_COVERAGE_NONE/ONE/MULTI for each pixel.
 *
 * FAST derives coverage from the normal grid result: undefined pixels become
 * NONE, overlap pixels become MULTI, and remaining valid pixels become ONE.
 * EXACT rechecks pixels with all-cell semantics, which detects total/nested
 * overlaps that have no winning-cell boundary.
 *
 * @param flags ALEA_GRID_COVERAGE_* flags
 * @param out_secondary_cell_ids Optional secondary cell IDs for MULTI pixels
 * @param out_coverage Optional coverage class array (size nu*nv)
 */
int alea_find_cells_grid_coverage(alea_system_t* sys,
                                  const alea_slice_view_t* view,
                                  int nu, int nv,
                                  int universe_depth,
                                  unsigned flags,
                                  int* out_cell_ids,
                                  int* out_material_ids,
                                  int* out_secondary_cell_ids,
                                  uint8_t* out_coverage,
                                  uint8_t* out_errors);

/**
 * @brief Find cells, coverage, and optional concrete path IDs on a slice grid
 *
 * Extends alea_find_cells_grid_coverage() with an interned path table. This is
 * currently intended for path-conditioned large-model error refinement. Path IDs
 * are filled when ALEA_GRID_PATH_IDS is set and out_path_ids/out_paths are
 * provided. Pixels without a stable path receive UINT32_MAX.
 */
int alea_find_cells_grid_coverage_paths(alea_system_t* sys,
                                        const alea_slice_view_t* view,
                                        int nu, int nv,
                                        int universe_depth,
                                        unsigned flags,
                                        int* out_cell_ids,
                                        int* out_material_ids,
                                        int* out_secondary_cell_ids,
                                        uint8_t* out_coverage,
                                        uint8_t* out_errors,
                                        uint32_t* out_path_ids,
                                        alea_slice_path_table_t* out_paths);

/** Free records owned by a slice path table */
void alea_slice_path_table_free(alea_slice_path_table_t* table);

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

/**
 * @brief Mark curve-adjacent overlap pixels from precomputed coverage
 *
 * Coverage-driven companion to alea_check_grid_overlaps_curves(). It performs
 * the same curve sampling, but never queries geometry. Pixels sampled near a
 * curve are marked as overlap when coverage[pixel] is ALEA_COVERAGE_MULTI.
 *
 * @return Number of newly marked overlap pixels, or -1 on error
 */
int alea_check_grid_overlaps_curves_coverage(
    const alea_slice_view_t* view,
    const alea_slice_curves_t* curves,
    int nu, int nv,
    const uint8_t* coverage,
    uint8_t* errors);

/**
 * @brief Refine coverage by exact queries near surface curves
 *
 * Samples surface curves, maps samples to nearby pixels, and rechecks only
 * those pixels with all-cell semantics. This is a selective exact pass for
 * total/nested overlaps near hidden boundaries, cheaper than exact coverage
 * over the full grid.
 *
 * @param out_secondary_cell_ids Optional secondary cell ID array, can be NULL
 * @param coverage Coverage array to update (size nu*nv)
 * @param errors Error array to update (size nu*nv)
 * @return Number of pixels newly or still classified as multi-coverage, or -1
 */
int alea_refine_grid_coverage_curves_exact(
    alea_system_t* sys,
    const alea_slice_view_t* view,
    const alea_slice_curves_t* curves,
    int nu, int nv,
    int universe_depth,
    int* out_secondary_cell_ids,
    uint8_t* coverage,
    uint8_t* errors);

/**
 * @brief Refine coverage by querying spatial candidates per tile
 *
 * Experimental tile-based exact coverage pass. For each tile, it queries the
 * spatial index once for candidate cells, evaluates those candidates over the
 * tile pixels, and updates coverage/errors. If a tile candidate list is
 * truncated, that tile falls back to exact per-pixel queries.
 *
 * @param tile_w,tile_h Tile dimensions; use positive values such as 16x16
 * @return Number of pixels newly classified as multi-coverage, or -1 on error
 */
int alea_refine_grid_coverage_tiles_exact(
    alea_system_t* sys,
    const alea_slice_view_t* view,
    int nu, int nv,
    int universe_depth,
    int tile_w, int tile_h,
    int* out_secondary_cell_ids,
    uint8_t* coverage,
    uint8_t* errors);

/**
 * @brief Refine coverage by querying spatial candidates per concrete path
 *
 * Experimental path-conditioned exact coverage pass. For each tile, pixels are
 * grouped by their concrete universe path, queried in that universe's local
 * coordinates, and evaluated against only sibling candidates. Pixels without a
 * stable path fall back to exact per-pixel queries.
 *
 * @param path_ids Path ID grid produced by alea_find_cells_grid_coverage_paths()
 * @param paths Interned path table produced with ALEA_GRID_PATH_IDS
 * @return Number of pixels newly classified as multi-coverage, or -1 on error
 */
int alea_refine_grid_coverage_paths_exact(
    alea_system_t* sys,
    const alea_slice_view_t* view,
    int nu, int nv,
    int universe_depth,
    int tile_w, int tile_h,
    const int* primary_cell_ids,
    const uint32_t* path_ids,
    const alea_slice_path_table_t* paths,
    int* out_secondary_cell_ids,
    uint8_t* coverage,
    uint8_t* errors);

/**
 * @brief Suppress exact-boundary overlap ambiguities in plot coverage grids
 *
 * Tests provisional MULTI pixels by sampling small offsets in the slice plane.
 * Boundary-only multi-hits are downgraded to ONE/no-error. Finite overlaps and
 * inconclusive cases are retained.
 *
 * @return Number of suppressed pixels, or -1 on invalid input
 */
int alea_filter_grid_boundary_ambiguities(
    alea_system_t* sys,
    const alea_slice_view_t* view,
    int nu, int nv,
    int universe_depth,
    const int* primary_cell_ids,
    int* secondary_cell_ids,
    uint8_t* coverage,
    uint8_t* errors,
    alea_boundary_filter_stats_t* out_stats);

/** Reset diagnostics for tile coverage refinement */
void alea_tile_coverage_stats_reset(void);

/** Return diagnostics from the most recent tile coverage refinement pass */
alea_tile_coverage_stats_t alea_tile_coverage_stats_get(void);

/** Reset diagnostics for exact point-coverage refinement */
void alea_point_coverage_stats_reset(void);

/** Return diagnostics from exact point-coverage refinement */
alea_point_coverage_stats_t alea_point_coverage_stats_get(void);

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
 * Works in both flat and hierarchical spatial modes.
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
 * @brief Evaluate a curve at parameter t, returning the 2D (u,v) point.
 *
 * Thin accessor over the internal curve evaluator so external consumers
 * (e.g. the geometry validator) can sample boundary curves without depending
 * on slice internals.
 *
 * @return 0 on success, -1 on error or if the point is undefined at t.
 */
int alea_slice_curve_eval(const alea_slice_curves_t* curves, size_t index,
                          double t, double* u, double* v);

/**
 * @brief Get the parameter range of a curve clipped to the view's viewport.
 *
 * Mirrors the clipping used by the slice error sampler. On error both outputs
 * are set to 0.
 */
void alea_slice_curve_param_range(const alea_slice_curves_t* curves, size_t index,
                                  const alea_slice_view_t* view,
                                  double* t_lo, double* t_hi);

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

/** Plot-level geometry error component kind */
typedef enum {
    ALEA_PLOT_ERR_UNDEFINED_REGION = 1, /**< Connected undefined/gap pixels */
    ALEA_PLOT_ERR_PARTIAL_OVERLAP  = 2, /**< Overlap with visible non-overlap for participants */
    ALEA_PLOT_ERR_TOTAL_OVERLAP    = 3  /**< Nested/coincident overlap with hidden participant */
} alea_plot_error_kind_t;

/** Connected component of plot-level error pixels */
typedef struct {
    alea_plot_error_kind_t kind;
    int primary_cell_id;
    int secondary_cell_id;
    int pixel_count;
    int min_i, min_j, max_i, max_j;
} alea_plot_error_component_t;

/** Result of plot-level component classification */
typedef struct {
    alea_plot_error_component_t* components;
    size_t component_count;
} alea_plot_error_component_result_t;

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
 * @deprecated Use alea_validate_geometry_slice() (alea_geo_validator.h), which
 * samples the same boundary curves but returns structured, primitive-keyed
 * events with five diagnostic types instead of OVERLAP/GAP segments, and avoids
 * deduplicated-surface false positives.
 *
 * @param sys CSG system
 * @param view Slice view
 * @param curves Previously computed curves from alea_get_slice_curves()
 * @param universe_depth Which level of universe hierarchy (-1 for innermost)
 * @return Error result (must be freed with alea_slice_errors_free), or NULL
 */
ALEA_DEPRECATED("use alea_validate_geometry_slice")
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
 *
 * @deprecated Use alea_validate_geometry_slice() (alea_geo_validator.h) for
 * structured, primitive-keyed boundary diagnostics. The grid per-pixel error
 * flags from alea_find_cells_grid() remain available for fast pixel overlays.
 */
ALEA_DEPRECATED("use alea_validate_geometry_slice")
alea_slice_error_result_t* alea_check_slice_errors_grid(
    alea_system_t* sys,
    const alea_slice_view_t* view,
    const alea_slice_curves_t* curves,
    const int* cell_ids,
    const uint8_t* grid_errors,
    int nu, int nv);

/**
 * @brief Check surface curves using pre-computed grid coverage
 *
 * Variant of alea_check_slice_errors_grid() that does not perform same-cell CSG
 * fallback queries. When both sides of a curve map to the same winning cell,
 * the supplied coverage/error grids decide whether that curve sample is an
 * overlap, gap, or clean region.
 *
 * @deprecated Use alea_validate_geometry_slice() (alea_geo_validator.h) for
 * structured, primitive-keyed boundary diagnostics.
 */
ALEA_DEPRECATED("use alea_validate_geometry_slice")
alea_slice_error_result_t* alea_check_slice_errors_grid_ex(
    const alea_slice_view_t* view,
    const alea_slice_curves_t* curves,
    const int* cell_ids,
    const uint8_t* coverage,
    const uint8_t* grid_errors,
    int nu, int nv);

/** Free error result */
void alea_slice_errors_free(alea_slice_error_result_t* result);

/**
 * @brief Classify connected error components from plot coverage arrays
 *
 * Undefined components are built from ALEA_COVERAGE_NONE. Overlap components are
 * built from ALEA_COVERAGE_MULTI and classified as partial when a participant
 * is visible in adjacent non-overlap pixels; otherwise they are classified as
 * total/nested overlap.
 *
 * @param cell_ids Primary cell IDs for the same grid
 * @param secondary_cell_ids Optional secondary IDs from coverage refinement
 * @param coverage Coverage class grid
 * @return Component result (free with alea_plot_error_components_free), or NULL
 */
alea_plot_error_component_result_t* alea_classify_plot_error_components(
    const int* cell_ids,
    const int* secondary_cell_ids,
    const uint8_t* coverage,
    int nu, int nv);

/** Free plot error component result */
void alea_plot_error_components_free(alea_plot_error_component_result_t* result);

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
