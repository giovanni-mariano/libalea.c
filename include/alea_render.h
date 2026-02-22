// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/**
 * @file alea_render.h
 * @brief 3D batch renderer for Alea CSG geometry
 *
 * Produces publication-quality 3D images of CSG models with cutaway views,
 * Phong shading, shadow rays, and material coloring.
 *
 * Uses the existing raycast pipeline (BVH + cell-aware traversal + lattice DDA).
 * Pure CPU, OpenMP parallelized via tile-based rendering.
 */

#ifndef ALEA_RENDER_H
#define ALEA_RENDER_H

#include "alea.h"
#include "alea_raycast.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ============================================================================
 * CONSTANTS
 * ============================================================================ */

#define RENDER_MAX_CLIPS       16
#define RENDER_MAX_AA          8
#define RENDER_DEFAULT_WIDTH   1920
#define RENDER_DEFAULT_HEIGHT  1080
#define RENDER_DEFAULT_FOV     45.0
#define RENDER_DEFAULT_TILE    32
#define RENDER_PALETTE_SIZE    32

/* ============================================================================
 * ENUMS
 * ============================================================================ */

/** Coloring mode */
typedef enum {
    RENDER_COLOR_MATERIAL = 0,  /**< Color by material ID (default) */
    RENDER_COLOR_CELL,          /**< Color by cell ID */
    RENDER_COLOR_UNIVERSE,      /**< Color by universe ID */
    RENDER_COLOR_DENSITY,       /**< Continuous density gradient */
} render_color_mode_t;

/** Rendering mode */
typedef enum {
    RENDER_MODE_SOLID = 0,      /**< Solid + cutaway (default) */
    RENDER_MODE_XRAY,           /**< Semi-transparent X-ray / radiograph */
    RENDER_MODE_DEPTH,          /**< Depth map (grayscale) */
    RENDER_MODE_CELLID,         /**< Cell ID map */
    RENDER_MODE_MATID,          /**< Material ID map */
} render_mode_t;

/* ============================================================================
 * DATA STRUCTURES
 * ============================================================================ */

/** Camera definition */
typedef struct {
    double eye[3];          /**< Camera position */
    double target[3];       /**< Look-at point */
    double up[3];           /**< Up vector */
    double fov_degrees;     /**< Vertical FOV (0 = orthographic) */
    double ortho_height;    /**< Ortho view height (if fov==0) */

    /* Derived (computed by render_camera_setup) */
    double forward[3];
    double right[3];
    double true_up[3];
    double half_height;
    double half_width;
    int auto_fit;           /**< 1 if eye/target were auto-computed */
} render_camera_t;

/** Clipping plane: visible where n.x + d > 0 */
typedef struct {
    double normal[3];
    double d;
} render_clip_plane_t;

/** Custom color entry (from palette file) */
typedef struct {
    int id;
    float r, g, b;
} render_color_entry_t;

/** Render configuration (parsed from CLI) */
typedef struct {
    /* Image */
    int width;
    int height;
    float background[3];

    /* Camera */
    double eye[3];
    double target[3];
    double up[3];
    double fov;
    double ortho_height;
    int eye_set;
    int target_set;

    /* Clipping */
    render_clip_plane_t clips[RENDER_MAX_CLIPS];
    int num_clips;

    /* Appearance */
    render_color_mode_t color_mode;
    render_mode_t render_mode;
    double light_dir[3];
    int light_set;
    int shadows;
    int edges;
    float ambient;
    float diffuse;
    float specular;
    float shininess;
    float cross_section_tint;

    /* Quality */
    int aa_samples;         /**< NxN supersampling (1=off) */
    int aa_adaptive;        /**< Adaptive AA */
    int universe_depth;     /**< -1=innermost, 0=root, N=level */

    /* Performance */
    int threads;
    int tile_size;
    int preview;

    /* Auxiliary output */
    int aux_output;

    /* Custom colors */
    render_color_entry_t* custom_colors;
    int num_custom_colors;

    /* X-ray mode */
    float xray_density_scale;

    /* Surface dedup */
    int dedup;

    /* Logging */
    int log_level;
} render_config_t;

/** Per-pixel hit result (internal) */
typedef struct {
    float color[3];         /**< Final shaded color */
    int cell_id;            /**< Cell ID at this pixel (-1 for background) */
    int material_id;        /**< Material ID */
    float depth;            /**< Distance from camera */
    float normal[3];        /**< World-space normal */
} render_pixel_t;

/** Framebuffer */
typedef struct {
    int width;
    int height;
    float* color;           /**< RGB float, width*height*3 */
    int* cell_id;           /**< Cell ID per pixel, width*height */
    int* material_id;       /**< Material ID per pixel, width*height */
    float* depth;           /**< Depth per pixel, width*height */
    float* normal;          /**< Normal per pixel, width*height*3 */
} render_framebuffer_t;

/* ============================================================================
 * API FUNCTIONS
 * ============================================================================ */

/**
 * @brief Initialize render config with defaults
 */
void render_config_init(render_config_t* cfg);

/**
 * @brief Free render config resources (custom color table)
 */
void render_config_free(render_config_t* cfg);

/**
 * @brief Set up camera from config (computes derived quantities)
 *
 * If eye/target not set, auto-fits from model bounding sphere.
 */
int render_camera_setup(render_camera_t* cam,
                        const render_config_t* cfg,
                        const alea_system_t* sys);

/**
 * @brief Generate camera ray for pixel (px, py) with sub-pixel offset (sx, sy)
 *
 * @param sx, sy Sub-pixel offset in [0,1] (0.5, 0.5 = pixel center)
 */
void render_camera_ray(const render_camera_t* cam,
                       int width, int height,
                       double px, double py,
                       double* ox, double* oy, double* oz,
                       double* dx, double* dy, double* dz);

/**
 * @brief Create framebuffer
 */
render_framebuffer_t* render_framebuffer_create(int width, int height, int aux);

/**
 * @brief Free framebuffer
 */
void render_framebuffer_free(render_framebuffer_t* fb);

/**
 * @brief Main render function — renders the model into the framebuffer
 *
 * @param sys CSG system (must be loaded and indexed)
 * @param cfg Render configuration
 * @param cam Camera (must be set up)
 * @param fb Output framebuffer
 * @return 0 on success, -1 on error
 */
int render_scene(const alea_system_t* sys,
                 const render_config_t* cfg,
                 const render_camera_t* cam,
                 render_framebuffer_t* fb);

/**
 * @brief Post-process: darken edges where cell ID changes
 */
void render_edge_darken(render_framebuffer_t* fb);

/**
 * @brief Convert float framebuffer to 8-bit RGB pixels
 *
 * @param fb Source framebuffer
 * @param pixels Output: width*height*3 bytes (caller-allocated)
 */
void render_tonemap(const render_framebuffer_t* fb, uint8_t* pixels);

/**
 * @brief Get material/cell/universe color from palette
 */
void render_get_color(int id, render_color_mode_t mode,
                      const render_config_t* cfg,
                      float* r, float* g, float* b);

/**
 * @brief Write PNG image (minimal, no dependencies)
 */
int render_write_png(const char* filename, const uint8_t* pixels,
                     int width, int height);

/**
 * @brief Write BMP image
 */
int render_write_bmp(const char* filename, const uint8_t* pixels,
                     int width, int height);

/**
 * @brief Write PPM image
 */
int render_write_ppm(const char* filename, const uint8_t* pixels,
                     int width, int height);

/**
 * @brief Write image, auto-detecting format from extension
 */
int render_write_image(const char* filename, const uint8_t* pixels,
                       int width, int height);

/**
 * @brief Write auxiliary maps (depth, cellid, matid, normal)
 *
 * @param base_filename Base filename (without extension)
 * @param fb Framebuffer with auxiliary data
 */
int render_write_aux(const char* base_filename, const render_framebuffer_t* fb);

/**
 * @brief Load custom color palette from file
 */
int render_load_colors(const char* filename, render_config_t* cfg);

#ifdef __cplusplus
}
#endif

#endif /* ALEA_RENDER_H */
