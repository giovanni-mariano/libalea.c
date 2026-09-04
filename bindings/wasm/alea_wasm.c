// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
// SPDX-License-Identifier: MPL-2.0

#include "alea_wasm.h"
#include "alea.h"
#include "alea_mcnp.h"
#include "alea_render.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef __EMSCRIPTEN__
#include <emscripten/emscripten.h>
#define ALEA_WASM_EXPORT EMSCRIPTEN_KEEPALIVE
#else
#define ALEA_WASM_EXPORT
#endif

typedef struct {
    mcnp_model_t* model;
    alea_system_t* system;
    render_config_t config;
    render_framebuffer_t* framebuffer;
    uint8_t* rgb;
    uint8_t* rgba;
    double center[3];
    double target[3];
    double radius;
    int width;
    int height;
    int config_initialized;
    int initialized;
} alea_wasm_state_t;

static alea_wasm_state_t state;
static char last_error[256];

static void set_error(const char* message) {
    snprintf(last_error, sizeof(last_error), "%s",
             message && message[0] ? message : "unknown libalea error");
}

static void configure_renderer(int width, int height) {
    render_config_init(&state.config);
    state.config_initialized = 1;
    state.config.width = width;
    state.config.height = height;
    state.config.background[0] = 0.025f;
    state.config.background[1] = 0.040f;
    state.config.background[2] = 0.075f;
    state.config.fov = 42.0;
    state.config.eye_set = 1;
    state.config.target_set = 1;
    state.config.up[2] = 1.0;
    state.config.color_mode = RENDER_COLOR_MATERIAL;
    state.config.render_mode = RENDER_MODE_SOLID;
    state.config.shadows = 0;
    state.config.edges = 1;
    state.config.aa_samples = 1;
    state.config.tile_size = 32;
    state.config.ambient = 0.58f;
    state.config.diffuse = 0.42f;
    state.config.specular = 0.08f;
    state.config.cross_section_tint = 0.80f;
    state.config.log_level = 0;
}

ALEA_WASM_EXPORT void alea_wasm_destroy(void) {
    free(state.rgb);
    free(state.rgba);
    render_framebuffer_free(state.framebuffer);
    if (state.config_initialized) render_config_free(&state.config);
    mcnp_model_destroy(state.model);
    memset(&state, 0, sizeof(state));
}

ALEA_WASM_EXPORT int alea_wasm_init(int width, int height) {
    alea_wasm_destroy();
    last_error[0] = '\0';
    if (width < 16 || height < 16 || width > 2048 || height > 2048) {
        set_error("frame dimensions must be between 16 and 2048 pixels");
        return -1;
    }
    state.width = width;
    state.height = height;
    configure_renderer(width, height);
    const size_t pixels = (size_t)width * (size_t)height;
    state.framebuffer = render_framebuffer_create(width, height, 0);
    state.rgb = malloc(pixels * 3);
    state.rgba = malloc(pixels * 4);
    if (!state.framebuffer || !state.rgb || !state.rgba) {
        set_error("failed to allocate the WASM framebuffer");
        alea_wasm_destroy();
        return -1;
    }
    state.initialized = 1;
    return 0;
}

ALEA_WASM_EXPORT int alea_wasm_resize(int width, int height) {
    if (!state.initialized) {
        set_error("renderer is not initialized");
        return -1;
    }
    if (width < 16 || height < 16 || width > 2048 || height > 2048) {
        set_error("frame dimensions must be between 16 and 2048 pixels");
        return -1;
    }
    if (width == state.width && height == state.height) return 0;

    const size_t pixels = (size_t)width * (size_t)height;
    render_framebuffer_t* framebuffer = render_framebuffer_create(width, height, 0);
    uint8_t* rgb = malloc(pixels * 3);
    uint8_t* rgba = malloc(pixels * 4);
    if (!framebuffer || !rgb || !rgba) {
        render_framebuffer_free(framebuffer);
        free(rgb);
        free(rgba);
        set_error("failed to resize the WASM framebuffer");
        return -1;
    }

    render_framebuffer_free(state.framebuffer);
    free(state.rgb);
    free(state.rgba);
    state.framebuffer = framebuffer;
    state.rgb = rgb;
    state.rgba = rgba;
    state.width = width;
    state.height = height;
    state.config.width = width;
    state.config.height = height;
    return 0;
}

ALEA_WASM_EXPORT int alea_wasm_load_mcnp(const char* input, int length) {
    if (!state.initialized) {
        set_error("renderer is not initialized");
        return -1;
    }
    if (!input || length <= 0) {
        set_error("MCNP input is empty");
        return -1;
    }
    mcnp_model_t* model = mcnp_load_string(input, (size_t)length);
    if (!model) {
        set_error(alea_error());
        return -1;
    }
    alea_system_t* system = mcnp_model_system(model);
    if (!system || alea_build_universe_index(system) != 0 ||
        alea_prepare_query_acceleration(system) != 0) {
        set_error(alea_error());
        mcnp_model_destroy(model);
        return -1;
    }
    double cx, cy, cz, radius;
    if (alea_compute_bounding_sphere(system, 0.05, &cx, &cy, &cz, &radius) != 0 ||
        !isfinite(radius) || radius <= 0.0) {
        set_error("MCNP geometry has no finite, bounded cells");
        mcnp_model_destroy(model);
        return -1;
    }
    mcnp_model_destroy(state.model);
    state.model = model;
    state.system = system;
    state.center[0] = cx;
    state.center[1] = cy;
    state.center[2] = cz;
    memcpy(state.target, state.center, sizeof(state.target));
    state.radius = radius;
    last_error[0] = '\0';
    return 0;
}

ALEA_WASM_EXPORT int alea_wasm_render(double azimuth, double elevation,
                                      double distance_scale,
                                      double clip_fraction, int edges) {
    if (!state.initialized || !state.system) {
        set_error("load an MCNP model before rendering");
        return -1;
    }
    if (distance_scale < 1.2) distance_scale = 1.2;
    if (distance_scale > 12.0) distance_scale = 12.0;
    if (elevation < -1.45) elevation = -1.45;
    if (elevation > 1.45) elevation = 1.45;
    const double distance = state.radius * distance_scale;
    const double horizontal = cos(elevation) * distance;
    state.config.eye[0] = state.target[0] + cos(azimuth) * horizontal;
    state.config.eye[1] = state.target[1] + sin(azimuth) * horizontal;
    state.config.eye[2] = state.target[2] + sin(elevation) * distance;
    memcpy(state.config.target, state.target, sizeof(state.target));
    state.config.edges = edges != 0;
    state.config.num_clips = 0;
    if (clip_fraction >= 0.0) {
        if (clip_fraction > 1.0) clip_fraction = 1.0;
        render_clip_plane_t* clip = &state.config.clips[0];
        clip->normal[0] = 1.0;
        clip->normal[1] = 0.0;
        clip->normal[2] = 0.0;
        clip->d = -(state.center[0] + (clip_fraction - 0.5) *
                    2.0 * state.radius);
        state.config.num_clips = 1;
    }
    render_camera_t camera;
    if (render_camera_setup(&camera, &state.config, state.system) != 0 ||
        render_scene(state.system, &state.config, &camera, state.framebuffer) != 0) {
        set_error(alea_error());
        return -1;
    }
    if (state.config.edges) render_edge_darken(state.framebuffer);
    render_tonemap(state.framebuffer, state.rgb);
    const size_t pixels = (size_t)state.width * (size_t)state.height;
    for (size_t i = 0; i < pixels; ++i) {
        state.rgba[i * 4 + 0] = state.rgb[i * 3 + 0];
        state.rgba[i * 4 + 1] = state.rgb[i * 3 + 1];
        state.rgba[i * 4 + 2] = state.rgb[i * 3 + 2];
        state.rgba[i * 4 + 3] = 255;
    }
    return 0;
}

ALEA_WASM_EXPORT const uint8_t* alea_wasm_pixels(void) { return state.rgba; }
ALEA_WASM_EXPORT int alea_wasm_pixel_bytes(void) {
    return state.initialized ? state.width * state.height * 4 : 0;
}
ALEA_WASM_EXPORT int alea_wasm_width(void) { return state.width; }
ALEA_WASM_EXPORT int alea_wasm_height(void) { return state.height; }
ALEA_WASM_EXPORT int alea_wasm_cell_count(void) {
    return state.system ? (int)alea_cell_count(state.system) : 0;
}
ALEA_WASM_EXPORT int alea_wasm_surface_count(void) {
    return state.system ? (int)alea_surface_count(state.system) : 0;
}
ALEA_WASM_EXPORT double alea_wasm_center_x(void) { return state.center[0]; }
ALEA_WASM_EXPORT double alea_wasm_center_y(void) { return state.center[1]; }
ALEA_WASM_EXPORT double alea_wasm_center_z(void) { return state.center[2]; }
ALEA_WASM_EXPORT double alea_wasm_radius(void) { return state.radius; }
ALEA_WASM_EXPORT int alea_wasm_set_target(double x, double y, double z) {
    if (!state.system) {
        set_error("load an MCNP model before setting the camera target");
        return -1;
    }
    if (!isfinite(x) || !isfinite(y) || !isfinite(z)) {
        set_error("camera target coordinates must be finite numbers");
        return -1;
    }
    state.target[0] = x;
    state.target[1] = y;
    state.target[2] = z;
    return 0;
}
ALEA_WASM_EXPORT int alea_wasm_parallel_enabled(void) {
    return alea_parallel_enabled();
}
ALEA_WASM_EXPORT int alea_wasm_parallel_max_threads(void) {
    return alea_parallel_max_threads();
}
ALEA_WASM_EXPORT int alea_wasm_openmp_enabled(void) { return alea_openmp_enabled(); }
ALEA_WASM_EXPORT int alea_wasm_openmp_max_threads(void) {
    return alea_openmp_max_threads();
}
ALEA_WASM_EXPORT const char* alea_wasm_last_error(void) {
    return last_error[0] ? last_error : alea_error();
}
