// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
// SPDX-License-Identifier: MPL-2.0

#ifndef ALEA_RENDER_VIEW_H
#define ALEA_RENDER_VIEW_H
#include "sdl_view.h"

typedef struct {
    alea_system_t* sys;
    render_config_t cfg, initial;
    render_camera_t cam;
    const char* output;
} render_viewer_t;

static int render_view_draw(void* state, uint8_t** pixels) {
    render_viewer_t* v = state;
    *pixels = NULL;
    render_framebuffer_t* fb = render_framebuffer_create(v->cfg.width, v->cfg.height, 1);
    uint8_t* rgb = malloc((size_t)v->cfg.width * v->cfg.height * 3);
    int rc = -1;
    if (fb && rgb && render_camera_setup(&v->cam, &v->cfg, v->sys) == 0 &&
        render_scene(v->sys, &v->cfg, &v->cam, fb) == 0) {
        postprocess_frame(&v->cfg, fb);
        render_tonemap(fb, rgb);
        *pixels = rgb;
        rgb = NULL;
        rc = 0;
    }
    render_framebuffer_free(fb);
    free(rgb);
    return rc;
}

static int render_view_save(const void* state, const uint8_t* pixels) {
    const render_viewer_t* v = state;
    int rc = render_write_image(v->output, pixels, v->cfg.width, v->cfg.height);
    fprintf(stderr, "%s %s\n", rc == 0 ? "Saved" : "Failed to save", v->output);
    return rc;
}

static int render_view_active(const void* state, SDL_Keycode key) {
    const render_viewer_t* v = state;
    return (key == SDLK_e && v->cfg.edges) || (key == SDLK_l && v->cfg.shadows);
}

/* Rotate a vector about a unit axis (Rodrigues' formula). */
static void view_rotate(double v[3], const double axis[3], double angle) {
    double c = cos(angle), s = sin(angle);
    double dot = v[0]*axis[0] + v[1]*axis[1] + v[2]*axis[2];
    double cross[3] = {axis[1]*v[2] - axis[2]*v[1],
                       axis[2]*v[0] - axis[0]*v[2],
                       axis[0]*v[1] - axis[1]*v[0]};
    for (int i = 0; i < 3; ++i)
        v[i] = v[i]*c + cross[i]*s + axis[i]*dot*(1-c);
}

static int render_view_input(void* state, const SDL_Event* event) {
    render_viewer_t* v = state;
    render_config_t* cfg = &v->cfg;
    /* Input can be coalesced between frames; derive axes from the latest state. */
    if (render_camera_setup(&v->cam, cfg, v->sys) != 0) return -1;
    double yaw = 0, pitch = 0, pan_x = 0, pan_y = 0, zoom = 1;
    if (event->type == SDL_MOUSEWHEEL) {
        int dy = event->wheel.y;
        if (event->wheel.direction == SDL_MOUSEWHEEL_FLIPPED) dy = -dy;
        zoom = dy > 0 ? 0.8 : dy < 0 ? 1.25 : 1;
    } else if (event->type == SDL_MOUSEMOTION) {
        if (event->motion.state & (SDL_BUTTON_RMASK | SDL_BUTTON_MMASK)) {
            pan_x = -2.0 * event->motion.xrel * v->cam.half_width / cfg->width;
            pan_y = 2.0 * event->motion.yrel * v->cam.half_height / cfg->height;
            if (cfg->fov != 0) {
                double distance = 0;
                for (int i = 0; i < 3; ++i)
                    distance += (cfg->eye[i]-cfg->target[i])*(cfg->eye[i]-cfg->target[i]);
                pan_x *= sqrt(distance); pan_y *= sqrt(distance);
            }
        } else if (event->motion.state & SDL_BUTTON_LMASK) {
            yaw = -event->motion.xrel * 0.008;
            pitch = -event->motion.yrel * 0.008;
        } else return 0;
    } else if (event->type == SDL_KEYDOWN) {
        switch (event->key.keysym.sym) {
        case SDLK_r: *cfg = v->initial; return 1;
        case SDLK_c: cfg->color_mode = (cfg->color_mode + 1) % 4; return 1;
        case SDLK_e: cfg->edges = !cfg->edges; return 1;
        case SDLK_l: cfg->shadows = !cfg->shadows; return 1;
        case SDLK_LEFT: yaw = 0.15; break;
        case SDLK_RIGHT: yaw = -0.15; break;
        case SDLK_UP: pitch = 0.15; break;
        case SDLK_DOWN: pitch = -0.15; break;
        case SDLK_EQUALS: case SDLK_PLUS: case SDLK_KP_PLUS: zoom = 0.8; break;
        case SDLK_MINUS: case SDLK_KP_MINUS: zoom = 1.25; break;
        default: return 0;
        }
    } else return 0;
    double offset[3];
    for (int i = 0; i < 3; ++i) offset[i] = cfg->eye[i] - cfg->target[i];
    if (cfg->fov == 0) {
        double height = cfg->ortho_height * zoom;
        if (height < 1e-10 || height > 1e20) return 0;
        cfg->ortho_height = height;
    } else {
        double distance = sqrt(offset[0]*offset[0]+offset[1]*offset[1]+offset[2]*offset[2]);
        if (distance * zoom < 1e-10 || distance * zoom > 1e20) return 0;
        for (int i = 0; i < 3; ++i) offset[i] *= zoom;
    }
    view_rotate(offset, v->cam.true_up, yaw);
    view_rotate(offset, v->cam.right, pitch);
    memcpy(cfg->up, v->cam.true_up, sizeof(cfg->up));
    view_rotate(cfg->up, v->cam.right, pitch);
    for (int i = 0; i < 3; ++i) {
        cfg->target[i] += pan_x*v->cam.right[i] + pan_y*v->cam.true_up[i];
        cfg->eye[i] = cfg->target[i] + offset[i];
    }
    return 1;
}

static int render_interactive(alea_system_t* sys, const render_config_t* cfg,
                              const render_camera_t* cam, const char* output) {
    if (cfg->width < 1 || cfg->height < 1 || cfg->width > 16384 || cfg->height > 16384)
        return -1;
    render_viewer_t v = {0};
    v.sys = sys; v.cfg = *cfg; v.cam = *cam; v.output = output;
    /* Freeze auto-fit so navigation is relative to the visible camera. */
    memcpy(v.cfg.eye, cam->eye, sizeof(v.cfg.eye));
    memcpy(v.cfg.target, cam->target, sizeof(v.cfg.target));
    memcpy(v.cfg.up, cam->true_up, sizeof(v.cfg.up));
    v.cfg.ortho_height = cam->ortho_height;
    v.cfg.eye_set = v.cfg.target_set = 1;
    v.initial = v.cfg; /* Borrowed config resources; freed only by the caller. */
    static const alea_view_button_t buttons[] = {
        {"Save", SDLK_s}, {"Reset", SDLK_r}, {"Zoom +", SDLK_PLUS},
        {"Zoom -", SDLK_MINUS}, {"Help", SDLK_h}, {"Quit", SDLK_q},
        {"Left", SDLK_LEFT}, {"Right", SDLK_RIGHT}, {"Up", SDLK_UP},
        {"Down", SDLK_DOWN}, {"Colors", SDLK_c}, {"Edges", SDLK_e},
        {"Shadows", SDLK_l}
    };
    return alea_sdl_view("Alea 3D",
        "Left drag / arrows: orbit\n"
        "Right / middle drag: pan   Wheel / +/-: zoom\n"
        "C: colors   E: edges   L: shadows\n"
        "R: reset   S: save   H: help   Q / Esc: quit",
        cfg->width, cfg->height, &v, sizeof(v), render_view_draw,
        render_view_save, render_view_input, render_view_active,
        buttons, (int)(sizeof(buttons)/sizeof(buttons[0])));
}
#endif
