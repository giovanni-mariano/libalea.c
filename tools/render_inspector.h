// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
// SPDX-License-Identifier: MPL-2.0
#ifndef ALEA_RENDER_INSPECTOR_H
#define ALEA_RENDER_INSPECTOR_H
#include <SDL.h>
#include "alea_render.h"
#include "alea_slice.h"

/* The three preset values are the fixed axes used by alea_slice_view_axis. */
typedef enum { SECTION_YZ = 0, SECTION_XZ, SECTION_XY, SECTION_TWO_PICKS } section_kind_t;

/* Immutable copy sent to the shared model worker for one section image. */
typedef struct {
    alea_slice_plane_t plane; /* Initial plane through the pick(s). */
    double u, v, offset, span; /* In-plane pan, normal displacement, horizontal span (cm). */
    double pick_du, pick_dv; /* Separation of the two picks in plane coordinates. */
    section_kind_t kind;
    int up_axis;
    int width, height, cell, revision;
} section_view_t;
typedef struct inspector inspector_t;
typedef enum { PROBE_CLOSED, PROBE_PICK, PROBE_MATERIALS } probe_view_t;
void inspector_toggle(inspector_t* p, probe_view_t view);
int inspector_active(const inspector_t* p, probe_view_t view);
int inspector_height(const inspector_t* p);
inspector_t* inspector_create(const alea_system_t* sys);
void inspector_destroy(inspector_t* p);
void inspector_close_graphics(inspector_t* p);
void inspector_select(inspector_t* p, const render_pick_t* pick);
void inspector_clear(inspector_t* p);
int inspector_event(inspector_t* p, SDL_Event* event, SDL_Rect rect);
/* Returns 1 for changed visibility, 2 to queue a slice, -1 on error. */
int inspector_draw(inspector_t* p, SDL_Window* window, SDL_Renderer* renderer,
                   SDL_Rect rect, render_config_t* cfg);
void inspector_section_event(inspector_t* p, Uint32 window_id, const SDL_Event* event);
section_view_t inspector_section_view(const inspector_t* p);
int inspector_slice(alea_system_t* sys, const render_config_t* cfg, inspector_t* p,
                    const section_view_t* view, unsigned char** rgb);
void inspector_show_slice(inspector_t* p, const section_view_t* view, unsigned char* rgb);
#endif
