// SPDX-License-Identifier: MPL-2.0
#include "render_inspector.h"
#include <float.h>
#include "sdl_text.h"
#if !SDL_VERSION_ATLEAST(2, 0, 18)
#error "The Nuklear SDL renderer requires SDL 2.0.18 or later"
#endif
#include <ctype.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define NK_INCLUDE_FIXED_TYPES
#define NK_INCLUDE_STANDARD_IO
#define NK_INCLUDE_STANDARD_VARARGS
#define NK_INCLUDE_DEFAULT_ALLOCATOR
#define NK_INCLUDE_VERTEX_BUFFER_OUTPUT
#define NK_INCLUDE_FONT_BAKING
#define NK_INCLUDE_DEFAULT_FONT
#define NK_IMPLEMENTATION
#include "../vendor/nuklear/nuklear.h"
#define NK_SDL_RENDERER_IMPLEMENTATION
#include "../vendor/nuklear/nuklear_sdl_renderer.h"

struct inspector {
    const alea_system_t* sys;
    struct nk_context* ctx;
    int captured, editing, unfocus, selected, material;
    render_pick_t pick;
    probe_view_t view;
    double previous_point[3];
    int has_previous;
    int* hidden;
    size_t hidden_count;
    int only_cell, only_material;
    char search[80], previous[80];
    int* matches;
    int match_count;
    SDL_Window* window;
    SDL_Renderer* renderer;
    SDL_Texture* slice;
    section_view_t section, initial, displayed;
    SDL_atomic_t revision;
    int request, preview, redraw;
    const char* section_error;
    const char* pick_error;
};

#include "render_section.h"

inspector_t* inspector_create(const alea_system_t* sys) {
    inspector_t* p = calloc(1, sizeof(*p));
    if (!p) return NULL;
    p->sys = sys;
    p->match_count = -1;
    return p;
}

void inspector_close_graphics(inspector_t* p) {
    close_section(p);
    if (p->ctx) nk_sdl_shutdown();
    p->ctx = NULL;
}

void inspector_destroy(inspector_t* p) {
    if (!p) return;
    free(p->matches);
    free(p->hidden);
    free(p);
}

int inspector_active(const inspector_t* p, probe_view_t view) {
    return p->view == view;
}

void inspector_toggle(inspector_t* p, probe_view_t view) {
    p->view = inspector_active(p, view) ? PROBE_CLOSED : view;
    p->editing = p->captured = 0;
    p->unfocus = 1;
}

int inspector_height(const inspector_t* p) {
    switch (p->view) {
    case PROBE_PICK: return p->selected ? 143 : 52;
    case PROBE_MATERIALS: return 166;
    default: return 0;
    }
}

void inspector_select(inspector_t* p, const render_pick_t* hit) {
    if (p->selected) {
        memcpy(p->previous_point, p->pick.point, sizeof(p->previous_point));
        p->has_previous = 1;
    }
    p->pick_error = NULL;
    p->pick = *hit;
    p->selected = 1;
    p->material = hit->material_id;
}

void inspector_clear(inspector_t* p) {
    p->hidden_count = 0;
    p->only_cell = p->only_material = 0;
    p->selected = p->material = p->has_previous = 0;
    p->pick_error = NULL;
}

int inspector_event(inspector_t* p, SDL_Event* e, SDL_Rect rect) {
    if (!p->ctx || rect.h <= 0) return 0;
    SDL_Point pos;
    switch (e->type) {
    case SDL_MOUSEBUTTONDOWN: case SDL_MOUSEBUTTONUP:
        pos = (SDL_Point){e->button.x, e->button.y};
        break;
    case SDL_MOUSEMOTION:
        pos = (SDL_Point){e->motion.x, e->motion.y};
        break;
    case SDL_MOUSEWHEEL:
        SDL_GetMouseState(&pos.x, &pos.y);
        break;
    case SDL_TEXTINPUT:
        // SDL can deliver a whole pasted or composed string in one event.
        if (!p->editing) return 0;
        for (int i = 0; e->text.text[i];) {
            nk_rune rune;
            int n = nk_utf_decode(e->text.text + i, &rune, (int)strlen(e->text.text + i));
            if (n <= 0) break;
            nk_input_unicode(p->ctx, rune);
            i += n;
        }
        return 1;
    case SDL_KEYDOWN: case SDL_KEYUP:
        if (!p->editing) return 0;
        if (e->type == SDL_KEYDOWN && (e->key.keysym.sym == SDLK_ESCAPE ||
            e->key.keysym.sym == SDLK_RETURN || e->key.keysym.sym == SDLK_KP_ENTER)) {
            p->editing = 0;
            p->unfocus = 1;
            return 1;
        }
        nk_sdl_handle_event(e);
        return 1;
    case SDL_WINDOWEVENT:
        if (e->window.event == SDL_WINDOWEVENT_FOCUS_LOST) {
            p->captured = p->editing = 0;
            p->unfocus = 1;
            nk_input_button(p->ctx, NK_BUTTON_LEFT, -1, -1, 0);
        }
        return 0;
    default: return 0;
    }
    int inside = SDL_PointInRect(&pos, &rect);
    int consume = inside || p->captured;
    if (e->type == SDL_MOUSEBUTTONDOWN) {
        p->captured = inside;
        if (!inside) { p->editing = 0; p->unfocus = 1; }
    }
    if (consume || e->type == SDL_MOUSEMOTION)
        nk_sdl_handle_event(e);
    if (e->type == SDL_MOUSEBUTTONUP) p->captured = 0;
    return consume;
}

static int contains_text(const char* text, const char* query) {
    if (!*query) return 1;
    for (; *text; ++text) {
        size_t i = 0;
        while (query[i] && text[i] &&
               tolower((unsigned char)text[i]) == tolower((unsigned char)query[i])) ++i;
        if (!query[i]) return 1;
    }
    return 0;
}

static int search_materials(inspector_t* p) {
    if (p->match_count >= 0 && !strcmp(p->previous, p->search)) return 0;
    size_t n = alea_material_count(p->sys);
    if (n && !p->matches) {
        p->matches = malloc(n * sizeof(*p->matches));
        if (!p->matches) return -1;
    }
    p->match_count = 0;
    for (size_t i = 0; i < n; ++i) {
        char id[24];
        snprintf(id, sizeof(id), "%d", alea_material_get_id(p->sys, (int)i));
        const char* name = alea_material_get_name(p->sys, (int)i);
        if (contains_text(id, p->search) || (name && contains_text(name, p->search)))
            p->matches[p->match_count++] = (int)i;
    }
    strcpy(p->previous, p->search);
    return 0;
}

static int hide_selected(inspector_t* p) {
    int id = p->pick.cell_id;
    size_t i = 0;
    while (i < p->hidden_count && p->hidden[i] < id) ++i;
    if (i == p->hidden_count || p->hidden[i] != id) {
        int* ids = realloc(p->hidden, (p->hidden_count + 1) * sizeof(*ids));
        if (!ids) return -1;
        p->hidden = ids;
        memmove(ids+i+1, ids+i, (p->hidden_count-i)*sizeof(*ids));
        ids[i] = id;
        ++p->hidden_count;
    }
    p->only_cell = 0;
    return 1;
}

static void apply_visibility(inspector_t* p, render_config_t* cfg) {
    cfg->cell_filter = (render_id_filter_t){0};
    cfg->material_filter = (render_id_filter_t){0};
    if (p->only_cell)
        cfg->cell_filter = (render_id_filter_t){RENDER_FILTER_INCLUDE, &p->only_cell, 1};
    else if (p->hidden_count)
        cfg->cell_filter = (render_id_filter_t){RENDER_FILTER_EXCLUDE, p->hidden, p->hidden_count};
    if (p->only_material)
        cfg->material_filter = (render_id_filter_t){RENDER_FILTER_INCLUDE, &p->only_material, 1};
}

static int show_all(inspector_t* p) {
    if (!p->only_cell && !p->only_material && !p->hidden_count) return 0;
    if (!nk_button_label(p->ctx, "Show all")) return 0;
    p->only_cell = p->only_material = 0;
    p->hidden_count = 0;
    return 1;
}

static void visibility_label(inspector_t* p) {
    if (p->only_cell) nk_labelf(p->ctx, NK_TEXT_RIGHT, "Showing cell %d", p->only_cell);
    else if (p->only_material && p->hidden_count)
        nk_labelf(p->ctx, NK_TEXT_RIGHT, "Material %d   %zu hidden", p->only_material, p->hidden_count);
    else if (p->only_material) nk_labelf(p->ctx, NK_TEXT_RIGHT, "Showing material %d", p->only_material);
    else if (p->hidden_count) nk_labelf(p->ctx, NK_TEXT_RIGHT, "%zu hidden", p->hidden_count);

}

static int pick_controls(inspector_t* p, const render_config_t* cfg) {
    struct nk_context* ctx = p->ctx;
    const render_pick_t* h = &p->pick;
    if (!p->selected) {
        nk_layout_row_dynamic(ctx, 24, 2);
        nk_label(ctx, "Click a surface to probe its boundary.", NK_TEXT_LEFT);
        return show_all(p);
    }
    nk_layout_row_dynamic(ctx, 18, 2);
    char boundary[40] = "";
    if (h->clip_plane >= 0) snprintf(boundary, sizeof(boundary), "   Clip %d", h->clip_plane + 1);
    else if (h->surface_id >= 0) snprintf(boundary, sizeof(boundary), "   Surface %d", h->surface_id);
    nk_labelf(ctx, NK_TEXT_LEFT, "Cell %d   Material %d%s", h->cell_id, h->material_id, boundary);
    visibility_label(p);
    nk_layout_row_dynamic(ctx, 18, 2);
    nk_labelf(ctx, NK_TEXT_LEFT, "XYZ (cm)   %.5g   %.5g   %.5g", h->point[0], h->point[1], h->point[2]);
    if (p->pick_error) nk_label(ctx, p->pick_error, NK_TEXT_RIGHT);
    else if (p->has_previous) {
        double dx = h->point[0] - p->previous_point[0], dy = h->point[1] - p->previous_point[1];
        double dz = h->point[2] - p->previous_point[2];
        nk_labelf(ctx, NK_TEXT_RIGHT, "From previous pick   %.5g cm", hypot(hypot(dx, dy), dz));
    }
    nk_layout_row_dynamic(ctx, 24, 3);
    if (nk_button_label(ctx, "Hide cell")) return hide_selected(p);
    if (nk_button_label(ctx, "Isolate cell")) {
        p->only_cell = h->cell_id; p->hidden_count = 0; p->only_material = 0;
        return 1;
    }
    if (show_all(p)) return 1;
    nk_layout_row_dynamic(ctx, 24, 3);
    if (nk_button_label(ctx, "XY section")) open_section(p, cfg, SECTION_XY, 1);
    if (nk_button_label(ctx, "XZ section")) open_section(p, cfg, SECTION_XZ, 2);
    if (nk_button_label(ctx, "YZ section")) open_section(p, cfg, SECTION_YZ, 2);
    if (p->has_previous) {
        nk_layout_row_dynamic(ctx, 24, 4);
        nk_label(ctx, "Last two picks", NK_TEXT_LEFT);
        const char* const labels[] = {"X up", "Y up", "Z up"};
        for (int axis = 0; axis < 3; ++axis)
            if (nk_button_label(ctx, labels[axis])) open_section(p, cfg, SECTION_TWO_PICKS, axis);
    } else {
        nk_layout_row_dynamic(ctx, 24, 1);
        nk_label(ctx, "Pick another point for a two-point section.", NK_TEXT_LEFT);
    }
    return 0;
}

static int material_controls(inspector_t* p, int height) {
    struct nk_context* ctx = p->ctx;
    char text[160];
    nk_layout_row_dynamic(ctx, 24, 3);
    nk_label(ctx, "Material ID or name", NK_TEXT_LEFT);
    p->editing = !!(nk_edit_string_zero_terminated(ctx, NK_EDIT_FIELD, p->search,
                                  sizeof(p->search), nk_filter_default) & NK_EDIT_ACTIVE);
    visibility_label(p);
    if (strcmp(p->previous, p->search)) nk_group_set_scroll(ctx, "Materials", 0, 0);
    if (search_materials(p)) return -1;
    struct nk_list_view list;
    nk_layout_row_dynamic(ctx, height - 76, 1);
    if (!p->match_count) nk_label(ctx, "No matching materials", NK_TEXT_LEFT);
    else if (nk_list_view_begin(ctx, &list, "Materials", NK_WINDOW_BORDER, 22, p->match_count)) {
        for (int row = list.begin; row < list.end; ++row) {
            int index = p->matches[row];
            int id = alea_material_get_id(p->sys, index);
            const char* name = alea_material_get_name(p->sys, index);
            snprintf(text, sizeof(text), "%d%s%s", id, name ? "  " : "", name ? name : "");
            nk_layout_row_dynamic(ctx, 22, 1);
            nk_bool selected = id == p->material;
            if (nk_selectable_label(ctx, text, NK_TEXT_LEFT, &selected)) p->material = id;
        }
        nk_list_view_end(&list);
    }
    nk_layout_row_dynamic(ctx, 24, 2);
    if (p->material) {
        snprintf(text, sizeof(text), "Isolate material %d", p->material);
        if (nk_button_label(ctx, text)) {
            p->only_material = p->material; p->only_cell = 0; p->hidden_count = 0;
            return 1;
        }
    } else nk_spacing(ctx, 1);
    return show_all(p);
}

int inspector_draw(inspector_t* p, SDL_Window* win, SDL_Renderer* renderer,
                   SDL_Rect rect, render_config_t* cfg) {
    if (!p->ctx) {
        if (rect.h <= 0) return 0;
        struct nk_font_atlas* atlas;
        p->ctx = nk_sdl_init(win, renderer);
        nk_sdl_font_stash_begin(&atlas);
        struct nk_font* font = nk_font_atlas_add_default(atlas, 16, NULL);
        nk_sdl_font_stash_end();
        if (!sdl.ogl.font_tex) return -1;
        if (font) nk_style_set_font(p->ctx, &font->handle);
        struct nk_color colors[NK_COLOR_COUNT];
        for (int i = 0; i < NK_COLOR_COUNT; ++i) colors[i] = nk_rgb(220, 220, 216);
        colors[NK_COLOR_TEXT] = colors[NK_COLOR_EDIT_CURSOR] = nk_rgb(25, 28, 30);
        colors[NK_COLOR_BORDER] = colors[NK_COLOR_SCROLLBAR_CURSOR] = nk_rgb(112, 116, 114);
        colors[NK_COLOR_BUTTON] = nk_rgb(200, 203, 200);
        colors[NK_COLOR_BUTTON_HOVER] = colors[NK_COLOR_EDIT] = nk_rgb(247, 248, 244);
        colors[NK_COLOR_SELECT_ACTIVE] = colors[NK_COLOR_BUTTON_ACTIVE] = nk_rgb(255, 201, 91);
        nk_style_from_table(p->ctx, colors);
        nk_input_begin(p->ctx);
    }
    nk_input_end(p->ctx);
    int action = 0;
    if (rect.h > 0) {
        if (nk_begin(p->ctx, "Probe tools", nk_rect(rect.x, rect.y, rect.w, rect.h),
                     NK_WINDOW_BORDER | NK_WINDOW_NO_SCROLLBAR)) {
            if (p->unfocus) { nk_edit_unfocus(p->ctx); p->unfocus = 0; }
            switch (p->view) {
            case PROBE_PICK: action = pick_controls(p, cfg); break;
            case PROBE_MATERIALS: action = material_controls(p, rect.h); break;
            default: break;
            }
        }
        nk_end(p->ctx);
        nk_sdl_render(NK_ANTI_ALIASING_ON);
    } else {
        p->editing = p->captured = 0;
        p->unfocus = 1;
        nk_clear(p->ctx);
    }
    nk_input_begin(p->ctx);
    if (action == 1) apply_visibility(p, cfg);
    draw_section(p);
    if (action == 1) return 1;
    if (!p->request) return 0;
    p->request = 0;
    return 2;
}
