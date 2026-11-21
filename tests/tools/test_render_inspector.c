// SPDX-License-Identifier: MPL-2.0
/* Drive the real probe and material controls with SDL mouse/text events. Locate controls in
 * Nuklear's text commands so font metrics and window layout can change. */
#include <SDL.h>
#include <assert.h>
static void locate_controls(void);
#define SDL_RenderGeometryRaw(...) (locate_controls(), SDL_RenderGeometryRaw(__VA_ARGS__))
#include "../../tools/render_inspector.c"
#undef SDL_RenderGeometryRaw

static const char* wanted;
static SDL_Point button;
static void locate_controls(void) {
    if (!wanted) return;
    const struct nk_command* cmd;
    nk_foreach(cmd, &sdl.ctx) {
        if (cmd->type != NK_COMMAND_TEXT) continue;
        const struct nk_command_text* t = (const struct nk_command_text*)cmd;
        if (strlen(wanted) == (size_t)t->length && !memcmp(t->string, wanted, t->length)) {
            button = (SDL_Point){t->x + 5, t->y + t->h/2};
            /* A numeric search string can also appear above the list row. */
        }
    }
}

/* The fixture picks the center of a radius-five sphere on an X cut face.
 * Even after another material is isolated in 3D, the XY section shows the
 * picked cell at the center and background in a corner. */
static int saw_reference, saw_empty_plane;
static void check_slice(inspector_t* p, const section_view_t* view, unsigned char* rgb) {
    if (!saw_reference && p->window && view->revision == p->section.revision) {
        assert(rgb && view->span > 10 && view->width == SECTION_PREVIEW_SIZE);
        unsigned char* center = rgb + 3 * ((view->height/2)*view->width + view->width/2);
        assert(center[0] == 255 && center[1] == 191 && center[2] == 50);
        assert(rgb[0] == rgb[1] && rgb[1] == rgb[2]);
        saw_reference = 1;
    }
    if (p->window && view->revision == p->section.revision && view->kind == SECTION_XY && view->offset > 5) {
        assert(rgb);
        for (int i = 0; i < view->width*view->height; ++i)
            assert(rgb[3*i] == 32 && rgb[3*i+1] == 32 && rgb[3*i+2] == 32);
        saw_empty_plane = 1;
    }
    /* Compare actual sampled pixels with both fixture spheres in world space,
     * including bands after the first and geometry away from the view center. */
    if (rgb && p->window && view->revision == p->section.revision) {
        double units = view->span / view->width;
        for (int row = 1; row < view->height; row += 7) {
            for (int col = 2; col < view->width; col += 7) {
                double u = view->u + (col + .5 - view->width*.5) * units;
                double v = view->v + (view->height*.5 - row - .5) * units;
                double xyz[3];
                for (int i = 0; i < 3; ++i)
                    xyz[i] = view->plane.origin[i] + view->offset * view->plane.normal[i]
                           + u * view->plane.u_axis[i] + v * view->plane.v_axis[i];
                double d1 = hypot(hypot(xyz[0], xyz[1]), xyz[2]) - 5;
                double d2 = hypot(hypot(xyz[0], xyz[1]-7), xyz[2]-1) - 1.6;
                if (fabs(d1) < 1e-6 || fabs(d2) < 1e-6) continue;
                unsigned char* pixel = rgb + 3*(row*view->width+col);
                assert((pixel[0]==255 && pixel[1]==191 && pixel[2]==50) == (d1<0));
                assert((pixel[0]==32 && pixel[1]==32 && pixel[2]==32) == (d1>0 && d2>0));
            }
        }
    }
    inspector_show_slice(p, view, rgb);
}

#define inspector_show_slice check_slice
static void advance(void*, SDL_Renderer*, SDL_Rect, int);
#define SDL_RenderPresent(r) advance(state, r, image, current_preview)
#define main viewer_main
#include "../../examples/c/render3d.c"
#undef main
#undef SDL_RenderPresent

static int step, release, finished, other_material;
static SDL_Point last;
static void mouse(Uint32 type, SDL_Point point) {
    SDL_Event event = {.type = SDL_MOUSEMOTION};
    event.motion.x = point.x;
    event.motion.y = point.y;
    assert(SDL_PushEvent(&event) == 1);
    event.type = type;
    event.button.button = SDL_BUTTON_LEFT;
    event.button.x = point.x;
    event.button.y = point.y;
    assert(SDL_PushEvent(&event) == 1);
}
static void click(SDL_Point p) {
    mouse(SDL_MOUSEBUTTONDOWN, p);
    last = p;
    release = 1;
    ++step;
}
static void key(SDL_Keycode k) {
    SDL_Event event = {.type = SDL_KEYDOWN};
    event.key.keysym.sym = k;
    assert(SDL_PushEvent(&event) == 1);
    ++step;
}
static void section_key(inspector_t* p, SDL_Keycode key) {
    SDL_Event event = {.type = SDL_KEYDOWN};
    event.key.windowID = SDL_GetWindowID(p->window);
    event.key.keysym.sym = key;
    assert(SDL_PushEvent(&event) == 1);
}
static int control(const char* label) {
    if (!wanted || strcmp(wanted, label)) {
        wanted = label; button.x = -1;
        return 0;
    }
    return button.x >= 0;
}

static void advance(void* state, SDL_Renderer* renderer, SDL_Rect image, int preview) {
    SDL_RenderPresent(renderer);
    render_viewer_t* v = state;
    inspector_t* p = v->inspector;
    static int previous = -1;
    static Uint32 entered;
    if (previous != step) { previous = step; entered = SDL_GetTicks(); }
    if (step && SDL_GetTicks() - entered > 12000) {
        fprintf(stderr, "Inspector stalled at step %d looking for '%s'\n", step, wanted ? wanted : "");
        abort();
    }
    if (release) { mouse(SDL_MOUSEBUTTONUP, last); release = 0; return; }
    switch (step) {
    case 0:
        assert(!inspector_height(p) && !p->ctx);
        if (preview) click((SDL_Point){image.x + image.w/2 + image.w/12, image.y + image.h/2 - image.h/12});
        break;
    case 1:
        if (v->has_pivot) {
            assert(!inspector_height(p));
            click((SDL_Point){image.x + image.w/2, image.y + image.h/2});
        }
        break;
    case 2:
        if (p->has_previous) {
            assert(memcmp(p->previous_point, p->pick.point, sizeof(p->previous_point)) != 0);
            key(SDLK_i);
        }
        break;
    case 3:
        if (control("Hide cell")) {
            SDL_Window* win = sdl.win;
            SDL_SetWindowSize(win, 640, 480);
            wanted = NULL;
            ++step;
        }
        break;
    case 4: {
        static char distance[80];
        double* prev = p->previous_point;
        double* hit = p->pick.point;
        snprintf(distance, sizeof(distance), "From previous pick   %.5g cm",
                 hypot(hypot(hit[0]-prev[0], hit[1]-prev[1]), hit[2]-prev[2]));
        if (control(distance)) ++step;
        break;
    }
    case 5:
        if (control("Hide cell")) click(button);
        break;
    case 6:
        if (!p->hidden_count) break;
        assert(v->cfg.cell_filter.mode == RENDER_FILTER_EXCLUDE);
        assert(v->cfg.cell_filter.ids[0] == v->pick.cell_id);
        if (control("Show all")) click(button);
        break;
    case 7:
        if (p->hidden_count) break;
        assert(v->cfg.cell_filter.mode == RENDER_FILTER_ALL);
        if (control("Isolate cell")) click(button);
        break;
    case 8:
        if (!p->only_cell) break;
        assert(v->cfg.cell_filter.mode == RENDER_FILTER_INCLUDE);
        key(SDLK_m);
        break;
    case 9:
        if (control("Material ID or name")) {
            int w, h;
            SDL_GetRendererOutputSize(renderer, &w, &h);
            click((SDL_Point){w/2, button.y});
        }
        break;
    case 10: {
        if (!p->editing) break;
        SDL_Event e = {.type=SDL_TEXTINPUT};
        for (size_t i = 0; i < alea_material_count(v->sys); ++i) {
            int id = alea_material_get_id(v->sys, (int)i);
            if (id > 0 && id != v->pick.material_id) { other_material = id; break; }
        }
        assert(other_material);
        snprintf(e.text.text, sizeof(e.text.text), "%d", other_material);
        assert(SDL_PushEvent(&e) == 1);
        ++step;
        break;
    }
    case 11: {
        if (!p->search[0]) break;
        assert(p->match_count > 0);
        int idx = p->matches[0];
        const char* name = alea_material_get_name(v->sys, idx);
        static char row[160];
        snprintf(row, sizeof(row), "%d%s%s", other_material, name ? "  " : "", name ? name : "");
        if (control(row)) click(button);
        break;
    }
    case 12: {
        if (p->material != other_material) break;
        static char label[48];
        snprintf(label, sizeof(label), "Isolate material %d", other_material);
        if (control(label)) click(button);
        break;
    }
    case 13:
        if (p->only_material != other_material) break;
        assert(v->cfg.material_filter.mode == RENDER_FILTER_INCLUDE && !p->only_cell);
        key(SDLK_i);
        break;
    case 14:
        if (control("XY section")) click(button);
        break;
    case 15:
        if (!p->slice) break;
        assert(saw_reference && p->view == PROBE_PICK && p->window != sdl.win);
        assert(!p->section_error);
        key(SDLK_i); /* Hide the 3D controls while the section stays open. */
        break;
    case 16: {
        assert(!inspector_height(p) && p->window);
        Uint32 id = SDL_GetWindowID(p->window);
        SDL_Event drag = {.type = SDL_MOUSEMOTION};
        drag.motion.windowID = id; drag.motion.state = SDL_BUTTON_LMASK;
        drag.motion.xrel = -10; drag.motion.yrel = 12;
        assert(SDL_PushEvent(&drag) == 1);
        SDL_Event wheel = {.type = SDL_MOUSEWHEEL};
        wheel.wheel.windowID = id; wheel.wheel.y = 1;
        assert(SDL_PushEvent(&wheel) == 1);
        for (int i = 0; i < 12; ++i) section_key(p, SDLK_PAGEUP);
        ++step;
        break;
    }
    case 17:
        if (p->displayed.revision != p->section.revision) break;
        assert(p->displayed.u > p->initial.u && p->displayed.v > p->initial.v);
        assert(saw_empty_plane && p->displayed.offset > 5 && p->displayed.span < p->initial.span);
        assert(!memcmp(v->cfg.eye, v->initial.eye, sizeof(v->cfg.eye)));
        SDL_SetWindowSize(p->window, 700, 375);
        ++step;
        break;
    case 18:
        if (p->displayed.width != 700) break;
        assert(p->displayed.height == 375 - SECTION_FOOTER);
        section_key(p, SDLK_r);
        ++step;
        break;
    case 19:
        if (p->displayed.revision != p->section.revision) break;
        assert(p->displayed.span == p->initial.span && p->displayed.offset == p->initial.offset);
        assert(p->displayed.u == p->initial.u && p->displayed.v == p->initial.v);
        section_key(p, SDLK_ESCAPE);
        ++step;
        break;
    case 20:
        assert(!p->window && !p->slice && v->has_pivot && p->only_material == other_material);
        key(SDLK_i);
        break;
    case 21: case 23: case 25: case 27: case 29: {
        const char* labels[] = {"XZ section", "YZ section", "X up", "Y up", "Z up"};
        if (control(labels[(step-21)/2])) click(button);
        break;
    }
    case 22: case 24: case 26: case 28: case 30: {
        if (!p->window || p->displayed.revision != p->section.revision) break;
        const section_view_t* s = &p->displayed;
        if (step == 22 || step == 24) {
            assert(s->kind == (step == 22 ? SECTION_XZ : SECTION_YZ));
            assert(s->plane.normal[s->kind] == 1 && s->plane.v_axis[2] == 1);
            assert(s->plane.u_axis[step == 22 ? 0 : 1] == 1);
            assert(!memcmp(s->plane.origin, p->pick.point, sizeof(p->pick.point)));
        } else {
            int up = (step-26)/2;
            assert(s->kind == SECTION_TWO_PICKS && s->up_axis == up);
            assert(fabs(s->plane.v_axis[up]-1) < 1e-12);
            for (int k = 0; k < 2; ++k) {
                const double* hit = k ? p->pick.point : p->previous_point;
                double normal = 0, horizontal = 0;
                for (int i = 0; i < 3; ++i) {
                    double d = hit[i]-s->plane.origin[i];
                    normal += d*s->plane.normal[i];
                    horizontal += d*s->plane.u_axis[i];
                }
                assert(fabs(normal) < 1e-12);
                assert(k ? horizontal > 0 : horizontal < 0);
                double side = k ? .5 : -.5;
                assert(fabs(horizontal - side*s->pick_du) < 1e-12);
                assert(fabs(hit[up]-s->plane.origin[up] - side*s->pick_dv) < 1e-12);
                assert(fabs(horizontal) < .5*s->span);
                assert(fabs(hit[up]-s->plane.origin[up]) < .5*s->span*s->height/s->width);
            }
        }
        if (step == 30) section_key(p, SDLK_PAGEUP);
        ++step;
        break;
    }
    case 31: {
        if (p->displayed.revision != p->section.revision || p->displayed.offset <= 0) break;
        SDL_Event close = {.type = SDL_WINDOWEVENT};
        close.window.windowID = SDL_GetWindowID(p->window);
        close.window.event = SDL_WINDOWEVENT_CLOSE;
        assert(SDL_PushEvent(&close) == 1);
        /* A queued event for a destroyed child must never close the 3D view. */
        section_key(p, SDLK_q);
        ++step;
        break;
    }
    case 32:
        assert(!p->window && v->has_pivot && p->only_material == other_material);
        finished = 1;
        SDL_Event e = {.type = SDL_WINDOWEVENT};
        e.window.windowID = SDL_GetWindowID(sdl.win);
        e.window.event = SDL_WINDOWEVENT_CLOSE;
        assert(SDL_PushEvent(&e) == 1);
        break;
    }
}

int main(int argc, char** argv) {
    double a[3] = {1,2,3}, b[3] = {1,2,7}, across, rise;
    alea_slice_plane_t plane;
    assert(section_from_picks(a, a, 0, &plane, &across, &rise));
    assert(section_from_picks(a, b, 2, &plane, &across, &rise));
    int rc = viewer_main(argc, argv);
    assert(!rc && finished);
    puts("SDL probe distance, material search, visibility and orthogonal/two-pick section navigation: PASS");
    return rc;
}
