// SPDX-License-Identifier: MPL-2.0
/* Exercise asynchronous regional publication, failure and quit with the
 * normal CLI, SDL event loop and render/save workers. */
#include <SDL.h>
#include <assert.h>
#include <stdatomic.h>
#include <string.h>
#include "alea_raycast.h"
#include "alea_render.h"

static atomic_int builds, prepared_frames, build_done, first_image;
static int rotated, queued, first_render = 1;
static Uint64 uploaded_revision;
static const char* scenario;
static int fixed_hierarchy;

static void key(SDL_Keycode code) {
    SDL_Event event = {.type=SDL_KEYDOWN};
    event.key.keysym.sym = code;
    assert(SDL_PushEvent(&event) == 1);
}

static int slow_first(alea_system_t* sys, const render_config_t* cfg,
                      const render_camera_t* cam, render_framebuffer_t* fb) {
    if (first_render) {
        first_render = 0;
        SDL_Delay(200);
    }
    return render_scene(sys, cfg, cam, fb);
}

static alea_prepared_queries_t* test_build(alea_system_t* sys, const alea_bbox_t* box) {
    if (!strcmp(scenario, "pick")) return alea_prepared_queries_create(sys, box);
    assert(atomic_load(&first_image));
    atomic_fetch_add(&builds, 1);
    if (!strcmp(scenario, "close")) {
        key(SDLK_q);
        while (!alea_interrupted()) SDL_Delay(1);
        atomic_store(&build_done, 1);
        return NULL;
    }
    alea_prepared_queries_t* idx = alea_prepared_queries_create(sys,
        strcmp(scenario, "fail") ? box : NULL);
    atomic_store(&build_done, 1);
    return idx;
}

static int test_progressive(alea_system_t* sys, alea_prepared_queries_t* index,
                            const render_config_t* cfg, const render_camera_t* cam,
                            render_framebuffer_t* fb, const render_updates_t* updates) {
    if (first_render) {
        first_render = 0;
        SDL_Delay(200);
    }
    if (index) {
        atomic_fetch_add(&prepared_frames, 1);
        if (!strcmp(scenario, "renderfail")) return -1;
    }
    return alea_render_scene_progressive(sys, index, cfg, cam, fb, updates);
}

static int upload(SDL_Texture* tex, const SDL_Rect* rect, const void* pixels,
                  int pitch, Uint64 revision) {
    int rc = SDL_UpdateTexture(tex, rect, pixels, pitch);
    if (!rc) uploaded_revision = revision;
    return rc;
}

static void present(SDL_Renderer* renderer) {
    SDL_RenderPresent(renderer);
    if (uploaded_revision) atomic_store(&first_image, 1);
    if (fixed_hierarchy && uploaded_revision && !queued) {
        queued = 1;
        key(SDLK_s);
        key(SDLK_q);
    }
    if (!queued && atomic_load(&build_done) && strcmp(scenario, "close")) {
        if (!rotated) {
            rotated = 1;
            key(SDLK_RIGHT);
        } else if (uploaded_revision == 2 &&
                   (!strcmp(scenario, "fail") || atomic_load(&prepared_frames))) {
            queued = 1;
            key(SDLK_s);
            key(SDLK_q);
        }
    }
}

static void pick_present(SDL_Renderer*, void*, SDL_Rect, int);
static SDL_Thread* start_thread(SDL_ThreadFunction, const char*, void*);
#define SDL_CreateThread start_thread
#define render_scene slow_first
#define alea_prepared_queries_create test_build
#define alea_render_scene_progressive test_progressive
#define SDL_UpdateTexture(t,r,p,s) upload(t,r,p,s,job.revision)
#define SDL_RenderPresent(r) pick_present(r, state, image, current_preview)
#define main viewer_main
#include "../../examples/c/render3d.c"
#undef main
#undef SDL_CreateThread
#undef SDL_RenderPresent

/* Click an off-center pixel on the first preview through the scaled SDL
 * image. Its world point must project back to that pixel after orbiting. */
static void pick_present(SDL_Renderer* renderer, void* state, SDL_Rect image, int preview) {
    if (strcmp(scenario, "pick")) { present(renderer); return; }
    static int clicked;
    SDL_RenderPresent(renderer);
    render_viewer_t* v = state;
    const int x = 36, y = 30;
    if (!clicked && preview) {
        SDL_Event e = {.type = SDL_MOUSEBUTTONDOWN};
        e.button.button = SDL_BUTTON_LEFT;
        e.button.x = image.x + (int)((x + 0.5) * image.w / v->cfg.width);
        e.button.y = image.y + (int)((y + 0.5) * image.h / v->cfg.height);
        assert(SDL_PushEvent(&e) == 1);
        e.type = SDL_MOUSEBUTTONUP;
        assert(SDL_PushEvent(&e) == 1);
        clicked = 1;
        atomic_store(&first_image, 1);
    } else if (v->has_pivot && !rotated) {
        assert(v->pick.cell_id == 1 && v->pick.material_id == 1);
        assert(v->pick.clip_plane == 0 && v->pick.surface_id == -1);
        assert(fabs(v->pick.point[0]) < 1e-8);
        key(SDLK_RIGHT);
        rotated = 1;
    } else if (rotated && uploaded_revision == 2 && !queued) {
        render_camera_t cam;
        assert(render_camera_setup(&cam, &v->cfg, v->sys) == 0);
        double o[3], d[3];
        render_camera_ray(&cam, v->cfg.width, v->cfg.height, x+0.5, y+0.5,
                          o, o+1, o+2, d, d+1, d+2);
        double t = (v->pick.point[0] - o[0]) / d[0];
        for (int i = 0; i < 3; ++i) assert(fabs(o[i] + t*d[i] - v->pick.point[i]) < 1e-7);
        queued = 1;
        key(SDLK_q);
    }
}

static SDL_Thread* start_thread(SDL_ThreadFunction fn, const char* name, void* data) {
    /* A slow preview must be visible before the scheduler commits to a full
     * CPU refinement. An input arriving on presentation can then supersede it. */
    if (!strcmp(name, "alea-render") && !atomic_load(&first_image))
        assert(((alea_view_job_t*)data)->preview);
    return SDL_CreateThread(fn, name, data);
}

int main(int argc, char** argv) {
    render_view_build_t build = {0};
    render_viewer_t view = {.build=&build};
    render_view_preview_shown(&view, ALEA_VIEW_REFINE_DELAY_MS - 1);
    assert(!build.attempted);
    scenario = getenv("ALEA_TEST_BUILD");
    assert(scenario);
    fixed_hierarchy = !strcmp(scenario, "cpu") || !strcmp(scenario, "xray");
    int rc = viewer_main(argc, argv);
    assert(rc == 0);
    if (!strcmp(scenario, "pick")) { assert(queued); return rc; }
    assert(atomic_load(&builds) == !fixed_hierarchy);
    assert(atomic_load(&build_done) == !fixed_hierarchy);
    int attempts = atomic_load(&prepared_frames);
    if (!strcmp(scenario, "renderfail")) assert(attempts == 1);
    else assert((attempts > 0) == !strcmp(scenario, "ready"));
    assert(!alea_interrupted());
    return rc;
}
