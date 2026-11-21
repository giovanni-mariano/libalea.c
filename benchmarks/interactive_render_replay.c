// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
// SPDX-License-Identifier: MPL-2.0
/* Drive the actual SDL viewer, including texture upload and presentation.
 * Run with SDL_VIDEODRIVER=dummy for a reproducible software-display replay. */
#include <SDL.h>
#include <stdatomic.h>
#include <stdio.h>
#include <stdlib.h>
#ifdef __linux__
#include <sys/resource.h>
#endif
#include "alea_raycast.h"
#include "alea_render.h"

static Uint64 started, issued;
static int uploaded, frame;
static int waiting;
static int cancel_zoom, zoom_partial, zoom_cancelled;
static Uint64 last_partial_revision;
static Uint64 revision = 1;
static const char* directory;
/* 0 before preparation, 1 during the build, 2 on success, -1 on error. */
static atomic_int build_state;

static double seconds_since(Uint64 tick) {
    return (double)(SDL_GetPerformanceCounter() - tick) / SDL_GetPerformanceFrequency();
}

static alea_prepared_queries_t* build_measured(alea_system_t* sys, const alea_bbox_t* bounds) {
    atomic_store(&build_state, 1);
    Uint64 start = SDL_GetPerformanceCounter();
    alea_prepared_queries_t* index = alea_prepared_queries_create(sys, bounds);
    fprintf(stderr, "build_s=%.6f\n", seconds_since(start));
#ifdef __linux__
    struct rusage usage;
    if (getrusage(RUSAGE_SELF, &usage) == 0)
        fprintf(stderr, "build_complete_peak_kib=%ld\n", usage.ru_maxrss);
#endif
    atomic_store(&build_state, index ? 2 : -1);
    return index;
}

static int record_texture(SDL_Texture* tex, const SDL_Rect* rect, const void* pixels,
                          int pitch, Uint64 rendered, int preview, int partial) {
    int rc = SDL_UpdateTexture(tex, rect, pixels, pitch);
    if (rc || rendered != revision) return rc;
    if (partial) {
        if (last_partial_revision != revision) {
            last_partial_revision = revision;
            fprintf(stderr, "partial revision=%llu latency_s=%.6f\n",
                    (unsigned long long)revision, seconds_since(issued));
        }
        if (frame == 2 && !zoom_cancelled) zoom_partial = 1;
        return rc;
    }
    if (preview != (frame != 5 && frame != 10)) return rc;
    int w, h;
    char path[1024];
    if (SDL_QueryTexture(tex, NULL, NULL, &w, &h) ||
        snprintf(path, sizeof(path), "%s/%02d.rgb", directory, frame) >= (int)sizeof(path))
        exit(2);
    FILE* out = fopen(path, "wb");
    if (!out) exit(2);
    for (int y = 0; y < h; ++y)
        if (fwrite((const Uint8*)pixels + y*pitch, 3, w, out) != (size_t)w) exit(2);
    if (fclose(out)) exit(2);
    uploaded = 1;
    return 0;
}

static void push(SDL_Event event) {
    if (SDL_PushEvent(&event) != 1) exit(2);
}

static void camera_input(int step) {
    issued = SDL_GetPerformanceCounter();
    ++revision;
    if (step == 1 || step == 3) {
        SDL_Event e = {.type = SDL_KEYDOWN};
        e.key.keysym.sym = step == 1 ? SDLK_MINUS : SDLK_PLUS;
        push(e);
        return;
    }
    SDL_Event e = {.type = SDL_MOUSEBUTTONDOWN};
    e.button.button = SDL_BUTTON_LEFT;
    e.button.x = 360;
    e.button.y = 300;
    push(e);
    e.type = SDL_MOUSEMOTION;
    e.motion.state = SDL_BUTTON_LMASK;
    e.motion.xrel = step == 0 ? 70 : -58;
    e.motion.yrel = step == 0 ? -25 : 32;
    push(e);
    e.type = SDL_MOUSEBUTTONUP;
    e.button.button = SDL_BUTTON_LEFT;
    push(e);
}

static void record_present(SDL_Renderer* renderer) {
    SDL_RenderPresent(renderer);
    if (atomic_load(&build_state) < 0) exit(2);
    if (waiting && atomic_load(&build_state) == 2) {
        waiting = 0;
        camera_input(0);
    }
    if (cancel_zoom && frame == 2 && zoom_partial && !zoom_cancelled &&
        seconds_since(issued) >= 0.5 && !uploaded) {
        zoom_cancelled = 1;
        camera_input(3); /* Undo the zoom while its preview is still running. */
    }
    if (!uploaded) return;
    uploaded = 0;
    fprintf(stderr, "present %d source_s=%.6f latency_s=%.6f\n", frame,
            seconds_since(started), seconds_since(issued));
    ++frame;
    /* Initial preview, four camera changes, full refinement, then four more
     * camera changes and another full refinement. */
    if (cancel_zoom && frame == 3) push((SDL_Event){.type = SDL_QUIT});
    else if (frame <= 4) camera_input(frame-1);
    else if (frame == 6 && atomic_load(&build_state) == 1) waiting = 1;
    else if (frame >= 6 && frame <= 9) camera_input(frame-6);
    else if (frame == 11) push((SDL_Event){.type = SDL_QUIT});
}

#define SDL_UpdateTexture(t,r,p,s) record_texture(t,r,p,s,job.revision,job.preview,worker != NULL)
#define SDL_RenderPresent record_present
#define alea_prepared_queries_create build_measured
#define main viewer_main
#include "../examples/c/render3d.c"
#undef main

int main(int argc, char** argv) {
    directory = getenv("ALEA_REPLAY_DIR");
    if (!directory) return 2;
    cancel_zoom = getenv("ALEA_REPLAY_CANCEL_ZOOM") != NULL;
    started = issued = SDL_GetPerformanceCounter();
    int rc = viewer_main(argc, argv);
    return rc || (cancel_zoom ? frame != 3 || !zoom_cancelled : frame != 11);
}
