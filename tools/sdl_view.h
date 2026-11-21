// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
// SPDX-License-Identifier: MPL-2.0

#ifndef ALEA_SDL_VIEW_H
#define ALEA_SDL_VIEW_H

#include "sdl_text.h"
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Only the SDL thread touches live state, windows, textures, or events. A job
 * owns a byte-copy of state and a newly allocated RGB image. Pointers inside a
 * snapshot must refer to immutable data, or worker-owned resources that only
 * the single worker accesses, kept alive until the worker joins.
 * One worker at a time also serializes access to the model's query caches. */
typedef int (*alea_view_draw_t)(void*, unsigned char**);
typedef int (*alea_view_save_t)(const void*, const unsigned char*);
typedef int (*alea_view_input_t)(void*, const SDL_Event*);
/* Optional quality hook runs on a job snapshot before dispatch. Saves always
 * receive preview=0; preview images must still have the requested dimensions. */
typedef void (*alea_view_prepare_t)(void*, int preview);
typedef int (*alea_view_active_t)(const void*, SDL_Keycode);
typedef struct { const char* label; SDL_Keycode key; } alea_view_button_t;
enum { ALEA_VIEW_REFINE_DELAY_MS = 150 };
/* Optional bottom strip for the 3D viewer. UI and texture callbacks run on the
 * SDL thread. task runs between image/pick jobs on a private snapshot. */
typedef struct {
    int (*height)(const void* state);
    int (*draw)(void*, SDL_Window*, SDL_Renderer*, SDL_Rect);
    int (*event)(void*, SDL_Event*, SDL_Rect);
    void (*other_window_event)(void*, Uint32, const SDL_Event*);
    int (*task)(void*, unsigned char**);
    void (*publish)(void*, const void*, unsigned char*);
    void (*shutdown)(void*);
    void* (*snapshot)(const void*);
} alea_view_panel_t;
typedef struct {
    void (*preview_shown)(void* state, Uint32 elapsed_ms);
    void (*finish)(void* state);
    /* A pick runs on a render-worker snapshot, between frames. Only the UI
     * thread applies a result, and only while its view is still current. */
    int (*pick)(void* snapshot, int x, int y);
    const char* (*picked)(void* state, const void* snapshot, int result);
    const alea_view_panel_t* panel;
    int delay_after_preview; /* Count the idle interval from presentation. */
} alea_view_lifecycle_t;

typedef struct {
    void* snapshot;
    unsigned char* pixels;
    alea_view_draw_t draw;
    alea_view_save_t save;
    SDL_atomic_t done, cancelled, progress;
    SDL_mutex* publication;
    unsigned char* partial;
    uint64_t revision;
    int saving, result, preview;
    int (*pick)(void*, int, int);
    int pick_x, pick_y, panel_task;
} alea_view_job_t;

/* Optional progressive jobs publish owned RGB snapshots through a single-slot
 * mailbox. The UI never reads the worker's mutable framebuffer. */
typedef void (*alea_view_bind_t)(void*, alea_view_job_t*);
static inline int alea_view_cancelled(alea_view_job_t* job) {
    return job && SDL_AtomicGet(&job->cancelled);
}
static inline void alea_view_publish(alea_view_job_t* job, unsigned char* pixels, int progress) {
    SDL_LockMutex(job->publication);
    free(job->partial);job->partial=pixels;
    SDL_AtomicSet(&job->progress,progress);
    SDL_UnlockMutex(job->publication);
}

static Uint32 alea_view_event_window(const SDL_Event* e) {
    switch (e->type) {
    case SDL_WINDOWEVENT: return e->window.windowID;
    case SDL_MOUSEBUTTONDOWN: case SDL_MOUSEBUTTONUP: return e->button.windowID;
    case SDL_MOUSEMOTION: return e->motion.windowID;
    case SDL_MOUSEWHEEL: return e->wheel.windowID;
    case SDL_KEYDOWN: case SDL_KEYUP: return e->key.windowID;
    case SDL_TEXTINPUT: return e->text.windowID;
    default: return 0;
    }
}

static int alea_view_worker(void* opaque) {
    alea_view_job_t* job = opaque;
    if (job->pick) {
        job->result = job->pick(job->snapshot, job->pick_x, job->pick_y);
        SDL_AtomicSet(&job->done, 1);
        return 0;
    }
    job->result = job->draw(job->snapshot, &job->pixels);
    if (job->result == 0 && !job->pixels) job->result = -1;
    if (job->result == 0 && job->saving)
        job->result = job->save(job->snapshot, job->pixels) == 0 ? 0 : -2;
    /* SDL atomics publish the image and result before the UI consumes them. */
    SDL_AtomicSet(&job->done, 1);
    return 0;
}

static void alea_ui_bevel(SDL_Renderer* r, SDL_Rect box, int down) {
    SDL_SetRenderDrawColor(r, down ? 168 : 192, down ? 168 : 192, down ? 168 : 192, 255);
    SDL_RenderFillRect(r, &box);
    for (int i=0;i<2;++i) {
        int light = down ? 96 : 255, dark = down ? 255 : 96;
        SDL_SetRenderDrawColor(r, light, light, light, 255);
        SDL_RenderDrawLine(r, box.x+i, box.y+i, box.x+box.w-i-1, box.y+i);
        SDL_RenderDrawLine(r, box.x+i, box.y+i, box.x+i, box.y+box.h-i-1);
        SDL_SetRenderDrawColor(r, dark, dark, dark, 255);
        SDL_RenderDrawLine(r, box.x+i, box.y+box.h-i-1, box.x+box.w-i-1, box.y+box.h-i-1);
        SDL_RenderDrawLine(r, box.x+box.w-i-1, box.y+i, box.x+box.w-i-1, box.y+box.h-i-1);
    }
}

static int alea_ui_layout(const alea_view_button_t* buttons, int count,
                          int width, SDL_Rect* rects) {
    int x = 8, y = 8;
    for (int i=0;i<count;++i) {
        int w = (int)strlen(buttons[i].label)*12 + 20;
        /* Keep the common file/view actions above the navigation controls. */
        if (i == 6 || x+w > width-8) { x=8; y+=36; }
        rects[i] = (SDL_Rect){x,y,w,30};
        x += w+6;
    }
    return y+38;
}

static int alea_ui_hit(const SDL_Rect* rects, int count, int x, int y) {
    SDL_Point point = {x,y};
    for (int i=0;i<count;++i) if (SDL_PointInRect(&point, &rects[i])) return i;
    return -1;
}

static int alea_sdl_view_progressive(const char* title, const char* help,
                         int width, int height, void* state, size_t state_size,
                         alea_view_draw_t draw, alea_view_save_t save,
                         alea_view_input_t input, alea_view_active_t active,
                         const alea_view_button_t* buttons, int button_count,
                         alea_view_prepare_t prepare, alea_view_bind_t bind,
                         const alea_view_lifecycle_t* lifecycle) {
    SDL_Window* window = NULL;
    SDL_Renderer* renderer = NULL;
    SDL_Texture* texture = NULL;
    SDL_Thread* worker = NULL;
    alea_view_job_t job = {0};
    SDL_Rect rects[32];
    void* pending_save = NULL;
    uint64_t revision = 1, displayed = 0, save_revision = 0, refined = 0;
    const alea_view_panel_t* panel = lifecycle ? lifecycle->panel : NULL;
    void* pending_task = NULL;
    uint64_t pending_pick = 0, image_revision = 0;
    int pick_x = 0, pick_y = 0, click_candidate = 0;
    SDL_Point down = {0,0};
    Uint32 idle_since = 0;
    Uint32 dispatched = 0;
    int rc = -1, closing = 0, have_image = 0, pressed = -1, focus = 0;
    int dragging = 0, show_help = 0;
    const char* notice = "Ready";
    if (width < 1 || height < 1 || width > 16384 || height > 16384 ||
        button_count < 1 || button_count > 32) {
        fprintf(stderr, "Invalid interactive image size or toolbar.\n");
        return -1;
    }
    if (SDL_Init(SDL_INIT_VIDEO) != 0) goto done;
    int window_w = width > 1200 ? 1200 : width < 720 ? 720 : width;
    int window_h = (int)((double)height * window_w / width) + 150;
    if (window_h > 960) window_h = 960;
    if (window_h < 480) window_h = 480;
    window = SDL_CreateWindow(title, SDL_WINDOWPOS_CENTERED,
        SDL_WINDOWPOS_CENTERED, window_w, window_h, SDL_WINDOW_RESIZABLE);
    if (!window) goto done;
    SDL_SetWindowMinimumSize(window, 640, 480);
    renderer = SDL_CreateRenderer(window, -1, SDL_RENDERER_ACCELERATED);
    if (!renderer) renderer = SDL_CreateRenderer(window, -1, SDL_RENDERER_SOFTWARE);
    if (!renderer) goto done;
    texture = SDL_CreateTexture(renderer, SDL_PIXELFORMAT_RGB24,
        SDL_TEXTUREACCESS_STREAMING, width, height);
    if (!texture) goto done;
    fprintf(stderr, "%s\n", help);
    idle_since = SDL_GetTicks();
    for (;;) {
        int current_preview = 0;
        Uint32 preview_ms = 0;
        /* Never join a running job on the event thread. The atomic completion
         * flag is also the ownership handoff for its image and snapshot. */
        if(worker && job.publication) {
            SDL_LockMutex(job.publication);
            unsigned char* partial=job.partial;job.partial=NULL;
            SDL_UnlockMutex(job.publication);
            if(partial && job.revision==revision && !closing) {
                int updated=SDL_UpdateTexture(texture,NULL,partial,width*3);
                free(partial);partial=NULL;
                if(updated!=0)goto done;
                have_image=1; image_revision=job.revision;
            }
            free(partial);
        }
        if (worker && SDL_AtomicGet(&job.done)) {
            SDL_WaitThread(worker, NULL);
            worker = NULL;
            int can_display = job.revision == revision ||
                (!bind && prepare && !job.saving && job.revision >= displayed);
            if (job.panel_task) {
                if (!closing) panel->publish(state, job.snapshot, job.pixels);
            } else if (job.pick) {
                if (!closing && job.revision == revision && !pending_pick)
                    notice = lifecycle->picked(state, job.snapshot, job.result);
            } else if ((job.result == 0 || job.result == -2) && can_display && !closing) {
                if (SDL_UpdateTexture(texture, NULL, job.pixels, width*3) != 0) goto done;
                have_image = 1; image_revision = job.revision;
                displayed = job.revision;
                current_preview = job.preview && job.revision == revision;
                preview_ms = SDL_GetTicks() - dispatched;
                /* Let input following the displayed preview precede refinement. */
                if (current_preview && (!bind || (lifecycle && lifecycle->delay_after_preview)))
                    idle_since = SDL_GetTicks();
                if (!job.preview) refined = job.revision;
            }
            if (!job.panel_task && !job.pick && job.result == -1) {
                notice = "Render failed - change view to retry";
                if (job.revision == revision) displayed = refined = revision;
            } else if (job.saving) notice = job.result == -2 ? "Save failed - check output path" : "Image saved";
            free(job.pixels); job.pixels = NULL;
            free(job.snapshot); job.snapshot = NULL;
            free(job.partial);job.partial=NULL;
            SDL_DestroyMutex(job.publication);job.publication=NULL;
        }
        if (closing && !worker && !pending_save) break;
        int refining = prepare && displayed == revision && refined != revision &&
            SDL_GetTicks()-idle_since >= ALEA_VIEW_REFINE_DELAY_MS;
        if (closing || pending_pick != revision) pending_pick = 0;
        if (closing) { free(pending_task); pending_task = NULL; }
        if (!worker && (pending_task || pending_pick || pending_save || (!closing && (displayed != revision || refining)))) {
            memset(&job, 0, sizeof(job));
            job.saving = pending_save != NULL;
            job.revision = job.saving ? save_revision : revision;
            if (job.saving) {
                job.snapshot = pending_save;
                pending_save = NULL;
            } else if (pending_task && !pending_pick) {
                job.snapshot = pending_task;
                pending_task = NULL;
                job.panel_task = 1;
            } else {
                if (panel) job.snapshot = panel->snapshot(state);
                else {
                    job.snapshot = malloc(state_size);
                    if (job.snapshot) memcpy(job.snapshot, state, state_size);
                }
                if (pending_pick) {
                    job.pick = lifecycle->pick;
                    job.pick_x = pick_x;
                    job.pick_y = pick_y;
                    pending_pick = 0;
                }
            }
            if (!job.snapshot) goto done;
            job.preview = prepare && !job.panel_task && !job.pick && !job.saving && !refining;
            if (prepare) prepare(job.snapshot, job.preview);
            if(bind && !job.pick && !job.panel_task) {
                job.publication=SDL_CreateMutex();if(!job.publication)goto done;
                bind(job.snapshot,&job);
            }
            job.draw = job.panel_task ? panel->task : draw; job.save = save;
            SDL_AtomicSet(&job.done, 0);
            dispatched = SDL_GetTicks();
            worker = SDL_CreateThread(alea_view_worker, "alea-render", &job);
            if (!worker) goto done;
        }
        SDL_GetWindowSize(window, &window_w, &window_h);
        int top = alea_ui_layout(buttons, button_count, window_w, rects);
        SDL_Rect area = {8,top,window_w-16,window_h-top-38};
        int strip_h = panel ? panel->height(state) : 0;
        if (strip_h > area.h / 2) strip_h = area.h / 2;
        SDL_Rect strip = {area.x, area.y + area.h - strip_h, area.w, strip_h};
        if (strip_h) area.h -= strip_h + 4;
        double scale_x = (double)area.w / width, scale_y = (double)area.h / height;
        double scale = scale_x < scale_y ? scale_x : scale_y;
        SDL_Rect image = {0,0,(int)(width*scale),(int)(height*scale)};
        if (image.w < 1) image.w = 1;
        if (image.h < 1) image.h = 1;
        image.x = area.x+(area.w-image.w)/2; image.y = area.y+(area.h-image.h)/2;
        if (!closing) {
            SDL_SetRenderDrawColor(renderer, 192,192,192,255);
            SDL_RenderClear(renderer);
            SDL_SetRenderDrawColor(renderer, 64,64,64,255);
            SDL_RenderFillRect(renderer, &area);
            if (have_image && SDL_RenderCopy(renderer, texture, NULL, &image) != 0) goto done;
            for (int i=0;i<button_count;++i) {
                int down = pressed == i || (active && active(state, buttons[i].key));
                alea_ui_bevel(renderer, rects[i], down);
                SDL_SetRenderDrawColor(renderer, 24,24,24,255);
                alea_ui_text(renderer, rects[i].x+10+down, rects[i].y+8+down, buttons[i].label);
                if (focus == i) {
                    SDL_Rect mark = {rects[i].x+4,rects[i].y+4,rects[i].w-8,rects[i].h-8};
                    SDL_RenderDrawRect(renderer, &mark);
                }
            }
            SDL_Rect status = {8,window_h-30,window_w-16,24};
            alea_ui_bevel(renderer, status, 1);
            SDL_SetRenderDrawColor(renderer, 24,24,24,255);
            const char* message = worker ? (pending_save ? "Rendering - save queued" : job.saving ? "Rendering and saving" : "Rendering - controls available") : notice;
            char progress_text[64];
            if(worker && job.panel_task) message = "Preparing section";
            else if(worker && job.pick) message = "Inspecting surface";
            else if(worker && bind && !job.saving) {
                snprintf(progress_text,sizeof(progress_text),"%s: %d / 100",job.preview?"Preview":"Refining",SDL_AtomicGet(&job.progress));
                message=progress_text;
            }
            alea_ui_text(renderer, 16,window_h-25,message);
            if (worker) {
                SDL_SetRenderDrawColor(renderer, 0,0,128,255);
                SDL_Rect pulse = {window_w-72+(int)((SDL_GetTicks()/100)%4)*12,window_h-23,8,10};
                SDL_RenderFillRect(renderer, &pulse);
            }
            if (panel) {
                int action = panel->draw(state, window, renderer, strip);
                if (action < 0) goto done;
                if (action == 2) {
                    free(pending_task);
                    pending_task = panel->snapshot(state);
                    if (!pending_task) goto done;
                    /* Let the section use the worker before finishing an idle 3D refinement. */
                    if (worker && bind && !job.preview && !job.saving && !job.pick && !job.panel_task)
                        SDL_AtomicSet(&job.cancelled, 1);
                }
                if (action == 1) {
                    ++revision; idle_since = SDL_GetTicks(); notice = "Ready";
                    if (worker && bind && !job.saving) SDL_AtomicSet(&job.cancelled, 1);
                }
            }
            if (show_help) {
                SDL_Rect panel = {area.x+4,area.y+4,area.w-8,160};
                alea_ui_bevel(renderer,panel,0);
                SDL_RenderSetClipRect(renderer,&panel);
                SDL_SetRenderDrawColor(renderer,0,0,128,255);
                alea_ui_text(renderer,panel.x+12,panel.y+12,"Controls - Help closes this panel");
                SDL_SetRenderDrawColor(renderer,24,24,24,255);
                alea_ui_text(renderer,panel.x+12,panel.y+40,help);
                SDL_RenderSetClipRect(renderer,NULL);
            }
            SDL_RenderPresent(renderer);
            if (lifecycle && lifecycle->preview_shown && current_preview) lifecycle->preview_shown(state, preview_ms);
        }
        SDL_Event event;
        if (!SDL_WaitEventTimeout(&event, 16)) continue;
        do {
            if (closing) continue;
            Uint32 target = alea_view_event_window(&event);
            if (target && target != SDL_GetWindowID(window)) {
                if (panel) panel->other_window_event(state, target, &event);
                continue;
            }
            if (panel && !dragging && pressed < 0 && panel->event(state, &event, strip)) continue;
            SDL_Keycode command = 0;
            if (event.type == SDL_QUIT || (event.type == SDL_WINDOWEVENT &&
                event.window.event == SDL_WINDOWEVENT_CLOSE)) command = SDLK_q;
            else if (event.type == SDL_KEYDOWN) {
                command = event.key.keysym.sym;
                if (command == SDLK_TAB) {
                    focus = (focus + ((event.key.keysym.mod & KMOD_SHIFT) ? button_count-1 : 1)) % button_count;
                    continue;
                }
                if (command == SDLK_RETURN || command == SDLK_SPACE) command = buttons[focus].key;
            } else if (event.type == SDL_MOUSEBUTTONDOWN) {
                SDL_Point point = {event.button.x,event.button.y};
                if (event.button.button == SDL_BUTTON_LEFT) {
                    pressed = alea_ui_hit(rects,button_count,point.x,point.y);
                    if (pressed >= 0) focus = pressed;
                }
                down = point;
                dragging = pressed < 0 && !show_help && SDL_PointInRect(&point,&image);
                click_candidate = dragging && lifecycle && lifecycle->pick &&
                    event.button.button == SDL_BUTTON_LEFT && image_revision == revision;
                continue;
            } else if (event.type == SDL_MOUSEBUTTONUP) {
                if (event.button.button == SDL_BUTTON_LEFT && pressed >= 0 &&
                    pressed == alea_ui_hit(rects,button_count,event.button.x,event.button.y))
                    command = buttons[pressed].key;
                SDL_Point point = {event.button.x, event.button.y};
                if (click_candidate && event.button.button == SDL_BUTTON_LEFT && !show_help &&
                    image_revision == revision && SDL_PointInRect(&point, &image) &&
                    abs(point.x-down.x) <= 3 && abs(point.y-down.y) <= 3) {
                    pick_x = (point.x-image.x) * width / image.w;
                    pick_y = (point.y-image.y) * height / image.h;
                    pending_pick = revision;
                    notice = "Inspection queued";
                }
                click_candidate = 0; pressed = -1; dragging = 0;
                if (!command) continue;
            } else if (event.type == SDL_MOUSEMOTION) {
                if (!dragging || show_help) continue;
                if (click_candidate) {
                    if (abs(event.motion.x-down.x) <= 3 && abs(event.motion.y-down.y) <= 3 &&
                        abs(event.motion.xrel) <= 3 && abs(event.motion.yrel) <= 3) continue;
                    click_candidate = 0;
                }
                event.motion.xrel = (int)(event.motion.xrel / scale);
                event.motion.yrel = (int)(event.motion.yrel / scale);
            } else if (event.type == SDL_MOUSEWHEEL) {
                SDL_Point point;
                SDL_GetMouseState(&point.x,&point.y);
                if (show_help || !SDL_PointInRect(&point,&image)) continue;
            } else if (event.type == SDL_WINDOWEVENT && event.window.event == SDL_WINDOWEVENT_FOCUS_LOST) {
                click_candidate = 0; pressed = -1; dragging = 0;
                continue;
            } else continue;
            if (command == SDLK_q || command == SDLK_ESCAPE) {
                closing = 1;
                if (panel) panel->shutdown(state);
                if(worker && bind && !job.saving)SDL_AtomicSet(&job.cancelled,1);
                /* Hide promptly, then pump events while the in-flight library
                 * call finishes. Model/config resources outlive that call. */
                SDL_HideWindow(window);
            } else if (command == SDLK_h) show_help = !show_help;
            else if (command == SDLK_s) {
                if (panel) {
                    free(pending_save);
                    pending_save = panel->snapshot(state);
                } else if (!pending_save) pending_save = malloc(state_size);
                if (!pending_save) goto done;
                if (!panel) memcpy(pending_save,state,state_size);
                save_revision = revision;
            } else {
                if (command) { event.type = SDL_KEYDOWN; event.key.keysym.sym = command; }
                int changed = input(state,&event);
                if (changed < 0) goto done;
                if (changed) {
                    ++revision; idle_since = SDL_GetTicks(); notice = "Ready";
                    if(worker && bind && !job.saving)SDL_AtomicSet(&job.cancelled,1);
                }
            }
        } while (SDL_PollEvent(&event));
    }
    rc = 0;
done:
    /* Exceptional SDL/allocation failures still cannot free a running job's
     * resources. Hide the failed window and continue pumping until it exits. */
    if (rc != 0) fprintf(stderr, "Viewer failed: %s\n", SDL_GetError());
    if (worker) {
        if(bind && !job.saving)SDL_AtomicSet(&job.cancelled,1);
        if (window) SDL_HideWindow(window);
        while (!SDL_AtomicGet(&job.done)) { SDL_PumpEvents(); SDL_Delay(10); }
        SDL_WaitThread(worker,NULL);
    }
    free(job.pixels); free(job.snapshot); free(pending_save); free(pending_task); free(job.partial);
    SDL_DestroyMutex(job.publication);
    if (panel && !closing) panel->shutdown(state);
    SDL_DestroyTexture(texture);
    SDL_DestroyRenderer(renderer);
    SDL_DestroyWindow(window);
    if (lifecycle && lifecycle->finish) lifecycle->finish(state);
    SDL_Quit();
    return rc;
}
static inline int alea_sdl_view_ex(const char* title, const char* help,
        int width, int height, void* state, size_t state_size,
        alea_view_draw_t draw, alea_view_save_t save,
        alea_view_input_t input, alea_view_active_t active,
        const alea_view_button_t* buttons, int button_count, alea_view_prepare_t prepare) {
    return alea_sdl_view_progressive(title,help,width,height,state,state_size,draw,save,
        input,active,buttons,button_count,prepare,NULL,NULL);
}
/* Preserve the full-resolution behavior of existing viewers. */
static inline int alea_sdl_view(const char* title, const char* help,
        int width, int height, void* state, size_t state_size,
        alea_view_draw_t draw, alea_view_save_t save,
        alea_view_input_t input, alea_view_active_t active,
        const alea_view_button_t* buttons, int button_count) {
    return alea_sdl_view_ex(title,help,width,height,state,state_size,draw,save,
        input,active,buttons,button_count,NULL);
}
#endif
