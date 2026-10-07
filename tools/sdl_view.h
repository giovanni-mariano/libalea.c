// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
// SPDX-License-Identifier: MPL-2.0

#ifndef ALEA_SDL_VIEW_H
#define ALEA_SDL_VIEW_H

#include <SDL.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Only the SDL thread touches live state, windows, textures, or events. A job
 * owns a byte-copy of state and a newly allocated RGB image. Pointers inside a
 * snapshot must refer to immutable data kept alive until the worker joins.
 * One worker at a time also serializes access to the model's query caches. */
typedef int (*alea_view_draw_t)(void*, unsigned char**);
typedef int (*alea_view_save_t)(const void*, const unsigned char*);
typedef int (*alea_view_input_t)(void*, const SDL_Event*);
typedef int (*alea_view_active_t)(const void*, SDL_Keycode);
typedef struct { const char* label; SDL_Keycode key; } alea_view_button_t;

typedef struct {
    void* snapshot;
    unsigned char* pixels;
    alea_view_draw_t draw;
    alea_view_save_t save;
    SDL_atomic_t done;
    uint64_t revision;
    int saving, result;
} alea_view_job_t;

static int alea_view_worker(void* opaque) {
    alea_view_job_t* job = opaque;
    job->result = job->draw(job->snapshot, &job->pixels);
    if (job->result == 0 && !job->pixels) job->result = -1;
    if (job->result == 0 && job->saving)
        job->result = job->save(job->snapshot, job->pixels) == 0 ? 0 : -2;
    /* SDL atomics publish the image and result before the UI consumes them. */
    SDL_AtomicSet(&job->done, 1);
    return 0;
}

/* Small built-in bitmap alphabet: no font installation or SDL_ttf required. */
static void alea_ui_text(SDL_Renderer* renderer, int x, int y, const char* text) {
    static const unsigned char letters[][7] = {
        {14,17,17,31,17,17,17},{30,17,17,30,17,17,30},
        {14,17,16,16,16,17,14},{30,17,17,17,17,17,30},
        {31,16,16,30,16,16,31},{31,16,16,30,16,16,16},
        {14,17,16,23,17,17,15},{17,17,17,31,17,17,17},
        {14,4,4,4,4,4,14},{7,2,2,2,18,18,12},
        {17,18,20,24,20,18,17},{16,16,16,16,16,16,31},
        {17,27,21,21,17,17,17},{17,25,25,21,19,19,17},
        {14,17,17,17,17,17,14},{30,17,17,30,16,16,16},
        {14,17,17,17,21,18,13},{30,17,17,30,20,18,17},
        {15,16,16,14,1,1,30},{31,4,4,4,4,4,4},
        {17,17,17,17,17,17,14},{17,17,17,17,17,10,4},
        {17,17,17,21,21,21,10},{17,17,10,4,10,17,17},
        {17,17,10,4,4,4,4},{31,1,2,4,8,16,31},
        {0,0,14,1,15,17,15},{16,16,30,17,17,17,30},
        {0,0,14,17,16,17,14},{1,1,15,17,17,17,15},
        {0,0,14,17,31,16,14},{6,8,8,28,8,8,8},
        {0,0,15,17,15,1,14},{16,16,30,17,17,17,17},
        {4,0,12,4,4,4,14},{2,0,6,2,2,18,12},
        {16,16,18,20,24,20,18},{12,4,4,4,4,4,14},
        {0,0,26,21,21,21,21},{0,0,30,17,17,17,17},
        {0,0,14,17,17,17,14},{0,0,30,17,30,16,16},
        {0,0,15,17,15,1,1},{0,0,22,25,16,16,16},
        {0,0,15,16,14,1,30},{8,8,28,8,8,9,6},
        {0,0,17,17,17,19,13},{0,0,17,17,17,10,4},
        {0,0,17,17,21,21,10},{0,0,17,10,4,10,17},
        {0,0,17,17,15,1,14},{0,0,31,2,4,8,31}
    };
    static const unsigned char digits[][7] = {
        {14,17,19,21,25,17,14},{4,12,4,4,4,4,14},
        {14,17,1,2,4,8,31},{30,1,1,14,1,1,30},
        {2,6,10,18,31,2,2},{31,16,16,30,1,1,30},
        {14,16,16,30,17,17,14},{31,1,2,4,8,8,8},
        {14,17,17,14,17,17,14},{14,17,17,15,1,1,14}
    };
    int start = x;
    for (; *text; ++text, x += 12) {
        unsigned char glyph[7] = {0};
        if (*text == '\n') { y += 24; x = start - 12; continue; }
        if (*text >= 'A' && *text <= 'Z') memcpy(glyph, letters[*text-'A'], 7);
        else if (*text >= 'a' && *text <= 'z') memcpy(glyph, letters[26+*text-'a'], 7);
        else if (*text >= '0' && *text <= '9') memcpy(glyph, digits[*text-'0'], 7);
        else if (*text == '-') glyph[3] = 14;
        else if (*text == '+') { glyph[1]=4; glyph[2]=4; glyph[3]=31; glyph[4]=4; glyph[5]=4; }
        else if (*text == '.') glyph[6] = 4;
        else if (*text == ':') { glyph[2]=4; glyph[5]=4; }
        else if (*text == '/') { for (int i=0;i<7;++i) glyph[i]=(unsigned char)(1 << (i*4/6)); }
        else if (*text == '|') { for (int i=0;i<7;++i) glyph[i]=4; }
        for (int row = 0; row < 7; ++row)
            for (int col = 0; col < 5; ++col)
                if (glyph[row] & (16 >> col)) {
                    SDL_Rect pixel = {x+col*2, y+row*2, 2, 2};
                    SDL_RenderFillRect(renderer, &pixel);
                }
    }
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

static int alea_sdl_view(const char* title, const char* help,
                         int width, int height, void* state, size_t state_size,
                         alea_view_draw_t draw, alea_view_save_t save,
                         alea_view_input_t input, alea_view_active_t active,
                         const alea_view_button_t* buttons, int button_count) {
    SDL_Window* window = NULL;
    SDL_Renderer* renderer = NULL;
    SDL_Texture* texture = NULL;
    SDL_Thread* worker = NULL;
    alea_view_job_t job = {0};
    SDL_Rect rects[32];
    void* pending_save = NULL;
    uint64_t revision = 1, displayed = 0, save_revision = 0;
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
    for (;;) {
        /* Never join a running job on the event thread. The atomic completion
         * flag is also the ownership handoff for its image and snapshot. */
        if (worker && SDL_AtomicGet(&job.done)) {
            SDL_WaitThread(worker, NULL);
            worker = NULL;
            if (job.result != -1 && job.revision == revision && !closing) {
                if (SDL_UpdateTexture(texture, NULL, job.pixels, width*3) != 0) goto done;
                have_image = 1;
                displayed = revision;
            }
            if (job.result == -1) {
                notice = "Render failed - change view to retry";
                if (job.revision == revision) displayed = revision;
            } else if (job.saving) notice = job.result == -2 ? "Save failed - check output path" : "Image saved";
            free(job.pixels); job.pixels = NULL;
            free(job.snapshot); job.snapshot = NULL;
        }
        if (closing && !worker && !pending_save) break;
        if (!worker && (pending_save || (!closing && displayed != revision))) {
            memset(&job, 0, sizeof(job));
            job.snapshot = pending_save ? pending_save : malloc(state_size);
            if (!job.snapshot) goto done;
            job.saving = pending_save != NULL;
            job.revision = pending_save ? save_revision : revision;
            if (!pending_save) memcpy(job.snapshot, state, state_size);
            pending_save = NULL;
            job.draw = draw; job.save = save;
            SDL_AtomicSet(&job.done, 0);
            worker = SDL_CreateThread(alea_view_worker, "alea-render", &job);
            if (!worker) goto done;
        }
        SDL_GetWindowSize(window, &window_w, &window_h);
        int top = alea_ui_layout(buttons, button_count, window_w, rects);
        SDL_Rect area = {8,top,window_w-16,window_h-top-38};
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
            alea_ui_text(renderer, 16,window_h-25,message);
            if (worker) {
                SDL_SetRenderDrawColor(renderer, 0,0,128,255);
                SDL_Rect pulse = {window_w-72+(int)((SDL_GetTicks()/100)%4)*12,window_h-23,8,10};
                SDL_RenderFillRect(renderer, &pulse);
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
        }
        SDL_Event event;
        if (!SDL_WaitEventTimeout(&event, 16)) continue;
        do {
            if (closing) continue;
            SDL_Keycode command = 0;
            if (event.type == SDL_QUIT) command = SDLK_q;
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
                dragging = pressed < 0 && !show_help && SDL_PointInRect(&point,&image);
                continue;
            } else if (event.type == SDL_MOUSEBUTTONUP) {
                if (event.button.button == SDL_BUTTON_LEFT && pressed >= 0 &&
                    pressed == alea_ui_hit(rects,button_count,event.button.x,event.button.y))
                    command = buttons[pressed].key;
                pressed = -1; dragging = 0;
                if (!command) continue;
            } else if (event.type == SDL_MOUSEMOTION) {
                if (!dragging || show_help) continue;
                event.motion.xrel = (int)(event.motion.xrel / scale);
                event.motion.yrel = (int)(event.motion.yrel / scale);
            } else if (event.type == SDL_MOUSEWHEEL) {
                SDL_Point point;
                SDL_GetMouseState(&point.x,&point.y);
                if (show_help || !SDL_PointInRect(&point,&image)) continue;
            } else if (event.type == SDL_WINDOWEVENT && event.window.event == SDL_WINDOWEVENT_FOCUS_LOST) {
                pressed = -1; dragging = 0;
                continue;
            } else continue;
            if (command == SDLK_q || command == SDLK_ESCAPE) {
                closing = 1;
                /* Hide promptly, then pump events while the in-flight library
                 * call finishes. Model/config resources outlive that call. */
                SDL_HideWindow(window);
            } else if (command == SDLK_h) show_help = !show_help;
            else if (command == SDLK_s) {
                if (!pending_save) pending_save = malloc(state_size);
                if (!pending_save) goto done;
                memcpy(pending_save,state,state_size);
                save_revision = revision;
            } else {
                if (command) { event.type = SDL_KEYDOWN; event.key.keysym.sym = command; }
                int changed = input(state,&event);
                if (changed < 0) goto done;
                if (changed) { ++revision; notice = "Ready"; }
            }
        } while (SDL_PollEvent(&event));
    }
    rc = 0;
done:
    /* Exceptional SDL/allocation failures still cannot free a running job's
     * resources. Hide the failed window and continue pumping until it exits. */
    if (rc != 0) fprintf(stderr, "Viewer failed: %s\n", SDL_GetError());
    if (worker) {
        if (window) SDL_HideWindow(window);
        while (!SDL_AtomicGet(&job.done)) { SDL_PumpEvents(); SDL_Delay(10); }
        SDL_WaitThread(worker,NULL);
    }
    free(job.pixels); free(job.snapshot); free(pending_save);
    SDL_DestroyTexture(texture);
    SDL_DestroyRenderer(renderer);
    SDL_DestroyWindow(window);
    SDL_Quit();
    return rc;
}
#endif
