// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
// SPDX-License-Identifier: MPL-2.0

/* Shared shell regression: a slow renderer, actual mouse events, rapid edits,
 * save snapshots, repaints/resizing during rendering, and close during a job. */
#include <assert.h>
#include <SDL.h>
static SDL_atomic_t painting, draws, saves, started;
static int busy_presents, busy_inputs, saved_value;
static SDL_threadID ui_thread;
static SDL_Window* test_window;
static void test_present(SDL_Renderer* renderer) {
    assert(SDL_ThreadID() == ui_thread);
    test_window = SDL_RenderGetWindow(renderer);
    if (SDL_AtomicGet(&painting)) ++busy_presents;
    SDL_RenderPresent(renderer);
}
#define SDL_RenderPresent test_present
#include "../../tools/sdl_view.h"
#undef SDL_RenderPresent

typedef struct { int value; } state_t;
static int slow_draw(void* opaque, unsigned char** pixels) {
    state_t* state = opaque;
    int value = state->value;
    assert(SDL_ThreadID() != ui_thread);
    assert(SDL_AtomicAdd(&painting,1) == 0);
    SDL_AtomicSet(&started,1);
    SDL_Delay(300);
    assert(state->value == value); /* UI edits never touch this snapshot. */
    *pixels = malloc(64*64*3);
    assert(*pixels);
    memset(*pixels,value,64*64*3);
    SDL_AtomicAdd(&draws,1);
    SDL_AtomicSet(&painting,0);
    return 0;
}
static int save_frame(const void* opaque, const unsigned char* pixels) {
    assert(SDL_ThreadID() != ui_thread);
    saved_value = ((const state_t*)opaque)->value;
    assert(pixels[0] == saved_value);
    SDL_AtomicAdd(&saves,1);
    return 0;
}
static int change_view(void* opaque, const SDL_Event* event) {
    assert(SDL_ThreadID() == ui_thread);
    if (SDL_AtomicGet(&painting)) ++busy_inputs;
    if (event->type == SDL_KEYDOWN && event->key.keysym.sym == SDLK_PLUS) {
        ++((state_t*)opaque)->value;
        return 1;
    }
    if (event->type == SDL_KEYDOWN && event->key.keysym.sym == SDLK_v) {
        assert(test_window);
        SDL_SetWindowSize(test_window,900,600);
    }
    return 0;
}
static void key(SDL_Keycode code) {
    SDL_Event event = {0};
    event.type = SDL_KEYDOWN;
    event.key.keysym.sym = code;
    assert(SDL_PushEvent(&event) == 1);
}
static const alea_view_button_t buttons[] = {
    {"Save", SDLK_s}, {"Reset", SDLK_r}, {"Zoom +", SDLK_PLUS},
    {"Zoom -", SDLK_MINUS}, {"Help", SDLK_h}, {"Quit", SDLK_q}
};
static void click(int index) {
    SDL_Rect rects[6];
    alea_ui_layout(buttons,6,720,rects);
    SDL_Event event = {0};
    event.type = SDL_MOUSEBUTTONDOWN;
    event.button.button = SDL_BUTTON_LEFT;
    event.button.x = rects[index].x+10;
    event.button.y = rects[index].y+10;
    assert(SDL_PushEvent(&event) == 1);
    event.type = SDL_MOUSEBUTTONUP;
    assert(SDL_PushEvent(&event) == 1);
}
static int feed_events(void* opaque) {
    int close_early = *(int*)opaque;
    while (!SDL_AtomicGet(&started)) SDL_Delay(1);
    if (close_early) { click(5); return 0; }
    click(2); /* value 2, while initial value 1 is rendering */
    click(0); /* save a snapshot of 2 */
    key(SDLK_PLUS); /* live view becomes 3 before save executes */
    key(SDLK_v); /* resize from the SDL thread */
    while (SDL_AtomicGet(&draws) < 3 || SDL_AtomicGet(&saves) < 1) SDL_Delay(1);
    SDL_Delay(100); /* Let the UI consume the completed frame. */
    key(SDLK_q);
    return 0;
}
int main(int argc, char** argv) {
    (void)argv;
    int close_early = argc > 1;
    ui_thread = SDL_ThreadID();
    assert(SDL_Init(SDL_INIT_VIDEO) == 0);
    SDL_Thread* feeder = SDL_CreateThread(feed_events,"test-input",&close_early);
    assert(feeder);
    state_t state = {1};
    assert(alea_sdl_view("Test", "Test controls",64,64,&state,sizeof(state),
        slow_draw,save_frame,change_view,NULL,buttons,6) == 0);
    SDL_WaitThread(feeder,NULL);
    assert(!SDL_AtomicGet(&painting));
    if (close_early) assert(SDL_AtomicGet(&draws) == 1);
    else {
        assert(state.value == 3 && saved_value == 2);
        assert(SDL_AtomicGet(&draws) == 3);
        assert(busy_inputs >= 2 && busy_presents >= 3);
    }
    puts(close_early ? "close during render: PASS" : "responsive controls and save snapshots: PASS");
    return 0;
}
