// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
// SPDX-License-Identifier: MPL-2.0

/* Shared shell regression: a slow renderer, actual mouse events, rapid edits,
 * save snapshots, repaints/resizing during rendering, and close during a job. */
#include <assert.h>
#include <SDL.h>
static SDL_atomic_t painting, draws, saves, started, previews, full_frames, latest_full;
static int busy_presents, busy_inputs, saved_value;
static SDL_atomic_t cancellations, publications, displayed_partial;
static int progressive_mode, live_value=1;
static SDL_threadID ui_thread;
static SDL_Window* test_window;
static void test_present(SDL_Renderer* renderer) {
    assert(SDL_ThreadID() == ui_thread);
    test_window = SDL_RenderGetWindow(renderer);
    if (SDL_AtomicGet(&painting)) ++busy_presents;
    SDL_RenderPresent(renderer);
}
static int test_update(SDL_Texture* texture,const SDL_Rect* rect,const void* pixels,int pitch) {
    assert(SDL_ThreadID()==ui_thread);
    if(progressive_mode){assert(((const unsigned char*)pixels)[0]==live_value);SDL_AtomicAdd(&displayed_partial,1);}
    return SDL_UpdateTexture(texture,rect,pixels,pitch);
}
#define SDL_UpdateTexture test_update
#define SDL_RenderPresent test_present
#include "../../tools/sdl_view.h"
#undef SDL_RenderPresent
#undef SDL_UpdateTexture

typedef struct { int value, preview; alea_view_job_t* job; } state_t;
static void bind_job(void* state,alea_view_job_t* job){((state_t*)state)->job=job;}
static void prepare_quality(void* state, int preview) {
    assert(SDL_ThreadID() == ui_thread);
    ((state_t*)state)->preview=preview;
}
static int slow_draw(void* opaque, unsigned char** pixels) {
    state_t* state = opaque;
    int value = state->value;
    assert(SDL_ThreadID() != ui_thread);
    assert(SDL_AtomicAdd(&painting,1) == 0);
    SDL_AtomicSet(&started,1);
    for(int i=0;i<100;++i) {
        SDL_Delay(3);
        if(state->job && alea_view_cancelled(state->job)) {
            SDL_AtomicAdd(&cancellations,1);SDL_AtomicSet(&painting,0);return 1;
        }
        if(state->job && !state->job->saving && i==10) {
            unsigned char* partial=malloc(64*64*3);assert(partial);memset(partial,value,64*64*3);
            alea_view_publish(state->job,partial,10);SDL_AtomicAdd(&publications,1);
        }
    }
    assert(state->value == value); /* UI edits never touch this snapshot. */
    *pixels = malloc(64*64*3);
    assert(*pixels);
    memset(*pixels,value,64*64*3);
    if (state->preview) SDL_AtomicAdd(&previews,1);
    else { SDL_AtomicAdd(&full_frames,1); SDL_AtomicSet(&latest_full,state->value); }
    SDL_AtomicAdd(&draws,1);
    SDL_AtomicSet(&painting,0);
    return 0;
}
static int save_frame(const void* opaque, const unsigned char* pixels) {
    assert(SDL_ThreadID() != ui_thread);
    assert(((const state_t*)opaque)->preview==0);
    saved_value = ((const state_t*)opaque)->value;
    assert(pixels[0] == saved_value);
    SDL_AtomicAdd(&saves,1);
    return 0;
}
static int change_view(void* opaque, const SDL_Event* event) {
    assert(SDL_ThreadID() == ui_thread);
    if (SDL_AtomicGet(&painting)) ++busy_inputs;
    if (event->type == SDL_KEYDOWN && event->key.keysym.sym == SDLK_PLUS) {
        ++((state_t*)opaque)->value;live_value=((state_t*)opaque)->value;
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
    int mode = *(int*)opaque;
    while (!SDL_AtomicGet(&started)) SDL_Delay(1);
    if (mode==1) { click(5); return 0; }
    if(mode==4){while(!SDL_AtomicGet(&displayed_partial))SDL_Delay(1);click(5);return 0;}
    if(mode==3)while(!SDL_AtomicGet(&displayed_partial))SDL_Delay(1);
    if (mode==2) while (SDL_AtomicGet(&draws)<2) SDL_Delay(1); /* Initial preview then idle refinement. */
    click(2); /* value 2, while initial value 1 is rendering */
    click(0); /* save a snapshot of 2 */
    key(SDLK_PLUS); /* live view becomes 3 before save executes */
    key(SDLK_v); /* resize from the SDL thread */
    while (SDL_AtomicGet(&draws) < 3 || SDL_AtomicGet(&saves) < 1) SDL_Delay(1);
    if (mode>=2) while (SDL_AtomicGet(&latest_full)!=3) SDL_Delay(1);
    SDL_Delay(100); /* Let the UI consume the completed frame. */
    key(SDLK_q);
    return 0;
}
int main(int argc, char** argv) {
    int mode = argc<2 ? 0 : strcmp(argv[1],"adaptive")==0 ? 2 : strcmp(argv[1],"progressive")==0 ? 3 : strcmp(argv[1],"progressive-close")==0 ? 4 : 1;
    progressive_mode=mode>=3;
    ui_thread = SDL_ThreadID();
    assert(SDL_Init(SDL_INIT_VIDEO) == 0);
    SDL_Thread* feeder = SDL_CreateThread(feed_events,"test-input",&mode);
    assert(feeder);
    state_t state = {.value=1};
    int rc = mode>=3 ? alea_sdl_view_progressive("Test","Test controls",64,64,&state,sizeof(state),
        slow_draw,save_frame,change_view,NULL,buttons,6,prepare_quality,bind_job) : mode==2 ? alea_sdl_view_ex("Test","Test controls",64,64,&state,sizeof(state),
        slow_draw,save_frame,change_view,NULL,buttons,6,prepare_quality) :
        alea_sdl_view("Test","Test controls",64,64,&state,sizeof(state),
        slow_draw,save_frame,change_view,NULL,buttons,6);
    assert(rc==0);
    SDL_WaitThread(feeder,NULL);
    assert(!SDL_AtomicGet(&painting));
    if(mode==4){assert(SDL_AtomicGet(&cancellations)==1);assert(SDL_AtomicGet(&draws)==0);}
    else if (mode==1) assert(SDL_AtomicGet(&draws) == 1);
    else {
        assert(state.value == 3 && saved_value == 2);
        if (mode==0) assert(SDL_AtomicGet(&draws) == 3);
        else if(mode==3){assert(SDL_AtomicGet(&cancellations)>=1);assert(SDL_AtomicGet(&publications)>=2);assert(SDL_AtomicGet(&displayed_partial)>=2);}
        else { assert(SDL_AtomicGet(&previews)>=1); assert(SDL_AtomicGet(&full_frames)>=3); }
        assert(busy_presents >= 3);
        if (mode==0) assert(busy_inputs >= 2);
    }
    puts(mode==4 ? "close cancels progressive work: PASS" : mode==3 ? "progressive publication, cancellation, stale frame rejection, and save snapshots: PASS" : mode==2 ? "preview, idle refinement, and full-quality save snapshots: PASS" : mode==1 ? "close during render: PASS" : "responsive controls and save snapshots: PASS");
    return 0;
}
