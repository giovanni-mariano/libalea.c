// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
// SPDX-License-Identifier: MPL-2.0

#ifndef ALEA_RENDER_VIEW_H
#define ALEA_RENDER_VIEW_H
#include <stdatomic.h>
#include "sdl_view.h"
#include "render_inspector.h"

typedef struct {
    alea_system_t* sys;
    alea_bbox_t domain;
    alea_prepared_queries_t* index;
    SDL_Thread* thread;
    atomic_int finished;
    int attempted;
} render_view_build_t;

typedef struct {
    render_framebuffer_t* frames[2]; /* One job at a time, tile workers own disjoint rectangles. */
    render_framebuffer_t* display[2];
    uint64_t preview_revision;
    int preview_edge;
} render_view_buffers_t;

typedef struct {
    alea_system_t* sys;
    alea_compute_t* compute;
    alea_prepared_queries_t* prepared;
    render_view_build_t* build;
    render_view_buffers_t* buffers;
    int preview, progressive, progress_error;
    alea_view_job_t* job;
    Uint32 last_publish;
    size_t published_pixels;
    render_config_t cfg, initial, working;
    render_camera_t cam;
    render_pick_t pick;
    int has_pivot;
    inspector_t* inspector; /* UI thread only; jobs own their filter IDs. */
    section_view_t section;
    char inspection[160];
    const char* output;
} render_viewer_t;

static int render_view_build_index(void* opaque) {
    render_view_build_t* build = opaque;
    SDL_SetThreadPriority(SDL_THREAD_PRIORITY_LOW);
    /* The viewer finished preparing hierarchy caches before opening. Build a
     * private handle on its read-only model and publish only the finished index. */
    build->index = alea_prepared_queries_create(build->sys, &build->domain);
    if (!build->index && !alea_interrupted())
        fprintf(stderr, "Regional preparation failed, continuing with hierarchy: %s\n", alea_error());
    atomic_store(&build->finished, 1);
    return 0;
}

static void render_view_preview_shown(void* state, Uint32 elapsed_ms) {
    render_view_build_t* build = ((render_viewer_t*)state)->build;
    /* A preview completed inside the viewer's refinement delay
     * already serves short interactions cheaply. Reconsider on later views. */
    if (build->attempted || elapsed_ms < ALEA_VIEW_REFINE_DELAY_MS) return;
    build->attempted = 1;
    build->thread = SDL_CreateThread(render_view_build_index, "alea-regional", build);
    if (!build->thread)
        fprintf(stderr, "Cannot start regional preparation: %s\n", SDL_GetError());
}

static void render_view_finish_build(void* state) {
    render_view_build_t* build = ((render_viewer_t*)state)->build;
    if (!build->thread) return;
    /* The standalone viewer has drained all render and save jobs by now.
     * The builder is the only remaining library operation in this process. */
    int cancel = !atomic_load(&build->finished) && !alea_interrupted();
    if (cancel) alea_interrupt();
    SDL_WaitThread(build->thread, NULL);
    if (cancel) alea_clear_interrupt();
    alea_prepared_queries_destroy(build->index);
}

/* Seed a refinement from the completed preview for exactly this camera. */
static void render_view_seed(render_viewer_t* v, render_framebuffer_t* fb) {
    const render_framebuffer_t* seed = !v->preview && v->buffers->preview_revision == v->job->revision
        ? v->buffers->frames[1] : NULL;
    for (int y = 0; y < fb->height; ++y) for (int x = 0; x < fb->width; ++x) {
        size_t i = (size_t)y * fb->width + x;
        size_t j = seed ? ((size_t)y * seed->height / fb->height) * seed->width +
                          (size_t)x * seed->width / fb->width : 0;
        fb->cell_id[i] = seed ? seed->cell_id[j] : -1;
        fb->material_id[i] = seed ? seed->material_id[j] : 0;
        fb->depth[i] = seed ? seed->depth[j] : INFINITY;
        for (int a = 0; a < 3; ++a) {
            fb->color[3*i+a] = seed ? seed->color[3*j+a] : v->cfg.background[a];
            fb->normal[3*i+a] = seed ? seed->normal[3*j+a] : 0;
        }
    }
}

static void render_view_copy(render_framebuffer_t* dst, const render_framebuffer_t* src,
                             size_t offset, size_t count) {
    memcpy(dst->color + 3*offset, src->color + 3*offset, count*3*sizeof(float));
    memcpy(dst->normal + 3*offset, src->normal + 3*offset, count*3*sizeof(float));
    memcpy(dst->depth + offset, src->depth + offset, count*sizeof(float));
    memcpy(dst->cell_id + offset, src->cell_id + offset, count*sizeof(int));
    memcpy(dst->material_id + offset, src->material_id + offset, count*sizeof(int));
}

/* Postprocessing is destructive. Keep raw progressive pixels separate so
 * repeatedly displaying a partial image cannot repeatedly darken its edges. */
static uint8_t* render_view_image(render_viewer_t* v, const render_config_t* cfg,
                                 const render_framebuffer_t* raw, const uint8_t* ready) {
    render_framebuffer_t** slot=&v->buffers->display[v->preview];
    if(*slot && ((*slot)->width!=raw->width || (*slot)->height!=raw->height)) {
        render_framebuffer_free(*slot);*slot=NULL;
    }
    if(!*slot)*slot=render_framebuffer_create(raw->width,raw->height,1);
    render_framebuffer_t* fb=*slot;if(!fb)return NULL;
    size_t n=(size_t)raw->width*raw->height;
    if (ready) {
        /* Other workers may still be writing any tile that is not marked ready. */
        render_view_seed(v, fb);
        int tile = cfg->tile_size;
        int cols = (raw->width + tile - 1) / tile;
        for (int y = 0; y < raw->height; ++y) {
            for (int x = 0; x < raw->width; x += tile) {
                if (ready[(y / tile) * cols + x / tile]) {
                    int width = raw->width - x < tile ? raw->width - x : tile;
                    render_view_copy(fb, raw, (size_t)y * raw->width + x, width);
                }
            }
        }
    } else {
        render_view_copy(fb, raw, 0, n);
    }
    postprocess_frame(cfg,fb);
    uint8_t* rgb=malloc((size_t)v->cfg.width*v->cfg.height*3);
    if(!rgb)return NULL;
    if(raw->width==v->cfg.width && raw->height==v->cfg.height)render_tonemap(fb,rgb);
    else {
        uint8_t* small=malloc(n*3);if(!small){free(rgb);return NULL;}
        render_tonemap(fb,small);
        for(int y=0;y<v->cfg.height;++y)for(int x=0;x<v->cfg.width;++x) {
            size_t src=((size_t)y*raw->height/v->cfg.height*raw->width+(size_t)x*raw->width/v->cfg.width)*3;
            memcpy(rgb+((size_t)y*v->cfg.width+x)*3,small+src,3);
        }
        free(small);
    }
    return rgb;
}
#ifdef ALEA_VIEW_METAL
static int render_view_progress(void* user,const render_framebuffer_t* fb,size_t completed,size_t total) {
    render_viewer_t* v=user;
    if(alea_view_cancelled(v->job))return 1;
    SDL_AtomicSet(&v->job->progress,(int)(100*completed/total));
    if(completed==v->published_pixels || completed==total || SDL_GetTicks()-v->last_publish<50)return 0;
    uint8_t* rgb=render_view_image(v,&v->working,fb,NULL);
    if(!rgb){v->progress_error=1;return 1;}
    alea_view_publish(v->job,rgb,(int)(100*completed/total));
    v->last_publish=SDL_GetTicks();v->published_pixels=completed;
    return alea_view_cancelled(v->job);
}
#endif
typedef struct {
    render_viewer_t* view;
    SDL_mutex* lock;
    SDL_atomic_t failed;
    uint8_t* ready;
    size_t completed;
    Uint32 last_publish;
} render_view_cpu_t;

static int render_view_cpu_cancelled(void* user) {
    render_view_cpu_t* cpu = user;
    return alea_view_cancelled(cpu->view->job) || SDL_AtomicGet(&cpu->failed);
}

static void render_view_tile_done(void* user, const render_framebuffer_t* fb,
                                  int x, int y, int width, int height) {
    render_view_cpu_t* cpu = user;
    render_viewer_t* v = cpu->view;
    int tile = v->working.tile_size;
    int cols = (fb->width + tile - 1) / tile;
    size_t total = (size_t)fb->width * fb->height;
    SDL_LockMutex(cpu->lock);
    cpu->ready[(y / tile) * cols + x / tile] = 1;
    cpu->completed += (size_t)width * height;
    int percent = (int)(100.0 * cpu->completed / total);
    SDL_AtomicSet(&v->job->progress, percent);
    if (cpu->completed < total && SDL_GetTicks() - cpu->last_publish >= 50 &&
        !render_view_cpu_cancelled(cpu)) {
        uint8_t* pixels = render_view_image(v, &v->working, fb, cpu->ready);
        if (pixels) alea_view_publish(v->job, pixels, percent);
        else SDL_AtomicSet(&cpu->failed, 1);
        cpu->last_publish = SDL_GetTicks();
    }
    SDL_UnlockMutex(cpu->lock);
}

static int render_view_cpu_frame(render_viewer_t* v, alea_prepared_queries_t* index,
                                 render_framebuffer_t* fb) {
    int tile = v->working.tile_size;
    size_t cols = (fb->width + tile - 1) / tile;
    size_t rows = (fb->height + tile - 1) / tile;
    render_view_cpu_t cpu = {.view = v, .lock = SDL_CreateMutex(),
                             .ready = calloc(cols * rows, 1)};
    const render_updates_t updates = {&cpu, render_view_cpu_cancelled, render_view_tile_done};
    int rc = -1;
    if (cpu.lock && cpu.ready) {
        cpu.last_publish = SDL_GetTicks();
        rc = alea_render_scene_progressive(v->sys, index, &v->working, &v->cam, fb, &updates);
    }
    SDL_DestroyMutex(cpu.lock);
    free(cpu.ready);
    return SDL_AtomicGet(&cpu.failed) ? -1 : rc;
}

static int render_view_draw(void* state, uint8_t** pixels) {
    render_viewer_t* v=state;*pixels=NULL;
    render_config_t cfg=v->cfg;
    if(v->preview) {
        if(v->progressive) {
            int edge=v->buffers->preview_edge?v->buffers->preview_edge:96;
            int longest=cfg.width>cfg.height?cfg.width:cfg.height;
            if(longest>edge) {
                cfg.width=(int)((double)cfg.width*edge/longest);if(cfg.width<1)cfg.width=1;
                cfg.height=(int)((double)cfg.height*edge/longest);if(cfg.height<1)cfg.height=1;
            }
        } else {cfg.width=(cfg.width+1)/2;cfg.height=(cfg.height+1)/2;}
        cfg.aa_samples=1;cfg.shadows=0;
    }
    render_framebuffer_t** slot=&v->buffers->frames[v->preview];
    if(*slot && ((*slot)->width!=cfg.width || (*slot)->height!=cfg.height)) {
        render_framebuffer_free(*slot);*slot=NULL;
    }
    if(!*slot)*slot=render_framebuffer_create(cfg.width,cfg.height,1);
    render_framebuffer_t* fb=*slot;
    if(!fb || render_camera_setup(&v->cam,&cfg,v->sys))return -1;
    alea_compute_t* compute = v->compute;
    if (compute && (cfg.cell_filter.mode != RENDER_FILTER_ALL || cfg.material_filter.mode != RENDER_FILTER_ALL)) {
        compute = NULL;
        if (alea_prepare_query_acceleration(v->sys)) return -1;
    }
    int rc;Uint32 start=SDL_GetTicks();v->working=cfg;
#ifdef ALEA_VIEW_METAL
    if(compute && v->progressive && v->job && !v->job->saving) {
        render_view_seed(v, fb);
        v->last_publish=start;v->published_pixels=0;v->progress_error=0;
        rc=alea_compute_render_progressive(v->compute,&cfg,&v->cam,fb,render_view_progress,v);
        if(v->progress_error)rc=-1;
        if(rc==0 && v->preview) {
            double elapsed=SDL_GetTicks()-start,edge=cfg.width>cfg.height?cfg.width:cfg.height;
            if(elapsed>250)edge*=sqrt(150/elapsed);else if(elapsed<50)edge*=1.5;
            v->buffers->preview_edge=(int)fmax(16,fmin(128,edge));
        }
    } else
#endif
    {
        alea_prepared_queries_t* index = v->prepared;
        if (v->build && atomic_load(&v->build->finished)) index = v->build->index;
        int cpu_updates = !compute && v->job && !v->job->saving && cfg.render_mode != RENDER_MODE_XRAY;
        if (cpu_updates) {
            int longest = cfg.width > cfg.height ? cfg.width : cfg.height;
            if (v->working.tile_size <= 0) v->working.tile_size = RENDER_DEFAULT_TILE;
            if (v->working.tile_size > longest) v->working.tile_size = longest;
            rc = render_view_cpu_frame(v, index, fb);
        } else {
            rc = render_selected(v->sys,compute,index,&cfg,&v->cam,fb);
        }
        if (rc < 0 && index && v->build && !alea_interrupted()) {
            /* Additional tile-worker storage is allocated on first use. The
             * render worker exclusively owns the published handle per frame. */
            fprintf(stderr, "Regional rendering failed, retrying on hierarchy: %s\n", alea_error());
            alea_prepared_queries_destroy(index);
            v->build->index = NULL;
            rc = cpu_updates ? render_view_cpu_frame(v, NULL, fb) : render_scene(v->sys, &cfg, &v->cam, fb);
        }
    }
    if (rc == 0 && v->preview && v->job) v->buffers->preview_revision = v->job->revision;
    (void)start;
    if(rc==0){*pixels=render_view_image(v,&cfg,fb,NULL);if(!*pixels)rc=-1;}
    if(rc<0)fprintf(stderr,"Rendering failed: %s\n",alea_error());
    return rc;
}

static void render_view_bind(void* state,alea_view_job_t* job) {
    ((render_viewer_t*)state)->job=job;
}

static void render_view_prepare(void* state, int preview) {
    ((render_viewer_t*)state)->preview = preview;
}

static int render_view_save(const void* state, const uint8_t* pixels) {
    const render_viewer_t* v = state;
    int rc = render_write_image(v->output, pixels, v->cfg.width, v->cfg.height);
    fprintf(stderr, "%s %s\n", rc == 0 ? "Saved" : "Failed to save", v->output);
    return rc;
}

static int render_view_active(const void* state, SDL_Keycode key) {
    const render_viewer_t* v = state;
    if (key == SDLK_i) return inspector_active(v->inspector, PROBE_PICK);
    if (key == SDLK_m) return inspector_active(v->inspector, PROBE_MATERIALS);
    return (key == SDLK_x && v->cfg.num_clips) || (key == SDLK_e && v->cfg.edges) || (key == SDLK_l && v->cfg.shadows);
}

static int render_view_pick(void* snapshot, int x, int y) {
    render_viewer_t* v = snapshot;
    alea_prepared_queries_t* index = v->prepared;
    if (v->build && atomic_load(&v->build->finished)) index = v->build->index;
    if (render_camera_setup(&v->cam, &v->cfg, v->sys)) return -1;
    int rc = alea_render_pick(v->sys, index, &v->cfg, &v->cam, x, y, &v->pick);
    if (rc < 0) fprintf(stderr, "Inspection failed: %s\n", alea_error());
    return rc;
}

static const char* render_view_picked(void* state, const void* snapshot, int result) {
    render_viewer_t* v = state;
    if (result < 0) return "Inspection failed";
    if (!result) return "Background";
    v->pick = ((const render_viewer_t*)snapshot)->pick;
    v->has_pivot = 1;
    inspector_select(v->inspector, &v->pick);
    char boundary[40] = "";
    if (v->pick.clip_plane >= 0)
        snprintf(boundary, sizeof(boundary), " Clip %d", v->pick.clip_plane + 1);
    else if (v->pick.surface_id >= 0)
        snprintf(boundary, sizeof(boundary), " Surface %d", v->pick.surface_id);
    snprintf(v->inspection, sizeof(v->inspection), "Cell %d Mat %d%s  %.5g %.5g %.5g",
             v->pick.cell_id, v->pick.material_id, boundary,
             v->pick.point[0], v->pick.point[1], v->pick.point[2]);
    fprintf(stderr, "%s\n", v->inspection);
    return v->inspection;
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
        case SDLK_i: inspector_toggle(v->inspector, PROBE_PICK); return 0;
        case SDLK_m: inspector_toggle(v->inspector, PROBE_MATERIALS); return 0;
        case SDLK_r:
            *cfg = v->initial; v->has_pivot = 0;
            inspector_clear(v->inspector);
            return 1;
        case SDLK_c: cfg->color_mode = (cfg->color_mode + 1) % (v->compute ? 2 : 4); return 1;
        case SDLK_e: cfg->edges = !cfg->edges; return 1;
        case SDLK_l:
            if (v->compute) { fprintf(stderr,"Shadows are not supported by Metal yet.\n"); return 0; }
            cfg->shadows = !cfg->shadows; return 1;
        case SDLK_x:
            if (cfg->num_clips) cfg->num_clips=0;
            else {
                cfg->num_clips=1;
                /* Retain a configured/moved plane when toggling it. */
                render_clip_plane_t* p=&cfg->clips[0];
                if (p->normal[0]==0 && p->normal[1]==0 && p->normal[2]==0) {
                    p->normal[2]=1; p->d=-cfg->target[2];
                }
            }
            return 1;
        case SDLK_LEFTBRACKET: case SDLK_RIGHTBRACKET: {
            if (!cfg->num_clips) return 0;
            double distance=0;
            for (int i=0;i<3;++i) distance+=pow(cfg->eye[i]-cfg->target[i],2);
            render_clip_plane_t* p=&cfg->clips[0];
            double scale=hypot(hypot(p->normal[0],p->normal[1]),p->normal[2]);
            double step=(cfg->fov==0?cfg->ortho_height:sqrt(distance))*.02*scale;
            p->d += event->key.keysym.sym==SDLK_LEFTBRACKET ? step : -step;
            return 1;
        }
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
    /* Rotate both ends of the camera around the world-space selection. This
     * keeps an off-center hit under its original pixel during an orbit. */
    if (v->has_pivot && (yaw || pitch)) {
        double target_offset[3];
        for (int i = 0; i < 3; ++i) target_offset[i] = cfg->target[i] - v->pick.point[i];
        view_rotate(target_offset, v->cam.true_up, yaw);
        view_rotate(target_offset, v->cam.right, pitch);
        for (int i = 0; i < 3; ++i) cfg->target[i] = v->pick.point[i] + target_offset[i];
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

static int render_panel_height(const void* state) {
    return inspector_height(((const render_viewer_t*)state)->inspector);
}

static int render_panel_draw(void* state, SDL_Window* win, SDL_Renderer* r, SDL_Rect rect) {
    render_viewer_t* v = state;
    return inspector_draw(v->inspector, win, r, rect, &v->cfg);
}

static int render_panel_event(void* state, SDL_Event* event, SDL_Rect rect) {
    return inspector_event(((render_viewer_t*)state)->inspector, event, rect);
}

static void render_section_event(void* state, Uint32 id, const SDL_Event* event) {
    inspector_section_event(((render_viewer_t*)state)->inspector, id, event);
}

static int render_panel_slice(void* snapshot, unsigned char** pixels) {
    render_viewer_t* v = snapshot;
    return inspector_slice(v->sys, &v->cfg, v->inspector, &v->section, pixels);
}

static void render_panel_publish(void* state, const void* snapshot, unsigned char* pixels) {
    const render_viewer_t* v = snapshot;
    inspector_show_slice(((render_viewer_t*)state)->inspector, &v->section, pixels);
}

static void render_panel_close(void* state) {
    inspector_close_graphics(((render_viewer_t*)state)->inspector);
}

/* Render, pick and queued save jobs may outlive the next visibility edit. */
static void* render_view_snapshot(const void* state) {
    const render_viewer_t* live = state;
    size_t cells = live->cfg.cell_filter.count, mats = live->cfg.material_filter.count;
    render_viewer_t* copy = malloc(sizeof(*copy) + (cells + mats) * sizeof(int));
    if (!copy) return NULL;
    *copy = *live;
    copy->section = inspector_section_view(live->inspector);
    int* ids = (int*)(copy + 1);
    if (cells) {
        memcpy(ids, live->cfg.cell_filter.ids, cells*sizeof(*ids));
        copy->cfg.cell_filter.ids = ids;
    }
    if (mats) {
        memcpy(ids+cells, live->cfg.material_filter.ids, mats*sizeof(*ids));
        copy->cfg.material_filter.ids = ids+cells;
    }
    return copy;
}

static int render_interactive(alea_system_t* sys, alea_compute_t* compute, alea_prepared_queries_t* prepared, const render_config_t* cfg,
                              const render_camera_t* cam, const char* output, int progressive, int automatic) {
    if (cfg->width < 1 || cfg->height < 1 || cfg->width > 16384 || cfg->height > 16384)
        return -1;
    render_viewer_t v = {0};
    render_view_build_t build = {.sys=sys, .domain=render_view_domain(cam)};
    atomic_init(&build.finished, 0);
    v.build = automatic ? &build : NULL;
    render_view_buffers_t buffers = {.preview_edge=alea_cell_count(sys)>10000?16:96};
    v.buffers = &buffers; v.compute = compute;v.progressive=progressive;
    v.prepared = prepared; v.sys = sys; v.cfg = *cfg; v.cam = *cam; v.output = output;
    /* Freeze auto-fit so navigation is relative to the visible camera. */
    memcpy(v.cfg.eye, cam->eye, sizeof(v.cfg.eye));
    memcpy(v.cfg.target, cam->target, sizeof(v.cfg.target));
    memcpy(v.cfg.up, cam->true_up, sizeof(v.cfg.up));
    v.cfg.ortho_height = cam->ortho_height;
    v.cfg.eye_set = v.cfg.target_set = 1;
    v.initial = v.cfg; /* Borrowed config resources; freed only by the caller. */
    v.inspector = inspector_create(sys);
    if (!v.inspector) return -1;
    static const alea_view_button_t buttons[] = {
        {"Save", SDLK_s}, {"Reset", SDLK_r}, {"Zoom +", SDLK_PLUS},
        {"Zoom -", SDLK_MINUS}, {"Help", SDLK_h}, {"Quit", SDLK_q},
        {"Left", SDLK_LEFT}, {"Right", SDLK_RIGHT}, {"Up", SDLK_UP},
        {"Down", SDLK_DOWN}, {"Colors", SDLK_c}, {"Edges", SDLK_e},
        {"Clip", SDLK_x}, {"Clip -", SDLK_LEFTBRACKET}, {"Clip +", SDLK_RIGHTBRACKET},
        {"Shadows", SDLK_l}, {"Probe", SDLK_i}, {"Materials", SDLK_m}
    };
    const char* title = compute ? "Alea 3D - Metal" :
                        prepared ? "Alea 3D - Regional CPU" : "Alea 3D - CPU";
    static const alea_view_panel_t panel = {
        render_panel_height, render_panel_draw, render_panel_event, render_section_event, render_panel_slice,
        render_panel_publish, render_panel_close, render_view_snapshot
    };
    const alea_view_lifecycle_t lifecycle = {
        .preview_shown = automatic ? render_view_preview_shown : NULL,
        .finish = automatic ? render_view_finish_build : NULL,
        .pick = render_view_pick, .picked = render_view_picked, .panel = &panel,
        .delay_after_preview = !progressive
    };
    int rc = alea_sdl_view_progressive(title,
        "Click: inspect and set orbit pivot   Left drag / arrows: orbit\n"
        "Right / middle drag: pan   Wheel / +/-: zoom\n"
        "C: colors   E: edges   L: shadows (CPU)\n"
        "I: probe   M: materials   X: clip   [ / ]: move clip\n"
        "R: reset   S: save   H: help   Q / Esc: quit",
        cfg->width, cfg->height, &v, sizeof(v), render_view_draw,
        render_view_save, render_view_input, render_view_active,
        buttons, (int)(sizeof(buttons)/sizeof(buttons[0])), render_view_prepare,
        (progressive || cfg->render_mode != RENDER_MODE_XRAY) ? render_view_bind : NULL, &lifecycle);
    inspector_destroy(v.inspector);
    render_framebuffer_free(buffers.frames[0]); render_framebuffer_free(buffers.frames[1]);
    render_framebuffer_free(buffers.display[0]);render_framebuffer_free(buffers.display[1]);
    return rc;
}
#endif
