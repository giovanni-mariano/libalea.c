// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
// SPDX-License-Identifier: MPL-2.0

#ifndef ALEA_RENDER_VIEW_H
#define ALEA_RENDER_VIEW_H
#include "sdl_view.h"

typedef struct {
    render_framebuffer_t* frames[2]; /* Only the serial render worker accesses these. */
    render_framebuffer_t* display[2];
    uint64_t preview_revision;
    int preview_edge;
} render_view_buffers_t;

typedef struct {
    alea_system_t* sys;
    alea_compute_t* compute;
    render_view_buffers_t* buffers;
    int preview, progressive, progress_error;
    alea_view_job_t* job;
    Uint32 last_publish;
    size_t published_pixels;
    render_config_t cfg, initial, working;
    render_camera_t cam;
    const char* output;
} render_viewer_t;

/* Postprocessing is destructive. Keep raw progressive pixels separate so
 * repeatedly displaying a partial image cannot repeatedly darken its edges. */
static uint8_t* render_view_image(render_viewer_t* v, const render_config_t* cfg,
                                 const render_framebuffer_t* raw) {
    render_framebuffer_t** slot=&v->buffers->display[v->preview];
    if(*slot && ((*slot)->width!=raw->width || (*slot)->height!=raw->height)) {
        render_framebuffer_free(*slot);*slot=NULL;
    }
    if(!*slot)*slot=render_framebuffer_create(raw->width,raw->height,1);
    render_framebuffer_t* fb=*slot;if(!fb)return NULL;
    size_t n=(size_t)raw->width*raw->height;
    memcpy(fb->color,raw->color,n*3*sizeof(float));
    memcpy(fb->normal,raw->normal,n*3*sizeof(float));
    memcpy(fb->depth,raw->depth,n*sizeof(float));
    memcpy(fb->cell_id,raw->cell_id,n*sizeof(int));
    memcpy(fb->material_id,raw->material_id,n*sizeof(int));
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
    uint8_t* rgb=render_view_image(v,&v->working,fb);
    if(!rgb){v->progress_error=1;return 1;}
    alea_view_publish(v->job,rgb,(int)(100*completed/total));
    v->last_publish=SDL_GetTicks();v->published_pixels=completed;
    return alea_view_cancelled(v->job);
}
#endif
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
    int rc;Uint32 start=SDL_GetTicks();v->working=cfg;
#ifdef ALEA_VIEW_METAL
    if(v->progressive && v->job && !v->job->saving) {
        {
            const render_framebuffer_t* seed=!v->preview && v->buffers->preview_revision==v->job->revision?v->buffers->frames[1]:NULL;
            for(int y=0;y<fb->height;++y)for(int x=0;x<fb->width;++x) {
                size_t i=(size_t)y*fb->width+x;
                size_t j=seed?((size_t)y*seed->height/fb->height)*seed->width+(size_t)x*seed->width/fb->width:0;
                fb->cell_id[i]=seed?seed->cell_id[j]:-1;fb->material_id[i]=seed?seed->material_id[j]:-1;
                fb->depth[i]=seed?seed->depth[j]:INFINITY;
                for(int a=0;a<3;++a){fb->color[3*i+a]=seed?seed->color[3*j+a]:cfg.background[a];fb->normal[3*i+a]=seed?seed->normal[3*j+a]:0;}
            }
        }
        v->last_publish=start;v->published_pixels=0;v->progress_error=0;
        rc=alea_compute_render_progressive(v->compute,&cfg,&v->cam,fb,render_view_progress,v);
        if(v->progress_error)rc=-1;
        if(rc==0 && v->preview) {
            v->buffers->preview_revision=v->job->revision;
            double elapsed=SDL_GetTicks()-start,edge=cfg.width>cfg.height?cfg.width:cfg.height;
            if(elapsed>250)edge*=sqrt(150/elapsed);else if(elapsed<50)edge*=1.5;
            v->buffers->preview_edge=(int)fmax(16,fmin(128,edge));
        }
    } else
#endif
        rc=render_selected(v->sys,v->compute,&cfg,&v->cam,fb);
    (void)start;
    if(rc==0){*pixels=render_view_image(v,&cfg,fb);if(!*pixels)rc=-1;}
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
    return (key == SDLK_x && v->cfg.num_clips) || (key == SDLK_e && v->cfg.edges) || (key == SDLK_l && v->cfg.shadows);
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
        case SDLK_r: *cfg = v->initial; return 1;
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

static int render_interactive(alea_system_t* sys, alea_compute_t* compute, const render_config_t* cfg,
                              const render_camera_t* cam, const char* output, int progressive) {
    if (cfg->width < 1 || cfg->height < 1 || cfg->width > 16384 || cfg->height > 16384)
        return -1;
    render_viewer_t v = {0};
    render_view_buffers_t buffers = {.preview_edge=alea_cell_count(sys)>10000?16:96};
    v.buffers = &buffers; v.compute = compute;v.progressive=progressive;
    v.sys = sys; v.cfg = *cfg; v.cam = *cam; v.output = output;
    /* Freeze auto-fit so navigation is relative to the visible camera. */
    memcpy(v.cfg.eye, cam->eye, sizeof(v.cfg.eye));
    memcpy(v.cfg.target, cam->target, sizeof(v.cfg.target));
    memcpy(v.cfg.up, cam->true_up, sizeof(v.cfg.up));
    v.cfg.ortho_height = cam->ortho_height;
    v.cfg.eye_set = v.cfg.target_set = 1;
    v.initial = v.cfg; /* Borrowed config resources; freed only by the caller. */
    static const alea_view_button_t buttons[] = {
        {"Save", SDLK_s}, {"Reset", SDLK_r}, {"Zoom +", SDLK_PLUS},
        {"Zoom -", SDLK_MINUS}, {"Help", SDLK_h}, {"Quit", SDLK_q},
        {"Left", SDLK_LEFT}, {"Right", SDLK_RIGHT}, {"Up", SDLK_UP},
        {"Down", SDLK_DOWN}, {"Colors", SDLK_c}, {"Edges", SDLK_e},
        {"Clip", SDLK_x}, {"Clip -", SDLK_LEFTBRACKET}, {"Clip +", SDLK_RIGHTBRACKET},
        {"Shadows", SDLK_l}
    };
    int rc = alea_sdl_view_progressive(compute ? "Alea 3D - Metal" : "Alea 3D - CPU",
        "Left drag / arrows: orbit\n"
        "Right / middle drag: pan   Wheel / +/-: zoom\n"
        "C: colors   E: edges   L: shadows (CPU)\n"
        "X: toggle clip   [ / ]: move clip\n"
        "R: reset   S: save   H: help   Q / Esc: quit",
        cfg->width, cfg->height, &v, sizeof(v), render_view_draw,
        render_view_save, render_view_input, render_view_active,
        buttons, (int)(sizeof(buttons)/sizeof(buttons[0])), render_view_prepare, progressive?render_view_bind:NULL);
    render_framebuffer_free(buffers.frames[0]); render_framebuffer_free(buffers.frames[1]);
    render_framebuffer_free(buffers.display[0]);render_framebuffer_free(buffers.display[1]);
    return rc;
}
#endif
