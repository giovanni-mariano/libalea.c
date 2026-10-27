// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
// SPDX-License-Identifier: MPL-2.0
#define _POSIX_C_SOURCE 200809L
#include "alea.h"
#include "alea_mcnp.h"
#include "alea_compute.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <errno.h>
#include <signal.h>
static double seconds(void){struct timespec t;clock_gettime(CLOCK_MONOTONIC,&t);return t.tv_sec+t.tv_nsec*1e-9;}
static void stop(int sig){(void)sig;alea_interrupt();}
static int integer(const char* s,int low,int high){char* e;errno=0;long v=strtol(s,&e,10);return errno || e==s || *e || v<low || v>high?-1:(int)v;}
typedef struct {double start,last,first,max_gap,cancel;size_t done,calls;} progress_state;
static int progress(void* user,const render_framebuffer_t* image,size_t done,size_t total) {
    (void)image;progress_state* p=user;double now=seconds();
    if(done){if(!p->done){p->first=now-p->start;printf("First %zu pixels after %.3f s\n",done,p->first);}}
    if(p->calls && now-p->last>p->max_gap)p->max_gap=now-p->last;
    p->last=now;p->done=done;++p->calls;
    if(done==total)printf("Completed %zu pixels\n",total);
    return p->cancel>0 && now-p->start>=p->cancel;
}
int main(int argc,char** argv) {
    if(argc!=2 && argc!=4 && argc!=5 && argc!=6){fprintf(stderr,"Usage: %s MODEL [WIDTH HEIGHT [CANCEL_MS [VIEW]]]\nITER VIEW 1=cutaway, 2=opposite cutaway, 3=exterior, 4=component close-up. Zero CANCEL_MS compares completed pixels against ordinary rendering.\n",argv[0]);return 1;}
    int w=argc>=4?integer(argv[2],1,4096):16,h=argc>=4?integer(argv[3],1,4096):12,ms=argc>=5?integer(argv[4],0,3600000):0,view=argc==6?integer(argv[5],1,4):1;
    if(w<0 || h<0 || ms<0 || view<0)return 1;
    setbuf(stdout,NULL);signal(SIGINT,stop);signal(SIGTERM,stop);
    unsetenv("ALEA_METAL_PROFILE");
    double start=seconds();mcnp_model_t* model=mcnp_load(argv[1]);if(!model){fprintf(stderr,"%s\n",alea_error());return 1;}
    printf("Load %.3f s\n",seconds()-start);
    render_config_t cfg;render_config_init(&cfg);cfg.width=w;cfg.height=h;cfg.log_level=0;
    cfg.eye_set=cfg.target_set=1;cfg.eye[0]=1800;cfg.eye[1]=-2200;cfg.eye[2]=1400;
    cfg.num_clips=1;cfg.clips[0]=(render_clip_plane_t){{0,1,0},0};
    if(view==2){cfg.eye[0]=-1800;cfg.eye[1]=2200;cfg.eye[2]=1000;cfg.clips[0].normal[1]=-1;}
    if(view==3)cfg.num_clips=0;
    if(view==4){cfg.eye[0]=300;cfg.eye[1]=80;cfg.eye[2]=-120;cfg.target[0]=300;cfg.target[1]=300;cfg.target[2]=-200;cfg.num_clips=0;}
    printf("View %d, %dx%d, AA1, cancel %d ms\n",view,w,h,ms);
    render_camera_t camera;alea_compute_t* gpu=NULL;render_framebuffer_t* image=NULL;render_framebuffer_t* baseline=NULL;int status=1;
    if(render_camera_setup(&camera,&cfg,model->sys))goto done;
    start=seconds();gpu=alea_compute_create(model->sys,ALEA_COMPUTE_METAL,ALEA_COMPUTE_FP32);if(!gpu)goto done;
    printf("Prepare %.3f s\n",seconds()-start);
    image=render_framebuffer_create(w,h,1);if(!image)goto done;
    progress_state p={.start=seconds(),.cancel=ms/1000.0};p.last=p.start;
    int rc=alea_compute_render_progressive(gpu,&cfg,&camera,image,progress,&p);if(rc<0)goto done;
    printf("%s after %.3f s: %zu/%d pixels, %zu callbacks, first pixels %.3f s, largest callback gap %.3f s\n",rc?"Cancelled":"Finished",seconds()-p.start,p.done,w*h,p.calls,p.first,p.max_gap);
    if(rc){status=0;goto done;}
    baseline=render_framebuffer_create(w,h,1);if(!baseline)goto done;
    start=seconds();if(alea_compute_render(gpu,&cfg,&camera,baseline))goto done;
    printf("Ordinary render %.3f s\n",seconds()-start);
    size_t differences=0,n=(size_t)w*h;
    for(size_t i=0;i<n;++i)differences+=image->cell_id[i]!=baseline->cell_id[i] || image->material_id[i]!=baseline->material_id[i] ||
        memcmp(image->color+3*i,baseline->color+3*i,3*sizeof(float)) || memcmp(image->normal+3*i,baseline->normal+3*i,3*sizeof(float)) ||
        memcmp(image->depth+i,baseline->depth+i,sizeof(float));
    printf("Pixels differing from ordinary Metal render: %zu/%zu\n",differences,n);status=differences?2:0;
done:
    if(status==1)fprintf(stderr,"%s\n",alea_error());
    render_framebuffer_free(image);render_framebuffer_free(baseline);alea_compute_destroy(gpu);render_config_free(&cfg);mcnp_model_destroy(model);return status;
}
