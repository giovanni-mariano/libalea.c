// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
// SPDX-License-Identifier: MPL-2.0
/* Native Metal diagnostic. Fixed queries/camera while dispatch size changes.
 * Private bridge access controls the existing resident pool without adding
 * experimental scheduling options to the library API. */
#include "../src/compute/compute_metal.m"
#include "gpu_geometry_fixture.h"

typedef struct {double seconds,maximum;size_t commands,timed;} gpu_times;
static gpu_times timing;
@interface TimedCommand : NSProxy { @public id<MTLCommandBuffer> target; }
@end
@implementation TimedCommand
- (NSMethodSignature*)methodSignatureForSelector:(SEL)s {return [(NSObject*)target methodSignatureForSelector:s];}
- (void)forwardInvocation:(NSInvocation*)i {[i invokeWithTarget:target];}
- (void)waitUntilCompleted {
    [target waitUntilCompleted];++timing.commands;
    double begin=target.GPUStartTime,end=target.GPUEndTime;
    if(begin>0 && end>begin){++timing.timed;timing.seconds+=end-begin;timing.maximum=fmax(timing.maximum,end-begin);}
}
- (void)dealloc {[target release];[super dealloc];}
@end
@interface TimedQueue : NSProxy { @public id<MTLCommandQueue> target; }
@end
@implementation TimedQueue
- (NSMethodSignature*)methodSignatureForSelector:(SEL)s {return [(NSObject*)target methodSignatureForSelector:s];}
- (void)forwardInvocation:(NSInvocation*)i {[i invokeWithTarget:target];}
- (id<MTLCommandBuffer>)commandBuffer {
    id<MTLCommandBuffer> command=[target commandBuffer];if(!command)return nil;
    TimedCommand* p=[TimedCommand alloc];p->target=[command retain];return (id<MTLCommandBuffer>)[p autorelease];
}
- (void)dealloc {[target release];[super dealloc];}
@end
static double now(void) {struct timespec t;clock_gettime(CLOCK_MONOTONIC,&t);return t.tv_sec+t.tv_nsec*1e-9;}
static void attach_timer(metal_state* s) {
    TimedQueue* q=[TimedQueue alloc];q->target=s->queue;s->queue=(id<MTLCommandQueue>)q;
}
static int same_frame(const render_framebuffer_t* a,const render_framebuffer_t* b) {
    size_t n=(size_t)a->width*a->height;
    return a->width==b->width && a->height==b->height &&
        !memcmp(a->cell_id,b->cell_id,n*sizeof(int)) && !memcmp(a->material_id,b->material_id,n*sizeof(int)) &&
        !memcmp(a->depth,b->depth,n*sizeof(float)) && !memcmp(a->normal,b->normal,3*n*sizeof(float)) &&
        !memcmp(a->color,b->color,3*n*sizeof(float));
}
static void print_time(const char* operation,size_t batch,int pass,size_t queries,double elapsed) {
    printf("operation=%s batch=%zu pass=%d queries=%zu wall_ms=%.6f gpu_ms=%.6f commands=%zu timed=%zu max_command_ms=%.6f queries_per_s=%.3f\n",
        operation,batch,pass,queries,elapsed*1000,timing.seconds*1000,timing.commands,timing.timed,timing.maximum*1000,queries/elapsed);
}
typedef struct {double start,first,max_gap,last;size_t completed;int cancel;} progress_state;
static int progress(void* user,const render_framebuffer_t* fb,size_t completed,size_t total) {
    (void)fb;(void)total;progress_state* p=user;double t=now();
    if(completed && !p->completed)p->first=t-p->start;
    p->max_gap=fmax(p->max_gap,t-p->last);p->last=t;p->completed=completed;
    return p->cancel && t-p->start>=.1;
}
static void camera(render_config_t* cfg,int width,int height) {
    render_config_init(cfg);cfg->width=width;cfg->height=height;cfg->log_level=0;
    cfg->eye_set=cfg->target_set=1;cfg->eye[0]=16;cfg->eye[1]=-20;cfg->eye[2]=18;
    cfg->target[0]=cfg->target[1]=cfg->target[2]=0;
}
static int points(metal_state* s,size_t batch,int passes,int selected) {
    const size_t n=128*128;double* xyz=malloc(n*3*sizeof(double));
    alea_compute_result_t* all=malloc(n*sizeof(*all));alea_compute_selection_t* one=malloc(n*sizeof(*one));
    if(!xyz || !all || !one){free(xyz);free(all);free(one);return -1;}
    int rc=-1;
    for(size_t i=0;i<n;++i){xyz[3*i]=-8+((i%128)+.5)/8;xyz[3*i+1]=-8+((i/128)+.5)/8;xyz[3*i+2]=0;}
    for(int pass=0;pass<passes;++pass) {
        timing=(gpu_times){0};double start=now();
        for(size_t i=0;i<n;i+=batch) {
            size_t count=MIN(batch,n-i);
            if(selected?alea_compute_metal_select(s,xyz+3*i,count,one+i):alea_compute_metal_classify(s,xyz+3*i,count,all+i))goto done;
        }
        double elapsed=now()-start;
        if(selected)for(size_t i=0;i<n;++i) {
            if(one[i].status!=ALEA_COMPUTE_FOUND || one[i].depth!=1)goto done;
            all[i]=(alea_compute_result_t){one[i].cell_id,one[i].material_id,1};
        }
        if(check_points(xyz,all,n))goto done;
        print_time(selected?"selected":"all-owner",batch,pass,n,elapsed);
    }
    rc=0;
done:free(xyz);free(all);free(one);return rc;
}
static int rays(alea_system_t* sys,metal_state* s,size_t batch,int passes,int resumable) {
    render_config_t cfg;camera(&cfg,192,128);render_camera_t cam;
    render_framebuffer_t* fb=render_framebuffer_create(cfg.width,cfg.height,1);
    render_framebuffer_t* reference=render_framebuffer_create(cfg.width,cfg.height,1);int rc=-1;
    if(!fb || !reference || render_camera_setup(&cam,&cfg,sys))goto done;
    /* Untimed complete-frame reference: identical ray set, no chunk-dependent oracle. */
    if(alea_compute_metal_render(s,&cfg,&cam,reference,NULL,NULL) || check_render(&cfg,&cam,reference))goto done;
    if(resumable){if(prepare_resume(s))goto done;s->resume_capacity=MIN(batch,s->resume_capacity);}
    else s->render_batch_limit=batch;
    for(int pass=0;pass<passes;++pass) {
        timing=(gpu_times){0};double start=now();progress_state p={.start=start,.last=start};
        if(alea_compute_metal_render(s,&cfg,&cam,fb,resumable?progress:NULL,&p))goto done;
        double elapsed=now()-start;
        if(!same_frame(fb,reference) || check_render(&cfg,&cam,fb))goto done;
        print_time(resumable?"resume":"rays",batch,pass,(size_t)cfg.width*cfg.height,elapsed);
        if(resumable)printf("first_pixels_ms=%.6f callback_gap_ms=%.6f completed=%zu resident_slots=%zu\n",p.first*1000,p.max_gap*1000,p.completed,s->resume_capacity);
    }
    rc=0;
done:render_framebuffer_free(fb);render_framebuffer_free(reference);render_config_free(&cfg);return rc;
}
static int interactive(alea_system_t* sys,metal_state* s,int passes) {
    render_config_t cfg;camera(&cfg,96,64);int rc=-1;
    render_framebuffer_t* preview=render_framebuffer_create(96,64,1);
    render_framebuffer_t* detail=render_framebuffer_create(384,256,1);
    render_framebuffer_t* reference=render_framebuffer_create(384,256,1);
    if(!preview || !detail || !reference)goto done;
    for(int pass=0;pass<passes;++pass) {
        cfg.width=96;cfg.height=64;cfg.eye[0]=16+.1*pass;render_camera_t cam;
        if(render_camera_setup(&cam,&cfg,sys))goto done;
        timing=(gpu_times){0};double start=now();progress_state p={.start=start,.last=start};
        if(alea_compute_metal_render(s,&cfg,&cam,preview,progress,&p))goto done;
        double preview_time=now()-start;
        if(check_render(&cfg,&cam,preview))goto done;
        print_time("preview",0,pass,96*64,preview_time);
        printf("first_pixels_ms=%.6f callback_gap_ms=%.6f\n",p.first*1000,p.max_gap*1000);
        cfg.width=384;cfg.height=256;if(render_camera_setup(&cam,&cfg,sys))goto done;
        timing=(gpu_times){0};start=now();p=(progress_state){.start=start,.last=start};
        if(alea_compute_metal_render(s,&cfg,&cam,detail,progress,&p))goto done;
        double elapsed=now()-start;
        print_time("detail",0,pass,384*256,elapsed);
        if(check_render(&cfg,&cam,detail) || alea_compute_metal_render(s,&cfg,&cam,reference,NULL,NULL) || !same_frame(detail,reference))goto done;
    }
    /* Force long work units to exercise cancellation before a detailed image
     * completes, then render a new camera to detect leaked continuation state. */
    s->resume_fixed=1;s->resume_work=1;
    render_camera_t cam;if(render_camera_setup(&cam,&cfg,sys))goto done;
    double start=now();progress_state p={.start=start,.last=start,.cancel=1};
    int cancelled=alea_compute_metal_render(s,&cfg,&cam,detail,progress,&p);
    printf("cancel_result=%d cancel_threshold_ms=100 callback_cancel_overrun_ms=%.6f completed=%zu\n",cancelled,(now()-start-.1)*1000,p.completed);
    if(cancelled!=1 || p.completed==(size_t)cfg.width*cfg.height)goto done;
    s->resume_fixed=0;s->resume_work=128;cfg.eye[0]+=1;
    if(render_camera_setup(&cam,&cfg,sys))goto done;
    start=now();p=(progress_state){.start=start,.last=start};
    if(alea_compute_metal_render(s,&cfg,&cam,detail,progress,&p) ||
        alea_compute_metal_render(s,&cfg,&cam,reference,NULL,NULL) || !same_frame(detail,reference))goto done;
    rc=0;
done:render_framebuffer_free(preview);render_framebuffer_free(detail);render_framebuffer_free(reference);render_config_free(&cfg);return rc;
}
int main(int argc,char** argv) {
    if(argc!=4){fprintf(stderr,"Usage: %s all-owner|selected|rays|resume|interactive BATCH PASSES\n",argv[0]);return 2;}
    char* end;errno=0;unsigned long batch=strtoul(argv[2],&end,10);
    if(errno || *end || end==argv[2] || batch<1 || batch>65536)return 2;
    errno=0;long passes=strtol(argv[3],&end,10);if(errno || *end || end==argv[3] || passes<1 || passes>10)return 2;
    if(strcmp(argv[1],"all-owner") && strcmp(argv[1],"selected") && strcmp(argv[1],"rays") && strcmp(argv[1],"resume") && strcmp(argv[1],"interactive"))return 2;
    @autoreleasepool {
        setbuf(stdout,NULL);alea_system_t* sys=scene();if(!sys)return 1;
        double start=now();metal_state* s=alea_compute_metal_create(sys);
        if(!s){fprintf(stderr,"%s\n",alea_error());alea_destroy(sys);return 1;}
        printf("prepare_ms=%.6f\n",(now()-start)*1000);attach_timer(s);
        int rc=!strcmp(argv[1],"interactive")?interactive(sys,s,(int)passes):
            (!strcmp(argv[1],"all-owner") || !strcmp(argv[1],"selected"))?points(s,batch,(int)passes,!strcmp(argv[1],"selected")):
            rays(sys,s,batch,(int)passes,!strcmp(argv[1],"resume"));
        printf("validation=%s buffers=%zu\n",rc?"FAIL":"PASS",context_buffer_bytes(s));
        if(rc)fprintf(stderr,"%s\n",alea_error());
        alea_compute_metal_destroy(s);alea_destroy(sys);return rc?1:0;
    }
}
