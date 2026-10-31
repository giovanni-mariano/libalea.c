// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
// SPDX-License-Identifier: MPL-2.0
/* One bounded diagnostic context per process; never linked into the package. */
#import <Foundation/Foundation.h>
#import <Metal/Metal.h>
static id<MTLBuffer> ray_audit_maps,ray_audit_rows;
static int ray_audit_enabled;
#define COMPUTE_RENDER_DIAGNOSTIC_BIND(e) do {if(ray_audit_enabled) { \
    [e setBuffer:ray_audit_maps offset:0 atIndex:17]; \
    [e setBuffer:ray_audit_rows offset:0 atIndex:18];     RAY_EXECUTE_BIND(e); }} while(0)
#ifndef RAY_EXECUTE_BIND
#define RAY_EXECUTE_BIND(e)
#endif
#include "../src/compute/compute_metal.m"
#include "alea_mcnp.h"
#include "gpu_geometry_fixture.h"
#ifdef ALEA_RAY_EXECUTION_EXPERIMENT
#include "ray_group_execute_metallib.h"
#define alea_ray_group_work_metallib alea_ray_group_execute_metallib
#define alea_ray_group_work_metallib_len alea_ray_group_execute_metallib_len
#else
#include "ray_group_work_metallib.h"
#endif
#include "ray_contract_oracle.h"
#include <inttypes.h>
#include <ctype.h>
#include <time.h>
#include <sys/resource.h>
#include <signal.h>
static double ray_seconds(void) {struct timespec t;clock_gettime(CLOCK_MONOTONIC,&t);return t.tv_sec+t.tv_nsec*1e-9;}
static void ray_interrupted(int sig) {(void)sig;alea_interrupt();}
#ifdef ALEA_RAY_EXECUTION_EXPERIMENT
#include "ray_lattice_fixture.h"
#include "ray_execute_driver.h"
#else
static int audit_render(metal_state* s,const render_config_t* cfg,const render_camera_t* camera,
                        render_framebuffer_t* fb,const char* path) {
    const size_t n=(size_t)fb->width*fb->height;int rc=-1;compute_render_pixel* reference=NULL;
    NSError* error=nil;id<MTLLibrary> library=nil;id<MTLComputePipelineState> pipeline=nil;
    dispatch_data_t data=dispatch_data_create(alea_ray_group_work_metallib,alea_ray_group_work_metallib_len,NULL,DISPATCH_DATA_DESTRUCTOR_DEFAULT);
    library=[s->device newLibraryWithData:data error:&error];dispatch_release(data);
    if(!library || ray_contract_check(s,library))goto done;
    if(alea_compute_metal_render(s,cfg,camera,fb,NULL,NULL))goto done;
    reference=malloc(n*sizeof(*reference));if(!reference)goto done;memcpy(reference,s->pixels.contents,n*sizeof(*reference));
    MTLFunctionConstantValues* constants=[[MTLFunctionConstantValues alloc] init];bool yes=true,no=false;
    [constants setConstantValue:&yes type:MTLDataTypeBool atIndex:0];
    [constants setConstantValue:&s->lazy_threshold type:MTLDataTypeUInt atIndex:2];
    [constants setConstantValue:&s->refill type:MTLDataTypeBool atIndex:4];
    [constants setConstantValue:&s->compact_code type:MTLDataTypeBool atIndex:7];
    [constants setConstantValue:&no type:MTLDataTypeBool atIndex:8];
    id<MTLFunction> f=[library newFunctionWithName:@"alea_render_image" constantValues:constants error:&error];
    if(f)pipeline=[s->device newComputePipelineStateWithFunction:f error:&error];[f release];[constants release];
    if(!pipeline)goto done;
    ray_audit_maps=[s->group_maps retain];
    if(!ray_audit_maps) {const uint32_t empty[5]={0};ray_audit_maps=[s->device newBufferWithBytes:empty length:sizeof(empty) options:MTLResourceStorageModeShared];}
    ray_audit_rows=[s->device newBufferWithLength:n*8*4 options:MTLResourceStorageModeShared];
    if(!ray_audit_maps || !ray_audit_rows)goto done;
    if(ray_audit_maps.length+s->group_refs.length+ray_audit_rows.length>64*1024*1024 ||
       s->device.currentAllocatedSize>1024ull*1024*1024)goto done;
    [s->render_pipeline release];s->render_pipeline=pipeline;pipeline=nil;ray_audit_enabled=1;
    if(alea_compute_metal_render(s,cfg,camera,fb,NULL,NULL))goto done;
    if(memcmp(reference,s->pixels.contents,n*sizeof(*reference))){fprintf(stderr,"Instrumented render changed outputs\n");goto done;}
    uint64_t totals[6]={0};const uint32_t* rows=ray_audit_rows.contents;
    const compute_profile* profiles=s->ray_costs.contents;uint64_t instructions=0,primitives=0;
    for(size_t i=0;i<n;++i) {
        for(unsigned k=0;k<6;++k){if(rows[8*i+k]==UINT32_MAX)goto done;totals[k]+=rows[8*i+k];}
        if(profiles[i].instructions==UINT32_MAX || profiles[i].primitives==UINT32_MAX)goto done;
        instructions+=profiles[i].instructions;primitives+=profiles[i].primitives;
    }
    if(totals[1]!=instructions || totals[2]!=primitives || totals[3]>totals[1] || totals[4]>totals[2]) {
        fprintf(stderr,"Ray work totals failed reconciliation\n");goto done;
    }
    FILE* out=fopen(path,"w");if(!out)goto done;
    fprintf(out,"pixels,aa,expression_calls,instructions,primitives,group_instructions,group_primitives,group_entries\n");
    fprintf(out,"%zu,%d",n,cfg->aa_samples);for(unsigned k=0;k<6;++k)fprintf(out,",%" PRIu64,totals[k]);fputc('\n',out);
    int failed=ferror(out);failed|=fclose(out)!=0;if(failed)goto done;
    printf("pixels=%zu group_instruction_percent=%.3f group_primitive_percent=%.3f outputs=identical totals=reconciled saturation=0\n",n,
        totals[1]?100.*totals[3]/totals[1]:0,totals[2]?100.*totals[4]/totals[2]:0);rc=0;
 done:
    if(rc)fprintf(stderr,"Ray audit failed: %s %s\n",alea_error(),error?error.localizedDescription.UTF8String:"");
    ray_audit_enabled=0;[ray_audit_maps release];[ray_audit_rows release];ray_audit_maps=ray_audit_rows=nil;
    [pipeline release];[library release];free(reference);return rc;
}
#endif
int main(int argc,char** argv) {
    if(argc!=2 && argc!=4){fprintf(stderr,"Usage: %s OUTPUT.csv [MODEL CAMERA.txt]\n",argv[0]);return 2;}
    @autoreleasepool {
        (void)check_render;(void)check_points;setbuf(stdout,NULL);
        signal(SIGINT,ray_interrupted);signal(SIGTERM,ray_interrupted);
        const char* threading=getenv("ALEA_METAL_JUMP_THREADING");
        if(threading && strcmp(threading,"0")){fprintf(stderr,"Ray coverage requires unthreaded code\n");return 2;}
        if(setenv("ALEA_METAL_CSG_GROUPS","1",1))return 2;
#ifdef ALEA_RAY_EXECUTION_EXPERIMENT
        unsetenv("ALEA_METAL_PROFILE");
        if(getenv("ALEA_BENCH_RAY_PROFILE"))setenv("ALEA_METAL_PROFILE",argv[1],1);
#else
        setenv("ALEA_METAL_PROFILE",argv[1],1);
#endif
        unsetenv("ALEA_METAL_ROOT_REUSE");
        if(argc==2)setenv("ALEA_BENCH_CSG_GROUPS","1",1);
        mcnp_model_t* model=argc==4?mcnp_load(argv[2]):NULL;
        alea_system_t* sys=argc==4?(model?model->sys:NULL):
#ifdef ALEA_RAY_EXECUTION_EXPERIMENT
            getenv("ALEA_BENCH_RAY_LATTICE")?ray_lattice_scene():
#endif
            scene();
        metal_state* s=NULL;render_framebuffer_t* fb=NULL;int rc=1;
        render_config_t cfg;render_config_init(&cfg);cfg.width=argc==4?4:32;cfg.height=argc==4?3:24;cfg.log_level=0;
        cfg.eye_set=cfg.target_set=1;cfg.eye[0]=0;cfg.eye[1]=-24;cfg.eye[2]=12;cfg.target[0]=cfg.target[1]=cfg.target[2]=0;
#ifdef ALEA_RAY_EXECUTION_EXPERIMENT
        if(argc==2 && getenv("ALEA_BENCH_RAY_CLIP")){cfg.num_clips=1;cfg.clips[0]=(render_clip_plane_t){{0,1,0},0};}
#endif
        if(argc==4) {
            double values[10];FILE* camera=fopen(argv[3],"r");if(!camera)goto done;
            bool valid=true;for(int j=0;j<10;++j)if(fscanf(camera,"%lf",values+j)!=1 || !isfinite(values[j]))valid=false;
            int trailing;do{trailing=fgetc(camera);}while(trailing!=EOF && isspace(trailing));
            fclose(camera);if(!valid || trailing!=EOF)goto done;
            for(int j=0;j<3;++j){cfg.eye[j]=values[j];cfg.target[j]=values[3+j];cfg.clips[0].normal[j]=values[6+j];}
            cfg.clips[0].d=values[9];cfg.num_clips=values[6]!=0 || values[7]!=0 || values[8]!=0;
        }
        render_camera_t camera;if(!sys || render_camera_setup(&camera,&cfg,sys))goto done;
        double preparation=ray_seconds();s=alea_compute_metal_create(sys);preparation=ray_seconds()-preparation;
        printf("preparation_seconds=%.6f\n",preparation);
        if(preparation>15)goto done;
        fb=render_framebuffer_create(cfg.width,cfg.height,1);
        if(!s || !fb || audit_render(s,&cfg,&camera,fb,argv[1]))goto done;
        struct rusage usage;if(getrusage(RUSAGE_SELF,&usage) || (uint64_t)usage.ru_maxrss>4ull*1024*1024*1024)goto done;
        printf("peak_process_bytes=%ld device_allocated_bytes=%lu\n",usage.ru_maxrss,(unsigned long)s->device.currentAllocatedSize);
        puts("result=PASS");rc=0;
 done:
        if(rc)fprintf(stderr,"Ray diagnostic failed: %s\n",alea_error());
        render_framebuffer_free(fb);alea_compute_metal_destroy(s);
        if(model)mcnp_model_destroy(model);else alea_destroy(sys);return rc;
    }
}
