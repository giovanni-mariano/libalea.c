// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
// SPDX-License-Identifier: MPL-2.0
/* Benchmark-only display approximation. Never supplies geometry query results. */
#import <Foundation/Foundation.h>
#import <Metal/Metal.h>
#include <time.h>
#include "gpu_geometry_fixture.h"
#include "reprojection_layout.h"
#include "reprojection_metallib.h"

_Static_assert(sizeof(preview_point)==28,"Point layout");
_Static_assert(sizeof(preview_params)==72,"Camera layout");
typedef struct {
    id<MTLDevice> device;
    id<MTLCommandQueue> queue;
    id<MTLComputePipelineState> pipeline[4];
    id<MTLBuffer> points,projected,depths,owners,color;
    size_t capacity;
} preview_gpu;
static double now(void) {struct timespec t;clock_gettime(CLOCK_MONOTONIC,&t);return t.tv_sec+t.tv_nsec*1e-9;}
static void destroy(preview_gpu* g) {
    [g->points release];[g->projected release];[g->depths release];[g->owners release];[g->color release];
    for(int i=0;i<4;++i)[g->pipeline[i] release];[g->queue release];[g->device release];
}
static int create(preview_gpu* g,size_t capacity) {
    /* Includes CPU readback and point staging, in addition to GPU buffers. */
    if(!capacity || capacity>64*1024*1024/(2*sizeof(preview_point)+sizeof(preview_projection)+40))return -1;
    g->capacity=capacity;g->device=MTLCreateSystemDefaultDevice();if(!g->device){fprintf(stderr,"Metal device unavailable\n");return -1;}
    g->queue=[g->device newCommandQueue];NSError* error=nil;
    dispatch_data_t data=dispatch_data_create(alea_reprojection_metallib,alea_reprojection_metallib_len,NULL,DISPATCH_DATA_DESTRUCTOR_DEFAULT);
    id<MTLLibrary> lib=[g->device newLibraryWithData:data error:&error];dispatch_release(data);
    if(!lib){fprintf(stderr,"Shader library: %s\n",error.localizedDescription.UTF8String);return -1;}
    const char* names[]={"preview_clear","preview_project","preview_choose","preview_resolve"};
    for(int i=0;i<4;++i) {
        id<MTLFunction> f=[lib newFunctionWithName:[NSString stringWithUTF8String:names[i]]];
        g->pipeline[i]=[g->device newComputePipelineStateWithFunction:f error:&error];[f release];
        if(!g->pipeline[i]){[lib release];return -1;}
    }
    [lib release];
    g->points=[g->device newBufferWithLength:capacity*sizeof(preview_point) options:MTLResourceStorageModeShared];
    g->projected=[g->device newBufferWithLength:capacity*sizeof(preview_projection) options:MTLResourceStorageModeShared];
    g->depths=[g->device newBufferWithLength:capacity*4 options:MTLResourceStorageModeShared];
    g->owners=[g->device newBufferWithLength:capacity*4 options:MTLResourceStorageModeShared];
    g->color=[g->device newBufferWithLength:capacity*16 options:MTLResourceStorageModeShared];
    return g->queue && g->points && g->projected && g->depths && g->owners && g->color?0:-1;
}
static int project(preview_gpu* g,const preview_params* p,float* rgba,double* gpu_ms) {
    size_t n=(size_t)p->width*p->height;
    if(!p->width || !p->height || n>g->capacity || p->count>g->capacity ||
       !isfinite(p->half_width) || !isfinite(p->half_height) || p->half_width<=0 || p->half_height<=0)return -1;
    id<MTLCommandBuffer> command=[g->queue commandBuffer];if(!command)return -1;
    for(int pass=0;pass<4;++pass) {
        size_t count=(pass==1 || pass==2)?p->count:n;if(!count)continue;
        id<MTLComputeCommandEncoder> e=[command computeCommandEncoder];if(!e)return -1;
        [e setComputePipelineState:g->pipeline[pass]];
        if(pass==0){[e setBuffer:g->depths offset:0 atIndex:0];[e setBuffer:g->owners offset:0 atIndex:1];}
        if(pass==1) {
            [e setBuffer:g->points offset:0 atIndex:0];[e setBytes:p length:sizeof(*p) atIndex:1];
            [e setBuffer:g->projected offset:0 atIndex:2];[e setBuffer:g->depths offset:0 atIndex:3];
        }
        if(pass==2) {
            [e setBuffer:g->projected offset:0 atIndex:0];[e setBuffer:g->depths offset:0 atIndex:1];
            [e setBuffer:g->owners offset:0 atIndex:2];
        }
        if(pass==3) {
            [e setBuffer:g->points offset:0 atIndex:0];[e setBuffer:g->owners offset:0 atIndex:1];
            [e setBuffer:g->color offset:0 atIndex:2];[e setBytes:p length:sizeof(*p) atIndex:3];
        }
        NSUInteger group=MIN((NSUInteger)128,g->pipeline[pass].maxTotalThreadsPerThreadgroup);
        [e dispatchThreads:MTLSizeMake(count,1,1) threadsPerThreadgroup:MTLSizeMake(group,1,1)];[e endEncoding];
    }
    [command commit];[command waitUntilCompleted];
    if(command.status!=MTLCommandBufferStatusCompleted){fprintf(stderr,"Projection: %s\n",command.error.localizedDescription.UTF8String);return -1;}
    *gpu_ms=(command.GPUEndTime-command.GPUStartTime)*1000;
    memcpy(rgba,g->color.contents,n*16);return 0;
}
static preview_params params(const render_camera_t* c,int w,int h,size_t count) {
    preview_params p={.width=w,.height=h,.count=(preview_uint)count,.perspective=c->fov_degrees>0,
        .half_width=c->half_width,.half_height=c->half_height};
    for(int k=0;k<3;++k){p.eye[k]=c->eye[k];p.forward[k]=c->forward[k];p.right[k]=c->right[k];p.up[k]=c->true_up[k];}
    return p;
}
static size_t cache(preview_gpu* g,const render_camera_t* cam,const render_framebuffer_t* fb) {
    preview_point* points=g->points.contents;size_t n=(size_t)fb->width*fb->height,hits=0;
    memset(points,0,n*sizeof(*points));
    for(size_t i=0;i<n;++i)if(fb->material_id[i]>0 && isfinite(fb->depth[i]) && fb->depth[i]>0) {
        double o[3],d[3];render_camera_ray(cam,fb->width,fb->height,i%fb->width+.5,i/fb->width+.5,
            &o[0],&o[1],&o[2],&d[0],&d[1],&d[2]);
        for(int k=0;k<3;++k){points[i].world[k]=o[k]+fb->depth[i]*d[k];points[i].color[k]=fb->color[3*i+k];}
        points[i].valid=1;++hits;
    }
    return hits;
}
static int unit_tests(preview_gpu* g,float* rgba) {
    preview_point a[]={{{0,0,3},{1,0,0},1},{{0,0,2},{0,1,0},1},{{0,0,2},{0,0,1},1},
        {{0,0,-1},{1,1,1},1},{{NAN,0,1},{1,1,1},1},{{INFINITY,0,1},{1,1,1},1},
        {{100,0,1},{1,1,1},1},{{0,0,1},{1,1,1},0},{{0,0,0},{1,1,1},1}};
    memcpy(g->points.contents,a,sizeof(a));
    preview_params p={.forward={0,0,1},.right={1,0,0},.up={0,1,0},.half_width=1,.half_height=1,
        .width=4,.height=4,.perspective=1,.count=sizeof(a)/sizeof(*a)};
    double ms;
    for(int repeat=0;repeat<5;++repeat) {
        if(project(g,&p,rgba,&ms))return -1;
        const uint32_t* owners=g->owners.contents;
        for(int i=0;i<16;++i)if(owners[i]!=(i==10?1:PREVIEW_NONE) || rgba[4*i+3]!=(i==10?1:0))return -1;
        if(rgba[41]!=1 || rgba[40]!=0 || rgba[42]!=0)return -1;
    }
    p.count=0;if(project(g,&p,rgba,&ms))return -1;
    for(int i=0;i<16;++i)if(rgba[4*i+3]!=0 || ((uint32_t*)g->owners.contents)[i]!=PREVIEW_NONE)return -1;
    /* Analytically constructed pixel centers in a translated, rotated camera.
     * Non-multiple dispatch dimensions exercise the exact thread-grid edges. */
    p=(preview_params){.eye={3,-2,1},.forward={0,1,0},.right={1,0,0},.up={0,0,1},
        .half_width=1.4f,.half_height=1,.width=7,.height=5,.count=35};
    for(int perspective=0;perspective<2;++perspective) {
        p.perspective=perspective;preview_point* points=g->points.contents;
        for(int i=0;i<35;++i) {
            float z=2+.125f*i,scale=perspective?z:1;
            points[i]=(preview_point){.world={3+(2*((i%7)+.5f)/7-1)*p.half_width*scale,
                -2+z,1+(1-2*((i/7)+.5f)/5)*scale},.color={i/35.f,.25f,.75f},.valid=1};
        }
        if(project(g,&p,rgba,&ms))return -1;
        for(int i=0;i<35;++i)if(((uint32_t*)g->owners.contents)[i]!=(uint32_t)i ||
            rgba[4*i+3]!=1 || rgba[4*i]!=i/35.f ||
            fabs(((float*)g->depths.contents)[i]-(2+.125f*i))>1e-5)return -1;
    }
    p.width=(preview_uint)g->capacity+1;if(project(g,&p,rgba,&ms)!=-1)return -1;
    puts("tests=nearest_depth,tie,invalid,behind,offscreen,reset,empty,bounds,camera_modes PASS");return 0;
}
static void camera(render_config_t* c) {
    render_config_init(c);c->width=256;c->height=192;c->aa_samples=1;c->shadows=0;c->log_level=0;
    c->eye_set=c->target_set=1;c->eye[0]=16;c->eye[1]=-20;c->eye[2]=18;
    c->target[0]=c->target[1]=c->target[2]=0;
}
static int compare(preview_gpu* g,const char* name,const render_camera_t* cam,const render_framebuffer_t* truth,const float* rgba) {
    size_t hits=0,covered=0,ghost=0,bad=0;double error=0;
    const float* z=g->depths.contents;size_t n=(size_t)truth->width*truth->height;
    for(size_t i=0;i<n;++i) {
        int actual=truth->material_id[i]>0,known=rgba[4*i+3]>0;
        hits+=actual;ghost+=known && !actual;
        if(!known || !actual)continue;++covered;
        double o[3],d[3],depth=0;render_camera_ray(cam,truth->width,truth->height,i%truth->width+.5,i/truth->width+.5,
            &o[0],&o[1],&o[2],&d[0],&d[1],&d[2]);
        for(int k=0;k<3;++k)depth+=(o[k]+truth->depth[i]*d[k]-cam->eye[k])*cam->forward[k];
        double tolerance=2*cam->half_height*(cam->fov_degrees>0?depth:1)/truth->height;
        bad+=fabs(z[i]-depth)>tolerance;
        for(int k=0;k<3;++k)error+=fabs(rgba[4*i+k]-truth->color[3*i+k]);
    }
    printf("case=%s hits=%zu covered=%zu missing=%zu coverage_pct=%.3f ghosts=%zu depth_bad=%zu rgb_mae=%.6f\n",
        name,hits,covered,hits-covered,hits?100.*covered/hits:0,ghost,bad,covered?error/(3*covered):0);
    return (strcmp(name,"identity")==0 || strcmp(name,"ortho-identity")==0) && (covered!=hits || ghost || bad || error>1e-4)?-1:0;
}
static int picture(const char* path,const render_framebuffer_t* src,const float* rgba,const render_framebuffer_t* truth) {
    int w=src->width,h=src->height;size_t n=(size_t)w*h;
    render_framebuffer_t* panel=render_framebuffer_create(3*w,h,0);uint8_t* rgb=malloc(9*n);
    if(!panel || !rgb){render_framebuffer_free(panel);free(rgb);return -1;}
    for(int y=0;y<h;++y)for(int x=0;x<w;++x)for(int k=0;k<3;++k) {
        size_t i=(size_t)y*w+x,j=(size_t)y*3*w+x;
        panel->color[3*j+k]=src->color[3*i+k];panel->color[3*(j+w)+k]=rgba[4*i+k];
        panel->color[3*(j+2*w)+k]=truth->color[3*i+k];
    }
    render_tonemap(panel,rgb);int rc=render_write_image(path,rgb,3*w,h);free(rgb);render_framebuffer_free(panel);return rc;
}
int main(int argc,char** argv) {
    if(argc>2){fprintf(stderr,"Usage: %s [public-comparison.png]\n",argv[0]);return 2;}
    @autoreleasepool {
        (void)check_points;setbuf(stdout,NULL);int rc=1;const size_t n=256*192;preview_gpu g={0};
        float* rgba=malloc(n*16);alea_system_t* sys=scene();alea_compute_t* metal=NULL;
        render_framebuffer_t* src=render_framebuffer_create(256,192,1);
        render_framebuffer_t* truth=render_framebuffer_create(256,192,1);
        render_framebuffer_t* sparse=render_framebuffer_create(16,12,1);
        render_config_t cfg;camera(&cfg);
        double start=now();
        if(!rgba || !sys || !src || !truth || !sparse || create(&g,n))goto done;
        printf("projection_setup_ms=%.3f scratch_bound_mib=%.3f\n",(now()-start)*1000,n*(2*sizeof(preview_point)+sizeof(preview_projection)+40)/1048576.);
        if(unit_tests(&g,rgba))goto done;
        metal=alea_compute_create(sys,ALEA_COMPUTE_METAL,ALEA_COMPUTE_FP32);if(!metal)goto done;
        const char* names[]={"identity","orbit-small","pan","zoom","orbit-large","ortho-identity","ortho-pan","sparse-identity"};
        for(int view=0;view<8;++view) {
            render_config_free(&cfg);camera(&cfg);if(view==5 || view==6){cfg.fov=0;cfg.ortho_height=24;}
            render_framebuffer_t* source=view==7?sparse:src;
            cfg.width=source->width;cfg.height=source->height;
            render_camera_t source_cam,target_cam;if(render_camera_setup(&source_cam,&cfg,sys))goto done;
            start=now();if(alea_compute_render(metal,&cfg,&source_cam,source))goto done;
            printf("case=%s source_render_ms=%.3f\n",names[view],(now()-start)*1000);
            if(check_render(&cfg,&source_cam,source))goto done;
            start=now();size_t count=cache(&g,&source_cam,source);
            printf("case=%s cache_ms=%.3f cached_hits=%zu\n",names[view],(now()-start)*1000,count);
            if(view==1 || view==4) {
                double a=(view==1?2:20)*M_PI/180,x=cfg.eye[0],y=cfg.eye[1];
                cfg.eye[0]=x*cos(a)-y*sin(a);cfg.eye[1]=x*sin(a)+y*cos(a);
            }
            if(view==2 || view==6){cfg.eye[0]+=.25;cfg.target[0]+=.25;}
            if(view==3)for(int k=0;k<3;++k)cfg.eye[k]*=.95;
            cfg.width=256;cfg.height=192;
            if(render_camera_setup(&target_cam,&cfg,sys))goto done;
            preview_params p=params(&target_cam,cfg.width,cfg.height,(size_t)source->width*source->height);
            for(int pass=0;pass<11;++pass) {
                double gpu_ms;start=now();if(project(&g,&p,rgba,&gpu_ms))goto done;
                printf("case=%s pass=%d projection_ms=%.6f gpu_ms=%.6f\n",names[view],pass,(now()-start)*1000,gpu_ms);
            }
            if(alea_compute_render(metal,&cfg,&target_cam,truth) || check_render(&cfg,&target_cam,truth) || compare(&g,names[view],&target_cam,truth,rgba))goto done;
            if(view==1 && argc==2 && picture(argv[1],src,rgba,truth))goto done;
        }
        rc=0;puts("result=PASS");
    done:
        if(rc)fprintf(stderr,"Reprojection benchmark failed\n");
        render_config_free(&cfg);render_framebuffer_free(src);render_framebuffer_free(truth);render_framebuffer_free(sparse);
        alea_compute_destroy(metal);alea_destroy(sys);free(rgba);destroy(&g);return rc;
    }
}
