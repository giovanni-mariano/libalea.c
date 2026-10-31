// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
// SPDX-License-Identifier: MPL-2.0
/* Display-only root bounds experiment; no CSG surface approximation is claimed. */
#include "../src/compute/compute_metal.m"
#include "gpu_geometry_fixture.h"
#include "bounds_preview_layout.h"
#include "bounds_preview_metallib.h"
_Static_assert(sizeof(preview_box)==24,"Box layout");
_Static_assert(sizeof(bounds_camera)==64,"Camera layout");

typedef struct {
    id<MTLRenderPipelineState> pipeline;
    id<MTLTexture> target;
    id<MTLBuffer> boxes,readback;
    size_t count,capacity;
    unsigned width,height;
} bounds_gpu;
static double seconds(void) {struct timespec t;clock_gettime(CLOCK_MONOTONIC,&t);return t.tv_sec+t.tv_nsec*1e-9;}
static void free_bounds(bounds_gpu* g) {
    [g->pipeline release];[g->target release];[g->boxes release];[g->readback release];memset(g,0,sizeof(*g));
}
static int init_bounds(bounds_gpu* g,id<MTLDevice> device,unsigned w,unsigned h,size_t capacity) {
    size_t bytes=(size_t)w*h*4;
    if(!w || w%64 || !h || w>1024 || h>1024 || capacity>BOUNDS_PREVIEW_LIMIT ||
        3*bytes+capacity*sizeof(preview_box)>64*1024*1024)return -1;
    g->width=w;g->height=h;g->capacity=capacity;
    NSError* error=nil;
    dispatch_data_t data=dispatch_data_create(alea_bounds_preview_metallib,alea_bounds_preview_metallib_len,NULL,DISPATCH_DATA_DESTRUCTOR_DEFAULT);
    id<MTLLibrary> lib=[device newLibraryWithData:data error:&error];dispatch_release(data);
    if(!lib)return -1;
    MTLRenderPipelineDescriptor* p=[[MTLRenderPipelineDescriptor alloc] init];
    p.vertexFunction=[lib newFunctionWithName:@"bounds_vertex_main"];
    p.fragmentFunction=[lib newFunctionWithName:@"bounds_fragment_main"];
    p.colorAttachments[0].pixelFormat=MTLPixelFormatRGBA8Unorm;
    g->pipeline=[device newRenderPipelineStateWithDescriptor:p error:&error];
    [p.vertexFunction release];[p.fragmentFunction release];[p release];[lib release];
    if(!g->pipeline){fprintf(stderr,"Bounds pipeline: %s\n",error.localizedDescription.UTF8String);return -1;}
    MTLTextureDescriptor* t=[MTLTextureDescriptor texture2DDescriptorWithPixelFormat:MTLPixelFormatRGBA8Unorm width:w height:h mipmapped:NO];
    t.usage=MTLTextureUsageRenderTarget;t.storageMode=MTLStorageModePrivate;
    g->target=[device newTextureWithDescriptor:t];
    g->boxes=[device newBufferWithLength:MAX(capacity,1)*sizeof(preview_box) options:MTLResourceStorageModeShared];
    g->readback=[device newBufferWithLength:bytes options:MTLResourceStorageModeShared];
    return g->target && g->boxes && g->readback?0:-1;
}
static int draw_bounds(bounds_gpu* g,id<MTLCommandQueue> queue,const render_camera_t* cam,int filled,uint8_t* rgba,double* gpu_ms) {
    if(g->count>g->capacity || !isfinite(cam->half_width) || !isfinite(cam->half_height) || cam->half_width<=0 || cam->half_height<=0)return -1;
    bounds_camera p={.half_width=cam->half_width,.half_height=cam->half_height,.perspective=cam->fov_degrees>0,.filled=filled};
    for(int k=0;k<3;++k){p.eye[k]=cam->eye[k];p.forward[k]=cam->forward[k];p.right[k]=cam->right[k];p.up[k]=cam->true_up[k];}
    id<MTLCommandBuffer> cmd=[queue commandBuffer];if(!cmd)return -1;
    MTLRenderPassDescriptor* pass=[MTLRenderPassDescriptor renderPassDescriptor];
    pass.colorAttachments[0].texture=g->target;pass.colorAttachments[0].loadAction=MTLLoadActionClear;
    pass.colorAttachments[0].storeAction=MTLStoreActionStore;pass.colorAttachments[0].clearColor=MTLClearColorMake(.035,.045,.06,0);
    id<MTLRenderCommandEncoder> e=[cmd renderCommandEncoderWithDescriptor:pass];if(!e)return -1;
    [e setRenderPipelineState:g->pipeline];[e setCullMode:MTLCullModeNone];
    [e setVertexBuffer:g->boxes offset:0 atIndex:0];[e setVertexBytes:&p length:sizeof(p) atIndex:1];
    if(g->count)[e drawPrimitives:filled?MTLPrimitiveTypeTriangle:MTLPrimitiveTypeLine vertexStart:0 vertexCount:filled?36:24 instanceCount:g->count];
    [e endEncoding];
    id<MTLBlitCommandEncoder> blit=[cmd blitCommandEncoder];if(!blit)return -1;
    [blit copyFromTexture:g->target sourceSlice:0 sourceLevel:0 sourceOrigin:MTLOriginMake(0,0,0)
        sourceSize:MTLSizeMake(g->width,g->height,1) toBuffer:g->readback destinationOffset:0
        destinationBytesPerRow:g->width*4 destinationBytesPerImage:g->width*g->height*4];[blit endEncoding];
    [cmd commit];[cmd waitUntilCompleted];
    if(cmd.status!=MTLCommandBufferStatusCompleted){fprintf(stderr,"Bounds draw: %s\n",cmd.error.localizedDescription.UTF8String);return -1;}
    *gpu_ms=(cmd.GPUEndTime-cmd.GPUStartTime)*1000;
    memcpy(rgba,g->readback.contents,(size_t)g->width*g->height*4);return 0;
}
/* Root BVH bounds are in world coordinates. Child-universe bounds are not. */
static int extract_bounds(const compute_render_scene* scene,bounds_gpu* g) {
    g->count=0;size_t leaves=0,omitted=0;
    if(scene->root==COMPUTE_NONE){puts("root_leaves=0 boxes=0 omitted=0");return 0;}
    if(scene->root>=scene->bvh_count)return -1;
    size_t end=scene->bvh[scene->root].escape;if(end>scene->bvh_count || end<=scene->root)return -1;
    preview_box* boxes=g->boxes.contents;
    for(size_t i=scene->root;i<end;++i)if(scene->bvh[i].cell!=COMPUTE_NONE) {
        ++leaves;const float* b=scene->bvh[i].bounds;int valid=1;
        for(int k=0;k<3;++k)if(!isfinite(b[2*k]) || !isfinite(b[2*k+1]) || b[2*k]>=b[2*k+1])valid=0;
        if(!valid){++omitted;continue;}
        if(g->count==g->capacity){fprintf(stderr,"Bounds count exceeds experiment capacity\n");return -1;}
        memcpy(boxes[g->count++].b,b,sizeof(boxes[0].b));
    }
    printf("root_leaves=%zu boxes=%zu omitted=%zu\n",leaves,g->count,omitted);return 0;
}
static int box_hit(const preview_box* box,const render_camera_t* cam,int w,int h,int x,int y) {
    double o[3],d[3];render_camera_ray(cam,w,h,x+.5,y+.5,&o[0],&o[1],&o[2],&d[0],&d[1],&d[2]);
    double cosine=0;for(int k=0;k<3;++k)cosine+=d[k]*cam->forward[k];
    double lo=BOUNDS_PREVIEW_NEAR/cosine,hi=cam->fov_degrees>0?INFINITY:BOUNDS_PREVIEW_FAR;
    for(int k=0;k<3;++k) {
        if(fabs(d[k])<1e-20){if(o[k]<box->b[2*k] || o[k]>box->b[2*k+1])return 0;continue;}
        double a=(box->b[2*k]-o[k])/d[k],b=(box->b[2*k+1]-o[k])/d[k];
        lo=fmax(lo,fmin(a,b));hi=fmin(hi,fmax(a,b));
    }
    return hi>=lo;
}
static int test_bounds(metal_state* s) {
    bounds_gpu g={0};int rc=-1;uint8_t rgba[64*48*4];double gpu;
    if(init_bounds(&g,s->device,64,48,4))goto done;
    render_config_t cfg;render_config_init(&cfg);cfg.width=64;cfg.height=48;cfg.eye_set=cfg.target_set=1;
    cfg.eye[0]=cfg.eye[1]=cfg.eye[2]=0;cfg.target[0]=cfg.target[1]=0;cfg.target[2]=1;cfg.up[0]=cfg.up[2]=0;cfg.up[1]=1;
    const preview_box cases[]={{{-.37f,.71f,-.43f,.68f,2,3}},{{-.4f,.4f,-.4f,.4f,-3,-1}},
        {{-.47f,.61f,-.38f,.57f,-1,2}},{{1,2,1,2,2,3}}};
    for(int mode=0;mode<2;++mode) {
        cfg.fov=mode?0:60;cfg.ortho_height=4;render_camera_t cam;if(render_camera_setup(&cam,&cfg,NULL))goto done;
        for(size_t c=0;c<sizeof(cases)/sizeof(*cases);++c) {
            memcpy(g.boxes.contents,&cases[c],sizeof(preview_box));g.count=1;
            if(draw_bounds(&g,s->queue,&cam,1,rgba,&gpu))goto done;
            for(int y=0;y<48;++y)for(int x=0;x<64;++x)
                if((rgba[4*(y*64+x)+3]!=0)!=box_hit(&cases[c],&cam,64,48,x,y)) {
                    fprintf(stderr,"Bounds raster oracle mismatch mode=%d case=%zu\n",mode,c);goto done;
                }
        }
        g.count=0;if(draw_bounds(&g,s->queue,&cam,1,rgba,&gpu))goto done;
        for(size_t i=0;i<64*48;++i)if(rgba[4*i+3])goto done;
        g.count=5;if(draw_bounds(&g,s->queue,&cam,1,rgba,&gpu)!=-1)goto done;
    }
    /* Extraction visits only the root tree, excludes invalid bounds and fails
     * rather than silently truncating a capacity overflow. */
    compute_render_bvh nodes[4]={0};
    nodes[0].cell=COMPUTE_NONE;nodes[0].escape=3;
    nodes[1].cell=0;nodes[1].escape=2;memcpy(nodes[1].bounds,cases[0].b,24);
    nodes[2].cell=1;nodes[2].escape=3;memcpy(nodes[2].bounds,cases[0].b,24);nodes[2].bounds[0]=-INFINITY;
    nodes[3].cell=2;nodes[3].escape=4;memcpy(nodes[3].bounds,cases[0].b,24);
    compute_render_scene scene={.bvh=nodes,.bvh_count=4,.root=0};
    if(extract_bounds(&scene,&g) || g.count!=1)goto done;
    scene.root=COMPUTE_NONE;if(extract_bounds(&scene,&g) || g.count)goto done;
    scene.root=0;g.capacity=0;if(extract_bounds(&scene,&g)!=-1)goto done;
    puts("tests=perspective,orthographic,behind,near_clip,empty,reset,capacity,root_scope,invalid_bounds PASS");rc=0;
 done:free_bounds(&g);return rc;
}
static int save_panel(const char* path,const uint8_t* wire,const uint8_t* filled,const render_framebuffer_t* truth,int w,int h) {
    size_t n=(size_t)w*h;uint8_t* rgb=malloc(n*9);uint8_t* exact=malloc(n*3);
    if(!rgb || !exact){free(rgb);free(exact);return -1;}
    if(truth)render_tonemap(truth,exact);else memset(exact,64,n*3);
    for(int y=0;y<h;++y)for(int x=0;x<w;++x)for(int k=0;k<3;++k) {
        size_t i=(size_t)y*w+x,j=(size_t)y*3*w+x;
        rgb[3*j+k]=wire[4*i+k];rgb[3*(j+w)+k]=filled[4*i+k];rgb[3*(j+2*w)+k]=exact[3*i+k];
    }
    int rc=render_write_image(path,rgb,w*3,h);free(rgb);free(exact);return rc;
}
static int evaluate_bounds(metal_state* s,alea_system_t* sys,render_config_t* cfg,const char* prefix,int analytic) {
    int rc=-1;bounds_gpu g={0};size_t n=(size_t)cfg->width*cfg->height;
    uint8_t* wire=malloc(n*4);uint8_t* filled=malloc(n*4);render_framebuffer_t* truth=NULL;
    double evaluation_start=seconds(),start=evaluation_start;if(!wire || !filled || init_bounds(&g,s->device,cfg->width,cfg->height,BOUNDS_PREVIEW_LIMIT))goto done;
    printf("raster_setup_ms=%.6f\n",(seconds()-start)*1000);start=seconds();
    if(extract_bounds(&s->scene.render,&g))goto done;
    printf("extract_upload_ms=%.6f box_bytes=%zu preview_storage_bound_mib=%.3f\n",(seconds()-start)*1000,g.count*sizeof(preview_box),
        (BOUNDS_PREVIEW_LIMIT*sizeof(preview_box)+4*n*4)/1048576.);
    if(analytic){truth=render_framebuffer_create(cfg->width,cfg->height,1);if(!truth)goto done;}
    const char* names[]={"exterior","orbit","orthographic"};double original[3];memcpy(original,cfg->eye,sizeof(original));
    for(int view=0;view<3;++view) {
        if(view==1) {
            double a=20*M_PI/180,x=original[0]-cfg->target[0],y=original[1]-cfg->target[1];
            cfg->eye[0]=cfg->target[0]+x*cos(a)-y*sin(a);cfg->eye[1]=cfg->target[1]+x*sin(a)+y*cos(a);
        }
        if(view==2){memcpy(cfg->eye,original,sizeof(original));double dist=0;for(int k=0;k<3;++k)dist+=pow(cfg->eye[k]-cfg->target[k],2);cfg->ortho_height=sqrt(dist)*.8;cfg->fov=0;}
        render_camera_t cam;if(render_camera_setup(&cam,cfg,sys))goto done;
        for(int pass=0;pass<11;++pass) {
            double gpu;start=seconds();if(draw_bounds(&g,s->queue,&cam,0,wire,&gpu))goto done;
            printf("case=%s pass=%d wire_ms=%.6f gpu_ms=%.6f\n",names[view],pass,(seconds()-start)*1000,gpu);
            if(view==0 && pass==0)printf("prepared_scene_to_first_readback_ms=%.6f\n",(seconds()-evaluation_start)*1000);
        }
        double gpu;if(draw_bounds(&g,s->queue,&cam,1,filled,&gpu))goto done;
        size_t pixels=0,lines=0,hits=0,covered=0;
        if(analytic && (alea_compute_metal_render(s,cfg,&cam,truth,NULL,NULL) || check_render(cfg,&cam,truth)))goto done;
        for(size_t i=0;i<n;++i){pixels+=filled[4*i+3]!=0;lines+=wire[4*i+3]!=0;if(truth && truth->material_id[i]>0){++hits;covered+=filled[4*i+3]!=0;}}
        printf("case=%s wire_pixels=%zu bounds_pixels=%zu total_pixels=%zu",names[view],lines,pixels,n);
        if(analytic)printf(" analytic_hits=%zu retained=%zu excess=%zu hit_coverage_pct=%.3f bounds_precision_pct=%.3f",hits,covered,pixels-covered,hits?100.*covered/hits:0,pixels?100.*covered/pixels:0);
        puts("");
        if(analytic && covered!=hits){fprintf(stderr,"Public bounds lost sphere hits\n");goto done;}
        if(prefix){char path[1024];int count=snprintf(path,sizeof(path),"%s-%s.png",prefix,names[view]);if(count<0 || count>=(int)sizeof(path) || save_panel(path,wire,filled,truth,cfg->width,cfg->height))goto done;}
    }
    rc=0;
 done:free_bounds(&g);free(wire);free(filled);render_framebuffer_free(truth);return rc;
}
int main(int argc,char** argv) {
    if(argc>2){fprintf(stderr,"Usage: %s [public-image-prefix]\n",argv[0]);return 2;}
    @autoreleasepool {
        (void)check_points;setbuf(stdout,NULL);int rc=1;alea_system_t* sys=scene();metal_state* s=NULL;
        render_config_t cfg;render_config_init(&cfg);cfg.width=512;cfg.height=384;cfg.aa_samples=1;cfg.log_level=0;cfg.shadows=0;
        cfg.eye_set=cfg.target_set=1;cfg.eye[0]=16;cfg.eye[1]=-20;cfg.eye[2]=18;cfg.target[0]=cfg.target[1]=cfg.target[2]=0;
        double start=seconds();if(!sys || !(s=alea_compute_metal_create(sys)))goto done;
        printf("full_scene_prepare_ms=%.6f\n",(seconds()-start)*1000);
        if(evaluate_bounds(s,sys,&cfg,argc==2?argv[1]:NULL,1) || test_bounds(s))goto done;
        puts("result=PASS");rc=0;
 done:if(rc)fprintf(stderr,"Bounds preview experiment failed\n");render_config_free(&cfg);alea_compute_metal_destroy(s);alea_destroy(sys);return rc;
    }
}
