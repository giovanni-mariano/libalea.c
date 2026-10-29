// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
// SPDX-License-Identifier: MPL-2.0

#import <Foundation/Foundation.h>
#import <Metal/Metal.h>
#include "compute_internal.h"
#include "alea.h"
#include "compute_metallib.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <errno.h>

/* Manual ownership keeps Objective-C entirely inside this bridge. */
typedef struct {
    id<MTLDevice> device;
    id<MTLCommandQueue> queue;
    id<MTLComputePipelineState> pipeline;
    id<MTLBuffer> nodes, primitives, cells, points, results;
    uint32_t cell_count;
    int hierarchical;
    id<MTLBuffer> point_cells;
    size_t capacity;
    compute_scene scene;
    id<MTLComputePipelineState> render_pipeline;
    id<MTLBuffer> bvh, transforms, lattices, fills, colors, pixels;
    uint32_t bvh_root;
    size_t pixel_capacity, render_batch_limit;
    id<MTLComputePipelineState> resume_pipeline;
    id<MTLBuffer> resume_states, resume_done, resume_pixels;
    size_t resume_stride,resume_capacity;
    int refill,trace;
    uint32_t resume_work;
    int resume_fixed;
    uint32_t lazy_threshold;
    char* profile_prefix;
    id<MTLBuffer> ray_costs, cell_costs, spatial;
} metal_state;

int alea_compute_metal_available(void) {
    @autoreleasepool {
        id<MTLDevice> device = MTLCreateSystemDefaultDevice();
        int available = device != nil;
        [device release];
        return available;
    }
}

void alea_compute_metal_destroy(void* state) {
    if (!state) return;
    metal_state* s = state;
    [s->results release]; [s->points release]; [s->point_cells release];
    [s->cells release]; [s->nodes release]; [s->primitives release];
    [s->pipeline release]; [s->queue release]; [s->device release];
    [s->resume_pipeline release];[s->resume_states release];[s->resume_done release];[s->resume_pixels release];
    [s->render_pipeline release];
    [s->bvh release];
    [s->colors release]; [s->pixels release];
    [s->transforms release]; [s->lattices release]; [s->fills release];
    alea_compute_scene_free(&s->scene);
    [s->ray_costs release]; [s->cell_costs release];
    [s->spatial release];
    free(s->profile_prefix);
    free(s);
}

void* alea_compute_metal_create(const alea_system_t* sys) {
    @autoreleasepool {
        compute_scene scene;
        if (alea_compute_pack(sys, &scene)) return NULL;
        metal_state* s = calloc(1, sizeof(*s));
        if (!s) {
            alea_compute_scene_free(&scene);
            alea_set_error_detail(ALEA_ERR_OUT_OF_MEMORY, "Allocating Metal context");
            return NULL;
        }
        const char* profile=getenv("ALEA_METAL_PROFILE");
        if(profile && *profile) {
            s->profile_prefix=strdup(profile);
            if(!s->profile_prefix) {alea_set_error_detail(ALEA_ERR_OUT_OF_MEMORY,"Copying profile path");goto fail;}
        }
        const char* batch_env=getenv("ALEA_METAL_BATCH_PIXELS");
        if(batch_env) {
            char* end;errno=0;unsigned long long value=strtoull(batch_env,&end,10);
            if(errno || end==batch_env || *end || value==0 || value>UINT32_MAX) {
                alea_set_error_detail(ALEA_ERR_INVALID_ARG,"ALEA_METAL_BATCH_PIXELS must be an integer in 1..4294967295");goto fail;
            }
            s->render_batch_limit=(size_t)value;
        } else if(scene.node_count>1000000)s->render_batch_limit=512;
        s->resume_work=128;s->lazy_threshold=32;s->refill=1;
        s->trace=getenv("ALEA_BENCH_TRACE")!=NULL;
        const char* refill_env=getenv("ALEA_METAL_REFILL");
        if(refill_env) {
            if(!strcmp(refill_env,"1"))s->refill=1;
            else if(!strcmp(refill_env,"0"))s->refill=0;
            else {alea_set_error_detail(ALEA_ERR_INVALID_ARG,"ALEA_METAL_REFILL must be 0 or 1");goto fail;}
        }
        const char* lazy_env=getenv("ALEA_METAL_LAZY");
        if(lazy_env) {
            if(!strcmp(lazy_env,"1"))s->lazy_threshold=32;
            else if(!strcmp(lazy_env,"all"))s->lazy_threshold=1;
            else if(!strcmp(lazy_env,"0"))s->lazy_threshold=0;
            else {alea_set_error_detail(ALEA_ERR_INVALID_ARG,"ALEA_METAL_LAZY must be 0, 1, or all");goto fail;}
        }
        const char* work_env=getenv("ALEA_METAL_RESUME_WORK");
        if(work_env) {
            char* end;errno=0;unsigned long value=strtoul(work_env,&end,10);
            if(errno || end==work_env || *end || value<1 || value>65536) {
                alea_set_error_detail(ALEA_ERR_INVALID_ARG,"ALEA_METAL_RESUME_WORK must be an integer in 1..65536");goto fail;
            }
            s->resume_work=(uint32_t)value;s->resume_fixed=1;
        }
        s->device = MTLCreateSystemDefaultDevice();
        s->queue = [s->device newCommandQueue];
        NSError* error = nil;
        dispatch_data_t data = dispatch_data_create(alea_compute_metallib,
            alea_compute_metallib_len, NULL, DISPATCH_DATA_DESTRUCTOR_DEFAULT);
        id<MTLLibrary> library = [s->device newLibraryWithData:data error:&error];
        dispatch_release(data);
        s->hierarchical = scene.point_cells != NULL;
        id<MTLFunction> function = [library newFunctionWithName:
            s->hierarchical ? @"alea_classify_hierarchy" : @"alea_classify_points"];
        if (function) s->pipeline = [s->device newComputePipelineStateWithFunction:function error:&error];
        [function release]; [library release];
        if (!s->device || !s->queue || !s->pipeline) {
            alea_set_error_detail(ALEA_ERR_INVALID_STATE, "Creating Metal pipeline: %s",
                error ? [[error localizedDescription] UTF8String] : "device or pipeline unavailable");
            goto fail;
        }
        size_t node_bytes = (scene.node_count ? scene.node_count : 1) * sizeof(compute_node);
        size_t primitive_bytes = (scene.primitive_count ? scene.primitive_count : 1) * sizeof(compute_primitive);
        size_t cell_bytes = (scene.cell_count ? scene.cell_count : 1) * sizeof(compute_cell);
        if (node_bytes > s->device.maxBufferLength || primitive_bytes > s->device.maxBufferLength || cell_bytes > s->device.maxBufferLength) {
            alea_set_error_detail(ALEA_ERR_OVERFLOW, "Metal geometry exceeds maximum buffer size");
            goto fail;
        }
        s->nodes = [s->device newBufferWithBytes:scene.nodes length:node_bytes options:MTLResourceStorageModeShared];
        s->primitives = [s->device newBufferWithBytes:scene.primitives length:primitive_bytes options:MTLResourceStorageModeShared];
        s->cells = [s->device newBufferWithBytes:scene.cells length:cell_bytes options:MTLResourceStorageModeShared];
        s->cell_count = (uint32_t)scene.cell_count;
        if (!s->nodes || !s->primitives || !s->cells) {
            alea_set_error_detail(ALEA_ERR_OUT_OF_MEMORY, "Allocating Metal geometry buffers");
            goto fail;
        }
        free(scene.nodes); scene.nodes = NULL;
        free(scene.primitives); scene.primitives = NULL;
        s->scene = scene;
        return s;
    fail:
        alea_compute_scene_free(&scene);
        alea_compute_metal_destroy(s);
        return NULL;
    }
}

/* Classification and rendering share the transform upload. Query links are
 * separate because the render BVH intentionally excludes terminal void. */
static int prepare_classify(metal_state* s) {
    if (!s->hierarchical || s->point_cells) return 0;
    size_t bytes = s->scene.cell_count * sizeof(metal_point_cell);
    size_t tb = (s->scene.render.transform_count ? s->scene.render.transform_count : 1)
        * sizeof(compute_transform);
    if (bytes > s->device.maxBufferLength || (!s->transforms && tb > s->device.maxBufferLength)) {
        alea_set_error_detail(ALEA_ERR_OVERFLOW, "Metal point hierarchy exceeds buffer limit");
        return -1;
    }
    id<MTLBuffer> links = [s->device newBufferWithBytes:s->scene.point_cells
        length:bytes options:MTLResourceStorageModeShared];
    id<MTLBuffer> transforms = s->transforms ? [s->transforms retain] :
        [s->device newBufferWithBytes:s->scene.render.transforms length:tb options:MTLResourceStorageModeShared];
    if (!links || !transforms) {
        [links release]; [transforms release];
        alea_set_error_detail(ALEA_ERR_OUT_OF_MEMORY, "Uploading Metal point hierarchy");
        return -1;
    }
    s->point_cells = links;
    [s->transforms release]; s->transforms = transforms;
    free(s->scene.point_cells); s->scene.point_cells = NULL;
    return 0;
}

int alea_compute_metal_classify(void* state, const double* xyz, size_t count,
                                alea_compute_result_t* results) {
    @autoreleasepool {
        metal_state* s = state;
        if (prepare_classify(s)) return -1;
        const size_t total = count;
        const size_t limit = s->hierarchical ? 4096 : count;
        for (size_t first = 0; first < total; first += count) {
            if (alea_interrupted()) {
                alea_set_error_detail(ALEA_ERR_INTERRUPTED, "Metal batch interrupted"); return -1;
            }
            count = total - first < limit ? total - first : limit;
            if (count > s->device.maxBufferLength / (3 * sizeof(float)) ||
                count > s->device.maxBufferLength / sizeof(compute_result)) {
                alea_set_error_detail(ALEA_ERR_OVERFLOW, "Metal batch exceeds maximum buffer size");
                return -1;
            }
            if (count > s->capacity) {
                id<MTLBuffer> points = [s->device newBufferWithLength:count * 3 * sizeof(float) options:MTLResourceStorageModeShared];
                id<MTLBuffer> output = [s->device newBufferWithLength:count * sizeof(compute_result) options:MTLResourceStorageModeShared];
                if (!points || !output) {
                    [points release]; [output release];
                    alea_set_error_detail(ALEA_ERR_OUT_OF_MEMORY, "Allocating Metal batch buffers");
                    return -1;
                }
                [s->points release]; [s->results release];
                s->points = points; s->results = output; s->capacity = count;
            }
            float* points = s->points.contents;
            for (size_t i = 0; i < count * 3; ++i) points[i] = (float)xyz[3 * first + i];
            id<MTLCommandBuffer> command = [s->queue commandBuffer];
            id<MTLComputeCommandEncoder> encoder = [command computeCommandEncoder];
            if (!command || !encoder) {
                alea_set_error_detail(ALEA_ERR_INVALID_STATE, "Creating Metal command buffer");
                return -1;
            }
            [encoder setComputePipelineState:s->pipeline];
            [encoder setBuffer:s->nodes offset:0 atIndex:0];
            [encoder setBuffer:s->cells offset:0 atIndex:1];
            [encoder setBuffer:s->points offset:0 atIndex:2];
            [encoder setBuffer:s->results offset:0 atIndex:3];
            uint32_t point_count = (uint32_t)count;
            uint32_t root_or_count = s->hierarchical ? s->scene.point_root : s->cell_count;
            [encoder setBytes:&root_or_count length:sizeof(root_or_count) atIndex:4];
            [encoder setBytes:&point_count length:sizeof(point_count) atIndex:5];
            [encoder setBuffer:s->primitives offset:0 atIndex:6];
            if (s->hierarchical) {
                [encoder setBuffer:s->transforms offset:0 atIndex:7];
                [encoder setBuffer:s->point_cells offset:0 atIndex:8];
            }
            NSUInteger width = MIN((NSUInteger)128, s->pipeline.maxTotalThreadsPerThreadgroup);
            [encoder dispatchThreads:MTLSizeMake(count, 1, 1) threadsPerThreadgroup:MTLSizeMake(width, 1, 1)];
            [encoder endEncoding];
            [command commit];
            [command waitUntilCompleted];
            if (command.status != MTLCommandBufferStatusCompleted) {
                alea_set_error_detail(ALEA_ERR_INVALID_STATE, "Metal compute failed: %s",
                    command.error ? [[command.error localizedDescription] UTF8String] : "unknown command error");
                return -1;
            }
            if (alea_interrupted()) {
                alea_set_error_detail(ALEA_ERR_INTERRUPTED, "Metal batch interrupted");
                return -1;
            }
            const compute_result* output = s->results.contents;
            for (size_t i = 0; i < count; ++i) {
                if (output[i].owner_count == UINT32_MAX) {
                    alea_set_error_detail(ALEA_ERR_OVERFLOW,
                        "FP32 evaluation or hierarchy limit failure at point %zu", first + i);
                    return -1;
                }
            }
            memcpy(results + first, output, count * sizeof(*results));
        }
        return 0;
    }
}

/* Upload render-only resources lazily; transform buffers may already be shared. */
static int prepare_render(metal_state* s) {
    if (s->render_pipeline) return 0;
    compute_render_scene r = s->scene.render;
    NSError* error = nil;
    dispatch_data_t data = dispatch_data_create(alea_compute_metallib,
        alea_compute_metallib_len, NULL, DISPATCH_DATA_DESTRUCTOR_DEFAULT);
    id<MTLLibrary> library = [s->device newLibraryWithData:data error:&error];
    dispatch_release(data);
    MTLFunctionConstantValues* constants=[[MTLFunctionConstantValues alloc] init];
    bool profile=s->profile_prefix!=NULL;
    [constants setConstantValue:&profile type:MTLDataTypeBool atIndex:0];
    [constants setConstantValue:&s->lazy_threshold type:MTLDataTypeUInt atIndex:2];
    bool refill=s->refill;
    [constants setConstantValue:&refill type:MTLDataTypeBool atIndex:4];
    id<MTLFunction> function = [library newFunctionWithName:@"alea_render_image" constantValues:constants error:&error];
    [constants release];
    id<MTLComputePipelineState> pipeline = function ?
        [s->device newComputePipelineStateWithFunction:function error:&error] : nil;
    [function release]; [library release];
    if (!pipeline) {
        alea_set_error_detail(ALEA_ERR_INVALID_STATE, "Creating Metal renderer: %s",
            error ? [[error localizedDescription] UTF8String] : "pipeline unavailable");
        return -1;
    }
    size_t bb = (r.bvh_count ? r.bvh_count : 1)*sizeof(compute_render_bvh);
    size_t tb = (r.transform_count ? r.transform_count : 1)*sizeof(compute_transform);
    size_t ab = (r.lattice_count ? r.lattice_count : 1)*sizeof(compute_lattice);
    size_t fb = (r.fill_count ? r.fill_count : 1)*sizeof(uint32_t);
    size_t color_bytes = (s->cell_count ? s->cell_count : 1)*3*sizeof(float);
    size_t spatial_bytes=s->scene.spatial_count*sizeof(compute_spatial_program);
    if (bb>s->device.maxBufferLength || tb>s->device.maxBufferLength || ab>s->device.maxBufferLength || fb>s->device.maxBufferLength || color_bytes>s->device.maxBufferLength ||
        spatial_bytes>s->device.maxBufferLength) {
        alea_set_error_detail(ALEA_ERR_OVERFLOW, "Metal render scene exceeds buffer limit");
        [pipeline release]; return -1;
    }
    id<MTLBuffer> bvh = [s->device newBufferWithBytes:r.bvh length:bb options:MTLResourceStorageModeShared];
    id<MTLBuffer> colors = [s->device newBufferWithLength:color_bytes options:MTLResourceStorageModeShared];
    id<MTLBuffer> transforms = s->transforms ? [s->transforms retain] :
        [s->device newBufferWithBytes:r.transforms length:tb options:MTLResourceStorageModeShared];
    id<MTLBuffer> lattices = [s->device newBufferWithBytes:r.lattices length:ab options:MTLResourceStorageModeShared];
    id<MTLBuffer> fills = [s->device newBufferWithBytes:r.fills length:fb options:MTLResourceStorageModeShared];
    /* Commit ownership only after every allocation succeeds. A failed render
     * may be retried on this context, with its original host scene intact. */
    id<MTLBuffer> spatial=nil;
    if (!bvh || !colors || !transforms || !lattices || !fills) {
        alea_set_error_detail(ALEA_ERR_OUT_OF_MEMORY, "Allocating Metal render scene");goto fail;
    }
    spatial=spatial_bytes?[s->device newBufferWithBytes:s->scene.spatial length:spatial_bytes options:MTLResourceStorageModeShared]:
        [s->device newBufferWithLength:sizeof(compute_spatial_program) options:MTLResourceStorageModeShared];
    if(!spatial) {
        alea_set_error_detail(ALEA_ERR_OUT_OF_MEMORY,"Uploading spatial Metal programs");goto fail;
    }
    s->spatial=spatial;
    s->bvh_root=r.root;
    free(s->scene.spatial);s->scene.spatial=NULL;
    alea_compute_render_free(&s->scene.render);
    [s->transforms release]; s->transforms=transforms; s->lattices=lattices; s->fills=fills;
    s->render_pipeline=pipeline; s->bvh=bvh; s->colors=colors;
    if(s->trace)fprintf(stderr,"Metal device %s: %.2f MiB persistent render buffers\n",s->device.name.UTF8String,
        (s->nodes.length+s->primitives.length+s->cells.length+bvh.length+transforms.length+lattices.length+fills.length+colors.length+s->spatial.length)/1048576.0);
    return 0;
fail:
    [spatial release];
    [transforms release];[lattices release];[fills release];[pipeline release];[bvh release];[colors release];
    return -1;
}

/* Query the device's private continuation layout, including SIMD padding.
 * Scratch storage is bounded independently of image size, with up to 4096 slots. */
static int prepare_resume(metal_state* s) {
    if(s->resume_pipeline)return 0;
    NSError* error=nil;
    dispatch_data_t data=dispatch_data_create(alea_compute_metallib,alea_compute_metallib_len,NULL,DISPATCH_DATA_DESTRUCTOR_DEFAULT);
    id<MTLLibrary> library=[s->device newLibraryWithData:data error:&error];dispatch_release(data);
    MTLFunctionConstantValues* constants=[[MTLFunctionConstantValues alloc] init];
    bool profile=s->profile_prefix!=NULL;
    [constants setConstantValue:&profile type:MTLDataTypeBool atIndex:0];
    [constants setConstantValue:&s->lazy_threshold type:MTLDataTypeUInt atIndex:2];
    bool refill=s->refill;
    [constants setConstantValue:&refill type:MTLDataTypeBool atIndex:4];
    id<MTLFunction> function=[library newFunctionWithName:@"alea_render_resume" constantValues:constants error:&error];
    [constants release];
    id<MTLComputePipelineState> pipeline=function?[s->device newComputePipelineStateWithFunction:function error:&error]:nil;
    [function release];
    function=[library newFunctionWithName:@"alea_resume_layout"];
    id<MTLComputePipelineState> layout=function?[s->device newComputePipelineStateWithFunction:function error:&error]:nil;
    [function release];[library release];
    id<MTLBuffer> size=[s->device newBufferWithLength:sizeof(uint32_t) options:MTLResourceStorageModeShared];
    id<MTLCommandBuffer> command=[s->queue commandBuffer];
    id<MTLComputeCommandEncoder> encoder=command?[command computeCommandEncoder]:nil;
    int rc=-1;
    if(!pipeline || !layout || !size || !encoder) {
        if(encoder)[encoder endEncoding];
        alea_set_error_detail(ALEA_ERR_INVALID_STATE,"Preparing resumable Metal renderer: %s",error?error.localizedDescription.UTF8String:"allocation failed");goto done;
    }
    [encoder setComputePipelineState:layout];[encoder setBuffer:size offset:0 atIndex:0];
    [encoder dispatchThreads:MTLSizeMake(1,1,1) threadsPerThreadgroup:MTLSizeMake(1,1,1)];
    [encoder endEncoding];[command commit];[command waitUntilCompleted];
    if(command.status!=MTLCommandBufferStatusCompleted){alea_set_error_detail(ALEA_ERR_INVALID_STATE,"Reading Metal continuation layout");goto done;}
    size_t stride=*(uint32_t*)size.contents;
    size_t available=MIN((size_t)64*1024*1024,(size_t)s->device.maxBufferLength);
    size_t capacity=stride?MIN((size_t)4096,available/(stride+2*sizeof(uint32_t))):0;
    if(!capacity){alea_set_error_detail(ALEA_ERR_OVERFLOW,"Metal ray continuations exceed the 64 MiB scratch budget");goto done;}
    s->resume_states=[s->device newBufferWithLength:stride*capacity options:MTLResourceStorageModePrivate];
    s->resume_done=[s->device newBufferWithLength:capacity*sizeof(uint32_t) options:MTLResourceStorageModeShared];
    s->resume_pixels=[s->device newBufferWithLength:capacity*sizeof(uint32_t) options:MTLResourceStorageModeShared];
    if(!s->resume_states || !s->resume_done || !s->resume_pixels) {
        [s->resume_states release];s->resume_states=nil;[s->resume_done release];s->resume_done=nil;
        [s->resume_pixels release];s->resume_pixels=nil;
        alea_set_error_detail(ALEA_ERR_OUT_OF_MEMORY,"Allocating resumable Metal rays");goto done;
    }
    s->resume_stride=stride;s->resume_capacity=capacity;
    if(s->trace)fprintf(stderr,"Metal device %s: continuation %zu bytes, %zu slots, %.2f MiB scratch, lazy threshold %u, refill %d\n",
        s->device.name.UTF8String,stride,capacity,capacity*(stride+2*sizeof(uint32_t))/1048576.0,s->lazy_threshold,s->refill);
    s->resume_pipeline=pipeline;pipeline=nil;rc=0;
done:
    [pipeline release];[layout release];[size release];return rc;
}

/* Enumerate the same centre-first 64x64 tiles without an image-sized queue. */
static uint32_t next_pixel(compute_render_params* p,size_t* order,size_t* local) {
    size_t nx=(p->width+63)/64;
    size_t center=(p->height/2/64)*nx+p->width/2/64;
    size_t id=*order==0?center:*order==center?0:*order;
    size_t x=(id%nx)*64,y=(id/nx)*64,w=MIN(64u,p->width-x),h=MIN(64u,p->height-y);
    uint32_t pixel=(uint32_t)((y+*local/w)*p->width+x+*local%w);
    if(++*local==w*h){*local=0;++*order;}
    return pixel;
}
static int render_resumable(metal_state* s,compute_render_params* p,render_framebuffer_t* fb,
        alea_compute_render_progress_t progress,void* user) {
    if(prepare_resume(s))return -1;
    size_t total=(size_t)p->width*p->height,completed=0;
    if(progress(user,fb,0,total))return 1;
    size_t edge=MIN((size_t)64,(size_t)sqrt((double)s->resume_capacity));
    size_t nx=(p->width+edge-1)/edge,ny=(p->height+edge-1)/edge,center=(p->height/2/edge)*nx+p->width/2/edge;
    size_t order=0,local=0,issued=0,submissions=0,active_samples=0;double gpu_seconds=0;
    uint32_t* done=s->resume_done.contents;uint32_t* slot_pixels=s->resume_pixels.contents;
    for(size_t t=0;t<(s->refill?1:nx*ny);++t) {
        size_t tile_id=t==0?center:t==center?0:t;
        compute_render_tile tile={.x=(tile_id%nx)*edge,.y=(tile_id/nx)*edge};tile.width=MIN(edge,p->width-tile.x);
        size_t count=s->refill?MIN(total,s->resume_capacity):tile.width*MIN(edge,p->height-tile.y);
        size_t remaining=s->refill?total:count,active=count;
        if(s->refill)for(size_t j=0;j<count;++j){slot_pixels[j]=next_pixel(p,&order,&local);done[j]=3;++issued;}
        compute_resume_params dispatch={.work=s->resume_fixed?s->resume_work:MIN(s->resume_work,128u),.reset=1};
        while(remaining) {
            @autoreleasepool {
            if(alea_interrupted()){alea_set_error_detail(ALEA_ERR_INTERRUPTED,"Metal render interrupted");return -1;}
            id<MTLCommandBuffer> command=[s->queue commandBuffer];
            id<MTLComputeCommandEncoder> encoder=[command computeCommandEncoder];
            if(!command || !encoder){alea_set_error_detail(ALEA_ERR_INVALID_STATE,"Creating resumable Metal command");return -1;}
            [encoder setComputePipelineState:s->resume_pipeline];
            id<MTLBuffer> buffers[]={s->nodes,s->cells,s->bvh,s->colors,s->pixels,nil,s->primitives,s->transforms,s->lattices,s->fills};
            for(int i=0;i<10;++i)if(i!=5)[encoder setBuffer:buffers[i] offset:0 atIndex:i];
            [encoder setBytes:p length:sizeof(*p) atIndex:5];
            [encoder setBytes:&tile length:sizeof(tile) atIndex:10];
            [encoder setBuffer:s->spatial offset:0 atIndex:15];
            [encoder setBuffer:s->resume_pixels offset:0 atIndex:20];
            if(s->profile_prefix){[encoder setBuffer:s->ray_costs offset:0 atIndex:11];[encoder setBuffer:s->cell_costs offset:0 atIndex:12];}
            [encoder setBuffer:s->resume_states offset:0 atIndex:16];[encoder setBuffer:s->resume_done offset:0 atIndex:17];
            [encoder setBytes:&dispatch length:sizeof(dispatch) atIndex:18];
            NSUInteger width=MIN((NSUInteger)64,s->resume_pipeline.maxTotalThreadsPerThreadgroup);
            [encoder dispatchThreads:MTLSizeMake(count,1,1) threadsPerThreadgroup:MTLSizeMake(width,1,1)];
            [encoder endEncoding];[command commit];[command waitUntilCompleted];
            if(command.status!=MTLCommandBufferStatusCompleted) {
                alea_set_error_detail(ALEA_ERR_INVALID_STATE,"Resumable Metal render failed: %s",command.error?command.error.localizedDescription.UTF8String:"unknown error");return -1;
            }
            if(alea_interrupted()){alea_set_error_detail(ALEA_ERR_INTERRUPTED,"Metal render interrupted");return -1;}
            ++submissions;active_samples+=active;gpu_seconds+=command.GPUEndTime-command.GPUStartTime;
            const compute_render_pixel* pixels=s->pixels.contents;
            for(size_t j=0;j<count;++j)if(done[j]==1) {
                size_t i=s->refill?slot_pixels[j]:(tile.y+j/tile.width)*p->width+tile.x+j%tile.width;
                if(pixels[i].error){alea_set_error_detail(ALEA_ERR_OVERFLOW,"Resumable Metal ray failed at pixel %zu (code %d)",i,pixels[i].error);return -1;}
                memcpy(fb->color+3*i,pixels[i].color,3*sizeof(float));fb->cell_id[i]=pixels[i].cell_id;
                if(fb->material_id)fb->material_id[i]=pixels[i].material_id;
                if(fb->depth)fb->depth[i]=pixels[i].depth;
                if(fb->normal)memcpy(fb->normal+3*i,pixels[i].normal,3*sizeof(float));
                done[j]=2;--remaining;++completed;
                if(s->refill && issued<total){slot_pixels[j]=next_pixel(p,&order,&local);done[j]=3;++issued;}
                else --active;
            }
            if(!s->resume_fixed) {
                double elapsed=command.GPUEndTime-command.GPUStartTime;
                if(elapsed>.025)dispatch.work=MAX(1u,(uint32_t)(dispatch.work*.020/elapsed));
                else if(elapsed<.010)dispatch.work=MIN(16384u,dispatch.work*2);
                s->resume_work=dispatch.work;
            }
            dispatch.reset=0;
            // Call even without finished pixels: all active rays may be inside CSG.
            if(progress(user,fb,completed,total))return 1;
            }
        }
    }
    if(s->trace)fprintf(stderr,"Metal schedule: %zu submissions, %.3f GPU seconds, %.1f active slots/submission, %.2f GiB continuation traffic upper bound\n",
        submissions,gpu_seconds,submissions?(double)active_samples/submissions:0,2.0*s->resume_stride*active_samples/(1024.0*1024*1024));
    if(s->profile_prefix && alea_compute_profile_write(s->profile_prefix,s->ray_costs.contents,
        s->cell_costs.contents,s->scene.cells,s->cell_count,p->width,p->height))return -1;
    return 0;
}

int alea_compute_metal_render(void* state, const render_config_t* cfg,
                              const render_camera_t* cam, render_framebuffer_t* fb,
                              alea_compute_render_progress_t progress, void* user) {
    @autoreleasepool {
        metal_state* s=state;
        if (cfg->num_clips<0 || cfg->num_clips>1 || cfg->shadows || cfg->aa_adaptive ||
            cfg->render_mode<RENDER_MODE_SOLID || cfg->render_mode>RENDER_MODE_MATID ||
            cfg->render_mode==RENDER_MODE_XRAY ||
            (cfg->color_mode!=RENDER_COLOR_CELL && cfg->color_mode!=RENDER_COLOR_MATERIAL) ||
            cfg->material_filter.mode!=RENDER_FILTER_ALL || cfg->cell_filter.mode!=RENDER_FILTER_ALL ||
            cfg->density_min!=0 || cfg->density_max!=INFINITY ||
            cfg->aa_samples<1 || cfg->aa_samples>4) {
            alea_set_error_detail(ALEA_ERR_UNSUPPORTED,
                "Metal rendering supports solid/depth/ID images, material/cell colors, one clip, AA 1..4, no shadows or filters");
            return -1;
        }
        compute_render_params p={0};
        for (int a=0;a<3;++a) {
            p.eye[a]=cam->eye[a]; p.forward[a]=cam->forward[a];
            p.right[a]=cam->right[a]; p.up[a]=cam->true_up[a];
            p.light[a]=cfg->light_dir[a]; p.background[a]=cfg->background[a];
            if (cfg->num_clips) p.clip[a]=cfg->clips[0].normal[a];
        }
        if (cfg->num_clips) p.clip[3]=cfg->clips[0].d;
        p.half_width=cam->half_width; p.half_height=cam->half_height;
        p.ambient=cfg->ambient; p.diffuse=cfg->diffuse; p.specular=cfg->specular;
        p.shininess=cfg->shininess; p.tint=cfg->cross_section_tint;
        double clip_len=hypot(hypot(p.clip[0],p.clip[1]),p.clip[2]);
        if (cfg->num_clips && clip_len==0) {
            alea_set_error_detail(ALEA_ERR_INVALID_ARG,"Clip plane normal must be nonzero"); return -1;
        }
        if (cfg->num_clips) for (int a=0;a<4;++a) p.clip[a]/=clip_len;
        /* Copy the scalar float prefix instead of aliasing across struct members. */
        float values[offsetof(compute_render_params,width)/sizeof(float)];
        memcpy(values,&p,sizeof(values));
        for (size_t i=0;i<sizeof(values)/sizeof(values[0]);++i)
            if (!isfinite(values[i])) {
                alea_set_error_detail(ALEA_ERR_INVALID_ARG,"Render parameters must fit finite FP32"); return -1;
            }
        if(progress && progress(user,fb,0,(size_t)fb->width*fb->height))return 1;
        if (prepare_render(s)) return -1;
        p.width=fb->width; p.height=fb->height; p.perspective=cam->fov_degrees!=0;
        p.light_set=cfg->light_set; p.clip_enabled=cfg->num_clips; p.bvh_root=s->bvh_root; p.aa=cfg->aa_samples;
        size_t count=(size_t)fb->width*fb->height;
        if (count>s->device.maxBufferLength/sizeof(compute_render_pixel)) {
            alea_set_error_detail(ALEA_ERR_OVERFLOW,"Metal image exceeds buffer limit"); return -1;
        }
        if (count>s->pixel_capacity) {
            id<MTLBuffer> pixels=[s->device newBufferWithLength:count*sizeof(compute_render_pixel) options:MTLResourceStorageModeShared];
            if (!pixels) { alea_set_error_detail(ALEA_ERR_OUT_OF_MEMORY,"Allocating Metal image"); return -1; }
            [s->pixels release]; s->pixels=pixels; s->pixel_capacity=count;
        }
        if(s->profile_prefix) {
            size_t cb=(s->cell_count?s->cell_count:1)*sizeof(compute_profile);
            size_t rb=count*sizeof(compute_profile);
            if(cb>s->device.maxBufferLength || rb>s->device.maxBufferLength) {
                alea_set_error_detail(ALEA_ERR_OVERFLOW,"Metal profile exceeds buffer limit");return -1;
            }
            if(!s->cell_costs)s->cell_costs=[s->device newBufferWithLength:cb options:MTLResourceStorageModeShared];
            if(!s->ray_costs || s->ray_costs.length<rb) {
                [s->ray_costs release];s->ray_costs=[s->device newBufferWithLength:rb options:MTLResourceStorageModeShared];
            }
            if(!s->cell_costs || !s->ray_costs) {
                alea_set_error_detail(ALEA_ERR_OUT_OF_MEMORY,"Allocating Metal profile");return -1;
            }
            memset(s->cell_costs.contents,0,cb);
        }
        float* colors=s->colors.contents;
        for (uint32_t i=0;i<s->cell_count;++i) {
            int id=cfg->color_mode==RENDER_COLOR_CELL?s->scene.cells[i].cell_id:s->scene.cells[i].material_id;
            render_get_color(id,cfg->color_mode,cfg,&colors[3*i],&colors[3*i+1],&colors[3*i+2]);
        }
        if(progress)return render_resumable(s,&p,fb,progress,user);
        /* Ordinary renders retain the existing complete-ray kernel. */
        size_t batch_limit=s->render_batch_limit?s->render_batch_limit:count;
        compute_render_tile tile={.width=p.width};
        for(size_t first=0;first<count;) {
            @autoreleasepool {
            if (alea_interrupted()) { alea_set_error_detail(ALEA_ERR_INTERRUPTED,"Metal render interrupted"); return -1; }
            size_t batch=MIN(batch_limit,count-first);
            id<MTLCommandBuffer> command=[s->queue commandBuffer];
            id<MTLComputeCommandEncoder> encoder=[command computeCommandEncoder];
            if (!command || !encoder) {
                alea_set_error_detail(ALEA_ERR_INVALID_STATE,"Creating Metal render command"); return -1;
            }
            [encoder setComputePipelineState:s->render_pipeline];
            [encoder setBuffer:s->nodes offset:0 atIndex:0];
            [encoder setBuffer:s->cells offset:0 atIndex:1];
            [encoder setBuffer:s->bvh offset:0 atIndex:2];
            [encoder setBuffer:s->colors offset:0 atIndex:3];
            [encoder setBuffer:s->pixels offset:0 atIndex:4];
            [encoder setBytes:&p length:sizeof(p) atIndex:5];
            [encoder setBuffer:s->primitives offset:0 atIndex:6];
            [encoder setBuffer:s->transforms offset:0 atIndex:7];
            [encoder setBuffer:s->lattices offset:0 atIndex:8];
            [encoder setBuffer:s->fills offset:0 atIndex:9];
            tile.offset=(uint32_t)first;
            [encoder setBytes:&tile length:sizeof(tile) atIndex:10];
            [encoder setBuffer:s->spatial offset:0 atIndex:15];
            if(s->profile_prefix) {
                [encoder setBuffer:s->ray_costs offset:0 atIndex:11];
                [encoder setBuffer:s->cell_costs offset:0 atIndex:12];
            }
            NSUInteger width=MIN((NSUInteger)64,s->render_pipeline.maxTotalThreadsPerThreadgroup);
            [encoder dispatchThreads:MTLSizeMake(batch,1,1) threadsPerThreadgroup:MTLSizeMake(width,1,1)];
            [encoder endEncoding]; [command commit]; [command waitUntilCompleted];
            if (command.status!=MTLCommandBufferStatusCompleted) {
                alea_set_error_detail(ALEA_ERR_INVALID_STATE,"Metal render failed: %s",
                    command.error?[[command.error localizedDescription] UTF8String]:"unknown command error"); return -1;
            }
            if (alea_interrupted()) { alea_set_error_detail(ALEA_ERR_INTERRUPTED,"Metal render interrupted"); return -1; }
            double elapsed=command.GPUEndTime-command.GPUStartTime;
            if(elapsed>1 && batch_limit>128) {
                batch_limit=MAX((size_t)128,batch_limit/2);s->render_batch_limit=batch_limit;
            }
            const compute_render_pixel* pixels=s->pixels.contents;
            for(size_t j=first;j<first+batch;++j) {
                size_t i=(tile.y+j/tile.width)*p.width+tile.x+j%tile.width;
                if(pixels[i].error) {
                    alea_set_error_detail(ALEA_ERR_OVERFLOW,"Metal ray failed at pixel %zu (code %d: 1=arithmetic or interval budget, 2=traversal budget, 3=lattice precision, 4=hierarchy depth, 5=BVH stack)",i,pixels[i].error);return -1;
                }
                memcpy(fb->color+3*i,pixels[i].color,3*sizeof(float));
                fb->cell_id[i]=pixels[i].cell_id;
                if(fb->material_id)fb->material_id[i]=pixels[i].material_id;
                if(fb->depth)fb->depth[i]=pixels[i].depth;
                if(fb->normal)memcpy(fb->normal+3*i,pixels[i].normal,3*sizeof(float));
            }
            first+=batch;
            }
        }
        if(s->profile_prefix && alea_compute_profile_write(s->profile_prefix,
            s->ray_costs.contents,s->cell_costs.contents,s->scene.cells,s->cell_count,p.width,p.height))return -1;
        return 0;
    }
}
