// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
// SPDX-License-Identifier: MPL-2.0
#include "alea_test.h"
#include <unistd.h>
/* Exercise private upload ownership without adding failure hooks to the API.
 * This translation unit supplies the bridge symbols instead of the archive. */
#include "../../src/compute/compute_metal.m"

static NSUInteger live_buffers;
@interface UploadBuffer : NSObject
@property(nonatomic) NSUInteger length;
@end
@implementation UploadBuffer
- (instancetype)init {if((self=[super init]))++live_buffers;return self;}
- (void)dealloc {--live_buffers;[super dealloc];}
@end

/* Pipeline creation uses the real device. Upload buffers are counted stand-ins;
 * most tests never submit them. The candidate retry test forwards real buffers. */
@interface UploadDevice : NSProxy {
@public
    id<MTLDevice> target;
    NSUInteger calls,fail_at;
    BOOL real_buffers;
}
@end
@implementation UploadDevice
- (NSMethodSignature*)methodSignatureForSelector:(SEL)selector {
    return [(NSObject*)target methodSignatureForSelector:selector];
}
- (void)forwardInvocation:(NSInvocation*)invocation {[invocation invokeWithTarget:target];}
- (id<MTLBuffer>)newBufferWithLength:(NSUInteger)length options:(MTLResourceOptions)options {
    if(++calls==fail_at)return nil;
    if(real_buffers)return [target newBufferWithLength:length options:options];
    UploadBuffer* buffer=[[UploadBuffer alloc] init];buffer.length=length;
    return (id<MTLBuffer>)buffer;
}
- (id<MTLBuffer>)newBufferWithBytes:(const void*)bytes length:(NSUInteger)length options:(MTLResourceOptions)options {
    (void)bytes;return [self newBufferWithLength:length options:options];
}
- (void)dealloc {[target release];[super dealloc];}
@end

TEST(metal_compact_storage_and_independent_context_formats) {
    @autoreleasepool {
        const char* flag=getenv("ALEA_METAL_COMPACT_CODE");char* saved=flag?strdup(flag):NULL;
        const char* candidates=getenv("ALEA_METAL_POINT_CANDIDATES");char* saved_candidates=candidates?strdup(candidates):NULL;
        ASSERT_EQ(alea_setenv("ALEA_METAL_POINT_CANDIDATES","1",1),0);
        alea_system_t* sys=alea_create();int mat=alea_add_material(sys,7);
        int box=alea_box_surface(sys,1,-2,2,-2,2,-2,2);
        int sphere=alea_sphere_surface(sys,2,0,0,0,1);
        alea_node_id_t root=alea_intersection(sys,alea_halfspace(sys,box,-1),alea_halfspace(sys,sphere,1));
        alea_add_cell(sys,42,root,mat,-1,0);
        metal_state* contexts[2];
        for(unsigned format=0;format<2;++format) {
            ASSERT_EQ(alea_setenv("ALEA_METAL_COMPACT_CODE",format?"1":"0",1),0);
            contexts[format]=alea_compute_metal_create(sys);ASSERT_NOT_NULL(contexts[format]);
            ASSERT_EQ(contexts[format]->compact_code,format);
        }
        ASSERT_EQ(contexts[0]->nodes.length,2*contexts[1]->nodes.length);
        /* Switch formats between contexts after creation and use the real host
         * proof as well as both GPU point pipelines on the same immutable code. */
        double xyz[]={0,0,0, 1.5,0,0, 3,0,0};
        for(unsigned pass=0;pass<2;++pass)for(unsigned format=0;format<2;++format) {
            metal_state* s=contexts[format];
            alea_compute_result_t classified[3];alea_compute_selection_t selected[3];
            ASSERT_EQ(alea_compute_metal_classify(s,xyz,3,classified),0);
            ASSERT_EQ(alea_compute_metal_select(s,xyz,3,selected),0);
            for(unsigned i=0;i<3;++i) {
                ASSERT_EQ(classified[i].owner_count,i==1);ASSERT_EQ(selected[i].status,i==1?ALEA_COMPUTE_FOUND:ALEA_COMPUTE_NOT_FOUND);
                ASSERT_EQ(selected[i].cell_id,i==1?42:-1);
            }
        }
        for(unsigned format=0;format<2;++format)alea_compute_metal_destroy(contexts[format]);
        const char* thread_value=getenv("ALEA_METAL_JUMP_THREADING");
        char* saved_thread=thread_value?strdup(thread_value):NULL;
        ASSERT_EQ(alea_setenv("ALEA_METAL_JUMP_THREADING","invalid",1),0);
        ASSERT_NULL(alea_compute_metal_create(sys));ASSERT_EQ(alea_error_code(),ALEA_ERR_INVALID_ARG);
        ASSERT_EQ(saved_thread?alea_setenv("ALEA_METAL_JUMP_THREADING",saved_thread,1):alea_unsetenv("ALEA_METAL_JUMP_THREADING"),0);
        free(saved_thread);
        ASSERT_EQ(alea_setenv("ALEA_METAL_COMPACT_CODE","invalid",1),0);
        ASSERT_NULL(alea_compute_metal_create(sys));ASSERT_EQ(alea_error_code(),ALEA_ERR_INVALID_ARG);
        ASSERT_EQ(alea_setenv("ALEA_METAL_COMPACT_CODE","1",1),0);
        metal_state* retry=alea_compute_metal_create(sys);ASSERT_NOT_NULL(retry);
        alea_compute_metal_destroy(retry);alea_destroy(sys);
        ASSERT_EQ(saved?alea_setenv("ALEA_METAL_COMPACT_CODE",saved,1):alea_unsetenv("ALEA_METAL_COMPACT_CODE"),0);
        ASSERT_EQ(saved_candidates?alea_setenv("ALEA_METAL_POINT_CANDIDATES",saved_candidates,1):alea_unsetenv("ALEA_METAL_POINT_CANDIDATES"),0);
        free(saved);free(saved_candidates);
    }
}

TEST(metal_upload_failure_releases_buffers_and_retries) {
    @autoreleasepool {
        alea_system_t* sys=alea_create();int mat=alea_add_material(sys,7);
        alea_node_id_t root=ALEA_NODE_ID_INVALID;
        for(int i=0;i<64;++i) {
            int sphere=alea_sphere_surface(sys,0,2*(i%4)-3,2*(i/4%4)-3,2*(i/16)-3,.5);
            alea_node_id_t leaf=alea_halfspace(sys,sphere,-1);
            root=i?alea_union(sys,root,leaf):leaf;
        }
        alea_add_cell(sys,42,root,mat,-1,0);
        /* Five main buffers, spatial programs. */
        for(NSUInteger failure=1;failure<=6;++failure) {
            metal_state* s=calloc(1,sizeof(*s));ASSERT_NOT_NULL(s);
            ASSERT_EQ(alea_compute_pack(sys,&s->scene),0);
            s->cell_count=(uint32_t)s->scene.cell_count;
            UploadDevice* device=[UploadDevice alloc];device->target=MTLCreateSystemDefaultDevice();
            ASSERT_NOT_NULL(device->target);device->fail_at=failure;s->device=(id<MTLDevice>)device;
            compute_render_bvh* original=s->scene.render.bvh;
            ASSERT_EQ(prepare_render(s),-1);ASSERT_EQ(alea_error_code(),ALEA_ERR_OUT_OF_MEMORY);
            ASSERT_EQ(live_buffers,0);
            ASSERT_NULL(s->render_pipeline);
            ASSERT_NULL(s->spatial);
            ASSERT_TRUE(s->scene.render.bvh==original);
            device->fail_at=0;
            ASSERT_EQ(prepare_render(s),0);ASSERT_EQ(live_buffers,6);
            ASSERT_NULL(s->scene.render.bvh);
            NSUInteger calls=device->calls;
            ASSERT_EQ(prepare_render(s),0);ASSERT_EQ(device->calls,calls);
            alea_compute_metal_destroy(s);ASSERT_EQ(live_buffers,0);
        }
        alea_destroy(sys);
    }
}
TEST(metal_point_upload_failure_and_shared_transform_ownership) {
    @autoreleasepool {
        alea_system_t* sys=alea_create();int mat=alea_add_material(sys,7);
        int sphere=alea_sphere_surface(sys,1,0,0,0,2);
        int parent=alea_add_cell(sys,10,alea_halfspace(sys,sphere,-1),ALEA_MATERIAL_VOID,0,0);
        ASSERT(alea_add_cell(sys,42,alea_halfspace(sys,sphere,-1),mat,-1,1)>=0);
        ASSERT_EQ(alea_set_fill(sys,parent,1,0),0);
        for(NSUInteger failure=1;failure<=9;++failure) {
            metal_state* s=calloc(1,sizeof(*s));ASSERT_NOT_NULL(s);
            ASSERT_EQ(alea_compute_pack(sys,&s->scene),0);
            s->hierarchical=1;s->cell_count=(uint32_t)s->scene.cell_count;
            UploadDevice* device=[UploadDevice alloc];device->target=MTLCreateSystemDefaultDevice();
            ASSERT_NOT_NULL(device->target);s->device=(id<MTLDevice>)device;
            metal_point_cell* original=s->scene.point_cells;
            if(failure<=5) {
                device->fail_at=failure;
                ASSERT_EQ(prepare_classify(s),-1);ASSERT_EQ(live_buffers,0);
                ASSERT_NULL(s->point_cells);ASSERT_NULL(s->transforms);
                ASSERT_TRUE(s->scene.point_cells==original);
                device->fail_at=0;
            }
            ASSERT_EQ(prepare_classify(s),0);ASSERT_EQ(live_buffers,5);
            ASSERT_NULL(s->scene.point_cells);
            id<MTLBuffer> transforms=s->transforms;
            NSUInteger calls=device->calls;
            ASSERT_EQ(prepare_classify(s),0);ASSERT_EQ(device->calls,calls);
            if(failure>5) {
                device->fail_at=device->calls+failure-5;
                ASSERT_EQ(prepare_render(s),-1);ASSERT_EQ(live_buffers,5);
                ASSERT_TRUE(s->transforms==transforms);ASSERT_NOT_NULL(s->scene.render.bvh);
                device->fail_at=0;
            }
            ASSERT_EQ(prepare_render(s),0);ASSERT_EQ(live_buffers,9);
            ASSERT_TRUE(s->transforms==transforms);
            calls=device->calls;ASSERT_EQ(prepare_classify(s),0);ASSERT_EQ(device->calls,calls);
            alea_compute_metal_destroy(s);ASSERT_EQ(live_buffers,0);
        }
        alea_destroy(sys);
    }
}

TEST(metal_point_upload_after_render_failure_preserves_shared_buffers) {
    @autoreleasepool {
        alea_system_t* sys=alea_create();int mat=alea_add_material(sys,7);
        int sphere=alea_sphere_surface(sys,1,0,0,0,2);
        int parent=alea_add_cell(sys,10,alea_halfspace(sys,sphere,-1),ALEA_MATERIAL_VOID,0,0);
        ASSERT(alea_add_cell(sys,42,alea_halfspace(sys,sphere,-1),mat,-1,1)>=0);
        ASSERT_EQ(alea_set_fill(sys,parent,1,0),0);
        for(NSUInteger failure=1;failure<=3;++failure) {
            metal_state* s=calloc(1,sizeof(*s));ASSERT_NOT_NULL(s);
            ASSERT_EQ(alea_compute_pack(sys,&s->scene),0);
            s->hierarchical=1;s->cell_count=(uint32_t)s->scene.cell_count;
            UploadDevice* device=[UploadDevice alloc];device->target=MTLCreateSystemDefaultDevice();
            ASSERT_NOT_NULL(device->target);s->device=(id<MTLDevice>)device;
            ASSERT_EQ(prepare_render(s),0);ASSERT_EQ(live_buffers,6);
            id<MTLBuffer> transforms=s->transforms;
            metal_point_cell* original=s->scene.point_cells;
            device->fail_at=device->calls+failure;
            ASSERT_EQ(prepare_classify(s),-1);ASSERT_EQ(live_buffers,6);
            ASSERT_TRUE(s->transforms==transforms);ASSERT_NULL(s->point_cells);
            ASSERT_TRUE(s->scene.point_cells==original);ASSERT_NOT_NULL(s->scene.point_lattices);
            device->fail_at=0;
            ASSERT_EQ(prepare_classify(s),0);ASSERT_EQ(live_buffers,9);
            ASSERT_TRUE(s->transforms==transforms);
            alea_compute_metal_destroy(s);ASSERT_EQ(live_buffers,0);
        }
        alea_destroy(sys);
    }
}

TEST(metal_point_candidate_allocation_failure_retries) {
    @autoreleasepool {
        alea_system_t* sys=alea_create();int mat=alea_add_material(sys,7);
        int sphere=alea_sphere_surface(sys,1,0,0,0,2);
        ASSERT(alea_add_cell(sys,42,alea_halfspace(sys,sphere,-1),mat,-1,0)>=0);
        ASSERT_EQ(alea_setenv("ALEA_METAL_POINT_CANDIDATES","1",1),0);
        metal_state* s=alea_compute_metal_create(sys);ASSERT_NOT_NULL(s);
        ASSERT_EQ(alea_unsetenv("ALEA_METAL_POINT_CANDIDATES"),0);alea_destroy(sys);
        UploadDevice* device=[UploadDevice alloc];device->target=s->device;
        device->real_buffers=YES;device->fail_at=1;s->device=(id<MTLDevice>)device;
        double point[]={0,0,0};alea_compute_result_t result;
        ASSERT_EQ(alea_compute_metal_classify(s,point,1,&result),-1);
        ASSERT_EQ(alea_error_code(),ALEA_ERR_OUT_OF_MEMORY);
        ASSERT_NULL(s->point_candidates);ASSERT_NULL(s->point_signs);
        device->fail_at=0;
        ASSERT_EQ(alea_compute_metal_classify(s,point,1,&result),0);ASSERT_EQ(result.cell_id,42);
        alea_compute_metal_destroy(s);
    }
}
TEST(metal_selected_counter_allocation_failure_retries) {
    @autoreleasepool {
        char path[160];snprintf(path,sizeof(path),"/tmp/alea-selected-resource-%ld.jsonl",(long)getpid());unlink(path);
        alea_system_t* sys=alea_create();int mat=alea_add_material(sys,7);
        int sphere=alea_sphere_surface(sys,1,0,0,0,2);
        ASSERT(alea_add_cell(sys,42,alea_halfspace(sys,sphere,-1),mat,-1,0)>=0);
        ASSERT_EQ(alea_setenv("ALEA_METAL_SELECT_PROFILE",path,1),0);
        metal_state* s=alea_compute_metal_create(sys);ASSERT_NOT_NULL(s);
        ASSERT_EQ(alea_unsetenv("ALEA_METAL_SELECT_PROFILE"),0);alea_destroy(sys);
        UploadDevice* device=[UploadDevice alloc];device->target=s->device;
        device->real_buffers=YES;device->fail_at=1;s->device=(id<MTLDevice>)device;
        double point[]={0,0,0};alea_compute_selection_t result;
        ASSERT_EQ(alea_compute_metal_select(s,point,1,&result),-1);
        ASSERT_EQ(alea_error_code(),ALEA_ERR_OUT_OF_MEMORY);ASSERT_NULL(s->select_costs);
        device->fail_at=0;ASSERT_EQ(alea_compute_metal_select(s,point,1,&result),0);
        ASSERT_EQ(result.cell_id,42);alea_compute_metal_destroy(s);unlink(path);
    }
}

TEST(root_reuse_exact_keys_errors_eviction_and_reset) {
    compute_root_reuse s={0};uint32_t key[7]={7,0,0,0,0,0x3f800000u,0};
    root_reuse_record(&s,key,COMPUTE_SPHERE,0);
    for(unsigned c=0;c<ROOT_REUSE_CAPACITIES;++c)ASSERT_EQ(s.counts[COMPUTE_SPHERE][ROOT_HITS+c],0);
    root_reuse_record(&s,key,COMPUTE_SPHERE,0);
    for(unsigned j=0;j<7;++j) {
        /* Includes +0/-0: equality is bitwise, not floating point. */
        key[j]^=0x80000000u;
        root_reuse_record(&s,key,COMPUTE_SPHERE,0);
        for(unsigned c=0;c<ROOT_REUSE_CAPACITIES;++c)ASSERT_EQ(s.counts[COMPUTE_SPHERE][ROOT_HITS+c],1+j);
        root_reuse_record(&s,key,COMPUTE_SPHERE,0);
    }
    root_reuse_record(&s,key,COMPUTE_SPHERE,1);
    ASSERT_EQ(s.counts[COMPUTE_SPHERE][ROOT_ERRORS],1);
    for(unsigned c=0;c<ROOT_REUSE_CAPACITIES;++c)ASSERT_EQ(s.counts[COMPUTE_SPHERE][ROOT_HITS+c],8);
    root_reuse_reset(&s);root_reuse_record(&s,key,COMPUTE_SPHERE,0);
    for(unsigned c=0;c<ROOT_REUSE_CAPACITIES;++c)ASSERT_EQ(s.counts[COMPUTE_SPHERE][ROOT_HITS+c],8);
    memset(&s,0,sizeof(s));
    for(unsigned i=0;i<513;++i){key[0]=i;root_reuse_record(&s,key,COMPUTE_SPHERE,0);}
    for(unsigned c=0;c<ROOT_REUSE_CAPACITIES;++c) {
        ASSERT_EQ(s.counts[COMPUTE_SPHERE][ROOT_HITS+c],0);
        ASSERT_TRUE(s.counts[COMPUTE_SPHERE][ROOT_EVICTIONS+c]>0);
    }
    root_reuse_record(&s,key,COMPUTE_SPHERE,0);
    for(unsigned c=0;c<ROOT_REUSE_CAPACITIES;++c)ASSERT_EQ(s.counts[COMPUTE_SPHERE][ROOT_HITS+c],1);
    memset(&s,0,sizeof(s));root_reuse_record(&s,key,COMPUTE_SPHERE,2);
    ASSERT_EQ(s.counts[COMPUTE_SPHERE][ROOT_NONFINITE],1);
    root_reuse_record(&s,key,COMPUTE_SPHERE,0);
    for(unsigned c=0;c<ROOT_REUSE_CAPACITIES;++c)ASSERT_EQ(s.counts[COMPUTE_SPHERE][ROOT_HITS+c],0);
    memset(&s,0,sizeof(s));root_reuse_record(&s,key,COMPUTE_SPHERE,1);
    root_reuse_record(&s,key,COMPUTE_SPHERE,0);
    for(unsigned c=0;c<ROOT_REUSE_CAPACITIES;++c)ASSERT_EQ(s.counts[COMPUTE_SPHERE][ROOT_HITS+c],0);
    s.counts[COMPUTE_SPHERE][ROOT_REQUESTS]=UINT32_MAX;
    root_reuse_record(&s,key,COMPUTE_SPHERE,0);
    ASSERT_EQ(s.counts[COMPUTE_SPHERE][ROOT_REQUESTS],UINT32_MAX);ASSERT_EQ(s.saturated,1);
}

static int reuse_progress(void* u,const render_framebuffer_t* fb,size_t done,size_t total) {
    (void)fb;(void)done;(void)total;++*(int*)u;return 0;
}
TEST(metal_root_reuse_preserves_output_samples_and_frames) {
    @autoreleasepool {
        char path[160];snprintf(path,sizeof(path),"/tmp/alea-root-reuse-%ld.csv",(long)getpid());
        const char* flag=getenv("ALEA_METAL_COMPACT_CODE");char* saved=flag?strdup(flag):NULL;
        alea_system_t* sys=alea_create();int mat=alea_add_material(sys,7);
        int box=alea_box_surface(sys,1,-5,5,-5,5,-5,5);
        int sphere=alea_sphere_surface(sys,2,0,0,0,4);
        alea_node_id_t shell=alea_intersection(sys,alea_halfspace(sys,box,-1),alea_halfspace(sys,sphere,1));
        ASSERT(alea_add_cell(sys,42,shell,mat,-1,1)>=0);
        for(int placed=0;placed<2;++placed) {
            int outer=alea_box_surface(sys,3+placed,20*placed-6,20*placed+6,-6,6,-6,6);
            int parent=alea_add_cell(sys,10+placed,alea_halfspace(sys,outer,-1),ALEA_MATERIAL_VOID,0,0);
            double transform[]={20*placed,0,0,0,1,0,-1,0,0,0,0,1};
            ASSERT_EQ(alea_add_transform(sys,1+placed,transform,placed?12:3,0),0);
            ASSERT_EQ(alea_set_fill(sys,parent,1,1+placed),0);
        }
        render_config_t cfg;render_config_init(&cfg);cfg.width=cfg.height=1;cfg.log_level=0;
        cfg.eye_set=cfg.target_set=1;cfg.target[1]=1;
        render_framebuffer_t* reference=render_framebuffer_create(1,1,1);
        render_framebuffer_t* measured=render_framebuffer_create(1,1,1);
        for(unsigned format=0;format<2;++format) {
            ASSERT_EQ(alea_setenv("ALEA_METAL_COMPACT_CODE",format?"1":"0",1),0);
            ASSERT_EQ(alea_unsetenv("ALEA_METAL_ROOT_REUSE"),0);
            metal_state* normal=alea_compute_metal_create(sys);ASSERT_NOT_NULL(normal);
            ASSERT_EQ(alea_setenv("ALEA_METAL_ROOT_REUSE",path,1),0);
            metal_state* diagnostic=alea_compute_metal_create(sys);ASSERT_NOT_NULL(diagnostic);
            ASSERT_EQ(alea_unsetenv("ALEA_METAL_ROOT_REUSE"),0);
            for(int placed=0;placed<2;++placed)for(int clip=0;clip<2;++clip) {
                cfg.eye[0]=cfg.target[0]=20*placed;cfg.eye[1]=cfg.eye[2]=cfg.target[2]=0;
                cfg.num_clips=clip;cfg.clips[0]=(render_clip_plane_t){{0,1,0},-1};
                render_camera_t cam;ASSERT_EQ(render_camera_setup(&cam,&cfg,sys),0);
                /* Identical AA rays make failure to reset keys directly visible. */
                cam.half_width=cam.half_height=0;
                compute_root_reuse single;
                for(unsigned aa=1;aa<=2;++aa)for(unsigned repeat=0;repeat<2;++repeat) {
                    cfg.aa_samples=aa;
                    ASSERT_EQ(alea_compute_metal_render(normal,&cfg,&cam,reference,NULL,NULL),0);
                    ASSERT_MSG(alea_compute_metal_render(diagnostic,&cfg,&cam,measured,NULL,NULL)==0,alea_error());
                    ASSERT_NULL(normal->root_reuse);
                    ASSERT_EQ(reference->cell_id[0],42);ASSERT_EQ(measured->cell_id[0],42);
                    ASSERT_EQ(measured->material_id[0],reference->material_id[0]);
                    ASSERT_EQ(memcmp(measured->depth,reference->depth,sizeof(float)),0);
                    ASSERT_EQ(memcmp(measured->normal,reference->normal,3*sizeof(float)),0);
                    ASSERT_EQ(memcmp(measured->color,reference->color,3*sizeof(float)),0);
                    const compute_root_reuse* row=diagnostic->root_reuse.contents;
                    ASSERT_EQ(row->saturated,0);
                    if(aa==1 && repeat==0)single=*row;
                    unsigned requests=0,hits=0;
                    for(unsigned f=0;f<ROOT_REUSE_FAMILIES;++f) {
                        requests+=row->counts[f][ROOT_REQUESTS];hits+=row->counts[f][ROOT_HITS+3];
                        ASSERT_EQ(row->counts[f][ROOT_ERRORS],0);
                        for(unsigned c=0;c<ROOT_REUSE_COUNTERS;++c)
                            ASSERT_EQ(row->counts[f][c],single.counts[f][c]*aa*aa);
                    }
                    ASSERT_TRUE(requests>0);ASSERT_TRUE(hits>0);
                    FILE* report=fopen(path,"r");ASSERT_NOT_NULL(report);char line[512];unsigned lines=0;
                    while(fgets(line,sizeof(line),report))++lines;
                    fclose(report);ASSERT_EQ(lines,1+ROOT_REUSE_FAMILIES*ROOT_REUSE_CAPACITIES);
                }
            }
            render_camera_t cam;ASSERT_EQ(render_camera_setup(&cam,&cfg,sys),0);int callbacks=0;
            ASSERT_EQ(alea_compute_metal_render(diagnostic,&cfg,&cam,measured,reuse_progress,&callbacks),-1);
            ASSERT_EQ(alea_error_code(),ALEA_ERR_UNSUPPORTED);ASSERT_EQ(callbacks,0);
            render_framebuffer_t large={.width=1025,.height=1};
            ASSERT_EQ(alea_compute_metal_render(diagnostic,&cfg,&cam,&large,NULL,NULL),-1);
            ASSERT_EQ(alea_error_code(),ALEA_ERR_OVERFLOW);
            alea_compute_metal_destroy(normal);alea_compute_metal_destroy(diagnostic);
        }
        render_framebuffer_free(reference);render_framebuffer_free(measured);render_config_free(&cfg);alea_destroy(sys);
        ASSERT_EQ(saved?alea_setenv("ALEA_METAL_COMPACT_CODE",saved,1):alea_unsetenv("ALEA_METAL_COMPACT_CODE"),0);free(saved);unlink(path);
    }
}
TEST(metal_root_reuse_allocation_failure_retries) {
    @autoreleasepool {
        metal_state* s=calloc(1,sizeof(*s));ASSERT_NOT_NULL(s);
        UploadDevice* device=[UploadDevice alloc];device->target=MTLCreateSystemDefaultDevice();
        device->real_buffers=YES;device->fail_at=1;s->device=(id<MTLDevice>)device;
        ASSERT_EQ(prepare_root_reuse(s,1),-1);ASSERT_EQ(alea_error_code(),ALEA_ERR_OUT_OF_MEMORY);ASSERT_NULL(s->root_reuse);
        device->fail_at=0;ASSERT_EQ(prepare_root_reuse(s,1),0);
        id<MTLBuffer> original=s->root_reuse;
        device->fail_at=device->calls+1;ASSERT_EQ(prepare_root_reuse(s,2),-1);ASSERT_TRUE(s->root_reuse==original);
        device->fail_at=0;ASSERT_EQ(prepare_root_reuse(s,2),0);
        ASSERT_EQ(prepare_root_reuse(s,ROOT_REUSE_MAX_PIXELS+1),-1);ASSERT_EQ(alea_error_code(),ALEA_ERR_OVERFLOW);
        alea_compute_metal_destroy(s);
    }
}


TEST(metal_root_reuse_preserves_solver_errors_and_report_failures) {
    @autoreleasepool {
        char path[160];snprintf(path,sizeof(path),"/tmp/alea-root-errors-%ld.csv",(long)getpid());unlink(path);
        alea_system_t* sys=alea_create();int mat=alea_add_material(sys,7);
        int sphere=alea_sphere_surface(sys,1,0,0,0,1e20);
        ASSERT(alea_add_cell(sys,42,alea_halfspace(sys,sphere,-1),mat,-1,0)>=0);
        ASSERT_EQ(alea_unsetenv("ALEA_METAL_ROOT_REUSE"),0);
        metal_state* normal=alea_compute_metal_create(sys);ASSERT_NOT_NULL(normal);
        ASSERT_EQ(alea_setenv("ALEA_METAL_ROOT_REUSE",path,1),0);
        metal_state* diagnostic=alea_compute_metal_create(sys);ASSERT_NOT_NULL(diagnostic);
        ASSERT_EQ(alea_unsetenv("ALEA_METAL_ROOT_REUSE"),0);
        render_config_t cfg;render_config_init(&cfg);cfg.width=cfg.height=1;
        cfg.eye_set=cfg.target_set=1;cfg.eye[0]=cfg.eye[2]=0;cfg.eye[1]=-10;
        cfg.target[0]=cfg.target[1]=cfg.target[2]=0;
        render_camera_t cam;ASSERT_EQ(render_camera_setup(&cam,&cfg,sys),0);
        render_framebuffer_t* fb=render_framebuffer_create(1,1,1);ASSERT_NOT_NULL(fb);
        ASSERT_EQ(alea_compute_metal_render(normal,&cfg,&cam,fb,NULL,NULL),-1);ASSERT_EQ(alea_error_code(),ALEA_ERR_OVERFLOW);
        ASSERT_EQ(alea_compute_metal_render(diagnostic,&cfg,&cam,fb,NULL,NULL),-1);ASSERT_EQ(alea_error_code(),ALEA_ERR_OVERFLOW);
        const compute_root_reuse* row=diagnostic->root_reuse.contents;
        ASSERT_EQ(row->counts[COMPUTE_SPHERE][ROOT_REQUESTS],1);
        ASSERT_EQ(row->counts[COMPUTE_SPHERE][ROOT_ERRORS],1);
        for(unsigned c=0;c<ROOT_REUSE_CAPACITIES;++c)ASSERT_EQ(row->counts[COMPUTE_SPHERE][ROOT_HITS+c],0);
        ASSERT_NE(access(path,F_OK),0); /* A failed frame must not look complete. */
        cfg.num_clips=1;cfg.clips[0]=(render_clip_plane_t){{0,-1,0},-20};
        free(diagnostic->root_reuse_path);diagnostic->root_reuse_path=strdup("/dev/null/unwritable");
        ASSERT_EQ(alea_compute_metal_render(diagnostic,&cfg,&cam,fb,NULL,NULL),-1);ASSERT_EQ(alea_error_code(),ALEA_ERR_INVALID_STATE);
        free(diagnostic->root_reuse_path);diagnostic->root_reuse_path=strdup(path);
        ASSERT_EQ(alea_compute_metal_render(diagnostic,&cfg,&cam,fb,NULL,NULL),0);
        ASSERT_EQ(fb->cell_id[0],-1);ASSERT_EQ(row->counts[COMPUTE_SPHERE][ROOT_REQUESTS],0);
        render_framebuffer_free(fb);render_config_free(&cfg);alea_destroy(sys);
        alea_compute_metal_destroy(normal);alea_compute_metal_destroy(diagnostic);unlink(path);
    }
}

#include "test_group_trace.h"
TEST(package_groups_trace_values_errors_and_aliases) {
    @autoreleasepool {
        metal_state state={0};state.device=MTLCreateSystemDefaultDevice();state.queue=[state.device newCommandQueue];
        NSError* error=nil;dispatch_data_t data=dispatch_data_create(alea_group_eval_probe_metallib,alea_group_eval_probe_metallib_len,NULL,DISPATCH_DATA_DESTRUCTOR_DEFAULT);
        id<MTLLibrary> library=[state.device newLibraryWithData:data error:&error];dispatch_release(data);ASSERT_NOT_NULL(library);
        for(unsigned format=0;format<2;++format){state.compact_code=format;ASSERT_EQ(package_trace_test(&state,library),0);}
        [library release];[state.queue release];[state.device release];
    }
}
TEST(package_groups_allocation_fallback_and_retry) {
    @autoreleasepool {
        compute_node code[]={{COMPUTE_LEAF,0},{COMPUTE_SKIP_INTERSECTION,2},{COMPUTE_LEAF,2},{COMPUTE_INTERSECTION,0}};
        compute_cell cell={.root=0,.count=4};compute_scene scene={.nodes=code,.node_count=4,.cells=&cell,.cell_count=1,.primitive_count=2};
        for(NSUInteger failure=1;failure<=2;++failure) {
            metal_state* s=calloc(1,sizeof(*s));ASSERT_NOT_NULL(s);
            UploadDevice* device=[UploadDevice alloc];device->target=MTLCreateSystemDefaultDevice();device->real_buffers=YES;device->fail_at=failure;s->device=(id<MTLDevice>)device;
            ASSERT_EQ(prepare_groups(s,&scene),0);ASSERT_FALSE(s->groups_enabled);ASSERT_NULL(s->group_maps);ASSERT_NULL(s->group_refs);
            device->fail_at=0;ASSERT_EQ(prepare_groups(s,&scene),1);ASSERT_TRUE(s->groups_enabled);
            ASSERT_NOT_NULL(s->group_maps);ASSERT_NOT_NULL(s->group_refs);alea_compute_metal_destroy(s);
        }
        metal_state* s=calloc(1,sizeof(*s));s->device=MTLCreateSystemDefaultDevice();
        alea_interrupt();int rc=prepare_groups(s,&scene);alea_clear_interrupt();ASSERT_EQ(rc,-1);ASSERT_EQ(alea_error_code(),ALEA_ERR_INTERRUPTED);
        ASSERT_NULL(s->group_maps);ASSERT_EQ(prepare_groups(s,&scene),1);alea_compute_metal_destroy(s);
    }
}
TEST(package_groups_configuration_and_profile_units) {
    @autoreleasepool {
        const char* keys[]={"ALEA_METAL_CSG_GROUPS","ALEA_METAL_POINT_PROFILE","ALEA_METAL_SELECT_PROFILE","ALEA_METAL_POINT_CANDIDATES"};char* saved[4]={0};
        for(unsigned k=0;k<4;++k){const char* v=getenv(keys[k]);saved[k]=v?strdup(v):NULL;}
        char point_path[160],select_path[160];snprintf(point_path,sizeof(point_path),"/tmp/alea-group-point-%ld.jsonl",(long)getpid());snprintf(select_path,sizeof(select_path),"/tmp/alea-group-select-%ld.jsonl",(long)getpid());unlink(point_path);unlink(select_path);
        alea_system_t* sys=alea_create();int mat=alea_add_material(sys,7);
        int sphere=alea_sphere_surface(sys,1,0,0,0,2),box=alea_box_surface(sys,2,-4,4,-4,4,-4,4);
        ASSERT(alea_add_cell(sys,42,alea_intersection(sys,alea_halfspace(sys,sphere,-1),alea_halfspace(sys,box,-1)),mat,-1,0)>=0);
        ASSERT_EQ(alea_setenv(keys[0],"invalid",1),0);ASSERT_NULL(alea_compute_metal_create(sys));ASSERT_EQ(alea_error_code(),ALEA_ERR_INVALID_ARG);
        ASSERT_EQ(alea_setenv(keys[0],"1",1),0);ASSERT_EQ(alea_setenv(keys[1],point_path,1),0);ASSERT_EQ(alea_setenv(keys[2],select_path,1),0);ASSERT_EQ(alea_setenv(keys[3],"0",1),0);
        metal_state* s=alea_compute_metal_create(sys);ASSERT_NOT_NULL(s);ASSERT_TRUE(s->groups_enabled);alea_destroy(sys);
        for(unsigned k=0;k<4;++k){ASSERT_EQ(saved[k]?alea_setenv(keys[k],saved[k],1):alea_unsetenv(keys[k]),0);free(saved[k]);}
        double xyz[]={0,0,0,3,0,0};alea_compute_result_t result[2];alea_compute_selection_t selected[2];
        ASSERT_EQ(alea_compute_metal_classify(s,xyz,2,result),0);ASSERT_EQ(result[0].cell_id,42);ASSERT_EQ(result[1].owner_count,0);
        ASSERT_EQ(alea_compute_metal_select(s,xyz,2,selected),0);ASSERT_EQ(selected[0].cell_id,42);ASSERT_EQ(selected[1].status,ALEA_COMPUTE_NOT_FOUND);
        char line[8192];const char* paths[]={point_path,select_path};
        for(unsigned k=0;k<2;++k){FILE* f=fopen(paths[k],"r");ASSERT_NOT_NULL(f);ASSERT_NOT_NULL(fgets(line,sizeof(line),f));fclose(f);ASSERT_NOT_NULL(strstr(line,"\"csg_groups\":true"));ASSERT_NOT_NULL(strstr(line,"\"instruction_unit\":\"bytecode_and_group_ops\""));}
        ASSERT_NOT_NULL(strstr(line,"\"primitives\":3"));
        alea_compute_metal_destroy(s);unlink(point_path);unlink(select_path);
    }
}

#include "ray_eval_probe_metallib.h"
#define ALEA_RAY_CONTRACT_GROUPS 1
#include "test_ray_contract.h"
TEST(package_ray_groups_analytic_contract) {
    @autoreleasepool {
        metal_state state={0};state.device=MTLCreateSystemDefaultDevice();state.queue=[state.device newCommandQueue];
        NSError* error=nil;dispatch_data_t data=dispatch_data_create(alea_ray_eval_probe_metallib,alea_ray_eval_probe_metallib_len,NULL,DISPATCH_DATA_DESTRUCTOR_DEFAULT);
        id<MTLLibrary> library=[state.device newLibraryWithData:data error:&error];dispatch_release(data);ASSERT_NOT_NULL(library);
        ASSERT_EQ(ray_contract_check(&state,library),0);
        [library release];[state.queue release];[state.device release];
    }
}
TEST(package_ray_groups_independent_settings_and_fallback) {
    const char* keys[]={"ALEA_METAL_CSG_GROUPS","ALEA_METAL_RAY_GROUPS"};char* saved[2]={0};
    for(unsigned k=0;k<2;++k){const char* v=getenv(keys[k]);saved[k]=v?strdup(v):NULL;}
    alea_system_t* sys=alea_create();int mat=alea_add_material(sys,7),sphere=alea_sphere_surface(sys,1,0,0,0,1);
    alea_node_id_t leaf=alea_halfspace(sys,sphere,-1);alea_add_cell(sys,42,alea_intersection(sys,leaf,leaf),mat,-1,0);
    ASSERT_EQ(alea_setenv(keys[1],"invalid",1),0);ASSERT_NULL(alea_compute_metal_create(sys));ASSERT_EQ(alea_error_code(),ALEA_ERR_INVALID_ARG);
    metal_state* contexts[4];
    for(unsigned mode=0;mode<4;++mode) {
        ASSERT_EQ(alea_setenv(keys[0],mode&1?"1":"0",1),0);ASSERT_EQ(alea_setenv(keys[1],mode&2?"1":"0",1),0);
        contexts[mode]=alea_compute_metal_create(sys);ASSERT_NOT_NULL(contexts[mode]);
        ASSERT_EQ(contexts[mode]->groups_enabled,!!(mode&1));ASSERT_EQ(contexts[mode]->ray_groups_enabled,!!(mode&2));
        ASSERT_EQ(contexts[mode]->group_maps!=nil,mode!=0);ASSERT_EQ(contexts[mode]->group_refs!=nil,mode!=0);
    }
    alea_destroy(sys);
    ASSERT_EQ(alea_setenv(keys[0],"0",1),0);ASSERT_EQ(alea_setenv(keys[1],"0",1),0);
    for(unsigned mode=0;mode<4;++mode){ASSERT_EQ(contexts[mode]->groups_enabled,!!(mode&1));ASSERT_EQ(contexts[mode]->ray_groups_enabled,!!(mode&2));alea_compute_metal_destroy(contexts[mode]);}
    sys=alea_create();sphere=alea_sphere_surface(sys,1,0,0,0,1);mat=alea_add_material(sys,7);alea_add_cell(sys,42,alea_halfspace(sys,sphere,-1),mat,-1,0);
    ASSERT_EQ(alea_setenv(keys[1],"1",1),0);metal_state* fallback=alea_compute_metal_create(sys);ASSERT_NOT_NULL(fallback);
    ASSERT_FALSE(fallback->ray_groups_enabled);ASSERT_NULL(fallback->group_maps);alea_compute_metal_destroy(fallback);alea_destroy(sys);
    for(unsigned k=0;k<2;++k){ASSERT_EQ(saved[k]?alea_setenv(keys[k],saved[k],1):alea_unsetenv(keys[k]),0);free(saved[k]);}
}

TEST(package_ray_groups_profile_units_and_generic_progressive) {
    const char* keys[]={"ALEA_METAL_CSG_GROUPS","ALEA_METAL_RAY_GROUPS","ALEA_METAL_PROFILE"};char* saved[3]={0};
    for(unsigned k=0;k<3;++k){const char* v=getenv(keys[k]);saved[k]=v?strdup(v):NULL;}
    char prefix[160],path[200];snprintf(prefix,sizeof(prefix),"/tmp/alea-ray-group-profile-%ld",(long)getpid());
    ASSERT_EQ(alea_setenv(keys[0],"0",1),0);ASSERT_EQ(alea_setenv(keys[1],"1",1),0);ASSERT_EQ(alea_setenv(keys[2],prefix,1),0);
    alea_system_t* sys=alea_create();int mat=alea_add_material(sys,7),sphere=alea_sphere_surface(sys,1,0,0,0,1);
    alea_node_id_t leaf=alea_halfspace(sys,sphere,-1);alea_add_cell(sys,42,alea_intersection(sys,leaf,leaf),mat,-1,0);
    render_config_t cfg;render_config_init(&cfg);cfg.width=cfg.height=1;cfg.eye_set=cfg.target_set=1;cfg.eye[0]=cfg.eye[2]=0;cfg.eye[1]=-5;cfg.target[0]=cfg.target[1]=cfg.target[2]=0;
    render_camera_t camera;ASSERT_EQ(render_camera_setup(&camera,&cfg,sys),0);metal_state* s=alea_compute_metal_create(sys);ASSERT_NOT_NULL(s);ASSERT_TRUE(s->ray_groups_enabled);ASSERT_FALSE(s->groups_enabled);alea_destroy(sys);
    for(unsigned k=0;k<3;++k){ASSERT_EQ(saved[k]?alea_setenv(keys[k],saved[k],1):alea_unsetenv(keys[k]),0);free(saved[k]);}
    render_framebuffer_t* fb=render_framebuffer_create(1,1,1);ASSERT_NOT_NULL(fb);
    for(unsigned progressive=0;progressive<2;++progressive) {
        int calls=0;ASSERT_EQ(alea_compute_metal_render(s,&cfg,&camera,fb,progressive?reuse_progress:NULL,&calls),0);
        ASSERT_EQ(fb->cell_id[0],42);snprintf(path,sizeof(path),"%s-rays.csv",prefix);FILE* f=fopen(path,"r");ASSERT_NOT_NULL(f);char line[1024];ASSERT_NOT_NULL(fgets(line,sizeof(line),f));ASSERT_NOT_NULL(strstr(line,"csg_groups,instruction_unit"));ASSERT_NOT_NULL(fgets(line,sizeof(line),f));fclose(f);
        ASSERT_NOT_NULL(strstr(line,progressive?",0,bytecode":",1,bytecode_and_group_ops"));
    }
    const char* suffixes[]={"-rays.csv","-cells.csv","-heatmap.png"};for(unsigned k=0;k<3;++k){snprintf(path,sizeof(path),"%s%s",prefix,suffixes[k]);ASSERT_EQ(unlink(path),0);}
    render_framebuffer_free(fb);render_config_free(&cfg);alea_compute_metal_destroy(s);
}
TEST_MAIN()
