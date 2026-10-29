// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
// SPDX-License-Identifier: MPL-2.0
#include "alea_test.h"
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
 * these tests never submit them to a command encoder. */
@interface UploadDevice : NSProxy {
@public
    id<MTLDevice> target;
    NSUInteger calls,fail_at;
}
@end
@implementation UploadDevice
- (NSMethodSignature*)methodSignatureForSelector:(SEL)selector {
    return [(NSObject*)target methodSignatureForSelector:selector];
}
- (void)forwardInvocation:(NSInvocation*)invocation {[invocation invokeWithTarget:target];}
- (id<MTLBuffer>)newBufferWithLength:(NSUInteger)length options:(MTLResourceOptions)options {
    (void)options;
    if(++calls==fail_at)return nil;
    UploadBuffer* buffer=[[UploadBuffer alloc] init];buffer.length=length;
    return (id<MTLBuffer>)buffer;
}
- (id<MTLBuffer>)newBufferWithBytes:(const void*)bytes length:(NSUInteger)length options:(MTLResourceOptions)options {
    (void)bytes;return [self newBufferWithLength:length options:options];
}
- (void)dealloc {[target release];[super dealloc];}
@end

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
            if(failure<=4) {
                device->fail_at=failure;
                ASSERT_EQ(prepare_classify(s),-1);ASSERT_EQ(live_buffers,0);
                ASSERT_NULL(s->point_cells);ASSERT_NULL(s->transforms);
                ASSERT_TRUE(s->scene.point_cells==original);
                device->fail_at=0;
            }
            ASSERT_EQ(prepare_classify(s),0);ASSERT_EQ(live_buffers,4);
            ASSERT_NULL(s->scene.point_cells);
            id<MTLBuffer> transforms=s->transforms;
            NSUInteger calls=device->calls;
            ASSERT_EQ(prepare_classify(s),0);ASSERT_EQ(device->calls,calls);
            if(failure>4) {
                device->fail_at=device->calls+failure-4;
                ASSERT_EQ(prepare_render(s),-1);ASSERT_EQ(live_buffers,4);
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

TEST_MAIN()
