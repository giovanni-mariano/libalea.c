// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
// SPDX-License-Identifier: MPL-2.0
#include "region_group_audit.h"
typedef struct {uint32_t* data;size_t count,capacity,limit;} group_storage;
static int group_append(void* user,region_group group) {
    group_storage* s=user;
    if(s->count==s->capacity) {
        size_t next=s->capacity?s->capacity*2:1024;
        if(next>(s->limit-4)/8)next=(s->limit-4)/8;
        if(next<=s->capacity){fprintf(stderr,"Group ranges exceed bounded allowance\n");return -1;}
        uint32_t* data=realloc(s->data,4+next*8);if(!data)return -1;
        s->data=data;s->capacity=next;
    }
    s->data[1+2*s->count]=group.begin;s->data[2+2*s->count]=group.end;++s->count;return 0;
}
static int group_compare(const void* x,const void* y) {
    const region_group* a=x,*b=y;
    if(a->begin!=b->begin)return a->begin<b->begin?-1:1;
    return a->end>b->end?-1:a->end!=b->end;
}
static id<MTLBuffer> group_buffer(metal_state* s,const uint32_t* metadata) {
    compute_scene host=s->scene;group_storage storage={0};
    if(s->compact_code)host.compact_nodes=s->nodes.contents;else if(s->nodes)host.nodes=s->nodes.contents;
    size_t records_bytes=(size_t)metadata[0]*16;
    size_t counter_bytes=4+records_bytes+256*8*AUDIT_REASONS*AUDIT_FIELDS*4;
    if(records_bytes+counter_bytes+16>64*1024*1024)return nil;
    storage.limit=(64*1024*1024-records_bytes-counter_bytes)/2;
    if(storage.limit>24*1024*1024)storage.limit=24*1024*1024;
    uint32_t* records=malloc(records_bytes?records_bytes:4);if(!records)return nil;
    for(size_t i=0;i<metadata[0];++i) {
        compute_cell c={.root=metadata[1+4*i],.count=metadata[2+4*i]};size_t first=storage.count;
        if(region_group_program(&host,c,group_append,&storage)) {
            fprintf(stderr,"Group recognition failed\n");free(records);free(storage.data);return nil;
        }
        if(storage.count>first)qsort(storage.data+1+2*first,storage.count-first,8,group_compare);
        for(size_t j=first+1;j<storage.count;++j)if(storage.data[1+2*j]<storage.data[2+2*(j-1)]) {
            fprintf(stderr,"Overlapping groups within a program\n");free(records);free(storage.data);return nil;
        }
        records[4*i]=c.root;records[4*i+1]=c.count;
        records[4*i+2]=(uint32_t)first;records[4*i+3]=(uint32_t)(storage.count-first);
    }
    size_t bytes=4+records_bytes+storage.count*8;
    if(bytes+records_bytes+4+storage.capacity*8+counter_bytes>64*1024*1024) {
        fprintf(stderr,"Group audit exceeds combined 64 MiB allowance\n");free(records);free(storage.data);return nil;
    }
    uint64_t covered=0;
    for(size_t i=0;i<storage.count;++i)covered+=storage.data[2+2*i]-storage.data[1+2*i];
    printf("groups=%zu group_static_instructions=%llu group_bytes=%zu group_capacity_bytes=%zu\n",storage.count,(unsigned long long)covered,bytes,4+8*storage.capacity);
    id<MTLBuffer> buffer=[s->device newBufferWithLength:bytes options:MTLResourceStorageModeShared];
    if(buffer) {
        uint32_t* target=buffer.contents;target[0]=metadata[0];memcpy(target+1,records,records_bytes);
        if(storage.count)memcpy(target+1+4*(size_t)metadata[0],storage.data+1,storage.count*8);
    }
    free(records);free(storage.data);return buffer;
}
static int group_map_self_test(metal_state* s,id<MTLLibrary> library) {
    const uint32_t map[]={3,0,1,0,0,0,7,0,1,2,4,1,1,0,7,2,6};
    const uint32_t expected[30]={0,0,0,0,0, 1,1,1,0,0, 1,1,1,0,0,
                                1,0,0,0,1, 1,0,0,0,1, 0,0,0,0,0};
    NSError* error=nil;
    id<MTLFunction> function=[library newFunctionWithName:@"alea_group_map_probe"];
    id<MTLComputePipelineState> pipeline=function?[s->device newComputePipelineStateWithFunction:function error:&error]:nil;
    id<MTLBuffer> input=[s->device newBufferWithBytes:map length:sizeof(map) options:MTLResourceStorageModeShared];
    id<MTLBuffer> output=[s->device newBufferWithLength:sizeof(expected) options:MTLResourceStorageModeShared];
    int rc=-1;
    if(pipeline && input && output) {
        id<MTLCommandBuffer> command=[s->queue commandBuffer];id<MTLComputeCommandEncoder> e=[command computeCommandEncoder];
        if(command && e) {
            [e setComputePipelineState:pipeline];[e setBuffer:input offset:0 atIndex:0];[e setBuffer:output offset:0 atIndex:1];
            [e dispatchThreads:MTLSizeMake(6,1,1) threadsPerThreadgroup:MTLSizeMake(6,1,1)];[e endEncoding];[command commit];[command waitUntilCompleted];
            if(command.status==MTLCommandBufferStatusCompleted && !memcmp(output.contents,expected,sizeof(expected)))rc=0;
        }
    }
    [function release];[pipeline release];[input release];[output release];
    if(!rc)puts("group_map_alias_oracle=passed");else fprintf(stderr,"Group map alias oracle failed\n");return rc;
}
