// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
// SPDX-License-Identifier: MPL-2.0
/* Included by the Objective-C bridge after its context definition. */
static int prepare_groups(metal_state* s,const compute_scene* scene) {
    const size_t budget=64u*1024u*1024u;compute_group_plan plan;
    int rc=alea_compute_groups_plan(scene,budget,&plan);if(rc<=0)return rc;
    id<MTLBuffer> maps=nil,refs=nil;size_t leaves=0,hash_capacity=1;
    while(hash_capacity<2*plan.program_count)hash_capacity*=2;
    size_t bytes=16+plan.program_count*16+plan.span_count*12+hash_capacity*4;
    size_t scratch=plan.program_capacity*sizeof(*plan.programs)+plan.span_capacity*sizeof(*plan.spans);
    if(bytes>budget || scratch>budget-bytes || bytes>s->device.maxBufferLength)goto fallback;
    for(size_t i=0;i<plan.span_count;++i) {
        uint32_t length=plan.spans[i].end-plan.spans[i].begin;
        if(length<4 || length%3!=1)goto invalid;
        leaves+=(length+2)/3;
    }
    size_t ref_bytes=MAX((size_t)4,leaves*4);
    if(ref_bytes>budget-bytes || ref_bytes>s->device.maxBufferLength)goto fallback;
    maps=[s->device newBufferWithLength:bytes options:MTLResourceStorageModeShared];if(!maps)goto fallback;
    uint32_t* target=maps.contents;
    target[0]=(uint32_t)plan.program_count;target[1]=(uint32_t)plan.span_count;
    target[2]=(uint32_t)(hash_capacity-1);target[3]=0;
    _Static_assert(sizeof(compute_group_program)==16,"Metal group program layout");
    memcpy(target+4,plan.programs,plan.program_count*16);
    size_t base=4+4*plan.program_count,offset=0;
    for(size_t i=0;i<plan.span_count;++i) {
        if(!(i&4095) && alea_interrupted())goto interrupted;
        target[base+3*i]=plan.spans[i].begin;target[base+3*i+1]=plan.spans[i].end;
        target[base+3*i+2]=(uint32_t)offset;offset+=(plan.spans[i].end-plan.spans[i].begin+2)/3;
    }
    uint32_t* slots=target+base+3*plan.span_count;memset(slots,0,hash_capacity*4);
    for(size_t i=0;i<plan.program_count;++i) {
        if(!(i&4095) && alea_interrupted())goto interrupted;
        uint32_t h=(target[4+4*i]*0x9e3779b9u)^(target[5+4*i]*0x85ebca6bu);h=(h^(h>>16))&(uint32_t)(hash_capacity-1);
        while(slots[h])h=(h+1)&(uint32_t)(hash_capacity-1);slots[h]=(uint32_t)i+1;
    }
    size_t groups=plan.span_count;
    alea_compute_groups_free(&plan);
    refs=[s->device newBufferWithLength:ref_bytes options:MTLResourceStorageModeShared];if(!refs)goto fallback;
    uint32_t* pool=refs.contents;size_t written=0;
    for(size_t i=0;i<groups;++i) {
        if(!(i&4095) && alea_interrupted())goto interrupted;
        uint32_t begin=target[base+3*i],end=target[base+3*i+1];
        if(written!=target[base+3*i+2])goto invalid;
        for(uint32_t pc=begin;pc<end;++pc) {
            compute_node n=alea_compute_node_at(scene,pc);
            if(n.op==COMPUTE_LEAF){if(written>=leaves)goto invalid;pool[written++]=n.arg;}
        }
    }
    if(written!=leaves)goto invalid;
    if(alea_interrupted())goto interrupted;
    s->group_maps=maps;s->group_refs=refs;s->groups_enabled=true;
    fprintf(stderr,"Metal CSG groups: groups=%zu bytes=%zu\n",groups,bytes+ref_bytes);return 1;
 interrupted:alea_set_error_detail(ALEA_ERR_INTERRUPTED,"CSG group preparation interrupted");goto fail;
 invalid:alea_set_error_detail(ALEA_ERR_INVALID_STATE,"Invalid CSG group list");
 fail:alea_compute_groups_free(&plan);[maps release];[refs release];return -1;
 fallback:alea_compute_groups_free(&plan);[maps release];[refs release];
    if(alea_interrupted()){alea_set_error_detail(ALEA_ERR_INTERRUPTED,"CSG group preparation interrupted");return -1;}
    return 0;
}
