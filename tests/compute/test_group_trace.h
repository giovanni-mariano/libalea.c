// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
// SPDX-License-Identifier: MPL-2.0
#include "group_eval_probe_metallib.h"
typedef struct {compute_group_span spans[8];size_t count;} trace_spans;
static int trace_emit(void* user,compute_group_span s) {trace_spans* p=user;if(p->count==8)return -1;p->spans[p->count++]=s;return 0;}
static int trace_compare(const void* a,const void* b) {const compute_group_span* x=a,*y=b;return x->begin<y->begin?-1:x->begin!=y->begin;}
static int trace_program_compare(const void* a,const void* b) {const uint32_t* x=a,*y=b;if(x[0]!=y[0])return x[0]<y[0]?-1:1;return x[1]<y[1]?-1:x[1]!=y[1];}
static id<MTLComputePipelineState> package_trace_pipeline(metal_state* s,id<MTLLibrary> library,const char* name,bool enabled,bool trace,bool profile) {
    (void)trace;NSError* error=nil;bool no=false;
    MTLFunctionConstantValues* c=[[MTLFunctionConstantValues alloc] init];
    [c setConstantValue:&profile type:MTLDataTypeBool atIndex:5];[c setConstantValue:&no type:MTLDataTypeBool atIndex:6];
    [c setConstantValue:&s->compact_code type:MTLDataTypeBool atIndex:7];[c setConstantValue:&enabled type:MTLDataTypeBool atIndex:9];
    id<MTLFunction> f=[library newFunctionWithName:[NSString stringWithUTF8String:name] constantValues:c error:&error];
    id<MTLComputePipelineState> pipeline=f?[s->device newComputePipelineStateWithFunction:f error:&error]:nil;
    [f release];[c release];return pipeline;
}
static int package_trace_test(metal_state* s,id<MTLLibrary> library) {
    const compute_node patterns[6][10]={
        {{COMPUTE_LEAF,0},{COMPUTE_SKIP_INTERSECTION,5},{COMPUTE_LEAF,2},{COMPUTE_SKIP_INTERSECTION,2},{COMPUTE_LEAF,4},{COMPUTE_INTERSECTION,0},{COMPUTE_INTERSECTION,0}},
        {{COMPUTE_LEAF,0},{COMPUTE_SKIP_INTERSECTION,5},{COMPUTE_LEAF,2},{COMPUTE_SKIP_INTERSECTION,2},{COMPUTE_LEAF,4},{COMPUTE_INTERSECTION,0},{COMPUTE_INTERSECTION,0},{COMPUTE_COMPLEMENT,0}},
        {{COMPUTE_LEAF,0},{COMPUTE_SKIP_INTERSECTION,2},{COMPUTE_LEAF,2},{COMPUTE_INTERSECTION,0},{COMPUTE_SKIP_UNION,2},{COMPUTE_LEAF,4},{COMPUTE_UNION,0}},
        {{COMPUTE_LEAF,4},{COMPUTE_SKIP_INTERSECTION,5},{COMPUTE_LEAF,0},{COMPUTE_SKIP_INTERSECTION,2},{COMPUTE_LEAF,2},{COMPUTE_INTERSECTION,0},{COMPUTE_DIFFERENCE,0}},
        {{COMPUTE_LEAF,0},{COMPUTE_SKIP_INTERSECTION,2},{COMPUTE_LEAF,2},{COMPUTE_INTERSECTION,0},{COMPUTE_SKIP_REVERSE_DIFFERENCE,2},{COMPUTE_LEAF,4},{COMPUTE_REVERSE_DIFFERENCE,0}},
        {{COMPUTE_LEAF,0},{COMPUTE_SKIP_INTERSECTION,2},{COMPUTE_LEAF,2},{COMPUTE_INTERSECTION,0},{COMPUTE_SKIP_UNION,5},{COMPUTE_LEAF,1},{COMPUTE_SKIP_INTERSECTION,2},{COMPUTE_LEAF,4},{COMPUTE_INTERSECTION,0},{COMPUTE_UNION,0}}};
    const uint32_t lengths[]={7,8,7,7,7,10};compute_node code[64];compute_cell cells[8]={0};
    uint32_t map[128]={8},refs[64];size_t used=0,nrefs=0,ngroups=0;
    trace_spans storage={0};compute_scene scene={.nodes=code,.primitive_count=3};
    for(unsigned k=0;k<6;++k){cells[k].root=(uint32_t)used;cells[k].count=lengths[k];memcpy(code+used,patterns[k],lengths[k]*sizeof(*code));used+=lengths[k];}
    cells[6]=(compute_cell){.root=0,.count=1};cells[7]=(compute_cell){.root=2,.count=4};
    scene.node_count=used;
    for(unsigned k=0;k<8;++k) {
        storage.count=0;if(alea_compute_group_program(&scene,cells[k],trace_emit,&storage)){return -1;}
        if(storage.count)qsort(storage.spans,storage.count,sizeof(compute_group_span),trace_compare);
        map[4+4*k]=cells[k].root;map[5+4*k]=cells[k].count;map[6+4*k]=(uint32_t)ngroups;map[7+4*k]=(uint32_t)storage.count;
        for(size_t j=0;j<storage.count;++j) {
            uint32_t begin=storage.spans[j].begin,end=storage.spans[j].end;
            map[36+3*ngroups]=begin;map[37+3*ngroups]=end;map[38+3*ngroups]=(uint32_t)nrefs;++ngroups;
            for(uint32_t pc=begin;pc<end;++pc)if(code[pc].op==COMPUTE_LEAF)refs[nrefs++]=code[pc].arg;
        }
    }
    map[1]=(uint32_t)ngroups;
    qsort(map+4,8,16,trace_program_compare);
    compute_primitive primitives[216*3]={0};float values[]={-1,-0.0f,0.0f,1,NAN,INFINITY};
    for(unsigned i=0;i<216;++i)for(unsigned j=0,divisor=1;j<3;++j,divisor*=6) {
        primitives[3*i+j].type=COMPUTE_PLANE;primitives[3*i+j].p[3]=values[i/divisor%6];
    }
    uint32_t compact[64];const void* encoded=code;size_t code_bytes=used*sizeof(*code);
    if(s->compact_code){if(alea_compute_code_compact(code,used,compact)!=1)return -1;encoded=compact;code_bytes=used*4;}
    id<MTLBuffer> inputs[5]={nil,nil,nil,nil,nil},outputs[3]={nil,nil,nil};id<MTLComputePipelineState> pipelines[3]={nil,nil,nil};
    const void* pointers[]={encoded,primitives,cells,map,refs};size_t sizes[]={code_bytes,sizeof(primitives),sizeof(cells),(36+3*ngroups+16)*4,nrefs*4};int rc=-1;
    for(unsigned i=0;i<5;++i){inputs[i]=[s->device newBufferWithBytes:pointers[i] length:sizes[i] options:MTLResourceStorageModeShared];if(!inputs[i])goto done;}
    for(unsigned i=0;i<3;++i) {
        if(i==2) {
            uint32_t* m=inputs[3].contents;m[2]=15;uint32_t* slots=m+36+3*ngroups;
            for(unsigned j=0;j<8;++j){uint32_t h=(m[4+4*j]*0x9e3779b9u)^(m[5+4*j]*0x85ebca6bu);h=(h^(h>>16))&15;while(slots[h])h=(h+1)&15;slots[h]=j+1;}
        }
        pipelines[i]=package_trace_pipeline(s,library,"alea_group_execute_probe",i!=0,true,false);
        outputs[i]=[s->device newBufferWithLength:1728*40 options:MTLResourceStorageModeShared];if(!pipelines[i] || !outputs[i])goto done;
        id<MTLCommandBuffer> command=[s->queue commandBuffer];id<MTLComputeCommandEncoder> e=[command computeCommandEncoder];if(!command || !e)goto done;
        [e setComputePipelineState:pipelines[i]];for(unsigned j=0;j<5;++j)[e setBuffer:inputs[j] offset:0 atIndex:j];[e setBuffer:outputs[i] offset:0 atIndex:5];
        [e dispatchThreads:MTLSizeMake(1728,1,1) threadsPerThreadgroup:MTLSizeMake(64,1,1)];[e endEncoding];[command commit];[command waitUntilCompleted];if(command.status!=MTLCommandBufferStatusCompleted)goto done;
    }
    if(memcmp(outputs[0].contents,outputs[1].contents,1728*40) || memcmp(outputs[0].contents,outputs[2].contents,1728*40)){fprintf(stderr,"Group scalar/trace/error mismatch\n");goto done;}
    for(unsigned program=0;program<8;++program)for(unsigned test=0;test<216;++test) {
        float stack[64];uint32_t ids[8]={0},visits=0;unsigned sp=0;bool error=false;
        for(unsigned pc=0;pc<cells[program].count;++pc) {
            compute_node n=code[cells[program].root+pc];
            if(n.op==COMPUTE_LEAF) {
                ids[visits++]=n.arg;unsigned divisor=(n.arg>>1)==0?1:(n.arg>>1)==1?6:36;
                float v=0.0f+values[test/divisor%6];if(!isfinite(v)){error=true;break;}
                stack[sp++]=(n.arg&1)?-v:v;
            } else if(n.op==COMPUTE_SKIP_UNION || n.op==COMPUTE_SKIP_INTERSECTION || n.op==COMPUTE_SKIP_REVERSE_DIFFERENCE) {
                float v=stack[sp-1];
                if((n.op==COMPUTE_SKIP_UNION && v<=0) || (n.op==COMPUTE_SKIP_INTERSECTION && v>0) || (n.op==COMPUTE_SKIP_REVERSE_DIFFERENCE && v<0)) {
                    if(n.op==COMPUTE_SKIP_REVERSE_DIFFERENCE)stack[sp-1]=-v;pc+=n.arg;
                }
            } else if(n.op==COMPUTE_COMPLEMENT)stack[sp-1]=-stack[sp-1];
            else {float b=stack[--sp],a=stack[sp-1];stack[sp-1]=n.op==COMPUTE_UNION?fminf(a,b):n.op==COMPUTE_INTERSECTION?fmaxf(a,b):n.op==COMPUTE_DIFFERENCE?fmaxf(a,-b):fmaxf(b,-a);}
        }
        const uint32_t* result=(const uint32_t*)outputs[2].contents+10*(program*216+test);float value;memcpy(&value,result,4);
        if(result[1]!=visits || memcmp(result+2,ids,sizeof(ids)) || error==isfinite(value) || (!error && ((value<=0)!=(stack[0]<=0)))) {
            fprintf(stderr,"Independent CPU trace/classification oracle failed\n");goto done;
        }
    }
    puts("group_trace_cases=1728 map_modes=fallback,hash results=bit-identical leaf_order=identical errors=identical cpu_oracle=passed");rc=0;
 done:for(unsigned i=0;i<5;++i)[inputs[i] release];for(unsigned i=0;i<3;++i){[pipelines[i] release];[outputs[i] release];}return rc;
}
