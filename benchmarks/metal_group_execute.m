// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
// SPDX-License-Identifier: MPL-2.0
/* Bounded native experiment. One context per process; no production switch. */
#import <Foundation/Foundation.h>
#import <Metal/Metal.h>
#include "alea.h"
#include "compute/compute_internal.h"
static int execute_pack(const alea_system_t*,compute_scene*);
#define alea_compute_pack execute_pack
#include "../src/compute/compute_metal.m"
#undef alea_compute_pack
#define ALEA_REGION_BRIDGE_INCLUDED 1
#define ALEA_REGION_HELPERS_ONLY 1
#include "metal_region_work.m"
#include "region_group_storage.h"
#include "group_execute_metallib.h"
static id<MTLBuffer> execute_maps,execute_leaves;
static int execute_enabled;
static double execute_prepare_ms;
static void execute_release(void) {[execute_maps release];[execute_leaves release];execute_maps=execute_leaves=nil;}
static int execute_pack(const alea_system_t* sys,compute_scene* scene) {
    (void)group_map_self_test;(void)check_render;
    const char* mode=getenv("ALEA_BENCH_GROUP_EXECUTE");
    if(mode && strcmp(mode,"0") && strcmp(mode,"1")) {
        alea_set_error_detail(ALEA_ERR_INVALID_ARG,"ALEA_BENCH_GROUP_EXECUTE must be 0 or 1");return -1;
    }
    const char* lookup=getenv("ALEA_BENCH_GROUP_HASH");
    if(lookup && strcmp(lookup,"0") && strcmp(lookup,"1")) {
        alea_set_error_detail(ALEA_ERR_INVALID_ARG,"ALEA_BENCH_GROUP_HASH must be 0 or 1");return -1;
    }
    int hashed=lookup && !strcmp(lookup,"1");
    execute_enabled=mode && !strcmp(mode,"1");execute_prepare_ms=0;
    int rc=alea_compute_pack(sys,scene);if(rc || !execute_enabled)return rc;
    double start=seconds();id<MTLDevice> device=MTLCreateSystemDefaultDevice();
    metal_state temporary={.device=device,.scene=*scene};
    size_t words=0,offset=0;uint32_t* metadata=region_metadata(&temporary,&words,&offset);
    id<MTLBuffer> old=metadata?group_buffer(&temporary,metadata):nil;
    id<MTLBuffer> maps=nil,leaves=nil;size_t count=0,leaf_count=0;bool unavailable=false;
    if(!old){unavailable=true;goto fail;}
    const uint32_t* source=old.contents;size_t programs=source[0],base=1+4*programs;
    count=programs?(size_t)source[3+4*(programs-1)]+source[4+4*(programs-1)]:0;
    size_t hash_capacity=0;
    if(hashed){hash_capacity=1;while(hash_capacity<2*programs)hash_capacity*=2;}
    size_t bytes=16+programs*16+count*12+hash_capacity*4;
    if(bytes+old.length+words*4>64*1024*1024){unavailable=true;goto fail;}
    maps=[device newBufferWithLength:bytes options:MTLResourceStorageModeShared];if(!maps)goto fail;
    uint32_t* target=maps.contents;target[0]=(uint32_t)programs;target[1]=(uint32_t)count;target[2]=hashed?(uint32_t)(hash_capacity-1):0;target[3]=0;
    memcpy(target+4,source+1,programs*16);
    for(size_t i=0;i<count;++i) {
        uint32_t begin=source[base+2*i],end=source[base+2*i+1],length=end-begin;
        if(length<4 || length%3!=1 || (i%4096==0 && alea_interrupted()))goto fail;
        target[4+4*programs+3*i]=begin;target[5+4*programs+3*i]=end;
        target[6+4*programs+3*i]=(uint32_t)leaf_count;leaf_count+=(length+2)/3;
    }
    if(hashed) {
        uint32_t* slots=target+4+4*programs+3*count;memset(slots,0,hash_capacity*4);
        for(size_t i=0;i<programs;++i) {
            uint32_t h=(target[4+4*i]*0x9e3779b9u)^(target[5+4*i]*0x85ebca6bu);h=(h^(h>>16))&(uint32_t)(hash_capacity-1);
            while(slots[h])h=(h+1)&(uint32_t)(hash_capacity-1);slots[h]=(uint32_t)i+1;
        }
    }
    [old release];old=nil;free(metadata);metadata=NULL;
    if(bytes+MAX((size_t)4,leaf_count*4)>64*1024*1024){unavailable=true;goto fail;}
    leaves=[device newBufferWithLength:MAX((size_t)4,leaf_count*4) options:MTLResourceStorageModeShared];if(!leaves)goto fail;
    uint32_t* pool=leaves.contents;size_t written=0;
    for(size_t i=0;i<count;++i) {
        if(!(i%4096) && alea_interrupted())goto fail;
        uint32_t begin=target[4+4*programs+3*i],end=target[5+4*programs+3*i];
        if(written!=target[6+4*programs+3*i])goto fail;
        for(uint32_t pc=begin;pc<end;++pc)if(scene->nodes[pc].op==COMPUTE_LEAF)pool[written++]=scene->nodes[pc].arg;
    }
    if(written!=leaf_count)goto fail;
    execute_maps=maps;execute_leaves=leaves;maps=leaves=nil;
    execute_prepare_ms=(seconds()-start)*1000;
    printf("group_prepare_ms=%.6f group_count=%zu group_bytes=%zu leaf_bytes=%zu hash=%d\n",execute_prepare_ms,count,execute_maps.length,execute_leaves.length,hashed);
    [device release];return 0;
 fail:free(metadata);[old release];[maps release];[leaves release];[device release];
    if(unavailable && !alea_interrupted()) {
        execute_enabled=0;puts("group_fallback=unavailable_or_budget generic=retained");return 0;
    }
    alea_compute_scene_free(scene);alea_set_error_detail(alea_interrupted()?ALEA_ERR_INTERRUPTED:ALEA_ERR_UNSUPPORTED,"Group preparation failed");return -1;
}
static id<MTLComputePipelineState> execute_pipeline(metal_state* s,id<MTLLibrary> library,const char* name,bool enabled,bool trace,bool profile) {
    NSError* error=nil;MTLFunctionConstantValues* c=[[MTLFunctionConstantValues alloc] init];bool no=false;
    [c setConstantValue:&profile type:MTLDataTypeBool atIndex:5];[c setConstantValue:&no type:MTLDataTypeBool atIndex:6];
    [c setConstantValue:&s->compact_code type:MTLDataTypeBool atIndex:7];
    [c setConstantValue:&enabled type:MTLDataTypeBool atIndex:9];[c setConstantValue:&trace type:MTLDataTypeBool atIndex:10];
    id<MTLFunction> f=[library newFunctionWithName:[NSString stringWithUTF8String:name] constantValues:c error:&error];
    id<MTLComputePipelineState> p=f?[s->device newComputePipelineStateWithFunction:f error:&error]:nil;
    if(!p)fprintf(stderr,"Group pipeline: %s\n",error.localizedDescription.UTF8String);
    [f release];[c release];return p;
}
static int execute_trace_test(metal_state* s,id<MTLLibrary> library) {
    const compute_node patterns[6][10]={
        {{COMPUTE_LEAF,0},{COMPUTE_SKIP_INTERSECTION,5},{COMPUTE_LEAF,2},{COMPUTE_SKIP_INTERSECTION,2},{COMPUTE_LEAF,4},{COMPUTE_INTERSECTION,0},{COMPUTE_INTERSECTION,0}},
        {{COMPUTE_LEAF,0},{COMPUTE_SKIP_INTERSECTION,5},{COMPUTE_LEAF,2},{COMPUTE_SKIP_INTERSECTION,2},{COMPUTE_LEAF,4},{COMPUTE_INTERSECTION,0},{COMPUTE_INTERSECTION,0},{COMPUTE_COMPLEMENT,0}},
        {{COMPUTE_LEAF,0},{COMPUTE_SKIP_INTERSECTION,2},{COMPUTE_LEAF,2},{COMPUTE_INTERSECTION,0},{COMPUTE_SKIP_UNION,2},{COMPUTE_LEAF,4},{COMPUTE_UNION,0}},
        {{COMPUTE_LEAF,4},{COMPUTE_SKIP_INTERSECTION,5},{COMPUTE_LEAF,0},{COMPUTE_SKIP_INTERSECTION,2},{COMPUTE_LEAF,2},{COMPUTE_INTERSECTION,0},{COMPUTE_DIFFERENCE,0}},
        {{COMPUTE_LEAF,0},{COMPUTE_SKIP_INTERSECTION,2},{COMPUTE_LEAF,2},{COMPUTE_INTERSECTION,0},{COMPUTE_SKIP_REVERSE_DIFFERENCE,2},{COMPUTE_LEAF,4},{COMPUTE_REVERSE_DIFFERENCE,0}},
        {{COMPUTE_LEAF,0},{COMPUTE_SKIP_INTERSECTION,2},{COMPUTE_LEAF,2},{COMPUTE_INTERSECTION,0},{COMPUTE_SKIP_UNION,5},{COMPUTE_LEAF,1},{COMPUTE_SKIP_INTERSECTION,2},{COMPUTE_LEAF,4},{COMPUTE_INTERSECTION,0},{COMPUTE_UNION,0}}};
    const uint32_t lengths[]={7,8,7,7,7,10};compute_node code[64];compute_cell cells[8]={0};
    uint32_t map[128]={8},refs[64];size_t used=0,nrefs=0,ngroups=0;
    group_storage storage={.limit=1024};compute_scene scene={.nodes=code,.primitive_count=3};
    for(unsigned k=0;k<6;++k){cells[k].root=(uint32_t)used;cells[k].count=lengths[k];memcpy(code+used,patterns[k],lengths[k]*sizeof(*code));used+=lengths[k];}
    cells[6]=(compute_cell){.root=0,.count=1};cells[7]=(compute_cell){.root=2,.count=4};
    scene.node_count=used;
    for(unsigned k=0;k<8;++k) {
        storage.count=0;if(region_group_program(&scene,cells[k],group_append,&storage)){free(storage.data);return -1;}
        if(storage.count)qsort(storage.data+1,storage.count,8,group_compare);
        map[4+4*k]=cells[k].root;map[5+4*k]=cells[k].count;map[6+4*k]=(uint32_t)ngroups;map[7+4*k]=(uint32_t)storage.count;
        for(size_t j=0;j<storage.count;++j) {
            uint32_t begin=storage.data[1+2*j],end=storage.data[2+2*j];
            map[36+3*ngroups]=begin;map[37+3*ngroups]=end;map[38+3*ngroups]=(uint32_t)nrefs;++ngroups;
            for(uint32_t pc=begin;pc<end;++pc)if(code[pc].op==COMPUTE_LEAF)refs[nrefs++]=code[pc].arg;
        }
    }
    free(storage.data);map[1]=(uint32_t)ngroups;
    qsort(map+4,8,16,region_compare);
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
        pipelines[i]=execute_pipeline(s,library,"alea_group_execute_probe",i!=0,true,false);
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
    puts("group_trace_cases=1728 map_modes=binary,hash results=bit-identical leaf_order=identical errors=identical cpu_oracle=passed");rc=0;
 done:for(unsigned i=0;i<5;++i)[inputs[i] release];for(unsigned i=0;i<3;++i){[pipelines[i] release];[outputs[i] release];}return rc;
}
static int execute_profile(metal_state* s,id<MTLLibrary> library,const double* xyz,size_t n,const void* reference) {
    id<MTLComputePipelineState> pipeline=execute_pipeline(s,library,"alea_classify_hierarchy",execute_enabled,false,true);
    id<MTLBuffer> points=[s->device newBufferWithLength:128*12 options:MTLResourceStorageModeShared];
    id<MTLBuffer> results=[s->device newBufferWithLength:128*sizeof(compute_result) options:MTLResourceStorageModeShared];
    id<MTLBuffer> profiles=[s->device newBufferWithLength:128*sizeof(compute_point_profile) options:MTLResourceStorageModeShared];
    uint64_t primitives=0;int rc=-1;
    if(!pipeline || !points || !results || !profiles || execute_maps.length+execute_leaves.length+profiles.length>64*1024*1024)goto done;
    for(size_t first=0;first<n;first+=128) @autoreleasepool {
        uint32_t count=(uint32_t)MIN((size_t)128,n-first),root=s->scene.point_root;float* p=points.contents;
        for(size_t i=0;i<3*count;++i)p[i]=(float)xyz[3*first+i];
        id<MTLCommandBuffer> command=[s->queue commandBuffer];id<MTLComputeCommandEncoder> e=[command computeCommandEncoder];if(!command || !e)goto done;
        [e setComputePipelineState:pipeline];id<MTLBuffer> buffers[]={s->nodes,s->cells,points,results,nil,nil,s->primitives,s->transforms,s->point_cells,s->point_lattices,s->point_fills,s->spatial,profiles};
        for(unsigned i=0;i<13;++i)if(buffers[i])[e setBuffer:buffers[i] offset:0 atIndex:i];
        [e setBytes:&root length:4 atIndex:4];[e setBytes:&count length:4 atIndex:5];
        [e setBuffer:execute_enabled?execute_maps:points offset:0 atIndex:17];[e setBuffer:execute_enabled?execute_leaves:points offset:0 atIndex:18];
        [e dispatchThreads:MTLSizeMake(count,1,1) threadsPerThreadgroup:MTLSizeMake(128,1,1)];[e endEncoding];[command commit];[command waitUntilCompleted];
        if(command.status!=MTLCommandBufferStatusCompleted || memcmp(results.contents,(const char*)reference+first*sizeof(compute_result),count*sizeof(compute_result)))goto done;
        const compute_point_profile* cost=profiles.contents;
        for(size_t i=0;i<count;++i) {
            for(unsigned k=0;k<POINT_COUNTER_COUNT;++k)if(cost[i].v[k]==UINT32_MAX)goto done;
            primitives+=(uint64_t)cost[i].v[POINT_ORIGINAL_PRIMITIVES]+cost[i].v[POINT_SPATIAL_PRIMITIVES];
        }
    }
    printf("diagnostic_points=%zu diagnostic_primitives=%llu outputs=identical saturation=0\n",n,(unsigned long long)primitives);rc=0;
 done:[pipeline release];[points release];[results release];[profiles release];return rc;
}
static int execute_queries(metal_state* s,const double* xyz,size_t n,size_t batch,const char* path,int public_oracle) {
    int rc=-1;bool selected=getenv("ALEA_BENCH_REGION_SELECT")!=NULL;
    size_t stride=selected?sizeof(compute_selection):sizeof(compute_result);
    void* reference=malloc(n*stride);id<MTLBuffer> points=nil,output=nil;id<MTLComputePipelineState> pipeline=nil;id<MTLLibrary> library=nil;
    if(s->point_profile_path || s->select_profile_path || s->profile_prefix || s->root_reuse_path)goto done;
    if(execute_enabled && execute_maps.device.registryID!=s->device.registryID)goto done;
    if(context_buffer_bytes(s)+execute_maps.length+execute_leaves.length+batch*(12+stride)>1024ull*1024*1024)goto done;
    if(!reference || !s->hierarchical || s->point_candidates_enabled || n>16384 || !n || !batch || batch>4096)goto done;
    if(selected?alea_compute_metal_select(s,xyz,n,reference):alea_compute_metal_classify(s,xyz,n,reference))goto done;
    if(public_oracle && !selected && check_points(xyz,reference,n))goto done;
    NSError* error=nil;dispatch_data_t data=dispatch_data_create(alea_group_execute_metallib,alea_group_execute_metallib_len,NULL,DISPATCH_DATA_DESTRUCTOR_DEFAULT);
    library=[s->device newLibraryWithData:data error:&error];dispatch_release(data);
    if(public_oracle && execute_trace_test(s,library))goto done;
    pipeline=execute_pipeline(s,library,selected?"alea_select_hierarchy":"alea_classify_hierarchy",execute_enabled,false,false);
    points=[s->device newBufferWithLength:batch*12 options:MTLResourceStorageModeShared];output=[s->device newBufferWithLength:batch*stride options:MTLResourceStorageModeShared];
    if(!pipeline || !points || !output)goto done;
    const char* profile=getenv("ALEA_BENCH_GROUP_PROFILE");
    bool diagnostic=profile && !strcmp(profile,"1");
    for(unsigned pass=0;pass<(diagnostic?1u:6u);++pass) {
        double start=seconds();
        for(size_t first=0;first<n;first+=batch) @autoreleasepool {
            uint32_t count=(uint32_t)MIN(batch,n-first),root=s->scene.point_root;float* p=points.contents;
            for(size_t i=0;i<3*count;++i)p[i]=(float)xyz[3*first+i];
            id<MTLCommandBuffer> command=[s->queue commandBuffer];id<MTLComputeCommandEncoder> e=[command computeCommandEncoder];if(!command || !e)goto done;
            [e setComputePipelineState:pipeline];
            id<MTLBuffer> buffers[]={s->nodes,s->cells,points,output,nil,nil,s->primitives,s->transforms,s->point_cells,s->point_lattices,s->point_fills,s->spatial};
            for(unsigned i=0;i<12;++i)if(buffers[i])[e setBuffer:buffers[i] offset:0 atIndex:i];
            [e setBytes:&root length:4 atIndex:4];[e setBytes:&count length:4 atIndex:5];
            [e setBuffer:execute_enabled?execute_maps:points offset:0 atIndex:17];[e setBuffer:execute_enabled?execute_leaves:points offset:0 atIndex:18];
            [e dispatchThreads:MTLSizeMake(count,1,1) threadsPerThreadgroup:MTLSizeMake(MIN((NSUInteger)128,pipeline.maxTotalThreadsPerThreadgroup),1,1)];
            [e endEncoding];[command commit];[command waitUntilCompleted];
            if(command.status!=MTLCommandBufferStatusCompleted){fprintf(stderr,"Group command: %s\n",command.error.localizedDescription.UTF8String);goto done;}
            if(memcmp(output.contents,(char*)reference+first*stride,count*stride)){fprintf(stderr,"Group query output mismatch\n");goto done;}
        }
        printf("pass=%u warm=%d execute=%d selected=%d wall_ms=%.6f points=%zu batch=%zu\n",pass,pass!=0,execute_enabled,selected,(seconds()-start)*1000,n,batch);
    }
    if(diagnostic){if(selected || execute_profile(s,library,xyz,MIN(n,(size_t)512),reference))goto done;}
    if(path){FILE* f=fopen(path,"wb");if(!f)goto done;int bad=fwrite(reference,stride,n,f)!=n;bad|=fclose(f)!=0;if(bad)goto done;}
    puts("outputs=identical profiling=disabled");rc=0;
 done:free(reference);[pipeline release];[points release];[output release];[library release];return rc;
}
static int execute_lifecycle_test(const alea_system_t* sys) {
    const char* mode=getenv("ALEA_BENCH_GROUP_EXECUTE"),*hash=getenv("ALEA_BENCH_GROUP_HASH");
    char* saved_mode=mode?strdup(mode):NULL;char* saved_hash=hash?strdup(hash):NULL;int rc=-1;
    if((mode && !saved_mode) || (hash && !saved_hash))goto done;
    setenv("ALEA_BENCH_GROUP_EXECUTE","invalid",1);
    void* context=alea_compute_metal_create(sys);if(context){alea_compute_metal_destroy(context);goto done;}
    if(alea_error_code()!=ALEA_ERR_INVALID_ARG)goto done;
    if(saved_mode)setenv("ALEA_BENCH_GROUP_EXECUTE",saved_mode,1);else unsetenv("ALEA_BENCH_GROUP_EXECUTE");
    setenv("ALEA_BENCH_GROUP_HASH","invalid",1);
    context=alea_compute_metal_create(sys);if(context){alea_compute_metal_destroy(context);goto done;}
    if(alea_error_code()!=ALEA_ERR_INVALID_ARG)goto done;
    if(saved_hash)setenv("ALEA_BENCH_GROUP_HASH",saved_hash,1);else unsetenv("ALEA_BENCH_GROUP_HASH");
    alea_interrupt();context=alea_compute_metal_create(sys);alea_clear_interrupt();
    if(context){alea_compute_metal_destroy(context);goto done;}
    if(alea_error_code()!=ALEA_ERR_INTERRUPTED)goto done;rc=0;
 done:
    if(saved_mode)setenv("ALEA_BENCH_GROUP_EXECUTE",saved_mode,1);else unsetenv("ALEA_BENCH_GROUP_EXECUTE");
    if(saved_hash)setenv("ALEA_BENCH_GROUP_HASH",saved_hash,1);else unsetenv("ALEA_BENCH_GROUP_HASH");
    free(saved_mode);free(saved_hash);return rc;
}
int main(int argc,char** argv) {
    if(argc!=2){fprintf(stderr,"Usage: %s OUTPUT.bin\n",argv[0]);return 2;}
    @autoreleasepool {
        setbuf(stdout,NULL);alea_system_t* sys=scene();metal_state* s=NULL;int rc=1;size_t n=16384;double* xyz=malloc(n*3*sizeof(*xyz));
        if(!sys || !xyz || execute_lifecycle_test(sys) || !(s=alea_compute_metal_create(sys)))goto done;
        puts("group_lifecycle=invalid-switch,interruption,retry-passed");
        for(size_t i=0;i<n;++i){xyz[3*i]=-8+((i%128)+.5)/8;xyz[3*i+1]=-8+((i/128)+.5)/8;xyz[3*i+2]=0;}
        if(execute_queries(s,xyz,n,4096,argv[1],1))goto done;rc=0;puts("result=PASS");
 done:free(xyz);alea_compute_metal_destroy(s);execute_release();alea_destroy(sys);return rc;
    }
}
