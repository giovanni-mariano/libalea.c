// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
// SPDX-License-Identifier: MPL-2.0
/* Diagnostic only. Attribute work to the actual selected program structure. */
#ifndef ALEA_REGION_BRIDGE_INCLUDED
#include "../src/compute/compute_metal.m"
#endif
#include "gpu_geometry_fixture.h"
#include "point_work_audit_layout.h"
#ifdef ALEA_GROUP_WORK_AUDIT
#include "group_work_audit_metallib.h"
#include "region_group_storage.h"
#define alea_region_work_audit_metallib alea_group_work_audit_metallib
#define alea_region_work_audit_metallib_len alea_group_work_audit_metallib_len
#endif
#if !defined(ALEA_GROUP_WORK_AUDIT) && !defined(ALEA_REGION_HELPERS_ONLY)
#include "region_work_audit_metallib.h"
#endif
#include "region_program_audit.h"
#ifndef ALEA_REGION_HELPERS_ONLY
static const char* audit_reasons[]={"specialized","no_program","uncertified","outside_domain","unchanged"};
#endif
static double seconds(void) {struct timespec t;clock_gettime(CLOCK_MONOTONIC,&t);return t.tv_sec+t.tv_nsec*1e-9;}
typedef struct {uint32_t root,count,kind,leaves;} region_range;
static int region_compare(const void* a,const void* b) {
    const region_range* x=a,*y=b;
    if(x->root!=y->root)return x->root<y->root?-1:1;
    return x->count<y->count?-1:x->count!=y->count;
}
static uint32_t* region_metadata(metal_state* s,size_t* words,size_t* offset) {
    compute_scene host=s->scene;
    if(s->compact_code)host.compact_nodes=s->nodes.contents;else if(s->nodes)host.nodes=s->nodes.contents;
    size_t capacity=host.cell_count+host.spatial_count;
    if(capacity>(64*1024*1024-256*8*AUDIT_REASONS*AUDIT_FIELDS*4)/32)return NULL;
    region_range* ranges=malloc((capacity?capacity:1)*sizeof(*ranges));if(!ranges)return NULL;
    size_t n=0;
    for(size_t i=0;i<host.cell_count+host.spatial_count;++i) {
        if((i&4095)==0 && alea_interrupted()){free(ranges);return NULL;}
        compute_cell c={0};
        if(i<host.cell_count)c=host.cells[i];else {
            const compute_spatial_program* p=host.spatial+i-host.cell_count;
            if(p->child!=COMPUTE_NONE)continue;
            c.root=p->root;c.count=p->count;
        }
        uint32_t leaves;int result=region_and_program(&host,c,&leaves);
        if(result<0){free(ranges);return NULL;}
        ranges[n++]=(region_range){c.root,c.count,result?(leaves>1?2:1):0,leaves};
    }
    qsort(ranges,n,sizeof(*ranges),region_compare);size_t unique=0;
    for(size_t i=0;i<n;++i)if(!unique || region_compare(ranges+i,ranges+unique-1))ranges[unique++]=ranges[i];
    *offset=1+unique*4;*words=*offset+256*8*AUDIT_REASONS*AUDIT_FIELDS;
    uint32_t* result=calloc(*words,4);if(!result){free(ranges);return NULL;}
    result[0]=(uint32_t)unique;memcpy(result+1,ranges,unique*sizeof(*ranges));
    uint64_t counts[3]={0},leaves[3]={0};
    for(size_t i=0;i<unique;++i){++counts[ranges[i].kind];leaves[ranges[i].kind]+=ranges[i].leaves;}
    printf("programs=%zu general=%llu single=%llu intersection=%llu intersection_leaf_bytes=%llu metadata_bytes=%zu\n",unique,
        (unsigned long long)counts[0],(unsigned long long)counts[1],(unsigned long long)counts[2],(unsigned long long)(4*leaves[2]),*words*4);
    free(ranges);return result;
}
#ifndef ALEA_REGION_HELPERS_ONLY
static int region_write(const uint32_t* table,const char* path) {
    FILE* f=fopen(path,"w");if(!f)return -1;
    const char* kinds[]={"general","single","intersection","unknown"};
    size_t offset=1+4*(size_t)table[0];int unknown=0;
    fprintf(f,"kind,depth,reason,calls,instructions,primitives,reject_instructions,reject_calls"
#ifdef ALEA_GROUP_WORK_AUDIT
        ",group_instructions,group_primitives,group_entries,group_skips,group_combines"
#endif
        "\n");
    for(unsigned k=0;k<4;++k)for(unsigned d=0;d<2;++d)for(unsigned r=0;r<AUDIT_REASONS;++r) {
        uint64_t row[AUDIT_FIELDS]={0};
        for(unsigned shard=0;shard<256;++shard)for(unsigned j=0;j<AUDIT_FIELDS;++j)
            row[j]+=table[offset+((shard*8+k*2+d)*AUDIT_REASONS+r)*AUDIT_FIELDS+j];
        if(!row[AUDIT_CALLS])continue;
#ifdef ALEA_GROUP_WORK_AUDIT
        if(row[AUDIT_GROUP_INSTRUCTIONS]!=row[AUDIT_GROUP_PRIMITIVES]+row[AUDIT_GROUP_SKIPS]+row[AUDIT_GROUP_COMBINES] ||
           row[AUDIT_GROUP_INSTRUCTIONS]>row[AUDIT_INSTRUCTIONS] || row[AUDIT_GROUP_PRIMITIVES]>row[AUDIT_PRIMITIVES] ||
           row[AUDIT_GROUP_ENTRIES]>row[AUDIT_GROUP_PRIMITIVES]){fclose(f);return -1;}
#endif
        if(k==3)unknown=1;
        fprintf(f,"%s,%s,%s",kinds[k],d?"child":"root",audit_reasons[r]);
        for(unsigned j=0;j<AUDIT_FIELDS;++j)fprintf(f,",%" PRIu64,row[j]);fputc('\n',f);
    }
    int bad=ferror(f);bad|=fclose(f)!=0;return bad || unknown?-1:0;
}
static int audit_work(metal_state* s,const double* xyz,size_t n,size_t batch,const char* path,int public_oracle) {
    if(!s->hierarchical || s->point_candidates_enabled || !n || n>16384 || !batch || batch>4096)return -1;
    const char* threading=getenv("ALEA_METAL_JUMP_THREADING");
    if(threading && strcmp(threading,"0")){fprintf(stderr,"Region audit requires unthreaded code\n");return -1;}
    bool selected=getenv("ALEA_BENCH_REGION_SELECT")!=NULL;
    double analysis_start=seconds();
    size_t words=0,offset=0;uint32_t* metadata=region_metadata(s,&words,&offset);
    if(!metadata)return -1;
    printf("analysis_ms=%.6f operation=%s\n",(seconds()-analysis_start)*1000,selected?"selected":"all-owner");

    int rc=-1;id<MTLComputePipelineState> pipeline=nil;
    id<MTLBuffer> counters=nil,profiles=nil,points=nil,output=nil,selection_costs=nil,selection_depths=nil;
#ifdef ALEA_GROUP_WORK_AUDIT
    id<MTLBuffer> groups=nil;
#endif
    uint32_t* totals=NULL;size_t result_bytes=selected?sizeof(compute_selection):sizeof(compute_result);
    void* reference=malloc(n*result_bytes);
    if(!reference || prepare_classify(s))goto done;
    double start=seconds();if(selected?alea_compute_metal_select(s,xyz,n,reference):alea_compute_metal_classify(s,xyz,n,reference))goto done;
    printf("ordinary_ms=%.6f queries=%zu batch=%zu\n",(seconds()-start)*1000,n,batch);
    if(public_oracle && !selected && check_points(xyz,reference,n))goto done;
    NSError* error=nil;dispatch_data_t data=dispatch_data_create(alea_region_work_audit_metallib,alea_region_work_audit_metallib_len,NULL,DISPATCH_DATA_DESTRUCTOR_DEFAULT);
    id<MTLLibrary> library=[s->device newLibraryWithData:data error:&error];dispatch_release(data);
#ifdef ALEA_GROUP_WORK_AUDIT
    if(public_oracle && group_map_self_test(s,library)){[library release];goto done;}
#endif
    MTLFunctionConstantValues* constants=[[MTLFunctionConstantValues alloc] init];bool yes=true,no=false;
    [constants setConstantValue:&yes type:MTLDataTypeBool atIndex:5];
    [constants setConstantValue:&no type:MTLDataTypeBool atIndex:6];
    [constants setConstantValue:&s->compact_code type:MTLDataTypeBool atIndex:7];
    id<MTLFunction> function=[library newFunctionWithName:selected?@"alea_select_hierarchy":@"alea_classify_hierarchy" constantValues:constants error:&error];
    if(function)pipeline=[s->device newComputePipelineStateWithFunction:function error:&error];
    [function release];[constants release];[library release];
    if(!pipeline)goto done;
    counters=[s->device newBufferWithLength:words*4 options:MTLResourceStorageModeShared];
    profiles=[s->device newBufferWithLength:batch*sizeof(compute_point_profile) options:MTLResourceStorageModeShared];
    points=[s->device newBufferWithLength:batch*3*sizeof(float) options:MTLResourceStorageModeShared];
    output=[s->device newBufferWithLength:batch*result_bytes options:MTLResourceStorageModeShared];
    if(!counters || !profiles || !points || !output)goto done;
    if(selected) {
        selection_costs=[s->device newBufferWithLength:s->scene.cell_count*SELECT_COUNTER_COUNT*4 options:MTLResourceStorageModeShared];
        selection_depths=[s->device newBufferWithLength:COMPUTE_HIER_STACK*SELECT_COUNTER_COUNT*4 options:MTLResourceStorageModeShared];
        if(!selection_costs || !selection_depths)goto done;
    }
    totals=counters.contents;memcpy(totals,metadata,words*4);free(metadata);metadata=NULL;
#ifdef ALEA_GROUP_WORK_AUDIT
    groups=group_buffer(s,totals);if(!groups)goto done;
#endif
    uint64_t expected[AUDIT_FIELDS]={0};start=seconds();
    for(size_t first=0;first<n;first+=batch) {
        @autoreleasepool {
            uint32_t count=(uint32_t)MIN(batch,n-first),root=s->scene.point_root;
            float* converted=points.contents;
            for(size_t i=0;i<3*count;++i)converted[i]=(float)xyz[3*first+i];
            id<MTLCommandBuffer> command=[s->queue commandBuffer];
            id<MTLComputeCommandEncoder> e=[command computeCommandEncoder];if(!command || !e)goto done;
            [e setComputePipelineState:pipeline];
#ifdef ALEA_GROUP_WORK_AUDIT
            [e setBuffer:groups offset:0 atIndex:21];
#endif
            id<MTLBuffer> buffers[]={s->nodes,s->cells,points,output,nil,nil,s->primitives,s->transforms,
                s->point_cells,s->point_lattices,s->point_fills,s->spatial,profiles};
            for(unsigned i=0;i<13;++i)if(buffers[i])[e setBuffer:buffers[i] offset:0 atIndex:i];
            [e setBytes:&root length:4 atIndex:4];[e setBytes:&count length:4 atIndex:5];
            if(selected) {
                memset(selection_costs.contents,0,selection_costs.length);memset(selection_depths.contents,0,selection_depths.length);
                [e setBuffer:selection_costs offset:0 atIndex:14];[e setBuffer:selection_depths offset:0 atIndex:15];
                [e setBuffer:counters offset:0 atIndex:16];
            } else [e setBuffer:counters offset:0 atIndex:14];
            [e dispatchThreads:MTLSizeMake(count,1,1) threadsPerThreadgroup:MTLSizeMake(MIN((NSUInteger)128,pipeline.maxTotalThreadsPerThreadgroup),1,1)];
            [e endEncoding];[command commit];[command waitUntilCompleted];
            if(command.status!=MTLCommandBufferStatusCompleted){fprintf(stderr,"Audit command: %s\n",command.error.localizedDescription.UTF8String);goto done;}
            if(memcmp(output.contents,(const char*)reference+first*result_bytes,count*result_bytes)){fprintf(stderr,"Audit changed query outputs\n");goto done;}
            const uint32_t* values=counters.contents;
            for(size_t i=offset;i<words;++i)if(values[i]==UINT32_MAX){fprintf(stderr,"Audit counter saturated\n");goto done;}
            if(selected) {
                const uint32_t* costs=selection_costs.contents;
                for(size_t i=0;i<s->scene.cell_count;++i) {
                    const uint32_t* row=costs+i*SELECT_COUNTER_COUNT;
                    for(unsigned k=0;k<SELECT_COUNTER_COUNT;++k)if(row[k]==UINT32_MAX)goto done;
                    expected[AUDIT_CALLS]+=row[SELECT_CALLS]-row[SELECT_LATTICES];
                    expected[AUDIT_INSTRUCTIONS]+=row[SELECT_INSTRUCTIONS];expected[AUDIT_PRIMITIVES]+=row[SELECT_PRIMITIVES];
                    expected[AUDIT_REJECT_CALLS]+=row[SELECT_REJECTED];
                }
            } else {
                const compute_point_profile* p=profiles.contents;
            for(size_t i=0;i<count;++i) {
                for(int j=0;j<POINT_COUNTER_COUNT;++j)if(p[i].v[j]==UINT32_MAX)goto done;
                expected[AUDIT_CALLS]+=p[i].v[POINT_CELLS]-p[i].v[POINT_LATTICES];
                expected[AUDIT_INSTRUCTIONS]+=(uint64_t)p[i].v[POINT_ORIGINAL_INSTRUCTIONS]+p[i].v[POINT_SPATIAL_INSTRUCTIONS];
                expected[AUDIT_PRIMITIVES]+=(uint64_t)p[i].v[POINT_ORIGINAL_PRIMITIVES]+p[i].v[POINT_SPATIAL_PRIMITIVES];
                expected[AUDIT_REJECT_INSTRUCTIONS]+=p[i].v[POINT_REJECTED_INSTRUCTIONS];
                expected[AUDIT_REJECT_CALLS]+=p[i].v[POINT_REJECTED];
            }
            }
        }
    }
    printf("instrumented_ms=%.6f counter_storage_mib=%.3f\n",(seconds()-start)*1000,words*4/1048576.);
    uint64_t actual[AUDIT_FIELDS]={0};
    for(unsigned r=0;r<AUDIT_REASONS;++r) {
        uint64_t sum[AUDIT_FIELDS]={0};
        for(size_t i=0;i<256*8;++i)for(unsigned k=0;k<AUDIT_FIELDS;++k)sum[k]+=totals[offset+(i*AUDIT_REASONS+r)*AUDIT_FIELDS+k];
        printf("reason=%s",audit_reasons[r]);for(unsigned k=0;k<AUDIT_FIELDS;++k){printf(" field%u=%" PRIu64,k,sum[k]);actual[k]+=sum[k];}puts("");
    }
    for(unsigned k=0;k<=AUDIT_REJECT_CALLS;++k) {
        /* Selection's existing profiler has no rejected-instruction subtotal. */
        if(selected && k==AUDIT_REJECT_INSTRUCTIONS)continue;
        if(actual[k]!=expected[k]){fprintf(stderr,"Independent profiles disagree\n");goto done;}
    }
#ifdef ALEA_GROUP_WORK_AUDIT
    if(public_oracle && getenv("ALEA_BENCH_CSG_GROUPS") && !getenv("ALEA_BENCH_CSG_CHAIN")) {
        uint64_t inside=0;
        for(size_t i=0;i<n;++i) {
            double x=xyz[3*i],y=xyz[3*i+1];
            x-=2*floor((x+8)/2)-7;y-=2*floor((y+8)/2)-7;
            inside+=x*x+y*y+xyz[3*i+2]*xyz[3*i+2]<=.75*.75;
        }
        uint64_t oracle[]={4*n,2*n,2*n-inside,2*n-inside,inside};
        for(unsigned k=0;k<5;++k)if(actual[AUDIT_GROUP_INSTRUCTIONS+k]!=oracle[k]) {
            fprintf(stderr,"Independent group-work oracle failed\n");goto done;
        }
        puts("group_work_oracle=passed");
    }
#endif
    if(region_write(totals,path))goto done;
    printf("outputs=identical totals=reconciled saturation=0 reconciled_fields=%s\n",selected?"calls,instructions,primitives,reject_calls":"all");rc=0;
 done:
#ifdef ALEA_GROUP_WORK_AUDIT
    [groups release];
#endif
    free(metadata);if(rc)fprintf(stderr,"Work audit failed\n");[pipeline release];[counters release];[profiles release];[points release];[output release];[selection_costs release];[selection_depths release];free(reference);return rc;
}
int main(int argc,char** argv) {
    if(argc!=2){fprintf(stderr,"Usage: %s OUTPUT.csv\n",argv[0]);return 2;}
    @autoreleasepool {
        (void)check_render;setbuf(stdout,NULL);alea_system_t* sys=scene();metal_state* s=NULL;int rc=1;
        const size_t n=128*128;double* xyz=malloc(3*n*sizeof(*xyz));
        if(!sys || !xyz || !(s=alea_compute_metal_create(sys)))goto done;
        for(size_t i=0;i<n;++i){xyz[3*i]=-8+((i%128)+.5)/8;xyz[3*i+1]=-8+((i/128)+.5)/8;xyz[3*i+2]=0;}
        if(audit_work(s,xyz,n,4096,argv[1],1))goto done;rc=0;puts("result=PASS");
 done:free(xyz);alea_compute_metal_destroy(s);alea_destroy(sys);return rc;
    }
}

#endif /* ALEA_REGION_HELPERS_ONLY */
