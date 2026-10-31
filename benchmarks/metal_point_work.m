// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
// SPDX-License-Identifier: MPL-2.0
/* Diagnostic only. Per-cell output can identify private model structure. */
#include "../src/compute/compute_metal.m"
#include "gpu_geometry_fixture.h"
#include "point_work_audit_layout.h"
#include "point_work_audit_metallib.h"
static const char* audit_reasons[]={"specialized","no_program","uncertified","outside_domain","unchanged"};
static double seconds(void) {struct timespec t;clock_gettime(CLOCK_MONOTONIC,&t);return t.tv_sec+t.tv_nsec*1e-9;}
static int write_work(metal_state* s,const uint32_t* table,const char* path) {
    FILE* f=fopen(path,"w");if(!f)return -1;
    fprintf(f,"cell_index,cell_id,material_id,program_instructions,kind,point_certificate,skip_ops,threadable_ops,reason,calls,instructions,primitives,reject_instructions,reject_calls\n");
    size_t budget=32u*1024u*1024u;const metal_point_cell* links=s->point_cells.contents;
    compute_scene host=s->scene;
    if(s->compact_code)host.compact_nodes=s->nodes.contents;else host.nodes=s->nodes.contents;
    const char* mode=getenv("ALEA_METAL_POINT_SPATIAL");
    for(size_t i=0;i<s->scene.cell_count;++i) {
        compute_cell c=s->scene.cells[i];const char* certificate="no_program";
        if(c.spatial!=COMPUTE_NONE) {
            if(mode && !strcmp(mode,"0"))certificate="disabled";
            else if(c.lattice!=COMPUTE_NONE)certificate="lattice";
            else if(c.count>budget)certificate="budget";
            else {budget-=c.count;certificate=links[i].spatial!=COMPUTE_NONE?"certified":"finite_unproven";}
        }
        const char* kind=c.lattice!=COMPUTE_NONE?"lattice":links[i].child!=COMPUTE_NONE?"container":"terminal";
        int visited=0;for(size_t r=0;r<AUDIT_REASONS;++r)visited|=table[(i*AUDIT_REASONS+r)*AUDIT_FIELDS+AUDIT_CALLS]!=0;
        if(!visited)continue;
        size_t skips=0,threadable=0;
        for(uint32_t pc=0;pc<c.count;++pc) {
            compute_node n=alea_compute_node_at(&host,(size_t)c.root+pc);
            skips+=n.op==COMPUTE_SKIP_UNION || n.op==COMPUTE_SKIP_INTERSECTION || n.op==COMPUTE_SKIP_REVERSE_DIFFERENCE;
            uint64_t target=(uint64_t)pc+n.arg+1;
            if((n.op==COMPUTE_SKIP_UNION || n.op==COMPUTE_SKIP_INTERSECTION) && target<c.count &&
               alea_compute_node_at(&host,(size_t)c.root+target).op==n.op)++threadable;
        }
        for(size_t r=0;r<AUDIT_REASONS;++r) {
            const uint32_t* row=table+(i*AUDIT_REASONS+r)*AUDIT_FIELDS;if(!row[AUDIT_CALLS])continue;
            fprintf(f,"%zu,%d,%d,%u,%s,%s,%zu,%zu,%s",i,c.cell_id,c.material_id,c.count,kind,certificate,skips,threadable,audit_reasons[r]);
            for(size_t k=0;k<AUDIT_FIELDS;++k)fprintf(f,",%" PRIu64,(uint64_t)row[k]);fputc('\n',f);
        }
    }
    int bad=ferror(f);bad|=fclose(f)!=0;return bad?-1:0;
}
static int audit_work(metal_state* s,const double* xyz,size_t n,size_t batch,const char* path,int public_oracle) {
    if(!s->hierarchical || s->point_candidates_enabled || !n || n>16384 || !batch || batch>4096)return -1;
    size_t words=s->scene.cell_count*AUDIT_REASONS*AUDIT_FIELDS;
    if(words>64*1024*1024/4){fprintf(stderr,"Audit counters exceed 64 MiB\n");return -1;}
    int rc=-1;id<MTLComputePipelineState> pipeline=nil;
    id<MTLBuffer> counters=nil,profiles=nil,points=nil,output=nil;
    uint32_t* totals=NULL;alea_compute_result_t* reference=malloc(n*sizeof(*reference));
    if(!reference || prepare_classify(s))goto done;
    double start=seconds();if(alea_compute_metal_classify(s,xyz,n,reference))goto done;
    printf("ordinary_ms=%.6f queries=%zu batch=%zu\n",(seconds()-start)*1000,n,batch);
    if(public_oracle && check_points(xyz,reference,n))goto done;
    NSError* error=nil;dispatch_data_t data=dispatch_data_create(alea_point_work_audit_metallib,alea_point_work_audit_metallib_len,NULL,DISPATCH_DATA_DESTRUCTOR_DEFAULT);
    id<MTLLibrary> library=[s->device newLibraryWithData:data error:&error];dispatch_release(data);
    MTLFunctionConstantValues* constants=[[MTLFunctionConstantValues alloc] init];bool yes=true,no=false;
    [constants setConstantValue:&yes type:MTLDataTypeBool atIndex:5];
    [constants setConstantValue:&no type:MTLDataTypeBool atIndex:6];
    [constants setConstantValue:&s->compact_code type:MTLDataTypeBool atIndex:7];
    id<MTLFunction> function=[library newFunctionWithName:@"alea_classify_hierarchy" constantValues:constants error:&error];
    if(function)pipeline=[s->device newComputePipelineStateWithFunction:function error:&error];
    [function release];[constants release];[library release];
    if(!pipeline)goto done;
    counters=[s->device newBufferWithLength:words*4 options:MTLResourceStorageModeShared];
    profiles=[s->device newBufferWithLength:batch*sizeof(compute_point_profile) options:MTLResourceStorageModeShared];
    points=[s->device newBufferWithLength:batch*3*sizeof(float) options:MTLResourceStorageModeShared];
    output=[s->device newBufferWithLength:batch*sizeof(compute_result) options:MTLResourceStorageModeShared];
    if(!counters || !profiles || !points || !output)goto done;
    totals=counters.contents;memset(totals,0,words*4);
    uint64_t expected[AUDIT_FIELDS]={0};start=seconds();
    for(size_t first=0;first<n;first+=batch) {
        @autoreleasepool {
            uint32_t count=(uint32_t)MIN(batch,n-first),root=s->scene.point_root;
            float* converted=points.contents;
            for(size_t i=0;i<3*count;++i)converted[i]=(float)xyz[3*first+i];
            id<MTLCommandBuffer> command=[s->queue commandBuffer];
            id<MTLComputeCommandEncoder> e=[command computeCommandEncoder];if(!command || !e)goto done;
            [e setComputePipelineState:pipeline];
            id<MTLBuffer> buffers[]={s->nodes,s->cells,points,output,nil,nil,s->primitives,s->transforms,
                s->point_cells,s->point_lattices,s->point_fills,s->spatial,profiles};
            for(unsigned i=0;i<13;++i)if(buffers[i])[e setBuffer:buffers[i] offset:0 atIndex:i];
            [e setBytes:&root length:4 atIndex:4];[e setBytes:&count length:4 atIndex:5];
            [e setBuffer:counters offset:0 atIndex:14];
            [e dispatchThreads:MTLSizeMake(count,1,1) threadsPerThreadgroup:MTLSizeMake(MIN((NSUInteger)128,pipeline.maxTotalThreadsPerThreadgroup),1,1)];
            [e endEncoding];[command commit];[command waitUntilCompleted];
            if(command.status!=MTLCommandBufferStatusCompleted){fprintf(stderr,"Audit command: %s\n",command.error.localizedDescription.UTF8String);goto done;}
            if(memcmp(output.contents,reference+first,count*sizeof(*reference))){fprintf(stderr,"Audit changed query outputs\n");goto done;}
            const uint32_t* values=counters.contents;
            for(size_t i=0;i<words;++i)if(values[i]==UINT32_MAX){fprintf(stderr,"Audit counter saturated\n");goto done;}
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
    printf("instrumented_ms=%.6f counter_storage_mib=%.3f\n",(seconds()-start)*1000,words*4/1048576.);
    uint64_t actual[AUDIT_FIELDS]={0};
    for(unsigned r=0;r<AUDIT_REASONS;++r) {
        uint64_t sum[AUDIT_FIELDS]={0};
        for(size_t i=0;i<s->scene.cell_count;++i)for(unsigned k=0;k<AUDIT_FIELDS;++k)sum[k]+=totals[(i*AUDIT_REASONS+r)*AUDIT_FIELDS+k];
        printf("reason=%s",audit_reasons[r]);for(unsigned k=0;k<AUDIT_FIELDS;++k){printf(" field%u=%" PRIu64,k,sum[k]);actual[k]+=sum[k];}puts("");
    }
    if(memcmp(actual,expected,sizeof(actual))){fprintf(stderr,"Cell and query profiles disagree\n");goto done;}
    if(write_work(s,totals,path))goto done;
    puts("outputs=identical totals=reconciled saturation=0");rc=0;
 done:if(rc)fprintf(stderr,"Work audit failed\n");[pipeline release];[counters release];[profiles release];[points release];[output release];free(reference);return rc;
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
