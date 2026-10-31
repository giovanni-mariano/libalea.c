// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
// SPDX-License-Identifier: MPL-2.0
/* Independent recursive fixture oracle. Roots are analytic constants on a
 * one-dimensional ray; this does not reuse the device primitive solver. */
typedef struct {int op,left,right;uint32_t arg;} ray_oracle_node;
static const ray_oracle_node ray_tree[]={
    {COMPUTE_LEAF,0,0,0},{COMPUTE_LEAF,0,0,2},{COMPUTE_LEAF,0,0,4},
    {COMPUTE_LEAF,0,0,5},{COMPUTE_INTERSECTION,0,1,0},
    {COMPUTE_INTERSECTION,4,3,0},{COMPUTE_UNION,1,2,0},
    {COMPUTE_INTERSECTION,0,6,0},{COMPUTE_DIFFERENCE,4,2,0},
    {COMPUTE_INTERSECTION,0,0,0},{COMPUTE_COMPLEMENT,5,0,0},
    {COMPUTE_LEAF,0,0,6},{COMPUTE_INTERSECTION,0,11,0}
};
static const int ray_tree_roots[]={4,5,7,8,9,10,12};
typedef struct {int sign;float next;uint32_t surface;} ray_oracle_value;
static ray_oracle_value ray_oracle(int index,float cursor,double y,uint32_t* at,uint32_t* count,uint32_t* trace) {
    ray_oracle_node n=ray_tree[index];
    if(n.op==COMPUTE_LEAF) {
        uint32_t id=n.arg>>1;float roots[2]={(id==0 || id==3)?5:id==1?3:3.5f,4.5f};if(id==2 && y!=0)roots[0]=roots[1]=4;
        float next=8;uint32_t surface=COMPUTE_NONE;
        for(int j=0;j<(id==2?2:1);++j) {
            if(roots[j]==cursor)*at=id;
            if(roots[j]>cursor && roots[j]<next){next=roots[j];surface=id;}
        }
        double x=-4+cursor+fmin((next-cursor)*.5, fmax(1,fabs(cursor)*.001));
        double value=(id==0 || id==3)?x-1:id==1?-x-1:x*x+y*y-.25;
        int sign=(value>0)-(value<0);if(n.arg&1)sign=-sign;
        ++*count;*trace=*trace*31+n.arg+1;
        return (ray_oracle_value){sign,next,surface};
    }
    ray_oracle_value a=ray_oracle(n.left,cursor,y,at,count,trace);
    if(n.op==COMPUTE_COMPLEMENT){a.sign=-a.sign;return a;}
    if((n.op==COMPUTE_UNION && a.sign<=0) || (n.op!=COMPUTE_UNION && a.sign>0))return a;
    ray_oracle_value b=ray_oracle(n.right,cursor,y,at,count,trace);
    a.sign=n.op==COMPUTE_UNION?(a.sign<b.sign?a.sign:b.sign):
        (a.sign>(n.op==COMPUTE_DIFFERENCE?-b.sign:b.sign)?a.sign:(n.op==COMPUTE_DIFFERENCE?-b.sign:b.sign));
    if(b.next<a.next){a.next=b.next;a.surface=b.surface;}return a;
}
static void ray_emit(int index,compute_node* code,uint32_t* size) {
    ray_oracle_node n=ray_tree[index];
    if(n.op==COMPUTE_LEAF){code[(*size)++]=(compute_node){COMPUTE_LEAF,n.arg};return;}
    ray_emit(n.left,code,size);
    if(n.op==COMPUTE_COMPLEMENT){code[(*size)++]=(compute_node){COMPUTE_COMPLEMENT,0};return;}
    uint32_t skip=(*size)++;
    ray_emit(n.right,code,size);code[(*size)++]=(compute_node){(uint32_t)n.op,0};
    code[skip]=(compute_node){n.op==COMPUTE_UNION?COMPUTE_SKIP_UNION:COMPUTE_SKIP_INTERSECTION,*size-skip-1};
}
static int ray_contract_check(metal_state* s,id<MTLLibrary> library) {
    compute_node code[64]={0};compute_cell cells[7]={0};uint32_t size=0;
    for(int i=0;i<7;++i){cells[i].spatial=COMPUTE_NONE;cells[i].root=size;ray_emit(ray_tree_roots[i],code,&size);cells[i].count=size-cells[i].root;}
    compute_primitive primitives[4]={{.type=COMPUTE_PLANE,.p={1,0,0,-1}},
        {.type=COMPUTE_PLANE,.p={-1,0,0,-1}},{.type=COMPUTE_SPHERE,.p={0,0,0,.5}},
        {.type=COMPUTE_PLANE,.p={1,0,0,-1}}};
    uint32_t compact[64];for(uint32_t i=0;i<size;++i)compact[i]=code[i].op|(code[i].arg<<4);
    const uint32_t empty[]={0,0,0,0,0};int result=-1;
#ifdef ALEA_RAY_CONTRACT_GROUPS
    compute_scene test_scene={.nodes=code,.primitives=primitives,.cells=cells,.node_count=size,.primitive_count=4,.cell_count=7};
    metal_state test_groups={.device=s->device};if(prepare_groups(&test_groups,&test_scene)!=1)return -1;
    const int modes=4;
#else
    const int modes=2;
#endif
    for(int mode=0;mode<modes;++mode) {
        NSError* error=nil;MTLFunctionConstantValues* constants=[[MTLFunctionConstantValues alloc] init];
        bool no=false,wide=(mode&1)!=0;uint32_t zero=0;
        [constants setConstantValue:&no type:MTLDataTypeBool atIndex:0];
        [constants setConstantValue:&zero type:MTLDataTypeUInt atIndex:2];
        [constants setConstantValue:&no type:MTLDataTypeBool atIndex:8];
        [constants setConstantValue:&wide type:MTLDataTypeBool atIndex:7];
#ifdef ALEA_RAY_CONTRACT_GROUPS
        bool grouped=mode>=2;[constants setConstantValue:&grouped type:MTLDataTypeBool atIndex:10];
#endif
        id<MTLFunction> f=[library newFunctionWithName:@"alea_ray_contract_probe" constantValues:constants error:&error];
        id<MTLComputePipelineState> pipeline=f?[s->device newComputePipelineStateWithFunction:f error:&error]:nil;
        id<MTLBuffer> buffers[]={
            [s->device newBufferWithBytes:(mode&1)?(void*)compact:(void*)code length:size*((mode&1)?4:8) options:MTLResourceStorageModeShared],
            [s->device newBufferWithBytes:primitives length:sizeof(primitives) options:MTLResourceStorageModeShared],
            [s->device newBufferWithBytes:cells length:sizeof(cells) options:MTLResourceStorageModeShared],
            [s->device newBufferWithLength:2*7*29*12*4 options:MTLResourceStorageModeShared],
            [s->device newBufferWithBytes:empty length:sizeof(empty) options:MTLResourceStorageModeShared]};
        bool ok=pipeline!=nil;for(int j=0;j<5;++j)ok=ok && buffers[j]!=nil;
        if(ok) {
            id<MTLCommandBuffer> command=[s->queue commandBuffer];id<MTLComputeCommandEncoder> e=[command computeCommandEncoder];
            ok=command && e;
            if(ok) {
                [e setComputePipelineState:pipeline];for(int j=0;j<5;++j)[e setBuffer:buffers[j] offset:0 atIndex:j];
#ifdef ALEA_RAY_CONTRACT_GROUPS
                [e setBuffer:test_groups.group_maps offset:0 atIndex:6];[e setBuffer:test_groups.group_refs offset:0 atIndex:7];
#endif
                [e dispatchThreads:MTLSizeMake(406,1,1) threadsPerThreadgroup:MTLSizeMake(32,1,1)];
                [e endEncoding];[command commit];[command waitUntilCompleted];ok=command.status==MTLCommandBufferStatusCompleted;
            }
        }
        if(ok)for(uint32_t i=0;i<406;++i) {
            uint32_t at=COMPUTE_NONE,count=0,trace=0;ray_oracle_value expected=ray_oracle(ray_tree_roots[(i%203)/29],(i%29)*.25f,i>=203?.5:0,&at,&count,&trace);
            uint32_t bits;memcpy(&bits,&expected.next,4);
            const uint32_t want[]={expected.sign<=0,bits,expected.surface,at,0,count,trace,expected.sign<=0,bits,expected.surface,at,0};
            if(memcmp((uint32_t*)buffers[3].contents+12*i,want,sizeof(want))){fprintf(stderr,"Ray contract failed format=%d case=%u\n",mode,i);ok=false;break;}
        }
        if(ok) {
            id<MTLFunction> ef=[library newFunctionWithName:@"alea_ray_error_probe" constantValues:constants error:&error];
            id<MTLComputePipelineState> ep=ef?[s->device newComputePipelineStateWithFunction:ef error:&error]:nil;
            if(!ep)ok=false;
            for(int bad=0;ok && bad<2;++bad) {
                compute_primitive changed[4];memcpy(changed,primitives,sizeof(changed));changed[1].p[3]=bad?INFINITY:NAN;
                id<MTLBuffer> invalid=[s->device newBufferWithBytes:changed length:sizeof(changed) options:MTLResourceStorageModeShared];
                id<MTLCommandBuffer> command=[s->queue commandBuffer];id<MTLComputeCommandEncoder> e=[command computeCommandEncoder];
                if(!invalid || !command || !e)ok=false;
                if(ok) {
                    [e setComputePipelineState:ep];for(int j=0;j<5;++j)[e setBuffer:j==1?invalid:buffers[j] offset:0 atIndex:j];
#ifdef ALEA_RAY_CONTRACT_GROUPS
                    [e setBuffer:test_groups.group_maps offset:0 atIndex:6];[e setBuffer:test_groups.group_refs offset:0 atIndex:7];
#endif
                    [e dispatchThreads:MTLSizeMake(2,1,1) threadsPerThreadgroup:MTLSizeMake(2,1,1)];[e endEncoding];[command commit];[command waitUntilCompleted];
                    const uint32_t want[]={1,2,34,0,1,1};
                    ok=command.status==MTLCommandBufferStatusCompleted && !memcmp(buffers[3].contents,want,sizeof(want));
                }
                [invalid release];
            }
            [ep release];[ef release];
            if(!ok)fprintf(stderr,"Nonfinite ray visitation failed format/mode=%d\n",mode);
        }
        if(ok) {
            id<MTLFunction> inf=[library newFunctionWithName:@"alea_ray_interval_probe" constantValues:constants error:&error];
            id<MTLComputePipelineState> inp=inf?[s->device newComputePipelineStateWithFunction:inf error:&error]:nil;
            id<MTLCommandBuffer> command=[s->queue commandBuffer];id<MTLComputeCommandEncoder> e=[command computeCommandEncoder];
            ok=inp && command && e;
            if(ok) {
                [e setComputePipelineState:inp];for(int j=0;j<5;++j)[e setBuffer:buffers[j] offset:0 atIndex:j];
#ifdef ALEA_RAY_CONTRACT_GROUPS
                [e setBuffer:test_groups.group_maps offset:0 atIndex:6];[e setBuffer:test_groups.group_refs offset:0 atIndex:7];
#endif
                [e dispatchThreads:MTLSizeMake(8,1,1) threadsPerThreadgroup:MTLSizeMake(8,1,1)];[e endEncoding];[command commit];[command waitUntilCompleted];
                ok=command.status==MTLCommandBufferStatusCompleted;
                const float entries[]={3,4.5f,0,3.25f},exits[]={3.5f,5,0,3.5f};
                const uint32_t surfaces[]={1,2,COMPUTE_NONE,COMPUTE_NONE};
                for(int i=0;ok && i<8;++i) {
                    uint32_t en,ex;memcpy(&en,&entries[i%4],4);memcpy(&ex,&exits[i%4],4);
                    const uint32_t want[]={i%4!=2,0,en,ex,surfaces[i%4]};
                    ok=!memcmp((uint32_t*)buffers[3].contents+5*i,want,sizeof(want));
                }
            }
            [inp release];[inf release];if(!ok)fprintf(stderr,"Analytic material intervals failed mode=%d\n",mode);
        }
        for(int j=0;j<5;++j)[buffers[j] release];[pipeline release];[f release];[constants release];
        if(!ok) {
#ifdef ALEA_RAY_CONTRACT_GROUPS
            [test_groups.group_maps release];[test_groups.group_refs release];
#endif
            return -1;
        }
    }
#ifdef ALEA_RAY_CONTRACT_GROUPS
    [test_groups.group_maps release];[test_groups.group_refs release];
    puts("ray_group_contract_modes=generic,grouped");
#endif
    puts("ray_material_intervals=passed nonfinite_visitation=passed");
    result=0;puts("ray_contract_cases=406 formats=wide,compact ordinary=resumed analytic_oracle=passed");return result;
}
