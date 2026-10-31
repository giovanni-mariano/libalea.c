// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
// SPDX-License-Identifier: MPL-2.0
static int ray_double_compare(const void* a,const void* b) {double x=*(const double*)a,y=*(const double*)b;return (x>y)-(x<y);}
static int audit_render(metal_state* s,const render_config_t* cfg,const render_camera_t* camera,
                        render_framebuffer_t* fb,const char* path) {
    const char* env=getenv("ALEA_BENCH_RAY_GROUPS");
    if(env && strcmp(env,"0") && strcmp(env,"1"))return -1;
    bool enabled=env && !strcmp(env,"1");int warm=5;
    const char* warm_env=getenv("ALEA_BENCH_RAY_WARM");
    if(warm_env){if(strlen(warm_env)!=1 || warm_env[0]<'1' || warm_env[0]>'5')return -1;warm=warm_env[0]-'0';}
    if(!s->groups_enabled || !s->group_maps || !s->group_refs)return -1;
    const size_t n=(size_t)fb->width*fb->height;int rc=-1;compute_render_pixel* reference=NULL;
    NSError* error=nil;id<MTLLibrary> library=nil;id<MTLComputePipelineState> pipeline=nil;
    double pipeline_start=ray_seconds();
    dispatch_data_t data=dispatch_data_create(alea_ray_group_work_metallib,alea_ray_group_work_metallib_len,NULL,DISPATCH_DATA_DESTRUCTOR_DEFAULT);
    library=[s->device newLibraryWithData:data error:&error];dispatch_release(data);
    if(!library)goto done;
    /* Run the independent probes once per public validation process. */
    if(getenv("ALEA_BENCH_RAY_VALIDATE") && ray_contract_check(s,library))goto done;
    if(alea_compute_metal_render(s,cfg,camera,fb,NULL,NULL))goto done;
    reference=malloc(n*sizeof(*reference));if(!reference)goto done;memcpy(reference,s->pixels.contents,n*sizeof(*reference));
    MTLFunctionConstantValues* constants=[[MTLFunctionConstantValues alloc] init];bool profile=s->profile_prefix!=NULL,no=false;
    [constants setConstantValue:&profile type:MTLDataTypeBool atIndex:0];
    [constants setConstantValue:&s->lazy_threshold type:MTLDataTypeUInt atIndex:2];
    [constants setConstantValue:&s->refill type:MTLDataTypeBool atIndex:4];
    [constants setConstantValue:&s->compact_code type:MTLDataTypeBool atIndex:7];
    [constants setConstantValue:&no type:MTLDataTypeBool atIndex:8];
    [constants setConstantValue:&enabled type:MTLDataTypeBool atIndex:10];
    [constants setConstantValue:&no type:MTLDataTypeBool atIndex:11];
    id<MTLFunction> f=[library newFunctionWithName:@"alea_render_image" constantValues:constants error:&error];
    if(f)pipeline=[s->device newComputePipelineStateWithFunction:f error:&error];[f release];[constants release];
    if(!pipeline)goto done;
    printf("experiment_setup_seconds=%.6f validation=%d\n",ray_seconds()-pipeline_start,getenv("ALEA_BENCH_RAY_VALIDATE")!=NULL);
    ray_audit_maps=[s->group_maps retain];ray_audit_rows=[s->device newBufferWithLength:n*8*4 options:MTLResourceStorageModeShared];
    if(!ray_audit_rows || s->group_maps.length+s->group_refs.length+ray_audit_rows.length>64*1024*1024 ||
       s->device.currentAllocatedSize>1024ull*1024*1024)goto done;
    [s->render_pipeline release];s->render_pipeline=pipeline;pipeline=nil;ray_audit_enabled=1;s->ray_groups_enabled=enabled;
    double times[5],first=0;
    for(int frame=-1;frame<warm;++frame) {
        double start=ray_seconds();if(alea_compute_metal_render(s,cfg,camera,fb,NULL,NULL))goto done;
        double elapsed=ray_seconds()-start;
        if(memcmp(reference,s->pixels.contents,n*sizeof(*reference))){fprintf(stderr,"Grouped/generic ray outputs differ\n");goto done;}
        if(getenv("ALEA_BENCH_RAY_VALIDATE") && (getenv("ALEA_BENCH_RAY_LATTICE")?
           ray_lattice_reference(cfg,camera,fb):check_render(cfg,camera,fb))){fprintf(stderr,"Analytic render reference failed\n");goto done;}
        if(frame<0)first=elapsed;else times[frame]=elapsed;
    }
    qsort(times,warm,sizeof(double),ray_double_compare);double median=(times[warm/2]+times[(warm-1)/2])*.5;
    FILE* out=fopen(path,"w");if(!out)goto done;
    fprintf(out,"groups,pixels,aa,profiled,first_seconds,median_seconds,min_seconds,max_seconds,warm_samples\n");
    fprintf(out,"%d,%zu,%d,%d,%.9f,%.9f,%.9f,%.9f,%d\n",enabled,n,cfg->aa_samples,profile,first,median,times[0],times[warm-1],warm);
    int failed=ferror(out);failed|=fclose(out)!=0;if(failed)goto done;
    printf("groups=%d pixels=%zu first_seconds=%.6f warm_median_seconds=%.6f outputs=identical instruction_unit=%s\n",enabled,n,first,median,enabled?"bytecode_and_group_ops":"bytecode");
    if(profile) {
        uint64_t primitives=0;const compute_profile* rows=s->ray_costs.contents;
        for(size_t i=0;i<n;++i){if(rows[i].primitives==UINT32_MAX || rows[i].instructions==UINT32_MAX)goto done;primitives+=rows[i].primitives;}
        printf("profile_primitives=%llu saturation=0\n",(unsigned long long)primitives);
    }
    rc=0;
 done:
    if(rc)fprintf(stderr,"Ray execution experiment failed: %s %s\n",alea_error(),error?error.localizedDescription.UTF8String:"");
    ray_audit_enabled=0;[ray_audit_maps release];[ray_audit_rows release];ray_audit_maps=ray_audit_rows=nil;
    [pipeline release];[library release];free(reference);return rc;
}
