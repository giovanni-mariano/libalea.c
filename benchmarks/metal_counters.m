// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
// SPDX-License-Identifier: MPL-2.0
#import <Foundation/Foundation.h>
#import <Metal/Metal.h>
#include <stdio.h>

/* Public synthetic probe: enumerate native capabilities and actually sample
 * available counter sets. Absence is reported, never replaced with zero work. */
int main(void) {@autoreleasepool {
    id<MTLDevice> device=MTLCreateSystemDefaultDevice();if(!device)return 1;
    printf("device=%s stage_boundary=%d dispatch_boundary=%d counter_sets=%lu\n",
        device.name.UTF8String,[device supportsCounterSampling:MTLCounterSamplingPointAtStageBoundary],
        [device supportsCounterSampling:MTLCounterSamplingPointAtDispatchBoundary],(unsigned long)device.counterSets.count);
    for(id<MTLCounterSet> set in device.counterSets) {
        printf("set=%s",set.name.UTF8String);
        for(id<MTLCounter> c in set.counters)printf(" counter=%s",c.name.UTF8String);
        puts("");
    }
    if(!device.counterSets.count || ![device supportsCounterSampling:MTLCounterSamplingPointAtStageBoundary]){[device release];return 0;}
    NSError* error=nil;
    NSString* source=@"#include <metal_stdlib>\nusing namespace metal; kernel void probe(device uint* out [[buffer(0)]],uint i [[thread_position_in_grid]]){uint x=i;for(uint j=0;j<1024;++j)x=x*1664525u+1013904223u;out[i]=x;}";
    id<MTLLibrary> lib=[device newLibraryWithSource:source options:nil error:&error];
    id<MTLFunction> fn=[lib newFunctionWithName:@"probe"];
    id<MTLComputePipelineState> pipeline=fn?[device newComputePipelineStateWithFunction:fn error:&error]:nil;
    id<MTLCommandQueue> queue=[device newCommandQueue];
    id<MTLBuffer> output=[device newBufferWithLength:4096*sizeof(uint32_t) options:MTLResourceStorageModeShared];
    int failed= !pipeline || !queue || !output;
    if(failed)fprintf(stderr,"Probe preparation failed: %s\n",error?error.localizedDescription.UTF8String:"unknown error");
    if(!failed)for(id<MTLCounterSet> set in device.counterSets) {
        MTLCounterSampleBufferDescriptor* descriptor=[[MTLCounterSampleBufferDescriptor alloc] init];
        descriptor.counterSet=set;descriptor.storageMode=MTLStorageModeShared;descriptor.sampleCount=2;
        id<MTLCounterSampleBuffer> samples=[device newCounterSampleBufferWithDescriptor:descriptor error:&error];[descriptor release];
        if(!samples){printf("set=%s sampling_error=%s\n",set.name.UTF8String,error?error.localizedDescription.UTF8String:"unknown error");failed=1;continue;}
        MTLComputePassDescriptor* pass=[MTLComputePassDescriptor computePassDescriptor];
        pass.sampleBufferAttachments[0].sampleBuffer=samples;
        pass.sampleBufferAttachments[0].startOfEncoderSampleIndex=0;pass.sampleBufferAttachments[0].endOfEncoderSampleIndex=1;
        id<MTLCommandBuffer> command=[queue commandBuffer];
        id<MTLComputeCommandEncoder> encoder=[command computeCommandEncoderWithDescriptor:pass];
        [encoder setComputePipelineState:pipeline];[encoder setBuffer:output offset:0 atIndex:0];
        [encoder dispatchThreads:MTLSizeMake(4096,1,1) threadsPerThreadgroup:MTLSizeMake(MIN(128,pipeline.maxTotalThreadsPerThreadgroup),1,1)];
        [encoder endEncoding];[command commit];[command waitUntilCompleted];
        NSData* data=command.status==MTLCommandBufferStatusCompleted?[samples resolveCounterRange:NSMakeRange(0,2)]:nil;
        printf("set=%s resolved_bytes=%lu gpu_seconds=%.9f\n",set.name.UTF8String,(unsigned long)data.length,command.GPUEndTime-command.GPUStartTime);
        if([set.name isEqualToString:MTLCommonCounterSetTimestamp] && data.length>=2*sizeof(MTLCounterResultTimestamp)) {
            const MTLCounterResultTimestamp* r=data.bytes;
            int valid=r[0].timestamp!=MTLCounterErrorValue && r[1].timestamp!=MTLCounterErrorValue && r[1].timestamp>r[0].timestamp;
            printf("timestamp_valid=%d timestamp_delta=%llu\n",valid,(unsigned long long)(r[1].timestamp-r[0].timestamp));
            failed|=!valid;
        }
        if(!data.length)failed=1;
        [samples release];
    }
    /* Verify execution, independently of whether any counters are exposed. */
    if(!failed)for(uint32_t i=0;i<4096;++i){uint32_t x=i;for(int j=0;j<1024;++j)x=x*1664525u+1013904223u;if(((const uint32_t*)output.contents)[i]!=x){failed=1;break;}}
    [output release];[queue release];[pipeline release];[fn release];[lib release];[device release];return failed;
}}
