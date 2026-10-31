// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
// SPDX-License-Identifier: MPL-2.0
#include <metal_stdlib>
using namespace metal;
#include "../../src/compute/compute_layout.h"
struct group_test_trace {uint count,ids[8];};
static void group_test_record(thread group_test_trace& trace,uint arg) {
    if(trace.count<8)trace.ids[trace.count]=arg;++trace.count;
}
#define COMPUTE_GROUP_PARAMETERS ,thread group_test_trace& group_work
#define COMPUTE_GROUP_ARGUMENTS ,group_work
#define COMPUTE_GROUP_LOCAL group_test_trace group_work={};
#define COMPUTE_GROUP_BUFFER
#define COMPUTE_GROUP_RECORD
#define COMPUTE_EXECUTE_SETUP
#define COMPUTE_EXECUTE_STEP
#define COMPUTE_EXECUTE_LEAF(arg) group_test_record(group_work,arg);
#define COMPUTE_NATIVE_TRACE_PARAMETERS ,thread group_test_trace& group_work
#define COMPUTE_NATIVE_TRACE_ARGUMENTS ,group_work
#define COMPUTE_NATIVE_TRACE_LEAF(arg) group_test_record(group_work,arg);
#include "../../src/compute/compute.metal"
kernel void alea_group_execute_probe(device const compute_node* code [[buffer(0)]],
    device const compute_primitive* primitives [[buffer(1)]],device const compute_cell* cells [[buffer(2)]],
    device const uint* groups [[buffer(3)]],device const uint* refs [[buffer(4)]],device uint* results [[buffer(5)]],
    uint id [[thread_position_in_grid]]) {
    if(id>=216*8)return;
    compute_point_profile cost={};group_test_trace trace={};
    float value=evaluate(code,primitives+(id%216)*3,cells[id/216],float3(0),cost,false,trace,groups,refs);
    results[10*id]=as_type<uint>(value);results[10*id+1]=trace.count;
    for(uint k=0;k<8;++k)results[10*id+2+k]=trace.ids[k];
}
