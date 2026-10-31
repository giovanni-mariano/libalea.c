// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
// SPDX-License-Identifier: MPL-2.0
#include <metal_stdlib>
using namespace metal;
#define ALEA_COMPUTE_NO_NATIVE_GROUPS 1
#include "../src/compute/compute_layout.h"
constant bool execute_groups [[function_constant(9)]];
constant bool trace_groups [[function_constant(10)]];
struct execute_trace {uint count,ids[8];};
static void execute_record(thread execute_trace& trace,uint arg) {
    if(trace_groups){if(trace.count<8)trace.ids[trace.count]=arg;++trace.count;}
}
static uint2 execute_find(device const uint* groups,uint root,uint count) {
    if(groups[2]) {
        uint mask=groups[2],at=((root*0x9e3779b9u)^(count*0x85ebca6bu));at=(at^(at>>16))&mask;
        uint base=4+4*groups[0]+3*groups[1];
        for(uint steps=0;steps<=mask;++steps,at=(at+1)&mask) {
            uint entry=groups[base+at];if(!entry)return uint2(0);
            uint p=4+4*(entry-1);
            if(groups[p]==root && groups[p+1]==count)return uint2(groups[p+2],groups[p+2]+groups[p+3]);
        }
        return uint2(0);
    }
    uint lo=0,hi=groups[0];
    while(lo<hi) {
        uint mid=lo+(hi-lo)/2,at=4+4*mid;
        if(groups[at]<root || (groups[at]==root && groups[at+1]<count))lo=mid+1;else hi=mid;
    }
    if(lo==groups[0] || groups[4+4*lo]!=root || groups[5+4*lo]!=count)return uint2(0);
    uint at=groups[6+4*lo];return uint2(at,at+groups[7+4*lo]);
}
static float primitive_value(device const compute_primitive& n,float3 p);
static float execute_list(device const uint* leaves,uint begin,uint count,
    device const compute_primitive* primitives,float3 p,thread execute_trace& trace,thread uint& visited) {
    float result=-INFINITY;
    for(uint j=0;j<count;++j) {
        uint arg=leaves[begin+j];++visited;execute_record(trace,arg);
        float v=primitive_value(primitives[arg>>1],p);
        if(!isfinite(v))return NAN;
        v=(arg&1)?-v:v;result=max(result,v);
        if(v>0)break;
    }
    return result;
}
#define COMPUTE_GROUP_PARAMETERS ,device const uint* group_ranges,device const uint* group_leaves,thread execute_trace& group_work
#define COMPUTE_GROUP_ARGUMENTS ,group_ranges,group_leaves,group_work
#define COMPUTE_GROUP_LOCAL execute_trace group_work={};
#define COMPUTE_GROUP_BUFFER device const uint* group_ranges [[buffer(17)]],device const uint* group_leaves [[buffer(18)]],
#define COMPUTE_GROUP_RECORD
#define COMPUTE_EXECUTE_SETUP uint2 execute_cursor=uint2(0);if(execute_groups && cell.count>=4)execute_cursor=execute_find(group_ranges,cell.root,cell.count);
#define COMPUTE_EXECUTE_LEAF(arg) execute_record(group_work,arg);
#define COMPUTE_EXECUTE_STEP \
    if(execute_groups && execute_cursor.x<execute_cursor.y) { \
        uint base=4+4*group_ranges[0],pc=cell.root+i; \
        while(execute_cursor.x<execute_cursor.y && group_ranges[base+3*execute_cursor.x+1]<=pc)++execute_cursor.x; \
        if(execute_cursor.x<execute_cursor.y && group_ranges[base+3*execute_cursor.x]==pc) { \
            uint at=base+3*execute_cursor.x,end=group_ranges[at+1],offset=group_ranges[at+2]; \
            uint count=(end-pc+2)/3; \
            uint visited=0; \
            float v=execute_list(group_leaves,offset,count,primitives,p,group_work,visited); \
            if(point_profile_enabled) { \
                point_add(cost.v[spatial?POINT_SPATIAL_PRIMITIVES:POINT_ORIGINAL_PRIMITIVES],visited); \
                point_add(cost.v[spatial?POINT_SPATIAL_INSTRUCTIONS:POINT_ORIGINAL_INSTRUCTIONS],visited+1); \
            } \
            if(!isfinite(v))return NAN; \
            values[sp++]=v;i=end-cell.root-1;continue; \
        } \
    }
#include "../src/compute/compute.metal"
/* Run identical bytecode with both pipeline settings, retaining the exact visited
 * signed-leaf sequence. Deliberately nonfinite primitive coefficients test errors. */
kernel void alea_group_execute_probe(device const compute_node* code [[buffer(0)]],
    device const compute_primitive* primitives [[buffer(1)]],device const compute_cell* cells [[buffer(2)]],
    device const uint* groups [[buffer(3)]],device const uint* leaves [[buffer(4)]],device uint* results [[buffer(5)]],
    uint id [[thread_position_in_grid]]) {
    if(id>=216*8)return;
    uint program=id/216,test=id%216;
    compute_cell c=cells[program];
    compute_point_profile cost={};execute_trace trace={};
    float v=evaluate(code,primitives+test*3,c,float3(0),cost,false,groups,leaves,trace);
    results[10*id]=as_type<uint>(v);results[10*id+1]=trace.count;
    for(uint j=0;j<8;++j)results[10*id+2+j]=trace.ids[j];
}
