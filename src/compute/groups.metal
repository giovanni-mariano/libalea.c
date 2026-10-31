// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
// SPDX-License-Identifier: MPL-2.0
#ifndef ALEA_COMPUTE_NO_NATIVE_GROUPS
constant bool native_groups_option [[function_constant(9)]];
constant bool native_groups = is_function_constant_defined(native_groups_option)?native_groups_option:false;
static float primitive_value(device const compute_primitive& n,float3 p);
static uint2 native_group_find(device const uint* map,uint root,uint count) {
    uint mask=map[2],h=(root*0x9e3779b9u)^(count*0x85ebca6bu);h=(h^(h>>16))&mask;
    uint base=4+4*map[0]+3*map[1];
    for(uint steps=0;steps<=mask;++steps,h=(h+1)&mask) {
        uint entry=map[base+h];if(!entry)return uint2(0);
        uint at=4+4*(entry-1);
        if(map[at]==root && map[at+1]==count)return uint2(map[at+2],map[at+2]+map[at+3]);
    }
    return uint2(0);
}
#ifndef COMPUTE_NATIVE_TRACE_PARAMETERS
#define COMPUTE_NATIVE_TRACE_PARAMETERS
#define COMPUTE_NATIVE_TRACE_ARGUMENTS
#define COMPUTE_NATIVE_TRACE_LEAF(arg)
#endif
static float native_group_value(device const uint* refs,uint offset,uint count,
    device const compute_primitive* primitives,float3 p,thread uint& visited COMPUTE_NATIVE_TRACE_PARAMETERS) {
    float result=-INFINITY;
    for(uint j=0;j<count;++j) {
        uint arg=refs[offset+j];++visited;COMPUTE_NATIVE_TRACE_LEAF(arg)
        float v=primitive_value(primitives[arg>>1],p);if(!isfinite(v))return NAN;
        v=(arg&1)?-v:v;result=max(result,v);if(v>0)break;
    }
    return result;
}
#define COMPUTE_NATIVE_PARAMETERS ,device const uint* native_maps,device const uint* native_refs
#define COMPUTE_NATIVE_ARGUMENTS ,native_maps,native_refs
#define COMPUTE_NATIVE_BUFFER device const uint* native_maps [[buffer(17)]],device const uint* native_refs [[buffer(18)]],
#define COMPUTE_NATIVE_SETUP uint2 native_cursor=uint2(0);if(native_groups && cell.count>=4)native_cursor=native_group_find(native_maps,cell.root,cell.count);
#define COMPUTE_NATIVE_STEP \
    if(native_groups && native_cursor.x<native_cursor.y) { \
        uint base=4+4*native_maps[0],pc=cell.root+i; \
        while(native_cursor.x<native_cursor.y && native_maps[base+3*native_cursor.x+1]<=pc)++native_cursor.x; \
        if(native_cursor.x<native_cursor.y && native_maps[base+3*native_cursor.x]==pc) { \
            uint at=base+3*native_cursor.x,end=native_maps[at+1],visited=0; \
            float v=native_group_value(native_refs,native_maps[at+2],(end-pc+2)/3,primitives,p,visited COMPUTE_NATIVE_TRACE_ARGUMENTS); \
            if(point_profile_enabled) { \
                point_add(cost.v[spatial?POINT_SPATIAL_PRIMITIVES:POINT_ORIGINAL_PRIMITIVES],visited); \
                point_add(cost.v[spatial?POINT_SPATIAL_INSTRUCTIONS:POINT_ORIGINAL_INSTRUCTIONS],visited+1); \
            } \
            if(!isfinite(v))return NAN;values[sp++]=v;i=end-cell.root-1;continue; \
        } \
    }
#else
#define COMPUTE_NATIVE_PARAMETERS
#define COMPUTE_NATIVE_ARGUMENTS
#define COMPUTE_NATIVE_BUFFER
#define COMPUTE_NATIVE_SETUP
#define COMPUTE_NATIVE_STEP
#endif
