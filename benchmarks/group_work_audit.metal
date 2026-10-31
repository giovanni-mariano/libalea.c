// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
// SPDX-License-Identifier: MPL-2.0
#include <metal_stdlib>
using namespace metal;
#include "../src/compute/compute_layout.h"
struct group_trace {uint v[5];};
/* Each exact program range owns its sorted disjoint groups. Partial aliases
 * must not borrow a larger program's group classification. */
static uint2 group_find(device const uint* groups,uint root,uint count) {
    uint lo=0,hi=groups[0];
    while(lo<hi) {
        uint mid=lo+(hi-lo)/2,at=1+4*mid;
        if(groups[at]<root || (groups[at]==root && groups[at+1]<count))lo=mid+1;else hi=mid;
    }
    if(lo==groups[0] || groups[1+4*lo]!=root || groups[2+4*lo]!=count)return uint2(0);
    uint at=groups[3+4*lo];return uint2(at,at+groups[4+4*lo]);
}
static void group_visit(device const uint* groups,thread uint2& cursor,uint pc,uint op,thread group_trace& trace) {
    uint base=1+4*groups[0];
    while(cursor.x<cursor.y && groups[base+1+2*cursor.x]<=pc)++cursor.x;
    if(cursor.x==cursor.y || pc<groups[base+2*cursor.x])return;
    ++trace.v[0];trace.v[1]+=op==COMPUTE_LEAF;
    trace.v[2]+=pc==groups[base+2*cursor.x];trace.v[3]+=op==COMPUTE_SKIP_INTERSECTION;
    trace.v[4]+=op==COMPUTE_INTERSECTION;
}
#define ALEA_GROUP_WORK_AUDIT 1
#define COMPUTE_GROUP_PARAMETERS ,device const uint* group_ranges,thread group_trace& group_work
#define COMPUTE_GROUP_ARGUMENTS ,group_ranges,group_work
#define COMPUTE_GROUP_LOCAL group_trace group_work={};
#define COMPUTE_GROUP_BUFFER device const uint* group_ranges [[buffer(21)]],
#define COMPUTE_GROUP_RECORD ,group_work
#include "region_work_audit.metal"
/* Public diagnostic regression: nested and partial code aliases have distinct
 * program-local maps, including a single-leaf view inside a larger AND. */
kernel void alea_group_map_probe(device const uint* map [[buffer(0)]],
                                device uint* result [[buffer(1)]],uint id [[thread_position_in_grid]]) {
    if(id>=6)return;
    const uint roots[]={0,0,2,2,0,0},counts[]={1,7,4,4,7,2},pcs[]={0,0,2,5,6,0};
    const uint ops[]={COMPUTE_LEAF,COMPUTE_LEAF,COMPUTE_LEAF,COMPUTE_INTERSECTION,COMPUTE_INTERSECTION,COMPUTE_LEAF};
    uint2 cursor=group_find(map,roots[id],counts[id]);group_trace trace={};
    group_visit(map,cursor,pcs[id],ops[id],trace);
    for(uint k=0;k<5;++k)result[5*id+k]=trace.v[k];
}
