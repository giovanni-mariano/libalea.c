// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
// SPDX-License-Identifier: MPL-2.0
#include <metal_stdlib>
using namespace metal;
#include "../src/compute/compute_layout.h"
#include "point_work_audit_layout.h"
static void region_add(device atomic_uint* p,uint v) {
    uint old=atomic_load_explicit(p,memory_order_relaxed);
    while(!atomic_compare_exchange_weak_explicit(p,&old,old+min(v,UINT_MAX-old),memory_order_relaxed,memory_order_relaxed)) {}
}
static void audit_record(device atomic_uint* table,uint id,uint reason,bool specialized,
    compute_cell program,uint depth,uint instructions,uint primitives,float value
#ifdef ALEA_GROUP_WORK_AUDIT
    ,thread const group_trace& trace
#endif
    ) {
    uint n=atomic_load_explicit(table,memory_order_relaxed),lo=0,hi=n;
    while(lo<hi) {
        uint mid=lo+(hi-lo)/2,at=1+4*mid;
        uint root=atomic_load_explicit(table+at,memory_order_relaxed);
        uint count=atomic_load_explicit(table+at+1,memory_order_relaxed);
        if(root<program.root || (root==program.root && count<program.count))lo=mid+1;else hi=mid;
    }
    uint kind=3;
    if(lo<n) {
        uint at=1+4*lo;
        if(atomic_load_explicit(table+at,memory_order_relaxed)==program.root &&
           atomic_load_explicit(table+at+1,memory_order_relaxed)==program.count)
            kind=atomic_load_explicit(table+at+2,memory_order_relaxed);
    }
    uint r=specialized?AUDIT_SPECIALIZED:reason==POINT_NO_PROGRAM?AUDIT_NO_PROGRAM:
        reason==POINT_UNCERTIFIED?AUDIT_UNCERTIFIED:reason==POINT_OUTSIDE_DOMAIN?AUDIT_OUTSIDE:AUDIT_UNCHANGED;
    uint values[AUDIT_FIELDS]={1,instructions,primitives,value>0?instructions:0,uint(value>0)};
#ifdef ALEA_GROUP_WORK_AUDIT
    for(uint k=0;k<5;++k)values[AUDIT_GROUP_INSTRUCTIONS+k]=trace.v[k];
#endif
    uint base=1+4*n+(((id%256)*8+kind*2+uint(depth!=0))*AUDIT_REASONS+r)*AUDIT_FIELDS;
    for(uint k=0;k<AUDIT_FIELDS;++k)if(values[k])region_add(table+base+k,values[k]);
}
#define ALEA_POINT_WORK_AUDIT 1
#define ALEA_REGION_WORK_AUDIT 1
#include "../src/compute/compute.metal"
