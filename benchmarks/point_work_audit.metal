// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
// SPDX-License-Identifier: MPL-2.0
#include <metal_stdlib>
using namespace metal;
#include "../src/compute/compute_layout.h"
#include "point_work_audit_layout.h"
static void audit_add(device atomic_uint* p,uint value) {
    uint old=atomic_load_explicit(p,memory_order_relaxed);
    while(!atomic_compare_exchange_weak_explicit(p,&old,old+min(value,UINT_MAX-old),
        memory_order_relaxed,memory_order_relaxed)) {}
}
static void audit_record(device atomic_uint* table,uint id,uint reason,bool specialized,compute_cell program,uint depth,
        uint instructions,uint primitives,float value) {
    uint r=specialized?AUDIT_SPECIALIZED:reason==POINT_NO_PROGRAM?AUDIT_NO_PROGRAM:
        reason==POINT_UNCERTIFIED?AUDIT_UNCERTIFIED:reason==POINT_OUTSIDE_DOMAIN?AUDIT_OUTSIDE:AUDIT_UNCHANGED;
    uint values[AUDIT_FIELDS]={1,instructions,primitives,value>0?instructions:0,uint(value>0)};
    for(uint k=0;k<AUDIT_FIELDS;++k)if(values[k])
        audit_add(table+(size_t(id)*AUDIT_REASONS+r)*AUDIT_FIELDS+k,values[k]);
}
#define ALEA_POINT_WORK_AUDIT 1
#include "../src/compute/compute.metal"
