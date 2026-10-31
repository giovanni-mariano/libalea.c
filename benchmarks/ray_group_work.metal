// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
// SPDX-License-Identifier: MPL-2.0
#include <metal_stdlib>
using namespace metal;
#include "../src/compute/compute_layout.h"
constant bool ray_audit_option [[function_constant(11)]];
constant bool ray_audit_active=is_function_constant_defined(ray_audit_option)?ray_audit_option:true;
struct ray_group_cost {uint v[7];uint trace;};
static void ray_add(thread uint& x,uint amount=1) {x+=min(amount,UINT_MAX-x);}
#define COMPUTE_RAY_AUDIT_PARAMETERS ,device const uint* audit_maps,thread ray_group_cost& audit_cost
#define COMPUTE_RAY_AUDIT_ARGUMENTS ,audit_maps,audit_cost
#define COMPUTE_RAY_AUDIT_BEGIN \
    uint2 audit_cursor=uint2(0);if(ray_audit_active){audit_cursor=native_group_find(audit_maps,cell.root,cell.count);ray_add(audit_cost.v[0]);}
#define COMPUTE_RAY_AUDIT_STEP \
    if(ray_audit_active) { ray_add(audit_cost.v[1]); \
    if(n.op==COMPUTE_LEAF){ray_add(audit_cost.v[2]);audit_cost.trace=audit_cost.trace*31+n.arg+1;} \
    while(audit_cursor.x<audit_cursor.y && audit_maps[5+4*audit_maps[0]+3*audit_cursor.x]<=cell.root+i)++audit_cursor.x; \
    if(audit_cursor.x<audit_cursor.y && audit_maps[4+4*audit_maps[0]+3*audit_cursor.x]<=cell.root+i) { \
        ray_add(audit_cost.v[3]);if(n.op==COMPUTE_LEAF)ray_add(audit_cost.v[4]); \
        if(audit_maps[4+4*audit_maps[0]+3*audit_cursor.x]==cell.root+i)ray_add(audit_cost.v[5]); \
    } \
    }
#define COMPUTE_RAY_AUDIT_BUFFER device const uint* audit_maps [[buffer(17)]],device ray_group_cost* audit_rows [[buffer(18)]],
#define COMPUTE_RAY_AUDIT_LOCAL ray_group_cost audit_cost={};
#define COMPUTE_RAY_AUDIT_STORE if(ray_audit_active)audit_rows[index]=audit_cost;
#define COMPUTE_RENDER_NO_RESUME 1
#include "../src/compute/compute.metal"
/* Test ordinary expression state against one-opcode resumptions. */
#include "../src/compute/resume.metal"
#ifndef RAY_PROBE_BUFFERS
#define RAY_PROBE_BUFFERS
#define RAY_PROBE_ARGUMENTS ,nullptr,nullptr
#endif
#include "../tests/compute/ray_contract_kernels.metal"
