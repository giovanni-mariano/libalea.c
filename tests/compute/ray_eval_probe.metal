// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
// SPDX-License-Identifier: MPL-2.0
#include <metal_stdlib>
using namespace metal;
#include "../../src/compute/compute_layout.h"
struct ray_group_cost {uint v[7];uint trace;};
static void ray_record(thread ray_group_cost& work,uint arg) {++work.v[2];work.trace=work.trace*31+arg+1;}
#define COMPUTE_RAY_GROUP_TRACE_LEAF(arg) ray_record(audit_cost,arg);
#define COMPUTE_RAY_AUDIT_PARAMETERS ,device const uint* audit_maps,thread ray_group_cost& audit_cost
#define COMPUTE_RAY_AUDIT_ARGUMENTS ,audit_maps,audit_cost
#define COMPUTE_RAY_AUDIT_BEGIN
#define COMPUTE_RAY_AUDIT_STEP if(n.op==COMPUTE_LEAF)ray_record(audit_cost,n.arg);
#define COMPUTE_RAY_AUDIT_BUFFER device const uint* audit_maps [[buffer(17)]],
#define COMPUTE_RAY_AUDIT_LOCAL ray_group_cost audit_cost={};
#define COMPUTE_RAY_AUDIT_STORE
#include "../../src/compute/compute.metal"
#define RAY_PROBE_BUFFERS device const uint* ray_maps [[buffer(6)]],device const uint* ray_refs [[buffer(7)]],
#define RAY_PROBE_ARGUMENTS ,ray_maps,ray_refs
#include "ray_contract_kernels.metal"
