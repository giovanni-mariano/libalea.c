// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
// SPDX-License-Identifier: MPL-2.0
kernel void alea_ray_contract_probe(device const compute_node* code [[buffer(0)]],
    device const compute_primitive* primitives [[buffer(1)]],device const compute_cell* cells [[buffer(2)]],
    device uint* output [[buffer(3)]],device const uint* maps [[buffer(4)]], RAY_PROBE_BUFFERS uint id [[thread_position_in_grid]]) {
    if(id>=2*7*29)return;
    uint model=(id%(7*29))/29;float y=id>=7*29?.5f:0;float cursor=float(id%29)*.25f;
    float next;uint surface,at=COMPUTE_NONE;bool error=false;compute_profile cost={};ray_group_cost work={};
    bool inside=ray_expression(code,primitives,cells[model],float3(-4,y,0),float3(1,0,0),cursor,8,next,surface,at,error,cost,nullptr,maps,work RAY_PROBE_ARGUMENTS);
    uint base=id*12;output[base]=inside;output[base+1]=as_type<uint>(next);output[base+2]=surface;
    output[base+3]=at;output[base+4]=error;output[base+5]=work.v[2];output[base+6]=work.trace;
    resume_expression e={};e.program=cells[model];e.end=8; e.at_cursor=COMPUTE_NONE;bool failed=false;
    while(e.pc<e.program.count)if(resume_expression_step(code,primitives,e,float3(-4,y,0),float3(1,0,0),cursor,failed,cost)<0)break;
    output[base+7]=e.values[0]<=0;output[base+8]=as_type<uint>(e.until[0]);
    output[base+9]=e.boundaries[0];output[base+10]=e.at_cursor;output[base+11]=failed;
}

/* Nonfinite plane coefficients must fail only if the leaf is visited. */
kernel void alea_ray_error_probe(device const compute_node* code [[buffer(0)]],
    device const compute_primitive* primitives [[buffer(1)]],device const compute_cell* cells [[buffer(2)]],
    device uint* output [[buffer(3)]],device const uint* maps [[buffer(4)]], RAY_PROBE_BUFFERS uint id [[thread_position_in_grid]]) {
    if(id>=2)return;float cursor=id?6:0,next;uint surface,at=COMPUTE_NONE;bool error=false;
    compute_profile cost={};ray_group_cost work={};
    ray_expression(code,primitives,cells[0],float3(-4,0,0),float3(1,0,0),cursor,8,next,surface,at,error,cost,nullptr,maps,work RAY_PROBE_ARGUMENTS);
    output[3*id]=error;output[3*id+1]=work.v[2];output[3*id+2]=work.trace;
}
/* A slab minus a sphere has two occupied intervals. Check interval walking as
 * well as expression tuples, with starts outside, inside and at a boundary. */
kernel void alea_ray_interval_probe(device const compute_node* code [[buffer(0)]],
    device const compute_primitive* primitives [[buffer(1)]],device const compute_cell* cells [[buffer(2)]],
    device uint* output [[buffer(3)]],device const uint* maps [[buffer(4)]], RAY_PROBE_BUFFERS uint id [[thread_position_in_grid]]) {
    if(id>=8)return;const float starts[]={0,3.5f,5,3.25f};
    float enter=0,exit=0;uint surface=COMPUTE_NONE;bool error=false;compute_profile cost={};ray_group_cost work={};
    bool hit=cell_interval_impl(code,primitives,cells[id<4?1:3],float3(-4,0,0),float3(1,0,0),starts[id%4],8,true,
        enter,exit,surface,error,cost,nullptr,nullptr,maps,work RAY_PROBE_ARGUMENTS);
    uint at=id*5;output[at]=hit;output[at+1]=error;output[at+2]=hit?as_type<uint>(enter):0;
    output[at+3]=hit?as_type<uint>(exit):0;output[at+4]=hit?surface:COMPUTE_NONE;
}
