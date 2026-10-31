// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
// SPDX-License-Identifier: MPL-2.0
#include <metal_stdlib>
using namespace metal;
#include "reprojection_layout.h"

kernel void preview_clear(device atomic_uint* depths [[buffer(0)]],
        device atomic_uint* owners [[buffer(1)]],uint i [[thread_position_in_grid]]) {
    atomic_store_explicit(depths+i,as_type<uint>(INFINITY),memory_order_relaxed);
    atomic_store_explicit(owners+i,PREVIEW_NONE,memory_order_relaxed);
}
kernel void preview_project(device const preview_point* points [[buffer(0)]],
        constant preview_params& p [[buffer(1)]],device preview_projection* projected [[buffer(2)]],
        device atomic_uint* depths [[buffer(3)]],uint i [[thread_position_in_grid]]) {
    projected[i]={PREVIEW_NONE,0};
    if(!points[i].valid)return;
    float3 world(points[i].world[0],points[i].world[1],points[i].world[2]);
    float3 q=world-float3(p.eye[0],p.eye[1],p.eye[2]);
    float z=dot(q,float3(p.forward[0],p.forward[1],p.forward[2]));
    if(!all(isfinite(world)) || !isfinite(z) || z<=0)return;
    float u=dot(q,float3(p.right[0],p.right[1],p.right[2]))/(p.half_width*(p.perspective?z:1));
    float v=dot(q,float3(p.up[0],p.up[1],p.up[2]))/(p.half_height*(p.perspective?z:1));
    float x=(.5f+.5f*u)*p.width,y=(.5f-.5f*v)*p.height;
    if(!isfinite(x) || !isfinite(y) || x<0 || y<0 || x>=p.width || y>=p.height)return;
    uint pixel=uint(y)*p.width+uint(x),depth=as_type<uint>(z);
    projected[i]={pixel,depth};
    atomic_fetch_min_explicit(depths+pixel,depth,memory_order_relaxed);
}
/* Separate passes avoid racing a nearest-depth update against its RGB write.
 * Positive finite float bits preserve depth order; equal depths use source order. */
kernel void preview_choose(device const preview_projection* projected [[buffer(0)]],
        device const atomic_uint* depths [[buffer(1)]],device atomic_uint* owners [[buffer(2)]],
        uint i [[thread_position_in_grid]]) {
    preview_projection q=projected[i];if(q.pixel==PREVIEW_NONE)return;
    if(q.depth==atomic_load_explicit(depths+q.pixel,memory_order_relaxed))
        atomic_fetch_min_explicit(owners+q.pixel,i,memory_order_relaxed);
}
kernel void preview_resolve(device const preview_point* points [[buffer(0)]],
        device const uint* owners [[buffer(1)]],device float4* color [[buffer(2)]],
        constant preview_params& p [[buffer(3)]],uint i [[thread_position_in_grid]]) {
    uint owner=owners[i];
    if(owner==PREVIEW_NONE) {
        float shade=(((i%p.width)/8+(i/p.width)/8)&1)?0.12f:0.22f;
        color[i]=float4(shade,shade,shade,0); // Alpha zero means unknown, including background.
    } else color[i]=float4(points[owner].color[0],points[owner].color[1],points[owner].color[2],1);
}
