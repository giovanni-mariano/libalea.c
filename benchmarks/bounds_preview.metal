// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
// SPDX-License-Identifier: MPL-2.0
#include <metal_stdlib>
using namespace metal;
#include "bounds_preview_layout.h"
struct bounds_vertex {float4 position [[position]];};
    constant uint edges[]={0,1,2,3,4,5,6,7,0,2,1,3,4,6,5,7,0,4,1,5,2,6,3,7};
    constant uint faces[]={0,1,3,0,3,2,4,6,7,4,7,5,0,4,5,0,5,1,
        2,3,7,2,7,6,0,2,6,0,6,4,1,5,7,1,7,3};
vertex bounds_vertex bounds_vertex_main(uint v [[vertex_id]],uint instance [[instance_id]],
        device const preview_box* boxes [[buffer(0)]],constant bounds_camera& p [[buffer(1)]]) {
    uint corner=p.filled?faces[v]:edges[v];preview_box b=boxes[instance];
    float3 world(b.b[corner&1],b.b[2+((corner>>1)&1)],b.b[4+((corner>>2)&1)]);
    float3 q=world-float3(p.eye[0],p.eye[1],p.eye[2]);
    float x=dot(q,float3(p.right[0],p.right[1],p.right[2]))/p.half_width;
    float y=dot(q,float3(p.up[0],p.up[1],p.up[2]))/p.half_height;
    float z=dot(q,float3(p.forward[0],p.forward[1],p.forward[2]));
    return {p.perspective?float4(x,y,z-BOUNDS_PREVIEW_NEAR,z):
        float4(x,y,(z-BOUNDS_PREVIEW_NEAR)/(BOUNDS_PREVIEW_FAR-BOUNDS_PREVIEW_NEAR),1)};
}
fragment float4 bounds_fragment_main() {return float4(.95f,.62f,.16f,1);}
