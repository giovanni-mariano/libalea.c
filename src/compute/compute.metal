// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
// SPDX-License-Identifier: MPL-2.0

#include <metal_stdlib>
using namespace metal;
#include "compute_layout.h"

constant bool compact_programs [[function_constant(7)]];
static compute_node code_at(device const compute_node* code,uint index) {
    if(compact_programs){uint word=((device const uint*)code)[index];return {word&15u,word>>4};}
    return code[index];
}

constant bool point_profile_enabled [[function_constant(5)]];
constant bool point_candidates_enabled [[function_constant(6)]];
static void point_add(thread uint& counter,uint amount=1) {
    counter=amount>UINT_MAX-counter?UINT_MAX:counter+amount;
}

static float primitive_value(device const compute_primitive& n, float3 p) {
    device const float* a = n.p;
    switch (n.type) {
    case COMPUTE_PLANE: return a[0]*p.x + a[1]*p.y + a[2]*p.z + a[3];
    case COMPUTE_SPHERE: {
        float x = p.x-a[0], y = p.y-a[1], z = p.z-a[2];
        return x*x + y*y + z*z - a[3]*a[3];
    }
    case COMPUTE_CYL_X: {
        float y = p.y-a[0], z = p.z-a[1];
        return y*y + z*z - a[2]*a[2];
    }
    case COMPUTE_CYL_Y: {
        float x = p.x-a[0], z = p.z-a[1];
        return x*x + z*z - a[2]*a[2];
    }
    case COMPUTE_CYL_Z: {
        float x = p.x-a[0], y = p.y-a[1];
        return x*x + y*y - a[2]*a[2];
    }
    case COMPUTE_BOX:
        return max(max(a[0]-p.x, p.x-a[1]),
               max(max(a[2]-p.y, p.y-a[3]), max(a[4]-p.z, p.z-a[5])));
    case COMPUTE_QUADRIC: {
        p-=float3(a[10],a[11],a[12]);
        return p.x*(a[0]*p.x+a[3]*p.y+a[5]*p.z+a[6])+
               p.y*(a[1]*p.y+a[4]*p.z+a[7])+p.z*(a[2]*p.z+a[8])+a[9];
    }
    case COMPUTE_CONE_X: case COMPUTE_CONE_Y: case COMPUTE_CONE_Z: {
        p-=float3(a[0],a[1],a[2]);
        uint axis=n.type-COMPUTE_CONE_X;
        if(a[4]*p[axis]<0)return 1;
        float axial=p[axis];p[axis]=0;
        return dot(p,p)-a[3]*axial*axial;
    }
    case COMPUTE_TORUS_X: case COMPUTE_TORUS_Y: case COMPUTE_TORUS_Z: {
        p-=float3(a[0],a[1],a[2]);
        uint axis=n.type-COMPUTE_TORUS_X;
        float axial=p[axis]*a[4]/a[5];p[axis]=0;
        float radial=length(p)-a[3];
        return radial*radial+axial*axial-a[4]*a[4];
    }
    case COMPUTE_RCC: {
        p-=float3(a[0],a[1],a[2]);float3 axis(a[3],a[4],a[5]);
        float h=dot(p,axis);p-=h*axis;
        return max(dot(p,p)-a[7]*a[7],max(-h,h-a[6]));
    }
    }
    return NAN;
}

/* Linear short-circuit bytecode. The compiler preserves operand order unless
 * it must reorder a deep branch to stay within the checked value-stack bound. */
static float evaluate(device const compute_node* code,
                      device const compute_primitive* primitives,
                      compute_cell cell,float3 p,thread compute_point_profile& cost,bool spatial) {
    float values[COMPUTE_STACK];uint sp=0;
    for(uint i=0;i<cell.count;++i) {
        compute_node n=code_at(code,cell.root+i);
        if(point_profile_enabled)point_add(cost.v[spatial?POINT_SPATIAL_INSTRUCTIONS:POINT_ORIGINAL_INSTRUCTIONS]);
        if(n.op==COMPUTE_CONSTANT)values[sp++]=n.arg?1:-1;
        else if(n.op==COMPUTE_LEAF) {
            if(point_profile_enabled)point_add(cost.v[spatial?POINT_SPATIAL_PRIMITIVES:POINT_ORIGINAL_PRIMITIVES]);
            float v=primitive_value(primitives[n.arg>>1],p);
            if(!isfinite(v))return NAN;
            values[sp++]=(n.arg&1)?-v:v;
        } else if(n.op==COMPUTE_SKIP_UNION || n.op==COMPUTE_SKIP_INTERSECTION || n.op==COMPUTE_SKIP_REVERSE_DIFFERENCE) {
            float v=values[sp-1];
            if((n.op==COMPUTE_SKIP_UNION && v<=0) || (n.op==COMPUTE_SKIP_INTERSECTION && v>0) ||
               (n.op==COMPUTE_SKIP_REVERSE_DIFFERENCE && v<0)) {
                if(n.op==COMPUTE_SKIP_REVERSE_DIFFERENCE)values[sp-1]=-v;
                i+=n.arg;
            }
        } else if(n.op==COMPUTE_COMPLEMENT)values[sp-1]=-values[sp-1];
        else {
            float b=values[--sp],a=values[sp-1];
            values[sp-1]=n.op==COMPUTE_UNION?min(a,b):n.op==COMPUTE_INTERSECTION?max(a,b):
                n.op==COMPUTE_DIFFERENCE?max(a,-b):max(b,-a);
        }
    }
    return values[0];
}

kernel void alea_classify_points(device const compute_node* nodes [[buffer(0)]],
                                 device const compute_cell* cells [[buffer(1)]],
                                 device const float* xyz [[buffer(2)]],
                                 device compute_result* results [[buffer(3)]],
                                 constant uint& cell_count [[buffer(4)]],
                                 constant uint& point_count [[buffer(5)]],
                                 device const compute_primitive* primitives [[buffer(6)]],
                                 device compute_point_profile* profiles [[buffer(12)]],
                                 device const uint* candidates [[buffer(13)]],
                                 uint index [[thread_position_in_grid]]) {
    if (index >= point_count) return;
    size_t offset = size_t(index) * 3;
    float3 p(xyz[offset], xyz[offset+1], xyz[offset+2]);
    compute_result r = {-1, -1, 0};
    compute_point_profile cost={};
    if(point_profile_enabled){cost.v[POINT_UNIVERSES]=1;cost.v[POINT_MAX_DEPTH]=1;}
    for (uint i = point_candidates_enabled?cell_count:0;
         point_candidates_enabled?i!=COMPUTE_NONE:i<cell_count;
         i=point_candidates_enabled?candidates[i]:i+1) {
        if(point_profile_enabled) {
            point_add(cost.v[POINT_CELLS]);point_add(cost.v[POINT_ORIGINAL_CALLS]);point_add(cost.v[POINT_NO_PROGRAM]);
        }
        uint before=cost.v[POINT_ORIGINAL_INSTRUCTIONS];
        float value = evaluate(nodes, primitives, cells[i], p, cost, false);
        if(point_profile_enabled)point_add(cost.v[value>0?POINT_REJECTED_INSTRUCTIONS:POINT_ACCEPTED_INSTRUCTIONS],cost.v[POINT_ORIGINAL_INSTRUCTIONS]-before);
        if (!isfinite(value)) {
            r.owner_count = 0xffffffffu;
            results[index] = r;
            if(point_profile_enabled)profiles[index]=cost;
            return;
        }
        if(point_profile_enabled)point_add(cost.v[value>0?POINT_REJECTED:POINT_OWNERS]);
        if (value <= 0.0f) {
            ++r.owner_count;
            r.cell_id = cells[i].cell_id;
            r.material_id = cells[i].material_id;
        }
    }
    if (r.owner_count != 1) r.cell_id = r.material_id = -1;
    results[index] = r;
    if(point_profile_enabled)profiles[index]=cost;
}

#include "classify.metal"
#include "select.metal"
#include "render.metal"
