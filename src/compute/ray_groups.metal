// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
// SPDX-License-Identifier: MPL-2.0
/* Ordinary-ray evaluation over the context's shared intersection-group data. */
#ifndef ALEA_COMPUTE_NO_NATIVE_GROUPS
#ifndef COMPUTE_RAY_GROUP_TRACE_LEAF
#define COMPUTE_RAY_GROUP_TRACE_LEAF(arg)
#endif
constant bool ray_groups_option [[function_constant(10)]];
constant bool ray_groups=is_function_constant_defined(ray_groups_option)?ray_groups_option:false;
#define COMPUTE_RAY_EXECUTE_PARAMETERS ,device const uint* ray_maps,device const uint* ray_refs
#define COMPUTE_RAY_EXECUTE_ARGUMENTS ,ray_maps,ray_refs
#define COMPUTE_RAY_EXECUTE_BUFFER device const uint* ray_maps [[buffer(19)]],device const uint* ray_refs [[buffer(20)]],
#define COMPUTE_RAY_EXECUTE_SETUP uint2 ray_cursor=uint2(0);if(ray_groups && cell.count>=4)ray_cursor=native_group_find(ray_maps,cell.root,cell.count);
/* The loop accumulates the earliest crossing among all visited leaves, even
 * when the last leaf rejects. Strict comparisons preserve the left boundary
 * at ties; at_cursor still follows the last matching visited root. */
#define COMPUTE_RAY_EXECUTE_STEP \
    if(ray_groups && ray_cursor.x<ray_cursor.y) { \
        uint base=4+4*ray_maps[0],pc=cell.root+i; \
        while(ray_cursor.x<ray_cursor.y && ray_maps[base+3*ray_cursor.x+1]<=pc)++ray_cursor.x; \
        if(ray_cursor.x<ray_cursor.y && ray_maps[base+3*ray_cursor.x]==pc) { \
            uint at=base+3*ray_cursor.x,end_pc=ray_maps[at+1],offset=ray_maps[at+2],length=(end_pc-pc+2)/3; \
            int value=-1;float valid_until=end;uint boundary=COMPUTE_NONE; \
            if(profile_enabled)cost.instructions=profile_sum(cost.instructions,1); \
            for(uint j=0;j<length;++j) { \
                uint arg=ray_refs[offset+j],id=arg>>1;float4 roots; \
                COMPUTE_RAY_GROUP_TRACE_LEAF(arg) \
                if(profile_enabled){cost.instructions=profile_sum(cost.instructions,1);cost.primitives=profile_sum(cost.primitives,1);} \
                uint count=primitive_roots(primitives[id],o,d,roots,error); \
                if(root_reuse_enabled) { \
                    uint key[7]={id,as_type<uint>(o.x),as_type<uint>(o.y),as_type<uint>(o.z),as_type<uint>(d.x),as_type<uint>(d.y),as_type<uint>(d.z)}; \
                    uint rejected=uint(error);if(!error)for(uint k=0;k<count;++k)if(!isfinite(roots[k]))rejected=2; \
                    root_reuse_record(reuse,key,primitives[id].type,rejected); \
                } \
                if(error)return false; \
                float crossing=end;uint leaf_boundary=COMPUTE_NONE; \
                for(uint k=0;k<count;++k) { \
                    if(roots[k]==cursor)at_cursor=id; \
                    if(roots[k]>cursor && roots[k]<crossing){crossing=roots[k];leaf_boundary=id;} \
                } \
                float probe=cursor+(crossing-cursor)*.5f; \
                probe=min(probe,cursor+max(1.0f,abs(cursor)*.001f));if(probe==cursor)probe=crossing; \
                float v=primitive_value(primitives[id],o+probe*d);if(!isfinite(v)){error=true;return false;} \
                int sign=(v>0)-(v<0);if(arg&1)sign=-sign;value=max(value,sign); \
                if(crossing<valid_until){valid_until=crossing;boundary=leaf_boundary;} \
                if(sign>0)break; \
            } \
            values[sp]=value;until[sp]=valid_until;boundaries[sp]=boundary;++sp;i=end_pc-cell.root-1;continue; \
        } \
    }
#else
#define COMPUTE_RAY_EXECUTE_PARAMETERS
#define COMPUTE_RAY_EXECUTE_ARGUMENTS
#define COMPUTE_RAY_EXECUTE_BUFFER
#define COMPUTE_RAY_EXECUTE_SETUP
#define COMPUTE_RAY_EXECUTE_STEP
#endif
