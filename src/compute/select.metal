// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
// SPDX-License-Identifier: MPL-2.0

static void select_add(device atomic_uint* p,uint value) {
    uint old=atomic_load_explicit(p,memory_order_relaxed);
    while(!atomic_compare_exchange_weak_explicit(p,&old,old+min(value,UINT_MAX-old),
        memory_order_relaxed,memory_order_relaxed)) {}
}
static void select_record(device atomic_uint* cells,device atomic_uint* depths,uint id,uint depth,
    thread const compute_point_profile& cost,bool rejected,bool spatial,bool lattice) {
    uint values[SELECT_COUNTER_COUNT]={1,cost.v[POINT_ORIGINAL_INSTRUCTIONS],
        cost.v[POINT_ORIGINAL_PRIMITIVES],uint(rejected),uint(spatial),uint(lattice)};
    for(uint k=0;k<SELECT_COUNTER_COUNT;++k)if(values[k]) {
        select_add(cells+size_t(id)*SELECT_COUNTER_COUNT+k,values[k]);
        select_add(depths+depth*SELECT_COUNTER_COUNT+k,values[k]);
    }
}

/* Ordered selection follows a single winning path. The last containing fill
 * remains the answer if no child matches; no sibling backtracking is allowed. */
kernel void alea_select_hierarchy(
    device const compute_node* nodes [[buffer(0)]],
    device const compute_cell* cells [[buffer(1)]],
    device const float* xyz [[buffer(2)]],
    device compute_selection* results [[buffer(3)]],
    constant uint& root [[buffer(4)]],constant uint& count [[buffer(5)]],
    device const compute_primitive* primitives [[buffer(6)]],
    device const compute_transform* transforms [[buffer(7)]],
    device const metal_point_cell* links [[buffer(8)]],
    device const compute_lattice* lattices [[buffer(9)]],
    device const uint* fills [[buffer(10)]],
    device const compute_spatial_program* spatial [[buffer(11)]],
    device const uint* candidates [[buffer(13)]],
    device atomic_uint* costs [[buffer(14)]],device atomic_uint* depths [[buffer(15)]],uint index [[thread_position_in_grid]]) {
    if(index>=count)return;
    float3 p(xyz[3*index],xyz[3*index+1],xyz[3*index+2]);
    uint id=root,depth=0;compute_selection result={-1,-1,-1,0};
    compute_point_profile unused={};
    while(id!=COMPUTE_NONE) {
        compute_cell cell=cells[id];
        uint next=point_candidates_enabled && depth==0?candidates[id]:links[id].next;
        if(cell.lattice!=COMPUTE_NONE) {
            unused={};
            bool error=false,located;
            float3 local=p;
            uint child=point_lattice(lattices[cell.lattice],fills,local,error,located);
            if(point_profile_enabled)select_record(costs,depths,id,depth,unused,false,false,true);
            if(error){result.status=COMPUTE_NONE;break;}
            if(!located){id=next;continue;}
            p=local;
            result={cell.cell_id,cell.material_id,int(depth),2};
            id=child;
        } else {
            uint selection=0;compute_cell program=point_program(cell,p,links[id].spatial,spatial,selection);
            unused={};
            float value=evaluate(nodes,primitives,program,p,unused,false);
            if(point_profile_enabled)select_record(costs,depths,id,depth,unused,value>0,
                program.root!=cell.root || program.count!=cell.count,false);
            if(!isfinite(value)){result.status=COMPUTE_NONE;break;}
            if(value>0){id=next;continue;}
            bool fill=links[id].child!=COMPUTE_NONE;
            result={cell.cell_id,cell.material_id,int(depth),fill?2u:1u};
            if(!fill)break;
            device const float* m=transforms[cell.transform].inverse;
            p=float3(m[0]*p.x+m[1]*p.y+m[2]*p.z+m[3],
                     m[4]*p.x+m[5]*p.y+m[6]*p.z+m[7],
                     m[8]*p.x+m[9]*p.y+m[10]*p.z+m[11]);
            if(!all(isfinite(p))){result.status=COMPUTE_NONE;break;}
            id=links[id].child;
        }
        if(id!=COMPUTE_NONE && ++depth>=COMPUTE_HIER_STACK){result.status=COMPUTE_NONE;break;}
    }
    results[index]=result;
}

kernel void alea_select_points(
    device const compute_node* nodes [[buffer(0)]],device const compute_cell* cells [[buffer(1)]],
    device const float* xyz [[buffer(2)]],device compute_selection* results [[buffer(3)]],
    constant uint& root_or_count [[buffer(4)]],constant uint& count [[buffer(5)]],
    device const compute_primitive* primitives [[buffer(6)]],
    device const uint* candidates [[buffer(13)]],
    device atomic_uint* costs [[buffer(14)]],device atomic_uint* depths [[buffer(15)]],uint index [[thread_position_in_grid]]) {
    if(index>=count)return;
    float3 p(xyz[3*index],xyz[3*index+1],xyz[3*index+2]);
    compute_selection result={-1,-1,-1,0};compute_point_profile unused={};
    for(uint id=point_candidates_enabled?root_or_count:0;
        point_candidates_enabled?id!=COMPUTE_NONE:id<root_or_count;
        id=point_candidates_enabled?candidates[id]:id+1) {
        unused={};
        float v=evaluate(nodes,primitives,cells[id],p,unused,false);
        if(point_profile_enabled)select_record(costs,depths,id,0,unused,v>0,false,false);
        if(!isfinite(v)){result.status=COMPUTE_NONE;break;}
        if(v<=0){result={cells[id].cell_id,cells[id].material_id,0,1};break;}
    }
    results[index]=result;
}
