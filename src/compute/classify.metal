// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
// SPDX-License-Identifier: MPL-2.0

/* Native Metal point traversal. These links include void; render BVHs do not.
 * One thread owns a point and counts every terminal occurrence under matching
 * ancestors. Select only point-certified spatial programs; retain the original
 * expression outside their domains. Render visibility bounds never reject owners. */
struct point_frame {
    float3 point;
    uint next;
};

static compute_cell point_program(compute_cell cell,float3 p,uint root,
                                   device const compute_spatial_program* spatial,thread uint& selection) {
    if(root==COMPUTE_NONE) {
        if(point_profile_enabled)selection=cell.spatial==COMPUTE_NONE?POINT_NO_PROGRAM:POINT_UNCERTIFIED;
        return cell;
    }
    for(uint a=0;a<3;++a)
        if(p[a]<spatial[root].bounds[2*a] || p[a]>spatial[root].bounds[2*a+1]) {
            if(point_profile_enabled)selection=POINT_OUTSIDE_DOMAIN;
            return cell;
        }
    uint id=root;
    while(spatial[id].child!=COMPUTE_NONE) {
        uint child=spatial[id].child,axis=spatial[id].axis;
        id=child+uint(p[axis]>=spatial[child+1].bounds[2*axis]);
    }
    if(point_profile_enabled)selection=cell.root==spatial[id].root && cell.count==spatial[id].count?
        POINT_UNCHANGED_REGION:POINT_SPATIAL_CALLS;
    cell.root=spatial[id].root;cell.count=spatial[id].count;
    return cell;
}

/* Point ownership uses floor, without the directional crossing bias used by
 * ray traversal. Select one element; its terminal cells still count separately. */
static uint point_lattice(device const compute_lattice& l, device const uint* fills,
                          thread float3& p, thread bool& error) {
    int3 relative;
    bool outside = false;
    for (uint a = 0; a < 3; ++a) {
        int length = l.dims[2*a+1] - l.dims[2*a] + 1;
        bool tiled = a == 0 || length > 1 || (l.repeating && (a != 2 || l.pitch[a] > 0));
        if (!l.repeating && l.outer == COMPUTE_NONE && length == 1 && l.pitch[a] > 0 &&
            (p[a] < l.lower[a] || p[a] > l.lower[a] + l.pitch[a])) return COMPUTE_NONE;
        float index = tiled ? floor((p[a] - l.lower[a]) / l.pitch[a]) : 0;
        /* FP32 must distinguish adjacent element indices. Reject rather than
         * silently wrap integers or collapse distant repeated placements. */
        if (!isfinite(index) || abs(index) > 16777215.0f) { error = true; return COMPUTE_NONE; }
        relative[a] = int(index);
        long logical = long(relative[a]) + long(l.dims[2*a]);
        if (logical < -16777215 || logical > 16777215) { error = true; return COMPUTE_NONE; }
        outside |= relative[a] < 0 || relative[a] >= length;
        float shift = l.zero_coords ? float(logical) * l.pitch[a] :
            a == 2 && !tiled ? 0 : l.lower[a] + (index + .5f) * l.pitch[a];
        p[a] -= shift;
        if (!isfinite(p[a])) { error = true; return COMPUTE_NONE; }
    }
    if (l.repeating) return fills[l.offset];
    if (outside) return l.outer;
    uint ny = uint(l.dims[3] - l.dims[2] + 1), nz = uint(l.dims[5] - l.dims[4] + 1);
    return fills[l.offset + (uint(relative.x) * ny + uint(relative.y)) * nz + uint(relative.z)];
}

kernel void alea_classify_hierarchy(
    device const compute_node* nodes [[buffer(0)]],
    device const compute_cell* cells [[buffer(1)]],
    device const float* xyz [[buffer(2)]],
    device compute_result* results [[buffer(3)]],
    constant uint& root [[buffer(4)]],
    constant uint& count [[buffer(5)]],
    device const compute_primitive* primitives [[buffer(6)]],
    device const compute_transform* transforms [[buffer(7)]],
    device const metal_point_cell* links [[buffer(8)]],
    device const compute_lattice* lattices [[buffer(9)]],
    device const uint* fills [[buffer(10)]],
    device const compute_spatial_program* spatial [[buffer(11)]],
    device compute_point_profile* profiles [[buffer(12)]],
    uint index [[thread_position_in_grid]]) {
    if (index >= count) return;
    point_frame stack[COMPUTE_HIER_STACK];
    uint top = 0;
    stack[0].point = float3(xyz[3 * index], xyz[3 * index + 1], xyz[3 * index + 2]);
    stack[0].next = root;
    compute_result result = {-1, -1, 0};
    compute_point_profile cost={};
    if(point_profile_enabled){cost.v[POINT_UNIVERSES]=1;cost.v[POINT_MAX_DEPTH]=1;}
    while (true) {
        uint id = stack[top].next;
        if (id == COMPUTE_NONE) {
            if (top == 0) break;
            --top;
            continue;
        }
        stack[top].next = links[id].next;
        compute_cell cell = cells[id];
        if(point_profile_enabled)point_add(cost.v[POINT_CELLS]);
        float3 p = stack[top].point;
        if (cell.lattice != COMPUTE_NONE) {
            if(point_profile_enabled)point_add(cost.v[POINT_LATTICES]);
            bool error = false;
            uint child = point_lattice(lattices[cell.lattice], fills, p, error);
            if (error) { result.owner_count = COMPUTE_NONE; break; }
            if (child == COMPUTE_NONE) continue;
            if (top + 1 == COMPUTE_HIER_STACK) { result.owner_count = COMPUTE_NONE; break; }
            ++top;
            stack[top].point = p;
            stack[top].next = child;
            if(point_profile_enabled) {
                point_add(cost.v[POINT_DESCENTS]);point_add(cost.v[POINT_UNIVERSES]);
                cost.v[POINT_MAX_DEPTH]=max(cost.v[POINT_MAX_DEPTH],top+1);
            }
            continue;
        }
        uint selection=POINT_NO_PROGRAM;
        compute_cell program=point_program(cell,p,links[id].spatial,spatial,selection);
        bool specialized=program.root!=cell.root || program.count!=cell.count;
        if(point_profile_enabled) {
            point_add(cost.v[specialized?POINT_SPATIAL_CALLS:POINT_ORIGINAL_CALLS]);
            if(!specialized)point_add(cost.v[selection]);
        }
        uint slot=specialized?POINT_SPATIAL_INSTRUCTIONS:POINT_ORIGINAL_INSTRUCTIONS;
        uint before=cost.v[slot];
        float value = evaluate(nodes, primitives, program, p, cost, specialized);
        if(point_profile_enabled) {
            point_add(cost.v[value>0?POINT_REJECTED_INSTRUCTIONS:POINT_ACCEPTED_INSTRUCTIONS],cost.v[slot]-before);
            if(value>0)point_add(cost.v[POINT_REJECTED]);
        }
        if (!isfinite(value)) { result.owner_count = COMPUTE_NONE; break; }
        if (value > 0) continue;
        if (links[id].child != COMPUTE_NONE) {
            if (top + 1 == COMPUTE_HIER_STACK) { result.owner_count = COMPUTE_NONE; break; }
            device const float* m = transforms[cell.transform].inverse;
            float3 local(m[0]*p.x + m[1]*p.y + m[2]*p.z + m[3],
                         m[4]*p.x + m[5]*p.y + m[6]*p.z + m[7],
                         m[8]*p.x + m[9]*p.y + m[10]*p.z + m[11]);
            if (!all(isfinite(local))) { result.owner_count = COMPUTE_NONE; break; }
            ++top;
            stack[top].point = local;
            stack[top].next = links[id].child;
            if(point_profile_enabled) {
                point_add(cost.v[POINT_DESCENTS]);point_add(cost.v[POINT_UNIVERSES]);
                cost.v[POINT_MAX_DEPTH]=max(cost.v[POINT_MAX_DEPTH],top+1);
            }
        } else {
            if (result.owner_count == COMPUTE_NONE - 1) { result.owner_count = COMPUTE_NONE; break; }
            ++result.owner_count;
            if(point_profile_enabled)point_add(cost.v[POINT_OWNERS]);
            result.cell_id = cell.cell_id;
            result.material_id = cell.material_id;
        }
    }
    if (result.owner_count != 1) result.cell_id = result.material_id = -1;
    results[index] = result;
    if(point_profile_enabled)profiles[index]=cost;
}
