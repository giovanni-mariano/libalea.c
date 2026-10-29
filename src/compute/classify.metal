// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
// SPDX-License-Identifier: MPL-2.0

/* Native Metal point traversal. These links include void; render BVHs do not.
 * One thread owns a point and counts every terminal occurrence under matching
 * ancestors. Keep the original CSG program, without render-only simplification. */
struct point_frame {
    float3 point;
    uint next;
};

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
    uint index [[thread_position_in_grid]]) {
    if (index >= count) return;
    point_frame stack[COMPUTE_HIER_STACK];
    uint top = 0;
    stack[0].point = float3(xyz[3 * index], xyz[3 * index + 1], xyz[3 * index + 2]);
    stack[0].next = root;
    compute_result result = {-1, -1, 0};
    while (true) {
        uint id = stack[top].next;
        if (id == COMPUTE_NONE) {
            if (top == 0) break;
            --top;
            continue;
        }
        stack[top].next = links[id].next;
        compute_cell cell = cells[id];
        float3 p = stack[top].point;
        float value = evaluate(nodes, primitives, cell, p);
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
        } else {
            if (result.owner_count == COMPUTE_NONE - 1) { result.owner_count = COMPUTE_NONE; break; }
            ++result.owner_count;
            result.cell_id = cell.cell_id;
            result.material_id = cell.material_id;
        }
    }
    if (result.owner_count != 1) result.cell_id = result.material_id = -1;
    results[index] = result;
}
