// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
// SPDX-License-Identifier: MPL-2.0
#include "compute_internal.h"
#include <stdlib.h>
#include <string.h>

/* Spatial CSG programs simplify the original expression, independently of
 * render BVH visibility. Point queries may use them only after proving every
 * original primitive finite throughout the tree's logical domain. Inspecting
 * even short-circuited operands is deliberately conservative: an outside
 * Boolean result must not hide a nonfinite left operand. Original programs and
 * all-owner cell lists remain the fallback outside certified regions. */
#define POINT_PROOF_BUDGET (32u*1024u*1024u)
static int program_finite(const compute_scene* s,compute_cell cell,const float b[6]) {
    for(uint32_t i=0;i<cell.count;++i) {
        compute_node n=s->nodes[cell.root+i];
        if(n.op==COMPUTE_LEAF) {
            int finite;
            alea_compute_primitive_sign_finite(&s->primitives[n.arg>>1],b,&finite);
            if(!finite)return 0;
        }
    }
    return 1;
}
int alea_compute_point_programs_pack(compute_scene* s) {
    const char* mode=getenv("ALEA_METAL_POINT_SPATIAL");
    if(mode && strcmp(mode,"0") && strcmp(mode,"1")) {
        alea_set_error_detail(ALEA_ERR_INVALID_ARG,"ALEA_METAL_POINT_SPATIAL must be 0 or 1");return -1;
    }
    if((mode && !strcmp(mode,"0")) || !s->point_cells)return 0;
    uint32_t* roots=malloc(s->cell_count*sizeof(*roots));
    if(!roots){alea_set_error_detail(ALEA_ERR_OUT_OF_MEMORY,"Certifying point programs");return -1;}
    size_t budget=POINT_PROOF_BUDGET;
    for(size_t i=0;i<s->cell_count;++i) {
        if(alea_interrupted()) {
            free(roots);alea_set_error_detail(ALEA_ERR_INTERRUPTED,"Point program certification interrupted");return -1;
        }
        roots[i]=COMPUTE_NONE;
        compute_cell cell=s->cells[i];
        if(cell.spatial==COMPUTE_NONE || cell.lattice!=COMPUTE_NONE || cell.count>budget)continue;
        budget-=cell.count;
        if(program_finite(s,cell,s->spatial[cell.spatial].bounds))roots[i]=cell.spatial;
    }
    for(size_t i=0;i<s->cell_count;++i)s->point_cells[i].spatial=roots[i];
    free(roots);return 0;
}
