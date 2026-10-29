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

/* Strict sign propagation over the original bytecode. Inspect all operands,
 * including skipped ones, before accepting a rejection proof: a Boolean proof
 * alone must never suppress a nonfinite error in the original evaluation.
 * Signs are -1/0/+1; the two cache sentinels are never propagated into CSG.
 * The source packer validates opcodes and the maximum value-stack depth. */
enum { POINT_NONFINITE=2, POINT_UNCHECKED=3 };
static int candidate_outside(const compute_scene* s,compute_cell cell,const float bounds[6],
    int8_t* signs,compute_candidate_stats* stats) {
    int stack[COMPUTE_STACK];unsigned top=0;
    for(uint32_t i=0;i<cell.count;++i) {
        if((i&1023)==0 && alea_interrupted())return -1;
        compute_node n=s->nodes[cell.root+i];++stats->instructions;
        if(n.op==COMPUTE_LEAF) {
            uint32_t primitive=n.arg>>1;
            if(signs[primitive]==POINT_UNCHECKED) {
                int finite,sign=alea_compute_primitive_sign_finite(&s->primitives[primitive],bounds,&finite);
                signs[primitive]=finite?sign:POINT_NONFINITE;++stats->primitives;
            }
            if(signs[primitive]==POINT_NONFINITE)return 0;
            stack[top++]=(n.arg&1)?-signs[primitive]:signs[primitive];
        } else if(n.op==COMPUTE_CONSTANT)stack[top++]=n.arg?1:-1;
        else if(n.op==COMPUTE_COMPLEMENT)stack[top-1]=-stack[top-1];
        else if(n.op==COMPUTE_SKIP_UNION || n.op==COMPUTE_SKIP_INTERSECTION || n.op==COMPUTE_SKIP_REVERSE_DIFFERENCE)continue;
        else {
            int right=stack[--top],left=stack[top-1];
            if(n.op==COMPUTE_DIFFERENCE)right=-right;
            if(n.op==COMPUTE_REVERSE_DIFFERENCE)left=-left;
            stack[top-1]=n.op==COMPUTE_UNION?(left<right?left:right):(left>right?left:right);
        }
    }
    return stack[0]>0;
}
int alea_compute_point_candidates(const compute_scene* s,const float bounds[6],
    int8_t* signs,uint32_t* next,uint32_t* root,compute_candidate_stats* stats) {
    memset(signs,POINT_UNCHECKED,s->primitive_count*sizeof(*signs));
    *root=COMPUTE_NONE;uint32_t previous=COMPUTE_NONE;
    size_t budget=POINT_PROOF_BUDGET;
    uint32_t id=s->point_cells?s->point_root:s->cell_count?0:COMPUTE_NONE;
    while(id!=COMPUTE_NONE) {
        if(alea_interrupted())goto interrupted;
        compute_cell cell=s->cells[id];int outside=0;++stats->cells;
        /* A lattice's element CSG is not its placement envelope. Keep it. */
        if(cell.lattice==COMPUTE_NONE) {
            if(cell.count<=budget) {
                budget-=cell.count;
                outside=candidate_outside(s,cell,bounds,signs,stats);
                if(outside<0)goto interrupted;
            } else ++stats->budget_fallbacks;
        }
        if(outside)++stats->rejected;
        else {
            if(previous==COMPUTE_NONE)*root=id;else next[previous]=id;
            next[id]=COMPUTE_NONE;previous=id;
        }
        id=s->point_cells?s->point_cells[id].next:id+1<s->cell_count?id+1:COMPUTE_NONE;
    }
    return 0;
interrupted:
    alea_set_error_detail(ALEA_ERR_INTERRUPTED,"Point candidate proof interrupted");return -1;
}
