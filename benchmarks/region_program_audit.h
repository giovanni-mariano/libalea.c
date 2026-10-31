// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
// SPDX-License-Identifier: MPL-2.0
#ifndef ALEA_REGION_PROGRAM_AUDIT_H
#define ALEA_REGION_PROGRAM_AUDIT_H
#include "compute/compute_internal.h"
/* Diagnostic recognizer for validated, unthreaded packed programs. Return 1 for
 * the strict signed-leaf AND family, 0 for unsupported structure, -1 for invalid
 * ranges/primitive references. No geometry is rewritten or evaluated. */
static int region_and_program(const compute_scene* s,compute_cell c,uint32_t* leaves) {
    struct {uint32_t begin,end,skip;} stack[COMPUTE_STACK];
    unsigned sp=0;uint32_t count=0;
    *leaves=0;
    if(!c.count || (uint64_t)c.root+c.count>s->node_count)return -1;
    for(uint32_t i=0;i<c.count;++i) {
        compute_node n=alea_compute_node_at(s,(size_t)c.root+i);
        if(n.op==COMPUTE_LEAF) {
            if((n.arg>>1)>=s->primitive_count)return -1;
            if(sp==COMPUTE_STACK)return 0;
            stack[sp].begin=stack[sp].end=i;stack[sp++].skip=COMPUTE_NONE;++count;
        } else if(n.op==COMPUTE_SKIP_INTERSECTION) {
            if(!sp || stack[sp-1].skip!=COMPUTE_NONE || stack[sp-1].end+1!=i)return 0;
            if((uint64_t)i+n.arg+1>c.count)return -1;
            stack[sp-1].skip=i;
        } else if(n.op==COMPUTE_INTERSECTION) {
            if(sp<2)return 0;
            unsigned a=sp-2,b=sp-1;uint32_t skip=stack[a].skip;
            if(skip==COMPUTE_NONE || stack[b].skip!=COMPUTE_NONE ||
                stack[b].begin!=skip+1 || stack[b].end+1!=i)return 0;
            compute_node branch=alea_compute_node_at(s,(size_t)c.root+skip);
            if((uint64_t)skip+branch.arg+1!=(uint64_t)i+1)return 0;
            stack[a].end=i;stack[a].skip=COMPUTE_NONE;--sp;
        } else return 0;
    }
    if(sp!=1 || stack[0].begin || stack[0].end+1!=c.count || stack[0].skip!=COMPUTE_NONE)return 0;
    *leaves=count;return 1;
}
#endif
