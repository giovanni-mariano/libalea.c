// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
// SPDX-License-Identifier: MPL-2.0
#include "compute_groups.h"
#include "alea.h"
#include <stdlib.h>
#include <string.h>
int alea_compute_group_program(const compute_scene* s,compute_cell c,
                                compute_group_emit emit,void* user) {
    struct frame {uint32_t begin,end,skip,leaves;int pure;} stack[COMPUTE_STACK];
    unsigned sp=0;
    if(!c.count || (uint64_t)c.root+c.count>s->node_count)goto invalid;
    for(uint32_t i=0;i<c.count;++i) {
        if(!(i&4095) && alea_interrupted()){alea_set_error_detail(ALEA_ERR_INTERRUPTED,"CSG group preparation interrupted");return -1;}
        compute_node n=alea_compute_node_at(s,(size_t)c.root+i);
        if(n.op==COMPUTE_LEAF || n.op==COMPUTE_CONSTANT) {
            if(sp==COMPUTE_STACK || (n.op==COMPUTE_LEAF && (n.arg>>1)>=s->primitive_count))goto invalid;
            stack[sp++]=(struct frame){i,i,COMPUTE_NONE,n.op==COMPUTE_LEAF, n.op==COMPUTE_LEAF};
        } else if(n.op==COMPUTE_SKIP_UNION || n.op==COMPUTE_SKIP_INTERSECTION || n.op==COMPUTE_SKIP_REVERSE_DIFFERENCE) {
            if(!sp || stack[sp-1].skip!=COMPUTE_NONE || stack[sp-1].end+1!=i || (uint64_t)i+n.arg+1>c.count)goto invalid;
            stack[sp-1].skip=i;
        } else if(n.op==COMPUTE_COMPLEMENT) {
            if(!sp || stack[sp-1].skip!=COMPUTE_NONE || stack[sp-1].end+1!=i)goto invalid;
            struct frame* a=stack+sp-1;
            if(a->pure && a->leaves>1 && emit(user,(compute_group_span){c.root+a->begin,c.root+a->end+1}))return -1;
            a->pure=0;a->leaves=0;a->end=i;
        } else {
            if(sp<2 || (n.op!=COMPUTE_UNION && n.op!=COMPUTE_INTERSECTION && n.op!=COMPUTE_DIFFERENCE && n.op!=COMPUTE_REVERSE_DIFFERENCE))goto invalid;
            struct frame a=stack[sp-2],b=stack[sp-1];
            if(a.skip==COMPUTE_NONE || b.skip!=COMPUTE_NONE || a.end+1!=a.skip || b.begin!=a.skip+1 || b.end+1!=i)goto invalid;
            compute_node branch=alea_compute_node_at(s,(size_t)c.root+a.skip);
            unsigned expected=n.op==COMPUTE_UNION?COMPUTE_SKIP_UNION:n.op==COMPUTE_REVERSE_DIFFERENCE?COMPUTE_SKIP_REVERSE_DIFFERENCE:COMPUTE_SKIP_INTERSECTION;
            if(branch.op!=expected || (uint64_t)a.skip+branch.arg+1!=(uint64_t)i+1)goto invalid;
            int pure=n.op==COMPUTE_INTERSECTION && a.pure && b.pure;
            if(!pure) {
                if(a.pure && a.leaves>1 && emit(user,(compute_group_span){c.root+a.begin,c.root+a.end+1}))return -1;
                if(b.pure && b.leaves>1 && emit(user,(compute_group_span){c.root+b.begin,c.root+b.end+1}))return -1;
            }
            stack[sp-2]=(struct frame){a.begin,i,COMPUTE_NONE,pure?a.leaves+b.leaves:0,pure};--sp;
        }
    }
    if(sp!=1 || stack[0].begin || stack[0].end+1!=c.count || stack[0].skip!=COMPUTE_NONE)goto invalid;
    if(stack[0].pure && stack[0].leaves>1 && emit(user,(compute_group_span){c.root,c.root+c.count}))return -1;
    return 0;
 invalid:alea_set_error_detail(ALEA_ERR_INVALID_STATE,"Invalid structured CSG group program");return -1;
}

void alea_compute_groups_free(compute_group_plan* p) {
    if(!p)return;free(p->programs);free(p->spans);memset(p,0,sizeof(*p));
}
static int compare_programs(const void* x,const void* y) {
    const compute_group_program* a=x,*b=y;
    if(a->root!=b->root)return a->root<b->root?-1:1;
    return a->count<b->count?-1:a->count!=b->count;
}
static int compare_spans(const void* x,const void* y) {
    const compute_group_span* a=x,*b=y;
    return a->begin<b->begin?-1:a->begin!=b->begin;
}
typedef struct {compute_group_plan* plan;size_t limit;int unavailable;} group_builder;
static int append_span(void* user,compute_group_span span) {
    group_builder* b=user;compute_group_plan* p=b->plan;
    if(p->span_count==p->span_capacity) {
        size_t next=p->span_capacity?p->span_capacity*2:1024;
        if(next>b->limit/sizeof(*p->spans))next=b->limit/sizeof(*p->spans);
        if(next<=p->span_capacity){b->unavailable=1;return -1;}
        compute_group_span* spans=realloc(p->spans,next*sizeof(*spans));
        if(!spans){b->unavailable=1;return -1;}p->spans=spans;p->span_capacity=next;
    }
    p->spans[p->span_count++]=span;return 0;
}
int alea_compute_groups_plan(const compute_scene* scene,size_t budget,compute_group_plan* out) {
    memset(out,0,sizeof(*out));compute_group_plan p={0};
    if(alea_interrupted())goto interrupted;
    if(scene->cell_count>SIZE_MAX-scene->spatial_count)return 0;
    size_t capacity=scene->cell_count+scene->spatial_count;
    if(!capacity || capacity>UINT32_MAX || capacity>budget/sizeof(*p.programs))return 0;
    size_t bytes=capacity*sizeof(*p.programs);
    p.programs=malloc(bytes);if(!p.programs)return 0;p.program_capacity=capacity;
    for(size_t i=0;i<capacity;++i) {
        if(!(i&4095) && alea_interrupted())goto interrupted;
        compute_cell c={0};
        if(i<scene->cell_count)c=scene->cells[i];else {
            const compute_spatial_program* s=scene->spatial+i-scene->cell_count;
            if(s->child!=COMPUTE_NONE)continue;c.root=s->root;c.count=s->count;
        }
        p.programs[p.program_count++]=(compute_group_program){c.root,c.count,0,0};
    }
    qsort(p.programs,p.program_count,sizeof(*p.programs),compare_programs);
    size_t n=0;
    for(size_t i=0;i<p.program_count;++i)if(!n || compare_programs(p.programs+i,p.programs+n-1))p.programs[n++]=p.programs[i];
    p.program_count=n;
    group_builder builder={&p,(budget-bytes)/2,0};
    if(builder.limit>24*1024*1024)builder.limit=24*1024*1024;
    for(size_t i=0;i<n;++i) {
        compute_group_program* record=p.programs+i;size_t first=p.span_count;
        compute_cell c={.root=record->root,.count=record->count};
        if(alea_compute_group_program(scene,c,append_span,&builder)) {
            if(builder.unavailable){alea_compute_groups_free(&p);return 0;}goto fail;
        }
        if(p.span_count>first)qsort(p.spans+first,p.span_count-first,sizeof(*p.spans),compare_spans);
        for(size_t j=first+1;j<p.span_count;++j)if(p.spans[j].begin<p.spans[j-1].end) {
            alea_set_error_detail(ALEA_ERR_INVALID_STATE,"Overlapping CSG groups within a program");goto fail;
        }
        record->offset=(uint32_t)first;record->groups=(uint32_t)(p.span_count-first);
    }
    if(!p.span_count){alea_compute_groups_free(&p);return 0;}
    if(alea_interrupted())goto interrupted;
    *out=p;return 1;
 interrupted:alea_set_error_detail(ALEA_ERR_INTERRUPTED,"CSG group preparation interrupted");
 fail:alea_compute_groups_free(&p);return -1;
}
