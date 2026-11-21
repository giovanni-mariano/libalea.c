// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
// SPDX-License-Identifier: MPL-2.0
/* Exact sign predicates within a certified spatial box. An unexpanded
 * subtree refers back to the original scalar evaluator. AND/OR operands
 * are shared between the prepared cells and leaves. */
static uint32_t guard_hash(uint32_t a,uint32_t b) {
    uint32_t x = a*UINT32_C(0x9e3779b1) ^ b*UINT32_C(0x85ebca77);
    x ^= x>>16; x *= UINT32_C(0x21f0aaad); return x^(x>>15);
}
static bool guard_grow(guard_pool_t* p) {
    size_t cap=p->capacity ? p->capacity*2 : 1024;
    uint32_t* table=calloc(cap,sizeof(*table));
    if (!table) return false;
    for (size_t i=0;i<p->pairs.count;i++) {
        guard_pair_t k=p->pairs.data[i];size_t s=guard_hash(k.a,k.b)&(cap-1);
        while (table[s]) s=(s+1)&(cap-1);
        table[s]=(uint32_t)i+1;
    }
    free(p->table);p->table=table;p->capacity=cap;
    return true;
}

static uint32_t guard_join(guard_pool_t* p,uint32_t a,uint32_t b,bool join) {
    uint32_t absorbing=join?GUARD_ANY:GUARD_EMPTY, neutral=join?GUARD_EMPTY:GUARD_ANY;
    if (a==absorbing || b==absorbing) return absorbing;
    if (a==neutral || a==b) return b;
    if (b==neutral) return a;
    if ((p->pairs.count+1)*2>=p->capacity && !guard_grow(p)) {
        p->failed=true;return GUARD_ANY;
    }
    size_t s=guard_hash(a,b)&(p->capacity-1);
    while (p->table[s]) {
        guard_pair_t k=p->pairs.data[p->table[s]-1];
        if (k.a==a && k.b==b) break;
        s=(s+1)&(p->capacity-1);
    }
    if (!p->table[s]) {
        guard_pair_t k={a,b};
        if (p->pairs.count>=GUARD_UNION || alea_vec_push(&p->pairs,k,guard_pair_t)) {
            p->failed=true;return GUARD_ANY;
        }
        p->table[s]=(uint32_t)p->pairs.count;
    }
    return GUARD_BRANCH | (join?GUARD_UNION:0) | (p->table[s]-1);
}

static uint32_t compile_guard(const alea_system_t* s,guard_pool_t* pool,
                              const guard_surface_t* surfaces,uint32_t id,interval_box_t box,
                              bool positive,unsigned depth,unsigned* budget) {
    if (id>=s->nodes.count) return positive?GUARD_ANY:GUARD_EMPTY;
    uint32_t ref=GUARD_REF|(2*id+positive);
    if (!*budget || depth>=768) return ref;
    --*budget;
    uint32_t known=pool->bound_ids?pool->bound_ids[id]:0;
    if(known) {
        interval_box_t tight=interval_box_from(alea_node_bbox_get(pool->certified+known-1));
        for(unsigned k=0;k<3;k++)
            if(box.hi[k]<tight.lo[k] || box.lo[k]>tight.hi[k])return positive?GUARD_ANY:GUARD_EMPTY;
    }
    const alea_node_t* n=&s->nodes.data[id];
    interval_box_t bounds=interval_box_from(alea_node_bbox_get(&n->bbox));
    bool partial=false;
    for (unsigned i=0;i<3;i++) {
        if (box.lo[i]>bounds.hi[i] || box.hi[i]<bounds.lo[i]) return positive?GUARD_ANY:GUARD_EMPTY;
        if (box.lo[i]<bounds.lo[i]) {box.lo[i]=bounds.lo[i];partial=true;}
        if (box.hi[i]>bounds.hi[i]) {box.hi[i]=bounds.hi[i];partial=true;}
    }
    alea_operation_t op=ALEA_GET_OPERATION(n);
    uint32_t result;
    if (op==ALEA_OP_PRIMITIVE) {
        uint32_t prim=n->primitive.primitive_id;
        if (prim>=s->primitives.count) return ref;
        bool invert=((n->primitive.sense>0)!=(n->primitive.inverted!=0)) != positive;
        int sign=surfaces[prim].kind?guard_box_sign(&surfaces[prim],box):0;
        result=sign ? (((sign<0)!=invert) ? GUARD_ANY : GUARD_EMPTY) : 2+2*prim+invert;
    } else if (op==ALEA_OP_COMPLEMENT) {
        result=compile_guard(s,pool,surfaces,n->operation.left,box,!positive,depth+1,budget);
    } else {
        bool join=(op==ALEA_OP_UNION)!=positive;
        /* In an intersection, evaluate the right operand first so a short
         * clipping expression can bound the usually left-deep remainder. */
        if (!positive && op==ALEA_OP_INTERSECTION) {
            uint32_t right=compile_guard(s,pool,surfaces,n->operation.right,box,false,depth+1,budget);
            result=right;
            if (right!=GUARD_EMPTY) {
                /* A right operand proven true covers the entire current box. */
                interval_box_t inner=box;
                if (right!=GUARD_ANY) {
                    unsigned steps=32;
                    inner=contract_cell(s,n->operation.right,box,false,depth+1,&steps);
                }
                uint32_t left=box_empty(inner)?GUARD_EMPTY:
                    compile_guard(s,pool,surfaces,n->operation.left,inner,false,depth+1,budget);
                result=guard_join(pool,right,left,false);
            }
        } else {
        uint32_t a=compile_guard(s,pool,surfaces,n->operation.left,box,positive,depth+1,budget);
        if (a==(join?GUARD_ANY:GUARD_EMPTY)) result=a;
        else {
            uint32_t b=compile_guard(s,pool,surfaces,n->operation.right,box,
                                     positive!=(op==ALEA_OP_DIFFERENCE),depth+1,budget);
            if (positive && op==ALEA_OP_UNION && b!=GUARD_ANY) {
                /* For a union, the scalar evaluator returns the left value
                 * immediately if it is zero. Accept that case as well. */
                uint32_t left_nonpositive=compile_guard(s,pool,surfaces,n->operation.left,
                    box,false,depth+1,budget);
                b=guard_join(pool,left_nonpositive,b,true);
            }
            result=guard_join(pool,a,b,join);
        }
        }
    }
    if (partial) result=guard_join(pool,GUARD_BOX|(2*id+positive),result,positive);
    return result;
}
