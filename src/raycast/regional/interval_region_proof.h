// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
// SPDX-License-Identifier: MPL-2.0
/* Conservative sets of possible scalar signs on one closed child box. A
 * union returns its left operand immediately on zero. Keep zero separate
 * so that proof reuse is valid in complemented expressions as well. */
enum { SIGN_NEG=1, SIGN_ZERO=2, SIGN_POS=4, SIGN_UNKNOWN=7 };
static unsigned reverse_signs(unsigned mask) {
    return ((mask&SIGN_NEG)<<2) | (mask&SIGN_ZERO) | ((mask&SIGN_POS)>>2);
}
static unsigned intersection_signs(unsigned a,unsigned b) {
    unsigned mask=(a&b&SIGN_NEG) | ((a|b)&SIGN_POS);
    if(((a&SIGN_ZERO) && (b&(SIGN_NEG|SIGN_ZERO))) ||
       ((b&SIGN_ZERO) && (a&(SIGN_NEG|SIGN_ZERO))))mask|=SIGN_ZERO;
    return mask;
}
static unsigned region_signs(interval_index_t* idx,uint32_t id,interval_box_t box,
                             unsigned depth,unsigned* remaining) {
    if(id>=idx->sys->nodes.count)return SIGN_POS;
    if(!*remaining || depth>=128)return SIGN_UNKNOWN;
    struct proof_entry* slot=idx->proof_table+(id&65535);
    if(slot->stamp==idx->proof_epoch && slot->node==id)return slot->signs;
    --*remaining;
    const alea_node_t* n=idx->sys->nodes.data+id;
    interval_box_t bound=interval_box_from(alea_node_bbox_get(&n->bbox));
    uint32_t certificate=idx->guards.bound_ids[id];
    if(certificate) {
        interval_box_t c=interval_box_from(alea_node_bbox_get(idx->guards.certified+certificate-1));
        for(unsigned k=0;k<3;k++)
            if(box.hi[k]<c.lo[k] || box.lo[k]>c.hi[k])return SIGN_POS;
    }
    bool partial=false;
    for(unsigned k=0;k<3;k++) {
        if(box.hi[k]<bound.lo[k] || box.lo[k]>bound.hi[k])return SIGN_POS;
        if(box.lo[k]<bound.lo[k] || box.hi[k]>bound.hi[k])partial=true;
    }
    unsigned result=SIGN_UNKNOWN;
    alea_operation_t op=ALEA_GET_OPERATION(n);
    if(op==ALEA_OP_PRIMITIVE) {
        uint32_t p=n->primitive.primitive_id;
        if(p<idx->sys->primitives.count && idx->surfaces[p].kind) {
            struct shape_proof* proof=idx->primitive_proofs+p;
            if(proof->stamp!=idx->proof_epoch) {
                int sign=guard_box_sign(idx->surfaces+p,box);
                proof->signs=sign>0?SIGN_POS:sign<0?SIGN_NEG:SIGN_UNKNOWN;
                proof->stamp=idx->proof_epoch;
            }
            result=proof->signs;
            if((n->primitive.sense>0)!=(n->primitive.inverted!=0))result=reverse_signs(result);
        }
    }else if(op==ALEA_OP_COMPLEMENT) {
        result=reverse_signs(region_signs(idx,n->operation.left,box,depth+1,remaining));
    }else if(op==ALEA_OP_UNION) {
        unsigned left=region_signs(idx,n->operation.left,box,depth+1,remaining);
        result=left&(SIGN_NEG|SIGN_ZERO);
        if(left&SIGN_POS) result|=region_signs(idx,n->operation.right,box,depth+1,remaining);
    }else if(op==ALEA_OP_INTERSECTION || op==ALEA_OP_DIFFERENCE) {
        unsigned right=region_signs(idx,n->operation.right,box,depth+1,remaining);
        if(op==ALEA_OP_DIFFERENCE)right=reverse_signs(right);
        result=right==SIGN_POS?SIGN_POS:
            intersection_signs(region_signs(idx,n->operation.left,box,depth+1,remaining),right);
    }
    if(partial)result|=SIGN_POS;
    *slot=(struct proof_entry){idx->proof_epoch,id,result};
    return result;
}
