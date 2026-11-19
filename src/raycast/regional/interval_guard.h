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

/* Polynomial in world or local coordinates Ax²+By²+Cz²+Dxy+Eyz+Fxz+
 * Gx+Hy+Iz+J. Unsupported primitives retain an unrestricted predicate. */
enum { GUARD_RADIAL=3, GUARD_TORUS=4 };
static void guard_surface_init(guard_surface_t* out, const alea_system_t* s, uint32_t id) {
    const void* p=alea_primitive_payload_const(s,id);
    if (!p) return;
    double* k=out->k;
    switch (s->primitives.data[id].type) {
    case ALEA_PRIMITIVE_PLANE: {
        const alea_plane_data_t* plane=p;
        k[6]=plane->a;k[7]=plane->b;k[8]=plane->c;k[9]=plane->d;out->kind=1;break;
    }
    case ALEA_PRIMITIVE_QUADRIC: memcpy(k,((const alea_quadric_data_t*)p)->coeffs,sizeof(out->k));out->kind=2;break;
    case ALEA_PRIMITIVE_SPHERE: case ALEA_PRIMITIVE_SPH: {
        const alea_sphere_data_t* v=p;
        k[0]=v->center_x;k[1]=v->center_y;k[2]=v->center_z;k[3]=v->radius;
        out->kind=GUARD_RADIAL;out->axis=3;break;
    }
    case ALEA_PRIMITIVE_CYLINDER_X: {
        const alea_cylinder_x_data_t* v=p;k[1]=v->center_y;k[2]=v->center_z;k[3]=v->radius;
        out->kind=GUARD_RADIAL;out->axis=0;break;
    }
    case ALEA_PRIMITIVE_CYLINDER_Y: {
        const alea_cylinder_y_data_t* v=p;k[0]=v->center_x;k[2]=v->center_z;k[3]=v->radius;
        out->kind=GUARD_RADIAL;out->axis=1;break;
    }
    case ALEA_PRIMITIVE_CYLINDER_Z: {
        const alea_cylinder_z_data_t* v=p;k[0]=v->center_x;k[1]=v->center_y;k[3]=v->radius;
        out->kind=GUARD_RADIAL;out->axis=2;break;
    }
    case ALEA_PRIMITIVE_TORUS_X: case ALEA_PRIMITIVE_TORUS_Y: case ALEA_PRIMITIVE_TORUS_Z: {
        const alea_torus_data_t* v=p;
        double radial2=v->minor_radius*v->minor_radius, axial2=v->axial_semiwidth_B*v->axial_semiwidth_B;
        if(!(v->major_radius>v->minor_radius && v->minor_radius>0 && v->axial_semiwidth_B>0) ||
           !isnormal(radial2) || !isnormal(axial2))break;
        k[0]=v->center_x;k[1]=v->center_y;k[2]=v->center_z;
        k[3]=v->major_radius;k[4]=v->minor_radius;k[5]=v->axial_semiwidth_B;
        out->kind=GUARD_TORUS;out->axis=s->primitives.data[id].type-ALEA_PRIMITIVE_TORUS_X;
        break;
    }
    default: break;
    }
}

/* Distance to an axis or a center is bounded by the nearest and farthest
 * box coordinates. For a ring torus, bound the radial and axial terms of
 * its native elliptical cross-section separately. */
static int radial_box_sign(const guard_surface_t* s,interval_box_t box) {
    double low2=0,high2=0,axial_low=0,axial_high=0;
    for(unsigned i=0;i<3;i++) {
        double center=s->k[i];
        double pad=64*DBL_EPSILON*(1+fabs(center)+fabs(box.lo[i])+fabs(box.hi[i]));
        double nearest=fmax(0,fmax(box.lo[i]-center,center-box.hi[i])-pad);
        double farthest=fmax(fabs(box.lo[i]-center),fabs(box.hi[i]-center))+pad;
        if(i==s->axis){axial_low=nearest;axial_high=farthest;}
        else {low2+=nearest*nearest;high2+=farthest*farthest;}
    }
    double threshold=s->k[3]*s->k[3];
    if(s->kind==GUARD_TORUS) {
        double lower=sqrt(low2),upper=sqrt(high2),major=s->k[3];
        double pad=64*DBL_EPSILON*(1+upper+major);
        double near=fmax(0,fmax(lower-major,major-upper)-pad)/s->k[4];
        double far=(fmax(fabs(lower-major),fabs(upper-major))+pad)/s->k[4];
        axial_low/=s->k[5];axial_high/=s->k[5];
        low2=near*near+axial_low*axial_low;
        high2=far*far+axial_high*axial_high;
        threshold=1;
    }
    double pad=256*(DBL_EPSILON*(high2+threshold)+DBL_MIN);
    if(!isfinite(pad))return 0;
    return low2>threshold+pad?1:high2<threshold-pad?-1:0;
}

/* Centered Taylor enclosure of the scalar polynomial on a closed box.
 * The absolute coefficient sum bounds floating point evaluation error,
 * including cancellation when shifting a quadric to the box center. */
static int guard_box_sign(const guard_surface_t* s,interval_box_t box) {
    if(s->kind>=GUARD_RADIAL)return radial_box_sign(s,box);
    const double* k=s->k;
    double c[3],h[3],m[3];
    for (unsigned i=0;i<3;i++) {
        c[i]=0.5*box.lo[i]+0.5*box.hi[i];h[i]=0.5*(box.hi[i]-box.lo[i]);
        m[i]=fmax(fabs(box.lo[i]),fabs(box.hi[i]));
    }
    double scale=fabs(k[9]), radius=0, value;
    for (unsigned i=0;i<3;i++) scale+=fabs(k[i+6])*m[i];
    if (s->kind==1) {
        value=k[6]*c[0]+k[7]*c[1]+k[8]*c[2]+k[9];
        for (unsigned i=0;i<3;i++) radius+=fabs(k[i+6])*h[i];
    } else {
        value=guard_value(k,c[0],c[1],c[2]);
        for (unsigned i=0;i<3;i++) {
            unsigned j=(i+1)%3,b=(i+2)%3;
            double grad=2*k[i]*c[i]+k[3+i]*c[j]+k[3+b]*c[b]+k[6+i];
            scale+=fabs(k[i])*m[i]*m[i]+fabs(k[3+i])*m[i]*m[j];
            radius+=fabs(grad)*h[i]+fabs(k[i])*h[i]*h[i]+fabs(k[3+i])*h[i]*h[j];
        }
    }
    double margin=radius+256*(scale*DBL_EPSILON+DBL_MIN);
    return value>margin ? 1 : value < -margin ? -1 : 0;
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
