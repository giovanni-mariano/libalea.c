// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
// SPDX-License-Identifier: MPL-2.0
/* A required predicate can narrow the interval before solving any physical
 * crossing. The chord minus max(A,0)*width²/4 bounds a quadratic below. */
static inline ray_gate_t* ray_quadric_gate(interval_level_t* l,uint32_t id) {
    const interval_index_t* idx=l->query->index;
    regional_worker_t* worker=l->query->worker;
    ray_gate_t* g=worker->quadric_gates+id;
    if(g->stamp!=worker->generation) {
        const query_quadric_t* q=idx->quadrics+id;
        ray_quadratic_t poly=quadric_polynomial(&l->local,q->data.coeffs,1);
        g->a=poly.a;g->b=poly.b;g->c=poly.c;
        double m=l->magnitude;
        g->error=1024*(DBL_EPSILON*((q->quadratic_scale*m+q->linear_scale)*m+fabs(q->data.coeffs[9]))+DBL_MIN);
        g->stamp=worker->generation;
    }
    return g;
}
/* 0 and 1 are certified false and true, 2 is undecided. If this atom
 * is required, discard only the part that is provably false. */
static int classify_ray_atom(interval_level_t* l,uint32_t code,double* lo,double* hi,bool required) {
    const interval_index_t* idx=l->query->index;
    regional_worker_t* worker=l->query->worker;
    uint32_t id=(code&QUERY_OPERAND)>>1;
    ray_gate_t* g;
    const alea_ray_t* ray=&l->local;
    if((code&~QUERY_OPERAND)==QUERY_PLANE) {
        g=worker->plane_gates+id;
        if(g->stamp!=worker->generation) {
            const alea_plane_data_t* p=&idx->planes[id].data;
            g->a=0;
            g->b=p->a*ray->dx+p->b*ray->dy+p->c*ray->dz;
            g->c=p->a*ray->ox+p->b*ray->oy+p->c*ray->oz+p->d;
            g->error=512*(DBL_EPSILON*((fabs(p->a)+fabs(p->b)+fabs(p->c))*l->magnitude+fabs(p->d))+DBL_MIN);
            g->stamp=worker->generation;
        }
    }else if((code&~QUERY_OPERAND)==QUERY_QUADRIC) {
        g=ray_quadric_gate(l,id);
    }else if((code&(GUARD_REF|GUARD_BOX))==GUARD_BOX) {
        const alea_node_bbox_t* b=idx->boxes.data+((code&~GUARD_BOX)>>1);
        double low[]={b->min_x,b->min_y,b->min_z},high[]={b->max_x,b->max_y,b->max_z};
        bool inside=true;
        for(unsigned k=0;k<3;k++) {
            double v=l->direction[k],o=l->origin[k],pad=l->pad[k];
            double x=o+v*(v<0?*hi:*lo),y=o+v*(v<0?*lo:*hi);
            if(y<low[k]-pad || x>high[k]+pad)return (int)(code&1);
            if(x<low[k]+pad || y>high[k]-pad)inside=false;
        }
        return inside?(int)!(code&1):2;
    }else return code<=GUARD_ANY?(int)code:2;
    double sign=(code&1)?-1:1, width=*hi-*lo;
    double a=g->a*sign, b=g->b*sign, c=g->c*sign;
    double left=(a*(*lo)+b)*(*lo)+c,right=(a*(*hi)+b)*(*hi)+c;
    double bow=a*width*width*0.25;
    if(!isfinite(left) || !isfinite(right) || !isfinite(bow) || !isfinite(g->error))return 2;
    double upper_margin=g->error+(bow<0?-bow:0);
    if(left < -upper_margin && right < -upper_margin)return 1;
    double lower_margin=g->error+(bow>0?bow:0);
    left-=lower_margin;right-=lower_margin;
    if(!isfinite(left) || !isfinite(right))return 2;
    if(left>0 && right>0)return 0;
    if(required && ((left>0)!=(right>0))) {
        double delta=left-right;
        if(!isfinite(delta))return 2;
        double crossing=*lo+width*(left/delta);
        double pad=64*DBL_EPSILON*(1+l->hi.t);
        if(left>0 && crossing-pad>*lo)*lo=crossing-pad;
        if(right>0 && crossing+pad<*hi)*hi=crossing+pad;
    }
    return 2;
}

static bool clip_program(interval_level_t* level,uint32_t code,double* lo,double* hi) {
    const interval_index_t* idx=level->query->index;
    if((code&(GUARD_BRANCH|QUERY_PROGRAM))!=QUERY_PROGRAM)
        return code&GUARD_BRANCH || classify_ray_atom(level,code,lo,hi,true)!=0;
    const uint32_t* pc=idx->decisions.data+(code&~QUERY_PROGRAM);
    for(unsigned step=0;step<18;step++) {
        query_decision_t ins=unpack_decision(idx,*pc);
        int value=classify_ray_atom(level,ins.operand,lo,hi,ins.no==0);
        if(value==2) {
            if(ins.no!=0 && ins.yes!=0)return true;
            value=ins.no==0;
        }
        unsigned jump=value?ins.yes:ins.no;
        if(jump==0 || jump==QUERY_ACCEPT)return jump!=0;
        pc-=jump;
    }
    return true;
}
