// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
// SPDX-License-Identifier: MPL-2.0
/* Conservative bounds in the (x+y,x-y), (y+z,y-z) and (z+x,z-x) frames. */
static interval_box_t rotate_box(interval_box_t b,unsigned i) {
    unsigned j=(i+1)%3;
    interval_box_t q=b;
    q.lo[i]=b.lo[i]+b.lo[j];q.hi[i]=b.hi[i]+b.hi[j];
    q.lo[j]=b.lo[i]-b.hi[j];q.hi[j]=b.hi[i]-b.lo[j];
    return q;
}
static interval_box_t oriented_bound(const interval_index_t* idx,uint32_t code,
                                     interval_box_t b,unsigned i,unsigned* budget) {
    if(!*budget || code<=GUARD_ANY)return b;
    --*budget;
    if(code&GUARD_BRANCH) {
        guard_pair_t p=idx->guards.pairs.data[code&~(GUARD_BRANCH|GUARD_UNION)];
        interval_box_t l=oriented_bound(idx,p.a,b,i,budget);
        if(!(code&GUARD_UNION))return box_empty(l)?l:oriented_bound(idx,p.b,l,i,budget);
        interval_box_t r=oriented_bound(idx,p.b,b,i,budget);
        if(box_empty(l))return r;
        if(box_empty(r))return l;
        for(unsigned a=0;a<3;a++) {l.lo[a]=fmin(l.lo[a],r.lo[a]);l.hi[a]=fmax(l.hi[a],r.hi[a]);}
        return l;
    }
    if(code&(GUARD_REF|GUARD_BOX))return b;
    const guard_surface_t* s=&idx->surfaces[(code-2)>>1];
    if(!s->kind || s->kind>=GUARD_RADIAL)return b;
    unsigned j=(i+1)%3,k=(i+2)%3;
    const double* p=s->k;
    double q[10];memcpy(q,p,sizeof(q));
    q[6+i]=0.5*(p[6+i]+p[6+j]);q[6+j]=0.5*(p[6+i]-p[6+j]);
    double m=1;
    for(unsigned a=0;a<3;a++)m=fmax(m,fmax(fabs(b.lo[a]),fabs(b.hi[a])));
    double error=64*DBL_EPSILON*((fabs(p[6+i])+fabs(p[6+j]))*m+fabs(p[9])+1);
    double sign=code&1?-1:1;
    if(s->kind==1) {
        alea_plane_data_t plane={q[6],q[7],q[8],q[9]-sign*error};
        return plane_contract(&plane,sign,b);
    }
    q[i]=0.25*(p[i]+p[j]+p[3+i]);q[j]=0.25*(p[i]+p[j]-p[3+i]);
    q[3+i]=0.5*(p[i]-p[j]);
    q[3+k]=0.5*(p[3+k]+p[3+j]);q[3+j]=0.5*(p[3+k]-p[3+j]);
    for(unsigned a=0;a<6;a++)error+=64*DBL_EPSILON*fabs(p[a])*m*m;
    q[9]-=sign*error;
    return quad_contract(q,sign,b);
}
/* Each projection supplies a second enclosing box for the same point. */
static int sign_in_views(const guard_surface_t* s,interval_box_t b,const interval_box_t views[3]) {
    int sign=guard_box_sign(s,b);
    if(sign || !s->kind || s->kind>=GUARD_RADIAL)return sign;
    const double* p=s->k;
    for(unsigned i=0;i<3;i++) {
        unsigned j=(i+1)%3,k=(i+2)%3;
        const interval_box_t* v=views+i;
        guard_surface_t g=*s;
        double* q=g.k;
        q[6+i]=0.5*(p[6+i]+p[6+j]);q[6+j]=0.5*(p[6+i]-p[6+j]);
        double m=1;
        for(unsigned a=0;a<3;a++)m=fmax(m,fmax(fabs(v->lo[a]),fabs(v->hi[a])));
        double error=128*DBL_EPSILON*((fabs(p[6+i])+fabs(p[6+j]))*m+fabs(p[9])+1);
        if(s->kind==2) {
            q[i]=0.25*(p[i]+p[j]+p[3+i]);q[j]=0.25*(p[i]+p[j]-p[3+i]);
            q[3+i]=0.5*(p[i]-p[j]);
            q[3+k]=0.5*(p[3+k]+p[3+j]);q[3+j]=0.5*(p[3+k]-p[3+j]);
            for(unsigned a=0;a<6;a++)error+=128*DBL_EPSILON*fabs(p[a])*m*m;
        }
        q[9]=p[9]-error;
        if(guard_box_sign(&g,*v)>0)return 1;
        q[9]=p[9]+error;
        if(guard_box_sign(&g,*v)<0)return -1;
    }
    return 0;
}
static uint32_t fold_views(interval_index_t* idx,uint32_t code,interval_box_t b,const interval_box_t views[3],unsigned* budget) {
    if(!*budget || code<2)return code;
    --*budget;
    if(code&GUARD_BRANCH) {
        bool join=!!(code&GUARD_UNION);
        guard_pair_t p=idx->guards.pairs.data[code&~(GUARD_BRANCH|GUARD_UNION)];
        uint32_t a=fold_views(idx,p.a,b,views,budget);
        if(a==(join?GUARD_ANY:GUARD_EMPTY))return a;
        uint32_t right=fold_views(idx,p.b,b,views,budget);
        return guard_join(&idx->guards,a,right,join);
    }
    if(code&(GUARD_REF|GUARD_BOX))return code;
    guard_surface_t* s=idx->surfaces+((code-2)>>1);
    if(!s->kind)return code;
    int value=sign_in_views(s,b,views);
    return value==0?code:((value>0)==(bool)(code&1)?GUARD_ANY:GUARD_EMPTY);
}
