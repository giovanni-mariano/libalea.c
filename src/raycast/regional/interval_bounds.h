// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
// SPDX-License-Identifier: MPL-2.0
/* Conservative bounds within the prepared universe domain. */
static bool box_empty(interval_box_t b) {
    return b.lo[0]>b.hi[0] || b.lo[1]>b.hi[1] || b.lo[2]>b.hi[2];
}

static interval_box_t plane_contract(const alea_plane_data_t* p, double sign, interval_box_t b) {
    double a[] = {sign*p->a, sign*p->b, sign*p->c};
    double scale = fabs(p->d), terms[3];
    for (unsigned i=0;i<3;i++) {
        terms[i] = a[i] * (a[i] >= 0 ? b.lo[i] : b.hi[i]);
        scale += fabs(a[i])*fmax(fabs(b.lo[i]),fabs(b.hi[i]));
    }
    if (!isfinite(scale)) return b;
    double slack = 64 * (scale * DBL_EPSILON + DBL_MIN);
    for (unsigned i=0;i<3;i++) {
        double rest = sign*p->d + terms[(i+1)%3] + terms[(i+2)%3];
        if (a[i]>0) b.hi[i] = fmin(b.hi[i], nextafter((slack-rest)/a[i], INFINITY));
        if (a[i]<0) b.lo[i] = fmax(b.lo[i], nextafter((slack-rest)/a[i], -INFINITY));
    }
    return b;
}

/* Minimize a one-variable quadratic on a closed interval. Callers use
 * a scale-based FP64 margin for coefficient shifts and all evaluations. */
static double min_quadratic(double a,double b,double c,double lo,double hi) {
    double v=fmin((a*lo+b)*lo+c,(a*hi+b)*hi+c);
    if (a>0 && 2*a*lo+b<0 && 2*a*hi+b>0) {
        double x=-0.5*b/a;
        v=fmin(v,(a*x+b)*x+c);
    }
    return v;
}

static interval_box_t quad_contract(const double* k,double sign,interval_box_t box) {
    double c[3],h[3],m[3],grad[3],scale=fabs(k[9]);
    for (unsigned i=0;i<3;i++) {
        c[i]=0.5*box.lo[i]+0.5*box.hi[i];
        h[i]=fmax(c[i]-box.lo[i],box.hi[i]-c[i]);
        m[i]=fmax(1,fmax(fabs(box.lo[i]),fabs(box.hi[i])));
    }
    double v=sign*(k[0]*c[0]*c[0]+k[1]*c[1]*c[1]+k[2]*c[2]*c[2]+
        k[3]*c[0]*c[1]+k[4]*c[1]*c[2]+k[5]*c[2]*c[0]+k[6]*c[0]+k[7]*c[1]+k[8]*c[2]+k[9]);
    for (unsigned i=0;i<3;i++) {
        unsigned j=(i+1)%3,b=(i+2)%3;
        grad[i]=sign*(2*k[i]*c[i]+k[i+3]*c[j]+k[b+3]*c[b]+k[i+6]);
        scale+=fabs(k[i])*m[i]*m[i]+fabs(k[i+3])*m[i]*m[j]+fabs(k[i+6])*m[i];
    }
    double eps=1024*(DBL_EPSILON*scale+DBL_MIN);
    if (!isfinite(eps)) return box;
    for (unsigned i=0;i<3;i++) {
        unsigned j=(i+1)%3,b=(i+2)%3;
        double a=sign*k[i];
        double beta=fabs(k[i+3])*h[j]+fabs(k[b+3])*h[b];
        double gamma=v-eps+min_quadratic(sign*k[j],grad[j],0,-h[j],h[j])
                          +min_quadratic(sign*k[b],grad[b],0,-h[b],h[b])-fabs(k[j+3])*h[j]*h[b];
        double low=INFINITY,high=-INFINITY;
        for (unsigned side=0;side<2;side++) {
            double left=side?0:-h[i],right=side?h[i]:0, slope=grad[i]+(side?-beta:beta);
            if (min_quadratic(a,slope,gamma,left,right)>0) continue;
            double l=left,r=right;
            if ((a*l+slope)*l+gamma>0) {
                double x=l,y=r;
                for (unsigned step=0;step<8;step++) {
                    double mid=0.5*(x+y);
                    if (min_quadratic(a,slope,gamma,l,mid)>0) x=mid; else y=mid;
                }
                left=x;
            }
            if ((a*r+slope)*r+gamma>0) {
                double x=l,y=r;
                for (unsigned step=0;step<8;step++) {
                    double mid=0.5*(x+y);
                    if (min_quadratic(a,slope,gamma,mid,r)>0) y=mid; else x=mid;
                }
                right=y;
            }
            low=fmin(low,left);high=fmax(high,right);
        }
        if (low>high) {box.lo[i]=INFINITY;return box;}
        double pad=32*DBL_EPSILON*m[i];
        box.lo[i]=fmax(box.lo[i],c[i]+low-pad);
        box.hi[i]=fmin(box.hi[i],c[i]+high+pad);
    }
    return box;
}

/* For a sphere or circular cylinder the squared distances of the other
 * active coordinates give a lower bound when contracting one coordinate. */
static interval_box_t radial_contract(const double center[3],unsigned axis,double radius,
                                      bool outside,interval_box_t box) {
    double nearest[3]={0},near2=0,far2=0,r2=radius*radius,scale=r2;
    for(unsigned i=0;i<3;i++) {
        if(i==axis)continue;
        double lo=box.lo[i]-center[i], hi=box.hi[i]-center[i];
        double pad=32*DBL_EPSILON*(1+fabs(center[i])+fabs(box.lo[i])+fabs(box.hi[i]));
        double near=fmax(0,fmax(lo,-hi)-pad);
        double far=fmax(fabs(lo),fabs(hi))+pad;
        nearest[i]=near*near;
        near2+=nearest[i];far2+=far*far;
        scale+=far*far;
    }
    double slack=128*(DBL_EPSILON*scale+DBL_MIN);
    if(!isfinite(slack))return box;
    if(outside ? far2<r2-slack : near2>r2+slack) {
        box.lo[0]=INFINITY;return box;
    }
    if(outside)return box;
    for(unsigned i=0;i<3;i++) {
        if(i==axis)continue;
        double reach=sqrt(fmax(0,r2+slack-nearest[(i+1)%3]-nearest[(i+2)%3]));
        double pad=32*DBL_EPSILON*(1+fabs(center[i])+reach);
        box.lo[i]=fmax(box.lo[i],center[i]-reach-pad);
        box.hi[i]=fmin(box.hi[i],center[i]+reach+pad);
    }
    return box;
}

static interval_box_t contract_cell(const alea_system_t* sys, uint32_t id, interval_box_t b,
                                    bool positive, unsigned depth, unsigned* budget) {
    if (depth>=192 || !*budget || id>=sys->nodes.count) return b;
    --*budget;
    const alea_node_t* n = &sys->nodes.data[id];
    interval_box_t node_box = interval_box_from(alea_node_bbox_get(&n->bbox));
    if (positive) {
        /* Outside this node's box the scalar evaluator returns +1. */
        for (unsigned i=0;i<3;i++) {
            if (b.lo[i]<node_box.lo[i] || b.hi[i]>node_box.hi[i]) return b;
        }
    } else {
        /* Most leaf boxes already fit. A strict comparison skips only a
         * known identity. Ties (including signed zero) and unordered
         * bounds retain the original floating-point min/max behavior. */
        for (unsigned i=0;i<3;i++) {
            if (!(node_box.lo[i] < b.lo[i]))
                b.lo[i]=fmax(b.lo[i],node_box.lo[i]);
            if (!(node_box.hi[i] > b.hi[i]))
                b.hi[i]=fmin(b.hi[i],node_box.hi[i]);
        }
    }
    if (box_empty(b)) return b;
    alea_operation_t op=ALEA_GET_OPERATION(n);
    if (op==ALEA_OP_PRIMITIVE) {
        if (n->primitive.prim_type==ALEA_PRIMITIVE_PLANE) {
            const alea_plane_data_t* p=alea_primitive_payload_const(sys,n->primitive.primitive_id);
            bool flip=(n->primitive.sense>0) != (n->primitive.inverted!=0);
            if (p) return plane_contract(p,flip!=positive ? -1.0 : 1.0,b);
        }
        if (n->primitive.prim_type==ALEA_PRIMITIVE_QUADRIC) {
            const alea_quadric_data_t* q=alea_primitive_payload_const(sys,n->primitive.primitive_id);
            bool flip=(n->primitive.sense>0) != (n->primitive.inverted!=0);
            if (q) return quad_contract(q->coeffs,flip!=positive ? -1.0 : 1.0,b);
        }
        unsigned axis=3;
        double center[3]={0},radius;
        const void* payload=alea_primitive_payload_const(sys,n->primitive.primitive_id);
        if(!payload)return b;
        switch(n->primitive.prim_type) {
        case ALEA_PRIMITIVE_SPHERE: case ALEA_PRIMITIVE_SPH: {
            const alea_sphere_data_t* s=payload;
            center[0]=s->center_x;center[1]=s->center_y;center[2]=s->center_z;radius=s->radius;break;
        }
        case ALEA_PRIMITIVE_CYLINDER_X: {
            const alea_cylinder_x_data_t* c=payload;
            center[1]=c->center_y;center[2]=c->center_z;radius=c->radius;axis=0;break;
        }
        case ALEA_PRIMITIVE_CYLINDER_Y: {
            const alea_cylinder_y_data_t* c=payload;
            center[0]=c->center_x;center[2]=c->center_z;radius=c->radius;axis=1;break;
        }
        case ALEA_PRIMITIVE_CYLINDER_Z: {
            const alea_cylinder_z_data_t* c=payload;
            center[0]=c->center_x;center[1]=c->center_y;radius=c->radius;axis=2;break;
        }
        default:return b;
        }
        bool flip=(n->primitive.sense>0)!=(n->primitive.inverted!=0);
        return radial_contract(center,axis,radius,flip!=positive,b);
    }
    if (op==ALEA_OP_COMPLEMENT)
        return contract_cell(sys,n->operation.left,b,!positive,depth+1,budget);
    /* Union returns its left operand immediately on zero. Its positive
     * side can therefore contain a zero even if the right side is negative. */
    if (positive && op==ALEA_OP_UNION)
        return contract_cell(sys,n->operation.left,b,true,depth+1,budget);
    bool right_positive = positive != (op==ALEA_OP_DIFFERENCE);
    interval_box_t left=contract_cell(sys,n->operation.left,b,positive,depth+1,budget);
    bool join=(op==ALEA_OP_UNION) != positive;
    if (!join && box_empty(left)) return left;
    interval_box_t right=contract_cell(sys,n->operation.right,join?b:left,right_positive,depth+1,budget);
    if (!join || box_empty(left)) return right;
    if (box_empty(right)) return left;
    for (unsigned i=0;i<3;i++) {
        left.lo[i]=fmin(left.lo[i],right.lo[i]);
        left.hi[i]=fmax(left.hi[i],right.hi[i]);
    }
    return left;
}

static interval_box_t cell_domain(const alea_system_t* sys, uint32_t id, alea_bbox_t domain) {
    interval_box_t b=interval_box_from(domain);
    const alea_cell_entry_t* c=&sys->cells.data[id];
    if (c->lat_type && c->lat_fill) return b;
    unsigned budget=32768;
    for (unsigned pass=0;pass<12 && budget;pass++) {
        interval_box_t before=b;
        b=contract_cell(sys,c->root_node_id,b,false,0,&budget);
        if (box_empty(b) || !memcmp(&before,&b,sizeof(b))) break;
    }
    return b;
}

static bool overlaps_box(interval_box_t b, const alea_bbox_t* r) {
    return b.lo[0]<=r->max_x && b.hi[0]>=r->min_x &&
           b.lo[1]<=r->max_y && b.hi[1]>=r->min_y &&
           b.lo[2]<=r->max_z && b.hi[2]>=r->min_z;
}
