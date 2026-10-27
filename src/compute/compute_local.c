// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
// SPDX-License-Identifier: MPL-2.0
#include "compute_internal.h"
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <math.h>

/* Bounds enclose FP32 arithmetic, including rounded intermediate results.
 * Unknown/overflow stays unknown. Strict signs preserve closed boundaries and
 * complements; no sampled point is used to prove a spatial predicate. */
typedef struct { double lo,hi; } interval;
static interval unknown(void) { return (interval){-INFINITY,INFINITY}; }
static interval outward(double lo,double hi) {
    if(!isfinite(lo) || !isfinite(hi) || lo < -FLT_MAX || hi > FLT_MAX)return unknown();
    float lower=nextafterf((float)lo,-INFINITY),upper=nextafterf((float)hi,INFINITY);
    /* Metal may flush subnormal intermediate results to zero. */
    if(fabsf(lower)<FLT_MIN)lower=fminf(lower,0);
    if(fabsf(upper)<FLT_MIN)upper=fmaxf(upper,0);
    return (interval){lower,upper};
}
static interval scalar(double x) { return (interval){x,x}; }
static interval add(interval a,interval b) { return outward(a.lo+b.lo,a.hi+b.hi); }
static interval neg(interval a) { return (interval){-a.hi,-a.lo}; }
static interval mul(interval a,interval b) {
    double v[]={a.lo*b.lo,a.lo*b.hi,a.hi*b.lo,a.hi*b.hi};
    double lo=v[0],hi=v[0];
    for(int i=0;i<4;++i){if(!isfinite(v[i]))return unknown();lo=fmin(lo,v[i]);hi=fmax(hi,v[i]);}
    return outward(lo,hi);
}
static interval square(interval a) {
    if(!isfinite(a.lo) || !isfinite(a.hi))return unknown();
    double x=a.lo*a.lo,y=a.hi*a.hi;
    return outward(a.lo<=0 && a.hi>=0?0:fmin(x,y),fmax(x,y));
}
static interval root(interval a) {
    if(!isfinite(a.lo) || !isfinite(a.hi) || a.hi<0)return unknown();
    interval r=outward(sqrt(fmax(0,a.lo)),sqrt(fmax(0,a.hi)));
    /* Reserve extra ULPs for the device length/sqrt implementation. */
    r.lo=nextafterf(nextafterf((float)r.lo,-INFINITY),-INFINITY);
    r.hi=nextafterf(nextafterf((float)r.hi,INFINITY),INFINITY);return r;
}
static interval maximum(interval a,interval b) {return (interval){fmax(a.lo,b.lo),fmax(a.hi,b.hi)};}
int alea_compute_primitive_sign(const compute_primitive* primitive,const float b[6]) {
    const float* a=primitive->p;interval p[3];
    for(int i=0;i<13;++i)if(a[i]!=0 && fabsf(a[i])<FLT_MIN)return 0;
    for(int i=0;i<6;++i)if(b[i]!=0 && fabsf(b[i])<FLT_MIN)return 0;
    for(int i=0;i<3;++i)p[i]=(interval){b[2*i],b[2*i+1]};
    interval v=unknown();
    switch(primitive->type) {
    case COMPUTE_PLANE:
        v=add(add(add(mul(scalar(a[0]),p[0]),mul(scalar(a[1]),p[1])),mul(scalar(a[2]),p[2])),scalar(a[3]));break;
    case COMPUTE_SPHERE:
        v=scalar(0);
        for(int i=0;i<3;++i)v=add(v,square(add(p[i],scalar(-a[i]))));
        v=add(v,neg(square(scalar(a[3]))));break;
    case COMPUTE_CYL_X:case COMPUTE_CYL_Y:case COMPUTE_CYL_Z: {
        int axis=primitive->type-COMPUTE_CYL_X,j=0;v=scalar(0);
        for(int i=0;i<3;++i)if(i!=axis)v=add(v,square(add(p[i],scalar(-a[j++]))));
        v=add(v,neg(square(scalar(a[2]))));break;
    }
    case COMPUTE_BOX:
        v=scalar(-INFINITY);
        for(int i=0;i<3;++i)v=maximum(v,maximum(add(scalar(a[2*i]),neg(p[i])),add(p[i],scalar(-a[2*i+1]))));
        break;
    case COMPUTE_QUADRIC: {
        for(int i=0;i<3;++i)p[i]=add(p[i],scalar(-a[10+i]));
        interval x=add(add(add(mul(scalar(a[0]),p[0]),mul(scalar(a[3]),p[1])),mul(scalar(a[5]),p[2])),scalar(a[6]));
        interval y=add(add(mul(scalar(a[1]),p[1]),mul(scalar(a[4]),p[2])),scalar(a[7]));
        interval z=add(mul(scalar(a[2]),p[2]),scalar(a[8]));
        v=add(add(add(mul(p[0],x),mul(p[1],y)),mul(p[2],z)),scalar(a[9]));break;
    }
    case COMPUTE_CONE_X:case COMPUTE_CONE_Y:case COMPUTE_CONE_Z: {
        int axis=primitive->type-COMPUTE_CONE_X;
        for(int i=0;i<3;++i)p[i]=add(p[i],scalar(-a[i]));
        interval sheet=a[4]==0?scalar(0):mul(scalar(a[4]),p[axis]);
        if(sheet.hi<0){v=scalar(1);break;}
        v=scalar(0);for(int i=0;i<3;++i)if(i!=axis)v=add(v,square(p[i]));
        v=add(v,neg(mul(mul(scalar(a[3]),p[axis]),p[axis])));
        if(sheet.lo<0){v.lo=fmin(v.lo,1);v.hi=fmax(v.hi,1);}break;
    }
    case COMPUTE_TORUS_X:case COMPUTE_TORUS_Y:case COMPUTE_TORUS_Z: {
        if(a[5]<=0)break;
        int axis=primitive->type-COMPUTE_TORUS_X;
        for(int i=0;i<3;++i)p[i]=add(p[i],scalar(-a[i]));
        interval radial=scalar(0);for(int i=0;i<3;++i)if(i!=axis)radial=add(radial,square(p[i]));
        radial=add(root(radial),scalar(-a[3]));
        interval axial=mul(p[axis],scalar(a[4]));
        axial=outward(axial.lo/a[5],axial.hi/a[5]);
        v=add(add(square(radial),square(axial)),neg(square(scalar(a[4]))));break;
    }
    case COMPUTE_RCC: {
        for(int i=0;i<3;++i)p[i]=add(p[i],scalar(-a[i]));
        interval h=scalar(0);for(int i=0;i<3;++i)h=add(h,mul(p[i],scalar(a[3+i])));
        v=scalar(0);for(int i=0;i<3;++i)v=add(v,square(add(p[i],neg(mul(h,scalar(a[3+i]))))));
        v=maximum(add(v,neg(square(scalar(a[7])))),maximum(neg(h),add(h,scalar(-a[6]))));break;
    }
    default:break;
    }
    return v.lo>0?1:v.hi<0?-1:0;
}

typedef struct { uint32_t op,arg,left,right,count; int sign; } expr;
static expr constant(int sign) { return (expr){COMPUTE_CONSTANT,sign>0,0,0,1,sign}; }
static expr complement(const expr* nodes,uint32_t child) {
    if(nodes[child].sign)return constant(-nodes[child].sign);
    if(nodes[child].op==COMPUTE_COMPLEMENT)return nodes[nodes[child].left];
    return (expr){COMPUTE_COMPLEMENT,0,child,0,nodes[child].count+1,0};
}
static int emit(const expr* nodes,uint32_t id,size_t count,compute_node* out,size_t* at) {
    /* Parent masks can combine several deep source programs. Iterative emission
     * keeps host stack use independent of the combined expression depth. */
    typedef struct {uint32_t id,phase;} frame;
    frame* stack=malloc(count*sizeof(*stack));
    if(!stack){alea_set_error_detail(ALEA_ERR_OUT_OF_MEMORY,"Emitting local CSG program");return -1;}
    size_t depth=1;stack[0]=(frame){id,0};
    while(depth) {
        frame* f=&stack[depth-1];expr n=nodes[f->id];
        if(n.op==COMPUTE_LEAF || n.op==COMPUTE_CONSTANT) {
            out[(*at)++]=(compute_node){n.op,n.arg};--depth;
        } else if(f->phase==0) {
            f->phase=1;stack[depth++]=(frame){n.left,0};
        } else if(f->phase==1 && n.op!=COMPUTE_COMPLEMENT) {
            f->phase=2;
            uint32_t skip=n.op==COMPUTE_UNION?COMPUTE_SKIP_UNION:
                n.op==COMPUTE_REVERSE_DIFFERENCE?COMPUTE_SKIP_REVERSE_DIFFERENCE:COMPUTE_SKIP_INTERSECTION;
            out[(*at)++]=(compute_node){skip,nodes[n.right].count+1};
            stack[depth++]=(frame){n.right,0};
        } else {out[(*at)++]=(compute_node){n.op,0};--depth;}
    }
    free(stack);return 0;
}
int alea_compute_local_simplify(const compute_node* code,size_t count,
    const compute_primitive* primitives,const float bounds[6],
    compute_node* output,size_t* output_count) {
    /* The source packer validates primitive indices and jump structure. Keep
     * local stack checks as well, so a compiler regression fails explicitly. */
    if(!code || !output || !output_count || !primitives || !bounds ||
        !count || count>UINT32_MAX || count>SIZE_MAX/sizeof(expr)) {
        alea_set_error_detail(ALEA_ERR_INVALID_STATE,"Invalid local CSG program size");return -1;
    }
    expr* nodes=calloc(count,sizeof(*nodes));
    if(!nodes){alea_set_error_detail(ALEA_ERR_OUT_OF_MEMORY,"Reducing local CSG program");return -1;}
    uint32_t stack[COMPUTE_STACK]={0};size_t sp=0;
    for(uint32_t i=0;i<count;++i) {
        compute_node n=code[i];expr result={n.op,n.arg,0,0,1,0};
        if(n.op==COMPUTE_SKIP_UNION || n.op==COMPUTE_SKIP_INTERSECTION || n.op==COMPUTE_SKIP_REVERSE_DIFFERENCE)continue;
        if(n.op==COMPUTE_LEAF) {
            int sign=alea_compute_primitive_sign(&primitives[n.arg>>1],bounds);if(n.arg&1)sign=-sign;
            if(sign)result=constant(sign);
        } else if(n.op==COMPUTE_CONSTANT)result=constant(n.arg?1:-1);
        else if(n.op==COMPUTE_COMPLEMENT) {
            if(!sp)goto invalid;
            result=complement(nodes,stack[--sp]);
        }
        else {
            if(sp<2 || (n.op!=COMPUTE_UNION && n.op!=COMPUTE_INTERSECTION &&
                n.op!=COMPUTE_DIFFERENCE && n.op!=COMPUTE_REVERSE_DIFFERENCE))goto invalid;
            uint32_t right=stack[--sp],left=stack[--sp];
            int a=nodes[left].sign,b=nodes[right].sign;
            result=(expr){n.op,0,left,right,nodes[left].count+nodes[right].count+2,0};
            if(n.op==COMPUTE_UNION) {
                if(a<0 || b<0)result=constant(-1);
                else if(a>0)result=nodes[right];else if(b>0)result=nodes[left];
            } else if(n.op==COMPUTE_INTERSECTION) {
                if(a>0 || b>0)result=constant(1);
                else if(a<0)result=nodes[right];else if(b<0)result=nodes[left];
            } else {
                if(n.op==COMPUTE_REVERSE_DIFFERENCE){uint32_t t=left;left=right;right=t;a=nodes[left].sign;b=nodes[right].sign;}
                if(a>0 || b<0)result=constant(1);
                else if(b>0)result=nodes[left];else if(a<0)result=complement(nodes,right);
            }
        }
        if(sp>=COMPUTE_STACK)goto invalid;
        nodes[i]=result;stack[sp++]=i;
    }
    if(sp!=1)goto invalid;
    size_t at=0;int rc=emit(nodes,stack[0],count,output,&at);*output_count=at;free(nodes);return rc;
invalid:
    free(nodes);alea_set_error_detail(ALEA_ERR_INVALID_STATE,"Invalid local CSG value stack");return -1;
}
