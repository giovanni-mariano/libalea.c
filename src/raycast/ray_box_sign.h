// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
// SPDX-License-Identifier: MPL-2.0
#ifndef ALEA_RAY_BOX_SIGN_H
#define ALEA_RAY_BOX_SIGN_H
#include "core/alea_system.h"
#include <float.h>
#include <math.h>
#include <string.h>

typedef struct { double lo[3], hi[3]; } interval_box_t;

static inline interval_box_t interval_box_from(alea_bbox_t b) {
    return (interval_box_t){{b.min_x,b.min_y,b.min_z},{b.max_x,b.max_y,b.max_z}};
}


typedef struct { double k[10]; int kind; unsigned axis; } guard_surface_t;
static inline double guard_value(const double* k,double x,double y,double z) {
    return k[0]*x*x+k[1]*y*y+k[2]*z*z+k[3]*x*y+k[4]*y*z+k[5]*x*z+
           k[6]*x+k[7]*y+k[8]*z+k[9];
}



/* Polynomial in world or local coordinates Ax²+By²+Cz²+Dxy+Eyz+Fxz+
 * Gx+Hy+Iz+J. Unsupported primitives retain an unrestricted predicate. */
enum { GUARD_RADIAL=3, GUARD_TORUS=4 };
static inline void guard_surface_init(guard_surface_t* out, const alea_system_t* s, uint32_t id, bool ring_only) {
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
        if ((ring_only && v->major_radius <= v->minor_radius) ||
           !(v->major_radius >= 0 && v->minor_radius > 0 && v->axial_semiwidth_B > 0) ||
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
static inline int radial_box_sign(const guard_surface_t* s,interval_box_t box) {
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
static inline int guard_box_sign(const guard_surface_t* s,interval_box_t box) {
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

#endif
