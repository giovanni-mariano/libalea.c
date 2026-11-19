// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
// SPDX-License-Identifier: MPL-2.0
#include "primitive_codec.h"
#include <string.h>
#include <math.h>
int alea_primitive_from_values(alea_primitive_type_t type, const double* v,
                            alea_primitive_data_t* d) {
    double scratch[50];
    const alea_primitive_data_t empty={0};
    size_t n=alea_primitive_to_values(type,&empty,scratch);
    if(!n)return -1;
    for(size_t i=0;i<n;i++)if(!isfinite(v[i]))return -1;
    memset(d, 0, sizeof(*d));
    switch (type) {
        case ALEA_PRIMITIVE_PLANE: memcpy(&d->plane, v, 4 * sizeof(double)); break;
        case ALEA_PRIMITIVE_SPHERE: memcpy(&d->sphere, v, 4 * sizeof(double)); break;
        case ALEA_PRIMITIVE_CYLINDER_X: memcpy(&d->cyl_x, v, 3 * sizeof(double)); break;
        case ALEA_PRIMITIVE_CYLINDER_Y: memcpy(&d->cyl_y, v, 3 * sizeof(double)); break;
        case ALEA_PRIMITIVE_CYLINDER_Z: memcpy(&d->cyl_z, v, 3 * sizeof(double)); break;
        case ALEA_PRIMITIVE_CONE_X:
            if (v[4] != trunc(v[4]) || v[4] < -1 || v[4] > 1) return -1;
            memcpy(&d->cone_x, v, 4 * sizeof(double)); d->cone_x.sheet_selection = (int)v[4]; break;
        case ALEA_PRIMITIVE_CONE_Y:
            if (v[4] != trunc(v[4]) || v[4] < -1 || v[4] > 1) return -1;
            memcpy(&d->cone_y, v, 4 * sizeof(double)); d->cone_y.sheet_selection = (int)v[4]; break;
        case ALEA_PRIMITIVE_CONE_Z:
            if (v[4] != trunc(v[4]) || v[4] < -1 || v[4] > 1) return -1;
            memcpy(&d->cone_z, v, 4 * sizeof(double)); d->cone_z.sheet_selection = (int)v[4]; break;
        case ALEA_PRIMITIVE_RPP: memcpy(&d->box, v, 6 * sizeof(double)); break;
        case ALEA_PRIMITIVE_QUADRIC: memcpy(d->quadric.coeffs, v, 10 * sizeof(double)); break;
        case ALEA_PRIMITIVE_TORUS_X: case ALEA_PRIMITIVE_TORUS_Y: case ALEA_PRIMITIVE_TORUS_Z:
            d->torus.axis = type == ALEA_PRIMITIVE_TORUS_X ? ALEA_AXIS_X :
                            type == ALEA_PRIMITIVE_TORUS_Y ? ALEA_AXIS_Y : ALEA_AXIS_Z;
            d->torus.center_x=v[0]; d->torus.center_y=v[1]; d->torus.center_z=v[2];
            d->torus.major_radius=v[3]; d->torus.minor_radius=v[4]; d->torus.axial_semiwidth_B=v[5]; break;
        case ALEA_PRIMITIVE_RCC: memcpy(&d->rcc, v, 7 * sizeof(double)); break;
        case ALEA_PRIMITIVE_BOX: memcpy(&d->box_general, v, 12 * sizeof(double)); break;
        case ALEA_PRIMITIVE_SPH: memcpy(&d->sph, v, 4 * sizeof(double)); break;
        case ALEA_PRIMITIVE_TRC: memcpy(&d->trc, v, 8 * sizeof(double)); break;
        case ALEA_PRIMITIVE_ELL: memcpy(&d->ell, v, 7 * sizeof(double)); break;
        case ALEA_PRIMITIVE_REC: memcpy(&d->rec, v, 12 * sizeof(double)); break;
        case ALEA_PRIMITIVE_WED: memcpy(&d->wed, v, 12 * sizeof(double)); break;
        case ALEA_PRIMITIVE_RHP: memcpy(&d->rhp, v, 15 * sizeof(double)); break;
        case ALEA_PRIMITIVE_ARB:
            if (v[0] != trunc(v[0]) || v[1] != trunc(v[1])) return -1;
            d->arb.num_corners=(int)v[0]; d->arb.num_faces=(int)v[1];
            if (d->arb.num_corners < 4 || d->arb.num_corners > 8 ||
                d->arb.num_faces < 4 || d->arb.num_faces > 6) return -1;
            memcpy(d->arb.corners, v + 2, 24 * sizeof(double));
            for (int i=0; i<24; i++) {
                if (v[26+i] != trunc(v[26+i]) || v[26+i] < 0 || v[26+i] > 8) return -1;
                ((int*)d->arb.faces)[i]=(int)v[26+i];
            }
            break;
        default: return -1;
    }
    return 0;
}

size_t alea_primitive_to_values(alea_primitive_type_t type,
                               const alea_primitive_data_t* d, double* v) {
    size_t n = 0;
    switch (type) {
        case ALEA_PRIMITIVE_PLANE: n=4; memcpy(v,&d->plane,n*sizeof(double)); break;
        case ALEA_PRIMITIVE_SPHERE: n=4; memcpy(v,&d->sphere,n*sizeof(double)); break;
        case ALEA_PRIMITIVE_CYLINDER_X: n=3; memcpy(v,&d->cyl_x,n*sizeof(double)); break;
        case ALEA_PRIMITIVE_CYLINDER_Y: n=3; memcpy(v,&d->cyl_y,n*sizeof(double)); break;
        case ALEA_PRIMITIVE_CYLINDER_Z: n=3; memcpy(v,&d->cyl_z,n*sizeof(double)); break;
        case ALEA_PRIMITIVE_CONE_X: n=5; memcpy(v,&d->cone_x,4*sizeof(double)); v[4]=d->cone_x.sheet_selection; break;
        case ALEA_PRIMITIVE_CONE_Y: n=5; memcpy(v,&d->cone_y,4*sizeof(double)); v[4]=d->cone_y.sheet_selection; break;
        case ALEA_PRIMITIVE_CONE_Z: n=5; memcpy(v,&d->cone_z,4*sizeof(double)); v[4]=d->cone_z.sheet_selection; break;
        case ALEA_PRIMITIVE_RPP: n=6; memcpy(v,&d->box,n*sizeof(double)); break;
        case ALEA_PRIMITIVE_QUADRIC: n=10; memcpy(v,d->quadric.coeffs,n*sizeof(double)); break;
        case ALEA_PRIMITIVE_TORUS_X: case ALEA_PRIMITIVE_TORUS_Y: case ALEA_PRIMITIVE_TORUS_Z:
            n=6; v[0]=d->torus.center_x; v[1]=d->torus.center_y; v[2]=d->torus.center_z;
            v[3]=d->torus.major_radius; v[4]=d->torus.minor_radius; v[5]=d->torus.axial_semiwidth_B; break;
        case ALEA_PRIMITIVE_RCC: n=7; memcpy(v,&d->rcc,n*sizeof(double)); break;
        case ALEA_PRIMITIVE_BOX: n=12; memcpy(v,&d->box_general,n*sizeof(double)); break;
        case ALEA_PRIMITIVE_SPH: n=4; memcpy(v,&d->sph,n*sizeof(double)); break;
        case ALEA_PRIMITIVE_TRC: n=8; memcpy(v,&d->trc,n*sizeof(double)); break;
        case ALEA_PRIMITIVE_ELL: n=7; memcpy(v,&d->ell,n*sizeof(double)); break;
        case ALEA_PRIMITIVE_REC: n=12; memcpy(v,&d->rec,n*sizeof(double)); break;
        case ALEA_PRIMITIVE_WED: n=12; memcpy(v,&d->wed,n*sizeof(double)); break;
        case ALEA_PRIMITIVE_RHP: n=15; memcpy(v,&d->rhp,n*sizeof(double)); break;
        case ALEA_PRIMITIVE_ARB:
            n=50; v[0]=d->arb.num_corners; v[1]=d->arb.num_faces;
            memcpy(v+2,d->arb.corners,24*sizeof(double));
            for (int i=0;i<24;i++) v[26+i]=((const int*)d->arb.faces)[i];
            break;
        default: break;
    }
    return n;
}
