// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
// SPDX-License-Identifier: MPL-2.0
/* Plane-bounds preparation benchmark. Optional baseline is linked separately. */
#include "alea.h"
#include "core/alea_system.h"
#include "primitives/bbox.h"
#include "util/compat.h"
#include <stdio.h>
#include <math.h>
#ifdef ALEA_BBOX_BASELINE
int baseline_alea_tighten_bbox_plane_constraints(const alea_system_t*,alea_node_id_t,double,alea_bbox_t*);
#endif
int main(void){
    int sizes[]={6,32,64,128,256,1024};
    for(size_t k=0;k<sizeof(sizes)/sizeof(*sizes);k++){
        alea_system_t* sys=alea_create();alea_node_id_t root=ALEA_NODE_ID_INVALID;
        for(int i=0;i<sizes[k];i++){
            double a=0,b=0,c=0,d=-10;
            if(i<6){if(i/2==0)a=i%2?1:-1;if(i/2==1)b=i%2?1:-1;if(i/2==2)c=i%2?1:-1;}
            else{a=1+i%7;b=1+i%11;c=1+i%13;d=-10*(a+b+c)-1;}
            int sid=alea_plane_surface(sys,0,a,b,c,d);alea_node_id_t n=alea_halfspace(sys,sid,-1);root=i?alea_intersection(sys,root,n):n;
        }
        alea_bbox_t out;double t=alea_monotonic_seconds();int status=0;
        for(int i=0;i<20;i++)status=alea_tighten_bbox_plane_constraints(sys,root,1e-8,&out);
        printf("planes=%d lp_ms=%.6f status=%d",sizes[k],(alea_monotonic_seconds()-t)*50,status);
#ifdef ALEA_BBOX_BASELINE
        t=alea_monotonic_seconds();for(int i=0;i<20;i++)status=baseline_alea_tighten_bbox_plane_constraints(sys,root,1e-8,&out);
        printf(" baseline_ms=%.6f baseline_status=%d",(alea_monotonic_seconds()-t)*50,status);
#endif
        puts("");alea_destroy(sys);
    }
    /* Cylinder interior intersected with z >= x and z <= x+3. Analytic
     * enclosure [-2,2] x [-2,2] x [-2,5], despite unbounded plane-only LP. */
    alea_system_t* sys=alea_create();
    int cylinder=alea_cylinder_z_surface(sys,0,0,0,2);
    int lower=alea_plane_surface(sys,0,-1,0,1,0);
    int upper=alea_plane_surface(sys,0,-1,0,1,-3);
    alea_node_id_t root=alea_intersection(sys,alea_halfspace(sys,cylinder,-1),
        alea_intersection(sys,alea_halfspace(sys,lower,1),alea_halfspace(sys,upper,-1)));
    alea_bbox_t out;
    double start=alea_monotonic_seconds();
    int status=0;
    for(int i=0;i<100;i++)status=alea_tighten_bbox_mixed_constraints(sys,root,&out);
    if(status || out.min_x> -2 || out.max_x<2 || out.min_z> -2 || out.max_z<5) {
        alea_destroy(sys);return 1;
    }
    printf("mixed=cylinder-oblique-caps lp_ms=%.6f status=%d volume=%.9g\n",
           (alea_monotonic_seconds()-start)*10,status,alea_bbox_volume(&out));
    alea_destroy(sys);
    return 0;
}
