// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
// SPDX-License-Identifier: MPL-2.0
#include "core/alea_system.h"
static alea_system_t* ray_lattice_scene(void) {
    alea_system_t* sys=alea_create();if(!sys)return NULL;int mat=alea_add_material(sys,7);
    int envelope=alea_box_surface(sys,0,-2,4,-3,3,-3,3);
    int parent=alea_add_cell(sys,10,alea_halfspace(sys,envelope,-1),0,0,0);
    double translation[]={1,0,0};
    if(mat<0 || parent<0 || alea_add_transform(sys,1,translation,3,0) || alea_set_fill(sys,parent,1,1))goto fail;
    int unit=alea_box_surface(sys,0,-1,1,-1,1,-1,1);
    int element=alea_add_cell(sys,20,alea_halfspace(sys,unit,-1),0,0,1);if(element<0)goto fail;
    alea_cell_entry_t* c=&sys->cells.data[element];c->lat_type=1;c->lat_fill_count=27;
    c->lat_fill=malloc(27*sizeof(int));if(!c->lat_fill)goto fail;
    for(int a=0;a<3;++a){c->lat_fill_dims[2*a]=-1;c->lat_fill_dims[2*a+1]=1;c->lat_pitch[a]=2;c->lat_lower_left[a]=-3;}
    for(int i=0;i<27;++i)c->lat_fill[i]=2;
    int sphere=alea_sphere_surface(sys,0,0,0,0,.6);alea_node_id_t leaf=alea_halfspace(sys,sphere,-1);
    if(alea_add_cell(sys,42,alea_intersection(sys,leaf,leaf),mat,-1,2)<0)goto fail;
    return sys;
 fail:alea_destroy(sys);return NULL;
}
/* Direct world-space spheres provide a reference independent of lattice walking. */
static int ray_lattice_reference(const render_config_t* cfg,const render_camera_t* cam,const render_framebuffer_t* fb) {
    for(int y=0;y<cfg->height;++y)for(int x=0;x<cfg->width;++x) {
        double o[3],d[3];render_camera_ray(cam,cfg->width,cfg->height,x+.5,y+.5,o,o+1,o+2,d,d+1,d+2);
        double begin=0,end=INFINITY,best=INFINITY,normal[3]={0};
        if(cfg->num_clips){if(d[1]>0)begin=fmax(begin,-o[1]/d[1]);else if(d[1]<0)end=-o[1]/d[1];else if(o[1]<0)end=-1;}
        for(int z=-1;z<=1;++z)for(int j=-1;j<=1;++j)for(int i=-1;i<=1;++i) {
            double q[]={o[0]-(1+2*i),o[1]-2*j,o[2]-2*z};
            double a=d[0]*d[0]+d[1]*d[1]+d[2]*d[2],b=q[0]*d[0]+q[1]*d[1]+q[2]*d[2];
            double c=q[0]*q[0]+q[1]*q[1]+q[2]*q[2]-.36,disc=b*b-a*c;if(disc<0)continue;
            double lo=(-b-sqrt(disc))/a,hi=(-b+sqrt(disc))/a,hit=fmax(begin,lo);
            if(hit>=fmin(end,hi) || hit>=best)continue;best=hit;
            for(int k=0;k<3;++k)normal[k]=(q[k]+hit*d[k])/.6;
            if(hit==begin && begin>0){normal[0]=normal[2]=0;normal[1]=-1;}
        }
        size_t at=(size_t)y*cfg->width+x;int material=isfinite(best)?7:0;
        if(fb->material_id[at]!=material || (material && (fb->cell_id[at]!=42 || !isfinite(fb->depth[at]) || fabs(fb->depth[at]-best)>.003)))return -1;
        if(material)for(int k=0;k<3;++k)if(!isfinite(fb->normal[3*at+k]) || fabs(fb->normal[3*at+k]-normal[k])>.003)return -1;
    }
    return 0;
}
