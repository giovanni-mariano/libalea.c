// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
// SPDX-License-Identifier: MPL-2.0
#define _POSIX_C_SOURCE 200809L
#include "alea.h"
#include "alea_compute.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

static double seconds(void) {
    struct timespec t; clock_gettime(CLOCK_MONOTONIC,&t);
    return t.tv_sec+t.tv_nsec*1e-9;
}

static double measure(alea_system_t* sys, alea_compute_t* gpu, render_framebuffer_t* fb) {
    render_config_t cfg; render_config_init(&cfg);
    cfg.width=fb->width; cfg.height=fb->height; cfg.log_level=0;
    cfg.eye_set=cfg.target_set=1;
    cfg.target[0]=cfg.target[1]=cfg.target[2]=0;
    double total=0;
    for (int frame=0;frame<7;++frame) {
        double angle=.6+frame*.04;
        cfg.eye[0]=18*cos(angle); cfg.eye[1]=-18*sin(angle); cfg.eye[2]=9;
        render_camera_t cam;
        if (render_camera_setup(&cam,&cfg,sys)) return -1;
        double t=seconds();
        int rc=gpu?alea_compute_render(gpu,&cfg,&cam,fb):render_scene(sys,&cfg,&cam,fb);
        if (rc) { fprintf(stderr,"%s\n",alea_error()); return -1; }
        t=seconds()-t;
        if (!frame) printf("%s first frame %.3f ms\n",gpu?"Metal":"CPU",t*1000);
        else total+=t;
    }
    render_config_free(&cfg);
    return total/6;
}

int main(void) {
    int status=1;
    alea_system_t* sys=alea_create();
    alea_compute_t* gpu=NULL;
    render_framebuffer_t* cpu=render_framebuffer_create(800,600,1);
    render_framebuffer_t* metal=render_framebuffer_create(800,600,1);
    if (!sys || !cpu || !metal) goto done;
    int mat=alea_add_material(sys,7);
    for (int x=0;x<4;++x) for (int y=0;y<4;++y) for (int z=0;z<4;++z) {
        int surf=alea_sphere_surface(sys,0,2*x-3,2*y-3,2*z-3,.85);
        if (surf<0 || alea_add_cell(sys,0,alea_halfspace(sys,surf,-1),mat,-1,0)<0) goto done;
    }
    double t=seconds();
    if (alea_prepare_query_acceleration(sys)) goto done;
    printf("CPU hierarchy preparation %.3f ms\n",(seconds()-t)*1000);
    t=seconds();
    gpu=alea_compute_create(sys,ALEA_COMPUTE_METAL,ALEA_COMPUTE_FP32);
    if (!gpu) { fprintf(stderr,"%s\n",alea_error()); goto done; }
    printf("Metal scene preparation %.3f ms\n",(seconds()-t)*1000);
    puts("64 spheres, 800x600, moving perspective camera, AA1, six warm frames.");
    puts("CPU: existing hierarchical raycaster; Metal: cell BVH + CSG intervals.");
    double a=measure(sys,NULL,cpu), b=measure(sys,gpu,metal);
    if (a<=0 || b<=0) goto done;
    size_t different=0, missing_normal=0; double max_depth=0, max_color=0, valid_color=0;
    for (int i=0;i<800*600;++i) {
        if (cpu->cell_id[i]!=metal->cell_id[i]) ++different;
        else if (cpu->cell_id[i]>=0) max_depth=fmax(max_depth,fabs(cpu->depth[i]-metal->depth[i]));
        int normal_missing=cpu->cell_id[i]>=0 && cpu->normal[3*i]==0 &&
            cpu->normal[3*i+1]==0 && cpu->normal[3*i+2]==0;
        missing_normal+=normal_missing;
        for (int k=0;k<3;++k) {
            double delta=fabs(cpu->color[3*i+k]-metal->color[3*i+k]);
            max_color=fmax(max_color,delta);
            if (!normal_missing) valid_color=fmax(valid_color,delta);
        }
    }
    printf("CPU %.3f ms/frame, Metal %.3f ms/frame, ratio %.2fx\n",a*1000,b*1000,a/b);
    printf("Last frame: ID differences %zu/480000, max depth delta %.6g, max RGB delta %.6g\n",different,max_depth,max_color);
    printf("CPU hit pixels without normals: %zu; max RGB delta on remaining pixels %.6g\n",missing_normal,valid_color);
    puts("Timings include GPU synchronization/readback, exclude SDL display and PNG output.");
    status=0;
done:
    alea_compute_destroy(gpu); alea_destroy(sys);
    render_framebuffer_free(cpu); render_framebuffer_free(metal);
    return status;
}
