// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
// SPDX-License-Identifier: MPL-2.0
#ifndef ALEA_GPU_GEOMETRY_FIXTURE_H
#define ALEA_GPU_GEOMETRY_FIXTURE_H
#include "alea.h"
#include "alea_compute.h"
#include <math.h>

static alea_system_t* scene(void) {
    alea_system_t* s = alea_create();
    if (!s) return NULL;
    int mat = alea_add_material(s, 7);
    int sphere = alea_sphere_surface(s, 1, 0, 0, 0, .75);
    alea_node_id_t solid=alea_halfspace(s,sphere,-1);
    /* Optional equivalent geometry with a long left-associated branch chain. */
    if(getenv("ALEA_BENCH_CSG_CHAIN"))for(int i=0;i<64;++i) {
        int outer=alea_sphere_surface(s,10000+i,0,0,0,1+i*.01);
        solid=alea_intersection(s,solid,alea_halfspace(s,outer,-1));
    }
    /* Equivalent union of intersection groups for structural diagnostics. */
    if(getenv("ALEA_BENCH_CSG_GROUPS")) {
        int a=alea_sphere_surface(s,20000,0,0,0,2);
        int b=alea_sphere_surface(s,20001,0,0,0,3);
        solid=alea_union(s,alea_intersection(s,solid,alea_halfspace(s,a,-1)),
                           alea_intersection(s,solid,alea_halfspace(s,b,-1)));
    }
    if (mat < 0 || sphere < 0 ||
        alea_add_cell(s, 42, solid, mat, -1, 1) < 0 ||
        alea_add_cell(s, 43, alea_halfspace(s, sphere, 1), ALEA_MATERIAL_VOID, 0, 1) < 0)
        goto fail;
    for (int y = 0; y < 8; ++y) for (int x = 0; x < 8; ++x) {
        int id = y * 8 + x;
        double tr[] = {2*x - 7, 2*y - 7, 0};
        int box = alea_box_surface(s, 100 + id, tr[0]-1, tr[0]+1,
                                   tr[1]-1, tr[1]+1, -1, 1);
        if (box < 0) goto fail;
        int parent = alea_add_cell(s, 1000 + id, alea_halfspace(s, box, -1),
                                   ALEA_MATERIAL_VOID, 0, 0);
        if (parent < 0 || alea_add_transform(s, id+1, tr, 3, 0) ||
            alea_set_fill(s, parent, 1, id+1)) goto fail;
    }
    return s;
fail:
    alea_destroy(s); return NULL;
}

static int check_points(const double* xyz, const alea_compute_result_t* r, size_t n) {
    for (size_t i = 0; i < n; ++i) {
        double x = xyz[3*i], y = xyz[3*i+1];
        double cx = 2*floor((x+8)/2)-7, cy = 2*floor((y+8)/2)-7;
        int solid = (x-cx)*(x-cx)+(y-cy)*(y-cy) <= .75*.75;
        if (r[i].owner_count != 1 || r[i].cell_id != (solid ? 42 : 43) ||
            r[i].material_id != (solid ? 7 : 0)) return -1;
    }
    return 0;
}

static int check_render(const render_config_t* cfg, const render_camera_t* cam,
                         const render_framebuffer_t* fb) {
    int mismatches=0;
    for (int y = 0; y < cfg->height; ++y) for (int x = 0; x < cfg->width; ++x) {
        double o[3], d[3];
        render_camera_ray(cam, cfg->width, cfg->height, x+.5, y+.5,
                           &o[0], &o[1], &o[2], &d[0], &d[1], &d[2]);
        double begin = 0, end = INFINITY, best = INFINITY, normal[3] = {0};
        if (cfg->num_clips) {
            if (d[1] > 0) begin = fmax(begin, -o[1]/d[1]);
            else if (d[1] < 0) end = -o[1]/d[1];
            else if (o[1] < 0) end = -1;
        }
        for (int j = 0; j < 8; ++j) for (int i = 0; i < 8; ++i) {
            double p[] = {o[0]-(2*i-7), o[1]-(2*j-7), o[2]};
            double a = d[0]*d[0]+d[1]*d[1]+d[2]*d[2];
            double b = p[0]*d[0]+p[1]*d[1]+p[2]*d[2];
            double c = p[0]*p[0]+p[1]*p[1]+p[2]*p[2]-.75*.75;
            double disc = b*b-a*c;
            if (disc < 0) continue;
            double lo = (-b-sqrt(disc))/a, hi = (-b+sqrt(disc))/a;
            double hit = fmax(begin, lo);
            if (hit >= fmin(end, hi) || hit >= best) continue;
            best = hit;
            for (int k = 0; k < 3; ++k) normal[k] = (p[k]+hit*d[k])/.75;
            if (hit == begin && begin > 0) { normal[0]=normal[2]=0; normal[1]=-1; }
        }
        size_t pixel = (size_t)y*cfg->width+x;
        int material = isfinite(best) ? 7 : 0;
        if (fb->material_id[pixel] != material ||
            (material && (fb->cell_id[pixel] != 42 || !isfinite(fb->depth[pixel]) || fabs(fb->depth[pixel]-best) > .003)))
            {++mismatches;continue;}
        if (material) for (int k = 0; k < 3; ++k)
            if (!isfinite(fb->normal[3*pixel+k]) || fabs(fb->normal[3*pixel+k]-normal[k]) > .003) {++mismatches;break;}
    }
    return mismatches;
}

#endif
