// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
// SPDX-License-Identifier: MPL-2.0
#ifndef ALEA_NATIVE_PROBE_KERNEL_H
#define ALEA_NATIVE_PROBE_KERNEL_H

#include <stdint.h>
#include <math.h>
#if defined(PROBE_CUDA) || defined(PROBE_HIP)
#define PROBE_DEVICE __host__ __device__
#else
#define PROBE_DEVICE
#endif

/* Prototype wire types, deliberately independent of Metal's private layout.
 * Only intersections of sphere interiors and negative plane halfspaces.
 * This is NOT libalea's general CSG evaluator or boundary contract. */
struct ProbePrimitive { uint32_t type; float x, y, z, w; };
struct ProbePoint { float x, y, z; };
struct ProbeRay { ProbePoint origin, direction; };
struct ProbeState { float enter, exit; uint32_t cursor; int32_t status; };
enum { PROBE_ACTIVE = 0, PROBE_HIT = 1, PROBE_MISS = 2, PROBE_INVALID = 3 };
static_assert(sizeof(ProbePrimitive) == 20, "primitive wire layout");
static_assert(sizeof(ProbeState) == 16, "continuation wire layout");

PROBE_DEVICE inline int probe_contains(const ProbePrimitive *p, uint32_t n,
                                       ProbePoint q) {
    for (uint32_t i = 0; i < n; ++i) {
        float v;
        if (p[i].type == 0) {
            float x = q.x-p[i].x, y = q.y-p[i].y, z = q.z-p[i].z;
            v = x*x+y*y+z*z-p[i].w*p[i].w;
        } else v = p[i].x*q.x+p[i].y*q.y+p[i].z*q.z+p[i].w;
        if (v > 0) return 0;
    }
    return 1;
}

PROBE_DEVICE inline void probe_slice(const ProbePrimitive *p, uint32_t n,
                                    ProbeRay ray, ProbeState *s, uint32_t budget) {
    if (s->status != PROBE_ACTIVE) return;
    float a = ray.direction.x*ray.direction.x + ray.direction.y*ray.direction.y
            + ray.direction.z*ray.direction.z;
    if (a == 0) { s->status = PROBE_INVALID; return; }
    while (budget-- && s->cursor < n) {
        const ProbePrimitive v = p[s->cursor++];
        float lo = -INFINITY, hi = INFINITY;
        if (v.type == 0) {
            float x=ray.origin.x-v.x, y=ray.origin.y-v.y, z=ray.origin.z-v.z;
            float b=x*ray.direction.x+y*ray.direction.y+z*ray.direction.z;
            float c=x*x+y*y+z*z-v.w*v.w, d=b*b-a*c;
            if (d < 0) { s->status=PROBE_MISS; return; }
            float q=-b-copysignf(sqrtf(d), b);
            if (q == 0) lo=hi=-b/a;
            else { lo=fminf(q/a,c/q); hi=fmaxf(q/a,c/q); }
        } else {
            float b=v.x*ray.direction.x+v.y*ray.direction.y+v.z*ray.direction.z;
            float c=v.x*ray.origin.x+v.y*ray.origin.y+v.z*ray.origin.z+v.w;
            if (b == 0) {
                if (c > 0) { s->status=PROBE_MISS; return; }
            } else if (b > 0) hi=-c/b;
            else lo=-c/b;
        }
        s->enter=fmaxf(s->enter,lo); s->exit=fminf(s->exit,hi);
        if (s->enter > s->exit) { s->status=PROBE_MISS; return; }
    }
    if (s->cursor == n) s->status=PROBE_HIT;
}
#endif
