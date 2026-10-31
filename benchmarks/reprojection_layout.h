// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
// SPDX-License-Identifier: MPL-2.0
#ifndef ALEA_BENCH_REPROJECTION_LAYOUT_H
#define ALEA_BENCH_REPROJECTION_LAYOUT_H
#ifdef __METAL_VERSION__
typedef uint preview_uint;
#else
#include <stdint.h>
typedef uint32_t preview_uint;
#endif
/* Display-only cache; no geometry ownership or authoritative query outputs. */
typedef struct {float world[3],color[3];preview_uint valid;} preview_point;
typedef struct {
    float eye[3],forward[3],right[3],up[3],half_width,half_height;
    preview_uint width,height,perspective,count;
} preview_params;
typedef struct {preview_uint pixel,depth;} preview_projection;
#define PREVIEW_NONE 0xffffffffu
#endif
