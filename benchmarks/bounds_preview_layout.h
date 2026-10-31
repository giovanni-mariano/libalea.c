// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
// SPDX-License-Identifier: MPL-2.0
#ifndef ALEA_BOUNDS_PREVIEW_LAYOUT_H
#define ALEA_BOUNDS_PREVIEW_LAYOUT_H
#ifdef __METAL_VERSION__
typedef uint bounds_uint;
#else
#include <stdint.h>
typedef uint32_t bounds_uint;
#endif
typedef struct {float b[6];} preview_box;
typedef struct {
    float eye[3],forward[3],right[3],up[3],half_width,half_height;
    bounds_uint perspective,filled;
} bounds_camera;
#define BOUNDS_PREVIEW_LIMIT 65536u
#define BOUNDS_PREVIEW_NEAR 0.001f
#define BOUNDS_PREVIEW_FAR 1.0e8f
#endif
