// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
// SPDX-License-Identifier: MPL-2.0
#ifndef ALEA_COMPUTE_LAYOUT_H
#define ALEA_COMPUTE_LAYOUT_H
#ifdef __METAL_VERSION__
typedef uint compute_uint;
typedef int compute_int;
#else
#include <stdint.h>
typedef uint32_t compute_uint;
typedef int32_t compute_int;
#endif
#define COMPUTE_NONE 0xffffffffu
#define COMPUTE_STACK 64
#define COMPUTE_HIER_STACK 32
#define COMPUTE_BVH_STACK 256
#define COMPUTE_CELL_LATTICE_SUBTREE 1u

enum { COMPUTE_LEAF, COMPUTE_UNION, COMPUTE_INTERSECTION,
       COMPUTE_DIFFERENCE, COMPUTE_COMPLEMENT, COMPUTE_REVERSE_DIFFERENCE,
       COMPUTE_SKIP_UNION, COMPUTE_SKIP_INTERSECTION, COMPUTE_SKIP_REVERSE_DIFFERENCE,
       COMPUTE_CONSTANT };
enum { COMPUTE_PLANE, COMPUTE_SPHERE, COMPUTE_CYL_X, COMPUTE_CYL_Y,
       COMPUTE_CYL_Z, COMPUTE_BOX, COMPUTE_QUADRIC,
       COMPUTE_CONE_X, COMPUTE_CONE_Y, COMPUTE_CONE_Z,
       COMPUTE_TORUS_X, COMPUTE_TORUS_Y, COMPUTE_TORUS_Z, COMPUTE_RCC };

/* Postfix CSG instructions decouple device stack size from source tree depth.
 * LEAF arg packs the primitive index and halfspace flip in its low bit. */
/* Metal-private point traversal links, separate from visibility-only BVHs.
 * Every cell, including terminal void, remains eligible for ownership. */
typedef struct { compute_uint next, child; } metal_point_cell;
typedef struct { compute_uint op, arg; } compute_node;
typedef struct { compute_uint type; float p[13]; } compute_primitive;
typedef struct {
    compute_uint root, count;
    compute_int cell_id, material_id;
    compute_uint child, transform, lattice, spatial;
    compute_uint flags; /* Low bit: lattice descendants; upper bits: child query cost. */
} compute_cell;
/* Adaptive render-only programs in the containing universe's coordinates.
 * Children are adjacent; the right child's lower bound is the split plane. */
typedef struct {
    float bounds[6];
    compute_uint root, count, child, axis;
} compute_spatial_program;
typedef struct {
    compute_int cell_id, material_id;
    compute_uint owner_count;
} compute_result;
/* Concatenated per-universe BVHs; child refers to the root of another BVH.
 * Escape links identify subtree ends and the right child during near-first traversal. */
typedef struct {
    float bounds[6];
    compute_uint escape, cell;
    /* Original bounds and a coordinate budget for conservative FP32 fallback. */
    float fallback[6], guard;
} compute_render_bvh;
typedef struct { float inverse[12]; } compute_transform;
typedef struct {
    compute_int dims[6];
    float pitch[3], lower[3];
    compute_uint offset, repeating, zero_coords;
    compute_uint outer; /* BVH root, or NONE */
} compute_lattice;
typedef struct {
    float eye[3], forward[3], right[3], up[3], light[3], background[3], clip[4];
    float half_width, half_height, ambient, diffuse, specular, shininess, tint;
    compute_uint width, height, perspective, light_set, clip_enabled, bvh_root, aa;
} compute_render_params;
typedef struct { compute_uint x, y, width, offset; } compute_render_tile;
typedef struct { compute_uint work, reset; } compute_resume_params;
typedef struct {
    float color[3], depth, normal[3];
    compute_int cell_id, material_id, error;
} compute_render_pixel;
/* Saturating diagnostic counters; field order is also used by the Metal
 * atomic cell accumulator. Normal pipelines compile instrumentation out. */
typedef struct {
    compute_uint bvh, cells, instructions, primitives, intervals, hierarchy;
    compute_uint parents, cache_hits, rejected;
} compute_profile;
#endif
