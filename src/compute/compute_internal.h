// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
// SPDX-License-Identifier: MPL-2.0
#ifndef ALEA_COMPUTE_INTERNAL_H
#define ALEA_COMPUTE_INTERNAL_H
#include "alea_compute.h"
#ifdef ALEA_USE_METAL
#include "compute_layout.h"
typedef struct {
    compute_render_bvh* bvh;
    compute_transform* transforms;
    compute_lattice* lattices;
    compute_uint* fills;
    size_t bvh_count, transform_count, lattice_count, fill_count;
    compute_uint root;
} compute_render_scene;
typedef struct {
    compute_node* nodes;
    compute_primitive* primitives;
    compute_cell* cells;
    size_t node_count, primitive_count, cell_count;
    compute_spatial_program* spatial;
    size_t spatial_count;
    compute_render_scene render;
} compute_scene;
int alea_compute_primitive_sign(const compute_primitive* primitive, const float bounds[6]);
int alea_compute_spatial_pack(compute_scene* scene);
int alea_compute_bvh_rebuild(compute_scene* scene);
/* Internal entry point also used to check conservative specialization. */
int alea_compute_local_simplify(const compute_node* code, size_t count,
    const compute_primitive* primitives, const float bounds[6],
    compute_node* output, size_t* output_count);
int alea_compute_pack(const alea_system_t* sys, compute_scene* scene);
int alea_compute_pack_geometry(const alea_system_t* sys, compute_scene* scene, int accelerate);
void alea_compute_scene_free(compute_scene* scene);
int alea_compute_render_pack(const alea_system_t* sys, compute_scene* scene);
void alea_compute_render_free(compute_render_scene* render);
int alea_compute_profile_write(const char* prefix, const compute_profile* rays,
    const compute_profile* costs, const compute_cell* cells, size_t cell_count,
    unsigned width, unsigned height);
int alea_compute_metal_available(void);
void* alea_compute_metal_create(const alea_system_t* sys);
void alea_compute_metal_destroy(void* state);
int alea_compute_metal_classify(void* state, const double* xyz, size_t count,
                                alea_compute_result_t* results);
int alea_compute_metal_render(void* state, const render_config_t* cfg,
                              const render_camera_t* cam, render_framebuffer_t* fb,
                              alea_compute_render_progress_t progress, void* user);
#endif
#endif
