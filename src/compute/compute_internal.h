// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
// SPDX-License-Identifier: MPL-2.0
#ifndef ALEA_COMPUTE_INTERNAL_H
#define ALEA_COMPUTE_INTERNAL_H
#include "compute_prepare.h"
#ifdef ALEA_USE_METAL
int alea_compute_metal_available(void);
void* alea_compute_metal_create(const alea_system_t* sys);
void alea_compute_metal_destroy(void* state);
int alea_compute_metal_classify(void* state, const double* xyz, size_t count,
                                alea_compute_result_t* results);
int alea_compute_metal_select(void* state,const double* xyz,size_t count,alea_compute_selection_t* results);
int alea_compute_metal_render(void* state, const render_config_t* cfg,
                              const render_camera_t* cam, render_framebuffer_t* fb,
                              alea_compute_render_progress_t progress, void* user);
#endif
#endif
