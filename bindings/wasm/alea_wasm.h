// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
// SPDX-License-Identifier: MPL-2.0

#ifndef ALEA_WASM_BINDING_H
#define ALEA_WASM_BINDING_H

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

int alea_wasm_init(int width, int height);
int alea_wasm_resize(int width, int height);
int alea_wasm_load_mcnp(const char* input, int length);
int alea_wasm_render(double azimuth, double elevation, double distance_scale,
                     double clip_fraction, int edges);
const uint8_t* alea_wasm_pixels(void);
int alea_wasm_pixel_bytes(void);
int alea_wasm_width(void);
int alea_wasm_height(void);
int alea_wasm_cell_count(void);
int alea_wasm_surface_count(void);
double alea_wasm_center_x(void);
double alea_wasm_center_y(void);
double alea_wasm_center_z(void);
double alea_wasm_radius(void);
int alea_wasm_set_target(double x, double y, double z);
int alea_wasm_parallel_enabled(void);
int alea_wasm_parallel_max_threads(void);
int alea_wasm_openmp_enabled(void);
int alea_wasm_openmp_max_threads(void);
const char* alea_wasm_last_error(void);
void alea_wasm_destroy(void);

#ifdef __cplusplus
}
#endif

#endif
