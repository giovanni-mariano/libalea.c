// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

#ifndef ALEA_RENDER_TILES_INTERNAL_H
#define ALEA_RENDER_TILES_INTERNAL_H

#include "alea_render.h"

/* Tile storage is compact: each tile occupies tile_size * tile_size slots,
 * including unused slots on cropped right/bottom boundary tiles. */
int alea_render_scene_tile_span(alea_system_t* sys,
    const render_config_t* cfg, const render_camera_t* cam,
    size_t first_tile, size_t tile_count, render_framebuffer_t* compact);

#endif
