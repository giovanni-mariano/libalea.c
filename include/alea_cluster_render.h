// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

#ifndef ALEA_CLUSTER_RENDER_H
#define ALEA_CLUSTER_RENDER_H

#include "alea_cluster_base.h"
#include "alea_render.h"

#ifdef __cplusplus
extern "C" {
#endif

/** Collectively render statically assigned image tiles. Every rank supplies
 * an equivalent system, render configuration, and prepared camera. Only rank
 * zero supplies a framebuffer with cfg width/height; workers may pass NULL.
 * The framebuffer receives color and cell IDs, plus material IDs, depth and
 * normals when all three auxiliary arrays are present. On success, edges are
 * darkened on root when cfg->edges is set. On failure the framebuffer may be
 * partially written. MPI calls stay on the initializing thread. */
alea_cluster_status_t alea_cluster_render_scene(
    alea_cluster_t* cluster, alea_system_t* sys,
    const render_config_t* cfg, const render_camera_t* cam,
    render_framebuffer_t* root_framebuffer);

#ifdef __cplusplus
}
#endif

#endif /* ALEA_CLUSTER_RENDER_H */
