// SPDX-License-Identifier: MPL-2.0
#ifndef ALEA_RENDER_BENCHMARK_VIEW_H
#define ALEA_RENDER_BENCHMARK_VIEW_H

#include "alea.h"
#include "alea_render.h"

/* Cutaway, reverse cutaway, exterior, and close-up of the bounded geometry. */
static int benchmark_set_view(render_config_t* cfg, alea_system_t* sys, int view) {
    double center[3], radius;
    if (alea_compute_bounding_sphere(sys, 1.0, &center[0], &center[1], &center[2],
                                     &radius) != 0 || !(radius > 0)) return -1;
    double side = view == 2 ? -1.0 : 1.0;
    double distance = view == 4 ? 0.6 * radius : 1.6 * radius;
    cfg->eye_set = cfg->target_set = 1;
    for (int i = 0; i < 3; ++i) cfg->target[i] = center[i];
    cfg->eye[0] = center[0] + side * distance;
    cfg->eye[1] = center[1] - side * distance;
    cfg->eye[2] = center[2] + distance * 0.7;
    cfg->num_clips = view == 1 || view == 2;
    if (cfg->num_clips)
        cfg->clips[0] = (render_clip_plane_t){{0, side, 0}, -side * center[1]};
    return 0;
}
#endif
