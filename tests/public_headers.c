// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/*
 * Compile-only check: every installed header must resolve with -Iinclude and
 * without src/ on the include path. Keep this list in sync with include/.
 */
#include "alea.h"
#include "alea_geo_validator.h"
#include "alea_log.h"
#include "alea_mcnp.h"
#include "alea_mesh.h"
#include "alea_nucdata.h"
#include "alea_nucdata_types.h"
#include "alea_openmc.h"
#include "alea_raycast.h"
#include "alea_render.h"
#include "alea_serpent.h"
#include "alea_slice.h"
#include "alea_types.h"

#ifdef alea_surface_at
#error "alea_surface_at must remain an internal-only helper"
#endif

void public_api_smoke(void) {
    alea_system_t* sys = alea_create();
    int surface = alea_sphere_surface(sys, 0, 0.0, 0.0, 0.0, 1.0);
    alea_node_id_t inside = alea_halfspace(sys, surface, -1);
    alea_node_id_t outside = alea_halfspace(sys, surface, +1);
    alea_raycast_result_t* result = alea_raycast_result_create();
    alea_ray_coverage_slice_result_t* coverage =
        alea_ray_coverage_slice_result_create();
    alea_ray_coverage_slice_options_t coverage_options;
    uint8_t flags = 0;
    double hit_t = 0.0;
    int hit_surface_id = -1;

    (void)inside;
    (void)outside;
    (void)alea_raycast_segment_resolution_flags(result, 0, &flags);
    (void)alea_raycast_hit_count(result);
    (void)alea_raycast_hit_get(result, 0, &hit_t, &hit_surface_id);
    alea_ray_coverage_slice_options_init(&coverage_options);
    (void)alea_ray_coverage_slice_row_count(coverage);
    (void)alea_ray_coverage_slice_row_offsets(coverage);
    alea_log_set_callback(NULL, NULL);
    alea_raycast_result_destroy(result);
    alea_ray_coverage_slice_result_destroy(coverage);
    alea_destroy(sys);
}
