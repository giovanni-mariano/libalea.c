// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
// SPDX-License-Identifier: MPL-2.0

#include "alea_test.h"
#include "alea.h"
#include "core/alea_system.h"
#include "raycast/raycast.h"
#include "raycast/ray_epsilon.h"
#include <float.h>

/* Child-universe cards exercise indexed scans; the parent has no crossings
 * in the tested windows. The rotation makes the global ray parallel to the
 * unplaced plane, so only its occurrence can produce the expected hit.
 * Payload edits below intentionally bypass ingestion
 * normalization to test the scan's existing floating-point contract. */
static alea_system_t* plane_model(void) {
    alea_system_t* sys = alea_create();
    if (!sys) abort();
    int plane = alea_plane_surface(sys, 10, 1, 0, 0, -2);
    int box = alea_box_surface(sys, 20, -100, 100, -100, 100, -100, 100);
    if (plane < 0 || box < 0 ||
        alea_add_cell(sys, 1, alea_halfspace(sys, plane, -1),
                      ALEA_MATERIAL_VOID, 0, 1) < 0) abort();
    int parent = alea_add_cell(sys, 2, alea_halfspace(sys, box, -1),
                               ALEA_MATERIAL_VOID, 0, 0);
    const double rotation[] = {0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, -1};
    if (parent < 0 || alea_add_transform(sys, 1, rotation, 12, 0) ||
        alea_set_fill(sys, parent, 1, 1) ||
        alea_prepare_query_acceleration(sys)) abort();
    return sys;
}

static int query(alea_system_t* sys, const alea_ray_t* ray,
                 double lo, double hi, int nearest, alea_raycast_result_t* out) {
    alea_ray_t initialized;
    alea_ray_init_normalized(&initialized, ray->ox, ray->oy, ray->oz,
                             ray->dx, ray->dy, ray->dz);
    alea_raycast_result_init(out);
    out->breakpoint_nearest_only = nearest;
    return alea_raycast_validation_breakpoints_reuse_nocache(
        sys, &initialized, lo, hi, 64, out);
}

TEST(nearest_plane_scan_preserves_endpoints_and_parallel_threshold) {
    alea_system_t* sys = plane_model();
    alea_ray_t ray = {.dy = 1};
    const double windows[][2] = {
        {0, 2}, {2, 3}, {0, nextafter(2, 0)},
        {nextafter(2, INFINITY), 3}, {nextafter(2, 0), nextafter(2, INFINITY)}
    };
    const size_t expected[] = {1, 1, 0, 0, 1};
    for (int reverse = 0; reverse < 2; reverse++) {
        ray.oy = reverse ? 4 : 0;
        ray.dy = reverse ? -1 : 1;
        for (size_t i = 0; i < sizeof(expected)/sizeof(expected[0]); i++) {
            for (int nearest = 0; nearest < 2; nearest++) {
                alea_raycast_result_t result;
                ASSERT_EQ(query(sys, &ray, windows[i][0], windows[i][1], nearest, &result), 0);
                ASSERT_EQ(result.hits.count, expected[i]);
                if (expected[i]) {
                    ASSERT(result.hits.data[0].t == 2);
                    ASSERT_EQ(result.hits.data[0].surface_id, 10);
                    ASSERT(result.hits.data[0].ny == 1);
                    ASSERT(result.hits.data[0].nx == 0 && result.hits.data[0].nz == 0);
                }
                alea_raycast_result_free(&result);
            }
        }
    }
    uint32_t primitive = sys->surfaces.data[0].primitive_id;
    alea_plane_data_t* plane = &sys->primitive_planes.data[sys->primitives.data[primitive].payload_index];
    const double directions[] = {0, -0.0, nextafter(RAY_EPSILON, 0),
                                 RAY_EPSILON, nextafter(RAY_EPSILON, INFINITY), -RAY_EPSILON};
    for (size_t i = 0; i < sizeof(directions)/sizeof(directions[0]); i++) {
        ray.oy = 0; ray.dy = directions[i];
        ray.dz = sqrt(1 - ray.dy * ray.dy);
        plane->d = -2 * ray.dy;
        for (int nearest = 0; nearest < 2; nearest++) {
            alea_raycast_result_t result;
            ASSERT_EQ(query(sys, &ray, 0, 4, nearest, &result), 0);
            ASSERT_EQ(result.hits.count, i < 3 ? 0 : 1);
            if (result.hits.count) ASSERT(result.hits.data[0].t == 2);
            alea_raycast_result_free(&result);
        }
    }
    alea_destroy(sys);
}

TEST(nearest_plane_scan_matches_generic_path_for_extreme_payloads) {
    alea_system_t* sys = plane_model();
    uint32_t primitive = sys->surfaces.data[0].primitive_id;
    alea_plane_data_t* plane = &sys->primitive_planes.data[sys->primitives.data[primitive].payload_index];
    const alea_plane_data_t cases[] = {
        {1, -0.0, 0, -2}, {1e100, 0, 0, -2e100}, {1e-100, 0, 0, -2e-100},
        {DBL_MAX/4, 0, 0, -DBL_MAX/2}, {1, 0, 0, INFINITY},
        {1, 0, 0, NAN}, {INFINITY, 0, 0, -2}, {NAN, 0, 0, -2}
    };
    const size_t expected[] = {1, 1, 0, 1, 0, 0, 0, 0};
    alea_ray_t ray = {.dy = 1};
    for (size_t i = 0; i < sizeof(cases)/sizeof(cases[0]); i++) {
        *plane = cases[i];
        alea_raycast_result_t direct, generic;
        ASSERT_EQ(query(sys, &ray, 0, 4, 1, &direct), 0);
        ASSERT_EQ(query(sys, &ray, 0, 4, 0, &generic), 0);
        ASSERT_EQ(direct.hits.count, expected[i]);
        ASSERT_EQ(direct.hits.count, generic.hits.count);
        if (direct.hits.count) {
            const alea_ray_hit_t* a = &direct.hits.data[0];
            const alea_ray_hit_t* b = &generic.hits.data[0];
            ASSERT(a->t == 2);
            ASSERT_EQ(a->surface_id, b->surface_id);
            ASSERT_EQ(a->primitive_id, b->primitive_id);
            ASSERT(memcmp(&a->t, &b->t, sizeof(double)) == 0);
            ASSERT(memcmp(&a->nx, &b->nx, sizeof(double)) == 0);
            ASSERT(memcmp(&a->ny, &b->ny, sizeof(double)) == 0);
            ASSERT(memcmp(&a->nz, &b->nz, sizeof(double)) == 0);
        }
        alea_raycast_result_free(&direct);
        alea_raycast_result_free(&generic);
    }
    alea_destroy(sys);
}

TEST(nearest_plane_scan_preserves_invalid_payload_failure_and_reuse) {
    alea_system_t* sys = plane_model();
    alea_cell_entry_t* child = &sys->cells.data[alea_cell_find(sys, 1)];
    ASSERT_NOT_NULL(child->surface_indices);
    uint32_t card = child->surface_indices[0];
    uint32_t primitive = sys->surfaces.data[card].primitive_id;
    uint32_t payload = sys->primitives.data[primitive].payload_index;
    alea_ray_t ray = {.dy = 1};
    for (int fault = 0; fault < 2; fault++) {
        if (fault == 0) child->surface_indices[0] = UINT32_MAX;
        if (fault == 1) sys->primitives.data[primitive].payload_index = UINT32_MAX;
        alea_raycast_result_t result;
        int rc = query(sys, &ray, 0, 4, 1, &result);
        child->surface_indices[0] = card;
        sys->surfaces.data[card].primitive_id = primitive;
        sys->primitives.data[primitive].payload_index = payload;
        ASSERT(rc != 0);
        ASSERT(result.breakpoint_failed);
        alea_raycast_result_free(&result);
        ASSERT_EQ(query(sys, &ray, 0, 4, 1, &result), 0);
        ASSERT_EQ(result.hits.count, 1);
        ASSERT(result.hits.data[0].t == 2);
        alea_raycast_result_free(&result);
    }
    alea_destroy(sys);
}

TEST(nearest_plane_scan_keeps_crossings_outside_terminal_cell_bounds) {
    alea_system_t* sys = plane_model();
    int child = alea_cell_find(sys, 1);
    int box = alea_box_surface(sys, 30, 10, 20, -1, 1, -1, 1);
    ASSERT(box >= 0);
    alea_node_id_t region = alea_intersection(sys,
        sys->cells.data[child].root_node_id, alea_halfspace(sys, box, -1));
    ASSERT_EQ(alea_cell_set_region(sys, child, region), 0);
    ASSERT_EQ(alea_prepare_query_acceleration(sys), 0);
    alea_ray_t ray = {.dy = 1};
    alea_raycast_result_t result;
    ASSERT_EQ(query(sys, &ray, 0, 4, 1, &result), 0);
    /* A referenced plane remains a diagnostic crossing even though the
     * terminal region cannot own any point in this window. */
    ASSERT_EQ(result.hits.count, 1);
    ASSERT(result.hits.data[0].t == 2);
    ASSERT_EQ(result.hits.data[0].surface_id, 10);
    alea_raycast_result_free(&result);
    alea_interrupt();
    int rc = query(sys, &ray, 0, 4, 1, &result);
    alea_clear_interrupt();
    ASSERT(rc != 0);
    alea_raycast_result_free(&result);
    ASSERT_EQ(query(sys, &ray, 0, 4, 1, &result), 0);
    ASSERT_EQ(result.hits.count, 1);
    ASSERT(result.hits.data[0].t == 2);
    alea_raycast_result_free(&result);
    alea_destroy(sys);
}

TEST_MAIN()
