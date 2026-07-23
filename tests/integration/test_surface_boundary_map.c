// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

#include "alea_test.h"
#include "alea.h"
#include "alea_slice.h"

static alea_system_t* build_sphere_graveyard(void) {
    alea_system_t* sys = alea_create();
    if (!sys) return NULL;

    int sphere = alea_sphere_surface(sys, 1, 0.0, 0.0, 0.0, 1.0);
    if (sphere < 0) goto fail;
    alea_node_id_t inside = alea_halfspace(sys, sphere, -1);
    alea_node_id_t outside = alea_halfspace(sys, sphere, 1);
    if (inside == ALEA_NODE_ID_INVALID || outside == ALEA_NODE_ID_INVALID)
        goto fail;
    int material = alea_add_material(sys, 1);
    if (material < 0 ||
        alea_add_cell(sys, 1, inside, material, 1.0, 0) < 0 ||
        alea_add_cell(sys, 99, outside, ALEA_MATERIAL_VOID, 0.0, 0) < 0)
        goto fail;
    return sys;
fail:
    alea_destroy(sys);
    return NULL;
}

TEST(surface_boundary_map_attributes_sphere_graveyard_edge) {
    alea_system_t* sys = build_sphere_graveyard();
    ASSERT_NOT_NULL(sys);
    ASSERT_EQ(alea_prepare_query_acceleration(sys), 0);

    const int width = 64, height = 64;
    int ids[width * height];
    alea_slice_view_t view;
    alea_slice_view_axis(&view, 2, 0.0, -2.0, 2.0, -2.0, 2.0);
    ASSERT_EQ(alea_find_cells_grid(sys, &view, width, height, -1,
                                   ids, NULL, NULL), 0);

    alea_slice_surface_boundary_map_t* map = NULL;
    ASSERT_EQ(alea_slice_surface_boundary_map_create(
                  sys, &view, width, height, ids,
                  alea_slice_classify_cell, NULL, &map), 0);
    ASSERT_NOT_NULL(map);

    int changed_edges = 0;
    int valid_edges = 0;
    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            if (x + 1 < width && ids[y * width + x] != ids[y * width + x + 1]) {
                changed_edges++;
                alea_slice_boundary_status_t status =
                    alea_slice_surface_boundary_status(
                        map, x, y, ALEA_SLICE_EDGE_RIGHT);
                if (status ==
                    ALEA_SLICE_BOUNDARY_VALID &&
                    alea_slice_surface_boundary_surface_count(
                        map, x, y, ALEA_SLICE_EDGE_RIGHT) == 1 &&
                    alea_slice_surface_boundary_surface_id(
                        map, x, y, ALEA_SLICE_EDGE_RIGHT, 0) == 1)
                    valid_edges++;
            }
            if (y > 0 && ids[y * width + x] != ids[(y - 1) * width + x]) {
                changed_edges++;
                alea_slice_boundary_status_t status =
                    alea_slice_surface_boundary_status(
                        map, x, y, ALEA_SLICE_EDGE_DOWN);
                if (status ==
                    ALEA_SLICE_BOUNDARY_VALID &&
                    alea_slice_surface_boundary_surface_count(
                        map, x, y, ALEA_SLICE_EDGE_DOWN) == 1 &&
                    alea_slice_surface_boundary_surface_id(
                        map, x, y, ALEA_SLICE_EDGE_DOWN, 0) == 1)
                    valid_edges++;
            }
        }
    }
    ASSERT(changed_edges > 0);
    ASSERT_EQ(valid_edges, changed_edges);

    alea_slice_surface_boundary_map_free(map);
    alea_destroy(sys);
}

TEST_MAIN()
