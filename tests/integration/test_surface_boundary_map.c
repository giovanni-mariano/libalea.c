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

TEST(surface_boundary_map_ignores_cell_only_event_for_material_grid) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);
    int first = alea_plane_surface(sys, 1, 1.0, 0.0, 0.0, 0.0);
    int second = alea_plane_surface(sys, 2, 1.0, 0.0, 0.0, -0.1);
    ASSERT(first >= 0 && second >= 0);
    int material_one = alea_add_material(sys, 1);
    int material_two = alea_add_material(sys, 2);
    ASSERT(material_one >= 0 && material_two >= 0);

    /* The sole coarse material edge crosses both surfaces. Surface 1 changes
     * cell ownership but leaves material 1 selected; surface 2 is the only
     * causal material transition and therefore the only label participant. */
    alea_node_id_t first_neg = alea_halfspace(sys, first, -1);
    alea_node_id_t first_pos = alea_halfspace(sys, first, 1);
    alea_node_id_t second_neg = alea_halfspace(sys, second, -1);
    alea_node_id_t second_pos = alea_halfspace(sys, second, 1);
    ASSERT(first_neg != ALEA_NODE_ID_INVALID &&
           first_pos != ALEA_NODE_ID_INVALID &&
           second_neg != ALEA_NODE_ID_INVALID &&
           second_pos != ALEA_NODE_ID_INVALID);
    ASSERT(alea_add_cell(sys, 1, first_neg,
                         material_one, 1.0, 0) >= 0);
    ASSERT(alea_add_cell(sys, 2,
                         alea_intersection(sys,
                             first_pos, second_neg),
                         material_one, 1.0, 0) >= 0);
    ASSERT(alea_add_cell(sys, 3, second_pos,
                         material_two, 1.0, 0) >= 0);
    ASSERT_EQ(alea_prepare_query_acceleration(sys), 0);

    const int width = 2, height = 1;
    int cell_ids[width * height], material_ids[width * height];
    alea_slice_view_t view;
    alea_slice_view_axis(&view, 2, 0.0, -1.0, 1.0, -0.5, 0.5);
    ASSERT_EQ(alea_find_cells_grid(sys, &view, width, height, -1,
                                   cell_ids, material_ids, NULL), 0);
    ASSERT_EQ(material_ids[0], 1);
    ASSERT_EQ(material_ids[1], 2);

    alea_slice_surface_boundary_map_t* map = NULL;
    ASSERT_EQ(alea_slice_surface_boundary_map_create(
                  sys, &view, width, height, material_ids,
                  alea_slice_classify_material, NULL, &map), 0);
    ASSERT_NOT_NULL(map);
    ASSERT_EQ(alea_slice_surface_boundary_status(
                  map, 0, 0, ALEA_SLICE_EDGE_RIGHT),
              ALEA_SLICE_BOUNDARY_VALID);
    ASSERT_EQ(alea_slice_surface_boundary_surface_count(
                  map, 0, 0, ALEA_SLICE_EDGE_RIGHT), 1);
    ASSERT_EQ(alea_slice_surface_boundary_surface_id(
                  map, 0, 0, ALEA_SLICE_EDGE_RIGHT, 0), 2);
    ASSERT_EQ(alea_slice_surface_boundary_group_count(
                  map, 0, 0, ALEA_SLICE_EDGE_RIGHT), 1);
    ASSERT(alea_slice_surface_boundary_group_fraction(
               map, 0, 0, ALEA_SLICE_EDGE_RIGHT, 0) > 0.5);
    ASSERT_EQ(alea_slice_surface_boundary_group_surface_count(
                  map, 0, 0, ALEA_SLICE_EDGE_RIGHT, 0), 1);
    ASSERT_EQ(alea_slice_surface_boundary_group_surface_id(
                  map, 0, 0, ALEA_SLICE_EDGE_RIGHT, 0, 0), 2);

    alea_slice_surface_boundary_map_free(map);

    /* The same causal filtering must survive connected-arc label extraction:
     * the long material contour has enough edges to earn a label, but its
     * same-material cell transition must not leak surface 1 into the result. */
    const int label_width = 64, label_height = 64;
    int* label_cell_ids = malloc(
        (size_t)label_width * label_height * sizeof(*label_cell_ids));
    int* label_material_ids = malloc(
        (size_t)label_width * label_height * sizeof(*label_material_ids));
    ASSERT_NOT_NULL(label_cell_ids);
    ASSERT_NOT_NULL(label_material_ids);
    alea_slice_view_axis(&view, 2, 0.0, -1.0, 1.0, -1.0, 1.0);
    ASSERT_EQ(alea_find_cells_grid(sys, &view, label_width, label_height, -1,
                                   label_cell_ids, label_material_ids, NULL), 0);
    map = NULL;
    ASSERT_EQ(alea_slice_surface_boundary_map_create(
                  sys, &view, label_width, label_height, label_material_ids,
                  alea_slice_classify_material, NULL, &map), 0);
    ASSERT_NOT_NULL(map);
    alea_label_position_t* labels = NULL;
    int label_count = 0;
    ASSERT_EQ(alea_find_surface_labels_on_boundary_map(map, 2, &labels,
                                                        &label_count), 0);
    ASSERT_EQ(label_count, 1);
    ASSERT_EQ(labels[0].id, 2);
    free(labels);
    alea_slice_surface_boundary_map_free(map);
    free(label_cell_ids);
    free(label_material_ids);
    alea_destroy(sys);
}

TEST(surface_boundary_map_keeps_coincident_ids_in_one_crossing_group) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);
    int owner = alea_plane_surface(sys, 1, 1.0, 0.0, 0.0, 0.0);
    int coincident = alea_plane_surface(sys, 2, 1.0, 0.0, 0.0, 0.0);
    ASSERT(owner >= 0 && coincident >= 0);
    int material = alea_add_material(sys, 1);
    ASSERT(material >= 0);
    alea_node_id_t owner_neg = alea_halfspace(sys, owner, -1);
    alea_node_id_t owner_pos = alea_halfspace(sys, owner, 1);
    ASSERT(owner_neg != ALEA_NODE_ID_INVALID && owner_pos != ALEA_NODE_ID_INVALID);
    ASSERT(alea_add_cell(sys, 1, owner_neg, material, 1.0, 0) >= 0);
    ASSERT(alea_add_cell(sys, 2, owner_pos, material, 1.0, 0) >= 0);
    ASSERT_EQ(alea_prepare_query_acceleration(sys), 0);

    const int width = 2, height = 1;
    int cell_ids[width * height];
    alea_slice_view_t view;
    alea_slice_view_axis(&view, 2, 0.0, -1.0, 1.0, -0.5, 0.5);
    ASSERT_EQ(alea_find_cells_grid(sys, &view, width, height, -1,
                                   cell_ids, NULL, NULL), 0);

    alea_slice_surface_boundary_map_t* map = NULL;
    ASSERT_EQ(alea_slice_surface_boundary_map_create(
                  sys, &view, width, height, cell_ids,
                  alea_slice_classify_cell, NULL, &map), 0);
    ASSERT_NOT_NULL(map);
    /* Both IDs are retained at one geometric crossing. They must not turn the
     * edge into MULTI_HIT, which is reserved for distinct positions. */
    ASSERT_EQ(alea_slice_surface_boundary_status(
                  map, 0, 0, ALEA_SLICE_EDGE_RIGHT),
              ALEA_SLICE_BOUNDARY_VALID);
    ASSERT_EQ(alea_slice_surface_boundary_surface_count(
                  map, 0, 0, ALEA_SLICE_EDGE_RIGHT), 2);
    ASSERT_EQ(alea_slice_surface_boundary_surface_id(
                  map, 0, 0, ALEA_SLICE_EDGE_RIGHT, 0), 1);
    ASSERT_EQ(alea_slice_surface_boundary_surface_id(
                  map, 0, 0, ALEA_SLICE_EDGE_RIGHT, 1), 2);
    ASSERT_EQ(alea_slice_surface_boundary_group_count(
                  map, 0, 0, ALEA_SLICE_EDGE_RIGHT), 1);
    ASSERT(alea_slice_surface_boundary_group_fraction(
               map, 0, 0, ALEA_SLICE_EDGE_RIGHT, 0) > 0.49 &&
           alea_slice_surface_boundary_group_fraction(
               map, 0, 0, ALEA_SLICE_EDGE_RIGHT, 0) < 0.51);
    ASSERT_EQ(alea_slice_surface_boundary_group_surface_count(
                  map, 0, 0, ALEA_SLICE_EDGE_RIGHT, 0), 2);
    ASSERT_EQ(alea_slice_surface_boundary_group_surface_id(
                  map, 0, 0, ALEA_SLICE_EDGE_RIGHT, 0, 0), 1);
    ASSERT_EQ(alea_slice_surface_boundary_group_surface_id(
                  map, 0, 0, ALEA_SLICE_EDGE_RIGHT, 0, 1), 2);

    alea_slice_surface_boundary_map_free(map);
    alea_destroy(sys);
}

TEST_MAIN()
