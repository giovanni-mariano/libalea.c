// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/* Same-universe component insertion primitive: old cells are carved before
 * the incoming component cells are added. */

#include "alea_test.h"
#include "alea.h"
#include "alea_mcnp.h"

static int material_at(alea_system_t* sys, double x, double y, double z) {
    alea_cell_hit_t hits[8];
    int count = alea_find_all_cells(sys, x, y, z, hits, 8);
    return count > 0 ? hits[count - 1].material_id : 0;
}

TEST(carve_universe_preserves_space_for_component) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);

    int host_mat = alea_add_material(sys, 1);
    int component_mat = alea_add_material(sys, 2);
    int host_surface = alea_box_surface(sys, 0, -5, 5, -5, 5, -5, 5);
    int component_surface = alea_sphere_surface(sys, 0, 0, 0, 0, 1);
    alea_node_id_t host = alea_halfspace(sys, host_surface, -1);
    alea_node_id_t component = alea_halfspace(sys, component_surface, -1);

    ASSERT(alea_add_cell(sys, 1, host, host_mat, -1.0, 0) >= 0);

    int modified = -1, removed = -1;
    ASSERT_EQ(alea_carve_universe(sys, 0, component, 1, -1,
                                  &modified, &removed), 0);
    ASSERT_EQ(modified, 1);
    ASSERT_EQ(removed, 0);
    ASSERT_EQ(material_at(sys, 0, 0, 0), 0);
    ASSERT_EQ(material_at(sys, 2, 0, 0), 1);

    ASSERT(alea_add_cell(sys, 2, component, component_mat, -2.0, 0) >= 0);
    ASSERT_EQ(alea_build_universe_index(sys), 0);
    ASSERT_EQ(material_at(sys, 0, 0, 0), 2);
    ASSERT_EQ(material_at(sys, 2, 0, 0), 1);

    alea_destroy(sys);
}

TEST(carve_universe_removes_empty_cell) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);

    int mat = alea_add_material(sys, 1);
    int surface = alea_sphere_surface(sys, 0, 0, 0, 0, 1);
    alea_node_id_t sphere = alea_halfspace(sys, surface, -1);
    ASSERT(alea_add_cell(sys, 10, sphere, mat, -1.0, 0) >= 0);

    int modified = -1, removed = -1;
    ASSERT_EQ(alea_carve_universe(sys, 0, sphere, 1, -1,
                                  &modified, &removed), 0);
    ASSERT_EQ(modified, 1);
    ASSERT_EQ(removed, 1);
    ASSERT_EQ(alea_cell_count(sys), 0);
    ASSERT_EQ(alea_cell_find(sys, 10), -1);

    alea_destroy(sys);
}

TEST(carve_universe_skips_disjoint_cells) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);

    int mat = alea_add_material(sys, 1);
    int near_surface = alea_sphere_surface(sys, 0, 0, 0, 0, 2);
    int far_surface = alea_sphere_surface(sys, 0, 100, 0, 0, 2);
    int carve_surface = alea_sphere_surface(sys, 0, 0, 0, 0, 1);
    alea_node_id_t near_root = alea_halfspace(sys, near_surface, -1);
    alea_node_id_t far_root = alea_halfspace(sys, far_surface, -1);
    alea_node_id_t carve_root = alea_halfspace(sys, carve_surface, -1);
    ASSERT(alea_add_cell(sys, 1, near_root, mat, -1.0, 0) >= 0);
    ASSERT(alea_add_cell(sys, 2, far_root, mat, -1.0, 0) >= 0);

    int modified = -1, removed = -1;
    ASSERT_EQ(alea_carve_universe(sys, 0, carve_root, 1, -1,
                                  &modified, &removed), 0);
    ASSERT_EQ(modified, 1);
    ASSERT_EQ(removed, 0);

    alea_cell_info_t far_info;
    ASSERT_EQ(alea_cell_find_info(sys, 2, &far_info), 0);
    ASSERT_EQ(far_info.root, far_root);
    alea_destroy(sys);
}

TEST(mcnp_params_follow_cell_compaction) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);
    int mat = alea_add_material(sys, 1);
    int first_surface = alea_sphere_surface(sys, 0, 0, 0, 0, 1);
    int second_surface = alea_sphere_surface(sys, 0, 4, 0, 0, 1);
    ASSERT(alea_add_cell(sys, 1, alea_halfspace(sys, first_surface, -1),
                         mat, -1.0, 0) >= 0);
    ASSERT(alea_add_cell(sys, 2, alea_halfspace(sys, second_surface, -1),
                         mat, -1.0, 0) >= 0);

    mcnp_model_t* model = mcnp_model_wrap(sys);
    ASSERT_NOT_NULL(model);
    mcnp_cell_params(model, 1)->imp_n = 7.0;
    ASSERT_EQ(alea_cell_remove(sys, 0), 0);
    ASSERT_EQ(model->cell_params_count, 1);
    ASSERT_NEAR(mcnp_cell_params(model, 0)->imp_n, 7.0, 1e-12);

    mcnp_model_destroy(model);
    alea_destroy(sys);
}

TEST_MAIN()
