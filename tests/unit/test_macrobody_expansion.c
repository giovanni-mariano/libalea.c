// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
// SPDX-License-Identifier: MPL-2.0
#include "alea_test.h"
#include "alea.h"
#include "core/alea_system.h"

TEST(expansion_survives_node_vector_growth) {
    for (int complement = 0; complement < 2; complement++) {
        alea_system_t* sys = alea_create();
        ASSERT_NOT_NULL(sys);
        int box = alea_box_surface(sys, 1, -1, 1, -1, 1, -1, 1);
        int sphere = alea_sphere_surface(sys, 2, 10, 0, 0, 1);
        ASSERT(box >= 0 && sphere >= 0);
        alea_node_id_t root = alea_union(sys, alea_halfspace(sys, box, -1),
                                       alea_halfspace(sys, sphere, -1));
        if (complement) root = alea_complement(sys, root);
        ASSERT_NE(root, ALEA_NODE_ID_INVALID);
        while (sys->nodes.count + 2 < sys->nodes.capacity) {
            ASSERT(alea_sphere_surface(sys, 0, 30 + (double)sys->nodes.count, 0, 0, 1) >= 0);
        }
        size_t capacity = sys->nodes.capacity;
        alea_node_id_t expanded = alea_expand_all_macrobodies(sys, root);
        ASSERT_NE(expanded, ALEA_NODE_ID_INVALID);
        ASSERT(sys->nodes.capacity > capacity);
        ASSERT_EQ(alea_point_inside(sys, expanded, 0, 0, 0), !complement);
        ASSERT_EQ(alea_point_inside(sys, expanded, 10, 0, 0), !complement);
        ASSERT_EQ(alea_point_inside(sys, expanded, 5, 0, 0), complement);
        alea_destroy(sys);
    }
}

TEST(expansion_propagates_invalid_child_before_bbox_access) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);
    int sphere = alea_sphere_surface(sys, 1, 0, 0, 0, 1);
    ASSERT(sphere >= 0);
    alea_node_id_t child = alea_halfspace(sys, sphere, -1);
    alea_node_id_t root = alea_complement(sys, child);
    sys->nodes.data[root].operation.left = ALEA_NODE_ID_INVALID;
    ASSERT_EQ(alea_expand_all_macrobodies(sys, root), ALEA_NODE_ID_INVALID);
    root = alea_union(sys, child, alea_halfspace(sys, sphere, 1));
    sys->nodes.data[root].operation.right = ALEA_NODE_ID_INVALID;
    ASSERT_EQ(alea_expand_all_macrobodies(sys, root), ALEA_NODE_ID_INVALID);
    sys->nodes.data[root].operation.left = ALEA_NODE_ID_INVALID;
    ASSERT_EQ(alea_expand_all_macrobodies(sys, root), ALEA_NODE_ID_INVALID);
    alea_destroy(sys);
}

TEST_MAIN()
