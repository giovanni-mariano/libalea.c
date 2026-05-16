// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/*
 * test_extract_region.c
 *
 * Tests for alea_extract_region — the region-restricted flatten path.
 * Covers:
 *   - parity with a "no-clip" bbox on a single-universe model
 *   - two-level fill with an offset transform (the world-vs-local bbox bug
 *     the previous flat sweep had)
 *   - rectangular lattice expansion + clip
 *   - rectangular lattice expansion without clip
 */

#include "alea_test.h"
#include "alea.h"
#include "core/alea_system.h"
#include "core/alea_cell.h"
#include "core/alea_universe.h"
#include <math.h>
#include <stdlib.h>

/* Helper: terminal material at a world-frame point. */
static int terminal_material(alea_system_t* sys, double x, double y, double z) {
    alea_cell_hit_t hits[16];
    int n = alea_find_all_cells(sys, x, y, z, hits, 16);
    if (n <= 0) return 0;
    return hits[n - 1].material_id;
}

/* ------------------------------------------------------------------- */
/* Single-universe parity: clip box covers everything, all cells       */
/* survive and materials at sample points match the source.            */
/* ------------------------------------------------------------------- */

TEST(extract_region_parity_no_clip) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);

    int s1 = alea_sphere_surface(sys, 0, -3, 0, 0, 1.0);
    int s2 = alea_sphere_surface(sys, 0,  3, 0, 0, 1.0);
    int sb = alea_box_surface(sys,    0, -10, 10, -10, 10, -10, 10);

    alea_node_id_t sphere1 = alea_halfspace(sys, s1, -1);
    alea_node_id_t sphere2 = alea_halfspace(sys, s2, -1);
    alea_node_id_t box     = alea_halfspace(sys, sb, -1);

    int m1 = alea_add_material(sys, 1);
    int m2 = alea_add_material(sys, 2);

    ASSERT(alea_add_cell(sys, 10, sphere1, m1, -1.0, 0) >= 0);
    ASSERT(alea_add_cell(sys, 11, sphere2, m2, -1.0, 0) >= 0);
    /* "Outside" cell — the box with the spheres carved out is messy to
     * build via complements here; just include a generous box and let the
     * point query land on a sphere cell. */
    (void)box;

    ASSERT_EQ(alea_build_universe_index(sys), 0);

    alea_bbox_t big = {-100, 100, -100, 100, -100, 100};
    alea_system_t* out = alea_extract_region(sys, &big);
    ASSERT_NOT_NULL(out);

    ASSERT_EQ(alea_build_universe_index(out), 0);

    /* Both spheres must survive. */
    ASSERT_EQ(terminal_material(out, -3, 0, 0), 1);
    ASSERT_EQ(terminal_material(out,  3, 0, 0), 2);

    alea_destroy(out);
    alea_destroy(sys);
}

/* ------------------------------------------------------------------- */
/* Two-level fill with non-identity offset transform.                  */
/*                                                                       */
/* Universe 0: cell 1 = big box, FILL=1 with translation TR=(10,0,0).  */
/* Universe 1: cell 10 = unit sphere at origin (material 1).           */
/*                                                                       */
/* In world coords, the sphere is at (10,0,0). The previous flat-sweep */
/* extract_region was wrong here because cell 10's stored bbox is in   */
/* universe-1 local coords, so a query at the world origin would       */
/* falsely include it. We confirm:                                     */
/*   - bbox around (10,0,0) → sphere material reachable                */
/*   - bbox around the origin → output has no cell with that material  */
/* ------------------------------------------------------------------- */

TEST(extract_region_two_level_fill_offset) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);

    /* Universe 1: sphere at origin (in universe-1 local frame), R=1.5 */
    int ss = alea_sphere_surface(sys, 0, 0, 0, 0, 1.5);
    alea_node_id_t sphere = alea_halfspace(sys, ss, -1);

    /* Universe 0: large box, FILL=1, translation transform. */
    int sb = alea_box_surface(sys, 0, -50, 50, -50, 50, -50, 50);
    alea_node_id_t box = alea_halfspace(sys, sb, -1);

    int m1 = alea_add_material(sys, 1);

    int sphere_idx = alea_add_cell(sys, 10, sphere, m1, -1.0, 1);
    ASSERT(sphere_idx >= 0);

    /* TR1: translate (10, 0, 0). value_count=3 (pure translation). */
    double tr_data[3] = {10.0, 0.0, 0.0};
    ASSERT_EQ(alea_add_transform(sys, 1, tr_data, 3, 0), 0);

    int container_idx = alea_add_cell(sys, 1, box, ALEA_MATERIAL_VOID, 0.0, 0);
    ASSERT(container_idx >= 0);
    ASSERT_EQ(alea_set_cell_fill(sys, container_idx, 1, 1), 0);

    ASSERT_EQ(alea_build_universe_index(sys), 0);

    /* Sanity: hierarchical traversal finds the sphere at (10, 0, 0). */
    ASSERT_EQ(terminal_material(sys, 10, 0, 0), 1);

    /* Extract a region centered on the sphere's world placement. */
    alea_bbox_t near_sphere = {8, 12, -2, 2, -2, 2};
    alea_system_t* hit = alea_extract_region(sys, &near_sphere);
    ASSERT_NOT_NULL(hit);
    ASSERT_EQ(alea_build_universe_index(hit), 0);
    ASSERT_EQ(terminal_material(hit, 10, 0, 0), 1);
    alea_destroy(hit);

    /* Extract a region centered on the world origin — the sphere is NOT
     * here (it's at world (10,0,0)). The output must not include the
     * sphere region's terminal material at the origin. The container's
     * fill universe still gets visited (the container bbox spans the
     * origin), but the sphere element is culled because its world bbox
     * is centered on (10,0,0). */
    alea_bbox_t at_origin = {-1, 1, -1, 1, -1, 1};
    alea_system_t* miss = alea_extract_region(sys, &at_origin);
    ASSERT_NOT_NULL(miss);
    ASSERT_EQ(alea_build_universe_index(miss), 0);

    /* The sphere material at the world origin must not be reachable —
     * the sphere cell wasn't extracted. */
    ASSERT_EQ(terminal_material(miss, 10, 0, 0), 0);

    alea_destroy(miss);
    alea_destroy(sys);
}

/* ------------------------------------------------------------------- */
/* Helper: build a rectangular lattice cell programmatically.          */
/* Returns the cell index in `sys`. Lattice fills universe `lat_univ`  */
/* in the parent universe `parent_univ`, with the given dims/pitch and */
/* element fills.                                                       */
/* ------------------------------------------------------------------- */

static int add_rect_lattice_cell(alea_system_t* sys,
                                 int cell_id,
                                 int parent_universe,
                                 const int dims[6],
                                 const double pitch[3],
                                 const double lower_left[3],
                                 const int* fills,
                                 size_t fill_count) {
    /* CSG region: lattice bounding box. */
    int ni = dims[1] - dims[0] + 1;
    int nj = dims[3] - dims[2] + 1;
    int nk = dims[5] - dims[4] + 1;
    int sb = alea_box_surface(sys, 0,
                              lower_left[0], lower_left[0] + ni * pitch[0],
                              lower_left[1], lower_left[1] + nj * pitch[1],
                              lower_left[2], lower_left[2] + nk * pitch[2]);
    alea_node_id_t box = alea_halfspace(sys, sb, -1);

    int idx = alea_add_cell(sys, cell_id, box, ALEA_MATERIAL_VOID,
                            0.0, parent_universe);
    if (idx < 0) return -1;

    alea_cell_entry_t* c = &sys->cells.data[idx];
    c->lat_type = 1;
    c->lat_fill_dims[0] = dims[0]; c->lat_fill_dims[1] = dims[1];
    c->lat_fill_dims[2] = dims[2]; c->lat_fill_dims[3] = dims[3];
    c->lat_fill_dims[4] = dims[4]; c->lat_fill_dims[5] = dims[5];
    c->lat_pitch[0] = pitch[0];
    c->lat_pitch[1] = pitch[1];
    c->lat_pitch[2] = pitch[2];
    c->lat_lower_left[0] = lower_left[0];
    c->lat_lower_left[1] = lower_left[1];
    c->lat_lower_left[2] = lower_left[2];
    c->lat_fill = malloc(fill_count * sizeof(int));
    if (!c->lat_fill) return -1;
    memcpy(c->lat_fill, fills, fill_count * sizeof(int));
    c->lat_fill_count = fill_count;
    return idx;
}

/* ------------------------------------------------------------------- */
/* 5x1x1 rectangular lattice. Each element is a unit-cube universe     */
/* with a single sphere of a distinct material.                        */
/*                                                                       */
/* Without clip: all 5 spheres survive. With a clip box around element */
/* index 2 only: only that sphere survives.                             */
/* ------------------------------------------------------------------- */

TEST(extract_region_rect_lattice_no_clip) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);

    /* Five element universes (1..5), each a sphere of radius 0.3 at
     * element-local origin with a distinct material. */
    int materials[5];
    for (int u = 0; u < 5; u++) {
        materials[u] = alea_add_material(sys, 100 + u);
        int sid = alea_sphere_surface(sys, 0, 0, 0, 0, 0.3);
        alea_node_id_t s = alea_halfspace(sys, sid, -1);
        ASSERT(alea_add_cell(sys, 100 + u, s, materials[u], -1.0, 1 + u) >= 0);
    }

    int dims[6] = {0, 4, 0, 0, 0, 0};
    double pitch[3] = {1.0, 1.0, 1.0};
    double lower_left[3] = {0.0, -0.5, -0.5};
    int fills[5] = {1, 2, 3, 4, 5};
    int lat_idx = add_rect_lattice_cell(sys, 1, 0, dims, pitch, lower_left,
                                        fills, 5);
    ASSERT(lat_idx >= 0);

    ASSERT_EQ(alea_build_universe_index(sys), 0);

    /* No-clip-equivalent bbox covering the whole lattice. */
    alea_bbox_t big = {-100, 100, -100, 100, -100, 100};
    alea_system_t* out = alea_extract_region(sys, &big);
    ASSERT_NOT_NULL(out);
    ASSERT_EQ(alea_build_universe_index(out), 0);

    /* All 5 sphere centers (at world (0.5,0,0), (1.5,0,0), ...) must hit
     * the corresponding material. */
    for (int u = 0; u < 5; u++) {
        double cx = lower_left[0] + (u + 0.5) * pitch[0];
        int got = terminal_material(out, cx, 0.0, 0.0);
        char msg[64];
        snprintf(msg, sizeof(msg), "element %d", u);
        ASSERT_MSG(got == 100 + u, msg);
    }

    alea_destroy(out);
    alea_destroy(sys);
}

TEST(extract_region_rect_lattice_clipped) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);

    int materials[5];
    for (int u = 0; u < 5; u++) {
        materials[u] = alea_add_material(sys, 100 + u);
        int sid = alea_sphere_surface(sys, 0, 0, 0, 0, 0.3);
        alea_node_id_t s = alea_halfspace(sys, sid, -1);
        ASSERT(alea_add_cell(sys, 100 + u, s, materials[u], -1.0, 1 + u) >= 0);
    }

    int dims[6] = {0, 4, 0, 0, 0, 0};
    double pitch[3] = {1.0, 1.0, 1.0};
    double lower_left[3] = {0.0, -0.5, -0.5};
    int fills[5] = {1, 2, 3, 4, 5};
    int lat_idx = add_rect_lattice_cell(sys, 1, 0, dims, pitch, lower_left,
                                        fills, 5);
    ASSERT(lat_idx >= 0);
    ASSERT_EQ(alea_build_universe_index(sys), 0);

    /* Clip a thin slab over element index 2 only (center at x=2.5).
     * Element 2's bounds in x are [2.0, 3.0]; clip [2.3, 2.7]. */
    alea_bbox_t clip = {2.3, 2.7, -0.4, 0.4, -0.4, 0.4};
    alea_system_t* out = alea_extract_region(sys, &clip);
    ASSERT_NOT_NULL(out);
    ASSERT_EQ(alea_build_universe_index(out), 0);

    /* Element 2's material is reachable at its center. */
    ASSERT_EQ(terminal_material(out, 2.5, 0.0, 0.0), 102);

    /* Elements 0 and 4 must NOT be reachable in the extracted system
     * (their centers shouldn't have a sphere cell anymore). */
    ASSERT_EQ(terminal_material(out, 0.5, 0.0, 0.0), 0);
    ASSERT_EQ(terminal_material(out, 4.5, 0.0, 0.0), 0);

    alea_destroy(out);
    alea_destroy(sys);
}

TEST_MAIN()
