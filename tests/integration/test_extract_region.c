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
#include "alea_mcnp.h"
#include "core/alea_system.h"
#include "core/alea_cell.h"
#include "core/alea_universe.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

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

/* ------------------------------------------------------------------- */
/* Surface integrity: every primitive type in the extracted system     */
/* must be a valid enum value. Catches memory corruption / stale       */
/* primitive_ids in surface entries.                                    */
/* ------------------------------------------------------------------- */

/* Nested lattices: outer 3x1x1 lattice of inner 2x1x1 lattices of spheres.
 * Exercises the recursive lattice expansion path that user reports breaks
 * with very-high "Unknown primitive type" numbers from the MCNP exporter.
 *
 * Every surface in the extracted system must have a primitive type in
 * range [PLANE=1, ARB=23] and primitive_id < primitive_count. */
TEST(extract_region_nested_lattices_surface_integrity) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);

    /* Universe 20: a sphere at origin (the leaf universe). */
    int sphere_surf = alea_sphere_surface(sys, 0, 0, 0, 0, 0.15);
    alea_node_id_t sphere = alea_halfspace(sys, sphere_surf, -1);
    int m_inner = alea_add_material(sys, 1);
    ASSERT(alea_add_cell(sys, 200, sphere, m_inner, -1.0, 20) >= 0);

    /* Universe 10: inner lattice cell, 2x1x1 of universe 20. */
    int inner_dims[6] = {0, 1, 0, 0, 0, 0};
    double inner_pitch[3] = {0.5, 1.0, 1.0};
    double inner_ll[3] = {0.0, -0.5, -0.5};
    int inner_fills[2] = {20, 20};
    /* Inner lattice cell has its CSG shell = the lattice bbox. */
    int inner_sb = alea_box_surface(sys, 0,
        inner_ll[0], inner_ll[0] + 2 * inner_pitch[0],
        inner_ll[1], inner_ll[1] + 1 * inner_pitch[1],
        inner_ll[2], inner_ll[2] + 1 * inner_pitch[2]);
    alea_node_id_t inner_box = alea_halfspace(sys, inner_sb, -1);
    int inner_lat_idx = alea_add_cell(sys, 100, inner_box,
                                      ALEA_MATERIAL_VOID, 0.0, 10);
    ASSERT(inner_lat_idx >= 0);
    alea_cell_entry_t* inner_lc = &sys->cells.data[inner_lat_idx];
    inner_lc->lat_type = 1;
    inner_lc->lat_fill_dims[0] = inner_dims[0]; inner_lc->lat_fill_dims[1] = inner_dims[1];
    inner_lc->lat_fill_dims[2] = inner_dims[2]; inner_lc->lat_fill_dims[3] = inner_dims[3];
    inner_lc->lat_fill_dims[4] = inner_dims[4]; inner_lc->lat_fill_dims[5] = inner_dims[5];
    inner_lc->lat_pitch[0] = inner_pitch[0];
    inner_lc->lat_pitch[1] = inner_pitch[1];
    inner_lc->lat_pitch[2] = inner_pitch[2];
    inner_lc->lat_lower_left[0] = inner_ll[0];
    inner_lc->lat_lower_left[1] = inner_ll[1];
    inner_lc->lat_lower_left[2] = inner_ll[2];
    inner_lc->lat_fill = malloc(2 * sizeof(int));
    ASSERT_NOT_NULL(inner_lc->lat_fill);
    memcpy(inner_lc->lat_fill, inner_fills, 2 * sizeof(int));
    inner_lc->lat_fill_count = 2;

    /* Universe 0: outer lattice 3x1x1 of universe 10. */
    int outer_dims[6] = {0, 2, 0, 0, 0, 0};
    double outer_pitch[3] = {2.0, 1.0, 1.0};
    double outer_ll[3] = {0.0, -0.5, -0.5};
    int outer_fills[3] = {10, 10, 10};
    int outer_sb = alea_box_surface(sys, 0,
        outer_ll[0], outer_ll[0] + 3 * outer_pitch[0],
        outer_ll[1], outer_ll[1] + 1 * outer_pitch[1],
        outer_ll[2], outer_ll[2] + 1 * outer_pitch[2]);
    alea_node_id_t outer_box = alea_halfspace(sys, outer_sb, -1);
    int outer_lat_idx = alea_add_cell(sys, 1, outer_box,
                                      ALEA_MATERIAL_VOID, 0.0, 0);
    ASSERT(outer_lat_idx >= 0);
    alea_cell_entry_t* outer_lc = &sys->cells.data[outer_lat_idx];
    outer_lc->lat_type = 1;
    outer_lc->lat_fill_dims[0] = outer_dims[0]; outer_lc->lat_fill_dims[1] = outer_dims[1];
    outer_lc->lat_fill_dims[2] = outer_dims[2]; outer_lc->lat_fill_dims[3] = outer_dims[3];
    outer_lc->lat_fill_dims[4] = outer_dims[4]; outer_lc->lat_fill_dims[5] = outer_dims[5];
    outer_lc->lat_pitch[0] = outer_pitch[0];
    outer_lc->lat_pitch[1] = outer_pitch[1];
    outer_lc->lat_pitch[2] = outer_pitch[2];
    outer_lc->lat_lower_left[0] = outer_ll[0];
    outer_lc->lat_lower_left[1] = outer_ll[1];
    outer_lc->lat_lower_left[2] = outer_ll[2];
    outer_lc->lat_fill = malloc(3 * sizeof(int));
    ASSERT_NOT_NULL(outer_lc->lat_fill);
    memcpy(outer_lc->lat_fill, outer_fills, 3 * sizeof(int));
    outer_lc->lat_fill_count = 3;

    ASSERT_EQ(alea_build_universe_index(sys), 0);

    /* Extract a region covering everything. */
    alea_bbox_t big = {-100, 100, -100, 100, -100, 100};
    alea_system_t* out = alea_extract_region(sys, &big);
    ASSERT_NOT_NULL(out);

    /* Validate every surface's primitive type is in range. */
    size_t nsurf = alea_surface_count(out);
    ASSERT(nsurf > 0);
    for (size_t i = 0; i < nsurf; i++) {
        int surf_id;
        alea_primitive_type_t type;
        alea_node_id_t pn, nn;
        alea_boundary_type_t bnd;
        ASSERT_EQ(alea_surface_get(out, i, &surf_id, &type, &pn, &nn, &bnd), 0);
        char msg[128];
        snprintf(msg, sizeof(msg),
                 "nested-lat surface idx=%zu id=%d has out-of-range type=%d (cell_count=%zu)",
                 i, surf_id, (int)type, alea_cell_count(out));
        ASSERT_MSG((int)type >= ALEA_PRIMITIVE_PLANE && (int)type <= ALEA_PRIMITIVE_ARB,
                   msg);
    }

    /* Export to MCNP. The export pass mutates the system (assigns missing
     * surface IDs, expands macrobodies) and is where the user reports
     * "Unknown primitive type %d" with very high numbers. Re-validate
     * types after export to catch any corruption introduced there. */
    FILE* tmpf = tmpfile();
    ASSERT_NOT_NULL(tmpf);
    ASSERT_EQ(mcnp_export_system_stream(out, tmpf), 0);
    fclose(tmpf);

    nsurf = alea_surface_count(out);
    for (size_t i = 0; i < nsurf; i++) {
        int surf_id;
        alea_primitive_type_t type;
        alea_node_id_t pn, nn;
        alea_boundary_type_t bnd;
        ASSERT_EQ(alea_surface_get(out, i, &surf_id, &type, &pn, &nn, &bnd), 0);
        char msg[160];
        snprintf(msg, sizeof(msg),
                 "post-export surface idx=%zu id=%d has out-of-range type=%d",
                 i, surf_id, (int)type);
        ASSERT_MSG((int)type >= ALEA_PRIMITIVE_PLANE && (int)type <= ALEA_PRIMITIVE_ARB,
                   msg);
    }

    alea_destroy(out);
    alea_destroy(sys);
}

TEST(extract_region_surface_types_in_range) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);

    /* Build a small universe with a mix of surface types. */
    int ss = alea_sphere_surface(sys, 0, 0, 0, 0, 1.5);
    int sp = alea_plane_surface(sys, 0, 1, 0, 0, -0.5);
    int sc = alea_cylinder_z_surface(sys, 0, 0, 0, 1.0);
    int sb = alea_box_surface(sys, 0, -5, 5, -5, 5, -5, 5);

    alea_node_id_t sphere   = alea_halfspace(sys, ss, -1);
    alea_node_id_t plane    = alea_halfspace(sys, sp, -1);
    alea_node_id_t cylinder = alea_halfspace(sys, sc, -1);
    alea_node_id_t box      = alea_halfspace(sys, sb, -1);

    int m1 = alea_add_material(sys, 1);
    int m2 = alea_add_material(sys, 2);
    int m3 = alea_add_material(sys, 3);

    /* Universe 1: three cells using different surface types. */
    ASSERT(alea_add_cell(sys, 100, sphere,   m1, -1.0, 1) >= 0);
    ASSERT(alea_add_cell(sys, 101, plane,    m2, -1.0, 1) >= 0);
    ASSERT(alea_add_cell(sys, 102, cylinder, m3, -1.0, 1) >= 0);

    /* TR1: offset universe 1 to (4, 0, 0) so cloning happens through a
     * transform (this is where mc_surface_id gets zeroed). */
    double tr_data[3] = {4.0, 0.0, 0.0};
    ASSERT_EQ(alea_add_transform(sys, 1, tr_data, 3, 0), 0);

    int container_idx = alea_add_cell(sys, 1, box, ALEA_MATERIAL_VOID, 0.0, 0);
    ASSERT(container_idx >= 0);
    ASSERT_EQ(alea_set_cell_fill(sys, container_idx, 1, 1), 0);

    ASSERT_EQ(alea_build_universe_index(sys), 0);

    /* Extract a region covering the placed universe. */
    alea_bbox_t clip = {2, 6, -2, 2, -2, 2};
    alea_system_t* out = alea_extract_region(sys, &clip);
    ASSERT_NOT_NULL(out);

    /* Every surface in the output must reference a primitive whose type
     * is in the valid enum range [PLANE=1, ARB=23]. A "very high" number
     * indicates a stale or out-of-bounds primitive_id. */
    size_t nsurf = alea_surface_count(out);
    for (size_t i = 0; i < nsurf; i++) {
        int surf_id;
        alea_primitive_type_t type;
        alea_node_id_t pn, nn;
        alea_boundary_type_t bnd;
        ASSERT_EQ(alea_surface_get(out, i, &surf_id, &type, &pn, &nn, &bnd), 0);
        char msg[96];
        snprintf(msg, sizeof(msg),
                 "surface idx=%zu id=%d has out-of-range type=%d",
                 i, surf_id, (int)type);
        ASSERT_MSG((int)type >= ALEA_PRIMITIVE_PLANE && (int)type <= ALEA_PRIMITIVE_ARB,
                   msg);
    }

    alea_destroy(out);
    alea_destroy(sys);
}

TEST_MAIN()
