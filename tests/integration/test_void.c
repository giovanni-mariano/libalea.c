// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

#define _POSIX_C_SOURCE 199309L

/*
 * test_void.c - Void generation integration tests using the public alea API
 */

#include "alea_test.h"
#include "alea.h"
#include "core/alea_void.h"
#include "core/alea_system.h"
#include "alea_mcnp.h"
#include <stdlib.h>
#include <string.h>
#include <time.h>

/* Reduce octree depth for fast tests (default 8 is too slow for CI) */
static void set_fast_void_config(alea_system_t* sys) {
    alea_config_t cfg = alea_get_config(sys);
    cfg.void_max_depth = 4;
    cfg.void_min_size = 0.5;
    alea_set_config(sys, &cfg);
}

/* Configure a deliberately coarse octree for tests that exercise probe
 * classification at the root node. */
static void set_coarse_void_config(alea_system_t* sys) {
    alea_config_t cfg = alea_get_config(sys);
    cfg.void_max_depth = 4;
    cfg.void_min_size = 0.1;
    cfg.void_probes_per_axis = 3;
    alea_set_config(sys, &cfg);
}

/* Helper: create a system with a single box cell at [-1,1]^3 in universe 0 */
static alea_system_t* make_box_system(void) {
    alea_system_t* sys = alea_create();
    if (!sys) return NULL;

    set_fast_void_config(sys);

    int si = alea_box_surface(sys, 0, -1, 1, -1, 1, -1, 1);
    alea_node_id_t box = alea_halfspace(sys, si, -1);

    int m1 = alea_add_material(sys, 1);
    alea_add_cell(sys, 1, box, m1, 1.0, 0);

    alea_build_universe_index(sys);
    return sys;
}

/* Helper: create a material box with a small internal void cavity that does
 * not lie on the default 3x3x3 probe grid for bounds [-4,4]^3. */
static alea_system_t* make_box_shell_with_off_grid_cavity(void) {
    alea_system_t* sys = alea_create();
    if (!sys) return NULL;

    set_coarse_void_config(sys);

    int outer_si = alea_box_surface(sys, 0, -5, 5, -5, 5, -5, 5);
    int cavity_si = alea_sphere_surface(sys, 0, 1.25, 0.75, 0.25, 0.25);
    if (outer_si < 0 || cavity_si < 0) {
        alea_destroy(sys);
        return NULL;
    }

    alea_node_id_t outer = alea_halfspace(sys, outer_si, -1);
    alea_node_id_t cavity = alea_halfspace(sys, cavity_si, -1);
    alea_node_id_t shell = alea_difference(sys, outer, cavity);
    if (shell == ALEA_NODE_ID_INVALID) {
        alea_destroy(sys);
        return NULL;
    }

    int mat = alea_add_material(sys, 1);
    if (alea_add_cell(sys, 1, shell, mat, 1.0, 0) < 0 ||
        alea_build_universe_index(sys) < 0) {
        alea_destroy(sys);
        return NULL;
    }

    return sys;
}

/* Helper: root universe contains only a filled container. Void generation over
 * bounds fully inside the container should conservatively produce no root-level
 * void cells. */
static alea_system_t* make_simple_fill_system(void) {
    alea_system_t* sys = alea_create();
    if (!sys) return NULL;

    set_fast_void_config(sys);

    int child_si = alea_sphere_surface(sys, 0, 0, 0, 0, 1.0);
    int parent_si = alea_box_surface(sys, 0, -5, 5, -5, 5, -5, 5);
    if (child_si < 0 || parent_si < 0) {
        alea_destroy(sys);
        return NULL;
    }

    alea_node_id_t child = alea_halfspace(sys, child_si, -1);
    alea_node_id_t parent = alea_halfspace(sys, parent_si, -1);

    int mat = alea_add_material(sys, 1);
    if (alea_add_cell(sys, 10, child, mat, 1.0, 1) < 0) {
        alea_destroy(sys);
        return NULL;
    }

    int fill_idx = alea_add_cell(sys, 1, parent, ALEA_MATERIAL_VOID, 0.0, 0);
    if (fill_idx < 0 || alea_set_fill(sys, fill_idx, 1, 0) < 0 ||
        alea_build_universe_index(sys) < 0) {
        alea_destroy(sys);
        return NULL;
    }

    return sys;
}

/* ------------------------------------------------------------------------- */
/* Core correctness                                                           */
/* ------------------------------------------------------------------------- */

TEST(void_single_box) {
    alea_system_t* sys = make_box_system();
    ASSERT_NOT_NULL(sys);

    alea_bbox_t bounds = {-5, 5, -5, 5, -5, 5};
    void_result_t* vr = alea_void_generate(sys, &bounds);
    ASSERT_NOT_NULL(vr);

    /* There IS void outside the box */
    ASSERT(alea_void_count(vr) > 0);

    /* Each returned bbox should be valid and within search bounds */
    for (size_t i = 0; i < alea_void_count(vr); i++) {
        alea_bbox_t box;
        ASSERT_EQ(alea_void_get(vr, i, &box), 0);

        /* Positive volume */
        double vol = (box.max_x - box.min_x) *
                     (box.max_y - box.min_y) *
                     (box.max_z - box.min_z);
        ASSERT(vol > 0);

        /* Within search bounds (with tolerance for float rounding) */
        ASSERT(box.min_x >= bounds.min_x - 1e-6);
        ASSERT(box.max_x <= bounds.max_x + 1e-6);
        ASSERT(box.min_y >= bounds.min_y - 1e-6);
        ASSERT(box.max_y <= bounds.max_y + 1e-6);
        ASSERT(box.min_z >= bounds.min_z - 1e-6);
        ASSERT(box.max_z <= bounds.max_z + 1e-6);
    }

    alea_void_free(vr);
    alea_destroy(sys);
}

/* ------------------------------------------------------------------------- */
/* CSG node creation                                                          */
/* ------------------------------------------------------------------------- */

TEST(void_to_node) {
    alea_system_t* sys = make_box_system();
    ASSERT_NOT_NULL(sys);

    alea_bbox_t bounds = {-5, 5, -5, 5, -5, 5};
    void_result_t* vr = alea_void_generate(sys, &bounds);
    ASSERT_NOT_NULL(vr);
    ASSERT(alea_void_count(vr) > 0);

    alea_node_id_t node = alea_void_to_node(sys, vr);
    ASSERT_NE(node, ALEA_NODE_ID_INVALID);

    alea_void_free(vr);
    alea_destroy(sys);
}

/* ------------------------------------------------------------------------- */
/* Cell registration                                                          */
/* ------------------------------------------------------------------------- */

TEST(void_add_cells) {
    alea_system_t* sys = make_box_system();
    ASSERT_NOT_NULL(sys);

    size_t before = alea_cell_count(sys);

    alea_bbox_t bounds = {-5, 5, -5, 5, -5, 5};
    void_result_t* vr = alea_void_generate(sys, &bounds);
    ASSERT_NOT_NULL(vr);

    size_t void_count = alea_void_count(vr);
    ASSERT(void_count > 0);

    int added = alea_void_add_cells(sys, vr);
    ASSERT(added > 0);
    ASSERT_EQ((size_t)added, void_count);
    ASSERT_EQ(alea_cell_count(sys), before + void_count);

    alea_void_free(vr);
    alea_destroy(sys);
}

/* ------------------------------------------------------------------------- */
/* Merge reduces count                                                        */
/* ------------------------------------------------------------------------- */

TEST(void_merge) {
    alea_system_t* sys = make_box_system();
    ASSERT_NOT_NULL(sys);

    alea_bbox_t bounds = {-5, 5, -5, 5, -5, 5};
    void_result_t* vr = alea_void_generate(sys, &bounds);
    ASSERT_NOT_NULL(vr);

    size_t count_before = alea_void_count(vr);
    ASSERT(count_before > 0);

    int merged = alea_void_merge(sys, vr);
    ASSERT(merged > 0);

    size_t count_after = alea_void_count(vr);
    ASSERT(count_after <= count_before);
    ASSERT(count_after > 0);

    alea_void_free(vr);
    alea_destroy(sys);
}

/* ------------------------------------------------------------------------- */
/* No void when cell covers entire bounds                                     */
/* ------------------------------------------------------------------------- */

TEST(void_no_void) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);
    set_fast_void_config(sys);

    /* Box [-10,10]^3 fully encloses bounds [-5,5]^3 */
    int si = alea_box_surface(sys, 0, -10, 10, -10, 10, -10, 10);
    alea_node_id_t box = alea_halfspace(sys, si, -1);

    int m1 = alea_add_material(sys, 1);
    alea_add_cell(sys, 1, box, m1, 1.0, 0);
    alea_build_universe_index(sys);

    alea_bbox_t bounds = {-5, 5, -5, 5, -5, 5};
    void_result_t* vr = alea_void_generate(sys, &bounds);
    ASSERT_NOT_NULL(vr);

    ASSERT_EQ(alea_void_count(vr), 0);

    alea_void_free(vr);
    alea_destroy(sys);
}

/* ------------------------------------------------------------------------- */
/* Non-box geometry (sphere)                                                  */
/* ------------------------------------------------------------------------- */

TEST(void_sphere) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);
    set_fast_void_config(sys);

    int si = alea_sphere_surface(sys, 0, 0, 0, 0, 3.0);
    alea_node_id_t s = alea_halfspace(sys, si, -1);

    int m1 = alea_add_material(sys, 1);
    alea_add_cell(sys, 1, s, m1, 1.0, 0);
    alea_build_universe_index(sys);

    alea_bbox_t bounds = {-5, 5, -5, 5, -5, 5};
    void_result_t* vr = alea_void_generate(sys, &bounds);
    ASSERT_NOT_NULL(vr);

    /* Void in corners outside sphere */
    ASSERT(alea_void_count(vr) > 0);

    /* All boxes within bounds */
    for (size_t i = 0; i < alea_void_count(vr); i++) {
        alea_bbox_t box;
        ASSERT_EQ(alea_void_get(vr, i, &box), 0);
        ASSERT(box.min_x >= bounds.min_x - 1e-6);
        ASSERT(box.max_x <= bounds.max_x + 1e-6);
        ASSERT(box.min_y >= bounds.min_y - 1e-6);
        ASSERT(box.max_y <= bounds.max_y + 1e-6);
        ASSERT(box.min_z >= bounds.min_z - 1e-6);
        ASSERT(box.max_z <= bounds.max_z + 1e-6);
    }

    alea_void_free(vr);
    alea_destroy(sys);
}

/* ------------------------------------------------------------------------- */
/* Point verification: CSG correctness after void cell registration           */
/* ------------------------------------------------------------------------- */

TEST(void_point_verification) {
    alea_system_t* sys = make_box_system();
    ASSERT_NOT_NULL(sys);

    int mat_cell = alea_find_cell(sys, 0, 0, 0);
    ASSERT(mat_cell >= 0);

    /* Before void: point outside box is unresolved */
    ASSERT_EQ(alea_find_cell(sys, 3, 3, 3), -1);

    alea_bbox_t bounds = {-5, 5, -5, 5, -5, 5};
    void_result_t* vr = alea_void_generate(sys, &bounds);
    ASSERT_NOT_NULL(vr);
    ASSERT(alea_void_count(vr) > 0);

    alea_void_add_cells(sys, vr);
    alea_build_universe_index(sys);

    /* Point inside material cell is still found */
    int found = alea_find_cell(sys, 0, 0, 0);
    ASSERT(found >= 0);

    /* Point that was void is now covered by a void cell */
    int void_cell = alea_find_cell(sys, 3, 3, 3);
    ASSERT(void_cell >= 0);

    /* Verify it's a void material cell (MCNP convention: material_id=0 for void) */
    int material_id;
    ASSERT_EQ(alea_cell_get(sys, (size_t)void_cell, NULL, &material_id,
                             NULL, NULL, NULL, NULL), 0);
    ASSERT_EQ(material_id, 0);

    alea_void_free(vr);
    alea_destroy(sys);
}

TEST(void_add_cells_invalidates_prepared_query_caches) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);
    set_fast_void_config(sys);

    int si = alea_sphere_surface(sys, 0, 0, 0, 0, 3.0);
    ASSERT(si >= 0);
    alea_node_id_t sphere = alea_halfspace(sys, si, -1);
    int m1 = alea_add_material(sys, 1);
    int sphere_cell = alea_add_cell(sys, 1, sphere, m1, 1.0, 0);
    ASSERT(sphere_cell >= 0);
    ASSERT_EQ(alea_build_universe_index(sys), 0);

    alea_bbox_t bounds = {-5, 5, -5, 5, -5, 5};
    void_result_t* vr = alea_void_generate(sys, &bounds);
    ASSERT_NOT_NULL(vr);
    ASSERT(alea_void_count(vr) > 0);

    ASSERT_EQ(alea_prepare_query_acceleration(sys), 0);
    alea_cell_hit_t hits[16];
    int count = alea_find_all_cells(sys, 0, 0, 0, hits, 16);
    ASSERT_EQ(count, 1);
    ASSERT_EQ(hits[0].cell_id, 1);

    int added = alea_void_add_cells(sys, vr);
    ASSERT(added > 0);

    ASSERT_EQ(alea_prepare_query_acceleration(sys), 0);
    ASSERT(sys->universe_index_built);

    count = alea_find_all_cells(sys, 0, 0, 0, hits, 16);
    ASSERT_EQ(count, 1);
    ASSERT_EQ(hits[0].cell_id, 1);

    count = alea_find_all_cells(sys, 4.5, 4.5, 4.5, hits, 16);
    ASSERT(count >= 1);
    ASSERT_EQ(hits[count - 1].material_id, 0);

    alea_void_free(vr);
    alea_destroy(sys);
}

TEST(void_detects_internal_cavity_missed_by_probes) {
    alea_system_t* sys = make_box_shell_with_off_grid_cavity();
    ASSERT_NOT_NULL(sys);

    /* The cavity point is not in the material shell before void generation. */
    ASSERT_EQ(alea_find_cell(sys, 1.25, 0.75, 0.25), -1);

    alea_bbox_t bounds = {-4, 4, -4, 4, -4, 4};
    void_result_t* vr = alea_void_generate(sys, &bounds);
    ASSERT_NOT_NULL(vr);
    ASSERT_MSG(alea_void_count(vr) > 0,
               "void generation must not classify the whole bounds as solid from probes alone");

    int added = alea_void_add_cells(sys, vr);
    ASSERT(added > 0);
    ASSERT_EQ(alea_build_universe_index(sys), 0);

    ASSERT(alea_find_cell(sys, 1.25, 0.75, 0.25) >= 0);
    ASSERT_EQ(alea_material_at(sys, 1.25, 0.75, 0.25), 0);
    ASSERT_EQ(alea_material_at(sys, 0, 0, 0), 1);

    alea_void_free(vr);
    alea_destroy(sys);
}

TEST(void_treats_filled_root_container_as_carveout) {
    alea_system_t* sys = make_simple_fill_system();
    ASSERT_NOT_NULL(sys);

    ASSERT_EQ(alea_material_at(sys, 0, 0, 0), 1);

    alea_bbox_t bounds = {-4, 4, -4, 4, -4, 4};
    void_result_t* vr = alea_void_generate(sys, &bounds);
    ASSERT_NOT_NULL(vr);

    ASSERT_MSG(alea_void_count(vr) == 0,
               "bounds fully inside a filled root container should be carved out, not filled with root-level void");

    alea_void_free(vr);
    alea_destroy(sys);
}

TEST(void_generate_without_add_does_not_commit_geometry) {
    alea_system_t* sys = make_box_system();
    ASSERT_NOT_NULL(sys);

    size_t surfaces_before = alea_vec_count(&sys->surfaces);
    size_t nodes_before = alea_vec_count(&sys->nodes);
    size_t cells_before = alea_cell_count(sys);

    alea_bbox_t bounds = {-5, 5, -5, 5, -5, 5};
    void_result_t* vr = alea_void_generate(sys, &bounds);
    ASSERT_NOT_NULL(vr);
    ASSERT(alea_void_count(vr) > 0);

    alea_void_free(vr);

    ASSERT_EQ(alea_cell_count(sys), cells_before);
    ASSERT_EQ(alea_vec_count(&sys->surfaces), surfaces_before);
    ASSERT_EQ(alea_vec_count(&sys->nodes), nodes_before);

    alea_destroy(sys);
}

TEST(void_add_cells_rejects_invalid_region_without_partial_commit) {
    alea_system_t* sys = make_box_system();
    ASSERT_NOT_NULL(sys);

    alea_bbox_t bounds = {-5, 5, -5, 5, -5, 5};
    void_result_t* vr = alea_void_generate(sys, &bounds);
    ASSERT_NOT_NULL(vr);

    /* Need >= 2 regions so invalidating the last entry actually exercises
     * atomicity: a non-atomic add would commit [0..n-2] before failing at n-1. */
    size_t n = alea_void_count(vr);
    ASSERT(n >= 2);

    size_t cells_before = alea_cell_count(sys);
    vr->void_regions.data[n - 1].node = ALEA_NODE_ID_INVALID;

    int added = alea_void_add_cells(sys, vr);
    ASSERT_EQ(added, -1);
    ASSERT_EQ(alea_cell_count(sys), cells_before);

    alea_void_free(vr);
    alea_destroy(sys);
}

/* Fault-injection counter: fail the Nth allocation, succeed on all others. */
typedef struct {
    int  trigger_at;  /* call index that should fail (0-based) */
    int  call_count;
} fail_at_counter_t;

static bool fail_at_nth(void* ud) {
    fail_at_counter_t* c = ud;
    return (c->call_count++ == c->trigger_at);
}

TEST(void_generation_aborts_on_partial_octree_allocation_failure) {
    /* Use a coarse octree so the partial subtree is small but real:
     * the root spawns 8 children, then each MIXED child spawns 8 more.
     * Failing somewhere mid-way exercises the partial-tree rollback path. */
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);
    set_coarse_void_config(sys);

    int outer_si  = alea_box_surface(sys, 0, -5, 5, -5, 5, -5, 5);
    int cavity_si = alea_sphere_surface(sys, 0, 1.25, 0.75, 0.25, 0.25);
    ASSERT(outer_si >= 0 && cavity_si >= 0);
    alea_node_id_t outer  = alea_halfspace(sys, outer_si, -1);
    alea_node_id_t cavity = alea_halfspace(sys, cavity_si, -1);
    alea_node_id_t shell  = alea_difference(sys, outer, cavity);
    int mat = alea_add_material(sys, 1);
    ASSERT(alea_add_cell(sys, 1, shell, mat, 1.0, 0) >= 0);
    ASSERT_EQ(alea_build_universe_index(sys), 0);

    size_t surfaces_before   = alea_vec_count(&sys->surfaces);
    size_t nodes_before      = alea_vec_count(&sys->nodes);
    size_t primitives_before = alea_vec_count(&sys->primitives);
    size_t cells_before      = alea_cell_count(sys);

    /* Fail the 5th octree-node allocation: root succeeds, several children
     * succeed, then one fails partway through subdivision. */
    fail_at_counter_t counter = { .trigger_at = 5, .call_count = 0 };
    alea_void_set_octree_alloc_failure(fail_at_nth, &counter);

    alea_bbox_t bounds = {-4, 4, -4, 4, -4, 4};
    void_result_t* vr = alea_void_generate(sys, &bounds);

    /* Disable injection before any further test work. */
    alea_void_set_octree_alloc_failure(NULL, NULL);

    /* Generation must report failure rather than return a partial result. */
    ASSERT_NULL(vr);

    /* And the system must look exactly as it did before generation —
     * any surfaces/nodes/primitives appended along the way are rolled back. */
    ASSERT_EQ(alea_cell_count(sys), cells_before);
    ASSERT_EQ(alea_vec_count(&sys->surfaces),   surfaces_before);
    ASSERT_EQ(alea_vec_count(&sys->nodes),      nodes_before);
    ASSERT_EQ(alea_vec_count(&sys->primitives), primitives_before);

    /* Sanity: a clean generate after the injection is removed should
     * succeed normally — proves rollback left sys in a usable state. */
    void_result_t* vr2 = alea_void_generate(sys, &bounds);
    ASSERT_NOT_NULL(vr2);
    alea_void_free(vr2);

    alea_destroy(sys);
}

/* ------------------------------------------------------------------------- */
/* Null safety                                                                */
/* ------------------------------------------------------------------------- */

TEST(void_null_safety) {
    /* NULL sys returns NULL */
    void_result_t* vr = alea_void_generate(NULL, NULL);
    ASSERT_NULL(vr);

    /* NULL bounds with valid sys uses auto-bounds (should not crash) */
    alea_system_t* sys = make_box_system();
    ASSERT_NOT_NULL(sys);

    vr = alea_void_generate(sys, NULL);
    ASSERT_NOT_NULL(vr);
    ASSERT(alea_void_count(vr) > 0);

    alea_void_free(vr);
    alea_destroy(sys);
}

/* ------------------------------------------------------------------------- */
/* Benchmark: face-sorted vs greedy merge                                    */
/* ------------------------------------------------------------------------- */

static double now_sec(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + ts.tv_nsec * 1e-9;
}

/* Helper: create sphere system for benchmarking */
static alea_system_t* make_sphere_system(void) {
    alea_system_t* sys = alea_create();
    if (!sys) return NULL;

    alea_config_t cfg = alea_get_config(sys);
    cfg.void_max_depth = 4;
    cfg.void_min_size = 0.5;
    alea_set_config(sys, &cfg);

    int si = alea_sphere_surface(sys, 0, 0, 0, 0, 3.0);
    alea_node_id_t s = alea_halfspace(sys, si, -1);

    int m1 = alea_add_material(sys, 1);
    alea_add_cell(sys, 1, s, m1, 1.0, 0);
    alea_build_universe_index(sys);
    return sys;
}

/* ------------------------------------------------------------------------- */
/* Per-region complement filtering                                            */
/* ------------------------------------------------------------------------- */

TEST(void_per_region_filtering) {
    /*
     * Two spheres far apart. Void regions near sphere A should NOT reference
     * sphere B's surfaces (and vice versa). With per-region filtering, each
     * void cell only includes overlapping cells, bounding complexity.
     *
     * We verify: (1) correct point classification, (2) void generation works,
     * (3) merge+add produces valid cells covering the void space.
     */
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);

    alea_config_t cfg = alea_get_config(sys);
    cfg.void_max_depth = 4;
    cfg.void_min_size = 0.5;
    alea_set_config(sys, &cfg);

    /* Sphere A: R=2 at (-10, 0, 0) */
    int s1 = alea_sphere_surface(sys, 0, -10, 0, 0, 2.0);
    alea_node_id_t n1 = alea_halfspace(sys, s1, -1);
    int m1 = alea_add_material(sys, 1);
    alea_add_cell(sys, 1, n1, m1, 1.0, 0);

    /* Sphere B: R=2 at (+10, 0, 0) */
    int s2 = alea_sphere_surface(sys, 0, 10, 0, 0, 2.0);
    alea_node_id_t n2 = alea_halfspace(sys, s2, -1);
    int m2 = alea_add_material(sys, 2);
    alea_add_cell(sys, 2, n2, m2, 1.0, 0);

    alea_build_universe_index(sys);

    /* Generate void with bounds covering both spheres and the gap */
    alea_bbox_t bounds = {-14, 14, -4, 4, -4, 4};
    void_result_t* vr = alea_void_generate(sys, &bounds);
    ASSERT_NOT_NULL(vr);
    ASSERT(alea_void_count(vr) > 0);

    /* Merge (no consolidation — we want to verify per-region cells) */
    alea_void_merge_config_t mcfg = ALEA_VOID_MERGE_DEFAULT;
    mcfg.consolidate_max_surfaces = 0;
    int merged = alea_merge_void_cells(sys, vr, &mcfg);
    ASSERT(merged > 0);

    /* Add void cells and verify point classification */
    alea_void_add_cells(sys, vr);
    alea_build_universe_index(sys);

    /* Inside sphere A = material 1 */
    ASSERT_EQ(alea_material_at(sys, -10, 0, 0), 1);

    /* Inside sphere B = material 2 */
    ASSERT_EQ(alea_material_at(sys, 10, 0, 0), 2);

    /* Midpoint between spheres = void */
    int mid_cell = alea_find_cell(sys, 0, 0, 0);
    ASSERT(mid_cell >= 0);
    ASSERT_EQ(alea_material_at(sys, 0, 0, 0), 0);

    /* Far corners = void */
    int corner_cell = alea_find_cell(sys, -13, 3, 3);
    ASSERT(corner_cell >= 0);
    ASSERT_EQ(alea_material_at(sys, -13, 3, 3), 0);

    alea_void_free(vr);
    alea_destroy(sys);
}

/* ------------------------------------------------------------------------- */
/* Void consolidation                                                         */
/* ------------------------------------------------------------------------- */

TEST(void_consolidate) {
    alea_system_t* sys = make_sphere_system();
    ASSERT_NOT_NULL(sys);

    alea_bbox_t bounds = {-5, 5, -5, 5, -5, 5};
    void_result_t* vr = alea_void_generate(sys, &bounds);
    ASSERT_NOT_NULL(vr);

    size_t before = alea_void_count(vr);
    ASSERT(before > 1);

    /* Current void generation creates box surfaces eagerly, which invalidates
     * the universe cache needed by consolidation's global-void path. */
    ASSERT_EQ(alea_build_universe_index(sys), 0);

    /* Merge with consolidation (default) */
    alea_void_merge_config_t cfg = ALEA_VOID_MERGE_DEFAULT;
    cfg.consolidate_max_surfaces = 100;
    int ret = alea_merge_void_cells(sys, vr, &cfg);
    ASSERT(ret > 0);
    ASSERT_EQ(alea_void_count(vr), 1);

    /* Add the single consolidated void cell */
    alea_void_add_cells(sys, vr);
    alea_build_universe_index(sys);

    /* Point inside sphere = material */
    int cell = alea_find_cell(sys, 0, 0, 0);
    ASSERT(cell >= 0);
    int mat = alea_material_at(sys, 0, 0, 0);
    ASSERT(mat > 0);

    /* Point outside sphere but inside bounds = void */
    int void_cell = alea_find_cell(sys, 4.5, 4.5, 4.5);
    ASSERT(void_cell >= 0);
    int void_mat = alea_material_at(sys, 4.5, 4.5, 4.5);
    ASSERT_EQ(void_mat, 0);

    alea_void_free(vr);
    alea_destroy(sys);
}

TEST(void_consolidate_replaces_single_merged_union) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);

    alea_config_t acfg = alea_get_config(sys);
    acfg.void_max_depth = 3;
    acfg.void_min_size = 1.0;
    alea_set_config(sys, &acfg);

    int sphere_si = alea_sphere_surface(sys, 0, 0, 0, 0, 3.0);
    int outer_si = alea_box_surface(sys, 0, -10, 10, -10, 10, -10, 10);
    ASSERT(sphere_si >= 0);
    ASSERT(outer_si >= 0);

    alea_node_id_t sphere = alea_halfspace(sys, sphere_si, -1);
    alea_node_id_t outside_outer = alea_halfspace(sys, outer_si, +1);
    ASSERT(sphere != ALEA_NODE_ID_INVALID);
    ASSERT(outside_outer != ALEA_NODE_ID_INVALID);

    int mat = alea_add_material(sys, 1);
    ASSERT(mat >= 0);
    ASSERT(alea_add_cell(sys, 1, sphere, mat, 10.0, 0) >= 0);
    ASSERT(alea_add_cell(sys, 2, outside_outer, ALEA_MATERIAL_VOID, 0.0, 0) >= 0);
    ASSERT_EQ(alea_build_universe_index(sys), 0);

    alea_bbox_t bounds = {-10, 10, -10, 10, -10, 10};
    void_result_t* vr = alea_void_generate(sys, &bounds);
    ASSERT_NOT_NULL(vr);
    ASSERT(alea_void_count(vr) > 1);

    alea_void_merge_config_t cfg = ALEA_VOID_MERGE_DEFAULT;
    cfg.max_surfaces_per_cell = 100000;
    cfg.consolidate_max_surfaces = 0;
    ASSERT(alea_merge_void_cells(sys, vr, &cfg) > 0);
    ASSERT_EQ(alea_void_count(vr), 1);

    cfg.consolidate_max_surfaces = 100;
    ASSERT(alea_merge_void_cells(sys, vr, &cfg) > 0);
    ASSERT_EQ(alea_void_count(vr), 1);
    ASSERT(alea_void_add_cells(sys, vr) == 1);

    FILE* f = tmpfile();
    ASSERT_NOT_NULL(f);
    ASSERT_EQ(mcnp_export_system_stream(sys, f), 0);
    ASSERT_EQ(fseek(f, 0, SEEK_END), 0);
    long len = ftell(f);
    ASSERT(len > 0);
    ASSERT_EQ(fseek(f, 0, SEEK_SET), 0);

    char* text = calloc((size_t)len + 1, 1);
    ASSERT_NOT_NULL(text);
    ASSERT_EQ(fread(text, 1, (size_t)len, f), (size_t)len);
    fclose(f);

    ASSERT_NULL(strstr(text, "):("));
    ASSERT(len < 4000);

    free(text);
    alea_void_free(vr);
    alea_destroy(sys);
}

TEST(void_consolidate_disabled) {
    alea_system_t* sys = make_sphere_system();
    ASSERT_NOT_NULL(sys);

    alea_bbox_t bounds = {-5, 5, -5, 5, -5, 5};
    void_result_t* vr = alea_void_generate(sys, &bounds);
    ASSERT_NOT_NULL(vr);

    size_t before = alea_void_count(vr);
    ASSERT(before > 1);

    /* Merge with consolidation disabled */
    alea_void_merge_config_t cfg = ALEA_VOID_MERGE_DEFAULT;
    cfg.consolidate_max_surfaces = 0;
    int ret = alea_merge_void_cells(sys, vr, &cfg);
    ASSERT(ret > 0);

    /* Without consolidation, count should still be > 1
     * (merge reduces but doesn't collapse to 1 for sphere) */
    ASSERT(alea_void_count(vr) > 1);

    alea_void_free(vr);
    alea_destroy(sys);
}

/* ------------------------------------------------------------------------- */
/* Benchmark                                                                  */
/* ------------------------------------------------------------------------- */

TEST(void_graveyard) {
    alea_system_t* sys = make_box_system();
    ASSERT_NOT_NULL(sys);

    alea_bbox_t bounds = {-5, 5, -5, 5, -5, 5};
    void_result_t* vr = alea_void_generate(sys, &bounds);
    ASSERT_NOT_NULL(vr);

    alea_void_add_cells(sys, vr);

    size_t cells_before = alea_cell_count(sys);
    int ret = alea_void_add_graveyard(sys, vr);
    ASSERT_EQ(ret, 1);
    /* Adds 2 cells: shell (bbox→sphere) + graveyard (outside sphere) */
    ASSERT_EQ(alea_cell_count(sys), cells_before + 2);

    alea_build_universe_index(sys);

    /* Point far outside sphere should be in graveyard cell (IMP:N=0) */
    int grav_cell = alea_find_cell(sys, 100, 100, 100);
    ASSERT(grav_cell >= 0);

    int material_id;
    ASSERT_EQ(alea_cell_get(sys, (size_t)grav_cell, NULL, &material_id,
                             NULL, NULL, NULL, NULL), 0);
    ASSERT_EQ(material_id, 0);

    /* Graveyard cell is identified by having vacuum boundary surfaces.
     * (MCNP importances like IMP:N=0 now live in mcnp_model_t, not cell entries) */
    alea_build_cell_surface_index(sys);
    alea_cell_entry_t* cell = &sys->cells.data[grav_cell];
    int has_vacuum_surf = 0;
    for (size_t si = 0; si < cell->surface_index_count; si++) {
        uint32_t surf_idx = cell->surface_indices[si];
        if (surf_idx < alea_vec_count(&sys->surfaces) &&
            sys->surfaces.data[surf_idx].boundary_type == ALEA_BOUNDARY_VACUUM) {
            has_vacuum_surf = 1;
            break;
        }
    }
    ASSERT(has_vacuum_surf);

    /* Point outside bbox but inside sphere should be in shell cell
     * Bounds [-5,5]^3, sphere R ≈ 8.83. Point (7,0,0): outside box, inside sphere. */
    int shell_cell = alea_find_cell(sys, 7, 0, 0);
    ASSERT(shell_cell >= 0);
    ASSERT_NE(shell_cell, grav_cell);

    ASSERT_EQ(alea_cell_get(sys, (size_t)shell_cell, NULL, &material_id,
                             NULL, NULL, NULL, NULL), 0);
    ASSERT_EQ(material_id, 0);

    /* Point inside material cell is still found correctly */
    int mat_cell = alea_find_cell(sys, 0, 0, 0);
    ASSERT(mat_cell >= 0);
    ASSERT_EQ(alea_material_at(sys, 0, 0, 0), 1);

    alea_void_free(vr);
    alea_destroy(sys);
}

/* ------------------------------------------------------------------------- */
/* Benchmark                                                                  */
/* ------------------------------------------------------------------------- */

TEST(void_merge_benchmark) {
    alea_bbox_t bounds = {-5, 5, -5, 5, -5, 5};
    alea_void_merge_config_t merge_cfg = ALEA_VOID_MERGE_DEFAULT;
    merge_cfg.consolidate_max_surfaces = 0;  /* Benchmark merge algorithms only */

    /* --- Face-sorted merge --- */
    alea_system_t* sys1 = make_sphere_system();
    ASSERT_NOT_NULL(sys1);
    void_result_t* vr1 = alea_void_generate(sys1, &bounds);
    ASSERT_NOT_NULL(vr1);

    size_t before1 = alea_void_count(vr1);
    ASSERT(before1 > 0);

    merge_cfg.use_greedy = false;
    double t0 = now_sec();
    int r1 = alea_merge_void_cells(sys1, vr1, &merge_cfg);
    double t_face = now_sec() - t0;
    ASSERT(r1 > 0);
    size_t after_face = alea_void_count(vr1);

    /* --- Greedy merge --- */
    alea_system_t* sys2 = make_sphere_system();
    ASSERT_NOT_NULL(sys2);
    void_result_t* vr2 = alea_void_generate(sys2, &bounds);
    ASSERT_NOT_NULL(vr2);

    size_t before2 = alea_void_count(vr2);
    ASSERT_EQ(before1, before2);

    merge_cfg.use_greedy = true;
    t0 = now_sec();
    int r2 = alea_merge_void_cells(sys2, vr2, &merge_cfg);
    double t_greedy = now_sec() - t0;
    ASSERT(r2 > 0);
    size_t after_greedy = alea_void_count(vr2);

    printf("  Void merge benchmark (sphere R=3, bounds [-5,5]^3, depth=4):\n");
    printf("    Regions before merge: %zu\n", before1);
    printf("    Face-sorted: %zu regions in %.3f ms\n", after_face, t_face * 1000);
    printf("    Greedy:      %zu regions in %.3f ms\n", after_greedy, t_greedy * 1000);
    printf("    Speedup:     %.1fx\n", t_greedy / (t_face > 1e-9 ? t_face : 1e-9));

    alea_void_free(vr1);
    alea_destroy(sys1);
    alea_void_free(vr2);
    alea_destroy(sys2);
}

TEST_MAIN()
