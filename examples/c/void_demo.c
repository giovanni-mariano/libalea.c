// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/*
 * void_demo.c - Void generation demonstration
 *
 * Demonstrates:
 *   - Building geometry with multiple cells
 *   - Generating void regions with explicit bounds
 *   - Inspecting void bounding boxes
 *   - Merging void regions
 *   - Adding void cells and verifying coverage
 *   - MCNP export before and after merge
 *
 * Build:
 *   make -C examples/c void_demo
 *
 * Run:
 *   examples/c/void_demo
 */

#include <stdio.h>
#include <alea.h>
#include <alea_mcnp.h>

static void print_separator(void) {
    printf("────────────────────────────────────────────────\n");
}

/* Build the demo geometry: box [-2,2]^3 + sphere R=5 at (8,0,0) */
static alea_system_t* build_geometry(void) {
    alea_system_t* sys = alea_create();
    if (!sys) return NULL;

    /* Cell 1: box [-2,2]^3 at origin, material 1 */
    int si_box = alea_box_surface(sys, 0, -2, 2, -2, 2, -2, 2);
    alea_node_id_t box = alea_halfspace(sys, si_box, -1);
    int m1 = alea_add_material(sys, 1);
    alea_add_cell(sys, 1, box, m1, 7.8, 0);

    /* Cell 2: sphere R=5 at (8,0,0), material 2 */
    int si_sph = alea_sphere_surface(sys, 0, 8, 0, 0, 5.0);
    alea_node_id_t sph = alea_halfspace(sys, si_sph, -1);
    int m2 = alea_add_material(sys, 2);
    alea_add_cell(sys, 2, sph, m2, 2.7, 0);

    /* Tune octree for fast demo */
    alea_config_t cfg = alea_get_config(sys);
    cfg.void_max_depth = 4;
    cfg.void_min_size = 2.0;
    alea_set_config(sys, &cfg);

    alea_build_universe_index(sys);
    return sys;
}

int main(void) {
    printf("Alea %s - Void Generation Demo\n\n", alea_version());

    alea_bbox_t bounds = {-6, 16, -8, 8, -8, 8};

    /* ══════════════════════════════════════════════════
     * Part 1: Generate void, export MCNP before merge
     * ══════════════════════════════════════════════════ */

    alea_system_t* sys = build_geometry();
    if (!sys) {
        fprintf(stderr, "Failed to create system\n");
        return 1;
    }

    printf("Geometry: box [-2,2]^3 + sphere R=5 at (8,0,0)\n");
    printf("Cells: %zu, Surfaces: %zu\n\n",
           alea_cell_count(sys), alea_surface_count(sys));

    /* ── Generate void ── */

    print_separator();
    printf("Void generation\n");
    print_separator();

    printf("Search bounds: [%.0f,%.0f] x [%.0f,%.0f] x [%.0f,%.0f]\n",
           bounds.min_x, bounds.max_x, bounds.min_y, bounds.max_y,
           bounds.min_z, bounds.max_z);

    double bounds_vol = (bounds.max_x - bounds.min_x) *
                        (bounds.max_y - bounds.min_y) *
                        (bounds.max_z - bounds.min_z);
    printf("Search volume: %.0f\n\n", bounds_vol);

    void_result_t* vr = alea_void_generate(sys, &bounds);
    if (!vr) {
        fprintf(stderr, "Void generation failed: %s\n", alea_error());
        alea_destroy(sys);
        return 1;
    }

    size_t count = alea_void_count(vr);
    printf("Void regions found: %zu\n\n", count);

    /* Print each region */
    double total_void_vol = 0;
    size_t display_limit = count < 10 ? count : 10;

    for (size_t i = 0; i < display_limit; i++) {
        alea_bbox_t b;
        alea_void_get(vr, i, &b);

        double dx = b.max_x - b.min_x;
        double dy = b.max_y - b.min_y;
        double dz = b.max_z - b.min_z;
        double vol = dx * dy * dz;
        total_void_vol += vol;

        printf("  [%2zu] [%6.2f,%6.2f] x [%6.2f,%6.2f] x [%6.2f,%6.2f]"
               "  dim=%.2f x %.2f x %.2f  vol=%.1f\n",
               i, b.min_x, b.max_x, b.min_y, b.max_y,
               b.min_z, b.max_z, dx, dy, dz, vol);
    }

    /* Sum remaining volumes */
    for (size_t i = display_limit; i < count; i++) {
        alea_bbox_t b;
        alea_void_get(vr, i, &b);
        total_void_vol += (b.max_x - b.min_x) *
                          (b.max_y - b.min_y) *
                          (b.max_z - b.min_z);
    }
    if (count > display_limit)
        printf("  ... and %zu more regions\n", count - display_limit);

    printf("\nTotal void bbox volume: %.1f (%.1f%% of search bounds)\n",
           total_void_vol, 100.0 * total_void_vol / bounds_vol);

    /* ── Export MCNP before merge ── */

    printf("\n");
    print_separator();
    printf("MCNP export (before merge)\n");
    print_separator();

    alea_void_add_cells(sys, vr);
    alea_build_universe_index(sys);

    {
        FILE* fout = fopen("void_before_merge.inp", "w");
        if (fout && mcnp_export_system_stream(sys, fout) == 0) {
            printf("Written: void_before_merge.inp (%zu cells)\n", alea_cell_count(sys));
        } else {
            fprintf(stderr, "MCNP export failed: %s\n", alea_error());
        }
        if (fout) fclose(fout);
    }

    alea_void_free(vr);
    alea_destroy(sys);

    /* ══════════════════════════════════════════════════
     * Part 2: Fresh system, merge, export MCNP after
     * ══════════════════════════════════════════════════ */

    sys = build_geometry();
    if (!sys) {
        fprintf(stderr, "Failed to recreate system\n");
        return 1;
    }

    vr = alea_void_generate(sys, &bounds);
    if (!vr) {
        fprintf(stderr, "Void generation failed: %s\n", alea_error());
        alea_destroy(sys);
        return 1;
    }

    /* ── Merge ── */

    printf("\n");
    print_separator();
    printf("Void merge\n");
    print_separator();

    size_t before_merge = alea_void_count(vr);
    int merged = alea_void_merge(sys, vr);
    if (merged >= 0) {
        printf("Before merge: %zu regions\n", before_merge);
        printf("After merge:  %zu regions\n", alea_void_count(vr));
        printf("Reduction:    %.1f%%\n",
               100.0 * (1.0 - (double)alea_void_count(vr) / before_merge));
    }

    /* ── Add void cells ── */

    printf("\n");
    print_separator();
    printf("Adding void cells\n");
    print_separator();

    size_t cells_before = alea_cell_count(sys);
    int added = alea_void_add_cells(sys, vr);
    printf("Cells before: %zu\n", cells_before);
    printf("Void cells added: %d\n", added);

    /* Add shell (bbox→sphere, IMP:N=1) + graveyard (outside sphere, IMP:N=0) */
    int grav = alea_void_add_graveyard(sys, vr);
    if (grav > 0)
        printf("Shell + graveyard cells added\n");

    printf("Total cells: %zu\n", alea_cell_count(sys));

    alea_build_universe_index(sys);

    /* ── Simplify cells ── */

    printf("\n");
    print_separator();
    printf("Cell simplification\n");
    print_separator();

    size_t cells_before_simplify = alea_cell_count(sys);
    alea_simplify_stats_t simp_stats = {0};
    alea_flatten_all_cells(sys, &simp_stats);
    alea_build_universe_index(sys);

    printf("Cells before: %zu\n", cells_before_simplify);
    printf("Cells after:  %zu\n", alea_cell_count(sys));
    printf("Nodes:        %zu -> %zu\n", simp_stats.nodes_before, simp_stats.nodes_after);
    if (simp_stats.empty_cells_removed > 0)
        printf("Empty cells removed: %zu\n", simp_stats.empty_cells_removed);
    if (simp_stats.contradictions_found > 0)
        printf("Contradictions found: %zu\n", simp_stats.contradictions_found);

    /* ── Export MCNP after merge + simplify ── */

    printf("\n");
    print_separator();
    printf("MCNP export (after merge + simplify)\n");
    print_separator();

    {
        FILE* fout = fopen("void_after_merge.inp", "w");
        if (fout && mcnp_export_system_stream(sys, fout) == 0) {
            printf("Written: void_after_merge.inp (%zu cells)\n", alea_cell_count(sys));
        } else {
            fprintf(stderr, "MCNP export failed: %s\n", alea_error());
        }
        if (fout) fclose(fout);
    }

    /* ── Point verification ── */

    printf("\n");
    print_separator();
    printf("Point verification\n");
    print_separator();

    struct {
        double x, y, z;
        const char* desc;
    } points[] = {
        {  0.0,  0.0,  0.0, "center of box (material)" },
        {  8.0,  0.0,  0.0, "center of sphere (material)" },
        {  4.0,  0.0,  0.0, "inside sphere (mat)" },
        { -5.0,  0.0,  0.0, "outside left (void)" },
        { 15.0,  0.0,  0.0, "outside right (void)" },
        {  0.0,  7.0,  7.0, "corner region (void)" },
    };

    printf("%-35s  %8s  %8s\n", "Point", "Cell", "Material");
    printf("%-35s  %8s  %8s\n", "─────", "────", "────────");

    for (size_t i = 0; i < sizeof(points)/sizeof(points[0]); i++) {
        double x = points[i].x, y = points[i].y, z = points[i].z;
        int cell = -1;
        int mat = 0;
        int cell_id = -1;
        if (alea_find_cell_at(sys, x, y, z, &cell_id, &mat) == 0) {
            cell = alea_cell_find(sys, cell_id);
        }

        const char* kind = (cell < 0) ? "MISS" :
                           (mat == 0) ? "void" : "mat";

        printf("(%5.1f,%5.1f,%5.1f) %-18s  %6d  %6d  [%s]\n",
               x, y, z, points[i].desc, cell, mat, kind);
    }

    /* ── Cleanup ── */

    alea_void_free(vr);
    alea_destroy(sys);

    printf("\nDone.\n");
    return 0;
}
