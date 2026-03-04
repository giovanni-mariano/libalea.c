// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/*
 * basic.c - Basic Alea usage example
 *
 * Demonstrates:
 *   - Creating geometry programmatically
 *   - Point-in-cell queries
 *   - Void region generation
 *   - MCNP export
 *
 * Build:
 *   gcc -o basic basic.c -I../include ../bin/libalea.a -lm
 *
 * Run:
 *   ./basic
 */

#include <stdio.h>
#include <alea.h>
#include <alea_mcnp.h>

int main(void) {
    printf("Alea %s - Basic Example\n\n", alea_version());

    /* Create a new CSG system */
    alea_system_t* sys = alea_create();
    if (!sys) {
        fprintf(stderr, "Failed to create system\n");
        return 1;
    }

    /*
     * Build a simple geometry: sphere with cylindrical hole
     *
     *   Outer sphere: radius 10 cm, centered at origin
     *   Inner cylinder: radius 3 cm along Z-axis
     *   Material: 1 (fuel)
     */

    /* Create outer sphere (sense=-1 means inside) */
    int si_outer = alea_sphere_surface(sys, 0, 0, 0, 0, 10.0);
    alea_node_id_t outer = alea_halfspace(sys, si_outer, -1);

    /* Create cylinder along Z (sense=-1 means inside) */
    int si_hole = alea_cylinder_z_surface(sys, 0, 0, 0, 3.0);
    alea_node_id_t hole = alea_halfspace(sys, si_hole, -1);

    /* Subtract cylinder from sphere: sphere AND NOT(cylinder) */
    alea_node_id_t region = alea_difference(sys, outer, hole);

    /* Register material 1 (fuel) */
    int m1 = alea_add_material(sys, 1);

    /* Add as cell 1 with material m1, density 10.0 g/cc, universe 0 */
    int cell_idx = alea_add_cell(sys, 1, region, m1, 10.0, 0);
    if (cell_idx < 0) {
        fprintf(stderr, "Failed to add cell\n");
        alea_destroy(sys);
        return 1;
    }

    /* Build universe index for queries */
    if (alea_build_universe_index(sys) != 0) {
        fprintf(stderr, "Failed to build universe index\n");
        alea_destroy(sys);
        return 1;
    }

    printf("Geometry created: sphere (r=10) with cylindrical hole (r=3)\n\n");

    /* Test some points */
    struct {
        double x, y, z;
        const char* desc;
    } test_points[] = {
        {  5.0,  0.0,  0.0, "inside sphere, outside hole" },
        {  0.0,  0.0,  0.0, "on axis (inside hole)"       },
        {  0.0,  2.0,  0.0, "inside hole"                 },
        { 15.0,  0.0,  0.0, "outside sphere"              },
        {  7.0,  7.0,  0.0, "corner (outside sphere)"     },
    };

    printf("Point queries:\n");
    printf("%-35s  %8s  %8s\n", "Point", "Cell", "Material");
    printf("%-35s  %8s  %8s\n", "-----", "----", "--------");

    for (size_t i = 0; i < sizeof(test_points)/sizeof(test_points[0]); i++) {
        double x = test_points[i].x;
        double y = test_points[i].y;
        double z = test_points[i].z;

        int cell = alea_find_cell(sys, x, y, z);
        int mat = alea_material_at(sys, x, y, z);

        printf("(%5.1f, %5.1f, %5.1f) %-15s  %8d  %8d\n",
               x, y, z, test_points[i].desc, cell, mat);
    }

    /*
     * Generate void regions
     *
     * Void generation finds all regions within a bounding box that are
     * not covered by any cell. This is useful for creating "graveyard"
     * or external void cells in MCNP.
     */
    printf("\nVoid generation:\n");

    /* Define bounding box: [-15, 15] in all dimensions */
    alea_bbox_t bounds = {
        .min_x = -15.0, .max_x = 15.0,
        .min_y = -15.0, .max_y = 15.0,
        .min_z = -15.0, .max_z = 15.0
    };

    void_result_t* voids = alea_void_generate(sys, &bounds);
    if (voids) {
        size_t void_count = alea_void_count(voids);
        printf("  Found %zu void region(s) within bounds [%.0f, %.0f]^3\n",
               void_count, bounds.min_x, bounds.max_x);

        /* Print first few void regions */
        for (size_t i = 0; i < void_count && i < 3; i++) {
            alea_bbox_t box;
            if (alea_void_get(voids, i, &box) == 0) {
                printf("  Void %zu: [%.1f,%.1f] x [%.1f,%.1f] x [%.1f,%.1f]\n",
                       i, box.min_x, box.max_x, box.min_y, box.max_y,
                       box.min_z, box.max_z);
            }
        }
        if (void_count > 3) {
            printf("  ... and %zu more\n", void_count - 3);
        }

        /* Add void cells to the system */
        int added = alea_void_add_cells(sys, voids);
        if (added > 0) {
            printf("  Added %d void cell(s) to geometry\n", added);
        }

        alea_void_free(voids);
    } else {
        printf("  Void generation failed\n");
    }

    /* Export to MCNP format (generic export, no MCNP model needed) */
    const char* output_file = "basic_output.inp";
    FILE* fout = fopen(output_file, "w");
    if (fout && mcnp_export_system_stream(sys, fout) == 0) {
        printf("\nExported to: %s\n", output_file);
    } else {
        fprintf(stderr, "\nFailed to export: %s\n", alea_error());
    }
    if (fout) fclose(fout);

    /* Cleanup */
    alea_destroy(sys);

    printf("\nDone.\n");
    return 0;
}
