// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/*
 * mcnp_mesh.c - Export MCNP geometry as a structured hex mesh
 *
 * Usage:
 *   ./mcnp_mesh <input.inp> -o <output> [options]
 *
 * Options:
 *   -o FILE          Output file (required; .msh or .vtk)
 *   --format FMT     gmsh (default) or vtk
 *   --nx N           Elements in X (default: 50)
 *   --ny N           Elements in Y (default: 50)
 *   --nz N           Elements in Z (default: 50)
 *   --bounds X0,X1,Y0,Y1,Z0,Z1   Manual bounding box
 *   --void-mat ID    Material ID for void (default: 0)
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <alea.h>
#include <alea_mesh.h>

static void print_usage(const char *prog) {
    fprintf(stderr, "Usage: %s <input.inp> -o <output> [options]\n\n", prog);
    fprintf(stderr, "Export MCNP geometry as a structured hex mesh.\n\n");
    fprintf(stderr, "Options:\n");
    fprintf(stderr, "  -o FILE          Output file (.msh or .vtk)\n");
    fprintf(stderr, "  --format FMT     gmsh (default) or vtk\n");
    fprintf(stderr, "  --nx N           Elements in X (default: 50)\n");
    fprintf(stderr, "  --ny N           Elements in Y (default: 50)\n");
    fprintf(stderr, "  --nz N           Elements in Z (default: 50)\n");
    fprintf(stderr, "  --bounds X0,X1,Y0,Y1,Z0,Z1   Bounding box\n");
    fprintf(stderr, "  --void-mat ID    Material ID for void (default: 0)\n");
}

int main(int argc, char **argv) {
    if (argc < 2 || strcmp(argv[1], "-h") == 0 || strcmp(argv[1], "--help") == 0) {
        print_usage(argv[0]);
        return argc < 2 ? 1 : 0;
    }

    const char *input_file = argv[1];
    const char *output_file = NULL;

    alea_mesh_config_t cfg;
    alea_mesh_config_init(&cfg);
    cfg.nx = cfg.ny = cfg.nz = 50;

    int explicit_format = 0;

    int i = 2;
    while (i < argc) {
        if (strcmp(argv[i], "-o") == 0) {
            if (i + 1 >= argc) { fprintf(stderr, "Error: -o requires argument\n"); return 1; }
            output_file = argv[++i];
        } else if (strcmp(argv[i], "--format") == 0) {
            if (i + 1 >= argc) { fprintf(stderr, "Error: --format requires argument\n"); return 1; }
            i++;
            if (strcmp(argv[i], "gmsh") == 0 || strcmp(argv[i], "msh") == 0) {
                cfg.format = ALEA_MESH_GMSH;
            } else if (strcmp(argv[i], "vtk") == 0) {
                cfg.format = ALEA_MESH_VTK;
            } else {
                fprintf(stderr, "Error: unknown format '%s' (use gmsh or vtk)\n", argv[i]);
                return 1;
            }
            explicit_format = 1;
        } else if (strcmp(argv[i], "--nx") == 0) {
            if (i + 1 >= argc) { fprintf(stderr, "Error: --nx requires argument\n"); return 1; }
            cfg.nx = atoi(argv[++i]);
        } else if (strcmp(argv[i], "--ny") == 0) {
            if (i + 1 >= argc) { fprintf(stderr, "Error: --ny requires argument\n"); return 1; }
            cfg.ny = atoi(argv[++i]);
        } else if (strcmp(argv[i], "--nz") == 0) {
            if (i + 1 >= argc) { fprintf(stderr, "Error: --nz requires argument\n"); return 1; }
            cfg.nz = atoi(argv[++i]);
        } else if (strcmp(argv[i], "--bounds") == 0) {
            if (i + 1 >= argc) { fprintf(stderr, "Error: --bounds requires argument\n"); return 1; }
            i++;
            if (sscanf(argv[i], "%lf,%lf,%lf,%lf,%lf,%lf",
                       &cfg.x_min, &cfg.x_max, &cfg.y_min, &cfg.y_max,
                       &cfg.z_min, &cfg.z_max) != 6) {
                fprintf(stderr, "Error: --bounds format: X0,X1,Y0,Y1,Z0,Z1\n");
                return 1;
            }
        } else if (strcmp(argv[i], "--void-mat") == 0) {
            if (i + 1 >= argc) { fprintf(stderr, "Error: --void-mat requires argument\n"); return 1; }
            cfg.void_material_id = atoi(argv[++i]);
        } else if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0) {
            print_usage(argv[0]);
            return 0;
        } else {
            fprintf(stderr, "Error: unexpected argument '%s'\n", argv[i]);
            print_usage(argv[0]);
            return 1;
        }
        i++;
    }

    if (!output_file) {
        fprintf(stderr, "Error: -o <output> is required\n");
        print_usage(argv[0]);
        return 1;
    }

    /* Auto-detect format from extension if not explicitly set */
    if (!explicit_format) {
        const char *ext = strrchr(output_file, '.');
        if (ext && strcmp(ext, ".vtk") == 0) {
            cfg.format = ALEA_MESH_VTK;
        }
    }

    printf("Alea %s - Mesh Export\n\n", alea_version());

    /* Load */
    printf("Loading: %s\n", input_file);
    alea_system_t *sys = alea_load_mcnp(input_file);
    if (!sys) {
        fprintf(stderr, "Error: %s\n", alea_error());
        return 1;
    }
    printf("Cells: %zu\n", alea_cell_count(sys));

    alea_build_universe_index(sys);

    /* Sample */
    printf("Grid: %d x %d x %d = %d elements\n",
           cfg.nx, cfg.ny, cfg.nz, cfg.nx * cfg.ny * cfg.nz);
    if (cfg.x_min == 0.0 && cfg.x_max == 0.0)
        printf("Bounds: auto-detect\n");
    else
        printf("Bounds: [%.2f,%.2f] x [%.2f,%.2f] x [%.2f,%.2f]\n",
               cfg.x_min, cfg.x_max, cfg.y_min, cfg.y_max, cfg.z_min, cfg.z_max);

    alea_mesh_result_t *mesh = alea_mesh_sample(sys, &cfg);
    if (!mesh) {
        fprintf(stderr, "Error: mesh sampling failed\n");
        alea_destroy(sys);
        return 1;
    }

    printf("Unique materials: %d\n", mesh->num_materials);
    printf("Bounds used: [%.4f,%.4f] x [%.4f,%.4f] x [%.4f,%.4f]\n",
           mesh->x_nodes[0], mesh->x_nodes[mesh->nx],
           mesh->y_nodes[0], mesh->y_nodes[mesh->ny],
           mesh->z_nodes[0], mesh->z_nodes[mesh->nz]);

    /* Export */
    const char *fmt_name = cfg.format == ALEA_MESH_VTK ? "VTK" : "Gmsh";
    printf("Writing %s: %s\n", fmt_name, output_file);

    int rc = alea_mesh_export(mesh, cfg.format, output_file);
    if (rc != 0) {
        fprintf(stderr, "Error: export failed\n");
        alea_mesh_result_free(mesh);
        alea_destroy(sys);
        return 1;
    }

    printf("Done.\n");

    alea_mesh_result_free(mesh);
    alea_destroy(sys);
    return 0;
}
