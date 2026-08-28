// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/*
 * mcnp_volume.c - Estimate cell volumes for an MCNP model
 *
 * Usage:
 *   ./mcnp_volume <input.inp> [options]
 *
 * Options:
 *   --rays N        Number of rays (default: auto-compute from bounding sphere)
 *   --density D     Ray density in rays/cm^2 (default: 1.0)
 *   --tol T         Bbox tightening tolerance (default: 1.0)
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <alea.h>
#include <alea_mcnp.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define DEFAULT_RAY_DENSITY 1.0
#define MIN_RAYS            1000

static void print_usage(const char* prog) {
    fprintf(stderr, "Usage: %s <input.inp> [options]\n\n", prog);
    fprintf(stderr, "Estimate cell volumes via Cauchy-Crofton ray tracing.\n\n");
    fprintf(stderr, "Options:\n");
    fprintf(stderr, "  --rays N        Number of rays (default: auto from bounding sphere)\n");
    fprintf(stderr, "  --density D     Ray density in rays/cm^2 (default: %.1f)\n", DEFAULT_RAY_DENSITY);
    fprintf(stderr, "  --tol T         Bbox tightening tolerance (default: 1.0)\n");
    fprintf(stderr, "\nThe estimator computes its physical-path bounding sphere automatically.\n");
    fprintf(stderr, "Auto rays: n_rays = max(%d, ceil(density * pi * R^2))\n", MIN_RAYS);
}

static int compute_n_rays(double radius, double density) {
    double n = density * M_PI * radius * radius;
    int result = (int)ceil(n);
    return result > MIN_RAYS ? result : MIN_RAYS;
}

int main(int argc, char** argv) {
    if (argc < 2 || strcmp(argv[1], "-h") == 0 || strcmp(argv[1], "--help") == 0) {
        print_usage(argv[0]);
        return argc < 2 ? 1 : 0;
    }

    /* Parse arguments */
    const char* input_file = argv[1];
    int n_rays = 0;             /* 0 = auto-compute */
    double ray_density = DEFAULT_RAY_DENSITY;
    double tol = 1.0;

    int i = 2;
    while (i < argc) {
        if (strcmp(argv[i], "--rays") == 0) {
            if (i + 1 >= argc) {
                fprintf(stderr, "Error: --rays requires 1 argument\n");
                return 1;
            }
            n_rays = atoi(argv[i + 1]);
            if (n_rays <= 0) {
                fprintf(stderr, "Error: ray count must be positive\n");
                return 1;
            }
            i += 2;
        } else if (strcmp(argv[i], "--density") == 0) {
            if (i + 1 >= argc) {
                fprintf(stderr, "Error: --density requires 1 argument\n");
                return 1;
            }
            ray_density = atof(argv[i + 1]);
            if (ray_density <= 0) {
                fprintf(stderr, "Error: ray density must be positive\n");
                return 1;
            }
            i += 2;
        } else if (strcmp(argv[i], "--tol") == 0) {
            if (i + 1 >= argc) {
                fprintf(stderr, "Error: --tol requires 1 argument\n");
                return 1;
            }
            tol = atof(argv[i + 1]);
            if (tol <= 0) {
                fprintf(stderr, "Error: tol must be positive\n");
                return 1;
            }
            i += 2;
        } else if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0) {
            print_usage(argv[0]);
            return 0;
        } else {
            fprintf(stderr, "Error: unexpected argument '%s'\n", argv[i]);
            print_usage(argv[0]);
            return 1;
        }
    }

    printf("Alea %s - Volume Estimation\n\n", alea_version());

    /* Load MCNP input */
    printf("Loading: %s\n", input_file);
    mcnp_model_t* model = mcnp_load(input_file);
    if (!model) {
        fprintf(stderr, "Error: %s\n", alea_error());
        return 1;
    }
    alea_system_t* sys = model->sys;

    size_t n_paths = alea_volume_path_count(sys);
    printf("Physical volume paths: %zu\n", n_paths);

    alea_build_universe_index(sys);

    /* Determine a representative sphere for automatic ray-count selection. */
    double cx, cy, cz, radius;
    printf("Computing bounding sphere (tol=%.2f)...\n", tol);
    if (alea_compute_bounding_sphere(sys, tol, &cx, &cy, &cz, &radius) != 0) {
        fprintf(stderr, "Error: could not compute bounding sphere\n");
        mcnp_model_destroy(model);
        return 1;
    }
    printf("Bounding sphere: center=(%.2f,%.2f,%.2f) R=%.2f\n",
           cx, cy, cz, radius);

    /* Compute ray count */
    int rays = (n_rays > 0) ? n_rays : compute_n_rays(radius, ray_density);
    printf("Rays: %d (density=%.1f rays/cm^2, pi*R^2=%.0f cm^2)\n\n",
           rays, ray_density, M_PI * radius * radius);

    /* Estimate volumes */
    double* volumes = calloc(n_paths, sizeof(double));
    double* rel_errors = calloc(n_paths, sizeof(double));
    alea_volume_path_t* paths = calloc(n_paths, sizeof(*paths));
    if (!volumes || !rel_errors || !paths) {
        fprintf(stderr, "Error: out of memory\n");
        free(volumes);
        free(rel_errors);
        free(paths);
        mcnp_model_destroy(model);
        return 1;
    }

    if (alea_volume_paths_get(sys, paths, n_paths) != n_paths ||
        alea_estimate_volumes(sys, rays, volumes, rel_errors) != 0) {
        fprintf(stderr, "Error: %s\n", alea_error());
        free(volumes);
        free(rel_errors);
        free(paths);
        mcnp_model_destroy(model);
        return 1;
    }

    /* Print table */
    printf("%5s  %7s  %10s  %14s  %7s\n",
           "Path", "Cell", "Material", "Volume", "R");
    for (size_t j = 0; j < n_paths; j++) {
        if (rel_errors[j] >= 0.0) {
            printf("%5llu  %7d  %10d  %14.6e  %7.4f\n",
                   (unsigned long long)paths[j].path_id,
                   paths[j].terminal_cell_id, paths[j].material_id,
                   volumes[j], rel_errors[j]);
        } else {
            printf("%5llu  %7d  %10d  %14.6e  %7s\n",
                   (unsigned long long)paths[j].path_id,
                   paths[j].terminal_cell_id, paths[j].material_id,
                   volumes[j], "-");
        }
    }

    free(volumes);
    free(rel_errors);
    free(paths);
    mcnp_model_destroy(model);
    return 0;
}
