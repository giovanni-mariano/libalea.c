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
 *   --origin X Y Z  Sphere center (default: auto from bounding sphere)
 *   --radius R      Sphere radius (default: auto from bounding sphere)
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
    fprintf(stderr, "  --origin X Y Z  Sphere center (default: auto from bounding sphere)\n");
    fprintf(stderr, "  --radius R      Sphere radius (default: auto from bounding sphere)\n");
    fprintf(stderr, "\nWhen --origin/--radius are not given, they are auto-computed\n");
    fprintf(stderr, "via alea_compute_bounding_sphere().\n");
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
    double user_cx = 0, user_cy = 0, user_cz = 0, user_radius = 0;
    int have_origin = 0, have_radius = 0;

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
        } else if (strcmp(argv[i], "--origin") == 0) {
            if (i + 3 >= argc) {
                fprintf(stderr, "Error: --origin requires 3 arguments (X Y Z)\n");
                return 1;
            }
            user_cx = atof(argv[i + 1]);
            user_cy = atof(argv[i + 2]);
            user_cz = atof(argv[i + 3]);
            have_origin = 1;
            i += 4;
        } else if (strcmp(argv[i], "--radius") == 0) {
            if (i + 1 >= argc) {
                fprintf(stderr, "Error: --radius requires 1 argument\n");
                return 1;
            }
            user_radius = atof(argv[i + 1]);
            if (user_radius <= 0) {
                fprintf(stderr, "Error: radius must be positive\n");
                return 1;
            }
            have_radius = 1;
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

    size_t nc = alea_cell_count(sys);
    printf("Cells: %zu\n", nc);

    alea_build_universe_index(sys);

    /* Determine sphere parameters */
    double cx, cy, cz, radius;
    if (have_origin && have_radius) {
        cx = user_cx;
        cy = user_cy;
        cz = user_cz;
        radius = user_radius;
        printf("User sphere: center=(%.2f,%.2f,%.2f) R=%.2f\n",
               cx, cy, cz, radius);
    } else {
        printf("Computing bounding sphere (tol=%.2f)...\n", tol);
        if (alea_compute_bounding_sphere(sys, tol, &cx, &cy, &cz, &radius) != 0) {
            fprintf(stderr, "Error: could not compute bounding sphere\n");
            mcnp_model_destroy(model);
            return 1;
        }
        printf("Bounding sphere: center=(%.2f,%.2f,%.2f) R=%.2f\n",
               cx, cy, cz, radius);
    }

    /* Compute ray count */
    int rays = (n_rays > 0) ? n_rays : compute_n_rays(radius, ray_density);
    printf("Rays: %d (density=%.1f rays/cm^2, pi*R^2=%.0f cm^2)\n\n",
           rays, ray_density, M_PI * radius * radius);

    /* Estimate volumes */
    double* volumes = calloc(nc, sizeof(double));
    double* rel_errors = calloc(nc, sizeof(double));
    if (!volumes || !rel_errors) {
        fprintf(stderr, "Error: out of memory\n");
        free(volumes);
        free(rel_errors);
        mcnp_model_destroy(model);
        return 1;
    }

    alea_estimate_cell_volumes(sys, cx, cy, cz, radius, rays,
                                   volumes, rel_errors);

    /* Print table */
    printf("%5s  %10s  %14s  %7s\n", "Cell", "Material", "Volume", "R");
    for (size_t j = 0; j < nc; j++) {
        int cell_id, material_id;
        alea_cell_get(sys, j, &cell_id, &material_id, NULL, NULL, NULL, NULL);
        if (rel_errors[j] >= 0.0) {
            printf("%5d  %10d  %14.6e  %7.4f\n",
                   cell_id, material_id, volumes[j], rel_errors[j]);
        } else {
            printf("%5d  %10d  %14.6e  %7s\n",
                   cell_id, material_id, volumes[j], "-");
        }
    }

    free(volumes);
    free(rel_errors);
    mcnp_model_destroy(model);
    return 0;
}
