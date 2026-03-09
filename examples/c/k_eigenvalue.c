// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/**
 * @file k_eigenvalue.c
 * @brief Simplified k-eigenvalue calculation for a bare sphere
 *
 * Estimates k_eff for a homogeneous sphere of fissile material using
 * the power iteration method (analog Monte Carlo).
 *
 * Usage:
 *   ./k_eigenvalue <path-to-fendl-ace-dir> [radius_cm]
 */

#include "alea_nucdata.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static uint64_t rng_state;
static double rng(void)
{
    rng_state = rng_state * 6364136223846793005ULL + 1442695040888963407ULL;
    return (double)(rng_state >> 11) / (double)(1ULL << 53);
}

/* Particle bank for fission source */
typedef struct { double x, y, z, energy; } source_point_t;

int main(int argc, char* argv[])
{
    const char* datadir = "test/fendl-FENDL-3.2c-neutron-ace/neutron/ace";
    double radius = 8.0; /* cm — approximate bare U-235 critical radius */

    if (argc > 1) datadir = argv[1];
    if (argc > 2) radius = atof(argv[2]);

    rng_state = (uint64_t)time(NULL) ^ 0xDEADBEEF;

    /* Load nuclear data */
    alea_nuc_xsdir_t* xsdir = alea_nuc_xsdir_load_dir(datadir);
    if (!xsdir) {
        fprintf(stderr, "Failed to load data from %s\n", datadir);
        return 1;
    }

    alea_nuc_nuclide_t* u235 = alea_nuc_load_nuclide(xsdir, "92235.32c");
    if (!u235) { fprintf(stderr, "U-235 not found\n"); alea_nuc_xsdir_free(xsdir); return 1; }

    /* U-235 metal: N = 0.04833 atoms/barn-cm (density 18.7 g/cm³) */
    alea_nuc_material_t* fuel = alea_nuc_material_create();
    alea_nuc_material_add(fuel, u235, 0.04833);

    int n_particles = 5000;
    int n_batches = 30;
    int n_inactive = 10;

    printf("=== k-eigenvalue: Bare U-235 Sphere ===\n");
    printf("Radius: %.1f cm\n", radius);
    printf("Batches: %d (%d inactive + %d active)\n",
           n_batches, n_inactive, n_batches - n_inactive);
    printf("Particles/batch: %d\n\n", n_particles);

    /* Initialize source: uniform in sphere, 1 MeV */
    source_point_t* source = malloc((size_t)n_particles * 4 * sizeof(source_point_t));
    source_point_t* fission_bank = source + n_particles;

    for (int i = 0; i < n_particles; i++) {
        /* Uniform sampling in sphere */
        double r = radius * cbrt(rng());
        double costh = 2.0 * rng() - 1.0;
        double phi = 2.0 * M_PI * rng();
        double sinth = sqrt(1.0 - costh * costh);
        source[i].x = r * sinth * cos(phi);
        source[i].y = r * sinth * sin(phi);
        source[i].z = r * costh;
        source[i].energy = 1.0;
    }

    double k_sum = 0.0, k_sq_sum = 0.0;
    int n_active = 0;

    for (int batch = 0; batch < n_batches; batch++) {
        int fission_count = 0;

        for (int p = 0; p < n_particles; p++) {
            double x = source[p].x, y = source[p].y, z = source[p].z;
            double energy = source[p].energy;

            /* Isotropic initial direction */
            double mu = 2.0 * rng() - 1.0;
            double phi = 2.0 * M_PI * rng();
            double sin_th = sqrt(1.0 - mu * mu);
            double ux = sin_th * cos(phi);
            double uy = sin_th * sin(phi);
            double uz = mu;

            int alive = 1;
            int max_collisions = 5000;

            while (alive && max_collisions-- > 0) {
                double sigma_t = alea_nuc_mat_xs_total(fuel, energy);
                if (sigma_t <= 0.0) break;
                double dist = -log(1.0 - rng()) / sigma_t;

                x += dist * ux;
                y += dist * uy;
                z += dist * uz;

                /* Check if escaped sphere */
                if (x * x + y * y + z * z > radius * radius)
                    break;

                /* Sample reaction */
                int mt;
                alea_nuc_sample_reaction(u235, energy, rng(), &mt);

                /* Fission: bank new source sites */
                if (mt == 18 || mt == 19) {
                    double nu = alea_nuc_nu_bar(u235, energy);
                    int n_new = (int)(nu + rng());
                    for (int s = 0; s < n_new && fission_count < n_particles * 3; s++) {
                        fission_bank[fission_count].x = x;
                        fission_bank[fission_count].y = y;
                        fission_bank[fission_count].z = z;
                        /* Sample fission energy from fission spectrum */
                        double xi3[3] = {rng(), rng(), rng()};
                        alea_nuc_interaction_t result;
                        alea_nuc_sample_collision(u235, mt, energy, xi3, &result);
                        fission_bank[fission_count].energy = result.energy_out;
                        fission_count++;
                    }
                    alive = 0; /* analog: neutron dies at fission */
                    continue;
                }

                /* Absorption (non-fission) */
                if (mt != 2 && (mt < 51 || mt > 91)) {
                    alive = 0;
                    continue;
                }

                /* Scattering */
                double xi[3] = {rng(), rng(), rng()};
                alea_nuc_interaction_t result;
                alea_nuc_sample_collision(u235, mt, energy, xi, &result);
                energy = result.energy_out;
                if (energy <= 1e-11) break;

                /* Rotate direction */
                double mu_s = result.mu;
                double cos_phi = cos(2.0 * M_PI * rng());
                double sin_mu = sqrt(1.0 - mu_s * mu_s);
                double sin_uz = sqrt(1.0 - uz * uz);

                if (sin_uz > 1e-10) {
                    double ux_new = mu_s * ux + (sin_mu / sin_uz) *
                                    (ux * uz * cos_phi - uy * sqrt(1.0 - cos_phi * cos_phi));
                    double uy_new = mu_s * uy + (sin_mu / sin_uz) *
                                    (uy * uz * cos_phi + ux * sqrt(1.0 - cos_phi * cos_phi));
                    double uz_new = mu_s * uz - sin_mu * sin_uz * cos_phi;
                    ux = ux_new; uy = uy_new; uz = uz_new;
                } else {
                    mu = 2.0 * rng() - 1.0;
                    phi = 2.0 * M_PI * rng();
                    sin_th = sqrt(1.0 - mu * mu);
                    ux = sin_th * cos(phi);
                    uy = sin_th * sin(phi);
                    uz = mu;
                }
                /* Normalize */
                double norm = sqrt(ux*ux + uy*uy + uz*uz);
                if (norm > 0) { ux /= norm; uy /= norm; uz /= norm; }
            }
        }

        double k_batch = (double)fission_count / n_particles;

        if (batch >= n_inactive) {
            k_sum += k_batch;
            k_sq_sum += k_batch * k_batch;
            n_active++;
        }

        printf("  Batch %2d: k = %.4f  (fissions: %d)", batch + 1, k_batch, fission_count);
        if (batch < n_inactive)
            printf("  [inactive]\n");
        else
            printf("\n");

        /* Resample source from fission bank */
        if (fission_count > 0) {
            for (int i = 0; i < n_particles; i++) {
                int idx = (int)(rng() * fission_count);
                if (idx >= fission_count) idx = fission_count - 1;
                source[i] = fission_bank[idx];
            }
        }
    }

    double k_avg = k_sum / n_active;
    double k_var = k_sq_sum / n_active - k_avg * k_avg;
    double k_std = (k_var > 0) ? sqrt(k_var / n_active) : 0.0;

    printf("\n=== Results ===\n");
    printf("k_eff = %.4f ± %.4f\n", k_avg, k_std);
    printf("(Godiva critical radius ≈ 8.7 cm → k_eff ≈ 1.0)\n");

    free(source);
    alea_nuc_material_destroy(fuel);
    alea_nuc_nuclide_free(u235);
    alea_nuc_xsdir_free(xsdir);
    return 0;
}
