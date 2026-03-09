// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/**
 * @file slab_transport.c
 * @brief Simple 1D slab neutron transport using libalea nucdata
 *
 * Simulates monoenergetic neutrons hitting a slab of material and
 * tallies transmission, reflection, and absorption probabilities.
 *
 * Usage:
 *   ./slab_transport <path-to-fendl-ace-dir> [thickness_cm] [n_histories]
 *
 * Example:
 *   ./slab_transport test/fendl-FENDL-3.2c-neutron-ace/neutron/ace 5.0 100000
 */

#include "alea_nucdata.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

/* LCG random number generator */
static uint64_t rng_state;

static double rng(void)
{
    rng_state = rng_state * 6364136223846793005ULL + 1442695040888963407ULL;
    return (double)(rng_state >> 11) / (double)(1ULL << 53);
}

int main(int argc, char* argv[])
{
    const char* datadir = "test/fendl-FENDL-3.2c-neutron-ace/neutron/ace";
    double thickness = 5.0;   /* cm */
    int n_histories = 100000;

    if (argc > 1) datadir = argv[1];
    if (argc > 2) thickness = atof(argv[2]);
    if (argc > 3) n_histories = atoi(argv[3]);

    rng_state = (uint64_t)time(NULL);

    /* Load nuclear data */
    alea_nuc_xsdir_t* xsdir = alea_nuc_xsdir_load_dir(datadir);
    if (!xsdir) {
        fprintf(stderr, "Failed to load xsdir from %s\n", datadir);
        return 1;
    }

    /* Load nuclides for water (H2O) */
    alea_nuc_nuclide_t* h1 = alea_nuc_load_nuclide(xsdir, "1001.32c");
    alea_nuc_nuclide_t* o16 = alea_nuc_load_nuclide(xsdir, "8016.32c");
    if (!h1 || !o16) {
        fprintf(stderr, "Failed to load H-1 or O-16\n");
        alea_nuc_nuclide_free(h1);
        alea_nuc_nuclide_free(o16);
        alea_nuc_xsdir_free(xsdir);
        return 1;
    }

    /* Create water material (density ~1 g/cm³) */
    alea_nuc_material_t* water = alea_nuc_material_create();
    alea_nuc_material_add(water, h1, 6.676e-2);  /* atoms/barn-cm */
    alea_nuc_material_add(water, o16, 3.338e-2);

    printf("=== 1D Slab Neutron Transport ===\n");
    printf("Material: Water (H2O)\n");
    printf("Slab thickness: %.1f cm\n", thickness);
    printf("Source: 1 MeV neutrons, normal incidence\n");
    printf("Histories: %d\n\n", n_histories);

    /* Tallies */
    int transmitted = 0, reflected = 0, absorbed = 0;
    double avg_collisions = 0.0;
    int energy_bins = 10;
    int* leak_spectrum = calloc((size_t)energy_bins, sizeof(int));

    for (int hist = 0; hist < n_histories; hist++) {
        double x = 0.0;       /* position along slab */
        double mu_dir = 1.0;  /* direction cosine (initially forward) */
        double energy = 1.0;  /* MeV */
        double weight = 1.0;
        int collisions = 0;
        int alive = 1;

        while (alive) {
            /* Sample distance to next collision */
            double sigma_t = alea_nuc_mat_xs_total(water, energy);
            if (sigma_t <= 0.0) { absorbed++; break; }
            double dist = -log(1.0 - rng()) / sigma_t;

            /* Move particle */
            x += dist * mu_dir;

            /* Check boundaries */
            if (x >= thickness) {
                transmitted++;
                /* Log energy spectrum of transmitted neutrons */
                int ebin = (int)(log10(energy + 1e-11) + 5.0);
                if (ebin < 0) ebin = 0;
                if (ebin >= energy_bins) ebin = energy_bins - 1;
                leak_spectrum[ebin]++;
                break;
            }
            if (x <= 0.0) {
                reflected++;
                break;
            }

            collisions++;

            /* Sample target nuclide */
            alea_nuc_nuclide_t* target = NULL;
            alea_nuc_sample_nuclide(water, energy, rng(), &target);
            if (!target) { absorbed++; break; }

            /* Sample reaction */
            int mt;
            alea_nuc_sample_reaction(target, energy, rng(), &mt);

            /* Check for absorption (capture, fission, etc.) */
            if (mt != 2 && mt < 51) {
                /* Non-elastic, non-level scattering — treat as absorption for simplicity */
                absorbed++;
                break;
            }

            /* Sample collision outcome */
            double xi[3] = {rng(), rng(), rng()};
            alea_nuc_interaction_t result;
            alea_nuc_sample_collision(target, mt, energy, xi, &result);

            /* Update energy and direction */
            energy = result.energy_out;
            if (energy <= 1e-11) { absorbed++; break; }

            /* Update direction (simplified 1D: rotate mu_dir by scattering cosine) */
            /* In 1D slab, new direction = old_mu * mu_cm + sin(theta)*cos(phi)*sin(old_theta) */
            /* Simplified: just use lab cosine directly for 1D */
            double mu_scat = result.mu;
            double sin_old = sqrt(1.0 - mu_dir * mu_dir);
            if (sin_old > 1e-10) {
                double phi = 2.0 * 3.14159265358979 * rng();
                mu_dir = mu_dir * mu_scat + sin_old * sqrt(1.0 - mu_scat * mu_scat) * cos(phi);
            } else {
                mu_dir = (mu_dir > 0.0) ? mu_scat : -mu_scat;
            }
            if (mu_dir > 1.0) mu_dir = 1.0;
            if (mu_dir < -1.0) mu_dir = -1.0;

            /* Russian roulette for low weight */
            (void)weight;

            /* Kill if too many collisions (thermalized) */
            if (collisions > 10000) { absorbed++; break; }
        }

        avg_collisions += collisions;
    }

    avg_collisions /= n_histories;

    printf("Results:\n");
    printf("  Transmitted:  %d (%.2f%%)\n", transmitted, 100.0 * transmitted / n_histories);
    printf("  Reflected:    %d (%.2f%%)\n", reflected, 100.0 * reflected / n_histories);
    printf("  Absorbed:     %d (%.2f%%)\n", absorbed, 100.0 * absorbed / n_histories);
    printf("  Avg collisions: %.1f\n", avg_collisions);

    printf("\nTransmitted energy spectrum:\n");
    const char* labels[] = {"<1e-4", "1e-4", "1e-3", "0.01", "0.1", "1", "10", "100", "1e3", ">1e4"};
    for (int i = 0; i < energy_bins; i++) {
        if (leak_spectrum[i] > 0)
            printf("  %s MeV: %d\n", labels[i], leak_spectrum[i]);
    }

    printf("\nCross sections at 1 MeV:\n");
    printf("  Σ_total = %.4f cm⁻¹\n", alea_nuc_mat_xs_total(water, 1.0));
    printf("  MFP     = %.3f cm\n", alea_nuc_mean_free_path(water, 1.0));

    free(leak_spectrum);
    alea_nuc_material_destroy(water);
    alea_nuc_nuclide_free(h1);
    alea_nuc_nuclide_free(o16);
    alea_nuc_xsdir_free(xsdir);

    return 0;
}
