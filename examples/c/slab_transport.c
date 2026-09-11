// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/**
 * @file slab_transport.c
 * @brief Restricted continuous-energy neutron collision example
 *
 * Tracks 1 MeV neutrons through a one-dimensional hydrogen slab. Geometry
 * traversal and history control remain in the application; libalea evaluates
 * the material, samples flights, and performs stationary-target elastic or
 * absorption collisions.
 *
 * Usage:
 *   ./slab_transport <path-to-ace-dir> [thickness_cm] [histories] [zaid]
 */

#include "alea_nucdata.h"

#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

static uint64_t rng_state;

static double random_uniform(void* context) {
    uint64_t* state = context;
    *state = *state * UINT64_C(6364136223846793005) +
             UINT64_C(1442695040888963407);
    return (double)(*state >> 11) / 9007199254740992.0;
}

int main(int argc, char* argv[]) {
    const char* datadir = "test/fendl-FENDL-3.2c-neutron-ace/neutron/ace";
    const char* zaid = "1001.32c";
    double thickness = 5.0;
    int histories = 100000;
    if (argc > 1) datadir = argv[1];
    if (argc > 2) thickness = atof(argv[2]);
    if (argc > 3) histories = atoi(argv[3]);
    if (argc > 4) zaid = argv[4];
    if (!(thickness > 0.0) || histories <= 0) {
        fprintf(stderr, "thickness and history count must be positive\n");
        return 2;
    }
    rng_state = (uint64_t)time(NULL);

    alea_nuc_xsdir_t* xsdir = alea_nuc_xsdir_load_dir(datadir);
    if (!xsdir) {
        fprintf(stderr, "failed to load nuclear-data directory: %s\n", datadir);
        return 1;
    }
    alea_nuc_nuclide_t* hydrogen = alea_nuc_load_nuclide(xsdir, zaid);
    if (!hydrogen) {
        fprintf(stderr, "failed to load %s\n", zaid);
        alea_nuc_xsdir_free(xsdir);
        return 1;
    }

    alea_nuc_material_t* material = alea_nuc_material_create();
    if (!material || alea_nuc_material_add(material, hydrogen, 6.676e-2) != ALEA_OK) {
        fprintf(stderr, "failed to build hydrogen material\n");
        alea_nuc_material_destroy(material);
        alea_nuc_nuclide_free(hydrogen);
        alea_nuc_xsdir_free(xsdir);
        return 1;
    }

    alea_nuc_prepare_requirements_t requirements = {
        .required_capabilities = ALEA_NUC_CAP_RESTRICTED_NEUTRON
    };
    alea_nuc_capability_report_t report;
    alea_nuc_prepared_material_t* prepared = NULL;
    alea_error_t err = alea_nuc_prepare_material(
        material, &requirements, &report, &prepared);
    if (err != ALEA_OK) {
        fprintf(stderr, "data are not supported by this example: %s",
                report.detail);
        if (report.mt) fprintf(stderr, " (MT %d)", report.mt);
        fputc('\n', stderr);
        alea_nuc_material_destroy(material);
        alea_nuc_nuclide_free(hydrogen);
        alea_nuc_xsdir_free(xsdir);
        return 1;
    }

    int transmitted = 0, reflected = 0, absorbed = 0;
    double collision_sum = 0.0;
    for (int history = 0; history < histories; history++) {
        double x = 0.0;
        int collisions = 0;
        alea_nuc_particle_state_t neutron = {
            .type = ALEA_NUC_PARTICLE_NEUTRON,
            .energy = 1.0,
            .direction = {1.0, 0.0, 0.0},
            .weight = 1.0,
            .time = 0.0
        };

        for (;;) {
            alea_nuc_evaluation_t evaluation;
            err = alea_nuc_evaluate(prepared, &neutron, &evaluation);
            if (err != ALEA_OK) {
                absorbed++;
                break;
            }
            double distance;
            err = alea_nuc_sample_flight(&evaluation, random_uniform,
                                         &rng_state, &distance);
            if (err != ALEA_OK) {
                fprintf(stderr, "flight sampling failed: %s\n",
                        alea_error_string(err));
                alea_nuc_prepared_material_free(prepared);
                alea_nuc_material_destroy(material);
                alea_nuc_nuclide_free(hydrogen);
                alea_nuc_xsdir_free(xsdir);
                return 1;
            }
            x += distance * neutron.direction[0];
            if (x >= thickness) {
                transmitted++;
                break;
            }
            if (x <= 0.0) {
                reflected++;
                break;
            }

            alea_nuc_collision_result_t collision;
            err = alea_nuc_collide(&evaluation, random_uniform,
                                   &rng_state, &collision);
            if (err != ALEA_OK) {
                fprintf(stderr, "collision sampling failed: %s\n",
                        alea_error_string(err));
                alea_nuc_prepared_material_free(prepared);
                alea_nuc_material_destroy(material);
                alea_nuc_nuclide_free(hydrogen);
                alea_nuc_xsdir_free(xsdir);
                return 1;
            }
            collisions++;
            if (collision.outcome == ALEA_NUC_OUTCOME_ABSORBED) {
                absorbed++;
                break;
            }
            neutron = collision.outgoing;
            if (neutron.energy <= hydrogen->energy[0] || collisions >= 10000) {
                absorbed++;
                break;
            }
        }
        collision_sum += collisions;
    }

    printf("restricted hydrogen slab, %s\n", zaid);
    printf("histories: %d, thickness: %.6g cm\n", histories, thickness);
    printf("transmitted: %d (%.3f%%)\n", transmitted,
           100.0 * transmitted / histories);
    printf("reflected:   %d (%.3f%%)\n", reflected,
           100.0 * reflected / histories);
    printf("absorbed:    %d (%.3f%%)\n", absorbed,
           100.0 * absorbed / histories);
    printf("mean collisions: %.6g\n", collision_sum / histories);

    alea_nuc_prepared_material_free(prepared);
    alea_nuc_material_destroy(material);
    alea_nuc_nuclide_free(hydrogen);
    alea_nuc_xsdir_free(xsdir);
    return 0;
}
