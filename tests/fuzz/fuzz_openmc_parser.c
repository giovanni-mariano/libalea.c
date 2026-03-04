// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/**
 * @file fuzz_openmc_parser.c
 * @brief Fuzz target for OpenMC XML parser
 *
 * Build: make fuzz-build
 * Run:   ./bin/fuzz_openmc corpus_openmc/ -max_len=65536
 *
 * Or manually:
 *   clang -g -O1 -fsanitize=fuzzer,address,undefined -I../include \
 *         fuzz_openmc_parser.c ../bin/libalea_full.a -lm -o fuzz_openmc
 */

#include <stdint.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <unistd.h>

#include "alea.h"
#include "alea_mcnp.h"
#include "alea_openmc.h"
#include "alea_raycast.h"

int LLVMFuzzerTestOneInput(const uint8_t *data, size_t size) {
    /* Skip empty or huge inputs */
    if (size == 0 || size > 1024 * 1024) return 0;

    /* OpenMC parser needs a file, so write to temp file */
    char tmpname[] = "/tmp/fuzz_openmc_XXXXXX";
    int fd = mkstemp(tmpname);
    if (fd < 0) return 0;

    if (write(fd, data, size) != (ssize_t)size) {
        close(fd);
        unlink(tmpname);
        return 0;
    }
    close(fd);

    /* Parse the XML */
    openmc_model_t *omc_model = openmc_load(tmpname);
    unlink(tmpname);

    if (!omc_model) return 0;
    alea_system_t *sys = omc_model->sys;

    /* ========================================================================
     * Phase 1: Indexing
     * ======================================================================== */
    alea_build_universe_index(sys);
    alea_build_spatial_index(sys);

    /* ========================================================================
     * Phase 2: Basic inspection
     * ======================================================================== */
    size_t n_cells = alea_cell_count(sys);
    size_t n_surfs = alea_surface_count(sys);
    (void)n_surfs;

    for (size_t i = 0; i < n_cells && i < 32; i++) {
        alea_cell_info_t info;
        alea_cell_get_info(sys, i, &info);
    }

    /* ========================================================================
     * Phase 3: Point queries
     * ======================================================================== */
    static const double test_points[][3] = {
        {0, 0, 0},
        {1, 2, 3},
        {-100, -100, -100},
        {100, 100, 100},
    };

    for (size_t i = 0; i < sizeof(test_points) / sizeof(test_points[0]); i++) {
        double x = test_points[i][0];
        double y = test_points[i][1];
        double z = test_points[i][2];

        alea_find_cell(sys, x, y, z);
        alea_material_at(sys, x, y, z);

        alea_cell_hit_t hits[16];
        alea_find_all_cells(sys, x, y, z, hits, 16);
    }

    /* ========================================================================
     * Phase 4: Overlap detection
     * ======================================================================== */
    int overlap_pairs[64];
    alea_find_overlaps(sys, overlap_pairs, 32);

    /* ========================================================================
     * Phase 5: Bounding computations
     * ======================================================================== */
    double cx, cy, cz, r;
    alea_compute_bounding_sphere(sys, 1.0, &cx, &cy, &cz, &r);

    /* ========================================================================
     * Phase 6: Raycast (if geometry is bounded)
     * ======================================================================== */
    if (r > 0 && r < 1e10) {
        alea_raycast_result_t *result = alea_raycast_result_create();
        if (result) {
            alea_raycast(sys, cx, cy, cz, 1, 0, 0, r * 3, result);

            size_t nseg = alea_raycast_segment_count(result);
            for (size_t j = 0; j < nseg && j < 8; j++) {
                double t_enter, t_exit;
                int cell_id, mat_id;
                double density;
                alea_raycast_segment_get(result, j, &t_enter, &t_exit,
                                         &cell_id, &mat_id, &density);
            }

            alea_raycast_result_destroy(result);
        }
    }

    /* ========================================================================
     * Phase 7: Flatten (on a copy)
     * ======================================================================== */
    if (n_cells > 0 && n_cells < 1000) {
        alea_system_t *copy = alea_clone(sys);
        if (copy) {
            alea_flatten(copy, 0);
            alea_destroy(copy);
        }
    }

    /* ========================================================================
     * Phase 8: Export roundtrip
     * ======================================================================== */
    FILE *devnull = fopen("/dev/null", "w");
    if (devnull) {
        openmc_export_system_stream(sys, devnull);
        fclose(devnull);
    }

    devnull = fopen("/dev/null", "w");
    if (devnull) {
        mcnp_export_system_stream(sys, devnull);
        fclose(devnull);
    }

    openmc_model_destroy(omc_model);
    return 0;
}
