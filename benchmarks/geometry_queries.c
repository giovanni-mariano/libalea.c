// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
// SPDX-License-Identifier: MPL-2.0

/* Build with the same flags/library for each comparison:
 * cc -std=c11 -O2 -Iinclude benchmarks/geometry_queries.c bin/libalea.a
 *    -pthread -lm -o build/geometry_queries
 * Run: ALEA_NUM_THREADS=1 build/geometry_queries [capture_paths: 0 or 1]
 * Times include synchronous query output, excluding loading and preparation.
 */
#define _POSIX_C_SOURCE 200809L
#include "alea.h"
#include "alea_raycast.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

static double seconds(void) {
    struct timespec t;
    clock_gettime(CLOCK_MONOTONIC, &t);
    return t.tv_sec + t.tv_nsec * 1e-9;
}

int main(int argc, char** argv) {
    enum { LAYERS = 32, RAYS = 512, FRAMES = 7 };
    double origins[3 * RAYS], directions[3 * RAYS] = {0};
    alea_system_t* sys = alea_create();
    alea_raycast_batch_result_t* result = alea_raycast_batch_result_create();
    int status = 1;
    if (!sys || !result) goto done;
    int material = alea_add_material(sys, 1);
    if (material < 0) goto done;
    for (int i = 0; i < LAYERS; ++i) {
        int surface = alea_box_surface(sys, i + 1, i, i + 1, -2, 2, -2, 2);
        if (surface < 0 || alea_add_cell(sys, i + 1,
                alea_halfspace(sys, surface, -1), material, -1, 0) < 0)
            goto done;
    }
    double start = seconds();
    if (alea_prepare_query_acceleration(sys) != 0) goto done;
    printf("preparation_ms %.3f\n", (seconds() - start) * 1000);
    for (int i = 0; i < RAYS; ++i) {
        origins[3*i] = 0.25;
        origins[3*i+1] = -1 + 2.0 * i / RAYS;
        origins[3*i+2] = 0;
        directions[3*i] = 1;
    }
    alea_raycast_batch_options_t options = {
        .struct_size = sizeof(options), .fields = ALEA_RAY_BATCH_MATERIAL,
        .projected_depth = -1
    };
    if (argc > 1 && atoi(argv[1])) options.fields |= ALEA_RAY_BATCH_FULL_PATHS;
    for (int frame = -1; frame < FRAMES; ++frame) {
        start = seconds();
        if (alea_raycast_hier_batch(sys, origins, directions, RAYS,
                                   LAYERS - 0.5, &options, result) != 0)
            goto done;
        double elapsed = seconds() - start;
        const uint64_t* rows = alea_raycast_batch_ray_offsets(result);
        const double* enter = alea_raycast_batch_t_enter(result);
        const double* exit = alea_raycast_batch_t_exit(result);
        const int32_t* cells = alea_raycast_batch_cell_ids(result);
        for (int ray = 0; ray < RAYS; ++ray) {
            if (rows[ray+1] - rows[ray] != LAYERS) goto done;
            for (int i = 0; i < LAYERS; ++i) {
                size_t j = rows[ray] + i;
                double expected_enter = i ? i - 0.25 : 0;
                double expected_exit = i == LAYERS-1 ? LAYERS-0.5 : i+0.75;
                if (cells[j] != i+1 || fabs(enter[j]-expected_enter) > 1e-9 ||
                    fabs(exit[j]-expected_exit) > 1e-9) goto done;
            }
        }
        printf("%s %d ms %.3f segments %zu\n", frame < 0 ? "warmup" : "frame",
               frame, elapsed * 1000, alea_raycast_batch_segment_count(result));
    }
    status = 0;
done:
    if (status) fprintf(stderr, "geometry query benchmark failed\n");
    alea_raycast_batch_result_destroy(result);
    alea_destroy(sys);
    return status;
}
