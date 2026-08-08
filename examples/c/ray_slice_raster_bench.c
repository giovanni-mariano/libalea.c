// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/*
 * Build (from the repository root):
 *   gcc -std=c11 -O3 -fopenmp -Iinclude -Isrc \
 *       examples/c/ray_slice_raster_bench.c bin/libalea.a -lm \
 *       -o bin/ray_slice_raster_bench
 *
 * Optional arguments: [nu] [nv] [iterations]
 */

#include "alea.h"
#include "alea_raycast.h"
#include "alea_slice.h"
#include "core/alea_system.h"

#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

static double now_seconds(void) {
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return (double)tv.tv_sec + 1e-6 * (double)tv.tv_usec;
}

static size_t parse_size(const char* text, size_t fallback) {
    if (!text) return fallback;
    char* end = NULL;
    unsigned long long value = strtoull(text, &end, 10);
    return end && *end == '\0' && value > 0 ? (size_t)value : fallback;
}

static alea_system_t* make_scene(void) {
    alea_system_t* sys = alea_create();
    if (!sys) return NULL;
    int material = alea_add_material(sys, 1);
    int sphere = alea_sphere_surface(sys, 1, 0.0, 0.0, 0.0, 4.0);
    if (material < 0 || sphere < 0 ||
        alea_add_cell(sys, 1, alea_surface_at(sys, sphere)->neg_node,
                      material, -1.0, 0) < 0) {
        alea_destroy(sys);
        return NULL;
    }
    return sys;
}

int main(int argc, char** argv) {
    const size_t nu = parse_size(argc > 1 ? argv[1] : NULL, 1024);
    const size_t nv = parse_size(argc > 2 ? argv[2] : NULL, 1024);
    const size_t iterations = parse_size(argc > 3 ? argv[3] : NULL, 5);
    if (nu > SIZE_MAX / nv ||
        nu * nv > SIZE_MAX / (sizeof(int32_t) + sizeof(int32_t) +
                              sizeof(double) + sizeof(uint8_t))) {
        fprintf(stderr, "raster dimensions overflow\n");
        return 2;
    }
    const size_t pixels = nu * nv;
    int32_t* cells = malloc(pixels * sizeof(*cells));
    int32_t* materials = malloc(pixels * sizeof(*materials));
    double* densities = malloc(pixels * sizeof(*densities));
    uint8_t* flags = malloc(pixels * sizeof(*flags));
    alea_system_t* sys = make_scene();
    if (!cells || !materials || !densities || !flags || !sys) {
        fprintf(stderr, "allocation or scene setup failed\n");
        free(cells); free(materials); free(densities); free(flags);
        alea_destroy(sys);
        return 1;
    }

    alea_slice_view_t view;
    alea_slice_view_axis(&view, 2, 0.0, -5.0, 5.0, -5.0, 5.0);
    alea_slice_raster_t raster;
    alea_slice_raster_init(&raster);
    raster.nu = nu;
    raster.nv = nv;
    raster.fields = ALEA_SLICE_RASTER_CELL_ID |
                    ALEA_SLICE_RASTER_MATERIAL_ID |
                    ALEA_SLICE_RASTER_DENSITY |
                    ALEA_SLICE_RASTER_RESOLUTION_FLAGS;
    raster.cell_ids = cells;
    raster.material_ids = materials;
    raster.densities = densities;
    raster.resolution_flags = flags;

    alea_raycast_batch_options_t trace_options = {
        .struct_size = sizeof(trace_options),
        .fields = ALEA_RAY_BATCH_MATERIAL | ALEA_RAY_BATCH_DENSITY |
                  ALEA_RAY_BATCH_RESOLUTION_FLAGS,
        .projected_depth = -1
    };
    double trace_seconds = 0.0;
    double raster_seconds = 0.0;
    size_t segment_count = 0;
    size_t compact_completed = 0;
    alea_raycast_batch_result_t* compact = NULL;
    for (size_t i = 0; i < iterations; i++) {
        compact = alea_raycast_batch_result_create();
        if (!compact) break;
        const double start = now_seconds();
        int rc = alea_trace_ray_slice_compact(sys, &view, nv, &trace_options, compact);
        trace_seconds += now_seconds() - start;
        if (rc != 0) break;
        segment_count = alea_raycast_batch_segment_count(compact);
        const double raster_start = now_seconds();
        rc = alea_rasterize_ray_slice_compact(&view, compact, &raster);
        raster_seconds += now_seconds() - raster_start;
        alea_raycast_batch_result_destroy(compact);
        compact = NULL;
        if (rc != 0) break;
        compact_completed++;
    }
    alea_raycast_batch_result_destroy(compact);
    if (compact_completed != iterations) {
        fprintf(stderr, "compact trace or rasterization failed: %s\n", alea_error());
        free(cells); free(materials); free(densities); free(flags);
        alea_destroy(sys);
        return 1;
    }

    alea_slice_raster_options_t fused_options;
    alea_slice_raster_options_init(&fused_options);
    double fused_seconds = 0.0;
    size_t fused_completed = 0;
    for (size_t i = 0; i < iterations; i++) {
        const double start = now_seconds();
        if (alea_trace_ray_slice_raster(sys, &view, &fused_options, &raster) != 0) break;
        fused_seconds += now_seconds() - start;
        fused_completed++;
    }
    if (fused_completed != iterations) {
        fprintf(stderr, "fused trace and rasterization failed: %s\n", alea_error());
        free(cells); free(materials); free(densities); free(flags);
        alea_destroy(sys);
        return 1;
    }

    printf("ray-slice raster benchmark: %zux%zu, %zu iterations\n", nu, nv, iterations);
    printf("  compact trace:       %.3f ms\n", 1e3 * trace_seconds / (double)iterations);
    printf("  raster from compact: %.3f ms\n", 1e3 * raster_seconds / (double)iterations);
    printf("  fused trace+raster:  %.3f ms\n", 1e3 * fused_seconds / (double)iterations);
    printf("  compact segments:    %zu\n", segment_count);
    printf("  raster field bytes:  %zu\n",
           pixels * (sizeof(*cells) + sizeof(*materials) + sizeof(*densities) + sizeof(*flags)));

    free(cells); free(materials); free(densities); free(flags);
    alea_destroy(sys);
    return 0;
}
