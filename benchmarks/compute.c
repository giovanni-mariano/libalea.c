// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
// SPDX-License-Identifier: MPL-2.0

#define _POSIX_C_SOURCE 200809L
#include "alea.h"
#include "alea_compute.h"
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#ifdef _WIN32
#include <windows.h>
#endif

static double seconds(void) {
#ifdef _WIN32
    LARGE_INTEGER counter, frequency;
    QueryPerformanceCounter(&counter);
    QueryPerformanceFrequency(&frequency);
    return (double)counter.QuadPart / (double)frequency.QuadPart;
#else
    struct timespec t;
    clock_gettime(CLOCK_MONOTONIC, &t);
    return (double)t.tv_sec + (double)t.tv_nsec * 1e-9;
#endif
}

static alea_system_t* scene_create(void) {
    alea_system_t* sys = alea_create();
    if (!sys) return NULL;
    int mat = alea_add_material(sys, 7);
    if (mat < 0) { alea_destroy(sys); return NULL; }
    for (int x = 0; x < 4; ++x)
        for (int y = 0; y < 4; ++y)
            for (int z = 0; z < 4; ++z) {
                int surface = alea_sphere_surface(sys, 0, 2*x-3, 2*y-3, 2*z-3, .85);
                if (surface < 0 || alea_add_cell(sys, 0,
                    alea_halfspace(sys, surface, -1), mat, -1, 0) < 0) {
                    alea_destroy(sys); return NULL;
                }
            }
    return sys;
}

static double measure(const alea_system_t* sys, alea_compute_backend_t backend,
                       const double* xyz, size_t count, alea_compute_result_t* output) {
    const char* name = backend == ALEA_COMPUTE_CPU ? "CPU FP64" : "Metal FP32";
    double start = seconds();
    alea_compute_t* ctx = alea_compute_create(sys, backend,
        backend == ALEA_COMPUTE_CPU ? ALEA_COMPUTE_FP64 : ALEA_COMPUTE_FP32);
    if (!ctx) { fprintf(stderr, "%s\n", alea_error()); return -1; }
    double setup = seconds() - start;
    double first = 0, total = 0;
    for (int run = 0; run < 6; ++run) {
        start = seconds();
        if (alea_compute_classify(ctx, xyz, count, output)) {
            fprintf(stderr, "%s\n", alea_error());
            alea_compute_destroy(ctx); return -1;
        }
        double elapsed = seconds() - start;
        if (run == 0) first = elapsed;
        else total += elapsed;
    }
    double mean = total / 5;
    printf("%-10s setup=%.6fs first=%.6fs mean_batch=%.6fs points/s=%.0f\n",
           name, setup, first, mean, count / mean);
    alea_compute_destroy(ctx);
    return mean;
}

int main(int argc, char** argv) {
    size_t count = 262144;
    if (argc > 1) {
        char* end;
        errno = 0;
        unsigned long long n = strtoull(argv[1], &end, 10);
        if (argv[1][0] == '-' || !argv[1][0] || *end || errno || !n ||
            n > UINT32_MAX || n > SIZE_MAX / (3 * sizeof(double))) {
            fprintf(stderr, "Usage: %s [positive point count]\n", argv[0]);
            return 1;
        }
        count = (size_t)n;
    }
    alea_system_t* sys = scene_create();
    double* xyz = malloc(count * 3 * sizeof(double));
    alea_compute_result_t* cpu = malloc(count * sizeof(*cpu));
    alea_compute_result_t* gpu = malloc(count * sizeof(*gpu));
    int status = 1;
    if (!sys || !xyz || !cpu || !gpu) { fprintf(stderr, "Allocation failed\n"); goto done; }
    uint32_t random = 12345;
    for (size_t i = 0; i < count * 3; ++i) {
        random = random * UINT32_C(1664525) + UINT32_C(1013904223);
        xyz[i] = 10.0 * ((double)random / 4294967296.0) - 5.0;
    }
    printf("64 disjoint spheres, %zu points; full cell scan on both backends\n", count);
    printf("Batch timings include validation, input conversion, dispatch and result copy.\n");
    double cpu_time = measure(sys, ALEA_COMPUTE_CPU, xyz, count, cpu);
    if (cpu_time < 0) goto done;
#ifdef ALEA_USE_METAL
    double gpu_time = measure(sys, ALEA_COMPUTE_METAL, xyz, count, gpu);
    if (gpu_time < 0) goto done;
    size_t differences = 0;
    for (size_t i = 0; i < count; ++i)
        if (cpu[i].owner_count != gpu[i].owner_count ||
            cpu[i].cell_id != gpu[i].cell_id || cpu[i].material_id != gpu[i].material_id)
            ++differences;
    printf("CPU/Metal batch time ratio: %.2fx; FP32/FP64 differing points: %zu/%zu\n",
           cpu_time / gpu_time, differences, count);
#endif
    status = 0;
done:
    free(gpu); free(cpu); free(xyz); alea_destroy(sys);
    return status;
}
