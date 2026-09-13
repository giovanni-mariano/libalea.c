#define _POSIX_C_SOURCE 200809L

// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/** @file bench_transport.c
 *  @brief Pointwise transport-kernel and history-parallel throughput benchmark
 */

#include "alea_nucdata.h"
#include "util/alea_parallel.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

typedef struct {
    double value;
    uint64_t calls;
} fixed_rng_t;

typedef struct {
    alea_nuc_urr_sample_t* urr;
    double* component_total;
    double* component_elastic;
    double* component_thermal;
    double* reaction_rates;
    alea_nuc_evaluation_workspace_t workspace;
} worker_t;

typedef struct {
    const alea_nuc_prepared_material_t* prepared;
    const alea_nuc_particle_state_t* incident;
    worker_t* workers;
    double* output;
} batch_t;

static double wall_time(void) {
    struct timespec value;
    clock_gettime(CLOCK_MONOTONIC, &value);
    return value.tv_sec + value.tv_nsec * 1e-9;
}

static double fixed_uniform(void* context) {
    fixed_rng_t* rng = context;
    rng->calls++;
    return rng->value;
}

static int run_histories(void* opaque, size_t worker_index,
                         size_t begin, size_t end) {
    batch_t* batch = opaque;
    worker_t* worker = &batch->workers[worker_index];
    for (size_t i = begin; i < end; i++) {
        fixed_rng_t urr_rng = {0.5, 0};
        alea_nuc_evaluation_t evaluation;
        alea_error_t status = alea_nuc_evaluate_urr(
            batch->prepared, batch->incident, fixed_uniform, &urr_rng,
            &worker->workspace, &evaluation);
        double distance = NAN;
        fixed_rng_t flight_rng = {0.5, 0};
        if (status == ALEA_OK)
            status = alea_nuc_sample_flight(
                &evaluation, fixed_uniform, &flight_rng, &distance);
        fixed_rng_t collision_rng = {0.01, 0};
        alea_nuc_collision_result_t collision;
        if (status == ALEA_OK)
            status = alea_nuc_collide(
                &evaluation, fixed_uniform, &collision_rng, &collision);
        if (status != ALEA_OK) return 1;
        batch->output[i] = distance + collision.outgoing.energy;
    }
    return 0;
}

static int worker_init(worker_t* worker, size_t components, size_t reactions) {
    worker->urr = calloc(components, sizeof(*worker->urr));
    worker->component_total = calloc(components, sizeof(double));
    worker->component_elastic = calloc(components, sizeof(double));
    worker->component_thermal = calloc(components, sizeof(double));
    worker->reaction_rates = calloc(reactions ? reactions : 1, sizeof(double));
    if (!worker->urr || !worker->component_total ||
        !worker->component_elastic || !worker->component_thermal ||
        !worker->reaction_rates)
        return 0;
    worker->workspace = (alea_nuc_evaluation_workspace_t){
        .components=worker->urr,
        .capacity=components,
        .component_total=worker->component_total,
        .component_elastic=worker->component_elastic,
        .component_thermal=worker->component_thermal,
        .component_rate_capacity=components,
        .reaction_rates=worker->reaction_rates,
        .reaction_rate_capacity=reactions
    };
    return 1;
}

static void worker_free(worker_t* worker) {
    free(worker->urr);
    free(worker->component_total);
    free(worker->component_elastic);
    free(worker->component_thermal);
    free(worker->reaction_rates);
}

static double benchmark_batch(batch_t* batch, size_t histories,
                              size_t workers, size_t* actual_workers) {
    double start = wall_time();
    alea_parallel_status_t status = alea_parallel_for(
        histories, 256, workers, ALEA_PARALLEL_DYNAMIC,
        run_histories, batch, actual_workers);
    double elapsed = wall_time() - start;
    return status == ALEA_PARALLEL_OK ? elapsed : -1.0;
}

static int benchmark_serial_phases(
    const alea_nuc_prepared_material_t* prepared,
    const alea_nuc_particle_state_t* incident, worker_t* worker,
    size_t iterations) {
    fixed_rng_t rng = {0.5, 0};
    alea_nuc_evaluation_t evaluation;
    volatile double checksum = 0.0;

    double start = wall_time();
    for (size_t i = 0; i < iterations; i++) {
        alea_error_t status = alea_nuc_evaluate_urr(
            prepared, incident, fixed_uniform, &rng,
            &worker->workspace, &evaluation);
        if (status != ALEA_OK) return 0;
        checksum += evaluation.macro_total;
    }
    double evaluate_elapsed = wall_time() - start;

    start = wall_time();
    for (size_t i = 0; i < iterations; i++) {
        double distance;
        alea_error_t status = alea_nuc_sample_flight(
            &evaluation, fixed_uniform, &rng, &distance);
        if (status != ALEA_OK) return 0;
        checksum += distance;
    }
    double flight_elapsed = wall_time() - start;

    rng.value = 0.01; /* Select elastic scattering for this U-238 case. */
    start = wall_time();
    for (size_t i = 0; i < iterations; i++) {
        alea_nuc_collision_result_t collision;
        alea_error_t status = alea_nuc_collide(
            &evaluation, fixed_uniform, &rng, &collision);
        if (status != ALEA_OK) return 0;
        checksum += collision.outgoing.energy;
    }
    double collision_elapsed = wall_time() - start;

    printf("serial evaluation: %.3f ns/op, %.3f M/s\n",
           evaluate_elapsed * 1e9 / iterations,
           iterations / evaluate_elapsed / 1e6);
    printf("serial flight:     %.3f ns/op, %.3f M/s\n",
           flight_elapsed * 1e9 / iterations,
           iterations / flight_elapsed / 1e6);
    printf("serial collision:  %.3f ns/op, %.3f M/s\n",
           collision_elapsed * 1e9 / iterations,
           iterations / collision_elapsed / 1e6);
    fprintf(stderr, "phase checksum=%.17g\n", checksum);
    return 1;
}

static int report_batch(batch_t* batch, size_t histories, size_t requested,
                        double serial, size_t* actual_workers) {
    double elapsed = benchmark_batch(
        batch, histories, requested, actual_workers);
    if (elapsed < 0.0) return 0;
    printf("workers=%zu time=%.6f s throughput=%.3f Mhistory/s speedup=%.2f\n",
           *actual_workers, elapsed, histories / elapsed / 1e6,
           serial / elapsed);
    return 1;
}

int main(int argc, char** argv) {
    const char* path = argc > 1 ? argv[1] : getenv("ALEA_ENDFB80_XSDIR");
    size_t histories = argc > 2 ? (size_t)strtoull(argv[2], NULL, 10) : 1000000;
    if (!path || !path[0] || histories == 0) {
        fprintf(stderr, "usage: %s XSDIR [histories [max-workers]]\n", argv[0]);
        return 2;
    }

    alea_nuc_xsdir_t* xsdir = alea_nuc_xsdir_load(path);
    alea_nuc_nuclide_t* uranium = xsdir
        ? alea_nuc_xsdir_get_nuclide(xsdir, "92238.00c") : NULL;
    alea_nuc_material_t* material = alea_nuc_material_create();
    if (!uranium || !material ||
        alea_nuc_material_add(material, uranium, 0.01) != ALEA_OK) return 3;
    alea_nuc_prepare_requirements_t requirements = {
        .required_capabilities = ALEA_NUC_CAP_CONTINUOUS_NEUTRON |
                                 ALEA_NUC_CAP_URR
    };
    alea_nuc_capability_report_t report;
    alea_nuc_prepared_material_t* prepared = NULL;
    if (alea_nuc_prepare_material(material, &requirements, &report,
                                  &prepared) != ALEA_OK) return 4;
    size_t components, reactions;
    if (alea_nuc_evaluation_workspace_sizes(
            prepared, &components, &reactions) != ALEA_OK) return 5;

    size_t runtime_workers = alea_parallel_max_workers();
    size_t maximum_workers = argc > 3
        ? (size_t)strtoull(argv[3], NULL, 10) : runtime_workers;
    if (maximum_workers == 0 || maximum_workers > runtime_workers)
        maximum_workers = runtime_workers;
    worker_t* workers = calloc(maximum_workers, sizeof(*workers));
    double* output = calloc(histories, sizeof(*output));
    if (!workers || !output) return 6;
    for (size_t i = 0; i < maximum_workers; i++)
        if (!worker_init(&workers[i], components, reactions)) return 7;

    alea_nuc_particle_state_t incident = {
        ALEA_NUC_PARTICLE_NEUTRON, 0.05, {0.0, 0.0, 1.0}, 1.0, 0.0
    };
    batch_t batch = {prepared, &incident, workers, output};
    size_t actual;
    if (benchmark_batch(&batch, 1, 1, &actual) < 0.0) return 8;

    printf("U-238 at 50 keV, %zu histories\n", histories);
    printf("workspace: %zu components, %zu reactions, %zu bytes/worker\n",
           components, reactions,
           components * (sizeof(alea_nuc_urr_sample_t) + 3 * sizeof(double)) +
               reactions * sizeof(double));
    if (!benchmark_serial_phases(
            prepared, &incident, &workers[0], histories)) return 9;
    double serial = benchmark_batch(&batch, histories, 1, &actual);
    if (serial < 0.0) return 10;
    printf("workers=%zu time=%.6f s throughput=%.3f Mhistory/s speedup=1.00\n",
           actual, serial, histories / serial / 1e6);
    size_t last_requested = 1;
    for (size_t requested = 2; requested <= maximum_workers;
         requested *= 2) {
        if (!report_batch(&batch, histories, requested, serial, &actual))
            return 11;
        last_requested = requested;
        if (requested > maximum_workers / 2) break;
    }
    if (last_requested != maximum_workers &&
        !report_batch(&batch, histories, maximum_workers, serial, &actual))
        return 12;

    volatile double checksum = 0.0;
    for (size_t i = 0; i < histories; i++) checksum += output[i];
    fprintf(stderr, "checksum=%.17g\n", checksum);
    for (size_t i = 0; i < maximum_workers; i++) worker_free(&workers[i]);
    free(workers);
    free(output);
    alea_nuc_prepared_material_free(prepared);
    alea_nuc_material_destroy(material);
    alea_nuc_xsdir_free(xsdir);
    return 0;
}
