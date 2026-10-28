// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
// SPDX-License-Identifier: MPL-2.0

/* Complete-slice benchmark; row profiling also reads internal work counters.
 * Usage: geometry_slice_queries model size frames xmin xmax ymin ymax z split|fused
 * Row profiling: replace the final mode with rows samples row_ms total_ms.
 * Compare identical builds at baseline and experimental commits. Hashes detect
 * regressions, not independent geometric truth. Wall times exclude hashing.
 */
#define _POSIX_C_SOURCE 200809L
#ifdef __APPLE__
#define _DARWIN_C_SOURCE
#endif
#include "alea.h"
#include "alea_mcnp.h"
#include "alea_raycast.h"
#include "alea_slice.h"
#include <inttypes.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/resource.h>
#include <time.h>
#include <signal.h>
#include <sys/time.h>
#include "raycast/raycast.h"

static double seconds(void) {
    struct timespec t;
    clock_gettime(CLOCK_MONOTONIC, &t);
    return t.tv_sec + t.tv_nsec * 1e-9;
}

static double peak_mib(void) {
    struct rusage r;
    if (getrusage(RUSAGE_SELF, &r)) return -1;
#ifdef __APPLE__
    return r.ru_maxrss / 1048576.0;
#else
    return r.ru_maxrss / 1024.0;
#endif
}

static uint64_t digest(uint64_t hash, const void* data, size_t bytes) {
    const unsigned char* p = data;
    for (size_t i = 0; i < bytes; ++i) hash = (hash ^ p[i]) * UINT64_C(1099511628211);
    return hash;
}



/* Standalone one-worker diagnostic owns SIGALRM. Cancellation is cooperative;
 * a parent process must impose a hard timeout as well. */
static volatile sig_atomic_t row_alarm;
static void row_deadline(int sig) {
    (void)sig;
    row_alarm = 1;
    alea_interrupt();
}

static int profile_rows(alea_system_t* sys, const alea_slice_view_t* view,
                        size_t rows, size_t samples, double row_ms, double total_ms,
                        const alea_raycast_batch_options_t* options,
                        alea_raycast_batch_result_t* tracks) {
    struct sigaction action = {0}, previous;
    action.sa_handler = row_deadline;
    sigemptyset(&action.sa_mask);
    if (sigaction(SIGALRM, &action, &previous)) return 1;
    sigset_t alarm_mask, old_mask;
    sigemptyset(&alarm_mask); sigaddset(&alarm_mask, SIGALRM);
    sigprocmask(SIG_UNBLOCK, &alarm_mask, NULL);
    double began = seconds();
    size_t attempted = 0, completed = 0, cancelled = 0, errors = 0;
    for (size_t i = 0; i < samples; ++i) {
        double remaining = total_ms - (seconds()-began)*1000;
        if (remaining <= 0) break;
        size_t row = (2*i+1)*rows/(2*samples);
        double v = view->v_min + ((double)row+0.5)*
                   (view->v_max-view->v_min)/(double)rows;
        double origin[3], direction[3];
        for (int axis = 0; axis < 3; ++axis) {
            origin[axis] = view->plane.origin[axis] + view->u_min*view->plane.u_axis[axis] + v*view->plane.v_axis[axis];
            direction[axis] = view->plane.u_axis[axis];
        }
        double allowance = fmin(row_ms, remaining);
        struct itimerval timer = {0}, off = {0};
        long long us = (long long)ceil(allowance*1000);
        timer.it_value.tv_sec = (time_t)(us/1000000);
        timer.it_value.tv_usec = (suseconds_t)(us%1000000);
        alea_clear_interrupt(); row_alarm = 0;
        printf("row_begin sample=%zu row=%zu allowance_ms=%.3f\n",i,row,allowance); fflush(stdout);
        if (setitimer(ITIMER_REAL, &timer, NULL)) { ++errors; break; }
        double start = seconds();
        int rc = alea_raycast_hier_batch(sys, origin, direction, 1,
                    view->u_max-view->u_min, options, tracks);
        double ms = (seconds()-start)*1000;
        /* Disarm while blocked, then deliver any already-pending alarm before
         * clearing the flag. A late alarm cannot cancel the following row. */
        sigprocmask(SIG_BLOCK, &alarm_mask, &old_mask);
        setitimer(ITIMER_REAL, &off, NULL);
        sigprocmask(SIG_SETMASK, &old_mask, NULL);
        int requested = row_alarm;
        int interrupted = rc != 0 && alea_get_last_error() == ALEA_ERR_INTERRUPTED;
        alea_clear_interrupt();
        ++attempted;
        if (rc == 0) {
            ++completed;
            size_t n = alea_raycast_batch_segment_count(tracks);
            uint64_t hash = UINT64_C(14695981039346656037);
#define HASH_ROW(accessor) hash = digest(hash, accessor(tracks), n*sizeof(*accessor(tracks)))
            HASH_ROW(alea_raycast_batch_t_enter);
            HASH_ROW(alea_raycast_batch_t_exit);
            HASH_ROW(alea_raycast_batch_cell_ids);
            HASH_ROW(alea_raycast_batch_material_ids);
            HASH_ROW(alea_raycast_batch_densities);
            HASH_ROW(alea_raycast_batch_resolution_flags);
#undef HASH_ROW
            alea_raycast_batch_work_stats_t work;
            if (alea_raycast_batch_result_get_work_stats_internal(tracks, &work)) { ++errors; break; }
            printf("row sample=%zu row=%zu status=complete ms=%.6f alarm=%d segments=%zu intervals=%llu snapshots=%llu root_queries=%llu full_queries=%llu hash=%016" PRIx64 " peak_mib=%.3f\n",
                i,row,ms,requested,n,(unsigned long long)work.max_selected_intervals_yielded,
                (unsigned long long)work.max_path_snapshot_copies,
                (unsigned long long)work.max_owner_root_queries,
                (unsigned long long)work.max_owner_full_queries,hash,peak_mib());
        } else {
            if (interrupted && requested) ++cancelled; else ++errors;
            printf("row sample=%zu row=%zu status=%s ms=%.6f alarm=%d peak_mib=%.3f\n",
                i,row,interrupted && requested ? "cancelled" : "error",ms,requested,peak_mib());
        }
        fflush(stdout);
        if (errors) break;
    }
    sigaction(SIGALRM, &previous, NULL);
    printf("rows_summary requested=%zu attempted=%zu completed=%zu cancelled=%zu errors=%zu elapsed_ms=%.6f\n",
           samples,attempted,completed,cancelled,errors,(seconds()-began)*1000);
    fflush(stdout);
    return errors ? 1 : cancelled || attempted < samples ? 3 : 0;
}

int main(int argc, char** argv) {
    if (argc != 10 && argc != 13) return 2;
    char* end;
    size_t size = strtoul(argv[2], &end, 10);
    if (end == argv[2] || *end || !size || size > 2048) return 2;
    long frames = strtol(argv[3], &end, 10);
    if (end == argv[3] || *end || frames < 1 || frames > 30) return 2;
    double bounds[5];
    for (int i = 0; i < 5; ++i) {
        bounds[i] = strtod(argv[4+i], &end);
        if (end == argv[4+i] || *end || !isfinite(bounds[i])) return 2;
    }
    if (!(bounds[0] < bounds[1] && bounds[2] < bounds[3])) return 2;
    int fused = !strcmp(argv[9], "fused");
    int rows_mode = !strcmp(argv[9], "rows");
    size_t samples = 0;
    double row_ms = 0, total_ms = 0;
    if (rows_mode) {
        if (argc != 13 || frames != 1) return 2;
        samples = strtoul(argv[10], &end, 10);
        if (end == argv[10] || *end || !samples || samples > size) return 2;
        row_ms = strtod(argv[11], &end);
        if (end == argv[11] || *end || !isfinite(row_ms) || row_ms < 0.001 || row_ms > 60000) return 2;
        total_ms = strtod(argv[12], &end);
        if (end == argv[12] || *end || !isfinite(total_ms) || total_ms < 0.001 || total_ms > 300000) return 2;
    } else if (argc != 10 || (!fused && strcmp(argv[9], "split"))) return 2;
    /* Workers created during preparation inherit the blocked alarm. Only the
     * profiling thread unmasks it after installing its handler. */
    if (rows_mode) {
        sigset_t mask;
        sigemptyset(&mask); sigaddset(&mask, SIGALRM);
        if (sigprocmask(SIG_BLOCK, &mask, NULL)) return 1;
    }
    double start = seconds();
    mcnp_model_t* model = mcnp_load(argv[1]);
    if (!model) return 1;
    printf("load ms=%.6f peak_mib=%.3f\n", (seconds()-start)*1000, peak_mib());
    fflush(stdout);
    int status = 1;
    size_t count = size*size;
    int32_t* cells = malloc(count*sizeof(*cells));
    int32_t* materials = malloc(count*sizeof(*materials));
    double* density = malloc(count*sizeof(*density));
    uint8_t* flags = malloc(count*sizeof(*flags));
    alea_raycast_batch_result_t* tracks = alea_raycast_batch_result_create();
    if (!cells || !materials || !density || !flags || !tracks) goto done;
    start = seconds();
    if (alea_prepare_query_acceleration(model->sys)) goto done;
    printf("prepare ms=%.6f peak_mib=%.3f pixel_bytes=%zu trace_budget=%u\n",
           (seconds()-start)*1000, peak_mib(), count*17, 256u*1024u*1024u);
    fflush(stdout);
    alea_slice_view_t view;
    alea_slice_view_axis(&view, 2, bounds[4], bounds[0], bounds[1], bounds[2], bounds[3]);
    alea_slice_raster_t raster;
    alea_slice_raster_init(&raster);
    raster.nu = raster.nv = size;
    raster.fields = ALEA_SLICE_RASTER_CELL_ID | ALEA_SLICE_RASTER_MATERIAL_ID |
                    ALEA_SLICE_RASTER_DENSITY | ALEA_SLICE_RASTER_RESOLUTION_FLAGS;
    raster.cell_ids = cells; raster.material_ids = materials;
    raster.densities = density; raster.resolution_flags = flags;
    alea_raycast_batch_options_t trace = {
        .struct_size = sizeof(trace), .projected_depth = -1,
        .fields = ALEA_RAY_BATCH_MATERIAL | ALEA_RAY_BATCH_DENSITY |
                  ALEA_RAY_BATCH_RESOLUTION_FLAGS,
        .max_segments = 1000000, .max_output_bytes = 256u*1024u*1024u
    };
    alea_slice_raster_options_t fused_options;
    alea_slice_raster_options_init(&fused_options);
    fused_options.max_segments = trace.max_segments;
    fused_options.max_trace_output_bytes = trace.max_output_bytes;
    if (rows_mode) {
        status = profile_rows(model->sys, &view, size, samples, row_ms, total_ms, &trace, tracks);
        goto done;
    }
    uint64_t previous = 0;
    for (int frame = 0; frame < frames; ++frame) {
        printf("begin frame=%d mode=%s\n", frame, argv[9]); fflush(stdout);
        double trace_ms = 0, raster_ms = 0;
        start = seconds();
        if (fused) {
            if (alea_trace_ray_slice_raster(model->sys, &view, &fused_options, &raster)) goto done;
        } else {
            if (alea_trace_ray_slice_compact(model->sys, &view, size, &trace, tracks)) goto done;
            trace_ms = (seconds()-start)*1000;
            printf("traced frame=%d ms=%.6f segments=%zu peak_mib=%.3f\n", frame,
                   trace_ms, alea_raycast_batch_segment_count(tracks), peak_mib());
            fflush(stdout);
            double raster_start = seconds();
            if (alea_rasterize_ray_slice_compact(&view, tracks, &raster)) goto done;
            raster_ms = (seconds()-raster_start)*1000;
        }
        double total_ms = fused ? (seconds()-start)*1000 : trace_ms + raster_ms;
        uint64_t hash = UINT64_C(14695981039346656037);
        hash = digest(hash, cells, count*sizeof(*cells));
        hash = digest(hash, materials, count*sizeof(*materials));
        hash = digest(hash, density, count*sizeof(*density));
        hash = digest(hash, flags, count*sizeof(*flags));
        size_t material_pixels = 0, flagged_pixels = 0;
        for (size_t i = 0; i < count; ++i) {
            material_pixels += materials[i] > 0;
            flagged_pixels += flags[i] != 0;
        }
        if (frame && previous != hash) goto done;
        previous = hash;
        printf("slice frame=%d size=%zu trace_ms=%.6f raster_ms=%.6f total_ms=%.6f material_pixels=%zu flagged_pixels=%zu hash=%016" PRIx64 " peak_mib=%.3f\n",
               frame, size, trace_ms, raster_ms, total_ms, material_pixels, flagged_pixels, hash, peak_mib());
        fflush(stdout);
    }
    status = 0;
done:
    if (status == 1) fprintf(stderr, "slice failed: %s\n", alea_error());
    free(cells); free(materials); free(density); free(flags);
    alea_raycast_batch_result_destroy(tracks);
    mcnp_model_destroy(model);
    return status;
}
