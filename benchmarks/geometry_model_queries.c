// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
// SPDX-License-Identifier: MPL-2.0

/* Representative-model validation, using public APIs only.
 * Usage: geometry_model_queries model rows frames xmin xmax ymin ymax z
 *        [navigation_sources navigation_steps]
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

static uint64_t track_digest(const alea_raycast_batch_result_t* result, int paths) {
    size_t n = alea_raycast_batch_segment_count(result);
    size_t rows = alea_raycast_batch_ray_count(result);
    uint64_t hash = UINT64_C(14695981039346656037);
#define HASH(accessor, count) hash = digest(hash, accessor(result), (count)*sizeof(*accessor(result)))
    HASH(alea_raycast_batch_ray_offsets, rows+1);
    HASH(alea_raycast_batch_t_enter, n);
    HASH(alea_raycast_batch_t_exit, n);
    HASH(alea_raycast_batch_cell_ids, n);
    HASH(alea_raycast_batch_material_ids, n);
    HASH(alea_raycast_batch_densities, n);
    HASH(alea_raycast_batch_resolution_flags, n);
    if (paths) {
        size_t entries = alea_raycast_batch_path_entry_count(result);
        HASH(alea_raycast_batch_segment_path_offsets, n+1);
        HASH(alea_raycast_batch_path_cell_ids, entries);
        HASH(alea_raycast_batch_path_material_ids, entries);
        HASH(alea_raycast_batch_path_universe_ids, entries);
        HASH(alea_raycast_batch_path_fill_universes, entries);
        HASH(alea_raycast_batch_path_depths, entries);
        HASH(alea_raycast_batch_path_is_lattice, entries);
        HASH(alea_raycast_batch_path_lattice_origins_xyz, 3*entries);
        HASH(alea_raycast_batch_path_occurrence_keys, entries);
    }
#undef HASH
    return hash;
}

int main(int argc, char** argv) {
    if (argc != 9 && argc != 11) return 2;
    char* end;
    size_t rows = strtoul(argv[2], &end, 10);
    if (*end || !rows || rows > 4096) return 2;
    int frames = (int)strtol(argv[3], &end, 10);
    if (*end || frames < 1 || frames > 100) return 2;
    size_t max_sources = 16;
    int nav_steps = 32;
    if (argc == 11) {
        max_sources = strtoul(argv[9], &end, 10);
        if (*end || !max_sources || max_sources > 1024) return 2;
        nav_steps = (int)strtol(argv[10], &end, 10);
        if (*end || nav_steps < 1 || nav_steps > 10000) return 2;
    }
    double box[5];
    for (int i = 0; i < 5; ++i) {
        box[i] = strtod(argv[4+i], &end);
        if (*end || !isfinite(box[i])) return 2;
    }
    if (!(box[1] > box[0]) || !(box[3] > box[2])) return 2;
    double started = seconds();
    mcnp_model_t* model = mcnp_load(argv[1]);
    if (!model) { fprintf(stderr, "%s\n", alea_error()); return 1; }
    printf("load ms=%.6f peak_mib=%.3f\n", (seconds()-started)*1000, peak_mib());
    fflush(stdout);
    int status = 1;
    alea_raycast_batch_result_t* result = alea_raycast_batch_result_create();
    alea_ray_navigator_t* nav = NULL;
    double* origins = calloc(3*rows, sizeof(double));
    double* directions = calloc(3*rows, sizeof(double));
    int32_t* pixels = NULL;
    if (!result || !origins || !directions) goto done;
    started = seconds();
    if (alea_prepare_query_acceleration(model->sys)) goto done;
    printf("prepare ms=%.6f peak_mib=%.3f\n", (seconds()-started)*1000, peak_mib());
    fflush(stdout);
    double length = box[1]-box[0];
    for (size_t i = 0; i < rows; ++i) {
        origins[3*i] = i%2 ? box[1] : box[0];
        origins[3*i+1] = box[2] + (i+0.5)*(box[3]-box[2])/rows;
        origins[3*i+2] = box[4];
        directions[3*i] = i%2 ? -1 : 1;
    }
    for (int paths = 0; paths < 2; ++paths) {
        uint64_t ownership_hash = 0;
        for (int surfaces = 0; surfaces < 2; ++surfaces) {
            alea_raycast_batch_options_t options = {
                .struct_size = sizeof(options), .projected_depth = -1,
                .fields = ALEA_RAY_BATCH_MATERIAL | ALEA_RAY_BATCH_DENSITY |
                          ALEA_RAY_BATCH_RESOLUTION_FLAGS,
                .max_segments = 1000000, .max_path_entries = 2000000,
                .max_output_bytes = 256u*1024u*1024u
            };
            if (paths) options.fields |= ALEA_RAY_BATCH_FULL_PATHS;
            if (surfaces) options.fields |= ALEA_RAY_BATCH_SURFACES;
            for (int f = -1; f < frames; ++f) {
                started = seconds();
                if (alea_raycast_hier_batch(model->sys, origins, directions,
                                            rows, length, &options, result)) goto done;
                double ms = (seconds()-started)*1000;
                uint64_t hash = track_digest(result, paths);
                if (!surfaces && f == -1) ownership_hash = hash;
                if (hash != ownership_hash) {
                    fprintf(stderr, "ownership changes with enrichment/repetition\n");
                    goto done;
                }
                uint64_t surface_hash = UINT64_C(14695981039346656037);
                if (surfaces) {
                    size_t n = alea_raycast_batch_segment_count(result);
                    surface_hash = digest(surface_hash, alea_raycast_batch_enter_surface_ids(result), n*sizeof(int32_t));
                    surface_hash = digest(surface_hash, alea_raycast_batch_exit_surface_ids(result), n*sizeof(int32_t));
                }
                printf("tracks paths=%d surfaces=%d frame=%d ms=%.6f segments=%zu hash=%016" PRIx64 " surface_hash=%016" PRIx64 " peak_mib=%.3f\n",
                       paths, surfaces, f, ms, alea_raycast_batch_segment_count(result),
                       hash, surface_hash, peak_mib());
                fflush(stdout);
            }
        }
    }
    nav = alea_ray_navigator_create(model->sys);
    if (!nav) goto done;
    const uint64_t* offsets = alea_raycast_batch_ray_offsets(result);
    const double* enter = alea_raycast_batch_t_enter(result);
    const double* exit = alea_raycast_batch_t_exit(result);
    const int32_t* materials = alea_raycast_batch_material_ids(result);
    /* Start inside selected material, then test strict navigation with direction
     * changes. Report ambiguous sources explicitly, rather than treating as void. */
    started = seconds();
    uint64_t nav_hash = UINT64_C(14695981039346656037);
    size_t events = 0, ambiguous = 0, failures = 0, sources = 0;
    for (size_t row = 0; row < rows && sources < max_sources; ++row) {
        for (size_t j = offsets[row]; j < offsets[row+1] && sources < max_sources; ++j) {
            if (materials[j] <= 0 || !(exit[j]-enter[j] > 1e-6)) continue;
            ++sources;
            double t = (enter[j]+exit[j])*0.5;
            double position[] = {origins[3*row]+directions[3*row]*t,
                                 origins[3*row+1], origins[3*row+2]};
            double direction[] = {1, 0.125, 0.0625};
            alea_nav_location_t location;
            if (alea_ray_navigator_restart(nav, position, direction, &location)) {
                fprintf(stderr, "navigation restart source=%zu: %s\n", sources, alea_error());
                ++failures; continue;
            }
            nav_hash = digest(nav_hash, &location.kind, sizeof(location.kind));
            if (location.kind != ALEA_NAV_MATERIAL) { ++ambiguous; continue; }
            double distance = fmin(0.05, (exit[j]-enter[j])*0.125);
            for (int step = 0; step < nav_steps; ++step) {
                alea_nav_event_t event;
                if (alea_ray_navigator_advance(nav, distance, 2*distance, &event)) {
                    fprintf(stderr, "navigation advance source=%zu step=%d position=%.17g,%.17g,%.17g direction=%.17g,%.17g,%.17g: %s\n",
                            sources, step, position[0], position[1], position[2],
                            direction[0], direction[1], direction[2], alea_error());
                    ++failures; break;
                }
                ++events;
                nav_hash = digest(nav_hash, &event.kind, sizeof(event.kind));
                nav_hash = digest(nav_hash, &event.distance, sizeof(event.distance));
                nav_hash = digest(nav_hash, &event.surface_id, sizeof(event.surface_id));
                nav_hash = digest(nav_hash, &event.before.cell_id, sizeof(event.before.cell_id));
                nav_hash = digest(nav_hash, &event.before.occurrence_key, sizeof(event.before.occurrence_key));
                nav_hash = digest(nav_hash, &event.after.cell_id, sizeof(event.after.cell_id));
                nav_hash = digest(nav_hash, &event.after.kind, sizeof(event.after.kind));
                memcpy(position, event.position, sizeof(position));
                if (event.kind == ALEA_NAV_VACUUM || event.kind == ALEA_NAV_BOUNDARY_ACTION ||
                    event.after.kind != ALEA_NAV_MATERIAL) break;
                direction[0] = -direction[0];
                direction[1] = step%3 ? -0.25 : 0.125;
                if (alea_ray_navigator_set_direction(nav, direction)) { ++failures; break; }
            }
        }
    }
    printf("navigation ms=%.6f sources=%zu events=%zu ambiguous=%zu failures=%zu hash=%016" PRIx64 " peak_mib=%.3f\n",
           (seconds()-started)*1000, sources, events, ambiguous, failures, nav_hash, peak_mib());
    fflush(stdout);
    size_t count = 1024*rows;
    pixels = malloc(2*count*sizeof(*pixels));
    if (!pixels) goto done;
    alea_slice_view_t view;
    alea_slice_view_axis(&view, 2, box[4], box[0], box[1], box[2], box[3]);
    alea_slice_raster_t raster;
    alea_slice_raster_init(&raster);
    raster.nu = 1024; raster.nv = rows;
    raster.fields = ALEA_SLICE_RASTER_CELL_ID | ALEA_SLICE_RASTER_MATERIAL_ID;
    raster.cell_ids = pixels; raster.material_ids = pixels+count;
    for (int f = -1; f < frames; ++f) {
        started = seconds();
        if (alea_trace_ray_slice_raster(model->sys, &view, NULL, &raster)) goto done;
        double ms = (seconds()-started)*1000;
        uint64_t hash = digest(UINT64_C(14695981039346656037), pixels, 2*count*sizeof(*pixels));
        printf("raster frame=%d ms=%.6f hash=%016" PRIx64 " peak_mib=%.3f\n", f, ms, hash, peak_mib());
        fflush(stdout);
    }
    status = 0;
done:
    if (status) fprintf(stderr, "query failed: %s\n", alea_error());
    free(pixels); free(origins); free(directions);
    alea_ray_navigator_destroy(nav);
    alea_raycast_batch_result_destroy(result);
    mcnp_model_destroy(model);
    return status;
}
