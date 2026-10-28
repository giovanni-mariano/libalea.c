// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
// SPDX-License-Identifier: MPL-2.0

/* Bounded characteristic tracks through three unit-width layers. Prescribed
 * removal coefficients demonstrate uncollided attenuation, without nuclear
 * data or a scattering solver. Both directions must give exp(-1.15).
 */
#include "alea.h"
#include "alea_raycast.h"
#include <math.h>
#include <stdio.h>

int main(void) {
    const double origins[] = {0.25, 0, 0, 2.75, 0, 0};
    const double directions[] = {1, 0, 0, -1, 0, 0};
    const double sigma[] = {0.2, 0.4, 0.8}; /* inverse model length */
    const double expected[] = {0.75, 1, 0.75};
    alea_system_t* sys = alea_create();
    alea_raycast_batch_result_t* tracks = alea_raycast_batch_result_create();
    alea_ray_navigator_t* nav = NULL;
    alea_ray_coverage_slice_result_t* coverage =
        alea_ray_coverage_slice_result_create();
    int status = 1;
    if (!sys || !tracks || !coverage) goto done;
    int material = alea_add_material(sys, 1);
    if (material < 0) goto done;
    /* Same material, distinct cells: never merge away the physical regions. */
    for (int i = 0; i < 3; ++i) {
        int surface = alea_box_surface(sys, i+1, i, i+1, -1, 1, -1, 1);
        if (surface < 0 || alea_add_cell(sys, i+1,
                alea_halfspace(sys, surface, -1), material, -1, 0) < 0)
            goto done;
    }
    const alea_raycast_batch_options_t options = {
        .struct_size = sizeof(options),
        .fields = ALEA_RAY_BATCH_MATERIAL | ALEA_RAY_BATCH_FULL_PATHS,
        .projected_depth = -1, .max_segments = 6, .max_path_entries = 6,
        .max_output_bytes = 65536
    };
    if (alea_raycast_hier_batch(sys, origins, directions, 2, 2.5,
                               &options, tracks) != 0) goto done;
    nav = alea_ray_navigator_create(sys);
    if (!nav) goto done;
    const uint64_t* offsets = alea_raycast_batch_ray_offsets(tracks);
    const double* enter = alea_raycast_batch_t_enter(tracks);
    const double* exit = alea_raycast_batch_t_exit(tracks);
    const int32_t* cells = alea_raycast_batch_cell_ids(tracks);
    const uint64_t* paths = alea_raycast_batch_segment_path_offsets(tracks);
    const uint64_t* keys = alea_raycast_batch_path_occurrence_keys(tracks);
    for (int ray = 0; ray < 2; ++ray) {
        /* Selected tracks alone do not establish absence of hidden overlaps. */
        alea_ray_coverage_slice_options_t check;
        alea_ray_coverage_slice_options_init(&check);
        check.t_max = 2.5;
        check.max_intervals = 16;
        check.max_owners = 16;
        check.max_output_bytes = 65536;
        if (alea_ray_coverage_query(sys, origins[3*ray], 0, 0,
                directions[3*ray], 0, 0, &check, coverage) != 0) goto done;
        const uint8_t* kinds = alea_ray_coverage_slice_kinds(coverage);
        for (size_t i = 0; i < alea_ray_coverage_slice_interval_count(coverage); ++i)
            if (kinds[i] != ALEA_RAY_COVERAGE_UNIQUE) goto done;
        alea_nav_location_t location;
        if (alea_ray_navigator_restart(nav, origins+3*ray, directions+3*ray,
                                       &location) != 0) goto done;
        if (offsets[ray+1] - offsets[ray] != 3) goto done;
        double optical_length = 0, total = 0;
        for (size_t j = offsets[ray]; j < offsets[ray+1]; ++j) {
            int cell = ray ? 3 - (int)(j-offsets[ray]) : 1 + (int)(j-offsets[ray]);
            double length = exit[j] - enter[j];
            if (cells[j] != cell || fabs(length-expected[cell-1]) > 1e-10 ||
                paths[j+1] - paths[j] != 1 || !keys[paths[j]]) goto done;
            alea_nav_event_t event;
            if (alea_ray_navigator_advance(nav, INFINITY, length, &event) != 0 ||
                event.before.kind != ALEA_NAV_MATERIAL ||
                event.before.cell_id != cell ||
                event.before.occurrence_key != keys[paths[j]] ||
                fabs(event.distance-length) > 1e-10) goto done;
            optical_length += sigma[cell-1] * length;
            total += length;
        }
        double attenuation = exp(-optical_length);
        if (fabs(total-2.5) > 1e-10 || fabs(attenuation-exp(-1.15)) > 1e-10)
            goto done;
        printf("track %d: length %.6f, attenuation %.12f\n", ray, total, attenuation);
    }
    status = 0;
done:
    if (status) fprintf(stderr, "characteristic track validation failed\n");
    alea_ray_navigator_destroy(nav);
    alea_ray_coverage_slice_result_destroy(coverage);
    alea_raycast_batch_result_destroy(tracks);
    alea_destroy(sys);
    return status;
}
