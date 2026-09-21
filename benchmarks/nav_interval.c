// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
// SPDX-License-Identifier: MPL-2.0

/* Build after make lib-core modules:
 * cc -O2 -Iinclude -Isrc benchmarks/nav_interval.c \
 *   bin/libalea_mcnp.a bin/libalea.a -lm -pthread -fopenmp \
 *   -o bin/nav_interval_bench
 * Run from the repository root for the repeating-lattice fixture. */
#define _POSIX_C_SOURCE 200809L
#include "alea.h"
#include "alea_mcnp.h"
#include "alea_raycast.h"
#include <math.h>
#include <stdio.h>
#include <time.h>

static double seconds(void) {
    struct timespec t;
    clock_gettime(CLOCK_MONOTONIC, &t);
    return (double)t.tv_sec + 1e-9 * t.tv_nsec;
}

static int run(alea_system_t* sys, const char* geometry, int histories,
               int flights, double distance) {
    const double source[3] = {geometry[0] == 's' ? 0.0 : 0.5,
                              geometry[0] == 's' ? 0.0 : 0.5, 0.0};
    const double direction[3] = {1, 0, 0};
    const alea_nav_validation_mode_t modes[] = {
        ALEA_NAV_VALIDATE_FAST, ALEA_NAV_VALIDATE_STRICT,
        ALEA_NAV_VALIDATE_INTERVAL
    };
    const char* names[] = {"FAST", "STRICT", "INTERVAL"};
    for (int mode = 0; mode < 3; mode++) {
        alea_ray_navigator_t* nav = alea_ray_navigator_create(sys);
        if (!nav || alea_ray_navigator_set_validation_mode(nav, modes[mode]))
            return 1;
        alea_nav_location_t location;
        if (alea_ray_navigator_restart(nav, source, direction, &location))
            return 1;
        double start = seconds();
        for (int h = 0; h < histories; h++) {
            if (alea_ray_navigator_restart(nav, source, direction, &location))
                return 1;
            for (int f = 0; f < flights; f++) {
                alea_nav_event_t event;
                if (alea_ray_navigator_advance(nav, INFINITY, distance,
                                                &event)) {
                    fprintf(stderr, "%s %s: %s\n", geometry, names[mode],
                            alea_get_error_detail());
                    return 1;
                }
            }
        }
        double elapsed = seconds() - start;
        printf("%s %s: %.3f us/flight (%d flights)\n", geometry,
               names[mode], 1e6 * elapsed / (histories * flights),
               histories * flights);
        alea_ray_navigator_destroy(nav);
    }
    return 0;
}

int main(void) {
    alea_system_t* sphere = alea_create();
    int surface = alea_sphere_surface(sphere, 1, 0, 0, 0, 100);
    int material = alea_add_material(sphere, 1);
    if (surface < 0 || material < 0 || alea_add_cell(sphere, 1,
          alea_halfspace(sphere, surface, -1), material, -1, 0) < 0)
        return 1;
    if (run(sphere, "sphere", 10000, 20, 0.1)) return 1;
    alea_destroy(sphere);
    mcnp_model_t* lattice = mcnp_load("tests/data/mcnp_lattice_repeating.mcnp");
    if (!lattice) return 1;
    if (run(lattice->sys, "lattice", 5000, 10, 0.02)) return 1;
    mcnp_model_destroy(lattice);
    return 0;
}
