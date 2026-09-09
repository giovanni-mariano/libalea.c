// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0
/* Build with -O3 -Iinclude -Isrc and link bin/libalea.a -lm -fopenmp.
 * Compare identical binaries against baseline and candidate libraries. */
#include "alea.h"
#include "core/alea_system.h"
#include "core/alea_eval.h"
#include <stdio.h>
#include <time.h>

int main(void) {
    alea_system_t *sys = alea_create();
    if (!sys) return 1;
    const alea_primitive_type_t types[] = {
        ALEA_PRIMITIVE_PLANE, ALEA_PRIMITIVE_SPHERE,
        ALEA_PRIMITIVE_QUADRIC, ALEA_PRIMITIVE_TORUS_Z};
    for (unsigned k = 0; k < sizeof(types)/sizeof(types[0]); ++k) {
        alea_primitive_data_t data = {0};
        if (k == 0) data.plane = (alea_plane_data_t){1, 2, 3, -4};
        if (k == 1) data.sphere.radius = 3;
        if (k == 2) {
            data.quadric.coeffs[0] = 1; data.quadric.coeffs[1] = 2;
            data.quadric.coeffs[2] = 3; data.quadric.coeffs[9] = -4;
        }
        if (k == 3) {
            data.torus.major_radius = 3; data.torus.minor_radius = 1;
            data.torus.axial_semiwidth_B = 2;
        }
        int8_t inverted = 0;
        uint32_t id = alea_get_or_create_primitive(sys, types[k], &data, &inverted);
        alea_node_id_t node = alea_add_primitive_node(sys, id, -1, inverted, 1);
        double sum = 0;
        clock_t begin = clock();
        for (unsigned i = 0; i < 10000000; ++i)
            sum += alea_evaluate_point(sys, node, (i % 101)*0.07, 0.2, 0.3);
        printf("type=%d seconds=%.6f checksum=%.17g\n", types[k],
               (double)(clock()-begin)/CLOCKS_PER_SEC, sum);
    }
    alea_destroy(sys);
    return 0;
}
