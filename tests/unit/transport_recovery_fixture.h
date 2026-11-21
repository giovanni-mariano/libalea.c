// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
// SPDX-License-Identifier: MPL-2.0
#ifndef ALEA_TEST_TRANSPORT_RECOVERY_FIXTURE_H
#define ALEA_TEST_TRANSPORT_RECOVERY_FIXTURE_H
#include "alea.h"

/* The reference extends the incident absorber to the outgoing void. */
static alea_system_t* recovery_absorber(bool gap) {
    alea_system_t* sys = alea_create();
    if (!sys) return NULL;
    const double outer_edge = 1.0 + 5e-11;
    int material = alea_add_material(sys, 1);
    int outer = alea_sphere_surface(sys, 3, 0, 0, 0, 3);
    int outgoing = alea_sphere_surface(sys, 2, 0, 0, 0, outer_edge);
    int incident = gap ? alea_sphere_surface(sys, 1, 0, 0, 0, 1) : outgoing;
    if (material < 0 || outer < 0 || incident < 0 || outgoing < 0 ||
        alea_material_add_nuclide(sys, material, 1001, ".80c", 1) ||
        alea_surface_set_boundary(sys, 3, ALEA_BOUNDARY_VACUUM) ||
        alea_add_cell(sys, 10, alea_halfspace(sys, outer, -1), ALEA_MATERIAL_VOID, 0, 0) != 0 ||
        alea_set_fill(sys, 0, 7, 0) ||
        alea_add_cell(sys, 11, alea_halfspace(sys, incident, -1), material, 0.02, 7) != 1 ||
        alea_add_cell(sys, 12, alea_halfspace(sys, outgoing, 1), ALEA_MATERIAL_VOID, 0, 7) != 2) {
        alea_destroy(sys);
        return NULL;
    }
    return sys;
}
#endif
