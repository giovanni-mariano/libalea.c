// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/* Small neutron adjoint example with manufactured one-group data.
 * A real calculation can replace the group constants with ACE data collapsed
 * by alea_nuc_mg_collapse(), using the same material-building call below. */
#include "alea_adjoint.h"
#include <math.h>
#include <stdio.h>

int main(void) {
    alea_system_t* sys = alea_create();
    if (!sys) return 1;
    int surface = alea_sphere_surface(sys, 1, 0, 0, 0, 1);
    int material = alea_add_material(sys, 1);
    if (surface < 0 || material < 0 ||
        alea_surface_set_boundary(sys, 1, ALEA_BOUNDARY_VACUUM) != 0 ||
        alea_add_cell(sys, 1, alea_halfspace(sys, surface, -1),
                      material, 1.0, 0) != 0) {
        alea_destroy(sys);
        return 1;
    }

    double bounds[] = {2.0, 0.0}; /* descending MeV */
    alea_nuc_multigroup_t* hydrogen = alea_nuc_mg_create(1, bounds);
    if (!hydrogen) { alea_destroy(sys); return 1; }
    hydrogen->sigma_t[0] = 10.0; /* illustrative microscopic barns */
    const alea_nuc_multigroup_t* nuclides[] = {hydrogen};
    double densities[] = {0.1}; /* atoms/(barn cm), so Sigma_t = 1/cm */
    double total[1], transfer[1];
    alea_error_t err = alea_mg_neutron_material_build(
        nuclides, densities, 1, 1, total, transfer);
    alea_nuc_mg_destroy(hydrogen);
    if (err != ALEA_OK) { alea_destroy(sys); return 1; }

    alea_mg_material_t cell_materials[] = {{total, transfer}};
    double physical_source[] = {2.0}; /* neutrons/(s cm^3), all directions */
    alea_adjoint_box_detector_t detector = {
        .lower = {-0.25, -0.25, -0.25},
        .upper = { 0.25,  0.25,  0.25},
        .group = 0, .response = 0.5 /* 1/cm */
    };
    alea_adjoint_problem_t problem = {
        .particle = ALEA_NUC_PARTICLE_NEUTRON,
        .n_groups = 1, .cell_count = 1, .cell_materials = cell_materials,
        .physical_source = physical_source,
        .detector_sampler = alea_adjoint_sample_box_detector,
        .detector_context = &detector
    };
    alea_adjoint_options_t options = {
        .histories = 20000, .seed = 12345,
        .max_events_per_history = 100, .max_segment_distance = 2.0
    };
    alea_adjoint_result_t result;
    alea_adjoint_failure_t failure;
    err = alea_adjoint_run(sys, &problem, &options, &result, &failure);
    if (err == ALEA_OK)
        printf("detector response: %.6g +/- %.2g per second\n",
               result.mean, result.standard_error);
    else
        fprintf(stderr, "adjoint failed at history %u: %s\n",
                failure.history_id, alea_get_error_detail());
    alea_destroy(sys);
    return err == ALEA_OK ? 0 : 1;
}
