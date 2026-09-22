// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

#include "alea_test.h"
#include "alea_adjoint.h"
#include <math.h>

static alea_error_t center_detector(void* context, uint64_t seed,
    uint32_t history, alea_adjoint_detector_particle_t* out) {
    (void)seed; (void)history;
    *out = (alea_adjoint_detector_particle_t){0};
    out->physical_direction[0] = 1.0;
    out->group = *(const size_t*)context;
    out->weight = 4.0 * acos(-1.0);
    return ALEA_OK;
}

TEST(adjoint_vacuum_and_absorber_match_path_integral_for_both_species) {
    for (int species = 0; species < 2; ++species) {
        alea_system_t* sys = alea_create();
        int sphere = alea_sphere_surface(sys, 1, 0, 0, 0, 1);
        ASSERT(sphere >= 0);
        ASSERT_EQ(alea_surface_set_boundary(sys, 1, ALEA_BOUNDARY_VACUUM), 0);
        int material = alea_add_material(sys, 5);
        ASSERT(material >= 0);
        ASSERT_EQ(alea_add_cell(sys, 1, alea_halfspace(sys, sphere, -1),
                                material, 1.0, 0), 0);
        double total[] = {0.0}, transfer[] = {0.0}, source[] = {2.0};
        alea_adjoint_material_t cells[] = {{total, transfer}};
        size_t detector_group = 0;
        alea_adjoint_problem_t p = {
            .particle = species ? ALEA_NUC_PARTICLE_PHOTON :
                                  ALEA_NUC_PARTICLE_NEUTRON,
            .n_groups = 1, .cell_count = 1, .cell_materials = cells,
            .physical_source = source, .detector_sampler = center_detector,
            .detector_context = &detector_group
        };
        alea_adjoint_options_t o = {
            .histories = 1000, .seed = 42, .max_events_per_history = 100,
            .max_segment_distance = 0.2
        };
        alea_adjoint_result_t result;
        alea_transport_failure_t failure;
        ASSERT_EQ(alea_adjoint_run(sys, &p, &o, &result, &failure), ALEA_OK);
        ASSERT_NEAR(result.mean, 2.0, 1e-8);
        ASSERT_NEAR(result.standard_error, 0.0, 1e-7);
        total[0] = 1.0;
        o.histories = 40000;
        ASSERT_EQ(alea_adjoint_run(sys, &p, &o, &result, &failure), ALEA_OK);
        ASSERT_NEAR(result.mean, 2.0*(1.0-exp(-1.0)), 0.015);
        ASSERT(result.standard_error > 0.0);
        alea_destroy(sys);
    }
}

TEST(adjoint_transposes_asymmetric_group_transfer_and_replays) {
    alea_system_t* sys = alea_create();
    int sphere = alea_sphere_surface(sys, 1, 0, 0, 0, 1);
    ASSERT(sphere >= 0);
    ASSERT_EQ(alea_surface_set_boundary(sys, 1, ALEA_BOUNDARY_VACUUM), 0);
    int material = alea_add_material(sys, 6);
    ASSERT(material >= 0);
    ASSERT_EQ(alea_add_cell(sys, 1, alea_halfspace(sys, sphere, -1),
                            material, 1.0, 0), 0);
    double total[] = {1.0, 1.0};
    double transfer[] = {0.0, 0.5, 0.0, 0.0};
    double source[] = {1.0, 0.0};
    alea_adjoint_material_t cells[] = {{total, transfer}};
    size_t detector_group = 1;
    alea_adjoint_problem_t p = {
        .particle = ALEA_NUC_PARTICLE_NEUTRON, .n_groups = 2,
        .cell_count = 1, .cell_materials = cells, .physical_source = source,
        .detector_sampler = center_detector,
        .detector_context = &detector_group
    };
    alea_adjoint_options_t o = {
        .histories = 20000, .seed = 123, .max_events_per_history = 100,
        .max_segment_distance = 0.25
    };
    alea_adjoint_result_t a, b;
    ASSERT_EQ(alea_adjoint_run(sys, &p, &o, &a, NULL), ALEA_OK);
    ASSERT(a.mean > 0.01);
    ASSERT_EQ(alea_adjoint_run(sys, &p, &o, &b, NULL), ALEA_OK);
    ASSERT_EQ(a.sum, b.sum);
    detector_group = 0;
    ASSERT_EQ(alea_adjoint_run(sys, &p, &o, &b, NULL), ALEA_OK);
    ASSERT(b.mean > a.mean);
    source[0] = 0.0;
    source[1] = 1.0;
    ASSERT_EQ(alea_adjoint_run(sys, &p, &o, &b, NULL), ALEA_OK);
    ASSERT_EQ(b.mean, 0.0);
    alea_destroy(sys);
}

TEST(adjoint_rejects_incomplete_histories_and_invalid_kernel) {
    alea_system_t* sys = alea_create();
    int sphere = alea_sphere_surface(sys, 1, 0, 0, 0, 1);
    ASSERT(sphere >= 0);
    ASSERT_EQ(alea_surface_set_boundary(sys, 1, ALEA_BOUNDARY_VACUUM), 0);
    int material = alea_add_material(sys, 7);
    ASSERT(material >= 0);
    ASSERT_EQ(alea_add_cell(sys, 1, alea_halfspace(sys, sphere, -1),
                            material, 1.0, 0), 0);
    double total[] = {0.0}, transfer[] = {0.0}, source[] = {1.0};
    alea_adjoint_material_t cells[] = {{total, transfer}};
    size_t detector_group = 0;
    alea_adjoint_problem_t p = {
        .particle = ALEA_NUC_PARTICLE_NEUTRON, .n_groups = 1,
        .cell_count = 1, .cell_materials = cells, .physical_source = source,
        .detector_sampler = center_detector,
        .detector_context = &detector_group
    };
    alea_adjoint_options_t o = {
        .histories = 1, .seed = 1, .max_events_per_history = 1,
        .max_segment_distance = 0.1
    };
    alea_adjoint_result_t result;
    alea_transport_failure_t failure;
    ASSERT_EQ(alea_adjoint_run(sys, &p, &o, &result, &failure),
              ALEA_ERR_OVERFLOW);
    ASSERT_EQ(result.histories, 0);
    ASSERT_EQ(failure.history_id, 0);
    transfer[0] = 1.0;
    ASSERT_EQ(alea_adjoint_run(sys, &p, &o, &result, NULL),
              ALEA_ERR_INVALID_ARG);
    alea_destroy(sys);
}

TEST(adjoint_box_detector_normalizes_and_rejects_zero_adjoint_rate) {
    alea_adjoint_box_detector_t box = {
        .lower = {-1, -2, -3}, .upper = {1, 2, 3},
        .group = 1, .response = 0.25
    };
    alea_adjoint_detector_particle_t sample;
    ASSERT_EQ(alea_adjoint_sample_box_detector(&box, 7, 0, &sample), ALEA_OK);
    ASSERT_EQ(sample.group, 1);
    ASSERT_NEAR(sample.weight, 12.0*4.0*acos(-1.0), 1e-12);
    for (int i = 0; i < 3; ++i) {
        ASSERT(sample.position[i] >= box.lower[i]);
        ASSERT(sample.position[i] <= box.upper[i]);
    }
    alea_adjoint_detector_particle_t repeat;
    ASSERT_EQ(alea_adjoint_sample_box_detector(&box, 7, 0, &repeat), ALEA_OK);
    ASSERT_EQ(sample.weight, repeat.weight);
    for (int i = 0; i < 3; ++i)
        ASSERT_EQ(sample.position[i], repeat.position[i]);

    alea_system_t* sys = alea_create();
    int sphere = alea_sphere_surface(sys, 1, 0, 0, 0, 1);
    ASSERT(sphere >= 0);
    int material = alea_add_material(sys, 8);
    ASSERT(material >= 0);
    ASSERT_EQ(alea_add_cell(sys, 1, alea_halfspace(sys, sphere, -1),
                            material, 1.0, 0), 0);
    double total[] = {1.0, 0.0};
    double transfer[] = {0.0, 0.1, 0.0, 0.0};
    double source[] = {1.0, 0.0};
    alea_adjoint_material_t cells[] = {{total, transfer}};
    size_t detector_group = 1;
    alea_adjoint_problem_t p = {
        .particle = ALEA_NUC_PARTICLE_NEUTRON, .n_groups = 2,
        .cell_count = 1, .cell_materials = cells, .physical_source = source,
        .detector_sampler = center_detector,
        .detector_context = &detector_group
    };
    alea_adjoint_options_t o = {
        .histories = 1, .max_events_per_history = 100,
        .max_segment_distance = 1
    };
    alea_adjoint_result_t result;
    ASSERT_EQ(alea_adjoint_run(sys, &p, &o, &result, NULL),
              ALEA_ERR_UNSUPPORTED);
    alea_destroy(sys);
}

TEST_MAIN()
