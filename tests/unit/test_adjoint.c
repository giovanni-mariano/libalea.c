// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

#include "alea_test.h"
#include "alea_adjoint.h"
#include "rng/alea_rng.h"
#include <math.h>
#include <stdio.h>

static alea_error_t center_detector(void* context, uint64_t seed,
    uint32_t history, alea_adjoint_detector_particle_t* out) {
    (void)seed; (void)history;
    *out = (alea_adjoint_detector_particle_t){0};
    out->physical_direction[0] = 1.0;
    out->group = *(const size_t*)context;
    out->weight = 4.0 * acos(-1.0);
    return ALEA_OK;
}

static alea_error_t offset_detector(void* context, uint64_t seed,
    uint32_t history, alea_adjoint_detector_particle_t* out) {
    alea_error_t err = center_detector(context, seed, history, out);
    out->position[0] = 0.5;
    return err;
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
        alea_adjoint_failure_t failure;
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
    alea_adjoint_failure_t failure;
    ASSERT_EQ(alea_adjoint_run(sys, &p, &o, &result, &failure),
              ALEA_ERR_OVERFLOW);
    ASSERT_EQ(result.histories, 0);
    ASSERT_EQ(failure.history_id, 0);
    ASSERT_EQ(failure.group, 0);
    transfer[0] = 1.0;
    ASSERT_EQ(alea_adjoint_run(sys, &p, &o, &result, NULL),
              ALEA_ERR_INVALID_ARG);
    alea_destroy(sys);
}

TEST(neutron_group_material_combines_nuclides_and_rejects_fission) {
    double bounds[] = {2.0, 1.0, 0.0};
    alea_nuc_multigroup_t* a = alea_nuc_mg_create(2, bounds);
    alea_nuc_multigroup_t* b = alea_nuc_mg_create(2, bounds);
    ASSERT_NOT_NULL(a);
    ASSERT_NOT_NULL(b);
    a->sigma_t[0] = 2.0; a->sigma_t[1] = 3.0;
    b->sigma_t[0] = 4.0; b->sigma_t[1] = 5.0;
    a->scatter[0] = 0.5; a->scatter[1] = 1.0;
    b->scatter[0] = 1.5; b->scatter[1] = 2.0;
    const alea_nuc_multigroup_t* nuclides[] = {a, b};
    double densities[] = {0.1, 0.2};
    double total[] = {-1.0, -1.0};
    double transfer[] = {-1.0, -1.0, -1.0, -1.0};
    ASSERT_EQ(alea_mg_neutron_material_build(nuclides, densities, 2, 2,
              total, transfer), ALEA_OK);
    ASSERT_NEAR(total[0], 1.0, 1e-14);
    ASSERT_NEAR(total[1], 1.3, 1e-14);
    ASSERT_NEAR(transfer[0], 0.35, 1e-14);
    ASSERT_NEAR(transfer[1], 0.5, 1e-14);
    b->sigma_f[0] = 0.1;
    total[0] = -1.0;
    ASSERT_EQ(alea_mg_neutron_material_build(nuclides, densities, 2, 2,
              total, transfer), ALEA_ERR_UNSUPPORTED);
    ASSERT_EQ(total[0], -1.0);
    b->sigma_f[0] = 0.0;
    b->bounds[1] = 0.5;
    ASSERT_EQ(alea_mg_neutron_material_build(nuclides, densities, 2, 2,
              total, transfer), ALEA_ERR_INVALID_ARG);
    ASSERT_EQ(total[0], -1.0);
    alea_nuc_mg_destroy(a);
    alea_nuc_mg_destroy(b);
}

TEST(single_adjoint_history_has_unknown_standard_error) {
    alea_system_t* sys = alea_create();
    int sphere = alea_sphere_surface(sys, 1, 0, 0, 0, 1);
    ASSERT(sphere >= 0);
    ASSERT_EQ(alea_surface_set_boundary(sys, 1, ALEA_BOUNDARY_VACUUM), 0);
    int material = alea_add_material(sys, 12);
    ASSERT(material >= 0);
    ASSERT_EQ(alea_add_cell(sys, 1, alea_halfspace(sys, sphere, -1),
                            material, 1.0, 0), 0);
    double total[] = {0.0}, transfer[] = {0.0}, source[] = {1.0};
    alea_mg_material_t cells[] = {{total, transfer}};
    size_t detector_group = 0;
    alea_adjoint_problem_t p = {
        .particle = ALEA_NUC_PARTICLE_NEUTRON, .n_groups = 1,
        .cell_count = 1, .cell_materials = cells, .physical_source = source,
        .detector_sampler = center_detector,
        .detector_context = &detector_group
    };
    alea_adjoint_options_t o = {
        .histories = 1, .seed = 1, .max_events_per_history = 100,
        .max_segment_distance = 0.25
    };
    alea_adjoint_result_t result;
    ASSERT_EQ(alea_adjoint_run(sys, &p, &o, &result, NULL), ALEA_OK);
    ASSERT_NEAR(result.mean, 1.0, 1e-10);
    ASSERT(isnan(result.standard_error));
    alea_destroy(sys);
}

/* Independent forward reference: sample an isotropic source in the shell and
 * integrate attenuated track length through the detector box analytically. */
TEST(adjoint_matches_forward_source_sampling_in_absorbing_sphere) {
    alea_system_t* sys = alea_create();
    int outer = alea_sphere_surface(sys, 2, 0, 0, 0, 2);
    ASSERT(outer >= 0);
    ASSERT_EQ(alea_surface_set_boundary(sys, 2, ALEA_BOUNDARY_VACUUM), 0);
    int material = alea_add_material(sys, 16);
    ASSERT(material >= 0);
    alea_node_id_t inside = alea_halfspace(sys, outer, -1);
    ASSERT_EQ(alea_add_cell(sys, 1, inside, material, 1.0, 0), 0);
    double total[] = {0.7}, transfer[] = {0.0}, source[] = {1.0};
    alea_mg_material_t cells[] = {{total, transfer}};
    alea_adjoint_box_detector_t detector = {
        .lower = {-0.8, -0.8, -0.8}, .upper = {0.8, 0.8, 0.8},
        .group = 0, .response = 0.5
    };
    alea_adjoint_problem_t p = {
        .particle = ALEA_NUC_PARTICLE_NEUTRON, .n_groups = 1,
        .cell_count = 1, .cell_materials = cells, .physical_source = source,
        .detector_sampler = alea_adjoint_sample_box_detector,
        .detector_context = &detector
    };
    const uint32_t histories = 150000;
    alea_adjoint_options_t o = {
        .histories = histories, .seed = 101, .max_events_per_history = 80,
        .max_segment_distance = 3.0
    };
    alea_adjoint_result_t adj;
    alea_adjoint_failure_t failure;
    alea_error_t adj_err = alea_adjoint_run(sys, &p, &o, &adj, &failure);
    if (adj_err != ALEA_OK)
        fprintf(stderr, "adjoint error %d history %u event %u cell %d group %zu: %s\n",
            adj_err, failure.history_id, failure.event_index,
            failure.cell_id, failure.group, alea_get_error_detail());
    ASSERT_EQ(adj_err, ALEA_OK);

    double forward_mean = 0.0, forward_m2 = 0.0;
    const double source_volume = (4.0/3.0)*acos(-1.0)*8.0;
    for (uint32_t h = 0; h < histories; ++h) {
        double u[5];
        for (uint32_t j = 0; j < 5; ++j)
            ASSERT_EQ(alea_rng_uniform53_at(ALEA_RNG_PHILOX4X32_10,
                202, ALEA_RNG_DOMAIN_TEST_ONLY, h, j, &u[j]), 0);
        double r = 2.0*cbrt(u[0]);
        double z0 = 2.0*u[1]-1.0, z1 = 2.0*u[3]-1.0;
        double a0 = 2.0*acos(-1.0)*u[2];
        double a1 = 2.0*acos(-1.0)*u[4];
        double position[3] = {r*sqrt(1.0-z0*z0)*cos(a0),
                              r*sqrt(1.0-z0*z0)*sin(a0), r*z0};
        double direction[3] = {sqrt(1.0-z1*z1)*cos(a1),
                               sqrt(1.0-z1*z1)*sin(a1), z1};
        double enter = 0.0, leave = INFINITY;
        for (int j = 0; j < 3; ++j) {
            if (fabs(direction[j]) < 1e-15) {
                if (fabs(position[j]) > 0.8) leave = -1.0;
                continue;
            }
            double a = (-0.8-position[j])/direction[j];
            double b = ( 0.8-position[j])/direction[j];
            if (a > b) { double tmp = a; a = b; b = tmp; }
            if (a > enter) enter = a;
            if (b < leave) leave = b;
        }
        double score = leave > enter ? source_volume*detector.response*
            (exp(-0.7*enter)-exp(-0.7*leave))/0.7 : 0.0;
        double delta = score-forward_mean;
        forward_mean += delta/(h+1.0);
        forward_m2 += delta*(score-forward_mean);
    }
    double forward_se = sqrt(forward_m2/(histories*(double)(histories-1)));
    ASSERT(fabs(adj.mean-forward_mean) <
           5.0*hypot(adj.standard_error, forward_se));
    alea_destroy(sys);
}

TEST(adjoint_preserves_optical_depth_across_adjacent_cells) {
    alea_system_t* sys = alea_create();
    int sphere = alea_sphere_surface(sys, 1, 0, 0, 0, 2);
    int plane = alea_plane_surface(sys, 2, 1, 0, 0, 0);
    ASSERT(sphere >= 0 && plane >= 0);
    ASSERT_EQ(alea_surface_set_boundary(sys, 1, ALEA_BOUNDARY_VACUUM), 0);
    int material = alea_add_material(sys, 19);
    ASSERT(material >= 0);
    alea_node_id_t inside = alea_halfspace(sys, sphere, -1);
    alea_node_id_t left = alea_intersection(sys, inside,
        alea_halfspace(sys, plane, -1));
    alea_node_id_t right = alea_intersection(sys, inside,
        alea_halfspace(sys, plane, +1));
    ASSERT_EQ(alea_add_cell(sys, 1, left, material, 1.0, 0), 0);
    ASSERT_EQ(alea_add_cell(sys, 2, right, material, 1.0, 0), 1);
    double total[] = {1.0}, transfer[] = {0.0};
    double source[] = {1.0, 0.0};
    alea_mg_material_t cells[] = {{total, transfer}, {total, transfer}};
    size_t detector_group = 0;
    alea_adjoint_problem_t p = {
        .particle = ALEA_NUC_PARTICLE_NEUTRON, .n_groups = 1,
        .cell_count = 2, .cell_materials = cells, .physical_source = source,
        .detector_sampler = offset_detector,
        .detector_context = &detector_group
    };
    alea_adjoint_options_t o = {
        .histories = 60000, .seed = 73, .max_events_per_history = 50,
        .max_segment_distance = 3.0
    };
    alea_adjoint_result_t result;
    alea_adjoint_failure_t failure;
    alea_error_t err = alea_adjoint_run(sys, &p, &o, &result, &failure);
    if (err != ALEA_OK)
        fprintf(stderr, "adjacent-cell adjoint error %d history %u event %u: %s\n",
            err, failure.history_id, failure.event_index,
            alea_get_error_detail());
    ASSERT_EQ(err, ALEA_OK);
    ASSERT_NEAR(result.mean, exp(-0.5)*(1.0-exp(-2.0)), 0.012);
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
