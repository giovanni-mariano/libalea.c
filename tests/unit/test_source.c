// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

#include "alea_test.h"
#include "alea_source.h"
#include <math.h>
#include <string.h>

TEST(point_and_box_components_replay_and_normalize) {
    alea_source_spec_t spec = {
        .particle = ALEA_NUC_PARTICLE_NEUTRON,
        .space = ALEA_SOURCE_POINT,
        .angle = ALEA_SOURCE_MONODIRECTIONAL,
        .position = {1, 2, 3}, .direction = {2, 0, 0},
        .energy = 14.1, .weight = 1, .time = 0.5
    };
    alea_source_t* source = NULL;
    ASSERT_EQ(alea_source_prepare(&spec, &source), ALEA_OK);
    ASSERT_NOT_NULL(source);
    ASSERT_EQ(alea_source_particle_mask(source), ALEA_NUC_BIND_NEUTRON);
    alea_transport_source_t sample;
    ASSERT_EQ(alea_source_sample(source, 11, 7, &sample), ALEA_OK);
    ASSERT_NEAR(sample.position[0], 1, 1e-14);
    ASSERT_NEAR(sample.particle.direction[0], 1, 1e-14);
    ASSERT_NEAR(sample.particle.direction[1], 0, 1e-14);
    ASSERT_NEAR(sample.particle.energy, 14.1, 1e-14);
    ASSERT_NEAR(sample.particle.time, 0.5, 1e-14);
    alea_source_free(source);

    spec.particle = ALEA_NUC_PARTICLE_PHOTON;
    spec.space = ALEA_SOURCE_BOX;
    spec.angle = ALEA_SOURCE_ISOTROPIC;
    spec.lower[0] = -2; spec.lower[1] = -1; spec.lower[2] = 0;
    spec.upper[0] = 3; spec.upper[1] = 4; spec.upper[2] = 5;
    ASSERT_EQ(alea_source_prepare(&spec, &source), ALEA_OK);
    ASSERT_EQ(alea_source_particle_mask(source), ALEA_NUC_BIND_PHOTON);
    alea_transport_source_t batch[5];
    ASSERT_EQ(alea_source_sample_batch(source, 99, 12, 5, batch), ALEA_OK);
    for (uint32_t i = 0; i < 5; ++i) {
        alea_transport_source_t replay;
        ASSERT_EQ(alea_source_sample(source, 99, 12 + i, &replay), ALEA_OK);
        ASSERT_EQ(memcmp(&batch[i], &replay, sizeof(replay)), 0);
        for (int axis = 0; axis < 3; ++axis) {
            ASSERT(batch[i].position[axis] >= spec.lower[axis]);
            ASSERT(batch[i].position[axis] < spec.upper[axis]);
        }
        double* d = batch[i].particle.direction;
        ASSERT_NEAR(d[0]*d[0] + d[1]*d[1] + d[2]*d[2], 1, 1e-12);
    }
    ASSERT_EQ(alea_source_sample_batch(source, 99, UINT32_MAX, 2, batch), ALEA_ERR_OVERFLOW);
    alea_source_free(source);
}

TEST(source_validation_rejects_bad_inputs) {
    alea_source_spec_t spec = {
        .space = ALEA_SOURCE_POINT, .angle = ALEA_SOURCE_MONODIRECTIONAL,
        .direction = {1, 0, 0}, .energy = 2, .weight = 1
    };
    alea_source_t* source = NULL;
    spec.energy = 0;
    ASSERT_EQ(alea_source_prepare(&spec, &source), ALEA_ERR_INVALID_ARG);
    ASSERT_EQ(source, NULL);
    spec.energy = 2;
    spec.direction[0] = 0;
    ASSERT_EQ(alea_source_prepare(&spec, &source), ALEA_ERR_INVALID_ARG);
    spec.angle = ALEA_SOURCE_ISOTROPIC;
    spec.space = ALEA_SOURCE_BOX;
    spec.lower[0] = 2; spec.upper[0] = 1;
    ASSERT_EQ(alea_source_prepare(&spec, &source), ALEA_ERR_INVALID_ARG);
}

TEST(analytic_spatial_sources_use_length_area_and_volume_measures) {
    alea_source_spec_t spec = {
        .particle = ALEA_NUC_PARTICLE_NEUTRON,
        .space = ALEA_SOURCE_LINE, .angle = ALEA_SOURCE_ISOTROPIC,
        .start = {0, 0, 0}, .end = {10, 0, 0},
        .energy = 14.1, .weight = 1
    };
    alea_source_t* source = NULL;
    ASSERT_EQ(alea_source_prepare(&spec, &source), ALEA_OK);
    double mean = 0;
    for (uint32_t i = 0; i < 10000; ++i) {
        alea_transport_source_t sample;
        ASSERT_EQ(alea_source_sample(source, 1, i, &sample), ALEA_OK);
        ASSERT(sample.position[0] >= 0 && sample.position[0] <= 10);
        ASSERT_EQ(sample.position[1], 0);
        mean += sample.position[0];
    }
    ASSERT_NEAR(mean / 10000, 5, 0.1);
    alea_source_free(source);

    spec.space = ALEA_SOURCE_SPHERE;
    spec.inner_radius = 1;
    spec.outer_radius = 3;
    ASSERT_EQ(alea_source_prepare(&spec, &source), ALEA_OK);
    mean = 0;
    for (uint32_t i = 0; i < 10000; ++i) {
        alea_transport_source_t sample;
        ASSERT_EQ(alea_source_sample(source, 2, i, &sample), ALEA_OK);
        double* p = sample.position;
        double r = sqrt(p[0]*p[0] + p[1]*p[1] + p[2]*p[2]);
        ASSERT(r >= 1 - 1e-12 && r <= 3 + 1e-12);
        mean += r*r*r;
    }
    ASSERT_NEAR(mean / 10000, 14, 0.3);
    alea_source_free(source);

    spec.space = ALEA_SOURCE_CYLINDER;
    spec.axis[0] = 4; spec.axis[1] = 0; spec.axis[2] = 0;
    spec.inner_radius = 1;
    spec.outer_radius = 3;
    ASSERT_EQ(alea_source_prepare(&spec, &source), ALEA_OK);
    double mean_r2 = 0, mean_x = 0;
    for (uint32_t i = 0; i < 10000; ++i) {
        alea_transport_source_t sample;
        ASSERT_EQ(alea_source_sample(source, 3, i, &sample), ALEA_OK);
        double* p = sample.position;
        double r2 = p[1]*p[1] + p[2]*p[2];
        ASSERT(p[0] >= 0 && p[0] <= 4);
        ASSERT(r2 >= 1 - 1e-12 && r2 <= 9 + 1e-12);
        mean_r2 += r2; mean_x += p[0];
    }
    ASSERT_NEAR(mean_r2 / 10000, 5, 0.1);
    ASSERT_NEAR(mean_x / 10000, 2, 0.1);
    alea_source_free(source);
}

TEST(discrete_energy_lines_copy_input_and_replay) {
    double energies[] = {2.45, 14.1};
    double weights[] = {1, 3};
    alea_source_spec_t spec = {
        .space = ALEA_SOURCE_POINT, .angle = ALEA_SOURCE_ISOTROPIC,
        .energy_type = ALEA_SOURCE_ENERGY_LINES,
        .energy_values = energies, .energy_weights = weights,
        .energy_count = 2, .weight = 1
    };
    alea_source_t* source = NULL;
    ASSERT_EQ(alea_source_prepare(&spec, &source), ALEA_OK);
    energies[0] = -100;
    weights[0] = 0;
    size_t high = 0;
    for (uint32_t i = 0; i < 10000; ++i) {
        alea_transport_source_t sample, replay;
        ASSERT_EQ(alea_source_sample(source, 37, i, &sample), ALEA_OK);
        ASSERT_EQ(alea_source_sample(source, 37, i, &replay), ALEA_OK);
        ASSERT_EQ(sample.particle.energy, replay.particle.energy);
        ASSERT(sample.particle.energy == 2.45 || sample.particle.energy == 14.1);
        high += sample.particle.energy == 14.1;
    }
    ASSERT_NEAR((double)high / 10000, 0.75, 0.02);
    alea_source_free(source);
    spec.energy_values = energies;
    ASSERT_EQ(alea_source_prepare(&spec, &source), ALEA_ERR_INVALID_ARG);
}

TEST(angular_laws_follow_solid_angle_and_polar_pdf) {
    alea_source_spec_t spec = {
        .space = ALEA_SOURCE_POINT, .position = {1, 0, 0},
        .angle = ALEA_SOURCE_CONE, .direction = {0, 0, 2},
        .cone_half_angle = 1.04719755119659774615,
        .energy = 14.1, .weight = 1
    };
    alea_source_t* source = NULL;
    ASSERT_EQ(alea_source_prepare(&spec, &source), ALEA_OK);
    double mean = 0;
    for (uint32_t i = 0; i < 10000; ++i) {
        alea_transport_source_t sample;
        ASSERT_EQ(alea_source_sample(source, 9, i, &sample), ALEA_OK);
        ASSERT(sample.particle.direction[2] >= 0.5 - 1e-14);
        mean += sample.particle.direction[2];
    }
    ASSERT_NEAR(mean / 10000, 0.75, 0.015);
    alea_source_free(source);

    spec.angle = ALEA_SOURCE_COSINE;
    ASSERT_EQ(alea_source_prepare(&spec, &source), ALEA_OK);
    mean = 0;
    for (uint32_t i = 0; i < 10000; ++i) {
        alea_transport_source_t sample;
        ASSERT_EQ(alea_source_sample(source, 10, i, &sample), ALEA_OK);
        ASSERT(sample.particle.direction[2] >= 0);
        mean += sample.particle.direction[2];
    }
    ASSERT_NEAR(mean / 10000, 2.0/3.0, 0.015);
    alea_source_free(source);

    double mu[] = {0, 1};
    double pdf[] = {0, 2};
    spec.angle = ALEA_SOURCE_TABULATED_MU;
    spec.angle_mu = mu; spec.angle_pdf = pdf; spec.angle_count = 2;
    spec.angle_interpolation = ALEA_SOURCE_PDF_LINEAR;
    ASSERT_EQ(alea_source_prepare(&spec, &source), ALEA_OK);
    mu[0] = -2; pdf[1] = -1;
    mean = 0;
    for (uint32_t i = 0; i < 10000; ++i) {
        alea_transport_source_t sample, replay;
        ASSERT_EQ(alea_source_sample(source, 11, i, &sample), ALEA_OK);
        ASSERT_EQ(alea_source_sample(source, 11, i, &replay), ALEA_OK);
        ASSERT_EQ(memcmp(&sample, &replay, sizeof(sample)), 0);
        ASSERT(sample.particle.direction[2] >= 0);
        mean += sample.particle.direction[2];
    }
    ASSERT_NEAR(mean / 10000, 2.0/3.0, 0.015);
    alea_source_free(source);
    ASSERT_EQ(alea_source_prepare(&spec, &source), ALEA_ERR_INVALID_ARG);

    mu[0] = 0; pdf[1] = 2;
    spec.angle_interpolation = ALEA_SOURCE_PDF_HISTOGRAM;
    pdf[0] = 1;
    ASSERT_EQ(alea_source_prepare(&spec, &source), ALEA_OK);
    mean = 0;
    for (uint32_t i = 0; i < 10000; ++i) {
        alea_transport_source_t draw;
        ASSERT_EQ(alea_source_sample(source, 13, i, &draw), ALEA_OK);
        mean += draw.particle.direction[2];
    }
    ASSERT_NEAR(mean / 10000, 0.5, 0.015);
    alea_source_free(source);

    spec.angle = ALEA_SOURCE_CONE;
    spec.cone_half_angle = 4;
    ASSERT_EQ(alea_source_prepare(&spec, &source), ALEA_ERR_INVALID_ARG);
    spec.cone_half_angle = 0.2;

    spec.angle = ALEA_SOURCE_RADIAL;
    spec.space = ALEA_SOURCE_LINE;
    spec.start[0] = 1; spec.end[0] = 2;
    spec.radial_inward = 1;
    ASSERT_EQ(alea_source_prepare(&spec, &source), ALEA_OK);
    alea_transport_source_t sample;
    ASSERT_EQ(alea_source_sample(source, 12, 0, &sample), ALEA_OK);
    ASSERT_NEAR(sample.particle.direction[0], -1, 1e-14);
    alea_source_free(source);
}

TEST(tokamak_rz_emissivity_uses_cylindrical_volume) {
    double r_edges[] = {0, 1, 2};
    double z_edges[] = {0, 1, 2};
    double emissivity[] = {1, 0, 1, 1}; /* R-major; inner, upper bin is dark */
    alea_source_spec_t spec = {
        .space = ALEA_SOURCE_TOKAMAK_RZ,
        .r_edges = r_edges, .r_edge_count = 3,
        .z_edges = z_edges, .z_edge_count = 3,
        .rz_emissivity = emissivity,
        .phi_min = 0, .phi_max = 1.57079632679489661923,
        .angle = ALEA_SOURCE_ISOTROPIC,
        .energy = 14.1, .weight = 1
    };
    alea_source_t* source = NULL;
    ASSERT_EQ(alea_source_prepare(&spec, &source), ALEA_OK);
    double integrated;
    ASSERT_EQ(alea_source_integrated_emissivity(source, &integrated), ALEA_OK);
    ASSERT_NEAR(integrated, 1.75*3.14159265358979323846, 1e-12);
    r_edges[1] = 100;
    z_edges[1] = 100;
    emissivity[0] = 0;
    int outer = 0, upper = 0;
    double mean_r2 = 0;
    for (uint32_t i = 0; i < 20000; ++i) {
        alea_transport_source_t sample;
        ASSERT_EQ(alea_source_sample(source, 73, i, &sample), ALEA_OK);
        double x = sample.position[0], y = sample.position[1];
        double r2 = x*x + y*y;
        ASSERT(x >= -1e-14 && y >= -1e-14);
        ASSERT(r2 >= 0 && r2 < 4 + 1e-12);
        ASSERT(sample.position[2] >= 0 && sample.position[2] < 2);
        ASSERT(!(r2 < 1 && sample.position[2] >= 1));
        outer += r2 >= 1;
        upper += sample.position[2] >= 1;
        mean_r2 += r2;
    }
    ASSERT_NEAR((double)outer / 20000, 6.0/7.0, 0.015);
    ASSERT_NEAR((double)upper / 20000, 3.0/7.0, 0.015);
    ASSERT_NEAR(mean_r2 / 20000, 15.5/7.0, 0.05);
    alea_transport_source_t first, replay;
    ASSERT_EQ(alea_source_sample(source, 73, 123, &first), ALEA_OK);
    ASSERT_EQ(alea_source_sample(source, 73, 123, &replay), ALEA_OK);
    ASSERT_EQ(memcmp(&first, &replay, sizeof(first)), 0);
    alea_source_free(source);

    spec.r_edges = r_edges; r_edges[1] = 1;
    spec.z_edges = z_edges; z_edges[1] = 1;
    emissivity[0] = 0;
    ASSERT_EQ(alea_source_prepare(&spec, &source), ALEA_OK);
    alea_source_free(source);
    spec.rz_emissivity = (double[]){0, 0, 0, 0};
    ASSERT_EQ(alea_source_prepare(&spec, &source), ALEA_ERR_INVALID_ARG);
    spec.rz_emissivity = emissivity;
    r_edges[0] = -1;
    ASSERT_EQ(alea_source_prepare(&spec, &source), ALEA_ERR_INVALID_ARG);
    r_edges[0] = 0;
    spec.phi_min = 0;
    spec.phi_max = 0;
    ASSERT_EQ(alea_source_prepare(&spec, &source), ALEA_OK);
    ASSERT_EQ(alea_source_integrated_emissivity(source, &integrated), ALEA_OK);
    ASSERT_NEAR(integrated, 6.0*3.14159265358979323846, 1e-12);
    alea_source_free(source);
}

TEST(tabulated_energy_samples_pdf_in_mev) {
    double values[] = {1, 3};
    double pdf[] = {1, 1};
    alea_source_spec_t spec = {
        .space = ALEA_SOURCE_POINT, .angle = ALEA_SOURCE_ISOTROPIC,
        .energy_type = ALEA_SOURCE_ENERGY_TABULATED,
        .energy_values = values, .energy_weights = pdf,
        .energy_count = 2, .energy_interpolation = ALEA_SOURCE_PDF_LINEAR,
        .weight = 1
    };
    alea_source_t* source = NULL;
    ASSERT_EQ(alea_source_prepare(&spec, &source), ALEA_OK);
    values[0] = -1; pdf[0] = -1;
    double mean = 0;
    for (uint32_t i = 0; i < 10000; ++i) {
        alea_transport_source_t sample;
        ASSERT_EQ(alea_source_sample(source, 42, i, &sample), ALEA_OK);
        ASSERT(sample.particle.energy >= 1 && sample.particle.energy <= 3);
        mean += sample.particle.energy;
    }
    ASSERT_NEAR(mean / 10000, 2, 0.03);
    alea_source_free(source);
    ASSERT_EQ(alea_source_prepare(&spec, &source), ALEA_ERR_INVALID_ARG);
}

TEST(source_mixture_selects_physical_strength_without_reweighting) {
    alea_source_spec_t spec = {
        .space = ALEA_SOURCE_POINT, .angle = ALEA_SOURCE_ISOTROPIC,
        .energy = 2, .weight = 1
    };
    alea_source_t* parts[2] = {NULL, NULL};
    ASSERT_EQ(alea_source_prepare(&spec, &parts[0]), ALEA_OK);
    spec.position[0] = 1;
    spec.particle = ALEA_NUC_PARTICLE_PHOTON;
    spec.energy = 3;
    ASSERT_EQ(alea_source_prepare(&spec, &parts[1]), ALEA_OK);
    double strengths[] = {1, 3};
    alea_source_t* mixture = NULL;
    ASSERT_EQ(alea_source_mixture_prepare(parts, strengths, 2, &mixture), ALEA_OK);
    ASSERT_EQ(alea_source_particle_mask(mixture),
              ALEA_NUC_BIND_NEUTRON | ALEA_NUC_BIND_PHOTON);
    int photons = 0;
    for (uint32_t i = 0; i < 10000; ++i) {
        alea_transport_source_t sample, replay;
        ASSERT_EQ(alea_source_sample(mixture, 17, i, &sample), ALEA_OK);
        ASSERT_EQ(alea_source_sample(mixture, 17, i, &replay), ALEA_OK);
        ASSERT_EQ(memcmp(&sample, &replay, sizeof(sample)), 0);
        ASSERT_EQ(sample.particle.weight, 1);
        if (sample.particle.type == ALEA_NUC_PARTICLE_PHOTON) {
            ASSERT_EQ(sample.position[0], 1);
            ASSERT_EQ(sample.particle.energy, 3);
            photons++;
        } else {
            ASSERT_EQ(sample.position[0], 0);
            ASSERT_EQ(sample.particle.energy, 2);
        }
    }
    ASSERT_NEAR((double)photons / 10000, 0.75, 0.02);
    alea_source_free(mixture);
    ASSERT_EQ(alea_source_mixture_prepare(NULL, strengths, 2, &mixture),
              ALEA_ERR_NULL_ARG);
    spec.particle = ALEA_NUC_PARTICLE_NEUTRON;
    spec.position[0] = 0;
    ASSERT_EQ(alea_source_prepare(&spec, &parts[0]), ALEA_OK);
    spec.particle = ALEA_NUC_PARTICLE_PHOTON;
    ASSERT_EQ(alea_source_prepare(&spec, &parts[1]), ALEA_OK);
    strengths[1] = 0;
    ASSERT_EQ(alea_source_mixture_prepare(parts, strengths, 2, &mixture), ALEA_OK);
    ASSERT_EQ(alea_source_particle_mask(mixture), ALEA_NUC_BIND_NEUTRON);
    for (uint32_t i = 0; i < 1000; ++i) {
        alea_transport_source_t sample;
        ASSERT_EQ(alea_source_sample(mixture, 19, i, &sample), ALEA_OK);
        ASSERT_EQ(sample.particle.type, ALEA_NUC_PARTICLE_NEUTRON);
    }
    alea_source_free(mixture);
}

TEST(nested_mixtures_use_independent_component_draws) {
    alea_source_spec_t spec = {
        .space = ALEA_SOURCE_POINT, .angle = ALEA_SOURCE_ISOTROPIC,
        .energy = 1, .weight = 1
    };
    alea_source_t* parts[3] = {NULL, NULL, NULL};
    for (int i = 0; i < 3; ++i) {
        spec.position[0] = i;
        ASSERT_EQ(alea_source_prepare(&spec, &parts[i]), ALEA_OK);
    }
    alea_source_t* inner = NULL;
    double equal[] = {1, 1};
    ASSERT_EQ(alea_source_mixture_prepare(&parts[1], equal, 2, &inner), ALEA_OK);
    alea_source_t* outer_parts[] = {parts[0], inner};
    alea_source_t* outer = NULL;
    ASSERT_EQ(alea_source_mixture_prepare(outer_parts, equal, 2, &outer), ALEA_OK);
    int counts[3] = {0};
    for (uint32_t i = 0; i < 20000; ++i) {
        alea_transport_source_t sample;
        ASSERT_EQ(alea_source_sample(outer, 71, i, &sample), ALEA_OK);
        int category = (int)sample.position[0];
        ASSERT(category >= 0 && category < 3);
        counts[category]++;
    }
    ASSERT_NEAR((double)counts[0] / 20000, 0.5, 0.02);
    ASSERT_NEAR((double)counts[1] / 20000, 0.25, 0.02);
    ASSERT_NEAR((double)counts[2] / 20000, 0.25, 0.02);
    alea_source_free(outer);
}

TEST_MAIN()
