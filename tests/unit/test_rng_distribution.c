// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

#define ALEA_TEST_IMPLEMENTATION
#include "alea_test.h"
#include "rng/alea_rng_distribution.h"

#include <math.h>
#include <stdint.h>

static void init_event(alea_rng_event_t* event, uint32_t domain,
                       uint64_t entity, uint32_t event_index) {
    ASSERT_EQ(alea_rng_event_init(event, ALEA_RNG_PHILOX4X32_10, 42,
                                  domain, entity, event_index), 0);
}

TEST(event_address_matches_indexed_rng_and_separates_events) {
    alea_rng_event_t event;
    init_event(&event, ALEA_RNG_DOMAIN_TRANSPORT_FREE_PATH, 19, 7);
    uint32_t first, second, indexed, other_event;
    ASSERT_EQ(alea_rng_event_next_u32(&event, &first), 0);
    ASSERT_EQ(alea_rng_event_next_u32(&event, &second), 0);
    ASSERT_EQ(alea_rng_u32_at(
        ALEA_RNG_PHILOX4X32_10, 42, ALEA_RNG_DOMAIN_TRANSPORT_FREE_PATH,
        19, (UINT64_C(7) << 32), &indexed), 0);
    ASSERT_EQ(first, indexed);
    ASSERT_EQ(first, UINT32_C(0xa0b5a6df));
    ASSERT_EQ(second, UINT32_C(0x3ff7df4e));
    ASSERT_EQ(alea_rng_event_position(&event), 2);
    ASSERT_EQ(alea_rng_event_seek(&event, 0), 0);
    ASSERT_EQ(alea_rng_event_next_u32(&event, &indexed), 0);
    ASSERT_EQ(first, indexed);

    alea_rng_event_t next;
    init_event(&next, ALEA_RNG_DOMAIN_TRANSPORT_FREE_PATH, 19, 8);
    ASSERT_EQ(alea_rng_event_next_u32(&next, &other_event), 0);
    ASSERT_NE(first, other_event);
    ASSERT_NE(first, second);
}

TEST(transport_entity_identity_has_stable_word_layout) {
    ASSERT_EQ(alea_rng_transport_entity_id(UINT32_C(0x12345678),
                                           UINT32_C(0x9abcdef0)),
              UINT64_C(0x123456789abcdef0));
}

TEST(per_entity_samples_are_invariant_to_execution_order) {
    double forward[256];
    for (uint32_t history = 0; history < 256; history++) {
        alea_rng_event_t event;
        init_event(&event, ALEA_RNG_DOMAIN_TRANSPORT_FREE_PATH,
                   alea_rng_transport_entity_id(history, 0), 5);
        ASSERT_EQ(alea_rng_sample_exponential(&event, 0.75,
                                              &forward[history]), 0);
    }
    for (uint32_t history = 256; history-- > 0;) {
        alea_rng_event_t event;
        double reordered;
        init_event(&event, ALEA_RNG_DOMAIN_TRANSPORT_FREE_PATH,
                   alea_rng_transport_entity_id(history, 0), 5);
        ASSERT_EQ(alea_rng_sample_exponential(&event, 0.75, &reordered), 0);
        ASSERT(reordered == forward[history]);
    }
}

TEST(event_rejects_legacy_and_uniform_failure_is_transactional) {
    alea_rng_event_t event;
    ASSERT(alea_rng_event_init(&event, ALEA_RNG_LEGACY_LCG, 42,
                               ALEA_RNG_DOMAIN_TRANSPORT_REACTION, 1, 0) < 0);
    init_event(&event, ALEA_RNG_DOMAIN_TRANSPORT_REACTION, 1, UINT32_MAX);
    ASSERT_EQ(alea_rng_event_seek(&event, (UINT64_C(1) << 32) - 1), 0);
    double output = -1.0;
    ASSERT(alea_rng_event_next_uniform(&event, &output) < 0);
    ASSERT_EQ(alea_rng_event_position(&event), (UINT64_C(1) << 32) - 1);
}

TEST(continuous_samplers_have_fixed_consumption_and_valid_ranges) {
    alea_rng_event_t event;
    init_event(&event, ALEA_RNG_DOMAIN_TRANSPORT_FREE_PATH, 3, 11);
    double value;
    ASSERT_EQ(alea_rng_sample_exponential(&event, 2.5, &value), 0);
    ASSERT(isfinite(value));
    ASSERT(value > 0.0);
    ASSERT_EQ(alea_rng_event_position(&event), 2);
    ASSERT_EQ(alea_rng_sample_uniform_range(&event, -4.0, 7.0, &value), 0);
    ASSERT(value >= -4.0 && value < 7.0);
    ASSERT_EQ(alea_rng_event_position(&event), 4);

    const uint64_t before = alea_rng_event_position(&event);
    ASSERT(alea_rng_sample_exponential(&event, 0.0, &value) < 0);
    ASSERT_EQ(alea_rng_event_position(&event), before);
    ASSERT(alea_rng_sample_uniform_range(&event, 2.0, 2.0, &value) < 0);
    ASSERT_EQ(alea_rng_event_position(&event), before);
}

TEST(bounded_integer_is_unbiased_smoke_and_within_bound) {
    alea_rng_event_t event;
    init_event(&event, ALEA_RNG_DOMAIN_TRANSPORT_REACTION, 81, 0);
    size_t counts[7] = {0};
    const size_t samples = 70000;
    for (size_t i = 0; i < samples; i++) {
        uint32_t value;
        ASSERT_EQ(alea_rng_sample_u32_bounded(&event, 7, 0, &value), 0);
        ASSERT(value < 7);
        counts[value]++;
    }
    for (size_t i = 0; i < 7; i++)
        ASSERT_NEAR((double)counts[i], 10000.0, 350.0);
}

TEST(angular_samplers_have_expected_norms_and_moments) {
    alea_rng_event_t isotropic, cosine;
    init_event(&isotropic, ALEA_RNG_DOMAIN_TRANSPORT_SCATTER_ANGLE, 5, 0);
    init_event(&cosine, ALEA_RNG_DOMAIN_TRANSPORT_SOURCE_DIRECTION, 5, 0);
    double iso_z = 0.0, cosine_z = 0.0;
    const size_t samples = 30000;
    for (size_t i = 0; i < samples; i++) {
        double a[3], b[3];
        ASSERT_EQ(alea_rng_sample_isotropic_direction(&isotropic, a), 0);
        ASSERT_EQ(alea_rng_sample_cosine_direction(&cosine, b), 0);
        ASSERT_NEAR(a[0] * a[0] + a[1] * a[1] + a[2] * a[2], 1.0, 1e-12);
        ASSERT_NEAR(b[0] * b[0] + b[1] * b[1] + b[2] * b[2], 1.0, 1e-12);
        ASSERT(b[2] >= 0.0);
        iso_z += a[2];
        cosine_z += b[2];
    }
    ASSERT_NEAR(iso_z / (double)samples, 0.0, 0.015);
    ASSERT_NEAR(cosine_z / (double)samples, 2.0 / 3.0, 0.015);
    ASSERT_EQ(alea_rng_event_position(&isotropic), samples * 4);
    ASSERT_EQ(alea_rng_event_position(&cosine), samples * 4);
}

TEST(discrete_cdf_and_alias_tables_follow_weights) {
    const double weights[3] = {1.0, 3.0, 6.0};
    alea_rng_discrete_table_t* table =
        alea_rng_discrete_table_create(weights, 3);
    ASSERT_NOT_NULL(table);
    ASSERT_EQ(alea_rng_discrete_table_count(table), 3);

    alea_rng_event_t cdf, alias;
    init_event(&cdf, ALEA_RNG_DOMAIN_TRANSPORT_NUCLIDE, 91, 0);
    init_event(&alias, ALEA_RNG_DOMAIN_TRANSPORT_REACTION, 91, 0);
    size_t cdf_counts[3] = {0}, alias_counts[3] = {0};
    const size_t samples = 40000;
    for (size_t i = 0; i < samples; i++) {
        size_t index;
        ASSERT_EQ(alea_rng_sample_discrete_cdf(&cdf, table, &index), 0);
        ASSERT(index < 3);
        cdf_counts[index]++;
        ASSERT_EQ(alea_rng_sample_discrete_alias(&alias, table, &index), 0);
        ASSERT(index < 3);
        alias_counts[index]++;
    }
    for (size_t i = 0; i < 3; i++) {
        const double expected = samples * weights[i] / 10.0;
        ASSERT_NEAR((double)cdf_counts[i], expected, 350.0);
        ASSERT_NEAR((double)alias_counts[i], expected, 350.0);
    }
    alea_rng_discrete_table_destroy(table);
}

TEST(distribution_tables_reject_malformed_inputs) {
    const double negative[2] = {1.0, -1.0};
    const double zero[2] = {0.0, 0.0};
    ASSERT_NULL(alea_rng_discrete_table_create(negative, 2));
    ASSERT_NULL(alea_rng_discrete_table_create(zero, 2));

    const double bad_x[3] = {0.0, 2.0, 1.0};
    const double pdf[3] = {1.0, 1.0, 1.0};
    ASSERT_NULL(alea_rng_tabular_table_create(
        bad_x, pdf, 3, ALEA_RNG_TABULAR_LIN_LIN));
}

TEST(histogram_and_linlin_tabular_sampling_match_analytic_means) {
    const double histogram_x[3] = {0.0, 1.0, 3.0};
    const double histogram_pdf[3] = {1.0, 2.0, 0.0};
    const double triangle_x[2] = {0.0, 1.0};
    const double triangle_pdf[2] = {0.0, 2.0};
    alea_rng_tabular_table_t* histogram = alea_rng_tabular_table_create(
        histogram_x, histogram_pdf, 3, ALEA_RNG_TABULAR_HISTOGRAM);
    alea_rng_tabular_table_t* triangle = alea_rng_tabular_table_create(
        triangle_x, triangle_pdf, 2, ALEA_RNG_TABULAR_LIN_LIN);
    ASSERT_NOT_NULL(histogram);
    ASSERT_NOT_NULL(triangle);

    alea_rng_event_t histogram_event, triangle_event;
    init_event(&histogram_event,
               ALEA_RNG_DOMAIN_TRANSPORT_SECONDARY_ENERGY, 13, 0);
    init_event(&triangle_event,
               ALEA_RNG_DOMAIN_TRANSPORT_SOURCE_ENERGY, 13, 0);
    double histogram_sum = 0.0, triangle_sum = 0.0;
    const size_t samples = 50000;
    for (size_t i = 0; i < samples; i++) {
        double value;
        ASSERT_EQ(alea_rng_sample_tabular(&histogram_event, histogram, &value), 0);
        ASSERT(value >= 0.0 && value <= 3.0);
        histogram_sum += value;
        ASSERT_EQ(alea_rng_sample_tabular(&triangle_event, triangle, &value), 0);
        ASSERT(value >= 0.0 && value <= 1.0);
        triangle_sum += value;
    }
    ASSERT_NEAR(histogram_sum / (double)samples, 1.7, 0.015);
    ASSERT_NEAR(triangle_sum / (double)samples, 2.0 / 3.0, 0.01);
    ASSERT_EQ(alea_rng_event_position(&histogram_event), samples * 2);
    ASSERT_EQ(alea_rng_event_position(&triangle_event), samples * 2);

    alea_rng_tabular_table_destroy(histogram);
    alea_rng_tabular_table_destroy(triangle);
}

TEST_MAIN()
