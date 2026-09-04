// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

#define ALEA_TEST_IMPLEMENTATION
#include "alea_test.h"
#include "rng/alea_rng.h"

#include <stdint.h>
#include <string.h>

TEST(philox4x32_10_zero_known_answer) {
    const uint32_t counter[4] = {0, 0, 0, 0};
    const uint32_t key[2] = {0, 0};
    const uint32_t expected[4] = {
        UINT32_C(0x6627e8d5), UINT32_C(0xe169c58d),
        UINT32_C(0xbc57ac4c), UINT32_C(0x9b00dbd8)
    };
    uint32_t actual[4];
    alea_philox4x32_10(counter, key, actual);
    ASSERT_EQ(memcmp(actual, expected, sizeof(expected)), 0);
}

TEST(philox4x32_10_random123_known_answer) {
    /* Random123's kat_vectors tuple is counter words followed by key words. */
    const uint32_t counter[4] = {
        UINT32_C(0x243f6a88), UINT32_C(0x85a308d3),
        UINT32_C(0x13198a2e), UINT32_C(0x03707344)
    };
    const uint32_t key[2] = {
        UINT32_C(0xa4093822), UINT32_C(0x299f31d0)
    };
    const uint32_t expected[4] = {
        UINT32_C(0xd16cfe09), UINT32_C(0x94fdcceb),
        UINT32_C(0x5001e420), UINT32_C(0x24126ea1)
    };
    uint32_t actual[4];
    alea_philox4x32_10(counter, key, actual);
    ASSERT_EQ(memcmp(actual, expected, sizeof(expected)), 0);
}

TEST(rng_indexed_and_stream_access_match) {
    alea_rng_stream_t stream;
    alea_rng_stream_init(&stream, ALEA_RNG_PHILOX4X32_10,
                         UINT64_C(123456789),
                         ALEA_RNG_DOMAIN_TEST_ONLY, UINT64_C(987654321));
    for (uint64_t draw = 0; draw < 20; draw++) {
        uint32_t indexed, sequential;
        ASSERT_EQ(alea_rng_u32_at(
            ALEA_RNG_PHILOX4X32_10, UINT64_C(123456789),
            ALEA_RNG_DOMAIN_TEST_ONLY, UINT64_C(987654321), draw,
            &indexed), 0);
        ASSERT_EQ(alea_rng_stream_next_u32(&stream, &sequential), 0);
        ASSERT_EQ(indexed, sequential);
    }
    ASSERT_EQ(alea_rng_stream_position(&stream), 20);
}

TEST(rng_block_and_u64_access_match_words) {
    uint32_t block[4];
    ASSERT_EQ(alea_rng_block_at(
        ALEA_RNG_PHILOX4X32_10, 99, ALEA_RNG_DOMAIN_TEST_ONLY, 123, 2,
        block), 0);
    for (uint64_t lane = 0; lane < 4; lane++) {
        uint32_t indexed;
        ASSERT_EQ(alea_rng_u32_at(
            ALEA_RNG_PHILOX4X32_10, 99, ALEA_RNG_DOMAIN_TEST_ONLY, 123,
            8 + lane, &indexed), 0);
        ASSERT_EQ(block[lane], indexed);
    }
    uint64_t value;
    ASSERT_EQ(alea_rng_u64_at(
        ALEA_RNG_PHILOX4X32_10, 99, ALEA_RNG_DOMAIN_TEST_ONLY, 123, 4,
        &value), 0);
    ASSERT_EQ(value, ((uint64_t)block[0] << 32) | block[1]);
}

TEST(rng_address_contract_known_answer) {
    const uint32_t expected[4] = {
        UINT32_C(0x3e3b8daa), UINT32_C(0xc52865e2),
        UINT32_C(0x458c381c), UINT32_C(0xccd6de94)
    };
    uint32_t actual[4];
    ASSERT_EQ(alea_rng_block_at(
        ALEA_RNG_PHILOX4X32_10, 42,
        ALEA_RNG_DOMAIN_VOLUME_RAY_DIRECTION, 17, 0, actual), 0);
    ASSERT_EQ(memcmp(actual, expected, sizeof(expected)), 0);
}

TEST(rng_seek_and_domains_are_stable) {
    alea_rng_stream_t stream;
    alea_rng_stream_init(&stream, ALEA_RNG_PHILOX4X32_10, 42,
                         ALEA_RNG_DOMAIN_VOLUME_RAY_DIRECTION, 17);
    alea_rng_stream_seek(&stream, UINT64_C(0x100000003));
    uint32_t sought, indexed, other_domain;
    ASSERT_EQ(alea_rng_stream_next_u32(&stream, &sought), 0);
    ASSERT_EQ(alea_rng_u32_at(
        ALEA_RNG_PHILOX4X32_10, 42,
        ALEA_RNG_DOMAIN_VOLUME_RAY_DIRECTION, 17,
        UINT64_C(0x100000003), &indexed), 0);
    ASSERT_EQ(sought, indexed);
    ASSERT_EQ(alea_rng_u32_at(
        ALEA_RNG_PHILOX4X32_10, 42,
        ALEA_RNG_DOMAIN_VOLUME_RAY_DISK, 17,
        UINT64_C(0x100000003), &other_domain), 0);
    ASSERT(sought != other_domain);
}

TEST(rng_uniform_is_closed_open) {
    for (uint64_t draw = 0; draw < 1000; draw++) {
        double uniform = -1.0;
        ASSERT_EQ(alea_rng_uniform_at(
            ALEA_RNG_PHILOX4X32_10, 7, ALEA_RNG_DOMAIN_TEST_ONLY,
            11, draw, &uniform), 0);
        ASSERT(uniform >= 0.0);
        ASSERT(uniform < 1.0);
    }
}

TEST(rng_uniform_conversion_endpoints_are_exact) {
    ASSERT_EQ(alea_rng_uniform32_from_u32(0), 0.0);
    ASSERT(alea_rng_uniform32_from_u32(UINT32_MAX) < 1.0);
    ASSERT_EQ(alea_rng_uniform53_from_u32(0, 0), 0.0);
    ASSERT(alea_rng_uniform53_from_u32(UINT32_MAX, UINT32_MAX) < 1.0);
    ASSERT(alea_rng_uniform53_open_from_u32(0, 0) > 0.0);
    ASSERT(alea_rng_uniform53_open_from_u32(UINT32_MAX, UINT32_MAX) < 1.0);

    double closed, open;
    ASSERT_EQ(alea_rng_uniform53_at(
        ALEA_RNG_PHILOX4X32_10, 7, ALEA_RNG_DOMAIN_TEST_ONLY, 11, 3,
        &closed), 0);
    ASSERT_EQ(alea_rng_uniform_open_at(
        ALEA_RNG_PHILOX4X32_10, 7, ALEA_RNG_DOMAIN_TEST_ONLY, 11, 3,
        &open), 0);
    ASSERT(open > closed);
}

TEST(rng_algorithm_names_are_receipt_safe) {
    ASSERT_EQ(strcmp(alea_rng_algorithm_name(ALEA_RNG_PHILOX4X32_10),
                     "philox4x32-10"), 0);
    ASSERT_EQ(strcmp(alea_rng_algorithm_name(ALEA_RNG_LEGACY_LCG),
                     "legacy-lcg32"), 0);
}

TEST_MAIN()
