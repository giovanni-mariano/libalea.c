// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/*
 * test_bitset.c - Unit tests for alea_bitset
 */

#include "alea_test.h"
#include "util/alea_bitset.h"

TEST(bitset_create_and_destroy) {
    alea_bitset_t bs = alea_bitset_create(100);
    ASSERT_NOT_NULL(bs.words);
    ASSERT_EQ(bs.n_words, 2);  /* ceil(100/64) = 2 */

    /* All bits should start cleared */
    for (size_t i = 0; i < 100; i++) {
        ASSERT_FALSE(alea_bitset_test(&bs, i));
    }

    alea_bitset_destroy(&bs);
    ASSERT_NULL(bs.words);
    ASSERT_EQ(bs.n_words, 0);
}

TEST(bitset_create_zero) {
    alea_bitset_t bs = alea_bitset_create(0);
    ASSERT_NOT_NULL(bs.words);
    ASSERT_EQ(bs.n_words, 1);  /* minimum 1 word */
    alea_bitset_destroy(&bs);
}

TEST(bitset_set_test_clear) {
    alea_bitset_t bs = alea_bitset_create(200);

    alea_bitset_set(&bs, 42);
    ASSERT_TRUE(alea_bitset_test(&bs, 42));
    ASSERT_FALSE(alea_bitset_test(&bs, 41));
    ASSERT_FALSE(alea_bitset_test(&bs, 43));

    alea_bitset_clear(&bs, 42);
    ASSERT_FALSE(alea_bitset_test(&bs, 42));

    alea_bitset_destroy(&bs);
}

TEST(bitset_edge_bit0) {
    alea_bitset_t bs = alea_bitset_create(64);

    alea_bitset_set(&bs, 0);
    ASSERT_TRUE(alea_bitset_test(&bs, 0));
    ASSERT_FALSE(alea_bitset_test(&bs, 1));

    alea_bitset_destroy(&bs);
}

TEST(bitset_edge_bit63) {
    alea_bitset_t bs = alea_bitset_create(64);

    alea_bitset_set(&bs, 63);
    ASSERT_TRUE(alea_bitset_test(&bs, 63));
    ASSERT_FALSE(alea_bitset_test(&bs, 62));

    alea_bitset_destroy(&bs);
}

TEST(bitset_edge_bit64) {
    alea_bitset_t bs = alea_bitset_create(128);

    /* Bit 64 is the first bit in the second word */
    alea_bitset_set(&bs, 64);
    ASSERT_TRUE(alea_bitset_test(&bs, 64));
    ASSERT_FALSE(alea_bitset_test(&bs, 63));
    ASSERT_FALSE(alea_bitset_test(&bs, 65));

    alea_bitset_destroy(&bs);
}

TEST(bitset_popcount) {
    alea_bitset_t bs = alea_bitset_create(256);

    ASSERT_EQ(alea_bitset_popcount(&bs), 0);

    alea_bitset_set(&bs, 0);
    alea_bitset_set(&bs, 63);
    alea_bitset_set(&bs, 64);
    alea_bitset_set(&bs, 127);
    alea_bitset_set(&bs, 200);
    ASSERT_EQ(alea_bitset_popcount(&bs), 5);

    alea_bitset_clear(&bs, 63);
    ASSERT_EQ(alea_bitset_popcount(&bs), 4);

    alea_bitset_destroy(&bs);
}

TEST(bitset_clear_all) {
    alea_bitset_t bs = alea_bitset_create(128);

    for (size_t i = 0; i < 128; i++) {
        alea_bitset_set(&bs, i);
    }
    ASSERT_EQ(alea_bitset_popcount(&bs), 128);

    alea_bitset_clear_all(&bs);
    ASSERT_EQ(alea_bitset_popcount(&bs), 0);

    for (size_t i = 0; i < 128; i++) {
        ASSERT_FALSE(alea_bitset_test(&bs, i));
    }

    alea_bitset_destroy(&bs);
}

TEST(bitset_large) {
    alea_bitset_t bs = alea_bitset_create(100000);
    ASSERT_NOT_NULL(bs.words);

    /* Set every 100th bit */
    for (size_t i = 0; i < 100000; i += 100) {
        alea_bitset_set(&bs, i);
    }
    ASSERT_EQ(alea_bitset_popcount(&bs), 1000);

    /* Verify */
    for (size_t i = 0; i < 100000; i++) {
        if (i % 100 == 0) {
            ASSERT_TRUE(alea_bitset_test(&bs, i));
        } else {
            ASSERT_FALSE(alea_bitset_test(&bs, i));
        }
    }

    alea_bitset_destroy(&bs);
}

TEST(bitset_exact_word_boundary) {
    /* Exactly 64 bits = 1 word */
    alea_bitset_t bs = alea_bitset_create(64);
    ASSERT_EQ(bs.n_words, 1);

    /* 65 bits = 2 words */
    alea_bitset_t bs2 = alea_bitset_create(65);
    ASSERT_EQ(bs2.n_words, 2);

    alea_bitset_destroy(&bs);
    alea_bitset_destroy(&bs2);
}

TEST_MAIN()
