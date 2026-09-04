// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/* Philox4x32-10 follows Salmon et al., "Parallel Random Numbers: As Easy as
 * 1, 2, 3" (SC11) and the Random123 reference constants. */

#include "rng/alea_rng.h"

#include <stddef.h>

#define PHILOX_M0 UINT32_C(0xd2511f53)
#define PHILOX_M1 UINT32_C(0xcd9e8d57)
#define PHILOX_W0 UINT32_C(0x9e3779b9)
#define PHILOX_W1 UINT32_C(0xbb67ae85)

static uint64_t rng_splitmix64(uint64_t value) {
    value += UINT64_C(0x9e3779b97f4a7c15);
    value = (value ^ (value >> 30)) * UINT64_C(0xbf58476d1ce4e5b9);
    value = (value ^ (value >> 27)) * UINT64_C(0x94d049bb133111eb);
    return value ^ (value >> 31);
}

static void philox_mulhilo(uint32_t multiplier, uint32_t value,
                           uint32_t* high, uint32_t* low) {
    const uint64_t product = (uint64_t)multiplier * (uint64_t)value;
    *low = (uint32_t)product;
    *high = (uint32_t)(product >> 32);
}

void alea_philox4x32_10(const uint32_t counter[4],
                        const uint32_t key[2],
                        uint32_t output[4]) {
    uint32_t c0 = counter[0], c1 = counter[1];
    uint32_t c2 = counter[2], c3 = counter[3];
    uint32_t k0 = key[0], k1 = key[1];

    for (int round = 0; round < 10; round++) {
        uint32_t hi0, lo0, hi1, lo1;
        philox_mulhilo(PHILOX_M0, c0, &hi0, &lo0);
        philox_mulhilo(PHILOX_M1, c2, &hi1, &lo1);
        const uint32_t next0 = hi1 ^ c1 ^ k0;
        const uint32_t next1 = lo1;
        const uint32_t next2 = hi0 ^ c3 ^ k1;
        const uint32_t next3 = lo0;
        c0 = next0;
        c1 = next1;
        c2 = next2;
        c3 = next3;
        if (round != 9) {
            k0 += PHILOX_W0;
            k1 += PHILOX_W1;
        }
    }
    output[0] = c0;
    output[1] = c1;
    output[2] = c2;
    output[3] = c3;
}

static void rng_domain_key(uint64_t seed, uint32_t domain, uint32_t key[2]) {
    const uint64_t mixed = rng_splitmix64(
        seed ^ (UINT64_C(0xd2b74407b1ce6e93) * (uint64_t)domain));
    key[0] = (uint32_t)mixed;
    key[1] = (uint32_t)(mixed >> 32);
}

static uint32_t legacy_lcg_at(uint64_t seed, uint64_t entity_id,
                              uint64_t draw_index) {
    const uint64_t mixed = rng_splitmix64(
        seed ^ rng_splitmix64(entity_id));
    uint32_t state = (uint32_t)(mixed ^ (mixed >> 32));
    for (uint64_t draw = 0; draw <= draw_index; draw++)
        state = state * UINT32_C(1664525) + UINT32_C(1013904223);
    return state;
}

int alea_rng_block_at(alea_rng_algorithm_t algorithm,
                      uint64_t seed,
                      uint32_t domain,
                      uint64_t entity_id,
                      uint64_t block_index,
                      uint32_t output[4]) {
    if (!output) return -1;
    if (algorithm == ALEA_RNG_LEGACY_LCG) {
        if (block_index > UINT64_MAX / UINT64_C(4)) return -1;
        const uint64_t first = block_index * UINT64_C(4);
        for (uint64_t lane = 0; lane < 4; lane++)
            output[lane] = legacy_lcg_at(seed, entity_id, first + lane);
        return 0;
    }
    if (algorithm != ALEA_RNG_PHILOX4X32_10) return -1;

    const uint32_t counter[4] = {
        (uint32_t)entity_id,
        (uint32_t)(entity_id >> 32),
        (uint32_t)block_index,
        (uint32_t)(block_index >> 32)
    };
    uint32_t key[2];
    rng_domain_key(seed, domain, key);
    alea_philox4x32_10(counter, key, output);
    return 0;
}

int alea_rng_u32_at(alea_rng_algorithm_t algorithm,
                    uint64_t seed,
                    uint32_t domain,
                    uint64_t entity_id,
                    uint64_t draw_index,
                    uint32_t* output) {
    uint32_t block[4];
    if (!output || alea_rng_block_at(algorithm, seed, domain, entity_id,
                                      draw_index >> 2, block) != 0) return -1;
    *output = block[draw_index & UINT64_C(3)];
    return 0;
}

int alea_rng_u64_at(alea_rng_algorithm_t algorithm,
                    uint64_t seed,
                    uint32_t domain,
                    uint64_t entity_id,
                    uint64_t draw_index,
                    uint64_t* output) {
    uint32_t high, low;
    if (!output || draw_index > UINT64_MAX / UINT64_C(2) ||
        alea_rng_u32_at(algorithm, seed, domain, entity_id,
                        draw_index * UINT64_C(2), &high) != 0 ||
        alea_rng_u32_at(algorithm, seed, domain, entity_id,
                        draw_index * UINT64_C(2) + 1, &low) != 0) return -1;
    *output = ((uint64_t)high << 32) | (uint64_t)low;
    return 0;
}

double alea_rng_uniform32_from_u32(uint32_t bits) {
    return (double)bits * 0x1p-32;
}

double alea_rng_uniform53_from_u32(uint32_t high, uint32_t low) {
    const uint64_t significand =
        ((uint64_t)(high >> 5) << 26) | (uint64_t)(low >> 6);
    return (double)significand * 0x1p-53;
}

double alea_rng_uniform53_open_from_u32(uint32_t high, uint32_t low) {
    const uint64_t significand =
        ((uint64_t)(high >> 5) << 26) | (uint64_t)(low >> 6);
    /* Avoid both endpoints without midpoint rounding turning the largest
     * significand into 1.0 under round-to-nearest. */
    if (significand == ((UINT64_C(1) << 53) - 1))
        return 0x1.fffffffffffffp-1;
    return (double)(significand + 1) * 0x1p-53;
}

int alea_rng_uniform_at(alea_rng_algorithm_t algorithm,
                        uint64_t seed,
                        uint32_t domain,
                        uint64_t entity_id,
                        uint64_t draw_index,
                        double* output) {
    uint32_t bits;
    if (!output || alea_rng_u32_at(
            algorithm, seed, domain, entity_id, draw_index, &bits) != 0) {
        return -1;
    }
    *output = alea_rng_uniform32_from_u32(bits);
    return 0;
}

static int rng_uniform_pair_at(alea_rng_algorithm_t algorithm,
                               uint64_t seed, uint32_t domain,
                               uint64_t entity_id, uint64_t draw_index,
                               uint32_t* high, uint32_t* low) {
    if (draw_index > UINT64_MAX / UINT64_C(2)) return -1;
    const uint64_t first = draw_index * UINT64_C(2);
    return alea_rng_u32_at(algorithm, seed, domain, entity_id, first, high) != 0 ||
           alea_rng_u32_at(algorithm, seed, domain, entity_id, first + 1, low) != 0
               ? -1 : 0;
}

int alea_rng_uniform53_at(alea_rng_algorithm_t algorithm,
                          uint64_t seed,
                          uint32_t domain,
                          uint64_t entity_id,
                          uint64_t draw_index,
                          double* output) {
    uint32_t high, low;
    if (!output || rng_uniform_pair_at(algorithm, seed, domain, entity_id,
                                        draw_index, &high, &low) != 0) return -1;
    *output = alea_rng_uniform53_from_u32(high, low);
    return 0;
}

int alea_rng_uniform_open_at(alea_rng_algorithm_t algorithm,
                             uint64_t seed,
                             uint32_t domain,
                             uint64_t entity_id,
                             uint64_t draw_index,
                             double* output) {
    uint32_t high, low;
    if (!output || rng_uniform_pair_at(algorithm, seed, domain, entity_id,
                                        draw_index, &high, &low) != 0) return -1;
    *output = alea_rng_uniform53_open_from_u32(high, low);
    return 0;
}

void alea_rng_stream_init(alea_rng_stream_t* stream,
                          alea_rng_algorithm_t algorithm,
                          uint64_t seed,
                          uint32_t domain,
                          uint64_t entity_id) {
    if (!stream) return;
    stream->seed = seed;
    stream->algorithm = algorithm;
    stream->entity_id = entity_id;
    stream->domain = domain;
    stream->draw_index = 0;
    stream->cached_block_index = 0;
    stream->cached_block_valid = 0;
}

int alea_rng_stream_next_u32(alea_rng_stream_t* stream, uint32_t* output) {
    if (!stream || !output || stream->draw_index == UINT64_MAX) return -1;
    const uint64_t block_index = stream->draw_index >> 2;
    if (!stream->cached_block_valid ||
        stream->cached_block_index != block_index) {
        if (alea_rng_block_at(stream->algorithm, stream->seed, stream->domain,
                              stream->entity_id, block_index,
                              stream->cached_block) != 0) return -1;
        stream->cached_block_index = block_index;
        stream->cached_block_valid = 1;
    }
    *output = stream->cached_block[stream->draw_index & UINT64_C(3)];
    stream->draw_index++;
    return 0;
}

int alea_rng_stream_next_u64(alea_rng_stream_t* stream, uint64_t* output) {
    uint32_t high, low;
    if (!output || alea_rng_stream_next_u32(stream, &high) != 0 ||
        alea_rng_stream_next_u32(stream, &low) != 0) return -1;
    *output = ((uint64_t)high << 32) | (uint64_t)low;
    return 0;
}

int alea_rng_stream_next_uniform(alea_rng_stream_t* stream, double* output) {
    uint32_t bits;
    if (!output || alea_rng_stream_next_u32(stream, &bits) != 0) return -1;
    *output = alea_rng_uniform32_from_u32(bits);
    return 0;
}

uint64_t alea_rng_stream_position(const alea_rng_stream_t* stream) {
    return stream ? stream->draw_index : 0;
}

void alea_rng_stream_seek(alea_rng_stream_t* stream, uint64_t draw_index) {
    if (stream) {
        stream->draw_index = draw_index;
        if (!stream->cached_block_valid ||
            stream->cached_block_index != (draw_index >> 2))
            stream->cached_block_valid = 0;
    }
}

const char* alea_rng_algorithm_name(alea_rng_algorithm_t algorithm) {
    switch (algorithm) {
        case ALEA_RNG_PHILOX4X32_10: return "philox4x32-10";
        case ALEA_RNG_LEGACY_LCG: return "legacy-lcg32";
        default: return "unknown";
    }
}
