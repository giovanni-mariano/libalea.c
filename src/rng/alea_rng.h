// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

#ifndef ALEA_RNG_INTERNAL_H
#define ALEA_RNG_INTERNAL_H

#include "alea.h"

#include <stdint.h>

typedef enum {
    ALEA_RNG_DOMAIN_VOLUME_RAY_DIRECTION = 0x00010001u,
    ALEA_RNG_DOMAIN_VOLUME_RAY_DISK = 0x00010002u,
    ALEA_RNG_DOMAIN_TRANSPORT_SOURCE_POSITION = 0x00020001u,
    ALEA_RNG_DOMAIN_TRANSPORT_SOURCE_DIRECTION = 0x00020002u,
    ALEA_RNG_DOMAIN_TRANSPORT_SOURCE_ENERGY = 0x00020003u,
    ALEA_RNG_DOMAIN_TRANSPORT_FREE_PATH = 0x00020004u,
    ALEA_RNG_DOMAIN_TRANSPORT_NUCLIDE = 0x00020005u,
    ALEA_RNG_DOMAIN_TRANSPORT_REACTION = 0x00020006u,
    ALEA_RNG_DOMAIN_TRANSPORT_SCATTER_ANGLE = 0x00020007u,
    ALEA_RNG_DOMAIN_TRANSPORT_SECONDARY_ENERGY = 0x00020008u,
    ALEA_RNG_DOMAIN_TRANSPORT_ROULETTE = 0x00020009u,
    ALEA_RNG_DOMAIN_TRANSPORT_WEIGHT_WINDOW = 0x0002000au,
    ALEA_RNG_DOMAIN_RENDER_STOCHASTIC_SAMPLE = 0x00030001u,
    ALEA_RNG_DOMAIN_TEST_ONLY = 0x7fff0001u
} alea_rng_domain_t;

typedef struct {
    uint64_t seed;
    alea_rng_algorithm_t algorithm;
    uint64_t entity_id;
    uint32_t domain;
    uint64_t draw_index;
    uint32_t cached_block[4];
    uint64_t cached_block_index;
    int cached_block_valid;
} alea_rng_stream_t;

void alea_philox4x32_10(const uint32_t counter[4],
                        const uint32_t key[2],
                        uint32_t output[4]);

int alea_rng_block_at(alea_rng_algorithm_t algorithm,
                      uint64_t seed,
                      uint32_t domain,
                      uint64_t entity_id,
                      uint64_t block_index,
                      uint32_t output[4]);

int alea_rng_u32_at(alea_rng_algorithm_t algorithm,
                    uint64_t seed,
                    uint32_t domain,
                    uint64_t entity_id,
                    uint64_t draw_index,
                    uint32_t* output);

int alea_rng_u64_at(alea_rng_algorithm_t algorithm,
                    uint64_t seed,
                    uint32_t domain,
                    uint64_t entity_id,
                    uint64_t draw_index,
                    uint64_t* output);

double alea_rng_uniform32_from_u32(uint32_t bits);
double alea_rng_uniform53_from_u32(uint32_t high, uint32_t low);
double alea_rng_uniform53_open_from_u32(uint32_t high, uint32_t low);

int alea_rng_uniform_at(alea_rng_algorithm_t algorithm,
                        uint64_t seed,
                        uint32_t domain,
                        uint64_t entity_id,
                        uint64_t draw_index,
                        double* output);

int alea_rng_uniform53_at(alea_rng_algorithm_t algorithm,
                          uint64_t seed,
                          uint32_t domain,
                          uint64_t entity_id,
                          uint64_t draw_index,
                          double* output);

int alea_rng_uniform_open_at(alea_rng_algorithm_t algorithm,
                             uint64_t seed,
                             uint32_t domain,
                             uint64_t entity_id,
                             uint64_t draw_index,
                             double* output);

void alea_rng_stream_init(alea_rng_stream_t* stream,
                          alea_rng_algorithm_t algorithm,
                          uint64_t seed,
                          uint32_t domain,
                          uint64_t entity_id);
int alea_rng_stream_next_u32(alea_rng_stream_t* stream, uint32_t* output);
int alea_rng_stream_next_u64(alea_rng_stream_t* stream, uint64_t* output);
int alea_rng_stream_next_uniform(alea_rng_stream_t* stream, double* output);
uint64_t alea_rng_stream_position(const alea_rng_stream_t* stream);
void alea_rng_stream_seek(alea_rng_stream_t* stream, uint64_t draw_index);

#endif
