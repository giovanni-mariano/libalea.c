// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

#ifndef ALEA_RNG_DISTRIBUTION_INTERNAL_H
#define ALEA_RNG_DISTRIBUTION_INTERNAL_H

#include "rng/alea_rng.h"

#include <stddef.h>
#include <stdint.h>

#define ALEA_RNG_EVENT_ADDRESS_VERSION 1u

/* An event owns 2^32 local 32-bit draws. The collision/event ordinal is kept
 * in the upper half of the RNG draw address, so retries cannot perturb another
 * event or domain. */
typedef struct {
    alea_rng_algorithm_t algorithm;
    uint64_t seed;
    uint32_t domain;
    uint64_t entity_id;
    uint32_t event_index;
    uint64_t local_draw;
    uint32_t cached_block[4];
    uint64_t cached_block_index;
    int cached_block_valid;
} alea_rng_event_t;

typedef struct alea_rng_discrete_table alea_rng_discrete_table_t;
typedef struct alea_rng_tabular_table alea_rng_tabular_table_t;

typedef enum {
    ALEA_RNG_TABULAR_HISTOGRAM = 1,
    ALEA_RNG_TABULAR_LIN_LIN = 2
} alea_rng_tabular_interpolation_t;

/* Phase-1 transport identity layout. Particle ordinals must be allocated in
 * deterministic creation order within one primary history. */
uint64_t alea_rng_transport_entity_id(uint32_t history_id,
                                      uint32_t particle_ordinal);

int alea_rng_event_init(alea_rng_event_t* event,
                        alea_rng_algorithm_t algorithm,
                        uint64_t seed,
                        uint32_t domain,
                        uint64_t entity_id,
                        uint32_t event_index);
uint64_t alea_rng_event_position(const alea_rng_event_t* event);
int alea_rng_event_seek(alea_rng_event_t* event, uint64_t local_draw);
int alea_rng_event_next_u32(alea_rng_event_t* event, uint32_t* output);
/* Both uniform operations consume exactly two local 32-bit draws. */
int alea_rng_event_next_uniform(alea_rng_event_t* event, double* output);
int alea_rng_event_next_uniform_open(alea_rng_event_t* event, double* output);

/* Consumes one draw per rejection attempt and changes no state on failure. */
int alea_rng_sample_u32_bounded(alea_rng_event_t* event,
                                uint32_t bound,
                                uint32_t maximum_attempts,
                                uint32_t* output);
int alea_rng_sample_uniform_range(alea_rng_event_t* event,
                                  double lower,
                                  double upper,
                                  double* output);
int alea_rng_sample_exponential(alea_rng_event_t* event,
                                double rate,
                                double* output);
/* Direction outputs are unit vectors in global coordinates. Isotropic and
 * +Z cosine-law sampling each consume exactly four local draws. */
int alea_rng_sample_isotropic_direction(alea_rng_event_t* event,
                                        double direction[3]);
int alea_rng_sample_cosine_direction(alea_rng_event_t* event,
                                     double direction[3]);

alea_rng_discrete_table_t* alea_rng_discrete_table_create(
    const double* weights, size_t count);
void alea_rng_discrete_table_destroy(alea_rng_discrete_table_t* table);
size_t alea_rng_discrete_table_count(const alea_rng_discrete_table_t* table);
int alea_rng_sample_discrete_cdf(alea_rng_event_t* event,
                                 const alea_rng_discrete_table_t* table,
                                 size_t* output_index);
/* Alias sampling consumes bounded-column attempts followed by one draw. */
int alea_rng_sample_discrete_alias(alea_rng_event_t* event,
                                   const alea_rng_discrete_table_t* table,
                                   size_t* output_index);

alea_rng_tabular_table_t* alea_rng_tabular_table_create(
    const double* x,
    const double* pdf,
    size_t point_count,
    alea_rng_tabular_interpolation_t interpolation);
void alea_rng_tabular_table_destroy(alea_rng_tabular_table_t* table);
/* Inverse-CDF tabular sampling consumes exactly two local draws. */
int alea_rng_sample_tabular(alea_rng_event_t* event,
                            const alea_rng_tabular_table_t* table,
                            double* output);

#endif
