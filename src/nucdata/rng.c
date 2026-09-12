// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

#include "alea_nucdata.h"
#include "rng/alea_rng_distribution.h"
#include <math.h>

alea_error_t alea_nuc_rng_init(alea_nuc_rng_t* rng, uint64_t seed,
    uint32_t history_id, uint32_t particle_ordinal, uint32_t event_index,
    alea_nuc_rng_domain_t domain) {
    if (!rng) return ALEA_ERR_NULL_ARG;
    uint32_t address_domain;
    switch (domain) {
    case ALEA_NUC_RNG_FLIGHT:
        address_domain = ALEA_RNG_DOMAIN_TRANSPORT_FREE_PATH; break;
    case ALEA_NUC_RNG_COLLISION:
        address_domain = ALEA_RNG_DOMAIN_TRANSPORT_REACTION; break;
    case ALEA_NUC_RNG_URR:
        address_domain = ALEA_RNG_DOMAIN_TRANSPORT_URR; break;
    default: return ALEA_ERR_INVALID_ARG;
    }
    *rng = (alea_nuc_rng_t){seed,
        alea_rng_transport_entity_id(history_id, particle_ordinal),
        0, event_index, address_domain};
    return ALEA_OK;
}

double alea_nuc_rng_uniform(void* context) {
    alea_nuc_rng_t* rng = context;
    if (!rng || rng->local_draw > UINT64_C(0x100000000) - 2 ||
        (rng->domain != ALEA_RNG_DOMAIN_TRANSPORT_FREE_PATH &&
         rng->domain != ALEA_RNG_DOMAIN_TRANSPORT_REACTION &&
         rng->domain != ALEA_RNG_DOMAIN_TRANSPORT_URR)) return NAN;
    double value;
    uint64_t address = ((uint64_t)rng->event_index << 32) | rng->local_draw;
    if (alea_rng_uniform53_at(ALEA_RNG_PHILOX4X32_10, rng->seed,
            rng->domain, rng->entity_id, address, &value) != 0) return NAN;
    rng->local_draw += 2;
    return value;
}
