// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

#ifndef ALEA_GEO_VALIDATOR_INTERNAL_H
#define ALEA_GEO_VALIDATOR_INTERNAL_H

#include "alea_geo_validator.h"

int alea_validator_cluster_source_init(alea_system_t* sys,
    const alea_geom_validator_options_t* supplied,
    alea_geom_validator_options_t* prepared,
    double bounds[6], double* t_max, uint64_t* rng);
void alea_validator_cluster_next_ray(uint64_t* rng,
    const double bounds[6], double origin[3], double direction[3]);
int alea_validator_cluster_merge_one(alea_system_t* sys,
    const double origin[3], const double direction[3], double t_max,
    const alea_geom_validator_options_t* options,
    alea_geom_validator_result_t* result,
    const alea_geom_validator_result_t* candidate);

#endif /* ALEA_GEO_VALIDATOR_INTERNAL_H */
