// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/** @file alea_tally.h
 * Spatial tally definitions and source-history moments for transport.
 */
#ifndef ALEA_TALLY_H
#define ALEA_TALLY_H

#include "alea.h"
#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct alea_tally_plan alea_tally_plan_t;
typedef struct alea_tally_results alea_tally_results_t;

typedef enum {
    ALEA_TALLY_CELL,
    ALEA_TALLY_UNIVERSE,
    ALEA_TALLY_CARTESIAN_MESH
} alea_tally_domain_t;

typedef enum {
    ALEA_TALLY_TRACK_LENGTH,
    ALEA_TALLY_COLLISION,
    ALEA_TALLY_REACTION_EVENT,
    ALEA_TALLY_REACTION_RATE,
    ALEA_TALLY_HEATING,
    ALEA_TALLY_LOCAL_DEPOSITION
} alea_tally_score_t;

#define ALEA_TALLY_NEUTRON (1u << 0)
#define ALEA_TALLY_PHOTON  (1u << 1)

typedef struct {
    alea_tally_score_t score;
    alea_tally_domain_t domain;
    /* Only for CARTESIAN_MESH. World coordinates, centimetres. Bins are
     * half-open and flattened with x fastest: (z*ny + y)*nx + x. */
    double lower[3];
    double upper[3];
    uint32_t dimensions[3];
    uint32_t particle_mask; /* zero: all particles */
    double energy_min; /* both energy bounds zero: all energies; else [min,max) */
    double energy_max;
    double time_min; /* both time bounds zero: all times; else [min,max), s */
    double time_max;
    int material_id; /* zero: all; -1: void; positive: specific material */
    int reaction_mt; /* event/rate/local deposition; zero: all reactions */
    int nuclide_zaid; /* event/rate/heating or collision; 1000*Z+A, zero: all */
    const double* energy_edges; /* optional, group_count+1 increasing MeV edges */
    size_t energy_group_count; /* zero: no energy axis */
} alea_tally_spec_t;

/** The geometry must outlive the plan and must not change during a run. */
alea_tally_plan_t* alea_tally_plan_create(const alea_system_t* sys);
void alea_tally_plan_free(alea_tally_plan_t* plan);
alea_error_t alea_tally_plan_add(alea_tally_plan_t* plan,
    const alea_tally_spec_t* spec, size_t* index);
size_t alea_tally_plan_count(const alea_tally_plan_t* plan);

typedef struct {
    alea_tally_score_t score;
    alea_tally_domain_t domain;
    size_t bin_count; /* spatial_bin_count * max(1, energy_group_count) */
    size_t spatial_bin_count;
    size_t energy_group_count;
    const double* energy_edges; /* copied into results; NULL without axis */
    uint32_t histories;
    const int* bin_ids; /* cell or terminal universe IDs; NULL for mesh */
    const double* sum; /* units depend on score; source-history sum */
    const double* sum_squared; /* squares of complete source-history scores */
    double lower[3];
    double upper[3];
    uint32_t dimensions[3];
    uint32_t particle_mask;
    double energy_min;
    double energy_max;
    double time_min;
    double time_max;
    int material_id;
    int reaction_mt;
    int nuclide_zaid;
} alea_tally_view_t;

size_t alea_tally_results_count(const alea_tally_results_t* results);
/** View and arrays are borrowed from results. Mean per source history is
 * sum[bin]/histories; sum_squared permits uncertainty calculation. */
alea_error_t alea_tally_results_view(const alea_tally_results_t* results,
    size_t index, alea_tally_view_t* view);
void alea_tally_results_free(alea_tally_results_t* results);

#ifdef __cplusplus
}
#endif

#endif
