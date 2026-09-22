// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/** @file alea_adjoint.h
 * Steady-state, isotropic multigroup adjoint transport. All cross sections
 * are macroscopic (1/cm); group sources are particles/(s cm^3), integrated
 * over physical direction. The adjoint source is the detector response.
 */
#ifndef ALEA_ADJOINT_H
#define ALEA_ADJOINT_H

#include "alea_transport.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    const double* total;    /* n_groups entries; macroscopic total, 1/cm */
    const double* transfer; /* n_groups^2; forward expected production,
                             * transfer[incoming * n_groups + outgoing], 1/cm */
} alea_adjoint_material_t;

typedef struct {
    double position[3];
    double physical_direction[3]; /* adjoint walker travels in its negative */
    size_t group;
    double weight; /* detector density / sampling PDF, including angle measure */
} alea_adjoint_detector_particle_t;

typedef alea_error_t (*alea_adjoint_detector_sampler_fn)(
    void* context, uint64_t seed, uint32_t history_id,
    alea_adjoint_detector_particle_t* output);

typedef struct {
    double lower[3], upper[3]; /* detector box, cm; must lie inside geometry */
    size_t group;
    double response; /* uniform response coefficient, 1/cm */
} alea_adjoint_box_detector_t;

/** Isotropic box detector sampler. Output weight includes box volume and 4pi.
 * The caller must ensure samples are in the intended detector region. */
alea_error_t alea_adjoint_sample_box_detector(
    void* context, uint64_t seed, uint32_t history_id,
    alea_adjoint_detector_particle_t* output);

typedef struct {
    alea_nuc_particle_t particle; /* neutron or photon; no coupling yet */
    size_t n_groups;
    size_t cell_count; /* must equal alea_cell_count(sys) */
    const alea_adjoint_material_t* cell_materials; /* cell_count; NULL in void */
    const double* physical_source; /* cell_count*n_groups, isotropic strength */
    alea_adjoint_detector_sampler_fn detector_sampler;
    void* detector_context;
} alea_adjoint_problem_t;

typedef struct {
    uint32_t histories;
    uint32_t history_offset;
    uint64_t seed;
    uint32_t max_events_per_history;
    double max_segment_distance; /* finite, positive cm */
    alea_nav_validation_mode_t navigation_validation;
    size_t max_navigation_breakpoints;
} alea_adjoint_options_t;

typedef struct {
    uint32_t histories;
    double sum;         /* sum of independent detector-history response scores */
    double sum_squared; /* sum of squared complete history scores */
    double mean;        /* physical detector response, per physical-source second */
    double standard_error;
} alea_adjoint_result_t;

/** Run independent detector-launched histories. On error, output is zeroed;
 * failure describes the first incomplete history. Unsupported boundaries
 * return ALEA_ERR_UNSUPPORTED. Materials and source arrays are borrowed.
 * The geometry must not change during the run. */
alea_error_t alea_adjoint_run(
    alea_system_t* sys, const alea_adjoint_problem_t* problem,
    const alea_adjoint_options_t* options, alea_adjoint_result_t* output,
    alea_transport_failure_t* failure);

#ifdef __cplusplus
}
#endif

#endif
