// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0
#ifndef ALEA_TALLY_INTERNAL_H
#define ALEA_TALLY_INTERNAL_H

#include "alea_tally.h"
#include "alea_raycast.h"
#include "alea_nucdata.h"

int alea_tally_plan_matches_geometry(const alea_tally_plan_t* plan,
    const alea_system_t* sys);
int alea_tally_plan_has_track_scores(const alea_tally_plan_t* plan);
int alea_tally_plan_has_collision_scores(const alea_tally_plan_t* plan);
int alea_tally_plan_has_deposition_scores(const alea_tally_plan_t* plan);
int alea_tally_plan_needs_target(const alea_tally_plan_t* plan);
alea_tally_results_t* alea_tally_results_create(
    const alea_tally_plan_t* plan);
alea_error_t alea_tally_record_track(alea_tally_results_t* results,
    const alea_tally_plan_t* plan, const alea_nav_location_t* location,
    const double start[3], const double direction[3], double distance,
    double weight, double energy, double time, double speed,
    uint32_t particle_mask, const alea_nuc_evaluation_t* evaluation);
alea_error_t alea_tally_record_collision(alea_tally_results_t* results,
    const alea_tally_plan_t* plan, const alea_nav_location_t* location,
    const double position[3], double weight, double energy, double time,
    uint32_t particle_mask, int reaction_mt, int nuclide_zaid);
alea_error_t alea_tally_record_deposition(alea_tally_results_t* results,
    const alea_tally_plan_t* plan, const alea_nav_location_t* location,
    const double position[3], double weight, double energy, double time,
    uint32_t particle_mask, int reaction_mt, int nuclide_zaid,
    double deposited_energy);
alea_error_t alea_tally_commit_history(alea_tally_results_t* results);
void alea_tally_results_finalize(alea_tally_results_t* results);

#endif
