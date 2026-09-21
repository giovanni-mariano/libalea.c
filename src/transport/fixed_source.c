// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

#include "alea_transport.h"
#include "tally_internal.h"
#include "alea.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

#define NEUTRON_REST_MEV 939.56542052
#define LIGHT_SPEED_CM_S 2.99792458e10
#define DEFAULT_BANK_CAPACITY 1024u

typedef struct {
    double position[3];
    alea_nuc_particle_state_t particle;
    uint32_t ordinal;
} pending_particle_t;

void alea_transport_result_free(alea_transport_result_t* result) {
    if (!result) return;
    free(result->track_length);
    free(result->track_length_squared);
    alea_tally_results_free(result->tallies);
    memset(result, 0, sizeof(*result));
}

static int normalize_direction(double direction[3]) {
    double scale = fmax(fabs(direction[0]),
                        fmax(fabs(direction[1]), fabs(direction[2])));
    if (!isfinite(scale) || scale <= 0.0) return -1;
    double x = direction[0] / scale;
    double y = direction[1] / scale;
    double z = direction[2] / scale;
    double length = sqrt(x*x + y*y + z*z);
    if (!isfinite(length) || length <= 0.0) return -1;
    direction[0] = x / length;
    direction[1] = y / length;
    direction[2] = z / length;
    return 0;
}

static double neutron_speed(double energy) {
    const double mass = NEUTRON_REST_MEV;
    return LIGHT_SPEED_CM_S * sqrt(energy * (energy + 2.0 * mass)) /
           (energy + mass);
}

static alea_error_t draw_optical_depth(uint64_t seed, uint32_t history,
                                      uint32_t particle_ordinal,
                                      uint32_t collision_index, double* tau) {
    alea_nuc_rng_t rng;
    alea_error_t err = alea_nuc_rng_init(&rng, seed, history, particle_ordinal,
                                        collision_index, ALEA_NUC_RNG_FLIGHT);
    if (err != ALEA_OK) return err;
    double u = alea_nuc_rng_uniform(&rng);
    if (!isfinite(u) || u < 0.0 || u >= 1.0) return ALEA_ERR_INVALID_STATE;
    *tau = -log1p(-u);
    return isfinite(*tau) && *tau >= 0.0 ? ALEA_OK : ALEA_ERR_INVALID_STATE;
}

static alea_error_t fail_transport(alea_transport_failure_t* failure,
    uint32_t history, uint32_t event_index, const alea_nav_location_t* location,
    const double position[3], double energy, alea_error_t error) {
    if (failure) {
        failure->history_id = history;
        failure->event_index = event_index;
        failure->cell_id = location ? location->cell_id : -1;
        failure->energy = energy;
        failure->error = error;
        memcpy(failure->position, position, sizeof(failure->position));
    }
    alea_set_error_detail(error,
        "transport history %u event %u cell %d at (%.9g, %.9g, %.9g), E=%.9g MeV",
        history, event_index, location ? location->cell_id : -1,
        position[0], position[1], position[2], energy);
    return error;
}

alea_error_t alea_transport_run_fixed_source(
    alea_system_t* sys, const alea_nuc_cell_bindings_t* bindings,
    const alea_transport_source_t* source,
    const alea_transport_options_t* options,
    alea_transport_result_t* output,
    alea_transport_failure_t* failure) {
    if (!output) return ALEA_ERR_NULL_ARG;
    memset(output, 0, sizeof(*output));
    if (failure) memset(failure, 0, sizeof(*failure));
    if (!sys || !bindings || !source || !options) return ALEA_ERR_NULL_ARG;
    if (!alea_nuc_cell_bindings_matches_geometry(bindings, sys))
        return ALEA_ERR_INVALID_STATE;
    if (options->tally_plan &&
        !alea_tally_plan_matches_geometry(options->tally_plan, sys))
        return ALEA_ERR_INVALID_STATE;
    if (options->histories == 0 || options->max_events_per_history == 0 ||
        !isfinite(options->max_segment_distance) ||
        options->max_segment_distance <= 0.0 ||
        (source->particle.type != ALEA_NUC_PARTICLE_NEUTRON &&
         source->particle.type != ALEA_NUC_PARTICLE_PHOTON) ||
        !isfinite(source->particle.energy) || source->particle.energy <= 0.0 ||
        !isfinite(source->particle.weight) || source->particle.weight <= 0.0 ||
        !isfinite(source->particle.time) ||
        !isfinite(source->position[0]) ||
        !isfinite(source->position[1]) ||
        !isfinite(source->position[2]))
        return ALEA_ERR_INVALID_ARG;

    const size_t cell_count = alea_cell_count(sys);
    const size_t bank_capacity = options->max_pending_particles ?
        options->max_pending_particles : DEFAULT_BANK_CAPACITY;
    const size_t scratch_capacity = bank_capacity > DEFAULT_BANK_CAPACITY ?
        bank_capacity : DEFAULT_BANK_CAPACITY;
    if (cell_count > SIZE_MAX / sizeof(double) ||
        bank_capacity > SIZE_MAX / sizeof(pending_particle_t) ||
        scratch_capacity > SIZE_MAX / sizeof(alea_nuc_particle_state_t))
        return ALEA_ERR_OVERFLOW;
    double* path_sum = calloc(cell_count ? cell_count : 1, sizeof(double));
    double* path_sum_squared = calloc(cell_count ? cell_count : 1, sizeof(double));
    double* history_path = calloc(cell_count ? cell_count : 1, sizeof(double));
    pending_particle_t* bank = calloc(bank_capacity, sizeof(*bank));
    alea_nuc_particle_state_t* emitted = calloc(scratch_capacity,
                                                sizeof(*emitted));
    alea_ray_navigator_t* navigator = alea_ray_navigator_create(sys);
    alea_tally_results_t* tallies = options->tally_plan ?
        alea_tally_results_create(options->tally_plan) : NULL;
    if (!path_sum || !path_sum_squared || !history_path || !bank || !emitted ||
        !navigator || (options->tally_plan && !tallies)) {
        free(path_sum); free(path_sum_squared); free(history_path);
        free(bank); free(emitted);
        alea_tally_results_free(tallies);
        alea_ray_navigator_destroy(navigator);
        return ALEA_ERR_OUT_OF_MEMORY;
    }
    alea_transport_result_t result = {
        .cell_count = cell_count, .histories = options->histories,
        .track_length = path_sum, .track_length_squared = path_sum_squared,
        .tallies = tallies
    };
    alea_error_t err = ALEA_OK;
    for (uint32_t h = 0; h < options->histories; ++h) {
        memset(history_path, 0, cell_count * sizeof(double));
        bank[0].particle = source->particle;
        if (normalize_direction(bank[0].particle.direction) != 0) {
            err = ALEA_ERR_INVALID_ARG;
            break;
        }
        memcpy(bank[0].position, source->position, sizeof(bank[0].position));
        bank[0].ordinal = 0;
        size_t bank_count = 1;
        uint64_t next_ordinal = 1;
        uint32_t events = 0;
        alea_nav_location_t last_location = {0};
        double last_position[3];
        memcpy(last_position, source->position, sizeof(last_position));
        double last_energy = source->particle.energy;
        uint32_t last_ordinal = 0;
        while (bank_count && err == ALEA_OK) {
            pending_particle_t current = bank[--bank_count];
            alea_nuc_particle_state_t particle = current.particle;
            const uint32_t particle_ordinal = current.ordinal;
            double position[3];
            memcpy(position, current.position, sizeof(position));
            alea_nav_location_t location;
            if (alea_ray_navigator_restart(navigator, position,
                                           particle.direction, &location) != 0) {
                err = fail_transport(failure, h, 0, NULL, position,
                                     particle.energy, ALEA_ERR_INVALID_STATE);
                if (failure) failure->particle_ordinal = particle_ordinal;
                break;
            }
            double remaining_tau = 0.0;
            bool new_flight = true;
            bool finished = false;
            uint32_t collisions = 0;
            while (!finished) {
                if (events >= options->max_events_per_history) {
                    err = fail_transport(failure, h, events, &location, position,
                                         particle.energy, ALEA_ERR_OVERFLOW);
                    break;
                }
                ++events;
                if (location.kind != ALEA_NAV_MATERIAL &&
                    location.kind != ALEA_NAV_VOID) {
                    err = fail_transport(failure, h, events, &location, position,
                                         particle.energy, ALEA_ERR_INVALID_STATE);
                    break;
                }
                if (new_flight) {
                    err = draw_optical_depth(options->seed, h, particle_ordinal,
                                             collisions,
                                             &remaining_tau);
                    if (err != ALEA_OK) {
                        err = fail_transport(failure, h, events, &location,
                                             position, particle.energy, err);
                        break;
                    }
                    new_flight = false;
                }
                alea_nuc_evaluation_t evaluation;
                double sigma = 0.0;
                if (location.kind == ALEA_NAV_MATERIAL) {
                    if (location.cell_index < 0 ||
                        (size_t)location.cell_index >= cell_count) {
                        err = ALEA_ERR_INVALID_STATE;
                    } else {
                        const alea_nuc_prepared_material_t* prepared =
                            alea_nuc_cell_bindings_get(bindings,
                                (size_t)location.cell_index,
                                particle.type);
                        err = prepared ? alea_nuc_evaluate(prepared, &particle,
                            &evaluation) : ALEA_ERR_NOT_FOUND;
                        if (err == ALEA_OK) sigma = evaluation.macro_total;
                    }
                    if (err != ALEA_OK || sigma < 0.0 || !isfinite(sigma)) {
                        err = fail_transport(failure, h, events, &location,
                            position, particle.energy,
                            err == ALEA_OK ? ALEA_ERR_INVALID_STATE : err);
                        break;
                    }
                }

                bool collision_now = location.kind == ALEA_NAV_MATERIAL &&
                                     sigma > 0.0 && remaining_tau == 0.0;
                alea_nav_event_t event = {0};
                if (collision_now) {
                    event.kind = ALEA_NAV_COLLISION;
                    event.before = location;
                    event.after = location;
                    memcpy(event.position, position, sizeof(position));
                } else {
                    double distance = sigma > 0.0 ? remaining_tau / sigma : INFINITY;
                    /* A subnormal distance below representable ray progress is
                     * an immediate collision in this material. */
                    if (sigma > 0.0 && distance == 0.0) {
                        event.kind = ALEA_NAV_COLLISION;
                        event.before = location;
                        event.after = location;
                        memcpy(event.position, position, sizeof(position));
                    } else if (!(distance > 0.0) && isfinite(distance)) {
                        err = fail_transport(failure, h, events, &location, position,
                                             particle.energy, ALEA_ERR_INVALID_STATE);
                        break;
                    } else if (alea_ray_navigator_advance(navigator, distance,
                            options->max_segment_distance, &event) != 0) {
                        err = fail_transport(failure, h, events, &location, position,
                                             particle.energy, ALEA_ERR_INVALID_STATE);
                        break;
                    }
                }
                if (!isfinite(event.distance) || event.distance < 0.0) {
                    err = fail_transport(failure, h, events, &location, position,
                                         particle.energy, ALEA_ERR_INVALID_STATE);
                    break;
                }
                if (event.distance > 0.0) {
                    double speed = particle.type == ALEA_NUC_PARTICLE_PHOTON ?
                        LIGHT_SPEED_CM_S : neutron_speed(particle.energy);
                    if (!(speed > 0.0) || !isfinite(speed)) {
                        err = fail_transport(failure, h, events, &location, position,
                                             particle.energy, ALEA_ERR_INVALID_STATE);
                        break;
                    }
                    if (tallies && alea_tally_plan_has_track_scores(
                            options->tally_plan)) {
                        err = alea_tally_record_track(tallies,
                            options->tally_plan, &location, position,
                            particle.direction, event.distance, particle.weight,
                            particle.energy, particle.time, speed,
                            particle.type == ALEA_NUC_PARTICLE_PHOTON ?
                                ALEA_TALLY_PHOTON : ALEA_TALLY_NEUTRON,
                            location.kind == ALEA_NAV_MATERIAL ?
                                &evaluation : NULL);
                        if (err != ALEA_OK) {
                            err = fail_transport(failure, h, events, &location,
                                position, particle.energy, err);
                            break;
                        }
                    }
                    if (particle.type == ALEA_NUC_PARTICLE_NEUTRON &&
                        location.cell_index >= 0 &&
                        (size_t)location.cell_index < cell_count) {
                        size_t cell = (size_t)location.cell_index;
                        history_path[cell] += particle.weight * event.distance;
                        if (!isfinite(history_path[cell])) {
                            err = fail_transport(failure, h, events, &location,
                                position, particle.energy, ALEA_ERR_OVERFLOW);
                            break;
                        }
                    }
                    particle.time += event.distance / speed;
                    if (!isfinite(particle.time)) {
                        err = fail_transport(failure, h, events, &location, position,
                                             particle.energy, ALEA_ERR_OVERFLOW);
                        break;
                    }
                    remaining_tau = fmax(0.0, remaining_tau - sigma * event.distance);
                }
                memcpy(position, event.position, sizeof(position));
                if (event.kind == ALEA_NAV_COLLISION) {
                    err = alea_nuc_evaluation_update_incident(&evaluation, &particle);
                    if (err != ALEA_OK) {
                        err = fail_transport(failure, h, events, &location,
                                             position, particle.energy, err);
                        break;
                    }
                    alea_nuc_rng_t rng;
                    err = alea_nuc_rng_init(&rng, options->seed, h,
                                            particle_ordinal, collisions,
                                            ALEA_NUC_RNG_COLLISION);
                    alea_nuc_collision_result_t collision;
                    alea_nuc_secondary_buffer_t secondaries = {
                        emitted, scratch_capacity, 0
                    };
                    if (err == ALEA_OK)
                        err = alea_nuc_collide_with_secondaries(&evaluation,
                            alea_nuc_rng_uniform, &rng, &secondaries, &collision);
                    if (err != ALEA_OK) {
                        /* The kernel reports insufficient secondary slots as
                         * OUT_OF_MEMORY; transport uses a bounded scratch. */
                        if (err == ALEA_ERR_OUT_OF_MEMORY) err = ALEA_ERR_OVERFLOW;
                        err = fail_transport(failure, h, events, &location,
                                             position, particle.energy, err);
                        break;
                    }
                    if (tallies && (alea_tally_plan_has_collision_scores(
                            options->tally_plan) ||
                        (particle.type == ALEA_NUC_PARTICLE_PHOTON &&
                         alea_tally_plan_has_deposition_scores(
                             options->tally_plan)))) {
                        int target_zaid = 0;
                        if (alea_tally_plan_needs_target(options->tally_plan)) {
                            err = alea_nuc_evaluation_component_zaid(&evaluation,
                                collision.component_index, &target_zaid);
                            if (err != ALEA_OK) {
                                err = fail_transport(failure, h, events, &location,
                                    position, particle.energy, err);
                                break;
                            }
                        }
                        uint32_t particle_mask =
                            particle.type == ALEA_NUC_PARTICLE_PHOTON ?
                                ALEA_TALLY_PHOTON : ALEA_TALLY_NEUTRON;
                        if (alea_tally_plan_has_collision_scores(
                                options->tally_plan))
                            err = alea_tally_record_collision(tallies,
                                options->tally_plan, &location, position,
                                particle.weight, particle.energy, particle.time,
                                particle_mask, collision.mt, target_zaid);
                        if (err == ALEA_OK &&
                            particle.type == ALEA_NUC_PARTICLE_PHOTON &&
                            alea_tally_plan_has_deposition_scores(
                                options->tally_plan)) {
                            err = collision.deposition_available ?
                                alea_tally_record_deposition(tallies,
                                    options->tally_plan, &location, position,
                                    particle.weight, particle.energy,
                                    particle.time, particle_mask,
                                    collision.mt, target_zaid,
                                    collision.local_energy_deposition) :
                                ALEA_ERR_UNSUPPORTED;
                        }
                        if (err != ALEA_OK) {
                            err = fail_transport(failure, h, events, &location,
                                position, particle.energy, err);
                            break;
                        }
                    }
                    result.collisions++;
                    if (particle.type == ALEA_NUC_PARTICLE_PHOTON)
                        result.photon_collisions++;
                    if (collisions == UINT32_MAX) {
                        err = fail_transport(failure, h, events, &location,
                                             position, particle.energy, ALEA_ERR_OVERFLOW);
                        break;
                    }
                    collisions++;
                    if (collision.n_emitted != secondaries.count ||
                        secondaries.count > bank_capacity - bank_count ||
                        next_ordinal + secondaries.count > UINT32_MAX) {
                        err = fail_transport(failure, h, events, &location,
                            position, particle.energy, ALEA_ERR_OVERFLOW);
                        break;
                    }
                    for (size_t i = 0; i < secondaries.count; ++i) {
                        if ((emitted[i].type != ALEA_NUC_PARTICLE_NEUTRON &&
                             emitted[i].type != ALEA_NUC_PARTICLE_PHOTON) ||
                            !isfinite(emitted[i].energy) ||
                            emitted[i].energy <= 0.0 ||
                            !isfinite(emitted[i].time) ||
                            !isfinite(emitted[i].weight) ||
                            emitted[i].weight <= 0.0 ||
                            normalize_direction(emitted[i].direction) != 0) {
                            err = fail_transport(failure, h, events, &location,
                                position, particle.energy, ALEA_ERR_UNSUPPORTED);
                            break;
                        }
                    }
                    if (err != ALEA_OK) break;
                    for (size_t i = 0; i < secondaries.count; ++i) {
                        pending_particle_t* child = &bank[bank_count++];
                        child->particle = emitted[i];
                        memcpy(child->position, position, sizeof(position));
                        child->ordinal = (uint32_t)next_ordinal++;
                        if (emitted[i].type == ALEA_NUC_PARTICLE_PHOTON)
                            result.emitted_photons++;
                        else result.emitted_neutrons++;
                    }
                    if (collision.outcome == ALEA_NUC_OUTCOME_ABSORBED) {
                        result.absorbed++;
                        if (particle.type == ALEA_NUC_PARTICLE_PHOTON)
                            result.photon_absorbed++;
                        finished = true;
                    } else if (collision.outcome == ALEA_NUC_OUTCOME_SCATTERED) {
                        particle = collision.outgoing;
                        if (!isfinite(particle.energy) || particle.energy <= 0.0 ||
                            !isfinite(particle.time) ||
                            normalize_direction(particle.direction) != 0) {
                            err = fail_transport(failure, h, events, &location,
                                position, particle.energy, ALEA_ERR_INVALID_STATE);
                            break;
                        }
                        if (alea_ray_navigator_set_direction(navigator,
                                particle.direction) != 0) {
                            err = fail_transport(failure, h, events, &location,
                                position, particle.energy, ALEA_ERR_INVALID_STATE);
                            break;
                        }
                        new_flight = true;
                    } else if (collision.outcome == ALEA_NUC_OUTCOME_REPLACED) {
                        result.replaced++;
                        if (particle.type == ALEA_NUC_PARTICLE_PHOTON)
                            result.photon_replaced++;
                        finished = true;
                    } else {
                        err = fail_transport(failure, h, events, &location,
                                             position, particle.energy,
                                             ALEA_ERR_UNSUPPORTED);
                        break;
                    }
                } else if (event.kind == ALEA_NAV_VACUUM) {
                    result.leaked++;
                    if (particle.type == ALEA_NUC_PARTICLE_PHOTON)
                        result.photon_leaked++;
                    finished = true;
                } else if (event.kind == ALEA_NAV_BOUNDARY) {
                    result.boundary_crossings++;
                    location = event.after;
                } else if (event.kind == ALEA_NAV_BOUNDARY_ACTION &&
                           event.boundary_type == ALEA_BOUNDARY_REFLECTIVE) {
                    double norm2 = event.normal[0]*event.normal[0] +
                                   event.normal[1]*event.normal[1] +
                                   event.normal[2]*event.normal[2];
                    if (!isfinite(norm2) || norm2 < 0.5 || norm2 > 1.5) {
                        err = fail_transport(failure, h, events, &location,
                                             position, particle.energy,
                                             ALEA_ERR_INVALID_STATE);
                        break;
                    }
                    double dot = particle.direction[0]*event.normal[0] +
                                 particle.direction[1]*event.normal[1] +
                                 particle.direction[2]*event.normal[2];
                    for (int j = 0; j < 3; ++j)
                        particle.direction[j] -= 2.0 * dot * event.normal[j] / norm2;
                    if (alea_ray_navigator_reflect_specular(navigator) != 0) {
                        err = fail_transport(failure, h, events, &location,
                                             position, particle.energy,
                                             ALEA_ERR_INVALID_STATE);
                        break;
                    }
                    result.reflections++;
                } else if (event.kind == ALEA_NAV_DISTANCE_LIMIT) {
                    /* Continue with the same optical-depth draw. */
                } else {
                    err = fail_transport(failure, h, events, &location,
                                         position, particle.energy,
                                         ALEA_ERR_UNSUPPORTED);
                    break;
                }
            }
            if (err != ALEA_OK && failure)
                failure->particle_ordinal = particle_ordinal;
            last_location = location;
            memcpy(last_position, position, sizeof(last_position));
            last_energy = particle.energy;
            last_ordinal = particle_ordinal;
        }
        if (err != ALEA_OK) break;
        for (size_t i = 0; i < cell_count; ++i) {
            path_sum[i] += history_path[i];
            path_sum_squared[i] += history_path[i] * history_path[i];
            if (!isfinite(path_sum[i]) || !isfinite(path_sum_squared[i])) {
                err = fail_transport(failure, h, events, &last_location,
                                     last_position, last_energy,
                                     ALEA_ERR_OVERFLOW);
                if (failure) failure->particle_ordinal = last_ordinal;
                break;
            }
        }
        if (err != ALEA_OK) break;
        if (tallies) {
            err = alea_tally_commit_history(tallies);
            if (err != ALEA_OK) {
                err = fail_transport(failure, h, events, &last_location,
                    last_position, last_energy, err);
                if (failure) failure->particle_ordinal = last_ordinal;
                break;
            }
        }
    }
    alea_ray_navigator_destroy(navigator);
    free(history_path);
    free(bank);
    free(emitted);
    if (err != ALEA_OK) {
        alea_transport_result_free(&result);
        return err;
    }
    alea_tally_results_finalize(tallies);
    *output = result;
    return ALEA_OK;
}

alea_error_t alea_transport_run_fixed_neutron(
    alea_system_t* sys, const alea_nuc_cell_bindings_t* bindings,
    const alea_transport_source_t* source,
    const alea_transport_options_t* options,
    alea_transport_result_t* output,
    alea_transport_failure_t* failure) {
    if (!output) return ALEA_ERR_NULL_ARG;
    if (!source || source->particle.type != ALEA_NUC_PARTICLE_NEUTRON) {
        memset(output, 0, sizeof(*output));
        if (failure) memset(failure, 0, sizeof(*failure));
        return source ? ALEA_ERR_INVALID_ARG : ALEA_ERR_NULL_ARG;
    }
    return alea_transport_run_fixed_source(sys, bindings, source, options,
                                           output, failure);
}
