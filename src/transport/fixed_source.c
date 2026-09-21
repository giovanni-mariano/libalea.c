// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

#include "alea_transport.h"
#include "alea.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

#define NEUTRON_REST_MEV 939.56542052
#define LIGHT_SPEED_CM_S 2.99792458e10

void alea_transport_result_free(alea_transport_result_t* result) {
    if (!result) return;
    free(result->track_length);
    free(result->track_length_squared);
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
                                      uint32_t collision_index, double* tau) {
    alea_nuc_rng_t rng;
    alea_error_t err = alea_nuc_rng_init(&rng, seed, history, 0,
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

alea_error_t alea_transport_run_fixed_neutron(
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
    if (options->histories == 0 || options->max_events_per_history == 0 ||
        !isfinite(options->max_segment_distance) ||
        options->max_segment_distance <= 0.0 ||
        source->particle.type != ALEA_NUC_PARTICLE_NEUTRON ||
        !isfinite(source->particle.energy) || source->particle.energy <= 0.0 ||
        !isfinite(source->particle.weight) || source->particle.weight <= 0.0 ||
        !isfinite(source->particle.time) ||
        !isfinite(source->position[0]) ||
        !isfinite(source->position[1]) ||
        !isfinite(source->position[2]))
        return ALEA_ERR_INVALID_ARG;

    const size_t cell_count = alea_cell_count(sys);
    if (cell_count > SIZE_MAX / sizeof(double)) return ALEA_ERR_OVERFLOW;
    double* path_sum = calloc(cell_count ? cell_count : 1, sizeof(double));
    double* path_sum_squared = calloc(cell_count ? cell_count : 1, sizeof(double));
    double* history_path = calloc(cell_count ? cell_count : 1, sizeof(double));
    alea_ray_navigator_t* navigator = alea_ray_navigator_create(sys);
    if (!path_sum || !path_sum_squared || !history_path || !navigator) {
        free(path_sum); free(path_sum_squared); free(history_path);
        alea_ray_navigator_destroy(navigator);
        return ALEA_ERR_OUT_OF_MEMORY;
    }
    alea_transport_result_t result = {
        .cell_count = cell_count, .histories = options->histories,
        .track_length = path_sum, .track_length_squared = path_sum_squared
    };
    alea_error_t err = ALEA_OK;
    for (uint32_t h = 0; h < options->histories; ++h) {
        memset(history_path, 0, cell_count * sizeof(double));
        alea_nuc_particle_state_t particle = source->particle;
        if (normalize_direction(particle.direction) != 0) {
            err = ALEA_ERR_INVALID_ARG;
            break;
        }
        double position[3];
        memcpy(position, source->position, sizeof(position));
        alea_nav_location_t location;
        if (alea_ray_navigator_restart(navigator, position,
                                       particle.direction, &location) != 0) {
            err = fail_transport(failure, h, 0, NULL, position,
                                 particle.energy, ALEA_ERR_INVALID_STATE);
            break;
        }
        double remaining_tau = 0.0;
        bool new_flight = true;
        bool finished = false;
        uint32_t collisions = 0;
        uint32_t events = 0;
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
                err = draw_optical_depth(options->seed, h, collisions,
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
                            ALEA_NUC_PARTICLE_NEUTRON);
                    err = prepared ? alea_nuc_evaluate(prepared, &particle,
                        &evaluation) : ALEA_ERR_INVALID_STATE;
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
                if (location.cell_index >= 0 &&
                    (size_t)location.cell_index < cell_count) {
                    size_t cell = (size_t)location.cell_index;
                    history_path[cell] += particle.weight * event.distance;
                    if (!isfinite(history_path[cell])) {
                        err = fail_transport(failure, h, events, &location,
                            position, particle.energy, ALEA_ERR_OVERFLOW);
                        break;
                    }
                }
                double speed = neutron_speed(particle.energy);
                if (!(speed > 0.0) || !isfinite(speed)) {
                    err = fail_transport(failure, h, events, &location, position,
                                         particle.energy, ALEA_ERR_INVALID_STATE);
                    break;
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
                err = alea_nuc_rng_init(&rng, options->seed, h, 0, collisions,
                                        ALEA_NUC_RNG_COLLISION);
                alea_nuc_collision_result_t collision;
                if (err == ALEA_OK)
                    err = alea_nuc_collide(&evaluation, alea_nuc_rng_uniform,
                                            &rng, &collision);
                if (err != ALEA_OK) {
                    err = fail_transport(failure, h, events, &location,
                                         position, particle.energy, err);
                    break;
                }
                result.collisions++;
                if (collisions == UINT32_MAX) {
                    err = fail_transport(failure, h, events, &location,
                                         position, particle.energy, ALEA_ERR_OVERFLOW);
                    break;
                }
                collisions++;
                if (collision.outcome == ALEA_NUC_OUTCOME_ABSORBED) {
                    result.absorbed++;
                    finished = true;
                } else if (collision.outcome == ALEA_NUC_OUTCOME_SCATTERED) {
                    particle = collision.outgoing;
                    if (alea_ray_navigator_set_direction(navigator,
                            particle.direction) != 0) {
                        err = fail_transport(failure, h, events, &location,
                            position, particle.energy, ALEA_ERR_INVALID_STATE);
                        break;
                    }
                    new_flight = true;
                } else {
                    err = fail_transport(failure, h, events, &location,
                                         position, particle.energy,
                                         ALEA_ERR_UNSUPPORTED);
                    break;
                }
            } else if (event.kind == ALEA_NAV_VACUUM) {
                result.leaked++;
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
        if (err != ALEA_OK) break;
        for (size_t i = 0; i < cell_count; ++i) {
            path_sum[i] += history_path[i];
            path_sum_squared[i] += history_path[i] * history_path[i];
            if (!isfinite(path_sum[i]) || !isfinite(path_sum_squared[i])) {
                err = fail_transport(failure, h, events, &location, position,
                                     particle.energy, ALEA_ERR_OVERFLOW);
                break;
            }
        }
        if (err != ALEA_OK) break;
    }
    alea_ray_navigator_destroy(navigator);
    free(history_path);
    if (err != ALEA_OK) {
        alea_transport_result_free(&result);
        return err;
    }
    *output = result;
    return ALEA_OK;
}
