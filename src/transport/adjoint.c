// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

#include "alea_adjoint.h"
#include "../rng/alea_rng.h"
#include <math.h>
#include <string.h>

#define ALEA_ADJOINT_TWO_PI 6.28318530717958647693
#define ALEA_ADJOINT_FOUR_PI 12.56637061435917295385

static int direction_normalize(double v[3]) {
    double m = fmax(fabs(v[0]), fmax(fabs(v[1]), fabs(v[2])));
    if (!(m > 0.0) || !isfinite(m)) return -1;
    double x = v[0]/m, y = v[1]/m, z = v[2]/m;
    double n = sqrt(x*x + y*y + z*z);
    if (!(n > 0.0) || !isfinite(n)) return -1;
    v[0] = x/n; v[1] = y/n; v[2] = z/n;
    return 0;
}

static alea_error_t draw(uint64_t seed, uint32_t history, uint32_t event,
                         alea_rng_domain_t domain, uint32_t index, double* u) {
    uint64_t address = ((uint64_t)event << 32) | index;
    if (alea_rng_uniform53_at(ALEA_RNG_PHILOX4X32_10, seed, domain,
            history, address, u) != 0 || !isfinite(*u) ||
            *u < 0.0 || *u >= 1.0) return ALEA_ERR_INVALID_STATE;
    return ALEA_OK;
}

alea_error_t alea_adjoint_sample_box_detector(
    void* context, uint64_t seed, uint32_t history_id,
    alea_adjoint_detector_particle_t* output) {
    if (!context || !output) return ALEA_ERR_NULL_ARG;
    const alea_adjoint_box_detector_t* box = context;
    double volume = 1.0, u[5];
    for (int i = 0; i < 3; ++i) {
        double width = box->upper[i] - box->lower[i];
        if (!isfinite(box->lower[i]) || !isfinite(box->upper[i]) ||
            !(width > 0.0)) return ALEA_ERR_INVALID_ARG;
        volume *= width;
    }
    if (!isfinite(box->response) || !(box->response > 0.0) ||
        !isfinite(volume * box->response * ALEA_ADJOINT_FOUR_PI))
        return ALEA_ERR_INVALID_ARG;
    for (uint32_t i = 0; i < 5; ++i) {
        alea_error_t err = draw(seed, history_id, 0,
            ALEA_RNG_DOMAIN_TRANSPORT_SOURCE_POSITION, i, &u[i]);
        if (err != ALEA_OK) return err;
    }
    memset(output, 0, sizeof(*output));
    for (int i = 0; i < 3; ++i)
        output->position[i] = box->lower[i] +
            (box->upper[i] - box->lower[i]) * u[i];
    double mu = 2.0*u[3] - 1.0;
    double phi = ALEA_ADJOINT_TWO_PI*u[4];
    double radial = sqrt(fmax(0.0, 1.0 - mu*mu));
    output->physical_direction[0] = radial*cos(phi);
    output->physical_direction[1] = radial*sin(phi);
    output->physical_direction[2] = mu;
    output->group = box->group;
    output->weight = volume * box->response * ALEA_ADJOINT_FOUR_PI;
    return ALEA_OK;
}

static alea_error_t fail(alea_transport_failure_t* failure, uint32_t history,
    uint32_t event, const alea_nav_location_t* location, const double pos[3],
    alea_error_t err) {
    if (failure) {
        failure->history_id = history;
        failure->event_index = event;
        failure->cell_id = location ? location->cell_id : -1;
        failure->error = err;
        memcpy(failure->position, pos, 3*sizeof(double));
    }
    alea_set_error_detail(err, "adjoint history %u event %u cell %d",
        history, event, location ? location->cell_id : -1);
    return err;
}

static alea_error_t validate(const alea_system_t* sys,
    const alea_adjoint_problem_t* p, const alea_adjoint_options_t* o) {
    if (!sys || !p || !o || !p->cell_materials || !p->physical_source ||
        !p->detector_sampler) return ALEA_ERR_NULL_ARG;
    if (p->particle != ALEA_NUC_PARTICLE_NEUTRON &&
        p->particle != ALEA_NUC_PARTICLE_PHOTON) return ALEA_ERR_INVALID_ARG;
    if (!p->n_groups || p->n_groups > SIZE_MAX/p->n_groups ||
        p->cell_count != alea_cell_count(sys) ||
        p->cell_count > SIZE_MAX/p->n_groups || !o->histories ||
        !o->max_events_per_history ||
        !(o->max_segment_distance > 0.0) ||
        !isfinite(o->max_segment_distance) ||
        (o->navigation_validation != ALEA_NAV_VALIDATE_STRICT &&
         o->navigation_validation != ALEA_NAV_VALIDATE_FAST &&
         o->navigation_validation != ALEA_NAV_VALIDATE_INTERVAL))
        return ALEA_ERR_INVALID_ARG;
    if (o->histories-1 > UINT32_MAX-o->history_offset)
        return ALEA_ERR_OVERFLOW;
    for (size_t cell = 0; cell < p->cell_count; ++cell) {
        const alea_adjoint_material_t* m = &p->cell_materials[cell];
        if ((m->total == NULL) != (m->transfer == NULL))
            return ALEA_ERR_INVALID_ARG;
        for (size_t g = 0; g < p->n_groups; ++g) {
            double q = p->physical_source[cell*p->n_groups + g];
            if (!isfinite(q) || q < 0.0) return ALEA_ERR_INVALID_ARG;
            if (!m->total) continue;
            double sigma = m->total[g];
            if (!isfinite(sigma) || sigma < 0.0) return ALEA_ERR_INVALID_ARG;
            for (size_t h = 0; h < p->n_groups; ++h) {
                double k = m->transfer[g*p->n_groups + h];
                if (!isfinite(k) || k < 0.0 || (sigma == 0.0 && k > 0.0))
                    return ALEA_ERR_INVALID_ARG;
            }
        }
        /* The collision estimator uses Sigma_t[g] as its event rate. A
         * transposed row with production but no event rate needs a different
         * tracking scheme and must not silently produce a zero response. */
        if (m->total) for (size_t g = 0; g < p->n_groups; ++g) {
            if (m->total[g] != 0.0) continue;
            for (size_t h = 0; h < p->n_groups; ++h)
                if (m->transfer[h*p->n_groups + g] > 0.0)
                    return ALEA_ERR_UNSUPPORTED;
        }
    }
    return ALEA_OK;
}

alea_error_t alea_adjoint_run(
    alea_system_t* sys, const alea_adjoint_problem_t* p,
    const alea_adjoint_options_t* o, alea_adjoint_result_t* output,
    alea_transport_failure_t* failure) {
    if (!output) return ALEA_ERR_NULL_ARG;
    memset(output, 0, sizeof(*output));
    if (failure) memset(failure, 0, sizeof(*failure));
    alea_error_t err = validate(sys, p, o);
    if (err != ALEA_OK) return err;
    alea_ray_navigator_t* nav = alea_ray_navigator_create(sys);
    if (!nav) return ALEA_ERR_OUT_OF_MEMORY;
    alea_ray_navigator_set_validation_mode(nav, o->navigation_validation);
    if (o->max_navigation_breakpoints)
        alea_ray_navigator_set_interval_budget(nav, o->max_navigation_breakpoints);

    double sum = 0.0, sum_squared = 0.0;
    for (uint32_t local = 0; local < o->histories; ++local) {
        uint32_t history = o->history_offset + local;
        alea_adjoint_detector_particle_t detector = {0};
        err = p->detector_sampler(p->detector_context, o->seed, history,
                                  &detector);
        double position[3] = {detector.position[0], detector.position[1],
                              detector.position[2]};
        double direction[3] = {-detector.physical_direction[0],
                               -detector.physical_direction[1],
                               -detector.physical_direction[2]};
        if (err != ALEA_OK || detector.group >= p->n_groups ||
            !isfinite(detector.weight) || !(detector.weight > 0.0) ||
            !isfinite(position[0]) || !isfinite(position[1]) ||
            !isfinite(position[2]) || direction_normalize(direction) != 0) {
            err = fail(failure, history, 0, NULL, position,
                       err == ALEA_OK ? ALEA_ERR_INVALID_ARG : err);
            break;
        }
        alea_nav_location_t location;
        if (alea_ray_navigator_restart(nav, position, direction,
                                        &location) != 0) {
            err = fail(failure, history, 0, NULL, position,
                       ALEA_ERR_INVALID_STATE);
            break;
        }
        size_t group = detector.group;
        double weight = detector.weight, score = 0.0, tau = 0.0;
        int new_flight = 1, done = 0;
        uint32_t collisions = 0, events = 0;
        while (!done) {
            if (events == o->max_events_per_history) {
                err = fail(failure, history, events, &location, position,
                           ALEA_ERR_OVERFLOW);
                break;
            }
            ++events;
            if (location.kind != ALEA_NAV_MATERIAL &&
                location.kind != ALEA_NAV_VOID) {
                err = fail(failure, history, events, &location, position,
                           ALEA_ERR_INVALID_STATE);
                break;
            }
            if (new_flight) {
                double u;
                err = draw(o->seed, history, collisions,
                    ALEA_RNG_DOMAIN_TRANSPORT_FREE_PATH, 0, &u);
                if (err != ALEA_OK) break;
                tau = -log1p(-u);
                new_flight = 0;
            }
            const alea_adjoint_material_t* m = NULL;
            double sigma = 0.0;
            if (location.cell_index >= 0 &&
                (size_t)location.cell_index < p->cell_count) {
                m = &p->cell_materials[location.cell_index];
                if (location.kind == ALEA_NAV_MATERIAL && !m->total) {
                    err = fail(failure, history, events, &location, position,
                               ALEA_ERR_INVALID_STATE);
                    break;
                }
                if (location.kind == ALEA_NAV_MATERIAL) sigma = m->total[group];
            } else if (location.kind == ALEA_NAV_MATERIAL) {
                err = fail(failure, history, events, &location, position,
                           ALEA_ERR_INVALID_STATE);
                break;
            }
            double distance = sigma > 0.0 ? tau/sigma : INFINITY;
            alea_nav_event_t event = {0};
            if (sigma > 0.0 && distance == 0.0) {
                event.kind = ALEA_NAV_COLLISION;
                event.before = event.after = location;
                memcpy(event.position, position, sizeof(position));
            } else if (alea_ray_navigator_advance(nav, distance,
                           o->max_segment_distance, &event) != 0) {
                err = fail(failure, history, events, &location, position,
                           ALEA_ERR_INVALID_STATE);
                break;
            }
            if (!(event.distance >= 0.0) || !isfinite(event.distance)) {
                err = fail(failure, history, events, &location, position,
                           ALEA_ERR_INVALID_STATE);
                break;
            }
            if (location.cell_index >= 0 &&
                (size_t)location.cell_index < p->cell_count) {
                double q = p->physical_source[
                    (size_t)location.cell_index*p->n_groups + group];
                score += weight * q / ALEA_ADJOINT_FOUR_PI * event.distance;
            }
            tau = fmax(0.0, tau - sigma*event.distance);
            memcpy(position, event.position, sizeof(position));
            if (!isfinite(score) || !isfinite(tau)) {
                err = fail(failure, history, events, &location, position,
                           ALEA_ERR_OVERFLOW);
                break;
            }
            if (event.kind == ALEA_NAV_COLLISION) {
                double row = 0.0;
                for (size_t h = 0; h < p->n_groups; ++h)
                    row += m->transfer[h*p->n_groups + group];
                if (!isfinite(row)) {
                    err = fail(failure, history, events, &location, position,
                               ALEA_ERR_OVERFLOW);
                    break;
                }
                if (row == 0.0) { done = 1; continue; }
                double u, mu, az;
                err = draw(o->seed, history, collisions,
                    ALEA_RNG_DOMAIN_TRANSPORT_REACTION, 0, &u);
                if (err == ALEA_OK) err = draw(o->seed, history, collisions,
                    ALEA_RNG_DOMAIN_TRANSPORT_SCATTER_ANGLE, 0, &mu);
                if (err == ALEA_OK) err = draw(o->seed, history, collisions,
                    ALEA_RNG_DOMAIN_TRANSPORT_SCATTER_ANGLE, 1, &az);
                if (err != ALEA_OK) break;
                double target = u*row, running = 0.0;
                size_t selected = group;
                for (size_t h = 0; h < p->n_groups; ++h) {
                    double k = m->transfer[h*p->n_groups + group];
                    if (k <= 0.0) continue;
                    running += k;
                    selected = h;
                    if (target < running) break;
                }
                weight *= row/sigma;
                if (!isfinite(weight) || !(weight > 0.0)) {
                    err = fail(failure, history, events, &location, position,
                               ALEA_ERR_OVERFLOW);
                    break;
                }
                group = selected;
                double z = 2.0*mu - 1.0;
                double r = sqrt(fmax(0.0, 1.0-z*z));
                direction[0] = r*cos(ALEA_ADJOINT_TWO_PI*az);
                direction[1] = r*sin(ALEA_ADJOINT_TWO_PI*az);
                direction[2] = z;
                if (alea_ray_navigator_set_direction(nav, direction) != 0) {
                    err = fail(failure, history, events, &location, position,
                               ALEA_ERR_INVALID_STATE);
                    break;
                }
                if (collisions == UINT32_MAX) {
                    err = fail(failure, history, events, &location, position,
                               ALEA_ERR_OVERFLOW);
                    break;
                }
                ++collisions;
                new_flight = 1;
            } else if (event.kind == ALEA_NAV_BOUNDARY) {
                location = event.after;
            } else if (event.kind == ALEA_NAV_VACUUM) {
                done = 1;
            } else if (event.kind == ALEA_NAV_DISTANCE_LIMIT) {
                /* Retain optical depth across artificial segments. */
            } else if (event.kind == ALEA_NAV_BOUNDARY_ACTION &&
                       event.boundary_type == ALEA_BOUNDARY_REFLECTIVE) {
                double nn = 0.0, dot = 0.0;
                for (int j = 0; j < 3; ++j) {
                    nn += event.normal[j]*event.normal[j];
                    dot += direction[j]*event.normal[j];
                }
                if (!(nn > 0.5 && nn < 1.5)) {
                    err = fail(failure, history, events, &location, position,
                               ALEA_ERR_INVALID_STATE);
                    break;
                }
                for (int j = 0; j < 3; ++j)
                    direction[j] -= 2.0*dot*event.normal[j]/nn;
                if (alea_ray_navigator_reflect_specular(nav) != 0) {
                    err = fail(failure, history, events, &location, position,
                               ALEA_ERR_INVALID_STATE);
                    break;
                }
            } else {
                err = fail(failure, history, events, &location, position,
                           ALEA_ERR_UNSUPPORTED);
                break;
            }
        }
        if (err != ALEA_OK) break;
        sum += score;
        sum_squared += score*score;
        if (!isfinite(sum) || !isfinite(sum_squared)) {
            err = fail(failure, history, events, &location, position,
                       ALEA_ERR_OVERFLOW);
            break;
        }
    }
    alea_ray_navigator_destroy(nav);
    if (err != ALEA_OK) return err;
    output->histories = o->histories;
    output->sum = sum;
    output->sum_squared = sum_squared;
    output->mean = sum / o->histories;
    if (o->histories > 1) {
        double centered = sum_squared - sum*output->mean;
        if (centered < 0.0) centered = 0.0;
        output->standard_error = sqrt(centered /
            (o->histories * (double)(o->histories-1)));
    }
    return ALEA_OK;
}
