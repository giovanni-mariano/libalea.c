// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

#include "alea_source.h"
#include "../rng/alea_rng.h"
#include "../rng/alea_rng_distribution.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

struct alea_source {
    alea_source_spec_t spec;
    double axis_unit[3];
    double radial_u[3];
    double radial_v[3];
    double axis_length;
    double angle_u[3];
    double angle_v[3];
    double* energy_values;
    alea_rng_discrete_table_t* energy_table;
    alea_rng_tabular_table_t* angle_table;
};

static double length3(const double v[3]) {
    double scale = fmax(fabs(v[0]), fmax(fabs(v[1]), fabs(v[2])));
    if (!isfinite(scale) || scale <= 0.0) return 0.0;
    return scale * sqrt((v[0]/scale)*(v[0]/scale) +
                        (v[1]/scale)*(v[1]/scale) +
                        (v[2]/scale)*(v[2]/scale));
}

static int normalize(double direction[3]) {
    double scale = fmax(fabs(direction[0]),
                        fmax(fabs(direction[1]), fabs(direction[2])));
    if (!isfinite(scale) || scale <= 0.0) return 0;
    double x = direction[0] / scale;
    double y = direction[1] / scale;
    double z = direction[2] / scale;
    double length = sqrt(x*x + y*y + z*z);
    if (!isfinite(length) || length <= 0.0) return 0;
    direction[0] = x / length;
    direction[1] = y / length;
    direction[2] = z / length;
    return 1;
}

static int angle_frame(alea_source_t* source) {
    double* n = source->spec.direction;
    double reference[3] = {0, 0, 1};
    if (fabs(n[2]) > 0.9) {
        reference[0] = 1; reference[2] = 0;
    }
    source->angle_u[0] = reference[1]*n[2] - reference[2]*n[1];
    source->angle_u[1] = reference[2]*n[0] - reference[0]*n[2];
    source->angle_u[2] = reference[0]*n[1] - reference[1]*n[0];
    if (!normalize(source->angle_u)) return 0;
    const double* u = source->angle_u;
    source->angle_v[0] = n[1]*u[2] - n[2]*u[1];
    source->angle_v[1] = n[2]*u[0] - n[0]*u[2];
    source->angle_v[2] = n[0]*u[1] - n[1]*u[0];
    return 1;
}

alea_error_t alea_source_prepare(const alea_source_spec_t* spec,
                                 alea_source_t** output) {
    if (!spec || !output) return ALEA_ERR_NULL_ARG;
    *output = NULL;
    if ((spec->particle != ALEA_NUC_PARTICLE_NEUTRON &&
         spec->particle != ALEA_NUC_PARTICLE_PHOTON) ||
        (spec->space < ALEA_SOURCE_POINT || spec->space > ALEA_SOURCE_CYLINDER) ||
        (spec->angle < ALEA_SOURCE_MONODIRECTIONAL ||
         spec->angle > ALEA_SOURCE_RADIAL) ||
        (spec->energy_type != ALEA_SOURCE_ENERGY_MONO &&
         spec->energy_type != ALEA_SOURCE_ENERGY_LINES) ||
        (spec->energy_type == ALEA_SOURCE_ENERGY_MONO &&
         (!isfinite(spec->energy) || spec->energy <= 0.0)) ||
        !isfinite(spec->time) || !isfinite(spec->weight) ||
        spec->weight <= 0.0) return ALEA_ERR_INVALID_ARG;
    if (spec->energy_type == ALEA_SOURCE_ENERGY_LINES) {
        if (!spec->energy_values || !spec->energy_weights ||
            spec->energy_count == 0 || spec->energy_count > UINT32_MAX ||
            spec->energy_count > SIZE_MAX / sizeof(double)) return ALEA_ERR_INVALID_ARG;
        for (size_t i = 0; i < spec->energy_count; ++i)
            if (!isfinite(spec->energy_values[i]) || spec->energy_values[i] <= 0.0)
                return ALEA_ERR_INVALID_ARG;
    }
    for (int i = 0; i < 3; ++i) {
        if (spec->space == ALEA_SOURCE_POINT && !isfinite(spec->position[i]))
            return ALEA_ERR_INVALID_ARG;
        if (spec->space == ALEA_SOURCE_BOX &&
            (!isfinite(spec->lower[i]) || !isfinite(spec->upper[i]) ||
             spec->upper[i] < spec->lower[i])) return ALEA_ERR_INVALID_ARG;
        if (spec->space == ALEA_SOURCE_LINE &&
            (!isfinite(spec->start[i]) || !isfinite(spec->end[i])))
            return ALEA_ERR_INVALID_ARG;
        if (spec->space == ALEA_SOURCE_SPHERE && !isfinite(spec->center[i]))
            return ALEA_ERR_INVALID_ARG;
        if (spec->space == ALEA_SOURCE_CYLINDER &&
            (!isfinite(spec->base[i]) || !isfinite(spec->axis[i])))
            return ALEA_ERR_INVALID_ARG;
    }
    if (spec->space == ALEA_SOURCE_SPHERE ||
        spec->space == ALEA_SOURCE_CYLINDER) {
        if (!isfinite(spec->inner_radius) || !isfinite(spec->outer_radius) ||
            spec->inner_radius < 0.0 ||
            spec->outer_radius <= spec->inner_radius ||
            spec->outer_radius > 1e100) return ALEA_ERR_INVALID_ARG;
    }
    if (spec->angle == ALEA_SOURCE_CONE &&
        (!isfinite(spec->cone_half_angle) || spec->cone_half_angle < 0.0 ||
         spec->cone_half_angle > 3.14159265358979323846))
        return ALEA_ERR_INVALID_ARG;
    if (spec->angle == ALEA_SOURCE_RADIAL) {
        for (int i = 0; i < 3; ++i)
            if (!isfinite(spec->angle_origin[i])) return ALEA_ERR_INVALID_ARG;
        if (spec->space == ALEA_SOURCE_POINT &&
            spec->position[0] == spec->angle_origin[0] &&
            spec->position[1] == spec->angle_origin[1] &&
            spec->position[2] == spec->angle_origin[2])
            return ALEA_ERR_INVALID_ARG;
    }
    if (spec->angle == ALEA_SOURCE_TABULATED_MU) {
        if (!spec->angle_mu || !spec->angle_pdf || spec->angle_count < 2 ||
            spec->angle_count > UINT32_MAX ||
            spec->angle_count > SIZE_MAX / sizeof(double) ||
            (spec->angle_interpolation != ALEA_SOURCE_PDF_HISTOGRAM &&
             spec->angle_interpolation != ALEA_SOURCE_PDF_LINEAR))
            return ALEA_ERR_INVALID_ARG;
        for (size_t i = 0; i < spec->angle_count; ++i)
            if (!isfinite(spec->angle_mu[i]) || spec->angle_mu[i] < -1.0 ||
                spec->angle_mu[i] > 1.0 ||
                (i && spec->angle_mu[i] <= spec->angle_mu[i-1]))
                return ALEA_ERR_INVALID_ARG;
    }
    alea_source_t* source = calloc(1, sizeof(*source));
    if (!source) return ALEA_ERR_OUT_OF_MEMORY;
    source->spec = *spec;
    if (spec->angle == ALEA_SOURCE_TABULATED_MU) {
        source->angle_table = alea_rng_tabular_table_create(
            spec->angle_mu, spec->angle_pdf, spec->angle_count,
            spec->angle_interpolation == ALEA_SOURCE_PDF_HISTOGRAM ?
                ALEA_RNG_TABULAR_HISTOGRAM : ALEA_RNG_TABULAR_LIN_LIN);
        if (!source->angle_table) {
            alea_error_t err = (alea_error_t)alea_error_code();
            alea_source_free(source);
            return err == ALEA_OK ? ALEA_ERR_INVALID_ARG : err;
        }
        source->spec.angle_mu = NULL;
        source->spec.angle_pdf = NULL;
    }
    if (spec->energy_type == ALEA_SOURCE_ENERGY_LINES) {
        source->energy_values = malloc(spec->energy_count * sizeof(double));
        if (!source->energy_values) { alea_source_free(source); return ALEA_ERR_OUT_OF_MEMORY; }
        memcpy(source->energy_values, spec->energy_values,
               spec->energy_count * sizeof(double));
        source->energy_table = alea_rng_discrete_table_create(
            spec->energy_weights, spec->energy_count);
        if (!source->energy_table) {
            alea_error_t err = (alea_error_t)alea_error_code();
            alea_source_free(source);
            return err == ALEA_OK ? ALEA_ERR_INVALID_ARG : err;
        }
        source->spec.energy_values = source->energy_values;
        source->spec.energy_weights = NULL;
    }
    if (spec->space == ALEA_SOURCE_LINE) {
        double displacement[3];
        for (int i = 0; i < 3; ++i)
            displacement[i] = spec->end[i] - spec->start[i];
        if (!isfinite(length3(displacement)) || length3(displacement) <= 0.0) {
            alea_source_free(source); return ALEA_ERR_INVALID_ARG;
        }
    }
    if (spec->space == ALEA_SOURCE_CYLINDER) {
        source->axis_length = length3(spec->axis);
        if (!isfinite(source->axis_length) || source->axis_length <= 0.0) {
            alea_source_free(source); return ALEA_ERR_INVALID_ARG;
        }
        for (int i = 0; i < 3; ++i)
            source->axis_unit[i] = spec->axis[i] / source->axis_length;
        double* n = source->axis_unit;
        double reference[3] = {0, 0, 1};
        if (fabs(n[2]) > 0.9) {
            reference[0] = 1; reference[2] = 0;
        }
        source->radial_u[0] = reference[1]*n[2] - reference[2]*n[1];
        source->radial_u[1] = reference[2]*n[0] - reference[0]*n[2];
        source->radial_u[2] = reference[0]*n[1] - reference[1]*n[0];
        if (!normalize(source->radial_u)) {
            alea_source_free(source); return ALEA_ERR_INVALID_ARG;
        }
        double* u = source->radial_u;
        source->radial_v[0] = n[1]*u[2] - n[2]*u[1];
        source->radial_v[1] = n[2]*u[0] - n[0]*u[2];
        source->radial_v[2] = n[0]*u[1] - n[1]*u[0];
    }
    if (spec->angle == ALEA_SOURCE_MONODIRECTIONAL ||
        spec->angle == ALEA_SOURCE_CONE || spec->angle == ALEA_SOURCE_COSINE ||
        spec->angle == ALEA_SOURCE_TABULATED_MU) {
        if (!normalize(source->spec.direction) ||
            (spec->angle != ALEA_SOURCE_MONODIRECTIONAL &&
             !angle_frame(source))) {
            alea_source_free(source);
            return ALEA_ERR_INVALID_ARG;
        }
    }
    *output = source;
    return ALEA_OK;
}

void alea_source_free(alea_source_t* source) {
    if (!source) return;
    alea_rng_discrete_table_destroy(source->energy_table);
    alea_rng_tabular_table_destroy(source->angle_table);
    free(source->energy_values);
    free(source);
}

uint32_t alea_source_particle_mask(const alea_source_t* source) {
    if (!source) return 0;
    return source->spec.particle == ALEA_NUC_PARTICLE_NEUTRON ?
        ALEA_NUC_BIND_NEUTRON : ALEA_NUC_BIND_PHOTON;
}

alea_error_t alea_source_sample(void* context, uint64_t seed,
                                uint32_t history_id,
                                alea_transport_source_t* output) {
    if (!context || !output) return ALEA_ERR_NULL_ARG;
    const alea_source_t* source = context;
    const alea_source_spec_t* spec = &source->spec;
    alea_transport_source_t sample = {0};
    sample.particle.type = spec->particle;
    sample.particle.energy = spec->energy;
    if (spec->energy_type == ALEA_SOURCE_ENERGY_LINES) {
        alea_rng_event_t event;
        size_t index;
        if (alea_rng_event_init(&event, ALEA_RNG_PHILOX4X32_10, seed,
            ALEA_RNG_DOMAIN_TRANSPORT_SOURCE_ENERGY, history_id, 0) != 0 ||
            alea_rng_sample_discrete_alias(&event, source->energy_table,
                                           &index) != 0) return ALEA_ERR_INVALID_STATE;
        sample.particle.energy = source->energy_values[index];
    }
    sample.particle.weight = spec->weight;
    sample.particle.time = spec->time;
    if (spec->space == ALEA_SOURCE_POINT) {
        memcpy(sample.position, spec->position, sizeof(sample.position));
    } else if (spec->space == ALEA_SOURCE_BOX) {
        for (int i = 0; i < 3; ++i) {
            double u;
            if (alea_rng_uniform53_at(ALEA_RNG_PHILOX4X32_10, seed,
                ALEA_RNG_DOMAIN_TRANSPORT_SOURCE_POSITION, history_id,
                (uint64_t)i, &u) != 0) return ALEA_ERR_INVALID_STATE;
            sample.position[i] = spec->lower[i] +
                u * (spec->upper[i] - spec->lower[i]);
            if (!isfinite(sample.position[i])) return ALEA_ERR_OVERFLOW;
        }
    } else {
        double u[3] = {0};
        int draws = spec->space == ALEA_SOURCE_LINE ? 1 : 3;
        for (int i = 0; i < draws; ++i)
            if (alea_rng_uniform53_at(ALEA_RNG_PHILOX4X32_10, seed,
                ALEA_RNG_DOMAIN_TRANSPORT_SOURCE_POSITION, history_id,
                (uint64_t)i, &u[i]) != 0) return ALEA_ERR_INVALID_STATE;
        if (spec->space == ALEA_SOURCE_LINE) {
            for (int i = 0; i < 3; ++i)
                sample.position[i] = spec->start[i] +
                    u[0]*(spec->end[i] - spec->start[i]);
        } else if (spec->space == ALEA_SOURCE_SPHERE) {
            double inner3 = spec->inner_radius*spec->inner_radius*spec->inner_radius;
            double outer3 = spec->outer_radius*spec->outer_radius*spec->outer_radius;
            double radius = cbrt(inner3 + u[0]*(outer3 - inner3));
            double mu = 2*u[1] - 1;
            double phi = 6.2831853071795864769*u[2];
            double transverse = sqrt(fmax(0.0, 1 - mu*mu));
            sample.position[0] = spec->center[0] + radius*transverse*cos(phi);
            sample.position[1] = spec->center[1] + radius*transverse*sin(phi);
            sample.position[2] = spec->center[2] + radius*mu;
        } else {
            double inner2 = spec->inner_radius*spec->inner_radius;
            double outer2 = spec->outer_radius*spec->outer_radius;
            double radius = sqrt(inner2 + u[0]*(outer2 - inner2));
            double phi = 6.2831853071795864769*u[1];
            double x = radius*cos(phi), y = radius*sin(phi);
            for (int i = 0; i < 3; ++i)
                sample.position[i] = spec->base[i] +
                    u[2]*spec->axis[i] +
                    x*source->radial_u[i] + y*source->radial_v[i];
        }
        for (int i = 0; i < 3; ++i)
            if (!isfinite(sample.position[i])) return ALEA_ERR_OVERFLOW;
    }
    if (spec->angle == ALEA_SOURCE_MONODIRECTIONAL) {
        memcpy(sample.particle.direction, spec->direction,
               sizeof(sample.particle.direction));
    } else if (spec->angle == ALEA_SOURCE_RADIAL) {
        for (int i = 0; i < 3; ++i)
            sample.particle.direction[i] =
                (spec->radial_inward ? -1.0 : 1.0) *
                (sample.position[i] - spec->angle_origin[i]);
        if (!normalize(sample.particle.direction)) return ALEA_ERR_INVALID_STATE;
    } else {
        alea_rng_event_t event;
        if (alea_rng_event_init(&event, ALEA_RNG_PHILOX4X32_10, seed,
            ALEA_RNG_DOMAIN_TRANSPORT_SOURCE_DIRECTION, history_id, 0) != 0)
            return ALEA_ERR_INVALID_STATE;
        double u = 0, v, mu;
        if (spec->angle == ALEA_SOURCE_TABULATED_MU) {
            if (alea_rng_sample_tabular(&event, source->angle_table, &mu) != 0)
                return ALEA_ERR_INVALID_STATE;
        } else {
            if (alea_rng_event_next_uniform(&event, &u) != 0)
                return ALEA_ERR_INVALID_STATE;
            mu = 2.0*u - 1.0;
        }
        if (alea_rng_event_next_uniform(&event, &v) != 0)
            return ALEA_ERR_INVALID_STATE;
        if (spec->angle == ALEA_SOURCE_CONE)
            mu = 1.0 - u*(1.0 - cos(spec->cone_half_angle));
        else if (spec->angle == ALEA_SOURCE_COSINE)
            mu = sqrt(u);
        double phi = 6.2831853071795864769*v;
        double transverse = sqrt(fmax(0.0, 1.0 - mu*mu));
        if (spec->angle == ALEA_SOURCE_ISOTROPIC) {
            sample.particle.direction[0] = transverse*cos(phi);
            sample.particle.direction[1] = transverse*sin(phi);
            sample.particle.direction[2] = mu;
        } else {
            for (int i = 0; i < 3; ++i)
                sample.particle.direction[i] = mu*spec->direction[i] +
                    transverse*(cos(phi)*source->angle_u[i] +
                                sin(phi)*source->angle_v[i]);
        }
    }
    *output = sample;
    return ALEA_OK;
}

alea_error_t alea_source_sample_batch(const alea_source_t* source,
    uint64_t seed, uint32_t history_offset, uint32_t count,
    alea_transport_source_t* output) {
    if (!source || (count && !output)) return ALEA_ERR_NULL_ARG;
    if (count && count - 1 > UINT32_MAX - history_offset)
        return ALEA_ERR_OVERFLOW;
    for (uint32_t i = 0; i < count; ++i) {
        alea_error_t err = alea_source_sample((void*)source, seed,
            history_offset + i, &output[i]);
        if (err != ALEA_OK) return err;
    }
    return ALEA_OK;
}
