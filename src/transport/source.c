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
    int is_mixture;
    alea_source_t** components;
    size_t component_count;
    uint32_t component_mask;
    alea_rng_discrete_table_t* component_table;
    alea_source_spec_t spec;
    double axis_unit[3];
    double radial_u[3];
    double radial_v[3];
    double axis_length;
    double angle_u[3];
    double angle_v[3];
    double* energy_values;
    alea_rng_discrete_table_t* energy_table;
    alea_rng_tabular_table_t* energy_pdf_table;
    alea_rng_tabular_table_t* angle_table;
    double* r_edges;
    double* z_edges;
    alea_rng_discrete_table_t* rz_table;
    double phi_min;
    double phi_span;
    double integrated_emissivity;
};

alea_error_t alea_source_mixture_prepare(alea_source_t* const* components,
    const double* strengths, size_t count, alea_source_t** output) {
    if (!output) return ALEA_ERR_NULL_ARG;
    *output = NULL;
    if (!components || !strengths) return ALEA_ERR_NULL_ARG;
    if (!count || count > UINT32_MAX || count > SIZE_MAX / sizeof(alea_source_t*))
        return ALEA_ERR_INVALID_ARG;
    uint32_t mask = 0;
    for (size_t i = 0; i < count; ++i) {
        if (!components[i] || !isfinite(strengths[i]) || strengths[i] < 0.0)
            return ALEA_ERR_INVALID_ARG;
        for (size_t j = 0; j < i; ++j)
            if (components[i] == components[j]) return ALEA_ERR_INVALID_ARG;
        if (strengths[i] > 0.0)
            mask |= alea_source_particle_mask(components[i]);
    }
    alea_rng_discrete_table_t* table =
        alea_rng_discrete_table_create(strengths, count);
    if (!table) {
        alea_error_t err = (alea_error_t)alea_error_code();
        return err == ALEA_OK ? ALEA_ERR_INVALID_ARG : err;
    }
    alea_source_t* source = calloc(1, sizeof(*source));
    alea_source_t** owned = malloc(count*sizeof(*owned));
    if (!source || !owned) {
        free(source); free(owned);
        alea_rng_discrete_table_destroy(table);
        return ALEA_ERR_OUT_OF_MEMORY;
    }
    memcpy(owned, components, count*sizeof(*owned));
    source->is_mixture = 1;
    source->components = owned;
    source->component_count = count;
    source->component_mask = mask;
    source->component_table = table;
    *output = source;
    return ALEA_OK;
}

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
        (spec->space < ALEA_SOURCE_POINT || spec->space > ALEA_SOURCE_TOKAMAK_RZ) ||
        (spec->angle < ALEA_SOURCE_MONODIRECTIONAL ||
         spec->angle > ALEA_SOURCE_RADIAL) ||
        (spec->energy_type < ALEA_SOURCE_ENERGY_MONO ||
         spec->energy_type > ALEA_SOURCE_ENERGY_TABULATED) ||
        (spec->energy_type == ALEA_SOURCE_ENERGY_MONO &&
         (!isfinite(spec->energy) || spec->energy <= 0.0)) ||
        !isfinite(spec->time) || !isfinite(spec->weight) ||
        spec->weight <= 0.0) return ALEA_ERR_INVALID_ARG;
    if (spec->energy_type != ALEA_SOURCE_ENERGY_MONO) {
        if (!spec->energy_values || !spec->energy_weights ||
            spec->energy_count < (spec->energy_type == ALEA_SOURCE_ENERGY_LINES ? 1u : 2u) ||
            spec->energy_count > UINT32_MAX ||
            spec->energy_count > SIZE_MAX / sizeof(double)) return ALEA_ERR_INVALID_ARG;
        for (size_t i = 0; i < spec->energy_count; ++i)
            if (!isfinite(spec->energy_values[i]) || spec->energy_values[i] <= 0.0 ||
                (spec->energy_type == ALEA_SOURCE_ENERGY_TABULATED && i &&
                 spec->energy_values[i] <= spec->energy_values[i-1]))
                return ALEA_ERR_INVALID_ARG;
        if (spec->energy_type == ALEA_SOURCE_ENERGY_TABULATED &&
            spec->energy_interpolation != ALEA_SOURCE_PDF_HISTOGRAM &&
            spec->energy_interpolation != ALEA_SOURCE_PDF_LINEAR)
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
    if (spec->space == ALEA_SOURCE_TOKAMAK_RZ) {
        if (!spec->r_edges || !spec->z_edges || !spec->rz_emissivity ||
            spec->r_edge_count < 2 || spec->z_edge_count < 2 ||
            spec->r_edge_count > UINT32_MAX || spec->z_edge_count > UINT32_MAX ||
            spec->r_edge_count > SIZE_MAX / sizeof(double) ||
            spec->z_edge_count > SIZE_MAX / sizeof(double) ||
            spec->r_edge_count - 1 > SIZE_MAX / (spec->z_edge_count - 1) ||
            (spec->r_edge_count - 1)*(spec->z_edge_count - 1) > UINT32_MAX)
            return ALEA_ERR_INVALID_ARG;
        double phi_min = spec->phi_min;
        double phi_max = spec->phi_max;
        if (phi_min == 0.0 && phi_max == 0.0)
            phi_max = 6.2831853071795864769;
        double span = phi_max - phi_min;
        if (!isfinite(phi_min) || !isfinite(phi_max) ||
            !isfinite(span) || span <= 0.0 ||
            span > 6.2831853071795864769) return ALEA_ERR_INVALID_ARG;
        for (size_t i = 0; i < spec->r_edge_count; ++i)
            if (!isfinite(spec->r_edges[i]) || spec->r_edges[i] < 0.0 ||
                spec->r_edges[i] > 1e100 ||
                (i && spec->r_edges[i] <= spec->r_edges[i-1]))
                return ALEA_ERR_INVALID_ARG;
        for (size_t i = 0; i < spec->z_edge_count; ++i)
            if (!isfinite(spec->z_edges[i]) || fabs(spec->z_edges[i]) > 1e100 ||
                (i && spec->z_edges[i] <= spec->z_edges[i-1]))
                return ALEA_ERR_INVALID_ARG;
        size_t count = (spec->r_edge_count - 1)*(spec->z_edge_count - 1);
        for (size_t i = 0; i < count; ++i)
            if (!isfinite(spec->rz_emissivity[i]) ||
                spec->rz_emissivity[i] < 0.0) return ALEA_ERR_INVALID_ARG;
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
    if (spec->space == ALEA_SOURCE_TOKAMAK_RZ) {
        size_t nr = spec->r_edge_count - 1, nz = spec->z_edge_count - 1;
        size_t count = nr*nz;
        double* masses = malloc(count*sizeof(double));
        source->r_edges = malloc(spec->r_edge_count*sizeof(double));
        source->z_edges = malloc(spec->z_edge_count*sizeof(double));
        if (!masses || !source->r_edges || !source->z_edges) {
            free(masses); alea_source_free(source); return ALEA_ERR_OUT_OF_MEMORY;
        }
        memcpy(source->r_edges, spec->r_edges, spec->r_edge_count*sizeof(double));
        memcpy(source->z_edges, spec->z_edges, spec->z_edge_count*sizeof(double));
        source->phi_min = spec->phi_min;
        source->phi_span = spec->phi_min == 0.0 && spec->phi_max == 0.0 ?
            6.2831853071795864769 : spec->phi_max - spec->phi_min;
        double total = 0.0, correction = 0.0;
        for (size_t ir = 0; ir < nr; ++ir) {
            double rlo = spec->r_edges[ir], rhi = spec->r_edges[ir+1];
            double radial_area = 0.5*(rhi-rlo)*(rhi+rlo);
            for (size_t iz = 0; iz < nz; ++iz) {
                size_t index = ir*nz + iz;
                masses[index] = spec->rz_emissivity[index] * radial_area *
                    (spec->z_edges[iz+1] - spec->z_edges[iz]) * source->phi_span;
                if (!isfinite(masses[index])) {
                    free(masses); alea_source_free(source);
                    return ALEA_ERR_INVALID_ARG;
                }
                double adjusted = masses[index] - correction;
                double next = total + adjusted;
                correction = (next - total) - adjusted;
                total = next;
            }
        }
        if (!isfinite(total) || total <= 0.0) {
            free(masses); alea_source_free(source);
            return ALEA_ERR_INVALID_ARG;
        }
        source->integrated_emissivity = total;
        source->rz_table = alea_rng_discrete_table_create(masses, count);
        free(masses);
        if (!source->rz_table) {
            alea_error_t err = (alea_error_t)alea_error_code();
            alea_source_free(source);
            return err == ALEA_OK ? ALEA_ERR_INVALID_ARG : err;
        }
        source->spec.r_edges = source->r_edges;
        source->spec.z_edges = source->z_edges;
        source->spec.rz_emissivity = NULL;
    }
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
    } else if (spec->energy_type == ALEA_SOURCE_ENERGY_TABULATED) {
        source->energy_pdf_table = alea_rng_tabular_table_create(
            spec->energy_values, spec->energy_weights, spec->energy_count,
            spec->energy_interpolation == ALEA_SOURCE_PDF_HISTOGRAM ?
                ALEA_RNG_TABULAR_HISTOGRAM : ALEA_RNG_TABULAR_LIN_LIN);
        if (!source->energy_pdf_table) {
            alea_error_t err = (alea_error_t)alea_error_code();
            alea_source_free(source);
            return err == ALEA_OK ? ALEA_ERR_INVALID_ARG : err;
        }
        source->spec.energy_values = NULL;
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
    for (size_t i = 0; i < source->component_count; ++i)
        alea_source_free(source->components[i]);
    alea_rng_discrete_table_destroy(source->component_table);
    free(source->components);
    alea_rng_discrete_table_destroy(source->energy_table);
    alea_rng_tabular_table_destroy(source->energy_pdf_table);
    alea_rng_tabular_table_destroy(source->angle_table);
    alea_rng_discrete_table_destroy(source->rz_table);
    free(source->energy_values);
    free(source->r_edges);
    free(source->z_edges);
    free(source);
}

uint32_t alea_source_particle_mask(const alea_source_t* source) {
    if (!source) return 0;
    if (source->is_mixture) return source->component_mask;
    return source->spec.particle == ALEA_NUC_PARTICLE_NEUTRON ?
        ALEA_NUC_BIND_NEUTRON : ALEA_NUC_BIND_PHOTON;
}

alea_error_t alea_source_integrated_emissivity(const alea_source_t* source,
                                               double* output) {
    if (!source || !output) return ALEA_ERR_NULL_ARG;
    if (source->is_mixture) return ALEA_ERR_INVALID_ARG;
    if (source->spec.space != ALEA_SOURCE_TOKAMAK_RZ)
        return ALEA_ERR_INVALID_ARG;
    *output = source->integrated_emissivity;
    return ALEA_OK;
}

static alea_error_t source_sample_inner(const alea_source_t* source,
    uint64_t seed, uint32_t history_id, uint32_t depth,
    alea_transport_source_t* output) {
    if (!source || !output) return ALEA_ERR_NULL_ARG;
    if (depth >= 64) return ALEA_ERR_INVALID_STATE;
    if (source->is_mixture) {
        alea_rng_event_t event;
        size_t index;
        if (alea_rng_event_init(&event, ALEA_RNG_PHILOX4X32_10, seed,
            ALEA_RNG_DOMAIN_TRANSPORT_SOURCE_COMPONENT, history_id, depth) != 0 ||
            alea_rng_sample_discrete_alias(&event, source->component_table,
                                           &index) != 0) return ALEA_ERR_INVALID_STATE;
        return source_sample_inner(source->components[index], seed,
                                   history_id, depth + 1, output);
    }
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
    } else if (spec->energy_type == ALEA_SOURCE_ENERGY_TABULATED) {
        alea_rng_event_t event;
        if (alea_rng_event_init(&event, ALEA_RNG_PHILOX4X32_10, seed,
            ALEA_RNG_DOMAIN_TRANSPORT_SOURCE_ENERGY, history_id, 0) != 0 ||
            alea_rng_sample_tabular(&event, source->energy_pdf_table,
                                    &sample.particle.energy) != 0)
            return ALEA_ERR_INVALID_STATE;
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
    } else if (spec->space == ALEA_SOURCE_TOKAMAK_RZ) {
        alea_rng_event_t event;
        size_t index;
        if (alea_rng_event_init(&event, ALEA_RNG_PHILOX4X32_10, seed,
            ALEA_RNG_DOMAIN_TRANSPORT_SOURCE_POSITION, history_id, 0) != 0 ||
            alea_rng_sample_discrete_alias(&event, source->rz_table,
                                           &index) != 0)
            return ALEA_ERR_INVALID_STATE;
        double u_r, u_z, u_phi;
        if (alea_rng_event_next_uniform(&event, &u_r) != 0 ||
            alea_rng_event_next_uniform(&event, &u_z) != 0 ||
            alea_rng_event_next_uniform(&event, &u_phi) != 0)
            return ALEA_ERR_INVALID_STATE;
        size_t nz = spec->z_edge_count - 1;
        size_t ir = index / nz, iz = index % nz;
        double rlo = spec->r_edges[ir], rhi = spec->r_edges[ir+1];
        double radius = sqrt(rlo*rlo + u_r*(rhi*rhi - rlo*rlo));
        double phi = source->phi_min + u_phi*source->phi_span;
        sample.position[0] = radius*cos(phi);
        sample.position[1] = radius*sin(phi);
        sample.position[2] = spec->z_edges[iz] +
            u_z*(spec->z_edges[iz+1] - spec->z_edges[iz]);
        for (int i = 0; i < 3; ++i)
            if (!isfinite(sample.position[i])) return ALEA_ERR_OVERFLOW;
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

alea_error_t alea_source_sample(void* context, uint64_t seed,
                                uint32_t history_id,
                                alea_transport_source_t* output) {
    return source_sample_inner(context, seed, history_id, 0, output);
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
