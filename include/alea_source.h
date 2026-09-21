// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/** @file alea_source.h
 * Prepared primary sources for fixed-source transport and data-free previews.
 */
#ifndef ALEA_SOURCE_H
#define ALEA_SOURCE_H

#include "alea_transport.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct alea_source alea_source_t;

typedef enum {
    ALEA_SOURCE_POINT, ALEA_SOURCE_BOX, ALEA_SOURCE_LINE,
    ALEA_SOURCE_SPHERE, ALEA_SOURCE_CYLINDER
} alea_source_space_t;
typedef enum {
    ALEA_SOURCE_MONODIRECTIONAL, ALEA_SOURCE_ISOTROPIC,
    ALEA_SOURCE_CONE, ALEA_SOURCE_COSINE, ALEA_SOURCE_TABULATED_MU,
    ALEA_SOURCE_RADIAL
}
    alea_source_angle_t;
typedef enum { ALEA_SOURCE_PDF_HISTOGRAM, ALEA_SOURCE_PDF_LINEAR }
    alea_source_pdf_t;
typedef enum { ALEA_SOURCE_ENERGY_MONO, ALEA_SOURCE_ENERGY_LINES }
    alea_source_energy_t;

typedef struct {
    alea_nuc_particle_t particle;
    alea_source_space_t space;
    alea_source_angle_t angle;
    double position[3]; /* point, cm */
    double lower[3];    /* box, cm */
    double upper[3];    /* box, cm */
    double start[3];    /* line endpoint, cm */
    double end[3];      /* line endpoint, cm */
    double center[3];   /* spherical volume center, cm */
    double base[3];     /* cylindrical volume base center, cm */
    double axis[3];     /* cylinder base-to-top vector, cm */
    double inner_radius; /* sphere/cylinder shell, cm; zero for full volume */
    double outer_radius; /* sphere/cylinder shell, cm */
    double direction[3]; /* mono direction or cone/cosine/tabulated axis */
    double cone_half_angle; /* cone half-angle in radians, [0, pi] */
    const double* angle_mu; /* tabulated polar cosine coordinates, [-1, 1] */
    const double* angle_pdf; /* density per unit mu; copied at preparation */
    size_t angle_count; /* tabulated coordinates and densities */
    alea_source_pdf_t angle_interpolation; /* histogram or linear */
    double angle_origin[3]; /* radial direction origin, cm */
    int radial_inward; /* radial: point toward origin when nonzero */
    double energy; /* constant MeV; positive */
    alea_source_energy_t energy_type; /* zero: mono */
    const double* energy_values; /* lines: positive MeV; copied at preparation */
    const double* energy_weights; /* lines: nonnegative masses; copied */
    size_t energy_count; /* lines only; at least one positive weight */
    double time;   /* constant seconds */
    double weight; /* positive; one for an unbiased physical source */
} alea_source_spec_t;

/** Validate and copy a source. No geometry or nuclear data is required.
 * An accepted source owns its description and may be reused across runs. */
alea_error_t alea_source_prepare(const alea_source_spec_t* spec,
                                 alea_source_t** output);
void alea_source_free(alea_source_t* source);
uint32_t alea_source_particle_mask(const alea_source_t* source);

/** Sample a global history ID. Seed/history identity is shared with transport.
 * This function also has the signature of alea_transport_source_sampler_fn. */
alea_error_t alea_source_sample(void* source, uint64_t seed,
                                uint32_t history_id,
                                alea_transport_source_t* output);

/** Sample consecutive global history IDs into caller-owned storage. */
alea_error_t alea_source_sample_batch(const alea_source_t* source,
    uint64_t seed, uint32_t history_offset, uint32_t count,
    alea_transport_source_t* output);

#ifdef __cplusplus
}
#endif

#endif
