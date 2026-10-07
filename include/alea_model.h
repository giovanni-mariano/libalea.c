// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/** @file alea_model.h @brief ALEA model and non-geometric metadata. */

#ifndef ALEA_MODEL_H
#define ALEA_MODEL_H

#include "alea.h"
#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct alea_model alea_model_t;

#define ALEA_CELL_PARAM_VOLUME (1u << 0)
#define ALEA_CELL_PARAM_PHOTON_PRODUCTION (1u << 1)
#define ALEA_CELL_PARAM_FISSION_MODE (1u << 2)
#define ALEA_CELL_PARAM_DETECTOR_PROBABILITY (1u << 3)
#define ALEA_CELL_PARAM_ENERGY_CUTOFF (1u << 4)
#define ALEA_CELL_PARAM_SECONDARY_STATE (1u << 5)
#define ALEA_CELL_PARAM_MAGNETIC_FIELD (1u << 6)
/* Compatibility aliases for the original flag names. */
#define ALEA_CELL_PARAM_PWT ALEA_CELL_PARAM_PHOTON_PRODUCTION
#define ALEA_CELL_PARAM_NONU ALEA_CELL_PARAM_FISSION_MODE
#define ALEA_CELL_PARAM_PD ALEA_CELL_PARAM_DETECTOR_PROBABILITY
#define ALEA_CELL_PARAM_ELPT ALEA_CELL_PARAM_ENERGY_CUTOFF
#define ALEA_CELL_PARAM_UNC ALEA_CELL_PARAM_SECONDARY_STATE
#define ALEA_CELL_PARAM_BFLCL ALEA_CELL_PARAM_MAGNETIC_FIELD

typedef enum {
    ALEA_PARTICLE_NEUTRON, ALEA_PARTICLE_PHOTON, ALEA_PARTICLE_ELECTRON,
    ALEA_PARTICLE_COUNT
} alea_particle_t;

typedef enum {
    ALEA_PHOTON_PRODUCTION_THRESHOLD,
    ALEA_PHOTON_PRODUCTION_ONE_PER_COLLISION,
    ALEA_PHOTON_PRODUCTION_OFF
} alea_photon_production_mode_t;

typedef enum {
    ALEA_PHOTON_WEIGHT_ABSOLUTE, ALEA_PHOTON_WEIGHT_SOURCE_RELATIVE
} alea_photon_weight_basis_t;

typedef struct {
    alea_photon_production_mode_t mode;
    alea_photon_weight_basis_t weight_basis;
    double weight_threshold; /* Positive magnitude, independent of MCNP signs. */
} alea_photon_production_t;

typedef enum {
    ALEA_FISSION_NORMAL,
    ALEA_FISSION_CAPTURE_WITH_PHOTONS,
    ALEA_FISSION_CAPTURE_WITHOUT_PHOTONS
} alea_fission_mode_t;

typedef enum {
    ALEA_SECONDARY_COLLIDED, ALEA_SECONDARY_UNCOLLIDED
} alea_secondary_collision_state_t;

typedef struct {
    int tally; /* Positive tally ID; the all-tally default is stored separately. */
    double probability;
} alea_detector_probability_t;

/** Optional metadata parallel to one system cell. */
typedef struct {
    char* name;
    double importance_neutron;
    double importance_photon;
    double importance_electron;
    unsigned has_importance_neutron : 1;
    unsigned has_importance_photon : 1;
    unsigned has_importance_electron : 1;
    double user_volume;
    alea_photon_production_t photon_production;
    double detector_contribution;
    double energy_cutoff; /* Default for all supported particles, in MeV. */
    alea_fission_mode_t fission_mode;
    alea_secondary_collision_state_t secondary_collision_state;
    int magnetic_field; /* Reference to a magnetic-field definition; zero = none. */
    double particle_energy_cutoff[ALEA_PARTICLE_COUNT];
    alea_secondary_collision_state_t particle_secondary_state[ALEA_PARTICLE_COUNT];
    uint32_t energy_cutoff_particles;
    uint32_t secondary_state_particles;
    /* Owned by the model; use the setter below to insert/update entries. */
    alea_detector_probability_t* detector_probabilities;
    size_t detector_probability_count;
    uint32_t parameter_flags;
} alea_model_cell_metadata_t;

/** Wrap a system. The caller retains ownership of the system. */
alea_model_t* alea_model_wrap(alea_system_t* sys);

/** Create a model that takes ownership of the supplied system. */
alea_model_t* alea_model_adopt(alea_system_t* sys);

void alea_model_destroy(alea_model_t* model);
alea_system_t* alea_model_system(alea_model_t* model);
const alea_system_t* alea_model_system_const(const alea_model_t* model);
alea_system_t* alea_model_take_system(alea_model_t* model);

const char* alea_model_name(const alea_model_t* model);
const char* alea_model_title(const alea_model_t* model);
const char* alea_model_comments(const alea_model_t* model);
int alea_model_set_name(alea_model_t* model, const char* value);
int alea_model_set_title(alea_model_t* model, const char* value);
int alea_model_set_comments(alea_model_t* model, const char* value);

size_t alea_model_cell_metadata_count(const alea_model_t* model);
const alea_model_cell_metadata_t* alea_model_cell_metadata(
    const alea_model_t* model, size_t cell_index);
alea_model_cell_metadata_t* alea_model_cell_metadata_mut(
    alea_model_t* model, size_t cell_index);
int alea_model_cell_set_name(alea_model_t* model, size_t cell_index,
                             const char* name);
/* tally=0 sets the all-tally default. Probability must be in [0,1]. */
int alea_model_cell_set_detector_probability(alea_model_t* model,
    size_t cell_index, int tally, double probability);

#ifdef __cplusplus
}
#endif

#endif
