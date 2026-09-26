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
#define ALEA_CELL_PARAM_PWT    (1u << 1)
#define ALEA_CELL_PARAM_NONU   (1u << 2)
#define ALEA_CELL_PARAM_PD     (1u << 3)
#define ALEA_CELL_PARAM_ELPT   (1u << 4)
#define ALEA_CELL_PARAM_UNC    (1u << 5)
#define ALEA_CELL_PARAM_BFLCL  (1u << 6)

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
    double photon_weight;
    double detector_contribution;
    double energy_cutoff;
    int fission_turnoff;
    int uncollided_secondaries;
    int magnetic_field;
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

#ifdef __cplusplus
}
#endif

#endif
