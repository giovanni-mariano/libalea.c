// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

#include "alea_model.h"
#include "core/alea_system.h"
#include "util/compat.h"
#include "util/alea_alloc.h"
#include "util/cell_parameter.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

struct alea_model {
    alea_system_t* sys;
    int owns_sys;
    char* name;
    char* title;
    char* comments;
    alea_model_cell_metadata_t* cells;
    size_t cell_count;
    size_t cell_capacity;
};

static void cell_metadata_free(alea_model_cell_metadata_t* meta) {
    if (!meta) return;
    free(meta->name);
    free(meta->detector_probabilities);
    memset(meta, 0, sizeof(*meta));
}

static void cell_metadata_init(alea_model_cell_metadata_t* meta) {
    memset(meta, 0, sizeof(*meta));
    meta->importance_neutron = 1.0;
    meta->importance_photon = 1.0;
    meta->importance_electron = 1.0;
    meta->photon_production.weight_basis = ALEA_PHOTON_WEIGHT_SOURCE_RELATIVE;
    meta->photon_production.weight_threshold = 1.0;
    meta->detector_contribution = 1.0;
    meta->fission_mode = ALEA_FISSION_NORMAL;
    meta->secondary_collision_state = ALEA_SECONDARY_UNCOLLIDED;
    for (int i = 0; i < ALEA_PARTICLE_COUNT; i++)
        meta->particle_secondary_state[i] = ALEA_SECONDARY_UNCOLLIDED;
}

static int reserve_cells(alea_model_t* model, size_t needed) {
    if (needed <= model->cell_capacity) return 0;
    const size_t max = SIZE_MAX / sizeof(*model->cells);
    if (needed > max) return -1;
    size_t capacity = model->cell_capacity ? model->cell_capacity : 8;
    while (capacity < needed) {
        if (capacity > max / 2) { capacity = needed; break; }
        capacity *= 2;
    }
    alea_model_cell_metadata_t* cells = alea_realloc(
        model->cells, capacity * sizeof(*cells));
    if (!cells) return -1;
    memset(cells + model->cell_capacity, 0,
           (capacity - model->cell_capacity) * sizeof(*cells));
    model->cells = cells;
    model->cell_capacity = capacity;
    return 0;
}

static int model_cells_reserve(void* userdata, size_t count) {
    return reserve_cells(userdata, count);
}

static void model_cell_added(void* userdata, size_t index) {
    alea_model_t* model = userdata;
    if (!model || index != model->cell_count || index >= model->cell_capacity) return;
    cell_metadata_init(&model->cells[index]);
    model->cell_count++;
}

static int model_cell_copied(void* userdata, size_t dst, size_t src) {
    alea_model_t* model = userdata;
    if (!model || src >= model->cell_count) return -1;
    if (dst == src) return 0;
    if (dst >= model->cell_count) model_cell_added(userdata, dst);
    if (dst >= model->cell_count) return -1;
    char* name = model->cells[src].name ? alea_alloc_strdup(model->cells[src].name) : NULL;
    alea_detector_probability_t* probabilities = alea_detector_probability_copy(
        model->cells[src].detector_probabilities,
        model->cells[src].detector_probability_count);
    if ((model->cells[src].name && !name) ||
        (model->cells[src].detector_probability_count && !probabilities)) {
        free(name);
        free(probabilities);
        alea_set_error_detail(ALEA_ERR_OUT_OF_MEMORY, "copying cell parameters");
        return -1;
    }
    cell_metadata_free(&model->cells[dst]);
    model->cells[dst] = model->cells[src];
    model->cells[dst].name = name;
    model->cells[dst].detector_probabilities = probabilities;
    return 0;
}

static void model_cell_removed(void* userdata, size_t index) {
    alea_model_t* model = userdata;
    if (!model || index >= model->cell_count) return;
    cell_metadata_free(&model->cells[index]);
    if (index + 1 < model->cell_count) {
        memmove(&model->cells[index], &model->cells[index + 1],
                (model->cell_count - index - 1) * sizeof(*model->cells));
    }
    model->cell_count--;
    memset(&model->cells[model->cell_count], 0, sizeof(*model->cells));
}

static void model_cells_pruned(void* userdata) {
    alea_model_t* model = userdata;
    size_t write = 0;
    for (size_t i = 0; i < model->cell_count; i++) {
        if (model->sys->cells.data[i].root_node_id == ALEA_NODE_ID_INVALID)
            cell_metadata_free(&model->cells[i]);
        else model->cells[write++] = model->cells[i];
    }
    if (write < model->cell_count)
        memset(model->cells + write, 0, (model->cell_count - write) * sizeof(*model->cells));
    model->cell_count = write;
}

static alea_model_t* model_create(alea_system_t* sys, int owns_sys) {
    if (!sys) return NULL;
    if (sys->on_cells_reserve || sys->on_cell_added || sys->on_cell_copied || sys->on_cell_removed || sys->on_cells_pruned) {
        alea_set_error_detail(ALEA_ERR_INVALID_STATE,
                              "system is already attached to a model wrapper");
        return NULL;
    }
    alea_model_t* model = calloc(1, sizeof(*model));
    if (!model) return NULL;
    model->sys = sys;
    model->owns_sys = owns_sys;
    size_t count = alea_vec_count(&sys->cells);
    if (count && reserve_cells(model, count)) {
        free(model);
        return NULL;
    }
    model->cell_count = count;
    for (size_t i = 0; i < count; i++) cell_metadata_init(&model->cells[i]);
    sys->cell_hook_userdata = model;
    sys->on_cells_reserve = model_cells_reserve;
    sys->on_cell_added = model_cell_added;
    sys->on_cell_copied = model_cell_copied;
    sys->on_cell_removed = model_cell_removed;
    sys->on_cells_pruned = model_cells_pruned;
    return model;
}

alea_model_t* alea_model_wrap(alea_system_t* sys) { return model_create(sys, 0); }
alea_model_t* alea_model_adopt(alea_system_t* sys) { return model_create(sys, 1); }

static void model_unhook(alea_model_t* model) {
    if (!model || !model->sys) return;
    if (model->sys->cell_hook_userdata == model) {
        model->sys->cell_hook_userdata = NULL;
        model->sys->on_cells_reserve = NULL;
        model->sys->on_cell_added = NULL;
        model->sys->on_cell_copied = NULL;
        model->sys->on_cell_removed = NULL;
        model->sys->on_cells_pruned = NULL;
    }
}

void alea_model_destroy(alea_model_t* model) {
    if (!model) return;
    model_unhook(model);
    for (size_t i = 0; i < model->cell_count; i++) cell_metadata_free(&model->cells[i]);
    free(model->cells);
    free(model->name);
    free(model->title);
    free(model->comments);
    if (model->owns_sys) alea_destroy(model->sys);
    free(model);
}

alea_system_t* alea_model_system(alea_model_t* model) { return model ? model->sys : NULL; }
const alea_system_t* alea_model_system_const(const alea_model_t* model) { return model ? model->sys : NULL; }

alea_system_t* alea_model_take_system(alea_model_t* model) {
    if (!model || !model->sys) return NULL;
    model_unhook(model);
    alea_system_t* sys = model->sys;
    model->sys = NULL;
    model->owns_sys = 0;
    return sys;
}

static int set_string(char** target, const char* value) {
    char* copy = value ? alea_strdup(value) : NULL;
    if (value && !copy) return -1;
    free(*target);
    *target = copy;
    return 0;
}

const char* alea_model_name(const alea_model_t* model) { return model ? model->name : NULL; }
const char* alea_model_title(const alea_model_t* model) { return model ? model->title : NULL; }
const char* alea_model_comments(const alea_model_t* model) { return model ? model->comments : NULL; }
int alea_model_set_name(alea_model_t* model, const char* value) { return model ? set_string(&model->name, value) : -1; }
int alea_model_set_title(alea_model_t* model, const char* value) { return model ? set_string(&model->title, value) : -1; }
int alea_model_set_comments(alea_model_t* model, const char* value) { return model ? set_string(&model->comments, value) : -1; }

size_t alea_model_cell_metadata_count(const alea_model_t* model) { return model ? model->cell_count : 0; }
const alea_model_cell_metadata_t* alea_model_cell_metadata(const alea_model_t* model, size_t index) {
    return model && index < model->cell_count ? &model->cells[index] : NULL;
}
alea_model_cell_metadata_t* alea_model_cell_metadata_mut(alea_model_t* model, size_t index) {
    return model && index < model->cell_count ? &model->cells[index] : NULL;
}
int alea_model_cell_set_name(alea_model_t* model, size_t index, const char* name) {
    return model && index < model->cell_count ? set_string(&model->cells[index].name, name) : -1;
}

int alea_model_cell_set_detector_probability(alea_model_t* model,
    size_t index, int tally, double probability) {
    alea_model_cell_metadata_t* meta = alea_model_cell_metadata_mut(model, index);
    if (!meta || tally < 0 || !isfinite(probability) ||
        probability < 0 || probability > 1) return -1;
    if (!tally) {
        meta->detector_contribution = probability;
        meta->parameter_flags |= ALEA_CELL_PARAM_DETECTOR_PROBABILITY;
        return 0;
    }
    return alea_detector_probability_store(&meta->detector_probabilities,
        &meta->detector_probability_count, tally, probability);
}
