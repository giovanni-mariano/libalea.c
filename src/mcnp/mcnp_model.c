// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/**
 * @file mcnp_model.c
 * @brief MCNP model implementation
 *
 * Manages the mcnp_model_t lifecycle and its parallel cell params array.
 * Cell event callbacks keep the params array in sync with sys->cells.
 */

#include "mcnp_model.h"
#include "alea.h"
#include "alea_model.h"
#include "core/alea_system.h"
#include "core/alea_export.h"
#include "util/alea_log.h"
#include <errno.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "util/compat.h"

static int load_profile_enabled(void) {
    const char* v = getenv("ALEA_PROFILE_LOAD");
    return v && *v && strcmp(v, "0") != 0;
}

const mcnp_export_config_t MCNP_EXPORT_CONFIG_DEFAULT = {
    .surface_policy = 0,    /* ALEA_EMIT_MACROBODY */
    .trcl_mode = 0,         /* preserve */
    .transform_mode = 0,    /* original */
    .mcnp_max_col = 80,
    .mcnp_cont_indent = 5,
};

/* ============================================================================
 * INTERNAL: Default cell params
 * ============================================================================ */

static void init_default_params(mcnp_cell_params_t* p) {
    memset(p, 0, sizeof(*p));
    p->imp_n = 1.0;
    p->imp_p = 1.0;
    p->imp_e = 1.0;
    p->trcl_inline_index = MCNP_INLINE_TRANSFORM_INVALID;
    p->fill_transform_index = MCNP_INLINE_TRANSFORM_INVALID;
}

static int set_string(char** target, const char* value) {
    char* copy = value ? alea_strdup(value) : NULL;
    if (value && !copy) return -1;
    free(*target);
    *target = copy;
    return 0;
}

/* ============================================================================
 * CELL EVENT CALLBACKS
 * ============================================================================ */

static void on_cell_added_cb(void* ud, size_t new_index) {
    mcnp_model_t* model = (mcnp_model_t*)ud;
    (void)new_index;
    mcnp_model_add_params(model);
}

static void on_cell_copied_cb(void* ud, size_t dst_index, size_t src_index) {
    mcnp_model_t* model = (mcnp_model_t*)ud;
    if (dst_index < model->cell_params_count && src_index < model->cell_params_count) {
        if (dst_index == src_index) return;
        mcnp_scoped_cell_params_t scoped = {0};
        if (mcnp_scoped_params_copy(&scoped, &model->cell_params[src_index].scoped)) {
            alea_set_error_detail(ALEA_ERR_OUT_OF_MEMORY, "copying MCNP cell parameters");
            return;
        }
        mcnp_scoped_params_free(&model->cell_params[dst_index].scoped);
        model->cell_params[dst_index] = model->cell_params[src_index];
        model->cell_params[dst_index].scoped = scoped;
        const char* source_name = model->cell_names[src_index];
        (void)set_string(&model->cell_names[dst_index], source_name);
    }
}

static void on_cell_removed_cb(void* ud, size_t index) {
    mcnp_model_t* model = (mcnp_model_t*)ud;
    if (!model || index >= model->cell_params_count) return;
    free(model->cell_names[index]);
    mcnp_scoped_params_free(&model->cell_params[index].scoped);
    const size_t trailing = model->cell_params_count - index - 1;
    if (trailing > 0) {
        memmove(&model->cell_params[index], &model->cell_params[index + 1],
                trailing * sizeof(*model->cell_params));
        memmove(&model->cell_names[index], &model->cell_names[index + 1],
                trailing * sizeof(*model->cell_names));
    }
    model->cell_params_count--;
    model->cell_names[model->cell_params_count] = NULL;
    memset(&model->cell_params[model->cell_params_count], 0, sizeof(*model->cell_params));
}

/* ============================================================================
 * PUBLIC API
 * ============================================================================ */

int mcnp_model_reserve_params(mcnp_model_t* model, size_t cap) {
    if (!model) return -1;
    if (cap <= model->cell_params_capacity) return 0;

    size_t new_cap = model->cell_params_capacity;
    if (new_cap == 0) new_cap = 64;
    while (new_cap < cap) new_cap *= 2;

    mcnp_cell_params_t* new_arr = realloc(model->cell_params,
                                           new_cap * sizeof(mcnp_cell_params_t));
    if (!new_arr) return -1;
    model->cell_params = new_arr;

    char** new_names = realloc(model->cell_names,
                               new_cap * sizeof(*model->cell_names));
    if (!new_names) return -1;
    memset(new_names + model->cell_params_capacity, 0,
           (new_cap - model->cell_params_capacity) * sizeof(*new_names));
    model->cell_names = new_names;
    model->cell_params_capacity = new_cap;
    return 0;
}

int mcnp_model_add_params(mcnp_model_t* model) {
    if (!model) return -1;

    if (model->cell_params_count >= model->cell_params_capacity) {
        if (mcnp_model_reserve_params(model, model->cell_params_count + 1) < 0)
            return -1;
    }

    size_t idx = model->cell_params_count++;
    init_default_params(&model->cell_params[idx]);
    model->cell_names[idx] = NULL;
    return (int)idx;
}

mcnp_cell_params_t* mcnp_cell_params(mcnp_model_t* m, size_t idx) {
    if (!m || idx >= m->cell_params_count) return NULL;
    return &m->cell_params[idx];
}

const mcnp_cell_params_t* mcnp_cell_params_const(const mcnp_model_t* m, size_t idx) {
    if (!m || idx >= m->cell_params_count) return NULL;
    return &m->cell_params[idx];
}

const char* mcnp_model_name(const mcnp_model_t* model) {
    return model ? model->name : NULL;
}

const char* mcnp_model_title(const mcnp_model_t* model) {
    return model ? model->title : NULL;
}

const char* mcnp_model_comments(const mcnp_model_t* model) {
    return model ? model->comments : NULL;
}

int mcnp_model_set_name(mcnp_model_t* model, const char* value) {
    return model ? set_string(&model->name, value) : -1;
}

int mcnp_model_set_title(mcnp_model_t* model, const char* value) {
    return model ? set_string(&model->title, value) : -1;
}

int mcnp_model_set_comments(mcnp_model_t* model, const char* value) {
    return model ? set_string(&model->comments, value) : -1;
}

const char* mcnp_model_cell_name(const mcnp_model_t* model, size_t index) {
    return model && index < model->cell_params_count
        ? model->cell_names[index] : NULL;
}

int mcnp_model_cell_set_name(mcnp_model_t* model, size_t index,
                             const char* value) {
    return model && index < model->cell_params_count
        ? set_string(&model->cell_names[index], value) : -1;
}

static int reserve_inline_transforms(mcnp_model_t* model, size_t cap) {
    if (cap <= model->inline_transform_capacity) return 0;
    const size_t max_elems = SIZE_MAX / sizeof(*model->inline_transforms);
    if (cap > max_elems) {
        alea_set_error_detail(ALEA_ERR_OUT_OF_MEMORY,
                              "Inline transform count %zu exceeds allocation limit",
                              cap);
        return -1;
    }

    size_t new_cap = model->inline_transform_capacity;
    if (new_cap == 0) new_cap = 16;
    while (new_cap < cap) {
        if (new_cap > max_elems / 2) {
            new_cap = cap;
            break;
        }
        new_cap *= 2;
    }

    mcnp_inline_transform_t* new_arr = realloc(
        model->inline_transforms, new_cap * sizeof(*new_arr));
    if (!new_arr) {
        alea_set_error_detail(ALEA_ERR_OUT_OF_MEMORY,
                              "Failed to grow inline transforms to %zu elements",
                              new_cap);
        return -1;
    }

    model->inline_transforms = new_arr;
    model->inline_transform_capacity = new_cap;
    return 0;
}

uint32_t mcnp_model_add_inline_transform(mcnp_model_t* model,
                                         const double* values,
                                         int count,
                                         int degrees) {
    if (!model || !values || count <= 0 || count > 13) {
        return MCNP_INLINE_TRANSFORM_INVALID;
    }

    size_t idx = model->inline_transform_count;
    if (idx >= UINT32_MAX) {
        return MCNP_INLINE_TRANSFORM_INVALID;
    }

    if (reserve_inline_transforms(model, idx + 1) < 0) {
        return MCNP_INLINE_TRANSFORM_INVALID;
    }
    mcnp_inline_transform_t* tr = &model->inline_transforms[idx];
    model->inline_transform_count++;

    memset(tr, 0, sizeof(*tr));
    tr->count = (uint8_t)count;
    tr->degrees = (uint8_t)(degrees != 0);
    memcpy(tr->values, values, (size_t)count * sizeof(double));
    return (uint32_t)idx;
}

const mcnp_inline_transform_t* mcnp_model_inline_transform_const(
    const mcnp_model_t* model,
    uint32_t index) {
    if (!model || index == MCNP_INLINE_TRANSFORM_INVALID ||
        (size_t)index >= model->inline_transform_count) {
        return NULL;
    }
    return &model->inline_transforms[index];
}

static mcnp_model_t* finalize_loaded_model(mcnp_model_t* model) {
    if (!model) return NULL;

    double t0 = alea_monotonic_seconds();
    if (alea_validate_cell_ids(model->sys) < 0) {
        alea_set_error_detail(ALEA_ERR_INVALID_ID, "Duplicate cell IDs");
        mcnp_model_destroy(model);
        return NULL;
    }
    double t1 = alea_monotonic_seconds();
    if (load_profile_enabled()) {
        fprintf(stderr, "[alea-load-profile] %-28s %.6f s\n",
                "validate_cell_ids", t1 - t0);
    }

    model->sys->source = ALEA_SOURCE_MCNP;

    return model;
}

void mcnp_model_register_hooks(mcnp_model_t* model) {
    if (!model || !model->sys) return;
    model->sys->cell_hook_userdata = model;
    model->sys->on_cell_added = on_cell_added_cb;
    model->sys->on_cell_copied = on_cell_copied_cb;
    model->sys->on_cell_removed = on_cell_removed_cb;
}

void mcnp_model_destroy(mcnp_model_t* model) {
    if (!model) return;

    /* Unhook before destroying */
    if (model->sys) {
        model->sys->cell_hook_userdata = NULL;
        model->sys->on_cell_added = NULL;
        model->sys->on_cell_copied = NULL;
        model->sys->on_cell_removed = NULL;
        if (model->owns_sys) {
            alea_system_destroy(model->sys);
        }
    }

    for (size_t i = 0; i < model->cell_params_count; i++) {
        free(model->cell_names[i]);
        mcnp_scoped_params_free(&model->cell_params[i].scoped);
    }
    free(model->cell_names);
    free(model->name);
    free(model->title);
    free(model->comments);
    free(model->cell_params);
    free(model->inline_transforms);
    free(model);
}

alea_system_t* mcnp_model_system(mcnp_model_t* model) {
    return model ? model->sys : NULL;
}

alea_system_t* mcnp_model_take_system(mcnp_model_t* model) {
    if (!model) return NULL;

    alea_system_t* sys = model->sys;
    if (!sys) return NULL;

    /* The caller takes only the generic system, not this MCNP sidecar.  Do
     * not leave mutation hooks pointing at the model after it is destroyed. */
    sys->cell_hook_userdata = NULL;
    sys->on_cell_added = NULL;
    sys->on_cell_copied = NULL;
    sys->on_cell_removed = NULL;

    model->sys = NULL;
    model->owns_sys = 0;
    return sys;
}

mcnp_model_t* mcnp_model_wrap(alea_system_t* sys) {
    if (!sys) return NULL;
    if (sys->on_cell_added || sys->on_cell_copied || sys->on_cell_removed) {
        alea_set_error_detail(ALEA_ERR_INVALID_STATE,
                              "system is already attached to a model wrapper");
        return NULL;
    }

    mcnp_model_t* model = calloc(1, sizeof(mcnp_model_t));
    if (!model) return NULL;

    model->sys = sys;
    model->owns_sys = 0;
    model->export_config = MCNP_EXPORT_CONFIG_DEFAULT;

    /* Populate cell params with defaults */
    size_t n = alea_vec_count(&sys->cells);
    if (n > 0) {
        if (mcnp_model_reserve_params(model, n) < 0) {
            free(model);
            return NULL;
        }
        for (size_t i = 0; i < n; i++) {
            mcnp_model_add_params(model);
        }
    }

    /* Keep params synchronized with future cell mutations on the wrapped system. */
    mcnp_model_register_hooks(model);

    return model;
}

static int params_to_metadata(const mcnp_cell_params_t* p,
                               alea_model_cell_metadata_t* m) {
    if (!p || !m) return -1;
    m->importance_neutron = p->has_imp_n ? p->imp_n : 1.0;
    m->importance_photon = p->has_imp_p ? p->imp_p : 1.0;
    m->importance_electron = p->has_imp_e ? p->imp_e : 1.0;
    m->has_importance_neutron = p->has_imp_n;
    m->has_importance_photon = p->has_imp_p;
    m->has_importance_electron = p->has_imp_e;
    if (p->has_vol) { m->user_volume = p->vol; m->parameter_flags |= ALEA_CELL_PARAM_VOLUME; }
    if (p->has_pwt) {
        if (alea_photon_production_from_mcnp(p->pwt, &m->photon_production)) return -1;
        m->parameter_flags |= ALEA_CELL_PARAM_PHOTON_PRODUCTION;
    }
    if (p->has_nonu) {
        if (alea_fission_mode_from_mcnp(p->nonu, &m->fission_mode)) return -1;
        m->parameter_flags |= ALEA_CELL_PARAM_FISSION_MODE;
    }
    if (p->has_pd) { m->detector_contribution = p->pd; m->parameter_flags |= ALEA_CELL_PARAM_PD; }
    if (p->has_elpt) { m->energy_cutoff = p->elpt; m->parameter_flags |= ALEA_CELL_PARAM_ELPT; }
    if (p->has_unc) { m->secondary_collision_state = (alea_secondary_collision_state_t)p->unc; m->parameter_flags |= ALEA_CELL_PARAM_UNC; }
    if (p->has_bflcl) { m->magnetic_field = p->bflcl; m->parameter_flags |= ALEA_CELL_PARAM_BFLCL; }
    m->energy_cutoff_particles = p->scoped.energy_cutoff_particles;
    m->secondary_state_particles = p->scoped.secondary_state_particles;
    for (int i = 0; i < ALEA_PARTICLE_COUNT; i++) {
        m->particle_energy_cutoff[i] = p->scoped.energy_cutoff[i];
        m->particle_secondary_state[i] = (alea_secondary_collision_state_t)p->scoped.secondary_state[i];
    }
    m->detector_probabilities = alea_detector_probability_copy(
        p->scoped.detector_probabilities, p->scoped.detector_probability_count);
    if (p->scoped.detector_probability_count && !m->detector_probabilities) return -1;
    m->detector_probability_count = p->scoped.detector_probability_count;
    if (!alea_cell_parameters_valid(m)) {
        alea_set_error_detail(ALEA_ERR_PARSE_ERROR, "MCNP cell parameters cannot be represented in ALEA");
        return -1;
    }
    return 0;
}

static int metadata_to_params(const alea_model_cell_metadata_t* m,
                               mcnp_cell_params_t* p) {
    if (!m || !p) return -1;
    if (!alea_cell_parameters_valid(m)) {
        alea_set_error_detail(ALEA_ERR_EXPORT_FAILED, "invalid ALEA cell parameters for MCNP conversion");
        return -1;
    }
    p->imp_n = m->importance_neutron;
    p->imp_p = m->importance_photon;
    p->imp_e = m->importance_electron;
    p->has_imp_n = m->has_importance_neutron;
    p->has_imp_p = m->has_importance_photon;
    p->has_imp_e = m->has_importance_electron;
    if (m->parameter_flags & ALEA_CELL_PARAM_VOLUME) { p->vol = m->user_volume; p->has_vol = 1; }
    if (m->parameter_flags & ALEA_CELL_PARAM_PHOTON_PRODUCTION) {
        if (alea_photon_production_to_mcnp(&m->photon_production, &p->pwt)) {
            alea_set_error_detail(ALEA_ERR_EXPORT_FAILED,
                "photon production threshold cannot be represented in MCNP");
            return -1;
        }
        p->has_pwt = 1;
    }
    if (m->parameter_flags & ALEA_CELL_PARAM_FISSION_MODE) {
        if (alea_fission_mode_to_mcnp(m->fission_mode, &p->nonu)) return -1;
        p->has_nonu = 1;
    }
    if (m->parameter_flags & ALEA_CELL_PARAM_PD) { p->pd = m->detector_contribution; p->has_pd = 1; }
    if (m->parameter_flags & ALEA_CELL_PARAM_ELPT) { p->elpt = m->energy_cutoff; p->has_elpt = 1; }
    if (m->parameter_flags & ALEA_CELL_PARAM_UNC) { p->unc = (int)m->secondary_collision_state; p->has_unc = 1; }
    if (m->parameter_flags & ALEA_CELL_PARAM_BFLCL) { p->bflcl = m->magnetic_field; p->has_bflcl = 1; }
    p->scoped.energy_cutoff_particles = m->energy_cutoff_particles;
    p->scoped.secondary_state_particles = m->secondary_state_particles;
    for (int i = 0; i < ALEA_PARTICLE_COUNT; i++) {
        p->scoped.energy_cutoff[i] = m->particle_energy_cutoff[i];
        p->scoped.secondary_state[i] = (int)m->particle_secondary_state[i];
    }
    p->scoped.detector_probabilities = alea_detector_probability_copy(
        m->detector_probabilities, m->detector_probability_count);
    if (m->detector_probability_count && !p->scoped.detector_probabilities) return -1;
    p->scoped.detector_probability_count = m->detector_probability_count;
    return 0;
}

alea_model_t* mcnp_model_to_alea_model(const mcnp_model_t* model) {
    if (!model || !model->sys) return NULL;
    alea_system_t* clone = alea_clone(model->sys);
    if (!clone) return NULL;
    alea_model_t* result = alea_model_adopt(clone);
    if (!result) { alea_destroy(clone); return NULL; }
    if (alea_model_set_name(result, model->name) ||
        alea_model_set_title(result, model->title) ||
        alea_model_set_comments(result, model->comments)) {
        alea_model_destroy(result);
        return NULL;
    }
    size_t count = alea_model_cell_metadata_count(result);
    if (count > model->cell_params_count) count = model->cell_params_count;
    for (size_t i = 0; i < count; i++) {
        if (params_to_metadata(&model->cell_params[i], alea_model_cell_metadata_mut(result, i))) {
            alea_model_destroy(result);
            return NULL;
        }
        if (alea_model_cell_set_name(result, i, model->cell_names[i])) {
            alea_model_destroy(result);
            return NULL;
        }
    }
    return result;
}

mcnp_model_t* mcnp_model_from_alea_model(const alea_model_t* model) {
    if (!model || !alea_model_system_const(model)) return NULL;
    alea_system_t* clone = alea_clone(alea_model_system_const(model));
    if (!clone) return NULL;
    mcnp_model_t* result = mcnp_model_wrap(clone);
    if (!result) { alea_destroy(clone); return NULL; }
    result->owns_sys = 1;
    if (mcnp_model_set_name(result, alea_model_name(model)) ||
        mcnp_model_set_title(result, alea_model_title(model)) ||
        mcnp_model_set_comments(result, alea_model_comments(model))) {
        mcnp_model_destroy(result);
        return NULL;
    }
    size_t count = result->cell_params_count;
    if (count > alea_model_cell_metadata_count(model))
        count = alea_model_cell_metadata_count(model);
    for (size_t i = 0; i < count; i++) {
        if (metadata_to_params(alea_model_cell_metadata(model, i), &result->cell_params[i])) {
            mcnp_model_destroy(result);
            return NULL;
        }
        const alea_model_cell_metadata_t* metadata = alea_model_cell_metadata(model, i);
        if (mcnp_model_cell_set_name(result, i, metadata ? metadata->name : NULL)) {
            mcnp_model_destroy(result);
            return NULL;
        }
    }
    return result;
}

mcnp_model_t* mcnp_load(const char* filename) {
    if (!filename) return NULL;

    double t_load0 = alea_monotonic_seconds();

    /* mcnp_convert_to_model creates the system, model, and populates
       cell params directly during conversion. */
    mcnp_model_t* model = mcnp_convert_to_model(filename);
    if (!model) return NULL;

    model = finalize_loaded_model(model);
    if (!model) {
        return NULL;
    }

    /* Build cell adjacency lazily on first use (raycasting, slicing, mesh
     * export). Eager adjacency construction is prohibitively expensive on
     * very large models and is not needed for pure conversion. */

    if (load_profile_enabled()) {
        double t_load1 = alea_monotonic_seconds();
        fprintf(stderr, "[alea-load-profile] %-28s %.6f s\n",
                "mcnp_load_total", t_load1 - t_load0);
    }

    return model;
}

mcnp_model_t* mcnp_load_string(const char* input, size_t len) {
    if (!input) return NULL;

    size_t actual_len = len > 0 ? len : strlen(input);
    mcnp_model_t* model = mcnp_convert_buffer_to_model(input, actual_len, "<memory>");
    return finalize_loaded_model(model);
}

int mcnp_export(const mcnp_model_t* model, const char* filename) {
    if (!model || !filename) return -1;

    FILE* f = fopen(filename, "w");
    if (!f) {
        alea_set_error_detail(ALEA_ERR_FILE_WRITE,
                              "failed to open MCNP export file '%s': %s",
                              filename, strerror(errno));
        return -1;
    }

    int ret = mcnp_export_stream(model, f);
    fclose(f);
    return ret;
}

/* Forward declaration for the MCNP export function */
extern int export_mcnp(alea_system_t* sys, export_context_t* ctx);

int mcnp_export_stream(const mcnp_model_t* model, FILE* out) {
    if (!model || !out || !model->sys) return -1;

    const mcnp_export_config_t* cfg = &model->export_config;
    alea_system_t* sys = model->sys;

    export_context_t* ctx = export_context_create(
        ALEA_EXPORT_FORMAT_MCNP,
        (alea_surface_emit_policy_t)cfg->surface_policy,
        out, true,
        alea_next_synthetic_surface_id(sys),
        sys->config.universe_depth,
        sys->config.fill_depth);
    if (!ctx) return -1;

    ctx->trcl_export_mode = (alea_trcl_export_mode_t)cfg->trcl_mode;
    ctx->transform_export_mode = (alea_transform_export_mode_t)cfg->transform_mode;
    ctx->mcnp_max_col = cfg->mcnp_max_col;
    ctx->mcnp_cont_indent = cfg->mcnp_cont_indent;
    ctx->module_data = model;

    int ret = export_mcnp(sys, ctx);
    export_context_destroy(ctx);
    return ret;
}

int mcnp_export_system(alea_system_t* sys, const char* filename) {
    if (!sys || !filename) return -1;

    FILE* f = fopen(filename, "w");
    if (!f) {
        alea_set_error_detail(ALEA_ERR_FILE_WRITE,
                              "failed to open MCNP export file '%s': %s",
                              filename, strerror(errno));
        return -1;
    }

    int ret = mcnp_export_system_stream(sys, f);
    fclose(f);
    return ret;
}

int mcnp_export_system_stream(alea_system_t* sys, FILE* out) {
    if (!sys || !out) return -1;

    export_context_t* ctx = export_context_create(
        ALEA_EXPORT_FORMAT_MCNP,
        ALEA_EMIT_MACROBODY,
        out, sys->config.dedup,
        alea_next_synthetic_surface_id(sys),
        sys->config.universe_depth,
        sys->config.fill_depth);
    if (!ctx) return -1;

    int ret = export_mcnp(sys, ctx);
    export_context_destroy(ctx);
    return ret;
}
