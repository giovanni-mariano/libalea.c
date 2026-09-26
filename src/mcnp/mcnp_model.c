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
        model->cell_params[dst_index] = model->cell_params[src_index];
    }
}

static void on_cell_removed_cb(void* ud, size_t index) {
    mcnp_model_t* model = (mcnp_model_t*)ud;
    if (!model || index >= model->cell_params_count) return;
    const size_t trailing = model->cell_params_count - index - 1;
    if (trailing > 0) {
        memmove(&model->cell_params[index], &model->cell_params[index + 1],
                trailing * sizeof(*model->cell_params));
    }
    model->cell_params_count--;
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

static void params_to_metadata(const mcnp_cell_params_t* p,
                               alea_model_cell_metadata_t* m) {
    if (!p || !m) return;
    m->importance_neutron = p->imp_n;
    m->importance_photon = p->imp_p;
    m->importance_electron = p->imp_e;
    m->has_importance_neutron = p->has_imp_n;
    m->has_importance_photon = p->has_imp_p;
    m->has_importance_electron = p->has_imp_e;
    if (p->has_vol) { m->user_volume = p->vol; m->parameter_flags |= ALEA_CELL_PARAM_VOLUME; }
    if (p->has_pwt) { m->photon_weight = p->pwt; m->parameter_flags |= ALEA_CELL_PARAM_PWT; }
    if (p->has_nonu) { m->fission_turnoff = p->nonu; m->parameter_flags |= ALEA_CELL_PARAM_NONU; }
    if (p->has_pd) { m->detector_contribution = p->pd; m->parameter_flags |= ALEA_CELL_PARAM_PD; }
    if (p->has_elpt) { m->energy_cutoff = p->elpt; m->parameter_flags |= ALEA_CELL_PARAM_ELPT; }
    if (p->has_unc) { m->uncollided_secondaries = p->unc; m->parameter_flags |= ALEA_CELL_PARAM_UNC; }
    if (p->has_bflcl) { m->magnetic_field = p->bflcl; m->parameter_flags |= ALEA_CELL_PARAM_BFLCL; }
}

static void metadata_to_params(const alea_model_cell_metadata_t* m,
                               mcnp_cell_params_t* p) {
    if (!m || !p) return;
    p->imp_n = m->importance_neutron;
    p->imp_p = m->importance_photon;
    p->imp_e = m->importance_electron;
    p->has_imp_n = m->has_importance_neutron;
    p->has_imp_p = m->has_importance_photon;
    p->has_imp_e = m->has_importance_electron;
    if (m->parameter_flags & ALEA_CELL_PARAM_VOLUME) { p->vol = m->user_volume; p->has_vol = 1; }
    if (m->parameter_flags & ALEA_CELL_PARAM_PWT) { p->pwt = m->photon_weight; p->has_pwt = 1; }
    if (m->parameter_flags & ALEA_CELL_PARAM_NONU) { p->nonu = m->fission_turnoff; p->has_nonu = 1; }
    if (m->parameter_flags & ALEA_CELL_PARAM_PD) { p->pd = m->detector_contribution; p->has_pd = 1; }
    if (m->parameter_flags & ALEA_CELL_PARAM_ELPT) { p->elpt = m->energy_cutoff; p->has_elpt = 1; }
    if (m->parameter_flags & ALEA_CELL_PARAM_UNC) { p->unc = m->uncollided_secondaries; p->has_unc = 1; }
    if (m->parameter_flags & ALEA_CELL_PARAM_BFLCL) { p->bflcl = m->magnetic_field; p->has_bflcl = 1; }
}

alea_model_t* mcnp_model_to_alea_model(const mcnp_model_t* model) {
    if (!model || !model->sys) return NULL;
    alea_system_t* clone = alea_clone(model->sys);
    if (!clone) return NULL;
    alea_model_t* result = alea_model_adopt(clone);
    if (!result) { alea_destroy(clone); return NULL; }
    size_t count = alea_model_cell_metadata_count(result);
    if (count > model->cell_params_count) count = model->cell_params_count;
    for (size_t i = 0; i < count; i++)
        params_to_metadata(&model->cell_params[i], alea_model_cell_metadata_mut(result, i));
    return result;
}

mcnp_model_t* mcnp_model_from_alea_model(const alea_model_t* model) {
    if (!model || !alea_model_system_const(model)) return NULL;
    alea_system_t* clone = alea_clone(alea_model_system_const(model));
    if (!clone) return NULL;
    mcnp_model_t* result = mcnp_model_wrap(clone);
    if (!result) { alea_destroy(clone); return NULL; }
    result->owns_sys = 1;
    size_t count = result->cell_params_count;
    if (count > alea_model_cell_metadata_count(model))
        count = alea_model_cell_metadata_count(model);
    for (size_t i = 0; i < count; i++)
        metadata_to_params(alea_model_cell_metadata(model, i), &result->cell_params[i]);
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
