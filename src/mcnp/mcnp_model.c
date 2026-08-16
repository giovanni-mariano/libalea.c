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
    (void)ud;
    (void)index;
    /* Cell removal is handled by compaction — the caller rebuilds the parallel
       array or we rely on count tracking.  For now, no-op. */
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

mcnp_model_t* mcnp_model_wrap(alea_system_t* sys) {
    if (!sys) return NULL;

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
