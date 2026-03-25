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
#include "core/alea_spatial.h"
#include "core/alea_export.h"
#include "util/alea_log.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "util/compat.h"

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

    return model;
}

mcnp_model_t* mcnp_load(const char* filename) {
    if (!filename) return NULL;

    /* mcnp_convert_to_model creates the system, model, and populates
       cell params directly during conversion. */
    mcnp_model_t* model = mcnp_convert_to_model(filename);
    if (!model) return NULL;

    if (alea_validate_cell_ids(model->sys) < 0) {
        alea_set_error_detail(ALEA_ERR_INVALID_ID, "Duplicate cell IDs");
        mcnp_model_destroy(model);
        return NULL;
    }

    model->sys->source = ALEA_SOURCE_MCNP;

    /* Build cell adjacency eagerly; spatial index is built lazily on first
     * query (raycasting, slicing, rendering) — not needed for pure conversion. */
    alea_build_cell_adjacency(model->sys);

    return model;
}

mcnp_model_t* mcnp_load_string(const char* input, size_t len) {
    if (!input) return NULL;

    size_t actual_len = len > 0 ? len : strlen(input);

    /* Write to temp file and use mcnp_load */
    char tmppath[256];
    FILE* f = alea_tmpfile(tmppath);
    if (!f) return NULL;

    fwrite(input, 1, actual_len, f);
    fclose(f);

    mcnp_model_t* model = mcnp_load(tmppath);
    remove(tmppath);

    return model;
}

int mcnp_export(const mcnp_model_t* model, const char* filename) {
    if (!model || !filename) return -1;

    FILE* f = fopen(filename, "w");
    if (!f) return -1;

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
    if (!f) return -1;

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
