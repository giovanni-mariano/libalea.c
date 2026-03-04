// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/**
 * @file openmc_model.c
 * @brief OpenMC model implementation
 *
 * Manages the openmc_model_t lifecycle and convenience export functions.
 */

#include "openmc_model.h"
#include "alea.h"
#include "core/alea_system.h"
#include "core/alea_export.h"
#include "util/alea_log.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

/* Forward declaration for the OpenMC export function */
extern int export_openmc(const alea_system_t* sys, export_context_t* ctx);

/* ============================================================================
 * DEDUP REPORT (shared utility, moved from alea_public_api.c)
 * ============================================================================ */

static void write_dedup_report(const alea_system_t* sys, const export_context_t* ctx,
                               const char* xml_filename) {
    if (!ctx->deduplicate || !ctx->prim_to_surface || ctx->dedup_hits == 0) return;

    /* Build report filename: replace .xml with .dedup, or append .dedup */
    char report_path[1024];
    size_t len = strlen(xml_filename);
    if (len > 4 && strcmp(xml_filename + len - 4, ".xml") == 0) {
        snprintf(report_path, sizeof(report_path), "%.*s.dedup", (int)(len - 4), xml_filename);
    } else {
        snprintf(report_path, sizeof(report_path), "%s.dedup", xml_filename);
    }

    FILE* rf = fopen(report_path, "w");
    if (!rf) return;

    fprintf(rf, "# Surface deduplication report\n");
    fprintf(rf, "# original_surface -> canonical_surface  sign\n");
    fprintf(rf, "#\n");

    size_t dedup_count = 0;
    for (size_t i = 0; i < alea_vec_count(&sys->surfaces); i++) {
        const alea_surface_entry_t* surf = &sys->surfaces.data[i];
        uint32_t prim_id = surf->primitive_id;
        int orig_id = surf->mcnp_surface_id;

        if (prim_id >= ctx->prim_to_surface_size) continue;
        int canon_id = ctx->prim_to_surface[prim_id];
        if (canon_id < 0 || canon_id == orig_id) continue;

        int8_t orig_inv = 0;
        if (surf->pos_node < alea_vec_count(&sys->nodes)) {
            orig_inv = sys->nodes.data[surf->pos_node].primitive.inverted;
        }
        int8_t canon_inv = ctx->prim_to_surface_inverted[prim_id];
        const char* sign = (orig_inv == canon_inv) ? "+" : "-";

        fprintf(rf, "%d -> %d  %s\n", orig_id, canon_id, sign);
        dedup_count++;
    }

    fprintf(rf, "#\n");
    fprintf(rf, "# %zu surfaces deduplicated\n", dedup_count);
    fclose(rf);

    ALEA_LOG_INFO("Deduplication report written to %s (%zu entries)", report_path, dedup_count);
}

/* ============================================================================
 * INTERNAL: Create export context for OpenMC
 * ============================================================================ */

static export_context_t* create_openmc_export_ctx(const alea_system_t* sys, FILE* out) {
    int next_id = alea_next_synthetic_surface_id(sys);
    /* OpenMC always uses ALEA_EMIT_SURFACES for primitive decomposition */
    export_context_t* ctx = export_context_create(
        ALEA_EXPORT_FORMAT_OPENMC,
        ALEA_EMIT_SURFACES,
        out,
        sys->config.dedup,
        next_id,
        sys->config.universe_depth,
        sys->config.fill_depth
    );
    return ctx;
}

/* ============================================================================
 * PUBLIC API
 * ============================================================================ */

void openmc_model_destroy(openmc_model_t* model) {
    if (!model) return;
    if (model->owns_sys && model->sys) {
        alea_destroy(model->sys);
    }
    free(model);
}

openmc_model_t* openmc_model_wrap(alea_system_t* sys) {
    if (!sys) return NULL;

    openmc_model_t* model = calloc(1, sizeof(openmc_model_t));
    if (!model) return NULL;

    model->sys = sys;
    model->owns_sys = 0;
    return model;
}

openmc_model_t* openmc_load(const char* filename) {
    if (!filename) return NULL;

    openmc_model_t* model = openmc_convert_to_model(filename);
    if (!model) return NULL;

    model->sys->source = ALEA_SOURCE_OPENMC;
    return model;
}

openmc_model_t* openmc_load_string(const char* input, size_t len) {
    if (!input) return NULL;

    openmc_model_t* model = openmc_convert_string_to_model(input, len);
    if (!model) return NULL;

    model->sys->source = ALEA_SOURCE_OPENMC;
    return model;
}

int openmc_export(const openmc_model_t* model, const char* filename) {
    if (!model || !filename || !model->sys) return -1;

    FILE* f = fopen(filename, "w");
    if (!f) return -1;

    export_context_t* ctx = create_openmc_export_ctx(model->sys, f);
    if (!ctx) {
        fclose(f);
        return -1;
    }

    int ret = export_openmc(model->sys, ctx);

    if (ret == 0) {
        write_dedup_report(model->sys, ctx, filename);
    }

    export_context_destroy(ctx);
    fclose(f);
    return ret;
}

int openmc_export_stream(const openmc_model_t* model, FILE* out) {
    if (!model || !out || !model->sys) return -1;

    export_context_t* ctx = create_openmc_export_ctx(model->sys, out);
    if (!ctx) return -1;

    int ret = export_openmc(model->sys, ctx);
    export_context_destroy(ctx);
    return ret;
}

int openmc_export_system(const alea_system_t* sys, const char* filename) {
    if (!sys || !filename) return -1;

    FILE* f = fopen(filename, "w");
    if (!f) return -1;

    export_context_t* ctx = create_openmc_export_ctx(sys, f);
    if (!ctx) {
        fclose(f);
        return -1;
    }

    int ret = export_openmc(sys, ctx);

    if (ret == 0) {
        write_dedup_report(sys, ctx, filename);
    }

    export_context_destroy(ctx);
    fclose(f);
    return ret;
}

int openmc_export_system_stream(const alea_system_t* sys, FILE* out) {
    if (!sys || !out) return -1;

    export_context_t* ctx = create_openmc_export_ctx(sys, out);
    if (!ctx) return -1;

    int ret = export_openmc(sys, ctx);
    export_context_destroy(ctx);
    return ret;
}
