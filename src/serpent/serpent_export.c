// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/**
 * @file serpent_export.c
 * @brief Export CSG geometry to Serpent input format
 */

#include "serpent_export.h"
#include "alea_serpent.h"
#include "core/alea_system.h"
#include "core/alea_export.h"
#include "core/alea_materials.h"
#include "util/str_builder.h"
#include "util/alea_bitset.h"
#include "util/alea_log.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define SERPENT_EXPR_BUF 4096

typedef struct {
    int source_material_id;
    double density;
    bool is_mass_density;
    int serpent_material_id;
} serpent_mat_entry_t;

ALEA_VEC_DEFINE(serpent_mat_vec, serpent_mat_entry_t);

typedef struct {
    serpent_mat_vec_t entries;
    int next_id;
} serpent_mat_map_t;

typedef enum {
    SERPENT_CTX_TOP,
    SERPENT_CTX_UNION,
    SERPENT_CTX_INTERSECTION
} serpent_expr_context_t;

static const alea_material_t* find_material_by_id(const alea_system_t* sys, int material_id) {
    for (size_t i = 0; i < alea_vec_count(&sys->materials); i++) {
        if (sys->materials.data[i].material_id == material_id) return &sys->materials.data[i];
    }
    return NULL;
}

static const alea_mixture_t* find_mixture_by_id(const alea_system_t* sys, int mixture_id) {
    for (size_t i = 0; i < alea_vec_count(&sys->mixtures); i++) {
        if (sys->mixtures.data[i].mixture_id == mixture_id) return &sys->mixtures.data[i];
    }
    return NULL;
}

static int max_material_like_id(const alea_system_t* sys) {
    int max_id = 0;
    for (size_t i = 0; i < alea_vec_count(&sys->materials); i++) {
        if (sys->materials.data[i].material_id > max_id) max_id = sys->materials.data[i].material_id;
    }
    for (size_t i = 0; i < alea_vec_count(&sys->mixtures); i++) {
        if (sys->mixtures.data[i].mixture_id > max_id) max_id = sys->mixtures.data[i].mixture_id;
    }
    return max_id;
}

static void serpent_mat_map_init(serpent_mat_map_t* map, int start_id) {
    alea_vec_init(&map->entries);
    map->next_id = start_id;
}

static void serpent_mat_map_free(serpent_mat_map_t* map) {
    alea_vec_free(&map->entries);
}

static int serpent_mat_map_get(serpent_mat_map_t* map, int material_id,
                               double density, bool is_mass_density) {
    if (material_id <= 0) return 0;

    double rho = fabs(density);
    for (size_t i = 0; i < map->entries.count; i++) {
        serpent_mat_entry_t* e = &map->entries.data[i];
        if (e->source_material_id == material_id &&
            e->is_mass_density == is_mass_density &&
            fabs(e->density - rho) < 1e-10) {
            return e->serpent_material_id;
        }
    }

    serpent_mat_entry_t* e = alea_vec_push_uninit(&map->entries, serpent_mat_entry_t);
    if (!e) return material_id;
    e->source_material_id = material_id;
    e->density = rho;
    e->is_mass_density = is_mass_density;

    bool first_for_material = true;
    for (size_t i = 0; i + 1 < map->entries.count; i++) {
        if (map->entries.data[i].source_material_id == material_id) {
            first_for_material = false;
            break;
        }
    }
    e->serpent_material_id = first_for_material ? material_id : map->next_id++;
    return e->serpent_material_id;
}

static void build_material_map(const alea_system_t* sys, serpent_mat_map_t* map) {
    serpent_mat_map_init(map, max_material_like_id(sys) + 1000);
    for (size_t i = 0; i < alea_vec_count(&sys->cells); i++) {
        const alea_cell_entry_t* cell = &sys->cells.data[i];
        if (cell->material_id > 0) {
            serpent_mat_map_get(map, cell->material_id, cell->density, cell->is_mass_density);
        }
    }
}

static int canonical_surface_id(export_context_t* ctx, const alea_node_t* node) {
    int surf_id = node->primitive.mc_surface_id;
    if (ctx && ctx->deduplicate && ctx->prim_to_surface) {
        uint32_t prim_id = node->primitive.primitive_id;
        if (prim_id < ctx->prim_to_surface_size && ctx->prim_to_surface[prim_id] >= 0) {
            surf_id = ctx->prim_to_surface[prim_id];
        }
    }
    return surf_id;
}

static int effective_sense(export_context_t* ctx, const alea_node_t* node) {
    int sense = node->primitive.sense;
    if (ctx && ctx->deduplicate && ctx->prim_to_surface_inverted) {
        uint32_t prim_id = node->primitive.primitive_id;
        int8_t surface_inverted = 0;
        if (prim_id < ctx->prim_to_surface_size) {
            surface_inverted = ctx->prim_to_surface_inverted[prim_id];
        }
        if (node->primitive.inverted != surface_inverted) sense = -sense;
    }
    return sense;
}

static bool tree_to_serpent_region(const alea_system_t* sys,
                                   export_context_t* ctx,
                                   uint32_t node_id,
                                   str_builder_t* sb,
                                   serpent_expr_context_t parent_ctx) {
    if (node_id == ALEA_NODE_ID_INVALID || node_id >= alea_vec_count(&sys->nodes)) return false;

    const alea_node_t* node = &sys->nodes.data[node_id];
    alea_operation_t op = ALEA_GET_OPERATION(node);

    switch (op) {
        case ALEA_OP_PRIMITIVE: {
            int surf_id = canonical_surface_id(ctx, node);
            int sense = effective_sense(ctx, node);
            if (sense < 0) str_builder_putc(sb, '-');
            str_builder_int(sb, surf_id);
            return !str_builder_error(sb);
        }

        case ALEA_OP_COMPLEMENT: {
            uint32_t left = node->operation.left;
            if (left >= alea_vec_count(&sys->nodes)) return false;
            const alea_node_t* child = &sys->nodes.data[left];
            if (ALEA_GET_OPERATION(child) == ALEA_OP_PRIMITIVE) {
                int surf_id = canonical_surface_id(ctx, child);
                int sense = effective_sense(ctx, child);
                if (sense > 0) str_builder_putc(sb, '-');
                str_builder_int(sb, surf_id);
                return !str_builder_error(sb);
            }
            str_builder_putc(sb, '#');
            str_builder_putc(sb, '(');
            if (!tree_to_serpent_region(sys, ctx, left, sb, SERPENT_CTX_TOP)) return false;
            str_builder_putc(sb, ')');
            return !str_builder_error(sb);
        }

        case ALEA_OP_UNION: {
            uint32_t left = node->operation.left;
            uint32_t right = node->operation.right;
            if (left >= alea_vec_count(&sys->nodes) || right >= alea_vec_count(&sys->nodes)) return false;
            bool need_parens = (parent_ctx == SERPENT_CTX_INTERSECTION);
            if (need_parens) str_builder_putc(sb, '(');
            if (!tree_to_serpent_region(sys, ctx, left, sb, SERPENT_CTX_UNION)) return false;
            str_builder_puts(sb, " : ");
            if (!tree_to_serpent_region(sys, ctx, right, sb, SERPENT_CTX_UNION)) return false;
            if (need_parens) str_builder_putc(sb, ')');
            return !str_builder_error(sb);
        }

        case ALEA_OP_INTERSECTION: {
            uint32_t left = node->operation.left;
            uint32_t right = node->operation.right;
            if (left >= alea_vec_count(&sys->nodes) || right >= alea_vec_count(&sys->nodes)) return false;
            bool need_parens = (parent_ctx == SERPENT_CTX_UNION);
            if (need_parens) str_builder_putc(sb, '(');
            if (!tree_to_serpent_region(sys, ctx, left, sb, SERPENT_CTX_INTERSECTION)) return false;
            str_builder_putc(sb, ' ');
            if (!tree_to_serpent_region(sys, ctx, right, sb, SERPENT_CTX_INTERSECTION)) return false;
            if (need_parens) str_builder_putc(sb, ')');
            return !str_builder_error(sb);
        }

        case ALEA_OP_DIFFERENCE: {
            uint32_t left = node->operation.left;
            uint32_t right = node->operation.right;
            if (left >= alea_vec_count(&sys->nodes) || right >= alea_vec_count(&sys->nodes)) return false;
            if (!tree_to_serpent_region(sys, ctx, left, sb, SERPENT_CTX_INTERSECTION)) return false;
            str_builder_puts(sb, " #(");
            if (!tree_to_serpent_region(sys, ctx, right, sb, SERPENT_CTX_TOP)) return false;
            str_builder_putc(sb, ')');
            return !str_builder_error(sb);
        }

        default:
            return false;
    }
}

static void write_comment_block(FILE* out, const char* text) {
    if (!text || !*text) return;
    const char* p = text;
    while (*p) {
        const char* nl = strchr(p, '\n');
        size_t len = nl ? (size_t)(nl - p) : strlen(p);
        if (len > 0) fprintf(out, "%% %.*s\n", (int)len, p);
        p = nl ? nl + 1 : p + len;
    }
}

static int write_serpent_surface(FILE* out, int surface_id,
                                 alea_primitive_type_t type,
                                 const alea_primitive_data_t* data,
                                 int8_t inverted) {
    fprintf(out, "surf %d ", surface_id);

    switch (type) {
        case ALEA_PRIMITIVE_PLANE: {
            double a = data->plane.a;
            double b = data->plane.b;
            double c = data->plane.c;
            double d = data->plane.d;
            if (inverted) { a = -a; b = -b; c = -c; d = -d; }
            if (fabs(b) < 1e-10 && fabs(c) < 1e-10 && fabs(a - 1.0) < 1e-10) {
                fprintf(out, "px %.16g\n", -d);
            } else if (fabs(a) < 1e-10 && fabs(c) < 1e-10 && fabs(b - 1.0) < 1e-10) {
                fprintf(out, "py %.16g\n", -d);
            } else if (fabs(a) < 1e-10 && fabs(b) < 1e-10 && fabs(c - 1.0) < 1e-10) {
                fprintf(out, "pz %.16g\n", -d);
            } else {
                fprintf(out, "plane %.16g %.16g %.16g %.16g\n", a, b, c, -d);
            }
            return 0;
        }
        case ALEA_PRIMITIVE_SPHERE:
            fprintf(out, "sph %.16g %.16g %.16g %.16g\n",
                    data->sphere.center_x, data->sphere.center_y,
                    data->sphere.center_z, data->sphere.radius);
            return 0;
        case ALEA_PRIMITIVE_CYLINDER_X:
            fprintf(out, "cylx %.16g %.16g %.16g\n",
                    data->cyl_x.center_y, data->cyl_x.center_z, data->cyl_x.radius);
            return 0;
        case ALEA_PRIMITIVE_CYLINDER_Y:
            fprintf(out, "cyly %.16g %.16g %.16g\n",
                    data->cyl_y.center_x, data->cyl_y.center_z, data->cyl_y.radius);
            return 0;
        case ALEA_PRIMITIVE_CYLINDER_Z:
            fprintf(out, "cylz %.16g %.16g %.16g\n",
                    data->cyl_z.center_x, data->cyl_z.center_y, data->cyl_z.radius);
            return 0;
        case ALEA_PRIMITIVE_CONE_X:
            if (data->cone_x.sheet_selection != 0) return -1;
            fprintf(out, "ckx %.16g %.16g %.16g %.16g\n",
                    data->cone_x.apex_x, data->cone_x.apex_y,
                    data->cone_x.apex_z, data->cone_x.tan_angle_sq);
            return 0;
        case ALEA_PRIMITIVE_CONE_Y:
            if (data->cone_y.sheet_selection != 0) return -1;
            fprintf(out, "cky %.16g %.16g %.16g %.16g\n",
                    data->cone_y.apex_x, data->cone_y.apex_y,
                    data->cone_y.apex_z, data->cone_y.tan_angle_sq);
            return 0;
        case ALEA_PRIMITIVE_CONE_Z:
            if (data->cone_z.sheet_selection != 0) return -1;
            fprintf(out, "ckz %.16g %.16g %.16g %.16g\n",
                    data->cone_z.apex_x, data->cone_z.apex_y,
                    data->cone_z.apex_z, data->cone_z.tan_angle_sq);
            return 0;
        case ALEA_PRIMITIVE_RPP:
            fprintf(out, "cuboid %.16g %.16g %.16g %.16g %.16g %.16g\n",
                    data->box.min_x, data->box.max_x, data->box.min_y,
                    data->box.max_y, data->box.min_z, data->box.max_z);
            return 0;
        case ALEA_PRIMITIVE_RCC:
            fprintf(out, "rcc %.16g %.16g %.16g %.16g %.16g %.16g %.16g\n",
                    data->rcc.base_x, data->rcc.base_y, data->rcc.base_z,
                    data->rcc.height_x, data->rcc.height_y, data->rcc.height_z,
                    data->rcc.radius);
            return 0;
        case ALEA_PRIMITIVE_QUADRIC:
            fprintf(out, "quadratic %.16g %.16g %.16g %.16g %.16g %.16g %.16g %.16g %.16g %.16g\n",
                    data->quadric.coeffs[0], data->quadric.coeffs[1],
                    data->quadric.coeffs[2], data->quadric.coeffs[3],
                    data->quadric.coeffs[4], data->quadric.coeffs[5],
                    data->quadric.coeffs[6], data->quadric.coeffs[7],
                    data->quadric.coeffs[8], data->quadric.coeffs[9]);
            return 0;
        case ALEA_PRIMITIVE_TORUS_X:
        case ALEA_PRIMITIVE_TORUS_Y:
        case ALEA_PRIMITIVE_TORUS_Z: {
            const char* kind = type == ALEA_PRIMITIVE_TORUS_X ? "torx" :
                               type == ALEA_PRIMITIVE_TORUS_Y ? "tory" : "torz";
            fprintf(out, "%s %.16g %.16g %.16g %.16g %.16g %.16g\n",
                    kind, data->torus.center_x, data->torus.center_y,
                    data->torus.center_z, data->torus.major_radius,
                    data->torus.minor_radius,
                    data->torus.axial_semiwidth_B > 0 ? data->torus.axial_semiwidth_B : data->torus.minor_radius);
            return 0;
        }
        default:
            return -1;
    }
}

static void write_material_nuclides(FILE* out, const alea_material_t* mat, double scale) {
    const double sign = mat->is_weight_fraction ? -1.0 : 1.0;
    for (size_t i = 0; i < alea_vec_count(&mat->nuclides); i++) {
        const alea_nuclide_t* nuc = &mat->nuclides.data[i];
        fprintf(out, "    %d", nuc->zaid);
        if (nuc->library && nuc->library[0]) fprintf(out, ".%s", nuc->library);
        fprintf(out, " %.12g\n", sign * nuc->fraction * scale);
    }
    for (size_t i = 0; i < alea_vec_count(&mat->elements); i++) {
        const alea_element_comp_t* elem_comp = &mat->elements.data[i];
        fprintf(out, "    %d", elem_comp->atomic_number * 1000);
        if (elem_comp->library && elem_comp->library[0]) fprintf(out, ".%s", elem_comp->library);
        fprintf(out, " %.12g\n", sign * elem_comp->fraction * scale);
    }
}

static void write_serpent_materials(FILE* out, const alea_system_t* sys,
                                    const serpent_mat_map_t* map) {
    fprintf(out, "\n%% Materials\n");
    for (size_t i = 0; i < map->entries.count; i++) {
        const serpent_mat_entry_t* entry = &map->entries.data[i];
        const alea_material_t* mat = find_material_by_id(sys, entry->source_material_id);
        const alea_mixture_t* mix = NULL;

        double density = entry->is_mass_density ? -entry->density : entry->density;

        if (mat) {
            write_comment_block(out, mat->comments);
            fprintf(out, "mat m%d %.12g\n", entry->serpent_material_id, density);
            write_material_nuclides(out, mat, 1.0);
            continue;
        }

        mix = find_mixture_by_id(sys, entry->source_material_id);
        if (mix) {
            write_comment_block(out, mix->comments);
            fprintf(out, "mat m%d %.12g\n", entry->serpent_material_id, density);
            for (size_t c = 0; c < alea_vec_count(&mix->components); c++) {
                const alea_mixture_comp_t* comp = &mix->components.data[c];
                const alea_material_t* comp_mat = find_material_by_id(sys, comp->material_id);
                if (comp_mat) write_material_nuclides(out, comp_mat, comp->fraction);
            }
        }
    }
}

static int write_serpent_cells(FILE* out, const alea_system_t* sys,
                               export_context_t* ctx,
                               serpent_mat_map_t* mat_map) {
    fprintf(out, "\n%% Cells\n");
    for (size_t i = 0; i < alea_vec_count(&sys->cells); i++) {
        ALEA_CHECK_INTERRUPTED(-1);
        const alea_cell_entry_t* cell = &sys->cells.data[i];

        write_comment_block(out, cell->comments);
        fprintf(out, "cell %d %d ", cell->mc_cell_id, cell->universe_id);

        if (cell->fill_universe > 0) {
            fprintf(out, "fill %d", cell->fill_universe);
        } else if (cell->material_id > 0) {
            int mat_id = serpent_mat_map_get(mat_map, cell->material_id,
                                             cell->density, cell->is_mass_density);
            fprintf(out, "m%d", mat_id);
        } else {
            fprintf(out, "void");
        }

        if (cell->root_node_id != ALEA_NODE_ID_INVALID) {
            str_builder_t sb;
            str_builder_init(&sb, &ctx->arena, SERPENT_EXPR_BUF);
            if (!tree_to_serpent_region(sys, ctx, cell->root_node_id, &sb, SERPENT_CTX_TOP)) {
                ALEA_LOG_ERROR("Failed to convert cell %d to Serpent region", cell->mc_cell_id);
                return -1;
            }
            fprintf(out, " %s", str_builder_get(&sb));
        }

        if (cell->inline_comment && *cell->inline_comment) {
            fprintf(out, " %% %s", cell->inline_comment);
        }
        fprintf(out, "\n");
        ctx->cells_written++;
    }
    return 0;
}

static int write_serpent_surfaces(FILE* out, alea_system_t* sys, export_context_t* ctx) {
    fprintf(out, "\n%% Surfaces\n");

    if (ctx->deduplicate && ctx->prim_to_surface) {
        alea_bitset_t prim_written = alea_bitset_create(ctx->prim_to_surface_size);
        for (size_t i = 0; i < alea_vec_count(&sys->surfaces); i++) {
            ALEA_CHECK_INTERRUPTED(-1);
            const alea_surface_entry_t* surface = &sys->surfaces.data[i];
            uint32_t prim_id = surface->primitive_id;
            if (prim_id < ctx->prim_to_surface_size &&
                ctx->prim_to_surface[prim_id] >= 0 &&
                surface->mc_surface_id != ctx->prim_to_surface[prim_id]) {
                continue;
            }
            if (alea_bitset_test(&prim_written, prim_id)) continue;
            const alea_primitive_entry_t* prim = &sys->primitives.data[prim_id];
            int8_t inverted = 0;
            if (surface->pos_node < alea_vec_count(&sys->nodes)) {
                inverted = sys->nodes.data[surface->pos_node].primitive.inverted;
            }
            if (write_serpent_surface(out, surface->mc_surface_id, prim->type, &prim->data, inverted) < 0) {
                alea_bitset_destroy(&prim_written);
                alea_set_error_detail(ALEA_ERR_UNSUPPORTED_SURFACE,
                                      "Serpent export does not support primitive type %d", prim->type);
                return -1;
            }
            ctx->surfaces_written++;
            alea_bitset_set(&prim_written, prim_id);
        }
        alea_bitset_destroy(&prim_written);
        return 0;
    }

    for (size_t i = 0; i < alea_vec_count(&sys->surfaces); i++) {
        ALEA_CHECK_INTERRUPTED(-1);
        const alea_surface_entry_t* surface = &sys->surfaces.data[i];
        const alea_primitive_entry_t* prim = &sys->primitives.data[surface->primitive_id];
        int8_t inverted = 0;
        if (surface->pos_node < alea_vec_count(&sys->nodes)) {
            inverted = sys->nodes.data[surface->pos_node].primitive.inverted;
        }
        if (write_serpent_surface(out, surface->mc_surface_id, prim->type, &prim->data, inverted) < 0) {
            alea_set_error_detail(ALEA_ERR_UNSUPPORTED_SURFACE,
                                  "Serpent export does not support primitive type %d", prim->type);
            return -1;
        }
        ctx->surfaces_written++;
    }
    return 0;
}

int export_serpent(alea_system_t* sys, export_context_t* ctx) {
    if (!sys || !ctx || !ctx->out) return -1;

    serpent_mat_map_t mat_map;
    build_material_map(sys, &mat_map);

    fprintf(ctx->out, "%% Serpent input exported by Alea\n");

    alea_assign_missing_surface_ids(sys, ctx);
    if (ctx->deduplicate) alea_build_canonical_surface_map(ctx, sys);

    if (write_serpent_cells(ctx->out, sys, ctx, &mat_map) < 0) {
        serpent_mat_map_free(&mat_map);
        return -1;
    }

    if (write_serpent_surfaces(ctx->out, sys, ctx) < 0) {
        serpent_mat_map_free(&mat_map);
        return -1;
    }

    write_serpent_materials(ctx->out, sys, &mat_map);
    serpent_mat_map_free(&mat_map);
    return 0;
}

int serpent_export_system(alea_system_t* sys, const char* filename) {
    if (!sys || !filename) return -1;
    FILE* f = fopen(filename, "w");
    if (!f) return -1;
    int ret = serpent_export_system_stream(sys, f);
    fclose(f);
    return ret;
}

int serpent_export_system_stream(alea_system_t* sys, FILE* out) {
    if (!sys || !out) return -1;
    export_context_t* ctx = export_context_create(
        ALEA_EXPORT_FORMAT_SERPENT,
        ALEA_EMIT_MACROBODY,
        out,
        sys->config.dedup,
        alea_next_synthetic_surface_id(sys),
        sys->config.universe_depth,
        sys->config.fill_depth);
    if (!ctx) return -1;

    int ret = export_serpent(sys, ctx);
    export_context_destroy(ctx);
    return ret;
}
