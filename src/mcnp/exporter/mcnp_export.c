// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/**
 * @file mcnp_export.c
 * @brief MCNP format export implementation
 *
 * MCNP-specific export functionality including:
 * - CSG tree to MCNP expression conversion
 * - Cell card writing with proper 80-column wrapping
 * - Surface card writing for all primitive types
 * - Material and transform card writing
 */

#include "mcnp_export.h"
#include "mcnp_str.h"
#include "mcnp/mcnp_model.h"
#include "core/alea_system.h"
#include "core/alea_export.h"
#include "core/alea_materials.h"
#include "core/alea_macrobody.h"
#include "primitives/primitive_desc.h"
#include "util/alea_log.h"
#include "util/alea_bitset.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include "util/math.h"

/* ============================================================================
 * INTERNAL CONSTANTS
 * ============================================================================ */

#define EXPR_BUF_SIZE 4096
#define SURFACE_BUF_SIZE 256

/* ============================================================================
 * HELPER FUNCTIONS
 * ============================================================================ */

/**
 * @brief Check if a node is a cell complement reference
 *
 * If the node is a ALEA_OP_COMPLEMENT that was created as a cell reference (#N),
 * returns the referenced cell ID. Otherwise returns -1.
 */
static int get_cell_complement_ref(const alea_system_t* sys, uint32_t node_id) {
    if (!sys) return -1;

    for (size_t i = 0; i < alea_vec_count(&sys->cell_refs); i++) {
        if (sys->cell_refs.data[i].placeholder_node == node_id) {
            return sys->cell_refs.data[i].referenced_cell_id;
        }
    }

    return -1;
}

/**
 * @brief Get the MCNP surface ID to use for a primitive node in cell expressions.
 *
 * When deduplication is enabled, returns the CANONICAL surface ID
 * (the one that was actually exported), not the original mc_surface_id
 * stored in the node.
 */
static int find_surface_id_for_node(export_context_t* ctx,
                                    const alea_system_t* sys,
                                    uint32_t node_id) {
    if (node_id >= alea_vec_count(&sys->nodes)) return -1;

    const alea_node_t* node = &sys->nodes.data[node_id];

    if (ALEA_GET_OPERATION(node) != ALEA_OP_PRIMITIVE) return -1;

    /* If deduplication enabled, use canonical surface ID based on primitive_id */
    if (ctx && ctx->deduplicate && ctx->prim_to_surface) {
        uint32_t prim_id = node->primitive.primitive_id;

        if (prim_id < ctx->prim_to_surface_size &&
            ctx->prim_to_surface[prim_id] >= 0) {
            return ctx->prim_to_surface[prim_id];
        }

        /* Primitive not in map (e.g. API-created or from macrobody expansion)
           — fall through to use the node's stored surface ID */
    }

    /* No dedup or not in map - use the ID stored in the node */
    return node->primitive.mc_surface_id;
}

/* ============================================================================
 * CSG TREE TO EXPRESSION
 * ============================================================================ */

/*
 * MCNP operator precedence (highest to lowest):
 *   1. Parentheses
 *   2. Complement (#)
 *   3. Intersection (space/juxtaposition)
 *   4. Union (:)
 *
 * We only need parens when a union appears inside an intersection context.
 */

typedef enum {
    CONTEXT_TOP,          /* Top level - no parens needed */
    CONTEXT_UNION,        /* Inside union - no parens needed for children */
    CONTEXT_INTERSECTION  /* Inside intersection - union children need parens */
} expr_context_t;

static bool tree_to_expr_recursive(const alea_system_t* sys,
                                   export_context_t* ctx,
                                   uint32_t node_id,
                                   mcnp_str_t* s,
                                   expr_context_t context) {
    /* Bounds checking */
    if (node_id == UINT32_MAX) return false;
    if (node_id >= alea_vec_count(&sys->nodes)) return false;

    const alea_node_t* node = &sys->nodes.data[node_id];
    alea_operation_t op = ALEA_GET_OPERATION(node);

    if (op == ALEA_OP_PRIMITIVE) {
        /* Look up by node_id or by primitive_id for new nodes */
        int surface_id = find_surface_id_for_node(ctx, sys, node_id);
        if (surface_id < 0) {
            ALEA_LOG_WARN("No surface ID for node %u (prim_id=%u)",
                    node_id, node->primitive.primitive_id);
            return false;
        }

        /* Compute effective sense relative to the EXPORTED surface coefficients.
           Surface export un-canonicalizes plane coefficients using the surface's
           pos_node inverted flag. Cell sense must be relative to that orientation.
           Formula: flip sense when node's inverted differs from the exported
           surface's inverted (XOR). In non-dedup mode, they always match
           (same surface), so effective_sense = sense. */
        int effective_sense = node->primitive.sense;
        if (ctx && ctx->deduplicate && ctx->prim_to_surface_inverted) {
            uint32_t prim_id = node->primitive.primitive_id;
            int8_t surface_inverted = 0;
            if (prim_id < ctx->prim_to_surface_size) {
                surface_inverted = ctx->prim_to_surface_inverted[prim_id];
            }
            if (node->primitive.inverted != surface_inverted) {
                effective_sense = -effective_sense;
            }
        }

        /* Check for decomposed macrobody */
        const macrobody_decomposition_t* decomp = NULL;
        if (ctx && ctx->surface_policy == ALEA_EMIT_SURFACES) {
            decomp = find_decomposition(ctx, surface_id);
        }

        if (decomp) {
            /* Emit decomposed expression */
            bool is_interior = (effective_sense < 0);

            if (!is_interior) {
                /* Exterior = complement of interior */
                mcnp_str_putc(s, '#');
            }

            mcnp_str_putc(s, '(');

            for (size_t i = 0; i < decomp->surface_count; i++) {
                if (i > 0) mcnp_str_putc(s, ' ');

                int surf = decomp->surfaces[i].surface_id;
                int8_t sense = decomp->surfaces[i].sense;

                if (sense < 0) surf = -surf;

                bool written = mcnp_str_int(s, surf);
                if (!written) return false;
            }

            mcnp_str_putc(s, ')');
            return true;
        }

        if (effective_sense < 0) {
            surface_id = -surface_id;
        }

        bool written = mcnp_str_int(s, surface_id);
        if (!written) return false;
        return true;
    }

    if (op == ALEA_OP_COMPLEMENT) {
        /* Check if this is a cell reference (#cell_id) */
        int cell_ref = get_cell_complement_ref(sys, node_id);
        if (cell_ref >= 0) {
            mcnp_str_putc(s, '#');
            bool written = mcnp_str_int(s, cell_ref);
            if (!written) return false;
            return true;
        }

        uint32_t left = node->operation.left;
        if (left >= alea_vec_count(&sys->nodes)) return false;

        const alea_node_t* child = &sys->nodes.data[left];
        alea_operation_t child_op = ALEA_GET_OPERATION(child);

        /* Optimization: if complement of single primitive, just flip the sign */
        if (child_op == ALEA_OP_PRIMITIVE) {
            int surface_id = find_surface_id_for_node(ctx, sys, left);
            if (surface_id < 0) {
                mcnp_str_putc(s, '#');
                if (!tree_to_expr_recursive(sys, ctx, left, s, CONTEXT_INTERSECTION))
                    return false;
                return true;
            }

            /* Compute effective sense relative to the exported surface
               (same logic as non-complement case, see comment there) */
            int effective_sense = child->primitive.sense;
            if (ctx && ctx->deduplicate && ctx->prim_to_surface_inverted) {
                uint32_t prim_id = child->primitive.primitive_id;
                int8_t surface_inverted = 0;
                if (prim_id < ctx->prim_to_surface_size) {
                    surface_inverted = ctx->prim_to_surface_inverted[prim_id];
                }
                if (child->primitive.inverted != surface_inverted) {
                    effective_sense = -effective_sense;
                }
            }

            /* Complement flips the sense */
            if (effective_sense > 0) {
                surface_id = -surface_id;
            }

            bool written = mcnp_str_int(s, surface_id);
            if (!written) return false;
            return true;
        }

        /* Complex expression: use #(...) */
        mcnp_str_putc(s, '#');
        mcnp_str_putc(s, '(');
        if (!tree_to_expr_recursive(sys, ctx, left, s, CONTEXT_TOP))
            return false;
        mcnp_str_putc(s, ')');
        return true;
    }

    if (op == ALEA_OP_UNION) {
        uint32_t left = node->operation.left;
        uint32_t right = node->operation.right;
        if (left >= alea_vec_count(&sys->nodes) || right >= alea_vec_count(&sys->nodes)) return false;

        /* Need parens if we're inside an intersection */
        bool need_parens = (context == CONTEXT_INTERSECTION);

        if (need_parens) mcnp_str_putc(s, '(');

        if (!tree_to_expr_recursive(sys, ctx, left, s, CONTEXT_UNION))
            return false;

        mcnp_str_putc(s, ':');

        if (!tree_to_expr_recursive(sys, ctx, right, s, CONTEXT_UNION))
            return false;

        if (need_parens) mcnp_str_putc(s, ')');
        return true;
    }

    if (op == ALEA_OP_INTERSECTION) {
        uint32_t left = node->operation.left;
        uint32_t right = node->operation.right;
        if (left >= alea_vec_count(&sys->nodes) || right >= alea_vec_count(&sys->nodes)) return false;

        /* Parenthesize intersection branches inside a union for readability */
        bool need_parens = (context == CONTEXT_UNION);

        if (need_parens) mcnp_str_putc(s, '(');

        if (!tree_to_expr_recursive(sys, ctx, left, s, CONTEXT_INTERSECTION))
            return false;

        mcnp_str_putc(s, ' ');

        if (!tree_to_expr_recursive(sys, ctx, right, s, CONTEXT_INTERSECTION))
            return false;

        if (need_parens) mcnp_str_putc(s, ')');
        return true;
    }

    if (op == ALEA_OP_DIFFERENCE) {
        uint32_t left = node->operation.left;
        uint32_t right = node->operation.right;
        if (left >= alea_vec_count(&sys->nodes) || right >= alea_vec_count(&sys->nodes)) return false;

        /* A - B = A #B */
        if (!tree_to_expr_recursive(sys, ctx, left, s, CONTEXT_INTERSECTION))
            return false;

        mcnp_str_putc(s, ' ');
        mcnp_str_putc(s, '#');

        const alea_node_t* right_node = &sys->nodes.data[right];
        alea_operation_t right_op = ALEA_GET_OPERATION(right_node);
        bool need_parens = (right_op != ALEA_OP_PRIMITIVE);

        if (need_parens) {
            mcnp_str_putc(s, '(');
            mcnp_str_putc(s, '(');
        }
        if (!tree_to_expr_recursive(sys, ctx, right, s, CONTEXT_INTERSECTION))
            return false;
        if (need_parens) {
            mcnp_str_putc(s, ')');
            mcnp_str_putc(s, ')');
        }

        return true;
    }

    ALEA_LOG_WARN("Unknown operation type %d at node %u", op, node_id);
    return false;
}

static int alea_tree_to_mcnp_expr(const alea_system_t* sys,
                                  export_context_t* ctx,
                                  uint32_t node_id,
                                  mcnp_str_t* s) {
    if (!sys) return -1;

    bool ok = tree_to_expr_recursive(sys, ctx, node_id, s, CONTEXT_TOP);
    if (!ok || mcnp_str_error(s)) return -1;

    return (int)s->sb.len;
}

/* ============================================================================
 * MCNP CELL WRITER
 * ============================================================================ */

/**
 * @brief Write a cell card with proper 80-column wrapping
 */
/* Per-cell flatten info computed once before the cell emit loop. See
 * build_flatten_info() for the detection rules. */
typedef struct {
    /* If skip is true, the entire cell is omitted from the output (it has been
     * folded into a LAT cell that fills directly into the cell's fill_universe). */
    bool skip;
    /* If override_fill_universe > 0, this LAT cell uses this universe instead of
     * its own lat_fill[0]. override_fill_transform (if non-zero) is the transform
     * to apply on the FILL. */
    int override_fill_universe;
    int override_fill_transform;
} flatten_cell_info_t;

/* For each surface id (sparse, keyed by mc_surface_id), an offset to apply to
 * the emitted primitive data. Used to position a flattened LAT cell's bounding
 * macrobody at [lower_left, lower_left + pitch] in the parent universe frame
 * without modifying the internal CSG. */
typedef struct {
    int surface_id;
    double off_x, off_y, off_z;
} surface_offset_t;

/* Shift primitive bounds by an offset when emitting a flattened LAT cell's
 * bounding macrobody. Only RPP is handled today — the synthetic LAT element
 * boxes that come out of build_lattice_element_tree are always RPPs. */
static void apply_flatten_surface_offset(int mc_surface_id,
                                         const surface_offset_t* offsets,
                                         size_t count,
                                         alea_primitive_type_t type,
                                         alea_primitive_data_t* data) {
    if (!offsets || count == 0 || type != ALEA_PRIMITIVE_RPP) return;
    for (size_t i = 0; i < count; i++) {
        if (offsets[i].surface_id != mc_surface_id) continue;
        data->box.min_x += offsets[i].off_x;
        data->box.max_x += offsets[i].off_x;
        data->box.min_y += offsets[i].off_y;
        data->box.max_y += offsets[i].off_y;
        data->box.min_z += offsets[i].off_z;
        data->box.max_z += offsets[i].off_z;
        return;
    }
}

static void write_mcnp_cell_line(FILE* out, const alea_cell_entry_t* cell,
                                  size_t cell_index, const char* expr,
                                  const alea_system_t* sys,
                                  export_context_t* ctx,
                                  const flatten_cell_info_t* flatten,
                                  const char* inline_comment) {

    mcnp_str_t s;
    mcnp_str_init(&s, &ctx->arena, EXPR_BUF_SIZE, ctx->mcnp_max_col, ctx->mcnp_cont_indent);

    /* Write cell ID and material */
    mcnp_str_int(&s, cell->mc_cell_id);
    mcnp_str_putc(&s, ' ');

    if (cell->material_id == 0) {
        mcnp_str_putc(&s, '0');
        mcnp_str_putc(&s, ' ');
    } else {
        mcnp_str_int(&s, cell->material_id);
        mcnp_str_putc(&s, ' ');
        mcnp_str_double(&s, cell->is_mass_density ? -cell->density : cell->density, 6);
        mcnp_str_putc(&s, ' ');
    }

    /* Write geometry expression */
    mcnp_str_puts(&s, expr);

    mcnp_str_newline(&s);

    /* Read MCNP-specific params from model if available */
    const mcnp_model_t* model = (const mcnp_model_t*)ctx->module_data;
    const mcnp_cell_params_t* mp = model ? mcnp_cell_params_const(model, cell_index) : NULL;

    /* IMP:N - neutron importance (default 1.0) */
    mcnp_str_puts(&s, " IMP:N=");
    mcnp_str_double(&s, (mp && mp->has_imp_n) ? mp->imp_n : 1.0, 4);

    /* IMP:P - photon importance (default 1.0) */
    mcnp_str_puts(&s, " IMP:P=");
    mcnp_str_double(&s, (mp && mp->has_imp_p) ? mp->imp_p : 1.0, 4);

    /* IMP:E - electron importance */
    if (mp && mp->has_imp_e) {
        mcnp_str_puts(&s, " IMP:E=");
        mcnp_str_double(&s, mp->imp_e, 4);
    }

    /* Simple keyword=value params (VOL, TMP, PWT, NONU, PD, ELPT, UNC, BFLCL) */
    #define MCNP_EXPORT_double(s, val, prec) mcnp_str_double(s, val, prec)
    #define MCNP_EXPORT_int(s, val, prec)    mcnp_str_int(s, val)
    #define X_EXPORT(name, type, kw, prec) \
        if (mp && mp->has_##name) { \
            mcnp_str_puts(&s, " " kw "="); \
            MCNP_EXPORT_##type(&s, mp->name, prec); \
        }
    MCNP_CELL_SIMPLE_PARAMS(X_EXPORT)
    #undef X_EXPORT
    #undef MCNP_EXPORT_double
    #undef MCNP_EXPORT_int

    /* U= - universe membership */
    if (cell->universe_id != 0) {
        mcnp_str_puts(&s, " U=");
        mcnp_str_int(&s, cell->universe_id);
    }

    /* FILL= or *FILL= - fill with universe */
    if (cell->fill_universe != 0) {
        const alea_transform_t* tr = NULL;
        if (cell->fill_transform != 0) {
            tr = alea_get_transform(sys, cell->fill_transform);
        }

        int use_degrees = tr ? tr->degrees : (mp ? mp->fill_transform_degrees : 0);

        if (use_degrees) {
            mcnp_str_puts(&s, " *FILL=");
        } else {
            mcnp_str_puts(&s, " FILL=");
        }
        mcnp_str_int(&s, cell->fill_universe);

        if (tr && tr->value_count > 0) {
            bool emit_inline = (ctx->transform_export_mode == ALEA_TR_EXPORT_INLINE) ||
                               (ctx->transform_export_mode == ALEA_TR_EXPORT_ORIGINAL && tr->from_inline);
            if (emit_inline) {
                mcnp_str_puts(&s, " (");
                for (int i = 0; i < tr->value_count; i++) {
                    if (i > 0) mcnp_str_putc(&s, ' ');
                    /* tr->data contains original values (degrees or cosines) */
                    mcnp_str_double(&s, tr->data[i], 10);
                }
                mcnp_str_putc(&s, ')');
            } else {
                mcnp_str_puts(&s, " (");
                mcnp_str_int(&s, cell->fill_transform);
                mcnp_str_putc(&s, ')');
            }
        } else if (cell->fill_transform != 0) {
            mcnp_str_puts(&s, " (");
            mcnp_str_int(&s, cell->fill_transform);
            mcnp_str_putc(&s, ')');
        }
    }

    /* LAT= - lattice type */
    if (cell->lat_type != 0) {
        mcnp_str_puts(&s, " LAT=");
        mcnp_str_int(&s, cell->lat_type);
    }

    /* Lattice FILL: emit simple "FILL=N" for 1x1x1 lattices (single element,
     * single universe), and the explicit "FILL=imin:imax ... U U U" array
     * form otherwise. Both are valid MCNP; the simple form is idiomatic for
     * degenerate 1x1x1 lattices that act as positioning wrappers. */
    if (cell->lat_fill && cell->lat_fill_count > 0) {
        int ni = cell->lat_fill_dims[1] - cell->lat_fill_dims[0] + 1;
        int nj = cell->lat_fill_dims[3] - cell->lat_fill_dims[2] + 1;
        int nk = cell->lat_fill_dims[5] - cell->lat_fill_dims[4] + 1;

        if (cell->lat_fill_count == 1 && ni == 1 && nj == 1 && nk == 1) {
            /* Flatten override: when the lat_fill[0] universe is a passthrough,
             * the pre-pass redirects FILL to its target. */
            int fill_univ = (flatten && flatten->override_fill_universe > 0)
                          ? flatten->override_fill_universe : cell->lat_fill[0];
            int fill_xform = (flatten && flatten->override_fill_universe > 0)
                           ? flatten->override_fill_transform : 0;

            mcnp_str_puts(&s, " FILL=");
            mcnp_str_int(&s, fill_univ);
            if (fill_xform != 0) {
                const alea_transform_t* tr = alea_get_transform(sys, fill_xform);
                if (tr && tr->value_count > 0 &&
                    ((ctx->transform_export_mode == ALEA_TR_EXPORT_INLINE) ||
                     (ctx->transform_export_mode == ALEA_TR_EXPORT_ORIGINAL && tr->from_inline))) {
                    mcnp_str_puts(&s, " (");
                    for (int i = 0; i < tr->value_count; i++) {
                        if (i > 0) mcnp_str_putc(&s, ' ');
                        mcnp_str_double(&s, tr->data[i], 10);
                    }
                    mcnp_str_putc(&s, ')');
                } else {
                    mcnp_str_puts(&s, " (");
                    mcnp_str_int(&s, fill_xform);
                    mcnp_str_putc(&s, ')');
                }
            }
        } else {
            mcnp_str_puts(&s, " FILL=");
            mcnp_str_int(&s, cell->lat_fill_dims[0]);
            mcnp_str_putc(&s, ':');
            mcnp_str_int(&s, cell->lat_fill_dims[1]);

            if (nj > 1 || nk > 1) {
                mcnp_str_putc(&s, ' ');
                mcnp_str_int(&s, cell->lat_fill_dims[2]);
                mcnp_str_putc(&s, ':');
                mcnp_str_int(&s, cell->lat_fill_dims[3]);
            }

            if (nk > 1) {
                mcnp_str_putc(&s, ' ');
                mcnp_str_int(&s, cell->lat_fill_dims[4]);
                mcnp_str_putc(&s, ':');
                mcnp_str_int(&s, cell->lat_fill_dims[5]);
            }

            for (size_t i = 0; i < cell->lat_fill_count; i++) {
                mcnp_str_putc(&s, ' ');
                mcnp_str_int(&s, cell->lat_fill[i]);
            }
        }
    }

    /* TRCL= or *TRCL= - cell transformation */
    if (mp && mp->has_trcl && mp->trcl != 0
        && ctx->trcl_export_mode != ALEA_TRCL_EXPORT_BAKE) {
        const alea_transform_t* tr = alea_get_transform(sys, mp->trcl);
        int use_degrees = tr ? tr->degrees : mp->trcl_degrees;

        if (use_degrees) {
            mcnp_str_puts(&s, " *TRCL=");
        } else {
            mcnp_str_puts(&s, " TRCL=");
        }

        bool trcl_inline = tr && tr->value_count > 0 &&
                           ((ctx->transform_export_mode == ALEA_TR_EXPORT_INLINE) ||
                            (ctx->transform_export_mode == ALEA_TR_EXPORT_ORIGINAL && tr->from_inline));
        if (trcl_inline) {
            mcnp_str_putc(&s, '(');
            for (int i = 0; i < tr->value_count; i++) {
                if (i > 0) mcnp_str_putc(&s, ' ');
                /* tr->data contains original values (degrees or cosines) */
                mcnp_str_double(&s, tr->data[i], 10);
            }
            mcnp_str_putc(&s, ')');
        } else {
            mcnp_str_int(&s, mp->trcl);
        }
    }

    if (inline_comment && *inline_comment) {
        mcnp_str_inline_comment(&s, inline_comment);
    }

    mcnp_str_write(&s, out);
}

/* ============================================================================
 * MCNP SURFACE WRITER
 * ============================================================================ */

static int write_mcnp_surface(FILE* out,
                              arena_t* arena,
                              int surface_id,
                              int transform_id,
                              alea_boundary_type_t boundary_type,
                              int periodic_surface_id,
                              alea_primitive_type_t type,
                              const alea_primitive_data_t* data,
                              int8_t inverted,
                              int max_col,
                              int cont_indent) {

    mcnp_str_t s;
    mcnp_str_init(&s, arena, SURFACE_BUF_SIZE, max_col, cont_indent);

    /* Boundary type prefix */
    switch (boundary_type) {
        case ALEA_BOUNDARY_REFLECTIVE:
            mcnp_str_putc(&s, '*');
            break;
        case ALEA_BOUNDARY_WHITE:
            mcnp_str_putc(&s, '+');
            break;
        default:
            break;
    }

    /* Surface ID */
    mcnp_str_int(&s, surface_id);
    mcnp_str_putc(&s, ' ');

    /* Transform ID or periodic surface reference */
    if (boundary_type == ALEA_BOUNDARY_PERIODIC && periodic_surface_id != 0) {
        mcnp_str_int(&s, -periodic_surface_id);
        mcnp_str_putc(&s, ' ');
    } else if (transform_id != 0) {
        mcnp_str_int(&s, transform_id);
        mcnp_str_putc(&s, ' ');
    }

    switch (type) {
        case ALEA_PRIMITIVE_PLANE: {
            double a = data->plane.a;
            double b = data->plane.b;
            double c = data->plane.c;
            double d = data->plane.d;
            if (inverted) {
                a = -a; b = -b; c = -c; d = -d;
            }
            if (fabs(b) < 1e-10 && fabs(c) < 1e-10 && fabs(a - 1.0) < 1e-10) {
                mcnp_str_puts(&s, "PX ");
                mcnp_str_double(&s, -d, 16);
            } else if (fabs(a) < 1e-10 && fabs(c) < 1e-10 && fabs(b - 1.0) < 1e-10) {
                mcnp_str_puts(&s, "PY ");
                mcnp_str_double(&s, -d, 16);
            } else if (fabs(a) < 1e-10 && fabs(b) < 1e-10 && fabs(c - 1.0) < 1e-10) {
                mcnp_str_puts(&s, "PZ ");
                mcnp_str_double(&s, -d, 16);
            } else {
                mcnp_str_puts(&s, "P ");
                mcnp_str_double(&s, a, 16);
                mcnp_str_putc(&s, ' ');
                mcnp_str_double(&s, b, 16);
                mcnp_str_putc(&s, ' ');
                mcnp_str_double(&s, c, 16);
                mcnp_str_putc(&s, ' ');
                mcnp_str_double(&s, -d, 16);
            }
            break;
        }

        case ALEA_PRIMITIVE_SPHERE: {
            const alea_sphere_data_t* sp = &data->sphere;
            if (fabs(sp->center_x) < 1e-10 && fabs(sp->center_y) < 1e-10 && fabs(sp->center_z) < 1e-10) {
                mcnp_str_puts(&s, "SO ");
                mcnp_str_double(&s, sp->radius, 16);
            } else {
                mcnp_str_puts(&s, "S ");
                mcnp_str_double(&s, sp->center_x, 16);
                mcnp_str_putc(&s, ' ');
                mcnp_str_double(&s, sp->center_y, 16);
                mcnp_str_putc(&s, ' ');
                mcnp_str_double(&s, sp->center_z, 16);
                mcnp_str_putc(&s, ' ');
                mcnp_str_double(&s, sp->radius, 16);
            }
            break;
        }

        case ALEA_PRIMITIVE_CYLINDER_X: {
            const alea_cylinder_x_data_t* c = &data->cyl_x;
            if (fabs(c->center_y) < 1e-10 && fabs(c->center_z) < 1e-10) {
                mcnp_str_puts(&s, "CX ");
                mcnp_str_double(&s, c->radius, 16);
            } else {
                mcnp_str_puts(&s, "C/X ");
                mcnp_str_double(&s, c->center_y, 16);
                mcnp_str_putc(&s, ' ');
                mcnp_str_double(&s, c->center_z, 16);
                mcnp_str_putc(&s, ' ');
                mcnp_str_double(&s, c->radius, 16);
            }
            break;
        }

        case ALEA_PRIMITIVE_CYLINDER_Y: {
            const alea_cylinder_y_data_t* c = &data->cyl_y;
            if (fabs(c->center_x) < 1e-10 && fabs(c->center_z) < 1e-10) {
                mcnp_str_puts(&s, "CY ");
                mcnp_str_double(&s, c->radius, 16);
            } else {
                mcnp_str_puts(&s, "C/Y ");
                mcnp_str_double(&s, c->center_x, 16);
                mcnp_str_putc(&s, ' ');
                mcnp_str_double(&s, c->center_z, 16);
                mcnp_str_putc(&s, ' ');
                mcnp_str_double(&s, c->radius, 16);
            }
            break;
        }

        case ALEA_PRIMITIVE_CYLINDER_Z: {
            const alea_cylinder_z_data_t* c = &data->cyl_z;
            if (fabs(c->center_x) < 1e-10 && fabs(c->center_y) < 1e-10) {
                mcnp_str_puts(&s, "CZ ");
                mcnp_str_double(&s, c->radius, 16);
            } else {
                mcnp_str_puts(&s, "C/Z ");
                mcnp_str_double(&s, c->center_x, 16);
                mcnp_str_putc(&s, ' ');
                mcnp_str_double(&s, c->center_y, 16);
                mcnp_str_putc(&s, ' ');
                mcnp_str_double(&s, c->radius, 16);
            }
            break;
        }

        case ALEA_PRIMITIVE_CONE_X: {
            const alea_cone_x_data_t* c = &data->cone_x;
            const int on_x_axis = (fabs(c->apex_y) < 1e-12 && fabs(c->apex_z) < 1e-12);
            if (on_x_axis) {
                mcnp_str_puts(&s, "KX ");
                mcnp_str_double(&s, c->apex_x, 16);
                mcnp_str_putc(&s, ' ');
            } else {
                mcnp_str_puts(&s, "K/X ");
                mcnp_str_double(&s, c->apex_x, 16);
                mcnp_str_putc(&s, ' ');
                mcnp_str_double(&s, c->apex_y, 16);
                mcnp_str_putc(&s, ' ');
                mcnp_str_double(&s, c->apex_z, 16);
                mcnp_str_putc(&s, ' ');
            }
            mcnp_str_double(&s, c->tan_angle_sq, 16);
            if (c->sheet_selection != 0) {
                mcnp_str_putc(&s, ' ');
                mcnp_str_int(&s, c->sheet_selection);
            }
            break;
        }

        case ALEA_PRIMITIVE_CONE_Y: {
            const alea_cone_y_data_t* c = &data->cone_y;
            const int on_y_axis = (fabs(c->apex_x) < 1e-12 && fabs(c->apex_z) < 1e-12);
            if (on_y_axis) {
                mcnp_str_puts(&s, "KY ");
                mcnp_str_double(&s, c->apex_y, 16);
                mcnp_str_putc(&s, ' ');
            } else {
                mcnp_str_puts(&s, "K/Y ");
                mcnp_str_double(&s, c->apex_x, 16);
                mcnp_str_putc(&s, ' ');
                mcnp_str_double(&s, c->apex_y, 16);
                mcnp_str_putc(&s, ' ');
                mcnp_str_double(&s, c->apex_z, 16);
                mcnp_str_putc(&s, ' ');
            }
            mcnp_str_double(&s, c->tan_angle_sq, 16);
            if (c->sheet_selection != 0) {
                mcnp_str_putc(&s, ' ');
                mcnp_str_int(&s, c->sheet_selection);
            }
            break;
        }

        case ALEA_PRIMITIVE_CONE_Z: {
            const alea_cone_z_data_t* c = &data->cone_z;
            const int on_z_axis = (fabs(c->apex_x) < 1e-12 && fabs(c->apex_y) < 1e-12);
            if (on_z_axis) {
                mcnp_str_puts(&s, "KZ ");
                mcnp_str_double(&s, c->apex_z, 16);
                mcnp_str_putc(&s, ' ');
            } else {
                mcnp_str_puts(&s, "K/Z ");
                mcnp_str_double(&s, c->apex_x, 16);
                mcnp_str_putc(&s, ' ');
                mcnp_str_double(&s, c->apex_y, 16);
                mcnp_str_putc(&s, ' ');
                mcnp_str_double(&s, c->apex_z, 16);
                mcnp_str_putc(&s, ' ');
            }
            mcnp_str_double(&s, c->tan_angle_sq, 16);
            if (c->sheet_selection != 0) {
                mcnp_str_putc(&s, ' ');
                mcnp_str_int(&s, c->sheet_selection);
            }
            break;
        }

        case ALEA_PRIMITIVE_RPP: {
            const alea_box_data_t* b = &data->box;
            mcnp_str_puts(&s, "RPP ");
            mcnp_str_double(&s, b->min_x, 16);
            mcnp_str_putc(&s, ' ');
            mcnp_str_double(&s, b->max_x, 16);
            mcnp_str_putc(&s, ' ');
            mcnp_str_double(&s, b->min_y, 16);
            mcnp_str_putc(&s, ' ');
            mcnp_str_double(&s, b->max_y, 16);
            mcnp_str_putc(&s, ' ');
            mcnp_str_double(&s, b->min_z, 16);
            mcnp_str_putc(&s, ' ');
            mcnp_str_double(&s, b->max_z, 16);
            break;
        }

        case ALEA_PRIMITIVE_RCC: {
            const alea_rcc_data_t* r = &data->rcc;
            mcnp_str_puts(&s, "RCC ");
            mcnp_str_double(&s, r->base_x, 16);
            mcnp_str_putc(&s, ' ');
            mcnp_str_double(&s, r->base_y, 16);
            mcnp_str_putc(&s, ' ');
            mcnp_str_double(&s, r->base_z, 16);
            mcnp_str_putc(&s, ' ');
            mcnp_str_double(&s, r->height_x, 16);
            mcnp_str_putc(&s, ' ');
            mcnp_str_double(&s, r->height_y, 16);
            mcnp_str_putc(&s, ' ');
            mcnp_str_double(&s, r->height_z, 16);
            mcnp_str_putc(&s, ' ');
            mcnp_str_double(&s, r->radius, 16);
            break;
        }

        case ALEA_PRIMITIVE_TORUS_X:
        case ALEA_PRIMITIVE_TORUS_Y:
        case ALEA_PRIMITIVE_TORUS_Z: {
            const alea_torus_data_t* t = &data->torus;
            const char* mnemonic = (type == ALEA_PRIMITIVE_TORUS_X) ? "TX " :
                                   (type == ALEA_PRIMITIVE_TORUS_Y) ? "TY " : "TZ ";
            mcnp_str_puts(&s, mnemonic);
            mcnp_str_double(&s, t->center_x, 16);
            mcnp_str_putc(&s, ' ');
            mcnp_str_double(&s, t->center_y, 16);
            mcnp_str_putc(&s, ' ');
            mcnp_str_double(&s, t->center_z, 16);
            mcnp_str_putc(&s, ' ');
            mcnp_str_double(&s, t->major_radius, 16);
            mcnp_str_putc(&s, ' ');
            mcnp_str_double(&s, t->axial_semiwidth_B, 16);
            mcnp_str_putc(&s, ' ');
            mcnp_str_double(&s, t->minor_radius, 16);
            break;
        }

        case ALEA_PRIMITIVE_QUADRIC: {
            const alea_quadric_data_t* q = &data->quadric;
            mcnp_str_puts(&s, "GQ ");
            for (int i = 0; i < 10; i++) {
                mcnp_str_double(&s, q->coeffs[i], 16);
                if (i < 9) mcnp_str_putc(&s, ' ');
            }
            break;
        }

        /* Macrobody surface cards. MCNP supports these natively at the
         * surface-card level; we just need to write the right mnemonic and
         * the data fields in the order MCNP expects. */

        case ALEA_PRIMITIVE_BOX: {
            const alea_box_general_data_t* b = &data->box_general;
            const double v[12] = {
                b->corner_x, b->corner_y, b->corner_z,
                b->v1_x, b->v1_y, b->v1_z,
                b->v2_x, b->v2_y, b->v2_z,
                b->v3_x, b->v3_y, b->v3_z,
            };
            mcnp_str_puts(&s, "BOX ");
            for (int i = 0; i < 12; i++) {
                mcnp_str_double(&s, v[i], 16);
                if (i < 11) mcnp_str_putc(&s, ' ');
            }
            break;
        }

        case ALEA_PRIMITIVE_SPH: {
            const alea_sph_data_t* sp = &data->sph;
            mcnp_str_puts(&s, "SPH ");
            mcnp_str_double(&s, sp->center_x, 16); mcnp_str_putc(&s, ' ');
            mcnp_str_double(&s, sp->center_y, 16); mcnp_str_putc(&s, ' ');
            mcnp_str_double(&s, sp->center_z, 16); mcnp_str_putc(&s, ' ');
            mcnp_str_double(&s, sp->radius, 16);
            break;
        }

        case ALEA_PRIMITIVE_TRC: {
            const alea_trc_data_t* t = &data->trc;
            const double v[8] = {
                t->base_x, t->base_y, t->base_z,
                t->height_x, t->height_y, t->height_z,
                t->base_radius, t->top_radius,
            };
            mcnp_str_puts(&s, "TRC ");
            for (int i = 0; i < 8; i++) {
                mcnp_str_double(&s, v[i], 16);
                if (i < 7) mcnp_str_putc(&s, ' ');
            }
            break;
        }

        case ALEA_PRIMITIVE_ELL: {
            /* ELL has two parametrizations:
             *   - Foci form: V1 (focus 1), V2 (focus 2), RM > 0 (major axis length)
             *   - Semi-axes form: V1 (center), V2 (semi-axis lengths), RM < 0
             * libalea stores both fields raw; the sign on major_axis_len
             * selects the form. We preserve that sign on export. */
            const alea_ell_data_t* e = &data->ell;
            mcnp_str_puts(&s, "ELL ");
            mcnp_str_double(&s, e->v1_x, 16); mcnp_str_putc(&s, ' ');
            mcnp_str_double(&s, e->v1_y, 16); mcnp_str_putc(&s, ' ');
            mcnp_str_double(&s, e->v1_z, 16); mcnp_str_putc(&s, ' ');
            mcnp_str_double(&s, e->v2_x, 16); mcnp_str_putc(&s, ' ');
            mcnp_str_double(&s, e->v2_y, 16); mcnp_str_putc(&s, ' ');
            mcnp_str_double(&s, e->v2_z, 16); mcnp_str_putc(&s, ' ');
            mcnp_str_double(&s, e->major_axis_len, 16);
            break;
        }

        case ALEA_PRIMITIVE_REC: {
            /* REC: base center + height vector + 2 perpendicular semi-axis
             * vectors. MCNP also accepts a short form (one axis + scalar
             * radius); we always emit the full 12-value form. */
            const alea_rec_data_t* r = &data->rec;
            const double v[12] = {
                r->base_x, r->base_y, r->base_z,
                r->height_x, r->height_y, r->height_z,
                r->axis1_x, r->axis1_y, r->axis1_z,
                r->axis2_x, r->axis2_y, r->axis2_z,
            };
            mcnp_str_puts(&s, "REC ");
            for (int i = 0; i < 12; i++) {
                mcnp_str_double(&s, v[i], 16);
                if (i < 11) mcnp_str_putc(&s, ' ');
            }
            break;
        }

        case ALEA_PRIMITIVE_WED: {
            const alea_wed_data_t* w = &data->wed;
            const double v[12] = {
                w->vertex_x, w->vertex_y, w->vertex_z,
                w->v1_x, w->v1_y, w->v1_z,
                w->v2_x, w->v2_y, w->v2_z,
                w->v3_x, w->v3_y, w->v3_z,
            };
            mcnp_str_puts(&s, "WED ");
            for (int i = 0; i < 12; i++) {
                mcnp_str_double(&s, v[i], 16);
                if (i < 11) mcnp_str_putc(&s, ' ');
            }
            break;
        }

        case ALEA_PRIMITIVE_RHP: {
            /* RHP: base center + height vector + three semi-vertex vectors
             * (r1, r2, r3). MCNP also accepts a short form with only r1
             * when the hex is regular; we always emit the full 15-value
             * form (matches HEX alias which is the same primitive). */
            const alea_rhp_data_t* h = &data->rhp;
            const double v[15] = {
                h->base_x, h->base_y, h->base_z,
                h->height_x, h->height_y, h->height_z,
                h->r1_x, h->r1_y, h->r1_z,
                h->r2_x, h->r2_y, h->r2_z,
                h->r3_x, h->r3_y, h->r3_z,
            };
            mcnp_str_puts(&s, "RHP ");
            for (int i = 0; i < 15; i++) {
                mcnp_str_double(&s, v[i], 16);
                if (i < 14) mcnp_str_putc(&s, ' ');
            }
            break;
        }

        case ALEA_PRIMITIVE_ARB: {
            /* ARB: 8 corner triples (24 values) + 6 face descriptors
             * (each is a 4-digit integer whose digits are 1..8 corner
             * indices, or 0-padded for triangles). MCNP requires exactly
             * 30 values regardless of num_corners / num_faces. */
            const alea_arb_data_t* a = &data->arb;
            mcnp_str_puts(&s, "ARB ");
            for (int i = 0; i < 8; i++) {
                for (int j = 0; j < 3; j++) {
                    mcnp_str_double(&s, a->corners[i][j], 16);
                    mcnp_str_putc(&s, ' ');
                }
            }
            for (int f = 0; f < 6; f++) {
                /* Encode face as concatenated digits: e.g. {1,2,3,4} -> 1234,
                 * {1,2,3,0} (triangle) -> 1230. */
                int code = 0;
                for (int j = 0; j < 4; j++) {
                    int d = a->faces[f][j];
                    if (d < 0 || d > 9) d = 0;
                    code = code * 10 + d;
                }
                mcnp_str_int(&s, code);
                if (f < 5) mcnp_str_putc(&s, ' ');
            }
            break;
        }

        default:
            ALEA_LOG_WARN("Unknown primitive type %d", type);
            return -1;
    }

    mcnp_str_write(&s, out);
    return 0;
}

/* ============================================================================
 * MCNP EXPORT - MAIN FUNCTION
 * ============================================================================ */

int export_mcnp(alea_system_t* sys, export_context_t* ctx) {
    /* Reusable builder for section comments */
    mcnp_str_t cs;
    mcnp_str_init(&cs, &ctx->arena, 256, ctx->mcnp_max_col, ctx->mcnp_cont_indent);

    /* Header */
    fprintf(ctx->out, "MCNP Input - Exported by CSG Library\n");
    mcnp_str_comment(&cs, "");
    mcnp_str_comment(&cs, "Cell Cards");
    mcnp_str_comment(&cs, "");
    mcnp_str_write(&cs, ctx->out);
    mcnp_str_reset(&cs);

    /* Assign surface IDs to primitives from expanded macrobodies */
    alea_assign_missing_surface_ids(sys, ctx);

    /* Build canonical surface map if deduplicating */
    if (ctx->deduplicate) {
        alea_build_canonical_surface_map(ctx, sys);
    }

    /* Macrobody expansion when ALEA_EMIT_SURFACES policy is set */
    if (ctx->surface_policy == ALEA_EMIT_SURFACES) {
        int expanded = alea_expand_macrobodies_tree_level(sys, ctx);
        if (expanded < 0) {
            ALEA_LOG_ERROR("Tree-level macrobody expansion failed");
            return -1;
        }
    }

    /* Build flatten info: collapse 1x1x1 LAT cells whose content universe is a
     * single passthrough cell (no geometry, just FILL=X with optional transform).
     * The LAT cell takes the passthrough's target universe and transform; the
     * passthrough cell itself is skipped on emission. This matches the cleaner
     * idiomatic MCNP form for OpenMC-style positioning-wrapper lattices.
     *
     * When the passthrough's translation T equals -(lower_left + pitch/2) — the
     * canonical OpenMC normalization pattern — we additionally drop the FILL
     * transform and shift the LAT cell's bounding RPP by (lower_left + pitch/2),
     * so the cell occupies its actual region in the parent universe frame. */
    size_t n_cells = alea_vec_count(&sys->cells);
    flatten_cell_info_t* flatten = NULL;
    surface_offset_t* surf_offsets = NULL;
    size_t surf_offset_count = 0;
    if (n_cells > 0) {
        flatten = calloc(n_cells, sizeof(*flatten));
        surf_offsets = calloc(n_cells, sizeof(*surf_offsets));
    }
    if (flatten && surf_offsets) {
        for (size_t i = 0; i < n_cells; i++) {
            const alea_cell_entry_t* lc = &sys->cells.data[i];
            if (lc->lat_type == 0 || !lc->lat_fill || lc->lat_fill_count != 1) continue;
            int ni = lc->lat_fill_dims[1] - lc->lat_fill_dims[0] + 1;
            int nj = lc->lat_fill_dims[3] - lc->lat_fill_dims[2] + 1;
            int nk = lc->lat_fill_dims[5] - lc->lat_fill_dims[4] + 1;
            if (ni != 1 || nj != 1 || nk != 1) continue;

            int u_int = lc->lat_fill[0];
            const alea_universe_t* u = alea_get_universe(sys, u_int);
            if (!u || u->cell_indices.count != 1) continue;

            size_t inner_idx = u->cell_indices.data[0];
            if (inner_idx >= n_cells) continue;
            const alea_cell_entry_t* inner = &sys->cells.data[inner_idx];
            if (inner->fill_universe <= 0) continue;
            if (inner->root_node_id != ALEA_NODE_ID_INVALID) continue;
            if (inner->lat_type != 0) continue;

            /* Element center C = lower_left + pitch/2 in the lattice frame. */
            double cx = lc->lat_lower_left[0] + 0.5 * lc->lat_pitch[0];
            double cy = lc->lat_lower_left[1] + 0.5 * lc->lat_pitch[1];
            double cz = lc->lat_lower_left[2] + 0.5 * lc->lat_pitch[2];

            /* Check if the passthrough's translation cancels C: T == -C. */
            int cancels = 0;
            if (inner->fill_transform != 0) {
                const alea_transform_t* tr = alea_get_transform(sys, inner->fill_transform);
                if (tr) {
                    const double eps = 1e-6;
                    if (fabs(tr->data[0] + cx) < eps &&
                        fabs(tr->data[1] + cy) < eps &&
                        fabs(tr->data[2] + cz) < eps) {
                        cancels = 1;
                    }
                }
            }

            flatten[i].override_fill_universe = inner->fill_universe;
            flatten[inner_idx].skip = true;

            if (cancels) {
                /* Drop the FILL transform and shift the bounding RPP by +C so the
                 * cell occupies [lower_left, lower_left + pitch] in the parent. */
                flatten[i].override_fill_transform = 0;
                const alea_node_t* root = (lc->root_node_id < alea_vec_count(&sys->nodes))
                    ? &sys->nodes.data[lc->root_node_id] : NULL;
                if (root && ALEA_GET_OPERATION(root) == ALEA_OP_PRIMITIVE &&
                    root->primitive.mc_surface_id != 0) {
                    surf_offsets[surf_offset_count].surface_id = root->primitive.mc_surface_id;
                    surf_offsets[surf_offset_count].off_x = cx;
                    surf_offsets[surf_offset_count].off_y = cy;
                    surf_offsets[surf_offset_count].off_z = cz;
                    surf_offset_count++;
                }
            } else {
                flatten[i].override_fill_transform = inner->fill_transform;
            }
        }
    }

    /* Prepare universe depth cache for filtering */
    int* depth_cache = NULL;
    size_t depth_cache_size = 0;
    if (ctx->universe_depth >= 0) {
        int max_univ = 0;
        for (size_t i = 0; i < alea_vec_count(&sys->cells); i++) {
            if (sys->cells.data[i].universe_id > max_univ) max_univ = sys->cells.data[i].universe_id;
            if (sys->cells.data[i].fill_universe > max_univ) max_univ = sys->cells.data[i].fill_universe;
        }
        depth_cache_size = (size_t)(max_univ + 1);
        depth_cache = malloc(depth_cache_size * sizeof(int));
        if (depth_cache) {
            for (size_t i = 0; i < depth_cache_size; i++) {
                depth_cache[i] = -2;
            }
        }
    }

    /* Write cells */
    for (size_t i = 0; i < alea_vec_count(&sys->cells); i++) {
        ALEA_CHECK_INTERRUPTED(-1);
        const alea_cell_entry_t* cell = &sys->cells.data[i];

        /* Skip passthrough cells folded into their parent LAT cell */
        if (flatten && flatten[i].skip) continue;

        /* Apply universe depth filter */
        if (ctx->universe_depth >= 0 &&
            !alea_should_export_cell(sys, cell, ctx->universe_depth, depth_cache, depth_cache_size)) {
            continue;
        }

        /* Use original (pre-TRCL) tree when preserving TRCL */
        uint32_t export_root = cell->root_node_id;
        const mcnp_model_t* exp_model = (const mcnp_model_t*)ctx->module_data;
        const mcnp_cell_params_t* exp_mp = exp_model ? mcnp_cell_params_const(exp_model, i) : NULL;
        if (exp_mp && exp_mp->has_trcl && cell->original_root_node_id != ALEA_NODE_ID_INVALID
            && ctx->trcl_export_mode == ALEA_TRCL_EXPORT_PRESERVE) {
            export_root = cell->original_root_node_id;
        }

        mcnp_str_t s;
        mcnp_str_init_unwrapped(&s, &ctx->arena, EXPR_BUF_SIZE);
        int len = alea_tree_to_mcnp_expr(sys, ctx, export_root, &s);
        if (len < 0) {
            ALEA_LOG_ERROR("Error converting cell %d to expression", cell->mc_cell_id);
            continue;
        }

        // Emit "C" comment lines before the cell card
        if (cell->comments && *cell->comments) {
            mcnp_str_t cs_cell;
            mcnp_str_init(&cs_cell, &ctx->arena, 1024, ctx->mcnp_max_col, ctx->mcnp_cont_indent);
            const char* cp = cell->comments;
            while (*cp) {
                const char* nl = strchr(cp, '\n');
                size_t line_len = nl ? (size_t)(nl - cp) : strlen(cp);
                // Add "c " prefix if line doesn't already have it
                bool has_prefix = (line_len >= 1 && (cp[0] == 'c' || cp[0] == 'C') &&
                                   (line_len == 1 || cp[1] == ' ' || cp[1] == '\n'));
                if (!has_prefix && line_len > 0) {
                    str_builder_write(&cs_cell.sb, "c ", 2);
                }
                str_builder_write(&cs_cell.sb, cp, line_len);
                str_builder_putc(&cs_cell.sb, '\n');
                cs_cell.col = 1;
                cp = nl ? nl + 1 : cp + line_len;
            }
            str_builder_finish(&cs_cell.sb);
            if (cs_cell.sb.len > 0) {
                fprintf(ctx->out, "%s", cs_cell.sb.buf);
            }
        }

        write_mcnp_cell_line(ctx->out, cell, i, s.sb.buf, sys, ctx,
                             flatten ? &flatten[i] : NULL, cell->inline_comment);
        ctx->cells_written++;
    }

    free(depth_cache);

    /* Blank line separator */
    fprintf(ctx->out, "\n");
    mcnp_str_comment(&cs, "Surface Cards");
    mcnp_str_comment(&cs, "");
    mcnp_str_write(&cs, ctx->out);
    mcnp_str_reset(&cs);

    /* Write surfaces */
    if (ctx->deduplicate && ctx->prim_to_surface) {
        alea_bitset_t prim_written = alea_bitset_create(ctx->prim_to_surface_size);

        for (size_t i = 0; i < alea_vec_count(&sys->surfaces); i++) {
            ALEA_CHECK_INTERRUPTED(-1);
            const alea_surface_entry_t* surface = &sys->surfaces.data[i];

            int mcnp_id = surface->mc_surface_id;
            uint32_t prim_id = surface->primitive_id;

            /* Only emit the canonical surface for each primitive (the one whose
               mc_surface_id matches prim_to_surface[prim_id]).  Cell expressions
               reference prim_to_surface[prim_id], so we must emit that ID. */
            if (prim_id < ctx->prim_to_surface_size &&
                ctx->prim_to_surface[prim_id] >= 0 &&
                mcnp_id != ctx->prim_to_surface[prim_id]) {
                continue;
            }

            if (alea_bitset_test(&prim_written, prim_id)) continue;

            if (prim_id >= alea_vec_count(&sys->primitives)) {
                ALEA_LOG_ERROR("Surface entry %zu (mc_surface_id=%d) has out-of-range "
                               "primitive_id=%u (primitive_count=%zu) in dedup path",
                               i, mcnp_id, prim_id, alea_vec_count(&sys->primitives));
                alea_bitset_destroy(&prim_written);
                return -1;
            }
            const alea_primitive_entry_t* prim = &sys->primitives.data[prim_id];

            if (ctx->surface_policy == ALEA_EMIT_SURFACES && is_macrobody(prim->type)) {
                alea_bitset_set(&prim_written, prim_id);
                continue;
            }

            int8_t inverted = 0;
            if (surface->pos_node < alea_vec_count(&sys->nodes)) {
                inverted = sys->nodes.data[surface->pos_node].primitive.inverted;
            }

            alea_primitive_type_t export_type = prim->type;
            alea_primitive_data_t source_data;
            alea_primitive_data_t export_data;
            if (!alea_primitive_copy_data(sys, prim_id, &source_data)) {
                alea_bitset_destroy(&prim_written);
                return -1;
            }
            export_data = source_data;
            int export_transform_id = 0;

            if (surface->transform_applied && surface->transform_id != 0) {
                const alea_transform_t* tr = alea_get_transform(sys, surface->transform_id);
                if (tr && alea_apply_inverse_transform_to_primitive(tr, prim->type, &source_data,
                                                                    &export_type, &export_data)) {
                    export_transform_id = surface->transform_id;
                }
            }

            apply_flatten_surface_offset(mcnp_id, surf_offsets, surf_offset_count,
                                         export_type, &export_data);

            write_mcnp_surface(ctx->out, &ctx->arena, mcnp_id, export_transform_id, surface->boundary_type,
                surface->periodic_surface_id, export_type, &export_data, inverted,
                ctx->mcnp_max_col, ctx->mcnp_cont_indent);
            ctx->surfaces_written++;
            alea_bitset_set(&prim_written, prim_id);
        }

        alea_bitset_destroy(&prim_written);
    } else {
        for (size_t i = 0; i < alea_vec_count(&sys->surfaces); i++) {
            ALEA_CHECK_INTERRUPTED(-1);
            const alea_surface_entry_t* surface = &sys->surfaces.data[i];

            /* Bounds-check primitive_id before dereferencing. If a surface
             * entry references an out-of-range primitive, reading past the
             * primitives vector yields garbage — most visibly a "very high"
             * prim->type printed by the "Unknown primitive type" branch
             * below. Fail loudly with diagnostics instead. */
            if (surface->primitive_id >= alea_vec_count(&sys->primitives)) {
                ALEA_LOG_ERROR("Surface entry %zu (mc_surface_id=%d) has out-of-range "
                               "primitive_id=%u (primitive_count=%zu) — corrupted surface table",
                               i, surface->mc_surface_id, surface->primitive_id,
                               alea_vec_count(&sys->primitives));
                return -1;
            }
            const alea_primitive_entry_t* prim = &sys->primitives.data[surface->primitive_id];
            if (ctx->surface_policy == ALEA_EMIT_SURFACES && is_macrobody(prim->type)) {
                continue;
            }

            int8_t inverted = 0;
            if (surface->pos_node < alea_vec_count(&sys->nodes)) {
                inverted = sys->nodes.data[surface->pos_node].primitive.inverted;
            }

            alea_primitive_type_t export_type = prim->type;
            alea_primitive_data_t source_data;
            alea_primitive_data_t export_data;
            if (!alea_primitive_copy_data(sys, surface->primitive_id, &source_data)) {
                return -1;
            }
            export_data = source_data;
            int export_transform_id = 0;

            if (surface->transform_applied && surface->transform_id != 0) {
                const alea_transform_t* tr = alea_get_transform(sys, surface->transform_id);
                if (tr && alea_apply_inverse_transform_to_primitive(tr, prim->type, &source_data,
                                                                    &export_type, &export_data)) {
                    export_transform_id = surface->transform_id;
                }
            }

            apply_flatten_surface_offset(surface->mc_surface_id, surf_offsets, surf_offset_count,
                                         export_type, &export_data);

            write_mcnp_surface(ctx->out, &ctx->arena, surface->mc_surface_id, export_transform_id,
                surface->boundary_type, surface->periodic_surface_id, export_type, &export_data, inverted,
                ctx->mcnp_max_col, ctx->mcnp_cont_indent);
            ctx->surfaces_written++;
        }
    }

    /* Write decomposed surfaces (legacy code path) */
    if (ctx->decompositions && ctx->decomposition_count > 0) {
        mcnp_str_comment(&cs, "");
        mcnp_str_comment(&cs, "Decomposed macrobody surfaces");
        mcnp_str_comment(&cs, "");
        mcnp_str_write(&cs, ctx->out);
        mcnp_str_reset(&cs);

        for (size_t d = 0; d < ctx->decomposition_count; d++) {
            const macrobody_decomposition_t* decomp = &ctx->decompositions[d];

            for (size_t si = 0; si < decomp->surface_count; si++) {
                int surf_id = decomp->surfaces[si].surface_id;

                for (size_t n = 0; n < alea_vec_count(&sys->nodes); n++) {
                    const alea_node_t* node = &sys->nodes.data[n];
                    if (ALEA_GET_OPERATION(node) == ALEA_OP_PRIMITIVE &&
                        node->primitive.mc_surface_id == surf_id) {
                        const alea_primitive_entry_t* prim = &sys->primitives.data[node->primitive.primitive_id];
                        alea_primitive_data_t prim_data;
                        if (!alea_primitive_copy_data(sys, node->primitive.primitive_id, &prim_data)) {
                            continue;
                        }
                        write_mcnp_surface(ctx->out, &ctx->arena, surf_id, 0, ALEA_BOUNDARY_TRANSMISSIVE, 0,
                                           prim->type, &prim_data, node->primitive.inverted,
                                           ctx->mcnp_max_col, ctx->mcnp_cont_indent);
                        ctx->surfaces_written++;
                        break;
                    }
                }
            }
        }
    }

    /* Write synthetic surfaces from immediate macrobody expansion */
    if (ctx->synthetic_surfaces_created > 0) {
        mcnp_str_comment(&cs, "");
        mcnp_str_comment(&cs, "Expanded macrobody component surfaces");
        mcnp_str_comment(&cs, "");
        mcnp_str_write(&cs, ctx->out);
        mcnp_str_reset(&cs);

        int max_registered_id = 0;
        for (size_t i = 0; i < alea_vec_count(&sys->surfaces); i++) {
            if (sys->surfaces.data[i].mc_surface_id > max_registered_id) {
                max_registered_id = sys->surfaces.data[i].mc_surface_id;
            }
        }

        bool* registered = NULL;
        if (max_registered_id > 0) {
            registered = calloc(max_registered_id + 1, sizeof(bool));
            if (registered) {
                for (size_t i = 0; i < alea_vec_count(&sys->surfaces); i++) {
                    registered[sys->surfaces.data[i].mc_surface_id] = true;
                }
            }
        }

        int* written_ids = calloc(ctx->next_synthetic_surface_id + 1, sizeof(int));
        if (written_ids) {
            for (size_t n = 0; n < alea_vec_count(&sys->nodes); n++) {
                const alea_node_t* node = &sys->nodes.data[n];
                if (ALEA_GET_OPERATION(node) != ALEA_OP_PRIMITIVE) continue;

                int surf_id = node->primitive.mc_surface_id;
                if (surf_id <= 0) continue;
                if (written_ids[surf_id]) continue;

                if (registered && surf_id <= max_registered_id && registered[surf_id]) continue;

                if (node->primitive.primitive_id >= alea_vec_count(&sys->primitives)) {
                    ALEA_LOG_ERROR("Node %zu (mc_surface_id=%d) has out-of-range "
                                   "primitive_id=%u (primitive_count=%zu) in synthetic-surface emission",
                                   n, surf_id, node->primitive.primitive_id,
                                   alea_vec_count(&sys->primitives));
                    continue;
                }
                const alea_primitive_entry_t* prim = &sys->primitives.data[node->primitive.primitive_id];
                alea_primitive_data_t prim_data;
                if (!alea_primitive_copy_data(sys, node->primitive.primitive_id, &prim_data)) {
                    continue;
                }
                write_mcnp_surface(ctx->out, &ctx->arena, surf_id, 0,
                                   ALEA_BOUNDARY_TRANSMISSIVE, 0, prim->type, &prim_data,
                                   node->primitive.inverted,
                                   ctx->mcnp_max_col, ctx->mcnp_cont_indent);
                ctx->surfaces_written++;
                written_ids[surf_id] = 1;
            }
            free(written_ids);
        }
        free(registered);
    }

    /* Materials */
    fprintf(ctx->out, "\n");
    mcnp_str_comment(&cs, "Data Cards");
    mcnp_str_comment(&cs, "");
    mcnp_str_comment(&cs, "Material Section");
    mcnp_str_comment(&cs, "");
    mcnp_str_write(&cs, ctx->out);
    mcnp_str_reset(&cs);

    for (size_t i = 0; i < alea_vec_count(&sys->materials); i++) {
        ALEA_CHECK_INTERRUPTED(-1);
        const alea_material_t* mat = &sys->materials.data[i];

        mcnp_str_t s;
        mcnp_str_init(&s, &ctx->arena, 1024, ctx->mcnp_max_col, ctx->mcnp_cont_indent);

        mcnp_str_putc(&s, 'M');
        mcnp_str_int(&s, mat->material_id);

        for (size_t j = 0; j < alea_vec_count(&mat->nuclides); j++) {
            const alea_nuclide_t* nuc = &mat->nuclides.data[j];

            mcnp_str_putc(&s, ' ');
            mcnp_str_int(&s, nuc->zaid);
            if (nuc->library) {
                mcnp_str_putc(&s, '.');
                mcnp_str_puts(&s, nuc->library);
            }
            mcnp_str_putc(&s, ' ');
            if (mat->is_weight_fraction) {
                mcnp_str_putc(&s, '-');
            }
            mcnp_str_double(&s, nuc->fraction, 6);

            int Z = alea_zaid_to_Z(nuc->zaid);
            int A = alea_zaid_to_A(nuc->zaid);
            const alea_element_t* elem = alea_get_element(Z);
            if (elem) {
                double elem_sum = 0.0;
                for (size_t k = 0; k < alea_vec_count(&mat->nuclides); k++) {
                    if (alea_zaid_to_Z(mat->nuclides.data[k].zaid) == Z) {
                        elem_sum += mat->nuclides.data[k].fraction;
                    }
                }

                double abundance = 0.0;
                for (size_t k = 0; k < elem->isotope_count; k++) {
                    if (elem->isotopes[k].mass_number == A) {
                        abundance = elem->isotopes[k].abundance * 100.0;
                        break;
                    }
                }
                char sym_upper[4];
                for (int k = 0; elem->symbol[k] && k < 3; k++) {
                    sym_upper[k] = (char)toupper((unsigned char)elem->symbol[k]);
                    sym_upper[k+1] = '\0';
                }
                double weight_pct = elem_sum > 0.0 ? (nuc->fraction / elem_sum) * 100.0 : 0.0;
                char comment[80];
                snprintf(comment, sizeof(comment), "%s %d WEIGHT(%%) %7.4f AB(%%) %6.2f",
                         sym_upper, A, weight_pct, abundance);
                mcnp_str_inline_comment(&s, comment);
            }

            if (j < alea_vec_count(&mat->nuclides)-1) mcnp_str_newline(&s);
        }

        mcnp_str_write(&s, ctx->out);
    }

    /* Mixtures — each mixture becomes a synthetic material card */
    for (size_t i = 0; i < alea_vec_count(&sys->mixtures); i++) {
        ALEA_CHECK_INTERRUPTED(-1);
        const alea_mixture_t* mix = &sys->mixtures.data[i];

        mcnp_str_t s;
        mcnp_str_init(&s, &ctx->arena, 1024, ctx->mcnp_max_col, ctx->mcnp_cont_indent);

        /* Use mixture_id as the material number */
        int mix_mat_id = mix->mixture_id;

        /* Header comment with mixture composition */
        {
            mcnp_str_t hdr;
            mcnp_str_init(&hdr, &ctx->arena, 256, ctx->mcnp_max_col, ctx->mcnp_cont_indent);
            if (mix->name && mix->name[0]) {
                mcnp_str_comment(&hdr, mix->name);
            }
            char comp_line[128];
            snprintf(comp_line, sizeof(comp_line), "Mixture M%d:", mix_mat_id);
            mcnp_str_comment(&hdr, comp_line);
            for (size_t c = 0; c < alea_vec_count(&mix->components); c++) {
                const alea_mixture_comp_t* comp = &mix->components.data[c];
                snprintf(comp_line, sizeof(comp_line), "  %.4g of M%d",
                         comp->fraction, comp->material_id);
                mcnp_str_comment(&hdr, comp_line);
            }
            mcnp_str_write(&hdr, ctx->out);
        }

        mcnp_str_putc(&s, 'M');
        mcnp_str_int(&s, mix_mat_id);

        /* Write nuclides from each component, scaled by component fraction */
        bool first_nuc = true;
        for (size_t c = 0; c < alea_vec_count(&mix->components); c++) {
            const alea_mixture_comp_t* comp = &mix->components.data[c];

            /* Find the component material */
            const alea_material_t* comp_mat = NULL;
            for (size_t m = 0; m < alea_vec_count(&sys->materials); m++) {
                if (sys->materials.data[m].material_id == comp->material_id) {
                    comp_mat = &sys->materials.data[m];
                    break;
                }
            }
            if (!comp_mat) {
                ALEA_LOG_WARN("Mixture M%d: component M%d not found, skipping",
                             mix_mat_id, comp->material_id);
                continue;
            }

            for (size_t j = 0; j < alea_vec_count(&comp_mat->nuclides); j++) {
                const alea_nuclide_t* nuc = &comp_mat->nuclides.data[j];

                if (!first_nuc) mcnp_str_newline(&s);
                first_nuc = false;

                mcnp_str_putc(&s, ' ');
                mcnp_str_int(&s, nuc->zaid);
                if (nuc->library) {
                    mcnp_str_putc(&s, '.');
                    mcnp_str_puts(&s, nuc->library);
                }
                mcnp_str_putc(&s, ' ');
                if (comp_mat->is_weight_fraction) {
                    mcnp_str_putc(&s, '-');
                }
                mcnp_str_double(&s, nuc->fraction * comp->fraction, 6);

                /* Inline comment: element info + source material */
                int Z = alea_zaid_to_Z(nuc->zaid);
                int A = alea_zaid_to_A(nuc->zaid);
                const alea_element_t* elem = alea_get_element(Z);
                if (elem) {
                    char sym_upper[4];
                    for (int k = 0; elem->symbol[k] && k < 3; k++) {
                        sym_upper[k] = (char)toupper((unsigned char)elem->symbol[k]);
                        sym_upper[k+1] = '\0';
                    }
                    char comment[80];
                    snprintf(comment, sizeof(comment), "%s %d from M%d",
                             sym_upper, A, comp->material_id);
                    mcnp_str_inline_comment(&s, comment);
                }
            }
        }

        mcnp_str_write(&s, ctx->out);
    }

    mcnp_str_comment(&cs, "END Material Section");
    mcnp_str_comment(&cs, "");
    mcnp_str_comment(&cs, "Transform Cards");
    mcnp_str_comment(&cs, "");
    mcnp_str_write(&cs, ctx->out);
    mcnp_str_reset(&cs);

    /* Transforms */
    for (size_t i = 0; i < alea_vec_count(&sys->transforms); i++) {
        const alea_transform_t* tr = &sys->transforms.data[i];

        /* ORIGINAL: skip inline-origin transforms (no TR card for them)
         * INLINE:    skip all transforms (everything emitted inline)
         * CARDS:     emit all transforms as TR cards */
        if (ctx->transform_export_mode == ALEA_TR_EXPORT_INLINE) {
            continue;
        }
        if (ctx->transform_export_mode == ALEA_TR_EXPORT_ORIGINAL && tr->from_inline) {
            continue;
        }

        mcnp_str_t s;
        mcnp_str_init(&s, &ctx->arena, 256, ctx->mcnp_max_col, ctx->mcnp_cont_indent);

        if (tr->degrees) {
            mcnp_str_puts(&s, "*TR");
        } else {
            mcnp_str_puts(&s, "TR");
        }
        mcnp_str_int(&s, tr->transform_id);
        mcnp_str_putc(&s, ' ');

        for (int j = 0; j < tr->value_count; j++) {
            if (j > 0) mcnp_str_putc(&s, ' ');
            /* tr->data contains original values (degrees or cosines) */
            mcnp_str_double(&s, tr->data[j], 16);
        }

        mcnp_str_write(&s, ctx->out);
        ctx->transforms_written++;
    }

    mcnp_str_comment(&cs, "END Transform Cards");
    mcnp_str_comment(&cs, "");
    mcnp_str_write(&cs, ctx->out);

    free(flatten);
    free(surf_offsets);

    return 0;
}
