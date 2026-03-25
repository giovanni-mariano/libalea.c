// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/**
 * @file mcnp_conversion.c
 * @brief MCNP to CSG system conversion functions
 *
 * Converts MCNP geometry, materials, and transforms into a alea_system_t
 * using automatic deduplication and dynamic array growth.
 */
#include "alea.h"
#include "mcnp/parser/mcnp_parser.h"
#include "mcnp/mcnp_model.h"
#include "surface_conv.h"
#include "cell_conv.h"
#include "core/alea_cell_complement.h"
#include "core/alea_universe.h"
#include "core/alea_materials.h"
#include "util/compat.h"
#include "util/alea_log.h"
#include "util/alea_bitset.h"
#include "core/alea_system.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>


// ============================================================================
// MATERIAL PARSING
// ============================================================================

/**
 * Parse material definition string into nuclides
 * Format: ZAID.LIB fraction ZAID.LIB fraction ...
 * e.g., "92235.80c -0.04 92238.80c -0.96"
 *
 * Uses alea_mat_add_nuclide() from the materials API.
 */
static int parse_material_definition(alea_material_t* mat, const char* definition) {
    if (!mat || !definition) return -1;

    const char* p = definition;

    /* Default to atom fractions */
    mat->is_weight_fraction = false;

    while (*p) {
        /* Skip whitespace */
        while (*p && isspace(*p)) p++;
        if (!*p) break;

        /* Parse ZAID */
        int zaid = 0;
        while (*p && isdigit(*p)) {
            zaid = zaid * 10 + (*p - '0');
            p++;
        }
        if (zaid == 0) break;

        /* Parse library suffix (optional) */
        char library[16] = {0};
        if (*p == '.') {
            p++;
            int li = 0;
            while (*p && !isspace(*p) && li < 15) {
                library[li++] = *p++;
            }
        }

        /* Skip whitespace */
        while (*p && isspace(*p)) p++;

        /* Parse fraction */
        double fraction = 0.0;
        bool negative = false;
        if (*p == '-') {
            negative = true;
            p++;
        } else if (*p == '+') {
            p++;
        }

        char* endp;
        fraction = strtod(p, &endp);
        if (endp == p) break;  /* No number found */
        p = endp;

        if (negative) {
            mat->is_weight_fraction = true;
        }

        /* Add nuclide using the materials API */
        if (alea_mat_add_nuclide(mat, zaid, library[0] ? library : NULL, fraction) < 0) {
            return -1;
        }
    }

    return alea_vec_count(&mat->nuclides) > 0 ? 0 : -1;
}

/**
 * Convert MCNP material to alea_material_t
 *
 * Materials are stored in a flat array in the system. We initialize each
 * material in place with proper capacity for nuclides.
 */
static int convert_material(alea_system_t* sys, const mcnp_material_t* mcnp_mat) {
    if (!sys || !mcnp_mat) return -1;

    ALEA_LOG_INFO("Converting material M%d", mcnp_mat->material_id);
    ALEA_LOG_INFO("Material count: %zu\n", alea_vec_count(&sys->materials));

    /* Register material via public API */
    int mat_idx = alea_add_material(sys, mcnp_mat->material_id);
    if (mat_idx < 0) return -1;
    alea_material_t* mat = &sys->materials.data[mat_idx];

    /* Parse the material definition */
    if (mcnp_mat->definition) {
        if (parse_material_definition(mat, mcnp_mat->definition) < 0) {
            ALEA_LOG_WARN("Warning: Failed to parse material M%d definition\n",
                    mcnp_mat->material_id);
        }
    }

    /* Copy comments */
    if (mcnp_mat->comments) {
        mat->comments = alea_strdup(mcnp_mat->comments);
    }

    ALEA_LOG_INFO("Converted material M%d: %zu nuclides (%s fractions)\n",
           mat->material_id, alea_vec_count(&mat->nuclides),
           mat->is_weight_fraction ? "weight" : "atom");

    return 0;
}

// ============================================================================
// TRANSFORM PARSING  
// ============================================================================

/**
 * Parse TRn card definition
 * Format: ox oy oz [b1 b2 b3 b4 b5 b6 b7 b8 b9] [M]
 * M=1 means degrees for rotation matrix
 */
static int parse_transform_definition(alea_transform_t* tr, const char* definition, int is_star) {
    if (!tr || !definition) return -1;
    
    const char* p = definition;
    double values[13];  /* Up to 12 values + possible M flag */
    int count = 0;
    
    while (*p && count < 13) {
        while (*p && isspace(*p)) p++;
        if (!*p) break;
        
        char* endp;
        double val = strtod(p, &endp);
        if (endp == p) break;
        
        values[count++] = val;
        p = endp;
    }
    
    if (count < 3) return -1;  /* Need at least ox oy oz */
    
    /* Check for M flag at end */
    int m_flag = 0;
    if (count == 4 || count == 13) {
        int last = (int)values[count-1];
        if (last == 1) {
            m_flag = 1;
            count--;
        }
    }
    
    tr->value_count = (count <= 3) ? 3 : 12;
    tr->degrees = is_star || m_flag;
    
    for (int i = 0; i < count && i < 12; i++) {
        tr->data[i] = values[i];
    }
    for (int i = count; i < 12; i++) {
        tr->data[i] = 0.0;
    }
    
    return 0;
}


// ============================================================================
// VACUUM BOUNDARY DETECTION
// ============================================================================

/**
 * @brief Detect and mark vacuum boundary surfaces
 *
 * Uses raycast from "infinity" to find the graveyard cell (IMP:N=0, material=0)
 * that represents the outside world. Surfaces bounding this cell are marked
 * as ALEA_BOUNDARY_VACUUM.
 *
 * @param sys CSG system with cells already converted
 * @return Number of surfaces marked as vacuum, or -1 on error
 */
static int detect_vacuum_boundaries(alea_system_t* sys, const mcnp_model_t* model) {
    if (!sys) return -1;

    /* Test points at "infinity" in 6 directions (9.9e5 cm = ~10 km) */
    static const double FAR = 9.9e5;
    static const double test_points[][3] = {
        { FAR, 0, 0}, {-FAR, 0, 0},
        {0,  FAR, 0}, {0, -FAR, 0},
        {0, 0,  FAR}, {0, 0, -FAR}
    };

    int graveyard_cell_id = -1;
    int graveyard_cell_idx = -1;

    /* Use alea_identify_cell_at_point which works correctly */
    for (int i = 0; i < 6 && graveyard_cell_idx < 0; i++) {
        int cell_idx = alea_identify_cell_at_point(sys,
                                                   test_points[i][0],
                                                   test_points[i][1],
                                                   test_points[i][2]);

        if (cell_idx >= 0 && (size_t)cell_idx < alea_vec_count(&sys->cells)) {
            alea_cell_entry_t* cell = &sys->cells.data[cell_idx];

            /* Check if this is graveyard: material=0 and IMP:N=0 */
            const mcnp_cell_params_t* mp = model ? mcnp_cell_params_const(model, (size_t)cell_idx) : NULL;
            if (cell->material_id == 0 && mp && mp->has_imp_n && mp->imp_n == 0.0) {
                graveyard_cell_id = cell->mc_cell_id;
                graveyard_cell_idx = cell_idx;
            }
        }
    }

    /* Second pass: if no explicit IMP:N=0, use any void cell at infinity */
    if (graveyard_cell_idx < 0) {
        for (int i = 0; i < 6 && graveyard_cell_idx < 0; i++) {
            int cell_idx = alea_identify_cell_at_point(sys,
                                                       test_points[i][0],
                                                       test_points[i][1],
                                                       test_points[i][2]);

            if (cell_idx >= 0 && (size_t)cell_idx < alea_vec_count(&sys->cells)) {
                alea_cell_entry_t* cell = &sys->cells.data[cell_idx];
                if (cell->material_id == 0) {
                    graveyard_cell_id = cell->mc_cell_id;
                    graveyard_cell_idx = cell_idx;
                    ALEA_LOG_DEBUG("Using void cell %d at infinity as graveyard (no explicit IMP:N=0)",
                                  cell->mc_cell_id);
                }
            }
        }
    }

    if (graveyard_cell_idx < 0) {
        ALEA_LOG_DEBUG("No graveyard cell found at infinity");
        return 0;
    }

    ALEA_LOG_INFO("Detected graveyard cell %d", graveyard_cell_id);

    /* Build surface index for the graveyard cell */
    alea_build_cell_surface_index(sys);

    alea_cell_entry_t* graveyard = &sys->cells.data[graveyard_cell_idx];

    if (!graveyard->surface_indices || graveyard->surface_index_count == 0) {
        ALEA_LOG_WARN("Graveyard cell has no surfaces");
        return 0;
    }

    /* Mark all surfaces of graveyard cell as vacuum boundaries */
    int marked = 0;
    for (size_t i = 0; i < graveyard->surface_index_count; i++) {
        uint32_t surf_idx = graveyard->surface_indices[i];
        if (surf_idx < alea_vec_count(&sys->surfaces)) {
            alea_surface_entry_t* surf = &sys->surfaces.data[surf_idx];
            /* Only mark if not already a special boundary */
            if (surf->boundary_type == ALEA_BOUNDARY_TRANSMISSIVE) {
                surf->boundary_type = ALEA_BOUNDARY_VACUUM;
                marked++;
                ALEA_LOG_DEBUG("Marked surface %d as vacuum boundary",
                             surf->mc_surface_id);
            }
        }
    }

    ALEA_LOG_INFO("Marked %d surfaces as vacuum boundaries", marked);
    return marked;
}


// ============================================================================
// RESERVE HINTS FOR VECTORS
// ============================================================================

/**
 * Count operand tokens in an MCNP geometry string.
 *
 * Only surface/cell references create CSG leaf nodes; the ':' union operator
 * is a pure delimiter and must NOT be counted.  For N operands the CSG binary
 * tree has exactly 2N-1 nodes (N leaves + N-1 internal), so the caller
 * multiplies this count by 2 to get a tight upper bound.
 *
 * Excluded from count: whitespace, '(', ')', ':'.
 */
static size_t count_tokens(const char* str) {
    if (!str || !*str) return 0;

    size_t tokens = 0;
    bool in_token = false;

    for (const char* p = str; *p; p++) {
        if (*p == ' ' || *p == '\t' || *p == '(' || *p == ')' || *p == ':') {
            in_token = false;
        } else if (!in_token) {
            in_token = true;
            tokens++;
        }
    }
    return tokens;
}

/**
 * Pre-allocate vectors based on MCNP context to avoid reallocations.
 */
static void reserve_from_mcnp(alea_system_t* sys, const mcnp_context_t* mcnp) {
    // Cells: tokens (excluding parens) become primitive refs + operation nodes
    // N tokens → N primitive refs + (N-1) operators ≈ 2N nodes
    size_t estimated_nodes = 0;
    for (size_t i = 0; i < mcnp->cell_count; i++) {
        size_t tokens = count_tokens(mcnp->cells[i]->geometry_definition);
        estimated_nodes += tokens * 2;
    }

    alea_vec_reserve(&sys->nodes, estimated_nodes, alea_node_t);
    alea_vec_reserve(&sys->surfaces, mcnp->surface_count, alea_surface_entry_t);
    alea_vec_reserve(&sys->primitives, mcnp->surface_count, alea_primitive_entry_t);
    alea_vec_reserve(&sys->cells, mcnp->cell_count, alea_cell_entry_t);
    alea_vec_reserve(&sys->transforms, mcnp->transform_count, alea_transform_t);
    alea_vec_reserve(&sys->materials, mcnp->material_count, alea_material_t);

    ALEA_LOG_DEBUG("Reserved: %zu nodes for %zu cells\n",
           estimated_nodes, mcnp->cell_count);
}

// ============================================================================
// COMPLETE FILE CONVERSION
// ============================================================================

mcnp_model_t* mcnp_convert_to_model(const char* filename) {
    if (!filename) return NULL;

    // Parse MCNP file
    mcnp_context_t* mcnp = NULL;
    if (!mcnp_parse_file(filename, &mcnp) || !mcnp) {
        ALEA_LOG_ERROR("Failed to parse MCNP file: %s\n", filename);
        return NULL;
    }

    ALEA_LOG_DEBUG("Parsed MCNP file: %zu surfaces, %zu cells, %zu materials\n",
           mcnp->surface_count, mcnp->cell_count, mcnp->material_count);

    // Create CSG system with automatic sizing
    alea_system_t* sys = alea_system_create();
    if (!sys) {
        ALEA_LOG_ERROR("Failed to create CSG system\n");
        mcnp_context_destroy(mcnp);
        return NULL;
    }

    // Create model early, register hooks so cell additions auto-grow params
    mcnp_model_t* model = calloc(1, sizeof(mcnp_model_t));
    if (!model) {
        alea_system_destroy(sys);
        mcnp_context_destroy(mcnp);
        return NULL;
    }
    model->sys = sys;
    model->owns_sys = 1;
    model->export_config = MCNP_EXPORT_CONFIG_DEFAULT;
    mcnp_model_reserve_params(model, mcnp->cell_count);
    mcnp_model_register_hooks(model);

    // Pre-allocate vectors based on parsed counts to avoid reallocations
    reserve_from_mcnp(sys, mcnp);

   // Convert all transforms (TRn cards)
    ALEA_LOG_INFO("\nConverting transforms...\n");
    for (size_t i = 0; i < mcnp->transform_count; i++) {
        mcnp_transform_t* tr = mcnp->transforms[i];
        alea_transform_t alea_tr;
        memset(&alea_tr, 0, sizeof(alea_tr));
        alea_tr.transform_id = tr->transform_id;
        alea_tr.from_inline = 0;

        if (tr->definition && parse_transform_definition(&alea_tr, tr->definition, tr->is_star) == 0) {
            alea_add_transform(sys, tr->transform_id, alea_tr.data, alea_tr.value_count, alea_tr.degrees);
        } else {
            ALEA_LOG_WARN("Warning: Failed to parse transform TR%d\n", tr->transform_id);
        }
    }

    // Set inline transform IDs to start above all TRn card IDs
    alea_finalize_transform_ids(sys);

    // Convert all surfaces (with automatic deduplication)
    ALEA_LOG_INFO("\nConverting surfaces...\n");

    // Build seen-array for O(1) duplicate surface ID detection
    int max_surf_id = 0;
    for (size_t i = 0; i < mcnp->surface_count; i++) {
        if (mcnp->surfaces[i]->surface_id > max_surf_id)
            max_surf_id = mcnp->surfaces[i]->surface_id;
    }
    alea_bitset_t surf_seen = alea_bitset_create((size_t)max_surf_id + 1);

    for (size_t i = 0; i < mcnp->surface_count; i++) {
        if (g_alea_interrupted) goto interrupted;

        int sid = mcnp->surfaces[i]->surface_id;
        if (surf_seen.words && alea_bitset_test(&surf_seen, sid)) {
            ALEA_LOG_WARN("Duplicate surface id=%d, skipping", sid);
            continue;
        }
        if (surf_seen.words) alea_bitset_set(&surf_seen, sid);

        if (alea_convert_surface(sys, mcnp->surfaces[i]) < 0) {
            ALEA_LOG_ERROR("Failed to convert surface %d",
                    mcnp->surfaces[i]->surface_id);
            alea_bitset_destroy(&surf_seen);
            mcnp_model_destroy(model);
            mcnp_context_destroy(mcnp);
            return NULL;
        }
    }
    alea_bitset_destroy(&surf_seen);

    // Build surface lookup table
    alea_build_surface_lookup(sys);

    // Convert all materials (before cells, so cell conversion can resolve material indices)
    ALEA_LOG_INFO("\nConverting materials...\n");
    for (size_t i = 0; i < mcnp->material_count; i++) {
        if (g_alea_interrupted) goto interrupted;
        if (convert_material(sys, mcnp->materials[i]) < 0) {
            ALEA_LOG_WARN("Warning: Failed to convert material M%d\n",
                    mcnp->materials[i]->material_id);
        }
    }

    // Convert all cells (hooks auto-grow model->cell_params)
    ALEA_LOG_INFO("\nConverting cells...\n");
    for (size_t i = 0; i < mcnp->cell_count; i++) {
        if (g_alea_interrupted) goto interrupted;
        if (alea_convert_cell(sys, mcnp->cells[i], model) == UINT32_MAX) {
            ALEA_LOG_WARN("Warning: Failed to convert cell %d\n",
                    mcnp->cells[i]->cell_id);
        }
    }

    if (!alea_vec_empty(&sys->cell_refs)) {
        ALEA_LOG_INFO("\nResolving cell complement references...\n");
        int resolved = alea_resolve_cell_complements(sys);
        if (resolved < 0) {
            ALEA_LOG_ERROR("Error resolving cell complements\n");
        }
    }

    // Resolve LIKE BUT cells (must happen before TRCL transforms)
    ALEA_LOG_INFO("\nResolving LIKE BUT cells...\n");
    int like_count = alea_resolve_like_cells(sys, model);
    if (like_count < 0) {
        ALEA_LOG_WARN("Warning: Error resolving LIKE cells\n");
    }

    // Apply TRCL transforms to cells
    ALEA_LOG_INFO("\nApplying TRCL transforms...\n");
    int trcl_count = alea_apply_trcl_transforms(sys, model);
    if (trcl_count < 0) {
        ALEA_LOG_WARN("Warning: Error applying TRCL transforms\n");
    }

    // Compute lattice pitch and lower_left from cell bounding boxes
    for (size_t i = 0; i < alea_vec_count(&sys->cells); i++) {
        alea_cell_entry_t* cell = &sys->cells.data[i];
        if (cell->lat_type == 0 || !cell->lat_fill) continue;

        alea_bbox_t bb = sys->nodes.data[cell->root_node_id].bbox;
        cell->lat_pitch[0] = bb.max_x - bb.min_x;
        cell->lat_pitch[1] = bb.max_y - bb.min_y;
        cell->lat_pitch[2] = bb.max_z - bb.min_z;
        cell->lat_lower_left[0] = bb.min_x + cell->lat_fill_dims[0] * cell->lat_pitch[0];
        cell->lat_lower_left[1] = bb.min_y + cell->lat_fill_dims[2] * cell->lat_pitch[1];
        cell->lat_lower_left[2] = bb.min_z + cell->lat_fill_dims[4] * cell->lat_pitch[2];

        ALEA_LOG_DEBUG("Lattice cell %d: pitch=(%.4f, %.4f, %.4f) lower_left=(%.4f, %.4f, %.4f)",
                     cell->mc_cell_id,
                     cell->lat_pitch[0], cell->lat_pitch[1], cell->lat_pitch[2],
                     cell->lat_lower_left[0], cell->lat_lower_left[1], cell->lat_lower_left[2]);
    }

    // Detect vacuum boundaries from graveyard cell
    ALEA_LOG_INFO("\nDetecting vacuum boundaries...\n");
    int vacuum_count = detect_vacuum_boundaries(sys, model);
    if (vacuum_count > 0) {
        ALEA_LOG_INFO("Marked %d vacuum boundary surfaces\n", vacuum_count);
    }

    // Print statistics
    alea_system_print_stats(sys);

    // Cleanup MCNP context
    mcnp_context_destroy(mcnp);

    return model;

interrupted:
    alea_set_error_detail(ALEA_ERR_INTERRUPTED, "MCNP conversion interrupted");
    mcnp_model_destroy(model);
    mcnp_context_destroy(mcnp);
    return NULL;
}


