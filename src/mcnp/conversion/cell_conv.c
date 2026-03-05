// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/**
 * @file cell_conv.c
 * @brief MCNP cell to CSG tree conversion - Updated for flattened architecture
 * 
 * This module converts MCNP cell geometry expressions to CSG trees.
 */


#include "alea.h"
#include "cell_conv.h"
#include "surface_conv.h"
#include "mcnp/parser/geom_parser.h"
#include "core/alea_universe.h"
#include "util/alea_log.h"
#include "util/compat.h"
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <stdio.h>


// ============================================================================
// CELL PARAMETERS PARSER
// ============================================================================

/**
 * @brief Case-insensitive string prefix match
 */
static int match_prefix(const char* str, const char* prefix) {
    while (*prefix) {
        if (tolower((unsigned char)*str) != tolower((unsigned char)*prefix)) {
            return 0;
        }
        str++;
        prefix++;
    }
    return 1;
}

/**
 * @brief Skip whitespace
 */
static const char* skip_ws(const char* p) {
    while (*p && isspace((unsigned char)*p)) p++;
    return p;
}

/**
 * @brief Parse a double value, handling optional '=' prefix
 */
static double parse_double_value(const char** cursor, int* ok) {
    const char* p = *cursor;
    p = skip_ws(p);
    if (*p == '=') p++;
    p = skip_ws(p);

    char* end;
    double val = strtod(p, &end);
    if (end > p) {
        *cursor = end;
        if (ok) *ok = 1;
        return val;
    }
    if (ok) *ok = 0;
    return 0.0;
}

/**
 * @brief Parse an integer value, handling optional '=' prefix
 */
static int parse_int_value(const char** cursor, int* ok) {
    const char* p = *cursor;
    p = skip_ws(p);
    if (*p == '=') p++;
    p = skip_ws(p);

    char* end;
    long val = strtol(p, &end, 10);
    if (end > p) {
        *cursor = end;
        if (ok) *ok = 1;
        return (int)val;
    }
    if (ok) *ok = 0;
    return 0;
}

/**
 * @brief Parse parenthesized transform values
 *
 * Helper to parse values like (t), (ox oy oz), or (ox oy oz b1..b9)
 * Returns number of values parsed, stores in values array.
 * Updates cursor past closing paren.
 */
static int parse_transform_parens(const char** cursor, double* values, int max_values,
                                   int* out_is_reference) {
    const char* p = *cursor;
    p = skip_ws(p);

    if (*p != '(') {
        return 0;
    }
    p++;
    p = skip_ws(p);

    int count = 0;
    int has_decimal = 0;
    char* end;

    while (count < max_values && *p && *p != ')') {
        const char* val_start = p;
        double val = strtod(p, &end);
        if (end == p) break;

        for (const char* c = val_start; c < end; c++) {
            if (*c == '.') {
                has_decimal = 1;
                break;
            }
        }

        values[count++] = val;
        p = end;
        p = skip_ws(p);
    }

    while (*p && *p != ')') p++;
    if (*p == ')') p++;

    *cursor = p;

    // Detect if single integer-like value (transform reference)
    *out_is_reference = (count == 1 && !has_decimal);

    return count;
}

/**
 * @brief Check if string looks like a lattice range (contains ':')
 */
static int is_lattice_fill(const char* p) {
    // Look ahead to see if we have "num:num" pattern
    // Skip leading whitespace
    while (*p && isspace(*p)) p++;

    // Check for optional negative sign and digits
    if (*p == '-' || *p == '+') p++;
    if (!isdigit(*p)) return 0;

    // Skip digits
    while (isdigit(*p)) p++;

    // Check for colon
    return (*p == ':');
}

/**
 * @brief Parse lattice FILL array format
 *
 * Format: FILL=i1:i2 j1:j2 [k1:k2] u1 u2 u3 ...
 * For 2D lattice: FILL=i1:i2 j1:j2 0:0 u1 u2 ...
 */
static void parse_lattice_fill(const char** cursor, alea_cell_params_t* params) {
    const char* p = *cursor;
    char* end;

    // Parse dimension ranges
    int dims[6] = {0, 0, 0, 0, 0, 0};  // imin, imax, jmin, jmax, kmin, kmax
    int dim_idx = 0;

    // Parse up to 3 ranges (i, j, k)
    while (dim_idx < 6) {
        const char* range_start = p;
        p = skip_ws(p);

        // Parse imin/jmin/kmin
        long val = strtol(p, &end, 10);
        if (end == p) break;
        p = end;

        p = skip_ws(p);
        if (*p != ':') {
            // Not a range - this is the first universe ID
            // Rewind to before this number so it gets parsed as a universe
            p = range_start;
            break;
        }
        // It's a valid range - store the min value
        dims[dim_idx++] = (int)val;
        p++;  // skip ':'

        // Parse imax/jmax/kmax
        val = strtol(p, &end, 10);
        if (end == p) break;
        dims[dim_idx++] = (int)val;
        p = end;
    }

    // Must have at least i range (2 values)
    if (dim_idx < 2) {
        *cursor = p;
        return;
    }

    // Store dimensions
    for (int i = 0; i < 6; i++) {
        params->lat_fill_dims[i] = dims[i];
    }

    // Calculate total count (use uint64_t to prevent int overflow)
    int ni = dims[1] - dims[0] + 1;
    int nj = (dim_idx >= 4) ? (dims[3] - dims[2] + 1) : 1;
    int nk = (dim_idx >= 6) ? (dims[5] - dims[4] + 1) : 1;
    if (ni <= 0 || nj <= 0 || nk <= 0) {
        *cursor = p;
        return;
    }
    uint64_t total = (uint64_t)ni * (uint64_t)nj * (uint64_t)nk;

    if (total == 0 || total > 1000000) {  // Sanity check
        *cursor = p;
        return;
    }

    // Allocate array for universes
    params->lat_fill = malloc(total * sizeof(int));
    if (!params->lat_fill) {
        *cursor = p;
        return;
    }

    // Parse universe IDs
    size_t count = 0;
    while (count < total) {
        p = skip_ws(p);
        if (!*p) break;

        // Stop at keywords
        if (isalpha(*p)) break;

        long univ = strtol(p, &end, 10);
        if (end == p) break;

        params->lat_fill[count++] = (int)univ;
        p = end;
    }

    params->lat_fill_count = count;
    params->lat_fill_is_array = 1;
    params->has_fill = 1;

    *cursor = p;
}

/**
 * @brief Parse FILL parameter which can have complex syntax
 *
 * Formats supported:
 *   FILL=n              - Simple fill with universe n
 *   FILL=n (t)          - Fill with universe n, apply transform TRt (reference)
 *   FILL=n (ox oy oz)   - Fill with inline translation
 *   FILL=n (ox oy oz b1 b2 b3 b4 b5 b6 b7 b8 b9) - Fill with inline full transform
 *   FILL=i1:i2 j1:j2 [k1:k2] u1 u2 ... - Lattice fill array
 *   *FILL=n (...)       - Same as above but angles in degrees
 *
 * Detection logic for parentheses content:
 *   - 1 value that looks like integer: transform reference
 *   - 3 values: inline translation
 *   - 12 values: inline full transform
 */
static void parse_fill_param(const char** cursor, alea_cell_params_t* params, int degrees) {
    const char* p = *cursor;
    p = skip_ws(p);
    if (*p == '=') p++;
    p = skip_ws(p);

    params->fill_transform_degrees = degrees;

    // Check for lattice fill format (contains ':')
    if (is_lattice_fill(p)) {
        parse_lattice_fill(&p, params);
        *cursor = p;
        return;
    }

    // Parse universe number
    char* end;
    long universe = strtol(p, &end, 10);
    if (end > p) {
        params->fill_universe = (int)universe;
        params->has_fill = 1;
        p = end;
    } else {
        *cursor = p;
        return;
    }

    // Check for optional transform in parentheses
    p = skip_ws(p);
    if (*p == '(') {
        p++;
        p = skip_ws(p);

        // Parse all values inside parentheses
        double values[12];
        int count = 0;
        int has_decimal = 0;

        while (count < 12 && *p && *p != ')') {
            const char* val_start = p;
            double val = strtod(p, &end);
            if (end == p) break;  // No more numbers

            // Check if this value has a decimal point
            for (const char* c = val_start; c < end; c++) {
                if (*c == '.') {
                    has_decimal = 1;
                    break;
                }
            }

            values[count++] = val;
            p = end;
            p = skip_ws(p);
        }

        // Skip to closing paren
        while (*p && *p != ')') p++;
        if (*p == ')') p++;

        // Interpret the values
        if (count == 1 && !has_decimal) {
            // Single integer-like value: transform reference
            params->fill_transform = (int)values[0];
            params->fill_transform_inline = 0;
            params->fill_transform_count = 0;
        } else if (count == 3 || count == 12) {
            // Inline transform
            params->fill_transform_inline = 1;
            params->fill_transform_count = count;
            for (int i = 0; i < count; i++) {
                params->fill_transform_data[i] = values[i];
            }
        } else if (count > 0) {
            // Unexpected count - store what we have
            params->fill_transform_inline = 1;
            params->fill_transform_count = count;
            for (int i = 0; i < count && i < 12; i++) {
                params->fill_transform_data[i] = values[i];
            }
        }
    }

    *cursor = p;
}

/**
 * @brief Parse TRCL parameter
 *
 * Formats supported:
 *   TRCL=n              - Apply transform TRn
 *   TRCL=(ox oy oz)     - Inline translation
 *   TRCL=(ox oy oz b1 b2 b3 b4 b5 b6 b7 b8 b9) - Inline full transform
 *   *TRCL=...           - Same but angles in degrees
 */
static void parse_trcl_param(const char** cursor, alea_cell_params_t* params, int degrees) {
    const char* p = *cursor;
    p = skip_ws(p);
    if (*p == '=') p++;
    p = skip_ws(p);

    params->trcl_degrees = degrees;
    params->has_trcl = 1;

    // Check if it starts with parenthesis (inline transform) or number (TR reference)
    if (*p == '(') {
        // Inline transform
        double values[12];
        int is_reference;
        int count = parse_transform_parens(&p, values, 12, &is_reference);

        if (count > 0) {
            if (is_reference) {
                // Single integer - transform reference in parens: TRCL=(n)
                params->trcl = (int)values[0];
                params->trcl_inline = 0;
                params->trcl_count = 0;
            } else {
                // Inline transform data
                params->trcl_inline = 1;
                params->trcl_count = count;
                for (int i = 0; i < count && i < 12; i++) {
                    params->trcl_data[i] = values[i];
                }
            }
        }
    } else {
        // Direct transform reference: TRCL=n
        char* end;
        long tr = strtol(p, &end, 10);
        if (end > p) {
            params->trcl = (int)tr;
            params->trcl_inline = 0;
            params->trcl_count = 0;
            p = end;
        }
    }

    *cursor = p;
}

int parse_cell_parameters(const char* params_str, alea_cell_params_t* out_params,
                          int cell_id) {
    if (!out_params) return -1;

    // Initialize output structure
    memset(out_params, 0, sizeof(alea_cell_params_t));
    out_params->imp_n = 1.0;  // Default importance
    out_params->imp_p = 1.0;
    out_params->imp_e = 1.0;

    if (!params_str || !*params_str) {
        return 0;  // No parameters to parse
    }

    const char* p = params_str;
    int ok;
    
    while (*p) {
        p = skip_ws(p);
        if (!*p) break;
        
        // IMP:X[,Y[,Z]]=value - Importance (supports IMP:N=1, IMP:N,P=1, IMP:N,P,E=1)
        if (match_prefix(p, "IMP:")) {
            p += 4;  // Skip "IMP:"

            // Parse particle types (N, P, E separated by commas)
            int set_n = 0, set_p = 0, set_e = 0;
            while (*p && *p != '=' && !isspace((unsigned char)*p)) {
                char c = toupper((unsigned char)*p);
                if (c == 'N') set_n = 1;
                else if (c == 'P') set_p = 1;
                else if (c == 'E') set_e = 1;
                p++;
                if (*p == ',') p++;  // Skip comma separator
            }

            // Parse the value
            double imp_val = parse_double_value(&p, &ok);
            if (!ok) {
                ALEA_LOG_WARN("Cell %d: invalid value for IMP parameter", cell_id);
                continue;
            }

            // Apply to selected particle types
            if (set_n) {
                if (out_params->has_imp_n)
                    ALEA_LOG_WARN("Cell %d: IMP:N specified multiple times, using last value", cell_id);
                out_params->imp_n = imp_val;
                out_params->has_imp_n = 1;
            }
            if (set_p) {
                if (out_params->has_imp_p)
                    ALEA_LOG_WARN("Cell %d: IMP:P specified multiple times, using last value", cell_id);
                out_params->imp_p = imp_val;
                out_params->has_imp_p = 1;
            }
            if (set_e) {
                if (out_params->has_imp_e)
                    ALEA_LOG_WARN("Cell %d: IMP:E specified multiple times, using last value", cell_id);
                out_params->imp_e = imp_val;
                out_params->has_imp_e = 1;
            }
        }
        // U= - Universe
        else if (match_prefix(p, "U=") || (match_prefix(p, "U") && (p[1] == '=' || isspace((unsigned char)p[1])))) {
            p++;  // Skip 'U'
            if (out_params->has_universe)
                ALEA_LOG_WARN("Cell %d: U specified multiple times, using last value", cell_id);
            out_params->universe_id = parse_int_value(&p, &ok);
            if (!ok) ALEA_LOG_WARN("Cell %d: invalid value for U parameter", cell_id);
            out_params->has_universe = 1;
        }
        // *FILL= - Fill with universe (angles in degrees)
        else if (match_prefix(p, "*FILL")) {
            p += 5;
            if (out_params->has_fill)
                ALEA_LOG_WARN("Cell %d: FILL specified multiple times, using last value", cell_id);
            parse_fill_param(&p, out_params, 1);  // degrees=1
        }
        // FILL= - Fill with universe
        else if (match_prefix(p, "FILL")) {
            p += 4;
            if (out_params->has_fill)
                ALEA_LOG_WARN("Cell %d: FILL specified multiple times, using last value", cell_id);
            parse_fill_param(&p, out_params, 0);  // degrees=0
        }
        // LAT= - Lattice type
        else if (match_prefix(p, "LAT")) {
            p += 3;
            if (out_params->has_lat)
                ALEA_LOG_WARN("Cell %d: LAT specified multiple times, using last value", cell_id);
            out_params->lat_type = parse_int_value(&p, &ok);
            if (!ok) ALEA_LOG_WARN("Cell %d: invalid value for LAT parameter", cell_id);
            out_params->has_lat = 1;
        }
        // Simple keyword=value params (VOL, TMP, PWT, NONU, PD, ELPT, UNC, BFLCL)
        #define MCNP_PARSE_double parse_double_value
        #define MCNP_PARSE_int    parse_int_value
        #define X_PARSE(name, type, kw, prec) \
            else if (match_prefix(p, kw)) { \
                p += sizeof(kw) - 1; \
                if (out_params->has_##name) \
                    ALEA_LOG_WARN("Cell %d: " kw " specified multiple times, using last value", cell_id); \
                out_params->name = (type)MCNP_PARSE_##type(&p, &ok); \
                if (!ok) ALEA_LOG_WARN("Cell %d: invalid value for " kw " parameter", cell_id); \
                out_params->has_##name = 1; \
            }
        MCNP_CELL_SIMPLE_PARAMS(X_PARSE)
        #undef X_PARSE
        #undef MCNP_PARSE_double
        #undef MCNP_PARSE_int
        // *TRCL= - Cell transformation (angles in degrees)
        else if (match_prefix(p, "*TRCL")) {
            p += 5;
            if (out_params->has_trcl)
                ALEA_LOG_WARN("Cell %d: TRCL specified multiple times, using last value", cell_id);
            parse_trcl_param(&p, out_params, 1);  // degrees=1
        }
        // TRCL= - Cell transformation
        else if (match_prefix(p, "TRCL")) {
            p += 4;
            if (out_params->has_trcl)
                ALEA_LOG_WARN("Cell %d: TRCL specified multiple times, using last value", cell_id);
            parse_trcl_param(&p, out_params, 0);  // degrees=0
        }
        // MAT= - Material override (for LIKE BUT)
        else if (match_prefix(p, "MAT")) {
            p += 3;
            if (out_params->has_mat)
                ALEA_LOG_WARN("Cell %d: MAT specified multiple times, using last value", cell_id);
            out_params->material_id = parse_int_value(&p, &ok);
            if (!ok) ALEA_LOG_WARN("Cell %d: invalid value for MAT parameter", cell_id);
            out_params->has_mat = 1;
        }
        // RHO= - Density override (for LIKE BUT)
        else if (match_prefix(p, "RHO")) {
            p += 3;
            if (out_params->has_rho)
                ALEA_LOG_WARN("Cell %d: RHO specified multiple times, using last value", cell_id);
            out_params->density = parse_double_value(&p, &ok);
            if (!ok) ALEA_LOG_WARN("Cell %d: invalid value for RHO parameter", cell_id);
            out_params->has_rho = 1;
        }
        // Skip unknown token
        else {
            // Advance to next whitespace or end
            while (*p && !isspace((unsigned char)*p)) p++;
        }
    }
    
    return 0;
}


// ============================================================================
// CELL CONVERSION

uint32_t alea_convert_cell(alea_system_t* sys, const mcnp_cell_t* cell,
                           mcnp_model_t* model) {
    if (!sys || !cell) return UINT32_MAX;

    // Parse the geometry expression using geometry parser
    // This converts the MCNP boolean expression (e.g., "-1 : (2 -3)") into CSG tree
    int like_cell_id = 0;
    const char* but_clause = NULL;
    alea_node_id_t root = parse_mcnp_geometry_with_like(
        sys,
        cell->geometry_definition,
        cell->cell_id,
        &like_cell_id,
        &but_clause
    );

    // Handle LIKE cells - defer geometry resolution until all cells parsed
    bool is_like_cell = (root == ALEA_NODE_ID_LIKE_PENDING && like_cell_id > 0);

    if (root == ALEA_NODE_ID_INVALID) {
        ALEA_LOG_ERROR("Failed to parse geometry for cell %d", cell->cell_id);
        return UINT32_MAX;
    }

    // For LIKE cells, set root to INVALID until resolution
    if (is_like_cell) {
        root = ALEA_NODE_ID_INVALID;
    }

    // Parse cell parameters (U=, FILL=, IMP:N=, etc.)
    alea_cell_params_t params;
    parse_cell_parameters(cell->parameters, &params, cell->cell_id);

    // For LIKE cells, parse BUT clause parameters (override template values)
    if (is_like_cell && but_clause) {
        alea_cell_params_t but_params;
        memset(&but_params, 0, sizeof(but_params));
        parse_cell_parameters(but_clause, &but_params, cell->cell_id);

        // Merge BUT parameters into params (BUT takes priority)
        if (but_params.has_mat) {
            params.material_id = but_params.material_id;
            params.has_mat = 1;
        }
        if (but_params.has_rho) {
            params.density = but_params.density;
            params.has_rho = 1;
        }
        if (but_params.has_universe) {
            params.universe_id = but_params.universe_id;
            params.has_universe = 1;
        }
        if (but_params.has_fill) {
            params.fill_universe = but_params.fill_universe;
            params.fill_transform = but_params.fill_transform;
            params.fill_transform_inline = but_params.fill_transform_inline;
            params.fill_transform_degrees = but_params.fill_transform_degrees;
            params.fill_transform_count = but_params.fill_transform_count;
            for (int i = 0; i < 12; i++) {
                params.fill_transform_data[i] = but_params.fill_transform_data[i];
            }
            params.has_fill = 1;
        }
        if (but_params.has_trcl) {
            params.trcl = but_params.trcl;
            params.trcl_inline = but_params.trcl_inline;
            params.trcl_degrees = but_params.trcl_degrees;
            params.trcl_count = but_params.trcl_count;
            for (int i = 0; i < 12; i++) {
                params.trcl_data[i] = but_params.trcl_data[i];
            }
            params.has_trcl = 1;
        }
        if (but_params.has_imp_n) {
            params.imp_n = but_params.imp_n;
            params.has_imp_n = 1;
        }
        if (but_params.has_imp_p) {
            params.imp_p = but_params.imp_p;
            params.has_imp_p = 1;
        }
        if (but_params.has_imp_e) {
            params.imp_e = but_params.imp_e;
            params.has_imp_e = 1;
        }
        #define X_MERGE(name, type, kw, prec) \
            if (but_params.has_##name) { \
                params.name = but_params.name; \
                params.has_##name = 1; \
            }
        MCNP_CELL_SIMPLE_PARAMS(X_MERGE)
        #undef X_MERGE
        if (but_params.has_lat) {
            params.lat_type = but_params.lat_type;
            params.has_lat = 1;
        }
    }

    // For LIKE cells, use material/density from BUT clause (or keep 0 for resolution)
    int material_id = is_like_cell && params.has_mat ? params.material_id : cell->material_id;
    double density = is_like_cell && params.has_rho ? params.density : cell->density;

    // Resolve MCNP material ID to material index (auto-register if missing)
    int mat_index = ALEA_MATERIAL_VOID;
    if (material_id != 0) {
        mat_index = alea_find_material_by_id(sys, material_id);
        if (mat_index < 0) {
            mat_index = alea_add_material(sys, material_id);
        }
    }

    // Add cell to system (use _with_id for file loading - validation done after all cells loaded)
    int cell_idx = alea_add_cell_with_id(sys, cell->cell_id, root,
                                        mat_index, density,
                                        params.universe_id);
    if (cell_idx < 0) {
        return UINT32_MAX;
    }

    // Write core fields to cell entry
    alea_cell_entry_t* entry = &sys->cells.data[cell_idx];
    entry->fill_universe = params.fill_universe;
    entry->fill_transform = params.fill_transform;
    entry->lat_type = params.lat_type;

    // Temperature: convert MCNP MeV to Kelvin
    if (params.has_tmp) {
        entry->temperature = params.tmp / 8.617333262e-11;
        entry->has_temperature = 1;
    }

    // Copy comments from parsed cell
    if (cell->comments)
        entry->comments = alea_strdup(cell->comments);
    if (cell->inline_comment)
        entry->inline_comment = alea_strdup(cell->inline_comment);

    // Handle inline FILL transforms (registered in sys, ID stored in entry)
    if (params.fill_transform_inline && params.fill_transform_count > 0) {
        int tr_id = alea_add_inline_transform(
            sys,
            params.fill_transform_data,
            params.fill_transform_count,
            params.fill_transform_degrees
        );
        if (tr_id > 0) {
            entry->fill_transform = tr_id;
        } else {
            entry->fill_transform = 0;
            ALEA_LOG_WARN("Failed to add inline FILL transform for cell %d",
                    cell->cell_id);
        }
    }

    // Lattice FILL array (transfer ownership of array)
    if (params.lat_fill_is_array && params.lat_fill && params.lat_fill_count > 0) {
        for (int i = 0; i < 6; i++) {
            entry->lat_fill_dims[i] = params.lat_fill_dims[i];
        }
        entry->lat_fill = params.lat_fill;
        entry->lat_fill_count = params.lat_fill_count;
        // Don't free params.lat_fill - ownership transferred
    }

    // For simple FILL=N with LAT, synthesize a 1x1x1 lattice fill
    if (entry->lat_type > 0 && entry->fill_universe > 0 && !entry->lat_fill) {
        int* fill = malloc(sizeof(int));
        if (fill) {
            fill[0] = entry->fill_universe;
            entry->lat_fill = fill;
            entry->lat_fill_count = 1;
            // lat_fill_dims defaults to [0,0,0,0,0,0] giving 1x1x1
        }
    }

    // Write MCNP-specific fields to model params
    if (model) {
        mcnp_cell_params_t* mp = mcnp_cell_params(model, (size_t)cell_idx);
        if (mp) {
            // Particle importances
            mp->imp_n = params.imp_n;
            mp->imp_p = params.imp_p;
            mp->imp_e = params.imp_e;
            mp->has_imp_n = params.has_imp_n;
            mp->has_imp_p = params.has_imp_p;
            mp->has_imp_e = params.has_imp_e;

            // Simple keyword=value params (VOL, TMP, PWT, etc.)
            #define X_COPY(name, type, kw, prec) \
                mp->name = params.name; mp->has_##name = params.has_##name;
            MCNP_CELL_SIMPLE_PARAMS(X_COPY)
            #undef X_COPY

            // FILL transform original data (for export)
            mp->fill_transform_inline = params.fill_transform_inline;
            mp->fill_transform_degrees = params.fill_transform_degrees;
            mp->fill_transform_count = params.fill_transform_count;
            for (int i = 0; i < params.fill_transform_count && i < 12; i++) {
                mp->fill_transform_data[i] = params.fill_transform_data[i];
            }

            // TRCL - cell transformation
            mp->trcl = params.trcl;
            mp->trcl_degrees = params.trcl_degrees;
            mp->has_trcl = params.has_trcl;
            if (params.trcl_inline && params.trcl_count > 0) {
                int trcl_id = alea_add_inline_transform(
                    sys, params.trcl_data, params.trcl_count, params.trcl_degrees);
                if (trcl_id > 0) {
                    mp->trcl = trcl_id;
                } else {
                    mp->trcl = 0;
                    ALEA_LOG_WARN("Failed to add inline TRCL transform for cell %d",
                            cell->cell_id);
                }
            }
            mp->trcl_inline = params.trcl_inline;
            mp->trcl_count = params.trcl_count;
            for (int i = 0; i < params.trcl_count && i < 12; i++) {
                mp->trcl_data[i] = params.trcl_data[i];
            }

            // LIKE BUT - store template cell ID for later resolution
            mp->like_cell_id = like_cell_id;
            if (is_like_cell) {
                mp->has_mat = params.has_mat;
                mp->has_rho = params.has_rho;
            }
        }
    }

    sys->stats.cells_converted++;

    {
        char detail[256] = "";
        int pos = 0;
        if (params.has_universe) {
            pos += snprintf(detail + pos, sizeof(detail) - pos, " U=%d", params.universe_id);
        }
        if (params.has_fill) {
            pos += snprintf(detail + pos, sizeof(detail) - pos,
                params.fill_transform_degrees ? " *FILL=%d" : " FILL=%d", params.fill_universe);
            if (params.fill_transform_inline && params.fill_transform_count > 0) {
                pos += snprintf(detail + pos, sizeof(detail) - pos, " (inline -> TR%d)", entry->fill_transform);
            } else if (params.fill_transform > 0) {
                pos += snprintf(detail + pos, sizeof(detail) - pos, " (TR%d)", params.fill_transform);
            }
        }
        if (is_like_cell) {
            ALEA_LOG_INFO("Registered cell %d LIKE %d (pending resolution)%s", cell->cell_id, like_cell_id, detail);
        } else {
            ALEA_LOG_INFO("Converted cell %d to CSG tree (root node %u)%s", cell->cell_id, root, detail);
        }
    }
    
    // For LIKE cells, root is ALEA_NODE_ID_INVALID (pending resolution) but the
    // cell was created successfully. Return a non-error value.
    if (is_like_cell) return ALEA_NODE_ID_LIKE_PENDING;
    return root;
}

// ============================================================================
// TRCL TRANSFORM APPLICATION
// ============================================================================

int alea_apply_trcl_transforms(alea_system_t* sys, mcnp_model_t* model) {
    if (!sys || !model) return -1;

    int count = 0;

    for (size_t i = 0; i < alea_vec_count(&sys->cells); i++) {
        alea_cell_entry_t* cell = &sys->cells.data[i];
        const mcnp_cell_params_t* mp = mcnp_cell_params_const(model, i);

        /* Skip cells without TRCL */
        if (!mp || !mp->has_trcl) continue;
        if (cell->root_node_id == ALEA_NODE_ID_INVALID) continue;

        /* Build transform matrix from TRCL data */
        alea_matrix_t mat;

        if (mp->trcl > 0) {
            /* TRCL references a TRn card */
            const alea_transform_t* tr = alea_get_transform(sys, mp->trcl);
            if (!tr) {
                ALEA_LOG_WARN("Cell %d: TRCL=%d references unknown transform, skipping",
                            cell->mcnp_cell_id, mp->trcl);
                continue;
            }
            /* Use tr->cosines which has pre-computed direction cosines */
            alea_matrix_from_mcnp(&mat, tr->cosines, tr->value_count, false);
        } else if (mp->trcl_inline && mp->trcl_count > 0) {
            /* Inline TRCL data */
            alea_matrix_from_mcnp(&mat, mp->trcl_data, mp->trcl_count,
                                mp->trcl_degrees);
        } else {
            /* has_trcl is set but no actual transform data - skip */
            continue;
        }

        /* Allocate primitive remap cache */
        uint32_t* prim_remap = NULL;
        int8_t* prim_inverted = NULL;
        size_t remap_size = 0;

        if (alea_vec_count(&sys->primitives) > 0) {
            remap_size = alea_vec_count(&sys->primitives);
            prim_remap = calloc(remap_size, sizeof(uint32_t));
            prim_inverted = calloc(remap_size, sizeof(int8_t));
            if (prim_remap) {
                for (size_t j = 0; j < remap_size; j++) {
                    prim_remap[j] = ALEA_PRIMITIVE_ID_INVALID;
                }
            }
        }

        /* Save original tree before cloning */
        cell->original_root_node_id = cell->root_node_id;

        /* Clone tree with transform applied */
        alea_node_id_t new_root = alea_clone_tree_transformed(
            sys, cell->root_node_id, &mat, prim_remap, prim_inverted, remap_size);

        free(prim_remap);
        free(prim_inverted);

        if (new_root == ALEA_NODE_ID_INVALID) {
            ALEA_LOG_WARN("Cell %d: Failed to apply TRCL transform",
                        cell->mcnp_cell_id);
            cell->original_root_node_id = ALEA_NODE_ID_INVALID;
            continue;
        }

        /* Update cell to use transformed tree (for evaluation) */
        cell->root_node_id = new_root;

        count++;
        ALEA_LOG_DEBUG("Cell %d: Applied TRCL transform (new root node %u)",
                     cell->mcnp_cell_id, new_root);
    }

    if (count > 0) {
        ALEA_LOG_INFO("Applied TRCL transforms to %d cells", count);
    }

    return count;
}

// ============================================================================
// LIKE BUT RESOLUTION
// ============================================================================

/**
 * @brief Clone a CSG tree within the same system (no transform)
 */
static alea_node_id_t clone_tree_simple(alea_system_t* sys, alea_node_id_t root_id) {
    return alea_clone_tree_transformed(sys, root_id, NULL, NULL, NULL, 0);
}

int alea_resolve_like_cells(alea_system_t* sys, mcnp_model_t* model) {
    if (!sys || !model) return -1;

    int count = 0;

    for (size_t i = 0; i < alea_vec_count(&sys->cells); i++) {
        alea_cell_entry_t* cell = &sys->cells.data[i];
        mcnp_cell_params_t* mp = mcnp_cell_params(model, i);

        /* Skip cells that aren't LIKE cells */
        if (!mp || mp->like_cell_id == 0) continue;

        /* Already resolved? */
        if (cell->root_node_id != ALEA_NODE_ID_INVALID) continue;

        /* Find template cell */
        int template_idx = alea_find_cell_by_id(sys, mp->like_cell_id);
        if (template_idx < 0) {
            ALEA_LOG_WARN("Cell %d: LIKE references unknown cell %d",
                        cell->mcnp_cell_id, mp->like_cell_id);
            continue;
        }
        alea_cell_entry_t* template_cell = &sys->cells.data[template_idx];

        /* Template must have valid geometry */
        if (template_cell->root_node_id == ALEA_NODE_ID_INVALID) {
            ALEA_LOG_WARN("Cell %d: LIKE references cell %d with no geometry",
                        cell->mcnp_cell_id, mp->like_cell_id);
            continue;
        }

        /* Clone the template's CSG tree */
        alea_node_id_t new_root = clone_tree_simple(sys, template_cell->root_node_id);

        if (new_root == ALEA_NODE_ID_INVALID) {
            ALEA_LOG_WARN("Cell %d: Failed to clone geometry from cell %d",
                        cell->mcnp_cell_id, mp->like_cell_id);
            continue;
        }

        /* Update cell with cloned tree */
        cell->root_node_id = new_root;

        /* Inherit core parameters from template if not overridden in BUT clause */
        if (!mp->has_mat && template_cell->material_id != 0) {
            cell->material_id = template_cell->material_id;
            cell->material_index = template_cell->material_index;
        }
        if (!mp->has_rho && template_cell->density != 0.0) {
            cell->density = template_cell->density;
        }

        /* Inherit MCNP-specific parameters from template */
        const mcnp_cell_params_t* template_mp = mcnp_cell_params_const(model, (size_t)template_idx);
        if (template_mp) {
            if (!mp->has_imp_n && template_mp->has_imp_n) {
                mp->imp_n = template_mp->imp_n;
                mp->has_imp_n = 1;
            }
            if (!mp->has_imp_p && template_mp->has_imp_p) {
                mp->imp_p = template_mp->imp_p;
                mp->has_imp_p = 1;
            }
            if (!mp->has_imp_e && template_mp->has_imp_e) {
                mp->imp_e = template_mp->imp_e;
                mp->has_imp_e = 1;
            }
        }

        /* Clear LIKE reference - resolved */
        int template_id = mp->like_cell_id;  // Save for logging
        mp->like_cell_id = 0;

        count++;
        ALEA_LOG_DEBUG("Cell %d: Resolved LIKE %d (new root node %u, mat=%d, dens=%.4g)",
                     cell->mcnp_cell_id, template_id, new_root,
                     cell->material_id, cell->density);
    }

    if (count > 0) {
        ALEA_LOG_INFO("Resolved %d LIKE BUT cells", count);
    }

    return count;
}