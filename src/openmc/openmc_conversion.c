// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/**
 * @file openmc_conversion.c
 * @brief Convert OpenMC XML geometry to CSG system
 *
 * Pipeline:
 *   1. Parse XML file -> openmc_xml_doc_t tree
 *   2. Extract <geometry> section -> surfaces[], cells[]
 *   3. Extract <materials> section -> materials[]
 *   4. Create alea_system_t
 *   5. Convert surfaces (with deduplication)
 *   6. Build surface lookup table
 *   7. Convert cells (parse region expressions -> CSG trees)
 *   8. Convert materials
 *   9. Cleanup intermediate data
 */

#include "openmc_conversion.h"
#include "openmc_parse.h"
#include "openmc_region.h"
#include "alea.h"
#include "core/alea_system.h"
#include "core/alea_materials.h"
#include "core/alea_ops.h"
#include "primitives/primitive_create.h"
#include "util/alea_log.h"
#include "util/compat.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>

/* ============================================================================
 * SURFACE CONVERSION
 * ============================================================================ */

/**
 * @brief Map OpenMC surface type name to CSG primitive type
 */
static alea_primitive_type_t get_primitive_type(const char* type_str) {
    if (!type_str) return 0;

    if (strcmp(type_str, "plane") == 0) return ALEA_PRIMITIVE_PLANE;
    if (strcmp(type_str, "x-plane") == 0) return ALEA_PRIMITIVE_PLANE;
    if (strcmp(type_str, "y-plane") == 0) return ALEA_PRIMITIVE_PLANE;
    if (strcmp(type_str, "z-plane") == 0) return ALEA_PRIMITIVE_PLANE;
    if (strcmp(type_str, "sphere") == 0) return ALEA_PRIMITIVE_SPHERE;
    if (strcmp(type_str, "x-cylinder") == 0) return ALEA_PRIMITIVE_CYLINDER_X;
    if (strcmp(type_str, "y-cylinder") == 0) return ALEA_PRIMITIVE_CYLINDER_Y;
    if (strcmp(type_str, "z-cylinder") == 0) return ALEA_PRIMITIVE_CYLINDER_Z;
    if (strcmp(type_str, "x-cone") == 0) return ALEA_PRIMITIVE_CONE_X;
    if (strcmp(type_str, "y-cone") == 0) return ALEA_PRIMITIVE_CONE_Y;
    if (strcmp(type_str, "z-cone") == 0) return ALEA_PRIMITIVE_CONE_Z;
    if (strcmp(type_str, "quadric") == 0) return ALEA_PRIMITIVE_QUADRIC;
    if (strcmp(type_str, "x-torus") == 0) return ALEA_PRIMITIVE_TORUS_X;
    if (strcmp(type_str, "y-torus") == 0) return ALEA_PRIMITIVE_TORUS_Y;
    if (strcmp(type_str, "z-torus") == 0) return ALEA_PRIMITIVE_TORUS_Z;

    return 0; /* Unknown */
}

/**
 * @brief Build primitive data from OpenMC surface attributes
 * @return true on success, false if surface type/coefficients are invalid
 */
static bool build_primitive_data(const char* type_str, const double* coeffs, size_t num_coeffs,
                                  alea_primitive_type_t prim_type, alea_primitive_data_t* out_data) {
    memset(out_data, 0, sizeof(*out_data));

    switch (prim_type) {
        case ALEA_PRIMITIVE_PLANE:
            if (strcmp(type_str, "x-plane") == 0 && num_coeffs >= 1) {
                *out_data = alea_make_plane(1.0, 0.0, 0.0, -coeffs[0]);
            } else if (strcmp(type_str, "y-plane") == 0 && num_coeffs >= 1) {
                *out_data = alea_make_plane(0.0, 1.0, 0.0, -coeffs[0]);
            } else if (strcmp(type_str, "z-plane") == 0 && num_coeffs >= 1) {
                *out_data = alea_make_plane(0.0, 0.0, 1.0, -coeffs[0]);
            } else if (num_coeffs >= 4) {
                *out_data = alea_make_plane(coeffs[0], coeffs[1], coeffs[2], -coeffs[3]);
            } else {
                return false;
            }
            break;

        case ALEA_PRIMITIVE_SPHERE:
            if (num_coeffs >= 4) {
                *out_data = alea_make_sphere(coeffs[0], coeffs[1], coeffs[2], coeffs[3]);
            } else {
                return false;
            }
            break;

        case ALEA_PRIMITIVE_CYLINDER_X:
            if (num_coeffs >= 3) {
                *out_data = alea_make_cylinder_x(coeffs[0], coeffs[1], coeffs[2]);
            } else {
                return false;
            }
            break;

        case ALEA_PRIMITIVE_CYLINDER_Y:
            if (num_coeffs >= 3) {
                *out_data = alea_make_cylinder_y(coeffs[0], coeffs[1], coeffs[2]);
            } else {
                return false;
            }
            break;

        case ALEA_PRIMITIVE_CYLINDER_Z:
            if (num_coeffs >= 3) {
                *out_data = alea_make_cylinder_z(coeffs[0], coeffs[1], coeffs[2]);
            } else {
                return false;
            }
            break;

        case ALEA_PRIMITIVE_CONE_X:
            if (num_coeffs >= 4) {
                double tan_angle = sqrt(coeffs[3]);
                *out_data = alea_make_cone_x(coeffs[0], coeffs[1], coeffs[2], tan_angle);
            } else {
                return false;
            }
            break;

        case ALEA_PRIMITIVE_CONE_Y:
            if (num_coeffs >= 4) {
                double tan_angle = sqrt(coeffs[3]);
                *out_data = alea_make_cone_y(coeffs[0], coeffs[1], coeffs[2], tan_angle);
            } else {
                return false;
            }
            break;

        case ALEA_PRIMITIVE_CONE_Z:
            if (num_coeffs >= 4) {
                double tan_angle = sqrt(coeffs[3]);
                *out_data = alea_make_cone_z(coeffs[0], coeffs[1], coeffs[2], tan_angle);
            } else {
                return false;
            }
            break;

        case ALEA_PRIMITIVE_QUADRIC:
            if (num_coeffs >= 10) {
                *out_data = alea_make_quadric(coeffs);
            } else {
                return false;
            }
            break;

        case ALEA_PRIMITIVE_TORUS_X:
            if (num_coeffs >= 6) {
                *out_data = alea_make_torus_x(coeffs[0], coeffs[1], coeffs[2],
                                              coeffs[3], coeffs[4], coeffs[5]);
            } else {
                return false;
            }
            break;

        case ALEA_PRIMITIVE_TORUS_Y:
            if (num_coeffs >= 6) {
                *out_data = alea_make_torus_y(coeffs[0], coeffs[1], coeffs[2],
                                              coeffs[3], coeffs[4], coeffs[5]);
            } else {
                return false;
            }
            break;

        case ALEA_PRIMITIVE_TORUS_Z:
            if (num_coeffs >= 6) {
                *out_data = alea_make_torus_z(coeffs[0], coeffs[1], coeffs[2],
                                              coeffs[3], coeffs[4], coeffs[5]);
            } else {
                return false;
            }
            break;

        default:
            return false;
    }

    return true;
}

/**
 * @brief Convert a single OpenMC surface element to CSG primitive
 *
 * Creates BOTH positive and negative sense nodes (like MCNP conversion),
 * referencing the same primitive. Returns the positive sense node ID.
 */
static alea_node_id_t convert_surface(alea_system_t* sys,
                                      openmc_xml_element_t* surf_elem,
                                      int* out_surface_id,
                                      alea_node_id_t* out_neg_node) {
    int id = openmc_xml_get_attr_int(surf_elem, "id", -1);
    const char* type_str = openmc_xml_get_attr(surf_elem, "type");
    const char* coeffs_str = openmc_xml_get_attr(surf_elem, "coeffs");
    const char* boundary = openmc_xml_get_attr(surf_elem, "boundary");

    *out_neg_node = ALEA_NODE_ID_INVALID;

    if (id < 0 || !type_str) {
        ALEA_LOG_WARN("Surface missing id or type attribute");
        return ALEA_NODE_ID_INVALID;
    }

    *out_surface_id = id;

    /* Keep auto-surface counter above explicit IDs */
    if (id >= sys->next_auto_surface_id) {
        sys->next_auto_surface_id = id + 1;
    }

    alea_primitive_type_t prim_type = get_primitive_type(type_str);
    if (prim_type == 0) {
        ALEA_LOG_WARN("Unknown surface type: %s (surface %d)", type_str, id);
        return ALEA_NODE_ID_INVALID;
    }

    /* Parse coefficients */
    double coeffs[16] = {0};
    size_t num_coeffs = 0;
    if (coeffs_str) {
        num_coeffs = openmc_parse_doubles(coeffs_str, coeffs, 16);
    }

    /* Build primitive data */
    alea_primitive_data_t prim_data;
    if (!build_primitive_data(type_str, coeffs, num_coeffs, prim_type, &prim_data)) {
        ALEA_LOG_WARN("Failed to build primitive for surface %d (type %s, %zu coeffs)",
                     id, type_str, num_coeffs);
        return ALEA_NODE_ID_INVALID;
    }

    /* Get or create primitive (with automatic deduplication) */
    int8_t inverted = 0;
    alea_primitive_id_t prim_id = alea_get_or_create_primitive(sys, prim_type, &prim_data, &inverted);
    if (prim_id == UINT32_MAX) {
        ALEA_LOG_WARN("Failed to create primitive for surface %d", id);
        return ALEA_NODE_ID_INVALID;
    }

    /* Create POSITIVE and NEGATIVE sense nodes */
    alea_node_id_t pos_node = alea_add_primitive_node(sys, prim_id, +1, inverted, id);
    if (pos_node == ALEA_NODE_ID_INVALID) {
        return ALEA_NODE_ID_INVALID;
    }

    alea_node_id_t neg_node = alea_add_primitive_node(sys, prim_id, -1, inverted, id);
    if (neg_node == ALEA_NODE_ID_INVALID) {
        return ALEA_NODE_ID_INVALID;
    }

    /* Register surface with both node IDs */
    alea_surface_entry_t* entry = alea_vec_push_uninit(&sys->surfaces, alea_surface_entry_t);
    if (!entry) {
        return ALEA_NODE_ID_INVALID;
    }
    memset(entry, 0, sizeof(*entry));
    entry->mcnp_surface_id = id;
    entry->primitive_id = prim_id;
    entry->pos_node = pos_node;
    entry->neg_node = neg_node;
    entry->expanded_pos_node = ALEA_NODE_ID_INVALID;
    entry->expanded_neg_node = ALEA_NODE_ID_INVALID;
    entry->boundary_type = ALEA_BOUNDARY_TRANSMISSIVE;

    /* Set boundary condition */
    if (boundary) {
        if (strcmp(boundary, "vacuum") == 0) entry->boundary_type = ALEA_BOUNDARY_VACUUM;
        else if (strcmp(boundary, "reflective") == 0) entry->boundary_type = ALEA_BOUNDARY_REFLECTIVE;
        else if (strcmp(boundary, "periodic") == 0) entry->boundary_type = ALEA_BOUNDARY_PERIODIC;
        else if (strcmp(boundary, "white") == 0) entry->boundary_type = ALEA_BOUNDARY_WHITE;
    }

    *out_neg_node = neg_node;
    return pos_node;
}

/* ============================================================================
 * CELL CONVERSION
 * ============================================================================ */

/**
 * @brief Convert a single OpenMC cell element
 */
static int convert_cell(alea_system_t* sys,
                         openmc_xml_element_t* cell_elem,
                         openmc_region_ctx_t* region_ctx) {
    int id = openmc_xml_get_attr_int(cell_elem, "id", -1);
    const char* region = openmc_xml_get_attr(cell_elem, "region");
    const char* material = openmc_xml_get_attr(cell_elem, "material");
    int universe = openmc_xml_get_attr_int(cell_elem, "universe", 0);
    int fill = openmc_xml_get_attr_int(cell_elem, "fill", 0);
    const char* translation = openmc_xml_get_attr(cell_elem, "translation");

    if (id < 0) {
        ALEA_LOG_WARN("Cell missing id attribute");
        return -1;
    }

    /* Parse region expression */
    alea_node_id_t root = ALEA_NODE_ID_INVALID;
    if (region && *region) {
        root = openmc_parse_region(region_ctx, region);
        if (root == ALEA_NODE_ID_INVALID) {
            ALEA_LOG_WARN("Failed to parse region for cell %d: %s",
                         id, openmc_region_get_error(region_ctx));
            return -1;
        }
    }

    /* Parse material */
    int material_id = 0;
    double density = 0.0;
    if (material) {
        if (strcmp(material, "void") != 0) {
            material_id = (int)strtol(material, NULL, 10);
        }
    }

    /* Add cell to system */
    int cell_idx = alea_add_cell(sys, id, root, material_id, density, universe);
    if (cell_idx < 0) {
        ALEA_LOG_WARN("Failed to add cell %d", id);
        return -1;
    }

    /* Handle fill */
    if (fill > 0) {
        int fill_transform = 0;

        /* Handle translation */
        if (translation) {
            double trans[3] = {0};
            if (openmc_parse_doubles(translation, trans, 3) >= 3) {
                /* Create an inline transform */
                int tr_id = alea_add_inline_transform(sys, trans, 3, 0);
                if (tr_id > 0) {
                    fill_transform = tr_id;
                }
            }
        }

        /* Set fill */
        sys->cells.data[cell_idx].fill_universe = fill;
        sys->cells.data[cell_idx].fill_transform = fill_transform;
    }

    return cell_idx;
}

/* ============================================================================
 * MATERIAL CONVERSION
 * ============================================================================ */

/**
 * @brief Parse OpenMC nuclide name to ZAID
 *
 * Converts names like "U235", "Fe56", "H1" to ZAID format (92235, 26056, 1001)
 * Uses the materials module element database for lookup.
 *
 * @param name OpenMC nuclide name (e.g., "U235")
 * @return ZAID or 0 if parsing failed
 */
static int openmc_name_to_zaid(const char* name) {
    if (!name || !*name) return 0;

    const char* p = name;
    char symbol[3] = {0};

    /* Extract element symbol (1-2 letters) */
    if (isalpha((unsigned char)p[0])) {
        symbol[0] = p[0];
        p++;
        if (isalpha((unsigned char)p[0]) && islower((unsigned char)p[0])) {
            symbol[1] = p[0];
            p++;
        }
    }

    /* Look up element by symbol */
    const alea_element_t* elem = alea_get_element_by_symbol(symbol);
    if (!elem) return 0;

    /* Parse mass number */
    int a = 0;
    while (*p && isdigit((unsigned char)*p)) {
        a = a * 10 + (*p++ - '0');
    }

    return alea_make_zaid(elem->atomic_number, a);
}

/**
 * @brief Convert a single OpenMC material element
 */
static int convert_material(alea_system_t* sys, openmc_xml_element_t* mat_elem) {
    int id = openmc_xml_get_attr_int(mat_elem, "id", -1);
    const char* name = openmc_xml_get_attr(mat_elem, "name");

    if (id < 0) {
        ALEA_LOG_WARN("Material missing id attribute");
        return -1;
    }

    /* Add material using vector (grows automatically) */
    size_t mat_idx = alea_vec_count(&sys->materials);
    alea_material_t* mat = alea_vec_push_uninit(&sys->materials, alea_material_t);
    if (!mat) return -1;

    memset(mat, 0, sizeof(*mat));
    mat->material_id = id;
    if (name) mat->name = alea_strdup(name);

    /* Initialize nuclide array */
    mat->nuclide_capacity = 8;
    mat->nuclides = (alea_nuclide_t*)calloc(mat->nuclide_capacity, sizeof(alea_nuclide_t));
    if (!mat->nuclides) return -1;

    /* Note: Density is stored per-cell in alea_system, not per-material.
       OpenMC materials have a density element, but we skip it here.
       Cells referencing this material will need to set density separately. */

    /* Determine fraction type by checking first nuclide with wo attr */
    mat->is_weight_fraction = false;

    /* Process nuclides */
    for (size_t i = 0; i < mat_elem->child_count; i++) {
        openmc_xml_element_t* child = mat_elem->children[i];
        if (strcmp(child->tag_name, "nuclide") != 0) continue;

        const char* nuc_name = openmc_xml_get_attr(child, "name");
        double ao = openmc_xml_get_attr_double(child, "ao", 0.0);
        double wo = openmc_xml_get_attr_double(child, "wo", 0.0);

        if (wo != 0.0) {
            mat->is_weight_fraction = true;
        }

        if (nuc_name) {
            int zaid = openmc_name_to_zaid(nuc_name);
            if (zaid > 0) {
                double fraction = (wo != 0.0) ? wo : ao;
                alea_material_add_nuclide(mat, zaid, NULL, fraction);
            } else {
                ALEA_LOG_WARN("Failed to parse nuclide name: %s", nuc_name);
            }
        }
    }

    /* Process S(a,b) thermal scattering */
    size_t sab_count = openmc_xml_count_children(mat_elem, "sab");
    if (sab_count > 0) {
        mat->thermal_laws = (alea_thermal_law_t*)calloc(sab_count, sizeof(alea_thermal_law_t));
        if (mat->thermal_laws) {
            for (size_t i = 0; i < mat_elem->child_count && mat->thermal_count < sab_count; i++) {
                openmc_xml_element_t* child = mat_elem->children[i];
                if (strcmp(child->tag_name, "sab") != 0) continue;

                const char* sab_name = openmc_xml_get_attr(child, "name");
                if (sab_name) {
                    alea_thermal_law_t* law = &mat->thermal_laws[mat->thermal_count++];
                    law->identifier = alea_strdup(sab_name);
                }
            }
        }
    }

    return (int)mat_idx;
}

/* ============================================================================
 * LATTICE PARSING
 * ============================================================================ */

/**
 * @brief Temporary storage for parsed lattice data
 */
typedef struct {
    int id;
    int lat_type;           /* 1=rect, 2=hex */
    int dims[6];            /* imin, imax, jmin, jmax, kmin, kmax */
    int* universes;         /* Array of universe IDs */
    size_t universe_count;
    double pitch[3];
    double lower_left[3];
} openmc_lattice_t;

/**
 * @brief Parse a single OpenMC lattice element
 */
static openmc_lattice_t* parse_lattice(openmc_xml_element_t* lat_elem, arena_t* arena) {
    int id = openmc_xml_get_attr_int(lat_elem, "id", -1);
    if (id < 0) {
        ALEA_LOG_WARN("Lattice missing id attribute");
        return NULL;
    }

    openmc_lattice_t* lat = (openmc_lattice_t*)arena_alloc(arena, sizeof(openmc_lattice_t));
    if (!lat) return NULL;
    memset(lat, 0, sizeof(*lat));
    lat->id = id;
    lat->lat_type = 1;  /* Default to rectangular */

    /* Check for hex lattice type */
    const char* type_str = openmc_xml_get_attr(lat_elem, "type");
    if (type_str && strcmp(type_str, "hexagonal") == 0) {
        lat->lat_type = 2;
    }

    /* Parse <dimension> */
    openmc_xml_element_t* dim_elem = openmc_xml_find_child(lat_elem, "dimension");
    if (dim_elem && dim_elem->text_content) {
        int dims[3] = {1, 1, 1};
        size_t ndims = openmc_parse_ints(dim_elem->text_content, dims, 3);
        /* Convert to min:max format (0-based to MCNP-style indices) */
        lat->dims[0] = 0;                   /* imin */
        lat->dims[1] = dims[0] - 1;         /* imax */
        lat->dims[2] = 0;                   /* jmin */
        lat->dims[3] = (ndims >= 2) ? dims[1] - 1 : 0;  /* jmax */
        lat->dims[4] = 0;                   /* kmin */
        lat->dims[5] = (ndims >= 3) ? dims[2] - 1 : 0;  /* kmax */
    }

    /* Parse <pitch> */
    openmc_xml_element_t* pitch_elem = openmc_xml_find_child(lat_elem, "pitch");
    if (pitch_elem && pitch_elem->text_content) {
        openmc_parse_doubles(pitch_elem->text_content, lat->pitch, 3);
    }

    /* Parse <lower_left> */
    openmc_xml_element_t* ll_elem = openmc_xml_find_child(lat_elem, "lower_left");
    if (ll_elem && ll_elem->text_content) {
        openmc_parse_doubles(ll_elem->text_content, lat->lower_left, 3);
    }

    /* Parse <universes> */
    openmc_xml_element_t* univ_elem = openmc_xml_find_child(lat_elem, "universes");
    if (univ_elem && univ_elem->text_content) {
        /* Count universes first */
        int ni = lat->dims[1] - lat->dims[0] + 1;
        int nj = lat->dims[3] - lat->dims[2] + 1;
        int nk = lat->dims[5] - lat->dims[4] + 1;
        size_t total = (size_t)(ni * nj * nk);
        if (total > 0 && total < 1000000) {
            lat->universes = (int*)malloc(total * sizeof(int));
            if (lat->universes) {
                lat->universe_count = openmc_parse_ints(univ_elem->text_content,
                                                         lat->universes, total);
            }
        }
    }

    ALEA_LOG_DEBUG("Parsed lattice %d: type=%d, dims=[%d:%d %d:%d %d:%d], %zu universes",
                 lat->id, lat->lat_type,
                 lat->dims[0], lat->dims[1], lat->dims[2], lat->dims[3],
                 lat->dims[4], lat->dims[5], lat->universe_count);

    return lat;
}

/**
 * @brief Convert OpenMC visual row layout to internal hex grid (y-orientation)
 *
 * Replicates fill_lattice_y() traversal from OpenMC C++ source.
 * Maps each input universe (in visual row order) to the internal flat grid
 * indexed as grid[(i - imin)*nj*nk + (j - jmin)*nk + k].
 *
 * @param input     Array of universe IDs in visual row order
 * @param n_input   Number of input universe IDs
 * @param n_rings   Number of hex rings
 * @param n_axial   Number of axial layers
 * @param dims      [imin, imax, jmin, jmax, kmin, kmax]
 * @param grid      Output grid (must be pre-allocated to ni*nj*nk)
 * @param grid_size Size of the output grid
 */
static void visual_rows_to_grid(const int* input, size_t n_input,
                                 int n_rings, int n_axial,
                                 const int* dims, int* grid, size_t grid_size) {
    int nj = dims[3] - dims[2] + 1;
    int nk = dims[5] - dims[4] + 1;
    int R = n_rings;
    size_t input_index = 0;

    for (int m = 0; m < n_axial; m++) {
        int i_x = 1;
        int i_a = R - 1;

        /* Upper triangular region: R-1 rows */
        for (int k = 0; k < R - 1; k++) {
            i_x -= 1;
            for (int j = 0; j <= k; j++) {
                if (input_index < n_input) {
                    int oi = i_x - dims[0];
                    int oj = i_a - dims[2];
                    size_t idx = (size_t)(oi * nj * nk + oj * nk + m);
                    if (idx < grid_size) {
                        grid[idx] = input[input_index];
                    }
                }
                input_index++;
                i_x += 2;
                i_a -= 1;
            }
            i_x -= 2 * (k + 1);
            i_a += (k + 1);
        }

        /* Middle square region: 2R-1 rows */
        for (int k = 0; k < 2 * R - 1; k++) {
            if ((k % 2) == 0) {
                i_x -= 1;
            } else {
                i_x += 1;
                i_a -= 1;
            }
            int row_len = R - (k % 2);
            for (int j = 0; j < row_len; j++) {
                if (input_index < n_input) {
                    int oi = i_x - dims[0];
                    int oj = i_a - dims[2];
                    size_t idx = (size_t)(oi * nj * nk + oj * nk + m);
                    if (idx < grid_size) {
                        grid[idx] = input[input_index];
                    }
                }
                input_index++;
                i_x += 2;
                i_a -= 1;
            }
            i_x -= 2 * row_len;
            i_a += row_len;
        }

        /* Lower triangular region: R-1 rows */
        for (int k = 0; k < R - 1; k++) {
            i_x += 1;
            i_a -= 1;
            int row_len = R - k - 1;
            for (int j = 0; j < row_len; j++) {
                if (input_index < n_input) {
                    int oi = i_x - dims[0];
                    int oj = i_a - dims[2];
                    size_t idx = (size_t)(oi * nj * nk + oj * nk + m);
                    if (idx < grid_size) {
                        grid[idx] = input[input_index];
                    }
                }
                input_index++;
                i_x += 2;
                i_a -= 1;
            }
            i_x -= 2 * row_len;
            i_a += row_len;
        }
    }
}

/**
 * @brief Parse an OpenMC <hex_lattice> element
 */
static openmc_lattice_t* parse_hex_lattice(openmc_xml_element_t* lat_elem, arena_t* arena) {
    int id = openmc_xml_get_attr_int(lat_elem, "id", -1);
    if (id < 0) {
        ALEA_LOG_WARN("Hex lattice missing id attribute");
        return NULL;
    }

    int n_rings = openmc_xml_get_attr_int(lat_elem, "n_rings", -1);
    if (n_rings < 1) {
        ALEA_LOG_WARN("Hex lattice %d: invalid n_rings", id);
        return NULL;
    }

    int n_axial = openmc_xml_get_attr_int(lat_elem, "n_axial", 1);

    openmc_lattice_t* lat = (openmc_lattice_t*)arena_alloc(arena, sizeof(openmc_lattice_t));
    if (!lat) return NULL;
    memset(lat, 0, sizeof(*lat));
    lat->id = id;
    lat->lat_type = 2;

    /* Grid dims: symmetric around origin */
    lat->dims[0] = -(n_rings - 1);  /* imin */
    lat->dims[1] = n_rings - 1;     /* imax */
    lat->dims[2] = -(n_rings - 1);  /* jmin */
    lat->dims[3] = n_rings - 1;     /* jmax */
    lat->dims[4] = 0;               /* kmin */
    lat->dims[5] = n_axial - 1;     /* kmax */

    /* Parse <pitch> — 1 value (2D) or 2 values (3D: xy-pitch + z-pitch) */
    openmc_xml_element_t* pitch_elem = openmc_xml_find_child(lat_elem, "pitch");
    if (pitch_elem && pitch_elem->text_content) {
        double pvals[2] = {0, 0};
        size_t np = openmc_parse_doubles(pitch_elem->text_content, pvals, 2);
        lat->pitch[0] = pvals[0];
        lat->pitch[1] = pvals[0];  /* hex uses same pitch for x and y */
        lat->pitch[2] = (np >= 2) ? pvals[1] : 0;
    }

    /* Parse <center> — used as lower_left for internal storage */
    openmc_xml_element_t* center_elem = openmc_xml_find_child(lat_elem, "center");
    if (center_elem && center_elem->text_content) {
        openmc_parse_doubles(center_elem->text_content, lat->lower_left, 3);
    }

    /* Parse <universes> — visual row layout, needs conversion to grid */
    openmc_xml_element_t* univ_elem = openmc_xml_find_child(lat_elem, "universes");
    if (univ_elem && univ_elem->text_content) {
        /* Expected element count: 3*R*(R-1)+1 per axial layer */
        size_t elements_per_layer = (size_t)(3 * n_rings * (n_rings - 1) + 1);
        size_t total_visual = elements_per_layer * (size_t)n_axial;

        /* Parse visual-order input */
        int* visual_input = (int*)arena_alloc(arena, total_visual * sizeof(int));
        if (!visual_input) return NULL;
        size_t n_parsed = openmc_parse_ints(univ_elem->text_content,
                                             visual_input, total_visual);

        /* Allocate full grid (ni * nj * nk) */
        int ni = 2 * n_rings - 1;
        int nj = 2 * n_rings - 1;
        int nk = n_axial;
        size_t grid_size = (size_t)(ni * nj * nk);
        lat->universes = (int*)malloc(grid_size * sizeof(int));
        if (!lat->universes) return NULL;

        /* Initialize grid with 0 (invalid positions) */
        memset(lat->universes, 0, grid_size * sizeof(int));

        /* Convert visual rows to grid */
        visual_rows_to_grid(visual_input, n_parsed, n_rings, n_axial,
                            lat->dims, lat->universes, grid_size);
        lat->universe_count = grid_size;
    }

    /* Parse outer universe */
    int outer = openmc_xml_get_attr_int(lat_elem, "outer", -1);
    (void)outer; /* Stored at cell level if needed */

    ALEA_LOG_DEBUG("Parsed hex_lattice %d: n_rings=%d, n_axial=%d, "
                 "dims=[%d:%d %d:%d %d:%d], %zu grid entries",
                 lat->id, n_rings, n_axial,
                 lat->dims[0], lat->dims[1], lat->dims[2], lat->dims[3],
                 lat->dims[4], lat->dims[5], lat->universe_count);

    return lat;
}

/**
 * @brief Find lattice by ID in array
 */
static openmc_lattice_t* find_lattice(openmc_lattice_t** lattices, size_t count, int id) {
    for (size_t i = 0; i < count; i++) {
        if (lattices[i] && lattices[i]->id == id) {
            return lattices[i];
        }
    }
    return NULL;
}

/**
 * @brief Apply lattice data to a cell
 */
static void apply_lattice_to_cell(alea_cell_entry_t* cell, const openmc_lattice_t* lat) {
    if (!cell || !lat) return;

    cell->lat_type = lat->lat_type;
    for (int i = 0; i < 6; i++) {
        cell->lat_fill_dims[i] = lat->dims[i];
    }

    if (lat->universes && lat->universe_count > 0) {
        cell->lat_fill = (int*)malloc(lat->universe_count * sizeof(int));
        if (cell->lat_fill) {
            memcpy(cell->lat_fill, lat->universes, lat->universe_count * sizeof(int));
            cell->lat_fill_count = lat->universe_count;
        }
    }

    for (int i = 0; i < 3; i++) {
        cell->lat_pitch[i] = lat->pitch[i];
        cell->lat_lower_left[i] = lat->lower_left[i];
    }

    /* Clear fill_universe since we're using lattice fill */
    cell->fill_universe = 0;
}

/* ============================================================================
 * LATTICE ELEMENT GEOMETRY
 * ============================================================================ */

/**
 * @brief Build CSG geometry for one lattice element
 *
 * Creates the element-shape CSG tree so that the MCNP exporter can emit
 * a valid cell card for the synthetic lattice cell.
 *
 * - Rectangular (lat_type=1): RPP macrobody from pitch
 * - Hex (lat_type=2): 6 oblique planes forming a regular hexagonal prism.
 *   Flat-top orientation matching lattice_hex_lookup():
 *     3 pairs of parallel faces with normals at 0, 60, 120 deg,
 *     each at distance p/2 from center (flat-to-flat = p).
 *   Within each pair, MCNP expects the positive-sense surface first
 *   (lower D), then negative-sense (higher D).
 */
static alea_node_id_t build_lattice_element_tree(alea_system_t* sys,
                                                  const openmc_lattice_t* lat) {
    if (!sys || !lat) return ALEA_NODE_ID_INVALID;

    if (lat->lat_type == 1) {
        /* Rectangular: RPP macrobody centered at origin */
        double px = lat->pitch[0];
        double py = lat->pitch[1];
        double pz = lat->pitch[2];
        if (px <= 0) return ALEA_NODE_ID_INVALID;
        if (py <= 0) py = px;
        double zmin = (pz > 0) ? -pz / 2.0 : -1e6;
        double zmax = (pz > 0) ?  pz / 2.0 :  1e6;

        int box_idx = alea_box_surface(sys, 0,
                                        -px / 2.0, px / 2.0,
                                        -py / 2.0, py / 2.0,
                                        zmin, zmax);
        if (box_idx < 0) return ALEA_NODE_ID_INVALID;
        return alea_halfspace(sys, box_idx, -1);
    }

    if (lat->lat_type == 2) {
        /* Hex: 6 bounding planes forming a regular hexagonal prism.
         * Flat-top orientation matching lattice_hex_lookup().
         * 3 pairs of parallel faces with outward normals at 0, 60, 120 deg.
         * Each face is at distance p/2 from center (flat-to-flat = p).
         *
         * MCNP surface ordering: within each pair, the surface with
         * positive sense (lower MCNP D) comes FIRST, then negative
         * sense (higher MCNP D).  This matches the (+s -s) pattern
         * seen in standard MCNP hex lattice inputs. */
        double p = lat->pitch[0];
        if (p <= 0) return ALEA_NODE_ID_INVALID;

        double hp = p / 2.0;

        /* Normal directions for the 3 face pairs */
        static const double normals[3][2] = {
            { 1.0,  0.0 },
            { 0.5,  0.86602540378443864676 },  /* (0.5, sqrt3/2) */
            {-0.5,  0.86602540378443864676 },  /* (-0.5, sqrt3/2) */
        };

        alea_node_id_t nodes[8];
        int n = 0;

        for (int i = 0; i < 3; i++) {
            double a = normals[i][0], b = normals[i][1];
            /* "left" face:  a*x + b*y + p/2 = 0  → sense +1 gives interior
             * MCNP D = -hp (lower value) → positive sense → listed first */
            int sl = alea_plane_surface(sys, 0, a, b, 0.0,  hp);
            /* "right" face: a*x + b*y - p/2 = 0  → sense -1 gives interior
             * MCNP D = +hp (higher value) → negative sense → listed second */
            int sr = alea_plane_surface(sys, 0, a, b, 0.0, -hp);
            if (sl < 0 || sr < 0) return ALEA_NODE_ID_INVALID;
            nodes[n++] = alea_halfspace(sys, sl, +1);
            nodes[n++] = alea_halfspace(sys, sr, -1);
        }

        /* Axial bounds if 3D lattice */
        if (lat->pitch[2] > 0) {
            double pz = lat->pitch[2];
            int sb = alea_plane_surface(sys, 0, 0.0, 0.0, 1.0,  pz / 2.0);
            int st = alea_plane_surface(sys, 0, 0.0, 0.0, 1.0, -pz / 2.0);
            if (sb < 0 || st < 0) return ALEA_NODE_ID_INVALID;
            nodes[n++] = alea_halfspace(sys, sb, +1);
            nodes[n++] = alea_halfspace(sys, st, -1);
        }

        return alea_intersection_n(sys, nodes, (size_t)n);
    }

    return ALEA_NODE_ID_INVALID;
}

/* ============================================================================
 * MAIN CONVERSION
 * ============================================================================ */

static alea_system_t* convert_document(openmc_xml_doc_t* doc) {
    if (!doc || !doc->root) {
        return NULL;
    }

    /* Find geometry and materials sections */
    openmc_xml_element_t* geometry = NULL;
    openmc_xml_element_t* materials = NULL;

    /* Check if root is <model> (combined file) or <geometry> (separate file) */
    if (strcmp(doc->root->tag_name, "model") == 0) {
        geometry = openmc_xml_find_child(doc->root, "geometry");
        materials = openmc_xml_find_child(doc->root, "materials");
    } else if (strcmp(doc->root->tag_name, "geometry") == 0) {
        geometry = doc->root;
    } else if (strcmp(doc->root->tag_name, "materials") == 0) {
        materials = doc->root;
    }

    if (!geometry) {
        ALEA_LOG_ERROR("No <geometry> section found in XML");
        return NULL;
    }

    /* Create CSG system (uses ALEA_CONFIG_DEFAULT) */
    alea_system_t* sys = alea_system_create();
    if (!sys) return NULL;

    /* Create arena for temporary allocations */
    arena_t arena;
    arena_init(&arena);

    /* Initialize region context */
    openmc_region_ctx_t region_ctx;
    openmc_region_ctx_init(&region_ctx, sys, &arena);

    /* ==================== CONVERT SURFACES ==================== */
    size_t surface_count = openmc_xml_count_children(geometry, "surface");
    ALEA_LOG_INFO("Converting %zu surfaces...", surface_count);

    openmc_xml_element_t** surfaces = (openmc_xml_element_t**)arena_alloc(
        &arena, surface_count * sizeof(openmc_xml_element_t*));
    openmc_xml_find_children(geometry, "surface", surfaces, surface_count);

    for (size_t i = 0; i < surface_count; i++) {
        /* Pre-check for duplicate surface IDs before full conversion */
        int peek_id = openmc_xml_get_attr_int(surfaces[i], "id", -1);
        if (peek_id >= 0 && (size_t)peek_id < region_ctx.surface_map_size &&
            region_ctx.surface_node_map[peek_id].pos_node != ALEA_NODE_ID_INVALID) {
            ALEA_LOG_WARN("Duplicate surface id=%d, skipping", peek_id);
            continue;
        }

        int surface_id = -1;
        alea_node_id_t neg_node = ALEA_NODE_ID_INVALID;
        alea_node_id_t pos_node = convert_surface(sys, surfaces[i], &surface_id, &neg_node);
        if (pos_node != ALEA_NODE_ID_INVALID && surface_id >= 0) {
            openmc_region_register_surface(&region_ctx, surface_id, pos_node, neg_node);
        }
    }

    /* ==================== CONVERT MATERIALS ==================== */
    if (materials) {
        size_t mat_count = openmc_xml_count_children(materials, "material");
        ALEA_LOG_INFO("Converting %zu materials...", mat_count);

        openmc_xml_element_t** mat_elems = (openmc_xml_element_t**)arena_alloc(
            &arena, mat_count * sizeof(openmc_xml_element_t*));
        openmc_xml_find_children(materials, "material", mat_elems, mat_count);

        for (size_t i = 0; i < mat_count; i++) {
            convert_material(sys, mat_elems[i]);
        }
    }

    /* ==================== PARSE LATTICES ==================== */
    size_t rect_lat_count = openmc_xml_count_children(geometry, "lattice");
    size_t hex_lat_count = openmc_xml_count_children(geometry, "hex_lattice");
    size_t lattice_count = rect_lat_count + hex_lat_count;
    openmc_lattice_t** lattices = NULL;
    if (lattice_count > 0) {
        ALEA_LOG_INFO("Parsing %zu lattices (%zu rect, %zu hex)...",
                     lattice_count, rect_lat_count, hex_lat_count);

        lattices = (openmc_lattice_t**)arena_alloc(
            &arena, lattice_count * sizeof(openmc_lattice_t*));
        size_t lat_idx = 0;

        /* Parse rectangular <lattice> elements */
        if (rect_lat_count > 0) {
            openmc_xml_element_t** lat_elems = (openmc_xml_element_t**)arena_alloc(
                &arena, rect_lat_count * sizeof(openmc_xml_element_t*));
            openmc_xml_find_children(geometry, "lattice", lat_elems, rect_lat_count);
            for (size_t i = 0; i < rect_lat_count; i++) {
                lattices[lat_idx++] = parse_lattice(lat_elems[i], &arena);
            }
        }

        /* Parse hex <hex_lattice> elements */
        if (hex_lat_count > 0) {
            openmc_xml_element_t** hex_elems = (openmc_xml_element_t**)arena_alloc(
                &arena, hex_lat_count * sizeof(openmc_xml_element_t*));
            openmc_xml_find_children(geometry, "hex_lattice", hex_elems, hex_lat_count);
            for (size_t i = 0; i < hex_lat_count; i++) {
                lattices[lat_idx++] = parse_hex_lattice(hex_elems[i], &arena);
            }
        }
    }

    /* ==================== CONVERT CELLS ==================== */
    size_t cell_count = openmc_xml_count_children(geometry, "cell");
    ALEA_LOG_INFO("Converting %zu cells...", cell_count);

    openmc_xml_element_t** cells = (openmc_xml_element_t**)arena_alloc(
        &arena, cell_count * sizeof(openmc_xml_element_t*));
    openmc_xml_find_children(geometry, "cell", cells, cell_count);

    for (size_t i = 0; i < cell_count; i++) {
        convert_cell(sys, cells[i], &region_ctx);
    }

    /* ==================== CREATE GRAVEYARD CELL IF NEEDED ==================== */
    /* In OpenMC, vacuum boundaries are explicit on surfaces but no graveyard cell
     * may exist. The graveyard is the "outside" of all vacuum surfaces - i.e.,
     * the opposite sense from how cells use those surfaces. */
    {
        /* Check if any vacuum boundaries exist */
        int has_vacuum = 0;
        for (size_t i = 0; i < alea_vec_count(&sys->surfaces); i++) {
            if (sys->surfaces.data[i].boundary_type == ALEA_BOUNDARY_VACUUM) {
                has_vacuum = 1;
                break;
            }
        }

        if (has_vacuum) {
            /* Check if graveyard already exists by testing point at infinity */
            int graveyard_exists = 0;
            int cell_idx = alea_identify_cell_at_point(sys, 9.9e5, 0, 0);
            if (cell_idx >= 0) {
                alea_cell_entry_t* cell = &sys->cells.data[cell_idx];
                if (cell->material_id == 0) {
                    graveyard_exists = 1;
                    if (!cell->has_imp_n) {
                        cell->imp_n = 0.0;
                        cell->imp_p = 0.0;
                        cell->has_imp_n = 1;
                        cell->has_imp_p = 1;
                    }
                }
            }

            if (!graveyard_exists) {
                /* For each vacuum surface, determine the "outside" sense.
                 * Cells use the "inside" sense, graveyard uses the opposite.
                 * The pos_node is +surface (outside for typical bounding surfaces),
                 * the neg_node is -surface (inside). */
                alea_node_id_t graveyard_root = ALEA_NODE_ID_INVALID;
                size_t vacuum_surface_count = 0;

                for (size_t i = 0; i < alea_vec_count(&sys->surfaces); i++) {
                    alea_surface_entry_t* surf = &sys->surfaces.data[i];
                    if (surf->boundary_type == ALEA_BOUNDARY_VACUUM) {
                        /* Use pos_node = outside of surface for graveyard */
                        if (graveyard_root == ALEA_NODE_ID_INVALID) {
                            graveyard_root = surf->pos_node;
                        } else {
                            /* Union of all "outside" regions */
                            graveyard_root = alea_create_union(sys,
                                                               graveyard_root,
                                                               surf->pos_node);
                        }
                        vacuum_surface_count++;
                    }
                }

                if (graveyard_root != ALEA_NODE_ID_INVALID) {
                    /* Find next available cell ID */
                    int max_cell_id = 0;
                    for (size_t i = 0; i < alea_vec_count(&sys->cells); i++) {
                        if (sys->cells.data[i].mcnp_cell_id > max_cell_id) {
                            max_cell_id = sys->cells.data[i].mcnp_cell_id;
                        }
                    }
                    int graveyard_id = max_cell_id + 1;

                    /* Add graveyard cell */
                    alea_cell_entry_t* graveyard = alea_vec_push_uninit(&sys->cells, alea_cell_entry_t);
                    if (graveyard) {
                        memset(graveyard, 0, sizeof(*graveyard));
                        graveyard->mcnp_cell_id = graveyard_id;
                        graveyard->root_node_id = graveyard_root;
                        graveyard->material_id = 0;
                        graveyard->density = 0.0;
                        graveyard->imp_n = 0.0;
                        graveyard->imp_p = 0.0;
                        graveyard->has_imp_n = 1;
                        graveyard->has_imp_p = 1;
                        graveyard->universe_id = 0;

                        ALEA_LOG_INFO("Created graveyard cell %d from %zu vacuum surfaces",
                                     graveyard_id, vacuum_surface_count);
                    }
                }
            }
        }
    }

    /* ==================== CREATE SYNTHETIC LATTICE CELLS ==================== */
    /* OpenMC has separate <cell fill="10"> + <lattice id="10"> entities.
     * Internally we follow the MCNP convention where lattice data lives on a
     * cell whose universe_id IS the lattice ID.  Create one synthetic cell per
     * parsed lattice so the exporter can emit <lattice>/<hex_lattice> and the
     * parent cells keep their fill_universe pointing at the lattice ID. */
    if (lattices && lattice_count > 0) {
        int synth_cell_id = alea_max_cell_id(sys);
        for (size_t li = 0; li < lattice_count; li++) {
            openmc_lattice_t* lat = lattices[li];
            if (!lat) continue;

            synth_cell_id++;
            int idx = alea_add_cell_with_id(sys, synth_cell_id,
                                             ALEA_NODE_ID_INVALID,
                                             0, 0.0, lat->id);
            if (idx >= 0) {
                alea_cell_entry_t* sc = &sys->cells.data[idx];
                sc->lat_type = lat->lat_type;
                memcpy(sc->lat_fill_dims, lat->dims, sizeof(lat->dims));
                if (lat->universes && lat->universe_count > 0) {
                    sc->lat_fill = malloc(lat->universe_count * sizeof(int));
                    if (sc->lat_fill) {
                        memcpy(sc->lat_fill, lat->universes,
                               lat->universe_count * sizeof(int));
                        sc->lat_fill_count = lat->universe_count;
                    }
                }
                memcpy(sc->lat_pitch, lat->pitch, sizeof(lat->pitch));
                memcpy(sc->lat_lower_left, lat->lower_left,
                       sizeof(lat->lower_left));

                /* Build element geometry so MCNP export can emit a valid cell */
                alea_node_id_t elem_root = build_lattice_element_tree(sys, lat);
                if (elem_root != ALEA_NODE_ID_INVALID) {
                    /* Re-fetch pointer (surface creation may grow vectors) */
                    sc = &sys->cells.data[idx];
                    sc->root_node_id = elem_root;
                }

                ALEA_LOG_DEBUG("Created synthetic lattice cell %d for lattice %d (universe %d)",
                             synth_cell_id, lat->id, lat->id);
            }
        }
        /* Parent cells keep fill_universe → lattice ID (no apply_lattice_to_cell) */

        /* Free lattice universe arrays (they were copied to synthetic cells) */
        for (size_t i = 0; i < lattice_count; i++) {
            if (lattices[i] && lattices[i]->universes) {
                free(lattices[i]->universes);
            }
        }
    }

    /* Cleanup */
    arena_free(&arena);

    /* Rebuild universe index to include any newly created cells (e.g., graveyard) */
    sys->universe_index_built = false;
    alea_build_universe_index(sys);

    ALEA_LOG_INFO("Conversion complete: %zu surfaces, %zu cells, %zu materials",
                 alea_vec_count(&sys->surfaces), alea_vec_count(&sys->cells), alea_vec_count(&sys->materials));

    return sys;
}

/* ============================================================================
 * PUBLIC API
 * ============================================================================ */

alea_system_t* openmc_convert_file(const char* filename) {
    if (!filename) {
        ALEA_LOG_ERROR("NULL filename");
        return NULL;
    }

    openmc_xml_doc_t* doc = openmc_xml_parse_file(filename);
    if (!doc) {
        ALEA_LOG_ERROR("Failed to parse XML file");
        return NULL;
    }

    if (!doc->root) {
        ALEA_LOG_ERROR("XML parse error: %s", openmc_xml_get_error(doc));
        openmc_xml_doc_free(doc);
        return NULL;
    }

    alea_system_t* sys = convert_document(doc);
    openmc_xml_doc_free(doc);

    return sys;
}

alea_system_t* openmc_convert_string(const char* xml_content, size_t length) {
    if (!xml_content) return NULL;

    openmc_xml_doc_t* doc = openmc_xml_parse_string(xml_content, length);
    if (!doc || !doc->root) {
        if (doc) openmc_xml_doc_free(doc);
        return NULL;
    }

    alea_system_t* sys = convert_document(doc);
    openmc_xml_doc_free(doc);

    return sys;
}
