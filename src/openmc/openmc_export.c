// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/**
 * @file openmc_export.c
 * @brief Export CSG geometry to OpenMC XML format
 */

#include "openmc_export.h"
#include "openmc_xml.h"
#include "core/alea_system.h"
#include "core/alea_export.h"
#include "core/alea_materials.h"
#include "util/str_builder.h"
#include "util/alea_log.h"
#include "util/alea_bitset.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

/* ============================================================================
 * MATERIAL-DENSITY DENORMALIZATION
 *
 * In MCNP, density is a cell property. In OpenMC, density is a material
 * property. When the same MCNP material is used with different densities
 * in different cells, we need to create separate OpenMC materials.
 * ============================================================================ */

/**
 * @brief Entry for unique (material_id, density) pair
 */
typedef struct {
    int mc_material_id;     /* Original MCNP material ID */
    double density;           /* Cell density (always positive) */
    bool is_mass_density;     /* true = g/cm³, false = atoms/b-cm */
    int openmc_material_id;   /* Synthetic OpenMC material ID */
} mat_density_entry_t;

/**
 * @brief Material-density mapping for OpenMC export
 */
typedef struct {
    mat_density_entry_t* entries;
    size_t count;
    size_t capacity;
    int next_synthetic_id;    /* Next ID for synthetic materials */
} mat_density_map_t;

#define DENSITY_TOLERANCE 1e-10

/** Snap near-zero values to exactly 0.0 using the system tolerance */
static inline double snap_zero(double v, double threshold) {
    return fabs(v) < threshold ? 0.0 : v;
}

static void mat_density_map_init(mat_density_map_t* map, int start_id) {
    map->entries = NULL;
    map->count = 0;
    map->capacity = 0;
    map->next_synthetic_id = start_id;
}

static void mat_density_map_free(mat_density_map_t* map) {
    free(map->entries);
    map->entries = NULL;
    map->count = 0;
    map->capacity = 0;
}

/**
 * @brief Get or create OpenMC material ID for (mcnp_mat_id, density) pair
 */
static int mat_density_map_get(mat_density_map_t* map, int mcnp_mat_id, double density, bool is_mass_density) {
    /* Void cells have no material */
    if (mcnp_mat_id == 0) return 0;

    /* Search for existing entry */
    for (size_t i = 0; i < map->count; i++) {
        if (map->entries[i].mc_material_id == mcnp_mat_id &&
            fabs(map->entries[i].density - density) < DENSITY_TOLERANCE) {
            return map->entries[i].openmc_material_id;
        }
    }

    /* Create new entry */
    if (map->count >= map->capacity) {
        size_t new_cap = map->capacity == 0 ? 16 : map->capacity * 2;
        mat_density_entry_t* new_entries = realloc(map->entries,
                                                    new_cap * sizeof(mat_density_entry_t));
        if (!new_entries) return mcnp_mat_id;  /* Fallback to original ID */
        map->entries = new_entries;
        map->capacity = new_cap;
    }

    mat_density_entry_t* entry = &map->entries[map->count++];
    entry->mc_material_id = mcnp_mat_id;
    entry->density = density;
    entry->is_mass_density = is_mass_density;

    /* Check if this is the first use of this material */
    int uses_count = 0;
    for (size_t i = 0; i < map->count - 1; i++) {
        if (map->entries[i].mc_material_id == mcnp_mat_id) {
            uses_count++;
            break;
        }
    }

    if (uses_count == 0) {
        /* First use - keep original ID */
        entry->openmc_material_id = mcnp_mat_id;
    } else {
        /* Additional use with different density - assign synthetic ID */
        entry->openmc_material_id = map->next_synthetic_id++;
    }

    return entry->openmc_material_id;
}

/**
 * @brief Build material-density map from all cells
 */
static void build_mat_density_map(const alea_system_t* sys, mat_density_map_t* map) {
    /* Find max material ID for synthetic ID starting point */
    int max_mat_id = 0;
    for (size_t i = 0; i < alea_vec_count(&sys->materials); i++) {
        if (sys->materials.data[i].material_id > max_mat_id) {
            max_mat_id = sys->materials.data[i].material_id;
        }
    }

    /* Start synthetic IDs after max material ID */
    mat_density_map_init(map, max_mat_id + 1000);

    /* Scan all cells to populate the map */
    for (size_t i = 0; i < alea_vec_count(&sys->cells); i++) {
        const alea_cell_entry_t* cell = &sys->cells.data[i];
        if (cell->material_id > 0) {
            mat_density_map_get(map, cell->material_id, cell->density, cell->is_mass_density);
        }
    }
}

/* ============================================================================
 * OPENMC SURFACE TYPE MAPPING
 * ============================================================================ */

#define PLANE_AXIS_TOL 1e-10

/**
 * @brief Plane alignment type
 */
typedef enum {
    PLANE_GENERAL,
    PLANE_X,
    PLANE_Y,
    PLANE_Z
} plane_alignment_t;

/**
 * @brief Detect if plane is axis-aligned
 *
 * A plane ax + by + cz + d = 0 is axis-aligned if only one normal component is non-zero.
 * We require the non-zero coefficient to be exactly 1.0 (normalized form from MCNP parsing).
 */
static plane_alignment_t detect_plane_alignment(const alea_plane_data_t* p) {
    int is_x = fabs(p->a - 1.0) < PLANE_AXIS_TOL && fabs(p->b) < PLANE_AXIS_TOL && fabs(p->c) < PLANE_AXIS_TOL;
    int is_y = fabs(p->a) < PLANE_AXIS_TOL && fabs(p->b - 1.0) < PLANE_AXIS_TOL && fabs(p->c) < PLANE_AXIS_TOL;
    int is_z = fabs(p->a) < PLANE_AXIS_TOL && fabs(p->b) < PLANE_AXIS_TOL && fabs(p->c - 1.0) < PLANE_AXIS_TOL;

    if (is_x) return PLANE_X;
    if (is_y) return PLANE_Y;
    if (is_z) return PLANE_Z;
    return PLANE_GENERAL;
}

/**
 * @brief Get OpenMC surface type name for a primitive
 *
 * For planes, detects axis-alignment and returns x-plane/y-plane/z-plane accordingly.
 */
static const char* get_openmc_surface_type(const alea_primitive_entry_t* prim, int8_t inverted) {
    switch (prim->type) {
        case ALEA_PRIMITIVE_PLANE: {
            if (inverted) {
                /* Un-canonicalized plane has flipped leading coeff — no longer axis-aligned */
                return "plane";
            }
            plane_alignment_t align = detect_plane_alignment(&prim->data.plane);
            switch (align) {
                case PLANE_X: return "x-plane";
                case PLANE_Y: return "y-plane";
                case PLANE_Z: return "z-plane";
                default:      return "plane";
            }
        }
        case ALEA_PRIMITIVE_SPHERE:     return "sphere";
        case ALEA_PRIMITIVE_CYLINDER_X: return "x-cylinder";
        case ALEA_PRIMITIVE_CYLINDER_Y: return "y-cylinder";
        case ALEA_PRIMITIVE_CYLINDER_Z: return "z-cylinder";
        case ALEA_PRIMITIVE_CONE_X:     return "x-cone";
        case ALEA_PRIMITIVE_CONE_Y:     return "y-cone";
        case ALEA_PRIMITIVE_CONE_Z:     return "z-cone";
        case ALEA_PRIMITIVE_QUADRIC:    return "quadric";
        case ALEA_PRIMITIVE_TORUS_X:    return "x-torus";
        case ALEA_PRIMITIVE_TORUS_Y:    return "y-torus";
        case ALEA_PRIMITIVE_TORUS_Z:    return "z-torus";
        default:                       return NULL; /* Macrobodies need decomposition */
    }
}

/**
 * @brief Check if a primitive type is a macrobody (expanded into constituent surfaces)
 */
static bool is_macrobody_type(alea_primitive_type_t type) {
    switch (type) {
        case ALEA_PRIMITIVE_RPP:
        case ALEA_PRIMITIVE_RCC:
        case ALEA_PRIMITIVE_BOX:
        case ALEA_PRIMITIVE_SPH:
        case ALEA_PRIMITIVE_TRC:
        case ALEA_PRIMITIVE_ELL:
        case ALEA_PRIMITIVE_REC:
        case ALEA_PRIMITIVE_WED:
        case ALEA_PRIMITIVE_RHP:
        case ALEA_PRIMITIVE_ARB:
            return true;
        default:
            return false;
    }
}

static const char* macrobody_type_name(alea_primitive_type_t type) {
    switch (type) {
        case ALEA_PRIMITIVE_RPP:         return "RPP";
        case ALEA_PRIMITIVE_RCC:         return "RCC";
        case ALEA_PRIMITIVE_BOX:         return "BOX";
        case ALEA_PRIMITIVE_SPH:         return "SPH";
        case ALEA_PRIMITIVE_TRC:         return "TRC";
        case ALEA_PRIMITIVE_ELL:         return "ELL";
        case ALEA_PRIMITIVE_REC:         return "REC";
        case ALEA_PRIMITIVE_WED:         return "WED";
        case ALEA_PRIMITIVE_RHP:         return "RHP";
        case ALEA_PRIMITIVE_ARB:         return "ARB";
        default:                         return "unknown";
    }
}

/**
 * @brief Write surface coefficients to XML builder
 *
 * Returns the coefficients as a space-separated string in the provided buffer.
 * For axis-aligned planes (a=1, b=0, c=0 form), outputs just the intercept value.
 *
 * Internal: ax + by + cz + d = 0
 * OpenMC x-plane: x = x0  (for a=1, d=-x0, so x0 = -d)
 * OpenMC plane: Ax + By + Cz = D  (D = -d in our convention)
 */
static int format_surface_coeffs(const alea_primitive_entry_t* prim, char* buf, size_t buf_size,
                                 int8_t inverted) {
    const alea_primitive_data_t* d = &prim->data;

    switch (prim->type) {
        case ALEA_PRIMITIVE_PLANE: {
            /* Un-canonicalize: restore original sign when inverted */
            double a = d->plane.a, b = d->plane.b, c = d->plane.c, pd = d->plane.d;
            if (inverted) { a = -a; b = -b; c = -c; pd = -pd; }

            /* Check for axis-aligned plane */
            alea_plane_data_t tmp_plane = { a, b, c, pd };
            plane_alignment_t align = detect_plane_alignment(&tmp_plane);
            switch (align) {
                case PLANE_X:
                    return snprintf(buf, buf_size, "%.16g", -pd);
                case PLANE_Y:
                    return snprintf(buf, buf_size, "%.16g", -pd);
                case PLANE_Z:
                    return snprintf(buf, buf_size, "%.16g", -pd);
                default:
                    /* OpenMC general plane: Ax + By + Cz = D, D = -d */
                    return snprintf(buf, buf_size, "%.16g %.16g %.16g %.16g",
                                   a, b, c, -pd);
            }
        }

        case ALEA_PRIMITIVE_SPHERE:
            /* OpenMC sphere: x0 y0 z0 R */
            return snprintf(buf, buf_size, "%.16g %.16g %.16g %.16g",
                           d->sphere.center_x, d->sphere.center_y, d->sphere.center_z,
                           d->sphere.radius);

        case ALEA_PRIMITIVE_CYLINDER_X:
            /* OpenMC x-cylinder: y0 z0 R */
            return snprintf(buf, buf_size, "%.16g %.16g %.16g",
                           d->cyl_x.center_y, d->cyl_x.center_z, d->cyl_x.radius);

        case ALEA_PRIMITIVE_CYLINDER_Y:
            /* OpenMC y-cylinder: x0 z0 R */
            return snprintf(buf, buf_size, "%.16g %.16g %.16g",
                           d->cyl_y.center_x, d->cyl_y.center_z, d->cyl_y.radius);

        case ALEA_PRIMITIVE_CYLINDER_Z:
            /* OpenMC z-cylinder: x0 y0 R */
            return snprintf(buf, buf_size, "%.16g %.16g %.16g",
                           d->cyl_z.center_x, d->cyl_z.center_y, d->cyl_z.radius);

        case ALEA_PRIMITIVE_CONE_X:
            /* OpenMC x-cone: x0 y0 z0 R^2 */
            return snprintf(buf, buf_size, "%.16g %.16g %.16g %.16g",
                           d->cone_x.apex_x, d->cone_x.apex_y, d->cone_x.apex_z,
                           d->cone_x.tan_angle_sq);

        case ALEA_PRIMITIVE_CONE_Y:
            /* OpenMC y-cone: x0 y0 z0 R^2 */
            return snprintf(buf, buf_size, "%.16g %.16g %.16g %.16g",
                           d->cone_y.apex_x, d->cone_y.apex_y, d->cone_y.apex_z,
                           d->cone_y.tan_angle_sq);

        case ALEA_PRIMITIVE_CONE_Z:
            /* OpenMC z-cone: x0 y0 z0 R^2 */
            return snprintf(buf, buf_size, "%.16g %.16g %.16g %.16g",
                           d->cone_z.apex_x, d->cone_z.apex_y, d->cone_z.apex_z,
                           d->cone_z.tan_angle_sq);

        case ALEA_PRIMITIVE_QUADRIC:
            /* OpenMC quadric: A B C D E F G H J K (10 coefficients) */
            return snprintf(buf, buf_size,
                           "%.16g %.16g %.16g %.16g %.16g %.16g %.16g %.16g %.16g %.16g",
                           d->quadric.coeffs[0], d->quadric.coeffs[1], d->quadric.coeffs[2],
                           d->quadric.coeffs[3], d->quadric.coeffs[4], d->quadric.coeffs[5],
                           d->quadric.coeffs[6], d->quadric.coeffs[7], d->quadric.coeffs[8],
                           d->quadric.coeffs[9]);

        case ALEA_PRIMITIVE_TORUS_X:
        case ALEA_PRIMITIVE_TORUS_Y:
        case ALEA_PRIMITIVE_TORUS_Z:
            /* OpenMC torus: x0 y0 z0 A B C
               A = major radius, B = C = minor radius (circular cross-section) */
            return snprintf(buf, buf_size, "%.16g %.16g %.16g %.16g %.16g %.16g",
                           d->torus.center_x, d->torus.center_y, d->torus.center_z,
                           d->torus.major_radius, d->torus.minor_radius,
                           d->torus.axial_semiwidth_B > 0 ? d->torus.axial_semiwidth_B : d->torus.minor_radius);

        default:
            return -1;
    }
}

/* ============================================================================
 * MATERIAL UTILITIES
 * ============================================================================ */

/**
 * @brief Convert ZAID to OpenMC nuclide name
 *
 * ZAID format: ZZAAA or ZZAAAI (with isomer)
 * OpenMC format: ElementMassNumber (e.g., "U235", "H1", "Pu239")
 *
 * Uses the materials module element database for symbol lookup.
 *
 * @param zaid ZAID number (e.g., 92235 for U-235)
 * @param buf Output buffer
 * @param buf_size Buffer size
 * @return Number of characters written, or -1 on error
 */
static int zaid_to_openmc_name(int zaid, char* buf, size_t buf_size) {
    if (!buf || buf_size < 8) return -1;

    int z = alea_zaid_to_Z(zaid);
    int a = alea_zaid_to_A(zaid);

    const alea_element_t* elem = alea_get_element(z);
    if (!elem) {
        return snprintf(buf, buf_size, "%d", zaid);  /* Fallback: raw number */
    }

    /* Natural element (A=0) uses just the symbol */
    if (a == 0) {
        return snprintf(buf, buf_size, "%s", elem->symbol);
    }

    return snprintf(buf, buf_size, "%s%d", elem->symbol, a);
}

/**
 * @brief Normalize material fractions so they sum to 1.0
 *
 * OpenMC requires fractions to sum to 1.0. This function normalizes
 * in-place without modifying the original.
 *
 * @param fractions Array of fractions
 * @param count Number of fractions
 * @param normalized Output array (must have 'count' elements)
 * @return Sum before normalization (useful for detecting bad input)
 */
static double normalize_fractions(const double* fractions, size_t count, double* normalized) {
    double sum = 0.0;
    for (size_t i = 0; i < count; i++) {
        sum += fractions[i];
    }

    if (sum > 0.0) {
        for (size_t i = 0; i < count; i++) {
            normalized[i] = fractions[i] / sum;
        }
    } else {
        /* Avoid division by zero - uniform distribution */
        double uniform = 1.0 / (double)count;
        for (size_t i = 0; i < count; i++) {
            normalized[i] = uniform;
        }
    }

    return sum;
}

/* ============================================================================
 * CSG TREE TO OPENMC REGION EXPRESSION
 * ============================================================================ */

/* Expression context for precedence handling */
typedef enum {
    OPENMC_CTX_TOP,
    OPENMC_CTX_UNION,
    OPENMC_CTX_INTERSECTION
} openmc_expr_context_t;

/**
 * @brief Recursively convert CSG tree to OpenMC region expression
 *
 * Uses str_builder_t for safe dynamic string building.
 */
static bool tree_to_openmc_region(const alea_system_t* sys,
                                   export_context_t* ctx,
                                   uint32_t node_id,
                                   str_builder_t* sb,
                                   openmc_expr_context_t parent_ctx) {
    if (node_id == ALEA_NODE_ID_INVALID || node_id >= alea_vec_count(&sys->nodes)) {
        return false;
    }

    const alea_node_t* node = &sys->nodes.data[node_id];
    alea_operation_t op = ALEA_GET_OPERATION(node);

    switch (op) {
        case ALEA_OP_PRIMITIVE: {
            /* In OpenMC: -N means inside (negative half-space), N means outside */
            int surf_id = node->primitive.mc_surface_id;

            /* Compute effective sense relative to the EXPORTED surface coefficients.
               Surface export un-canonicalizes plane coefficients using the surface's
               pos_node inverted flag. Cell sense must be relative to that orientation.
               Formula: flip sense when node's inverted differs from the exported
               surface's inverted (XOR). In non-dedup mode, they always match
               (same surface), so effective_sense = sense. */
            int effective_sense = node->primitive.sense;

            /* Use canonical surface ID when deduplication is enabled */
            if (ctx && ctx->deduplicate && ctx->prim_to_surface) {
                uint32_t prim_id = node->primitive.primitive_id;
                if (prim_id < ctx->prim_to_surface_size &&
                    ctx->prim_to_surface[prim_id] >= 0) {
                    surf_id = ctx->prim_to_surface[prim_id];
                }
                /* XOR: flip sense when node's inverted differs from canonical surface's inverted */
                int8_t surface_inverted = 0;
                if (prim_id < ctx->prim_to_surface_size && ctx->prim_to_surface_inverted) {
                    surface_inverted = ctx->prim_to_surface_inverted[prim_id];
                }
                if (node->primitive.inverted != surface_inverted) {
                    effective_sense = -effective_sense;
                }
            }

            if (effective_sense < 0) {
                str_builder_putc(sb, '-');
            }
            str_builder_int(sb, surf_id);
            return !str_builder_error(sb);
        }

        case ALEA_OP_COMPLEMENT: {
            /* OpenMC complement: ~(...) */
            uint32_t left = node->operation.left;
            if (left >= alea_vec_count(&sys->nodes)) return false;

            const alea_node_t* child = &sys->nodes.data[left];
            alea_operation_t child_op = ALEA_GET_OPERATION(child);

            /* Optimization: complement of single primitive just flips sign */
            if (child_op == ALEA_OP_PRIMITIVE) {
                int surf_id = child->primitive.mc_surface_id;
                int effective_sense = child->primitive.sense;

                /* Use canonical surface ID when deduplication is enabled */
                if (ctx && ctx->deduplicate && ctx->prim_to_surface) {
                    uint32_t prim_id = child->primitive.primitive_id;
                    if (prim_id < ctx->prim_to_surface_size &&
                        ctx->prim_to_surface[prim_id] >= 0) {
                        surf_id = ctx->prim_to_surface[prim_id];
                    }
                    /* XOR: flip sense when node's inverted differs from canonical surface's inverted */
                    int8_t surface_inverted = 0;
                    if (prim_id < ctx->prim_to_surface_size && ctx->prim_to_surface_inverted) {
                        surface_inverted = ctx->prim_to_surface_inverted[prim_id];
                    }
                    if (child->primitive.inverted != surface_inverted) {
                        effective_sense = -effective_sense;
                    }
                }

                /* Complement flips the sense */
                if (effective_sense > 0) {
                    /* Was outside, complement is inside */
                    str_builder_putc(sb, '-');
                }
                str_builder_int(sb, surf_id);
                return !str_builder_error(sb);
            }

            /* Complex expression: use ~(...) */
            str_builder_putc(sb, '~');
            str_builder_putc(sb, '(');
            if (!tree_to_openmc_region(sys, ctx, left, sb, OPENMC_CTX_TOP)) {
                return false;
            }
            str_builder_putc(sb, ')');
            return !str_builder_error(sb);
        }

        case ALEA_OP_UNION: {
            /* OpenMC union: (A | B) */
            uint32_t left = node->operation.left;
            uint32_t right = node->operation.right;
            if (left >= alea_vec_count(&sys->nodes) || right >= alea_vec_count(&sys->nodes)) return false;

            bool need_parens = (parent_ctx == OPENMC_CTX_INTERSECTION);

            if (need_parens) str_builder_putc(sb, '(');

            if (!tree_to_openmc_region(sys, ctx, left, sb, OPENMC_CTX_UNION)) {
                return false;
            }

            str_builder_puts(sb, " | ");

            if (!tree_to_openmc_region(sys, ctx, right, sb, OPENMC_CTX_UNION)) {
                return false;
            }

            if (need_parens) str_builder_putc(sb, ')');
            return !str_builder_error(sb);
        }

        case ALEA_OP_INTERSECTION: {
            /* OpenMC intersection: A B (space-separated) */
            uint32_t left = node->operation.left;
            uint32_t right = node->operation.right;
            if (left >= alea_vec_count(&sys->nodes) || right >= alea_vec_count(&sys->nodes)) return false;

            /* Parenthesize intersection branches inside a union so OpenMC
               doesn't misparse "A B | C D" — it needs "(A B) | (C D)" */
            bool need_parens = (parent_ctx == OPENMC_CTX_UNION);

            if (need_parens) str_builder_putc(sb, '(');

            if (!tree_to_openmc_region(sys, ctx, left, sb, OPENMC_CTX_INTERSECTION)) {
                return false;
            }

            str_builder_putc(sb, ' ');

            if (!tree_to_openmc_region(sys, ctx, right, sb, OPENMC_CTX_INTERSECTION)) {
                return false;
            }

            if (need_parens) str_builder_putc(sb, ')');
            return !str_builder_error(sb);
        }

        case ALEA_OP_DIFFERENCE: {
            /* A - B = A & ~B in OpenMC */
            uint32_t left = node->operation.left;
            uint32_t right = node->operation.right;
            if (left >= alea_vec_count(&sys->nodes) || right >= alea_vec_count(&sys->nodes)) return false;

            if (!tree_to_openmc_region(sys, ctx, left, sb, OPENMC_CTX_INTERSECTION)) {
                return false;
            }

            str_builder_puts(sb, " ~(");

            if (!tree_to_openmc_region(sys, ctx, right, sb, OPENMC_CTX_TOP)) {
                return false;
            }

            str_builder_putc(sb, ')');
            return !str_builder_error(sb);
        }

        default:
            return false;
    }
}

/* ============================================================================
 * OPENMC MODEL.XML EXPORT (Combined geometry, materials, settings)
 * ============================================================================ */

/**
 * @brief Find material by MCNP ID
 */
static const alea_material_t* find_material_by_id(const alea_system_t* sys, int material_id) {
    for (size_t i = 0; i < alea_vec_count(&sys->materials); i++) {
        if (sys->materials.data[i].material_id == material_id) {
            return &sys->materials.data[i];
        }
    }
    return NULL;
}

/**
 * @brief Write a single material entry to XML
 *
 * Used for each unique (material_id, density) pair.
 */
static const alea_mixture_t* find_mixture_by_id(const alea_system_t* sys, int mixture_id) {
    for (size_t i = 0; i < alea_vec_count(&sys->mixtures); i++) {
        if (sys->mixtures.data[i].mixture_id == mixture_id) {
            return &sys->mixtures.data[i];
        }
    }
    return NULL;
}

/**
 * @brief Write a mixture as a material entry to XML
 *
 * Combines nuclides from all component materials, scaling fractions.
 */
static bool write_mixture_entry(const alea_system_t* sys, openmc_xml_t* xml,
                                arena_t* arena, const mat_density_entry_t* entry,
                                const alea_mixture_t* mix) {
    /* <material id="N"> */
    if (!openmc_xml_start_element(xml, "material")) return false;
    if (!openmc_xml_attribute_i(xml, "id", entry->openmc_material_id)) return false;

    char mat_name[128];
    if (mix->name && mix->name[0]) {
        snprintf(mat_name, sizeof(mat_name), "%s", mix->name);
    } else {
        snprintf(mat_name, sizeof(mat_name), "Mixture_M%d", mix->mixture_id);
    }
    if (!openmc_xml_attribute(xml, "name", mat_name)) return false;
    if (!openmc_xml_end_start_tag(xml, false)) return false;

    /* <density> */
    if (entry->density != 0.0) {
        if (!openmc_xml_start_element(xml, "density")) return false;
        if (!openmc_xml_attribute_f(xml, "value", entry->density, 10, false)) return false;
        if (entry->is_mass_density) {
            if (!openmc_xml_attribute(xml, "units", "g/cc")) return false;
        } else {
            if (!openmc_xml_attribute(xml, "units", "atom/b-cm")) return false;
        }
        if (!openmc_xml_end_start_tag(xml, true)) return false;
    }

    /* Count total nuclides across all components */
    size_t total_nuclides = 0;
    for (size_t c = 0; c < alea_vec_count(&mix->components); c++) {
        const alea_material_t* comp_mat = find_material_by_id(sys, mix->components.data[c].material_id);
        if (comp_mat) total_nuclides += alea_vec_count(&comp_mat->nuclides);
    }
    if (total_nuclides == 0) {
        if (!openmc_xml_end_element(xml, "material")) return false;
        return true;
    }

    /* Collect scaled fractions */
    double* fractions = (double*)arena_alloc(arena, total_nuclides * sizeof(double));
    if (!fractions) return false;

    size_t idx = 0;
    bool is_weight = false;
    for (size_t c = 0; c < alea_vec_count(&mix->components); c++) {
        const alea_mixture_comp_t* comp = &mix->components.data[c];
        const alea_material_t* comp_mat = find_material_by_id(sys, comp->material_id);
        if (!comp_mat) continue;
        if (comp_mat->is_weight_fraction) is_weight = true;
        for (size_t j = 0; j < alea_vec_count(&comp_mat->nuclides); j++) {
            fractions[idx++] = comp_mat->nuclides.data[j].fraction * comp->fraction;
        }
    }

    /* Normalize */
    double* normalized = (double*)arena_alloc(arena, total_nuclides * sizeof(double));
    if (!normalized) return false;
    normalize_fractions(fractions, total_nuclides, normalized);

    /* Write nuclides */
    const char* fraction_attr = is_weight ? "wo" : "ao";
    idx = 0;
    for (size_t c = 0; c < alea_vec_count(&mix->components); c++) {
        const alea_mixture_comp_t* comp = &mix->components.data[c];
        const alea_material_t* comp_mat = find_material_by_id(sys, comp->material_id);
        if (!comp_mat) continue;
        for (size_t j = 0; j < alea_vec_count(&comp_mat->nuclides); j++) {
            const alea_nuclide_t* nuc = &comp_mat->nuclides.data[j];
            char nuclide_name[32];
            if (zaid_to_openmc_name(nuc->zaid, nuclide_name, sizeof(nuclide_name)) < 0) {
                ALEA_LOG_WARN("Failed to convert ZAID %d", nuc->zaid);
                idx++;
                continue;
            }
            if (!openmc_xml_start_element(xml, "nuclide")) return false;
            if (!openmc_xml_attribute(xml, "name", nuclide_name)) return false;
            if (!openmc_xml_attribute_f(xml, fraction_attr, normalized[idx], 10, false)) return false;
            if (!openmc_xml_end_start_tag(xml, true)) return false;
            idx++;
        }
    }

    /* Thermal scattering laws from all components */
    for (size_t c = 0; c < alea_vec_count(&mix->components); c++) {
        const alea_material_t* comp_mat = find_material_by_id(sys, mix->components.data[c].material_id);
        if (!comp_mat) continue;
        for (size_t t = 0; t < alea_vec_count(&comp_mat->thermal_laws); t++) {
            const alea_thermal_law_t* law = &comp_mat->thermal_laws.data[t];
            if (law->identifier && law->identifier[0]) {
                if (!openmc_xml_start_element(xml, "sab")) return false;
                if (!openmc_xml_attribute(xml, "name", law->identifier)) return false;
                if (!openmc_xml_end_start_tag(xml, true)) return false;
            }
        }
    }

    if (!openmc_xml_end_element(xml, "material")) return false;
    return true;
}

static bool write_material_entry(const alea_system_t* sys, openmc_xml_t* xml,
                                  arena_t* arena, const mat_density_entry_t* entry) {
    const alea_material_t* mat = find_material_by_id(sys, entry->mc_material_id);
    if (!mat) {
        /* Check if it's a mixture */
        const alea_mixture_t* mix = find_mixture_by_id(sys, entry->mc_material_id);
        if (mix) return write_mixture_entry(sys, xml, arena, entry, mix);
        return true;  /* Unknown material, skip */
    }
    if (alea_vec_count(&mat->nuclides) == 0) return true;  /* Skip empty materials */

    /* <material id="N"> */
    if (!openmc_xml_start_element(xml, "material")) return false;
    if (!openmc_xml_attribute_i(xml, "id", entry->openmc_material_id)) return false;

    /* Name: include original material ID and density for synthetic materials */
    char mat_name[128];
    if (entry->openmc_material_id != entry->mc_material_id) {
        /* Synthetic material - note the original */
        snprintf(mat_name, sizeof(mat_name), "M%d_rho%.4g",
                 entry->mc_material_id, entry->density);
    } else if (mat->name && mat->name[0]) {
        snprintf(mat_name, sizeof(mat_name), "%s", mat->name);
    } else {
        snprintf(mat_name, sizeof(mat_name), "M%d", entry->mc_material_id);
    }
    if (!openmc_xml_attribute(xml, "name", mat_name)) return false;

    if (!openmc_xml_end_start_tag(xml, false)) return false;

    /* <density value="X" units="g/cc" /> */
    if (entry->density != 0.0) {
        if (!openmc_xml_start_element(xml, "density")) return false;
        if (!openmc_xml_attribute_f(xml, "value", entry->density, 10, false)) return false;
        if (entry->is_mass_density) {
            if (!openmc_xml_attribute(xml, "units", "g/cc")) return false;
        } else {
            if (!openmc_xml_attribute(xml, "units", "atom/b-cm")) return false;
        }
        if (!openmc_xml_end_start_tag(xml, true)) return false;
    }

    /* Normalize fractions */
    size_t nuc_count = alea_vec_count(&mat->nuclides);
    double* normalized = (double*)arena_alloc(arena, nuc_count * sizeof(double));
    if (!normalized) return false;

    double* orig_fractions = (double*)arena_alloc(arena, nuc_count * sizeof(double));
    if (!orig_fractions) return false;

    for (size_t i = 0; i < nuc_count; i++) {
        orig_fractions[i] = mat->nuclides.data[i].fraction;
    }

    double orig_sum = normalize_fractions(orig_fractions, nuc_count, normalized);

    /* Warn if fractions were significantly off from 1.0 */
    if (fabs(orig_sum - 1.0) > 0.01) {
        ALEA_LOG_DEBUG("Note: Material %d fractions summed to %.6f, normalized to 1.0",
                     mat->material_id, orig_sum);
    }

    /* Export nuclides */
    const char* fraction_attr = mat->is_weight_fraction ? "wo" : "ao";

    for (size_t i = 0; i < nuc_count; i++) {
        const alea_nuclide_t* nuc = &mat->nuclides.data[i];

        char nuclide_name[32];
        if (zaid_to_openmc_name(nuc->zaid, nuclide_name, sizeof(nuclide_name)) < 0) {
            ALEA_LOG_WARN("Failed to convert ZAID %d", nuc->zaid);
            continue;
        }

        /* <nuclide name="U235" wo="0.05" /> */
        if (!openmc_xml_start_element(xml, "nuclide")) return false;
        if (!openmc_xml_attribute(xml, "name", nuclide_name)) return false;
        if (!openmc_xml_attribute_f(xml, fraction_attr, normalized[i], 10, false)) return false;
        if (!openmc_xml_end_start_tag(xml, true)) return false;
    }

    /* Thermal scattering laws (S(a,b) tables) */
    for (size_t t = 0; t < alea_vec_count(&mat->thermal_laws); t++) {
        const alea_thermal_law_t* law = &mat->thermal_laws.data[t];
        if (law->identifier && law->identifier[0]) {
            /* <sab name="c_H_in_H2O" /> */
            if (!openmc_xml_start_element(xml, "sab")) return false;
            if (!openmc_xml_attribute(xml, "name", law->identifier)) return false;
            if (!openmc_xml_end_start_tag(xml, true)) return false;
        }
    }

    /* </material> */
    if (!openmc_xml_end_element(xml, "material")) return false;

    return true;
}

/**
 * @brief Write materials section to XML builder
 *
 * Uses the material-density map to create separate materials for each
 * unique (material_id, density) pair. This handles the MCNP->OpenMC
 * model difference where MCNP has cell-level density but OpenMC has
 * material-level density.
 */
static bool write_materials_section(const alea_system_t* sys, openmc_xml_t* xml,
                                     arena_t* arena, const mat_density_map_t* map) {
    if (map->count == 0) return true;  /* No materials to export */

    /* <materials> */
    if (!openmc_xml_start_element(xml, "materials")) return false;
    if (!openmc_xml_end_start_tag(xml, false)) return false;

    /* Export each unique (material, density) pair */
    for (size_t i = 0; i < map->count; i++) {
        if (!write_material_entry(sys, xml, arena, &map->entries[i])) {
            return false;
        }
    }

    /* </materials> */
    if (!openmc_xml_end_element(xml, "materials")) return false;

    return true;
}

/* ============================================================================
 * HEX LATTICE HELPERS
 * ============================================================================ */

/**
 * @brief Compute n_rings from hex lattice fill_dims
 *
 * For a hex lattice centered at the origin, n_rings = max cube-distance + 1
 * over all positions in the axial bounding rect.
 */
static int compute_n_rings(const int* fill_dims) {
    int max_ring = 0;
    for (int i = fill_dims[0]; i <= fill_dims[1]; i++) {
        for (int j = fill_dims[2]; j <= fill_dims[3]; j++) {
            int k_cube = -i - j;
            int ring = (abs(i) + abs(j) + abs(k_cube)) / 2;
            if (ring > max_ring) max_ring = ring;
        }
    }
    return max_ring + 1;
}

/**
 * @brief Convert internal hex grid to OpenMC visual row order (y-orientation)
 *
 * Replicates fill_lattice_y() traversal from OpenMC, but reads from our
 * internal grid instead of writing. Our axial coords (i,j) map directly
 * to OpenMC's (i_x, i_a).
 *
 * @param lat_fill     Internal universe array indexed as [oi*nj*nk + oj*nk + ok]
 * @param fill_dims    [imin, imax, jmin, jmax, kmin, kmax]
 * @param n_rings      Number of hex rings
 * @param map_func     Callback to map original universe ID through wrapper table (or NULL)
 * @param map_ctx      Opaque context for map_func
 * @param sb           String builder to append universe IDs to
 */
typedef int (*univ_map_func_t)(int uid, void* ctx);

static void hex_grid_to_visual_rows(const int* lat_fill, const int* fill_dims,
                                     size_t fill_count,
                                     int n_rings, int n_axial,
                                     univ_map_func_t map_func, void* map_ctx,
                                     str_builder_t* sb) {
    int nj = fill_dims[3] - fill_dims[2] + 1;
    int nk = fill_dims[5] - fill_dims[4] + 1;
    int R = n_rings;
    bool first = true;

    for (int m = 0; m < n_axial; m++) {
        int i_x = 1;
        int i_a = R - 1;

        /* Upper triangular region: R-1 rows */
        for (int k = 0; k < R - 1; k++) {
            i_x -= 1;
            for (int j = 0; j <= k; j++) {
                int oi = i_x - fill_dims[0];
                int oj = i_a - fill_dims[2];
                size_t idx = (size_t)(oi * nj * nk + oj * nk + m);
                int uid = (idx < fill_count) ? lat_fill[idx] : 0;
                if (map_func) uid = map_func(uid, map_ctx);
                if (!first) str_builder_putc(sb, ' ');
                str_builder_int(sb, uid);
                first = false;
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
                int oi = i_x - fill_dims[0];
                int oj = i_a - fill_dims[2];
                size_t idx = (size_t)(oi * nj * nk + oj * nk + m);
                int uid = (idx < fill_count) ? lat_fill[idx] : 0;
                if (map_func) uid = map_func(uid, map_ctx);
                if (!first) str_builder_putc(sb, ' ');
                str_builder_int(sb, uid);
                first = false;
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
                int oi = i_x - fill_dims[0];
                int oj = i_a - fill_dims[2];
                size_t idx = (size_t)(oi * nj * nk + oj * nk + m);
                int uid = (idx < fill_count) ? lat_fill[idx] : 0;
                if (map_func) uid = map_func(uid, map_ctx);
                if (!first) str_builder_putc(sb, ' ');
                str_builder_int(sb, uid);
                first = false;
                i_x += 2;
                i_a -= 1;
            }
            i_x -= 2 * row_len;
            i_a += row_len;
        }
    }
}

/**
 * @brief Write geometry section to XML builder
 */
static bool write_geometry_section(const alea_system_t* sys, export_context_t* ctx,
                                    openmc_xml_t* xml, arena_t* arena) {
    /* <geometry> */
    if (!openmc_xml_start_element(xml, "geometry")) return false;
    if (!openmc_xml_end_start_tag(xml, false)) return false;

    /* ========== UNIVERSE REACHABILITY ========== */
    /* Walk the fill graph from universe 0 (root) to find reachable universes.
     * Orphaned universes (defined in MCNP but never filled) are skipped —
     * OpenMC requires exactly one root universe. */
    int max_id = 0;
    for (size_t i = 0; i < alea_vec_count(&sys->cells); i++) {
        const alea_cell_entry_t* c = &sys->cells.data[i];
        if (c->universe_id > max_id) max_id = c->universe_id;
        if (c->fill_universe > max_id) max_id = c->fill_universe;
        if (c->mc_cell_id > max_id) max_id = c->mc_cell_id;
        for (size_t u = 0; u < c->lat_fill_count; u++) {
            if (c->lat_fill[u] > max_id) max_id = c->lat_fill[u];
        }
    }

    alea_bitset_t reachable = alea_bitset_create(max_id + 1);
    if (!reachable.words) return false;
    alea_bitset_set(&reachable, 0);  /* root universe */

    for (bool changed = true; changed; ) {
        changed = false;
        for (size_t i = 0; i < alea_vec_count(&sys->cells); i++) {
            const alea_cell_entry_t* c = &sys->cells.data[i];
            if (!alea_bitset_test(&reachable, c->universe_id)) continue;

            /* Regular fill */
            if (c->fill_universe > 0 && !alea_bitset_test(&reachable, c->fill_universe)) {
                alea_bitset_set(&reachable, c->fill_universe);
                changed = true;
            }
            /* Lattice fill — all universes in the fill array */
            if (c->lat_type > 0 && c->lat_fill) {
                for (size_t u = 0; u < c->lat_fill_count; u++) {
                    int uid = c->lat_fill[u];
                    if (uid >= 0 && uid <= max_id && !alea_bitset_test(&reachable, uid)) {
                        alea_bitset_set(&reachable, uid);
                        changed = true;
                    }
                }
            }
        }
    }

    {
        int n_reachable = 0, n_total = 0;
        for (size_t i = 0; i < alea_vec_count(&sys->cells); i++) {
            int uid = sys->cells.data[i].universe_id;
            n_total++;
            if (alea_bitset_test(&reachable, uid)) n_reachable++;
        }
        if (n_reachable < n_total) {
            ALEA_LOG_WARN("Pruning %d/%d cells in orphaned universes (not reachable from root)",
                    n_total - n_reachable, n_total);
            /* List each unreachable universe with its cell count */
            alea_bitset_t warned = alea_bitset_create(max_id + 1);
            for (size_t i = 0; i < alea_vec_count(&sys->cells); i++) {
                int uid = sys->cells.data[i].universe_id;
                if (!alea_bitset_test(&reachable, uid) && (!warned.words || !alea_bitset_test(&warned, uid))) {
                    int count = 0;
                    for (size_t j = i; j < alea_vec_count(&sys->cells); j++) {
                        if (sys->cells.data[j].universe_id == uid) count++;
                    }
                    ALEA_LOG_WARN("  Universe %d: %d cells (not filled by any cell)", uid, count);
                    if (warned.words) alea_bitset_set(&warned, uid);
                }
            }
            alea_bitset_destroy(&warned);
        }
    }

    /* ========== LATTICES ========== */
    /* In OpenMC, a lattice IS a universe: lattice ID = cell's universe_id.
     * Cells referencing fill=U where U is the lattice cell's universe will
     * directly fill with the lattice.  The lattice cell itself is NOT exported
     * as an OpenMC cell (the reference openmc_mcnp_adapter does the same).
     *
     * When the lattice element center is not at (0,0,0), we create wrapper
     * universes with translated cells, matching openmc_mcnp_adapter behavior.
     */
    int next_id = max_id + 1;

    /* Wrapper cells for lattice center-offset translation */
    typedef struct {
        int wrapper_univ_id;
        int original_univ_id;
        int wrapper_cell_id;
        double translation[3];
    } wrapper_cell_t;

    /* Deferred lattice emission data */
    typedef struct {
        int lattice_id;
        int mc_cell_id;
        double pitch[3];
        int dims[3];
        double lower_left[3];
        int is_3d;
        int outer_universe;
        const char* universes_str;  /* pre-built, arena-allocated */
        int lat_type;               /* 1=rect, 2=hex */
        int fill_dims[6];          /* imin, imax, jmin, jmax, kmin, kmax */
        double center[3];          /* lattice center (for hex) */
    } lattice_emit_t;

    size_t wrapper_count = 0;
    size_t wrapper_cap = 64;
    wrapper_cell_t* wrappers = (wrapper_cell_t*)arena_alloc(arena,
            wrapper_cap * sizeof(wrapper_cell_t));

    size_t lattice_count = 0;
    size_t lattice_cap = 16;
    lattice_emit_t* lattices = (lattice_emit_t*)arena_alloc(arena,
            lattice_cap * sizeof(lattice_emit_t));

    for (size_t i = 0; i < alea_vec_count(&sys->cells); i++) {
        const alea_cell_entry_t* cell = &sys->cells.data[i];

        /* Skip non-lattice cells */
        if (cell->lat_type == 0 || !cell->lat_fill || cell->lat_fill_count == 0) {
            continue;
        }

        /* Lattice ID = the cell's universe (the lattice replaces this universe) */
        int lattice_id = cell->universe_id;
        if (lattice_id <= 0) {
            ALEA_LOG_WARN("Lattice cell %d has no universe assignment, skipping",
                    cell->mc_cell_id);
            continue;
        }

        /* Skip lattices in unreachable universes */
        if (!alea_bitset_test(&reachable, lattice_id)) continue;

        /* Calculate dimensions from lat_fill_dims */
        int ni = cell->lat_fill_dims[1] - cell->lat_fill_dims[0] + 1;
        int nj = cell->lat_fill_dims[3] - cell->lat_fill_dims[2] + 1;
        int nk = cell->lat_fill_dims[5] - cell->lat_fill_dims[4] + 1;
        if (ni <= 0) ni = 1;
        if (nj <= 0) nj = 1;
        if (nk <= 0) nk = 1;

        /* Use stored pitch from MCNP conversion (computed from cell bbox) */
        double pitch_x = cell->lat_pitch[0];
        double pitch_y = cell->lat_pitch[1];
        double pitch_z = cell->lat_pitch[2];
        if (pitch_x <= 0) pitch_x = 1.0;
        if (pitch_y <= 0) pitch_y = 1.0;
        if (pitch_z <= 0) pitch_z = 1.0;

        double ll_x = cell->lat_lower_left[0];
        double ll_y = cell->lat_lower_left[1];
        double ll_z = cell->lat_lower_left[2];

        /* 3D if z-pitch is finite (cell has z-bounding surfaces) */
        int is_3d = (pitch_z > 0 && pitch_z < 1e8);

        /* Compute element center from bbox (midpoint of bounding surfaces).
         * If non-zero, fill universes need translation = -center so that
         * their local coordinates align with the element frame. */
        double center[3] = {0, 0, 0};
        if (cell->root_node_id != ALEA_NODE_ID_INVALID &&
            cell->root_node_id < alea_vec_count(&sys->nodes)) {
            alea_bbox_t bbox = sys->nodes.data[cell->root_node_id].bbox;
            center[0] = (bbox.min_x + bbox.max_x) / 2.0;
            center[1] = (bbox.min_y + bbox.max_y) / 2.0;
            if (is_3d)
                center[2] = (bbox.min_z + bbox.max_z) / 2.0;
        }
        int needs_translation = (fabs(center[0]) > 1e-10 ||
                                 fabs(center[1]) > 1e-10 ||
                                 fabs(center[2]) > 1e-10);

        /* Build universe mapping: original_uid -> wrapper_uid (if translating).
         * Small inline lookup — lattice fill arrays rarely exceed ~100 unique IDs. */
        int univ_map_orig[256];
        int univ_map_wrap[256];
        int univ_map_count = 0;

        if (needs_translation) {
            /* Collect unique fill universe IDs and assign wrapper IDs */
            for (size_t u = 0; u < cell->lat_fill_count; u++) {
                int uid = cell->lat_fill[u];
                int found = 0;
                for (int m = 0; m < univ_map_count; m++) {
                    if (univ_map_orig[m] == uid) { found = 1; break; }
                }
                if (!found && univ_map_count < 256) {
                    int wrapper_univ = next_id++;
                    int wrapper_cell = next_id++;
                    univ_map_orig[univ_map_count] = uid;
                    univ_map_wrap[univ_map_count] = wrapper_univ;
                    univ_map_count++;

                    /* Record wrapper cell for emission later */
                    if (wrapper_count >= wrapper_cap) {
                        wrapper_cap *= 2;
                        wrapper_cell_t* new_w = (wrapper_cell_t*)arena_alloc(arena,
                                wrapper_cap * sizeof(wrapper_cell_t));
                        memcpy(new_w, wrappers, wrapper_count * sizeof(wrapper_cell_t));
                        wrappers = new_w;
                    }
                    wrappers[wrapper_count].wrapper_univ_id = wrapper_univ;
                    wrappers[wrapper_count].original_univ_id = uid;
                    wrappers[wrapper_count].wrapper_cell_id = wrapper_cell;
                    wrappers[wrapper_count].translation[0] = -center[0];
                    wrappers[wrapper_count].translation[1] = -center[1];
                    wrappers[wrapper_count].translation[2] = -center[2];
                    wrapper_count++;
                }
            }
        }

        /* Helper: map a fill universe through the wrapper table */
        #define MAP_UNIV(uid) \
            ({ int _r = (uid); \
               if (needs_translation) { \
                   for (int _m = 0; _m < univ_map_count; _m++) { \
                       if (univ_map_orig[_m] == _r) { _r = univ_map_wrap[_m]; break; } \
                   } \
               } \
               _r; })

        /* Determine outer universe (for infinite / simple-fill lattices) */
        int outer_universe = -1;
        if (cell->lat_fill_count == 1) {
            outer_universe = MAP_UNIV(cell->lat_fill[0]);
        }

        /* Build mapped copy of lat_fill (applying wrapper universe mapping) */
        int* mapped_fill = (int*)arena_alloc(arena,
                cell->lat_fill_count * sizeof(int));
        for (size_t u = 0; u < cell->lat_fill_count; u++) {
            mapped_fill[u] = MAP_UNIV(cell->lat_fill[u]);
        }

        #undef MAP_UNIV

        /* Pre-build universes string (arena-allocated, survives this scope) */
        str_builder_t univ_sb;
        str_builder_init(&univ_sb, arena, 4096);

        if (cell->lat_type == 2) {
            /* Hex: emit in OpenMC visual row order */
            int n_rings = compute_n_rings(cell->lat_fill_dims);
            hex_grid_to_visual_rows(mapped_fill, cell->lat_fill_dims,
                                    cell->lat_fill_count,
                                    n_rings, nk, NULL, NULL, &univ_sb);
        } else {
            /* Rect: flat order */
            for (size_t u = 0; u < cell->lat_fill_count; u++) {
                if (u > 0) str_builder_putc(&univ_sb, ' ');
                str_builder_int(&univ_sb, mapped_fill[u]);
            }
        }
        const char* universes_str = str_builder_get(&univ_sb);

        /* Store lattice data for deferred emission (after regular cells) */
        if (lattice_count >= lattice_cap) {
            lattice_cap *= 2;
            lattice_emit_t* new_l = (lattice_emit_t*)arena_alloc(arena,
                    lattice_cap * sizeof(lattice_emit_t));
            memcpy(new_l, lattices, lattice_count * sizeof(lattice_emit_t));
            lattices = new_l;
        }
        lattice_emit_t* le = &lattices[lattice_count++];
        le->lattice_id = lattice_id;
        le->mc_cell_id = cell->mc_cell_id;
        le->pitch[0] = pitch_x; le->pitch[1] = pitch_y; le->pitch[2] = pitch_z;
        le->dims[0] = ni; le->dims[1] = nj; le->dims[2] = nk;
        le->lower_left[0] = ll_x; le->lower_left[1] = ll_y; le->lower_left[2] = ll_z;
        le->is_3d = is_3d;
        le->outer_universe = outer_universe;
        le->universes_str = universes_str;
        le->lat_type = cell->lat_type;
        memcpy(le->fill_dims, cell->lat_fill_dims, sizeof(le->fill_dims));
        le->center[0] = center[0];
        le->center[1] = center[1];
        le->center[2] = center[2];

        ALEA_LOG_DEBUG("Prepared lattice %d from cell %d: %dx%dx%d%s",
                     lattice_id, cell->mc_cell_id, ni, nj, nk,
                     needs_translation ? " (with center translation)" : "");
    }

    /* ========== CELLS ========== */
    for (size_t i = 0; i < alea_vec_count(&sys->cells); i++) {
        if (g_alea_interrupted) return false;
        const alea_cell_entry_t* cell = &sys->cells.data[i];

        /* Skip cells in unreachable universes (orphaned in MCNP model) */
        if (!alea_bitset_test(&reachable, cell->universe_id)) continue;

        /* Skip lattice-defining cells — replaced by <lattice> element above */
        if (cell->lat_type > 0 && cell->lat_fill && cell->lat_fill_count > 0) {
            continue;
        }

        /* Skip FILL cells without geometry */
        if (cell->fill_universe > 0 && cell->root_node_id == ALEA_NODE_ID_INVALID) {
            continue;
        }

        /* Build region expression */
        str_builder_t region_sb;
        str_builder_init(&region_sb, arena, 4096);

        if (cell->root_node_id != ALEA_NODE_ID_INVALID) {
            if (!tree_to_openmc_region(sys, ctx, cell->root_node_id, &region_sb, OPENMC_CTX_TOP)) {
                ALEA_LOG_WARN("Failed to convert CSG tree for cell %d",
                        cell->mc_cell_id);
                continue;
            }
        }

        const char* region = str_builder_get(&region_sb);
        size_t region_len = region_sb.len;

        /* Pre-compute translation and rotation strings for fill cells */
        char trans_str[128] = {0};
        char rot_str[256] = {0};
        if (cell->fill_universe > 0 && cell->fill_transform > 0) {
            const alea_transform_t* tr = alea_get_transform(sys, cell->fill_transform);
            double zt = sys->config.zero_threshold;
            if (tr && tr->value_count >= 3) {
                snprintf(trans_str, sizeof(trans_str), "%.10g %.10g %.10g",
                         snap_zero(tr->data[0], zt),
                         snap_zero(tr->data[1], zt),
                         snap_zero(tr->data[2], zt));
                if (tr->value_count >= 12) {
                    snprintf(rot_str, sizeof(rot_str),
                             "%.10g %.10g %.10g %.10g %.10g %.10g %.10g %.10g %.10g",
                             snap_zero(tr->cosines[3], zt),
                             snap_zero(tr->cosines[4], zt),
                             snap_zero(tr->cosines[5], zt),
                             snap_zero(tr->cosines[6], zt),
                             snap_zero(tr->cosines[7], zt),
                             snap_zero(tr->cosines[8], zt),
                             snap_zero(tr->cosines[9], zt),
                             snap_zero(tr->cosines[10], zt),
                             snap_zero(tr->cosines[11], zt));
                }
            }
        }

        /* <cell ...> — attribute order: fill, id, material, region, rotation, translation, universe */
        xml->no_wrap = true;
        if (!openmc_xml_start_element(xml, "cell")) return false;

        /* fill */
        if (cell->fill_universe > 0) {
            if (!openmc_xml_attribute_i(xml, "fill", cell->fill_universe)) return false;
        }

        /* id */
        if (!openmc_xml_attribute_i(xml, "id", cell->mc_cell_id)) return false;

        /* material — skip for fill cells (OpenMC: material and fill are exclusive) */
        if (cell->fill_universe == 0) {
            if (cell->material_id == 0) {
                if (!openmc_xml_attribute(xml, "material", "void")) return false;
            } else {
                int openmc_mat_id = cell->material_id;
                if (ctx->mat_map) {
                    const mat_density_map_t* map = ctx->mat_map;
                    for (size_t m = 0; m < map->count; m++) {
                        if (map->entries[m].mc_material_id == cell->material_id &&
                            fabs(map->entries[m].density - cell->density) < DENSITY_TOLERANCE) {
                            openmc_mat_id = map->entries[m].openmc_material_id;
                            break;
                        }
                    }
                }
                char mat_str[32];
                snprintf(mat_str, sizeof(mat_str), "%d", openmc_mat_id);
                if (!openmc_xml_attribute(xml, "material", mat_str)) return false;
            }
        }

        /* region */
        if (region_len > 0) {
            if (!openmc_xml_attribute(xml, "region", region)) return false;
        }

        /* rotation */
        if (rot_str[0]) {
            if (!openmc_xml_attribute(xml, "rotation", rot_str)) return false;
        }

        /* translation */
        if (trans_str[0]) {
            if (!openmc_xml_attribute(xml, "translation", trans_str)) return false;
        }

        /* universe */
        if (!openmc_xml_attribute_i(xml, "universe", cell->universe_id)) return false;

        if (!openmc_xml_end_start_tag(xml, true)) return false;
        xml->no_wrap = false;
        ctx->cells_written++;
    }

    alea_bitset_destroy(&reachable);

    /* ========== WRAPPER CELLS for lattice center-offset ========== */
    for (size_t w = 0; w < wrapper_count; w++) {
        xml->no_wrap = true;
        if (!openmc_xml_start_element(xml, "cell")) return false;
        if (!openmc_xml_attribute_i(xml, "fill", wrappers[w].original_univ_id)) return false;
        if (!openmc_xml_attribute_i(xml, "id", wrappers[w].wrapper_cell_id)) return false;
        {
            char buf[128];
            snprintf(buf, sizeof(buf), "%.10g %.10g %.10g",
                     wrappers[w].translation[0],
                     wrappers[w].translation[1],
                     wrappers[w].translation[2]);
            if (!openmc_xml_attribute(xml, "translation", buf)) return false;
        }
        if (!openmc_xml_attribute_i(xml, "universe", wrappers[w].wrapper_univ_id)) return false;
        if (!openmc_xml_end_start_tag(xml, true)) return false;
        xml->no_wrap = false;
        ctx->cells_written++;
    }

    /* ========== LATTICES (deferred from computation pass) ========== */
    for (size_t li = 0; li < lattice_count; li++) {
        const lattice_emit_t* le = &lattices[li];

        if (le->lat_type == 2) {
            /* ---- HEX LATTICE ---- */
            int n_rings = compute_n_rings(le->fill_dims);

            if (!openmc_xml_start_element(xml, "hex_lattice")) return false;
            if (!openmc_xml_attribute_i(xml, "id", le->lattice_id)) return false;
            {
                char lat_name[64];
                snprintf(lat_name, sizeof(lat_name), "lattice_cell_%d", le->mc_cell_id);
                if (!openmc_xml_attribute(xml, "name", lat_name)) return false;
            }
            if (!openmc_xml_attribute_i(xml, "n_rings", n_rings)) return false;
            if (le->is_3d) {
                if (!openmc_xml_attribute_i(xml, "n_axial", le->dims[2])) return false;
            }
            if (le->outer_universe >= 0) {
                if (!openmc_xml_attribute_i(xml, "outer", le->outer_universe)) return false;
            }
            if (!openmc_xml_end_start_tag(xml, false)) return false;

            /* <center> */
            if (!openmc_xml_start_element(xml, "center")) return false;
            {
                char buf[128];
                if (le->is_3d)
                    snprintf(buf, sizeof(buf), "%.10g %.10g %.10g",
                             le->center[0], le->center[1], le->center[2]);
                else
                    snprintf(buf, sizeof(buf), "%.10g %.10g",
                             le->center[0], le->center[1]);
                if (!openmc_xml_text(xml, buf)) return false;
            }
            if (!openmc_xml_end_element(xml, "center")) return false;

            /* <pitch> — hex uses single pitch (2D) or two values (3D) */
            if (!openmc_xml_start_element(xml, "pitch")) return false;
            {
                char buf[128];
                if (le->is_3d)
                    snprintf(buf, sizeof(buf), "%.10g %.10g",
                             le->pitch[0], le->pitch[2]);
                else
                    snprintf(buf, sizeof(buf), "%.10g", le->pitch[0]);
                if (!openmc_xml_text(xml, buf)) return false;
            }
            if (!openmc_xml_end_element(xml, "pitch")) return false;

            /* <universes> — visual row order */
            if (!openmc_xml_start_element(xml, "universes")) return false;
            if (!openmc_xml_text(xml, le->universes_str)) return false;
            if (!openmc_xml_end_element(xml, "universes")) return false;

            if (!openmc_xml_end_element(xml, "hex_lattice")) return false;
        } else {
            /* ---- RECTANGULAR LATTICE ---- */
            if (!openmc_xml_start_element(xml, "lattice")) return false;
            if (!openmc_xml_attribute_i(xml, "id", le->lattice_id)) return false;
            {
                char lat_name[64];
                snprintf(lat_name, sizeof(lat_name), "lattice_cell_%d", le->mc_cell_id);
                if (!openmc_xml_attribute(xml, "name", lat_name)) return false;
            }
            if (!openmc_xml_end_start_tag(xml, false)) return false;

            /* <pitch> */
            if (!openmc_xml_start_element(xml, "pitch")) return false;
            {
                char buf[128];
                if (le->is_3d)
                    snprintf(buf, sizeof(buf), "%.10g %.10g %.10g",
                             le->pitch[0], le->pitch[1], le->pitch[2]);
                else
                    snprintf(buf, sizeof(buf), "%.10g %.10g",
                             le->pitch[0], le->pitch[1]);
                if (!openmc_xml_text(xml, buf)) return false;
            }
            if (!openmc_xml_end_element(xml, "pitch")) return false;

            /* <dimension> */
            if (!openmc_xml_start_element(xml, "dimension")) return false;
            {
                char buf[64];
                if (le->is_3d)
                    snprintf(buf, sizeof(buf), "%d %d %d",
                             le->dims[0], le->dims[1], le->dims[2]);
                else
                    snprintf(buf, sizeof(buf), "%d %d",
                             le->dims[0], le->dims[1]);
                if (!openmc_xml_text(xml, buf)) return false;
            }
            if (!openmc_xml_end_element(xml, "dimension")) return false;

            /* <lower_left> */
            if (!openmc_xml_start_element(xml, "lower_left")) return false;
            {
                char buf[128];
                if (le->is_3d)
                    snprintf(buf, sizeof(buf), "%.10g %.10g %.10g",
                             le->lower_left[0], le->lower_left[1], le->lower_left[2]);
                else
                    snprintf(buf, sizeof(buf), "%.10g %.10g",
                             le->lower_left[0], le->lower_left[1]);
                if (!openmc_xml_text(xml, buf)) return false;
            }
            if (!openmc_xml_end_element(xml, "lower_left")) return false;

            /* <outer> */
            if (le->outer_universe >= 0) {
                if (!openmc_xml_start_element(xml, "outer")) return false;
                {
                    char buf[32];
                    snprintf(buf, sizeof(buf), "%d", le->outer_universe);
                    if (!openmc_xml_text(xml, buf)) return false;
                }
                if (!openmc_xml_end_element(xml, "outer")) return false;
            }

            /* <universes> */
            if (!openmc_xml_start_element(xml, "universes")) return false;
            if (!openmc_xml_text(xml, le->universes_str)) return false;
            if (!openmc_xml_end_element(xml, "universes")) return false;

            if (!openmc_xml_end_element(xml, "lattice")) return false;
        }
    }

    /* ========== SURFACES ========== */
    /* Track emitted primitives to avoid duplicates (surfaces can share primitives) */
    alea_bitset_t prim_emitted = alea_bitset_create(alea_vec_count(&sys->primitives) + 1);
    if (!prim_emitted.words) return false;

    for (size_t i = 0; i < alea_vec_count(&sys->surfaces); i++) {
        if (g_alea_interrupted) { alea_bitset_destroy(&prim_emitted); return false; }
        const alea_surface_entry_t* surf = &sys->surfaces.data[i];
        int surf_id = surf->mc_surface_id;

        uint32_t prim_id = surf->primitive_id;
        if (prim_id == ALEA_PRIMITIVE_ID_INVALID || prim_id >= alea_vec_count(&sys->primitives)) continue;

        /* When dedup is enabled, only emit the canonical surface for each primitive
           (the one whose mc_surface_id matches prim_to_surface[prim_id]).
           Cell expressions use prim_to_surface[prim_id], so we must emit that ID. */
        if (ctx->deduplicate && ctx->prim_to_surface &&
            prim_id < ctx->prim_to_surface_size &&
            ctx->prim_to_surface[prim_id] >= 0) {
            if (surf_id != ctx->prim_to_surface[prim_id]) continue;
        }

        /* Only skip duplicate primitives when dedup is enabled — without dedup,
           each surface entry must be emitted with its own mc_surface_id */
        if (ctx->deduplicate) {
            if (alea_bitset_test(&prim_emitted, prim_id)) continue;
            alea_bitset_set(&prim_emitted, prim_id);
        }

        const alea_primitive_entry_t* prim = &sys->primitives.data[prim_id];

        /* Un-canonicalize: restore original surface orientation */
        int8_t inverted = 0;
        if (surf->pos_node < alea_vec_count(&sys->nodes)) {
            inverted = sys->nodes.data[surf->pos_node].primitive.inverted;
        }

        /* Skip 1-sheet cones — they are expanded into (double-sheet cone ∩ plane)
           and the original surface ID is not referenced in cell expressions */
        if ((prim->type == ALEA_PRIMITIVE_CONE_X && prim->data.cone_x.sheet_selection != 0) ||
            (prim->type == ALEA_PRIMITIVE_CONE_Y && prim->data.cone_y.sheet_selection != 0) ||
            (prim->type == ALEA_PRIMITIVE_CONE_Z && prim->data.cone_z.sheet_selection != 0)) {
            continue;
        }

        const char* type_name = get_openmc_surface_type(prim, inverted);

        if (!type_name) {
            if (is_macrobody_type(prim->type)) {
                ALEA_LOG_INFO("Macrobody %s surface %d expanded into constituent surfaces",
                        macrobody_type_name(prim->type), surf_id);
            } else {
                ALEA_LOG_WARN("Skipping surface %d with unknown type %d",
                        surf_id, prim->type);
            }
            continue;
        }

        char coeffs[512];
        if (format_surface_coeffs(prim, coeffs, sizeof(coeffs), inverted) < 0) {
            ALEA_LOG_WARN("Failed to format coefficients for surface %d", surf_id);
            continue;
        }

        /* <surface boundary="..." coeffs="..." id="N" type="TYPE"/> */
        xml->no_wrap = true;
        if (!openmc_xml_start_element(xml, "surface")) { alea_bitset_destroy(&prim_emitted); return false; }

        /* Boundary condition (first, if present) */
        if (surf->boundary_type == ALEA_BOUNDARY_VACUUM) {
            if (!openmc_xml_attribute(xml, "boundary", "vacuum")) { alea_bitset_destroy(&prim_emitted); return false; }
        } else if (surf->boundary_type == ALEA_BOUNDARY_WHITE) {
            if (!openmc_xml_attribute(xml, "boundary", "white")) { alea_bitset_destroy(&prim_emitted); return false; }
        } else if (surf->boundary_type == ALEA_BOUNDARY_REFLECTIVE) {
            if (!openmc_xml_attribute(xml, "boundary", "reflective")) { alea_bitset_destroy(&prim_emitted); return false; }
        } else if (surf->boundary_type == ALEA_BOUNDARY_PERIODIC) {
            if (!openmc_xml_attribute(xml, "boundary", "periodic")) { alea_bitset_destroy(&prim_emitted); return false; }
        }

        if (!openmc_xml_attribute(xml, "coeffs", coeffs)) { alea_bitset_destroy(&prim_emitted); return false; }
        if (!openmc_xml_attribute_i(xml, "id", surf_id)) { alea_bitset_destroy(&prim_emitted); return false; }
        if (!openmc_xml_attribute(xml, "type", type_name)) { alea_bitset_destroy(&prim_emitted); return false; }

        if (!openmc_xml_end_start_tag(xml, true)) { alea_bitset_destroy(&prim_emitted); return false; }
        xml->no_wrap = false;
        ctx->surfaces_written++;
    }

    alea_bitset_destroy(&prim_emitted);

    /* ========== SYNTHETIC SURFACES (from macrobody expansion) ========== */
    /* Build set of registered surface IDs for quick lookup */
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

    /* Track written to avoid duplicates */
    int* written_synth = calloc(ctx->next_synthetic_surface_id + 1, sizeof(int));
    if (written_synth) {
        for (size_t n = 0; n < alea_vec_count(&sys->nodes); n++) {
            const alea_node_t* node = &sys->nodes.data[n];
            if (ALEA_GET_OPERATION(node) != ALEA_OP_PRIMITIVE) continue;

            int surf_id = node->primitive.mc_surface_id;
            if (surf_id <= 0) continue;
            if (written_synth[surf_id]) continue;
            if (registered && surf_id <= max_registered_id && registered[surf_id]) continue;

            uint32_t prim_id = node->primitive.primitive_id;
            if (prim_id >= alea_vec_count(&sys->primitives)) continue;

            const alea_primitive_entry_t* prim = &sys->primitives.data[prim_id];
            const char* type_name = get_openmc_surface_type(prim, node->primitive.inverted);

            if (!type_name) {
                if (is_macrobody_type(prim->type)) {
                    ALEA_LOG_INFO("Macrobody %s surface %d expanded into constituent surfaces",
                            macrobody_type_name(prim->type), surf_id);
                } else {
                    ALEA_LOG_WARN("Skipping synthetic surface %d with unknown type %d",
                            surf_id, prim->type);
                }
                continue;
            }

            char coeffs[512];
            if (format_surface_coeffs(prim, coeffs, sizeof(coeffs), node->primitive.inverted) < 0) {
                ALEA_LOG_WARN("Failed to format coefficients for synthetic surface %d", surf_id);
                continue;
            }

            xml->no_wrap = true;
            if (!openmc_xml_start_element(xml, "surface")) goto synth_fail;
            if (!openmc_xml_attribute(xml, "coeffs", coeffs)) goto synth_fail;
            if (!openmc_xml_attribute_i(xml, "id", surf_id)) goto synth_fail;
            if (!openmc_xml_attribute(xml, "type", type_name)) goto synth_fail;
            if (!openmc_xml_end_start_tag(xml, true)) goto synth_fail;
            xml->no_wrap = false;
            ctx->surfaces_written++;
            written_synth[surf_id] = 1;
        }
        free(written_synth);
    }
    free(registered);
    goto synth_done;

synth_fail:
    free(written_synth);
    free(registered);
    return false;

synth_done:

    /* </geometry> */
    if (!openmc_xml_end_element(xml, "geometry")) return false;

    return true;
}

/**
 * @brief Write settings section to XML builder
 */
static bool write_settings_section(const alea_system_t* sys, openmc_xml_t* xml) {
    (void)sys;  /* May use system info in future */

    /* <settings> */
    if (!openmc_xml_start_element(xml, "settings")) return false;
    if (!openmc_xml_end_start_tag(xml, false)) return false;

    /* <run_mode>eigenvalue</run_mode> - default mode */
    if (!openmc_xml_start_element(xml, "run_mode")) return false;
    if (!openmc_xml_text(xml, "eigenvalue")) return false;
    if (!openmc_xml_end_element(xml, "run_mode")) return false;

    /* <particles>10000</particles> - reasonable default */
    if (!openmc_xml_start_element(xml, "particles")) return false;
    if (!openmc_xml_text(xml, "10000")) return false;
    if (!openmc_xml_end_element(xml, "particles")) return false;

    /* <batches>100</batches> */
    if (!openmc_xml_start_element(xml, "batches")) return false;
    if (!openmc_xml_text(xml, "100")) return false;
    if (!openmc_xml_end_element(xml, "batches")) return false;

    /* <inactive>10</inactive> */
    if (!openmc_xml_start_element(xml, "inactive")) return false;
    if (!openmc_xml_text(xml, "10")) return false;
    if (!openmc_xml_end_element(xml, "inactive")) return false;

    /* </settings> */
    if (!openmc_xml_end_element(xml, "settings")) return false;

    return true;
}

/**
 * @brief Assign surface IDs to primitives that don't have them (from macrobody expansion)
 *
 * Same logic as MCNP export - uses primitive-to-surface-id map for consistency.
 */
static void assign_surface_ids_recursive_openmc(alea_system_t* sys, int* next_id,
                                                 alea_node_id_t node_id,
                                                 int** prim_surf_map, size_t* map_size) {
    if (node_id == ALEA_NODE_ID_INVALID || node_id >= alea_vec_count(&sys->nodes)) return;

    alea_node_t* node = &sys->nodes.data[node_id];
    alea_operation_t op = ALEA_GET_OPERATION(node);

    if (op == ALEA_OP_PRIMITIVE) {
        if (node->primitive.mc_surface_id == 0) {
            alea_primitive_id_t prim_id = node->primitive.primitive_id;

            /* Grow map if needed */
            if (prim_id >= *map_size) {
                size_t new_size = (prim_id + 1) * 2;
                int* new_map = realloc(*prim_surf_map, new_size * sizeof(int));
                if (!new_map) return;
                for (size_t i = *map_size; i < new_size; i++) {
                    new_map[i] = 0;
                }
                *prim_surf_map = new_map;
                *map_size = new_size;
            }

            /* Check if we already assigned an ID to this primitive */
            if ((*prim_surf_map)[prim_id] != 0) {
                node->primitive.mc_surface_id = (*prim_surf_map)[prim_id];
            } else {
                int new_id = (*next_id)++;
                node->primitive.mc_surface_id = new_id;
                (*prim_surf_map)[prim_id] = new_id;

                /* Register a surface entry so the canonical surface map
                   and surface emission can find this primitive */
                alea_surface_entry_t* entry = alea_vec_push_uninit(&sys->surfaces, alea_surface_entry_t);
                if (entry) {
                    memset(entry, 0, sizeof(*entry));
                    entry->mc_surface_id = new_id;
                    entry->primitive_id = prim_id;
                    entry->pos_node = node_id;
                    entry->neg_node = ALEA_NODE_ID_INVALID;
                    entry->expanded_pos_node = ALEA_NODE_ID_INVALID;
                    entry->expanded_neg_node = ALEA_NODE_ID_INVALID;
                }
            }
        }
        return;
    }

    /* Recurse into children */
    if (op == ALEA_OP_COMPLEMENT) {
        assign_surface_ids_recursive_openmc(sys, next_id, node->operation.left, prim_surf_map, map_size);
    } else if (op == ALEA_OP_UNION || op == ALEA_OP_INTERSECTION || op == ALEA_OP_DIFFERENCE) {
        assign_surface_ids_recursive_openmc(sys, next_id, node->operation.left, prim_surf_map, map_size);
        assign_surface_ids_recursive_openmc(sys, next_id, node->operation.right, prim_surf_map, map_size);
    }
}

static void assign_missing_surface_ids_openmc(alea_system_t* sys, int* next_id) {
    int* prim_surf_map = NULL;
    size_t map_size = 0;

    for (size_t i = 0; i < alea_vec_count(&sys->cells); i++) {
        alea_cell_entry_t* cell = &sys->cells.data[i];
        if (cell->root_node_id == ALEA_NODE_ID_INVALID) continue;

        assign_surface_ids_recursive_openmc(sys, next_id, cell->root_node_id, &prim_surf_map, &map_size);
    }

    free(prim_surf_map);
}

int export_openmc(const alea_system_t* sys, export_context_t* ctx) {
    if (!sys || !ctx || !ctx->out) return -1;

    arena_t arena;
    arena_init(&arena);

    /* OpenMC doesn't support macrobodies - replace them with expanded primitives */
    /* This uses pre-computed expanded nodes from surface entries */
    int expanded = alea_expand_macrobodies_in_cells((alea_system_t*)sys);
    if (expanded > 0) {
        ALEA_LOG_DEBUG("Expanded %d cells with macrobodies for OpenMC export", expanded);
    }

    /* Assign surface IDs to primitives from immediate macrobody expansion */
    /* Uses ctx->next_synthetic_surface_id which is set from max surface ID in alea_api.c */
    assign_missing_surface_ids_openmc((alea_system_t*)sys, &ctx->next_synthetic_surface_id);

    /* Build canonical surface map so cell expressions use deduplicated surface IDs */
    if (ctx->deduplicate) {
        alea_build_canonical_surface_map(ctx, sys);
    }

    /* Build material-density map for MCNP->OpenMC density conversion */
    /* In MCNP, density is per-cell; in OpenMC, density is per-material */
    mat_density_map_t mat_map;
    build_mat_density_map(sys, &mat_map);
    ctx->mat_map = &mat_map;

    if (mat_map.count > alea_vec_count(&sys->materials)) {
        ALEA_LOG_INFO("Created %zu OpenMC materials from %zu MCNP materials "
                    "(some materials used with multiple densities)",
                    mat_map.count, alea_vec_count(&sys->materials));
    }

    openmc_xml_t xml;
    openmc_xml_init_ex(&xml, &arena, 16384, 100, 2, true);

    /* XML declaration */
    if (!openmc_xml_start_document(&xml, "1.0", "utf-8")) {
        mat_density_map_free(&mat_map);
        arena_free(&arena);
        return -1;
    }

    /* <model> root element */
    if (!openmc_xml_start_element(&xml, "model")) goto fail;
    if (!openmc_xml_end_start_tag(&xml, false)) goto fail;

    /* Materials section */
    if (!write_materials_section(sys, &xml, &arena, &mat_map)) goto fail;

    /* Geometry section */
    if (!write_geometry_section(sys, ctx, &xml, &arena)) goto fail;

    /* Settings section */
    if (!write_settings_section(sys, &xml)) goto fail;

    /* </model> */
    if (!openmc_xml_end_element(&xml, "model")) goto fail;

    /* Write to file */
    if (!openmc_xml_write(&xml, ctx->out)) goto fail;

    mat_density_map_free(&mat_map);
    arena_free(&arena);
    return 0;

fail:
    mat_density_map_free(&mat_map);
    arena_free(&arena);
    return -1;
}
