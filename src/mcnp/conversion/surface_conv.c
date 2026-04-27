// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/**
 * @file surface_conv.c
 * @brief MCNP surface to CSG primitive conversion - Updated for flattened architecture
 * 
 * Uses alea_system_t with automatic deduplication and growth.
 */


#include "surface_conv.h"
#include "mcnp/parser/mcnp_parser.h"
#include "core/alea_macrobody.h"
#include "primitives/bbox.h"
#include "util/alea_log.h"
#include "util/compat.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>

/* Error function declared in alea_system.h (thread-local) */

// ============================================================================
// HELPER FUNCTIONS
// ============================================================================

/**
 * @brief Parse coefficient from string
 */
static double parse_coefficient(const char** cursor) {
    while (**cursor && isspace(**cursor)) (*cursor)++;
    if (!**cursor) return 0.0;
    
    char* endptr;
    double value = strtod(*cursor, &endptr);
    *cursor = endptr;
    
    return value;
}

static bool torus_is_sphere(double major_radius, double axial_semiwidth,
                            double minor_radius) {
    if (axial_semiwidth <= 0.0 || minor_radius <= 0.0) return false;

    double scale = fmax(fabs(axial_semiwidth), fabs(minor_radius));
    if (scale < 1.0) scale = 1.0;
    const double eps = 1e-12 * scale;

    return fabs(major_radius) <= eps &&
           fabs(axial_semiwidth - minor_radius) <= eps;
}

static void set_sphere_from_torus(alea_primitive_type_t* out_type,
                                  alea_primitive_data_t* out_data,
                                  double cx, double cy, double cz,
                                  double radius) {
    *out_type = ALEA_PRIMITIVE_SPHERE;
    out_data->sphere.center_x = cx;
    out_data->sphere.center_y = cy;
    out_data->sphere.center_z = cz;
    out_data->sphere.radius = radius;
}

/**
 * @brief Count how many valid numeric coefficients exist in a string
 */
static int count_coefficients(const char* str) {
    if (!str) return 0;
    int count = 0;
    const char* p = str;
    while (*p) {
        while (*p && isspace((unsigned char)*p)) p++;
        if (!*p) break;
        char* endptr;
        strtod(p, &endptr);
        if (endptr == p) break;  // not a number
        p = endptr;
        count++;
    }
    return count;
}

/**
 * @brief Get valid coefficient count range for an MCNP surface mnemonic
 * @param mnemonic  Surface type string
 * @param out_min   Minimum required coefficients (-1 if variable)
 * @param out_max   Maximum accepted coefficients (-1 if variable)
 */
static void get_expected_coeff_range(const char* mnemonic, int* out_min, int* out_max) {
    // min/max = valid coefficient range; -1 = variable (validated by handler)
    static const struct { const char* name; int min; int max; } specs[] = {
        {"px",1,1}, {"py",1,1}, {"pz",1,1}, {"p",-1,-1},
        {"so",1,1}, {"sx",2,2}, {"sy",2,2}, {"sz",2,2}, {"s",4,4},
        {"cx",1,1}, {"cy",1,1}, {"cz",1,1},
        {"c/x",3,3}, {"c/y",3,3}, {"c/z",3,3},
        {"kx",2,3}, {"ky",2,3}, {"kz",2,3},
        {"k/x",4,5}, {"k/y",4,5}, {"k/z",4,5},
        {"rpp",6,6}, {"rcc",7,7}, {"sph",4,4}, {"trc",8,8},
        {"ell",7,7}, {"rec",12,12}, {"box",12,12}, {"wed",12,12},
        {"rhp",-1,-1}, {"hex",-1,-1},
        {"tx",6,6}, {"ty",6,6}, {"tz",6,6},
        {"gq",10,10}, {"sq",7,10},
        {"x",-1,-1}, {"y",-1,-1}, {"z",-1,-1}, {"arb",-1,-1},
        {NULL, 0, 0}
    };
    for (int i = 0; specs[i].name; i++) {
        if (alea_strcasecmp(mnemonic, specs[i].name) == 0) {
            *out_min = specs[i].min;
            *out_max = specs[i].max;
            return;
        }
    }
    *out_min = -1;
    *out_max = -1;
}

/**
 * @brief Validate coefficient count before parsing
 * @return false if wrong number of coefficients (error), true otherwise
 */
static bool validate_coefficient_count(int surface_id, const char* mnemonic,
                                        const char* coefficients) {
    int min_expected, max_expected;
    get_expected_coeff_range(mnemonic, &min_expected, &max_expected);
    if (min_expected < 0) return true;  // variable-count or unknown — skip

    int actual = count_coefficients(coefficients);
    if (actual < min_expected) {
        ALEA_LOG_ERROR("Surface %d (%s): expected %d coefficients, got %d",
                       surface_id, mnemonic, min_expected, actual);
        return false;
    }
    if (actual > max_expected) {
        ALEA_LOG_WARN("Surface %d (%s): expected at most %d coefficients, got %d",
                      surface_id, mnemonic, max_expected, actual);
    }
    return true;
}

// ============================================================================
// AXISYMMETRIC SURFACE HELPER
// ============================================================================

/**
 * @brief Build plane/cylinder/cone/quadric from point pairs
 * 
 * Used by X, Y, Z cards which define surfaces by (axis_coord, radius) pairs.
 * - 1 point  -> plane perpendicular to axis
 * - 2 points -> cylinder (if radii equal) or cone (if radii differ)
 * - 3 points -> general quadric (conic section of revolution)
 */
static bool build_axisymmetric_surface(
    alea_axis_t axis,
    const char* coeffs,
    alea_primitive_type_t* out_type,
    alea_primitive_data_t* out_data
) {
    /* Parse up to 3 (coord, radius) pairs */
    double pts[6] = {0};
    int count = 0;
    const char* cursor = coeffs;

    while (*cursor && count < 6) {
        while (*cursor && isspace(*cursor)) cursor++;
        if (!*cursor) break;

        char* endptr;
        pts[count] = strtod(cursor, &endptr);
        if (endptr == cursor) break;
        cursor = endptr;
        count++;
    }

    int n_points = count / 2;
    if (n_points < 1 || n_points > 3)
        return false;

    double a1 = pts[0], r1 = pts[1];
    double a2 = pts[2], r2 = pts[3];
    double a3 = pts[4], r3 = pts[5];

    /* 1 POINT -> PLANE */
    if (n_points == 1) {
        *out_type = ALEA_PRIMITIVE_PLANE;
        out_data->plane.a = (axis == ALEA_AXIS_X) ? 1.0 : 0.0;
        out_data->plane.b = (axis == ALEA_AXIS_Y) ? 1.0 : 0.0;
        out_data->plane.c = (axis == ALEA_AXIS_Z) ? 1.0 : 0.0;
        out_data->plane.d = -a1;
        return true;
    }

    /* 2 POINTS -> CYLINDER OR CONE */
    if (n_points == 2) {
        if (fabs(r1 - r2) < 1e-10) {
            /* Cylinder */
            if (axis == ALEA_AXIS_X) {
                *out_type = ALEA_PRIMITIVE_CYLINDER_X;
                out_data->cyl_x.center_y = 0.0;
                out_data->cyl_x.center_z = 0.0;
                out_data->cyl_x.radius = r1;
            } else if (axis == ALEA_AXIS_Y) {
                *out_type = ALEA_PRIMITIVE_CYLINDER_Y;
                out_data->cyl_y.center_x = 0.0;
                out_data->cyl_y.center_z = 0.0;
                out_data->cyl_y.radius = r1;
            } else {
                *out_type = ALEA_PRIMITIVE_CYLINDER_Z;
                out_data->cyl_z.center_x = 0.0;
                out_data->cyl_z.center_y = 0.0;
                out_data->cyl_z.radius = r1;
            }
            return true;
        }

        /* Cone */
        double slope = (r2 - r1) / (a2 - a1);
        double apex = a1 - r1 / slope;
        double tan_sq = slope * slope;

        /* Determine sheet selection: both defining points should be on same side of apex */
        int sheet = 0;
        if (a1 > apex && a2 > apex) {
            sheet = 1;   /* Both points on positive side of apex */
        } else if (a1 < apex && a2 < apex) {
            sheet = -1;  /* Both points on negative side of apex */
        }
        /* If points straddle apex (shouldn't happen for valid input), use both sheets */

        if (axis == ALEA_AXIS_X) {
            *out_type = ALEA_PRIMITIVE_CONE_X;
            out_data->cone_x.apex_x = apex;
            out_data->cone_x.apex_y = 0.0;
            out_data->cone_x.apex_z = 0.0;
            out_data->cone_x.tan_angle_sq = tan_sq;
            out_data->cone_x.sheet_selection = sheet;
        } else if (axis == ALEA_AXIS_Y) {
            *out_type = ALEA_PRIMITIVE_CONE_Y;
            out_data->cone_y.apex_x = 0.0;
            out_data->cone_y.apex_y = apex;
            out_data->cone_y.apex_z = 0.0;
            out_data->cone_y.tan_angle_sq = tan_sq;
            out_data->cone_y.sheet_selection = sheet;
        } else {
            *out_type = ALEA_PRIMITIVE_CONE_Z;
            out_data->cone_z.apex_x = 0.0;
            out_data->cone_z.apex_y = 0.0;
            out_data->cone_z.apex_z = apex;
            out_data->cone_z.tan_angle_sq = tan_sq;
            out_data->cone_z.sheet_selection = sheet;
        }
        return true;
    }

    /* 3 POINTS -> QUADRIC */
    double denom = (a1-a2)*(a1-a3)*(a2-a3);
    if (fabs(denom) < 1e-20)
        return false;

    double r1sq = r1*r1, r2sq = r2*r2, r3sq = r3*r3;

    /* Lagrange interpolation: r^2 = A*a^2 + B*a + C */
    double A = (r1sq*(a2-a3) + r2sq*(a3-a1) + r3sq*(a1-a2)) / denom;
    double B = (r1sq*(a3*a3-a2*a2) + r2sq*(a1*a1-a3*a3) + r3sq*(a2*a2-a1*a1)) / denom;
    double C = (r1sq*a2*a3*(a2-a3) + r2sq*a1*a3*(a3-a1) + r3sq*a1*a2*(a1-a2)) / denom;

    /* Convert to quadric: perpendicular^2 - A*axis^2 - B*axis - C = 0 */
    *out_type = ALEA_PRIMITIVE_QUADRIC;
    double* q = out_data->quadric.coeffs;
    q[0] = (axis == ALEA_AXIS_X) ? -A : 1.0;  /* x^2 */
    q[1] = (axis == ALEA_AXIS_Y) ? -A : 1.0;  /* y^2 */
    q[2] = (axis == ALEA_AXIS_Z) ? -A : 1.0;  /* z^2 */
    q[3] = q[4] = q[5] = 0.0;            /* xy, yz, xz */
    q[6] = (axis == ALEA_AXIS_X) ? -B : 0.0;  /* x */
    q[7] = (axis == ALEA_AXIS_Y) ? -B : 0.0;  /* y */
    q[8] = (axis == ALEA_AXIS_Z) ? -B : 0.0;  /* z */
    q[9] = -C;                            /* constant */

    return true;
}

/**
 * @brief Build surface lookup table for fast MCNP surface ID to index mapping
 * 
 * @param sys CSG system
 * @return 0 on success, -1 on error
 */

int alea_build_surface_lookup(alea_system_t* sys) {
    if (!sys) return -1;

    // Find maximum surface ID
    int max_id = 0;
    for (size_t i = 0; i < alea_vec_count(&sys->surfaces); i++) {
        if (sys->surfaces.data[i].mc_surface_id > max_id) {
            max_id = sys->surfaces.data[i].mc_surface_id;
        }
    }

    // Allocate lookup table
    size_t table_size = (size_t)(max_id + 1);
    sys->surface_lookup = malloc(table_size * sizeof(size_t));
    if (!sys->surface_lookup) {
        alea_set_error(ALEA_ERR_OUT_OF_MEMORY, "Failed to allocate surface lookup");
        return -1;
    }
    sys->surface_lookup_size = table_size;

    // Initialize to "not found"
    for (size_t i = 0; i < table_size; i++) {
        sys->surface_lookup[i] = UINT32_MAX;
    }

    // Fill lookup table
    for (size_t i = 0; i < alea_vec_count(&sys->surfaces); i++) {
        int id = sys->surfaces.data[i].mc_surface_id;
        sys->surface_lookup[id] = i;
    }

    ALEA_LOG_INFO("Built surface lookup: %zu entries, max_id=%d",
           alea_vec_count(&sys->surfaces), max_id);

    return 0;
}

/**
 * Get the appropriate sense node for a signed surface reference.
 * 
 * @param sys System
 * @param mc_surface_id Signed surface ID:
 *        +5 → positive sense (outside)
 *        -5 → negative sense (inside)
 * @return Node ID, or ALEA_NODE_ID_INVALID if not found
 */
alea_node_id_t alea_surface_node(const alea_system_t* sys, int mc_surface_id, bool negative) {
    if (!sys || !sys->surface_lookup) {
        return ALEA_NODE_ID_INVALID;
    }
    
    ALEA_LOG_INFO("Looking up surface node for MCNP ID %d (negative=%d)",
           mc_surface_id, negative);
    //int abs_id = negative ? mc_surface_id : -mc_surface_id;
    
    // Bounds check
    if ((size_t)mc_surface_id >= sys->surface_lookup_size) {
        return ALEA_NODE_ID_INVALID;
    }
    
    // Lookup
    size_t idx = sys->surface_lookup[mc_surface_id];
    if (idx == UINT32_MAX) {
        return ALEA_NODE_ID_INVALID;
    }
    const alea_surface_entry_t* surf = &sys->surfaces.data[idx];

    return negative ? surf->neg_node : surf->pos_node;
}

// ============================================================================
// SURFACE CONVERSION 
// ============================================================================

static bool surface_to_primitive(
    const mcnp_surface_t* surf,
    alea_primitive_type_t* out_type,
    alea_primitive_data_t* out_data
) {
    if (!surf || !surf->mnemonic || !out_type || !out_data) {
        return false;
    }
    
    const char* mnemonic = surf->mnemonic;
    const char* coeffs = surf->coefficients_str;

    if (!coeffs) {
        return false;
    }

    // Validate coefficient count before parsing
    if (!validate_coefficient_count(surf->surface_id, mnemonic, coeffs))
        return false;

    const char* cursor = coeffs;
    
    // ========================================================================
    // PLANE SURFACES
    // ========================================================================
    
    if (alea_strcasecmp(mnemonic, "px") == 0) {
        double d = parse_coefficient(&cursor);
        *out_type = ALEA_PRIMITIVE_PLANE;
        out_data->plane.a = 1.0;
        out_data->plane.b = 0.0;
        out_data->plane.c = 0.0;
        out_data->plane.d = -d;
        
    } else if (alea_strcasecmp(mnemonic, "py") == 0) {
        double d = parse_coefficient(&cursor);
        *out_type = ALEA_PRIMITIVE_PLANE;
        out_data->plane.a = 0.0;
        out_data->plane.b = 1.0;
        out_data->plane.c = 0.0;
        out_data->plane.d = -d;
        
    } else if (alea_strcasecmp(mnemonic, "pz") == 0) {
        double d = parse_coefficient(&cursor);
        *out_type = ALEA_PRIMITIVE_PLANE;
        out_data->plane.a = 0.0;
        out_data->plane.b = 0.0;
        out_data->plane.c = 1.0;
        out_data->plane.d = -d;
        
    } else if (alea_strcasecmp(mnemonic, "p") == 0) {
        *out_type = ALEA_PRIMITIVE_PLANE;
        int ncoeffs = count_coefficients(coeffs);
        if (ncoeffs == 4) {
            // Standard form: A B C D
            out_data->plane.a = parse_coefficient(&cursor);
            out_data->plane.b = parse_coefficient(&cursor);
            out_data->plane.c = parse_coefficient(&cursor);
            out_data->plane.d = -parse_coefficient(&cursor);
        } else if (ncoeffs == 9) {
            // Three-point form: x1 y1 z1 x2 y2 z2 x3 y3 z3
            double x1 = parse_coefficient(&cursor), y1 = parse_coefficient(&cursor), z1 = parse_coefficient(&cursor);
            double x2 = parse_coefficient(&cursor), y2 = parse_coefficient(&cursor), z2 = parse_coefficient(&cursor);
            double x3 = parse_coefficient(&cursor), y3 = parse_coefficient(&cursor), z3 = parse_coefficient(&cursor);
            // Normal = (P2-P1) x (P3-P1)
            double ux = x2 - x1, uy = y2 - y1, uz = z2 - z1;
            double vx = x3 - x1, vy = y3 - y1, vz = z3 - z1;
            double a = uy * vz - uz * vy;
            double b = uz * vx - ux * vz;
            double c = ux * vy - uy * vx;
            double norm = sqrt(a*a + b*b + c*c);
            if (norm < 1e-20) {
                ALEA_LOG_ERROR("Surface %d (P): three points are collinear", surf->surface_id);
                return false;
            }
            a /= norm; b /= norm; c /= norm;
            double d = a * x1 + b * y1 + c * z1;
            out_data->plane.a = a;
            out_data->plane.b = b;
            out_data->plane.c = c;
            out_data->plane.d = -d;
        } else {
            ALEA_LOG_ERROR("Surface %d (P): expected 4 or 9 coefficients, got %d", surf->surface_id, ncoeffs);
            return false;
        }
    
    // ========================================================================
    // SPHERE SURFACES
    // ========================================================================
    
    } else if (alea_strcasecmp(mnemonic, "so") == 0) {
        double r = parse_coefficient(&cursor);
        *out_type = ALEA_PRIMITIVE_SPHERE;
        out_data->sphere.center_x = 0.0;
        out_data->sphere.center_y = 0.0;
        out_data->sphere.center_z = 0.0;
        out_data->sphere.radius = r;

    } else if (alea_strcasecmp(mnemonic, "sx") == 0) {
        *out_type = ALEA_PRIMITIVE_SPHERE;
        out_data->sphere.center_x = parse_coefficient(&cursor);
        out_data->sphere.center_y = 0.0;
        out_data->sphere.center_z = 0.0;
        out_data->sphere.radius = parse_coefficient(&cursor);
    
    } else if (alea_strcasecmp(mnemonic, "sy") == 0) {
        *out_type = ALEA_PRIMITIVE_SPHERE;
        out_data->sphere.center_x = 0.0;
        out_data->sphere.center_y = parse_coefficient(&cursor);
        out_data->sphere.center_z = 0.0;
        out_data->sphere.radius = parse_coefficient(&cursor);
    
    } else if (alea_strcasecmp(mnemonic, "sz") == 0) {
        *out_type = ALEA_PRIMITIVE_SPHERE;
        out_data->sphere.center_x = 0.0;
        out_data->sphere.center_y = 0.0;
        out_data->sphere.center_z = parse_coefficient(&cursor);
        out_data->sphere.radius = parse_coefficient(&cursor);
        
    } else if (alea_strcasecmp(mnemonic, "s") == 0) {
        *out_type = ALEA_PRIMITIVE_SPHERE;
        out_data->sphere.center_x = parse_coefficient(&cursor);
        out_data->sphere.center_y = parse_coefficient(&cursor);
        out_data->sphere.center_z = parse_coefficient(&cursor);
        out_data->sphere.radius = parse_coefficient(&cursor);
    
    // ========================================================================
    // CYLINDER SURFACES (axis-aligned)
    // ========================================================================
    
    } else if (alea_strcasecmp(mnemonic, "cx") == 0) {
        *out_type = ALEA_PRIMITIVE_CYLINDER_X;
        out_data->cyl_x.center_y = 0.0;
        out_data->cyl_x.center_z = 0.0;
        out_data->cyl_x.radius = parse_coefficient(&cursor);
        
    } else if (alea_strcasecmp(mnemonic, "cy") == 0) {
        *out_type = ALEA_PRIMITIVE_CYLINDER_Y;
        out_data->cyl_y.center_x = 0.0;
        out_data->cyl_y.center_z = 0.0;
        out_data->cyl_y.radius = parse_coefficient(&cursor);
        
    } else if (alea_strcasecmp(mnemonic, "cz") == 0) {
        *out_type = ALEA_PRIMITIVE_CYLINDER_Z;
        out_data->cyl_z.center_x = 0.0;
        out_data->cyl_z.center_y = 0.0;
        out_data->cyl_z.radius = parse_coefficient(&cursor);
    
    // ========================================================================
    // CYLINDER SURFACES (off-axis)
    // ========================================================================
    
    } else if (alea_strcasecmp(mnemonic, "c/x") == 0) {
        *out_type = ALEA_PRIMITIVE_CYLINDER_X;
        out_data->cyl_x.center_y = parse_coefficient(&cursor);
        out_data->cyl_x.center_z = parse_coefficient(&cursor);
        out_data->cyl_x.radius = parse_coefficient(&cursor);
        
    } else if (alea_strcasecmp(mnemonic, "c/y") == 0) {
        *out_type = ALEA_PRIMITIVE_CYLINDER_Y;
        out_data->cyl_y.center_x = parse_coefficient(&cursor);
        out_data->cyl_y.center_z = parse_coefficient(&cursor);
        out_data->cyl_y.radius = parse_coefficient(&cursor);
        
    } else if (alea_strcasecmp(mnemonic, "c/z") == 0) {
        *out_type = ALEA_PRIMITIVE_CYLINDER_Z;
        out_data->cyl_z.center_x = parse_coefficient(&cursor);
        out_data->cyl_z.center_y = parse_coefficient(&cursor);
        out_data->cyl_z.radius = parse_coefficient(&cursor);
    
    // ========================================================================
    // CONE SURFACES (on-axis)
    // ========================================================================
    
    } else if (alea_strcasecmp(mnemonic, "kx") == 0) {
        *out_type = ALEA_PRIMITIVE_CONE_X;
        out_data->cone_x.apex_x = parse_coefficient(&cursor);
        out_data->cone_x.apex_y = 0.0;
        out_data->cone_x.apex_z = 0.0;
        double t_sq = parse_coefficient(&cursor);
        if (t_sq < 0) {
            ALEA_LOG_WARN("Negative t² in KX surface %d", surf->surface_id);
            return false;
        }
        out_data->cone_x.tan_angle_sq = t_sq;
        out_data->cone_x.sheet_selection = parse_coefficient(&cursor);
        
    } else if (alea_strcasecmp(mnemonic, "ky") == 0) {
        *out_type = ALEA_PRIMITIVE_CONE_Y;
        out_data->cone_y.apex_x = 0.0;
        out_data->cone_y.apex_y = parse_coefficient(&cursor);
        out_data->cone_y.apex_z = 0.0;
        double t_sq = parse_coefficient(&cursor);
        if (t_sq < 0) {
            ALEA_LOG_WARN("Negative t² in KY surface %d", surf->surface_id);
            return false;
        }
        out_data->cone_y.tan_angle_sq = t_sq;
        out_data->cone_y.sheet_selection = parse_coefficient(&cursor);
        
    } else if (alea_strcasecmp(mnemonic, "kz") == 0) {
        *out_type = ALEA_PRIMITIVE_CONE_Z;
        out_data->cone_z.apex_x = 0.0;
        out_data->cone_z.apex_y = 0.0;
        out_data->cone_z.apex_z = parse_coefficient(&cursor);
        double t_sq = parse_coefficient(&cursor);
        if (t_sq < 0) {
            ALEA_LOG_WARN("Negative t² in KZ surface %d", surf->surface_id);
            return false;
        }
        out_data->cone_z.tan_angle_sq = t_sq;
        out_data->cone_z.sheet_selection = parse_coefficient(&cursor);
    
    // ========================================================================
    // CONE SURFACES (off-axis)
    // ========================================================================
    
    } else if (alea_strcasecmp(mnemonic, "k/x") == 0) {
        *out_type = ALEA_PRIMITIVE_CONE_X;
        out_data->cone_x.apex_x = parse_coefficient(&cursor);
        out_data->cone_x.apex_y = parse_coefficient(&cursor);
        out_data->cone_x.apex_z = parse_coefficient(&cursor);
        double t_sq = parse_coefficient(&cursor);
        if (t_sq < 0) {
            ALEA_LOG_WARN("Negative t² in K/X surface %d", surf->surface_id);
            return false;
        }
        out_data->cone_x.tan_angle_sq = t_sq;
        out_data->cone_x.sheet_selection = parse_coefficient(&cursor);
        
    } else if (alea_strcasecmp(mnemonic, "k/y") == 0) {
        *out_type = ALEA_PRIMITIVE_CONE_Y;
        out_data->cone_y.apex_x = parse_coefficient(&cursor);
        out_data->cone_y.apex_y = parse_coefficient(&cursor);
        out_data->cone_y.apex_z = parse_coefficient(&cursor);
        double t_sq = parse_coefficient(&cursor);
        if (t_sq < 0) {
            ALEA_LOG_WARN("Negative t² in K/Y surface %d", surf->surface_id);
            return false;
        }
        out_data->cone_y.tan_angle_sq = t_sq;
        out_data->cone_y.sheet_selection = parse_coefficient(&cursor);
        
    } else if (alea_strcasecmp(mnemonic, "k/z") == 0) {
        *out_type = ALEA_PRIMITIVE_CONE_Z;
        out_data->cone_z.apex_x = parse_coefficient(&cursor);
        out_data->cone_z.apex_y = parse_coefficient(&cursor);
        out_data->cone_z.apex_z = parse_coefficient(&cursor);
        double t_sq = parse_coefficient(&cursor);
        if (t_sq < 0) {
            ALEA_LOG_WARN("Negative t² in K/Z surface %d", surf->surface_id);
            return false;
        }
        out_data->cone_z.tan_angle_sq = t_sq;
        out_data->cone_z.sheet_selection = parse_coefficient(&cursor);
    
    // ========================================================================
    // BOX SURFACE (RPP)
    // ========================================================================
    
    } else if (alea_strcasecmp(mnemonic, "rpp") == 0) {
        *out_type = ALEA_PRIMITIVE_RPP;
        out_data->box.min_x = parse_coefficient(&cursor);
        out_data->box.max_x = parse_coefficient(&cursor);
        out_data->box.min_y = parse_coefficient(&cursor);
        out_data->box.max_y = parse_coefficient(&cursor);
        out_data->box.min_z = parse_coefficient(&cursor);
        out_data->box.max_z = parse_coefficient(&cursor);
    
    // ========================================================================
    // RIGHT CIRCULAR CYLINDER (RCC)
    // ========================================================================
    
    } else if (alea_strcasecmp(mnemonic, "rcc") == 0) {
        *out_type = ALEA_PRIMITIVE_RCC;
        out_data->rcc.base_x = parse_coefficient(&cursor);
        out_data->rcc.base_y = parse_coefficient(&cursor);
        out_data->rcc.base_z = parse_coefficient(&cursor);
        out_data->rcc.height_x = parse_coefficient(&cursor);
        out_data->rcc.height_y = parse_coefficient(&cursor);
        out_data->rcc.height_z = parse_coefficient(&cursor);
        out_data->rcc.radius = parse_coefficient(&cursor);

    // ========================================================================
    // SPHERE MACROBODY (SPH)
    // ========================================================================

    } else if (alea_strcasecmp(mnemonic, "sph") == 0) {
        *out_type = ALEA_PRIMITIVE_SPH;
        out_data->sph.center_x = parse_coefficient(&cursor);
        out_data->sph.center_y = parse_coefficient(&cursor);
        out_data->sph.center_z = parse_coefficient(&cursor);
        out_data->sph.radius = parse_coefficient(&cursor);

    // ========================================================================
    // TRUNCATED RIGHT CONE (TRC)
    // ========================================================================

    } else if (alea_strcasecmp(mnemonic, "trc") == 0) {
        *out_type = ALEA_PRIMITIVE_TRC;
        out_data->trc.base_x = parse_coefficient(&cursor);
        out_data->trc.base_y = parse_coefficient(&cursor);
        out_data->trc.base_z = parse_coefficient(&cursor);
        out_data->trc.height_x = parse_coefficient(&cursor);
        out_data->trc.height_y = parse_coefficient(&cursor);
        out_data->trc.height_z = parse_coefficient(&cursor);
        out_data->trc.base_radius = parse_coefficient(&cursor);
        out_data->trc.top_radius = parse_coefficient(&cursor);

    // ========================================================================
    // ELLIPSOID (ELL)
    // ========================================================================

    } else if (alea_strcasecmp(mnemonic, "ell") == 0) {
        *out_type = ALEA_PRIMITIVE_ELL;
        out_data->ell.v1_x = parse_coefficient(&cursor);
        out_data->ell.v1_y = parse_coefficient(&cursor);
        out_data->ell.v1_z = parse_coefficient(&cursor);
        out_data->ell.v2_x = parse_coefficient(&cursor);
        out_data->ell.v2_y = parse_coefficient(&cursor);
        out_data->ell.v2_z = parse_coefficient(&cursor);
        out_data->ell.major_axis_len = parse_coefficient(&cursor);

    // ========================================================================
    // RIGHT ELLIPTICAL CYLINDER (REC)
    // ========================================================================

    } else if (alea_strcasecmp(mnemonic, "rec") == 0) {
        *out_type = ALEA_PRIMITIVE_REC;
        out_data->rec.base_x = parse_coefficient(&cursor);
        out_data->rec.base_y = parse_coefficient(&cursor);
        out_data->rec.base_z = parse_coefficient(&cursor);
        out_data->rec.height_x = parse_coefficient(&cursor);
        out_data->rec.height_y = parse_coefficient(&cursor);
        out_data->rec.height_z = parse_coefficient(&cursor);
        out_data->rec.axis1_x = parse_coefficient(&cursor);
        out_data->rec.axis1_y = parse_coefficient(&cursor);
        out_data->rec.axis1_z = parse_coefficient(&cursor);
        out_data->rec.axis2_x = parse_coefficient(&cursor);
        out_data->rec.axis2_y = parse_coefficient(&cursor);
        out_data->rec.axis2_z = parse_coefficient(&cursor);

    // ========================================================================
    // GENERAL ORIENTED BOX (BOX)
    // ========================================================================

    } else if (alea_strcasecmp(mnemonic, "box") == 0) {
        *out_type = ALEA_PRIMITIVE_BOX;
        out_data->box_general.corner_x = parse_coefficient(&cursor);
        out_data->box_general.corner_y = parse_coefficient(&cursor);
        out_data->box_general.corner_z = parse_coefficient(&cursor);
        out_data->box_general.v1_x = parse_coefficient(&cursor);
        out_data->box_general.v1_y = parse_coefficient(&cursor);
        out_data->box_general.v1_z = parse_coefficient(&cursor);
        out_data->box_general.v2_x = parse_coefficient(&cursor);
        out_data->box_general.v2_y = parse_coefficient(&cursor);
        out_data->box_general.v2_z = parse_coefficient(&cursor);
        out_data->box_general.v3_x = parse_coefficient(&cursor);
        out_data->box_general.v3_y = parse_coefficient(&cursor);
        out_data->box_general.v3_z = parse_coefficient(&cursor);

    // ========================================================================
    // WEDGE (WED)
    // ========================================================================

    } else if (alea_strcasecmp(mnemonic, "wed") == 0) {
        *out_type = ALEA_PRIMITIVE_WED;
        out_data->wed.vertex_x = parse_coefficient(&cursor);
        out_data->wed.vertex_y = parse_coefficient(&cursor);
        out_data->wed.vertex_z = parse_coefficient(&cursor);
        out_data->wed.v1_x = parse_coefficient(&cursor);
        out_data->wed.v1_y = parse_coefficient(&cursor);
        out_data->wed.v1_z = parse_coefficient(&cursor);
        out_data->wed.v2_x = parse_coefficient(&cursor);
        out_data->wed.v2_y = parse_coefficient(&cursor);
        out_data->wed.v2_z = parse_coefficient(&cursor);
        out_data->wed.v3_x = parse_coefficient(&cursor);
        out_data->wed.v3_y = parse_coefficient(&cursor);
        out_data->wed.v3_z = parse_coefficient(&cursor);

    // ========================================================================
    // RIGHT HEXAGONAL PRISM (RHP/HEX)
    // ========================================================================

    } else if (alea_strcasecmp(mnemonic, "rhp") == 0 || alea_strcasecmp(mnemonic, "hex") == 0) {
        *out_type = ALEA_PRIMITIVE_RHP;
        out_data->rhp.base_x = parse_coefficient(&cursor);
        out_data->rhp.base_y = parse_coefficient(&cursor);
        out_data->rhp.base_z = parse_coefficient(&cursor);
        out_data->rhp.height_x = parse_coefficient(&cursor);
        out_data->rhp.height_y = parse_coefficient(&cursor);
        out_data->rhp.height_z = parse_coefficient(&cursor);
        out_data->rhp.r1_x = parse_coefficient(&cursor);
        out_data->rhp.r1_y = parse_coefficient(&cursor);
        out_data->rhp.r1_z = parse_coefficient(&cursor);
        out_data->rhp.r2_x = parse_coefficient(&cursor);
        out_data->rhp.r2_y = parse_coefficient(&cursor);
        out_data->rhp.r2_z = parse_coefficient(&cursor);
        out_data->rhp.r3_x = parse_coefficient(&cursor);
        out_data->rhp.r3_y = parse_coefficient(&cursor);
        out_data->rhp.r3_z = parse_coefficient(&cursor);

    // ========================================================================
    // ARBITRARY POLYHEDRON (ARB)
    // ========================================================================

    } else if (alea_strcasecmp(mnemonic, "arb") == 0) {
        *out_type = ALEA_PRIMITIVE_ARB;
        /* ARB format: 8 corners (x y z) then face definitions as 4-digit numbers */
        /* Parse all 8 corners */
        for (int i = 0; i < 8; i++) {
            out_data->arb.corners[i][0] = parse_coefficient(&cursor);
            out_data->arb.corners[i][1] = parse_coefficient(&cursor);
            out_data->arb.corners[i][2] = parse_coefficient(&cursor);
        }
        out_data->arb.num_corners = 8;

        /* Parse face definitions - each is a 4-digit number where digits are vertex indices */
        out_data->arb.num_faces = 0;
        for (int f = 0; f < 6; f++) {
            double face_val = parse_coefficient(&cursor);
            if (face_val == 0.0) break;

            int face_int = (int)face_val;
            out_data->arb.faces[f][0] = (face_int / 1000) % 10;
            out_data->arb.faces[f][1] = (face_int / 100) % 10;
            out_data->arb.faces[f][2] = (face_int / 10) % 10;
            out_data->arb.faces[f][3] = face_int % 10;
            out_data->arb.num_faces++;
        }

    // ========================================================================
    // TORUS SURFACES
    // ========================================================================
    
    } else if (alea_strcasecmp(mnemonic, "tx") == 0) {
        double cx = parse_coefficient(&cursor);
        double cy = parse_coefficient(&cursor);
        double cz = parse_coefficient(&cursor);
        double major_radius = parse_coefficient(&cursor);
        double axial_semiwidth = parse_coefficient(&cursor);
        double minor_radius = parse_coefficient(&cursor);
        if (torus_is_sphere(major_radius, axial_semiwidth, minor_radius)) {
            set_sphere_from_torus(out_type, out_data, cx, cy, cz, axial_semiwidth);
        } else {
            *out_type = ALEA_PRIMITIVE_TORUS_X;
            out_data->torus.axis = ALEA_AXIS_X;
            out_data->torus.center_x = cx;
            out_data->torus.center_y = cy;
            out_data->torus.center_z = cz;
            out_data->torus.major_radius = major_radius;
            out_data->torus.axial_semiwidth_B = axial_semiwidth;
            out_data->torus.minor_radius = minor_radius;
        }
    
    } else if (alea_strcasecmp(mnemonic, "ty") == 0) {
        double cx = parse_coefficient(&cursor);
        double cy = parse_coefficient(&cursor);
        double cz = parse_coefficient(&cursor);
        double major_radius = parse_coefficient(&cursor);
        double axial_semiwidth = parse_coefficient(&cursor);
        double minor_radius = parse_coefficient(&cursor);
        if (torus_is_sphere(major_radius, axial_semiwidth, minor_radius)) {
            set_sphere_from_torus(out_type, out_data, cx, cy, cz, axial_semiwidth);
        } else {
            *out_type = ALEA_PRIMITIVE_TORUS_Y;
            out_data->torus.axis = ALEA_AXIS_Y;
            out_data->torus.center_x = cx;
            out_data->torus.center_y = cy;
            out_data->torus.center_z = cz;
            out_data->torus.major_radius = major_radius;
            out_data->torus.axial_semiwidth_B = axial_semiwidth;
            out_data->torus.minor_radius = minor_radius;
        }
    
    } else if (alea_strcasecmp(mnemonic, "tz") == 0) {
        double cx = parse_coefficient(&cursor);
        double cy = parse_coefficient(&cursor);
        double cz = parse_coefficient(&cursor);
        double major_radius = parse_coefficient(&cursor);
        double axial_semiwidth = parse_coefficient(&cursor);
        double minor_radius = parse_coefficient(&cursor);
        if (torus_is_sphere(major_radius, axial_semiwidth, minor_radius)) {
            set_sphere_from_torus(out_type, out_data, cx, cy, cz, axial_semiwidth);
        } else {
            *out_type = ALEA_PRIMITIVE_TORUS_Z;
            out_data->torus.axis = ALEA_AXIS_Z;
            out_data->torus.center_x = cx;
            out_data->torus.center_y = cy;
            out_data->torus.center_z = cz;
            out_data->torus.major_radius = major_radius;
            out_data->torus.axial_semiwidth_B = axial_semiwidth;
            out_data->torus.minor_radius = minor_radius;
        }
    
    // ========================================================================
    // QUADRIC SURFACES
    // ========================================================================
    
    } else if (alea_strcasecmp(mnemonic, "gq") == 0) {
        *out_type = ALEA_PRIMITIVE_QUADRIC;
        for (int i = 0; i < 10; i++) {
            out_data->quadric.coeffs[i] = parse_coefficient(&cursor);
        }
    
    } else if (alea_strcasecmp(mnemonic, "sq") == 0) {
        *out_type = ALEA_PRIMITIVE_QUADRIC;
        // Special quadric - convert to general quadric
        // MCNP SQ: A(x-x₀)² + B(y-y₀)² + C(z-z₀)² + 2D(x-x₀) + 2E(y-y₀) + 2F(z-z₀) + G = 0
        double a = parse_coefficient(&cursor);
        double b = parse_coefficient(&cursor);
        double c = parse_coefficient(&cursor);
        double d = parse_coefficient(&cursor);
        double e = parse_coefficient(&cursor);
        double f = parse_coefficient(&cursor);
        double g = parse_coefficient(&cursor);
        double x0 = parse_coefficient(&cursor);
        double y0 = parse_coefficient(&cursor);
        double z0 = parse_coefficient(&cursor);
        
        // Convert to general quadric form: Ax² + By² + Cz² + Dxy + Eyz + Fxz + Gx + Hy + Iz + J = 0
        out_data->quadric.coeffs[0] = a;  // x² coefficient
        out_data->quadric.coeffs[1] = b;  // y² coefficient  
        out_data->quadric.coeffs[2] = c;  // z² coefficient
        out_data->quadric.coeffs[3] = 0;  // xy coefficient
        out_data->quadric.coeffs[4] = 0;  // yz coefficient
        out_data->quadric.coeffs[5] = 0;  // xz coefficient
        out_data->quadric.coeffs[6] = 2*d - 2*a*x0;  // x coefficient
        out_data->quadric.coeffs[7] = 2*e - 2*b*y0;  // y coefficient
        out_data->quadric.coeffs[8] = 2*f - 2*c*z0;  // z coefficient
        out_data->quadric.coeffs[9] = g + a*x0*x0 + b*y0*y0 + c*z0*z0 - 2*d*x0 - 2*e*y0 - 2*f*z0;  // constant
    
    // ========================================================================
    // AXISYMMETRIC SURFACES DEFINED BY POINTS (X, Y, Z cards)
    // ========================================================================
    
    } else if (alea_strcasecmp(mnemonic, "x") == 0) {
        if (!build_axisymmetric_surface(ALEA_AXIS_X, coeffs, out_type, out_data))
            return false;
            
    } else if (alea_strcasecmp(mnemonic, "y") == 0) {
        if (!build_axisymmetric_surface(ALEA_AXIS_Y, coeffs, out_type, out_data))
            return false;
            
    } else if (alea_strcasecmp(mnemonic, "z") == 0) {
        if (!build_axisymmetric_surface(ALEA_AXIS_Z, coeffs, out_type, out_data))
            return false;

    } else {
        ALEA_LOG_WARN("Unsupported surface type '%s' for surface %d",
               mnemonic, surf->surface_id);
        return false;
    }
    
    return true;
}


int alea_convert_surface(alea_system_t* sys, const mcnp_surface_t* surf) {
    if (!sys || !surf) {
        alea_set_error(ALEA_ERR_NULL_ARG, "NULL argument");
        return -1;
    }
    
    // Parse MCNP surface into primitive (in local/auxiliary coordinates)
    alea_primitive_type_t local_type;
    alea_primitive_data_t local_data;
    
    if (!surface_to_primitive(surf, &local_type, &local_data)) {
        sys->stats.failed_surfaces++;
        return -1;
    }

    // Apply transform if present (convert to global coordinates)
    alea_primitive_type_t global_type = local_type;
    alea_primitive_data_t global_data = local_data;
    bool transform_applied = false;

    if (surf->transform_id != 0) {
        
        const alea_transform_t* tr = alea_get_transform(sys, surf->transform_id);
        
        if (!tr) {
            sys->stats.failed_surfaces++;
            alea_set_error_detail(ALEA_ERR_INVALID_ID,
                    "Surface %d references unknown transform TR%d",
                    surf->surface_id, surf->transform_id);
            return -1;
        }

        if (!alea_apply_transform_to_primitive(tr, local_type, &local_data,
                                               &global_type, &global_data)) {
            sys->stats.failed_surfaces++;
            alea_set_error_detail(ALEA_ERR_UNSUPPORTED,
                    "Surface %d uses transform TR%d, but primitive type %d cannot be transformed exactly",
                    surf->surface_id, surf->transform_id, local_type);
            return -1;
        }

        transform_applied = true;
        if (global_type != local_type) {
            ALEA_LOG_INFO("Surface %d: transformed from type %d to type %d",
                   surf->surface_id, local_type, global_type);
        }
    }
    
    alea_node_id_t pos_node, neg_node;
    alea_node_id_t expanded_pos_node = ALEA_NODE_ID_INVALID;
    alea_node_id_t expanded_neg_node = ALEA_NODE_ID_INVALID;
    alea_primitive_id_t prim_id;

    // Get or create primitive (with automatic deduplication)
    int8_t inverted = 0;
    prim_id = alea_get_or_create_primitive(sys, global_type, &global_data, &inverted);
    if (prim_id == UINT32_MAX) {
        sys->stats.failed_surfaces++;
        return -1;
    }

    // Log dedup: ref_count > 1 means this primitive already existed
    if (sys->primitives.data[prim_id].ref_count > 1) {
        ALEA_LOG_INFO("Surface %d: deduplicated to primitive %u (ref_count=%u)",
                     surf->surface_id, prim_id,
                     sys->primitives.data[prim_id].ref_count);
    }

    // === Create POSITIVE and NEGATIVE sense nodes ===
    pos_node = alea_add_primitive_node(sys, prim_id, +1, inverted, surf->surface_id);
    if (pos_node == ALEA_NODE_ID_INVALID) {
        return -1;
    }

    neg_node = alea_add_primitive_node(sys, prim_id, -1, inverted, surf->surface_id);
    if (neg_node == ALEA_NODE_ID_INVALID) {
        return -1;
    }

    // === For macrobodies: also create expanded representation ===
    if (alea_is_macrobody(global_type)) {
        if (alea_expand_macrobody_immediate(sys, global_type, &global_data,
                                            &expanded_neg_node, &expanded_pos_node) == 0) {
            ALEA_LOG_DEBUG("Surface %d: macrobody with expanded form (exp_neg=%u, exp_pos=%u)",
                          surf->surface_id, expanded_neg_node, expanded_pos_node);
        } else {
            ALEA_LOG_WARN("Macrobody expansion failed for surface %d", surf->surface_id);
            // Keep expanded nodes as INVALID - expansion will happen at export if needed
        }
    }

    // === For 1-sheet cones: also create expanded representation ===
    // Check for sheet_selection != 0 which indicates a 1-sheet cone
    bool is_1sheet = false;
    switch (global_type) {
        case ALEA_PRIMITIVE_CONE_X:
            is_1sheet = (global_data.cone_x.sheet_selection != 0);
            break;
        case ALEA_PRIMITIVE_CONE_Y:
            is_1sheet = (global_data.cone_y.sheet_selection != 0);
            break;
        case ALEA_PRIMITIVE_CONE_Z:
            is_1sheet = (global_data.cone_z.sheet_selection != 0);
            break;
        default:
            break;
    }
    if (is_1sheet) {
        if (alea_expand_1sheet_cone_immediate(sys, global_type, &global_data,
                                              &expanded_neg_node, &expanded_pos_node) == 0) {
            ALEA_LOG_DEBUG("Surface %d: 1-sheet cone with expanded form (exp_neg=%u, exp_pos=%u)",
                          surf->surface_id, expanded_neg_node, expanded_pos_node);
        } else {
            ALEA_LOG_WARN("1-sheet cone expansion failed for surface %d", surf->surface_id);
        }
    }

    // === Map MCNP boundary type to CSG boundary type ===
    alea_boundary_type_t boundary = ALEA_BOUNDARY_TRANSMISSIVE;
    switch (surf->boundary_type) {
        case SURFACE_REFLECTIVE:
            boundary = ALEA_BOUNDARY_REFLECTIVE;
            break;
        case SURFACE_WHITE:
            boundary = ALEA_BOUNDARY_WHITE;
            break;
        case SURFACE_PERIODIC:
            boundary = ALEA_BOUNDARY_PERIODIC;
            break;
        case SURFACE_VACUUM:
            boundary = ALEA_BOUNDARY_VACUUM;
            break;
        default:
            boundary = ALEA_BOUNDARY_TRANSMISSIVE;
            break;
    }
    
    // === Store surface entry using vector API ===
    alea_surface_entry_t* surf_entry = alea_vec_push_uninit(&sys->surfaces, alea_surface_entry_t);
    if (!surf_entry) {
        return -1;
    }
    *surf_entry = (alea_surface_entry_t){
        .mc_surface_id = surf->surface_id,
        .primitive_id = prim_id,
        .pos_node = pos_node,
        .neg_node = neg_node,
        .transform_id = surf->transform_id,  // Keep original for roundtrip export
        .transform_applied = transform_applied,
        .boundary_type = boundary,
        .periodic_surface_id = surf->periodic_surface_id,
        .expanded_pos_node = expanded_pos_node,
        .expanded_neg_node = expanded_neg_node
    };

    sys->stats.surfaces_converted++;
    alea_system_invalidate_query_caches(sys, ALEA_CACHE_ALL);

    ALEA_LOG_DEBUG("Surface %d: prim=%u (type=%d), pos=%u, neg=%u, boundary=%d%s",
                  surf->surface_id, prim_id, global_type, pos_node, neg_node, boundary,
                  transform_applied ? " (transformed)" : "");
    
    return 0;
}



// ============================================================================
// LOOKUP FUNCTIONS
// ============================================================================

uint32_t alea_find_surface_primitive(const alea_system_t* sys, int mc_surface_id) {
    if (!sys) return UINT32_MAX;

    int abs_id = abs(mc_surface_id);

    for (size_t i = 0; i < alea_vec_count(&sys->surfaces); i++) {
        if (sys->surfaces.data[i].mc_surface_id == abs_id) {
            return sys->surfaces.data[i].primitive_id;
        }
    }

    return UINT32_MAX;
}

bool alea_surface_sense(int mc_surface_id) {
    // In MCNP: positive ID means outside  (positive sense)
    //          negative ID means inside (negative sense)
    return mc_surface_id > 0;
}
