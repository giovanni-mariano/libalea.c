// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/**
 * @file alea_macrobody.c
 * @brief Macrobody expansion - convert macrobodies to primitive surface CSG trees
 */

#include "alea_macrobody.h"
#include "alea_simplify.h"
#include "primitives/bbox.h"
#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

/* ============================================================================
 * EXPANSION CACHE
 * ============================================================================ */

/**
 * @brief Cache entry for expanded macrobody/1-sheet cone
 *
 * Stores both positive and negative sense representations of an expanded
 * primitive. This allows reuse when the same macrobody is referenced
 * multiple times in the geometry.
 */
typedef struct {
    alea_primitive_id_t prim_id;     /* Original primitive ID */
    alea_node_id_t neg_node;         /* Interior (negative sense) representation */
    alea_node_id_t pos_node;         /* Exterior (positive sense) representation */
} expansion_cache_entry_t;

/**
 * @brief Expansion cache - maps primitive IDs to expanded node pairs
 */
typedef struct {
    expansion_cache_entry_t* entries;
    size_t count;
    size_t capacity;
} expansion_cache_t;

static void cache_init(expansion_cache_t* cache) {
    cache->entries = NULL;
    cache->count = 0;
    cache->capacity = 0;
}

static void cache_free(expansion_cache_t* cache) {
    free(cache->entries);
    cache->entries = NULL;
    cache->count = 0;
    cache->capacity = 0;
}

static expansion_cache_entry_t* cache_find(expansion_cache_t* cache, alea_primitive_id_t prim_id) {
    for (size_t i = 0; i < cache->count; i++) {
        if (cache->entries[i].prim_id == prim_id) {
            return &cache->entries[i];
        }
    }
    return NULL;
}

static expansion_cache_entry_t* cache_add(expansion_cache_t* cache, alea_primitive_id_t prim_id,
                                           alea_node_id_t neg_node, alea_node_id_t pos_node) {
    if (cache->count >= cache->capacity) {
        size_t new_cap = cache->capacity ? cache->capacity * 2 : 16;
        expansion_cache_entry_t* new_entries = realloc(cache->entries,
                                                        new_cap * sizeof(expansion_cache_entry_t));
        if (!new_entries) return NULL;
        cache->entries = new_entries;
        cache->capacity = new_cap;
    }

    expansion_cache_entry_t* entry = &cache->entries[cache->count++];
    entry->prim_id = prim_id;
    entry->neg_node = neg_node;
    entry->pos_node = pos_node;
    return entry;
}

/* ============================================================================
 * HELPERS
 * ============================================================================ */

static double vec3_len(double x, double y, double z) {
    return sqrt(x*x + y*y + z*z);
}

/**
 * @brief Create a primitive node WITHOUT assigning MCNP surface ID
 *
 * Used during macrobody expansion at creation time. Surface IDs are
 * assigned later at export.
 */
static alea_node_id_t create_primitive_node_internal(alea_system_t* sys,
                                                     alea_primitive_type_t type,
                                                     alea_primitive_data_t* data,
                                                     int8_t sense) {
    int8_t inverted = 0;
    alea_primitive_id_t prim_id = alea_get_or_create_primitive(sys, type, data, &inverted);
    if (prim_id == UINT32_MAX) return ALEA_NODE_ID_INVALID;

    return alea_add_primitive_node(sys, prim_id, sense, inverted, 0);
}

static void vec3_normalize(double* x, double* y, double* z) {
    double len = vec3_len(*x, *y, *z);
    if (len > 1e-20) {
        *x /= len; *y /= len; *z /= len;
    }
}

static void vec3_cross(double ax, double ay, double az,
                       double bx, double by, double bz,
                       double* rx, double* ry, double* rz) {
    *rx = ay*bz - az*by;
    *ry = az*bx - ax*bz;
    *rz = ax*by - ay*bx;
}


/* Create an intersection of two nodes */
static alea_node_id_t create_intersection(alea_system_t* sys,
                                          alea_node_id_t left,
                                          alea_node_id_t right) {
    if (left == ALEA_NODE_ID_INVALID) return right;
    if (right == ALEA_NODE_ID_INVALID) return left;

    alea_node_id_t node_id = alea_alloc_node(sys);
    if (node_id == ALEA_NODE_ID_INVALID) return ALEA_NODE_ID_INVALID;

    alea_node_t* node = &sys->nodes.data[node_id];
    ALEA_SET_OPERATION(node, ALEA_OP_INTERSECTION);
    node->operation.left = left;
    node->operation.right = right;

    /* Compute bbox as intersection of children */
    alea_bbox_t* bl = &sys->nodes.data[left].bbox;
    alea_bbox_t* br = &sys->nodes.data[right].bbox;
    node->bbox.min_x = fmax(bl->min_x, br->min_x);
    node->bbox.max_x = fmin(bl->max_x, br->max_x);
    node->bbox.min_y = fmax(bl->min_y, br->min_y);
    node->bbox.max_y = fmin(bl->max_y, br->max_y);
    node->bbox.min_z = fmax(bl->min_z, br->min_z);
    node->bbox.max_z = fmin(bl->max_z, br->max_z);

    return node_id;
}

/* Create a complement of a node */
static alea_node_id_t create_complement(alea_system_t* sys, alea_node_id_t child) {
    if (child == ALEA_NODE_ID_INVALID) return ALEA_NODE_ID_INVALID;

    alea_node_id_t node_id = alea_alloc_node(sys);
    if (node_id == ALEA_NODE_ID_INVALID) return ALEA_NODE_ID_INVALID;

    alea_node_t* node = &sys->nodes.data[node_id];
    ALEA_SET_OPERATION(node, ALEA_OP_COMPLEMENT);
    node->operation.left = child;
    node->operation.right = ALEA_NODE_ID_INVALID;

    /* Complement has same bbox (conservative) */
    node->bbox = sys->nodes.data[child].bbox;

    return node_id;
}

/* ============================================================================
 * TYPE CHECK
 * ============================================================================ */

bool alea_is_macrobody(alea_primitive_type_t type) {
    switch (type) {
        case ALEA_PRIMITIVE_RCC:
        case ALEA_PRIMITIVE_RPP:           /* RPP - axis-aligned box */
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

bool alea_is_1sheet_cone(const alea_system_t* sys, alea_node_id_t node_id) {
    if (!sys || node_id >= alea_vec_count(&sys->nodes)) return false;

    const alea_node_t* node = &sys->nodes.data[node_id];
    if (ALEA_GET_OPERATION(node) != ALEA_OP_PRIMITIVE) return false;

    alea_primitive_id_t prim_id = node->primitive.primitive_id;
    if (prim_id >= alea_vec_count(&sys->primitives)) return false;

    const alea_primitive_entry_t* prim = &sys->primitives.data[prim_id];

    switch (prim->type) {
        case ALEA_PRIMITIVE_CONE_X:
            return prim->data.cone_x.sheet_selection != 0;
        case ALEA_PRIMITIVE_CONE_Y:
            return prim->data.cone_y.sheet_selection != 0;
        case ALEA_PRIMITIVE_CONE_Z:
            return prim->data.cone_z.sheet_selection != 0;
        default:
            return false;
    }
}

/* ============================================================================
 * EXPANSION FUNCTIONS FOR EACH MACROBODY TYPE
 * ============================================================================ */

/* BOX/RPP: Axis-aligned box -> 6 planes */
static alea_node_id_t expand_box(alea_system_t* sys, const alea_box_data_t* box, int8_t sense) {
    alea_node_id_t result = ALEA_NODE_ID_INVALID;
    alea_primitive_data_t data;

    /* Create 6 axis-aligned planes */

    /* Plane 1: x = min_x (normal +x), interior has x > min_x */
    memset(&data, 0, sizeof(data));
    data.plane.a = 1.0;
    data.plane.b = 0.0;
    data.plane.c = 0.0;
    data.plane.d = -box->min_x;
    alea_node_id_t p1 = create_primitive_node_internal(sys, ALEA_PRIMITIVE_PLANE, &data, sense);
    result = p1;

    /* Plane 2: x = max_x (normal -x), interior has x < max_x */
    data.plane.a = -1.0;
    data.plane.b = 0.0;
    data.plane.c = 0.0;
    data.plane.d = box->max_x;
    alea_node_id_t p2 = create_primitive_node_internal(sys, ALEA_PRIMITIVE_PLANE, &data, sense);
    result = create_intersection(sys, result, p2);

    /* Plane 3: y = min_y (normal +y), interior has y > min_y */
    data.plane.a = 0.0;
    data.plane.b = 1.0;
    data.plane.c = 0.0;
    data.plane.d = -box->min_y;
    alea_node_id_t p3 = create_primitive_node_internal(sys, ALEA_PRIMITIVE_PLANE, &data, sense);
    result = create_intersection(sys, result, p3);

    /* Plane 4: y = max_y (normal -y), interior has y < max_y */
    data.plane.a = 0.0;
    data.plane.b = -1.0;
    data.plane.c = 0.0;
    data.plane.d = box->max_y;
    alea_node_id_t p4 = create_primitive_node_internal(sys, ALEA_PRIMITIVE_PLANE, &data, sense);
    result = create_intersection(sys, result, p4);

    /* Plane 5: z = min_z (normal +z), interior has z > min_z */
    data.plane.a = 0.0;
    data.plane.b = 0.0;
    data.plane.c = 1.0;
    data.plane.d = -box->min_z;
    alea_node_id_t p5 = create_primitive_node_internal(sys, ALEA_PRIMITIVE_PLANE, &data, sense);
    result = create_intersection(sys, result, p5);

    /* Plane 6: z = max_z (normal -z), interior has z < max_z */
    data.plane.a = 0.0;
    data.plane.b = 0.0;
    data.plane.c = -1.0;
    data.plane.d = box->max_z;
    alea_node_id_t p6 = create_primitive_node_internal(sys, ALEA_PRIMITIVE_PLANE, &data, sense);
    result = create_intersection(sys, result, p6);

    return result;
}

/* RCC: Right Circular Cylinder -> cylinder + 2 planes */
static alea_node_id_t expand_rcc(alea_system_t* sys, const alea_rcc_data_t* rcc, int8_t sense) {
    double hx = rcc->height_x, hy = rcc->height_y, hz = rcc->height_z;
    double h_len = vec3_len(hx, hy, hz);
    if (h_len < 1e-20) return ALEA_NODE_ID_INVALID;

    /* Axis unit vector */
    double ax = hx/h_len, ay = hy/h_len, az = hz/h_len;

    /* Check if axis-aligned */
    bool is_x = (fabs(ax) > 0.999 && fabs(ay) < 0.001 && fabs(az) < 0.001);
    bool is_y = (fabs(ay) > 0.999 && fabs(ax) < 0.001 && fabs(az) < 0.001);
    bool is_z = (fabs(az) > 0.999 && fabs(ax) < 0.001 && fabs(ay) < 0.001);

    alea_node_id_t cyl_node = ALEA_NODE_ID_INVALID;
    alea_primitive_data_t data;
    memset(&data, 0, sizeof(data));

    if (is_z) {
        /* Z-aligned cylinder */
        data.cyl_z.center_x = rcc->base_x;
        data.cyl_z.center_y = rcc->base_y;
        data.cyl_z.radius = rcc->radius;
        cyl_node = create_primitive_node_internal(sys, ALEA_PRIMITIVE_CYLINDER_Z, &data, sense);
    } else if (is_y) {
        /* Y-aligned cylinder */
        data.cyl_y.center_x = rcc->base_x;
        data.cyl_y.center_z = rcc->base_z;
        data.cyl_y.radius = rcc->radius;
        cyl_node = create_primitive_node_internal(sys, ALEA_PRIMITIVE_CYLINDER_Y, &data, sense);
    } else if (is_x) {
        /* X-aligned cylinder */
        data.cyl_x.center_y = rcc->base_y;
        data.cyl_x.center_z = rcc->base_z;
        data.cyl_x.radius = rcc->radius;
        cyl_node = create_primitive_node_internal(sys, ALEA_PRIMITIVE_CYLINDER_X, &data, sense);
    } else {
        /* General orientation - use quadric surface for infinite cylinder
         *
         * For cylinder with axis n=(ax,ay,az) through point P=(px,py,pz), radius r:
         * The surface equation is ||X - P - ((X-P)·n)n||² = r²
         *
         * This expands to standard quadric form:
         * Ax² + By² + Cz² + Dxy + Eyz + Fxz + Gx + Hy + Iz + J = 0
         */
        double px = rcc->base_x, py = rcc->base_y, pz = rcc->base_z;
        double r = rcc->radius;

        /* Quadric coefficients for infinite cylinder */
        double A = 1.0 - ax*ax;
        double B = 1.0 - ay*ay;
        double C = 1.0 - az*az;
        double D = -2.0 * ax * ay;
        double E = -2.0 * ay * az;
        double F = -2.0 * ax * az;

        double dot = ax*px + ay*py + az*pz;
        double G = 2.0 * (ax*dot - px);
        double H = 2.0 * (ay*dot - py);
        double I = 2.0 * (az*dot - pz);
        double J = px*px + py*py + pz*pz - dot*dot - r*r;

        data.quadric.coeffs[0] = A;
        data.quadric.coeffs[1] = B;
        data.quadric.coeffs[2] = C;
        data.quadric.coeffs[3] = D;
        data.quadric.coeffs[4] = E;
        data.quadric.coeffs[5] = F;
        data.quadric.coeffs[6] = G;
        data.quadric.coeffs[7] = H;
        data.quadric.coeffs[8] = I;
        data.quadric.coeffs[9] = J;

        cyl_node = create_primitive_node_internal(sys, ALEA_PRIMITIVE_QUADRIC, &data, sense);
    }

    if (cyl_node == ALEA_NODE_ID_INVALID) return ALEA_NODE_ID_INVALID;

    /* Base plane: ax*x + ay*y + az*z = ax*bx + ay*by + az*bz */
    double d_base = ax*rcc->base_x + ay*rcc->base_y + az*rcc->base_z;
    memset(&data, 0, sizeof(data));
    data.plane.a = ax;
    data.plane.b = ay;
    data.plane.c = az;
    data.plane.d = -d_base;
    alea_node_id_t plane_base = create_primitive_node_internal(sys, ALEA_PRIMITIVE_PLANE, &data, sense);

    /* Top plane */
    double tx = rcc->base_x + hx, ty = rcc->base_y + hy, tz = rcc->base_z + hz;
    double d_top = ax*tx + ay*ty + az*tz;
    data.plane.a = -ax;  /* Opposite normal */
    data.plane.b = -ay;
    data.plane.c = -az;
    data.plane.d = d_top;
    alea_node_id_t plane_top = create_primitive_node_internal(sys, ALEA_PRIMITIVE_PLANE, &data, sense);

    /* Build intersection: cylinder ∩ base_plane ∩ top_plane */
    alea_node_id_t result = create_intersection(sys, cyl_node, plane_base);
    result = create_intersection(sys, result, plane_top);

    return result;
}

/* BOX_GENERAL: 6 planes */
static alea_node_id_t expand_box_general(alea_system_t* sys, const alea_box_general_data_t* box, int8_t sense) {
    alea_node_id_t result = ALEA_NODE_ID_INVALID;
    alea_primitive_data_t data;

    /* Corner and edge vectors */
    double cx = box->corner_x, cy = box->corner_y, cz = box->corner_z;
    double v1x = box->v1_x, v1y = box->v1_y, v1z = box->v1_z;
    double v2x = box->v2_x, v2y = box->v2_y, v2z = box->v2_z;
    double v3x = box->v3_x, v3y = box->v3_y, v3z = box->v3_z;

    /* Compute face normals (cross products) */
    double n1x, n1y, n1z;  /* v2 x v3 */
    double n2x, n2y, n2z;  /* v3 x v1 */
    double n3x, n3y, n3z;  /* v1 x v2 */

    vec3_cross(v2x, v2y, v2z, v3x, v3y, v3z, &n1x, &n1y, &n1z);
    vec3_cross(v3x, v3y, v3z, v1x, v1y, v1z, &n2x, &n2y, &n2z);
    vec3_cross(v1x, v1y, v1z, v2x, v2y, v2z, &n3x, &n3y, &n3z);

    vec3_normalize(&n1x, &n1y, &n1z);
    vec3_normalize(&n2x, &n2y, &n2z);
    vec3_normalize(&n3x, &n3y, &n3z);

    /* Six planes, two for each normal direction */
    /* Plane 1: n1 . (p - corner) >= 0 */
    memset(&data, 0, sizeof(data));
    data.plane.a = n1x; data.plane.b = n1y; data.plane.c = n1z;
    data.plane.d = -(n1x*cx + n1y*cy + n1z*cz);
    alea_node_id_t p1 = create_primitive_node_internal(sys, ALEA_PRIMITIVE_PLANE, &data, sense);
    result = p1;

    /* Plane 2: -n1 . (p - corner - v1) >= 0 */
    double ox = cx + v1x, oy = cy + v1y, oz = cz + v1z;
    data.plane.a = -n1x; data.plane.b = -n1y; data.plane.c = -n1z;
    data.plane.d = (n1x*ox + n1y*oy + n1z*oz);
    alea_node_id_t p2 = create_primitive_node_internal(sys, ALEA_PRIMITIVE_PLANE, &data, sense);
    result = create_intersection(sys, result, p2);

    /* Planes 3,4 for n2 direction */
    data.plane.a = n2x; data.plane.b = n2y; data.plane.c = n2z;
    data.plane.d = -(n2x*cx + n2y*cy + n2z*cz);
    alea_node_id_t p3 = create_primitive_node_internal(sys, ALEA_PRIMITIVE_PLANE, &data, sense);
    result = create_intersection(sys, result, p3);

    ox = cx + v2x; oy = cy + v2y; oz = cz + v2z;
    data.plane.a = -n2x; data.plane.b = -n2y; data.plane.c = -n2z;
    data.plane.d = (n2x*ox + n2y*oy + n2z*oz);
    alea_node_id_t p4 = create_primitive_node_internal(sys, ALEA_PRIMITIVE_PLANE, &data, sense);
    result = create_intersection(sys, result, p4);

    /* Planes 5,6 for n3 direction */
    data.plane.a = n3x; data.plane.b = n3y; data.plane.c = n3z;
    data.plane.d = -(n3x*cx + n3y*cy + n3z*cz);
    alea_node_id_t p5 = create_primitive_node_internal(sys, ALEA_PRIMITIVE_PLANE, &data, sense);
    result = create_intersection(sys, result, p5);

    ox = cx + v3x; oy = cy + v3y; oz = cz + v3z;
    data.plane.a = -n3x; data.plane.b = -n3y; data.plane.c = -n3z;
    data.plane.d = (n3x*ox + n3y*oy + n3z*oz);
    alea_node_id_t p6 = create_primitive_node_internal(sys, ALEA_PRIMITIVE_PLANE, &data, sense);
    result = create_intersection(sys, result, p6);

    return result;
}

/* TRC: Truncated Right Cone -> cone/cylinder + 2 planes */
static alea_node_id_t expand_trc(alea_system_t* sys, const alea_trc_data_t* trc, int8_t sense) {
    double hx = trc->height_x, hy = trc->height_y, hz = trc->height_z;
    double h_len = vec3_len(hx, hy, hz);
    if (h_len < 1e-20) return ALEA_NODE_ID_INVALID;

    double ax = hx/h_len, ay = hy/h_len, az = hz/h_len;

    /* Check if radii are equal (cylinder) or different (cone) */
    bool is_cylinder = fabs(trc->base_radius - trc->top_radius) < 1e-10;

    /* Check if axis-aligned */
    bool is_x = (fabs(ax) > 0.999 && fabs(ay) < 0.001 && fabs(az) < 0.001);
    bool is_y = (fabs(ay) > 0.999 && fabs(ax) < 0.001 && fabs(az) < 0.001);
    bool is_z = (fabs(az) > 0.999 && fabs(ax) < 0.001 && fabs(ay) < 0.001);

    alea_node_id_t surf_node = ALEA_NODE_ID_INVALID;
    alea_primitive_data_t data;
    memset(&data, 0, sizeof(data));

    if (is_cylinder) {
        /* Cylinder case: equal radii */
        if (is_z) {
            data.cyl_z.center_x = trc->base_x;
            data.cyl_z.center_y = trc->base_y;
            data.cyl_z.radius = trc->base_radius;
            surf_node = create_primitive_node_internal(sys, ALEA_PRIMITIVE_CYLINDER_Z, &data, sense);
        } else if (is_y) {
            data.cyl_y.center_x = trc->base_x;
            data.cyl_y.center_z = trc->base_z;
            data.cyl_y.radius = trc->base_radius;
            surf_node = create_primitive_node_internal(sys, ALEA_PRIMITIVE_CYLINDER_Y, &data, sense);
        } else if (is_x) {
            data.cyl_x.center_y = trc->base_y;
            data.cyl_x.center_z = trc->base_z;
            data.cyl_x.radius = trc->base_radius;
            surf_node = create_primitive_node_internal(sys, ALEA_PRIMITIVE_CYLINDER_X, &data, sense);
        } else {
            /* Non-axis-aligned cylinder - use quadric */
            double px = trc->base_x, py = trc->base_y, pz = trc->base_z;
            double r = trc->base_radius;

            double A = 1.0 - ax*ax;
            double B = 1.0 - ay*ay;
            double C = 1.0 - az*az;
            double D = -2.0 * ax * ay;
            double E = -2.0 * ay * az;
            double F = -2.0 * ax * az;

            double dot = ax*px + ay*py + az*pz;
            double G = 2.0 * (ax*dot - px);
            double H = 2.0 * (ay*dot - py);
            double I = 2.0 * (az*dot - pz);
            double J = px*px + py*py + pz*pz - dot*dot - r*r;

            data.quadric.coeffs[0] = A;
            data.quadric.coeffs[1] = B;
            data.quadric.coeffs[2] = C;
            data.quadric.coeffs[3] = D;
            data.quadric.coeffs[4] = E;
            data.quadric.coeffs[5] = F;
            data.quadric.coeffs[6] = G;
            data.quadric.coeffs[7] = H;
            data.quadric.coeffs[8] = I;
            data.quadric.coeffs[9] = J;

            surf_node = create_primitive_node_internal(sys, ALEA_PRIMITIVE_QUADRIC, &data, sense);
        }
    } else {
        /* Cone case: different radii */
        double r1 = trc->base_radius;
        double r2 = trc->top_radius;
        double slope = (r2 - r1) / h_len;
        double t2 = slope * slope;  /* tan²(half-angle) */

        /* Compute apex position along axis */
        double apex_dist = -r1 / slope;  /* distance from base to apex along axis */
        double apex_x = trc->base_x + apex_dist * ax;
        double apex_y = trc->base_y + apex_dist * ay;
        double apex_z = trc->base_z + apex_dist * az;

        if (is_z) {
            data.cone_z.apex_x = apex_x;
            data.cone_z.apex_y = apex_y;
            data.cone_z.apex_z = apex_z;
            data.cone_z.tan_angle_sq = t2;
            data.cone_z.sheet_selection = (slope > 0) ? 1 : -1;
            surf_node = create_primitive_node_internal(sys, ALEA_PRIMITIVE_CONE_Z, &data, sense);
        } else if (is_y) {
            data.cone_y.apex_x = apex_x;
            data.cone_y.apex_y = apex_y;
            data.cone_y.apex_z = apex_z;
            data.cone_y.tan_angle_sq = t2;
            data.cone_y.sheet_selection = (slope > 0) ? 1 : -1;
            surf_node = create_primitive_node_internal(sys, ALEA_PRIMITIVE_CONE_Y, &data, sense);
        } else if (is_x) {
            data.cone_x.apex_x = apex_x;
            data.cone_x.apex_y = apex_y;
            data.cone_x.apex_z = apex_z;
            data.cone_x.tan_angle_sq = t2;
            data.cone_x.sheet_selection = (slope > 0) ? 1 : -1;
            surf_node = create_primitive_node_internal(sys, ALEA_PRIMITIVE_CONE_X, &data, sense);
        } else {
            /* Non-axis-aligned cone - use quadric
             * Cone equation: (v·n)² = t² * ||v||²
             * where v = X - apex, n = axis direction, t² = tan²(half-angle)
             *
             * Expanding to quadric form:
             * (ax² - t²)x² + (ay² - t²)y² + (az² - t²)z² + 2axay·xy + 2ayaz·yz + 2axaz·xz + G·x + H·y + I·z + J = 0
             */
            double px = apex_x, py = apex_y, pz = apex_z;

            double A = ax*ax - t2;
            double B = ay*ay - t2;
            double C = az*az - t2;
            double D = 2.0 * ax * ay;
            double E = 2.0 * ay * az;
            double F = 2.0 * ax * az;

            /* Linear terms from expanding (x-px)², etc */
            double G = -2.0*A*px - D*py - F*pz;
            double H = -D*px - 2.0*B*py - E*pz;
            double I = -F*px - E*py - 2.0*C*pz;

            /* Constant term */
            double J = A*px*px + B*py*py + C*pz*pz + D*px*py + E*py*pz + F*px*pz;

            data.quadric.coeffs[0] = A;
            data.quadric.coeffs[1] = B;
            data.quadric.coeffs[2] = C;
            data.quadric.coeffs[3] = D;
            data.quadric.coeffs[4] = E;
            data.quadric.coeffs[5] = F;
            data.quadric.coeffs[6] = G;
            data.quadric.coeffs[7] = H;
            data.quadric.coeffs[8] = I;
            data.quadric.coeffs[9] = J;

            surf_node = create_primitive_node_internal(sys, ALEA_PRIMITIVE_QUADRIC, &data, sense);
        }
    }

    if (surf_node == ALEA_NODE_ID_INVALID) return ALEA_NODE_ID_INVALID;

    /* Base plane */
    double d_base = ax*trc->base_x + ay*trc->base_y + az*trc->base_z;
    memset(&data, 0, sizeof(data));
    data.plane.a = ax;
    data.plane.b = ay;
    data.plane.c = az;
    data.plane.d = -d_base;
    alea_node_id_t plane_base = create_primitive_node_internal(sys, ALEA_PRIMITIVE_PLANE, &data, sense);

    /* Top plane */
    double tx = trc->base_x + hx, ty = trc->base_y + hy, tz = trc->base_z + hz;
    double d_top = ax*tx + ay*ty + az*tz;
    data.plane.a = -ax;
    data.plane.b = -ay;
    data.plane.c = -az;
    data.plane.d = d_top;
    alea_node_id_t plane_top = create_primitive_node_internal(sys, ALEA_PRIMITIVE_PLANE, &data, sense);

    /* Build intersection */
    alea_node_id_t result = create_intersection(sys, surf_node, plane_base);
    result = create_intersection(sys, result, plane_top);

    return result;
}

/* WED: Wedge -> 5 planes */
static alea_node_id_t expand_wed(alea_system_t* sys, const alea_wed_data_t* wed, int8_t sense) {
    alea_node_id_t result = ALEA_NODE_ID_INVALID;
    alea_primitive_data_t data;

    double vx = wed->vertex_x, vy = wed->vertex_y, vz = wed->vertex_z;
    double v1x = wed->v1_x, v1y = wed->v1_y, v1z = wed->v1_z;
    double v2x = wed->v2_x, v2y = wed->v2_y, v2z = wed->v2_z;
    double v3x = wed->v3_x, v3y = wed->v3_y, v3z = wed->v3_z;

    /* Base normal (v1 x v2) - points into the solid */
    double nb_x, nb_y, nb_z;
    vec3_cross(v1x, v1y, v1z, v2x, v2y, v2z, &nb_x, &nb_y, &nb_z);
    vec3_normalize(&nb_x, &nb_y, &nb_z);

    /* Plane 1: Base plane */
    memset(&data, 0, sizeof(data));
    data.plane.a = nb_x; data.plane.b = nb_y; data.plane.c = nb_z;
    data.plane.d = -(nb_x*vx + nb_y*vy + nb_z*vz);
    alea_node_id_t p1 = create_primitive_node_internal(sys, ALEA_PRIMITIVE_PLANE, &data, sense);
    result = p1;

    /* Plane 2: Top plane (at vertex + v3) */
    double ox = vx + v3x, oy = vy + v3y, oz = vz + v3z;
    data.plane.a = -nb_x; data.plane.b = -nb_y; data.plane.c = -nb_z;
    data.plane.d = (nb_x*ox + nb_y*oy + nb_z*oz);
    alea_node_id_t p2 = create_primitive_node_internal(sys, ALEA_PRIMITIVE_PLANE, &data, sense);
    result = create_intersection(sys, result, p2);

    /* Plane 3: v1 side (normal = v3 x v1) */
    double n3x, n3y, n3z;
    vec3_cross(v3x, v3y, v3z, v1x, v1y, v1z, &n3x, &n3y, &n3z);
    vec3_normalize(&n3x, &n3y, &n3z);
    data.plane.a = n3x; data.plane.b = n3y; data.plane.c = n3z;
    data.plane.d = -(n3x*vx + n3y*vy + n3z*vz);
    alea_node_id_t p3 = create_primitive_node_internal(sys, ALEA_PRIMITIVE_PLANE, &data, sense);
    result = create_intersection(sys, result, p3);

    /* Plane 4: v2 side (normal = v2 x v3) */
    double n4x, n4y, n4z;
    vec3_cross(v2x, v2y, v2z, v3x, v3y, v3z, &n4x, &n4y, &n4z);
    vec3_normalize(&n4x, &n4y, &n4z);
    data.plane.a = n4x; data.plane.b = n4y; data.plane.c = n4z;
    data.plane.d = -(n4x*vx + n4y*vy + n4z*vz);
    alea_node_id_t p4 = create_primitive_node_internal(sys, ALEA_PRIMITIVE_PLANE, &data, sense);
    result = create_intersection(sys, result, p4);

    /* Plane 5: Hypotenuse (from vertex+v1 to vertex+v2) */
    /* Normal = (v2-v1) x v3 */
    double dx = v2x - v1x, dy = v2y - v1y, dz = v2z - v1z;
    double n5x, n5y, n5z;
    vec3_cross(dx, dy, dz, v3x, v3y, v3z, &n5x, &n5y, &n5z);
    vec3_normalize(&n5x, &n5y, &n5z);
    ox = vx + v1x; oy = vy + v1y; oz = vz + v1z;
    data.plane.a = n5x; data.plane.b = n5y; data.plane.c = n5z;
    data.plane.d = -(n5x*ox + n5y*oy + n5z*oz);
    alea_node_id_t p5 = create_primitive_node_internal(sys, ALEA_PRIMITIVE_PLANE, &data, sense);
    result = create_intersection(sys, result, p5);

    return result;
}

/* RHP: Right Hexagonal Prism -> 8 planes (6 hex + 2 caps) */
static alea_node_id_t expand_rhp(alea_system_t* sys, const alea_rhp_data_t* rhp, int8_t sense) {
    alea_node_id_t result = ALEA_NODE_ID_INVALID;
    alea_primitive_data_t data;

    double bx = rhp->base_x, by = rhp->base_y, bz = rhp->base_z;
    double hx = rhp->height_x, hy = rhp->height_y, hz = rhp->height_z;
    double h_len = vec3_len(hx, hy, hz);
    if (h_len < 1e-20) return ALEA_NODE_ID_INVALID;

    double ax = hx/h_len, ay = hy/h_len, az = hz/h_len;

    /* Base cap plane */
    memset(&data, 0, sizeof(data));
    data.plane.a = ax; data.plane.b = ay; data.plane.c = az;
    data.plane.d = -(ax*bx + ay*by + az*bz);
    alea_node_id_t cap_base = create_primitive_node_internal(sys, ALEA_PRIMITIVE_PLANE, &data, sense);
    result = cap_base;

    /* Top cap plane */
    double tx = bx + hx, ty = by + hy, tz = bz + hz;
    data.plane.a = -ax; data.plane.b = -ay; data.plane.c = -az;
    data.plane.d = (ax*tx + ay*ty + az*tz);
    alea_node_id_t cap_top = create_primitive_node_internal(sys, ALEA_PRIMITIVE_PLANE, &data, sense);
    result = create_intersection(sys, result, cap_top);

    /* 6 hex faces - each r_i direction gives 2 parallel planes */
    double r1x = rhp->r1_x, r1y = rhp->r1_y, r1z = rhp->r1_z;
    double r2x = rhp->r2_x, r2y = rhp->r2_y, r2z = rhp->r2_z;
    double r3x = rhp->r3_x, r3y = rhp->r3_y, r3z = rhp->r3_z;

    double r1_len = vec3_len(r1x, r1y, r1z);
    double r2_len = vec3_len(r2x, r2y, r2z);
    double r3_len = vec3_len(r3x, r3y, r3z);

    /* r1 direction planes */
    double n1x = r1x/r1_len, n1y = r1y/r1_len, n1z = r1z/r1_len;
    data.plane.a = n1x; data.plane.b = n1y; data.plane.c = n1z;
    data.plane.d = -(n1x*(bx+r1x) + n1y*(by+r1y) + n1z*(bz+r1z));
    alea_node_id_t h1p = create_primitive_node_internal(sys, ALEA_PRIMITIVE_PLANE, &data, -sense);
    result = create_intersection(sys, result, h1p);

    data.plane.a = -n1x; data.plane.b = -n1y; data.plane.c = -n1z;
    data.plane.d = (n1x*(bx-r1x) + n1y*(by-r1y) + n1z*(bz-r1z));
    alea_node_id_t h1n = create_primitive_node_internal(sys, ALEA_PRIMITIVE_PLANE, &data, -sense);
    result = create_intersection(sys, result, h1n);

    /* r2 direction planes */
    double n2x = r2x/r2_len, n2y = r2y/r2_len, n2z = r2z/r2_len;
    data.plane.a = n2x; data.plane.b = n2y; data.plane.c = n2z;
    data.plane.d = -(n2x*(bx+r2x) + n2y*(by+r2y) + n2z*(bz+r2z));
    alea_node_id_t h2p = create_primitive_node_internal(sys, ALEA_PRIMITIVE_PLANE, &data, -sense);
    result = create_intersection(sys, result, h2p);

    data.plane.a = -n2x; data.plane.b = -n2y; data.plane.c = -n2z;
    data.plane.d = (n2x*(bx-r2x) + n2y*(by-r2y) + n2z*(bz-r2z));
    alea_node_id_t h2n = create_primitive_node_internal(sys, ALEA_PRIMITIVE_PLANE, &data, -sense);
    result = create_intersection(sys, result, h2n);

    /* r3 direction planes */
    double n3x = r3x/r3_len, n3y = r3y/r3_len, n3z = r3z/r3_len;
    data.plane.a = n3x; data.plane.b = n3y; data.plane.c = n3z;
    data.plane.d = -(n3x*(bx+r3x) + n3y*(by+r3y) + n3z*(bz+r3z));
    alea_node_id_t h3p = create_primitive_node_internal(sys, ALEA_PRIMITIVE_PLANE, &data, -sense);
    result = create_intersection(sys, result, h3p);

    data.plane.a = -n3x; data.plane.b = -n3y; data.plane.c = -n3z;
    data.plane.d = (n3x*(bx-r3x) + n3y*(by-r3y) + n3z*(bz-r3z));
    alea_node_id_t h3n = create_primitive_node_internal(sys, ALEA_PRIMITIVE_PLANE, &data, -sense);
    result = create_intersection(sys, result, h3n);

    return result;
}

/* ARB: Arbitrary Polyhedron -> N face planes */
static alea_node_id_t expand_arb(alea_system_t* sys, const alea_arb_data_t* arb, int8_t sense) {
    if (arb->num_faces < 4) return ALEA_NODE_ID_INVALID;

    alea_node_id_t result = ALEA_NODE_ID_INVALID;
    alea_primitive_data_t data;

    for (int f = 0; f < arb->num_faces; f++) {
        int i0 = arb->faces[f][0] - 1;
        int i1 = arb->faces[f][1] - 1;
        int i2 = arb->faces[f][2] - 1;

        if (i0 < 0 || i1 < 0 || i2 < 0 ||
            i0 >= arb->num_corners || i1 >= arb->num_corners || i2 >= arb->num_corners) {
            continue;
        }

        double x0 = arb->corners[i0][0], y0 = arb->corners[i0][1], z0 = arb->corners[i0][2];
        double x1 = arb->corners[i1][0], y1 = arb->corners[i1][1], z1 = arb->corners[i1][2];
        double x2 = arb->corners[i2][0], y2 = arb->corners[i2][1], z2 = arb->corners[i2][2];

        double e1x = x1 - x0, e1y = y1 - y0, e1z = z1 - z0;
        double e2x = x2 - x0, e2y = y2 - y0, e2z = z2 - z0;

        double nx, ny, nz;
        vec3_cross(e1x, e1y, e1z, e2x, e2y, e2z, &nx, &ny, &nz);
        double nlen = vec3_len(nx, ny, nz);
        if (nlen < 1e-20) continue;
        nx /= nlen; ny /= nlen; nz /= nlen;

        memset(&data, 0, sizeof(data));
        data.plane.a = nx;
        data.plane.b = ny;
        data.plane.c = nz;
        data.plane.d = -(nx*x0 + ny*y0 + nz*z0);

        alea_node_id_t plane_node = create_primitive_node_internal(sys, ALEA_PRIMITIVE_PLANE, &data, sense);
        result = create_intersection(sys, result, plane_node);
    }

    return result;
}

/* ============================================================================
 * 1-SHEET CONE EXPANSION
 * ============================================================================ */

/**
 * @brief Expand a 1-sheet cone into double-sheet cone + plane
 *
 * A 1-sheet cone with sheet_selection != 0 is expanded to:
 * - A double-sheet cone (same geometry, sheet_selection = 0)
 * - A plane at the apex that selects the appropriate sheet
 *
 * For sheet_selection = +1: interior = cone ∩ (axis > apex)
 * For sheet_selection = -1: interior = cone ∩ (axis < apex)
 */
static alea_node_id_t expand_1sheet_cone(alea_system_t* sys, alea_node_id_t node_id) {
    if (!sys || node_id >= alea_vec_count(&sys->nodes)) return ALEA_NODE_ID_INVALID;

    alea_node_t* node = &sys->nodes.data[node_id];
    if (ALEA_GET_OPERATION(node) != ALEA_OP_PRIMITIVE) return ALEA_NODE_ID_INVALID;

    alea_primitive_id_t prim_id = node->primitive.primitive_id;
    if (prim_id >= alea_vec_count(&sys->primitives)) return ALEA_NODE_ID_INVALID;

    const alea_primitive_entry_t* prim = &sys->primitives.data[prim_id];
    /* Always create interior (negative sense) representation - exterior handled by complement */
    const int8_t sense = -1;
    alea_primitive_data_t data;
    memset(&data, 0, sizeof(data));

    alea_node_id_t cone_node = ALEA_NODE_ID_INVALID;
    alea_node_id_t plane_node = ALEA_NODE_ID_INVALID;
    int sheet_sel = 0;
    double apex_pos = 0.0;

    switch (prim->type) {
        case ALEA_PRIMITIVE_CONE_Z: {
            const alea_cone_z_data_t* c = &prim->data.cone_z;
            sheet_sel = c->sheet_selection;
            if (sheet_sel == 0) return node_id;  /* Not a 1-sheet cone */

            /* Create double-sheet cone (sheet_selection = 0) */
            data.cone_z = *c;
            data.cone_z.sheet_selection = 0;
            cone_node = create_primitive_node_internal(sys, ALEA_PRIMITIVE_CONE_Z, &data, sense);

            /* Create plane at apex to select the correct nappe.
               With sense=-1, the plane half-space is "inside" when eval <= 0.
               For sheet_sel>0 we want z >= apex (positive nappe):
                 plane c=-1, d=apex -> eval = -z+apex, inside when -z+apex<=0 -> z>=apex
               For sheet_sel<0 we want z <= apex (negative nappe):
                 plane c=1, d=-apex -> eval = z-apex, inside when z-apex<=0 -> z<=apex */
            apex_pos = c->apex_z;
            memset(&data, 0, sizeof(data));
            data.plane.a = 0;
            data.plane.b = 0;
            data.plane.c = (sheet_sel > 0) ? -1.0 : 1.0;
            data.plane.d = (sheet_sel > 0) ? apex_pos : -apex_pos;
            plane_node = create_primitive_node_internal(sys, ALEA_PRIMITIVE_PLANE, &data, sense);
            break;
        }

        case ALEA_PRIMITIVE_CONE_Y: {
            const alea_cone_y_data_t* c = &prim->data.cone_y;
            sheet_sel = c->sheet_selection;
            if (sheet_sel == 0) return node_id;

            data.cone_y = *c;
            data.cone_y.sheet_selection = 0;
            cone_node = create_primitive_node_internal(sys, ALEA_PRIMITIVE_CONE_Y, &data, sense);

            apex_pos = c->apex_y;
            memset(&data, 0, sizeof(data));
            data.plane.a = 0;
            data.plane.b = (sheet_sel > 0) ? -1.0 : 1.0;
            data.plane.c = 0;
            data.plane.d = (sheet_sel > 0) ? apex_pos : -apex_pos;
            plane_node = create_primitive_node_internal(sys, ALEA_PRIMITIVE_PLANE, &data, sense);
            break;
        }

        case ALEA_PRIMITIVE_CONE_X: {
            const alea_cone_x_data_t* c = &prim->data.cone_x;
            sheet_sel = c->sheet_selection;
            if (sheet_sel == 0) return node_id;

            data.cone_x = *c;
            data.cone_x.sheet_selection = 0;
            cone_node = create_primitive_node_internal(sys, ALEA_PRIMITIVE_CONE_X, &data, sense);

            apex_pos = c->apex_x;
            memset(&data, 0, sizeof(data));
            data.plane.a = (sheet_sel > 0) ? -1.0 : 1.0;
            data.plane.b = 0;
            data.plane.c = 0;
            data.plane.d = (sheet_sel > 0) ? apex_pos : -apex_pos;
            plane_node = create_primitive_node_internal(sys, ALEA_PRIMITIVE_PLANE, &data, sense);
            break;
        }

        default:
            return node_id;  /* Not a cone */
    }

    if (cone_node == ALEA_NODE_ID_INVALID || plane_node == ALEA_NODE_ID_INVALID) {
        return ALEA_NODE_ID_INVALID;
    }

    /* Interior of 1-sheet cone = intersection of double-sheet cone and plane halfspace */
    return create_intersection(sys, cone_node, plane_node);
}

/* ============================================================================
 * MAIN EXPANSION FUNCTION (WITH CACHE)
 * ============================================================================ */

/**
 * @brief Internal expansion with cache support
 *
 * This is the core expansion function that uses a cache to avoid
 * creating duplicate expansions for the same primitive.
 */
static alea_node_id_t expand_macrobody_cached(alea_system_t* sys, alea_node_id_t node_id,
                                              expansion_cache_t* cache) {
    if (!sys || node_id >= alea_vec_count(&sys->nodes)) return ALEA_NODE_ID_INVALID;

    alea_node_t* node = &sys->nodes.data[node_id];
    if (ALEA_GET_OPERATION(node) != ALEA_OP_PRIMITIVE) {
        return ALEA_NODE_ID_INVALID;  /* Not a primitive node */
    }

    alea_primitive_type_t type = node->primitive.prim_type;
    alea_primitive_id_t prim_id = node->primitive.primitive_id;
    int8_t sense = node->primitive.sense;

    /* Check cache first */
    expansion_cache_entry_t* cached = cache_find(cache, prim_id);
    if (cached) {
        /* Return cached result based on sense */
        return (sense < 0) ? cached->neg_node : cached->pos_node;
    }

    /* Not cached - need to expand */
    alea_node_id_t neg_node = ALEA_NODE_ID_INVALID;

    /* Check for 1-sheet cones first */
    if (alea_is_1sheet_cone(sys, node_id)) {
        neg_node = expand_1sheet_cone(sys, node_id);
    } else if (alea_is_macrobody(type)) {
        if (prim_id >= alea_vec_count(&sys->primitives)) return ALEA_NODE_ID_INVALID;
        const alea_primitive_data_t* data = &sys->primitives.data[prim_id].data;

        /* Always expand with sense=-1 (interior representation) */
        switch (type) {
            case ALEA_PRIMITIVE_RCC:
                neg_node = expand_rcc(sys, &data->rcc, -1);
                break;

            case ALEA_PRIMITIVE_RPP:
                neg_node = expand_box(sys, &data->box, -1);
                break;

            case ALEA_PRIMITIVE_BOX:
                neg_node = expand_box_general(sys, &data->box_general, -1);
                break;

            case ALEA_PRIMITIVE_TRC:
                neg_node = expand_trc(sys, &data->trc, -1);
                break;

            case ALEA_PRIMITIVE_WED:
                neg_node = expand_wed(sys, &data->wed, -1);
                break;

            case ALEA_PRIMITIVE_RHP:
                neg_node = expand_rhp(sys, &data->rhp, -1);
                break;

            case ALEA_PRIMITIVE_ARB:
                neg_node = expand_arb(sys, &data->arb, -1);
                break;

            case ALEA_PRIMITIVE_SPH: {
                /* SPH is just a sphere - create equivalent sphere primitive */
                alea_primitive_data_t sphere_data;
                memset(&sphere_data, 0, sizeof(sphere_data));
                sphere_data.sphere.center_x = data->sph.center_x;
                sphere_data.sphere.center_y = data->sph.center_y;
                sphere_data.sphere.center_z = data->sph.center_z;
                sphere_data.sphere.radius = data->sph.radius;
                neg_node = create_primitive_node_internal(sys, ALEA_PRIMITIVE_SPHERE, &sphere_data, -1);
                break;
            }

            case ALEA_PRIMITIVE_ELL:
            case ALEA_PRIMITIVE_REC:
                /* These stay as-is - no expansion */
                return node_id;

            default:
                return ALEA_NODE_ID_INVALID;
        }
    } else {
        /* Not a macrobody or 1-sheet cone */
        return node_id;
    }

    if (neg_node == ALEA_NODE_ID_INVALID) {
        return ALEA_NODE_ID_INVALID;
    }

    /* Create complement for exterior (positive sense) representation */
    alea_node_id_t pos_node = create_complement(sys, neg_node);
    if (pos_node == ALEA_NODE_ID_INVALID) {
        return ALEA_NODE_ID_INVALID;
    }

    /* Cache both representations */
    cache_add(cache, prim_id, neg_node, pos_node);

    /* Return the appropriate one based on requested sense */
    return (sense < 0) ? neg_node : pos_node;
}

/**
 * @brief Public expansion function (creates temporary cache)
 */
alea_node_id_t alea_expand_macrobody(alea_system_t* sys, alea_node_id_t node_id) {
    expansion_cache_t cache;
    cache_init(&cache);

    alea_node_id_t result = expand_macrobody_cached(sys, node_id, &cache);

    cache_free(&cache);
    return result;
}

/* ============================================================================
 * RECURSIVE TREE EXPANSION (WITH CACHE)
 * ============================================================================ */

/**
 * @brief Internal recursive expansion with cache
 */
static alea_node_id_t expand_all_cached(alea_system_t* sys, alea_node_id_t root_id,
                                        expansion_cache_t* cache) {
    if (!sys || root_id >= alea_vec_count(&sys->nodes)) return ALEA_NODE_ID_INVALID;

    alea_node_t* node = &sys->nodes.data[root_id];
    alea_operation_t op = ALEA_GET_OPERATION(node);

    if (op == ALEA_OP_PRIMITIVE) {
        /* Try to expand if it's a macrobody or 1-sheet cone */
        return expand_macrobody_cached(sys, root_id, cache);
    }

    if (op == ALEA_OP_COMPLEMENT) {
        /* Complement has single child in left */
        alea_node_id_t new_child = expand_all_cached(sys, node->operation.left, cache);
        if (new_child != node->operation.left) {
            /* Create new complement node with expanded child */
            alea_node_id_t new_node = alea_alloc_node(sys);
            if (new_node == ALEA_NODE_ID_INVALID) return ALEA_NODE_ID_INVALID;
            alea_node_t* nn = &sys->nodes.data[new_node];
            ALEA_SET_OPERATION(nn, ALEA_OP_COMPLEMENT);
            nn->operation.left = new_child;
            nn->operation.right = ALEA_NODE_ID_INVALID;
            nn->bbox = sys->nodes.data[new_child].bbox;
            return new_node;
        }
        return root_id;
    }

    /* Binary operation */
    alea_node_id_t new_left = expand_all_cached(sys, node->operation.left, cache);
    alea_node_id_t new_right = expand_all_cached(sys, node->operation.right, cache);

    if (new_left != node->operation.left || new_right != node->operation.right) {
        /* Create new operation node with expanded children */
        alea_node_id_t new_node = alea_alloc_node(sys);
        if (new_node == ALEA_NODE_ID_INVALID) return ALEA_NODE_ID_INVALID;

        alea_node_t* nn = &sys->nodes.data[new_node];
        ALEA_SET_OPERATION(nn, op);
        nn->operation.left = new_left;
        nn->operation.right = new_right;

        /* Recompute bbox based on operation */
        alea_bbox_t* bl = &sys->nodes.data[new_left].bbox;
        alea_bbox_t* br = &sys->nodes.data[new_right].bbox;
        if (op == ALEA_OP_UNION) {
            nn->bbox.min_x = fmin(bl->min_x, br->min_x);
            nn->bbox.max_x = fmax(bl->max_x, br->max_x);
            nn->bbox.min_y = fmin(bl->min_y, br->min_y);
            nn->bbox.max_y = fmax(bl->max_y, br->max_y);
            nn->bbox.min_z = fmin(bl->min_z, br->min_z);
            nn->bbox.max_z = fmax(bl->max_z, br->max_z);
        } else {
            nn->bbox.min_x = fmax(bl->min_x, br->min_x);
            nn->bbox.max_x = fmin(bl->max_x, br->max_x);
            nn->bbox.min_y = fmax(bl->min_y, br->min_y);
            nn->bbox.max_y = fmin(bl->max_y, br->max_y);
            nn->bbox.min_z = fmax(bl->min_z, br->min_z);
            nn->bbox.max_z = fmin(bl->max_z, br->max_z);
        }
        return new_node;
    }

    return root_id;
}

/**
 * @brief Public function to expand all macrobodies in a tree
 */
alea_node_id_t alea_expand_all_macrobodies(alea_system_t* sys, alea_node_id_t root_id) {
    expansion_cache_t cache;
    cache_init(&cache);

    alea_node_id_t result = expand_all_cached(sys, root_id, &cache);

    cache_free(&cache);
    return result;
}

/* ============================================================================
 * IMMEDIATE EXPANSION AT CREATION TIME
 * ============================================================================ */

/**
 * @brief Expand a macrobody at creation time
 *
 * Called from surface_conv.c when a macrobody is parsed. Creates the
 * expanded CSG tree (component primitives) immediately and returns
 * both sense nodes.
 *
 * The component primitives do NOT get MCNP surface IDs assigned here -
 * that happens at export time.
 *
 * The pos_node is created using De Morgan transformation, so:
 *   neg_node = -cyl ∩ -plane_base ∩ -plane_top  (intersection)
 *   pos_node = +cyl ∪ +plane_base ∪ +plane_top  (union, no COMPLEMENT wrapper)
 *
 * @param sys CSG system
 * @param type Macrobody primitive type
 * @param data Macrobody primitive data
 * @param out_neg_node Output: negative sense node (interior)
 * @param out_pos_node Output: positive sense node (exterior, De Morgan form)
 * @return 0 on success, -1 on error
 */
int alea_expand_macrobody_immediate(alea_system_t* sys,
                                    alea_primitive_type_t type,
                                    const alea_primitive_data_t* data,
                                    alea_node_id_t* out_neg_node,
                                    alea_node_id_t* out_pos_node) {
    if (!sys || !data || !out_neg_node || !out_pos_node) return -1;
    if (!alea_is_macrobody(type)) return -1;

    alea_node_id_t neg_node = ALEA_NODE_ID_INVALID;

    /* Expand based on macrobody type - always create interior (sense=-1) first */
    switch (type) {
        case ALEA_PRIMITIVE_RCC:
            neg_node = expand_rcc(sys, &data->rcc, -1);
            break;

        case ALEA_PRIMITIVE_RPP:
            neg_node = expand_box(sys, &data->box, -1);
            break;

        case ALEA_PRIMITIVE_BOX:
            neg_node = expand_box_general(sys, &data->box_general, -1);
            break;

        case ALEA_PRIMITIVE_TRC:
            neg_node = expand_trc(sys, &data->trc, -1);
            break;

        case ALEA_PRIMITIVE_WED:
            neg_node = expand_wed(sys, &data->wed, -1);
            break;

        case ALEA_PRIMITIVE_RHP:
            neg_node = expand_rhp(sys, &data->rhp, -1);
            break;

        case ALEA_PRIMITIVE_ARB:
            neg_node = expand_arb(sys, &data->arb, -1);
            break;

        case ALEA_PRIMITIVE_SPH: {
            /* SPH is just a sphere - create equivalent sphere primitive */
            alea_primitive_data_t sphere_data;
            memset(&sphere_data, 0, sizeof(sphere_data));
            sphere_data.sphere.center_x = data->sph.center_x;
            sphere_data.sphere.center_y = data->sph.center_y;
            sphere_data.sphere.center_z = data->sph.center_z;
            sphere_data.sphere.radius = data->sph.radius;
            neg_node = create_primitive_node_internal(sys, ALEA_PRIMITIVE_SPHERE, &sphere_data, -1);
            break;
        }

        case ALEA_PRIMITIVE_ELL:
        case ALEA_PRIMITIVE_REC:
            /* These don't expand - return error, caller should handle directly */
            return -1;

        default:
            return -1;
    }

    if (neg_node == ALEA_NODE_ID_INVALID) {
        return -1;
    }

    /* Create exterior (positive sense) using De Morgan transformation */
    /* This gives us a union of positive-sense primitives instead of #(intersection) */
    alea_node_id_t pos_node = alea_apply_demorgan(sys, neg_node);
    if (pos_node == ALEA_NODE_ID_INVALID) {
        return -1;
    }

    *out_neg_node = neg_node;
    *out_pos_node = pos_node;
    return 0;
}

/**
 * @brief Check if cone data represents a 1-sheet cone
 */
static bool is_1sheet_cone_data(alea_primitive_type_t type, const alea_primitive_data_t* data) {
    switch (type) {
        case ALEA_PRIMITIVE_CONE_X:
            return data->cone_x.sheet_selection != 0;
        case ALEA_PRIMITIVE_CONE_Y:
            return data->cone_y.sheet_selection != 0;
        case ALEA_PRIMITIVE_CONE_Z:
            return data->cone_z.sheet_selection != 0;
        default:
            return false;
    }
}

int alea_expand_1sheet_cone_immediate(alea_system_t* sys,
                                      alea_primitive_type_t type,
                                      const alea_primitive_data_t* data,
                                      alea_node_id_t* out_neg_node,
                                      alea_node_id_t* out_pos_node) {
    if (!sys || !data || !out_neg_node || !out_pos_node) return -1;

    /* Verify this is a 1-sheet cone */
    if (!is_1sheet_cone_data(type, data)) return -1;

    alea_node_id_t cone_node = ALEA_NODE_ID_INVALID;
    alea_node_id_t plane_node = ALEA_NODE_ID_INVALID;
    int sheet_sel = 0;
    double apex_pos = 0.0;
    alea_primitive_data_t pdata;
    memset(&pdata, 0, sizeof(pdata));

    switch (type) {
        case ALEA_PRIMITIVE_CONE_Z: {
            const alea_cone_z_data_t* c = &data->cone_z;
            sheet_sel = c->sheet_selection;
            apex_pos = c->apex_z;

            /* Create double-sheet cone (sheet_selection = 0) */
            pdata.cone_z = *c;
            pdata.cone_z.sheet_selection = 0;
            cone_node = create_primitive_node_internal(sys, ALEA_PRIMITIVE_CONE_Z, &pdata, -1);

            /* Create plane at apex to select the correct nappe.
               See expand_1sheet_cone() for the sign derivation. */
            memset(&pdata, 0, sizeof(pdata));
            pdata.plane.a = 0;
            pdata.plane.b = 0;
            pdata.plane.c = (sheet_sel > 0) ? -1.0 : 1.0;
            pdata.plane.d = (sheet_sel > 0) ? apex_pos : -apex_pos;
            plane_node = create_primitive_node_internal(sys, ALEA_PRIMITIVE_PLANE, &pdata, -1);
            break;
        }

        case ALEA_PRIMITIVE_CONE_Y: {
            const alea_cone_y_data_t* c = &data->cone_y;
            sheet_sel = c->sheet_selection;
            apex_pos = c->apex_y;

            pdata.cone_y = *c;
            pdata.cone_y.sheet_selection = 0;
            cone_node = create_primitive_node_internal(sys, ALEA_PRIMITIVE_CONE_Y, &pdata, -1);

            memset(&pdata, 0, sizeof(pdata));
            pdata.plane.a = 0;
            pdata.plane.b = (sheet_sel > 0) ? -1.0 : 1.0;
            pdata.plane.c = 0;
            pdata.plane.d = (sheet_sel > 0) ? apex_pos : -apex_pos;
            plane_node = create_primitive_node_internal(sys, ALEA_PRIMITIVE_PLANE, &pdata, -1);
            break;
        }

        case ALEA_PRIMITIVE_CONE_X: {
            const alea_cone_x_data_t* c = &data->cone_x;
            sheet_sel = c->sheet_selection;
            apex_pos = c->apex_x;

            pdata.cone_x = *c;
            pdata.cone_x.sheet_selection = 0;
            cone_node = create_primitive_node_internal(sys, ALEA_PRIMITIVE_CONE_X, &pdata, -1);

            memset(&pdata, 0, sizeof(pdata));
            pdata.plane.a = (sheet_sel > 0) ? -1.0 : 1.0;
            pdata.plane.b = 0;
            pdata.plane.c = 0;
            pdata.plane.d = (sheet_sel > 0) ? apex_pos : -apex_pos;
            plane_node = create_primitive_node_internal(sys, ALEA_PRIMITIVE_PLANE, &pdata, -1);
            break;
        }

        default:
            return -1;
    }

    if (cone_node == ALEA_NODE_ID_INVALID || plane_node == ALEA_NODE_ID_INVALID) {
        return -1;
    }

    /* Interior of 1-sheet cone = intersection of double-sheet cone and plane halfspace */
    alea_node_id_t neg_node = create_intersection(sys, cone_node, plane_node);
    if (neg_node == ALEA_NODE_ID_INVALID) {
        return -1;
    }

    /* Create exterior using De Morgan transformation */
    alea_node_id_t pos_node = alea_apply_demorgan(sys, neg_node);
    if (pos_node == ALEA_NODE_ID_INVALID) {
        return -1;
    }

    *out_neg_node = neg_node;
    *out_pos_node = pos_node;
    return 0;
}
