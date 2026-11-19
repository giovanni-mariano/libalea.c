// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

#include "bbox.h"
#include "alea.h"
#include "primitive_desc.h"
#include "core/alea_eval.h"
#include "core/alea_universe.h"  /* for alea_matrix_t (struct alea_matrix) */
#include "util/math.h"
#include <math.h>
#include <float.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

/**
 * @file bbox.c
 * @brief Implementation of bounding box operations
 *
 * Uses the primitive descriptor table for primitive bbox computation.
 */

#define BBOX_LARGE 1.0e10

alea_bbox_t alea_primitive_bbox(
    alea_primitive_type_t type,
    const alea_primitive_data_t* data
) {
    if (!data) return alea_bbox_empty();

    const alea_primitive_desc_t* desc = alea_primitive_get_desc(type);
    if (!desc || !desc->bbox) {
        return alea_bbox_infinite();
    }

    return desc->bbox(data);
}

// ============================================================================
// BOUNDING BOX OPERATIONS
// ============================================================================

alea_bbox_t alea_bbox_union(const alea_bbox_t* a, const alea_bbox_t* b) {
    if (!a || !b) {
        return alea_bbox_empty();
    }
    
    alea_bbox_t result = {
        MIN(a->min_x, b->min_x), MAX(a->max_x, b->max_x),
        MIN(a->min_y, b->min_y), MAX(a->max_y, b->max_y),
        MIN(a->min_z, b->min_z), MAX(a->max_z, b->max_z)
    };
    
    return result;
}

alea_bbox_t alea_bbox_intersection(const alea_bbox_t* a, const alea_bbox_t* b) {
    if (!a || !b) {
        return alea_bbox_empty();
    }
    
    alea_bbox_t result = {
        MAX(a->min_x, b->min_x), MIN(a->max_x, b->max_x),
        MAX(a->min_y, b->min_y), MIN(a->max_y, b->max_y),
        MAX(a->min_z, b->min_z), MIN(a->max_z, b->max_z)
    };
    
    // Check if result is valid (intersection exists)
    if (result.min_x > result.max_x ||
        result.min_y > result.max_y ||
        result.min_z > result.max_z) {
        return alea_bbox_empty();
    }
    
    return result;
}

bool alea_bbox_intersects(const alea_bbox_t* a, const alea_bbox_t* b) {
    if (!a || !b) {
        return false;
    }
    
    // Separating axis test - boxes intersect if NOT separated on any axis
    return !(a->max_x < b->min_x || a->min_x > b->max_x ||
             a->max_y < b->min_y || a->min_y > b->max_y ||
             a->max_z < b->min_z || a->min_z > b->max_z);
}

bool alea_bbox_contains_point(const alea_bbox_t* bbox, double x, double y, double z) {
    if (!bbox) {
        return false;
    }

    return x >= bbox->min_x && x <= bbox->max_x &&
           y >= bbox->min_y && y <= bbox->max_y &&
           z >= bbox->min_z && z <= bbox->max_z;
}

bool alea_bbox_contains(const alea_bbox_t* outer, const alea_bbox_t* inner) {
    if (!outer || !inner) {
        return false;
    }

    return inner->min_x >= outer->min_x && inner->max_x <= outer->max_x &&
           inner->min_y >= outer->min_y && inner->max_y <= outer->max_y &&
           inner->min_z >= outer->min_z && inner->max_z <= outer->max_z;
}

void alea_bbox_center(const alea_bbox_t* bbox, double* out_x, double* out_y, double* out_z) {
    if (!bbox || !out_x || !out_y || !out_z) {
        return;
    }
    
    *out_x = (bbox->min_x + bbox->max_x) * 0.5;
    *out_y = (bbox->min_y + bbox->max_y) * 0.5;
    *out_z = (bbox->min_z + bbox->max_z) * 0.5;
}

void alea_bbox_dimensions(const alea_bbox_t* bbox, double* out_width, double* out_height, double* out_depth) {
    if (!bbox || !out_width || !out_height || !out_depth) {
        return;
    }
    
    *out_width = bbox->max_x - bbox->min_x;
    *out_height = bbox->max_y - bbox->min_y;
    *out_depth = bbox->max_z - bbox->min_z;
}

double alea_bbox_volume(const alea_bbox_t* bbox) {
    if (!bbox || !alea_bbox_is_valid(bbox)) {
        return 0.0;
    }
    
    double width = bbox->max_x - bbox->min_x;
    double height = bbox->max_y - bbox->min_y;
    double depth = bbox->max_z - bbox->min_z;
    
    return width * height * depth;
}

double alea_bbox_surface_area(const alea_bbox_t* bbox) {
    if (!bbox || !alea_bbox_is_valid(bbox)) {
        return 0.0;
    }
    
    double width = bbox->max_x - bbox->min_x;
    double height = bbox->max_y - bbox->min_y;
    double depth = bbox->max_z - bbox->min_z;
    
    return 2.0 * (width*height + height*depth + depth*width);
}

bool alea_bbox_is_valid(const alea_bbox_t* bbox) {
    if (!bbox) {
        return false;
    }
    
    return bbox->min_x <= bbox->max_x &&
           bbox->min_y <= bbox->max_y &&
           bbox->min_z <= bbox->max_z;
}

alea_bbox_t alea_bbox_empty(void) {
    // Empty bbox has min > max
    return (alea_bbox_t){
        DBL_MAX, -DBL_MAX,
        DBL_MAX, -DBL_MAX,
        DBL_MAX, -DBL_MAX
    };
}

alea_bbox_t alea_bbox_infinite(void) {
    return (alea_bbox_t){
        -BBOX_LARGE, BBOX_LARGE,
        -BBOX_LARGE, BBOX_LARGE,
        -BBOX_LARGE, BBOX_LARGE
    };
}

alea_bbox_t alea_bbox_transform(const alea_bbox_t* bbox,
                                const struct alea_matrix* mat) {
    if (!bbox || !mat) return alea_bbox_empty();

    /* Decompose rotation matrix into min/max contributions to avoid
     * enumerating all 8 corners (18 mul + 18 add vs 192 mul-adds).
     * Matrix layout: m[row*4 + col], row-major 3x4 with translation in col 3. */
    const double bmin[3] = {bbox->min_x, bbox->min_y, bbox->min_z};
    const double bmax[3] = {bbox->max_x, bbox->max_y, bbox->max_z};
    double rmin[3], rmax[3];

    for (int i = 0; i < 3; i++) {
        rmin[i] = rmax[i] = mat->m[i * 4 + 3]; /* translation */
        for (int j = 0; j < 3; j++) {
            double e = mat->m[i * 4 + j] * bmin[j];
            double f = mat->m[i * 4 + j] * bmax[j];
            if (e < f) { rmin[i] += e; rmax[i] += f; }
            else       { rmin[i] += f; rmax[i] += e; }
        }
    }

    return (alea_bbox_t){rmin[0], rmax[0], rmin[1], rmax[1], rmin[2], rmax[2]};
}

alea_bbox_t alea_get_bbox(const alea_system_t* sys, alea_node_id_t id) {
    if (id == ALEA_NODE_ID_INVALID || !sys) {
        return alea_bbox_empty();
    }
    
    if (id >= alea_vec_count(&sys->nodes)) {
        return alea_bbox_empty();
    }
    const alea_node_t* node = &sys->nodes.data[id];
    
    alea_operation_t op = ALEA_GET_OPERATION(node);
    
    if (op == ALEA_OP_PRIMITIVE) {
        /* Use the primitive_id from the node, not the node_id */
        uint32_t prim_id = node->primitive.primitive_id;
        if (prim_id >= alea_vec_count(&sys->primitives)) {
            return alea_bbox_empty();
        }
        alea_primitive_data_t data;
        if (!alea_primitive_copy_data(sys, prim_id, &data)) {
            return alea_bbox_empty();
        }
        /* Compute effective sense: inverted flips the sense */
        int8_t effective_sense = node->primitive.inverted ? -node->primitive.sense : node->primitive.sense;
        return alea_halfspace_bbox(node->primitive.prim_type, &data, effective_sense);
    }
    
    // Recursively compute for operations
    alea_bbox_t left_bbox = alea_get_bbox(sys, node->operation.left);
    alea_bbox_t right_bbox = alea_get_bbox(sys, node->operation.right);
    
    switch (op) {
        case ALEA_OP_UNION:
            return alea_bbox_union(&left_bbox, &right_bbox);
        case ALEA_OP_INTERSECTION:
            return alea_bbox_intersection(&left_bbox, &right_bbox);
        case ALEA_OP_DIFFERENCE:
            return left_bbox; // Conservative
        default:
            return left_bbox;
    }
}

alea_bbox_t alea_halfspace_bbox(
    alea_primitive_type_t type,
    const alea_primitive_data_t* data,
    int8_t sense
) {
    if (!data) return alea_bbox_empty();

    alea_bbox_t surf_bbox = alea_primitive_bbox(type, data);

    /* For bounded primitives (sphere, box, rcc, etc.):
     * - sense < 0 (inside): return the primitive's bbox
     * - sense > 0 (outside): return infinite bbox
     */
    switch (type) {
        case ALEA_PRIMITIVE_SPHERE:
        case ALEA_PRIMITIVE_RPP:
        case ALEA_PRIMITIVE_RCC:
        case ALEA_PRIMITIVE_TRC:
        case ALEA_PRIMITIVE_SPH:
        case ALEA_PRIMITIVE_BOX:
        case ALEA_PRIMITIVE_RHP:
        case ALEA_PRIMITIVE_ELL:
        case ALEA_PRIMITIVE_WED:
        case ALEA_PRIMITIVE_ARB:
        case ALEA_PRIMITIVE_REC:
            /* Bounded primitives: inside has finite bbox, outside is infinite */
            if (sense < 0) {
                return surf_bbox;  /* Inside: bounded */
            } else {
                return alea_bbox_infinite();  /* Outside: infinite */
            }

        case ALEA_PRIMITIVE_PLANE: {
            /* For axis-aligned planes, compute proper half-space bbox */
            double a = data->plane.a;
            double b = data->plane.b;
            double c = data->plane.c;
            double d = data->plane.d;

            /* Z-aligned plane: a ≈ 0, b ≈ 0, |c| ≈ 1 */
            if (fabs(a) < 1e-9 && fabs(b) < 1e-9 && fabs(c) > 1e-9) {
                double z0 = -d / c;
                /* Plane equation: c*z + d = 0 → z = -d/c
                 * Positive sense: c*z + d > 0
                 * If c > 0: z > z0; if c < 0: z < z0 */
                if ((c > 0 && sense > 0) || (c < 0 && sense < 0)) {
                    /* z > z0 */
                    return (alea_bbox_t){-BBOX_LARGE, BBOX_LARGE, -BBOX_LARGE, BBOX_LARGE, z0, BBOX_LARGE};
                } else {
                    /* z < z0 */
                    return (alea_bbox_t){-BBOX_LARGE, BBOX_LARGE, -BBOX_LARGE, BBOX_LARGE, -BBOX_LARGE, z0};
                }
            }

            /* Y-aligned plane: a ≈ 0, |b| > 0, c ≈ 0 */
            if (fabs(a) < 1e-9 && fabs(b) > 1e-9 && fabs(c) < 1e-9) {
                double y0 = -d / b;
                if ((b > 0 && sense > 0) || (b < 0 && sense < 0)) {
                    return (alea_bbox_t){-BBOX_LARGE, BBOX_LARGE, y0, BBOX_LARGE, -BBOX_LARGE, BBOX_LARGE};
                } else {
                    return (alea_bbox_t){-BBOX_LARGE, BBOX_LARGE, -BBOX_LARGE, y0, -BBOX_LARGE, BBOX_LARGE};
                }
            }

            /* X-aligned plane: |a| > 0, b ≈ 0, c ≈ 0 */
            if (fabs(a) > 1e-9 && fabs(b) < 1e-9 && fabs(c) < 1e-9) {
                double x0 = -d / a;
                if ((a > 0 && sense > 0) || (a < 0 && sense < 0)) {
                    return (alea_bbox_t){x0, BBOX_LARGE, -BBOX_LARGE, BBOX_LARGE, -BBOX_LARGE, BBOX_LARGE};
                } else {
                    return (alea_bbox_t){-BBOX_LARGE, x0, -BBOX_LARGE, BBOX_LARGE, -BBOX_LARGE, BBOX_LARGE};
                }
            }

            /* General plane: infinite in all directions */
            return alea_bbox_infinite();
        }

        case ALEA_PRIMITIVE_CYLINDER_X:
        case ALEA_PRIMITIVE_CYLINDER_Y:
        case ALEA_PRIMITIVE_CYLINDER_Z:
            /* Infinite cylinders: inside has radial bounds, outside is infinite */
            if (sense < 0) {
                return surf_bbox;  /* Inside: bounded radially */
            } else {
                return alea_bbox_infinite();  /* Outside: infinite */
            }

        case ALEA_PRIMITIVE_CONE_X:
        case ALEA_PRIMITIVE_CONE_Y:
        case ALEA_PRIMITIVE_CONE_Z:
            /* Cones are infinite along the axis; inside is radially unbounded too */
            return alea_bbox_infinite();

        case ALEA_PRIMITIVE_TORUS_X: {
            /* Tori are bounded if they have finite axial extent */
            const alea_torus_data_t* t = &data->torus;
            /* Check if torus has finite extent (axial_semiwidth_B > 0 means bounded) */
            bool bounded = (t->axial_semiwidth_B > 0.0) ||
                           (surf_bbox.min_x > -BBOX_LARGE * 0.9 && surf_bbox.max_x < BBOX_LARGE * 0.9);
            if (bounded) {
                if (sense < 0) {
                    return surf_bbox;  /* Inside: bounded */
                } else {
                    return alea_bbox_infinite();  /* Outside: infinite */
                }
            }
            return alea_bbox_infinite();
        }

        case ALEA_PRIMITIVE_TORUS_Y: {
            const alea_torus_data_t* t = &data->torus;
            bool bounded = (t->axial_semiwidth_B > 0.0) ||
                           (surf_bbox.min_y > -BBOX_LARGE * 0.9 && surf_bbox.max_y < BBOX_LARGE * 0.9);
            if (bounded) {
                if (sense < 0) {
                    return surf_bbox;
                } else {
                    return alea_bbox_infinite();
                }
            }
            return alea_bbox_infinite();
        }

        case ALEA_PRIMITIVE_TORUS_Z: {
            const alea_torus_data_t* t = &data->torus;
            bool bounded = (t->axial_semiwidth_B > 0.0) ||
                           (surf_bbox.min_z > -BBOX_LARGE * 0.9 && surf_bbox.max_z < BBOX_LARGE * 0.9);
            if (bounded) {
                if (sense < 0) {
                    return surf_bbox;
                } else {
                    return alea_bbox_infinite();
                }
            }
            return alea_bbox_infinite();
        }

        case ALEA_PRIMITIVE_QUADRIC: {
            /* Quadrics: check if the computed bbox is finite */
            bool is_bounded = (surf_bbox.min_x > -BBOX_LARGE * 0.9 && surf_bbox.max_x < BBOX_LARGE * 0.9 &&
                               surf_bbox.min_y > -BBOX_LARGE * 0.9 && surf_bbox.max_y < BBOX_LARGE * 0.9 &&
                               surf_bbox.min_z > -BBOX_LARGE * 0.9 && surf_bbox.max_z < BBOX_LARGE * 0.9);
            if (is_bounded) {
                /* Bounded quadric (ellipsoid, etc.) */
                if (sense < 0) {
                    return surf_bbox;  /* Inside: bounded */
                } else {
                    return alea_bbox_infinite();  /* Outside: infinite */
                }
            }
            /* Unbounded quadric (hyperboloid, paraboloid, etc.) */
            return alea_bbox_infinite();
        }

        default:
            /* Unknown or complex: use conservative infinite bbox */
            return alea_bbox_infinite();
    }
}

/* Forward declaration for numerical tightening (defined below) */
static int tighten_bbox_plane_lp(const alea_system_t* sys,
                                    alea_node_id_t root,
                                    alea_bbox_t* out);

/* ============================================================================
 * BOUNDING BOX TIGHTENING
 *
 * Binary-search each axis to find the true geometric extent of a cell.
 * The initial bbox comes from the CSG tree node (conservative).
 * Each iteration splits the bbox and tests whether one half is entirely
 * outside the cell.  Converges in O(log(extent/tol)) steps per axis.
 * ============================================================================ */

void alea_tighten_tree_bbox(const alea_system_t* sys,
                            alea_node_id_t root,
                            const alea_bbox_t* start,
                            double tol,
                            alea_bbox_t* out) {
    *out = *start;

    /* For each axis: 0=X, 1=Y, 2=Z */
    for (int axis = 0; axis < 3; axis++) {
        double* lo;
        double* hi;
        switch (axis) {
            case 0: lo = &out->min_x; hi = &out->max_x; break;
            case 1: lo = &out->min_y; hi = &out->max_y; break;
            default: lo = &out->min_z; hi = &out->max_z; break;
        }

        /* Tighten lower bound (increase min) */
        {
            double search_lo = *lo;
            double search_hi = *hi;
            while (search_hi - search_lo > tol) {
                double mid = (search_lo + search_hi) * 0.5;

                /* Test lower half [lo, mid] */
                alea_bbox_t test = *out;
                switch (axis) {
                    case 0: test.max_x = mid; break;
                    case 1: test.max_y = mid; break;
                    default: test.max_z = mid; break;
                }

                alea_box_relation_t rel = alea_tree_box_relation(sys, root, &test);
                if (rel == ALEA_RELATION_POSITIVE) {
                    /* Lower half is empty -> tighten */
                    *lo = mid;
                    search_lo = mid;
                } else {
                    /* Cell has content in lower half -> stop search */
                    search_hi = mid;
                }
            }
        }

        /* Tighten upper bound (decrease max) */
        {
            double search_lo = *lo;
            double search_hi = *hi;
            while (search_hi - search_lo > tol) {
                double mid = (search_lo + search_hi) * 0.5;

                /* Test upper half [mid, hi] */
                alea_bbox_t test = *out;
                switch (axis) {
                    case 0: test.min_x = mid; break;
                    case 1: test.min_y = mid; break;
                    default: test.min_z = mid; break;
                }

                alea_box_relation_t rel = alea_tree_box_relation(sys, root, &test);
                if (rel == ALEA_RELATION_POSITIVE) {
                    /* Upper half is empty -> tighten */
                    *hi = mid;
                    search_hi = mid;
                } else {
                    /* Cell has content in upper half -> stop search */
                    search_lo = mid;
                }
            }
        }
    }
}

int alea_compute_bounding_sphere(alea_system_t* sys,
                                 double tol,
                                 double* cx, double* cy, double* cz,
                                 double* radius) {
    if (!sys || tol <= 0.0 || !cx || !cy || !cz || !radius) return -1;

    size_t n_cells = alea_vec_count(&sys->cells);
    if (n_cells == 0) return -1;

    /* Union of all tightened cell bboxes */
    double model_min_x = 1e30, model_min_y = 1e30, model_min_z = 1e30;
    double model_max_x = -1e30, model_max_y = -1e30, model_max_z = -1e30;
    int n_bounded = 0;

    for (size_t i = 0; i < n_cells; i++) {
        const alea_cell_entry_t* cell = &sys->cells.data[i];
        if (cell->root_node_id == ALEA_NODE_ID_INVALID) continue;

        const alea_bbox_t box_v = alea_node_bbox_get(&sys->nodes.data[cell->root_node_id].bbox);
        const alea_bbox_t* box = &box_v;
        if (box->min_x > box->max_x) continue;
        double dx = box->max_x - box->min_x;
        double dy = box->max_y - box->min_y;
        double dz = box->max_z - box->min_z;
        if (dx > 9e5 || dy > 9e5 || dz > 9e5) continue;

        /* First check: is this cell entirely empty? */
        alea_box_relation_t rel = alea_tree_box_relation(sys, cell->root_node_id, box);
        if (rel == ALEA_RELATION_POSITIVE) continue;  /* empty cell */

        alea_bbox_t tight;
        alea_tighten_tree_bbox(sys, cell->root_node_id, box, tol, &tight);

        /* Skip cells that tightened to nothing */
        if (tight.min_x >= tight.max_x ||
            tight.min_y >= tight.max_y ||
            tight.min_z >= tight.max_z) continue;

        if (tight.min_x < model_min_x) model_min_x = tight.min_x;
        if (tight.min_y < model_min_y) model_min_y = tight.min_y;
        if (tight.min_z < model_min_z) model_min_z = tight.min_z;
        if (tight.max_x > model_max_x) model_max_x = tight.max_x;
        if (tight.max_y > model_max_y) model_max_y = tight.max_y;
        if (tight.max_z > model_max_z) model_max_z = tight.max_z;
        n_bounded++;
    }

    /* Pass 2: try numerical tightening for cells skipped due to infinite bbox */
    for (size_t i = 0; i < n_cells; i++) {
        const alea_cell_entry_t* cell = &sys->cells.data[i];
        if (cell->root_node_id == ALEA_NODE_ID_INVALID) continue;

        const alea_bbox_t box_v = alea_node_bbox_get(&sys->nodes.data[cell->root_node_id].bbox);
        const alea_bbox_t* box = &box_v;
        if (box->min_x > box->max_x) continue;
        double dx2 = box->max_x - box->min_x;
        double dy2 = box->max_y - box->min_y;
        double dz2 = box->max_z - box->min_z;
        if (dx2 <= 9e5 && dy2 <= 9e5 && dz2 <= 9e5) continue; /* already handled */

        alea_bbox_t tight;
        if (alea_tighten_bbox_numerical(sys, cell->root_node_id, tol, &tight) != 0)
            continue;

        if (tight.min_x >= tight.max_x ||
            tight.min_y >= tight.max_y ||
            tight.min_z >= tight.max_z) continue;

        /* Skip results that are still too large */
        double tdx = tight.max_x - tight.min_x;
        double tdy = tight.max_y - tight.min_y;
        double tdz = tight.max_z - tight.min_z;
        if (tdx > 9e5 || tdy > 9e5 || tdz > 9e5) continue;

        if (tight.min_x < model_min_x) model_min_x = tight.min_x;
        if (tight.min_y < model_min_y) model_min_y = tight.min_y;
        if (tight.min_z < model_min_z) model_min_z = tight.min_z;
        if (tight.max_x > model_max_x) model_max_x = tight.max_x;
        if (tight.max_y > model_max_y) model_max_y = tight.max_y;
        if (tight.max_z > model_max_z) model_max_z = tight.max_z;
        n_bounded++;
    }

    if (n_bounded == 0) return -1;

    /* Center of model bbox */
    *cx = (model_min_x + model_max_x) * 0.5;
    *cy = (model_min_y + model_max_y) * 0.5;
    *cz = (model_min_z + model_max_z) * 0.5;

    /* Radius: half-diagonal with small margin */
    double dx = model_max_x - model_min_x;
    double dy = model_max_y - model_min_y;
    double dz = model_max_z - model_min_z;
    *radius = 0.5 * sqrt(dx * dx + dy * dy + dz * dz) * 1.01;

    return 0;
}

int alea_tighten_all_bboxes(alea_system_t* sys, double tol) {
    if (!sys || tol <= 0.0) return -1;

    size_t n_cells = alea_vec_count(&sys->cells);
    int tightened = 0;

    for (size_t i = 0; i < n_cells; i++) {
        alea_cell_entry_t* cell = &sys->cells.data[i];
        if (cell->root_node_id == ALEA_NODE_ID_INVALID) continue;

        alea_node_bbox_t* node_box = &sys->nodes.data[cell->root_node_id].bbox;
        alea_bbox_t box_v = alea_node_bbox_get(node_box);
        alea_bbox_t* box = &box_v;
        if (box->min_x > box->max_x) continue;
        double dx = box->max_x - box->min_x;
        double dy = box->max_y - box->min_y;
        double dz = box->max_z - box->min_z;
        if (dx > 9e5 || dy > 9e5 || dz > 9e5) {
            /* Use globally certified mixed constraints, never a clipped search domain. */
            alea_bbox_t lp_box;
            if (alea_tighten_bbox_mixed_constraints(sys, cell->root_node_id, &lp_box) == 0) {
                alea_bbox_t tight;
                alea_tighten_tree_bbox(sys, cell->root_node_id, &lp_box, tol, &tight);
                if (tight.min_x < tight.max_x &&
                    tight.min_y < tight.max_y &&
                    tight.min_z < tight.max_z) {
                    alea_node_bbox_set(node_box, &tight);
                    tightened++;
                }
            }
            continue;
        }

        alea_box_relation_t rel = alea_tree_box_relation(sys, cell->root_node_id, box);
        if (rel == ALEA_RELATION_POSITIVE) continue;

        alea_bbox_t tight;
        alea_tighten_tree_bbox(sys, cell->root_node_id, box, tol, &tight);

        if (tight.min_x < tight.max_x &&
            tight.min_y < tight.max_y &&
            tight.min_z < tight.max_z) {
            alea_node_bbox_set(node_box, &tight);
            tightened++;
        }
    }

    if (tightened) alea_system_invalidate_query_caches(sys, ALEA_CACHE_ALL & ~ALEA_CACHE_CELL_SURFACES);
    return tightened;
}

/* ============================================================================
 * NUMERICAL BOUNDING BOX TIGHTENING
 *
 * Two strategies for cells with infinite analytical bboxes:
 * 1. certified plane LP (fast, exact for pure plane intersections)
 * 2. Octree recursive bisection (general-purpose fallback)
 * ============================================================================ */

/* Deterministic incremental linear programs in at most three dimensions.
 * The optimizer proposes extrema; nonnegative dual combinations independently
 * certify enclosing bounds. Artificial box faces are never accepted as global
 * bounds. Uncertain, empty and unbounded systems leave the caller's box alone. */
#define LP_MAX_PLANES 65536u
#define LP_MAX_WORK 16000000u
#define LP_BIG 1e12

typedef struct { double a[3], b; } lp_constraint;
typedef struct { size_t work; lp_constraint* scratch[2]; } lp_workspace;

static int lp_solve(int d, const lp_constraint* h, size_t n,
                    const double c[3], const double bounds[3][2],
                    double x[3], lp_workspace* w) {
    if (alea_interrupted()) return -1;
    if (d == 1) {
        double lo=bounds[0][0], hi=bounds[0][1];
        for(size_t i=0;i<n;i++) {
            if (++w->work>LP_MAX_WORK) return -1;
            double a=h[i].a[0], b=h[i].b;
            if(a==0) { if(b<0) return -1; }
            else if(a>0) hi=fmin(hi,b/a);
            else lo=fmax(lo,b/a);
        }
        if(lo>hi || !isfinite(lo) || !isfinite(hi)) return -1;
        x[0]=c[0]<0?hi:lo; return 0;
    }
    for(int k=0;k<d;k++) x[k]=bounds[k][c[k]<0];
    for(size_t i=0;i<n;i++) {
        if(++w->work>LP_MAX_WORK) return -1;
        double v=0,mag=fabs(h[i].b);
        for(int k=0;k<d;k++){v+=h[i].a[k]*x[k];mag+=fabs(h[i].a[k]*x[k]);}
        if(v<=h[i].b+32*DBL_EPSILON*(1+mag)) continue;
        int pivot=0, keep[2],nk=0;
        for(int k=1;k<d;k++) if(fabs(h[i].a[k])>fabs(h[i].a[pivot]))pivot=k;
        if(h[i].a[pivot]==0) return -1;
        for(int k=0;k<d;k++) if(k!=pivot) keep[nk++]=k;
        double r[2], origin=h[i].b/h[i].a[pivot];
        for(int k=0;k<nk;k++) r[k]=h[i].a[keep[k]]/h[i].a[pivot];
        lp_constraint* sub=w->scratch[d-2];
        for(size_t j=0;j<i;j++) {
            sub[j]=(lp_constraint){{0,0,0}, h[j].b-h[j].a[pivot]*origin};
            for(int k=0;k<nk;k++)sub[j].a[k]=h[j].a[keep[k]]-h[j].a[pivot]*r[k];
        }
        sub[i]=(lp_constraint){{0,0,0},bounds[pivot][1]-origin};
        sub[i+1]=(lp_constraint){{0,0,0},origin-bounds[pivot][0]};
        double nc[3]={0},nb[3][2]={{0}},y[3]={0};
        for(int k=0;k<nk;k++) {
            sub[i].a[k]=-r[k];sub[i+1].a[k]=r[k];
            nc[k]=c[keep[k]]-c[pivot]*r[k];
            nb[k][0]=bounds[keep[k]][0];nb[k][1]=bounds[keep[k]][1];
        }
        if(lp_solve(d-1,sub,i+2,nc,nb,y,w))return -1;
        x[pivot]=origin;
        for(int k=0;k<nk;k++){x[keep[k]]=y[k];x[pivot]-=r[k]*y[k];}
    }
    for(int k=0;k<d;k++)if(!isfinite(x[k]))return -1;
    return 0;
}

/* Directed interval operations used only for the independent certificate.
 * Each elementary result is widened by one representable number. */
typedef struct { double lo,hi; } lp_interval;
static lp_interval lp_add(lp_interval a,lp_interval b) {
    return (lp_interval){nextafter(a.lo+b.lo,-INFINITY),nextafter(a.hi+b.hi,INFINITY)};
}
static lp_interval lp_product(double a,double b) {
    double p=a*b;return (lp_interval){nextafter(p,-INFINITY),nextafter(p,INFINITY)};
}

static lp_interval lp_neg(lp_interval a) { return (lp_interval){-a.hi,-a.lo}; }
static lp_interval lp_mul(lp_interval a,lp_interval b) {
    double p[4]={a.lo*b.lo,a.lo*b.hi,a.hi*b.lo,a.hi*b.hi};
    double lo=p[0],hi=p[0];
    for(int i=0;i<4;i++) {
        if(!isfinite(p[i]))return (lp_interval){-INFINITY,INFINITY};
        lo=fmin(lo,p[i]);hi=fmax(hi,p[i]);
    }
    return (lp_interval){nextafter(lo,-INFINITY),nextafter(hi,INFINITY)};
}
static lp_interval lp_sub(lp_interval a,lp_interval b) { return lp_add(a,lp_neg(b)); }
static lp_interval lp_scale(lp_interval a,double b) {
    return lp_mul(a,(lp_interval){b,b});
}

/* A certified slab l <= a.x <= u can close an otherwise unbounded quadric.
 * For q(x)<=0, q(x)+mu*(a.x-mid)^2-mu*radius^2<=0 is necessary.
 * If that new quadratic is provably positive definite, its ellipsoid is a
 * GLOBAL enclosure of the conjunction. No small eigenvalue is erased and no
 * artificial LP search box or user window participates in this proof.
 * Interval Sylvester minors certify definiteness; an interval inverse and
 * completing the square enclose the center and coordinate support radii. */
static int lp_quadric_enclosure(const double coeff[10],double sense,
                                const lp_constraint* slab,double lower,
                                alea_bbox_t* out) {
    double q[10],scale=0;
    for(int i=0;i<10;i++) {
        q[i]=-sense*coeff[i];
        if(!isfinite(q[i]))return -1;
        if(i<6)scale=fmax(scale,fabs(q[i]));
    }
    lp_interval a[3][3],b[3],constant={q[9],q[9]};
    for(int i=0;i<3;i++) { a[i][i]=(lp_interval){q[i],q[i]};b[i]=(lp_interval){q[6+i],q[6+i]}; }
    a[0][1]=a[1][0]=lp_product(q[3],0.5);
    a[1][2]=a[2][1]=lp_product(q[4],0.5);
    a[0][2]=a[2][0]=lp_product(q[5],0.5);
    if(slab) {
        double norm=0;for(int i=0;i<3;i++)norm+=slab->a[i]*slab->a[i];
        double mu=scale/norm,upper=slab->b;
        if(!isfinite(mu) || mu<=0 || !isfinite(lower) || lower>upper)return -1;
        double mid=lower/2+upper/2;
        lp_interval dl=lp_sub((lp_interval){lower,lower},(lp_interval){mid,mid});
        lp_interval du=lp_sub((lp_interval){upper,upper},(lp_interval){mid,mid});
        double radius=fmax(fmax(fabs(dl.lo),fabs(dl.hi)),fmax(fabs(du.lo),fabs(du.hi)));
        /* Subtract an upper radius squared: enlarging the ellipsoid is safe. */
        double r2=lp_product(radius,radius).hi;
        constant=lp_add(constant,lp_scale(lp_sub(lp_product(mid,mid),(lp_interval){r2,r2}),mu));
        for(int i=0;i<3;i++) {
            b[i]=lp_sub(b[i],lp_scale(lp_scale(lp_product(mid,slab->a[i]),mu),2));
            for(int j=i;j<3;j++)
                a[j][i]=a[i][j]=lp_add(a[i][j],lp_scale(lp_product(slab->a[i],slab->a[j]),mu));
        }
    }
    lp_interval cof[3][3];
    cof[0][0]=lp_sub(lp_mul(a[1][1],a[2][2]),lp_mul(a[1][2],a[1][2]));
    cof[1][1]=lp_sub(lp_mul(a[0][0],a[2][2]),lp_mul(a[0][2],a[0][2]));
    cof[2][2]=lp_sub(lp_mul(a[0][0],a[1][1]),lp_mul(a[0][1],a[0][1]));
    cof[0][1]=cof[1][0]=lp_sub(lp_mul(a[0][2],a[1][2]),lp_mul(a[0][1],a[2][2]));
    cof[0][2]=cof[2][0]=lp_sub(lp_mul(a[0][1],a[1][2]),lp_mul(a[0][2],a[1][1]));
    cof[1][2]=cof[2][1]=lp_sub(lp_mul(a[0][1],a[0][2]),lp_mul(a[0][0],a[1][2]));
    lp_interval det={0,0};
    for(int i=0;i<3;i++)det=lp_add(det,lp_mul(a[0][i],cof[0][i]));
    if(!(a[0][0].lo>0 && cof[2][2].lo>0 && det.lo>0) || !isfinite(det.hi))return -1;
    lp_interval reciprocal={nextafter(1/det.hi,-INFINITY),nextafter(1/det.lo,INFINITY)};
    lp_interval inverse[3][3],center[3],radius2=lp_neg(constant);
    for(int i=0;i<3;i++) {
        lp_interval row={0,0};
        for(int j=0;j<3;j++) {
            inverse[i][j]=lp_mul(cof[i][j],reciprocal);
            row=lp_add(row,lp_mul(inverse[i][j],b[j]));
        }
        center[i]=lp_scale(row,-0.5);
        radius2=lp_add(radius2,lp_scale(lp_mul(b[i],row),0.25));
    }
    if(!isfinite(radius2.hi) || radius2.hi<0)return -1; /* No empty-set claim. */
    double lo[3],hi[3];
    for(int i=0;i<3;i++) {
        double support2=lp_product(radius2.hi,inverse[i][i].hi).hi;
        if(!isfinite(support2) || support2<0)return -1;
        double support=nextafter(sqrt(support2),INFINITY);
        lo[i]=nextafter(center[i].lo-support,-INFINITY);
        hi[i]=nextafter(center[i].hi+support,INFINITY);
        if(!isfinite(lo[i]) || !isfinite(hi[i]) || lo[i]>hi[i])return -1;
    }
    *out=(alea_bbox_t){lo[0],hi[0],lo[1],hi[1],lo[2],hi[2]};return 0;
}

static int lp_certificate(const lp_constraint* h,size_t n,const double x[3],
                          int axis,int sign,double* upper,double* residual) {
    size_t active[32],count=0;
    for(size_t i=0;i<n && count<32;i++) {
        double v=0,mag=fabs(h[i].b);
        for(int k=0;k<3;k++){v+=h[i].a[k]*x[k];mag+=fabs(h[i].a[k]*x[k]);}
        if(fabs(v-h[i].b)<=1e-7*(1+mag))active[count++]=i;
    }
    /* Zero columns allow certificates supported on one or two constraints. */
    for(size_t i=0;i<count;i++)for(size_t j=i;j<count;j++)for(size_t k=j;k<count;k++) {
        size_t ids[3]={active[i],active[j],active[k]};
        double lambda[3]={0};
        int m=(i==j&&j==k)?1:(i==j||j==k)?2:3;
        if(m==2){ids[0]=active[i];ids[1]=active[k];}
        if(m==1) {
            double a=h[ids[0]].a[axis];if(a==0)continue;
            lambda[0]=sign/a;
        } else if(m==2) {
            double best=0;int p=0,q=1;
            for(int a=0;a<3;a++)for(int b=a+1;b<3;b++) {
                double det=h[ids[0]].a[a]*h[ids[1]].a[b]-h[ids[1]].a[a]*h[ids[0]].a[b];
                if(fabs(det)>fabs(best)){best=det;p=a;q=b;}
            }
            if(best==0)continue;
            double cp=p==axis?sign:0,cq=q==axis?sign:0;
            lambda[0]=(cp*h[ids[1]].a[q]-cq*h[ids[1]].a[p])/best;
            lambda[1]=(cq*h[ids[0]].a[p]-cp*h[ids[0]].a[q])/best;
        } else {
            mat3 A;
            for(int col=0;col<3;col++)for(int row=0;row<3;row++)A.m[col][row]=h[ids[col]].a[row];
            vec3 c={axis==0?sign:0,axis==1?sign:0,axis==2?sign:0},v;
            if(mat3_solve_cramer(A,c,&v))continue;
            lambda[0]=v.x;lambda[1]=v.y;lambda[2]=v.z;
        }
        int valid=1;for(int t=0;t<m;t++)if(lambda[t]<0 || !isfinite(lambda[t]))valid=0;
        if(!valid)continue;
        lp_interval sum={0,0};
        for(int t=0;t<m;t++)sum=lp_add(sum,lp_product(lambda[t],h[ids[t]].b));
        double error=0;
        for(int a=0;a<3;a++) {
            lp_interval r={a==axis?sign:0,a==axis?sign:0};
            for(int t=0;t<m;t++)r=lp_add(r,lp_product(-lambda[t],h[ids[t]].a[a]));
            double err=fmax(fabs(r.lo),fabs(r.hi));
            error=nextafter(error+err,INFINITY);
        }
        if(isfinite(sum.hi) && error<0.25){*upper=sum.hi;*residual=error;return 0;}
    }
    return -1;
}

/* Necessary constraints only: an enclosing box for a curved primitive is a
 * relaxation, never an approximation of its membership predicate. */
static size_t lp_primitive_constraints(const alea_primitive_data_t* data,
        alea_primitive_type_t type, double sense, lp_constraint h[6]) {
    if(type==ALEA_PRIMITIVE_PLANE) {
        h[0]=(lp_constraint){{-sense*data->plane.a,-sense*data->plane.b,
                              -sense*data->plane.c},sense*data->plane.d};
        return 1;
    }
    if(sense>=0)return 0; /* Exterior curved half-spaces have no such bounds. */
    double center[3]={0},radius=0,lo[3],hi[3];
    int skip=-1;
    switch(type) {
    case ALEA_PRIMITIVE_SPHERE:
        center[0]=data->sphere.center_x;center[1]=data->sphere.center_y;
        center[2]=data->sphere.center_z;radius=data->sphere.radius;break;
    case ALEA_PRIMITIVE_CYLINDER_X:
        skip=0;center[1]=data->cyl_x.center_y;center[2]=data->cyl_x.center_z;
        radius=data->cyl_x.radius;break;
    case ALEA_PRIMITIVE_CYLINDER_Y:
        skip=1;center[0]=data->cyl_y.center_x;center[2]=data->cyl_y.center_z;
        radius=data->cyl_y.radius;break;
    case ALEA_PRIMITIVE_CYLINDER_Z:
        skip=2;center[0]=data->cyl_z.center_x;center[1]=data->cyl_z.center_y;
        radius=data->cyl_z.radius;break;
    case ALEA_PRIMITIVE_RPP:
        lo[0]=data->box.min_x;hi[0]=data->box.max_x;
        lo[1]=data->box.min_y;hi[1]=data->box.max_y;
        lo[2]=data->box.min_z;hi[2]=data->box.max_z;break;
    default:return 0;
    }
    if(type!=ALEA_PRIMITIVE_RPP) {
        if(!isfinite(radius) || radius<=0)return 0;
        for(int k=0;k<3;k++) {
            lo[k]=nextafter(center[k]-radius,-INFINITY);
            hi[k]=nextafter(center[k]+radius,INFINITY);
        }
    }
    size_t n=0;
    for(int k=0;k<3;k++)if(k!=skip) {
        if(!isfinite(lo[k]) || !isfinite(hi[k]) || lo[k]>hi[k])return 0;
        h[n]=(lp_constraint){{0},hi[k]};h[n++].a[k]=1;
        h[n]=(lp_constraint){{0},-lo[k]};h[n++].a[k]=-1;
    }
    return n;
}

static int tighten_bbox_lp(const alea_system_t* sys,
                           alea_node_id_t root,alea_bbox_t* out,int mixed) {
    /* Iterative collection bounds both memory and pathological DAG expansion. */
    size_t capacity=64,n=0,top=0,visits=0;
    /* Bounded extra work: at most four quadrics and four exact opposite slabs.
     * Unsupported or uncertain candidates simply contribute no constraints. */
    alea_node_id_t quadrics[4];size_t quadric_count=0;
    lp_constraint* h=malloc(capacity*sizeof(*h));
    alea_node_id_t* stack=malloc(LP_MAX_PLANES*sizeof(*stack));
    if(!h || !stack){free(h);free(stack);return -1;}
    stack[top++]=root;int status=-1;
    while(top) {
        if(++visits>2*LP_MAX_PLANES || alea_interrupted())goto done;
        alea_node_id_t id=stack[--top];if(id>=sys->nodes.count)goto done;
        const alea_node_t* node=&sys->nodes.data[id];
        alea_operation_t op=ALEA_GET_OPERATION(node);
        if(op==ALEA_OP_INTERSECTION) {
            if(top+2>LP_MAX_PLANES)goto done;
            stack[top++]=node->operation.right;stack[top++]=node->operation.left;continue;
        }
        if(mixed && op==ALEA_OP_DIFFERENCE) {
            /* A minus B is contained in A. No constraints from B are valid. */
            stack[top++]=node->operation.left;continue;
        }
        if(op!=ALEA_OP_PRIMITIVE) {
            if(mixed)continue; /* Do not descend into unions or complements. */
            goto done;
        }
        if(!mixed && node->primitive.prim_type!=ALEA_PRIMITIVE_PLANE)goto done;
        alea_primitive_data_t data;
        if(!alea_primitive_copy_data(sys,node->primitive.primitive_id,&data))goto done;
        double sense=node->primitive.sense*(node->primitive.inverted?-1:1);
        if(mixed && node->primitive.prim_type==ALEA_PRIMITIVE_QUADRIC && quadric_count<4) {
            int negative=0;
            for(int k=0;k<3;k++)if(-sense*data.quadric.coeffs[k]<0)negative++;
            if(negative<2)quadrics[quadric_count++]=id;
        }
        lp_constraint leaf[6];
        size_t count=lp_primitive_constraints(&data,node->primitive.prim_type,sense,leaf);
        if(n+count>LP_MAX_PLANES)goto done;
        if(n+count>capacity){capacity*=2;void* grown=realloc(h,capacity*sizeof(*h));if(!grown)goto done;h=grown;}
        for(size_t j=0;j<count;j++) {
            for(int k=0;k<3;k++)if(!isfinite(leaf[j].a[k]))goto done;
            if(!isfinite(leaf[j].b))goto done;
            h[n++]=leaf[j];
        }
    }
    if(quadric_count) {
        lp_constraint slabs[4];double lowers[4];size_t slab_count=0;
        size_t search=n<64?n:64;
        for(size_t i=0;i<search && slab_count<4;i++)for(size_t j=i+1;j<search && slab_count<4;j++) {
            int opposite=1,nonzero=0;
            for(int k=0;k<3;k++) {
                if(h[i].a[k]!=-h[j].a[k])opposite=0;
                if(h[i].a[k]!=0)nonzero=1;
            }
            if(opposite && nonzero && -h[j].b<=h[i].b) {
                slabs[slab_count]=h[i];lowers[slab_count++]=-h[j].b;
            }
        }
        for(size_t q=0;q<quadric_count;q++) {
            const alea_node_t* node=&sys->nodes.data[quadrics[q]];
            alea_primitive_data_t data;
            if(!alea_primitive_copy_data(sys,node->primitive.primitive_id,&data))goto done;
            double sense=node->primitive.sense*(node->primitive.inverted?-1:1);
            alea_bbox_t enclosure;int bounded=0;
            for(size_t s=0;s<=slab_count;s++) {
                alea_bbox_t box;
                if(lp_quadric_enclosure(data.quadric.coeffs,sense,s?&slabs[s-1]:NULL,
                                        s?lowers[s-1]:0,&box))continue;
                enclosure=bounded?alea_bbox_intersection(&enclosure,&box):box;bounded=1;
            }
            if(!bounded || !alea_bbox_is_valid(&enclosure))continue;
            if(n+6>LP_MAX_PLANES)goto done;
            if(n+6>capacity) {
                capacity*=2;void* grown=realloc(h,capacity*sizeof(*h));
                if(!grown)goto done;
                h=grown;
            }
            double lo[3]={enclosure.min_x,enclosure.min_y,enclosure.min_z};
            double hi[3]={enclosure.max_x,enclosure.max_y,enclosure.max_z};
            for(int k=0;k<3;k++) {
                h[n]=(lp_constraint){{0},hi[k]};h[n++].a[k]=1;
                h[n]=(lp_constraint){{0},-lo[k]};h[n++].a[k]=-1;
            }
        }
    }
    if(n<4)goto done;
    uint64_t seed=UINT64_C(0x2545f4914f6cdd1d);
    for(size_t i=n-1;i>0;i--){seed^=seed<<13;seed^=seed>>7;seed^=seed<<17;size_t j=seed%(i+1);lp_constraint tmp=h[i];h[i]=h[j];h[j]=tmp;}
    lp_workspace w={0};
    w.scratch[0]=malloc((n+4)*sizeof(*h));w.scratch[1]=malloc((n+4)*sizeof(*h));
    if(!w.scratch[0] || !w.scratch[1]){free(w.scratch[0]);free(w.scratch[1]);goto done;}
    double bounds[3][2]={{-LP_BIG,LP_BIG},{-LP_BIG,LP_BIG},{-LP_BIG,LP_BIG}},v[6],errors[6];
    int ok=1;
    for(int axis=0;axis<3 && ok;axis++)for(int side=0;side<2;side++) {
        int sign=side?1:-1;double c[3]={0},x[3]={0},u;
        c[axis]=-sign;
        if(lp_solve(3,h,n,c,bounds,x,&w) || lp_certificate(h,n,x,axis,sign,&u,&errors[2*axis+side])){ok=0;break;}
        v[2*axis+side]=u;
    }
    free(w.scratch[0]);free(w.scratch[1]);
    if(ok) {
        /* For every feasible point, M=max |x_i| satisfies
         * M <= max(rhs) + max(residual)*M. This closes the bound
         * without assuming that the artificial search box encloses it. */
        double rhs=0,error=0;
        for(int i=0;i<6;i++){rhs=fmax(rhs,fabs(v[i]));error=fmax(error,errors[i]);}
        double denominator=nextafter(1-error,-INFINITY);
        double extent=nextafter(rhs/denominator,INFINITY);
        if(!isfinite(extent))goto done;
        for(int i=0;i<6;i++)v[i]=nextafter(v[i]+nextafter(errors[i]*extent,INFINITY),INFINITY);
        alea_bbox_t result={-v[0],v[1],-v[2],v[3],-v[4],v[5]};
        if(alea_bbox_is_valid(&result)){*out=result;status=0;}
    }
done:
    free(h);free(stack);return status;
}

static int tighten_bbox_plane_lp(const alea_system_t* sys,
                                 alea_node_id_t root,alea_bbox_t* out) {
    return tighten_bbox_lp(sys,root,out,0);
}

int alea_tighten_bbox_mixed_constraints(const alea_system_t* sys,
                                       alea_node_id_t root,alea_bbox_t* out) {
    if(!sys || !out)return -1;
    return tighten_bbox_lp(sys,root,out,1);
}

int alea_tighten_bbox_plane_constraints(const alea_system_t* sys,
                                       alea_node_id_t root,double tol,
                                       alea_bbox_t* out) {
    if(!sys || !out || !isfinite(tol) || tol<=0)return -1;
    /* The LP certificates already establish global enclosure. */
    return tighten_bbox_plane_lp(sys,root,out);
}

/**
 * Octree recursive bisection: start from a large search box, recursively
 * subdivide, and accumulate bbox from non-empty leaf boxes.
 */
static void octree_collect_bbox(
    const alea_system_t* sys,
    alea_node_id_t root,
    const alea_bbox_t* box,
    int depth,
    alea_bbox_t* accum
) {
    alea_box_relation_t rel = alea_tree_box_relation(sys, root, box);
    if (rel == ALEA_RELATION_POSITIVE) return;  /* entirely outside — skip */

    if (rel == ALEA_RELATION_NEGATIVE || depth <= 0) {
        /* This box has content — union into accumulator */
        *accum = alea_bbox_union(accum, box);
        return;
    }

    /* Bisect along longest axis */
    double dx = box->max_x - box->min_x;
    double dy = box->max_y - box->min_y;
    double dz = box->max_z - box->min_z;

    alea_bbox_t a = *box, b = *box;
    if (dx >= dy && dx >= dz) {
        double mid = (box->min_x + box->max_x) * 0.5;
        a.max_x = mid;  b.min_x = mid;
    } else if (dy >= dz) {
        double mid = (box->min_y + box->max_y) * 0.5;
        a.max_y = mid;  b.min_y = mid;
    } else {
        double mid = (box->min_z + box->max_z) * 0.5;
        a.max_z = mid;  b.min_z = mid;
    }

    octree_collect_bbox(sys, root, &a, depth - 1, accum);
    octree_collect_bbox(sys, root, &b, depth - 1, accum);
}

/**
 * Octree-based bbox tightening.
 * Starts from [-1e6, 1e6]^3, recursively bisects to depth 10.
 * Returns 0 on success, -1 if no content found.
 */
static int tighten_bbox_octree(
    const alea_system_t* sys,
    alea_node_id_t root,
    alea_bbox_t* out
) {
    alea_bbox_t search = {
        -1e6, 1e6,
        -1e6, 1e6,
        -1e6, 1e6
    };

    alea_bbox_t accum = alea_bbox_empty();
    octree_collect_bbox(sys, root, &search, 10, &accum);

    if (accum.min_x > accum.max_x) return -1;  /* nothing found */

    *out = accum;
    return 0;
}

int alea_tighten_bbox_numerical(
    const alea_system_t* sys,
    alea_node_id_t root,
    double tol,
    alea_bbox_t* out
) {
    if (!sys || !out || root == ALEA_NODE_ID_INVALID) return -1;

    alea_bbox_t coarse;

    /* Strategy 1: globally certified plane and mixed-primitive bounds */
    if (alea_tighten_bbox_mixed_constraints(sys, root, &coarse) == 0) {
        /* Refine with binary search */
        alea_tighten_tree_bbox(sys, root, &coarse, tol, out);
        return 0;
    }

    /* Strategy 2: Octree recursive bisection (general fallback) */
    if (tighten_bbox_octree(sys, root, &coarse) == 0) {
        /* Refine with binary search */
        alea_tighten_tree_bbox(sys, root, &coarse, tol, out);
        return 0;
    }

    return -1;
}
