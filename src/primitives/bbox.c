// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

#include "bbox.h"
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

alea_bbox_t alea_get_bbox(alea_system_t* sys, alea_node_id_t id) {
    if (id == ALEA_NODE_ID_INVALID || !sys) {
        return alea_bbox_empty();
    }
    
    alea_node_t* node = alea_get_node(sys, id);
    if (!node) {
        return alea_bbox_empty();
    }
    
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
static int tighten_bbox_vertex_enum(const alea_system_t* sys,
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
            /* Try LP vertex enumeration only (fast, no octree) */
            alea_bbox_t lp_box;
            if (tighten_bbox_vertex_enum(sys, cell->root_node_id, &lp_box) == 0) {
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

    return tightened;
}

/* ============================================================================
 * NUMERICAL BOUNDING BOX TIGHTENING
 *
 * Two strategies for cells with infinite analytical bboxes:
 * 1. LP vertex enumeration (fast, exact for pure plane intersections)
 * 2. Octree recursive bisection (general-purpose fallback)
 * ============================================================================ */

/* Maximum planes for vertex enumeration */
#define MAX_LP_PLANES 64

/* Half-space constraint: n·x <= d  (i.e. eff_sense*(ax+by+cz+d) <= 0) */
typedef struct {
    double a, b, c, d;  /* n·x <= d  means  a*x + b*y + c*z <= d */
} halfspace_t;

/**
 * Collect plane half-space constraints from a pure intersection tree.
 * Returns the number of planes collected, or -1 if the tree contains
 * non-plane primitives or union/difference nodes.
 */
static int collect_plane_constraints(
    const alea_system_t* sys,
    alea_node_id_t node_id,
    halfspace_t* planes,
    int max_planes,
    int count
) {
    if (node_id == ALEA_NODE_ID_INVALID || node_id >= alea_vec_count(&sys->nodes))
        return -1;

    const alea_node_t* node = &sys->nodes.data[node_id];
    alea_operation_t op = ALEA_GET_OPERATION(node);

    if (op == ALEA_OP_INTERSECTION) {
        count = collect_plane_constraints(sys, node->operation.left, planes, max_planes, count);
        if (count < 0) return -1;
        count = collect_plane_constraints(sys, node->operation.right, planes, max_planes, count);
        return count;
    }

    if (op != ALEA_OP_PRIMITIVE) return -1;  /* union/difference/complement */
    if (node->primitive.prim_type != ALEA_PRIMITIVE_PLANE) return -1;

    if (count >= max_planes) return -1;

    uint32_t prim_id = node->primitive.primitive_id;
    if (prim_id >= alea_vec_count(&sys->primitives)) return -1;

    alea_primitive_data_t data;
    if (!alea_primitive_copy_data(sys, prim_id, &data)) return -1;
    const alea_plane_data_t* p = &data.plane;
    int8_t eff_sense = node->primitive.sense;
    if (node->primitive.inverted) eff_sense = -eff_sense;

    /* Half-space convention:
     * eff_sense = -1 (negative side): ax+by+cz+d <= 0 → ax+by+cz <= -d
     * eff_sense = +1 (positive side): ax+by+cz+d >= 0 → -ax-by-cz <= d
     * Unified: (-eff_sense*a)*x + (-eff_sense*b)*y + (-eff_sense*c)*z
     *          <= -(-eff_sense*d) = eff_sense*d (wait, let's be explicit):
     * Store as: n·x <= rhs where n = -eff_sense*(a,b,c), rhs = eff_sense*d */
    double ns = -(double)eff_sense;
    planes[count].a = ns * p->a;
    planes[count].b = ns * p->b;
    planes[count].c = ns * p->c;
    planes[count].d = (double)eff_sense * p->d;  /* RHS */
    count++;
    return count;
}

/**
 * Check if a point satisfies all half-space constraints.
 */
static bool point_feasible(const halfspace_t* planes, int n, vec3 pt) {
    for (int i = 0; i < n; i++) {
        double val = planes[i].a * pt.x + planes[i].b * pt.y + planes[i].c * pt.z;
        if (val > planes[i].d + 1e-9) return false;
    }
    return true;
}

/**
 * LP vertex enumeration: enumerate all triple-plane intersection vertices,
 * keep feasible ones, compute their bounding box.
 *
 * Returns 0 on success, -1 on failure.
 */
static int tighten_bbox_vertex_enum(
    const alea_system_t* sys,
    alea_node_id_t root,
    alea_bbox_t* out
) {
    halfspace_t planes[MAX_LP_PLANES];
    int n = collect_plane_constraints(sys, root, planes, MAX_LP_PLANES, 0);
    if (n < 3) return -1;  /* Need at least 3 planes for a vertex */

    alea_bbox_t result = alea_bbox_empty();
    int n_feasible = 0;

    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
            for (int k = j + 1; k < n; k++) {
                /* Build 3x3 system: row r = plane normal, RHS = d */
                mat3 A;
                /* Column-major: A.m[col][row] */
                A.m[0][0] = planes[i].a; A.m[1][0] = planes[i].b; A.m[2][0] = planes[i].c;
                A.m[0][1] = planes[j].a; A.m[1][1] = planes[j].b; A.m[2][1] = planes[j].c;
                A.m[0][2] = planes[k].a; A.m[1][2] = planes[k].b; A.m[2][2] = planes[k].c;

                vec3 rhs = { planes[i].d, planes[j].d, planes[k].d };
                vec3 pt;
                if (mat3_solve_cramer(A, rhs, &pt) != 0) continue;

                if (point_feasible(planes, n, pt)) {
                    /* Expand result bbox */
                    if (pt.x < result.min_x) result.min_x = pt.x;
                    if (pt.x > result.max_x) result.max_x = pt.x;
                    if (pt.y < result.min_y) result.min_y = pt.y;
                    if (pt.y > result.max_y) result.max_y = pt.y;
                    if (pt.z < result.min_z) result.min_z = pt.z;
                    if (pt.z > result.max_z) result.max_z = pt.z;
                    n_feasible++;
                }
            }
        }
    }

    if (n_feasible == 0) return -1;

    /* Add small margin for numerical safety */
    double margin = 1.0;
    result.min_x -= margin;
    result.max_x += margin;
    result.min_y -= margin;
    result.max_y += margin;
    result.min_z -= margin;
    result.max_z += margin;

    *out = result;
    return 0;
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

    /* Strategy 1: LP vertex enumeration (fast, exact for plane-only) */
    if (tighten_bbox_vertex_enum(sys, root, &coarse) == 0) {
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
