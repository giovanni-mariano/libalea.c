// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

#include "alea_system.h"
#include <math.h>
#include <string.h>

// ============================================================================
// FUZZY NUMBER EQUALITY
// ============================================================================

bool alea_numbers_equal(double a, double b, const alea_config_t* tol) {
    // Exact equality (includes NaN, Inf)
    if (a == b) return true;
    
    // Both near zero
    if (fabs(a) < tol->zero_threshold && fabs(b) < tol->zero_threshold) {
        return true;
    }
    
    // Absolute difference
    double diff = fabs(a - b);
    if (diff <= tol->abs_tol) {
        return true;
    }
    
    // Relative difference
    double max_abs = fmax(fabs(a), fabs(b));
    if (diff <= max_abs * tol->rel_tol) {
        return true;
    }
    
    return false;
}

// ============================================================================
// PRIMITIVE CANONICALIZATION
// ============================================================================

void alea_canonicalize_primitive(alea_primitive_type_t type, 
                                alea_primitive_data_t* data,
                                int8_t* inverted) {
    
    *inverted = 0;
    switch (type) {
        case ALEA_PRIMITIVE_PLANE: {
            alea_plane_data_t* plane = &data->plane;
            
            // Normalize the normal vector
            double norm = sqrt(plane->a * plane->a + 
                             plane->b * plane->b + 
                             plane->c * plane->c);
            
            if (norm > 1e-10) {
                plane->a /= norm;
                plane->b /= norm;
                plane->c /= norm;
                plane->d /= norm;
            }
            
            // Ensure first non-zero coefficient is positive
            if (plane->a < -1e-10 || 
                (fabs(plane->a) < 1e-10 && plane->b < -1e-10) ||
                (fabs(plane->a) < 1e-10 && fabs(plane->b) < 1e-10 && plane->c < -1e-10)) {
                plane->a = -plane->a;
                plane->b = -plane->b;
                plane->c = -plane->c;
                plane->d = -plane->d;
                *inverted = 1;
            }
            break;
        }
        
        case ALEA_PRIMITIVE_SPHERE: {
            break;
        }
        
        case ALEA_PRIMITIVE_CYLINDER_X:
        case ALEA_PRIMITIVE_CYLINDER_Y:
        case ALEA_PRIMITIVE_CYLINDER_Z: {
            break;
        }
        
        case ALEA_PRIMITIVE_CONE_X:
        case ALEA_PRIMITIVE_CONE_Y:
        case ALEA_PRIMITIVE_CONE_Z: {
            break;
        }
        
        case ALEA_PRIMITIVE_TORUS_X:
        case ALEA_PRIMITIVE_TORUS_Y:
        case ALEA_PRIMITIVE_TORUS_Z: {
            break;
        }

        case ALEA_PRIMITIVE_QUADRIC: {
            break;
        }
        
        case ALEA_PRIMITIVE_RPP: {
            alea_box_data_t* box = &data->box;
            if (box->min_x > box->max_x) {
                double tmp = box->min_x;
                box->min_x = box->max_x;
                box->max_x = tmp;
            }
            if (box->min_y > box->max_y) {
                double tmp = box->min_y;
                box->min_y = box->max_y;
                box->max_y = tmp;
            }
            if (box->min_z > box->max_z) {
                double tmp = box->min_z;
                box->min_z = box->max_z;
                box->max_z = tmp;
            }
            break;
        }

        case ALEA_PRIMITIVE_RCC: {
            // RCC is already canonical: base + height vector + radius
            // No normalization needed - the height vector defines both direction and length
            break;
        }

        case ALEA_PRIMITIVE_BOX: {
            // General box: corner + 3 edge vectors, no canonicalization needed
            break;
        }

        case ALEA_PRIMITIVE_SPH: {
            // SPH macrobody: same as sphere, no canonicalization needed
            break;
        }

        case ALEA_PRIMITIVE_TRC: {
            // TRC: base + height + two radii, no canonicalization needed
            break;
        }

        case ALEA_PRIMITIVE_ELL: {
            // Ellipsoid: two foci + major axis, no canonicalization needed
            break;
        }

        case ALEA_PRIMITIVE_REC: {
            // Right elliptical cylinder: base + height + 2 semi-axis vectors
            break;
        }

        case ALEA_PRIMITIVE_WED: {
            // Wedge: vertex + 3 edge vectors
            break;
        }

        case ALEA_PRIMITIVE_RHP: {
            // Right hexagonal prism: base + height + 3 vertex vectors
            break;
        }

        case ALEA_PRIMITIVE_ARB: {
            // Arbitrary polyhedron: up to 8 corners + 6 faces
            break;
        }
    }
}

// ============================================================================
// PRIMITIVE EQUALITY
// ============================================================================

bool alea_primitives_equal(alea_primitive_type_t type1, const alea_primitive_data_t* data1,
                         alea_primitive_type_t type2, const alea_primitive_data_t* data2,
                         const alea_config_t* tol,
                         int8_t* match_inverted) {
    if (match_inverted) *match_inverted = 0;
    if (type1 != type2) return false;

    switch (type1) {
        case ALEA_PRIMITIVE_PLANE: {
            const alea_plane_data_t* p1 = &data1->plane;
            const alea_plane_data_t* p2 = &data2->plane;

            // Check if normals are parallel
            double dot = p1->a * p2->a + p1->b * p2->b + p1->c * p2->c;

            if (alea_numbers_equal(fabs(dot), 1.0, tol)) {
                if (dot > 0) {
                    // Same direction
                    return alea_numbers_equal(p1->d, p2->d, tol);
                } else {
                    // Opposite direction — signal inversion needed
                    if (alea_numbers_equal(p1->d, -p2->d, tol)) {
                        if (match_inverted) *match_inverted = 1;
                        return true;
                    }
                    return false;
                }
            }
            return false;
        }
        
        case ALEA_PRIMITIVE_SPHERE: {
            const alea_sphere_data_t* s1 = &data1->sphere;
            const alea_sphere_data_t* s2 = &data2->sphere;
            
            return alea_numbers_equal(s1->center_x, s2->center_x, tol) &&
                   alea_numbers_equal(s1->center_y, s2->center_y, tol) &&
                   alea_numbers_equal(s1->center_z, s2->center_z, tol) &&
                   alea_numbers_equal(s1->radius, s2->radius, tol);
        }
        
        case ALEA_PRIMITIVE_CYLINDER_Z: {
            const alea_cylinder_z_data_t* c1 = &data1->cyl_z;
            const alea_cylinder_z_data_t* c2 = &data2->cyl_z;
            
            return alea_numbers_equal(c1->center_x, c2->center_x, tol) &&
                   alea_numbers_equal(c1->center_y, c2->center_y, tol) &&
                   alea_numbers_equal(c1->radius, c2->radius, tol);
        }

        case ALEA_PRIMITIVE_CYLINDER_X: {
            const alea_cylinder_x_data_t* c1 = &data1->cyl_x;
            const alea_cylinder_x_data_t* c2 = &data2->cyl_x;

            return alea_numbers_equal(c1->center_y, c2->center_y, tol) &&
                   alea_numbers_equal(c1->center_z, c2->center_z, tol) &&
                   alea_numbers_equal(c1->radius, c2->radius, tol);
        }

        case ALEA_PRIMITIVE_CYLINDER_Y: {
            const alea_cylinder_y_data_t* c1 = &data1->cyl_y;
            const alea_cylinder_y_data_t* c2 = &data2->cyl_y;

            return alea_numbers_equal(c1->center_x, c2->center_x, tol) &&
                   alea_numbers_equal(c1->center_z, c2->center_z, tol) &&
                   alea_numbers_equal(c1->radius, c2->radius, tol);
        }
        
        case ALEA_PRIMITIVE_CONE_X: {
            const alea_cone_x_data_t* c1 = &data1->cone_x;
            const alea_cone_x_data_t* c2 = &data2->cone_x;

            return c1->sheet_selection == c2->sheet_selection &&
                   alea_numbers_equal(c1->apex_x, c2->apex_x, tol) &&
                   alea_numbers_equal(c1->apex_y, c2->apex_y, tol) &&
                   alea_numbers_equal(c1->apex_z, c2->apex_z, tol) &&
                   alea_numbers_equal(c1->tan_angle_sq, c2->tan_angle_sq, tol);
        }

        case ALEA_PRIMITIVE_CONE_Y: {
            const alea_cone_y_data_t* c1 = &data1->cone_y;
            const alea_cone_y_data_t* c2 = &data2->cone_y;

            return c1->sheet_selection == c2->sheet_selection &&
                   alea_numbers_equal(c1->apex_x, c2->apex_x, tol) &&
                   alea_numbers_equal(c1->apex_y, c2->apex_y, tol) &&
                   alea_numbers_equal(c1->apex_z, c2->apex_z, tol) &&
                   alea_numbers_equal(c1->tan_angle_sq, c2->tan_angle_sq, tol);
        }

        case ALEA_PRIMITIVE_CONE_Z: {
            const alea_cone_z_data_t* c1 = &data1->cone_z;
            const alea_cone_z_data_t* c2 = &data2->cone_z;

            return c1->sheet_selection == c2->sheet_selection &&
                   alea_numbers_equal(c1->apex_x, c2->apex_x, tol) &&
                   alea_numbers_equal(c1->apex_y, c2->apex_y, tol) &&
                   alea_numbers_equal(c1->apex_z, c2->apex_z, tol) &&
                   alea_numbers_equal(c1->tan_angle_sq, c2->tan_angle_sq, tol);
        }
        
        case ALEA_PRIMITIVE_TORUS_X:
        case ALEA_PRIMITIVE_TORUS_Y:
        case ALEA_PRIMITIVE_TORUS_Z: {
            const alea_torus_data_t* t1 = &data1->torus;
            const alea_torus_data_t* t2 = &data2->torus;
            
            return t1->axis == t2->axis &&
                   alea_numbers_equal(t1->center_x, t2->center_x, tol) &&
                   alea_numbers_equal(t1->center_y, t2->center_y, tol) &&
                   alea_numbers_equal(t1->center_z, t2->center_z, tol) &&
                   alea_numbers_equal(t1->major_radius, t2->major_radius, tol) &&
                   alea_numbers_equal(t1->minor_radius, t2->minor_radius, tol) &&
                   alea_numbers_equal(t1->axial_semiwidth_B, t2->axial_semiwidth_B, tol);
        }
        
        case ALEA_PRIMITIVE_RPP: {
            const alea_box_data_t* b1 = &data1->box;
            const alea_box_data_t* b2 = &data2->box;
            
            return alea_numbers_equal(b1->min_x, b2->min_x, tol) &&
                   alea_numbers_equal(b1->max_x, b2->max_x, tol) &&
                   alea_numbers_equal(b1->min_y, b2->min_y, tol) &&
                   alea_numbers_equal(b1->max_y, b2->max_y, tol) &&
                   alea_numbers_equal(b1->min_z, b2->min_z, tol) &&
                   alea_numbers_equal(b1->max_z, b2->max_z, tol);
        }
        
        case ALEA_PRIMITIVE_QUADRIC: {
            const alea_quadric_data_t* q1 = &data1->quadric;
            const alea_quadric_data_t* q2 = &data2->quadric;
            
            // Compare all 10 coefficients
            for (int i = 0; i < 10; i++) {
                if (!alea_numbers_equal(q1->coeffs[i], q2->coeffs[i], tol)) {
                    return false;
                }
            }
            return true;
        }

        case ALEA_PRIMITIVE_RCC: {
            const alea_rcc_data_t* r1 = &data1->rcc;
            const alea_rcc_data_t* r2 = &data2->rcc;

            return alea_numbers_equal(r1->base_x, r2->base_x, tol) &&
                   alea_numbers_equal(r1->base_y, r2->base_y, tol) &&
                   alea_numbers_equal(r1->base_z, r2->base_z, tol) &&
                   alea_numbers_equal(r1->height_x, r2->height_x, tol) &&
                   alea_numbers_equal(r1->height_y, r2->height_y, tol) &&
                   alea_numbers_equal(r1->height_z, r2->height_z, tol) &&
                   alea_numbers_equal(r1->radius, r2->radius, tol);
        }

        case ALEA_PRIMITIVE_BOX: {
            const alea_box_general_data_t* b1 = &data1->box_general;
            const alea_box_general_data_t* b2 = &data2->box_general;

            return alea_numbers_equal(b1->corner_x, b2->corner_x, tol) &&
                   alea_numbers_equal(b1->corner_y, b2->corner_y, tol) &&
                   alea_numbers_equal(b1->corner_z, b2->corner_z, tol) &&
                   alea_numbers_equal(b1->v1_x, b2->v1_x, tol) &&
                   alea_numbers_equal(b1->v1_y, b2->v1_y, tol) &&
                   alea_numbers_equal(b1->v1_z, b2->v1_z, tol) &&
                   alea_numbers_equal(b1->v2_x, b2->v2_x, tol) &&
                   alea_numbers_equal(b1->v2_y, b2->v2_y, tol) &&
                   alea_numbers_equal(b1->v2_z, b2->v2_z, tol) &&
                   alea_numbers_equal(b1->v3_x, b2->v3_x, tol) &&
                   alea_numbers_equal(b1->v3_y, b2->v3_y, tol) &&
                   alea_numbers_equal(b1->v3_z, b2->v3_z, tol);
        }

        case ALEA_PRIMITIVE_SPH: {
            const alea_sph_data_t* s1 = &data1->sph;
            const alea_sph_data_t* s2 = &data2->sph;

            return alea_numbers_equal(s1->center_x, s2->center_x, tol) &&
                   alea_numbers_equal(s1->center_y, s2->center_y, tol) &&
                   alea_numbers_equal(s1->center_z, s2->center_z, tol) &&
                   alea_numbers_equal(s1->radius, s2->radius, tol);
        }

        case ALEA_PRIMITIVE_TRC: {
            const alea_trc_data_t* t1 = &data1->trc;
            const alea_trc_data_t* t2 = &data2->trc;

            return alea_numbers_equal(t1->base_x, t2->base_x, tol) &&
                   alea_numbers_equal(t1->base_y, t2->base_y, tol) &&
                   alea_numbers_equal(t1->base_z, t2->base_z, tol) &&
                   alea_numbers_equal(t1->height_x, t2->height_x, tol) &&
                   alea_numbers_equal(t1->height_y, t2->height_y, tol) &&
                   alea_numbers_equal(t1->height_z, t2->height_z, tol) &&
                   alea_numbers_equal(t1->base_radius, t2->base_radius, tol) &&
                   alea_numbers_equal(t1->top_radius, t2->top_radius, tol);
        }

        case ALEA_PRIMITIVE_ELL: {
            const alea_ell_data_t* e1 = &data1->ell;
            const alea_ell_data_t* e2 = &data2->ell;

            return alea_numbers_equal(e1->v1_x, e2->v1_x, tol) &&
                   alea_numbers_equal(e1->v1_y, e2->v1_y, tol) &&
                   alea_numbers_equal(e1->v1_z, e2->v1_z, tol) &&
                   alea_numbers_equal(e1->v2_x, e2->v2_x, tol) &&
                   alea_numbers_equal(e1->v2_y, e2->v2_y, tol) &&
                   alea_numbers_equal(e1->v2_z, e2->v2_z, tol) &&
                   alea_numbers_equal(e1->major_axis_len, e2->major_axis_len, tol);
        }

        case ALEA_PRIMITIVE_REC: {
            const alea_rec_data_t* r1 = &data1->rec;
            const alea_rec_data_t* r2 = &data2->rec;

            return alea_numbers_equal(r1->base_x, r2->base_x, tol) &&
                   alea_numbers_equal(r1->base_y, r2->base_y, tol) &&
                   alea_numbers_equal(r1->base_z, r2->base_z, tol) &&
                   alea_numbers_equal(r1->height_x, r2->height_x, tol) &&
                   alea_numbers_equal(r1->height_y, r2->height_y, tol) &&
                   alea_numbers_equal(r1->height_z, r2->height_z, tol) &&
                   alea_numbers_equal(r1->axis1_x, r2->axis1_x, tol) &&
                   alea_numbers_equal(r1->axis1_y, r2->axis1_y, tol) &&
                   alea_numbers_equal(r1->axis1_z, r2->axis1_z, tol) &&
                   alea_numbers_equal(r1->axis2_x, r2->axis2_x, tol) &&
                   alea_numbers_equal(r1->axis2_y, r2->axis2_y, tol) &&
                   alea_numbers_equal(r1->axis2_z, r2->axis2_z, tol);
        }

        case ALEA_PRIMITIVE_WED: {
            const alea_wed_data_t* w1 = &data1->wed;
            const alea_wed_data_t* w2 = &data2->wed;

            return alea_numbers_equal(w1->vertex_x, w2->vertex_x, tol) &&
                   alea_numbers_equal(w1->vertex_y, w2->vertex_y, tol) &&
                   alea_numbers_equal(w1->vertex_z, w2->vertex_z, tol) &&
                   alea_numbers_equal(w1->v1_x, w2->v1_x, tol) &&
                   alea_numbers_equal(w1->v1_y, w2->v1_y, tol) &&
                   alea_numbers_equal(w1->v1_z, w2->v1_z, tol) &&
                   alea_numbers_equal(w1->v2_x, w2->v2_x, tol) &&
                   alea_numbers_equal(w1->v2_y, w2->v2_y, tol) &&
                   alea_numbers_equal(w1->v2_z, w2->v2_z, tol) &&
                   alea_numbers_equal(w1->v3_x, w2->v3_x, tol) &&
                   alea_numbers_equal(w1->v3_y, w2->v3_y, tol) &&
                   alea_numbers_equal(w1->v3_z, w2->v3_z, tol);
        }

        case ALEA_PRIMITIVE_RHP: {
            const alea_rhp_data_t* h1 = &data1->rhp;
            const alea_rhp_data_t* h2 = &data2->rhp;

            return alea_numbers_equal(h1->base_x, h2->base_x, tol) &&
                   alea_numbers_equal(h1->base_y, h2->base_y, tol) &&
                   alea_numbers_equal(h1->base_z, h2->base_z, tol) &&
                   alea_numbers_equal(h1->height_x, h2->height_x, tol) &&
                   alea_numbers_equal(h1->height_y, h2->height_y, tol) &&
                   alea_numbers_equal(h1->height_z, h2->height_z, tol) &&
                   alea_numbers_equal(h1->r1_x, h2->r1_x, tol) &&
                   alea_numbers_equal(h1->r1_y, h2->r1_y, tol) &&
                   alea_numbers_equal(h1->r1_z, h2->r1_z, tol) &&
                   alea_numbers_equal(h1->r2_x, h2->r2_x, tol) &&
                   alea_numbers_equal(h1->r2_y, h2->r2_y, tol) &&
                   alea_numbers_equal(h1->r2_z, h2->r2_z, tol) &&
                   alea_numbers_equal(h1->r3_x, h2->r3_x, tol) &&
                   alea_numbers_equal(h1->r3_y, h2->r3_y, tol) &&
                   alea_numbers_equal(h1->r3_z, h2->r3_z, tol);
        }

        case ALEA_PRIMITIVE_ARB: {
            const alea_arb_data_t* a1 = &data1->arb;
            const alea_arb_data_t* a2 = &data2->arb;

            if (a1->num_corners != a2->num_corners ||
                a1->num_faces != a2->num_faces) {
                return false;
            }

            // Compare corners
            for (int i = 0; i < a1->num_corners; i++) {
                if (!alea_numbers_equal(a1->corners[i][0], a2->corners[i][0], tol) ||
                    !alea_numbers_equal(a1->corners[i][1], a2->corners[i][1], tol) ||
                    !alea_numbers_equal(a1->corners[i][2], a2->corners[i][2], tol)) {
                    return false;
                }
            }

            // Compare faces (exact integer match)
            for (int i = 0; i < a1->num_faces; i++) {
                for (int j = 0; j < 4; j++) {
                    if (a1->faces[i][j] != a2->faces[i][j]) {
                        return false;
                    }
                }
            }
            return true;
        }
    }

    return false;
}