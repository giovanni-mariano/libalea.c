// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/**
 * @file primitive_desc.c
 * @brief Primitive descriptor table and operations
 *
 * All primitive-specific code lives here. Each primitive type has:
 * - An eval function (signed distance)
 * - A bbox function (bounding box)
 *
 * The descriptor table maps types to these functions.
 */

#include "primitive_desc.h"
#include "util/math.h"
#include <math.h>
#include <float.h>
#include <stddef.h>

#define BBOX_LARGE 1.0e10

/* Forward declarations for bbox functions used by interval functions */
static alea_bbox_t bbox_box_general(const alea_primitive_data_t* data);
static alea_bbox_t bbox_wed(const alea_primitive_data_t* data);
static alea_bbox_t bbox_arb(const alea_primitive_data_t* data);

/* ============================================================================
 * EVAL FUNCTIONS
 * ============================================================================ */

static double eval_plane(const alea_primitive_data_t* data, double x, double y, double z) {
    const alea_plane_data_t* p = &data->plane;
    return p->a * x + p->b * y + p->c * z + p->d;
}

static double eval_sphere(const alea_primitive_data_t* data, double x, double y, double z) {
    const alea_sphere_data_t* s = &data->sphere;
    double dx = x - s->center_x;
    double dy = y - s->center_y;
    double dz = z - s->center_z;
    return dx*dx + dy*dy + dz*dz - s->radius*s->radius;
}

static double eval_cylinder_x(const alea_primitive_data_t* data, double x, double y, double z) {
    const alea_cylinder_x_data_t* c = &data->cyl_x;
    (void)x;
    double dy = y - c->center_y;
    double dz = z - c->center_z;
    return dy*dy + dz*dz - c->radius*c->radius;
}

static double eval_cylinder_y(const alea_primitive_data_t* data, double x, double y, double z) {
    const alea_cylinder_y_data_t* c = &data->cyl_y;
    (void)y;
    double dx = x - c->center_x;
    double dz = z - c->center_z;
    return dx*dx + dz*dz - c->radius*c->radius;
}

static double eval_cylinder_z(const alea_primitive_data_t* data, double x, double y, double z) {
    const alea_cylinder_z_data_t* c = &data->cyl_z;
    (void)z;
    double dx = x - c->center_x;
    double dy = y - c->center_y;
    return dx*dx + dy*dy - c->radius*c->radius;
}

static double eval_cone_x(const alea_primitive_data_t* data, double x, double y, double z) {
    const alea_cone_x_data_t* c = &data->cone_x;
    double dx = x - c->apex_x;
    double dy = y - c->apex_y;
    double dz = z - c->apex_z;
    double rho_sq = dy*dy + dz*dz;
    double f = rho_sq - c->tan_angle_sq * dx*dx;
    if (c->sheet_selection != 0) {
        if ((c->sheet_selection > 0 && dx < 0) || (c->sheet_selection < 0 && dx > 0))
            return 1.0;
    }
    return f;
}

static double eval_cone_y(const alea_primitive_data_t* data, double x, double y, double z) {
    const alea_cone_y_data_t* c = &data->cone_y;
    double dx = x - c->apex_x;
    double dy = y - c->apex_y;
    double dz = z - c->apex_z;
    double rho_sq = dx*dx + dz*dz;
    double f = rho_sq - c->tan_angle_sq * dy*dy;
    if (c->sheet_selection != 0) {
        if ((c->sheet_selection > 0 && dy < 0) || (c->sheet_selection < 0 && dy > 0))
            return 1.0;
    }
    return f;
}

static double eval_cone_z(const alea_primitive_data_t* data, double x, double y, double z) {
    const alea_cone_z_data_t* c = &data->cone_z;
    double dx = x - c->apex_x;
    double dy = y - c->apex_y;
    double dz = z - c->apex_z;
    double rho_sq = dx*dx + dy*dy;
    double f = rho_sq - c->tan_angle_sq * dz*dz;
    if (c->sheet_selection != 0) {
        if ((c->sheet_selection > 0 && dz < 0) || (c->sheet_selection < 0 && dz > 0))
            return 1.0;
    }
    return f;
}

static double eval_box(const alea_primitive_data_t* data, double x, double y, double z) {
    const alea_box_data_t* b = &data->box;
    double dx = MAX(b->min_x - x, x - b->max_x);
    double dy = MAX(b->min_y - y, y - b->max_y);
    double dz = MAX(b->min_z - z, z - b->max_z);
    return MAX(dx, MAX(dy, dz));
}

static double eval_quadric(const alea_primitive_data_t* data, double x, double y, double z) {
    const double* c = data->quadric.coeffs;
    return c[0]*x*x + c[1]*y*y + c[2]*z*z +
           c[3]*x*y + c[4]*y*z + c[5]*x*z +
           c[6]*x + c[7]*y + c[8]*z + c[9];
}

static double eval_torus_x(const alea_primitive_data_t* data, double x, double y, double z) {
    const alea_torus_data_t* t = &data->torus;
    double dx = x - t->center_x;
    double dy = y - t->center_y;
    double dz = z - t->center_z;
    double rho = sqrt(dy*dy + dz*dz);
    double C = t->minor_radius;
    if (C == 0.0) return 1.0;
    double radial = (rho - t->major_radius);
    double term_radial = (radial * radial) / (C * C);
    double term_axial = 0.0;
    if (t->axial_semiwidth_B > 0.0)
        term_axial = (dx * dx) / (t->axial_semiwidth_B * t->axial_semiwidth_B);
    return term_axial + term_radial - 1.0;
}

static double eval_torus_y(const alea_primitive_data_t* data, double x, double y, double z) {
    const alea_torus_data_t* t = &data->torus;
    double dx = x - t->center_x;
    double dy = y - t->center_y;
    double dz = z - t->center_z;
    double rho = sqrt(dx*dx + dz*dz);
    double C = t->minor_radius;
    if (C == 0.0) return 1.0;
    double radial = (rho - t->major_radius);
    double term_radial = (radial * radial) / (C * C);
    double term_axial = 0.0;
    if (t->axial_semiwidth_B > 0.0)
        term_axial = (dy * dy) / (t->axial_semiwidth_B * t->axial_semiwidth_B);
    return term_axial + term_radial - 1.0;
}

static double eval_torus_z(const alea_primitive_data_t* data, double x, double y, double z) {
    const alea_torus_data_t* t = &data->torus;
    double dx = x - t->center_x;
    double dy = y - t->center_y;
    double dz = z - t->center_z;
    double rho = sqrt(dx*dx + dy*dy);
    double C = t->minor_radius;
    if (C == 0.0) return 1.0;
    double radial = (rho - t->major_radius);
    double term_radial = (radial * radial) / (C * C);
    double term_axial = 0.0;
    if (t->axial_semiwidth_B > 0.0)
        term_axial = (dz * dz) / (t->axial_semiwidth_B * t->axial_semiwidth_B);
    return term_axial + term_radial - 1.0;
}


static double eval_rcc(const alea_primitive_data_t* data, double x, double y, double z) {
      const alea_rcc_data_t* rcc = &data->rcc;
      double px = x - rcc->base_x;
      double py = y - rcc->base_y;
      double pz = z - rcc->base_z;
      double h_len_sq = rcc->height_x*rcc->height_x + rcc->height_y*rcc->height_y + rcc->height_z*rcc->height_z;
      if (h_len_sq < 1e-20) return 1.0;

      double h_len = sqrt(h_len_sq);
      double dot = px*rcc->height_x + py*rcc->height_y + pz*rcc->height_z;
      double t = dot / h_len_sq;

      // Clamp t to [0,1] for finding closest point on finite cylinder axis
      double t_clamped = (t < 0.0) ? 0.0 : (t > 1.0) ? 1.0 : t;

      double ax = rcc->base_x + t_clamped * rcc->height_x;
      double ay = rcc->base_y + t_clamped * rcc->height_y;
      double az = rcc->base_z + t_clamped * rcc->height_z;
      double dx = x - ax, dy = y - ay, dz = z - az;
      double radial_dist_sq = dx*dx + dy*dy + dz*dz;
      double radial_dist = sqrt(radial_dist_sq);

      // Axial distance (negative inside, positive outside)
      double axial_dist = (t < 0.0) ? -t * h_len :
                          (t > 1.0) ? (t - 1.0) * h_len :
                          -fmin(t, 1.0 - t) * h_len;  // negative when inside

      // Radial distance (negative inside, positive outside)
      double radial_sdf = radial_dist - rcc->radius;

      // For a capped cylinder: max of radial and axial distances
      return fmax(radial_sdf, axial_dist);
  }

/* SPH macrobody - same as sphere but different storage */
static double eval_sph(const alea_primitive_data_t* data, double x, double y, double z) {
    const alea_sph_data_t* s = &data->sph;
    double dx = x - s->center_x;
    double dy = y - s->center_y;
    double dz = z - s->center_z;
    return dx*dx + dy*dy + dz*dz - s->radius*s->radius;
}

/* TRC - Truncated Right Cone (frustum) */
static double eval_trc(const alea_primitive_data_t* data, double x, double y, double z) {
    const alea_trc_data_t* trc = &data->trc;

    /* Transform point to local coords where base is at origin, axis is along +Z */
    double px = x - trc->base_x;
    double py = y - trc->base_y;
    double pz = z - trc->base_z;

    double h_len_sq = trc->height_x*trc->height_x + trc->height_y*trc->height_y + trc->height_z*trc->height_z;
    if (h_len_sq < 1e-20) return 1.0;
    double h_len = sqrt(h_len_sq);

    /* Project onto axis: t in [0, h_len] when inside */
    double t = (px*trc->height_x + py*trc->height_y + pz*trc->height_z) / h_len;

    /* Radius varies linearly from base_radius to top_radius */
    double r_base = trc->base_radius;
    double r_top = trc->top_radius;

    /* Clamp t for finding closest point */
    double t_clamped = (t < 0.0) ? 0.0 : (t > h_len) ? h_len : t;
    double t_norm = t_clamped / h_len;
    double r_at_t = r_base + (r_top - r_base) * t_norm;

    /* Find closest point on axis */
    double ax = trc->base_x + (t_clamped / h_len) * trc->height_x;
    double ay = trc->base_y + (t_clamped / h_len) * trc->height_y;
    double az = trc->base_z + (t_clamped / h_len) * trc->height_z;

    double dx = x - ax, dy = y - ay, dz = z - az;
    double radial_dist = sqrt(dx*dx + dy*dy + dz*dz);

    /* Signed radial distance */
    double radial_sdf = radial_dist - r_at_t;

    /* Axial distance */
    double axial_dist;
    if (t < 0.0) {
        axial_dist = -t;
    } else if (t > h_len) {
        axial_dist = t - h_len;
    } else {
        axial_dist = -fmin(t, h_len - t);
    }

    return fmax(radial_sdf, axial_dist);
}

/* ELL - Ellipsoid defined by two foci and major axis length */
static double eval_ell(const alea_primitive_data_t* data, double x, double y, double z) {
    const alea_ell_data_t* ell = &data->ell;

    /* Distance to both foci */
    double d1 = sqrt((x - ell->v1_x)*(x - ell->v1_x) +
                     (y - ell->v1_y)*(y - ell->v1_y) +
                     (z - ell->v1_z)*(z - ell->v1_z));
    double d2 = sqrt((x - ell->v2_x)*(x - ell->v2_x) +
                     (y - ell->v2_y)*(y - ell->v2_y) +
                     (z - ell->v2_z)*(z - ell->v2_z));

    /* Ellipsoid: sum of distances to foci = major axis length */
    return d1 + d2 - ell->major_axis_len;
}

/* BOX_GENERAL - Oriented box defined by corner and 3 edge vectors */
static double eval_box_general(const alea_primitive_data_t* data, double x, double y, double z) {
    const alea_box_general_data_t* b = &data->box_general;

    /* Transform to local coordinates */
    double px = x - b->corner_x;
    double py = y - b->corner_y;
    double pz = z - b->corner_z;

    /* Edge lengths */
    double len1 = sqrt(b->v1_x*b->v1_x + b->v1_y*b->v1_y + b->v1_z*b->v1_z);
    double len2 = sqrt(b->v2_x*b->v2_x + b->v2_y*b->v2_y + b->v2_z*b->v2_z);
    double len3 = sqrt(b->v3_x*b->v3_x + b->v3_y*b->v3_y + b->v3_z*b->v3_z);

    if (len1 < 1e-20 || len2 < 1e-20 || len3 < 1e-20) return 1.0;

    /* Unit vectors */
    double u1x = b->v1_x/len1, u1y = b->v1_y/len1, u1z = b->v1_z/len1;
    double u2x = b->v2_x/len2, u2y = b->v2_y/len2, u2z = b->v2_z/len2;
    double u3x = b->v3_x/len3, u3y = b->v3_y/len3, u3z = b->v3_z/len3;

    /* Project onto each axis */
    double t1 = px*u1x + py*u1y + pz*u1z;
    double t2 = px*u2x + py*u2y + pz*u2z;
    double t3 = px*u3x + py*u3y + pz*u3z;

    /* SDF for each axis */
    double d1 = MAX(-t1, t1 - len1);
    double d2 = MAX(-t2, t2 - len2);
    double d3 = MAX(-t3, t3 - len3);

    return MAX(d1, MAX(d2, d3));
}

/* REC - Right Elliptical Cylinder */
static double eval_rec(const alea_primitive_data_t* data, double x, double y, double z) {
    const alea_rec_data_t* rec = &data->rec;

    /* Transform to local coordinates */
    double px = x - rec->base_x;
    double py = y - rec->base_y;
    double pz = z - rec->base_z;

    /* Height vector length */
    double h_len_sq = rec->height_x*rec->height_x + rec->height_y*rec->height_y + rec->height_z*rec->height_z;
    if (h_len_sq < 1e-20) return 1.0;
    double h_len = sqrt(h_len_sq);

    /* Semi-axis lengths */
    double a_len = sqrt(rec->axis1_x*rec->axis1_x + rec->axis1_y*rec->axis1_y + rec->axis1_z*rec->axis1_z);
    double b_len = sqrt(rec->axis2_x*rec->axis2_x + rec->axis2_y*rec->axis2_y + rec->axis2_z*rec->axis2_z);
    if (a_len < 1e-20 || b_len < 1e-20) return 1.0;

    /* Project onto height axis */
    double t = (px*rec->height_x + py*rec->height_y + pz*rec->height_z) / h_len_sq;

    /* Clamp for closest point */
    double t_clamped = (t < 0.0) ? 0.0 : (t > 1.0) ? 1.0 : t;

    /* Find point on axis */
    double ax = rec->base_x + t_clamped * rec->height_x;
    double ay = rec->base_y + t_clamped * rec->height_y;
    double az = rec->base_z + t_clamped * rec->height_z;

    /* Vector from axis to point */
    double dx = x - ax, dy = y - ay, dz = z - az;

    /* Project onto ellipse axes (normalized) */
    double e1 = (dx*rec->axis1_x + dy*rec->axis1_y + dz*rec->axis1_z) / (a_len * a_len);
    double e2 = (dx*rec->axis2_x + dy*rec->axis2_y + dz*rec->axis2_z) / (b_len * b_len);

    /* Ellipse distance: (e1/a)^2 + (e2/b)^2 - 1 */
    double ellipse_dist = e1*e1 + e2*e2 - 1.0;

    /* Axial distance */
    double axial_dist;
    if (t < 0.0) {
        axial_dist = -t * h_len;
    } else if (t > 1.0) {
        axial_dist = (t - 1.0) * h_len;
    } else {
        axial_dist = -fmin(t, 1.0 - t) * h_len;
    }

    return fmax(ellipse_dist, axial_dist);
}

/* WED - Wedge (triangular prism) */
static double eval_wed(const alea_primitive_data_t* data, double x, double y, double z) {
    const alea_wed_data_t* wed = &data->wed;

    /* Transform to local coordinates */
    double px = x - wed->vertex_x;
    double py = y - wed->vertex_y;
    double pz = z - wed->vertex_z;

    /* Vector lengths */
    double len1 = sqrt(wed->v1_x*wed->v1_x + wed->v1_y*wed->v1_y + wed->v1_z*wed->v1_z);
    double len2 = sqrt(wed->v2_x*wed->v2_x + wed->v2_y*wed->v2_y + wed->v2_z*wed->v2_z);
    double len3 = sqrt(wed->v3_x*wed->v3_x + wed->v3_y*wed->v3_y + wed->v3_z*wed->v3_z);

    if (len1 < 1e-20 || len2 < 1e-20 || len3 < 1e-20) return 1.0;

    /* Unit vectors */
    double u1x = wed->v1_x/len1, u1y = wed->v1_y/len1, u1z = wed->v1_z/len1;
    double u2x = wed->v2_x/len2, u2y = wed->v2_y/len2, u2z = wed->v2_z/len2;
    double u3x = wed->v3_x/len3, u3y = wed->v3_y/len3, u3z = wed->v3_z/len3;

    /* Project onto each axis */
    double t1 = px*u1x + py*u1y + pz*u1z;
    double t2 = px*u2x + py*u2y + pz*u2z;
    double t3 = px*u3x + py*u3y + pz*u3z;

    /* Box constraints for v1, v2 directions */
    double d1 = MAX(-t1, t1 - len1);
    double d2 = MAX(-t2, t2 - len2);

    /* Height constraint (along v3) */
    double d3 = MAX(-t3, t3 - len3);

    /* Triangular constraint: t1/len1 + t2/len2 <= 1 */
    double diag = (t1/len1 + t2/len2 - 1.0);
    /* Scale to distance (approximate) */
    double diag_dist = diag * fmin(len1, len2) / sqrt(2.0);

    return MAX(MAX(d1, d2), MAX(d3, diag_dist));
}

/* RHP - Right Hexagonal Prism */
static double eval_rhp(const alea_primitive_data_t* data, double x, double y, double z) {
    const alea_rhp_data_t* rhp = &data->rhp;

    /* Transform to local coordinates relative to base center */
    double px = x - rhp->base_x;
    double py = y - rhp->base_y;
    double pz = z - rhp->base_z;

    /* Height vector */
    double h_len_sq = rhp->height_x*rhp->height_x + rhp->height_y*rhp->height_y + rhp->height_z*rhp->height_z;
    if (h_len_sq < 1e-20) return 1.0;
    double h_len = sqrt(h_len_sq);

    /* Project onto height axis */
    double t = (px*rhp->height_x + py*rhp->height_y + pz*rhp->height_z) / h_len;
    double axial_dist = MAX(-t, t - h_len);

    /* Remove axial component to get position in hexagonal cross-section */
    double hx = rhp->height_x / h_len;
    double hy = rhp->height_y / h_len;
    double hz = rhp->height_z / h_len;
    double qx = px - t * hx;
    double qy = py - t * hy;
    double qz = pz - t * hz;

    /* Distances to each pair of hex faces (normal is along r_i direction) */
    double r1_len = sqrt(rhp->r1_x*rhp->r1_x + rhp->r1_y*rhp->r1_y + rhp->r1_z*rhp->r1_z);
    double r2_len = sqrt(rhp->r2_x*rhp->r2_x + rhp->r2_y*rhp->r2_y + rhp->r2_z*rhp->r2_z);
    double r3_len = sqrt(rhp->r3_x*rhp->r3_x + rhp->r3_y*rhp->r3_y + rhp->r3_z*rhp->r3_z);

    if (r1_len < 1e-20 || r2_len < 1e-20 || r3_len < 1e-20) return 1.0;

    /* Project onto each direction and check bounds [-len, +len] */
    double p1 = qx*rhp->r1_x/r1_len + qy*rhp->r1_y/r1_len + qz*rhp->r1_z/r1_len;
    double p2 = qx*rhp->r2_x/r2_len + qy*rhp->r2_y/r2_len + qz*rhp->r2_z/r2_len;
    double p3 = qx*rhp->r3_x/r3_len + qy*rhp->r3_y/r3_len + qz*rhp->r3_z/r3_len;

    double d1 = fabs(p1) - r1_len;
    double d2 = fabs(p2) - r2_len;
    double d3 = fabs(p3) - r3_len;

    double hex_dist = MAX(d1, MAX(d2, d3));

    return MAX(axial_dist, hex_dist);
}

/* ARB - Arbitrary Polyhedron */
static double eval_arb(const alea_primitive_data_t* data, double x, double y, double z) {
    const alea_arb_data_t* arb = &data->arb;

    /* For each face, compute signed distance to plane */
    double max_dist = -1e30;

    for (int f = 0; f < arb->num_faces; f++) {
        int i0 = arb->faces[f][0] - 1;  /* MCNP uses 1-based indices */
        int i1 = arb->faces[f][1] - 1;
        int i2 = arb->faces[f][2] - 1;

        if (i0 < 0 || i1 < 0 || i2 < 0 ||
            i0 >= arb->num_corners || i1 >= arb->num_corners || i2 >= arb->num_corners) {
            continue;
        }

        /* Three vertices of face */
        double x0 = arb->corners[i0][0], y0 = arb->corners[i0][1], z0 = arb->corners[i0][2];
        double x1 = arb->corners[i1][0], y1 = arb->corners[i1][1], z1 = arb->corners[i1][2];
        double x2 = arb->corners[i2][0], y2 = arb->corners[i2][1], z2 = arb->corners[i2][2];

        /* Edge vectors */
        double e1x = x1 - x0, e1y = y1 - y0, e1z = z1 - z0;
        double e2x = x2 - x0, e2y = y2 - y0, e2z = z2 - z0;

        /* Normal (cross product, pointing inward for MCNP convention) */
        double nx = e1y*e2z - e1z*e2y;
        double ny = e1z*e2x - e1x*e2z;
        double nz = e1x*e2y - e1y*e2x;

        double nlen = sqrt(nx*nx + ny*ny + nz*nz);
        if (nlen < 1e-20) continue;

        nx /= nlen; ny /= nlen; nz /= nlen;

        /* Signed distance to plane (negative inside) */
        double dist = (x - x0)*nx + (y - y0)*ny + (z - z0)*nz;
        if (dist > max_dist) max_dist = dist;
    }

    return max_dist;
}

/* ============================================================================
 * INTERVAL EVAL FUNCTIONS (for geometric pruning)
 * ============================================================================ */

static alea_interval_t interval_plane(const alea_primitive_data_t* data, const alea_bbox_t* box) {
    const alea_plane_data_t* p = &data->plane;
    /* f(x,y,z) = ax + by + cz + d */
    alea_interval_t ix = alea_iv_make(box->min_x, box->max_x);
    alea_interval_t iy = alea_iv_make(box->min_y, box->max_y);
    alea_interval_t iz = alea_iv_make(box->min_z, box->max_z);

    alea_interval_t result = alea_iv_make(p->d, p->d);
    result = alea_iv_add(result, alea_iv_mul_scalar(ix, p->a));
    result = alea_iv_add(result, alea_iv_mul_scalar(iy, p->b));
    result = alea_iv_add(result, alea_iv_mul_scalar(iz, p->c));
    return result;
}

static alea_interval_t interval_sphere(const alea_primitive_data_t* data, const alea_bbox_t* box) {
    const alea_sphere_data_t* s = &data->sphere;
    /* f = (x-cx)^2 + (y-cy)^2 + (z-cz)^2 - R^2 */
    alea_interval_t dx = alea_iv_sub(alea_iv_make(box->min_x, box->max_x),
                                    alea_iv_make(s->center_x, s->center_x));
    alea_interval_t dy = alea_iv_sub(alea_iv_make(box->min_y, box->max_y),
                                    alea_iv_make(s->center_y, s->center_y));
    alea_interval_t dz = alea_iv_sub(alea_iv_make(box->min_z, box->max_z),
                                    alea_iv_make(s->center_z, s->center_z));

    alea_interval_t d2 = alea_iv_add(alea_iv_add(alea_iv_sqr(dx), alea_iv_sqr(dy)), alea_iv_sqr(dz));
    double R2 = s->radius * s->radius;
    return alea_iv_sub(d2, alea_iv_make(R2, R2));
}

static alea_interval_t interval_cylinder_x(const alea_primitive_data_t* data, const alea_bbox_t* box) {
    const alea_cylinder_x_data_t* c = &data->cyl_x;
    /* For infinite cylinder along X: f = (y-cy)^2 + (z-cz)^2 - R^2 */
    alea_interval_t dy = alea_iv_sub(alea_iv_make(box->min_y, box->max_y),
                                    alea_iv_make(c->center_y, c->center_y));
    alea_interval_t dz = alea_iv_sub(alea_iv_make(box->min_z, box->max_z),
                                    alea_iv_make(c->center_z, c->center_z));

    alea_interval_t d2 = alea_iv_add(alea_iv_sqr(dy), alea_iv_sqr(dz));
    double R2 = c->radius * c->radius;
    return alea_iv_sub(d2, alea_iv_make(R2, R2));
}

static alea_interval_t interval_cylinder_y(const alea_primitive_data_t* data, const alea_bbox_t* box) {
    const alea_cylinder_y_data_t* c = &data->cyl_y;
    /* For infinite cylinder along Y: f = (x-cx)^2 + (z-cz)^2 - R^2 */
    alea_interval_t dx = alea_iv_sub(alea_iv_make(box->min_x, box->max_x),
                                    alea_iv_make(c->center_x, c->center_x));
    alea_interval_t dz = alea_iv_sub(alea_iv_make(box->min_z, box->max_z),
                                    alea_iv_make(c->center_z, c->center_z));

    alea_interval_t d2 = alea_iv_add(alea_iv_sqr(dx), alea_iv_sqr(dz));
    double R2 = c->radius * c->radius;
    return alea_iv_sub(d2, alea_iv_make(R2, R2));
}

static alea_interval_t interval_cylinder_z(const alea_primitive_data_t* data, const alea_bbox_t* box) {
    const alea_cylinder_z_data_t* c = &data->cyl_z;
    /* For infinite cylinder along Z: f = (x-cx)^2 + (y-cy)^2 - R^2 */
    alea_interval_t dx = alea_iv_sub(alea_iv_make(box->min_x, box->max_x),
                                    alea_iv_make(c->center_x, c->center_x));
    alea_interval_t dy = alea_iv_sub(alea_iv_make(box->min_y, box->max_y),
                                    alea_iv_make(c->center_y, c->center_y));

    alea_interval_t d2 = alea_iv_add(alea_iv_sqr(dx), alea_iv_sqr(dy));
    double R2 = c->radius * c->radius;
    return alea_iv_sub(d2, alea_iv_make(R2, R2));
}

static alea_interval_t interval_cone_x(const alea_primitive_data_t* data, const alea_bbox_t* box) {
    const alea_cone_x_data_t* c = &data->cone_x;
    /* f = (y-ay)^2 + (z-az)^2 - t^2*(x-ax)^2 */
    alea_interval_t dx = alea_iv_sub(alea_iv_make(box->min_x, box->max_x),
                                    alea_iv_make(c->apex_x, c->apex_x));
    alea_interval_t dy = alea_iv_sub(alea_iv_make(box->min_y, box->max_y),
                                    alea_iv_make(c->apex_y, c->apex_y));
    alea_interval_t dz = alea_iv_sub(alea_iv_make(box->min_z, box->max_z),
                                    alea_iv_make(c->apex_z, c->apex_z));

    alea_interval_t rho_sq = alea_iv_add(alea_iv_sqr(dy), alea_iv_sqr(dz));
    alea_interval_t axial_sq = alea_iv_mul_scalar(alea_iv_sqr(dx), c->tan_angle_sq);
    return alea_iv_sub(rho_sq, axial_sq);
}

static alea_interval_t interval_cone_y(const alea_primitive_data_t* data, const alea_bbox_t* box) {
    const alea_cone_y_data_t* c = &data->cone_y;
    /* f = (x-ax)^2 + (z-az)^2 - t^2*(y-ay)^2 */
    alea_interval_t dx = alea_iv_sub(alea_iv_make(box->min_x, box->max_x),
                                    alea_iv_make(c->apex_x, c->apex_x));
    alea_interval_t dy = alea_iv_sub(alea_iv_make(box->min_y, box->max_y),
                                    alea_iv_make(c->apex_y, c->apex_y));
    alea_interval_t dz = alea_iv_sub(alea_iv_make(box->min_z, box->max_z),
                                    alea_iv_make(c->apex_z, c->apex_z));

    alea_interval_t rho_sq = alea_iv_add(alea_iv_sqr(dx), alea_iv_sqr(dz));
    alea_interval_t axial_sq = alea_iv_mul_scalar(alea_iv_sqr(dy), c->tan_angle_sq);
    return alea_iv_sub(rho_sq, axial_sq);
}

static alea_interval_t interval_cone_z(const alea_primitive_data_t* data, const alea_bbox_t* box) {
    const alea_cone_z_data_t* c = &data->cone_z;
    /* f = (x-ax)^2 + (y-ay)^2 - t^2*(z-az)^2 */
    alea_interval_t dx = alea_iv_sub(alea_iv_make(box->min_x, box->max_x),
                                    alea_iv_make(c->apex_x, c->apex_x));
    alea_interval_t dy = alea_iv_sub(alea_iv_make(box->min_y, box->max_y),
                                    alea_iv_make(c->apex_y, c->apex_y));
    alea_interval_t dz = alea_iv_sub(alea_iv_make(box->min_z, box->max_z),
                                    alea_iv_make(c->apex_z, c->apex_z));

    alea_interval_t rho_sq = alea_iv_add(alea_iv_sqr(dx), alea_iv_sqr(dy));
    alea_interval_t axial_sq = alea_iv_mul_scalar(alea_iv_sqr(dz), c->tan_angle_sq);
    return alea_iv_sub(rho_sq, axial_sq);
}

static alea_interval_t interval_box(const alea_primitive_data_t* data, const alea_bbox_t* box) {
    const alea_box_data_t* b = &data->box;
    /* Box SDF: max(dx, dy, dz) where d = max(min - pos, pos - max) */

    /* Distance from query box to primitive box boundaries */
    double dx_min = MAX(b->min_x - box->max_x, box->min_x - b->max_x);
    double dx_max = MAX(fabs(box->max_x - b->min_x), fabs(box->min_x - b->max_x));
    double dy_min = MAX(b->min_y - box->max_y, box->min_y - b->max_y);
    double dy_max = MAX(fabs(box->max_y - b->min_y), fabs(box->min_y - b->max_y));
    double dz_min = MAX(b->min_z - box->max_z, box->min_z - b->max_z);
    double dz_max = MAX(fabs(box->max_z - b->min_z), fabs(box->min_z - b->max_z));

    /* Check if boxes overlap */
    if (box->min_x <= b->max_x && box->max_x >= b->min_x &&
        box->min_y <= b->max_y && box->max_y >= b->min_y &&
        box->min_z <= b->max_z && box->max_z >= b->min_z) {
        /* Overlapping: result includes negative values */
        double min_dist = -MIN(MIN(b->max_x - box->min_x, box->max_x - b->min_x),
                               MIN(MIN(b->max_y - box->min_y, box->max_y - b->min_y),
                                   MIN(b->max_z - box->min_z, box->max_z - b->min_z)));
        return alea_iv_make(min_dist, MAX(dx_max, MAX(dy_max, dz_max)));
    }

    /* Non-overlapping: positive distance */
    double pos_dist = MAX(dx_min, MAX(dy_min, dz_min));
    return alea_iv_make(pos_dist, MAX(dx_max, MAX(dy_max, dz_max)));
}

static alea_interval_t interval_quadric(const alea_primitive_data_t* data, const alea_bbox_t* box) {
    const double* c = data->quadric.coeffs;
    /* f = Ax^2 + By^2 + Cz^2 + Dxy + Eyz + Fxz + Gx + Hy + Iz + J */
    alea_interval_t ix = alea_iv_make(box->min_x, box->max_x);
    alea_interval_t iy = alea_iv_make(box->min_y, box->max_y);
    alea_interval_t iz = alea_iv_make(box->min_z, box->max_z);

    alea_interval_t result = alea_iv_make(c[9], c[9]); /* J */

    result = alea_iv_add(result, alea_iv_mul_scalar(alea_iv_sqr(ix), c[0])); /* Ax^2 */
    result = alea_iv_add(result, alea_iv_mul_scalar(alea_iv_sqr(iy), c[1])); /* By^2 */
    result = alea_iv_add(result, alea_iv_mul_scalar(alea_iv_sqr(iz), c[2])); /* Cz^2 */

    result = alea_iv_add(result, alea_iv_mul_scalar(ix, c[6])); /* Gx */
    result = alea_iv_add(result, alea_iv_mul_scalar(iy, c[7])); /* Hy */
    result = alea_iv_add(result, alea_iv_mul_scalar(iz, c[8])); /* Iz */

    /* Cross terms */
    if (c[3] != 0) result = alea_iv_add(result, alea_iv_mul_scalar(alea_iv_mul(ix, iy), c[3])); /* Dxy */
    if (c[4] != 0) result = alea_iv_add(result, alea_iv_mul_scalar(alea_iv_mul(iy, iz), c[4])); /* Eyz */
    if (c[5] != 0) result = alea_iv_add(result, alea_iv_mul_scalar(alea_iv_mul(iz, ix), c[5])); /* Fxz */

    return result;
}

static alea_interval_t interval_torus_x(const alea_primitive_data_t* data, const alea_bbox_t* box) {
    const alea_torus_data_t* t = &data->torus;
    /* Torus X: f = (radial/C)^2 + (axial/B)^2 - 1
     * where radial = sqrt(dy^2 + dz^2) - R, axial = dx
     * C = minor_radius, B = axial_semiwidth_B (or infinite if 0)
     */
    double R = t->major_radius;
    double C = t->minor_radius;
    if (C == 0.0) return alea_iv_make(1.0, 1.0);  /* Degenerate */

    alea_interval_t dy = alea_iv_sub(alea_iv_make(box->min_y, box->max_y),
                                    alea_iv_make(t->center_y, t->center_y));
    alea_interval_t dz = alea_iv_sub(alea_iv_make(box->min_z, box->max_z),
                                    alea_iv_make(t->center_z, t->center_z));
    alea_interval_t dx = alea_iv_sub(alea_iv_make(box->min_x, box->max_x),
                                    alea_iv_make(t->center_x, t->center_x));

    /* Radial distance in YZ plane: rho = sqrt(dy^2 + dz^2) */
    alea_interval_t rho_sq = alea_iv_add(alea_iv_sqr(dy), alea_iv_sqr(dz));
    double rho_min = sqrt(rho_sq.min);
    double rho_max = sqrt(rho_sq.max);

    /* radial = rho - R */
    double radial_min, radial_max;
    if (rho_min > R) {
        radial_min = rho_min - R;
        radial_max = rho_max - R;
    } else if (rho_max < R) {
        radial_min = R - rho_max;
        radial_max = R - rho_min;
    } else {
        radial_min = 0;
        radial_max = MAX(rho_max - R, R - rho_min);
    }

    /* Radial term: (radial/C)^2 */
    double C2 = C * C;
    alea_interval_t radial_term = alea_iv_make(radial_min * radial_min / C2,
                                              radial_max * radial_max / C2);

    /* Axial term: (dx/B)^2, or 0 if B = 0 (infinite in axial direction) */
    alea_interval_t axial_term;
    if (t->axial_semiwidth_B > 0.0) {
        double B2 = t->axial_semiwidth_B * t->axial_semiwidth_B;
        alea_interval_t dx_sq = alea_iv_sqr(dx);
        axial_term = alea_iv_make(dx_sq.min / B2, dx_sq.max / B2);
    } else {
        axial_term = alea_iv_make(0.0, 0.0);
    }

    /* f = radial_term + axial_term - 1 */
    return alea_iv_sub(alea_iv_add(radial_term, axial_term), alea_iv_make(1.0, 1.0));
}

static alea_interval_t interval_torus_y(const alea_primitive_data_t* data, const alea_bbox_t* box) {
    const alea_torus_data_t* t = &data->torus;

    /* Torus Y: f = (radial/C)^2 + (axial/B)^2 - 1
     * where radial = sqrt(dx^2 + dz^2) - R, axial = dy */
    double R = t->major_radius;
    double C = t->minor_radius;  /* radial minor radius */

    /* Compute intervals for displacements from center */
    alea_interval_t dx = alea_iv_sub(alea_iv_make(box->min_x, box->max_x),
                                    alea_iv_make(t->center_x, t->center_x));
    alea_interval_t dy = alea_iv_sub(alea_iv_make(box->min_y, box->max_y),
                                    alea_iv_make(t->center_y, t->center_y));
    alea_interval_t dz = alea_iv_sub(alea_iv_make(box->min_z, box->max_z),
                                    alea_iv_make(t->center_z, t->center_z));

    /* rho = sqrt(dx^2 + dz^2) - distance from Y axis */
    alea_interval_t rho_sq = alea_iv_add(alea_iv_sqr(dx), alea_iv_sqr(dz));
    double rho_min = sqrt(rho_sq.min);
    double rho_max = sqrt(rho_sq.max);

    /* radial = rho - R, handling the case where R is inside [rho_min, rho_max] */
    double radial_min, radial_max;
    if (rho_min > R) {
        radial_min = rho_min - R;
        radial_max = rho_max - R;
    } else if (rho_max < R) {
        radial_min = R - rho_max;
        radial_max = R - rho_min;
    } else {
        radial_min = 0;
        radial_max = MAX(rho_max - R, R - rho_min);
    }

    /* Radial term: (radial/C)^2 */
    double C2 = C * C;
    alea_interval_t radial_term = alea_iv_make(radial_min * radial_min / C2,
                                              radial_max * radial_max / C2);

    /* Axial term: (dy/B)^2, or 0 if B = 0 (infinite in axial direction) */
    alea_interval_t axial_term;
    if (t->axial_semiwidth_B > 0.0) {
        double B2 = t->axial_semiwidth_B * t->axial_semiwidth_B;
        alea_interval_t dy_sq = alea_iv_sqr(dy);
        axial_term = alea_iv_make(dy_sq.min / B2, dy_sq.max / B2);
    } else {
        axial_term = alea_iv_make(0.0, 0.0);
    }

    /* f = radial_term + axial_term - 1 */
    return alea_iv_sub(alea_iv_add(radial_term, axial_term), alea_iv_make(1.0, 1.0));
}

static alea_interval_t interval_torus_z(const alea_primitive_data_t* data, const alea_bbox_t* box) {
    const alea_torus_data_t* t = &data->torus;

    /* Torus Z: f = (radial/C)^2 + (axial/B)^2 - 1
     * where radial = sqrt(dx^2 + dy^2) - R, axial = dz */
    double R = t->major_radius;
    double C = t->minor_radius;  /* radial minor radius */

    /* Compute intervals for displacements from center */
    alea_interval_t dx = alea_iv_sub(alea_iv_make(box->min_x, box->max_x),
                                    alea_iv_make(t->center_x, t->center_x));
    alea_interval_t dy = alea_iv_sub(alea_iv_make(box->min_y, box->max_y),
                                    alea_iv_make(t->center_y, t->center_y));
    alea_interval_t dz = alea_iv_sub(alea_iv_make(box->min_z, box->max_z),
                                    alea_iv_make(t->center_z, t->center_z));

    /* rho = sqrt(dx^2 + dy^2) - distance from Z axis */
    alea_interval_t rho_sq = alea_iv_add(alea_iv_sqr(dx), alea_iv_sqr(dy));
    double rho_min = sqrt(rho_sq.min);
    double rho_max = sqrt(rho_sq.max);

    /* radial = rho - R, handling the case where R is inside [rho_min, rho_max] */
    double radial_min, radial_max;
    if (rho_min > R) {
        radial_min = rho_min - R;
        radial_max = rho_max - R;
    } else if (rho_max < R) {
        radial_min = R - rho_max;
        radial_max = R - rho_min;
    } else {
        radial_min = 0;
        radial_max = MAX(rho_max - R, R - rho_min);
    }

    /* Radial term: (radial/C)^2 */
    double C2 = C * C;
    alea_interval_t radial_term = alea_iv_make(radial_min * radial_min / C2,
                                              radial_max * radial_max / C2);

    /* Axial term: (dz/B)^2, or 0 if B = 0 (infinite in axial direction) */
    alea_interval_t axial_term;
    if (t->axial_semiwidth_B > 0.0) {
        double B2 = t->axial_semiwidth_B * t->axial_semiwidth_B;
        alea_interval_t dz_sq = alea_iv_sqr(dz);
        axial_term = alea_iv_make(dz_sq.min / B2, dz_sq.max / B2);
    } else {
        axial_term = alea_iv_make(0.0, 0.0);
    }

    /* f = radial_term + axial_term - 1 */
    return alea_iv_sub(alea_iv_add(radial_term, axial_term), alea_iv_make(1.0, 1.0));
}

static alea_interval_t interval_rcc(const alea_primitive_data_t* data, const alea_bbox_t* box) {
    const alea_rcc_data_t* rcc = &data->rcc;
    /* RCC is complex (finite cylinder with arbitrary axis) */
    /* Use conservative spherical bound: sphere at center with max extent */
    double cx = rcc->base_x + rcc->height_x * 0.5;
    double cy = rcc->base_y + rcc->height_y * 0.5;
    double cz = rcc->base_z + rcc->height_z * 0.5;
    double h_len = sqrt(rcc->height_x*rcc->height_x + rcc->height_y*rcc->height_y + rcc->height_z*rcc->height_z);
    double bounding_r = sqrt(rcc->radius*rcc->radius + (h_len*0.5)*(h_len*0.5));

    alea_interval_t dx = alea_iv_sub(alea_iv_make(box->min_x, box->max_x), alea_iv_make(cx, cx));
    alea_interval_t dy = alea_iv_sub(alea_iv_make(box->min_y, box->max_y), alea_iv_make(cy, cy));
    alea_interval_t dz = alea_iv_sub(alea_iv_make(box->min_z, box->max_z), alea_iv_make(cz, cz));

    alea_interval_t d2 = alea_iv_add(alea_iv_add(alea_iv_sqr(dx), alea_iv_sqr(dy)), alea_iv_sqr(dz));
    return alea_iv_sub(d2, alea_iv_make(bounding_r*bounding_r, bounding_r*bounding_r));
}

/* SPH interval - same as sphere */
static alea_interval_t interval_sph(const alea_primitive_data_t* data, const alea_bbox_t* box) {
    const alea_sph_data_t* s = &data->sph;
    alea_interval_t dx = alea_iv_sub(alea_iv_make(box->min_x, box->max_x), alea_iv_make(s->center_x, s->center_x));
    alea_interval_t dy = alea_iv_sub(alea_iv_make(box->min_y, box->max_y), alea_iv_make(s->center_y, s->center_y));
    alea_interval_t dz = alea_iv_sub(alea_iv_make(box->min_z, box->max_z), alea_iv_make(s->center_z, s->center_z));
    alea_interval_t d2 = alea_iv_add(alea_iv_add(alea_iv_sqr(dx), alea_iv_sqr(dy)), alea_iv_sqr(dz));
    double R2 = s->radius * s->radius;
    return alea_iv_sub(d2, alea_iv_make(R2, R2));
}

/* TRC interval - use conservative spherical bound */
static alea_interval_t interval_trc(const alea_primitive_data_t* data, const alea_bbox_t* box) {
    const alea_trc_data_t* trc = &data->trc;
    double cx = trc->base_x + trc->height_x * 0.5;
    double cy = trc->base_y + trc->height_y * 0.5;
    double cz = trc->base_z + trc->height_z * 0.5;
    double h_len = sqrt(trc->height_x*trc->height_x + trc->height_y*trc->height_y + trc->height_z*trc->height_z);
    double r_max = MAX(trc->base_radius, trc->top_radius);
    double bounding_r = sqrt(r_max*r_max + (h_len*0.5)*(h_len*0.5));

    alea_interval_t dx = alea_iv_sub(alea_iv_make(box->min_x, box->max_x), alea_iv_make(cx, cx));
    alea_interval_t dy = alea_iv_sub(alea_iv_make(box->min_y, box->max_y), alea_iv_make(cy, cy));
    alea_interval_t dz = alea_iv_sub(alea_iv_make(box->min_z, box->max_z), alea_iv_make(cz, cz));
    alea_interval_t d2 = alea_iv_add(alea_iv_add(alea_iv_sqr(dx), alea_iv_sqr(dy)), alea_iv_sqr(dz));
    return alea_iv_sub(d2, alea_iv_make(bounding_r*bounding_r, bounding_r*bounding_r));
}

/* ELL interval - use bounding sphere */
static alea_interval_t interval_ell(const alea_primitive_data_t* data, const alea_bbox_t* box) {
    const alea_ell_data_t* ell = &data->ell;
    double cx = (ell->v1_x + ell->v2_x) * 0.5;
    double cy = (ell->v1_y + ell->v2_y) * 0.5;
    double cz = (ell->v1_z + ell->v2_z) * 0.5;
    double a = ell->major_axis_len * 0.5;

    alea_interval_t dx = alea_iv_sub(alea_iv_make(box->min_x, box->max_x), alea_iv_make(cx, cx));
    alea_interval_t dy = alea_iv_sub(alea_iv_make(box->min_y, box->max_y), alea_iv_make(cy, cy));
    alea_interval_t dz = alea_iv_sub(alea_iv_make(box->min_z, box->max_z), alea_iv_make(cz, cz));
    alea_interval_t d2 = alea_iv_add(alea_iv_add(alea_iv_sqr(dx), alea_iv_sqr(dy)), alea_iv_sqr(dz));
    return alea_iv_sub(d2, alea_iv_make(a*a, a*a));
}

/* BOX_GENERAL interval - conservative box bound */
static alea_interval_t interval_box_general(const alea_primitive_data_t* data, const alea_bbox_t* box) {
    alea_bbox_t prim_bbox = bbox_box_general(data);
    /* Same logic as interval_box */
    double dx_min = MAX(prim_bbox.min_x - box->max_x, box->min_x - prim_bbox.max_x);
    double dx_max = MAX(fabs(box->max_x - prim_bbox.min_x), fabs(box->min_x - prim_bbox.max_x));
    double dy_min = MAX(prim_bbox.min_y - box->max_y, box->min_y - prim_bbox.max_y);
    double dy_max = MAX(fabs(box->max_y - prim_bbox.min_y), fabs(box->min_y - prim_bbox.max_y));
    double dz_min = MAX(prim_bbox.min_z - box->max_z, box->min_z - prim_bbox.max_z);
    double dz_max = MAX(fabs(box->max_z - prim_bbox.min_z), fabs(box->min_z - prim_bbox.max_z));

    if (box->min_x <= prim_bbox.max_x && box->max_x >= prim_bbox.min_x &&
        box->min_y <= prim_bbox.max_y && box->max_y >= prim_bbox.min_y &&
        box->min_z <= prim_bbox.max_z && box->max_z >= prim_bbox.min_z) {
        double min_dist = -MIN(MIN(prim_bbox.max_x - box->min_x, box->max_x - prim_bbox.min_x),
                               MIN(MIN(prim_bbox.max_y - box->min_y, box->max_y - prim_bbox.min_y),
                                   MIN(prim_bbox.max_z - box->min_z, box->max_z - prim_bbox.min_z)));
        return alea_iv_make(min_dist, MAX(dx_max, MAX(dy_max, dz_max)));
    }
    double pos_dist = MAX(dx_min, MAX(dy_min, dz_min));
    return alea_iv_make(pos_dist, MAX(dx_max, MAX(dy_max, dz_max)));
}

/* REC interval - conservative cylindrical bound */
static alea_interval_t interval_rec(const alea_primitive_data_t* data, const alea_bbox_t* box) {
    const alea_rec_data_t* rec = &data->rec;
    double a_len = sqrt(rec->axis1_x*rec->axis1_x + rec->axis1_y*rec->axis1_y + rec->axis1_z*rec->axis1_z);
    double b_len = sqrt(rec->axis2_x*rec->axis2_x + rec->axis2_y*rec->axis2_y + rec->axis2_z*rec->axis2_z);
    double r_max = MAX(a_len, b_len);
    double cx = rec->base_x + rec->height_x * 0.5;
    double cy = rec->base_y + rec->height_y * 0.5;
    double cz = rec->base_z + rec->height_z * 0.5;
    double h_len = sqrt(rec->height_x*rec->height_x + rec->height_y*rec->height_y + rec->height_z*rec->height_z);
    double bounding_r = sqrt(r_max*r_max + (h_len*0.5)*(h_len*0.5));

    alea_interval_t dx = alea_iv_sub(alea_iv_make(box->min_x, box->max_x), alea_iv_make(cx, cx));
    alea_interval_t dy = alea_iv_sub(alea_iv_make(box->min_y, box->max_y), alea_iv_make(cy, cy));
    alea_interval_t dz = alea_iv_sub(alea_iv_make(box->min_z, box->max_z), alea_iv_make(cz, cz));
    alea_interval_t d2 = alea_iv_add(alea_iv_add(alea_iv_sqr(dx), alea_iv_sqr(dy)), alea_iv_sqr(dz));
    return alea_iv_sub(d2, alea_iv_make(bounding_r*bounding_r, bounding_r*bounding_r));
}

/* WED interval - use bbox-based bound */
static alea_interval_t interval_wed(const alea_primitive_data_t* data, const alea_bbox_t* box) {
    alea_bbox_t prim_bbox = bbox_wed(data);
    double dx_min = MAX(prim_bbox.min_x - box->max_x, box->min_x - prim_bbox.max_x);
    double dx_max = MAX(fabs(box->max_x - prim_bbox.min_x), fabs(box->min_x - prim_bbox.max_x));
    double dy_min = MAX(prim_bbox.min_y - box->max_y, box->min_y - prim_bbox.max_y);
    double dy_max = MAX(fabs(box->max_y - prim_bbox.min_y), fabs(box->min_y - prim_bbox.max_y));
    double dz_min = MAX(prim_bbox.min_z - box->max_z, box->min_z - prim_bbox.max_z);
    double dz_max = MAX(fabs(box->max_z - prim_bbox.min_z), fabs(box->min_z - prim_bbox.max_z));

    if (box->min_x <= prim_bbox.max_x && box->max_x >= prim_bbox.min_x &&
        box->min_y <= prim_bbox.max_y && box->max_y >= prim_bbox.min_y &&
        box->min_z <= prim_bbox.max_z && box->max_z >= prim_bbox.min_z) {
        double min_dist = -MIN(MIN(prim_bbox.max_x - box->min_x, box->max_x - prim_bbox.min_x),
                               MIN(MIN(prim_bbox.max_y - box->min_y, box->max_y - prim_bbox.min_y),
                                   MIN(prim_bbox.max_z - box->min_z, box->max_z - prim_bbox.min_z)));
        return alea_iv_make(min_dist, MAX(dx_max, MAX(dy_max, dz_max)));
    }
    double pos_dist = MAX(dx_min, MAX(dy_min, dz_min));
    return alea_iv_make(pos_dist, MAX(dx_max, MAX(dy_max, dz_max)));
}

/* RHP interval - use bounding cylinder */
static alea_interval_t interval_rhp(const alea_primitive_data_t* data, const alea_bbox_t* box) {
    const alea_rhp_data_t* rhp = &data->rhp;
    double r1 = sqrt(rhp->r1_x*rhp->r1_x + rhp->r1_y*rhp->r1_y + rhp->r1_z*rhp->r1_z);
    double r2 = sqrt(rhp->r2_x*rhp->r2_x + rhp->r2_y*rhp->r2_y + rhp->r2_z*rhp->r2_z);
    double r3 = sqrt(rhp->r3_x*rhp->r3_x + rhp->r3_y*rhp->r3_y + rhp->r3_z*rhp->r3_z);
    double r_max = MAX(r1, MAX(r2, r3));
    double cx = rhp->base_x + rhp->height_x * 0.5;
    double cy = rhp->base_y + rhp->height_y * 0.5;
    double cz = rhp->base_z + rhp->height_z * 0.5;
    double h_len = sqrt(rhp->height_x*rhp->height_x + rhp->height_y*rhp->height_y + rhp->height_z*rhp->height_z);
    double bounding_r = sqrt(r_max*r_max + (h_len*0.5)*(h_len*0.5));

    alea_interval_t dx = alea_iv_sub(alea_iv_make(box->min_x, box->max_x), alea_iv_make(cx, cx));
    alea_interval_t dy = alea_iv_sub(alea_iv_make(box->min_y, box->max_y), alea_iv_make(cy, cy));
    alea_interval_t dz = alea_iv_sub(alea_iv_make(box->min_z, box->max_z), alea_iv_make(cz, cz));
    alea_interval_t d2 = alea_iv_add(alea_iv_add(alea_iv_sqr(dx), alea_iv_sqr(dy)), alea_iv_sqr(dz));
    return alea_iv_sub(d2, alea_iv_make(bounding_r*bounding_r, bounding_r*bounding_r));
}

/* ARB interval - use bbox-based bound */
static alea_interval_t interval_arb(const alea_primitive_data_t* data, const alea_bbox_t* box) {
    alea_bbox_t prim_bbox = bbox_arb(data);
    double dx_min = MAX(prim_bbox.min_x - box->max_x, box->min_x - prim_bbox.max_x);
    double dx_max = MAX(fabs(box->max_x - prim_bbox.min_x), fabs(box->min_x - prim_bbox.max_x));
    double dy_min = MAX(prim_bbox.min_y - box->max_y, box->min_y - prim_bbox.max_y);
    double dy_max = MAX(fabs(box->max_y - prim_bbox.min_y), fabs(box->min_y - prim_bbox.max_y));
    double dz_min = MAX(prim_bbox.min_z - box->max_z, box->min_z - prim_bbox.max_z);
    double dz_max = MAX(fabs(box->max_z - prim_bbox.min_z), fabs(box->min_z - prim_bbox.max_z));

    if (box->min_x <= prim_bbox.max_x && box->max_x >= prim_bbox.min_x &&
        box->min_y <= prim_bbox.max_y && box->max_y >= prim_bbox.min_y &&
        box->min_z <= prim_bbox.max_z && box->max_z >= prim_bbox.min_z) {
        double min_dist = -MIN(MIN(prim_bbox.max_x - box->min_x, box->max_x - prim_bbox.min_x),
                               MIN(MIN(prim_bbox.max_y - box->min_y, box->max_y - prim_bbox.min_y),
                                   MIN(prim_bbox.max_z - box->min_z, box->max_z - prim_bbox.min_z)));
        return alea_iv_make(min_dist, MAX(dx_max, MAX(dy_max, dz_max)));
    }
    double pos_dist = MAX(dx_min, MAX(dy_min, dz_min));
    return alea_iv_make(pos_dist, MAX(dx_max, MAX(dy_max, dz_max)));
}

/* ============================================================================
 * BBOX FUNCTIONS
 * ============================================================================ */

static alea_bbox_t bbox_plane(const alea_primitive_data_t* data) {
    (void)data;  /* unused */
    /* A plane is a 2D surface with no volume - return infinite bbox.
     *
     * Note: Returning thin slabs at the plane position causes problems
     * when intersecting planes at different positions (results in empty bbox
     * because min > max after intersection).
     *
     * With infinite bbox:
     * - Cells made only of planes = infinite bbox (correct, less efficient)
     * - Cells with planes + bounded primitives = bounded primitive's bbox
     */
    return (alea_bbox_t){-BBOX_LARGE, BBOX_LARGE, -BBOX_LARGE, BBOX_LARGE, -BBOX_LARGE, BBOX_LARGE};
}

static alea_bbox_t bbox_sphere(const alea_primitive_data_t* data) {
    const alea_sphere_data_t* s = &data->sphere;
    return (alea_bbox_t){
        s->center_x - s->radius, s->center_x + s->radius,
        s->center_y - s->radius, s->center_y + s->radius,
        s->center_z - s->radius, s->center_z + s->radius
    };
}

static alea_bbox_t bbox_cylinder_x(const alea_primitive_data_t* data) {
    const alea_cylinder_x_data_t* c = &data->cyl_x;
    return (alea_bbox_t){
        -INFINITY, INFINITY,
        c->center_y - c->radius, c->center_y + c->radius,
        c->center_z - c->radius, c->center_z + c->radius
    };
}

static alea_bbox_t bbox_cylinder_y(const alea_primitive_data_t* data) {
    const alea_cylinder_y_data_t* c = &data->cyl_y;
    return (alea_bbox_t){
        c->center_x - c->radius, c->center_x + c->radius,
        -INFINITY, INFINITY,
        c->center_z - c->radius, c->center_z + c->radius
    };
}

static alea_bbox_t bbox_cylinder_z(const alea_primitive_data_t* data) {
    const alea_cylinder_z_data_t* c = &data->cyl_z;
    return (alea_bbox_t){
        c->center_x - c->radius, c->center_x + c->radius,
        c->center_y - c->radius, c->center_y + c->radius,
        -INFINITY, INFINITY
    };
}

static alea_bbox_t bbox_cone_x(const alea_primitive_data_t* data) {
    (void)data;
    return (alea_bbox_t){ -INFINITY, INFINITY, -INFINITY, INFINITY, -INFINITY, INFINITY };
}

static alea_bbox_t bbox_cone_y(const alea_primitive_data_t* data) {
    (void)data;
    return (alea_bbox_t){ -INFINITY, INFINITY, -INFINITY, INFINITY, -INFINITY, INFINITY };
}

static alea_bbox_t bbox_cone_z(const alea_primitive_data_t* data) {
    (void)data;
    return (alea_bbox_t){ -INFINITY, INFINITY, -INFINITY, INFINITY, -INFINITY, INFINITY };
}

static alea_bbox_t bbox_box(const alea_primitive_data_t* data) {
    const alea_box_data_t* b = &data->box;
    return (alea_bbox_t){b->min_x, b->max_x, b->min_y, b->max_y, b->min_z, b->max_z};
}

static alea_bbox_t bbox_quadric(const alea_primitive_data_t* data) {
    /*
     * General quadric: Ax² + By² + Cz² + Dxy + Eyz + Fxz + Gx + Hy + Iz + J = 0
     * coeffs: [A, B, C, D, E, F, G, H, I, J]
     *
     * Strategy:
     * 1. For bounded quadrics (ellipsoids), compute actual bounds
     * 2. For semi-bounded (cylinders, cones), bound in finite directions
     * 3. Fall back to BBOX_LARGE for unbounded directions
     *
     * At an extreme point in x-direction: ∂f/∂y = 0, ∂f/∂z = 0, f = 0
     * This gives: 2By + Dx + Ez + H = 0 and 2Cz + Ey + Fx + I = 0
     * Solve for y, z in terms of x, substitute into f = 0.
     */
    const double* c = data->quadric.coeffs;
    double A = c[0], B = c[1], C = c[2];
    double D = c[3], E = c[4], F = c[5];
    double G = c[6], H = c[7], I = c[8], J = c[9];

    alea_bbox_t bbox = {-BBOX_LARGE, BBOX_LARGE, -BBOX_LARGE, BBOX_LARGE, -BBOX_LARGE, BBOX_LARGE};

    /* Determinant of the 2x2 system for y,z given x */
    double det_yz = 4.0*B*C - E*E;
    double det_xz = 4.0*A*C - F*F;
    double det_xy = 4.0*A*B - D*D;

    /* Find x-extent: solve for y,z from ∂f/∂y = ∂f/∂z = 0, then substitute */
    if (fabs(det_yz) > 1e-12 && A > 1e-12) {
        /* Can solve for y,z: y = (E*I - 2*C*H - (E*F - 2*C*D)*x) / det_yz
         *                    z = (E*H - 2*B*I - (E*D - 2*B*F)*x) / det_yz
         * Substitute into f = 0 to get quadratic in x: a*x² + b*x + c = 0
         */
        double coef_y_x = (E*F - 2.0*C*D) / det_yz;
        double coef_y_c = (E*I - 2.0*C*H) / det_yz;
        double coef_z_x = (E*D - 2.0*B*F) / det_yz;
        double coef_z_c = (E*H - 2.0*B*I) / det_yz;

        /* Substitute y = coef_y_x*x + coef_y_c, z = coef_z_x*x + coef_z_c into f */
        double qa = A + B*coef_y_x*coef_y_x + C*coef_z_x*coef_z_x
                  + D*coef_y_x + E*coef_y_x*coef_z_x + F*coef_z_x;
        double qb = 2.0*B*coef_y_x*coef_y_c + 2.0*C*coef_z_x*coef_z_c
                  + D*coef_y_c + E*(coef_y_x*coef_z_c + coef_y_c*coef_z_x) + F*coef_z_c
                  + G + H*coef_y_x + I*coef_z_x;
        double qc = B*coef_y_c*coef_y_c + C*coef_z_c*coef_z_c
                  + E*coef_y_c*coef_z_c + H*coef_y_c + I*coef_z_c + J;

        if (qa > 1e-12) {
            double disc = qb*qb - 4.0*qa*qc;
            if (disc >= 0) {
                double sqrt_disc = sqrt(disc);
                bbox.min_x = (-qb - sqrt_disc) / (2.0*qa);
                bbox.max_x = (-qb + sqrt_disc) / (2.0*qa);
            }
        }
    }

    /* Find y-extent similarly */
    if (fabs(det_xz) > 1e-12 && B > 1e-12) {
        double coef_x_y = (E*F - 2.0*C*D) / det_xz;
        double coef_x_c = (F*I - 2.0*C*G) / det_xz;
        double coef_z_y = (F*D - 2.0*A*E) / det_xz;
        double coef_z_c = (F*G - 2.0*A*I) / det_xz;

        double qa = B + A*coef_x_y*coef_x_y + C*coef_z_y*coef_z_y
                  + D*coef_x_y + E*coef_z_y + F*coef_x_y*coef_z_y;
        double qb = 2.0*A*coef_x_y*coef_x_c + 2.0*C*coef_z_y*coef_z_c
                  + D*coef_x_c + E*coef_z_c + F*(coef_x_y*coef_z_c + coef_x_c*coef_z_y)
                  + H + G*coef_x_y + I*coef_z_y;
        double qc = A*coef_x_c*coef_x_c + C*coef_z_c*coef_z_c
                  + F*coef_x_c*coef_z_c + G*coef_x_c + I*coef_z_c + J;

        if (qa > 1e-12) {
            double disc = qb*qb - 4.0*qa*qc;
            if (disc >= 0) {
                double sqrt_disc = sqrt(disc);
                bbox.min_y = (-qb - sqrt_disc) / (2.0*qa);
                bbox.max_y = (-qb + sqrt_disc) / (2.0*qa);
            }
        }
    }

    /* Find z-extent similarly */
    if (fabs(det_xy) > 1e-12 && C > 1e-12) {
        double coef_x_z = (D*E - 2.0*B*F) / det_xy;
        double coef_x_c = (D*H - 2.0*B*G) / det_xy;
        double coef_y_z = (D*F - 2.0*A*E) / det_xy;
        double coef_y_c = (D*G - 2.0*A*H) / det_xy;

        double qa = C + A*coef_x_z*coef_x_z + B*coef_y_z*coef_y_z
                  + D*coef_x_z*coef_y_z + E*coef_y_z + F*coef_x_z;
        double qb = 2.0*A*coef_x_z*coef_x_c + 2.0*B*coef_y_z*coef_y_c
                  + D*(coef_x_z*coef_y_c + coef_x_c*coef_y_z) + E*coef_y_c + F*coef_x_c
                  + I + G*coef_x_z + H*coef_y_z;
        double qc = A*coef_x_c*coef_x_c + B*coef_y_c*coef_y_c
                  + D*coef_x_c*coef_y_c + G*coef_x_c + H*coef_y_c + J;

        if (qa > 1e-12) {
            double disc = qb*qb - 4.0*qa*qc;
            if (disc >= 0) {
                double sqrt_disc = sqrt(disc);
                bbox.min_z = (-qb - sqrt_disc) / (2.0*qa);
                bbox.max_z = (-qb + sqrt_disc) / (2.0*qa);
            }
        }
    }

    return bbox;
}

static alea_bbox_t bbox_torus_x(const alea_primitive_data_t* data) {
    const alea_torus_data_t* t = &data->torus;
    double R = t->major_radius + t->minor_radius;
    double min_x = (t->axial_semiwidth_B > 0.0) ? t->center_x - t->axial_semiwidth_B : -BBOX_LARGE;
    double max_x = (t->axial_semiwidth_B > 0.0) ? t->center_x + t->axial_semiwidth_B : BBOX_LARGE;
    return (alea_bbox_t){
        min_x, max_x,
        t->center_y - R, t->center_y + R,
        t->center_z - R, t->center_z + R
    };
}

static alea_bbox_t bbox_torus_y(const alea_primitive_data_t* data) {
    const alea_torus_data_t* t = &data->torus;
    double R = t->major_radius + t->minor_radius;
    double min_y = (t->axial_semiwidth_B > 0.0) ? t->center_y - t->axial_semiwidth_B : -BBOX_LARGE;
    double max_y = (t->axial_semiwidth_B > 0.0) ? t->center_y + t->axial_semiwidth_B : BBOX_LARGE;
    return (alea_bbox_t){
        t->center_x - R, t->center_x + R,
        min_y, max_y,
        t->center_z - R, t->center_z + R
    };
}

static alea_bbox_t bbox_torus_z(const alea_primitive_data_t* data) {
    const alea_torus_data_t* t = &data->torus;
    double R = t->major_radius + t->minor_radius;
    double min_z = (t->axial_semiwidth_B > 0.0) ? t->center_z - t->axial_semiwidth_B : -BBOX_LARGE;
    double max_z = (t->axial_semiwidth_B > 0.0) ? t->center_z + t->axial_semiwidth_B : BBOX_LARGE;
    return (alea_bbox_t){
        t->center_x - R, t->center_x + R,
        t->center_y - R, t->center_y + R,
        min_z, max_z
    };
}

static alea_bbox_t bbox_rcc(const alea_primitive_data_t* data) {
    const alea_rcc_data_t* rcc = &data->rcc;
    double bx = rcc->base_x, by = rcc->base_y, bz = rcc->base_z;
    double tx = bx + rcc->height_x, ty = by + rcc->height_y, tz = bz + rcc->height_z;
    double r = rcc->radius;
    double h_len = sqrt(rcc->height_x*rcc->height_x + rcc->height_y*rcc->height_y + rcc->height_z*rcc->height_z);
    if (h_len < 1e-20) {
        return (alea_bbox_t){bx - r, bx + r, by - r, by + r, bz - r, bz + r};
    }
    double ax = rcc->height_x / h_len, ay = rcc->height_y / h_len, az = rcc->height_z / h_len;
    double rx = r * sqrt(1.0 - ax*ax), ry = r * sqrt(1.0 - ay*ay), rz = r * sqrt(1.0 - az*az);
    return (alea_bbox_t){
        MIN(bx, tx) - rx, MAX(bx, tx) + rx,
        MIN(by, ty) - ry, MAX(by, ty) + ry,
        MIN(bz, tz) - rz, MAX(bz, tz) + rz
    };
}

static alea_bbox_t bbox_sph(const alea_primitive_data_t* data) {
    const alea_sph_data_t* s = &data->sph;
    return (alea_bbox_t){
        s->center_x - s->radius, s->center_x + s->radius,
        s->center_y - s->radius, s->center_y + s->radius,
        s->center_z - s->radius, s->center_z + s->radius
    };
}

static alea_bbox_t bbox_trc(const alea_primitive_data_t* data) {
    const alea_trc_data_t* trc = &data->trc;
    double bx = trc->base_x, by = trc->base_y, bz = trc->base_z;
    double tx = bx + trc->height_x, ty = by + trc->height_y, tz = bz + trc->height_z;
    double r_max = MAX(trc->base_radius, trc->top_radius);

    double h_len = sqrt(trc->height_x*trc->height_x + trc->height_y*trc->height_y + trc->height_z*trc->height_z);
    if (h_len < 1e-20) {
        return (alea_bbox_t){bx - r_max, bx + r_max, by - r_max, by + r_max, bz - r_max, bz + r_max};
    }

    double ax = trc->height_x / h_len, ay = trc->height_y / h_len, az = trc->height_z / h_len;
    double rx = r_max * sqrt(1.0 - ax*ax), ry = r_max * sqrt(1.0 - ay*ay), rz = r_max * sqrt(1.0 - az*az);
    return (alea_bbox_t){
        MIN(bx, tx) - rx, MAX(bx, tx) + rx,
        MIN(by, ty) - ry, MAX(by, ty) + ry,
        MIN(bz, tz) - rz, MAX(bz, tz) + rz
    };
}

static alea_bbox_t bbox_ell(const alea_primitive_data_t* data) {
    const alea_ell_data_t* ell = &data->ell;
    /* Center is midpoint of foci */
    double cx = (ell->v1_x + ell->v2_x) * 0.5;
    double cy = (ell->v1_y + ell->v2_y) * 0.5;
    double cz = (ell->v1_z + ell->v2_z) * 0.5;

    /* Semi-major axis is half the major axis length */
    double a = ell->major_axis_len * 0.5;

    /* Conservative: sphere with radius = semi-major axis */
    return (alea_bbox_t){
        cx - a, cx + a,
        cy - a, cy + a,
        cz - a, cz + a
    };
}

static alea_bbox_t bbox_box_general(const alea_primitive_data_t* data) {
    const alea_box_general_data_t* b = &data->box_general;

    /* 8 corners of the box */
    double corners[8][3];
    for (int i = 0; i < 8; i++) {
        corners[i][0] = b->corner_x;
        corners[i][1] = b->corner_y;
        corners[i][2] = b->corner_z;
        if (i & 1) { corners[i][0] += b->v1_x; corners[i][1] += b->v1_y; corners[i][2] += b->v1_z; }
        if (i & 2) { corners[i][0] += b->v2_x; corners[i][1] += b->v2_y; corners[i][2] += b->v2_z; }
        if (i & 4) { corners[i][0] += b->v3_x; corners[i][1] += b->v3_y; corners[i][2] += b->v3_z; }
    }

    alea_bbox_t bbox = {corners[0][0], corners[0][0], corners[0][1], corners[0][1], corners[0][2], corners[0][2]};
    for (int i = 1; i < 8; i++) {
        bbox.min_x = MIN(bbox.min_x, corners[i][0]);
        bbox.max_x = MAX(bbox.max_x, corners[i][0]);
        bbox.min_y = MIN(bbox.min_y, corners[i][1]);
        bbox.max_y = MAX(bbox.max_y, corners[i][1]);
        bbox.min_z = MIN(bbox.min_z, corners[i][2]);
        bbox.max_z = MAX(bbox.max_z, corners[i][2]);
    }
    return bbox;
}

static alea_bbox_t bbox_rec(const alea_primitive_data_t* data) {
    const alea_rec_data_t* rec = &data->rec;

    /* Semi-axis lengths */
    double a_len = sqrt(rec->axis1_x*rec->axis1_x + rec->axis1_y*rec->axis1_y + rec->axis1_z*rec->axis1_z);
    double b_len = sqrt(rec->axis2_x*rec->axis2_x + rec->axis2_y*rec->axis2_y + rec->axis2_z*rec->axis2_z);
    double r_max = MAX(a_len, b_len);

    /* Top center */
    double tx = rec->base_x + rec->height_x;
    double ty = rec->base_y + rec->height_y;
    double tz = rec->base_z + rec->height_z;

    return (alea_bbox_t){
        MIN(rec->base_x, tx) - r_max, MAX(rec->base_x, tx) + r_max,
        MIN(rec->base_y, ty) - r_max, MAX(rec->base_y, ty) + r_max,
        MIN(rec->base_z, tz) - r_max, MAX(rec->base_z, tz) + r_max
    };
}

static alea_bbox_t bbox_wed(const alea_primitive_data_t* data) {
    const alea_wed_data_t* wed = &data->wed;

    /* 6 vertices of the wedge */
    double verts[6][3];
    /* Base triangle */
    verts[0][0] = wed->vertex_x; verts[0][1] = wed->vertex_y; verts[0][2] = wed->vertex_z;
    verts[1][0] = wed->vertex_x + wed->v1_x; verts[1][1] = wed->vertex_y + wed->v1_y; verts[1][2] = wed->vertex_z + wed->v1_z;
    verts[2][0] = wed->vertex_x + wed->v2_x; verts[2][1] = wed->vertex_y + wed->v2_y; verts[2][2] = wed->vertex_z + wed->v2_z;
    /* Top triangle (+ v3) */
    for (int i = 0; i < 3; i++) {
        verts[i+3][0] = verts[i][0] + wed->v3_x;
        verts[i+3][1] = verts[i][1] + wed->v3_y;
        verts[i+3][2] = verts[i][2] + wed->v3_z;
    }

    alea_bbox_t bbox = {verts[0][0], verts[0][0], verts[0][1], verts[0][1], verts[0][2], verts[0][2]};
    for (int i = 1; i < 6; i++) {
        bbox.min_x = MIN(bbox.min_x, verts[i][0]);
        bbox.max_x = MAX(bbox.max_x, verts[i][0]);
        bbox.min_y = MIN(bbox.min_y, verts[i][1]);
        bbox.max_y = MAX(bbox.max_y, verts[i][1]);
        bbox.min_z = MIN(bbox.min_z, verts[i][2]);
        bbox.max_z = MAX(bbox.max_z, verts[i][2]);
    }
    return bbox;
}

static alea_bbox_t bbox_rhp(const alea_primitive_data_t* data) {
    const alea_rhp_data_t* rhp = &data->rhp;

    /* Maximum radial extent */
    double r1 = sqrt(rhp->r1_x*rhp->r1_x + rhp->r1_y*rhp->r1_y + rhp->r1_z*rhp->r1_z);
    double r2 = sqrt(rhp->r2_x*rhp->r2_x + rhp->r2_y*rhp->r2_y + rhp->r2_z*rhp->r2_z);
    double r3 = sqrt(rhp->r3_x*rhp->r3_x + rhp->r3_y*rhp->r3_y + rhp->r3_z*rhp->r3_z);
    double r_max = MAX(r1, MAX(r2, r3));

    /* Top center */
    double tx = rhp->base_x + rhp->height_x;
    double ty = rhp->base_y + rhp->height_y;
    double tz = rhp->base_z + rhp->height_z;

    return (alea_bbox_t){
        MIN(rhp->base_x, tx) - r_max, MAX(rhp->base_x, tx) + r_max,
        MIN(rhp->base_y, ty) - r_max, MAX(rhp->base_y, ty) + r_max,
        MIN(rhp->base_z, tz) - r_max, MAX(rhp->base_z, tz) + r_max
    };
}

static alea_bbox_t bbox_arb(const alea_primitive_data_t* data) {
    const alea_arb_data_t* arb = &data->arb;

    if (arb->num_corners < 1) {
        return (alea_bbox_t){-BBOX_LARGE, BBOX_LARGE, -BBOX_LARGE, BBOX_LARGE, -BBOX_LARGE, BBOX_LARGE};
    }

    alea_bbox_t bbox = {
        arb->corners[0][0], arb->corners[0][0],
        arb->corners[0][1], arb->corners[0][1],
        arb->corners[0][2], arb->corners[0][2]
    };

    for (int i = 1; i < arb->num_corners; i++) {
        bbox.min_x = MIN(bbox.min_x, arb->corners[i][0]);
        bbox.max_x = MAX(bbox.max_x, arb->corners[i][0]);
        bbox.min_y = MIN(bbox.min_y, arb->corners[i][1]);
        bbox.max_y = MAX(bbox.max_y, arb->corners[i][1]);
        bbox.min_z = MIN(bbox.min_z, arb->corners[i][2]);
        bbox.max_z = MAX(bbox.max_z, arb->corners[i][2]);
    }
    return bbox;
}

/* ============================================================================
 * TRANSFORM HELPER FUNCTIONS
 * ============================================================================ */

/*
 * Matrix layout: 3x4 row-major affine
 * mat[0..2]  = row 0: R00, R01, R02
 * mat[3]     = Tx
 * mat[4..6]  = row 1: R10, R11, R12
 * mat[7]     = Ty
 * mat[8..10] = row 2: R20, R21, R22
 * mat[11]    = Tz
 */

/* Transform a point: p' = R*p + t */
static void xform_point(const double* mat, double* x, double* y, double* z) {
    double px = *x, py = *y, pz = *z;
    *x = mat[0]*px + mat[1]*py + mat[2]*pz  + mat[3];
    *y = mat[4]*px + mat[5]*py + mat[6]*pz  + mat[7];
    *z = mat[8]*px + mat[9]*py + mat[10]*pz + mat[11];
}

/* Transform a vector (direction, no translation): v' = R*v */
static void xform_vector(const double* mat, double* vx, double* vy, double* vz) {
    double x = *vx, y = *vy, z = *vz;
    *vx = mat[0]*x + mat[1]*y + mat[2]*z;
    *vy = mat[4]*x + mat[5]*y + mat[6]*z;
    *vz = mat[8]*x + mat[9]*y + mat[10]*z;
}

/* Transform a normal: n' = (R^-T)*n, for orthogonal R this is just R*n then normalize */
static void xform_normal(const double* mat, double* nx, double* ny, double* nz) {
    xform_vector(mat, nx, ny, nz);
    double len = sqrt((*nx)*(*nx) + (*ny)*(*ny) + (*nz)*(*nz));
    if (len > 1e-15) {
        *nx /= len;
        *ny /= len;
        *nz /= len;
    }
}

/* Check if matrix is translation only (rotation part is identity) */
static bool is_translation_only(const double* mat) {
    double eps = 1e-10;
    return fabs(mat[0] - 1.0) < eps && fabs(mat[1]) < eps && fabs(mat[2]) < eps &&
           fabs(mat[4]) < eps && fabs(mat[5] - 1.0) < eps && fabs(mat[6]) < eps &&
           fabs(mat[8]) < eps && fabs(mat[9]) < eps && fabs(mat[10] - 1.0) < eps;
}

/* ============================================================================
 * TRANSFORM FUNCTIONS
 * ============================================================================ */

static bool xform_plane(const alea_primitive_data_t* in, const double* mat,
                        alea_primitive_type_t* out_type, alea_primitive_data_t* out) {
    *out_type = ALEA_PRIMITIVE_PLANE;
    out->plane = in->plane;

    /* Transform normal */
    double a = in->plane.a, b = in->plane.b, c = in->plane.c;
    xform_normal(mat, &a, &b, &c);

    /* Transform d: for plane n·x + d = 0, d_new = d_old - n_new · translation */
    double d_new = in->plane.d - (a * mat[3] + b * mat[7] + c * mat[11]);

    out->plane.a = a;
    out->plane.b = b;
    out->plane.c = c;
    out->plane.d = d_new;
    return true;
}

static bool xform_sphere(const alea_primitive_data_t* in, const double* mat,
                         alea_primitive_type_t* out_type, alea_primitive_data_t* out) {
    *out_type = ALEA_PRIMITIVE_SPHERE;
    out->sphere = in->sphere;
    xform_point(mat, &out->sphere.center_x, &out->sphere.center_y, &out->sphere.center_z);
    return true;
}

/* Helper for cylinder transforms - returns new axis-aligned type or QUADRIC */
static bool xform_cylinder_common(const double* mat,
                                  double ax, double ay, double az,
                                  double px, double py, double pz,
                                  double radius,
                                  alea_primitive_type_t* out_type,
                                  alea_primitive_data_t* out) {
    /* Transform axis direction and point on axis */
    xform_vector(mat, &ax, &ay, &az);
    xform_point(mat, &px, &py, &pz);

    /* Normalize axis - required for quadric formula */
    double axis_len = sqrt(ax*ax + ay*ay + az*az);
    if (axis_len > 1e-15) {
        ax /= axis_len;
        ay /= axis_len;
        az /= axis_len;
    }

    double eps = 1e-10;

    /* Check if still axis-aligned */
    if (fabs(ay) < eps && fabs(az) < eps && fabs(ax) > eps) {
        *out_type = ALEA_PRIMITIVE_CYLINDER_X;
        out->cyl_x.center_y = py;
        out->cyl_x.center_z = pz;
        out->cyl_x.radius = radius;
        return true;
    } else if (fabs(ax) < eps && fabs(az) < eps && fabs(ay) > eps) {
        *out_type = ALEA_PRIMITIVE_CYLINDER_Y;
        out->cyl_y.center_x = px;
        out->cyl_y.center_z = pz;
        out->cyl_y.radius = radius;
        return true;
    } else if (fabs(ax) < eps && fabs(ay) < eps && fabs(az) > eps) {
        *out_type = ALEA_PRIMITIVE_CYLINDER_Z;
        out->cyl_z.center_x = px;
        out->cyl_z.center_y = py;
        out->cyl_z.radius = radius;
        return true;
    }

    /* General orientation - convert to quadric */
    *out_type = ALEA_PRIMITIVE_QUADRIC;
    double* c = out->quadric.coeffs;
    double r = radius;

    c[0] = 1.0 - ax*ax;
    c[1] = 1.0 - ay*ay;
    c[2] = 1.0 - az*az;
    c[3] = -2.0 * ax * ay;
    c[4] = -2.0 * ay * az;
    c[5] = -2.0 * ax * az;
    c[6] = 2.0 * ((ax*ax - 1.0)*px + ax*ay*py + ax*az*pz);
    c[7] = 2.0 * (ax*ay*px + (ay*ay - 1.0)*py + ay*az*pz);
    c[8] = 2.0 * (ax*az*px + ay*az*py + (az*az - 1.0)*pz);
    c[9] = (1.0 - ax*ax)*px*px + (1.0 - ay*ay)*py*py + (1.0 - az*az)*pz*pz
         - 2.0*ax*ay*px*py - 2.0*ay*az*py*pz - 2.0*ax*az*px*pz - r*r;
    return true;
}

static bool xform_cylinder_x(const alea_primitive_data_t* in, const double* mat,
                             alea_primitive_type_t* out_type, alea_primitive_data_t* out) {
    const alea_cylinder_x_data_t* cyl = &in->cyl_x;
    if (is_translation_only(mat)) {
        *out_type = ALEA_PRIMITIVE_CYLINDER_X;
        out->cyl_x = *cyl;
        out->cyl_x.center_y += mat[7];
        out->cyl_x.center_z += mat[11];
        return true;
    }
    return xform_cylinder_common(mat, 1, 0, 0, 0, cyl->center_y, cyl->center_z,
                                 cyl->radius, out_type, out);
}

static bool xform_cylinder_y(const alea_primitive_data_t* in, const double* mat,
                             alea_primitive_type_t* out_type, alea_primitive_data_t* out) {
    const alea_cylinder_y_data_t* cyl = &in->cyl_y;
    if (is_translation_only(mat)) {
        *out_type = ALEA_PRIMITIVE_CYLINDER_Y;
        out->cyl_y = *cyl;
        out->cyl_y.center_x += mat[3];
        out->cyl_y.center_z += mat[11];
        return true;
    }
    return xform_cylinder_common(mat, 0, 1, 0, cyl->center_x, 0, cyl->center_z,
                                 cyl->radius, out_type, out);
}

static bool xform_cylinder_z(const alea_primitive_data_t* in, const double* mat,
                             alea_primitive_type_t* out_type, alea_primitive_data_t* out) {
    const alea_cylinder_z_data_t* cyl = &in->cyl_z;
    if (is_translation_only(mat)) {
        *out_type = ALEA_PRIMITIVE_CYLINDER_Z;
        out->cyl_z = *cyl;
        out->cyl_z.center_x += mat[3];
        out->cyl_z.center_y += mat[7];
        return true;
    }
    return xform_cylinder_common(mat, 0, 0, 1, cyl->center_x, cyl->center_y, 0,
                                 cyl->radius, out_type, out);
}

/* Helper for cone transforms */
static bool xform_cone_common(const double* mat,
                              double ax, double ay, double az,
                              double apex_x, double apex_y, double apex_z,
                              double tan_sq, int sheet,
                              alea_primitive_type_t* out_type,
                              alea_primitive_data_t* out) {
    xform_vector(mat, &ax, &ay, &az);
    xform_point(mat, &apex_x, &apex_y, &apex_z);

    /* Normalize axis - required for quadric formula */
    double axis_len = sqrt(ax*ax + ay*ay + az*az);
    if (axis_len > 1e-15) {
        ax /= axis_len;
        ay /= axis_len;
        az /= axis_len;
    }

    double eps = 1e-10;

    if (fabs(ay) < eps && fabs(az) < eps && fabs(ax) > eps) {
        *out_type = ALEA_PRIMITIVE_CONE_X;
        out->cone_x.apex_x = apex_x;
        out->cone_x.apex_y = apex_y;
        out->cone_x.apex_z = apex_z;
        out->cone_x.tan_angle_sq = tan_sq;
        out->cone_x.sheet_selection = sheet;
        return true;
    } else if (fabs(ax) < eps && fabs(az) < eps && fabs(ay) > eps) {
        *out_type = ALEA_PRIMITIVE_CONE_Y;
        out->cone_y.apex_x = apex_x;
        out->cone_y.apex_y = apex_y;
        out->cone_y.apex_z = apex_z;
        out->cone_y.tan_angle_sq = tan_sq;
        out->cone_y.sheet_selection = sheet;
        return true;
    } else if (fabs(ax) < eps && fabs(ay) < eps && fabs(az) > eps) {
        *out_type = ALEA_PRIMITIVE_CONE_Z;
        out->cone_z.apex_x = apex_x;
        out->cone_z.apex_y = apex_y;
        out->cone_z.apex_z = apex_z;
        out->cone_z.tan_angle_sq = tan_sq;
        out->cone_z.sheet_selection = sheet;
        return true;
    }

    /* General cone -> quadric */
    *out_type = ALEA_PRIMITIVE_QUADRIC;
    double* c = out->quadric.coeffs;
    double t2 = tan_sq;

    c[0] = 1.0 - ax*ax*(1.0 + t2);
    c[1] = 1.0 - ay*ay*(1.0 + t2);
    c[2] = 1.0 - az*az*(1.0 + t2);
    c[3] = -2.0*ax*ay*(1.0 + t2);
    c[4] = -2.0*ay*az*(1.0 + t2);
    c[5] = -2.0*ax*az*(1.0 + t2);
    c[6] = 2.0*((ax*ax*(1.0+t2) - 1.0)*apex_x + ax*ay*(1.0+t2)*apex_y + ax*az*(1.0+t2)*apex_z);
    c[7] = 2.0*(ax*ay*(1.0+t2)*apex_x + (ay*ay*(1.0+t2) - 1.0)*apex_y + ay*az*(1.0+t2)*apex_z);
    c[8] = 2.0*(ax*az*(1.0+t2)*apex_x + ay*az*(1.0+t2)*apex_y + (az*az*(1.0+t2) - 1.0)*apex_z);
    c[9] = (1.0 - ax*ax*(1.0+t2))*apex_x*apex_x
         + (1.0 - ay*ay*(1.0+t2))*apex_y*apex_y
         + (1.0 - az*az*(1.0+t2))*apex_z*apex_z
         - 2.0*ax*ay*(1.0+t2)*apex_x*apex_y
         - 2.0*ay*az*(1.0+t2)*apex_y*apex_z
         - 2.0*ax*az*(1.0+t2)*apex_x*apex_z;
    return true;
}

static bool xform_cone_x(const alea_primitive_data_t* in, const double* mat,
                         alea_primitive_type_t* out_type, alea_primitive_data_t* out) {
    const alea_cone_x_data_t* cone = &in->cone_x;
    if (is_translation_only(mat)) {
        *out_type = ALEA_PRIMITIVE_CONE_X;
        out->cone_x = *cone;
        out->cone_x.apex_x += mat[3];
        out->cone_x.apex_y += mat[7];
        out->cone_x.apex_z += mat[11];
        return true;
    }
    return xform_cone_common(mat, 1, 0, 0, cone->apex_x, cone->apex_y, cone->apex_z,
                             cone->tan_angle_sq, cone->sheet_selection, out_type, out);
}

static bool xform_cone_y(const alea_primitive_data_t* in, const double* mat,
                         alea_primitive_type_t* out_type, alea_primitive_data_t* out) {
    const alea_cone_y_data_t* cone = &in->cone_y;
    if (is_translation_only(mat)) {
        *out_type = ALEA_PRIMITIVE_CONE_Y;
        out->cone_y = *cone;
        out->cone_y.apex_x += mat[3];
        out->cone_y.apex_y += mat[7];
        out->cone_y.apex_z += mat[11];
        return true;
    }
    return xform_cone_common(mat, 0, 1, 0, cone->apex_x, cone->apex_y, cone->apex_z,
                             cone->tan_angle_sq, cone->sheet_selection, out_type, out);
}

static bool xform_cone_z(const alea_primitive_data_t* in, const double* mat,
                         alea_primitive_type_t* out_type, alea_primitive_data_t* out) {
    const alea_cone_z_data_t* cone = &in->cone_z;
    if (is_translation_only(mat)) {
        *out_type = ALEA_PRIMITIVE_CONE_Z;
        out->cone_z = *cone;
        out->cone_z.apex_x += mat[3];
        out->cone_z.apex_y += mat[7];
        out->cone_z.apex_z += mat[11];
        return true;
    }
    return xform_cone_common(mat, 0, 0, 1, cone->apex_x, cone->apex_y, cone->apex_z,
                             cone->tan_angle_sq, cone->sheet_selection, out_type, out);
}

static bool xform_box(const alea_primitive_data_t* in, const double* mat,
                      alea_primitive_type_t* out_type, alea_primitive_data_t* out) {
    if (is_translation_only(mat)) {
        *out_type = ALEA_PRIMITIVE_RPP;
        out->box.min_x = in->box.min_x + mat[3];
        out->box.max_x = in->box.max_x + mat[3];
        out->box.min_y = in->box.min_y + mat[7];
        out->box.max_y = in->box.max_y + mat[7];
        out->box.min_z = in->box.min_z + mat[11];
        out->box.max_z = in->box.max_z + mat[11];
        return true;
    }
    /* Rotated box -> convert to BOX_GENERAL */
    *out_type = ALEA_PRIMITIVE_BOX;
    double cx = in->box.min_x, cy = in->box.min_y, cz = in->box.min_z;
    xform_point(mat, &cx, &cy, &cz);
    out->box_general.corner_x = cx;
    out->box_general.corner_y = cy;
    out->box_general.corner_z = cz;

    double dx = in->box.max_x - in->box.min_x;
    double dy = in->box.max_y - in->box.min_y;
    double dz = in->box.max_z - in->box.min_z;

    out->box_general.v1_x = mat[0]*dx;
    out->box_general.v1_y = mat[4]*dx;
    out->box_general.v1_z = mat[8]*dx;
    out->box_general.v2_x = mat[1]*dy;
    out->box_general.v2_y = mat[5]*dy;
    out->box_general.v2_z = mat[9]*dy;
    out->box_general.v3_x = mat[2]*dz;
    out->box_general.v3_y = mat[6]*dz;
    out->box_general.v3_z = mat[10]*dz;
    return true;
}

static bool xform_quadric(const alea_primitive_data_t* in, const double* mat,
                          alea_primitive_type_t* out_type, alea_primitive_data_t* out) {
    /*
     * Full quadric transform: Q(x) = x^T M x + b^T x + c = 0
     * Where M is symmetric 3x3, b is linear coeffs, c is constant.
     *
     * Under affine transform x' = R*x + t (forward), we have x = R^T*(x' - t)
     * Substituting into Q(x) = 0:
     *   (R^T(x'-t))^T M (R^T(x'-t)) + b^T (R^T(x'-t)) + c = 0
     *   (x'-t)^T R M R^T (x'-t) + b^T R^T (x'-t) + c = 0
     *
     * Let M' = R M R^T, b' = R b, t' = t
     * Expanding: x'^T M' x' - 2 x'^T M' t' + t'^T M' t' + b'^T x' - b'^T t' + c = 0
     *
     * New coefficients:
     *   M_new = M' = R M R^T
     *   b_new = b' - 2 M' t' = R b - 2 R M R^T t
     *   c_new = c + t'^T M' t' - b'^T t' = c + t^T R M R^T t - (R b)^T t
     *
     * Mat layout: [R00,R01,R02,Tx, R10,R11,R12,Ty, R20,R21,R22,Tz]
     */
    *out_type = ALEA_PRIMITIVE_QUADRIC;

    const double* ci = in->quadric.coeffs;
    double* co = out->quadric.coeffs;

    /* Extract R and t from mat */
    double R[3][3] = {
        {mat[0], mat[1], mat[2]},
        {mat[4], mat[5], mat[6]},
        {mat[8], mat[9], mat[10]}
    };
    double t[3] = {mat[3], mat[7], mat[11]};

    /* Build symmetric matrix M from quadric coeffs */
    /* Q = Ax² + By² + Cz² + Dxy + Eyz + Fxz + Gx + Hy + Iz + J */
    /* M = [A, D/2, F/2; D/2, B, E/2; F/2, E/2, C] */
    double M[3][3] = {
        {ci[0],      ci[3]/2.0, ci[5]/2.0},
        {ci[3]/2.0,  ci[1],     ci[4]/2.0},
        {ci[5]/2.0,  ci[4]/2.0, ci[2]}
    };
    double b[3] = {ci[6], ci[7], ci[8]};
    double c = ci[9];

    /* Compute R * M */
    double RM[3][3] = {{0}};
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            for (int k = 0; k < 3; k++) {
                RM[i][j] += R[i][k] * M[k][j];
            }
        }
    }

    /* Compute M' = R * M * R^T */
    double Mp[3][3] = {{0}};
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            for (int k = 0; k < 3; k++) {
                Mp[i][j] += RM[i][k] * R[j][k]; /* R^T[k][j] = R[j][k] */
            }
        }
    }

    /* Compute b' = R * b */
    double bp[3] = {0};
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            bp[i] += R[i][j] * b[j];
        }
    }

    /* Compute M' * t */
    double Mpt[3] = {0};
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            Mpt[i] += Mp[i][j] * t[j];
        }
    }

    /* New linear coeffs: b_new = b' - 2*M'*t */
    double bn[3];
    for (int i = 0; i < 3; i++) {
        bn[i] = bp[i] - 2.0 * Mpt[i];
    }

    /* New constant: c_new = c + t^T * M' * t - b'^T * t */
    double cn = c;
    for (int i = 0; i < 3; i++) {
        cn += t[i] * Mpt[i];  /* t^T * M' * t */
        cn -= bp[i] * t[i];   /* b'^T * t */
    }

    /* Write output coefficients */
    co[0] = Mp[0][0];           /* A */
    co[1] = Mp[1][1];           /* B */
    co[2] = Mp[2][2];           /* C */
    co[3] = 2.0 * Mp[0][1];     /* D */
    co[4] = 2.0 * Mp[1][2];     /* E */
    co[5] = 2.0 * Mp[0][2];     /* F */
    co[6] = bn[0];              /* G */
    co[7] = bn[1];              /* H */
    co[8] = bn[2];              /* I */
    co[9] = cn;                 /* J */

    return true;
}

static bool xform_torus(const alea_primitive_data_t* in, const double* mat,
                        alea_primitive_type_t in_type,
                        alea_primitive_type_t* out_type, alea_primitive_data_t* out) {
    const alea_torus_data_t* t = &in->torus;

    /* Transform center */
    double cx = t->center_x, cy = t->center_y, cz = t->center_z;
    xform_point(mat, &cx, &cy, &cz);

    /* Transform axis */
    double ax = 0, ay = 0, az = 0;
    switch (in_type) {
        case ALEA_PRIMITIVE_TORUS_X: ax = 1; break;
        case ALEA_PRIMITIVE_TORUS_Y: ay = 1; break;
        case ALEA_PRIMITIVE_TORUS_Z: az = 1; break;
        default: return false;
    }
    xform_vector(mat, &ax, &ay, &az);

    double eps = 1e-10;
    alea_primitive_type_t new_type;
    alea_axis_t new_axis;

    if (fabs(ay) < eps && fabs(az) < eps) {
        new_type = ALEA_PRIMITIVE_TORUS_X;
        new_axis = ALEA_AXIS_X;
    } else if (fabs(ax) < eps && fabs(az) < eps) {
        new_type = ALEA_PRIMITIVE_TORUS_Y;
        new_axis = ALEA_AXIS_Y;
    } else if (fabs(ax) < eps && fabs(ay) < eps) {
        new_type = ALEA_PRIMITIVE_TORUS_Z;
        new_axis = ALEA_AXIS_Z;
    } else {
        /* Tori are quartic - can't represent as quadric */
        return false;
    }

    *out_type = new_type;
    out->torus.axis = new_axis;
    out->torus.center_x = cx;
    out->torus.center_y = cy;
    out->torus.center_z = cz;
    out->torus.major_radius = t->major_radius;
    out->torus.minor_radius = t->minor_radius;
    out->torus.axial_semiwidth_B = t->axial_semiwidth_B;
    return true;
}

static bool xform_torus_x(const alea_primitive_data_t* in, const double* mat,
                          alea_primitive_type_t* out_type, alea_primitive_data_t* out) {
    return xform_torus(in, mat, ALEA_PRIMITIVE_TORUS_X, out_type, out);
}

static bool xform_torus_y(const alea_primitive_data_t* in, const double* mat,
                          alea_primitive_type_t* out_type, alea_primitive_data_t* out) {
    return xform_torus(in, mat, ALEA_PRIMITIVE_TORUS_Y, out_type, out);
}

static bool xform_torus_z(const alea_primitive_data_t* in, const double* mat,
                          alea_primitive_type_t* out_type, alea_primitive_data_t* out) {
    return xform_torus(in, mat, ALEA_PRIMITIVE_TORUS_Z, out_type, out);
}

static bool xform_rcc(const alea_primitive_data_t* in, const double* mat,
                      alea_primitive_type_t* out_type, alea_primitive_data_t* out) {
    *out_type = ALEA_PRIMITIVE_RCC;
    out->rcc = in->rcc;
    xform_point(mat, &out->rcc.base_x, &out->rcc.base_y, &out->rcc.base_z);
    xform_vector(mat, &out->rcc.height_x, &out->rcc.height_y, &out->rcc.height_z);
    return true;
}

static bool xform_box_general(const alea_primitive_data_t* in, const double* mat,
                              alea_primitive_type_t* out_type, alea_primitive_data_t* out) {
    *out_type = ALEA_PRIMITIVE_BOX;
    out->box_general = in->box_general;
    xform_point(mat, &out->box_general.corner_x, &out->box_general.corner_y, &out->box_general.corner_z);
    xform_vector(mat, &out->box_general.v1_x, &out->box_general.v1_y, &out->box_general.v1_z);
    xform_vector(mat, &out->box_general.v2_x, &out->box_general.v2_y, &out->box_general.v2_z);
    xform_vector(mat, &out->box_general.v3_x, &out->box_general.v3_y, &out->box_general.v3_z);
    return true;
}

static bool xform_sph(const alea_primitive_data_t* in, const double* mat,
                      alea_primitive_type_t* out_type, alea_primitive_data_t* out) {
    *out_type = ALEA_PRIMITIVE_SPH;
    out->sph = in->sph;
    xform_point(mat, &out->sph.center_x, &out->sph.center_y, &out->sph.center_z);
    return true;
}

static bool xform_trc(const alea_primitive_data_t* in, const double* mat,
                      alea_primitive_type_t* out_type, alea_primitive_data_t* out) {
    *out_type = ALEA_PRIMITIVE_TRC;
    out->trc = in->trc;
    xform_point(mat, &out->trc.base_x, &out->trc.base_y, &out->trc.base_z);
    xform_vector(mat, &out->trc.height_x, &out->trc.height_y, &out->trc.height_z);
    return true;
}

static bool xform_ell(const alea_primitive_data_t* in, const double* mat,
                      alea_primitive_type_t* out_type, alea_primitive_data_t* out) {
    *out_type = ALEA_PRIMITIVE_ELL;
    out->ell = in->ell;
    xform_point(mat, &out->ell.v1_x, &out->ell.v1_y, &out->ell.v1_z);
    xform_point(mat, &out->ell.v2_x, &out->ell.v2_y, &out->ell.v2_z);
    return true;
}

static bool xform_rec(const alea_primitive_data_t* in, const double* mat,
                      alea_primitive_type_t* out_type, alea_primitive_data_t* out) {
    *out_type = ALEA_PRIMITIVE_REC;
    out->rec = in->rec;
    xform_point(mat, &out->rec.base_x, &out->rec.base_y, &out->rec.base_z);
    xform_vector(mat, &out->rec.height_x, &out->rec.height_y, &out->rec.height_z);
    xform_vector(mat, &out->rec.axis1_x, &out->rec.axis1_y, &out->rec.axis1_z);
    xform_vector(mat, &out->rec.axis2_x, &out->rec.axis2_y, &out->rec.axis2_z);
    return true;
}

static bool xform_wed(const alea_primitive_data_t* in, const double* mat,
                      alea_primitive_type_t* out_type, alea_primitive_data_t* out) {
    *out_type = ALEA_PRIMITIVE_WED;
    out->wed = in->wed;
    xform_point(mat, &out->wed.vertex_x, &out->wed.vertex_y, &out->wed.vertex_z);
    xform_vector(mat, &out->wed.v1_x, &out->wed.v1_y, &out->wed.v1_z);
    xform_vector(mat, &out->wed.v2_x, &out->wed.v2_y, &out->wed.v2_z);
    xform_vector(mat, &out->wed.v3_x, &out->wed.v3_y, &out->wed.v3_z);
    return true;
}

static bool xform_rhp(const alea_primitive_data_t* in, const double* mat,
                      alea_primitive_type_t* out_type, alea_primitive_data_t* out) {
    *out_type = ALEA_PRIMITIVE_RHP;
    out->rhp = in->rhp;
    xform_point(mat, &out->rhp.base_x, &out->rhp.base_y, &out->rhp.base_z);
    xform_vector(mat, &out->rhp.height_x, &out->rhp.height_y, &out->rhp.height_z);
    xform_vector(mat, &out->rhp.r1_x, &out->rhp.r1_y, &out->rhp.r1_z);
    xform_vector(mat, &out->rhp.r2_x, &out->rhp.r2_y, &out->rhp.r2_z);
    xform_vector(mat, &out->rhp.r3_x, &out->rhp.r3_y, &out->rhp.r3_z);
    return true;
}

static bool xform_arb(const alea_primitive_data_t* in, const double* mat,
                      alea_primitive_type_t* out_type, alea_primitive_data_t* out) {
    *out_type = ALEA_PRIMITIVE_ARB;
    out->arb = in->arb;
    for (int i = 0; i < in->arb.num_corners; i++) {
        double x = in->arb.corners[i][0];
        double y = in->arb.corners[i][1];
        double z = in->arb.corners[i][2];
        xform_point(mat, &x, &y, &z);
        out->arb.corners[i][0] = x;
        out->arb.corners[i][1] = y;
        out->arb.corners[i][2] = z;
    }
    return true;
}

/* ============================================================================
 * DESCRIPTOR TABLE
 * ============================================================================ */

/* Table indexed by primitive type. Index 0 is unused (types start at 1). */
static const alea_primitive_desc_t g_prim_desc[] = {
    /*                           name          eval             bbox              interval             transform */
    [0]                       = { "invalid",    NULL,            NULL,             NULL,                NULL },
    [ALEA_PRIMITIVE_PLANE]     = { "plane",      eval_plane,      bbox_plane,       interval_plane,      xform_plane },
    [ALEA_PRIMITIVE_SPHERE]    = { "sphere",     eval_sphere,     bbox_sphere,      interval_sphere,     xform_sphere },
    [ALEA_PRIMITIVE_CYLINDER_X]= { "cylinder_x", eval_cylinder_x, bbox_cylinder_x,  interval_cylinder_x, xform_cylinder_x },
    [ALEA_PRIMITIVE_CYLINDER_Y]= { "cylinder_y", eval_cylinder_y, bbox_cylinder_y,  interval_cylinder_y, xform_cylinder_y },
    [ALEA_PRIMITIVE_CYLINDER_Z]= { "cylinder_z", eval_cylinder_z, bbox_cylinder_z,  interval_cylinder_z, xform_cylinder_z },
    [ALEA_PRIMITIVE_CONE_X]    = { "cone_x",     eval_cone_x,     bbox_cone_x,      interval_cone_x,     xform_cone_x },
    [ALEA_PRIMITIVE_CONE_Y]    = { "cone_y",     eval_cone_y,     bbox_cone_y,      interval_cone_y,     xform_cone_y },
    [ALEA_PRIMITIVE_CONE_Z]    = { "cone_z",     eval_cone_z,     bbox_cone_z,      interval_cone_z,     xform_cone_z },
    [ALEA_PRIMITIVE_RPP]       = { "rpp",        eval_box,        bbox_box,         interval_box,        xform_box },
    [ALEA_PRIMITIVE_QUADRIC]   = { "quadric",    eval_quadric,    bbox_quadric,     interval_quadric,    xform_quadric },
    [ALEA_PRIMITIVE_TORUS_X]   = { "torus_x",    eval_torus_x,    bbox_torus_x,     interval_torus_x,    xform_torus_x },
    [ALEA_PRIMITIVE_TORUS_Y]   = { "torus_y",    eval_torus_y,    bbox_torus_y,     interval_torus_y,    xform_torus_y },
    [ALEA_PRIMITIVE_TORUS_Z]   = { "torus_z",    eval_torus_z,    bbox_torus_z,     interval_torus_z,    xform_torus_z },
    [ALEA_PRIMITIVE_RCC]       = { "rcc",        eval_rcc,        bbox_rcc,         interval_rcc,        xform_rcc },
    /* Macrobodies */
    [ALEA_PRIMITIVE_BOX]       = { "box",        eval_box_general, bbox_box_general, interval_box_general, xform_box_general },
    [ALEA_PRIMITIVE_SPH]       = { "sph",        eval_sph,        bbox_sph,         interval_sph,        xform_sph },
    [ALEA_PRIMITIVE_TRC]       = { "trc",        eval_trc,        bbox_trc,         interval_trc,        xform_trc },
    [ALEA_PRIMITIVE_ELL]       = { "ell",        eval_ell,        bbox_ell,         interval_ell,        xform_ell },
    [ALEA_PRIMITIVE_REC]       = { "rec",        eval_rec,        bbox_rec,         interval_rec,        xform_rec },
    [ALEA_PRIMITIVE_WED]       = { "wed",        eval_wed,        bbox_wed,         interval_wed,        xform_wed },
    [ALEA_PRIMITIVE_RHP]       = { "rhp",        eval_rhp,        bbox_rhp,         interval_rhp,        xform_rhp },
    [ALEA_PRIMITIVE_ARB]       = { "arb",        eval_arb,        bbox_arb,         interval_arb,        xform_arb },
};

#define PRIM_DESC_COUNT (sizeof(g_prim_desc) / sizeof(g_prim_desc[0]))

const alea_primitive_desc_t* alea_primitive_get_desc(alea_primitive_type_t type) {
    if ((unsigned)type >= PRIM_DESC_COUNT) return NULL;
    if (!g_prim_desc[type].eval) return NULL;
    return &g_prim_desc[type];
}

const char* alea_primitive_type_name(alea_primitive_type_t type) {
    if ((unsigned)type >= PRIM_DESC_COUNT) return "unknown";
    return g_prim_desc[type].name ? g_prim_desc[type].name : "unknown";
}

bool alea_primitive_transform(
    alea_primitive_type_t in_type,
    const alea_primitive_data_t* in_data,
    const double* mat,
    alea_primitive_type_t* out_type,
    alea_primitive_data_t* out_data)
{
    if (!in_data || !mat || !out_type || !out_data) {
        return false;
    }

    const alea_primitive_desc_t* desc = alea_primitive_get_desc(in_type);
    if (!desc || !desc->transform) {
        /* No transform function - just copy */
        *out_type = in_type;
        *out_data = *in_data;
        return false;
    }

    return desc->transform(in_data, mat, out_type, out_data);
}

/* ============================================================================
 * INTERVAL ARITHMETIC DISPATCH
 * ============================================================================ */

alea_interval_t alea_primitive_interval_eval(
    alea_primitive_type_t type,
    const alea_primitive_data_t* data,
    const alea_bbox_t* box)
{
    const alea_primitive_desc_t* desc = alea_primitive_get_desc(type);
    if (!desc || !desc->interval_eval) {
        /* No interval function - return infinite interval (conservative) */
        return alea_iv_make(-1e30, 1e30);
    }
    return desc->interval_eval(data, box);
}

