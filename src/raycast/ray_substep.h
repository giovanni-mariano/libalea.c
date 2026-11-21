// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
// SPDX-License-Identifier: MPL-2.0
#ifndef ALEA_RAY_SUBSTEP_H
#define ALEA_RAY_SUBSTEP_H
#include "raycast.h"
#include "core/alea_spatial_hier.h"

bool alea_substep_plane_cuts(const alea_ray_t* r, const alea_plane_data_t* p,
    const alea_plane_data_t* q, double t, double* enter, double* exit, double* width);

/* Resolve an overrun interval between two physical crossings. Updates
 * the path only after establishing the owner in its placement. */
bool alea_ray_resolve_short_interval(alea_system_t* sys, const alea_ray_t* world,
    uint32_t first, const alea_matrix_t* first_frame, uint32_t second, const alea_matrix_t* second_frame,
    double lo, double hi, bool adjacent, bool left_lattice, alea_hier_ray_path_t* path);

/* Both adjacent endpoint primitives must be expressed in this cell's frame. */
int alea_ray_adjacent_plane_cell_sign(const alea_system_t* sys, alea_node_id_t node,
    const alea_ray_t* local, uint32_t first, uint32_t second, double lo, double hi);

/* Original CSG sign at a narrow interval midpoint. Zero is inconclusive. */
int alea_ray_interval_midpoint_sign(const alea_system_t* sys, alea_node_id_t root,
    const alea_ray_t* local, uint32_t first, uint32_t second, double t);

/* Check the same interior sample with extra precision, or use the exact
 * midpoint if no FP64 sample exists. A false result is inconclusive. */
bool alea_ray_verify_short_owner(const alea_system_t* sys, const alea_ray_t* ray,
    double lo, double hi, double sample_t, const alea_hier_ray_path_t* path);

/* The caller has matched the parent occurrence of the incident cell to
 * the path of the outgoing cell. surface <= 0 searches all incident cards. */
int alea_substep_fill_transition(alea_system_t* sys, const alea_ray_t* world,
    double t, int incident, int surface, const alea_hier_ray_path_t* path,
    alea_ray_substep_fill_t* out);
int alea_substep_fill_at(alea_system_t* sys, const alea_ray_t* ray, double t,
                         alea_ray_substep_fill_t* out);
#endif
