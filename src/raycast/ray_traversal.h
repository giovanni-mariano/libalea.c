// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
// SPDX-License-Identifier: MPL-2.0
#ifndef ALEA_RAY_TRAVERSAL_H
#define ALEA_RAY_TRAVERSAL_H

#include "raycast.h"
#include "core/alea_system.h"
#include "core/alea_transform.h"

/* Express O + t D in a placement's local frame, preserving the world ray's t. */
void alea_ray_transform_inverse(const alea_matrix_t* frame, const alea_ray_t* world,
                                alea_ray_t* local);

/* Select an incident cell card at exactly t using the scalar walker's card order. */
int alea_ray_cell_surface_at(alea_system_t* sys, const alea_ray_t* local,
                             uint32_t cell, double t, alea_raycast_result_t* result,
                             uint32_t* primitive);

#define MAX_FILL_RAYCAST_DEPTH 32

/* Optional candidate source for root recovery in a scalar segment walk.
 * Return a verified root cell index or -1 to use the hierarchy. */
typedef struct {
    int (*find)(void* context, double x, double y, double z, int crossed_surface);
    void* context;
} alea_ray_root_lookup_t;

int alea_raycast_hier_segments_with_root_lookup(
    alea_system_t* sys, double ox, double oy, double oz,
    double dx, double dy, double dz, double t_max, alea_raycast_result_t* result,
    const alea_ray_root_lookup_t* lookup);

void alea_raycast_tree_primitives(alea_system_t* sys,
                                    const alea_ray_t* ray,
                                    alea_node_id_t node_id,
                                    double t_min, double t_max,
                                    alea_raycast_result_t* result,
                                    const alea_matrix_t* local_to_world);
int alea_lattice_raycast_interval(const alea_ray_t* ray,
                                    const alea_cell_entry_t* cell,
                                    double t_min, double t_max,
                                    double* out_enter, double* out_exit);
int alea_lattice_raycast_step_limit(const alea_ray_t* ray,
                                      const alea_cell_entry_t* cell,
                                      double t_enter, double t_exit);
double alea_lattice_next_boundary(const alea_ray_t* ray,
                                    const alea_cell_entry_t* lat_cell,
                                    double t_min,
                                    double t_max);
int alea_ray_append_selected_segment(alea_raycast_result_t* result, int* previous_cell,
                                     int cell, const alea_ray_segment_t* segment);

void alea_raycast_visible_normal(alea_system_t* sys, const alea_ray_t* ray, double t,
                                 const alea_matrix_t* frame, alea_ray_first_visible_result_t* visible);
#endif
