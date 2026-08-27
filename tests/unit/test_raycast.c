// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/*
 * test_raycast.c - Unit tests for ray casting module
 */

#include "alea_test.h"
#include "alea.h"
#include "alea_geo_validator.h"
#include "alea_mcnp.h"
#include "alea_raycast.h"
#include "alea_slice.h"
#include "raycast/raycast.h"
#include "raycast/ray_intersect.h"
#include "core/alea_system.h"
#include "core/alea_spatial_hier.h"
#include "geo_validator/transition_slice_critical.h"

#define EPS 1e-6

/* ------------------------------------------------------------------------- */
/* Ray-Primitive Intersection Tests                                           */
/* ------------------------------------------------------------------------- */

TEST(ray_sphere_intersection) {
    alea_ray_t ray;
    alea_sphere_data_t sphere = {0, 0, 0, 5.0};
    double t[2];
    int count;

    /* Ray from outside, hitting sphere */
    alea_ray_init(&ray, -10, 0, 0, 1, 0, 0);
    count = ray_intersect_sphere(&ray, &sphere, t);
    ASSERT_EQ(count, 2);
    ASSERT_NEAR(t[0], 5.0, EPS);   /* Enter at x=-5 */
    ASSERT_NEAR(t[1], 15.0, EPS);  /* Exit at x=5 */

    /* Ray from inside sphere */
    alea_ray_init(&ray, 0, 0, 0, 1, 0, 0);
    count = ray_intersect_sphere(&ray, &sphere, t);
    ASSERT_EQ(count, 2);
    ASSERT_NEAR(t[0], -5.0, EPS);  /* Behind origin */
    ASSERT_NEAR(t[1], 5.0, EPS);   /* Forward */

    /* Ray missing sphere */
    alea_ray_init(&ray, -10, 10, 0, 1, 0, 0);
    count = ray_intersect_sphere(&ray, &sphere, t);
    ASSERT_EQ(count, 0);

    /* Ray tangent to sphere (grazing) */
    alea_ray_init(&ray, -10, 5, 0, 1, 0, 0);
    count = ray_intersect_sphere(&ray, &sphere, t);
    ASSERT_EQ(count, 2);
    ASSERT_NEAR(t[0], t[1], EPS);  /* Single touch point */
}

TEST(ray_plane_intersection) {
    alea_ray_t ray;
    alea_plane_data_t plane;
    double t[2], nx, ny, nz;
    int count;

    /* Plane at z=5 (normal pointing +z) */
    plane.a = 0; plane.b = 0; plane.c = 1; plane.d = -5;

    /* Ray pointing towards plane */
    alea_ray_init(&ray, 0, 0, 0, 0, 0, 1);
    count = ray_intersect_plane(&ray, &plane, t, &nx, &ny, &nz);
    ASSERT_EQ(count, 1);
    ASSERT_NEAR(t[0], 5.0, EPS);
    ASSERT_NEAR(nz, 1.0, EPS);

    /* Ray parallel to plane */
    alea_ray_init(&ray, 0, 0, 0, 1, 0, 0);
    count = ray_intersect_plane(&ray, &plane, t, NULL, NULL, NULL);
    ASSERT_EQ(count, 0);

    /* Ray pointing away from plane */
    alea_ray_init(&ray, 0, 0, 0, 0, 0, -1);
    count = ray_intersect_plane(&ray, &plane, t, NULL, NULL, NULL);
    ASSERT_EQ(count, 1);
    ASSERT_NEAR(t[0], -5.0, EPS);  /* Intersection behind origin */
}

TEST(surface_project_along_uses_nearest_signed_native_intersection) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);
    ASSERT(alea_sphere_surface(sys, 101, 0.0, 0.0, 0.0, 1.0) >= 0);

    const double point[3] = {1.00001, 0.0, 0.0};
    const double direction[3] = {2.0, 0.0, 0.0};
    double parameter = 0.0;
    double projected[3] = {0.0, 0.0, 0.0};
    alea_primitive_type_t type = 0;
    ASSERT_EQ(alea_surface_project_along(
                  sys, 101, point, direction, &parameter, projected, &type), 1);
    ASSERT_EQ(type, ALEA_PRIMITIVE_SPHERE);
    ASSERT_NEAR(parameter, -1.0e-5, 1.0e-10);
    ASSERT_NEAR(projected[0], 1.0, 1.0e-10);
    ASSERT_NEAR(projected[1], 0.0, 1.0e-12);
    ASSERT_NEAR(projected[2], 0.0, 1.0e-12);

    ASSERT(alea_torus_z_surface(sys, 102, 0.0, 0.0, 0.0, 2.0, 1.0) >= 0);
    const double torus_point[3] = {3.00001, 0.0, 0.0};
    ASSERT_EQ(alea_surface_project_along(
                  sys, 102, torus_point, direction, &parameter, projected,
                  &type), 1);
    ASSERT_EQ(type, ALEA_PRIMITIVE_TORUS_Z);
    ASSERT_NEAR(parameter, -1.0e-5, 1.0e-9);
    ASSERT_NEAR(projected[0], 3.0, 1.0e-9);

    ASSERT_EQ(alea_surface_project_along(
                  sys, 999, point, direction, &parameter, projected, NULL), -1);
    alea_destroy(sys);
}

TEST(ray_cylinder_z_intersection) {
    alea_ray_t ray;
    alea_cylinder_z_data_t cyl = {0, 0, 3.0};
    double t[2];
    int count;

    /* Ray perpendicular to cylinder axis */
    alea_ray_init(&ray, -10, 0, 5, 1, 0, 0);
    count = ray_intersect_cylinder_z(&ray, &cyl, t);
    ASSERT_EQ(count, 2);
    ASSERT_NEAR(t[0], 7.0, EPS);   /* Enter at x=-3 */
    ASSERT_NEAR(t[1], 13.0, EPS);  /* Exit at x=3 */

    /* Ray parallel to cylinder axis (inside) */
    alea_ray_init(&ray, 1, 1, 0, 0, 0, 1);
    count = ray_intersect_cylinder_z(&ray, &cyl, t);
    ASSERT_EQ(count, 0);  /* No intersection with infinite cylinder wall */

    /* Ray missing cylinder */
    alea_ray_init(&ray, -10, 5, 0, 1, 0, 0);
    count = ray_intersect_cylinder_z(&ray, &cyl, t);
    ASSERT_EQ(count, 0);
}

TEST(ray_box_intersection) {
    alea_ray_t ray;
    alea_box_data_t box = {-1, 1, -2, 2, -3, 3};
    double t[2];
    int count;

    /* Ray through center */
    alea_ray_init(&ray, -10, 0, 0, 1, 0, 0);
    count = ray_intersect_box(&ray, &box, t);
    ASSERT_EQ(count, 2);
    ASSERT_NEAR(t[0], 9.0, EPS);   /* Enter at x=-1 */
    ASSERT_NEAR(t[1], 11.0, EPS);  /* Exit at x=1 */

    /* Ray from inside */
    alea_ray_init(&ray, 0, 0, 0, 1, 0, 0);
    count = ray_intersect_box(&ray, &box, t);
    ASSERT_EQ(count, 2);
    ASSERT_NEAR(t[0], -1.0, EPS);  /* Behind */
    ASSERT_NEAR(t[1], 1.0, EPS);   /* Forward */

    /* Ray missing box */
    alea_ray_init(&ray, -10, 10, 0, 1, 0, 0);
    count = ray_intersect_box(&ray, &box, t);
    ASSERT_EQ(count, 0);

    /* Diagonal ray */
    alea_ray_init(&ray, -10, -10, -10, 1, 1, 1);
    count = ray_intersect_box(&ray, &box, t);
    ASSERT_EQ(count, 2);
    ASSERT(t[0] < t[1]);
}

TEST(ray_cone_z_intersection) {
    alea_ray_t ray;
    alea_cone_z_data_t cone = {0, 0, 0, 1.0, 0};  /* 45-degree cone */
    double t[2];
    int count;

    /* Ray perpendicular to axis */
    alea_ray_init(&ray, -10, 0, 5, 1, 0, 0);  /* At z=5, cone radius is 5 */
    count = ray_intersect_cone_z(&ray, &cone, t);
    ASSERT_EQ(count, 2);
    ASSERT_NEAR(t[0], 5.0, EPS);   /* Enter at x=-5 */
    ASSERT_NEAR(t[1], 15.0, EPS);  /* Exit at x=5 */
}

TEST(ray_rcc_intersection) {
    alea_ray_t ray;
    /* RCC: base at (0,0,0), height vector (0,0,10), radius 2 */
    alea_rcc_data_t rcc = {0, 0, 0, 0, 0, 10, 2.0};
    double t[2];
    int count;

    /* Ray through side */
    alea_ray_init(&ray, -5, 0, 5, 1, 0, 0);
    count = ray_intersect_rcc(&ray, &rcc, t);
    ASSERT_EQ(count, 2);
    ASSERT_NEAR(t[0], 3.0, EPS);   /* Enter at x=-2 */
    ASSERT_NEAR(t[1], 7.0, EPS);   /* Exit at x=2 */

    /* Ray through caps */
    alea_ray_init(&ray, 0, 0, -5, 0, 0, 1);
    count = ray_intersect_rcc(&ray, &rcc, t);
    ASSERT_EQ(count, 2);
    ASSERT_NEAR(t[0], 5.0, EPS);   /* Enter bottom cap */
    ASSERT_NEAR(t[1], 15.0, EPS);  /* Exit top cap */

    /* Ray missing */
    alea_ray_init(&ray, -5, 5, 5, 1, 0, 0);
    count = ray_intersect_rcc(&ray, &rcc, t);
    ASSERT_EQ(count, 0);
}

/* ------------------------------------------------------------------------- */
/* Ray Utilities Tests                                                        */
/* ------------------------------------------------------------------------- */

TEST(ray_normalization) {
    alea_ray_t ray;

    /* Non-unit direction should be normalized */
    alea_ray_init(&ray, 0, 0, 0, 3, 4, 0);
    ASSERT_NEAR(ray.dx, 0.6, EPS);
    ASSERT_NEAR(ray.dy, 0.8, EPS);
    ASSERT_NEAR(ray.dz, 0.0, EPS);

    /* Check magnitude is 1 */
    double mag = sqrt(ray.dx*ray.dx + ray.dy*ray.dy + ray.dz*ray.dz);
    ASSERT_NEAR(mag, 1.0, EPS);
}

TEST(cauchy_crofton_rays_are_deterministic_and_bound_to_sphere) {
    const double center[3] = {3.0, -2.0, 7.0};
    const double radius = 11.0;
    double origins_a[12], directions_a[12];
    double origins_b[12], directions_b[12];
    uint32_t state_a = 12345;
    uint32_t state_b = 12345;

    ASSERT_EQ(alea_generate_cauchy_crofton_rays(
                  center[0], center[1], center[2], radius, &state_a, 4,
                  origins_a, directions_a), 0);
    ASSERT_EQ(alea_generate_cauchy_crofton_rays(
                  center[0], center[1], center[2], radius, &state_b, 4,
                  origins_b, directions_b), 0);
    ASSERT_EQ(state_a, state_b);
    for (size_t i = 0; i < 4; i++) {
        double dx = directions_a[i * 3];
        double dy = directions_a[i * 3 + 1];
        double dz = directions_a[i * 3 + 2];
        double disk_x = origins_a[i * 3] + 2.0 * radius * dx - center[0];
        double disk_y = origins_a[i * 3 + 1] + 2.0 * radius * dy - center[1];
        double disk_z = origins_a[i * 3 + 2] + 2.0 * radius * dz - center[2];
        ASSERT_NEAR(dx * dx + dy * dy + dz * dz, 1.0, EPS);
        ASSERT_NEAR(disk_x * dx + disk_y * dy + disk_z * dz, 0.0, EPS);
        ASSERT(disk_x * disk_x + disk_y * disk_y + disk_z * disk_z <=
               radius * radius + EPS);
        for (size_t axis = 0; axis < 3; axis++) {
            ASSERT_NEAR(origins_a[i * 3 + axis], origins_b[i * 3 + axis], EPS);
            ASSERT_NEAR(directions_a[i * 3 + axis], directions_b[i * 3 + axis], EPS);
        }
    }
    ASSERT_EQ(alea_generate_cauchy_crofton_rays(
                  center[0], center[1], center[2], 0.0, &state_a, 0,
                  NULL, NULL), -1);
}

/* ------------------------------------------------------------------------- */
/* Full Raycast Tests                                                         */
/* ------------------------------------------------------------------------- */

typedef struct {
    size_t count;
    double material_length;
} selected_segment_probe_t;

static int probe_selected_segment(void* context,
                                  const alea_ray_segment_t* segment) {
    selected_segment_probe_t* probe = context;
    probe->count++;
    if (segment->cell_id >= 0 && segment->material_id != 0)
        probe->material_length += segment->t_exit - segment->t_enter;
    return 0;
}

typedef struct {
    size_t counts[2];
    double material_lengths[2];
} batch_selected_segment_probe_t;

static int probe_batch_selected_segment(void* context, size_t ray_index,
                                        const alea_ray_segment_t* segment) {
    batch_selected_segment_probe_t* probe = context;
    if (ray_index >= 2) return -1;
    probe->counts[ray_index]++;
    if (segment->cell_id >= 0 && segment->material_id != 0)
        probe->material_lengths[ray_index] += segment->t_exit - segment->t_enter;
    return 0;
}

TEST(raycast_simple_geometry) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);

    /* Create sphere surface and get interior node */
    int surf_idx = alea_sphere_surface(sys, 1, 0, 0, 0, 5.0);
    alea_node_id_t sphere = alea_surface_at(sys, surf_idx)->neg_node;

    /* Add cell */
    int m1 = alea_add_material(sys, 1);

    alea_add_cell(sys, 1, sphere, m1, -2.7, 0);
    ASSERT_EQ(alea_prepare_query_acceleration(sys), 0);

    /* Cast ray through sphere */
    alea_raycast_result_t result;
    alea_raycast_result_init(&result);

    int rc = alea_raycast(sys, -10, 0, 0, 1, 0, 0, 100, &result);
    ASSERT_EQ(rc, 0);

    /* Should have 2 hits (enter and exit) */
    ASSERT_EQ(result.hits.count, 2);
    ASSERT_NEAR(result.hits.data[0].t, 5.0, EPS);
    ASSERT_NEAR(result.hits.data[1].t, 15.0, EPS);

    /* Should have segments */
    ASSERT(result.segments.count >= 1);

    alea_raycast_result_free(&result);
    alea_destroy(sys);
}

TEST(raycast_selected_segment_visitor_streams_without_publication) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);
    const int surface = alea_sphere_surface(sys, 1, 0, 0, 0, 5.0);
    ASSERT(surface >= 0);
    const alea_node_id_t inside = alea_surface_at(sys, surface)->neg_node;
    const int material = alea_add_material(sys, 1);
    ASSERT(material >= 0);
    ASSERT(alea_add_cell(sys, 1, inside, material, -2.7, 0) >= 0);
    ASSERT_EQ(alea_prepare_query_acceleration(sys), 0);

    alea_ray_t ray;
    ASSERT_EQ(alea_ray_init(&ray, -10, 0, 0, 1, 0, 0), 0);
    alea_raycast_result_t scratch;
    alea_raycast_result_init(&scratch);
    selected_segment_probe_t probe = {0};
    ASSERT_EQ(alea_raycast_hier_visit_segments_nocache(
                  sys, &ray, 100, &scratch,
                  probe_selected_segment, &probe), 0);
    ASSERT_EQ(probe.count, 3);
    ASSERT_NEAR(probe.material_length, 10.0, EPS);
    ASSERT_EQ(scratch.segments.count, 0);
    ASSERT_EQ(scratch.hits.count, 0);

    alea_raycast_result_free(&scratch);
    alea_destroy(sys);
}

TEST(raycast_batch_selected_segment_visitor_preserves_ray_slots) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);
    const int surface = alea_sphere_surface(sys, 1, 0, 0, 0, 5.0);
    ASSERT(surface >= 0);
    const alea_node_id_t inside = alea_surface_at(sys, surface)->neg_node;
    const int material = alea_add_material(sys, 1);
    ASSERT(material >= 0);
    ASSERT(alea_add_cell(sys, 1, inside, material, -2.7, 0) >= 0);

    const double origins[6] = { -10, 0, 0, 20, 0, 0 };
    const double directions[6] = { 1, 0, 0, 1, 0, 0 };
    batch_selected_segment_probe_t probe = {0};
    ASSERT_EQ(alea_raycast_hier_visit_segments_batch_nocache(
                  sys, origins, directions, 2, 100,
                  probe_batch_selected_segment, &probe), 0);
    ASSERT_EQ(probe.counts[0], 3);
    ASSERT_NEAR(probe.material_lengths[0], 10.0, EPS);
    ASSERT_EQ(probe.counts[1], 1);
    ASSERT_NEAR(probe.material_lengths[1], 0.0, EPS);

    alea_destroy(sys);
}

TEST(raycast_path_length) {
    alea_raycast_result_t result;
    alea_raycast_result_init(&result);

    /* Manually create some segments */
    alea_ray_segment_t seg1 = {0, 5, 1, 1, -2.7, -1, -1, -1};
    alea_ray_segment_t seg2 = {5, 10, 2, 2, -8.0, -1, -1, -1};
    alea_ray_segment_t seg3 = {10, 15, 3, 1, -2.7, -1, -1, -1};

    result.segments.data = malloc(3 * sizeof(alea_ray_segment_t));
    ASSERT_NOT_NULL(result.segments.data);
    result.segments.data[0] = seg1;
    result.segments.data[1] = seg2;
    result.segments.data[2] = seg3;
    result.segments.count = 3;
    result.segments.capacity = 3;

    /* Total path length through material 1 */
    double len1 = alea_raycast_path_length(&result, 1);
    ASSERT_NEAR(len1, 10.0, EPS);  /* 5 + 5 */

    /* Total path length through material 2 */
    double len2 = alea_raycast_path_length(&result, 2);
    ASSERT_NEAR(len2, 5.0, EPS);

    /* Total path length through all materials */
    double total = alea_raycast_path_length(&result, -1);
    ASSERT_NEAR(total, 15.0, EPS);

    alea_raycast_result_free(&result);
}

TEST(remove_cells_by_volume_rebuilds_structural_indexes) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);

    int m1 = alea_add_material(sys, 1);
    int m2 = alea_add_material(sys, 2);
    ASSERT(m1 >= 0);
    ASSERT(m2 >= 0);

    int s1 = alea_sphere_surface(sys, 1, 0, 0, 0, 5.0);
    int s2 = alea_sphere_surface(sys, 2, 20, 0, 0, 5.0);
    ASSERT(s1 >= 0);
    ASSERT(s2 >= 0);

    alea_node_id_t c1_root = alea_surface_at(sys, s1)->neg_node;
    alea_node_id_t c2_root = alea_surface_at(sys, s2)->neg_node;
    ASSERT(alea_add_cell(sys, 10, c1_root, m1, -1.0, 0) >= 0);
    ASSERT(alea_add_cell(sys, 20, c2_root, m2, -1.0, 0) >= 0);

    ASSERT_EQ(alea_prepare_query_acceleration(sys), 0);
    ASSERT_NOT_NULL(sys->hier_spatial_index);
    ASSERT_EQ(alea_build_cell_adjacency(sys), 0);
    ASSERT(sys->cells.data[0].surface_index_count > 0);
    ASSERT(sys->cells.data[1].surface_index_count > 0);
    ASSERT(sys->cell_adjacency_built);
    ASSERT_NOT_NULL(sys->surface_cell_offsets);
    ASSERT_NOT_NULL(sys->surface_cell_refs);

    double volumes[2] = {0.0, 1.0};
    int removed = alea_remove_cells_by_volume(sys, volumes, 0.1);
    ASSERT_EQ(removed, 1);
    ASSERT_EQ((int)alea_cell_count(sys), 1);

    alea_cell_info_t info;
    ASSERT_EQ(alea_cell_find_info(sys, 10, &info), -1);
    ASSERT_EQ(alea_cell_find_info(sys, 20, &info), 0);
    ASSERT_EQ(info.cell_id, 20);

    ASSERT(!sys->universe_index_built);
    ASSERT_NULL(sys->hier_spatial_index);
    ASSERT(!sys->cell_adjacency_built);
    ASSERT_NULL(sys->surface_cell_offsets);
    ASSERT_NULL(sys->surface_cell_refs);
    ASSERT_EQ((int)sys->surface_cell_ref_count, 0);
    ASSERT_EQ((int)sys->cells.data[0].surface_index_count, 0);

    ASSERT_EQ(alea_prepare_query_acceleration(sys), 0);
    ASSERT_NOT_NULL(sys->hier_spatial_index);

    int cell_id = -1;
    int material_id = -1;
    ASSERT_EQ(alea_find_cell_at(sys, 20, 0, 0, &cell_id, &material_id), 0);
    ASSERT_EQ(cell_id, 20);
    ASSERT_EQ(material_id, 2);

    alea_destroy(sys);
}

TEST(surface_reference_csr_preserves_exact_cards_and_all_references) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);

    /* Geometrically identical planes deduplicate to one primitive, but their
     * MCNP card identities must remain separate in the reverse index. */
    int exact_a = alea_plane_surface(sys, 101, 1.0, 0.0, 0.0, 0.0);
    int exact_b = alea_plane_surface(sys, 102, 1.0, 0.0, 0.0, 0.0);
    ASSERT(exact_a >= 0);
    ASSERT(exact_b >= 0);
    ASSERT_EQ(alea_add_cell(sys, 11, alea_surface_at(sys, exact_a)->neg_node,
                            ALEA_MATERIAL_VOID, 0.0, 7), 0);
    ASSERT_EQ(alea_add_cell(sys, 22, alea_surface_at(sys, exact_b)->pos_node,
                            ALEA_MATERIAL_VOID, 0.0, 7), 1);
    ASSERT_EQ(alea_build_cell_adjacency(sys), 0);

    alea_surface_reference_stats_t stats;
    ASSERT_EQ(alea_surface_reference_stats(sys, &stats), 0);
    ASSERT(stats.built);
    ASSERT_EQ(stats.surface_count, (size_t)2);
    ASSERT_EQ(stats.reference_count, (size_t)2);
    ASSERT_EQ(stats.max_references_per_surface, (size_t)1);
    ASSERT(stats.memory_bytes >= 2 * sizeof(alea_surface_cell_reference_t));

    alea_surface_cell_reference_t refs[2];
    size_t count = 0;
    ASSERT_EQ(alea_surface_cell_references(sys, 101, refs, 2, &count), 0);
    ASSERT_EQ(count, (size_t)1);
    ASSERT_EQ(refs[0].cell_id, 11);
    ASSERT_EQ(refs[0].universe_id, 7);
    ASSERT_EQ(refs[0].sense, -1);

    ASSERT_EQ(alea_surface_cell_references(sys, 102, refs, 2, &count), 0);
    ASSERT_EQ(count, (size_t)1);
    ASSERT_EQ(refs[0].cell_id, 22);
    ASSERT_EQ(refs[0].sense, 1);
    alea_destroy(sys);

    sys = alea_create();
    ASSERT_NOT_NULL(sys);
    int shared = alea_plane_surface(sys, 201, 1.0, 0.0, 0.0, 0.0);
    ASSERT(shared >= 0);
    for (int i = 0; i < 129; i++) {
        ASSERT_EQ(alea_add_cell(sys, 1000 + i,
                               alea_surface_at(sys, shared)->neg_node,
                               ALEA_MATERIAL_VOID, 0.0, 9), i);
    }
    ASSERT_EQ(alea_build_cell_adjacency(sys), 0);

    /* The pairwise neighbor pool intentionally skips surfaces above 128
     * references. The authoritative CSR must retain every reference. */
    ASSERT_EQ(alea_surface_cell_references(sys, 201, NULL, 0, &count), 0);
    ASSERT_EQ(count, (size_t)129);
    ASSERT_EQ(alea_surface_reference_stats(sys, &stats), 0);
    ASSERT_EQ(stats.reference_count, (size_t)129);
    ASSERT_EQ(stats.max_references_per_surface, (size_t)129);
    alea_destroy(sys);

    sys = alea_create();
    ASSERT_NOT_NULL(sys);
    int sx = alea_plane_surface(sys, 301, 1.0, 0.0, 0.0, 0.0);
    int sy = alea_plane_surface(sys, 302, 0.0, 1.0, 0.0, 0.0);
    ASSERT(sx >= 0);
    ASSERT(sy >= 0);
    alea_node_id_t both_neg = alea_intersection(
        sys, alea_surface_at(sys, sx)->neg_node,
        alea_surface_at(sys, sy)->neg_node);
    alea_node_id_t both_pos = alea_intersection(
        sys, alea_surface_at(sys, sx)->pos_node,
        alea_surface_at(sys, sy)->pos_node);
    ASSERT_EQ(alea_add_cell(sys, 31, both_neg, ALEA_MATERIAL_VOID, 0.0, 4), 0);
    ASSERT_EQ(alea_add_cell(sys, 32, both_pos, ALEA_MATERIAL_VOID, 0.0, 4), 1);
    ASSERT_EQ(alea_build_cell_adjacency(sys), 0);
    ASSERT_EQ(alea_surface_cell_references(sys, 301, refs, 2, &count), 0);
    ASSERT_EQ(count, (size_t)2);
    ASSERT_EQ(alea_surface_cell_references(sys, 302, refs, 2, &count), 0);
    ASSERT_EQ(count, (size_t)2);
    alea_destroy(sys);
}

TEST(point_coverage_classifies_competing_unequal_depth_leaves) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);
    int sphere = alea_sphere_surface(sys, 401, 0.0, 0.0, 0.0, 2.0);
    ASSERT(sphere >= 0);
    alea_node_id_t inside = alea_surface_at(sys, sphere)->neg_node;

    int terminal_root = alea_add_cell(
        sys, 41, inside, ALEA_MATERIAL_VOID, 0.0, 0);
    int fill_root = alea_add_cell(
        sys, 42, inside, ALEA_MATERIAL_VOID, 0.0, 0);
    int terminal_child = alea_add_cell(
        sys, 43, inside, ALEA_MATERIAL_VOID, 0.0, 8);
    ASSERT(terminal_root >= 0);
    ASSERT(fill_root >= 0);
    ASSERT(terminal_child >= 0);
    ASSERT_EQ(alea_set_cell_fill(sys, fill_root, 8, 0), 0);
    ASSERT_EQ(alea_prepare_query_acceleration(sys), 0);

    alea_cell_hit_t hits[8];
    uint64_t keys[8], parents[8];
    uint8_t owners[8];
    int count = alea_find_all_cells_coverage_chain(
        sys, 0.0, 0.0, 0.0, hits, keys, parents, 8);
    ASSERT_EQ(count, 3);

    alea_point_coverage_classification_t classification;
    ASSERT_EQ(alea_classify_point_coverage_chain(
                  hits, keys, parents, (size_t)count, -1,
                  owners, &classification), 0);
    ASSERT_EQ(classification.kind, ALEA_POINT_COVERAGE_OVERLAP);
    ASSERT_EQ(classification.owner_count, (size_t)2);
    ASSERT_EQ(owners[terminal_root], 1);
    ASSERT_EQ(owners[fill_root], 0);
    ASSERT_EQ(owners[terminal_child], 1);

    /* Entering universe 8 directly treats its coordinates as local and
     * starts a fresh occurrence tree at depth zero. */
    count = alea_find_all_cells_in_universe_coverage_chain(
        sys, 8, 0.0, 0.0, 0.0, hits, keys, parents, 8);
    ASSERT_EQ(count, 1);
    ASSERT_EQ(hits[0].cell_id, 43);
    ASSERT_EQ(hits[0].depth, 0);
    ASSERT_EQ(parents[0], UINT64_C(0));
    ASSERT_EQ(alea_classify_point_coverage_chain(
                  hits, keys, parents, (size_t)count, -1,
                  owners, &classification), 0);
    ASSERT_EQ(classification.kind, ALEA_POINT_COVERAGE_UNIQUE);
    ASSERT_EQ(classification.owner_count, (size_t)1);

    alea_destroy(sys);
}

static int check_x_plane_transition(
    alea_system_t* sys, int current_cell_id, int primary_surface_id,
    const int* tied_surface_ids, size_t tied_surface_count,
    alea_transition_result_t* result) {
    const double point[3] = {0.0, 0.0, 0.0};
    const double direction[3] = {1.0, 0.0, 0.0};
    if (alea_prepare_query_acceleration(sys) != 0) return -1;
    return alea_check_transition_local(
        sys, 0, current_cell_id, primary_surface_id,
        tied_surface_ids, tied_surface_count,
        point, direction, NULL, result);
}

TEST(supplied_transition_kernel_classifies_adjacency_and_anomalies) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);
    int surface = alea_plane_surface(sys, 501, 1.0, 0.0, 0.0, 0.0);
    ASSERT(surface >= 0);
    ASSERT_EQ(alea_add_cell(sys, 51, alea_surface_at(sys, surface)->neg_node,
                            ALEA_MATERIAL_VOID, 0.0, 0), 0);
    ASSERT_EQ(alea_add_cell(sys, 52, alea_surface_at(sys, surface)->pos_node,
                            ALEA_MATERIAL_VOID, 0.0, 0), 1);
    alea_transition_result_t result;
    ASSERT_EQ(check_x_plane_transition(sys, 51, 501, NULL, 0, &result), 0);
    ASSERT_EQ(result.kind, ALEA_TRANSITION_VALID);
    ASSERT_EQ(result.after_cell_id, 52);
    ASSERT_EQ(result.primary_candidate_count, (size_t)1);
    ASSERT_EQ(result.primary_containing_count, (size_t)1);
    ASSERT_EQ(result.coverage_fallbacks, (size_t)0);
    ASSERT(result.flags & ALEA_TRANSITION_FLAG_OFFSET_STABLE);
    alea_destroy(sys);

    sys = alea_create();
    ASSERT_NOT_NULL(sys);
    surface = alea_plane_surface(sys, 511, 1.0, 0.0, 0.0, 0.0);
    ASSERT(surface >= 0);
    ASSERT_EQ(alea_add_cell(sys, 61, alea_surface_at(sys, surface)->neg_node,
                            ALEA_MATERIAL_VOID, 0.0, 0), 0);
    ASSERT_EQ(check_x_plane_transition(sys, 61, 511, NULL, 0, &result), 0);
    ASSERT_EQ(result.kind, ALEA_TRANSITION_GAP);
    ASSERT(result.flags & ALEA_TRANSITION_FLAG_PRIMARY_MISSING);
    ASSERT(result.flags & ALEA_TRANSITION_FLAG_COVERAGE_FALLBACK);
    ASSERT(result.coverage_fallbacks > 0);
    alea_destroy(sys);

    sys = alea_create();
    ASSERT_NOT_NULL(sys);
    surface = alea_plane_surface(sys, 521, 1.0, 0.0, 0.0, 0.0);
    ASSERT(surface >= 0);
    ASSERT_EQ(alea_add_cell(sys, 71, alea_surface_at(sys, surface)->neg_node,
                            ALEA_MATERIAL_VOID, 0.0, 0), 0);
    ASSERT_EQ(alea_add_cell(sys, 72, alea_surface_at(sys, surface)->pos_node,
                            ALEA_MATERIAL_VOID, 0.0, 0), 1);
    ASSERT_EQ(alea_add_cell(sys, 73, alea_surface_at(sys, surface)->pos_node,
                            ALEA_MATERIAL_VOID, 0.0, 0), 2);
    ASSERT_EQ(check_x_plane_transition(sys, 71, 521, NULL, 0, &result), 0);
    ASSERT_EQ(result.kind, ALEA_TRANSITION_OVERLAP);
    ASSERT_EQ(result.primary_containing_count, (size_t)2);
    ASSERT_EQ(result.after_owner_count, (size_t)2);
    alea_destroy(sys);

    /* Exact-card mismatch: geometrically coincident cards remain distinct.
     * Without tied evidence the unique after owner is non-adjacent. */
    sys = alea_create();
    ASSERT_NOT_NULL(sys);
    int primary = alea_plane_surface(sys, 531, 1.0, 0.0, 0.0, 0.0);
    int secondary = alea_plane_surface(sys, 532, 1.0, 0.0, 0.0, 0.0);
    ASSERT(primary >= 0);
    ASSERT(secondary >= 0);
    ASSERT_EQ(alea_add_cell(sys, 81, alea_surface_at(sys, primary)->neg_node,
                            ALEA_MATERIAL_VOID, 0.0, 0), 0);
    ASSERT_EQ(alea_add_cell(sys, 82, alea_surface_at(sys, secondary)->pos_node,
                            ALEA_MATERIAL_VOID, 0.0, 0), 1);
    ASSERT_EQ(check_x_plane_transition(sys, 81, 531, NULL, 0, &result), 0);
    ASSERT_EQ(result.kind, ALEA_TRANSITION_NON_ADJACENT);
    ASSERT_EQ(result.after_cell_id, 82);
    alea_destroy(sys);

    /* A cell may use the primary card in both orientations in disjoint
     * Boolean branches.  Concrete before/after containment still proves an
     * exit, so the checker must continue to coverage rather than stop at the
     * card-level sense ambiguity. */
    sys = alea_create();
    ASSERT_NOT_NULL(sys);
    primary = alea_plane_surface(sys, 536, 1.0, 0.0, 0.0, 0.0);
    secondary = alea_plane_surface(sys, 537, 1.0, 0.0, 0.0, 0.0);
    int gate = alea_plane_surface(sys, 538, 0.0, 1.0, 0.0, 0.0);
    ASSERT(primary >= 0);
    ASSERT(secondary >= 0);
    ASSERT(gate >= 0);
    alea_node_id_t negative_branch = alea_intersection(
        sys, alea_surface_at(sys, primary)->neg_node,
        alea_surface_at(sys, gate)->neg_node);
    alea_node_id_t positive_branch = alea_intersection(
        sys, alea_surface_at(sys, primary)->pos_node,
        alea_surface_at(sys, gate)->pos_node);
    alea_node_id_t both_senses = alea_union(
        sys, negative_branch, positive_branch);
    alea_node_id_t after_region = alea_intersection(
        sys, alea_surface_at(sys, secondary)->pos_node,
        alea_surface_at(sys, gate)->neg_node);
    ASSERT_EQ(alea_add_cell(sys, 86, both_senses,
                            ALEA_MATERIAL_VOID, 0.0, 0), 0);
    ASSERT_EQ(alea_add_cell(sys, 87, after_region,
                            ALEA_MATERIAL_VOID, 0.0, 0), 1);
    const double gated_point[3] = {0.0, -1.0, 0.0};
    const double positive_x[3] = {1.0, 0.0, 0.0};
    ASSERT_EQ(alea_prepare_query_acceleration(sys), 0);
    ASSERT_EQ(alea_check_transition_local(
        sys, 0, 86, 536, NULL, 0, gated_point, positive_x, NULL, &result), 0);
    ASSERT_EQ(result.current_sense, 2);
    ASSERT(result.flags & ALEA_TRANSITION_FLAG_CURRENT_BEFORE_CONTAINS);
    ASSERT(!(result.flags & ALEA_TRANSITION_FLAG_CURRENT_AFTER_CONTAINS));
    ASSERT(result.flags & ALEA_TRANSITION_FLAG_COVERAGE_FALLBACK);
    ASSERT_EQ(result.kind, ALEA_TRANSITION_NON_ADJACENT);
    ASSERT_EQ(result.after_cell_id, 87);
    alea_destroy(sys);

    /* If the current cell also carries the coincident card and the event
     * supplies it as tied evidence, retain the surface-chain/corner cause. */
    sys = alea_create();
    ASSERT_NOT_NULL(sys);
    primary = alea_plane_surface(sys, 541, 1.0, 0.0, 0.0, 0.0);
    secondary = alea_plane_surface(sys, 542, 1.0, 0.0, 0.0, 0.0);
    alea_node_id_t current = alea_intersection(
        sys, alea_surface_at(sys, primary)->neg_node,
        alea_surface_at(sys, secondary)->neg_node);
    ASSERT_EQ(alea_add_cell(sys, 91, current,
                            ALEA_MATERIAL_VOID, 0.0, 0), 0);
    ASSERT_EQ(alea_add_cell(sys, 92, alea_surface_at(sys, secondary)->pos_node,
                            ALEA_MATERIAL_VOID, 0.0, 0), 1);
    ASSERT_EQ(check_x_plane_transition(sys, 91, 541, NULL, 0, &result), 0);
    ASSERT_EQ(result.kind, ALEA_TRANSITION_SURFACE_CHAIN_CORNER);
    ASSERT_EQ(result.after_cell_id, 92);
    ASSERT_EQ(result.connecting_surface_id, 542);
    ASSERT(result.flags & ALEA_TRANSITION_FLAG_TIED_SURFACE_CONNECTS);
    alea_destroy(sys);
}

TEST(transition_slice_screen_streams_only_bounded_findings) {
    alea_slice_view_t view;
    alea_slice_view_init(&view, 0.0, 0.0, 0.0,
                         0.0, 0.0, 1.0,
                         0.0, 1.0, 0.0,
                         -2.0, 2.0, -1.0, 1.0);
    alea_transition_slice_options_t options;
    alea_transition_slice_options_init(&options);
    options.horizontal_rays = 4;
    options.vertical_rays = 3;

    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);
    int shared = alea_plane_surface(sys, 551, 1.0, 0.0, 0.0, 0.0);
    ASSERT(shared >= 0);
    ASSERT(alea_add_cell(sys, 101, alea_surface_at(sys, shared)->neg_node,
                         ALEA_MATERIAL_VOID, 0.0, 0) >= 0);
    ASSERT(alea_add_cell(sys, 102, alea_surface_at(sys, shared)->pos_node,
                         ALEA_MATERIAL_VOID, 0.0, 0) >= 0);
    alea_transition_slice_result_t* result =
        alea_transition_slice_result_create();
    ASSERT_NOT_NULL(result);
    ASSERT_EQ(alea_transition_slice_screen(sys, &view, &options, result), 0);
    alea_transition_slice_stats_t stats;
    ASSERT_EQ(alea_transition_slice_stats(result, &stats), 0);
    ASSERT(stats.complete);
    ASSERT_EQ(stats.requested_rays, (size_t)7);
    ASSERT_EQ(stats.executed_rays, (size_t)7);
    ASSERT_EQ(stats.valid_transitions, (size_t)4);
    ASSERT_EQ(stats.coverage_fallbacks, (size_t)0);
    ASSERT_EQ(alea_transition_slice_finding_count(result), (size_t)0);
    ASSERT(stats.peak_live_event_bytes <= options.max_scratch_bytes);
    alea_transition_slice_result_destroy(result);
    alea_destroy(sys);

    /* Geometrically coincident but card-distinct ownership is retained as a
     * finding on each horizontal row; vertical rows have no crossing. */
    sys = alea_create();
    ASSERT_NOT_NULL(sys);
    int primary = alea_plane_surface(sys, 561, 1.0, 0.0, 0.0, 0.0);
    int other = alea_plane_surface(sys, 562, 1.0, 0.0, 0.0, 0.0);
    ASSERT(primary >= 0 && other >= 0);
    ASSERT(alea_add_cell(sys, 111, alea_surface_at(sys, primary)->neg_node,
                         ALEA_MATERIAL_VOID, 0.0, 0) >= 0);
    ASSERT(alea_add_cell(sys, 112, alea_surface_at(sys, other)->pos_node,
                         ALEA_MATERIAL_VOID, 0.0, 0) >= 0);
    result = alea_transition_slice_result_create();
    ASSERT_NOT_NULL(result);
    ASSERT_EQ(alea_transition_slice_screen(sys, &view, &options, result), 0);
    ASSERT_EQ(alea_transition_slice_stats(result, &stats), 0);
    ASSERT(stats.complete);
    ASSERT_EQ(alea_transition_slice_finding_count(result), (size_t)4);
    ASSERT_EQ(alea_transition_slice_component_count(result), (size_t)1);
    alea_transition_slice_component_t component;
    ASSERT_EQ(alea_transition_slice_component_get(result, 0, &component), 0);
    ASSERT_EQ(component.kind, ALEA_TRANSITION_NON_ADJACENT);
    ASSERT_EQ(component.orientation, ALEA_TRANSITION_SLICE_HORIZONTAL);
    ASSERT_EQ(component.finding_count, (size_t)4);
    ASSERT_NEAR(component.uv_min[0], 0.0, EPS);
    ASSERT_NEAR(component.uv_max[0], 0.0, EPS);
    ASSERT_NEAR(component.uv_min[1], -0.75, EPS);
    ASSERT_NEAR(component.uv_max[1], 0.75, EPS);
    ASSERT_EQ(stats.valid_transitions, (size_t)0);
    ASSERT(stats.coverage_fallbacks > 0);
    for (size_t i = 0; i < alea_transition_slice_finding_count(result); i++) {
        alea_transition_slice_finding_t finding;
        ASSERT_EQ(alea_transition_slice_finding_get(result, i, &finding), 0);
        ASSERT_EQ(finding.transition.kind, ALEA_TRANSITION_NON_ADJACENT);
        ASSERT_EQ(finding.orientation, ALEA_TRANSITION_SLICE_HORIZONTAL);
        ASSERT_EQ(finding.transition.current_cell_id, 111);
        ASSERT_EQ(finding.transition.after_cell_id, 112);
        ASSERT_EQ(finding.base_ray_index, i);
        ASSERT_EQ(finding.refinement_depth, (uint32_t)0);
        ASSERT_NEAR(finding.transverse_coordinate,
                    -0.75 + 0.5 * (double)i, EPS);
        ASSERT_NEAR(finding.uv[0], 0.0, EPS);
        ASSERT_NEAR(finding.world_point[0], 0.0, EPS);
    }

    options.vertical_rays = 0;
    options.max_refinement_depth = 1;
    ASSERT_EQ(alea_transition_slice_screen(sys, &view, &options, result), 0);
    ASSERT_EQ(alea_transition_slice_stats(result, &stats), 0);
    ASSERT(stats.complete);
    ASSERT(stats.converged);
    ASSERT_EQ(stats.refinement_status,
              ALEA_TRANSITION_SLICE_REFINEMENT_CONVERGED);
    ASSERT_EQ(stats.executed_rays, (size_t)4);
    ASSERT_EQ(stats.refined_rays_executed, (size_t)0);
    options.vertical_rays = 3;
    options.max_refinement_depth = 0;

    options.max_findings = 1;
    ASSERT_EQ(alea_transition_slice_screen(sys, &view, &options, result), 0);
    ASSERT_EQ(alea_transition_slice_stats(result, &stats), 0);
    ASSERT(!stats.complete);
    ASSERT_EQ(stats.stop_reason, ALEA_TRANSITION_SLICE_STOP_MAX_FINDINGS);
    ASSERT_EQ(alea_transition_slice_finding_count(result), (size_t)1);

    options.max_findings = 10;
    options.max_coverage_fallbacks = 1;
    ASSERT_EQ(alea_transition_slice_screen(sys, &view, &options, result), 0);
    ASSERT_EQ(alea_transition_slice_stats(result, &stats), 0);
    ASSERT(!stats.complete);
    ASSERT_EQ(stats.stop_reason,
              ALEA_TRANSITION_SLICE_STOP_MAX_COVERAGE_FALLBACKS);
    ASSERT_EQ(stats.coverage_fallbacks, (size_t)1);
    ASSERT_EQ(alea_transition_slice_finding_count(result), (size_t)1);
    alea_transition_slice_finding_t truncated;
    ASSERT_EQ(alea_transition_slice_finding_get(result, 0, &truncated), 0);
    ASSERT_EQ(truncated.transition.kind, ALEA_TRANSITION_TRUNCATED);
    alea_transition_slice_result_destroy(result);
    alea_destroy(sys);
}

TEST(transition_slice_screen_batch_respects_scratch_budget) {
    enum { PAGE_COUNT = 4 };
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);
    int primary = alea_plane_surface(sys, 563, 1.0, 0.0, 0.0, 0.0);
    int other = alea_plane_surface(sys, 564, 1.0, 0.0, 0.0, 0.0);
    ASSERT(primary >= 0 && other >= 0);
    ASSERT(alea_add_cell(sys, 113, alea_halfspace(sys, primary, -1),
                         ALEA_MATERIAL_VOID, 0.0, 0) >= 0);
    ASSERT(alea_add_cell(sys, 114, alea_halfspace(sys, other, 1),
                         ALEA_MATERIAL_VOID, 0.0, 0) >= 0);

    alea_slice_view_t views[PAGE_COUNT];
    alea_transition_slice_result_t* results[PAGE_COUNT];
    for (size_t page = 0; page < PAGE_COUNT; page++) {
        double centre = -0.75 + 0.5 * (double)page;
        alea_slice_view_init(&views[page], 0.0, 0.0, 0.0,
                             0.0, 0.0, 1.0,
                             0.0, 1.0, 0.0,
                             -0.25, 0.25, centre - 0.1, centre + 0.1);
        results[page] = alea_transition_slice_result_create();
        ASSERT_NOT_NULL(results[page]);
    }
    alea_transition_slice_options_t options;
    alea_transition_slice_options_init(&options);
    options.horizontal_rays = 1;
    options.vertical_rays = 0;
    options.max_scratch_bytes = 1024;
    options.max_critical_scratch_bytes = 1024;

    alea_transition_slice_batch_stats_t batch;
    ASSERT_EQ(alea_transition_slice_screen_batch(
        sys, views, PAGE_COUNT, &options, 4, 4096, results, &batch), 0);
    ASSERT_EQ(batch.page_count, (size_t)PAGE_COUNT);
    ASSERT_EQ(batch.completed_page_count, (size_t)PAGE_COUNT);
    ASSERT_EQ(batch.requested_workers, (size_t)4);
    ASSERT(batch.actual_workers >= (size_t)1);
    ASSERT(batch.actual_workers <= (size_t)2);
    ASSERT_EQ(batch.reserved_scratch_bytes_per_worker, (uint64_t)2048);
    ASSERT(batch.reserved_parallel_scratch_bytes >=
           (uint64_t)batch.actual_workers * 2048);
    ASSERT(batch.reserved_parallel_scratch_bytes <= (uint64_t)4096);
    for (size_t page = 0; page < PAGE_COUNT; page++) {
        alea_transition_slice_stats_t stats;
        ASSERT_EQ(alea_transition_slice_stats(results[page], &stats), 0);
        ASSERT(stats.complete);
        ASSERT_EQ(stats.executed_rays, (size_t)1);
        ASSERT_EQ(alea_transition_slice_finding_count(results[page]),
                  (size_t)1);
    }

    ASSERT_EQ(alea_transition_slice_screen_batch(
        sys, views, PAGE_COUNT, &options, 4, 0, results, &batch), 0);
    ASSERT_EQ(batch.actual_workers, (size_t)1);
    ASSERT_EQ(batch.completed_page_count, (size_t)PAGE_COUNT);
    for (size_t page = 0; page < PAGE_COUNT; page++)
        alea_transition_slice_result_destroy(results[page]);
    alea_destroy(sys);
}

TEST(transition_slice_screen_refines_changed_event_signatures) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);
    int split = alea_plane_surface(sys, 571, 0.0, 1.0, 0.0, -0.25);
    int shared = alea_plane_surface(sys, 572, 1.0, 0.0, 0.0, 0.0);
    int top_left = alea_plane_surface(sys, 573, 1.0, 0.0, 0.0, 0.0);
    int top_right = alea_plane_surface(sys, 574, 1.0, 0.0, 0.0, 0.0);
    ASSERT(split >= 0 && shared >= 0 && top_left >= 0 && top_right >= 0);
    alea_node_id_t bottom_left = alea_intersection(
        sys, alea_surface_at(sys, split)->neg_node,
        alea_surface_at(sys, shared)->neg_node);
    alea_node_id_t bottom_right = alea_intersection(
        sys, alea_surface_at(sys, split)->neg_node,
        alea_surface_at(sys, shared)->pos_node);
    alea_node_id_t upper_left = alea_intersection(
        sys, alea_surface_at(sys, split)->pos_node,
        alea_surface_at(sys, top_left)->neg_node);
    alea_node_id_t upper_right = alea_intersection(
        sys, alea_surface_at(sys, split)->pos_node,
        alea_surface_at(sys, top_right)->pos_node);
    ASSERT(alea_add_cell(sys, 121, bottom_left,
                         ALEA_MATERIAL_VOID, 0.0, 0) >= 0);
    ASSERT(alea_add_cell(sys, 122, bottom_right,
                         ALEA_MATERIAL_VOID, 0.0, 0) >= 0);
    ASSERT(alea_add_cell(sys, 123, upper_left,
                         ALEA_MATERIAL_VOID, 0.0, 0) >= 0);
    ASSERT(alea_add_cell(sys, 124, upper_right,
                         ALEA_MATERIAL_VOID, 0.0, 0) >= 0);

    alea_slice_view_t view;
    alea_slice_view_init(&view, 0.0, 0.0, 0.0,
                         0.0, 0.0, 1.0, 0.0, 1.0, 0.0,
                         -2.0, 2.0, -1.0, 1.0);
    alea_transition_slice_options_t options;
    alea_transition_slice_options_init(&options);
    options.horizontal_rays = 2;
    options.vertical_rays = 0;
    options.max_refinement_depth = 1;
    options.max_rays = 8;
    alea_transition_slice_result_t* result =
        alea_transition_slice_result_create();
    ASSERT_NOT_NULL(result);
    ASSERT_EQ(alea_transition_slice_screen(sys, &view, &options, result), 0);
    alea_transition_slice_stats_t stats;
    ASSERT_EQ(alea_transition_slice_stats(result, &stats), 0);
    ASSERT(stats.complete);
    ASSERT(!stats.converged);
    ASSERT_EQ(stats.refinement_status,
              ALEA_TRANSITION_SLICE_REFINEMENT_MAX_DEPTH);
    ASSERT_EQ(stats.requested_rays, (size_t)2);
    ASSERT_EQ(stats.executed_rays, (size_t)3);
    ASSERT_EQ(stats.refined_rays_executed, (size_t)1);
    ASSERT_EQ(stats.max_refinement_depth_reached, (uint32_t)1);
    ASSERT(stats.peak_row_scratch_bytes <= options.max_row_scratch_bytes);
    options.enable_critical_refinement = 1;
    ASSERT_EQ(alea_transition_slice_screen(sys, &view, &options, result), 0);
    ASSERT_EQ(alea_transition_slice_stats(result, &stats), 0);
    ASSERT(stats.complete);
    ASSERT(stats.critical_enabled);
    ASSERT(stats.critical_complete);
    ASSERT_EQ(stats.critical_stop_reason,
              ALEA_TRANSITION_SLICE_CRITICAL_NONE);
    ASSERT(stats.critical_tiles_processed > 0);
    ASSERT(stats.critical_region_hits > 0);
    ASSERT(stats.critical_occurrence_seed_points >= (size_t)9);
    ASSERT(stats.critical_occurrence_paths > 0);
    ASSERT(stats.critical_occurrence_universe_queries > 0);
    ASSERT_EQ(stats.critical_root_region_fallbacks, (size_t)0);
    /* Surface cards 572, 573, and 574 share primitive geometry.  Critical
     * enumeration must retain their exact card identities while deduplicating
     * repeated references to card 572 in the same root occurrence. */
    ASSERT(stats.critical_curves >= (size_t)4);
    ASSERT(stats.critical_duplicate_surface_occurrences > 0);
    ASSERT(stats.peak_critical_curves <= options.max_curves_per_tile);
    ASSERT(stats.critical_points > 0);
    ASSERT(stats.critical_probes > 0);
    ASSERT_EQ(alea_transition_slice_critical_finding_count(result),
              stats.critical_findings);
    ASSERT_EQ(stats.critical_findings, stats.critical_probe_findings);
    if (stats.critical_findings > 0) {
        alea_transition_slice_critical_finding_t critical_finding;
        ASSERT_EQ(alea_transition_slice_critical_finding_get(
                      result, 0, &critical_finding), 0);
        ASSERT(critical_finding.source_cell_id > 0);
        ASSERT(critical_finding.source_surface_id > 0);
        ASSERT(critical_finding.radius > 0.0);
        ASSERT(critical_finding.transition.kind != ALEA_TRANSITION_VALID);
        ASSERT(critical_finding.boundary_piece_count > 0);
        const alea_transition_slice_boundary_piece_t* piece =
            &critical_finding.boundary_pieces[0];
        ASSERT_EQ(piece->surface_id, critical_finding.source_surface_id);
        ASSERT(piece->role_flags &
               ALEA_TRANSITION_SLICE_BOUNDARY_ROLE_SOURCE);
        ASSERT_EQ(piece->point_count,
                  (size_t)ALEA_TRANSITION_SLICE_BOUNDARY_POINT_CAPACITY);
        for (size_t point = 0; point < piece->point_count; point++) {
            ASSERT(isfinite(piece->uv[point][0]));
            ASSERT(isfinite(piece->uv[point][1]));
        }
    }
    ASSERT(stats.critical_boundary_evidence > 0);
    ASSERT_EQ(stats.critical_boundary_evidence, stats.critical_findings);
    ASSERT_EQ(stats.omitted_critical_boundary_evidence, (size_t)0);
    ASSERT(stats.critical_curve_pair_candidates > 0);
    ASSERT(stats.critical_curve_pairs_tested > 0);
    ASSERT(stats.critical_active_segments > 0);
    ASSERT(stats.critical_sector_witnesses > 0);
    ASSERT(stats.peak_critical_scratch_bytes <=
           options.max_critical_scratch_bytes);
    ASSERT(alea_transition_slice_refinement_frontier_count(result) > 0);
    ASSERT(alea_transition_slice_critical_tile_count(result) > 0);
    alea_transition_slice_refinement_frontier_t frontier;
    ASSERT_EQ(alea_transition_slice_refinement_frontier_get(
                  result, 0, &frontier), 0);
    ASSERT_EQ(frontier.orientation, ALEA_TRANSITION_SLICE_HORIZONTAL);
    ASSERT(frontier.uv_max[0] > frontier.uv_min[0]);
    alea_transition_slice_critical_tile_t tile;
    ASSERT_EQ(alea_transition_slice_critical_tile_get(result, 0, &tile), 0);
    ASSERT(tile.source_flags &
           (1u << ALEA_TRANSITION_SLICE_TILE_SOURCE_REFINEMENT_FRONTIER));
    options.horizontal_rays = 0;
    options.max_refinement_depth = 0;
    options.critical_full_view = 1;
    ASSERT_EQ(alea_transition_slice_screen(sys, &view, &options, result), 0);
    ASSERT_EQ(alea_transition_slice_stats(result, &stats), 0);
    ASSERT(stats.complete);
    ASSERT(stats.critical_complete);
    ASSERT_EQ(stats.executed_rays, (size_t)0);
    ASSERT_EQ(alea_transition_slice_critical_tile_count(result), (size_t)1);
    ASSERT_EQ(alea_transition_slice_critical_tile_get(result, 0, &tile), 0);
    ASSERT(tile.source_flags &
           (1u << ALEA_TRANSITION_SLICE_TILE_SOURCE_FULL_VIEW));
    ASSERT(stats.critical_findings > 0);
    options.horizontal_rays = 2;
    options.max_refinement_depth = 1;
    options.critical_full_view = 0;
    options.max_critical_boundary_evidence = 1;
    ASSERT_EQ(alea_transition_slice_screen(sys, &view, &options, result), 0);
    ASSERT_EQ(alea_transition_slice_stats(result, &stats), 0);
    ASSERT(stats.critical_complete);
    ASSERT(stats.critical_boundary_evidence <= (size_t)1);
    if (stats.critical_findings > stats.critical_boundary_evidence)
        ASSERT(stats.omitted_critical_boundary_evidence > 0);
    options.max_critical_boundary_evidence = 1024;
    options.max_critical_points = 1;
    ASSERT_EQ(alea_transition_slice_screen(sys, &view, &options, result), 0);
    ASSERT_EQ(alea_transition_slice_stats(result, &stats), 0);
    ASSERT(!stats.critical_complete);
    ASSERT_EQ(stats.critical_stop_reason,
              ALEA_TRANSITION_SLICE_CRITICAL_MAX_POINTS);
    options.max_critical_points = 2048;
    options.max_critical_probes = 1;
    ASSERT_EQ(alea_transition_slice_screen(sys, &view, &options, result), 0);
    ASSERT_EQ(alea_transition_slice_stats(result, &stats), 0);
    ASSERT(!stats.critical_complete);
    ASSERT_EQ(stats.critical_stop_reason,
              ALEA_TRANSITION_SLICE_CRITICAL_MAX_PROBES);
    ASSERT_EQ(stats.critical_probes, (size_t)1);
    options.max_critical_probes = 4096;
    options.max_curve_pairs = 1;
    ASSERT_EQ(alea_transition_slice_screen(sys, &view, &options, result), 0);
    ASSERT_EQ(alea_transition_slice_stats(result, &stats), 0);
    ASSERT(!stats.critical_complete);
    ASSERT_EQ(stats.critical_stop_reason,
              ALEA_TRANSITION_SLICE_CRITICAL_MAX_CURVE_PAIRS);
    ASSERT_EQ(stats.critical_curve_pairs_tested, (size_t)1);
    options.max_curve_pairs = 100000;
    options.max_critical_sector_witnesses = 1;
    ASSERT_EQ(alea_transition_slice_screen(sys, &view, &options, result), 0);
    ASSERT_EQ(alea_transition_slice_stats(result, &stats), 0);
    ASSERT(!stats.critical_complete);
    ASSERT_EQ(stats.critical_stop_reason,
              ALEA_TRANSITION_SLICE_CRITICAL_MAX_SECTOR_WITNESSES);
    ASSERT_EQ(stats.critical_sector_witnesses, (size_t)1);
    options.max_critical_sector_witnesses = 8192;
    options.max_curves_per_tile = 1;
    ASSERT_EQ(alea_transition_slice_screen(sys, &view, &options, result), 0);
    ASSERT_EQ(alea_transition_slice_stats(result, &stats), 0);
    ASSERT(!stats.complete);
    ASSERT(!stats.critical_complete);
    ASSERT_EQ(stats.critical_stop_reason,
              ALEA_TRANSITION_SLICE_CRITICAL_MAX_CURVES);
    ASSERT(stats.peak_critical_curves <= (size_t)1);
    options.max_curves_per_tile = 1024;
    options.enable_critical_refinement = 0;
    options.max_rays = 2;
    ASSERT_EQ(alea_transition_slice_screen(sys, &view, &options, result), 0);
    ASSERT_EQ(alea_transition_slice_stats(result, &stats), 0);
    ASSERT(!stats.complete);
    ASSERT_EQ(stats.stop_reason, ALEA_TRANSITION_SLICE_STOP_MAX_RAYS);
    ASSERT_EQ(stats.executed_rays, (size_t)2);
    ASSERT_EQ(stats.refinement_status,
              ALEA_TRANSITION_SLICE_REFINEMENT_STOPPED);

    options.max_rays = 8;
    options.min_transverse_spacing = 0.6;
    ASSERT_EQ(alea_transition_slice_screen(sys, &view, &options, result), 0);
    ASSERT_EQ(alea_transition_slice_stats(result, &stats), 0);
    ASSERT(stats.complete);
    ASSERT(!stats.converged);
    ASSERT_EQ(stats.refinement_status,
              ALEA_TRANSITION_SLICE_REFINEMENT_MIN_SPACING);
    ASSERT_EQ(stats.executed_rays, (size_t)2);
    alea_transition_slice_result_destroy(result);
    alea_destroy(sys);
}

TEST(transition_slice_critical_enumeration_preserves_transformed_occurrences) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);
    int parent_material = alea_add_material(sys, 1);
    int child_material = alea_add_material(sys, 2);
    int left_boundary = alea_sphere_surface(sys, 701, -5.0, 0.0, 0.0, 3.0);
    int right_boundary = alea_sphere_surface(sys, 702, 5.0, 0.0, 0.0, 3.0);
    int child_left = alea_plane_surface(sys, 703, 1.0, 0.0, 0.0, 0.0);
    int child_right = alea_plane_surface(sys, 704, 1.0, 0.0, 0.0, 0.0);
    ASSERT(parent_material >= 0 && child_material >= 0);
    ASSERT(left_boundary >= 0 && right_boundary >= 0);
    ASSERT(child_left >= 0 && child_right >= 0);
    int left = alea_add_cell(
        sys, 201, alea_surface_at(sys, left_boundary)->neg_node,
        parent_material, -1.0, 0);
    int right = alea_add_cell(
        sys, 202, alea_surface_at(sys, right_boundary)->neg_node,
        parent_material, -1.0, 0);
    ASSERT(left >= 0 && right >= 0);
    ASSERT(alea_add_cell(
        sys, 203, alea_surface_at(sys, child_left)->neg_node,
        child_material, -2.0, 10) >= 0);
    ASSERT(alea_add_cell(
        sys, 204, alea_surface_at(sys, child_right)->pos_node,
        child_material, -2.0, 10) >= 0);
    const double left_translation[3] = {-5.0, 0.0, 0.0};
    const double right_translation[3] = {5.0, 0.0, 0.0};
    ASSERT_EQ(alea_add_transform(sys, 201, left_translation, 3, 0), 0);
    ASSERT_EQ(alea_add_transform(sys, 202, right_translation, 3, 0), 0);
    ASSERT_EQ(alea_set_fill(sys, left, 10, 201), 0);
    ASSERT_EQ(alea_set_fill(sys, right, 10, 202), 0);

    alea_slice_view_t view;
    alea_slice_view_init(&view, 0.0, 0.0, 0.0,
                         0.0, 0.0, 1.0, 0.0, 1.0, 0.0,
                         -9.0, 9.0, -1.0, 1.0);
    alea_transition_slice_options_t options;
    alea_transition_slice_options_init(&options);
    options.horizontal_rays = 3;
    options.vertical_rays = 0;
    options.enable_critical_refinement = 1;
    options.critical_tile_padding = 0.1;
    alea_transition_slice_result_t* result =
        alea_transition_slice_result_create();
    ASSERT_NOT_NULL(result);
    ASSERT_EQ(alea_transition_slice_screen(sys, &view, &options, result), 0);
    alea_transition_slice_stats_t stats;
    ASSERT_EQ(alea_transition_slice_stats(result, &stats), 0);
    ASSERT(stats.critical_complete);
    ASSERT_EQ(stats.critical_root_region_fallbacks, (size_t)0);
    ASSERT(stats.critical_occurrence_paths > 0);
    /* Both distinct placements of universe 10 must be queried.  The root is
     * represented by seed-path neighborhoods rather than a broad BLAS query,
     * and cards 703/704 must survive independently in both occurrences. */
    ASSERT(stats.critical_occurrence_universe_queries >= (size_t)2);
    ASSERT(stats.critical_curves >= (size_t)4);
    ASSERT(alea_transition_slice_component_count(result) >= (size_t)2);
    alea_transition_slice_result_destroy(result);
    alea_destroy(sys);
}

TEST(transition_slice_critical_enumeration_preserves_lattice_occurrences) {
    mcnp_model_t* model = mcnp_load("tests/data/mcnp_lattice_eval.mcnp");
    if (!model) SKIP("Test data file not found");
    alea_system_t* sys = model->sys;
    int duplicate = alea_cylinder_z_surface(sys, 705, 0.0, 0.0, 0.3);
    ASSERT(duplicate >= 0);
    int duplicate_cell = alea_add_cell(
        sys, 205, alea_surface_at(sys, duplicate)->neg_node,
        ALEA_MATERIAL_VOID, 0.0, 1);
    ASSERT(duplicate_cell >= 0);

    alea_slice_view_t view;
    alea_slice_view_init(&view, 0.0, 0.0, 0.0,
                         0.0, 0.0, 1.0, 0.0, 1.0, 0.0,
                         -1.6, 1.6, -1.6, 1.6);
    alea_transition_slice_options_t options;
    alea_transition_slice_options_init(&options);
    options.horizontal_rays = 8;
    options.vertical_rays = 8;
    options.coverage_uniform_probes_per_ray = 16;
    options.coverage_probe_selected_intervals = 1;
    options.enable_critical_refinement = 1;
    options.critical_tile_padding = 0.05;
    alea_transition_slice_result_t* result =
        alea_transition_slice_result_create();
    ASSERT_NOT_NULL(result);
    ASSERT_EQ(alea_transition_slice_screen(sys, &view, &options, result), 0);
    alea_transition_slice_stats_t stats;
    ASSERT_EQ(alea_transition_slice_stats(result, &stats), 0);
    ASSERT(alea_transition_slice_coverage_component_count(result) > 0);
    ASSERT(stats.critical_tiles > 0);
    ASSERT(stats.critical_occurrence_paths > 0);
    ASSERT(stats.critical_occurrence_universe_queries > (size_t)2);
    ASSERT_EQ(stats.critical_root_region_fallbacks, (size_t)0);
    /* Surface 10 and its duplicate card 705 occur in more than one concrete
     * lattice element and must not collapse to one universe-definition key. */
    ASSERT(stats.critical_curves >= (size_t)4);
    alea_transition_slice_result_destroy(result);
    mcnp_model_destroy(model);
}

TEST(transition_slice_screen_bounds_component_publication) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);
    int split = alea_plane_surface(sys, 581, 0.0, 1.0, 0.0, 0.0);
    int bottom_left = alea_plane_surface(sys, 582, 1.0, 0.0, 0.0, 0.0);
    int bottom_right = alea_plane_surface(sys, 583, 1.0, 0.0, 0.0, 0.0);
    int top_left = alea_plane_surface(sys, 584, 1.0, 0.0, 0.0, 0.0);
    int top_right = alea_plane_surface(sys, 585, 1.0, 0.0, 0.0, 0.0);
    ASSERT(split >= 0 && bottom_left >= 0 && bottom_right >= 0 &&
           top_left >= 0 && top_right >= 0);
    alea_node_id_t bottom_left_region = alea_intersection(
        sys, alea_surface_at(sys, split)->neg_node,
        alea_surface_at(sys, bottom_left)->neg_node);
    alea_node_id_t bottom_right_region = alea_intersection(
        sys, alea_surface_at(sys, split)->neg_node,
        alea_surface_at(sys, bottom_right)->pos_node);
    alea_node_id_t top_left_region = alea_intersection(
        sys, alea_surface_at(sys, split)->pos_node,
        alea_surface_at(sys, top_left)->neg_node);
    alea_node_id_t top_right_region = alea_intersection(
        sys, alea_surface_at(sys, split)->pos_node,
        alea_surface_at(sys, top_right)->pos_node);
    ASSERT(alea_add_cell(sys, 131, bottom_left_region,
                         ALEA_MATERIAL_VOID, 0.0, 0) >= 0);
    ASSERT(alea_add_cell(sys, 132, bottom_right_region,
                         ALEA_MATERIAL_VOID, 0.0, 0) >= 0);
    ASSERT(alea_add_cell(sys, 133, top_left_region,
                         ALEA_MATERIAL_VOID, 0.0, 0) >= 0);
    ASSERT(alea_add_cell(sys, 134, top_right_region,
                         ALEA_MATERIAL_VOID, 0.0, 0) >= 0);
    alea_slice_view_t view;
    alea_slice_view_init(&view, 0.0, 0.0, 0.0,
                         0.0, 0.0, 1.0, 0.0, 1.0, 0.0,
                         -2.0, 2.0, -1.0, 1.0);
    alea_transition_slice_options_t options;
    alea_transition_slice_options_init(&options);
    options.horizontal_rays = 2;
    options.vertical_rays = 0;
    alea_transition_slice_result_t* result =
        alea_transition_slice_result_create();
    ASSERT_NOT_NULL(result);
    ASSERT_EQ(alea_transition_slice_screen(sys, &view, &options, result), 0);
    alea_transition_slice_stats_t stats;
    ASSERT_EQ(alea_transition_slice_stats(result, &stats), 0);
    ASSERT(stats.complete);
    ASSERT_EQ(alea_transition_slice_finding_count(result), (size_t)2);
    ASSERT_EQ(alea_transition_slice_component_count(result), (size_t)2);

    options.enable_critical_refinement = 1;
    options.horizontal_rays = 3;
    options.refine_signals = ALEA_TRANSITION_SLICE_REFINE_FINDING;
    options.max_refinement_frontiers = 1;
    ASSERT_EQ(alea_transition_slice_screen(sys, &view, &options, result), 0);
    ASSERT_EQ(alea_transition_slice_stats(result, &stats), 0);
    ASSERT(!stats.critical_complete);
    ASSERT_EQ(stats.critical_stop_reason,
              ALEA_TRANSITION_SLICE_CRITICAL_MAX_FRONTIERS);
    ASSERT_EQ(alea_transition_slice_refinement_frontier_count(result),
              (size_t)1);
    ASSERT(stats.omitted_refinement_frontiers > 0);
    options.horizontal_rays = 2;
    options.refine_signals = 0;
    options.max_refinement_frontiers = 1024;
    ASSERT_EQ(alea_transition_slice_screen(sys, &view, &options, result), 0);
    ASSERT_EQ(alea_transition_slice_stats(result, &stats), 0);
    ASSERT(stats.critical_complete);
    ASSERT_EQ(alea_transition_slice_critical_tile_count(result), (size_t)2);
    ASSERT_EQ(alea_transition_slice_critical_tile_source_count(result),
              (size_t)2);
    options.max_critical_tiles = 1;
    ASSERT_EQ(alea_transition_slice_screen(sys, &view, &options, result), 0);
    ASSERT_EQ(alea_transition_slice_stats(result, &stats), 0);
    ASSERT(!stats.complete);
    ASSERT(!stats.critical_complete);
    ASSERT_EQ(stats.critical_stop_reason,
              ALEA_TRANSITION_SLICE_CRITICAL_MAX_TILES);
    ASSERT_EQ(alea_transition_slice_critical_tile_count(result), (size_t)1);
    options.max_critical_tiles = 256;
    options.max_critical_tile_sources = 1;
    ASSERT_EQ(alea_transition_slice_screen(sys, &view, &options, result), 0);
    ASSERT_EQ(alea_transition_slice_stats(result, &stats), 0);
    ASSERT(!stats.critical_complete);
    ASSERT_EQ(stats.critical_stop_reason,
              ALEA_TRANSITION_SLICE_CRITICAL_MAX_TILE_SOURCES);
    ASSERT_EQ(alea_transition_slice_critical_tile_source_count(result),
              (size_t)1);

    options.max_components = 1;
    ASSERT_EQ(alea_transition_slice_screen(sys, &view, &options, result), 0);
    ASSERT_EQ(alea_transition_slice_stats(result, &stats), 0);
    ASSERT(!stats.complete);
    ASSERT_EQ(stats.stop_reason, ALEA_TRANSITION_SLICE_STOP_MAX_COMPONENTS);
    ASSERT_EQ(alea_transition_slice_finding_count(result), (size_t)2);
    ASSERT_EQ(alea_transition_slice_component_count(result), (size_t)1);
    alea_transition_slice_result_destroy(result);
    alea_destroy(sys);
}

TEST(transition_slice_screen_streams_bounded_point_coverage_findings) {
    alea_slice_view_t view;
    alea_slice_view_init(&view, 0.0, 0.0, 0.0,
                         0.0, 0.0, 1.0, 0.0, 1.0, 0.0,
                         -4.0, 4.0, -1.0, 1.0);
    alea_transition_slice_options_t options;
    alea_transition_slice_options_init(&options);
    options.horizontal_rays = 1;
    options.vertical_rays = 0;
    options.coverage_uniform_probes_per_ray = 17;
    options.max_coverage_probes = 100;

    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);
    int outer = alea_sphere_surface(sys, 591, 0.0, 0.0, 0.0, 3.0);
    int inner = alea_sphere_surface(sys, 592, 0.0, 0.0, 0.0, 1.0);
    ASSERT(outer >= 0 && inner >= 0);
    ASSERT(alea_add_cell(sys, 141, alea_surface_at(sys, outer)->neg_node,
                         ALEA_MATERIAL_VOID, 0.0, 0) >= 0);
    ASSERT(alea_add_cell(sys, 142, alea_surface_at(sys, inner)->neg_node,
                         ALEA_MATERIAL_VOID, 0.0, 0) >= 0);
    alea_transition_slice_result_t* result =
        alea_transition_slice_result_create();
    ASSERT_NOT_NULL(result);
    ASSERT_EQ(alea_transition_slice_screen(sys, &view, &options, result), 0);
    alea_transition_slice_stats_t stats;
    ASSERT_EQ(alea_transition_slice_stats(result, &stats), 0);
    ASSERT(stats.complete);
    ASSERT_EQ(stats.coverage_probes, (size_t)17);
    ASSERT(stats.unique_coverage_probes > 0);
    ASSERT(stats.skipped_unowned_coverage_probes > 0);
    ASSERT(alea_transition_slice_coverage_finding_count(result) > 0);
    ASSERT_EQ(alea_transition_slice_coverage_component_count(result), (size_t)1);
    alea_transition_slice_coverage_component_t component;
    ASSERT_EQ(alea_transition_slice_coverage_component_get(
                  result, 0, &component), 0);
    ASSERT_EQ(component.kind, ALEA_POINT_COVERAGE_OVERLAP);
    ASSERT(component.finding_count > 0);
    for (size_t i = 0;
         i < alea_transition_slice_coverage_finding_count(result); i++) {
        alea_transition_slice_coverage_finding_t finding;
        ASSERT_EQ(alea_transition_slice_coverage_finding_get(
                      result, i, &finding), 0);
        ASSERT_EQ(finding.kind, ALEA_POINT_COVERAGE_OVERLAP);
        ASSERT_EQ(finding.owner_count, (size_t)2);
        ASSERT_EQ(finding.owner_count_lower_bound, (size_t)2);
        ASSERT(finding.bracket_t_exit > finding.bracket_t_enter);
    }

    options.coverage_uniform_probes_per_ray = 0;
    options.coverage_probe_selected_intervals = 1;
    ASSERT_EQ(alea_transition_slice_screen(sys, &view, &options, result), 0);
    ASSERT_EQ(alea_transition_slice_stats(result, &stats), 0);
    ASSERT(stats.complete);
    ASSERT(stats.coverage_probes > 0);
    ASSERT(alea_transition_slice_coverage_finding_count(result) > 0);

    options.coverage_uniform_probes_per_ray = 17;
    options.coverage_probe_selected_intervals = 0;
    options.max_coverage_probes = 1;
    ASSERT_EQ(alea_transition_slice_screen(sys, &view, &options, result), 0);
    ASSERT_EQ(alea_transition_slice_stats(result, &stats), 0);
    ASSERT(!stats.complete);
    ASSERT_EQ(stats.stop_reason,
              ALEA_TRANSITION_SLICE_STOP_MAX_COVERAGE_PROBES);
    ASSERT_EQ(stats.coverage_probes, (size_t)1);
    alea_transition_slice_result_destroy(result);
    alea_destroy(sys);

    /* Unowned samples are opt-in because an arbitrary viewport normally also
     * contains legitimate exterior void.  With the policy enabled, a strip
     * between two halfspaces is retained as sampled gap evidence. */
    sys = alea_create();
    ASSERT_NOT_NULL(sys);
    int left = alea_plane_surface(sys, 593, 1.0, 0.0, 0.0, 0.5);
    int right = alea_plane_surface(sys, 594, 1.0, 0.0, 0.0, -0.5);
    ASSERT(left >= 0 && right >= 0);
    ASSERT(alea_add_cell(sys, 143, alea_surface_at(sys, left)->neg_node,
                         ALEA_MATERIAL_VOID, 0.0, 0) >= 0);
    ASSERT(alea_add_cell(sys, 144, alea_surface_at(sys, right)->pos_node,
                         ALEA_MATERIAL_VOID, 0.0, 0) >= 0);
    alea_transition_slice_options_init(&options);
    options.horizontal_rays = 1;
    options.vertical_rays = 0;
    options.coverage_uniform_probes_per_ray = 17;
    options.report_unowned_coverage = 1;
    result = alea_transition_slice_result_create();
    ASSERT_NOT_NULL(result);
    ASSERT_EQ(alea_transition_slice_screen(sys, &view, &options, result), 0);
    ASSERT(alea_transition_slice_coverage_finding_count(result) > 0);
    int saw_gap = 0;
    for (size_t i = 0;
         i < alea_transition_slice_coverage_finding_count(result); i++) {
        alea_transition_slice_coverage_finding_t finding;
        ASSERT_EQ(alea_transition_slice_coverage_finding_get(
                      result, i, &finding), 0);
        if (finding.kind == ALEA_POINT_COVERAGE_GAP &&
            finding.world_point[0] > -0.5 && finding.world_point[0] < 0.5)
            saw_gap = 1;
    }
    ASSERT(saw_gap);
    alea_transition_slice_result_destroy(result);
    alea_destroy(sys);

    /* A card-distinct transition into an overlapped interval produces both
     * verdicts.  The link identifies the shared selected-event boundary
     * without a geometric nearest-neighbor search. */
    sys = alea_create();
    ASSERT_NOT_NULL(sys);
    int primary = alea_plane_surface(sys, 595, 1.0, 0.0, 0.0, 0.0);
    int after_a = alea_plane_surface(sys, 596, 1.0, 0.0, 0.0, 0.0);
    int after_b = alea_plane_surface(sys, 597, 1.0, 0.0, 0.0, 0.0);
    ASSERT(primary >= 0 && after_a >= 0 && after_b >= 0);
    ASSERT(alea_add_cell(sys, 145, alea_surface_at(sys, primary)->neg_node,
                         ALEA_MATERIAL_VOID, 0.0, 0) >= 0);
    ASSERT(alea_add_cell(sys, 146, alea_surface_at(sys, after_a)->pos_node,
                         ALEA_MATERIAL_VOID, 0.0, 0) >= 0);
    ASSERT(alea_add_cell(sys, 147, alea_surface_at(sys, after_b)->pos_node,
                         ALEA_MATERIAL_VOID, 0.0, 0) >= 0);
    alea_transition_slice_options_init(&options);
    options.horizontal_rays = 1;
    options.vertical_rays = 0;
    options.coverage_probe_selected_intervals = 1;
    result = alea_transition_slice_result_create();
    ASSERT_NOT_NULL(result);
    ASSERT_EQ(alea_transition_slice_screen(sys, &view, &options, result), 0);
    ASSERT_EQ(alea_transition_slice_stats(result, &stats), 0);
    ASSERT(stats.complete);
    ASSERT_EQ(alea_transition_slice_component_count(result), (size_t)1);
    ASSERT_EQ(alea_transition_slice_coverage_component_count(result),
              (size_t)1);
    ASSERT_EQ(alea_transition_slice_component_link_count(result), (size_t)1);
    ASSERT_EQ(stats.component_links, (size_t)1);
    alea_transition_slice_component_link_t link;
    ASSERT_EQ(alea_transition_slice_component_link_get(result, 0, &link), 0);
    ASSERT_EQ(link.transition_component_index, (size_t)0);
    ASSERT_EQ(link.coverage_component_index, (size_t)0);
    ASSERT(link.boundary_sides & ALEA_TRANSITION_SLICE_LINK_ENTER);
    ASSERT_EQ(link.witness_pair_count, (size_t)1);
    alea_transition_slice_result_destroy(result);
    alea_destroy(sys);

    sys = alea_create();
    ASSERT_NOT_NULL(sys);
    int split = alea_plane_surface(sys, 598, 0.0, 1.0, 0.0, 0.0);
    int x_planes[6];
    for (int i = 0; i < 6; i++)
        x_planes[i] = alea_plane_surface(
            sys, 599 + i, 1.0, 0.0, 0.0, 0.0);
    ASSERT(split >= 0);
    for (int i = 0; i < 6; i++) ASSERT(x_planes[i] >= 0);
    for (int half = 0; half < 2; half++) {
        alea_node_id_t split_node = half == 0
            ? alea_surface_at(sys, split)->neg_node
            : alea_surface_at(sys, split)->pos_node;
        for (int region = 0; region < 3; region++) {
            int plane = x_planes[half * 3 + region];
            alea_node_id_t x_node = region == 0
                ? alea_surface_at(sys, plane)->neg_node
                : alea_surface_at(sys, plane)->pos_node;
            ASSERT(alea_add_cell(
                sys, 148 + half * 3 + region,
                alea_intersection(sys, split_node, x_node),
                ALEA_MATERIAL_VOID, 0.0, 0) >= 0);
        }
    }
    alea_transition_slice_options_init(&options);
    options.horizontal_rays = 2;
    options.vertical_rays = 0;
    options.coverage_probe_selected_intervals = 1;
    result = alea_transition_slice_result_create();
    ASSERT_NOT_NULL(result);
    ASSERT_EQ(alea_transition_slice_screen(sys, &view, &options, result), 0);
    ASSERT_EQ(alea_transition_slice_component_link_count(result), (size_t)2);
    options.max_component_links = 1;
    ASSERT_EQ(alea_transition_slice_screen(sys, &view, &options, result), 0);
    ASSERT_EQ(alea_transition_slice_stats(result, &stats), 0);
    ASSERT(!stats.complete);
    ASSERT_EQ(stats.stop_reason,
              ALEA_TRANSITION_SLICE_STOP_MAX_COMPONENT_LINKS);
    ASSERT_EQ(alea_transition_slice_component_link_count(result), (size_t)1);
    alea_transition_slice_result_destroy(result);
    alea_destroy(sys);
}

/* Regression: terminal segment must clamp t_exit to max_distance when the
 * ray's endpoint lies inside a cell. Previously t_exit was DBL_MAX (or the
 * geometric exit beyond max_distance). */
TEST(raycast_clamps_terminal_segment_to_t_max) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);

    int surf_idx = alea_sphere_surface(sys, 1, 0, 0, 0, 5.0);
    alea_node_id_t sphere = alea_surface_at(sys, surf_idx)->neg_node;
    int m1 = alea_add_material(sys, 1);
    alea_add_cell(sys, 1, sphere, m1, -2.7, 0);
    ASSERT_EQ(alea_prepare_query_acceleration(sys), 0);

    /* Endpoint (t=10) is at the sphere center, well inside the cell.
     * Geometric exit at t=15 lies past max_distance=10. */
    alea_raycast_result_t result;
    alea_raycast_result_init(&result);

    ASSERT_EQ(alea_raycast(sys, -10, 0, 0, 1, 0, 0, 10.0, &result), 0);
    ASSERT(result.segments.count >= 2);
    const alea_ray_segment_t* last = &result.segments.data[result.segments.count - 1];
    ASSERT_EQ(last->cell_id, 1);
    ASSERT_NEAR(last->t_enter, 5.0, EPS);
    ASSERT_NEAR(last->t_exit, 10.0, EPS);

    /* Same expectation for cell-aware path. */
    alea_raycast_result_clear(&result);
    ASSERT_EQ(alea_raycast_cell_aware(sys, -10, 0, 0, 1, 0, 0, 10.0, &result), 0);
    ASSERT(result.segments.count >= 2);
    last = &result.segments.data[result.segments.count - 1];
    ASSERT_EQ(last->cell_id, 1);
    ASSERT_NEAR(last->t_enter, 5.0, EPS);
    ASSERT_NEAR(last->t_exit, 10.0, EPS);

    alea_raycast_result_free(&result);
    alea_destroy(sys);
}

TEST(compact_hierarchical_batch_matches_single_ray_segments) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);

    int surf_idx = alea_sphere_surface(sys, 1, 0, 0, 0, 5.0);
    ASSERT(surf_idx >= 0);
    alea_node_id_t sphere = alea_surface_at(sys, surf_idx)->neg_node;
    int material = alea_add_material(sys, 7);
    ASSERT(material >= 0);
    ASSERT(alea_add_cell(sys, 10, sphere, material, -1.5, 0) >= 0);

    const double origins[] = {
        -10.0, 0.0, 0.0,
         10.0, 0.0, 0.0,
          0.0, 0.0, 0.0
    };
    const double directions[] = {
        1.0, 0.0, 0.0,
       -1.0, 0.0, 0.0,
        1.0, 0.0, 0.0
    };
    alea_raycast_batch_options_t options = {
        .struct_size = sizeof(options),
        .fields = ALEA_RAY_BATCH_MATERIAL | ALEA_RAY_BATCH_DENSITY |
                  ALEA_RAY_BATCH_SURFACES | ALEA_RAY_BATCH_RESOLUTION_FLAGS |
                  ALEA_RAY_BATCH_PROJECTED_OWNER | ALEA_RAY_BATCH_FULL_PATHS,
        .projected_depth = -1,
        .max_segments = 0,
        .max_path_entries = 0
    };
    alea_raycast_batch_result_t* batch = alea_raycast_batch_result_create();
    ASSERT_NOT_NULL(batch);
    ASSERT_EQ(alea_raycast_hier_batch(sys, origins, directions, 3, 30.0,
                                      &options, batch), 0);
    ASSERT_EQ(alea_raycast_batch_ray_count(batch), 3);
    ASSERT_NOT_NULL(alea_raycast_batch_ray_offsets(batch));
    ASSERT_EQ(alea_raycast_batch_ray_offsets(batch)[0], 0);
    ASSERT_EQ(alea_raycast_batch_ray_offsets(batch)[3],
              alea_raycast_batch_segment_count(batch));
    ASSERT_NOT_NULL(alea_raycast_batch_material_ids(batch));
    ASSERT_NOT_NULL(alea_raycast_batch_densities(batch));
    ASSERT_NOT_NULL(alea_raycast_batch_enter_surface_ids(batch));
    ASSERT_NOT_NULL(alea_raycast_batch_exit_surface_ids(batch));
    ASSERT_NOT_NULL(alea_raycast_batch_resolution_flags(batch));
    ASSERT_NOT_NULL(alea_raycast_batch_projected_cell_ids(batch));
    ASSERT_NOT_NULL(alea_raycast_batch_projected_material_ids(batch));
    ASSERT_NOT_NULL(alea_raycast_batch_segment_path_offsets(batch));
    ASSERT_NOT_NULL(alea_raycast_batch_path_cell_ids(batch));

    const uint64_t* offsets = alea_raycast_batch_ray_offsets(batch);
    const double* t_enter = alea_raycast_batch_t_enter(batch);
    const double* t_exit = alea_raycast_batch_t_exit(batch);
    const int32_t* cell_ids = alea_raycast_batch_cell_ids(batch);
    const int32_t* material_ids = alea_raycast_batch_material_ids(batch);
    const double* densities = alea_raycast_batch_densities(batch);
    const int32_t* enter_surfaces = alea_raycast_batch_enter_surface_ids(batch);
    const int32_t* exit_surfaces = alea_raycast_batch_exit_surface_ids(batch);
    const uint8_t* flags = alea_raycast_batch_resolution_flags(batch);
    const int32_t* projected_cells = alea_raycast_batch_projected_cell_ids(batch);
    const int32_t* projected_materials = alea_raycast_batch_projected_material_ids(batch);
    const uint64_t* path_offsets = alea_raycast_batch_segment_path_offsets(batch);
    const int32_t* path_cells = alea_raycast_batch_path_cell_ids(batch);

    for (size_t i = 0; i < 3; i++) {
        alea_raycast_result_t single;
        alea_raycast_result_init(&single);
        ASSERT_EQ(alea_raycast_hier_fast_segments(
                      sys, origins[i * 3], origins[i * 3 + 1], origins[i * 3 + 2],
                      directions[i * 3], directions[i * 3 + 1], directions[i * 3 + 2],
                      30.0, &single), 0);
        ASSERT_EQ(offsets[i + 1] - offsets[i], single.segments.count);
        for (size_t j = 0; j < single.segments.count; j++) {
            const alea_ray_segment_t* seg = &single.segments.data[j];
            size_t k = (size_t)offsets[i] + j;
            ASSERT_NEAR(t_enter[k], seg->t_enter, EPS);
            ASSERT_NEAR(t_exit[k], seg->t_exit, EPS);
            ASSERT_EQ(cell_ids[k], seg->cell_id);
            ASSERT_EQ(material_ids[k], seg->material_id);
            ASSERT_NEAR(densities[k], seg->density, EPS);
            ASSERT_EQ(enter_surfaces[k], seg->enter_surface_id);
            ASSERT_EQ(exit_surfaces[k], seg->exit_surface_id);
            ASSERT_EQ(flags[k], seg->resolution_flags);
            ASSERT_EQ(projected_cells[k], seg->cell_id);
            ASSERT_EQ(projected_materials[k], seg->material_id);
            if (path_offsets[k + 1] > path_offsets[k]) {
                ASSERT_EQ(path_cells[path_offsets[k + 1] - 1], projected_cells[k]);
            } else {
                ASSERT_EQ(projected_cells[k], -1);
            }
        }
        alea_raycast_result_free(&single);
    }

    size_t old_count = alea_raycast_batch_segment_count(batch);
    options.max_segments = 1;
    ASSERT_EQ(alea_raycast_hier_batch(sys, origins, directions, 3, 30.0,
                                      &options, batch), -1);
    ASSERT_EQ(alea_raycast_batch_segment_count(batch), old_count);

    options.max_segments = 0;
    options.max_path_entries = 1;
    ASSERT_EQ(alea_raycast_hier_batch(sys, origins, directions, 3, 30.0,
                                      &options, batch), -1);
    ASSERT_EQ(alea_raycast_batch_segment_count(batch), old_count);

    options.max_path_entries = 0;
    options.max_output_bytes = 1;
    ASSERT_EQ(alea_raycast_hier_batch(sys, origins, directions, 3, 30.0,
                                      &options, batch), -1);
    ASSERT_EQ(alea_raycast_batch_segment_count(batch), old_count);

    alea_raycast_batch_result_destroy(batch);
    alea_destroy(sys);
}

TEST(compact_ray_slice_returns_view_u_coordinates) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);
    int surf_idx = alea_sphere_surface(sys, 1, 0, 0, 0, 5.0);
    ASSERT(surf_idx >= 0);
    int material = alea_add_material(sys, 1);
    ASSERT(material >= 0);
    ASSERT(alea_add_cell(sys, 10, alea_surface_at(sys, surf_idx)->neg_node,
                         material, -1.0, 0) >= 0);

    alea_slice_view_t view;
    alea_slice_view_init(&view, 0.0, 0.0, 0.0,
                         1.0, 1.0, 1.0,
                         0.0, 0.0, 1.0,
                         -10.0, 10.0, -1.0, 1.0);
    alea_raycast_batch_options_t options = {
        .struct_size = sizeof(options),
        .fields = ALEA_RAY_BATCH_MATERIAL,
        .projected_depth = -1,
        .max_segments = 0,
        .max_path_entries = 0,
        .max_output_bytes = 0
    };
    alea_raycast_batch_result_t* batch = alea_raycast_batch_result_create();
    ASSERT_NOT_NULL(batch);
    ASSERT_EQ(alea_trace_ray_slice_compact(sys, &view, 2, &options, batch), 0);
    ASSERT_EQ(alea_raycast_batch_ray_count(batch), 2);

    const uint64_t* offsets = alea_raycast_batch_ray_offsets(batch);
    const double* u_enter = alea_raycast_batch_t_enter(batch);
    const double* u_exit = alea_raycast_batch_t_exit(batch);
    const int32_t* cell_ids = alea_raycast_batch_cell_ids(batch);
    for (size_t row = 0; row < 2; row++) {
        int found = 0;
        for (size_t i = (size_t)offsets[row]; i < (size_t)offsets[row + 1]; i++) {
            ASSERT(u_enter[i] >= view.u_min - EPS);
            ASSERT(u_exit[i] <= view.u_max + EPS);
            if (cell_ids[i] == 10) {
                found = 1;
                ASSERT(u_enter[i] < -4.8);
                ASSERT(u_exit[i] > 4.8);
            }
        }
        ASSERT(found);
    }

    alea_raycast_batch_result_destroy(batch);
    alea_destroy(sys);
}

TEST(compact_batch_surface_fields_preserve_selected_segments) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);
    int surface = alea_sphere_surface(sys, 1, 0, 0, 0, 5.0);
    int material = alea_add_material(sys, 7);
    ASSERT(surface >= 0 && material >= 0);
    ASSERT(alea_add_cell(sys, 10, alea_surface_at(sys, surface)->neg_node,
                         material, -1.5, 0) >= 0);

    const double origins[] = {-10, 0, 0, 10, 0, 0};
    const double directions[] = {1, 0, 0, -1, 0, 0};
    const alea_raycast_batch_options_t surface_off = {
        .struct_size = sizeof(surface_off),
        .fields = ALEA_RAY_BATCH_MATERIAL
    };
    const alea_raycast_batch_options_t surface_on = {
        .struct_size = sizeof(surface_on),
        .fields = ALEA_RAY_BATCH_MATERIAL | ALEA_RAY_BATCH_SURFACES
    };
    alea_raycast_batch_result_t* off = alea_raycast_batch_result_create();
    alea_raycast_batch_result_t* on = alea_raycast_batch_result_create();
    ASSERT_NOT_NULL(off);
    ASSERT_NOT_NULL(on);
    ASSERT_EQ(alea_raycast_hier_batch(sys, origins, directions, 2, 30,
                                      &surface_off, off), 0);
    ASSERT_EQ(alea_raycast_hier_batch(sys, origins, directions, 2, 30,
                                      &surface_on, on), 0);

    alea_raycast_batch_work_stats_t work_stats;
    ASSERT_EQ(alea_raycast_batch_result_get_work_stats_internal(on, &work_stats), 0);
    ASSERT_EQ(work_stats.peak_trace_staging_bytes, 0);
    ASSERT_EQ(work_stats.total_result_buffer_growths, 0);

    ASSERT_NULL(alea_raycast_batch_enter_surface_ids(off));
    ASSERT_NOT_NULL(alea_raycast_batch_enter_surface_ids(on));
    ASSERT_EQ(alea_raycast_batch_ray_count(off), alea_raycast_batch_ray_count(on));
    ASSERT_EQ(alea_raycast_batch_segment_count(off),
              alea_raycast_batch_segment_count(on));
    const uint64_t* off_offsets = alea_raycast_batch_ray_offsets(off);
    const uint64_t* on_offsets = alea_raycast_batch_ray_offsets(on);
    const double* off_enter = alea_raycast_batch_t_enter(off);
    const double* on_enter = alea_raycast_batch_t_enter(on);
    const double* off_exit = alea_raycast_batch_t_exit(off);
    const double* on_exit = alea_raycast_batch_t_exit(on);
    const int32_t* off_cells = alea_raycast_batch_cell_ids(off);
    const int32_t* on_cells = alea_raycast_batch_cell_ids(on);
    for (size_t i = 0; i <= 2; i++)
        ASSERT_EQ(off_offsets[i], on_offsets[i]);
    for (size_t i = 0; i < alea_raycast_batch_segment_count(off); i++) {
        ASSERT_NEAR(off_enter[i], on_enter[i], EPS);
        ASSERT_NEAR(off_exit[i], on_exit[i], EPS);
        ASSERT_EQ(off_cells[i], on_cells[i]);
    }

    /* The common worker-arena route has the same operation-wide limits and
     * transactional lifetime contract as the legacy rich-trace route. */
    const size_t segment_count = alea_raycast_batch_segment_count(on);
    ASSERT(segment_count > 0);
    alea_raycast_batch_options_t limited = surface_on;
    limited.max_segments = segment_count;
    ASSERT_EQ(alea_raycast_hier_batch(sys, origins, directions, 2, 30,
                                      &limited, on), 0);
    const uint64_t* previous_offsets = alea_raycast_batch_ray_offsets(on);
    limited.max_segments = segment_count - 1;
    ASSERT_EQ(alea_raycast_hier_batch(sys, origins, directions, 2, 30,
                                      &limited, on), -1);
    ASSERT_EQ(alea_raycast_batch_ray_offsets(on), previous_offsets);
    ASSERT_EQ(alea_raycast_batch_segment_count(on), segment_count);

    const size_t exact_bytes = 3 * sizeof(uint64_t) + segment_count *
        (2 * sizeof(double) + 4 * sizeof(int32_t));
    limited.max_segments = 0;
    limited.max_output_bytes = exact_bytes;
    ASSERT_EQ(alea_raycast_hier_batch(sys, origins, directions, 2, 30,
                                      &limited, on), 0);
    previous_offsets = alea_raycast_batch_ray_offsets(on);
    limited.max_output_bytes = exact_bytes - 1;
    ASSERT_EQ(alea_raycast_hier_batch(sys, origins, directions, 2, 30,
                                      &limited, on), -1);
    ASSERT_EQ(alea_raycast_batch_ray_offsets(on), previous_offsets);
    ASSERT_EQ(alea_raycast_batch_segment_count(on), segment_count);

    alea_raycast_batch_result_destroy(on);
    alea_raycast_batch_result_destroy(off);
    alea_destroy(sys);
}

/* The same child universe occurs twice through different transformed fills.
 * This exercises hierarchy-path CSR and transformed placement metadata, not
 * just flat segment values. */
TEST(compact_hierarchical_batch_matches_transformed_fill_paths) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);
    int container_material = alea_add_material(sys, 1);
    int child_material = alea_add_material(sys, 2);
    int left_surface = alea_sphere_surface(sys, 1, -5.0, 0.0, 0.0, 3.0);
    int right_surface = alea_sphere_surface(sys, 2, 5.0, 0.0, 0.0, 3.0);
    int child_surface = alea_sphere_surface(sys, 3, 0.0, 0.0, 0.0, 1.0);
    ASSERT(left_surface >= 0 && right_surface >= 0 && child_surface >= 0);
    int left = alea_add_cell(sys, 10, alea_surface_at(sys, left_surface)->neg_node,
                             container_material, -1.0, 0);
    int right = alea_add_cell(sys, 20, alea_surface_at(sys, right_surface)->neg_node,
                              container_material, -1.0, 0);
    ASSERT(left >= 0 && right >= 0);
    ASSERT(alea_add_cell(sys, 30, alea_surface_at(sys, child_surface)->neg_node,
                         child_material, -2.0, 10) >= 0);
    /* A filled universe must cover the parent placement. The explicit
     * complement makes its non-fuel region a resolved void instead of an
     * undefined-fill fallback, so the ray steps across the child boundaries. */
    ASSERT(alea_add_cell(sys, 31, alea_surface_at(sys, child_surface)->pos_node,
                         ALEA_MATERIAL_VOID, 0.0, 10) >= 0);
    const double left_translation[3] = {-5.0, 0.0, 0.0};
    const double right_translation[3] = {5.0, 0.0, 0.0};
    ASSERT_EQ(alea_add_transform(sys, 101, left_translation, 3, 0), 0);
    ASSERT_EQ(alea_add_transform(sys, 102, right_translation, 3, 0), 0);
    ASSERT_EQ(alea_set_fill(sys, left, 10, 101), 0);
    ASSERT_EQ(alea_set_fill(sys, right, 10, 102), 0);

    const double origins[] = {-10.0, 0.0, 0.0, -10.0, 2.0, 0.0};
    const double directions[] = {1.0, 0.0, 0.0, 1.0, 0.0, 0.0};
    alea_raycast_batch_options_t options = {
        .struct_size = sizeof(options),
        .fields = ALEA_RAY_BATCH_MATERIAL | ALEA_RAY_BATCH_DENSITY |
                  ALEA_RAY_BATCH_SURFACES | ALEA_RAY_BATCH_RESOLUTION_FLAGS |
                  ALEA_RAY_BATCH_PROJECTED_OWNER | ALEA_RAY_BATCH_FULL_PATHS,
        .projected_depth = -1
    };
    alea_raycast_batch_result_t* batch = alea_raycast_batch_result_create();
    ASSERT_NOT_NULL(batch);
    ASSERT_EQ(alea_raycast_hier_batch(sys, origins, directions, 2, 25.0,
                                      &options, batch), 0);

    const uint64_t* offsets = alea_raycast_batch_ray_offsets(batch);
    const uint64_t* path_offsets = alea_raycast_batch_segment_path_offsets(batch);
    const int32_t* projected_cells = alea_raycast_batch_projected_cell_ids(batch);
    const uint64_t* projected_keys = alea_raycast_batch_projected_occurrence_keys(batch);
    const int32_t* path_cells = alea_raycast_batch_path_cell_ids(batch);
    const int32_t* path_materials = alea_raycast_batch_path_material_ids(batch);
    const int32_t* path_universes = alea_raycast_batch_path_universe_ids(batch);
    const int32_t* path_fills = alea_raycast_batch_path_fill_universes(batch);
    const int32_t* path_depths = alea_raycast_batch_path_depths(batch);
    const uint8_t* path_lattice = alea_raycast_batch_path_is_lattice(batch);
    const double* path_origins = alea_raycast_batch_path_lattice_origins_xyz(batch);
    const uint64_t* path_keys = alea_raycast_batch_path_occurrence_keys(batch);
    uint64_t left_child_key = 0;
    uint64_t right_child_key = 0;

    for (size_t ray_index = 0; ray_index < 2; ray_index++) {
        alea_raycast_result_t single;
        alea_raycast_result_init(&single);
        alea_raycast_result_set_path_capture(&single, 1);
        ASSERT_EQ(alea_raycast_hier_fast_segments(
                      sys, origins[ray_index * 3], origins[ray_index * 3 + 1],
                      origins[ray_index * 3 + 2], directions[ray_index * 3],
                      directions[ray_index * 3 + 1], directions[ray_index * 3 + 2],
                      25.0, &single), 0);
        ASSERT_EQ(offsets[ray_index + 1] - offsets[ray_index], single.segments.count);
        for (size_t j = 0; j < single.segments.count; j++) {
            size_t segment = (size_t)offsets[ray_index] + j;
            const alea_ray_segment_t* expected = &single.segments.data[j];
            size_t count = alea_raycast_segment_path_count(&single, j);
            ASSERT_EQ(path_offsets[segment + 1] - path_offsets[segment], count);
            ASSERT_EQ(projected_cells[segment], expected->cell_id);
            for (size_t p = 0; p < count; p++) {
                alea_raycast_path_entry_t entry;
                size_t q = (size_t)path_offsets[segment] + p;
                ASSERT_EQ(alea_raycast_segment_path_get(&single, j, p, &entry), 0);
                ASSERT_EQ(path_cells[q], entry.cell_id);
                ASSERT_EQ(path_materials[q], entry.material_id);
                ASSERT_EQ(path_universes[q], entry.universe_id);
                ASSERT_EQ(path_fills[q], entry.fill_universe);
                ASSERT_EQ(path_depths[q], entry.depth);
                ASSERT_EQ(path_lattice[q], entry.is_lattice);
                ASSERT_NEAR(path_origins[q * 3], entry.lattice_origin[0], EPS);
                ASSERT_NEAR(path_origins[q * 3 + 1], entry.lattice_origin[1], EPS);
                ASSERT_NEAR(path_origins[q * 3 + 2], entry.lattice_origin[2], EPS);
                ASSERT_EQ(path_keys[q], entry.occurrence_key);
            }
            if (count > 0) {
                size_t leaf = (size_t)path_offsets[segment + 1] - 1;
                ASSERT_EQ(projected_keys[segment], path_keys[leaf]);
                if (ray_index == 0 && expected->cell_id == 30) {
                    if (expected->t_enter < 10.0)
                        left_child_key = projected_keys[segment];
                    else
                        right_child_key = projected_keys[segment];
                }
            }
        }
        alea_raycast_result_free(&single);
    }
    ASSERT(left_child_key != 0);
    ASSERT(right_child_key != 0);
    ASSERT_NE(left_child_key, right_child_key);
    alea_raycast_batch_result_destroy(batch);
    alea_destroy(sys);
}

TEST(selected_boundary_events_retain_transformed_occurrence_receipts) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);
    int container_material = alea_add_material(sys, 1);
    int child_material = alea_add_material(sys, 2);
    int left_surface = alea_sphere_surface(sys, 11, -5.0, 0.0, 0.0, 3.0);
    int right_surface = alea_sphere_surface(sys, 12, 5.0, 0.0, 0.0, 3.0);
    int child_surface = alea_sphere_surface(sys, 13, 0.0, 0.0, 0.0, 1.0);
    ASSERT(left_surface >= 0 && right_surface >= 0 && child_surface >= 0);
    int left = alea_add_cell(sys, 110,
                             alea_surface_at(sys, left_surface)->neg_node,
                             container_material, -1.0, 0);
    int right = alea_add_cell(sys, 120,
                              alea_surface_at(sys, right_surface)->neg_node,
                              container_material, -1.0, 0);
    ASSERT(left >= 0 && right >= 0);
    ASSERT(alea_add_cell(sys, 130,
                         alea_surface_at(sys, child_surface)->neg_node,
                         child_material, -2.0, 10) >= 0);
    ASSERT(alea_add_cell(sys, 131,
                         alea_surface_at(sys, child_surface)->pos_node,
                         ALEA_MATERIAL_VOID, 0.0, 10) >= 0);
    const double left_translation[3] = {-5.0, 0.0, 0.0};
    const double right_translation[3] = {5.0, 0.0, 0.0};
    ASSERT_EQ(alea_add_transform(sys, 111, left_translation, 3, 0), 0);
    ASSERT_EQ(alea_add_transform(sys, 112, right_translation, 3, 0), 0);
    ASSERT_EQ(alea_set_fill(sys, left, 10, 111), 0);
    ASSERT_EQ(alea_set_fill(sys, right, 10, 112), 0);
    ASSERT_EQ(alea_prepare_query_acceleration(sys), 0);

    alea_ray_t ray;
    ASSERT_EQ(alea_ray_init(&ray, -10.0, 0.0, 0.0, 1.0, 0.0, 0.0), 0);
    alea_raycast_result_t scratch;
    alea_ray_boundary_event_result_t events;
    alea_raycast_result_init(&scratch);
    alea_ray_boundary_event_result_init(&events);
    ASSERT_EQ(alea_raycast_selected_boundary_events_nocache(
                  sys, &ray, 20.0, &scratch, &events), 0);

    uint64_t child_occurrences[2] = {0, 0};
    size_t child_occurrence_count = 0;
    size_t child_event_count = 0;
    for (size_t i = 0; i < events.events.count; i++) {
        const alea_ray_boundary_event_t* event = &events.events.data[i];
        if (event->surface_id != 13 || event->active_universe_id != 10)
            continue;
        child_event_count++;
        ASSERT(event->provenance_flags & ALEA_BOUNDARY_PROVENANCE_ACTIVE_FRAME);
        ASSERT(event->provenance_flags & ALEA_BOUNDARY_PROVENANCE_BEFORE_OWNER);
        ASSERT(event->provenance_flags & ALEA_BOUNDARY_PROVENANCE_AFTER_OWNER);
        ASSERT(event->active_cell_id == 130 || event->active_cell_id == 131);
        ASSERT_EQ(event->active_depth, 1);
        ASSERT(event->active_occurrence_key != 0);
        ASSERT(event->active_parent_occurrence_key != 0);
        ASSERT(event->before_occurrence_key != 0);
        ASSERT(event->after_occurrence_key != 0);
        ASSERT_NE(event->before_occurrence_key, event->after_occurrence_key);
        ASSERT_NEAR(fabs(event->local_point[0]), 1.0, 1e-7);
        ASSERT_NEAR(event->local_point[1], 0.0, EPS);
        ASSERT_NEAR(event->local_point[2], 0.0, EPS);
        ASSERT_NEAR(event->local_direction[0], 1.0, EPS);

        alea_transition_result_t transition;
        ASSERT_EQ(alea_check_selected_boundary_event_transition_nocache(
                      sys, event, NULL, &transition), 0);
        ASSERT_EQ(transition.kind, ALEA_TRANSITION_VALID);
        ASSERT_EQ(transition.universe_id, 10);
        ASSERT_EQ(transition.current_cell_id, event->active_cell_id);
        ASSERT_EQ(transition.occurrence_depth, 1);
        ASSERT_EQ(transition.current_occurrence_key,
                  event->active_occurrence_key);
        ASSERT_EQ(transition.before_occurrence_key,
                  event->before_occurrence_key);
        ASSERT_EQ(transition.selected_after_occurrence_key,
                  event->after_occurrence_key);
        ASSERT_EQ(transition.coverage_fallbacks, (size_t)0);

        if (event->active_cell_id == 130) {
            int seen = 0;
            for (size_t j = 0; j < child_occurrence_count; j++)
                if (child_occurrences[j] == event->active_occurrence_key)
                    seen = 1;
            if (!seen && child_occurrence_count < 2)
                child_occurrences[child_occurrence_count++] =
                    event->active_occurrence_key;
        }
    }
    ASSERT_EQ(child_event_count, (size_t)4);
    ASSERT_EQ(child_occurrence_count, (size_t)2);
    ASSERT_NE(child_occurrences[0], child_occurrences[1]);

    alea_ray_boundary_event_query_result_t* public_events =
        alea_ray_boundary_event_query_result_create();
    ASSERT_NOT_NULL(public_events);
    alea_ray_boundary_event_options_t public_options;
    alea_ray_boundary_event_options_init(&public_options);
    public_options.t_max = 20.0;
    public_options.include_occurrence_provenance = 1;
    ASSERT_EQ(alea_ray_boundary_event_query(
                  sys, -10.0, 0.0, 0.0, 1.0, 0.0, 0.0,
                  &public_options, public_events), 0);
    size_t public_child_events = 0;
    for (size_t i = 0; i < alea_ray_boundary_event_count(public_events); i++) {
        double t = 0.0;
        int kind = -1, surface_id = -1;
        ASSERT_EQ(alea_ray_boundary_event_get(
                      public_events, i, &t, &kind, &surface_id,
                      NULL, NULL, NULL, NULL, NULL, NULL,
                      NULL, NULL, NULL), 0);
        if (kind != ALEA_RAY_EVENT_PHYSICAL || surface_id != 13) continue;
        alea_ray_boundary_event_provenance_t provenance;
        ASSERT_EQ(alea_ray_boundary_event_provenance_get(
                      public_events, i, &provenance), 0);
        ASSERT(provenance.flags &
               ALEA_RAY_BOUNDARY_PROVENANCE_ACTIVE_FRAME);
        ASSERT_EQ(provenance.active_universe_id, 10);
        ASSERT_NEAR(fabs(provenance.local_point[0]), 1.0, 1e-7);
        public_child_events++;
    }
    ASSERT_EQ(public_child_events, (size_t)4);
    alea_ray_boundary_event_query_result_destroy(public_events);

    alea_ray_boundary_event_result_free(&events);
    alea_raycast_result_free(&scratch);
    alea_destroy(sys);
}

TEST(compact_ray_slice_matches_explicit_generic_batch) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);
    int surface = alea_sphere_surface(sys, 1, 0.0, 0.0, 0.0, 5.0);
    int material = alea_add_material(sys, 1);
    ASSERT(surface >= 0 && material >= 0);
    ASSERT(alea_add_cell(sys, 10, alea_surface_at(sys, surface)->neg_node,
                         material, -1.0, 0) >= 0);
    alea_slice_view_t view;
    alea_slice_view_init(&view, 1.0, -2.0, 3.0, 1.0, 2.0, 3.0, 0.0, 0.0, 1.0,
                         -8.0, 9.0, -3.0, 4.0);
    enum { rows = 3 };
    double origins[rows * 3];
    double directions[rows * 3];
    for (size_t row = 0; row < rows; row++) {
        double v = view.v_min + ((double)row + 0.5) *
                   (view.v_max - view.v_min) / (double)rows;
        for (size_t axis = 0; axis < 3; axis++) {
            origins[row * 3 + axis] = view.plane.origin[axis] +
                view.u_min * view.plane.u_axis[axis] + v * view.plane.v_axis[axis];
            directions[row * 3 + axis] = view.plane.u_axis[axis];
        }
    }
    alea_raycast_batch_options_t options = {
        .struct_size = sizeof(options),
        .fields = ALEA_RAY_BATCH_MATERIAL | ALEA_RAY_BATCH_DENSITY |
                  ALEA_RAY_BATCH_SURFACES | ALEA_RAY_BATCH_RESOLUTION_FLAGS
    };
    alea_raycast_batch_result_t* slice = alea_raycast_batch_result_create();
    alea_raycast_batch_result_t* generic = alea_raycast_batch_result_create();
    ASSERT_NOT_NULL(slice);
    ASSERT_NOT_NULL(generic);
    ASSERT_EQ(alea_trace_ray_slice_compact(sys, &view, rows, &options, slice), 0);
    ASSERT_EQ(alea_raycast_hier_batch(sys, origins, directions, rows,
                                      view.u_max - view.u_min, &options, generic), 0);
    ASSERT_EQ(alea_raycast_batch_ray_count(slice), alea_raycast_batch_ray_count(generic));
    ASSERT_EQ(alea_raycast_batch_segment_count(slice), alea_raycast_batch_segment_count(generic));
    const uint64_t* slice_offsets = alea_raycast_batch_ray_offsets(slice);
    const uint64_t* generic_offsets = alea_raycast_batch_ray_offsets(generic);
    const double* slice_enter = alea_raycast_batch_t_enter(slice);
    const double* slice_exit = alea_raycast_batch_t_exit(slice);
    const double* generic_enter = alea_raycast_batch_t_enter(generic);
    const double* generic_exit = alea_raycast_batch_t_exit(generic);
    const int32_t* slice_cells = alea_raycast_batch_cell_ids(slice);
    const int32_t* generic_cells = alea_raycast_batch_cell_ids(generic);
    for (size_t row = 0; row < rows; row++) ASSERT_EQ(slice_offsets[row], generic_offsets[row]);
    for (size_t i = 0; i < alea_raycast_batch_segment_count(slice); i++) {
        ASSERT_NEAR(slice_enter[i], generic_enter[i] + view.u_min, EPS);
        ASSERT_NEAR(slice_exit[i], generic_exit[i] + view.u_min, EPS);
        ASSERT_EQ(slice_cells[i], generic_cells[i]);
    }
    alea_raycast_batch_result_destroy(generic);
    alea_raycast_batch_result_destroy(slice);
    alea_destroy(sys);
}

TEST(compact_batch_internal_segments_support_per_ray_ranges) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);
    int surface = alea_sphere_surface(sys, 1, 0.0, 0.0, 0.0, 2.0);
    int material = alea_add_material(sys, 1);
    ASSERT(surface >= 0 && material >= 0);
    ASSERT(alea_add_cell(sys, 10, alea_surface_at(sys, surface)->neg_node,
                         material, -1.0, 0) >= 0);
    const double origins[] = {-5, 0, 0, -5, 0, 0};
    const double directions[] = {1, 0, 0, 1, 0, 0};
    const double t_mins[] = {4.0, 0.0};
    const double t_maxs[] = {7.0, 4.0};
    const alea_ray_batch_query_t query = {
        .kind = ALEA_RAY_QUERY_SEGMENTS,
        .t_mins = t_mins,
        .t_maxs = t_maxs
    };
    const alea_raycast_batch_options_t options = {
        .struct_size = sizeof(options),
        .fields = ALEA_RAY_BATCH_MATERIAL | ALEA_RAY_BATCH_SURFACES
    };
    alea_raycast_batch_result_t* batch = alea_raycast_batch_result_create();
    ASSERT_NOT_NULL(batch);
    ASSERT_EQ(alea_raycast_hier_batch_query_nocache(
                  sys, origins, directions, 2, &query, &options, batch), 0);
    const uint64_t* offsets = alea_raycast_batch_ray_offsets(batch);
    const double* enter = alea_raycast_batch_t_enter(batch);
    const double* exit = alea_raycast_batch_t_exit(batch);
    const int32_t* cells = alea_raycast_batch_cell_ids(batch);
    const int32_t* enter_surfaces = alea_raycast_batch_enter_surface_ids(batch);
    ASSERT_EQ(offsets[0], 0);
    ASSERT_EQ(offsets[1], 1);
    ASSERT_EQ(offsets[2], 3);
    ASSERT_NEAR(enter[0], 4.0, EPS);
    ASSERT_NEAR(exit[0], 7.0, EPS);
    ASSERT_EQ(cells[0], 10);
    ASSERT_EQ(enter_surfaces[0], -1);
    ASSERT_NEAR(enter[1], 0.0, EPS);
    ASSERT_NEAR(exit[1], 3.0, EPS);
    ASSERT_EQ(cells[1], -1);
    ASSERT_NEAR(enter[2], 3.0, EPS);
    ASSERT_NEAR(exit[2], 4.0, EPS);
    ASSERT_EQ(cells[2], 10);
    alea_raycast_batch_result_destroy(batch);
    alea_destroy(sys);
}

TEST(compact_batch_internal_first_visible_preserves_order_and_ranges) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);
    int surface = alea_sphere_surface(sys, 1, 0.0, 0.0, 0.0, 2.0);
    int material = alea_add_material(sys, 1);
    ASSERT(surface >= 0 && material >= 0);
    ASSERT(alea_add_cell(sys, 10, alea_surface_at(sys, surface)->neg_node,
                         material, -1.0, 0) >= 0);
    const double origins[] = {
        -5, 0, 0, -5, 0, 0, -5, 5, 0, -5, 0, 0, -5, 5, 0,
        -5, 0, 0, -5, 0, 0, -5, 5, 0, -5, 0, 0
    };
    const double directions[] = {
        1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0,
        1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0
    };
    const double t_mins[] = {0.0, 4.0, 0.0, 0.0, 0.0, 0.0, 4.0, 0.0, 0.0};
    const double t_maxs[] = {8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0};
    const alea_ray_batch_query_t query = {
        .kind = ALEA_RAY_QUERY_FIRST_VISIBLE,
        .fields = ALEA_RAY_QUERY_FIELD_DENSITY |
                  ALEA_RAY_QUERY_FIELD_SURFACE_ID |
                  ALEA_RAY_QUERY_FIELD_SURFACE_NORMAL |
                  ALEA_RAY_QUERY_FIELD_RESOLUTION_FLAGS,
        .material_filter = -1,
        .t_mins = t_mins,
        .t_maxs = t_maxs
    };
    alea_ray_first_visible_batch_result_t batch;
    alea_ray_first_visible_batch_result_init(&batch);
    ASSERT_EQ(alea_raycast_hier_first_visible_batch_nocache(
                  sys, origins, directions, 5, &query, &batch), 0);
    ASSERT_EQ(batch.ray_count, 5);
    ASSERT_EQ(batch.found[0], 1);
    ASSERT_EQ(batch.found[1], 1);
    ASSERT_EQ(batch.found[2], 0);
    ASSERT_NEAR(batch.t[0], 3.0, EPS);
    ASSERT_NEAR(batch.t[1], 4.0, EPS);
    ASSERT_EQ(batch.cell_ids[0], 10);
    ASSERT_EQ(batch.cell_ids[1], 10);
    ASSERT_EQ(batch.cell_ids[2], -1);
    ASSERT_EQ(batch.found[3], 1);
    ASSERT_EQ(batch.found[4], 0);
    ASSERT_NEAR(batch.t[3], 3.0, EPS);
    ASSERT_EQ(batch.cell_ids[3], 10);
    ASSERT_EQ(batch.cell_ids[4], -1);
    ASSERT_EQ(batch.surface_ids[0], 1);
    ASSERT_EQ(batch.surface_ids[1], -1);
    ASSERT_NEAR(fabs(batch.normals_xyz[0]), 1.0, EPS);
    ASSERT_NEAR(batch.normals_xyz[3], 0.0, EPS);
    const size_t packet_counts[] = {0, 1, 3, 4, 5, 9};
    const uint8_t expected_found[] = {1, 1, 0, 1, 0, 1, 1, 0, 1};
    for (size_t count_index = 0;
         count_index < sizeof(packet_counts) / sizeof(packet_counts[0]);
         count_index++) {
        const size_t count = packet_counts[count_index];
        ASSERT_EQ(alea_raycast_hier_first_visible_batch_nocache(
                      sys, origins, directions, count, &query, &batch), 0);
        ASSERT_EQ(batch.ray_count, count);
        for (size_t i = 0; i < count; i++)
            ASSERT_EQ(batch.found[i], expected_found[i]);
    }
    alea_ray_batch_query_t limited = query;
    limited.max_output_bytes = 1;
    ASSERT_EQ(alea_raycast_hier_first_visible_batch_nocache(
                  sys, origins, directions, 9, &limited, &batch), -1);
    ASSERT_EQ(batch.ray_count, 9);
    ASSERT_EQ(batch.found[0], 1);
    alea_interrupt();
    ASSERT_EQ(alea_raycast_hier_first_visible_batch_nocache(
                  sys, origins, directions, 9, &query, &batch), -1);
    ASSERT_EQ(batch.ray_count, 9);
    ASSERT_EQ(batch.found[0], 1);
    alea_clear_interrupt();
    alea_ray_first_visible_batch_result_free(&batch);
    alea_destroy(sys);
}

TEST(compact_batch_internal_any_hit_preserves_order_and_filters) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);
    int surface = alea_sphere_surface(sys, 1, 0.0, 0.0, 0.0, 2.0);
    int material = alea_add_material(sys, 1);
    ASSERT(surface >= 0 && material >= 0);
    ASSERT(alea_add_cell(sys, 10, alea_surface_at(sys, surface)->neg_node,
                         material, -1.0, 0) >= 0);
    const double origins[] = {
        -5, 0, 0, -5, 0, 0, -5, 5, 0, -5, 0, 0, -5, 5, 0,
        -5, 0, 0, -5, 0, 0, -5, 5, 0, -5, 0, 0
    };
    const double directions[] = {
        1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0,
        1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0
    };
    const double t_mins[] = {0.0, 4.0, 0.0, 0.0, 0.0, 0.0, 4.0, 0.0, 0.0};
    const double t_maxs[] = {8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0};
    const alea_ray_batch_query_t query = {
        .kind = ALEA_RAY_QUERY_ANY_HIT,
        .material_filter = -1,
        .t_mins = t_mins,
        .t_maxs = t_maxs
    };
    alea_ray_any_hit_batch_result_t batch;
    alea_ray_any_hit_batch_result_init(&batch);
    ASSERT_EQ(alea_raycast_hier_any_hit_batch_nocache(
                  sys, origins, directions, 5, &query, &batch), 0);
    ASSERT_EQ(batch.ray_count, 5);
    ASSERT_EQ(batch.hits[0], 1);
    ASSERT_EQ(batch.hits[1], 1);
    ASSERT_EQ(batch.hits[2], 0);
    ASSERT_EQ(batch.hits[3], 1);
    ASSERT_EQ(batch.hits[4], 0);
    alea_ray_batch_query_t filtered = query;
    filtered.material_filter = 99;
    ASSERT_EQ(alea_raycast_hier_any_hit_batch_nocache(
                  sys, origins, directions, 5, &filtered, &batch), 0);
    ASSERT_EQ(batch.hits[0], 0);
    ASSERT_EQ(batch.hits[1], 0);
    ASSERT_EQ(batch.hits[2], 0);
    ASSERT_EQ(batch.hits[3], 0);
    ASSERT_EQ(batch.hits[4], 0);
    const size_t packet_counts[] = {0, 1, 3, 4, 5, 9};
    const uint8_t expected_hits[] = {1, 1, 0, 1, 0, 1, 1, 0, 1};
    for (size_t count_index = 0;
         count_index < sizeof(packet_counts) / sizeof(packet_counts[0]);
         count_index++) {
        const size_t count = packet_counts[count_index];
        ASSERT_EQ(alea_raycast_hier_any_hit_batch_nocache(
                      sys, origins, directions, count, &query, &batch), 0);
        ASSERT_EQ(batch.ray_count, count);
        for (size_t i = 0; i < count; i++)
            ASSERT_EQ(batch.hits[i], expected_hits[i]);
    }
    alea_interrupt();
    ASSERT_EQ(alea_raycast_hier_any_hit_batch_nocache(
                  sys, origins, directions, 9, &query, &batch), -1);
    ASSERT_EQ(batch.ray_count, 9);
    ASSERT_EQ(batch.hits[0], 1);
    alea_clear_interrupt();
    alea_ray_any_hit_batch_result_free(&batch);
    alea_destroy(sys);
}

TEST(compact_batch_internal_boundary_events_preserve_canonical_contract) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);
    int surface = alea_sphere_surface(sys, 1, 0.0, 0.0, 0.0, 2.0);
    int material = alea_add_material(sys, 1);
    ASSERT(surface >= 0 && material >= 0);
    ASSERT(alea_add_cell(sys, 10, alea_surface_at(sys, surface)->neg_node,
                         material, -1.0, 0) >= 0);
    const double origins[] = {-5, 0, 0, -5, 0, 0, -5, 5, 0};
    const double directions[] = {1, 0, 0, 1, 0, 0, 1, 0, 0};
    const double t_mins[] = {0.0, 4.0, 0.0};
    const double t_maxs[] = {8.0, 8.0, 8.0};
    const alea_ray_batch_query_t query = {
        .kind = ALEA_RAY_QUERY_BOUNDARY_EVENTS,
        .fields = ALEA_RAY_QUERY_FIELD_PRIMITIVE_ID |
                  ALEA_RAY_QUERY_FIELD_SURFACE_NORMAL,
        .material_filter = -1,
        .t_mins = t_mins,
        .t_maxs = t_maxs
    };
    alea_ray_boundary_event_batch_result_t batch;
    alea_ray_boundary_event_batch_result_init(&batch);
    ASSERT_EQ(alea_raycast_boundary_events_batch_nocache(
                  sys, origins, directions, 3, &query, &batch), 0);
    ASSERT_EQ(batch.ray_count, 3);
    ASSERT_EQ(batch.event_count, 3);
    ASSERT_EQ(batch.ray_offsets[0], 0);
    ASSERT_EQ(batch.ray_offsets[1], 2);
    ASSERT_EQ(batch.ray_offsets[2], 3);
    ASSERT_EQ(batch.ray_offsets[3], 3);
    ASSERT_NEAR(batch.t[0], 3.0, EPS);
    ASSERT_NEAR(batch.t[1], 7.0, EPS);
    ASSERT_NEAR(batch.t[2], 7.0, EPS);
    ASSERT_EQ(batch.kinds[0], ALEA_RAY_BOUNDARY_EVENT_PHYSICAL);
    ASSERT_EQ(batch.surface_ids[0], 1);
    ASSERT_EQ(batch.cell_before[0], -1);
    ASSERT_EQ(batch.cell_after[0], 10);
    ASSERT_EQ(batch.cell_before[1], 10);
    ASSERT_EQ(batch.cell_after[1], -1);
    ASSERT_NEAR(fabs(batch.normals_xyz[0]), 1.0, EPS);
    ASSERT_NE(batch.primitive_ids[0], UINT32_MAX);
    alea_ray_boundary_event_batch_result_free(&batch);
    alea_destroy(sys);
}

TEST(selected_boundary_events_collect_complete_local_surface_group) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);
    int surface = alea_sphere_surface(sys, 1, 0.0, 0.0, 0.0, 2.0);
    int material = alea_add_material(sys, 1);
    ASSERT(surface >= 0 && material >= 0);
    ASSERT(alea_add_cell(sys, 10, alea_surface_at(sys, surface)->neg_node,
                         material, -1.0, 0) >= 0);
    ASSERT_EQ(alea_prepare_query_acceleration(sys), 0);
    alea_ray_t ray;
    ASSERT_EQ(alea_ray_init(&ray, -5.0, 0.0, 0.0, 1.0, 0.0, 0.0), 0);
    alea_raycast_result_t scratch;
    alea_ray_boundary_event_result_t events;
    alea_raycast_result_init(&scratch);
    alea_ray_boundary_event_result_init(&events);
    ASSERT_EQ(alea_raycast_selected_boundary_events_nocache(
                  sys, &ray, 8.0, &scratch, &events), 0);
    ASSERT_EQ(events.events.count, 2);
    ASSERT_EQ(events.events.data[0].kind, ALEA_RAY_BOUNDARY_EVENT_PHYSICAL);
    ASSERT_EQ(events.events.data[0].local_surface_complete, 0);
    ASSERT_EQ(events.events.data[1].kind, ALEA_RAY_BOUNDARY_EVENT_PHYSICAL);
    ASSERT_EQ(events.events.data[1].local_surface_complete, 1);
    ASSERT_EQ(events.events.data[1].local_surface_count, 1);
    ASSERT_EQ(events.events.data[1].local_surface_ids[0], 1);
    ASSERT_EQ(alea_raycast_selected_boundary_events_bidirectional_nocache(
                  sys, &ray, 8.0, &scratch, &scratch, &events), 0);
    ASSERT_EQ(events.events.count, 2);
    for (size_t i = 0; i < events.events.count; i++) {
        ASSERT_EQ(events.events.data[i].local_surface_complete, 1);
        ASSERT_EQ(events.events.data[i].local_surface_count, 1);
        ASSERT_EQ(events.events.data[i].local_surface_ids[0], 1);
    }
    alea_ray_boundary_event_result_free(&events);
    alea_raycast_result_free(&scratch);
    alea_destroy(sys);
}

TEST(compact_batch_boundary_events_preserve_coincident_physical_surfaces) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);
    int first = alea_sphere_surface(sys, 1, 0.0, 0.0, 0.0, 2.0);
    int second = alea_sphere_surface(sys, 2, 0.0, 0.0, 0.0, 2.0);
    int first_material = alea_add_material(sys, 1);
    int second_material = alea_add_material(sys, 2);
    ASSERT(first >= 0 && second >= 0);
    ASSERT(first_material >= 0 && second_material >= 0);
    ASSERT(alea_add_cell(sys, 10, alea_surface_at(sys, first)->neg_node,
                         first_material, -1.0, 0) >= 0);
    ASSERT(alea_add_cell(sys, 20, alea_surface_at(sys, second)->neg_node,
                         second_material, -1.0, 0) >= 0);
    ASSERT_EQ(alea_prepare_query_acceleration(sys), 0);

    alea_ray_t ray;
    ASSERT_EQ(alea_ray_init(&ray, -5.0, 0.0, 0.0, 1.0, 0.0, 0.0), 0);
    alea_raycast_result_t scalar_trace;
    alea_ray_boundary_event_result_t scalar_events;
    alea_raycast_result_init(&scalar_trace);
    alea_ray_boundary_event_result_init(&scalar_events);
    const alea_ray_boundary_event_options_internal_t options = {
        .include_all_coincident_physical = true
    };
    ASSERT_EQ(alea_raycast_boundary_events_with_options(
                  sys, &ray, 8.0, &options, &scalar_trace, &scalar_events), 0);
    ASSERT(scalar_events.events.count >= 4);

    const double origins[] = {-5.0, 0.0, 0.0};
    const double directions[] = {1.0, 0.0, 0.0};
    const double t_maxs[] = {8.0};
    const alea_ray_batch_query_t query = {
        .kind = ALEA_RAY_QUERY_BOUNDARY_EVENTS,
        .material_filter = -1,
        .t_maxs = t_maxs,
        .include_all_coincident_physical = true
    };
    alea_ray_boundary_event_batch_result_t batch;
    alea_ray_boundary_event_batch_result_init(&batch);
    ASSERT_EQ(alea_raycast_boundary_events_batch_nocache(
                  sys, origins, directions, 1, &query, &batch), 0);
    ASSERT_EQ(batch.event_count, scalar_events.events.count);
    for (size_t i = 0; i < scalar_events.events.count; i++) {
        ASSERT_NEAR(batch.t[i], scalar_events.events.data[i].t, EPS);
        ASSERT_EQ(batch.kinds[i], scalar_events.events.data[i].kind);
        ASSERT_EQ(batch.surface_ids[i], scalar_events.events.data[i].surface_id);
    }
    alea_ray_boundary_event_batch_result_free(&batch);
    alea_ray_boundary_event_batch_result_init(&batch);
    alea_ray_batch_query_t blas_query = query;
    blas_query.use_hier_blas = true;
    ASSERT_EQ(alea_raycast_boundary_events_batch_nocache(
                  sys, origins, directions, 1, &blas_query, &batch), 0);
    ASSERT_EQ(batch.event_count, scalar_events.events.count);
    for (size_t i = 0; i < scalar_events.events.count; i++) {
        ASSERT_NEAR(batch.t[i], scalar_events.events.data[i].t, EPS);
        ASSERT_EQ(batch.kinds[i], scalar_events.events.data[i].kind);
        ASSERT_EQ(batch.surface_ids[i], scalar_events.events.data[i].surface_id);
    }

    alea_ray_t selected_ray;
    ASSERT_EQ(alea_ray_init(&selected_ray, -5.0, 0.0, 0.0, 1.0, 0.0, 0.0), 0);
    alea_raycast_result_t selected_scratch;
    alea_ray_boundary_event_result_t selected_events;
    alea_raycast_result_init(&selected_scratch);
    alea_ray_boundary_event_result_init(&selected_events);
    ASSERT_EQ(alea_raycast_selected_boundary_events_bidirectional_nocache(
                  sys, &selected_ray, 8.0, &selected_scratch,
                  &selected_scratch, &selected_events), 1);
    alea_ray_boundary_event_result_free(&selected_events);
    alea_raycast_result_free(&selected_scratch);

    alea_ray_boundary_event_batch_result_free(&batch);
    alea_ray_boundary_event_result_free(&scalar_events);
    alea_raycast_result_free(&scalar_trace);
    alea_destroy(sys);
}

TEST(compact_ray_slice_fast_bidirectional_reuses_plot_cache) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);
    int surface = alea_sphere_surface(sys, 1, 0.0, 0.0, 0.0, 3.0);
    int material = alea_add_material(sys, 1);
    ASSERT(surface >= 0 && material >= 0);
    ASSERT(alea_add_cell(sys, 10, alea_surface_at(sys, surface)->neg_node,
                         material, -1.0, 0) >= 0);

    alea_slice_view_t view;
    alea_slice_view_init(&view, 0.0, 0.0, 0.0,
                         0.0, 0.0, 1.0,
                         0.0, 1.0, 0.0,
                         -5.0, 5.0, -2.0, 2.0);
    alea_raycast_batch_result_t* plot = alea_raycast_batch_result_create();
    alea_ray_slice_validation_result_t* validation =
        alea_ray_slice_validation_result_create();
    ASSERT_NOT_NULL(plot);
    ASSERT_NOT_NULL(validation);
    ASSERT_EQ(alea_trace_ray_slice_compact(sys, &view, 3, NULL, plot), 0);

    alea_ray_slice_validation_options_t options;
    alea_ray_slice_validation_options_init(&options);
    options.checks = ALEA_RAY_SLICE_VALIDATE_FAST_BIDIRECTIONAL;
    ASSERT_EQ(alea_validate_ray_slice_compact(sys, &view, 3, &options, NULL,
                                               plot, validation), 0);
    ASSERT_EQ(alea_ray_slice_validation_reused_trace_mask(validation),
              ALEA_RAY_SLICE_TRACE_FAST_FORWARD);
    ASSERT_EQ(alea_ray_slice_validation_executed_trace_mask(validation),
              ALEA_RAY_SLICE_TRACE_FAST_REVERSE);
    ASSERT_EQ(alea_ray_slice_validation_interval_count(validation), 0);

    ASSERT_EQ(alea_validate_ray_slice_compact(sys, &view, 3, &options, NULL,
                                               NULL, validation), 0);
    ASSERT_EQ(alea_ray_slice_validation_reused_trace_mask(validation), 0);
    ASSERT_EQ(alea_ray_slice_validation_executed_trace_mask(validation),
              ALEA_RAY_SLICE_TRACE_FAST_FORWARD |
              ALEA_RAY_SLICE_TRACE_FAST_REVERSE);
    ASSERT_EQ(alea_ray_slice_validation_interval_count(validation), 0);

    /* Any geometry change makes the plotted trace stale: validation must
     * rebuild forward rather than silently comparing a previous geometry. */
    ASSERT(alea_sphere_surface(sys, 99, 20.0, 0.0, 0.0, 1.0) >= 0);
    ASSERT_EQ(alea_validate_ray_slice_compact(sys, &view, 3, &options, NULL,
                                               plot, validation), 0);
    ASSERT_EQ(alea_ray_slice_validation_reused_trace_mask(validation), 0);
    ASSERT_EQ(alea_ray_slice_validation_executed_trace_mask(validation),
              ALEA_RAY_SLICE_TRACE_FAST_FORWARD |
              ALEA_RAY_SLICE_TRACE_FAST_REVERSE);

    alea_ray_slice_validation_result_destroy(validation);
    alea_raycast_batch_result_destroy(plot);
    alea_destroy(sys);
}

/* Two concentric cells are a total overlap that no selected-owner trace can
 * see: forward and reverse both pick the same claimant and agree.  Only the
 * complete-coverage sweep reports it. */
TEST(public_compact_coverage_slice_publishes_borrowed_csr) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);
    const int surface = alea_sphere_surface(sys, 1, 0, 0, 0, 1.0);
    const int material = alea_add_material(sys, 1);
    ASSERT(surface >= 0 && material >= 0);
    ASSERT(alea_add_cell(sys, 1, alea_surface_at(sys, surface)->neg_node,
                         material, 1.0, 0) >= 0);
    const double origins[] = {-2, 0, 0, -2, 2, 0};
    const double directions[] = {1, 0, 0, 1, 0, 0};
    const uint8_t tags[] = {4, 4};
    const double coordinates[] = {0, 2};
    alea_ray_coverage_slice_options_t options;
    alea_ray_coverage_slice_options_init(&options);
    options.t_max = 4.0;
    options.flags = ALEA_RAY_COVERAGE_DOMAIN;
    options.domain_t_min = 0.5;
    options.domain_t_max = 3.5;
    alea_ray_coverage_slice_result_t* result =
        alea_ray_coverage_slice_result_create();
    ASSERT_NOT_NULL(result);
    ASSERT_EQ(alea_ray_coverage_slice_query(
                  sys, origins, directions, 2, tags, coordinates, &options,
                  result), 0);
    ASSERT_EQ(alea_ray_coverage_slice_row_count(result), (size_t)2);
    ASSERT_EQ(alea_ray_coverage_slice_interval_count(result), (size_t)4);
    ASSERT_EQ(alea_ray_coverage_slice_owner_count(result), (size_t)1);
    const size_t* offsets = alea_ray_coverage_slice_row_offsets(result);
    const uint8_t* kinds = alea_ray_coverage_slice_kinds(result);
    const int* owner_cells = alea_ray_coverage_slice_owner_cell_ids(result);
    ASSERT_NOT_NULL(offsets);
    ASSERT_NOT_NULL(kinds);
    ASSERT_NOT_NULL(owner_cells);
    ASSERT_EQ(offsets[0], (size_t)0);
    ASSERT_EQ(offsets[1], (size_t)3);
    ASSERT_EQ(offsets[2], (size_t)4);
    ASSERT_EQ(kinds[1], (uint8_t)ALEA_RAY_COVERAGE_UNIQUE);
    ASSERT_EQ(owner_cells[0], 1);
    alea_ray_coverage_slice_result_t* scalar =
        alea_ray_coverage_slice_result_create();
    ASSERT_NOT_NULL(scalar);
    ASSERT_EQ(alea_ray_coverage_query(sys, -2, 0, 0, 1, 0, 0,
                                      &options, scalar), 0);
    ASSERT_EQ(alea_ray_coverage_slice_row_count(scalar), (size_t)1);
    ASSERT_EQ(alea_ray_coverage_slice_interval_count(scalar), (size_t)3);
    alea_ray_coverage_slice_result_destroy(scalar);
    alea_ray_coverage_slice_result_t* exterior =
        alea_ray_coverage_slice_result_create();
    ASSERT_NOT_NULL(exterior);
    options.flags = ALEA_RAY_COVERAGE_DOMAIN |
                    ALEA_RAY_COVERAGE_REPORT_EXTERIOR;
    options.domain_t_min = 1.0;
    options.domain_t_max = 3.0;
    ASSERT_EQ(alea_ray_coverage_query(sys, -2, 0, 0, 1, 0, 0,
                                      &options, exterior), 0);
    ASSERT_EQ(alea_ray_coverage_slice_interval_count(exterior), (size_t)3);
    const uint8_t* exterior_kinds = alea_ray_coverage_slice_kinds(exterior);
    ASSERT_EQ(exterior_kinds[0], (uint8_t)ALEA_RAY_COVERAGE_ALLOWED_EXTERIOR);
    ASSERT_EQ(exterior_kinds[1], (uint8_t)ALEA_RAY_COVERAGE_UNIQUE);
    ASSERT_EQ(exterior_kinds[2], (uint8_t)ALEA_RAY_COVERAGE_ALLOWED_EXTERIOR);
    alea_ray_coverage_slice_result_destroy(exterior);
    options.flags = ALEA_RAY_COVERAGE_DOMAIN;
    options.domain_t_min = 0.5;
    options.domain_t_max = 3.5;
    alea_ray_coverage_slice_result_t* adaptive =
        alea_ray_coverage_slice_result_create();
    ASSERT_NOT_NULL(adaptive);
    options.max_refinement_depth = 1;
    ASSERT_EQ(alea_ray_coverage_slice_query(
                  sys, origins, directions, 2, tags, coordinates, &options,
                  adaptive), 0);
    ASSERT_EQ(alea_ray_coverage_slice_row_count(adaptive), (size_t)3);
    ASSERT_EQ(alea_ray_coverage_slice_refinement_status(adaptive),
              ALEA_RAY_COVERAGE_REFINEMENT_MAX_DEPTH);
    options.max_refinement_depth = 0;
    alea_ray_coverage_slice_result_destroy(adaptive);
    const size_t* previous_offsets = offsets;
    options.max_output_bytes = 1;
    ASSERT_EQ(alea_ray_coverage_slice_query(
                  sys, origins, directions, 2, tags, coordinates, &options,
                  result), -1);
    ASSERT_EQ(alea_ray_coverage_slice_row_offsets(result), previous_offsets);
    ASSERT_EQ(alea_ray_coverage_slice_interval_count(result), (size_t)4);
    alea_ray_coverage_slice_result_destroy(result);
    alea_destroy(sys);
}

TEST(compact_ray_slice_coverage_reports_overlap_and_gaps) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);
    const int outer = alea_sphere_surface(sys, 1, 0.0, 0.0, 0.0, 3.0);
    const int inner = alea_sphere_surface(sys, 2, 0.0, 0.0, 0.0, 1.0);
    const int material = alea_add_material(sys, 1);
    ASSERT(outer >= 0 && inner >= 0 && material >= 0);
    ASSERT(alea_add_cell(sys, 10, alea_surface_at(sys, outer)->neg_node,
                         material, -1.0, 0) >= 0);
    ASSERT(alea_add_cell(sys, 20, alea_surface_at(sys, inner)->neg_node,
                         material, -1.0, 0) >= 0);

    alea_slice_view_t view;
    alea_slice_view_init(&view, 0.0, 0.0, 0.0,
                         0.0, 0.0, 1.0,
                         0.0, 1.0, 0.0,
                         -5.0, 5.0, -2.0, 2.0);
    alea_ray_slice_validation_result_t* validation =
        alea_ray_slice_validation_result_create();
    ASSERT_NOT_NULL(validation);

    alea_ray_slice_validation_options_t options;
    alea_ray_slice_validation_options_init(&options);
    options.checks = ALEA_RAY_SLICE_VALIDATE_COVERAGE;

    /* Neither an explicit domain nor a uniform policy: lowering must reject
     * the request rather than guess where ownership is required. */
    ASSERT_EQ(alea_validate_ray_slice_compact(sys, &view, 3, &options, NULL,
                                              NULL, validation), -1);
    options.coverage_flags = ALEA_RAY_SLICE_COVERAGE_HAS_DOMAIN |
                             ALEA_RAY_SLICE_COVERAGE_UNOWNED_IS_EXTERIOR;
    ASSERT_EQ(alea_validate_ray_slice_compact(sys, &view, 3, &options, NULL,
                                              NULL, validation), -1);
    options.coverage_flags = ALEA_RAY_SLICE_COVERAGE_HAS_DOMAIN;
    options.coverage_domain_u_min = 2.0;
    options.coverage_domain_u_max = 2.0;
    ASSERT_EQ(alea_validate_ray_slice_compact(sys, &view, 3, &options, NULL,
                                              NULL, validation), -1);

    /* An explicit domain spanning the view: unowned space inside it is an
     * interior gap. */
    options.coverage_domain_u_min = -5.0;
    options.coverage_domain_u_max = 5.0;
    ASSERT_EQ(alea_validate_ray_slice_compact(sys, &view, 3, &options, NULL,
                                              NULL, validation), 0);
    ASSERT(alea_ray_slice_validation_fields(validation) &
           ALEA_RAY_SLICE_VALIDATION_FIELD_COVERAGE);
    const uint64_t* offsets = alea_ray_slice_validation_row_offsets(validation);
    const double* u_enter = alea_ray_slice_validation_u_enter(validation);
    const double* u_exit = alea_ray_slice_validation_u_exit(validation);
    const uint32_t* flags =
        alea_ray_slice_validation_diagnostic_flags(validation);
    const uint32_t* owners =
        alea_ray_slice_validation_coverage_owner_counts(validation);
    ASSERT_NOT_NULL(offsets);
    ASSERT_NOT_NULL(owners);
    ASSERT_EQ(alea_ray_slice_validation_row_count(validation), (size_t)3);

    /* The centre row crosses both spheres: gap, overlap, gap.  The uniquely
     * owned shell intervals are the expected answer and are not reported. */
    const size_t begin = (size_t)offsets[1], end = (size_t)offsets[2];
    ASSERT_EQ(end - begin, (size_t)3);
    ASSERT_EQ(flags[begin], ALEA_RAY_SLICE_DIAG_COVERAGE_GAP);
    ASSERT_NEAR(u_enter[begin], -5.0, 1e-9);
    ASSERT_NEAR(u_exit[begin], -3.0, 1e-9);
    ASSERT_EQ(owners[begin], (uint32_t)0);
    ASSERT_EQ(flags[begin + 1], ALEA_RAY_SLICE_DIAG_COVERAGE_OVERLAP);
    ASSERT_NEAR(u_enter[begin + 1], -1.0, 1e-9);
    ASSERT_NEAR(u_exit[begin + 1], 1.0, 1e-9);
    ASSERT_EQ(owners[begin + 1], (uint32_t)2);
    ASSERT_EQ(flags[begin + 2], ALEA_RAY_SLICE_DIAG_COVERAGE_GAP);
    ASSERT_NEAR(u_enter[begin + 2], 3.0, 1e-9);
    ASSERT_NEAR(u_exit[begin + 2], 5.0, 1e-9);

    /* An off-centre row misses the inner sphere and reports gaps only. */
    for (size_t i = (size_t)offsets[0]; i < (size_t)offsets[1]; i++)
        ASSERT_EQ(flags[i], ALEA_RAY_SLICE_DIAG_COVERAGE_GAP);

    /* The uniform policy declares all unowned space exterior, so the same
     * geometry yields the overlap alone.  Nothing about the viewport changed. */
    options.coverage_flags = ALEA_RAY_SLICE_COVERAGE_UNOWNED_IS_EXTERIOR;
    ASSERT_EQ(alea_validate_ray_slice_compact(sys, &view, 3, &options, NULL,
                                              NULL, validation), 0);
    offsets = alea_ray_slice_validation_row_offsets(validation);
    flags = alea_ray_slice_validation_diagnostic_flags(validation);
    ASSERT_EQ(alea_ray_slice_validation_interval_count(validation), (size_t)1);
    ASSERT_EQ(offsets[1], (uint64_t)0);
    ASSERT_EQ(offsets[2], (uint64_t)1);
    ASSERT_EQ(flags[0], ALEA_RAY_SLICE_DIAG_COVERAGE_OVERLAP);

    options.coverage_flags = ALEA_RAY_SLICE_COVERAGE_UNOWNED_IS_EXTERIOR |
                             ALEA_RAY_SLICE_COVERAGE_REPORT_EXTERIOR;
    ASSERT_EQ(alea_validate_ray_slice_compact(sys, &view, 3, &options, NULL,
                                              NULL, validation), 0);
    offsets = alea_ray_slice_validation_row_offsets(validation);
    flags = alea_ray_slice_validation_diagnostic_flags(validation);
    ASSERT_EQ((size_t)offsets[2] - (size_t)offsets[1], (size_t)3);
    ASSERT_EQ(flags[(size_t)offsets[1]],
              ALEA_RAY_SLICE_DIAG_COVERAGE_ALLOWED_EXTERIOR);
    ASSERT_EQ(flags[(size_t)offsets[1] + 1],
              ALEA_RAY_SLICE_DIAG_COVERAGE_OVERLAP);

    /* Directional agreement is not evidence of validity: the same rows the
     * coverage sweep condemns produce no bidirectional mismatch at all. */
    options.checks = ALEA_RAY_SLICE_VALIDATE_FAST_BIDIRECTIONAL;
    options.coverage_flags = 0;
    ASSERT_EQ(alea_validate_ray_slice_compact(sys, &view, 3, &options, NULL,
                                              NULL, validation), 0);
    ASSERT_EQ(alea_ray_slice_validation_interval_count(validation), (size_t)0);
    ASSERT_EQ(alea_ray_slice_validation_fields(validation) &
              ALEA_RAY_SLICE_VALIDATION_FIELD_COVERAGE, (uint32_t)0);

    /* Both checks together keep the two evidence classes in separate flags. */
    options.checks = ALEA_RAY_SLICE_VALIDATE_FAST_BIDIRECTIONAL |
                     ALEA_RAY_SLICE_VALIDATE_COVERAGE;
    options.coverage_flags = ALEA_RAY_SLICE_COVERAGE_UNOWNED_IS_EXTERIOR;
    ASSERT_EQ(alea_validate_ray_slice_compact(sys, &view, 3, &options, NULL,
                                              NULL, validation), 0);
    flags = alea_ray_slice_validation_diagnostic_flags(validation);
    ASSERT_EQ(alea_ray_slice_validation_interval_count(validation), (size_t)1);
    ASSERT_EQ(flags[0], ALEA_RAY_SLICE_DIAG_COVERAGE_OVERLAP);

    alea_ray_slice_validation_result_destroy(validation);
    alea_destroy(sys);
}

/* The inner cell only shows up on rows close to the centre, so adjacent
 * coverage signatures differ and adaptive refinement has something to find
 * between the requested rows. */
TEST(compact_ray_slice_coverage_refines_rows_adaptively) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);
    const int outer = alea_sphere_surface(sys, 1, 0.0, 0.0, 0.0, 3.0);
    const int inner = alea_sphere_surface(sys, 2, 0.0, 0.0, 0.0, 1.0);
    const int material = alea_add_material(sys, 1);
    ASSERT(outer >= 0 && inner >= 0 && material >= 0);
    ASSERT(alea_add_cell(sys, 10, alea_surface_at(sys, outer)->neg_node,
                         material, -1.0, 0) >= 0);
    ASSERT(alea_add_cell(sys, 20, alea_surface_at(sys, inner)->neg_node,
                         material, -1.0, 0) >= 0);

    alea_slice_view_t view;
    alea_slice_view_init(&view, 0.0, 0.0, 0.0,
                         0.0, 0.0, 1.0,
                         0.0, 1.0, 0.0,
                         -5.0, 5.0, -2.0, 2.0);
    alea_ray_slice_validation_result_t* validation =
        alea_ray_slice_validation_result_create();
    ASSERT_NOT_NULL(validation);

    alea_ray_slice_validation_options_t options;
    alea_ray_slice_validation_options_init(&options);
    options.checks = ALEA_RAY_SLICE_VALIDATE_COVERAGE;
    options.coverage_flags = ALEA_RAY_SLICE_COVERAGE_UNOWNED_IS_EXTERIOR;

    /* No refinement requested: the published rows stay one-to-one with the
     * requested rows and no row provenance is published. */
    ASSERT_EQ(alea_validate_ray_slice_compact(sys, &view, 3, &options, NULL,
                                              NULL, validation), 0);
    ASSERT_EQ(alea_ray_slice_validation_row_count(validation), (size_t)3);
    ASSERT_EQ(alea_ray_slice_validation_fields(validation) &
              ALEA_RAY_SLICE_VALIDATION_FIELD_ADAPTIVE_ROWS, (uint32_t)0);
    ASSERT_EQ(alea_ray_slice_validation_refinement_status(validation),
              ALEA_RAY_SLICE_REFINEMENT_NOT_REQUESTED);
    ASSERT_NULL(alea_ray_slice_validation_row_base_indices(validation));

    /* Refinement without the coverage check has nothing to measure. */
    options.checks = ALEA_RAY_SLICE_VALIDATE_FAST_BIDIRECTIONAL;
    options.coverage_max_refinement_depth = 1;
    ASSERT_EQ(alea_validate_ray_slice_compact(sys, &view, 3, &options, NULL,
                                              NULL, validation), -1);

    /* One wave splits both differing pairs, publishing five ordered rows. */
    options.checks = ALEA_RAY_SLICE_VALIDATE_COVERAGE |
                     ALEA_RAY_SLICE_VALIDATE_FAST_BIDIRECTIONAL;
    ASSERT_EQ(alea_validate_ray_slice_compact(sys, &view, 3, &options, NULL,
                                              NULL, validation), 0);
    ASSERT_EQ(alea_ray_slice_validation_row_count(validation), (size_t)5);
    ASSERT(alea_ray_slice_validation_fields(validation) &
           ALEA_RAY_SLICE_VALIDATION_FIELD_ADAPTIVE_ROWS);
    ASSERT_EQ(alea_ray_slice_validation_refinement_status(validation),
              ALEA_RAY_SLICE_REFINEMENT_MAX_DEPTH);

    const size_t* base = alea_ray_slice_validation_row_base_indices(validation);
    const double* coordinates =
        alea_ray_slice_validation_row_transverse_coordinates(validation);
    const uint8_t* tags = alea_ray_slice_validation_row_direction_tags(validation);
    ASSERT_NOT_NULL(base);
    ASSERT_NOT_NULL(coordinates);
    ASSERT_NOT_NULL(tags);
    /* Requested rows keep their index; refined rows sit between them. */
    ASSERT_EQ(base[0], (size_t)0);
    ASSERT_EQ(base[1], SIZE_MAX);
    ASSERT_EQ(base[2], (size_t)1);
    ASSERT_EQ(base[3], SIZE_MAX);
    ASSERT_EQ(base[4], (size_t)2);
    ASSERT_NEAR(coordinates[0], -4.0 / 3.0, 1e-12);
    ASSERT_NEAR(coordinates[1], -2.0 / 3.0, 1e-12);
    ASSERT_NEAR(coordinates[2], 0.0, 1e-12);
    ASSERT_NEAR(coordinates[3], 2.0 / 3.0, 1e-12);
    ASSERT_NEAR(coordinates[4], 4.0 / 3.0, 1e-12);
    for (size_t row = 0; row + 1 < 5; row++) {
        ASSERT_EQ(tags[row], (uint8_t)0);
        ASSERT(coordinates[row] < coordinates[row + 1]);
    }

    /* A refined row still reports the overlap it was generated to find, and
     * reports no directional evidence: the fast traces cover the requested
     * viewport rows only. */
    const uint64_t* offsets = alea_ray_slice_validation_row_offsets(validation);
    const uint32_t* flags =
        alea_ray_slice_validation_diagnostic_flags(validation);
    const int32_t* forward_cells =
        alea_ray_slice_validation_fast_forward_cell_ids(validation);
    const uint32_t* owners =
        alea_ray_slice_validation_coverage_owner_counts(validation);
    ASSERT_EQ((size_t)offsets[4] - (size_t)offsets[3], (size_t)1);
    const size_t refined = (size_t)offsets[3];
    ASSERT_EQ(flags[refined], ALEA_RAY_SLICE_DIAG_COVERAGE_OVERLAP);
    ASSERT_EQ(owners[refined], (uint32_t)2);
    ASSERT_EQ(forward_cells[refined], -1);

    /* With an event cache the requested rows keep boundary provenance, while a
     * refined row publishes none rather than borrowing a neighbour's line. */
    alea_slice_directional_event_cache_t* cache =
        alea_slice_directional_event_cache_create(sys, &view, 64, 3);
    ASSERT_NOT_NULL(cache);
    options.flags = ALEA_RAY_SLICE_VALIDATION_INCLUDE_AGREEMENTS;
    ASSERT_EQ(alea_validate_ray_slice_compact_with_event_cache(
                  sys, &view, 3, &options, NULL, NULL, cache, validation), 0);
    ASSERT_EQ(alea_ray_slice_validation_row_count(validation), (size_t)5);
    offsets = alea_ray_slice_validation_row_offsets(validation);
    base = alea_ray_slice_validation_row_base_indices(validation);
    const int32_t* enter_surfaces =
        alea_ray_slice_validation_u_enter_forward_surface_ids(validation);
    ASSERT_NOT_NULL(enter_surfaces);
    int requested_surface_seen = 0;
    for (size_t row = 0; row < 5; row++) {
        for (size_t i = (size_t)offsets[row]; i < (size_t)offsets[row + 1]; i++) {
            if (base[row] == SIZE_MAX)
                ASSERT_EQ(enter_surfaces[i], -1);
            else if (enter_surfaces[i] > 0)
                requested_surface_seen = 1;
        }
    }
    ASSERT_EQ(requested_surface_seen, 1);
    alea_slice_directional_event_cache_destroy(cache);
    options.flags = 0;

    /* Row budget and spacing limits publish completed rows with an explicit
     * limited status rather than claiming the criteria converged. */
    options.checks = ALEA_RAY_SLICE_VALIDATE_COVERAGE;
    options.coverage_max_refinement_depth = 4;
    options.max_coverage_rows = 3;
    ASSERT_EQ(alea_validate_ray_slice_compact(sys, &view, 3, &options, NULL,
                                              NULL, validation), 0);
    ASSERT_EQ(alea_ray_slice_validation_row_count(validation), (size_t)3);
    ASSERT_EQ(alea_ray_slice_validation_refinement_status(validation),
              ALEA_RAY_SLICE_REFINEMENT_MAX_ROWS);

    options.max_coverage_rows = 0;
    options.coverage_min_transverse_spacing = 1.0;
    ASSERT_EQ(alea_validate_ray_slice_compact(sys, &view, 3, &options, NULL,
                                              NULL, validation), 0);
    ASSERT_EQ(alea_ray_slice_validation_row_count(validation), (size_t)3);
    ASSERT_EQ(alea_ray_slice_validation_refinement_status(validation),
              ALEA_RAY_SLICE_REFINEMENT_MIN_SPACING);

    /* A signal that needs a tolerance must not silently disable itself. */
    options.coverage_min_transverse_spacing = 0.0;
    options.coverage_refine_signals = ALEA_RAY_SLICE_REFINE_DISPLACEMENT;
    ASSERT_EQ(alea_validate_ray_slice_compact(sys, &view, 3, &options, NULL,
                                              NULL, validation), -1);
    options.coverage_refine_signals = 1u << 20;
    ASSERT_EQ(alea_validate_ray_slice_compact(sys, &view, 3, &options, NULL,
                                              NULL, validation), -1);

    alea_ray_slice_validation_result_destroy(validation);
    alea_destroy(sys);
}

TEST(transition_slice_active_boundary_filter_partitions_narrow_regions) {
    alea_slice_view_t view;
    alea_slice_view_init(&view, 0.0, 0.0, 0.0,
                         0.0, 0.0, 1.0, 0.0, 1.0, 0.0,
                         -2.0, 2.0, -1.0, 1.0);
    alea_transition_slice_critical_tile_t tile;
    memset(&tile, 0, sizeof(tile));
    tile.uv_min[0] = -2.0; tile.uv_max[0] = 2.0;
    tile.uv_min[1] = -1.0; tile.uv_max[1] = 1.0;
    alea_transition_slice_options_t options;
    alea_transition_slice_options_init(&options);

    /* The first branch is wholly contained by x < 1.  All three of its card
     * curves cross the tile but none is part of the cell boundary. */
    alea_system_t* redundant = alea_create();
    ASSERT_NOT_NULL(redundant);
    const int x0 = alea_plane_surface(redundant, 801, 1, 0, 0, 0);
    const int y0 = alea_plane_surface(redundant, 802, 0, 1, 0, 0);
    const int yd = alea_plane_surface(redundant, 803, 0, 1, 0, -1e-5);
    const int x1 = alea_plane_surface(redundant, 804, 1, 0, 0, -1);
    const int rxlo = alea_plane_surface(redundant, 805, 1, 0, 0, 3);
    const int rylo = alea_plane_surface(redundant, 806, 0, 1, 0, 2);
    const int ryhi = alea_plane_surface(redundant, 807, 0, 1, 0, -2);
    const int rzlo = alea_plane_surface(redundant, 808, 0, 0, 1, 1);
    const int rzhi = alea_plane_surface(redundant, 809, 0, 0, 1, -1);
    ASSERT(x0 >= 0 && y0 >= 0 && yd >= 0 && x1 >= 0 && rxlo >= 0 &&
           rylo >= 0 && ryhi >= 0 && rzlo >= 0 && rzhi >= 0);
    alea_node_id_t narrow = alea_intersection(
        redundant, alea_halfspace(redundant, x0, -1),
        alea_intersection(redundant,
            alea_halfspace(redundant, y0, 1),
            alea_halfspace(redundant, yd, -1)));
    alea_node_id_t redundant_union = alea_union(
        redundant, narrow, alea_halfspace(redundant, x1, -1));
    const alea_node_id_t redundant_nodes[] = {
        redundant_union,
        alea_halfspace(redundant, rxlo, 1),
        alea_halfspace(redundant, rylo, 1),
        alea_halfspace(redundant, ryhi, -1),
        alea_halfspace(redundant, rzlo, 1),
        alea_halfspace(redundant, rzhi, -1)
    };
    alea_node_id_t region = alea_intersection_n(
        redundant, redundant_nodes,
        sizeof(redundant_nodes) / sizeof(redundant_nodes[0]));
    ASSERT(alea_add_cell(redundant, 801, region,
                         ALEA_MATERIAL_VOID, 0.0, 0) >= 0);
    ASSERT_EQ(alea_prepare_query_acceleration(redundant), 0);
    alea_transition_slice_stats_t stats;
    memset(&stats, 0, sizeof(stats));
    ASSERT_EQ(alea_transition_slice_enumerate_critical_tiles(
                  redundant, &view, &options, &tile, 1, NULL, NULL, &stats), 0);
    ASSERT_EQ(stats.critical_curves, (size_t)1);
    ASSERT(stats.critical_curves_culled >= (size_t)3);
    ASSERT(stats.critical_active_boundary_tests > 0);
    ASSERT_EQ(stats.critical_active_segments, (size_t)1);
    alea_destroy(redundant);

    /* Exact arrangement side tests also prove inactivity when a card occurs
     * with both senses.  Here the union of its two halfspaces is a Boolean
     * tautology inside a finite box: retaining the whole x=0 line would add
     * conservative work without representing a cell boundary. */
    alea_system_t* both_sense_redundant = alea_create();
    ASSERT_NOT_NULL(both_sense_redundant);
    const int bsx = alea_plane_surface(
        both_sense_redundant, 8101, 1, 0, 0, 0);
    const int bsxlo = alea_plane_surface(
        both_sense_redundant, 8102, 1, 0, 0, 3);
    const int bsxhi = alea_plane_surface(
        both_sense_redundant, 8103, 1, 0, 0, -3);
    const int bsylo = alea_plane_surface(
        both_sense_redundant, 8104, 0, 1, 0, 2);
    const int bsyhi = alea_plane_surface(
        both_sense_redundant, 8105, 0, 1, 0, -2);
    const int bszlo = alea_plane_surface(
        both_sense_redundant, 8106, 0, 0, 1, 1);
    const int bszhi = alea_plane_surface(
        both_sense_redundant, 8107, 0, 0, 1, -1);
    ASSERT(bsx >= 0 && bsxlo >= 0 && bsxhi >= 0 && bsylo >= 0 &&
           bsyhi >= 0 && bszlo >= 0 && bszhi >= 0);
    const alea_node_id_t both_sense_tautology = alea_union(
        both_sense_redundant,
        alea_halfspace(both_sense_redundant, bsx, -1),
        alea_halfspace(both_sense_redundant, bsx, 1));
    const alea_node_id_t both_sense_nodes[] = {
        both_sense_tautology,
        alea_halfspace(both_sense_redundant, bsxlo, 1),
        alea_halfspace(both_sense_redundant, bsxhi, -1),
        alea_halfspace(both_sense_redundant, bsylo, 1),
        alea_halfspace(both_sense_redundant, bsyhi, -1),
        alea_halfspace(both_sense_redundant, bszlo, 1),
        alea_halfspace(both_sense_redundant, bszhi, -1)
    };
    region = alea_intersection_n(
        both_sense_redundant, both_sense_nodes,
        sizeof(both_sense_nodes) / sizeof(both_sense_nodes[0]));
    ASSERT(alea_add_cell(both_sense_redundant, 8101, region,
                         ALEA_MATERIAL_VOID, 0.0, 0) >= 0);
    ASSERT_EQ(alea_prepare_query_acceleration(both_sense_redundant), 0);
    memset(&stats, 0, sizeof(stats));
    ASSERT_EQ(alea_transition_slice_enumerate_critical_tiles(
                  both_sense_redundant, &view, &options, &tile, 1,
                  NULL, NULL, &stats), 0);
    ASSERT(stats.critical_curves_culled >= (size_t)1);
    ASSERT_EQ(stats.critical_whole_curve_fallbacks, (size_t)0);
    alea_destroy(both_sense_redundant);

    /* The same 1e-5 interval is much narrower than any fixed 3x3 seed grid,
     * but curve intersections partition x=0 at both endpoints. */
    alea_system_t* active = alea_create();
    ASSERT_NOT_NULL(active);
    const int ax = alea_plane_surface(active, 811, 1, 0, 0, 0);
    const int ay0 = alea_plane_surface(active, 812, 0, 1, 0, 0);
    const int ayd = alea_plane_surface(active, 813, 0, 1, 0, -1e-5);
    const int axlo = alea_plane_surface(active, 814, 1, 0, 0, 3);
    const int azlo = alea_plane_surface(active, 815, 0, 0, 1, 1);
    const int azhi = alea_plane_surface(active, 816, 0, 0, 1, -1);
    ASSERT(ax >= 0 && ay0 >= 0 && ayd >= 0 &&
           axlo >= 0 && azlo >= 0 && azhi >= 0);
    alea_node_id_t active_narrow = alea_intersection(
        active, alea_halfspace(active, ax, -1),
        alea_intersection(active,
            alea_halfspace(active, ay0, 1),
            alea_halfspace(active, ayd, -1)));
    const alea_node_id_t active_nodes[] = {
        active_narrow,
        alea_halfspace(active, axlo, 1),
        alea_halfspace(active, azlo, 1),
        alea_halfspace(active, azhi, -1)
    };
    region = alea_intersection_n(
        active, active_nodes,
        sizeof(active_nodes) / sizeof(active_nodes[0]));
    ASSERT(alea_add_cell(active, 811, region,
                         ALEA_MATERIAL_VOID, 0.0, 0) >= 0);
    ASSERT_EQ(alea_prepare_query_acceleration(active), 0);
    memset(&stats, 0, sizeof(stats));
    ASSERT_EQ(alea_transition_slice_enumerate_critical_tiles(
                  active, &view, &options, &tile, 1, NULL, NULL, &stats), 0);
    ASSERT_EQ(stats.critical_curves, (size_t)3);
    ASSERT_EQ(stats.critical_curves_culled, (size_t)0);
    ASSERT_EQ(stats.critical_active_segments, (size_t)3);

    /* One exact card can contribute multiple disconnected active patches.  It
     * must be charged as two segments, not collapsed back to one whole line. */
    alea_system_t* disjoint = alea_create();
    ASSERT_NOT_NULL(disjoint);
    const int dx = alea_plane_surface(disjoint, 821, 1, 0, 0, 0);
    const int dy0 = alea_plane_surface(disjoint, 822, 0, 1, 0, 0);
    const int dy1 = alea_plane_surface(disjoint, 823, 0, 1, 0, -1e-5);
    const int dy2 = alea_plane_surface(disjoint, 824, 0, 1, 0, -0.5);
    const int dy3 = alea_plane_surface(disjoint, 825, 0, 1, 0, -0.50001);
    const int dxlo = alea_plane_surface(disjoint, 826, 1, 0, 0, 3);
    const int dzlo = alea_plane_surface(disjoint, 827, 0, 0, 1, 1);
    const int dzhi = alea_plane_surface(disjoint, 828, 0, 0, 1, -1);
    ASSERT(dx >= 0 && dy0 >= 0 && dy1 >= 0 && dy2 >= 0 && dy3 >= 0 &&
           dxlo >= 0 && dzlo >= 0 && dzhi >= 0);
    const alea_node_id_t strip0 = alea_intersection(
        disjoint, alea_halfspace(disjoint, dx, -1),
        alea_intersection(disjoint,
            alea_halfspace(disjoint, dy0, 1),
            alea_halfspace(disjoint, dy1, -1)));
    const alea_node_id_t strip1 = alea_intersection(
        disjoint, alea_halfspace(disjoint, dx, -1),
        alea_intersection(disjoint,
            alea_halfspace(disjoint, dy2, 1),
            alea_halfspace(disjoint, dy3, -1)));
    const alea_node_id_t disjoint_nodes[] = {
        alea_union(disjoint, strip0, strip1),
        alea_halfspace(disjoint, dxlo, 1),
        alea_halfspace(disjoint, dzlo, 1),
        alea_halfspace(disjoint, dzhi, -1)
    };
    region = alea_intersection_n(
        disjoint, disjoint_nodes,
        sizeof(disjoint_nodes) / sizeof(disjoint_nodes[0]));
    ASSERT(alea_add_cell(disjoint, 821, region,
                         ALEA_MATERIAL_VOID, 0.0, 0) >= 0);
    ASSERT_EQ(alea_prepare_query_acceleration(disjoint), 0);
    memset(&stats, 0, sizeof(stats));
    ASSERT_EQ(alea_transition_slice_enumerate_critical_tiles(
                  disjoint, &view, &options, &tile, 1,
                  NULL, NULL, &stats), 0);
    ASSERT_EQ(stats.critical_active_segments, (size_t)6);
    ASSERT_EQ(stats.critical_curves, (size_t)6);
    alea_destroy(disjoint);

    /* A cylinder sliced parallel to its axis yields one analytical
     * PARALLEL_LINES curve.  It must become two independently clipped and
     * CSG-tested active segments. */
    alea_system_t* parallel = alea_create();
    ASSERT_NOT_NULL(parallel);
    const int cyl = alea_cylinder_z_surface(parallel, 831, 0, 0, 1);
    const int pzlo = alea_plane_surface(parallel, 832, 0, 0, 1, 2);
    const int pzhi = alea_plane_surface(parallel, 833, 0, 0, 1, -2);
    ASSERT(cyl >= 0 && pzlo >= 0 && pzhi >= 0);
    const alea_node_id_t cylinder_nodes[] = {
        alea_halfspace(parallel, cyl, -1),
        alea_halfspace(parallel, pzlo, 1),
        alea_halfspace(parallel, pzhi, -1)
    };
    region = alea_intersection_n(
        parallel, cylinder_nodes,
        sizeof(cylinder_nodes) / sizeof(cylinder_nodes[0]));
    ASSERT(alea_add_cell(parallel, 831, region,
                         ALEA_MATERIAL_VOID, 0.0, 0) >= 0);
    ASSERT_EQ(alea_prepare_query_acceleration(parallel), 0);
    alea_slice_view_t parallel_view;
    alea_slice_view_init(&parallel_view, 0.0, 0.0, 0.0,
                         1.0, 0.0, 0.0, 0.0, 0.0, 1.0,
                         -2.0, 2.0, -1.0, 1.0);
    memset(&stats, 0, sizeof(stats));
    ASSERT_EQ(alea_transition_slice_enumerate_critical_tiles(
                  parallel, &parallel_view, &options, &tile, 1,
                  NULL, NULL, &stats), 0);
    ASSERT_EQ(stats.critical_active_segments, (size_t)2);
    ASSERT_EQ(stats.critical_curves, (size_t)2);
    ASSERT_EQ(stats.critical_whole_curve_fallbacks, (size_t)0);
    alea_destroy(parallel);

    /* A full circular cell boundary is retained as one proven active arc,
     * rather than consuming a conservative whole-circle fallback. */
    alea_system_t* circular = alea_create();
    ASSERT_NOT_NULL(circular);
    const int sphere = alea_sphere_surface(circular, 841, 0, 0, 0, 0.75);
    ASSERT(sphere >= 0);
    ASSERT(alea_add_cell(circular, 841,
                         alea_halfspace(circular, sphere, -1),
                         ALEA_MATERIAL_VOID, 0.0, 0) >= 0);
    ASSERT_EQ(alea_prepare_query_acceleration(circular), 0);
    memset(&stats, 0, sizeof(stats));
    ASSERT_EQ(alea_transition_slice_enumerate_critical_tiles(
                  circular, &view, &options, &tile, 1,
                  NULL, NULL, &stats), 0);
    ASSERT_EQ(stats.critical_active_segments, (size_t)1);
    ASSERT_EQ(stats.critical_curves, (size_t)1);
    ASSERT_EQ(stats.critical_whole_curve_fallbacks, (size_t)0);
    ASSERT_EQ(stats.critical_unsupported_curves, (size_t)0);
    alea_destroy(circular);

    /* A tangent slice produces one exact point, not an unbounded source that
     * needs conservative active-boundary fallback. */
    alea_system_t* tangent = alea_create();
    ASSERT_NOT_NULL(tangent);
    const int tangent_sphere = alea_sphere_surface(
        tangent, 843, 0, 0, 0, 0.75);
    ASSERT(tangent_sphere >= 0);
    ASSERT(alea_add_cell(
        tangent, 843, alea_halfspace(tangent, tangent_sphere, -1),
        ALEA_MATERIAL_VOID, 0.0, 0) >= 0);
    ASSERT_EQ(alea_prepare_query_acceleration(tangent), 0);
    alea_slice_view_t tangent_view;
    alea_slice_view_init(&tangent_view, 0.0, 0.0, 0.75,
                         0.0, 0.0, 1.0, 0.0, 1.0, 0.0,
                         -2.0, 2.0, -1.0, 1.0);
    memset(&stats, 0, sizeof(stats));
    ASSERT_EQ(alea_transition_slice_enumerate_critical_tiles(
                  tangent, &tangent_view, &options, &tile, 1,
                  NULL, NULL, &stats), 0);
    ASSERT_EQ(stats.critical_active_unsupported_point_fallbacks,
              (size_t)0);
    ASSERT_EQ(stats.critical_whole_curve_fallbacks, (size_t)0);
    ASSERT_EQ(stats.critical_curves, (size_t)1);
    ASSERT(stats.critical_points >= (size_t)1);
    alea_destroy(tangent);

    /* A box slice is a polygon source.  Its four finite edges are filtered
     * and retained independently instead of falling back to the complete
     * polygon item. */
    alea_system_t* polygonal = alea_create();
    ASSERT_NOT_NULL(polygonal);
    const int box = alea_box_surface(
        polygonal, 845, -0.75, 0.75, -0.75, 0.75, -0.5, 0.5);
    ASSERT(box >= 0);
    ASSERT(alea_add_cell(
        polygonal, 845, alea_halfspace(polygonal, box, -1),
        ALEA_MATERIAL_VOID, 0.0, 0) >= 0);
    ASSERT_EQ(alea_prepare_query_acceleration(polygonal), 0);
    memset(&stats, 0, sizeof(stats));
    ASSERT_EQ(alea_transition_slice_enumerate_critical_tiles(
                  polygonal, &view, &options, &tile, 1,
                  NULL, NULL, &stats), 0);
    ASSERT_EQ(stats.critical_active_unsupported_polygon_fallbacks,
              (size_t)0);
    ASSERT_EQ(stats.critical_whole_curve_fallbacks, (size_t)0);
    ASSERT_EQ(stats.critical_active_segments, (size_t)4);
    ASSERT_EQ(stats.critical_curves, (size_t)4);
    alea_destroy(polygonal);

    /* A same-cell plane partitions the circle into a finite active half-arc;
     * the retained line/arc pair remains supported by the pair solver. */
    alea_system_t* half_circle = alea_create();
    ASSERT_NOT_NULL(half_circle);
    const int half_sphere = alea_sphere_surface(
        half_circle, 851, 0, 0, 0, 0.75);
    const int half_plane = alea_plane_surface(
        half_circle, 852, 1, 0, 0, 0);
    ASSERT(half_sphere >= 0 && half_plane >= 0);
    region = alea_intersection(
        half_circle, alea_halfspace(half_circle, half_sphere, -1),
        alea_halfspace(half_circle, half_plane, -1));
    ASSERT(alea_add_cell(half_circle, 851, region,
                         ALEA_MATERIAL_VOID, 0.0, 0) >= 0);
    ASSERT_EQ(alea_prepare_query_acceleration(half_circle), 0);
    memset(&stats, 0, sizeof(stats));
    ASSERT_EQ(alea_transition_slice_enumerate_critical_tiles(
                  half_circle, &view, &options, &tile, 1,
                  NULL, NULL, &stats), 0);
    ASSERT_EQ(stats.critical_active_segments, (size_t)2);
    ASSERT_EQ(stats.critical_curves, (size_t)2);
    ASSERT_EQ(stats.critical_whole_curve_fallbacks, (size_t)0);
    ASSERT_EQ(stats.critical_unsupported_curves, (size_t)0);
    ASSERT(stats.critical_curve_pairs_tested > 0);
    ASSERT_EQ(stats.critical_unsupported_curve_pairs, (size_t)0);
    options.max_curves_per_tile = 2;
    memset(&stats, 0, sizeof(stats));
    ASSERT_EQ(alea_transition_slice_enumerate_critical_tiles(
                  half_circle, &view, &options, &tile, 1,
                  NULL, NULL, &stats), 0);
    ASSERT_EQ(stats.critical_active_capacity_fallbacks, (size_t)0);
    ASSERT_EQ(stats.critical_whole_curve_fallbacks, (size_t)0);
    ASSERT_EQ(stats.critical_active_segments, (size_t)2);
    options.max_curves_per_tile = 512;
    alea_destroy(half_circle);

    /* General-quadric ellipse slices use the same bounded active-arc contract. */
    alea_system_t* half_ellipse = alea_create();
    ASSERT_NOT_NULL(half_ellipse);
    const int ellipsoid = alea_quadric_surface(
        half_ellipse, 861, 0.25, 1.0, 1.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.0);
    const int ellipse_plane = alea_plane_surface(
        half_ellipse, 862, 1, 0, 0, 0);
    ASSERT(ellipsoid >= 0 && ellipse_plane >= 0);
    region = alea_intersection(
        half_ellipse, alea_halfspace(half_ellipse, ellipsoid, -1),
        alea_halfspace(half_ellipse, ellipse_plane, -1));
    ASSERT(alea_add_cell(half_ellipse, 861, region,
                         ALEA_MATERIAL_VOID, 0.0, 0) >= 0);
    ASSERT_EQ(alea_prepare_query_acceleration(half_ellipse), 0);
    memset(&stats, 0, sizeof(stats));
    ASSERT_EQ(alea_transition_slice_enumerate_critical_tiles(
                  half_ellipse, &view, &options, &tile, 1,
                  NULL, NULL, &stats), 0);
    ASSERT_EQ(stats.critical_active_segments, (size_t)2);
    ASSERT_EQ(stats.critical_curves, (size_t)2);
    ASSERT_EQ(stats.critical_whole_curve_fallbacks, (size_t)0);
    ASSERT_EQ(stats.critical_unsupported_curves, (size_t)0);
    ASSERT_EQ(stats.critical_unsupported_curve_pairs, (size_t)0);
    alea_destroy(half_ellipse);

    /* Closed conics are partitioned at their exact quartic intersections;
     * neither circle/ellipse nor ellipse/ellipse needs whole-curve retention. */
    alea_system_t* circle_ellipse = alea_create();
    ASSERT_NOT_NULL(circle_ellipse);
    const int ce_sphere = alea_sphere_surface(
        circle_ellipse, 871, 0, 0, 0, 0.9);
    const int ce_ellipsoid = alea_quadric_surface(
        circle_ellipse, 872, 0.25, 4.0, 1.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.0);
    ASSERT(ce_sphere >= 0 && ce_ellipsoid >= 0);
    region = alea_intersection(
        circle_ellipse, alea_halfspace(circle_ellipse, ce_sphere, -1),
        alea_halfspace(circle_ellipse, ce_ellipsoid, -1));
    ASSERT(alea_add_cell(circle_ellipse, 871, region,
                         ALEA_MATERIAL_VOID, 0.0, 0) >= 0);
    ASSERT_EQ(alea_prepare_query_acceleration(circle_ellipse), 0);
    memset(&stats, 0, sizeof(stats));
    ASSERT_EQ(alea_transition_slice_enumerate_critical_tiles(
                  circle_ellipse, &view, &options, &tile, 1,
                  NULL, NULL, &stats), 0);
    ASSERT(stats.critical_active_segments >= (size_t)2);
    ASSERT_EQ(stats.critical_whole_curve_fallbacks, (size_t)0);
    ASSERT_EQ(stats.critical_unsupported_curves, (size_t)0);
    ASSERT_EQ(stats.critical_unsupported_curve_pairs, (size_t)0);
    alea_destroy(circle_ellipse);

    alea_system_t* ellipse_ellipse = alea_create();
    ASSERT_NOT_NULL(ellipse_ellipse);
    const int ee_horizontal = alea_quadric_surface(
        ellipse_ellipse, 881, 0.25, 2.7777777777777778, 1.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.0);
    const int ee_vertical = alea_quadric_surface(
        ellipse_ellipse, 882, 2.7777777777777778, 0.25, 1.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.0);
    ASSERT(ee_horizontal >= 0 && ee_vertical >= 0);
    region = alea_intersection(
        ellipse_ellipse,
        alea_halfspace(ellipse_ellipse, ee_horizontal, -1),
        alea_halfspace(ellipse_ellipse, ee_vertical, -1));
    ASSERT(alea_add_cell(ellipse_ellipse, 881, region,
                         ALEA_MATERIAL_VOID, 0.0, 0) >= 0);
    ASSERT_EQ(alea_prepare_query_acceleration(ellipse_ellipse), 0);
    memset(&stats, 0, sizeof(stats));
    ASSERT_EQ(alea_transition_slice_enumerate_critical_tiles(
                  ellipse_ellipse, &view, &options, &tile, 1,
                  NULL, NULL, &stats), 0);
    ASSERT(stats.critical_active_segments >= (size_t)2);
    ASSERT_EQ(stats.critical_whole_curve_fallbacks, (size_t)0);
    ASSERT_EQ(stats.critical_unsupported_curves, (size_t)0);
    ASSERT_EQ(stats.critical_unsupported_curve_pairs, (size_t)0);
    alea_destroy(ellipse_ellipse);

    /* General-conic point refinement solves parabola/parabola intersections
     * analytically through their quartic resultant. */
    alea_system_t* parabola_pair = alea_create();
    ASSERT_NOT_NULL(parabola_pair);
    const int parabola_up = alea_quadric_surface(
        parabola_pair, 891, -1.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0);
    const int parabola_down = alea_quadric_surface(
        parabola_pair, 892, 1.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 1.0, 0.0, -1.0);
    ASSERT(parabola_up >= 0 && parabola_down >= 0);
    ASSERT(alea_add_cell(parabola_pair, 891,
                         alea_halfspace(parabola_pair, parabola_up, -1),
                         ALEA_MATERIAL_VOID, 0.0, 0) >= 0);
    ASSERT(alea_add_cell(parabola_pair, 892,
                         alea_halfspace(parabola_pair, parabola_down, -1),
                         ALEA_MATERIAL_VOID, 0.0, 0) >= 0);
    ASSERT_EQ(alea_prepare_query_acceleration(parabola_pair), 0);
    memset(&stats, 0, sizeof(stats));
    ASSERT_EQ(alea_transition_slice_enumerate_critical_tiles(
                  parabola_pair, &view, &options, &tile, 1,
                  NULL, NULL, &stats), 0);
    ASSERT_EQ(stats.critical_unsupported_curves, (size_t)0);
    ASSERT_EQ(stats.critical_unsupported_curve_pairs, (size_t)0);
    ASSERT_EQ(stats.critical_active_unsupported_quartic_fallbacks, (size_t)0);
    ASSERT_EQ(stats.critical_whole_curve_fallbacks, (size_t)0);
    ASSERT(stats.critical_active_segments >= (size_t)1);
    ASSERT(stats.critical_curve_pairs_tested > (size_t)0);
    ASSERT(stats.critical_pair_intersection_points >= (size_t)2);
    ASSERT_EQ(stats.critical_active_unsupported_parabola_fallbacks,
              (size_t)0);
    ASSERT_EQ(stats.critical_whole_curve_fallbacks, (size_t)0);
    ASSERT(stats.critical_active_segments >= (size_t)2);
    alea_destroy(parabola_pair);

    /* Canonicalization is orientation-independent. */
    alea_system_t* rotated_parabola = alea_create();
    ASSERT_NOT_NULL(rotated_parabola);
    const double inv_sqrt_two = 0.70710678118654752440;
    const int rotated_parabola_surface = alea_quadric_surface(
        rotated_parabola, 911, -0.5, -0.5, 0.0,
        -1.0, 0.0, 0.0, -inv_sqrt_two, inv_sqrt_two, 0.0, 0.0);
    ASSERT(rotated_parabola_surface >= 0);
    ASSERT(alea_add_cell(rotated_parabola, 911,
                         alea_halfspace(
                             rotated_parabola, rotated_parabola_surface, -1),
                         ALEA_MATERIAL_VOID, 0.0, 0) >= 0);
    ASSERT_EQ(alea_prepare_query_acceleration(rotated_parabola), 0);
    memset(&stats, 0, sizeof(stats));
    ASSERT_EQ(alea_transition_slice_enumerate_critical_tiles(
                  rotated_parabola, &view, &options, &tile, 1,
                  NULL, NULL, &stats), 0);
    ASSERT_EQ(stats.critical_active_unsupported_parabola_fallbacks,
              (size_t)0);
    ASSERT_EQ(stats.critical_whole_curve_fallbacks, (size_t)0);
    ASSERT(stats.critical_active_segments >= (size_t)1);
    alea_destroy(rotated_parabola);

    /* Translation can make a conic's linear and constant coefficients many
     * orders larger than its quadratic coefficients.  Canonical rank tests
     * remain degree-scaled, so this distant parabola is still extracted. */
    alea_system_t* distant_parabola = alea_create();
    ASSERT_NOT_NULL(distant_parabola);
    const double distant_x = 1.0e6;
    const int distant_parabola_surface = alea_quadric_surface(
        distant_parabola, 912, 1.0, 0.0, 0.0,
        0.0, 0.0, 0.0, -2.0*distant_x, -1.0, 0.0,
        distant_x*distant_x);
    const int distant_xlo = alea_plane_surface(
        distant_parabola, 925, 1, 0, 0, -(distant_x-2.0));
    const int distant_xhi = alea_plane_surface(
        distant_parabola, 926, 1, 0, 0, -(distant_x+2.0));
    const int distant_ylo = alea_plane_surface(
        distant_parabola, 927, 0, 1, 0, 1.0);
    const int distant_yhi = alea_plane_surface(
        distant_parabola, 928, 0, 1, 0, -1.0);
    const int distant_zlo = alea_plane_surface(
        distant_parabola, 929, 0, 0, 1, 1.0);
    const int distant_zhi = alea_plane_surface(
        distant_parabola, 930, 0, 0, 1, -1.0);
    ASSERT(distant_parabola_surface >= 0 && distant_xlo >= 0 &&
           distant_xhi >= 0 && distant_ylo >= 0 && distant_yhi >= 0 &&
           distant_zlo >= 0 && distant_zhi >= 0);
    const alea_node_id_t distant_nodes[] = {
        alea_halfspace(distant_parabola, distant_parabola_surface, -1),
        alea_halfspace(distant_parabola, distant_xlo, 1),
        alea_halfspace(distant_parabola, distant_xhi, -1),
        alea_halfspace(distant_parabola, distant_ylo, 1),
        alea_halfspace(distant_parabola, distant_yhi, -1),
        alea_halfspace(distant_parabola, distant_zlo, 1),
        alea_halfspace(distant_parabola, distant_zhi, -1)
    };
    region = alea_intersection_n(
        distant_parabola, distant_nodes,
        sizeof(distant_nodes)/sizeof(distant_nodes[0]));
    ASSERT(alea_add_cell(distant_parabola, 912, region,
                         ALEA_MATERIAL_VOID, 0.0, 0) >= 0);
    ASSERT_EQ(alea_prepare_query_acceleration(distant_parabola), 0);
    alea_slice_view_t distant_view;
    alea_slice_view_init(&distant_view, 0.0, 0.0, 0.0,
                         0.0, 0.0, 1.0, 0.0, 1.0, 0.0,
                         distant_x-2.0, distant_x+2.0, -1.0, 1.0);
    alea_transition_slice_critical_tile_t distant_tile = tile;
    distant_tile.uv_min[0] = distant_x-2.0;
    distant_tile.uv_max[0] = distant_x+2.0;
    memset(&stats, 0, sizeof(stats));
    ASSERT_EQ(alea_transition_slice_enumerate_critical_tiles(
                  distant_parabola, &distant_view, &options,
                  &distant_tile, 1, NULL, NULL, &stats), 0);
    ASSERT_EQ(stats.critical_active_evaluation_general_conic_fallbacks,
              (size_t)0);
    ASSERT_EQ(stats.critical_whole_curve_fallbacks, (size_t)0);
    ASSERT(stats.critical_active_segments >= (size_t)1);
    alea_destroy(distant_parabola);

    /* A narrow off-center parabola can fail the canonical parameterization:
     * its axial linear coefficient is significant at machine precision but
     * smaller than the canonical relative tolerance.  It is neither empty
     * nor a degenerate pair of lines, so the implicit scanline path must
     * partition and retain its active branches without a whole-curve
     * fallback. */
    alea_system_t* implicit_parabola = alea_create();
    ASSERT_NOT_NULL(implicit_parabola);
    const double implicit_tangential = 200.0;
    const double implicit_axial = 1.0e-10;
    const int implicit_parabola_surface = alea_quadric_surface(
        implicit_parabola, 913, 1.0, 0.0, 0.0,
        0.0, 0.0, 0.0, implicit_tangential, -implicit_axial, 0.0,
        0.0);
    const int implicit_parabola_box = alea_box_surface(
        implicit_parabola, 914, -2.0, 2.0, -1.0, 1.0, -0.5, 0.5);
    ASSERT(implicit_parabola_surface >= 0 && implicit_parabola_box >= 0);
    region = alea_intersection(
        implicit_parabola,
        alea_halfspace(
            implicit_parabola, implicit_parabola_surface, -1),
        alea_halfspace(implicit_parabola, implicit_parabola_box, -1));
    ASSERT(alea_add_cell(implicit_parabola, 913, region,
                         ALEA_MATERIAL_VOID, 0.0, 0) >= 0);
    ASSERT_EQ(alea_prepare_query_acceleration(implicit_parabola), 0);
    alea_slice_view_t implicit_view;
    alea_slice_view_init(&implicit_view, 0.0, 0.0, 0.0,
                         0.0, 0.0, 1.0, 0.0, 1.0, 0.0,
                         -2.0, 2.0, -1.0, 1.0);
    alea_transition_slice_critical_tile_t implicit_tile = tile;
    memset(&stats, 0, sizeof(stats));
    ASSERT_EQ(alea_transition_slice_enumerate_critical_tiles(
                  implicit_parabola, &implicit_view, &options,
                  &implicit_tile, 1, NULL, NULL, &stats), 0);
    ASSERT_EQ(stats.critical_active_open_conic_canonical_fallbacks,
              (size_t)0);
    ASSERT_EQ(stats.critical_whole_curve_fallbacks, (size_t)0);
    ASSERT(stats.critical_active_segments >= (size_t)1);
    alea_destroy(implicit_parabola);

    /* The same implicit fallback must partition at intersections with a
     * general (tilted-slice) torus.  This exercises the quartic/quadratic
     * Sylvester resultant rather than the special circle-union path. */
    alea_system_t* implicit_torus = alea_create();
    ASSERT_NOT_NULL(implicit_torus);
    const int implicit_torus_surface = alea_torus_z_surface(
        implicit_torus, 915, 0.0, 0.0, 0.0, 1.5, 0.5);
    const int implicit_torus_parabola = alea_quadric_surface(
        implicit_torus, 916, 1.0, 0.0, 0.0,
        0.0, 0.0, 0.0, implicit_tangential, -implicit_axial, 0.0, 0.0);
    ASSERT(implicit_torus_surface >= 0 && implicit_torus_parabola >= 0);
    region = alea_intersection(
        implicit_torus,
        alea_halfspace(implicit_torus, implicit_torus_surface, -1),
        alea_halfspace(implicit_torus, implicit_torus_parabola, -1));
    ASSERT(alea_add_cell(implicit_torus, 915, region,
                         ALEA_MATERIAL_VOID, 0.0, 0) >= 0);
    ASSERT_EQ(alea_prepare_query_acceleration(implicit_torus), 0);
    alea_slice_view_t implicit_torus_view;
    alea_slice_view_init(&implicit_torus_view, 0.0, 0.0, 0.0,
                         0.0, 0.2, 1.0, 0.0, 0.0, 1.0,
                         -3.0, 3.0, -3.0, 3.0);
    alea_transition_slice_critical_tile_t implicit_torus_tile;
    memset(&implicit_torus_tile, 0, sizeof(implicit_torus_tile));
    implicit_torus_tile.uv_min[0] = -3.0;
    implicit_torus_tile.uv_max[0] = 3.0;
    implicit_torus_tile.uv_min[1] = -3.0;
    implicit_torus_tile.uv_max[1] = 3.0;
    memset(&stats, 0, sizeof(stats));
    ASSERT_EQ(alea_transition_slice_enumerate_critical_tiles(
                  implicit_torus, &implicit_torus_view, &options,
                  &implicit_torus_tile, 1, NULL, NULL, &stats), 0);
    ASSERT_EQ(stats.critical_active_open_conic_canonical_fallbacks,
              (size_t)0);
    ASSERT_EQ(stats.critical_whole_curve_fallbacks, (size_t)0);
    ASSERT(stats.critical_active_segments >= (size_t)1);
    alea_destroy(implicit_torus);

    /* A degenerate hyperbola is an exact pair of intersecting lines.  It is
     * factored and filtered as two components instead of failing canonical
     * positive-scale construction at its zero centered constant. */
    alea_system_t* degenerate_hyperbola = alea_create();
    ASSERT_NOT_NULL(degenerate_hyperbola);
    const int degenerate_hyperbola_surface = alea_quadric_surface(
        degenerate_hyperbola, 931, 1.0, -1.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
    const int degenerate_hyperbola_box = alea_box_surface(
        degenerate_hyperbola, 932, -1.5, 1.5, -0.75, 0.75, -0.5, 0.5);
    ASSERT(degenerate_hyperbola_surface >= 0 &&
           degenerate_hyperbola_box >= 0);
    region = alea_intersection(
        degenerate_hyperbola,
        alea_halfspace(
            degenerate_hyperbola, degenerate_hyperbola_surface, -1),
        alea_halfspace(
            degenerate_hyperbola, degenerate_hyperbola_box, -1));
    ASSERT(alea_add_cell(
        degenerate_hyperbola, 931, region,
        ALEA_MATERIAL_VOID, 0.0, 0) >= 0);
    ASSERT_EQ(alea_prepare_query_acceleration(degenerate_hyperbola), 0);
    memset(&stats, 0, sizeof(stats));
    ASSERT_EQ(alea_transition_slice_enumerate_critical_tiles(
                  degenerate_hyperbola, &view, &options, &tile, 1,
                  NULL, NULL, &stats), 0);
    ASSERT_EQ(stats.critical_active_open_conic_canonical_fallbacks,
              (size_t)0);
    ASSERT_EQ(stats.critical_whole_curve_fallbacks, (size_t)0);
    ASSERT(stats.critical_active_segments >= (size_t)2);
    alea_destroy(degenerate_hyperbola);

    /* Rectangular hyperbolas exercise the linear-in-v resultant path. */
    alea_system_t* hyperbola_pair = alea_create();
    ASSERT_NOT_NULL(hyperbola_pair);
    const int hyperbola_one = alea_quadric_surface(
        hyperbola_pair, 893, 0.0, 0.0, 0.0,
        1.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.0);
    const int hyperbola_two = alea_quadric_surface(
        hyperbola_pair, 894, 0.0, 0.0, 0.0,
        1.0, 0.0, 0.0, 1.0, 1.0, 0.0, -3.25);
    ASSERT(hyperbola_one >= 0 && hyperbola_two >= 0);
    ASSERT(alea_add_cell(hyperbola_pair, 893,
                         alea_halfspace(hyperbola_pair, hyperbola_one, -1),
                         ALEA_MATERIAL_VOID, 0.0, 0) >= 0);
    ASSERT(alea_add_cell(hyperbola_pair, 894,
                         alea_halfspace(hyperbola_pair, hyperbola_two, -1),
                         ALEA_MATERIAL_VOID, 0.0, 0) >= 0);
    ASSERT_EQ(alea_prepare_query_acceleration(hyperbola_pair), 0);
    memset(&stats, 0, sizeof(stats));
    ASSERT_EQ(alea_transition_slice_enumerate_critical_tiles(
                  hyperbola_pair, &view, &options, &tile, 1,
                  NULL, NULL, &stats), 0);
    ASSERT_EQ(stats.critical_unsupported_curves, (size_t)0);
    ASSERT_EQ(stats.critical_unsupported_curve_pairs, (size_t)0);
    ASSERT(stats.critical_pair_intersection_points >= (size_t)1);
    ASSERT_EQ(stats.critical_active_unsupported_hyperbola_fallbacks,
              (size_t)0);
    ASSERT_EQ(stats.critical_whole_curve_fallbacks, (size_t)0);
    ASSERT(stats.critical_active_segments >= (size_t)2);
    alea_destroy(hyperbola_pair);

    /* Line/general and closed/general combinations share the same exact
     * domain-filtered critical-point path. */
    alea_system_t* mixed_conics = alea_create();
    ASSERT_NOT_NULL(mixed_conics);
    const int mixed_parabola = alea_quadric_surface(
        mixed_conics, 895, -1.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0);
    const int mixed_sphere = alea_sphere_surface(
        mixed_conics, 896, 0.0, 0.0, 0.0, 0.9);
    const int mixed_plane = alea_plane_surface(
        mixed_conics, 897, 1.0, 0.0, 0.0, 0.5);
    ASSERT(mixed_parabola >= 0 && mixed_sphere >= 0 && mixed_plane >= 0);
    ASSERT(alea_add_cell(mixed_conics, 895,
                         alea_halfspace(mixed_conics, mixed_parabola, -1),
                         ALEA_MATERIAL_VOID, 0.0, 0) >= 0);
    ASSERT(alea_add_cell(mixed_conics, 896,
                         alea_halfspace(mixed_conics, mixed_sphere, -1),
                         ALEA_MATERIAL_VOID, 0.0, 0) >= 0);
    ASSERT(alea_add_cell(mixed_conics, 897,
                         alea_halfspace(mixed_conics, mixed_plane, -1),
                         ALEA_MATERIAL_VOID, 0.0, 0) >= 0);
    ASSERT_EQ(alea_prepare_query_acceleration(mixed_conics), 0);
    memset(&stats, 0, sizeof(stats));
    ASSERT_EQ(alea_transition_slice_enumerate_critical_tiles(
                  mixed_conics, &view, &options, &tile, 1,
                  NULL, NULL, &stats), 0);
    ASSERT_EQ(stats.critical_unsupported_curves, (size_t)0);
    ASSERT_EQ(stats.critical_unsupported_curve_pairs, (size_t)0);
    ASSERT(stats.critical_pair_intersection_points >= (size_t)3);
    alea_destroy(mixed_conics);

    /* Supported general-conic intersections partition both the open conic
     * and the other parameterized boundaries into finite active pieces. */
    alea_system_t* mixed_active = alea_create();
    ASSERT_NOT_NULL(mixed_active);
    const int active_parabola = alea_quadric_surface(
        mixed_active, 902, -1.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0);
    const int active_sphere = alea_sphere_surface(
        mixed_active, 903, 0.0, 0.0, 0.0, 0.9);
    const int active_plane = alea_plane_surface(
        mixed_active, 904, 1.0, 0.0, 0.0, 0.5);
    ASSERT(active_parabola >= 0 && active_sphere >= 0 && active_plane >= 0);
    region = alea_intersection(
        mixed_active,
        alea_intersection(
            mixed_active, alea_halfspace(mixed_active, active_parabola, -1),
            alea_halfspace(mixed_active, active_sphere, -1)),
        alea_halfspace(mixed_active, active_plane, -1));
    ASSERT(alea_add_cell(mixed_active, 902, region,
                         ALEA_MATERIAL_VOID, 0.0, 0) >= 0);
    ASSERT_EQ(alea_prepare_query_acceleration(mixed_active), 0);
    memset(&stats, 0, sizeof(stats));
    ASSERT_EQ(alea_transition_slice_enumerate_critical_tiles(
                  mixed_active, &view, &options, &tile, 1,
                  NULL, NULL, &stats), 0);
    ASSERT_EQ(stats.critical_active_unsupported_parabola_fallbacks,
              (size_t)0);
    ASSERT_EQ(stats.critical_whole_curve_fallbacks, (size_t)0);
    ASSERT(stats.critical_active_segments >= (size_t)2);
    alea_destroy(mixed_active);

    /* A line crossing a torus slice uses exact quartic substitution.  The
     * central slice exercises the two-concentric-circle torus mode. */
    alea_system_t* line_torus = alea_create();
    ASSERT_NOT_NULL(line_torus);
    const int torus_surface = alea_torus_z_surface(
        line_torus, 898, 0.0, 0.0, 0.0, 1.5, 0.5);
    const int torus_plane = alea_plane_surface(
        line_torus, 899, 1.0, 0.0, 0.0, 0.5);
    ASSERT(torus_surface >= 0 && torus_plane >= 0);
    ASSERT(alea_add_cell(line_torus, 898,
                         alea_halfspace(line_torus, torus_surface, -1),
                         ALEA_MATERIAL_VOID, 0.0, 0) >= 0);
    ASSERT(alea_add_cell(line_torus, 899,
                         alea_halfspace(line_torus, torus_plane, -1),
                         ALEA_MATERIAL_VOID, 0.0, 0) >= 0);
    ASSERT_EQ(alea_prepare_query_acceleration(line_torus), 0);
    memset(&stats, 0, sizeof(stats));
    ASSERT_EQ(alea_transition_slice_enumerate_critical_tiles(
                  line_torus, &view, &options, &tile, 1,
                  NULL, NULL, &stats), 0);
    ASSERT_EQ(stats.critical_unsupported_curves, (size_t)0);
    ASSERT_EQ(stats.critical_unsupported_curve_pairs, (size_t)0);
    ASSERT(stats.critical_pair_intersection_points >= (size_t)2);
    alea_destroy(line_torus);

    alea_system_t* active_line_torus = alea_create();
    ASSERT_NOT_NULL(active_line_torus);
    const int active_torus_surface = alea_torus_z_surface(
        active_line_torus, 905, 0.0, 0.0, 0.0, 1.5, 0.5);
    const int active_torus_plane = alea_plane_surface(
        active_line_torus, 906, 1.0, 0.0, 0.0, 0.5);
    ASSERT(active_torus_surface >= 0 && active_torus_plane >= 0);
    region = alea_intersection(
        active_line_torus,
        alea_halfspace(active_line_torus, active_torus_surface, -1),
        alea_halfspace(active_line_torus, active_torus_plane, -1));
    ASSERT(alea_add_cell(active_line_torus, 905, region,
                         ALEA_MATERIAL_VOID, 0.0, 0) >= 0);
    ASSERT_EQ(alea_prepare_query_acceleration(active_line_torus), 0);
    memset(&stats, 0, sizeof(stats));
    ASSERT_EQ(alea_transition_slice_enumerate_critical_tiles(
                  active_line_torus, &view, &options, &tile, 1,
                  NULL, NULL, &stats), 0);
    ASSERT_EQ(stats.critical_active_unsupported_quartic_fallbacks,
              (size_t)0);
    ASSERT_EQ(stats.critical_whole_curve_fallbacks, (size_t)0);
    ASSERT(stats.critical_active_segments >= (size_t)2);
    alea_destroy(active_line_torus);

    /* Special torus slices are exact circle unions, so their nonlinear pairs
     * do not require a higher-degree polynomial solver. */
    alea_system_t* special_torus_pairs = alea_create();
    ASSERT_NOT_NULL(special_torus_pairs);
    const int pair_torus = alea_torus_z_surface(
        special_torus_pairs, 907, 0.0, 0.0, 0.0, 1.5, 0.5);
    const int pair_ellipse = alea_quadric_surface(
        special_torus_pairs, 908, 0.6944444444444444, 0.25, 1.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.0);
    const int pair_parabola = alea_quadric_surface(
        special_torus_pairs, 909, -1.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0);
    const int pair_torus_shifted = alea_torus_z_surface(
        special_torus_pairs, 910, 0.5, 0.0, 0.0, 1.5, 0.5);
    ASSERT(pair_torus >= 0 && pair_ellipse >= 0 &&
           pair_parabola >= 0 && pair_torus_shifted >= 0);
    ASSERT(alea_add_cell(special_torus_pairs, 907,
                         alea_halfspace(
                             special_torus_pairs, pair_torus, -1),
                         ALEA_MATERIAL_VOID, 0.0, 0) >= 0);
    ASSERT(alea_add_cell(special_torus_pairs, 908,
                         alea_halfspace(
                             special_torus_pairs, pair_ellipse, -1),
                         ALEA_MATERIAL_VOID, 0.0, 0) >= 0);
    ASSERT(alea_add_cell(special_torus_pairs, 909,
                         alea_halfspace(
                             special_torus_pairs, pair_parabola, -1),
                         ALEA_MATERIAL_VOID, 0.0, 0) >= 0);
    ASSERT(alea_add_cell(special_torus_pairs, 910,
                         alea_halfspace(
                             special_torus_pairs, pair_torus_shifted, -1),
                         ALEA_MATERIAL_VOID, 0.0, 0) >= 0);
    ASSERT_EQ(alea_prepare_query_acceleration(special_torus_pairs), 0);
    memset(&stats, 0, sizeof(stats));
    ASSERT_EQ(alea_transition_slice_enumerate_critical_tiles(
                  special_torus_pairs, &view, &options, &tile, 1,
                  NULL, NULL, &stats), 0);
    ASSERT_EQ(stats.critical_unsupported_curve_pairs, (size_t)0);
    ASSERT_EQ(stats.critical_unsupported_quartic_pairs, (size_t)0);
    ASSERT(stats.critical_pair_intersection_points >= (size_t)4);
    alea_destroy(special_torus_pairs);

    alea_system_t* tilted_line_torus = alea_create();
    ASSERT_NOT_NULL(tilted_line_torus);
    const int tilted_torus_surface = alea_torus_z_surface(
        tilted_line_torus, 900, 0.0, 0.0, 0.0, 1.5, 0.5);
    const int tilted_torus_plane = alea_plane_surface(
        tilted_line_torus, 901, 1.0, 0.0, 0.0, 0.0);
    ASSERT(tilted_torus_surface >= 0 && tilted_torus_plane >= 0);
    ASSERT(alea_add_cell(tilted_line_torus, 900,
                         alea_halfspace(
                             tilted_line_torus, tilted_torus_surface, -1),
                         ALEA_MATERIAL_VOID, 0.0, 0) >= 0);
    ASSERT(alea_add_cell(tilted_line_torus, 901,
                         alea_halfspace(
                             tilted_line_torus, tilted_torus_plane, -1),
                         ALEA_MATERIAL_VOID, 0.0, 0) >= 0);
    ASSERT_EQ(alea_prepare_query_acceleration(tilted_line_torus), 0);
    alea_slice_view_t tilted_view;
    alea_slice_view_init(&tilted_view, 0.0, 0.0, 0.0,
                         0.0, 0.2, 1.0, 0.0, 0.0, 1.0,
                         -3.0, 3.0, -3.0, 3.0);
    alea_transition_slice_critical_tile_t tilted_tile;
    memset(&tilted_tile, 0, sizeof(tilted_tile));
    tilted_tile.uv_min[0] = -3.0; tilted_tile.uv_max[0] = 3.0;
    tilted_tile.uv_min[1] = -3.0; tilted_tile.uv_max[1] = 3.0;
    memset(&stats, 0, sizeof(stats));
    ASSERT_EQ(alea_transition_slice_enumerate_critical_tiles(
                  tilted_line_torus, &tilted_view, &options, &tilted_tile, 1,
                  NULL, NULL, &stats), 0);
    ASSERT_EQ(stats.critical_unsupported_curves, (size_t)0);
    ASSERT_EQ(stats.critical_unsupported_curve_pairs, (size_t)0);
    ASSERT(stats.critical_curve_pairs_tested > (size_t)0);
    ASSERT(stats.critical_pair_intersection_points >= (size_t)1);
    alea_destroy(tilted_line_torus);

    alea_slice_view_t degree_eight_view;
    alea_slice_view_init(&degree_eight_view, 0.0, 0.0, 0.0,
                         0.0, 1.0, 1.0, 0.0, 0.0, 1.0,
                         -3.0, 3.0, -3.0, 3.0);

    alea_system_t* tilted_torus_closed = alea_create();
    ASSERT_NOT_NULL(tilted_torus_closed);
    const int degree_eight_torus = alea_torus_z_surface(
        tilted_torus_closed, 913, 0.0, 0.0, 0.0, 1.5, 0.5);
    const int degree_eight_sphere = alea_sphere_surface(
        tilted_torus_closed, 914, 0.0, 0.0, 0.0, sqrt(2.5));
    ASSERT(degree_eight_torus >= 0 && degree_eight_sphere >= 0);
    ASSERT(alea_add_cell(tilted_torus_closed, 913,
                         alea_halfspace(
                             tilted_torus_closed, degree_eight_torus, -1),
                         ALEA_MATERIAL_VOID, 0.0, 0) >= 0);
    ASSERT(alea_add_cell(tilted_torus_closed, 914,
                         alea_halfspace(
                             tilted_torus_closed, degree_eight_sphere, -1),
                         ALEA_MATERIAL_VOID, 0.0, 0) >= 0);
    ASSERT_EQ(alea_prepare_query_acceleration(tilted_torus_closed), 0);
    memset(&stats, 0, sizeof(stats));
    ASSERT_EQ(alea_transition_slice_enumerate_critical_tiles(
                  tilted_torus_closed, &degree_eight_view, &options,
                  &tilted_tile, 1, NULL, NULL, &stats), 0);
    ASSERT_EQ(stats.critical_unsupported_curve_pairs, (size_t)0);
    ASSERT(stats.critical_curve_pair_candidates > (size_t)0);
    ASSERT(stats.critical_curve_pairs_tested > (size_t)0);
    ASSERT(stats.critical_pair_algebraic_points >= (size_t)4);
    ASSERT(stats.critical_pair_domain_rejections <
           stats.critical_pair_algebraic_points);
    ASSERT(stats.critical_points >= (size_t)4);
    alea_destroy(tilted_torus_closed);

    /* The same torus/sphere pair in one CSG region must be partitioned into
     * active pieces.  In particular, the general torus quartic must not fall
     * back to publishing its complete implicit curve. */
    alea_system_t* active_tilted_torus_closed = alea_create();
    ASSERT_NOT_NULL(active_tilted_torus_closed);
    const int active_degree_eight_torus = alea_torus_z_surface(
        active_tilted_torus_closed, 923, 0.0, 0.0, 0.0, 1.5, 0.5);
    const int active_degree_eight_sphere = alea_sphere_surface(
        active_tilted_torus_closed, 924, 0.0, 0.0, 0.0, sqrt(2.5));
    ASSERT(active_degree_eight_torus >= 0 &&
           active_degree_eight_sphere >= 0);
    region = alea_intersection(
        active_tilted_torus_closed,
        alea_halfspace(active_tilted_torus_closed,
                       active_degree_eight_torus, -1),
        alea_halfspace(active_tilted_torus_closed,
                       active_degree_eight_sphere, -1));
    ASSERT(alea_add_cell(active_tilted_torus_closed, 923, region,
                         ALEA_MATERIAL_VOID, 0.0, 0) >= 0);
    ASSERT_EQ(alea_prepare_query_acceleration(active_tilted_torus_closed), 0);
    memset(&stats, 0, sizeof(stats));
    ASSERT_EQ(alea_transition_slice_enumerate_critical_tiles(
                  active_tilted_torus_closed, &degree_eight_view, &options,
                  &tilted_tile, 1, NULL, NULL, &stats), 0);
    ASSERT_EQ(stats.critical_active_unsupported_quartic_fallbacks,
              (size_t)0);
    ASSERT_EQ(stats.critical_whole_curve_fallbacks, (size_t)0);
    ASSERT(stats.critical_active_segments >= (size_t)2);
    alea_destroy(active_tilted_torus_closed);

    alea_system_t* tilted_torus_open = alea_create();
    ASSERT_NOT_NULL(tilted_torus_open);
    const int degree_eight_open_torus = alea_torus_z_surface(
        tilted_torus_open, 915, 0.0, 0.0, 0.0, 1.5, 0.5);
    const int degree_eight_parabola = alea_quadric_surface(
        tilted_torus_open, 916, -1.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 2.25);
    ASSERT(degree_eight_open_torus >= 0 && degree_eight_parabola >= 0);
    ASSERT(alea_add_cell(tilted_torus_open, 915,
                         alea_halfspace(
                             tilted_torus_open, degree_eight_open_torus, -1),
                         ALEA_MATERIAL_VOID, 0.0, 0) >= 0);
    ASSERT(alea_add_cell(tilted_torus_open, 916,
                         alea_halfspace(
                             tilted_torus_open, degree_eight_parabola, -1),
                         ALEA_MATERIAL_VOID, 0.0, 0) >= 0);
    ASSERT_EQ(alea_prepare_query_acceleration(tilted_torus_open), 0);
    memset(&stats, 0, sizeof(stats));
    ASSERT_EQ(alea_transition_slice_enumerate_critical_tiles(
                  tilted_torus_open, &degree_eight_view, &options,
                  &tilted_tile, 1, NULL, NULL, &stats), 0);
    ASSERT_EQ(stats.critical_unsupported_curve_pairs, (size_t)0);
    ASSERT(stats.critical_pair_intersection_points >= (size_t)1);
    alea_destroy(tilted_torus_open);

    alea_system_t* tilted_torus_hyperbola = alea_create();
    ASSERT_NOT_NULL(tilted_torus_hyperbola);
    const double known_rho = 1.5 + sqrt(3.0) / 4.0;
    const double known_x = sqrt(known_rho * known_rho - 0.25 * 0.25);
    const int degree_eight_hyperbola_torus = alea_torus_z_surface(
        tilted_torus_hyperbola, 917, 0.0, 0.0, 0.0, 1.5, 0.5);
    const int degree_eight_hyperbola = alea_quadric_surface(
        tilted_torus_hyperbola, 918, 0.0, 0.0, 0.0,
        1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.25 * known_x);
    ASSERT(degree_eight_hyperbola_torus >= 0 && degree_eight_hyperbola >= 0);
    ASSERT(alea_add_cell(tilted_torus_hyperbola, 917,
                         alea_halfspace(tilted_torus_hyperbola,
                                        degree_eight_hyperbola_torus, -1),
                         ALEA_MATERIAL_VOID, 0.0, 0) >= 0);
    ASSERT(alea_add_cell(tilted_torus_hyperbola, 918,
                         alea_halfspace(tilted_torus_hyperbola,
                                        degree_eight_hyperbola, -1),
                         ALEA_MATERIAL_VOID, 0.0, 0) >= 0);
    ASSERT_EQ(alea_prepare_query_acceleration(tilted_torus_hyperbola), 0);
    memset(&stats, 0, sizeof(stats));
    ASSERT_EQ(alea_transition_slice_enumerate_critical_tiles(
                  tilted_torus_hyperbola, &degree_eight_view, &options,
                  &tilted_tile, 1, NULL, NULL, &stats), 0);
    ASSERT_EQ(stats.critical_unsupported_curve_pairs, (size_t)0);
    ASSERT(stats.critical_pair_intersection_points >= (size_t)2);
    alea_destroy(tilted_torus_hyperbola);

    alea_system_t* tilted_torus_pair = alea_create();
    ASSERT_NOT_NULL(tilted_torus_pair);
    const int first_degree_sixteen_torus = alea_torus_z_surface(
        tilted_torus_pair, 919, 0.0, 0.0, 0.0, 1.5, 0.5);
    const int second_degree_sixteen_torus = alea_torus_z_surface(
        tilted_torus_pair, 920, 0.0, 0.0, 0.0, 1.7, 0.5);
    ASSERT(first_degree_sixteen_torus >= 0 &&
           second_degree_sixteen_torus >= 0);
    ASSERT(alea_add_cell(tilted_torus_pair, 919,
                         alea_halfspace(tilted_torus_pair,
                                        first_degree_sixteen_torus, -1),
                         ALEA_MATERIAL_VOID, 0.0, 0) >= 0);
    ASSERT(alea_add_cell(tilted_torus_pair, 920,
                         alea_halfspace(tilted_torus_pair,
                                        second_degree_sixteen_torus, -1),
                         ALEA_MATERIAL_VOID, 0.0, 0) >= 0);
    ASSERT_EQ(alea_prepare_query_acceleration(tilted_torus_pair), 0);
    memset(&stats, 0, sizeof(stats));
    ASSERT_EQ(alea_transition_slice_enumerate_critical_tiles(
                  tilted_torus_pair, &degree_eight_view, &options,
                  &tilted_tile, 1, NULL, NULL, &stats), 0);
    ASSERT_EQ(stats.critical_unsupported_curve_pairs, (size_t)0);
    ASSERT(stats.critical_pair_intersection_points >= (size_t)4);
    alea_destroy(tilted_torus_pair);

    alea_system_t* mixed_torus_pair = alea_create();
    ASSERT_NOT_NULL(mixed_torus_pair);
    const int mixed_general_torus = alea_torus_z_surface(
        mixed_torus_pair, 921, 0.0, 0.0, 0.0, 1.5, 0.5);
    const int mixed_special_torus = alea_torus_x_surface(
        mixed_torus_pair, 922, 0.0, 0.0, 0.0, 1.5, 1.5);
    ASSERT(mixed_general_torus >= 0 && mixed_special_torus >= 0);
    ASSERT(alea_add_cell(mixed_torus_pair, 921,
                         alea_halfspace(mixed_torus_pair,
                                        mixed_general_torus, -1),
                         ALEA_MATERIAL_VOID, 0.0, 0) >= 0);
    ASSERT(alea_add_cell(mixed_torus_pair, 922,
                         alea_halfspace(mixed_torus_pair,
                                        mixed_special_torus, -1),
                         ALEA_MATERIAL_VOID, 0.0, 0) >= 0);
    ASSERT_EQ(alea_prepare_query_acceleration(mixed_torus_pair), 0);
    memset(&stats, 0, sizeof(stats));
    ASSERT_EQ(alea_transition_slice_enumerate_critical_tiles(
                  mixed_torus_pair, &degree_eight_view, &options,
                  &tilted_tile, 1, NULL, NULL, &stats), 0);
    ASSERT(stats.critical_curve_pairs_tested >= (size_t)1);
    ASSERT_EQ(stats.critical_unsupported_curve_pairs, (size_t)0);
    alea_destroy(mixed_torus_pair);

    options.max_active_boundary_tests = 0;
    memset(&stats, 0, sizeof(stats));
    ASSERT_EQ(alea_transition_slice_enumerate_critical_tiles(
                  active, &view, &options, &tile, 1, NULL, NULL, &stats), 0);
    ASSERT(stats.critical_active_boundary_fallbacks > 0);
    ASSERT(stats.critical_whole_curve_fallbacks > 0);
    ASSERT_EQ(stats.critical_curves, (size_t)3);
    alea_destroy(active);
}

TEST_MAIN()
