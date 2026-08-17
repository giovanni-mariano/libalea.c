// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/*
 * test_raycast.c - Unit tests for ray casting module
 */

#include "alea_test.h"
#include "alea.h"
#include "alea_geo_validator.h"
#include "alea_raycast.h"
#include "alea_slice.h"
#include "raycast/raycast.h"
#include "raycast/ray_intersect.h"
#include "core/alea_system.h"

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

TEST_MAIN()
