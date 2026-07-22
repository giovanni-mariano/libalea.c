// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/*
 * test_mcnp.c - MCNP parsing, conversion, and export tests
 */

#include "alea_test.h"
#include "alea.h"
#include "core/alea_system.h"
#include "core/alea_export.h"
#include "core/alea_universe.h"
#include "core/alea_spatial_hier.h"
#include "alea_mcnp.h"
#include "alea_openmc.h"
#include "raycast/raycast.h"
#include "alea_raycast.h"

/* ------------------------------------------------------------------------- */
/* TRCL (Cell Transformation) Tests                                          */
/* ------------------------------------------------------------------------- */

static const char* TRCL_INLINE_INPUT =
"Test TRCL\n"
"c Cell at origin with TRCL to translate to (10, 0, 0)\n"
"1 1 -1.0 -1 TRCL=(10 0 0)\n"
"2 0 1\n"
"\n"
"1 SO 5\n"
"\n"
"M1 92235.80c 1.0\n";

static const char* TRCL_TR_INPUT =
"Test TRCL with TR card\n"
"c Cell at origin with TRCL=1 referencing TR1\n"
"1 1 -1.0 -1 TRCL=1\n"
"2 0 1\n"
"\n"
"1 SO 5\n"
"\n"
"TR1 20 0 0\n"
"M1 92235.80c 1.0\n";

TEST(trcl_inline_translation) {
    FILE* f = fopen("test_trcl_inline_tmp.mcnp", "w");
    ASSERT_NOT_NULL(f);
    fprintf(f, "%s", TRCL_INLINE_INPUT);
    fclose(f);

    mcnp_model_t* model = mcnp_load("test_trcl_inline_tmp.mcnp");
    ASSERT_NOT_NULL(model);
    alea_system_t* sys = model->sys;
    ASSERT_EQ(alea_prepare_query_acceleration(sys), 0);

    /* Origin should be void (sphere translated away) */
    int origin_mat = alea_material_at(sys, 0, 0, 0);
    ASSERT_EQ(origin_mat, 0);

    /* New center (10,0,0) should have material 1 */
    int center_mat = alea_material_at(sys, 10, 0, 0);
    ASSERT_EQ(center_mat, 1);

    /* Point inside translated sphere should have material 1 */
    int near_mat = alea_material_at(sys, 10, 4, 0);
    ASSERT_EQ(near_mat, 1);

    /* Point outside translated sphere should be void */
    int far_mat = alea_material_at(sys, 10, 6, 0);
    ASSERT_EQ(far_mat, 0);

    mcnp_model_destroy(model);
    remove("test_trcl_inline_tmp.mcnp");
}

TEST(trcl_with_tr_card) {
    FILE* f = fopen("test_trcl_tr_tmp.mcnp", "w");
    ASSERT_NOT_NULL(f);
    fprintf(f, "%s", TRCL_TR_INPUT);
    fclose(f);

    mcnp_model_t* model = mcnp_load("test_trcl_tr_tmp.mcnp");
    ASSERT_NOT_NULL(model);
    alea_system_t* sys = model->sys;
    ASSERT_EQ(alea_prepare_query_acceleration(sys), 0);

    /* Origin should be void (sphere translated away) */
    int origin_mat = alea_material_at(sys, 0, 0, 0);
    ASSERT_EQ(origin_mat, 0);

    /* (20,0,0) should have material 1 */
    int center_mat = alea_material_at(sys, 20, 0, 0);
    ASSERT_EQ(center_mat, 1);

    mcnp_model_destroy(model);
    remove("test_trcl_tr_tmp.mcnp");
}

/* ------------------------------------------------------------------------- */
/* MCNP File Parsing Tests                                                    */
/* ------------------------------------------------------------------------- */

static void assert_raycast_results_equivalent(alea_system_t* sys,
                                              double ox, double oy, double oz,
                                              double dx, double dy, double dz,
                                              double t_max) {
    alea_raycast_result_t canonical;
    alea_raycast_result_t cell_aware;
    alea_raycast_result_init(&canonical);
    alea_raycast_result_init(&cell_aware);

    ASSERT_EQ(alea_raycast(sys, ox, oy, oz, dx, dy, dz, t_max, &canonical), 0);
    ASSERT_EQ(alea_raycast_cell_aware(sys, ox, oy, oz, dx, dy, dz, t_max,
                                      &cell_aware), 0);

    ASSERT_EQ(canonical.hits.count, cell_aware.hits.count);
    for (size_t i = 0; i < canonical.hits.count; i++) {
        ASSERT_NEAR(canonical.hits.data[i].t, cell_aware.hits.data[i].t, 1e-9);
        ASSERT_EQ(canonical.hits.data[i].surface_id,
                  cell_aware.hits.data[i].surface_id);
    }

    ASSERT_EQ(canonical.segments.count, cell_aware.segments.count);
    for (size_t i = 0; i < canonical.segments.count; i++) {
        ASSERT_NEAR(canonical.segments.data[i].t_enter,
                    cell_aware.segments.data[i].t_enter, 1e-9);
        ASSERT_NEAR(canonical.segments.data[i].t_exit,
                    cell_aware.segments.data[i].t_exit, 1e-9);
        ASSERT_EQ(canonical.segments.data[i].cell_id,
                  cell_aware.segments.data[i].cell_id);
        ASSERT_EQ(canonical.segments.data[i].material_id,
                  cell_aware.segments.data[i].material_id);
        ASSERT_NEAR(canonical.segments.data[i].density,
                    cell_aware.segments.data[i].density, 1e-12);
    }

    alea_raycast_result_free(&canonical);
    alea_raycast_result_free(&cell_aware);
}

/* The reusable internal path is deliberately a semantic twin of the public
 * hit-producing raycast.  Keep this comparison here because lattice crossing
 * synthesis is the part most likely to diverge during future optimizations. */
static void assert_reusable_global_raycast_equivalent(alea_system_t* sys,
                                                      double ox, double oy, double oz,
                                                      double dx, double dy, double dz,
                                                      double t_max) {
    alea_raycast_result_t canonical;
    alea_raycast_result_t reusable;
    alea_raycast_result_init(&canonical);
    alea_raycast_result_init(&reusable);

    ASSERT_EQ(alea_raycast(sys, ox, oy, oz, dx, dy, dz, t_max, &canonical), 0);

    alea_ray_t ray;
    ASSERT_EQ(alea_ray_init(&ray, ox, oy, oz, dx, dy, dz), 0);
    ASSERT_EQ(alea_raycast_global_reuse_nocache(sys, &ray, t_max, &reusable), 0);

    ASSERT_EQ(canonical.hits.count, reusable.hits.count);
    for (size_t i = 0; i < canonical.hits.count; i++) {
        ASSERT_NEAR(canonical.hits.data[i].t, reusable.hits.data[i].t, 1e-9);
        ASSERT_EQ(canonical.hits.data[i].surface_id,
                  reusable.hits.data[i].surface_id);
        ASSERT_EQ(canonical.hits.data[i].primitive_id,
                  reusable.hits.data[i].primitive_id);
    }

    ASSERT_EQ(canonical.segments.count, reusable.segments.count);
    for (size_t i = 0; i < canonical.segments.count; i++) {
        ASSERT_NEAR(canonical.segments.data[i].t_enter,
                    reusable.segments.data[i].t_enter, 1e-9);
        ASSERT_NEAR(canonical.segments.data[i].t_exit,
                    reusable.segments.data[i].t_exit, 1e-9);
        ASSERT_EQ(canonical.segments.data[i].cell_id,
                  reusable.segments.data[i].cell_id);
        ASSERT_EQ(canonical.segments.data[i].material_id,
                  reusable.segments.data[i].material_id);
        ASSERT_NEAR(canonical.segments.data[i].density,
                    reusable.segments.data[i].density, 1e-12);
        ASSERT_EQ(canonical.segments.data[i].enter_surface_id,
                  reusable.segments.data[i].enter_surface_id);
        ASSERT_EQ(canonical.segments.data[i].exit_surface_id,
                  reusable.segments.data[i].exit_surface_id);
    }

    /* Repeat with the same buffer: capacity may be retained, contents may not. */
    ASSERT_EQ(alea_raycast_global_reuse_nocache(sys, &ray, t_max, &reusable), 0);
    ASSERT_EQ(canonical.hits.count, reusable.hits.count);
    ASSERT_EQ(canonical.segments.count, reusable.segments.count);

    alea_raycast_result_free(&canonical);
    alea_raycast_result_free(&reusable);
}

static void assert_lattice_boundary_event_contract(alea_system_t* sys) {
    alea_raycast_result_t trace;
    alea_ray_boundary_event_result_t events;
    alea_raycast_result_init(&trace);
    alea_ray_boundary_event_result_init(&events);

    alea_ray_t ray;
    ASSERT_EQ(alea_ray_init(&ray, -1.5, 0, 0, 1, 0, 0), 0);
    ASSERT_EQ(alea_raycast_boundary_events_reuse_nocache(sys, &ray, 7.0,
                                                          &trace, &events), 0);

    /* Nine ownership boundaries: seven physical intersections and two
     * synthetic lattice/DDA transitions. */
    ASSERT_EQ(events.events.count, 9);
    ASSERT_EQ(events.events.data[0].kind, ALEA_RAY_BOUNDARY_EVENT_PHYSICAL);
    ASSERT_EQ(events.events.data[0].surface_id, 1);
    ASSERT_EQ(events.events.data[0].cell_before, -1);
    ASSERT_EQ(events.events.data[0].cell_after, 2);
    ASSERT_NEAR(events.events.data[0].t, 0.5, 1e-9);

    size_t synthetic_count = 0;
    for (size_t i = 0; i < events.events.count; i++) {
        if (events.events.data[i].kind == ALEA_RAY_BOUNDARY_EVENT_SYNTHETIC_LATTICE) {
            synthetic_count++;
            ASSERT_EQ(events.events.data[i].surface_id, 0);
        }
    }
    ASSERT_EQ(synthetic_count, 2);

    /* Result storage is reusable and must not retain stale events. */
    ASSERT_EQ(alea_raycast_boundary_events_reuse_nocache(sys, &ray, 7.0,
                                                          &trace, &events), 0);
    ASSERT_EQ(events.events.count, 9);

    alea_ray_boundary_event_result_free(&events);
    alea_raycast_result_free(&trace);
}

static void assert_hier_raycast_equivalent(alea_system_t* sys,
                                           double ox, double oy, double oz,
                                           double dx, double dy, double dz,
                                           double t_max) {
    alea_raycast_result_t canonical;
    alea_raycast_result_t hier;
    alea_raycast_result_init(&canonical);
    alea_raycast_result_init(&hier);

    ASSERT_EQ(alea_raycast(sys, ox, oy, oz, dx, dy, dz, t_max, &canonical), 0);
    ASSERT_EQ(alea_raycast_hier(sys, ox, oy, oz, dx, dy, dz, t_max, &hier), 0);

    ASSERT_EQ(hier.segments.count, canonical.segments.count);
    for (size_t i = 0; i < canonical.segments.count; i++) {
        ASSERT_NEAR(hier.segments.data[i].t_enter,
                    canonical.segments.data[i].t_enter, 1e-9);
        ASSERT_NEAR(hier.segments.data[i].t_exit,
                    canonical.segments.data[i].t_exit, 1e-9);
        ASSERT_EQ(hier.segments.data[i].cell_id,
                  canonical.segments.data[i].cell_id);
        ASSERT_EQ(hier.segments.data[i].material_id,
                  canonical.segments.data[i].material_id);
        ASSERT_NEAR(hier.segments.data[i].density,
                    canonical.segments.data[i].density, 1e-12);
    }

    alea_raycast_result_free(&canonical);
    alea_raycast_result_free(&hier);
}

static void assert_hier_blas_raycast_equivalent(alea_system_t* sys,
                                                double ox, double oy, double oz,
                                                double dx, double dy, double dz,
                                                double t_max) {
    alea_raycast_result_t hier;
    alea_raycast_result_t blas;
    alea_raycast_result_init(&hier);
    alea_raycast_result_init(&blas);

    ASSERT_EQ(alea_raycast_hier(sys, ox, oy, oz, dx, dy, dz, t_max, &hier), 0);
    ASSERT_EQ(alea_raycast_hier_blas_experimental(sys, ox, oy, oz,
                                                  dx, dy, dz, t_max,
                                                  &blas), 0);

    ASSERT_EQ(blas.segments.count, hier.segments.count);
    for (size_t i = 0; i < hier.segments.count; i++) {
        ASSERT_NEAR(blas.segments.data[i].t_enter,
                    hier.segments.data[i].t_enter, 1e-9);
        ASSERT_NEAR(blas.segments.data[i].t_exit,
                    hier.segments.data[i].t_exit, 1e-9);
        ASSERT_EQ(blas.segments.data[i].cell_id,
                  hier.segments.data[i].cell_id);
        ASSERT_EQ(blas.segments.data[i].material_id,
                  hier.segments.data[i].material_id);
        ASSERT_NEAR(blas.segments.data[i].density,
                    hier.segments.data[i].density, 1e-12);
    }

    alea_raycast_result_free(&hier);
    alea_raycast_result_free(&blas);
}

static void assert_hier_cell_raycast_segments_equivalent(alea_system_t* sys,
                                                         double ox, double oy, double oz,
                                                         double dx, double dy, double dz,
                                                         double t_max) {
    alea_raycast_result_t canonical;
    alea_raycast_result_t hier;
    alea_raycast_result_init(&canonical);
    alea_raycast_result_init(&hier);

    ASSERT_EQ(alea_raycast(sys, ox, oy, oz, dx, dy, dz, t_max, &canonical), 0);
    ASSERT_EQ(alea_raycast_hier_cell_aware(sys, ox, oy, oz, dx, dy, dz,
                                           t_max, &hier), 0);

    ASSERT_EQ(hier.segments.count, canonical.segments.count);
    for (size_t i = 0; i < canonical.segments.count; i++) {
        ASSERT_NEAR(hier.segments.data[i].t_enter,
                    canonical.segments.data[i].t_enter, 1e-9);
        ASSERT_NEAR(hier.segments.data[i].t_exit,
                    canonical.segments.data[i].t_exit, 1e-9);
        ASSERT_EQ(hier.segments.data[i].cell_id,
                  canonical.segments.data[i].cell_id);
        ASSERT_EQ(hier.segments.data[i].material_id,
                  canonical.segments.data[i].material_id);
        ASSERT_NEAR(hier.segments.data[i].density,
                    canonical.segments.data[i].density, 1e-12);
        ASSERT_EQ(hier.segments.data[i].enter_surface_id,
                  canonical.segments.data[i].enter_surface_id);
        ASSERT_EQ(hier.segments.data[i].exit_surface_id,
                  canonical.segments.data[i].exit_surface_id);
    }

    alea_raycast_result_free(&canonical);
    alea_raycast_result_free(&hier);
}

TEST(parse_simple_sphere) {
    mcnp_model_t* model = mcnp_load("tests/data/simple_sphere.mcnp");
    ASSERT_NOT_NULL(model);
    alea_system_t* sys = model->sys;
    ASSERT(alea_cell_count(sys) > 0);
    mcnp_model_destroy(model);
}

TEST(parse_simple_box) {
    mcnp_model_t* model = mcnp_load("tests/data/simple_box.mcnp");
    ASSERT_NOT_NULL(model);
    alea_system_t* sys = model->sys;
    ASSERT(alea_cell_count(sys) > 0);
    mcnp_model_destroy(model);
}

TEST(raycast_boundary_events_group_coincident_surfaces) {
    const char* input =
        "Coincident boundary-event surfaces\n"
        "1 1 -1.0 -1\n"
        "2 2 -1.0 1\n"
        "3 0 2\n"
        "\n"
        "1 PX 0\n"
        "2 PX 0\n"
        "\n"
        "M1 1001.80c 1.0\n"
        "M2 1001.80c 1.0\n";
    mcnp_model_t* model = mcnp_load_string(input, strlen(input));
    ASSERT_NOT_NULL(model);
    alea_system_t* sys = model->sys;
    ASSERT_EQ(alea_prepare_query_acceleration(sys), 0);

    alea_ray_t ray;
    alea_raycast_result_t trace;
    alea_ray_boundary_event_result_t events;
    ASSERT_EQ(alea_ray_init(&ray, -1, 0, 0, 1, 0, 0), 0);
    alea_raycast_result_init(&trace);
    alea_ray_boundary_event_result_init(&events);

    ASSERT_EQ(alea_raycast_boundary_events_reuse_nocache(sys, &ray, 2.0,
                                                          &trace, &events), 0);
    ASSERT_EQ(events.events.count, 1);
    ASSERT_EQ(events.events.data[0].kind, ALEA_RAY_BOUNDARY_EVENT_PHYSICAL);
    ASSERT_EQ(events.events.data[0].surface_id, 1);
    ASSERT_NEAR(events.events.data[0].t, 1.0, 1e-9);

    const alea_ray_boundary_event_options_t all_physical = {
        .include_all_coincident_physical = true
    };
    ASSERT_EQ(alea_raycast_boundary_events_with_options(
                  sys, &ray, 2.0, &all_physical, &trace, &events), 0);
    ASSERT_EQ(events.events.count, 2);
    ASSERT_EQ(events.events.data[0].surface_id, 1);
    ASSERT_EQ(events.events.data[1].surface_id, 2);
    ASSERT_NEAR(events.events.data[0].t, events.events.data[1].t, 1e-12);

    alea_ray_boundary_event_result_free(&events);
    alea_raycast_result_free(&trace);
    mcnp_model_destroy(model);
}

TEST(raycast_boundary_events_are_bidirectionally_normalized) {
    const char* input =
        "Bidirectional boundary events\n"
        "1 1 -1.0 -1\n"
        "2 0 1\n"
        "\n"
        "1 SO 1\n"
        "\n"
        "M1 1001.80c 1.0\n";
    mcnp_model_t* model = mcnp_load_string(input, strlen(input));
    ASSERT_NOT_NULL(model);
    alea_system_t* sys = model->sys;
    ASSERT_EQ(alea_prepare_query_acceleration(sys), 0);

    alea_raycast_result_t forward_trace, reverse_trace;
    alea_ray_boundary_event_result_t forward, reverse;
    alea_raycast_result_init(&forward_trace);
    alea_raycast_result_init(&reverse_trace);
    alea_ray_boundary_event_result_init(&forward);
    alea_ray_boundary_event_result_init(&reverse);
    alea_ray_t forward_ray, reverse_ray;
    ASSERT_EQ(alea_ray_init(&forward_ray, -2, 0, 0, 1, 0, 0), 0);
    ASSERT_EQ(alea_ray_init(&reverse_ray, 2, 0, 0, -1, 0, 0), 0);
    ASSERT_EQ(alea_raycast_boundary_events_reuse_nocache(
                  sys, &forward_ray, 4.0, &forward_trace, &forward), 0);
    ASSERT_EQ(alea_raycast_boundary_events_reuse_nocache(
                  sys, &reverse_ray, 4.0, &reverse_trace, &reverse), 0);

    ASSERT_EQ(forward.events.count, 2);
    ASSERT_EQ(reverse.events.count, forward.events.count);
    for (size_t i = 0; i < forward.events.count; i++) {
        const alea_ray_boundary_event_t* f = &forward.events.data[i];
        const alea_ray_boundary_event_t* r =
            &reverse.events.data[reverse.events.count - 1 - i];
        ASSERT_EQ(f->kind, ALEA_RAY_BOUNDARY_EVENT_PHYSICAL);
        ASSERT_EQ(r->kind, ALEA_RAY_BOUNDARY_EVENT_PHYSICAL);
        ASSERT_NEAR(f->t, 4.0 - r->t, 1e-9);
        ASSERT_EQ(f->surface_id, r->surface_id);
        ASSERT_EQ(f->cell_before, r->cell_after);
        ASSERT_EQ(f->cell_after, r->cell_before);
        ASSERT_EQ(f->material_before, r->material_after);
        ASSERT_EQ(f->material_after, r->material_before);
    }

    alea_ray_boundary_event_result_free(&forward);
    alea_ray_boundary_event_result_free(&reverse);
    alea_raycast_result_free(&forward_trace);
    alea_raycast_result_free(&reverse_trace);
    mcnp_model_destroy(model);
}

TEST(raycast_descends_into_fill_universe) {
    const char* input =
        "Fill ray test\n"
        "1 0 1 -2 3 -4 5 -6 FILL=10\n"
        "2 0 -1 : 2 : -3 : 4 : -5 : 6\n"
        "10 1 -1.0 -7 U=10\n"
        "11 2 -2.0 7 U=10\n"
        "\n"
        "1 PX -10\n"
        "2 PX 10\n"
        "3 PY -1\n"
        "4 PY 1\n"
        "5 PZ -1\n"
        "6 PZ 1\n"
        "7 PX 0\n"
        "\n"
        "M1 1001.80c 1.0\n"
        "M2 1001.80c 1.0\n";

    mcnp_model_t* model = mcnp_load_string(input, strlen(input));
    ASSERT_NOT_NULL(model);
    alea_system_t* sys = model->sys;
    ASSERT_EQ(alea_prepare_query_acceleration(sys), 0);

    alea_cell_hit_t hits[4];
    int n = alea_find_all_cells_at_point(sys, -5.0, 0.0, 0.0, hits, 4);
    ASSERT_EQ(n, 2);
    ASSERT_EQ(hits[1].cell_id, 10);
    n = alea_find_all_cells_at_point(sys, 5.0, 0.0, 0.0, hits, 4);
    ASSERT_EQ(n, 2);
    ASSERT_EQ(hits[1].cell_id, 11);

    alea_raycast_result_t result;
    alea_raycast_result_init(&result);
    ASSERT_EQ(alea_raycast(sys, -5.0, 0.0, 0.0, 1.0, 0.0, 0.0, 10.0, &result), 0);
    ASSERT(result.segments.count >= 2);
    ASSERT_EQ(result.segments.data[0].cell_id, 10);
    ASSERT_EQ(result.segments.data[0].material_id, 1);
    ASSERT_NEAR(result.segments.data[0].t_exit, 5.0, 1e-9);
    ASSERT_EQ(result.segments.data[1].cell_id, 11);
    ASSERT_EQ(result.segments.data[1].material_id, 2);

    assert_hier_raycast_equivalent(sys, -5.0, 0.0, 0.0, 1.0, 0.0, 0.0,
                                   10.0);
    assert_hier_cell_raycast_segments_equivalent(sys, -5.0, 0.0, 0.0,
                                                 1.0, 0.0, 0.0, 10.0);

    alea_raycast_result_free(&result);
    mcnp_model_destroy(model);
}

TEST(raycast_hier_transformed_fill_matches_flat) {
    const char* input =
        "Transformed fill ray test\n"
        "1 0 -1 FILL=3 (5 0 0)\n"
        "10 1 -1.0 -2 U=3\n"
        "11 2 -2.0 2 U=3\n"
        "\n"
        "1 SO 20.0\n"
        "2 SO 3.0\n"
        "\n"
        "M1 1001.80c 1.0\n"
        "M2 1001.80c 1.0\n";

    mcnp_model_t* model = mcnp_load_string(input, strlen(input));
    ASSERT_NOT_NULL(model);
    alea_system_t* sys = model->sys;
    ASSERT_EQ(alea_prepare_query_acceleration(sys), 0);

    assert_hier_raycast_equivalent(sys, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
                                   10.0);
    assert_hier_blas_raycast_equivalent(sys, 0.0, 0.0, 0.0,
                                        1.0, 0.0, 0.0, 10.0);
    assert_hier_cell_raycast_segments_equivalent(sys, 0.0, 0.0, 0.0,
                                                 1.0, 0.0, 0.0, 10.0);

    mcnp_model_destroy(model);
}

TEST(raycast_hier_fill_container_exit_precedes_terminal_surface) {
    const char* input =
        "Fill container boundary ray test\n"
        "1 0 -1 FILL=3\n"
        "10 1 -1.0 -2 U=3\n"
        "\n"
        "1 SO 5.0\n"
        "2 SO 50.0\n"
        "\n"
        "M1 1001.80c 1.0\n";

    mcnp_model_t* model = mcnp_load_string(input, strlen(input));
    ASSERT_NOT_NULL(model);
    alea_system_t* sys = model->sys;
    ASSERT_EQ(alea_prepare_query_acceleration(sys), 0);

    assert_hier_raycast_equivalent(sys, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
                                   10.0);
    assert_hier_cell_raycast_segments_equivalent(sys, 0.0, 0.0, 0.0,
                                                 1.0, 0.0, 0.0, 10.0);

    alea_raycast_result_t hier;
    alea_raycast_result_init(&hier);
    ASSERT_EQ(alea_raycast_hier(sys, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
                                10.0, &hier), 0);
    ASSERT(hier.segments.count >= 2);
    ASSERT_NEAR(hier.segments.data[0].t_exit, 5.0, 1e-9);
    ASSERT_EQ(hier.segments.data[0].cell_id, 10);
    ASSERT_EQ(hier.segments.data[1].cell_id, -1);

    alea_raycast_result_free(&hier);
    mcnp_model_destroy(model);
}

/* ------------------------------------------------------------------------- */
/* Deduplication Export Tests                                                 */
/* ------------------------------------------------------------------------- */

TEST(dedup_identical_spheres) {
    mcnp_model_t* model = mcnp_load("tests/data/dedup_test.mcnp");
    if (!model) {
        SKIP("Test data file not found");
    }
    alea_system_t* sys = model->sys;

    alea_build_universe_index(sys);

    /* Verify geometry before export */
    int cell_lower = alea_identify_cell_at_point(sys, 0, 0, -5);
    int cell_upper = alea_identify_cell_at_point(sys, 0, 0, 5);
    ASSERT(cell_lower >= 0);
    ASSERT(cell_upper >= 0);

    /* Export with dedup enabled */
    alea_config_t cfg = alea_get_config(sys);
    cfg.dedup = true;
    alea_set_config(sys, &cfg);
    int rc = mcnp_export_system(sys, "test_dedup_tmp.mcnp");
    ASSERT_EQ(rc, 0);

    mcnp_model_destroy(model);

    /* Roundtrip: load exported file and verify */
    mcnp_model_t* model2 = mcnp_load("test_dedup_tmp.mcnp");
    ASSERT_NOT_NULL(model2);
    alea_system_t* sys2 = model2->sys;

    alea_build_universe_index(sys2);
    cell_lower = alea_identify_cell_at_point(sys2, 0, 0, -5);
    cell_upper = alea_identify_cell_at_point(sys2, 0, 0, 5);
    ASSERT(cell_lower >= 0);
    ASSERT(cell_upper >= 0);

    mcnp_model_destroy(model2);
    remove("test_dedup_tmp.mcnp");
}

TEST(dedup_opposite_planes) {
    mcnp_model_t* model = mcnp_load("tests/data/dedup_opposite_signs.mcnp");
    if (!model) {
        SKIP("Test data file not found");
    }
    alea_system_t* sys = model->sys;

    alea_build_universe_index(sys);

    int cell_left = alea_identify_cell_at_point(sys, 0, 0, 0);
    int cell_right = alea_identify_cell_at_point(sys, 10, 0, 0);
    ASSERT(cell_left >= 0);
    ASSERT(cell_right >= 0);

    alea_config_t cfg = alea_get_config(sys);
    cfg.dedup = true;
    alea_set_config(sys, &cfg);
    int rc = mcnp_export_system(sys, "test_dedup_opposite_tmp.mcnp");
    ASSERT_EQ(rc, 0);

    mcnp_model_destroy(model);

    /* Roundtrip verification */
    mcnp_model_t* model2 = mcnp_load("test_dedup_opposite_tmp.mcnp");
    ASSERT_NOT_NULL(model2);
    alea_system_t* sys2 = model2->sys;

    alea_build_universe_index(sys2);
    cell_left = alea_identify_cell_at_point(sys2, 0, 0, 0);
    cell_right = alea_identify_cell_at_point(sys2, 10, 0, 0);
    ASSERT(cell_left >= 0);
    ASSERT(cell_right >= 0);

    mcnp_model_destroy(model2);
    remove("test_dedup_opposite_tmp.mcnp");
}

/* ------------------------------------------------------------------------- */
/* Negated plane roundtrip test                                               */
/*                                                                            */
/* Surface 1 is "P -1 0 0 -5" (x = 5 with negated normal).                  */
/* Canonicalization flips to (1,0,0,-5) with inverted=1.                     */
/* The MCNP exporter un-canonicalizes surface coefficients, so the cell      */
/* sense must be computed relative to the un-canonicalized surface.           */
/* This test verifies materials are correct after roundtrip (non-dedup).     */
/* ------------------------------------------------------------------------- */

TEST(negated_plane_roundtrip) {
    mcnp_model_t* model = mcnp_load("tests/data/mcnp_negated_plane.mcnp");
    if (!model) SKIP("Test data file not found");
    alea_system_t* sys = model->sys;

    alea_build_universe_index(sys);

    /* Verify original parse: x=0 should be in cell 1 (mat 1, x < 5) */
    int mat_left = alea_material_at(sys, 0, 0, 0);
    ASSERT_EQ(mat_left, 1);

    /* x=10 should be in cell 2 (mat 2, x > 5) */
    int mat_right = alea_material_at(sys, 10, 0, 0);
    ASSERT_EQ(mat_right, 2);

    /* Export without dedup (exercises un-canonicalization of negated plane) */
    int rc = mcnp_export_system(sys, "test_negated_plane_tmp.mcnp");
    ASSERT_EQ(rc, 0);
    mcnp_model_destroy(model);

    /* Roundtrip: re-parse the exported file */
    mcnp_model_t* model2 = mcnp_load("test_negated_plane_tmp.mcnp");
    ASSERT_NOT_NULL(model2);
    alea_system_t* sys2 = model2->sys;

    alea_build_universe_index(sys2);

    /* Materials must be preserved after roundtrip */
    int mat_left2 = alea_material_at(sys2, 0, 0, 0);
    ASSERT_MSG(mat_left2 == 1,
        "x=0 should have mat 1 (x < 5) after roundtrip");

    int mat_right2 = alea_material_at(sys2, 10, 0, 0);
    ASSERT_MSG(mat_right2 == 2,
        "x=10 should have mat 2 (x > 5) after roundtrip");

    mcnp_model_destroy(model2);
    remove("test_negated_plane_tmp.mcnp");
}

/* Also test with dedup enabled */
TEST(negated_plane_roundtrip_dedup) {
    mcnp_model_t* model = mcnp_load("tests/data/mcnp_negated_plane.mcnp");
    if (!model) SKIP("Test data file not found");
    alea_system_t* sys = model->sys;

    alea_build_universe_index(sys);

    /* Export with dedup */
    alea_config_t cfg = alea_get_config(sys);
    cfg.dedup = true;
    alea_set_config(sys, &cfg);
    int rc = mcnp_export_system(sys, "test_negated_plane_dedup_tmp.mcnp");
    ASSERT_EQ(rc, 0);
    mcnp_model_destroy(model);

    /* Roundtrip verification */
    mcnp_model_t* model2 = mcnp_load("test_negated_plane_dedup_tmp.mcnp");
    ASSERT_NOT_NULL(model2);
    alea_system_t* sys2 = model2->sys;

    alea_build_universe_index(sys2);

    int mat_left = alea_material_at(sys2, 0, 0, 0);
    ASSERT_MSG(mat_left == 1,
        "x=0 should have mat 1 (x < 5) after dedup roundtrip");

    int mat_right = alea_material_at(sys2, 10, 0, 0);
    ASSERT_MSG(mat_right == 2,
        "x=10 should have mat 2 (x > 5) after dedup roundtrip");

    mcnp_model_destroy(model2);
    remove("test_negated_plane_dedup_tmp.mcnp");
}

/* ------------------------------------------------------------------------- */
/* Lattice Tests                                                              */
/* ------------------------------------------------------------------------- */

TEST(lattice_mcnp_parse) {
    mcnp_model_t* model = mcnp_load("tests/data/mcnp_lattice.mcnp");
    if (!model) {
        SKIP("Test data file not found");
    }
    alea_system_t* sys = model->sys;

    ASSERT(alea_cell_count(sys) > 0);
    mcnp_model_destroy(model);
}

/*
 * Lattice point query test.
 *
 * 3x3 rectangular lattice, pitch 2x2x2 (cell bounded by ±1 planes).
 * Universe 1: mat 1 inside cz 0.3, mat 2 outside.
 * Universe 3: mat 3 inside cz 0.3, mat 4 outside.
 * Checkerboard fill: 1 3 1 / 3 1 3 / 1 3 1
 */
TEST(lattice_eval_point_query) {
    mcnp_model_t* model = mcnp_load("tests/data/mcnp_lattice_eval.mcnp");
    if (!model) SKIP("Test data file not found");
    alea_system_t* sys = model->sys;

    alea_build_universe_index(sys);

    int cell_id, material;

    /* Element (0,0,0) center = origin, universe 1 */
    ASSERT_EQ(alea_find_cell_lazy(sys, 0, 0, 0, &cell_id, &material, NULL), 0);
    ASSERT_EQ(material, 1);

    /* Element (0,0,0), outside cylinder → material 2 */
    ASSERT_EQ(alea_find_cell_lazy(sys, 0.5, 0.5, 0, &cell_id, &material, NULL), 0);
    ASSERT_EQ(material, 2);

    /* Element (1,0,0) center = (2,0,0), universe 3 */
    ASSERT_EQ(alea_find_cell_lazy(sys, 2, 0, 0, &cell_id, &material, NULL), 0);
    ASSERT_EQ(material, 3);

    /* Element (1,0,0), outside cylinder → material 4 */
    ASSERT_EQ(alea_find_cell_lazy(sys, 2.5, 0.5, 0, &cell_id, &material, NULL), 0);
    ASSERT_EQ(material, 4);

    /* Element (2,0,0) center = (4,0,0), universe 1 again */
    ASSERT_EQ(alea_find_cell_lazy(sys, 4, 0, 0, &cell_id, &material, NULL), 0);
    ASSERT_EQ(material, 1);

    /* Element (0,1,0) center = (0,2,0), universe 3 */
    ASSERT_EQ(alea_find_cell_lazy(sys, 0, 2, 0, &cell_id, &material, NULL), 0);
    ASSERT_EQ(material, 3);

    /* Outside lattice bounds → void */
    ASSERT_EQ(alea_find_cell_lazy(sys, 10, 0, 0, &cell_id, &material, NULL), -1);

    mcnp_model_destroy(model);
}

TEST(lattice_public_point_query_matches_deepest_hit) {
    mcnp_model_t* model = mcnp_load("tests/data/mcnp_lattice_eval.mcnp");
    if (!model) SKIP("Test data file not found");
    alea_system_t* sys = model->sys;

    alea_build_universe_index(sys);

    int cell_id = -1;
    int material = -1;

    ASSERT_EQ(alea_material_at(sys, 2, 0, 0), 3);
    ASSERT_EQ(alea_find_cell_at(sys, 2, 0, 0, &cell_id, &material), 0);
    ASSERT_EQ(cell_id, 3);
    ASSERT_EQ(material, 3);

    ASSERT_EQ(alea_material_at(sys, 2.5, 0.5, 0), 4);
    ASSERT_EQ(alea_find_cell_at(sys, 2.5, 0.5, 0, &cell_id, &material), 0);
    ASSERT_EQ(cell_id, 4);
    ASSERT_EQ(material, 4);

    mcnp_model_destroy(model);
}

/*
 * Hex lattice point query via OpenMC.
 *
 * 3x3 hex lattice, pitch 2.0, flat-top orientation.
 * Basis: e1=(2,0), e2=(1, sqrt(3)).
 * Element (0,0) at origin → univ 1, (1,0) at (2,0) → univ 3, etc.
 */
TEST(lattice_hex_eval) {
    openmc_model_t* omc = openmc_load("tests/data/openmc_hex_lattice.xml");
    if (!omc) SKIP("Test data file not found");
    alea_system_t* sys = omc->sys;

    alea_build_universe_index(sys);

    int cell_id, material;

    /* Element (0,0): center at origin, universe 1 → mat 1 inside cyl */
    ASSERT_EQ(alea_find_cell_lazy(sys, 0, 0, 0, &cell_id, &material, NULL), 0);
    ASSERT_EQ(material, 1);

    /* Element (1,0): center at (2, 0), universe 3 → mat 3 inside cyl */
    ASSERT_EQ(alea_find_cell_lazy(sys, 2, 0, 0, &cell_id, &material, NULL), 0);
    ASSERT_EQ(material, 3);

    /* Element (0,1): center at (1, sqrt(3)≈1.732), universe 3 → mat 3 */
    ASSERT_EQ(alea_find_cell_lazy(sys, 1.0, 1.732, 0, &cell_id, &material, NULL), 0);
    ASSERT_EQ(material, 3);

    /* Element (1,1): center at (3, sqrt(3)), universe 1 → mat 1 */
    ASSERT_EQ(alea_find_cell_lazy(sys, 3.0, 1.732, 0, &cell_id, &material, NULL), 0);
    ASSERT_EQ(material, 1);

    /* Outside cylinder in element (0,0) → mat 2 */
    ASSERT_EQ(alea_find_cell_lazy(sys, 0.5, 0.5, 0, &cell_id, &material, NULL), 0);
    ASSERT_EQ(material, 2);

    openmc_model_destroy(omc);
}

/*
 * Raycast through a rectangular lattice.
 *
 * Ray along x-axis at y=0, z=0 through a 3x3 lattice (pitch=2).
 * Elements (0,0), (1,0), (2,0) have universes 1, 3, 1.
 * Each has a cylinder r=0.3 at the center.
 * Expected: mat2, mat1, mat2, mat4, mat3, mat4, mat2, mat1, mat2
 */
TEST(lattice_raycast) {
    mcnp_model_t* model = mcnp_load("tests/data/mcnp_lattice_eval.mcnp");
    if (!model) SKIP("Test data file not found");
    alea_system_t* sys = model->sys;

    alea_build_universe_index(sys);
    ASSERT_EQ(alea_prepare_query_acceleration(sys), 0);
    assert_raycast_results_equivalent(sys, -1.5, 0, 0, 1, 0, 0, 7.0);
    assert_reusable_global_raycast_equivalent(sys, -1.5, 0, 0, 1, 0, 0, 7.0);
    assert_lattice_boundary_event_contract(sys);
    assert_hier_raycast_equivalent(sys, -1.5, 0, 0, 1, 0, 0, 7.0);
    assert_hier_cell_raycast_segments_equivalent(sys, -1.5, 0, 0,
                                                 1, 0, 0, 7.0);

    alea_raycast_result_t result;
    alea_raycast_result_init(&result);

    /* Ray along +x through all 3 lattice elements, y=z=0 */
    int rc = alea_raycast(sys, -1.5, 0, 0, 1, 0, 0, 7.0, &result);
    ASSERT_EQ(rc, 0);

    /* Debug: print what we got */
    printf("  hits=%zu segments=%zu\n", result.hits.count, result.segments.count);
    for (size_t i = 0; i < result.segments.count; i++) {
        printf("  seg[%zu]: t=[%.4f, %.4f] cell=%d mat=%d\n",
               i, result.segments.data[i].t_enter, result.segments.data[i].t_exit,
               result.segments.data[i].cell_id, result.segments.data[i].material_id);
    }
    for (size_t i = 0; i < result.hits.count; i++) {
        printf("  hit[%zu]: t=%.6f surf=%d\n",
               i, result.hits.data[i].t, result.hits.data[i].surface_id);
    }

    /* We should have segments with correct materials */
    ASSERT(result.segments.count >= 5);

    /* Verify path through cylinder intersections:
     * Cylinder r=0.3 at each element center (0, 2, 4).
     * Entry at x=-0.3 (t=1.2), exit at x=0.3 (t=1.8) for element 0
     * Entry at x=1.7 (t=3.2), exit at x=2.3 (t=3.8) for element 1
     * Entry at x=3.7 (t=5.2), exit at x=4.3 (t=5.8) for element 2
     */
    int found_mat1 = 0, found_mat3 = 0;
    for (size_t i = 0; i < result.segments.count; i++) {
        if (result.segments.data[i].material_id == 1) found_mat1++;
        if (result.segments.data[i].material_id == 3) found_mat3++;
    }

    /* Should find material 1 (univ 1 cylinder) and material 3 (univ 3 cylinder) */
    ASSERT(found_mat1 >= 2);  /* Elements (0,0) and (2,0) */
    ASSERT(found_mat3 >= 1);  /* Element (1,0) */

    alea_raycast_result_free(&result);
    mcnp_model_destroy(model);
}

/*
 * Raycast through a hex lattice (via OpenMC input).
 *
 * Ray along x-axis at y=0, z=0 through a 3x3 hex lattice (pitch=2).
 * Element (0,0) at origin → univ 1, element (1,0) at (2,0) → univ 3.
 * Each has a cylinder r=0.3.
 * Expected: outside cyl (univ 1), inside cyl (univ 1), outside,
 *           boundary, outside cyl (univ 3), inside cyl (univ 3), ...
 */
TEST(lattice_hex_raycast) {
    openmc_model_t* omc = openmc_load("tests/data/openmc_hex_lattice.xml");
    if (!omc) SKIP("Test data file not found");
    alea_system_t* sys = omc->sys;

    alea_build_universe_index(sys);
    ASSERT_EQ(alea_prepare_query_acceleration(sys), 0);
    assert_raycast_results_equivalent(sys, -0.5, 0, 0, 1, 0, 0, 5.0);
    assert_hier_raycast_equivalent(sys, -0.5, 0, 0, 1, 0, 0, 5.0);
    assert_hier_cell_raycast_segments_equivalent(sys, -0.5, 0, 0,
                                                 1, 0, 0, 5.0);

    alea_raycast_result_t result;
    alea_raycast_result_init(&result);

    /* Ray along +x at y=0, z=0 */
    int rc = alea_raycast(sys, -0.5, 0, 0, 1, 0, 0, 5.0, &result);
    ASSERT_EQ(rc, 0);

    printf("  hex hits=%zu segments=%zu\n",
           result.hits.count, result.segments.count);
    for (size_t i = 0; i < result.segments.count; i++) {
        printf("  seg[%zu]: t=[%.4f, %.4f] cell=%d mat=%d\n",
               i, result.segments.data[i].t_enter, result.segments.data[i].t_exit,
               result.segments.data[i].cell_id, result.segments.data[i].material_id);
    }

    /* Should find both material 1 (univ 1 cyl) and material 3 (univ 3 cyl) */
    int found_mat1 = 0, found_mat3 = 0;
    for (size_t i = 0; i < result.segments.count; i++) {
        if (result.segments.data[i].material_id == 1) found_mat1++;
        if (result.segments.data[i].material_id == 3) found_mat3++;
    }
    ASSERT(found_mat1 >= 1);
    ASSERT(found_mat3 >= 1);

    alea_raycast_result_free(&result);
    openmc_model_destroy(omc);
}

/*
 * Nested lattice (lattice of lattices) point query.
 *
 * Outer: 2x1 rectangular lattice, pitch 2 (cell 100, fill=0:1 0:0 0:0).
 *   Element (0,0,0) center at (0,0,0) → universe 10
 *   Element (1,0,0) center at (2,0,0) → universe 10
 *
 * Inner (universe 10): 2x2 rectangular lattice, pitch 1 (cell 200,
 *   fill=-1:0 -1:0 0:0).  Covers [-1,1]x[-1,1] in element-local coords.
 *   fill[0]=1 (-1,-1) center (-0.5,-0.5)
 *   fill[1]=3 (-1, 0) center (-0.5, 0.5)
 *   fill[2]=3 ( 0,-1) center ( 0.5,-0.5)
 *   fill[3]=1 ( 0, 0) center ( 0.5, 0.5)
 *
 * Universe 1: mat 1 inside cz 0.2, mat 2 outside.
 * Universe 3: mat 3 inside cz 0.2, mat 4 outside.
 */
TEST(lattice_nested) {
    mcnp_model_t* model = mcnp_load("tests/data/mcnp_nested_lattice.mcnp");
    if (!model) SKIP("Test data file not found");
    alea_system_t* sys = model->sys;

    alea_build_universe_index(sys);

    int cell_id, material;

    /* ---- Outer element (0,0,0), inner element (-1,-1) → univ 1 ---- */
    /* Center of inner element at local (-0.5,-0.5) = global (-0.5,-0.5) */
    ASSERT_EQ(alea_find_cell_lazy(sys, -0.5, -0.5, 0, &cell_id, &material, NULL), 0);
    ASSERT_EQ(material, 1);   /* inside cylinder */

    /* Same element, outside cylinder */
    ASSERT_EQ(alea_find_cell_lazy(sys, -0.8, -0.8, 0, &cell_id, &material, NULL), 0);
    ASSERT_EQ(material, 2);

    /* ---- Outer element (0,0,0), inner element (0,-1) → univ 3 ---- */
    /* Center at local (0.5,-0.5) = global (0.5,-0.5) */
    ASSERT_EQ(alea_find_cell_lazy(sys, 0.5, -0.5, 0, &cell_id, &material, NULL), 0);
    ASSERT_EQ(material, 3);

    /* ---- Outer element (0,0,0), inner element (-1,0) → univ 3 ---- */
    /* Center at local (-0.5,0.5) = global (-0.5,0.5) */
    ASSERT_EQ(alea_find_cell_lazy(sys, -0.5, 0.5, 0, &cell_id, &material, NULL), 0);
    ASSERT_EQ(material, 3);

    /* ---- Outer element (0,0,0), inner element (0,0) → univ 1 ---- */
    /* Center at local (0.5,0.5) = global (0.5,0.5) */
    ASSERT_EQ(alea_find_cell_lazy(sys, 0.5, 0.5, 0, &cell_id, &material, NULL), 0);
    ASSERT_EQ(material, 1);

    /* ---- Outer element (1,0,0), inner element (-1,-1) → univ 1 ---- */
    /* Outer center at (2,0,0). Inner center at local (-0.5,-0.5) = global (1.5,-0.5) */
    ASSERT_EQ(alea_find_cell_lazy(sys, 1.5, -0.5, 0, &cell_id, &material, NULL), 0);
    ASSERT_EQ(material, 1);

    /* ---- Outer element (1,0,0), inner element (0,-1) → univ 3 ---- */
    /* Global (2.5,-0.5) */
    ASSERT_EQ(alea_find_cell_lazy(sys, 2.5, -0.5, 0, &cell_id, &material, NULL), 0);
    ASSERT_EQ(material, 3);

    /* ---- Outside lattice bounds → void ---- */
    ASSERT_EQ(alea_find_cell_lazy(sys, 10, 0, 0, &cell_id, &material, NULL), -1);

    mcnp_model_destroy(model);
}

/*
 * Undefined-fill regions (Option C): a container whose filling universe has
 * no cell at the point is returned as the answer, flagged
 * ALEA_RESOLVE_UNDEFINED_FILL. Universe 10 covers only x < 0 of container
 * cell 1; x > 0 inside the container is an MCNP undefined region.
 */
TEST(undefined_fill_flagged) {
    const char* input =
        "Undefined fill flag test\n"
        "1 0 1 -2 3 -4 5 -6 FILL=10\n"
        "2 0 -1 : 2 : -3 : 4 : -5 : 6\n"
        "10 1 -1.0 -7 1 3 -4 5 -6 U=10\n"
        "\n"
        "1 PX -10\n"
        "2 PX 10\n"
        "3 PY -1\n"
        "4 PY 1\n"
        "5 PZ -1\n"
        "6 PZ 1\n"
        "7 PX 0\n"
        "\n"
        "M1 1001.80c 1.0\n";

    mcnp_model_t* model = mcnp_load_string(input, strlen(input));
    ASSERT_NOT_NULL(model);
    alea_system_t* sys = model->sys;
    ASSERT_EQ(alea_prepare_query_acceleration(sys), 0);

    /* Covered half: terminal cell, no flag */
    alea_cell_hit_t hit;
    ASSERT_EQ(alea_find_deepest_cell_hit_at_point(sys, -5, 0, 0, &hit), 0);
    ASSERT_EQ(hit.cell_id, 10);
    ASSERT_EQ(hit.resolution_flags, 0);

    /* Uncovered half: the container is the answer, flagged */
    ASSERT_EQ(alea_find_deepest_cell_hit_at_point(sys, 5, 0, 0, &hit), 0);
    ASSERT_EQ(hit.cell_id, 1);
    ASSERT_EQ(hit.resolution_flags & ALEA_RESOLVE_UNDEFINED_FILL,
              ALEA_RESOLVE_UNDEFINED_FILL);

    /* Hier resolver agrees, including the flag */
    alea_hier_cell_hit_t hhit;
    ASSERT_EQ(alea_hier_spatial_find_deepest_cell_at_point(sys, 5, 0, 0,
                                                           &hhit), 1);
    ASSERT_EQ(hhit.hit.cell_id, 1);
    ASSERT_EQ(hhit.hit.resolution_flags & ALEA_RESOLVE_UNDEFINED_FILL,
              ALEA_RESOLVE_UNDEFINED_FILL);

    /* Cache-poisoning order: an undefined-fill query caches a chain ending
     * in the container; a following query where the fill HAS content must
     * not be served that truncated chain (the container contains both
     * points, so containment validation alone cannot reject it). */
    alea_cell_hit_t chain[8];
    int nh = alea_find_all_cells(sys, 5, 0, 0, chain, 8);   /* undefined   */
    ASSERT_EQ(nh, 1);
    nh = alea_find_all_cells(sys, -5, 0, 0, chain, 8);       /* covered     */
    ASSERT_EQ(nh, 2);
    ASSERT_EQ(chain[1].cell_id, 10);
    ASSERT_EQ(chain[1].resolution_flags, 0);

    /* Trace segments carry the flag and agree with the point answers */
    alea_raycast_result_t r;
    alea_raycast_result_init(&r);
    ASSERT_EQ(alea_raycast_hier(sys, -5, 0, 0, 1, 0, 0, 15.0, &r), 0);
    ASSERT(r.segments.count >= 2);
    ASSERT_EQ(r.segments.data[0].cell_id, 10);
    ASSERT_EQ(r.segments.data[0].resolution_flags, 0);
    ASSERT_EQ(r.segments.data[1].cell_id, 1);
    ASSERT_EQ(r.segments.data[1].resolution_flags &
              ALEA_RESOLVE_UNDEFINED_FILL, ALEA_RESOLVE_UNDEFINED_FILL);
    alea_raycast_result_free(&r);

    /* Flat pipeline segments too */
    alea_raycast_result_init(&r);
    ASSERT_EQ(alea_raycast(sys, -5, 0, 0, 1, 0, 0, 15.0, &r), 0);
    int saw_flagged_container = 0;
    for (size_t i = 0; i < r.segments.count; i++) {
        if (r.segments.data[i].cell_id == 1) {
            saw_flagged_container = 1;
            ASSERT_EQ(r.segments.data[i].resolution_flags &
                      ALEA_RESOLVE_UNDEFINED_FILL,
                      ALEA_RESOLVE_UNDEFINED_FILL);
        }
    }
    ASSERT(saw_flagged_container);
    alea_raycast_result_free(&r);

    mcnp_model_destroy(model);
}

/*
 * Interval defect classification: an isolated overlap (cell 2 fully inside
 * cell 1, no complement) produces no ownership transition and is invisible
 * to trace segments; the owner-set classifier must still report it with
 * exact extent, and a coverage hole must classify as a gap.
 */
TEST(ray_classify_intervals) {
    const char* input =
        "Isolated overlap + gap\n"
        "1 1 -1.0 -1 U=0\n"
        "2 2 -1.0 -2 U=0\n"
        "3 3 -1.0 -4 1 3 U=0\n"
        "99 0 4\n"
        "\n"
        "1 SO 5\n"
        "2 SO 3\n"
        "3 SO 6\n"
        "4 SO 10\n"
        "\n"
        "M1 1001.80c 1.0\n"
        "M2 1001.80c 1.0\n"
        "M3 1001.80c 1.0\n";
    /* Radial structure: [0,3] cells 1+2 (overlap), (3,5] cell 1,
     * (5,6) NOTHING (gap between SO 5 and SO 6), [6,10) cell 3. */

    mcnp_model_t* model = mcnp_load_string(input, strlen(input));
    ASSERT_NOT_NULL(model);
    alea_system_t* sys = model->sys;

    alea_ray_interval_finding_t f[16];
    int n = alea_ray_classify_intervals(sys, -10, 0, 0, 1, 0, 0, 20, f, 16);
    ASSERT_EQ(n, 7);

    /* -10..-6 cell 3; -6..-5 GAP; -5..-3 cell 1; -3..3 OVERLAP(1,2);
       3..5 cell 1; 5..6 GAP; 6..10 cell 3 */
    ASSERT_EQ(f[0].kind, ALEA_INTERVAL_OK);
    ASSERT_EQ(f[0].cell_id, 3);
    ASSERT_EQ(f[1].kind, ALEA_INTERVAL_GAP);
    ASSERT_NEAR(f[1].t_enter, 4.0, 1e-6);
    ASSERT_NEAR(f[1].t_exit, 5.0, 1e-6);
    ASSERT_EQ(f[2].kind, ALEA_INTERVAL_OK);
    ASSERT_EQ(f[2].cell_id, 1);
    ASSERT_EQ(f[3].kind, ALEA_INTERVAL_OVERLAP);
    ASSERT_EQ(f[3].cell_id, 1);
    ASSERT_EQ(f[3].overlap_cell_id, 2);
    ASSERT_NEAR(f[3].t_enter, 7.0, 1e-6);
    ASSERT_NEAR(f[3].t_exit, 13.0, 1e-6);
    ASSERT_EQ(f[4].kind, ALEA_INTERVAL_OK);
    ASSERT_EQ(f[5].kind, ALEA_INTERVAL_GAP);
    ASSERT_EQ(f[6].kind, ALEA_INTERVAL_OK);
    ASSERT_EQ(f[6].cell_id, 3);

    mcnp_model_destroy(model);
}

TEST_MAIN()
