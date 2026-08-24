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
#include "alea_slice.h"
#include "alea_geo_validator.h"
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

/* FIRST_VISIBLE and ANY_HIT must retain the global trace's answer even when
 * reaching it requires entering a finite lattice from void.  Exercise the
 * public batch wrapper too: today it orchestrates scalar queries, so this is
 * both a policy and an API-contract regression test. */
static void assert_lattice_query_policy_equivalent(
    alea_system_t* sys,
    double ox, double oy, double oz,
    double dx, double dy, double dz,
    double t_min, double t_max,
    int material_filter) {
    alea_ray_t ray;
    alea_raycast_result_t global_scratch;
    alea_raycast_result_t hier_scratch;
    alea_ray_query_output_t global_output;
    alea_ray_first_visible_result_t hier_visible;
    alea_raycast_result_init(&global_scratch);
    alea_raycast_result_init(&hier_scratch);
    ASSERT_EQ(alea_ray_init(&ray, ox, oy, oz, dx, dy, dz), 0);

    const alea_ray_query_t query = {
        .kind = ALEA_RAY_QUERY_FIRST_VISIBLE,
        .fields = ALEA_RAY_QUERY_FIELD_SURFACE_ID |
                  ALEA_RAY_QUERY_FIELD_SURFACE_NORMAL,
        .t_min = t_min, .t_max = t_max, .material_filter = material_filter,
        .backend = ALEA_RAY_QUERY_BACKEND_GLOBAL
    };
    ASSERT_EQ(alea_raycast_query_reuse_nocache(
                  sys, &ray, &query, &global_scratch, NULL, &global_output), 0);

    const alea_ray_query_t first_cell_global = {
        .kind = ALEA_RAY_QUERY_FIRST_CELL,
        .backend = ALEA_RAY_QUERY_BACKEND_GLOBAL,
        .t_min = t_min, .t_max = t_max, .material_filter = material_filter
    };
    alea_ray_query_output_t global_first_cell;
    ASSERT_EQ(alea_raycast_query_reuse_nocache(
                  sys, &ray, &first_cell_global, &global_scratch, NULL,
                  &global_first_cell), 0);
    alea_ray_query_t first_cell_auto = first_cell_global;
    first_cell_auto.backend = ALEA_RAY_QUERY_BACKEND_AUTO;
    alea_ray_query_output_t auto_first_cell;
    ASSERT_EQ(alea_raycast_query_reuse_nocache(
                  sys, &ray, &first_cell_auto, &hier_scratch, NULL,
                  &auto_first_cell), 0);
    ASSERT_EQ(auto_first_cell.first_cell_id, global_first_cell.first_cell_id);
    ASSERT_NEAR(auto_first_cell.first_cell_t, global_first_cell.first_cell_t,
                1e-9);
    ASSERT_EQ(hier_scratch.hits.count, 0);
    ASSERT_EQ(hier_scratch.segments.count, 0);

    ASSERT_EQ(alea_raycast_hier_first_visible_nocache(
                  sys, &ray, t_min, t_max, material_filter, 1,
                  &hier_scratch, &hier_visible), 0);
    ASSERT_EQ(hier_visible.found, global_output.first_visible.found);
    if (hier_visible.found) {
        ASSERT_NEAR(hier_visible.t, global_output.first_visible.t, 1e-9);
        ASSERT_EQ(hier_visible.cell_id, global_output.first_visible.cell_id);
        ASSERT_EQ(hier_visible.material_id, global_output.first_visible.material_id);
        ASSERT_EQ(hier_visible.surface_id, global_output.first_visible.surface_id);
        ASSERT_NEAR(hier_visible.nx, global_output.first_visible.nx, 1e-9);
        ASSERT_NEAR(hier_visible.ny, global_output.first_visible.ny, 1e-9);
        ASSERT_NEAR(hier_visible.nz, global_output.first_visible.nz, 1e-9);
    }
    ASSERT_EQ(hier_scratch.segments.count, 0);

    int any_hit = 0;
    ASSERT_EQ(alea_raycast_hier_any_hit_nocache(
                  sys, &ray, t_min, t_max, material_filter,
                  &hier_scratch, &any_hit), 0);
    ASSERT_EQ(any_hit, global_output.first_visible.found ? 1 : 0);

    const alea_ray_query_t any_query = {
        .kind = ALEA_RAY_QUERY_ANY_HIT,
        .backend = ALEA_RAY_QUERY_BACKEND_AUTO,
        .t_min = t_min, .t_max = t_max, .material_filter = material_filter
    };
    alea_ray_query_output_t any_output;
    ASSERT_EQ(alea_raycast_query_reuse_nocache(
                  sys, &ray, &any_query, &hier_scratch, NULL, &any_output), 0);
    ASSERT_EQ(any_output.any_hit, any_hit != 0);
    ASSERT_EQ(hier_scratch.hits.count, 0);
    ASSERT_EQ(hier_scratch.segments.count, 0);

    /* AUTO selected intervals must use the same coherent lattice walker as
     * the explicit forward backend, rather than silently changing engines
     * when a model contains lattice placements. */
    const alea_ray_query_t segments_auto = {
        .kind = ALEA_RAY_QUERY_SEGMENTS,
        .backend = ALEA_RAY_QUERY_BACKEND_AUTO,
        .t_min = t_min, .t_max = t_max, .material_filter = material_filter
    };
    alea_ray_query_t segments_fast = segments_auto;
    segments_fast.backend = ALEA_RAY_QUERY_BACKEND_FAST_FORWARD;
    alea_raycast_result_t auto_segments;
    alea_raycast_result_t fast_segments;
    alea_raycast_result_init(&auto_segments);
    alea_raycast_result_init(&fast_segments);
    ASSERT_EQ(alea_raycast_query_reuse_nocache(
                  sys, &ray, &segments_auto, &auto_segments, NULL, NULL), 0);
    ASSERT_EQ(alea_raycast_query_reuse_nocache(
                  sys, &ray, &segments_fast, &fast_segments, NULL, NULL), 0);
    ASSERT_EQ(auto_segments.hits.count, 0);
    ASSERT_EQ(auto_segments.segments.count, fast_segments.segments.count);
    for (size_t i = 0; i < auto_segments.segments.count; i++) {
        const alea_ray_segment_t* a = &auto_segments.segments.data[i];
        const alea_ray_segment_t* b = &fast_segments.segments.data[i];
        ASSERT_NEAR(a->t_enter, b->t_enter, 1e-9);
        ASSERT_NEAR(a->t_exit, b->t_exit, 1e-9);
        ASSERT_EQ(a->cell_id, b->cell_id);
        ASSERT_EQ(a->material_id, b->material_id);
    }
    alea_raycast_result_free(&fast_segments);
    alea_raycast_result_free(&auto_segments);

    const double origins[] = {ox, oy, oz};
    const double directions[] = {dx, dy, dz};
    const double t_mins[] = {t_min};
    const double t_maxs[] = {t_max};
    const alea_ray_batch_query_t batch_query = {
        .kind = ALEA_RAY_QUERY_FIRST_VISIBLE,
        .fields = query.fields,
        .material_filter = material_filter,
        .t_mins = t_mins,
        .t_maxs = t_maxs
    };
    alea_ray_first_visible_batch_result_t batch;
    alea_ray_first_visible_batch_result_init(&batch);
    ASSERT_EQ(alea_raycast_hier_first_visible_batch_nocache(
                  sys, origins, directions, 1, &batch_query, &batch), 0);
    ASSERT_EQ(batch.found[0], global_output.first_visible.found ? 1 : 0);
    if (batch.found[0]) {
        ASSERT_NEAR(batch.t[0], global_output.first_visible.t, 1e-9);
        ASSERT_EQ(batch.cell_ids[0], global_output.first_visible.cell_id);
        ASSERT_EQ(batch.material_ids[0], global_output.first_visible.material_id);
        ASSERT_EQ(batch.surface_ids[0], global_output.first_visible.surface_id);
    }
    alea_ray_first_visible_batch_result_free(&batch);
    alea_raycast_result_free(&hier_scratch);
    alea_raycast_result_free(&global_scratch);
}

static void assert_hier_packet_batch_matches_scalar(
    alea_system_t* sys, const double* origins, const double* directions,
    size_t ray_count, const double* t_mins, const double* t_maxs,
    int material_filter) {
    const alea_ray_batch_query_t visible_query = {
        .kind = ALEA_RAY_QUERY_FIRST_VISIBLE,
        .fields = ALEA_RAY_QUERY_FIELD_SURFACE_ID |
                  ALEA_RAY_QUERY_FIELD_SURFACE_NORMAL,
        .material_filter = material_filter,
        .t_mins = t_mins, .t_maxs = t_maxs
    };
    const alea_ray_batch_query_t any_query = {
        .kind = ALEA_RAY_QUERY_ANY_HIT,
        .material_filter = material_filter,
        .t_mins = t_mins, .t_maxs = t_maxs
    };
    alea_ray_first_visible_batch_result_t visible;
    alea_ray_any_hit_batch_result_t any;
    alea_raycast_result_t scratch;
    alea_ray_first_visible_batch_result_init(&visible);
    alea_ray_any_hit_batch_result_init(&any);
    alea_raycast_result_init(&scratch);
    ASSERT_EQ(alea_raycast_hier_first_visible_batch_nocache(
                  sys, origins, directions, ray_count, &visible_query, &visible), 0);
    ASSERT_EQ(alea_raycast_hier_any_hit_batch_nocache(
                  sys, origins, directions, ray_count, &any_query, &any), 0);
    for (size_t i = 0; i < ray_count; i++) {
        alea_ray_t ray;
        alea_ray_first_visible_result_t scalar;
        ASSERT_EQ(alea_ray_init(&ray, origins[i * 3], origins[i * 3 + 1],
                                origins[i * 3 + 2], directions[i * 3],
                                directions[i * 3 + 1], directions[i * 3 + 2]), 0);
        ASSERT_EQ(alea_raycast_hier_first_visible_nocache(
                      sys, &ray, t_mins ? t_mins[i] : 0.0,
                      t_maxs ? t_maxs[i] : 0.0, material_filter, 1,
                      &scratch, &scalar), 0);
        ASSERT_EQ(visible.found[i], scalar.found ? 1 : 0);
        ASSERT_EQ(any.hits[i], scalar.found ? 1 : 0);
        if (scalar.found) {
            ASSERT_NEAR(visible.t[i], scalar.t, 1e-9);
            ASSERT_EQ(visible.cell_ids[i], scalar.cell_id);
            ASSERT_EQ(visible.material_ids[i], scalar.material_id);
            ASSERT_EQ(visible.surface_ids[i], scalar.surface_id);
        }
    }

    /* Fixed-output execution is transactional too: a preflight resource
     * rejection must leave the prior multi-ray publication usable. */
    if (ray_count > 1) {
        const uint8_t* visible_found = visible.found;
        const uint8_t* any_hits = any.hits;
        alea_ray_batch_query_t limited_visible = visible_query;
        alea_ray_batch_query_t limited_any = any_query;
        limited_visible.max_output_bytes = 1;
        limited_any.max_output_bytes = ray_count - 1;
        ASSERT_EQ(alea_raycast_hier_first_visible_batch_nocache(
                      sys, origins, directions, ray_count, &limited_visible,
                      &visible), -1);
        ASSERT_EQ(visible.found, visible_found);
        ASSERT_EQ(alea_raycast_hier_any_hit_batch_nocache(
                      sys, origins, directions, ray_count, &limited_any,
                      &any), -1);
        ASSERT_EQ(any.hits, any_hits);
    }
    alea_raycast_result_free(&scratch);
    alea_ray_any_hit_batch_result_free(&any);
    alea_ray_first_visible_batch_result_free(&visible);
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

    /* Ten ownership boundaries: seven physical intersections, two internal
     * DDA transitions, and the finite lattice-array exit. */
    ASSERT_EQ(events.events.count, 10);
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
    ASSERT_EQ(synthetic_count, 3);

    /* Result storage is reusable and must not retain stale events. */
    ASSERT_EQ(alea_raycast_boundary_events_reuse_nocache(sys, &ray, 7.0,
                                                          &trace, &events), 0);
    ASSERT_EQ(events.events.count, 10);

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

    const alea_ray_boundary_event_options_internal_t all_physical = {
        .include_all_coincident_physical = true
    };
    ASSERT_EQ(alea_raycast_boundary_events_with_options(
                  sys, &ray, 2.0, &all_physical, &trace, &events), 0);
    ASSERT_EQ(events.events.count, 2);
    ASSERT_EQ(events.events.data[0].surface_id, 1);
    ASSERT_EQ(events.events.data[1].surface_id, 2);
    ASSERT_NEAR(events.events.data[0].t, events.events.data[1].t, 1e-12);

    /* Boundary-event batches stage compact records in worker arenas but must
     * retain input order and t-min clipping through transactional compaction. */
    const double batch_origins[] = {
        -1, 0, 0, -1, 0, 0, -1, 0, 0
    };
    const double batch_directions[] = {
        1, 0, 0, 1, 0, 0, 1, 0, 0
    };
    const double batch_t_mins[] = {0.0, 1.1, 0.0};
    const double batch_t_maxs[] = {2.0, 2.0, 2.0};
    const alea_ray_batch_query_t batch_query = {
        .kind = ALEA_RAY_QUERY_BOUNDARY_EVENTS,
        .material_filter = -1,
        .t_mins = batch_t_mins, .t_maxs = batch_t_maxs
    };
    alea_ray_boundary_event_batch_result_t batch = {0};
    ASSERT_EQ(alea_raycast_boundary_events_batch_nocache(
                  sys, batch_origins, batch_directions, 3, &batch_query,
                  &batch), 0);
    ASSERT_EQ(batch.event_count, (size_t)2);
    ASSERT_EQ(batch.ray_offsets[0], (uint64_t)0);
    ASSERT_EQ(batch.ray_offsets[1], (uint64_t)1);
    ASSERT_EQ(batch.ray_offsets[2], (uint64_t)1);
    ASSERT_EQ(batch.ray_offsets[3], (uint64_t)2);
    ASSERT_NEAR(batch.t[0], 1.0, 1e-9);
    ASSERT_NEAR(batch.t[1], 1.0, 1e-9);
    const uint64_t* batch_offsets = batch.ray_offsets;
    alea_ray_batch_query_t limited_batch_query = batch_query;
    limited_batch_query.max_events = 1;
    ASSERT_EQ(alea_raycast_boundary_events_batch_nocache(
                  sys, batch_origins, batch_directions, 3,
                  &limited_batch_query, &batch), -1);
    ASSERT_EQ(batch.ray_offsets, batch_offsets);
    ASSERT_EQ(batch.event_count, (size_t)2);
    alea_ray_boundary_event_batch_result_free(&batch);

    alea_ray_boundary_event_result_free(&events);
    alea_raycast_result_free(&trace);
    mcnp_model_destroy(model);
}

TEST(raycast_boundary_events_support_many_coincident_surfaces) {
    char input[8192];
    size_t used = (size_t)snprintf(input, sizeof(input),
        "Many coincident boundary-event surfaces\n"
        "1 1 -1.0 -1\n"
        "2 2 -1.0 1\n"
        "3 0 2\n\n");
    ASSERT(used < sizeof(input));
    for (int surface_id = 1; surface_id <= 65; surface_id++) {
        int written = snprintf(input + used, sizeof(input) - used,
                               "%d PX 0\n", surface_id);
        ASSERT(written > 0 && (size_t)written < sizeof(input) - used);
        used += (size_t)written;
    }
    int written = snprintf(input + used, sizeof(input) - used,
                           "\nM1 1001.80c 1.0\nM2 1001.80c 1.0\n");
    ASSERT(written > 0 && (size_t)written < sizeof(input) - used);

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

    const alea_ray_boundary_event_options_internal_t all_physical = {
        .include_all_coincident_physical = true
    };
    ASSERT_EQ(alea_raycast_boundary_events_with_options(
                  sys, &ray, 2.0, &all_physical, &trace, &events), 0);
    ASSERT_EQ(events.events.count, 65);
    ASSERT_EQ(events.events.data[0].surface_id, 1);
    ASSERT_EQ(events.events.data[64].surface_id, 65);

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

TEST(raycast_internal_query_policies_match_canonical_trace) {
    const char* input =
        "Internal ray query policies\n"
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

    alea_ray_t ray;
    alea_raycast_result_t trace;
    alea_ray_boundary_event_result_t events;
    alea_ray_query_output_t output;
    ASSERT_EQ(alea_ray_init(&ray, -2, 0, 0, 1, 0, 0), 0);
    alea_raycast_result_init(&trace);
    alea_ray_boundary_event_result_init(&events);

    const alea_ray_query_t visible = {
        .kind = ALEA_RAY_QUERY_FIRST_VISIBLE,
        .fields = ALEA_RAY_QUERY_FIELD_SURFACE_NORMAL |
                  ALEA_RAY_QUERY_FIELD_SURFACE_ID,
        .t_min = 0,
        .t_max = 4,
        .material_filter = -1
    };
    ASSERT_EQ(alea_raycast_query_reuse_nocache(
                  sys, &ray, &visible, &trace, NULL, &output), 0);
    ASSERT(output.first_visible.found);
    ASSERT_EQ(output.first_visible.cell_id, 1);
    ASSERT_EQ(output.first_visible.material_id, 1);
    ASSERT_EQ(output.first_visible.surface_id, 1);
    ASSERT_NEAR(output.first_visible.t, 1.0, 1e-9);
    ASSERT_NEAR(fabs(output.first_visible.nx), 1.0, 1e-9);
    ASSERT_EQ(trace.hits.count, 0);
    ASSERT_EQ(trace.segments.count, 0);
    const int first_visible_steps = trace.step_iterations;

    alea_ray_query_t fast_visible = visible;
    fast_visible.backend = ALEA_RAY_QUERY_BACKEND_FAST_FORWARD;
    ASSERT_EQ(alea_raycast_query_reuse_nocache(
                  sys, &ray, &fast_visible, &trace, NULL, &output), 0);
    ASSERT(output.first_visible.found);
    ASSERT_EQ(output.first_visible.cell_id, 1);
    ASSERT_EQ(output.first_visible.surface_id, 1);
    ASSERT_NEAR(fabs(output.first_visible.nx), 1.0, 1e-9);
    ASSERT_EQ(trace.hits.count, 0);
    ASSERT_EQ(trace.segments.count, 0);

    alea_raycast_result_t full_trace;
    alea_raycast_result_init(&full_trace);
    ASSERT_EQ(alea_raycast_hier_with_hits_nocache(sys, &ray, 4.0,
                                                   &full_trace), 0);
    ASSERT(full_trace.segments.count > 0);
    ASSERT_NEAR(full_trace.segments.data[1].t_enter,
                output.first_visible.t, 1e-9);
    ASSERT_EQ(full_trace.segments.data[1].cell_id,
              output.first_visible.cell_id);
    ASSERT(full_trace.step_iterations > first_visible_steps);
    alea_raycast_result_free(&full_trace);

    int first_cell_id = -1;
    double first_cell_t = -1;
    ASSERT_EQ(alea_raycast_hier_first_cell_nocache(
                  sys, &ray, 0, 4.0, -1,
                  &trace, &first_cell_id, &first_cell_t), 0);
    ASSERT_EQ(first_cell_id, 2);
    ASSERT_NEAR(first_cell_t, 0.0, 1e-9);
    ASSERT_EQ(trace.hits.count, 0);
    ASSERT_EQ(trace.segments.count, 0);
    const alea_ray_query_t first_cell = {
        .kind = ALEA_RAY_QUERY_FIRST_CELL,
        .t_min = 0, .t_max = 4, .material_filter = -1
    };
    ASSERT_EQ(alea_raycast_query_reuse_nocache(
                  sys, &ray, &first_cell, &trace, NULL, &output), 0);
    ASSERT_EQ(output.first_cell_id, 2);
    ASSERT_NEAR(output.first_cell_t, 0.0, 1e-9);
    ASSERT_EQ(trace.segments.count, 0);
    alea_ray_query_t fast_first_cell = first_cell;
    fast_first_cell.backend = ALEA_RAY_QUERY_BACKEND_FAST_FORWARD;
    ASSERT_EQ(alea_raycast_query_reuse_nocache(
                  sys, &ray, &fast_first_cell, &trace, NULL, &output), 0);
    ASSERT_EQ(output.first_cell_id, 2);
    ASSERT_NEAR(output.first_cell_t, 0.0, 1e-9);
    ASSERT_EQ(trace.segments.count, 0);
    const alea_ray_query_t filtered_first_cell = {
        .kind = ALEA_RAY_QUERY_FIRST_CELL,
        .t_min = 1.5, .t_max = 4, .material_filter = 1
    };
    ASSERT_EQ(alea_raycast_query_reuse_nocache(
                  sys, &ray, &filtered_first_cell, &trace, NULL, &output), 0);
    ASSERT_EQ(output.first_cell_id, 1);
    ASSERT_NEAR(output.first_cell_t, 1.5, 1e-9);
    ASSERT_EQ(trace.segments.count, 0);

    alea_ray_first_visible_options_t public_visible;
    alea_ray_first_visible_options_init(&public_visible);
    public_visible.fields = ALEA_RAY_FIRST_VISIBLE_SURFACE_ID |
                            ALEA_RAY_FIRST_VISIBLE_SURFACE_NORMAL;
    public_visible.t_max = 4;
    alea_ray_first_visible_query_result_t* public_result =
        alea_ray_first_visible_query_result_create();
    ASSERT_NOT_NULL(public_result);
    ASSERT_EQ(alea_ray_first_visible_query(sys, -2, 0, 0, 1, 0, 0,
                                           &public_visible, public_result), 0);
    ASSERT(alea_ray_first_visible_found(public_result));
    ASSERT_EQ(alea_ray_first_visible_cell_id(public_result), 1);
    ASSERT_NEAR(alea_ray_first_visible_t(public_result), 1.0, 1e-9);
    ASSERT_EQ(alea_ray_first_visible_surface_id(public_result), 1);
    double public_nx, public_ny, public_nz;
    ASSERT_EQ(alea_ray_first_visible_normal(public_result,
                                             &public_nx, &public_ny, &public_nz), 0);
    ASSERT_NEAR(fabs(public_nx), 1.0, 1e-9);
    /* Older callers may provide only the stable options prefix; omitted
     * material_filter must retain its documented default. */
    public_visible.material_filter = 99;
    public_visible.struct_size = offsetof(alea_ray_first_visible_options_t,
                                          material_filter);
    ASSERT_EQ(alea_ray_first_visible_query(sys, -2, 0, 0, 1, 0, 0,
                                           &public_visible, public_result), 0);
    ASSERT(alea_ray_first_visible_found(public_result));
    alea_ray_first_visible_query_result_destroy(public_result);

    alea_ray_query_t no_normal = visible;
    no_normal.fields = ALEA_RAY_QUERY_FIELD_SURFACE_ID;
    ASSERT_EQ(alea_raycast_query_reuse_nocache(
                  sys, &ray, &no_normal, &trace, NULL, &output), 0);
    ASSERT(output.first_visible.found);
    ASSERT_EQ(output.first_visible.surface_id, 1);
    ASSERT_NEAR(output.first_visible.nx, 0.0, 1e-12);
    ASSERT_NEAR(output.first_visible.ny, 0.0, 1e-12);
    ASSERT_NEAR(output.first_visible.nz, 0.0, 1e-12);

    int occluded = 0;
    ASSERT_EQ(alea_raycast_hier_any_hit_nocache(
                  sys, &ray, 0, 4, -1, &trace, &occluded), 0);
    ASSERT_EQ(occluded, 1);
    ASSERT_EQ(trace.hits.count, 0);
    ASSERT_EQ(trace.segments.count, 0);

    alea_ray_query_t clipped = visible;
    clipped.t_min = 1.5;
    ASSERT_EQ(alea_raycast_query_reuse_nocache(
                  sys, &ray, &clipped, &trace, NULL, &output), 0);
    ASSERT(output.first_visible.found);
    ASSERT_NEAR(output.first_visible.t, 1.5, 1e-9);
    ASSERT_EQ(output.first_visible.surface_id, -1);

    const alea_ray_query_t any_hit = {
        .kind = ALEA_RAY_QUERY_ANY_HIT,
        .t_min = 0, .t_max = 4, .material_filter = 1
    };
    ASSERT_EQ(alea_raycast_query_reuse_nocache(
                  sys, &ray, &any_hit, &trace, NULL, &output), 0);
    ASSERT(output.any_hit);
    ASSERT_EQ(trace.hits.count, 0);
    ASSERT_EQ(trace.segments.count, 0);

    /* Scalar SEGMENTS must have the same range/cross-section semantics as
     * the compact per-ray executor. */
    const alea_ray_query_t segments = {
        .kind = ALEA_RAY_QUERY_SEGMENTS,
        .t_min = 1.5, .t_max = 4, .material_filter = -1
    };
    ASSERT_EQ(alea_raycast_query_reuse_nocache(
                  sys, &ray, &segments, &trace, NULL, &output), 0);
    ASSERT_EQ(trace.hits.count, 0);
    ASSERT_EQ(trace.segments.count, 2);
    ASSERT_NEAR(trace.segments.data[0].t_enter, 1.5, 1e-9);
    ASSERT_NEAR(trace.segments.data[0].t_exit, 3.0, 1e-9);
    ASSERT_EQ(trace.segments.data[0].cell_id, 1);
    ASSERT_EQ(trace.segments.data[0].enter_surface_id, -1);
    ASSERT_EQ(trace.segments.data[0].enter_hit_index, -1);
    ASSERT_NEAR(trace.segments.data[1].t_enter, 3.0, 1e-9);
    ASSERT_NEAR(trace.segments.data[1].t_exit, 4.0, 1e-9);
    ASSERT_EQ(trace.segments.data[1].cell_id, 2);

    alea_ray_query_t fast_segments = segments;
    fast_segments.backend = ALEA_RAY_QUERY_BACKEND_FAST_FORWARD;
    ASSERT_EQ(alea_raycast_query_reuse_nocache(
                  sys, &ray, &fast_segments, &trace, NULL, &output), 0);
    ASSERT_EQ(trace.segments.count, 2);
    ASSERT_NEAR(trace.segments.data[0].t_enter, 1.5, 1e-9);
    ASSERT_NEAR(trace.segments.data[0].t_exit, 3.0, 1e-9);
    ASSERT_EQ(trace.segments.data[0].cell_id, 1);

    fast_segments.backend = ALEA_RAY_QUERY_BACKEND_FAST_REVERSE;
    ASSERT_EQ(alea_raycast_query_reuse_nocache(
                  sys, &ray, &fast_segments, &trace, NULL, &output), 0);
    ASSERT_EQ(trace.segments.count, 2);
    ASSERT_NEAR(trace.segments.data[0].t_enter, 1.5, 1e-9);
    ASSERT_NEAR(trace.segments.data[0].t_exit, 3.0, 1e-9);
    ASSERT_EQ(trace.segments.data[0].cell_id, 1);

    fast_segments.backend = ALEA_RAY_QUERY_BACKEND_FAST_FORWARD_REVERSE;
    ASSERT_EQ(alea_raycast_query_reuse_nocache(
                  sys, &ray, &fast_segments, &trace, NULL, &output), 0);
    ASSERT(!output.directional_mismatch);
    ASSERT_EQ(trace.segments.count, 2);
    ASSERT_EQ(trace.segments.data[0].cell_id, 1);

    const alea_ray_query_t boundaries = {
        .kind = ALEA_RAY_QUERY_BOUNDARY_EVENTS,
        .t_min = 0, .t_max = 4, .material_filter = -1
    };
    ASSERT_EQ(alea_raycast_query_reuse_nocache(
                  sys, &ray, &boundaries, &trace, &events, &output), 0);
    ASSERT_EQ(events.events.count, 2);
    ASSERT_EQ(events.events.data[0].surface_id, 1);

    /* Event limits apply to the requested ray range, not crossings clipped
     * away by t_min. */
    alea_ray_query_t clipped_boundaries = boundaries;
    clipped_boundaries.t_min = 2.0;
    clipped_boundaries.max_events = 1;
    ASSERT_EQ(alea_raycast_query_reuse_nocache(
                  sys, &ray, &clipped_boundaries, &trace, &events, &output), 0);
    ASSERT_EQ(events.events.count, 1);
    ASSERT(events.events.data[0].t >= clipped_boundaries.t_min);

    alea_ray_query_t fast_boundaries = boundaries;
    fast_boundaries.backend = ALEA_RAY_QUERY_BACKEND_FAST_FORWARD_REVERSE;
    ASSERT_EQ(alea_raycast_query_reuse_nocache(
                  sys, &ray, &fast_boundaries, &trace, &events, &output), 0);
    ASSERT_EQ(events.events.count, 2);
    ASSERT_EQ(events.events.data[0].surface_id, 1);

    alea_ray_query_t clipped_fast_boundaries = fast_boundaries;
    clipped_fast_boundaries.t_min = 2.0;
    clipped_fast_boundaries.max_events = 1;
    ASSERT_EQ(alea_raycast_query_reuse_nocache(
                  sys, &ray, &clipped_fast_boundaries,
                  &trace, &events, &output), 0);
    ASSERT_EQ(events.events.count, 1);
    ASSERT(events.events.data[0].t >= clipped_fast_boundaries.t_min);

    alea_ray_boundary_event_options_t public_events;
    alea_ray_boundary_event_options_init(&public_events);
    public_events.fields = ALEA_RAY_BOUNDARY_EVENT_PRIMITIVE_ID |
                           ALEA_RAY_BOUNDARY_EVENT_NORMAL;
    public_events.t_max = 4;
    alea_ray_boundary_event_query_result_t* public_event_result =
        alea_ray_boundary_event_query_result_create();
    ASSERT_NOT_NULL(public_event_result);
    ASSERT_EQ(alea_ray_boundary_event_query(sys, -2, 0, 0, 1, 0, 0,
                                            &public_events, public_event_result), 0);
    ASSERT_EQ(alea_ray_boundary_event_count(public_event_result), 2);
    int public_surface, public_kind;
    ASSERT_EQ(alea_ray_boundary_event_get(public_event_result, 0, NULL,
                                          &public_kind, &public_surface, NULL, NULL,
                                          NULL, NULL, NULL, NULL, NULL, NULL, NULL), 0);
    ASSERT_EQ(public_kind, ALEA_RAY_EVENT_PHYSICAL);
    ASSERT_EQ(public_surface, 1);
    public_events.max_events = 1;
    ASSERT_EQ(alea_ray_boundary_event_query(sys, -2, 0, 0, 1, 0, 0,
                                            &public_events, public_event_result), -1);
    ASSERT_EQ(alea_ray_boundary_event_count(public_event_result), 0);
    /* An older options prefix omits max_events and must retain its default. */
    public_events.struct_size = offsetof(alea_ray_boundary_event_options_t,
                                         max_events);
    ASSERT_EQ(alea_ray_boundary_event_query(sys, -2, 0, 0, 1, 0, 0,
                                            &public_events, public_event_result), 0);
    ASSERT_EQ(alea_ray_boundary_event_count(public_event_result), 2);
    alea_ray_boundary_event_query_result_destroy(public_event_result);

    alea_ray_query_t bounded = boundaries;
    bounded.max_events = 1;
    ASSERT_EQ(alea_raycast_query_reuse_nocache(
                  sys, &ray, &bounded, &trace, &events, &output), -1);
    ASSERT_EQ(events.events.count, 0);
    ASSERT_EQ(trace.segments.count, 0);

    alea_ray_boundary_event_result_free(&events);
    alea_raycast_result_free(&trace);
    mcnp_model_destroy(model);
}

TEST(surface_boundary_map_keeps_coincident_surface_labels) {
    const char* input =
        "Coincident surface labels\n"
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

    enum { width = 8, height = 1 };
    int ids[width * height];
    alea_slice_view_t view;
    alea_slice_view_axis(&view, 2, 0.0, -1.0, 1.0, -0.1, 0.1);
    ASSERT_EQ(alea_find_cells_grid(sys, &view, width, height, -1,
                                   ids, NULL, NULL), 0);

    alea_slice_surface_boundary_map_t* map = NULL;
    ASSERT_EQ(alea_slice_surface_boundary_map_create(
                  sys, &view, width, height, ids,
                  alea_slice_classify_cell, NULL, &map), 0);
    ASSERT_NOT_NULL(map);

    int found = 0;
    for (int x = 0; x + 1 < width; x++) {
        if (ids[x] == ids[x + 1]) continue;
        /* Coincident boundary definitions are a geometric overlap. The map
         * must still retain both candidate surface IDs for diagnostics/labels. */
        ASSERT_EQ(alea_slice_surface_boundary_status(
                      map, x, 0, ALEA_SLICE_EDGE_RIGHT),
                  ALEA_SLICE_BOUNDARY_OVERLAP);
        ASSERT_EQ(alea_slice_surface_boundary_surface_count(
                      map, x, 0, ALEA_SLICE_EDGE_RIGHT), 2);
        ASSERT_EQ(alea_slice_surface_boundary_surface_id(
                      map, x, 0, ALEA_SLICE_EDGE_RIGHT, 0), 1);
        ASSERT_EQ(alea_slice_surface_boundary_surface_id(
                      map, x, 0, ALEA_SLICE_EDGE_RIGHT, 1), 2);
        found++;
    }
    ASSERT_EQ(found, 1);

    alea_slice_surface_boundary_map_free(map);
    alea_slice_directional_trace_cache_t* cache =
        alea_slice_directional_trace_cache_create(sys, &view, width, height);
    ASSERT_NOT_NULL(cache);
    ASSERT_EQ(alea_slice_directional_trace_cache_matches(
                  cache, sys, &view, width, height), 1);
    alea_slice_view_t shifted_view = view;
    shifted_view.u_max += 0.01;
    ASSERT_EQ(alea_slice_directional_trace_cache_matches(
                  cache, sys, &shifted_view, width, height), 0);
    const alea_ray_boundary_event_t* events = NULL;
    size_t event_count = 0;
    ASSERT_EQ(alea_slice_directional_event_cache_line_events(
                  cache, ALEA_SLICE_EDGE_RIGHT, 0, 0, &events, &event_count), 0);
    ASSERT(event_count >= 2);
    ASSERT_NOT_NULL(events);
    ASSERT_EQ(events[0].surface_id, 1);
    ASSERT_EQ(events[1].surface_id, 2);
    alea_slice_surface_boundary_map_t* shared_map = NULL;
    ASSERT_EQ(alea_slice_surface_boundary_map_create_with_directional_cache(
                  sys, &view, width, height, ids, alea_slice_classify_cell,
                  NULL, cache, &shared_map), 0);
    ASSERT_NOT_NULL(shared_map);
    alea_slice_surface_boundary_map_free(shared_map);
    alea_ray_slice_validation_options_t validation_options;
    alea_ray_slice_validation_options_init(&validation_options);
    validation_options.checks = ALEA_RAY_SLICE_VALIDATE_FAST_BIDIRECTIONAL;
    validation_options.flags = ALEA_RAY_SLICE_VALIDATION_INCLUDE_AGREEMENTS;
    alea_ray_slice_validation_result_t* validation =
        alea_ray_slice_validation_result_create();
    ASSERT_NOT_NULL(validation);
    ASSERT_EQ(alea_validate_ray_slice_compact_with_directional_cache(
                  sys, &view, height, &validation_options, NULL, NULL,
                  cache, validation), 0);
    ASSERT_EQ(alea_ray_slice_validation_reused_trace_mask(validation),
              ALEA_RAY_SLICE_TRACE_FAST_FORWARD |
              ALEA_RAY_SLICE_TRACE_FAST_REVERSE);
    ASSERT_EQ(alea_ray_slice_validation_executed_trace_mask(validation), 0);
    ASSERT(alea_ray_slice_validation_interval_count(validation) > 0);
    ASSERT_NOT_NULL(
        alea_ray_slice_validation_u_enter_forward_surface_ids(validation));
    ASSERT_NOT_NULL(
        alea_ray_slice_validation_u_enter_reverse_surface_ids(validation));
    ASSERT_NOT_NULL(
        alea_ray_slice_validation_u_enter_provenance_flags(validation));
    const int32_t* provenance_enter =
        alea_ray_slice_validation_u_enter_forward_surface_ids(validation);
    const int32_t* provenance_exit =
        alea_ray_slice_validation_u_exit_forward_surface_ids(validation);
    const uint32_t* provenance_enter_flags =
        alea_ray_slice_validation_u_enter_provenance_flags(validation);
    const uint32_t* provenance_exit_flags =
        alea_ray_slice_validation_u_exit_provenance_flags(validation);
    int found_coincident_boundary = 0;
    for (size_t i = 0; i < alea_ray_slice_validation_interval_count(validation); i++) {
        if ((provenance_enter[i] == 1 &&
             (provenance_enter_flags[i] & ALEA_RAY_SLICE_BOUNDARY_PROVENANCE_COINCIDENT)) ||
            (provenance_exit[i] == 1 &&
             (provenance_exit_flags[i] & ALEA_RAY_SLICE_BOUNDARY_PROVENANCE_COINCIDENT))) {
            found_coincident_boundary = 1;
            break;
        }
    }
    ASSERT(found_coincident_boundary);
    /* Explicit provenance reuse rejects a cache with different sampling
     * dimensions and preserves the last published diagnostic result. */
    alea_slice_directional_trace_cache_t* incompatible_cache =
        alea_slice_directional_trace_cache_create(sys, &view, width, height + 1);
    ASSERT_NOT_NULL(incompatible_cache);
    const size_t published_interval_count =
        alea_ray_slice_validation_interval_count(validation);
    ASSERT_EQ(alea_validate_ray_slice_compact_with_directional_cache(
                  sys, &view, height, &validation_options, NULL, NULL,
                  incompatible_cache, validation), -1);
    ASSERT_EQ(alea_ray_slice_validation_interval_count(validation),
              published_interval_count);
    alea_slice_directional_trace_cache_destroy(incompatible_cache);
    ASSERT(alea_sphere_surface(sys, 99, 10, 0, 0, 1) >= 0);
    ASSERT_EQ(alea_slice_directional_trace_cache_matches(
                  cache, sys, &view, width, height), 0);
    /* A stale public cache is rejected transactionally: the last diagnostic
     * publication remains available to the caller. */
    ASSERT_EQ(alea_validate_ray_slice_compact_with_directional_cache(
                  sys, &view, height, &validation_options, NULL, NULL,
                  cache, validation), -1);
    ASSERT_EQ(alea_ray_slice_validation_interval_count(validation),
              published_interval_count);
    alea_ray_slice_validation_result_destroy(validation);
    alea_slice_directional_trace_cache_destroy(cache);
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

TEST(lattice_entry_respects_transformed_fill_ancestors) {
    const char* input =
        "Inactive transformed lattice placement\n"
        /* The fill is translated to x=5, but its enclosing cell only exists
         * at x<0.  Its lattice BVH box still overlaps x=[4,6], so DDA entry
         * must reject that placement through its enclosing-chain check. */
        "1 0 -1 FILL=10 (5 0 0)\n"
        "2 0 1\n"
        "100 0 2 -3 4 -5 6 -7 LAT=1 U=10 FILL=1:1 0:0 0:0\n"
        "     1\n"
        "10 1 -1.0 -8 U=1\n"
        "11 2 -1.0 8 U=1\n"
        "\n"
        "1 PX 0\n"
        "2 PX -1\n"
        "3 PX 1\n"
        "4 PY -1\n"
        "5 PY 1\n"
        "6 PZ -1\n"
        "7 PZ 1\n"
        "8 CZ 0.3\n"
        "\n"
        "M1 1001.80c 1.0\n"
        "M2 1001.80c 1.0\n";

    mcnp_model_t* model = mcnp_load_string(input, strlen(input));
    ASSERT_NOT_NULL(model);
    alea_system_t* sys = model->sys;
    ASSERT_EQ(alea_prepare_query_acceleration(sys), 0);
    ASSERT_EQ(alea_material_at(sys, 5.0, 0.0, 0.0), 0);
    assert_lattice_query_policy_equivalent(
        sys, 3.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 4.0, -1);

    mcnp_model_destroy(model);
}

TEST(repeating_lattice_support_bounds_prune_inactive_transformed_fill) {
    const char* input =
        "Inactive transformed repeating lattice placement\n"
        /* Cell 1 owns a finite sphere at the origin, but its filled universe is translated to x=5.
         * The simple LAT=1 FILL=1 form is infinite in that child universe.
         * Its occurrence must still be pruned for a ray wholly in x>0 from
         * the trusted support bounds of cell 1, before lattice DDA begins. */
        "1 0 -1 FILL=10 (5 0 0)\n"
        "100 0 -8 9 -10 11 -12 13 LAT=1 U=10 FILL=1\n"
        "10 1 -1.0 -14 U=1\n"
        "11 2 -1.0 14 U=1\n"
        "\n"
        "1 SO 1\n"
        "2 PX -1\n"
        "3 PY -1\n"
        "4 PY 1\n"
        "5 PZ -1\n"
        "6 PZ 1\n"
        "8 PX -1\n"
        "9 PX 1\n"
        "10 PY -1\n"
        "11 PY 1\n"
        "12 PZ -1\n"
        "13 PZ 1\n"
        "14 CZ 0.3\n"
        "\n"
        "M1 1001.80c 1.0\n"
        "M2 1001.80c 1.0\n";

    mcnp_model_t* model = mcnp_load_string(input, strlen(input));
    ASSERT_NOT_NULL(model);
    alea_system_t* sys = model->sys;
    ASSERT_EQ(alea_prepare_query_acceleration(sys), 0);

    alea_raycast_result_t trace;
    alea_raycast_result_init(&trace);
    ASSERT_EQ(alea_raycast_hier_fast_segments(
                  sys, 3.0, 0.0, 0.0, 1.0, 0.0, 0.0, 10000.0, &trace), 0);
    ASSERT(trace.lattice_entry_calls > 0);
    ASSERT_EQ(trace.lattice_entry_candidates, 0);
    ASSERT_EQ(trace.lattice_entry_dda_steps, 0);
    ASSERT_EQ(trace.lattice_entry_future_entry_results, 0);

    alea_raycast_result_free(&trace);
    mcnp_model_destroy(model);
}

TEST(repeating_lattice_jumps_to_exact_ancestor_support_entry) {
    const char* input =
        "Repeating lattice must not DDA through inactive support bounds\n"
        /* The support AABB spans x=[-100,100], but cell 1 is only the thin
         * spherical shell 99 < r < 100.  A ray from the origin reaches the
         * actual occurrence at x=99.  The lattice pitch is one, so the old
         * entry search walked roughly one hundred irrelevant pitches first. */
        "1 0 -1 2 FILL=10\n"
        "100 0 -3 LAT=1 U=10 FILL=1\n"
        "10 1 -1.0 -9 U=1\n"
        "\n"
        "1 SO 100\n"
        "2 SO 99\n"
        "3 RPP -0.5 0.5 -0.5 0.5 -0.5 0.5\n"
        "9 SO 0.25\n"
        "\n"
        "M1 1001.80c 1.0\n";

    mcnp_model_t* model = mcnp_load_string(input, strlen(input));
    ASSERT_NOT_NULL(model);
    alea_system_t* sys = model->sys;
    ASSERT_EQ(alea_prepare_query_acceleration(sys), 0);

    alea_raycast_result_t trace;
    alea_raycast_result_init(&trace);
    ASSERT_EQ(alea_raycast_hier_fast_segments(
                  sys, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 110.0, &trace), 0);
    ASSERT(trace.lattice_entry_calls > 0);
    ASSERT(trace.lattice_entry_candidates > 0);
    ASSERT(trace.lattice_entry_ancestor_surface_tests > 0);
    ASSERT(trace.lattice_entry_ancestor_events > 0);
    ASSERT(trace.lattice_entry_dda_steps < 10);

    alea_raycast_result_free(&trace);
    mcnp_model_destroy(model);
}

TEST(repeating_lattice_entry_respects_canonical_competing_owner) {
    const char* input =
        "Overlapping root fill must shadow later lattice occurrence\n"
        /* Cell 1 is deck-first and contains the ray, but its child universe
         * has no matching cell. Cell 2 has the same CSG and fills a repeating
         * lattice. The entry helper must not manufacture a transition into
         * cell 2 merely because its enclosing CSG also contains the point. */
        "1 0 -1 FILL=20\n"
        "2 0 -1 FILL=10\n"
        "100 0 -3 LAT=1 U=10 FILL=1\n"
        "10 1 -1.0 -4 U=1\n"
        "200 0 -5 U=20\n"
        "\n"
        "1 SO 10\n"
        "3 RPP -0.5 0.5 -0.5 0.5 -0.5 0.5\n"
        "4 SO 0.25\n"
        "5 S 1000 0 0 1\n"
        "\n"
        "M1 1001.80c 1.0\n";

    mcnp_model_t* model = mcnp_load_string(input, strlen(input));
    ASSERT_NOT_NULL(model);
    alea_system_t* sys = model->sys;
    ASSERT_EQ(alea_prepare_query_acceleration(sys), 0);

    alea_raycast_result_t trace;
    alea_raycast_result_init(&trace);
    ASSERT_EQ(alea_raycast_hier_fast_segments(
                  sys, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 5.0, &trace), 0);
    ASSERT(trace.lattice_entry_calls > 0);
    ASSERT(trace.lattice_entry_candidates > 0);
    ASSERT(trace.lattice_entry_canonical_rejections > 0);
    ASSERT_EQ(trace.lattice_entry_future_entry_results, 0);

    for (size_t i = 0; i < trace.segments.count; i++) {
        ASSERT_NE(trace.segments.data[i].cell_id, 10);
    }

    alea_raycast_result_free(&trace);
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

TEST(lattice_compact_first_visible_and_any_hit_match_canonical_trace) {
    mcnp_model_t* model = mcnp_load("tests/data/mcnp_lattice_eval.mcnp");
    if (!model) SKIP("Test data file not found");
    alea_system_t* sys = model->sys;
    const double origins[] = {-5.0, 0.0, 0.0};
    const double directions[] = {1.0, 0.0, 0.0};
    const double t_mins[] = {0.0};
    const double t_maxs[] = {12.0};
    const alea_ray_batch_query_t visible_query = {
        .kind = ALEA_RAY_QUERY_FIRST_VISIBLE,
        .fields = ALEA_RAY_QUERY_FIELD_SURFACE_NORMAL |
                  ALEA_RAY_QUERY_FIELD_SURFACE_ID,
        .material_filter = -1,
        .t_mins = t_mins,
        .t_maxs = t_maxs
    };
    alea_ray_first_visible_batch_result_t visible;
    alea_ray_first_visible_batch_result_init(&visible);
    ASSERT_EQ(alea_raycast_hier_first_visible_batch_nocache(
                  sys, origins, directions, 1, &visible_query, &visible), 0);
    ASSERT_EQ(visible.ray_count, 1);

    alea_ray_t ray;
    alea_raycast_result_t trace;
    alea_ray_query_output_t scalar;
    alea_raycast_result_init(&trace);
    ASSERT_EQ(alea_ray_init(&ray, origins[0], origins[1], origins[2],
                            directions[0], directions[1], directions[2]), 0);
    const alea_ray_query_t scalar_query = {
        .kind = ALEA_RAY_QUERY_FIRST_VISIBLE,
        .fields = visible_query.fields,
        .t_min = 0.0, .t_max = 12.0, .material_filter = -1,
        .backend = ALEA_RAY_QUERY_BACKEND_GLOBAL
    };
    ASSERT_EQ(alea_raycast_query_reuse_nocache(
                  sys, &ray, &scalar_query, &trace, NULL, &scalar), 0);
    ASSERT_EQ(visible.found[0], scalar.first_visible.found ? 1 : 0);
    if (visible.found[0]) {
        ASSERT_NEAR(visible.t[0], scalar.first_visible.t, 1e-9);
        ASSERT_EQ(visible.cell_ids[0], scalar.first_visible.cell_id);
        ASSERT_EQ(visible.material_ids[0], scalar.first_visible.material_id);
    }

    /* The hierarchical early-stop path must be equivalent to the canonical
     * full lattice trace, without retaining segments behind the answer. */
    alea_raycast_result_t hier_scratch;
    alea_ray_first_visible_result_t hier_visible;
    alea_raycast_result_init(&hier_scratch);
    ASSERT_EQ(alea_raycast_hier_first_visible_nocache(
                  sys, &ray, 0.0, 12.0, -1, 1,
                  &hier_scratch, &hier_visible), 0);
    ASSERT_EQ(hier_visible.found, scalar.first_visible.found);
    if (hier_visible.found) {
        ASSERT_NEAR(hier_visible.t, scalar.first_visible.t, 1e-9);
        ASSERT_EQ(hier_visible.cell_id, scalar.first_visible.cell_id);
        ASSERT_EQ(hier_visible.material_id, scalar.first_visible.material_id);
        ASSERT_EQ(hier_visible.surface_id, scalar.first_visible.surface_id);
        ASSERT_NEAR(hier_visible.nx, scalar.first_visible.nx, 1e-9);
        ASSERT_NEAR(hier_visible.ny, scalar.first_visible.ny, 1e-9);
        ASSERT_NEAR(hier_visible.nz, scalar.first_visible.nz, 1e-9);
    }
    ASSERT_EQ(hier_scratch.segments.count, 0);

    const alea_ray_batch_query_t any_query = {
        .kind = ALEA_RAY_QUERY_ANY_HIT,
        .material_filter = -1,
        .t_mins = t_mins,
        .t_maxs = t_maxs
    };
    alea_ray_any_hit_batch_result_t any;
    alea_ray_any_hit_batch_result_init(&any);
    ASSERT_EQ(alea_raycast_hier_any_hit_batch_nocache(
                  sys, origins, directions, 1, &any_query, &any), 0);
    ASSERT_EQ(any.hits[0], visible.found[0]);
    alea_ray_any_hit_batch_result_free(&any);
    alea_raycast_result_free(&hier_scratch);
    alea_raycast_result_free(&trace);
    alea_ray_first_visible_batch_result_free(&visible);
    mcnp_model_destroy(model);
}

TEST(lattice_rect_query_policies_match_global_forward_and_reverse) {
    mcnp_model_t* model = mcnp_load("tests/data/mcnp_lattice_eval.mcnp");
    if (!model) SKIP("Test data file not found");
    alea_system_t* sys = model->sys;
    ASSERT_EQ(alea_prepare_query_acceleration(sys), 0);

    assert_lattice_query_policy_equivalent(
        sys, -5.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 12.0, -1);
    assert_lattice_query_policy_equivalent(
        sys, 5.5, 0.0, 0.0, -1.0, 0.0, 0.0, 0.0, 7.0, -1);
    assert_lattice_query_policy_equivalent(
        sys, -5.0, 0.0, 0.0, 1.0, 0.0, 0.0, 2.2, 10.0, 3);
    assert_lattice_query_policy_equivalent(
        sys, 5.5, 0.0, 0.0, -1.0, 0.0, 0.0, 1.2, 7.0, 1);

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
    /* Exercise negative DDA stepping through the same synthetic boundaries. */
    assert_raycast_results_equivalent(sys, 5.5, 0, 0, -1, 0, 0, 7.0);
    assert_hier_raycast_equivalent(sys, 5.5, 0, 0, -1, 0, 0, 7.0);
    assert_hier_cell_raycast_segments_equivalent(sys, 5.5, 0, 0,
                                                 -1, 0, 0, 7.0);

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
    assert_raycast_results_equivalent(sys, 4.5, 0, 0, -1, 0, 0, 5.0);
    assert_hier_raycast_equivalent(sys, 4.5, 0, 0, -1, 0, 0, 5.0);
    assert_hier_cell_raycast_segments_equivalent(sys, 4.5, 0, 0,
                                                 -1, 0, 0, 5.0);

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

    /* Filter past the first element so this parity case crosses a synthetic
     * hex-DDA transition before reaching the visible material. */
    alea_ray_t ray;
    alea_raycast_result_t global_trace;
    alea_raycast_result_t hier_scratch;
    alea_ray_query_output_t global_output;
    alea_ray_first_visible_result_t hier_visible;
    alea_raycast_result_init(&global_trace);
    alea_raycast_result_init(&hier_scratch);
    ASSERT_EQ(alea_ray_init(&ray, -0.5, 0.0, 0.0, 1.0, 0.0, 0.0), 0);
    const alea_ray_query_t global_query = {
        .kind = ALEA_RAY_QUERY_FIRST_VISIBLE,
        .fields = ALEA_RAY_QUERY_FIELD_SURFACE_ID |
                  ALEA_RAY_QUERY_FIELD_SURFACE_NORMAL,
        .t_min = 0.0, .t_max = 5.0, .material_filter = 3,
        .backend = ALEA_RAY_QUERY_BACKEND_GLOBAL
    };
    ASSERT_EQ(alea_raycast_query_reuse_nocache(
                  sys, &ray, &global_query, &global_trace,
                  NULL, &global_output), 0);
    ASSERT_EQ(alea_raycast_hier_first_visible_nocache(
                  sys, &ray, 0.0, 5.0, 3, 1,
                  &hier_scratch, &hier_visible), 0);
    ASSERT_EQ(hier_visible.found, global_output.first_visible.found);
    ASSERT(hier_visible.found);
    ASSERT_NEAR(hier_visible.t, global_output.first_visible.t, 1e-9);
    ASSERT_EQ(hier_visible.cell_id, global_output.first_visible.cell_id);
    ASSERT_EQ(hier_visible.material_id, global_output.first_visible.material_id);
    ASSERT_EQ(hier_visible.surface_id, global_output.first_visible.surface_id);
    ASSERT_EQ(hier_scratch.segments.count, 0);
    alea_raycast_result_free(&hier_scratch);
    alea_raycast_result_free(&global_trace);

    alea_raycast_result_free(&result);
    openmc_model_destroy(omc);
}

TEST(lattice_hex_query_policies_match_global_forward_and_reverse) {
    openmc_model_t* omc = openmc_load("tests/data/openmc_hex_lattice.xml");
    if (!omc) SKIP("Test data file not found");
    alea_system_t* sys = omc->sys;
    ASSERT_EQ(alea_prepare_query_acceleration(sys), 0);

    assert_lattice_query_policy_equivalent(
        sys, -0.5, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 5.0, -1);
    assert_lattice_query_policy_equivalent(
        sys, 4.5, 0.0, 0.0, -1.0, 0.0, 0.0, 0.0, 5.0, -1);
    assert_lattice_query_policy_equivalent(
        sys, -0.5, 0.0, 0.0, 1.0, 0.0, 0.0, 1.2, 5.0, 3);
    assert_lattice_query_policy_equivalent(
        sys, 4.5, 0.0, 0.0, -1.0, 0.0, 0.0, 0.8, 5.0, 1);

    const double origins[] = {
        -0.5, 0.0, 0.0, 4.5, 0.0, 0.0, -0.5, 0.0, 0.0,
        4.5, 0.0, 0.0, -0.5, 3.5, 0.0
    };
    const double directions[] = {
        1.0, 0.0, 0.0, -1.0, 0.0, 0.0, 1.0, 0.0, 0.0,
        -1.0, 0.0, 0.0, 1.0, 0.0, 0.0
    };
    const double t_mins[] = {0.0, 0.0, 1.2, 0.8, 0.0};
    const double t_maxs[] = {5.0, 5.0, 5.0, 5.0, 5.0};
    assert_hier_packet_batch_matches_scalar(
        sys, origins, directions, 5, t_mins, t_maxs, 3);

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
    ASSERT_EQ(alea_prepare_query_acceleration(sys), 0);

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

    /* Exercise the scalar early-stop path before packetizing it. The complete
     * trace reference is covered separately: nested lattice ray stepping must
     * first prove it can descend and accept a terminal interval. */
    alea_ray_t nested_ray;
    alea_raycast_result_t nested_scratch;
    alea_ray_first_visible_result_t nested_visible;
    alea_raycast_result_init(&nested_scratch);
    ASSERT_EQ(alea_ray_init(&nested_ray, -2.0, -0.5, 0.0,
                            1.0, 0.0, 0.0), 0);
    ASSERT_EQ(alea_raycast_hier_first_visible_nocache(
                  sys, &nested_ray, 0.0, 6.0, -1, 1,
                  &nested_scratch, &nested_visible), 0);
    ASSERT(nested_visible.found);
    ASSERT_EQ(nested_scratch.segments.count, 0);
    const alea_ray_first_visible_result_t scalar_visible = nested_visible;

    const alea_ray_query_t nested_fast_query = {
        .kind = ALEA_RAY_QUERY_FIRST_VISIBLE,
        .fields = ALEA_RAY_QUERY_FIELD_SURFACE_ID |
                  ALEA_RAY_QUERY_FIELD_SURFACE_NORMAL,
        .t_min = 0.0, .t_max = 6.0, .material_filter = -1,
        .backend = ALEA_RAY_QUERY_BACKEND_FAST_FORWARD
    };
    alea_ray_query_output_t nested_fast_output;
    ASSERT_EQ(alea_raycast_query_reuse_nocache(
                  sys, &nested_ray, &nested_fast_query, &nested_scratch,
                  NULL, &nested_fast_output), 0);
    ASSERT_EQ(nested_fast_output.first_visible.found, scalar_visible.found);
    ASSERT_NEAR(nested_fast_output.first_visible.t, scalar_visible.t, 1e-9);
    ASSERT_EQ(nested_fast_output.first_visible.cell_id, scalar_visible.cell_id);
    ASSERT_EQ(nested_fast_output.first_visible.material_id,
              scalar_visible.material_id);
    ASSERT_EQ(nested_fast_output.first_visible.surface_id,
              scalar_visible.surface_id);
    alea_raycast_result_free(&nested_scratch);

    const double nested_origins[] = {-2.0, -0.5, 0.0};
    const double nested_directions[] = {1.0, 0.0, 0.0};
    const double nested_t_mins[] = {0.0};
    const double nested_t_maxs[] = {6.0};
    const alea_ray_batch_query_t nested_batch_query = {
        .kind = ALEA_RAY_QUERY_FIRST_VISIBLE,
        .fields = ALEA_RAY_QUERY_FIELD_SURFACE_ID |
                  ALEA_RAY_QUERY_FIELD_SURFACE_NORMAL,
        .material_filter = -1,
        .t_mins = nested_t_mins,
        .t_maxs = nested_t_maxs
    };
    alea_ray_first_visible_batch_result_t nested_batch_visible;
    alea_ray_first_visible_batch_result_init(&nested_batch_visible);
    ASSERT_EQ(alea_raycast_hier_first_visible_batch_nocache(
                  sys, nested_origins, nested_directions, 1,
                  &nested_batch_query, &nested_batch_visible), 0);
    ASSERT_EQ(nested_batch_visible.found[0], scalar_visible.found ? 1 : 0);
    ASSERT_NEAR(nested_batch_visible.t[0], scalar_visible.t, 1e-9);
    ASSERT_EQ(nested_batch_visible.cell_ids[0], scalar_visible.cell_id);
    ASSERT_EQ(nested_batch_visible.material_ids[0], scalar_visible.material_id);
    ASSERT_EQ(nested_batch_visible.surface_ids[0], scalar_visible.surface_id);
    alea_ray_first_visible_batch_result_free(&nested_batch_visible);

    const alea_ray_batch_query_t nested_any_query = {
        .kind = ALEA_RAY_QUERY_ANY_HIT,
        .material_filter = -1,
        .t_mins = nested_t_mins,
        .t_maxs = nested_t_maxs
    };
    alea_ray_any_hit_batch_result_t nested_batch_any;
    alea_ray_any_hit_batch_result_init(&nested_batch_any);
    ASSERT_EQ(alea_raycast_hier_any_hit_batch_nocache(
                  sys, nested_origins, nested_directions, 1,
                  &nested_any_query, &nested_batch_any), 0);
    ASSERT_EQ(nested_batch_any.hits[0], scalar_visible.found ? 1 : 0);
    alea_ray_any_hit_batch_result_free(&nested_batch_any);

    const double packet_origins[] = {
        -2.0, -0.5, 0.0, -2.0, 0.5, 0.0, 4.0, -0.5, 0.0,
        4.0, 0.5, 0.0, -2.0, 3.0, 0.0
    };
    const double packet_directions[] = {
        1.0, 0.0, 0.0, 1.0, 0.0, 0.0, -1.0, 0.0, 0.0,
        -1.0, 0.0, 0.0, 1.0, 0.0, 0.0
    };
    const double packet_t_mins[] = {0.0, 0.5, 0.0, 0.5, 0.0};
    const double packet_t_maxs[] = {6.0, 6.0, 6.0, 6.0, 6.0};
    assert_hier_packet_batch_matches_scalar(
        sys, packet_origins, packet_directions, 5, packet_t_mins,
        packet_t_maxs, -1);

    alea_raycast_result_init(&nested_scratch);
    ASSERT_EQ(alea_raycast_hier_fast_segments(
                  sys, -2.0, -0.5, 0.0, 1.0, 0.0, 0.0, 6.0,
                  &nested_scratch), 0);
    ASSERT(nested_scratch.segments.count > 0);
    alea_raycast_result_free(&nested_scratch);

    alea_raycast_result_init(&nested_scratch);
    ASSERT_EQ(alea_raycast_hier_with_hits(
                  sys, -2.0, -0.5, 0.0, 1.0, 0.0, 0.0, 6.0,
                  &nested_scratch), 0);
    ASSERT(nested_scratch.segments.count > 0);
    alea_raycast_result_free(&nested_scratch);

    mcnp_model_destroy(model);
}

TEST(lattice_nested_transformed_packet_parity) {
    mcnp_model_t* model =
        mcnp_load("tests/data/mcnp_nested_lattice_transformed.mcnp");
    if (!model) SKIP("Test data file not found");
    alea_system_t* sys = model->sys;
    ASSERT_EQ(alea_prepare_query_acceleration(sys), 0);

    int cell_id, material;
    ASSERT_EQ(alea_find_cell_lazy(sys, 4.5, -0.5, 0.0,
                                  &cell_id, &material, NULL), 0);
    ASSERT_EQ(material, 1);
    const double origins[] = {
        3.0, -0.5, 0.0, 3.0, 0.5, 0.0, 9.0, -0.5, 0.0,
        9.0, 0.5, 0.0, 3.0, 3.0, 0.0
    };
    const double directions[] = {
        1.0, 0.0, 0.0, 1.0, 0.0, 0.0, -1.0, 0.0, 0.0,
        -1.0, 0.0, 0.0, 1.0, 0.0, 0.0
    };
    const double t_mins[] = {0.0, 0.5, 0.0, 0.5, 0.0};
    const double t_maxs[] = {8.0, 8.0, 8.0, 8.0, 8.0};
    assert_hier_packet_batch_matches_scalar(
        sys, origins, directions, 5, t_mins, t_maxs, -1);
    assert_hier_packet_batch_matches_scalar(
        sys, origins, directions, 5, t_mins, t_maxs, 3);

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

typedef struct {
    size_t count;
    size_t overlap_owner_count;
    int saw_occurrence_keys;
    int saw_truncated;
    int stop_after_first;
} coverage_stream_probe_t;

typedef struct {
    alea_ray_coverage_kind_t kinds[8];
    double enters[8];
    double exits[8];
    size_t count;
} coverage_kind_trace_t;

typedef struct {
    size_t row_indices[8];
    alea_ray_coverage_kind_t kinds[8];
    double enters[8];
    double exits[8];
    size_t count;
} coverage_row_trace_t;

static int probe_coverage_interval(
    void* context, const alea_ray_coverage_interval_t* interval) {
    coverage_stream_probe_t* probe = context;
    probe->count++;
    if (interval->kind == ALEA_RAY_COVERAGE_OVERLAP)
        probe->overlap_owner_count = interval->owner_count;
    if (interval->kind == ALEA_RAY_COVERAGE_TRUNCATED)
        probe->saw_truncated = 1;
    for (size_t i = 0; i < interval->owner_count; i++)
        if (interval->owners[i].occurrence_key != 0)
            probe->saw_occurrence_keys = 1;
    return probe->stop_after_first && probe->count == 1;
}

static int collect_coverage_kind(void* context,
                                 const alea_ray_coverage_interval_t* interval) {
    coverage_kind_trace_t* trace = context;
    if (trace->count >= 8) return -1;
    trace->kinds[trace->count] = interval->kind;
    trace->enters[trace->count] = interval->t_enter;
    trace->exits[trace->count] = interval->t_exit;
    trace->count++;
    return 0;
}

static int collect_coverage_row(void* context, size_t row_index,
                                const alea_ray_coverage_interval_t* interval) {
    coverage_row_trace_t* trace = context;
    if (trace->count >= 8) return -1;
    trace->row_indices[trace->count] = row_index;
    trace->kinds[trace->count] = interval->kind;
    trace->enters[trace->count] = interval->t_enter;
    trace->exits[trace->count] = interval->t_exit;
    trace->count++;
    return 0;
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

    /* The coverage oracle consumes global breakpoints only.  It must not
     * build deck-precedence segments as an implementation side effect. */
    ASSERT_EQ(alea_prepare_query_acceleration(sys), 0);
    alea_ray_t diagnostic_ray;
    ASSERT_EQ(alea_ray_init(&diagnostic_ray, -10, 0, 0, 1, 0, 0), 0);
    alea_raycast_result_t breakpoint_scratch;
    alea_raycast_result_init(&breakpoint_scratch);
    ASSERT_EQ(alea_raycast_global_breakpoints_reuse_nocache(
                  sys, &diagnostic_ray, 20.0, &breakpoint_scratch), 0);
    ASSERT(breakpoint_scratch.hits.count > 0);
    ASSERT_EQ(breakpoint_scratch.segments.count, 0);
    alea_raycast_result_free(&breakpoint_scratch);

    alea_ray_interval_finding_t f[16];
    int n = alea_ray_classify_intervals(sys, -10, 0, 0, 1, 0, 0, 20, f, 16);
    ASSERT_EQ(n, 7);

    alea_raycast_result_t coverage_scratch;
    alea_raycast_result_init(&coverage_scratch);
    coverage_stream_probe_t probe = {0};
    ASSERT_EQ(alea_ray_coverage_sweep_reuse_nocache(
                  sys, &diagnostic_ray, 20.0, &coverage_scratch,
                  probe_coverage_interval, &probe), n);
    ASSERT_EQ(probe.count, (size_t)n);
    ASSERT_EQ(probe.overlap_owner_count, 2);
    ASSERT(probe.saw_occurrence_keys);
    alea_raycast_result_free(&coverage_scratch);

    alea_raycast_result_init(&coverage_scratch);
    coverage_stream_probe_t stopped_probe = { .stop_after_first = 1 };
    ASSERT_EQ(alea_ray_coverage_sweep_reuse_nocache(
                  sys, &diagnostic_ray, 20.0, &coverage_scratch,
                  probe_coverage_interval, &stopped_probe), 1);
    ASSERT_EQ(stopped_probe.count, 1);
    alea_raycast_result_free(&coverage_scratch);

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

TEST(ray_coverage_reports_owner_budget_truncation) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);
    const int sphere = alea_sphere_surface(sys, 1, 0, 0, 0, 1.0);
    ASSERT(sphere >= 0);
    const alea_node_id_t inside = alea_halfspace(sys, sphere, -1);
    ASSERT_NE(inside, ALEA_NODE_ID_INVALID);
    const int material = alea_add_material(sys, 1);
    ASSERT(material >= 0);
    /* The coverage API retains 32 owners and reserves one sentinel slot to
     * detect saturation.  All cells intentionally claim the same volume. */
    for (int i = 0; i < 33; i++)
        ASSERT(alea_add_cell(sys, i + 1, inside, material, 1.0, 0) >= 0);
    ASSERT_EQ(alea_prepare_query_acceleration(sys), 0);

    alea_ray_t ray;
    ASSERT_EQ(alea_ray_init(&ray, -2, 0, 0, 1, 0, 0), 0);
    alea_raycast_result_t scratch;
    alea_raycast_result_init(&scratch);
    coverage_stream_probe_t probe = {0};
    ASSERT_EQ(alea_ray_coverage_sweep_reuse_nocache(
                  sys, &ray, 4.0, &scratch,
                  probe_coverage_interval, &probe), 3);
    ASSERT_EQ(probe.count, 3);
    ASSERT(probe.saw_truncated);
    alea_raycast_result_free(&scratch);

    /* The legacy interval adapter must preserve incomplete coverage rather
     * than deriving an overlap/unique answer from the retained owner prefix. */
    alea_ray_interval_finding_t findings[4];
    const int finding_count = alea_ray_classify_intervals(
        sys, -2, 0, 0, 1, 0, 0, 4.0, findings, 4);
    ASSERT_EQ(finding_count, 3);
    ASSERT_EQ(findings[1].kind, ALEA_INTERVAL_TRUNCATED);

    /* The executor must preserve owner saturation as published data, while
     * enforcing max_owners across all worker arenas rather than per worker. */
    alea_ray_coverage_row_t rows[2] = {
        { .ray = ray, .t_max = 4.0, .direction_tag = 1,
          .transverse_coordinate = 0.0 },
        { .ray = ray, .t_max = 4.0, .direction_tag = 1,
          .transverse_coordinate = 1.0 }
    };
    alea_ray_coverage_slice_limits_t limits;
    alea_ray_coverage_slice_limits_init(&limits);
    alea_ray_coverage_executor_t executor;
    alea_ray_coverage_executor_init(&executor);
    ASSERT_EQ(alea_ray_coverage_executor_prepare(&executor, 2), 0);
    alea_ray_coverage_slice_result_t slice;
    alea_ray_coverage_slice_result_init(&slice);
    ASSERT_EQ(alea_ray_coverage_slice_build_executor_nocache(
                  sys, rows, 2, &limits, &executor, &slice), 0);
    ASSERT_EQ(slice.row_count, (size_t)2);
    ASSERT_EQ(slice.interval_count, (size_t)6);
    ASSERT_EQ(slice.owner_count, (size_t)64);
    ASSERT_EQ(slice.kinds[1], (uint8_t)ALEA_RAY_COVERAGE_TRUNCATED);
    ASSERT_EQ(slice.kinds[4], (uint8_t)ALEA_RAY_COVERAGE_TRUNCATED);
    ASSERT_EQ(slice.owner_count_lower_bounds[1], (size_t)33);
    ASSERT_EQ(slice.owner_count_lower_bounds[4], (size_t)33);
    const size_t* prior_offsets = slice.row_offsets;
    limits.max_owners = 63;
    ASSERT_EQ(alea_ray_coverage_slice_build_executor_nocache(
                  sys, rows, 2, &limits, &executor, &slice), -1);
    ASSERT_EQ(slice.row_offsets, prior_offsets);
    ASSERT_EQ(slice.owner_count, (size_t)64);
    alea_ray_coverage_slice_result_free(&slice);
    alea_ray_coverage_executor_free(&executor);

    alea_destroy(sys);
}

TEST(ray_coverage_domain_partitions_exterior_from_interior_gap) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);
    const int sphere = alea_sphere_surface(sys, 1, 0, 0, 0, 1.0);
    ASSERT(sphere >= 0);
    const alea_node_id_t inside = alea_halfspace(sys, sphere, -1);
    ASSERT_NE(inside, ALEA_NODE_ID_INVALID);
    const int material = alea_add_material(sys, 1);
    ASSERT(material >= 0);
    ASSERT(alea_add_cell(sys, 1, inside, material, 1.0, 0) >= 0);
    ASSERT_EQ(alea_prepare_query_acceleration(sys), 0);

    alea_ray_t ray;
    ASSERT_EQ(alea_ray_init(&ray, -2, 0, 0, 1, 0, 0), 0);
    const alea_ray_coverage_domain_t domain = {
        .t_min = 0.5, .t_max = 3.5,
        .has_domain = 1, .report_allowed_exterior = 1
    };
    alea_raycast_result_t scratch;
    alea_raycast_result_init(&scratch);
    coverage_kind_trace_t trace = {0};
    ASSERT_EQ(alea_ray_coverage_sweep_domain_reuse_nocache(
                  sys, &ray, 4.0, &domain, &scratch,
                  collect_coverage_kind, &trace), 5);
    ASSERT_EQ(trace.count, (size_t)5);
    ASSERT_EQ(trace.kinds[0], ALEA_RAY_COVERAGE_ALLOWED_EXTERIOR);
    ASSERT_NEAR(trace.exits[0], 0.5, 1e-9);
    ASSERT_EQ(trace.kinds[1], ALEA_RAY_COVERAGE_GAP);
    ASSERT_NEAR(trace.enters[1], 0.5, 1e-9);
    ASSERT_NEAR(trace.exits[1], 1.0, 1e-9);
    ASSERT_EQ(trace.kinds[2], ALEA_RAY_COVERAGE_UNIQUE);
    ASSERT_EQ(trace.kinds[3], ALEA_RAY_COVERAGE_GAP);
    ASSERT_EQ(trace.kinds[4], ALEA_RAY_COVERAGE_ALLOWED_EXTERIOR);
    ASSERT_NEAR(trace.enters[4], 3.5, 1e-9);

    alea_raycast_result_free(&scratch);
    alea_destroy(sys);
}

TEST(ray_coverage_rows_stream_in_input_order) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);
    const int sphere = alea_sphere_surface(sys, 1, 0, 0, 0, 1.0);
    ASSERT(sphere >= 0);
    const alea_node_id_t inside = alea_halfspace(sys, sphere, -1);
    ASSERT_NE(inside, ALEA_NODE_ID_INVALID);
    const int material = alea_add_material(sys, 1);
    ASSERT(material >= 0);
    ASSERT(alea_add_cell(sys, 1, inside, material, 1.0, 0) >= 0);
    ASSERT_EQ(alea_prepare_query_acceleration(sys), 0);

    const alea_ray_coverage_domain_t domain = {
        .t_min = 0.5, .t_max = 3.5, .has_domain = 1
    };
    alea_ray_coverage_row_t rows[2] = {
        { .t_max = 4.0, .domain = domain, .use_domain = 1,
          .direction_tag = 0, .transverse_coordinate = 0.0 },
        { .t_max = 4.0, .domain = domain, .use_domain = 1,
          .direction_tag = 0, .transverse_coordinate = 2.0 }
    };
    ASSERT_EQ(alea_ray_init(&rows[0].ray, -2, 0, 0, 1, 0, 0), 0);
    ASSERT_EQ(alea_ray_init(&rows[1].ray, -2, 2, 0, 1, 0, 0), 0);

    alea_raycast_result_t scratch;
    alea_raycast_result_init(&scratch);
    alea_ray_coverage_executor_t executor;
    alea_ray_coverage_executor_init(&executor);
    ASSERT_EQ(alea_ray_coverage_executor_prepare(&executor, 2), 0);
    ASSERT_EQ(executor.worker_count, (size_t)2);
    ASSERT_NOT_NULL(executor.workers);
    ASSERT_EQ(alea_ray_coverage_executor_worker_for_row(&executor, 0),
              &executor.workers[0]);
    ASSERT_EQ(alea_ray_coverage_executor_worker_for_row(&executor, 3),
              &executor.workers[1]);
    ASSERT_EQ(alea_ray_coverage_executor_prepare(&executor, 1), 0);
    ASSERT_EQ(executor.worker_count, (size_t)1);
    alea_ray_coverage_executor_free(&executor);
    coverage_row_trace_t trace = {0};
    ASSERT_EQ(alea_ray_coverage_rows_serial_reuse_nocache(
                  sys, rows, 2, &scratch, collect_coverage_row, &trace), 0);
    ASSERT_EQ(trace.count, (size_t)4);
    ASSERT_EQ(trace.row_indices[0], (size_t)0);
    ASSERT_EQ(trace.row_indices[1], (size_t)0);
    ASSERT_EQ(trace.row_indices[2], (size_t)0);
    ASSERT_EQ(trace.row_indices[3], (size_t)1);
    ASSERT_EQ(trace.kinds[0], ALEA_RAY_COVERAGE_GAP);
    ASSERT_EQ(trace.kinds[1], ALEA_RAY_COVERAGE_UNIQUE);
    ASSERT_EQ(trace.kinds[2], ALEA_RAY_COVERAGE_GAP);
    ASSERT_EQ(trace.kinds[3], ALEA_RAY_COVERAGE_GAP);
    ASSERT_NEAR(trace.enters[3], 0.5, 1e-9);
    ASSERT_NEAR(trace.exits[3], 3.5, 1e-9);

    alea_raycast_result_free(&scratch);
    alea_destroy(sys);
}

TEST(ray_coverage_slice_builds_transactional_csr) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);
    const int sphere = alea_sphere_surface(sys, 1, 0, 0, 0, 1.0);
    ASSERT(sphere >= 0);
    const alea_node_id_t inside = alea_halfspace(sys, sphere, -1);
    ASSERT_NE(inside, ALEA_NODE_ID_INVALID);
    const int material = alea_add_material(sys, 1);
    ASSERT(material >= 0);
    ASSERT(alea_add_cell(sys, 1, inside, material, 1.0, 0) >= 0);
    ASSERT_EQ(alea_prepare_query_acceleration(sys), 0);

    const alea_ray_coverage_domain_t domain = {
        .t_min = 0.5, .t_max = 3.5, .has_domain = 1
    };
    alea_ray_coverage_row_t rows[2] = {
        { .t_max = 4.0, .domain = domain, .use_domain = 1,
          .direction_tag = 4, .transverse_coordinate = 0.0 },
        { .t_max = 4.0, .domain = domain, .use_domain = 1,
          .direction_tag = 4, .transverse_coordinate = 2.0 }
    };
    ASSERT_EQ(alea_ray_init(&rows[0].ray, -2, 0, 0, 1, 0, 0), 0);
    ASSERT_EQ(alea_ray_init(&rows[1].ray, -2, 2, 0, 1, 0, 0), 0);

    alea_raycast_result_t scratch;
    alea_raycast_result_init(&scratch);
    alea_ray_coverage_slice_result_t result;
    alea_ray_coverage_slice_result_init(&result);
    alea_ray_coverage_slice_limits_t limits;
    alea_ray_coverage_slice_limits_init(&limits);
    ASSERT_EQ(alea_ray_coverage_slice_build_serial_nocache(
                  sys, rows, 2, &limits, &scratch, &result), 0);
    ASSERT_EQ(result.row_count, (size_t)2);
    ASSERT_EQ(result.interval_count, (size_t)4);
    ASSERT_EQ(result.owner_count, (size_t)1);
    ASSERT_EQ(result.row_offsets[0], (size_t)0);
    ASSERT_EQ(result.row_offsets[1], (size_t)3);
    ASSERT_EQ(result.row_offsets[2], (size_t)4);
    ASSERT_EQ(result.row_direction_tags[0], (uint8_t)4);
    ASSERT_NEAR(result.row_transverse_coordinates[1], 2.0, 1e-12);
    ASSERT_EQ(result.kinds[0], (uint8_t)ALEA_RAY_COVERAGE_GAP);
    ASSERT_EQ(result.kinds[1], (uint8_t)ALEA_RAY_COVERAGE_UNIQUE);
    ASSERT_EQ(result.kinds[3], (uint8_t)ALEA_RAY_COVERAGE_GAP);
    ASSERT_NEAR(result.t_enter[1], 1.0, 1e-9);
    ASSERT_NEAR(result.t_exit[1], 3.0, 1e-9);
    ASSERT_EQ(result.owner_offsets[0], (size_t)0);
    ASSERT_EQ(result.owner_offsets[1], (size_t)0);
    ASSERT_EQ(result.owner_offsets[2], (size_t)1);
    ASSERT_EQ(result.owner_offsets[4], (size_t)1);
    ASSERT_EQ(result.owner_count_lower_bounds[1], (size_t)1);
    ASSERT_EQ(result.owner_cell_ids[0], 1);
    ASSERT_EQ(result.owner_parent_occurrence_keys[0], (uint64_t)0);
    ASSERT_EQ(alea_ray_coverage_slice_rows_same_signature(&result, 0, 0), 1);
    ASSERT_EQ(alea_ray_coverage_slice_rows_same_signature(&result, 0, 1), 0);
    uint8_t refine_between[1] = {0};
    ASSERT_EQ(alea_ray_coverage_slice_mark_refinement_boundaries(
                  &result, refine_between), 1);
    ASSERT_EQ(refine_between[0], (uint8_t)1);

    /* Phase 10 must change only scheduling and staging: the worker-arena
     * executor compacts back to exactly the serial CSR order. */
    alea_ray_coverage_executor_t executor;
    alea_ray_coverage_executor_init(&executor);
    /* More workers than rows exercises deterministic empty worker arenas. */
    ASSERT_EQ(alea_ray_coverage_executor_prepare(&executor, 3), 0);
    alea_ray_coverage_slice_result_t executor_result;
    alea_ray_coverage_slice_result_init(&executor_result);
    ASSERT_EQ(alea_ray_coverage_slice_build_executor_nocache(
                  sys, rows, 2, &limits, &executor, &executor_result), 0);
    ASSERT_EQ(executor_result.row_count, result.row_count);
    ASSERT_EQ(executor_result.interval_count, result.interval_count);
    ASSERT_EQ(executor_result.owner_count, result.owner_count);
    ASSERT_EQ(memcmp(executor_result.row_offsets, result.row_offsets,
                     3 * sizeof(*result.row_offsets)), 0);
    ASSERT_EQ(memcmp(executor_result.t_enter, result.t_enter,
                     4 * sizeof(*result.t_enter)), 0);
    ASSERT_EQ(memcmp(executor_result.t_exit, result.t_exit,
                     4 * sizeof(*result.t_exit)), 0);
    ASSERT_EQ(memcmp(executor_result.kinds, result.kinds,
                     4 * sizeof(*result.kinds)), 0);
    ASSERT_EQ(memcmp(executor_result.owner_offsets, result.owner_offsets,
                     5 * sizeof(*result.owner_offsets)), 0);
    ASSERT_EQ(memcmp(executor_result.owner_occurrence_keys,
                     result.owner_occurrence_keys,
                     sizeof(*result.owner_occurrence_keys)), 0);
    const size_t* executor_previous_offsets = executor_result.row_offsets;
    limits.max_intervals = 3;
    ASSERT_EQ(alea_ray_coverage_slice_build_executor_nocache(
                  sys, rows, 2, &limits, &executor, &executor_result), -1);
    ASSERT_EQ(executor_result.row_offsets, executor_previous_offsets);
    ASSERT_EQ(executor_result.interval_count, (size_t)4);
    limits.max_intervals = 0;
    /* A malformed row owned by a different worker aborts the complete
     * operation; already-staged rows never leak into publication. */
    alea_ray_coverage_row_t invalid_rows[2] = { rows[0], rows[1] };
    invalid_rows[1].t_max = 0.0;
    ASSERT_EQ(alea_ray_coverage_slice_build_executor_nocache(
                  sys, invalid_rows, 2, &limits, &executor,
                  &executor_result), -1);
    ASSERT_EQ(executor_result.row_offsets, executor_previous_offsets);
    ASSERT_EQ(executor_result.interval_count, (size_t)4);
    /* Interruption is an operation failure too and must retain the last
     * publication without depending on a scheduling race. */
    alea_interrupt();
    ASSERT_EQ(alea_ray_coverage_slice_build_executor_nocache(
                  sys, rows, 2, &limits, &executor, &executor_result), -1);
    ASSERT_EQ(executor_result.row_offsets, executor_previous_offsets);
    ASSERT_EQ(executor_result.interval_count, (size_t)4);
    alea_clear_interrupt();
    alea_ray_coverage_slice_result_free(&executor_result);
    alea_ray_coverage_executor_free(&executor);

    alea_ray_coverage_row_t refined_rows[3] = {0};
    size_t refined_count = 0;
    ASSERT_EQ(alea_ray_coverage_rows_refine_midpoints(
                  rows, 2, refine_between, 3, refined_rows, 3,
                  &refined_count), 0);
    ASSERT_EQ(refined_count, (size_t)3);
    ASSERT_NEAR(refined_rows[0].transverse_coordinate, 0.0, 1e-12);
    ASSERT_NEAR(refined_rows[1].transverse_coordinate, 1.0, 1e-12);
    ASSERT_NEAR(refined_rows[2].transverse_coordinate, 2.0, 1e-12);
    ASSERT_NEAR(refined_rows[1].ray.oy, 1.0, 1e-12);
    const alea_ray_coverage_row_t untouched_row = {
        .direction_tag = 99, .transverse_coordinate = -1.0
    };
    alea_ray_coverage_row_t failed_refinement[3] = { untouched_row };
    size_t failed_count = 123;
    ASSERT_EQ(alea_ray_coverage_rows_refine_midpoints(
                  rows, 2, refine_between, 2, failed_refinement, 3,
                  &failed_count), -1);
    ASSERT_EQ(failed_refinement[0].direction_tag, (uint8_t)99);
    ASSERT_EQ(failed_count, (size_t)123);

    alea_ray_coverage_slice_result_t adaptive;
    alea_ray_coverage_slice_result_init(&adaptive);
    ASSERT_EQ(alea_ray_coverage_slice_build_adaptive_serial_nocache(
                  sys, rows, 2, 0, NULL, &scratch, &adaptive), 0);
    ASSERT_EQ(adaptive.row_count, (size_t)2);
    ASSERT_EQ(adaptive.refinement_status,
              ALEA_RAY_COVERAGE_REFINEMENT_MAX_DEPTH);
    alea_ray_coverage_slice_result_free(&adaptive);
    limits = (alea_ray_coverage_slice_limits_t){ .max_rows = 2 };
    ASSERT_EQ(alea_ray_coverage_slice_build_adaptive_serial_nocache(
                  sys, rows, 2, 2, &limits, &scratch, &adaptive), 0);
    ASSERT_EQ(adaptive.row_count, (size_t)2);
    ASSERT_EQ(adaptive.refinement_status,
              ALEA_RAY_COVERAGE_REFINEMENT_MAX_ROWS);
    alea_ray_coverage_slice_result_free(&adaptive);
    ASSERT_EQ(alea_ray_coverage_slice_build_adaptive_serial_nocache(
                  sys, rows, 2, 1, NULL, &scratch, &adaptive), 0);
    ASSERT_EQ(adaptive.row_count, (size_t)3);
    ASSERT_EQ(adaptive.refinement_status,
              ALEA_RAY_COVERAGE_REFINEMENT_MAX_DEPTH);
    alea_ray_coverage_slice_result_free(&adaptive);

    alea_ray_coverage_slice_limits_init(&limits);
    alea_ray_coverage_executor_init(&executor);
    ASSERT_EQ(alea_ray_coverage_executor_prepare(&executor, 2), 0);
    ASSERT_EQ(alea_ray_coverage_slice_build_adaptive_policy_executor_nocache(
                  sys, rows, 2, 1, NULL, &limits, &executor, &adaptive), 0);
    alea_ray_coverage_slice_result_t serial_adaptive;
    alea_ray_coverage_slice_result_init(&serial_adaptive);
    ASSERT_EQ(alea_ray_coverage_slice_build_adaptive_serial_nocache(
                  sys, rows, 2, 1, &limits, &scratch, &serial_adaptive), 0);
    ASSERT_EQ(adaptive.row_count, (size_t)3);
    ASSERT_EQ(adaptive.row_count, serial_adaptive.row_count);
    ASSERT_EQ(adaptive.interval_count, serial_adaptive.interval_count);
    ASSERT_EQ(adaptive.owner_count, serial_adaptive.owner_count);
    ASSERT_EQ(adaptive.refinement_status,
              ALEA_RAY_COVERAGE_REFINEMENT_MAX_DEPTH);
    ASSERT_EQ(adaptive.refinement_status, serial_adaptive.refinement_status);
    ASSERT_NEAR(adaptive.row_transverse_coordinates[1], 1.0, 1e-12);
    ASSERT_EQ(memcmp(adaptive.row_offsets, serial_adaptive.row_offsets,
                     (adaptive.row_count + 1) * sizeof(*adaptive.row_offsets)), 0);
    ASSERT_EQ(memcmp(adaptive.t_enter, serial_adaptive.t_enter,
                     adaptive.interval_count * sizeof(*adaptive.t_enter)), 0);
    ASSERT_EQ(memcmp(adaptive.t_exit, serial_adaptive.t_exit,
                     adaptive.interval_count * sizeof(*adaptive.t_exit)), 0);
    ASSERT_EQ(memcmp(adaptive.kinds, serial_adaptive.kinds,
                     adaptive.interval_count * sizeof(*adaptive.kinds)), 0);
    ASSERT_EQ(memcmp(adaptive.owner_offsets, serial_adaptive.owner_offsets,
                     (adaptive.interval_count + 1) *
                     sizeof(*adaptive.owner_offsets)), 0);
    ASSERT_EQ(memcmp(adaptive.owner_occurrence_keys,
                     serial_adaptive.owner_occurrence_keys,
                     adaptive.owner_count *
                     sizeof(*adaptive.owner_occurrence_keys)), 0);
    alea_ray_coverage_slice_result_free(&serial_adaptive);
    alea_ray_coverage_slice_result_free(&adaptive);
    alea_ray_coverage_executor_free(&executor);

    const size_t* previous_offsets = result.row_offsets;
    limits.max_intervals = 3;
    ASSERT_EQ(alea_ray_coverage_slice_build_serial_nocache(
                  sys, rows, 2, &limits, &scratch, &result), -1);
    ASSERT_EQ(result.row_offsets, previous_offsets);
    ASSERT_EQ(result.interval_count, (size_t)4);
    ASSERT_EQ(result.owner_count, (size_t)1);

    /* The byte limit measures the published CSR arrays exactly, rather than
     * allocator capacity.  One byte below this known result must also retain
     * the prior publication. */
    const size_t published_bytes =
        3 * sizeof(size_t) + 2 * sizeof(uint8_t) + 2 * sizeof(double) +
        4 * (2 * sizeof(double) + sizeof(uint8_t) + sizeof(size_t)) +
        5 * sizeof(size_t) +
        6 * sizeof(int) + 2 * sizeof(uint64_t) + sizeof(uint8_t);
    limits = (alea_ray_coverage_slice_limits_t){
        .max_bytes = published_bytes - 1
    };
    ASSERT_EQ(alea_ray_coverage_slice_build_serial_nocache(
                  sys, rows, 2, &limits, &scratch, &result), -1);
    ASSERT_EQ(result.row_offsets, previous_offsets);
    ASSERT_EQ(result.interval_count, (size_t)4);

    /* The executor reserves this limit across all worker arenas, rather than
     * giving every worker an independent copy of the operation budget. */
    alea_ray_coverage_executor_init(&executor);
    ASSERT_EQ(alea_ray_coverage_executor_prepare(&executor, 3), 0);
    alea_ray_coverage_slice_result_init(&executor_result);
    limits.max_bytes = published_bytes;
    ASSERT_EQ(alea_ray_coverage_slice_build_executor_nocache(
                  sys, rows, 2, &limits, &executor, &executor_result), 0);
    ASSERT_EQ(executor_result.interval_count, (size_t)4);
    const size_t* executor_byte_limited_offsets = executor_result.row_offsets;
    limits.max_bytes = published_bytes - 1;
    ASSERT_EQ(alea_ray_coverage_slice_build_executor_nocache(
                  sys, rows, 2, &limits, &executor, &executor_result), -1);
    ASSERT_EQ(executor_result.row_offsets, executor_byte_limited_offsets);
    ASSERT_EQ(executor_result.interval_count, (size_t)4);
    alea_ray_coverage_slice_result_free(&executor_result);
    alea_ray_coverage_executor_free(&executor);

    alea_ray_coverage_slice_result_free(&result);
    alea_raycast_result_free(&scratch);
    alea_destroy(sys);
}

/* Signature equality deliberately ignores endpoints, so two chords through one
 * sphere look identical to it.  The displacement, density, and finding signals
 * exist precisely to keep such a pair refinable. */
TEST(ray_coverage_refinement_signals_select_matching_rows) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);
    const int sphere = alea_sphere_surface(sys, 1, 0, 0, 0, 1.0);
    ASSERT(sphere >= 0);
    const alea_node_id_t inside = alea_halfspace(sys, sphere, -1);
    ASSERT_NE(inside, ALEA_NODE_ID_INVALID);
    const int material = alea_add_material(sys, 1);
    ASSERT(material >= 0);
    ASSERT(alea_add_cell(sys, 1, inside, material, 1.0, 0) >= 0);
    ASSERT_EQ(alea_prepare_query_acceleration(sys), 0);

    const alea_ray_coverage_domain_t domain = {
        .t_min = 0.5, .t_max = 3.5, .has_domain = 1
    };
    alea_ray_coverage_row_t rows[2] = {
        { .t_max = 4.0, .domain = domain, .use_domain = 1,
          .direction_tag = 0, .transverse_coordinate = 0.0 },
        { .t_max = 4.0, .domain = domain, .use_domain = 1,
          .direction_tag = 0, .transverse_coordinate = 0.5 }
    };
    ASSERT_EQ(alea_ray_init(&rows[0].ray, -2, 0.0, 0, 1, 0, 0), 0);
    ASSERT_EQ(alea_ray_init(&rows[1].ray, -2, 0.5, 0, 1, 0, 0), 0);

    alea_raycast_result_t scratch;
    alea_raycast_result_init(&scratch);
    alea_ray_coverage_slice_result_t result;
    alea_ray_coverage_slice_result_init(&result);
    alea_ray_coverage_slice_limits_t limits;
    alea_ray_coverage_slice_limits_init(&limits);
    ASSERT_EQ(alea_ray_coverage_slice_build_serial_nocache(
                  sys, rows, 2, &limits, &scratch, &result), 0);
    /* Both rows: GAP, UNIQUE, GAP -- identical owners, displaced boundaries. */
    ASSERT_EQ(result.row_count, (size_t)2);
    ASSERT_EQ(result.interval_count, (size_t)6);
    ASSERT_EQ(alea_ray_coverage_slice_rows_same_signature(&result, 0, 1), 1);
    ASSERT_NEAR(result.t_enter[1], 1.0, 1e-9);
    ASSERT_NEAR(result.t_enter[4], 2.0 - sqrt(0.75), 1e-9);

    uint8_t refine_between[1] = {0};
    size_t spacing_limited = 123;
    alea_ray_coverage_refinement_policy_t policy;
    alea_ray_coverage_refinement_policy_init(&policy);
    ASSERT_EQ(policy.signals, (uint32_t)ALEA_RAY_COVERAGE_REFINE_SIGNATURE);
    ASSERT_EQ(alea_ray_coverage_slice_mark_refinement_boundaries_policy(
                  &result, &policy, refine_between, &spacing_limited), 0);
    ASSERT_EQ(refine_between[0], (uint8_t)0);
    ASSERT_EQ(spacing_limited, (size_t)0);

    /* Endpoints move by about 0.134; a tolerance either side of that decides. */
    policy.signals = ALEA_RAY_COVERAGE_REFINE_DISPLACEMENT;
    policy.endpoint_displacement = 0.01;
    ASSERT_EQ(alea_ray_coverage_slice_mark_refinement_boundaries_policy(
                  &result, &policy, refine_between, NULL), 1);
    ASSERT_EQ(refine_between[0], (uint8_t)1);
    policy.endpoint_displacement = 0.5;
    ASSERT_EQ(alea_ray_coverage_slice_mark_refinement_boundaries_policy(
                  &result, &policy, refine_between, NULL), 0);
    ASSERT_EQ(refine_between[0], (uint8_t)0);

    policy.signals = ALEA_RAY_COVERAGE_REFINE_DENSITY;
    policy.crossing_density = 3;
    ASSERT_EQ(alea_ray_coverage_slice_mark_refinement_boundaries_policy(
                  &result, &policy, refine_between, NULL), 1);
    policy.crossing_density = 4;
    ASSERT_EQ(alea_ray_coverage_slice_mark_refinement_boundaries_policy(
                  &result, &policy, refine_between, NULL), 0);

    /* Both rows carry GAP intervals, so the finding signal selects the pair. */
    policy.signals = ALEA_RAY_COVERAGE_REFINE_FINDING;
    ASSERT_EQ(alea_ray_coverage_slice_mark_refinement_boundaries_policy(
                  &result, &policy, refine_between, NULL), 1);

    /* Splitting a 0.5 gap yields 0.25 spacing, below the configured minimum. */
    policy.min_transverse_spacing = 0.4;
    ASSERT_EQ(alea_ray_coverage_slice_mark_refinement_boundaries_policy(
                  &result, &policy, refine_between, &spacing_limited), 0);
    ASSERT_EQ(refine_between[0], (uint8_t)0);
    ASSERT_EQ(spacing_limited, (size_t)1);
    policy.min_transverse_spacing = 0.2;
    ASSERT_EQ(alea_ray_coverage_slice_mark_refinement_boundaries_policy(
                  &result, &policy, refine_between, &spacing_limited), 1);
    ASSERT_EQ(spacing_limited, (size_t)0);

    /* A selected signal without its tolerance is a policy error, not a
     * silently disabled probe. */
    policy = (alea_ray_coverage_refinement_policy_t){
        .signals = ALEA_RAY_COVERAGE_REFINE_DISPLACEMENT
    };
    ASSERT_EQ(alea_ray_coverage_slice_mark_refinement_boundaries_policy(
                  &result, &policy, refine_between, NULL), -1);
    policy = (alea_ray_coverage_refinement_policy_t){
        .signals = ALEA_RAY_COVERAGE_REFINE_DENSITY
    };
    ASSERT_EQ(alea_ray_coverage_slice_mark_refinement_boundaries_policy(
                  &result, &policy, refine_between, NULL), -1);
    policy = (alea_ray_coverage_refinement_policy_t){ .signals = 1u << 16 };
    ASSERT_EQ(alea_ray_coverage_slice_mark_refinement_boundaries_policy(
                  &result, &policy, refine_between, NULL), -1);

    /* A spacing-limited wave is a successful, explicitly limited result. */
    alea_ray_coverage_slice_result_t adaptive;
    alea_ray_coverage_slice_result_init(&adaptive);
    policy = (alea_ray_coverage_refinement_policy_t){
        .signals = ALEA_RAY_COVERAGE_REFINE_FINDING,
        .min_transverse_spacing = 0.4
    };
    ASSERT_EQ(alea_ray_coverage_slice_build_adaptive_policy_serial_nocache(
                  sys, rows, 2, 4, &policy, NULL, &scratch, &adaptive), 0);
    ASSERT_EQ(adaptive.row_count, (size_t)2);
    ASSERT_EQ(adaptive.refinement_status,
              ALEA_RAY_COVERAGE_REFINEMENT_MIN_SPACING);

    /* Relaxing the spacing lets the same signal refine until it is reached. */
    policy.min_transverse_spacing = 0.1;
    ASSERT_EQ(alea_ray_coverage_slice_build_adaptive_policy_serial_nocache(
                  sys, rows, 2, 4, &policy, NULL, &scratch, &adaptive), 0);
    ASSERT_EQ(adaptive.row_count, (size_t)5);
    ASSERT_EQ(adaptive.refinement_status,
              ALEA_RAY_COVERAGE_REFINEMENT_MIN_SPACING);
    ASSERT_NEAR(adaptive.row_transverse_coordinates[1], 0.125, 1e-12);
    alea_ray_coverage_slice_result_free(&adaptive);

    alea_ray_coverage_slice_result_free(&result);
    alea_raycast_result_free(&scratch);
    alea_destroy(sys);
}

TEST(ray_query_lowering_declares_semantics) {
    const alea_ray_query_t visible = {
        .kind = ALEA_RAY_QUERY_FIRST_VISIBLE,
        .backend = ALEA_RAY_QUERY_BACKEND_AUTO,
        .fields = ALEA_RAY_QUERY_FIELD_SURFACE_NORMAL,
        .t_min = 0.25,
        .t_max = 8.0,
        .material_filter = -1
    };
    alea_ray_plan_t plan;
    ASSERT_EQ(alea_ray_query_lower(&visible, &plan), 0);
    ASSERT_EQ(plan.engine, ALEA_RAY_ENGINE_SELECTED_WALKER);
    ASSERT_EQ(plan.ownership, ALEA_RAY_OWNERSHIP_TRACK_COHERENT);
    ASSERT(plan.requirements.need_selected_owner);
    ASSERT(plan.requirements.need_surface_identity);
    ASSERT(plan.requirements.need_normal);

    const alea_ray_query_t segments = {
        .kind = ALEA_RAY_QUERY_SEGMENTS,
        .backend = ALEA_RAY_QUERY_BACKEND_AUTO,
        .t_max = 8.0,
        .material_filter = -1
    };
    ASSERT_EQ(alea_ray_query_lower(&segments, &plan), 0);
    ASSERT_EQ(plan.engine, ALEA_RAY_ENGINE_SELECTED_WALKER);
    ASSERT_EQ(plan.ownership, ALEA_RAY_OWNERSHIP_TRACK_COHERENT);
    ASSERT(plan.requirements.need_selected_owner);

    const alea_ray_query_t boundaries = {
        .kind = ALEA_RAY_QUERY_BOUNDARY_EVENTS,
        .backend = ALEA_RAY_QUERY_BACKEND_AUTO,
        .t_max = 8.0,
        .material_filter = -1
    };
    ASSERT_EQ(alea_ray_query_lower(&boundaries, &plan), 0);
    ASSERT_EQ(plan.engine, ALEA_RAY_ENGINE_GLOBAL_BREAKPOINTS);
    ASSERT_EQ(plan.ownership, ALEA_RAY_OWNERSHIP_SELECT_CANONICAL);
    ASSERT(plan.requirements.need_all_coincident_primitives);

    alea_ray_query_t fast_boundaries = boundaries;
    fast_boundaries.backend = ALEA_RAY_QUERY_BACKEND_FAST_FORWARD;
    ASSERT_EQ(alea_ray_query_lower(&fast_boundaries, &plan), 0);
    ASSERT_EQ(plan.engine, ALEA_RAY_ENGINE_SELECTED_WALKER);
    ASSERT_EQ(plan.ownership, ALEA_RAY_OWNERSHIP_TRACK_COHERENT);
    ASSERT(!plan.requirements.need_all_coincident_primitives);

    const alea_ray_query_t invalid = {
        .kind = ALEA_RAY_QUERY_ANY_HIT,
        .backend = ALEA_RAY_QUERY_BACKEND_AUTO,
        .t_min = 2.0,
        .t_max = 1.0
    };
    ASSERT_EQ(alea_ray_query_lower(&invalid, &plan), -1);
}

TEST(coverage_occurrence_keys_distinguish_lattice_elements) {
    mcnp_model_t* model = mcnp_load("tests/data/mcnp_lattice_eval.mcnp");
    ASSERT_NOT_NULL(model);
    ASSERT_EQ(alea_prepare_query_acceleration(model->sys), 0);

    alea_cell_hit_t origin_hits[8], repeated_hits[8];
    uint64_t origin_keys[8], repeated_keys[8];
    const int origin_count = alea_find_all_cells_at_point_coverage_recursive(
        model->sys, 0.0, 0.0, 0.0, origin_hits, origin_keys, 8);
    const int repeated_count = alea_find_all_cells_at_point_coverage_recursive(
        model->sys, 4.0, 0.0, 0.0, repeated_hits, repeated_keys, 8);
    ASSERT(origin_count > 1);
    ASSERT_EQ(origin_count, repeated_count);

    int origin_child = -1, repeated_child = -1;
    for (int i = 0; i < origin_count; i++)
        if (origin_hits[i].depth > 0) { origin_child = i; break; }
    for (int i = 0; i < repeated_count; i++)
        if (repeated_hits[i].depth > 0) { repeated_child = i; break; }
    ASSERT(origin_child >= 0);
    ASSERT(repeated_child >= 0);
    ASSERT_EQ(origin_hits[origin_child].cell_index,
              repeated_hits[repeated_child].cell_index);
    ASSERT(origin_keys[origin_child] != repeated_keys[repeated_child]);
    mcnp_model_destroy(model);
}

TEST_MAIN()
