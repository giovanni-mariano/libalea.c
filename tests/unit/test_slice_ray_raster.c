// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

#include "alea_test.h"
#include "alea.h"
#include "alea_raycast.h"
#include "alea_slice.h"
#include "core/alea_system.h"

#include <stdlib.h>
#include <string.h>

#ifdef _OPENMP
#include <omp.h>
#endif

enum { RASTER_NU = 16, RASTER_NV = 3, RASTER_PIXELS = RASTER_NU * RASTER_NV };

static alea_system_t* raster_test_system(void) {
    alea_system_t* sys = alea_create();
    if (!sys) return NULL;
    int sphere = alea_sphere_surface(sys, 1, 0.0, 0.0, 0.0, 1.0);
    int material = alea_add_material(sys, 7);
    if (sphere < 0 || material < 0 ||
        alea_add_cell(sys, 10, alea_surface_at(sys, sphere)->neg_node,
                      material, -2.0, 0) < 0) {
        alea_destroy(sys);
        return NULL;
    }
    return sys;
}

static void raster_test_view(alea_slice_view_t* view) {
    alea_slice_view_axis(view, 2, 0.0, -2.0, 2.0, -0.5, 0.5);
}

static void raster_init_output(alea_slice_raster_t* output,
                               int32_t* cells, int32_t* materials,
                               double* densities, uint8_t* flags) {
    alea_slice_raster_init(output);
    output->nu = RASTER_NU;
    output->nv = RASTER_NV;
    output->fields = ALEA_SLICE_RASTER_CELL_ID |
                     ALEA_SLICE_RASTER_MATERIAL_ID |
                     ALEA_SLICE_RASTER_DENSITY |
                     ALEA_SLICE_RASTER_RESOLUTION_FLAGS;
    output->cell_ids = cells;
    output->material_ids = materials;
    output->densities = densities;
    output->resolution_flags = flags;
}

TEST(ray_slice_raster_compact_pixel_centers_and_defaults) {
    alea_system_t* sys = raster_test_system();
    ASSERT_NOT_NULL(sys);
    alea_slice_view_t view;
    raster_test_view(&view);

    alea_raycast_batch_options_t options = {
        .struct_size = sizeof(options),
        .fields = ALEA_RAY_BATCH_MATERIAL | ALEA_RAY_BATCH_DENSITY |
                  ALEA_RAY_BATCH_RESOLUTION_FLAGS,
        .projected_depth = -1
    };
    alea_raycast_batch_result_t* compact = alea_raycast_batch_result_create();
    ASSERT_NOT_NULL(compact);
    ASSERT_EQ(alea_trace_ray_slice_compact(sys, &view, RASTER_NV,
                                            &options, compact), 0);

    int32_t cells[RASTER_PIXELS], materials[RASTER_PIXELS];
    double densities[RASTER_PIXELS];
    uint8_t flags[RASTER_PIXELS];
    alea_slice_raster_t output;
    raster_init_output(&output, cells, materials, densities, flags);
    ASSERT_EQ(alea_rasterize_ray_slice_compact(&view, compact, &output), 0);

    for (size_t row = 0; row < RASTER_NV; row++) {
        for (size_t column = 0; column < RASTER_NU; column++) {
            const size_t pixel = row * RASTER_NU + column;
            if (column >= 4 && column < 12) {
                ASSERT_EQ(cells[pixel], 10);
                ASSERT_NE(materials[pixel], 0);
                ASSERT_NEAR(densities[pixel], 2.0, 1e-12);
            } else {
                ASSERT_EQ(cells[pixel], -1);
                ASSERT_EQ(materials[pixel], 0);
                ASSERT_NEAR(densities[pixel], 0.0, 1e-12);
                ASSERT_EQ(flags[pixel], 0);
            }
        }
    }
    alea_raycast_batch_result_destroy(compact);
    alea_destroy(sys);
}

TEST(ray_slice_raster_fused_matches_compact) {
    alea_system_t* sys = raster_test_system();
    ASSERT_NOT_NULL(sys);
    alea_slice_view_t view;
    raster_test_view(&view);

    int32_t compact_cells[RASTER_PIXELS], compact_materials[RASTER_PIXELS];
    int32_t fused_cells[RASTER_PIXELS], fused_materials[RASTER_PIXELS];
    int32_t compact_universes[RASTER_PIXELS], compact_fills[RASTER_PIXELS];
    int32_t fused_universes[RASTER_PIXELS], fused_fills[RASTER_PIXELS];
    double compact_densities[RASTER_PIXELS], fused_densities[RASTER_PIXELS];
    uint8_t compact_flags[RASTER_PIXELS], fused_flags[RASTER_PIXELS];
    alea_slice_raster_t compact_output, fused_output;
    raster_init_output(&compact_output, compact_cells, compact_materials,
                       compact_densities, compact_flags);
    raster_init_output(&fused_output, fused_cells, fused_materials,
                       fused_densities, fused_flags);
    compact_output.fields |= ALEA_SLICE_RASTER_UNIVERSE_ID |
                             ALEA_SLICE_RASTER_FILL_UNIVERSE;
    compact_output.universe_ids = compact_universes;
    compact_output.fill_universe_ids = compact_fills;
    fused_output.fields |= ALEA_SLICE_RASTER_UNIVERSE_ID |
                           ALEA_SLICE_RASTER_FILL_UNIVERSE;
    fused_output.universe_ids = fused_universes;
    fused_output.fill_universe_ids = fused_fills;

    alea_raycast_batch_options_t batch_options = {
        .struct_size = sizeof(batch_options),
        .fields = ALEA_RAY_BATCH_MATERIAL | ALEA_RAY_BATCH_DENSITY |
                  ALEA_RAY_BATCH_RESOLUTION_FLAGS |
                  ALEA_RAY_BATCH_PROJECTED_OWNER,
        .projected_depth = -1
    };
    alea_raycast_batch_result_t* compact = alea_raycast_batch_result_create();
    ASSERT_NOT_NULL(compact);
    ASSERT_EQ(alea_trace_ray_slice_compact(sys, &view, RASTER_NV,
                                            &batch_options, compact), 0);
    ASSERT_EQ(alea_rasterize_ray_slice_compact(&view, compact, &compact_output), 0);

    alea_slice_raster_options_t options;
    alea_slice_raster_options_init(&options);
    ASSERT_EQ(alea_trace_ray_slice_raster(sys, &view, &options, &fused_output), 0);
    ASSERT_EQ(memcmp(compact_cells, fused_cells, sizeof(compact_cells)), 0);
    ASSERT_EQ(memcmp(compact_materials, fused_materials, sizeof(compact_materials)), 0);
    ASSERT_EQ(memcmp(compact_universes, fused_universes, sizeof(compact_universes)), 0);
    ASSERT_EQ(memcmp(compact_fills, fused_fills, sizeof(compact_fills)), 0);
    ASSERT_EQ(memcmp(compact_densities, fused_densities, sizeof(compact_densities)), 0);
    ASSERT_EQ(memcmp(compact_flags, fused_flags, sizeof(compact_flags)), 0);

    alea_raycast_batch_result_destroy(compact);
    alea_destroy(sys);
}

TEST(ray_slice_raster_rejects_generic_batch_without_writing) {
    alea_system_t* sys = raster_test_system();
    ASSERT_NOT_NULL(sys);
    alea_slice_view_t view;
    raster_test_view(&view);
    const double origins[] = {-2.0, 0.0, 0.0};
    const double directions[] = {1.0, 0.0, 0.0};
    alea_raycast_batch_result_t* generic = alea_raycast_batch_result_create();
    ASSERT_NOT_NULL(generic);
    ASSERT_EQ(alea_raycast_hier_batch(sys, origins, directions, 1, 4.0,
                                      NULL, generic), 0);

    int32_t cells[RASTER_PIXELS];
    for (size_t i = 0; i < RASTER_PIXELS; i++) cells[i] = 12345;
    alea_slice_raster_t output;
    alea_slice_raster_init(&output);
    output.nu = RASTER_NU;
    output.nv = RASTER_NV;
    output.fields = ALEA_SLICE_RASTER_CELL_ID;
    output.cell_ids = cells;
    ASSERT_EQ(alea_rasterize_ray_slice_compact(&view, generic, &output), -1);
    for (size_t i = 0; i < RASTER_PIXELS; i++) ASSERT_EQ(cells[i], 12345);

    alea_raycast_batch_result_destroy(generic);
    alea_destroy(sys);
}

TEST(ray_slice_raster_parallel_rows_are_deterministic) {
    enum { nu = 256, nv = 256 };
    const size_t pixels = (size_t)nu * nv;
    alea_system_t* sys = raster_test_system();
    ASSERT_NOT_NULL(sys);
    alea_slice_view_t view;
    raster_test_view(&view);
    alea_raycast_batch_result_t* compact = alea_raycast_batch_result_create();
    ASSERT_NOT_NULL(compact);
    ASSERT_EQ(alea_trace_ray_slice_compact(sys, &view, nv, NULL, compact), 0);

    int32_t* one_thread = malloc(pixels * sizeof(*one_thread));
    int32_t* many_threads = malloc(pixels * sizeof(*many_threads));
    ASSERT_NOT_NULL(one_thread);
    ASSERT_NOT_NULL(many_threads);
    alea_slice_raster_t output;
    alea_slice_raster_init(&output);
    output.nu = nu;
    output.nv = nv;
    output.fields = ALEA_SLICE_RASTER_CELL_ID;
    output.cell_ids = one_thread;
#ifdef _OPENMP
    omp_set_num_threads(1);
#endif
    ASSERT_EQ(alea_rasterize_ray_slice_compact(&view, compact, &output), 0);
    output.cell_ids = many_threads;
#ifdef _OPENMP
    omp_set_num_threads(4);
#endif
    ASSERT_EQ(alea_rasterize_ray_slice_compact(&view, compact, &output), 0);
    ASSERT_EQ(memcmp(one_thread, many_threads, pixels * sizeof(*one_thread)), 0);

    free(many_threads);
    free(one_thread);
    alea_raycast_batch_result_destroy(compact);
    alea_destroy(sys);
}

TEST_MAIN()
