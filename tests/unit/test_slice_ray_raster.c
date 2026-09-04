// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

#include "alea_test.h"
#include "alea.h"
#include "alea_mcnp.h"
#include "alea_raycast.h"
#include "alea_slice.h"
#include "core/alea_system.h"

#include <stdlib.h>
#include <string.h>

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
    ASSERT_EQ(alea_rasterize_ray_slice_compact(&view, compact, &output), 0);
    output.cell_ids = many_threads;
    ASSERT_EQ(alea_rasterize_ray_slice_compact(&view, compact, &output), 0);
    ASSERT_EQ(memcmp(one_thread, many_threads, pixels * sizeof(*one_thread)), 0);

    free(many_threads);
    free(one_thread);
    alea_raycast_batch_result_destroy(compact);
    alea_destroy(sys);
}

TEST(ray_slice_raster_rejects_invalid_descriptors_and_trace_budget) {
    alea_system_t* sys = raster_test_system();
    ASSERT_NOT_NULL(sys);
    alea_slice_view_t view;
    raster_test_view(&view);
    int32_t cells[RASTER_PIXELS];
    for (size_t i = 0; i < RASTER_PIXELS; i++) cells[i] = 7654;
    alea_slice_raster_t output;
    alea_slice_raster_init(&output);
    output.nu = RASTER_NU;
    output.nv = RASTER_NV;
    output.fields = ALEA_SLICE_RASTER_CELL_ID;
    output.cell_ids = cells;

    output.fields |= 1u << 31;
    ASSERT_EQ(alea_trace_ray_slice_raster(sys, &view, NULL, &output), -1);
    output.fields = ALEA_SLICE_RASTER_CELL_ID | ALEA_SLICE_RASTER_MATERIAL_ID;
    output.material_ids = cells;
    ASSERT_EQ(alea_trace_ray_slice_raster(sys, &view, NULL, &output), -1);
    output.fields = ALEA_SLICE_RASTER_CELL_ID;
    output.material_ids = NULL;
    alea_slice_raster_options_t options;
    alea_slice_raster_options_init(&options);
    options.max_trace_output_bytes = 1;
    ASSERT_EQ(alea_trace_ray_slice_raster(sys, &view, &options, &output), -1);
    for (size_t i = 0; i < RASTER_PIXELS; i++) ASSERT_EQ(cells[i], 7654);
    alea_destroy(sys);
}

TEST(ray_slice_raster_lattice_projected_owner_matches_fused) {
    enum { nu = 24, nv = 24 };
    const size_t pixels = (size_t)nu * nv;
    mcnp_model_t* model = mcnp_load("tests/data/mcnp_lattice_eval.mcnp");
    if (!model) SKIP("lattice fixture not found");
    alea_system_t* sys = model->sys;
    alea_slice_view_t view;
    alea_slice_view_axis(&view, 2, 0.0, -1.0, 5.0, -1.0, 5.0);
    int32_t *compact_cells = calloc(pixels, sizeof(*compact_cells));
    int32_t *fused_cells = calloc(pixels, sizeof(*fused_cells));
    int32_t *compact_materials = calloc(pixels, sizeof(*compact_materials));
    int32_t *fused_materials = calloc(pixels, sizeof(*fused_materials));
    int32_t *compact_universes = calloc(pixels, sizeof(*compact_universes));
    int32_t *fused_universes = calloc(pixels, sizeof(*fused_universes));
    int32_t *compact_fills = calloc(pixels, sizeof(*compact_fills));
    int32_t *fused_fills = calloc(pixels, sizeof(*fused_fills));
    double *compact_densities = calloc(pixels, sizeof(*compact_densities));
    double *fused_densities = calloc(pixels, sizeof(*fused_densities));
    uint8_t *compact_flags = calloc(pixels, sizeof(*compact_flags));
    uint8_t *fused_flags = calloc(pixels, sizeof(*fused_flags));
    ASSERT_NOT_NULL(compact_cells); ASSERT_NOT_NULL(fused_cells);
    ASSERT_NOT_NULL(compact_materials); ASSERT_NOT_NULL(fused_materials);
    ASSERT_NOT_NULL(compact_universes); ASSERT_NOT_NULL(fused_universes);
    ASSERT_NOT_NULL(compact_fills); ASSERT_NOT_NULL(fused_fills);
    ASSERT_NOT_NULL(compact_densities); ASSERT_NOT_NULL(fused_densities);
    ASSERT_NOT_NULL(compact_flags); ASSERT_NOT_NULL(fused_flags);

    alea_slice_raster_t compact_output, fused_output;
    alea_slice_raster_init(&compact_output);
    alea_slice_raster_init(&fused_output);
    const alea_slice_raster_fields_t fields =
        ALEA_SLICE_RASTER_CELL_ID | ALEA_SLICE_RASTER_MATERIAL_ID |
        ALEA_SLICE_RASTER_UNIVERSE_ID | ALEA_SLICE_RASTER_FILL_UNIVERSE |
        ALEA_SLICE_RASTER_DENSITY | ALEA_SLICE_RASTER_RESOLUTION_FLAGS;
    compact_output = (alea_slice_raster_t){
        .struct_size = sizeof(compact_output), .nu = nu, .nv = nv, .fields = fields,
        .cell_ids = compact_cells, .material_ids = compact_materials,
        .universe_ids = compact_universes, .fill_universe_ids = compact_fills,
        .densities = compact_densities, .resolution_flags = compact_flags
    };
    fused_output = (alea_slice_raster_t){
        .struct_size = sizeof(fused_output), .nu = nu, .nv = nv, .fields = fields,
        .cell_ids = fused_cells, .material_ids = fused_materials,
        .universe_ids = fused_universes, .fill_universe_ids = fused_fills,
        .densities = fused_densities, .resolution_flags = fused_flags
    };
    alea_raycast_batch_options_t batch_options = {
        .struct_size = sizeof(batch_options),
        .fields = ALEA_RAY_BATCH_MATERIAL | ALEA_RAY_BATCH_DENSITY |
                  ALEA_RAY_BATCH_RESOLUTION_FLAGS | ALEA_RAY_BATCH_PROJECTED_OWNER,
        .projected_depth = 0
    };
    alea_raycast_batch_result_t* compact = alea_raycast_batch_result_create();
    ASSERT_NOT_NULL(compact);
    ASSERT_EQ(alea_trace_ray_slice_compact(sys, &view, nv, &batch_options, compact), 0);
    ASSERT_EQ(alea_rasterize_ray_slice_compact(&view, compact, &compact_output), 0);
    alea_slice_raster_options_t options;
    alea_slice_raster_options_init(&options);
    options.projected_depth = 0;
    ASSERT_EQ(alea_trace_ray_slice_raster(sys, &view, &options, &fused_output), 0);
    ASSERT_EQ(memcmp(compact_cells, fused_cells, pixels * sizeof(*compact_cells)), 0);
    ASSERT_EQ(memcmp(compact_materials, fused_materials, pixels * sizeof(*compact_materials)), 0);
    ASSERT_EQ(memcmp(compact_universes, fused_universes, pixels * sizeof(*compact_universes)), 0);
    ASSERT_EQ(memcmp(compact_fills, fused_fills, pixels * sizeof(*compact_fills)), 0);
    ASSERT_EQ(memcmp(compact_densities, fused_densities, pixels * sizeof(*compact_densities)), 0);
    ASSERT_EQ(memcmp(compact_flags, fused_flags, pixels * sizeof(*compact_flags)), 0);

    alea_raycast_batch_result_destroy(compact);
    free(fused_flags); free(compact_flags); free(fused_densities); free(compact_densities);
    free(fused_fills); free(compact_fills); free(fused_universes); free(compact_universes);
    free(fused_materials); free(compact_materials); free(fused_cells); free(compact_cells);
    mcnp_model_destroy(model);
}

TEST(ray_slice_raster_nested_lattice_leaf_matches_fused) {
    enum { nu = 32, nv = 16 };
    const size_t pixels = (size_t)nu * nv;
    mcnp_model_t* model = mcnp_load("tests/data/mcnp_nested_lattice.mcnp");
    if (!model) SKIP("nested lattice fixture not found");
    alea_slice_view_t view;
    alea_slice_view_axis(&view, 2, 0.0, -1.0, 3.0, -1.0, 1.0);
    int32_t* compact_cells = malloc(pixels * sizeof(*compact_cells));
    int32_t* fused_cells = malloc(pixels * sizeof(*fused_cells));
    int32_t* compact_materials = malloc(pixels * sizeof(*compact_materials));
    int32_t* fused_materials = malloc(pixels * sizeof(*fused_materials));
    ASSERT_NOT_NULL(compact_cells); ASSERT_NOT_NULL(fused_cells);
    ASSERT_NOT_NULL(compact_materials); ASSERT_NOT_NULL(fused_materials);
    alea_slice_raster_t compact_output, fused_output;
    alea_slice_raster_init(&compact_output);
    alea_slice_raster_init(&fused_output);
    compact_output.nu = fused_output.nu = nu;
    compact_output.nv = fused_output.nv = nv;
    compact_output.fields = fused_output.fields =
        ALEA_SLICE_RASTER_CELL_ID | ALEA_SLICE_RASTER_MATERIAL_ID;
    compact_output.cell_ids = compact_cells;
    compact_output.material_ids = compact_materials;
    fused_output.cell_ids = fused_cells;
    fused_output.material_ids = fused_materials;
    alea_raycast_batch_options_t batch_options = {
        .struct_size = sizeof(batch_options),
        .fields = ALEA_RAY_BATCH_MATERIAL,
        .projected_depth = -1
    };
    alea_raycast_batch_result_t* compact = alea_raycast_batch_result_create();
    ASSERT_NOT_NULL(compact);
    ASSERT_EQ(alea_trace_ray_slice_compact(model->sys, &view, nv,
                                            &batch_options, compact), 0);
    ASSERT_EQ(alea_rasterize_ray_slice_compact(&view, compact, &compact_output), 0);
    ASSERT_EQ(alea_trace_ray_slice_raster(model->sys, &view, NULL, &fused_output), 0);
    ASSERT_EQ(memcmp(compact_cells, fused_cells, pixels * sizeof(*compact_cells)), 0);
    ASSERT_EQ(memcmp(compact_materials, fused_materials,
                     pixels * sizeof(*compact_materials)), 0);
    int material_count = 0;
    for (size_t i = 0; i < pixels; i++) if (fused_materials[i] > 0) material_count++;
    ASSERT(material_count > 0);
    alea_raycast_batch_result_destroy(compact);
    free(fused_materials); free(compact_materials); free(fused_cells); free(compact_cells);
    mcnp_model_destroy(model);
}

TEST(ray_slice_raster_preserves_undefined_fill_flags) {
    const char* input =
        "Undefined fill raster test\n"
        "1 0 1 -2 3 -4 5 -6 FILL=10\n"
        "2 0 -1 : 2 : -3 : 4 : -5 : 6\n"
        "10 1 -1.0 -7 1 3 -4 5 -6 U=10\n"
        "\n"
        "1 PX -10\n2 PX 10\n3 PY -1\n4 PY 1\n5 PZ -1\n6 PZ 1\n7 PX 0\n"
        "\nM1 1001.80c 1.0\n";
    mcnp_model_t* model = mcnp_load_string(input, strlen(input));
    ASSERT_NOT_NULL(model);
    enum { nu = 32, nv = 4 };
    int32_t cells[nu * nv];
    uint8_t flags[nu * nv];
    alea_slice_view_t view;
    alea_slice_view_axis(&view, 2, 0.0, -10.0, 10.0, -1.0, 1.0);
    alea_slice_raster_t output;
    alea_slice_raster_init(&output);
    output.nu = nu;
    output.nv = nv;
    output.fields = ALEA_SLICE_RASTER_CELL_ID | ALEA_SLICE_RASTER_RESOLUTION_FLAGS;
    output.cell_ids = cells;
    output.resolution_flags = flags;
    ASSERT_EQ(alea_trace_ray_slice_raster(model->sys, &view, NULL, &output), 0);
    int saw_resolved = 0;
    int saw_undefined = 0;
    for (size_t i = 0; i < (size_t)nu * nv; i++) {
        if (cells[i] == 10 && flags[i] == 0) saw_resolved = 1;
        if (cells[i] == 1 && (flags[i] & ALEA_RESOLVE_UNDEFINED_FILL))
            saw_undefined = 1;
    }
    ASSERT(saw_resolved);
    ASSERT(saw_undefined);
    mcnp_model_destroy(model);
}

TEST_MAIN()
