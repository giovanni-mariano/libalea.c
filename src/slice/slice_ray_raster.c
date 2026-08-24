// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

#include "alea.h"
#include "alea_slice.h"
#include "raycast/raycast.h"

#include <math.h>
#include <stdint.h>
#include <stddef.h>
#include <string.h>

#define ALEA_SLICE_RASTER_ALL_FIELDS \
    (ALEA_SLICE_RASTER_CELL_ID | ALEA_SLICE_RASTER_MATERIAL_ID | \
     ALEA_SLICE_RASTER_UNIVERSE_ID | ALEA_SLICE_RASTER_FILL_UNIVERSE | \
     ALEA_SLICE_RASTER_DENSITY | ALEA_SLICE_RASTER_RESOLUTION_FLAGS)

/* A parallel region costs more than it saves on small viewports. */
#define ALEA_SLICE_RASTER_OMP_MIN_PIXELS ((size_t)65536)

typedef struct {
    uintptr_t begin;
    uintptr_t end;
} raster_range_t;

typedef struct {
    size_t pixels;
    size_t bytes;
} raster_layout_t;

static int raster_add_bytes(size_t* total, size_t count, size_t element_size) {
    if (count != 0 && element_size > SIZE_MAX / count) return -1;
    const size_t bytes = count * element_size;
    if (*total > SIZE_MAX - bytes) return -1;
    *total += bytes;
    return 0;
}

static int raster_add_range(raster_range_t* ranges, size_t* count,
                            const void* pointer, size_t bytes) {
    if (!pointer) return -1;
    const uintptr_t begin = (uintptr_t)pointer;
    if (bytes > UINTPTR_MAX - begin) return -1;
    ranges[*count].begin = begin;
    ranges[*count].end = begin + bytes;
    (*count)++;
    return 0;
}

static int raster_validate_descriptor(const alea_slice_raster_t* output,
                                      raster_layout_t* layout) {
    if (!output || !layout) {
        alea_set_error_detail(ALEA_ERR_NULL_ARG,
                              "ray slice raster requires an output descriptor");
        return -1;
    }
    if (output->struct_size < sizeof(*output)) {
        alea_set_error_detail(ALEA_ERR_INVALID_ARG,
                              "ray slice raster descriptor is too small");
        return -1;
    }
    if (output->nu == 0 || output->nv == 0 ||
        output->nu > SIZE_MAX / output->nv) {
        alea_set_error_detail(ALEA_ERR_OVERFLOW,
                              "ray slice raster dimensions are invalid or overflow");
        return -1;
    }
    if (output->fields == 0 || (output->fields & ~ALEA_SLICE_RASTER_ALL_FIELDS)) {
        alea_set_error_detail(ALEA_ERR_INVALID_ARG,
                              "ray slice raster field mask is empty or unsupported");
        return -1;
    }

    layout->pixels = output->nu * output->nv;
    layout->bytes = 0;
    raster_range_t ranges[6];
    size_t range_count = 0;
    if ((output->fields & ALEA_SLICE_RASTER_CELL_ID) &&
        (raster_add_bytes(&layout->bytes, layout->pixels, sizeof(*output->cell_ids)) ||
         raster_add_range(ranges, &range_count, output->cell_ids,
                          layout->pixels * sizeof(*output->cell_ids)))) goto invalid;
    if ((output->fields & ALEA_SLICE_RASTER_MATERIAL_ID) &&
        (raster_add_bytes(&layout->bytes, layout->pixels, sizeof(*output->material_ids)) ||
         raster_add_range(ranges, &range_count, output->material_ids,
                          layout->pixels * sizeof(*output->material_ids)))) goto invalid;
    if ((output->fields & ALEA_SLICE_RASTER_UNIVERSE_ID) &&
        (raster_add_bytes(&layout->bytes, layout->pixels, sizeof(*output->universe_ids)) ||
         raster_add_range(ranges, &range_count, output->universe_ids,
                          layout->pixels * sizeof(*output->universe_ids)))) goto invalid;
    if ((output->fields & ALEA_SLICE_RASTER_FILL_UNIVERSE) &&
        (raster_add_bytes(&layout->bytes, layout->pixels, sizeof(*output->fill_universe_ids)) ||
         raster_add_range(ranges, &range_count, output->fill_universe_ids,
                          layout->pixels * sizeof(*output->fill_universe_ids)))) goto invalid;
    if ((output->fields & ALEA_SLICE_RASTER_DENSITY) &&
        (raster_add_bytes(&layout->bytes, layout->pixels, sizeof(*output->densities)) ||
         raster_add_range(ranges, &range_count, output->densities,
                          layout->pixels * sizeof(*output->densities)))) goto invalid;
    if ((output->fields & ALEA_SLICE_RASTER_RESOLUTION_FLAGS) &&
        (raster_add_bytes(&layout->bytes, layout->pixels, sizeof(*output->resolution_flags)) ||
         raster_add_range(ranges, &range_count, output->resolution_flags,
                          layout->pixels * sizeof(*output->resolution_flags)))) goto invalid;

    for (size_t i = 0; i < range_count; i++)
        for (size_t j = i + 1; j < range_count; j++)
            if (ranges[i].begin < ranges[j].end && ranges[j].begin < ranges[i].end) {
                alea_set_error_detail(ALEA_ERR_INVALID_ARG,
                                      "requested ray slice raster buffers overlap");
                return -1;
            }
    return 0;

invalid:
    alea_set_error_detail(ALEA_ERR_INVALID_ARG,
                          "ray slice raster field buffer is missing or overflows");
    return -1;
}

static int raster_validate_view(const alea_slice_view_t* view) {
    if (!view) {
        alea_set_error_detail(ALEA_ERR_NULL_ARG, "ray slice raster requires a view");
        return -1;
    }
    if (!isfinite(view->u_min) || !isfinite(view->u_max) ||
        !isfinite(view->v_min) || !isfinite(view->v_max) ||
        !(view->u_max > view->u_min) || !(view->v_max > view->v_min)) {
        alea_set_error_detail(ALEA_ERR_INVALID_ARG, "ray slice raster view bounds are invalid");
        return -1;
    }
    return 0;
}

static size_t raster_pixel_index(double u, const alea_slice_view_t* view,
                                 size_t nu) {
    const double du = (view->u_max - view->u_min) / (double)nu;
    const double index = (u - view->u_min) / du - 0.5;
    if (index <= 0.0) return 0;
    if (index >= (double)nu) return nu;
    return (size_t)ceil(index);
}

static void raster_fill_i32(int32_t* out, size_t begin, size_t end, int32_t value) {
    for (size_t i = begin; i < end; i++) out[i] = value;
}

static void raster_fill_f64(double* out, size_t begin, size_t end, double value) {
    for (size_t i = begin; i < end; i++) out[i] = value;
}

static void raster_fill_u8(uint8_t* out, size_t begin, size_t end, uint8_t value) {
    for (size_t i = begin; i < end; i++) out[i] = value;
}

void alea_slice_raster_init(alea_slice_raster_t* raster) {
    if (!raster) return;
    memset(raster, 0, sizeof(*raster));
    raster->struct_size = sizeof(*raster);
}

void alea_slice_raster_options_init(alea_slice_raster_options_t* options) {
    if (!options) return;
    memset(options, 0, sizeof(*options));
    options->struct_size = sizeof(*options);
    options->projected_depth = -1;
}

int alea_rasterize_ray_slice_compact(
    const alea_slice_view_t* view,
    const alea_raycast_batch_result_t* segments,
    alea_slice_raster_t* output) {
    raster_layout_t layout;
    int projected_depth = -1;
    const uint64_t* offsets;
    const double* enters;
    const double* exits;
    const int32_t* cells;
    const int32_t* materials = NULL;
    const int32_t* universes = NULL;
    const int32_t* fills = NULL;
    const double* densities = NULL;
    const uint8_t* flags = NULL;
    const size_t segment_count = alea_raycast_batch_segment_count(segments);

    if (raster_validate_view(view) != 0 || raster_validate_descriptor(output, &layout) != 0)
        return -1;
    if (!segments || !alea_raycast_batch_result_get_compact_slice_provenance(
                         segments, view, output->nv, &projected_depth)) {
        alea_set_error_detail(ALEA_ERR_INVALID_ARG,
                              "ray slice raster requires a matching compact slice result");
        return -1;
    }

    offsets = alea_raycast_batch_ray_offsets(segments);
    enters = alea_raycast_batch_t_enter(segments);
    exits = alea_raycast_batch_t_exit(segments);
    cells = alea_raycast_batch_cell_ids(segments);
    const uint32_t fields = alea_raycast_batch_fields(segments);
    const int projected_owner_required =
        (output->fields & (ALEA_SLICE_RASTER_UNIVERSE_ID |
                           ALEA_SLICE_RASTER_FILL_UNIVERSE)) ||
        (projected_depth >= 0 &&
         (output->fields & (ALEA_SLICE_RASTER_CELL_ID |
                            ALEA_SLICE_RASTER_MATERIAL_ID)));
    if (!offsets || !enters || !exits || !cells ||
        (projected_owner_required && !(fields & ALEA_RAY_BATCH_PROJECTED_OWNER)) ||
        ((output->fields & ALEA_SLICE_RASTER_MATERIAL_ID) && projected_depth < 0 &&
         !(fields & ALEA_RAY_BATCH_MATERIAL)) ||
        ((output->fields & ALEA_SLICE_RASTER_DENSITY) &&
         !(fields & ALEA_RAY_BATCH_DENSITY)) ||
        ((output->fields & ALEA_SLICE_RASTER_RESOLUTION_FLAGS) &&
         !(fields & ALEA_RAY_BATCH_RESOLUTION_FLAGS))) {
        alea_set_error_detail(ALEA_ERR_INVALID_ARG,
                              "compact ray slice lacks requested raster source fields");
        return -1;
    }

    if (output->fields & ALEA_SLICE_RASTER_CELL_ID)
        cells = projected_depth >= 0 ? alea_raycast_batch_projected_cell_ids(segments) : cells;
    if (output->fields & ALEA_SLICE_RASTER_MATERIAL_ID)
        materials = projected_depth >= 0 ?
            alea_raycast_batch_projected_material_ids(segments) :
            alea_raycast_batch_material_ids(segments);
    if (output->fields & ALEA_SLICE_RASTER_UNIVERSE_ID)
        universes = alea_raycast_batch_projected_universe_ids(segments);
    if (output->fields & ALEA_SLICE_RASTER_FILL_UNIVERSE)
        fills = alea_raycast_batch_projected_fill_universes(segments);
    if (output->fields & ALEA_SLICE_RASTER_DENSITY)
        densities = alea_raycast_batch_densities(segments);
    if (output->fields & ALEA_SLICE_RASTER_RESOLUTION_FLAGS)
        flags = alea_raycast_batch_resolution_flags(segments);
    if (((output->fields & ALEA_SLICE_RASTER_CELL_ID) && !cells) ||
        ((output->fields & ALEA_SLICE_RASTER_MATERIAL_ID) && !materials) ||
        ((output->fields & ALEA_SLICE_RASTER_UNIVERSE_ID) && !universes) ||
        ((output->fields & ALEA_SLICE_RASTER_FILL_UNIVERSE) && !fills) ||
        ((output->fields & ALEA_SLICE_RASTER_DENSITY) && !densities) ||
        ((output->fields & ALEA_SLICE_RASTER_RESOLUTION_FLAGS) && !flags)) {
        alea_set_error_detail(ALEA_ERR_INVALID_ARG,
                              "compact ray slice source field pointer is missing");
        return -1;
    }

    if (offsets[0] != 0 || offsets[output->nv] != segment_count) {
        alea_set_error_detail(ALEA_ERR_INVALID_ARG, "compact ray slice CSR offsets are invalid");
        return -1;
    }
    for (size_t row = 0; row < output->nv; row++) {
        if (offsets[row] > offsets[row + 1] || offsets[row + 1] > segment_count) {
            alea_set_error_detail(ALEA_ERR_INVALID_ARG, "compact ray slice CSR offsets are invalid");
            return -1;
        }
        double previous_enter = -INFINITY;
        for (size_t i = (size_t)offsets[row]; i < (size_t)offsets[row + 1]; i++) {
            if (!isfinite(enters[i]) || !isfinite(exits[i]) ||
                enters[i] < view->u_min || exits[i] > view->u_max ||
                enters[i] > exits[i] || enters[i] < previous_enter) {
                alea_set_error_detail(ALEA_ERR_INVALID_ARG,
                                      "compact ray slice interval is invalid");
                return -1;
            }
            previous_enter = enters[i];
        }
    }

#ifdef _OPENMP
    #pragma omp parallel for if(layout.pixels >= ALEA_SLICE_RASTER_OMP_MIN_PIXELS) schedule(static)
#endif
    for (size_t row = 0; row < output->nv; row++) {
        const size_t row_offset = row * output->nu;
        const size_t row_end = row_offset + output->nu;
        if (output->fields & ALEA_SLICE_RASTER_CELL_ID)
            raster_fill_i32(output->cell_ids, row_offset, row_end, -1);
        if (output->fields & ALEA_SLICE_RASTER_MATERIAL_ID)
            raster_fill_i32(output->material_ids, row_offset, row_end, 0);
        if (output->fields & ALEA_SLICE_RASTER_UNIVERSE_ID)
            raster_fill_i32(output->universe_ids, row_offset, row_end, -1);
        if (output->fields & ALEA_SLICE_RASTER_FILL_UNIVERSE)
            raster_fill_i32(output->fill_universe_ids, row_offset, row_end, -1);
        if (output->fields & ALEA_SLICE_RASTER_DENSITY)
            raster_fill_f64(output->densities, row_offset, row_end, 0.0);
        if (output->fields & ALEA_SLICE_RASTER_RESOLUTION_FLAGS)
            raster_fill_u8(output->resolution_flags, row_offset, row_end, 0);
        for (size_t i = (size_t)offsets[row]; i < (size_t)offsets[row + 1]; i++) {
            const size_t lo = raster_pixel_index(enters[i], view, output->nu);
            const size_t hi = raster_pixel_index(exits[i], view, output->nu);
            if (lo >= hi) continue;
            const size_t begin = row_offset + lo;
            const size_t end = row_offset + hi;
            if (output->fields & ALEA_SLICE_RASTER_CELL_ID)
                raster_fill_i32(output->cell_ids, begin, end, cells[i]);
            if (output->fields & ALEA_SLICE_RASTER_MATERIAL_ID)
                raster_fill_i32(output->material_ids, begin, end, materials[i]);
            if (output->fields & ALEA_SLICE_RASTER_UNIVERSE_ID)
                raster_fill_i32(output->universe_ids, begin, end, universes[i]);
            if (output->fields & ALEA_SLICE_RASTER_FILL_UNIVERSE)
                raster_fill_i32(output->fill_universe_ids, begin, end, fills[i]);
            if (output->fields & ALEA_SLICE_RASTER_DENSITY)
                raster_fill_f64(output->densities, begin, end, densities[i]);
            if (output->fields & ALEA_SLICE_RASTER_RESOLUTION_FLAGS)
                raster_fill_u8(output->resolution_flags, begin, end, flags[i]);
        }
    }
    return 0;
}

int alea_trace_ray_slice_raster(
    alea_system_t* sys,
    const alea_slice_view_t* view,
    const alea_slice_raster_options_t* input_options,
    alea_slice_raster_t* output) {
    raster_layout_t layout;
    alea_slice_raster_options_t defaults;
    alea_slice_raster_options_t options;
    alea_raycast_batch_options_t trace_options;
    alea_raycast_batch_result_t* compact;
    uint32_t fields = 0;

    if (!sys || !view) {
        alea_set_error_detail(ALEA_ERR_NULL_ARG,
                              "ray slice raster trace requires a system and view");
        return -1;
    }
    if (raster_validate_view(view) != 0 || raster_validate_descriptor(output, &layout) != 0)
        return -1;
    alea_slice_raster_options_init(&defaults);
    options = defaults;
    if (input_options) {
        if (input_options->struct_size < sizeof(*input_options)) {
            alea_set_error_detail(ALEA_ERR_INVALID_ARG,
                                  "ray slice raster options are too small");
            return -1;
        }
        options = *input_options;
    }
    if (options.projected_depth < -1) {
        alea_set_error_detail(ALEA_ERR_INVALID_ARG,
                              "ray slice raster projected depth is invalid");
        return -1;
    }

    if ((output->fields & (ALEA_SLICE_RASTER_UNIVERSE_ID |
                           ALEA_SLICE_RASTER_FILL_UNIVERSE)) ||
        (options.projected_depth >= 0 &&
         (output->fields & (ALEA_SLICE_RASTER_CELL_ID |
                            ALEA_SLICE_RASTER_MATERIAL_ID))))
        fields |= ALEA_RAY_BATCH_PROJECTED_OWNER;
    if ((output->fields & ALEA_SLICE_RASTER_MATERIAL_ID) && options.projected_depth < 0)
        fields |= ALEA_RAY_BATCH_MATERIAL;
    if (output->fields & ALEA_SLICE_RASTER_DENSITY)
        fields |= ALEA_RAY_BATCH_DENSITY;
    if (output->fields & ALEA_SLICE_RASTER_RESOLUTION_FLAGS)
        fields |= ALEA_RAY_BATCH_RESOLUTION_FLAGS;

    trace_options = (alea_raycast_batch_options_t){
        .struct_size = sizeof(trace_options),
        .fields = fields,
        .projected_depth = options.projected_depth,
        .max_segments = options.max_segments,
        .max_path_entries = 0,
        .max_output_bytes = options.max_trace_output_bytes
    };
    compact = alea_raycast_batch_result_create();
    if (!compact) {
        alea_set_error_detail(ALEA_ERR_OUT_OF_MEMORY,
                              "failed to allocate compact ray slice result");
        return -1;
    }
    const int trace_rc = alea_trace_ray_slice_compact(sys, view, output->nv,
                                                       &trace_options, compact);
    const int raster_rc = trace_rc == 0 ?
        alea_rasterize_ray_slice_compact(view, compact, output) : -1;
    alea_raycast_batch_result_destroy(compact);
    return trace_rc != 0 ? trace_rc : raster_rc;
}
