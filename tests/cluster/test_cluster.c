// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

#include "alea.h"
#include "alea_cluster.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int fail(int rank, const char* message) {
    fprintf(stderr, "rank %d: %s\n", rank, message);
    return 1;
}

static int stop_volume(size_t completed, size_t maximum, double error, void* data) {
    (void)completed; (void)maximum; (void)error;
    int* calls = data;
    ++*calls;
    return 1;
}

static int same_batch(const alea_raycast_batch_result_t* a,
                      const alea_raycast_batch_result_t* b) {
    size_t rays = alea_raycast_batch_ray_count(a);
    size_t segments = alea_raycast_batch_segment_count(a);
    size_t paths = alea_raycast_batch_path_entry_count(a);
    if (rays != alea_raycast_batch_ray_count(b) ||
        segments != alea_raycast_batch_segment_count(b) ||
        paths != alea_raycast_batch_path_entry_count(b) ||
        alea_raycast_batch_fields(a) != alea_raycast_batch_fields(b)) return 0;
#define SAME(field, count, type) do { \
    const type* left = alea_raycast_batch_##field(a); \
    const type* right = alea_raycast_batch_##field(b); \
    if ((count) && (!left || !right || \
        memcmp(left, right, (count) * sizeof(type)) != 0)) return 0; \
} while (0)
    SAME(ray_offsets, rays + 1, uint64_t);
    SAME(t_enter, segments, double);
    SAME(t_exit, segments, double);
    SAME(cell_ids, segments, int32_t);
    uint32_t fields = alea_raycast_batch_fields(a);
    if (fields & ALEA_RAY_BATCH_MATERIAL)
        SAME(material_ids, segments, int32_t);
    if (fields & ALEA_RAY_BATCH_DENSITY)
        SAME(densities, segments, double);
    if (fields & ALEA_RAY_BATCH_SURFACES) {
        SAME(enter_surface_ids, segments, int32_t);
        SAME(exit_surface_ids, segments, int32_t);
    }
    if (fields & ALEA_RAY_BATCH_RESOLUTION_FLAGS)
        SAME(resolution_flags, segments, uint8_t);
    if (fields & ALEA_RAY_BATCH_PROJECTED_OWNER) {
        SAME(projected_cell_ids, segments, int32_t);
        SAME(projected_material_ids, segments, int32_t);
        SAME(projected_universe_ids, segments, int32_t);
        SAME(projected_fill_universes, segments, int32_t);
        SAME(projected_depths, segments, int32_t);
        SAME(projected_is_lattice, segments, uint8_t);
        SAME(projected_occurrence_keys, segments, uint64_t);
    }
    if (fields & ALEA_RAY_BATCH_FULL_PATHS) {
        SAME(segment_path_offsets, segments + 1, uint64_t);
        SAME(path_cell_ids, paths, int32_t);
        SAME(path_material_ids, paths, int32_t);
        SAME(path_universe_ids, paths, int32_t);
        SAME(path_fill_universes, paths, int32_t);
        SAME(path_depths, paths, int32_t);
        SAME(path_is_lattice, paths, uint8_t);
        size_t path_xyz_count = paths * 3;
        SAME(path_lattice_origins_xyz, path_xyz_count, double);
        SAME(path_occurrence_keys, paths, uint64_t);
    }
#undef SAME
    return 1;
}

typedef struct {
    alea_system_t* sys;
    const double* origins;
    const double* directions;
    const alea_raycast_batch_options_t* options;
    size_t next_ray;
    size_t calls;
    size_t stop_after;
    double t_max;
    int mismatch;
} stream_check_t;

static int check_stream_batch(size_t first_ray,
                              const alea_raycast_batch_result_t* batch,
                              void* user_data) {
    stream_check_t* check = user_data;
    size_t count = alea_raycast_batch_ray_count(batch);
    alea_raycast_batch_result_t* serial = alea_raycast_batch_result_create();
    if (!serial || first_ray != check->next_ray || count == 0 ||
        alea_raycast_hier_batch(check->sys,
            check->origins + 3 * first_ray,
            check->directions + 3 * first_ray,
            count, check->t_max, check->options, serial) != 0 ||
        !same_batch(batch, serial)) check->mismatch = 1;
    alea_raycast_batch_result_destroy(serial);
    check->next_ray += count;
    ++check->calls;
    return check->mismatch ||
        (check->stop_after && check->calls == check->stop_after);
}

typedef struct {
    alea_system_t* sys;
    const alea_raycast_batch_options_t* options;
    size_t rays_seen;
    size_t calls;
    size_t stop_after;
    int mismatch;
} shard_check_t;

static int check_shard(size_t first_ray,
                       const alea_raycast_batch_result_t* shard,
                       void* user_data) {
    shard_check_t* check = user_data;
    size_t count = alea_raycast_batch_ray_count(shard);
    double* origins = calloc(count * 3, sizeof(double));
    double* directions = calloc(count * 3, sizeof(double));
    alea_raycast_batch_result_t* serial = alea_raycast_batch_result_create();
    if (!origins || !directions || !serial || !count)
        check->mismatch = 1;
    else {
        for (size_t i = 0; i < count; ++i) {
            origins[3 * i + 2] = -2.0;
            directions[3 * i + 2] = 1.0;
            if (first_ray + i == 1024) origins[3 * i] = 2.0;
        }
        if (alea_raycast_hier_batch(check->sys, origins, directions,
                count, 4.0, check->options, serial) != 0 ||
            !same_batch(shard, serial)) check->mismatch = 1;
    }
    free(origins); free(directions);
    alea_raycast_batch_result_destroy(serial);
    check->rays_seen += count;
    ++check->calls;
    return check->mismatch ||
        (check->stop_after && check->calls == check->stop_after);
}

static int same_coverage(const alea_ray_coverage_slice_result_t* a,
                         const alea_ray_coverage_slice_result_t* b) {
    size_t rows = alea_ray_coverage_slice_row_count(a);
    size_t intervals = alea_ray_coverage_slice_interval_count(a);
    size_t owners = alea_ray_coverage_slice_owner_count(a);
    if (rows != alea_ray_coverage_slice_row_count(b) ||
        intervals != alea_ray_coverage_slice_interval_count(b) ||
        owners != alea_ray_coverage_slice_owner_count(b)) return 0;
#define SAME(getter, count, type) do { \
    if ((count) && memcmp(getter(a), getter(b), \
                          (count) * sizeof(type))) return 0; \
} while (0)
    SAME(alea_ray_coverage_slice_row_offsets, rows + 1, size_t);
    SAME(alea_ray_coverage_slice_row_direction_tags, rows, uint8_t);
    SAME(alea_ray_coverage_slice_row_transverse_coordinates, rows, double);
    SAME(alea_ray_coverage_slice_t_enter, intervals, double);
    SAME(alea_ray_coverage_slice_t_exit, intervals, double);
    SAME(alea_ray_coverage_slice_kinds, intervals, uint8_t);
    SAME(alea_ray_coverage_slice_owner_offsets, intervals + 1, size_t);
    SAME(alea_ray_coverage_slice_owner_count_lower_bounds, intervals, size_t);
    SAME(alea_ray_coverage_slice_owner_cell_ids, owners, int);
    SAME(alea_ray_coverage_slice_owner_material_ids, owners, int);
    SAME(alea_ray_coverage_slice_owner_universe_ids, owners, int);
    SAME(alea_ray_coverage_slice_owner_fill_universes, owners, int);
    SAME(alea_ray_coverage_slice_owner_depths, owners, int);
    SAME(alea_ray_coverage_slice_owner_occurrence_keys, owners, uint64_t);
    SAME(alea_ray_coverage_slice_owner_parent_occurrence_keys, owners, uint64_t);
    SAME(alea_ray_coverage_slice_owner_resolution_flags, owners, uint8_t);
#undef SAME
    return 1;
}

typedef struct {
    alea_system_t* sys;
    alea_ray_coverage_slice_options_t options;
    size_t rays_seen;
    size_t calls;
    size_t stop_after;
    int mismatch;
} coverage_check_t;

static int check_coverage_shard(size_t first_row,
        const alea_ray_coverage_slice_result_t* shard, void* user_data) {
    coverage_check_t* check = user_data;
    size_t count = alea_ray_coverage_slice_row_count(shard);
    double* origins = calloc(count * 3, sizeof(double));
    double* directions = calloc(count * 3, sizeof(double));
    uint8_t* tags = calloc(count, sizeof(uint8_t));
    double* coordinates = calloc(count, sizeof(double));
    alea_ray_coverage_slice_result_t* serial =
        alea_ray_coverage_slice_result_create();
    if (!origins || !directions || !tags || !coordinates || !serial || !count)
        check->mismatch = 1;
    else {
        for (size_t i = 0; i < count; ++i) {
            origins[3 * i + 2] = -2.0;
            directions[3 * i + 2] = 1.0;
            tags[i] = (uint8_t)((first_row + i) & 1u);
            coordinates[i] = (double)(first_row + i);
            if (first_row + i == 256) origins[3 * i] = 2.0;
        }
        if (alea_ray_coverage_slice_query(check->sys, origins, directions,
                count, tags, coordinates, &check->options, serial) != 0 ||
            !same_coverage(shard, serial)) check->mismatch = 1;
    }
    alea_ray_coverage_slice_result_destroy(serial);
    free(origins); free(directions); free(tags); free(coordinates);
    check->rays_seen += count;
    ++check->calls;
    return check->mismatch ||
        (check->stop_after && check->calls == check->stop_after);
}

static int test_nested_batch(alea_cluster_t* cluster, int rank,
                             const alea_raycast_batch_options_t* options) {
    alea_system_t* nested = alea_create();
    if (!nested) return fail(rank, "nested system creation failed");
    int outer_surface = alea_sphere_surface(nested, 10, 0, 0, 0, 2);
    int inner_surface = alea_sphere_surface(nested, 20, 0, 0, 0, 1);
    int material = alea_add_material(nested, 1);
    if (outer_surface < 0 || inner_surface < 0 || material < 0)
        return fail(rank, "nested primitive creation failed");
    alea_node_id_t outer = alea_halfspace(nested, outer_surface, -1);
    alea_node_id_t inner = alea_halfspace(nested, inner_surface, -1);
    int container = alea_add_cell(nested, 10, outer,
                                  ALEA_MATERIAL_VOID, 0, 0);
    if (container < 0 || alea_set_fill(nested, container, 1, 0) != 0 ||
        alea_add_cell(nested, 20, inner, material, 1.0, 1) < 0 ||
        alea_build_universe_index(nested) != 0)
        return fail(rank, "nested fill setup failed");
    const double origins[9] = {0, 0, -3, 3, 0, -3, 0, 0, 0};
    const double directions[9] = {0, 0, 1, 0, 0, 1, 0, 0, 1};
    alea_raycast_batch_result_t* gathered = rank == 0
        ? alea_raycast_batch_result_create() : NULL;
    alea_raycast_batch_result_t* serial = rank == 0
        ? alea_raycast_batch_result_create() : NULL;
    alea_cluster_status_t status = alea_cluster_raycast_batch(
        cluster, nested, rank == 0 ? origins : NULL,
        rank == 0 ? directions : NULL, rank == 0 ? 3 : 0,
        rank == 0 ? 6.0 : 0.0, rank == 0 ? options : NULL, gathered);
    if (status != ALEA_CLUSTER_OK)
        return fail(rank, "nested packed batch failed");
    if (rank == 0) {
        if (alea_raycast_hier_batch(nested, origins, directions, 3, 6.0,
                                    options, serial) != 0 ||
            !same_batch(gathered, serial))
            return fail(rank, "nested packed batch differs from local");
        const uint64_t* offsets =
            alea_raycast_batch_segment_path_offsets(gathered);
        size_t segments = alea_raycast_batch_segment_count(gathered);
        int has_nested_path = 0;
        for (size_t i = 0; i < segments; ++i)
            if (offsets[i + 1] - offsets[i] >= 2) has_nested_path = 1;
        if (!has_nested_path)
            return fail(rank, "nested packed batch has no nested path");
    }
    alea_raycast_batch_result_destroy(gathered);
    alea_raycast_batch_result_destroy(serial);
    alea_destroy(nested);
    return 0;
}

static int test_render(alea_cluster_t* cluster, alea_system_t* sys,
                       int rank, render_mode_t mode, int aa, int edges,
                       int clipping, int width, int height, int auxiliary) {
    render_config_t cfg;
    render_config_init(&cfg);
    cfg.width = width;
    cfg.height = height;
    cfg.tile_size = 8;
    cfg.render_mode = mode;
    cfg.aa_samples = aa;
    cfg.edges = edges;
    cfg.shadows = mode == RENDER_MODE_SOLID;
    cfg.log_level = 0;
    cfg.threads = rank == 0 ? 1 : 2;
    if (clipping) {
        cfg.num_clips = 1;
        cfg.clips[0].normal[0] = 1.0;
        cfg.clips[0].d = 0.0;
    }
    render_camera_t cam;
    int local_error = render_camera_setup(&cam, &cfg, sys) != 0;
    render_framebuffer_t* gathered = rank == 0
        ? render_framebuffer_create(cfg.width, cfg.height, auxiliary) : NULL;
    render_framebuffer_t* serial = rank == 0
        ? render_framebuffer_create(cfg.width, cfg.height, auxiliary) : NULL;
    if (rank == 0 && (!gathered || !serial)) local_error = 1;
    alea_cluster_status_t status = alea_cluster_agree(cluster,
        local_error ? ALEA_CLUSTER_OUT_OF_MEMORY : ALEA_CLUSTER_OK);
    if (status != ALEA_CLUSTER_OK) return fail(rank, "render setup failed");
    status = alea_cluster_render_scene(cluster, sys, &cfg, &cam, gathered);
    if (status != ALEA_CLUSTER_OK)
        return fail(rank, "cluster render failed");
    if (rank == 0) {
        const size_t pixels = (size_t)cfg.width * cfg.height;
        if (render_scene(sys, &cfg, &cam, serial) != 0)
            local_error = 1;
        if (edges) render_edge_darken(serial);
        if (memcmp(gathered->color, serial->color,
                   pixels * 3 * sizeof(float)) ||
            memcmp(gathered->cell_id, serial->cell_id,
                   pixels * sizeof(int))) local_error = 1;
        if (auxiliary &&
            (memcmp(gathered->material_id, serial->material_id,
                    pixels * sizeof(int)) ||
             memcmp(gathered->depth, serial->depth,
                    pixels * sizeof(float)) ||
             memcmp(gathered->normal, serial->normal,
                    pixels * 3 * sizeof(float)))) local_error = 1;
    }
    status = alea_cluster_agree(cluster,
        local_error ? ALEA_CLUSTER_COMPUTE_ERROR : ALEA_CLUSTER_OK);
    render_framebuffer_free(gathered);
    render_framebuffer_free(serial);
    render_config_free(&cfg);
    return status == ALEA_CLUSTER_OK ? 0
        : fail(rank, "cluster render differs from serial");
}

static void free_slice_arrays(alea_slice_raster_t* raster) {
    free(raster->cell_ids); free(raster->material_ids);
    free(raster->universe_ids); free(raster->fill_universe_ids);
    free(raster->densities); free(raster->resolution_flags);
}

typedef struct {
    alea_system_t* sys;
    const alea_slice_view_t* views;
    const alea_slice_raster_options_t* options;
    alea_slice_raster_t* serial;
    size_t calls;
    size_t stop_after;
    int mismatch;
} slice_stack_check_t;

static int check_slice_plane(size_t index,
        const alea_slice_raster_t* raster, void* user_data) {
    slice_stack_check_t* check = user_data;
    const size_t pixels = raster->nu * raster->nv;
    if (index != check->calls ||
        alea_trace_ray_slice_raster(check->sys, &check->views[index],
            check->options, check->serial) != 0 ||
        memcmp(raster->cell_ids, check->serial->cell_ids,
               pixels * sizeof(int32_t)) ||
        memcmp(raster->densities, check->serial->densities,
               pixels * sizeof(double))) check->mismatch = 1;
    ++check->calls;
    return check->mismatch ||
        (check->stop_after && check->calls == check->stop_after);
}

static int test_slice_raster(alea_cluster_t* cluster, alea_system_t* sys,
                             int rank, size_t rows) {
    alea_slice_view_t view = {0};
    view.plane.normal[2] = 1.0;
    view.plane.u_axis[0] = 1.0;
    view.plane.v_axis[1] = 1.0;
    view.u_min = view.v_min = -1.5;
    view.u_max = view.v_max = 1.5;
    alea_slice_raster_t gathered, serial;
    alea_slice_raster_init(&gathered);
    alea_slice_raster_init(&serial);
    const size_t columns = 17, pixels = columns * rows;
    if (rank == 0) {
        gathered.nu = serial.nu = columns;
        gathered.nv = serial.nv = rows;
        gathered.fields = serial.fields =
            ALEA_SLICE_RASTER_CELL_ID | ALEA_SLICE_RASTER_MATERIAL_ID |
            ALEA_SLICE_RASTER_UNIVERSE_ID |
            ALEA_SLICE_RASTER_FILL_UNIVERSE |
            ALEA_SLICE_RASTER_DENSITY |
            ALEA_SLICE_RASTER_RESOLUTION_FLAGS;
#define ALLOC(member, type) do { \
    gathered.member = malloc(pixels * sizeof(type)); \
    serial.member = malloc(pixels * sizeof(type)); \
} while (0)
        ALLOC(cell_ids, int32_t); ALLOC(material_ids, int32_t);
        ALLOC(universe_ids, int32_t); ALLOC(fill_universe_ids, int32_t);
        ALLOC(densities, double); ALLOC(resolution_flags, uint8_t);
#undef ALLOC
    }
    int bad = rank == 0 &&
        (!gathered.cell_ids || !gathered.material_ids ||
         !gathered.universe_ids || !gathered.fill_universe_ids ||
         !gathered.densities || !gathered.resolution_flags ||
         !serial.cell_ids || !serial.material_ids ||
         !serial.universe_ids || !serial.fill_universe_ids ||
         !serial.densities || !serial.resolution_flags);
    alea_cluster_status_t status = alea_cluster_agree(cluster,
        bad ? ALEA_CLUSTER_OUT_OF_MEMORY : ALEA_CLUSTER_OK);
    if (status != ALEA_CLUSTER_OK) return fail(rank, "slice allocation failed");
    alea_slice_raster_options_t options;
    alea_slice_raster_options_init(&options);
    status = alea_cluster_slice_raster(cluster, sys, &view,
        rank == 0 ? &options : NULL, rank == 0 ? &gathered : NULL);
    if (status != ALEA_CLUSTER_OK)
        return fail(rank, "cluster slice raster failed");
    if (rank == 0) {
        if (alea_trace_ray_slice_raster(sys, &view, &options, &serial))
            bad = 1;
#define SAME(member, type) \
    if (memcmp(gathered.member, serial.member, \
               pixels * sizeof(type))) bad = 1
        SAME(cell_ids, int32_t); SAME(material_ids, int32_t);
        SAME(universe_ids, int32_t); SAME(fill_universe_ids, int32_t);
        SAME(densities, double); SAME(resolution_flags, uint8_t);
#undef SAME
    }
    status = alea_cluster_agree(cluster,
        bad ? ALEA_CLUSTER_COMPUTE_ERROR : ALEA_CLUSTER_OK);
    if (status != ALEA_CLUSTER_OK)
        return fail(rank, "cluster slice differs from serial");
    options.max_segments = 1;
    status = alea_cluster_slice_raster(cluster, sys, &view,
        rank == 0 ? &options : NULL, rank == 0 ? &gathered : NULL);
    if (status != ALEA_CLUSTER_OUTPUT_LIMIT)
        return fail(rank, "cluster slice global segment limit failed");
    options.max_segments = 0;
    alea_slice_view_t views[3] = {view, view, view};
    views[0].plane.origin[2] = -0.5;
    views[2].plane.origin[2] = 0.5;
    slice_stack_check_t stack = {sys, views, &options, &serial, 0, 0, 0};
    status = alea_cluster_slice_stack_stream(cluster, sys, views, 3,
        rank == 0 ? &options : NULL,
        rank == 0 ? &gathered : NULL,
        rank == 0 ? check_slice_plane : NULL,
        rank == 0 ? &stack : NULL);
    if (status != ALEA_CLUSTER_OK ||
        (rank == 0 && (stack.mismatch || stack.calls != 3)))
        return fail(rank, "cluster slice stack differs from serial");
    if (alea_cluster_size(cluster) > 1) {
        status = alea_cluster_slice_stack_stream(cluster, sys, views,
            rank == 0 ? 3 : 2, rank == 0 ? &options : NULL,
            rank == 0 ? &gathered : NULL,
            rank == 0 ? check_slice_plane : NULL,
            rank == 0 ? &stack : NULL);
        if (status != ALEA_CLUSTER_INVALID_ARGUMENT)
            return fail(rank, "slice stack count mismatch was not rejected");
    }
    stack.calls = 0;
    stack.stop_after = 2;
    status = alea_cluster_slice_stack_stream(cluster, sys, views, 3,
        rank == 0 ? &options : NULL,
        rank == 0 ? &gathered : NULL,
        rank == 0 ? check_slice_plane : NULL,
        rank == 0 ? &stack : NULL);
    if (status != ALEA_CLUSTER_INTERRUPTED ||
        (rank == 0 && (stack.mismatch || stack.calls != 2)))
        return fail(rank, "slice stack callback cancellation failed");
    free_slice_arrays(&gathered);
    free_slice_arrays(&serial);
    return 0;
}

static int same_mesh(const alea_mesh_result_t* a,
                     const alea_mesh_result_t* b) {
    if (!a || !b || a->nx != b->nx || a->ny != b->ny ||
        a->nz != b->nz || a->fields != b->fields ||
        a->bounds_source != b->bounds_source ||
        a->bounds_padding != b->bounds_padding ||
        a->num_materials != b->num_materials ||
        a->mixed_count != b->mixed_count ||
        a->fraction_count != b->fraction_count ||
        a->cell_fraction_count != b->cell_fraction_count)
        return 0;
    size_t cells = (size_t)a->nx * a->ny * a->nz;
#define SAME(member, count, type) do { \
    if ((count) && memcmp(a->member, b->member, \
                          (count) * sizeof(type))) return 0; \
} while (0)
    SAME(x_nodes, (size_t)a->nx + 1, double);
    SAME(y_nodes, (size_t)a->ny + 1, double);
    SAME(z_nodes, (size_t)a->nz + 1, double);
    SAME(unique_materials, (size_t)a->num_materials, int);
    if (a->fields & ALEA_MESH_FIELD_MATERIAL_ID)
        SAME(material_ids, cells, int);
    if (a->fields & ALEA_MESH_FIELD_CELL_ID)
        SAME(cell_ids, cells, int);
    if (a->fields & ALEA_MESH_FIELD_MIXED_FLAG)
        SAME(mixed_flags, cells, unsigned char);
    if (a->fields & ALEA_MESH_FIELD_DOMINANT_FRACTION)
        SAME(dominant_fractions, cells, double);
    if (a->fields & ALEA_MESH_FIELD_ESTIMATED_ERROR)
        SAME(estimated_errors, cells, double);
    if (a->fields & ALEA_MESH_FIELD_SAMPLE_COUNT)
        SAME(sample_counts, cells, uint32_t);
    if (a->fields & ALEA_MESH_FIELD_TIE_FLAG)
        SAME(tie_flags, cells, uint8_t);
    if (a->fields & ALEA_MESH_FIELD_REFINEMENT_FLAG)
        SAME(refinement_flags, cells, uint8_t);
    if (a->fields & ALEA_MESH_FIELD_SAMPLED_FRACTIONS) {
        SAME(fraction_spans, cells, alea_mesh_fraction_span_t);
        for (size_t i = 0; i < a->fraction_count; ++i)
            if (a->fractions[i].material_id != b->fractions[i].material_id ||
                a->fractions[i].fraction != b->fractions[i].fraction)
                return 0;
    }
    if (a->fields & ALEA_MESH_FIELD_CELL_FRACTIONS) {
        SAME(cell_fraction_spans, cells, alea_mesh_fraction_span_t);
        for (size_t i = 0; i < a->cell_fraction_count; ++i)
            if (a->cell_fractions[i].cell_id != b->cell_fractions[i].cell_id ||
                a->cell_fractions[i].material_id !=
                    b->cell_fractions[i].material_id ||
                a->cell_fractions[i].fraction !=
                    b->cell_fractions[i].fraction)
                return 0;
    }
#undef SAME
    return 1;
}

typedef struct {
    const alea_mesh_result_t* full;
    size_t expected_first;
    size_t expected_slices;
    size_t calls;
    int cancel;
    int mismatch;
} mesh_slab_check_t;

static int check_mesh_slab(size_t first,
        const alea_mesh_result_t* slab, void* user_data) {
    mesh_slab_check_t* check = user_data;
    const alea_mesh_result_t* full = check->full;
    const size_t nxy = (size_t)full->nx * full->ny;
    if (first != check->expected_first ||
        (size_t)slab->nz != check->expected_slices ||
        slab->bounds_source != full->bounds_source ||
        slab->bounds_padding != full->bounds_padding ||
        memcmp(slab->x_nodes, full->x_nodes,
               ((size_t)full->nx + 1) * sizeof(double)) ||
        memcmp(slab->y_nodes, full->y_nodes,
               ((size_t)full->ny + 1) * sizeof(double)) ||
        memcmp(slab->z_nodes, full->z_nodes + first,
               ((size_t)slab->nz + 1) * sizeof(double)))
        check->mismatch = 1;
    for (size_t i = 0; i < nxy * (size_t)slab->nz &&
                       !check->mismatch; ++i) {
        const size_t global = first * nxy + i;
        if (slab->material_ids[i] != full->material_ids[global] ||
            slab->cell_ids[i] != full->cell_ids[global] ||
            slab->sample_counts[i] != full->sample_counts[global] ||
            slab->dominant_fractions[i] != full->dominant_fractions[global] ||
            slab->estimated_errors[i] != full->estimated_errors[global] ||
            slab->refinement_flags[i] != full->refinement_flags[global])
            check->mismatch = 1;
        if (slab->fraction_spans && full->fraction_spans) {
            const alea_mesh_fraction_span_t x = slab->fraction_spans[i];
            const alea_mesh_fraction_span_t y = full->fraction_spans[global];
            if (x.count != y.count) check->mismatch = 1;
            for (uint32_t j = 0; j < x.count && !check->mismatch; ++j)
                if (slab->fractions[x.offset + j].material_id !=
                        full->fractions[y.offset + j].material_id ||
                    slab->fractions[x.offset + j].fraction !=
                        full->fractions[y.offset + j].fraction)
                    check->mismatch = 1;
        }
    }
    ++check->calls;
    return check->mismatch || check->cancel;
}

typedef struct {
    const alea_mesh_result_t* full;
    size_t next;
    size_t stop_after;
    int mismatch;
} mesh_visit_check_t;

static int check_mesh_visit(const alea_mesh_voxel_sample_t* sample,
                            void* user_data) {
    mesh_visit_check_t* check = user_data;
    const alea_mesh_result_t* full = check->full;
    const size_t nx = (size_t)full->nx, ny = (size_t)full->ny;
    const size_t index = check->next;
    const size_t i = index % nx;
    const size_t j = (index / nx) % ny;
    const size_t k = index / (nx * ny);
    if (sample->i != (int)i || sample->j != (int)j ||
        sample->k != (int)k ||
        sample->x_min != full->x_nodes[i] ||
        sample->x_max != full->x_nodes[i + 1] ||
        sample->y_min != full->y_nodes[j] ||
        sample->y_max != full->y_nodes[j + 1] ||
        sample->z_min != full->z_nodes[k] ||
        sample->z_max != full->z_nodes[k + 1] ||
        sample->material_id != full->material_ids[index] ||
        sample->cell_id != full->cell_ids[index] ||
        sample->mixed != full->mixed_flags[index] ||
        sample->tie_flags != full->tie_flags[index] ||
        sample->dominant_fraction != full->dominant_fractions[index] ||
        sample->sample_count != full->sample_counts[index] ||
        sample->estimated_error != full->estimated_errors[index] ||
        sample->refinement_flags != full->refinement_flags[index] ||
        sample->fraction_count != full->fraction_spans[index].count ||
        sample->cell_fraction_count != full->cell_fraction_spans[index].count)
        check->mismatch = 1;
    const alea_mesh_fraction_span_t materials = full->fraction_spans[index];
    for (uint32_t p = 0; p < sample->fraction_count && !check->mismatch; ++p)
        if (sample->fractions[p].material_id !=
                full->fractions[materials.offset + p].material_id ||
            sample->fractions[p].fraction !=
                full->fractions[materials.offset + p].fraction)
            check->mismatch = 1;
    const alea_mesh_fraction_span_t cells = full->cell_fraction_spans[index];
    for (uint32_t p = 0; p < sample->cell_fraction_count && !check->mismatch; ++p)
        if (sample->cell_fractions[p].cell_id !=
                full->cell_fractions[cells.offset + p].cell_id ||
            sample->cell_fractions[p].material_id !=
                full->cell_fractions[cells.offset + p].material_id ||
            sample->cell_fractions[p].fraction !=
                full->cell_fractions[cells.offset + p].fraction)
            check->mismatch = 1;
    ++check->next;
    return check->mismatch ||
        (check->stop_after && check->next == check->stop_after);
}

static int test_mesh(alea_cluster_t* cluster, alea_system_t* sys,
                     int rank, int custom_nodes,
                     alea_mesh_sampling_mode_t mode) {
    alea_mesh_config_t cfg;
    alea_mesh_config_init(&cfg);
    cfg.nx = 5; cfg.ny = 4; cfg.nz = 3;
    cfg.x_min = cfg.y_min = cfg.z_min = -1.5;
    cfg.x_max = cfg.y_max = cfg.z_max = 1.5;
    cfg.bounds_mode = ALEA_MESH_BOUNDS_EXPLICIT;
    cfg.sampling_mode = mode;
    cfg.subsamples_per_axis = 2;
    cfg.max_total_samples =
        mode == ALEA_MESH_SAMPLE_ADAPTIVE || mode == ALEA_MESH_SAMPLE_RAY ? 0 :
        (uint64_t)cfg.nx * cfg.ny * cfg.nz * 8;
    if (mode == ALEA_MESH_SAMPLE_ADAPTIVE) cfg.max_refine_depth = 2;
    if (mode == ALEA_MESH_SAMPLE_ADAPTIVE &&
        (custom_nodes == 9 || custom_nodes == 10))
        cfg.max_total_samples = (uint64_t)cfg.nx * cfg.ny * cfg.nz * 8 + 64;
    if (mode == ALEA_MESH_SAMPLE_RAY) {
        cfg.ray_directions = custom_nodes == 3 ? ALEA_MESH_RAY_X :
            custom_nodes == 4 ? ALEA_MESH_RAY_Y :
            (custom_nodes == 5 || custom_nodes == 7) ? ALEA_MESH_RAY_Z :
            (custom_nodes == 6 || custom_nodes == 8) ? ALEA_MESH_RAY_XYZ :
            ALEA_MESH_RAY_X | ALEA_MESH_RAY_Y;
        cfg.ray_origin_mode = ALEA_MESH_RAY_ORIGINS_SOBOL;
        cfg.ray_samples = 4;
        if (custom_nodes == 4) {
            cfg.ray_origin_mode = ALEA_MESH_RAY_ORIGINS_GRID;
            cfg.ray_grid_u = cfg.ray_grid_v = 2;
        }
    }
    cfg.workers = rank == 0 ? 1 : 2;
    const double xn[6] = {-1.5, -1.0, -0.2, 0.2, 1.0, 1.5};
    const double yn[5] = {-1.5, -0.5, 0.0, 0.7, 1.5};
    const double zn[4] = {-1.5, -0.3, 0.6, 1.5};
    if (custom_nodes == 1 || custom_nodes == 7) {
        cfg.x_nodes = xn; cfg.y_nodes = yn; cfg.z_nodes = zn;
    } else if (custom_nodes == 2 || custom_nodes == 8 ||
               custom_nodes == 10) {
        cfg.bounds_mode = ALEA_MESH_BOUNDS_AUTO;
        cfg.x_min = cfg.x_max = cfg.y_min = cfg.y_max = 0.0;
        cfg.z_min = cfg.z_max = 0.0;
    }
    alea_mesh_result_t* gathered = NULL;
    alea_cluster_status_t status = alea_cluster_mesh_sample(
        cluster, sys, &cfg, rank == 0 ? &gathered : NULL);
    if (status != ALEA_CLUSTER_OK)
        return fail(rank, "cluster mesh sampling failed");
    int bad = 0;
    if (rank == 0) {
        alea_mesh_result_t* serial = alea_mesh_sample(sys, &cfg);
        bad = !same_mesh(gathered, serial);
        alea_mesh_result_free(serial);
    }
    status = alea_cluster_agree(cluster,
        bad ? ALEA_CLUSTER_COMPUTE_ERROR : ALEA_CLUSTER_OK);
    alea_mesh_result_free(gathered);
    if (status == ALEA_CLUSTER_OK &&
        (custom_nodes == 2 || custom_nodes == 8 ||
         custom_nodes == 10) &&
        (mode == ALEA_MESH_SAMPLE_STRATIFIED ||
         mode == ALEA_MESH_SAMPLE_ADAPTIVE ||
         mode == ALEA_MESH_SAMPLE_RAY)) {
        alea_mesh_result_t* full = alea_mesh_sample(sys, &cfg);
        alea_cluster_status_t prepared = alea_cluster_agree(cluster,
            full ? ALEA_CLUSTER_OK : ALEA_CLUSTER_COMPUTE_ERROR);
        if (prepared != ALEA_CLUSTER_OK)
            return fail(rank, "mesh slab reference failed");
        const size_t ranks = (size_t)alea_cluster_size(cluster);
        const size_t mine = (size_t)cfg.nz / ranks +
            ((size_t)rank < (size_t)cfg.nz % ranks);
        const size_t first = (size_t)rank * ((size_t)cfg.nz / ranks) +
            ((size_t)rank < (size_t)cfg.nz % ranks
                ? (size_t)rank : (size_t)cfg.nz % ranks);
        mesh_slab_check_t check = {full, first, mine, 0, 0, 0};
        prepared = alea_cluster_mesh_sample_shards(cluster, sys, &cfg,
            check_mesh_slab, &check);
        if (prepared != ALEA_CLUSTER_OK || check.mismatch ||
            check.calls != (mine ? 1u : 0u))
            return fail(rank, "rank-owned mesh slab differs from serial");
        if (ranks > 1) {
            prepared = alea_cluster_mesh_sample_shards(cluster, sys, &cfg,
                rank == 1 ? NULL : check_mesh_slab, &check);
            if (prepared != ALEA_CLUSTER_INVALID_ARGUMENT)
                return fail(rank, "one-rank mesh callback failure was not agreed");
        }
        check.calls = 0;
        check.cancel = rank == 0;
        prepared = alea_cluster_mesh_sample_shards(cluster, sys, &cfg,
            check_mesh_slab, &check);
        if (prepared != ALEA_CLUSTER_INTERRUPTED)
            return fail(rank, "mesh slab callback cancellation failed");
        mesh_visit_check_t visit = {full, 0, 0, 0};
        prepared = alea_cluster_mesh_visit_root(cluster, sys, &cfg,
            rank == 0 ? check_mesh_visit : NULL,
            rank == 0 ? &visit : NULL);
        if (prepared != ALEA_CLUSTER_OK ||
            (rank == 0 && (visit.mismatch ||
                           visit.next != (size_t)cfg.nx * cfg.ny * cfg.nz)))
            return fail(rank, "ordered root mesh visitor differs from serial");
        visit.next = 0;
        prepared = alea_cluster_mesh_stream_root(cluster, sys, &cfg,
            rank == 0 ? check_mesh_visit : NULL,
            rank == 0 ? &visit : NULL);
        if (prepared != ALEA_CLUSTER_OK ||
            (rank == 0 && (visit.mismatch ||
                           visit.next != (size_t)cfg.nx * cfg.ny * cfg.nz)))
            return fail(rank, "streamed root mesh differs from serial");
        if (ranks > 1) {
            prepared = alea_cluster_mesh_stream_root(cluster, sys, &cfg,
                NULL, rank == 0 ? &visit : NULL);
            if (prepared != ALEA_CLUSTER_INVALID_ARGUMENT)
                return fail(rank, "root mesh stream callback failure was not agreed");
        }
        visit.next = 0;
        visit.stop_after = ranks > 1
            ? ((size_t)cfg.nz / ranks + (0u < (size_t)cfg.nz % ranks)) *
                (size_t)cfg.nx * cfg.ny + 1
            : 3;
        prepared = alea_cluster_mesh_stream_root(cluster, sys, &cfg,
            rank == 0 ? check_mesh_visit : NULL,
            rank == 0 ? &visit : NULL);
        if (prepared != ALEA_CLUSTER_INTERRUPTED ||
            (rank == 0 && visit.next != visit.stop_after))
            return fail(rank, "root mesh stream cancellation failed");
        visit.next = 0;
        visit.stop_after = 3;
        prepared = alea_cluster_mesh_visit_root(cluster, sys, &cfg,
            rank == 0 ? check_mesh_visit : NULL,
            rank == 0 ? &visit : NULL);
        alea_mesh_result_free(full);
        if (prepared != ALEA_CLUSTER_INTERRUPTED ||
            (rank == 0 && visit.next != 3))
            return fail(rank, "root mesh visitor cancellation failed");
    }
    return status == ALEA_CLUSTER_OK ? 0
        : fail(rank, "cluster mesh differs from serial");
}

static int same_validation(const alea_geom_validator_result_t* a,
                           const alea_geom_validator_result_t* b) {
    if (a->error_count != b->error_count ||
        a->crossings_checked != b->crossings_checked ||
        a->adjacency_hits != b->adjacency_hits ||
        a->exact_queries != b->exact_queries ||
        a->ambiguous_crossings != b->ambiguous_crossings ||
        a->suppressed_samples != b->suppressed_samples ||
        a->sample_limited_curves != b->sample_limited_curves ||
        a->incomplete_rays != b->incomplete_rays ||
        a->incomplete_slice_samples != b->incomplete_slice_samples ||
        a->truncated != b->truncated) return 0;
    for (size_t i = 0; i < a->error_count; ++i) {
        const alea_geom_error_t* x = &a->errors[i];
        const alea_geom_error_t* y = &b->errors[i];
        if (x->type != y->type || x->source != y->source ||
            x->previous_cell_id != y->previous_cell_id ||
            x->found_cell_id != y->found_cell_id ||
            x->surface_id != y->surface_id || x->flags != y->flags ||
            x->cause != y->cause ||
            x->t != y->t || x->curve_index != y->curve_index ||
            x->component_index != y->component_index ||
            memcmp(x->uv, y->uv, sizeof(x->uv)) != 0 ||
            memcmp(x->crossing_point, y->crossing_point,
                   sizeof(x->crossing_point)) != 0 ||
            memcmp(x->sample_point, y->sample_point,
                   sizeof(x->sample_point)) != 0) return 0;
    }
    return 1;
}

static int test_slice_validator(alea_cluster_t* cluster, alea_system_t* sys,
                                 int rank, size_t max_crossings) {
    alea_slice_view_t view = {0};
    view.plane.normal[2] = 1.0;
    view.plane.u_axis[0] = 1.0;
    view.plane.v_axis[1] = 1.0;
    view.u_min = view.v_min = -1.5;
    view.u_max = view.v_max = 1.5;
    alea_slice_curves_t* curves = alea_get_slice_curves(sys, &view);
    alea_cluster_status_t status = alea_cluster_agree(cluster,
        curves && alea_slice_curves_count(curves) >
            (max_crossings == 17 ? 4u : 0u)
            ? ALEA_CLUSTER_OK : ALEA_CLUSTER_COMPUTE_ERROR);
    if (status != ALEA_CLUSTER_OK)
        return fail(rank, "slice validation curve setup failed");
    alea_geom_validator_options_t options;
    alea_geom_validator_options_init(&options);
    options.max_samples_per_curve = 16;
    options.max_samples_per_signature = 2;
    options.max_crossings = max_crossings;
    options.max_errors = 5;
    alea_geom_validator_result_t gathered;
    alea_geom_validator_result_init(&gathered);
    status = alea_cluster_validate_slice_curves(cluster, sys, &view, curves,
        &options, rank == 0 ? &gathered : NULL);
    if (status != ALEA_CLUSTER_OK)
        return fail(rank, "cluster slice validation failed");
    int bad = 0;
    if (rank == 0) {
        alea_geom_validator_result_t serial;
        alea_geom_validator_result_init(&serial);
        bad = alea_validate_geometry_slice(sys, &view, curves,
            &options, &serial) != 0 || !same_validation(&gathered, &serial);
        alea_geom_validator_result_free(&serial);
    }
    status = alea_cluster_agree(cluster,
        bad ? ALEA_CLUSTER_COMPUTE_ERROR : ALEA_CLUSTER_OK);
    if (status == ALEA_CLUSTER_OK && alea_cluster_size(cluster) > 1 &&
        max_crossings == 0) {
        alea_geom_validator_options_t changed = options;
        if (rank == 1) changed.max_samples_per_curve++;
        alea_cluster_status_t mismatch = alea_cluster_validate_slice_curves(
            cluster, sys, &view, curves, &changed,
            rank == 0 ? &gathered : NULL);
        if (mismatch != ALEA_CLUSTER_INVALID_ARGUMENT)
            status = ALEA_CLUSTER_COMPUTE_ERROR;
        alea_slice_view_t changed_view = view;
        if (rank == 1) changed_view.u_max += 0.1;
        mismatch = alea_cluster_validate_slice_curves(cluster, sys,
            &changed_view, curves, &options,
            rank == 0 ? &gathered : NULL);
        if (mismatch != ALEA_CLUSTER_INVALID_ARGUMENT)
            status = ALEA_CLUSTER_COMPUTE_ERROR;
    }
    alea_geom_validator_result_free(&gathered);
    alea_slice_curves_free(curves);
    return status == ALEA_CLUSTER_OK ? 0
        : fail(rank, "cluster slice validation differs from serial");
}

static int test_validator(alea_cluster_t* cluster, alea_system_t* sys,
                          int rank, size_t max_crossings, int bounded_domain) {
    alea_geom_validator_options_t options;
    alea_geom_validator_options_init(&options);
    options.ray_count = 130;
    options.max_crossings = max_crossings;
    options.max_samples_per_signature = 2;
    if (bounded_domain) {
        options.flags |= ALEA_GEOM_VALIDATE_DOMAIN_BOUNDS;
        options.validation_bounds[0] = options.validation_bounds[2] =
            options.validation_bounds[4] = -2.0;
        options.validation_bounds[1] = options.validation_bounds[3] =
            options.validation_bounds[5] = 2.0;
        options.max_errors = 5;
    }
    alea_geom_validator_result_t gathered;
    alea_geom_validator_result_init(&gathered);
    alea_cluster_status_t status = alea_cluster_validate_geometry(
        cluster, sys, &options, rank == 0 ? &gathered : NULL);
    if (status != ALEA_CLUSTER_OK) {
        alea_geom_validator_result_free(&gathered);
        return fail(rank, "cluster validation failed");
    }
    int bad = 0;
    if (rank == 0) {
        alea_geom_validator_result_t serial;
        alea_geom_validator_result_init(&serial);
        bad = alea_validate_geometry(sys, &options, &serial) != 0 ||
              !same_validation(&gathered, &serial) ||
              (bounded_domain && serial.error_count == 0);
        alea_geom_validator_result_free(&serial);
    }
    status = alea_cluster_agree(cluster,
        bad ? ALEA_CLUSTER_COMPUTE_ERROR : ALEA_CLUSTER_OK);
    if (status == ALEA_CLUSTER_OK) {
        if (rank == 0) gathered.truncated = 1;
        alea_cluster_status_t skipped = alea_cluster_validate_geometry(
            cluster, sys, &options, rank == 0 ? &gathered : NULL);
        if (skipped != ALEA_CLUSTER_OK ||
            (rank == 0 && gathered.truncated != 1))
            status = ALEA_CLUSTER_COMPUTE_ERROR;
    }
    if (status == ALEA_CLUSTER_OK && alea_cluster_size(cluster) > 1 &&
        !max_crossings && !bounded_domain) {
        alea_geom_validator_options_t changed = options;
        if (rank == 1) changed.seed++;
        alea_cluster_status_t mismatch = alea_cluster_validate_geometry(
            cluster, sys, &changed, rank == 0 ? &gathered : NULL);
        if (mismatch != ALEA_CLUSTER_INVALID_ARGUMENT)
            status = ALEA_CLUSTER_COMPUTE_ERROR;
    }
    alea_geom_validator_result_free(&gathered);
    return status == ALEA_CLUSTER_OK ? 0
        : fail(rank, "cluster validation differs from serial");
}

static int test_validator_incomplete(alea_cluster_t* cluster,
                                     alea_system_t* sys, int rank) {
    alea_geom_validator_options_t options;
    alea_geom_validator_options_init(&options);
    options.ray_count = 32;
    options.max_errors = 64;
    options.max_breakpoints = 1;
    alea_geom_validator_result_t gathered;
    alea_geom_validator_result_init(&gathered);
    alea_cluster_status_t status = alea_cluster_validate_geometry(
        cluster, sys, &options, rank == 0 ? &gathered : NULL);
    int bad = status != ALEA_CLUSTER_OK;
    if (rank == 0 && !bad) {
        alea_geom_validator_result_t serial;
        alea_geom_validator_result_init(&serial);
        bad = alea_validate_geometry(sys, &options, &serial) != 0 ||
              serial.incomplete_rays == 0 ||
              !same_validation(&gathered, &serial);
        alea_geom_validator_result_free(&serial);
    }
    alea_geom_validator_result_free(&gathered);
    status = alea_cluster_agree(cluster,
        bad ? ALEA_CLUSTER_COMPUTE_ERROR : ALEA_CLUSTER_OK);
    return status == ALEA_CLUSTER_OK ? 0
        : fail(rank, "cluster incomplete-ray receipt differs from serial");
}

int main(int argc, char** argv) {
    alea_cluster_status_t status = alea_cluster_initialize(&argc, &argv);
    if (status != ALEA_CLUSTER_OK) return fail(-1, "cluster initialization failed");
    alea_cluster_t* cluster = alea_cluster_create();
    if (!cluster) return fail(-1, "cluster context creation failed");
    const int rank = alea_cluster_rank(cluster);

    char* file_data = NULL;
    size_t file_length = 0;
    status = alea_cluster_read_file(cluster,
        rank == 0 ? "README.md" : NULL, &file_data, &file_length);
    if (status != ALEA_CLUSTER_OK || !file_data || file_length < 10 ||
        file_data[file_length] != '\0' ||
        strstr(file_data, "# libalea.c") == NULL)
        return fail(rank, "collective file read failed");
    free(file_data);
    file_data = NULL;
    file_length = 42;
    status = alea_cluster_read_file(cluster,
        rank == 0 ? "/definitely/not/a/libalea/input/file" : NULL,
        &file_data, &file_length);
    if (status != ALEA_CLUSTER_IO_ERROR || file_data || file_length != 0)
        return fail(rank, "collective file error was not propagated");
    status = alea_cluster_read_file(cluster, rank == 0 ? "README.md" : NULL,
        rank == 1 ? NULL : &file_data, &file_length);
    if (alea_cluster_size(cluster) > 1 &&
        (status != ALEA_CLUSTER_INVALID_ARGUMENT || file_data))
        return fail(rank, "asymmetric argument error was not propagated");
    if (alea_cluster_size(cluster) == 1 && status != ALEA_CLUSTER_OK)
        return fail(rank, "single-rank file read failed");
    free(file_data);

    status = alea_cluster_read_mcnp_input(cluster,
        rank == 0 ? "tests/cluster/mcnp_deps_main.inp" : NULL,
        &file_data, &file_length);
    if (status != ALEA_CLUSTER_OK || !file_data ||
        file_length != strlen(file_data) ||
        !strstr(file_data, "1 so 1\n") ||
        !strstr(file_data, "m1 1001 1\n") ||
        strstr(file_data, "READ FILE=") ||
        strstr(file_data, "ReAd FiLe="))
        return fail(rank, "MCNP dependency expansion failed");
    free(file_data);
    file_data = NULL;
    status = alea_cluster_read_mcnp_input(cluster,
        rank == 0 ? "tests/cluster/mcnp_deps_missing.inp" : NULL,
        &file_data, &file_length);
    if (status != ALEA_CLUSTER_IO_ERROR || file_data || file_length != 0)
        return fail(rank, "missing MCNP dependency was not propagated");
    status = alea_cluster_read_mcnp_input(cluster,
        rank == 0 ? "tests/cluster/mcnp_deps_bad.inp" : NULL,
        &file_data, &file_length);
    if (status != ALEA_CLUSTER_COMPUTE_ERROR || file_data || file_length != 0)
        return fail(rank, "malformed MCNP dependency was not rejected");
    status = alea_cluster_read_mcnp_input(cluster,
        rank == 0 ? "tests/cluster/mcnp_deps_cycle.inp" : NULL,
        &file_data, &file_length);
    if (status != ALEA_CLUSTER_COMPUTE_ERROR || file_data || file_length != 0)
        return fail(rank, "cyclic MCNP dependency was not rejected");

    alea_system_t* sys = alea_create();
    if (!sys) return fail(rank, "system creation failed");
    int material = alea_add_material(sys, 1);
    int surface = alea_sphere_surface(sys, 1, 0.0, 0.0, 0.0, 1.0);
    alea_node_id_t interior = alea_halfspace(sys, surface, -1);
    if (material < 0 || surface < 0 ||
        alea_add_cell(sys, 1, interior, material, 1.0, 0) < 0)
        return fail(rank, "sphere model creation failed");

    if (test_render(cluster, sys, rank, RENDER_MODE_SOLID, 2, 1, 1,
                    35, 27, 1) ||
        test_render(cluster, sys, rank, RENDER_MODE_XRAY, 1, 0, 0,
                    35, 27, 1) ||
        test_render(cluster, sys, rank, RENDER_MODE_DEPTH, 1, 0, 0,
                    5, 5, 0))
        return 1;
    if (test_slice_raster(cluster, sys, rank, 13) ||
        test_slice_raster(cluster, sys, rank, 2))
        return 1;
    if (test_mesh(cluster, sys, rank, 0, ALEA_MESH_SAMPLE_CENTER) ||
        test_mesh(cluster, sys, rank, 1, ALEA_MESH_SAMPLE_CORNERS) ||
        test_mesh(cluster, sys, rank, 0, ALEA_MESH_SAMPLE_SUBCELL) ||
        test_mesh(cluster, sys, rank, 1, ALEA_MESH_SAMPLE_SUBCELL) ||
        test_mesh(cluster, sys, rank, 2, ALEA_MESH_SAMPLE_SUBCELL) ||
        test_mesh(cluster, sys, rank, 0, ALEA_MESH_SAMPLE_STRATIFIED) ||
        test_mesh(cluster, sys, rank, 1, ALEA_MESH_SAMPLE_STRATIFIED) ||
        test_mesh(cluster, sys, rank, 2, ALEA_MESH_SAMPLE_STRATIFIED) ||
        test_mesh(cluster, sys, rank, 0, ALEA_MESH_SAMPLE_ADAPTIVE) ||
        test_mesh(cluster, sys, rank, 2, ALEA_MESH_SAMPLE_ADAPTIVE) ||
        test_mesh(cluster, sys, rank, 9, ALEA_MESH_SAMPLE_ADAPTIVE) ||
        test_mesh(cluster, sys, rank, 10, ALEA_MESH_SAMPLE_ADAPTIVE) ||
        test_mesh(cluster, sys, rank, 0, ALEA_MESH_SAMPLE_RAY) ||
        test_mesh(cluster, sys, rank, 2, ALEA_MESH_SAMPLE_RAY) ||
        test_mesh(cluster, sys, rank, 3, ALEA_MESH_SAMPLE_RAY) ||
        test_mesh(cluster, sys, rank, 4, ALEA_MESH_SAMPLE_RAY) ||
        test_mesh(cluster, sys, rank, 5, ALEA_MESH_SAMPLE_RAY) ||
        test_mesh(cluster, sys, rank, 6, ALEA_MESH_SAMPLE_RAY) ||
        test_mesh(cluster, sys, rank, 7, ALEA_MESH_SAMPLE_RAY) ||
        test_mesh(cluster, sys, rank, 8, ALEA_MESH_SAMPLE_RAY))
        return 1;
    if (test_validator(cluster, sys, rank, 0, 0) ||
        test_validator_incomplete(cluster, sys, rank) ||
        test_validator(cluster, sys, rank, 8, 0) ||
        test_validator(cluster, sys, rank, 0, 1) ||
        test_slice_validator(cluster, sys, rank, 0) ||
        test_slice_validator(cluster, sys, rank, 8))
        return 1;
    alea_system_t* overlap = alea_create();
    if (!overlap) return fail(rank, "overlap model creation failed");
    int overlap_material = alea_add_material(overlap, 1);
    int overlap_a = alea_sphere_surface(overlap, 1, 0, 0, 0, 1);
    int overlap_b = alea_sphere_surface(overlap, 2, 0.4, 0, 0, 1);
    if (overlap_material < 0 || overlap_a < 0 || overlap_b < 0 ||
        alea_add_cell(overlap, 1,
            alea_halfspace(overlap, overlap_a, -1),
            overlap_material, 1.0, 0) < 0 ||
        alea_add_cell(overlap, 2,
            alea_halfspace(overlap, overlap_b, -1),
            overlap_material, 1.0, 0) < 0)
        return fail(rank, "overlap model setup failed");
    int overlap_bad = test_validator(cluster, overlap, rank, 0, 1) ||
        test_slice_validator(cluster, overlap, rank, 0) ||
        test_slice_validator(cluster, overlap, rank, 8);
    alea_destroy(overlap);
    if (overlap_bad) return 1;
    alea_system_t* many_curves = alea_create();
    if (!many_curves) return fail(rank, "multi-curve model creation failed");
    int curve_material = alea_add_material(many_curves, 1);
    for (int i = 0; i < 6; ++i) {
        int sid = alea_sphere_surface(many_curves, 20 + i,
            -1.0 + 0.4 * i, 0.0, 0.0, 0.35);
        if (curve_material < 0 || sid < 0 ||
            alea_add_cell(many_curves, 20 + i,
                alea_halfspace(many_curves, sid, -1),
                curve_material, 1.0, 0) < 0)
            return fail(rank, "multi-curve model setup failed");
    }
    int curves_bad = test_slice_validator(cluster, many_curves, rank, 0) ||
        test_slice_validator(cluster, many_curves, rank, 17);
    alea_destroy(many_curves);
    if (curves_bad) return 1;
    if (alea_cluster_size(cluster) > 1) {
        render_config_t mismatched_cfg;
        render_config_init(&mismatched_cfg);
        mismatched_cfg.width = rank == 0 ? 16 : 17;
        mismatched_cfg.height = 16;
        mismatched_cfg.log_level = 0;
        render_camera_t mismatched_cam;
        if (render_camera_setup(&mismatched_cam, &mismatched_cfg, sys))
            return fail(rank, "render mismatch setup failed");
        render_framebuffer_t* mismatch_frame = rank == 0
            ? render_framebuffer_create(16, 16, 0) : NULL;
        status = alea_cluster_agree(cluster,
            rank == 0 && !mismatch_frame
                ? ALEA_CLUSTER_OUT_OF_MEMORY : ALEA_CLUSTER_OK);
        if (status != ALEA_CLUSTER_OK)
            return fail(rank, "render mismatch allocation failed");
        status = alea_cluster_render_scene(cluster, sys, &mismatched_cfg,
                                           &mismatched_cam, mismatch_frame);
        render_framebuffer_free(mismatch_frame);
        render_config_free(&mismatched_cfg);
        if (status != ALEA_CLUSTER_INVALID_ARGUMENT)
            return fail(rank, "render config mismatch was not rejected");
    }

    const double ray_origins[9] = {0, 0, -2, 2, 0, -2, 0, 0, 0};
    const double ray_directions[9] = {0, 0, 1, 0, 0, 1, 0, 0, 1};
    unsigned char ray_hits[3] = {0};
    int32_t ray_cells[3] = {0};
    double ray_enters[3] = {0}, ray_exits[3] = {0};
    status = alea_cluster_raycast_first_segments(cluster, sys,
        rank == 0 ? ray_origins : NULL,
        rank == 0 ? ray_directions : NULL,
        rank == 0 ? 3 : 0, rank == 0 ? 4.0 : 0.0,
        rank == 0 ? ray_hits : NULL, rank == 0 ? ray_cells : NULL,
        rank == 0 ? ray_enters : NULL, rank == 0 ? ray_exits : NULL);
    if (status != ALEA_CLUSTER_OK)
        return fail(rank, "cluster raycast failed");
    if (rank == 0 && (!ray_hits[0] || ray_hits[1] || !ray_hits[2] ||
                      ray_cells[0] != 1 || ray_cells[2] != 1 ||
                      fabs(ray_enters[0] - 1.0) > 1e-8 ||
                      fabs(ray_exits[0] - 3.0) > 1e-8)) {
        fprintf(stderr, "raycast receipts: %u/%d %.12g %.12g; %u; %u/%d %.12g %.12g\n",
            ray_hits[0], ray_cells[0], ray_enters[0], ray_exits[0],
            ray_hits[1], ray_hits[2], ray_cells[2], ray_enters[2], ray_exits[2]);
        return fail(rank, "cluster raycast returned incorrect first segments");
    }

    alea_raycast_batch_options_t batch_options = {
        sizeof(batch_options),
        ALEA_RAY_BATCH_MATERIAL | ALEA_RAY_BATCH_DENSITY |
        ALEA_RAY_BATCH_SURFACES | ALEA_RAY_BATCH_RESOLUTION_FLAGS |
        ALEA_RAY_BATCH_PROJECTED_OWNER | ALEA_RAY_BATCH_FULL_PATHS,
        -1, 0, 0, 0
    };
    alea_raycast_batch_result_t* gathered = rank == 0
        ? alea_raycast_batch_result_create() : NULL;
    alea_raycast_batch_result_t* serial = rank == 0
        ? alea_raycast_batch_result_create() : NULL;
    status = alea_cluster_raycast_batch(cluster, sys,
        rank == 0 ? ray_origins : NULL,
        rank == 0 ? ray_directions : NULL,
        rank == 0 ? 3 : 0, rank == 0 ? 4.0 : 0.0,
        rank == 0 ? &batch_options : NULL, gathered);
    if (status != ALEA_CLUSTER_OK)
        return fail(rank, "cluster packed batch failed");
    if (rank == 0 && (alea_raycast_hier_batch(sys, ray_origins,
            ray_directions, 3, 4.0, &batch_options, serial) != 0 ||
            !same_batch(gathered, serial)))
        return fail(rank, "cluster packed batch differs from local batch");
    batch_options.max_output_bytes = 8;
    status = alea_cluster_raycast_batch(cluster, sys,
        rank == 0 ? ray_origins : NULL,
        rank == 0 ? ray_directions : NULL,
        rank == 0 ? 3 : 0, rank == 0 ? 4.0 : 0.0,
        rank == 0 ? &batch_options : NULL, gathered);
    if (status != ALEA_CLUSTER_OUTPUT_LIMIT ||
        (rank == 0 && !same_batch(gathered, serial)))
        return fail(rank, "cluster output limit changed prior result");
    batch_options.max_output_bytes = 0;
    batch_options.max_segments = 1;
    status = alea_cluster_raycast_batch(cluster, sys,
        rank == 0 ? ray_origins : NULL,
        rank == 0 ? ray_directions : NULL,
        rank == 0 ? 3 : 0, rank == 0 ? 4.0 : 0.0,
        rank == 0 ? &batch_options : NULL, gathered);
    if (status != ALEA_CLUSTER_OUTPUT_LIMIT ||
        (rank == 0 && !same_batch(gathered, serial)))
        return fail(rank, "cluster segment limit changed prior result");
    batch_options.max_segments = 0;
    alea_raycast_batch_result_destroy(gathered);
    alea_raycast_batch_result_destroy(serial);

    if (test_nested_batch(cluster, rank, &batch_options) != 0) return 1;

    gathered = rank == 0 ? alea_raycast_batch_result_create() : NULL;
    serial = rank == 0 ? alea_raycast_batch_result_create() : NULL;
    status = alea_cluster_raycast_batch(cluster, sys,
        rank == 0 ? ray_origins : NULL,
        rank == 0 ? ray_directions : NULL,
        rank == 0 ? 3 : 0, rank == 0 ? 4.0 : 0.0, NULL, gathered);
    if (status != ALEA_CLUSTER_OK)
        return fail(rank, "basic cluster packed batch failed");
    if (rank == 0 && (alea_raycast_hier_batch(sys, ray_origins,
            ray_directions, 3, 4.0, NULL, serial) != 0 ||
            !same_batch(gathered, serial) ||
            alea_raycast_batch_material_ids(gathered) != NULL ||
            alea_raycast_batch_segment_path_offsets(gathered) != NULL))
        return fail(rank, "basic cluster packed batch differs from local");
    alea_raycast_batch_result_destroy(gathered);
    alea_raycast_batch_result_destroy(serial);

    gathered = rank == 0 ? alea_raycast_batch_result_create() : NULL;
    serial = rank == 0 ? alea_raycast_batch_result_create() : NULL;
    status = alea_cluster_raycast_batch(cluster, sys, NULL, NULL, 0,
        rank == 0 ? 4.0 : 0.0,
        rank == 0 ? &batch_options : NULL, gathered);
    if (status != ALEA_CLUSTER_OK)
        return fail(rank, "empty cluster packed batch failed");
    if (rank == 0 && (alea_raycast_hier_batch(sys, NULL, NULL, 0,
            4.0, &batch_options, serial) != 0 ||
            !same_batch(gathered, serial)))
        return fail(rank, "empty cluster packed batch differs from local");
    alea_raycast_batch_result_destroy(gathered);
    alea_raycast_batch_result_destroy(serial);

    enum { MANY_RAYS = 1025 };
    double* many_origins = rank == 0
        ? calloc(MANY_RAYS * 3, sizeof(double)) : NULL;
    double* many_directions = rank == 0
        ? calloc(MANY_RAYS * 3, sizeof(double)) : NULL;
    unsigned char* many_hits = rank == 0
        ? calloc(MANY_RAYS, sizeof(unsigned char)) : NULL;
    int32_t* many_cells = rank == 0
        ? calloc(MANY_RAYS, sizeof(int32_t)) : NULL;
    double* many_enters = rank == 0
        ? calloc(MANY_RAYS, sizeof(double)) : NULL;
    double* many_exits = rank == 0
        ? calloc(MANY_RAYS, sizeof(double)) : NULL;
    if (rank == 0 && many_origins && many_directions) {
        for (size_t i = 0; i < MANY_RAYS; ++i) {
            many_origins[3 * i + 2] = -2.0;
            many_directions[3 * i + 2] = 1.0;
        }
        many_origins[3 * (MANY_RAYS - 1)] = 2.0;
    }
    status = alea_cluster_raycast_first_segments(cluster, sys,
        many_origins, many_directions, rank == 0 ? MANY_RAYS : 0,
        rank == 0 ? 4.0 : 0.0, many_hits, many_cells,
        many_enters, many_exits);
    if (status != ALEA_CLUSTER_OK)
        return fail(rank, "multi-batch cluster raycast failed");
    if (rank == 0 && (!many_hits[1023] || many_cells[1023] != 1 ||
                      many_hits[1024]))
        return fail(rank, "cluster raycast batch boundary was incorrect");
    gathered = rank == 0 ? alea_raycast_batch_result_create() : NULL;
    serial = rank == 0 ? alea_raycast_batch_result_create() : NULL;
    status = alea_cluster_raycast_batch(cluster, sys, many_origins,
        many_directions, rank == 0 ? MANY_RAYS : 0,
        rank == 0 ? 4.0 : 0.0, rank == 0 ? &batch_options : NULL, gathered);
    if (status != ALEA_CLUSTER_OK)
        return fail(rank, "multi-batch packed raycast failed");
    if (rank == 0 && (alea_raycast_hier_batch(sys, many_origins,
            many_directions, MANY_RAYS, 4.0, &batch_options, serial) != 0 ||
            !same_batch(gathered, serial)))
        return fail(rank, "multi-batch packed raycast differs from local");
    stream_check_t stream = {
        sys, many_origins, many_directions, &batch_options,
        0, 0, 0, 4.0, 0
    };
    status = alea_cluster_raycast_batch_stream(cluster, sys,
        many_origins, many_directions, rank == 0 ? MANY_RAYS : 0,
        rank == 0 ? 4.0 : 0.0,
        rank == 0 ? &batch_options : NULL,
        rank == 0 ? check_stream_batch : NULL,
        rank == 0 ? &stream : NULL);
    if (status != ALEA_CLUSTER_OK ||
        (rank == 0 && (stream.mismatch || stream.calls != 2 ||
                       stream.next_ray != MANY_RAYS)))
        return fail(rank, "streamed packed batches differ from local");
    stream.next_ray = stream.calls = 0;
    stream.stop_after = 1;
    status = alea_cluster_raycast_batch_stream(cluster, sys,
        many_origins, many_directions, rank == 0 ? MANY_RAYS : 0,
        rank == 0 ? 4.0 : 0.0,
        rank == 0 ? &batch_options : NULL,
        rank == 0 ? check_stream_batch : NULL,
        rank == 0 ? &stream : NULL);
    if (status != ALEA_CLUSTER_INTERRUPTED ||
        (rank == 0 && (stream.mismatch || stream.calls != 1 ||
                       stream.next_ray != 1024)))
        return fail(rank, "stream callback cancellation was not coordinated");
    stream.next_ray = stream.calls = stream.stop_after = 0;
    status = alea_cluster_raycast_batch_stream(cluster, sys,
        many_origins, many_directions, rank == 0 ? MANY_RAYS : 0,
        rank == 0 ? 4.0 : 0.0,
        rank == 0 ? &batch_options : NULL,
        rank == 0 ? check_stream_batch : NULL,
        rank == 0 ? &stream : NULL);
    if (status != ALEA_CLUSTER_OK ||
        (rank == 0 && (stream.mismatch || stream.calls != 2 ||
                       stream.next_ray != MANY_RAYS)))
        return fail(rank, "cluster did not recover after stream cancellation");
    stream.next_ray = stream.calls = 0;
    status = alea_cluster_raycast_batch_stream(cluster, sys,
        NULL, NULL, 0, rank == 0 ? 4.0 : 0.0,
        rank == 0 ? &batch_options : NULL,
        rank == 0 ? check_stream_batch : NULL,
        rank == 0 ? &stream : NULL);
    if (status != ALEA_CLUSTER_OK || (rank == 0 && stream.calls != 0))
        return fail(rank, "empty stream invoked callback");
    alea_raycast_batch_options_t bounded_options = {
        sizeof(bounded_options), 0, -1, 0, 0, 0
    };
    if (rank == 0) {
        size_t first_segments =
            (size_t)alea_raycast_batch_ray_offsets(gathered)[1024];
        bounded_options.max_output_bytes =
            (1024u + 1u) * sizeof(uint64_t) +
            first_segments * (2u * sizeof(double) + sizeof(int32_t));
    }
    status = alea_cluster_raycast_batch(cluster, sys,
        many_origins, many_directions, rank == 0 ? MANY_RAYS : 0,
        rank == 0 ? 4.0 : 0.0,
        rank == 0 ? &bounded_options : NULL, gathered);
    if (status != ALEA_CLUSTER_OUTPUT_LIMIT)
        return fail(rank, "full-result byte limit was not enforced");
    stream.options = &bounded_options;
    stream.next_ray = stream.calls = 0;
    status = alea_cluster_raycast_batch_stream(cluster, sys,
        many_origins, many_directions, rank == 0 ? MANY_RAYS : 0,
        rank == 0 ? 4.0 : 0.0,
        rank == 0 ? &bounded_options : NULL,
        rank == 0 ? check_stream_batch : NULL,
        rank == 0 ? &stream : NULL);
    if (status != ALEA_CLUSTER_OK ||
        (rank == 0 && (stream.mismatch || stream.calls != 2 ||
                       stream.next_ray != MANY_RAYS)))
        return fail(rank, "streaming did not honor per-batch byte limit");
    shard_check_t shard_check = {sys, &batch_options, 0, 0, 0, 0};
    status = alea_cluster_raycast_batch_shards(cluster, sys,
        many_origins, many_directions, rank == 0 ? MANY_RAYS : 0,
        rank == 0 ? 4.0 : 0.0,
        rank == 0 ? &batch_options : NULL,
        check_shard, &shard_check);
    if (status != ALEA_CLUSTER_OK || shard_check.mismatch ||
        shard_check.rays_seen == 0)
        return fail(rank, "rank-owned ray shard failed");
    if (alea_cluster_size(cluster) > 1) {
        shard_check.rays_seen = shard_check.calls = 0;
        shard_check.stop_after = rank == 1 ? 1 : 0;
        status = alea_cluster_raycast_batch_shards(cluster, sys,
            many_origins, many_directions, rank == 0 ? MANY_RAYS : 0,
            rank == 0 ? 4.0 : 0.0,
            rank == 0 ? &batch_options : NULL,
            check_shard, &shard_check);
        if (status != ALEA_CLUSTER_INTERRUPTED)
            return fail(rank, "non-root shard cancellation failed");
    }
    enum { COVERAGE_ROWS = 257 };
    double* coverage_origins = rank == 0
        ? calloc(COVERAGE_ROWS * 3, sizeof(double)) : NULL;
    double* coverage_directions = rank == 0
        ? calloc(COVERAGE_ROWS * 3, sizeof(double)) : NULL;
    uint8_t* coverage_tags = rank == 0
        ? calloc(COVERAGE_ROWS, sizeof(uint8_t)) : NULL;
    double* coverage_coordinates = rank == 0
        ? calloc(COVERAGE_ROWS, sizeof(double)) : NULL;
    if (rank == 0) {
        if (!coverage_origins || !coverage_directions ||
            !coverage_tags || !coverage_coordinates)
            return fail(rank, "coverage input allocation failed");
        for (size_t i = 0; i < COVERAGE_ROWS; ++i) {
            coverage_origins[3 * i + 2] = -2.0;
            coverage_directions[3 * i + 2] = 1.0;
            coverage_tags[i] = (uint8_t)(i & 1u);
            coverage_coordinates[i] = (double)i;
        }
        coverage_origins[3 * (COVERAGE_ROWS - 1)] = 2.0;
    }
    coverage_check_t coverage_check = {.sys = sys};
    alea_ray_coverage_slice_options_init(&coverage_check.options);
    coverage_check.options.t_max = 4.0;
    status = alea_cluster_coverage_shards(cluster, sys,
        coverage_origins, coverage_directions,
        rank == 0 ? COVERAGE_ROWS : 0,
        coverage_tags, coverage_coordinates,
        rank == 0 ? &coverage_check.options : NULL,
        check_coverage_shard, &coverage_check);
    if (status != ALEA_CLUSTER_OK || coverage_check.mismatch ||
        coverage_check.rays_seen == 0)
        return fail(rank, "rank-owned coverage shard failed");
    if (alea_cluster_size(cluster) > 1) {
        status = alea_cluster_coverage_shards(cluster, sys,
            coverage_origins, coverage_directions,
            rank == 0 ? COVERAGE_ROWS : 0,
            coverage_tags, coverage_coordinates,
            rank == 0 ? &coverage_check.options : NULL,
            rank == 1 ? NULL : check_coverage_shard, &coverage_check);
        if (status != ALEA_CLUSTER_INVALID_ARGUMENT)
            return fail(rank, "one-rank coverage callback failure was not agreed");
    }
    coverage_check.options.max_rows = 1;
    status = alea_cluster_coverage_shards(cluster, sys,
        coverage_origins, coverage_directions,
        rank == 0 ? COVERAGE_ROWS : 0,
        coverage_tags, coverage_coordinates,
        rank == 0 ? &coverage_check.options : NULL,
        check_coverage_shard, &coverage_check);
    if (status != ALEA_CLUSTER_OUTPUT_LIMIT)
        return fail(rank, "coverage shard row limit failed");
    coverage_check.options.max_rows = 0;
    if (alea_cluster_size(cluster) > 1) {
        coverage_check.calls = coverage_check.rays_seen = 0;
        coverage_check.stop_after = rank == 1 ? 1 : 0;
        status = alea_cluster_coverage_shards(cluster, sys,
            coverage_origins, coverage_directions,
            rank == 0 ? COVERAGE_ROWS : 0,
            coverage_tags, coverage_coordinates,
            rank == 0 ? &coverage_check.options : NULL,
            check_coverage_shard, &coverage_check);
        if (status != ALEA_CLUSTER_INTERRUPTED)
            return fail(rank, "non-root coverage cancellation failed");
    }
    coverage_check.calls = coverage_check.rays_seen = 0;
    coverage_check.stop_after = 0;
    status = alea_cluster_coverage_stream(cluster, sys,
        coverage_origins, coverage_directions,
        rank == 0 ? COVERAGE_ROWS : 0,
        coverage_tags, coverage_coordinates,
        rank == 0 ? &coverage_check.options : NULL,
        rank == 0 ? check_coverage_shard : NULL,
        rank == 0 ? &coverage_check : NULL);
    if (status != ALEA_CLUSTER_OK ||
        (rank == 0 && (coverage_check.mismatch ||
                       coverage_check.calls != 2 ||
                       coverage_check.rays_seen != COVERAGE_ROWS)))
        return fail(rank, "root coverage stream differs from serial");
    alea_ray_coverage_slice_result_t* whole = rank == 0
        ? alea_ray_coverage_slice_result_create() : NULL;
    alea_ray_coverage_slice_result_t* expected = rank == 0
        ? alea_ray_coverage_slice_result_create() : NULL;
    if (rank == 0 && (!whole || !expected))
        return fail(rank, "whole coverage result allocation failed");
    status = alea_cluster_coverage(cluster, sys,
        coverage_origins, coverage_directions,
        rank == 0 ? COVERAGE_ROWS : 0,
        coverage_tags, coverage_coordinates,
        rank == 0 ? &coverage_check.options : NULL, whole);
    if (status != ALEA_CLUSTER_OK ||
        (rank == 0 &&
         (alea_ray_coverage_slice_query(sys, coverage_origins,
             coverage_directions, COVERAGE_ROWS, coverage_tags,
             coverage_coordinates, &coverage_check.options, expected) != 0 ||
          !same_coverage(whole, expected))))
        return fail(rank, "whole coverage differs from serial");
    coverage_check.options.max_rows = COVERAGE_ROWS - 1;
    status = alea_cluster_coverage(cluster, sys,
        coverage_origins, coverage_directions,
        rank == 0 ? COVERAGE_ROWS : 0,
        coverage_tags, coverage_coordinates,
        rank == 0 ? &coverage_check.options : NULL, whole);
    if (status != ALEA_CLUSTER_OUTPUT_LIMIT ||
        (rank == 0 && !same_coverage(whole, expected)))
        return fail(rank, "whole coverage global limit failed");
    coverage_check.options.max_rows = 0;
    status = alea_cluster_coverage(cluster, sys, NULL, NULL, 0,
        NULL, NULL, rank == 0 ? &coverage_check.options : NULL, whole);
    if (status != ALEA_CLUSTER_OK ||
        (rank == 0 && (alea_ray_coverage_slice_row_count(whole) != 0 ||
                       !alea_ray_coverage_slice_row_offsets(whole) ||
                       !alea_ray_coverage_slice_owner_offsets(whole) ||
                       alea_ray_coverage_slice_row_offsets(whole)[0] != 0 ||
                       alea_ray_coverage_slice_owner_offsets(whole)[0] != 0)))
        return fail(rank, "empty whole coverage failed");
    coverage_check.options.max_output_bytes = 1;
    status = alea_cluster_coverage(cluster, sys, NULL, NULL, 0,
        NULL, NULL, rank == 0 ? &coverage_check.options : NULL, whole);
    if (status != ALEA_CLUSTER_OUTPUT_LIMIT ||
        (rank == 0 && !alea_ray_coverage_slice_row_offsets(whole)))
        return fail(rank, "empty whole coverage byte limit failed");
    coverage_check.options.max_output_bytes = 0;
    alea_ray_coverage_slice_result_destroy(whole);
    alea_ray_coverage_slice_result_destroy(expected);
    coverage_check.options.max_rows = 130;
    status = alea_cluster_coverage_stream(cluster, sys,
        coverage_origins, coverage_directions,
        rank == 0 ? COVERAGE_ROWS : 0,
        coverage_tags, coverage_coordinates,
        rank == 0 ? &coverage_check.options : NULL,
        rank == 0 ? check_coverage_shard : NULL,
        rank == 0 ? &coverage_check : NULL);
    if (status != ALEA_CLUSTER_OUTPUT_LIMIT)
        return fail(rank, "root coverage batch limit failed");
    coverage_check.options.max_rows = 0;
    coverage_check.calls = coverage_check.rays_seen = 0;
    coverage_check.stop_after = 1;
    status = alea_cluster_coverage_stream(cluster, sys,
        coverage_origins, coverage_directions,
        rank == 0 ? COVERAGE_ROWS : 0,
        coverage_tags, coverage_coordinates,
        rank == 0 ? &coverage_check.options : NULL,
        rank == 0 ? check_coverage_shard : NULL,
        rank == 0 ? &coverage_check : NULL);
    if (status != ALEA_CLUSTER_INTERRUPTED ||
        (rank == 0 && coverage_check.calls != 1))
        return fail(rank, "root coverage stream cancellation failed");
    free(coverage_origins); free(coverage_directions);
    free(coverage_tags); free(coverage_coordinates);
    alea_raycast_batch_result_destroy(gathered);
    alea_raycast_batch_result_destroy(serial);
    free(many_origins);
    free(many_directions);
    free(many_hits);
    free(many_cells);
    free(many_enters);
    free(many_exits);

    alea_volume_estimate_options_t options;
    alea_volume_estimate_options_init(&options);
    options.max_rays = 1003;
    options.batch_size = 211;
    options.seed = UINT64_C(987654321);
    options.requested_workers = rank == 0 ? 1 : 2;
    options.use_sampling_sphere = true;
    options.sampling_radius = 1.1;

    double volume = 0.0, error = 0.0;
    alea_cluster_volume_stats_t stats;
    status = alea_cluster_estimate_volumes(
        cluster, sys, &options, &volume, &error, &stats);
    if (status != ALEA_CLUSTER_OK)
        return fail(rank, alea_cluster_status_string(status));
    const double exact = 4.0 * 3.14159265358979323846 / 3.0;
    if (!(volume > exact * 0.7 && volume < exact * 1.3))
        return fail(rank, "volume estimate outside tolerance");
    if (!(error >= 0.0) || stats.volume.rays_completed != options.max_rays ||
        stats.rank_count != alea_cluster_size(cluster))
        return fail(rank, "invalid execution statistics");

    if (alea_cluster_is_root(cluster))
        printf("cluster=%s ranks=%d volume=%.8g rel_error=%.5g\n",
               alea_cluster_backend(cluster), alea_cluster_size(cluster),
               volume, error);

    int progress_calls = 0;
    options.progress = stop_volume;
    options.progress_user_data = &progress_calls;
    status = alea_cluster_estimate_volumes(
        cluster, sys, &options, &volume, &error, &stats);
    if (status != ALEA_CLUSTER_OK || !stats.volume.cancelled ||
        stats.volume.rays_completed != options.batch_size ||
        progress_calls != (rank == 0 ? 1 : 0))
        return fail(rank, "volume progress cancellation was not coordinated");
    options.progress = NULL;
    options.max_rays = options.batch_size = 1;
    options.target_rel_error = 1.0;
    options.sampling_radius = 1.0;
    status = alea_cluster_estimate_volumes(
        cluster, sys, &options, &volume, &error, &stats);
    if (status != ALEA_CLUSTER_OK || stats.volume.converged ||
        !(volume > 0.0) || !isinf(error))
        return fail(rank, "single ray falsely established volume uncertainty");
    if (rank == alea_cluster_size(cluster) - 1) alea_interrupt();
    status = alea_cluster_estimate_volumes(
        cluster, sys, &options, &volume, &error, &stats);
    alea_clear_interrupt();
    if (status != ALEA_CLUSTER_INTERRUPTED || !stats.volume.cancelled ||
        stats.volume.rays_completed != 0)
        return fail(rank, "volume interruption was not coordinated");

    if (alea_cluster_size(cluster) > 1) {
        alea_system_t* mismatch = alea_create();
        if (!mismatch) return fail(rank, "mismatch system creation failed");
        int mismatch_material = alea_add_material(mismatch, 1);
        int mismatch_surface = alea_sphere_surface(
            mismatch, 1, 0.0, 0.0, 0.0, 1.0);
        alea_node_id_t mismatch_interior =
            alea_halfspace(mismatch, mismatch_surface, -1);
        if (mismatch_material < 0 || mismatch_surface < 0 ||
            alea_add_cell(mismatch, rank == 0 ? 1 : 2, mismatch_interior,
                          mismatch_material,
                          1.0, 0) < 0)
            return fail(rank, "mismatch model creation failed");
        status = alea_cluster_estimate_volumes(
            cluster, mismatch, &options, &volume, &error, NULL);
        alea_destroy(mismatch);
        if (status != ALEA_CLUSTER_MODEL_MISMATCH)
            return fail(rank, "model mismatch was not detected");

        mismatch = alea_create();
        if (!mismatch) return fail(rank, "geometry mismatch setup failed");
        mismatch_material = alea_add_material(mismatch, 1);
        mismatch_surface = alea_sphere_surface(
            mismatch, 1, 0.0, 0.0, 0.0, rank == 0 ? 1.0 : 0.8);
        mismatch_interior = alea_halfspace(mismatch, mismatch_surface, -1);
        if (mismatch_material < 0 || mismatch_surface < 0 ||
            alea_add_cell(mismatch, 1, mismatch_interior,
                          mismatch_material, 1.0, 0) < 0)
            return fail(rank, "geometry mismatch model creation failed");
        status = alea_cluster_estimate_volumes(
            cluster, mismatch, &options, &volume, &error, NULL);
        if (status != ALEA_CLUSTER_MODEL_MISMATCH)
            return fail(rank, "primitive mismatch was not detected");
        status = alea_cluster_raycast_first_segments(cluster, mismatch,
            rank == 0 ? ray_origins : NULL,
            rank == 0 ? ray_directions : NULL,
            rank == 0 ? 3 : 0, rank == 0 ? 4.0 : 0.0,
            rank == 0 ? ray_hits : NULL, rank == 0 ? ray_cells : NULL,
            rank == 0 ? ray_enters : NULL, rank == 0 ? ray_exits : NULL);
        alea_destroy(mismatch);
        if (status != ALEA_CLUSTER_MODEL_MISMATCH)
            return fail(rank, "cluster raycast model mismatch was not detected");
    }

    alea_destroy(sys);
    alea_cluster_destroy(cluster);
    status = alea_cluster_finalize();
    return status == ALEA_CLUSTER_OK ? 0 : fail(rank, "cluster finalize failed");
}
