// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

#include "cluster_internal.h"
#include "mesh/mesh_cluster_internal.h"

#include <limits.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#define MESH_FIELDS (ALEA_MESH_FIELD_MATERIAL_ID | ALEA_MESH_FIELD_CELL_ID | \
    ALEA_MESH_FIELD_MIXED_FLAG | ALEA_MESH_FIELD_DOMINANT_FRACTION | \
    ALEA_MESH_FIELD_SAMPLED_FRACTIONS | ALEA_MESH_FIELD_SAMPLE_COUNT | \
    ALEA_MESH_FIELD_TIE_FLAG | ALEA_MESH_FIELD_ESTIMATED_ERROR | \
    ALEA_MESH_FIELD_REFINEMENT_FLAG | ALEA_MESH_FIELD_CELL_FRACTIONS)

static size_t slab_count(size_t total, size_t rank, size_t ranks) {
    return total / ranks + (rank < total % ranks);
}

static size_t slab_first(size_t total, size_t rank, size_t ranks) {
    return rank * (total / ranks) + (rank < total % ranks ? rank : total % ranks);
}

static uint64_t mesh_hash(uint64_t hash, const void* bytes, size_t count) {
    const unsigned char* data = bytes;
    for (size_t i = 0; i < count; ++i) {
        hash ^= data[i];
        hash *= UINT64_C(1099511628211);
    }
    return hash;
}

static uint64_t config_fingerprint(const alea_mesh_config_t* cfg) {
    uint64_t hash = UINT64_C(1469598103934665603);
#define HASH(value) (hash = mesh_hash(hash, &(value), sizeof(value)))
    HASH(cfg->x_min); HASH(cfg->x_max); HASH(cfg->y_min); HASH(cfg->y_max);
    HASH(cfg->z_min); HASH(cfg->z_max);
    HASH(cfg->nx); HASH(cfg->ny); HASH(cfg->nz);
    HASH(cfg->format); HASH(cfg->void_material_id); HASH(cfg->auto_pad);
    HASH(cfg->sampling_mode); HASH(cfg->subsamples_per_axis);
    HASH(cfg->mixed_threshold); HASH(cfg->target_error);
    HASH(cfg->max_refine_depth); HASH(cfg->max_samples_per_voxel);
    HASH(cfg->max_total_samples); HASH(cfg->sampling_seed);
    HASH(cfg->bounds_mode); HASH(cfg->fields);
    HASH(cfg->ray_grid_u); HASH(cfg->ray_grid_v);
    HASH(cfg->ray_origin_mode); HASH(cfg->ray_samples);
    HASH(cfg->ray_point_count); HASH(cfg->ray_directions);
    if (cfg->sampling_mode == ALEA_MESH_SAMPLE_RAY &&
        cfg->ray_origin_mode == ALEA_MESH_RAY_ORIGINS_CUSTOM &&
        cfg->ray_points)
        hash = mesh_hash(hash, cfg->ray_points,
            (size_t)cfg->ray_point_count * 2 * sizeof(double));
    int has_x = cfg->x_nodes != NULL;
    int has_y = cfg->y_nodes != NULL;
    int has_z = cfg->z_nodes != NULL;
    HASH(has_x); HASH(has_y); HASH(has_z);
    if (has_x) hash = mesh_hash(hash, cfg->x_nodes,
                                ((size_t)cfg->nx + 1) * sizeof(double));
    if (has_y) hash = mesh_hash(hash, cfg->y_nodes,
                                ((size_t)cfg->ny + 1) * sizeof(double));
    if (has_z) hash = mesh_hash(hash, cfg->z_nodes,
                                ((size_t)cfg->nz + 1) * sizeof(double));
#undef HASH
    return hash ? hash : 1;
}

static int compare_ints(const void* left, const void* right) {
    const int a = *(const int*)left, b = *(const int*)right;
    return (a > b) - (a < b);
}

static alea_cluster_status_t sample_local_slab(
        alea_cluster_t* cluster, alea_system_t* sys,
        const alea_mesh_config_t* part_config, size_t slices, size_t first,
        double global_z_min, double global_z_max, int ordered_budget,
        alea_mesh_result_t** output) {
    *output = NULL;
    const size_t ranks = (size_t)cluster->size;
    const size_t rank = (size_t)cluster->rank;
    uint64_t cumulative = 0;
    const size_t passes = ordered_budget ? ranks : 1;
    for (size_t peer = 0; peer < passes; ++peer) {
        const int active = slices && (!ordered_budget || rank == peer);
        uint64_t next = cumulative;
        if (active)
            *output = alea_mesh_sample_with_z_offset(sys, part_config,
                (int)first, global_z_min, global_z_max,
                cumulative, &next);
        alea_cluster_status_t local = !active || *output
            ? ALEA_CLUSTER_OK
            : alea_interrupted() ? ALEA_CLUSTER_INTERRUPTED
                                 : ALEA_CLUSTER_COMPUTE_ERROR;
        alea_cluster_status_t status = alea_cluster_agree(cluster, local);
        if (status != ALEA_CLUSTER_OK) return status;
        if (ordered_budget) {
            if (alea_cluster_backend_broadcast_u64(cluster->backend,
                    &next, (int)peer)) {
                cluster->usable = 0;
                return ALEA_CLUSTER_BACKEND_ERROR;
            }
            cumulative = next;
        }
    }
    return ALEA_CLUSTER_OK;
}

alea_cluster_status_t alea_cluster_mesh_sample(
        alea_cluster_t* cluster, alea_system_t* sys,
        const alea_mesh_config_t* cfg,
        alea_mesh_result_t** root_result) {
    if (!cluster) return ALEA_CLUSTER_INVALID_ARGUMENT;
    if (root_result) *root_result = NULL;
    alea_cluster_status_t local = ALEA_CLUSTER_OK;
    if (!sys || !cfg || cfg->nx <= 0 || cfg->ny <= 0 || cfg->nz <= 0 ||
        (cfg->sampling_mode != ALEA_MESH_SAMPLE_CENTER &&
         cfg->sampling_mode != ALEA_MESH_SAMPLE_CORNERS &&
         cfg->sampling_mode != ALEA_MESH_SAMPLE_SUBCELL &&
         cfg->sampling_mode != ALEA_MESH_SAMPLE_STRATIFIED &&
         cfg->sampling_mode != ALEA_MESH_SAMPLE_ADAPTIVE &&
         cfg->sampling_mode != ALEA_MESH_SAMPLE_RAY) ||
        (cfg->sampling_mode == ALEA_MESH_SAMPLE_RAY &&
         (!(cfg->ray_directions & ALEA_MESH_RAY_XYZ) ||
          cfg->max_total_samples != 0 ||
          (cfg->ray_origin_mode == ALEA_MESH_RAY_ORIGINS_CUSTOM &&
           (uint64_t)cfg->ray_point_count * 2u * sizeof(double) >
               SIZE_MAX))) ||
        cfg->visit || cfg->progress || (cfg->fields & ~MESH_FIELDS))
        local = ALEA_CLUSTER_INVALID_ARGUMENT;
    if (cluster->rank == 0 && !root_result)
        local = ALEA_CLUSTER_INVALID_ARGUMENT;
    alea_cluster_status_t status = alea_cluster_agree(cluster, local);
    if (status != ALEA_CLUSTER_OK) return status;
    const size_t nx = (size_t)cfg->nx, ny = (size_t)cfg->ny;
    const size_t nz = (size_t)cfg->nz;
    local = nx > SIZE_MAX / ny || nx * ny > SIZE_MAX / nz ||
        nx * ny * nz > INT_MAX ||
        nx >= SIZE_MAX / sizeof(double) ||
        ny >= SIZE_MAX / sizeof(double) ||
        nz >= SIZE_MAX / sizeof(double)
        ? ALEA_CLUSTER_OUTPUT_LIMIT : ALEA_CLUSTER_OK;
    status = alea_cluster_agree(cluster, local);
    if (status != ALEA_CLUSTER_OK) return status;
    uint64_t samples_per_voxel = cfg->sampling_mode == ALEA_MESH_SAMPLE_CENTER
        ? 1u : cfg->sampling_mode == ALEA_MESH_SAMPLE_CORNERS ? 8u : 0u;
    if (cfg->sampling_mode == ALEA_MESH_SAMPLE_SUBCELL ||
        cfg->sampling_mode == ALEA_MESH_SAMPLE_STRATIFIED ||
        cfg->sampling_mode == ALEA_MESH_SAMPLE_ADAPTIVE) {
        uint64_t n = cfg->subsamples_per_axis > 0
            ? (uint64_t)cfg->subsamples_per_axis : 0;
        if (n && n <= UINT64_MAX / n && n * n <= UINT64_MAX / n)
            samples_per_voxel = n * n * n;
    }
    if (cfg->sampling_mode == ALEA_MESH_SAMPLE_RAY) samples_per_voxel = 1;
    local = !samples_per_voxel ||
        samples_per_voxel > cfg->max_samples_per_voxel ||
        (cfg->max_total_samples &&
         (samples_per_voxel > UINT64_MAX / (nx * ny * nz) ||
          samples_per_voxel * (nx * ny * nz) > cfg->max_total_samples))
        ? ALEA_CLUSTER_INVALID_ARGUMENT : ALEA_CLUSTER_OK;
    status = alea_cluster_agree(cluster, local);
    if (status != ALEA_CLUSTER_OK) return status;
    uint64_t model = alea_cluster_system_fingerprint(sys);
    local = model ? ALEA_CLUSTER_OK : ALEA_CLUSTER_COMPUTE_ERROR;
    status = alea_cluster_agree(cluster, local);
    if (status != ALEA_CLUSTER_OK) return status;
    int match = alea_cluster_fingerprints_match(cluster, model);
    if (match < 0) return ALEA_CLUSTER_BACKEND_ERROR;
    if (!match) return ALEA_CLUSTER_MODEL_MISMATCH;
    match = alea_cluster_fingerprints_match(cluster, config_fingerprint(cfg));
    if (match < 0) return ALEA_CLUSTER_BACKEND_ERROR;
    if (!match) return ALEA_CLUSTER_INVALID_ARGUMENT;

    const int all_bounds_zero = cfg->x_min == 0 && cfg->x_max == 0 &&
        cfg->y_min == 0 && cfg->y_max == 0 &&
        cfg->z_min == 0 && cfg->z_max == 0;
    const int auto_bounds = cfg->bounds_mode == ALEA_MESH_BOUNDS_AUTO ||
        (cfg->bounds_mode == ALEA_MESH_BOUNDS_LEGACY && all_bounds_zero);
    alea_mesh_config_t resolved = *cfg;
    if (auto_bounds && (!cfg->x_nodes || !cfg->y_nodes || !cfg->z_nodes)) {
        double bounds[6];
        local = alea_mesh_cluster_auto_bounds(sys, cfg, bounds) == 0
            ? ALEA_CLUSTER_OK : ALEA_CLUSTER_COMPUTE_ERROR;
        status = alea_cluster_agree(cluster, local);
        if (status != ALEA_CLUSTER_OK) return status;
        if (!cfg->x_nodes) { resolved.x_min = bounds[0]; resolved.x_max = bounds[1]; }
        if (!cfg->y_nodes) { resolved.y_min = bounds[2]; resolved.y_max = bounds[3]; }
        if (!cfg->z_nodes) { resolved.z_min = bounds[4]; resolved.z_max = bounds[5]; }
    }
    local = (!cfg->x_nodes && !(resolved.x_max > resolved.x_min)) ||
            (!cfg->y_nodes && !(resolved.y_max > resolved.y_min)) ||
            (!cfg->z_nodes && !(resolved.z_max > resolved.z_min))
        ? ALEA_CLUSTER_INVALID_ARGUMENT : ALEA_CLUSTER_OK;
    status = alea_cluster_agree(cluster, local);
    if (status != ALEA_CLUSTER_OK) return status;

    const size_t ranks = (size_t)cluster->size;
    const size_t slices = slab_count(nz, (size_t)cluster->rank, ranks);
    const size_t first = slab_first(nz, (size_t)cluster->rank, ranks);
    const size_t nxy = nx * ny, cells = nxy * nz;
    double* z_nodes = malloc((nz + 1) * sizeof(double));
    local = z_nodes ? ALEA_CLUSTER_OK : ALEA_CLUSTER_OUT_OF_MEMORY;
    status = alea_cluster_agree(cluster, local);
    if (status != ALEA_CLUSTER_OK) { free(z_nodes); return status; }
    if (cfg->z_nodes)
        memcpy(z_nodes, cfg->z_nodes, (nz + 1) * sizeof(double));
    else {
        double step = (resolved.z_max - resolved.z_min) / cfg->nz;
        for (size_t k = 0; k <= nz; ++k)
            z_nodes[k] = resolved.z_min + (double)k * step;
    }
    alea_mesh_result_t* partial = NULL;
    alea_mesh_config_t part_cfg = resolved;
    if (slices) {
        part_cfg.nz = (int)slices;
        part_cfg.z_nodes = z_nodes + first;
        part_cfg.z_min = z_nodes[first];
        part_cfg.z_max = z_nodes[first + slices];
        part_cfg.bounds_mode = ALEA_MESH_BOUNDS_EXPLICIT;
    }
    status = sample_local_slab(cluster, sys, &part_cfg, slices, first,
        z_nodes[0], z_nodes[nz],
        cfg->sampling_mode == ALEA_MESH_SAMPLE_ADAPTIVE &&
            cfg->max_total_samples != 0,
        &partial);
    if (status != ALEA_CLUSTER_OK) goto cleanup;
    if (partial) {
        for (size_t i = 0; i < partial->fraction_count; ++i) {
            int id = partial->fractions[i].material_id;
            double fraction = partial->fractions[i].fraction;
            memset(&partial->fractions[i], 0, sizeof(partial->fractions[i]));
            partial->fractions[i].material_id = id;
            partial->fractions[i].fraction = fraction;
        }
        for (size_t i = 0; i < partial->cell_fraction_count; ++i) {
            int cell = partial->cell_fractions[i].cell_id;
            int material = partial->cell_fractions[i].material_id;
            double fraction = partial->cell_fractions[i].fraction;
            memset(&partial->cell_fractions[i], 0,
                   sizeof(partial->cell_fractions[i]));
            partial->cell_fractions[i].cell_id = cell;
            partial->cell_fractions[i].material_id = material;
            partial->cell_fractions[i].fraction = fraction;
        }
    }

    uint64_t* fractions_by_rank = cluster->rank == 0
        ? calloc(ranks, sizeof(uint64_t)) : NULL;
    uint64_t* cells_by_rank = cluster->rank == 0
        ? calloc(ranks, sizeof(uint64_t)) : NULL;
    uint64_t* materials_by_rank = cluster->rank == 0
        ? calloc(ranks, sizeof(uint64_t)) : NULL;
    uint64_t* mixed_by_rank = cluster->rank == 0
        ? calloc(ranks, sizeof(uint64_t)) : NULL;
    size_t* rank_bytes = cluster->rank == 0
        ? calloc(ranks, sizeof(size_t)) : NULL;
    alea_mesh_result_t* result = cluster->rank == 0
        ? calloc(1, sizeof(*result)) : NULL;
    local = cluster->rank != 0 ||
        (fractions_by_rank && cells_by_rank && materials_by_rank &&
         mixed_by_rank && rank_bytes && result)
        ? ALEA_CLUSTER_OK : ALEA_CLUSTER_OUT_OF_MEMORY;
    status = alea_cluster_agree(cluster, local);
    if (status != ALEA_CLUSTER_OK) goto cleanup_result;

#define GATHER_COUNT(value, target) do { \
    if (alea_cluster_backend_gather_u64(cluster->backend, (uint64_t)(value), \
                                        (target), 0)) { \
        cluster->usable = 0; status = ALEA_CLUSTER_BACKEND_ERROR; \
        goto cleanup_result; \
    } \
} while (0)
    GATHER_COUNT(partial ? partial->fraction_count : 0, fractions_by_rank);
    GATHER_COUNT(partial ? partial->cell_fraction_count : 0, cells_by_rank);
    GATHER_COUNT(partial ? partial->num_materials : 0, materials_by_rank);
    GATHER_COUNT(partial ? partial->mixed_count : 0, mixed_by_rank);
#undef GATHER_COUNT
    size_t material_count = 0;
    local = ALEA_CLUSTER_OK;
    if (cluster->rank == 0) {
        result->nx = cfg->nx; result->ny = cfg->ny; result->nz = cfg->nz;
        result->fields = cfg->fields;
        result->bounds_source = cfg->x_nodes && cfg->y_nodes && cfg->z_nodes
            ? ALEA_MESH_BOUNDS_SOURCE_CUSTOM_NODES
            : auto_bounds ? ALEA_MESH_BOUNDS_SOURCE_INFERRED_ROOT_AABB
                          : ALEA_MESH_BOUNDS_SOURCE_EXPLICIT;
        result->bounds_padding = auto_bounds ? cfg->auto_pad : 0.0;
        result->sampling_mode = cfg->sampling_mode;
        result->sampling_seed = cfg->sampling_seed;
        result->target_error = cfg->target_error;
        for (size_t rank = 0; rank < ranks; ++rank) {
            if (fractions_by_rank[rank] > UINT32_MAX - result->fraction_count ||
                cells_by_rank[rank] > UINT32_MAX - result->cell_fraction_count ||
                materials_by_rank[rank] > SIZE_MAX - material_count ||
                mixed_by_rank[rank] >
                    (uint64_t)(INT_MAX - result->mixed_count)) {
                local = ALEA_CLUSTER_OUTPUT_LIMIT;
                break;
            }
            result->fraction_count += (size_t)fractions_by_rank[rank];
            result->cell_fraction_count += (size_t)cells_by_rank[rank];
            material_count += (size_t)materials_by_rank[rank];
            result->mixed_count += (int)mixed_by_rank[rank];
        }
        if (local == ALEA_CLUSTER_OK) {
            result->x_nodes = malloc((nx + 1) * sizeof(double));
            result->y_nodes = malloc((ny + 1) * sizeof(double));
            if (!result->x_nodes || !result->y_nodes)
                local = ALEA_CLUSTER_OUT_OF_MEMORY;
            else {
                memcpy(result->x_nodes, partial->x_nodes,
                       (nx + 1) * sizeof(double));
                memcpy(result->y_nodes, partial->y_nodes,
                       (ny + 1) * sizeof(double));
            }
        }
#define ALLOC_FIXED(flag, member, type) \
    if (local == ALEA_CLUSTER_OK && (cfg->fields & (flag)) && \
        !(result->member = malloc(cells * sizeof(type)))) \
        local = ALEA_CLUSTER_OUT_OF_MEMORY
        ALLOC_FIXED(ALEA_MESH_FIELD_MATERIAL_ID, material_ids, int);
        ALLOC_FIXED(ALEA_MESH_FIELD_CELL_ID, cell_ids, int);
        ALLOC_FIXED(ALEA_MESH_FIELD_MIXED_FLAG, mixed_flags, unsigned char);
        ALLOC_FIXED(ALEA_MESH_FIELD_DOMINANT_FRACTION, dominant_fractions, double);
        ALLOC_FIXED(ALEA_MESH_FIELD_ESTIMATED_ERROR, estimated_errors, double);
        ALLOC_FIXED(ALEA_MESH_FIELD_SAMPLE_COUNT, sample_counts, uint32_t);
        ALLOC_FIXED(ALEA_MESH_FIELD_TIE_FLAG, tie_flags, uint8_t);
        ALLOC_FIXED(ALEA_MESH_FIELD_REFINEMENT_FLAG, refinement_flags, uint8_t);
        ALLOC_FIXED(ALEA_MESH_FIELD_SAMPLED_FRACTIONS,
                    fraction_spans, alea_mesh_fraction_span_t);
        ALLOC_FIXED(ALEA_MESH_FIELD_CELL_FRACTIONS,
                    cell_fraction_spans, alea_mesh_fraction_span_t);
#undef ALLOC_FIXED
        if (local == ALEA_CLUSTER_OK) {
            result->fractions = malloc((result->fraction_count
                ? result->fraction_count : 1) * sizeof(*result->fractions));
            result->cell_fractions = malloc((result->cell_fraction_count
                ? result->cell_fraction_count : 1) *
                sizeof(*result->cell_fractions));
            result->unique_materials = malloc((material_count
                ? material_count : 1) * sizeof(int));
            if (!result->fractions || !result->cell_fractions ||
                !result->unique_materials)
                local = ALEA_CLUSTER_OUT_OF_MEMORY;
        }
    }
    status = alea_cluster_agree(cluster, local);
    if (status != ALEA_CLUSTER_OK) goto cleanup_result;

#define GATHER_FIELD(member, type, local_count, count_for_rank, target) do { \
    if (cluster->rank == 0) \
        for (size_t rank = 0; rank < ranks; ++rank) \
            rank_bytes[rank] = (count_for_rank) * sizeof(type); \
    if (alea_cluster_backend_gather_bytes(cluster->backend, \
            partial ? partial->member : NULL, (local_count) * sizeof(type), \
            cluster->rank == 0 ? (target) : NULL, rank_bytes, 0)) { \
        cluster->usable = 0; status = ALEA_CLUSTER_BACKEND_ERROR; \
        goto cleanup_result; \
    } \
} while (0)
#define FIXED(flag, member, type) do { \
    if (cfg->fields & (flag)) \
        GATHER_FIELD(member, type, slices * nxy, \
            slab_count(nz, rank, ranks) * nxy, result->member); \
} while (0)
    FIXED(ALEA_MESH_FIELD_MATERIAL_ID, material_ids, int);
    FIXED(ALEA_MESH_FIELD_CELL_ID, cell_ids, int);
    FIXED(ALEA_MESH_FIELD_MIXED_FLAG, mixed_flags, unsigned char);
    FIXED(ALEA_MESH_FIELD_DOMINANT_FRACTION, dominant_fractions, double);
    FIXED(ALEA_MESH_FIELD_ESTIMATED_ERROR, estimated_errors, double);
    FIXED(ALEA_MESH_FIELD_SAMPLE_COUNT, sample_counts, uint32_t);
    FIXED(ALEA_MESH_FIELD_TIE_FLAG, tie_flags, uint8_t);
    FIXED(ALEA_MESH_FIELD_REFINEMENT_FLAG, refinement_flags, uint8_t);
    FIXED(ALEA_MESH_FIELD_SAMPLED_FRACTIONS,
          fraction_spans, alea_mesh_fraction_span_t);
    FIXED(ALEA_MESH_FIELD_CELL_FRACTIONS,
          cell_fraction_spans, alea_mesh_fraction_span_t);
#undef FIXED
    if (cfg->fields & ALEA_MESH_FIELD_SAMPLED_FRACTIONS) {
        GATHER_FIELD(fractions, alea_mesh_material_fraction_t,
            partial ? partial->fraction_count : 0,
            (size_t)fractions_by_rank[rank], result->fractions);
        if (cluster->rank == 0) {
            size_t offset = 0;
            for (size_t rank = 0; rank < ranks; ++rank) {
                size_t first_cell = slab_first(nz, rank, ranks) * nxy;
                size_t count = slab_count(nz, rank, ranks) * nxy;
                for (size_t v = 0; v < count; ++v)
                    result->fraction_spans[first_cell + v].offset +=
                        (uint32_t)offset;
                offset += (size_t)fractions_by_rank[rank];
            }
        }
    }
    if (cfg->fields & ALEA_MESH_FIELD_CELL_FRACTIONS) {
        GATHER_FIELD(cell_fractions, alea_mesh_cell_fraction_t,
            partial ? partial->cell_fraction_count : 0,
            (size_t)cells_by_rank[rank], result->cell_fractions);
        if (cluster->rank == 0) {
            size_t offset = 0;
            for (size_t rank = 0; rank < ranks; ++rank) {
                size_t first_cell = slab_first(nz, rank, ranks) * nxy;
                size_t count = slab_count(nz, rank, ranks) * nxy;
                for (size_t v = 0; v < count; ++v)
                    result->cell_fraction_spans[first_cell + v].offset +=
                        (uint32_t)offset;
                offset += (size_t)cells_by_rank[rank];
            }
        }
    }
    GATHER_FIELD(unique_materials, int,
        partial ? (size_t)partial->num_materials : 0,
        (size_t)materials_by_rank[rank], result->unique_materials);
#undef GATHER_FIELD
    if (cluster->rank == 0) {
        qsort(result->unique_materials, material_count, sizeof(int),
              compare_ints);
        size_t unique = 0;
        for (size_t i = 0; i < material_count; ++i)
            if (!unique || result->unique_materials[i] !=
                           result->unique_materials[unique - 1])
                result->unique_materials[unique++] = result->unique_materials[i];
        result->num_materials = (int)unique;
        result->z_nodes = z_nodes;
        z_nodes = NULL;
        *root_result = result;
        result = NULL;
    }
    status = ALEA_CLUSTER_OK;

cleanup_result:
    free(fractions_by_rank); free(cells_by_rank);
    free(materials_by_rank); free(mixed_by_rank); free(rank_bytes);
    alea_mesh_result_free(result);
cleanup:
    alea_mesh_result_free(partial);
    free(z_nodes);
    return status;
}

static int visit_mesh_slab(const alea_mesh_result_t* mesh, size_t first_z,
        alea_mesh_voxel_visit_fn callback, void* user_data) {
    const size_t nx = (size_t)mesh->nx, ny = (size_t)mesh->ny;
    for (size_t k = 0; k < (size_t)mesh->nz; ++k)
        for (size_t j = 0; j < ny; ++j)
            for (size_t i = 0; i < nx; ++i) {
                const size_t index = (k * ny + j) * nx + i;
                const alea_mesh_fraction_span_t material_span =
                    mesh->fraction_spans[index];
                const alea_mesh_fraction_span_t cell_span =
                    mesh->cell_fraction_spans[index];
                alea_mesh_voxel_sample_t sample = {
                    .i = (int)i, .j = (int)j, .k = (int)(first_z + k),
                    .x_min = mesh->x_nodes[i], .x_max = mesh->x_nodes[i + 1],
                    .y_min = mesh->y_nodes[j], .y_max = mesh->y_nodes[j + 1],
                    .z_min = mesh->z_nodes[k], .z_max = mesh->z_nodes[k + 1],
                    .material_id = mesh->material_ids[index],
                    .cell_id = mesh->cell_ids[index],
                    .mixed = mesh->mixed_flags[index],
                    .tie_flags = mesh->tie_flags[index],
                    .dominant_fraction = mesh->dominant_fractions[index],
                    .estimated_error = mesh->estimated_errors[index],
                    .sample_count = mesh->sample_counts[index],
                    .refinement_flags = mesh->refinement_flags[index],
                    .fractions = material_span.count
                        ? mesh->fractions + material_span.offset : NULL,
                    .fraction_count = material_span.count,
                    .cell_fractions = cell_span.count
                        ? mesh->cell_fractions + cell_span.offset : NULL,
                    .cell_fraction_count = cell_span.count
                };
                if (callback(&sample, user_data)) return -1;
            }
    return 0;
}

static alea_cluster_status_t mesh_stream_root_slabs(
        alea_cluster_t* cluster, const alea_mesh_result_t* local_slab,
        const alea_mesh_config_t* cfg,
        alea_mesh_voxel_visit_fn callback, void* user_data) {
    const size_t ranks = (size_t)cluster->size;
    const size_t nz = (size_t)cfg->nz;
    const size_t nx = (size_t)cfg->nx, ny = (size_t)cfg->ny;
    size_t* rank_bytes = cluster->rank == 0
        ? calloc(ranks, sizeof(size_t)) : NULL;
    alea_cluster_status_t local = cluster->rank != 0 || rank_bytes
        ? ALEA_CLUSTER_OK : ALEA_CLUSTER_OUT_OF_MEMORY;
    alea_cluster_status_t status = alea_cluster_agree(cluster, local);
    if (status != ALEA_CLUSTER_OK) { free(rank_bytes); return status; }
    alea_mesh_result_t* received = NULL;
    for (size_t peer = 0; peer < ranks; ++peer) {
        const size_t slices = slab_count(nz, peer, ranks);
        if (!slices) continue;
        const size_t first = slab_first(nz, peer, ranks);
        if (peer == 0) {
            local = cluster->rank == 0 &&
                visit_mesh_slab(local_slab, first, callback, user_data)
                    ? ALEA_CLUSTER_INTERRUPTED : ALEA_CLUSTER_OK;
            status = alea_cluster_agree(cluster, local);
            if (status != ALEA_CLUSTER_OK) goto cleanup;
            continue;
        }
        uint64_t counts[2] = {0, 0};
        if ((size_t)cluster->rank == peer) {
            counts[0] = (uint64_t)local_slab->fraction_count;
            counts[1] = (uint64_t)local_slab->cell_fraction_count;
        }
        if (alea_cluster_backend_broadcast_bytes(cluster->backend, counts,
                sizeof(counts), (int)peer)) {
            cluster->usable = 0;
            status = ALEA_CLUSTER_BACKEND_ERROR;
            goto cleanup;
        }
        const size_t cells = nx * ny * slices;
        local = ALEA_CLUSTER_OK;
        if (cluster->rank == 0) {
            if (counts[0] > UINT32_MAX || counts[1] > UINT32_MAX) {
                local = ALEA_CLUSTER_OUTPUT_LIMIT;
            } else {
                received = calloc(1, sizeof(*received));
                if (!received) local = ALEA_CLUSTER_OUT_OF_MEMORY;
            }
            if (received) {
                received->nx = cfg->nx;
                received->ny = cfg->ny;
                received->nz = (int)slices;
                received->fields = MESH_FIELDS;
                received->fraction_count = (size_t)counts[0];
                received->cell_fraction_count = (size_t)counts[1];
#define ALLOC(member, count, type) do { \
    size_t n = (count); \
    if (n > SIZE_MAX / sizeof(type)) local = ALEA_CLUSTER_OUTPUT_LIMIT; \
    else if (!(received->member = malloc((n ? n : 1) * sizeof(type)))) \
        local = ALEA_CLUSTER_OUT_OF_MEMORY; \
} while (0)
                ALLOC(x_nodes, nx + 1, double);
                ALLOC(y_nodes, ny + 1, double);
                ALLOC(z_nodes, slices + 1, double);
                ALLOC(material_ids, cells, int);
                ALLOC(cell_ids, cells, int);
                ALLOC(mixed_flags, cells, unsigned char);
                ALLOC(dominant_fractions, cells, double);
                ALLOC(estimated_errors, cells, double);
                ALLOC(sample_counts, cells, uint32_t);
                ALLOC(tie_flags, cells, uint8_t);
                ALLOC(refinement_flags, cells, uint8_t);
                ALLOC(fraction_spans, cells, alea_mesh_fraction_span_t);
                ALLOC(cell_fraction_spans, cells, alea_mesh_fraction_span_t);
                ALLOC(fractions, received->fraction_count,
                    alea_mesh_material_fraction_t);
                ALLOC(cell_fractions, received->cell_fraction_count,
                    alea_mesh_cell_fraction_t);
#undef ALLOC
            }
        }
        status = alea_cluster_agree(cluster, local);
        if (status != ALEA_CLUSTER_OK) goto cleanup;
#define TRANSFER(member, count, type) do { \
    const size_t bytes = (count) * sizeof(type); \
    if (cluster->rank == 0) rank_bytes[peer] = bytes; \
    const void* source = (size_t)cluster->rank == peer \
        ? local_slab->member : NULL; \
    if (alea_cluster_backend_gather_bytes(cluster->backend, source, \
            (size_t)cluster->rank == peer ? bytes : 0, \
            cluster->rank == 0 ? received->member : NULL, rank_bytes, 0)) { \
        cluster->usable = 0; status = ALEA_CLUSTER_BACKEND_ERROR; \
        goto cleanup; \
    } \
    if (cluster->rank == 0) rank_bytes[peer] = 0; \
} while (0)
        TRANSFER(x_nodes, nx + 1, double);
        TRANSFER(y_nodes, ny + 1, double);
        TRANSFER(z_nodes, slices + 1, double);
        TRANSFER(material_ids, cells, int);
        TRANSFER(cell_ids, cells, int);
        TRANSFER(mixed_flags, cells, unsigned char);
        TRANSFER(dominant_fractions, cells, double);
        TRANSFER(estimated_errors, cells, double);
        TRANSFER(sample_counts, cells, uint32_t);
        TRANSFER(tie_flags, cells, uint8_t);
        TRANSFER(refinement_flags, cells, uint8_t);
        TRANSFER(fraction_spans, cells, alea_mesh_fraction_span_t);
        TRANSFER(cell_fraction_spans, cells, alea_mesh_fraction_span_t);
        TRANSFER(fractions, (size_t)counts[0],
            alea_mesh_material_fraction_t);
        TRANSFER(cell_fractions, (size_t)counts[1],
            alea_mesh_cell_fraction_t);
#undef TRANSFER
        local = cluster->rank == 0 &&
            visit_mesh_slab(received, first, callback, user_data)
                ? ALEA_CLUSTER_INTERRUPTED : ALEA_CLUSTER_OK;
        status = alea_cluster_agree(cluster, local);
        alea_mesh_result_free(received);
        received = NULL;
        if (status != ALEA_CLUSTER_OK) goto cleanup;
    }
cleanup:
    alea_mesh_result_free(received);
    free(rank_bytes);
    return status;
}

static alea_cluster_status_t mesh_sample_shards_impl(
        alea_cluster_t* cluster, alea_system_t* sys,
        const alea_mesh_config_t* cfg,
        alea_cluster_mesh_slab_callback_t callback, void* user_data,
        alea_mesh_voxel_visit_fn root_visitor, int root_stream) {
    if (!cluster) return ALEA_CLUSTER_INVALID_ARGUMENT;
    alea_cluster_status_t local = ALEA_CLUSTER_OK;
    if (!sys || !cfg ||
        (root_stream ? (cluster->rank == 0 && !root_visitor) : !callback) ||
        cfg->nx <= 0 || cfg->ny <= 0 ||
        cfg->nz <= 0 || cfg->visit || cfg->progress ||
        (cfg->fields & ~MESH_FIELDS) ||
        (cfg->sampling_mode != ALEA_MESH_SAMPLE_CENTER &&
         cfg->sampling_mode != ALEA_MESH_SAMPLE_CORNERS &&
         cfg->sampling_mode != ALEA_MESH_SAMPLE_SUBCELL &&
         cfg->sampling_mode != ALEA_MESH_SAMPLE_STRATIFIED &&
         cfg->sampling_mode != ALEA_MESH_SAMPLE_ADAPTIVE &&
         cfg->sampling_mode != ALEA_MESH_SAMPLE_RAY) ||
        (cfg->sampling_mode == ALEA_MESH_SAMPLE_RAY &&
         (!(cfg->ray_directions & ALEA_MESH_RAY_XYZ) ||
          cfg->max_total_samples != 0 ||
          (cfg->ray_origin_mode == ALEA_MESH_RAY_ORIGINS_CUSTOM &&
           (uint64_t)cfg->ray_point_count * 2u * sizeof(double) >
               SIZE_MAX))))
        local = ALEA_CLUSTER_INVALID_ARGUMENT;
    alea_cluster_status_t status = alea_cluster_agree(cluster, local);
    if (status != ALEA_CLUSTER_OK) return status;
    const size_t nx = (size_t)cfg->nx, ny = (size_t)cfg->ny;
    const size_t nz = (size_t)cfg->nz;
    local = nx > SIZE_MAX / ny || nx * ny > SIZE_MAX / nz ||
        nx >= SIZE_MAX / sizeof(double) ||
        ny >= SIZE_MAX / sizeof(double) ||
        nz >= SIZE_MAX / sizeof(double)
        ? ALEA_CLUSTER_OUTPUT_LIMIT : ALEA_CLUSTER_OK;
    status = alea_cluster_agree(cluster, local);
    if (status != ALEA_CLUSTER_OK) return status;
    uint64_t per_voxel = cfg->sampling_mode == ALEA_MESH_SAMPLE_CENTER
        ? 1u : cfg->sampling_mode == ALEA_MESH_SAMPLE_CORNERS ? 8u : 0u;
    if (!per_voxel) {
        uint64_t n = cfg->subsamples_per_axis > 0
            ? (uint64_t)cfg->subsamples_per_axis : 0;
        if (n && n <= UINT64_MAX / n && n * n <= UINT64_MAX / n)
            per_voxel = n * n * n;
    }
    if (cfg->sampling_mode == ALEA_MESH_SAMPLE_RAY) per_voxel = 1;
    local = !per_voxel || per_voxel > cfg->max_samples_per_voxel ||
        (cfg->max_total_samples &&
         (per_voxel > UINT64_MAX / (nx * ny * nz) ||
          per_voxel * (nx * ny * nz) > cfg->max_total_samples))
        ? ALEA_CLUSTER_INVALID_ARGUMENT : ALEA_CLUSTER_OK;
    status = alea_cluster_agree(cluster, local);
    if (status != ALEA_CLUSTER_OK) return status;
    uint64_t model = alea_cluster_system_fingerprint(sys);
    local = model ? ALEA_CLUSTER_OK : ALEA_CLUSTER_COMPUTE_ERROR;
    status = alea_cluster_agree(cluster, local);
    if (status != ALEA_CLUSTER_OK) return status;
    int match = alea_cluster_fingerprints_match(cluster, model);
    if (match < 0) return ALEA_CLUSTER_BACKEND_ERROR;
    if (!match) return ALEA_CLUSTER_MODEL_MISMATCH;
    match = alea_cluster_fingerprints_match(cluster, config_fingerprint(cfg));
    if (match < 0) return ALEA_CLUSTER_BACKEND_ERROR;
    if (!match) return ALEA_CLUSTER_INVALID_ARGUMENT;

    const int all_zero = cfg->x_min == 0 && cfg->x_max == 0 &&
        cfg->y_min == 0 && cfg->y_max == 0 &&
        cfg->z_min == 0 && cfg->z_max == 0;
    const int inferred = cfg->bounds_mode == ALEA_MESH_BOUNDS_AUTO ||
        (cfg->bounds_mode == ALEA_MESH_BOUNDS_LEGACY && all_zero);
    alea_mesh_config_t resolved = *cfg;
    if (inferred && (!cfg->x_nodes || !cfg->y_nodes || !cfg->z_nodes)) {
        double bounds[6];
        local = alea_mesh_cluster_auto_bounds(sys, cfg, bounds) == 0
            ? ALEA_CLUSTER_OK : ALEA_CLUSTER_COMPUTE_ERROR;
        status = alea_cluster_agree(cluster, local);
        if (status != ALEA_CLUSTER_OK) return status;
        if (!cfg->x_nodes) { resolved.x_min = bounds[0]; resolved.x_max = bounds[1]; }
        if (!cfg->y_nodes) { resolved.y_min = bounds[2]; resolved.y_max = bounds[3]; }
        if (!cfg->z_nodes) { resolved.z_min = bounds[4]; resolved.z_max = bounds[5]; }
    }
    local = (!cfg->x_nodes && !(resolved.x_max > resolved.x_min)) ||
            (!cfg->y_nodes && !(resolved.y_max > resolved.y_min)) ||
            (!cfg->z_nodes && !(resolved.z_max > resolved.z_min))
        ? ALEA_CLUSTER_INVALID_ARGUMENT : ALEA_CLUSTER_OK;
    status = alea_cluster_agree(cluster, local);
    if (status != ALEA_CLUSTER_OK) return status;

    const size_t ranks = (size_t)cluster->size;
    const size_t rank = (size_t)cluster->rank;
    const size_t slices = slab_count(nz, rank, ranks);
    const size_t first = slab_first(nz, rank, ranks);
    double* local_nodes = slices
        ? malloc((slices + 1) * sizeof(double)) : NULL;
    local = !slices || local_nodes ? ALEA_CLUSTER_OK
                                   : ALEA_CLUSTER_OUT_OF_MEMORY;
    status = alea_cluster_agree(cluster, local);
    if (status != ALEA_CLUSTER_OK) { free(local_nodes); return status; }
    alea_mesh_result_t* slab = NULL;
    alea_mesh_config_t part = resolved;
    double global_z_min = cfg->z_nodes ? cfg->z_nodes[0] : resolved.z_min;
    double global_z_max = cfg->z_nodes
        ? cfg->z_nodes[nz]
        : resolved.z_min + (double)nz *
            ((resolved.z_max - resolved.z_min) / cfg->nz);
    if (slices) {
        if (cfg->z_nodes)
            memcpy(local_nodes, cfg->z_nodes + first,
                   (slices + 1) * sizeof(double));
        else {
            const double step = (resolved.z_max - resolved.z_min) / cfg->nz;
            for (size_t k = 0; k <= slices; ++k)
                local_nodes[k] = resolved.z_min + (double)(first + k) * step;
        }
        part.nz = (int)slices;
        part.z_nodes = local_nodes;
        part.z_min = local_nodes[0];
        part.z_max = local_nodes[slices];
        part.bounds_mode = ALEA_MESH_BOUNDS_EXPLICIT;
    }
    status = sample_local_slab(cluster, sys, &part, slices, first,
        global_z_min, global_z_max,
        cfg->sampling_mode == ALEA_MESH_SAMPLE_ADAPTIVE &&
            cfg->max_total_samples != 0,
        &slab);
    if (slab) {
        slab->bounds_source = cfg->x_nodes && cfg->y_nodes && cfg->z_nodes
            ? ALEA_MESH_BOUNDS_SOURCE_CUSTOM_NODES
            : inferred ? ALEA_MESH_BOUNDS_SOURCE_INFERRED_ROOT_AABB
                       : ALEA_MESH_BOUNDS_SOURCE_EXPLICIT;
        slab->bounds_padding = inferred ? cfg->auto_pad : 0.0;
    }
    if (status == ALEA_CLUSTER_OK) {
        if (root_stream)
            status = mesh_stream_root_slabs(cluster, slab, cfg,
                root_visitor, user_data);
        else {
            local = slices && callback(first, slab, user_data)
                ? ALEA_CLUSTER_INTERRUPTED : ALEA_CLUSTER_OK;
            status = alea_cluster_agree(cluster, local);
        }
    }
    alea_mesh_result_free(slab);
    free(local_nodes);
    return status;
}

alea_cluster_status_t alea_cluster_mesh_sample_shards(
        alea_cluster_t* cluster, alea_system_t* sys,
        const alea_mesh_config_t* cfg,
        alea_cluster_mesh_slab_callback_t callback, void* user_data) {
    return mesh_sample_shards_impl(cluster, sys, cfg, callback, user_data,
        NULL, 0);
}

alea_cluster_status_t alea_cluster_mesh_stream_root(
        alea_cluster_t* cluster, alea_system_t* sys,
        const alea_mesh_config_t* cfg,
        alea_mesh_voxel_visit_fn callback, void* user_data) {
    if (!cluster) return ALEA_CLUSTER_INVALID_ARGUMENT;
    alea_cluster_status_t local = cfg
        ? ALEA_CLUSTER_OK : ALEA_CLUSTER_INVALID_ARGUMENT;
    alea_cluster_status_t status = alea_cluster_agree(cluster, local);
    if (status != ALEA_CLUSTER_OK) return status;
    alea_mesh_config_t full_config = *cfg;
    full_config.fields = MESH_FIELDS;
    return mesh_sample_shards_impl(cluster, sys, &full_config,
        NULL, user_data, callback, 1);
}

alea_cluster_status_t alea_cluster_mesh_visit_root(
        alea_cluster_t* cluster, alea_system_t* sys,
        const alea_mesh_config_t* cfg,
        alea_mesh_voxel_visit_fn callback, void* user_data) {
    if (!cluster) return ALEA_CLUSTER_INVALID_ARGUMENT;
    alea_cluster_status_t local = cfg &&
        (cluster->rank != 0 || callback)
        ? ALEA_CLUSTER_OK : ALEA_CLUSTER_INVALID_ARGUMENT;
    alea_cluster_status_t status = alea_cluster_agree(cluster, local);
    if (status != ALEA_CLUSTER_OK) return status;
    alea_mesh_config_t full_config = *cfg;
    full_config.fields = MESH_FIELDS;
    alea_mesh_result_t* result = NULL;
    status = alea_cluster_mesh_sample(cluster, sys, &full_config,
        cluster->rank == 0 ? &result : NULL);
    if (status != ALEA_CLUSTER_OK) return status;
    local = ALEA_CLUSTER_OK;
    if (cluster->rank == 0 &&
        visit_mesh_slab(result, 0, callback, user_data))
        local = ALEA_CLUSTER_INTERRUPTED;
    status = alea_cluster_agree(cluster, local);
    alea_mesh_result_free(result);
    return status;
}
