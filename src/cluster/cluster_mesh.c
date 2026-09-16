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
        (cfg->sampling_mode == ALEA_MESH_SAMPLE_ADAPTIVE &&
         cfg->max_total_samples != 0) ||
        (cfg->sampling_mode == ALEA_MESH_SAMPLE_RAY &&
         ((cfg->ray_directions & ALEA_MESH_RAY_Z) ||
          !(cfg->ray_directions & (ALEA_MESH_RAY_X | ALEA_MESH_RAY_Y)) ||
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
    if (slices) {
        alea_mesh_config_t part_cfg = resolved;
        part_cfg.nz = (int)slices;
        part_cfg.z_nodes = z_nodes + first;
        part_cfg.z_min = z_nodes[first];
        part_cfg.z_max = z_nodes[first + slices];
        part_cfg.bounds_mode = ALEA_MESH_BOUNDS_EXPLICIT;
        partial = alea_mesh_sample_with_z_offset(sys, &part_cfg, (int)first);
    }
    local = !slices || partial ? ALEA_CLUSTER_OK :
        alea_interrupted() ? ALEA_CLUSTER_INTERRUPTED : ALEA_CLUSTER_COMPUTE_ERROR;
    status = alea_cluster_agree(cluster, local);
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

alea_cluster_status_t alea_cluster_mesh_sample_shards(
        alea_cluster_t* cluster, alea_system_t* sys,
        const alea_mesh_config_t* cfg,
        alea_cluster_mesh_slab_callback_t callback, void* user_data) {
    if (!cluster) return ALEA_CLUSTER_INVALID_ARGUMENT;
    alea_cluster_status_t local = ALEA_CLUSTER_OK;
    if (!sys || !cfg || !callback || cfg->nx <= 0 || cfg->ny <= 0 ||
        cfg->nz <= 0 || cfg->visit || cfg->progress ||
        (cfg->fields & ~MESH_FIELDS) ||
        (cfg->sampling_mode != ALEA_MESH_SAMPLE_CENTER &&
         cfg->sampling_mode != ALEA_MESH_SAMPLE_CORNERS &&
         cfg->sampling_mode != ALEA_MESH_SAMPLE_SUBCELL &&
         cfg->sampling_mode != ALEA_MESH_SAMPLE_STRATIFIED &&
         cfg->sampling_mode != ALEA_MESH_SAMPLE_ADAPTIVE &&
         cfg->sampling_mode != ALEA_MESH_SAMPLE_RAY) ||
        (cfg->sampling_mode == ALEA_MESH_SAMPLE_ADAPTIVE &&
         cfg->max_total_samples != 0) ||
        (cfg->sampling_mode == ALEA_MESH_SAMPLE_RAY &&
         ((cfg->ray_directions & ALEA_MESH_RAY_Z) ||
          !(cfg->ray_directions & (ALEA_MESH_RAY_X | ALEA_MESH_RAY_Y)) ||
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
    if (slices) {
        if (cfg->z_nodes)
            memcpy(local_nodes, cfg->z_nodes + first,
                   (slices + 1) * sizeof(double));
        else {
            const double step = (resolved.z_max - resolved.z_min) / cfg->nz;
            for (size_t k = 0; k <= slices; ++k)
                local_nodes[k] = resolved.z_min + (double)(first + k) * step;
        }
        alea_mesh_config_t part = resolved;
        part.nz = (int)slices;
        part.z_nodes = local_nodes;
        part.z_min = local_nodes[0];
        part.z_max = local_nodes[slices];
        part.bounds_mode = ALEA_MESH_BOUNDS_EXPLICIT;
        slab = alea_mesh_sample_with_z_offset(sys, &part, (int)first);
        if (slab) {
            slab->bounds_source = cfg->x_nodes && cfg->y_nodes && cfg->z_nodes
                ? ALEA_MESH_BOUNDS_SOURCE_CUSTOM_NODES
                : inferred ? ALEA_MESH_BOUNDS_SOURCE_INFERRED_ROOT_AABB
                           : ALEA_MESH_BOUNDS_SOURCE_EXPLICIT;
            slab->bounds_padding = inferred ? cfg->auto_pad : 0.0;
        }
    }
    local = !slices || slab ? ALEA_CLUSTER_OK :
        alea_interrupted() ? ALEA_CLUSTER_INTERRUPTED : ALEA_CLUSTER_COMPUTE_ERROR;
    status = alea_cluster_agree(cluster, local);
    if (status == ALEA_CLUSTER_OK) {
        local = slices && callback(first, slab, user_data)
            ? ALEA_CLUSTER_INTERRUPTED : ALEA_CLUSTER_OK;
        status = alea_cluster_agree(cluster, local);
    }
    alea_mesh_result_free(slab);
    free(local_nodes);
    return status;
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
    if (cluster->rank == 0) {
        const size_t nx = (size_t)result->nx;
        const size_t ny = (size_t)result->ny;
        const size_t nz = (size_t)result->nz;
        for (size_t k = 0; k < nz && local == ALEA_CLUSTER_OK; ++k)
            for (size_t j = 0; j < ny && local == ALEA_CLUSTER_OK; ++j)
                for (size_t i = 0; i < nx; ++i) {
                    const size_t index = (k * ny + j) * nx + i;
                    const alea_mesh_fraction_span_t material_span =
                        result->fraction_spans[index];
                    const alea_mesh_fraction_span_t cell_span =
                        result->cell_fraction_spans[index];
                    alea_mesh_voxel_sample_t sample = {
                        .i = (int)i, .j = (int)j, .k = (int)k,
                        .x_min = result->x_nodes[i],
                        .x_max = result->x_nodes[i + 1],
                        .y_min = result->y_nodes[j],
                        .y_max = result->y_nodes[j + 1],
                        .z_min = result->z_nodes[k],
                        .z_max = result->z_nodes[k + 1],
                        .material_id = result->material_ids[index],
                        .cell_id = result->cell_ids[index],
                        .mixed = result->mixed_flags[index],
                        .tie_flags = result->tie_flags[index],
                        .dominant_fraction =
                            result->dominant_fractions[index],
                        .estimated_error = result->estimated_errors[index],
                        .sample_count = result->sample_counts[index],
                        .refinement_flags = result->refinement_flags[index],
                        .fractions = result->fractions + material_span.offset,
                        .fraction_count = material_span.count,
                        .cell_fractions =
                            result->cell_fractions + cell_span.offset,
                        .cell_fraction_count = cell_span.count
                    };
                    if (callback(&sample, user_data)) {
                        local = ALEA_CLUSTER_INTERRUPTED;
                        break;
                    }
                }
    }
    status = alea_cluster_agree(cluster, local);
    alea_mesh_result_free(result);
    return status;
}
