// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

#include "cluster_internal.h"
#include "render/render_tiles_internal.h"

#include <limits.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

static uint64_t render_hash(uint64_t value, const void* data, size_t length) {
    const unsigned char* bytes = data;
    for (size_t i = 0; i < length; ++i) {
        value ^= bytes[i];
        value *= UINT64_C(1099511628211);
    }
    return value;
}

static uint64_t render_fingerprint(const render_config_t* cfg,
                                   const render_camera_t* cam) {
    uint64_t hash = UINT64_C(1469598103934665603);
#define HASH(value) (hash = render_hash(hash, &(value), sizeof(value)))
#define HASH_ARRAY(value) (hash = render_hash(hash, (value), sizeof(value)))
    HASH(cfg->width); HASH(cfg->height); HASH_ARRAY(cfg->background);
    HASH_ARRAY(cfg->eye); HASH_ARRAY(cfg->target); HASH_ARRAY(cfg->up);
    HASH(cfg->fov); HASH(cfg->ortho_height);
    HASH(cfg->eye_set); HASH(cfg->target_set);
    HASH(cfg->num_clips); HASH(cfg->clip_mode);
    for (int i = 0; i < cfg->num_clips; ++i) {
        HASH_ARRAY(cfg->clips[i].normal); HASH(cfg->clips[i].d);
    }
    const render_id_filter_t* filters[2] = {
        &cfg->material_filter, &cfg->cell_filter
    };
    for (int i = 0; i < 2; ++i) {
        HASH(filters[i]->mode); HASH(filters[i]->count);
        hash = render_hash(hash, filters[i]->ids,
                           filters[i]->count * sizeof(int));
    }
    HASH(cfg->color_mode); HASH(cfg->render_mode);
    HASH_ARRAY(cfg->light_dir); HASH(cfg->light_set);
    HASH(cfg->shadows); HASH(cfg->edges);
    HASH(cfg->ambient); HASH(cfg->diffuse); HASH(cfg->specular);
    HASH(cfg->shininess); HASH(cfg->cross_section_tint);
    HASH(cfg->aa_samples); HASH(cfg->aa_adaptive); HASH(cfg->tile_size);
    HASH(cfg->xray_density_scale); HASH(cfg->num_custom_colors);
    HASH(cfg->custom_colors_sorted);
    for (int i = 0; i < cfg->num_custom_colors; ++i) {
        HASH(cfg->custom_colors[i].id);
        HASH(cfg->custom_colors[i].r);
        HASH(cfg->custom_colors[i].g);
        HASH(cfg->custom_colors[i].b);
    }
    HASH_ARRAY(cam->eye); HASH_ARRAY(cam->target); HASH_ARRAY(cam->up);
    HASH(cam->fov_degrees); HASH(cam->ortho_height);
    HASH_ARRAY(cam->forward); HASH_ARRAY(cam->right);
    HASH_ARRAY(cam->true_up); HASH(cam->half_height); HASH(cam->half_width);
    HASH(cam->auto_fit);
#undef HASH_ARRAY
#undef HASH
    return hash ? hash : 1;
}

static size_t rank_tiles(size_t total, size_t rank, size_t ranks) {
    return total / ranks + (rank < total % ranks);
}

static size_t rank_first_tile(size_t total, size_t rank, size_t ranks) {
    return rank * (total / ranks) + (rank < total % ranks ? rank : total % ranks);
}

static void copy_tiles(void* target, const void* source,
                       int width, int height, size_t tile,
                       size_t tiles_x, size_t tiles,
                       size_t component_size, size_t components) {
    unsigned char* dst = target;
    const unsigned char* src = source;
    const size_t pixel_bytes = component_size * components;
    for (size_t index = 0; index < tiles; ++index) {
        size_t x0 = index % tiles_x * tile;
        size_t y0 = index / tiles_x * tile;
        size_t tw = tile < (size_t)width - x0 ? tile : (size_t)width - x0;
        size_t th = tile < (size_t)height - y0 ? tile : (size_t)height - y0;
        for (size_t y = 0; y < th; ++y)
            memcpy(dst + ((y0 + y) * (size_t)width + x0) * pixel_bytes,
                   src + (index * tile * tile + y * tile) * pixel_bytes,
                   tw * pixel_bytes);
    }
}

alea_cluster_status_t alea_cluster_render_scene(
        alea_cluster_t* cluster, alea_system_t* sys,
        const render_config_t* cfg, const render_camera_t* cam,
        render_framebuffer_t* root_framebuffer) {
    if (!cluster) return ALEA_CLUSTER_INVALID_ARGUMENT;
    alea_cluster_status_t local = ALEA_CLUSTER_OK;
    if (!sys || !cfg || !cam || cfg->width <= 0 || cfg->height <= 0 ||
        cfg->num_clips < 0 || cfg->num_clips > RENDER_MAX_CLIPS ||
        cfg->num_custom_colors < 0 ||
        (cfg->num_custom_colors && !cfg->custom_colors) ||
        (cfg->material_filter.count && !cfg->material_filter.ids) ||
        (cfg->cell_filter.count && !cfg->cell_filter.ids) ||
        cfg->material_filter.count > SIZE_MAX / sizeof(int) ||
        cfg->cell_filter.count > SIZE_MAX / sizeof(int) ||
        cfg->aa_samples > RENDER_MAX_AA)
        local = ALEA_CLUSTER_INVALID_ARGUMENT;
    int aux = 0;
    if (cluster->rank == 0) {
        if (!cfg || !root_framebuffer || !root_framebuffer->color ||
            !root_framebuffer->cell_id ||
            (cfg && (root_framebuffer->width != cfg->width ||
                     root_framebuffer->height != cfg->height)))
            local = ALEA_CLUSTER_INVALID_ARGUMENT;
        if (root_framebuffer) {
            int count = !!root_framebuffer->material_id +
                        !!root_framebuffer->depth + !!root_framebuffer->normal;
            if (count != 0 && count != 3)
                local = ALEA_CLUSTER_INVALID_ARGUMENT;
            aux = count == 3;
        }
    }
    alea_cluster_status_t status = alea_cluster_agree(cluster, local);
    if (status != ALEA_CLUSTER_OK) return status;
    if (alea_cluster_backend_broadcast_int(cluster->backend, &aux, 0)) {
        cluster->usable = 0;
        return ALEA_CLUSTER_BACKEND_ERROR;
    }
    const size_t tile = cfg->tile_size > 0
        ? (size_t)cfg->tile_size : RENDER_DEFAULT_TILE;
    const size_t tiles_x = ((size_t)cfg->width + tile - 1) / tile;
    const size_t tiles_y = ((size_t)cfg->height + tile - 1) / tile;
    local = tiles_x > INT_MAX || tiles_y > SIZE_MAX / tiles_x ||
        tile > SIZE_MAX / tile ||
        tiles_x * tiles_y > SIZE_MAX / (tile * tile) ||
        tiles_x * tiles_y * tile * tile > PTRDIFF_MAX / (3 * sizeof(float))
        ? ALEA_CLUSTER_OUTPUT_LIMIT : ALEA_CLUSTER_OK;
    status = alea_cluster_agree(cluster, local);
    if (status != ALEA_CLUSTER_OK) return status;
    const size_t total_tiles = tiles_x * tiles_y;
    const size_t tile_pixels = tile * tile;
    const size_t total_pixels = total_tiles * tile_pixels;
    uint64_t model = alea_cluster_system_fingerprint(sys);
    local = model ? ALEA_CLUSTER_OK : ALEA_CLUSTER_COMPUTE_ERROR;
    status = alea_cluster_agree(cluster, local);
    if (status != ALEA_CLUSTER_OK) return status;
    int match = alea_cluster_fingerprints_match(cluster, model);
    if (match < 0) return ALEA_CLUSTER_BACKEND_ERROR;
    if (!match) return ALEA_CLUSTER_MODEL_MISMATCH;
    match = alea_cluster_fingerprints_match(
        cluster, render_fingerprint(cfg, cam));
    if (match < 0) return ALEA_CLUSTER_BACKEND_ERROR;
    if (!match) return ALEA_CLUSTER_INVALID_ARGUMENT;

    const size_t ranks = (size_t)cluster->size;
    const size_t local_tiles = rank_tiles(total_tiles,
        (size_t)cluster->rank, ranks);
    const size_t local_pixels = local_tiles * tile_pixels;
    render_framebuffer_t compact = {0};
    compact.color = calloc((local_pixels ? local_pixels : 1) * 3,
                           sizeof(float));
    compact.cell_id = calloc((local_pixels ? local_pixels : 1), sizeof(int));
    if (aux) {
        compact.material_id = calloc((local_pixels ? local_pixels : 1),
                                     sizeof(int));
        compact.depth = calloc((local_pixels ? local_pixels : 1),
                               sizeof(float));
        compact.normal = calloc((local_pixels ? local_pixels : 1) * 3,
                                sizeof(float));
    }
    unsigned char* received = cluster->rank == 0
        ? malloc(total_pixels * 3 * sizeof(float)) : NULL;
    size_t* bytes_by_rank = cluster->rank == 0
        ? malloc(ranks * sizeof(size_t)) : NULL;
    local = compact.color && compact.cell_id &&
        (!aux || (compact.material_id && compact.depth && compact.normal)) &&
        (cluster->rank != 0 || (received && bytes_by_rank))
        ? ALEA_CLUSTER_OK : ALEA_CLUSTER_OUT_OF_MEMORY;
    status = alea_cluster_agree(cluster, local);
    if (status != ALEA_CLUSTER_OK) goto cleanup;

    local = local_tiles && alea_render_scene_tile_span(sys, cfg, cam,
        rank_first_tile(total_tiles, (size_t)cluster->rank, ranks),
        local_tiles, &compact) != 0
        ? ALEA_CLUSTER_COMPUTE_ERROR : ALEA_CLUSTER_OK;
    if (alea_interrupted()) local = ALEA_CLUSTER_INTERRUPTED;
    status = alea_cluster_agree(cluster, local);
    if (status != ALEA_CLUSTER_OK) goto cleanup;

#define GATHER_RENDER(member, type, components) do { \
    if (cluster->rank == 0) \
        for (size_t rank = 0; rank < ranks; ++rank) \
            bytes_by_rank[rank] = rank_tiles(total_tiles, rank, ranks) * \
                tile_pixels * (components) * sizeof(type); \
    if (alea_cluster_backend_gather_bytes(cluster->backend, compact.member, \
            local_pixels * (components) * sizeof(type), received, \
            bytes_by_rank, 0)) { \
        cluster->usable = 0; \
        status = ALEA_CLUSTER_BACKEND_ERROR; goto cleanup; \
    } \
    if (cluster->rank == 0) copy_tiles(root_framebuffer->member, received, \
        cfg->width, cfg->height, tile, tiles_x, total_tiles, \
        sizeof(type), (components)); \
} while (0)
    GATHER_RENDER(color, float, 3);
    GATHER_RENDER(cell_id, int, 1);
    if (aux) {
        GATHER_RENDER(material_id, int, 1);
        GATHER_RENDER(depth, float, 1);
        GATHER_RENDER(normal, float, 3);
    }
#undef GATHER_RENDER
    if (cluster->rank == 0 && cfg->edges)
        render_edge_darken(root_framebuffer);
    status = ALEA_CLUSTER_OK;

cleanup:
    free(compact.color); free(compact.cell_id);
    free(compact.material_id); free(compact.depth); free(compact.normal);
    free(received); free(bytes_by_rank);
    return status;
}
