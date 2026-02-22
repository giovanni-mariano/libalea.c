// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

#include "alea_system.h"
#include "alea_types.h"
#include "alea.h"
#include "core/alea_eval.h"
#include "core/alea_ops.h"
#include "core/alea_universe.h"
#include "core/alea_spatial.h"
#include "util/alea_log.h"
#include "util/compat.h"
#include "raycast/bvh.h"
#include "primitives/primitive_desc.h"
#include "primitives/bbox.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include "util/math.h"

#define INITIAL_NODE_CAPACITY 1024
#define INITIAL_PRIMITIVE_CAPACITY 256
#define INITIAL_SURFACE_CAPACITY 256
#define INITIAL_MATERIAL_CAPACITY 256
#define INITIAL_CELL_CAPACITY 256
#define INITIAL_TRANSFORM_CAPACITY 32

// ============================================================================
// DEFAULT CONFIGURATION
// ============================================================================

const alea_config_t ALEA_CONFIG_DEFAULT = {
    /* Tolerance */
    .abs_tol = 1e-14,
    .rel_tol = 1e-12,
    .zero_threshold = 1e-14,

    /* Behavior */
    .dedup = true,
    .log_level = 2, /* ALEA_LOG_WARN */

    /* Export */
    .surface_policy = ALEA_EMIT_MACROBODY,
    .export_materials = true,
    .export_transforms = true,
    .universe_depth = -1,
    .fill_depth = 0,
    .trcl_mode = 0, /* ALEA_TRCL_EXPORT_PRESERVE */
    .transform_mode = 0, /* ALEA_TR_EXPORT_ORIGINAL */
    .mcnp_max_col = 80,
    .mcnp_cont_indent = 5,

    /* Void generation */
    .void_max_depth = 8,
    .void_min_size = 0.1,
    .void_samples = 27,

    /* Void merge */
    .merge_cell_weight = 1.0,
    .merge_surface_weight = 0.1,
    .merge_max_surfaces = 24,
    .merge_min_cells = 1,

    /* Flatten */
    .flatten_max_depth = 0,
};

// ============================================================================
// CELL HASH TABLE (cell_id -> cell_index) — now uses alea_hashmap
// ============================================================================

// ============================================================================
// SYSTEM CREATION AND DESTRUCTION
// ============================================================================

alea_system_t* alea_system_create(void) {
    alea_system_t* sys = calloc(1, sizeof(alea_system_t));
    if (!sys) return NULL;

    /* All vectors are zero-initialized by calloc and grow automatically */
    sys->next_inline_transform_id = 1;
    sys->next_auto_surface_id = 1;  // Start at 1 for auto-assigned surface IDs
    sys->next_auto_cell_id = 1;     // Start at 1 for auto-assigned cell IDs
    sys->source = ALEA_SOURCE_EMPTY;
    /* cell_refs is zero-initialized by calloc (equivalent to ALEA_VEC_INIT) */

    // Create hash table for deduplication
    sys->primitive_index = primitive_hash_table_create();
    if (!sys->primitive_index) {
        alea_system_destroy(sys);
        return NULL;
    }

    // Create hash maps for cell and universe ID lookup
    sys->cell_index = cell_hashmap_create(256);
    sys->universe_index = universe_hashmap_create(64);

    // Set default configuration
    sys->config = ALEA_CONFIG_DEFAULT;

    // BVH acceleration (lazy-built)
    sys->surface_bvh = NULL;
    sys->bvh_dirty = true;

    // Track initial memory (all vectors start empty)
    sys->stats.current_memory = 0;
    sys->stats.peak_memory = 0;

    return sys;
}

void alea_system_destroy_internals(alea_system_t* sys) {
    if (!sys) return;

    // Invalidate thread-local spatial cache to prevent stale hits
    // if a new system is allocated at the same address (ABA problem)
    alea_spatial_reset_cache();

    // Free per-cell data
    for (size_t i = 0; i < alea_vec_count(&sys->cells); i++) {
        free(sys->cells.data[i].surface_indices);
        free(sys->cells.data[i].neighbors);
        free(sys->cells.data[i].lat_fill);
    }
    alea_vec_free(&sys->cells);

    alea_vec_free(&sys->nodes);
    alea_vec_free(&sys->primitives);
    alea_vec_free(&sys->surfaces);
    /* Free internal arrays of each material before freeing the vector */
    for (size_t i = 0; i < alea_vec_count(&sys->materials); i++) {
        for (size_t j = 0; j < sys->materials.data[i].nuclide_count; j++) {
            free(sys->materials.data[i].nuclides[j].library);
        }
        free(sys->materials.data[i].nuclides);
        for (size_t j = 0; j < sys->materials.data[i].element_count; j++) {
            free(sys->materials.data[i].elements[j].library);
        }
        free(sys->materials.data[i].elements);
        for (size_t j = 0; j < sys->materials.data[i].thermal_count; j++) {
            free(sys->materials.data[i].thermal_laws[j].identifier);
        }
        free(sys->materials.data[i].thermal_laws);
        free(sys->materials.data[i].name);
        free(sys->materials.data[i].comments);
    }
    alea_vec_free(&sys->materials);
    alea_vec_free(&sys->transforms);
    free(sys->surface_lookup);
    free(sys->prim_to_surface);
    free(sys->mcnp_id_to_surface);

    if (sys->primitive_index) {
        primitive_hash_table_destroy(sys->primitive_index);
    }

    cell_hashmap_destroy(&sys->cell_index);
    universe_hashmap_destroy(&sys->universe_index);

    // Free BVH
    if (sys->surface_bvh) {
        alea_bvh_free(sys->surface_bvh);
    }

    // Free universe index (each universe has internal arrays)
    for (size_t i = 0; i < alea_vec_count(&sys->universes); i++) {
        free(sys->universes.data[i].cell_indices);
    }
    alea_vec_free(&sys->universes);

    alea_free_cell_refs(sys);

    // Free spatial index
    if (sys->spatial_index) {
        alea_spatial_index_free(sys->spatial_index);
    }

    // Free mixtures
    for (size_t i = 0; i < alea_vec_count(&sys->mixtures); i++) {
        free(sys->mixtures.data[i].components);
        free(sys->mixtures.data[i].name);
        free(sys->mixtures.data[i].comments);
    }
    alea_vec_free(&sys->mixtures);
}

void alea_system_destroy(alea_system_t* sys) {
    if (!sys) return;
    alea_system_destroy_internals(sys);
    free(sys);
}

void alea_system_reset(alea_system_t* sys) {
    if (!sys) return;

    // Free per-cell dynamic arrays before reset
    for (size_t i = 0; i < alea_vec_count(&sys->cells); i++) {
        free(sys->cells.data[i].surface_indices);
        sys->cells.data[i].surface_indices = NULL;
        sys->cells.data[i].surface_index_count = 0;
        free(sys->cells.data[i].neighbors);
        sys->cells.data[i].neighbors = NULL;
        sys->cells.data[i].neighbor_count = 0;
        free(sys->cells.data[i].lat_fill);
        sys->cells.data[i].lat_fill = NULL;
        sys->cells.data[i].lat_fill_count = 0;
    }
    alea_vec_clear(&sys->cells);

    alea_vec_clear(&sys->nodes);
    alea_vec_clear(&sys->primitives);
    alea_vec_clear(&sys->surfaces);
    /* Free internal arrays of each material before clearing vector */
    for (size_t i = 0; i < alea_vec_count(&sys->materials); i++) {
        for (size_t j = 0; j < sys->materials.data[i].nuclide_count; j++) {
            free(sys->materials.data[i].nuclides[j].library);
        }
        free(sys->materials.data[i].nuclides);
        for (size_t j = 0; j < sys->materials.data[i].element_count; j++) {
            free(sys->materials.data[i].elements[j].library);
        }
        free(sys->materials.data[i].elements);
        for (size_t j = 0; j < sys->materials.data[i].thermal_count; j++) {
            free(sys->materials.data[i].thermal_laws[j].identifier);
        }
        free(sys->materials.data[i].thermal_laws);
        free(sys->materials.data[i].name);
        free(sys->materials.data[i].comments);
    }
    alea_vec_clear(&sys->materials);
    alea_vec_clear(&sys->transforms);
    sys->next_inline_transform_id = 1;
    
    // Free universe index internal arrays and clear vector
    for (size_t i = 0; i < alea_vec_count(&sys->universes); i++) {
        free(sys->universes.data[i].cell_indices);
        sys->universes.data[i].cell_indices = NULL;
    }
    alea_vec_clear(&sys->universes);
    universe_hashmap_clear(&sys->universe_index);
    sys->universe_index_built = false;

    // Free BVH (will be rebuilt on next use)
    if (sys->surface_bvh) {
        alea_bvh_free(sys->surface_bvh);
        sys->surface_bvh = NULL;
    }
    sys->bvh_dirty = true;

    // Free primitive-to-surface map
    free(sys->prim_to_surface);
    sys->prim_to_surface = NULL;
    sys->prim_to_surface_size = 0;

    // Free mcnp-id-to-surface map
    free(sys->mcnp_id_to_surface);
    sys->mcnp_id_to_surface = NULL;
    sys->mcnp_id_to_surface_size = 0;

    // Clear hash tables
    if (sys->primitive_index) {
        primitive_hash_table_destroy(sys->primitive_index);
        sys->primitive_index = primitive_hash_table_create();
    }
    cell_hashmap_clear(&sys->cell_index);
    
    memset(&sys->stats, 0, sizeof(alea_stats_t));
    alea_free_cell_refs(sys);
    /* cell_refs is zero-initialized by calloc (equivalent to ALEA_VEC_INIT) */                                                                     

}

/* grow_nodes removed - nodes now use vector API */
/* grow_primitives removed - primitives now use vector API */

// ============================================================================
// NODE OPERATIONS
// ============================================================================

alea_node_id_t alea_alloc_node(alea_system_t* sys) {
    if (!sys) {
        alea_set_error_detail(ALEA_ERR_NULL_ARG, "alea_alloc_node: system is NULL");
        return ALEA_NODE_ID_INVALID;
    }

    uint32_t id = (uint32_t)alea_vec_count(&sys->nodes);
    alea_node_t* node = alea_vec_push_uninit(&sys->nodes, alea_node_t);
    if (!node) {
        alea_set_error_detail(ALEA_ERR_OUT_OF_MEMORY, "alea_alloc_node: failed to allocate node %u", id);
        return ALEA_NODE_ID_INVALID;
    }

    memset(node, 0, sizeof(alea_node_t));
    return id;
}

alea_node_t* alea_get_node(alea_system_t* sys, uint32_t id) {
    if (!sys || id >= alea_vec_count(&sys->nodes)) return NULL;
    return &sys->nodes.data[id];
}

alea_node_id_t alea_clone_primitive(alea_system_t* sys, alea_node_id_t node_id, int8_t sense) {
    if (!sys || node_id >= alea_vec_count(&sys->nodes)) return ALEA_NODE_ID_INVALID;

    alea_node_t* src = &sys->nodes.data[node_id];
    if (ALEA_GET_OPERATION(src) != ALEA_OP_PRIMITIVE) return ALEA_NODE_ID_INVALID;

    alea_node_id_t new_id = alea_add_primitive_node(sys, src->primitive.primitive_id, sense,
                                                   src->primitive.inverted,
                                                   src->primitive.mcnp_surface_id);
    if (new_id == ALEA_NODE_ID_INVALID) return ALEA_NODE_ID_INVALID;

    /* Copy material_id from source (alea_add_primitive_node doesn't set this) */
    sys->nodes.data[new_id].material_id = src->material_id;

    return new_id;
}

// ============================================================================
// PRIMITIVE OPERATIONS WITH DEDUPLICATION
// ============================================================================

alea_primitive_id_t alea_get_or_create_primitive(alea_system_t* sys,
                                               alea_primitive_type_t type,
                                               alea_primitive_data_t* data,
                                               int8_t* inverted) {
    if (!sys) {
        alea_set_error_detail(ALEA_ERR_NULL_ARG, "alea_get_or_create_primitive: system is NULL");
        return UINT32_MAX;
    }

    // Canonicalize the primitive
    alea_canonicalize_primitive(type, data, inverted);

    // Compute hash
    uint64_t hash = alea_compute_primitive_hash(type, data, &sys->config);

    // Try to find existing primitive
    int8_t match_inverted = 0;
    uint32_t existing = primitive_hash_table_find(sys->primitive_index, sys,
                                                  type, data, hash, &match_inverted);

    if (existing != UINT32_MAX) {
        // Found duplicate!
        sys->primitives.data[existing].ref_count++;
        sys->stats.dedup_hits++;
        sys->stats.dedup_saved_bytes += sizeof(alea_primitive_entry_t);

        // If the match was via opposite normals, toggle the inverted flag
        if (match_inverted) {
            *inverted = !(*inverted);
            ALEA_LOG_INFO("Dedup hit: primitive %u (%s) reused with OPPOSITE normal (inverted=%d)",
                         existing, alea_primitive_type_name(type), *inverted);
        }

        if (alea_log_enabled(ALEA_LOG_LEVEL_DEBUG)) {
            int orig_surf_id = -1;
            for (size_t i = 0; i < alea_vec_count(&sys->surfaces); i++) {
                if (sys->surfaces.data[i].primitive_id == existing) {
                    orig_surf_id = sys->surfaces.data[i].mcnp_surface_id;
                    break;
                }
            }
            ALEA_LOG_DEBUG("Dedup hit: primitive %u (%s) reused, original surface %d (ref_count=%u)",
                          existing, alea_primitive_type_name(type),
                          orig_surf_id, sys->primitives.data[existing].ref_count);
        }

        return existing;
    }

    // Create new primitive using vector API
    uint32_t id = (uint32_t)alea_vec_count(&sys->primitives);
    alea_primitive_entry_t* prim = alea_vec_push_uninit(&sys->primitives, alea_primitive_entry_t);
    if (!prim) {
        alea_set_error_detail(ALEA_ERR_OUT_OF_MEMORY, "alea_get_or_create_primitive: failed to allocate primitive %u (type %d)", id, type);
        return UINT32_MAX;
    }

    prim->type = type;
    prim->data = *data;
    prim->ref_count = 1;

    // Add to hash table
    primitive_hash_table_insert(sys->primitive_index, id, hash);

    sys->stats.unique_primitives++;

    ALEA_LOG_INFO("New primitive %u (%s) created", id, alea_primitive_type_name(type));

    return id;
}

const alea_primitive_entry_t* alea_get_primitive(const alea_system_t* sys, uint32_t id) {
    if (!sys || id >= alea_vec_count(&sys->primitives)) return NULL;
    return &sys->primitives.data[id];
}

alea_node_id_t alea_add_primitive_node(alea_system_t* sys, uint32_t primitive_id,
                                      int8_t sense, int8_t inverted,
                                      int32_t mcnp_surface_id) {
    if (!sys || primitive_id >= alea_vec_count(&sys->primitives)) return ALEA_NODE_ID_INVALID;

    alea_node_id_t node_id = alea_alloc_node(sys);
    if (node_id == ALEA_NODE_ID_INVALID) return ALEA_NODE_ID_INVALID;

    alea_node_t* node = &sys->nodes.data[node_id];
    ALEA_SET_OPERATION(node, ALEA_OP_PRIMITIVE);
    node->primitive.primitive_id = primitive_id;
    node->primitive.prim_type = sys->primitives.data[primitive_id].type;
    node->primitive.sense = sense;
    node->primitive.inverted = inverted;
    node->primitive.mcnp_surface_id = mcnp_surface_id;

    /* Compute bounding box (sense-aware for proper halfspace bounds) */
    /* Note: inverted flips the effective sense */
    int8_t effective_sense = inverted ? -sense : sense;
    node->bbox = alea_halfspace_bbox(sys->primitives.data[primitive_id].type,
                                     &sys->primitives.data[primitive_id].data, effective_sense);

    sys->primitives.data[primitive_id].ref_count++;

    return node_id;
}

// ============================================================================
// STATISTICS
// ============================================================================

/**
 * @brief Compute current memory usage of all dynamic arrays in the system
 * @return Total allocated bytes
 */
size_t alea_system_memory_usage(const alea_system_t* sys) {
    if (!sys) return 0;

    size_t total = 0;

    /* Vector capacities (allocated memory, not just used) */
    total += sys->nodes.capacity * sizeof(alea_node_t);
    total += sys->primitives.capacity * sizeof(alea_primitive_entry_t);
    total += sys->surfaces.capacity * sizeof(alea_surface_entry_t);
    total += sys->materials.capacity * sizeof(alea_material_t);
    total += sys->cells.capacity * sizeof(alea_cell_entry_t);
    total += sys->transforms.capacity * sizeof(alea_transform_t);
    total += sys->universes.capacity * sizeof(alea_universe_t);
    total += sys->mixtures.capacity * sizeof(alea_mixture_t);
    total += sys->cell_refs.capacity * sizeof(alea_cell_ref_t);

    /* Additional allocations */
    if (sys->surface_lookup) {
        total += sys->surface_lookup_size * sizeof(alea_node_id_t);
    }
    if (sys->prim_to_surface) {
        total += sys->prim_to_surface_size * sizeof(uint32_t);
    }

    /* Per-material internal arrays */
    for (size_t i = 0; i < alea_vec_count(&sys->materials); i++) {
        const alea_material_t* m = &sys->materials.data[i];
        total += m->nuclide_count * sizeof(alea_nuclide_t);
        total += m->element_count * sizeof(alea_element_comp_t);
        total += m->thermal_count * sizeof(alea_thermal_law_t);
    }

    /* Per-cell internal arrays */
    for (size_t i = 0; i < alea_vec_count(&sys->cells); i++) {
        const alea_cell_entry_t* c = &sys->cells.data[i];
        if (c->surface_indices) {
            total += c->surface_index_count * sizeof(uint32_t);
        }
        if (c->neighbors) {
            total += c->neighbor_count * sizeof(alea_cell_neighbor_t);
        }
    }

    /* Per-universe internal arrays */
    for (size_t i = 0; i < alea_vec_count(&sys->universes); i++) {
        const alea_universe_t* u = &sys->universes.data[i];
        total += u->cell_capacity * sizeof(size_t);
    }

    return total;
}

void alea_system_print_stats(const alea_system_t* sys) {
    if (!sys) return;

    size_t memory = alea_system_memory_usage(sys);

    ALEA_LOG_DEBUG("=== CSG System Statistics ===\n");
    ALEA_LOG_DEBUG("Nodes: %zu\n", alea_vec_count(&sys->nodes));

    ALEA_LOG_DEBUG("Primitives: %zu unique, %zu total references\n",
           sys->stats.unique_primitives,
           sys->stats.dedup_hits + sys->stats.unique_primitives);

    if (sys->stats.dedup_hits > 0) {
        ALEA_LOG_DEBUG("Deduplication: %zu hits (%.1f%% reuse), saved %.1f MB\n",
               sys->stats.dedup_hits,
               100.0 * sys->stats.dedup_hits /
               (sys->stats.dedup_hits + sys->stats.unique_primitives),
               sys->stats.dedup_saved_bytes / 1048576.0);
    }

    ALEA_LOG_DEBUG("Memory: %.2f MB\n", memory / 1048576.0);

    if (sys->stats.surfaces_converted > 0) {
        ALEA_LOG_DEBUG("Surface conversion: %zu/%zu succeeded (%.1f%%)\n",
               sys->stats.surfaces_converted,
               sys->stats.surfaces_converted + sys->stats.failed_surfaces,
               100.0 * sys->stats.surfaces_converted /
               (sys->stats.surfaces_converted + sys->stats.failed_surfaces));
    }

    if (sys->stats.cells_converted > 0) {
        ALEA_LOG_DEBUG("Cells converted: %zu\n", sys->stats.cells_converted);
    }

    ALEA_LOG_DEBUG("Surfaces: %zu, Materials: %zu, Cells: %zu\n",
           alea_vec_count(&sys->surfaces), alea_vec_count(&sys->materials), alea_vec_count(&sys->cells));
}


// ============================================================================
// TRANSFORM OPERATIONS
// ============================================================================

/**
 * @brief Helper to populate a transform's data and cosines arrays
 */
static void populate_transform_data(alea_transform_t* tr, const double* data,
                                    int value_count, int degrees) {
    // Copy translation (first 3 values) - same for both arrays
    for (int i = 0; i < 3 && i < value_count; i++) {
        tr->data[i] = data[i];
        tr->cosines[i] = data[i];  // Translation is unchanged
    }

    // Handle rotation matrix (values 3-11)
    if (value_count == 12) {
        for (int i = 3; i < 12; i++) {
            tr->data[i] = data[i];  // Store original value
            if (degrees) {
                // Convert degrees to cosines
                tr->cosines[i] = cos(DEG_TO_RAD(data[i]));
            } else {
                // Already cosines
                tr->cosines[i] = data[i];
            }
        }
    }

    // Zero out unused slots
    for (int i = value_count; i < 12; i++) {
        tr->data[i] = 0.0;
        tr->cosines[i] = 0.0;
    }
}

int alea_add_transform(alea_system_t* sys, int transform_id,
                      const double* data, int value_count, int degrees) {
    if (!sys || !data) {
        alea_set_error_detail(ALEA_ERR_NULL_ARG, "alea_add_transform: NULL argument");
        return -1;
    }
    if (value_count != 3 && value_count != 12) {
        alea_set_error_detail(ALEA_ERR_INVALID_ARG, "alea_add_transform: value_count must be 3 or 12, got %d", value_count);
        return -1;
    }

    // Check if transform already exists
    for (size_t i = 0; i < alea_vec_count(&sys->transforms); i++) {
        if (sys->transforms.data[i].transform_id == transform_id) {
            // Update existing transform
            alea_transform_t* tr = &sys->transforms.data[i];
            populate_transform_data(tr, data, value_count, degrees);
            tr->value_count = value_count;
            tr->degrees = degrees;
            return 0;
        }
    }

    // Add new transform (vector grows automatically)
    alea_transform_t* tr = alea_vec_push_uninit(&sys->transforms, alea_transform_t);
    if (!tr) {
        alea_set_error_detail(ALEA_ERR_OUT_OF_MEMORY, "alea_add_transform: failed to allocate transform TR%d", transform_id);
        return -1;
    }

    tr->transform_id = transform_id;
    tr->value_count = value_count;
    tr->degrees = degrees;
    tr->from_inline = 0;

    populate_transform_data(tr, data, value_count, degrees);

    ALEA_LOG_INFO("Added transform TR%d (%d values, %s)",
           transform_id, value_count,
           degrees ? "degrees" : "cosines");

    return 0;
}

int alea_add_inline_transform(alea_system_t* sys, const double* data,
                             int value_count, int degrees) {
    if (!sys || !data) {
        alea_set_error_detail(ALEA_ERR_NULL_ARG, "alea_add_inline_transform: NULL argument");
        return -1;
    }
    if (value_count < 1 || value_count > 12) {
        alea_set_error_detail(ALEA_ERR_INVALID_ARG, "alea_add_inline_transform: value_count must be 1-12, got %d", value_count);
        return -1;
    }

    // Build candidate cosines for comparison
    alea_transform_t candidate;
    memset(&candidate, 0, sizeof(candidate));
    candidate.value_count = value_count;
    candidate.degrees = degrees;
    populate_transform_data(&candidate, data, value_count, degrees);

    // Deduplicate: scan existing inline transforms for a match on cosines
    static const double TRANSFORM_TOL = 1e-10;
    for (size_t i = 0; i < alea_vec_count(&sys->transforms); i++) {
        const alea_transform_t* existing = &sys->transforms.data[i];
        if (!existing->from_inline) continue;
        if (existing->value_count != value_count) continue;

        bool match = true;
        for (int j = 0; j < value_count; j++) {
            if (fabs(existing->cosines[j] - candidate.cosines[j]) > TRANSFORM_TOL) {
                match = false;
                break;
            }
        }
        if (match) {
            ALEA_LOG_INFO("Dedup inline transform -> existing ID %d",
                   existing->transform_id);
            return existing->transform_id;
        }
    }

    int assigned_id = sys->next_inline_transform_id++;

    // Add new transform (vector grows automatically)
    alea_transform_t* tr = alea_vec_push_uninit(&sys->transforms, alea_transform_t);
    if (!tr) {
        alea_set_error_detail(ALEA_ERR_OUT_OF_MEMORY, "alea_add_inline_transform: failed to allocate inline transform");
        return -1;
    }

    tr->transform_id = assigned_id;
    tr->value_count = candidate.value_count;
    tr->degrees = candidate.degrees;
    tr->from_inline = 1;
    memcpy(tr->data, candidate.data, sizeof(tr->data));
    memcpy(tr->cosines, candidate.cosines, sizeof(tr->cosines));

    ALEA_LOG_INFO("Added inline transform (assigned ID %d, %d values, %s)",
           assigned_id, value_count,
           degrees ? "degrees" : "cosines");

    return assigned_id;
}

void alea_finalize_transform_ids(alea_system_t* sys) {
    if (!sys) return;

    int max_id = 0;
    for (size_t i = 0; i < alea_vec_count(&sys->transforms); i++) {
        int id = sys->transforms.data[i].transform_id;
        if (id > max_id) max_id = id;
    }
    sys->next_inline_transform_id = max_id + 1;
}

// ============================================================================
// UNIVERSE INDEX
// ============================================================================

#define INITIAL_UNIVERSE_CELLS_CAPACITY 64

static alea_universe_t* find_or_create_universe(alea_system_t* sys, int universe_id) {
    // O(1) lookup via hashmap
    int* idx_ptr = universe_hashmap_get(&sys->universe_index, universe_id);
    if (idx_ptr) {
        return &sys->universes.data[*idx_ptr];
    }

    // Create new (vector grows automatically)
    int new_index = (int)alea_vec_count(&sys->universes);
    alea_universe_t* u = alea_vec_push_uninit(&sys->universes, alea_universe_t);
    if (!u) return NULL;

    memset(u, 0, sizeof(alea_universe_t));
    u->universe_id = universe_id;
    u->cell_indices = malloc(INITIAL_UNIVERSE_CELLS_CAPACITY * sizeof(size_t));
    if (!u->cell_indices) {
        alea_vec_pop_discard(&sys->universes);  // Rollback
        return NULL;
    }
    u->cell_capacity = INITIAL_UNIVERSE_CELLS_CAPACITY;
    u->bbox = (alea_bbox_t){1e30, -1e30, 1e30, -1e30, 1e30, -1e30};  // Empty bbox

    universe_hashmap_put(&sys->universe_index, universe_id, new_index);
    return u;
}

static int add_cell_to_universe(alea_universe_t* u, size_t cell_index) {
    if (u->cell_count >= u->cell_capacity) {
        size_t new_capacity = u->cell_capacity * 2;
        size_t* new_indices = realloc(u->cell_indices, new_capacity * sizeof(size_t));
        if (!new_indices) return -1;
        u->cell_indices = new_indices;
        u->cell_capacity = new_capacity;
    }
    u->cell_indices[u->cell_count++] = cell_index;
    return 0;
}

/**
 * @brief Find cell index by cell ID
 * @return Cell index (0-based) or -1 if not found
 */
int alea_find_cell_by_id(const alea_system_t* sys, int cell_id) {
    if (!sys) return -1;
    if (sys->cell_index.entries) {
        int* val = cell_hashmap_get(&sys->cell_index, cell_id);
        return val ? *val : -1;
    }
    /* Fallback linear scan if hash map not available */
    for (size_t i = 0; i < alea_vec_count(&sys->cells); i++) {
        if (sys->cells.data[i].mcnp_cell_id == cell_id) {
            return (int)i;
        }
    }
    return -1;
}

/**
 * @brief Internal: add cell with explicit ID (used by file loaders)
 *
 * Fast path - no duplicate checking. Call alea_validate_cell_ids() after
 * loading all cells to check for duplicates.
 */
int alea_add_cell_with_id(alea_system_t* sys, int cell_id, alea_node_id_t root_node,
                         int material_id, double density, int universe_id) {
    if (!sys) {
        alea_set_error_detail(ALEA_ERR_NULL_ARG, "alea_add_cell_with_id: system is NULL");
        return -1;
    }

    int idx = (int)alea_vec_count(&sys->cells);
    alea_cell_entry_t* cell = alea_vec_push_uninit(&sys->cells, alea_cell_entry_t);
    if (!cell) {
        alea_set_error_detail(ALEA_ERR_OUT_OF_MEMORY, "alea_add_cell_with_id: failed to allocate cell %d", cell_id);
        return -1;
    }
    memset(cell, 0, sizeof(*cell));

    cell->mcnp_cell_id = cell_id;
    cell->root_node_id = root_node;
    cell->original_root_node_id = ALEA_NODE_ID_INVALID;
    cell->material_id = material_id;
    cell->is_mass_density = (density < 0) ? 1 : 0;
    cell->density = fabs(density);
    cell->universe_id = universe_id;
    cell->imp_n = 1.0;
    cell->imp_p = 1.0;
    cell->imp_e = 1.0;

    sys->universe_index_built = false;

    /* Update next_auto_cell_id if this ID is >= it */
    if (cell_id >= sys->next_auto_cell_id) {
        sys->next_auto_cell_id = cell_id + 1;
    }

    /* Insert into cell hash map */
    cell_hashmap_put(&sys->cell_index, cell_id, idx);

    return idx;
}

/**
 * @brief Add cell with auto-assigned or validated ID
 *
 * If cell_id <= 0, auto-assigns a new ID.
 * If cell_id > 0 but already exists, logs warning and auto-assigns.
 */
int alea_add_cell(alea_system_t* sys, int cell_id, alea_node_id_t root_node,
                 int material_id, double density, int universe_id) {
    if (!sys) {
        alea_set_error_detail(ALEA_ERR_NULL_ARG, "alea_add_cell: system is NULL");
        return -1;
    }

    int final_cell_id;

    if (cell_id <= 0) {
        /* Auto-assign */
        final_cell_id = sys->next_auto_cell_id++;
    } else {
        /* Check for duplicate */
        if (alea_find_cell_by_id(sys, cell_id) >= 0) {
            ALEA_LOG_WARN("Cell ID %d already exists, auto-assigning new ID %d",
                         cell_id, sys->next_auto_cell_id);
            final_cell_id = sys->next_auto_cell_id++;
        } else {
            final_cell_id = cell_id;
            /* Update next_auto_cell_id if needed */
            if (cell_id >= sys->next_auto_cell_id) {
                sys->next_auto_cell_id = cell_id + 1;
            }
        }
    }

    /* Mark as programmatic if adding cells manually */
    if (sys->source == ALEA_SOURCE_EMPTY) {
        sys->source = ALEA_SOURCE_PROGRAMMATIC;
    }

    int idx = (int)alea_vec_count(&sys->cells);
    alea_cell_entry_t* cell = alea_vec_push_uninit(&sys->cells, alea_cell_entry_t);
    if (!cell) {
        alea_set_error_detail(ALEA_ERR_OUT_OF_MEMORY, "alea_add_cell: failed to allocate cell %d", final_cell_id);
        return -1;
    }
    memset(cell, 0, sizeof(*cell));

    cell->mcnp_cell_id = final_cell_id;
    cell->root_node_id = root_node;
    cell->original_root_node_id = ALEA_NODE_ID_INVALID;
    cell->material_id = material_id;
    cell->is_mass_density = (density < 0) ? 1 : 0;
    cell->density = fabs(density);
    cell->universe_id = universe_id;
    cell->imp_n = 1.0;
    cell->imp_p = 1.0;
    cell->imp_e = 1.0;

    sys->universe_index_built = false;

    /* Insert into cell hash map */
    cell_hashmap_put(&sys->cell_index, final_cell_id, idx);

    return idx;
}

int alea_max_cell_id(const alea_system_t* sys) {
    if (!sys || alea_vec_count(&sys->cells) == 0) return 0;

    int max_id = 0;
    for (size_t i = 0; i < alea_vec_count(&sys->cells); i++) {
        if (sys->cells.data[i].mcnp_cell_id > max_id) {
            max_id = sys->cells.data[i].mcnp_cell_id;
        }
    }
    return max_id;
}

/* alea_cell_get (individual pointer args) is in alea_public_api.c */
/* alea_cell_get_info (struct-based) is also in alea_public_api.c */

int alea_set_cell_fill(alea_system_t* sys, int cell_index, int fill_universe, int fill_transform) {
    if (!sys) return -1;
    if (cell_index < 0 || (size_t)cell_index >= alea_vec_count(&sys->cells)) return -1;

    sys->cells.data[cell_index].fill_universe = fill_universe;
    sys->cells.data[cell_index].fill_transform = fill_transform;
    return 0;
}

/* alea_cell_count, alea_universe_count defined in alea_public_api.c */

alea_material_id_t alea_material_at_point(const alea_system_t* sys, double x, double y, double z) {
    if (!sys) return ALEA_MATERIAL_NONE;

    int cell_idx = alea_identify_cell_at_point(sys, x, y, z);
    if (cell_idx < 0) return ALEA_MATERIAL_NONE;

    return sys->cells.data[cell_idx].material_id;
}

/* ============================================================================
 * MIXTURE OPERATIONS
 * ============================================================================ */

int alea_add_mixture(alea_system_t* sys, const alea_mixture_t* mixture) {
    if (!sys || !mixture) {
        alea_set_error_detail(ALEA_ERR_NULL_ARG, "alea_add_mixture: NULL argument");
        return -1;
    }

    /* Get index before push */
    size_t idx = alea_vec_count(&sys->mixtures);

    /* Allocate new slot (grows automatically if needed) */
    alea_mixture_t* dst = alea_vec_push_uninit(&sys->mixtures, alea_mixture_t);
    if (!dst) {
        alea_set_error_detail(ALEA_ERR_OUT_OF_MEMORY, "alea_add_mixture: failed to allocate mixture slot");
        return -1;
    }

    /* Deep copy the mixture */
    dst->mixture_id = mixture->mixture_id;
    dst->mcnp_material_id = mixture->mcnp_material_id;
    dst->is_weight_fraction = mixture->is_weight_fraction;
    dst->component_count = mixture->component_count;
    dst->component_capacity = mixture->component_count;

    /* Copy components array */
    if (mixture->component_count > 0) {
        dst->components = malloc(mixture->component_count * sizeof(alea_mixture_comp_t));
        if (!dst->components) {
            alea_set_error_detail(ALEA_ERR_OUT_OF_MEMORY, "alea_add_mixture: failed to allocate mixture components");
            /* Rollback: remove the element we just added */
            alea_vec_pop_discard(&sys->mixtures);
            return -1;
        }
        memcpy(dst->components, mixture->components,
               mixture->component_count * sizeof(alea_mixture_comp_t));
    } else {
        dst->components = NULL;
    }

    /* Copy strings if present */
    dst->name = mixture->name ? alea_strdup(mixture->name) : NULL;
    dst->comments = mixture->comments ? alea_strdup(mixture->comments) : NULL;

    return (int)idx;
}

/* Comparison function for qsort */
static int compare_cell_ids(const void* a, const void* b) {
    return *(const int*)a - *(const int*)b;
}

int alea_validate_cell_ids(alea_system_t* sys) {
    if (!sys) {
        alea_set_error_detail(ALEA_ERR_NULL_ARG, "alea_validate_cell_ids: system is NULL");
        return -1;
    }
    if (alea_vec_count(&sys->cells) <= 1) {
        return 0;  /* No duplicates possible */
    }

    /* Extract cell IDs into temporary array */
    int* ids = malloc(alea_vec_count(&sys->cells) * sizeof(int));
    if (!ids) {
        alea_set_error_detail(ALEA_ERR_OUT_OF_MEMORY, "alea_validate_cell_ids: failed to allocate memory for %zu IDs", alea_vec_count(&sys->cells));
        return -1;
    }

    for (size_t i = 0; i < alea_vec_count(&sys->cells); i++) {
        ids[i] = sys->cells.data[i].mcnp_cell_id;
    }

    /* Sort and check for adjacent duplicates */
    qsort(ids, alea_vec_count(&sys->cells), sizeof(int), compare_cell_ids);

    int duplicate_found = 0;
    for (size_t i = 1; i < alea_vec_count(&sys->cells); i++) {
        if (ids[i] == ids[i - 1]) {
            ALEA_LOG_ERROR("Duplicate cell ID: %d", ids[i]);
            duplicate_found = 1;
            /* Continue to report all duplicates */
        }
    }

    free(ids);

    if (duplicate_found) {
        alea_set_error_detail(ALEA_ERR_INVALID_ID, "alea_validate_cell_ids: duplicate cell IDs found");
        return -1;
    }

    return 0;
}

int alea_build_universe_index(alea_system_t* sys) {
    if (!sys) {
        alea_set_error_detail(ALEA_ERR_NULL_ARG, "alea_build_universe_index: system is NULL");
        return -1;
    }

    // Clear existing index (free internal arrays first)
    for (size_t i = 0; i < alea_vec_count(&sys->universes); i++) {
        free(sys->universes.data[i].cell_indices);
    }
    alea_vec_clear(&sys->universes);
    universe_hashmap_clear(&sys->universe_index);

    // Build index
    for (size_t i = 0; i < alea_vec_count(&sys->cells); i++) {
        ALEA_CHECK_INTERRUPTED(-1);
        const alea_cell_entry_t* cell = &sys->cells.data[i];

        alea_universe_t* u = find_or_create_universe(sys, cell->universe_id);
        if (!u) {
            alea_set_error_detail(ALEA_ERR_OUT_OF_MEMORY, "alea_build_universe_index: failed to create universe %d", cell->universe_id);
            return -1;
        }

        if (add_cell_to_universe(u, i) < 0) {
            alea_set_error_detail(ALEA_ERR_OUT_OF_MEMORY, "alea_build_universe_index: failed to add cell to universe %d", cell->universe_id);
            return -1;
        }

        // Expand universe bbox with cell bbox
        if (cell->root_node_id < alea_vec_count(&sys->nodes)) {
            const alea_node_t* node = &sys->nodes.data[cell->root_node_id];
            if (node->bbox.min_x < u->bbox.min_x) u->bbox.min_x = node->bbox.min_x;
            if (node->bbox.max_x > u->bbox.max_x) u->bbox.max_x = node->bbox.max_x;
            if (node->bbox.min_y < u->bbox.min_y) u->bbox.min_y = node->bbox.min_y;
            if (node->bbox.max_y > u->bbox.max_y) u->bbox.max_y = node->bbox.max_y;
            if (node->bbox.min_z < u->bbox.min_z) u->bbox.min_z = node->bbox.min_z;
            if (node->bbox.max_z > u->bbox.max_z) u->bbox.max_z = node->bbox.max_z;
        }
    }

    /* Check if any cell has a lattice */
    sys->has_lattice = false;
    for (size_t i = 0; i < alea_vec_count(&sys->cells); i++) {
        if (sys->cells.data[i].lat_type != 0 && sys->cells.data[i].lat_fill) {
            sys->has_lattice = true;
            break;
        }
    }

    sys->universe_index_built = true;

    ALEA_LOG_INFO("Built universe index: %zu universes", alea_vec_count(&sys->universes));
    for (size_t i = 0; i < alea_vec_count(&sys->universes); i++) {
        ALEA_LOG_INFO("  Universe %d: %zu cells",
               sys->universes.data[i].universe_id, sys->universes.data[i].cell_count);
    }

    return 0;
}

const alea_universe_t* alea_get_universe(const alea_system_t* sys, int universe_id) {
    if (!sys || alea_vec_empty(&sys->universes)) return NULL;

    /* O(1) lookup via hashmap */
    int* idx_ptr = universe_hashmap_get(&sys->universe_index, universe_id);
    if (idx_ptr) {
        return &sys->universes.data[*idx_ptr];
    }
    /* Fallback: linear scan (hashmap may not be populated yet) */
    for (size_t i = 0; i < alea_vec_count(&sys->universes); i++) {
        if (sys->universes.data[i].universe_id == universe_id) {
            return &sys->universes.data[i];
        }
    }
    return NULL;
}

int alea_get_universe_cells(const alea_system_t* sys, int universe_id,
                           const alea_cell_entry_t** out_cells, size_t max_cells) {
    if (!sys || !out_cells) return -1;
    
    const alea_universe_t* u = alea_get_universe(sys, universe_id);
    if (!u) return 0;
    
    size_t count = (u->cell_count < max_cells) ? u->cell_count : max_cells;
    for (size_t i = 0; i < count; i++) {
        out_cells[i] = &sys->cells.data[u->cell_indices[i]];
    }
    return (int)count;
}

int alea_identify_cell_at_point(const alea_system_t* sys, double x, double y, double z) {
    if (!sys) return -1;

    // Test all cells in base universe (universe 0)
    // If universe index not built, test all cells with universe_id == 0
    if (!sys->universe_index_built) {
        if (alea_build_universe_index((alea_system_t*)sys) < 0) {
            return -1;
        }
    }
    
    if (sys->universe_index_built) {
        const alea_universe_t* base = alea_get_universe(sys, 0);
        
        if (!base) return -1;

        
        for (size_t i = 0; i < base->cell_count; i++) {
            size_t cell_idx = base->cell_indices[i];
            const alea_cell_entry_t* cell = &sys->cells.data[cell_idx];
            if (alea_contains_point(sys, cell->root_node_id, x, y, z)) {
                return (int)cell_idx;
            }
        }
    } else {
        // Linear search through all cells
        for (size_t i = 0; i < alea_vec_count(&sys->cells); i++) {
            const alea_cell_entry_t* cell = &sys->cells.data[i];
            if (cell->universe_id != 0) continue;
            
            if (alea_contains_point(sys, cell->root_node_id, x, y, z)) {
                return (int)i;
            }
        }
    }
    
    return -1;  // In void
}

int alea_find_overlaps(const alea_system_t* sys, int* out_pairs, size_t max_pairs) {
    if (!sys || !out_pairs || max_pairs == 0) return 0;
    
    size_t found = 0;
    
    // Only check base universe cells
    for (size_t i = 0; i < alea_vec_count(&sys->cells) && found < max_pairs; i++) {
        ALEA_CHECK_INTERRUPTED(-1);
        const alea_cell_entry_t* cell_i = &sys->cells.data[i];
        if (cell_i->universe_id != 0) continue;
        if (cell_i->root_node_id >= alea_vec_count(&sys->nodes)) continue;

        const alea_bbox_t* bbox_i = &sys->nodes.data[cell_i->root_node_id].bbox;

        for (size_t j = i + 1; j < alea_vec_count(&sys->cells) && found < max_pairs; j++) {
            const alea_cell_entry_t* cell_j = &sys->cells.data[j];
            if (cell_j->universe_id != 0) continue;
            if (cell_j->root_node_id >= alea_vec_count(&sys->nodes)) continue;
            
            const alea_bbox_t* bbox_j = &sys->nodes.data[cell_j->root_node_id].bbox;
            
            // Check bbox overlap first (fast rejection)
            if (bbox_i->max_x < bbox_j->min_x || bbox_j->max_x < bbox_i->min_x ||
                bbox_i->max_y < bbox_j->min_y || bbox_j->max_y < bbox_i->min_y ||
                bbox_i->max_z < bbox_j->min_z || bbox_j->max_z < bbox_i->min_z) {
                continue;  // No bbox overlap
            }
            
            // Sample points in intersection of bboxes
            double min_x = (bbox_i->min_x > bbox_j->min_x) ? bbox_i->min_x : bbox_j->min_x;
            double max_x = (bbox_i->max_x < bbox_j->max_x) ? bbox_i->max_x : bbox_j->max_x;
            double min_y = (bbox_i->min_y > bbox_j->min_y) ? bbox_i->min_y : bbox_j->min_y;
            double max_y = (bbox_i->max_y < bbox_j->max_y) ? bbox_i->max_y : bbox_j->max_y;
            double min_z = (bbox_i->min_z > bbox_j->min_z) ? bbox_i->min_z : bbox_j->min_z;
            double max_z = (bbox_i->max_z < bbox_j->max_z) ? bbox_i->max_z : bbox_j->max_z;
            
            // Sample 8 points (corners) + center
            double test_points[9][3] = {
                {min_x, min_y, min_z}, {max_x, min_y, min_z},
                {min_x, max_y, min_z}, {max_x, max_y, min_z},
                {min_x, min_y, max_z}, {max_x, min_y, max_z},
                {min_x, max_y, max_z}, {max_x, max_y, max_z},
                {(min_x+max_x)/2, (min_y+max_y)/2, (min_z+max_z)/2}
            };
            
            bool overlap_found = false;
            for (int p = 0; p < 9 && !overlap_found; p++) {
                if (alea_contains_point(sys, cell_i->root_node_id, 
                                       test_points[p][0], test_points[p][1], test_points[p][2]) &&
                    alea_contains_point(sys, cell_j->root_node_id,
                                       test_points[p][0], test_points[p][1], test_points[p][2])) {
                    overlap_found = true;
                }
            }
            
            if (overlap_found) {
                out_pairs[found * 2] = (int)i;
                out_pairs[found * 2 + 1] = (int)j;
                found++;
            }
        }
    }
    
    return (int)found;
}

// ============================================================================
// TRANSFORM OPERATIONS
// ============================================================================

const alea_transform_t* alea_get_transform(const alea_system_t* sys, int transform_id) {
    if (!sys || transform_id == 0) return NULL;
    for (size_t i = 0; i < alea_vec_count(&sys->transforms); i++) {
        if (sys->transforms.data[i].transform_id == transform_id) {
            return &sys->transforms.data[i];
        }
    }
    return NULL;
}

bool alea_transform_is_translation_only(const alea_transform_t* tr) {
    if (!tr) return true;  // No transform = identity = translation only
    return tr->value_count <= 3;
}

void alea_transform_point(const alea_transform_t* tr, double* x, double* y, double* z) {
    if (!tr || !x || !y || !z) return;

    /* Use tr->cosines which contains pre-computed direction cosines */
    double ox = tr->cosines[0];
    double oy = tr->cosines[1];
    double oz = tr->cosines[2];

    if (tr->value_count <= 3) {
        // Translation only
        *x += ox;
        *y += oy;
        *z += oz;
    } else {
        // Full rotation + translation
        // MCNP convention: R has COLUMNS [b1,b2,b3], [b4,b5,b6], [b7,b8,b9]
        // p_main = o + R * p_aux where R = [b1 b4 b7; b2 b5 b8; b3 b6 b9]
        double b1 = tr->cosines[3], b2 = tr->cosines[4], b3 = tr->cosines[5];
        double b4 = tr->cosines[6], b5 = tr->cosines[7], b6 = tr->cosines[8];
        double b7 = tr->cosines[9], b8 = tr->cosines[10], b9 = tr->cosines[11];

        double px = *x, py = *y, pz = *z;

        // p_main = o + R * p_aux (columns are [b1,b2,b3], [b4,b5,b6], [b7,b8,b9])
        *x = ox + b1*px + b4*py + b7*pz;
        *y = oy + b2*px + b5*py + b8*pz;
        *z = oz + b3*px + b6*py + b9*pz;
    }
}

void alea_transform_vector(const alea_transform_t* tr, double* vx, double* vy, double* vz) {
    if (!tr || !vx || !vy || !vz) return;

    if (tr->value_count <= 3) {
        // Translation only - vectors unchanged
        return;
    }

    /* Use tr->cosines which contains pre-computed direction cosines */
    // Rotation only (no translation for vectors)
    // MCNP convention: R has COLUMNS [b1,b2,b3], [b4,b5,b6], [b7,b8,b9]
    double b1 = tr->cosines[3], b2 = tr->cosines[4], b3 = tr->cosines[5];
    double b4 = tr->cosines[6], b5 = tr->cosines[7], b6 = tr->cosines[8];
    double b7 = tr->cosines[9], b8 = tr->cosines[10], b9 = tr->cosines[11];

    double px = *vx, py = *vy, pz = *vz;

    // v_main = R * v_aux (columns are [b1,b2,b3], [b4,b5,b6], [b7,b8,b9])
    *vx = b1*px + b4*py + b7*pz;
    *vy = b2*px + b5*py + b8*pz;
    *vz = b3*px + b6*py + b9*pz;
}


bool alea_apply_transform_to_primitive(const alea_transform_t* tr,
                                      alea_primitive_type_t in_type,
                                      const alea_primitive_data_t* in_data,
                                      alea_primitive_type_t* out_type,
                                      alea_primitive_data_t* out_data) {
    if (!tr || !in_data || !out_type || !out_data) {
        /* No transform - just copy */
        if (out_type) *out_type = in_type;
        if (out_data && in_data) *out_data = *in_data;
        return true;
    }

    /* Convert alea_transform_t to 3x4 row-major matrix
     * alea_transform_t.cosines layout: [Ox,Oy,Oz, b1,b2,b3, b4,b5,b6, b7,b8,b9]
     * where R = [b1 b4 b7; b2 b5 b8; b3 b6 b9] (column vectors)
     * Our format: [R00,R01,R02,Tx, R10,R11,R12,Ty, R20,R21,R22,Tz] (row-major)
     *
     * Note: Use tr->cosines which contains pre-computed direction cosines,
     * regardless of whether the original input was in degrees.
     */
    double mat[12];

    if (tr->value_count <= 3) {
        /* Translation only - identity rotation */
        mat[0] = 1.0;  mat[1] = 0.0;  mat[2] = 0.0;  mat[3] = tr->cosines[0];
        mat[4] = 0.0;  mat[5] = 1.0;  mat[6] = 0.0;  mat[7] = tr->cosines[1];
        mat[8] = 0.0;  mat[9] = 0.0;  mat[10] = 1.0; mat[11] = tr->cosines[2];
    } else {
        /* Full rotation + translation */
        mat[0] = tr->cosines[3];   /* b1 = R00 */
        mat[1] = tr->cosines[6];   /* b4 = R01 */
        mat[2] = tr->cosines[9];   /* b7 = R02 */
        mat[3] = tr->cosines[0];   /* Ox = Tx */
        mat[4] = tr->cosines[4];   /* b2 = R10 */
        mat[5] = tr->cosines[7];   /* b5 = R11 */
        mat[6] = tr->cosines[10];  /* b8 = R12 */
        mat[7] = tr->cosines[1];   /* Oy = Ty */
        mat[8] = tr->cosines[5];   /* b3 = R20 */
        mat[9] = tr->cosines[8];   /* b6 = R21 */
        mat[10] = tr->cosines[11]; /* b9 = R22 */
        mat[11] = tr->cosines[2];  /* Oz = Tz */
    }

    return alea_primitive_transform(in_type, in_data, mat, out_type, out_data);
}

bool alea_apply_inverse_transform_to_primitive(const alea_transform_t* tr,
                                              alea_primitive_type_t in_type,
                                              const alea_primitive_data_t* in_data,
                                              alea_primitive_type_t* out_type,
                                              alea_primitive_data_t* out_data) {
    if (!tr || !in_data || !out_type || !out_data) {
        /* No transform - just copy */
        if (out_type) *out_type = in_type;
        if (out_data && in_data) *out_data = *in_data;
        return true;
    }

    /* Build inverse matrix from alea_transform_t
     * Forward: p' = R * p + T
     * Inverse: p = R^T * (p' - T) = R^T * p' - R^T * T
     *
     * For the inverse matrix in row-major [R00,R01,R02,Tx, R10,R11,R12,Ty, R20,R21,R22,Tz]:
     * R_inv = R^T, T_inv = -R^T * T
     *
     * Note: Use tr->cosines which contains pre-computed direction cosines.
     */
    double mat[12];

    if (tr->value_count <= 3) {
        /* Translation only - just negate translation */
        mat[0] = 1.0;  mat[1] = 0.0;  mat[2] = 0.0;  mat[3] = -tr->cosines[0];
        mat[4] = 0.0;  mat[5] = 1.0;  mat[6] = 0.0;  mat[7] = -tr->cosines[1];
        mat[8] = 0.0;  mat[9] = 0.0;  mat[10] = 1.0; mat[11] = -tr->cosines[2];
    } else {
        /* Full rotation + translation
         * Original R (column vectors): [b1 b4 b7; b2 b5 b8; b3 b6 b9]
         * R^T (inverse rotation) = [b1 b2 b3; b4 b5 b6; b7 b8 b9]
         * T_inv = -R^T * T
         */
        double b1 = tr->cosines[3], b2 = tr->cosines[4], b3 = tr->cosines[5];
        double b4 = tr->cosines[6], b5 = tr->cosines[7], b6 = tr->cosines[8];
        double b7 = tr->cosines[9], b8 = tr->cosines[10], b9 = tr->cosines[11];
        double ox = tr->cosines[0], oy = tr->cosines[1], oz = tr->cosines[2];

        /* R^T (transpose of R) */
        mat[0] = b1;  mat[1] = b2;  mat[2] = b3;   /* First row of R^T */
        mat[4] = b4;  mat[5] = b5;  mat[6] = b6;   /* Second row of R^T */
        mat[8] = b7;  mat[9] = b8;  mat[10] = b9;  /* Third row of R^T */

        /* T_inv = -R^T * T */
        mat[3]  = -(b1*ox + b2*oy + b3*oz);
        mat[7]  = -(b4*ox + b5*oy + b6*oz);
        mat[11] = -(b7*ox + b8*oy + b9*oz);
    }

    return alea_primitive_transform(in_type, in_data, mat, out_type, out_data);
}

// ============================================================================
// CELL SURFACE INDEX
// ============================================================================

/**
 * Build/rebuild the primitive-to-surface mapping stored in sys
 */
static int ensure_prim_to_surface_map(alea_system_t* sys) {
    size_t prim_count = alea_vec_count(&sys->primitives);

    /* Check if map needs rebuilding */
    if (sys->prim_to_surface && sys->prim_to_surface_size >= prim_count) {
        return 0;  /* Already built and big enough */
    }

    /* Free old map */
    free(sys->prim_to_surface);
    sys->prim_to_surface = NULL;
    sys->prim_to_surface_size = 0;

    if (prim_count == 0) {
        return 0;
    }

    /* Allocate new map */
    sys->prim_to_surface = malloc(prim_count * sizeof(uint32_t));
    if (!sys->prim_to_surface) {
        return -1;
    }
    sys->prim_to_surface_size = prim_count;

    /* Initialize to "no surface" */
    for (size_t i = 0; i < sys->prim_to_surface_size; i++) {
        sys->prim_to_surface[i] = UINT32_MAX;
    }

    /* Build mapping from surfaces */
    for (size_t i = 0; i < alea_vec_count(&sys->surfaces); i++) {
        uint32_t prim_id = sys->surfaces.data[i].primitive_id;
        if (prim_id < sys->prim_to_surface_size) {
            sys->prim_to_surface[prim_id] = (uint32_t)i;
        }
    }

    return 0;
}

/**
 * Build/rebuild the mcnp-surface-id-to-surface mapping stored in sys
 * This provides O(1) lookup by MCNP surface ID instead of O(n) linear search.
 */
static int ensure_mcnp_id_to_surface_map(alea_system_t* sys) {
    size_t surf_count = alea_vec_count(&sys->surfaces);

    if (surf_count == 0) {
        return 0;
    }

    /* Find max MCNP surface ID to size the array */
    int max_mcnp_id = 0;
    for (size_t i = 0; i < surf_count; i++) {
        int id = sys->surfaces.data[i].mcnp_surface_id;
        if (id > max_mcnp_id) {
            max_mcnp_id = id;
        }
    }

    if (max_mcnp_id <= 0) {
        return 0;  /* No valid MCNP surface IDs */
    }

    size_t needed_size = (size_t)(max_mcnp_id + 1);

    /* Check if map needs rebuilding */
    if (sys->mcnp_id_to_surface && sys->mcnp_id_to_surface_size >= needed_size) {
        return 0;  /* Already built and big enough */
    }

    /* Free old map */
    free(sys->mcnp_id_to_surface);
    sys->mcnp_id_to_surface = NULL;
    sys->mcnp_id_to_surface_size = 0;

    /* Allocate new map */
    sys->mcnp_id_to_surface = malloc(needed_size * sizeof(uint32_t));
    if (!sys->mcnp_id_to_surface) {
        return -1;
    }
    sys->mcnp_id_to_surface_size = needed_size;

    /* Initialize to "no surface" */
    for (size_t i = 0; i < needed_size; i++) {
        sys->mcnp_id_to_surface[i] = UINT32_MAX;
    }

    /* Build mapping from surfaces */
    for (size_t i = 0; i < surf_count; i++) {
        int mcnp_id = sys->surfaces.data[i].mcnp_surface_id;
        if (mcnp_id > 0 && (size_t)mcnp_id < needed_size) {
            sys->mcnp_id_to_surface[mcnp_id] = (uint32_t)i;
        }
    }

    return 0;
}

/**
 * Helper: find surface index by MCNP surface ID
 * Returns UINT32_MAX if not found
 * Uses O(1) lookup table if available, falls back to O(n) search otherwise.
 */
static uint32_t find_surface_by_mcnp_id(alea_system_t* sys, int mcnp_surface_id) {
    if (mcnp_surface_id <= 0) {
        return UINT32_MAX;
    }

    /* Use lookup table if available */
    if (sys->mcnp_id_to_surface && (size_t)mcnp_surface_id < sys->mcnp_id_to_surface_size) {
        return sys->mcnp_id_to_surface[mcnp_surface_id];
    }

    /* Fallback: linear search (should not happen if map is built) */
    for (size_t i = 0; i < alea_vec_count(&sys->surfaces); i++) {
        if (sys->surfaces.data[i].mcnp_surface_id == mcnp_surface_id) {
            return (uint32_t)i;
        }
    }
    return UINT32_MAX;
}

/**
 * Helper: recursively collect surface indices from a CSG tree
 *
 * Uses node->primitive.mcnp_surface_id to find the correct surface,
 * which handles cases where multiple surfaces share the same primitive
 * due to primitive deduplication. Falls back to prim_to_surface if
 * mcnp_surface_id is not set.
 */
static void collect_surfaces_recursive(alea_system_t* sys,
                                       alea_node_id_t node_id,
                                       uint32_t** surface_indices,
                                       size_t* count,
                                       size_t* capacity) {
    if (node_id == ALEA_NODE_ID_INVALID || node_id >= alea_vec_count(&sys->nodes))
        return;

    const alea_node_t* node = &sys->nodes.data[node_id];
    alea_operation_t op = ALEA_GET_OPERATION(node);

    if (op == ALEA_OP_PRIMITIVE) {
        uint32_t surf_idx = UINT32_MAX;

        /* First try: use node's mcnp_surface_id for correct surface lookup
         * This handles primitive deduplication correctly */
        int mcnp_id = node->primitive.mcnp_surface_id;
        if (mcnp_id > 0) {
            surf_idx = find_surface_by_mcnp_id(sys, mcnp_id);
        }

        /* Fallback: use prim_to_surface if mcnp_id not set or not found */
        if (surf_idx == UINT32_MAX) {
            alea_primitive_id_t prim_id = node->primitive.primitive_id;
            if (prim_id < sys->prim_to_surface_size) {
                surf_idx = sys->prim_to_surface[prim_id];
            }
        }

        if (surf_idx != UINT32_MAX) {
            /* Check if already in array */
            bool found = false;
            for (size_t j = 0; j < *count; j++) {
                if ((*surface_indices)[j] == surf_idx) {
                    found = true;
                    break;
                }
            }
            if (!found) {
                /* Add to array */
                if (*count >= *capacity) {
                    size_t new_cap = *capacity ? *capacity * 2 : 16;
                    uint32_t* new_arr = realloc(*surface_indices, new_cap * sizeof(uint32_t));
                    if (!new_arr) return;
                    *surface_indices = new_arr;
                    *capacity = new_cap;
                }
                (*surface_indices)[(*count)++] = surf_idx;
            }
        }
    } else if (op == ALEA_OP_COMPLEMENT) {
        /* Complement has one child in left */
        collect_surfaces_recursive(sys, node->operation.left, surface_indices, count, capacity);
    } else {
        /* Binary operation: union, intersection, difference */
        collect_surfaces_recursive(sys, node->operation.left, surface_indices, count, capacity);
        collect_surfaces_recursive(sys, node->operation.right, surface_indices, count, capacity);
    }
}

/**
 * @brief Build surface index for all cells
 *
 * Traverses each cell's CSG tree and collects all surface indices
 * that the cell references. This enables faster raycasting by only
 * testing surfaces relevant to each cell.
 *
 * @param sys The CSG system
 * @return 0 on success, -1 on failure
 */
int alea_build_cell_surface_index(alea_system_t* sys) {
    if (!sys) return -1;

    /* Build/update primitive->surface map */
    if (ensure_prim_to_surface_map(sys) != 0) {
        return -1;
    }

    /* Build/update MCNP surface ID -> surface index map for O(1) lookup */
    if (ensure_mcnp_id_to_surface_map(sys) != 0) {
        return -1;
    }

    for (size_t i = 0; i < alea_vec_count(&sys->cells); i++) {
        alea_cell_entry_t* cell = &sys->cells.data[i];

        /* Free old index if present */
        free(cell->surface_indices);
        cell->surface_indices = NULL;
        cell->surface_index_count = 0;

        /* Skip void cells (no geometry) */
        if (cell->root_node_id == ALEA_NODE_ID_INVALID)
            continue;

        /* Collect surface indices */
        uint32_t* indices = NULL;
        size_t count = 0;
        size_t capacity = 0;

        collect_surfaces_recursive(sys, cell->root_node_id, &indices, &count, &capacity);

        /* Store result (shrink to exact size) */
        if (count > 0) {
            cell->surface_indices = realloc(indices, count * sizeof(uint32_t));
            if (!cell->surface_indices) {
                cell->surface_indices = indices;  /* Keep original if shrink fails */
            }
            cell->surface_index_count = count;
        } else {
            free(indices);
        }
    }

    return 0;
}

/* ============================================================================
 * CELL ADJACENCY
 * ============================================================================ */

/**
 * @brief Entry in the surface-to-cells map
 */
typedef struct {
    uint32_t cell_index;
    int8_t sense;  /* +1 or -1 */
} surface_cell_ref_t;

typedef struct {
    surface_cell_ref_t* refs;
    size_t count;
    size_t capacity;
} surface_cell_list_t;

/**
 * @brief Helper: collect primitive senses from a cell's CSG tree
 */
typedef struct {
    uint32_t surface_idx;
    int8_t sense;
} cell_surface_sense_t;

static void collect_cell_surface_senses(const alea_system_t* sys,
                                        alea_node_id_t node_id,
                                        int current_sense,
                                        cell_surface_sense_t** out_senses,
                                        size_t* count,
                                        size_t* capacity) {
    if (node_id == ALEA_NODE_ID_INVALID || node_id >= alea_vec_count(&sys->nodes))
        return;

    const alea_node_t* node = &sys->nodes.data[node_id];
    alea_operation_t op = ALEA_GET_OPERATION(node);

    if (op == ALEA_OP_PRIMITIVE) {
        alea_primitive_id_t prim_id = node->primitive.primitive_id;
        int8_t node_sense = node->primitive.sense;
        int8_t inverted = node->primitive.inverted;
        /* Combine: node sense * dedup inversion * complement context */
        int8_t final_sense = (int8_t)(current_sense * node_sense * (inverted ? -1 : 1));

        if (prim_id < sys->prim_to_surface_size) {
            uint32_t surf_idx = sys->prim_to_surface[prim_id];
            if (surf_idx != UINT32_MAX) {
                /* Check if already collected */
                for (size_t i = 0; i < *count; i++) {
                    if ((*out_senses)[i].surface_idx == surf_idx) {
                        return;  /* Already have this surface */
                    }
                }
                /* Add to collection */
                if (*count >= *capacity) {
                    size_t new_cap = *capacity ? *capacity * 2 : 16;
                    cell_surface_sense_t* new_arr = realloc(*out_senses,
                        new_cap * sizeof(cell_surface_sense_t));
                    if (!new_arr) return;
                    *out_senses = new_arr;
                    *capacity = new_cap;
                }
                (*out_senses)[*count].surface_idx = surf_idx;
                (*out_senses)[*count].sense = final_sense;
                (*count)++;
            }
        }
    } else if (op == ALEA_OP_COMPLEMENT) {
        collect_cell_surface_senses(sys, node->operation.left, -current_sense,
                                   out_senses, count, capacity);
    } else {
        collect_cell_surface_senses(sys, node->operation.left, current_sense,
                                   out_senses, count, capacity);
        collect_cell_surface_senses(sys, node->operation.right, current_sense,
                                   out_senses, count, capacity);
    }
}

int alea_build_cell_adjacency(alea_system_t* sys) {
    if (!sys) return -1;

    /* Early exit if already built */
    if (sys->cell_adjacency_built) return 0;

    /* Ensure surface index is built */
    if (alea_build_cell_surface_index(sys) != 0) {
        return -1;
    }

    /* Free existing adjacency data */
    for (size_t i = 0; i < alea_vec_count(&sys->cells); i++) {
        free(sys->cells.data[i].neighbors);
        sys->cells.data[i].neighbors = NULL;
        sys->cells.data[i].neighbor_count = 0;
    }
    sys->cell_adjacency_built = false;

    /* Create surface->cells map */
    surface_cell_list_t* surf_map = calloc(alea_vec_count(&sys->surfaces), sizeof(surface_cell_list_t));
    if (!surf_map) return -1;

    /* For each cell, collect its surfaces with senses and add to map */
    for (size_t ci = 0; ci < alea_vec_count(&sys->cells); ci++) {
        alea_cell_entry_t* cell = &sys->cells.data[ci];
        if (cell->root_node_id == ALEA_NODE_ID_INVALID) continue;

        /* Collect surfaces with senses */
        cell_surface_sense_t* senses = NULL;
        size_t sense_count = 0;
        size_t sense_cap = 0;

        collect_cell_surface_senses(sys, cell->root_node_id, 1,
                                   &senses, &sense_count, &sense_cap);

        /* Add cell references to surface map */
        for (size_t si = 0; si < sense_count; si++) {
            uint32_t surf_idx = senses[si].surface_idx;
            if (surf_idx >= alea_vec_count(&sys->surfaces)) continue;

            surface_cell_list_t* list = &surf_map[surf_idx];

            /* Grow if needed */
            if (list->count >= list->capacity) {
                size_t new_cap = list->capacity ? list->capacity * 2 : 4;
                surface_cell_ref_t* new_arr = realloc(list->refs,
                    new_cap * sizeof(surface_cell_ref_t));
                if (!new_arr) {
                    free(senses);
                    goto cleanup;
                }
                list->refs = new_arr;
                list->capacity = new_cap;
            }

            list->refs[list->count].cell_index = (uint32_t)ci;
            list->refs[list->count].sense = senses[si].sense;
            list->count++;
        }

        free(senses);
    }

    /* Build neighbor lists from surface map */
    /* For each surface, if it has cells on both sides (opposite senses),
     * those cells are neighbors */
    for (size_t si = 0; si < alea_vec_count(&sys->surfaces); si++) {
        surface_cell_list_t* list = &surf_map[si];

        /* Need at least 2 cells to have neighbors */
        if (list->count < 2) continue;

        /* Find pairs with opposite senses */
        for (size_t a = 0; a < list->count; a++) {
            for (size_t b = a + 1; b < list->count; b++) {
                /* Check if opposite senses */
                if (list->refs[a].sense == list->refs[b].sense) continue;

                uint32_t cell_a = list->refs[a].cell_index;
                uint32_t cell_b = list->refs[b].cell_index;

                /* Only pair cells in the same universe */
                if (sys->cells.data[cell_a].universe_id !=
                    sys->cells.data[cell_b].universe_id) continue;
                int8_t sense_a = list->refs[a].sense;
                int8_t sense_b = list->refs[b].sense;

                /* Add B as neighbor of A */
                alea_cell_entry_t* ca = &sys->cells.data[cell_a];
                size_t na = ca->neighbor_count;
                alea_cell_neighbor_t* new_neighbors_a = realloc(ca->neighbors,
                    (na + 1) * sizeof(alea_cell_neighbor_t));
                if (new_neighbors_a) {
                    ca->neighbors = new_neighbors_a;
                    ca->neighbors[na].surface_id = sys->surfaces.data[si].mcnp_surface_id;
                    ca->neighbors[na].surface_index = (uint32_t)si;
                    ca->neighbors[na].neighbor_cell_id = sys->cells.data[cell_b].mcnp_cell_id;
                    ca->neighbors[na].neighbor_index = cell_b;
                    ca->neighbors[na].our_sense = sense_a;
                    ca->neighbor_count++;
                }

                /* Add A as neighbor of B */
                alea_cell_entry_t* cb = &sys->cells.data[cell_b];
                size_t nb = cb->neighbor_count;
                alea_cell_neighbor_t* new_neighbors_b = realloc(cb->neighbors,
                    (nb + 1) * sizeof(alea_cell_neighbor_t));
                if (new_neighbors_b) {
                    cb->neighbors = new_neighbors_b;
                    cb->neighbors[nb].surface_id = sys->surfaces.data[si].mcnp_surface_id;
                    cb->neighbors[nb].surface_index = (uint32_t)si;
                    cb->neighbors[nb].neighbor_cell_id = sys->cells.data[cell_a].mcnp_cell_id;
                    cb->neighbors[nb].neighbor_index = cell_a;
                    cb->neighbors[nb].our_sense = sense_b;
                    cb->neighbor_count++;
                }
            }
        }
    }

    sys->cell_adjacency_built = true;

cleanup:
    /* Free surface map */
    for (size_t i = 0; i < alea_vec_count(&sys->surfaces); i++) {
        free(surf_map[i].refs);
    }
    free(surf_map);

    return sys->cell_adjacency_built ? 0 : -1;
}

int alea_find_neighbor_cell(const alea_system_t* sys,
                           uint32_t cell_index,
                           int surface_id) {
    if (!sys || cell_index >= alea_vec_count(&sys->cells)) return -1;
    if (!sys->cell_adjacency_built) return -1;

    const alea_cell_entry_t* cell = &sys->cells.data[cell_index];

    for (size_t i = 0; i < cell->neighbor_count; i++) {
        if (cell->neighbors[i].surface_id == surface_id) {
            return (int)cell->neighbors[i].neighbor_index;
        }
    }

    return -1;  /* No neighbor found (exterior/void) */
}