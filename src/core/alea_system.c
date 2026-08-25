// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

#include "alea_system.h"
#include "alea_types.h"
#include "alea.h"
#include "core/alea_eval.h"
#include "core/alea_ops.h"
#include "core/alea_universe.h"
#include "core/alea_spatial_hier.h"
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

static atomic_uint_fast64_t g_next_system_id = 1;

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
    .export_materials = true,
    .export_transforms = true,
    .universe_depth = -1,
    .fill_depth = 0,

    /* Void generation */
    .void_max_depth = 8,
    .void_min_size = 0.1,
    .void_probes_per_axis = 3,

    /* Void merge */
    .merge_cell_weight = 1.0,
    .merge_surface_weight = 0.1,
    .merge_max_surfaces = 24,
    .merge_min_cells = 1,
    .void_consolidate = 100,

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
    sys->next_auto_material_id = 1; // Start at 1 for auto-assigned material IDs
    sys->source = ALEA_SOURCE_EMPTY;
    sys->system_id = atomic_fetch_add(&g_next_system_id, 1);
    atomic_init(&sys->geometry_generation, 1);
    atomic_init(&sys->query_cache_state, 0);
    atomic_flag_clear(&sys->query_cache_build_lock);
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

    // BVH acceleration (prepared on demand by query acceleration setup)
    sys->surface_bvh = NULL;
    sys->bvh_dirty = true;

    // Track initial memory (all vectors start empty)
    sys->stats.current_memory = 0;
    sys->stats.peak_memory = 0;

    return sys;
}

int alea_system_query_cache_ready(const alea_system_t* sys, unsigned flags) {
    if (!sys) return 0;
    if (flags == ALEA_CACHE_ALL) flags = ALEA_CACHE_RAYCAST | ALEA_CACHE_UNIVERSE;
    unsigned state = atomic_load(&sys->query_cache_state);
    return (state & flags) == flags;
}

uint64_t alea_system_geometry_generation(const alea_system_t* sys) {
    if (!sys) return 0;
    return atomic_load(&sys->geometry_generation);
}

static void alea_free_surface_lookup(alea_system_t* sys) {
    free(sys->surface_lookup);
    sys->surface_lookup = NULL;
    sys->surface_lookup_size = 0;
}

static void alea_free_cell_dynamic_fields(alea_system_t* sys) {
    for (size_t i = 0; i < alea_vec_count(&sys->cells); i++) {
        free(sys->cells.data[i].surface_indices);
        sys->cells.data[i].surface_indices = NULL;
        sys->cells.data[i].surface_index_count = 0;
        if (!sys->neighbor_pool)
            free(sys->cells.data[i].neighbors);
        sys->cells.data[i].neighbors = NULL;
        sys->cells.data[i].neighbor_count = 0;
        free(sys->cells.data[i].lat_fill);
        sys->cells.data[i].lat_fill = NULL;
        sys->cells.data[i].lat_fill_count = 0;
        free(sys->cells.data[i].comments);
        sys->cells.data[i].comments = NULL;
        free(sys->cells.data[i].inline_comment);
        sys->cells.data[i].inline_comment = NULL;
    }
    free(sys->neighbor_pool);
    sys->neighbor_pool = NULL;
    free(sys->surface_cell_offsets);
    sys->surface_cell_offsets = NULL;
    free(sys->surface_cell_refs);
    sys->surface_cell_refs = NULL;
    sys->surface_cell_ref_count = 0;
    sys->cell_adjacency_built = false;
}

static void alea_free_material_contents(alea_material_t* m) {
    for (size_t j = 0; j < alea_vec_count(&m->nuclides); j++) {
        free(m->nuclides.data[j].library);
    }
    alea_vec_free(&m->nuclides);
    for (size_t j = 0; j < alea_vec_count(&m->elements); j++) {
        free(m->elements.data[j].library);
    }
    alea_vec_free(&m->elements);
    for (size_t j = 0; j < alea_vec_count(&m->thermal_laws); j++) {
        free(m->thermal_laws.data[j].identifier);
    }
    alea_vec_free(&m->thermal_laws);
    free(m->name);
    free(m->comments);
}

static void alea_free_all_material_contents(alea_system_t* sys) {
    for (size_t i = 0; i < alea_vec_count(&sys->materials); i++) {
        alea_free_material_contents(&sys->materials.data[i]);
    }
}

static void alea_free_mixture_contents(alea_mixture_t* m) {
    alea_vec_free(&m->components);
    free(m->name);
    free(m->comments);
}

static void alea_free_all_mixture_contents(alea_system_t* sys) {
    for (size_t i = 0; i < alea_vec_count(&sys->mixtures); i++) {
        alea_free_mixture_contents(&sys->mixtures.data[i]);
    }
}

static void alea_clear_universe_cache(alea_system_t* sys, bool free_vector) {
    for (size_t i = 0; i < alea_vec_count(&sys->universes); i++) {
        alea_vec_free(&sys->universes.data[i].cell_indices);
        alea_vec_free(&sys->universes.data[i].point_bvh_nodes);
        free(sys->universes.data[i].point_bvh_indices);
        sys->universes.data[i].point_bvh_indices = NULL;
        sys->universes.data[i].point_bvh_built = false;
        sys->universes.data[i].point_bvh_disabled = false;
    }
    if (free_vector)
        alea_vec_free(&sys->universes);
    else
        alea_vec_clear(&sys->universes);
    universe_hashmap_clear(&sys->universe_index);
    sys->universe_index_built = false;
}

static void alea_free_query_cache_storage(alea_system_t* sys, unsigned flags,
                                          bool free_universe_vector) {
    if ((flags & ALEA_CACHE_HIER_SPATIAL) &&
        sys->volume_path_index) {
        alea_volume_path_index_free(sys->volume_path_index);
        sys->volume_path_index = NULL;
    }

    if (flags & ALEA_CACHE_UNIVERSE) {
        alea_clear_universe_cache(sys, free_universe_vector);
    }

    if (flags & ALEA_CACHE_CELL_SURFACES) {
        for (size_t i = 0; i < alea_vec_count(&sys->cells); i++) {
            free(sys->cells.data[i].surface_indices);
            sys->cells.data[i].surface_indices = NULL;
            sys->cells.data[i].surface_index_count = 0;
        }
        free(sys->prim_to_surface);
        sys->prim_to_surface = NULL;
        sys->prim_to_surface_size = 0;
        free(sys->mc_id_to_surface);
        sys->mc_id_to_surface = NULL;
        sys->mc_id_to_surface_size = 0;
    }

    if (flags & ALEA_CACHE_HIER_SPATIAL) {
        if (sys->hier_spatial_index) {
            alea_hier_spatial_index_free(sys->hier_spatial_index);
            sys->hier_spatial_index = NULL;
        }
    }

    if (flags & ALEA_CACHE_SURFACE_BVH) {
        if (sys->surface_bvh) {
            alea_bvh_free(sys->surface_bvh);
            sys->surface_bvh = NULL;
        }
        sys->bvh_dirty = true;
    }

    if (flags & ALEA_CACHE_ADJACENCY) {
        free(sys->neighbor_pool);
        sys->neighbor_pool = NULL;
        free(sys->surface_cell_offsets);
        sys->surface_cell_offsets = NULL;
        free(sys->surface_cell_refs);
        sys->surface_cell_refs = NULL;
        sys->surface_cell_ref_count = 0;
        for (size_t i = 0; i < alea_vec_count(&sys->cells); i++) {
            sys->cells.data[i].neighbors = NULL;
            sys->cells.data[i].neighbor_count = 0;
        }
        sys->cell_adjacency_built = false;
    }
}

void alea_system_invalidate_query_caches(alea_system_t* sys, unsigned flags) {
    if (!sys || flags == 0) return;

    unsigned invalidated_flags = flags;

    unsigned prev_state = atomic_fetch_and(&sys->query_cache_state, ~flags);
    atomic_fetch_add(&sys->geometry_generation, 1);

    if ((invalidated_flags & ALEA_CACHE_HIER_SPATIAL) &&
        sys->volume_path_index) {
        alea_volume_path_index_free(sys->volume_path_index);
        sys->volume_path_index = NULL;
    }

    /* Restrict expensive teardown to caches that were actually built.
     * During bulk load, callers invalidate on every add; without this,
     * each add walks all prior cells/universes (O(N^2) load time). */
    flags &= prev_state;
    if (flags == 0) return;

    alea_free_query_cache_storage(sys, flags, false);
}

int alea_system_prepare_query_caches(alea_system_t* sys, unsigned flags) {
    if (!sys) {
        alea_set_error_detail(ALEA_ERR_NULL_ARG,
                              "alea_system_prepare_query_caches: system is NULL");
        return -1;
    }

    if (flags == ALEA_CACHE_ALL) {
        flags = ALEA_CACHE_UNIVERSE | ALEA_CACHE_RAYCAST;
    }
    if (flags & ALEA_CACHE_HIER_SPATIAL)
        flags |= ALEA_CACHE_UNIVERSE | ALEA_CACHE_CELL_SURFACES;
    if (flags & ALEA_CACHE_ADJACENCY)
        flags |= ALEA_CACHE_CELL_SURFACES;

    if (alea_system_query_cache_ready(sys, flags))
        return 0;

    while (atomic_flag_test_and_set(&sys->query_cache_build_lock)) {
        /* Another thread is preparing shared caches. */
    }

    if (alea_system_query_cache_ready(sys, flags)) {
        atomic_flag_clear(&sys->query_cache_build_lock);
        return 0;
    }

    if ((flags & ALEA_CACHE_UNIVERSE) &&
        !alea_system_query_cache_ready(sys, ALEA_CACHE_UNIVERSE)) {
        if (alea_build_universe_index(sys) != 0) goto fail;
        atomic_fetch_or(&sys->query_cache_state, ALEA_CACHE_UNIVERSE);
    }
    if ((flags & ALEA_CACHE_CELL_SURFACES) &&
        !alea_system_query_cache_ready(sys, ALEA_CACHE_CELL_SURFACES)) {
        if (alea_build_cell_surface_index(sys) != 0) goto fail;
        atomic_fetch_or(&sys->query_cache_state, ALEA_CACHE_CELL_SURFACES);
    }

    if ((flags & ALEA_CACHE_HIER_SPATIAL) &&
        !alea_system_query_cache_ready(sys, ALEA_CACHE_HIER_SPATIAL)) {
        if (alea_hier_spatial_index_build(sys) != 0) goto fail;
        atomic_fetch_or(&sys->query_cache_state, ALEA_CACHE_HIER_SPATIAL);
    }

    if ((flags & ALEA_CACHE_SURFACE_BVH) &&
        !alea_system_query_cache_ready(sys, ALEA_CACHE_SURFACE_BVH)) {
        if (sys->surface_bvh) {
            alea_bvh_free(sys->surface_bvh);
            sys->surface_bvh = NULL;
        }
        if (alea_vec_count(&sys->surfaces) > 0) {
            sys->surface_bvh = alea_bvh_build(sys);
            if (!sys->surface_bvh) goto fail;
        }
        sys->bvh_dirty = false;
        atomic_fetch_or(&sys->query_cache_state, ALEA_CACHE_SURFACE_BVH);
    }

    if ((flags & ALEA_CACHE_ADJACENCY) &&
        !alea_system_query_cache_ready(sys, ALEA_CACHE_ADJACENCY)) {
        if (alea_build_cell_adjacency(sys) != 0) goto fail;
        atomic_fetch_or(&sys->query_cache_state, ALEA_CACHE_ADJACENCY);
    }

    atomic_flag_clear(&sys->query_cache_build_lock);
    return 0;

fail:
    atomic_flag_clear(&sys->query_cache_build_lock);
    return -1;
}

int alea_prepare_query_acceleration(alea_system_t* sys) {
    return alea_system_prepare_query_caches(sys, ALEA_CACHE_ALL);
}

void alea_system_destroy_internals(alea_system_t* sys) {
    if (!sys) return;

    atomic_store(&sys->query_cache_state, 0);
    atomic_fetch_add(&sys->geometry_generation, 1);
    alea_free_query_cache_storage(sys, ALEA_CACHE_ALL, true);
    alea_free_cell_dynamic_fields(sys);
    alea_vec_free(&sys->cells);

    alea_vec_free(&sys->nodes);
    alea_vec_free(&sys->primitives);
    alea_vec_free(&sys->primitive_planes);
    alea_vec_free(&sys->primitive_spheres);
    alea_vec_free(&sys->primitive_cyl_x);
    alea_vec_free(&sys->primitive_cyl_y);
    alea_vec_free(&sys->primitive_cyl_z);
    alea_vec_free(&sys->primitive_cone_x);
    alea_vec_free(&sys->primitive_cone_y);
    alea_vec_free(&sys->primitive_cone_z);
    alea_vec_free(&sys->primitive_boxes);
    alea_vec_free(&sys->primitive_quadrics);
    alea_vec_free(&sys->primitive_toruses);
    alea_vec_free(&sys->primitive_rccs);
    alea_vec_free(&sys->primitive_box_generals);
    alea_vec_free(&sys->primitive_sphs);
    alea_vec_free(&sys->primitive_trcs);
    alea_vec_free(&sys->primitive_ells);
    alea_vec_free(&sys->primitive_recs);
    alea_vec_free(&sys->primitive_weds);
    alea_vec_free(&sys->primitive_rhps);
    alea_vec_free(&sys->primitive_arbs);
    alea_vec_free(&sys->surfaces);
    alea_free_all_material_contents(sys);
    alea_vec_free(&sys->materials);
    alea_vec_free(&sys->transforms);
    alea_free_surface_lookup(sys);

    if (sys->primitive_index) {
        primitive_hash_table_destroy(sys->primitive_index);
    }

    cell_hashmap_destroy(&sys->cell_index);
    universe_hashmap_destroy(&sys->universe_index);

    alea_free_cell_refs(sys);

    alea_free_all_mixture_contents(sys);
    alea_vec_free(&sys->mixtures);
}

void alea_system_destroy(alea_system_t* sys) {
    if (!sys) return;
    alea_system_destroy_internals(sys);
    free(sys);
}

void alea_system_reset(alea_system_t* sys) {
    if (!sys) return;

    atomic_store(&sys->query_cache_state, 0);
    atomic_fetch_add(&sys->geometry_generation, 1);
    alea_free_query_cache_storage(sys, ALEA_CACHE_ALL, false);
    alea_free_cell_dynamic_fields(sys);
    alea_vec_clear(&sys->cells);

    alea_vec_clear(&sys->nodes);
    alea_vec_clear(&sys->primitives);
    alea_vec_clear(&sys->primitive_planes);
    alea_vec_clear(&sys->primitive_spheres);
    alea_vec_clear(&sys->primitive_cyl_x);
    alea_vec_clear(&sys->primitive_cyl_y);
    alea_vec_clear(&sys->primitive_cyl_z);
    alea_vec_clear(&sys->primitive_cone_x);
    alea_vec_clear(&sys->primitive_cone_y);
    alea_vec_clear(&sys->primitive_cone_z);
    alea_vec_clear(&sys->primitive_boxes);
    alea_vec_clear(&sys->primitive_quadrics);
    alea_vec_clear(&sys->primitive_toruses);
    alea_vec_clear(&sys->primitive_rccs);
    alea_vec_clear(&sys->primitive_box_generals);
    alea_vec_clear(&sys->primitive_sphs);
    alea_vec_clear(&sys->primitive_trcs);
    alea_vec_clear(&sys->primitive_ells);
    alea_vec_clear(&sys->primitive_recs);
    alea_vec_clear(&sys->primitive_weds);
    alea_vec_clear(&sys->primitive_rhps);
    alea_vec_clear(&sys->primitive_arbs);
    alea_vec_clear(&sys->surfaces);
    alea_free_all_material_contents(sys);
    alea_vec_clear(&sys->materials);
    alea_vec_clear(&sys->transforms);
    sys->next_inline_transform_id = 1;
    sys->next_auto_surface_id = 1;
    sys->next_auto_cell_id = 1;
    sys->next_auto_material_id = 1;

    alea_free_surface_lookup(sys);
    alea_free_all_mixture_contents(sys);
    alea_vec_clear(&sys->mixtures);

    // Clear hash tables
    if (sys->primitive_index) {
        primitive_hash_table_destroy(sys->primitive_index);
        sys->primitive_index = primitive_hash_table_create();
    }
    cell_hashmap_clear(&sys->cell_index);
    universe_hashmap_clear(&sys->universe_index);
    
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
                                                   src->primitive.mc_surface_id);
    if (new_id == ALEA_NODE_ID_INVALID) return ALEA_NODE_ID_INVALID;

    /* Copy material_id from source (alea_add_primitive_node doesn't set this) */
    sys->nodes.data[new_id].material_id = src->material_id;

    return new_id;
}

// ============================================================================
// PRIMITIVE OPERATIONS WITH DEDUPLICATION
// ============================================================================

static uint32_t alea_store_primitive_payload(alea_system_t* sys,
                                             alea_primitive_type_t type,
                                             const alea_primitive_data_t* data) {
#define STORE_PAYLOAD(kind, vec, field) do {                                  \
        uint32_t idx = (uint32_t)alea_vec_count(&(sys)->vec);                  \
        kind* dst = alea_vec_push_uninit(&(sys)->vec, kind);                   \
        if (!dst) return UINT32_MAX;                                           \
        *dst = data->field;                                                    \
        return idx;                                                            \
    } while (0)

    switch (type) {
        case ALEA_PRIMITIVE_PLANE:      STORE_PAYLOAD(alea_plane_data_t,       primitive_planes,       plane);
        case ALEA_PRIMITIVE_SPHERE:     STORE_PAYLOAD(alea_sphere_data_t,      primitive_spheres,      sphere);
        case ALEA_PRIMITIVE_CYLINDER_X: STORE_PAYLOAD(alea_cylinder_x_data_t,  primitive_cyl_x,        cyl_x);
        case ALEA_PRIMITIVE_CYLINDER_Y: STORE_PAYLOAD(alea_cylinder_y_data_t,  primitive_cyl_y,        cyl_y);
        case ALEA_PRIMITIVE_CYLINDER_Z: STORE_PAYLOAD(alea_cylinder_z_data_t,  primitive_cyl_z,        cyl_z);
        case ALEA_PRIMITIVE_CONE_X:     STORE_PAYLOAD(alea_cone_x_data_t,      primitive_cone_x,       cone_x);
        case ALEA_PRIMITIVE_CONE_Y:     STORE_PAYLOAD(alea_cone_y_data_t,      primitive_cone_y,       cone_y);
        case ALEA_PRIMITIVE_CONE_Z:     STORE_PAYLOAD(alea_cone_z_data_t,      primitive_cone_z,       cone_z);
        case ALEA_PRIMITIVE_RPP:        STORE_PAYLOAD(alea_box_data_t,         primitive_boxes,        box);
        case ALEA_PRIMITIVE_QUADRIC:    STORE_PAYLOAD(alea_quadric_data_t,     primitive_quadrics,     quadric);
        case ALEA_PRIMITIVE_TORUS_X:
        case ALEA_PRIMITIVE_TORUS_Y:
        case ALEA_PRIMITIVE_TORUS_Z:    STORE_PAYLOAD(alea_torus_data_t,       primitive_toruses,      torus);
        case ALEA_PRIMITIVE_RCC:        STORE_PAYLOAD(alea_rcc_data_t,         primitive_rccs,         rcc);
        case ALEA_PRIMITIVE_BOX:        STORE_PAYLOAD(alea_box_general_data_t, primitive_box_generals, box_general);
        case ALEA_PRIMITIVE_SPH:        STORE_PAYLOAD(alea_sph_data_t,         primitive_sphs,         sph);
        case ALEA_PRIMITIVE_TRC:        STORE_PAYLOAD(alea_trc_data_t,         primitive_trcs,         trc);
        case ALEA_PRIMITIVE_ELL:        STORE_PAYLOAD(alea_ell_data_t,         primitive_ells,         ell);
        case ALEA_PRIMITIVE_REC:        STORE_PAYLOAD(alea_rec_data_t,         primitive_recs,         rec);
        case ALEA_PRIMITIVE_WED:        STORE_PAYLOAD(alea_wed_data_t,         primitive_weds,         wed);
        case ALEA_PRIMITIVE_RHP:        STORE_PAYLOAD(alea_rhp_data_t,         primitive_rhps,         rhp);
        case ALEA_PRIMITIVE_ARB:        STORE_PAYLOAD(alea_arb_data_t,         primitive_arbs,         arb);
        default: return UINT32_MAX;
    }

#undef STORE_PAYLOAD
}

const void* alea_primitive_payload_const(const alea_system_t* sys,
                                         uint32_t id) {
    if (!sys || id >= alea_vec_count(&sys->primitives)) return NULL;
    const alea_primitive_entry_t* prim = &sys->primitives.data[id];
    uint32_t idx = prim->payload_index;

#define GET_PAYLOAD(vec) \
    ((idx < alea_vec_count(&(sys)->vec)) ? (const void*)&(sys)->vec.data[idx] : NULL)

    switch (prim->type) {
        case ALEA_PRIMITIVE_PLANE:       return GET_PAYLOAD(primitive_planes);
        case ALEA_PRIMITIVE_SPHERE:      return GET_PAYLOAD(primitive_spheres);
        case ALEA_PRIMITIVE_CYLINDER_X:  return GET_PAYLOAD(primitive_cyl_x);
        case ALEA_PRIMITIVE_CYLINDER_Y:  return GET_PAYLOAD(primitive_cyl_y);
        case ALEA_PRIMITIVE_CYLINDER_Z:  return GET_PAYLOAD(primitive_cyl_z);
        case ALEA_PRIMITIVE_CONE_X:      return GET_PAYLOAD(primitive_cone_x);
        case ALEA_PRIMITIVE_CONE_Y:      return GET_PAYLOAD(primitive_cone_y);
        case ALEA_PRIMITIVE_CONE_Z:      return GET_PAYLOAD(primitive_cone_z);
        case ALEA_PRIMITIVE_RPP:         return GET_PAYLOAD(primitive_boxes);
        case ALEA_PRIMITIVE_QUADRIC:     return GET_PAYLOAD(primitive_quadrics);
        case ALEA_PRIMITIVE_TORUS_X:
        case ALEA_PRIMITIVE_TORUS_Y:
        case ALEA_PRIMITIVE_TORUS_Z:     return GET_PAYLOAD(primitive_toruses);
        case ALEA_PRIMITIVE_RCC:         return GET_PAYLOAD(primitive_rccs);
        case ALEA_PRIMITIVE_BOX:         return GET_PAYLOAD(primitive_box_generals);
        case ALEA_PRIMITIVE_SPH:         return GET_PAYLOAD(primitive_sphs);
        case ALEA_PRIMITIVE_TRC:         return GET_PAYLOAD(primitive_trcs);
        case ALEA_PRIMITIVE_ELL:         return GET_PAYLOAD(primitive_ells);
        case ALEA_PRIMITIVE_REC:         return GET_PAYLOAD(primitive_recs);
        case ALEA_PRIMITIVE_WED:         return GET_PAYLOAD(primitive_weds);
        case ALEA_PRIMITIVE_RHP:         return GET_PAYLOAD(primitive_rhps);
        case ALEA_PRIMITIVE_ARB:         return GET_PAYLOAD(primitive_arbs);
        default: return NULL;
    }

#undef GET_PAYLOAD
}

bool alea_primitive_copy_data(const alea_system_t* sys,
                              uint32_t id,
                              alea_primitive_data_t* out) {
    if (!out) return false;
    memset(out, 0, sizeof(*out));
    if (!sys || id >= alea_vec_count(&sys->primitives)) return false;

    const alea_primitive_entry_t* prim = &sys->primitives.data[id];
    const void* payload = alea_primitive_payload_const(sys, id);
    if (!payload) return false;

    switch (prim->type) {
        case ALEA_PRIMITIVE_PLANE:       out->plane = *(const alea_plane_data_t*)payload; return true;
        case ALEA_PRIMITIVE_SPHERE:      out->sphere = *(const alea_sphere_data_t*)payload; return true;
        case ALEA_PRIMITIVE_CYLINDER_X:  out->cyl_x = *(const alea_cylinder_x_data_t*)payload; return true;
        case ALEA_PRIMITIVE_CYLINDER_Y:  out->cyl_y = *(const alea_cylinder_y_data_t*)payload; return true;
        case ALEA_PRIMITIVE_CYLINDER_Z:  out->cyl_z = *(const alea_cylinder_z_data_t*)payload; return true;
        case ALEA_PRIMITIVE_CONE_X:      out->cone_x = *(const alea_cone_x_data_t*)payload; return true;
        case ALEA_PRIMITIVE_CONE_Y:      out->cone_y = *(const alea_cone_y_data_t*)payload; return true;
        case ALEA_PRIMITIVE_CONE_Z:      out->cone_z = *(const alea_cone_z_data_t*)payload; return true;
        case ALEA_PRIMITIVE_RPP:         out->box = *(const alea_box_data_t*)payload; return true;
        case ALEA_PRIMITIVE_QUADRIC:     out->quadric = *(const alea_quadric_data_t*)payload; return true;
        case ALEA_PRIMITIVE_TORUS_X:
        case ALEA_PRIMITIVE_TORUS_Y:
        case ALEA_PRIMITIVE_TORUS_Z:     out->torus = *(const alea_torus_data_t*)payload; return true;
        case ALEA_PRIMITIVE_RCC:         out->rcc = *(const alea_rcc_data_t*)payload; return true;
        case ALEA_PRIMITIVE_BOX:         out->box_general = *(const alea_box_general_data_t*)payload; return true;
        case ALEA_PRIMITIVE_SPH:         out->sph = *(const alea_sph_data_t*)payload; return true;
        case ALEA_PRIMITIVE_TRC:         out->trc = *(const alea_trc_data_t*)payload; return true;
        case ALEA_PRIMITIVE_ELL:         out->ell = *(const alea_ell_data_t*)payload; return true;
        case ALEA_PRIMITIVE_REC:         out->rec = *(const alea_rec_data_t*)payload; return true;
        case ALEA_PRIMITIVE_WED:         out->wed = *(const alea_wed_data_t*)payload; return true;
        case ALEA_PRIMITIVE_RHP:         out->rhp = *(const alea_rhp_data_t*)payload; return true;
        case ALEA_PRIMITIVE_ARB:         out->arb = *(const alea_arb_data_t*)payload; return true;
        default: return false;
    }
}

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
        // Found duplicate. Reusing a primitive does not change node ownership;
        // ref_count tracks only primitive nodes that actually reference it.
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
                    orig_surf_id = sys->surfaces.data[i].mc_surface_id;
                    break;
                }
            }
            ALEA_LOG_DEBUG("Dedup hit: primitive %u (%s) reused, original surface %d (ref_count=%u)",
                          existing, alea_primitive_type_name(type),
                          orig_surf_id, sys->primitives.data[existing].ref_count);
        }

        return existing;
    }

    uint32_t payload_index = alea_store_primitive_payload(sys, type, data);
    if (payload_index == UINT32_MAX) {
        alea_set_error_detail(ALEA_ERR_OUT_OF_MEMORY,
                              "alea_get_or_create_primitive: failed to allocate primitive payload %u (type %d)",
                              (uint32_t)alea_vec_count(&sys->primitives), type);
        return UINT32_MAX;
    }

    // Create new primitive using vector API
    uint32_t id = (uint32_t)alea_vec_count(&sys->primitives);
    alea_primitive_entry_t* prim = alea_vec_push_uninit(&sys->primitives, alea_primitive_entry_t);
    if (!prim) {
        alea_set_error_detail(ALEA_ERR_OUT_OF_MEMORY, "alea_get_or_create_primitive: failed to allocate primitive %u (type %d)", id, type);
        return UINT32_MAX;
    }

    prim->type = type;
    prim->payload_index = payload_index;
    prim->ref_count = 0;

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
                                      int32_t mc_surface_id) {
    if (!sys || primitive_id >= alea_vec_count(&sys->primitives)) return ALEA_NODE_ID_INVALID;

    alea_node_id_t node_id = alea_alloc_node(sys);
    if (node_id == ALEA_NODE_ID_INVALID) return ALEA_NODE_ID_INVALID;

    alea_node_t* node = &sys->nodes.data[node_id];
    ALEA_SET_OPERATION(node, ALEA_OP_PRIMITIVE);
    node->primitive.primitive_id = primitive_id;
    node->primitive.prim_type = sys->primitives.data[primitive_id].type;
    node->primitive.sense = sense;
    node->primitive.inverted = inverted;
    node->primitive.mc_surface_id = mc_surface_id;

    /* Compute bounding box (sense-aware for proper halfspace bounds) */
    /* Note: inverted flips the effective sense */
    int8_t effective_sense = inverted ? -sense : sense;
    alea_primitive_data_t data;
    if (!alea_primitive_copy_data(sys, primitive_id, &data)) {
        return ALEA_NODE_ID_INVALID;
    }
    alea_bbox_t prim_bbox = alea_halfspace_bbox(sys->primitives.data[primitive_id].type,
                                     &data, effective_sense);
    alea_node_bbox_set(&node->bbox, &prim_bbox);

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
    total += sys->primitive_planes.capacity * sizeof(alea_plane_data_t);
    total += sys->primitive_spheres.capacity * sizeof(alea_sphere_data_t);
    total += sys->primitive_cyl_x.capacity * sizeof(alea_cylinder_x_data_t);
    total += sys->primitive_cyl_y.capacity * sizeof(alea_cylinder_y_data_t);
    total += sys->primitive_cyl_z.capacity * sizeof(alea_cylinder_z_data_t);
    total += sys->primitive_cone_x.capacity * sizeof(alea_cone_x_data_t);
    total += sys->primitive_cone_y.capacity * sizeof(alea_cone_y_data_t);
    total += sys->primitive_cone_z.capacity * sizeof(alea_cone_z_data_t);
    total += sys->primitive_boxes.capacity * sizeof(alea_box_data_t);
    total += sys->primitive_quadrics.capacity * sizeof(alea_quadric_data_t);
    total += sys->primitive_toruses.capacity * sizeof(alea_torus_data_t);
    total += sys->primitive_rccs.capacity * sizeof(alea_rcc_data_t);
    total += sys->primitive_box_generals.capacity * sizeof(alea_box_general_data_t);
    total += sys->primitive_sphs.capacity * sizeof(alea_sph_data_t);
    total += sys->primitive_trcs.capacity * sizeof(alea_trc_data_t);
    total += sys->primitive_ells.capacity * sizeof(alea_ell_data_t);
    total += sys->primitive_recs.capacity * sizeof(alea_rec_data_t);
    total += sys->primitive_weds.capacity * sizeof(alea_wed_data_t);
    total += sys->primitive_rhps.capacity * sizeof(alea_rhp_data_t);
    total += sys->primitive_arbs.capacity * sizeof(alea_arb_data_t);
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
    if (sys->mc_id_to_surface) {
        total += sys->mc_id_to_surface_size * sizeof(uint32_t);
    }
    if (sys->surface_cell_offsets) {
        total += (alea_vec_count(&sys->surfaces) + 1) * sizeof(size_t);
    }
    if (sys->surface_cell_refs) {
        total += sys->surface_cell_ref_count * sizeof(alea_surface_cell_ref_t);
    }

    /* Per-material internal arrays */
    for (size_t i = 0; i < alea_vec_count(&sys->materials); i++) {
        const alea_material_t* m = &sys->materials.data[i];
        total += m->nuclides.capacity * sizeof(alea_nuclide_t);
        total += m->elements.capacity * sizeof(alea_element_comp_t);
        total += m->thermal_laws.capacity * sizeof(alea_thermal_law_t);
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
        total += u->cell_indices.capacity * sizeof(size_t);
    }

    return total;
}

void alea_system_shrink_to_fit(alea_system_t* sys) {
    if (!sys) return;
    /* Release the doubling-growth slack from the big arrays after loading.
     * The node array dominates (2x growth can leave ~50% unused tail). */
    alea_vec_shrink_to_fit(&sys->nodes, alea_node_t);
    alea_vec_shrink_to_fit(&sys->primitives, alea_primitive_entry_t);
    alea_vec_shrink_to_fit(&sys->primitive_planes, alea_plane_data_t);
    alea_vec_shrink_to_fit(&sys->primitive_spheres, alea_sphere_data_t);
    alea_vec_shrink_to_fit(&sys->primitive_cyl_x, alea_cylinder_x_data_t);
    alea_vec_shrink_to_fit(&sys->primitive_cyl_y, alea_cylinder_y_data_t);
    alea_vec_shrink_to_fit(&sys->primitive_cyl_z, alea_cylinder_z_data_t);
    alea_vec_shrink_to_fit(&sys->primitive_cone_x, alea_cone_x_data_t);
    alea_vec_shrink_to_fit(&sys->primitive_cone_y, alea_cone_y_data_t);
    alea_vec_shrink_to_fit(&sys->primitive_cone_z, alea_cone_z_data_t);
    alea_vec_shrink_to_fit(&sys->primitive_boxes, alea_box_data_t);
    alea_vec_shrink_to_fit(&sys->primitive_quadrics, alea_quadric_data_t);
    alea_vec_shrink_to_fit(&sys->primitive_toruses, alea_torus_data_t);
    alea_vec_shrink_to_fit(&sys->primitive_rccs, alea_rcc_data_t);
    alea_vec_shrink_to_fit(&sys->primitive_box_generals, alea_box_general_data_t);
    alea_vec_shrink_to_fit(&sys->primitive_sphs, alea_sph_data_t);
    alea_vec_shrink_to_fit(&sys->primitive_trcs, alea_trc_data_t);
    alea_vec_shrink_to_fit(&sys->primitive_ells, alea_ell_data_t);
    alea_vec_shrink_to_fit(&sys->primitive_recs, alea_rec_data_t);
    alea_vec_shrink_to_fit(&sys->primitive_weds, alea_wed_data_t);
    alea_vec_shrink_to_fit(&sys->primitive_rhps, alea_rhp_data_t);
    alea_vec_shrink_to_fit(&sys->primitive_arbs, alea_arb_data_t);
    alea_vec_shrink_to_fit(&sys->surfaces, alea_surface_entry_t);
    alea_vec_shrink_to_fit(&sys->cells, alea_cell_entry_t);
    alea_vec_shrink_to_fit(&sys->transforms, alea_transform_t);
    alea_vec_shrink_to_fit(&sys->materials, alea_material_t);
    alea_vec_shrink_to_fit(&sys->mixtures, alea_mixture_t);
    alea_vec_shrink_to_fit(&sys->cell_refs, alea_cell_ref_t);
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

static double local_vec3_dot(const double a[3], const double b[3]) {
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

static void local_vec3_cross(const double a[3], const double b[3], double out[3]) {
    out[0] = a[1]*b[2] - a[2]*b[1];
    out[1] = a[2]*b[0] - a[0]*b[2];
    out[2] = a[0]*b[1] - a[1]*b[0];
}

static bool local_vec3_normalize(double v[3]) {
    double n2 = local_vec3_dot(v, v);
    if (n2 <= 1e-24 || !isfinite(n2)) return false;
    double inv_n = 1.0 / sqrt(n2);
    v[0] *= inv_n;
    v[1] *= inv_n;
    v[2] *= inv_n;
    return true;
}

static bool complete_one_vector(const double* r, double R[9]) {
    double c0[3] = {r[0], r[1], r[2]};
    if (!local_vec3_normalize(c0)) return false;

    double helper[3] = {1.0, 0.0, 0.0};
    if (fabs(c0[0]) < fabs(c0[1]) && fabs(c0[0]) < fabs(c0[2])) {
        helper[0] = 1.0; helper[1] = 0.0; helper[2] = 0.0;
    } else if (fabs(c0[1]) < fabs(c0[2])) {
        helper[0] = 0.0; helper[1] = 1.0; helper[2] = 0.0;
    } else {
        helper[0] = 0.0; helper[1] = 0.0; helper[2] = 1.0;
    }

    double c1[3], c2[3];
    local_vec3_cross(helper, c0, c1);
    if (!local_vec3_normalize(c1)) return false;
    local_vec3_cross(c0, c1, c2);
    if (!local_vec3_normalize(c2)) return false;

    R[0] = c0[0]; R[1] = c0[1]; R[2] = c0[2];
    R[3] = c1[0]; R[4] = c1[1]; R[5] = c1[2];
    R[6] = c2[0]; R[7] = c2[1]; R[8] = c2[2];
    return true;
}

static bool complete_two_vectors(const double* r, double R[9]) {
    double c0[3] = {r[0], r[1], r[2]};
    double c1[3] = {r[3], r[4], r[5]};
    if (!local_vec3_normalize(c0) || !local_vec3_normalize(c1)) return false;
    if (fabs(local_vec3_dot(c0, c1)) > 1e-8) return false;

    double c2[3];
    local_vec3_cross(c0, c1, c2);
    if (!local_vec3_normalize(c2)) return false;

    R[0] = c0[0]; R[1] = c0[1]; R[2] = c0[2];
    R[3] = c1[0]; R[4] = c1[1]; R[5] = c1[2];
    R[6] = c2[0]; R[7] = c2[1]; R[8] = c2[2];
    return true;
}

static bool complete_five_values(const double* r, double R[9]) {
    double c0[3] = {r[0], r[1], r[2]};
    if (!local_vec3_normalize(c0)) return false;
    if (fabs(c0[2]) < 1e-12) return false;

    double c1[3] = {r[3], r[4], 0.0};
    c1[2] = -(c0[0]*c1[0] + c0[1]*c1[1]) / c0[2];
    if (!local_vec3_normalize(c1)) return false;
    if (fabs(local_vec3_dot(c0, c1)) > 1e-8) return false;

    double c2[3];
    local_vec3_cross(c0, c1, c2);
    if (!local_vec3_normalize(c2)) return false;

    R[0] = c0[0]; R[1] = c0[1]; R[2] = c0[2];
    R[3] = c1[0]; R[4] = c1[1]; R[5] = c1[2];
    R[6] = c2[0]; R[7] = c2[1]; R[8] = c2[2];
    return true;
}

int alea_normalize_mcnp_transform_values(const double* data, int value_count,
                                        int degrees, double out[12],
                                        int* out_count) {
    if (!data || !out || !out_count || value_count < 3 || value_count > 13) {
        alea_set_error_detail(ALEA_ERR_INVALID_ARG,
                "transform: value_count must be between 3 and 13, got %d",
                value_count);
        return -1;
    }

    for (int i = 0; i < value_count; i++) {
        if (!isfinite(data[i])) {
            alea_set_error_detail(ALEA_ERR_INVALID_ARG,
                    "transform: value %d is not finite", i);
            return -1;
        }
    }

    int m = 1;
    int count = value_count;
    int rot_count = count - 3;
    if ((count == 4 || count == 7 || count == 9 || count == 10 || count == 13) &&
        (data[count - 1] == 1.0 || data[count - 1] == -1.0)) {
        m = (int)data[count - 1];
        count--;
        rot_count = count - 3;
    }

    if (rot_count != 0 && rot_count != 3 && rot_count != 5 &&
        rot_count != 6 && rot_count != 9) {
        alea_set_error_detail(ALEA_ERR_INVALID_ARG,
                "transform: expected 0, 3, 5, 6, or 9 rotation entries, got %d",
                rot_count);
        return -1;
    }

    double R[9] = {1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0};
    double r[9] = {0};
    for (int i = 0; i < rot_count; i++) {
        r[i] = degrees ? cos(DEG_TO_RAD(data[3 + i])) : data[3 + i];
        if (!isfinite(r[i])) {
            alea_set_error_detail(ALEA_ERR_INVALID_ARG,
                    "transform: rotation entry %d is not finite", i);
            return -1;
        }
    }

    bool ok = true;
    if (rot_count == 9) {
        memcpy(R, r, sizeof(R));
    } else if (rot_count == 6) {
        ok = complete_two_vectors(r, R);
    } else if (rot_count == 5) {
        ok = complete_five_values(r, R);
    } else if (rot_count == 3) {
        ok = complete_one_vector(r, R);
    }
    if (!ok) {
        alea_set_error_detail(ALEA_ERR_INVALID_ARG,
                "transform: could not complete partial rotation matrix");
        return -1;
    }

    double d[3] = {data[0], data[1], data[2]};
    double T[3];
    if (m == -1) {
        T[0] = -(R[0]*d[0] + R[3]*d[1] + R[6]*d[2]);
        T[1] = -(R[1]*d[0] + R[4]*d[1] + R[7]*d[2]);
        T[2] = -(R[2]*d[0] + R[5]*d[1] + R[8]*d[2]);
    } else {
        T[0] = d[0];
        T[1] = d[1];
        T[2] = d[2];
    }

    out[0] = T[0]; out[1] = T[1]; out[2] = T[2];
    if (rot_count == 0) {
        *out_count = 3;
        for (int i = 3; i < 12; i++) out[i] = 0.0;
        return 0;
    }

    for (int i = 0; i < 9; i++) out[3 + i] = R[i];
    *out_count = 12;
    return 0;
}

static bool validate_transform_values(const char* caller, const char* context,
                                      const double* data, int value_count,
                                      int degrees) {
    if (value_count < 1 || value_count > 12 ||
        (value_count > 3 && value_count != 12)) {
        alea_set_error_detail(ALEA_ERR_INVALID_ARG,
                "%s: value_count must be 1-3 for translation or 12 for rotation, got %d",
                caller, value_count);
        return false;
    }

    for (int i = 0; i < value_count; i++) {
        if (!isfinite(data[i])) {
            alea_set_error_detail(ALEA_ERR_INVALID_ARG,
                    "%s: transform value %d is not finite", caller, i);
            return false;
        }
    }

    if (value_count <= 3) return true;

    double b[9];
    for (int i = 0; i < 9; i++) {
        b[i] = degrees ? cos(DEG_TO_RAD(data[3 + i])) : data[3 + i];
        if (!isfinite(b[i])) {
            alea_set_error_detail(ALEA_ERR_INVALID_ARG,
                    "%s: transform cosine %d is not finite", caller, i);
            return false;
        }
    }

    /* MCNP stores the three auxiliary basis vectors as columns:
     * [b1,b2,b3], [b4,b5,b6], [b7,b8,b9]. */
    double c0[3] = {b[0], b[1], b[2]};
    double c1[3] = {b[3], b[4], b[5]};
    double c2[3] = {b[6], b[7], b[8]};
    double n0 = c0[0]*c0[0] + c0[1]*c0[1] + c0[2]*c0[2];
    double n1 = c1[0]*c1[0] + c1[1]*c1[1] + c1[2]*c1[2];
    double n2 = c2[0]*c2[0] + c2[1]*c2[1] + c2[2]*c2[2];
    double d01 = c0[0]*c1[0] + c0[1]*c1[1] + c0[2]*c1[2];
    double d02 = c0[0]*c2[0] + c0[1]*c2[1] + c0[2]*c2[2];
    double d12 = c1[0]*c2[0] + c1[1]*c2[1] + c1[2]*c2[2];
    double det = c0[0]*(c1[1]*c2[2] - c1[2]*c2[1])
               - c1[0]*(c0[1]*c2[2] - c0[2]*c2[1])
               + c2[0]*(c0[1]*c1[2] - c0[2]*c1[1]);
    const double tol = 1e-5;

    if (fabs(n0 - 1.0) > tol || fabs(n1 - 1.0) > tol ||
        fabs(n2 - 1.0) > tol || fabs(d01) > tol ||
        fabs(d02) > tol || fabs(d12) > tol ||
        fabs(fabs(det) - 1.0) > tol) {
        alea_set_error_detail(ALEA_ERR_INVALID_ARG,
                "%s: transform rotation must be finite, orthonormal, and invertible "
                "(norms %.12g %.12g %.12g, dots %.12g %.12g %.12g, det %.12g)",
                caller, n0, n1, n2, d01, d02, d12, det);
        return false;
    }

    if (det < 0.0) {
        if (context && context[0] != '\0') {
            ALEA_LOG_INFO("%s: transform rotation for %s has negative determinant %.12g (reflection + rotation)",
                          caller, context, det);
        } else {
            ALEA_LOG_INFO("%s: transform rotation has negative determinant %.12g (reflection + rotation)",
                          caller, det);
        }
    }

    return true;
}

int alea_add_transform(alea_system_t* sys, int transform_id,
                      const double* data, int value_count, int degrees) {
    if (!sys || !data) {
        alea_set_error_detail(ALEA_ERR_NULL_ARG, "alea_add_transform: NULL argument");
        return -1;
    }
    double normalized[12];
    int normalized_count = 0;
    if (alea_normalize_mcnp_transform_values(data, value_count, degrees,
                                             normalized, &normalized_count) != 0) {
        return -1;
    }
    char context[32];
    snprintf(context, sizeof(context), "TR%d", transform_id);
    if (!validate_transform_values("alea_add_transform", context, normalized,
                                   normalized_count, 0)) {
        return -1;
    }

    // Check if transform already exists
    for (size_t i = 0; i < alea_vec_count(&sys->transforms); i++) {
        if (sys->transforms.data[i].transform_id == transform_id) {
            // Update existing transform
            alea_transform_t* tr = &sys->transforms.data[i];
            populate_transform_data(tr, normalized, normalized_count, 0);
            tr->value_count = normalized_count;
            tr->degrees = 0;
            alea_system_invalidate_query_caches(sys, ALEA_CACHE_ALL);
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
    tr->value_count = normalized_count;
    tr->degrees = 0;
    tr->from_inline = 0;

    populate_transform_data(tr, normalized, normalized_count, 0);

    ALEA_LOG_INFO("Added transform TR%d (%d values, %s)",
           transform_id, value_count,
           "normalized cosines");

    alea_system_invalidate_query_caches(sys, ALEA_CACHE_ALL);
    return 0;
}

int alea_add_inline_transform(alea_system_t* sys, const double* data,
                             int value_count, int degrees,
                             int cell_id, const char* role) {
    if (!sys || !data) {
        alea_set_error_detail(ALEA_ERR_NULL_ARG, "alea_add_inline_transform: NULL argument");
        return -1;
    }
    double normalized[12];
    int normalized_count = 0;
    if (alea_normalize_mcnp_transform_values(data, value_count, degrees,
                                             normalized, &normalized_count) != 0) {
        return -1;
    }
    char context[64];
    const char* context_ptr = NULL;
    if (cell_id != 0 && role && role[0] != '\0') {
        snprintf(context, sizeof(context), "cell %d inline %s", cell_id, role);
        context_ptr = context;
    } else if (cell_id != 0) {
        snprintf(context, sizeof(context), "cell %d inline transform", cell_id);
        context_ptr = context;
    }
    if (!validate_transform_values("alea_add_inline_transform", context_ptr, normalized,
                                   normalized_count, 0)) {
        return -1;
    }

    // Build candidate cosines for comparison
    alea_transform_t candidate;
    memset(&candidate, 0, sizeof(candidate));
    candidate.value_count = normalized_count;
    candidate.degrees = 0;
    populate_transform_data(&candidate, normalized, normalized_count, 0);

    // Deduplicate: scan existing inline transforms for a match on cosines
    static const double TRANSFORM_TOL = 1e-10;
    for (size_t i = 0; i < alea_vec_count(&sys->transforms); i++) {
        const alea_transform_t* existing = &sys->transforms.data[i];
        if (!existing->from_inline) continue;
        if (existing->value_count != normalized_count) continue;

        bool match = true;
        for (int j = 0; j < normalized_count; j++) {
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
           assigned_id, normalized_count,
           "normalized cosines");

    alea_system_invalidate_query_caches(sys, ALEA_CACHE_ALL);
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
    alea_vec_init(&u->cell_indices);
    int vres = alea_vec_reserve(&u->cell_indices, INITIAL_UNIVERSE_CELLS_CAPACITY, size_t);
    if (vres != 0) {
        alea_vec_pop_discard(&sys->universes);  // Rollback
        return NULL;
    }
    u->bbox = (alea_bbox_t){1e30, -1e30, 1e30, -1e30, 1e30, -1e30};  // Empty bbox

    universe_hashmap_put(&sys->universe_index, universe_id, new_index);
    return u;
}

static int add_cell_to_universe(alea_universe_t* u, size_t cell_index) {
    int res = alea_vec_push(&u->cell_indices, cell_index, size_t);
    return res != 0 ? -1 : 0;
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
        if (sys->cells.data[i].mc_cell_id == cell_id) {
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
                         int material_index, double density, int universe_id) {
    if (!sys) {
        alea_set_error_detail(ALEA_ERR_NULL_ARG, "alea_add_cell_with_id: system is NULL");
        return -1;
    }

    /* Resolve material index to MCNP material ID */
    int mat_id = 0;
    int mat_idx = -1;
    if (material_index == ALEA_MATERIAL_VOID || material_index < 0) {
        mat_id = 0;
        mat_idx = -1;
    } else {
        if ((size_t)material_index >= alea_vec_count(&sys->materials)) {
            alea_set_error_detail(ALEA_ERR_INVALID_ARG,
                "alea_add_cell_with_id: material index %d out of range (have %zu materials)",
                material_index, alea_vec_count(&sys->materials));
            return -1;
        }
        mat_id = sys->materials.data[material_index].material_id;
        mat_idx = material_index;
    }

    int idx = (int)alea_vec_count(&sys->cells);
    alea_cell_entry_t* cell = alea_vec_push_uninit(&sys->cells, alea_cell_entry_t);
    if (!cell) {
        alea_set_error_detail(ALEA_ERR_OUT_OF_MEMORY, "alea_add_cell_with_id: failed to allocate cell %d", cell_id);
        return -1;
    }
    memset(cell, 0, sizeof(*cell));

    cell->mc_cell_id = cell_id;
    cell->root_node_id = root_node;
    cell->original_root_node_id = ALEA_NODE_ID_INVALID;
    cell->material_id = mat_id;
    cell->material_index = mat_idx;
    cell->is_mass_density = (density < 0) ? 1 : 0;
    cell->density = fabs(density);
    cell->universe_id = universe_id;

    alea_system_invalidate_query_caches(sys, ALEA_CACHE_ALL);

    /* Update next_auto_cell_id if this ID is >= it */
    if (cell_id >= sys->next_auto_cell_id) {
        sys->next_auto_cell_id = cell_id + 1;
    }

    /* Insert into cell hash map */
    cell_hashmap_put(&sys->cell_index, cell_id, idx);

    /* Notify module hooks */
    if (sys->on_cell_added)
        sys->on_cell_added(sys->cell_hook_userdata, (size_t)idx);

    return idx;
}

/**
 * @brief Add cell with auto-assigned or validated ID
 *
 * If cell_id <= 0, auto-assigns a new ID.
 * If cell_id > 0 but already exists, logs warning and auto-assigns.
 */
int alea_add_cell(alea_system_t* sys, int cell_id, alea_node_id_t root_node,
                 int material_index, double density, int universe_id) {
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

    /* Resolve material index to MCNP material ID */
    int mat_id = 0;
    int mat_idx = -1;
    if (material_index == ALEA_MATERIAL_VOID || material_index < 0) {
        mat_id = 0;
        mat_idx = -1;
    } else {
        if ((size_t)material_index >= alea_vec_count(&sys->materials)) {
            alea_set_error_detail(ALEA_ERR_INVALID_ARG,
                "alea_add_cell: material index %d out of range (have %zu materials)",
                material_index, alea_vec_count(&sys->materials));
            return -1;
        }
        mat_id = sys->materials.data[material_index].material_id;
        mat_idx = material_index;
    }

    int idx = (int)alea_vec_count(&sys->cells);
    alea_cell_entry_t* cell = alea_vec_push_uninit(&sys->cells, alea_cell_entry_t);
    if (!cell) {
        alea_set_error_detail(ALEA_ERR_OUT_OF_MEMORY, "alea_add_cell: failed to allocate cell %d", final_cell_id);
        return -1;
    }
    memset(cell, 0, sizeof(*cell));

    cell->mc_cell_id = final_cell_id;
    cell->root_node_id = root_node;
    cell->original_root_node_id = ALEA_NODE_ID_INVALID;
    cell->material_id = mat_id;
    cell->material_index = mat_idx;
    cell->is_mass_density = (density < 0) ? 1 : 0;
    cell->density = fabs(density);
    cell->universe_id = universe_id;

    alea_system_invalidate_query_caches(sys, ALEA_CACHE_ALL);

    /* Insert into cell hash map */
    cell_hashmap_put(&sys->cell_index, final_cell_id, idx);

    /* Notify module hooks */
    if (sys->on_cell_added)
        sys->on_cell_added(sys->cell_hook_userdata, (size_t)idx);

    return idx;
}

int alea_max_cell_id(const alea_system_t* sys) {
    if (!sys || alea_vec_count(&sys->cells) == 0) return 0;

    int max_id = 0;
    for (size_t i = 0; i < alea_vec_count(&sys->cells); i++) {
        if (sys->cells.data[i].mc_cell_id > max_id) {
            max_id = sys->cells.data[i].mc_cell_id;
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

alea_material_id_t alea_material_at_point(alea_system_t* sys, double x, double y, double z) {
    if (!sys) return ALEA_MATERIAL_NONE;

    alea_cell_hit_t hit;
    if (alea_find_deepest_cell_hit_at_point(sys, x, y, z, &hit) != 0) {
        return ALEA_MATERIAL_NONE;
    }

    return hit.material_id;
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
    dst->mc_material_id = mixture->mc_material_id;
    dst->is_weight_fraction = mixture->is_weight_fraction;

    /* Copy components array via vec */
    alea_vec_init(&dst->components);
    size_t comp_count = alea_vec_count(&mixture->components);
    if (comp_count > 0) {
        int r = alea_vec_reserve(&dst->components, comp_count, alea_mixture_comp_t);
        if (r != 0) {
            alea_set_error_detail(ALEA_ERR_OUT_OF_MEMORY, "alea_add_mixture: failed to allocate mixture components");
            alea_vec_pop_discard(&sys->mixtures);
            return -1;
        }
        memcpy(dst->components.data, mixture->components.data,
               comp_count * sizeof(alea_mixture_comp_t));
        dst->components.count = comp_count;
    }

    /* Copy strings if present */
    dst->name = mixture->name ? alea_strdup(mixture->name) : NULL;
    dst->comments = mixture->comments ? alea_strdup(mixture->comments) : NULL;

    return (int)idx;
}

/* Comparison function for qsort */
static int compare_cell_ids(const void* a, const void* b) {
    int ia = *(const int*)a, ib = *(const int*)b;
    return (ia > ib) - (ia < ib);
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
        ids[i] = sys->cells.data[i].mc_cell_id;
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
        alea_vec_free(&sys->universes.data[i].cell_indices);
        alea_vec_free(&sys->universes.data[i].point_bvh_nodes);
        free(sys->universes.data[i].point_bvh_indices);
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
    atomic_fetch_or(&sys->query_cache_state, ALEA_CACHE_UNIVERSE);

    ALEA_LOG_INFO("Built universe index: %zu universes", alea_vec_count(&sys->universes));
    for (size_t i = 0; i < alea_vec_count(&sys->universes); i++) {
        ALEA_LOG_INFO("  Universe %d: %zu cells",
               sys->universes.data[i].universe_id, sys->universes.data[i].cell_indices.count);
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
    
    size_t count = (u->cell_indices.count < max_cells) ? u->cell_indices.count : max_cells;
    for (size_t i = 0; i < count; i++) {
        out_cells[i] = &sys->cells.data[u->cell_indices.data[i]];
    }
    return (int)count;
}

int alea_identify_cell_at_point(alea_system_t* sys, double x, double y, double z) {
    if (!sys) return -1;

    alea_cell_hit_t hit;
    if (alea_find_deepest_cell_hit_at_point(sys, x, y, z, &hit) != 0) {
        return -1;
    }
    return hit.cell_index;
}

int alea_find_overlaps(alea_system_t* sys, int* out_pairs, size_t max_pairs) {
    if (!sys || !out_pairs || max_pairs == 0) return 0;
    
    size_t found = 0;
    
    // Only check base universe cells
    for (size_t i = 0; i < alea_vec_count(&sys->cells) && found < max_pairs; i++) {
        ALEA_CHECK_INTERRUPTED(-1);
        const alea_cell_entry_t* cell_i = &sys->cells.data[i];
        if (cell_i->universe_id != 0) continue;
        if (cell_i->root_node_id >= alea_vec_count(&sys->nodes)) continue;

        const alea_bbox_t bbox_i_v = alea_node_bbox_get(&sys->nodes.data[cell_i->root_node_id].bbox);
        const alea_bbox_t* bbox_i = &bbox_i_v;

        for (size_t j = i + 1; j < alea_vec_count(&sys->cells) && found < max_pairs; j++) {
            const alea_cell_entry_t* cell_j = &sys->cells.data[j];
            if (cell_j->universe_id != 0) continue;
            if (cell_j->root_node_id >= alea_vec_count(&sys->nodes)) continue;

            const alea_bbox_t bbox_j_v = alea_node_bbox_get(&sys->nodes.data[cell_j->root_node_id].bbox);
            const alea_bbox_t* bbox_j = &bbox_j_v;

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
static int ensure_mc_id_to_surface_map(alea_system_t* sys) {
    size_t surf_count = alea_vec_count(&sys->surfaces);

    if (surf_count == 0) {
        return 0;
    }

    /* Find max MCNP surface ID to size the array */
    int max_mc_id = 0;
    for (size_t i = 0; i < surf_count; i++) {
        int id = sys->surfaces.data[i].mc_surface_id;
        if (id > max_mc_id) {
            max_mc_id = id;
        }
    }

    if (max_mc_id <= 0) {
        return 0;  /* No valid MCNP surface IDs */
    }

    size_t needed_size = (size_t)(max_mc_id + 1);

    /* Check if map needs rebuilding */
    if (sys->mc_id_to_surface && sys->mc_id_to_surface_size >= needed_size) {
        return 0;  /* Already built and big enough */
    }

    /* Free old map */
    free(sys->mc_id_to_surface);
    sys->mc_id_to_surface = NULL;
    sys->mc_id_to_surface_size = 0;

    /* Allocate new map */
    sys->mc_id_to_surface = malloc(needed_size * sizeof(uint32_t));
    if (!sys->mc_id_to_surface) {
        return -1;
    }
    sys->mc_id_to_surface_size = needed_size;

    /* Initialize to "no surface" */
    for (size_t i = 0; i < needed_size; i++) {
        sys->mc_id_to_surface[i] = UINT32_MAX;
    }

    /* Build mapping from surfaces */
    for (size_t i = 0; i < surf_count; i++) {
        int mcnp_id = sys->surfaces.data[i].mc_surface_id;
        if (mcnp_id > 0 && (size_t)mcnp_id < needed_size) {
            sys->mc_id_to_surface[mcnp_id] = (uint32_t)i;
        }
    }

    return 0;
}

/**
 * Helper: find surface index by MCNP surface ID
 * Returns UINT32_MAX if not found
 * Uses O(1) lookup table if available, falls back to O(n) search otherwise.
 */
static uint32_t find_surface_by_mc_id(const alea_system_t* sys, int mc_surface_id) {
    if (mc_surface_id <= 0) {
        return UINT32_MAX;
    }

    /* Use lookup table if available */
    if (sys->mc_id_to_surface && (size_t)mc_surface_id < sys->mc_id_to_surface_size) {
        return sys->mc_id_to_surface[mc_surface_id];
    }

    /* Fallback: linear search (should not happen if map is built) */
    for (size_t i = 0; i < alea_vec_count(&sys->surfaces); i++) {
        if (sys->surfaces.data[i].mc_surface_id == mc_surface_id) {
            return (uint32_t)i;
        }
    }
    return UINT32_MAX;
}

/**
 * Helper: recursively collect surface indices from a CSG tree
 *
 * Uses node->primitive.mc_surface_id to find the correct surface,
 * which handles cases where multiple surfaces share the same primitive
 * due to primitive deduplication. Falls back to prim_to_surface if
 * mc_surface_id is not set.
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

        /* First try: use node's mc_surface_id for correct surface lookup
         * This handles primitive deduplication correctly */
        int mcnp_id = node->primitive.mc_surface_id;
        if (mcnp_id > 0) {
            surf_idx = find_surface_by_mc_id(sys, mcnp_id);
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
    if (ensure_mc_id_to_surface_map(sys) != 0) {
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

    atomic_fetch_or(&sys->query_cache_state, ALEA_CACHE_CELL_SURFACES);
    return 0;
}

/* ============================================================================
 * CELL ADJACENCY
 * ============================================================================ */

/**
 * @brief Entry in the surface-to-cells map
 */
/**
 * @brief Helper: collect primitive senses from a cell's CSG tree
 */
typedef struct {
    uint32_t surface_idx;
    int8_t sense;
} cell_surface_sense_t;

static int collect_cell_surface_senses(const alea_system_t* sys,
                                       alea_node_id_t node_id,
                                       int current_sense,
                                       cell_surface_sense_t** out_senses,
                                       size_t* count,
                                       size_t* capacity) {
    if (node_id == ALEA_NODE_ID_INVALID || node_id >= alea_vec_count(&sys->nodes))
        return 0;

    const alea_node_t* node = &sys->nodes.data[node_id];
    alea_operation_t op = ALEA_GET_OPERATION(node);

    if (op == ALEA_OP_PRIMITIVE) {
        alea_primitive_id_t prim_id = node->primitive.primitive_id;
        int8_t node_sense = node->primitive.sense;
        int8_t inverted = node->primitive.inverted;
        /* Combine: node sense * dedup inversion * complement context */
        int8_t final_sense = (int8_t)(current_sense * node_sense * (inverted ? -1 : 1));

        uint32_t surf_idx = UINT32_MAX;
        int mc_surface_id = node->primitive.mc_surface_id;
        if (mc_surface_id > 0)
            surf_idx = find_surface_by_mc_id(sys, mc_surface_id);
        if (surf_idx == UINT32_MAX && prim_id < sys->prim_to_surface_size)
            surf_idx = sys->prim_to_surface[prim_id];

        if (surf_idx != UINT32_MAX) {
            /* Preserve both senses if a card appears in both branches, while
             * suppressing repeated uses of the same oriented halfspace. */
            for (size_t i = 0; i < *count; i++) {
                if ((*out_senses)[i].surface_idx == surf_idx &&
                    (*out_senses)[i].sense == final_sense) {
                    return 0;
                }
            }
            if (*count >= *capacity) {
                size_t new_cap = *capacity ? *capacity * 2 : 16;
                cell_surface_sense_t* new_arr = realloc(*out_senses,
                    new_cap * sizeof(cell_surface_sense_t));
                if (!new_arr) return -1;
                *out_senses = new_arr;
                *capacity = new_cap;
            }
            (*out_senses)[*count].surface_idx = surf_idx;
            (*out_senses)[*count].sense = final_sense;
            (*count)++;
        }
    } else if (op == ALEA_OP_COMPLEMENT) {
        if (collect_cell_surface_senses(sys, node->operation.left, -current_sense,
                                        out_senses, count, capacity) != 0)
            return -1;
    } else {
        if (collect_cell_surface_senses(sys, node->operation.left, current_sense,
                                        out_senses, count, capacity) != 0 ||
            collect_cell_surface_senses(sys, node->operation.right, current_sense,
                                        out_senses, count, capacity) != 0)
            return -1;
    }
    return 0;
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
    if (sys->neighbor_pool) {
        free(sys->neighbor_pool);
        sys->neighbor_pool = NULL;
    }
    free(sys->surface_cell_offsets);
    sys->surface_cell_offsets = NULL;
    free(sys->surface_cell_refs);
    sys->surface_cell_refs = NULL;
    sys->surface_cell_ref_count = 0;
    for (size_t i = 0; i < alea_vec_count(&sys->cells); i++) {
        sys->cells.data[i].neighbors = NULL;
        sys->cells.data[i].neighbor_count = 0;
    }
    sys->cell_adjacency_built = false;

    size_t num_surfaces = alea_vec_count(&sys->surfaces);
    size_t num_cells = alea_vec_count(&sys->cells);

    /* ------------------------------------------------------------------ *
     * Phase 1: Build flat surface→cell map (two-pass to avoid per-surface
     * reallocs that fragment the heap)
     * ------------------------------------------------------------------ */

    /* Reusable senses buffer across all cells */
    cell_surface_sense_t* senses_buf = NULL;
    size_t senses_cap = 0;

    /* Pass 1a: Count how many cell refs each surface will have */
    size_t* surf_counts = calloc(num_surfaces, sizeof(size_t));
    if (!surf_counts) { free(senses_buf); return -1; }

    for (size_t ci = 0; ci < num_cells; ci++) {
        alea_cell_entry_t* cell = &sys->cells.data[ci];
        if (cell->root_node_id == ALEA_NODE_ID_INVALID) continue;

        size_t sense_count = 0;
        if (collect_cell_surface_senses(sys, cell->root_node_id, 1,
                                        &senses_buf, &sense_count,
                                        &senses_cap) != 0) {
            free(surf_counts); free(senses_buf); return -1;
        }

        for (size_t si = 0; si < sense_count; si++) {
            uint32_t surf_idx = senses_buf[si].surface_idx;
            if (surf_idx < num_surfaces)
                surf_counts[surf_idx]++;
        }
    }

    /* Compute prefix-sum offsets and single allocation for all refs */
    size_t* surf_offsets = malloc((num_surfaces + 1) * sizeof(size_t));
    if (!surf_offsets) { free(surf_counts); free(senses_buf); return -1; }

    size_t total_refs = 0;
    for (size_t i = 0; i < num_surfaces; i++) {
        surf_offsets[i] = total_refs;
        total_refs += surf_counts[i];
    }
    surf_offsets[num_surfaces] = total_refs;

    alea_surface_cell_ref_t* all_refs = NULL;
    if (total_refs > 0) {
        all_refs = malloc(total_refs * sizeof(alea_surface_cell_ref_t));
        if (!all_refs) {
            free(surf_offsets); free(surf_counts); free(senses_buf);
            return -1;
        }
    }

    /* Reset counts for the filling pass */
    memset(surf_counts, 0, num_surfaces * sizeof(size_t));

    /* Pass 1b: Fill the flat refs array */
    for (size_t ci = 0; ci < num_cells; ci++) {
        alea_cell_entry_t* cell = &sys->cells.data[ci];
        if (cell->root_node_id == ALEA_NODE_ID_INVALID) continue;

        size_t sense_count = 0;
        if (collect_cell_surface_senses(sys, cell->root_node_id, 1,
                                        &senses_buf, &sense_count,
                                        &senses_cap) != 0)
            goto cleanup;

        for (size_t si = 0; si < sense_count; si++) {
            uint32_t surf_idx = senses_buf[si].surface_idx;
            if (surf_idx >= num_surfaces) continue;

            size_t idx = surf_offsets[surf_idx] + surf_counts[surf_idx];
            all_refs[idx].cell_index = (uint32_t)ci;
            all_refs[idx].sense = senses_buf[si].sense;
            surf_counts[surf_idx]++;
        }
    }

    free(senses_buf);
    senses_buf = NULL;

    /* ------------------------------------------------------------------ *
     * Phase 2: Build neighbor lists — two passes, zero pair materialization.
     *
     * The naive approach (enumerate all cell pairs per surface) creates
     * O(n²) pairs for surfaces shared by many cells (up to 2638 cells on
     * one surface → millions of pairs).  Instead we iterate per-cell:
     *
     * For each cell, walk its surfaces in the flat surf_map.  For each
     * surface, scan cells on the opposite sense in the same universe.
     * A reusable bitset (16 KB for 116K cells) deduplicates neighbors
     * within a single cell in O(1).
     *
     * Pass 1: count unique neighbors per cell.
     * Pass 2: allocate one pool, fill directly.
     * ------------------------------------------------------------------ */

    /* Build cell→surface-with-sense index so we can iterate per cell.
     * Reuse the senses_buf approach: for each cell, collect_cell_surface_senses
     * gives us the list of (surface_idx, sense).  We already did this twice
     * in Phase 1.  To avoid a third traversal, build a flat cell→surface map
     * from the existing surf_map by inverting it. */

    /* cell_surf_counts[ci] = number of surface refs for cell ci */
    size_t* cell_surf_counts = calloc(num_cells, sizeof(size_t));
    if (!cell_surf_counts) goto cleanup;

    for (size_t r = 0; r < total_refs; r++)
        cell_surf_counts[all_refs[r].cell_index]++;

    /* Prefix-sum for offsets */
    size_t* cell_surf_offsets = malloc(num_cells * sizeof(size_t));
    if (!cell_surf_offsets) { free(cell_surf_counts); goto cleanup; }
    {
        size_t off = 0;
        for (size_t i = 0; i < num_cells; i++) {
            cell_surf_offsets[i] = off;
            off += cell_surf_counts[i];
        }
    }

    /* Invert: for each surface ref, store (surface_idx, sense) keyed by cell */
    typedef struct { uint32_t surf_idx; int8_t sense; } cell_surf_entry_t;
    cell_surf_entry_t* cell_surf_refs = malloc(total_refs * sizeof(cell_surf_entry_t));
    if (!cell_surf_refs) { free(cell_surf_offsets); free(cell_surf_counts); goto cleanup; }

    /* Reset counts for filling */
    memset(cell_surf_counts, 0, num_cells * sizeof(size_t));
    for (size_t si = 0; si < num_surfaces; si++) {
        alea_surface_cell_ref_t* refs = &all_refs[surf_offsets[si]];
        size_t count = surf_counts[si];
        for (size_t r = 0; r < count; r++) {
            uint32_t ci = refs[r].cell_index;
            size_t idx = cell_surf_offsets[ci] + cell_surf_counts[ci];
            cell_surf_refs[idx].surf_idx = (uint32_t)si;
            cell_surf_refs[idx].sense = refs[r].sense;
            cell_surf_counts[ci]++;
        }
    }

    /* Surfaces shared by very many cells (e.g. a plane cutting through
     * the entire reactor) create O(n²) neighbor cliques that are useless
     * for local adjacency walking and dominate memory.  Skip them. */
    #define ADJACENCY_MAX_CELLS_PER_SURFACE 128

    /* Bitset for dedup: one bit per cell index */
    size_t bitset_words = (num_cells + 63) / 64;
    uint64_t* seen = calloc(bitset_words, sizeof(uint64_t));
    if (!seen) { free(cell_surf_refs); free(cell_surf_offsets); free(cell_surf_counts); goto cleanup; }

    /* Track which words we touched so we can reset in O(neighbors) not O(num_cells) */
    size_t* dirty_words = malloc(bitset_words * sizeof(size_t));
    if (!dirty_words) { free(seen); free(cell_surf_refs); free(cell_surf_offsets); free(cell_surf_counts); goto cleanup; }

    /* ---- Pass 1: Count unique neighbors per cell ---- */
    size_t* neighbor_count_arr = calloc(num_cells, sizeof(size_t));
    if (!neighbor_count_arr) { free(dirty_words); free(seen); free(cell_surf_refs); free(cell_surf_offsets); free(cell_surf_counts); goto cleanup; }

    for (size_t ci = 0; ci < num_cells; ci++) {
        int my_universe = sys->cells.data[ci].universe_id;
        size_t n_dirty = 0;
        size_t cnt = 0;

        /* Walk this cell's surfaces */
        size_t cs_off = cell_surf_offsets[ci];
        size_t cs_cnt = cell_surf_counts[ci];
        for (size_t s = 0; s < cs_cnt; s++) {
            uint32_t si = cell_surf_refs[cs_off + s].surf_idx;
            int8_t my_sense = cell_surf_refs[cs_off + s].sense;

            /* Skip surfaces shared by too many cells (global planes) */
            size_t ref_cnt = surf_counts[si];
            if (ref_cnt > ADJACENCY_MAX_CELLS_PER_SURFACE) continue;

            /* Walk cells on this surface with opposite sense */
            alea_surface_cell_ref_t* refs = &all_refs[surf_offsets[si]];
            for (size_t r = 0; r < ref_cnt; r++) {
                if (refs[r].sense == my_sense) continue;
                uint32_t other = refs[r].cell_index;
                if (other == (uint32_t)ci) continue;
                if (sys->cells.data[other].universe_id != my_universe) continue;

                /* Bitset dedup */
                size_t word = other / 64;
                uint64_t bit = (uint64_t)1 << (other % 64);
                if (seen[word] & bit) continue;
                if (seen[word] == 0) dirty_words[n_dirty++] = word;
                seen[word] |= bit;
                cnt++;
            }
        }

        neighbor_count_arr[ci] = cnt;

        /* Reset only touched bitset words */
        for (size_t d = 0; d < n_dirty; d++)
            seen[dirty_words[d]] = 0;
    }

    /* Compute pool size and offsets */
    size_t total_neighbors = 0;
    size_t* neighbor_offsets = malloc(num_cells * sizeof(size_t));
    if (!neighbor_offsets) { free(neighbor_count_arr); free(dirty_words); free(seen); free(cell_surf_refs); free(cell_surf_offsets); free(cell_surf_counts); goto cleanup; }

    for (size_t i = 0; i < num_cells; i++) {
        neighbor_offsets[i] = total_neighbors;
        total_neighbors += neighbor_count_arr[i];
    }

    ALEA_LOG_INFO("Cell adjacency: %zu unique neighbors, pool %.1f MB",
                  total_neighbors,
                  total_neighbors * sizeof(alea_cell_neighbor_t) / (1024.0*1024.0));

    /* Allocate pool */
    alea_cell_neighbor_t* pool = NULL;
    if (total_neighbors > 0) {
        pool = malloc(total_neighbors * sizeof(alea_cell_neighbor_t));
        if (!pool) { free(neighbor_offsets); free(neighbor_count_arr); free(dirty_words); free(seen); free(cell_surf_refs); free(cell_surf_offsets); free(cell_surf_counts); goto cleanup; }
    }

    /* Point each cell into its slice */
    for (size_t i = 0; i < num_cells; i++) {
        alea_cell_entry_t* c = &sys->cells.data[i];
        c->neighbors = (neighbor_count_arr[i] > 0) ? pool + neighbor_offsets[i] : NULL;
        c->neighbor_count = 0;
    }

    free(neighbor_count_arr);
    neighbor_count_arr = NULL;

    /* ---- Pass 2: Fill pool ---- */
    for (size_t ci = 0; ci < num_cells; ci++) {
        int my_universe = sys->cells.data[ci].universe_id;
        alea_cell_entry_t* ca = &sys->cells.data[ci];
        size_t n_dirty = 0;

        size_t cs_off = cell_surf_offsets[ci];
        size_t cs_cnt = cell_surf_counts[ci];
        for (size_t s = 0; s < cs_cnt; s++) {
            uint32_t si = cell_surf_refs[cs_off + s].surf_idx;
            int8_t my_sense = cell_surf_refs[cs_off + s].sense;

            size_t ref_cnt = surf_counts[si];
            if (ref_cnt > ADJACENCY_MAX_CELLS_PER_SURFACE) continue;

            alea_surface_cell_ref_t* refs = &all_refs[surf_offsets[si]];
            for (size_t r = 0; r < ref_cnt; r++) {
                if (refs[r].sense == my_sense) continue;
                uint32_t other = refs[r].cell_index;
                if (other == (uint32_t)ci) continue;
                if (sys->cells.data[other].universe_id != my_universe) continue;

                size_t word = other / 64;
                uint64_t bit = (uint64_t)1 << (other % 64);
                if (seen[word] & bit) continue;
                if (seen[word] == 0) dirty_words[n_dirty++] = word;
                seen[word] |= bit;

                size_t n = ca->neighbor_count;
                ca->neighbors[n].surface_id =
                    sys->surfaces.data[si].mc_surface_id;
                ca->neighbors[n].surface_index = si;
                ca->neighbors[n].neighbor_cell_id =
                    sys->cells.data[other].mc_cell_id;
                ca->neighbors[n].neighbor_index = other;
                ca->neighbors[n].our_sense = my_sense;
                ca->neighbor_count++;
            }
        }

        for (size_t d = 0; d < n_dirty; d++)
            seen[dirty_words[d]] = 0;
    }

    sys->neighbor_pool = pool;
    sys->surface_cell_offsets = surf_offsets;
    sys->surface_cell_refs = all_refs;
    sys->surface_cell_ref_count = total_refs;
    surf_offsets = NULL;
    all_refs = NULL;

    free(neighbor_offsets);
    free(dirty_words);
    free(seen);
    free(cell_surf_refs);
    free(cell_surf_offsets);
    free(cell_surf_counts);

    /* surf_offsets/all_refs are retained as the exact reverse CSR. */
    free(surf_counts);   surf_counts = NULL;

    sys->cell_adjacency_built = true;
    atomic_fetch_or(&sys->query_cache_state, ALEA_CACHE_ADJACENCY);

cleanup:
    free(all_refs);
    free(surf_offsets);
    free(surf_counts);
    free(senses_buf);

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
