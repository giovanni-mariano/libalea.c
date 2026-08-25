// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/**
 * @file alea_public_api.c
 * @brief Alea Public API Implementation
 *
 * Implementation of the public alea_* API functions.
 * Calls internal functions directly.
 */

#include "alea.h"
#include "core/alea_system.h"
#include "core/alea_materials.h"
#include "core/alea_ops.h"
#include "core/alea_eval.h"
#include "core/alea_void.h"
#include "core/alea_universe.h"
#include "core/alea_spatial_hier.h"
#include "core/alea_simplify.h"
#include "primitives/primitive_create.h"
#include "primitives/bbox.h"
#include "util/alea_log.h"

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "util/arena.h"
#include "util/str_builder.h"
#include "util/compat.h"

#ifdef _OPENMP
#include <omp.h>
#endif

#define STRINGIFY_(x) #x
#define STRINGIFY(x) STRINGIFY_(x)
#define ALEA_VERSION_STRING \
    STRINGIFY(ALEA_VERSION_MAJOR) "." \
    STRINGIFY(ALEA_VERSION_MINOR) "." \
    STRINGIFY(ALEA_VERSION_PATCH)

/* ============================================================================
 * VERSION
 * ============================================================================ */

const char* alea_version(void) {
    return ALEA_VERSION_STRING;
}

int alea_openmp_enabled(void) {
#ifdef _OPENMP
    return 1;
#else
    return 0;
#endif
}

int alea_openmp_max_threads(void) {
#ifdef _OPENMP
    return omp_get_max_threads();
#else
    return 1;
#endif
}

/* ============================================================================
 * ERROR HANDLING
 * ============================================================================ */

const char* alea_error(void) {
    return alea_get_error_detail();
}

int alea_error_code(void) {
    return (int)alea_get_last_error();
}

void alea_error_clear(void) {
    alea_clear_error_detail();
}

/* ============================================================================
 * INTERRUPT SUPPORT
 * ============================================================================ */

void alea_interrupt(void) {
    g_alea_interrupted = 1;
}

void alea_clear_interrupt(void) {
    g_alea_interrupted = 0;
}

bool alea_interrupted(void) {
    return g_alea_interrupted != 0;
}

/* ============================================================================
 * SYSTEM LIFECYCLE
 * ============================================================================ */

alea_system_t* alea_create(void) {
    alea_system_t* sys = alea_system_create();
    if (!sys) {
        alea_set_error_detail(ALEA_ERR_OUT_OF_MEMORY, "Failed to allocate CSG system");
    }
    return sys;
}

void alea_destroy(alea_system_t* sys) {
    if (sys) {
        alea_system_destroy(sys);
    }
}

/* Helper macro to clone a vector (shallow copy — safe only for POD types) */
#define CLONE_VEC(dst, src, type) do { \
    size_t count = alea_vec_count(&(src)); \
    if (count > 0) { \
        int _r = alea_vec_reserve(&(dst), count, type); \
        if (_r != 0) goto clone_error; \
        memcpy((dst).data, (src).data, count * sizeof(type)); \
        (dst).count = count; \
    } \
} while(0)

alea_system_t* alea_clone(const alea_system_t* sys) {
    if (!sys) return NULL;

    alea_system_t* clone = alea_system_create();
    if (!clone) return NULL;

    /* Clone vectors without internal pointers (safe shallow copy) */
    CLONE_VEC(clone->nodes, sys->nodes, alea_node_t);
    CLONE_VEC(clone->primitives, sys->primitives, alea_primitive_entry_t);
    CLONE_VEC(clone->primitive_planes, sys->primitive_planes, alea_plane_data_t);
    CLONE_VEC(clone->primitive_spheres, sys->primitive_spheres, alea_sphere_data_t);
    CLONE_VEC(clone->primitive_cyl_x, sys->primitive_cyl_x, alea_cylinder_x_data_t);
    CLONE_VEC(clone->primitive_cyl_y, sys->primitive_cyl_y, alea_cylinder_y_data_t);
    CLONE_VEC(clone->primitive_cyl_z, sys->primitive_cyl_z, alea_cylinder_z_data_t);
    CLONE_VEC(clone->primitive_cone_x, sys->primitive_cone_x, alea_cone_x_data_t);
    CLONE_VEC(clone->primitive_cone_y, sys->primitive_cone_y, alea_cone_y_data_t);
    CLONE_VEC(clone->primitive_cone_z, sys->primitive_cone_z, alea_cone_z_data_t);
    CLONE_VEC(clone->primitive_boxes, sys->primitive_boxes, alea_box_data_t);
    CLONE_VEC(clone->primitive_quadrics, sys->primitive_quadrics, alea_quadric_data_t);
    CLONE_VEC(clone->primitive_toruses, sys->primitive_toruses, alea_torus_data_t);
    CLONE_VEC(clone->primitive_rccs, sys->primitive_rccs, alea_rcc_data_t);
    CLONE_VEC(clone->primitive_box_generals, sys->primitive_box_generals, alea_box_general_data_t);
    CLONE_VEC(clone->primitive_sphs, sys->primitive_sphs, alea_sph_data_t);
    CLONE_VEC(clone->primitive_trcs, sys->primitive_trcs, alea_trc_data_t);
    CLONE_VEC(clone->primitive_ells, sys->primitive_ells, alea_ell_data_t);
    CLONE_VEC(clone->primitive_recs, sys->primitive_recs, alea_rec_data_t);
    CLONE_VEC(clone->primitive_weds, sys->primitive_weds, alea_wed_data_t);
    CLONE_VEC(clone->primitive_rhps, sys->primitive_rhps, alea_rhp_data_t);
    CLONE_VEC(clone->primitive_arbs, sys->primitive_arbs, alea_arb_data_t);
    CLONE_VEC(clone->surfaces, sys->surfaces, alea_surface_entry_t);
    CLONE_VEC(clone->transforms, sys->transforms, alea_transform_t);
    CLONE_VEC(clone->cell_refs, sys->cell_refs, alea_cell_ref_t);

    /* Clone cells: shallow copy then fix up internal pointers */
    CLONE_VEC(clone->cells, sys->cells, alea_cell_entry_t);
    for (size_t i = 0; i < alea_vec_count(&clone->cells); i++) {
        alea_cell_entry_t* c = &clone->cells.data[i];
        const alea_cell_entry_t* s = &sys->cells.data[i];
        /* Index data — rebuilt on demand, NULL out to prevent double-free */
        c->surface_indices = NULL;
        c->surface_index_count = 0;
        c->neighbors = NULL;
        c->neighbor_count = 0;
        /* Deep-copy owned lattice fill array */
        if (s->lat_fill && s->lat_fill_count > 0) {
            c->lat_fill = malloc(s->lat_fill_count * sizeof(int));
            if (!c->lat_fill) goto clone_error;
            memcpy(c->lat_fill, s->lat_fill, s->lat_fill_count * sizeof(int));
        } else {
            c->lat_fill = NULL;
        }
        /* Deep-copy comment strings */
        c->comments = s->comments ? alea_strdup(s->comments) : NULL;
        c->inline_comment = s->inline_comment ? alea_strdup(s->inline_comment) : NULL;
    }

    /* Clone materials: shallow copy then deep-copy internal arrays */
    CLONE_VEC(clone->materials, sys->materials, alea_material_t);
    for (size_t i = 0; i < alea_vec_count(&clone->materials); i++) {
        alea_material_t* m = &clone->materials.data[i];
        const alea_material_t* s = &sys->materials.data[i];
        /* NULL everything first so destroy is safe on partial failure */
        alea_vec_init(&m->nuclides);
        alea_vec_init(&m->elements);
        alea_vec_init(&m->thermal_laws);
        m->name = NULL;
        m->comments = NULL;
        /* Deep-copy nuclides */
        if (alea_vec_count(&s->nuclides) > 0) {
            int r = alea_vec_reserve(&m->nuclides, alea_vec_count(&s->nuclides), alea_nuclide_t);
            if (r != 0) goto clone_error;
            memcpy(m->nuclides.data, s->nuclides.data, alea_vec_count(&s->nuclides) * sizeof(alea_nuclide_t));
            m->nuclides.count = alea_vec_count(&s->nuclides);
            for (size_t j = 0; j < alea_vec_count(&s->nuclides); j++) {
                m->nuclides.data[j].library = alea_strdup(s->nuclides.data[j].library);
            }
        }
        /* Deep-copy elements */
        if (alea_vec_count(&s->elements) > 0) {
            int r = alea_vec_reserve(&m->elements, alea_vec_count(&s->elements), alea_element_comp_t);
            if (r != 0) goto clone_error;
            memcpy(m->elements.data, s->elements.data, alea_vec_count(&s->elements) * sizeof(alea_element_comp_t));
            m->elements.count = alea_vec_count(&s->elements);
            for (size_t j = 0; j < alea_vec_count(&s->elements); j++) {
                m->elements.data[j].library = alea_strdup(s->elements.data[j].library);
            }
        }
        /* Deep-copy thermal laws */
        if (alea_vec_count(&s->thermal_laws) > 0) {
            int r = alea_vec_reserve(&m->thermal_laws, alea_vec_count(&s->thermal_laws), alea_thermal_law_t);
            if (r != 0) goto clone_error;
            memcpy(m->thermal_laws.data, s->thermal_laws.data, alea_vec_count(&s->thermal_laws) * sizeof(alea_thermal_law_t));
            m->thermal_laws.count = alea_vec_count(&s->thermal_laws);
            for (size_t j = 0; j < alea_vec_count(&s->thermal_laws); j++) {
                m->thermal_laws.data[j].identifier = alea_strdup(s->thermal_laws.data[j].identifier);
            }
        }
        m->name = alea_strdup(s->name);
        m->comments = alea_strdup(s->comments);
    }

    /* Clone universes: shallow copy then NULL out index data */
    CLONE_VEC(clone->universes, sys->universes, alea_universe_t);
    for (size_t i = 0; i < alea_vec_count(&clone->universes); i++) {
        alea_vec_init(&clone->universes.data[i].cell_indices);
        alea_vec_init(&clone->universes.data[i].point_bvh_nodes);
        clone->universes.data[i].point_bvh_indices = NULL;
        clone->universes.data[i].point_bvh_built = false;
        clone->universes.data[i].point_bvh_disabled = false;
    }

    /* Clone mixtures: shallow copy then deep-copy internal arrays */
    CLONE_VEC(clone->mixtures, sys->mixtures, alea_mixture_t);
    for (size_t i = 0; i < alea_vec_count(&clone->mixtures); i++) {
        alea_mixture_t* m = &clone->mixtures.data[i];
        const alea_mixture_t* s = &sys->mixtures.data[i];
        alea_vec_init(&m->components);
        m->name = NULL;
        m->comments = NULL;
        if (alea_vec_count(&s->components) > 0) {
            int r = alea_vec_reserve(&m->components, alea_vec_count(&s->components), alea_mixture_comp_t);
            if (r != 0) goto clone_error;
            memcpy(m->components.data, s->components.data, alea_vec_count(&s->components) * sizeof(alea_mixture_comp_t));
            m->components.count = alea_vec_count(&s->components);
        }
        m->name = alea_strdup(s->name);
        m->comments = alea_strdup(s->comments);
    }

    /* Clone surface lookup if present */
    if (sys->surface_lookup && sys->surface_lookup_size > 0) {
        clone->surface_lookup = malloc(sys->surface_lookup_size * sizeof(alea_node_id_t));
        if (!clone->surface_lookup) goto clone_error;
        memcpy(clone->surface_lookup, sys->surface_lookup,
               sys->surface_lookup_size * sizeof(alea_node_id_t));
        clone->surface_lookup_size = sys->surface_lookup_size;
    }

    /* Clone prim_to_surface if present */
    if (sys->prim_to_surface && sys->prim_to_surface_size > 0) {
        clone->prim_to_surface = malloc(sys->prim_to_surface_size * sizeof(uint32_t));
        if (!clone->prim_to_surface) goto clone_error;
        memcpy(clone->prim_to_surface, sys->prim_to_surface,
               sys->prim_to_surface_size * sizeof(uint32_t));
        clone->prim_to_surface_size = sys->prim_to_surface_size;
    }

    /* Copy configuration */
    clone->config = sys->config;
    clone->source = sys->source;
    clone->next_inline_transform_id = sys->next_inline_transform_id;
    clone->next_auto_surface_id = sys->next_auto_surface_id;
    clone->stats = sys->stats;
    /* Universe index was invalidated (cell_indices NULLed), must be rebuilt */
    clone->universe_index_built = false;

    /* Rebuild cell hash map for cloned system */
    cell_hashmap_clear(&clone->cell_index);
    for (size_t i = 0; i < alea_vec_count(&clone->cells); i++) {
        cell_hashmap_put(&clone->cell_index, clone->cells.data[i].mc_cell_id, (int)i);
    }

    /* Rebuild the primitive dedup index. The primitives vector was copied above,
       but alea_system_create() gave the clone an empty index. Unlike the query
       caches below, the dedup index has no build-on-demand path, so it must be
       repopulated here or mutating the clone would miss dedup and skew stats. */
    for (size_t i = 0; i < alea_vec_count(&clone->primitives); i++) {
        const alea_primitive_entry_t* p = &clone->primitives.data[i];
        alea_primitive_data_t data;
        if (!alea_primitive_copy_data(clone, (uint32_t)i, &data)) goto clone_error;
        uint64_t hash = alea_compute_primitive_hash(p->type, &data, &clone->config);
        primitive_hash_table_insert(clone->primitive_index, (uint32_t)i, hash);
    }

    /* Note: instance_cache, surface_bvh, spatial_index, hier_spatial_index
       are NOT cloned - they are query caches rebuilt on demand */

    return clone;

clone_error:
    alea_system_destroy(clone);
    return NULL;
}

#undef CLONE_VEC

void alea_reset(alea_system_t* sys) {
    if (sys) {
        alea_system_reset(sys);
    }
}

/* ============================================================================
 * CONFIGURATION
 * ============================================================================ */

alea_config_t alea_get_config(const alea_system_t* sys) {
    if (!sys) return ALEA_CONFIG_DEFAULT;
    return sys->config;
}

void alea_set_config(alea_system_t* sys, const alea_config_t* config) {
    if (!sys || !config) return;
    sys->config = *config;
}

/* ============================================================================
 * LOADING
 *
 * MCNP loading: mcnp_load() in alea_mcnp.h
 * OpenMC loading: openmc_load() in alea_openmc.h
 * ============================================================================ */

/* ============================================================================
 * INDEXING
 * ============================================================================ */

int alea_query_acceleration_stats(const alea_system_t* sys,
                                  alea_query_acceleration_stats_t* out_stats) {
    if (!sys || !out_stats) return -1;

    memset(out_stats, 0, sizeof(*out_stats));

    const alea_hier_spatial_stats_t* stats =
        alea_hier_spatial_index_stats(sys->hier_spatial_index);
    out_stats->built = stats != NULL;
    if (!stats) return 0;

    out_stats->hier_universe_count = stats->universe_count;
    out_stats->hier_blas_count = stats->blas_count;
    out_stats->hier_linear_universe_count = stats->linear_universe_count;
    out_stats->hier_blas_cell_count = stats->blas_cell_count;
    out_stats->hier_blas_node_count = stats->blas_node_count;
    out_stats->hier_fill_cell_count = stats->fill_cell_count;
    out_stats->hier_lattice_cell_count = stats->lattice_cell_count;
    out_stats->hier_transform_count = stats->transform_count;
    out_stats->hier_placement_count = stats->placement_count;
    out_stats->hier_root_placement_count = stats->root_placement_count;
    out_stats->hier_fill_placement_count = stats->fill_placement_count;
    out_stats->hier_lattice_placement_count = stats->lattice_placement_count;
    out_stats->hier_max_placement_depth = stats->max_placement_depth;
    out_stats->hier_max_universe_cells = stats->max_universe_cells;
    out_stats->hier_largest_universe_id = stats->largest_universe_id;
    out_stats->memory_bytes = stats->memory_bytes;
    return 0;
}

int alea_surface_reference_stats(const alea_system_t* sys,
                                 alea_surface_reference_stats_t* out_stats) {
    if (!sys || !out_stats) return -1;

    memset(out_stats, 0, sizeof(*out_stats));
    out_stats->built = sys->cell_adjacency_built &&
                       sys->surface_cell_offsets != NULL;
    if (!out_stats->built) return 0;

    out_stats->surface_count = alea_vec_count(&sys->surfaces);
    out_stats->reference_count = sys->surface_cell_ref_count;
    out_stats->memory_bytes =
        (out_stats->surface_count + 1) * sizeof(size_t) +
        out_stats->reference_count * sizeof(alea_surface_cell_ref_t);
    for (size_t i = 0; i < out_stats->surface_count; i++) {
        size_t count = sys->surface_cell_offsets[i + 1] -
                       sys->surface_cell_offsets[i];
        if (count > out_stats->max_references_per_surface)
            out_stats->max_references_per_surface = count;
    }
    return 0;
}

int alea_surface_cell_references(const alea_system_t* sys,
                                 int surface_id,
                                 alea_surface_cell_reference_t* out_refs,
                                 size_t capacity,
                                 size_t* out_count) {
    if (!sys || !out_count || (capacity > 0 && !out_refs)) {
        alea_set_error_detail(ALEA_ERR_NULL_ARG,
                              "alea_surface_cell_references: invalid output");
        return -1;
    }
    *out_count = 0;
    if (!sys->cell_adjacency_built || !sys->surface_cell_offsets) {
        alea_set_error_detail(ALEA_ERR_INVALID_ARG,
                              "alea_surface_cell_references: adjacency cache is not prepared");
        return -1;
    }

    uint32_t surface_index = UINT32_MAX;
    if (surface_id > 0 && sys->mc_id_to_surface &&
        (size_t)surface_id < sys->mc_id_to_surface_size) {
        surface_index = sys->mc_id_to_surface[surface_id];
    }
    if (surface_index == UINT32_MAX) {
        alea_set_error_detail(ALEA_ERR_NOT_FOUND,
                              "alea_surface_cell_references: surface %d not found",
                              surface_id);
        return -1;
    }

    size_t begin = sys->surface_cell_offsets[surface_index];
    size_t end = sys->surface_cell_offsets[surface_index + 1];
    size_t count = end - begin;
    *out_count = count;
    size_t copied = capacity < count ? capacity : count;
    for (size_t i = 0; i < copied; i++) {
        const alea_surface_cell_ref_t* ref =
            &sys->surface_cell_refs[begin + i];
        const alea_cell_entry_t* cell = &sys->cells.data[ref->cell_index];
        out_refs[i].cell_index = ref->cell_index;
        out_refs[i].cell_id = cell->mc_cell_id;
        out_refs[i].universe_id = cell->universe_id;
        out_refs[i].sense = ref->sense;
    }
    return 0;
}

/* ============================================================================
 * GEOMETRY QUERIES
 * ============================================================================ */

int alea_find_cell(alea_system_t* sys, double x, double y, double z) {
    if (!sys) return -1;

    alea_cell_hit_t hit;
    if (alea_find_deepest_cell_hit_at_point(sys, x, y, z, &hit) != 0) {
        return -1;
    }
    return hit.cell_index;
}

int alea_find_all_cells(alea_system_t* sys, double x, double y, double z,
                            alea_cell_hit_t* hits, size_t max_hits) {
    if (!sys || !hits || max_hits == 0) return -1;
    /* Internal alea_cell_hit_t has same layout as public alea_cell_hit_t */
    return alea_find_all_cells_at_point(sys, x, y, z, (alea_cell_hit_t*)hits, max_hits);
}

int alea_find_all_cells_coverage_chain(alea_system_t* sys,
                                       double x, double y, double z,
                                       alea_cell_hit_t* hits,
                                       uint64_t* occurrence_keys,
                                       uint64_t* parent_occurrence_keys,
                                       size_t max_hits) {
    if (!sys || !hits || !occurrence_keys || !parent_occurrence_keys ||
        max_hits == 0) return -1;
    return alea_find_all_cells_at_point_coverage_chain_recursive(
        sys, x, y, z, hits, occurrence_keys, parent_occurrence_keys, max_hits);
}

int alea_find_all_cells_in_universe_coverage_chain(
    alea_system_t* sys, int universe_id,
    double x, double y, double z,
    alea_cell_hit_t* hits,
    uint64_t* occurrence_keys,
    uint64_t* parent_occurrence_keys,
    size_t max_hits) {
    if (!sys || !hits || !occurrence_keys || !parent_occurrence_keys ||
        max_hits == 0) return -1;
    return alea_find_all_cells_in_universe_at_point_coverage_chain_recursive(
        sys, universe_id, x, y, z, hits, occurrence_keys,
        parent_occurrence_keys, max_hits);
}

int alea_classify_point_coverage_chain(
    const alea_cell_hit_t* hits,
    const uint64_t* occurrence_keys,
    const uint64_t* parent_occurrence_keys,
    size_t hit_count,
    int universe_depth,
    uint8_t* out_owner_mask,
    alea_point_coverage_classification_t* out_classification) {
    if (!out_classification ||
        (hit_count > 0 && (!hits || !occurrence_keys || !parent_occurrence_keys)))
        return -1;

    *out_classification = (alea_point_coverage_classification_t){
        .kind = ALEA_POINT_COVERAGE_GAP,
        .target_depth = universe_depth,
        .owner_count = 0
    };
    if (out_owner_mask && hit_count > 0)
        memset(out_owner_mask, 0, hit_count * sizeof(*out_owner_mask));
    if (hit_count == 0) return 0;

    size_t* child_counts = calloc(hit_count, sizeof(*child_counts));
    if (!child_counts) return -1;
    size_t root_count = 0;
    int malformed = 0;
    for (size_t i = 0; i < hit_count && !malformed; i++) {
        if (occurrence_keys[i] == 0) {
            malformed = 1;
            break;
        }
        for (size_t j = i + 1; j < hit_count; j++) {
            if (occurrence_keys[i] == occurrence_keys[j]) {
                malformed = 1;
                break;
            }
        }
        if (malformed) break;
        if (parent_occurrence_keys[i] == 0) {
            root_count++;
            continue;
        }
        size_t parent = 0;
        while (parent < hit_count &&
               occurrence_keys[parent] != parent_occurrence_keys[i])
            parent++;
        if (parent == hit_count || hits[i].depth != hits[parent].depth + 1) {
            malformed = 1;
            break;
        }
        child_counts[parent]++;
    }
    if (malformed || root_count == 0) {
        out_classification->kind = ALEA_POINT_COVERAGE_UNRESOLVED;
        free(child_counts);
        return 0;
    }

    size_t owner_count = 0;
    int one_owner = -1;
    int common_depth = -1;
    int mixed_depths = 0;
    for (size_t i = 0; i < hit_count; i++) {
        int selected = universe_depth >= 0
            ? hits[i].depth == universe_depth
            : child_counts[i] == 0;
        if (!selected) continue;
        if (out_owner_mask) out_owner_mask[i] = 1;
        one_owner = (int)i;
        owner_count++;
        if (common_depth < 0) common_depth = hits[i].depth;
        else if (common_depth != hits[i].depth) mixed_depths = 1;
    }
    out_classification->owner_count = owner_count;
    out_classification->target_depth = universe_depth >= 0
        ? universe_depth : (mixed_depths ? -1 : common_depth);
    if (owner_count == 0) {
        out_classification->kind = ALEA_POINT_COVERAGE_GAP;
    } else if (owner_count > 1) {
        out_classification->kind = ALEA_POINT_COVERAGE_OVERLAP;
    } else if (hits[one_owner].resolution_flags & ALEA_RESOLVE_UNDEFINED_FILL) {
        out_classification->kind = ALEA_POINT_COVERAGE_UNDEFINED_FILL;
    } else {
        out_classification->kind = ALEA_POINT_COVERAGE_UNIQUE;
    }
    free(child_counts);
    return 0;
}

bool alea_point_inside(const alea_system_t* sys, alea_node_id_t node,
                           double x, double y, double z) {
    if (!sys) return false;
    return alea_contains_point(sys, node, x, y, z);
}

int alea_material_at(alea_system_t* sys, double x, double y, double z) {
    if (!sys) return -1;

    alea_cell_hit_t hit;
    if (alea_find_deepest_cell_hit_at_point(sys, x, y, z, &hit) != 0) {
        return 0;  /* void */
    }
    return hit.material_id;
}


/* ============================================================================
 * CSG CONSTRUCTION - SURFACES
 * ============================================================================ */

/**
 * @brief Helper to create a complete surface entry with both sense nodes
 *
 * Follows the MCNP conversion pattern:
 * 1. Get or create primitive (with deduplication)
 * 2. Create positive sense node
 * 3. Create negative sense node
 * 4. Register surface entry
 *
 * @return Index into sys->surfaces array, or -1 on error
 */
static int create_surface_entry(alea_system_t* sys,
                                int surface_id,
                                alea_primitive_type_t type,
                                alea_primitive_data_t* data) {
    if (!sys || !data) return -1;

    if (surface_id <= 0) {
        surface_id = sys->next_auto_surface_id++;
    } else if (surface_id >= sys->next_auto_surface_id) {
        sys->next_auto_surface_id = surface_id + 1;
    }

    /* Get or create primitive (with automatic deduplication) */
    int8_t inverted = 0;
    alea_primitive_id_t prim_id = alea_get_or_create_primitive(sys, type, data, &inverted);
    if (prim_id == UINT32_MAX) return -1;

    /* Create POSITIVE and NEGATIVE sense nodes */
    alea_node_id_t pos_node = alea_add_primitive_node(sys, prim_id, +1, inverted, surface_id);
    if (pos_node == ALEA_NODE_ID_INVALID) return -1;

    alea_node_id_t neg_node = alea_add_primitive_node(sys, prim_id, -1, inverted, surface_id);
    if (neg_node == ALEA_NODE_ID_INVALID) return -1;

    /* Register surface entry - get index before push (count - 1 after push) */
    size_t index = alea_vec_count(&sys->surfaces);
    alea_surface_entry_t* surf = alea_vec_push_uninit(&sys->surfaces, alea_surface_entry_t);
    if (!surf) return -1;

    *surf = (alea_surface_entry_t){
        .mc_surface_id = surface_id,
        .primitive_id = prim_id,
        .pos_node = pos_node,
        .neg_node = neg_node,
        .transform_id = 0,
        .transform_applied = false,
        .boundary_type = ALEA_BOUNDARY_TRANSMISSIVE,
        .periodic_surface_id = 0,
        .expanded_pos_node = ALEA_NODE_ID_INVALID,
        .expanded_neg_node = ALEA_NODE_ID_INVALID
    };

    alea_system_invalidate_query_caches(sys, ALEA_CACHE_ALL);
    return (int)index;
}

int alea_plane_surface(alea_system_t* sys, int surface_id,
                           double a, double b, double c, double d) {
    alea_primitive_data_t data = {0};
    data.plane.a = a;
    data.plane.b = b;
    data.plane.c = c;
    data.plane.d = d;
    return create_surface_entry(sys, surface_id, ALEA_PRIMITIVE_PLANE, &data);
}

int alea_sphere_surface(alea_system_t* sys, int surface_id,
                                              double cx, double cy, double cz, double r) {
    alea_primitive_data_t data = {0};
    data.sphere.center_x = cx;
    data.sphere.center_y = cy;
    data.sphere.center_z = cz;
    data.sphere.radius = r;
    return create_surface_entry(sys, surface_id, ALEA_PRIMITIVE_SPHERE, &data);
}

int alea_cylinder_z_surface(alea_system_t* sys, int surface_id,
                                                  double cx, double cy, double r) {
    alea_primitive_data_t data = {0};
    data.cyl_z.center_x = cx;
    data.cyl_z.center_y = cy;
    data.cyl_z.radius = r;
    return create_surface_entry(sys, surface_id, ALEA_PRIMITIVE_CYLINDER_Z, &data);
}

int alea_cylinder_x_surface(alea_system_t* sys, int surface_id,
                                                  double cy, double cz, double r) {
    alea_primitive_data_t data = {0};
    data.cyl_x.center_y = cy;
    data.cyl_x.center_z = cz;
    data.cyl_x.radius = r;
    return create_surface_entry(sys, surface_id, ALEA_PRIMITIVE_CYLINDER_X, &data);
}

int alea_cylinder_y_surface(alea_system_t* sys, int surface_id,
                                                  double cx, double cz, double r) {
    alea_primitive_data_t data = {0};
    data.cyl_y.center_x = cx;
    data.cyl_y.center_z = cz;
    data.cyl_y.radius = r;
    return create_surface_entry(sys, surface_id, ALEA_PRIMITIVE_CYLINDER_Y, &data);
}

int alea_box_surface(alea_system_t* sys, int surface_id,
                                           double xmin, double xmax,
                                           double ymin, double ymax,
                                           double zmin, double zmax) {
    alea_primitive_data_t data = {0};
    data.box.min_x = xmin;
    data.box.max_x = xmax;
    data.box.min_y = ymin;
    data.box.max_y = ymax;
    data.box.min_z = zmin;
    data.box.max_z = zmax;
    return create_surface_entry(sys, surface_id, ALEA_PRIMITIVE_RPP, &data);
}

int alea_cone_z_surface(alea_system_t* sys, int surface_id,
                                              double cx, double cy, double cz,
                                              double t_squared) {
    alea_primitive_data_t data = {0};
    data.cone_z.apex_x = cx;
    data.cone_z.apex_y = cy;
    data.cone_z.apex_z = cz;
    data.cone_z.tan_angle_sq = t_squared;
    return create_surface_entry(sys, surface_id, ALEA_PRIMITIVE_CONE_Z, &data);
}

int alea_cone_x_surface(alea_system_t* sys, int surface_id,
                                              double cx, double cy, double cz,
                                              double t_squared) {
    alea_primitive_data_t data = {0};
    data.cone_x.apex_x = cx;
    data.cone_x.apex_y = cy;
    data.cone_x.apex_z = cz;
    data.cone_x.tan_angle_sq = t_squared;
    return create_surface_entry(sys, surface_id, ALEA_PRIMITIVE_CONE_X, &data);
}

int alea_cone_y_surface(alea_system_t* sys, int surface_id,
                                              double cx, double cy, double cz,
                                              double t_squared) {
    alea_primitive_data_t data = {0};
    data.cone_y.apex_x = cx;
    data.cone_y.apex_y = cy;
    data.cone_y.apex_z = cz;
    data.cone_y.tan_angle_sq = t_squared;
    return create_surface_entry(sys, surface_id, ALEA_PRIMITIVE_CONE_Y, &data);
}

int alea_torus_z_surface(alea_system_t* sys, int surface_id,
                                               double cx, double cy, double cz,
                                               double major_radius, double minor_radius) {
    alea_primitive_data_t data = {0};
    data.torus.axis = ALEA_AXIS_Z;
    data.torus.center_x = cx;
    data.torus.center_y = cy;
    data.torus.center_z = cz;
    data.torus.major_radius = major_radius;
    data.torus.minor_radius = minor_radius;
    data.torus.axial_semiwidth_B = minor_radius;  /* Circular cross-section */
    return create_surface_entry(sys, surface_id, ALEA_PRIMITIVE_TORUS_Z, &data);
}

int alea_torus_x_surface(alea_system_t* sys, int surface_id,
                                               double cx, double cy, double cz,
                                               double major_radius, double minor_radius) {
    alea_primitive_data_t data = {0};
    data.torus.axis = ALEA_AXIS_X;
    data.torus.center_x = cx;
    data.torus.center_y = cy;
    data.torus.center_z = cz;
    data.torus.major_radius = major_radius;
    data.torus.minor_radius = minor_radius;
    data.torus.axial_semiwidth_B = minor_radius;
    return create_surface_entry(sys, surface_id, ALEA_PRIMITIVE_TORUS_X, &data);
}

int alea_torus_y_surface(alea_system_t* sys, int surface_id,
                                               double cx, double cy, double cz,
                                               double major_radius, double minor_radius) {
    alea_primitive_data_t data = {0};
    data.torus.axis = ALEA_AXIS_Y;
    data.torus.center_x = cx;
    data.torus.center_y = cy;
    data.torus.center_z = cz;
    data.torus.major_radius = major_radius;
    data.torus.minor_radius = minor_radius;
    data.torus.axial_semiwidth_B = minor_radius;
    return create_surface_entry(sys, surface_id, ALEA_PRIMITIVE_TORUS_Y, &data);
}

int alea_quadric_surface(alea_system_t* sys, int surface_id,
                                               double A, double B, double C,
                                               double D, double E, double F,
                                               double G, double H, double I, double J) {
    alea_primitive_data_t data = {0};
    data.quadric.coeffs[0] = A;
    data.quadric.coeffs[1] = B;
    data.quadric.coeffs[2] = C;
    data.quadric.coeffs[3] = D;
    data.quadric.coeffs[4] = E;
    data.quadric.coeffs[5] = F;
    data.quadric.coeffs[6] = G;
    data.quadric.coeffs[7] = H;
    data.quadric.coeffs[8] = I;
    data.quadric.coeffs[9] = J;
    return create_surface_entry(sys, surface_id, ALEA_PRIMITIVE_QUADRIC, &data);
}

/* --- Macrobody surfaces --- */

int alea_rcc_surface(alea_system_t* sys, int surface_id,
                                           double base_x, double base_y, double base_z,
                                           double height_x, double height_y, double height_z,
                                           double radius) {
    alea_primitive_data_t data = {0};
    data.rcc.base_x = base_x;
    data.rcc.base_y = base_y;
    data.rcc.base_z = base_z;
    data.rcc.height_x = height_x;
    data.rcc.height_y = height_y;
    data.rcc.height_z = height_z;
    data.rcc.radius = radius;
    return create_surface_entry(sys, surface_id, ALEA_PRIMITIVE_RCC, &data);
}

int alea_box_general_surface(alea_system_t* sys, int surface_id,
                                                   double corner_x, double corner_y, double corner_z,
                                                   double v1_x, double v1_y, double v1_z,
                                                   double v2_x, double v2_y, double v2_z,
                                                   double v3_x, double v3_y, double v3_z) {
    alea_primitive_data_t data = {0};
    data.box_general.corner_x = corner_x;
    data.box_general.corner_y = corner_y;
    data.box_general.corner_z = corner_z;
    data.box_general.v1_x = v1_x;
    data.box_general.v1_y = v1_y;
    data.box_general.v1_z = v1_z;
    data.box_general.v2_x = v2_x;
    data.box_general.v2_y = v2_y;
    data.box_general.v2_z = v2_z;
    data.box_general.v3_x = v3_x;
    data.box_general.v3_y = v3_y;
    data.box_general.v3_z = v3_z;
    return create_surface_entry(sys, surface_id, ALEA_PRIMITIVE_BOX, &data);
}

int alea_sph_surface(alea_system_t* sys, int surface_id,
                                           double cx, double cy, double cz, double r) {
    alea_primitive_data_t data = {0};
    data.sph.center_x = cx;
    data.sph.center_y = cy;
    data.sph.center_z = cz;
    data.sph.radius = r;
    return create_surface_entry(sys, surface_id, ALEA_PRIMITIVE_SPH, &data);
}

int alea_trc_surface(alea_system_t* sys, int surface_id,
                                           double base_x, double base_y, double base_z,
                                           double height_x, double height_y, double height_z,
                                           double base_radius, double top_radius) {
    alea_primitive_data_t data = {0};
    data.trc.base_x = base_x;
    data.trc.base_y = base_y;
    data.trc.base_z = base_z;
    data.trc.height_x = height_x;
    data.trc.height_y = height_y;
    data.trc.height_z = height_z;
    data.trc.base_radius = base_radius;
    data.trc.top_radius = top_radius;
    return create_surface_entry(sys, surface_id, ALEA_PRIMITIVE_TRC, &data);
}

int alea_ell_surface(alea_system_t* sys, int surface_id,
                                           double v1_x, double v1_y, double v1_z,
                                           double v2_x, double v2_y, double v2_z,
                                           double major_axis_len) {
    alea_primitive_data_t data = {0};
    data.ell.v1_x = v1_x;
    data.ell.v1_y = v1_y;
    data.ell.v1_z = v1_z;
    data.ell.v2_x = v2_x;
    data.ell.v2_y = v2_y;
    data.ell.v2_z = v2_z;
    data.ell.major_axis_len = major_axis_len;
    return create_surface_entry(sys, surface_id, ALEA_PRIMITIVE_ELL, &data);
}

int alea_rec_surface(alea_system_t* sys, int surface_id,
                                           double base_x, double base_y, double base_z,
                                           double height_x, double height_y, double height_z,
                                           double axis1_x, double axis1_y, double axis1_z,
                                           double axis2_x, double axis2_y, double axis2_z) {
    alea_primitive_data_t data = {0};
    data.rec.base_x = base_x;
    data.rec.base_y = base_y;
    data.rec.base_z = base_z;
    data.rec.height_x = height_x;
    data.rec.height_y = height_y;
    data.rec.height_z = height_z;
    data.rec.axis1_x = axis1_x;
    data.rec.axis1_y = axis1_y;
    data.rec.axis1_z = axis1_z;
    data.rec.axis2_x = axis2_x;
    data.rec.axis2_y = axis2_y;
    data.rec.axis2_z = axis2_z;
    return create_surface_entry(sys, surface_id, ALEA_PRIMITIVE_REC, &data);
}

int alea_wed_surface(alea_system_t* sys, int surface_id,
                                           double vertex_x, double vertex_y, double vertex_z,
                                           double v1_x, double v1_y, double v1_z,
                                           double v2_x, double v2_y, double v2_z,
                                           double v3_x, double v3_y, double v3_z) {
    alea_primitive_data_t data = {0};
    data.wed.vertex_x = vertex_x;
    data.wed.vertex_y = vertex_y;
    data.wed.vertex_z = vertex_z;
    data.wed.v1_x = v1_x;
    data.wed.v1_y = v1_y;
    data.wed.v1_z = v1_z;
    data.wed.v2_x = v2_x;
    data.wed.v2_y = v2_y;
    data.wed.v2_z = v2_z;
    data.wed.v3_x = v3_x;
    data.wed.v3_y = v3_y;
    data.wed.v3_z = v3_z;
    return create_surface_entry(sys, surface_id, ALEA_PRIMITIVE_WED, &data);
}

int alea_rhp_surface(alea_system_t* sys, int surface_id,
                                           double base_x, double base_y, double base_z,
                                           double height_x, double height_y, double height_z,
                                           double r1_x, double r1_y, double r1_z,
                                           double r2_x, double r2_y, double r2_z,
                                           double r3_x, double r3_y, double r3_z) {
    alea_primitive_data_t data = {0};
    data.rhp.base_x = base_x;
    data.rhp.base_y = base_y;
    data.rhp.base_z = base_z;
    data.rhp.height_x = height_x;
    data.rhp.height_y = height_y;
    data.rhp.height_z = height_z;
    data.rhp.r1_x = r1_x;
    data.rhp.r1_y = r1_y;
    data.rhp.r1_z = r1_z;
    data.rhp.r2_x = r2_x;
    data.rhp.r2_y = r2_y;
    data.rhp.r2_z = r2_z;
    data.rhp.r3_x = r3_x;
    data.rhp.r3_y = r3_y;
    data.rhp.r3_z = r3_z;
    return create_surface_entry(sys, surface_id, ALEA_PRIMITIVE_RHP, &data);
}

alea_node_id_t alea_halfspace(const alea_system_t* sys, int surface_index, int sense) {
    if (!sys || surface_index < 0 || (size_t)surface_index >= alea_vec_count(&sys->surfaces))
        return ALEA_NODE_ID_INVALID;
    const alea_surface_entry_t* s = &sys->surfaces.data[surface_index];
    return (sense < 0) ? s->neg_node : s->pos_node;
}

/* ============================================================================
 * CSG CONSTRUCTION - OPERATIONS
 * ============================================================================ */

alea_node_id_t alea_union(alea_system_t* sys, alea_node_id_t a, alea_node_id_t b) {
    if (!sys) return ALEA_NODE_ID_INVALID;
    return alea_create_union(sys, a, b);
}

alea_node_id_t alea_intersection(alea_system_t* sys, alea_node_id_t a, alea_node_id_t b) {
    if (!sys) return ALEA_NODE_ID_INVALID;
    return alea_create_intersection(sys, a, b);
}

alea_node_id_t alea_difference(alea_system_t* sys, alea_node_id_t a, alea_node_id_t b) {
    if (!sys) return ALEA_NODE_ID_INVALID;
    return alea_create_difference(sys, a, b);
}

alea_node_id_t alea_complement(alea_system_t* sys, alea_node_id_t a) {
    if (!sys) return ALEA_NODE_ID_INVALID;
    return alea_create_complement(sys, a);
}

alea_node_id_t alea_union_n(alea_system_t* sys, const alea_node_id_t* nodes, size_t count) {
    if (!sys || !nodes || count == 0) return ALEA_NODE_ID_INVALID;
    if (count == 1) return nodes[0];

    alea_node_id_t result = nodes[0];
    for (size_t i = 1; i < count; i++) {
        result = alea_create_union(sys, result, nodes[i]);
        if (result == ALEA_NODE_ID_INVALID) return ALEA_NODE_ID_INVALID;
    }
    return result;
}

alea_node_id_t alea_intersection_n(alea_system_t* sys, const alea_node_id_t* nodes, size_t count) {
    if (!sys || !nodes || count == 0) return ALEA_NODE_ID_INVALID;
    if (count == 1) return nodes[0];

    alea_node_id_t result = nodes[0];
    for (size_t i = 1; i < count; i++) {
        result = alea_create_intersection(sys, result, nodes[i]);
        if (result == ALEA_NODE_ID_INVALID) return ALEA_NODE_ID_INVALID;
    }
    return result;
}

/* ============================================================================
 * CSG CONSTRUCTION - CELLS
 * ============================================================================ */


int alea_set_fill(alea_system_t* sys, int cell_index, int fill_universe, int transform) {
    if (!sys || cell_index < 0 || (size_t)cell_index >= alea_vec_count(&sys->cells))
        return -1;

    sys->cells.data[cell_index].fill_universe = fill_universe;
    sys->cells.data[cell_index].fill_transform = transform;
    return 0;
}

int alea_cell_set_comment(alea_system_t* sys, int cell_index, const char* comment) {
    if (!sys || cell_index < 0 || (size_t)cell_index >= alea_vec_count(&sys->cells))
        return -1;
    alea_cell_entry_t* cell = &sys->cells.data[cell_index];
    free(cell->comments);
    cell->comments = comment ? alea_strdup(comment) : NULL;
    return 0;
}

int alea_cell_set_inline_comment(alea_system_t* sys, int cell_index, const char* comment) {
    if (!sys || cell_index < 0 || (size_t)cell_index >= alea_vec_count(&sys->cells))
        return -1;
    alea_cell_entry_t* cell = &sys->cells.data[cell_index];
    free(cell->inline_comment);
    cell->inline_comment = comment ? alea_strdup(comment) : NULL;
    return 0;
}

/* ============================================================================
 * CELL PROPERTY SETTERS
 * ============================================================================ */

int alea_cell_set_material(alea_system_t* sys, int cell_index, int material_index) {
    if (!sys || cell_index < 0 || (size_t)cell_index >= alea_vec_count(&sys->cells))
        return -1;

    alea_cell_entry_t* cell = &sys->cells.data[cell_index];

    if (material_index == ALEA_MATERIAL_VOID || material_index < 0) {
        cell->material_id = 0;
        cell->material_index = -1;
    } else {
        if ((size_t)material_index >= alea_vec_count(&sys->materials)) {
            alea_set_error_detail(ALEA_ERR_INVALID_ARG,
                "alea_cell_set_material: material index %d out of range (have %zu materials)",
                material_index, alea_vec_count(&sys->materials));
            return -1;
        }
        cell->material_id = sys->materials.data[material_index].material_id;
        cell->material_index = material_index;
    }
    return 0;
}

int alea_cell_set_mixture(alea_system_t* sys, int cell_index, int mixture_id) {
    if (!sys || cell_index < 0 || (size_t)cell_index >= alea_vec_count(&sys->cells))
        return -1;

    if (alea_find_mixture_by_id(sys, mixture_id) < 0) {
        alea_set_error_detail(ALEA_ERR_INVALID_ARG,
            "alea_cell_set_mixture: mixture ID %d not found", mixture_id);
        return -1;
    }

    alea_cell_entry_t* cell = &sys->cells.data[cell_index];
    cell->material_id = mixture_id;
    cell->material_index = -1;  /* Not a direct material */
    return 0;
}

int alea_cell_set_density(alea_system_t* sys, int cell_index, double density) {
    if (!sys || cell_index < 0 || (size_t)cell_index >= alea_vec_count(&sys->cells))
        return -1;

    alea_cell_entry_t* cell = &sys->cells.data[cell_index];
    cell->is_mass_density = (density < 0) ? 1 : 0;
    cell->density = fabs(density);
    return 0;
}

int alea_cell_set_universe(alea_system_t* sys, int cell_index, int universe_id) {
    if (!sys || cell_index < 0 || (size_t)cell_index >= alea_vec_count(&sys->cells))
        return -1;

    sys->cells.data[cell_index].universe_id = universe_id;
    alea_system_invalidate_query_caches(sys, ALEA_CACHE_ALL);
    return 0;
}

/* ============================================================================
 * CELL REMOVAL
 * ============================================================================ */

int alea_cell_remove(alea_system_t* sys, int cell_index) {
    if (!sys || cell_index < 0 || (size_t)cell_index >= alea_vec_count(&sys->cells))
        return -1;

    size_t n_cells = alea_vec_count(&sys->cells);
    alea_cell_entry_t* cell = &sys->cells.data[cell_index];

    /* Notify hook before removal */
    if (sys->on_cell_removed)
        sys->on_cell_removed(sys->cell_hook_userdata, (size_t)cell_index);

    /* Remove from cell ID hashmap */
    cell_hashmap_remove(&sys->cell_index, cell->mc_cell_id);

    /* Free per-cell allocations (skip neighbors if pool-allocated) */
    free(cell->surface_indices);
    if (!sys->neighbor_pool)
        free(cell->neighbors);
    free(cell->lat_fill);
    free(cell->comments);
    free(cell->inline_comment);

    /* Compact: shift remaining cells down */
    size_t last = n_cells - 1;
    if ((size_t)cell_index < last) {
        memmove(&sys->cells.data[cell_index],
                &sys->cells.data[cell_index + 1],
                (last - (size_t)cell_index) * sizeof(alea_cell_entry_t));
    }
    sys->cells.count = last;

    /* Rebuild hashmap — indices shifted */
    cell_hashmap_clear(&sys->cell_index);
    for (size_t i = 0; i < last; i++) {
        cell_hashmap_put(&sys->cell_index, sys->cells.data[i].mc_cell_id, (int)i);
    }

    alea_system_invalidate_query_caches(sys, ALEA_CACHE_ALL);
    return 0;
}

/* ============================================================================
 * UNIVERSE OPERATIONS
 * ============================================================================ */

int alea_flatten(alea_system_t* sys, int universe_id) {
    if (!sys) return -1;
    alea_flatten_config_t config = ALEA_FLATTEN_DEFAULT;
    config.starting_universe_id = universe_id;
    return alea_flatten_in_place(sys, &config);
}


alea_system_t* alea_extract_universe(const alea_system_t* sys, int universe_id) {
    if (!sys) return NULL;

    alea_system_t* extracted = alea_system_create();
    if (!extracted) return NULL;

    primitive_remap_t* remap = alea_create_remap_table(alea_vec_count(&sys->primitives));
    if (!remap) {
        alea_system_destroy(extracted);
        return NULL;
    }

    /* Find cells in the target universe */
    for (size_t i = 0; i < alea_vec_count(&sys->cells); i++) {
        const alea_cell_entry_t* cell = &sys->cells.data[i];
        if (cell->universe_id != universe_id) continue;

        /* Deep copy node tree */
        alea_node_id_t new_root = ALEA_NODE_ID_INVALID;
        if (cell->root_node_id != ALEA_NODE_ID_INVALID) {
            new_root = alea_clone_tree_to_system(extracted, sys, cell->root_node_id, remap);
        }

        int idx = alea_add_cell(extracted, cell->mc_cell_id, new_root,
                               ALEA_MATERIAL_VOID, cell->density, cell->universe_id);
        if (idx >= 0) {
            /* Copy core cell fields */
            alea_cell_entry_t* dst_cell = &extracted->cells.data[idx];
            dst_cell->material_id = cell->material_id;
            dst_cell->material_index = cell->material_index;
            dst_cell->is_mass_density = cell->is_mass_density;
            dst_cell->original_root_node_id = cell->original_root_node_id;
            dst_cell->temperature = cell->temperature;
            dst_cell->has_temperature = cell->has_temperature;
            dst_cell->fill_universe = cell->fill_universe;
            dst_cell->fill_transform = cell->fill_transform;
            dst_cell->comments = cell->comments ? alea_strdup(cell->comments) : NULL;
            dst_cell->inline_comment = cell->inline_comment ? alea_strdup(cell->inline_comment) : NULL;
        }
    }

    /* Copy referenced auxiliary data */
    alea_copy_surfaces_with_remap(extracted, sys, remap);
    alea_destroy_remap_table(remap);
    alea_copy_referenced_materials(extracted, sys);
    alea_copy_referenced_mixtures(extracted, sys);
    alea_copy_referenced_transforms(extracted, sys);
    alea_copy_referenced_cell_refs(extracted, sys);

    extracted->source = sys->source;
    return extracted;
}

int alea_merge(alea_system_t* target, const alea_system_t* source, int id_offset) {
    if (!target || !source) return -1;

    size_t node_offset = alea_vec_count(&target->nodes);
    size_t source_prim_count = alea_vec_count(&source->primitives);
    uint32_t* prim_map = NULL;
    int8_t* prim_inverted = NULL;

    if (source_prim_count > 0) {
        prim_map = malloc(source_prim_count * sizeof(uint32_t));
        prim_inverted = calloc(source_prim_count, sizeof(int8_t));
        if (!prim_map || !prim_inverted) {
            free(prim_map);
            free(prim_inverted);
            return -1;
        }
    }

    /* Copy primitives */
    for (size_t i = 0; i < source_prim_count; i++) {
        const alea_primitive_entry_t* prim = &source->primitives.data[i];
        alea_primitive_data_t data;
        int8_t inverted = 0;
        if (!alea_primitive_copy_data(source, (uint32_t)i, &data)) {
            free(prim_map);
            free(prim_inverted);
            return -1;
        }
        uint32_t new_id = alea_get_or_create_primitive(
            target, prim->type, &data, &inverted);
        if (new_id == UINT32_MAX) {
            free(prim_map);
            free(prim_inverted);
            return -1;
        }
        prim_map[i] = new_id;
        prim_inverted[i] = inverted;
    }

    /* Copy nodes with adjusted references */
    for (size_t i = 0; i < alea_vec_count(&source->nodes); i++) {
        alea_node_t node = source->nodes.data[i];
        alea_operation_t op = ALEA_GET_OPERATION(&node);

        if (op == ALEA_OP_PRIMITIVE) {
            uint32_t old_id = node.primitive.primitive_id;
            if (old_id >= source_prim_count) {
                free(prim_map);
                free(prim_inverted);
                return -1;
            }
            node.primitive.primitive_id = prim_map[old_id];
            if (prim_inverted[old_id]) {
                node.primitive.inverted = !node.primitive.inverted;
            }
        } else {
            node.operation.left += (uint32_t)node_offset;
            if (op != ALEA_OP_COMPLEMENT) {
                node.operation.right += (uint32_t)node_offset;
            }
        }
        alea_vec_push(&target->nodes, node, alea_node_t);
    }

    /* Copy cells with adjusted IDs and node references */
    int cells_added = 0;
    for (size_t i = 0; i < alea_vec_count(&source->cells); i++) {
        alea_cell_entry_t cell = source->cells.data[i];
        cell.mc_cell_id += id_offset;
        if (cell.root_node_id != ALEA_NODE_ID_INVALID) {
            cell.root_node_id += (uint32_t)node_offset;
        }
        if (cell.material_id != 0) {
            cell.material_id += id_offset;
        }
        /* NULL out cached index data to prevent double-free */
        cell.surface_indices = NULL;
        cell.surface_index_count = 0;
        cell.neighbors = NULL;
        cell.neighbor_count = 0;
        /* Deep-copy owned lattice fill array */
        if (cell.lat_fill && cell.lat_fill_count > 0) {
            const alea_cell_entry_t* src = &source->cells.data[i];
            int* copy = malloc(src->lat_fill_count * sizeof(int));
            if (copy) {
                memcpy(copy, src->lat_fill, src->lat_fill_count * sizeof(int));
            }
            cell.lat_fill = copy;
        } else {
            cell.lat_fill = NULL;
        }
        /* Deep-copy comment strings */
        cell.comments = cell.comments ? alea_strdup(cell.comments) : NULL;
        cell.inline_comment = cell.inline_comment ? alea_strdup(cell.inline_comment) : NULL;
        int idx = (int)alea_vec_count(&target->cells);
        alea_vec_push(&target->cells, cell, alea_cell_entry_t);
        cell_hashmap_put(&target->cell_index, cell.mc_cell_id, idx);
        cells_added++;
    }

    /* Copy surfaces with adjusted IDs */
    for (size_t i = 0; i < alea_vec_count(&source->surfaces); i++) {
        alea_surface_entry_t surf = source->surfaces.data[i];
        surf.mc_surface_id += id_offset;
        surf.pos_node += (uint32_t)node_offset;
        surf.neg_node += (uint32_t)node_offset;
        if (surf.primitive_id >= source_prim_count) {
            free(prim_map);
            free(prim_inverted);
            return -1;
        }
        surf.primitive_id = prim_map[surf.primitive_id];
        alea_vec_push(&target->surfaces, surf, alea_surface_entry_t);
    }

    /* Mark universe index as needing rebuild */
    target->universe_index_built = false;

    free(prim_map);
    free(prim_inverted);
    return cells_added;
}

/* ============================================================================
 * INFORMATION
 * ============================================================================ */

size_t alea_cell_count(const alea_system_t* sys) {
    return sys ? alea_vec_count(&sys->cells) : 0;
}

size_t alea_surface_count(const alea_system_t* sys) {
    return sys ? alea_vec_count(&sys->surfaces) : 0;
}

size_t alea_universe_count(const alea_system_t* sys) {
    return sys ? alea_vec_count(&sys->universes) : 0;
}

void alea_stats(const alea_system_t* sys, alea_stats_t* stats) {
    if (!sys || !stats) return;
    *stats = sys->stats;
}

void alea_print_summary(const alea_system_t* sys) {
    if (!sys) return;

    ALEA_LOG_INFO("Alea System Summary:");
    ALEA_LOG_INFO("  Cells:      %zu", alea_vec_count(&sys->cells));
    ALEA_LOG_INFO("  Surfaces:   %zu", alea_vec_count(&sys->surfaces));
    ALEA_LOG_INFO("  Primitives: %zu", alea_vec_count(&sys->primitives));
    ALEA_LOG_INFO("  Nodes:      %zu", alea_vec_count(&sys->nodes));
    ALEA_LOG_INFO("  Universes:  %zu", alea_vec_count(&sys->universes));
    ALEA_LOG_INFO("  Transforms: %zu", alea_vec_count(&sys->transforms));
}

/* ============================================================================
 * TREE DEBUG PRINTING
 * ============================================================================ */

static const char* tree_op_name(alea_operation_t op) {
    switch (op) {
        case ALEA_OP_PRIMITIVE:    return "PRIMITIVE";
        case ALEA_OP_UNION:        return "UNION";
        case ALEA_OP_INTERSECTION: return "INTERSECTION";
        case ALEA_OP_DIFFERENCE:   return "DIFFERENCE";
        default:                   return "UNKNOWN";
    }
}

static const char* tree_prim_name(alea_primitive_type_t type) {
    switch (type) {
        case ALEA_PRIMITIVE_PLANE:      return "PLANE";
        case ALEA_PRIMITIVE_SPHERE:     return "SPHERE";
        case ALEA_PRIMITIVE_CYLINDER_X: return "CYLINDER_X";
        case ALEA_PRIMITIVE_CYLINDER_Y: return "CYLINDER_Y";
        case ALEA_PRIMITIVE_CYLINDER_Z: return "CYLINDER_Z";
        case ALEA_PRIMITIVE_CONE_X:     return "CONE_X";
        case ALEA_PRIMITIVE_CONE_Y:     return "CONE_Y";
        case ALEA_PRIMITIVE_CONE_Z:     return "CONE_Z";
        case ALEA_PRIMITIVE_RPP:        return "RPP";
        case ALEA_PRIMITIVE_QUADRIC:    return "QUADRIC";
        case ALEA_PRIMITIVE_TORUS_X:    return "TORUS_X";
        case ALEA_PRIMITIVE_TORUS_Y:    return "TORUS_Y";
        case ALEA_PRIMITIVE_TORUS_Z:    return "TORUS_Z";
        default:                        return "UNKNOWN";
    }
}

static void tree_print_impl(const alea_system_t* sys, alea_node_id_t node_id,
                             const char* prefix, bool is_last) {
    if (node_id == ALEA_NODE_ID_INVALID || node_id >= alea_vec_count(&sys->nodes))
        return;

    const alea_node_t* node = &sys->nodes.data[node_id];
    alea_operation_t op = ALEA_GET_OPERATION(node);

    if (op == ALEA_OP_PRIMITIVE) {
        ALEA_LOG_INFO("%s%s%s", prefix, is_last ? "└── " : "├── ",
                      tree_prim_name(node->primitive.prim_type));
    } else {
        ALEA_LOG_INFO("%s%s%s", prefix, is_last ? "└── " : "├── ", tree_op_name(op));

        char new_prefix[256];
        snprintf(new_prefix, sizeof(new_prefix), "%s%s",
                 prefix, is_last ? "    " : "│   ");

        tree_print_impl(sys, node->operation.left, new_prefix, false);
        tree_print_impl(sys, node->operation.right, new_prefix, true);
    }
}

void alea_tree_print(const alea_system_t* sys, alea_node_id_t node_id) {
    if (!sys) return;
    ALEA_LOG_INFO("CSG Tree:");
    tree_print_impl(sys, node_id, "", true);
}

/* ============================================================================
 * VOLUME ESTIMATION & BOUNDING BOX TIGHTENING
 * ============================================================================ */

/* alea_compute_bounding_sphere, alea_tighten_cell_bbox, alea_tighten_all_bboxes
 * are in primitives/bbox.c (already included via primitives/bbox.h) */

/* Weak stubs in alea_module_stubs.c, strong overrides in raycast_api.c */
int alea_estimate_cell_volumes(alea_system_t* sys,
                               double ox, double oy, double oz,
                               double radius, int n_rays,
                               double* volumes, double* rel_errors);
int alea_estimate_path_volumes(alea_system_t* sys,
                               int n_rays,
                               double* volumes,
                               double* rel_errors);
int alea_remove_cells_by_volume(alea_system_t* sys,
                                const double* volumes,
                                double threshold);

typedef struct {
    uint64_t hash;
    uint64_t path_id;
    uint8_t occupied;
} volume_path_hash_slot_t;

typedef struct alea_volume_path_index {
    alea_volume_path_t* paths;
    size_t count;
    size_t capacity;
    size_t max_path_count;
    uint64_t geometry_generation;
    volume_path_hash_slot_t* slots;
    size_t slot_capacity;
} alea_volume_path_index_t;

typedef struct {
    alea_volume_path_t* out;
    size_t max;
    size_t count;
    alea_volume_path_index_t* index;
    int failed;
} volume_path_enum_ctx_t;

typedef struct {
    int ancestor_count;
    int ancestor_cell_indices[ALEA_VOLUME_PATH_MAX_DEPTH];
    int ancestor_universe_ids[ALEA_VOLUME_PATH_MAX_DEPTH];
    int lattice_step_count;
    alea_volume_lattice_step_t lattice_steps[ALEA_VOLUME_PATH_MAX_DEPTH];
    alea_matrix_t transform;
    bool truncated;
} volume_path_builder_t;

static size_t volume_path_max_count_from_env(void) {
    const char* env = getenv("ALEA_VOLUME_PATH_MAX_COUNT");
    if (!env || !env[0]) return SIZE_MAX;
    char* end = NULL;
    unsigned long long value = alea_strtoull(env, &end, 10);
    if (end == env || value == 0) return SIZE_MAX;
    if (value > (unsigned long long)SIZE_MAX) return SIZE_MAX;
    return (size_t)value;
}

static void volume_path_copy_transform(double dst[12], const alea_matrix_t* mat) {
    if (!dst) return;
    if (!mat) {
        memset(dst, 0, 12 * sizeof(double));
        dst[0] = dst[5] = dst[10] = 1.0;
        return;
    }
    memcpy(dst, mat->m, 12 * sizeof(double));
}

static int volume_path_fill_matrix(const alea_system_t* sys,
                                   const alea_cell_entry_t* cell,
                                   alea_matrix_t* out) {
    if (!cell || !out) return -1;
    if (cell->fill_transform != 0) {
        const alea_transform_t* tr = alea_get_transform(sys, cell->fill_transform);
        if (!tr) return -1;
        if (!alea_matrix_from_mcnp(out, tr->cosines, tr->value_count, false)) return -1;
    } else {
        alea_matrix_identity(out);
    }
    return alea_matrix_invert(out) ? 0 : -1;
}

static void volume_path_push_ancestor(volume_path_builder_t* b,
                                      int cell_index,
                                      int universe_id) {
    if (!b) return;
    if (b->ancestor_count >= ALEA_VOLUME_PATH_MAX_DEPTH) {
        b->truncated = true;
        return;
    }
    b->ancestor_cell_indices[b->ancestor_count] = cell_index;
    b->ancestor_universe_ids[b->ancestor_count] = universe_id;
    b->ancestor_count++;
}

static void volume_path_push_lattice(volume_path_builder_t* b,
                                     const alea_volume_lattice_step_t* step) {
    if (!b || !step) return;
    if (b->lattice_step_count >= ALEA_VOLUME_PATH_MAX_DEPTH) {
        b->truncated = true;
        return;
    }
    b->lattice_steps[b->lattice_step_count++] = *step;
}

static int volume_path_index_append(alea_volume_path_index_t* idx,
                                    const alea_volume_path_t* path) {
    if (!idx || !path) return -1;
    if (idx->count >= idx->max_path_count) {
        alea_set_error_detail(ALEA_ERR_INVALID_STATE,
                              "volume path count exceeds ALEA_VOLUME_PATH_MAX_COUNT=%zu",
                              idx->max_path_count);
        return -1;
    }
    if (idx->count == idx->capacity) {
        size_t next = idx->capacity ? idx->capacity * 2 : 16;
        alea_volume_path_t* grown = realloc(idx->paths, next * sizeof(*grown));
        if (!grown) return -1;
        idx->paths = grown;
        idx->capacity = next;
    }
    idx->paths[idx->count] = *path;
    idx->paths[idx->count].path_id = (uint64_t)idx->count;
    idx->count++;
    return 0;
}

static void volume_path_emit(volume_path_enum_ctx_t* ctx,
                             const alea_system_t* sys,
                             int cell_index,
                             int depth,
                             const volume_path_builder_t* b) {
    if (!ctx || !sys || cell_index < 0 ||
        (size_t)cell_index >= alea_vec_count(&sys->cells)) {
        return;
    }

    size_t id = ctx->count++;

    const alea_cell_entry_t* cell = &sys->cells.data[cell_index];
    alea_volume_path_t local;
    alea_volume_path_t* path = &local;
    memset(path, 0, sizeof(*path));
    path->path_id = (uint64_t)id;
    path->terminal_cell_index = cell_index;
    path->terminal_cell_id = cell->mc_cell_id;
    path->material_id = cell->material_id;
    path->universe_id = cell->universe_id;
    path->depth = depth;

    int ancestor_count = b ? b->ancestor_count : 0;
    if (ancestor_count > ALEA_VOLUME_PATH_MAX_DEPTH) ancestor_count = ALEA_VOLUME_PATH_MAX_DEPTH;
    path->ancestor_count = (uint8_t)ancestor_count;
    for (int i = 0; i < ancestor_count; i++) {
        path->ancestor_cell_indices[i] = b->ancestor_cell_indices[i];
        path->ancestor_universe_ids[i] = b->ancestor_universe_ids[i];
    }

    int lattice_count = b ? b->lattice_step_count : 0;
    if (lattice_count > ALEA_VOLUME_PATH_MAX_DEPTH) lattice_count = ALEA_VOLUME_PATH_MAX_DEPTH;
    path->lattice_step_count = (uint8_t)lattice_count;
    for (int i = 0; i < lattice_count; i++) {
        path->lattice_steps[i] = b->lattice_steps[i];
    }

    path->truncated = (uint8_t)(b && b->truncated);
    volume_path_copy_transform(path->world_to_local, b ? &b->transform : NULL);

    if (ctx->index && volume_path_index_append(ctx->index, path) != 0) {
        ctx->failed = 1;
    }
    if (ctx->out && id < ctx->max) {
        ctx->out[id] = *path;
    }
}

static int volume_paths_enumerate_universe(alea_system_t* sys,
                                           int universe_id,
                                           int depth,
                                           volume_path_builder_t* builder,
                                           volume_path_enum_ctx_t* ctx);

static int volume_paths_enumerate_lattice(alea_system_t* sys,
                                          const alea_cell_entry_t* cell,
                                          int cell_index,
                                          int depth,
                                          volume_path_builder_t* builder,
                                          volume_path_enum_ctx_t* ctx) {
    int ni = cell->lat_fill_dims[1] - cell->lat_fill_dims[0] + 1;
    int nj = cell->lat_fill_dims[3] - cell->lat_fill_dims[2] + 1;
    int nk = cell->lat_fill_dims[5] - cell->lat_fill_dims[4] + 1;
    if (ni <= 0 || nj <= 0 || nk <= 0 || !cell->lat_fill) return -1;

    for (int oi = 0; oi < ni; oi++) {
        for (int oj = 0; oj < nj; oj++) {
            for (int ok = 0; ok < nk; ok++) {
                size_t linear = (size_t)(oi * nj * nk + oj * nk + ok);
                if (linear >= cell->lat_fill_count) continue;
                int fill_universe = cell->lat_fill[linear];
                if (fill_universe <= 0) continue;

                volume_path_builder_t child = *builder;
                volume_path_push_ancestor(&child, cell_index, cell->universe_id);
                alea_volume_lattice_step_t step = {
                    .lattice_cell_index = cell_index,
                    .fill_universe = fill_universe,
                    .i = cell->lat_fill_dims[0] + oi,
                    .j = cell->lat_fill_dims[2] + oj,
                    .k = cell->lat_fill_dims[4] + ok,
                    .linear_index = (int)linear
                };
                volume_path_push_lattice(&child, &step);

                alea_matrix_t translation;
                alea_matrix_identity(&translation);
                if (cell->lat_type == 2) {
                    double p = cell->lat_pitch[0] > 0.0 ? cell->lat_pitch[0] : 1.0;
                    translation.m[3] = step.i * p + step.j * p * 0.5;
                    translation.m[7] = step.j * p * sqrt(3.0) * 0.5;
                    translation.m[11] = cell->lat_fill_zero_element_coords
                        ? step.k * cell->lat_pitch[2]
                        : ((nk == 1) ? 0.0
                           : cell->lat_lower_left[2] + (ok + 0.5) * cell->lat_pitch[2]);
                } else {
                    if (cell->lat_fill_zero_element_coords) {
                        translation.m[3] = step.i * cell->lat_pitch[0];
                        translation.m[7] = step.j * cell->lat_pitch[1];
                        translation.m[11] = step.k * cell->lat_pitch[2];
                    } else {
                        translation.m[3] = cell->lat_lower_left[0] + (oi + 0.5) * cell->lat_pitch[0];
                        translation.m[7] = cell->lat_lower_left[1] + (oj + 0.5) * cell->lat_pitch[1];
                        translation.m[11] = cell->lat_lower_left[2] + (ok + 0.5) * cell->lat_pitch[2];
                    }
                }
                translation.has_inverse = false;

                alea_matrix_t composed;
                alea_matrix_multiply(&composed, &builder->transform, &translation);
                if (!alea_matrix_invert(&composed)) return -1;
                child.transform = composed;

                if (volume_paths_enumerate_universe(sys, fill_universe,
                                                    depth + 1, &child, ctx) != 0) {
                    return -1;
                }
            }
        }
    }
    return 0;
}

static int volume_paths_enumerate_universe(alea_system_t* sys,
                                           int universe_id,
                                           int depth,
                                           volume_path_builder_t* builder,
                                           volume_path_enum_ctx_t* ctx) {
    if (!sys || !builder || !ctx) return -1;
    if (depth >= ALEA_VOLUME_PATH_MAX_DEPTH) {
        builder->truncated = true;
        return 0;
    }

    if (!sys->universe_index_built && alea_build_universe_index(sys) != 0) return -1;
    const alea_universe_t* univ = alea_get_universe(sys, universe_id);
    if (!univ) return 0;

    for (size_t i = 0; i < univ->cell_indices.count; i++) {
        int cell_index = univ->cell_indices.data[i];
        if (cell_index < 0 || (size_t)cell_index >= alea_vec_count(&sys->cells)) continue;
        const alea_cell_entry_t* cell = &sys->cells.data[cell_index];

        if (cell->lat_type != 0 && cell->lat_fill) {
            if (volume_paths_enumerate_lattice(sys, cell, cell_index, depth,
                                               builder, ctx) != 0) {
                return -1;
            }
            continue;
        }

        if (cell->fill_universe > 0) {
            alea_matrix_t fill;
            if (volume_path_fill_matrix(sys, cell, &fill) != 0) return -1;
            volume_path_builder_t child = *builder;
            volume_path_push_ancestor(&child, cell_index, cell->universe_id);
            alea_matrix_t composed;
            alea_matrix_multiply(&composed, &builder->transform, &fill);
            if (!alea_matrix_invert(&composed)) return -1;
            child.transform = composed;
            if (volume_paths_enumerate_universe(sys, cell->fill_universe,
                                                depth + 1, &child, ctx) != 0) {
                return -1;
            }
            continue;
        }

        volume_path_emit(ctx, sys, cell_index, depth, builder);
    }
    return 0;
}

static int alea_volume_path_index_ensure(alea_system_t* sys);
static const alea_volume_path_t*
alea_volume_path_index_find(const alea_volume_path_index_t* idx,
                            const alea_volume_path_t* structural_path);

size_t alea_volume_path_count(alea_system_t* sys) {
    if (!sys || alea_volume_path_index_ensure(sys) != 0 || !sys->volume_path_index)
        return 0;
    return sys->volume_path_index->count;
}

size_t alea_volume_paths_get(alea_system_t* sys,
                             alea_volume_path_t* out_paths,
                             size_t max_paths) {
    if (!out_paths && max_paths > 0) return 0;
    if (!sys || alea_volume_path_index_ensure(sys) != 0 || !sys->volume_path_index)
        return 0;
    size_t count = sys->volume_path_index->count;
    size_t n = count < max_paths ? count : max_paths;
    if (out_paths && n > 0)
        memcpy(out_paths, sys->volume_path_index->paths, n * sizeof(*out_paths));
    return count;
}

static int volume_lattice_lookup_step(const alea_cell_entry_t* cell,
                                      double px,
                                      double py,
                                      double pz,
                                      alea_volume_lattice_step_t* step,
                                      double* ox,
                                      double* oy,
                                      double* oz) {
    if (!cell || !step || !ox || !oy || !oz || !cell->lat_fill) return -1;

    int ni = cell->lat_fill_dims[1] - cell->lat_fill_dims[0] + 1;
    int nj = cell->lat_fill_dims[3] - cell->lat_fill_dims[2] + 1;
    int nk = cell->lat_fill_dims[5] - cell->lat_fill_dims[4] + 1;
    if (ni <= 0 || nj <= 0 || nk <= 0) return -1;

    int oi, oj, ok;
    int i, j, k;
    if (cell->lat_type == 2) {
        double p = cell->lat_pitch[0];
        if (p <= 0.0) return -1;

        double fj = py / (p * sqrt(3.0) * 0.5);
        double fi = px / p - fj * 0.5;
        double cx = fi;
        double cz = fj;
        double cy = -fi - fj;

        int ri = (int)round(cx);
        int rj = (int)round(cy);
        int rk = (int)round(cz);

        double dx = fabs(ri - cx);
        double dy = fabs(rj - cy);
        double dz = fabs(rk - cz);
        if (dx > dy && dx > dz) {
            ri = -rj - rk;
        } else if (dy > dz) {
            rj = -ri - rk;
        } else {
            rk = -ri - rj;
        }

        i = ri;
        j = rk;
        k = (nk == 1 && !cell->lat_fill_repeating) ? cell->lat_fill_dims[4]
                      : cell->lat_fill_dims[4] +
                        (int)floor((pz - cell->lat_lower_left[2]) / cell->lat_pitch[2]);
        oi = i - cell->lat_fill_dims[0];
        oj = j - cell->lat_fill_dims[2];
        ok = k - cell->lat_fill_dims[4];
        *ox = ri * p + rk * p * 0.5;
        *oy = rk * p * sqrt(3.0) * 0.5;
        *oz = cell->lat_fill_zero_element_coords
            ? k * cell->lat_pitch[2]
            : ((nk == 1 && !cell->lat_fill_repeating) ? 0.0
               : cell->lat_lower_left[2] + (ok + 0.5) * cell->lat_pitch[2]);
    } else {
        oi = (int)floor((px - cell->lat_lower_left[0]) / cell->lat_pitch[0]);
        oj = (cell->lat_fill_repeating || nj > 1)
            ? (int)floor((py - cell->lat_lower_left[1]) / cell->lat_pitch[1]) : 0;
        ok = (cell->lat_fill_repeating || nk > 1)
            ? (int)floor((pz - cell->lat_lower_left[2]) / cell->lat_pitch[2]) : 0;
        i = cell->lat_fill_dims[0] + oi;
        j = cell->lat_fill_dims[2] + oj;
        k = cell->lat_fill_dims[4] + ok;
        if (cell->lat_fill_zero_element_coords) {
            *ox = i * cell->lat_pitch[0];
            *oy = j * cell->lat_pitch[1];
            *oz = k * cell->lat_pitch[2];
        } else {
            *ox = cell->lat_lower_left[0] + (oi + 0.5) * cell->lat_pitch[0];
            *oy = cell->lat_lower_left[1] + (oj + 0.5) * cell->lat_pitch[1];
            *oz = cell->lat_lower_left[2] + (ok + 0.5) * cell->lat_pitch[2];
        }
    }

    if (!cell->lat_fill_repeating &&
        (oi < 0 || oi >= ni || oj < 0 || oj >= nj || ok < 0 || ok >= nk)) return -1;
    size_t linear = cell->lat_fill_repeating
        ? 0 : (size_t)(oi * nj * nk + oj * nk + ok);
    if (linear >= cell->lat_fill_count) return -1;

    step->lattice_cell_index = -1;
    step->fill_universe = cell->lat_fill[linear];
    step->i = i;
    step->j = j;
    step->k = k;
    step->linear_index = (int)linear;
    return step->fill_universe;
}

static bool volume_paths_same_identity(const alea_volume_path_t* a,
                                       const alea_volume_path_t* b) {
    if (!a || !b) return false;
    if (a->terminal_cell_index != b->terminal_cell_index ||
        a->universe_id != b->universe_id ||
        a->depth != b->depth ||
        a->ancestor_count != b->ancestor_count ||
        a->lattice_step_count != b->lattice_step_count) {
        return false;
    }
    for (uint8_t i = 0; i < a->ancestor_count; i++) {
        if (a->ancestor_cell_indices[i] != b->ancestor_cell_indices[i] ||
            a->ancestor_universe_ids[i] != b->ancestor_universe_ids[i]) {
            return false;
        }
    }
    for (uint8_t i = 0; i < a->lattice_step_count; i++) {
        const alea_volume_lattice_step_t* x = &a->lattice_steps[i];
        const alea_volume_lattice_step_t* y = &b->lattice_steps[i];
        if (x->lattice_cell_index != y->lattice_cell_index ||
            x->fill_universe != y->fill_universe ||
            x->i != y->i || x->j != y->j || x->k != y->k ||
            x->linear_index != y->linear_index) {
            return false;
        }
    }
    return true;
}

void alea_volume_path_index_free(struct alea_volume_path_index* idx) {
    if (!idx) return;
    free(idx->paths);
    free(idx->slots);
    free(idx);
}

static uint64_t volume_path_hash_mix_u64(uint64_t h, uint64_t v) {
    h ^= v;
    h *= 1099511628211ULL;
    return h;
}

static uint64_t volume_path_hash_mix_i32(uint64_t h, int v) {
    return volume_path_hash_mix_u64(h, (uint64_t)(uint32_t)v);
}

static uint64_t volume_path_hash_identity(const alea_volume_path_t* path) {
    uint64_t h = 1469598103934665603ULL;
    h = volume_path_hash_mix_i32(h, path->terminal_cell_index);
    h = volume_path_hash_mix_i32(h, path->universe_id);
    h = volume_path_hash_mix_i32(h, path->depth);
    h = volume_path_hash_mix_u64(h, path->ancestor_count);
    for (uint8_t i = 0; i < path->ancestor_count; i++) {
        h = volume_path_hash_mix_i32(h, path->ancestor_cell_indices[i]);
        h = volume_path_hash_mix_i32(h, path->ancestor_universe_ids[i]);
    }
    h = volume_path_hash_mix_u64(h, path->lattice_step_count);
    for (uint8_t i = 0; i < path->lattice_step_count; i++) {
        const alea_volume_lattice_step_t* step = &path->lattice_steps[i];
        h = volume_path_hash_mix_i32(h, step->lattice_cell_index);
        h = volume_path_hash_mix_i32(h, step->fill_universe);
        h = volume_path_hash_mix_i32(h, step->i);
        h = volume_path_hash_mix_i32(h, step->j);
        h = volume_path_hash_mix_i32(h, step->k);
        h = volume_path_hash_mix_i32(h, step->linear_index);
    }
    return h ? h : 1;
}

static size_t volume_path_hash_capacity(size_t count) {
    size_t cap = 16;
    while (cap < count * 2) cap *= 2;
    return cap;
}

static int volume_path_index_build_hash(alea_volume_path_index_t* idx) {
    if (!idx) return -1;
    free(idx->slots);
    idx->slots = NULL;
    idx->slot_capacity = 0;
    if (idx->count == 0) return 0;

    idx->slot_capacity = volume_path_hash_capacity(idx->count);
    idx->slots = calloc(idx->slot_capacity, sizeof(*idx->slots));
    if (!idx->slots) {
        idx->slot_capacity = 0;
        return -1;
    }

    for (size_t i = 0; i < idx->count; i++) {
        uint64_t hash = volume_path_hash_identity(&idx->paths[i]);
        size_t mask = idx->slot_capacity - 1;
        size_t pos = (size_t)hash & mask;
        for (;;) {
            volume_path_hash_slot_t* slot = &idx->slots[pos];
            if (!slot->occupied) {
                slot->occupied = 1;
                slot->hash = hash;
                slot->path_id = idx->paths[i].path_id;
                break;
            }
            if (slot->hash == hash &&
                slot->path_id < idx->count &&
                volume_paths_same_identity(&idx->paths[slot->path_id],
                                           &idx->paths[i])) {
                alea_set_error_detail(ALEA_ERR_INVALID_STATE,
                                      "duplicate volume path identity");
                return -1;
            }
            pos = (pos + 1) & mask;
        }
    }
    return 0;
}

static const alea_volume_path_t*
alea_volume_path_index_find(const alea_volume_path_index_t* idx,
                            const alea_volume_path_t* structural_path) {
    if (!idx || !structural_path || !idx->slots || idx->slot_capacity == 0)
        return NULL;

    uint64_t hash = volume_path_hash_identity(structural_path);
    size_t mask = idx->slot_capacity - 1;
    size_t pos = (size_t)hash & mask;
    for (;;) {
        const volume_path_hash_slot_t* slot = &idx->slots[pos];
        if (!slot->occupied) return NULL;
        if (slot->hash == hash && slot->path_id < idx->count) {
            const alea_volume_path_t* candidate = &idx->paths[slot->path_id];
            if (volume_paths_same_identity(candidate, structural_path))
                return candidate;
        }
        pos = (pos + 1) & mask;
    }
}

static int volume_path_index_build(alea_system_t* sys,
                                   alea_volume_path_index_t** out_idx) {
    if (!sys || !out_idx) return -1;
    *out_idx = NULL;

    alea_volume_path_index_t* idx = calloc(1, sizeof(*idx));
    if (!idx) return -1;
    idx->geometry_generation = alea_system_geometry_generation(sys);
    idx->max_path_count = volume_path_max_count_from_env();

    if (alea_system_prepare_query_caches(sys, ALEA_CACHE_HIER_SPATIAL) != 0) {
        alea_volume_path_index_free(idx);
        return -1;
    }

    volume_path_enum_ctx_t ctx = {.index = idx, .count = 0};
    volume_path_builder_t builder;
    memset(&builder, 0, sizeof(builder));
    alea_matrix_identity(&builder.transform);
    if (volume_paths_enumerate_universe(sys, 0, 0, &builder, &ctx) != 0 ||
        ctx.failed) {
        alea_volume_path_index_free(idx);
        return -1;
    }

    if (volume_path_index_build_hash(idx) != 0) {
        alea_volume_path_index_free(idx);
        return -1;
    }

    *out_idx = idx;
    return 0;
}

static bool volume_path_index_fresh(const alea_system_t* sys,
                                    const alea_volume_path_index_t* idx) {
    if (!sys || !idx) return false;
    return idx->geometry_generation == alea_system_geometry_generation(sys);
}

static int alea_volume_path_index_ensure(alea_system_t* sys) {
    if (!sys) return -1;
    if (volume_path_index_fresh(sys, sys->volume_path_index)) return 0;

    alea_volume_path_index_t* rebuilt = NULL;
    if (volume_path_index_build(sys, &rebuilt) != 0) return -1;

    alea_volume_path_index_free(sys->volume_path_index);
    sys->volume_path_index = rebuilt;
    return 0;
}

static int volume_path_lookup_universe(alea_system_t* sys,
                                       int universe_id,
                                       double x,
                                       double y,
                                       double z,
                                       int depth,
                                       volume_path_builder_t* builder,
                                       alea_volume_path_t* out,
                                       bool* found) {
    if (!sys || !builder || !out || !found) return -1;
    if (depth >= ALEA_VOLUME_PATH_MAX_DEPTH) {
        builder->truncated = true;
        return 0;
    }

    const alea_universe_t* univ = alea_get_universe(sys, universe_id);
    if (!univ) return -1;

    for (size_t i = 0; i < univ->cell_indices.count; i++) {
        int cell_index = univ->cell_indices.data[i];
        if (cell_index < 0 || (size_t)cell_index >= alea_vec_count(&sys->cells)) continue;
        const alea_cell_entry_t* cell = &sys->cells.data[cell_index];
        if (cell->lat_type != 0 && cell->lat_fill) {
            alea_volume_lattice_step_t step;
            double ox = 0.0, oy = 0.0, oz = 0.0;
            int fill_universe = volume_lattice_lookup_step(cell, x, y, z,
                                                           &step, &ox, &oy, &oz);
            if (fill_universe <= 0) continue;
            step.lattice_cell_index = cell_index;

            volume_path_builder_t child = *builder;
            volume_path_push_ancestor(&child, cell_index, cell->universe_id);
            volume_path_push_lattice(&child, &step);

            alea_matrix_t translation;
            alea_matrix_identity(&translation);
            translation.m[3] = ox;
            translation.m[7] = oy;
            translation.m[11] = oz;
            translation.has_inverse = false;

            alea_matrix_t composed;
            alea_matrix_multiply(&composed, &builder->transform, &translation);
            if (!alea_matrix_invert(&composed)) return -1;
            child.transform = composed;

            if (volume_path_lookup_universe(sys, fill_universe,
                                            x - ox, y - oy, z - oz,
                                            depth + 1, &child, out, found) != 0) {
                return -1;
            }
            if (*found) return 0;
            continue;
        }

        if (!alea_contains_point(sys, cell->root_node_id, x, y, z)) continue;

        if (cell->fill_universe > 0) {
            alea_matrix_t fill;
            if (volume_path_fill_matrix(sys, cell, &fill) != 0) return -1;

            double child_x = x;
            double child_y = y;
            double child_z = z;
            alea_matrix_transform_point_inverse(&fill, &child_x, &child_y, &child_z);

            volume_path_builder_t child = *builder;
            volume_path_push_ancestor(&child, cell_index, cell->universe_id);
            alea_matrix_t composed;
            alea_matrix_multiply(&composed, &builder->transform, &fill);
            if (!alea_matrix_invert(&composed)) return -1;
            child.transform = composed;

            if (volume_path_lookup_universe(sys, cell->fill_universe,
                                            child_x, child_y, child_z,
                                            depth + 1, &child, out, found) != 0) {
                return -1;
            }
            if (*found) return 0;
            continue;
        }

        volume_path_enum_ctx_t ctx = {.out = out, .max = 1, .count = 0};
        volume_path_emit(&ctx, sys, cell_index, depth, builder);
        *found = true;
        return 0;
    }
    return 0;
}

int alea_volume_path_at_point(alea_system_t* sys,
                              double x, double y, double z,
                              alea_volume_path_t* out_path) {
    if (!sys || !out_path) return -1;
    if (alea_volume_path_index_ensure(sys) != 0) return -1;
    if (alea_system_prepare_query_caches(sys, ALEA_CACHE_HIER_SPATIAL) != 0) return -1;

    volume_path_builder_t builder;
    memset(&builder, 0, sizeof(builder));
    alea_matrix_identity(&builder.transform);

    bool found = false;
    if (volume_path_lookup_universe(sys, 0, x, y, z, 0,
                                    &builder, out_path, &found) != 0) {
        return -1;
    }
    if (!found) return 0;

    const alea_volume_path_t* canonical =
        alea_volume_path_index_find(sys->volume_path_index, out_path);
    if (canonical) *out_path = *canonical;
    return 1;
}





int alea_tighten_cell_bbox(const alea_system_t* sys,
                                size_t cell_index,
                                double tol,
                                alea_bbox_t* out) {
    if (!sys || !out || tol <= 0.0) return -1;
    size_t n_cells = alea_vec_count(&sys->cells);
    if (cell_index >= n_cells) return -1;

    const alea_cell_entry_t* cell = &sys->cells.data[cell_index];
    if (cell->root_node_id == ALEA_NODE_ID_INVALID) return -1;

    const alea_bbox_t box_v = alea_node_bbox_get(&sys->nodes.data[cell->root_node_id].bbox);
    const alea_bbox_t* box = &box_v;
    if (box->min_x > box->max_x) return -1;

    alea_tighten_tree_bbox(sys, cell->root_node_id, box, tol, out);
    return 0;
}


int alea_tighten_cell_bbox_numerical(alea_system_t* sys, int cell_index) {
    if (!sys || cell_index < 0) return -1;
    size_t n_cells = alea_vec_count(&sys->cells);
    if ((size_t)cell_index >= n_cells) return -1;

    alea_cell_entry_t* cell = &sys->cells.data[cell_index];
    if (cell->root_node_id == ALEA_NODE_ID_INVALID) return -1;

    alea_bbox_t tight;
    if (alea_tighten_bbox_numerical(sys, cell->root_node_id, 1.0, &tight) != 0)
        return -1;

    alea_node_bbox_set(&sys->nodes.data[cell->root_node_id].bbox, &tight);
    return 0;
}

/* ============================================================================
 * CONFIGURATION
 * ============================================================================ */

/* ============================================================================
 * EXPORT
 *
 * MCNP export: mcnp_export() in alea_mcnp.h
 * OpenMC export: openmc_export() in alea_openmc.h
 * ============================================================================ */

/* ============================================================================
 * VOID GENERATION
 * ============================================================================ */

void_result_t* alea_void_generate_in_region(alea_system_t* sys,
                                            alea_node_id_t bounds_region) {
    if (!sys) return NULL;

    /* Build local octree config from sys->config */
    octree_config_t local_config = {
        .max_depth = sys->config.void_max_depth,
        .min_size = sys->config.void_min_size,
        .probes_per_axis = sys->config.void_probes_per_axis,
    };
    return alea_generate_void_in_region(sys, bounds_region, &local_config);
}

void_result_t* alea_void_generate_in_bbox(alea_system_t* sys,
                                          const alea_bbox_t* bounds) {
    if (!sys) return NULL;

    /* Build local octree config from sys->config */
    octree_config_t local_config = {
        .max_depth = sys->config.void_max_depth,
        .min_size = sys->config.void_min_size,
        .probes_per_axis = sys->config.void_probes_per_axis,
    };
    return alea_generate_void_in_bbox(sys, bounds, &local_config);
}

size_t alea_void_count(const void_result_t* result) {
    if (!result) return 0;
    return result->void_regions.count;
}

int alea_void_get(const void_result_t* result, size_t index, alea_bbox_t* box) {
    if (!result || !box) return -1;
    if (index >= result->void_regions.count) return -1;
    *box = result->void_regions.data[index].bbox;
    return 0;
}

void alea_void_free(void_result_t* result) {
    if (result) {
        alea_void_result_destroy(result);
    }
}

/* alea_void_to_node is defined in alea_void.c */

int alea_void_merge(alea_system_t* sys, void_result_t* result) {
    if (!sys || !result) return -1;

    /* Build local merge config from sys->config */
    alea_void_merge_config_t internal_config = {
        .cell_weight = sys->config.merge_cell_weight,
        .surface_weight = sys->config.merge_surface_weight,
        .max_surfaces_per_cell = sys->config.merge_max_surfaces,
        .min_cells = sys->config.merge_min_cells,
        .use_greedy = sys->config.merge_use_greedy,
        .consolidate_max_surfaces = sys->config.void_consolidate
    };

    return alea_merge_void_cells(sys, result, &internal_config);
}

/* ============================================================================
 * LOGGING
 * ============================================================================ */




/* ============================================================================
 * MATERIAL OPERATIONS
 * ============================================================================ */

int alea_add_material(alea_system_t* sys, int material_id) {
    if (!sys) return -1;

    int final_id;
    if (material_id <= 0) {
        final_id = sys->next_auto_material_id++;
    } else {
        /* Check for duplicate MCNP ID */
        if (alea_find_material_by_id(sys, material_id) >= 0) {
            ALEA_LOG_WARN("Material ID %d already exists, auto-assigning new ID %d",
                         material_id, sys->next_auto_material_id);
            final_id = sys->next_auto_material_id++;
        } else {
            final_id = material_id;
            if (material_id >= sys->next_auto_material_id) {
                sys->next_auto_material_id = material_id + 1;
            }
        }
    }

    int idx = (int)alea_vec_count(&sys->materials);
    alea_material_t* mat = alea_vec_push_uninit(&sys->materials, alea_material_t);
    if (!mat) return -1;

    memset(mat, 0, sizeof(*mat));
    mat->material_id = final_id;
    return idx;
}

int alea_find_material_by_id(const alea_system_t* sys, int material_id) {
    if (!sys) return -1;
    for (size_t i = 0; i < alea_vec_count(&sys->materials); i++) {
        if (sys->materials.data[i].material_id == material_id) {
            return (int)i;
        }
    }
    return -1;
}

int alea_create_mixture(alea_system_t* sys, const int* mat_ids,
                            const double* fractions, size_t count, int new_mat_id) {
    if (!sys || !mat_ids || !fractions || count == 0) return -1;

    /* Auto-assign ID if not specified */
    if (new_mat_id <= 0) {
        new_mat_id = 1;
        for (size_t i = 0; i < alea_vec_count(&sys->mixtures); i++) {
            if (sys->mixtures.data[i].mixture_id >= new_mat_id) {
                new_mat_id = sys->mixtures.data[i].mixture_id + 1;
            }
        }
    }

    /* Create mixture and add components */
    alea_mixture_t* mix = alea_mixture_create(new_mat_id);
    if (!mix) return -1;

    for (size_t i = 0; i < count; i++) {
        if (alea_mixture_add_component(mix, mat_ids[i], fractions[i]) < 0) {
            alea_mixture_destroy(mix);
            return -1;
        }
    }

    int result = alea_add_mixture(sys, mix);
    alea_mixture_destroy(mix);
    return (result >= 0) ? new_mat_id : -1;
}

size_t alea_material_count(const alea_system_t* sys) {
    return sys ? alea_vec_count(&sys->materials) : 0;
}

int alea_material_get_id(const alea_system_t* sys, int mat_index) {
    if (!sys || mat_index < 0 || (size_t)mat_index >= alea_vec_count(&sys->materials))
        return -1;
    return sys->materials.data[mat_index].material_id;
}

int alea_material_add_nuclide(alea_system_t* sys, int mat_index,
                              int zaid, const char* library, double fraction) {
    if (!sys || mat_index < 0 || (size_t)mat_index >= alea_vec_count(&sys->materials))
        return -1;
    return alea_mat_add_nuclide(&sys->materials.data[mat_index], zaid, library, fraction);
}

int alea_material_add_element(alea_system_t* sys, int mat_index,
                              int Z, const char* library, double fraction) {
    if (!sys || mat_index < 0 || (size_t)mat_index >= alea_vec_count(&sys->materials))
        return -1;
    return alea_mat_add_element(&sys->materials.data[mat_index], Z, library, fraction);
}

int alea_material_set_density(alea_system_t* sys, int mat_index, double density) {
    if (!sys || mat_index < 0 || (size_t)mat_index >= alea_vec_count(&sys->materials))
        return -1;
    alea_mat_set_density(&sys->materials.data[mat_index], density);
    return 0;
}

int alea_material_set_weight_fraction(alea_system_t* sys, int mat_index, bool is_weight) {
    if (!sys || mat_index < 0 || (size_t)mat_index >= alea_vec_count(&sys->materials))
        return -1;
    sys->materials.data[mat_index].is_weight_fraction = is_weight;
    return 0;
}

int alea_material_expand_elements(alea_system_t* sys, int mat_index) {
    if (!sys || mat_index < 0 || (size_t)mat_index >= alea_vec_count(&sys->materials))
        return -1;
    return alea_mat_expand_elements(&sys->materials.data[mat_index]);
}

size_t alea_material_nuclide_count(const alea_system_t* sys, int mat_index) {
    if (!sys || mat_index < 0 || (size_t)mat_index >= alea_vec_count(&sys->materials))
        return 0;
    return alea_vec_count(&sys->materials.data[mat_index].nuclides);
}

int alea_material_nuclide_get(const alea_system_t* sys, int mat_index,
                              size_t nuc_index, int* zaid,
                              const char** library, double* fraction) {
    if (!sys || mat_index < 0 || (size_t)mat_index >= alea_vec_count(&sys->materials))
        return -1;
    const alea_material_t* mat = &sys->materials.data[mat_index];
    if (nuc_index >= alea_vec_count(&mat->nuclides)) return -1;
    const alea_nuclide_t* nuc = &mat->nuclides.data[nuc_index];
    if (zaid) *zaid = nuc->zaid;
    if (library) *library = nuc->library;
    if (fraction) *fraction = nuc->fraction;
    return 0;
}

size_t alea_material_element_count(const alea_system_t* sys, int mat_index) {
    if (!sys || mat_index < 0 || (size_t)mat_index >= alea_vec_count(&sys->materials))
        return 0;
    return alea_vec_count(&sys->materials.data[mat_index].elements);
}

int alea_material_element_get(const alea_system_t* sys, int mat_index,
                              size_t elem_index, int* Z,
                              const char** library, double* fraction) {
    if (!sys || mat_index < 0 || (size_t)mat_index >= alea_vec_count(&sys->materials))
        return -1;
    const alea_material_t* mat = &sys->materials.data[mat_index];
    if (elem_index >= alea_vec_count(&mat->elements)) return -1;
    const alea_element_comp_t* elem = &mat->elements.data[elem_index];
    if (Z) *Z = elem->atomic_number;
    if (library) *library = elem->library;
    if (fraction) *fraction = elem->fraction;
    return 0;
}

int alea_material_get_density(const alea_system_t* sys, int mat_index,
                              double* density, bool* has_density) {
    if (!sys || mat_index < 0 || (size_t)mat_index >= alea_vec_count(&sys->materials))
        return -1;
    const alea_material_t* mat = &sys->materials.data[mat_index];
    if (density) *density = mat->standard_density;
    if (has_density) *has_density = mat->has_standard_density;
    return 0;
}

bool alea_material_is_weight_fraction(const alea_system_t* sys, int mat_index) {
    if (!sys || mat_index < 0 || (size_t)mat_index >= alea_vec_count(&sys->materials))
        return false;
    return sys->materials.data[mat_index].is_weight_fraction;
}

/* ============================================================================
 * MIXTURE QUERY OPERATIONS
 * ============================================================================ */

int alea_find_mixture_by_id(const alea_system_t* sys, int mixture_id) {
    if (!sys) return -1;
    for (size_t i = 0; i < alea_vec_count(&sys->mixtures); i++) {
        if (sys->mixtures.data[i].mixture_id == mixture_id) {
            return (int)i;
        }
    }
    return -1;
}

size_t alea_mixture_count(const alea_system_t* sys) {
    return sys ? alea_vec_count(&sys->mixtures) : 0;
}

int alea_mixture_get_id(const alea_system_t* sys, int mix_index) {
    if (!sys || mix_index < 0 || (size_t)mix_index >= alea_vec_count(&sys->mixtures))
        return -1;
    return sys->mixtures.data[mix_index].mixture_id;
}

size_t alea_mixture_component_count(const alea_system_t* sys, int mix_index) {
    if (!sys || mix_index < 0 || (size_t)mix_index >= alea_vec_count(&sys->mixtures))
        return 0;
    return alea_vec_count(&sys->mixtures.data[mix_index].components);
}

int alea_mixture_component_get(const alea_system_t* sys, int mix_index,
                               size_t comp_index, int* material_id,
                               double* fraction) {
    if (!sys || mix_index < 0 || (size_t)mix_index >= alea_vec_count(&sys->mixtures))
        return -1;
    const alea_mixture_t* mix = &sys->mixtures.data[mix_index];
    if (comp_index >= alea_vec_count(&mix->components)) return -1;
    const alea_mixture_comp_t* comp = &mix->components.data[comp_index];
    if (material_id) *material_id = comp->material_id;
    if (fraction) *fraction = comp->fraction;
    return 0;
}

/* ============================================================================
 * EXTRACT / FILTER
 * ============================================================================ */

int alea_carve_universe(alea_system_t* sys, int universe_id,
                        alea_node_id_t carve_root, int simplify, int cell_limit,
                        int* out_modified, int* out_removed) {
    if (out_modified) *out_modified = 0;
    if (out_removed) *out_removed = 0;
    if (!sys || carve_root == ALEA_NODE_ID_INVALID ||
        carve_root >= alea_vec_count(&sys->nodes)) {
        alea_set_error_detail(ALEA_ERR_INVALID_ARG,
            "alea_carve_universe: invalid system or carve region");
        return -1;
    }

    size_t scan_count = alea_vec_count(&sys->cells);
    if (cell_limit >= 0 && (size_t)cell_limit < scan_count)
        scan_count = (size_t)cell_limit;

    size_t candidate_count = 0;
    for (size_t i = 0; i < scan_count; i++) {
        const alea_cell_entry_t* cell = &sys->cells.data[i];
        if (cell->universe_id != universe_id ||
            cell->root_node_id == ALEA_NODE_ID_INVALID) continue;
        if (cell->lat_type != 0 && cell->lat_fill) {
            alea_set_error_detail(ALEA_ERR_INVALID_ARG,
                "alea_carve_universe: universe %d contains lattice cell %d; "
                "lattice masking is not supported", universe_id,
                cell->mc_cell_id);
            return -1;
        }
        candidate_count++;
    }

    if (candidate_count == 0) return 0;

    const alea_node_t* carve_node = &sys->nodes.data[carve_root];
    const alea_bbox_t carve_bbox = alea_node_bbox_get(&carve_node->bbox);

    size_t* indices = malloc(candidate_count * sizeof(*indices));
    alea_node_id_t* roots = malloc(candidate_count * sizeof(*roots));
    if (!indices || !roots) {
        free(indices);
        free(roots);
        alea_set_error_detail(ALEA_ERR_OUT_OF_MEMORY,
            "alea_carve_universe: failed to allocate rewrite table");
        return -1;
    }

    size_t n = 0;
    for (size_t i = 0; i < scan_count; i++) {
        const alea_cell_entry_t* cell = &sys->cells.data[i];
        if (cell->universe_id != universe_id ||
            cell->root_node_id == ALEA_NODE_ID_INVALID) continue;

        const alea_node_t* cell_node = &sys->nodes.data[cell->root_node_id];
        const alea_bbox_t cell_bbox = alea_node_bbox_get(&cell_node->bbox);
        if (!alea_bbox_intersects(&cell_bbox, &carve_bbox)) continue;

        alea_node_id_t root = alea_create_difference(sys, cell->root_node_id,
                                                      carve_root);
        if (root == ALEA_NODE_ID_INVALID) {
            free(indices);
            free(roots);
            alea_set_error_detail(ALEA_ERR_OUT_OF_MEMORY,
                "alea_carve_universe: failed to build difference for cell %d",
                cell->mc_cell_id);
            return -1;
        }
        if (simplify) {
            alea_error_clear();
            root = alea_tree_simplify(sys, root, NULL);
            if (root == ALEA_NODE_ID_INVALID && alea_error_code() != ALEA_OK) {
                free(indices);
                free(roots);
                return -1;
            }
        }
        indices[n] = i;
        roots[n] = root;
        n++;
    }

    int modified = 0;
    for (size_t i = 0; i < n; i++) {
        alea_cell_entry_t* cell = &sys->cells.data[indices[i]];
        cell->root_node_id = roots[i];
        /* A modified region must never export an old pre-TRCL expression. */
        cell->original_root_node_id = ALEA_NODE_ID_INVALID;
        modified++;
    }

    int removed = 0;
    /* Remove from highest index so earlier recorded indices stay valid. */
    for (size_t i = n; i > 0; i--) {
        if (roots[i - 1] == ALEA_NODE_ID_INVALID) {
            if (alea_cell_remove(sys, (int)indices[i - 1]) < 0) {
                free(indices);
                free(roots);
                return -1;
            }
            removed++;
        }
    }

    free(indices);
    free(roots);
    alea_system_invalidate_query_caches(sys, ALEA_CACHE_ALL);

    if (out_modified) *out_modified = modified;
    if (out_removed) *out_removed = removed;
    return 0;
}

alea_system_t* alea_extract_region(const alea_system_t* sys, const alea_bbox_t* bbox) {
    if (!sys || !bbox) return NULL;

    /* Region-restricted flatten: descend the universe tree from the root,
     * propagate transforms, and only emit cells whose world-frame bbox
     * overlaps the query region. This prunes entire fill sub-universes
     * (and out-of-range lattice elements) at the parent level, so cost
     * is proportional to geometry inside the region — not to the total
     * cell count. Necessary for large lattice models where a full flatten
     * would not fit in memory. */
    alea_flatten_config_t cfg = ALEA_FLATTEN_DEFAULT;
    cfg.starting_universe_id = 0;
    cfg.clip_active = true;
    cfg.clip_bbox = *bbox;

    /* The source isn't semantically mutated; alea_flatten_to_new_system may
     * lazily build the universe index, which is the usual reason for the
     * const cast on read-only entry points. */
    return alea_flatten_to_new_system((alea_system_t*)sys, &cfg);
}

size_t alea_get_cells_in_bbox(const alea_system_t* sys, const alea_bbox_t* bbox,
                                  int* out_indices, size_t max_count) {
    if (!sys || !bbox) return 0;

    size_t count = 0;
    for (size_t i = 0; i < alea_vec_count(&sys->cells); i++) {
        const alea_cell_entry_t* cell = &sys->cells.data[i];
        /* Check if cell bbox intersects query bbox */
        if (cell->root_node_id == ALEA_NODE_ID_INVALID) continue;

        const alea_bbox_t cb_v = alea_node_bbox_get(&sys->nodes.data[cell->root_node_id].bbox);
        const alea_bbox_t* cb = &cb_v;
        bool intersects = (cb->min_x <= bbox->max_x && cb->max_x >= bbox->min_x &&
                          cb->min_y <= bbox->max_y && cb->max_y >= bbox->min_y &&
                          cb->min_z <= bbox->max_z && cb->max_z >= bbox->min_z);
        if (intersects) {
            if (out_indices && count < max_count) {
                out_indices[count] = (int)i;
            }
            count++;
        }
    }
    return count;
}

/* ============================================================================
 * CELL INFO
 * ============================================================================ */

int alea_cell_get(const alea_system_t* sys, size_t index,
                      int* cell_id, int* material_id, double* density,
                      int* universe_id, int* fill_universe, alea_node_id_t* root) {
    if (!sys) return -1;
    if (index >= alea_vec_count(&sys->cells)) return -1;

    const alea_cell_entry_t* c = &sys->cells.data[index];
    if (cell_id) *cell_id = c->mc_cell_id;
    if (material_id) *material_id = c->material_id;
    if (density) *density = c->density;
    if (universe_id) *universe_id = c->universe_id;
    if (fill_universe) *fill_universe = c->fill_universe;
    if (root) *root = c->root_node_id;
    return 0;
}

int alea_cell_find(const alea_system_t* sys, int cell_id) {
    if (!sys) return -1;
    return alea_find_cell_by_id(sys, cell_id);
}

int alea_cell_get_info(const alea_system_t* sys, size_t index, alea_cell_info_t* info) {
    if (!sys || !info) return -1;
    if (index >= alea_vec_count(&sys->cells)) return -1;

    const alea_cell_entry_t* c = &sys->cells.data[index];
    info->cell_id = c->mc_cell_id;
    info->material_id = c->material_id;
    info->density = c->density;
    info->is_mass_density = c->is_mass_density;
    info->universe_id = c->universe_id;
    info->fill_universe = c->fill_universe;
    info->fill_transform = c->fill_transform;
    info->temperature = c->temperature;
    info->has_temperature = c->has_temperature;
    info->root = c->root_node_id;

    /* Get bbox from root node if available */
    if (c->root_node_id != ALEA_NODE_ID_INVALID &&
        c->root_node_id < alea_vec_count(&sys->nodes)) {
        info->bbox = alea_node_bbox_get(&sys->nodes.data[c->root_node_id].bbox);
    } else {
        info->bbox = (alea_bbox_t){-1e30, 1e30, -1e30, 1e30, -1e30, 1e30};
    }

    /* Lattice info */
    info->lat_type = c->lat_type;
    memcpy(info->lat_fill_dims, c->lat_fill_dims, sizeof(info->lat_fill_dims));
    info->lat_fill = c->lat_fill;
    info->lat_fill_count = c->lat_fill_count;
    info->lat_fill_repeating = (int)c->lat_fill_repeating;
    info->lat_fill_zero_element_coords = (int)c->lat_fill_zero_element_coords;
    memcpy(info->lat_pitch, c->lat_pitch, sizeof(info->lat_pitch));
    memcpy(info->lat_lower_left, c->lat_lower_left, sizeof(info->lat_lower_left));

    /* Comments */
    info->comments = c->comments;
    info->inline_comment = c->inline_comment;
    return 0;
}

int alea_cell_find_info(const alea_system_t* sys, int cell_id, alea_cell_info_t* info) {
    int index = alea_cell_find(sys, cell_id);
    if (index < 0) return -1;
    return alea_cell_get_info(sys, (size_t)index, info);
}

int alea_cells_in_universe(const alea_system_t* sys, int universe_id,
                                int* out_indices, size_t max_count) {
    if (!sys) return 0;

    size_t found = 0;
    for (size_t i = 0; i < alea_vec_count(&sys->cells); i++) {
        if (sys->cells.data[i].universe_id == universe_id) {
            if (out_indices && found < max_count) {
                out_indices[found] = (int)i;
            }
            found++;
        }
    }
    return (int)found;
}

/* ============================================================================
 * SURFACE INFO
 * ============================================================================ */

int alea_surface_get(const alea_system_t* sys, size_t index,
                          int* surface_id, alea_primitive_type_t* type,
                          alea_node_id_t* pos_node, alea_node_id_t* neg_node,
                          alea_boundary_type_t* boundary_type) {
    if (!sys) return -1;
    if (index >= alea_vec_count(&sys->surfaces)) return -1;

    const alea_surface_entry_t* s = &sys->surfaces.data[index];
    if (surface_id) *surface_id = s->mc_surface_id;
    if (pos_node) *pos_node = s->pos_node;
    if (neg_node) *neg_node = s->neg_node;
    if (boundary_type) *boundary_type = s->boundary_type;
    if (type && s->primitive_id < alea_vec_count(&sys->primitives)) {
        *type = sys->primitives.data[s->primitive_id].type;
    }
    return 0;
}

int alea_surface_find(const alea_system_t* sys, int surface_id) {
    if (!sys) return -1;
    /* Fast path: O(1) direct-address table (built after surface conversion) */
    if (surface_id > 0 && sys->mc_id_to_surface &&
        (size_t)surface_id < sys->mc_id_to_surface_size) {
        uint32_t idx = sys->mc_id_to_surface[surface_id];
        return (idx != UINT32_MAX) ? (int)idx : -1;
    }
    /* Fallback: linear scan (table not yet built or ID out of range) */
    for (size_t i = 0; i < alea_vec_count(&sys->surfaces); i++) {
        if (sys->surfaces.data[i].mc_surface_id == surface_id)
            return (int)i;
    }
    return -1;
}

int alea_surface_id_at(const alea_system_t* sys, size_t idx) {
    if (!sys) return -1;
    if (idx >= alea_vec_count(&sys->surfaces)) return -1;
    return sys->surfaces.data[idx].mc_surface_id;
}

size_t alea_get_surface_ids(const alea_system_t* sys, int* out_ids) {
    if (!sys || !out_ids) return 0;
    size_t n = alea_vec_count(&sys->surfaces);
    for (size_t i = 0; i < n; i++) {
        out_ids[i] = sys->surfaces.data[i].mc_surface_id;
    }
    return n;
}

/* ============================================================================
 * UNIVERSE INFO
 * ============================================================================ */

int alea_universe_get(const alea_system_t* sys, size_t index,
                           int* universe_id, size_t* cell_count, alea_bbox_t* bbox) {
    if (!sys) return -1;
    if (index >= alea_vec_count(&sys->universes)) return -1;

    const alea_universe_t* u = &sys->universes.data[index];
    if (universe_id) *universe_id = u->universe_id;
    if (cell_count) *cell_count = u->cell_indices.count;
    if (bbox) *bbox = u->bbox;
    return 0;
}

int alea_universe_find(const alea_system_t* sys, int universe_id) {
    if (!sys) return -1;
    /* O(1) lookup via hashmap */
    int* idx_ptr = universe_hashmap_get(&sys->universe_index, universe_id);
    if (idx_ptr) return *idx_ptr;
    /* Fallback: linear scan (hashmap may not be populated yet) */
    for (size_t i = 0; i < alea_vec_count(&sys->universes); i++) {
        if (sys->universes.data[i].universe_id == universe_id) {
            return (int)i;
        }
    }
    return -1;
}

/* ============================================================================
 * EXTENDED GEOMETRY QUERIES
 * ============================================================================ */

int alea_find_cell_at(alea_system_t* sys, double x, double y, double z,
                           int* out_cell_id, int* out_material) {
    if (!sys) return -1;

    alea_cell_hit_t hit;
    if (alea_find_deepest_cell_hit_at_point(sys, x, y, z, &hit) != 0) {
        return -1;
    }

    if (out_cell_id) *out_cell_id = hit.cell_id;
    if (out_material) *out_material = hit.material_id;
    return 0;
}

void alea_set_debug_trace(int enable) {
    alea_set_debug_point_trace(enable);
}

/* ============================================================================
 * CSG TREE CONSTRUCTION - EXTENDED
 * ============================================================================ */

int alea_get_cell_id(const alea_system_t* sys, int cell_index) {
    if (!sys) return -1;
    if (cell_index < 0 || (size_t)cell_index >= alea_vec_count(&sys->cells)) return -1;
    return sys->cells.data[cell_index].mc_cell_id;
}

/* ============================================================================
 * RENUMBERING
 * ============================================================================ */

int alea_renumber_cells(alea_system_t* sys, int start_id) {
    if (!sys) return -1;
    for (size_t i = 0; i < alea_vec_count(&sys->cells); i++) {
        sys->cells.data[i].mc_cell_id = start_id + (int)i;
    }
    sys->next_auto_cell_id = start_id + (int)alea_vec_count(&sys->cells);

    /* Rebuild cell hash map */
    cell_hashmap_clear(&sys->cell_index);
    for (size_t i = 0; i < alea_vec_count(&sys->cells); i++) {
        cell_hashmap_put(&sys->cell_index, sys->cells.data[i].mc_cell_id, (int)i);
    }
    return 0;
}

int alea_renumber_surfaces(alea_system_t* sys, int start_id) {
    if (!sys) return -1;
    for (size_t i = 0; i < alea_vec_count(&sys->surfaces); i++) {
        sys->surfaces.data[i].mc_surface_id = start_id + (int)i;
    }
    sys->next_auto_surface_id = start_id + (int)alea_vec_count(&sys->surfaces);
    return 0;
}

int alea_offset_cell_ids(alea_system_t* sys, int offset) {
    if (!sys) return -1;
    for (size_t i = 0; i < alea_vec_count(&sys->cells); i++) {
        sys->cells.data[i].mc_cell_id += offset;
    }
    sys->next_auto_cell_id += offset;

    /* Rebuild cell hash map */
    cell_hashmap_clear(&sys->cell_index);
    for (size_t i = 0; i < alea_vec_count(&sys->cells); i++) {
        cell_hashmap_put(&sys->cell_index, sys->cells.data[i].mc_cell_id, (int)i);
    }
    return 0;
}

int alea_offset_surface_ids(alea_system_t* sys, int offset) {
    if (!sys) return -1;
    for (size_t i = 0; i < alea_vec_count(&sys->surfaces); i++) {
        sys->surfaces.data[i].mc_surface_id += offset;
    }
    sys->next_auto_surface_id += offset;
    return 0;
}

int alea_offset_material_ids(alea_system_t* sys, int offset) {
    if (!sys) return -1;
    for (size_t i = 0; i < alea_vec_count(&sys->cells); i++) {
        if (sys->cells.data[i].material_id != 0) {
            sys->cells.data[i].material_id += offset;
        }
    }
    for (size_t i = 0; i < alea_vec_count(&sys->materials); i++) {
        if (sys->materials.data[i].material_id != 0) {
            sys->materials.data[i].material_id += offset;
        }
    }
    return 0;
}

/* ============================================================================
 * FILTER OPERATIONS
 * ============================================================================ */

size_t alea_get_cells_by_material(const alea_system_t* sys, int material_id,
                                       int* out_indices, size_t max_count) {
    if (!sys) return 0;
    size_t count = 0;
    for (size_t i = 0; i < alea_vec_count(&sys->cells); i++) {
        if (sys->cells.data[i].material_id == material_id) {
            if (out_indices && count < max_count) {
                out_indices[count] = (int)i;
            }
            count++;
        }
    }
    return count;
}

size_t alea_get_cells_by_universe(const alea_system_t* sys, int universe_id,
                                       int* out_indices, size_t max_count) {
    if (!sys) return 0;
    size_t count = 0;
    for (size_t i = 0; i < alea_vec_count(&sys->cells); i++) {
        if (sys->cells.data[i].universe_id == universe_id) {
            if (out_indices && count < max_count) {
                out_indices[count] = (int)i;
            }
            count++;
        }
    }
    return count;
}

size_t alea_get_cells_filling_universe(const alea_system_t* sys, int universe_id,
                                            int* out_indices, size_t max_count) {
    if (!sys) return 0;
    size_t count = 0;
    for (size_t i = 0; i < alea_vec_count(&sys->cells); i++) {
        if (sys->cells.data[i].fill_universe == universe_id) {
            if (out_indices && count < max_count) {
                out_indices[count] = (int)i;
            }
            count++;
        }
    }
    return count;
}

/* ============================================================================
 * VALIDATION
 * ============================================================================ */

int alea_validate(const alea_system_t* sys) {
    if (!sys) return 1;

    int issues = 0;
    for (size_t i = 0; i < alea_vec_count(&sys->nodes); i++) {
        const alea_node_t* n = &sys->nodes.data[i];
        alea_operation_t op = ALEA_GET_OPERATION(n);

        if (op != ALEA_OP_PRIMITIVE) {
            if (n->operation.left >= alea_vec_count(&sys->nodes)) issues++;
            if (op != ALEA_OP_COMPLEMENT && n->operation.right >= alea_vec_count(&sys->nodes)) issues++;
        }
    }

    for (size_t i = 0; i < alea_vec_count(&sys->cells); i++) {
        if (sys->cells.data[i].root_node_id >= alea_vec_count(&sys->nodes)) issues++;
    }

    return issues;
}

/* ============================================================================
 * MACROBODY EXPANSION
 * ============================================================================ */

#include "core/alea_macrobody.h"




int alea_expand_macrobodies_in_system(alea_system_t* sys) {
    if (!sys) return -1;

    int count = 0;
    for (size_t i = 0; i < alea_vec_count(&sys->cells); i++) {
        alea_cell_entry_t* cell = &sys->cells.data[i];
        if (cell->root_node_id == ALEA_NODE_ID_INVALID) continue;

        alea_node_id_t new_root = alea_expand_all_macrobodies(sys, cell->root_node_id);
        if (new_root != cell->root_node_id) {
            cell->root_node_id = new_root;
            count++;
        }
    }
    return count;
}

/* ============================================================================
 * CSG NODE INSPECTION
 * ============================================================================ */

alea_operation_t alea_node_operation(const alea_system_t* sys, alea_node_id_t node) {
    if (!sys || node >= alea_vec_count(&sys->nodes)) return ALEA_OP_PRIMITIVE;
    return ALEA_GET_OPERATION(&sys->nodes.data[node]);
}

alea_node_id_t alea_node_left(const alea_system_t* sys, alea_node_id_t node) {
    if (!sys || node >= alea_vec_count(&sys->nodes)) return ALEA_NODE_ID_INVALID;
    const alea_node_t* n = &sys->nodes.data[node];
    if (ALEA_GET_OPERATION(n) == ALEA_OP_PRIMITIVE) return ALEA_NODE_ID_INVALID;
    return n->operation.left;
}

alea_node_id_t alea_node_right(const alea_system_t* sys, alea_node_id_t node) {
    if (!sys || node >= alea_vec_count(&sys->nodes)) return ALEA_NODE_ID_INVALID;
    const alea_node_t* n = &sys->nodes.data[node];
    alea_operation_t op = ALEA_GET_OPERATION(n);
    if (op == ALEA_OP_PRIMITIVE || op == ALEA_OP_COMPLEMENT) return ALEA_NODE_ID_INVALID;
    return n->operation.right;
}

alea_primitive_type_t alea_node_primitive_type(const alea_system_t* sys, alea_node_id_t node) {
    if (!sys || node >= alea_vec_count(&sys->nodes)) return 0;
    const alea_node_t* n = &sys->nodes.data[node];
    if (ALEA_GET_OPERATION(n) != ALEA_OP_PRIMITIVE) return 0;
    return n->primitive.prim_type;
}

alea_primitive_id_t alea_node_primitive_id(const alea_system_t* sys, alea_node_id_t node) {
    if (!sys || node >= alea_vec_count(&sys->nodes)) return ALEA_PRIMITIVE_ID_INVALID;
    const alea_node_t* n = &sys->nodes.data[node];
    if (ALEA_GET_OPERATION(n) != ALEA_OP_PRIMITIVE) return ALEA_PRIMITIVE_ID_INVALID;
    return n->primitive.primitive_id;
}

int alea_node_primitive_data(const alea_system_t* sys, alea_node_id_t node,
                                 alea_primitive_data_t* out) {
    if (!sys || !out || node >= alea_vec_count(&sys->nodes)) return -1;
    const alea_node_t* n = &sys->nodes.data[node];
    if (ALEA_GET_OPERATION(n) != ALEA_OP_PRIMITIVE) return -1;

    uint32_t prim_id = n->primitive.primitive_id;
    if (prim_id >= alea_vec_count(&sys->primitives)) return -1;
    return alea_primitive_copy_data(sys, prim_id, out) ? 0 : -1;
}

int alea_node_sense(const alea_system_t* sys, alea_node_id_t node) {
    if (!sys || node >= alea_vec_count(&sys->nodes)) return 0;
    const alea_node_t* n = &sys->nodes.data[node];
    if (ALEA_GET_OPERATION(n) != ALEA_OP_PRIMITIVE) return 0;
    return n->primitive.sense;
}

int alea_node_surface_id(const alea_system_t* sys, alea_node_id_t node) {
    if (!sys || node >= alea_vec_count(&sys->nodes)) return 0;
    const alea_node_t* n = &sys->nodes.data[node];
    if (ALEA_GET_OPERATION(n) != ALEA_OP_PRIMITIVE) return 0;
    return n->primitive.mc_surface_id;
}

/* ============================================================================
 * CSG EXPRESSION STRINGIFIER
 * ============================================================================ */

typedef enum { CTX_TOP, CTX_UNION, CTX_INTERSECTION } expr_ctx_t;

static bool cell_expr_recursive(const alea_system_t* sys, uint32_t node_id,
                                const char* union_op, const char* inter_op,
                                const char* compl_op,
                                str_builder_t* sb, expr_ctx_t ctx) {
    if (node_id == ALEA_NODE_ID_INVALID) return false;
    if (node_id >= alea_vec_count(&sys->nodes)) return false;

    const alea_node_t* node = &sys->nodes.data[node_id];
    alea_operation_t op = ALEA_GET_OPERATION(node);

    if (op == ALEA_OP_PRIMITIVE) {
        int surface_id = node->primitive.mc_surface_id;
        if (node->primitive.sense < 0)
            surface_id = -surface_id;
        str_builder_int(sb, surface_id);
        return true;
    }

    if (op == ALEA_OP_COMPLEMENT) {
        uint32_t left = node->operation.left;
        if (left >= alea_vec_count(&sys->nodes)) return false;

        const alea_node_t* child = &sys->nodes.data[left];
        alea_operation_t child_op = ALEA_GET_OPERATION(child);

        /* Complement of a primitive: just flip the sign */
        if (child_op == ALEA_OP_PRIMITIVE) {
            int surface_id = child->primitive.mc_surface_id;
            /* Complement flips: positive sense -> negative surface_id, and vice versa */
            if (child->primitive.sense > 0)
                surface_id = -surface_id;
            str_builder_int(sb, surface_id);
            return true;
        }

        /* Complex expression: compl_op(inner) */
        str_builder_puts(sb, compl_op);
        str_builder_putc(sb, '(');
        if (!cell_expr_recursive(sys, left, union_op, inter_op, compl_op, sb, CTX_TOP))
            return false;
        str_builder_putc(sb, ')');
        return true;
    }

    if (op == ALEA_OP_UNION) {
        uint32_t left = node->operation.left;
        uint32_t right = node->operation.right;

        bool need_parens = (ctx == CTX_INTERSECTION);
        if (need_parens) str_builder_putc(sb, '(');

        if (!cell_expr_recursive(sys, left, union_op, inter_op, compl_op, sb, CTX_UNION))
            return false;
        str_builder_puts(sb, union_op);
        if (!cell_expr_recursive(sys, right, union_op, inter_op, compl_op, sb, CTX_UNION))
            return false;

        if (need_parens) str_builder_putc(sb, ')');
        return true;
    }

    if (op == ALEA_OP_INTERSECTION) {
        uint32_t left = node->operation.left;
        uint32_t right = node->operation.right;

        bool need_parens = (ctx == CTX_UNION);
        if (need_parens) str_builder_putc(sb, '(');

        if (!cell_expr_recursive(sys, left, union_op, inter_op, compl_op, sb, CTX_INTERSECTION))
            return false;
        str_builder_puts(sb, inter_op);
        if (!cell_expr_recursive(sys, right, union_op, inter_op, compl_op, sb, CTX_INTERSECTION))
            return false;

        if (need_parens) str_builder_putc(sb, ')');
        return true;
    }

    if (op == ALEA_OP_DIFFERENCE) {
        uint32_t left = node->operation.left;
        uint32_t right = node->operation.right;

        /* A - B = A inter compl(B) */
        if (!cell_expr_recursive(sys, left, union_op, inter_op, compl_op, sb, CTX_INTERSECTION))
            return false;

        str_builder_puts(sb, inter_op);
        str_builder_puts(sb, compl_op);

        /* Always parenthesize the right side unless it's a primitive */
        if (right < alea_vec_count(&sys->nodes)) {
            const alea_node_t* rn = &sys->nodes.data[right];
            alea_operation_t rop = ALEA_GET_OPERATION(rn);
            if (rop != ALEA_OP_PRIMITIVE) {
                str_builder_putc(sb, '(');
                if (!cell_expr_recursive(sys, right, union_op, inter_op, compl_op, sb, CTX_TOP))
                    return false;
                str_builder_putc(sb, ')');
            } else {
                if (!cell_expr_recursive(sys, right, union_op, inter_op, compl_op, sb, CTX_TOP))
                    return false;
            }
        } else {
            return false;
        }
        return true;
    }

    return false;
}

char* alea_cell_expr(const alea_system_t* sys, size_t cell_index,
                     const char* union_op, const char* inter_op,
                     const char* compl_op) {
    if (!sys || !union_op || !inter_op || !compl_op) return NULL;
    if (cell_index >= alea_vec_count(&sys->cells)) return NULL;

    const alea_cell_entry_t* cell = &sys->cells.data[cell_index];
    if (cell->root_node_id == ALEA_NODE_ID_INVALID) return NULL;

    arena_t arena;
    if (!arena_init_with_size(&arena, 4096)) return NULL;

    str_builder_t sb;
    str_builder_init(&sb, &arena, 256);

    bool ok = cell_expr_recursive(sys, cell->root_node_id,
                                  union_op, inter_op, compl_op,
                                  &sb, CTX_TOP);

    char* result = NULL;
    if (ok && !str_builder_error(&sb)) {
        str_builder_finish(&sb);
        const char* s = str_builder_get(&sb);
        if (s) result = alea_strdup(s);
    }

    arena_free(&arena);
    return result;
}
