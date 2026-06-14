// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

#ifndef ALEA_SYSTEM_H
#define ALEA_SYSTEM_H

#include <stddef.h>
#include <stdint.h>
#include <stdbool.h>
#include <signal.h>
#include <math.h>    /* nextafterf for alea_node_bbox_set */
#include <float.h>   /* FLT_MAX for empty-bbox sentinel */
#include "util/alea_atomic.h"
#include "alea_types.h"  // All primitive and node types
#include "alea_materials.h"  // Material, element, mixture types
#include "util/alea_vec.h"  // Dynamic array types
#include "util/alea_hashmap.h"  // Generic hash map
#include "core/alea_cell_complement.h"
#include "core/alea_cell.h"
#include "core/alea_surface.h"
#include "core/alea_transform.h"
#include "core/alea_universe.h"

/* Global interrupt flag (defined in alea_error.c) */
extern volatile sig_atomic_t g_alea_interrupted;

#define ALEA_CACHE_UNIVERSE      (1u << 0)
#define ALEA_CACHE_CELL_SURFACES (1u << 1)
#define ALEA_CACHE_SURFACE_BVH   (1u << 3)
#define ALEA_CACHE_ADJACENCY     (1u << 4)
#define ALEA_CACHE_HIER_SPATIAL  (1u << 5)
/* The hierarchical spatial index is the only spatial backend. */
#define ALEA_CACHE_RAYCAST \
    (ALEA_CACHE_UNIVERSE | ALEA_CACHE_CELL_SURFACES | \
     ALEA_CACHE_HIER_SPATIAL | ALEA_CACHE_SURFACE_BVH | \
     ALEA_CACHE_ADJACENCY)
#define ALEA_CACHE_RAYCAST_HIER ALEA_CACHE_RAYCAST
#define ALEA_CACHE_ALL           0xFFFFFFFFu

/**
 * @brief Check interrupt flag and return early if set
 *
 * Use in long-running loops. Sets error detail and returns the given value.
 */
#define ALEA_CHECK_INTERRUPTED(retval) do { \
    if (g_alea_interrupted) { \
        alea_set_error_detail(ALEA_ERR_INTERRUPTED, "Operation interrupted"); \
        return (retval); \
    } \
} while (0)

/* Forward declaration for BVH */
struct alea_bvh;
struct alea_hier_spatial_index;
struct alea_volume_path_index;

/* ============================================================================
 * VECTOR TYPES FOR CSG SYSTEM ARRAYS
 * ============================================================================ */

ALEA_VEC_DEFINE(alea_mixture_vec, alea_mixture_t);
ALEA_VEC_DEFINE(alea_cell_ref_vec, alea_cell_ref_t);
ALEA_VEC_DEFINE(alea_material_vec, alea_material_t);


/**
 * @file alea_system.h
 * @brief Core CSG system - Single flattened structure
 *
 * This is the main public API for the CSG system. It provides a single
 * alea_system_t structure that manages all primitives, nodes, and mappings
 * with automatic growth and deduplication.
 *
 * Key features:
 * - Automatic array growth (no manual sizing)
 * - Automatic primitive deduplication (fuzzy matching)
 * - Direct array access (no getters/setters)
 * - Built-in statistics tracking
 */

// Forward declarations (from MCNP parser - internal)
typedef struct mcnp_context mcnp_context_t;

// ============================================================================
// PRIMITIVE STORAGE
// ============================================================================

/**
 * @brief Entry for a stored primitive
 *
 * Primitives are deduplicated - identical primitives (within tolerance)
 * share the same entry and primitive_id.
 */
typedef struct {
    alea_primitive_type_t type;
    uint32_t payload_index;
    uint32_t ref_count;  // How many nodes reference this primitive
} alea_primitive_entry_t;

ALEA_VEC_DEFINE(alea_primitive_vec, alea_primitive_entry_t);
ALEA_VEC_DEFINE(alea_plane_data_vec, alea_plane_data_t);
ALEA_VEC_DEFINE(alea_sphere_data_vec, alea_sphere_data_t);
ALEA_VEC_DEFINE(alea_cylinder_x_data_vec, alea_cylinder_x_data_t);
ALEA_VEC_DEFINE(alea_cylinder_y_data_vec, alea_cylinder_y_data_t);
ALEA_VEC_DEFINE(alea_cylinder_z_data_vec, alea_cylinder_z_data_t);
ALEA_VEC_DEFINE(alea_cone_x_data_vec, alea_cone_x_data_t);
ALEA_VEC_DEFINE(alea_cone_y_data_vec, alea_cone_y_data_t);
ALEA_VEC_DEFINE(alea_cone_z_data_vec, alea_cone_z_data_t);
ALEA_VEC_DEFINE(alea_box_data_vec, alea_box_data_t);
ALEA_VEC_DEFINE(alea_quadric_data_vec, alea_quadric_data_t);
ALEA_VEC_DEFINE(alea_torus_data_vec, alea_torus_data_t);
ALEA_VEC_DEFINE(alea_rcc_data_vec, alea_rcc_data_t);
ALEA_VEC_DEFINE(alea_box_general_data_vec, alea_box_general_data_t);
ALEA_VEC_DEFINE(alea_sph_data_vec, alea_sph_data_t);
ALEA_VEC_DEFINE(alea_trc_data_vec, alea_trc_data_t);
ALEA_VEC_DEFINE(alea_ell_data_vec, alea_ell_data_t);
ALEA_VEC_DEFINE(alea_rec_data_vec, alea_rec_data_t);
ALEA_VEC_DEFINE(alea_wed_data_vec, alea_wed_data_t);
ALEA_VEC_DEFINE(alea_rhp_data_vec, alea_rhp_data_t);
ALEA_VEC_DEFINE(alea_arb_data_vec, alea_arb_data_t);

// ============================================================================
// CSG TREE NODES
// ============================================================================

// Macros for accessing packed operation field
#define ALEA_GET_OPERATION(node) ((alea_operation_t)((node)->type_and_flags & 0xFF))
#define ALEA_SET_OPERATION(node, op) ((node)->type_and_flags = ((node)->type_and_flags & 0xFFFFFF00) | (op))



/**
 * @brief Compact float bounding box stored per CSG node.
 *
 * The per-node bbox is only a conservative spatial cull cache: containment is
 * always confirmed by full primitive evaluation, so float storage cannot affect
 * correctness, only how many subtrees are visited. Storing 6 floats instead of
 * 6 doubles halves the per-node bbox footprint (48B -> 24B), which dominates the
 * node array for large models. Field names mirror alea_bbox_t so existing
 * `node->bbox.min_x` reads still compile (float -> double promotion).
 *
 * Conversions go through alea_node_bbox_set/get which round outward (set) and
 * widen (get), guaranteeing the float box always encloses the true double box.
 */
typedef struct alea_node_bbox {
    float min_x, max_x;
    float min_y, max_y;
    float min_z, max_z;
} alea_node_bbox_t;

/* Point-in-box test directly on the float bbox (hot path; avoids a double copy). */
static inline bool alea_node_bbox_contains_point(const alea_node_bbox_t* b,
                                                 double x, double y, double z) {
    return x >= (double)b->min_x && x <= (double)b->max_x &&
           y >= (double)b->min_y && y <= (double)b->max_y &&
           z >= (double)b->min_z && z <= (double)b->max_z;
}

/* Widen a stored float bbox back to a double bbox (exact, enclosing). */
static inline alea_bbox_t alea_node_bbox_get(const alea_node_bbox_t* b) {
    alea_bbox_t r;
    r.min_x = (double)b->min_x; r.max_x = (double)b->max_x;
    r.min_y = (double)b->min_y; r.max_y = (double)b->max_y;
    r.min_z = (double)b->min_z; r.max_z = (double)b->max_z;
    return r;
}

/* Store a double bbox into float, rounding each bound OUTWARD so the float box
 * always contains the double box. Preserves the empty (min>max) sentinel. */
static inline void alea_node_bbox_set(alea_node_bbox_t* dst, const alea_bbox_t* s) {
    if (!(s->min_x <= s->max_x && s->min_y <= s->max_y && s->min_z <= s->max_z)) {
        dst->min_x = dst->min_y = dst->min_z =  FLT_MAX;
        dst->max_x = dst->max_y = dst->max_z = -FLT_MAX;
        return;
    }
    dst->min_x = nextafterf((float)s->min_x, -INFINITY);
    dst->min_y = nextafterf((float)s->min_y, -INFINITY);
    dst->min_z = nextafterf((float)s->min_z, -INFINITY);
    dst->max_x = nextafterf((float)s->max_x,  INFINITY);
    dst->max_y = nextafterf((float)s->max_y,  INFINITY);
    dst->max_z = nextafterf((float)s->max_z,  INFINITY);
}

/**
 * @brief CSG tree node
 *
 * Nodes are stored in a flat array. Child nodes are referenced by ID,
 * not by pointer, allowing the array to grow without invalidating references.
 *
 * The node uses a packed type_and_flags field to save space:
 * - Lower 8 bits: operation type
 * - Upper 24 bits: flags (for future use)
 */
typedef struct alea_node {
    // Node type and flags (packed)
    uint32_t type_and_flags;  // Operation type (8 bits) + flags (24 bits)

    // Material ID
    alea_material_id_t material_id;

    // Bounding box (compact float cull cache; see alea_node_bbox_t)
    alea_node_bbox_t bbox;

    // Union: either primitive data or operation children
    union {
        // For primitive nodes (ALEA_OP_PRIMITIVE)
        struct {
            alea_primitive_id_t primitive_id;
            alea_primitive_type_t prim_type;
            int8_t sense;               /* +1 outside, -1 inside */
            int8_t inverted;            /* 1 if normal was flipped during canonicalization */
            int mc_surface_id;        /* MCNP surface ID for export */
        } primitive;

        // For operation nodes (union/intersection/difference)
        // NOTE: Uses node IDs instead of pointers to allow pool reallocation
        struct {
            alea_node_id_t left;
            alea_node_id_t right;
        } operation;
    };
} alea_node_t;

ALEA_VEC_DEFINE(alea_node_vec, alea_node_t);

// ============================================================================
// HASH TABLE FOR PRIMITIVE DEDUPLICATION
// ============================================================================

#define PRIMITIVE_HASH_INITIAL_SIZE 256

/**
 * @brief Hash table entry for primitive lookup
 */
typedef struct primitive_hash_entry {
    uint32_t primitive_id;
    uint64_t hash;
    struct primitive_hash_entry* next;
} primitive_hash_entry_t;

/**
 * @brief Hash table for O(1) primitive deduplication
 */
typedef struct {
    primitive_hash_entry_t** buckets;
    size_t bucket_count;
    size_t entry_count;
} primitive_hash_table_t;

// ============================================================================
// HASH MAP FOR CELL ID LOOKUP (open-addressing, generated by alea_hashmap.h)
// ============================================================================

static inline uint64_t cell_id_hash(int id) { return (uint64_t)(unsigned int)id * 0x9E3779B97F4A7C15ULL; }
static inline bool cell_id_eq(int a, int b) { return a == b; }
#include <limits.h>
ALEA_HASHMAP_DEFINE(cell_hashmap, int, int, cell_id_hash, cell_id_eq, INT_MIN);

// ============================================================================
// HASH MAP FOR UNIVERSE ID LOOKUP (open-addressing, generated by alea_hashmap.h)
// ============================================================================

static inline uint64_t universe_id_hash(int id) { return (uint64_t)(unsigned int)id * 0x9E3779B97F4A7C15ULL; }
static inline bool universe_id_eq(int a, int b) { return a == b; }
ALEA_HASHMAP_DEFINE(universe_hashmap, int, int, universe_id_hash, universe_id_eq, INT_MIN);

// ============================================================================
// STATISTICS
// ============================================================================

/**
 * @brief Statistics about the CSG system
 */
typedef struct alea_stats {
    size_t unique_primitives;
    size_t dedup_hits;
    size_t dedup_saved_bytes;
    size_t current_memory;
    size_t peak_memory;
    size_t reallocations;
    size_t surfaces_converted;
    size_t failed_surfaces;
    size_t cells_converted;
} alea_stats_t;

// ============================================================================
// THE MAIN CSG SYSTEM
// ============================================================================

/**
 * @brief Main CSG system structure
 *
 * This is the ONLY structure needed for the entire CSG system.
 * It replaces the old multi-layer architecture (node_pool, primitive_pool,
 * world, parse_context) with a single flattened structure.
 *
 * Key design principles:
 * - Everything in one place
 * - Automatic growth (no manual sizing)
 * - Direct array access (users can read sys->nodes.data[id] directly)
 * - Built-in deduplication
 * - Built-in statistics
 */
typedef struct alea_system {
    // Core data arrays
    alea_node_vec_t nodes;

    // Primitives with deduplication
    alea_primitive_vec_t primitives;
    alea_plane_data_vec_t primitive_planes;
    alea_sphere_data_vec_t primitive_spheres;
    alea_cylinder_x_data_vec_t primitive_cyl_x;
    alea_cylinder_y_data_vec_t primitive_cyl_y;
    alea_cylinder_z_data_vec_t primitive_cyl_z;
    alea_cone_x_data_vec_t primitive_cone_x;
    alea_cone_y_data_vec_t primitive_cone_y;
    alea_cone_z_data_vec_t primitive_cone_z;
    alea_box_data_vec_t primitive_boxes;
    alea_quadric_data_vec_t primitive_quadrics;
    alea_torus_data_vec_t primitive_toruses;
    alea_rcc_data_vec_t primitive_rccs;
    alea_box_general_data_vec_t primitive_box_generals;
    alea_sph_data_vec_t primitive_sphs;
    alea_trc_data_vec_t primitive_trcs;
    alea_ell_data_vec_t primitive_ells;
    alea_rec_data_vec_t primitive_recs;
    alea_wed_data_vec_t primitive_weds;
    alea_rhp_data_vec_t primitive_rhps;
    alea_arb_data_vec_t primitive_arbs;

    // Primitive lookup for deduplication
    primitive_hash_table_t* primitive_index;

    // Mappings
    alea_surface_vec_t surfaces;

    // Fast lookup: MCNP surface ID -> CSG node ID (built after surface conversion)
    alea_node_id_t* surface_lookup;
    size_t surface_lookup_size;

    // Reverse mapping: primitive_id -> surface index (for cell surface index)
    uint32_t* prim_to_surface;      // prim_to_surface[prim_id] = surface index (UINT32_MAX if none)
    size_t prim_to_surface_size;

    // Reverse mapping: mc_surface_id -> surface index (for fast lookup)
    uint32_t* mc_id_to_surface;   // mc_id_to_surface[mcnp_id] = surface index (UINT32_MAX if none)
    size_t mc_id_to_surface_size;

    alea_material_vec_t materials;

    /* Mixtures (references to pure materials with fractions) */
    alea_mixture_vec_t mixtures;

    // Cells
    alea_cell_vec_t cells;
    cell_hashmap_t cell_index;  // O(1) cell ID -> cell index lookup

    /* Cell complement reference tracking */
    alea_cell_ref_vec_t cell_refs;


    // Transforms (TRn cards and inline from FILL)
    alea_transform_vec_t transforms;
    int next_inline_transform_id;  // Next auto-assigned ID for inline transforms

    // Auto-assigned surface ID counter (for programmatic surface creation)
    int next_auto_surface_id;

    // Auto-assigned material ID counter (for programmatic material creation)
    int next_auto_material_id;

    // Universe index (built after cell conversion)
    alea_universe_vec_t universes;
    universe_hashmap_t  universe_index;  // O(1) universe_id -> vec index
    bool universe_index_built;
    bool has_lattice;  // True if any cell has lat_type != 0

    // Configuration
    alea_config_t config;

    // Statistics
    alea_stats_t stats;

    /* Cell event callbacks (NULL = no-op, set by module models) */
    void* cell_hook_userdata;
    void (*on_cell_added)(void* ud, size_t new_index);
    void (*on_cell_copied)(void* ud, size_t dst_index, size_t src_index);
    void (*on_cell_removed)(void* ud, size_t index);

    /* Surface BVH acceleration (prepared query cache) */
    struct alea_bvh* surface_bvh;
    bool bvh_dirty;  /* True if surfaces changed since BVH build */

    /* Cell adjacency (built by alea_build_cell_adjacency) */
    bool cell_adjacency_built;
    struct alea_cell_neighbor* neighbor_pool;  /* Single allocation for all neighbor data */

    /* Persistent surface->cell map (surface index -> cells referencing that
     * surface, with sense). Unlike the per-cell neighbor lists, this is NOT
     * pruned for high-degree surfaces, so it supports neighbor-walking across
     * surfaces shared by many cells. Built alongside the neighbor lists. */
    struct alea_surface_cell_ref* surf_cell_pool;  /* all refs, grouped by surface */
    size_t* surf_cell_offset;  /* [num_surfaces] offset into surf_cell_pool */
    size_t* surf_cell_count;   /* [num_surfaces] ref count per surface */
    size_t surf_cell_num_surfaces;

    /* Hierarchical spatial index (prepared query cache) */
    struct alea_hier_spatial_index* hier_spatial_index;
    struct alea_volume_path_index* volume_path_index;
    atomic_int spatial_build_state;  /* 0=pending, 1=building, 2=done */
    atomic_uint query_cache_state;
    atomic_flag query_cache_build_lock;
    uint64_t system_id;
    atomic_uint_fast64_t geometry_generation;

    /* Source tracking - how this system was created */
    enum {
        ALEA_SOURCE_EMPTY = 0,       /* Created empty via alea_system_create() */
        ALEA_SOURCE_MCNP,            /* Loaded from MCNP file */
        ALEA_SOURCE_OPENMC,          /* Loaded from OpenMC file */
        ALEA_SOURCE_PROGRAMMATIC     /* Cells added programmatically */
    } source;

    /* Next auto-assigned cell ID (for programmatic creation) */
    int next_auto_cell_id;
} alea_system_t;

int alea_system_query_cache_ready(const alea_system_t* sys, unsigned flags);
uint64_t alea_system_geometry_generation(const alea_system_t* sys);
int alea_system_prepare_query_caches(alea_system_t* sys, unsigned flags);
void alea_system_invalidate_query_caches(alea_system_t* sys, unsigned flags);
void alea_volume_path_index_free(struct alea_volume_path_index* idx);

// ============================================================================
// API - SYSTEM MANAGEMENT
// ============================================================================

/**
 * @brief Create a new CSG system
 *
 * Creates a system with default initial capacities. Arrays grow automatically
 * as needed, so you don't need to guess sizes.
 *
 * @return Pointer to new system, or NULL on failure
 */
alea_system_t* alea_system_create(void);

/**
 * @brief Free all internal resources of a CSG system without freeing the struct itself
 *
 * @param sys System whose internals to destroy
 */
void alea_system_destroy_internals(alea_system_t* sys);

/**
 * @brief Destroy a CSG system and free all memory
 *
 * @param sys System to destroy
 */
void alea_system_destroy(alea_system_t* sys);

/**
 * @brief Reset system to initial state (keep allocated memory)
 *
 * @param sys System to reset
 */
void alea_system_reset(alea_system_t* sys);

// ============================================================================
// API - UNIVERSE OPERATIONS
// ============================================================================

/**
 * @brief Build universe index from cells
 *
 * Groups all cells by their universe_id. Must be called after all cells
 * are converted. Computes bounding boxes for each universe.
 *
 * @param sys CSG system
 * @return 0 on success, -1 on error
 */
int alea_build_universe_index(alea_system_t* sys);

/**
 * @brief Build cell adjacency graph
 *
 * Computes which cells neighbor each other across shared surfaces.
 * This enables O(1) cell lookup when crossing a surface.
 *
 * @param sys CSG system
 * @return 0 on success, -1 on failure
 */
int alea_build_cell_adjacency(alea_system_t* sys);

/**
 * @brief Find neighbor cell across a surface
 *
 * Given a cell and a surface it touches, find the cell on the other side.
 *
 * @param sys CSG system
 * @param cell_index Current cell index
 * @param surface_id MCNP surface ID being crossed
 * @return Neighbor cell index, or -1 if not found
 */
int alea_find_neighbor_cell(const alea_system_t* sys,
                           uint32_t cell_index,
                           int surface_id);

/* Entry in the persistent surface->cell map: a cell referencing a surface,
 * with the sense in which it does so. */
struct alea_surface_cell_ref {
    uint32_t cell_index;
    int8_t sense;  /* +1 or -1 */
};

/**
 * @brief Get all cells referencing a surface (by mc_surface_id).
 *
 * Returns the unpruned surface->cell list built by alea_build_cell_adjacency,
 * enabling neighbor-walking across high-degree surfaces. Returns 0 with
 * *out_count == 0 when adjacency is not built or the surface is unknown.
 *
 * @param sys         CSG system
 * @param mc_surface_id MCNP surface ID
 * @param out_refs    Receives pointer into the internal pool (do not free)
 * @param out_count   Receives the number of refs
 * @return 0 on success (including the empty case), -1 on bad arguments
 */
int alea_surface_cells(const alea_system_t* sys, int mc_surface_id,
                       const struct alea_surface_cell_ref** out_refs,
                       size_t* out_count);

/**
 * @brief Get universe by ID
 *
 * @param sys CSG system
 * @param universe_id Universe ID to find
 * @return Pointer to universe, or NULL if not found
 */
const alea_universe_t* alea_get_universe(const alea_system_t* sys, int universe_id);

/**
 * @brief Get all cells in a universe
 *
 * @param sys CSG system
 * @param universe_id Universe ID
 * @param out_cells Output array of cell pointers (caller provides)
 * @param max_cells Maximum cells to return
 * @return Number of cells found, or -1 on error
 */
int alea_get_universe_cells(const alea_system_t* sys, int universe_id,
                           const alea_cell_entry_t** out_cells, size_t max_cells);

/**
 * @brief Get number of universes in the system
 */
size_t alea_universe_count(const alea_system_t* sys);

/**
 * @brief Get material ID at a point
 * @return Material ID, or ALEA_MATERIAL_NONE if in void
 */
alea_material_id_t alea_material_at_point(alea_system_t* sys, double x, double y, double z);

/**
 * @brief Add a mixture to the system
 *
 * Mixtures store references to pure materials with fractions, allowing
 * export with comments explaining the composition.
 *
 * @param sys CSG system
 * @param mixture Mixture to add (copied into system)
 * @return Result with mixture index, or error code
 */
int alea_add_mixture(alea_system_t* sys, const alea_mixture_t* mixture);


// ============================================================================
// API - VOID GENERATION
// ============================================================================

/**
 * @brief Identify which cell contains a point
 *
 * Tests point against all cells in base universe (universe 0).
 *
 * @param sys CSG system
 * @param x X coordinate
 * @param y Y coordinate
 * @param z Z coordinate
 * @return Cell index (into sys->cells), or -1 if in void
 */
int alea_identify_cell_at_point(alea_system_t* sys, double x, double y, double z);

/**
 * @brief Find overlapping cells
 *
 * Detects cells that occupy the same space (geometry errors).
 * Uses bounding box pre-filtering then point sampling.
 *
 * @param sys CSG system
 * @param out_pairs Output array of cell index pairs [i1,j1,i2,j2,...]
 * @param max_pairs Maximum pairs to return
 * @return Number of overlapping pairs found
 */
int alea_find_overlaps(alea_system_t* sys, int* out_pairs, size_t max_pairs);

// ============================================================================
// API - NODE OPERATIONS
// ============================================================================

/**
 * @brief Allocate a new node
 *
 * Automatically grows the node array if needed.
 *
 * @param sys CSG system
 * @return Result with node ID, or error code on failure
 */
alea_node_id_t alea_alloc_node(alea_system_t* sys);

/**
 * @brief Get pointer to a node
 *
 * @param sys CSG system
 * @param id Node ID
 * @return Pointer to node, or NULL if invalid ID
 */
alea_node_t* alea_get_node(alea_system_t* sys, uint32_t id);

/**
 * @brief Clone a primitive node with specified sense
 *
 * Creates a new node referencing the same primitive but with the given sense.
 * Use this when you need the same surface with opposite side selection.
 *
 * Example:
 *   int idx = alea_sphere_surface(sys, 1, 0, 0, 0, 5.0);
 *   alea_node_id_t sphere_in = alea_surface_at(sys, idx)->neg_node;   // interior
 *   alea_node_id_t sphere_out = alea_surface_at(sys, idx)->pos_node;  // exterior
 *
 * @param sys CSG system
 * @param node_id Source primitive node ID
 * @param sense -1 for interior, +1 for exterior
 * @return New node ID, or ALEA_NODE_ID_INVALID if not a primitive node
 */
alea_node_id_t alea_clone_primitive(alea_system_t* sys, alea_node_id_t node_id, int8_t sense);

// ============================================================================
// API - PRIMITIVE OPERATIONS (with automatic deduplication)
// ============================================================================

/**
 * @brief Get or create a primitive (with deduplication)
 *
 * If an identical primitive already exists (within tolerance), returns
 * its ID. Otherwise, creates a new primitive.
 *
 * This is the core deduplication function.
 *
 * @param sys CSG system
 * @param type Primitive type
 * @param data Primitive data
 * @return Result with primitive ID, or error code
 */
alea_primitive_id_t alea_get_or_create_primitive(alea_system_t* sys,
                                               alea_primitive_type_t type,
                                               alea_primitive_data_t* data,
                                               int8_t* inverted);

/**
 * @brief Create a primitive node
 *
 * Creates a CSG node referencing a primitive with given sense.
 *
 * @param sys CSG system
 * @param primitive_id Primitive ID to reference
 * @param sense -1 for interior, +1 for exterior
 * @param inverted Whether the normal was flipped during canonicalization
 * @param mc_surface_id MCNP surface ID for tracking (0 if none)
 * @return Node ID
 */
alea_node_id_t alea_add_primitive_node(alea_system_t* sys, uint32_t primitive_id,
                                      int8_t sense, int8_t inverted,
                                      int32_t mc_surface_id);

/**
 * @brief Get pointer to a primitive entry
 *
 * @param sys CSG system
 * @param id Primitive ID
 * @return Pointer to primitive entry (read-only)
 */
const alea_primitive_entry_t* alea_get_primitive(const alea_system_t* sys, uint32_t id);
bool alea_primitive_copy_data(const alea_system_t* sys,
                              uint32_t id,
                              alea_primitive_data_t* out);
const void* alea_primitive_payload_const(const alea_system_t* sys,
                                         uint32_t id);

// ============================================================================
// API - UTILITIES
// ============================================================================

/**
 * @brief Print statistics about the CSG system
 *
 * @param sys CSG system
 */
void alea_system_print_stats(const alea_system_t* sys);

/**
 * @brief Compute current memory usage of all dynamic arrays in the system
 * @param sys CSG system
 * @return Total allocated bytes
 */
size_t alea_system_memory_usage(const alea_system_t* sys);

/**
 * @brief Release unused capacity from the system's large arrays.
 *
 * After loading, the 2x-growth vectors (especially nodes) can hold up to ~50%
 * unused tail capacity. Reallocs them down to their current count. Safe to call
 * once geometry is final; do not call mid-construction.
 */
void alea_system_shrink_to_fit(alea_system_t* sys);

/**
 * @brief Get statistics
 *
 * @param sys CSG system
 * @return Pointer to stats (read-only)
 */
const alea_stats_t* alea_system_get_stats(const alea_system_t* sys);

/**
 * @brief Reset statistics
 *
 * @param sys CSG system
 */
void alea_system_reset_stats(alea_system_t* sys);

// ============================================================================
// TOLERANCE AND COMPARISON FUNCTIONS
// ============================================================================

/**
 * @brief Compare two numbers with tolerance
 *
 * @param a First number
 * @param b Second number
 * @param tol Tolerance configuration
 * @return True if numbers are equal within tolerance
 */
bool alea_numbers_equal(double a, double b, const alea_config_t* tol);

/**
 * @brief Canonicalize a primitive to standard form
 *
 * Modifies the primitive data to ensure consistent representation.
 * For example, normalizes plane normals and orders bounding box coordinates.
 *
 * @param type Primitive type
 * @param data Primitive data (modified in place)
 */
void alea_canonicalize_primitive(alea_primitive_type_t type,
                                alea_primitive_data_t* data,
                                int8_t* inverted);

/**
 * @brief Compare two primitives for equality
 *
 * @param type1 First primitive type
 * @param data1 First primitive data
 * @param type2 Second primitive type
 * @param data2 Second primitive data
 * @param tol Tolerance configuration
 * @param match_inverted If non-NULL, set to 1 when the match requires an
 *                       inversion (e.g. planes with opposite normals)
 * @return True if primitives are equal within tolerance
 */
bool alea_primitives_equal(alea_primitive_type_t type1, const alea_primitive_data_t* data1,
                         alea_primitive_type_t type2, const alea_primitive_data_t* data2,
                         const alea_config_t* tol,
                         int8_t* match_inverted);

// ============================================================================
// HASH TABLE OPERATIONS (internal, but exposed for implementation)
// ============================================================================

/**
 * @brief Create a new hash table for primitive deduplication
 *
 * @return Pointer to new hash table, or NULL on failure
 */
primitive_hash_table_t* primitive_hash_table_create(void);

/**
 * @brief Destroy a hash table
 *
 * @param table Hash table to destroy
 */
void primitive_hash_table_destroy(primitive_hash_table_t* table);

/**
 * @brief Compute hash for a primitive
 *
 * @param type Primitive type
 * @param data Primitive data
 * @return Hash value
 */
uint64_t alea_compute_primitive_hash(alea_primitive_type_t type, const alea_primitive_data_t* data,
                                     const alea_config_t* config);

/**
 * @brief Find primitive in hash table
 *
 * @param table Hash table
 * @param sys CSG system (for primitive comparison)
 * @param type Primitive type to find
 * @param data Primitive data to find
 * @param hash Pre-computed hash value
 * @return Primitive ID if found, or ALEA_PRIMITIVE_ID_INVALID if not found
 */
uint32_t primitive_hash_table_find(primitive_hash_table_t* table, alea_system_t* sys,
                                   alea_primitive_type_t type, alea_primitive_data_t* data,
                                   uint64_t hash, int8_t* match_inverted);

/**
 * @brief Insert primitive into hash table
 *
 * @param table Hash table
 * @param primitive_id Primitive ID to insert
 * @param hash Hash value for the primitive
 */
void primitive_hash_table_insert(primitive_hash_table_t* table, uint32_t primitive_id, uint64_t hash);


// ============================================================================
// ERROR HANDLING API (thread-local storage - safe for parallel use)
// ============================================================================

void alea_set_error(alea_error_t code, const char* fmt, ...);

const char* alea_get_error(void);


#endif // ALEA_SYSTEM_H
