// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
// SPDX-License-Identifier: MPL-2.0
#ifndef ALEA_INTERVAL_INTERNAL_H
#define ALEA_INTERVAL_INTERNAL_H

#include "interval_index.h"
#include "core/alea_system.h"
#include "util/compat.h"

typedef struct {
    uint32_t first,count,right,axis,packet;double split;
} spatial_node_t;
ALEA_VEC_DEFINE(spatial_nodes, spatial_node_t);
ALEA_VEC_DEFINE(cell_ids, uint32_t);

#include "raycast/ray_box_sign.h"

enum { GUARD_EMPTY=0, GUARD_ANY=1 };
#define GUARD_BOX UINT32_C(0x20000000)
#define GUARD_REF UINT32_C(0x40000000)
#define GUARD_BRANCH UINT32_C(0x80000000)
#define GUARD_UNION UINT32_C(0x40000000)
typedef struct { uint32_t a,b; } guard_pair_t;
ALEA_VEC_DEFINE(guard_pairs,guard_pair_t);
typedef struct {
    guard_pairs_t pairs;
    uint32_t* bound_ids;
    alea_node_bbox_t* certified;
    uint32_t* table;
    size_t capacity;
    bool failed;
} guard_pool_t;

typedef struct {
    float lo[3], hi[3];
    uint32_t first, count;float dl[6],dh[6];
} patch_bounds_t;
typedef struct {uint32_t first,count,cell,program,active;} spatial_patch_t;
enum { BOUND_LANES=4, BOUND_AXES=9 };
typedef struct {float lo[BOUND_LANES],hi[BOUND_LANES];} packet_axis_t;
typedef struct {packet_axis_t xyz[3];uint32_t first,mask;} bound_packet_t;
ALEA_VEC_DEFINE(packet_axes,packet_axis_t);
ALEA_VEC_DEFINE(bound_packets,bound_packet_t);
ALEA_VEC_DEFINE(spatial_patches,spatial_patch_t);
typedef struct {double center[3];uint32_t first[8],count[8],solid[8];} point_tiles_t;
enum { POINT_AXES=511, POINT_LATTICE=512, POINT_OFFSET_BITS=10 };
typedef struct {uint32_t program,bounds;} point_item_t;
ALEA_VEC_DEFINE(point_tile_sets,point_tiles_t);
ALEA_VEC_DEFINE(point_items,point_item_t);

typedef struct {alea_plane_data_t data;} query_plane_t;
typedef struct {alea_quadric_data_t data;double linear_scale,quadratic_scale;} query_quadric_t;
typedef struct {uint32_t shape;int number;} query_card_t;
typedef struct {
    uint32_t count;
    int mc_cell_id,material_id,fill_universe;
    bool lattice,transformed,container;
    double density;
} query_cell_t;

typedef struct {
    double low[3],scale[3];
    uint32_t first;
    bool ready;
} point_grid_t;
ALEA_VEC_DEFINE(query_cards,query_card_t);
ALEA_VEC_DEFINE(query_boxes,alea_node_bbox_t);
#define CARD_PLANE UINT32_C(0x40000000)
#define CARD_QUADRIC UINT32_C(0x80000000)
#define CARD_ID UINT32_C(0x3fffffff)

#include "interval_program.h"

typedef struct {uint32_t key,index;} point_order_t;
typedef struct {double a,b,c,error;uint32_t stamp;} ray_gate_t;
/* Root slots are direct-mapped by primitive ID. A collision only causes a
 * fresh intersection, and never changes the ordered list of candidate cells. */
typedef struct { uint32_t generation, primitive; int count; double t[4]; } cached_roots_t;
typedef struct { uint32_t generation, slot; } visited_cell_t;
struct regional_worker {
    const interval_index_t* index;
    cached_roots_t* roots;
    size_t root_mask;
    uint32_t generation;
    visited_cell_t* visited;
    ray_gate_t* plane_gates;
    ray_gate_t* quadric_gates;
    struct interval_level* levels;
    interval_fallback_stats_t fallback;
};
enum { CHILD_BOUND_BLOCK_SIZE = 256 };
enum { CHILD_REGION_ABSENT, CHILD_REGION_BUILDING, CHILD_REGION_READY };
typedef struct child_bound_block {
    struct child_bound_block* next;
    alea_node_bbox_t envelope;
    alea_node_bbox_t cells[];
} child_bound_block_t;

typedef struct {
    alea_node_bbox_t box;
    uint32_t branch, cell_root;
} render_bound_t;
ALEA_VEC_DEFINE(render_bounds, render_bound_t);

struct interval_index {
    atomic_uint references;
    alea_mapped_file_t mapping;
    void (*detach_mapping)(struct interval_index*);
    size_t leaf_size;
    bool selective;
    bool points_only;
    bool root_only; /* Points independent of hierarchy, with optional child indexes. */
    interval_child_mode_t child_mode;
    bool time_child_points;
    size_t child_count;
    size_t scan_count;
    size_t child_bounds_count;
    struct universe_scans {
        atomic_uint_fast64_t visits;
        atomic_uint_fast64_t tests;
        atomic_uint_fast64_t local_ns;
        atomic_uint_fast64_t unfiltered_ns;
        atomic_uint_fast64_t build_ns;
        _Atomic size_t covered;
        atomic_int building;
        atomic_int region_state;
        child_bound_block_t* first;
        child_bound_block_t* last;
        interval_index_t* region;
        bool indexed;
        alea_node_bbox_t* bounds;
    }* scans;
    double build_seconds;

    struct proof_entry {uint32_t stamp,node,signs;}* proof_table;
    struct shape_proof {uint32_t stamp,signs;}* primitive_proofs;
    uint32_t proof_epoch;
    regional_worker_t* workers;
    size_t worker_count;
    point_tile_sets_t tiles;point_items_t point_items;
    cell_ids_t point_nodes;
    point_grid_t* point_grids;
    interval_point_t* centers;
    size_t center_capacity;
    point_order_t* order;
    size_t order_capacity;
    uint32_t* active_shapes;
    uint32_t active_epoch;
    query_decisions_t decisions;
    wide_decisions_t wide;
    query_plane_t* planes;
    query_quadric_t* quadrics;
    size_t plane_count,quadric_count;
    uint32_t* shape_id;
    query_boxes_t boxes;
    cell_ids_t patch_cards;
    query_cards_t surface_cards;
    query_cell_t* cell_cards;
    spatial_patches_t patches;
    bound_packets_t packets;
    packet_axes_t extra_axes;
    cell_ids_t cards;
    guard_pool_t guards;
    guard_pairs_t residual;
    uint32_t* box_ids;
    size_t finished_patches, finished_items;
    guard_surface_t* surfaces;
    interval_box_t* bounds;
    alea_bbox_t domain,envelope;alea_bbox_t* domains;
    struct exterior_region {alea_bbox_t box;uint32_t root;} exterior[6];
    bool exterior_ready;
    bool render_bounds_attempted, render_bounds_ready;
    size_t render_void_prefix_count; /* Leading terminal voids in world definition order. */
    uint32_t render_exterior_void_root;
    alea_bbox_t render_core;
    render_bounds_t render_boxes;
    alea_size_vec_t exterior_point_cells;
    double proof_margin;
    alea_system_t* sys;
    uint64_t geometry_generation, material_generation;
    uint32_t* universe_roots;
    spatial_nodes_t nodes;
    struct interval_card { uint32_t primitive; int number; } *card_info;
};
/* The pointer is written once before publication. Edits and teardown require
 * exclusive model access, so a reader can retain it for the whole lookup. */
static inline const interval_index_t* published_child_region(const struct universe_scans* scans) {
    return atomic_load_explicit(&scans->region_state, memory_order_acquire) == CHILD_REGION_READY
        ? scans->region : NULL;
}
void interval_extend_child_bounds(const interval_index_t* idx, const alea_universe_t* u,
                                  size_t covered);

/* A compact instruction has a branch byte and 24 operand bits.
 * Boxes use two adjacent operand tags for their larger node indices. */
enum {
    SHORT_BRANCH_LIMIT = 120,
    BRANCH_ACCEPT_REJECT = 240,
    BRANCH_REJECT_ACCEPT,
    BRANCH_NEXT_ACCEPT,
    BRANCH_ACCEPT_NEXT,
    BRANCH_WIDE = 255
};
static inline query_decision_t unpack_decision(const interval_index_t* idx, uint32_t word) {
    unsigned branch = word >> 24;
    if (branch == BRANCH_WIDE) return idx->wide.data[word & 0xffffff];
    static const uint32_t bases[] = {
        0, QUERY_PLANE, QUERY_QUADRIC, GUARD_BOX, GUARD_BOX + (1u << 21), GUARD_REF
    };
    unsigned tag = (word >> 21) & 7;
    uint32_t operand = bases[tag] + (word & 0x1fffff);
    if (branch < SHORT_BRANCH_LIMIT) return (query_decision_t){operand, 1, branch};
    if (branch < BRANCH_ACCEPT_REJECT)
        return (query_decision_t){operand, branch - SHORT_BRANCH_LIMIT, 1};
    switch (branch) {
    case BRANCH_ACCEPT_REJECT: return (query_decision_t){operand, QUERY_ACCEPT, 0};
    case BRANCH_REJECT_ACCEPT: return (query_decision_t){operand, 0, QUERY_ACCEPT};
    case BRANCH_NEXT_ACCEPT:   return (query_decision_t){operand, 1, QUERY_ACCEPT};
    default:                   return (query_decision_t){operand, QUERY_ACCEPT, 1};
    }
}
void alea_regional_free_scratch(struct interval_level* scratch);
#endif
