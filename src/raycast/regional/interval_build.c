// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
// SPDX-License-Identifier: MPL-2.0
/* Build per-universe spatial and expression indexes for selected-owner queries. */
#include "interval_internal.h"
#include "raycast/ray_traversal.h"
#include "util/compat.h"
#include "util/load_profile.h"
#include "util/alea_parallel.h"
#include "util/alea_alloc.h"

#include "core/alea_regional.h"
#include "core/alea_universe.h"
#include "primitives/bbox.h"
#include <float.h>
#include <math.h>
#include <stdlib.h>

#include "interval_bounds.h"
#include "interval_guard.h"
#include "interval_split.h"
#include "interval_region_proof.h"
#include "interval_point_tiles.h"

#include "interval_projections.h"
/* Cell-root certificates can also eliminate complemented cells wherever
 * their original CSG subtree is shared by another expression. */
static int prepare_root_bounds(interval_index_t* idx, const alea_universe_t* child) {
    const alea_system_t* s=idx->sys;
    size_t n = child ? child->cell_indices.count : s->cells.count;
    idx->guards.bound_ids=calloc(s->nodes.count+1,sizeof(uint32_t));
    idx->guards.certified=calloc(n+1,sizeof(alea_node_bbox_t));
    if(!idx->guards.bound_ids || !idx->guards.certified)return -1;
    size_t count=0;
    for(size_t j=0;j<n;j++) {
        size_t i = child ? child->cell_indices.data[j] : j;
        const alea_cell_entry_t* c=s->cells.data+i;
        uint32_t root=c->root_node_id;
        if (!child && idx->root_only && c->universe_id != 0) {
            if (!idx->child_count) continue;
            const alea_universe_t* u = alea_get_universe(s, c->universe_id);
            if (!u || !idx->scans[u - s->universes.data].indexed) continue;
        }
        if(root>=s->nodes.count || (c->lat_type && c->lat_fill) ||
           idx->guards.bound_ids[root] || c->surface_index_count<4)continue;
        interval_box_t current=interval_box_from(alea_node_bbox_get(&s->nodes.data[root].bbox));
        if(current.hi[0]-current.lo[0]<1e5 && current.hi[1]-current.lo[1]<1e5 &&
           current.hi[2]-current.lo[2]<1e5)continue;
        alea_bbox_t box;
        if(alea_tighten_bbox_mixed_constraints(s,root,&box))continue;
        interval_box_t certified=interval_box_from(box);
        for(unsigned k=0;k<3;k++) {
            double pad=1e-6*(1+fabs(certified.lo[k])+fabs(certified.hi[k]));
            certified.lo[k]-=pad;certified.hi[k]+=pad;
        }
        box=(alea_bbox_t){certified.lo[0],certified.hi[0],certified.lo[1],certified.hi[1],
                          certified.lo[2],certified.hi[2]};
        alea_node_bbox_set(idx->guards.certified+count,&box);
        idx->guards.bound_ids[root]=(uint32_t)++count;
    }
    return 0;
}
static interval_box_t bounded_cell_domain(const interval_index_t* idx,uint32_t id,alea_bbox_t box) {
    const alea_cell_entry_t* cell=idx->sys->cells.data+id;
    if(cell->lat_type && cell->lat_fill)return interval_box_from(box);
    uint32_t root=cell->root_node_id;
    if(root<idx->sys->nodes.count && idx->guards.bound_ids[root]) {
        alea_bbox_t bound=alea_node_bbox_get(idx->guards.certified+idx->guards.bound_ids[root]-1);
        box=alea_bbox_intersection(&box,&bound);
    }
    return cell_domain(idx->sys,id,box);
}
typedef struct {
    interval_index_t* index;
    const alea_universe_t* universe;
    alea_bbox_t domain;
} cell_bounds_work_t;

static int prepare_cell_bounds(void* context, size_t worker, size_t begin, size_t end) {
    (void)worker;
    const cell_bounds_work_t* work = context;
    for (size_t j = begin; j < end; ++j) {
        uint32_t id = (uint32_t)work->universe->cell_indices.data[j];
        work->index->bounds[id] = bounded_cell_domain(work->index, id, work->domain);
    }
    return 0;
}

#include "interval_domains.h"

/* One reader extends the immutable prefix. Readers that see an older prefix
 * can finish their lookup by testing the remaining cells in definition order. */
void interval_extend_child_bounds(const interval_index_t* idx, const alea_universe_t* u,
                                  size_t covered) {
    struct universe_scans* scans = &idx->scans[u - idx->sys->universes.data];
    int idle = 0;
    if (!atomic_compare_exchange_strong(&scans->building, &idle, 1)) return;
    if (atomic_load(&scans->covered) != covered || alea_interrupted()) goto done;

    atomic_store(&scans->unfiltered_ns, 0);
    double start = alea_monotonic_seconds();
    size_t count = u->cell_indices.count - covered;
    if (count > CHILD_BOUND_BLOCK_SIZE) count = CHILD_BOUND_BLOCK_SIZE;
    child_bound_block_t* block = alea_malloc(sizeof(*block) + count * sizeof(*block->cells));
    if (!block) goto done;
    block->next = NULL;
    alea_bbox_t envelope = alea_bbox_empty();
    for (size_t i = 0; i < count; ++i) {
        const alea_cell_entry_t* cell = &idx->sys->cells.data[u->cell_indices.data[covered + i]];
        if (cell->lat_type && cell->lat_fill) {
            alea_node_bbox_set(&block->cells[i], &idx->domain);
        } else if (cell->root_node_id < idx->sys->nodes.count) {
            block->cells[i] = idx->sys->nodes.data[cell->root_node_id].bbox;
        } else {
            alea_bbox_t empty = alea_bbox_empty();
            alea_node_bbox_set(&block->cells[i], &empty);
        }
        alea_bbox_t box = alea_node_bbox_get(&block->cells[i]);
        envelope = alea_bbox_union(&envelope, &box);
    }
    alea_node_bbox_set(&block->envelope, &envelope);
    if (covered) scans->last->next = block;
    else scans->first = block;
    scans->last = block;
    uint64_t ns = (uint64_t)((alea_monotonic_seconds() - start) * 1e9);
    atomic_fetch_add_explicit(&scans->build_ns, ns, memory_order_relaxed);
    atomic_store(&scans->covered, covered + count);
done:
    atomic_store(&scans->building, 0);
}

/* Bounds and unspecialized references retain their original subtree. Its
 * physical primitives (including aliases) can still supply a crossing. */
static bool mark_original_shapes(interval_index_t* idx,uint32_t node,unsigned* budget,unsigned depth) {
    if(!*budget || depth>=768 || node>=idx->sys->nodes.count)return false;
    --*budget;
    const alea_node_t* n=idx->sys->nodes.data+node;
    alea_operation_t op=ALEA_GET_OPERATION(n);
    if(op==ALEA_OP_PRIMITIVE) {
        uint32_t id=n->primitive.primitive_id;
        if(id>=idx->sys->primitives.count)return false;
        idx->active_shapes[id]=idx->active_epoch;
        return true;
    }
    return mark_original_shapes(idx,n->operation.left,budget,depth+1) &&
           (op==ALEA_OP_COMPLEMENT || mark_original_shapes(idx,n->operation.right,budget,depth+1));
}
static bool mark_active_shapes(interval_index_t* idx,uint32_t code,unsigned* budget) {
    if(code<=GUARD_ANY)return true;
    if(!*budget)return false;
    --*budget;
    if(code&GUARD_BRANCH) {
        guard_pair_t pair=idx->guards.pairs.data[code&~(GUARD_BRANCH|GUARD_UNION)];
        return mark_active_shapes(idx,pair.a,budget) && mark_active_shapes(idx,pair.b,budget);
    }
    if(code&(GUARD_REF|GUARD_BOX))
        return mark_original_shapes(idx,(code&~(GUARD_REF|GUARD_BOX))>>1,budget,0);
    idx->active_shapes[(code-2)>>1]=idx->active_epoch;
    return true;
}

/* The projected box of the three primary axes already certifies many
 * diagonal slabs. Store and test only directions that tighten a lane. */
static int store_packet(interval_index_t* idx,const packet_axis_t axes[BOUND_AXES],unsigned mask) {
    if(idx->extra_axes.count>UINT32_MAX-6)return -1;
    bound_packet_t p={.first=(uint32_t)idx->extra_axes.count,.mask=mask};
    memcpy(p.xyz,axes,sizeof(p.xyz));
    for(unsigned k=0;k<6;k++)if(mask&(1u<<k)) {
        if(alea_vec_push(&idx->extra_axes,axes[3+k],packet_axis_t))return -1;
    }
    return alea_vec_push(&idx->packets,p,bound_packet_t);
}
static int spatial_leaf(interval_index_t* idx, alea_bbox_t box, uint32_t node, const uint32_t* ids,
                        size_t count) {
    if (count > UINT32_MAX - idx->patches.count) return -1;
    size_t start=idx->patches.count;
    uint32_t packet_start=idx->packets.count;
    packet_axis_t axes[BOUND_AXES]={0};
    unsigned mask=0;
    for (size_t i=0;i<count;i++) {
        uint32_t id=ids[i];
        const alea_cell_entry_t* cell=&idx->sys->cells.data[id];
        interval_box_t b=interval_box_from(box);
        for (unsigned a=0;a<3;a++) {
            b.lo[a]=fmax(b.lo[a],idx->bounds[id].lo[a]);
            b.hi[a]=fmin(b.hi[a],idx->bounds[id].hi[a]);
        }
        if (box_empty(b)) continue;
        if (!(cell->lat_type && cell->lat_fill)) {
            unsigned steps=2048;
            for (unsigned pass=0;pass<3 && steps;pass++) {
                interval_box_t before=b;
                b=contract_cell(idx->sys,cell->root_node_id,b,false,0,&steps);
                if (box_empty(b) || !memcmp(&b,&before,sizeof(b))) break;
            }
            if (box_empty(b)) continue;
        }
        /* Leave room for the rounded spatial cuts and local ray evaluation. */
        for (unsigned a=0;a<3;a++) {
            double margin=512*DBL_EPSILON*(1+fabs(b.lo[a])+fabs(b.hi[a]));
            b.lo[a]-=margin; b.hi[a]+=margin;
        }
        if (idx->cards.count>UINT32_MAX) return -1;
        patch_bounds_t patch={.first=(uint32_t)idx->cards.count};
        for (unsigned a=0;a<3;a++) {
            patch.lo[a]=nextafterf((float)b.lo[a],-INFINITY);
            patch.hi[a]=nextafterf((float)b.hi[a],INFINITY);
        }
        /* The ray slabs include an FP64 arithmetic pad. Certify signs on
         * the full outward-rounded box that a query can actually visit. */
        for (unsigned a=0;a<3;a++) {
            b.lo[a]=(double)patch.lo[a]-idx->proof_margin;
            b.hi[a]=(double)patch.hi[a]+idx->proof_margin;
        }
        unsigned budget=32768;
        uint32_t test=cell->lat_type && cell->lat_fill ? GUARD_ANY :
            compile_guard(idx->sys,&idx->guards,idx->surfaces,cell->root_node_id,b,false,0,&budget);
        if(idx->guards.failed)return -1;
        if(test==GUARD_EMPTY)continue;
        interval_box_t views[3];
        for (unsigned a=0;a<3;a++) {
            unsigned j=(a+1)%3,remaining=256;
            interval_box_t q=rotate_box(b,a),o=q;
            for(unsigned pass=0;pass<2 && remaining;pass++) {
                q=oriented_bound(idx,test,q,a,&remaining);
                if (box_empty(q)) {q=o;break;}
            }
            patch.dl[2*a]=nextafterf((float)q.lo[a],-INFINITY);
            patch.dh[2*a]=nextafterf((float)q.hi[a],INFINITY);
            patch.dl[2*a+1]=nextafterf((float)q.lo[j],-INFINITY);
            patch.dh[2*a+1]=nextafterf((float)q.hi[j],INFINITY);
            views[a]=o;
            views[a].lo[a]=(double)patch.dl[2*a]-4*idx->proof_margin;
            views[a].hi[a]=(double)patch.dh[2*a]+4*idx->proof_margin;
            views[a].lo[j]=(double)patch.dl[2*a+1]-4*idx->proof_margin;
            views[a].hi[j]=(double)patch.dh[2*a+1]+4*idx->proof_margin;
        }
        if (idx->guards.failed) return -1;
        if (test==GUARD_EMPTY) continue;
        unsigned extra=1024;
        test=fold_views(idx,test,b,views,&extra);
        if(idx->guards.failed)return -1;
        if(test==GUARD_EMPTY)continue;
        uint32_t active = 0;
        if (!idx->points_only) {
            for (size_t k=0;k<cell->surface_index_count;k++) {
                uint32_t si=cell->surface_indices[k];
                uint32_t prim=idx->sys->surfaces.data[si].primitive_id;
                if (idx->surfaces[prim].kind && sign_in_views(&idx->surfaces[prim],b,views)) continue;
                if (alea_vec_push(&idx->cards,si,uint32_t)) return -1;
            }
            patch.count=(uint32_t)(idx->cards.count-patch.first);
            if(++idx->active_epoch==0) {
                memset(idx->active_shapes,0,idx->sys->primitives.count*sizeof(*idx->active_shapes));
                idx->active_epoch=1;
            }
            unsigned visits=4096;
            bool complete=mark_active_shapes(idx,test,&visits);
            for(uint32_t j=0;j<patch.count;j++) {
                uint32_t* cards=idx->cards.data+patch.first;
                uint32_t si=cards[j];
                struct interval_card info=idx->card_info[si];
                if(!complete || idx->active_shapes[info.primitive]==idx->active_epoch) {
                    cards[j]=cards[active];cards[active++]=si;
                }
            }
            if (active==patch.count && patch.count==cell->surface_index_count) {
                idx->cards.count=patch.first;
                patch.count=UINT32_MAX;
            }
        }
        unsigned lane=(idx->patches.count-start)%BOUND_LANES;
        if(lane==0){memset(axes,0,sizeof(axes));mask=0;}
        for(unsigned k=0;k<BOUND_AXES;k++){
            axes[k].lo[lane]=k<3?patch.lo[k]:patch.dl[k-3];
            axes[k].hi[lane]=k<3?patch.hi[k]:patch.dh[k-3];
        }
        for(unsigned i=0;i<3;i++) {
            unsigned j=(i+1)%3;
            double low=patch.lo[i],high=patch.hi[i];
            if(patch.dl[2*i]>low+patch.lo[j] || patch.dh[2*i]<high+patch.hi[j])
                mask|=1u<<(2*i);
            if(patch.dl[2*i+1]>low-patch.hi[j] || patch.dh[2*i+1]<high-patch.lo[j])
                mask|=2u<<(2*i);
        }
        if(lane==BOUND_LANES-1 && store_packet(idx,axes,mask))return -1;
        spatial_patch_t meta={.first=patch.first,.count=patch.count,.cell=id,.program=test,.active=active};
        if (alea_vec_push(&idx->patches,meta,spatial_patch_t)) return -1;
    }
    if((idx->patches.count-start)%BOUND_LANES && store_packet(idx,axes,mask))return -1;
    idx->nodes.data[node]=(spatial_node_t){.packet=packet_start,.first=(uint32_t)start,.count=(uint32_t)(idx->patches.count-start)};
    return prepare_point_tiles(idx,idx->nodes.data+node,box);
}


static int build_region(interval_index_t* idx, alea_bbox_t box, const uint32_t* ids,
                        size_t count, unsigned depth, size_t budget, uint32_t* root) {
    if (alea_interrupted() || idx->nodes.count >= UINT32_MAX) return -1;
    uint32_t node = (uint32_t)idx->nodes.count;
    spatial_node_t blank = {0};
    if (alea_vec_push(&idx->nodes, blank, spatial_node_t)) return -1;
    *root = node;
    if (count <= idx->leaf_size || depth == 24 || budget < 2 * count)
        return spatial_leaf(idx, box, node, ids, count);
    alea_bbox_t halves[2]={box,box};
    unsigned axis;
    double split=select_split(idx,box,ids,count,&axis);
    if(axis==0)halves[0].max_x=halves[1].min_x=split;
    else if(axis==1)halves[0].max_y=halves[1].min_y=split;
    else halves[0].max_z=halves[1].min_z=split;
    if(count>SIZE_MAX/(2*sizeof(uint32_t)))return -1;
    uint32_t* children = malloc(2 * count * sizeof(*children));
    if (!children) return -1;
    size_t sizes[2] = {0};
    for(unsigned side=0;side<2;side++) {
        if(++idx->proof_epoch==0) {
            memset(idx->proof_table,0,65536*sizeof(*idx->proof_table));
            memset(idx->primitive_proofs,0,idx->sys->primitives.count*sizeof(*idx->primitive_proofs));
            idx->proof_epoch=1;
        }
        interval_box_t region=interval_box_from(halves[side]);
        for(size_t i=0;i<count;i++) {
            uint32_t id=ids[i];
            if(!overlaps_box(idx->bounds[id],halves+side))continue;
            const alea_cell_entry_t* cell=idx->sys->cells.data+id;
            unsigned visits=256;
            if((cell->lat_type && cell->lat_fill) ||
               region_signs(idx,cell->root_node_id,region,0,&visits)!=SIGN_POS)
                children[side*count+sizes[side]++]=id;
        }
    }
    size_t total = sizes[0] + sizes[1];
    size_t left_budget = total ? (size_t)((double)budget * sizes[0] / total) : 0;
    uint32_t left, right;
    int rc = build_region(idx, halves[0], children, sizes[0], depth + 1, left_budget, &left);
    if (!rc) rc = build_region(idx, halves[1], children + count, sizes[1], depth + 1,
                               budget - left_budget, &right);
    if (!rc) idx->nodes.data[node] = (spatial_node_t){.first=left,.right=right,.axis=axis,.split=split};
    free(children);
    return rc;
}

static uint32_t place_box(interval_index_t* idx,uint32_t code,uint32_t* ids) {
    if((code&(GUARD_BRANCH|GUARD_REF|GUARD_BOX))!=GUARD_BOX)return code;
    uint32_t node=(code&~GUARD_BOX)>>1;
    if(!ids[node]) {
        if(alea_vec_push(&idx->boxes,idx->sys->nodes.data[node].bbox,alea_node_bbox_t)) {
            idx->guards.failed=true;return code;
        }
        ids[node]=idx->boxes.count;
    }
    return GUARD_BOX|((ids[node]-1)*2+(code&1));
}
static uint32_t runtime_operand(interval_index_t* idx,uint32_t code,uint32_t* boxes) {
    if(code<2 || code&(GUARD_BRANCH|GUARD_REF|GUARD_BOX))return place_box(idx,code,boxes);
    if (idx->root_only) {
        const alea_primitive_entry_t* p = &idx->sys->primitives.data[(code - 2) >> 1];
        if (p->payload_index > (QUERY_OPERAND >> 1)) return code;
        uint32_t operand = 2 * p->payload_index + (code & 1);
        if (p->type == ALEA_PRIMITIVE_PLANE) return QUERY_PLANE | operand;
        if (p->type == ALEA_PRIMITIVE_QUADRIC) return QUERY_QUADRIC | operand;
        return code;
    }
    uint32_t id=(code-2)>>1,shape=idx->shape_id[id],sense=code&1;
    if(shape&CARD_PLANE)return QUERY_PLANE|((shape&CARD_ID)*2+sense);
    if(shape&CARD_QUADRIC)return QUERY_QUADRIC|((shape&CARD_ID)*2+sense);
    return code;
}

/* During emission -1 and -2 are false and true. Nonnegative targets
 * are instruction offsets in this program's own short instruction list. */
static int emit_decision(const interval_index_t* idx,uint32_t code,int yes,int no,
                         query_decision_t* words,unsigned* count) {
    if(yes==no)return yes;
    if(code==GUARD_ANY)return yes;
    if(code==GUARD_EMPTY)return no;
    if(*count>=QUERY_DECISION_LIMIT)return -3;
    if(code&GUARD_BRANCH) {
        guard_pair_t pair=idx->guards.pairs.data[code&~(GUARD_BRANCH|GUARD_UNION)];
        int right=emit_decision(idx,pair.b,yes,no,words,count);
        if(right==-3)return right;
        return (code&GUARD_UNION)?emit_decision(idx,pair.a,yes,right,words,count):
            emit_decision(idx,pair.a,right,no,words,count);
    }
    unsigned pc=(*count)++;
    query_decision_t instruction={code,yes==-1?0:yes==-2?QUERY_ACCEPT:pc-(unsigned)yes,
                                      no==-1?0:no==-2?QUERY_ACCEPT:pc-(unsigned)no};
    words[pc]=instruction;
    return (int)pc;
}
/* Both ray patches and point octants own regional expression roots. */
static uint32_t* regional_program(interval_index_t* idx,size_t i) {
    size_t patches = idx->patches.count - idx->finished_patches;
    if (i < patches) return &idx->patches.data[idx->finished_patches + i].program;
    return &idx->point_items.data[idx->finished_items + i - patches].program;
}
/* Most decision nodes fall through on one branch and take a short jump on the other. */
static unsigned branch_code(uint16_t yes, uint16_t no) {
    if (yes == 1 && no < SHORT_BRANCH_LIMIT) return no;
    if (no == 1 && yes < SHORT_BRANCH_LIMIT) return SHORT_BRANCH_LIMIT + yes;
    if (yes == QUERY_ACCEPT && no == 0) return BRANCH_ACCEPT_REJECT;
    if (yes == 0 && no == QUERY_ACCEPT) return BRANCH_REJECT_ACCEPT;
    if (yes == 1 && no == QUERY_ACCEPT) return BRANCH_NEXT_ACCEPT;
    if (yes == QUERY_ACCEPT && no == 1) return BRANCH_ACCEPT_NEXT;
    return BRANCH_WIDE;
}
static int pack_decision(interval_index_t* idx, query_decision_t word, uint32_t* output) {
    unsigned branch = branch_code(word.yes, word.no);
    uint32_t payload = word.operand & QUERY_OPERAND;
    unsigned tag;
    switch (word.operand & ~QUERY_OPERAND) {
    case 0:             tag = 0; break;
    case QUERY_PLANE:   tag = 1; break;
    case QUERY_QUADRIC: tag = 2; break;
    case GUARD_BOX:     tag = 3; break;
    case GUARD_REF:     tag = 5; break;
    default:            tag = 7; break;
    }
    unsigned bits = tag == 3 ? 22 : 21;
    if (branch != BRANCH_WIDE && tag < 7 && payload < (1u << bits)) {
        *output = (branch << 24) + (tag << 21) + payload;
        return 0;
    }
    if (idx->wide.count > 0xffffff) return 1;
    *output = 0xff000000u | (uint32_t)idx->wide.count;
    return alea_vec_push(&idx->wide, word, query_decision_t) ? -1 : 0;
}
static size_t pending_programs(const interval_index_t* idx) {
    return (idx->patches.count - idx->finished_patches) +
           (idx->point_items.count - idx->finished_items);
}
static int compile_decisions(interval_index_t* idx) {
    size_t n=idx->guards.pairs.count;
    uint32_t* map=calloc(2*n+1,sizeof(*map));
    if(!map)return -1;
    for(size_t i=0;i<pending_programs(idx);i++) {
        uint32_t* slot=regional_program(idx,i);
        uint32_t code=*slot;
        if(!(code&GUARD_BRANCH))continue;
        uint32_t key=(code&~(GUARD_BRANCH|GUARD_UNION))*2+!!(code&GUARD_UNION);
        if(!map[key]) {
            query_decision_t buffer[QUERY_DECISION_LIMIT];unsigned count=0;
            int root=emit_decision(idx,code,-2,-1,buffer,&count);
            if(root<0 || idx->decisions.count+count>=GUARD_BOX)map[key]=UINT32_MAX;
            else{
                size_t begin=idx->decisions.count;
                if(alea_vec_reserve(&idx->decisions,begin+count,uint32_t)) {free(map);return -1;}
                size_t wide=idx->wide.count;int rc=0;
                for(unsigned i=0;i<count && !rc;i++)rc=pack_decision(idx,buffer[i],idx->decisions.data+begin+i);
                if(rc<0){free(map);return -1;}
                if(rc){idx->wide.count=wide;map[key]=UINT32_MAX;}
                else {idx->decisions.count+=count;map[key]=(uint32_t)begin+root+1;}
            }
        }
        if(map[key]!=UINT32_MAX)*slot=QUERY_PROGRAM|(map[key]-1);
    }
    free(map);return 0;
}
/* A physical-card list is reused by many neighbouring spatial leaves. Intern
 * the lists before translating them to the compact query geometry. */
typedef struct {size_t begin,count;uint64_t hash;} card_list_t;
typedef struct { card_list_t* entries; size_t capacity,used; } card_table_t;

static int grow_card_table(card_table_t* table) {
    if (table->capacity > SIZE_MAX / (2 * sizeof(card_list_t))) return -1;
    size_t capacity=table->capacity?2*table->capacity:1024;
    card_list_t* entries=calloc(capacity,sizeof(*entries));
    if(!entries)return -1;
    for(size_t i=0;i<table->capacity;i++) {
        card_list_t entry=table->entries[i];
        if(!entry.count)continue;
        size_t slot=(entry.hash^(entry.hash>>32))&(capacity-1);
        while(entries[slot].count)slot=(slot+1)&(capacity-1);
        entries[slot]=entry;
    }
    free(table->entries);table->entries=entries;table->capacity=capacity;
    return 0;
}
static int intern_cards(interval_index_t* idx,card_table_t* pool,
                        const uint32_t* ids,size_t n,uint32_t* first) {
    if(!n){*first=0;return 0;}
    if(pool->used>=pool->capacity/2 && grow_card_table(pool))return -1;
    size_t mask=pool->capacity-1;
    card_list_t* table=pool->entries;
    uint64_t hash=UINT64_C(14695981039346656037);
    for(size_t i=0;i<n;i++){hash^=ids[i];hash*=UINT64_C(1099511628211);}
    size_t slot=(hash^(hash>>32))&mask;
    while(table[slot].count) {
        card_list_t entry=table[slot];
        if(entry.hash==hash && entry.count==n){
            if (!memcmp(idx->patch_cards.data + entry.begin, ids, n * sizeof(*ids))) {
                *first = (uint32_t)entry.begin;
                return 0;
            }
        }
        slot=(slot+1)&mask;
    }
    if(n>UINT32_MAX-idx->patch_cards.count ||
       alea_vec_reserve(&idx->patch_cards,idx->patch_cards.count+n,uint32_t))return -1;
    *first=(uint32_t)idx->patch_cards.count;
    memcpy(idx->patch_cards.data + *first, ids, n * sizeof(*ids));
    idx->patch_cards.count += n;
    pool->used++;
    table[slot]=(card_list_t){*first,n,hash};return 0;
}
/* Full physical surface lists already belong to the editable model. Only
 * intern reduced lists, using a table scoped to this universe's patches. */
static int finish_patch_cards(interval_index_t* idx) {
    card_table_t table={0};
    int rc=0;
    for(size_t i=idx->finished_patches;!rc && i<idx->patches.count;i++) {
        spatial_patch_t* p=idx->patches.data+i;
        if(p->count==UINT32_MAX) {
            p->first = UINT32_MAX;
            p->count = idx->cell_cards[p->cell].count;
        }else if(p->count)rc=intern_cards(idx,&table,idx->cards.data+p->first,p->count,&p->first);
    }
    free(table.entries);
    return rc;
}
/* Branch operands always refer to an earlier interned pair. Mark the
 * residual DAG from roots that outgrew the short decision programs. */
static uint32_t remap_pair(uint32_t code,const uint32_t* map) {
    if(!(code&GUARD_BRANCH))return code;
    uint32_t flags=code&(GUARD_BRANCH|GUARD_UNION);
    return flags|map[code&~flags];
}
static int compact_pairs(interval_index_t* idx) {
    guard_pairs_t* pairs=&idx->guards.pairs;
    uint32_t* map=calloc(pairs->count+1,sizeof(*map));
    if(!map)return -1;
    for(size_t i=0;i<pending_programs(idx);i++) {
        uint32_t* slot=regional_program(idx,i);
        uint32_t code=*slot;
        if(code&GUARD_BRANCH)map[code&~(GUARD_BRANCH|GUARD_UNION)]=1;
    }
    for(size_t i=pairs->count;i-->0;)if(map[i]) {
        guard_pair_t pair=pairs->data[i];
        if(pair.a&GUARD_BRANCH)map[pair.a&~(GUARD_BRANCH|GUARD_UNION)]=1;
        if(pair.b&GUARD_BRANCH)map[pair.b&~(GUARD_BRANCH|GUARD_UNION)]=1;
    }
    for(size_t i=0;i<pairs->count;i++)if(map[i]) {
        guard_pair_t p=pairs->data[i];
        map[i]=(uint32_t)idx->residual.count;
        guard_pair_t result={remap_pair(p.a,map),remap_pair(p.b,map)};
        if (idx->residual.count >= GUARD_UNION || alea_vec_push(&idx->residual,result,guard_pair_t)) {
            free(map);
            return -1;
        }
    }
    for(size_t i=0;i<pending_programs(idx);i++) {
        uint32_t* slot=regional_program(idx,i);
        *slot=remap_pair(*slot,map);
    }
    free(map);
    alea_vec_free(pairs);
    return 0;
}
static int pack_shapes(interval_index_t* idx) {
    const alea_system_t* sys=idx->sys;
    if(sys->primitives.count>CARD_ID)return -1;
    for(size_t i=0;i<sys->primitives.count;i++) {
        if(sys->primitives.data[i].type==ALEA_PRIMITIVE_PLANE)idx->plane_count++;
        if(sys->primitives.data[i].type==ALEA_PRIMITIVE_QUADRIC)idx->quadric_count++;
    }
    idx->planes=calloc(idx->plane_count+1,sizeof(*idx->planes));
    idx->quadrics=calloc(idx->quadric_count+1,sizeof(*idx->quadrics));
    idx->shape_id=calloc(sys->primitives.count+1,sizeof(*idx->shape_id));
    idx->cell_cards=calloc(sys->cells.count+1,sizeof(*idx->cell_cards));
    if(!idx->planes || !idx->quadrics || !idx->shape_id || !idx->cell_cards)return -1;
    size_t planes=0,quadrics=0;
    for(size_t i=0;i<sys->primitives.count;i++) {
        alea_primitive_type_t type = sys->primitives.data[i].type;
        const void* payload=alea_primitive_payload_const(sys,(uint32_t)i);
        if (!payload) return -1;
        if(type==ALEA_PRIMITIVE_PLANE) {
            idx->planes[planes].data=*(const alea_plane_data_t*)payload;
            idx->shape_id[i]=CARD_PLANE|(uint32_t)planes++;
        }else if(type==ALEA_PRIMITIVE_QUADRIC) {
            idx->quadrics[quadrics].data=*(const alea_quadric_data_t*)payload;
            query_quadric_t* q=idx->quadrics+quadrics;
            for(unsigned j=0;j<6;j++)q->quadratic_scale+=fabs(q->data.coeffs[j]);
            for(unsigned j=6;j<9;j++)q->linear_scale+=fabs(q->data.coeffs[j]);
            idx->shape_id[i]=CARD_QUADRIC|(uint32_t)quadrics++;
        }else idx->shape_id[i]=(uint32_t)i;
    }
    if (!idx->points_only) {
        if (alea_vec_reserve(&idx->surface_cards, sys->surfaces.count, query_card_t)) return -1;
        for (size_t i=0; i<sys->surfaces.count; ++i) {
            struct interval_card card = idx->card_info[i];
            idx->surface_cards.data[i] = (query_card_t){idx->shape_id[card.primitive], card.number};
        }
        idx->surface_cards.count = sys->surfaces.count;
    }
    for (size_t i=0; i<sys->cells.count; ++i) {
        const alea_cell_entry_t* c=sys->cells.data+i;
        if (c->surface_index_count>UINT32_MAX) return -1;
        idx->cell_cards[i]=(query_cell_t){
            .count=(uint32_t)c->surface_index_count, .mc_cell_id=c->mc_cell_id,
            .material_id=c->material_id, .fill_universe=c->fill_universe,
            .lattice=c->lat_type && c->lat_fill, .transformed=c->fill_transform>0,
            .container=alea_cell_entry_is_container(c), .density=c->density
        };
    }
    return 0;
}
static int finish_region_programs(interval_index_t* idx) {
    uint32_t* ids=idx->box_ids;
    for(size_t i=0;i<idx->guards.pairs.count;i++) {
        guard_pair_t* p=idx->guards.pairs.data+i;
        p->a=runtime_operand(idx,p->a,ids);p->b=runtime_operand(idx,p->b,ids);
    }
    for(size_t i=0;i<pending_programs(idx);i++) {
        uint32_t* slot=regional_program(idx,i);
        *slot=runtime_operand(idx,*slot,ids);
    }
    free(idx->guards.table);
    idx->guards.table=NULL;
    idx->guards.capacity=0;
    if(idx->guards.failed || compile_decisions(idx) || compact_pairs(idx))return -1;
    if (!idx->points_only && finish_patch_cards(idx)) return -1;
    alea_vec_free(&idx->cards);
    idx->finished_patches = idx->patches.count;
    idx->finished_items = idx->point_items.count;
    return 0;
}
static int finish_index(interval_index_t* idx) {

    idx->guards.pairs=idx->residual;
    idx->residual=(guard_pairs_t){0};
    free(idx->box_ids);
    idx->box_ids=NULL;
    alea_vec_shrink_to_fit(&idx->tiles,point_tiles_t);
    alea_vec_shrink_to_fit(&idx->point_items,point_item_t);
    alea_vec_shrink_to_fit(&idx->decisions,uint32_t);
    alea_vec_shrink_to_fit(&idx->wide,query_decision_t);
    alea_vec_shrink_to_fit(&idx->patch_cards,uint32_t);
    alea_vec_shrink_to_fit(&idx->boxes,alea_node_bbox_t);
    alea_vec_shrink_to_fit(&idx->patches,spatial_patch_t);
    alea_vec_shrink_to_fit(&idx->packets,bound_packet_t);
    alea_vec_shrink_to_fit(&idx->extra_axes,packet_axis_t);
    alea_vec_shrink_to_fit(&idx->nodes,spatial_node_t);
    free(idx->surfaces);idx->surfaces=NULL;free(idx->card_info);idx->card_info=NULL;
    free(idx->shape_id);idx->shape_id=NULL;
    return 0;
}

/* Select a required union from a product of constraints. Its alternatives
 * cover the cell; the other constraints are retained in the inherited box. */
static uint32_t render_union(const alea_system_t* sys, uint32_t id, unsigned depth, unsigned* size) {
    const alea_node_t* n = sys->nodes.data + id;
    if (ALEA_GET_OPERATION(n) == ALEA_OP_UNION) {
        *size = 1;
        while (ALEA_GET_OPERATION(n) == ALEA_OP_UNION && *size < 512) {
            ++*size;
            n = sys->nodes.data + n->operation.left;
        }
        return id;
    }
    if (!depth || ALEA_GET_OPERATION(n) != ALEA_OP_INTERSECTION) return UINT32_MAX;
    unsigned left_size = 0, right_size = 0;
    uint32_t left = render_union(sys, n->operation.left, depth-1, &left_size);
    uint32_t right = render_union(sys, n->operation.right, depth-1, &right_size);
    *size = left_size > right_size ? left_size : right_size;
    return left_size > right_size ? left : right;
}

static bool add_render_bounds(interval_index_t* idx, uint32_t root, uint32_t id, interval_box_t bounds,
                              unsigned depth, unsigned* remaining) {
    for (;;) {
        if (box_empty(bounds)) return true;
        alea_bbox_t box = hierarchy_box(bounds);
        if (alea_bbox_contains(&idx->domain, &box)) {
            idx->render_core = alea_bbox_union(&idx->render_core, &box);
            return true;
        }
        unsigned size = 0;
        uint32_t split = *remaining && depth < 32 ? render_union(idx->sys, id, 12, &size) : UINT32_MAX;
        if (split == UINT32_MAX) {
            alea_bbox_t lp;
            if (ALEA_GET_OPERATION(idx->sys->nodes.data + id) == ALEA_OP_INTERSECTION &&
                !alea_tighten_bbox_mixed_constraints(idx->sys, id, &lp)) {
                box = alea_bbox_intersection(&box, &lp);
                if (box_empty(interval_box_from(box))) return true;
                if (alea_bbox_contains(&idx->domain, &box)) {
                    idx->render_core = alea_bbox_union(&idx->render_core, &box);
                    return true;
                }
            }
            render_bound_t stored = {.branch=id, .cell_root=root};
            alea_node_bbox_set(&stored.box, &box);
            return !alea_vec_push(&idx->render_boxes, stored, render_bound_t);
        }
        --*remaining;
        const alea_node_t* n = idx->sys->nodes.data + split;
        unsigned budget = 256;
        interval_box_t right = contract_cell(idx->sys, n->operation.right, bounds, false, 0, &budget);
        if (!add_render_bounds(idx, root, n->operation.right, right, depth+1, remaining)) return false;
        id = n->operation.left;
        /* Continue down a left-deep union with its original enclosing box. */
        if (ALEA_GET_OPERATION(idx->sys->nodes.data + id) != ALEA_OP_UNION) {
            budget = 256;
            bounds = contract_cell(idx->sys, id, bounds, false, 0, &budget);
        }
    }
}

/* These cells take precedence over every potentially visible world cell.
 * A filled or lattice cell can contain material even when its own ID is void. */
static void prepare_render_void_prefix(interval_index_t* idx, const alea_universe_t* world) {
    idx->render_void_prefix_count = 0;
    idx->render_exterior_void_root = UINT32_MAX;
    for (size_t i = 0; i < world->cell_indices.count; ++i) {
        const alea_cell_entry_t* cell = idx->sys->cells.data + world->cell_indices.data[i];
        if (cell->material_id || cell->lat_type || alea_cell_entry_is_container(cell) ||
            cell->root_node_id >= idx->sys->nodes.count) break;
        ++idx->render_void_prefix_count;
        const alea_node_t* root = idx->sys->nodes.data + cell->root_node_id;
        if (ALEA_GET_OPERATION(root) != ALEA_OP_PRIMITIVE ||
            root->primitive.prim_type != ALEA_PRIMITIVE_SPHERE) continue;
        bool outside = (root->primitive.sense > 0) != (root->primitive.inverted != 0);
        const alea_sphere_data_t* s = alea_primitive_payload_const(idx->sys, root->primitive.primitive_id);
        if (outside && s && isfinite(s->center_x) && isfinite(s->center_y) &&
            isfinite(s->center_z) && isfinite(s->radius * s->radius))
            idx->render_exterior_void_root = cell->root_node_id;
    }
}

/* A world cell can contain visible material directly or through a fill. A
 * lattice can extend past its card's CSG bounds, so leave it to the walker. */
int alea_regional_prepare_render_bounds(interval_index_t* idx) {
    /* Called by the exclusive render consumer, before dispatching tile workers. */
    if (idx->render_bounds_attempted || !idx->domains) return 0;
    const alea_universe_t* world = alea_get_universe(idx->sys, 0);
    if (!world) return 0;
    idx->render_core = alea_bbox_empty();
    prepare_render_void_prefix(idx, world);
    for (size_t j = 0; j < world->cell_indices.count; ++j) {
        uint32_t id = world->cell_indices.data[j];
        const alea_cell_entry_t* cell = idx->sys->cells.data + id;
        if (!cell->material_id && !alea_cell_entry_is_container(cell)) continue;
        if (alea_interrupted()) {
            alea_set_error_detail(ALEA_ERR_INTERRUPTED, "Rendering interrupted");
            goto fail;
        }
        if (cell->lat_type || cell->root_node_id >= idx->sys->nodes.count) {
            alea_vec_free(&idx->render_boxes);
            idx->render_bounds_attempted = true;
            return 0;
        }
        alea_bbox_t raw = alea_node_bbox_get(&idx->sys->nodes.data[cell->root_node_id].bbox);
        interval_box_t bounds = interval_box_from(raw);
        if (!alea_bbox_contains(&idx->domain, &raw))
            bounds = cell_domain(idx->sys, id, raw);
        unsigned remaining = 4096;
        if (!add_render_bounds(idx, cell->root_node_id, cell->root_node_id, bounds, 0, &remaining)) goto fail;
    }
    idx->render_bounds_attempted = true;
    idx->render_bounds_ready = true;
    return 0;
fail:
    alea_vec_free(&idx->render_boxes);
    return -1;
}

/* Keep the inner resolution when a ray enters from outside the requested
 * box. Six disjoint slabs extend the world index by its width on each side. */
static int build_world_exterior(interval_index_t* idx,const alea_universe_t* world) {
    interval_box_t outer=interval_box_from(idx->domain);
    size_t n=world->cell_indices.count;
    uint32_t* ids=malloc(n*sizeof(*ids));
    if(!ids)return -1;
    for(unsigned a=0;a<3;a++) {
        double width=outer.hi[a]-outer.lo[a];
        if(!isfinite(outer.lo[a]-width) || !isfinite(outer.hi[a]+width)){free(ids);return 0;}
        for(unsigned side=0;side<2;side++) {
            interval_box_t slab=outer;
            if(side){slab.lo[a]=outer.hi[a];slab.hi[a]+=width;}
            else {slab.hi[a]=outer.lo[a];slab.lo[a]-=width;}
            struct exterior_region* region=&idx->exterior[2*a+side];
            region->box=(alea_bbox_t){.min_x=slab.lo[0],.max_x=slab.hi[0],
                .min_y=slab.lo[1],.max_y=slab.hi[1],.min_z=slab.lo[2],.max_z=slab.hi[2]};
            size_t count=0;
            for(size_t i=0;i<n;i++) {
                uint32_t id=world->cell_indices.data[i];
                idx->bounds[id]=bounded_cell_domain(idx,id,region->box);
                if(!box_empty(idx->bounds[id]))ids[count++]=id;
            }
            size_t budget = idx->points_only ? 64 + count * 4 : 1024 + count * 16;
            if(build_region(idx,region->box,ids,count,0,budget,&region->root)){
                free(ids);return -1;
            }
        }
        outer.lo[a]-=width;outer.hi[a]+=width;
    }
    idx->envelope=(alea_bbox_t){.min_x=outer.lo[0],.max_x=outer.hi[0],.min_y=outer.lo[1],
        .max_y=outer.hi[1],.min_z=outer.lo[2],.max_z=outer.hi[2]};
    idx->exterior_ready=true;free(ids);return 0;
}

/* Start point traversal below the upper split planes. Each table entry
 * covers its entire tile, including a pad for rounded grid coordinates. */
enum { POINT_GRID_WIDTH=16 };
static int build_point_starts(interval_index_t* idx) {
    size_t n=idx->sys->universes.count;
    idx->point_grids=calloc(n+1,sizeof(*idx->point_grids));
    if(!idx->point_grids)return -1;
    for(size_t i=0;i<n;i++) {
        uint32_t root=idx->universe_roots[i];
        if(root==UINT32_MAX || !idx->nodes.data[root].right ||
           idx->sys->universes.data[i].cell_indices.count<64)continue;
        point_grid_t* grid=idx->point_grids+i;
        interval_box_t domain=interval_box_from(idx->domains[i]);
        double step[3],margin[3];
        bool valid=true;
        for(unsigned k=0;k<3;k++) {
            grid->low[k]=domain.lo[k];
            step[k]=(domain.hi[k]-domain.lo[k])/POINT_GRID_WIDTH;
            grid->scale[k]=1.0/step[k];
            margin[k]=128*DBL_EPSILON*(1+fabs(domain.lo[k])+fabs(domain.hi[k]));
            if(!isfinite(grid->scale[k]) || step[k]<=margin[k])valid=false;
        }
        if(!valid)continue;
        size_t tiles=POINT_GRID_WIDTH*POINT_GRID_WIDTH*POINT_GRID_WIDTH;
        if(idx->point_nodes.count>UINT32_MAX-tiles ||
           alea_vec_reserve(&idx->point_nodes,idx->point_nodes.count+tiles,uint32_t))return -1;
        grid->first=(uint32_t)idx->point_nodes.count;
        for(unsigned tile=0;tile<tiles;tile++) {
            interval_box_t box;
            for(unsigned k=0;k<3;k++) {
                unsigned axis=(tile>>(4*k))&15u;
                box.lo[k]=domain.lo[k]+axis*step[k]-margin[k];
                box.hi[k]=domain.lo[k]+(axis+1)*step[k]+margin[k];
            }
            uint32_t id=root;
            while(idx->nodes.data[id].right) {
                spatial_node_t node=idx->nodes.data[id];
                if(box.hi[node.axis]<=node.split)id=node.first;
                else if(box.lo[node.axis]>node.split)id=node.right;
                else break;
            }
            idx->point_nodes.data[idx->point_nodes.count++]=id;
        }
        grid->ready=true;
    }
    return 0;
}


static void free_worker(regional_worker_t* worker) {
    free(worker->roots);
    free(worker->visited);
    free(worker->plane_gates);
    free(worker->quadric_gates);
    alea_regional_free_scratch(worker->levels);
}

static int init_worker(regional_worker_t* w, const interval_index_t* idx) {
    size_t roots = 1;
    while (roots < idx->sys->primitives.count && roots < 4096) roots *= 2;
    *w = (regional_worker_t){.index = idx, .root_mask = roots - 1};
    w->roots = calloc(roots, sizeof(*w->roots));
    if (idx->domains)
        w->visited = calloc(idx->sys->cells.count + 1, sizeof(*w->visited));
    w->plane_gates = calloc(idx->plane_count + 1, sizeof(*w->plane_gates));
    w->quadric_gates = calloc(idx->quadric_count + 1, sizeof(*w->quadric_gates));
    return w->roots && (!idx->domains || w->visited) && w->plane_gates && w->quadric_gates ? 0 : -1;
}

regional_worker_t* alea_regional_worker_create(const interval_index_t* idx) {
    regional_worker_t* w = malloc(sizeof(*w));
    if (w && init_worker(w, idx) == 0) return w;
    alea_regional_worker_destroy(w);
    alea_set_error_detail(ALEA_ERR_OUT_OF_MEMORY, "Allocating regional batch worker");
    return NULL;
}

void alea_regional_worker_destroy(regional_worker_t* w) {
    if (w) free_worker(w);
    free(w);
}

int alea_regional_reserve_workers(interval_index_t* idx, size_t count) {
    if (count <= idx->worker_count) return 0;
    if (count > SIZE_MAX / sizeof(*idx->workers)) {
        alea_set_error_detail(ALEA_ERR_OVERFLOW, "Regional worker count exceeds allocation limit");
        return -1;
    }
    regional_worker_t* workers = realloc(idx->workers, count * sizeof(*workers));
    if (!workers) goto fail;
    idx->workers = workers;
    for (size_t i = idx->worker_count; i < count; ++i) {
        if (init_worker(&workers[i], idx)) {
            for (size_t j = idx->worker_count; j <= i; ++j) free_worker(&workers[j]);
            goto fail;
        }
    }
    idx->worker_count = count;
    return 0;
fail:
    alea_set_error_detail(ALEA_ERR_OUT_OF_MEMORY, "Allocating regional ray worker storage");
    return -1;
}

regional_worker_t* alea_regional_worker(interval_index_t* idx, size_t slot) {
    return &idx->workers[slot];
}

void alea_prepared_queries_destroy(interval_index_t* idx) {
    if (!idx || atomic_fetch_sub(&idx->references, 1) != 1) return;
    if (idx->detach_mapping) idx->detach_mapping(idx);
    free(idx->box_ids);
    alea_vec_free(&idx->residual);
    free(idx->guards.bound_ids);free(idx->guards.certified);
    free(idx->proof_table);free(idx->primitive_proofs);
    alea_vec_free(&idx->tiles);alea_vec_free(&idx->point_items);
    free(idx->point_grids);alea_vec_free(&idx->point_nodes);
    for (size_t i = 0; i < idx->scan_count; ++i) {
        free(idx->scans[i].bounds);
        alea_prepared_queries_destroy(idx->scans[i].region);
        child_bound_block_t* block = idx->scans[i].first;
        while (block) {
            child_bound_block_t* next = block->next;
            free(block);
            block = next;
        }
    }
    alea_vec_free(&idx->render_boxes);
    free(idx->scans);
    alea_vec_free(&idx->exterior_point_cells);
    free(idx->centers);free(idx->order);
    for (size_t i = 0; i < idx->worker_count; ++i) free_worker(&idx->workers[i]);
    free(idx->workers);
    free(idx->active_shapes);
    free(idx->planes);
    free(idx->quadrics);
    free(idx->shape_id);
    free(idx->cell_cards);
    alea_vec_free(&idx->boxes);
    alea_vec_free(&idx->surface_cards);
    alea_vec_free(&idx->patch_cards);
    alea_vec_free(&idx->decisions);alea_vec_free(&idx->wide);
    alea_vec_free(&idx->extra_axes);
    alea_vec_free(&idx->packets);free(idx->guards.table);
    alea_vec_free(&idx->guards.pairs);
    alea_vec_free(&idx->cards);alea_vec_free(&idx->patches);
    free(idx->surfaces);
    free(idx->bounds);
    free(idx->card_info);
    free(idx->domains);free(idx->universe_roots);
    alea_vec_free(&idx->nodes);
    free(idx);
}

interval_fallback_stats_t interval_index_fallback_stats(const interval_index_t* idx) {
    return idx->workers[0].fallback;
}

bool interval_index_is_regional(const interval_index_t* idx) {
    return idx->nodes.count!=0;
}

/* Allocated model index and per-primitive working tables. Per-depth temporary
 * ray buffers grow with the query and are included in the process peak RSS. */
size_t interval_index_bytes(const interval_index_t* idx) {
    size_t cells = idx->sys->cells.count + 1;
    size_t bytes = sizeof(*idx) + (idx->sys->universes.count+1)*sizeof(*idx->universe_roots);
    if (idx->child_mode == CHILD_POINTS_ADAPTIVE) {
        for (size_t i = 0; i < idx->scan_count; ++i) {
            size_t cells = atomic_load(&idx->scans[i].covered);
            size_t blocks = cells / CHILD_BOUND_BLOCK_SIZE + (cells % CHILD_BOUND_BLOCK_SIZE != 0);
            bytes += cells * sizeof(alea_node_bbox_t) + blocks * sizeof(child_bound_block_t);
            const interval_index_t* child = published_child_region(&idx->scans[i]);
            if (child) bytes += interval_index_bytes(child);
        }
    }
    bytes += idx->child_bounds_count * sizeof(alea_node_bbox_t);
    bytes += idx->exterior_point_cells.capacity * sizeof(*idx->exterior_point_cells.data);
    if (idx->scans) bytes += idx->sys->universes.count * sizeof(*idx->scans);
    for (size_t i = 0; i < idx->worker_count; ++i) {
        const regional_worker_t* w = &idx->workers[i];
        bytes += sizeof(*w) + (w->root_mask + 1) * sizeof(*w->roots)
            + (idx->plane_count + idx->quadric_count + 2) * sizeof(ray_gate_t);
        if (w->visited) bytes += cells * sizeof(*w->visited);
    }
    bytes+=idx->tiles.capacity*sizeof(point_tiles_t)+idx->point_items.capacity*sizeof(point_item_t);
    if(idx->point_grids)bytes+=(idx->sys->universes.count+1)*sizeof(*idx->point_grids)
        +idx->point_nodes.capacity*sizeof(*idx->point_nodes.data);
    bytes += idx->render_boxes.capacity * sizeof(render_bound_t);
    if(idx->domains)bytes+=(idx->sys->universes.count+1)*sizeof(*idx->domains);
    bytes+=idx->extra_axes.capacity*sizeof(packet_axis_t);
    bytes+=idx->decisions.capacity*sizeof(*idx->decisions.data);
    if (idx->planes) bytes+=(idx->plane_count+1)*sizeof(*idx->planes);
    if (idx->quadrics) bytes+=(idx->quadric_count+1)*sizeof(*idx->quadrics);
    if (idx->cell_cards) bytes+=cells*sizeof(*idx->cell_cards);
    bytes+=idx->boxes.capacity*sizeof(*idx->boxes.data)+
        idx->patch_cards.capacity*sizeof(uint32_t)+idx->surface_cards.capacity*sizeof(query_card_t)+idx->wide.capacity*sizeof(query_decision_t);
    return bytes+idx->packets.capacity*sizeof(*idx->packets.data)+idx->patches.capacity*sizeof(*idx->patches.data)+
        idx->cards.capacity*sizeof(*idx->cards.data)+
        idx->nodes.capacity*sizeof(*idx->nodes.data)+
        idx->guards.pairs.capacity*sizeof(*idx->guards.pairs.data);
}

/* Outside the indexed world box, only lattices and cells whose original
 * root bounds cross the box can own a point. Preserve definition order. */
static int collect_exterior_point_cells(interval_index_t* idx, const alea_universe_t* world) {
    for (size_t i = 0; i < world->cell_indices.count; ++i) {
        size_t id = world->cell_indices.data[i];
        const alea_cell_entry_t* cell = &idx->sys->cells.data[id];
        if (!(cell->lat_type && cell->lat_fill)) {
            if (cell->root_node_id >= idx->sys->nodes.count) continue;
            alea_bbox_t box = alea_node_bbox_get(&idx->sys->nodes.data[cell->root_node_id].bbox);
            if (alea_bbox_contains(&idx->domain, &box)) continue;
        }
        if (alea_vec_push(&idx->exterior_point_cells, id, size_t)) return -1;
    }
    return 0;
}

typedef enum {
    INDEX_POINTS_AND_RAYS,
    INDEX_ALL_UNIVERSE_POINTS,
    INDEX_ROOT_POINTS,
    INDEX_CHILD_POINTS
} index_scope_t;

static interval_index_t* build_index(alea_system_t* sys, const alea_bbox_t* domain,
                                     bool selective, index_scope_t scope,
                                     const interval_point_children_t* children) {
    const bool points_only = scope != INDEX_POINTS_AND_RAYS;
    const alea_universe_t* only_child = scope == INDEX_CHILD_POINTS
        ? alea_get_universe(sys, children->ids[0]) : NULL;
    const bool root_only = scope == INDEX_ROOT_POINTS || scope == INDEX_CHILD_POINTS;
    interval_index_t* idx = alea_calloc(1, sizeof(*idx));
    if (!idx) {
        alea_set_error_detail(ALEA_ERR_OUT_OF_MEMORY, "Allocating regional index");
        return NULL;
    }
    double stage = alea_load_profile_start();
    idx->points_only = points_only;
    idx->root_only = root_only;
    atomic_init(&idx->references, 1);
    idx->sys = sys;
    idx->geometry_generation = alea_system_geometry_generation(sys);
    idx->material_generation = alea_system_material_generation(sys);
    idx->selective = selective;
    if (root_only && sys->universes.count) {
        idx->scans = calloc(sys->universes.count, sizeof(*idx->scans));
        if (!idx->scans) goto fail;
        idx->scan_count = sys->universes.count;
        const char* timing = getenv("ALEA_PROFILE_CHILD_POINTS");
        idx->time_child_points = timing && *timing && strcmp(timing, "0");
        for (size_t i = 0; i < idx->scan_count; ++i) {
            atomic_init(&idx->scans[i].local_ns, 0);
            atomic_init(&idx->scans[i].unfiltered_ns, 0);
            atomic_init(&idx->scans[i].build_ns, 0);
            atomic_init(&idx->scans[i].covered, 0);
            atomic_init(&idx->scans[i].building, 0);
            atomic_init(&idx->scans[i].region_state, CHILD_REGION_ABSENT);
            atomic_init(&idx->scans[i].tests, 0);
            atomic_init(&idx->scans[i].visits, 0);
        }
        if (children) {
            idx->child_count = children->count;
            idx->child_mode = children->mode;
            for (size_t i = 0; i < children->count; ++i) {
                const alea_universe_t* u = alea_get_universe(sys, children->ids[i]);
                idx->scans[u - sys->universes.data].indexed = true;
            }
        }
    }
    if (!points_only) {
        idx->card_info = calloc(sys->surfaces.count + 1, sizeof(*idx->card_info));
        if (!idx->card_info) goto fail;
    }
    idx->universe_roots = calloc(sys->universes.count + 1, sizeof(uint32_t));
    if (!idx->universe_roots) goto fail;
    /* A promoted child shares the model already validated by its parent. */
    for (size_t i=0;!only_child && i<sys->surfaces.count;i++) {
        const alea_surface_entry_t* s=&sys->surfaces.data[i];
        if(s->primitive_id>=sys->primitives.count)goto fail;
        if (!points_only) idx->card_info[i]=(struct interval_card){s->primitive_id,s->mc_surface_id};
    }
    for(size_t i=0;!only_child && i<sys->cells.count;i++) {
        const alea_cell_entry_t* c=sys->cells.data+i;
        for(size_t j=0;j<c->surface_index_count;j++)
            if(c->surface_indices[j]>=sys->surfaces.count)goto fail;
    }
    if(sys->primitives.count >= (QUERY_PLANE-2)/2 || sys->nodes.count>=GUARD_BOX/2)domain=NULL;
    if(domain) {
        interval_box_t b=interval_box_from(*domain);
        for(unsigned k=0;k<3;k++)
            if(!(b.lo[k]<b.hi[k]) || !isfinite(b.hi[k]-b.lo[k]))goto fail;
    }
    if (domain) {
        idx->proof_table=calloc(65536,sizeof(*idx->proof_table));
        idx->primitive_proofs=calloc(sys->primitives.count+1,sizeof(*idx->primitive_proofs));
        if(!idx->proof_table || !idx->primitive_proofs)goto fail;
        if (!points_only) {
            idx->active_shapes=calloc(sys->primitives.count+1,sizeof(*idx->active_shapes));
            if(!idx->active_shapes)goto fail;
        }
        idx->domain = *domain;
        idx->domains=calloc(sys->universes.count+1,sizeof(*idx->domains));
        if(!idx->domains || prepare_root_bounds(idx, only_child))goto fail;
        if (root_only) {
            /* Use the same finite numeric box in each indexed local frame.
             * Points outside it retain the ordered local scan. */
            for (size_t i = 0; i < sys->universes.count; ++i)
                idx->domains[i] = *domain;
        } else if (choose_domains(idx)) goto fail;
        idx->proof_margin=2048*DBL_EPSILON*(1+fabs(domain->min_x)+fabs(domain->max_x)+
            fabs(domain->min_y)+fabs(domain->max_y)+fabs(domain->min_z)+fabs(domain->max_z));
        idx->surfaces=calloc(sys->primitives.count+1,sizeof(guard_surface_t));
        if (!idx->surfaces) goto fail;
        for (size_t i=0;i<sys->primitives.count;i++) guard_surface_init(&idx->surfaces[i],sys,(uint32_t)i,true);
        idx->bounds = calloc(sys->cells.count + 1, sizeof(interval_box_t));
        if (!idx->bounds) goto fail;
    }
    alea_load_profile_end("regional_setup_bounds", stage);
    stage = alea_load_profile_start();
    if (!root_only && pack_shapes(idx)) goto fail;
    alea_load_profile_end("regional_pack", stage);
    stage = alea_load_profile_start();
    idx->box_ids=calloc(sys->nodes.count+1,sizeof(*idx->box_ids));
    if (!idx->box_ids) goto fail;
    for (size_t i = 0; i < sys->universes.count; i++) {
        idx->universe_roots[i] = UINT32_MAX;
        const alea_universe_t* u = &sys->universes.data[i];
        size_t n = u->cell_indices.count;
        if (only_child && u != only_child) continue;
        if (!domain || !n || (root_only && u->universe_id != 0 && !idx->scans[i].indexed)) continue;
        double region_start = root_only && u->universe_id != 0 ? alea_load_profile_start() : -1;
        cell_bounds_work_t work = {idx, u, idx->domains[i]};
        double bounds_start = root_only ? alea_load_profile_start() : -1;
        if (root_only && n >= 128) {
            if (alea_parallel_for_independent(n, 32, prepare_cell_bounds, &work) != ALEA_PARALLEL_OK) {
                goto fail;
            }
        } else {
            prepare_cell_bounds(&work, 0, 0, n);
        }
        if (root_only) alea_load_profile_end("regional_cell_bounds", bounds_start);
        if (root_only && u->universe_id != 0 && idx->child_mode == CHILD_POINTS_BOUNDS) {
            alea_node_bbox_t* bounds = calloc(n, sizeof(*bounds));
            if (!bounds) goto fail;
            idx->scans[i].bounds = bounds;
            idx->child_bounds_count += n;
            for (size_t j = 0; j < n; ++j) {
                uint32_t id = (uint32_t)u->cell_indices.data[j];
                alea_bbox_t box = hierarchy_box(idx->bounds[id]);
                alea_node_bbox_set(&bounds[j], &box);
            }
            alea_load_profile_end("regional_point_child", region_start);
            continue;
        }
        uint32_t* ids = malloc(n * sizeof(*ids));
        if (!ids) goto fail;
        size_t count = 0;
        for (size_t j = 0; j < n; j++) {
            uint32_t id = (uint32_t)u->cell_indices.data[j];
            if (!box_empty(idx->bounds[id])) ids[count++] = id;
        }
        /* Every global query starts in universe zero. Dense collections of
         * short cell expressions also repay finer subdivision cheaply. */
        uint64_t cards = 0;
        for (size_t j = 0; j < count; ++j)
            cards += sys->cells.data[ids[j]].surface_index_count;
        bool dense = count >= 512 && cards / count < 16;
        size_t budget = 512 + count * 24;
        idx->leaf_size = 12;
        if (u->universe_id == 0) {
            budget = 8192 + count * 384;
            idx->leaf_size = 8;
        } else if (dense) {
            budget = 1024 + count * 100;
            idx->leaf_size = 8;
        }
        if (points_only) {
            /* Every lookup pays for the world search. Spatial duplication in
             * child universes only benefits queries that descend into them. */
            bool world = u->universe_id == 0;
            bool fine_root = world && !root_only;
            budget = 64 + count * (fine_root ? 32 : 4);
            idx->leaf_size = fine_root ? 12 : 24;
        }
        double subdivision_start = region_start >= 0 ? alea_monotonic_seconds() : -1;
        int rc = build_region(idx, idx->domains[i], ids, count, 0, budget, &idx->universe_roots[i]);
        alea_load_profile_end("regional_child_subdivision", subdivision_start);
        free(ids);
        double programs_start = region_start >= 0 ? alea_monotonic_seconds() : -1;
        if (rc || finish_region_programs(idx)) goto fail;
        alea_load_profile_end("regional_child_finish_programs", programs_start);
        alea_load_profile_end("regional_point_child", region_start);
    }
    alea_load_profile_end("regional_subdivide_compile", stage);
    stage = alea_load_profile_start();
    const alea_universe_t* world=alea_get_universe(sys,0);
    if (domain && world && idx->universe_roots[world - sys->universes.data] != UINT32_MAX) {
        if (root_only) {
            if (collect_exterior_point_cells(idx, world)) goto fail;
        } else if (build_world_exterior(idx, world) || finish_region_programs(idx)) goto fail;
    }
    alea_load_profile_end("regional_exterior", stage);
    stage = alea_load_profile_start();
    free(idx->bounds);idx->bounds=NULL;
    free(idx->surfaces);idx->surfaces=NULL;
    free(idx->guards.table);idx->guards.table=NULL;
    free(idx->active_shapes);idx->active_shapes=NULL;
    free(idx->proof_table);idx->proof_table=NULL;
    free(idx->primitive_proofs);idx->primitive_proofs=NULL;
    free(idx->guards.bound_ids);idx->guards.bound_ids=NULL;
    free(idx->guards.certified);idx->guards.certified=NULL;
    if(finish_index(idx) || (domain && build_point_starts(idx)))goto fail;
    if (!points_only && alea_regional_reserve_workers(idx, 1)) goto fail;
    alea_load_profile_end("regional_finalize", stage);
    return idx;
fail:
    alea_prepared_queries_destroy(idx);
    return NULL;
}

interval_index_t* interval_index_build(alea_system_t* sys, const alea_bbox_t* domain, bool selective) {
    return build_index(sys, domain, selective, INDEX_POINTS_AND_RAYS, NULL);
}

static bool selected_world_domain(alea_system_t* sys, alea_bbox_t* box) {
    const alea_universe_t* world = alea_get_universe(sys, 0);
    const alea_bbox_t infinite = alea_bbox_infinite();
    *box = alea_bbox_empty();
    for (size_t i = 0; world && i < world->cell_indices.count; ++i) {
        const alea_cell_entry_t* cell = &sys->cells.data[world->cell_indices.data[i]];
        if (cell->root_node_id >= sys->nodes.count) continue;
        alea_bbox_t b = alea_node_bbox_get(&sys->nodes.data[cell->root_node_id].bbox);
        if (alea_bbox_is_valid(&b) && b.min_x > infinite.min_x && b.max_x < infinite.max_x &&
            b.min_y > infinite.min_y && b.max_y < infinite.max_y &&
            b.min_z > infinite.min_z && b.max_z < infinite.max_z)
            *box = alea_bbox_union(box, &b);
    }
    return box->min_x < box->max_x && box->min_y < box->max_y && box->min_z < box->max_z;
}

static bool same_domain(alea_bbox_t a, alea_bbox_t b) {
    return a.min_x == b.min_x && a.max_x == b.max_x &&
           a.min_y == b.min_y && a.max_y == b.max_y &&
           a.min_z == b.min_z && a.max_z == b.max_z;
}

int alea_prepare_ray_acceleration(alea_system_t* sys, const alea_bbox_t* requested) {
    if (alea_system_prepare_query_caches(sys, ALEA_CACHE_HIER_SPATIAL)) return -1;
    alea_bbox_t domain;
    if (sys->ray_batch_index && (!requested || same_domain(*requested, sys->ray_batch_index->domain)))
        return 0;
    if (!requested && !selected_world_domain(sys, &domain)) return 0;
    interval_index_t* next = alea_prepared_queries_create(sys, requested ? requested : &domain);
    if (!next) return -1;
    int rc = alea_set_ray_batch_accelerator(sys, next);
    alea_prepared_queries_destroy(next);
    return rc;
}

static int prepare_selected_points(alea_system_t* sys, const alea_bbox_t* requested,
                                   bool root_only, const interval_point_children_t* children, bool automatic) {
    interval_child_mode_t mode = children ? children->mode : CHILD_POINTS_BOUNDS;
    if (children && children->count && !children->ids) {
        alea_set_error_detail(ALEA_ERR_INVALID_ARG, "Missing point child universe IDs");
        return -1;
    }
    if (requested) {
        interval_box_t bounds = interval_box_from(*requested);
        for (unsigned i = 0; i < 3; ++i) {
            if (!(bounds.lo[i] < bounds.hi[i]) || !isfinite(bounds.hi[i] - bounds.lo[i])) {
                alea_set_error_detail(ALEA_ERR_INVALID_ARG, "Invalid selected query domain");
                return -1;
            }
        }
    }
    alea_clear_error_detail();
    unsigned prerequisites = root_only
        ? ALEA_CACHE_UNIVERSE | ALEA_CACHE_CELL_SURFACES | ALEA_CACHE_FILL_MATRICES
        : ALEA_CACHE_HIER_SPATIAL;
    if (alea_system_prepare_query_caches(sys, prerequisites)) return -1;
    while (atomic_flag_test_and_set(&sys->query_cache_build_lock)) {
        if (alea_interrupted()) {
            alea_set_error_detail(ALEA_ERR_INTERRUPTED, "Selected query preparation interrupted");
            return -1;
        }
    }
    /* Explicit publication can supersede the automatic index. Keep the old
     * allocation alive for existing readers until the next exclusive clear. */
    if (!sys->selected_index && sys->automatic_point_index && root_only && children &&
        children->mode == CHILD_POINTS_ADAPTIVE && !children->count && !automatic) {
        /* Explicit preparation claims the same allocation and scope. */
        atomic_fetch_add(&sys->automatic_point_index->references, 1);
        sys->selected_index = sys->automatic_point_index;
        atomic_fetch_or(&sys->query_cache_state, ALEA_CACHE_SELECTED_POINTS);
    }
    if (automatic && (sys->automatic_point_index || sys->selected_index)) {
        atomic_flag_clear(&sys->query_cache_build_lock);
        return 0;
    }
    int rc = 0;
    double start = alea_monotonic_seconds();
    alea_bbox_t domain;
    size_t child_count = children ? children->count : 0;
    bool same_children = !sys->selected_index || sys->selected_index->child_count == child_count;
    for (size_t i = 0; i < child_count; ++i) {
        int id = children->ids[i];
        const alea_universe_t* u = alea_get_universe(sys, id);
        if (id <= 0 || !u || (i && id <= children->ids[i - 1])) {
            alea_set_error_detail(ALEA_ERR_INVALID_ARG, "Expected distinct existing child universe IDs in increasing order (got %d)", id);
            atomic_flag_clear(&sys->query_cache_build_lock);
            return -1;
        }
        if (sys->selected_index && (!sys->selected_index->scans ||
            !sys->selected_index->scans[u - sys->universes.data].indexed)) same_children = false;
    }
    if (sys->selected_index) {
        if (!same_children || mode != sys->selected_index->child_mode ||
            root_only != sys->selected_index->root_only ||
            (requested && !same_domain(*requested, sys->selected_index->domain))) {
            alea_set_error_detail(ALEA_ERR_INVALID_STATE, "Clear the selected query index before changing its domain or scope");
            rc = -1;
        }
    } else if (requested || selected_world_domain(sys, &domain)) {
        interval_index_t* index = build_index(sys, requested ? requested : &domain, true,
            root_only ? INDEX_ROOT_POINTS : INDEX_ALL_UNIVERSE_POINTS, children);
        if (!index || alea_interrupted()) {
            alea_prepared_queries_destroy(index);
            if (alea_interrupted())
                alea_set_error_detail(ALEA_ERR_INTERRUPTED, "Selected query preparation interrupted");
            else if (alea_get_last_error() == ALEA_OK)
                alea_set_error_detail(ALEA_ERR_INVALID_STATE, "Could not build selected query index");
            rc = -1;
        } else if (interval_index_is_regional(index)) {
            index->build_seconds = alea_monotonic_seconds() - start;
            if (automatic) sys->automatic_point_index = index;
            else sys->selected_index = index;
            atomic_store(&sys->selected_fallbacks, 0);
            atomic_fetch_or(&sys->query_cache_state, automatic ? ALEA_CACHE_AUTO_POINTS : ALEA_CACHE_SELECTED_POINTS);
        } else {
            alea_prepared_queries_destroy(index);
        }
    }
    if (child_count && !rc && !sys->selected_index) {
        alea_set_error_detail(ALEA_ERR_INVALID_STATE, "Child point indexes require an indexed world domain");
        rc = -1;
    }
    atomic_flag_clear(&sys->query_cache_build_lock);
    return rc;
}

int alea_prepare_selected_queries(alea_system_t* sys, const alea_bbox_t* domain) {
    return prepare_selected_points(sys, domain, false, NULL, false);
}

int alea_prepare_point_acceleration(alea_system_t* sys) {
    const interval_point_children_t children = {.mode = CHILD_POINTS_ADAPTIVE};
    return prepare_selected_points(sys, NULL, true, &children, false);
}

int alea_prepare_automatic_points(alea_system_t* sys) {
    const interval_point_children_t children = {.mode = CHILD_POINTS_ADAPTIVE};
    return prepare_selected_points(sys, NULL, true, &children, true);
}

/* Internal first-build experiment for the lifecycle benchmark. */
int alea_prepare_root_points(alea_system_t* sys, const interval_point_children_t* children) {
    return prepare_selected_points(sys, NULL, true, children, false);
}

int alea_promote_point_child(alea_system_t* sys, int universe_id) {
    interval_index_t* idx = sys ? alea_system_selected_index(sys) : NULL;
    const alea_universe_t* u = idx ? alea_get_universe(sys, universe_id) : NULL;
    if (!idx || idx->child_mode != CHILD_POINTS_ADAPTIVE || universe_id <= 0 || !u) {
        alea_set_error_detail(ALEA_ERR_INVALID_ARG, "Expected an adaptive point index and existing child universe");
        return -1;
    }
    size_t i = (size_t)(u - sys->universes.data);
    if (idx->universe_roots[i] != UINT32_MAX) return 0;
    struct universe_scans* scans = &idx->scans[i];
    int expected = CHILD_REGION_ABSENT;
    if (!atomic_compare_exchange_strong(&scans->region_state, &expected, CHILD_REGION_BUILDING))
        return expected == CHILD_REGION_READY ? 0 : 1;
    const interval_point_children_t child = {&universe_id, 1, CHILD_POINTS_REGIONS};
    alea_clear_error_detail();
    interval_index_t* region = build_index(sys, &idx->domain, true, INDEX_CHILD_POINTS, &child);
    if (!region || region->universe_roots[i] == UINT32_MAX || alea_interrupted()) {
        alea_prepared_queries_destroy(region);
        if (alea_interrupted()) alea_set_error_detail(ALEA_ERR_INTERRUPTED, "Child point preparation interrupted");
        else if (alea_get_last_error() == ALEA_OK)
            alea_set_error_detail(ALEA_ERR_INVALID_STATE, "Could not build a region for the selected child");
        atomic_store_explicit(&scans->region_state, CHILD_REGION_ABSENT, memory_order_release);
        return -1;
    }
    scans->region = region;
    atomic_store_explicit(&scans->region_state, CHILD_REGION_READY, memory_order_release);
    return 0;
}

uint64_t alea_root_point_scan_tests(const alea_system_t* sys, size_t universe,
                                   uint64_t* visits, double* local_seconds) {
    if (local_seconds) *local_seconds = 0;
    const interval_index_t* idx = alea_system_selected_index(sys);
    if (!idx || !idx->scans || universe >= sys->universes.count) {
        *visits = 0;
        return 0;
    }
    if (local_seconds) *local_seconds = 1e-9 * atomic_load(&idx->scans[universe].local_ns);
    *visits = atomic_load(&idx->scans[universe].visits);
    return atomic_load(&idx->scans[universe].tests);
}

size_t alea_adaptive_child_progress(const alea_system_t* sys, size_t universe, double* seconds) {
    const interval_index_t* idx = alea_system_selected_index(sys);
    if (seconds) *seconds = 0;
    if (!idx || idx->child_mode != CHILD_POINTS_ADAPTIVE || universe >= idx->scan_count) return 0;
    size_t covered = atomic_load(&idx->scans[universe].covered);
    if (seconds) *seconds = atomic_load(&idx->scans[universe].build_ns) * 1e-9;
    return covered;
}

void alea_clear_selected_queries(alea_system_t* sys) {
    if (!sys) return;
    if (sys->selected_index || sys->automatic_point_index) {
        atomic_fetch_and(&sys->query_cache_state, ~(ALEA_CACHE_SELECTED_POINTS | ALEA_CACHE_AUTO_POINTS));
        alea_prepared_queries_destroy(sys->selected_index);
        alea_prepared_queries_destroy(sys->automatic_point_index);
        sys->selected_index = NULL;
        sys->automatic_point_index = NULL;
    }
    atomic_store(&sys->point_lookup_count, 0);
    atomic_store(&sys->point_sample_ns, 0);
    atomic_store(&sys->point_promotion_attempted, 0);
}

int alea_selected_query_stats(const alea_system_t* sys, alea_selected_query_stats_t* stats) {
    const interval_index_t* idx = sys ? alea_system_selected_index(sys) : NULL;
    if (stats) {
        *stats = (alea_selected_query_stats_t){0};
        if (idx) {
            stats->domain = idx->domain;
            stats->bytes = interval_index_bytes(idx);
            if (idx != sys->automatic_point_index && sys->automatic_point_index)
                stats->bytes += interval_index_bytes(sys->automatic_point_index);
            stats->build_seconds = idx->build_seconds;
            stats->fallback_visits = atomic_load(&sys->selected_fallbacks);
        }
    }
    return idx != NULL;
}
