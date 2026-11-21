// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
// SPDX-License-Identifier: MPL-2.0

/* Selected ownership along finite rays through nested CSG universes. */
#include "raycast/ray_traversal.h"
#include "raycast/ray_substep.h"
#include "raycast/ray_quadratic.h"
#include "raycast/ray_epsilon.h"
#include "raycast/ray_intersect.h"
#include "core/alea_eval.h"
#include "core/alea_universe.h"
#include "alea.h"
#include <string.h>
#include "interval_internal.h"
#include <float.h>
#include <limits.h>
#include <math.h>
#include "primitives/bbox.h"

typedef struct {
    double t;
    int surface;
    uint32_t primitive;
} interval_end_t;

typedef struct {
    interval_end_t lo, hi;
    uint32_t cell;
    alea_lattice_location_t element;
} cell_run_t;
ALEA_VEC_DEFINE(cell_runs, cell_run_t);
ALEA_VEC_DEFINE(interval_ends, interval_end_t);

typedef struct { double original, refined; uint32_t primitive; } root_correction_t;
ALEA_VEC_DEFINE(root_corrections, root_correction_t);

typedef struct {
    alea_raycast_render_accept_t accept;
    void* context;
    bool normal;
    alea_ray_first_visible_result_t* hit;
} regional_visibility_t;

typedef struct {
    regional_visibility_t* visibility;
    regional_nav_span_t* nav;
    alea_system_t* sys;
    const alea_ray_t* world;
    alea_raycast_result_t* out;
    int previous_cell;
    interval_end_t incident_exit;
    alea_matrix_t incident_frame; /* The child frame can unwind before the next entry. */
    const interval_index_t* index;
    regional_worker_t* worker;
} interval_query_t;

typedef struct { uint32_t cell,head; double lo,hi; } possible_cell_t;
ALEA_VEC_DEFINE(possible_cells,possible_cell_t);

typedef struct { uint32_t item,next; double lo,hi; } ray_patch_t;
ALEA_VEC_DEFINE(ray_patches,ray_patch_t);

typedef struct interval_level {
    double magnitude;
    double local_direction_length;
    double os[9],ds[9],iv[9],err[9];
    double origin[3],direction[3],inverse[3],pad[3];
    uint32_t head;
    ray_patches_t patches;
    possible_cells_t possible;
    bool specialized;
    root_corrections_t corrections;
    bool recovering_boundaries, has_inactive;
    bool shares_root_frame;
    interval_query_t* query;
    alea_ray_t local;
    interval_end_t lo, hi;
    cell_runs_t runs;
    interval_ends_t cuts;
    cell_ids_t heap;
} interval_level_t;

#include "interval_packets.h"
#include "interval_ray_clip.h"
#include "interval_quadric_refine.h"

static int compare_ends(const void* a, const void* b) {
    const interval_end_t *x = a, *y = b;
    if (x->t != y->t) return x->t < y->t ? -1 : 1;
    return (x->surface > y->surface) - (x->surface < y->surface);
}

/* Keep comparisons in the caller for the many small lists of roots and
 * runs. A depth limit bounds the cost for adversarial orderings. */
static void sort_end_range(interval_end_t* data,size_t n,unsigned depth) {
    while(n>20) {
        if(!depth--) {
            qsort(data,n,sizeof(*data),compare_ends);
            return;
        }
        interval_end_t x=data[0],y=data[n/2],z=data[n-1];
        if(compare_ends(&x,&y)>0){interval_end_t t=x;x=y;y=t;}
        if(compare_ends(&y,&z)>0)y=compare_ends(&x,&z)>0?x:z;
        size_t low=0,i=0,high=n;
        while(i<high) {
            int side=compare_ends(data+i,&y);
            if(side<0){interval_end_t t=data[i];data[i++]=data[low];data[low++]=t;}
            else if(side>0){interval_end_t t=data[i];data[i]=data[--high];data[high]=t;}
            else i++;
        }
        if(low<n-high) {
            sort_end_range(data,low,depth);
            data+=high;n-=high;
        }else {
            sort_end_range(data+high,n-high,depth);
            n=low;
        }
    }
    for(size_t i=1;i<n;i++) {
        interval_end_t end=data[i];
        size_t j=i;
        while(j && compare_ends(&end,data+j-1)<0){data[j]=data[j-1];j--;}
        data[j]=end;
    }
}
static void sort_ends(interval_ends_t* cuts) {
    sort_end_range(cuts->data,cuts->count,32);
}

void alea_regional_free_scratch(interval_level_t* levels) {
    if (!levels) return;
    for (int i=0;i<MAX_FILL_RAYCAST_DEPTH;i++) {
        alea_vec_free(&levels[i].patches);
        alea_vec_free(&levels[i].possible);
        alea_vec_free(&levels[i].runs);
        alea_vec_free(&levels[i].cuts);
        alea_vec_free(&levels[i].heap);
    }
    free(levels);
}
static double corrected_root(const interval_level_t* l, uint32_t primitive, double t) {
    size_t lo = 0, hi = l->corrections.count;
    while (lo < hi) {
        size_t m = lo + (hi - lo) / 2;
        if (l->corrections.data[m].original < t) lo = m + 1;
        else hi = m;
    }
    for (; lo < l->corrections.count && l->corrections.data[lo].original == t; ++lo)
        if (l->corrections.data[lo].primitive == primitive) return l->corrections.data[lo].refined;
    return t;
}

static int append_card_roots(interval_level_t* level,const uint32_t* cards,size_t count,
                             double active_lo,double active_hi) {
    interval_query_t* q=level->query;
    for (size_t i = 0; i < count; i++) {
        const interval_index_t* idx=q->index;
        query_card_t card=idx->surface_cards.data[cards[i]];
        uint32_t primitive = q->sys->surfaces.data[cards[i]].primitive_id;
        uint32_t id = card.shape & CARD_ID;
        regional_worker_t* worker = q->worker;
        cached_roots_t* g = &worker->roots[primitive & worker->root_mask];
        if (g->generation != worker->generation || g->primitive != primitive) {
            if (card.shape & CARD_PLANE) {
                g->count = ray_intersect_plane(&level->local, &idx->planes[id].data,
                                                g->t, NULL, NULL, NULL);
            } else if (card.shape & CARD_QUADRIC) {
                ray_gate_t* poly = ray_quadric_gate(level, id);
                g->count = solve_quadratic_stable(poly->a, poly->b, poly->c, g->t);
            } else {
                alea_primitive_type_t type = q->sys->primitives.data[id].type;
                const void* payload=alea_primitive_payload_const(q->sys,id);
                if (!payload) return -1;
                switch(type) {
                case ALEA_PRIMITIVE_CYLINDER_X:g->count=ray_intersect_cylinder_x(&level->local,payload,g->t);break;
                case ALEA_PRIMITIVE_CYLINDER_Y:g->count=ray_intersect_cylinder_y(&level->local,payload,g->t);break;
                case ALEA_PRIMITIVE_CYLINDER_Z:g->count=ray_intersect_cylinder_z(&level->local,payload,g->t);break;
                case ALEA_PRIMITIVE_SPHERE:case ALEA_PRIMITIVE_SPH:
                    g->count=ray_intersect_sphere(&level->local,payload,g->t);break;
                case ALEA_PRIMITIVE_TORUS_X:case ALEA_PRIMITIVE_TORUS_Y:case ALEA_PRIMITIVE_TORUS_Z:
                    if(level->specialized && !ray_torus_window_may_hit(&level->local,payload,active_lo,active_hi))
                        continue;
                    g->count=ray_intersect_torus(&level->local,payload,g->t);break;
                default: {
                    alea_primitive_data_t data;
                    if(!alea_primitive_copy_data(q->sys,id,&data))return -1;
                    g->count=ray_intersect_primitive(&level->local,type,&data,g->t);
                    break;
                }
                }
            }
            g->primitive = primitive;
            g->generation = worker->generation;
            q->out->surfaces_tested++;
        }
        for(int k=0;k<g->count;k++) {
            double t=g->t[k];
            if (level->corrections.count) t = corrected_root(level, primitive, t);
            if (!(t > level->lo.t && t < level->hi.t && t >= active_lo && t <= active_hi)) continue;
            interval_end_t end = {t, card.number, primitive};
            if (alea_vec_push(&level->cuts, end, interval_end_t)) return -1;
        }
    }
    return 0;
}

static int cell_roots(interval_level_t* level,uint32_t id) {
    interval_query_t* q = level->query;
    const interval_index_t* idx=q->index;
    const query_cell_t* span=idx->cell_cards+id;
    if (!span->count) {
        const alea_cell_entry_t* cell=q->sys->cells.data+id;
        alea_raycast_result_t scratch = {0};
        scratch.breakpoint_hit_limit = SIZE_MAX;
        alea_raycast_tree_primitives(q->sys, &level->local, cell->root_node_id,
            level->lo.t, level->hi.t, &scratch, NULL);
        int rc = scratch.breakpoint_failed ? -1 : 0;
        for (size_t i = 0; rc == 0 && i < scratch.hits.count; i++) {
            const alea_ray_hit_t* h = &scratch.hits.data[i];
            double t = h->t;
            if (level->corrections.count) t = corrected_root(level, h->primitive_id, t);
            if (!(t > level->lo.t && t < level->hi.t)) continue;
            interval_end_t end = {t, h->surface_id, h->primitive_id};
            rc = alea_vec_push(&level->cuts, end, interval_end_t);
        }
        q->out->surfaces_tested += scratch.surfaces_tested;
        alea_raycast_result_free(&scratch);
        return rc;
    }
    if (!level->specialized) return append_card_roots(level,q->sys->cells.data[id].surface_indices,
        span->count,level->lo.t,level->hi.t);
    for (uint32_t k=level->head;k!=UINT32_MAX;k=level->patches.data[k].next) {
        ray_patch_t ref=level->patches.data[k];
        spatial_patch_t p=idx->patches.data[ref.item];
        const uint32_t* cards = NULL;
        if (p.count) cards = p.first == UINT32_MAX ? q->sys->cells.data[id].surface_indices
                                                    : idx->patch_cards.data+p.first;
        if (p.count > p.active) level->has_inactive = true;
        if (level->recovering_boundaries) {
            if(p.count>p.active &&
               append_card_roots(level,cards+p.active,p.count-p.active,ref.lo,ref.hi))return -1;
        }else if(append_card_roots(level,cards,p.active,ref.lo,ref.hi))return -1;
    }
    return 0;
}

static int lattice_runs(interval_level_t* l, uint32_t index) {
    const alea_cell_entry_t* cell = &l->query->sys->cells.data[index];
    double a, b;
    int rc = alea_lattice_raycast_interval(&l->local, cell, l->lo.t, l->hi.t, &a, &b);
    if (rc != 0) return rc < 0 ? -1 : 0;
    if (!(a < b)) return 0;
    int limit = alea_lattice_raycast_step_limit(&l->local, cell, a, b);
    for (int step = 0; step < limit && a < b; step++) {
        double end = alea_lattice_next_boundary(&l->local, cell, a, b);
        if (!(end > a) || end > b) return -1;
        double x, y, z;
        alea_ray_point_at(&l->local, a + (end - a) * 0.5, &x, &y, &z);
        cell_run_t run = {.cell = index,
            .lo = a == l->lo.t ? l->lo : (interval_end_t){a, 0, UINT32_MAX},
            .hi = end == l->hi.t ? l->hi : (interval_end_t){end, 0, UINT32_MAX}};
        rc = alea_lattice_locate_point(l->query->sys, cell, x, y, z, &run.element);
        if (rc < 0 || (rc && alea_vec_push(&l->runs, run, cell_run_t))) return -1;
        a = end;
    }
    return a < b ? -1 : 0;
}

static inline bool atom_contains(const interval_index_t* idx,uint32_t code,double x,double y,double z) {
    if(code<=GUARD_ANY)return code==GUARD_ANY;
    double value;
    if(code&GUARD_REF) {
        value=alea_evaluate_point(idx->sys,(code&~GUARD_REF)>>1,x,y,z);
    }else if(code&GUARD_BOX) {
        bool inside=alea_node_bbox_contains_point(&idx->boxes.data[(code&~GUARD_BOX)>>1],x,y,z);
        return inside!=(bool)(code&1);
    }else if(code&QUERY_PLANE) {
        uint32_t id = (code & QUERY_OPERAND) >> 1;
        const alea_plane_data_t* p = idx->root_only ? &idx->sys->primitive_planes.data[id] : &idx->planes[id].data;
        value=p->a*x+p->b*y+p->c*z+p->d;
    }else if(code&QUERY_QUADRIC) {
        uint32_t id = (code & QUERY_OPERAND) >> 1;
        const double* k = idx->root_only ? idx->sys->primitive_quadrics.data[id].coeffs : idx->quadrics[id].data.coeffs;
        value=guard_value(k,x,y,z);
    }else {
        uint32_t id=(code-2)>>1;
        value=alea_primitive_eval_payload(idx->sys->primitives.data[id].type,alea_primitive_payload_const(idx->sys,id),x,y,z);
    }
    return (code&1)?value>=0:value<=0;
}
static bool guard_contains(const interval_index_t* idx,uint32_t code,double x,double y,double z) {
    if((code&(GUARD_BRANCH|QUERY_PROGRAM))==QUERY_PROGRAM) {
        const uint32_t* pc=idx->decisions.data+(code&~QUERY_PROGRAM);
        for(;;) {
            query_decision_t ins=unpack_decision(idx,*pc);
            unsigned jump=atom_contains(idx,ins.operand,x,y,z)?ins.yes:ins.no;
            if(jump==0 || jump==QUERY_ACCEPT)return jump!=0;
            pc-=jump;
        }
    }
    if(code&GUARD_BRANCH) {
        guard_pair_t p=idx->guards.pairs.data[code&~(GUARD_BRANCH|GUARD_UNION)];
        bool a=guard_contains(idx,p.a,x,y,z);
        if(a==!!(code&GUARD_UNION))return a;
        return guard_contains(idx,p.b,x,y,z);
    }
    return atom_contains(idx,code,x,y,z);
}

/* A midpoint well away from zero has the same sign in the 1-D ray
 * polynomial and the FP64 primitive. Resolve uncertain signs in 3-D. */
static inline bool ray_atom(const interval_level_t* l,uint32_t code,double t,
                            double x,double y,double z) {
    const interval_index_t* idx=l->query->index;
    if((code&~QUERY_OPERAND)==QUERY_QUADRIC) {
        ray_gate_t* g=l->query->worker->quadric_gates+((code&QUERY_OPERAND)>>1);
        if(g->stamp==l->query->worker->generation) {
            double v=(g->a*t+g->b)*t+g->c;
            if(isfinite(v) && fabs(v)>g->error)return (v>0)==(bool)(code&1);
        }
    }
    return atom_contains(idx,code,x,y,z);
}
static bool ray_guard_contains(const interval_level_t* l,uint32_t code,double t,
                               double x,double y,double z) {
    const interval_index_t* idx=l->query->index;
    if((code&(GUARD_BRANCH|QUERY_PROGRAM))==QUERY_PROGRAM) {
        const uint32_t* pc=idx->decisions.data+(code&~QUERY_PROGRAM);
        for(;;) {
            query_decision_t ins=unpack_decision(idx,*pc);
            unsigned jump=ray_atom(l,ins.operand,t,x,y,z)?ins.yes:ins.no;
            if(!jump || jump==QUERY_ACCEPT)return jump!=0;
            pc-=jump;
        }
    }
    if(code&GUARD_BRANCH) {
        guard_pair_t pair=idx->guards.pairs.data[code&~(GUARD_BRANCH|GUARD_UNION)];
        bool left=ray_guard_contains(l,pair.a,t,x,y,z);
        if(left==!!(code&GUARD_UNION))return left;
        return ray_guard_contains(l,pair.b,t,x,y,z);
    }
    return ray_atom(l,code,t,x,y,z);
}

static bool contains_sample(const interval_level_t* l,uint32_t id,double t,double x,double y,double z) {
    const interval_index_t* idx=l->query->index;
    if (!l->specialized) return alea_contains_point(idx->sys,idx->sys->cells.data[id].root_node_id,x,y,z);
    for (uint32_t k=l->head;k!=UINT32_MAX;k=l->patches.data[k].next) {
        ray_patch_t ref=l->patches.data[k];
        if (t>=ref.lo && t<=ref.hi)
            return ray_guard_contains(l,idx->patches.data[ref.item].program,t,x,y,z);
    }
    return false;
}

static int collect_runs(uint32_t index, double box_enter, double box_exit, void* ctx) {
    interval_level_t* l = ctx;
    l->query->out->blas_cell_candidates++;
    if (l->query->index->cell_cards[index].lattice) return lattice_runs(l,index);
    l->cuts.count = 0;
    l->has_inactive = false;
    l->recovering_boundaries = false;
    if (cell_roots(l,index)) return -1;
    if (l->cuts.count > 1)
        sort_ends(&l->cuts);
    interval_end_t start = l->lo;
    /* Between successive roots the full cell expression has constant sign.
     * Merge adjacent true intervals of this cell before selecting an owner. */
    size_t first = l->runs.count;
    for (size_t i = 0; i <= l->cuts.count; i++) {
        interval_end_t end = i == l->cuts.count ? l->hi : l->cuts.data[i];
        if (end.t <= start.t) continue;
        double t = start.t + (end.t - start.t) * 0.5;
        if (t >= box_enter && t <= box_exit) {
            double x, y, z;
            bool narrow = end.t-start.t <= 64*DBL_EPSILON*l->magnitude;
            double length = l->local_direction_length;
            /* Keep the established local-distance midpoint for the original
             * CSG test near a boundary. Rounding the point in the other
             * order can put it on an overlapping face. Only this containment
             * sample uses local distances. Every crossing keeps world t. */
            if (narrow && length != 1.0 && isfinite(end.t * length)) {
                double lo = start.t * length, hi = end.t * length;
                double mid = lo + 0.5 * (hi - lo);
                x = l->local.ox + mid * (l->local.dx / length);
                y = l->local.oy + mid * (l->local.dy / length);
                z = l->local.oz + mid * (l->local.dz / length);
            } else {
                alea_ray_point_at(&l->local, t, &x, &y, &z);
            }
            bool inside;
            if (narrow) {
                const alea_system_t* sys = l->query->sys;
                alea_node_id_t node = sys->cells.data[index].root_node_id;
                int sign = 0;
                if (t > start.t && t < end.t) {
                    sign = alea_ray_interval_midpoint_sign(sys, node, &l->local, start.primitive, end.primitive, t);
                } else if ((start.t > l->lo.t && end.t < l->hi.t) || l->shares_root_frame) {
                    /* Internal cuts both belong to this local ray. A clip
                     * inherited from a fill can use another placement. */
                    sign = alea_ray_adjacent_plane_cell_sign(sys, node, &l->local, start.primitive,
                                                            end.primitive, start.t, end.t);
                }
                inside = sign ? sign < 0 : alea_contains_point(sys, node, x, y, z);
            } else {
                inside = contains_sample(l, index, t, x, y, z);
            }
            if (inside) {
                if (l->runs.count > first && alea_vec_last(&l->runs).hi.t == start.t) {
                    alea_vec_last(&l->runs).hi = end;
                } else {
                    cell_run_t run = {.cell = index, .lo = start, .hi = end};
                    if (alea_vec_push(&l->runs, run, cell_run_t)) return -1;
                }
            }
        }
        start = end;
    }
    /* Inactive cards can name an exact boundary or introduce an adjacent
     * floating-point interval. Check both before publishing the cell runs. */
    if (l->has_inactive && l->runs.count > first) {
        l->recovering_boundaries = true;
        l->cuts.count=0;
        if(cell_roots(l,index))return -1;
        for(size_t k=0;k<l->cuts.count;k++) {
            interval_end_t end=l->cuts.data[k];
            for(size_t i=first;i<l->runs.count;i++) {
                cell_run_t* run=l->runs.data+i;
                double pad=64*DBL_EPSILON*l->magnitude;
                if ((end.t != run->lo.t && fabs(end.t-run->lo.t)<=pad) ||
                    (end.t != run->hi.t && fabs(end.t-run->hi.t)<=pad)) {
                    /* A redundant card splits off a representable sliver.
                     * Replay just this cell with its original root ordering. */
                    l->runs.count=first;
                    l->specialized=false;
                    int rc=collect_runs(index,l->lo.t,l->hi.t,l);
                    l->specialized=true;
                    return rc;
                }
                if(end.t==run->lo.t && end.surface<run->lo.surface)run->lo=end;
                if(end.t==run->hi.t && end.surface<run->hi.surface)run->hi=end;
            }
        }
    }
    return 0;
}

/* Traverse a closed pair of child boxes. Pad the split in coordinate
 * space for FP64 ray evaluation and round the crossing distance outward. */
static int visit_region(interval_level_t* l, uint32_t node,double lo,double hi) {
    const interval_index_t* idx = l->query->index;
    spatial_node_t n = idx->nodes.data[node];
    if (!n.right) {
        for(uint32_t offset=0;offset<n.count;offset+=BOUND_LANES){
            const bound_packet_t* pack=&idx->packets.data[n.packet+offset/BOUND_LANES];
            double aa[BOUND_LANES],bb[BOUND_LANES];
            unsigned valid=packet_ray(l,pack,lo,hi,aa,bb);
            for(unsigned lane=0;valid && lane<BOUND_LANES && offset+lane<n.count;lane++) {
                if(!(valid&(1u<<lane)))continue;
                uint32_t i=offset+lane;
                double a=aa[lane],z=bb[lane];if(a>z)continue;
                const spatial_patch_t* patch=idx->patches.data+n.first+i;
                if(patch->active>=3 && !clip_program(l,patch->program,&a,&z))continue;
                uint32_t id=patch->cell;
                visited_cell_t* visit = &l->query->worker->visited[id];
                if (visit->generation != l->query->worker->generation) {
                    visit->generation=l->query->worker->generation;
                    visit->slot=(uint32_t)l->possible.count;
                    possible_cell_t candidate={id,UINT32_MAX,a,z};
                    if (alea_vec_push(&l->possible,candidate,possible_cell_t)) return -1;
                } else {
                    possible_cell_t* c=&l->possible.data[visit->slot];
                    if (a<c->lo) c->lo=a;
                    if (z>c->hi) c->hi=z;
                }
                possible_cell_t* c=&l->possible.data[visit->slot];
                if(l->patches.count>=UINT32_MAX)return -1;
                ray_patch_t ref={n.first+i,c->head,a,z};
                c->head=(uint32_t)l->patches.count;
                if (alea_vec_push(&l->patches,ref,ray_patch_t)) return -1;
            }
        }
        return 0;
    }
    unsigned axis=n.axis;
    double origin=l->origin[axis],dir=l->direction[axis],pad=l->pad[axis]+32*DBL_EPSILON*fabs(n.split);
    if (dir==0) {
        if(origin<=n.split+pad && visit_region(l,n.first,lo,hi))return -1;
        if(origin>=n.split-pad && visit_region(l,n.right,lo,hi))return -1;
        return 0;
    }
    uint32_t first=dir>0?n.first:n.right,second=dir>0?n.right:n.first;
    double center=(n.split-origin)/dir,delta=pad/fabs(dir)+64*DBL_EPSILON*(1+fabs(center));
    double upper=center+delta,lower=center-delta;
    if (lo<=upper && visit_region(l,first,lo,upper<hi?upper:hi))return -1;
    if (hi>=lower && visit_region(l,second,lower>lo?lower:lo,hi))return -1;
    return 0;
}

static int visit_sector(interval_level_t* l,uint32_t node,const alea_bbox_t* box) {
    interval_box_t b=interval_box_from(*box);
    double a=l->lo.t,z=l->hi.t;
    for(unsigned k=0;k<3;k++){
        double o=l->origin[k],v=l->direction[k],pad=l->pad[k];
        if(!v) {
            if(o<b.lo[k]-pad || o>b.hi[k]+pad)return 0;
        }else{
            double x=(b.lo[k]-pad-o)/v,y=(b.hi[k]+pad-o)/v;
            if(v<0){double t=x;x=y;y=t;}
            if(x>a)a=x;
            if(y<z)z=y;
        }
    }
    return a<=z?visit_region(l,node,a,z):0;
}

static bool inside_domain(const alea_ray_t* ray, double t, const alea_bbox_t* box) {
    double x,y,z;
    alea_ray_point_at(ray,t,&x,&y,&z);
    double pad = 32 * DBL_EPSILON * (1 + fabs(ray->ox) + fabs(ray->oy) + fabs(ray->oz) + t);
    return x > box->min_x + pad && x < box->max_x - pad &&
           y > box->min_y + pad && y < box->max_y - pad &&
           z > box->min_z + pad && z < box->max_z - pad;
}

static int trace_level(interval_query_t*, int, const alea_matrix_t*, int,
                       interval_end_t, interval_end_t, int, const alea_matrix_t*, const alea_matrix_t*, bool);

static int trace_owner(interval_query_t* q, const cell_run_t* run,
                       const alea_matrix_t* frame, int parent,
                       interval_end_t lo, interval_end_t hi, int depth,
                       const alea_matrix_t* entry_frame, const alea_matrix_t* exit_frame, bool owner_exit) {
    int index = run ? (int)run->cell : parent;
    const query_cell_t* c=index<0?NULL:q->index->cell_cards+index;
    if (q->nav) {
        alea_hier_ray_path_t* path = q->nav->path;
        path->count = depth;
        if (run) {
            if (depth >= ALEA_HIER_RAY_PATH_MAX) return -1;
            const alea_cell_entry_t* cell = &q->sys->cells.data[index];
            path->entries[depth] = (alea_hier_ray_path_entry_t){
                .cell_index = run->cell, .cell_id = c->mc_cell_id,
                .material_id = c->material_id, .universe_id = cell->universe_id,
                .fill_universe = c->lattice ? run->element.fill_universe : c->fill_universe,
                .is_lattice = c->lattice, .lat_fill_universe = run->element.fill_universe,
                .lat_i = run->element.i, .lat_j = run->element.j, .lat_k = run->element.k,
                .lat_ox = run->element.ox, .lat_oy = run->element.oy, .lat_oz = run->element.oz,
                .depth = depth, .transform = *frame
            };
            ++path->count;
        }
    }
    if (run && c->container) {
        alea_matrix_t translation, child;
        const alea_matrix_t* transform;
        int universe = c->fill_universe;
        if (c->lattice) {
            universe = run->element.fill_universe;
            alea_matrix_identity(&translation);
            translation.m[3] = run->element.ox;
            translation.m[7] = run->element.oy;
            translation.m[11] = run->element.oz;
            translation.has_inverse = false;
            transform = &translation;
        } else if(!c->transformed && universe>0) {
            return trace_level(q,universe,frame,index,lo,hi,depth+1,entry_frame,exit_frame,owner_exit);
        } else {
            transform = alea_cell_fill_matrix(q->sys, run->cell);
        }
        if (!transform) return -1;
        alea_matrix_multiply(&child, frame, transform);
        if (!alea_matrix_invert(&child)) return -1;
        if (universe > 0)
            return trace_level(q, universe, &child, index, lo, hi, depth + 1, entry_frame, exit_frame, owner_exit);
    }
    if (q->nav) {
        q->nav->cell = index;
        q->nav->end = hi.t;
        q->nav->surface = hi.surface;
        q->nav->primitive = hi.primitive;
        q->nav->exit_frame = *exit_frame;
        q->nav->resolution_flags = c && c->container ? ALEA_RESOLVE_UNDEFINED_FILL : 0;
        return 1;
    }
    if (q->incident_exit.surface > 0 && q->incident_exit.t == lo.t) {
        lo = q->incident_exit;
        if (q->visibility && q->visibility->normal) entry_frame = &q->incident_frame;
    }
    if (q->visibility) {
        regional_visibility_t* v = q->visibility;
        if (!c || !c->material_id || hi.t <= RAY_EPSILON ||
            (v->accept && !v->accept(v->context, c->mc_cell_id, c->material_id, c->density))) return 0;
        /* Global recovery from unowned space and synthetic lattice faces
         * have renderer-specific rules for reporting a physical normal. */
        if (v->normal && (!entry_frame || lo.surface == 0)) return 1;
        alea_ray_first_visible_result_t* hit = v->hit;
        hit->t = lo.t;
        hit->cell_id = c->mc_cell_id;
        hit->material_id = c->material_id;
        hit->density = c->density;
        hit->resolution_flags = c->container ? ALEA_RESOLVE_UNDEFINED_FILL : 0;
        hit->surface_id = lo.surface;
        hit->primitive_id = lo.surface > 0 ? lo.primitive : UINT32_MAX;
        if (v->normal && hit->primitive_id < q->sys->primitives.count)
            alea_raycast_visible_normal(q->sys, q->world, lo.t, entry_frame, hit);
        hit->found = true;
        return 1; /* Stop the ordered sweep, including all enclosing fills. */
    }
    alea_ray_segment_t segment = {.path_index = UINT32_MAX, .enter_hit_index = -1, .t_enter = lo.t, .t_exit = hi.t,
        .cell_id = c ? c->mc_cell_id : -1, .material_id = c ? c->material_id : 0,
        .density = c ? c->density : 0, .enter_surface_id = lo.surface,
        .exit_surface_id = hi.surface,
        .resolution_flags = c && c->container ? ALEA_RESOLVE_UNDEFINED_FILL : 0};
    return alea_ray_append_selected_segment(q->out, &q->previous_cell, index, &segment);
}

/* Higher priority is the earlier cell in model order. Break ties between
 * runs from the same cell in their original order. */
static bool earlier_run(const cell_runs_t* runs,uint32_t a,uint32_t b) {
    uint32_t left=runs->data[a].cell,right=runs->data[b].cell;
    return left<right || (left==right && a<b);
}
static void activate_run(interval_level_t* l,uint32_t run) {
    size_t i=l->heap.count++;
    while(i) {
        size_t parent=(i-1)/2;
        if(!earlier_run(&l->runs,run,l->heap.data[parent]))break;
        l->heap.data[i]=l->heap.data[parent];i=parent;
    }
    l->heap.data[i]=run;
}
static void expire_run(interval_level_t* l) {
    uint32_t tail=l->heap.data[--l->heap.count];
    if(!l->heap.count)return;
    size_t i=0;
    while(2*i+1<l->heap.count) {
        size_t child=2*i+1;
        if(child+1<l->heap.count && earlier_run(&l->runs,l->heap.data[child+1],l->heap.data[child]))child++;
        if(!earlier_run(&l->runs,l->heap.data[child],tail))break;
        l->heap.data[i]=l->heap.data[child];i=child;
    }
    l->heap.data[i]=tail;
}

/* Sweep the union of the candidate runs before publishing any children.
 * Recompute nearly coincident quadric roots and sub-step plane cuts in
 * ordinary fills. The width only triggers investigation of the equations. */
static int collect_root_corrections(interval_level_t* l, bool fill) {
    const alea_system_t* sys = l->query->sys;
    interval_end_t covered = l->lo;
    for (size_t i = 0; i < l->cuts.count; ++i) {
        const cell_run_t* run = &l->runs.data[l->cuts.data[i].surface];
        double width = run->lo.t - covered.t;
        if (width > 0 && width <= 64 * DBL_EPSILON * l->magnitude &&
            covered.t > l->lo.t && covered.primitive != run->lo.primitive &&
            covered.primitive < sys->primitives.count && run->lo.primitive < sys->primitives.count) {
            const interval_end_t ends[] = {covered, run->lo};
            double cuts[] = {covered.t, run->lo.t};
            alea_primitive_type_t type = sys->primitives.data[covered.primitive].type;
            if (type == sys->primitives.data[run->lo.primitive].type) {
                if (type == ALEA_PRIMITIVE_QUADRIC) {
                    for (size_t j = 0; j < 2; ++j) {
                        const alea_quadric_data_t* q = alea_primitive_payload_const(sys,ends[j].primitive);
                        if (!q) return -1;
                        cuts[j] = refine_quadric_crossing(&l->local,q->coeffs,ends[j].t);
                    }
                } else if (fill && type == ALEA_PRIMITIVE_PLANE &&
                           run->lo.t <= nextafter(covered.t,INFINITY)) {
                    double enter, exit, physical_width;
                    if (alea_substep_plane_cuts(&l->local,
                            alea_primitive_payload_const(sys,covered.primitive),
                            alea_primitive_payload_const(sys,run->lo.primitive),
                            covered.t,&enter,&exit,&physical_width)) {
                        cuts[0] = enter;
                        cuts[1] = exit;
                    }
                }
            }
            for (size_t j = 0; j < 2; ++j) {
                if (cuts[j] == ends[j].t) continue;
                root_correction_t correction = {ends[j].t, cuts[j], ends[j].primitive};
                if (alea_vec_push(&l->corrections, correction, root_correction_t)) return -1;
            }
        }
        if (run->hi.t > covered.t) covered = run->hi;
    }
    return 0;
}

static int trace_level(interval_query_t* q, int universe, const alea_matrix_t* frame,
                       int parent, interval_end_t lo, interval_end_t hi, int depth,
                       const alea_matrix_t* entry_frame, const alea_matrix_t* exit_frame, bool clipped_by_exit) {
    if (depth >= MAX_FILL_RAYCAST_DEPTH || alea_ray_interrupted(q->out)) return -1;
    interval_level_t l = {.query = q, .lo = lo, .hi = hi};
    l.shares_root_frame = depth == 1 && parent >= 0 &&
        !q->index->cell_cards[parent].transformed && !q->index->cell_cards[parent].lattice;
    l.local_direction_length = 1.0;
    if (depth) {
        alea_ray_transform_inverse(frame, q->world, &l.local);
        double dx = l.local.dx, dy = l.local.dy, dz = l.local.dz;
        l.local_direction_length = sqrt(dx * dx + dy * dy + dz * dz);
        if (!isfinite(l.local_direction_length) || l.local_direction_length <= RAY_EPSILON) return -1;
    } else {
        l.local = *q->world;
    }
    regional_worker_t* worker = q->worker;
    interval_level_t* scratch = &worker->levels[depth];
    l.patches=scratch->patches;l.possible=scratch->possible;l.runs=scratch->runs;l.cuts=scratch->cuts;l.heap=scratch->heap;
    q->out->blas_universe_queries++;
    const alea_universe_t* u = alea_get_universe(q->sys, universe);
    const interval_index_t* idx = q->index;
    uint32_t root = u ? idx->universe_roots[u - q->sys->universes.data] : UINT32_MAX;
    int rc;
collect_level:
    l.patches.count=l.possible.count=l.runs.count=l.cuts.count=l.heap.count=0;
    if (++worker->generation == 0) {
        if (worker->visited)
            memset(worker->visited, 0, q->sys->cells.count * sizeof(*worker->visited));
        memset(worker->roots, 0, (worker->root_mask + 1) * sizeof(*worker->roots));
        memset(worker->plane_gates, 0, idx->plane_count * sizeof(ray_gate_t));
        memset(worker->quadric_gates, 0, idx->quadric_count * sizeof(ray_gate_t));
        worker->generation++;
    }
    /* Bound the local coordinates swept through the world-distance interval. */
    double ray_magnitude = 1 + fabs(l.local.ox) + fabs(l.local.oy) + fabs(l.local.oz) +
                           fmax(1.0, l.local_direction_length) * hi.t;
    l.magnitude=ray_magnitude;
    const alea_bbox_t* domain=idx->domains && u?idx->domains+(u-q->sys->universes.data):&idx->domain;
    bool core=inside_domain(&l.local,l.lo.t,domain) && inside_domain(&l.local,l.hi.t,domain);
    bool outer=!core && universe==0 && idx->exterior_ready &&
        inside_domain(&l.local,l.lo.t,&idx->envelope) && inside_domain(&l.local,l.hi.t,&idx->envelope);
    if (root!=UINT32_MAX && 256*DBL_EPSILON*ray_magnitude<idx->proof_margin && (core || outer)) {
        l.specialized=true;
        l.origin[0]=l.local.ox;l.origin[1]=l.local.oy;l.origin[2]=l.local.oz;
        l.direction[0]=l.local.dx;l.direction[1]=l.local.dy;l.direction[2]=l.local.dz;
        l.inverse[0]=l.local.inv_dx;l.inverse[1]=l.local.inv_dy;l.inverse[2]=l.local.inv_dz;
        for(unsigned k=0;k<3;k++)
            l.pad[k]=64*DBL_EPSILON*(1+fabs(l.origin[k])+fabs(l.direction[k])*l.hi.t);
        for(unsigned i=0;i<3;i++) {
            l.os[i]=l.origin[i];l.ds[i]=l.direction[i];l.iv[i]=l.inverse[i];l.err[i]=l.pad[i];
            unsigned j=(i+1)%3;
            for(unsigned k=0;k<2;k++) {
                unsigned axis=3+2*i+k;
                double sign=k?-1:1;
                l.os[axis]=l.origin[i]+sign*l.origin[j];
                l.ds[axis]=l.direction[i]+sign*l.direction[j];
                l.iv[axis]=l.ds[axis]?1/l.ds[axis]:0;
                l.err[axis]=4*(l.pad[i]+l.pad[j]);
            }
        }
        if(core)rc=visit_region(&l,root,l.lo.t,l.hi.t);
        else {
            rc=visit_sector(&l,root,domain);
            for(unsigned i=0;!rc && i<6;i++)
                rc=visit_sector(&l,idx->exterior[i].root,&idx->exterior[i].box);
        }
        for (size_t i=0; !rc && i<l.possible.count; i++) {
            possible_cell_t c=l.possible.data[i];
            l.head=c.head;
            rc=collect_runs(c.cell,c.lo,c.hi,&l);
        }
    } else {
        size_t before=q->out->blas_cell_candidates,roots=q->out->surfaces_tested;
        rc=alea_hier_spatial_visit_universe_ray(q->sys,universe,&l.local,l.lo.t,l.hi.t,collect_runs,&l);
        worker->fallback.visits++;
        worker->fallback.candidates+=q->out->blas_cell_candidates-before;
        worker->fallback.intersections+=q->out->surfaces_tested-roots;
    }
    if (rc < 0) goto done;
    rc = -1;
    l.cuts.count = 0;
    if(l.runs.count>INT_MAX || alea_vec_reserve(&l.heap,l.runs.count,uint32_t))goto done;
    for(size_t i=0;i<l.runs.count;i++) {
        cell_run_t* run=l.runs.data+i;
        if(run->lo.t>=run->hi.t)continue;
        interval_end_t event={run->lo.t,(int)i, 0};
        if(alea_vec_push(&l.cuts,event,interval_end_t))goto done;
    }
    sort_ends(&l.cuts);
    if (!l.corrections.count && (idx->quadric_count || parent >= 0) && l.cuts.count > 1) {
        if (collect_root_corrections(&l, parent >= 0)) goto done;
        /* Reclassify every candidate using the same corrected cuts. No
         * segments or renderer callbacks have been published at this level. */
        if (l.corrections.count) goto collect_level;
    }
    /* Only a new run or the current winner's exit can change the owner. */
    size_t cursor=0;
    const cell_run_t* owner=NULL;
    interval_end_t start=lo;
    const alea_matrix_t* start_frame = entry_frame;
    double t=lo.t;
    while(t<hi.t) {
        while(cursor<l.cuts.count && l.cuts.data[cursor].t<=t)
            activate_run(&l,(uint32_t)l.cuts.data[cursor++].surface);
        while(l.heap.count && l.runs.data[l.heap.data[0]].hi.t<=t)expire_run(&l);
        const cell_run_t* next=l.heap.count?l.runs.data+l.heap.data[0]:NULL;
        if(next!=owner) {
            interval_end_t boundary={t,lo.surface,lo.primitive};
            if(owner && owner->hi.t==t)boundary=owner->hi;
            else if(next && next->lo.t==t)boundary=next->lo;
            if (start.t < t) {
                rc = trace_owner(q,owner,frame,parent,start,boundary,depth,start_frame,frame,
                                 owner && owner->hi.t == t && boundary.surface > 0);
                if (rc) goto done;
            }
            /* A child clipped by its parent's entry inherits that boundary's
             * transform. Other crossings were computed in this universe. */
            start_frame = t == lo.t ? entry_frame :
                          (!owner && parent < 0 ? NULL : frame);
            start=boundary;owner=next;
        }
        t=next && next->hi.t<hi.t?next->hi.t:hi.t;
        if(cursor<l.cuts.count && l.cuts.data[cursor].t<t)t=l.cuts.data[cursor].t;
    }
    /* At a physical parent exit the incident terminal may cross its own
     * surface at exactly the clip. Its card and placement belong to both
     * sides of this crossing. An entry into a fill has no such preference. */
    bool terminal_exit = false;
    if (clipped_by_exit && owner && !idx->cell_cards[owner->cell].container) {
        uint32_t primitive;
        int surface = alea_ray_cell_surface_at(q->sys, &l.local, owner->cell,
                                                hi.t, q->out, &primitive);
        if (surface > 0) {
            hi = (interval_end_t){hi.t, surface, primitive};
            terminal_exit = true;
            exit_frame = frame;
        }
    }
    rc=trace_owner(q,owner,frame,parent,start,hi,depth,start_frame,exit_frame,clipped_by_exit);
    if (rc == 0 && terminal_exit) {
        q->incident_exit = hi;
        if (q->visibility && q->visibility->normal) q->incident_frame = *frame;
    }
done:
    alea_vec_free(&l.corrections);
    scratch->patches=l.patches;scratch->possible=l.possible;scratch->runs=l.runs;scratch->cuts=l.cuts;scratch->heap=l.heap;
    return rc;
}

#include "interval_render_tail.h"

/* Only trace the finite, indexed prefix. A miss cannot discard geometry
 * farther along an unbounded camera or shadow ray. */
int alea_regional_render_visible(regional_worker_t* worker, const alea_ray_t* ray, double limit,
                                int normal, alea_raycast_render_accept_t accept, void* context,
                                alea_raycast_result_t* scratch, alea_ray_first_visible_result_t* hit) {
    const interval_index_t* idx = worker->index;
    const alea_bbox_t* box = &idx->domain;
    const alea_bbox_t* eligible = idx->exterior_ready ? &idx->envelope : box;
    double span = fmax(box->max_x-box->min_x, fmax(box->max_y-box->min_y, box->max_z-box->min_z));
    double prefix = fmin(limit, span * 0.5);
    if (!idx->domains || !(prefix > 0)) goto native;
    if (!inside_domain(ray, 0, eligible) || !inside_domain(ray, prefix, eligible)) {
        /* Clipping moves the ray origin to the cut plane. Near the horizon
         * that can be far outside the index, even with the camera inside it.
         * The model-wide bounds can still certify the entire retained ray. */
        if (!render_tail_empty(idx, ray, 0, limit)) goto native;
        alea_raycast_result_clear(scratch);
        scratch->ray = *ray;
        *hit = (alea_ray_first_visible_result_t){.cell_id=-1, .surface_id=-1, .primitive_id=UINT32_MAX};
        return 0;
    }
    do {
        if (!worker->levels) worker->levels = calloc(MAX_FILL_RAYCAST_DEPTH, sizeof(*worker->levels));
        if (!worker->levels) goto native;
        alea_raycast_result_clear(scratch);
        scratch->ray = *ray;
        *hit = (alea_ray_first_visible_result_t){.cell_id=-1, .surface_id=-1, .primitive_id=UINT32_MAX};
        regional_visibility_t v = {.accept=accept, .context=context, .normal=normal, .hit=hit};
        interval_query_t q = {.visibility=&v, .sys=idx->sys, .world=ray, .out=scratch, .index=idx, .worker=worker};
        alea_matrix_t identity;
        alea_matrix_identity(&identity);
        int rc = trace_level(&q, 0, &identity, -1, (interval_end_t){0,-1,UINT32_MAX},
                             (interval_end_t){prefix,-1,UINT32_MAX}, 0, &identity, &identity, false);
        if (rc < 0 && alea_ray_interrupted(scratch)) return -1;
        if (hit->found) return 0;
        if (rc != 0) break;
        if (prefix == limit || render_tail_empty(idx, ray, prefix, limit)) return 0;

        /* Keep the short search for foreground hits. On a miss, extend it
         * once past the last possible central material, if still indexed. */
        double core_exit;
        if (!idx->render_bounds_ready ||
            !render_box_interval(ray, idx->render_core, prefix, limit, NULL, &core_exit)) break;
        double extended = fmin(limit, core_exit + 1e-7 * (1 + core_exit));
        if (!(extended > prefix)) break;
        prefix = extended;
    } while (inside_domain(ray, prefix, eligible));
native:
    return alea_raycast_render_first_visible_nocache(idx->sys, ray, limit, normal, accept, context, scratch, hit);
}

int alea_regional_trace_segments(regional_worker_t* worker, const alea_ray_t* ray, double limit,
                                 alea_raycast_result_t* out) {
    const interval_index_t* index = worker->index;
    if (!worker->levels) {
        worker->levels=calloc(MAX_FILL_RAYCAST_DEPTH,sizeof(*worker->levels));
        if (!worker->levels)return -1;
    }
    alea_raycast_result_clear(out);
    out->ray = *ray;
    interval_query_t q = {.sys = index->sys, .world = ray, .out = out, .previous_cell = -2, .index = index, .worker = worker};
    alea_matrix_t identity;
    alea_matrix_identity(&identity);
    int rc = trace_level(&q, 0, &identity, -1,
        (interval_end_t){0, -1, UINT32_MAX}, (interval_end_t){limit, -1, UINT32_MAX}, 0, &identity, &identity, false);
    if (rc) alea_raycast_result_clear(out);
    return rc;
}

int alea_regional_first_span(regional_worker_t* worker, const alea_ray_t* ray,
                            double start, double limit, alea_raycast_result_t* scratch,
                            regional_nav_span_t* span) {
    if (!(limit > start) || !isfinite(limit)) return -1;
    if (!worker->levels) {
        worker->levels = calloc(MAX_FILL_RAYCAST_DEPTH, sizeof(*worker->levels));
        if (!worker->levels) return -1;
    }
    alea_raycast_result_clear(scratch);
    scratch->ray = *ray;
    interval_query_t q = {.nav=span, .sys=worker->index->sys, .world=ray,
        .out=scratch, .index=worker->index, .worker=worker};
    alea_matrix_t identity;
    alea_matrix_identity(&identity);
    int rc = trace_level(&q, 0, &identity, -1, (interval_end_t){start,-1,UINT32_MAX},
                         (interval_end_t){limit,-1,UINT32_MAX}, 0, &identity, &identity, false);
    return rc == 1 && span->end > start && span->end <= limit ? 0 : -1;
}

size_t alea_regional_worker_candidate_tree_visits(const regional_worker_t* worker) {
    return worker->fallback.visits;
}

void alea_regional_index_retain(interval_index_t* idx) {
    atomic_fetch_add(&idx->references, 1);
}

static int trace_ray(interval_index_t* index, double ox, double oy, double oz,
                     double dx, double dy, double dz, double limit, alea_raycast_result_t* out) {
    /* Keep native rounding for ordinary directions and scale extreme inputs
     * before squaring them in the shared ray initializer. */
    double largest = fmax(fabs(dx), fmax(fabs(dy), fabs(dz)));
    if (largest > 0 && (largest <= RAY_EPSILON || largest >= sqrt(DBL_MAX / 3))) {
        dx /= largest;
        dy /= largest;
        dz /= largest;
    }
    alea_ray_t ray;
    if (!isfinite(limit) || !(limit > 0) || alea_ray_init(&ray, ox, oy, oz, dx, dy, dz))
        return -1;
    return alea_regional_trace_segments(&index->workers[0], &ray, limit, out);
}

/* Test only the slab directions that cut this particular octant. */
static bool point_item_covers(const interval_index_t* idx,const spatial_node_t* leaf,
                              point_item_t item,const double p[BOUND_AXES]) {
    uint32_t offset = item.bounds >> POINT_OFFSET_BITS;
    unsigned axes = item.bounds & POINT_AXES;
    const bound_packet_t* packet=idx->packets.data+leaf->packet+offset/BOUND_LANES;
    unsigned lane=offset%BOUND_LANES;
    for(unsigned k=0;k<3;k++)if(axes&(1u<<k)) {
        if(p[k]<packet->xyz[k].lo[lane] || p[k]>packet->xyz[k].hi[lane])return false;
    }
    if(axes>>3) {
        size_t pos=packet->first;
        for(unsigned k=0;k<6;k++)if(packet->mask&(1u<<k)) {
            const packet_axis_t* axis=idx->extra_axes.data+pos++;
            if((axes&(8u<<k)) && (p[3+k]<axis->lo[lane] || p[3+k]>axis->hi[lane]))return false;
        }
    }
    return true;
}

static int adaptive_child_point(const interval_index_t* idx, const alea_universe_t* u,
                                double x, double y, double z, size_t* tested) {
    struct universe_scans* scans = &idx->scans[u - idx->sys->universes.data];
    const size_t covered = atomic_load(&scans->covered);
    *tested = 0;
    if (covered) {
        const child_bound_block_t* block = scans->first;
        for (size_t pos = 0; pos < covered;) {
            size_t count = covered - pos;
            if (count > CHILD_BOUND_BLOCK_SIZE) count = CHILD_BOUND_BLOCK_SIZE;
            if (alea_node_bbox_contains_point(&block->envelope, x, y, z)) {
                alea_size_vec_t cells = {.data = u->cell_indices.data + pos, .count = count};
                size_t checks;
                int id = alea_scan_ordered_cells(idx->sys, &cells, block->cells, x, y, z, &checks);
                *tested += checks;
                if (id >= 0) return id;
            }
            pos += count;
            /* A concurrent builder may be writing the last published block's link. */
            if (pos < covered) block = block->next;
        }
        if (covered == u->cell_indices.count) return -1;
    }
    alea_size_vec_t tail = {.data = u->cell_indices.data + covered,
                            .count = u->cell_indices.count - covered};
    size_t scanned;
    double start = alea_monotonic_seconds();
    int id = alea_scan_ordered_cells(idx->sys, &tail, NULL, x, y, z, &scanned);
    uint64_t elapsed = (uint64_t)((alea_monotonic_seconds() - start) * 1e9);
    *tested += scanned;
    if (scanned < 8 && scanned < tail.count) return id;
    uint64_t previous = atomic_fetch_add_explicit(&scans->unfiltered_ns, elapsed, memory_order_relaxed);
    /* Require repeated work beyond the current bounds. Limit the construction
     * charged to any one lookup, even for a universe containing millions of cells. */
    if (previous && previous + elapsed >= 500000 && scanned)
        interval_extend_child_bounds(idx, u, covered);
    return id;
}

static spatial_node_t point_leaf(const interval_index_t* idx, size_t ui, uint32_t root,
                                 const double p[BOUND_AXES]) {
    const point_grid_t* grid=idx->point_grids+ui;
    if(grid->ready && root==idx->universe_roots[ui]) {
        unsigned tile=0;
        for(unsigned k=0;k<3;k++) {
            double q=(p[k]-grid->low[k])*grid->scale[k];
            unsigned j=q<=0?0:q>=15?15:(unsigned)q;
            tile|=j<<(4*k);
        }
        root=idx->point_nodes.data[grid->first+tile];
    }
    spatial_node_t node=idx->nodes.data[root];
    while(node.right)
        node=idx->nodes.data[p[node.axis]<=node.split?node.first:node.right];
    return node;
}

int alea_regional_ray_root_candidate(const interval_index_t* idx, double x, double y,
                                     double z, int crossed_surface, size_t* candidates) {
    if (!idx || !idx->domains) return -1;
    /* The reference's linear search ignores the surface hint, including
     * when a double coordinate exceeds the compact cell tree's range. */
    if (crossed_surface > 0 && (fabs(x) > FLT_MAX || fabs(y) > FLT_MAX || fabs(z) > FLT_MAX ||
                                !alea_hier_spatial_has_cell_tree(idx->sys, 0))) return -1;
    const alea_universe_t* u = alea_get_universe(idx->sys, 0);
    if (!u) return -1;
    size_t ui = (size_t)(u - idx->sys->universes.data);
    uint32_t root = idx->universe_roots[ui];
    if (root == UINT32_MAX || !alea_bbox_contains_point(&idx->domains[ui], x, y, z))
        return -1;

    double p[BOUND_AXES] = {x, y, z, x+y, x-y, y+z, y-z, z+x, z-x};
    spatial_node_t leaf = point_leaf(idx, ui, root, p);
    /* Point tiles may stop at a solid cell. Ray recovery can prefer a later
     * overlapping cell sharing the crossed surface. Use the full leaf in
     * deck order and verify each candidate with the original predicate. */
    for (uint32_t i = 0; i < leaf.count; i += BOUND_LANES) {
        const bound_packet_t* packet = &idx->packets.data[leaf.packet + i / BOUND_LANES];
        unsigned mask = packet_point(idx, packet, p);
        for (unsigned lane = 0; lane < BOUND_LANES && i + lane < leaf.count; ++lane) {
            if (!(mask & (1u << lane))) continue;
            uint32_t cell = idx->patches.data[leaf.first + i + lane].cell;
            ++*candidates;
            if (alea_ray_owner_candidate_contains(idx->sys, cell, x, y, z, crossed_surface))
                return (int)cell;
        }
    }
    return -1;
}

static int select_point_cell_impl(const interval_index_t* idx,alea_system_t* sys,
                             const alea_universe_t* u,double lx,double ly,double lz,
                             alea_lattice_location_t* element) {
    size_t ui = (size_t)(u - sys->universes.data);
    bool in_domain = alea_bbox_contains_point(idx->domains + ui,lx,ly,lz);
    /* The child has the same local working box as its parent point index. */
    if (in_domain && idx->child_mode == CHILD_POINTS_ADAPTIVE && idx->universe_roots[ui] == UINT32_MAX) {
        const interval_index_t* child = published_child_region(&idx->scans[ui]);
        if (child) idx = child;
    }
    uint32_t root=idx->universe_roots[ui];
    if (!in_domain) {
        root=UINT32_MAX;
        if(u->universe_id==0 && idx->exterior_ready)for(unsigned i=0;i<6;i++)
            if(alea_bbox_contains_point(&idx->exterior[i].box,lx,ly,lz)){
                root=idx->exterior[i].root;break;
            }
    }
    int selected=-1;
    if(root==UINT32_MAX) {
        if (idx->points_only) atomic_fetch_add(&sys->selected_fallbacks, 1);
        if (idx->root_only) {
            size_t tested;
            const alea_size_vec_t* cells = u->universe_id == 0
                ? &idx->exterior_point_cells : &u->cell_indices;
            struct universe_scans* scans = &idx->scans[u - sys->universes.data];
            /* In-domain children with a prepared region already have a root. */
            if (idx->child_mode == CHILD_POINTS_ADAPTIVE && u->universe_id != 0 && cells->count >= 8 && in_domain) {
                selected = adaptive_child_point(idx, u, lx, ly, lz, &tested);
            } else {
                const alea_node_bbox_t* bounds = in_domain ? scans->bounds : NULL;
                selected = alea_scan_ordered_cells(sys, cells, bounds, lx, ly, lz, &tested);
            }
            atomic_fetch_add_explicit(&scans->tests, tested, memory_order_relaxed);
            atomic_fetch_add_explicit(&scans->visits, 1, memory_order_relaxed);
        } else {
            selected=alea_hier_spatial_find_ordered_cell_in_universe(sys,u->universe_id,lx,ly,lz,-1);
        }
        if(selected < -1)return -2;
        if(selected>=0 && sys->cells.data[selected].lat_type && sys->cells.data[selected].lat_fill &&
           alea_lattice_locate_point(sys,&sys->cells.data[selected],lx,ly,lz,element)!=1)
            return -2;
    }else {
    double p[BOUND_AXES]={lx,ly,lz,lx+ly,lx-ly,ly+lz,ly-lz,lz+lx,lz-lx};
    spatial_node_t node = point_leaf(idx, ui, root, p);
    if (node.axis) {
        const point_tiles_t* tiles = idx->tiles.data + node.axis - 1;
        unsigned octant = (lx > tiles->center[0]) | ((ly > tiles->center[1]) << 1) | ((lz > tiles->center[2]) << 2);
        for (uint32_t i = tiles->first[octant]; i < tiles->first[octant] + tiles->count[octant]; i++) {
            point_item_t item = idx->point_items.data[i];
            if ((item.bounds & POINT_AXES) && !point_item_covers(idx, &node, item, p)) continue;
            if (!guard_contains(idx, item.program, lx, ly, lz)) continue;
            uint32_t cell = idx->patches.data[node.first + (item.bounds >> POINT_OFFSET_BITS)].cell;
            if ((item.bounds & POINT_LATTICE) &&
                alea_lattice_locate_point(sys, sys->cells.data + cell, lx, ly, lz, element) != 1) continue;
            return (int)cell;
        }
        return tiles->solid[octant] == UINT32_MAX ? -1 : (int)tiles->solid[octant];
    }
    for(uint32_t offset=0;offset<node.count && selected<0;offset+=BOUND_LANES) {
        const bound_packet_t* packet=&idx->packets.data[node.packet+offset/BOUND_LANES];
        unsigned possible=packet_point(idx,packet,p);
        for(unsigned lane=0;lane<BOUND_LANES && offset+lane<node.count;lane++) {
            if(!(possible&(1u<<lane)))continue;
            uint32_t item=node.first+offset+lane,id=idx->patches.data[item].cell;
            const alea_cell_entry_t* c=sys->cells.data+id;
            if(idx->patches.data[item].program==GUARD_ANY && c->lat_type && c->lat_fill) {
                if(alea_lattice_locate_point(sys,c,lx,ly,lz,element)==1){selected=(int)id;break;}
                continue;
            }
            if(guard_contains(idx,idx->patches.data[item].program,lx,ly,lz)){selected=(int)id;break;}
        }
    }
    }
    return selected;
}

static int select_point_cell(const interval_index_t* idx, alea_system_t* sys,
                             const alea_universe_t* u, double x, double y, double z,
                             alea_lattice_location_t* element) {
    if (!idx->time_child_points || u->universe_id == 0)
        return select_point_cell_impl(idx, sys, u, x, y, z, element);
    double start = alea_monotonic_seconds();
    int owner = select_point_cell_impl(idx, sys, u, x, y, z, element);
    uint64_t ns = (uint64_t)((alea_monotonic_seconds() - start) * 1e9);
    atomic_fetch_add_explicit(&idx->scans[u - sys->universes.data].local_ns, ns, memory_order_relaxed);
    return owner;
}

/* Only ordinary fills accumulate transforms. A lattice establishes a new
 * local origin, as in the public ordered point query. */
static int descend_point(const interval_index_t* idx,alea_system_t* sys,interval_point_t world,
                         const alea_universe_t* u,int owner,const alea_matrix_t* frame,
                         int depth,int* cell_id,int* material,int* index) {
    double x=world.x,y=world.y,z=world.z,wx=x,wy=y,wz=z,lx=x,ly=y,lz=z;
    alea_matrix_t accumulated;
    if(frame) {
        accumulated=*frame;
        alea_matrix_transform_point_inverse(frame,&lx,&ly,&lz);
    }else alea_matrix_identity(&accumulated);
    for(;u && depth<MAX_FILL_RAYCAST_DEPTH;depth++) {
        alea_lattice_location_t element={0};
        int selected=select_point_cell(idx,sys,u,lx,ly,lz,&element);
        if(selected < -1)goto hierarchy;
        if(selected<0)break;
        owner=selected;
        const alea_cell_entry_t* c=sys->cells.data+owner;
        if(c->lat_type && c->lat_fill) {
            if(element.fill_universe<=0)break;
            wx=lx-element.ox;wy=ly-element.oy;wz=lz-element.oz;
            lx=wx;ly=wy;lz=wz;
            alea_matrix_identity(&accumulated);
            u=alea_get_universe(sys,element.fill_universe);
            if(depth+1==MAX_FILL_RAYCAST_DEPTH)goto hierarchy;
            continue;
        }
        if(c->fill_universe<=0)break;
        if(depth+1==MAX_FILL_RAYCAST_DEPTH)goto hierarchy;
        u=alea_get_universe(sys,c->fill_universe);
        if(c->fill_transform<=0)continue;
        const alea_matrix_t* fill=alea_cell_fill_matrix(sys,owner);
        if(!fill)goto hierarchy;
        alea_matrix_t next;
        alea_matrix_multiply(&next,&accumulated,fill);
        if(!alea_matrix_invert(&next))goto hierarchy;
        accumulated=next;
        lx=wx;ly=wy;lz=wz;
        alea_matrix_transform_point_inverse(&accumulated,&lx,&ly,&lz);
    }
    if(owner<0)return -1;
    if(index)*index=owner;
    if(cell_id)*cell_id=sys->cells.data[owner].mc_cell_id;
    if(material)*material=sys->cells.data[owner].material_id;
    return 0;

hierarchy:
    if (idx->points_only) atomic_fetch_add(&sys->selected_fallbacks, 1);
    return alea_hierarchy_point_owner(sys, x, y, z, cell_id, material, index);
}

static int query_error(alea_error_t code, const char* detail) {
    alea_set_error_detail(code, "%s", detail);
    return -1;
}

static int validate_batch(const interval_index_t* idx, size_t count,
                          const void* inputs, const void* outputs) {
    if (!idx || (count && (!inputs || !outputs)))
        return query_error(ALEA_ERR_NULL_ARG, "Prepared queries require a handle and buffers for nonempty batches");
    if (count > UINT32_MAX)
        return query_error(ALEA_ERR_OVERFLOW, "Prepared query batch exceeds the index limit");
    if (idx->geometry_generation != alea_system_geometry_generation(idx->sys) ||
        idx->material_generation != alea_system_material_generation(idx->sys))
        return query_error(ALEA_ERR_INVALID_STATE, "Prepared queries require their original model generation");
    return 0;
}

int alea_regional_point_owner(const interval_index_t* idx, double x, double y, double z,
                              int* cell_id, int* material, int* index) {
    alea_system_t* sys = idx->sys;
    if(!idx->domains)return alea_hierarchy_point_owner(sys,x,y,z,cell_id,material,index);
    return descend_point(idx,sys,(interval_point_t){x,y,z},alea_get_universe(sys,0),
                         -1,NULL,0,cell_id,material,index);
}

/* Morton ordering improves locality between independent queries. Choosing
 * the order and scattering answers to their input positions are timed. */
static unsigned spread4(unsigned x) {
    return (x&1u) | ((x&2u)<<2) | ((x&4u)<<4) | ((x&8u)<<6);
}
static unsigned point_bin(double x,double low,double scale) {
    double v=(x-low)*scale;
    return v>0 ? (v<15 ? (unsigned)v : 15u) : 0;
}
static int order_queries(interval_index_t* idx,const interval_point_t* points,size_t count) {
    if(count>UINT32_MAX || count>SIZE_MAX/(2*sizeof(point_order_t)))
        return query_error(ALEA_ERR_OVERFLOW, "Prepared query order buffer exceeds the size limit");
    bool sort=count>=64 && idx->domains;
    if(sort && idx->order_capacity<count) {
        point_order_t* order=realloc(idx->order,2*count*sizeof(*order));
        if (!order) return query_error(ALEA_ERR_OUT_OF_MEMORY, "Allocating prepared query order buffer");
        idx->order=order;idx->order_capacity=count;
    }
    if(sort) {
        double low[3]={INFINITY,INFINITY,INFINITY},high[3]={-INFINITY,-INFINITY,-INFINITY};
        for(size_t i=0;i<count;i++) {
            double v[]={points[i].x,points[i].y,points[i].z};
            for(unsigned k=0;k<3;k++) {
                if(v[k]<low[k])low[k]=v[k];
                if(v[k]>high[k])high[k]=v[k];
            }
        }
        double scale[3];
        for(unsigned k=0;k<3;k++)scale[k]=high[k]>low[k]?16/(high[k]-low[k]):0;
        size_t bins[64];
        for(size_t i=0;i<count;i++) {
            interval_point_t p=points[i];
            unsigned key=spread4(point_bin(p.x,low[0],scale[0])) |
                        (spread4(point_bin(p.y,low[1],scale[1]))<<1) |
                        (spread4(point_bin(p.z,low[2],scale[2]))<<2);
            idx->order[i]=(point_order_t){key,(uint32_t)i};
        }
        for(unsigned pass=0;pass<2;pass++) {
            point_order_t* from=idx->order+(pass?count:0);
            point_order_t* to=idx->order+(pass?0:count);
            memset(bins,0,sizeof(bins));
            for(size_t i=0;i<count;i++)bins[(from[i].key>>(6*pass))&63]++;
            size_t offset=0;
            for(unsigned k=0;k<64;k++){size_t n=bins[k];bins[k]=offset;offset+=n;}
            for(size_t i=0;i<count;i++)to[bins[(from[i].key>>(6*pass))&63]++]=from[i];
        }
    }
    return sort;
}

int alea_prepared_locate_points(interval_index_t* idx,
                             const interval_point_t* points,size_t count,interval_owner_t* out) {
    if (validate_batch(idx, count, points, out)) return -1;
    alea_system_t* sys=idx->sys;
    for(size_t i=0;i<count;i++)
        if(!isfinite(points[i].x) || !isfinite(points[i].y) || !isfinite(points[i].z))
            return query_error(ALEA_ERR_INVALID_ARG, "Prepared point coordinates must be finite");
    int sort=order_queries(idx,points,count);
    if(sort<0)return -1;
    const alea_universe_t* world=alea_get_universe(sys,0);
    size_t filled=0;
    point_order_t* pending=sort?idx->order+count:NULL;
    uint32_t largest=0;
    for(size_t i=0;i<count;i++) {
        size_t k=sort?idx->order[i].index:i;
        interval_point_t p=points[k];
        out[k]=(interval_owner_t){-1,-1,-1};
        if(!sort || !world) {
            out[k].status=alea_regional_point_owner(idx,p.x,p.y,p.z,&out[k].cell_id,&out[k].material_id,NULL);
            continue;
        }
        alea_lattice_location_t element;
        int id=select_point_cell(idx,sys,world,p.x,p.y,p.z,&element);
        if(id==-1)continue;
        if(id<0 || idx->cell_cards[id].lattice) {
            out[k].status=alea_regional_point_owner(idx,p.x,p.y,p.z,&out[k].cell_id,&out[k].material_id,NULL);
            continue;
        }
        const query_cell_t* c=idx->cell_cards+id;
        out[k]=(interval_owner_t){0,c->mc_cell_id,c->material_id};
        if(c->fill_universe>0) {
            pending[filled++]=(point_order_t){(uint32_t)id,(uint32_t)k};
            if((uint32_t)id>largest)largest=(uint32_t)id;
        }
    }
    /* Share the child universe and its transform across points in the same
     * parent. Stable bucketing preserves spatial order within each group. */
    point_order_t* tmp=idx->order;
    for(unsigned shift=0;shift<32 && (largest>>shift);shift+=8) {
        size_t bins[256]={0},offset=0;
        for(size_t i=0;i<filled;i++)bins[(pending[i].key>>shift)&255]++;
        for(unsigned j=0;j<256;j++){size_t n=bins[j];bins[j]=offset;offset+=n;}
        for(size_t i=0;i<filled;i++)tmp[bins[(pending[i].key>>shift)&255]++]=pending[i];
        point_order_t* swap=tmp;tmp=pending;pending=swap;
    }
    size_t i=0;
    while(i<filled) {
        uint32_t parent=pending[i].key;
        const query_cell_t* c=idx->cell_cards+parent;
        const alea_universe_t* child=alea_get_universe(sys,c->fill_universe);
        const alea_matrix_t* frame=c->transformed?alea_cell_fill_matrix(sys,parent):NULL;
        do {
            size_t k=pending[i++].index;
            interval_point_t p=points[k];
            if(c->transformed && (!frame || !frame->has_inverse))
                out[k].status=alea_hierarchy_point_owner(sys,p.x,p.y,p.z,&out[k].cell_id,&out[k].material_id,NULL);
            else
                out[k].status=descend_point(idx,sys,p,child,(int)parent,frame,1,&out[k].cell_id,&out[k].material_id,NULL);
            if(out[k].status)out[k].cell_id=out[k].material_id=-1;
        }while(i<filled && pending[i].key==parent);
    }
    return 0;
}


int alea_prepared_trace_rays(interval_index_t* idx,
                           const interval_line_t* rays,size_t count,alea_raycast_result_t** out) {
    if (validate_batch(idx, count, rays, out)) return -1;
    for (size_t i=0;i<count;i++) {
        const interval_line_t* ray=rays+i;
        if (!out[i]) return query_error(ALEA_ERR_NULL_ARG, "Prepared ray requires a result buffer");
        if (!isfinite(ray->distance) || ray->distance <= 0)
            return query_error(ALEA_ERR_INVALID_ARG, "Prepared ray distance must be finite and positive");
        for(unsigned k=0;k<3;k++)
            if(!isfinite(ray->origin[k]) || !isfinite(ray->direction[k]))
                return query_error(ALEA_ERR_INVALID_ARG, "Prepared ray origin and direction must be finite");
        if(!ray->direction[0] && !ray->direction[1] && !ray->direction[2])
            return query_error(ALEA_ERR_INVALID_ARG, "Prepared ray direction must be nonzero");
    }
    alea_clear_error_detail();
    if (count && alea_system_prepare_query_caches(idx->sys, ALEA_CACHE_HIER_SPATIAL)) return -1;
    idx->workers[0].fallback=(interval_fallback_stats_t){0};
    int sort=0;
    if(count>=64 && idx->domains) {
        if (count > SIZE_MAX/sizeof(*idx->centers))
            return query_error(ALEA_ERR_OVERFLOW, "Prepared ray center buffer exceeds the size limit");
        if(idx->center_capacity<count) {
            interval_point_t* centers=realloc(idx->centers,count*sizeof(*centers));
            if (!centers) return query_error(ALEA_ERR_OUT_OF_MEMORY, "Allocating prepared ray center buffer");
            idx->centers=centers;idx->center_capacity=count;
        }
        for(size_t i=0;i<count;i++) {
            const double* origin=rays[i].origin;
            const double* dir=rays[i].direction;
            double n=sqrt(dir[0]*dir[0]+dir[1]*dir[1]+dir[2]*dir[2]);
            double mid=n>0 && isfinite(n)?0.5*rays[i].distance/n:0;
            idx->centers[i]=(interval_point_t){origin[0]+mid*dir[0],origin[1]+mid*dir[1],origin[2]+mid*dir[2]};
        }
        sort=order_queries(idx,idx->centers,count);
        if(sort<0)return -1;
    }
    for(size_t i=0;i<count;i++) {
        size_t k=sort?idx->order[i].index:i;
        const double* o=rays[k].origin;
        const double* v=rays[k].direction;
        if (trace_ray(idx,o[0],o[1],o[2],v[0],v[1],v[2],rays[k].distance,out[k])) {
            if (alea_interrupted()) return query_error(ALEA_ERR_INTERRUPTED, "Prepared ray traversal interrupted");
            if (alea_get_last_error() == ALEA_OK)
                return query_error(ALEA_ERR_INVALID_STATE, "Could not complete prepared ray traversal");
            return -1;
        }
    }
    return 0;
}

alea_prepared_queries_t* alea_prepared_queries_create(alea_system_t* sys, const alea_bbox_t* domain) {
    if (!sys || !domain) {
        query_error(ALEA_ERR_NULL_ARG, "Prepared queries require a model and a domain");
        return NULL;
    }
    interval_box_t box = interval_box_from(*domain);
    for (unsigned axis = 0; axis < 3; ++axis) {
        if (!(box.lo[axis] < box.hi[axis]) || !isfinite(box.hi[axis] - box.lo[axis])) {
            query_error(ALEA_ERR_INVALID_ARG, "Prepared query domain requires finite positive extents");
            return NULL;
        }
    }
    alea_clear_error_detail();
    if (alea_system_prepare_query_caches(sys, ALEA_CACHE_HIER_SPATIAL)) return NULL;
    interval_index_t* idx = interval_index_build(sys, domain, true);
    if (!idx) {
        if (alea_interrupted()) query_error(ALEA_ERR_INTERRUPTED, "Prepared query construction interrupted");
        else if (alea_get_last_error() == ALEA_OK)
            query_error(ALEA_ERR_INVALID_STATE, "Could not build the regional query index");
    }
    return idx;
}

int alea_set_ray_batch_accelerator(alea_system_t* sys, alea_prepared_queries_t* idx) {
    if (!sys) return query_error(ALEA_ERR_NULL_ARG, "Ray batch acceleration requires a model");
    if (idx && (idx->sys != sys || idx->points_only))
        return query_error(ALEA_ERR_INVALID_ARG, "Ray batch accelerator must belong to this model and support rays");
    if (idx && validate_batch(idx, 0, NULL, NULL)) return -1;
    if (idx != sys->ray_batch_index) {
        if (idx) atomic_fetch_add(&idx->references, 1);
        alea_prepared_queries_destroy(sys->ray_batch_index);
        sys->ray_batch_index = idx;
        if (idx) atomic_fetch_or(&sys->query_cache_state, ALEA_CACHE_BATCH_RAYS);
        else atomic_fetch_and(&sys->query_cache_state, ~ALEA_CACHE_BATCH_RAYS);
    }
    return 0;
}

bool alea_regional_nav_ready(const interval_index_t* idx, const alea_system_t* sys) {
    return idx && idx->sys == sys && !idx->points_only && validate_batch(idx, 0, NULL, NULL) == 0;
}

bool alea_regional_render_ready(const interval_index_t* idx, const alea_system_t* sys) {
    if (idx && idx->sys == sys) return validate_batch(idx, 0, NULL, NULL) == 0;
    query_error(ALEA_ERR_INVALID_ARG, "Regional rendering requires a prepared handle for this model");
    return false;
}
