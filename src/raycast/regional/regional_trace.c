// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
// SPDX-License-Identifier: MPL-2.0

/* Selected ownership along finite rays through nested CSG universes. */
#include "raycast/ray_traversal.h"
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

typedef struct {
    alea_raycast_render_accept_t accept;
    void* context;
    bool normal;
    alea_ray_first_visible_result_t* hit;
} regional_visibility_t;

typedef struct {
    regional_visibility_t* visibility;
    alea_system_t* sys;
    const alea_ray_t* world;
    alea_raycast_result_t* out;
    int previous_cell;
    const interval_index_t* index;
    regional_worker_t* worker;
} interval_query_t;

typedef struct { uint32_t cell,head; double lo,hi; } possible_cell_t;
ALEA_VEC_DEFINE(possible_cells,possible_cell_t);

typedef struct { uint32_t item,next; double lo,hi; } ray_patch_t;
ALEA_VEC_DEFINE(ray_patches,ray_patch_t);

typedef struct interval_level {
    double magnitude;
    double os[9],ds[9],iv[9],err[9];
    double origin[3],direction[3],inverse[3],pad[3];
    uint32_t head;
    ray_patches_t patches;
    possible_cells_t possible;
    bool specialized;
    bool recovering_boundaries, has_inactive;
    interval_query_t* query;
    alea_ray_t local;
    interval_end_t lo, hi;
    cell_runs_t runs;
    interval_ends_t cuts;
    cell_ids_t heap;
} interval_level_t;

#include "interval_packets.h"
#include "interval_ray_clip.h"

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
            if (!(h->t > level->lo.t && h->t < level->hi.t)) continue;
            interval_end_t end = {h->t, h->surface_id, h->primitive_id};
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
        const alea_plane_data_t* p=&idx->planes[(code&QUERY_OPERAND)>>1].data;
        value=p->a*x+p->b*y+p->c*z+p->d;
    }else if(code&QUERY_QUADRIC) {
        const double* k=idx->quadrics[(code&QUERY_OPERAND)>>1].data.coeffs;
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
            alea_ray_point_at(&l->local, t, &x, &y, &z);
            /* A midpoint can round onto either nearby surface. Use the
             * original FP64 expression to classify that narrow interval. */
            bool narrow = l->specialized && end.t-start.t <= 64*DBL_EPSILON*l->magnitude;
            bool inside = narrow ? alea_contains_point(l->query->sys,
                l->query->sys->cells.data[index].root_node_id,x,y,z) : contains_sample(l,index,t,x,y,z);
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
                       interval_end_t, interval_end_t, int, const alea_matrix_t*);

static int trace_owner(interval_query_t* q, const cell_run_t* run,
                       const alea_matrix_t* frame, int parent,
                       interval_end_t lo, interval_end_t hi, int depth,
                       const alea_matrix_t* entry_frame) {
    int index = run ? (int)run->cell : parent;
    const query_cell_t* c=index<0?NULL:q->index->cell_cards+index;
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
            return trace_level(q,universe,frame,index,lo,hi,depth+1,entry_frame);
        } else {
            transform = alea_hier_spatial_get_cell_fill_matrix(q->sys, run->cell);
        }
        if (!transform) return -1;
        alea_matrix_multiply(&child, frame, transform);
        if (!alea_matrix_invert(&child)) return -1;
        if (universe > 0)
            return trace_level(q, universe, &child, index, lo, hi, depth + 1, entry_frame);
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
static int trace_level(interval_query_t* q, int universe, const alea_matrix_t* frame,
                       int parent, interval_end_t lo, interval_end_t hi, int depth,
                       const alea_matrix_t* entry_frame) {
    if (depth >= MAX_FILL_RAYCAST_DEPTH || alea_interrupted()) return -1;
    interval_level_t l = {.query = q, .lo = lo, .hi = hi};
    double scale = 1.0;
    if (depth) {
        /* A rounded MCNP rotation may change the direction's norm. Kernels
         * take unit local directions, so their roots use local distances. */
        const alea_ray_t* w = q->world;
        double dx = frame->inv[0] * w->dx + frame->inv[1] * w->dy + frame->inv[2] * w->dz;
        double dy = frame->inv[4] * w->dx + frame->inv[5] * w->dy + frame->inv[6] * w->dz;
        double dz = frame->inv[8] * w->dx + frame->inv[9] * w->dy + frame->inv[10] * w->dz;
        scale = sqrt(dx * dx + dy * dy + dz * dz);
        double x = w->ox, y = w->oy, z = w->oz;
        alea_matrix_transform_point_inverse(frame, &x, &y, &z);
        if (!isfinite(scale) || alea_ray_init(&l.local, x, y, z, dx, dy, dz)) return -1;
        l.lo.t *= scale;
        l.hi.t *= scale;
        if (!isfinite(l.hi.t) || !(l.lo.t < l.hi.t)) return -1;
    } else {
        l.local = *q->world;
    }
    regional_worker_t* worker = q->worker;
    interval_level_t* scratch = &worker->levels[depth];
    l.patches=scratch->patches;l.possible=scratch->possible;l.runs=scratch->runs;l.cuts=scratch->cuts;l.heap=scratch->heap;
    l.patches.count=l.possible.count=l.runs.count=l.cuts.count=l.heap.count=0;
    q->out->blas_universe_queries++;
    const alea_universe_t* u = alea_get_universe(q->sys, universe);
    const interval_index_t* idx = q->index;
    uint32_t root = u ? idx->universe_roots[u - q->sys->universes.data] : UINT32_MAX;
    int rc;
    if (++worker->generation == 0) {
        if (worker->visited)
            memset(worker->visited, 0, q->sys->cells.count * sizeof(*worker->visited));
        memset(worker->roots, 0, (worker->root_mask + 1) * sizeof(*worker->roots));
        memset(worker->plane_gates, 0, idx->plane_count * sizeof(ray_gate_t));
        memset(worker->quadric_gates, 0, idx->quadric_count * sizeof(ray_gate_t));
        worker->generation++;
    }
    double ray_magnitude=1+fabs(l.local.ox)+fabs(l.local.oy)+fabs(l.local.oz)+l.hi.t;
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
        run->lo.t=run->lo.t==l.lo.t?lo.t:run->lo.t/scale;
        run->hi.t=run->hi.t==l.hi.t?hi.t:run->hi.t/scale;
        if(run->lo.t>=run->hi.t)continue;
        interval_end_t event={run->lo.t,(int)i, 0};
        if(alea_vec_push(&l.cuts,event,interval_end_t))goto done;
    }
    sort_ends(&l.cuts);
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
                rc = trace_owner(q,owner,frame,parent,start,boundary,depth,start_frame);
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
    rc=trace_owner(q,owner,frame,parent,start,hi,depth,start_frame);
done:
    scratch->patches=l.patches;scratch->possible=l.possible;scratch->runs=l.runs;scratch->cuts=l.cuts;scratch->heap=l.heap;
    return rc;
}

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
    if (idx->domains && prefix > 0 && inside_domain(ray, 0, eligible) &&
        inside_domain(ray, prefix, eligible)) {
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
                             (interval_end_t){prefix,-1,UINT32_MAX}, 0, &identity);
        if (rc < 0 && alea_interrupted()) return -1;
        if (hit->found || (rc == 0 && prefix == limit)) return 0;
    }
native:
    return alea_raycast_render_first_visible_nocache(idx->sys, ray, limit, normal, accept, context, scratch, hit);
}

static int trace_ray(interval_index_t* index, alea_system_t* sys, double ox, double oy, double oz,
                               double dx, double dy, double dz, double limit,
                               alea_raycast_result_t* out) {
    regional_worker_t* worker = &index->workers[0];
    if (!worker->levels) {
        worker->levels=calloc(MAX_FILL_RAYCAST_DEPTH,sizeof(*worker->levels));
        if (!worker->levels)return -1;
    }
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
    if (alea_system_prepare_query_caches(sys, ALEA_CACHE_HIER_SPATIAL | ALEA_CACHE_CELL_SURFACES))
        return -1;
    alea_raycast_result_clear(out);
    out->ray = ray;
    interval_query_t q = {.sys = sys, .world = &ray, .out = out, .previous_cell = -2, .index = index, .worker = worker};
    alea_matrix_t identity;
    alea_matrix_identity(&identity);
    int rc = trace_level(&q, 0, &identity, -1,
        (interval_end_t){0, -1, UINT32_MAX}, (interval_end_t){limit, -1, UINT32_MAX}, 0, &identity);
    if (rc) alea_raycast_result_clear(out);
    return rc;
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

static int select_point_cell(const interval_index_t* idx,alea_system_t* sys,
                             const alea_universe_t* u,double lx,double ly,double lz,
                             alea_lattice_location_t* element) {
    uint32_t root=idx->universe_roots[u-sys->universes.data];
    if(!alea_bbox_contains_point(idx->domains+(u-sys->universes.data),lx,ly,lz)) {
        root=UINT32_MAX;
        if(u->universe_id==0 && idx->exterior_ready)for(unsigned i=0;i<6;i++)
            if(alea_bbox_contains_point(&idx->exterior[i].box,lx,ly,lz)){
                root=idx->exterior[i].root;break;
            }
    }
    int selected=-1;
    if(root==UINT32_MAX) {
        selected=alea_hier_spatial_find_ordered_cell_in_universe(sys,u->universe_id,lx,ly,lz,-1);
        if(selected < -1)return -2;
        if(selected>=0 && sys->cells.data[selected].lat_type && sys->cells.data[selected].lat_fill &&
           alea_lattice_locate_point(sys,&sys->cells.data[selected],lx,ly,lz,element)!=1)
            return -2;
    }else {
    double p[BOUND_AXES]={lx,ly,lz,lx+ly,lx-ly,ly+lz,ly-lz,lz+lx,lz-lx};
    const point_grid_t* grid=idx->point_grids+(u-sys->universes.data);
    if(grid->ready && root==idx->universe_roots[u-sys->universes.data]) {
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
            if(idx->patches.data[item].program==GUARD_ANY && idx->cell_cards[id].lattice) {
                if(alea_lattice_locate_point(sys,c,lx,ly,lz,element)==1){selected=(int)id;break;}
                continue;
            }
            if(guard_contains(idx,idx->patches.data[item].program,lx,ly,lz)){selected=(int)id;break;}
        }
    }
    }
    return selected;
}

/* Only ordinary fills accumulate transforms. A lattice establishes a new
 * local origin, as in the public ordered point query. */
static int descend_point(interval_index_t* idx,alea_system_t* sys,interval_point_t world,
                         const alea_universe_t* u,int owner,const alea_matrix_t* frame,
                         int depth,int* cell_id,int* material) {
    double x=world.x,y=world.y,z=world.z,wx=x,wy=y,wz=z,lx=x,ly=y,lz=z;
    alea_matrix_t accumulated;
    if(frame) {
        accumulated=*frame;
        alea_matrix_transform_point_inverse(frame,&lx,&ly,&lz);
    }else alea_matrix_identity(&accumulated);
    for(;u && depth<MAX_FILL_RAYCAST_DEPTH;depth++) {
        alea_lattice_location_t element={0};
        int selected=select_point_cell(idx,sys,u,lx,ly,lz,&element);
        if(selected < -1)return alea_find_cell_at(sys,x,y,z,cell_id,material);
        if(selected<0)break;
        owner=selected;
        const query_cell_t* c=idx->cell_cards+owner;
        if(c->lattice) {
            if(element.fill_universe<=0)break;
            wx=lx-element.ox;wy=ly-element.oy;wz=lz-element.oz;
            lx=wx;ly=wy;lz=wz;
            alea_matrix_identity(&accumulated);
            u=alea_get_universe(sys,element.fill_universe);
            if(depth+1==MAX_FILL_RAYCAST_DEPTH)return alea_find_cell_at(sys,x,y,z,cell_id,material);
            continue;
        }
        if(c->fill_universe<=0)break;
        if(depth+1==MAX_FILL_RAYCAST_DEPTH)return alea_find_cell_at(sys,x,y,z,cell_id,material);
        u=alea_get_universe(sys,c->fill_universe);
        if(!c->transformed)continue;
        const alea_matrix_t* fill=alea_hier_spatial_get_cell_fill_matrix(sys,owner);
        if(!fill)return alea_find_cell_at(sys,x,y,z,cell_id,material);
        alea_matrix_t next;
        alea_matrix_multiply(&next,&accumulated,fill);
        if(!alea_matrix_invert(&next))return alea_find_cell_at(sys,x,y,z,cell_id,material);
        accumulated=next;
        lx=wx;ly=wy;lz=wz;
        alea_matrix_transform_point_inverse(&accumulated,&lx,&ly,&lz);
    }
    if(owner<0)return -1;
    if(cell_id)*cell_id=idx->cell_cards[owner].mc_cell_id;
    if(material)*material=idx->cell_cards[owner].material_id;
    return 0;
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

static int locate_point(interval_index_t* idx,alea_system_t* sys,
                             double x,double y,double z,int* cell_id,int* material) {
    if(!idx->domains)return alea_find_cell_at(sys,x,y,z,cell_id,material);
    return descend_point(idx,sys,(interval_point_t){x,y,z},alea_get_universe(sys,0),
                         -1,NULL,0,cell_id,material);
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
            out[k].status=locate_point(idx,sys,p.x,p.y,p.z,&out[k].cell_id,&out[k].material_id);
            continue;
        }
        alea_lattice_location_t element;
        int id=select_point_cell(idx,sys,world,p.x,p.y,p.z,&element);
        if(id==-1)continue;
        if(id<0 || idx->cell_cards[id].lattice) {
            out[k].status=locate_point(idx,sys,p.x,p.y,p.z,&out[k].cell_id,&out[k].material_id);
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
        const alea_matrix_t* frame=c->transformed?alea_hier_spatial_get_cell_fill_matrix(sys,parent):NULL;
        do {
            size_t k=pending[i++].index;
            interval_point_t p=points[k];
            if(c->transformed && (!frame || !frame->has_inverse))
                out[k].status=alea_find_cell_at(sys,p.x,p.y,p.z,&out[k].cell_id,&out[k].material_id);
            else
                out[k].status=descend_point(idx,sys,p,child,(int)parent,frame,1,&out[k].cell_id,&out[k].material_id);
            if(out[k].status)out[k].cell_id=out[k].material_id=-1;
        }while(i<filled && pending[i].key==parent);
    }
    return 0;
}


int alea_prepared_trace_rays(interval_index_t* idx,
                           const interval_line_t* rays,size_t count,alea_raycast_result_t** out) {
    if (validate_batch(idx, count, rays, out)) return -1;
    alea_system_t* sys=idx->sys;
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
        if (trace_ray(idx,sys,o[0],o[1],o[2],v[0],v[1],v[2],rays[k].distance,out[k])) {
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

bool alea_regional_render_ready(const interval_index_t* idx, const alea_system_t* sys) {
    if (idx && idx->sys == sys) return validate_batch(idx, 0, NULL, NULL) == 0;
    query_error(ALEA_ERR_INVALID_ARG, "Regional rendering requires a prepared handle for this model");
    return false;
}
