// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
// SPDX-License-Identifier: MPL-2.0
#include "compute_internal.h"
#include "core/alea_system.h"
#include "core/alea_transform.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

_Static_assert(sizeof(compute_render_bvh)==60,"Render BVH layout");
_Static_assert(sizeof(compute_render_params)==144,"Render parameters layout");
_Static_assert(sizeof(compute_render_pixel)==40,"Render pixel layout");
_Static_assert(sizeof(compute_lattice)==64,"Lattice layout");

typedef struct { float bounds[6], center[3]; uint32_t cell; } bvh_item;
/* Low bit marks procedural descendants. Upper bits estimate the complete
 * child query cost, so lazy clipping is selected only for sparse cheap contents. */
static uint32_t mark_cost(compute_scene* s,uint32_t id,uint32_t* memo) {
    if(id==COMPUTE_NONE)return 0;
    if(memo[id]!=UINT32_MAX)return memo[id];
    compute_render_bvh* b=&s->render.bvh[id];uint32_t left,right,flag;uint64_t cost;
    if(b->cell==COMPUTE_NONE) {
        left=mark_cost(s,id+1,memo);right=mark_cost(s,s->render.bvh[id+1].escape,memo);
        flag=(left|right)&COMPUTE_CELL_LATTICE_SUBTREE;cost=(uint64_t)(left>>1)+(right>>1);
    } else {
        compute_cell* c=&s->cells[b->cell];
        left=c->lattice!=COMPUTE_NONE?COMPUTE_CELL_LATTICE_SUBTREE:mark_cost(s,c->child,memo);
        c->flags=left;flag=left&COMPUTE_CELL_LATTICE_SUBTREE;cost=(uint64_t)c->count+(left>>1);
    }
    if(cost>UINT32_MAX/4)cost=UINT32_MAX/4;
    return memo[id]=((uint32_t)cost<<1)|flag;
}
static int compare_x(const void* a, const void* b) {
    float x = ((const bvh_item*)a)->center[0], y = ((const bvh_item*)b)->center[0];
    return (x > y) - (x < y);
}
static int compare_y(const void* a, const void* b) {
    float x = ((const bvh_item*)a)->center[1], y = ((const bvh_item*)b)->center[1];
    return (x > y) - (x < y);
}
static int compare_z(const void* a, const void* b) {
    float x = ((const bvh_item*)a)->center[2], y = ((const bvh_item*)b)->center[2];
    return (x > y) - (x < y);
}

static void build_bvh(bvh_item* items, size_t count, compute_render_scene* out) {
    size_t index = out->bvh_count++;
    compute_render_bvh* node = &out->bvh[index];
    memcpy(node->bounds, items[0].bounds, sizeof(node->bounds));
    double lo[3], hi[3];
    for (int a = 0; a < 3; ++a) lo[a] = hi[a] = items[0].center[a];
    for (size_t i = 1; i < count; ++i) for (int a = 0; a < 3; ++a) {
        node->bounds[2*a] = fminf(node->bounds[2*a], items[i].bounds[2*a]);
        node->bounds[2*a+1] = fmaxf(node->bounds[2*a+1], items[i].bounds[2*a+1]);
        lo[a] = fmin(lo[a], items[i].center[a]);
        hi[a] = fmax(hi[a], items[i].center[a]);
    }
    node->cell = count == 1 ? items[0].cell : UINT32_MAX;
    if (count > 1) {
        int axis = 0;
        if (hi[1]-lo[1] > hi[axis]-lo[axis]) axis = 1;
        if (hi[2]-lo[2] > hi[axis]-lo[axis]) axis = 2;
        qsort(items, count, sizeof(*items), axis == 0 ? compare_x : axis == 1 ? compare_y : compare_z);
        size_t half = count / 2;
        build_bvh(items, half, out);
        build_bvh(items + half, count - half, out);
    }
    node->escape = (uint32_t)out->bvh_count;
    memcpy(node->fallback,node->bounds,sizeof(node->bounds));node->guard=INFINITY;
}


typedef struct { int universe; uint32_t cell; } cell_order;
typedef struct { int id; size_t first,count; uint32_t root; unsigned state,height; } universe;
static int order_cells(const void* a,const void* b) {
    const cell_order* x=a;const cell_order* y=b;
    if(x->universe!=y->universe)return (x->universe>y->universe)-(x->universe<y->universe);
    return (x->cell>y->cell)-(x->cell<y->cell);
}
static int find_universe(const universe* u,size_t count,int id) {
    size_t a=0,b=count;
    while(a<b){size_t m=a+(b-a)/2;if(u[m].id<id)a=m+1;else b=m;}
    return a<count && u[a].id==id?(int)a:-1;
}
static int check_universe(const alea_system_t* sys,universe* u,size_t count,
                           const cell_order* order,int at,unsigned depth) {
    if(depth>=COMPUTE_HIER_STACK || depth+u[at].height>COMPUTE_HIER_STACK) {
        alea_set_error_detail(ALEA_ERR_UNSUPPORTED,"Metal universe depth exceeds %d",COMPUTE_HIER_STACK);return -1;
    }
    if(u[at].state==2)return 0;
    if(u[at].state==1){alea_set_error_detail(ALEA_ERR_INVALID_STATE,"Cycle in fill hierarchy at universe %d",u[at].id);return -1;}
    u[at].state=1;u[at].height=1;
    for(size_t k=u[at].first;k<u[at].first+u[at].count;++k) {
        const alea_cell_entry_t* c=&sys->cells.data[order[k].cell];
        size_t children=c->lat_type?c->lat_fill_count+(c->lat_outer_universe>0):c->fill_universe>0?1:0;
        for(size_t j=0;j<children;++j) {
            int id=c->lat_type?(j<c->lat_fill_count?c->lat_fill[j]:c->lat_outer_universe):c->fill_universe;
            if(id<0 || (c->lat_type && id==0))continue;
            int child=find_universe(u,count,id);
            if(child<0){alea_set_error_detail(ALEA_ERR_INVALID_STATE,"Missing fill universe %d in cell %d",id,c->mc_cell_id);return -1;}
            if(check_universe(sys,u,count,order,child,depth+1))return -1;
            if(u[child].height>=u[at].height)u[at].height=u[child].height+1;
        }
    }
    u[at].state=2;return 0;
}

void alea_compute_render_free(compute_render_scene* s) {
    free(s->bvh);free(s->transforms);free(s->lattices);free(s->fills);
    memset(s,0,sizeof(*s));
}

static void copy_bounds(float out[6],const double in[6]) {
    for(int a=0;a<3;++a) {
        if(isnan(in[2*a]) || isnan(in[2*a+1]) || in[2*a]>in[2*a+1]) {
            out[2*a]=-INFINITY;out[2*a+1]=INFINITY;
        } else {
            out[2*a]=nextafterf((float)in[2*a],-INFINITY);
            out[2*a+1]=nextafterf((float)in[2*a+1],INFINITY);
        }
    }
}

int alea_compute_render_pack(const alea_system_t* sys,compute_scene* scene) {
    compute_render_scene* out=&scene->render; memset(out,0,sizeof(*out));out->root=COMPUTE_NONE;
    size_t n=scene->cell_count,nu=0,nt=0,nl=0,nf=0;
    if(n>UINT32_MAX/2 || n>SIZE_MAX/(2*sizeof(compute_render_bvh))) {
        alea_set_error_detail(ALEA_ERR_OVERFLOW,"Too many render cells");return -1;
    }
    float (*bounds)[6]=calloc(n?n:1,sizeof(*bounds));
    out->bvh=calloc(n?2*n:1,sizeof(*out->bvh));
    bvh_item* items=calloc(n?n:1,sizeof(*items));
    cell_order* order=calloc(n?n:1,sizeof(*order));
    universe* univ=calloc(n?n:1,sizeof(*univ));
    if(!bounds || !out->bvh || !items || !order || !univ)goto memory;
    for(size_t i=0;i<n;++i) {
        const alea_cell_entry_t* c=&sys->cells.data[i];
        order[i]=(cell_order){c->universe_id,(uint32_t)i};
        if(c->fill_universe>0)++nt;
        if(c->lat_type) {
            if(c->lat_type!=1 || !c->lat_fill || !c->lat_fill_count || c->fill_transform>0) {
                alea_set_error_detail(ALEA_ERR_UNSUPPORTED,"Metal requires a defined rectangular lattice without a lattice-local fill transform (cell %d)",c->mc_cell_id);goto fail;
            }
            for(int a=0;a<3;++a)if(c->lat_fill_dims[2*a]>c->lat_fill_dims[2*a+1] ||
                !isfinite(c->lat_pitch[a]) || c->lat_pitch[a]<0 || !isfinite(c->lat_lower_left[a])) {
                alea_set_error_detail(ALEA_ERR_INVALID_STATE,"Invalid lattice bounds (cell %d)",c->mc_cell_id);goto fail;
            }
            if(c->lat_pitch[0]<=0 || c->lat_pitch[1]<=0 ||
                (c->lat_fill_dims[4]!=c->lat_fill_dims[5] && c->lat_pitch[2]<=0)) {
                alea_set_error_detail(ALEA_ERR_INVALID_STATE,"Invalid lattice pitch (cell %d)",c->mc_cell_id);goto fail;
            }
            uint64_t elements=1;
            for(int a=0;a<3;++a) {
                int64_t length=(int64_t)c->lat_fill_dims[2*a+1]-c->lat_fill_dims[2*a]+1;
                if(length>INT32_MAX || elements>UINT32_MAX/(uint64_t)length) {
                    alea_set_error_detail(ALEA_ERR_OVERFLOW,"Lattice dimensions exceed device indices");goto fail;
                }
                elements*=length;
            }
            if((!c->lat_fill_repeating && elements!=c->lat_fill_count) ||
               (c->lat_fill_repeating && c->lat_fill_count!=1)) {
                alea_set_error_detail(ALEA_ERR_INVALID_STATE,"Lattice fill dimensions and table disagree (cell %d)",c->mc_cell_id);goto fail;
            }
            if(c->lat_fill_count>UINT32_MAX-nf){alea_set_error_detail(ALEA_ERR_OVERFLOW,"Lattice fill table too large");goto fail;}
            nf+=c->lat_fill_count;++nl;
        }
        const alea_node_t* root=&sys->nodes.data[c->root_node_id];
        double b[6]={root->bbox.min_x,root->bbox.max_x,root->bbox.min_y,root->bbox.max_y,root->bbox.min_z,root->bbox.max_z};
        if(c->lat_type)for(int a=0;a<3;++a) {
            if(c->lat_fill_repeating || c->lat_outer_universe>0 || c->lat_pitch[a]==0) {b[2*a]=-INFINITY;b[2*a+1]=INFINITY;}
            else {b[2*a]=c->lat_lower_left[a];b[2*a+1]=b[2*a]+((double)c->lat_fill_dims[2*a+1]-c->lat_fill_dims[2*a]+1)*c->lat_pitch[a];}
        }
        copy_bounds(bounds[i],b);
    }
    qsort(order,n,sizeof(*order),order_cells);
    for(size_t i=0;i<n;) {
        size_t first=i;int id=order[i].universe;
        while(i<n && order[i].universe==id)++i;
        univ[nu++]=(universe){id,first,i-first,COMPUTE_NONE,0,0};
    }
    int root=find_universe(univ,nu,0);
    if(root<0 && n){alea_set_error_detail(ALEA_ERR_INVALID_STATE,"No root universe 0");goto fail;}
    if(root>=0 && check_universe(sys,univ,nu,order,root,0))goto fail;
    /* Point queries share programs and transforms with rendering, but must
     * enumerate terminal void and all overlapping occurrences. No spatial
     * pruning is applied until its FP32 ownership contract is validated. */
    scene->point_root = COMPUTE_NONE;
    if (nl || nt || nu > 1) {
        scene->point_cells = calloc(n, sizeof(*scene->point_cells));
        if (!scene->point_cells) goto memory;
        if (root >= 0) scene->point_root = order[univ[root].first].cell;
        for (size_t u = 0; u < nu; ++u) {
            for (size_t j = univ[u].first; j < univ[u].first + univ[u].count; ++j) {
                uint32_t i = order[j].cell;
                metal_point_cell* point = &scene->point_cells[i];
                point->next = j + 1 < univ[u].first + univ[u].count
                    ? order[j + 1].cell : COMPUTE_NONE;
                point->child = COMPUTE_NONE;
                if (univ[u].state == 2 && !sys->cells.data[i].lat_type && sys->cells.data[i].fill_universe > 0) {
                    int child = find_universe(univ, nu, sys->cells.data[i].fill_universe);
                    point->child = order[univ[child].first].cell;
                }
            }
        }
    }
    /* Build only root-reachable universes. Unused input definitions are not
     * instances and must not appear in a rendered frame. */
    for(size_t u=0;u<nu;++u) {
        if(univ[u].state!=2)continue;
        size_t visible=0;
        for(size_t j=univ[u].first;j<univ[u].first+univ[u].count;++j) {
            uint32_t i=order[j].cell;const alea_cell_entry_t* c=&sys->cells.data[i];
            if(!c->material_id && !alea_cell_entry_is_container(c))continue;
            bvh_item* item=&items[visible++];item->cell=i;
            memcpy(item->bounds,bounds[i],sizeof(item->bounds));
            for(int a=0;a<3;++a){double center=.5*((double)item->bounds[2*a]+item->bounds[2*a+1]);item->center[a]=isfinite(center)?center:0;}
        }
        if(visible){univ[u].root=(uint32_t)out->bvh_count;build_bvh(items,visible,out);}
    }
    if(root>=0)out->root=univ[root].root;
    out->transforms=calloc(nt?nt:1,sizeof(*out->transforms));
    out->lattices=calloc(nl?nl:1,sizeof(*out->lattices));
    out->fills=calloc(nf?nf:1,sizeof(*out->fills));
    scene->point_lattices=calloc(nl?nl:1,sizeof(*scene->point_lattices));
    scene->point_fills=calloc(nf?nf:1,sizeof(*scene->point_fills));
    if(!out->transforms || !out->lattices || !out->fills || !scene->point_lattices || !scene->point_fills)goto memory;
    for(size_t i=0;i<n;++i) {
        const alea_cell_entry_t* c=&sys->cells.data[i];
        int u=find_universe(univ,nu,c->universe_id);if(u<0 || univ[u].state!=2)continue;
        compute_cell* dst=&scene->cells[i];
        if(c->lat_type) {
            compute_lattice* l=&out->lattices[out->lattice_count];
            dst->lattice=(uint32_t)out->lattice_count++;
            memcpy(l->dims,c->lat_fill_dims,sizeof(l->dims));
            for(int a=0;a<3;++a){l->pitch[a]=c->lat_pitch[a];l->lower[a]=c->lat_lower_left[a];
                if(!isfinite(l->pitch[a]) || !isfinite(l->lower[a]) || (c->lat_pitch[a]>0 && l->pitch[a]==0)) {
                    alea_set_error_detail(ALEA_ERR_UNSUPPORTED,"Lattice coordinates do not fit FP32");goto fail;
                }}
            l->offset=(uint32_t)out->fill_count;l->repeating=c->lat_fill_repeating;l->zero_coords=c->lat_fill_zero_element_coords;l->outer=COMPUTE_NONE;
            for(size_t j=0;j<c->lat_fill_count;++j) {
                int child=c->lat_fill[j]>0?find_universe(univ,nu,c->lat_fill[j]):-1;
                scene->point_fills[out->fill_count]=child>=0?order[univ[child].first].cell:COMPUTE_NONE;
                out->fills[out->fill_count++]=child>=0?univ[child].root:COMPUTE_NONE;
            }
            int child=c->lat_outer_universe>0?find_universe(univ,nu,c->lat_outer_universe):-1;if(child>=0)l->outer=univ[child].root;
            scene->point_lattices[dst->lattice]=*l;
            scene->point_lattices[dst->lattice].outer=child>=0?order[univ[child].first].cell:COMPUTE_NONE;
        } else if(c->fill_universe>0) {
            int child=find_universe(univ,nu,c->fill_universe);
            dst->child=univ[child].root;
            alea_matrix_t matrix;alea_matrix_identity(&matrix);
            if(c->fill_transform>0) {
                const alea_transform_t* tr=alea_get_transform(sys,c->fill_transform);
                if(!tr || !alea_matrix_from_mcnp(&matrix,tr->cosines,tr->value_count,false) || !alea_matrix_invert(&matrix)) {
                    alea_set_error_detail(ALEA_ERR_INVALID_STATE,"Invalid fill transform in cell %d",c->mc_cell_id);goto fail;
                }
            }
            dst->transform=(uint32_t)out->transform_count;
            compute_transform* t=&out->transforms[out->transform_count++];
            for(int j=0;j<12;++j){t->inverse[j]=matrix.inv[j];if(!isfinite(t->inverse[j])) {
                alea_set_error_detail(ALEA_ERR_UNSUPPORTED,"Fill transform does not fit FP32");goto fail;
            }}
        }
    }
    scene->point_lattice_count=out->lattice_count;
    scene->point_fill_count=out->fill_count;
    if(out->bvh_count) {
        uint32_t* memo=malloc(out->bvh_count*sizeof(*memo));if(!memo)goto memory;
        memset(memo,255,out->bvh_count*sizeof(*memo));
        for(size_t i=0;i<out->bvh_count;++i)mark_cost(scene,(uint32_t)i,memo);
        free(memo);
    }
    free(items);free(order);free(univ);free(bounds);return 0;
memory:
    alea_set_error_detail(ALEA_ERR_OUT_OF_MEMORY,"Preparing Metal hierarchy");
fail:
    free(items);free(order);free(univ);free(bounds);alea_compute_render_free(out);return -1;
}
