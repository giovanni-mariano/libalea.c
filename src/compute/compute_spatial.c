// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
// SPDX-License-Identifier: MPL-2.0
#include "compute_internal.h"
#include "core/alea_universe.h"
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Render-only acceleration. Original point-classification programs are never
 * modified. Shared universe BVHs stay shared; repeating lattices stay procedural. */
#define SPATIAL_NODES (256u*1024u)
#define SPATIAL_CODE (8u*1024u*1024u)
#define SPATIAL_DEPTH 6
_Static_assert(sizeof(compute_spatial_program)==40,"Spatial program layout");
typedef struct {
    compute_scene* s;
    unsigned char* visited;
    uint32_t* universe;
    float (*hints)[6];
    size_t original,code_capacity,tree_capacity,tightened,compiled,empty,tree_limit;
} compiler;
static void empty_bounds(float b[6]) {for(int a=0;a<3;++a){b[2*a]=INFINITY;b[2*a+1]=-INFINITY;}}
static int empty(const float b[6]) {return b[0]>b[1] || b[2]>b[3] || b[4]>b[5];}
static void unite(float b[6],const float a[6]) {
    for(int k=0;k<3;++k){b[2*k]=fminf(b[2*k],a[2*k]);b[2*k+1]=fmaxf(b[2*k+1],a[2*k+1]);}
}
static int finite_bounds(const float b[6]) {
    for(int a=0;a<3;++a)if(!isfinite(b[2*a]) || !isfinite(b[2*a+1]) ||
        fabsf(b[2*a])>=1e9f || fabsf(b[2*a+1])>=1e9f || b[2*a]>=b[2*a+1])return 0;
    return 1;
}
/* Invert the actual FP32 device matrix, not the source FP64 transform. Unknown
 * or ill-conditioned transforms cannot tighten a box. Leave substantial slack
 * for transformed ray arithmetic as well as outward rounding of the bounds. */
static int placed_bounds(const compute_transform* t,const float child[6],float out[6],float* guard) {
    if(empty(child)){empty_bounds(out);if(guard)*guard=0;return 1;}
    if(!finite_bounds(child))return 0;
    alea_matrix_t m={0};for(int j=0;j<12;++j)m.m[j]=t->inverse[j];
    if(!alea_matrix_invert(&m))return 0;
    double norm=0,inorm=0;
    for(int a=0;a<3;++a) {
        double row=0,ir=0;
        for(int j=0;j<3;++j){row+=fabs(m.m[4*a+j]);ir+=fabs(m.inv[4*a+j]);}
        norm=fmax(norm,row);inorm=fmax(inorm,ir);
    }
    if(norm*inorm>16)return 0;
    double min_pad=INFINITY;
    for(int a=0;a<3;++a) {
        double lo=m.inv[4*a+3],hi=lo,scale=fabs(lo);
        for(int j=0;j<3;++j) {
            double x=m.inv[4*a+j]*child[2*j],y=m.inv[4*a+j]*child[2*j+1];
            lo+=fmin(x,y);hi+=fmax(x,y);scale+=fmax(fabs(x),fabs(y));
        }
        double pad=fmax(.01,scale*.001);min_pad=fmin(min_pad,pad);
        out[2*a]=nextafterf((float)(lo-pad),-INFINITY);
        out[2*a+1]=nextafterf((float)(hi+pad),INFINITY);
    }
    if(guard)*guard=(float)(min_pad/(64*FLT_EPSILON*norm*inorm));
    return finite_bounds(out);
}
static int refit(compiler* c,uint32_t id,int contents,uint32_t universe) {
    compute_render_scene* r=&c->s->render;
    if(id==COMPUTE_NONE || c->visited[id])return 0;
    if(alea_interrupted()){alea_set_error_detail(ALEA_ERR_INTERRUPTED,"Spatial hierarchy compilation interrupted");return -1;}
    c->visited[id]=1;c->universe[id]=universe;compute_render_bvh* n=&r->bvh[id];
    if(n->cell==COMPUTE_NONE) {
        uint32_t left=id+1,right=r->bvh[left].escape;
        if(refit(c,left,contents,universe) || refit(c,right,contents,universe))return -1;
        memcpy(n->bounds,r->bvh[left].bounds,sizeof(n->bounds));unite(n->bounds,r->bvh[right].bounds);
        n->guard=fminf(r->bvh[left].guard,r->bvh[right].guard);
    } else {
        compute_cell* cell=&c->s->cells[n->cell];
        if(cell->lattice!=COMPUTE_NONE) {
            const compute_lattice* l=&r->lattices[cell->lattice];
            size_t count=l->repeating?1:(size_t)(l->dims[1]-l->dims[0]+1)*(l->dims[3]-l->dims[2]+1)*(l->dims[5]-l->dims[4]+1);
            for(size_t i=0;i<count;++i)if(refit(c,r->fills[l->offset+i],contents,r->fills[l->offset+i]))return -1;
            if(refit(c,l->outer,contents,l->outer))return -1;
        } else if(cell->child!=COMPUTE_NONE) {
            if(refit(c,cell->child,contents,cell->child))return -1;
            float b[6],guard;
            if(contents && placed_bounds(&r->transforms[cell->transform],r->bvh[cell->child].bounds,b,&guard)) {
                int changed=0;
                for(int a=0;a<3;++a) {
                    float lo=fmaxf(n->bounds[2*a],b[2*a]),hi=fminf(n->bounds[2*a+1],b[2*a+1]);
                    changed|=lo>n->bounds[2*a] || hi<n->bounds[2*a+1];
                    n->bounds[2*a]=lo;n->bounds[2*a+1]=hi;
                }
                if(changed || memcmp(n->bounds,n->fallback,sizeof(n->bounds))) {
                    const float* m=r->transforms[cell->transform].inverse;float norm=0,shift=0,child_scale=0;
                    for(int a=0;a<3;++a) {
                        norm=fmaxf(norm,fabsf(m[4*a])+fabsf(m[4*a+1])+fabsf(m[4*a+2]));
                        shift=fmaxf(shift,fabsf(m[4*a+3]));
                        child_scale=fmaxf(child_scale,fmaxf(fabsf(r->bvh[cell->child].bounds[2*a]),fabsf(r->bvh[cell->child].bounds[2*a+1])));
                    }
                    float child_guard=r->bvh[cell->child].guard;
                    if(isfinite(child_guard))guard=fminf(guard,fmaxf(0,(child_guard-8*(shift+child_scale+1))/fmaxf(1,norm)));
                    n->guard=fminf(n->guard,guard);c->tightened+=changed;
                }
            }
        }
    }
    if(empty(n->bounds))empty_bounds(n->bounds);
    return 0;
}
/* These hints choose where to cache proofs; they are NOT geometry bounds.
 * A universe with unbounded cells can still optimize its finite components.
 * The evaluator uses the complete original program outside the chosen domain. */
static void domain_hint(compiler* c,uint32_t id) {
    if(id==COMPUTE_NONE || c->visited[id])return;
    c->visited[id]=1;float* hint=c->hints[id];empty_bounds(hint);
    compute_render_bvh* n=&c->s->render.bvh[id];
    if(finite_bounds(n->bounds)){memcpy(hint,n->bounds,6*sizeof(float));return;}
    if(n->cell==COMPUTE_NONE) {
        uint32_t left=id+1,right=c->s->render.bvh[left].escape;
        domain_hint(c,left);domain_hint(c,right);unite(hint,c->hints[left]);unite(hint,c->hints[right]);
    } else {
        compute_cell cell=c->s->cells[n->cell];
        if(cell.child!=COMPUTE_NONE) {
            domain_hint(c,cell.child);
            placed_bounds(&c->s->render.transforms[cell.transform],c->hints[cell.child],hint,NULL);
        }
    }
}
typedef struct {uint32_t universe;float bounds[6];} domain_sample;
static int sample_order(const void* a,const void* b) {
    uint32_t x=((const domain_sample*)a)->universe,y=((const domain_sample*)b)->universe;return (x>y)-(x<y);
}
static int float_order(const void* a,const void* b) {float x=*(const float*)a,y=*(const float*)b;return (x>y)-(x<y);}
/* Very large construction solids are common in MCNP. A trimmed distribution
 * of finite cell bounds selects a useful cache domain without treating those
 * outliers as the model's physical extent. Outside it, original CSG still runs. */
static int typical_domains(compiler* c) {
    size_t n=0;compute_render_scene* r=&c->s->render;
    domain_sample* samples=malloc(r->bvh_count*sizeof(*samples));
    float* values=malloc(r->bvh_count*sizeof(*values));
    if(!samples || !values){free(samples);free(values);return -1;}
    for(size_t i=0;i<r->bvh_count;++i)if(r->bvh[i].cell!=COMPUTE_NONE && finite_bounds(r->bvh[i].bounds)) {
        samples[n].universe=c->universe[i];memcpy(samples[n++].bounds,r->bvh[i].bounds,6*sizeof(float));
    }
    qsort(samples,n,sizeof(*samples),sample_order);
    for(size_t first=0;first<n;) {
        size_t end=first+1;while(end<n && samples[end].universe==samples[first].universe)++end;
        size_t count=end-first;
        if(count>=4) {
            float b[6];
            for(int a=0;a<6;++a) {
                for(size_t j=0;j<count;++j)values[j]=samples[first+j].bounds[a];
                qsort(values,count,sizeof(*values),float_order);
                size_t k=(a&1)?(count-1)*9/10:(count-1)/10;b[a]=values[k];
            }
            if(finite_bounds(b)) {
                for(int a=0;a<3;++a){float pad=(b[2*a+1]-b[2*a])*.5f;b[2*a]-=pad;b[2*a+1]+=pad;}
                memcpy(c->hints[samples[first].universe],b,sizeof(b));
            }
        }
        first=end;
    }
    free(samples);free(values);return 0;
}
static int reserve_tree(compiler* c,size_t n) {
    compute_scene* s=c->s;
    if(s->spatial_count+n>c->tree_limit)return 1;
    if(s->spatial_count+n<=c->tree_capacity)return 0;
    size_t capacity=s->spatial_count+n+1024;if(capacity>c->tree_limit)capacity=c->tree_limit;
    void* p=realloc(s->spatial,capacity*sizeof(*s->spatial));
    if(!p){alea_set_error_detail(ALEA_ERR_OUT_OF_MEMORY,"Allocating spatial regions");return -1;}
    s->spatial=p;c->tree_capacity=capacity;return 0;
}
static int reduce(compiler* c,compute_cell cell,const float b[6],compute_node* code,size_t* count) {
    float proof[6];
    for(int a=0;a<3;++a) {
        float margin=fmaxf((b[2*a+1]-b[2*a])*.001f,fmaxf(fabsf(b[2*a]),fabsf(b[2*a+1]))*1e-5f);
        proof[2*a]=nextafterf(b[2*a]-margin,-INFINITY);proof[2*a+1]=nextafterf(b[2*a+1]+margin,INFINITY);
    }
    return alea_compute_local_simplify(c->s->nodes+cell.root,cell.count,c->s->primitives,proof,code,count);
}
static int save_program(compiler* c,uint32_t id,compute_cell cell,const compute_node* code,size_t count) {
    compute_scene* s=c->s;compute_spatial_program* p=&s->spatial[id];p->root=cell.root;p->count=cell.count;
    if((count>=cell.count && !(count==1 && code[0].op==COMPUTE_CONSTANT)) || count>SPATIAL_CODE-(s->node_count-c->original) || count>UINT32_MAX-s->node_count)return 0;
    if(s->node_count+count>c->code_capacity) {
        size_t cap=s->node_count+count+65536;
        compute_node* grown=realloc(s->nodes,cap*sizeof(*grown));
        if(!grown){alea_set_error_detail(ALEA_ERR_OUT_OF_MEMORY,"Allocating spatial programs");return -1;}
        s->nodes=grown;c->code_capacity=cap;
    }
    p->root=(uint32_t)s->node_count;p->count=(uint32_t)count;
    memcpy(s->nodes+s->node_count,code,count*sizeof(*code));s->node_count+=count;return 0;
}
/* The first three levels expose spatial locality even when one split alone
 * proves nothing. Further splits must reduce mean residual size by 10%. Uniform
 * regions stop immediately. Depth and two independent memory budgets bound
 * compilation even for pathological or deliberately enormous expressions. */
static int build_region(compiler* c,compute_cell cell,uint32_t id,unsigned depth,compute_node* code,size_t count) {
    compute_scene* s=c->s;
    if(alea_interrupted()){alea_set_error_detail(ALEA_ERR_INTERRUPTED,"Spatial CSG compilation interrupted");return -1;}
    if(count<32 || depth==SPATIAL_DEPTH || s->node_count-c->original>=SPATIAL_CODE || s->spatial_count+2>c->tree_limit)
        return save_program(c,id,cell,code,count);
    float box[6],b[2][6];memcpy(box,s->spatial[id].bounds,sizeof(box));
    unsigned axis=0;for(unsigned a=1;a<3;++a)if(box[2*a+1]-box[2*a]>box[2*axis+1]-box[2*axis])axis=a;
    compute_node* scratch=malloc(4*(size_t)cell.count*sizeof(*scratch));
    if(!scratch){alea_set_error_detail(ALEA_ERR_OUT_OF_MEMORY,"Splitting spatial CSG region");return -1;}
    size_t n[2]={0,0},best=SIZE_MAX;int rc=-1;
    unsigned first_axis=axis;
    for(unsigned trial=0;trial<3;++trial) {
        unsigned a=(first_axis+trial)%3;
        float mid=box[2*a]+(box[2*a+1]-box[2*a])*.5f;
        if(mid<=box[2*a] || mid>=box[2*a+1])continue;
        float trial_bounds[2][6];memcpy(trial_bounds[0],box,sizeof(box));memcpy(trial_bounds[1],box,sizeof(box));
        trial_bounds[0][2*a+1]=trial_bounds[1][2*a]=mid;
        size_t sizes[2];compute_node* trial_code=scratch+2*cell.count;
        if(reduce(c,cell,trial_bounds[0],trial_code,&sizes[0]) ||
            reduce(c,cell,trial_bounds[1],trial_code+cell.count,&sizes[1]))goto done;
        if(sizes[0]+sizes[1]<best) {
            best=sizes[0]+sizes[1];axis=a;n[0]=sizes[0];n[1]=sizes[1];memcpy(b,trial_bounds,sizeof(b));
            memcpy(scratch,trial_code,n[0]*sizeof(*scratch));
            memcpy(scratch+cell.count,trial_code+cell.count,n[1]*sizeof(*scratch));
        }
    }
    if(best==SIZE_MAX || (depth>=3 && best>=count*1.8)) {rc=save_program(c,id,cell,code,count);goto done;}
    if(reserve_tree(c,2))goto done;
    uint32_t child=(uint32_t)s->spatial_count;s->spatial_count+=2;
    s->spatial[id].child=child;s->spatial[id].axis=axis;
    for(int j=0;j<2;++j) {
        s->spatial[child+j]=(compute_spatial_program){.child=COMPUTE_NONE};
        memcpy(s->spatial[child+j].bounds,b[j],sizeof(b[j]));
        if(build_region(c,cell,child+j,depth+1,scratch+j*cell.count,n[j]))goto done;
    }
    rc=0;
done:free(scratch);return rc;
}
static void occupied_bounds(const compute_scene* s,uint32_t id,float b[6],float* margin) {
    const compute_spatial_program* p=&s->spatial[id];
    if(p->child!=COMPUTE_NONE){occupied_bounds(s,p->child,b,margin);occupied_bounds(s,p->child+1,b,margin);}
    else {
        for(int a=0;a<3;++a)*margin=fminf(*margin,fmaxf((p->bounds[2*a+1]-p->bounds[2*a])*.001f,
            fmaxf(fabsf(p->bounds[2*a]),fabsf(p->bounds[2*a+1]))*1e-5f));
        if(!(p->count==1 && s->nodes[p->root].op==COMPUTE_CONSTANT && s->nodes[p->root].arg))unite(b,p->bounds);
    }
}
int alea_compute_spatial_pack(compute_scene* s) {
    const char* env=getenv("ALEA_METAL_SPATIAL");int mode=2;
    if(env){if(!strcmp(env,"0"))mode=0;else if(!strcmp(env,"bounds"))mode=1;else if(strcmp(env,"1")) {
        alea_set_error_detail(ALEA_ERR_INVALID_ARG,"ALEA_METAL_SPATIAL must be 0, bounds, or 1");return -1;
    }}
    if(!mode || s->render.root==COMPUTE_NONE)return 0;
    if(s->render.root>=s->render.bvh_count){alea_set_error_detail(ALEA_ERR_INVALID_STATE,"Invalid spatial BVH root");return -1;}
    compiler c={.s=s,.original=s->node_count,.code_capacity=s->node_count,.tree_limit=SPATIAL_NODES};
    c.visited=calloc(s->render.bvh_count,1);
    c.universe=calloc(s->render.bvh_count,sizeof(*c.universe));
    c.hints=calloc(s->render.bvh_count,sizeof(*c.hints));
    if(!c.visited || !c.universe || !c.hints)goto memory;
    if(refit(&c,s->render.root,1,s->render.root))goto fail;
    memset(c.visited,0,s->render.bvh_count);
    for(size_t i=0;i<s->render.bvh_count;++i)domain_hint(&c,(uint32_t)i);
    if(mode==2 && typical_domains(&c))goto memory;
    if(mode==2)for(size_t i=0;i<s->render.bvh_count;++i) {
        compute_render_bvh* b=&s->render.bvh[i];if(b->cell==COMPUTE_NONE || empty(b->bounds))continue;
        compute_cell cell=s->cells[b->cell];
        if(cell.child==COMPUTE_NONE || cell.lattice!=COMPUTE_NONE || cell.count>262144)continue;
        int complete=1;float domain[6];
        const float* hint=c.hints[c.universe[i]];
        if(!finite_bounds(hint))hint=c.hints[i];
        memcpy(domain,b->bounds,sizeof(domain));
        if(finite_bounds(hint))for(int a=0;a<3;++a) {
            domain[2*a]=fmaxf(domain[2*a],hint[2*a]);domain[2*a+1]=fminf(domain[2*a+1],hint[2*a+1]);
        }
        if(!finite_bounds(domain))memcpy(domain,b->bounds,sizeof(domain));
        if(!finite_bounds(domain))continue;
        for(int a=0;a<6;++a)if(domain[a]!=b->bounds[a])complete=0;
        if(s->spatial_count+1>c.tree_limit || s->node_count-c.original>=SPATIAL_CODE)break;
        size_t old_tree=s->spatial_count,old_code=s->node_count;
        if(reserve_tree(&c,1))goto memory;
        uint32_t root=(uint32_t)s->spatial_count++;
        s->spatial[root]=(compute_spatial_program){.child=COMPUTE_NONE};memcpy(s->spatial[root].bounds,domain,sizeof(domain));
        compute_node* temp=malloc((size_t)cell.count*sizeof(*temp));size_t count;
        if(!temp)goto memory;
        int rc=reduce(&c,cell,domain,temp,&count);
        if(!rc)rc=build_region(&c,cell,root,0,temp,count);
        free(temp);if(rc)goto fail;
        if(s->node_count==old_code){s->spatial_count=old_tree;continue;}
        s->cells[b->cell].spatial=root;++c.compiled;
        if(complete) {
            float occupied[6],margin=INFINITY;empty_bounds(occupied);occupied_bounds(s,root,occupied,&margin);
            if(empty(occupied))++c.empty;
            if(memcmp(b->bounds,occupied,sizeof(b->bounds))) {
                b->guard=fminf(b->guard,margin/(64*FLT_EPSILON));
                memcpy(b->bounds,occupied,sizeof(b->bounds));
            }
        }
    }
    memset(c.visited,0,s->render.bvh_count);
    if(refit(&c,s->render.root,1,s->render.root))goto fail;
    if(c.tightened || c.compiled)fprintf(stderr,"Metal spatial CSG: %zu content bounds tightened, %zu containers specialized (%zu empty), %zu regions, %.2f MiB code, %.2f MiB tree\n",
        c.tightened,c.compiled,c.empty,s->spatial_count,(s->node_count-c.original)*sizeof(compute_node)/1048576.0,s->spatial_count*sizeof(*s->spatial)/1048576.0);
    free(c.visited);free(c.universe);free(c.hints);return 0;
memory:alea_set_error_detail(ALEA_ERR_OUT_OF_MEMORY,"Compiling spatial render hierarchy");
fail:free(c.visited);free(c.universe);free(c.hints);return -1;
}
