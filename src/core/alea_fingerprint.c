// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
// SPDX-License-Identifier: MPL-2.0
#include "alea_fingerprint.h"
#include "alea_system.h"
#include "alea_transform.h"
#include "primitives/primitive_codec.h"
#include "util/alea_digest.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
typedef struct {
    const alea_system_t* sys;
    alea_fingerprint_profile_t profile;
    alea_fingerprints_t* out;
    alea_digest_t* nodes;
    unsigned char *node_state,*cell_state,*universe_state;
    size_t work;
    size_t *first,*next;
    alea_digest_t *materials,*mixtures;
    unsigned char *material_ready,*mixture_ready;
}
context;
static int cmp_digest(const void* a,const void* b) {
    return memcmp(a,b,32);
}
static int cmp_fp(const void* a,const void* b) {
    int x=((const alea_fingerprint_t*)a)->id,y=((const alea_fingerprint_t*)b)->id;
    return (x>y)-(x<y);
}
static void begin(context* c,alea_sha256_t* h,unsigned tag) {
    alea_sha256_init(h);
    alea_digest_string(h,"ALEA fingerprint");
    alea_digest_u64(h,c->profile==ALEA_FINGERPRINT_LEGACY?ALEA_FINGERPRINT_VERSION:3);
    alea_digest_u64(h,tag);
}
static size_t universe_index(const alea_fingerprints_t* f,int id) {
    size_t lo=0,hi=f->universe_count;
    while(lo<hi) {
        size_t mid=lo+(hi-lo)/2;
        if(f->universes[mid].id<id)lo=mid+1;
        else hi=mid;
    }
    return lo<f->universe_count && f->universes[lo].id==id?lo:f->universe_count;
}
static int node_hash(context* c,uint32_t id,unsigned depth);
static int gather(context* c,uint32_t id,alea_operation_t op,alea_digest_t** list,size_t* n,size_t* cap,unsigned depth) {
    if(id>=c->sys->nodes.count || depth>1024 || ++c->work>100000000 || alea_interrupted())return -1;
    const alea_node_t* node=&c->sys->nodes.data[id];
    if(ALEA_GET_OPERATION(node)==op) {
        if(gather(c,node->operation.left,op,list,n,cap,depth+1))return -1;
        return gather(c,node->operation.right,op,list,n,cap,depth+1);
    }
    if(node_hash(c,id,depth+1))return -1;
    if(*n==*cap) {
        size_t size=*cap?*cap*2:8;
        void* p=realloc(*list,size*sizeof(**list));
        if(!p)return -1;
        *list=p;
        *cap=size;
    }
    (*list)[(*n)++]=c->nodes[id];
    return 0;
}
static int node_hash(context* c,uint32_t id,unsigned depth) {
    if(id>=c->sys->nodes.count || depth>1024 || alea_interrupted())return -1;
    if(c->node_state[id]==2)return 0;
    if(c->node_state[id])return -1;
    c->node_state[id]=1;
    const alea_node_t* node=&c->sys->nodes.data[id];
    alea_operation_t op=ALEA_GET_OPERATION(node);
    alea_sha256_t h;
    begin(c,&h,10+op);
    if(op==ALEA_OP_PRIMITIVE) {
        alea_primitive_data_t data;
        double values[50];
        if(!alea_primitive_copy_data(c->sys,node->primitive.primitive_id,&data))return -1;
        size_t n=alea_primitive_to_values(node->primitive.prim_type,&data,values);
        if(!n)return -1;
        int sense=node->primitive.sense*(node->primitive.inverted?-1:1);
        if(node->primitive.prim_type==ALEA_PRIMITIVE_PLANE || node->primitive.prim_type==ALEA_PRIMITIVE_QUADRIC) {
            double scale=0;
            for(size_t i=0;i<n;i++)if(fabs(values[i])>fabs(scale))scale=values[i];
            if(scale==0)return -1;
            if(scale<0)sense=-sense;
            for(size_t i=0;i<n;i++)values[i]/=scale;
        }
        alea_digest_u64(&h,node->primitive.prim_type);
        alea_digest_u64(&h,(uint64_t)(int64_t)sense);
        for(size_t i=0;i<n;i++)if(alea_digest_double(&h,values[i]))return -1;
    }
    else if(op==ALEA_OP_UNION || op==ALEA_OP_INTERSECTION) {
        alea_digest_t* list=NULL;
        size_t n=0,cap=0;
        if(gather(c,node->operation.left,op,&list,&n,&cap,depth+1)||gather(c,node->operation.right,op,&list,&n,&cap,depth+1)) {
            free(list);
            return -1;
        }
        qsort(list,n,sizeof(*list),cmp_digest);
        size_t unique=0;
        for(size_t i=0;i<n;i++)if(i==0 || memcmp(&list[i],&list[i-1],32)) {
            alea_sha256_update(&h,&list[i],32);
            unique++;
        }
        /* Idempotence also removes the operator when only one term remains. */
        if(unique==1) {
            c->nodes[id]=list[0];
            c->node_state[id]=2;
            free(list);
            return 0;
        }
        alea_digest_u64(&h,unique);
        free(list);
    }
    else if(op==ALEA_OP_DIFFERENCE || op==ALEA_OP_COMPLEMENT) {
        if(node_hash(c,node->operation.left,depth+1))return -1;
        alea_sha256_update(&h,&c->nodes[node->operation.left],32);
        if(op==ALEA_OP_DIFFERENCE) {
            if(node_hash(c,node->operation.right,depth+1))return -1;
            alea_sha256_update(&h,&c->nodes[node->operation.right],32);
        }
    }
    else return -1;
    alea_sha256_final(&h,c->nodes[id].bytes);
    c->node_state[id]=2;
    return 0;
}
static int component_hash(context* c,alea_digest_t* digest,int kind,int id,const char* library,double fraction) {
    alea_sha256_t h;
    begin(c,&h,40+kind);
    alea_digest_u64(&h,id);
    alea_digest_string(&h,library);
    if(alea_digest_double(&h,fraction))return -1;
    alea_sha256_final(&h,digest->bytes);
    return 0;
}
typedef struct {
    alea_digest_t identity;
    double fraction;
} composition_part;
static int cmp_part(const void* a,const void* b) {
    const composition_part* x=a;const composition_part* y=b;
    int cmp=memcmp(&x->identity,&y->identity,32);
    return cmp?cmp:(x->fraction>y->fraction)-(x->fraction<y->fraction);
}
/* Scaling before summation avoids overflow even for DBL_MAX fractions.
 * Sorting makes duplicate accumulation independent of source order. */
static int normalized_parts(context* c,composition_part* parts,size_t count,
                            alea_digest_t* out,size_t* written) {
    double scale=0;
    *written=0;
    for(size_t i=0;i<count;i++) {
        if(!isfinite(parts[i].fraction) || parts[i].fraction<0)return -1;
        if(parts[i].fraction>scale)scale=parts[i].fraction;
    }
    if(!count)return 0;
    if(scale==0)return -1;
    qsort(parts,count,sizeof(*parts),cmp_part);
    long double total=0;
    for(size_t i=0;i<count;i++)total+=(long double)parts[i].fraction/scale;
    for(size_t i=0;i<count;) {
        size_t end=i;long double sum=0;
        while(end<count && !memcmp(&parts[i].identity,&parts[end].identity,32)) {
            sum+=(long double)parts[end].fraction/scale;end++;
        }
        if(sum>0) {
            double fraction=(double)(sum/total);
            /* Never silently drop a positive component through underflow. */
            if(fraction==0 || !isfinite(fraction))return -1;
            alea_sha256_t h;begin(c,&h,72);
            alea_sha256_update(&h,&parts[i].identity,32);
            if(alea_digest_double(&h,fraction))return -1;
            alea_sha256_final(&h,out[(*written)++].bytes);
        }
        i=end;
    }
    return 0;
}
static int semantic_material_parts(context* c,const alea_material_t* m,
                                   alea_digest_t* out,size_t* written) {
    size_t count=m->nuclides.count+m->elements.count;
    composition_part* parts=calloc(count?count:1,sizeof(*parts));
    if(!parts)return -1;
    size_t n=0;int rc=-1;
    for(size_t i=0;i<m->nuclides.count;i++) {
        const alea_nuclide_t* v=&m->nuclides.data[i];
        if(component_hash(c,&parts[n].identity,0,v->zaid,v->library,0))goto done;
        parts[n++].fraction=v->fraction;
    }
    for(size_t i=0;i<m->elements.count;i++) {
        const alea_element_comp_t* v=&m->elements.data[i];
        if(component_hash(c,&parts[n].identity,1,v->atomic_number,v->library,0))goto done;
        parts[n++].fraction=v->fraction;
    }
    rc=normalized_parts(c,parts,n,out,written);
done:
    free(parts);return rc;
}
static int material_hash(context* c,const alea_material_t* m,alea_digest_t* digest) {
    alea_sha256_t h;
    begin(c,&h,30);
    alea_digest_u64(&h,m->is_weight_fraction);
    alea_digest_u64(&h,m->has_standard_density);
    if(m->has_standard_density && alea_digest_double(&h,m->standard_density))return -1;
    size_t count=m->nuclides.count+m->elements.count+m->thermal_laws.count;
    alea_digest_t* parts=calloc(count?count:1,sizeof(*parts));
    if(!parts)return -1;
    size_t n=0;
    int rc=-1;
    if(c->profile==ALEA_FINGERPRINT_SEMANTIC) {
        if(semantic_material_parts(c,m,parts,&n))goto done;
    } else {
    /* Legacy fractions are intentionally preserved. */
    for(size_t i=0;i<m->nuclides.count;i++) {
        const alea_nuclide_t* v=&m->nuclides.data[i];
        if(component_hash(c,&parts[n++],0,v->zaid,v->library,v->fraction))goto done;
    }
    for(size_t i=0;i<m->elements.count;i++) {
        const alea_element_comp_t* v=&m->elements.data[i];
        if(component_hash(c,&parts[n++],1,v->atomic_number,v->library,v->fraction))goto done;
    }
    }
    for(size_t i=0;i<m->thermal_laws.count;i++) {
        const alea_thermal_law_t* v=&m->thermal_laws.data[i];
        if(component_hash(c,&parts[n++],2,v->zaid_match,v->identifier,0))goto done;
    }
    qsort(parts,n,sizeof(*parts),cmp_digest);
    alea_digest_u64(&h,n);
    alea_sha256_update(&h,parts,n*sizeof(*parts));
    alea_sha256_final(&h,digest->bytes);
    rc=0;
    done:free(parts);
    return rc;
}
static int cached_material(context* c,size_t index,alea_digest_t* out) {
    if(index>=c->sys->materials.count)return -1;
    if(!c->material_ready[index]) {
        if(material_hash(c,&c->sys->materials.data[index],&c->materials[index]))return -1;
        c->material_ready[index]=1;
    }
    *out=c->materials[index];
    return 0;
}
static int mixture_hash(context* c,int id,alea_digest_t* out) {
    size_t index=0;
    while(index<c->sys->mixtures.count && c->sys->mixtures.data[index].mixture_id!=id)index++;
    if(index==c->sys->mixtures.count)return -1;
    if(c->mixture_ready[index]) {
        *out=c->mixtures[index];
        return 0;
    }
    const alea_mixture_t* mix=&c->sys->mixtures.data[index];
    alea_digest_t* parts=calloc(mix->components.count?mix->components.count:1,sizeof(*parts));
    if(!parts)return -1;
    int status=-1;
    composition_part* normalized=calloc(mix->components.count?mix->components.count:1,sizeof(*normalized));
    if(!normalized){free(parts);return -1;}
    for(size_t i=0;i<mix->components.count;i++) {
        int material=alea_find_material_by_id(c->sys,mix->components.data[i].material_id);
        alea_digest_t base;
        alea_sha256_t hash;
        begin(c,&hash,70);
        if(material<0 || cached_material(c,(size_t)material,&base))goto done;
        normalized[i].identity=base;
        normalized[i].fraction=mix->components.data[i].fraction;
        alea_sha256_update(&hash,&base,32);
        if(alea_digest_double(&hash,mix->components.data[i].fraction))goto done;
        alea_sha256_final(&hash,parts[i].bytes);
    }
    size_t count=mix->components.count;
    if(c->profile==ALEA_FINGERPRINT_SEMANTIC && normalized_parts(c,normalized,count,parts,&count))goto done;
    qsort(parts,count,sizeof(*parts),cmp_digest);
    alea_sha256_t hash;
    begin(c,&hash,71);
    alea_digest_u64(&hash,mix->is_weight_fraction);
    alea_digest_u64(&hash,count);
    alea_sha256_update(&hash,parts,count*sizeof(*parts));
    alea_sha256_final(&hash,c->mixtures[index].bytes);
    c->mixture_ready[index]=1;
    *out=c->mixtures[index];
    status=0;
    done:free(normalized);free(parts);
    return status;
}
static int universe_hash(context* c,int id,unsigned depth);
static int target(context* c,alea_sha256_t* g,alea_sha256_t* h,int id,int self,unsigned depth) {
    if(id==self) {
        alea_digest_u64(g,1);
        alea_digest_u64(h,1);
        return 0;
    }
    if(id<=0) {
        alea_digest_u64(g,0);
        alea_digest_u64(h,0);
        return 0;
    }
    if(universe_hash(c,id,depth+1))return -1;
    size_t u=universe_index(c->out,id);
    if(u==c->out->universe_count)return -1;
    alea_digest_u64(g,2);
    alea_digest_u64(h,2);
    alea_sha256_update(g,&c->out->universes[u].geometry,32);
    alea_sha256_update(h,&c->out->universes[u].content,32);
    return 0;
}
static int cell_hash(context* c,size_t i,unsigned depth) {
    if(depth>128)return -1;
    if(c->cell_state[i]==2)return 0;
    if(c->cell_state[i])return -1;
    c->cell_state[i]=1;
    const alea_cell_entry_t* cell=&c->sys->cells.data[i];
    if(node_hash(c,cell->root_node_id,0))return -1;
    alea_sha256_t g,h;
    begin(c,&g,50);
    begin(c,&h,51);
    alea_sha256_update(&g,&c->nodes[cell->root_node_id],32);
    alea_sha256_update(&h,&c->nodes[cell->root_node_id],32);
    const alea_transform_t* tr=cell->fill_transform>0?alea_get_transform(c->sys,cell->fill_transform):NULL;
    if(cell->fill_transform>0 && !tr)return -1;
    double identity[12]= {
        0,0,0,1,0,0,0,1,0,0,0,1
    };
    for(int k=0;k<12;k++) {
        double v=tr?tr->cosines[k]:identity[k];
        if(alea_digest_double(&g,v)||alea_digest_double(&h,v))return -1;
    }
#define BOTH(v) do{alea_digest_u64(&g,(uint64_t)(int64_t)(v));alea_digest_u64(&h,(uint64_t)(int64_t)(v));}while(0)
    BOTH(cell->lat_type);
    BOTH(cell->lat_fill_repeating);
    BOTH(cell->lat_fill_zero_element_coords);
    if(cell->lat_type) {
        for(int k=0;k<6;k++)BOTH(cell->lat_fill_dims[k]);
        for(int k=0;k<3;k++) {
            if(alea_digest_double(&g,cell->lat_pitch[k])||alea_digest_double(&h,cell->lat_pitch[k])||alea_digest_double(&g,cell->lat_lower_left[k])||alea_digest_double(&h,cell->lat_lower_left[k]))return -1;
        }
        BOTH(cell->lat_fill_count);
        if(cell->lat_fill_count && !cell->lat_fill)return -1;
        for(size_t k=0;k<cell->lat_fill_count;k++)if(target(c,&g,&h,cell->lat_fill[k],cell->universe_id,depth))return -1;
        if(target(c,&g,&h,cell->lat_outer_universe,cell->universe_id,depth))return -1;
    }
    else if(target(c,&g,&h,cell->fill_universe,-1,depth))return -1;
#undef BOTH
    if(alea_digest_double(&h,cell->density)||alea_digest_double(&h,cell->temperature))return -1;
    alea_digest_u64(&h,cell->is_mass_density);
    alea_digest_u64(&h,cell->has_temperature);
    alea_digest_t material= {
        {
            0
        }
    };
    if(cell->material_index>=0) {
        if((size_t)cell->material_index>=c->sys->materials.count || cached_material(c,(size_t)cell->material_index,&material))return -1;
    }
    else if(cell->material_id!=0 && mixture_hash(c,cell->material_id,&material))return -1;
    alea_sha256_update(&h,&material,32);
    alea_sha256_final(&g,c->out->cells[i].geometry.bytes);
    alea_sha256_final(&h,c->out->cells[i].content.bytes);
    c->cell_state[i]=2;
    return 0;
}
static int universe_hash(context* c,int id,unsigned depth) {
    if(depth>128 || alea_interrupted())return -1;
    size_t u=universe_index(c->out,id);
    if(u==c->out->universe_count)return -1;
    if(c->universe_state[u]==2)return 0;
    if(c->universe_state[u])return -1;
    c->universe_state[u]=1;
    size_t n=0;
    for(size_t i=c->first[u];i!=SIZE_MAX;i=c->next[i])n++;
    alea_digest_t* parts=calloc(n?n:1,sizeof(*parts));
    if(!parts)return -1;
    alea_sha256_t g,h;
    begin(c,&g,60);
    begin(c,&h,61);
    for(size_t i=c->first[u];i!=SIZE_MAX;i=c->next[i])if(cell_hash(c,i,depth+1)) {
        free(parts);
        return -1;
    }
    for(int content=0;content<2;content++) {
        size_t k=0;
        for(size_t i=c->first[u];i!=SIZE_MAX;i=c->next[i])parts[k++]=content?c->out->cells[i].content:c->out->cells[i].geometry;
        qsort(parts,n,sizeof(*parts),cmp_digest);
        alea_sha256_t* hash=content?&h:&g;
        alea_digest_u64(hash,n);
        alea_sha256_update(hash,parts,n*sizeof(*parts));
    }
    free(parts);
    alea_sha256_final(&g,c->out->universes[u].geometry.bytes);
    alea_sha256_final(&h,c->out->universes[u].content.bytes);
    c->universe_state[u]=2;
    return 0;
}
static alea_fingerprints_t* fingerprints_profile(const alea_system_t* sys,alea_fingerprint_profile_t profile) {
    if(!sys)return NULL;
    context c= {
        0
    };
    c.sys=sys;
    c.profile=profile;
    c.out=calloc(1,sizeof(*c.out));
    if(!c.out)return NULL;
    c.out->version=profile==ALEA_FINGERPRINT_LEGACY?ALEA_FINGERPRINT_VERSION:3;
    c.out->cell_count=sys->cells.count;
    c.out->cells=calloc(sys->cells.count+1,sizeof(*c.out->cells));
    c.out->universes=calloc(sys->cells.count+1,sizeof(*c.out->universes));
    c.nodes=calloc(sys->nodes.count+1,sizeof(*c.nodes));
    c.node_state=calloc(sys->nodes.count+1,1);
    c.cell_state=calloc(sys->cells.count+1,1);
    c.universe_state=calloc(sys->cells.count+1,1);
    c.first=malloc((sys->cells.count+1)*sizeof(size_t));
    c.next=malloc((sys->cells.count+1)*sizeof(size_t));
    c.materials=calloc(sys->materials.count+1,sizeof(*c.materials));
    c.material_ready=calloc(sys->materials.count+1,1);
    c.mixtures=calloc(sys->mixtures.count+1,sizeof(*c.mixtures));
    c.mixture_ready=calloc(sys->mixtures.count+1,1);
    int ok=c.materials && c.material_ready && c.mixtures && c.mixture_ready && c.first && c.next && c.out->cells && c.out->universes && c.nodes && c.node_state && c.cell_state && c.universe_state;
    if(ok) {
        for(size_t i=0;i<sys->cells.count;i++) {
            c.out->cells[i].id=sys->cells.data[i].mc_cell_id;
            c.out->universes[i].id=sys->cells.data[i].universe_id;
        }
        qsort(c.out->universes,sys->cells.count,sizeof(*c.out->universes),cmp_fp);
        for(size_t i=0;i<sys->cells.count;i++)if(!i || c.out->universes[i].id!=c.out->universes[i-1].id)c.out->universes[c.out->universe_count++].id=c.out->universes[i].id;
        for(size_t u=0;u<c.out->universe_count;u++)c.first[u]=SIZE_MAX;
        for(size_t i=0;i<sys->cells.count;i++) {
            size_t u=universe_index(c.out,sys->cells.data[i].universe_id);
            c.next[i]=c.first[u];
            c.first[u]=i;
        }
        for(size_t u=0;u<c.out->universe_count;u++)if(universe_hash(&c,c.out->universes[u].id,0)) {
            ok=0;
            break;
        }
    }
    free(c.materials);
    free(c.mixtures);
    free(c.material_ready);
    free(c.mixture_ready);
    free(c.first);
    free(c.next);
    free(c.nodes);
    free(c.node_state);
    free(c.cell_state);
    free(c.universe_state);
    if(!ok) {
        alea_fingerprints_destroy(c.out);
        alea_set_error_detail(ALEA_ERR_INVALID_STATE,"Cannot fingerprint invalid, cyclic, excessive or interrupted geometry");
        return NULL;
    }
    return c.out;
}
alea_fingerprints_t* alea_fingerprints(const alea_system_t* sys) {
    return fingerprints_profile(sys,ALEA_FINGERPRINT_LEGACY);
}
void alea_fingerprints_destroy(alea_fingerprints_t* f) {
    if(f) {
        free(f->cells);
        free(f->universes);
        free(f);
    }
}
void alea_digest_hex(const alea_digest_t* d,char hex[65]) {
    static const char chars[]="0123456789abcdef";
    for(int i=0;i<32;i++) {
        hex[2*i]=chars[d->bytes[i]>>4];
        hex[2*i+1]=chars[d->bytes[i]&15];
    }
    hex[64]=0;
}
static int cmp_content(const void* a,const void* b) {
    return memcmp(&((const alea_fingerprint_t*)a)->content,&((const alea_fingerprint_t*)b)->content,32);
}
static int diff_entries(alea_model_diff_t* d,const alea_fingerprint_t* a,size_t na,const alea_fingerprint_t* b,size_t nb,int universe) {
    alea_fingerprint_t* left=malloc((na+1)*sizeof(*left));
    alea_fingerprint_t* right=malloc((nb+1)*sizeof(*right));
    if(!left || !right) {
        free(left);
        free(right);
        return -1;
    }
    if(na)memcpy(left,a,na*sizeof(*a));
    if(nb)memcpy(right,b,nb*sizeof(*b));
    qsort(left,na,sizeof(*left),cmp_fp);
    qsort(right,nb,sizeof(*right),cmp_fp);
    size_t i=0,j=0,nl=0,nr=0;
    while(i<na || j<nb) {
        if(j==nb || (i<na && left[i].id<right[j].id)) {
            left[nl++]=left[i++];
            continue;
        }
        if(i==na || right[j].id<left[i].id) {
            right[nr++]=right[j++];
            continue;
        }
        alea_diff_kind_t kind=memcmp(&left[i].geometry,&right[j].geometry,32)?ALEA_DIFF_GEOMETRY:memcmp(&left[i].content,&right[j].content,32)?ALEA_DIFF_CONTENT:ALEA_DIFF_UNCHANGED;
        d->changes[d->count++]=(alea_model_change_t) {
            universe,left[i].id,right[j].id,kind
        };
        i++;
        j++;
    }
    qsort(left,nl,sizeof(*left),cmp_content);
    qsort(right,nr,sizeof(*right),cmp_content);
    i=j=0;
    while(i<nl || j<nr) {
        int cmp=i==nl?1:j==nr?-1:cmp_content(&left[i],&right[j]);
        if(cmp<0) {
            d->changes[d->count++]=(alea_model_change_t) {
                universe,left[i++].id,-1,ALEA_DIFF_REMOVED
            };
            continue;
        }
        if(cmp>0) {
            d->changes[d->count++]=(alea_model_change_t) {
                universe,-1,right[j++].id,ALEA_DIFF_ADDED
            };
            continue;
        }
        size_t ie=i+1,je=j+1;
        while(ie<nl && !cmp_content(&left[i],&left[ie]))ie++;
        while(je<nr && !cmp_content(&right[j],&right[je]))je++;
        if(ie==i+1 && je==j+1)d->changes[d->count++]=(alea_model_change_t) {
            universe,left[i].id,right[j].id,ALEA_DIFF_RENUMBERED
        };
        else {
            for(size_t k=i;k<ie;k++)d->changes[d->count++]=(alea_model_change_t) {
                universe,left[k].id,-1,ALEA_DIFF_AMBIGUOUS
            };
            for(size_t k=j;k<je;k++)d->changes[d->count++]=(alea_model_change_t) {
                universe,-1,right[k].id,ALEA_DIFF_ADDED
            };
        }
        i=ie;
        j=je;
    }
    free(left);
    free(right);
    return 0;
}
static alea_model_diff_t* fingerprint_diff(const alea_fingerprints_t* a,const alea_fingerprints_t* b) {
    if(!a || !b || b->version!=a->version)return NULL;
    size_t n=a->cell_count+b->cell_count+a->universe_count+b->universe_count;
    if(n>SIZE_MAX/sizeof(alea_model_change_t))return NULL;
    alea_model_diff_t* d=calloc(1,sizeof(*d));
    if(!d)return NULL;
    d->changes=calloc(n?n:1,sizeof(*d->changes));
    if(!d->changes || diff_entries(d,a->cells,a->cell_count,b->cells,b->cell_count,0)||diff_entries(d,a->universes,a->universe_count,b->universes,b->universe_count,1)) {
        alea_model_diff_destroy(d);
        return NULL;
    }
    return d;
}
void alea_model_diff_destroy(alea_model_diff_t* d) {
    if(d) {
        free(d->changes);
        free(d);
    }
}

alea_model_diff_t* alea_fingerprint_diff(const alea_fingerprints_t* a,const alea_fingerprints_t* b) {
    if(!a || a->version!=ALEA_FINGERPRINT_VERSION)return NULL;
    return fingerprint_diff(a,b);
}
struct alea_fingerprint_set {
    alea_fingerprint_profile_t profile;
    alea_fingerprints_t* values;
};
alea_fingerprint_set_t* alea_fingerprint_set_create(const alea_system_t* sys,
    const alea_fingerprint_options_t* options) {
    if(!options || options->struct_size!=sizeof(*options) ||
       (options->profile!=ALEA_FINGERPRINT_LEGACY && options->profile!=ALEA_FINGERPRINT_SEMANTIC)) {
        alea_set_error_detail(ALEA_ERR_INVALID_ARG,"Invalid fingerprint options");
        return NULL;
    }
    alea_fingerprint_set_t* set=calloc(1,sizeof(*set));
    if(!set)return NULL;
    set->profile=options->profile;
    set->values=fingerprints_profile(sys,options->profile);
    if(!set->values){free(set);return NULL;}
    return set;
}
const alea_fingerprints_t* alea_fingerprint_set_values(const alea_fingerprint_set_t* set) {
    return set?set->values:NULL;
}
alea_fingerprint_profile_t alea_fingerprint_set_profile(const alea_fingerprint_set_t* set) {
    return set?set->profile:ALEA_FINGERPRINT_LEGACY;
}
void alea_fingerprint_set_destroy(alea_fingerprint_set_t* set) {
    if(set){alea_fingerprints_destroy(set->values);free(set);}
}
alea_model_diff_t* alea_fingerprint_set_diff(const alea_fingerprint_set_t* a,
    const alea_fingerprint_set_t* b) {
    if(!a || !b || a->profile!=b->profile){
        alea_set_error_detail(ALEA_ERR_INVALID_ARG,"Incompatible fingerprint profiles");
        return NULL;
    }
    return fingerprint_diff(a->values,b->values);
}
