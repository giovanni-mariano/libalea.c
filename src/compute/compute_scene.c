// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
// SPDX-License-Identifier: MPL-2.0
#include "compute_internal.h"
#include "core/alea_system.h"
#include <stdlib.h>
#include <string.h>

_Static_assert(sizeof(metal_point_cell)==12,"Metal point links layout");
_Static_assert(sizeof(compute_node)==8,"GPU instruction layout");
_Static_assert(sizeof(compute_primitive)==56,"GPU primitive layout");
_Static_assert(sizeof(compute_cell)==36,"GPU cell layout");
_Static_assert(sizeof(compute_result)==sizeof(alea_compute_result_t),"GPU result layout");

/* Translate quadrics near their stationary point in double precision before
 * rounding coefficients. A small cylinder far from the origin otherwise loses
 * its radius in the constant coefficient. Pivoted elimination also handles
 * rank-deficient quadrics; unconstrained coordinates remain zero. */
static void shift_quadric(double p[13]) {
    double a[3][4]={{p[0],p[3]*.5,p[5]*.5,-p[6]*.5},
                    {p[3]*.5,p[1],p[4]*.5,-p[7]*.5},
                    {p[5]*.5,p[4]*.5,p[2],-p[8]*.5}};
    int cols[3]={0,1,2}; double scale=0;
    for(int i=0;i<3;++i)for(int j=0;j<3;++j)scale=fmax(scale,fabs(a[i][j]));
    int rank=0;
    for(int k=0;k<3;++k) {
        int row=k,col=k;
        for(int i=k;i<3;++i)for(int j=k;j<3;++j)if(fabs(a[i][j])>fabs(a[row][col])){row=i;col=j;}
        /* Near-null directions of rounded MCNP cylinders can have a formal
         * stationary point hundreds of millions of units away. Translating
         * there destroys local detail in FP32. Stop solving those directions;
         * all original quadratic terms and residual linear terms are retained,
         * so this changes the coordinate origin, not the represented surface. */
        if(fabs(a[row][col])<=scale*(32.0*FLT_EPSILON))break;
        for(int j=0;j<4;++j){double t=a[k][j];a[k][j]=a[row][j];a[row][j]=t;}
        for(int i=0;i<3;++i){double t=a[i][k];a[i][k]=a[i][col];a[i][col]=t;}
        int c=cols[k];cols[k]=cols[col];cols[col]=c;
        double pivot=a[k][k];for(int j=k;j<4;++j)a[k][j]/=pivot;
        for(int i=0;i<3;++i)if(i!=k){double f=a[i][k];for(int j=k;j<4;++j)a[i][j]-=f*a[k][j];}
        rank=k+1;
    }
    double o[3]={0}; for(int i=0;i<rank;++i)o[cols[i]]=a[i][3];
    double x=o[0],y=o[1],z=o[2];
    p[9]+=x*(p[0]*x+p[3]*y+p[5]*z+p[6])+y*(p[1]*y+p[4]*z+p[7])+z*(p[2]*z+p[8]);
    p[6]+=2*p[0]*x+p[3]*y+p[5]*z;
    p[7]+=2*p[1]*y+p[3]*x+p[4]*z;
    p[8]+=2*p[2]*z+p[5]*x+p[4]*y;
    double norm=0;for(int i=0;i<10;++i)norm=fmax(norm,fabs(p[i]));
    if(norm>0)for(int i=0;i<10;++i)p[i]/=norm;
    p[10]=x;p[11]=y;p[12]=z;
}

static int pack_primitive(const alea_system_t* sys, uint32_t id,
                          compute_primitive* out) {
    const void* data = id < sys->primitives.count ? alea_primitive_payload_const(sys, id) : NULL;
    if (!data) {
        alea_set_error_detail(ALEA_ERR_INVALID_ID, "Invalid compute primitive %u", id);
        return -1;
    }
    double p[13] = {0};
    switch (sys->primitives.data[id].type) {
    case ALEA_PRIMITIVE_PLANE: {
        const alea_plane_data_t* v = data;
        out->type = COMPUTE_PLANE;
        p[0] = v->a; p[1] = v->b; p[2] = v->c; p[3] = v->d;
        break;
    }
    case ALEA_PRIMITIVE_SPHERE: {
        const alea_sphere_data_t* v = data;
        out->type = COMPUTE_SPHERE;
        p[0] = v->center_x; p[1] = v->center_y; p[2] = v->center_z; p[3] = v->radius;
        break;
    }
    case ALEA_PRIMITIVE_CYLINDER_X: {
        const alea_cylinder_x_data_t* v = data;
        out->type = COMPUTE_CYL_X;
        p[0] = v->center_y; p[1] = v->center_z; p[2] = v->radius;
        break;
    }
    case ALEA_PRIMITIVE_CYLINDER_Y: {
        const alea_cylinder_y_data_t* v = data;
        out->type = COMPUTE_CYL_Y;
        p[0] = v->center_x; p[1] = v->center_z; p[2] = v->radius;
        break;
    }
    case ALEA_PRIMITIVE_CYLINDER_Z: {
        const alea_cylinder_z_data_t* v = data;
        out->type = COMPUTE_CYL_Z;
        p[0] = v->center_x; p[1] = v->center_y; p[2] = v->radius;
        break;
    }
    case ALEA_PRIMITIVE_RPP: {
        const alea_box_data_t* v = data;
        out->type = COMPUTE_BOX;
        p[0] = v->min_x; p[1] = v->max_x; p[2] = v->min_y;
        p[3] = v->max_y; p[4] = v->min_z; p[5] = v->max_z;
        break;
    }
    case ALEA_PRIMITIVE_QUADRIC: {
        const alea_quadric_data_t* v=data;
        out->type=COMPUTE_QUADRIC;
        memcpy(p,v->coeffs,10*sizeof(double));
        shift_quadric(p);
        break;
    }
    case ALEA_PRIMITIVE_CONE_X:
    case ALEA_PRIMITIVE_CONE_Y:
    case ALEA_PRIMITIVE_CONE_Z: {
        /* The three public cone structures have the same scalar layout. */
        double apex[4]; int sheet;
        if (sys->primitives.data[id].type==ALEA_PRIMITIVE_CONE_X) {
            const alea_cone_x_data_t* v=data;
            apex[0]=v->apex_x; apex[1]=v->apex_y; apex[2]=v->apex_z; apex[3]=v->tan_angle_sq; sheet=v->sheet_selection;
        } else if (sys->primitives.data[id].type==ALEA_PRIMITIVE_CONE_Y) {
            const alea_cone_y_data_t* v=data;
            apex[0]=v->apex_x; apex[1]=v->apex_y; apex[2]=v->apex_z; apex[3]=v->tan_angle_sq; sheet=v->sheet_selection;
        } else {
            const alea_cone_z_data_t* v=data;
            apex[0]=v->apex_x; apex[1]=v->apex_y; apex[2]=v->apex_z; apex[3]=v->tan_angle_sq; sheet=v->sheet_selection;
        }
        out->type=COMPUTE_CONE_X+sys->primitives.data[id].type-ALEA_PRIMITIVE_CONE_X;
        memcpy(p,apex,sizeof(apex)); p[4]=sheet;
        break;
    }
    case ALEA_PRIMITIVE_TORUS_X:
    case ALEA_PRIMITIVE_TORUS_Y:
    case ALEA_PRIMITIVE_TORUS_Z: {
        const alea_torus_data_t* v=data;
        out->type=COMPUTE_TORUS_X+sys->primitives.data[id].type-ALEA_PRIMITIVE_TORUS_X;
        p[0]=v->center_x; p[1]=v->center_y; p[2]=v->center_z;
        p[3]=v->major_radius; p[4]=v->minor_radius; p[5]=v->axial_semiwidth_B;
        if (p[3]<0 || p[4]<=0 || p[5]<=0) {
            alea_set_error_detail(ALEA_ERR_UNSUPPORTED,"Metal requires positive torus tube radii (primitive %u)",id); return -1;
        }
        break;
    }
    case ALEA_PRIMITIVE_RCC: {
        const alea_rcc_data_t* v=data;
        double h=hypot(hypot(v->height_x,v->height_y),v->height_z);
        if (h<=0 || v->radius<=0) {
            alea_set_error_detail(ALEA_ERR_INVALID_STATE,"Degenerate RCC primitive %u",id); return -1;
        }
        out->type=COMPUTE_RCC;
        p[0]=v->base_x; p[1]=v->base_y; p[2]=v->base_z;
        p[3]=v->height_x/h; p[4]=v->height_y/h; p[5]=v->height_z/h;
        p[6]=h; p[7]=v->radius;
        break;
    }
    default:
        alea_set_error_detail(ALEA_ERR_UNSUPPORTED,
            "Metal does not yet support primitive type %d", sys->primitives.data[id].type);
        return -1;
    }
    for (int i = 0; i < 13; ++i) {
        if (!isfinite(p[i]) || fabs(p[i]) > FLT_MAX) {
            alea_set_error_detail(ALEA_ERR_UNSUPPORTED, "Primitive %u does not fit FP32", id);
            return -1;
        }
        out->p[i] = (float)p[i];
    }
    return 0;
}


typedef struct { uint32_t count; uint16_t height; uint8_t need,state; } node_info;
typedef struct {
    const alea_system_t* sys;
    compute_scene* scene;
    node_info* info;
    uint32_t* primitive_map;
} pack_context;

#define COMPUTE_SOURCE_DEPTH 4096
static int inspect_node(pack_context* ctx,uint32_t id,unsigned depth) {
    if(id>=ctx->sys->nodes.count) {
        alea_set_error_detail(ALEA_ERR_INVALID_ID,"Invalid compute node %u",id); return -1;
    }
    node_info* m=&ctx->info[id];
    if(depth>=COMPUTE_SOURCE_DEPTH || depth+m->height>COMPUTE_SOURCE_DEPTH) {
        alea_set_error_detail(ALEA_ERR_UNSUPPORTED,"Metal source CSG depth exceeds %d",COMPUTE_SOURCE_DEPTH); return -1;
    }
    if((m->state&3)==2)return 0;
    if((m->state&3)==1) { alea_set_error_detail(ALEA_ERR_INVALID_STATE,"Cycle in compute CSG"); return -1; }
    m->state=1;
    const alea_node_t* n=&ctx->sys->nodes.data[id];
    alea_operation_t op=ALEA_GET_OPERATION(n);
    m->height=m->need=m->count=1;
    if(op==ALEA_OP_PRIMITIVE) {
        uint32_t p=n->primitive.primitive_id;
        if(p>=ctx->sys->primitives.count) {alea_set_error_detail(ALEA_ERR_INVALID_ID,"Invalid primitive %u",p);return -1;}
        if(ctx->primitive_map[p]==COMPUTE_NONE) {
            uint32_t at=(uint32_t)ctx->scene->primitive_count;
            if(pack_primitive(ctx->sys,p,&ctx->scene->primitives[at]))return -1;
            ctx->primitive_map[p]=at; ++ctx->scene->primitive_count;
        }
    } else {
        if(op!=ALEA_OP_UNION && op!=ALEA_OP_INTERSECTION && op!=ALEA_OP_DIFFERENCE && op!=ALEA_OP_COMPLEMENT) {
            alea_set_error_detail(ALEA_ERR_UNSUPPORTED,"Unsupported compute CSG operation");return -1;
        }
        if(inspect_node(ctx,n->operation.left,depth+1))return -1;
        const node_info* a=&ctx->info[n->operation.left];
        uint64_t count=(uint64_t)a->count+1;
        m->height=a->height+1;m->need=a->need;
        if(op!=ALEA_OP_COMPLEMENT) {
            if(inspect_node(ctx,n->operation.right,depth+1))return -1;
            const node_info* b=&ctx->info[n->operation.right];
            count+=(uint64_t)b->count+1; /* conditional jump before the second operand */
            if(b->height>=m->height)m->height=b->height+1;
            unsigned ordered=a->need>b->need+1?a->need:b->need+1;
            if(ordered>COMPUTE_STACK && b->need>a->need) {
                m->state|=4; m->need=b->need>a->need+1?b->need:a->need+1;
            } else m->need=ordered;
        }
        if(count>UINT32_MAX || m->need>COMPUTE_STACK) {
            alea_set_error_detail(ALEA_ERR_OVERFLOW,"CSG program exceeds device limits");return -1;
        }
        m->count=(uint32_t)count;
    }
    m->state=(m->state&4)|2;return 0;
}

static void emit_node(pack_context* ctx,uint32_t id,size_t* at) {
    const alea_node_t* n=&ctx->sys->nodes.data[id];
    alea_operation_t op=ALEA_GET_OPERATION(n);
    compute_node instruction={0};
    if(op==ALEA_OP_PRIMITIVE) {
        uint32_t p=ctx->primitive_map[n->primitive.primitive_id];
        unsigned flip=(n->primitive.sense>0)!=(n->primitive.inverted!=0);
        instruction.op=COMPUTE_LEAF;instruction.arg=(p<<1)|flip;
    } else {
        uint32_t left=n->operation.left,right=n->operation.right;
        int swap=(ctx->info[id].state&4)!=0;
        emit_node(ctx,swap?right:left,at);
        if(op!=ALEA_OP_COMPLEMENT) {
            uint32_t second=swap?left:right;
            compute_node jump={op==ALEA_OP_UNION?COMPUTE_SKIP_UNION:
                (op==ALEA_OP_DIFFERENCE && swap)?COMPUTE_SKIP_REVERSE_DIFFERENCE:COMPUTE_SKIP_INTERSECTION,
                ctx->info[second].count+1};
            ctx->scene->nodes[(*at)++]=jump;
            emit_node(ctx,second,at);
        }
        instruction.op=op==ALEA_OP_COMPLEMENT?COMPUTE_COMPLEMENT:
            op==ALEA_OP_UNION?COMPUTE_UNION:op==ALEA_OP_INTERSECTION?COMPUTE_INTERSECTION:
            swap?COMPUTE_REVERSE_DIFFERENCE:COMPUTE_DIFFERENCE;
    }
    ctx->scene->nodes[(*at)++]=instruction;
}

void alea_compute_scene_free(compute_scene* s) {
    free(s->nodes);free(s->primitives);free(s->cells);free(s->spatial);free(s->point_cells);
    free(s->point_lattices);free(s->point_fills);
    alea_compute_render_free(&s->render);
    memset(s,0,sizeof(*s));
}

int alea_compute_pack_geometry(const alea_system_t* sys,compute_scene* s,int accelerate) {
    memset(s,0,sizeof(*s));
    if(alea_compute_validate_rect_lattices(sys))return -1;
    size_t nn=sys->nodes.count,nc=sys->cells.count,np=sys->primitives.count;
    if(nn>UINT32_MAX || nc>=UINT32_MAX || np>UINT32_MAX/2 ||
       nn>SIZE_MAX/sizeof(node_info) || np>SIZE_MAX/sizeof(compute_primitive) || nc>SIZE_MAX/sizeof(compute_cell)) {
        alea_set_error_detail(ALEA_ERR_OVERFLOW,"Compute geometry exceeds index limits");return -1;
    }
    node_info* info=calloc(nn?nn:1,sizeof(*info));
    uint32_t* map=malloc((np?np:1)*sizeof(*map));
    s->primitives=calloc(np?np:1,sizeof(*s->primitives));
    s->cells=calloc(nc?nc:1,sizeof(*s->cells));
    if(!info || !map || !s->primitives || !s->cells)goto memory;
    for(size_t i=0;i<np;++i)map[i]=COMPUTE_NONE;
    pack_context ctx={sys,s,info,map};
    size_t total=0;
    for(size_t i=0;i<nc;++i) {
        if(alea_interrupted()){alea_set_error_detail(ALEA_ERR_INTERRUPTED,"Compute packing interrupted");goto fail;}
        uint32_t root=sys->cells.data[i].root_node_id;
        if(inspect_node(&ctx,root,0))goto fail;
        if(info[root].count>UINT32_MAX-total) {alea_set_error_detail(ALEA_ERR_OVERFLOW,"Compute instruction stream exceeds index limit");goto fail;}
        total+=info[root].count;
    }
    if(total>SIZE_MAX/sizeof(compute_node))goto memory;
    s->nodes=calloc(total?total:1,sizeof(*s->nodes));
    if(!s->nodes)goto memory;
    size_t at=0;
    for(size_t i=0;i<nc;++i) {
        const alea_cell_entry_t* c=&sys->cells.data[i];
        s->cells[i]=(compute_cell){(uint32_t)at,info[c->root_node_id].count,
            c->mc_cell_id,c->material_id,COMPUTE_NONE,COMPUTE_NONE,COMPUTE_NONE,COMPUTE_NONE,0};
        emit_node(&ctx,c->root_node_id,&at);
    }
    s->node_count=total;s->cell_count=nc;
    free(info);free(map);info=NULL;map=NULL;
    if(alea_compute_render_pack(sys,s))goto fail;
    if(accelerate && alea_compute_spatial_pack(s))goto fail;
    if(alea_compute_point_programs_pack(s))goto fail;
    if(accelerate && alea_compute_bvh_rebuild(s))goto fail;
    return 0;
memory:
    alea_set_error_detail(ALEA_ERR_OUT_OF_MEMORY,"Packing Metal geometry");
fail:
    free(info);free(map);alea_compute_scene_free(s);return -1;
}

int alea_compute_pack(const alea_system_t* sys,compute_scene* s) {
    return alea_compute_pack_geometry(sys,s,1);
}

/* Exact encoding; all arguments retain their bits. A wider argument or future
 * opcode keeps the existing format instead of truncating or rejecting input. */
int alea_compute_code_compact(const compute_node* code,size_t count,uint32_t* output) {
    for(size_t i=0;i<count;++i) {
        if((i&65535)==0 && alea_interrupted()) {
            alea_set_error_detail(ALEA_ERR_INTERRUPTED,"Code encoding interrupted");return -1;
        }
        if(code[i].op>15 || code[i].arg>0x0fffffffu)return 0;
    }
    for(size_t i=0;i<count;++i) {
        if((i&65535)==0 && alea_interrupted()) {
            alea_set_error_detail(ALEA_ERR_INTERRUPTED,"Code encoding interrupted");return -1;
        }
        compute_node node=code[i];
        output[i]=(node.arg<<4)|node.op;
    }
    return 1;
}

/* Called only after all structured-program consumers have finished. */
typedef struct { uint32_t root,count; } jump_range;
static int jump_range_order(const void* a,const void* b) {
    const jump_range* x=a,*y=b;
    if(x->root!=y->root)return x->root<y->root?-1:1;
    return (x->count>y->count)-(x->count<y->count);
}
int alea_compute_thread_jumps(compute_scene* s,size_t* rewritten) {
    *rewritten=0;
    size_t capacity=s->cell_count+s->spatial_count;
    if(capacity<s->cell_count || capacity>64*1024*1024/sizeof(jump_range)) {
        alea_set_error_detail(ALEA_ERR_OVERFLOW,"Jump range scratch exceeds 64 MiB");return -1;
    }
    jump_range* ranges=malloc((capacity?capacity:1)*sizeof(*ranges));
    if(!ranges){alea_set_error_detail(ALEA_ERR_OUT_OF_MEMORY,"Allocating jump ranges");return -1;}
    size_t n=0;
    for(size_t i=0;i<s->cell_count;++i)ranges[n++]=(jump_range){s->cells[i].root,s->cells[i].count};
    for(size_t i=0;i<s->spatial_count;++i)if(s->spatial[i].child==COMPUTE_NONE)
        ranges[n++]=(jump_range){s->spatial[i].root,s->spatial[i].count};
    qsort(ranges,n,sizeof(*ranges),jump_range_order);
    size_t unique=0;
    for(size_t i=0;i<n;++i) {
        jump_range r=ranges[i];
        if((uint64_t)r.root+r.count>s->node_count)goto invalid;
        if(unique && r.root==ranges[unique-1].root && r.count==ranges[unique-1].count)continue;
        ranges[unique++]=r;
    }
    /* Validate every branch before changing any instruction. */
    for(size_t i=0;i<unique;++i)for(uint32_t pc=0;pc<ranges[i].count;++pc) {
        if((pc&65535)==0 && alea_interrupted())goto interrupted;
        compute_node node=s->nodes[(size_t)ranges[i].root+pc];
        if(node.op==COMPUTE_SKIP_UNION || node.op==COMPUTE_SKIP_INTERSECTION || node.op==COMPUTE_SKIP_REVERSE_DIFFERENCE)
            if((uint64_t)pc+node.arg+1>ranges[i].count)goto invalid;
    }
    for(size_t first=0;first<unique;) {
        size_t last=first+1;uint64_t end=(uint64_t)ranges[first].root+ranges[first].count;
        while(last<unique && ranges[last].root<end) {
            uint64_t next=(uint64_t)ranges[last].root+ranges[last].count;if(next>end)end=next;++last;
        }
        /* Nonidentical overlapping programs retain their original branches. */
        if(last==first+1) {
            jump_range r=ranges[first];
            for(uint32_t pc=r.count;pc-->0;) {
                if((pc&65535)==0 && alea_interrupted())goto interrupted;
                compute_node* node=&s->nodes[(size_t)r.root+pc];
                if(node->op!=COMPUTE_SKIP_UNION && node->op!=COMPUTE_SKIP_INTERSECTION)continue;
                uint64_t target=(uint64_t)pc+node->arg+1;
                if(target>=r.count)continue;
                compute_node next=s->nodes[(size_t)r.root+target];
                uint64_t arg=(uint64_t)node->arg+next.arg+1;
                if(next.op==node->op && arg<=0x0fffffffu){node->arg=(uint32_t)arg;++*rewritten;}
            }
        }
        first=last;
    }
    free(ranges);return 0;
interrupted:
    free(ranges);alea_set_error_detail(ALEA_ERR_INTERRUPTED,"Jump threading interrupted; discard snapshot");return -1;
invalid:
    free(ranges);alea_set_error_detail(ALEA_ERR_INVALID_STATE,"Invalid jump program range");return -1;
}
