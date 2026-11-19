// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
// SPDX-License-Identifier: MPL-2.0
#include "alea_snapshot.h"
#include "snapshot_internal.h"
#include "primitives/primitive_codec.h"
#include "alea_engine_hash.h"
#include "util/compat.h"
#include "util/load_profile.h"
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <math.h>
const char* alea_snapshot_engine_id(void){return ALEA_ENGINE_HASH;}
int alea_snapshot_read_limits(const alea_snapshot_read_options_t* options,
                              size_t* file_bytes,size_t* decoded_bytes) {
    if(options && options->struct_size<sizeof(*options)) {
        alea_set_error_detail(ALEA_ERR_INVALID_ARG,"Invalid snapshot read options size");
        return -1;
    }
    *file_bytes=options?options->max_file_bytes:0;
    *decoded_bytes=options?options->max_decoded_bytes:0;
    if(!*decoded_bytes) {
        uint64_t default_limit=UINT64_C(8)<<30;
        *decoded_bytes=default_limit>SIZE_MAX?SIZE_MAX:(size_t)default_limit;
    }
    return 0;
}
void alea_snapshot_profile_allocations(const alea_alloc_budget_t* budget) {
    if(alea_load_profile_start()>=0)
        fprintf(stderr,"[alea-snapshot-allocation] requested=%zu limit=%zu exceeded=%d\n",
                budget->requested,budget->limit,budget->exceeded);
}
static void flush_snapshot(alea_snapshot_io* io){
    if(io->failed || !io->buffered)return;
    if(fwrite(io->buffer,1,io->buffered,io->file)!=io->buffered){io->failed=1;return;}
    alea_sha256_update(&io->digest,io->buffer,io->buffered);io->buffered=0;
}
void alea_snapshot_bytes(alea_snapshot_io* io,void* data,size_t n){
    if(io->failed)return;
    if(!n)return;
    if(n>io->remaining){io->failed=1;return;}
    if(io->mapped) {
        memcpy(data,io->mapping.data+io->position,n);
        io->position+=n;io->remaining-=n;return;
    }
    unsigned char* bytes=data;
    while(n){
        size_t take;
        if(io->read){
            if(io->position==io->buffered){
                size_t amount=io->remaining<sizeof(io->buffer)?(size_t)io->remaining:sizeof(io->buffer);
                if(fread(io->buffer,1,amount,io->file)!=amount){io->failed=1;return;}
                double start=io->profile?alea_monotonic_seconds():0;
                alea_sha256_update(&io->digest,io->buffer,amount);
                if(io->profile)io->decode_hash_seconds+=alea_monotonic_seconds()-start;
                io->position=0;io->buffered=amount;
            }
            take=io->buffered-io->position;if(take>n)take=n;
            memcpy(bytes,io->buffer+io->position,take);io->position+=take;
        }else{
            if(io->buffered==sizeof(io->buffer))flush_snapshot(io);
            if(io->failed)return;
            take=sizeof(io->buffer)-io->buffered;if(take>n)take=n;
            memcpy(io->buffer+io->buffered,bytes,take);io->buffered+=take;
        }
        io->remaining-=take;bytes+=take;n-=take;
    }
}
/* Consume a contiguous field without copying through a temporary byte buffer.
 * On a buffer boundary, leave the cursor untouched for the checked slow path. */
static inline const unsigned char* snapshot_read_span(alea_snapshot_io* io,size_t n) {
    if(io->failed)return NULL;
    if(n>io->remaining){io->failed=1;return NULL;}
    const unsigned char* p;
    if(io->mapped)p=(const unsigned char*)io->mapping.data+io->position;
    else {
        if(n>io->buffered-io->position)return NULL;
        p=io->buffer+io->position;
    }
    io->position+=n;io->remaining-=n;
    return p;
}

uint64_t alea_snapshot_uint(alea_snapshot_io* io,uint64_t value){
    /* Zigzag keeps signed sentinels and small positive IDs compact. */
    uint64_t wire=(value<<1)-(value>>63);
    if(value>>63)wire=~(value<<1);
    if(io->read){
        wire=0;
        for(unsigned shift=0;shift<70;shift+=7){
            const unsigned char* p=snapshot_read_span(io,1);
            unsigned char byte=0;
            if(p)byte=*p;
            else alea_snapshot_bytes(io,&byte,1);
            if(io->failed || (shift==63 && byte>1)){io->failed=1;return 0;}
            wire|=(uint64_t)(byte&127)<<shift;
            if(!(byte&128)){if(shift && byte==0)io->failed=1;return (wire>>1)^(UINT64_C(0)-(wire&1));}
        }
        io->failed=1;return 0;
    }
    do{unsigned char byte=(unsigned char)(wire&127);wire>>=7;if(wire)byte|=128;alea_snapshot_bytes(io,&byte,1);}while(wire && !io->failed);
    return value;
}
double alea_snapshot_double(alea_snapshot_io* io,double value){
    if(io->read) {
        unsigned char fallback[8]={0};
        const unsigned char* p=snapshot_read_span(io,8);
        if(!p){alea_snapshot_bytes(io,fallback,8);p=fallback;}
        uint64_t bits=0;
        for(int i=0;i<8;i++)bits|=(uint64_t)p[i]<<(8*i);
        memcpy(&value,&bits,8);
        if(isnan(value))io->failed=1;
        return value;
    }
    uint64_t bits;memcpy(&bits,&value,8);unsigned char bytes[8];
    for(int i=0;i<8;i++)bytes[i]=(unsigned char)(bits>>(8*i));
    alea_snapshot_bytes(io,bytes,8);bits=0;
    for(int i=0;i<8;i++)bits|=(uint64_t)bytes[i]<<(8*i);
    memcpy(&value,&bits,8);if(isnan(value))io->failed=1;return value;
}
float alea_snapshot_float(alea_snapshot_io* io,float value){
    if(io->read) {
        unsigned char fallback[4]={0};
        const unsigned char* p=snapshot_read_span(io,4);
        if(!p){alea_snapshot_bytes(io,fallback,4);p=fallback;}
        uint32_t bits=0;
        for(int i=0;i<4;i++)bits|=(uint32_t)p[i]<<(8*i);
        memcpy(&value,&bits,4);
        if(isnan(value))io->failed=1;
        return value;
    }
    uint32_t bits;memcpy(&bits,&value,4);unsigned char bytes[4];
    for(int i=0;i<4;i++)bytes[i]=(unsigned char)(bits>>(8*i));
    alea_snapshot_bytes(io,bytes,4);bits=0;
    for(int i=0;i<4;i++)bits|=(uint32_t)bytes[i]<<(8*i);
    memcpy(&value,&bits,4);if(isnan(value))io->failed=1;return value;
}
size_t alea_snapshot_count(alea_snapshot_io* io,size_t n,size_t size) {
    uint64_t value=alea_snapshot_uint(io,n);
    if(io->failed || value>UINT32_MAX || value>SIZE_MAX/(size?size:1) || value>io->remaining) {
        io->failed=1;
        return 0;
    }
    return (size_t)value;
}
void alea_snapshot_string(alea_snapshot_io* io,char** text) {
    uint64_t n=alea_snapshot_uint(io,*text?strlen(*text)+1:0);
    if(io->failed || n>io->remaining || n>SIZE_MAX) {
        io->failed=1;
        return;
    }
    if(io->read) {
        *text=n?alea_malloc((size_t)n):NULL;
        if(n && !*text) {
            io->failed=1;
            return;
        }
    }
    if(n) {
        alea_snapshot_bytes(io,*text,(size_t)n);
        if((*text)[n-1] || memchr(*text,0,(size_t)n-1))io->failed=1;
    }
}
static int snapshot_header(alea_snapshot_io* io,unsigned kind) {
    unsigned char magic[8]={'A','L','E','A','B','I','N',3},expected[8];
    memcpy(expected,magic,8);
    char engine[65]=ALEA_ENGINE_HASH;
    char wanted[65]=ALEA_ENGINE_HASH;
    alea_snapshot_bytes(io,magic,8);
    uint64_t actual=alea_snapshot_uint(io,kind);
    alea_snapshot_bytes(io,engine,65);
    if(memcmp(magic,expected,8)||actual!=kind||memcmp(engine,wanted,65))io->failed=1;
    return io->failed?-1:0;
}

static int mapped_checksum(alea_snapshot_io* io) {
    alea_sha256_t hash;
    unsigned char digest[32];
    alea_sha256_init(&hash);
    alea_sha256_update(&hash,io->mapping.data,io->mapping.size-32);
    alea_sha256_final(&hash,digest);
    return memcmp(digest,io->mapping.data+io->mapping.size-32,32)?-1:0;
}

int alea_snapshot_open_read(alea_snapshot_io* io,const char* path,size_t max_bytes,unsigned kind) {
    *io=(alea_snapshot_io){0};
    const char* mode=getenv("ALEA_SNAPSHOT_IO");
    if(mode && !strcmp(mode,"mapped") && alea_file_map(path,&io->mapping)) {
        io->mapped=io->read=1;
        io->failed=1;
        size_t length=io->mapping.size;
        uint64_t limit=max_bytes?max_bytes:(UINT64_C(2)<<30);
        if(length<113 || (uint64_t)length>limit)return -1;
        double start=alea_load_profile_start();
        io->profile=start>=0;
        if(mapped_checksum(io))return -1;
        alea_load_profile_end("snapshot_checksum",start);
        io->remaining=length-32;io->failed=0;
        if(start>=0)fprintf(stderr,"[alea-snapshot-io] mapped\n");
        return snapshot_header(io,kind);
    }
    /* Failed map handles are already closed by the compatibility layer. */
    FILE* file=fopen(path,"rb");
    if(!file){*io=(alea_snapshot_io){.failed=1};return -1;}
    if(alea_load_profile_start()>=0)fprintf(stderr,"[alea-snapshot-io] buffered\n");
    return alea_snapshot_begin(io,file,1,max_bytes,kind);
}

void alea_snapshot_close_read(alea_snapshot_io* io) {
    if(io->mapped)alea_file_unmap(&io->mapping);
    else if(io->file)fclose(io->file);
    io->mapped=0;io->file=NULL;
}

int alea_snapshot_begin(alea_snapshot_io* io,FILE* file,int read,size_t max_bytes,unsigned kind) {
    *io=(alea_snapshot_io){.file=file,.read=read,.failed=1};

    if(read) {
        double checksum_start=alea_load_profile_start();
        if(fseek(file,0,SEEK_END))return -1;
        long length=ftell(file);
        uint64_t limit=max_bytes?max_bytes:(UINT64_C(2)<<30);
        if(length<113 || (uint64_t)length>limit || fseek(file,0,SEEK_SET))return -1;
        alea_sha256_t hash;
        alea_sha256_init(&hash);
        unsigned char buffer[65536],actual[32],expected[32];
        size_t left=(size_t)length-32;
        while(left) {
            size_t n=left<sizeof(buffer)?left:sizeof(buffer);
            if(fread(buffer,1,n,file)!=n)return -1;
            alea_sha256_update(&hash,buffer,n);
            left-=n;
        }
        alea_sha256_final(&hash,actual);
        if(fread(expected,1,32,file)!=32 || memcmp(expected,actual,32) || fseek(file,0,SEEK_SET))return -1;
        max_bytes=(size_t)length-32;
        alea_load_profile_end("snapshot_checksum",checksum_start);
    }
    *io=(alea_snapshot_io) {
        .file=file,.read=read,.remaining=max_bytes?max_bytes:(UINT64_C(2)<<30),
        .profile=read && alea_load_profile_start()>=0
    };
    alea_sha256_init(&io->digest);
    return snapshot_header(io,kind);
}
int alea_snapshot_end(alea_snapshot_io* io) {
    if(io->failed)return -1;
    if(io->read && io->remaining)return -1;
    if(io->mapped) {
        double start=io->profile?alea_monotonic_seconds():0;
        int status=mapped_checksum(io);
        if(io->profile)alea_load_profile_end("snapshot_decode_checksum",start);
        return status;
    }
    if(!io->read)flush_snapshot(io);
    if(io->failed)return -1;
    unsigned char digest[32],stored[32];
    alea_sha256_final(&io->digest,digest);
    if(io->read) {
        if(io->profile)fprintf(stderr,"[alea-load-profile] %-28s %.9f s\n",
                              "snapshot_decode_checksum",io->decode_hash_seconds);
        if(fread(stored,1,32,io->file)!=32 || memcmp(digest,stored,32) || fgetc(io->file)!=EOF)io->failed=1;
    }
    else if(fwrite(digest,1,32,io->file)!=32 || fflush(io->file))io->failed=1;
    return io->failed?-1:0;
}
static void box(alea_snapshot_io* io,alea_node_bbox_t* b) {
    SF(b->min_x);
    SF(b->max_x);
    SF(b->min_y);
    SF(b->max_y);
    SF(b->min_z);
    SF(b->max_z);
}
static int validate_node_graph(const alea_system_t* sys,uint32_t id,unsigned char* state,unsigned depth){
    if(id>=sys->nodes.count || depth>1024 || alea_interrupted())return -1;
    if(state[id]==2)return 0;
    if(state[id])return -1;
    state[id]=1;
    const alea_node_t* node=&sys->nodes.data[id];
    alea_operation_t op=ALEA_GET_OPERATION(node);
    if(op!=ALEA_OP_PRIMITIVE){
        if(validate_node_graph(sys,node->operation.left,state,depth+1))return -1;
        if(op!=ALEA_OP_COMPLEMENT && validate_node_graph(sys,node->operation.right,state,depth+1))return -1;
    }
    state[id]=2;
    return 0;
}
int alea_snapshot_system(alea_snapshot_io* io,alea_system_t* sys) {
    double stage_start=io->read?alea_load_profile_start():-1;
    SD(sys->config.abs_tol);
    SD(sys->config.rel_tol);
    SD(sys->config.zero_threshold);
    SU(sys->config.dedup);
    SU(sys->config.log_level);
    SU(sys->config.export_materials);
    SU(sys->config.export_transforms);
    SU(sys->config.universe_depth);
    SU(sys->config.fill_depth);
    SU(sys->config.void_max_depth);
    SD(sys->config.void_min_size);
    SU(sys->config.void_probes_per_axis);
    SD(sys->config.merge_cell_weight);
    SD(sys->config.merge_surface_weight);
    SU(sys->config.merge_max_surfaces);
    SU(sys->config.merge_min_cells);
    SU(sys->config.merge_use_greedy);
    SU(sys->config.void_consolidate);
    SU(sys->config.flatten_max_depth);
    SU(sys->next_inline_transform_id);
    SU(sys->next_auto_surface_id);
    SU(sys->next_auto_material_id);
    SU(sys->next_auto_cell_id);
    SU(sys->source);
    /* Preserve the counters reported by alea_stats(), including load history. */
    SU(sys->stats.unique_primitives);
    SU(sys->stats.dedup_hits);
    SU(sys->stats.dedup_saved_bytes);
    SU(sys->stats.current_memory);
    SU(sys->stats.peak_memory);
    SU(sys->stats.reallocations);
    SU(sys->stats.surfaces_converted);
    SU(sys->stats.failed_surfaces);
    SU(sys->stats.cells_converted);
    size_t n=alea_snapshot_count(io,sys->primitives.count,sizeof(alea_primitive_entry_t));
    for(size_t i=0;i<n && !io->failed;i++) {
        alea_primitive_data_t data= {
            0
        };
        double values[50]= {
            0
        };
        alea_primitive_type_t type=io->read?0:sys->primitives.data[i].type;
        SU(type);
        size_t count=0;
        if(!io->read) {
            if(!alea_primitive_copy_data(sys,(uint32_t)i,&data))io->failed=1;
            count=alea_primitive_to_values(type,&data,values);
        }
        size_t encoded=alea_snapshot_count(io,count,sizeof(double));
        if(encoded>50 || !encoded) {
            io->failed=1;
            break;
        }
        for(size_t k=0;k<encoded;k++)SD(values[k]);
        if(io->read) {
            if(alea_primitive_to_values(type,&data,(double[50]) {
                0
            })!=encoded || alea_primitive_from_values(type,values,&data)) {
                io->failed=1;
                break;
            }
            alea_primitive_entry_t entry= {
                .type=type
            };
            switch(type) {
                case ALEA_PRIMITIVE_PLANE: entry.payload_index=(uint32_t)sys->primitive_planes.count;
                if(alea_vec_push(&sys->primitive_planes,data.plane,alea_plane_data_t)<0)io->failed=1;
                break;
                case ALEA_PRIMITIVE_SPHERE: entry.payload_index=(uint32_t)sys->primitive_spheres.count;
                if(alea_vec_push(&sys->primitive_spheres,data.sphere,alea_sphere_data_t)<0)io->failed=1;
                break;
                case ALEA_PRIMITIVE_CYLINDER_X: entry.payload_index=(uint32_t)sys->primitive_cyl_x.count;
                if(alea_vec_push(&sys->primitive_cyl_x,data.cyl_x,alea_cylinder_x_data_t)<0)io->failed=1;
                break;
                case ALEA_PRIMITIVE_CYLINDER_Y: entry.payload_index=(uint32_t)sys->primitive_cyl_y.count;
                if(alea_vec_push(&sys->primitive_cyl_y,data.cyl_y,alea_cylinder_y_data_t)<0)io->failed=1;
                break;
                case ALEA_PRIMITIVE_CYLINDER_Z: entry.payload_index=(uint32_t)sys->primitive_cyl_z.count;
                if(alea_vec_push(&sys->primitive_cyl_z,data.cyl_z,alea_cylinder_z_data_t)<0)io->failed=1;
                break;
                case ALEA_PRIMITIVE_CONE_X: entry.payload_index=(uint32_t)sys->primitive_cone_x.count;
                if(alea_vec_push(&sys->primitive_cone_x,data.cone_x,alea_cone_x_data_t)<0)io->failed=1;
                break;
                case ALEA_PRIMITIVE_CONE_Y: entry.payload_index=(uint32_t)sys->primitive_cone_y.count;
                if(alea_vec_push(&sys->primitive_cone_y,data.cone_y,alea_cone_y_data_t)<0)io->failed=1;
                break;
                case ALEA_PRIMITIVE_CONE_Z: entry.payload_index=(uint32_t)sys->primitive_cone_z.count;
                if(alea_vec_push(&sys->primitive_cone_z,data.cone_z,alea_cone_z_data_t)<0)io->failed=1;
                break;
                case ALEA_PRIMITIVE_RPP: entry.payload_index=(uint32_t)sys->primitive_boxes.count;
                if(alea_vec_push(&sys->primitive_boxes,data.box,alea_box_data_t)<0)io->failed=1;
                break;
                case ALEA_PRIMITIVE_QUADRIC: entry.payload_index=(uint32_t)sys->primitive_quadrics.count;
                if(alea_vec_push(&sys->primitive_quadrics,data.quadric,alea_quadric_data_t)<0)io->failed=1;
                break;
                case ALEA_PRIMITIVE_TORUS_X: entry.payload_index=(uint32_t)sys->primitive_toruses.count;
                if(alea_vec_push(&sys->primitive_toruses,data.torus,alea_torus_data_t)<0)io->failed=1;
                break;
                case ALEA_PRIMITIVE_TORUS_Y: entry.payload_index=(uint32_t)sys->primitive_toruses.count;
                if(alea_vec_push(&sys->primitive_toruses,data.torus,alea_torus_data_t)<0)io->failed=1;
                break;
                case ALEA_PRIMITIVE_TORUS_Z: entry.payload_index=(uint32_t)sys->primitive_toruses.count;
                if(alea_vec_push(&sys->primitive_toruses,data.torus,alea_torus_data_t)<0)io->failed=1;
                break;
                case ALEA_PRIMITIVE_RCC: entry.payload_index=(uint32_t)sys->primitive_rccs.count;
                if(alea_vec_push(&sys->primitive_rccs,data.rcc,alea_rcc_data_t)<0)io->failed=1;
                break;
                case ALEA_PRIMITIVE_BOX: entry.payload_index=(uint32_t)sys->primitive_box_generals.count;
                if(alea_vec_push(&sys->primitive_box_generals,data.box_general,alea_box_general_data_t)<0)io->failed=1;
                break;
                case ALEA_PRIMITIVE_SPH: entry.payload_index=(uint32_t)sys->primitive_sphs.count;
                if(alea_vec_push(&sys->primitive_sphs,data.sph,alea_sph_data_t)<0)io->failed=1;
                break;
                case ALEA_PRIMITIVE_TRC: entry.payload_index=(uint32_t)sys->primitive_trcs.count;
                if(alea_vec_push(&sys->primitive_trcs,data.trc,alea_trc_data_t)<0)io->failed=1;
                break;
                case ALEA_PRIMITIVE_ELL: entry.payload_index=(uint32_t)sys->primitive_ells.count;
                if(alea_vec_push(&sys->primitive_ells,data.ell,alea_ell_data_t)<0)io->failed=1;
                break;
                case ALEA_PRIMITIVE_REC: entry.payload_index=(uint32_t)sys->primitive_recs.count;
                if(alea_vec_push(&sys->primitive_recs,data.rec,alea_rec_data_t)<0)io->failed=1;
                break;
                case ALEA_PRIMITIVE_WED: entry.payload_index=(uint32_t)sys->primitive_weds.count;
                if(alea_vec_push(&sys->primitive_weds,data.wed,alea_wed_data_t)<0)io->failed=1;
                break;
                case ALEA_PRIMITIVE_RHP: entry.payload_index=(uint32_t)sys->primitive_rhps.count;
                if(alea_vec_push(&sys->primitive_rhps,data.rhp,alea_rhp_data_t)<0)io->failed=1;
                break;
                case ALEA_PRIMITIVE_ARB: entry.payload_index=(uint32_t)sys->primitive_arbs.count;
                if(alea_vec_push(&sys->primitive_arbs,data.arb,alea_arb_data_t)<0)io->failed=1;
                break;
                default:io->failed=1;
                break;
            }
            if(alea_vec_push(&sys->primitives,entry,alea_primitive_entry_t)<0)io->failed=1;
        }
    }
    alea_load_profile_end("snapshot_primitives",stage_start);
    stage_start=io->read?alea_load_profile_start():-1;
    n=alea_snapshot_count(io,sys->nodes.count,sizeof(alea_node_t));
    if(io->read && !io->failed) {
        /* Every node needs at least two header bytes, six float bounds and
         * two child IDs (28 wire bytes). Check before reserving to reject an
         * impossible count without a large speculative allocation. */
        if(n>io->remaining/28 || alea_vec_reserve(&sys->nodes,n,alea_node_t))
            io->failed=1;
    }
    for(size_t i=0;i<n && !io->failed;i++) {
        alea_node_t node= {
            0
        };
        if(!io->read)node=sys->nodes.data[i];
        SU(node.operation.type);
        SU(node.material_id);
        box(io,&node.bbox);
        alea_operation_t op=ALEA_GET_OPERATION(&node);
        if(op==ALEA_OP_PRIMITIVE) {
            SU(node.primitive.primitive_id);
            SU(node.primitive.prim_type);
            SU(node.primitive.sense);
            SU(node.primitive.inverted);
            SU(node.primitive.mc_surface_id);
            if(node.primitive.primitive_id>=sys->primitives.count || sys->primitives.data[node.primitive.primitive_id].type!=node.primitive.prim_type || (node.primitive.sense!=-1 && node.primitive.sense!=1))io->failed=1;
            if(io->read && !io->failed)sys->primitives.data[node.primitive.primitive_id].ref_count++;
        }
        else {
            SU(node.operation.left);
            SU(node.operation.right);
            if(op>ALEA_OP_COMPLEMENT || node.operation.left>=n || (op!=ALEA_OP_COMPLEMENT && node.operation.right>=n))io->failed=1;
        }
        if(io->read && alea_vec_push(&sys->nodes,node,alea_node_t)<0)io->failed=1;
    }
    alea_load_profile_end("snapshot_nodes",stage_start);
    stage_start=io->read?alea_load_profile_start():-1;
    n=alea_snapshot_count(io,sys->surfaces.count,sizeof(alea_surface_entry_t));
    for(size_t i=0;i<n && !io->failed;i++) {
        alea_surface_entry_t v= {
            0
        };
        if(!io->read)v=sys->surfaces.data[i];
        SU(v.mc_surface_id);
        SU(v.primitive_id);
        SU(v.pos_node);
        SU(v.neg_node);
        SU(v.transform_id);
        SU(v.transform_applied);
        SU(v.boundary_type);
        SU(v.periodic_surface_id);
        SU(v.expanded_pos_node);
        SU(v.expanded_neg_node);
        if(v.primitive_id>=sys->primitives.count || v.pos_node>=sys->nodes.count || v.neg_node>=sys->nodes.count)io->failed=1;
        if(io->read && alea_vec_push(&sys->surfaces,v,alea_surface_entry_t)<0)io->failed=1;
    }
    n=alea_snapshot_count(io,sys->transforms.count,sizeof(alea_transform_t));
    for(size_t i=0;i<n && !io->failed;i++) {
        alea_transform_t v= {
            0
        };
        if(!io->read)v=sys->transforms.data[i];
        SU(v.transform_id);
        SU(v.value_count);
        /* Normalized transforms have either translation or all 12 values.
         * Exporters use this count to index the fixed-size data array. */
        if(v.value_count!=3 && v.value_count!=12) {
            io->failed=1;
            break;
        }
        SU(v.degrees);
        SU(v.from_inline);
        for(int k=0;k<12;k++) {
            SD(v.data[k]);
            SD(v.cosines[k]);
        }
        if(io->read && alea_vec_push(&sys->transforms,v,alea_transform_t)<0)io->failed=1;
    }
    alea_load_profile_end("snapshot_surfaces_transforms",stage_start);
    stage_start=io->read?alea_load_profile_start():-1;
    n=alea_snapshot_count(io,sys->materials.count,sizeof(alea_material_t));
    for(size_t i=0;i<n && !io->failed;i++) {
        alea_material_t* v;
        if(io->read) {
            v=alea_vec_push_uninit(&sys->materials,alea_material_t);
            if(!v) {
                io->failed=1;
                break;
            }
            memset(v,0,sizeof(*v));
        }
        else v=&sys->materials.data[i];
        SU(v->material_id);
        SU(v->is_weight_fraction);
        SU(v->has_standard_density);
        SD(v->standard_density);
        SS(v->name);
        SS(v->comments);
        size_t count=alea_snapshot_count(io,v->nuclides.count,sizeof(alea_nuclide_t));
        for(size_t k=0;k<count && !io->failed;k++) {
            alea_nuclide_t* q=io->read?alea_vec_push_uninit(&v->nuclides,alea_nuclide_t):&v->nuclides.data[k];
            if(!q) {
                io->failed=1;
                break;
            }
            if(io->read)memset(q,0,sizeof(*q));
            SU(q->zaid);
            SD(q->fraction);
            SS(q->library);
        }
        count=alea_snapshot_count(io,v->elements.count,sizeof(alea_element_comp_t));
        for(size_t k=0;k<count && !io->failed;k++) {
            alea_element_comp_t* q=io->read?alea_vec_push_uninit(&v->elements,alea_element_comp_t):&v->elements.data[k];
            if(!q) {
                io->failed=1;
                break;
            }
            if(io->read)memset(q,0,sizeof(*q));
            SU(q->atomic_number);
            SD(q->fraction);
            SS(q->library);
        }
        count=alea_snapshot_count(io,v->thermal_laws.count,sizeof(alea_thermal_law_t));
        for(size_t k=0;k<count && !io->failed;k++) {
            alea_thermal_law_t* q=io->read?alea_vec_push_uninit(&v->thermal_laws,alea_thermal_law_t):&v->thermal_laws.data[k];
            if(!q) {
                io->failed=1;
                break;
            }
            if(io->read)memset(q,0,sizeof(*q));
            SU(q->zaid_match);
            SS(q->identifier);
        }
    }
    n=alea_snapshot_count(io,sys->mixtures.count,sizeof(alea_mixture_t));
    for(size_t i=0;i<n && !io->failed;i++) {
        alea_mixture_t* v=io->read?alea_vec_push_uninit(&sys->mixtures,alea_mixture_t):&sys->mixtures.data[i];
        if(!v) {
            io->failed=1;
            break;
        }
        if(io->read)memset(v,0,sizeof(*v));
        SU(v->mixture_id);
        SU(v->mc_material_id);
        SU(v->is_weight_fraction);
        SS(v->name);
        SS(v->comments);
        size_t count=alea_snapshot_count(io,v->components.count,sizeof(alea_mixture_comp_t));
        for(size_t k=0;k<count && !io->failed;k++) {
            alea_mixture_comp_t q= {
                0
            };
            if(!io->read)q=v->components.data[k];
            SU(q.material_id);
            SD(q.fraction);
            if(io->read && alea_vec_push(&v->components,q,alea_mixture_comp_t)<0)io->failed=1;
        }
    }
    alea_load_profile_end("snapshot_materials",stage_start);
    stage_start=io->read?alea_load_profile_start():-1;
    n=alea_snapshot_count(io,sys->cells.count,sizeof(alea_cell_entry_t));
    for(size_t i=0;i<n && !io->failed;i++) {
        alea_cell_entry_t* v=io->read?alea_vec_push_uninit(&sys->cells,alea_cell_entry_t):&sys->cells.data[i];
        if(!v) {
            io->failed=1;
            break;
        }
        if(io->read)alea_cell_entry_init(v);
        SU(v->mc_cell_id);
        SU(v->root_node_id);
        SU(v->original_root_node_id);
        SU(v->material_id);
        SU(v->material_index);
        SU(v->is_mass_density);
        SU(v->has_temperature);
        SU(v->universe_id);
        SU(v->fill_universe);
        SU(v->fill_transform);
        SU(v->lat_type);
        SU(v->lat_outer_universe);
        SU(v->lat_fill_repeating);
        SU(v->lat_fill_zero_element_coords);
        SD(v->density);
        SD(v->temperature);
        SS(v->comments);
        SS(v->inline_comment);
        for(int k=0;k<6;k++)SU(v->lat_fill_dims[k]);
        for(int k=0;k<3;k++) {
            SD(v->lat_pitch[k]);
            SD(v->lat_lower_left[k]);
        }
        SC(v->lat_fill_count,v->lat_fill);
        for(size_t k=0;k<v->lat_fill_count && !io->failed;k++)SU(v->lat_fill[k]);
        if(v->root_node_id>=sys->nodes.count || (v->original_root_node_id!=ALEA_NODE_ID_INVALID && v->original_root_node_id>=sys->nodes.count) || v->material_index < -1 || (v->material_index>=0 && (size_t)v->material_index>=sys->materials.count))io->failed=1;
    }
    n=alea_snapshot_count(io,sys->cell_refs.count,sizeof(alea_cell_ref_t));
    for(size_t i=0;i<n && !io->failed;i++) {
        alea_cell_ref_t v= {
            0
        };
        if(!io->read)v=sys->cell_refs.data[i];
        SU(v.placeholder_node);
        SU(v.referenced_cell_id);
        SU(v.referencing_cell_id);
        if(io->read && alea_vec_push(&sys->cell_refs,v,alea_cell_ref_t)<0)io->failed=1;
    }
    alea_load_profile_end("snapshot_cells",stage_start);
    stage_start=io->read?alea_load_profile_start():-1;
    if(io->read && !io->failed) {
        for(size_t i=0;i<sys->cells.count;i++)if(!cell_hashmap_put(&sys->cell_index,sys->cells.data[i].mc_cell_id,(int)i))io->failed=1;
        for(size_t i=0;i<sys->primitives.count;i++) {
            alea_primitive_data_t data;
            if(!alea_primitive_copy_data(sys,(uint32_t)i,&data)) {
                io->failed=1;
                break;
            }
            primitive_hash_table_insert(sys->primitive_index,(uint32_t)i,alea_compute_primitive_hash(sys->primitives.data[i].type,&data,&sys->config));
            if(sys->primitive_index->entry_count!=i+1) { io->failed=1; break; }
        }
    }
    alea_load_profile_end("snapshot_lookup_indexes",stage_start);
    stage_start=io->read?alea_load_profile_start():-1;
    if(io->read && !io->failed){
        unsigned char* state=alea_calloc(sys->nodes.count?sys->nodes.count:1,1);
        if(!state)io->failed=1;
        for(size_t i=0;i<sys->nodes.count && !io->failed;i++)
            if(validate_node_graph(sys,(uint32_t)i,state,0))io->failed=1;
        free(state);
    }
    alea_load_profile_end("snapshot_graph_validation",stage_start);
    stage_start=io->read?alea_load_profile_start():-1;
    if(!io->failed && alea_snapshot_hierarchy(io,sys))io->failed=1;
    alea_load_profile_end("snapshot_hierarchy",stage_start);
    return io->failed?-1:0;
}
int alea_snapshot_write(const alea_system_t* sys,const char* path) {
    if(!sys || !path)return -1;
    char temporary[4096];
    FILE* file=alea_sibling_tmpfile(path,temporary,sizeof(temporary));
    if(!file)return -1;
    alea_snapshot_io io;
    int rc=alea_snapshot_begin(&io,file,0,0,0);
    if(!rc)rc=alea_snapshot_system(&io,(alea_system_t*)sys);
    if(!rc)rc=alea_snapshot_end(&io);
    if(fclose(file))rc=-1;
    if(!rc)rc=alea_replace_file(temporary,path);
    if(rc)remove(temporary);
    return rc;
}
alea_system_t* alea_snapshot_read(const char* path,size_t max_bytes) {
    alea_snapshot_read_options_t options=ALEA_SNAPSHOT_READ_OPTIONS_INIT;
    options.max_file_bytes=max_bytes;
    return alea_snapshot_read_with_options(path,&options);
}
alea_system_t* alea_snapshot_read_with_options(
    const char* path,const alea_snapshot_read_options_t* options) {
    if(!path)return NULL;
    size_t max_bytes,decoded_bytes;
    if(alea_snapshot_read_limits(options,&max_bytes,&decoded_bytes))return NULL;
    alea_snapshot_io io;
    int rc=alea_snapshot_open_read(&io,path,max_bytes,0);
    alea_alloc_budget_t budget;
    alea_alloc_budget_begin(&budget,decoded_bytes);
    alea_system_t* sys=rc?NULL:alea_create();
    if(!sys || budget.failed)rc=-1;
    double decode_start=rc?-1:alea_load_profile_start();
    if(!rc)rc=alea_snapshot_system(&io,sys);
    if(!rc)rc=alea_snapshot_end(&io);
    if(budget.failed)rc=-1;
    alea_alloc_budget_end(&budget);
    alea_snapshot_profile_allocations(&budget);
    alea_load_profile_end("snapshot_decode",decode_start);
    alea_snapshot_close_read(&io);
    if(rc) {
        alea_destroy(sys);
        if(budget.exceeded)alea_set_error_detail(ALEA_ERR_OUT_OF_MEMORY,
            "Snapshot decoded-allocation budget exceeded (%zu bytes)",decoded_bytes);
        else if(budget.failed)alea_set_error_detail(ALEA_ERR_OUT_OF_MEMORY,
            "Snapshot allocation failed");
        else alea_set_error_detail(ALEA_ERR_FILE_READ,"Invalid or incompatible ALEA snapshot");
        return NULL;
    }
    return sys;
}
