// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
// SPDX-License-Identifier: MPL-2.0
#include "alea_test.h"
#include "alea_mcnp.h"
#include "alea_fingerprint.h"
#include "alea_snapshot.h"
#include "alea_raycast.h"
#include "util/compat.h"
#include "core/alea_system.h"
#include "core/snapshot_internal.h"
#include "util/alea_digest.h"
#include "util/alea_alloc.h"
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>
#if !defined(_WIN32) && !defined(__EMSCRIPTEN__)
#include <pthread.h>
#endif
static const char first[]="identity\n1 1 -1 -10\n2 0 10\n\n10 so 5\n\nm1 1001 1\n";
static const char renumbered[]="identity\n9 7 -1 -20\n8 0 20\n\n20 so 5\n\nm7 1001 1\n";
static void compare_stats(const alea_system_t* a,const alea_system_t* b) {
    alea_stats_t x,y;
    alea_stats(a,&x);
    alea_stats(b,&y);
    ASSERT_EQ(x.unique_primitives,y.unique_primitives);
    ASSERT_EQ(x.dedup_hits,y.dedup_hits);
    ASSERT_EQ(x.dedup_saved_bytes,y.dedup_saved_bytes);
    ASSERT_EQ(x.current_memory,y.current_memory);
    ASSERT_EQ(x.peak_memory,y.peak_memory);
    ASSERT_EQ(x.reallocations,y.reallocations);
    ASSERT_EQ(x.surfaces_converted,y.surfaces_converted);
    ASSERT_EQ(x.failed_surfaces,y.failed_surfaces);
    ASSERT_EQ(x.cells_converted,y.cells_converted);
}
TEST(sha256_standard_vectors) {
    alea_sha256_t h;
    alea_digest_t d;
    char hex[65];
    alea_sha256_init(&h);
    alea_sha256_update(&h,"abc",3);
    alea_sha256_final(&h,d.bytes);
    alea_digest_hex(&d,hex);
    ASSERT_STR_EQ(hex,"ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad");
    alea_sha256_init(&h);
    char block[1000];
    memset(block,'a',sizeof(block));
    for(int i=0;i<1000;i++)alea_sha256_update(&h,block,sizeof(block));
    alea_sha256_final(&h,d.bytes);
    alea_digest_hex(&d,hex);
    ASSERT_STR_EQ(hex,"cdc76e5c9914fb9281a1c7e284d73e67f1809a48a497200e046d39ccc7112cd0");
}
TEST(fingerprint_renumbering_and_content) {
    mcnp_model_t* a=mcnp_load_string(first,0);
    mcnp_model_t* b=mcnp_load_string(renumbered,0);
    ASSERT_NOT_NULL(a);
    ASSERT_NOT_NULL(b);
    alea_fingerprints_t* fa=alea_fingerprints(a->sys);
    alea_fingerprints_t* fb=alea_fingerprints(b->sys);
    ASSERT_NOT_NULL(fa);
    ASSERT_NOT_NULL(fb);
    ASSERT_EQ(memcmp(&fa->cells[0].geometry,&fb->cells[0].geometry,32),0);
    ASSERT_EQ(memcmp(&fa->cells[0].content,&fb->cells[0].content,32),0);
    alea_model_diff_t* diff=alea_fingerprint_diff(fa,fb);
    ASSERT_NOT_NULL(diff);
    ASSERT_EQ(diff->changes[0].kind,ALEA_DIFF_RENUMBERED);
    alea_model_diff_destroy(diff);
    b->sys->cells.data[0].density=2;
    alea_fingerprints_destroy(fb);
    fb=alea_fingerprints(b->sys);
    ASSERT_NOT_NULL(fb);
    ASSERT_EQ(memcmp(&fa->cells[0].geometry,&fb->cells[0].geometry,32),0);
    ASSERT_TRUE(memcmp(&fa->cells[0].content,&fb->cells[0].content,32)!=0);
    alea_fingerprints_destroy(fa);
    alea_fingerprints_destroy(fb);
    mcnp_model_destroy(a);
    mcnp_model_destroy(b);
}
TEST(snapshot_roundtrip_and_mutation) {
    mcnp_model_t* a=mcnp_load_string(first,0);
    ASSERT_NOT_NULL(a);
    char path[4096];
    FILE* temporary=alea_tmpfile(path);
    ASSERT_NOT_NULL(temporary);
    fclose(temporary);
    ASSERT_EQ(alea_prepare_query_acceleration(a->sys),0);
    ASSERT_EQ(alea_snapshot_write(a->sys,path),0);
    ASSERT_NULL(alea_snapshot_read(path,16));
    alea_system_t* b=alea_snapshot_read(path,0);
    ASSERT_NOT_NULL(b);
    compare_stats(a->sys,b);
    alea_query_acceleration_stats_t stats;
    ASSERT_EQ(alea_query_acceleration_stats(b,&stats),0);
    ASSERT_TRUE(stats.built);
    alea_fingerprints_t* fa=alea_fingerprints(a->sys);
    alea_fingerprints_t* fb=alea_fingerprints(b);
    ASSERT_NOT_NULL(fa);
    ASSERT_NOT_NULL(fb);
    ASSERT_EQ(memcmp(&fa->cells[0],&fb->cells[0],sizeof(fa->cells[0])),0);
    ASSERT_EQ(alea_prepare_query_acceleration(b),0);
    int cell=0,mat=0;
    ASSERT_EQ(alea_find_cell_at(b,0,0,0,&cell,&mat),0);
    ASSERT_EQ(cell,1);
    ASSERT_EQ(mat,1);
    int s=alea_sphere_surface(b,99,20,0,0,1);
    ASSERT_TRUE(s>0);
    alea_stats_t before,after;
    alea_stats(a->sys,&before);
    alea_stats(b,&after);
    ASSERT_EQ(after.unique_primitives,before.unique_primitives+1);
    alea_fingerprints_destroy(fa);
    alea_fingerprints_destroy(fb);
    alea_destroy(b);
    FILE* f=fopen(path,"r+b");
    ASSERT_NOT_NULL(f);
    ASSERT_EQ(fseek(f,100,SEEK_SET),0);
    fputc(0xff,f);
    fclose(f);
    b=alea_snapshot_read(path,0);
    ASSERT_NULL(b);
    remove(path);
    mcnp_model_destroy(a);
}
static void compare_export(const mcnp_model_t* a,const mcnp_model_t* b) {
    FILE* x=tmpfile();
    FILE* y=tmpfile();
    ASSERT_NOT_NULL(x);
    ASSERT_NOT_NULL(y);
    ASSERT_EQ(mcnp_export_stream(a,x),0);
    ASSERT_EQ(mcnp_export_stream(b,y),0);
    rewind(x);
    rewind(y);
    int c,d;
    do {
        c=fgetc(x);
        d=fgetc(y);
        ASSERT_EQ(c,d);
    }
    while(c!=EOF);
    fclose(x);
    fclose(y);
}
TEST(mcnp_snapshot_preserves_metadata_and_export) {
    const char* text="metadata\n1 1 -1 -10 imp:n=0.5 tmp=2.5e-8 vol=42 pwt=-1\n2 0 10\n\n10 so 5\n\nm1 1001.80c 1\n";
    mcnp_model_t* a=mcnp_load_string(text,0);
    ASSERT_NOT_NULL(a);
    char path[4096];
    FILE* temp=alea_tmpfile(path);
    ASSERT_NOT_NULL(temp);
    fclose(temp);
    ASSERT_EQ(mcnp_snapshot_write(a,path),0);
    mcnp_model_t* b=mcnp_snapshot_read(path,0);
    ASSERT_NOT_NULL(b);
    compare_export(a,b);
    compare_stats(a->sys,b->sys);
    ASSERT_NEAR(b->cell_params[0].imp_n,.5,0);
    ASSERT_NEAR(b->cell_params[0].vol,42,0);
    int sid=alea_sphere_surface(b->sys,99,20,0,0,1);
    alea_node_id_t root=alea_halfspace(b->sys,sid,-1);
    ASSERT_TRUE(alea_add_cell(b->sys,99,root,-1,0,0)>=0);
    ASSERT_EQ(b->cell_params_count,b->sys->cells.count);
    remove(path);
    mcnp_model_destroy(a);
    mcnp_model_destroy(b);
}
TEST(cache_hit_dependency_change_and_corruption) {
    char input[4096],dependency[4096],cache[4096];
    FILE* file=alea_tmpfile(input);
    ASSERT_NOT_NULL(file);
    FILE* inc=alea_tmpfile(dependency);
    ASSERT_NOT_NULL(inc);
    fputs("10 so 5\n",inc);
    fclose(inc);
    FILE* entry=alea_tmpfile(cache);
    ASSERT_NOT_NULL(entry);
    fclose(entry);
    fprintf(file,"dependencies\n1 1 -1 -10\n2 0 10\n\nREAD FILE=\"%s\"\n\nm1 1001 1\n",dependency);
    fclose(file);
    int hit=-1;
    mcnp_model_t* a=mcnp_load_cached(input,cache,&hit);
    ASSERT_NOT_NULL(a);
    ASSERT_EQ(hit,0);
    mcnp_model_t* b=mcnp_load_cached(input,cache,&hit);
    ASSERT_NOT_NULL(b);
    ASSERT_EQ(hit,1);
    compare_export(a,b);
    compare_stats(a->sys,b->sys);
    mcnp_model_destroy(b);
    inc=fopen(dependency,"w");
    ASSERT_NOT_NULL(inc);
    fputs("10 so 6\n",inc);
    fclose(inc);
    b=mcnp_load_cached(input,cache,&hit);
    ASSERT_NOT_NULL(b);
    ASSERT_EQ(hit,0);
    alea_fingerprints_t* fa=alea_fingerprints(a->sys);
    alea_fingerprints_t* fb=alea_fingerprints(b->sys);
    ASSERT_NOT_NULL(fa);
    ASSERT_NOT_NULL(fb);
    ASSERT_TRUE(memcmp(&fa->cells[0].geometry,&fb->cells[0].geometry,32)!=0);
    alea_fingerprints_destroy(fa);
    alea_fingerprints_destroy(fb);
    mcnp_model_destroy(b);
    entry=fopen(cache,"wb");
    ASSERT_NOT_NULL(entry);
    fputs("broken",entry);
    fclose(entry);
    b=mcnp_load_cached(input,cache,&hit);
    ASSERT_NOT_NULL(b);
    ASSERT_EQ(hit,0);
    mcnp_model_destroy(b);
    /* Failed cache writes still return the freshly loaded model. */
    b=mcnp_load_cached(input,"/nonexistent/alea-cache.bin",&hit);
    ASSERT_NOT_NULL(b);
    ASSERT_EQ(hit,0);
    mcnp_model_destroy(b);
    b=mcnp_load_cached(input,input,&hit);
    ASSERT_NOT_NULL(b);
    ASSERT_EQ(hit,0);
    mcnp_model_destroy(b);
    b=mcnp_load_cached(input,dependency,&hit);
    ASSERT_NOT_NULL(b);
    ASSERT_EQ(hit,0);
    mcnp_model_destroy(b);
    /* A different spelling of the include path must also protect the source. */
    const char* separator=strrchr(dependency,'/');
    if(!separator)separator=strrchr(dependency,'\\');
    ASSERT_NOT_NULL(separator);
    char alias[sizeof(dependency)+2];
    snprintf(alias,sizeof(alias),"%.*s.%c%s",(int)(separator-dependency+1),
             dependency,*separator,separator+1);
    b=mcnp_load_cached(input,alias,&hit);
    ASSERT_NOT_NULL(b);
    ASSERT_EQ(hit,0);
    mcnp_model_destroy(b);
    b=mcnp_load_cached(input,cache,&hit);
    ASSERT_NOT_NULL(b);
    ASSERT_EQ(hit,1);
    mcnp_model_destroy(b);
    remove(input);
    remove(dependency);
    remove(cache);
    mcnp_model_destroy(a);
}
TEST(fingerprint_associativity_order_and_ambiguous_duplicates) {
    const char* x="x\n1 0 -1 -2 -3\n\n1 so 5\n2 px 1\n3 py 1\n\n";
    const char* y="y\n1 0 -3 (-2 -1)\n\n1 so 5\n2 px 1\n3 py 1\n\n";
    mcnp_model_t* a=mcnp_load_string(x,0);
    mcnp_model_t* b=mcnp_load_string(y,0);
    ASSERT_NOT_NULL(a);
    ASSERT_NOT_NULL(b);
    alea_fingerprints_t* fa=alea_fingerprints(a->sys);
    alea_fingerprints_t* fb=alea_fingerprints(b->sys);
    ASSERT_NOT_NULL(fa);
    ASSERT_NOT_NULL(fb);
    ASSERT_EQ(memcmp(&fa->cells[0].geometry,&fb->cells[0].geometry,32),0);
    alea_fingerprints_destroy(fa);
    alea_fingerprints_destroy(fb);
    mcnp_model_destroy(a);
    mcnp_model_destroy(b);
    alea_fingerprint_t left[2]= {
        {
            .id=1
        }, {
            .id=2
        }
    },right[2]= {
        {
            .id=3
        }, {
            .id=4
        }
    };
    alea_fingerprints_t l= {
        .version=ALEA_FINGERPRINT_VERSION,.cell_count=2,.cells=left
    },r= {
        .version=ALEA_FINGERPRINT_VERSION,.cell_count=2,.cells=right
    };
    alea_model_diff_t* diff=alea_fingerprint_diff(&l,&r);
    ASSERT_NOT_NULL(diff);
    ASSERT_EQ(diff->count,4);
    ASSERT_EQ(diff->changes[0].kind,ALEA_DIFF_AMBIGUOUS);
    alea_model_diff_destroy(diff);
}
TEST(snapshot_nested_fill_retains_paths) {
    const char* text="fills\n1 0 -10 fill=5\n2 0 10\n3 1 -1 -20 u=5\n4 0 20 u=5\n\n10 so 10\n20 so 5\n\nm1 1001 1\n";
    mcnp_model_t* a=mcnp_load_string(text,0);
    ASSERT_NOT_NULL(a);
    ASSERT_EQ(alea_prepare_query_acceleration(a->sys),0);
    char path[4096];
    FILE* temp=alea_tmpfile(path);
    ASSERT_NOT_NULL(temp);
    fclose(temp);
    ASSERT_EQ(mcnp_snapshot_write(a,path),0);
    mcnp_model_t* b=mcnp_snapshot_read(path,0);
    ASSERT_NOT_NULL(b);
    ASSERT_EQ(alea_prepare_query_acceleration(b->sys),0);
    for(int k=-12;k<=12;k++) {
        int ca=0,cb=0,ma=0,mb=0;
        ASSERT_EQ(alea_find_cell_at(a->sys,k,.1,0,&ca,&ma),alea_find_cell_at(b->sys,k,.1,0,&cb,&mb));
        ASSERT_EQ(ca,cb);
        ASSERT_EQ(ma,mb);
    }
    compare_export(a,b);
    remove(path);
    mcnp_model_destroy(a);
    mcnp_model_destroy(b);
}
TEST(mixture_fingerprints_and_snapshot) {
    alea_system_t* sys=alea_create();
    int h=alea_add_material(sys,1),o=alea_add_material(sys,2);
    ASSERT_EQ(alea_material_add_nuclide(sys,h,1001,NULL,1),0);
    ASSERT_EQ(alea_material_add_nuclide(sys,o,8016,NULL,1),0);
    int materials[2]= {
        1,2
    };
    double weights[2]= {
        .2,.8
    };
    int mixture=alea_create_mixture(sys,materials,weights,2,7);
    ASSERT_EQ(mixture,7);
    int sphere=alea_sphere_surface(sys,1,0,0,0,5);
    int cell=alea_add_cell(sys,1,alea_halfspace(sys,sphere,-1),h,-1,0);
    ASSERT_TRUE(cell>=0);
    ASSERT_EQ(alea_cell_set_mixture(sys,cell,mixture),0);
    alea_fingerprints_t* before=alea_fingerprints(sys);
    ASSERT_NOT_NULL(before);
    char path[4096];
    FILE* file=alea_tmpfile(path);
    ASSERT_NOT_NULL(file);
    fclose(file);
    ASSERT_EQ(alea_snapshot_write(sys,path),0);
    alea_system_t* restored=alea_snapshot_read(path,0);
    ASSERT_NOT_NULL(restored);
    alea_fingerprints_t* after=alea_fingerprints(restored);
    ASSERT_NOT_NULL(after);
    ASSERT_EQ(memcmp(&before->cells[0].content,&after->cells[0].content,32),0);
    restored->mixtures.data[0].components.data[0].fraction=.3;
    alea_fingerprints_destroy(after);
    after=alea_fingerprints(restored);
    ASSERT_NOT_NULL(after);
    ASSERT_TRUE(memcmp(&before->cells[0].content,&after->cells[0].content,32)!=0);
    alea_fingerprints_destroy(before);
    alea_fingerprints_destroy(after);
    alea_destroy(sys);
    alea_destroy(restored);
    remove(path);
}
TEST(snapshot_integer_wire_boundaries){
    FILE* file=tmpfile();ASSERT_NOT_NULL(file);alea_snapshot_io io;
    uint64_t values[]={0,1,63,64,127,128,UINT32_MAX,UINT64_MAX,UINT64_C(1)<<63};
    ASSERT_EQ(alea_snapshot_begin(&io,file,0,0,7),0);
    for(size_t i=0;i<sizeof(values)/sizeof(*values);i++)alea_snapshot_uint(&io,values[i]);
    ASSERT_EQ(alea_snapshot_end(&io),0);rewind(file);
    ASSERT_EQ(alea_snapshot_begin(&io,file,1,0,7),0);
    for(size_t i=0;i<sizeof(values)/sizeof(*values);i++)ASSERT_TRUE(alea_snapshot_uint(&io,0)==values[i]);
    ASSERT_EQ(alea_snapshot_end(&io),0);fclose(file);
}
TEST(fingerprint_duplicate_terms_collapse_to_child) {
    const char* regions[]={"-1", "-1 -1", "-1:-1", "(-1:-1) (-1 -1)",
                           "-1 -1 -1", "-1:-1:-1"};
    alea_fingerprints_t* base=NULL;
    for(size_t i=0;i<sizeof(regions)/sizeof(*regions);i++) {
        char input[256];
        snprintf(input,sizeof(input),"duplicates\n1 0 %s\n\n1 so 5\n\n",regions[i]);
        mcnp_model_t* model=mcnp_load_string(input,0);
        ASSERT_NOT_NULL(model);
        alea_fingerprints_t* fp=alea_fingerprints(model->sys);
        ASSERT_NOT_NULL(fp);
        if(!base)base=fp;
        else {
            alea_model_diff_t* diff=alea_fingerprint_diff(base,fp);
            ASSERT_NOT_NULL(diff);
            ASSERT_EQ(diff->count,2); /* Cell and enclosing universe. */
            for(size_t j=0;j<diff->count;j++)ASSERT_EQ(diff->changes[j].kind,ALEA_DIFF_UNCHANGED);
            alea_model_diff_destroy(diff);
            alea_fingerprints_destroy(fp);
        }
        mcnp_model_destroy(model);
    }
    alea_fingerprints_destroy(base);
}

/* Change the length preceding a distinctive translation payload, then resign
 * the file. The reader must reject the field even with a correct checksum. */
static void corrupt_transform_count(const char* path,unsigned char wire_count) {
    FILE* file=fopen(path,"rb");
    ASSERT_NOT_NULL(file);
    ASSERT_EQ(fseek(file,0,SEEK_END),0);
    long length=ftell(file);
    ASSERT_TRUE(length>32);
    rewind(file);
    unsigned char* bytes=malloc((size_t)length);
    ASSERT_NOT_NULL(bytes);
    ASSERT_EQ(fread(bytes,1,(size_t)length,file),(size_t)length);
    fclose(file);
    /* Little-endian 1.25 in data[0] and cosines[0], then 2.5 in both arrays. */
    const unsigned char marker[]={
        0,0,0,0,0,0,0xf4,0x3f, 0,0,0,0,0,0,0xf4,0x3f,
        0,0,0,0,0,0,0x04,0x40, 0,0,0,0,0,0,0x04,0x40
    };
    size_t found=0;
    for(size_t i=3;i+sizeof(marker)<=(size_t)length-32;i++) {
        if(memcmp(bytes+i,marker,sizeof(marker)))continue;
        ASSERT_EQ(bytes[i-3],6); /* Zigzag-encoded value_count=3. */
        ASSERT_EQ(bytes[i-2],0); /* degrees */
        ASSERT_EQ(bytes[i-1],0); /* from_inline */
        bytes[i-3]=wire_count;
        found++;
    }
    ASSERT_EQ(found,1);
    alea_sha256_t hash;
    alea_sha256_init(&hash);
    alea_sha256_update(&hash,bytes,(size_t)length-32);
    alea_sha256_final(&hash,bytes+length-32);
    file=fopen(path,"wb");
    ASSERT_NOT_NULL(file);
    ASSERT_EQ(fwrite(bytes,1,(size_t)length,file),(size_t)length);
    ASSERT_EQ(fclose(file),0);
    free(bytes);
}

TEST(snapshot_rejects_invalid_transform_counts) {
    const char* input="transforms\n1 0 -1\n\n1 so 5\n\n"
        "tr1 1.25 2.5 3.75\ntr2 4 5 6 1 0 0 0 1 0 0 0 1\n";
    mcnp_model_t* model=mcnp_load_string(input,0);
    ASSERT_NOT_NULL(model);
    ASSERT_EQ(model->sys->transforms.count,2);
    ASSERT_EQ(model->sys->transforms.data[0].value_count,3);
    ASSERT_EQ(model->sys->transforms.data[1].value_count,12);
    char path[4096];
    FILE* file=alea_tmpfile(path);
    ASSERT_NOT_NULL(file);
    fclose(file);
    /* Both legal lengths must still round-trip and export. */
    ASSERT_EQ(mcnp_snapshot_write(model,path),0);
    mcnp_model_t* restored=mcnp_snapshot_read(path,0);
    ASSERT_NOT_NULL(restored);
    compare_export(model,restored);
    mcnp_model_destroy(restored);
    const int invalid[]={-1,0,1,11,13,63};
    for(size_t i=0;i<sizeof(invalid)/sizeof(*invalid);i++) {
        unsigned char wire=invalid[i]<0?1:(unsigned char)(2*invalid[i]);
        ASSERT_EQ(mcnp_snapshot_write(model,path),0);
        corrupt_transform_count(path,wire);
        ASSERT_NULL(mcnp_snapshot_read(path,0));
        ASSERT_EQ(alea_snapshot_write(model->sys,path),0);
        corrupt_transform_count(path,wire);
        ASSERT_NULL(alea_snapshot_read(path,0));
        model->sys->transforms.data[0].value_count=invalid[i];
        ASSERT_EQ(mcnp_snapshot_write(model,path),-1);
        ASSERT_EQ(alea_snapshot_write(model->sys,path),-1);
        model->sys->transforms.data[0].value_count=3;
    }
    /* A malformed but checksummed cache must fall back to loading the source. */
    char source[4096];
    file=alea_tmpfile(source);
    ASSERT_NOT_NULL(file);
    fputs(input,file);
    ASSERT_EQ(fclose(file),0);
    int hit=-1;
    restored=mcnp_load_cached(source,path,&hit);
    ASSERT_NOT_NULL(restored);
    ASSERT_EQ(hit,0);
    mcnp_model_destroy(restored);
    corrupt_transform_count(path,26); /* value_count=13 */
    restored=mcnp_load_cached(source,path,&hit);
    ASSERT_NOT_NULL(restored);
    ASSERT_EQ(hit,0);
    compare_export(model,restored);
    compare_stats(model->sys,restored->sys);
    mcnp_model_destroy(restored);
    restored=mcnp_load_cached(source,path,&hit);
    ASSERT_NOT_NULL(restored);
    ASSERT_EQ(hit,1);
    compare_stats(model->sys,restored->sys);
    mcnp_model_destroy(restored);
    remove(source);
    remove(path);
    mcnp_model_destroy(model);
}

TEST(snapshot_preserves_dedup_statistics) {
    const char* input="dedup\n1 0 -1\n2 0 2\n\n1 so 5\n2 so 5\n\n";
    mcnp_model_t* model=mcnp_load_string(input,0);
    ASSERT_NOT_NULL(model);
    alea_stats_t stats;
    alea_stats(model->sys,&stats);
    ASSERT_TRUE(stats.unique_primitives>0);
    ASSERT_TRUE(stats.dedup_hits>0);
    ASSERT_TRUE(stats.dedup_saved_bytes>0);
    ASSERT_TRUE(stats.surfaces_converted>0);
    ASSERT_TRUE(stats.cells_converted>0);
    char path[4096];
    FILE* file=alea_tmpfile(path);
    ASSERT_NOT_NULL(file);
    fclose(file);
    ASSERT_EQ(alea_snapshot_write(model->sys,path),0);
    alea_system_t* restored=alea_snapshot_read(path,0);
    ASSERT_NOT_NULL(restored);
    compare_stats(model->sys,restored);
    alea_destroy(restored);
    mcnp_model_destroy(model);
    remove(path);
}

TEST(semantic_composition_profiles) {
    const char* compositions[]={"1001.80c 1 8016.80c 2",
        "8016.80c 4 1001.80c 2", "1001.80c .5 8016.80c 2 1001.80c .5"};
    alea_fingerprint_options_t options={sizeof(options),ALEA_FINGERPRINT_SEMANTIC};
    alea_fingerprint_set_t* base=NULL;
    alea_fingerprints_t* legacy=NULL;
    for(size_t i=0;i<3;i++) {
        char input[256];
        snprintf(input,sizeof(input),"profiles\n1 1 -1 -1\n\n1 so 5\n\nm1 %s\n",compositions[i]);
        mcnp_model_t* m=mcnp_load_string(input,0);ASSERT_NOT_NULL(m);
        double original=m->sys->materials.data[0].nuclides.data[0].fraction;
        alea_fingerprint_set_t* set=alea_fingerprint_set_create(m->sys,&options);ASSERT_NOT_NULL(set);
        ASSERT_EQ(alea_fingerprint_set_values(set)->version,3);
        ASSERT_EQ(alea_fingerprint_set_profile(set),ALEA_FINGERPRINT_SEMANTIC);
        ASSERT_NEAR(original,m->sys->materials.data[0].nuclides.data[0].fraction,0);
        alea_fingerprints_t* exact=alea_fingerprints(m->sys);ASSERT_NOT_NULL(exact);
        if(i==0){base=set;legacy=exact;}
        else {
            alea_model_diff_t* diff=alea_fingerprint_set_diff(base,set);ASSERT_NOT_NULL(diff);
            for(size_t j=0;j<diff->count;j++)ASSERT_EQ(diff->changes[j].kind,ALEA_DIFF_UNCHANGED);
            ASSERT_TRUE(memcmp(&legacy->cells[0].content,&exact->cells[0].content,32)!=0);
            alea_model_diff_destroy(diff);alea_fingerprints_destroy(exact);alea_fingerprint_set_destroy(set);
        }
        mcnp_model_destroy(m);
    }
    const char* changed[]={"1001.81c 1 8016.80c 2", "1001.80c 1 8016.80c 3",
                          "1001.80c -1 8016.80c -2"};
    for(size_t i=0;i<3;i++) {
        char input[256];snprintf(input,sizeof(input),"profiles\n1 1 -1 -1\n\n1 so 5\n\nm1 %s\n",changed[i]);
        mcnp_model_t* m=mcnp_load_string(input,0);ASSERT_NOT_NULL(m);
        alea_fingerprint_set_t* set=alea_fingerprint_set_create(m->sys,&options);ASSERT_NOT_NULL(set);
        alea_model_diff_t* diff=alea_fingerprint_set_diff(base,set);ASSERT_NOT_NULL(diff);
        ASSERT_EQ(diff->changes[0].kind,ALEA_DIFF_CONTENT);
        alea_model_diff_destroy(diff);alea_fingerprint_set_destroy(set);
        options.profile=ALEA_FINGERPRINT_LEGACY;
        set=alea_fingerprint_set_create(m->sys,&options);ASSERT_NOT_NULL(set);
        ASSERT_NULL(alea_fingerprint_set_diff(base,set));
        alea_fingerprint_set_destroy(set);
        options.profile=ALEA_FINGERPRINT_SEMANTIC;
        double invalid[]={0,-1,INFINITY,NAN};
        for(size_t j=0;j<4;j++) {
            for(size_t n=0;n<m->sys->materials.data[0].nuclides.count;n++)
                m->sys->materials.data[0].nuclides.data[n].fraction=invalid[j];
            ASSERT_NULL(alea_fingerprint_set_create(m->sys,&options));
        }
        mcnp_model_destroy(m);
    }
    options.struct_size=0;
    ASSERT_NULL(alea_fingerprint_set_create(NULL,&options));
    alea_fingerprints_destroy(legacy);alea_fingerprint_set_destroy(base);
}
TEST(semantic_mixture_scaling_and_duplicates) {
    alea_system_t* a=alea_create();ASSERT_NOT_NULL(a);
    int h=alea_add_material(a,1),o=alea_add_material(a,2);
    ASSERT_EQ(alea_material_add_nuclide(a,h,1001,NULL,1),0);
    ASSERT_EQ(alea_material_add_nuclide(a,o,8016,NULL,1),0);
    int ids[]={1,2};double fractions[]={1,2};
    ASSERT_EQ(alea_create_mixture(a,ids,fractions,2,7),7);
    int surf=alea_sphere_surface(a,1,0,0,0,5);
    int cell=alea_add_cell(a,1,alea_halfspace(a,surf,-1),h,-1,0);
    ASSERT_TRUE(cell>=0);ASSERT_EQ(alea_cell_set_mixture(a,cell,7),0);
    alea_fingerprint_options_t options={sizeof(options),ALEA_FINGERPRINT_SEMANTIC};
    alea_fingerprint_set_t* before=alea_fingerprint_set_create(a,&options);ASSERT_NOT_NULL(before);
    a->mixtures.data[0].components.data[0].fraction=DBL_MAX/2;
    a->mixtures.data[0].components.data[1].fraction=DBL_MAX;
    alea_fingerprint_set_t* after=alea_fingerprint_set_create(a,&options);ASSERT_NOT_NULL(after);
    alea_model_diff_t* diff=alea_fingerprint_set_diff(before,after);ASSERT_NOT_NULL(diff);
    ASSERT_EQ(diff->changes[0].kind,ALEA_DIFF_UNCHANGED);
    alea_model_diff_destroy(diff);alea_fingerprint_set_destroy(after);
    a->mixtures.data[0].components.data[0].fraction=.5;
    a->mixtures.data[0].components.data[1].fraction=2;
    alea_mixture_comp_t duplicate={1,.5};
    ASSERT_TRUE(alea_vec_push(&a->mixtures.data[0].components,duplicate,alea_mixture_comp_t)>=0);
    after=alea_fingerprint_set_create(a,&options);ASSERT_NOT_NULL(after);
    diff=alea_fingerprint_set_diff(before,after);ASSERT_NOT_NULL(diff);
    ASSERT_EQ(diff->changes[0].kind,ALEA_DIFF_UNCHANGED);
    alea_model_diff_destroy(diff);alea_fingerprint_set_destroy(before);alea_fingerprint_set_destroy(after);alea_destroy(a);
}
TEST(snapshot_short_inputs_and_owned_restore) {
    mcnp_model_t* model=mcnp_load_string(first,0);
    ASSERT_NOT_NULL(model);
    char path[4096];
    FILE* file=alea_tmpfile(path);
    ASSERT_NOT_NULL(file);fclose(file);
    ASSERT_EQ(alea_snapshot_write(model->sys,path),0);
    file=fopen(path,"rb");ASSERT_NOT_NULL(file);
    ASSERT_EQ(fseek(file,0,SEEK_END),0);
    long length=ftell(file);ASSERT(length>113);
    unsigned char* bytes=malloc((size_t)length);ASSERT_NOT_NULL(bytes);
    rewind(file);ASSERT_EQ(fread(bytes,1,(size_t)length,file),(size_t)length);fclose(file);
    size_t sizes[]={0,1,112,(size_t)length-1};
    for(size_t i=0;i<sizeof(sizes)/sizeof(*sizes);i++) {
        file=fopen(path,"wb");ASSERT_NOT_NULL(file);
        ASSERT_EQ(fwrite(bytes,1,sizes[i],file),sizes[i]);fclose(file);
        ASSERT_NULL(alea_snapshot_read(path,0));
        ASSERT_NULL(mcnp_snapshot_read(path,0));
    }
    file=fopen(path,"wb");ASSERT_NOT_NULL(file);
    ASSERT_EQ(fwrite(bytes,1,(size_t)length,file),(size_t)length);fclose(file);free(bytes);
    alea_system_t* restored=alea_snapshot_read(path,0);ASSERT_NOT_NULL(restored);
    /* Windows cannot delete an open mapped file: this also checks handle release.
     * Subsequent mutation/queries must use owned arrays after unmapping. */
    ASSERT_EQ(remove(path),0);
    ASSERT(alea_sphere_surface(restored,99,20,0,0,1)>=0);
    ASSERT_EQ(alea_prepare_query_acceleration(restored),0);
    int cell,material;
    ASSERT_EQ(alea_find_cell_at(restored,0,0,0,&cell,&material),0);
    ASSERT_EQ(cell,1);
    alea_destroy(restored);mcnp_model_destroy(model);
}

TEST(snapshot_numeric_fields_cross_buffer_boundaries) {
    /* Literal wire values provide an independent decoder oracle. */
    const unsigned char wires[3][10]={
        {0,0,0xa0,0xbf}, /* float -1.25 */
        {0,0,0,0,0,0,0xf4,0xbf}, /* double -1.25 */
        {255,255,255,255,255,255,255,255,255,1} /* zigzag INT64_MIN */
    };
    const size_t lengths[]={4,8,10};
    unsigned char padding[65536]={0x0a}; /* Must survive Windows binary I/O. */
    for(int type=0;type<3;type++)for(size_t gap=0;gap<=10;gap++) {
        char path[4096];FILE* file=alea_tmpfile(path);ASSERT_NOT_NULL(file);
        alea_snapshot_io io;
        ASSERT_EQ(alea_snapshot_begin(&io,file,0,0,7),0);
        size_t n=sizeof(io.buffer)-io.buffered-gap;
        alea_snapshot_bytes(&io,padding,n);
        alea_snapshot_bytes(&io,(void*)wires[type],lengths[type]);
        ASSERT_EQ(alea_snapshot_end(&io),0);fclose(file);
        ASSERT_EQ(alea_snapshot_open_read(&io,path,0,7),0);
        alea_snapshot_bytes(&io,padding,n);
        if(type==0)ASSERT(alea_snapshot_float(&io,0)==-1.25f);
        else if(type==1)ASSERT(alea_snapshot_double(&io,0)==-1.25);
        else ASSERT(alea_snapshot_uint(&io,0)==(UINT64_C(1)<<63));
        ASSERT_EQ(alea_snapshot_end(&io),0);
        alea_snapshot_close_read(&io);ASSERT_EQ(remove(path),0);
    }
}

TEST(snapshot_numeric_fast_paths_reject_bad_fields) {
    const unsigned char invalid[5][10]={
        {0x80,0}, /* noncanonical integer */
        {255,255,255,255,255,255,255,255,255,2}, /* integer overflow */
        {0,0,0xa0,0xbf}, /* incomplete double */
        {0,0,0xc0,0x7f}, /* float NaN */
        {0,0,0,0,0,0,0xf8,0x7f} /* double NaN */
    };
    const size_t lengths[]={2,10,4,4,8};
    for(int type=0;type<5;type++) {
        char path[4096];FILE* file=alea_tmpfile(path);ASSERT_NOT_NULL(file);
        alea_snapshot_io io;unsigned char padding[64]={0};
        ASSERT_EQ(alea_snapshot_begin(&io,file,0,0,7),0);
        alea_snapshot_bytes(&io,padding,sizeof(padding));
        alea_snapshot_bytes(&io,(void*)invalid[type],lengths[type]);
        ASSERT_EQ(alea_snapshot_end(&io),0);fclose(file);
        ASSERT_EQ(alea_snapshot_open_read(&io,path,0,7),0);
        alea_snapshot_bytes(&io,padding,sizeof(padding));
        if(type==2 || type==4)(void)alea_snapshot_double(&io,0);
        else if(type==3)(void)alea_snapshot_float(&io,0);
        else (void)alea_snapshot_uint(&io,0);
        ASSERT(io.failed);ASSERT_EQ(alea_snapshot_end(&io),-1);
        alea_snapshot_close_read(&io);ASSERT_EQ(remove(path),0);
    }
}

TEST(snapshot_decoded_budget_covers_core_wrapper_and_hierarchy) {
    mcnp_model_t* model=mcnp_load_string(first,0);
    ASSERT_NOT_NULL(model);
    ASSERT_EQ(alea_prepare_query_acceleration(model->sys),0);
    char path[4096];FILE* file=alea_tmpfile(path);ASSERT_NOT_NULL(file);fclose(file);
    size_t core_bytes=0;
    for(int wrapper=0;wrapper<2;wrapper++) {
        ASSERT_EQ(wrapper?mcnp_snapshot_write(model,path):alea_snapshot_write(model->sys,path),0);
        alea_alloc_budget_t observed;
        alea_alloc_budget_begin(&observed,SIZE_MAX);
        mcnp_model_t* m=wrapper?mcnp_snapshot_read(path,0):NULL;
        alea_system_t* s=wrapper?NULL:alea_snapshot_read(path,0);
        alea_alloc_budget_end(&observed);
        ASSERT(!observed.failed);ASSERT(observed.requested>0);
        if(wrapper) { ASSERT_NOT_NULL(m);ASSERT(observed.requested>core_bytes); }
        else { ASSERT_NOT_NULL(s);core_bytes=observed.requested; }
        mcnp_model_destroy(m);alea_destroy(s);
        alea_snapshot_read_options_t options=ALEA_SNAPSHOT_READ_OPTIONS_INIT;
        options.max_decoded_bytes=observed.requested;
        m=wrapper?mcnp_snapshot_read_with_options(path,&options):NULL;
        s=wrapper?NULL:alea_snapshot_read_with_options(path,&options);
        if(wrapper) { ASSERT_NOT_NULL(m);compare_export(model,m); }
        else { ASSERT_NOT_NULL(s);compare_stats(model->sys,s); }
        mcnp_model_destroy(m);alea_destroy(s);
        /* One byte below total demand must fail, even with a large file limit. */
        options.max_file_bytes=SIZE_MAX;
        options.max_decoded_bytes--;
        ASSERT_NULL(wrapper?(void*)mcnp_snapshot_read_with_options(path,&options):
                            (void*)alea_snapshot_read_with_options(path,&options));
        ASSERT_NOT_NULL(strstr(alea_error(),"decoded-allocation budget exceeded"));
        options.max_decoded_bytes=SIZE_MAX;options.max_file_bytes=1;
        ASSERT_NULL(wrapper?(void*)mcnp_snapshot_read_with_options(path,&options):
                            (void*)alea_snapshot_read_with_options(path,&options));
        options.struct_size=0;
        ASSERT_NULL(wrapper?(void*)mcnp_snapshot_read_with_options(path,&options):
                            (void*)alea_snapshot_read_with_options(path,&options));
        /* Failure must not leave a thread-local limit active on the next read. */
        m=wrapper?mcnp_snapshot_read(path,0):NULL;
        s=wrapper?NULL:alea_snapshot_read(path,0);
        if(wrapper)ASSERT_NOT_NULL(m);else ASSERT_NOT_NULL(s);
        mcnp_model_destroy(m);alea_destroy(s);
    }
    ASSERT_EQ(remove(path),0);mcnp_model_destroy(model);
}

TEST(snapshot_allocation_failure_cleanup_and_retry) {
    mcnp_model_t* model=mcnp_load_string(first,0);ASSERT_NOT_NULL(model);
    ASSERT_EQ(alea_prepare_query_acceleration(model->sys),0);
    char path[4096];FILE* file=alea_tmpfile(path);ASSERT_NOT_NULL(file);fclose(file);
    ASSERT_EQ(mcnp_snapshot_write(model,path),0);
    int succeeded=0;
    for(long allocation=0;allocation<512;allocation++) {
        alea_alloc_fail_after(allocation);
        mcnp_model_t* copy=mcnp_snapshot_read(path,0);
        alea_alloc_fail_after(-1);
        if(copy) { mcnp_model_destroy(copy);succeeded=1;break; }
        copy=mcnp_snapshot_read(path,0);ASSERT_NOT_NULL(copy);
        compare_export(model,copy);mcnp_model_destroy(copy);
    }
    ASSERT(succeeded);
    ASSERT_EQ(remove(path),0);mcnp_model_destroy(model);
}

TEST(snapshot_budget_rejects_checksummed_expansion_and_aggregate_demand) {
    /* A compact count passes the wire-size checks but would expand 64 bytes
     * into 64 transform structures. Sign the payload so checksum rejection
     * cannot mask the allocation gate. */
    char path[4096];FILE* file=alea_tmpfile(path);ASSERT_NOT_NULL(file);
    alea_snapshot_io stream;alea_snapshot_io* io=&stream;unsigned char padding[128]={0};
    ASSERT_EQ(alea_snapshot_begin(io,file,0,0,7),0);
    (void)alea_snapshot_uint(io,64);
    alea_snapshot_bytes(io,padding,sizeof(padding));
    ASSERT_EQ(alea_snapshot_end(io),0);fclose(file);
    ASSERT_EQ(alea_snapshot_open_read(io,path,0,7),0);
    alea_alloc_budget_t budget;
    alea_alloc_budget_begin(&budget,128);
    size_t count=0;alea_transform_t* transforms=NULL;
    SC(count,transforms);
    alea_alloc_budget_end(&budget);
    ASSERT(io->failed);ASSERT(budget.exceeded);ASSERT_EQ(budget.requested,0);
    ASSERT_NULL(transforms);ASSERT_EQ(count,0);
    alea_snapshot_close_read(io);ASSERT_EQ(remove(path),0);

    /* Free does not refund the cumulative limit; full realloc capacity counts.
     * Nested scopes cannot bypass the caller's enclosing budget. */
    alea_alloc_budget_t inner;
    alea_alloc_budget_begin(&budget,100);
    void* p=alea_malloc(40);ASSERT_NOT_NULL(p);
    p=alea_realloc(p,50);ASSERT_NOT_NULL(p);free(p);
    alea_alloc_budget_begin(&inner,1000);
    ASSERT_NULL(alea_malloc(11));
    alea_alloc_budget_end(&inner);alea_alloc_budget_end(&budget);
    ASSERT(budget.exceeded);ASSERT(inner.failed);ASSERT_EQ(budget.requested,90);
    ASSERT_EQ(inner.requested,0);
}

TEST(snapshot_budget_cache_miss_falls_back_to_source) {
    char source[4096],cache[4096];
    FILE* file=alea_tmpfile(source);ASSERT_NOT_NULL(file);fputs(first,file);fclose(file);
    file=alea_tmpfile(cache);ASSERT_NOT_NULL(file);fclose(file);
    int hit=-1;
    mcnp_model_t* model=mcnp_load_cached(source,cache,&hit);
    ASSERT_NOT_NULL(model);ASSERT_EQ(hit,0);mcnp_model_destroy(model);
    model=mcnp_load_cached(source,cache,&hit);
    ASSERT_NOT_NULL(model);ASSERT_EQ(hit,1);mcnp_model_destroy(model);
    alea_snapshot_read_options_t options=ALEA_SNAPSHOT_READ_OPTIONS_INIT;
    options.max_decoded_bytes=1;
    model=mcnp_load_cached_with_options(source,cache,&hit,&options);
    ASSERT_NOT_NULL(model);ASSERT_EQ(hit,0);mcnp_model_destroy(model);
    options.struct_size=0;
    ASSERT_NULL(mcnp_load_cached_with_options(source,cache,&hit,&options));
    ASSERT_EQ(hit,0);
    model=mcnp_load_cached(source,cache,&hit);
    ASSERT_NOT_NULL(model);ASSERT_EQ(hit,1);mcnp_model_destroy(model);
    ASSERT_EQ(remove(source),0);ASSERT_EQ(remove(cache),0);
}

TEST(snapshot_checksummed_mutations_remain_allocation_bounded) {
    mcnp_model_t* model=mcnp_load_string(first,0);ASSERT_NOT_NULL(model);
    ASSERT_EQ(alea_prepare_query_acceleration(model->sys),0);
    char path[4096];FILE* file=alea_tmpfile(path);ASSERT_NOT_NULL(file);fclose(file);
    ASSERT_EQ(mcnp_snapshot_write(model,path),0);mcnp_model_destroy(model);
    file=fopen(path,"rb");ASSERT_NOT_NULL(file);
    ASSERT_EQ(fseek(file,0,SEEK_END),0);long length=ftell(file);ASSERT(length>106);
    rewind(file);
    unsigned char* original=malloc((size_t)length);
    unsigned char* changed=malloc((size_t)length);
    ASSERT_NOT_NULL(original);ASSERT_NOT_NULL(changed);
    ASSERT_EQ(fread(original,1,(size_t)length,file),(size_t)length);fclose(file);
    alea_snapshot_read_options_t options=ALEA_SNAPSHOT_READ_OPTIONS_INIT;
    options.max_decoded_bytes=64*1024;
    /* Exercise counts/indices and payload fields beyond the valid header.
     * Re-sign every mutation so it reaches decoding and allocation checks. */
    for(size_t sample=0;sample<128;sample++) {
        memcpy(changed,original,(size_t)length);
        size_t offset=74+sample*((size_t)length-32-74)/128;
        changed[offset]=(sample%3==0)?255:(sample%3==1)?126:0;
        alea_sha256_t hash;alea_sha256_init(&hash);
        alea_sha256_update(&hash,changed,(size_t)length-32);
        alea_sha256_final(&hash,changed+length-32);
        file=fopen(path,"wb");ASSERT_NOT_NULL(file);
        ASSERT_EQ(fwrite(changed,1,(size_t)length,file),(size_t)length);fclose(file);
        alea_alloc_budget_t observed;alea_alloc_budget_begin(&observed,SIZE_MAX);
        model=mcnp_snapshot_read_with_options(path,&options);
        alea_alloc_budget_end(&observed);
        ASSERT(observed.requested<=options.max_decoded_bytes);
        mcnp_model_destroy(model);
    }
    free(changed);free(original);ASSERT_EQ(remove(path),0);
}

#if !defined(_WIN32) && !defined(__EMSCRIPTEN__)
static void* budget_thread_read(void* path) { return mcnp_snapshot_read(path,0); }
TEST(snapshot_allocation_budget_is_thread_local) {
    mcnp_model_t* model=mcnp_load_string(first,0);ASSERT_NOT_NULL(model);
    char path[4096];FILE* file=alea_tmpfile(path);ASSERT_NOT_NULL(file);fclose(file);
    ASSERT_EQ(mcnp_snapshot_write(model,path),0);mcnp_model_destroy(model);
    alea_alloc_budget_t budget;alea_alloc_budget_begin(&budget,1);
    pthread_t thread;void* result=NULL;
    int created=pthread_create(&thread,NULL,budget_thread_read,path);
    int joined=created?created:pthread_join(thread,&result);
    alea_alloc_budget_end(&budget);
    ASSERT_EQ(created,0);ASSERT_EQ(joined,0);ASSERT_NOT_NULL(result);
    ASSERT_EQ(budget.requested,0);ASSERT(!budget.failed);
    mcnp_model_destroy(result);ASSERT_EQ(remove(path),0);
}
#endif

TEST_MAIN()
