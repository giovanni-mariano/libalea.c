// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
// SPDX-License-Identifier: MPL-2.0
#include "mcnp_model.h"
#include "core/snapshot_internal.h"
#include "util/compat.h"
#include "util/alea_mcnp_input.h"
#include "util/load_profile.h"
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <sys/stat.h>
#include <errno.h>
#include <math.h>
#ifdef _WIN32
#include <direct.h>
#include <sys/utime.h>
#else
#include <utime.h>
#endif
static int metadata(alea_snapshot_io* io,mcnp_model_t* model) {
    SS(model->name);
    SS(model->title);
    SS(model->comments);
    size_t n=alea_snapshot_count(io,model->cell_params_count,sizeof(mcnp_cell_params_t));
    if(n!=model->sys->cells.count || n!=model->cell_params_count) {
        io->failed=1;
        return -1;
    }
    for(size_t i=0;i<n && !io->failed;i++) {
        SS(model->cell_names[i]);
        mcnp_cell_params_t* p=&model->cell_params[i];
        SD(p->imp_n);
        SD(p->imp_p);
        SD(p->imp_e);
#define FIELD_double SD
#define FIELD_int SU
#define FIELD(name,type,keyword,precision) do { FIELD_##type(p->name); SU(p->has_##name); } while(0);
        MCNP_CELL_SIMPLE_PARAMS(FIELD)
#undef FIELD
#undef FIELD_int
#undef FIELD_double
        SU(p->trcl);
        SU(p->trcl_inline);
        SU(p->trcl_degrees);
        SU(p->trcl_count);
        SU(p->trcl_inline_index);
        SU(p->fill_transform_inline);
        SU(p->fill_transform_degrees);
        SU(p->fill_transform_count);
        SU(p->fill_transform_index);
        SU(p->like_cell_id);
        SU(p->has_imp_n);
        SU(p->has_imp_p);
        SU(p->has_imp_e);
        SU(p->has_trcl);
        SU(p->has_mat);
        SU(p->has_rho);
        SU(p->scoped.energy_cutoff_particles);
        SU(p->scoped.secondary_state_particles);
        for(int k=0;k<ALEA_PARTICLE_COUNT;k++) {
            SD(p->scoped.energy_cutoff[k]);
            SU(p->scoped.secondary_state[k]);
        }
        SC(p->scoped.detector_probability_count,p->scoped.detector_probabilities);
        for(size_t k=0;k<p->scoped.detector_probability_count && !io->failed;k++) {
            SU(p->scoped.detector_probabilities[k].tally);
            SD(p->scoped.detector_probabilities[k].probability);
        }
    }
    SC(model->inline_transform_count,model->inline_transforms);
    if(io->read)model->inline_transform_capacity=model->inline_transform_count;
    for(size_t i=0;i<model->inline_transform_count && !io->failed;i++) {
        mcnp_inline_transform_t* p=&model->inline_transforms[i];
        SU(p->count);
        SU(p->degrees);
        for(int k=0;k<13;k++)SD(p->values[k]);
        if(p->count>13)io->failed=1;
    }
    SU(model->export_config.surface_policy);
    SU(model->export_config.trcl_mode);
    SU(model->export_config.transform_mode);
    SU(model->export_config.mcnp_max_col);
    SU(model->export_config.mcnp_cont_indent);
    return io->failed?-1:0;
}
static int write_model(const mcnp_model_t* model,const char* path,const unsigned char* key) {
    if(!model || !path)return -1;
    char temporary[4096];
    FILE* file=alea_sibling_tmpfile(path,temporary,sizeof(temporary));
    if(!file)return -1;
    alea_snapshot_io io;
    int rc=alea_snapshot_begin(&io,file,0,0,key?2:1);
    if(key)alea_snapshot_bytes(&io,(void*)key,32);
    if(!rc)rc=alea_snapshot_system(&io,model->sys);
    if(!rc)rc=metadata(&io,(mcnp_model_t*)model);
    if(!rc)rc=alea_snapshot_end(&io);
    if(fclose(file))rc=-1;
    if(!rc)rc=alea_replace_file(temporary,path);
    if(rc)remove(temporary);
    return rc;
}
static mcnp_model_t* read_model(const char* path,const alea_snapshot_read_options_t* options,const unsigned char* key) {
    if(!path)return NULL;
    size_t max_bytes,decoded_bytes;
    if(alea_snapshot_read_limits(options,&max_bytes,&decoded_bytes))return NULL;
    alea_snapshot_io io;
    int rc=alea_snapshot_open_read(&io,path,max_bytes,key?2:1);
    unsigned char stored[32];
    if(key && !rc) {
        alea_snapshot_bytes(&io,stored,32);
        if(memcmp(stored,key,32))rc=-1;
    }
    alea_alloc_budget_t budget;
    alea_alloc_budget_begin(&budget,decoded_bytes);
    alea_system_t* sys=rc?NULL:alea_create();
    double decode_start=rc?-1:alea_load_profile_start();
    if(!sys || budget.failed)rc=-1;
    if(!rc)rc=alea_snapshot_system(&io,sys);
    mcnp_model_t* model=rc?NULL:mcnp_model_wrap(sys);
    if(!model)rc=-1;
    if(model)model->owns_sys=1;
    if(!rc)rc=metadata(&io,model);
    if(!rc)rc=alea_snapshot_end(&io);
    if(budget.failed)rc=-1;
    alea_alloc_budget_end(&budget);
    alea_snapshot_profile_allocations(&budget);
    alea_load_profile_end("snapshot_decode",decode_start);
    alea_snapshot_close_read(&io);
    if(rc) {
        if(model)mcnp_model_destroy(model);
        else alea_destroy(sys);
        if(budget.exceeded)alea_set_error_detail(ALEA_ERR_OUT_OF_MEMORY,
            "Snapshot decoded-allocation budget exceeded (%zu bytes)",decoded_bytes);
        else if(budget.failed)alea_set_error_detail(ALEA_ERR_OUT_OF_MEMORY,
            "Snapshot allocation failed");
        else alea_set_error_detail(ALEA_ERR_FILE_READ,"Invalid or incompatible MCNP snapshot");
        return NULL;
    }
    return model;
}
int mcnp_snapshot_write(const mcnp_model_t* model,const char* path) {
    return write_model(model,path,NULL);
}
mcnp_model_t* mcnp_snapshot_read(const char* path,size_t max_bytes) {
    alea_snapshot_read_options_t options=ALEA_SNAPSHOT_READ_OPTIONS_INIT;
    options.max_file_bytes=max_bytes;
    return mcnp_snapshot_read_with_options(path,&options);
}
mcnp_model_t* mcnp_snapshot_read_with_options(
    const char* path,const alea_snapshot_read_options_t* options) {
    return read_model(path,options,NULL);
}
typedef struct {
    char* path;
    uint64_t size;
    time_t touched;
}
cache_entry;
static int oldest_cache_entry(const void* a,const void* b) {
    const cache_entry* x=a;
    const cache_entry* y=b;
    return x->touched!=y->touched?(x->touched>y->touched?1:-1):strcmp(x->path,y->path);
}
static void trim_cache_directory(const char* directory,const char* current) {
    double gb=8;
    const char* setting=getenv("ALEA_CACHE_GB");
    if(setting) {
        char* end;
        double value=strtod(setting,&end);
        if(end!=setting && !*end && isfinite(value) && value>=0 && value<1000000)gb=value;
    }
    uint64_t budget=(uint64_t)(gb*1073741824.0),total=0;
    alea_dir_t* dir=alea_dir_open(directory);
    if(!dir)return;
    cache_entry* entries=NULL;
    size_t count=0;
    const char* name;
    while((name=alea_dir_next(dir))!=NULL) {
        size_t n=strlen(name);
        if(n!=80 || strncmp(name,"alea-cache-",11) || strcmp(name+n-5,".abin"))continue;
        int valid=1;
        for(size_t i=11;i<75;i++)if(!isxdigit((unsigned char)name[i]))valid=0;
        if(!valid || count==100000)continue;
        char path[4096];
        int length=snprintf(path,sizeof(path),"%s/%s",directory,name);
        struct stat info;
        if(length<0 || (size_t)length>=sizeof(path) || stat(path,&info) || info.st_size<0)continue;
        if(!strcmp(path,current))continue;
        cache_entry* grown=realloc(entries,(count+1)*sizeof(*entries));
        if(!grown)break;
        entries=grown;
        char* copy=malloc(strlen(path)+1);
        if(!copy)break;
        strcpy(copy,path);
        entries[count++]=(cache_entry) {
            copy,(uint64_t)info.st_size,info.st_mtime
        };
        total+=(uint64_t)info.st_size;
    }
    alea_dir_close(dir);
    struct stat info;
    uint64_t current_size=!stat(current,&info) && info.st_size>0?(uint64_t)info.st_size:0;
    total+=current_size;
    if(count)qsort(entries,count,sizeof(*entries),oldest_cache_entry);
    for(size_t i=0;i<count;i++) {
        if(total>budget && !remove(entries[i].path))total-=entries[i].size;
        free(entries[i].path);
    }
    free(entries);
    if(total>budget && current_size)remove(current);
}
static int cache_directory_path(const char* directory,const unsigned char key[32],char path[4096]) {
#ifdef _WIN32
    if(_mkdir(directory) && errno!=EEXIST)return -1;
#else
    if(mkdir(directory,0700) && errno!=EEXIST)return -1;
#endif
    char hex[65];
    static const char digits[]="0123456789abcdef";
    for(int i=0;i<32;i++) {
        hex[2*i]=digits[key[i]>>4];
        hex[2*i+1]=digits[key[i]&15];
    }
    hex[64]=0;
    int n=snprintf(path,4096,"%s/alea-cache-%s.abin",directory,hex);
    return n<0 || n>=4096?-1:0;
}
static void touch_cache_entry(const char* path) {
#ifdef _WIN32
    (void)_utime(path,NULL);
#else
    (void)utime(path,NULL);
#endif
}
mcnp_model_t* mcnp_load_cached(const char* filename,const char* cache_path,int* cache_hit) {
    return mcnp_load_cached_with_options(filename,cache_path,cache_hit,NULL);
}
mcnp_model_t* mcnp_load_cached_with_options(const char* filename,const char* cache_path,
    int* cache_hit,const alea_snapshot_read_options_t* options) {
    if(cache_hit)*cache_hit=0;
    if(!filename)return NULL;
    size_t file_limit,decoded_limit;
    if(alea_snapshot_read_limits(options,&file_limit,&decoded_limit))return NULL;
    double expansion_start=alea_load_profile_start();
    char* input=NULL;
    size_t length=0;
    int conflict=0;
    if(alea_read_mcnp_input(filename,cache_path,&input,&length,&conflict)) {
        alea_set_error_detail(ALEA_ERR_FILE_READ,"Cannot expand MCNP input dependencies");
        return NULL;
    }
    alea_sha256_t h;
    alea_load_profile_end("cache_source_expansion",expansion_start);
    double hash_start=alea_load_profile_start();
    unsigned char key[32];
    alea_sha256_init(&h);
    alea_digest_string(&h,"ALEA MCNP source cache v1");
    alea_digest_string(&h,filename);
    alea_sha256_update(&h,input,length);
    alea_digest_string(&h,getenv("ALEA_HIER_BLAS_THRESHOLD"));
    alea_digest_string(&h,getenv("ALEA_UNIVERSE_POINT_BVH_THRESHOLD"));
    alea_digest_string(&h,alea_snapshot_engine_id());
    alea_sha256_final(&h,key);
    alea_load_profile_end("cache_source_hash",hash_start);
    const char* directory=cache_path?NULL:getenv("ALEA_CACHE_DIR");
    char directory_path[4096];
    const char* enabled=getenv("ALEA_CACHE");
    int allowed=!enabled || strcmp(enabled,"off");
    if(directory && allowed && !cache_directory_path(directory,key,directory_path))cache_path=directory_path;
    int use_cache=cache_path && !conflict && allowed;
    mcnp_model_t* model=use_cache?read_model(cache_path,options,key):NULL;
    if(model) {
        if(cache_hit)*cache_hit=1;
        if(directory) {
            touch_cache_entry(cache_path);
            trim_cache_directory(directory,cache_path);
        }
        free(input);
        return model;
    }
    if(use_cache)alea_clear_error_detail(); /* A rejected entry is an ordinary miss. */
    model=mcnp_load_string(input,length);
    free(input);
    if(model) {
        if(mcnp_model_set_name(model,filename)) {
            mcnp_model_destroy(model);
            return NULL;
        }
        if(use_cache && !alea_prepare_query_acceleration(model->sys)) {
            if(!write_model(model,cache_path,key) && directory)trim_cache_directory(directory,cache_path);
        }
    }
    return model;
}
