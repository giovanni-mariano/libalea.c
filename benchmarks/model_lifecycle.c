// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
// SPDX-License-Identifier: MPL-2.0
/* One process per scenario: high-water RSS must not combine unrelated loads. */
#include "alea_mcnp.h"
#include "util/ray_profile.h"
#include "alea_raycast.h"
#include "alea_fingerprint.h"
#include "core/snapshot_internal.h"
#include "core/alea_universe.h"
#include "core/alea_spatial_hier.h"
#include "primitives/bbox.h"
#include "raycast/recovery_planes.h"
#include "util/compat.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#ifndef _WIN32
#include <sys/resource.h>
#endif

typedef struct { char kind; double v[7]; } query_t;
static int hierarchical_rays;
static query_t* read_workload(const char* path,size_t* count) {
    FILE* f=fopen(path,"r");
    if(!f)return NULL;
    query_t* queries=NULL;
    char line[1024],extra;
    size_t capacity=0;
    *count=0;
    while(fgets(line,sizeof(line),f)) {
        query_t q={0};
        if(line[0]=='#' || line[0]=='\n')continue;
        int n;
        if(sscanf(line," %c",&q.kind)!=1)goto fail;
        if(q.kind=='p') {
            n=sscanf(line," %c %lf %lf %lf %c",&q.kind,&q.v[0],&q.v[1],&q.v[2],&extra);
            if(n!=4)goto fail;
        } else if(q.kind=='r') {
            n=sscanf(line," %c %lf %lf %lf %lf %lf %lf %lf %c",&q.kind,
                     &q.v[0],&q.v[1],&q.v[2],&q.v[3],&q.v[4],&q.v[5],&q.v[6],&extra);
            if(n!=8)goto fail;
        } else goto fail;
        for(int k=0;k<n-1;k++)if(!isfinite(q.v[k]))goto fail;
        if(q.kind=='r' && (q.v[6]<=0 ||
           (q.v[3]==0 && q.v[4]==0 && q.v[5]==0)))goto fail;
        if(*count==1000000)goto fail;
        if(*count==capacity) {
            capacity=capacity?2*capacity:256;
            void* grown=realloc(queries,capacity*sizeof(*queries));
            if(!grown)goto fail;
            queries=grown;
        }
        queries[(*count)++]=q;
    }
    if(ferror(f) || !*count)goto fail;
    fclose(f);return queries;
fail:
    fclose(f);free(queries);return NULL;
}

static int point_query(alea_system_t* sys,const query_t* q,alea_sha256_t* hash) {
    int cell=-1,material=-1;
    int rc=alea_find_cell_at(sys,q->v[0],q->v[1],q->v[2],&cell,&material);
    alea_digest_u64(hash,(uint64_t)(int64_t)rc);
    alea_digest_u64(hash,(uint64_t)(int64_t)cell);
    alea_digest_u64(hash,(uint64_t)(int64_t)material);
    return cell>=0;
}

static int ray_query(alea_system_t* sys,const query_t* q,
                     alea_raycast_result_t* result,alea_sha256_t* hash,size_t* segments) {
    const double* v=q->v;
    int rc=hierarchical_rays
        ? alea_raycast_hier_fast_segments(sys,v[0],v[1],v[2],v[3],v[4],v[5],v[6],result)
        : alea_raycast(sys,v[0],v[1],v[2],v[3],v[4],v[5],v[6],result);
    if(rc)return -1;
    size_t count=alea_raycast_segment_count(result);
    alea_digest_u64(hash,count);*segments+=count;
    for(size_t i=0;i<count;i++) {
        double enter,exit,density;int cell,material,s0,s1;uint8_t flags;
        if(alea_raycast_segment_get(result,i,&enter,&exit,&cell,&material,&density,&s0,&s1) ||
           alea_raycast_segment_resolution_flags(result,i,&flags))return -1;
        if(alea_digest_double(hash,enter) || alea_digest_double(hash,exit) ||
           alea_digest_double(hash,density))return -1;
        alea_digest_u64(hash,(uint64_t)(int64_t)cell);
        alea_digest_u64(hash,(uint64_t)(int64_t)material);
        alea_digest_u64(hash,(uint64_t)(int64_t)s0);
        alea_digest_u64(hash,(uint64_t)(int64_t)s1);
        alea_digest_u64(hash,flags);
    }
    return 0;
}

/* Experimental opt-in. Only replace by an intersected, certified enclosure;
 * leave unsupported roots untouched and invalidate all dependent indexes. */
static size_t mixed_bounds(alea_system_t* sys) {
    size_t changed=0;
    for(size_t i=0;i<sys->cells.count;i++) {
        alea_node_id_t root=sys->cells.data[i].root_node_id;
        if(root>=sys->nodes.count)continue;
        alea_bbox_t certified;
        if(alea_tighten_bbox_mixed_constraints(sys,root,&certified))continue;
        alea_bbox_t old=alea_node_bbox_get(&sys->nodes.data[root].bbox);
        if(!alea_bbox_is_valid(&old))continue;
        alea_bbox_t box=alea_bbox_intersection(&old,&certified);
        alea_node_bbox_t packed;
        alea_node_bbox_set(&packed,&box);
        if(alea_bbox_is_valid(&box) && memcmp(&packed,&sys->nodes.data[root].bbox,sizeof(packed))) {
            sys->nodes.data[root].bbox=packed;changed++;
        }
    }
    if(changed)alea_system_invalidate_query_caches(sys,ALEA_CACHE_ALL & ~ALEA_CACHE_CELL_SURFACES);
    return changed;
}

static void peak_rss(void) {
#ifndef _WIN32
    struct rusage usage;
    if(!getrusage(RUSAGE_SELF,&usage)) {
#ifdef __APPLE__
        printf("%llu",(unsigned long long)usage.ru_maxrss);
#else
        printf("%llu",(unsigned long long)usage.ru_maxrss*1024);
#endif
        return;
    }
#endif
    printf("null");
}

/* Exact ordered acceleration metadata, outside the timed replay. Hash fields,
 * not struct padding; aliases and opposite senses must remain distinguishable. */
static void acceleration_digest(const alea_system_t* sys,char hex[65]) {
    alea_sha256_t hash;alea_digest_t result;
    alea_sha256_init(&hash);
    alea_digest_u64(&hash,sys->cells.count);
    for(size_t i=0;i<sys->cells.count;i++) {
        const alea_cell_entry_t* cell=&sys->cells.data[i];
        alea_digest_u64(&hash,(uint64_t)(int64_t)cell->mc_cell_id);
        alea_digest_u64(&hash,cell->surface_index_count);
        for(size_t k=0;k<cell->surface_index_count;k++)
            alea_digest_u64(&hash,cell->surface_indices[k]);
        alea_digest_u64(&hash,cell->neighbor_count);
        for(size_t k=0;k<cell->neighbor_count;k++) {
            const alea_cell_neighbor_t* n=&cell->neighbors[k];
            alea_digest_u64(&hash,(uint64_t)(int64_t)n->surface_id);
            alea_digest_u64(&hash,n->surface_index);
            alea_digest_u64(&hash,(uint64_t)(int64_t)n->neighbor_cell_id);
            alea_digest_u64(&hash,n->neighbor_index);
            alea_digest_u64(&hash,(uint64_t)(int64_t)n->our_sense);
        }
    }
    alea_digest_u64(&hash,sys->surfaces.count);
    for(size_t i=0;i<sys->surfaces.count;i++) {
        size_t begin=sys->surface_cell_offsets[i],end=sys->surface_cell_offsets[i+1];
        alea_digest_u64(&hash,end-begin);
        for(size_t k=begin;k<end;k++) {
            alea_digest_u64(&hash,sys->surface_cell_refs[k].cell_index);
            alea_digest_u64(&hash,(uint64_t)(int64_t)sys->surface_cell_refs[k].sense);
        }
    }
    alea_sha256_final(&hash,result.bytes);alea_digest_hex(&result,hex);
}

int main(int argc,char** argv) {
    if(argc<5 || argc>9) {
        fprintf(stderr,"Usage: model_lifecycle MODE MODEL WORKLOAD AUX [--mixed-bounds] [--hierarchical-rays] [--candidate-profile] [--ray-profile]\n"
                "MODE: fresh, cache, snapshot, write, edit; AUX: cache/snapshot path (fresh: -)\n");
        return 2;
    }
    int use_mixed=0,profile_candidates=0,profile_rays=0;
    for(int i=5;i<argc;i++) {
        if(!strcmp(argv[i],"--mixed-bounds"))use_mixed=1;
        else if(!strcmp(argv[i],"--hierarchical-rays"))hierarchical_rays=1;
        else if(!strcmp(argv[i],"--candidate-profile"))profile_candidates=1;
        else if(!strcmp(argv[i],"--ray-profile"))profile_rays=1;
        else return 2;
    }
#ifndef ALEA_RAY_PROFILE
    if(profile_rays) { fprintf(stderr,"--ray-profile requires an ALEA_RAY_PROFILE build\n");return 2; }
#endif
    const char* mode=argv[1];
    if(strcmp(mode,"fresh") && strcmp(mode,"cache") && strcmp(mode,"snapshot") &&
       strcmp(mode,"write") && strcmp(mode,"edit"))return 2;
    size_t count=0;
    query_t* queries=read_workload(argv[3],&count);
    if(!queries){fprintf(stderr,"Invalid workload: %s\n",argv[3]);return 2;}
    mcnp_model_t *model=NULL,*previous=NULL;
    alea_raycast_result_t* result=NULL;
    uint64_t *point_cells=NULL,*ray_cells=NULL;
    int status=1,cache_hit=0;
    double start=alea_monotonic_seconds();
    if(!strcmp(mode,"edit")) {
        previous=mcnp_snapshot_read(argv[4],0);
        if(!previous)goto done;
    }
    double previous_loaded=alea_monotonic_seconds();
    if(!strcmp(mode,"snapshot"))model=mcnp_snapshot_read(argv[4],0);
    else if(!strcmp(mode,"cache"))model=mcnp_load_cached(argv[2],argv[4],&cache_hit);
    else model=mcnp_load(argv[2]);
    double loaded=alea_monotonic_seconds();
    if(!model)goto done;
    size_t changed=use_mixed?mixed_bounds(model->sys):0;
    double bounded=alea_monotonic_seconds();
    unsigned flags[]={ALEA_CACHE_UNIVERSE,ALEA_CACHE_CELL_SURFACES,ALEA_CACHE_HIER_SPATIAL,
                      ALEA_CACHE_SURFACE_BVH,ALEA_CACHE_ADJACENCY,
                      ALEA_CACHE_RECOVERY_CHILDREN,ALEA_CACHE_RECOVERY_SURFACES};
    const char* names[]={"universe","cell_surfaces","hierarchy","surface_bvh","adjacency",
                         "recovery_children","recovery_surfaces"};
    double stages[7];
    for(size_t k=0;k<7;k++) {
        double t=alea_monotonic_seconds();
        if(alea_system_prepare_query_caches(model->sys,flags[k]))goto done;
        stages[k]=alea_monotonic_seconds()-t;
    }
    double ready=alea_monotonic_seconds();
    fprintf(stderr,"[alea-load-profile] ready %.9f s\n",ready-start);
    if(!strcmp(mode,"write") && mcnp_snapshot_write(model,argv[4]))goto done;
    double written=alea_monotonic_seconds();
    result=alea_raycast_result_create();if(!result)goto done;
    alea_sha256_t points,rays,warmup;alea_digest_t pd,rd;
    alea_sha256_init(&warmup);
    size_t warm_segments=0;
    double first_start=alea_monotonic_seconds();
    if(queries[0].kind=='p')point_query(model->sys,&queries[0],&warmup);
    else if(ray_query(model->sys,&queries[0],result,&warmup,&warm_segments))goto done;
    double first_s=alea_monotonic_seconds()-first_start;
    if(profile_candidates) {
        point_cells=calloc(model->sys->cells.count?model->sys->cells.count:1,sizeof(*point_cells));
        ray_cells=calloc(model->sys->cells.count?model->sys->cells.count:1,sizeof(*ray_cells));
        if(!point_cells || !ray_cells)goto done;
    }
    alea_sha256_init(&points);alea_sha256_init(&rays);
    size_t np=0,nr=0,hits=0,segments=0;
    alea_universe_point_bvh_stats_reset();
    alea_hier_debug_candidates_reset();
    alea_hier_debug_candidate_cells_set(model->sys,point_cells,model->sys->cells.count);
    double query_start=alea_monotonic_seconds();
    for(size_t i=0;i<count;i++)if(queries[i].kind=='p') {
        hits+=point_query(model->sys,&queries[i],&points);np++;
    }
    double point_s=alea_monotonic_seconds()-query_start;
    fprintf(stderr,"[alea-load-profile] point_replay %.9f s\n",point_s);
    alea_universe_point_bvh_stats_t point_stats;
    alea_universe_point_bvh_stats_get(&point_stats);
    size_t point_candidates=alea_hier_debug_candidates_get();
    alea_hier_debug_candidates_reset();
    alea_hier_debug_candidate_cells_set(model->sys,ray_cells,model->sys->cells.count);
#ifdef ALEA_RAY_PROFILE
    if(profile_rays && alea_ray_profile_start())goto done;
#endif
    query_start=alea_monotonic_seconds();
    for(size_t i=0;i<count;i++)if(queries[i].kind=='r') {
#ifdef ALEA_RAY_PROFILE
        alea_ray_profile.ray_ordinal = nr;
#endif
        if(ray_query(model->sys,&queries[i],result,&rays,&segments))goto done;
        nr++;
    }
    double ray_s=alea_monotonic_seconds()-query_start;
#ifdef ALEA_RAY_PROFILE
    if(profile_rays)alea_ray_profile_stop();
#endif
    alea_hier_debug_candidate_cells_set(NULL,NULL,0);
    alea_sha256_final(&points,pd.bytes);alea_sha256_final(&rays,rd.bytes);
    char ph[65],rh[65],ah[65];alea_digest_hex(&pd,ph);alea_digest_hex(&rd,rh);
    acceleration_digest(model->sys,ah);
    printf("{\"schema\":2,\"engine\":\"%s\",\"mode\":\"%s\",\"cache_hit\":%d,"
           "\"cells\":%zu,\"nodes\":%zu,\"bounds_changed\":%zu,\"previous_load_s\":%.9g,"
           "\"load_s\":%.9g,\"bounds_s\":%.9g,\"prepare_s\":%.9g,\"ready_s\":%.9g,"
           "\"write_s\":%.9g,\"first_query_s\":%.9g,\"points\":%zu,\"point_hits\":%zu,"
           "\"rays\":%zu,\"segments\":%zu,\"point_digest\":\"%s\",\"ray_digest\":\"%s\","
           "\"point_s\":%.9g,\"ray_s\":%.9g,\"legacy_point_bbox_tests\":%zu,"
           "\"legacy_point_exact_tests\":%zu,\"point_candidates\":%zu,"
           "\"ray_point_candidates\":%zu,\"peak_rss_bytes\":",
           alea_snapshot_engine_id(),mode,cache_hit,alea_cell_count(model->sys),
           model->sys->nodes.count,changed,previous_loaded-start,loaded-previous_loaded,
           bounded-loaded,ready-bounded,ready-start,written-ready,first_s,np,hits,nr,segments,
           ph,rh,point_s,ray_s,point_stats.bvh_bbox_tests,point_stats.exact_cell_tests,point_candidates,
           alea_hier_debug_candidates_get());
    peak_rss();printf(",\"ray_api\":\"%s\",\"acceleration_digest\":\"%s\",\"preparation_stages_s\":{",
                     hierarchical_rays?"hierarchical-segments":"global",ah);
    for(size_t k=0;k<7;k++)printf("%s\"%s\":%.9g",k?",":"",names[k],stages[k]);
    printf("}");
    if(profile_candidates) {
        printf(",\"candidate_profile\":[");
        int separator=0;
        for(size_t i=0;i<model->sys->cells.count;i++) {
            if(!point_cells[i] && !ray_cells[i])continue;
            const alea_cell_entry_t* cell=&model->sys->cells.data[i];
            uint64_t mask=0;
            for(size_t k=0;k<cell->surface_index_count;k++) {
                uint32_t si=cell->surface_indices[k];
                if(si>=model->sys->surfaces.count)continue;
                uint32_t pi=model->sys->surfaces.data[si].primitive_id;
                if(pi<model->sys->primitives.count) {
                    unsigned type=model->sys->primitives.data[pi].type;
                    if(type<64)mask|=UINT64_C(1)<<type;
                }
            }
            printf("%s{\"cell_index\":%zu,\"cell_id\":%d,\"universe\":%d,"
                   "\"point_candidates\":%llu,\"ray_candidates\":%llu,"
                   "\"surfaces\":%zu,\"primitive_mask\":%llu}",separator?",":"",
                   i,cell->mc_cell_id,cell->universe_id,
                   (unsigned long long)point_cells[i],(unsigned long long)ray_cells[i],
                   cell->surface_index_count,(unsigned long long)mask);
            separator=1;
        }
        printf("]");
    }
#ifdef ALEA_RAY_PROFILE
    if(profile_rays) {
        static const char* names[]={"other","hierarchy","membership","intersection","recovery","recovery_scan","recovery_cache"};
        printf(",\"ray_profile\":{\"timing\":\"exclusive-instrumented-wall\",\"stages\":[");
        for(int i=0;i<RAY_STAGES;++i)printf("%s{\"name\":\"%s\",\"seconds\":%.9g,\"calls\":%llu}",
            i?",":"",names[i],alea_ray_profile.seconds[i],(unsigned long long)alea_ray_profile.calls[i]);
        printf("],\"csg_nodes\":%llu,\"csg_primitives\":%llu,\"dropped_root_samples\":%llu,\"top_roots\":[",
            (unsigned long long)alea_ray_profile.csg_nodes,(unsigned long long)alea_ray_profile.csg_primitives,
            (unsigned long long)alea_ray_profile.dropped_roots);
        for(unsigned top=0;top<32;++top) {
            ray_profile_root *best=NULL;
            for(unsigned j=0;j<RAY_PROFILE_ROOTS;++j) {
                ray_profile_root *r=&alea_ray_profile.roots[j];
                if(r->used && (!best || r->seconds>best->seconds))best=r;
            }
            if(!best)break;
            printf("%s{\"root_node\":%u,\"calls\":%llu,\"seconds\":%.9g,\"cell_ids\":[",
                   top?",":"",best->root,(unsigned long long)best->calls,best->seconds);
            int sep=0;
            for(size_t j=0;j<model->sys->cells.count;++j)if(model->sys->cells.data[j].root_node_id==best->root)
                {printf("%s%d",sep?",":"",model->sys->cells.data[j].mc_cell_id);sep=1;}
            printf("]}");best->used=0;
        }
        printf("],\"recovery_root_builds\":%llu,\"recovery_root_tests\":%llu,\"recovery_root_queries\":%llu,\"plane_index_groups\":%llu,\"plane_index_pruned\":%llu,\"recovery_cache_hits\":%llu,\"dropped_list_samples\":%llu,\"scan_types\":[",
            (unsigned long long)alea_ray_profile.recovery_root_builds,
            (unsigned long long)alea_ray_profile.recovery_root_tests,
            (unsigned long long)alea_ray_profile.recovery_root_queries,
            (unsigned long long)alea_ray_profile.plane_index_groups,
            (unsigned long long)alea_ray_profile.plane_index_pruned,
            (unsigned long long)alea_ray_profile.recovery_cache_hits,
            (unsigned long long)alea_ray_profile.dropped_lists);
        for (unsigned i=0;i<RAY_PROFILE_TYPES;++i)
            printf("%s{\"type\":%u,\"tests\":%llu,\"window_roots\":%llu}", i?",":"",i,
                (unsigned long long)alea_ray_profile.scan_tests[i],
                (unsigned long long)alea_ray_profile.scan_window_roots[i]);
        printf("],\"fill_traversal\":{\"nearest_passes\":%llu,\"nearer_hits\":%llu,"
               "\"child_visits\":%llu,\"duplicate_skips\":%llu,"
               "\"bbox_rejections\":%llu,\"descents\":%llu},",
               (unsigned long long)alea_ray_profile.fill_passes,
               (unsigned long long)alea_ray_profile.fill_improved,
               (unsigned long long)alea_ray_profile.fill_children,
               (unsigned long long)alea_ray_profile.fill_duplicates,
               (unsigned long long)alea_ray_profile.fill_bbox_rejected,
               (unsigned long long)alea_ray_profile.fill_descents);
        static const char* owners[] = {
            [RAY_OWNER_NEIGHBOR] = "neighbor", [RAY_OWNER_CACHED_PATH] = "cached_path",
            [RAY_OWNER_PATH_DESCENT] = "path_descent", [RAY_OWNER_ROOT] = "root",
            [RAY_OWNER_FULL] = "full_lookup"
        };
        printf("\"dropped_recoveries\":%llu,\"owner_probe_recoveries\":[",
               (unsigned long long)alea_ray_profile.dropped_recoveries);
        for (size_t i=0; i<alea_ray_profile.recovery_count; i++) {
            const ray_profile_recovery_event* e = &alea_ray_profile.recoveries[i];
            printf("%s{\"ray\":%zu,\"attempt\":%d,\"source\":\"%s\",\"path_length\":%d,"
                   "\"root_cell\":%d,\"owner_cell\":%d,\"owner_universe\":%d,\"occurrence\":\"%016llx\","
                   "\"current_t\":%.17g,\"proposed_t\":%.17g,\"proposed_surface\":%d,"
                   "\"recovered_t\":%.17g,\"recovered_surface\":%d}",
                   i ? "," : "", e->ray, e->attempt, owners[e->source], e->path_length,
                   e->root_cell, e->owner_cell, e->owner_universe, (unsigned long long)e->occurrence,
                   e->current_t, e->proposed_t, e->proposed_surface, e->recovered_t, e->recovered_surface);
        }
        printf("],");
        printf("\"indexed_root_lists\":%llu,\"roots_sorted\":%llu,"
               "\"prepared_card_checks_reused\":%llu,\"cached_non_torus_candidates\":%llu,"
               "\"root_build_types\":[",
               (unsigned long long)alea_ray_profile.indexed_root_lists,
               (unsigned long long)alea_ray_profile.roots_sorted,
               (unsigned long long)alea_ray_profile.prepared_card_checks_reused,
               (unsigned long long)alea_ray_profile.cached_root_candidates);
        for (unsigned i=0;i<RAY_PROFILE_TYPES;++i)
            printf("%s{\"type\":%u,\"tests\":%llu,\"window_roots\":%llu}", i?",":"",i,
                   (unsigned long long)alea_ray_profile.root_tests_by_type[i],
                   (unsigned long long)alea_ray_profile.roots_retained_by_type[i]);
        printf("],\"top_scan_lists\":[");
        for (unsigned top=0;top<32;++top) {
            ray_profile_list *best=NULL;
            for (unsigned j=0;j<RAY_PROFILE_LISTS;++j) {
                ray_profile_list *r=&alea_ray_profile.lists[j];
                if (r->calls && (!best || r->tests>best->tests)) best=r;
            }
            if (!best) break;
            printf("%s{\"surface_count\":%zu,\"calls\":%llu,\"tests\":%llu,\"universe_ids\":[",
                top?",":"",best->count,(unsigned long long)best->calls,(unsigned long long)best->tests);
            int sep=0;
            for (size_t j=0;j<model->sys->universes.count;++j) {
                const alea_universe_t *u=&model->sys->universes.data[j];
                if (u->recovery_surfaces==best->indices) {printf("%s%d",sep?",":"",u->universe_id);sep=1;}
            }
            printf("],\"cell_ids\":[");sep=0;
            for (size_t j=0;j<model->sys->cells.count;++j) {
                const alea_cell_entry_t *c=&model->sys->cells.data[j];
                if (c->surface_indices==best->indices) {printf("%s%d",sep?",":"",c->mc_cell_id);sep=1;}
            }
            printf("]}");best->calls=0;
        }
        printf("]}");
    }
#endif
    size_t plane_bytes=0, plane_groups=0, plane_members=0;
    size_t plane_blocks=0, block_members=0, quadric_cards=0, quadric_bytes=0;
    for (size_t i=0;i<model->sys->universes.count;++i) {
        const recovery_plane_index_t* index=model->sys->universes.data[i].recovery_planes;
        if (!index) continue;
        size_t plane_storage = index->quadric_count ? index->quadric_offset : index->bytes;
        plane_bytes += plane_storage;
        quadric_bytes += index->bytes - plane_storage;
        quadric_cards += index->quadric_count;
        plane_groups += index->groups;
        plane_members += index->members;
        plane_blocks += index->blocks;
        block_members += index->block_members;
    }
    printf(",\"recovery_plane_index\":{\"bytes\":%zu,\"groups\":%zu,\"members\":%zu,"
           "\"blocks\":%zu,\"block_members\":%zu}",
           plane_bytes,plane_groups,plane_members,plane_blocks,block_members);
    printf(",\"recovery_quadric_packets\":{\"cards\":%zu,\"bytes\":%zu}",quadric_cards,quadric_bytes);
    puts("}");status=0;
done:
#ifdef ALEA_RAY_PROFILE
    alea_ray_profile_free();
#endif
    alea_hier_debug_candidate_cells_set(NULL,NULL,0);
    free(point_cells);free(ray_cells);
    if(status)fprintf(stderr,"Lifecycle scenario failed: %s\n",alea_error());
    alea_raycast_result_destroy(result);mcnp_model_destroy(previous);
    mcnp_model_destroy(model);free(queries);return status;
}
