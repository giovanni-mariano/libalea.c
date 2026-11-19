// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
// SPDX-License-Identifier: MPL-2.0
/* Benchmark prepared parent-interval point and ray queries. */
#include "alea.h"
#include "alea_mcnp.h"
#include "alea_raycast.h"
#include "core/alea_system.h"
#include "core/alea_eval.h"
#include "core/alea_spatial_hier.h"
#include "core/snapshot_internal.h"
#include "util/compat.h"
#include "interval_image.h"
#include "interval_memory.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct {char kind; double v[7];} query;

int main(int argc,char** argv) {
    bool regional = argc == 6 && !strcmp(argv[5], "--regional");
    if (regional) {
        fprintf(stderr, "The parent-interval experiment requires full cell expressions.\n");
        return 2;
    }
    bool grid = argc == 12 && !strcmp(argv[5], "--point-grid");
    alea_bbox_t domain = {0};
    if (grid) {
        double* limits[] = {&domain.min_x, &domain.max_x, &domain.min_y,
                            &domain.max_y, &domain.min_z, &domain.max_z};
        for (int i = 0; i < 6; i++) {
            char* end;
            *limits[i] = strtod(argv[i+6], &end);
            if (end == argv[i+6] || *end || !isfinite(*limits[i])) return 2;
        }
    }
    if((argc!=5 && !grid) || (strcmp(argv[1],"fresh") && strcmp(argv[1],"write") && strcmp(argv[1],"snapshot")))return 2;
    FILE* f=fopen(argv[3],"r");if(!f)return 2;
    query* qs=NULL; size_t nq=0,np=0,nr=0; char line[1024],extra;
    while(fgets(line,sizeof(line),f)) {
        query q={0}; int n;
        if(sscanf(line," %c",&q.kind)!=1 || q.kind=='#')continue;
        if(q.kind=='p')n=sscanf(line," %c %lf %lf %lf %c",&q.kind,&q.v[0],&q.v[1],&q.v[2],&extra);
        else if(q.kind=='r')n=sscanf(line," %c %lf %lf %lf %lf %lf %lf %lf %c",&q.kind,&q.v[0],&q.v[1],&q.v[2],&q.v[3],&q.v[4],&q.v[5],&q.v[6],&extra);
        else return 2;
        if(n!=(q.kind=='p'?4:8))return 2;
        for(int j=0;j<n-1;j++)if(!isfinite(q.v[j]))return 2;
        if(q.kind=='r' && (q.v[6]<=0 || (q.v[3]==0 && q.v[4]==0 && q.v[5]==0)))return 2;
        if(nq>=1000000)return 2;
        void* grown=realloc(qs,(nq+1)*sizeof(*qs));if(!grown)return 1;qs=grown;qs[nq++]=q;
        if(q.kind=='p')np++;else nr++;
    }
    if(ferror(f) || !np || !nr)return 2;
    fclose(f);
    double start=alea_monotonic_seconds();
    mcnp_model_t* model=!strcmp(argv[1],"snapshot")?mcnp_snapshot_read(argv[4],0):mcnp_load(argv[2]);
    if(!model){fprintf(stderr,"load: %s\n",alea_error());return 1;}
    double load_s=alea_monotonic_seconds()-start;
    interval_memory("model decoded");
    /* The private tracer and its full-expression fallback both use the
     * universe hierarchy. The global walker additionally needs recovery,
     * adjacency and surface BVHs. Keep snapshot writes on the public path. */
    bool private_grid = grid && strcmp(argv[1],"write") != 0;
    if (private_grid) {
        if (alea_system_prepare_query_caches(model->sys,ALEA_CACHE_HIER_SPATIAL)) return 1;
    } else if (alea_prepare_query_acceleration(model->sys)) return 1;
    double grid_s = 0;
    alea_point_grid_stats_t grid_stats = {0};
    if (grid && !private_grid) {
        double begin = alea_monotonic_seconds();
        if (alea_prepare_world_point_grid(model->sys, &domain, &grid_stats)) return 1;
        grid_s = alea_monotonic_seconds() - begin;
    }
    double idx_start = alea_monotonic_seconds();
    interval_memory("base acceleration ready");
    const char* restore = getenv("ALEA_INTERVAL_IMAGE_READ");
    const char* save = getenv("ALEA_INTERVAL_IMAGE_WRITE");
    if ((restore || save) && (!grid || strcmp(argv[1],"snapshot") || (restore && save))) return 2;
    interval_index_t* index = restore ? interval_index_restore(model->sys,&domain,restore)
                                      : grid && getenv("ALEA_INTERVAL_SELECTIVE") ? alea_prepared_queries_create(model->sys, &domain)
                                      : interval_index_build(model->sys, grid ? &domain : NULL, false);
    if (!index) return 1;
    double idx_s = alea_monotonic_seconds() - idx_start;
    interval_memory("regional index ready");
    if (save && interval_index_save(index,save)) return 1;
    double ready_s=alea_monotonic_seconds()-start;
    start=alea_monotonic_seconds();
    if(!strcmp(argv[1],"write") && mcnp_snapshot_write(model,argv[4]))return 1;
    double write_s=alea_monotonic_seconds()-start;
    interval_owner_t* ps=calloc(np,sizeof(*ps));
    alea_raycast_result_t** rs=calloc(nr,sizeof(*rs));
    interval_point_t* point_queries=calloc(np,sizeof(*point_queries));
    if(!ps || !rs || !point_queries)return 1;
    for(size_t i=0,k=0;i<nq;i++)if(qs[i].kind=='p')
        point_queries[k++]=(interval_point_t){qs[i].v[0],qs[i].v[1],qs[i].v[2]};
    for(size_t i=0;i<nr;i++)if(!(rs[i]=alea_raycast_result_create()))return 1;
    interval_line_t* ray_queries=calloc(nr,sizeof(*ray_queries));
    if(!ray_queries)return 1;
    for(size_t i=0,k=0;i<nq;i++)if(qs[i].kind=='r')
        {
            const double* v=qs[i].v;
            ray_queries[k++]=(interval_line_t){{v[0],v[1],v[2]}, {v[3],v[4],v[5]}, v[6]};
        }
    start=alea_monotonic_seconds();
    interval_owner_t first;
    if (alea_prepared_locate_points(index,point_queries,1,&first)) return 1;
    double first_s=alea_monotonic_seconds()-start;
    alea_perf_reset();
    alea_hier_debug_candidates_reset();
    start=alea_monotonic_seconds();
    if(alea_prepared_locate_points(index,point_queries,np,ps))return 1;
    double point_s=alea_monotonic_seconds()-start;
    interval_memory("first point batch");
    size_t point_candidates = alea_hier_debug_candidates_get();
    size_t point_primitives = alea_perf_get().primitive_evaluations;
    start=alea_monotonic_seconds();
    if(alea_prepared_trace_rays(index,ray_queries,nr,rs))return 1;
    double ray_s=alea_monotonic_seconds()-start;
    interval_memory("first ray batch");
    size_t visits = 0, candidates = 0, intersections = 0;
    for (size_t i = 0; i < nr; i++) {
        visits += rs[i]->blas_universe_queries;
        candidates += rs[i]->blas_cell_candidates;
        intersections += rs[i]->surfaces_tested;
    }
    unsigned long long peak_kb=0;
    FILE* status=fopen("/proc/self/status","r");
    if(!status)return 1;
    while(fgets(line,sizeof(line),status))
        if(sscanf(line,"VmHWM: %llu kB",&peak_kb)==1)break;
    fclose(status);
    putchar('{');
    printf("\"ray_index_s\":%.9f,\"ray_index_bytes\":%zu,", idx_s, interval_index_bytes(index));
    const char* algorithm=interval_index_is_regional(index)?
        "parent-interval-regional":"parent-interval-full-expression";
    printf("\"point_work_scope\":\"original evaluator and BVH fallback only\",");
    interval_fallback_stats_t fallback=interval_index_fallback_stats(index);
    printf("\"ray_fallback_work\":{\"visits\":%zu,\"candidates\":%zu,\"intersections\":%zu},",
           fallback.visits,fallback.candidates,fallback.intersections);
    free(ray_queries);
    alea_prepared_queries_destroy(index);
    printf("\"ray_algorithm\":\"%s\","
           "\"ray_work\":{\"universe_visits\":%zu,\"candidates\":%zu,"
           "\"intersections\":%zu},", algorithm, visits, candidates, intersections);
    printf("\"schema\":1,\"engine\":\"libalea\",\"engine_hash\":\"%s\","
        "\"contract\":\"prepared native batches returning the selected deepest owner\","
        "\"threads\":1,\"mode\":\"%s\",\"cells\":%zu,\"load_s\":%.17g,\"ready_s\":%.17g,"
        "\"write_s\":%.17g,\"first_query_s\":%.17g,\"point_s\":%.17g,\"ray_s\":%.17g,"
        "\"regional_s\":0,\"regional_roots\":0,\"regional_bytes\":0,"
        "\"grid_s\":%.17g,\"grid_tiles\":%zu,\"grid_entries\":%zu,\"grid_bytes\":%zu,"
        "\"grid_node_visits\":%zu,\"point_candidates\":%zu,\"point_primitives\":%zu,"
        "\"peak_rss_bytes\":%llu,\"points\":[",
        alea_snapshot_engine_id(),argv[1],model->sys->cells.count,load_s,ready_s,write_s,first_s,point_s,ray_s,
        grid_s,grid_stats.tiles,grid_stats.candidate_entries,grid_stats.bytes,
        grid_stats.node_visits,point_candidates,point_primitives,peak_kb*1024);
    for(size_t i=0;i<np;i++)printf("%s{\"status\":\"selected\",\"cell\":%d,\"material\":%d,\"rc\":%d}",
        i?",":"",ps[i].cell_id,ps[i].material_id,ps[i].status);
    printf("],\"rays\":[");
    for(size_t i=0;i<nr;i++) {
        printf("%s[",i?",":"");
        for(size_t k=0;k<alea_raycast_segment_count(rs[i]);k++) {
            double lo,hi,density;int cell,mat,s0,s1;uint8_t flags;
            if(alea_raycast_segment_get(rs[i],k,&lo,&hi,&cell,&mat,&density,&s0,&s1) ||
               alea_raycast_segment_resolution_flags(rs[i],k,&flags))return 1;
            printf("%s{\"enter\":%.17g,\"exit\":%.17g,\"label\":{\"status\":\"selected\",\"cell\":%d,"
                "\"material\":%d,\"density\":%.17g,\"flags\":%u},\"surface_in\":%d,\"surface_out\":%d}",
                k?",":"",lo,hi,cell,mat,density,flags,s0,s1);
        }
        printf("]");alea_raycast_result_destroy(rs[i]);
    }
    puts("]}");free(point_queries);free(qs);free(ps);free(rs);mcnp_model_destroy(model);return 0;
}
