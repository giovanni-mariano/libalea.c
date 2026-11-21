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

#include "interval_overlap.h"

typedef struct {char kind; double v[7];} query;

typedef struct {
    double distance;
    size_t ordinal;
} replay_row_t;

static int by_distance(const void* lhs, const void* rhs) {
    const replay_row_t* a = lhs;
    const replay_row_t* b = rhs;
    if (a->distance != b->distance) return a->distance < b->distance ? -1 : 1;
    return (a->ordinal > b->ordinal) - (a->ordinal < b->ordinal);
}

static size_t matching_segments(const alea_raycast_batch_result_t* batch, size_t offset,
                                size_t count, const alea_raycast_result_t* candidate,
                                bool compare_surfaces, bool compare_flags, double tolerance) {
    size_t m = alea_raycast_segment_count(candidate);
    size_t j = 0;
    for (; j < count && j < m; ++j) {
        double lo, hi, density;
        int cell, material, enter_surface, exit_surface;
        size_t k = offset + j;
        if (alea_raycast_segment_get(candidate, j, &lo, &hi, &cell, &material,
                                     &density, &enter_surface, &exit_surface) ||
            cell != alea_raycast_batch_cell_ids(batch)[k] ||
            material != alea_raycast_batch_material_ids(batch)[k] ||
            density != alea_raycast_batch_densities(batch)[k] ||
            (compare_flags && candidate->segments.data[j].resolution_flags !=
                              alea_raycast_batch_resolution_flags(batch)[k]) ||
            !(fabs(lo - alea_raycast_batch_t_enter(batch)[k]) <= tolerance) ||
            !(fabs(hi - alea_raycast_batch_t_exit(batch)[k]) <= tolerance) ||
            (compare_surfaces && (enter_surface != alea_raycast_batch_enter_surface_ids(batch)[k] ||
                                  exit_surface != alea_raycast_batch_exit_surface_ids(batch)[k]))) break;
    }
    return j;
}

/* Walk the common partition from the first failed segment. A long interval
 * can contain an extra crossing on the other path, so comparing the two
 * entries at the same ordinal would report the wrong boundary. */
static void report_boundary(alea_system_t* sys, size_t ray,
                            const alea_raycast_batch_result_t* batch, size_t begin, size_t end,
                            const alea_raycast_result_t* candidate, size_t match) {
    const double* enter = alea_raycast_batch_t_enter(batch);
    const double* exit = alea_raycast_batch_t_exit(batch);
    const int32_t* cells = alea_raycast_batch_cell_ids(batch);
    size_t m = candidate->segments.count;
    double start = INFINITY;
    if (match < end - begin) start = enter[begin + match];
    if (match < m) start = fmin(start, candidate->segments.data[match].t_enter);
    size_t i = begin, j = 0;
    while (i < end && j < m) {
        const alea_ray_segment_t* seg = &candidate->segments.data[j];
        double lo = fmax(start, fmax(enter[i], seg->t_enter));
        double hi = fmin(exit[i], seg->t_exit);
        if (lo < hi && cells[i] != seg->cell_id) {
            fprintf(stderr, "boundary ray %zu public %d %d %a %d %a prepared %d %d %a %d %a overlap %a %a\n",
                    ray, cells[i], alea_raycast_batch_enter_surface_ids(batch)[i], enter[i],
                    alea_raycast_batch_exit_surface_ids(batch)[i], exit[i],
                    seg->cell_id, seg->enter_surface_id, seg->t_enter,
                    seg->exit_surface_id, seg->t_exit, lo, hi);
            const double probes[] = {nextafter(lo, hi), lo + (hi-lo)*0.5, nextafter(hi, lo),
                                     nextafter(lo, -INFINITY), nextafter(hi, INFINITY)};
            for (size_t p = 0; p < 5; ++p) {
                double t = probes[p];
                if (p < 3 && !(lo < t && t < hi)) continue;
                if (t < 0 || t > candidate->segments.data[m-1].t_exit) continue;
                double x, y, z;
                int cell = -1, material = -1;
                alea_ray_point_at(&candidate->ray, t, &x, &y, &z);
                int status = alea_find_cell_at(sys, x, y, z, &cell, &material);
                fprintf(stderr, "probe ray %zu sample %zu t %a xyz %a %a %a status %d cell %d material %d\n",
                        ray, p, t, x, y, z, status, cell, material);
            }
            /* Exhaust short intervals in parameter space. Several adjacent
             * parameters may still round to the same world coordinate. */
            size_t samples = 0, public_owner = 0, prepared_owner = 0, missed = 0;
            double t = nextafter(lo, hi);
            while (t < hi && samples < 65536) {
                double x, y, z;
                int cell;
                alea_ray_point_at(&candidate->ray, t, &x, &y, &z);
                if (alea_find_cell_at(sys, x, y, z, &cell, NULL)) ++missed;
                else if (cell == cells[i]) ++public_owner;
                else if (cell == seg->cell_id) ++prepared_owner;
                ++samples;
                t = nextafter(t, hi);
            }
            fprintf(stderr, "interior ray %zu count %zu public %zu prepared %zu other %zu missed %zu complete %d\n",
                    ray, samples, public_owner, prepared_owner,
                    samples - public_owner - prepared_owner - missed, missed, t == hi);
            return;
        }
        double public_end = exit[i];
        if (public_end <= seg->t_exit) ++i;
        if (seg->t_exit <= public_end) ++j;
    }
    fprintf(stderr, "boundary ray %zu no overlapping ownership difference after segment %zu\n", ray, match);
}

/* Preserve the geometry failure. This reports the ordered physical cuts even
 * when the two rounded distances in the diagnostic are identical. */
static int report_substep_fills(alea_system_t* sys, size_t ordinal, const interval_line_t* ray,
                               const alea_raycast_result_t* result) {
    for (size_t j = 0; j < result->segments.count; ++j) {
        const alea_ray_segment_t* s = &result->segments.data[j];
        if (!(s->resolution_flags & ALEA_RESOLVE_UNDEFINED_FILL) ||
            s->t_exit > nextafter(s->t_enter, INFINITY)) continue;
        alea_ray_substep_fill_t gap;
        int rc = alea_ray_diagnose_substep_fill(sys,ray->origin,ray->direction,s->t_enter,&gap);
        if (rc < 0) return -1;
        if (rc && gap.parent_cell_id == s->cell_id)
            fprintf(stderr, "ray %zu substep undefined fill, parent %d universe %d, "
                    "surfaces %d -> %d, rounded distances %a %a, width %.17g cm\n",
                    ordinal,gap.parent_cell_id,gap.fill_universe,gap.enter_surface_id,
                    gap.exit_surface_id,gap.t_enter,gap.t_exit,gap.width);
    }
    return 0;
}

/* Public batches take a single distance. Keep original ray ordinals while
 * grouping equal limits, then compare each CSR row with its prepared result. */
static int compare_public_rays(alea_system_t* sys, interval_index_t* index, const interval_line_t* rays,
                               alea_raycast_result_t** prepared, size_t count) {
    const char* mode = getenv("ALEA_INTERVAL_COMPARE_PUBLIC");
    bool shared = mode && !strcmp(mode, "shared");
    if (shared) {
        if (alea_set_ray_batch_accelerator(sys, index)) return -1;
    } else if (alea_prepare_query_acceleration(sys)) return -1;
    alea_raycast_batch_options_t options = {
        .struct_size = sizeof(options),
        .fields = ALEA_RAY_BATCH_MATERIAL | ALEA_RAY_BATCH_DENSITY
    };
    alea_raycast_batch_result_t* batch = alea_raycast_batch_result_create();
    alea_raycast_result_t scratch;
    alea_raycast_result_init(&scratch);
    replay_row_t* rows = calloc(count, sizeof(*rows));
    double (*origins)[3] = calloc(count, sizeof(*origins));
    double (*directions)[3] = calloc(count, sizeof(*directions));
    int rc = -1;
    if (!batch || !rows || !origins || !directions) goto done;
    for (size_t i = 0; i < count; ++i) rows[i] = (replay_row_t){rays[i].distance, i};
    qsort(rows, count, sizeof(*rows), by_distance);
    bool boundaries = mode && strcmp(mode, "boundaries") == 0;
    bool compare_surfaces = mode && strcmp(mode, "surfaces") == 0;
    bool score_overlaps = mode && (!strcmp(mode, "overlaps") || !strcmp(mode, "1"));
    bool details = boundaries || (mode && strcmp(mode, "details") == 0);
    if (boundaries || compare_surfaces) options.fields |= ALEA_RAY_BATCH_SURFACES;
    if (score_overlaps || shared) options.fields |= ALEA_RAY_BATCH_RESOLUTION_FLAGS;
    size_t mismatches = 0;
    size_t overlap_only = 0;
    for (size_t begin = 0, end; begin < count; begin = end) {
        end = begin;
        do {
            const interval_line_t* r = &rays[rows[end].ordinal];
            memcpy(origins[end-begin], r->origin, sizeof(*origins));
            memcpy(directions[end-begin], r->direction, sizeof(*directions));
            ++end;
        } while (end < count && rows[end].distance == rows[begin].distance);
        if (alea_raycast_hier_batch(sys, origins[0], directions[0], end - begin,
                                    rows[begin].distance, &options, batch)) goto done;
        if (shared) {
            alea_raycast_batch_work_stats_t stats;
            if (alea_raycast_batch_result_get_work_stats_internal(batch, &stats) ||
                stats.regional_ray_count != end - begin) goto done;
        }
        const uint64_t* offsets = alea_raycast_batch_ray_offsets(batch);
        for (size_t i = begin; i < end; ++i) {
            size_t row = i - begin, ordinal = rows[i].ordinal;
            size_t n = offsets[row + 1] - offsets[row];
            size_t m = alea_raycast_segment_count(prepared[ordinal]);
            size_t match = matching_segments(batch, offsets[row], n, prepared[ordinal],
                                             compare_surfaces, score_overlaps || shared, shared ? 0.0 : 1e-6);
            if (match == n && match == m) continue;
            if (score_overlaps && overlap_only_ray(sys, batch, offsets[row], offsets[row+1],
                                                   prepared[ordinal], rows[i].distance, &scratch)) {
                fprintf(stderr, "ray %zu differs only within verified overlaps\n", ordinal);
                ++overlap_only;
                continue;
            }
            if (boundaries)
                report_boundary(sys, ordinal, batch, offsets[row], offsets[row + 1], prepared[ordinal], match);
            if ((boundaries || score_overlaps) &&
                report_substep_fills(sys, ordinal, &rays[ordinal], prepared[ordinal])) goto done;
            if (details || score_overlaps || mismatches < 3)
                fprintf(stderr, "ray %zu, first differing segment %zu, public/prepared counts %zu/%zu\n",
                        ordinal, match, n, m);
            if (compare_surfaces && mismatches < 3 && match < n && match < m) {
                size_t k = offsets[row] + match;
                const alea_ray_segment_t* segment = &prepared[ordinal]->segments.data[match];
                fprintf(stderr, "  entry/exit surface IDs: public %d/%d, prepared %d/%d\n",
                        alea_raycast_batch_enter_surface_ids(batch)[k],
                        alea_raycast_batch_exit_surface_ids(batch)[k],
                        segment->enter_surface_id, segment->exit_surface_id);
            }
            if (details && match < n) {
                size_t k = offsets[row] + match;
                fprintf(stderr, "  public   cell %d material %d density %.17g [%.17g, %.17g]\n",
                        alea_raycast_batch_cell_ids(batch)[k], alea_raycast_batch_material_ids(batch)[k],
                        alea_raycast_batch_densities(batch)[k], alea_raycast_batch_t_enter(batch)[k],
                        alea_raycast_batch_t_exit(batch)[k]);
            }
            if (details && match < m) {
                double lo, hi, density;
                int cell, material;
                if (alea_raycast_segment_get(prepared[ordinal], match, &lo, &hi,
                                             &cell, &material, &density, NULL, NULL)) goto done;
                fprintf(stderr, "  prepared cell %d material %d density %.17g [%.17g, %.17g]\n",
                        cell, material, density, lo, hi);
            }
            ++mismatches;
        }
    }
    if (score_overlaps) {
        printf("\"geometry_comparison\":{\"rays\":%zu,\"matched\":%zu,\"overlap_only\":%zu,"
               "\"remaining\":%zu},", count, count - mismatches - overlap_only, overlap_only, mismatches);
        fprintf(stderr, "raw disagreements %zu, verified overlap only %zu, remaining %zu\n",
                mismatches + overlap_only, overlap_only, mismatches);
    }
    fprintf(stderr, "public/prepared ray mismatches %zu/%zu%s\n", mismatches, count,
            compare_surfaces ? " including surface IDs" : "");
    rc = mismatches != 0;
done:
    alea_raycast_result_free(&scratch);
    free(rows);
    free(origins);
    free(directions);
    alea_raycast_batch_result_destroy(batch);
    return rc;
}

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
    int comparison = 0;
    if (getenv("ALEA_INTERVAL_COMPARE_PUBLIC"))
        comparison = compare_public_rays(model->sys, index, ray_queries, rs, nr);
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
    puts("]}");free(point_queries);free(qs);free(ps);free(rs);mcnp_model_destroy(model);return comparison != 0;
}
