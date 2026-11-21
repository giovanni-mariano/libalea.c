// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
// SPDX-License-Identifier: MPL-2.0
/* One process per scenario: high-water RSS must not combine unrelated loads. */
#include "alea_mcnp.h"
#include "util/ray_profile.h"
#include "alea_raycast.h"
#include "alea_fingerprint.h"
#include "core/snapshot_internal.h"
#include "core/alea_universe.h"
#include "raycast/regional/interval_index.h"
#include "core/alea_spatial_hier.h"
#include "primitives/bbox.h"
#include "raycast/recovery_planes.h"
#include "raycast/ray_traversal.h"
#include "util/compat.h"
#include "util/alea_parallel.h"
#include "interval_memory.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <errno.h>
#include <limits.h>
#ifndef _WIN32
#include <sys/resource.h>
#endif

typedef struct { char kind; double v[7]; } query_t;
static int hierarchical_rays;
static bool study_ray_root;
typedef struct {
    const interval_index_t* index;
    size_t queries, hits, candidates;
} root_replay_t;

static int replay_root_candidate(void* context, double x, double y, double z, int surface) {
    root_replay_t* replay = context;
    ++replay->queries;
    int cell = alea_regional_ray_root_candidate(replay->index, x, y, z, surface, &replay->candidates);
    if (cell >= 0) ++replay->hits;
    return cell;
}

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

static int timed_point_query(alea_system_t* sys,const query_t* q,alea_sha256_t* hash, double* elapsed) {
    int cell=-1,material=-1;
    double begin = elapsed ? alea_monotonic_seconds() : 0;
    int rc=alea_find_cell_at(sys,q->v[0],q->v[1],q->v[2],&cell,&material);
    if (elapsed) *elapsed += alea_monotonic_seconds() - begin;
    alea_digest_u64(hash,(uint64_t)(int64_t)rc);
    alea_digest_u64(hash,(uint64_t)(int64_t)cell);
    alea_digest_u64(hash,(uint64_t)(int64_t)material);
    return cell>=0;
}

static int point_query(alea_system_t* sys,const query_t* q,alea_sha256_t* hash) {
    return timed_point_query(sys, q, hash, NULL);
}

typedef struct {
    alea_system_t* sys;
    const query_t* queries;
    alea_point_owner_t* owners;
    bool compare;
} point_replay_t;

static int replay_owner_range(void* context, size_t worker, size_t begin, size_t end) {
    (void)worker;
    point_replay_t* replay = context;
    for (size_t i = begin; i < end; ++i) {
        if (replay->queries[i].kind != 'p') continue;
        const double* v = replay->queries[i].v;
        alea_point_owner_t owner = {.cell_id = -1, .material_id = -1};
        owner.status = alea_find_cell_at(replay->sys, v[0], v[1], v[2],
                                         &owner.cell_id, &owner.material_id);
        if (replay->compare) {
            alea_point_owner_t ref = replay->owners[i];
            if (ref.status != owner.status || ref.cell_id != owner.cell_id ||
                ref.material_id != owner.material_id) return -1;
        } else {
            replay->owners[i] = owner;
        }
    }
    return 0;
}

static int measure_point_scaling(alea_system_t* sys, const query_t* queries, size_t count) {
    alea_point_owner_t* owners = calloc(count, sizeof(*owners));
    if (!owners) return -1;
    point_replay_t replay = {sys, queries, owners, false};
    printf(",\"point_scaling\":[");
    int rc = 0;
    for (size_t cap = 1; cap <= 8; cap *= 2) {
        size_t workers;
        double start = alea_monotonic_seconds();
        if (alea_parallel_for(count, 1, cap, ALEA_PARALLEL_STATIC_BLOCK,
                              replay_owner_range, &replay, &workers) != ALEA_PARALLEL_OK) {
            rc = -1;
            break;
        }
        printf("%s{\"workers\":%zu,\"seconds\":%.9g}", cap == 1 ? "" : ",",
               workers, alea_monotonic_seconds() - start);
        replay.compare = true;
    }
    putchar(']');
    free(owners);
    return rc;
}

static int ray_query(alea_system_t* sys,const query_t* q,
                     alea_raycast_result_t* result,alea_sha256_t* hash,size_t* segments,
                     const alea_ray_root_lookup_t* lookup) {
    const double* v=q->v;
    int rc;
    if (lookup)
        rc = alea_raycast_hier_segments_with_root_lookup(sys, v[0], v[1], v[2],
                    v[3], v[4], v[5], v[6], result, lookup);
    else rc=hierarchical_rays
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

/* Each ray runs once per method, separated by one whole pass. Alternating
 * the assignment within each pass spreads load and cache drift across both. */
static int compare_ray_roots(alea_system_t* sys, const query_t* queries, size_t count,
                             alea_raycast_result_t* result, const alea_ray_root_lookup_t* regional) {
    alea_digest_t* first = calloc(count, sizeof(*first));
    if (!first) return -1;
    double seconds[2][2] = {{0}};
    size_t segments[2] = {0};
    int rc = -1;
    for (size_t pass = 0; pass < 2; ++pass) {
        for (size_t i = 0; i < count; ++i) {
            size_t method = (i + pass) % 2;
            alea_sha256_t hash;
            alea_digest_t digest;
            alea_sha256_init(&hash);
            double start = alea_monotonic_seconds();
            if (ray_query(sys, &queries[i], result, &hash, &segments[method],
                          method ? regional : NULL)) goto done;
            seconds[pass][method] += alea_monotonic_seconds() - start;
            alea_sha256_final(&hash, digest.bytes);
            if (!pass) first[i] = digest;
            else if (memcmp(first[i].bytes, digest.bytes, sizeof(digest.bytes))) {
                fprintf(stderr, "Root lookup changed the segment digest for ray %zu\n", i);
                goto done;
            }
        }
    }
    printf("{\"rays_per_method\":%zu,\"segments_per_method\":%zu,\"hierarchy_s\":%.9g,\"regional_s\":%.9g,"
           "\"passes\":[{\"hierarchy_s\":%.9g,\"regional_s\":%.9g},{\"hierarchy_s\":%.9g,\"regional_s\":%.9g}]}\n",
           count, segments[0], seconds[0][0] + seconds[1][0], seconds[0][1] + seconds[1][1],
           seconds[0][0], seconds[0][1], seconds[1][0], seconds[1][1]);
    rc = 0;
done:
    free(first);
    return rc;
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

static int compare_universe_ids(const void* a, const void* b) {
    int left = *(const int*)a, right = *(const int*)b;
    return (left > right) - (left < right);
}

static int parse_child_ids(const char* list, alea_int_vec_t* ids) {
    const char* p = list;
    for (;;) {
        if (*p < '0' || *p > '9') return -1;
        char* end;
        errno = 0;
        long id = strtol(p, &end, 10);
        if (errno || id <= 0 || id > INT_MAX || alea_vec_push(ids, (int)id, int)) return -1;
        if (!*end) break;
        if (*end != ',') return -1;
        p = end + 1;
    }
    qsort(ids->data, ids->count, sizeof(*ids->data), compare_universe_ids);
    return 0;
}

typedef struct { double total_s, window_s; } child_sample_t;
enum { PROMOTION_WINDOW = 2048 };

typedef struct {
    alea_system_t* sys;
    const query_t* queries;
    size_t begin, count;
    int universe;
    alea_sha256_t* hash;
    double replay_start, decision_start, build_s, next_result_s;
    _Atomic size_t completed;
    size_t completions_during_build;
} promotion_tail_t;

/* Two TinyPar tasks share the read-only model. The reader alone updates the
 * ordered digest; both tasks finish before the caller releases the model. */
static int build_while_reading(void* context, size_t worker, size_t begin, size_t end) {
    (void)worker;
    promotion_tail_t* tail = context;
    for (size_t task = begin; task < end; ++task) {
        if (task == 0) {
            for (size_t i = tail->begin; i < tail->count; ++i) {
                point_query(tail->sys, &tail->queries[i], tail->hash);
                atomic_store_explicit(&tail->completed, i + 1 - tail->begin, memory_order_relaxed);
                if (i == tail->begin) tail->next_result_s = alea_monotonic_seconds() - tail->decision_start;
                if ((i + 1) % PROMOTION_WINDOW == 0 || i + 1 == tail->count)
                    printf(",%.9g", alea_monotonic_seconds() - tail->replay_start);
            }
        } else {
            size_t prior = atomic_load_explicit(&tail->completed, memory_order_relaxed);
            double start = alea_monotonic_seconds();
            if (alea_promote_point_child(tail->sys, tail->universe)) return -1;
            tail->build_s = alea_monotonic_seconds() - start;
            tail->completions_during_build = atomic_load_explicit(&tail->completed, memory_order_relaxed) - prior;
        }
    }
    return 0;
}

/* Single-reader experiment. Selection only sees completed query windows.
 * Background measurement includes dispatch and the final join. */
static int promotion_replay(alea_system_t* sys, const query_t* queries, size_t count,
                            long threshold_ms, bool child_only, bool background, double ready_s, double prepare_s) {
    const size_t window = PROMOTION_WINDOW;
    double threshold = threshold_ms * 1e-3;
    child_sample_t* samples = NULL;
    if (threshold_ms) {
        samples = calloc(sys->universes.count + 1, sizeof(*samples));
        if (!samples) return -1;
    }
    int promoted = 0;
    size_t after = 0;
    double evidence_s = 0, rebuild_s = 0;
    size_t during_build = 0;
    double decision_start = 0, next_result_s = 0;
    alea_sha256_t hash;
    alea_sha256_init(&hash);
    alea_selected_query_stats_t stats;
    alea_selected_query_stats(sys, &stats);
    const char* build_mode = "none";
    if (threshold_ms) build_mode = child_only ? "child" : "replacement";
    if (background) build_mode = "background-child";
    printf("{\"points\":%zu,\"threshold_ms\":%ld,\"promotion_build\":\"%s\","
           "\"ready_s\":%.9g,\"prepare_s\":%.9g,\"ready_index_bytes\":%zu,\"windows_s\":[",
           count, threshold_ms, build_mode, ready_s, prepare_s, stats.bytes);
    double start = alea_monotonic_seconds();
    for (size_t begin = 0; begin < count; begin += window) {
        size_t end = begin + window < count ? begin + window : count;
        size_t i = begin;
        if (after && after == begin) {
            point_query(sys, &queries[i++], &hash);
            next_result_s = alea_monotonic_seconds() - decision_start;
        }
        for (; i < end; ++i) point_query(sys, &queries[i], &hash);
        if (threshold_ms && !promoted && end - begin == window) {
            for (size_t i = 0; i < sys->universes.count; ++i) {
                uint64_t visits;
                double total;
                alea_root_point_scan_tests(sys, i, &visits, &total);
                double recent = total - samples[i].total_s;
                if (total >= threshold && recent >= threshold / 32 &&
                    samples[i].window_s >= threshold / 32 && total > evidence_s &&
                    alea_adaptive_child_progress(sys, i, NULL)) {
                    promoted = sys->universes.data[i].universe_id;
                    evidence_s = total;
                }
                samples[i] = (child_sample_t){total, recent};
            }
            if (promoted) {
                interval_point_children_t choice = {&promoted, 1, CHILD_POINTS_ADAPTIVE};
                after = end;
                decision_start = alea_monotonic_seconds();
                if (background) {
                    printf("%s%.9g", begin ? "," : "", decision_start - start);
                    promotion_tail_t tail = {.sys=sys, .queries=queries, .begin=end, .count=count,
                        .universe=promoted, .hash=&hash, .replay_start=start, .decision_start=decision_start};
                    atomic_init(&tail.completed, 0);
                    size_t workers;
                    if (alea_parallel_for(2, 1, 2, ALEA_PARALLEL_STATIC_BLOCK, build_while_reading,
                                          &tail, &workers) != ALEA_PARALLEL_OK || workers != 2) {
                        free(samples);
                        return -1;
                    }
                    rebuild_s = tail.build_s;
                    during_build = tail.completions_during_build;
                    next_result_s = tail.next_result_s;
                    break;
                }
                int rc;
                if (child_only) rc = alea_promote_point_child(sys, promoted);
                else {
                    alea_clear_selected_queries(sys);
                    rc = alea_prepare_root_points(sys, &choice);
                }
                if (rc) {
                    free(samples);
                    return -1;
                }
                rebuild_s = alea_monotonic_seconds() - decision_start;
            }
        }
        printf("%s%.9g", begin ? "," : "", alea_monotonic_seconds() - start);
    }
    double replay_s = alea_monotonic_seconds() - start;
    free(samples);
    if (sys->hier_spatial_index) return -1;
    alea_digest_t digest;
    char hex[65];
    alea_sha256_final(&hash, digest.bytes);
    alea_digest_hex(&digest, hex);
    alea_selected_query_stats(sys, &stats);
    printf("],\"point_digest\":\"%s\",\"point_s\":%.9g,\"promoted_universe\":%d,"
           "\"after_points\":%zu,\"local_s_at_promotion\":%.9g,\"rebuild_s\":%.9g,"
           "\"index_bytes\":%zu,\"completions_during_build\":%zu,\"next_result_s\":", hex, replay_s, promoted, after,
           evidence_s, rebuild_s, stats.bytes, during_build);
    if (after && after < count) printf("%.9g", next_result_s);
    else printf("null");
    printf(",\"peak_rss_bytes\":");
    peak_rss();
    puts("}");
    return 0;
}

/* Record geometry after the timers so an independent replay can verify
 * each discrete field and compare the physical intervals. */
static int print_public_ray_records(size_t count, alea_raycast_result_t** results,
                                    const alea_raycast_batch_result_t* batch) {
    printf(",\"rays\":[");
    const uint64_t* offsets = batch ? alea_raycast_batch_ray_offsets(batch) : NULL;
    for (size_t i = 0; i < count; ++i) {
        printf("%s[", i ? "," : "");
        size_t segments = results ? alea_raycast_segment_count(results[i]) : offsets[i+1]-offsets[i];
        for (size_t j = 0; j < segments; ++j) {
            int cell, material;
            double lo, hi, density;
            uint8_t flags;
            if (results) {
                if (alea_raycast_segment_get(results[i], j, &lo, &hi, &cell, &material, &density, NULL, NULL) ||
                    alea_raycast_segment_resolution_flags(results[i], j, &flags)) return -1;
            } else {
                size_t k = offsets[i]+j;
                lo = alea_raycast_batch_t_enter(batch)[k];
                hi = alea_raycast_batch_t_exit(batch)[k];
                cell = alea_raycast_batch_cell_ids(batch)[k];
                material = alea_raycast_batch_material_ids(batch)[k];
                density = alea_raycast_batch_densities(batch)[k];
                flags = alea_raycast_batch_resolution_flags(batch)[k];
            }
            printf("%s[%.17g,%.17g,%d,%d,%.17g,%u]", j ? "," : "", lo, hi, cell, material, density, (unsigned)flags);
        }
        printf("]");
    }
    printf("]");
    return 0;
}

/* Prototype for scalar callers whose eventual query count is unknown. Sample
 * one lookup per 16 calls, excluding caller time and digest construction.
 * The 64-call minimum keeps isolated queries out of index construction. */
typedef struct {
    size_t observed, promoted_after;
    double threshold_s, sampled_s, build_s;
    bool attempted;
} point_promotion_study_t;

static int study_point_promotion(alea_system_t* sys, const query_t* row,
                                 alea_sha256_t* hash, point_promotion_study_t* study) {
    if (!study->attempted && study->observed >= 64 &&
        study->sampled_s * 16 >= study->threshold_s) {
        study->attempted = true;
        double start = alea_monotonic_seconds();
        if (alea_prepare_point_acceleration(sys)) return -1;
        study->build_s = alea_monotonic_seconds() - start;
        if (alea_system_selected_index(sys)) study->promoted_after = study->observed;
    }
    if (study->attempted) {
        point_query(sys, row, hash);
    } else {
        ++study->observed;
        timed_point_query(sys, row, hash, study->observed % 16 == 0 ? &study->sampled_s : NULL);
    }
    return 0;
}

/* Benchmark normal entry points from a fresh source load. The hierarchy
 * ray batch builds its caches lazily. Each process measures one mode. */
static int public_query_study(int argc, char** argv) {
    const bool automatic = !strcmp(argv[4], "auto-points");
    const bool lookup = automatic || !strcmp(argv[4], "lookup") || !strcmp(argv[4], "hierarchy-lookup");
    const bool points = lookup || !strcmp(argv[4], "points") || !strcmp(argv[4], "adaptive");
    const bool prepared = !strcmp(argv[4], "prepared");
    const bool shared = !strcmp(argv[4], "shared-rays");
    const bool ready_rays = !strcmp(argv[4], "ready-rays");
    if (!points && !prepared && !shared && !ready_rays && strcmp(argv[4], "rays")) return 2;
    point_promotion_study_t study = {.threshold_s = 0.5};
    if (automatic && argc == 6) {
        char* end;
        long ms;
        errno = 0;
        ms = strtol(argv[5], &end, 10);
        if (argv[5][0] < '0' || argv[5][0] > '9' || *end || errno || ms > INT_MAX) return 2;
        study.threshold_s = ms / 1000.;
    } else if (argc != (prepared || shared ? 11 : 5)) return 2;
    alea_bbox_t domain = {0};
    if (prepared || shared) {
        double* limits[] = {&domain.min_x, &domain.max_x, &domain.min_y,
                            &domain.max_y, &domain.min_z, &domain.max_z};
        for (int k = 0; k < 6; ++k) {
            char* end;
            *limits[k] = strtod(argv[k + 5], &end);
            if (end == argv[k + 5] || *end || !isfinite(*limits[k])) return 2;
        }
    }
    size_t count;
    query_t* rows = read_workload(argv[3], &count);
    if (!rows) return 2;
    int status = 1;
    mcnp_model_t* model = NULL;
    alea_prepared_queries_t* index = NULL;
    alea_raycast_batch_result_t* batch = NULL;
    alea_raycast_result_t** results = NULL;
    alea_segment_ray_t* rays = NULL;
    double* origins = NULL;
    double* directions = NULL;
    for (size_t i = 0; i < count; ++i) {
        if (rows[i].kind != (points ? 'p' : 'r') || (!points && rows[i].v[6] != rows[0].v[6])) {
            fprintf(stderr, "Public replay requires points only or rays of one shared length\n");
            goto done;
        }
    }
    if (prepared) {
        rays = calloc(count, sizeof(*rays));
        results = calloc(count, sizeof(*results));
        if (!rays || !results) goto done;
        for (size_t i = 0; i < count; ++i) {
            memcpy(rays[i].origin, rows[i].v, 3 * sizeof(double));
            memcpy(rays[i].direction, rows[i].v + 3, 3 * sizeof(double));
            rays[i].distance = rows[i].v[6];
            results[i] = alea_raycast_result_create();
            if (!results[i]) goto done;
        }
    } else if (!points) {
        origins = calloc(count, 3 * sizeof(double));
        directions = calloc(count, 3 * sizeof(double));
        batch = alea_raycast_batch_result_create();
        if (!origins || !directions || !batch) goto done;
        for (size_t i = 0; i < count; ++i) {
            memcpy(origins + 3*i, rows[i].v, 3 * sizeof(double));
            memcpy(directions + 3*i, rows[i].v + 3, 3 * sizeof(double));
        }
    }

    double start = alea_monotonic_seconds();
    model = mcnp_load(argv[2]);
    if (!model) goto done;
    double loaded = alea_monotonic_seconds();
    if (lookup && alea_prepare_point_lookup(model->sys)) goto done;
    if ((!strcmp(argv[4], "points") || ready_rays) && alea_prepare_query_acceleration(model->sys)) goto done;
    if (!strcmp(argv[4], "adaptive") && alea_prepare_point_acceleration(model->sys)) goto done;
    /* Keep the independent trigger experiment and fixed hierarchy controls isolated. */
    if (automatic || !strcmp(argv[4], "points") || !strcmp(argv[4], "hierarchy-lookup"))
        atomic_store(&model->sys->point_promotion_attempted, 1);
    if (shared && alea_prepare_ray_acceleration(model->sys, &domain)) goto done;
    if (prepared && !(index = alea_prepared_queries_create(model->sys, &domain))) goto done;
    double ready = alea_monotonic_seconds();
    size_t sizes[] = {1, count < 16 ? count : 16, count};
    const alea_raycast_batch_options_t opts = {.struct_size = sizeof(opts),
        .fields = ALEA_RAY_BATCH_MATERIAL | ALEA_RAY_BATCH_DENSITY | ALEA_RAY_BATCH_RESOLUTION_FLAGS};
    printf("{\"load_s\":%.9g,\"prepare_s\":%.9g,\"batches\":[", loaded-start, ready-loaded);
    for (int pass = 0; pass < 3; ++pass) {
        size_t n = sizes[pass];
        alea_sha256_t hash;
        alea_sha256_init(&hash);
        double begin = alea_monotonic_seconds();
        if (automatic) {
            for (size_t i = 0; i < n; ++i)
                if (study_point_promotion(model->sys, rows+i, &hash, &study)) goto done;
        } else if (points) {
            for (size_t i = 0; i < n; ++i) point_query(model->sys, rows+i, &hash);
        } else if (prepared) {
            if (alea_prepared_trace_rays(index, rays, n, results)) goto done;
        } else if (alea_raycast_hier_batch(model->sys, origins, directions, n, rows[0].v[6], &opts, batch)) goto done;
        double end = alea_monotonic_seconds();
        printf("%s{\"n\":%zu,\"query_s\":%.9g,\"source_elapsed_s\":%.9g,\"peak_rss_bytes\":",
               pass ? "," : "", n, end-begin, end-start);
        peak_rss();
        if (points) {
            alea_digest_t digest;
            char hex[65];
            alea_sha256_final(&hash, digest.bytes);
            alea_digest_hex(&digest, hex);
            printf(",\"digest\":\"%s\"", hex);
        }
        printf("}");
    }
    printf("]");
    if (automatic) printf(",\"promoted_after\":%zu,\"policy_build_s\":%.9g,\"estimated_lookup_s\":%.9g",
                          study.promoted_after, study.build_s, study.sampled_s * 16);
    if (lookup && !automatic) printf(",\"regional_active\":%s",
        alea_selected_query_stats(model->sys, NULL) ? "true" : "false");
    if (!points && print_public_ray_records(count, results, batch)) goto done;
    puts("}");
    status = 0;
done:
    alea_prepared_queries_destroy(index);
    if (results) for (size_t i = 0; i < count; ++i) alea_raycast_result_destroy(results[i]);
    alea_raycast_batch_result_destroy(batch);
    if (model) mcnp_model_destroy(model);
    free(rows); free(rays); free(results); free(origins); free(directions);
    return status;
}

#include "navigator_replay.h"

int main(int argc,char** argv) {
    if (argc >= 5 && (!strcmp(argv[1], "navigator-script") || !strcmp(argv[1], "navigator"))) {
        bool regional = argc > 5 && !strcmp(argv[5], "regional");
        if (argc != 5 && !(regional && (argc == 6 || argc == 12))) return 2;
        alea_bbox_t box;
        if (argc == 12) {
            double* bounds[] = {&box.min_x, &box.max_x, &box.min_y, &box.max_y, &box.min_z, &box.max_z};
            for (int i = 0; i < 6; ++i) {
                char* end;
                *bounds[i] = strtod(argv[6+i], &end);
                if (*end || end == argv[6+i] || !isfinite(*bounds[i])) return 2;
            }
        }
        const alea_bbox_t* domain = argc == 12 ? &box : NULL;
        if (!strcmp(argv[1], "navigator-script"))
            return navigator_script(argv[2], argv[3], argv[4], regional, domain);
        return navigator_replay(argv[2], argv[3], argv[4], regional, domain);
    }
    if (argc >= 5 && !strcmp(argv[1], "public")) return public_query_study(argc, argv);
    if(argc<5) {
        fprintf(stderr,"Usage: model_lifecycle MODE MODEL WORKLOAD AUX [--mixed-bounds] [--selected-points|--root-points [[--adaptive-child-bounds] [--child-regions IDS]|--child-bounds IDS]|--hierarchy-points] [--point-scaling] [--repeat-points] [--hierarchical-rays|--study-ray-root|--compare-ray-root] [--candidate-profile] [--ray-profile]\n"
                "Navigator study: navigator MODEL WORKLOAD strict|interval|fast [regional [xmin xmax ymin ymax zmin zmax]]\n"
                "Particle operations: navigator-script MODEL SCRIPT strict|interval|fast [regional [xmin xmax ymin ymax zmin zmax]]\n"
                "Public API study: public MODEL WORKLOAD points|lookup|hierarchy-lookup|auto-points [ms]|adaptive|rays|ready-rays|prepared|shared-rays [xmin xmax ymin ymax zmin zmax]\n"
                "MODE: fresh, cache, snapshot, write, edit; AUX: cache/snapshot path (fresh: -)\n"
                "Point-only promotion study: fresh MODEL WORKLOAD - --root-points --study-promotion-ms N [--study-child-only|--study-background-child] (0 = control)\n");
        return 2;
    }
    int use_mixed=0,profile_candidates=0,profile_rays=0;
    bool selected_points = false, point_scaling = false, repeat_points = false;
    bool root_points = false, hierarchy_points = false;
    const char* child_list = NULL;
    bool adaptive = false, child_regions = false;
    long promotion_ms = -1;
    bool study_child_only = false;
    bool background = false, compare_roots = false;
    for(int i=5;i<argc;i++) {
        if(!strcmp(argv[i],"--mixed-bounds"))use_mixed=1;
        else if(!strcmp(argv[i],"--repeat-points"))repeat_points=true;
        else if(!strcmp(argv[i],"--point-scaling"))point_scaling=true;
        else if(!strcmp(argv[i],"--selected-points"))selected_points=true;
        else if(!strcmp(argv[i],"--root-points"))root_points=true;
        else if (!strcmp(argv[i], "--adaptive-child-bounds")) adaptive = true;
        else if (!strcmp(argv[i], "--study-child-only")) study_child_only = true;
        else if (!strcmp(argv[i], "--study-background-child")) { background = true; study_child_only = true; }
        else if (!strcmp(argv[i], "--study-promotion-ms")) {
            if (++i == argc || promotion_ms >= 0 || argv[i][0] < '0' || argv[i][0] > '9') return 2;
            char* end;
            errno = 0;
            promotion_ms = strtol(argv[i], &end, 10);
            if (*end || errno || promotion_ms > INT_MAX) return 2;
            adaptive = true;
        }
        else if (!strcmp(argv[i], "--child-regions") || !strcmp(argv[i], "--child-bounds")) {
            child_regions = !strcmp(argv[i], "--child-regions");
            if (child_list || ++i == argc) return 2;
            child_list = argv[i];
        }
        else if(!strcmp(argv[i],"--hierarchy-points"))hierarchy_points=true;
        else if (!strcmp(argv[i], "--compare-ray-root")) {
            compare_roots = true;
            hierarchical_rays = 1;
            study_ray_root = true;
        }
        else if (!strcmp(argv[i], "--study-ray-root")) { study_ray_root = true; hierarchical_rays = 1; }
        else if(!strcmp(argv[i],"--hierarchical-rays"))hierarchical_rays=1;
        else if(!strcmp(argv[i],"--candidate-profile"))profile_candidates=1;
        else if(!strcmp(argv[i],"--ray-profile"))profile_rays=1;
        else return 2;
    }
    if (compare_roots && argc != 6) return 2;
    if (adaptive && child_list && !child_regions) return 2;
    if (study_child_only && promotion_ms <= 0) return 2;
    if (background && alea_parallel_effective_workers(2, 1, 2) < 2) {
        fprintf(stderr, "The background promotion study needs two workers\n");
        return 2;
    }
    if (promotion_ms >= 0) {
        if (child_list || repeat_points || point_scaling || use_mixed || profile_candidates ||
            profile_rays || hierarchical_rays) return 2;
        const char* timing = getenv("ALEA_PROFILE_CHILD_POINTS");
        if (promotion_ms && (!timing || !*timing || !strcmp(timing, "0"))) {
            fprintf(stderr, "Promotion measurement requires ALEA_PROFILE_CHILD_POINTS=1\n");
            return 2;
        }
    }
    if ((adaptive || child_list) && !root_points) return 2;
    if (repeat_points && ((!root_points && !hierarchy_points) || point_scaling)) return 2;
#ifndef ALEA_RAY_PROFILE
    if(profile_rays) { fprintf(stderr,"--ray-profile requires an ALEA_RAY_PROFILE build\n");return 2; }
#endif
    const char* mode=argv[1];
    if(strcmp(mode,"fresh") && strcmp(mode,"cache") && strcmp(mode,"snapshot") &&
       strcmp(mode,"write") && strcmp(mode,"edit"))return 2;
    size_t count=0;
    query_t* queries=read_workload(argv[3],&count);
    if(!queries){fprintf(stderr,"Invalid workload: %s\n",argv[3]);return 2;}
    if (study_ray_root) {
        if (strcmp(mode, "fresh") || root_points || selected_points || hierarchy_points) {
            free(queries);
            return 2;
        }
        for (size_t i = 0; i < count; ++i) if (queries[i].kind != 'r') {
            free(queries);
            return 2;
        }
    }
    if (root_points || hierarchy_points) {
        if (strcmp(mode,"fresh") || (root_points && hierarchy_points) || selected_points) {
            free(queries);
            return 2;
        }
        for (size_t i=0;i<count;i++) if (queries[i].kind != 'p') {
            free(queries);
            return 2;
        }
    }
    alea_int_vec_t ids = {0};
    if (child_list && parse_child_ids(child_list, &ids)) {
        fprintf(stderr, "Expected comma-separated positive child universe IDs\n");
        alea_vec_free(&ids);
        free(queries);
        return 2;
    }
    interval_child_mode_t child_mode = adaptive ? CHILD_POINTS_ADAPTIVE :
                                       child_regions ? CHILD_POINTS_REGIONS : CHILD_POINTS_BOUNDS;
    interval_point_children_t children = {ids.data, ids.count, child_mode};
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
    size_t stage_count = root_points ? 0 : hierarchy_points ? 3 : 7;
    double stages[7];
    for(size_t k=0;k<stage_count;k++) {
        double t=alea_monotonic_seconds();
        if(alea_system_prepare_query_caches(model->sys,flags[k]))goto done;
        stages[k]=alea_monotonic_seconds()-t;
    }
    if (selected_points && alea_prepare_selected_queries(model->sys, NULL)) goto done;
    if (root_points && alea_prepare_root_points(model->sys, &children)) goto done;
    if (study_ray_root && alea_prepare_root_points(model->sys, NULL)) goto done;
    root_replay_t root_replay = {.index = alea_system_selected_index(model->sys)};
    alea_ray_root_lookup_t root_lookup = {replay_root_candidate, &root_replay};
    const alea_ray_root_lookup_t* lookup = study_ray_root ? &root_lookup : NULL;
    double ready=alea_monotonic_seconds();
    if (promotion_ms >= 0) {
        status = promotion_replay(model->sys, queries, count, promotion_ms, study_child_only, background,
                                  ready-start, ready-bounded) != 0;
        goto done;
    }
    alea_selected_query_stats_t ready_index = {0};
    if (root_points) alea_selected_query_stats(model->sys, &ready_index);
    interval_memory("ready");
    fprintf(stderr,"[alea-load-profile] ready %.9f s\n",ready-start);
    if(!strcmp(mode,"write") && mcnp_snapshot_write(model,argv[4]))goto done;
    double written=alea_monotonic_seconds();
    result=alea_raycast_result_create();if(!result)goto done;
    if (compare_roots) {
        status = compare_ray_roots(model->sys, queries, count, result, lookup) != 0;
        goto done;
    }
    alea_sha256_t points,rays,warmup;alea_digest_t pd,rd;
    alea_sha256_init(&warmup);
    size_t warm_segments=0;
    double first_start=alea_monotonic_seconds();
    if (hierarchy_points) atomic_store(&model->sys->point_promotion_attempted, 1);
    if(queries[0].kind=='p')point_query(model->sys,&queries[0],&warmup);
    else if(ray_query(model->sys,&queries[0],result,&warmup,&warm_segments,lookup))goto done;
    root_replay.queries = root_replay.hits = root_replay.candidates = 0;
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
    interval_memory("scalar points complete");
    fprintf(stderr,"[alea-load-profile] point_replay %.9f s\n",point_s);
    alea_universe_point_bvh_stats_t point_stats;
    alea_universe_point_bvh_stats_get(&point_stats);
    size_t point_candidates=alea_hier_debug_candidates_get();
    alea_hier_debug_candidates_reset();
    alea_hier_debug_candidate_cells_set(model->sys,ray_cells,model->sys->cells.count);
#ifdef ALEA_RAY_PROFILE
    if(profile_rays && alea_ray_profile_start(model->sys->universes.count))goto done;
#endif
    query_start=alea_monotonic_seconds();
    for(size_t i=0;i<count;i++)if(queries[i].kind=='r') {
#ifdef ALEA_RAY_PROFILE
        alea_ray_profile.ray_ordinal = nr;
#endif
        if(ray_query(model->sys,&queries[i],result,&rays,&segments,lookup))goto done;
        nr++;
    }
    double ray_s=alea_monotonic_seconds()-query_start;
#ifdef ALEA_RAY_PROFILE
    if(profile_rays)alea_ray_profile_stop();
#endif
    alea_hier_debug_candidate_cells_set(NULL,NULL,0);
    alea_sha256_final(&points,pd.bytes);alea_sha256_final(&rays,rd.bytes);
    char ph[65],rh[65],ah[65];alea_digest_hex(&pd,ph);alea_digest_hex(&rd,rh);
    if (root_points && model->sys->hier_spatial_index) {
        fprintf(stderr, "Root-only replay unexpectedly built the hierarchy\n");
        goto done;
    }
    if (!root_points && !hierarchy_points) acceleration_digest(model->sys,ah);
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
    peak_rss();
    if (point_scaling) {
        if (measure_point_scaling(model->sys, queries, count)) goto done;
        interval_memory("concurrent points complete");
        printf(",\"point_scaling_peak_rss_bytes\":");
        peak_rss();
    }
    if (root_points) {
        printf(",\"point_index_bytes_at_ready\":%zu", ready_index.bytes);
        if (ids.count) {
            const char* kind = child_regions ? "regions" : "bounds";
            if (adaptive) kind = "regions-and-adaptive-bounds";
            printf(",\"child_index\":\"%s\",\"child_universes\":[", kind);
            for (size_t i = 0; i < ids.count; ++i) printf("%s%d", i ? "," : "", ids.data[i]);
            putchar(']');
        } else if (adaptive) printf(",\"child_index\":\"adaptive-bounds\"");
        uint64_t total = 0;
        bool first = true;
        printf(",\"root_scans_by_universe\":[");
        for (size_t i = 0; i < model->sys->universes.count; ++i) {
            uint64_t visits;
            double seconds;
            uint64_t tests = alea_root_point_scan_tests(model->sys, i, &visits, &seconds);
            if (!visits && !seconds) continue;
            printf("%s{\"universe\":%d,\"visits\":%llu,\"cell_tests\":%llu",
                   first ? "" : ",", model->sys->universes.data[i].universe_id,
                   (unsigned long long)visits, (unsigned long long)tests);
            if (adaptive) {
                double build_seconds;
                size_t cells = alea_adaptive_child_progress(model->sys, i, &build_seconds);
                if (cells) printf(",\"bounded_cells\":%zu,\"build_s\":%.9g", cells, build_seconds);
            }
            if (seconds) printf(",\"local_point_s\":%.9g", seconds);
            putchar('}');
            first = false;
            total += tests;
        }
        printf("],\"root_scan_cell_tests\":%llu", (unsigned long long)total);
    }
    if (study_ray_root) {
        printf(",\"regional_root_queries\":%zu,\"regional_root_hits\":%zu,\"regional_root_candidates\":%zu",
               root_replay.queries, root_replay.hits, root_replay.candidates);
    }
    if (selected_points || root_points || study_ray_root) {
        alea_selected_query_stats_t stats;
        int active = alea_selected_query_stats(model->sys, &stats);
        alea_bbox_t b = stats.domain;
        printf(",\"selected_points\":{\"active\":%d,\"bytes\":%zu,\"build_s\":%.9g,"
               "\"fallback_visits\":%llu,\"domain\":[%.17g,%.17g,%.17g,%.17g,%.17g,%.17g]}",
               active, stats.bytes, stats.build_seconds, (unsigned long long)stats.fallback_visits,
               b.min_x, b.max_x, b.min_y, b.max_y, b.min_z, b.max_z);
    }
    if (repeat_points) {
        alea_digest_t repeated;
        alea_sha256_init(&points);
        double begin = alea_monotonic_seconds();
        for (size_t i = 0; i < count; ++i)
            point_query(model->sys, &queries[i], &points);
        double seconds = alea_monotonic_seconds() - begin;
        alea_sha256_final(&points, repeated.bytes);
        if (memcmp(repeated.bytes, pd.bytes, sizeof(pd.bytes)) ||
            (root_points && model->sys->hier_spatial_index)) {
            fprintf(stderr, "Repeated point replay changed ownership or built the hierarchy\n");
            goto done;
        }
        printf(",\"repeat_point_s\":%.9g", seconds);
        if (root_points) {
            alea_selected_query_stats_t stats;
            alea_selected_query_stats(model->sys, &stats);
            printf(",\"repeat_point_index_bytes\":%zu", stats.bytes);
        }
    }
    if (root_points || hierarchy_points)
        printf(",\"point_backend\":\"%s\"", (ids.count || adaptive) ? "regional-root-child" : root_points ? "regional-root" : "hierarchy");
    else
        printf(",\"ray_api\":\"%s\",\"acceleration_digest\":\"%s\"",
               hierarchical_rays ? "hierarchical-segments" : "global", ah);
    printf(",\"preparation_stages_s\":{");
    for(size_t k=0;k<stage_count;k++)printf("%s\"%s\":%.9g",k?",":"",names[k],stages[k]);
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
        printf("],\"universe_costs\":[");
        bool first_universe = true;
        for (size_t i = 0; i < model->sys->universes.count; ++i) {
            const ray_profile_universe* u = &alea_ray_profile.universes[i];
            if (!u->lookups && !u->recoveries) continue;
            printf("%s{\"universe\":%d,\"lookups\":%llu,\"lookup_self_s\":%.9g,"
                   "\"membership_calls\":%llu,\"membership_s\":%.9g,"
                   "\"enumerations\":%llu,\"candidates\":%llu,\"enumeration_s\":%.9g,"
                   "\"boundary_retries\":%llu,\"boundary_retry_s\":%.9g}",
                   first_universe ? "" : ",", model->sys->universes.data[i].universe_id,
                   (unsigned long long)u->lookups, u->lookup_seconds,
                   (unsigned long long)u->membership_calls, u->membership_seconds,
                   (unsigned long long)u->enumerations, (unsigned long long)u->candidates,
                   u->enumeration_seconds, (unsigned long long)u->recoveries, u->recovery_seconds);
            first_universe = false;
        }
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
               "\"cached_empty_windows\":%llu,\"peak_root_bytes\":%zu,\"root_build_types\":[",
               (unsigned long long)alea_ray_profile.indexed_root_lists,
               (unsigned long long)alea_ray_profile.roots_sorted,
               (unsigned long long)alea_ray_profile.prepared_card_checks_reused,
               (unsigned long long)alea_ray_profile.cached_root_candidates,
               (unsigned long long)alea_ray_profile.cached_empty_windows, alea_ray_profile.peak_root_bytes);
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
    mcnp_model_destroy(model);alea_vec_free(&ids);free(queries);return status;
}
