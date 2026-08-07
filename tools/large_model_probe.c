// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

#include "alea.h"
#include "alea_mcnp.h"
#include "core/alea_system.h"
#include "core/alea_spatial_hier.h"
#include "core/alea_universe.h"
#include "core/alea_eval.h"
#include "raycast/raycast.h"
#include "primitives/bbox.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#if defined(_WIN32)
#  include <windows.h>
#  include <psapi.h>
#else
#  include <sys/resource.h>
#  include <sys/time.h>
#endif

#define POINT_SAMPLE_THRESHOLD 16
#define POINT_BVH_THRESHOLD 8192

/* Short names for the primitive-type histogram, indexed by alea_primitive_type_t. */
static const char* primitive_type_name(int type) {
    switch (type) {
        case ALEA_PRIMITIVE_PLANE:      return "plane";
        case ALEA_PRIMITIVE_SPHERE:     return "sphere";
        case ALEA_PRIMITIVE_CYLINDER_X: return "cyl_x";
        case ALEA_PRIMITIVE_CYLINDER_Y: return "cyl_y";
        case ALEA_PRIMITIVE_CYLINDER_Z: return "cyl_z";
        case ALEA_PRIMITIVE_CONE_X:     return "cone_x";
        case ALEA_PRIMITIVE_CONE_Y:     return "cone_y";
        case ALEA_PRIMITIVE_CONE_Z:     return "cone_z";
        case ALEA_PRIMITIVE_RPP:        return "rpp";
        case ALEA_PRIMITIVE_QUADRIC:    return "quadric";
        case ALEA_PRIMITIVE_TORUS_X:    return "torus_x";
        case ALEA_PRIMITIVE_TORUS_Y:    return "torus_y";
        case ALEA_PRIMITIVE_TORUS_Z:    return "torus_z";
        case ALEA_PRIMITIVE_RCC:        return "rcc";
        case ALEA_PRIMITIVE_BOX:        return "box";
        case ALEA_PRIMITIVE_SPH:        return "sph";
        case ALEA_PRIMITIVE_TRC:        return "trc";
        case ALEA_PRIMITIVE_ELL:        return "ell";
        case ALEA_PRIMITIVE_REC:        return "rec";
        case ALEA_PRIMITIVE_WED:        return "wed";
        case ALEA_PRIMITIVE_RHP:        return "rhp";
        case ALEA_PRIMITIVE_ARB:        return "arb";
        default:                        return "?";
    }
}

/* Print the Phase 2 per-ray surface-test instrumentation for a hier raycast. */
static void print_surface_instrumentation(const alea_raycast_result_t* result) {
    double avg_surf = result->crossed_cell_count
        ? (double)result->sum_cell_surface_count / (double)result->crossed_cell_count
        : 0.0;
    printf("  surf_tests: terminal=%d lattice=%d ancestor=%d\n",
           result->terminal_surfaces_tested,
           result->lattice_surfaces_tested,
           result->ancestor_surfaces_tested);
    printf("  crossed_cells=%zu cell_surfaces avg=%.1f max=%u\n",
           result->crossed_cell_count, avg_surf,
           result->max_cell_surface_count);
    printf("  prim_type_tests:");
    for (int t = 0; t < ALEA_RAYCAST_PRIM_TYPE_BINS; t++) {
        if (result->prim_type_tests[t])
            printf(" %s=%u", primitive_type_name(t), result->prim_type_tests[t]);
    }
    printf("\n");
}

static size_t point_bvh_threshold(void) {
    const char* env = getenv("ALEA_UNIVERSE_POINT_BVH_THRESHOLD");
    if (env && env[0]) {
        char* end = NULL;
        unsigned long value = strtoul(env, &end, 10);
        if (end != env && value > 0) return (size_t)value;
    }
    return POINT_BVH_THRESHOLD;
}

static double now_seconds(void) {
#if defined(_WIN32)
    LARGE_INTEGER freq, counter;
    QueryPerformanceFrequency(&freq);
    QueryPerformanceCounter(&counter);
    return (double)counter.QuadPart / (double)freq.QuadPart;
#else
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return (double)tv.tv_sec + (double)tv.tv_usec * 1e-6;
#endif
}

static long peak_rss_kib(void) {
#if defined(_WIN32)
    PROCESS_MEMORY_COUNTERS pmc;
    if (!GetProcessMemoryInfo(GetCurrentProcess(), &pmc, sizeof(pmc))) return -1;
    return (long)(pmc.PeakWorkingSetSize / 1024);
#else
    struct rusage ru;
    if (getrusage(RUSAGE_SELF, &ru) != 0) return -1;
#  if defined(__APPLE__)
    /* macOS reports ru_maxrss in bytes; Linux/BSD in KiB. */
    return ru.ru_maxrss / 1024;
#  else
    return ru.ru_maxrss;
#  endif
#endif
}

static alea_bbox_t cell_bbox(alea_system_t* sys, size_t cell_index) {
    const alea_cell_entry_t* cell = &sys->cells.data[cell_index];
    if (cell->root_node_id == ALEA_NODE_ID_INVALID) {
        return alea_bbox_empty();
    }
    const alea_node_t* root = &sys->nodes.data[cell->root_node_id];
    if (root->bbox.min_x <= root->bbox.max_x) {
        return alea_node_bbox_get(&root->bbox);
    }
    return alea_get_bbox(sys, cell->root_node_id);
}

static int universe_has_lattice(const alea_system_t* sys, const alea_universe_t* univ,
                                size_t* out_lattice_cells) {
    size_t n = 0;
    for (size_t i = 0; i < univ->cell_indices.count; i++) {
        const alea_cell_entry_t* cell = &sys->cells.data[univ->cell_indices.data[i]];
        if (cell->lat_type != 0 && cell->lat_fill) n++;
    }
    if (out_lattice_cells) *out_lattice_cells = n;
    return n > 0;
}

static void print_id_stats(const alea_system_t* sys) {
    int min_cell = 0, max_cell = 0;
    int min_surf = 0, max_surf = 0;

    for (size_t i = 0; i < sys->cells.count; i++) {
        int id = sys->cells.data[i].mc_cell_id;
        if (i == 0 || id < min_cell) min_cell = id;
        if (i == 0 || id > max_cell) max_cell = id;
    }
    for (size_t i = 0; i < sys->surfaces.count; i++) {
        int id = sys->surfaces.data[i].mc_surface_id;
        if (i == 0 || id < min_surf) min_surf = id;
        if (i == 0 || id > max_surf) max_surf = id;
    }

    double cell_density = max_cell > 0 ? (double)sys->cells.count / (double)max_cell : 0.0;
    double surf_density = max_surf > 0 ? (double)sys->surfaces.count / (double)max_surf : 0.0;
    double dense_cell_mib = max_cell > 0 ? ((double)(max_cell + 1) * sizeof(int)) / 1048576.0 : 0.0;
    double dense_surf_mib = max_surf > 0 ? ((double)(max_surf + 1) * sizeof(uint32_t)) / 1048576.0 : 0.0;

    printf("ID stats:\n");
    printf("  cells:    count=%zu min=%d max=%d density=%.6f dense-int-map=%.1f MiB\n",
           sys->cells.count, min_cell, max_cell, cell_density, dense_cell_mib);
    printf("  surfaces: count=%zu min=%d max=%d density=%.6f dense-u32-map=%.1f MiB\n",
           sys->surfaces.count, min_surf, max_surf, surf_density, dense_surf_mib);
}

static void print_universe_stats(const alea_system_t* sys) {
    size_t large = 0, large_lattice = 0, large_non_lattice = 0;
    size_t max_cells = 0, max_idx = 0;
    size_t total_indexed_cells = 0;
    size_t threshold = point_bvh_threshold();

    for (size_t i = 0; i < sys->universes.count; i++) {
        const alea_universe_t* univ = &sys->universes.data[i];
        size_t lattice_cells = 0;
        int has_lattice = universe_has_lattice(sys, univ, &lattice_cells);
        if (univ->cell_indices.count > max_cells) {
            max_cells = univ->cell_indices.count;
            max_idx = i;
        }
        if (univ->cell_indices.count > threshold) {
            large++;
            if (has_lattice) large_lattice++;
            else {
                large_non_lattice++;
                total_indexed_cells += univ->cell_indices.count;
            }
        }
    }

    printf("Universe stats:\n");
    printf("  universes=%zu threshold=%zu\n", sys->universes.count, threshold);
    printf("  large=%zu large_non_lattice=%zu large_lattice_skipped=%zu\n",
           large, large_non_lattice, large_lattice);
    printf("  current-BLAS-eligible cells=%zu\n", total_indexed_cells);
    if (sys->universes.count > 0) {
        const alea_universe_t* u = &sys->universes.data[max_idx];
        size_t lattice_cells = 0;
        universe_has_lattice(sys, u, &lattice_cells);
        printf("  largest universe: id=%d cells=%zu lattice_cells=%zu\n",
               u->universe_id, u->cell_indices.count, lattice_cells);
    }
}

static void run_center_queries(alea_system_t* sys, size_t max_queries) {
    size_t done = 0;
    alea_universe_point_bvh_stats_reset();
    double t0 = now_seconds();
    size_t total_hits = 0;
    size_t errors = 0;

    for (size_t ui = 0; ui < sys->universes.count && done < max_queries; ui++) {
        const alea_universe_t* univ = &sys->universes.data[ui];
        size_t lattice_cells = 0;
        if (univ->cell_indices.count <= POINT_SAMPLE_THRESHOLD ||
            universe_has_lattice(sys, univ, &lattice_cells)) {
            continue;
        }

        for (size_t i = 0; i < univ->cell_indices.count && done < max_queries; i++) {
            size_t cell_idx = univ->cell_indices.data[i];
            alea_bbox_t bbox = cell_bbox(sys, cell_idx);
            if (!alea_bbox_is_valid(&bbox)) continue;
            double x = 0.5 * (bbox.min_x + bbox.max_x);
            double y = 0.5 * (bbox.min_y + bbox.max_y);
            double z = 0.5 * (bbox.min_z + bbox.max_z);
            alea_cell_hit_t hits[64];
            int n = alea_find_all_cells_at_point_recursive(sys, x, y, z, hits, 64);
            if (n < 0) errors++;
            else total_hits += (size_t)n;
            done++;
        }
    }

    double t1 = now_seconds();
    alea_universe_point_bvh_stats_t stats;
    alea_universe_point_bvh_stats_get(&stats);
    printf("Recursive center-query probe:\n");
    printf("  queries=%zu errors=%zu total_hits=%zu time=%.3f s avg=%.3f ms/query\n",
           done, errors, total_hits, t1 - t0,
           done ? (t1 - t0) * 1000.0 / (double)done : 0.0);
    printf("  bvh: builds=%zu queries=%zu node_visits=%zu bbox_tests=%zu candidates=%zu\n",
           stats.bvh_builds, stats.bvh_queries, stats.bvh_node_visits,
           stats.bvh_bbox_tests, stats.bvh_candidates);
    printf("  linear: universe_scans=%zu cell_tests=%zu exact_cell_tests=%zu\n",
           stats.linear_universe_scans, stats.linear_cell_tests,
           stats.exact_cell_tests);
    printf("  averages: bvh_candidates/query=%.2f exact_tests/query=%.2f linear_cell_tests/query=%.2f\n",
           stats.bvh_queries ? (double)stats.bvh_candidates / (double)stats.bvh_queries : 0.0,
           done ? (double)stats.exact_cell_tests / (double)done : 0.0,
           done ? (double)stats.linear_cell_tests / (double)done : 0.0);
}

typedef struct {
    double ms;
    size_t prim;
    size_t bool_ops;
    size_t candidates;
    size_t universe_id;
    size_t cell_idx;
    double x, y, z;
} hier_query_sample_t;

static int sample_cmp_desc(const void* a, const void* b) {
    double da = ((const hier_query_sample_t*)a)->ms;
    double db = ((const hier_query_sample_t*)b)->ms;
    return (da < db) - (da > db);
}

static void run_hier_center_queries(alea_system_t* sys, size_t max_queries) {
    if (alea_hier_spatial_index_build(sys) != 0) {
        fprintf(stderr, "hierarchical spatial build failed: %s\n", alea_error());
        return;
    }

    hier_query_sample_t* samples = calloc(max_queries, sizeof(*samples));
    if (!samples) { fprintf(stderr, "alloc failed\n"); return; }

    size_t done = 0;
    size_t total_hits = 0;
    size_t errors = 0;
    double t0 = now_seconds();

    for (size_t ui = 0; ui < sys->universes.count && done < max_queries; ui++) {
        const alea_universe_t* univ = &sys->universes.data[ui];
        size_t lattice_cells = 0;
        if (univ->cell_indices.count <= POINT_SAMPLE_THRESHOLD ||
            universe_has_lattice(sys, univ, &lattice_cells)) {
            continue;
        }

        for (size_t i = 0; i < univ->cell_indices.count && done < max_queries; i++) {
            size_t cell_idx = univ->cell_indices.data[i];
            alea_bbox_t bbox = cell_bbox(sys, cell_idx);
            if (!alea_bbox_is_valid(&bbox)) continue;
            double x = 0.5 * (bbox.min_x + bbox.max_x);
            double y = 0.5 * (bbox.min_y + bbox.max_y);
            double z = 0.5 * (bbox.min_z + bbox.max_z);

            alea_perf_reset();
            alea_hier_debug_candidates_reset();
            double q0 = now_seconds();
            int n;
            alea_cell_hit_t hits[64];
            if (getenv("PROBE_SINGLE_CELL")) {
                int rc = alea_hier_spatial_find_cell_in_universe(sys, 0, x, y, z);
                hits[0].cell_index = rc;
                n = (rc >= 0) ? 1 : 0;
            } else {
                n = alea_hier_spatial_find_cells_at_point(sys, x, y, z, hits, 64);
            }
            double q1 = now_seconds();
            alea_perf_counters_t pc = alea_perf_get();

            samples[done].ms = (q1 - q0) * 1000.0;
            samples[done].prim = pc.primitive_evaluations;
            samples[done].bool_ops = pc.boolean_operations;
            samples[done].candidates = alea_hier_debug_candidates_get();
            samples[done].universe_id = ui;
            samples[done].cell_idx = cell_idx;
            samples[done].x = x; samples[done].y = y; samples[done].z = z;

            if (n < 0) errors++;
            else total_hits += (size_t)n;
            done++;
        }
    }

    double t1 = now_seconds();
    printf("Hierarchical center-query probe:\n");
    printf("  queries=%zu errors=%zu total_hits=%zu time=%.3f s avg=%.3f ms/query\n",
           done, errors, total_hits, t1 - t0,
           done ? (t1 - t0) * 1000.0 / (double)done : 0.0);
    /* Per-query distribution */
    if (done > 0) {
        double total_ms = 0;
        size_t total_prim = 0, total_bool = 0, total_cands = 0;
        for (size_t i = 0; i < done; i++) {
            total_ms += samples[i].ms;
            total_prim += samples[i].prim;
            total_bool += samples[i].bool_ops;
            total_cands += samples[i].candidates;
        }
        qsort(samples, done, sizeof(*samples), sample_cmp_desc);

        /* Buckets: <0.1, <1, <10, <100, >=100 ms */
        size_t b[5] = {0};
        for (size_t i = 0; i < done; i++) {
            double m = samples[i].ms;
            if      (m < 0.1)   b[0]++;
            else if (m < 1.0)   b[1]++;
            else if (m < 10.0)  b[2]++;
            else if (m < 100.0) b[3]++;
            else                b[4]++;
        }

        printf("  histogram (ms/query):\n");
        printf("    <0.1:  %zu (%.1f%%)\n", b[0], 100.0*b[0]/done);
        printf("    <1:    %zu (%.1f%%)\n", b[1], 100.0*b[1]/done);
        printf("    <10:   %zu (%.1f%%)\n", b[2], 100.0*b[2]/done);
        printf("    <100:  %zu (%.1f%%)\n", b[3], 100.0*b[3]/done);
        printf("    >=100: %zu (%.1f%%)\n", b[4], 100.0*b[4]/done);

        printf("  totals: prim_evals=%zu bool_ops=%zu csg_ops/query=%.0f cands/query=%.0f\n",
               total_prim, total_bool, (double)(total_prim + total_bool) / done,
               (double)total_cands / done);

        /* Diagnose universe 0 bbox quality at the origin */
        const alea_universe_t* u0 = alea_get_universe(sys, 0);
        if (u0) {
            size_t cells_in_u0 = u0->cell_indices.count;
            size_t containing = 0;
            for (size_t i = 0; i < cells_in_u0; i++) {
                size_t ci = u0->cell_indices.data[i];
                const alea_cell_entry_t* c = &sys->cells.data[ci];
                if (c->root_node_id == ALEA_NODE_ID_INVALID) continue;
                const alea_bbox_t bb_v = alea_node_bbox_get(&sys->nodes.data[c->root_node_id].bbox);
                const alea_bbox_t* bb = &bb_v;
                if (bb->min_x <= 0 && 0 <= bb->max_x &&
                    bb->min_y <= 0 && 0 <= bb->max_y &&
                    bb->min_z <= 0 && 0 <= bb->max_z) {
                    containing++;
                }
            }
            printf("  universe 0: cells=%zu  bboxes_containing_origin=%zu (%.1f%%)\n",
                   cells_in_u0, containing,
                   cells_in_u0 ? 100.0 * containing / cells_in_u0 : 0.0);
        }

        printf("  per-cand: csg_ops=%.0f prim=%.0f bool=%.0f\n",
               total_cands ? (double)(total_prim + total_bool) / total_cands : 0,
               total_cands ? (double)total_prim / total_cands : 0,
               total_cands ? (double)total_bool / total_cands : 0);

        /* Top 5 slowest */
        size_t show = done < 5 ? done : 5;
        printf("  top-%zu slowest:\n", show);
        for (size_t i = 0; i < show; i++) {
            const alea_universe_t* u = &sys->universes.data[samples[i].universe_id];
            const alea_cell_entry_t* c = &sys->cells.data[samples[i].cell_idx];
            printf("    %.2f ms  univ=%d  cell=%d  cands=%zu prim=%zu bool=%zu  pt=(%.1f,%.1f,%.1f)\n",
                   samples[i].ms, u->universe_id, c->mc_cell_id,
                   samples[i].candidates, samples[i].prim, samples[i].bool_ops,
                   samples[i].x, samples[i].y, samples[i].z);
        }
    }
    free(samples);
}

static int run_hier_slice_query(alea_system_t* sys,
                                double z,
                                double x_min,
                                double x_max,
                                double y_min,
                                double y_max,
                                size_t max_hits) {
    if (alea_hier_spatial_index_build(sys) != 0) {
        fprintf(stderr, "hierarchical spatial build failed: %s\n", alea_error());
        return -1;
    }

    alea_spatial_hit_t* hits = calloc(max_hits ? max_hits : 1, sizeof(*hits));
    if (!hits) {
        fprintf(stderr, "failed to allocate slice hit buffer\n");
        return -1;
    }

    double t0 = now_seconds();
    int n = alea_hier_spatial_query_slice_z(sys, z, x_min, x_max, y_min, y_max,
                                            hits, max_hits);
    double t1 = now_seconds();
    free(hits);

    if (n < 0) {
        fprintf(stderr, "hierarchical slice query failed: %s\n", alea_error());
        return -1;
    }

    printf("Hierarchical slice probe:\n");
    printf("  z=%.6g x=[%.6g,%.6g] y=[%.6g,%.6g] max_hits=%zu hits=%d time=%.3f s\n",
           z, x_min, x_max, y_min, y_max, max_hits, n, t1 - t0);
    return 0;
}

static int run_hier_raycast(alea_system_t* sys,
                            double ox, double oy, double oz,
                            double dx, double dy, double dz,
                            double t_max,
                            int fast_segments,
                            int blas_experimental,
                            size_t repeat_count) {
    alea_raycast_result_t result;
    alea_raycast_result_init(&result);

    double t0 = now_seconds();
    int rc = blas_experimental
        ? alea_raycast_hier_blas_experimental(sys, ox, oy, oz, dx, dy, dz,
                                              t_max, &result)
        : (fast_segments
            ? alea_raycast_hier_fast_segments(sys, ox, oy, oz, dx, dy, dz,
                                              t_max, &result)
            : alea_raycast_hier(sys, ox, oy, oz, dx, dy, dz, t_max,
                                &result));
    double t1 = now_seconds();

    if (rc != 0) {
        alea_raycast_result_free(&result);
        fprintf(stderr, "hierarchical raycast failed: %s\n", alea_error());
        return -1;
    }

    printf("Hierarchical %s%sraycast probe:\n",
           blas_experimental ? "BLAS experimental " : "",
           fast_segments ? "fast-segments " : "");
    if (fast_segments) {
        printf("  note=segments-first; carries fill/lattice path state and does not reconstruct full hit lists\n");
    }
    if (blas_experimental) {
        printf("  note=experimental; TLAS/BLAS-pruned surface hit generation\n");
    }
    printf("  origin=(%.6g,%.6g,%.6g) dir=(%.6g,%.6g,%.6g) t_max=%.6g\n",
           ox, oy, oz, dx, dy, dz, t_max);
    printf("  hits=%zu segments=%zu surfaces_tested=%d point_lookups=%d steps=%d time=%.3f s\n",
           result.hits.count, result.segments.count, result.surfaces_tested,
           result.point_lookups, result.step_iterations, t1 - t0);
    print_surface_instrumentation(&result);
    if (blas_experimental) {
        printf("  blas_counts: placements=%zu pruned=%zu universe_queries=%zu cell_candidates=%zu cells_tested=%zu hits_pre_dedup=%zu\n",
               result.blas_placement_candidates,
               result.blas_placements_pruned,
               result.blas_universe_queries,
               result.blas_cell_candidates,
               result.blas_cells_tested,
               result.blas_hits_before_dedup);
    }
    if (repeat_count > 1) {
        double warm_total = 0.0;
        double warm_min = 1e300;
        double warm_max = 0.0;
        size_t warm_runs = repeat_count - 1;
        for (size_t i = 0; i < warm_runs; i++) {
            alea_raycast_result_clear(&result);
            double tw0 = now_seconds();
            rc = blas_experimental
                ? alea_raycast_hier_blas_experimental(sys, ox, oy, oz,
                                                       dx, dy, dz, t_max,
                                                       &result)
                : (fast_segments
                    ? alea_raycast_hier_fast_segments(sys, ox, oy, oz,
                                                      dx, dy, dz, t_max,
                                                      &result)
                    : alea_raycast_hier(sys, ox, oy, oz, dx, dy, dz,
                                        t_max, &result));
            double tw1 = now_seconds();
            if (rc != 0) {
                alea_raycast_result_free(&result);
                fprintf(stderr, "hierarchical repeated raycast failed: %s\n",
                        alea_error());
                return -1;
            }
            double dt = tw1 - tw0;
            warm_total += dt;
            if (dt < warm_min) warm_min = dt;
            if (dt > warm_max) warm_max = dt;
        }
        printf("  repeat=%zu warm_avg=%.6f s warm_min=%.6f s warm_max=%.6f s\n",
               repeat_count, warm_total / warm_runs, warm_min, warm_max);
        printf("  warm_last: hits=%zu segments=%zu surfaces_tested=%d point_lookups=%d steps=%d\n",
               result.hits.count, result.segments.count,
               result.surfaces_tested, result.point_lookups,
               result.step_iterations);
        print_surface_instrumentation(&result);
        if (blas_experimental) {
            printf("  warm_last_blas_counts: placements=%zu pruned=%zu universe_queries=%zu cell_candidates=%zu cells_tested=%zu hits_pre_dedup=%zu\n",
                   result.blas_placement_candidates,
                   result.blas_placements_pruned,
                   result.blas_universe_queries,
                   result.blas_cell_candidates,
                   result.blas_cells_tested,
                   result.blas_hits_before_dedup);
        }
    }
    if (result.segments.count > 0) {
        const alea_ray_segment_t* first = &result.segments.data[0];
        const alea_ray_segment_t* last =
            &result.segments.data[result.segments.count - 1];
        printf("  first_segment: t=[%.6g,%.6g] cell=%d material=%d\n",
               first->t_enter, first->t_exit, first->cell_id,
               first->material_id);
        printf("  last_segment:  t=[%.6g,%.6g] cell=%d material=%d\n",
               last->t_enter, last->t_exit, last->cell_id,
               last->material_id);
    }

    alea_raycast_result_free(&result);
    return 0;
}

int main(int argc, char** argv) {
    if (argc < 2) {
        fprintf(stderr, "usage: %s MODEL.mcnp [--queries N] [--hier-build] [--hier-queries N] [--hier-slice Z XMIN XMAX YMIN YMAX MAX_HITS] [--hier-raycast OX OY OZ DX DY DZ TMAX] [--hier-fast-segments OX OY OZ DX DY DZ TMAX] [--hier-cell-raycast OX OY OZ DX DY DZ TMAX] [--hier-blas-raycast OX OY OZ DX DY DZ TMAX] [--raycast-repeat N]\n", argv[0]);
        fprintf(stderr, "set ALEA_HIER_BLAS_THRESHOLD=N to skip BLAS builds for universes below N cells\n");
        return 2;
    }

    size_t queries = 0;
    size_t hier_queries = 0;
    int hier_slice = 0;
    double slice_z = 0.0, slice_x_min = 0.0, slice_x_max = 0.0;
    double slice_y_min = 0.0, slice_y_max = 0.0;
    size_t slice_max_hits = 0;
    int hier_build = 0;
    int hier_raycast = 0;
    int hier_fast_segments = 0;
    int hier_blas_raycast = 0;
    double ray_ox = 0.0, ray_oy = 0.0, ray_oz = 0.0;
    double ray_dx = 0.0, ray_dy = 0.0, ray_dz = 1.0, ray_t_max = 0.0;
    size_t raycast_repeat = 1;
    for (int i = 2; i < argc; i++) {
        if (strcmp(argv[i], "--queries") == 0 && i + 1 < argc) {
            queries = (size_t)strtoull(argv[++i], NULL, 10);
        } else if (strcmp(argv[i], "--hier-build") == 0) {
            hier_build = 1;
        } else if (strcmp(argv[i], "--hier-queries") == 0 && i + 1 < argc) {
            hier_queries = (size_t)strtoull(argv[++i], NULL, 10);
        } else if (strcmp(argv[i], "--hier-slice") == 0 && i + 6 < argc) {
            hier_slice = 1;
            slice_z = strtod(argv[++i], NULL);
            slice_x_min = strtod(argv[++i], NULL);
            slice_x_max = strtod(argv[++i], NULL);
            slice_y_min = strtod(argv[++i], NULL);
            slice_y_max = strtod(argv[++i], NULL);
            slice_max_hits = (size_t)strtoull(argv[++i], NULL, 10);
        } else if (strcmp(argv[i], "--hier-raycast") == 0 && i + 7 < argc) {
            hier_raycast = 1;
            ray_ox = strtod(argv[++i], NULL);
            ray_oy = strtod(argv[++i], NULL);
            ray_oz = strtod(argv[++i], NULL);
            ray_dx = strtod(argv[++i], NULL);
            ray_dy = strtod(argv[++i], NULL);
            ray_dz = strtod(argv[++i], NULL);
            ray_t_max = strtod(argv[++i], NULL);
        } else if ((strcmp(argv[i], "--hier-fast-segments") == 0 ||
                    strcmp(argv[i], "--hier-cell-raycast") == 0) &&
                   i + 7 < argc) {
            hier_fast_segments = 1;
            ray_ox = strtod(argv[++i], NULL);
            ray_oy = strtod(argv[++i], NULL);
            ray_oz = strtod(argv[++i], NULL);
            ray_dx = strtod(argv[++i], NULL);
            ray_dy = strtod(argv[++i], NULL);
            ray_dz = strtod(argv[++i], NULL);
            ray_t_max = strtod(argv[++i], NULL);
        } else if (strcmp(argv[i], "--hier-blas-raycast") == 0 && i + 7 < argc) {
            hier_blas_raycast = 1;
            ray_ox = strtod(argv[++i], NULL);
            ray_oy = strtod(argv[++i], NULL);
            ray_oz = strtod(argv[++i], NULL);
            ray_dx = strtod(argv[++i], NULL);
            ray_dy = strtod(argv[++i], NULL);
            ray_dz = strtod(argv[++i], NULL);
            ray_t_max = strtod(argv[++i], NULL);
        } else if (strcmp(argv[i], "--raycast-repeat") == 0 && i + 1 < argc) {
            raycast_repeat = (size_t)strtoull(argv[++i], NULL, 10);
            if (raycast_repeat == 0) raycast_repeat = 1;
        }
    }

    alea_log_set_level(1);

    double t0 = now_seconds();
    mcnp_model_t* model = mcnp_load(argv[1]);
    double t1 = now_seconds();
    if (!model) {
        fprintf(stderr, "load failed: %s\n", alea_error());
        return 1;
    }

    alea_system_t* sys = model->sys;
    printf("Loaded %s\n", argv[1]);
    printf("  load_time=%.3f s peak_rss=%ld KiB\n", t1 - t0, peak_rss_kib());
    print_id_stats(sys);

    double tu0 = now_seconds();
    if (alea_build_universe_index(sys) != 0) {
        fprintf(stderr, "universe index failed: %s\n", alea_error());
        mcnp_model_destroy(model);
        return 1;
    }
    double tu1 = now_seconds();
    printf("Universe index:\n");
    printf("  build_time=%.3f s peak_rss=%ld KiB\n", tu1 - tu0, peak_rss_kib());
    print_universe_stats(sys);

    if (hier_build) {
        double th0 = now_seconds();
        if (alea_hier_spatial_index_build(sys) != 0) {
            fprintf(stderr, "hierarchical spatial build failed: %s\n", alea_error());
            mcnp_model_destroy(model);
            return 1;
        }
        double th1 = now_seconds();
        const alea_hier_spatial_stats_t* stats =
            alea_hier_spatial_index_stats(sys->hier_spatial_index);
        printf("Hierarchical spatial BLAS probe:\n");
        printf("  build_time=%.3f s peak_rss=%ld KiB\n", th1 - th0, peak_rss_kib());
        if (stats) {
            printf("  universes=%zu blas=%zu linear=%zu largest_universe=%d cells=%d\n",
                   stats->universe_count, stats->blas_count,
                   stats->linear_universe_count, stats->largest_universe_id,
                   stats->max_universe_cells);
            printf("  blas_cells=%zu blas_nodes=%zu memory_estimate=%.1f MiB\n",
                   stats->blas_cell_count, stats->blas_node_count,
                   (double)stats->memory_bytes / 1048576.0);
            printf("  placements=%zu root=%zu fill=%zu lattice=%zu max_depth=%d transforms=%zu\n",
                   stats->placement_count, stats->root_placement_count,
                   stats->fill_placement_count, stats->lattice_placement_count,
                   stats->max_placement_depth, stats->transform_count);
            printf("  placement_sources: fill_cells=%zu lattice_cells=%zu\n",
                   stats->fill_cell_count, stats->lattice_cell_count);
        }
    }

    if (queries > 0) {
        run_center_queries(sys, queries);
        printf("  peak_rss_after_queries=%ld KiB\n", peak_rss_kib());
    }

    if (hier_queries > 0) {
        run_hier_center_queries(sys, hier_queries);
        printf("  peak_rss_after_hier_queries=%ld KiB\n", peak_rss_kib());
    }

    if (hier_slice) {
        if (run_hier_slice_query(sys, slice_z, slice_x_min, slice_x_max,
                                 slice_y_min, slice_y_max, slice_max_hits) != 0) {
            mcnp_model_destroy(model);
            return 1;
        }
        printf("  peak_rss_after_hier_slice=%ld KiB\n", peak_rss_kib());
    }

    if (hier_raycast || hier_fast_segments || hier_blas_raycast) {
        if (run_hier_raycast(sys, ray_ox, ray_oy, ray_oz,
                             ray_dx, ray_dy, ray_dz, ray_t_max,
                             hier_fast_segments, hier_blas_raycast,
                             raycast_repeat) != 0) {
            mcnp_model_destroy(model);
            return 1;
        }
        printf("  peak_rss_after_hier_%s%sraycast=%ld KiB\n",
               hier_blas_raycast ? "blas_" : "",
               hier_fast_segments ? "fast_segments_" : "", peak_rss_kib());
    }

    mcnp_model_destroy(model);
    return 0;
}
