// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

#include "alea.h"
#include "alea_mcnp.h"
#include "core/alea_system.h"
#include "core/alea_spatial_hier.h"
#include "core/alea_universe.h"
#include "primitives/bbox.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/resource.h>
#include <sys/time.h>

#define POINT_SAMPLE_THRESHOLD 16
#define POINT_BVH_THRESHOLD 8192

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
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return (double)tv.tv_sec + (double)tv.tv_usec * 1e-6;
}

static long peak_rss_kib(void) {
    struct rusage ru;
    if (getrusage(RUSAGE_SELF, &ru) != 0) return -1;
    return ru.ru_maxrss;
}

static alea_bbox_t cell_bbox(alea_system_t* sys, size_t cell_index) {
    const alea_cell_entry_t* cell = &sys->cells.data[cell_index];
    if (cell->root_node_id == ALEA_NODE_ID_INVALID) {
        return alea_bbox_empty();
    }
    const alea_node_t* root = &sys->nodes.data[cell->root_node_id];
    if (root->bbox.min_x <= root->bbox.max_x) {
        return root->bbox;
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

int main(int argc, char** argv) {
    if (argc < 2) {
        fprintf(stderr, "usage: %s MODEL.mcnp [--queries N] [--hier-build]\n", argv[0]);
        fprintf(stderr, "set ALEA_DISABLE_UNIVERSE_POINT_BVH=1 to benchmark the old linear recursive path\n");
        fprintf(stderr, "set ALEA_HIER_BLAS_THRESHOLD=N to skip BLAS builds for universes below N cells\n");
        return 2;
    }

    size_t queries = 0;
    int hier_build = 0;
    for (int i = 2; i < argc; i++) {
        if (strcmp(argv[i], "--queries") == 0 && i + 1 < argc) {
            queries = (size_t)strtoull(argv[++i], NULL, 10);
        } else if (strcmp(argv[i], "--hier-build") == 0) {
            hier_build = 1;
        }
    }

    alea_log_set_level(1);
    const char* disabled = getenv("ALEA_DISABLE_UNIVERSE_POINT_BVH");
    printf("Universe point BVH: %s\n",
           (disabled && disabled[0] && strcmp(disabled, "0") != 0) ? "disabled" : "enabled");

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

    mcnp_model_destroy(model);
    return 0;
}
