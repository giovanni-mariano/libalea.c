// SPDX-License-Identifier: MPL-2.0
/*
 * slice_raycast_bench — compare the scanline-raytrace slice fill
 * (alea_find_cells_grid_raycast) against the point-query fill
 * (alea_find_cells_grid): timing + per-pixel cell agreement.
 *
 * Usage:
 *   slice_raycast_bench model.i axis value umin umax vmin vmax [res ...]
 *     axis: 0=X-normal, 1=Y-normal, 2=Z-normal
 *   defaults: a few square resolutions if none given.
 *
 * Build:
 *   cc -O2 -fopenmp -o slice_raycast_bench slice_raycast_bench.c \
 *      -Iinclude -Isrc -Lbin -lalea_full -lm
 */
#include "alea.h"
#include "alea_mcnp.h"
#include "alea_slice.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

static double now(void) {
    struct timespec t; clock_gettime(CLOCK_MONOTONIC, &t);
    return t.tv_sec + t.tv_nsec * 1e-9;
}

/* Fraction of disagreements that touch a cell boundary in the point-query
 * result (4-connected) — these are the expected sub-pixel boundary mismatches. */
static void compare(alea_system_t* sys, const alea_slice_view_t* view,
                    const int* a, const int* b, int nu, int nv) {
    long total = (long)nu * nv, diff = 0, diff_at_boundary = 0;
    int shown = 0;
    const alea_slice_plane_t* p = &view->plane;
    double du = (view->u_max - view->u_min) / nu;
    double dv = (view->v_max - view->v_min) / nv;
    for (int j = 0; j < nv; j++) {
        for (int i = 0; i < nu; i++) {
            long idx = (long)j * nu + i;
            if (a[idx] == b[idx]) continue;
            diff++;
            int bnd = 0;
            if (i > 0      && a[idx-1]  != a[idx]) bnd = 1;
            if (i < nu-1   && a[idx+1]  != a[idx]) bnd = 1;
            if (j > 0      && a[idx-nu] != a[idx]) bnd = 1;
            if (j < nv-1   && a[idx+nu] != a[idx]) bnd = 1;
            if (bnd) diff_at_boundary++;
            /* For a few NON-boundary diffs, print ground truth (find_cell_at). */
            else if (shown < 6) {
                double u = view->u_min + (i + 0.5) * du;
                double v = view->v_min + (j + 0.5) * dv;
                double x = p->origin[0] + u*p->u_axis[0] + v*p->v_axis[0];
                double y = p->origin[1] + u*p->u_axis[1] + v*p->v_axis[1];
                double z = p->origin[2] + u*p->u_axis[2] + v*p->v_axis[2];
                int gc = -1, gm = 0;
                alea_find_cell_at(sys, x, y, z, &gc, &gm);
                printf("      diff @(%d,%d) pt=(%.2f,%.2f,%.2f): pointquery=%d raytrace=%d  find_cell_at=%d\n",
                       i, j, x, y, z, a[idx], b[idx], gc);
                shown++;
            }
        }
    }
    double agree = 100.0 * (double)(total - diff) / total;
    printf("    agreement=%.4f%%  diffs=%ld (%.3f%%)  of-which-at-boundary=%ld (%.1f%%)\n",
           agree, diff, 100.0*diff/total, diff_at_boundary,
           diff ? 100.0*diff_at_boundary/diff : 0.0);
}

int main(int argc, char** argv) {
    if (argc < 8) {
        fprintf(stderr, "usage: %s model axis value umin umax vmin vmax [res ...]\n", argv[0]);
        return 1;
    }
    const char* path = argv[1];
    int axis = atoi(argv[2]);
    double value = atof(argv[3]);
    double umin = atof(argv[4]), umax = atof(argv[5]);
    double vmin = atof(argv[6]), vmax = atof(argv[7]);

    int res_default[] = {512, 1024, 2048, 4096};
    int n_res; int* res;
    if (argc > 8) {
        n_res = argc - 8; res = malloc(n_res * sizeof(int));
        for (int k = 0; k < n_res; k++) res[k] = atoi(argv[8 + k]);
    } else {
        n_res = 4; res = res_default;
    }

    double t0 = now();
    mcnp_model_t* m = mcnp_load(path);
    if (!m) { fprintf(stderr, "load failed: %s\n", alea_error()); return 1; }
    alea_system_t* s = m->sys;
    alea_build_universe_index(s);
    alea_prepare_query_acceleration(s);
    printf("loaded %s: %zu cells, %.1fs\n", path, alea_cell_count(s), now()-t0);
    printf("slice axis=%d value=%g  u=[%g,%g] v=[%g,%g]\n",
           axis, value, umin, umax, vmin, vmax);

    alea_slice_view_t view;
    alea_slice_view_axis(&view, axis, value, umin, umax, vmin, vmax);

    for (int k = 0; k < n_res; k++) {
        int n = res[k];
        size_t npx = (size_t)n * n;
        int* cg = malloc(npx * sizeof(int));
        int* mg = malloc(npx * sizeof(int));
        int* cr = malloc(npx * sizeof(int));
        int* mr = malloc(npx * sizeof(int));
        if (!cg || !mg || !cr || !mr) { fprintf(stderr, "oom at %d\n", n); return 1; }

        double tg0 = now();
        alea_find_cells_grid(s, &view, n, n, -1, cg, mg, NULL);
        double tg = now() - tg0;

        double tr0 = now();
        alea_find_cells_grid_raycast(s, &view, n, n, -1, cr, mr, NULL);
        double tr = now() - tr0;

        printf("  %dx%d:  point-query=%.3fs  raytrace=%.3fs  speedup=%.1fx\n",
               n, n, tg, tr, tr > 0 ? tg/tr : 0.0);
        compare(s, &view, cg, cr, n, n);

        free(cg); free(mg); free(cr); free(mr);
    }

    mcnp_model_destroy(m);
    if (argc > 8) free(res);
    return 0;
}
