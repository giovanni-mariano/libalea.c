// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
// SPDX-License-Identifier: MPL-2.0
/* E-lite viewer workload. Compare complete warm surface frames at two camera
 * positions using one prepared handle. RSS fields are KiB on Linux. */
#include "alea.h"
#include "alea_mcnp.h"
#include "alea_render.h"
#include "util/alea_parallel.h"
#include "util/compat.h"
#include <math.h>
#include <stdio.h>
#include <string.h>

enum { WIDTH = 96, HEIGHT = 72, VIEWS = 2, REPEATS = 3 };

static int same_frame(const render_framebuffer_t* a, const render_framebuffer_t* b) {
    size_t n = (size_t)a->width * a->height;
    return !memcmp(a->color, b->color, 3*n*sizeof(float)) &&
           !memcmp(a->normal, b->normal, 3*n*sizeof(float)) &&
           !memcmp(a->depth, b->depth, n*sizeof(float)) &&
           !memcmp(a->cell_id, b->cell_id, n*sizeof(int)) &&
           !memcmp(a->material_id, b->material_id, n*sizeof(int));
}

static void print_memory(void) {
    FILE* f = fopen("/proc/self/status", "r");
    if (!f) return;
    char line[256];
    unsigned long kib;
    while (fgets(line, sizeof(line), f)) {
        if (sscanf(line, "VmRSS: %lu", &kib) == 1) printf(" rss_kib=%lu", kib);
        if (sscanf(line, "VmHWM: %lu", &kib) == 1) printf(" peak_kib=%lu", kib);
    }
    fclose(f);
}

static double median3(const double t[REPEATS]) {
    return fmax(fmin(t[0], t[1]), fmin(fmax(t[0], t[1]), t[2]));
}

int main(int argc, char** argv) {
    if (argc != 2) { fprintf(stderr, "usage: %s E-lite_R200430_mod\n", argv[0]); return 2; }
    render_config_t cfg;
    render_config_init(&cfg);
    cfg.width = WIDTH;
    cfg.height = HEIGHT;
    cfg.tile_size = 8;
    cfg.log_level = 0;
    cfg.eye_set = cfg.target_set = 1;
    cfg.target[0] = 800;
    const double eyes[VIEWS][3] = {{1800,-600,350}, {1900,120,420}};
    double span = 2 * hypot(hypot(1000, 600), 350);
    alea_bbox_t domain = {800-span,800+span,-span,span,-span,span};
    mcnp_model_t* m = mcnp_load(argv[1]);
    if (!m) return 1;
    int rc = 1;
    alea_prepared_queries_t* prepared = NULL;
    render_framebuffer_t* reference[VIEWS] = {0};
    render_framebuffer_t* frame = NULL;
    render_camera_t cameras[VIEWS];
    if (alea_prepare_query_acceleration(m->sys)) goto done;
    prepared = alea_prepared_queries_create(m->sys, &domain);
    if (!prepared) goto done;
    printf("ready"); print_memory(); puts(""); fflush(stdout);
    frame = render_framebuffer_create(WIDTH, HEIGHT, 1);
    if (!frame) goto done;
    for (int view = 0; view < VIEWS; ++view) {
        memcpy(cfg.eye, eyes[view], sizeof(cfg.eye));
        if (render_camera_setup(&cameras[view], &cfg, m->sys)) goto done;
        reference[view] = render_framebuffer_create(WIDTH, HEIGHT, 1);
        cfg.threads = 1;
        if (!reference[view] || alea_render_scene_prepared(m->sys, prepared, &cfg,
                &cameras[view], reference[view])) goto done;
    }
    for (cfg.threads = 1; cfg.threads <= 8; cfg.threads *= 2) {
        size_t tiles = ((WIDTH+cfg.tile_size-1)/cfg.tile_size) * ((HEIGHT+cfg.tile_size-1)/cfg.tile_size);
        size_t workers = alea_parallel_effective_workers(tiles, 1, cfg.threads);
        double samples[VIEWS][REPEATS];
        for (int pass = -1; pass < REPEATS; ++pass) {
            for (int view = 0; view < VIEWS; ++view) {
                double start = alea_monotonic_seconds();
                if (alea_render_scene_prepared(m->sys, prepared, &cfg, &cameras[view], frame)) goto done;
                double elapsed = alea_monotonic_seconds() - start;
                if (!same_frame(frame, reference[view])) {
                    fprintf(stderr, "Frame mismatch with %zu workers at view %d\n", workers, view);
                    goto done;
                }
                if (pass >= 0) samples[view][pass] = elapsed;
            }
        }
        printf("requested=%d workers=%zu view1_s=%.4f view2_s=%.4f",
               cfg.threads, workers, median3(samples[0]), median3(samples[1]));
        print_memory(); puts(""); fflush(stdout);
    }
    rc = 0;
done:
    for (int view = 0; view < VIEWS; ++view) render_framebuffer_free(reference[view]);
    render_framebuffer_free(frame);
    alea_prepared_queries_destroy(prepared);
    mcnp_model_destroy(m);
    render_config_free(&cfg);
    return rc;
}
