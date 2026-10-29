// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
// SPDX-License-Identifier: MPL-2.0

/* Exercise cancellation at the actual plotter boundary, including reuse. */
#define main plotter_main
#include "../../tools/mc_plotter.c"
#undef main
#include <assert.h>

typedef struct { int calls, stop; } cancellation_t;
static int cancel_after(void* data) {
    cancellation_t* c = data;
    return ++c->calls >= c->stop;
}

int main(int argc, char** argv) {
    assert(argc == 2);
    loaded_model_t model;
    assert(load_model(argv[1], &model) == 0);
    plot_params_t p = {0};
    p.type = SLICE_Z;
    p.u_min = p.v_min = -10;
    p.u_max = p.v_max = 10;
    p.width = p.height = 96; /* Three GPU submissions. */
    p.metal_backend = 1;
    p.compute = alea_compute_create(model.sys, ALEA_COMPUTE_METAL, ALEA_COMPUTE_FP32);
    assert(p.compute);
    uint8_t* baseline = NULL;
    assert(render_plot(model.sys, &p, 0, &baseline) == 0 && baseline);
    for (int stop = 1; stop <= 4; ++stop) {
        cancellation_t cancel = {0, stop};
        p.cancelled = cancel_after;
        p.cancel_data = &cancel;
        uint8_t* image = baseline; /* Failure must explicitly clear the output. */
        assert(render_plot(model.sys, &p, 0, &image) != 0);
        assert(image == NULL && cancel.calls == stop);
        p.cancelled = NULL;
        assert(render_plot(model.sys, &p, 0, &image) == 0 && image);
        assert(memcmp(baseline, image, (size_t)p.width * p.height * 3) == 0);
        free(image);
    }
    free(baseline);
    alea_compute_destroy(p.compute);
    destroy_model(&model);
    puts("Metal slice: cancellation before/between submissions and before publication, reuse PASS");
    return 0;
}
