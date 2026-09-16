// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/* Render an MCNP or OpenMC geometry on MPI ranks; write one PPM on root.
 * Usage: cluster_render INPUT OUTPUT.ppm [WIDTH HEIGHT] */

#include "alea_cluster.h"
#include "alea_mcnp.h"
#include "alea_openmc.h"

#include <errno.h>
#include <limits.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int ends_with(const char* path, const char* suffix) {
    size_t a = strlen(path), b = strlen(suffix);
    return a >= b && strcmp(path + a - b, suffix) == 0;
}

static int dimension(const char* input, int* output) {
    char* end;
    errno = 0;
    long value = strtol(input, &end, 10);
    if (errno || *end || value <= 0 || value > INT_MAX) return -1;
    *output = (int)value;
    return 0;
}

int main(int argc, char** argv) {
    int result = 1;
    alea_cluster_status_t status = alea_cluster_initialize(&argc, &argv);
    if (status != ALEA_CLUSTER_OK) return 1;
    alea_cluster_t* cluster = alea_cluster_create();
    if (!cluster) { alea_cluster_finalize(); return 1; }
    int rank = alea_cluster_rank(cluster);
    render_config_t cfg;
    render_config_init(&cfg);
    cfg.width = 640;
    cfg.height = 480;
    cfg.log_level = rank == 0 ? 1 : 0;
    char* input = NULL;
    size_t input_length = 0;
    mcnp_model_t* mcnp = NULL;
    openmc_model_t* openmc = NULL;
    render_framebuffer_t* frame = NULL;
    uint8_t* pixels = NULL;

    int valid = (argc == 3 || argc == 5) &&
        (argc != 5 || (dimension(argv[3], &cfg.width) == 0 &&
                       dimension(argv[4], &cfg.height) == 0));
    status = alea_cluster_agree(cluster, valid ? ALEA_CLUSTER_OK
                                         : ALEA_CLUSTER_INVALID_ARGUMENT);
    if (status != ALEA_CLUSTER_OK) {
        if (rank == 0)
            fprintf(stderr, "usage: %s INPUT OUTPUT.ppm [WIDTH HEIGHT]\n",
                    argv[0]);
        goto done;
    }
    int is_openmc = ends_with(argv[1], ".xml");
    status = is_openmc
        ? alea_cluster_read_file(cluster, rank == 0 ? argv[1] : NULL,
                                 &input, &input_length)
        : alea_cluster_read_mcnp_input(cluster,
            rank == 0 ? argv[1] : NULL, &input, &input_length);
    if (status != ALEA_CLUSTER_OK) {
        if (rank == 0) fprintf(stderr, "input: %s\n",
                               alea_cluster_status_string(status));
        goto done;
    }
    alea_system_t* sys = NULL;
    if (is_openmc) {
        openmc = openmc_load_string(input, input_length);
        if (openmc) sys = openmc_model_system(openmc);
    } else {
        mcnp = mcnp_load_string(input, input_length);
        if (mcnp) sys = mcnp_model_system(mcnp);
    }
    free(input);
    input = NULL;
    status = alea_cluster_agree(cluster, sys ? ALEA_CLUSTER_OK
                                           : ALEA_CLUSTER_COMPUTE_ERROR);
    if (status != ALEA_CLUSTER_OK) {
        if (rank == 0) fprintf(stderr, "model parsing: %s\n",
                               alea_cluster_status_string(status));
        goto done;
    }
    render_camera_t camera;
    status = alea_cluster_agree(cluster,
        render_camera_setup(&camera, &cfg, sys) == 0
            ? ALEA_CLUSTER_OK : ALEA_CLUSTER_COMPUTE_ERROR);
    if (status != ALEA_CLUSTER_OK) goto done;
    if (rank == 0) frame = render_framebuffer_create(
        cfg.width, cfg.height, 1);
    status = alea_cluster_agree(cluster,
        rank == 0 && !frame ? ALEA_CLUSTER_OUT_OF_MEMORY : ALEA_CLUSTER_OK);
    if (status != ALEA_CLUSTER_OK) goto done;
    status = alea_cluster_render_scene(cluster, sys, &cfg, &camera, frame);
    if (status != ALEA_CLUSTER_OK) {
        if (rank == 0) fprintf(stderr, "render: %s\n",
                               alea_cluster_status_string(status));
        goto done;
    }
    if (rank == 0) {
        size_t count = (size_t)cfg.width * cfg.height;
        if (count <= SIZE_MAX / 3)
            pixels = malloc(count * 3);
        if (!pixels) { fprintf(stderr, "output allocation failed\n"); goto done; }
        render_tonemap(frame, pixels);
        if (render_write_ppm(argv[2], pixels, cfg.width, cfg.height) != 0) {
            fprintf(stderr, "cannot write %s\n", argv[2]);
            goto done;
        }
    }
    result = 0;

done:
    free(pixels);
    render_framebuffer_free(frame);
    openmc_model_destroy(openmc);
    mcnp_model_destroy(mcnp);
    free(input);
    render_config_free(&cfg);
    alea_cluster_destroy(cluster);
    if (alea_cluster_finalize() != ALEA_CLUSTER_OK) result = 1;
    return result;
}
