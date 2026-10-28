// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
// SPDX-License-Identifier: MPL-2.0
#include "alea_test.h"
#include "alea.h"
#include "core/alea_system.h"
#include "raycast/raycast.h"
#include "render/render3d.h"
#include "util/alea_alloc.h"
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#ifdef _WIN32
#include <direct.h>
#define test_mkdir(path) _mkdir(path)
#define test_rmdir(path) _rmdir(path)
#else
#include <unistd.h>
#include <signal.h>
#include <sys/resource.h>
#include <sys/wait.h>
#define test_mkdir(path) mkdir(path, 0700)
#define test_rmdir(path) rmdir(path)
#endif

static alea_system_t* sphere_system(double x) {
    alea_system_t* sys = alea_create();
    if (!sys) return NULL;
    int surface = alea_sphere_surface(sys, 1, x, 0, 0, 1);
    if (surface < 0 || alea_add_cell(sys, 1, alea_halfspace(sys, surface, -1),
                                    ALEA_MATERIAL_VOID, 0, 0) < 0) {
        alea_destroy(sys); return NULL;
    }
    return sys;
}

TEST(merge_invalidates_prepared_point_and_ray_queries) {
    alea_system_t* target = sphere_system(0);
    alea_system_t* source = sphere_system(10);
    ASSERT_NOT_NULL(target); ASSERT_NOT_NULL(source);
    ASSERT_EQ(alea_prepare_query_acceleration(target), 0);
    uint64_t generation = alea_system_geometry_generation(target);
    ASSERT_EQ(alea_find_cell(target, 10, 0, 0), -1);
    ASSERT_EQ(alea_merge(target, source, 100), 1);
    ASSERT(alea_system_geometry_generation(target) > generation);
    ASSERT_FALSE(alea_system_query_cache_ready(target, ALEA_CACHE_ALL));
    ASSERT_EQ(alea_prepare_query_acceleration(target), 0);
    ASSERT_EQ(alea_find_cell(target, 10, 0, 0), 1);
    ASSERT_EQ(alea_find_cell(target, 0, 0, 0), 0);
    alea_raycast_result_t ray;
    alea_raycast_result_init(&ray);
    ASSERT_EQ(alea_raycast_cell_aware(target, 8, 0, 0, 1, 0, 0, 4, &ray), 0);
    int found = 0;
    for (size_t i = 0; i < ray.segments.count; i++) {
        const alea_ray_segment_t* segment = &ray.segments.data[i];
        if (segment->cell_id != 101) continue;
        found = 1;
        ASSERT_NEAR(segment->t_enter, 1, 1e-8);
        ASSERT_NEAR(segment->t_exit, 3, 1e-8);
        ASSERT_EQ(segment->enter_surface_id, 101);
        ASSERT_EQ(segment->exit_surface_id, 101);
    }
    ASSERT(found);
    alea_raycast_result_free(&ray);
    alea_destroy(source); alea_destroy(target);
}

TEST(merge_allocation_failure_preserves_existing_geometry) {
    for (long failure = 0; failure < 3; failure++) {
        alea_system_t* target = sphere_system(0);
        alea_system_t* source = sphere_system(10);
        ASSERT_NOT_NULL(target); ASSERT_NOT_NULL(source);
        source->cells.data[0].lat_fill = malloc(sizeof(int));
        ASSERT_NOT_NULL(source->cells.data[0].lat_fill);
        source->cells.data[0].lat_fill[0] = 7;
        source->cells.data[0].lat_fill_count = 1;
        ASSERT_EQ(alea_prepare_query_acceleration(target), 0);
        size_t nodes = target->nodes.count, surfaces = target->surfaces.count;
        alea_error_clear();
        alea_alloc_fail_after(failure);
        int result = alea_merge(target, source, 100);
        alea_alloc_fail_after(-1);
        ASSERT_EQ(result, -1);
        ASSERT_EQ(target->cells.count, 1);
        ASSERT_EQ(target->nodes.count, nodes);
        ASSERT_EQ(target->surfaces.count, surfaces);
        ASSERT_EQ(alea_prepare_query_acceleration(target), 0);
        ASSERT_EQ(alea_find_cell(target, 0, 0, 0), 0);
        ASSERT_EQ(alea_find_cell(target, 10, 0, 0), -1);
        ASSERT_EQ(alea_merge(target, source, 100), 1);
        ASSERT_EQ(alea_prepare_query_acceleration(target), 0);
        ASSERT_EQ(alea_find_cell(target, 10, 0, 0), 1);
        ASSERT(target->cells.data[1].lat_fill != source->cells.data[0].lat_fill);
        alea_destroy(source); alea_destroy(target);
    }
}

TEST(surface_index_failure_is_reported_and_retryable) {
    alea_system_t* sys = sphere_system(0);
    ASSERT_NOT_NULL(sys);
    ASSERT_EQ(alea_build_cell_surface_index(sys), 0);
    alea_alloc_fail_after(0);
    int result = alea_build_cell_surface_index(sys);
    alea_alloc_fail_after(-1);
    ASSERT_EQ(result, -1);
    ASSERT_EQ(alea_error_code(), ALEA_ERR_OUT_OF_MEMORY);
    ASSERT_FALSE(alea_system_query_cache_ready(sys, ALEA_CACHE_CELL_SURFACES));
    ASSERT_EQ(alea_build_cell_surface_index(sys), 0);
    ASSERT_EQ(sys->cells.data[0].surface_index_count, 1);
    ASSERT(alea_system_query_cache_ready(sys, ALEA_CACHE_CELL_SURFACES));
    alea_destroy(sys);
}

TEST(color_load_failure_preserves_previous_configuration) {
    const char* path = "test_resource_colors.tmp";
    FILE* file = fopen(path, "w");
    ASSERT_NOT_NULL(file);
    for (int i = 64; i >= 0; i--) fprintf(file, "%d 1 0 0\n", i);
    ASSERT_EQ(fclose(file), 0);
    render_config_t cfg;
    render_config_init(&cfg);
    ASSERT_EQ(render_load_colors(path, &cfg), 0);
    ASSERT_EQ(cfg.num_custom_colors, 65);
    ASSERT_EQ(cfg.custom_colors[0].id, 0);
    render_color_entry_t* previous = cfg.custom_colors;
    for (long failure = 0; failure < 2; failure++) {
        alea_alloc_fail_after(failure);
        int result = render_load_colors(path, &cfg);
        alea_alloc_fail_after(-1);
        ASSERT_EQ(result, -1);
        ASSERT_EQ(alea_error_code(), ALEA_ERR_OUT_OF_MEMORY);
        ASSERT(cfg.custom_colors == previous);
        ASSERT_EQ(cfg.num_custom_colors, 65);
        ASSERT_EQ(cfg.custom_colors[64].id, 64);
    }
    ASSERT_EQ(render_load_colors(path, &cfg), 0);
    ASSERT_EQ(cfg.num_custom_colors, 65);
    render_config_free(&cfg);
    remove(path);
}

TEST(auxiliary_output_propagates_each_write_failure) {
    int cell = 1, material = 2;
    float depth = 1, normal[3] = {0,0,1};
    const char* paths[] = {"test_resource_aux_depth.png", "test_resource_aux_cellid.png",
                           "test_resource_aux_matid.png", "test_resource_aux_normal.png"};
    for (int buffer = 0; buffer < 4; buffer++) {
        render_framebuffer_t fb = {0};
        fb.width = fb.height = 1;
        fb.cell_id = &cell; fb.depth = &depth;
        fb.material_id = &material; fb.normal = normal;
        /* A directory at this output filename makes only this write fail.
         * Earlier maps must be written successfully before reaching it. */
        ASSERT_EQ(test_mkdir(paths[buffer]), 0);
        int result = render_write_aux("test_resource_aux.png", &fb);
        ASSERT_EQ(test_rmdir(paths[buffer]), 0);
        for (int i = 0; i < 4; i++) if (i != buffer) remove(paths[i]);
        ASSERT_EQ(result, -1);
        ASSERT_EQ(alea_error_code(), ALEA_ERR_FILE_WRITE);
    }
}

TEST(image_writer_allocation_failures_are_reported) {
    uint8_t pixels[12] = {0};
    for (int format = 0; format < 2; format++) {
        for (long failure = 0; failure < (format ? 1 : 2); failure++) {
            alea_alloc_fail_after(failure);
            int rc = format ? render_write_bmp("test_resource_image.bmp", pixels, 2, 2)
                            : render_write_png("test_resource_image.png", pixels, 2, 2);
            alea_alloc_fail_after(-1);
            ASSERT_EQ(rc, -1);
            ASSERT_EQ(alea_error_code(), ALEA_ERR_OUT_OF_MEMORY);
        }
    }
}

#ifndef _WIN32
/* Small images fail at buffered close, large images during fwrite. Keep the
 * process-wide file limit in a child so logs and other tests stay writable. */
TEST(image_writers_report_truncated_output_and_auxiliary_failures) {
    for (int format = 0; format < 7; format++) {
        for (int large = 0; large < 2; large++) {
            pid_t child = fork();
            ASSERT(child >= 0);
            if (child == 0) {
                struct rlimit limit = {16, 16};
                signal(SIGXFSZ, SIG_IGN);
                if (setrlimit(RLIMIT_FSIZE, &limit)) _exit(2);
                uint8_t pixels[128 * 128 * 3] = {0};
                int cells[128 * 128] = {0};
                int dim = large ? 128 : 2;
                const char* ext[] = {"png", "bmp", "ppm"};
                char path[128], output[160];
                snprintf(path, sizeof(path), "test_resource_image_%ld.%s",
                         (long)getpid(), ext[format % 3]);
                int rc;
                if (format == 0) rc = render_write_png(path, pixels, dim, dim);
                else if (format == 1) rc = render_write_bmp(path, pixels, dim, dim);
                else if (format == 2) rc = render_write_ppm(path, pixels, dim, dim);
                else if (format < 6) rc = render_write_image(path, pixels, dim, dim);
                else {
                    render_framebuffer_t fb = {0};
                    fb.width = fb.height = dim; fb.cell_id = cells;
                    rc = render_write_aux(path, &fb);
                }
                if (format == 6) {
                    char* dot = strrchr(path, '.');
                    if (dot) *dot = '\0';
                    snprintf(output, sizeof(output), "%s_cellid.png", path);
                } else snprintf(output, sizeof(output), "%s", path);
                struct stat st;
                int correct = rc == -1 && alea_error_code() == ALEA_ERR_FILE_WRITE &&
                              stat(output, &st) == 0 && st.st_size == 16;
                remove(output);
                _exit(correct ? 0 : 3);
            }
            int status;
            ASSERT_EQ(waitpid(child, &status, 0), child);
            ASSERT(WIFEXITED(status));
            ASSERT_EQ(WEXITSTATUS(status), 0);
        }
    }
}
#endif

TEST(png_writer_produces_complete_chunks) {
    const char* path = "test_resource_image.png";
    uint8_t pixels[12] = {0};
    ASSERT_EQ(render_write_png(path, pixels, 2, 2), 0);
    FILE* f = fopen(path, "rb");
    ASSERT_NOT_NULL(f);
    unsigned char signature[8], end[12];
    ASSERT_EQ(fread(signature, 1, sizeof(signature), f), sizeof(signature));
    ASSERT_EQ(memcmp(signature, "\211PNG\r\n\032\n", 8), 0);
    ASSERT_EQ(fseek(f, -12, SEEK_END), 0);
    ASSERT_EQ(fread(end, 1, sizeof(end), f), sizeof(end));
    ASSERT_EQ(fclose(f), 0);
    remove(path);
    const unsigned char expected[12] = {0,0,0,0,'I','E','N','D',0xae,0x42,0x60,0x82};
    ASSERT_EQ(memcmp(end, expected, sizeof(end)), 0);
}

TEST_MAIN()
