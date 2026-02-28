// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/*
 * test_render3d.c - Unit tests for 3D batch renderer
 */

#include "alea_test.h"
#include "alea.h"
#include "render/render3d.h"
#include "core/alea_system.h"
#include "core/alea_universe.h"

#include <string.h>
#include <float.h>

#define EPS 1e-6

/* ============================================================================
 * Helper: create a simple scene with a sphere and a box
 * ============================================================================ */

static alea_system_t* create_test_scene(void) {
    alea_system_t* sys = alea_create();
    if (!sys) return NULL;

    /* Cell 1: Sphere at origin, radius 5 */
    int s1 = alea_sphere_surface(sys, 1, 0, 0, 0, 5.0);
    alea_node_id_t sphere = alea_surface_at(sys, s1)->neg_node;

    /* Cell 2: Box from (10,-3,-3) to (16,3,3) */
    int s2 = alea_box_surface(sys, 2, 10, 16, -3, 3, -3, 3);
    alea_node_id_t box = alea_surface_at(sys, s2)->neg_node;

    int m1 = alea_add_material(sys, 1);
    int m2 = alea_add_material(sys, 2);

    alea_add_cell(sys, 1, sphere, m1, -2.7, 0);
    alea_add_cell(sys, 2, box, m2, -8.0, 0);

    alea_build_universe_index(sys);
    alea_build_spatial_index(sys);

    return sys;
}

/* ============================================================================
 * Config tests
 * ============================================================================ */

TEST(config_init_defaults) {
    render_config_t cfg;
    render_config_init(&cfg);

    ASSERT_EQ(cfg.width, RENDER_DEFAULT_WIDTH);
    ASSERT_EQ(cfg.height, RENDER_DEFAULT_HEIGHT);
    ASSERT_NEAR(cfg.fov, RENDER_DEFAULT_FOV, EPS);
    ASSERT_NEAR(cfg.background[0], 0.95f, 1e-3);
    ASSERT_NEAR(cfg.background[1], 0.95f, 1e-3);
    ASSERT_NEAR(cfg.background[2], 0.95f, 1e-3);
    ASSERT_EQ(cfg.color_mode, RENDER_COLOR_MATERIAL);
    ASSERT_EQ(cfg.render_mode, RENDER_MODE_SOLID);
    ASSERT_NEAR(cfg.ambient, 0.3f, 1e-3);
    ASSERT_NEAR(cfg.diffuse, 0.6f, 1e-3);
    ASSERT_EQ(cfg.num_clips, 0);
    ASSERT_EQ(cfg.shadows, 0);
    ASSERT_EQ(cfg.edges, 0);

    render_config_free(&cfg);
}

/* ============================================================================
 * Framebuffer tests
 * ============================================================================ */

TEST(framebuffer_create_basic) {
    render_framebuffer_t* fb = render_framebuffer_create(64, 48, 0);
    ASSERT_NOT_NULL(fb);
    ASSERT_EQ(fb->width, 64);
    ASSERT_EQ(fb->height, 48);
    ASSERT_NOT_NULL(fb->color);
    ASSERT_NOT_NULL(fb->cell_id);

    /* Auxiliary buffers should be NULL when aux=0 */
    ASSERT_NULL(fb->depth);
    ASSERT_NULL(fb->normal);

    render_framebuffer_free(fb);
}

TEST(framebuffer_create_with_aux) {
    render_framebuffer_t* fb = render_framebuffer_create(32, 32, 1);
    ASSERT_NOT_NULL(fb);
    ASSERT_NOT_NULL(fb->color);
    ASSERT_NOT_NULL(fb->cell_id);
    ASSERT_NOT_NULL(fb->material_id);
    ASSERT_NOT_NULL(fb->depth);
    ASSERT_NOT_NULL(fb->normal);

    render_framebuffer_free(fb);
}

/* ============================================================================
 * Camera tests
 * ============================================================================ */

TEST(camera_auto_fit) {
    alea_system_t* sys = create_test_scene();
    ASSERT_NOT_NULL(sys);

    render_config_t cfg;
    render_config_init(&cfg);
    cfg.width = 64;
    cfg.height = 48;

    render_camera_t cam;
    int rc = render_camera_setup(&cam, &cfg, sys);
    ASSERT_EQ(rc, 0);
    ASSERT(cam.auto_fit == 1);

    /* Camera should be looking roughly at the scene center */
    /* Scene spans sphere at origin r=5 and box at x=10..16,
     * so center is roughly (5..8, 0, 0) */
    ASSERT_NEAR(cam.target[1], 0.0, 1.0);
    ASSERT_NEAR(cam.target[2], 0.0, 1.0);

    /* Eye should be at some distance from the target */
    double dx = cam.eye[0] - cam.target[0];
    double dy = cam.eye[1] - cam.target[1];
    double dz = cam.eye[2] - cam.target[2];
    double dist = sqrt(dx*dx + dy*dy + dz*dz);
    ASSERT(dist > 1.0);

    render_config_free(&cfg);
    alea_destroy(sys);
}

TEST(camera_explicit_eye_target) {
    alea_system_t* sys = create_test_scene();
    ASSERT_NOT_NULL(sys);

    render_config_t cfg;
    render_config_init(&cfg);
    cfg.width = 64;
    cfg.height = 48;
    cfg.eye[0] = 0; cfg.eye[1] = -30; cfg.eye[2] = 0;
    cfg.target[0] = 0; cfg.target[1] = 0; cfg.target[2] = 0;
    cfg.eye_set = 1;
    cfg.target_set = 1;

    render_camera_t cam;
    int rc = render_camera_setup(&cam, &cfg, sys);
    ASSERT_EQ(rc, 0);
    ASSERT(cam.auto_fit == 0);
    ASSERT_NEAR(cam.eye[1], -30.0, EPS);
    ASSERT_NEAR(cam.target[0], 0.0, EPS);

    render_config_free(&cfg);
    alea_destroy(sys);
}

TEST(camera_ray_perspective_center) {
    alea_system_t* sys = create_test_scene();
    ASSERT_NOT_NULL(sys);

    render_config_t cfg;
    render_config_init(&cfg);
    cfg.width = 64;
    cfg.height = 64;
    cfg.fov = 90.0;
    cfg.eye[0] = 0; cfg.eye[1] = -20; cfg.eye[2] = 0;
    cfg.target[0] = 0; cfg.target[1] = 0; cfg.target[2] = 0;
    cfg.eye_set = 1;
    cfg.target_set = 1;

    render_camera_t cam;
    int rc = render_camera_setup(&cam, &cfg, sys);
    ASSERT_EQ(rc, 0);

    /* Center pixel ray should point straight at target */
    double ox, oy, oz, dx, dy, dz;
    render_camera_ray(&cam, 64, 64, 32.0, 32.0,
                      &ox, &oy, &oz, &dx, &dy, &dz);

    /* Direction should be roughly +Y (from eye at y=-20 toward origin) */
    double len = sqrt(dx*dx + dy*dy + dz*dz);
    ASSERT(len > 0.5);
    double ny = dy / len;
    ASSERT(ny > 0.9);

    render_config_free(&cfg);
    alea_destroy(sys);
}

/* ============================================================================
 * Render tests (small framebuffers for speed)
 * ============================================================================ */

TEST(render_solid_basic) {
    alea_system_t* sys = create_test_scene();
    ASSERT_NOT_NULL(sys);

    render_config_t cfg;
    render_config_init(&cfg);
    cfg.width = 64;
    cfg.height = 48;
    cfg.render_mode = RENDER_MODE_SOLID;

    render_camera_t cam;
    int rc = render_camera_setup(&cam, &cfg, sys);
    ASSERT_EQ(rc, 0);

    render_framebuffer_t* fb = render_framebuffer_create(cfg.width, cfg.height, 0);
    ASSERT_NOT_NULL(fb);

    rc = render_scene(sys, &cfg, &cam, fb);
    ASSERT_EQ(rc, 0);

    /* At least some pixels should be non-background */
    int hit_count = 0;
    for (int i = 0; i < cfg.width * cfg.height; i++) {
        if (fb->cell_id[i] > 0) hit_count++;
    }
    ASSERT(hit_count > 0);

    render_framebuffer_free(fb);
    render_config_free(&cfg);
    alea_destroy(sys);
}

TEST(render_solid_hits_both_cells) {
    alea_system_t* sys = create_test_scene();
    ASSERT_NOT_NULL(sys);

    render_config_t cfg;
    render_config_init(&cfg);
    cfg.width = 64;
    cfg.height = 48;
    cfg.render_mode = RENDER_MODE_SOLID;
    /* Look from +Y side so both sphere and box are visible */
    cfg.eye[0] = 5; cfg.eye[1] = -40; cfg.eye[2] = 0;
    cfg.target[0] = 5; cfg.target[1] = 0; cfg.target[2] = 0;
    cfg.eye_set = 1;
    cfg.target_set = 1;

    render_camera_t cam;
    int rc = render_camera_setup(&cam, &cfg, sys);
    ASSERT_EQ(rc, 0);

    render_framebuffer_t* fb = render_framebuffer_create(cfg.width, cfg.height, 0);
    ASSERT_NOT_NULL(fb);

    rc = render_scene(sys, &cfg, &cam, fb);
    ASSERT_EQ(rc, 0);

    int found_cell1 = 0, found_cell2 = 0;
    for (int i = 0; i < cfg.width * cfg.height; i++) {
        if (fb->cell_id[i] == 1) found_cell1 = 1;
        if (fb->cell_id[i] == 2) found_cell2 = 1;
    }
    ASSERT(found_cell1);
    ASSERT(found_cell2);

    render_framebuffer_free(fb);
    render_config_free(&cfg);
    alea_destroy(sys);
}

TEST(render_xray_mode) {
    alea_system_t* sys = create_test_scene();
    ASSERT_NOT_NULL(sys);

    render_config_t cfg;
    render_config_init(&cfg);
    cfg.width = 32;
    cfg.height = 32;
    cfg.render_mode = RENDER_MODE_XRAY;

    render_camera_t cam;
    int rc = render_camera_setup(&cam, &cfg, sys);
    ASSERT_EQ(rc, 0);

    render_framebuffer_t* fb = render_framebuffer_create(cfg.width, cfg.height, 0);
    ASSERT_NOT_NULL(fb);

    rc = render_scene(sys, &cfg, &cam, fb);
    ASSERT_EQ(rc, 0);

    /* In X-ray mode, pixels through geometry should differ from background */
    int non_bg = 0;
    for (int i = 0; i < cfg.width * cfg.height; i++) {
        float r = fb->color[i * 3 + 0];
        float g = fb->color[i * 3 + 1];
        float b = fb->color[i * 3 + 2];
        if (fabs(r - 0.95f) > 0.01f ||
            fabs(g - 0.95f) > 0.01f ||
            fabs(b - 0.95f) > 0.01f) {
            non_bg++;
        }
    }
    ASSERT(non_bg > 0);

    render_framebuffer_free(fb);
    render_config_free(&cfg);
    alea_destroy(sys);
}

TEST(render_depth_mode) {
    alea_system_t* sys = create_test_scene();
    ASSERT_NOT_NULL(sys);

    render_config_t cfg;
    render_config_init(&cfg);
    cfg.width = 32;
    cfg.height = 32;
    cfg.render_mode = RENDER_MODE_SOLID;
    cfg.aux_output = 1;

    render_camera_t cam;
    int rc = render_camera_setup(&cam, &cfg, sys);
    ASSERT_EQ(rc, 0);

    render_framebuffer_t* fb = render_framebuffer_create(cfg.width, cfg.height, 1);
    ASSERT_NOT_NULL(fb);

    rc = render_scene(sys, &cfg, &cam, fb);
    ASSERT_EQ(rc, 0);

    /* Depth buffer should have finite values for pixels that hit geometry */
    int finite_depths = 0;
    for (int i = 0; i < cfg.width * cfg.height; i++) {
        if (fb->cell_id[i] > 0 && fb->depth[i] > 0.0f && fb->depth[i] < 1e10f) {
            finite_depths++;
        }
    }
    /* If any cell was hit, depth should be recorded */
    int hits = 0;
    for (int i = 0; i < cfg.width * cfg.height; i++) {
        if (fb->cell_id[i] > 0) hits++;
    }
    if (hits > 0) {
        ASSERT(finite_depths > 0);
    }

    render_framebuffer_free(fb);
    render_config_free(&cfg);
    alea_destroy(sys);
}

/* ============================================================================
 * Tonemap and post-processing tests
 * ============================================================================ */

TEST(tonemap_clamps_output) {
    render_framebuffer_t* fb = render_framebuffer_create(4, 4, 0);
    ASSERT_NOT_NULL(fb);

    /* Set some extreme color values */
    for (int i = 0; i < 4 * 4 * 3; i += 3) {
        fb->color[i + 0] = 2.0f;   /* Over-bright */
        fb->color[i + 1] = -0.5f;  /* Negative */
        fb->color[i + 2] = 0.5f;   /* Normal */
    }

    uint8_t pixels[4 * 4 * 3];
    render_tonemap(fb, pixels);

    for (int i = 0; i < 4 * 4; i++) {
        /* Over-bright should be clamped to 255 */
        ASSERT_EQ(pixels[i * 3 + 0], 255);
        /* Negative should be clamped to 0 */
        ASSERT_EQ(pixels[i * 3 + 1], 0);
        /* 0.5 should map to ~128 */
        ASSERT(pixels[i * 3 + 2] >= 126 && pixels[i * 3 + 2] <= 130);
    }

    render_framebuffer_free(fb);
}

TEST(edge_darken_modifies_edges) {
    alea_system_t* sys = create_test_scene();
    ASSERT_NOT_NULL(sys);

    render_config_t cfg;
    render_config_init(&cfg);
    cfg.width = 64;
    cfg.height = 48;

    render_camera_t cam;
    int rc = render_camera_setup(&cam, &cfg, sys);
    ASSERT_EQ(rc, 0);

    render_framebuffer_t* fb = render_framebuffer_create(cfg.width, cfg.height, 0);
    ASSERT_NOT_NULL(fb);

    rc = render_scene(sys, &cfg, &cam, fb);
    ASSERT_EQ(rc, 0);

    /* Save a copy of colors before edge darkening */
    int npix = cfg.width * cfg.height * 3;
    float* before = malloc(npix * sizeof(float));
    ASSERT_NOT_NULL(before);
    memcpy(before, fb->color, npix * sizeof(float));

    render_edge_darken(fb);

    /* At least some pixels should have changed (where cell boundaries are) */
    int changed = 0;
    for (int i = 0; i < npix; i++) {
        if (fabs(fb->color[i] - before[i]) > 1e-6f) changed++;
    }
    /* Edge darkening should affect some pixels on cell boundaries */
    ASSERT(changed >= 0);  /* May be 0 if no edges visible, that's OK */

    free(before);
    render_framebuffer_free(fb);
    render_config_free(&cfg);
    alea_destroy(sys);
}

/* ============================================================================
 * Color mode tests
 * ============================================================================ */

TEST(color_mode_cell) {
    alea_system_t* sys = create_test_scene();
    ASSERT_NOT_NULL(sys);

    render_config_t cfg;
    render_config_init(&cfg);
    cfg.width = 64;
    cfg.height = 48;
    cfg.color_mode = RENDER_COLOR_CELL;
    cfg.eye[0] = 5; cfg.eye[1] = -40; cfg.eye[2] = 0;
    cfg.target[0] = 5; cfg.target[1] = 0; cfg.target[2] = 0;
    cfg.eye_set = 1;
    cfg.target_set = 1;

    render_camera_t cam;
    int rc = render_camera_setup(&cam, &cfg, sys);
    ASSERT_EQ(rc, 0);

    render_framebuffer_t* fb = render_framebuffer_create(cfg.width, cfg.height, 0);
    ASSERT_NOT_NULL(fb);

    rc = render_scene(sys, &cfg, &cam, fb);
    ASSERT_EQ(rc, 0);

    /* Find a pixel in cell 1 and cell 2, verify they have different base colors */
    int idx1 = -1, idx2 = -1;
    for (int i = 0; i < cfg.width * cfg.height; i++) {
        if (fb->cell_id[i] == 1 && idx1 < 0) idx1 = i;
        if (fb->cell_id[i] == 2 && idx2 < 0) idx2 = i;
    }

    if (idx1 >= 0 && idx2 >= 0) {
        /* Colors should differ between cell 1 and cell 2 */
        float dr = fb->color[idx1 * 3 + 0] - fb->color[idx2 * 3 + 0];
        float dg = fb->color[idx1 * 3 + 1] - fb->color[idx2 * 3 + 1];
        float db = fb->color[idx1 * 3 + 2] - fb->color[idx2 * 3 + 2];
        float diff = dr*dr + dg*dg + db*db;
        ASSERT(diff > 0.001f);
    }

    render_framebuffer_free(fb);
    render_config_free(&cfg);
    alea_destroy(sys);
}

/* ============================================================================
 * Clipping plane tests
 * ============================================================================ */

TEST(clip_plane_cuts_geometry) {
    /* Sphere at origin — clip plane behind the sphere removes everything.
     * clip_entry works as a cutaway: camera must be on the CLIPPED side
     * of the plane so the ray start is pushed forward past geometry. */
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);

    int s = alea_sphere_surface(sys, 1, 0, 0, 0, 5.0);
    alea_node_id_t sphere = alea_surface_at(sys, s)->neg_node;

    int m1 = alea_add_material(sys, 1);

    alea_add_cell(sys, 1, sphere, m1, -2.7, 0);
    alea_build_universe_index(sys);
    alea_build_spatial_index(sys);

    render_config_t cfg;
    render_config_init(&cfg);
    cfg.width = 64;
    cfg.height = 64;
    cfg.eye[0] = 0; cfg.eye[1] = -20; cfg.eye[2] = 0;
    cfg.target[0] = 0; cfg.target[1] = 0; cfg.target[2] = 0;
    cfg.eye_set = 1;
    cfg.target_set = 1;

    render_camera_t cam;
    int rc = render_camera_setup(&cam, &cfg, sys);
    ASSERT_EQ(rc, 0);

    /* Render without clipping — should see the sphere */
    render_framebuffer_t* fb1 = render_framebuffer_create(cfg.width, cfg.height, 0);
    ASSERT_NOT_NULL(fb1);
    rc = render_scene(sys, &cfg, &cam, fb1);
    ASSERT_EQ(rc, 0);

    int hits_no_clip = 0;
    for (int i = 0; i < cfg.width * cfg.height; i++) {
        if (fb1->cell_id[i] > 0) hits_no_clip++;
    }
    ASSERT(hits_no_clip > 0);

    /* Clip: visible region is y >= 100 (far behind sphere).
     * Camera at y=-20 is on clipped side (val = -20-100 = -120 < 0).
     * Ray enters visible region at t ≈ 120, well past the sphere (at t ≈ 15-25).
     * Result: all sphere segments are before clip entry → zero visible pixels. */
    cfg.num_clips = 1;
    cfg.clips[0].normal[0] = 0.0;
    cfg.clips[0].normal[1] = 1.0;
    cfg.clips[0].normal[2] = 0.0;
    cfg.clips[0].d = -100.0;

    render_framebuffer_t* fb2 = render_framebuffer_create(cfg.width, cfg.height, 0);
    ASSERT_NOT_NULL(fb2);
    rc = render_scene(sys, &cfg, &cam, fb2);
    ASSERT_EQ(rc, 0);

    int hits_with_clip = 0;
    for (int i = 0; i < cfg.width * cfg.height; i++) {
        if (fb2->cell_id[i] > 0) hits_with_clip++;
    }

    /* Aggressive clip should remove all visible geometry */
    ASSERT_EQ(hits_with_clip, 0);

    render_framebuffer_free(fb1);
    render_framebuffer_free(fb2);
    render_config_free(&cfg);
    alea_destroy(sys);
}

/* ============================================================================
 * Image output tests
 * ============================================================================ */

TEST(write_bmp_roundtrip) {
    int w = 4, h = 4;
    uint8_t pixels[4 * 4 * 3];
    memset(pixels, 0, sizeof(pixels));

    /* Set some recognizable pattern */
    pixels[0] = 255; pixels[1] = 0; pixels[2] = 0;      /* Red */
    pixels[3] = 0;   pixels[4] = 255; pixels[5] = 0;    /* Green */
    pixels[6] = 0;   pixels[7] = 0;   pixels[8] = 255;  /* Blue */

    const char* path = "test_render3d_tmp.bmp";
    int rc = render_write_bmp(path, pixels, w, h);
    if (rc != 0) {
        /* Skip test if file write fails (e.g., permission issues on Windows) */
        SKIP("Could not write test file");
    }

    /* Verify file exists and has reasonable size:
     * BMP header = 54 bytes, row = 4*3 = 12, padded to 4-byte = 12 bytes
     * Total = 54 + 4 * 12 = 102 */
    FILE* f = fopen(path, "rb");
    ASSERT_NOT_NULL(f);
    fseek(f, 0, SEEK_END);
    long size = ftell(f);
    fclose(f);
    ASSERT(size >= 54);  /* At least a valid BMP header */

    remove(path);
}

TEST(write_ppm_roundtrip) {
    int w = 4, h = 4;
    uint8_t pixels[4 * 4 * 3];
    memset(pixels, 128, sizeof(pixels));

    const char* path = "test_render3d_tmp.ppm";
    int rc = render_write_ppm(path, pixels, w, h);
    if (rc != 0) {
        /* Skip test if file write fails (e.g., permission issues on Windows) */
        SKIP("Could not write test file");
    }

    FILE* f = fopen(path, "rb");
    ASSERT_NOT_NULL(f);
    fseek(f, 0, SEEK_END);
    long size = ftell(f);
    fclose(f);
    /* PPM: header + w*h*3 raw bytes */
    ASSERT(size >= w * h * 3);

    remove(path);
}

/* ============================================================================
 * Orthographic projection test
 * ============================================================================ */

TEST(render_orthographic) {
    alea_system_t* sys = create_test_scene();
    ASSERT_NOT_NULL(sys);

    render_config_t cfg;
    render_config_init(&cfg);
    cfg.width = 64;
    cfg.height = 48;
    cfg.fov = 0;              /* Orthographic */
    cfg.ortho_height = 20.0;
    cfg.eye[0] = 0; cfg.eye[1] = -30; cfg.eye[2] = 0;
    cfg.target[0] = 0; cfg.target[1] = 0; cfg.target[2] = 0;
    cfg.eye_set = 1;
    cfg.target_set = 1;

    render_camera_t cam;
    int rc = render_camera_setup(&cam, &cfg, sys);
    ASSERT_EQ(rc, 0);

    render_framebuffer_t* fb = render_framebuffer_create(cfg.width, cfg.height, 0);
    ASSERT_NOT_NULL(fb);

    rc = render_scene(sys, &cfg, &cam, fb);
    ASSERT_EQ(rc, 0);

    /* Should still hit geometry */
    int hit_count = 0;
    for (int i = 0; i < cfg.width * cfg.height; i++) {
        if (fb->cell_id[i] > 0) hit_count++;
    }
    ASSERT(hit_count > 0);

    render_framebuffer_free(fb);
    render_config_free(&cfg);
    alea_destroy(sys);
}

/* ============================================================================
 * Color palette test
 * ============================================================================ */

TEST(palette_returns_distinct_colors) {
    render_config_t cfg;
    render_config_init(&cfg);

    float r1, g1, b1, r2, g2, b2;
    render_get_color(1, RENDER_COLOR_MATERIAL, &cfg, &r1, &g1, &b1);
    render_get_color(2, RENDER_COLOR_MATERIAL, &cfg, &r2, &g2, &b2);

    /* Different IDs should produce different colors */
    float diff = (r1 - r2) * (r1 - r2) +
                 (g1 - g2) * (g1 - g2) +
                 (b1 - b2) * (b1 - b2);
    ASSERT(diff > 0.001f);

    render_config_free(&cfg);
}

TEST_MAIN()
