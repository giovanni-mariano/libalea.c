// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/**
 * @file render3d.c
 * @brief 3D batch renderer implementation
 *
 * Camera, clipping, shading, tile-based rendering, image output.
 * Uses internal raycast API for buffer reuse (avoids malloc per pixel).
 */

#define _USE_MATH_DEFINES
#include "render3d.h"
#include "raycast/raycast.h"
#include "raycast/ray_intersect.h"
#include "core/alea_system.h"
#include "core/alea_universe.h"

#include <errno.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <stdio.h>
#include "util/math.h"

#include "util/compat.h"

#ifdef _OPENMP
#include <omp.h>
#endif

/* ============================================================================
 * COLOR PALETTE (32 qualitative colors)
 * ============================================================================ */

static const float palette[RENDER_PALETTE_SIZE][3] = {
    {0.122f, 0.467f, 0.706f},  /* Blue */
    {1.000f, 0.498f, 0.055f},  /* Orange */
    {0.173f, 0.627f, 0.173f},  /* Green */
    {0.839f, 0.153f, 0.157f},  /* Red */
    {0.580f, 0.404f, 0.741f},  /* Purple */
    {0.549f, 0.337f, 0.294f},  /* Brown */
    {0.890f, 0.467f, 0.761f},  /* Pink */
    {0.498f, 0.498f, 0.498f},  /* Gray */
    {0.737f, 0.741f, 0.133f},  /* Olive */
    {0.090f, 0.745f, 0.812f},  /* Cyan */
    {0.682f, 0.780f, 0.910f},  /* Light Blue */
    {1.000f, 0.733f, 0.471f},  /* Light Orange */
    {0.596f, 0.875f, 0.541f},  /* Light Green */
    {1.000f, 0.596f, 0.588f},  /* Light Red */
    {0.773f, 0.690f, 0.835f},  /* Light Purple */
    {0.769f, 0.612f, 0.580f},  /* Light Brown */
    {0.400f, 0.761f, 0.647f},  /* Teal */
    {0.988f, 0.553f, 0.384f},  /* Salmon */
    {0.553f, 0.627f, 0.796f},  /* Steel Blue */
    {0.906f, 0.541f, 0.765f},  /* Orchid */
    {0.651f, 0.847f, 0.329f},  /* Yellow Green */
    {1.000f, 0.851f, 0.184f},  /* Gold */
    {0.898f, 0.769f, 0.580f},  /* Tan */
    {0.702f, 0.702f, 0.702f},  /* Silver */
    {0.200f, 0.627f, 0.173f},  /* Forest */
    {0.984f, 0.604f, 0.600f},  /* Coral */
    {0.369f, 0.310f, 0.635f},  /* Indigo */
    {0.600f, 0.600f, 0.200f},  /* Dark Olive */
    {0.400f, 0.200f, 0.000f},  /* Chocolate */
    {0.000f, 0.502f, 0.502f},  /* Dark Teal */
    {0.863f, 0.078f, 0.235f},  /* Crimson */
    {0.255f, 0.412f, 0.882f},  /* Royal Blue */
};

/* ============================================================================
 * CONFIG INITIALIZATION
 * ============================================================================ */

void render_config_init(render_config_t* cfg) {
    memset(cfg, 0, sizeof(*cfg));

    cfg->width = RENDER_DEFAULT_WIDTH;
    cfg->height = RENDER_DEFAULT_HEIGHT;
    cfg->background[0] = 0.95f;
    cfg->background[1] = 0.95f;
    cfg->background[2] = 0.95f;

    cfg->up[2] = 1.0;  /* Z-up */
    cfg->fov = RENDER_DEFAULT_FOV;
    cfg->ortho_height = 0;

    cfg->color_mode = RENDER_COLOR_MATERIAL;
    cfg->render_mode = RENDER_MODE_SOLID;

    cfg->ambient = 0.3f;
    cfg->diffuse = 0.6f;
    cfg->specular = 0.2f;
    cfg->shininess = 30.0f;
    cfg->cross_section_tint = 0.85f;

    cfg->aa_samples = 1;
    cfg->universe_depth = -1;
    cfg->tile_size = RENDER_DEFAULT_TILE;
    cfg->xray_density_scale = 0.1f;
    cfg->log_level = 2;
}

void render_config_free(render_config_t* cfg) {
    free(cfg->custom_colors);
    cfg->custom_colors = NULL;
    cfg->num_custom_colors = 0;
}

/* ============================================================================
 * CAMERA
 * ============================================================================ */

int render_camera_setup(render_camera_t* cam,
                        const render_config_t* cfg,
                        alea_system_t* sys) {
    memset(cam, 0, sizeof(*cam));

    cam->up[0] = cfg->up[0];
    cam->up[1] = cfg->up[1];
    cam->up[2] = cfg->up[2];
    if (v3arr_length(cam->up) < 1e-10) {
        cam->up[2] = 1.0;
    }

    cam->fov_degrees = cfg->fov;
    cam->ortho_height = cfg->ortho_height;

    /* Auto-fit if eye/target not set */
    if (cfg->eye_set && cfg->target_set) {
        memcpy(cam->eye, cfg->eye, sizeof(cam->eye));
        memcpy(cam->target, cfg->target, sizeof(cam->target));
        cam->auto_fit = 0;
    } else {
        double cx, cy, cz, radius;
        int rc = alea_compute_bounding_sphere(sys, 1.0, &cx, &cy, &cz, &radius);
        if (rc != 0 || radius <= 0) {
            fprintf(stderr, "render: failed to compute bounding sphere\n");
            return -1;
        }

        if (cfg->target_set) {
            memcpy(cam->target, cfg->target, sizeof(cam->target));
        } else {
            cam->target[0] = cx;
            cam->target[1] = cy;
            cam->target[2] = cz;
        }

        if (cfg->eye_set) {
            memcpy(cam->eye, cfg->eye, sizeof(cam->eye));
        } else {
            /* Isometric direction, 2.5x radius */
            double dist = radius * 2.5;
            double iso[3] = {1.0, 1.0, 0.8};
            v3arr_normalize(iso);
            cam->eye[0] = cam->target[0] + iso[0] * dist;
            cam->eye[1] = cam->target[1] + iso[1] * dist;
            cam->eye[2] = cam->target[2] + iso[2] * dist;
        }

        cam->auto_fit = 1;

        /* Auto ortho_height if orthographic and not set */
        if (cam->fov_degrees <= 0 && cam->ortho_height <= 0) {
            cam->ortho_height = radius * 2.2;
        }
    }

    /* Compute derived quantities */
    v3arr_sub(cam->target, cam->eye, cam->forward);
    v3arr_normalize(cam->forward);

    v3arr_cross(cam->forward, cam->up, cam->right);
    if (v3arr_length(cam->right) < 1e-10) {
        /* Forward parallel to up — use fallback */
        double fallback_up[3] = {0, 1, 0};
        if (fabs(cam->forward[1]) > 0.9) {
            fallback_up[0] = 1; fallback_up[1] = 0; fallback_up[2] = 0;
        }
        v3arr_cross(cam->forward, fallback_up, cam->right);
    }
    v3arr_normalize(cam->right);

    v3arr_cross(cam->right, cam->forward, cam->true_up);
    v3arr_normalize(cam->true_up);

    /* Perspective half-extents */
    double aspect = (double)cfg->width / (double)cfg->height;
    if (cam->fov_degrees > 0) {
        cam->half_height = tan(cam->fov_degrees * M_PI / 360.0);
        cam->half_width = cam->half_height * aspect;
    } else {
        cam->half_height = cam->ortho_height * 0.5;
        cam->half_width = cam->half_height * aspect;
    }

    return 0;
}

void render_camera_ray(const render_camera_t* cam,
                       int width, int height,
                       double px, double py,
                       double* ox, double* oy, double* oz,
                       double* dx, double* dy, double* dz) {
    double u = (2.0 * px / width - 1.0) * cam->half_width;
    double v = (1.0 - 2.0 * py / height) * cam->half_height;

    if (cam->fov_degrees > 0) {
        /* Perspective */
        *ox = cam->eye[0];
        *oy = cam->eye[1];
        *oz = cam->eye[2];
        double d[3];
        d[0] = cam->forward[0] + u * cam->right[0] + v * cam->true_up[0];
        d[1] = cam->forward[1] + u * cam->right[1] + v * cam->true_up[1];
        d[2] = cam->forward[2] + u * cam->right[2] + v * cam->true_up[2];
        double len = sqrt(d[0]*d[0] + d[1]*d[1] + d[2]*d[2]);
        *dx = d[0] / len;
        *dy = d[1] / len;
        *dz = d[2] / len;
    } else {
        /* Orthographic */
        *ox = cam->eye[0] + u * cam->right[0] + v * cam->true_up[0];
        *oy = cam->eye[1] + u * cam->right[1] + v * cam->true_up[1];
        *oz = cam->eye[2] + u * cam->right[2] + v * cam->true_up[2];
        *dx = cam->forward[0];
        *dy = cam->forward[1];
        *dz = cam->forward[2];
    }
}

/* ============================================================================
 * CLIPPING PLANES
 * ============================================================================ */

/**
 * Compute t_min where the ray enters the visible region (intersection of
 * all clip half-spaces). Returns 0 if no clipping needed.
 * Returns -1 if ray is entirely clipped (never enters visible region).
 */
static double clip_entry(double ox, double oy, double oz,
                         double dx, double dy, double dz,
                         const render_clip_plane_t* clips, int num_clips) {
    double t_min = 0;

    for (int c = 0; c < num_clips; c++) {
        const double* n = clips[c].normal;
        double d = clips[c].d;

        /* Evaluate plane at ray origin: val = n.origin + d */
        double val_at_origin = n[0]*ox + n[1]*oy + n[2]*oz + d;
        /* Rate of change along ray: dval/dt = n.dir */
        double dval_dt = n[0]*dx + n[1]*dy + n[2]*dz;

        if (val_at_origin >= 0) {
            /* Origin is on visible side */
            if (dval_dt >= 0) continue;  /* Stays visible */
            /* Ray exits visible region at some t_exit — not relevant for entry */
            continue;
        } else {
            /* Origin is on clipped side */
            if (dval_dt <= 1e-15) return -1;  /* Never enters visible side */
            /* Ray enters visible side at t = -val/dval */
            double t_enter = -val_at_origin / dval_dt;
            if (t_enter > t_min) t_min = t_enter;
        }
    }

    return t_min;
}

/**
 * Check if a point is on the visible side of all clip planes.
 */
static int clip_visible(double px, double py, double pz,
                        const render_clip_plane_t* clips, int num_clips) {
    for (int c = 0; c < num_clips; c++) {
        const double* n = clips[c].normal;
        double val = n[0]*px + n[1]*py + n[2]*pz + clips[c].d;
        if (val < -1e-8) return 0;
    }
    return 1;
}

/* ============================================================================
 * COLOR MAPPING
 * ============================================================================ */

void render_get_color(int id, render_color_mode_t mode,
                      const render_config_t* cfg,
                      float* r, float* g, float* b) {
    /* Check custom colors first */
    if (cfg->custom_colors) {
        for (int i = 0; i < cfg->num_custom_colors; i++) {
            if (cfg->custom_colors[i].id == id) {
                *r = cfg->custom_colors[i].r;
                *g = cfg->custom_colors[i].g;
                *b = cfg->custom_colors[i].b;
                return;
            }
        }
    }

    /* Void = white */
    if (id <= 0) {
        *r = 1.0f; *g = 1.0f; *b = 1.0f;
        return;
    }

    if (mode == RENDER_COLOR_DENSITY) {
        /* Density mode is handled separately in shade_density() */
        *r = 0.5f; *g = 0.5f; *b = 0.5f;
        return;
    }

    /* Hash-based palette lookup */
    unsigned int h = (unsigned int)id;
    h = ((h >> 16) ^ h) * 0x45d9f3b;
    h = ((h >> 16) ^ h) * 0x45d9f3b;
    h = (h >> 16) ^ h;
    int idx = (int)(h % RENDER_PALETTE_SIZE);

    *r = palette[idx][0];
    *g = palette[idx][1];
    *b = palette[idx][2];
}

static void shade_density(double density, float* r, float* g, float* b) {
    /* white (void) -> yellow -> orange -> red -> dark red */
    double d = fabs(density);
    if (d <= 0) {
        *r = 1.0f; *g = 1.0f; *b = 1.0f;
        return;
    }
    /* Log scale, typical range 0.001 to 20 g/cc */
    double t = (log10(d) + 3.0) / 4.3;  /* maps 0.001..20 to 0..1 */
    t = clamp01(t);

    if (t < 0.33) {
        double s = t / 0.33;
        *r = (float)(1.0);
        *g = (float)(1.0 - 0.3 * s);
        *b = (float)(1.0 - s);
    } else if (t < 0.66) {
        double s = (t - 0.33) / 0.33;
        *r = (float)(1.0);
        *g = (float)(0.7 - 0.5 * s);
        *b = (float)(0.0);
    } else {
        double s = (t - 0.66) / 0.34;
        *r = (float)(1.0 - 0.4 * s);
        *g = (float)(0.2 - 0.2 * s);
        *b = (float)(0.0);
    }
}

/* ============================================================================
 * SHADING
 * ============================================================================ */

static void shade_surface(float base_r, float base_g, float base_b,
                          double nx, double ny, double nz,
                          const double* light_dir,
                          const double* view_dir,
                          const render_config_t* cfg,
                          float* out_r, float* out_g, float* out_b) {
    /* Ensure normal faces toward camera */
    double ndv = nx*view_dir[0] + ny*view_dir[1] + nz*view_dir[2];
    if (ndv > 0) {
        nx = -nx; ny = -ny; nz = -nz;
    }

    /* Lambertian diffuse */
    double ndl = -(nx*light_dir[0] + ny*light_dir[1] + nz*light_dir[2]);
    if (ndl < 0) ndl = 0;

    /* Phong specular */
    double spec = 0;
    if (ndl > 0 && cfg->specular > 0) {
        /* Reflect light about normal: R = 2(N.L)N - L */
        double rx = 2.0 * ndl * nx + light_dir[0];
        double ry = 2.0 * ndl * ny + light_dir[1];
        double rz = 2.0 * ndl * nz + light_dir[2];
        double rdv = -(rx*view_dir[0] + ry*view_dir[1] + rz*view_dir[2]);
        if (rdv > 0) {
            spec = pow(rdv, (double)cfg->shininess);
        }
    }

    float ka = cfg->ambient;
    float kd = cfg->diffuse * (float)ndl;
    float ks = cfg->specular * (float)spec;

    *out_r = clamp01f(base_r * (ka + kd) + ks);
    *out_g = clamp01f(base_g * (ka + kd) + ks);
    *out_b = clamp01f(base_b * (ka + kd) + ks);
}

/* ============================================================================
 * FRAMEBUFFER
 * ============================================================================ */

render_framebuffer_t* render_framebuffer_create(int width, int height, int aux) {
    render_framebuffer_t* fb = calloc(1, sizeof(render_framebuffer_t));
    if (!fb) return NULL;

    fb->width = width;
    fb->height = height;
    size_t npix = (size_t)width * height;

    fb->color = calloc(npix * 3, sizeof(float));
    fb->cell_id = calloc(npix, sizeof(int));
    if (!fb->color || !fb->cell_id) {
        render_framebuffer_free(fb);
        return NULL;
    }

    /* Initialize cell_id to -1 */
    for (size_t i = 0; i < npix; i++) fb->cell_id[i] = -1;

    if (aux) {
        fb->material_id = calloc(npix, sizeof(int));
        fb->depth = calloc(npix, sizeof(float));
        fb->normal = calloc(npix * 3, sizeof(float));
        if (!fb->material_id || !fb->depth || !fb->normal) {
            render_framebuffer_free(fb);
            return NULL;
        }
    }

    return fb;
}

void render_framebuffer_free(render_framebuffer_t* fb) {
    if (!fb) return;
    free(fb->color);
    free(fb->cell_id);
    free(fb->material_id);
    free(fb->depth);
    free(fb->normal);
    free(fb);
}

/* ============================================================================
 * PER-PIXEL RENDERING (SOLID MODE)
 * ============================================================================ */

/**
 * Render a single pixel in solid/cutaway mode.
 * Uses a thread-local alea_raycast_result_t for buffer reuse.
 */
static void render_pixel_solid(alea_system_t* sys,
                               const render_config_t* cfg,
                               const render_camera_t* cam,
                               double px, double py,
                               int img_width, int img_height,
                               alea_raycast_result_t* result,
                               float* out_color,
                               int* out_cell_id,
                               int* out_material_id,
                               float* out_depth,
                               float* out_normal) {
    double ox, oy, oz, dx, dy, dz;
    render_camera_ray(cam, img_width, img_height, px, py,
                      &ox, &oy, &oz, &dx, &dy, &dz);

    /* Default: background */
    out_color[0] = cfg->background[0];
    out_color[1] = cfg->background[1];
    out_color[2] = cfg->background[2];
    *out_cell_id = -1;
    *out_material_id = 0;
    *out_depth = 1e30f;
    if (out_normal) { out_normal[0] = 0; out_normal[1] = 0; out_normal[2] = 0; }

    /* Compute clip entry */
    double t_min = clip_entry(ox, oy, oz, dx, dy, dz,
                              cfg->clips, cfg->num_clips);
    if (t_min < 0) return;  /* Entirely clipped */

    double t_max = 1e15;

    /* Light and view direction */
    double light_dir[3];
    if (cfg->light_set) {
        light_dir[0] = cfg->light_dir[0];
        light_dir[1] = cfg->light_dir[1];
        light_dir[2] = cfg->light_dir[2];
    } else {
        /* Headlight: light from camera direction */
        light_dir[0] = dx;
        light_dir[1] = dy;
        light_dir[2] = dz;
    }
    double view_dir[3] = {dx, dy, dz};

    /* Check cross-section at clip entry */
    if (t_min > 1e-8) {
        double cpx = ox + (t_min + 1e-6) * dx;
        double cpy = oy + (t_min + 1e-6) * dy;
        double cpz = oz + (t_min + 1e-6) * dz;

        int cs_cell_id = -1, cs_material_id = 0;
        alea_find_cell_at(sys, cpx, cpy, cpz, &cs_cell_id, &cs_material_id);

        if (cs_cell_id >= 0 && cs_material_id != 0) {
            /* Cross-section face: flat shade (skip void material) */
            float br, bg, bb;
            if (cfg->color_mode == RENDER_COLOR_DENSITY) {
                int ci = alea_cell_find(sys, cs_cell_id);
                double density = 0;
                if (ci >= 0) alea_cell_get(sys, ci, NULL, NULL, &density, NULL, NULL, NULL);
                shade_density(density, &br, &bg, &bb);
            } else {
                int color_id = (cfg->color_mode == RENDER_COLOR_CELL) ? cs_cell_id : cs_material_id;
                render_get_color(color_id, cfg->color_mode, cfg, &br, &bg, &bb);
            }
            float tint = cfg->cross_section_tint;
            out_color[0] = br * tint;
            out_color[1] = bg * tint;
            out_color[2] = bb * tint;
            *out_cell_id = cs_cell_id;
            *out_material_id = cs_material_id;
            *out_depth = (float)t_min;
            return;
        }
    }

    /* Full raycast from origin */
    alea_raycast_result_clear(result);

    /* Solid rendering needs exactly one visible interval.  The hierarchical
     * stepper handles lattice DDA and fill expansion, so it can stop as soon
     * as that interval is found instead of materializing a full global trace. */
    alea_ray_t ray;
    alea_ray_init_normalized(&ray, ox, oy, oz, dx, dy, dz);
    alea_ray_first_visible_result_t visible;
    if (alea_raycast_hier_first_visible_nocache(
            sys, &ray, t_min, t_max, -1, 1, result, &visible) != 0 ||
        !visible.found)
        return;

    alea_ray_segment_t seg = {
        .t_enter = visible.t,
        .t_exit = t_max,
        .cell_id = visible.cell_id,
        .material_id = visible.material_id,
        .density = visible.density,
        .enter_surface_id = visible.surface_id,
        .exit_surface_id = -1,
        .enter_hit_index = -1,
        .resolution_flags = visible.resolution_flags,
        .path_index = UINT32_MAX
    };
    if (visible.surface_id > 0) {
        const alea_ray_hit_t hit = {
            .t = visible.t,
            .surface_id = visible.surface_id,
            .primitive_id = visible.primitive_id,
            .nx = visible.nx, .ny = visible.ny, .nz = visible.nz
        };
        if (alea_vec_push(&result->hits, hit, alea_ray_hit_t) != 0)
            return;
        seg.enter_hit_index = 0;
    }
    if (alea_vec_push(&result->segments, seg, alea_ray_segment_t) != 0)
        return;
    /* Find first non-void segment in visible region (past clips) */
    for (size_t i = 0; i < result->segments.count; i++) {
        const alea_ray_segment_t* seg = &result->segments.data[i];

        /* Skip void segments and void-material cells (graveyard etc.) */
        if (seg->cell_id < 0 || seg->material_id == 0) continue;

        /* Skip segments entirely before clip */
        if (seg->t_exit <= t_min + 1e-8) continue;

        /* Determine hit point */
        double t_hit;
        int is_cross_section = 0;

        if (seg->t_enter < t_min + 1e-6) {
            /* Segment straddles clip boundary — cross-section face */
            t_hit = t_min;
            is_cross_section = 1;
        } else {
            t_hit = seg->t_enter;
        }

        /* Get base color */
        float br, bg, bb;
        if (cfg->color_mode == RENDER_COLOR_DENSITY) {
            shade_density(seg->density, &br, &bg, &bb);
        } else {
            int color_id;
            if (cfg->color_mode == RENDER_COLOR_CELL)
                color_id = seg->cell_id;
            else
                color_id = seg->material_id;
            render_get_color(color_id, cfg->color_mode, cfg, &br, &bg, &bb);
        }

        if (is_cross_section) {
            /* Flat shade for cross-section */
            float tint = cfg->cross_section_tint;
            out_color[0] = br * tint;
            out_color[1] = bg * tint;
            out_color[2] = bb * tint;
        } else {
            /* Find surface normal from hit list (O(1) via enter_hit_index) */
            double nx = 0, ny = 0, nz = 1;
            if (seg->enter_hit_index >= 0 &&
                (size_t)seg->enter_hit_index < result->hits.count) {
                nx = result->hits.data[seg->enter_hit_index].nx;
                ny = result->hits.data[seg->enter_hit_index].ny;
                nz = result->hits.data[seg->enter_hit_index].nz;
            } else {
                /* Fallback: linear search (cross-section or cell-aware path) */
                for (size_t h = 0; h < result->hits.count; h++) {
                    if (fabs(result->hits.data[h].t - t_hit) < 1e-6) {
                        nx = result->hits.data[h].nx;
                        ny = result->hits.data[h].ny;
                        nz = result->hits.data[h].nz;
                        break;
                    }
                }
            }

            /* Shade with lighting */
            shade_surface(br, bg, bb, nx, ny, nz,
                         light_dir, view_dir, cfg,
                         &out_color[0], &out_color[1], &out_color[2]);

            /* Shadow ray */
            if (cfg->shadows) {
                double hp[3] = {
                    ox + t_hit * dx,
                    oy + t_hit * dy,
                    oz + t_hit * dz
                };

                /* Offset hit point along normal to avoid self-intersection */
                double bias = 1e-4;
                double sp_ox = hp[0] - nx * bias;
                double sp_oy = hp[1] - ny * bias;
                double sp_oz = hp[2] - nz * bias;

                /* Shadow ray toward light */
                double sl_dx = -light_dir[0];
                double sl_dy = -light_dir[1];
                double sl_dz = -light_dir[2];

                /* Check if shadow point is in visible region (clip planes) */
                int shadow_visible = clip_visible(sp_ox, sp_oy, sp_oz,
                                                   cfg->clips, cfg->num_clips);

                if (shadow_visible) {
                    /* Quick occlusion test (early-out, no segment building) */
                    int occluded = alea_ray_is_occluded(sys,
                        sp_ox, sp_oy, sp_oz, sl_dx, sl_dy, sl_dz,
                        1e6);

                    if (occluded) {
                        /* In shadow — darken */
                        float shadow_factor = 0.5f;
                        out_color[0] *= shadow_factor;
                        out_color[1] *= shadow_factor;
                        out_color[2] *= shadow_factor;
                    }
                }
            }

            /* Store normal */
            if (out_normal) {
                out_normal[0] = (float)nx;
                out_normal[1] = (float)ny;
                out_normal[2] = (float)nz;
            }
        }

        *out_cell_id = seg->cell_id;
        *out_material_id = seg->material_id;
        *out_depth = (float)t_hit;
        return;
    }
}

/* ============================================================================
 * PER-PIXEL RENDERING (X-RAY MODE)
 * ============================================================================ */

typedef struct {
    const render_config_t* cfg;
    float accum_r;
    float accum_g;
    float accum_b;
    float accum_alpha;
    int cell_id;
} render_xray_accumulator_t;

static int render_xray_accumulate_segment(
    void* context, const alea_ray_segment_t* seg) {
    render_xray_accumulator_t* accum = context;
    if (seg->cell_id < 0 || seg->material_id == 0) return 0;
    const double len = seg->t_exit - seg->t_enter;
    if (len <= 0 || len > 1e10) return 0;

    float alpha = (float)(fabs(seg->density) * len *
                          accum->cfg->xray_density_scale);
    if (alpha > 1.0f) alpha = 1.0f;

    float cr, cg, cb;
    const int color_id = accum->cfg->color_mode == RENDER_COLOR_CELL
        ? seg->cell_id : seg->material_id;
    render_get_color(color_id, accum->cfg->color_mode, accum->cfg,
                     &cr, &cg, &cb);
    const float factor = alpha * (1.0f - accum->accum_alpha);
    accum->accum_r += cr * factor;
    accum->accum_g += cg * factor;
    accum->accum_b += cb * factor;
    accum->accum_alpha += factor;
    if (accum->accum_alpha > 0.99f) return 1;
    if (accum->cell_id < 0) accum->cell_id = seg->cell_id;
    return 0;
}

static void render_pixel_xray(alea_system_t* sys,
                               const render_config_t* cfg,
                               const render_camera_t* cam,
                               double px, double py,
                               int img_width, int img_height,
                               alea_raycast_result_t* result,
                               float* out_color,
                               int* out_cell_id) {
    double ox, oy, oz, dx, dy, dz;
    render_camera_ray(cam, img_width, img_height, px, py,
                      &ox, &oy, &oz, &dx, &dy, &dz);

    out_color[0] = cfg->background[0];
    out_color[1] = cfg->background[1];
    out_color[2] = cfg->background[2];
    *out_cell_id = -1;

    /* X-ray consumes selected material intervals directly.  The callback
     * adapter handles fills and lattice DDA without publishing a per-pixel
     * segment vector. */
    render_xray_accumulator_t accum = { .cfg = cfg, .cell_id = -1 };
    alea_ray_t ray;
    if (alea_ray_init(&ray, ox, oy, oz, dx, dy, dz) != 0 ||
        alea_raycast_hier_visit_segments_nocache(
            sys, &ray, 0, result, render_xray_accumulate_segment, &accum) != 0)
        return;

    /* Composite over background */
    const float bg_factor = 1.0f - accum.accum_alpha;
    out_color[0] = accum.accum_r + cfg->background[0] * bg_factor;
    out_color[1] = accum.accum_g + cfg->background[1] * bg_factor;
    out_color[2] = accum.accum_b + cfg->background[2] * bg_factor;
    *out_cell_id = accum.cell_id;
}

/* Compact X-ray tiles are scheduled outside the batch executor's OpenMP
 * region.  That keeps one parallel level while avoiding per-pixel result
 * vectors and nested parallel tracing.  Lattice and AA paths intentionally
 * retain the canonical scalar renderer below. */
static int render_scene_xray_batched(alea_system_t* sys,
                                     const render_config_t* cfg,
                                     const render_camera_t* cam,
                                     render_framebuffer_t* fb,
                                     int tile) {
    const int w = fb->width, h = fb->height;
    const int tiles_x = (w + tile - 1) / tile;
    const int tiles_y = (h + tile - 1) / tile;
    const int n_tiles = tiles_x * tiles_y;
    const alea_raycast_batch_options_t options = {
        .struct_size = sizeof(options),
        .fields = ALEA_RAY_BATCH_MATERIAL | ALEA_RAY_BATCH_DENSITY
    };

    for (int tile_index = 0; tile_index < n_tiles; tile_index++) {
        const int tx = tile_index % tiles_x, ty = tile_index / tiles_x;
        const int x0 = tx * tile, y0 = ty * tile;
        const int x1 = (x0 + tile < w) ? x0 + tile : w;
        const int y1 = (y0 + tile < h) ? y0 + tile : h;
        const size_t ray_count = (size_t)(x1 - x0) * (size_t)(y1 - y0);
        double* origins = malloc(ray_count * 3 * sizeof(*origins));
        double* directions = malloc(ray_count * 3 * sizeof(*directions));
        alea_raycast_batch_result_t* batch = alea_raycast_batch_result_create();
        if (!origins || !directions || !batch) {
            free(origins); free(directions);
            alea_raycast_batch_result_destroy(batch);
            return -1;
        }
        size_t ray_index = 0;
        for (int y = y0; y < y1; y++) {
            for (int x = x0; x < x1; x++, ray_index++) {
                render_camera_ray(cam, w, h, x + 0.5, y + 0.5,
                                  &origins[ray_index * 3],
                                  &origins[ray_index * 3 + 1],
                                  &origins[ray_index * 3 + 2],
                                  &directions[ray_index * 3],
                                  &directions[ray_index * 3 + 1],
                                  &directions[ray_index * 3 + 2]);
            }
        }
        const int batch_rc = alea_raycast_hier_batch(
            sys, origins, directions, ray_count, 0, &options, batch);
        free(origins); free(directions);
        if (batch_rc != 0) {
            /* Retain the established best-effort renderer behavior if a
             * compact allocation or trace fails for this tile. */
            alea_raycast_result_t fallback;
            alea_raycast_result_init(&fallback);
            ray_index = 0;
            for (int y = y0; y < y1; y++) for (int x = x0; x < x1; x++, ray_index++) {
                const size_t pixel = (size_t)y * w + x;
                float color[3]; int cell_id;
                render_pixel_xray(sys, cfg, cam, x + 0.5, y + 0.5, w, h,
                                  &fallback, color, &cell_id);
                fb->color[pixel * 3] = color[0];
                fb->color[pixel * 3 + 1] = color[1];
                fb->color[pixel * 3 + 2] = color[2];
                fb->cell_id[pixel] = cell_id;
            }
            alea_raycast_result_free(&fallback);
            alea_raycast_batch_result_destroy(batch);
            continue;
        }
        const uint64_t* offsets = alea_raycast_batch_ray_offsets(batch);
        const double* enters = alea_raycast_batch_t_enter(batch);
        const double* exits = alea_raycast_batch_t_exit(batch);
        const int32_t* cells = alea_raycast_batch_cell_ids(batch);
        const int32_t* materials = alea_raycast_batch_material_ids(batch);
        const double* densities = alea_raycast_batch_densities(batch);
        if (!offsets || !enters || !exits || !cells || !materials || !densities) {
            alea_raycast_batch_result_destroy(batch);
            return -1;
        }
        ray_index = 0;
        for (int y = y0; y < y1; y++) for (int x = x0; x < x1; x++, ray_index++) {
            float accum_r = 0, accum_g = 0, accum_b = 0, accum_alpha = 0;
            int cell_id = -1;
            for (size_t i = offsets[ray_index]; i < offsets[ray_index + 1]; i++) {
                if (cells[i] < 0 || materials[i] == 0) continue;
                const double len = exits[i] - enters[i];
                if (len <= 0 || len > 1e10) continue;
                float alpha = (float)(fabs(densities[i]) * len * cfg->xray_density_scale);
                if (alpha > 1.0f) alpha = 1.0f;
                float cr, cg, cb;
                const int color_id = cfg->color_mode == RENDER_COLOR_CELL ?
                    cells[i] : materials[i];
                render_get_color(color_id, cfg->color_mode, cfg, &cr, &cg, &cb);
                const float factor = alpha * (1.0f - accum_alpha);
                accum_r += cr * factor; accum_g += cg * factor; accum_b += cb * factor;
                accum_alpha += factor;
                if (accum_alpha > 0.99f) break;
                if (cell_id < 0) cell_id = cells[i];
            }
            const size_t pixel = (size_t)y * w + x;
            const float bg_factor = 1.0f - accum_alpha;
            fb->color[pixel * 3] = accum_r + cfg->background[0] * bg_factor;
            fb->color[pixel * 3 + 1] = accum_g + cfg->background[1] * bg_factor;
            fb->color[pixel * 3 + 2] = accum_b + cfg->background[2] * bg_factor;
            fb->cell_id[pixel] = cell_id;
            if (fb->material_id) fb->material_id[pixel] = 0;
            if (fb->depth) fb->depth[pixel] = 1e30f;
            if (fb->normal) {
                fb->normal[pixel * 3] = 0;
                fb->normal[pixel * 3 + 1] = 0;
                fb->normal[pixel * 3 + 2] = 0;
            }
        }
        alea_raycast_batch_result_destroy(batch);
        fprintf(stderr, "\rrender: %d/%d tiles (%.0f%%)", tile_index + 1, n_tiles,
                100.0 * (tile_index + 1) / n_tiles);
    }
    fprintf(stderr, "\rrender: %d/%d tiles (100%%)\n", n_tiles, n_tiles);
    return 0;
}

/* ============================================================================
 * TILE-BASED PARALLEL RENDERING
 * ============================================================================ */

int render_scene(alea_system_t* sys,
                 const render_config_t* cfg,
                 const render_camera_t* cam,
                 render_framebuffer_t* fb) {
    if (!sys || !cfg || !cam || !fb) return -1;

    int w = fb->width;
    int h = fb->height;
    int tile = cfg->tile_size > 0 ? cfg->tile_size : RENDER_DEFAULT_TILE;
    int tiles_x = (w + tile - 1) / tile;
    int tiles_y = (h + tile - 1) / tile;
    int n_tiles = tiles_x * tiles_y;
    int aa = cfg->aa_samples > 0 ? cfg->aa_samples : 1;

    /* Ensure caches are built before parallel section */
    if (alea_raycast_ensure_caches(sys) != 0) return -1;

    int num_threads = cfg->threads;
#ifdef _OPENMP
    if (num_threads <= 0) num_threads = omp_get_max_threads();
    omp_set_num_threads(num_threads);
#else
    num_threads = 1;
#endif

    fprintf(stderr, "render: %dx%d, %d tiles (%dx%d), %d thread%s, aa=%dx%d\n",
            w, h, n_tiles, tile, tile, num_threads,
            num_threads > 1 ? "s" : "", aa, aa);

    if (cfg->render_mode == RENDER_MODE_XRAY && !sys->has_lattice && aa == 1)
        return render_scene_xray_batched(sys, cfg, cam, fb, tile);

    int progress_done = 0;

    #pragma omp parallel
    {
        /* Thread-local raycast result (reused across pixels) */
        alea_raycast_result_t result;
        alea_raycast_result_init(&result);
        alea_raycast_result_reserve(&result, 64, 32);

        #pragma omp for schedule(dynamic, 1)
        for (int t = 0; t < n_tiles; t++) {
            int tx = t % tiles_x;
            int ty = t / tiles_x;
            int x0 = tx * tile;
            int y0 = ty * tile;
            int x1 = x0 + tile; if (x1 > w) x1 = w;
            int y1 = y0 + tile; if (y1 > h) y1 = h;

            for (int y = y0; y < y1; y++) {
                for (int x = x0; x < x1; x++) {
                    size_t pidx = (size_t)y * w + x;

                    if (aa <= 1) {
                        /* Single sample at pixel center */
                        float color[3];
                        int cell_id = -1, mat_id = 0;
                        float depth = 1e30f;
                        float normal[3] = {0};

                        switch (cfg->render_mode) {
                        case RENDER_MODE_XRAY:
                            render_pixel_xray(sys, cfg, cam,
                                x + 0.5, y + 0.5, w, h,
                                &result, color, &cell_id);
                            break;
                        case RENDER_MODE_DEPTH:
                        case RENDER_MODE_CELLID:
                        case RENDER_MODE_MATID:
                        case RENDER_MODE_SOLID:
                        default:
                            render_pixel_solid(sys, cfg, cam,
                                x + 0.5, y + 0.5, w, h,
                                &result, color, &cell_id, &mat_id,
                                &depth, normal);
                            break;
                        }

                        fb->color[pidx * 3 + 0] = color[0];
                        fb->color[pidx * 3 + 1] = color[1];
                        fb->color[pidx * 3 + 2] = color[2];
                        fb->cell_id[pidx] = cell_id;
                        if (fb->material_id) fb->material_id[pidx] = mat_id;
                        if (fb->depth) fb->depth[pidx] = depth;
                        if (fb->normal) {
                            fb->normal[pidx * 3 + 0] = normal[0];
                            fb->normal[pidx * 3 + 1] = normal[1];
                            fb->normal[pidx * 3 + 2] = normal[2];
                        }
                    } else {
                        /* Supersampling NxN */
                        float acc_r = 0, acc_g = 0, acc_b = 0;
                        int center_cell = -1, center_mat = 0;
                        float center_depth = 1e30f;
                        float center_normal[3] = {0};
                        float inv_samples = 1.0f / (float)(aa * aa);

                        for (int sy = 0; sy < aa; sy++) {
                            for (int sx = 0; sx < aa; sx++) {
                                double spx = x + (sx + 0.5) / aa;
                                double spy = y + (sy + 0.5) / aa;

                                float color[3];
                                int cell_id = -1, mat_id = 0;
                                float depth = 1e30f;
                                float normal[3] = {0};

                                switch (cfg->render_mode) {
                                case RENDER_MODE_XRAY:
                                    render_pixel_xray(sys, cfg, cam,
                                        spx, spy, w, h,
                                        &result, color, &cell_id);
                                    break;
                                default:
                                    render_pixel_solid(sys, cfg, cam,
                                        spx, spy, w, h,
                                        &result, color, &cell_id, &mat_id,
                                        &depth, normal);
                                    break;
                                }

                                acc_r += color[0];
                                acc_g += color[1];
                                acc_b += color[2];

                                /* Use center sample for ID/depth/normal */
                                if (sx == aa/2 && sy == aa/2) {
                                    center_cell = cell_id;
                                    center_mat = mat_id;
                                    center_depth = depth;
                                    center_normal[0] = normal[0];
                                    center_normal[1] = normal[1];
                                    center_normal[2] = normal[2];
                                }
                            }
                        }

                        fb->color[pidx * 3 + 0] = acc_r * inv_samples;
                        fb->color[pidx * 3 + 1] = acc_g * inv_samples;
                        fb->color[pidx * 3 + 2] = acc_b * inv_samples;
                        fb->cell_id[pidx] = center_cell;
                        if (fb->material_id) fb->material_id[pidx] = center_mat;
                        if (fb->depth) fb->depth[pidx] = center_depth;
                        if (fb->normal) {
                            fb->normal[pidx * 3 + 0] = center_normal[0];
                            fb->normal[pidx * 3 + 1] = center_normal[1];
                            fb->normal[pidx * 3 + 2] = center_normal[2];
                        }
                    }
                }
            }

            /* Progress reporting */
            #pragma omp atomic
            progress_done++;

            if (progress_done % (n_tiles / 20 + 1) == 0) {
                #pragma omp critical
                {
                    fprintf(stderr, "\rrender: %d/%d tiles (%.0f%%)",
                            progress_done, n_tiles,
                            100.0 * progress_done / n_tiles);
                }
            }
        }

        alea_raycast_result_free(&result);
    }

    fprintf(stderr, "\rrender: %d/%d tiles (100%%)\n", n_tiles, n_tiles);

    return 0;
}

/* ============================================================================
 * POST-PROCESSING: EDGE DARKENING
 * ============================================================================ */

void render_edge_darken(render_framebuffer_t* fb) {
    if (!fb || !fb->cell_id) return;
    int w = fb->width;
    int h = fb->height;

    /* Detect cell ID discontinuities in 3x3 kernel and darken */
    for (int y = 1; y < h - 1; y++) {
        for (int x = 1; x < w - 1; x++) {
            int center = fb->cell_id[y * w + x];
            int is_edge = 0;

            for (int dy = -1; dy <= 1 && !is_edge; dy++) {
                for (int dx = -1; dx <= 1 && !is_edge; dx++) {
                    if (dx == 0 && dy == 0) continue;
                    int neighbor = fb->cell_id[(y + dy) * w + (x + dx)];
                    if (neighbor != center) is_edge = 1;
                }
            }

            if (is_edge) {
                size_t idx = (size_t)(y * w + x) * 3;
                fb->color[idx + 0] *= 0.3f;
                fb->color[idx + 1] *= 0.3f;
                fb->color[idx + 2] *= 0.3f;
            }
        }
    }
}

/* ============================================================================
 * TONEMAP (float -> uint8)
 * ============================================================================ */

void render_tonemap(const render_framebuffer_t* fb, uint8_t* pixels) {
    size_t npix = (size_t)fb->width * fb->height;
    for (size_t i = 0; i < npix; i++) {
        pixels[i * 3 + 0] = (uint8_t)(clamp01f(fb->color[i * 3 + 0]) * 255.0f + 0.5f);
        pixels[i * 3 + 1] = (uint8_t)(clamp01f(fb->color[i * 3 + 1]) * 255.0f + 0.5f);
        pixels[i * 3 + 2] = (uint8_t)(clamp01f(fb->color[i * 3 + 2]) * 255.0f + 0.5f);
    }
}

/* ============================================================================
 * IMAGE WRITERS
 * ============================================================================ */

/* PNG CRC table */
static uint32_t png_crc_table[256];
static int png_crc_table_ready = 0;

static void png_init_crc_table(void) {
    if (png_crc_table_ready) return;
    for (int n = 0; n < 256; n++) {
        uint32_t c = (uint32_t)n;
        for (int k = 0; k < 8; k++) {
            c = (c & 1) ? (0xEDB88320 ^ (c >> 1)) : (c >> 1);
        }
        png_crc_table[n] = c;
    }
    png_crc_table_ready = 1;
}

static uint32_t png_crc32(const uint8_t* data, size_t len) {
    png_init_crc_table();
    uint32_t crc = 0xFFFFFFFF;
    for (size_t i = 0; i < len; i++) {
        crc = png_crc_table[(crc ^ data[i]) & 0xFF] ^ (crc >> 8);
    }
    return crc ^ 0xFFFFFFFF;
}

static uint32_t png_adler32(const uint8_t* data, size_t len) {
    uint32_t a = 1, b = 0;
    for (size_t i = 0; i < len; i++) {
        a = (a + data[i]) % 65521;
        b = (b + a) % 65521;
    }
    return (b << 16) | a;
}

static void png_write_be32(uint8_t* p, uint32_t v) {
    p[0] = (v >> 24) & 0xFF;
    p[1] = (v >> 16) & 0xFF;
    p[2] = (v >> 8) & 0xFF;
    p[3] = v & 0xFF;
}

static void png_write_chunk(FILE* f, const char* type, const uint8_t* data, size_t len) {
    uint8_t header[8];
    png_write_be32(header, (uint32_t)len);
    memcpy(header + 4, type, 4);
    fwrite(header, 1, 8, f);

    if (len > 0 && data) fwrite(data, 1, len, f);

    uint8_t* crc_buf = malloc(4 + len);
    if (crc_buf) {
        memcpy(crc_buf, type, 4);
        if (len > 0 && data) memcpy(crc_buf + 4, data, len);
        uint32_t crc = png_crc32(crc_buf, 4 + len);
        free(crc_buf);
        uint8_t crc_bytes[4];
        png_write_be32(crc_bytes, crc);
        fwrite(crc_bytes, 1, 4, f);
    }
}

int render_write_png(const char* filename, const uint8_t* pixels,
                     int width, int height) {
    FILE* f = fopen(filename, "wb");
    if (!f) {
        alea_set_error_detail(ALEA_ERR_FILE_WRITE,
                              "failed to open render output file '%s': %s",
                              filename, strerror(errno));
        return -1;
    }

    static const uint8_t png_sig[8] = {137, 80, 78, 71, 13, 10, 26, 10};
    fwrite(png_sig, 1, 8, f);

    uint8_t ihdr[13];
    png_write_be32(ihdr, width);
    png_write_be32(ihdr + 4, height);
    ihdr[8] = 8; ihdr[9] = 2; ihdr[10] = 0; ihdr[11] = 0; ihdr[12] = 0;
    png_write_chunk(f, "IHDR", ihdr, 13);

    size_t raw_row = 1 + (size_t)width * 3;
    size_t raw_size = raw_row * height;
    uint8_t* raw = malloc(raw_size);
    if (!raw) { fclose(f); return -1; }

    for (int y = 0; y < height; y++) {
        raw[y * raw_row] = 0;
        memcpy(raw + y * raw_row + 1, pixels + y * width * 3, width * 3);
    }

    size_t max_block = 65535;
    size_t num_blocks = (raw_size + max_block - 1) / max_block;
    size_t zlib_size = 2 + raw_size + num_blocks * 5 + 4;
    uint8_t* zlib = malloc(zlib_size);
    if (!zlib) { free(raw); fclose(f); return -1; }

    size_t zp = 0;
    zlib[zp++] = 0x78; zlib[zp++] = 0x01;

    size_t remaining = raw_size, src_pos = 0;
    while (remaining > 0) {
        size_t block_size = (remaining > max_block) ? max_block : remaining;
        int is_final = (remaining <= max_block) ? 1 : 0;
        zlib[zp++] = is_final ? 0x01 : 0x00;
        zlib[zp++] = block_size & 0xFF;
        zlib[zp++] = (block_size >> 8) & 0xFF;
        zlib[zp++] = (~block_size) & 0xFF;
        zlib[zp++] = ((~block_size) >> 8) & 0xFF;
        memcpy(zlib + zp, raw + src_pos, block_size);
        zp += block_size;
        src_pos += block_size;
        remaining -= block_size;
    }

    uint32_t adler = png_adler32(raw, raw_size);
    png_write_be32(zlib + zp, adler);
    zp += 4;

    png_write_chunk(f, "IDAT", zlib, zp);
    free(zlib);
    free(raw);

    png_write_chunk(f, "IEND", NULL, 0);
    fclose(f);
    return 0;
}

int render_write_bmp(const char* filename, const uint8_t* pixels,
                     int width, int height) {
    FILE* f = fopen(filename, "wb");
    if (!f) {
        alea_set_error_detail(ALEA_ERR_FILE_WRITE,
                              "failed to open render output file '%s': %s",
                              filename, strerror(errno));
        return -1;
    }

    int row_size = ((width * 3 + 3) / 4) * 4;
    int data_size = row_size * height;
    int file_size = 54 + data_size;

    uint8_t header[54] = {0};
    header[0] = 'B'; header[1] = 'M';
    header[2] = file_size & 0xFF;
    header[3] = (file_size >> 8) & 0xFF;
    header[4] = (file_size >> 16) & 0xFF;
    header[5] = (file_size >> 24) & 0xFF;
    header[10] = 54; header[14] = 40;
    header[18] = width & 0xFF; header[19] = (width >> 8) & 0xFF;
    header[20] = (width >> 16) & 0xFF; header[21] = (width >> 24) & 0xFF;
    header[22] = height & 0xFF; header[23] = (height >> 8) & 0xFF;
    header[24] = (height >> 16) & 0xFF; header[25] = (height >> 24) & 0xFF;
    header[26] = 1; header[28] = 24;

    fwrite(header, 1, 54, f);

    uint8_t* row = malloc(row_size);
    for (int y = height - 1; y >= 0; y--) {
        memset(row, 0, row_size);
        for (int x = 0; x < width; x++) {
            int src = (y * width + x) * 3;
            int dst = x * 3;
            row[dst + 0] = pixels[src + 2];  /* BGR */
            row[dst + 1] = pixels[src + 1];
            row[dst + 2] = pixels[src + 0];
        }
        fwrite(row, 1, row_size, f);
    }

    free(row);
    fclose(f);
    return 0;
}

int render_write_ppm(const char* filename, const uint8_t* pixels,
                     int width, int height) {
    FILE* f = fopen(filename, "wb");
    if (!f) {
        alea_set_error_detail(ALEA_ERR_FILE_WRITE,
                              "failed to open render output file '%s': %s",
                              filename, strerror(errno));
        return -1;
    }

    fprintf(f, "P6\n%d %d\n255\n", width, height);
    fwrite(pixels, 1, (size_t)width * height * 3, f);
    fclose(f);
    return 0;
}

int render_write_image(const char* filename, const uint8_t* pixels,
                       int width, int height) {
    const char* ext = strrchr(filename, '.');
    if (!ext) return render_write_png(filename, pixels, width, height);

    if (alea_strcasecmp(ext, ".bmp") == 0)
        return render_write_bmp(filename, pixels, width, height);
    if (alea_strcasecmp(ext, ".ppm") == 0)
        return render_write_ppm(filename, pixels, width, height);
    return render_write_png(filename, pixels, width, height);
}

/* ============================================================================
 * AUXILIARY OUTPUT
 * ============================================================================ */

int render_write_aux(const char* base_filename, const render_framebuffer_t* fb) {
    if (!fb) return -1;
    int w = fb->width, h = fb->height;
    size_t npix = (size_t)w * h;

    /* Allocate pixel buffer */
    uint8_t* pixels = malloc(npix * 3);
    if (!pixels) return -1;

    /* Strip extension from base filename */
    char base[512];
    strncpy(base, base_filename, sizeof(base) - 1);
    base[sizeof(base) - 1] = '\0';
    char* dot = strrchr(base, '.');
    if (dot) *dot = '\0';

    char path[600];

    /* Depth map */
    if (fb->depth) {
        float d_min = 1e30f, d_max = 0;
        for (size_t i = 0; i < npix; i++) {
            if (fb->cell_id[i] >= 0) {
                if (fb->depth[i] < d_min) d_min = fb->depth[i];
                if (fb->depth[i] > d_max) d_max = fb->depth[i];
            }
        }
        float range = (d_max > d_min) ? (d_max - d_min) : 1.0f;

        for (size_t i = 0; i < npix; i++) {
            uint8_t v = 0;
            if (fb->cell_id[i] >= 0) {
                float t = 1.0f - (fb->depth[i] - d_min) / range;
                v = (uint8_t)(clamp01f(t) * 255.0f + 0.5f);
            }
            pixels[i * 3 + 0] = v;
            pixels[i * 3 + 1] = v;
            pixels[i * 3 + 2] = v;
        }
        snprintf(path, sizeof(path), "%s_depth.png", base);
        render_write_png(path, pixels, w, h);
    }

    /* Cell ID map */
    if (fb->cell_id) {
        for (size_t i = 0; i < npix; i++) {
            int id = fb->cell_id[i];
            if (id < 0) {
                pixels[i * 3 + 0] = 0;
                pixels[i * 3 + 1] = 0;
                pixels[i * 3 + 2] = 0;
            } else {
                /* RGB-encode integer */
                pixels[i * 3 + 0] = id & 0xFF;
                pixels[i * 3 + 1] = (id >> 8) & 0xFF;
                pixels[i * 3 + 2] = (id >> 16) & 0xFF;
            }
        }
        snprintf(path, sizeof(path), "%s_cellid.png", base);
        render_write_png(path, pixels, w, h);
    }

    /* Material ID map */
    if (fb->material_id) {
        for (size_t i = 0; i < npix; i++) {
            int id = fb->material_id[i];
            if (fb->cell_id[i] < 0) {
                pixels[i * 3 + 0] = 0;
                pixels[i * 3 + 1] = 0;
                pixels[i * 3 + 2] = 0;
            } else {
                pixels[i * 3 + 0] = id & 0xFF;
                pixels[i * 3 + 1] = (id >> 8) & 0xFF;
                pixels[i * 3 + 2] = (id >> 16) & 0xFF;
            }
        }
        snprintf(path, sizeof(path), "%s_matid.png", base);
        render_write_png(path, pixels, w, h);
    }

    /* Normal map */
    if (fb->normal) {
        for (size_t i = 0; i < npix; i++) {
            if (fb->cell_id[i] < 0) {
                pixels[i * 3 + 0] = 128;
                pixels[i * 3 + 1] = 128;
                pixels[i * 3 + 2] = 255;
            } else {
                pixels[i * 3 + 0] = (uint8_t)((fb->normal[i * 3 + 0] * 0.5f + 0.5f) * 255.0f);
                pixels[i * 3 + 1] = (uint8_t)((fb->normal[i * 3 + 1] * 0.5f + 0.5f) * 255.0f);
                pixels[i * 3 + 2] = (uint8_t)((fb->normal[i * 3 + 2] * 0.5f + 0.5f) * 255.0f);
            }
        }
        snprintf(path, sizeof(path), "%s_normal.png", base);
        render_write_png(path, pixels, w, h);
    }

    free(pixels);
    return 0;
}

/* ============================================================================
 * CUSTOM COLOR LOADING
 * ============================================================================ */

int render_load_colors(const char* filename, render_config_t* cfg) {
    FILE* f = fopen(filename, "r");
    if (!f) return -1;

    int capacity = 64;
    cfg->custom_colors = malloc(capacity * sizeof(render_color_entry_t));
    cfg->num_custom_colors = 0;

    char line[256];
    while (fgets(line, sizeof(line), f)) {
        if (line[0] == '#' || line[0] == '\n') continue;

        int id;
        float r, g, b;
        if (sscanf(line, "%d %f %f %f", &id, &r, &g, &b) == 4) {
            if (cfg->num_custom_colors >= capacity) {
                capacity *= 2;
                cfg->custom_colors = realloc(cfg->custom_colors,
                                             capacity * sizeof(render_color_entry_t));
            }
            render_color_entry_t* e = &cfg->custom_colors[cfg->num_custom_colors++];
            e->id = id;
            e->r = r;
            e->g = g;
            e->b = b;
        }
    }

    fclose(f);
    return 0;
}
