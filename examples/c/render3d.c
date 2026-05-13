// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/**
 * @file render3d.c
 * @brief 3D batch renderer for Alea CSG geometry
 *
 * Produces publication-quality 3D images of CSG models.
 *
 * Usage:
 *   render3d model.i -o render.png --width 1920 --height 1080 \
 *       --eye 100,50,200 --target 0,0,0 --fov 45 \
 *       --clip-x=0 --clip-z=-50 \
 *       --color material --light 1,1,2 --shadows \
 *       --aa 2 --threads 8
 *
 * Build:
 *   gcc -O2 -fopenmp -o render3d render3d.c -I../include -I../src \
 *       -L../bin -lalea_render -lalea_raycast \
 *       -lalea_mcnp -lalea_openmc -lalea -lm
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>

#ifdef _WIN32
#include <windows.h>
#define strcasecmp _stricmp
static double get_time_ms(void) {
    LARGE_INTEGER freq, count;
    QueryPerformanceFrequency(&freq);
    QueryPerformanceCounter(&count);
    return (double)count.QuadPart / freq.QuadPart * 1000.0;
}
#else
#include <sys/time.h>
#include <strings.h>
static double get_time_ms(void) {
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return tv.tv_sec * 1000.0 + tv.tv_usec / 1000.0;
}
#endif

#include "alea.h"
#include "alea_mcnp.h"
#include "alea_openmc.h"
#include "alea_raycast.h"
#include "alea_render.h"

/* ============================================================================
 * ARGUMENT PARSING HELPERS
 * ============================================================================ */

static int parse_vec3(const char* s, double* x, double* y, double* z) {
    return sscanf(s, "%lf,%lf,%lf", x, y, z) == 3 ? 0 : -1;
}

static int parse_color3(const char* s, float* r, float* g, float* b) {
    return sscanf(s, "%f,%f,%f", r, g, b) == 3 ? 0 : -1;
}

static int ends_with(const char* str, const char* suffix) {
    size_t sl = strlen(str), xl = strlen(suffix);
    if (sl < xl) return 0;
    return strcasecmp(str + sl - xl, suffix) == 0;
}

static void print_usage(const char* prog) {
    fprintf(stderr,
"Usage: %s <input> [options]\n"
"\n"
"INPUT:\n"
"  input.i / input.inp          MCNP geometry file\n"
"  input.xml                    OpenMC geometry file\n"
"\n"
"CAMERA:\n"
"  --eye X,Y,Z                  Camera position (default: auto-fit)\n"
"  --target X,Y,Z               Look-at point (default: model center)\n"
"  --up X,Y,Z                   Up vector (default: 0,0,1)\n"
"  --fov DEGREES                Vertical FOV (default: 45; 0=orthographic)\n"
"  --ortho-height H             Ortho view height in model units\n"
"\n"
"IMAGE:\n"
"  -o, --output FILE            Output file (default: render.png)\n"
"  -w, --width N                Image width (default: 1920)\n"
"  --height N                   Image height (default: 1080)\n"
"  --background R,G,B           Background color (default: 0.95,0.95,0.95)\n"
"\n"
"CUTAWAY:\n"
"  --clip-x V                   Show only X > V\n"
"  --clip-y V                   Show only Y > V\n"
"  --clip-z V                   Show only Z > V\n"
"  --clip-nx V                  Show only X < V (reversed)\n"
"  --clip-ny V                  Show only Y < V\n"
"  --clip-nz V                  Show only Z < V\n"
"  --clip A,B,C,D               Arbitrary clip plane (visible where ax+by+cz+d > 0)\n"
"\n"
"APPEARANCE:\n"
"  --color MODE                 Coloring: material (default), cell, universe, density\n"
"  --colors FILE                Custom color palette file\n"
"  --light X,Y,Z                Light direction (default: from camera)\n"
"  --shadows                    Enable shadow rays\n"
"  --edges                      Enable cell-boundary edge darkening\n"
"  --ambient F                  Ambient light factor (default: 0.3)\n"
"  --specular F                 Specular intensity (default: 0.2)\n"
"  --cross-section-tint F       Darken cross-section faces (default: 0.85)\n"
"\n"
"RENDERING MODE:\n"
"  --mode MODE                  solid (default), xray, depth, cellid, matid\n"
"  --xray-scale F               X-ray density scale factor (default: 0.1)\n"
"\n"
"QUALITY:\n"
"  --aa N                       Supersampling NxN (default: 1)\n"
"  --universe-depth N           -1=innermost (default), 0=root, N=level\n"
"\n"
"PERFORMANCE:\n"
"  --threads N                  Thread count (default: all cores)\n"
"  --tile-size N                Tile size in pixels (default: 32)\n"
"  --preview                    Write quick 1/4-res preview first\n"
"\n"
"AUXILIARY:\n"
"  --aux                        Write depth/cellid/matid/normal maps\n"
"  --dedup                      Enable surface deduplication\n"
"  --log-level N                Verbosity (0=quiet, 3=debug)\n"
"\n", prog);
}

static void add_clip_axis(render_config_t* cfg, int axis, double value, int negative) {
    if (cfg->num_clips >= RENDER_MAX_CLIPS) return;
    render_clip_plane_t* cp = &cfg->clips[cfg->num_clips++];
    memset(cp, 0, sizeof(*cp));
    cp->normal[axis] = negative ? -1.0 : 1.0;
    cp->d = negative ? value : -value;
}

/* ============================================================================
 * MAIN
 * ============================================================================ */

int main(int argc, char** argv) {
    if (argc < 2) {
        print_usage(argv[0]);
        return 1;
    }

    const char* input_file = NULL;
    const char* output_file = "render.png";
    render_config_t cfg;
    render_config_init(&cfg);

    /* Parse arguments */
    for (int i = 1; i < argc; i++) {
        const char* arg = argv[i];

        if (arg[0] != '-' && !input_file) {
            input_file = arg;
            continue;
        }

        if (strcmp(arg, "-o") == 0 || strcmp(arg, "--output") == 0) {
            if (++i < argc) output_file = argv[i];
        } else if (strcmp(arg, "-w") == 0 || strcmp(arg, "--width") == 0) {
            if (++i < argc) cfg.width = atoi(argv[i]);
        } else if (strcmp(arg, "--height") == 0) {
            if (++i < argc) cfg.height = atoi(argv[i]);
        } else if (strcmp(arg, "--eye") == 0) {
            if (++i < argc) {
                if (parse_vec3(argv[i], &cfg.eye[0], &cfg.eye[1], &cfg.eye[2]) == 0)
                    cfg.eye_set = 1;
            }
        } else if (strcmp(arg, "--target") == 0) {
            if (++i < argc) {
                if (parse_vec3(argv[i], &cfg.target[0], &cfg.target[1], &cfg.target[2]) == 0)
                    cfg.target_set = 1;
            }
        } else if (strcmp(arg, "--up") == 0) {
            if (++i < argc) parse_vec3(argv[i], &cfg.up[0], &cfg.up[1], &cfg.up[2]);
        } else if (strcmp(arg, "--fov") == 0) {
            if (++i < argc) cfg.fov = atof(argv[i]);
        } else if (strcmp(arg, "--ortho-height") == 0) {
            if (++i < argc) { cfg.ortho_height = atof(argv[i]); cfg.fov = 0; }
        } else if (strcmp(arg, "--background") == 0) {
            if (++i < argc) parse_color3(argv[i], &cfg.background[0], &cfg.background[1], &cfg.background[2]);
        }
        /* Clipping */
        else if (strncmp(arg, "--clip-x=", 9) == 0) {
            add_clip_axis(&cfg, 0, atof(arg + 9), 0);
        } else if (strncmp(arg, "--clip-y=", 9) == 0) {
            add_clip_axis(&cfg, 1, atof(arg + 9), 0);
        } else if (strncmp(arg, "--clip-z=", 9) == 0) {
            add_clip_axis(&cfg, 2, atof(arg + 9), 0);
        } else if (strncmp(arg, "--clip-nx=", 10) == 0) {
            add_clip_axis(&cfg, 0, atof(arg + 10), 1);
        } else if (strncmp(arg, "--clip-ny=", 10) == 0) {
            add_clip_axis(&cfg, 1, atof(arg + 10), 1);
        } else if (strncmp(arg, "--clip-nz=", 10) == 0) {
            add_clip_axis(&cfg, 2, atof(arg + 10), 1);
        } else if (strcmp(arg, "--clip") == 0) {
            if (++i < argc) {
                double a, b, c, d;
                if (sscanf(argv[i], "%lf,%lf,%lf,%lf", &a, &b, &c, &d) == 4) {
                    if (cfg.num_clips < RENDER_MAX_CLIPS) {
                        render_clip_plane_t* cp = &cfg.clips[cfg.num_clips++];
                        cp->normal[0] = a; cp->normal[1] = b; cp->normal[2] = c;
                        cp->d = d;
                    }
                }
            }
        }
        /* Appearance */
        else if (strcmp(arg, "--color") == 0) {
            if (++i < argc) {
                if (strcmp(argv[i], "cell") == 0) cfg.color_mode = RENDER_COLOR_CELL;
                else if (strcmp(argv[i], "universe") == 0) cfg.color_mode = RENDER_COLOR_UNIVERSE;
                else if (strcmp(argv[i], "density") == 0) cfg.color_mode = RENDER_COLOR_DENSITY;
                else cfg.color_mode = RENDER_COLOR_MATERIAL;
            }
        } else if (strcmp(arg, "--colors") == 0) {
            if (++i < argc) render_load_colors(argv[i], &cfg);
        } else if (strcmp(arg, "--light") == 0) {
            if (++i < argc) {
                if (parse_vec3(argv[i], &cfg.light_dir[0], &cfg.light_dir[1], &cfg.light_dir[2]) == 0) {
                    /* Normalize light direction */
                    double len = sqrt(cfg.light_dir[0]*cfg.light_dir[0] +
                                      cfg.light_dir[1]*cfg.light_dir[1] +
                                      cfg.light_dir[2]*cfg.light_dir[2]);
                    if (len > 1e-10) {
                        cfg.light_dir[0] /= len;
                        cfg.light_dir[1] /= len;
                        cfg.light_dir[2] /= len;
                    }
                    cfg.light_set = 1;
                }
            }
        } else if (strcmp(arg, "--shadows") == 0) {
            cfg.shadows = 1;
        } else if (strcmp(arg, "--edges") == 0) {
            cfg.edges = 1;
        } else if (strcmp(arg, "--ambient") == 0) {
            if (++i < argc) cfg.ambient = (float)atof(argv[i]);
        } else if (strcmp(arg, "--specular") == 0) {
            if (++i < argc) cfg.specular = (float)atof(argv[i]);
        } else if (strcmp(arg, "--cross-section-tint") == 0) {
            if (++i < argc) cfg.cross_section_tint = (float)atof(argv[i]);
        }
        /* Rendering mode */
        else if (strcmp(arg, "--mode") == 0) {
            if (++i < argc) {
                if (strcmp(argv[i], "xray") == 0) cfg.render_mode = RENDER_MODE_XRAY;
                else if (strcmp(argv[i], "depth") == 0) cfg.render_mode = RENDER_MODE_DEPTH;
                else if (strcmp(argv[i], "cellid") == 0) cfg.render_mode = RENDER_MODE_CELLID;
                else if (strcmp(argv[i], "matid") == 0) cfg.render_mode = RENDER_MODE_MATID;
                else cfg.render_mode = RENDER_MODE_SOLID;
            }
        } else if (strcmp(arg, "--xray-scale") == 0) {
            if (++i < argc) cfg.xray_density_scale = (float)atof(argv[i]);
        }
        /* Quality */
        else if (strcmp(arg, "--aa") == 0) {
            if (++i < argc) cfg.aa_samples = atoi(argv[i]);
        } else if (strcmp(arg, "--universe-depth") == 0) {
            if (++i < argc) cfg.universe_depth = atoi(argv[i]);
        }
        /* Performance */
        else if (strcmp(arg, "--threads") == 0) {
            if (++i < argc) cfg.threads = atoi(argv[i]);
        } else if (strcmp(arg, "--tile-size") == 0) {
            if (++i < argc) cfg.tile_size = atoi(argv[i]);
        } else if (strcmp(arg, "--preview") == 0) {
            cfg.preview = 1;
        }
        /* Auxiliary */
        else if (strcmp(arg, "--aux") == 0) {
            cfg.aux_output = 1;
        } else if (strcmp(arg, "--dedup") == 0) {
            cfg.dedup = 1;
        } else if (strcmp(arg, "--log-level") == 0) {
            if (++i < argc) cfg.log_level = atoi(argv[i]);
        }
        /* Help */
        else if (strcmp(arg, "--help") == 0 || strcmp(arg, "-h") == 0) {
            print_usage(argv[0]);
            return 0;
        } else {
            fprintf(stderr, "Unknown option: %s\n", arg);
            print_usage(argv[0]);
            return 1;
        }
    }

    if (!input_file) {
        fprintf(stderr, "Error: no input file specified\n");
        print_usage(argv[0]);
        return 1;
    }

    /* Validate */
    if (cfg.width <= 0 || cfg.height <= 0) {
        fprintf(stderr, "Error: invalid image dimensions %dx%d\n", cfg.width, cfg.height);
        return 1;
    }
    if (cfg.aa_samples < 1) cfg.aa_samples = 1;
    if (cfg.aa_samples > RENDER_MAX_AA) cfg.aa_samples = RENDER_MAX_AA;

    /* Set log level */
    alea_log_set_level(cfg.log_level);

    /* ====================================================================
     * LOAD MODEL
     * ==================================================================== */

    fprintf(stderr, "Loading %s...\n", input_file);
    double t0 = get_time_ms();

    alea_system_t* sys = NULL;
    mcnp_model_t* model = NULL;
    openmc_model_t* omc_model = NULL;
    if (ends_with(input_file, ".xml")) {
        omc_model = openmc_load(input_file);
        if (omc_model) sys = omc_model->sys;
    } else {
        model = mcnp_load(input_file);
        if (model) sys = model->sys;
    }

    if (!sys) {
        fprintf(stderr, "Error loading %s: %s\n", input_file, alea_error());
        return 1;
    }

    double t1 = get_time_ms();
    fprintf(stderr, "Loaded: %zu cells, %zu surfaces (%.1f ms)\n",
            alea_cell_count(sys), alea_surface_count(sys), t1 - t0);

    /* Configure dedup */
    if (cfg.dedup) {
        alea_config_t scfg = alea_get_config(sys);
        scfg.dedup = true;
        alea_set_config(sys, &scfg);
    }

    /* Build mode-aware query caches */
    alea_build_universe_index(sys);
    alea_prepare_query_acceleration(sys);

    /* Print summary */
    alea_print_summary(sys);

    /* ====================================================================
     * SET UP CAMERA
     * ==================================================================== */

    render_camera_t cam;
    if (render_camera_setup(&cam, &cfg, sys) != 0) {
        fprintf(stderr, "Error: camera setup failed\n");
        if (model) mcnp_model_destroy(model); else if (omc_model) openmc_model_destroy(omc_model); else alea_destroy(sys);
        return 1;
    }

    fprintf(stderr, "Camera: eye=(%.1f, %.1f, %.1f) target=(%.1f, %.1f, %.1f) fov=%.0f%s\n",
            cam.eye[0], cam.eye[1], cam.eye[2],
            cam.target[0], cam.target[1], cam.target[2],
            cam.fov_degrees, cam.auto_fit ? " (auto)" : "");

    if (cfg.num_clips > 0) {
        fprintf(stderr, "Clip planes: %d\n", cfg.num_clips);
    }

    /* ====================================================================
     * PREVIEW PASS (optional)
     * ==================================================================== */

    if (cfg.preview) {
        fprintf(stderr, "\n--- Preview pass (1/4 resolution) ---\n");
        render_config_t preview_cfg = cfg;
        preview_cfg.width /= 4;
        preview_cfg.height /= 4;
        preview_cfg.aa_samples = 1;
        preview_cfg.shadows = 0;

        render_camera_t preview_cam;
        render_camera_setup(&preview_cam, &preview_cfg, sys);

        render_framebuffer_t* preview_fb =
            render_framebuffer_create(preview_cfg.width, preview_cfg.height, 0);
        if (preview_fb) {
            double tp0 = get_time_ms();
            render_scene(sys, &preview_cfg, &preview_cam, preview_fb);
            double tp1 = get_time_ms();

            uint8_t* preview_pixels = malloc((size_t)preview_cfg.width * preview_cfg.height * 3);
            if (preview_pixels) {
                render_tonemap(preview_fb, preview_pixels);

                char preview_path[600];
                const char* dot = strrchr(output_file, '.');
                if (dot) {
                    int base_len = (int)(dot - output_file);
                    snprintf(preview_path, sizeof(preview_path),
                             "%.*s_preview%s", base_len, output_file, dot);
                } else {
                    snprintf(preview_path, sizeof(preview_path),
                             "%s_preview.png", output_file);
                }

                render_write_image(preview_path, preview_pixels,
                                   preview_cfg.width, preview_cfg.height);
                fprintf(stderr, "Preview: %s (%dx%d, %.1f ms)\n",
                        preview_path, preview_cfg.width, preview_cfg.height, tp1 - tp0);
                free(preview_pixels);
            }
            render_framebuffer_free(preview_fb);
        }
    }

    /* ====================================================================
     * FULL RENDER
     * ==================================================================== */

    fprintf(stderr, "\n--- Full render ---\n");

    render_framebuffer_t* fb = render_framebuffer_create(cfg.width, cfg.height,
                                                          cfg.aux_output);
    if (!fb) {
        fprintf(stderr, "Error: failed to allocate framebuffer (%dx%d)\n",
                cfg.width, cfg.height);
        if (model) mcnp_model_destroy(model); else if (omc_model) openmc_model_destroy(omc_model); else alea_destroy(sys);
        return 1;
    }

    double tr0 = get_time_ms();
    int rc = render_scene(sys, &cfg, &cam, fb);
    double tr1 = get_time_ms();

    if (rc != 0) {
        fprintf(stderr, "Error: rendering failed\n");
        render_framebuffer_free(fb);
        if (model) mcnp_model_destroy(model); else if (omc_model) openmc_model_destroy(omc_model); else alea_destroy(sys);
        return 1;
    }

    fprintf(stderr, "Render time: %.1f ms (%.1f s)\n", tr1 - tr0, (tr1 - tr0) / 1000.0);

    /* Post-processing */
    if (cfg.edges) {
        fprintf(stderr, "Post-processing: edge darkening...\n");
        render_edge_darken(fb);
    }

    if (cfg.render_mode == RENDER_MODE_DEPTH) {
        /* Depth mode renders with solid shading then overwrites with depth grayscale */
        /* The depth data was collected during rendering; now convert to visual */
        if (fb->depth) {
            int w = fb->width, h = fb->height;
            float d_min = 1e30f, d_max = 0;
            for (int i = 0; i < w * h; i++) {
                if (fb->cell_id[i] >= 0) {
                    if (fb->depth[i] < d_min) d_min = fb->depth[i];
                    if (fb->depth[i] > d_max) d_max = fb->depth[i];
                }
            }
            float range = (d_max > d_min) ? (d_max - d_min) : 1.0f;
            for (int i = 0; i < w * h; i++) {
                if (fb->cell_id[i] >= 0) {
                    float t = 1.0f - (fb->depth[i] - d_min) / range;
                    fb->color[i * 3 + 0] = t;
                    fb->color[i * 3 + 1] = t;
                    fb->color[i * 3 + 2] = t;
                }
            }
        }
    } else if (cfg.render_mode == RENDER_MODE_CELLID) {
        /* Overwrite colors with palette by cell ID */
        int w = fb->width, h = fb->height;
        for (int i = 0; i < w * h; i++) {
            if (fb->cell_id[i] >= 0) {
                float r, g, b;
                render_get_color(fb->cell_id[i], RENDER_COLOR_CELL, &cfg, &r, &g, &b);
                fb->color[i * 3 + 0] = r;
                fb->color[i * 3 + 1] = g;
                fb->color[i * 3 + 2] = b;
            }
        }
    } else if (cfg.render_mode == RENDER_MODE_MATID) {
        int w = fb->width, h = fb->height;
        for (int i = 0; i < w * h; i++) {
            if (fb->cell_id[i] >= 0) {
                int mid = fb->material_id ? fb->material_id[i] : 0;
                float r, g, b;
                render_get_color(mid, RENDER_COLOR_MATERIAL, &cfg, &r, &g, &b);
                fb->color[i * 3 + 0] = r;
                fb->color[i * 3 + 1] = g;
                fb->color[i * 3 + 2] = b;
            }
        }
    }

    /* ====================================================================
     * WRITE OUTPUT
     * ==================================================================== */

    uint8_t* pixels = malloc((size_t)cfg.width * cfg.height * 3);
    if (!pixels) {
        fprintf(stderr, "Error: failed to allocate pixel buffer\n");
        render_framebuffer_free(fb);
        if (model) mcnp_model_destroy(model); else if (omc_model) openmc_model_destroy(omc_model); else alea_destroy(sys);
        return 1;
    }

    render_tonemap(fb, pixels);

    double tw0 = get_time_ms();
    rc = render_write_image(output_file, pixels, cfg.width, cfg.height);
    double tw1 = get_time_ms();

    if (rc != 0) {
        fprintf(stderr, "Error: failed to write %s\n", output_file);
    } else {
        fprintf(stderr, "Output: %s (%dx%d, %.0f ms)\n",
                output_file, cfg.width, cfg.height, tw1 - tw0);
    }

    /* Auxiliary outputs */
    if (cfg.aux_output) {
        fprintf(stderr, "Writing auxiliary maps...\n");
        render_write_aux(output_file, fb);
    }

    /* Print stats */
    double total = tw1 - t0;
    fprintf(stderr, "\nTotal time: %.1f s\n", total / 1000.0);

    /* Cleanup */
    free(pixels);
    render_framebuffer_free(fb);
    render_config_free(&cfg);
    if (model) mcnp_model_destroy(model); else if (omc_model) openmc_model_destroy(omc_model); else alea_destroy(sys);

    return 0;
}
