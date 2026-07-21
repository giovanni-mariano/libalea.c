// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/**
 * @file plotter.c
 * @brief CSG geometry slice plotter
 *
 * Renders 2D slices of CSG geometry using:
 * - alea_find_cells_grid_*() for cell/material queries
 * - alea_get_slice_curves_*() for surface boundaries
 * - alea_find_label_positions() for label placement
 *
 * Features:
 * - PNG and BMP output (auto-detected from extension)
 * - Axis-aligned slices (X, Y, Z) and arbitrary planes
 * - Color by cell ID or material ID
 * - Optional tick labels showing coordinates
 * - Batch mode for multiple plots
 *
 * Build:
 *   gcc -O2 -fopenmp -o plotter plotter.c -I../include \
 *       -L../bin -lalea_slice -lalea -lm
 *
 * Usage (single plot):
 *   ./plotter input.i <axis> <value> <u_min> <u_max> <v_min> <v_max> <WxH> [output.png] [options]
 *
 *   axis = Z  -> XY plane at z=value
 *   axis = Y  -> XZ plane at y=value
 *   axis = X  -> YZ plane at x=value
 *
 *   WxH = resolution (e.g., 800 for square, 800x600 for non-square)
 *
 * Usage (batch mode):
 *   ./plotter input.i --batch=plots.txt
 *
 * Options:
 *   --labels=cells          Show cell ID labels at region centroids
 *   --labels=materials      Show material ID labels at region centroids
 *   --labels=surfaces       Show surface ID labels along boundaries
 *   --labels=all            Show all labels
 *   --color=cells           Color by cell ID (default)
 *   --color=materials       Color by material ID
 *   --contours=cells        Draw contours at cell boundaries (default)
 *   --contours=materials    Draw contours at material boundaries
 *   --ticks                 Show axis tick labels
 *   --errors                Show analytical error lines (overlaps/gaps)
 *
 * Batch file format (one plot per line):
 *   Z value u_min u_max v_min v_max WxH output.png [options]
 *   Y value u_min u_max v_min v_max WxH output.png [options]
 *   X value u_min u_max v_min v_max WxH output.png [options]
 *   PLANE  ox oy oz nx ny nz ux uy uz u_min u_max v_min v_max WxH output [options]
 *   PLANE3 x1 y1 z1 x2 y2 z2 x3 y3 z3 u_min u_max v_min v_max WxH output [options]
 *
 *   Batch options: labels=... color=cells|materials ticks
 *   Lines starting with # are comments
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
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
#include "alea_slice.h"  /* Includes alea_slice_curve_set_debug() */
#include "alea_geo_validator.h"

typedef struct {
    alea_system_t* sys;
    mcnp_model_t* mcnp;
    openmc_model_t* openmc;
} loaded_model_t;

static void count_plot_error_state(const int* cell_ids,
                                   const uint8_t* coverage,
                                   const uint8_t* errors,
                                   int n,
                                   int* cov_none,
                                   int* cov_one,
                                   int* cov_multi,
                                   int* err_overlap,
                                   int* err_undefined) {
    int none = 0, one = 0, multi = 0, overlap = 0, undefined = 0;
    for (int i = 0; i < n; i++) {
        if (coverage) {
            if (coverage[i] == ALEA_COVERAGE_MULTI) multi++;
            else if (coverage[i] == ALEA_COVERAGE_ONE) one++;
            else none++;
        } else if (cell_ids && cell_ids[i] < 0) {
            none++;
        } else {
            one++;
        }

        if (errors) {
            if (errors[i] == ALEA_GRID_OVERLAP) overlap++;
            else if (errors[i] == ALEA_GRID_UNDEFINED) undefined++;
        }
    }
    if (cov_none) *cov_none = none;
    if (cov_one) *cov_one = one;
    if (cov_multi) *cov_multi = multi;
    if (err_overlap) *err_overlap = overlap;
    if (err_undefined) *err_undefined = undefined;
}

static void print_point_coverage_stats(const char* label) {
    alea_point_coverage_stats_t ps = alea_point_coverage_stats_get();
    if (ps.queries == 0) return;

    double avg_candidates = ps.spatial_queries > 0
        ? (double)ps.candidate_total / (double)ps.spatial_queries : 0.0;
    double tests_per_query = ps.queries > 0
        ? (double)ps.contains_tests / (double)ps.queries : 0.0;
    printf("    [PLOT ERROR STATS] %s point coverage: queries=%zu spatial=%zu recursive_fallback=%zu lattice_fallback=%zu truncated_fallback=%zu early_multi=%zu candidates avg=%.1f max=%zu contains_per_query=%.1f query_errors=%zu\n",
           label, ps.queries, ps.spatial_queries, ps.recursive_fallbacks,
           ps.lattice_fallbacks, ps.truncated_fallbacks,
           ps.spatial_multi_early_exit, avg_candidates,
           ps.candidate_max, tests_per_query, ps.query_errors);
}

static void apply_boundary_ambiguity_filter(
    alea_system_t* sys,
    const alea_slice_view_t* view,
    int width,
    int height,
    const int* cell_ids,
    int* secondary_ids,
    uint8_t* coverage,
    uint8_t* errors,
    int plot_error_stats,
    const char* label) {
    alea_boundary_filter_stats_t bfs;
    int suppressed = alea_filter_grid_boundary_ambiguities(
        sys, view, width, height, -1, cell_ids, secondary_ids,
        coverage, errors, &bfs);
    if (plot_error_stats && suppressed >= 0 &&
        (bfs.checked > 0 || bfs.inconclusive > 0)) {
        printf("    [PLOT ERROR STATS] %s boundary filter: checked=%zu suppressed=%zu retained=%zu inconclusive=%zu\n",
               label, bfs.checked, bfs.suppressed,
               bfs.retained, bfs.inconclusive);
    }
}

static void read_tile_coverage_size(int* tile_w, int* tile_h) {
    int w = 16;
    int h = 16;
    const char* size_env = getenv("ALEA_PLOT_TILE_SIZE");
    if (size_env && size_env[0]) {
        char* end = NULL;
        long parsed_w = strtol(size_env, &end, 10);
        if (parsed_w > 0) {
            w = (int)parsed_w;
            if (*end == 'x' || *end == 'X' || *end == ',') {
                char* end_h = NULL;
                long parsed_h = strtol(end + 1, &end_h, 10);
                if (parsed_h > 0) h = (int)parsed_h;
            } else {
                h = w;
            }
        }
    }

    const char* w_env = getenv("ALEA_PLOT_TILE_W");
    const char* h_env = getenv("ALEA_PLOT_TILE_H");
    if (w_env && atoi(w_env) > 0) w = atoi(w_env);
    if (h_env && atoi(h_env) > 0) h = atoi(h_env);

    if (w < 1) w = 1;
    if (h < 1) h = 1;
    if (w > 256) w = 256;
    if (h > 256) h = 256;
    *tile_w = w;
    *tile_h = h;
}

static uint8_t exact_coverage_at_point(alea_system_t* sys,
                                       double x, double y, double z,
                                       int universe_depth) {
    alea_cell_hit_t hits[64];
    int num_hits = alea_find_all_cells(sys, x, y, z, hits, 64);
    if (num_hits <= 0) return ALEA_COVERAGE_NONE;

    int target_depth = universe_depth;
    if (target_depth < 0)
        target_depth = hits[num_hits - 1].depth;

    int count = 0;
    for (int h = 0; h < num_hits; h++) {
        if (hits[h].depth != target_depth) continue;
        count++;
        if (count > 1) return ALEA_COVERAGE_MULTI;
    }
    return count == 1 ? ALEA_COVERAGE_ONE : ALEA_COVERAGE_NONE;
}

static int verify_plot_coverage_sample(alea_system_t* sys,
                                       const alea_slice_view_t* view,
                                       int nu, int nv,
                                       int universe_depth,
                                       const uint8_t* coverage,
                                       int interval,
                                       int* out_checked,
                                       int* out_mismatches) {
    if (!sys || !view || !coverage || nu <= 0 || nv <= 0 || interval <= 0)
        return -1;

    const alea_slice_plane_t* plane = &view->plane;
    double du = (view->u_max - view->u_min) / nu;
    double dv = (view->v_max - view->v_min) / nv;
    int checked = 0;
    int mismatches = 0;

    for (int j = 0; j < nv; j++) {
        double v = view->v_min + (j + 0.5) * dv;
        for (int i = 0; i < nu; i++) {
            int idx = j * nu + i;
            if (idx % interval != 0) continue;

            double u = view->u_min + (i + 0.5) * du;
            double x = plane->origin[0] + u * plane->u_axis[0] + v * plane->v_axis[0];
            double y = plane->origin[1] + u * plane->u_axis[1] + v * plane->v_axis[1];
            double z = plane->origin[2] + u * plane->u_axis[2] + v * plane->v_axis[2];
            uint8_t exact = exact_coverage_at_point(sys, x, y, z, universe_depth);
            checked++;
            if (coverage[idx] != exact) mismatches++;
        }
    }

    if (out_checked) *out_checked = checked;
    if (out_mismatches) *out_mismatches = mismatches;
    return mismatches == 0 ? 0 : 1;
}

static int ends_with(const char* s, const char* suffix) {
    size_t slen = strlen(s), xlen = strlen(suffix);
    return slen >= xlen && strcmp(s + slen - xlen, suffix) == 0;
}

static int load_model(const char* path, loaded_model_t* m) {
    memset(m, 0, sizeof(*m));
    if (ends_with(path, ".xml")) {
        m->openmc = openmc_load(path);
        if (!m->openmc) return -1;
        m->sys = m->openmc->sys;
    } else {
        m->mcnp = mcnp_load(path);
        if (!m->mcnp) return -1;
        m->sys = m->mcnp->sys;
    }
    return 0;
}

static void destroy_model(loaded_model_t* m) {
    if (m->mcnp)   mcnp_model_destroy(m->mcnp);
    if (m->openmc)  openmc_model_destroy(m->openmc);
}

/* Get resident set size in MB (Linux only) */
static double get_rss_mb(void) {
#ifdef _WIN32
    return 0;
#else
    FILE* f = fopen("/proc/self/status", "r");
    if (!f) return 0;
    char line[256];
    double rss = 0;
    while (fgets(line, sizeof(line), f)) {
        if (strncmp(line, "VmRSS:", 6) == 0) {
            long kb = 0;
            sscanf(line + 6, "%ld", &kb);
            rss = kb / 1024.0;
            break;
        }
    }
    fclose(f);
    return rss;
#endif
}

/* ============================================================================
 * LABEL RENDERING - Simple 5x7 bitmap font
 * ============================================================================ */

/* Label display options */
typedef struct {
    int show_cells;
    int show_materials;
    int show_surfaces;
} label_options_t;

/* 5x7 bitmap font for digits 0-9 and some characters */
static const uint8_t font_5x7[][7] = {
    /* '0' */ {0x0E, 0x11, 0x13, 0x15, 0x19, 0x11, 0x0E},
    /* '1' */ {0x04, 0x0C, 0x04, 0x04, 0x04, 0x04, 0x0E},
    /* '2' */ {0x0E, 0x11, 0x01, 0x02, 0x04, 0x08, 0x1F},
    /* '3' */ {0x1F, 0x02, 0x04, 0x02, 0x01, 0x11, 0x0E},
    /* '4' */ {0x02, 0x06, 0x0A, 0x12, 0x1F, 0x02, 0x02},
    /* '5' */ {0x1F, 0x10, 0x1E, 0x01, 0x01, 0x11, 0x0E},
    /* '6' */ {0x06, 0x08, 0x10, 0x1E, 0x11, 0x11, 0x0E},
    /* '7' */ {0x1F, 0x01, 0x02, 0x04, 0x08, 0x08, 0x08},
    /* '8' */ {0x0E, 0x11, 0x11, 0x0E, 0x11, 0x11, 0x0E},
    /* '9' */ {0x0E, 0x11, 0x11, 0x0F, 0x01, 0x02, 0x0C},
    /* 'C' */ {0x0E, 0x11, 0x10, 0x10, 0x10, 0x11, 0x0E},
    /* 'M' */ {0x11, 0x1B, 0x15, 0x15, 0x11, 0x11, 0x11},
    /* 'S' */ {0x0E, 0x11, 0x10, 0x0E, 0x01, 0x11, 0x0E},
    /* '-' */ {0x00, 0x00, 0x00, 0x1F, 0x00, 0x00, 0x00},
    /* ' ' */ {0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00},
};

#define FONT_WIDTH 5
#define FONT_HEIGHT 7
#define CHAR_SPACING 1

/* Get font index for a character */
static int font_index(char c) {
    if (c >= '0' && c <= '9') return c - '0';
    if (c == 'C') return 10;
    if (c == 'M') return 11;
    if (c == 'S') return 12;
    if (c == '-') return 13;
    return 14; /* space */
}

/* Draw a single character at (x, y) with given color */
static void draw_char(uint8_t* pixels, int width, int height,
                      int x, int y, char c,
                      uint8_t r, uint8_t g, uint8_t b) {
    int idx = font_index(c);
    for (int row = 0; row < FONT_HEIGHT; row++) {
        int py = y + row;
        if (py < 0 || py >= height) continue;
        uint8_t bits = font_5x7[idx][row];
        for (int col = 0; col < FONT_WIDTH; col++) {
            int px = x + col;
            if (px < 0 || px >= width) continue;
            if (bits & (0x10 >> col)) {
                int pidx = (py * width + px) * 3;
                pixels[pidx + 0] = r;
                pixels[pidx + 1] = g;
                pixels[pidx + 2] = b;
            }
        }
    }
}

/* Draw a string at (x, y) - centered horizontally */
static void draw_string_centered(uint8_t* pixels, int width, int height,
                                  int cx, int cy, const char* str,
                                  uint8_t r, uint8_t g, uint8_t b) {
    int len = (int)strlen(str);
    int text_width = len * (FONT_WIDTH + CHAR_SPACING) - CHAR_SPACING;
    int x = cx - text_width / 2;
    int y = cy - FONT_HEIGHT / 2;

    /* Draw outline (black background for readability) */
    for (int dy = -1; dy <= 1; dy++) {
        for (int dx = -1; dx <= 1; dx++) {
            if (dx == 0 && dy == 0) continue;
            int ox = x + dx;
            int oy = y + dy;
            for (int i = 0; str[i]; i++) {
                draw_char(pixels, width, height, ox, oy, str[i], 0, 0, 0);
                ox += FONT_WIDTH + CHAR_SPACING;
            }
        }
    }

    /* Draw text */
    for (int i = 0; str[i]; i++) {
        draw_char(pixels, width, height, x, y, str[i], r, g, b);
        x += FONT_WIDTH + CHAR_SPACING;
    }
}

/* ============================================================================
 * LABEL DRAWING (uses library functions for position computation)
 * ============================================================================ */

/* Draw cell labels using library's interior point algorithm */
static void draw_cell_labels(uint8_t* pixels, int width, int height,
                              const int* cell_ids) {
    alea_label_position_t* labels = NULL;
    int count = 0;

    if (alea_find_label_positions(cell_ids, width, height, 100, &labels, &count) != 0) {
        return;
    }

    for (int i = 0; i < count; i++) {
        char label[32];
        snprintf(label, sizeof(label), "C%d", labels[i].id);
        int py = height - 1 - labels[i].py;  /* grid order -> image order */
        draw_string_centered(pixels, width, height, labels[i].px, py,
                            label, 255, 255, 255);
    }

    free(labels);
}

/* Draw material labels using library's interior point algorithm */
static void draw_material_labels(uint8_t* pixels, int width, int height,
                                  const int* cell_ids, const int* material_ids) {
    /* Get cell label positions */
    alea_label_position_t* cell_labels = NULL;
    int cell_count = 0;

    if (alea_find_label_positions(cell_ids, width, height, 100, &cell_labels, &cell_count) != 0) {
        return;
    }

    /* Build cell-to-material map */
    int cell_to_mat[10000];
    memset(cell_to_mat, 0, sizeof(cell_to_mat));

    for (int i = 0; i < width * height; i++) {
        int cell = cell_ids[i];
        if (cell >= 0 && cell < 10000) {
            cell_to_mat[cell] = material_ids[i];
        }
    }

    /* Draw material labels offset from cell positions */
    for (int i = 0; i < cell_count; i++) {
        int mat = (cell_labels[i].id < 10000) ? cell_to_mat[cell_labels[i].id] : 0;
        char label[32];
        snprintf(label, sizeof(label), "M%d", mat);
        int py = height - 1 - cell_labels[i].py;  /* grid order -> image order */
        /* Offset below cell label to avoid overlap */
        draw_string_centered(pixels, width, height,
                            cell_labels[i].px, py + FONT_HEIGHT + 2,
                            label, 255, 200, 100);
    }

    free(cell_labels);
}

/* ============================================================================
 * COLOR PALETTE AND CONTOURS
 * ============================================================================ */

/* Default color palette (same as slice_analytical.c) */
static const uint8_t palette[][3] = {
    {0x1f, 0x77, 0xb4}, {0xff, 0x7f, 0x0e}, {0x2c, 0xa0, 0x2c}, {0xd6, 0x27, 0x28},
    {0x94, 0x67, 0xbd}, {0x8c, 0x56, 0x4b}, {0xe3, 0x77, 0xc2}, {0x7f, 0x7f, 0x7f},
    {0xbc, 0xbd, 0x22}, {0x17, 0xbe, 0xcf}, {0xae, 0xc7, 0xe8}, {0xff, 0xbb, 0x78},
    {0x98, 0xdf, 0x8a}, {0xff, 0x98, 0x96}, {0xc5, 0xb0, 0xd5}, {0xc4, 0x9c, 0x94}
};
#define PALETTE_SIZE 16

/* Error codes (must match alea_grid_error_t) */
#define GRID_ERR_NONE      0
#define GRID_ERR_OVERLAP   1
#define GRID_ERR_UNDEFINED 2

/* Write BMP file */
static int write_bmp(const char* filename, const uint8_t* pixels, int width, int height) {
    FILE* f = fopen(filename, "wb");
    if (!f) return -1;

    int row_size = ((width * 3 + 3) / 4) * 4;
    int data_size = row_size * height;
    int file_size = 54 + data_size;

    uint8_t header[54] = {0};
    header[0] = 'B'; header[1] = 'M';
    header[2] = file_size & 0xFF;
    header[3] = (file_size >> 8) & 0xFF;
    header[4] = (file_size >> 16) & 0xFF;
    header[5] = (file_size >> 24) & 0xFF;
    header[10] = 54;
    header[14] = 40;
    header[18] = width & 0xFF;
    header[19] = (width >> 8) & 0xFF;
    header[20] = (width >> 16) & 0xFF;
    header[21] = (width >> 24) & 0xFF;
    header[22] = height & 0xFF;
    header[23] = (height >> 8) & 0xFF;
    header[24] = (height >> 16) & 0xFF;
    header[25] = (height >> 24) & 0xFF;
    header[26] = 1;
    header[28] = 24;

    fwrite(header, 1, 54, f);

    uint8_t* row = malloc(row_size);
    for (int y = height - 1; y >= 0; y--) {
        memset(row, 0, row_size);
        for (int x = 0; x < width; x++) {
            int src_idx = (y * width + x) * 3;
            int dst_idx = x * 3;
            row[dst_idx + 0] = pixels[src_idx + 2];  /* BGR */
            row[dst_idx + 1] = pixels[src_idx + 1];
            row[dst_idx + 2] = pixels[src_idx + 0];
        }
        fwrite(row, 1, row_size, f);
    }

    free(row);
    fclose(f);
    return 0;
}

/* ============================================================================
 * PNG WRITING - Minimal implementation (no dependencies)
 * ============================================================================ */

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

    if (len > 0 && data) {
        fwrite(data, 1, len, f);
    }

    /* CRC covers type + data */
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

/* Write PNG file (RGB pixels, no compression) */
static int write_png(const char* filename, const uint8_t* pixels, int width, int height) {
    FILE* f = fopen(filename, "wb");
    if (!f) return -1;

    /* PNG signature */
    static const uint8_t png_sig[8] = {137, 80, 78, 71, 13, 10, 26, 10};
    fwrite(png_sig, 1, 8, f);

    /* IHDR chunk */
    uint8_t ihdr[13];
    png_write_be32(ihdr, width);
    png_write_be32(ihdr + 4, height);
    ihdr[8] = 8;   /* bit depth */
    ihdr[9] = 2;   /* color type: RGB */
    ihdr[10] = 0;  /* compression */
    ihdr[11] = 0;  /* filter */
    ihdr[12] = 0;  /* interlace */
    png_write_chunk(f, "IHDR", ihdr, 13);

    /* Prepare raw image data with filter bytes */
    size_t raw_row = 1 + width * 3;  /* filter byte + RGB */
    size_t raw_size = raw_row * height;
    uint8_t* raw = malloc(raw_size);
    if (!raw) {
        fclose(f);
        return -1;
    }

    for (int y = 0; y < height; y++) {
        raw[y * raw_row] = 0;  /* filter: none */
        memcpy(raw + y * raw_row + 1, pixels + y * width * 3, width * 3);
    }

    /* Build zlib stream with uncompressed DEFLATE blocks */
    size_t max_block = 65535;
    size_t num_blocks = (raw_size + max_block - 1) / max_block;
    size_t zlib_size = 2 + raw_size + num_blocks * 5 + 4;

    uint8_t* zlib = malloc(zlib_size);
    if (!zlib) {
        free(raw);
        fclose(f);
        return -1;
    }

    size_t zp = 0;
    zlib[zp++] = 0x78;  /* CMF */
    zlib[zp++] = 0x01;  /* FLG */

    size_t remaining = raw_size;
    size_t src_pos = 0;
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

    /* IEND chunk */
    png_write_chunk(f, "IEND", NULL, 0);

    fclose(f);
    return 0;
}

/* Write image file (auto-detect format from extension) */
static int write_image(const char* filename, const uint8_t* pixels, int width, int height) {
    const char* ext = strrchr(filename, '.');
    if (ext && (strcasecmp(ext, ".png") == 0)) {
        return write_png(filename, pixels, width, height);
    }
    return write_bmp(filename, pixels, width, height);
}

/**
 * @brief Draw contours with error-aware coloring
 *
 * Normal boundaries: solid black
 * Overlap boundaries: dashed red
 * Undefined boundaries: dashed orange
 */
static void draw_contours_ex(uint8_t* pixels, const int* boundary_ids,
                              const uint8_t* errors, int width, int height) {
    /* Dash pattern: 3 pixels on, 3 pixels off */
    const int dash_period = 6;
    const int dash_on = 3;

    /* cell_ids is in grid order (row 0 = v_min), pixels in image order
     * (row 0 = v_max). Map pixel row y to grid row (height-1-y). */
    for (int y = 0; y < height; y++) {
        int gy = height - 1 - y;
        for (int x = 0; x < width; x++) {
            int gidx = gy * width + x;
            int cell = boundary_ids[gidx];
            int is_boundary = 0;
            int error_type = GRID_ERR_NONE;

            /* Check right neighbor (same grid row) */
            if (x + 1 < width && boundary_ids[gidx + 1] != cell) {
                is_boundary = 1;
                if (errors) {
                    uint8_t err_here = errors[gidx];
                    uint8_t err_right = errors[gidx + 1];
                    if (err_here == GRID_ERR_OVERLAP || err_right == GRID_ERR_OVERLAP) {
                        error_type = GRID_ERR_OVERLAP;
                    } else if (err_here == GRID_ERR_UNDEFINED || err_right == GRID_ERR_UNDEFINED) {
                        error_type = GRID_ERR_UNDEFINED;
                    }
                }
            }

            /* Check bottom pixel neighbor = grid row gy-1 */
            if (gy > 0 && boundary_ids[(gy - 1) * width + x] != cell) {
                is_boundary = 1;
                if (errors) {
                    uint8_t err_here = errors[gidx];
                    uint8_t err_below = errors[(gy - 1) * width + x];
                    if (err_here == GRID_ERR_OVERLAP || err_below == GRID_ERR_OVERLAP) {
                        error_type = GRID_ERR_OVERLAP;
                    } else if (error_type != GRID_ERR_OVERLAP &&
                               (err_here == GRID_ERR_UNDEFINED || err_below == GRID_ERR_UNDEFINED)) {
                        error_type = GRID_ERR_UNDEFINED;
                    }
                }
            }

            if (is_boundary) {
                /* For error boundaries, use dashed pattern */
                if (error_type != GRID_ERR_NONE) {
                    int dash_pos = (x + y) % dash_period;
                    if (dash_pos >= dash_on) continue;  /* Skip "off" pixels */
                }

                int pidx = (y * width + x) * 3;
                if (error_type == GRID_ERR_OVERLAP) {
                    /* Red for overlap */
                    pixels[pidx + 0] = 255;
                    pixels[pidx + 1] = 0;
                    pixels[pidx + 2] = 0;
                } else if (error_type == GRID_ERR_UNDEFINED) {
                    /* Orange for undefined */
                    pixels[pidx + 0] = 255;
                    pixels[pidx + 1] = 165;
                    pixels[pidx + 2] = 0;
                } else {
                    /* Black for normal boundary */
                    pixels[pidx + 0] = 0;
                    pixels[pidx + 1] = 0;
                    pixels[pidx + 2] = 0;
                }
            }
        }
    }
}


/* Parse label options from string */
static label_options_t parse_label_string(const char* str) {
    label_options_t opts = {0, 0, 0};
    if (!str) return opts;

    char buf[256];
    strncpy(buf, str, sizeof(buf) - 1);
    buf[sizeof(buf) - 1] = '\0';

    char* token = strtok(buf, ",");
    while (token) {
        while (*token == ' ') token++;  /* skip leading spaces */
        if (strcmp(token, "cells") == 0) opts.show_cells = 1;
        else if (strcmp(token, "materials") == 0) opts.show_materials = 1;
        else if (strcmp(token, "surfaces") == 0) opts.show_surfaces = 1;
        else if (strcmp(token, "all") == 0) {
            opts.show_cells = 1;
            opts.show_materials = 1;
            opts.show_surfaces = 1;
        }
        token = strtok(NULL, ",");
    }
    return opts;
}

/* ============================================================================
 * PLOT PARAMETERS
 * ============================================================================ */

typedef enum {
    SLICE_Z = 0,    /* XY plane at z=constant */
    SLICE_Y = 1,    /* XZ plane at y=constant */
    SLICE_X = 2,    /* YZ plane at x=constant */
    SLICE_PLANE = 3 /* Arbitrary plane */
} slice_type_t;

typedef enum {
    COLOR_BY_CELL = 0,
    COLOR_BY_MATERIAL = 1
} color_mode_t;

typedef enum {
    CONTOUR_BY_CELL = 0,
    CONTOUR_BY_MATERIAL = 1
} contour_mode_t;

typedef struct {
    slice_type_t type;
    double value;           /* Axis value for Z/Y/X slices */
    double u_min, u_max;    /* Horizontal axis range */
    double v_min, v_max;    /* Vertical axis range */
    int width, height;      /* Image dimensions */
    char output[512];
    label_options_t labels;
    color_mode_t color_mode;
    contour_mode_t contour_mode;
    int show_ticks;         /* Draw axis tick labels */
    int show_errors;        /* Draw analytical error lines (overlaps/gaps) */

    /* For arbitrary planes */
    double origin[3];
    double normal[3];
    double up[3];
} plot_params_t;

/* Build an alea_slice_view_t from plot parameters */
static void build_view_from_params(const plot_params_t* p, alea_slice_view_t* view) {
    switch (p->type) {
        case SLICE_Z:
            alea_slice_view_axis(view, 2, p->value,
                p->u_min, p->u_max, p->v_min, p->v_max);
            break;
        case SLICE_Y:
            alea_slice_view_axis(view, 1, p->value,
                p->u_min, p->u_max, p->v_min, p->v_max);
            break;
        case SLICE_X:
            alea_slice_view_axis(view, 0, p->value,
                p->u_min, p->u_max, p->v_min, p->v_max);
            break;
        case SLICE_PLANE:
            alea_slice_view_init(view,
                p->origin[0], p->origin[1], p->origin[2],
                p->normal[0], p->normal[1], p->normal[2],
                p->up[0], p->up[1], p->up[2],
                p->u_min, p->u_max, p->v_min, p->v_max);
            break;
    }
}

/* ============================================================================
 * TICK LABEL RENDERING
 * ============================================================================ */

/* Margin for tick labels (pixels) */
#define TICK_MARGIN 40
#define TICK_LENGTH 5

/* Format a number nicely for tick labels */
static void format_tick_value(char* buf, size_t size, double value) {
    if (fabs(value) < 1e-10) {
        snprintf(buf, size, "0");
    } else if (fabs(value) >= 1000 || fabs(value) < 0.01) {
        snprintf(buf, size, "%.2g", value);
    } else if (fabs(value - (int)value) < 1e-10) {
        snprintf(buf, size, "%d", (int)value);
    } else {
        snprintf(buf, size, "%.2f", value);
    }
}

/* Draw tick marks and labels on the plot */
static void draw_ticks(uint8_t* pixels, int width, int height,
                       double u_min, double u_max, double v_min, double v_max,
                       const char* u_label, const char* v_label) {
    /* Calculate nice tick intervals */
    double u_range = u_max - u_min;
    double v_range = v_max - v_min;

    /* Target ~5-8 ticks per axis */
    double u_step = pow(10, floor(log10(u_range / 5)));
    double v_step = pow(10, floor(log10(v_range / 5)));

    /* Adjust to 1, 2, or 5 multiples */
    if (u_range / u_step > 8) u_step *= 2;
    if (u_range / u_step > 8) u_step *= 2.5;
    if (v_range / v_step > 8) v_step *= 2;
    if (v_range / v_step > 8) v_step *= 2.5;

    /* Draw bottom edge ticks (U axis) */
    double u_start = ceil(u_min / u_step) * u_step;
    for (double u = u_start; u <= u_max; u += u_step) {
        int px = (int)((u - u_min) / u_range * (width - 1));
        if (px < 0 || px >= width) continue;

        /* Draw tick mark */
        for (int y = height - TICK_LENGTH; y < height; y++) {
            int idx = (y * width + px) * 3;
            pixels[idx + 0] = 0;
            pixels[idx + 1] = 0;
            pixels[idx + 2] = 0;
        }

        /* Draw tick label */
        char label[32];
        format_tick_value(label, sizeof(label), u);
        int label_y = height - TICK_LENGTH - 2;
        draw_string_centered(pixels, width, height, px, label_y, label, 0, 0, 0);
    }

    /* Draw left edge ticks (V axis) */
    double v_start = ceil(v_min / v_step) * v_step;
    for (double v = v_start; v <= v_max; v += v_step) {
        int py = height - 1 - (int)((v - v_min) / v_range * (height - 1));
        if (py < 0 || py >= height) continue;

        /* Draw tick mark */
        for (int x = 0; x < TICK_LENGTH; x++) {
            int idx = (py * width + x) * 3;
            pixels[idx + 0] = 0;
            pixels[idx + 1] = 0;
            pixels[idx + 2] = 0;
        }

        /* Draw tick label */
        char label[32];
        format_tick_value(label, sizeof(label), v);
        int label_x = TICK_LENGTH + 2 + (int)strlen(label) * 3;
        draw_string_centered(pixels, width, height, label_x, py, label, 0, 0, 0);
    }

    /* Draw axis labels */
    if (u_label && *u_label) {
        int label_x = width / 2;
        int label_y = height - TICK_LENGTH - 12;
        draw_string_centered(pixels, width, height, label_x, label_y, u_label, 0, 0, 0);
    }
    if (v_label && *v_label) {
        int label_x = TICK_LENGTH + 12;
        int label_y = height / 2;
        draw_string_centered(pixels, width, height, label_x, label_y, v_label, 0, 0, 0);
    }
}

/* Get slice curves for surface labels */
static alea_slice_curves_t* get_curves_for_plot(alea_system_t* sys, const plot_params_t* p) {
    alea_slice_view_t view;
    build_view_from_params(p, &view);
    return alea_get_slice_curves(sys, &view);
}


/* Get axis labels for a slice type */
static void get_axis_labels(slice_type_t type, const char** u_label, const char** v_label) {
    switch (type) {
        case SLICE_Z: *u_label = "X"; *v_label = "Y"; break;
        case SLICE_Y: *u_label = "X"; *v_label = "Z"; break;
        case SLICE_X: *u_label = "Y"; *v_label = "Z"; break;
        default:      *u_label = "U"; *v_label = "V"; break;
    }
}

/* Draw geometry-validator boundary events (overlaps/gaps/etc.) onto the pixel
 * buffer as dots. Faithful two-colour scheme: red for overlaps, orange for every
 * other event type (undefined / non-adjacent / missing-neighbor / ambiguous). */
static void draw_validator_error_dots(uint8_t* pixels, int width, int height,
                                      const plot_params_t* p,
                                      const alea_slice_view_t* view,
                                      const alea_geom_validator_result_t* res) {
    double du = (p->u_max - p->u_min) / width;
    double dv = (p->v_max - p->v_min) / height;
    const double* origin = view->plane.origin;
    const double* u_axis = view->plane.u_axis;
    const double* v_axis = view->plane.v_axis;

    for (size_t i = 0; i < res->error_count; i++) {
        const alea_geom_error_t* e = &res->errors[i];

        /* Project the world crossing point onto plane (u,v) coordinates. */
        double dx = e->crossing_point[0] - origin[0];
        double dy = e->crossing_point[1] - origin[1];
        double dz = e->crossing_point[2] - origin[2];
        double u = dx * u_axis[0] + dy * u_axis[1] + dz * u_axis[2];
        double v = dx * v_axis[0] + dy * v_axis[1] + dz * v_axis[2];

        /* Choose color: red for overlap, orange for everything else. */
        uint8_t r, g, b;
        if (e->type == ALEA_GEOM_ERR_OVERLAP_AFTER_CROSSING) {
            r = 255; g = 0; b = 0;
        } else {
            r = 255; g = 165; b = 0;
        }

        int px = (int)((u - p->u_min) / du);
        int py_grid = (int)((v - p->v_min) / dv);
        int py = height - 1 - py_grid;  /* flip to image coords */

        /* Draw 3x3 dot */
        for (int ddy = -1; ddy <= 1; ddy++) {
            for (int ddx = -1; ddx <= 1; ddx++) {
                int x = px + ddx;
                int y = py + ddy;
                if (x < 0 || x >= width || y < 0 || y >= height)
                    continue;
                int pidx = (y * width + x) * 3;
                pixels[pidx + 0] = r;
                pixels[pidx + 1] = g;
                pixels[pidx + 2] = b;
            }
        }
    }
}

/* Render a single plot and save to file */
static int render_plot(alea_system_t* sys, const plot_params_t* p, int verbose) {
    double t0, t1;
    const char* axis_names[] = {"Z", "Y", "X", "Plane"};

    if (verbose) {
        if (p->type == SLICE_PLANE) {
            printf("  Plane origin=(%.4g,%.4g,%.4g) normal=(%.4g,%.4g,%.4g), %dx%d -> %s\n",
                   p->origin[0], p->origin[1], p->origin[2],
                   p->normal[0], p->normal[1], p->normal[2],
                   p->width, p->height, p->output);
        } else {
            printf("  %s=%.4g, bounds=[%.4g,%.4g]x[%.4g,%.4g], %dx%d -> %s\n",
                   axis_names[p->type], p->value,
                   p->u_min, p->u_max, p->v_min, p->v_max,
                   p->width, p->height, p->output);
        }
    }

    int num_pixels = p->width * p->height;
    int* cell_ids = malloc(num_pixels * sizeof(int));
    int* material_ids = malloc(num_pixels * sizeof(int));
    int* secondary_ids = NULL;
    uint32_t* path_ids = NULL;
    alea_slice_path_table_t path_table = {0};
    uint8_t* errors = malloc(num_pixels);
    uint8_t* coverage = malloc(num_pixels);
    uint8_t* pixels = malloc(num_pixels * 3);

    if (!cell_ids || !material_ids || !errors || !coverage || !pixels) {
        fprintf(stderr, "Error: Failed to allocate buffers\n");
        free(cell_ids); free(material_ids); free(secondary_ids); free(path_ids);
        alea_slice_path_table_free(&path_table);
        free(errors); free(coverage); free(pixels);
        return -1;
    }

    if (verbose) {
        printf("    [RSS before grid query: %.1f MB]\n", get_rss_mb());
    }

    t0 = get_time_ms();

    alea_slice_view_t view;
    build_view_from_params(p, &view);
    int selective_exact_coverage = getenv("ALEA_PLOT_EXACT_COVERAGE") != NULL;
    int full_exact_coverage = getenv("ALEA_PLOT_FULL_EXACT_COVERAGE") != NULL;
    int tile_exact_coverage = getenv("ALEA_PLOT_TILE_COVERAGE") != NULL;
    int path_exact_coverage = getenv("ALEA_PLOT_PATH_COVERAGE") != NULL;
    int plot_error_stats = getenv("ALEA_PLOT_ERROR_STATS") != NULL;
    const char* verify_env = getenv("ALEA_PLOT_ERROR_VERIFY");
    int verify_interval = (verify_env && atoi(verify_env) > 0)
        ? atoi(verify_env) : 0;
    if (plot_error_stats) {
        secondary_ids = malloc(num_pixels * sizeof(int));
        if (!secondary_ids) {
            fprintf(stderr, "Error: Failed to allocate secondary ID buffer\n");
            free(cell_ids); free(material_ids); free(errors);
            free(path_ids); alea_slice_path_table_free(&path_table);
            free(coverage); free(pixels);
            return -1;
        }
    }
    if (path_exact_coverage) {
        path_ids = malloc((size_t)num_pixels * sizeof(uint32_t));
        if (!path_ids) {
            fprintf(stderr, "Error: Failed to allocate path ID buffer\n");
            free(cell_ids); free(material_ids); free(secondary_ids);
            free(errors); free(coverage); free(pixels);
            return -1;
        }
    }
    unsigned coverage_flags = full_exact_coverage
        ? ALEA_GRID_COVERAGE_EXACT
        : ALEA_GRID_COVERAGE_FAST;
    if (secondary_ids) coverage_flags |= ALEA_GRID_SECONDARY_CELL_IDS;
    if (path_exact_coverage) coverage_flags |= ALEA_GRID_PATH_IDS;
    int rc = path_exact_coverage
        ? alea_find_cells_grid_coverage_paths(
              sys, &view, p->width, p->height, -1, coverage_flags,
              cell_ids, material_ids, secondary_ids, coverage, errors,
              path_ids, &path_table)
        : alea_find_cells_grid_coverage(
              sys, &view, p->width, p->height, -1, coverage_flags,
              cell_ids, material_ids, secondary_ids, coverage, errors);

    t1 = get_time_ms();

    if (rc != 0) {
        fprintf(stderr, "Error: Grid query failed\n");
        free(cell_ids); free(material_ids); free(secondary_ids); free(path_ids);
        alea_slice_path_table_free(&path_table);
        free(errors); free(coverage); free(pixels);
        return -1;
    }

    if (verbose) {
        printf("    Grid query: %.1f ms (%.2f Mpx/s)\n",
               t1 - t0, num_pixels / (t1 - t0) / 1000.0);
        if (full_exact_coverage) {
            printf("    Full exact coverage enabled via ALEA_PLOT_FULL_EXACT_COVERAGE\n");
        } else if (path_exact_coverage) {
            printf("    Path exact coverage enabled via ALEA_PLOT_PATH_COVERAGE (%zu paths)\n",
                   path_table.count);
        } else if (tile_exact_coverage) {
            printf("    Tile exact coverage enabled via ALEA_PLOT_TILE_COVERAGE\n");
        } else if (selective_exact_coverage) {
            printf("    Selective exact coverage enabled via ALEA_PLOT_EXACT_COVERAGE\n");
        }
        printf("    [RSS after grid query: %.1f MB]\n", get_rss_mb());
    }
    if (plot_error_stats && full_exact_coverage) {
        print_point_coverage_stats("full-grid exact");
    }

    if (path_exact_coverage && !full_exact_coverage) {
        int tile_w = 16, tile_h = 16;
        read_tile_coverage_size(&tile_w, &tile_h);
        t0 = get_time_ms();
        int refined = alea_refine_grid_coverage_paths_exact(
            sys, &view, p->width, p->height, -1, tile_w, tile_h,
            cell_ids, path_ids, &path_table, secondary_ids, coverage, errors);
        int used_path_refinement = refined >= 0;
        if (refined < 0) {
            if (verbose || plot_error_stats) {
                printf("    Path exact coverage unavailable; falling back to tile exact coverage\n");
            }
            refined = alea_refine_grid_coverage_tiles_exact(
                sys, &view, p->width, p->height, -1, tile_w, tile_h,
                secondary_ids, coverage, errors);
        }
        t1 = get_time_ms();
        if (verbose || plot_error_stats) {
            printf("    %s exact coverage: %d pixels refined, tile=%dx%d (%.1f ms)\n",
                   used_path_refinement ? "Path" : "Tile fallback",
                   refined, tile_w, tile_h, t1 - t0);
        }
        if (plot_error_stats) {
            alea_tile_coverage_stats_t ts = alea_tile_coverage_stats_get();
            double avg_raw = ts.tiles > 0
                ? (double)ts.candidate_total / (double)ts.tiles : 0.0;
            double avg_dedup = ts.tiles > 0
                ? (double)ts.dedup_candidate_total / (double)ts.tiles : 0.0;
            double tests_per_pixel = ts.pixels > 0
                ? (double)ts.candidate_pixel_tests / (double)ts.pixels : 0.0;
            double bbox_reject_pct = ts.bbox_pixel_tests > 0
                ? 100.0 * (double)ts.bbox_pixel_rejects /
                  (double)ts.bbox_pixel_tests : 0.0;
            printf("    [PLOT ERROR STATS] %s coverage: tiles=%zu fallback=%zu query_errors=%zu raw_candidates avg=%.1f max=%zu dedup avg=%.1f max=%zu tests_per_pixel=%.1f exact_fallback_pixels=%zu paths=%zu\n",
                   used_path_refinement ? "path" : "tile fallback",
                   ts.tiles, ts.fallback_tiles, ts.query_errors,
                   avg_raw, ts.candidate_max, avg_dedup,
                   ts.dedup_candidate_max, tests_per_pixel,
                   ts.exact_fallback_pixels, path_table.count);
            if (used_path_refinement) {
                printf("    [PLOT ERROR STATS] path groups=%zu max_pixels=%zu max_candidates=%zu bbox_tests=%zu bbox_rejects=%zu (%.1f%%) contains=%zu primary_skips=%zu early_multi_skips=%zu\n",
                       ts.path_groups, ts.path_group_pixels_max,
                       ts.path_group_candidates_max,
                       ts.bbox_pixel_tests, ts.bbox_pixel_rejects,
                       bbox_reject_pct, ts.candidate_pixel_tests,
                       ts.primary_cell_skips,
                       ts.early_multi_skips);
                if (ts.path_2d_verify_queries > 0 ||
                    ts.path_2d_missing_candidates > 0) {
                    printf("    [PLOT ERROR STATS] path 2d verify: queries=%zu missing_groups=%zu missing_candidates=%zu\n",
                           ts.path_2d_verify_queries,
                           ts.path_2d_missing_tiles,
                           ts.path_2d_missing_candidates);
                }
            }
            print_point_coverage_stats("path fallback");
        }
    } else if (tile_exact_coverage && !full_exact_coverage) {
        int tile_w = 16, tile_h = 16;
        read_tile_coverage_size(&tile_w, &tile_h);
        t0 = get_time_ms();
        int refined = alea_refine_grid_coverage_tiles_exact(
            sys, &view, p->width, p->height, -1, tile_w, tile_h,
            secondary_ids, coverage, errors);
        t1 = get_time_ms();
        if (verbose || plot_error_stats) {
            printf("    Tile exact coverage: %d pixels refined, tile=%dx%d (%.1f ms)\n",
                   refined, tile_w, tile_h, t1 - t0);
        }
        if (plot_error_stats) {
            alea_tile_coverage_stats_t ts = alea_tile_coverage_stats_get();
            double avg_raw = ts.tiles > 0
                ? (double)ts.candidate_total / (double)ts.tiles : 0.0;
            double avg_dedup = ts.tiles > 0
                ? (double)ts.dedup_candidate_total / (double)ts.tiles : 0.0;
            double tests_per_pixel = ts.pixels > 0
                ? (double)ts.candidate_pixel_tests / (double)ts.pixels : 0.0;
            printf("    [PLOT ERROR STATS] tile coverage: tiles=%zu fallback=%zu query_errors=%zu raw_candidates avg=%.1f max=%zu dedup avg=%.1f max=%zu tests_per_pixel=%.1f exact_fallback_pixels=%zu\n",
                   ts.tiles, ts.fallback_tiles, ts.query_errors,
                   avg_raw, ts.candidate_max, avg_dedup,
                   ts.dedup_candidate_max, tests_per_pixel,
                   ts.exact_fallback_pixels);
            print_point_coverage_stats("tile fallback");
        }
    }

    if (p->show_errors &&
        (full_exact_coverage || path_exact_coverage || tile_exact_coverage)) {
        apply_boundary_ambiguity_filter(
            sys, &view, p->width, p->height, cell_ids, secondary_ids,
            coverage, errors, plot_error_stats, "after exact coverage");
    }

    if (plot_error_stats) {
        int c0, c1, c2, eo, eu;
        count_plot_error_state(cell_ids, coverage, errors, num_pixels,
                               &c0, &c1, &c2, &eo, &eu);
        printf("    [PLOT ERROR STATS] after grid: coverage none=%d one=%d multi=%d, errors overlap=%d undefined=%d\n",
               c0, c1, c2, eo, eu);
    }

    /* Select ID array based on color mode */
    const int* color_ids = (p->color_mode == COLOR_BY_MATERIAL) ? material_ids : cell_ids;

    /* Color pixels (grid row 0 = v_min = bottom of image) */
    for (int y = 0; y < p->height; y++) {
        int gy = p->height - 1 - y;  /* pixel top -> grid bottom */
        for (int x = 0; x < p->width; x++) {
            int id = color_ids[gy * p->width + x];
            int mat = material_ids[gy * p->width + x];
            int pidx = (y * p->width + x) * 3;
            if (id < 0 || mat == 0) {
                pixels[pidx + 0] = 255;
                pixels[pidx + 1] = 255;
                pixels[pidx + 2] = 255;
            } else {
                int color_idx = id % PALETTE_SIZE;
                pixels[pidx + 0] = palette[color_idx][0];
                pixels[pidx + 1] = palette[color_idx][1];
                pixels[pidx + 2] = palette[color_idx][2];
            }
        }
    }

    /* Detect nested overlaps using curves if error display is on */
    alea_slice_curves_t* err_curves = NULL;
    if (p->show_errors) {
        t0 = get_time_ms();
        err_curves = get_curves_for_plot(sys, p);
        if (err_curves) {
            int refined = 0;
            if (selective_exact_coverage && !full_exact_coverage &&
                !tile_exact_coverage && !path_exact_coverage) {
                refined = alea_refine_grid_coverage_curves_exact(
                    sys, &view, err_curves, p->width, p->height,
                    -1, secondary_ids, coverage, errors);
                if (verbose) {
                    printf("    Selective exact coverage: %d pixels refined\n",
                           refined);
                }
                if (plot_error_stats) {
                    print_point_coverage_stats("selective exact");
                }
            }
            int new_overlaps = 0;
            if (selective_exact_coverage && !full_exact_coverage &&
                !tile_exact_coverage && !path_exact_coverage) {
                new_overlaps = alea_check_grid_overlaps_curves_coverage(
                    &view, err_curves, p->width, p->height,
                    coverage, errors);
            } else if (!full_exact_coverage && !tile_exact_coverage &&
                       !path_exact_coverage) {
                new_overlaps = alea_check_grid_overlaps_curves(
                    sys, &view, err_curves, p->width, p->height,
                    -1, cell_ids, errors);
            }
            t1 = get_time_ms();
            if (verbose) {
                printf("    Nested overlap check: %d new overlap pixels (%.1f ms)\n",
                       new_overlaps, t1 - t0);
            }
            for (int i = 0; i < num_pixels; i++) {
                if (errors[i] == GRID_ERR_OVERLAP) {
                    coverage[i] = ALEA_COVERAGE_MULTI;
                } else if (errors[i] == GRID_ERR_UNDEFINED || cell_ids[i] < 0) {
                    coverage[i] = ALEA_COVERAGE_NONE;
                }
            }
            if (selective_exact_coverage || full_exact_coverage ||
                tile_exact_coverage || path_exact_coverage) {
                apply_boundary_ambiguity_filter(
                    sys, &view, p->width, p->height, cell_ids,
                    secondary_ids, coverage, errors, plot_error_stats,
                    "after curve check");
            }
            if (plot_error_stats) {
                int c0, c1, c2, eo, eu;
                count_plot_error_state(cell_ids, coverage, errors, num_pixels,
                                       &c0, &c1, &c2, &eo, &eu);
                printf("    [PLOT ERROR STATS] after curve check: coverage none=%d one=%d multi=%d, errors overlap=%d undefined=%d, refined=%d, curve_new=%d\n",
                       c0, c1, c2, eo, eu, refined, new_overlaps);
            }
        }
    }

    if (plot_error_stats) {
        alea_plot_error_component_result_t* comps =
            alea_classify_plot_error_components(
                cell_ids, secondary_ids, coverage, p->width, p->height);
        if (comps) {
            size_t undefined = 0, partial = 0, total = 0;
            for (size_t i = 0; i < comps->component_count; i++) {
                const alea_plot_error_component_t* c = &comps->components[i];
                switch (c->kind) {
                    case ALEA_PLOT_ERR_UNDEFINED_REGION: undefined++; break;
                    case ALEA_PLOT_ERR_PARTIAL_OVERLAP: partial++; break;
                    case ALEA_PLOT_ERR_TOTAL_OVERLAP: total++; break;
                }
            }
            printf("    [PLOT ERROR STATS] components undefined=%zu partial_overlap=%zu total_overlap=%zu\n",
                   undefined, partial, total);
            alea_plot_error_components_free(comps);
        }
    }

    /* Draw contours (grid is in math order, pixels in image order) */
    const int* boundary_ids = (p->contour_mode == CONTOUR_BY_MATERIAL) ? material_ids : cell_ids;
    draw_contours_ex(pixels, boundary_ids, errors, p->width, p->height);

    /* Draw geometry-validator boundary diagnostics if requested */
    if (p->show_errors) {
        t0 = get_time_ms();
        size_t error_event_count = 0;

        alea_geom_validator_options_t vopts;
        alea_geom_validator_options_init(&vopts);
        /* Preserve current --errors behaviour: outer boundaries to void are
         * still flagged (do not allow exterior void). Depth -1 matches the
         * grid coverage query above. */
        vopts.universe_depth = -1;
        vopts.max_errors    = (size_t)p->width * (size_t)p->height; /* don't truncate the overlay */
        vopts.max_crossings = 0;  /* 0 -> validator default cap */

        alea_geom_validator_result_t vres;
        alea_geom_validator_result_init(&vres);
        if (err_curves &&
            alea_validate_geometry_slice(sys, &view, err_curves,
                                         &vopts, &vres) == 0) {
            error_event_count = vres.error_count;
            if (verbose) {
                printf("    Error check: %zu boundary events found\n",
                       vres.error_count);
            }
            if (vres.error_count > 0) {
                draw_validator_error_dots(pixels, p->width, p->height, p,
                                          &view, &vres);
            }
        }
        alea_geom_validator_result_free(&vres);

        t1 = get_time_ms();
        if (verbose) {
            printf("    Error overlay: %.1f ms\n", t1 - t0);
        }
        if (plot_error_stats) {
            printf("    [PLOT ERROR STATS] boundary events=%zu\n",
                   error_event_count);
        }
    }
    if (err_curves) alea_slice_curves_free(err_curves);

    if (verify_interval > 0) {
        t0 = get_time_ms();
        int checked = 0, mismatches = 0;
        int verify_rc = verify_plot_coverage_sample(
            sys, &view, p->width, p->height, -1, coverage,
            verify_interval, &checked, &mismatches);
        t1 = get_time_ms();
        printf("    [PLOT ERROR VERIFY] interval=%d checked=%d mismatches=%d (%.1f ms)\n",
               verify_interval, checked, mismatches, t1 - t0);
        if (verify_rc != 0) {
            fprintf(stderr, "Error: Coverage verifier found %d mismatches\n",
                    mismatches);
            free(cell_ids);
            free(material_ids);
            free(secondary_ids);
            free(path_ids);
            alea_slice_path_table_free(&path_table);
            free(errors);
            free(coverage);
            free(pixels);
            return -1;
        }
    }

    /* Draw labels if requested */
    if (p->labels.show_cells) {
        draw_cell_labels(pixels, p->width, p->height, cell_ids);
    }
    if (p->labels.show_materials) {
        draw_material_labels(pixels, p->width, p->height, cell_ids, material_ids);
    }
    if (p->labels.show_surfaces) {
        alea_slice_surface_boundary_map_t* boundary_map = NULL;
        alea_slice_classify_point_fn classify =
            (p->contour_mode == CONTOUR_BY_MATERIAL)
                ? alea_slice_classify_material : alea_slice_classify_cell;
        if (alea_slice_surface_boundary_map_create(
                sys, &view, p->width, p->height, boundary_ids,
                classify, NULL, &boundary_map) == 0 && boundary_map) {
            alea_label_position_t* surf_labels = NULL;
            int surf_count = 0;
            if (alea_find_surface_labels_on_boundary_map(
                    boundary_map, 20, &surf_labels, &surf_count) == 0 &&
                surf_labels) {
                for (int i = 0; i < surf_count; i++) {
                    char label[32];
                    snprintf(label, sizeof(label), "S%d", surf_labels[i].id);
                    int py = p->height - 1 - surf_labels[i].py;
                    draw_string_centered(pixels, p->width, p->height,
                                        surf_labels[i].px, py,
                                        label, 100, 150, 255);
                }
                free(surf_labels);
            }
            alea_slice_surface_boundary_map_free(boundary_map);
        }
    }

    /* Draw tick labels if requested */
    if (p->show_ticks) {
        const char* u_label;
        const char* v_label;
        get_axis_labels(p->type, &u_label, &v_label);
        draw_ticks(pixels, p->width, p->height,
                   p->u_min, p->u_max, p->v_min, p->v_max,
                   u_label, v_label);
    }

    /* Write output */
    rc = write_image(p->output, pixels, p->width, p->height);

    free(cell_ids);
    free(material_ids);
    free(secondary_ids);
    free(path_ids);
    alea_slice_path_table_free(&path_table);
    free(errors);
    free(coverage);
    free(pixels);

    if (rc != 0) {
        fprintf(stderr, "Error: Failed to write %s\n", p->output);
        return -1;
    }

    return 0;
}

/* ============================================================================
 * BATCH FILE PARSING
 * ============================================================================
 *
 * Batch file format (one plot per line):
 *   Z value u_min u_max v_min v_max WxH output.png [options...]
 *   Y value u_min u_max v_min v_max WxH output.png [options...]
 *   X value u_min u_max v_min v_max WxH output.png [options...]
 *   PLANE  ox oy oz nx ny nz ux uy uz u_min u_max v_min v_max WxH output [options...]
 *   PLANE3 x1 y1 z1 x2 y2 z2 x3 y3 z3 u_min u_max v_min v_max WxH output [options...]
 *
 * Resolution can be: 800 (square) or 800x600 (width x height)
 *
 * Options:
 *   labels=cells,materials,surfaces  - Show labels
 *   color=cells|materials            - Color by cell or material ID
 *   ticks                            - Show axis tick labels
 *   errors                           - Show analytical error lines
 *
 * Lines starting with # are comments
 * Empty lines are skipped
 */

/* Parse resolution string (e.g., "800" or "800x600") */
static void parse_resolution(const char* str, int* width, int* height) {
    int w, h;
    if (sscanf(str, "%dx%d", &w, &h) == 2) {
        *width = w;
        *height = h;
    } else if (sscanf(str, "%d", &w) == 1) {
        *width = w;
        *height = w;  /* Square */
    } else {
        *width = 800;
        *height = 800;
    }
}

/* Parse batch line options (labels, color, ticks) */
static void parse_options(const char* line, plot_params_t* p) {
    /* Look for labels= */
    const char* labels_ptr = strstr(line, "labels=");
    if (labels_ptr) {
        char buf[256];
        if (sscanf(labels_ptr + 7, "%255[^, \t\n]", buf) == 1) {
            p->labels = parse_label_string(buf);
        }
    }

    /* Look for color= */
    const char* color_ptr = strstr(line, "color=");
    if (color_ptr) {
        char buf[32];
        if (sscanf(color_ptr + 6, "%31[^, \t\n]", buf) == 1) {
            if (strcasecmp(buf, "materials") == 0 || strcasecmp(buf, "material") == 0) {
                p->color_mode = COLOR_BY_MATERIAL;
            } else {
                p->color_mode = COLOR_BY_CELL;
            }
        }
    }

    /* Look for contours= */
    const char* contour_ptr = strstr(line, "contours=");
    if (contour_ptr) {
        char buf[32];
        if (sscanf(contour_ptr + 9, "%31[^, \t\n]", buf) == 1) {
            if (strcasecmp(buf, "materials") == 0 || strcasecmp(buf, "material") == 0) {
                p->contour_mode = CONTOUR_BY_MATERIAL;
            } else {
                p->contour_mode = CONTOUR_BY_CELL;
            }
        }
    }

    /* Look for ticks */
    if (strstr(line, "ticks")) {
        p->show_ticks = 1;
    }

    /* Look for errors */
    if (strstr(line, "errors")) {
        p->show_errors = 1;
    }
}

/* Compute plane parameters from 3 points */
static void plane_from_3_points(plot_params_t* p,
                                 double x1, double y1, double z1,
                                 double x2, double y2, double z2,
                                 double x3, double y3, double z3) {
    /* Origin = P1 */
    p->origin[0] = x1;
    p->origin[1] = y1;
    p->origin[2] = z1;

    /* U = P2 - P1 (will be normalized by slice_plane_init) */
    double ux = x2 - x1;
    double uy = y2 - y1;
    double uz = z2 - z1;

    /* V' = P3 - P1 */
    double vx = x3 - x1;
    double vy = y3 - y1;
    double vz = z3 - z1;

    /* Normal = U x V' */
    p->normal[0] = uy * vz - uz * vy;
    p->normal[1] = uz * vx - ux * vz;
    p->normal[2] = ux * vy - uy * vx;

    /* Up = use V' direction (will be orthogonalized by slice_plane_init) */
    p->up[0] = vx;
    p->up[1] = vy;
    p->up[2] = vz;
}

static int parse_batch_line(const char* line, plot_params_t* p) {
    /* Skip leading whitespace */
    while (*line == ' ' || *line == '\t') line++;

    /* Skip comments and empty lines */
    if (*line == '#' || *line == '\n' || *line == '\0') return 0;

    /* Initialize defaults */
    memset(p, 0, sizeof(*p));
    p->width = 800;
    p->height = 800;

    char type_str[16] = {0};
    char res_str[32] = {0};
    char output[512] = {0};
    int n;

    /* Try to parse the type prefix */
    n = sscanf(line, "%15s", type_str);
    if (n < 1) return -1;

    /* Determine slice type */
    if (strcasecmp(type_str, "Z") == 0) {
        p->type = SLICE_Z;
        n = sscanf(line, "%*s %lf %lf %lf %lf %lf %31s %511s",
                   &p->value, &p->u_min, &p->u_max, &p->v_min, &p->v_max,
                   res_str, output);
        if (n < 7) return -1;
    }
    else if (strcasecmp(type_str, "Y") == 0) {
        p->type = SLICE_Y;
        n = sscanf(line, "%*s %lf %lf %lf %lf %lf %31s %511s",
                   &p->value, &p->u_min, &p->u_max, &p->v_min, &p->v_max,
                   res_str, output);
        if (n < 7) return -1;
    }
    else if (strcasecmp(type_str, "X") == 0) {
        p->type = SLICE_X;
        n = sscanf(line, "%*s %lf %lf %lf %lf %lf %31s %511s",
                   &p->value, &p->u_min, &p->u_max, &p->v_min, &p->v_max,
                   res_str, output);
        if (n < 7) return -1;
    }
    else if (strcasecmp(type_str, "PLANE") == 0) {
        /* Plane defined by origin + normal + up */
        p->type = SLICE_PLANE;
        n = sscanf(line, "%*s %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %31s %511s",
                   &p->origin[0], &p->origin[1], &p->origin[2],
                   &p->normal[0], &p->normal[1], &p->normal[2],
                   &p->up[0], &p->up[1], &p->up[2],
                   &p->u_min, &p->u_max, &p->v_min, &p->v_max,
                   res_str, output);
        if (n < 15) return -1;
    }
    else if (strcasecmp(type_str, "PLANE3") == 0) {
        /* Plane defined by 3 points */
        p->type = SLICE_PLANE;
        double x1, y1, z1, x2, y2, z2, x3, y3, z3;
        n = sscanf(line, "%*s %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %31s %511s",
                   &x1, &y1, &z1, &x2, &y2, &z2, &x3, &y3, &z3,
                   &p->u_min, &p->u_max, &p->v_min, &p->v_max,
                   res_str, output);
        if (n < 15) return -1;
        plane_from_3_points(p, x1, y1, z1, x2, y2, z2, x3, y3, z3);
    }
    else {
        /* Legacy format: assume Z slice with just numbers */
        p->type = SLICE_Z;
        n = sscanf(line, "%lf %lf %lf %lf %lf %31s %511s",
                   &p->value, &p->u_min, &p->u_max, &p->v_min, &p->v_max,
                   res_str, output);
        if (n < 7) return -1;
    }

    /* Parse resolution */
    parse_resolution(res_str, &p->width, &p->height);

    strncpy(p->output, output, sizeof(p->output) - 1);
    p->output[sizeof(p->output) - 1] = '\0';

    /* Parse options (labels, color, ticks) from rest of line */
    parse_options(line, p);

    return 1;  /* Valid plot line */
}

static int run_batch(const char* batch_file,
                     const char* input_file) {
    FILE* f = fopen(batch_file, "r");
    if (!f) {
        fprintf(stderr, "Error: Cannot open batch file: %s\n", batch_file);
        return 1;
    }

    printf("CSG Slice Plotter - Batch Mode\n");
    printf("==========================================\n");
    printf("Input:      %s\n", input_file);
    printf("Batch file: %s\n", batch_file);
    printf("\n");

    /* Load geometry */
    printf("Loading geometry... ");
    fflush(stdout);
    double t0 = get_time_ms();

    loaded_model_t model;
    if (load_model(input_file, &model) != 0) {
        printf("FAILED: %s\n", alea_error());
        fclose(f);
        return 1;
    }
    alea_system_t* sys = model.sys;

    double t1 = get_time_ms();
    printf("OK (%.1f ms)\n", t1 - t0);
    printf("  Cells: %zu, Surfaces: %zu\n",
           alea_cell_count(sys), alea_surface_count(sys));

    /* Prepare query acceleration */
    printf("Preparing query acceleration... ");
    fflush(stdout);
    t0 = get_time_ms();

    if (alea_prepare_query_acceleration(sys) != 0) {
        printf("FAILED: %s\n", alea_error());
        destroy_model(&model);
        fclose(f);
        return 1;
    }

    t1 = get_time_ms();
    printf("OK (%.1f ms)\n\n", t1 - t0);

    /* Process batch file */
    printf("Processing plots:\n");
    char line[1024];
    int line_num = 0;
    int plot_count = 0;
    int error_count = 0;

    while (fgets(line, sizeof(line), f)) {
        line_num++;
        plot_params_t p;

        int rc = parse_batch_line(line, &p);
        if (rc == 0) continue;  /* Comment or empty */
        if (rc < 0) {
            fprintf(stderr, "  Line %d: Parse error\n", line_num);
            error_count++;
            continue;
        }

        plot_count++;
        if (render_plot(sys, &p, 1) != 0) {
            error_count++;
        }
    }

    fclose(f);
    destroy_model(&model);

    printf("\nBatch complete: %d plots, %d errors\n", plot_count, error_count);
    return error_count > 0 ? 1 : 0;
}

static void print_usage(const char* prog) {
    fprintf(stderr, "Usage: %s <input.i> <axis> <value> <u_min> <u_max> <v_min> <v_max> <WxH> [output.png] [options]\n", prog);
    fprintf(stderr, "       %s <input.i> --batch=<file.txt>\n", prog);
    fprintf(stderr, "\nAxis-aligned slices:\n");
    fprintf(stderr, "  axis = Z  -> XY plane at z=value, u=X, v=Y\n");
    fprintf(stderr, "  axis = Y  -> XZ plane at y=value, u=X, v=Z\n");
    fprintf(stderr, "  axis = X  -> YZ plane at x=value, u=Y, v=Z\n");
    fprintf(stderr, "\nResolution: 800 (square) or 800x600 (width x height)\n");
    fprintf(stderr, "\nOptions:\n");
    fprintf(stderr, "  --labels=cells          Show cell ID labels\n");
    fprintf(stderr, "  --labels=materials      Show material ID labels\n");
    fprintf(stderr, "  --labels=surfaces       Show surface ID labels\n");
    fprintf(stderr, "  --labels=all            Show all labels\n");
    fprintf(stderr, "  --color=cells           Color by cell ID (default)\n");
    fprintf(stderr, "  --color=materials       Color by material ID\n");
    fprintf(stderr, "  --contours=cells        Contours at cell boundaries (default)\n");
    fprintf(stderr, "  --contours=materials    Contours at material boundaries\n");
    fprintf(stderr, "  --ticks                 Show axis tick labels\n");
    fprintf(stderr, "  --debug                 Print verbose curve generation debug info\n");
    fprintf(stderr, "  --trace=px,py           Trace cell lookup at pixel (px,py) for debugging\n");
    fprintf(stderr, "\nBatch file format (one plot per line):\n");
    fprintf(stderr, "  Z value u_min u_max v_min v_max WxH output.png [options]\n");
    fprintf(stderr, "  Y value u_min u_max v_min v_max WxH output.png [options]\n");
    fprintf(stderr, "  X value u_min u_max v_min v_max WxH output.png [options]\n");
    fprintf(stderr, "  PLANE  ox oy oz nx ny nz ux uy uz u_min u_max v_min v_max WxH output [options]\n");
    fprintf(stderr, "  PLANE3 x1 y1 z1 x2 y2 z2 x3 y3 z3 u_min u_max v_min v_max WxH output [options]\n");
    fprintf(stderr, "\n  Batch options: labels=... color=cells|materials ticks\n");
    fprintf(stderr, "\nExamples:\n");
    fprintf(stderr, "  %s model.i Z 0 -50 50 -50 50 1000 slice.png\n", prog);
    fprintf(stderr, "  %s model.i Z 0 -50 50 -50 50 800x600 slice.png --ticks\n", prog);
    fprintf(stderr, "  %s model.i Y 0 -50 50 -50 50 1000 slice.png --color=materials --labels=materials\n", prog);
    fprintf(stderr, "  %s model.i --batch=plots.txt\n", prog);
}

int main(int argc, char** argv) {
    if (argc < 3) {
        print_usage(argv[0]);
        return 1;
    }

    const char* input_file = argv[1];

    /* Check for batch mode */
    for (int i = 2; i < argc; i++) {
        if (strncmp(argv[i], "--batch=", 8) == 0) {
            return run_batch(argv[i] + 8, input_file);
        }
    }

    /* Single plot mode - need at least 9 args (prog input axis value u_min u_max v_min v_max res) */
    if (argc < 9) {
        print_usage(argv[0]);
        return 1;
    }

    /* Parse axis type */
    plot_params_t plot;
    memset(&plot, 0, sizeof(plot));

    const char* axis_str = argv[2];
    if (strcasecmp(axis_str, "Z") == 0) {
        plot.type = SLICE_Z;
    } else if (strcasecmp(axis_str, "Y") == 0) {
        plot.type = SLICE_Y;
    } else if (strcasecmp(axis_str, "X") == 0) {
        plot.type = SLICE_X;
    } else {
        fprintf(stderr, "Error: Invalid axis '%s'. Use Z, Y, or X.\n", axis_str);
        return 1;
    }

    plot.value = atof(argv[3]);
    plot.u_min = atof(argv[4]);
    plot.u_max = atof(argv[5]);
    plot.v_min = atof(argv[6]);
    plot.v_max = atof(argv[7]);
    parse_resolution(argv[8], &plot.width, &plot.height);

    /* Find output file (first non-option argument after resolution) */
    const char* output_file = "plot.png";
    for (int i = 9; i < argc; i++) {
        if (strncmp(argv[i], "--", 2) != 0) {
            output_file = argv[i];
            break;
        }
    }
    strncpy(plot.output, output_file, sizeof(plot.output) - 1);

    /* Parse options */
    int debug_curves = 0;
    int trace_px = -1, trace_py = -1;  /* Pixel to trace (-1 = none) */
    for (int i = 9; i < argc; i++) {
        if (strncmp(argv[i], "--labels=", 9) == 0) {
            plot.labels = parse_label_string(argv[i] + 9);
        } else if (strncmp(argv[i], "--color=", 8) == 0) {
            const char* mode = argv[i] + 8;
            if (strcasecmp(mode, "materials") == 0 || strcasecmp(mode, "material") == 0) {
                plot.color_mode = COLOR_BY_MATERIAL;
            }
        } else if (strncmp(argv[i], "--contours=", 11) == 0) {
            const char* mode = argv[i] + 11;
            if (strcasecmp(mode, "materials") == 0 || strcasecmp(mode, "material") == 0) {
                plot.contour_mode = CONTOUR_BY_MATERIAL;
            } else {
                plot.contour_mode = CONTOUR_BY_CELL;
            }
        } else if (strcmp(argv[i], "--ticks") == 0) {
            plot.show_ticks = 1;
        } else if (strcmp(argv[i], "--errors") == 0) {
            plot.show_errors = 1;
        } else if (strcmp(argv[i], "--debug") == 0) {
            debug_curves = 1;
        } else if (strncmp(argv[i], "--trace=", 8) == 0) {
            sscanf(argv[i] + 8, "%d,%d", &trace_px, &trace_py);
        }
    }

    const char* axis_names[] = {"Z", "Y", "X"};
    printf("CSG Slice Plotter\n");
    printf("==========================================\n");
    printf("Input:      %s\n", input_file);
    printf("%s plane:   %.4g\n", axis_names[plot.type], plot.value);
    printf("Bounds:     [%.4g, %.4g] x [%.4g, %.4g]\n", plot.u_min, plot.u_max, plot.v_min, plot.v_max);
    printf("Resolution: %d x %d\n", plot.width, plot.height);
    printf("Output:     %s\n", output_file);
    printf("Color by:   %s\n", plot.color_mode == COLOR_BY_MATERIAL ? "materials" : "cells");
    printf("Contours:   %s\n", plot.contour_mode == CONTOUR_BY_MATERIAL ? "materials" : "cells");
    if (plot.show_ticks) printf("Ticks:      enabled\n");
    if (plot.show_errors) printf("Errors:     enabled\n");
    if (plot.labels.show_cells || plot.labels.show_materials || plot.labels.show_surfaces) {
        printf("Labels:     %s%s%s\n",
               plot.labels.show_cells ? "cells " : "",
               plot.labels.show_materials ? "materials " : "",
               plot.labels.show_surfaces ? "surfaces " : "");
    }
    printf("\n");

    /* Load geometry */
    printf("Loading geometry... ");
    fflush(stdout);
    double t0 = get_time_ms();

    loaded_model_t model;
    if (load_model(input_file, &model) != 0) {
        printf("FAILED: %s\n", alea_error());
        return 1;
    }
    alea_system_t* sys = model.sys;

    double t1 = get_time_ms();
    printf("OK (%.1f ms)\n", t1 - t0);
    printf("  Cells: %zu, Surfaces: %zu\n",
           alea_cell_count(sys), alea_surface_count(sys));
    printf("  [RSS after load: %.1f MB]\n", get_rss_mb());
    fflush(stdout);

    /* Prepare query acceleration */
    printf("Preparing query acceleration... ");
    fflush(stdout);
    t0 = get_time_ms();

    if (alea_prepare_query_acceleration(sys) != 0) {
        printf("FAILED: %s\n", alea_error());
        destroy_model(&model);
        return 1;
    }

    t1 = get_time_ms();
    printf("OK (%.1f ms)\n", t1 - t0);
    printf("  [RSS after query acceleration: %.1f MB]\n\n", get_rss_mb());

    /* Enable curve debug output if requested */
    if (debug_curves) {
        printf("Debug mode enabled - verbose curve output:\n");
        printf("==========================================\n");
        alea_slice_curve_set_debug(1);
    }

    /* Trace a specific pixel if requested */
    if (trace_px >= 0 && trace_py >= 0) {
        printf("Tracing pixel (%d, %d):\n", trace_px, trace_py);
        printf("==========================================\n");

        /* Compute world coordinates for the pixel */
        double du = (plot.u_max - plot.u_min) / plot.width;
        double dv = (plot.v_max - plot.v_min) / plot.height;
        double u = plot.u_min + (trace_px + 0.5) * du;
        double v = plot.v_max - (trace_py + 0.5) * dv;

        printf("  UV coordinates: (%.6f, %.6f)\n", u, v);

        /* Enable point trace and do a single lookup */
        alea_slice_point_trace_set_debug(1);

        int cell_id = -1, material_id = 0;
        alea_slice_view_t trace_view;
        build_view_from_params(&plot, &trace_view);
        trace_view.u_min = u - du/2;
        trace_view.u_max = u + du/2;
        trace_view.v_min = v - dv/2;
        trace_view.v_max = v + dv/2;
        alea_find_cells_grid(sys, &trace_view, 1, 1,
                                  -1, &cell_id, &material_id, NULL);

        alea_slice_point_trace_set_debug(0);
        printf("  Result: cell_id=%d, material_id=%d\n\n", cell_id, material_id);
    }

    /* Render the plot */
    printf("Rendering...\n");
    int rc = render_plot(sys, &plot, 1);

    destroy_model(&model);

    if (rc == 0) {
        printf("\nDone: %s\n", output_file);
    }
    return rc == 0 ? 0 : 1;
}
