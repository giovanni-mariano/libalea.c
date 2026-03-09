// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/**
 * @file alea_svg.h
 * @brief Minimal SVG plot generator for scientific data
 *
 * Generates log-log and lin-lin plots as standalone SVG files.
 * Zero external dependencies — outputs plain XML text.
 */

#ifndef ALEA_SVG_H
#define ALEA_SVG_H

#include <stdio.h>

/* Maximum curves per plot */
#define ALEA_SVG_MAX_CURVES 16

/* Plot axis scale */
typedef enum {
    ALEA_SVG_SCALE_LINEAR,
    ALEA_SVG_SCALE_LOG,
} alea_svg_scale_t;

/* One data curve */
typedef struct {
    const double* x;
    const double* y;
    int n;
    const char* label;
    const char* color;
    double width;       /* stroke width, 0 = default (1.5) */
    int dashed;         /* 1 = dashed line */
} alea_svg_curve_t;

/* Plot configuration */
typedef struct {
    /* Canvas */
    int width;          /* SVG width in px (0 = 900) */
    int height;         /* SVG height in px (0 = 500) */

    /* Axes */
    alea_svg_scale_t x_scale;
    alea_svg_scale_t y_scale;
    const char* x_label;
    const char* y_label;
    const char* title;

    /* Axis range (0,0 = auto) */
    double x_min, x_max;
    double y_min, y_max;

    /* Curves */
    alea_svg_curve_t curves[ALEA_SVG_MAX_CURVES];
    int n_curves;
} alea_svg_plot_t;

/**
 * @brief Initialize a plot with defaults
 */
static inline void alea_svg_plot_init(alea_svg_plot_t* p) {
    *p = (alea_svg_plot_t){
        .width = 900, .height = 500,
        .x_scale = ALEA_SVG_SCALE_LOG, .y_scale = ALEA_SVG_SCALE_LOG,
    };
}

/**
 * @brief Add a curve to the plot
 * @return 0 on success, -1 if max curves reached
 */
static inline int alea_svg_plot_add(alea_svg_plot_t* p, const double* x,
                                     const double* y, int n,
                                     const char* label, const char* color) {
    if (p->n_curves >= ALEA_SVG_MAX_CURVES) return -1;
    p->curves[p->n_curves++] = (alea_svg_curve_t){
        .x = x, .y = y, .n = n,
        .label = label, .color = color,
    };
    return 0;
}

/**
 * @brief Write the plot as SVG to a file
 * @return 0 on success, -1 on error
 */
int alea_svg_plot_write(const alea_svg_plot_t* p, const char* path);

/**
 * @brief Write the plot as SVG to an open FILE*
 */
int alea_svg_plot_fwrite(const alea_svg_plot_t* p, FILE* fp);

#endif /* ALEA_SVG_H */
