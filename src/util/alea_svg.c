// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/**
 * @file alea_svg.c
 * @brief SVG plot generator — log-log and lin-lin scientific plots
 */

#include "alea_svg.h"
#include <math.h>
#include <string.h>
#include <float.h>

/* Layout constants */
#define MARGIN_LEFT   75
#define MARGIN_RIGHT  25
#define MARGIN_TOP    40
#define MARGIN_BOTTOM 50
#define LEGEND_LINE   18

/* Default colors for curves */
static const char* default_colors[] = {
    "#2266cc", "#cc2222", "#22aa44", "#cc8800",
    "#8844cc", "#00aaaa", "#cc44aa", "#666666",
    "#4488ff", "#ff6644", "#44cc88", "#ffaa22",
    "#aa66ff", "#44cccc", "#ff66cc", "#999999",
};

static double log10_safe(double x) {
    return (x > 0.0) ? log10(x) : -30.0;
}

/* Round down to nearest power of 10 */
static double floor_pow10(double x) {
    if (x <= 0.0) return 1e-30;
    return pow(10.0, floor(log10(x)));
}

/* Round up to nearest power of 10 */
static double ceil_pow10(double x) {
    if (x <= 0.0) return 1e-30;
    return pow(10.0, ceil(log10(x)));
}

/* Compute auto-range from all curves */
static void compute_range(const alea_svg_plot_t* p,
                           double* xmin, double* xmax,
                           double* ymin, double* ymax) {
    double x0 = DBL_MAX, x1 = -DBL_MAX;
    double y0 = DBL_MAX, y1 = -DBL_MAX;

    for (int c = 0; c < p->n_curves; c++) {
        const alea_svg_curve_t* cv = &p->curves[c];
        for (int i = 0; i < cv->n; i++) {
            double vx = cv->x[i], vy = cv->y[i];
            if (p->x_scale == ALEA_SVG_SCALE_LOG && vx <= 0.0) continue;
            if (p->y_scale == ALEA_SVG_SCALE_LOG && vy <= 0.0) continue;
            if (vx < x0) x0 = vx;
            if (vx > x1) x1 = vx;
            if (vy < y0) y0 = vy;
            if (vy > y1) y1 = vy;
        }
    }

    if (p->x_scale == ALEA_SVG_SCALE_LOG) {
        *xmin = (p->x_min > 0.0) ? p->x_min : floor_pow10(x0);
        *xmax = (p->x_max > 0.0) ? p->x_max : ceil_pow10(x1);
    } else {
        *xmin = (p->x_min != 0.0 || p->x_max != 0.0) ? p->x_min : x0;
        *xmax = (p->x_min != 0.0 || p->x_max != 0.0) ? p->x_max : x1;
    }

    if (p->y_scale == ALEA_SVG_SCALE_LOG) {
        *ymin = (p->y_min > 0.0) ? p->y_min : floor_pow10(y0);
        *ymax = (p->y_max > 0.0) ? p->y_max : ceil_pow10(y1);
    } else {
        *ymin = (p->y_min != 0.0 || p->y_max != 0.0) ? p->y_min : y0;
        *ymax = (p->y_min != 0.0 || p->y_max != 0.0) ? p->y_max : y1;
        /* Add 5% padding so curves don't touch the edges */
        double margin = (*ymax - *ymin) * 0.05;
        if (margin > 0.0) {
            *ymin -= margin;
            *ymax += margin;
        }
    }
}

/* Map data value to pixel coordinate */
static double map_x(double val, double xmin, double xmax,
                     int plot_w, alea_svg_scale_t scale) {
    double t;
    if (scale == ALEA_SVG_SCALE_LOG) {
        t = (log10_safe(val) - log10_safe(xmin)) /
            (log10_safe(xmax) - log10_safe(xmin));
    } else {
        t = (val - xmin) / (xmax - xmin);
    }
    return MARGIN_LEFT + t * plot_w;
}

static double map_y(double val, double ymin, double ymax,
                     int plot_h, alea_svg_scale_t scale) {
    double t;
    if (scale == ALEA_SVG_SCALE_LOG) {
        t = (log10_safe(val) - log10_safe(ymin)) /
            (log10_safe(ymax) - log10_safe(ymin));
    } else {
        t = (val - ymin) / (ymax - ymin);
    }
    return MARGIN_TOP + plot_h - t * plot_h;
}

/* Format tick label for log scale */
static void format_log_label(char* buf, size_t sz, double val) {
    double e = log10(val);
    int ie = (int)round(e);
    if (fabs(e - ie) < 0.01) {
        if (ie == 0) snprintf(buf, sz, "1");
        else if (ie == 1) snprintf(buf, sz, "10");
        else if (ie == -1) snprintf(buf, sz, "0.1");
        else snprintf(buf, sz, "10<tspan font-size=\"9\" baseline-shift=\"super\">%d</tspan>", ie);
    } else {
        snprintf(buf, sz, "%.2g", val);
    }
}

/* Format tick label for linear scale */
static void format_lin_label(char* buf, size_t sz, double val) {
    if (fabs(val) < 1e-15) snprintf(buf, sz, "0");
    else if (fabs(val) >= 1000.0 || fabs(val) < 0.01)
        snprintf(buf, sz, "%.2g", val);
    else
        snprintf(buf, sz, "%.4g", val);
}

int alea_svg_plot_fwrite(const alea_svg_plot_t* p, FILE* fp) {
    if (!p || !fp || p->n_curves == 0) return -1;

    int W = p->width > 0 ? p->width : 900;
    int H = p->height > 0 ? p->height : 500;
    int legend_h = p->n_curves > 1 ? p->n_curves * LEGEND_LINE + 10 : 0;
    int plot_w = W - MARGIN_LEFT - MARGIN_RIGHT;
    int plot_h = H - MARGIN_TOP - MARGIN_BOTTOM - (legend_h > 0 ? 0 : 0);

    double xmin, xmax, ymin, ymax;
    compute_range(p, &xmin, &xmax, &ymin, &ymax);

    /* SVG header */
    fprintf(fp, "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n");
    fprintf(fp, "<svg xmlns=\"http://www.w3.org/2000/svg\" "
            "width=\"%d\" height=\"%d\" viewBox=\"0 0 %d %d\">\n",
            W, H + legend_h, W, H + legend_h);

    /* Background */
    fprintf(fp, "<rect width=\"%d\" height=\"%d\" fill=\"white\"/>\n",
            W, H + legend_h);

    /* Style */
    fprintf(fp, "<style>\n"
            "  text { font-family: 'Helvetica Neue', Arial, sans-serif; }\n"
            "  .tick { font-size: 11px; fill: #333; }\n"
            "  .label { font-size: 13px; fill: #222; font-weight: 500; }\n"
            "  .title { font-size: 15px; fill: #111; font-weight: 600; }\n"
            "  .legend { font-size: 12px; fill: #333; }\n"
            "  .grid { stroke: #e0e0e0; stroke-width: 0.5; }\n"
            "</style>\n");

    /* Plot area clip path */
    fprintf(fp, "<defs><clipPath id=\"plot-area\">"
            "<rect x=\"%d\" y=\"%d\" width=\"%d\" height=\"%d\"/>"
            "</clipPath></defs>\n",
            MARGIN_LEFT, MARGIN_TOP, plot_w, plot_h);

    /* Plot area border */
    fprintf(fp, "<rect x=\"%d\" y=\"%d\" width=\"%d\" height=\"%d\" "
            "fill=\"none\" stroke=\"#888\" stroke-width=\"1\"/>\n",
            MARGIN_LEFT, MARGIN_TOP, plot_w, plot_h);

    /* Grid lines and tick labels — X axis */
    if (p->x_scale == ALEA_SVG_SCALE_LOG) {
        double decade = floor_pow10(xmin);
        while (decade <= xmax * 1.01) {
            for (int m = 1; m <= 9; m++) {
                double val = decade * m;
                if (val < xmin * 0.999 || val > xmax * 1.001) continue;
                double px = map_x(val, xmin, xmax, plot_w, ALEA_SVG_SCALE_LOG);
                int is_major = (m == 1);
                fprintf(fp, "<line x1=\"%.1f\" y1=\"%d\" x2=\"%.1f\" y2=\"%d\" "
                        "class=\"grid\" %s/>\n",
                        px, MARGIN_TOP, px, MARGIN_TOP + plot_h,
                        is_major ? "stroke-width=\"0.8\"" : "stroke-dasharray=\"2,4\"");
                if (is_major) {
                    char lbl[128];
                    format_log_label(lbl, sizeof(lbl), val);
                    fprintf(fp, "<text x=\"%.1f\" y=\"%d\" text-anchor=\"middle\" "
                            "class=\"tick\">%s</text>\n",
                            px, MARGIN_TOP + plot_h + 16, lbl);
                }
            }
            decade *= 10.0;
        }
    } else {
        /* Linear: ~8 ticks */
        double range = xmax - xmin;
        double step = pow(10.0, floor(log10(range / 6.0)));
        if (range / step > 12) step *= 2;
        if (range / step < 4) step *= 0.5;
        double v = ceil(xmin / step) * step;
        while (v <= xmax * 1.001) {
            double px = map_x(v, xmin, xmax, plot_w, ALEA_SVG_SCALE_LINEAR);
            fprintf(fp, "<line x1=\"%.1f\" y1=\"%d\" x2=\"%.1f\" y2=\"%d\" class=\"grid\"/>\n",
                    px, MARGIN_TOP, px, MARGIN_TOP + plot_h);
            char lbl[32];
            format_lin_label(lbl, sizeof(lbl), v);
            fprintf(fp, "<text x=\"%.1f\" y=\"%d\" text-anchor=\"middle\" "
                    "class=\"tick\">%s</text>\n",
                    px, MARGIN_TOP + plot_h + 16, lbl);
            v += step;
        }
    }

    /* Grid lines and tick labels — Y axis */
    if (p->y_scale == ALEA_SVG_SCALE_LOG) {
        double decade = floor_pow10(ymin);
        while (decade <= ymax * 1.01) {
            for (int m = 1; m <= 9; m++) {
                double val = decade * m;
                if (val < ymin * 0.999 || val > ymax * 1.001) continue;
                double py = map_y(val, ymin, ymax, plot_h, ALEA_SVG_SCALE_LOG);
                int is_major = (m == 1);
                fprintf(fp, "<line x1=\"%d\" y1=\"%.1f\" x2=\"%d\" y2=\"%.1f\" "
                        "class=\"grid\" %s/>\n",
                        MARGIN_LEFT, py, MARGIN_LEFT + plot_w, py,
                        is_major ? "stroke-width=\"0.8\"" : "stroke-dasharray=\"2,4\"");
                if (is_major) {
                    char lbl[128];
                    format_log_label(lbl, sizeof(lbl), val);
                    fprintf(fp, "<text x=\"%d\" y=\"%.1f\" text-anchor=\"end\" "
                            "dominant-baseline=\"middle\" class=\"tick\">%s</text>\n",
                            MARGIN_LEFT - 6, py, lbl);
                }
            }
            decade *= 10.0;
        }
    } else {
        double range = ymax - ymin;
        double step = pow(10.0, floor(log10(range / 6.0)));
        if (range / step > 12) step *= 2;
        if (range / step < 4) step *= 0.5;
        double v = ceil(ymin / step) * step;
        while (v <= ymax * 1.001) {
            double py = map_y(v, ymin, ymax, plot_h, ALEA_SVG_SCALE_LINEAR);
            fprintf(fp, "<line x1=\"%d\" y1=\"%.1f\" x2=\"%d\" y2=\"%.1f\" class=\"grid\"/>\n",
                    MARGIN_LEFT, py, MARGIN_LEFT + plot_w, py);
            char lbl[32];
            format_lin_label(lbl, sizeof(lbl), v);
            fprintf(fp, "<text x=\"%d\" y=\"%.1f\" text-anchor=\"end\" "
                    "dominant-baseline=\"middle\" class=\"tick\">%s</text>\n",
                    MARGIN_LEFT - 6, py, lbl);
            v += step;
        }
    }

    /* X axis label */
    if (p->x_label) {
        fprintf(fp, "<text x=\"%d\" y=\"%d\" text-anchor=\"middle\" "
                "class=\"label\">%s</text>\n",
                MARGIN_LEFT + plot_w / 2, MARGIN_TOP + plot_h + 38, p->x_label);
    }

    /* Y axis label (rotated) */
    if (p->y_label) {
        fprintf(fp, "<text x=\"%d\" y=\"%d\" text-anchor=\"middle\" "
                "transform=\"rotate(-90,%d,%d)\" "
                "class=\"label\">%s</text>\n",
                16, MARGIN_TOP + plot_h / 2,
                16, MARGIN_TOP + plot_h / 2, p->y_label);
    }

    /* Title */
    if (p->title) {
        fprintf(fp, "<text x=\"%d\" y=\"%d\" text-anchor=\"middle\" "
                "class=\"title\">%s</text>\n",
                MARGIN_LEFT + plot_w / 2, MARGIN_TOP - 12, p->title);
    }

    /* Draw curves */
    for (int c = 0; c < p->n_curves; c++) {
        const alea_svg_curve_t* cv = &p->curves[c];
        const char* color = cv->color ? cv->color : default_colors[c % 16];
        double sw = cv->width > 0.0 ? cv->width : 1.5;

        fprintf(fp, "<polyline fill=\"none\" stroke=\"%s\" stroke-width=\"%.1f\" "
                "clip-path=\"url(#plot-area)\" ",
                color, sw);
        if (cv->dashed)
            fprintf(fp, "stroke-dasharray=\"6,3\" ");
        fprintf(fp, "points=\"");

        int first = 1;
        for (int i = 0; i < cv->n; i++) {
            double vx = cv->x[i], vy = cv->y[i];
            if (p->x_scale == ALEA_SVG_SCALE_LOG && vx <= 0.0) continue;
            if (p->y_scale == ALEA_SVG_SCALE_LOG && vy <= 0.0) continue;

            double px = map_x(vx, xmin, xmax, plot_w, p->x_scale);
            double py = map_y(vy, ymin, ymax, plot_h, p->y_scale);

            if (!first) fprintf(fp, " ");
            fprintf(fp, "%.1f,%.1f", px, py);
            first = 0;
        }
        fprintf(fp, "\"/>\n");
    }

    /* Legend */
    if (p->n_curves > 1) {
        int lx = MARGIN_LEFT + 12;
        int ly = H + 4;

        for (int c = 0; c < p->n_curves; c++) {
            const alea_svg_curve_t* cv = &p->curves[c];
            if (!cv->label) continue;
            const char* color = cv->color ? cv->color : default_colors[c % 16];

            fprintf(fp, "<line x1=\"%d\" y1=\"%d\" x2=\"%d\" y2=\"%d\" "
                    "stroke=\"%s\" stroke-width=\"2\" %s/>\n",
                    lx, ly + c * LEGEND_LINE + 9,
                    lx + 24, ly + c * LEGEND_LINE + 9,
                    color, cv->dashed ? "stroke-dasharray=\"6,3\"" : "");
            fprintf(fp, "<text x=\"%d\" y=\"%d\" dominant-baseline=\"middle\" "
                    "class=\"legend\">%s</text>\n",
                    lx + 30, ly + c * LEGEND_LINE + 9, cv->label);
        }
    }

    fprintf(fp, "</svg>\n");
    return 0;
}

int alea_svg_plot_write(const alea_svg_plot_t* p, const char* path) {
    FILE* fp = fopen(path, "w");
    if (!fp) return -1;
    int rc = alea_svg_plot_fwrite(p, fp);
    fclose(fp);
    return rc;
}
