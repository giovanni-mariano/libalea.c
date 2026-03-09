// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/**
 * @file lua_svg.c
 * @brief Lua bindings for the SVG plot generator
 *
 * Provides the alea.SvgPlot userdata type for creating scientific plots.
 *
 * Example:
 *   local plot = alea.svg_plot()
 *   plot:set_title("My Plot")
 *   plot:set_x_label("Energy (MeV)")
 *   plot:set_y_label("Cross section (b)")
 *   plot:set_scale("log", "log")
 *   plot:add({1e-5, 1e-3, 1.0, 10.0},
 *            {100,   20,  5.0, 1.0},
 *            "Total", "#2266cc")
 *   plot:write("output.svg")
 */

#include "alea_lua.h"
#include "util/alea_svg.h"
#include "util/compat.h"
#include <stdlib.h>
#include <string.h>

#define ALEA_SVGPLOT_MT "alea.SvgPlot"

/* Each plot owns copies of x/y data and label/color strings */
typedef struct {
    alea_svg_plot_t plot;
    /* Owned copies of curve data */
    struct {
        double* x;
        double* y;
        char* label;
        char* color;
    } owned[ALEA_SVG_MAX_CURVES];
    /* Owned copies of label strings */
    char* x_label;
    char* y_label;
    char* title;
} lua_svg_plot_t;

static lua_svg_plot_t* check_svgplot(lua_State* L, int idx) {
    return (lua_svg_plot_t*)luaL_checkudata(L, idx, ALEA_SVGPLOT_MT);
}

/* alea.svg_plot() -> SvgPlot */
static int l_svg_plot_create(lua_State* L) {
    lua_svg_plot_t* ud = lua_newuserdata(L, sizeof(lua_svg_plot_t));
    memset(ud, 0, sizeof(*ud));
    alea_svg_plot_init(&ud->plot);
    luaL_setmetatable(L, ALEA_SVGPLOT_MT);
    return 1;
}

/* plot:set_title(title) */
static int l_svgplot_set_title(lua_State* L) {
    lua_svg_plot_t* ud = check_svgplot(L, 1);
    free(ud->title);
    ud->title = alea_strdup(luaL_checkstring(L, 2));
    ud->plot.title = ud->title;
    return 0;
}

/* plot:set_x_label(label) */
static int l_svgplot_set_x_label(lua_State* L) {
    lua_svg_plot_t* ud = check_svgplot(L, 1);
    free(ud->x_label);
    ud->x_label = alea_strdup(luaL_checkstring(L, 2));
    ud->plot.x_label = ud->x_label;
    return 0;
}

/* plot:set_y_label(label) */
static int l_svgplot_set_y_label(lua_State* L) {
    lua_svg_plot_t* ud = check_svgplot(L, 1);
    free(ud->y_label);
    ud->y_label = alea_strdup(luaL_checkstring(L, 2));
    ud->plot.y_label = ud->y_label;
    return 0;
}

/* plot:set_scale(x_scale, y_scale)  -- "log" or "linear" */
static int l_svgplot_set_scale(lua_State* L) {
    lua_svg_plot_t* ud = check_svgplot(L, 1);
    const char* xs = luaL_checkstring(L, 2);
    const char* ys = luaL_checkstring(L, 3);
    ud->plot.x_scale = (strcmp(xs, "log") == 0) ? ALEA_SVG_SCALE_LOG : ALEA_SVG_SCALE_LINEAR;
    ud->plot.y_scale = (strcmp(ys, "log") == 0) ? ALEA_SVG_SCALE_LOG : ALEA_SVG_SCALE_LINEAR;
    return 0;
}

/* plot:set_size(width, height) */
static int l_svgplot_set_size(lua_State* L) {
    lua_svg_plot_t* ud = check_svgplot(L, 1);
    ud->plot.width = (int)luaL_checkinteger(L, 2);
    ud->plot.height = (int)luaL_checkinteger(L, 3);
    return 0;
}

/* plot:set_range(x_min, x_max, y_min, y_max) -- 0,0 for auto */
static int l_svgplot_set_range(lua_State* L) {
    lua_svg_plot_t* ud = check_svgplot(L, 1);
    ud->plot.x_min = luaL_checknumber(L, 2);
    ud->plot.x_max = luaL_checknumber(L, 3);
    ud->plot.y_min = luaL_checknumber(L, 4);
    ud->plot.y_max = luaL_checknumber(L, 5);
    return 0;
}

/* plot:add(x_table, y_table, label, color [, opts])
 * opts = {width=2.0, dashed=true} */
static int l_svgplot_add(lua_State* L) {
    lua_svg_plot_t* ud = check_svgplot(L, 1);
    luaL_checktype(L, 2, LUA_TTABLE);
    luaL_checktype(L, 3, LUA_TTABLE);
    const char* label = luaL_checkstring(L, 4);
    const char* color = luaL_optstring(L, 5, NULL);

    int nx = (int)lua_rawlen(L, 2);
    int ny = (int)lua_rawlen(L, 3);
    if (nx != ny)
        return luaL_error(L, "x and y arrays must have same length (%d vs %d)", nx, ny);
    if (nx == 0)
        return luaL_error(L, "arrays must not be empty");

    int ci = ud->plot.n_curves;
    if (ci >= ALEA_SVG_MAX_CURVES)
        return luaL_error(L, "maximum %d curves per plot", ALEA_SVG_MAX_CURVES);

    /* Copy data arrays */
    double* x = malloc((size_t)nx * sizeof(double));
    double* y = malloc((size_t)nx * sizeof(double));
    if (!x || !y) { free(x); free(y); return luaL_error(L, "out of memory"); }

    for (int i = 0; i < nx; i++) {
        lua_rawgeti(L, 2, i + 1);
        x[i] = luaL_checknumber(L, -1);
        lua_pop(L, 1);
        lua_rawgeti(L, 3, i + 1);
        y[i] = luaL_checknumber(L, -1);
        lua_pop(L, 1);
    }

    ud->owned[ci].x = x;
    ud->owned[ci].y = y;
    ud->owned[ci].label = alea_strdup(label);
    ud->owned[ci].color = color ? alea_strdup(color) : NULL;

    alea_svg_curve_t* c = &ud->plot.curves[ci];
    c->x = x;
    c->y = y;
    c->n = nx;
    c->label = ud->owned[ci].label;
    c->color = ud->owned[ci].color;
    c->width = 0;
    c->dashed = 0;

    /* Optional opts table */
    if (lua_istable(L, 6)) {
        lua_getfield(L, 6, "width");
        if (lua_isnumber(L, -1)) c->width = lua_tonumber(L, -1);
        lua_pop(L, 1);
        lua_getfield(L, 6, "dashed");
        if (lua_isboolean(L, -1)) c->dashed = lua_toboolean(L, -1);
        lua_pop(L, 1);
    }

    ud->plot.n_curves++;
    return 0;
}

/* plot:write(path) */
static int l_svgplot_write(lua_State* L) {
    lua_svg_plot_t* ud = check_svgplot(L, 1);
    const char* path = luaL_checkstring(L, 2);
    if (alea_svg_plot_write(&ud->plot, path) != 0)
        return luaL_error(L, "failed to write SVG to: %s", path);
    return 0;
}

/* plot:n_curves() -> int */
static int l_svgplot_n_curves(lua_State* L) {
    lua_svg_plot_t* ud = check_svgplot(L, 1);
    lua_pushinteger(L, ud->plot.n_curves);
    return 1;
}

static int l_svgplot_gc(lua_State* L) {
    lua_svg_plot_t* ud = check_svgplot(L, 1);
    for (int i = 0; i < ud->plot.n_curves; i++) {
        free(ud->owned[i].x);
        free(ud->owned[i].y);
        free(ud->owned[i].label);
        free(ud->owned[i].color);
    }
    free(ud->x_label);
    free(ud->y_label);
    free(ud->title);
    return 0;
}

/* ============================================================================
 * Registration
 * ============================================================================ */

static const luaL_Reg svgplot_methods[] = {
    {"set_title",   l_svgplot_set_title},
    {"set_x_label", l_svgplot_set_x_label},
    {"set_y_label", l_svgplot_set_y_label},
    {"set_scale",   l_svgplot_set_scale},
    {"set_size",    l_svgplot_set_size},
    {"set_range",   l_svgplot_set_range},
    {"add",         l_svgplot_add},
    {"write",       l_svgplot_write},
    {"n_curves",    l_svgplot_n_curves},
    {NULL, NULL}
};

static const luaL_Reg svgplot_meta[] = {
    {"__gc", l_svgplot_gc},
    {NULL, NULL}
};

int luaopen_alea_svg(lua_State* L) {
    /* Register SvgPlot metatable */
    luaL_newmetatable(L, ALEA_SVGPLOT_MT);
    luaL_setfuncs(L, svgplot_meta, 0);
    lua_newtable(L);
    luaL_setfuncs(L, svgplot_methods, 0);
    lua_setfield(L, -2, "__index");
    lua_pop(L, 1);

    /* Add factory function to alea table */
    lua_pushcfunction(L, l_svg_plot_create);
    lua_setfield(L, -2, "svg_plot");

    return 0;
}
