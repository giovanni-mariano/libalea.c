// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

#include "alea_lua.h"
#include "alea_slice.h"
#include <stdlib.h>
#include <string.h>

/* ============================================================================
 * Userdata for curves
 * ============================================================================ */

typedef struct {
    alea_slice_curves_t* ptr;
} alea_lua_slice_curves_t;

static alea_lua_slice_curves_t* check_curves(lua_State* L, int idx) {
    return (alea_lua_slice_curves_t*)luaL_checkudata(L, idx, ALEA_CURVES_MT);
}

/* ============================================================================
 * Slice view helpers
 * ============================================================================ */

/* Push a slice_view_t as a Lua table */
static void push_slice_view(lua_State* L, const alea_slice_view_t* view) {
    lua_createtable(L, 0, 16);

    /* plane.origin */
    lua_pushnumber(L, view->plane.origin[0]); lua_setfield(L, -2, "ox");
    lua_pushnumber(L, view->plane.origin[1]); lua_setfield(L, -2, "oy");
    lua_pushnumber(L, view->plane.origin[2]); lua_setfield(L, -2, "oz");

    /* plane.normal */
    lua_pushnumber(L, view->plane.normal[0]); lua_setfield(L, -2, "nx");
    lua_pushnumber(L, view->plane.normal[1]); lua_setfield(L, -2, "ny");
    lua_pushnumber(L, view->plane.normal[2]); lua_setfield(L, -2, "nz");

    /* plane.u_axis */
    lua_pushnumber(L, view->plane.u_axis[0]); lua_setfield(L, -2, "ux");
    lua_pushnumber(L, view->plane.u_axis[1]); lua_setfield(L, -2, "uy");
    lua_pushnumber(L, view->plane.u_axis[2]); lua_setfield(L, -2, "uz");

    /* plane.v_axis */
    lua_pushnumber(L, view->plane.v_axis[0]); lua_setfield(L, -2, "vx");
    lua_pushnumber(L, view->plane.v_axis[1]); lua_setfield(L, -2, "vy");
    lua_pushnumber(L, view->plane.v_axis[2]); lua_setfield(L, -2, "vz");

    /* viewport */
    lua_pushnumber(L, view->u_min); lua_setfield(L, -2, "u_min");
    lua_pushnumber(L, view->u_max); lua_setfield(L, -2, "u_max");
    lua_pushnumber(L, view->v_min); lua_setfield(L, -2, "v_min");
    lua_pushnumber(L, view->v_max); lua_setfield(L, -2, "v_max");
}

/* Read a Lua table into a slice_view_t */
static void lua_to_slice_view(lua_State* L, int idx, alea_slice_view_t* view) {
    luaL_checktype(L, idx, LUA_TTABLE);

    lua_getfield(L, idx, "ox"); view->plane.origin[0] = luaL_checknumber(L, -1); lua_pop(L, 1);
    lua_getfield(L, idx, "oy"); view->plane.origin[1] = luaL_checknumber(L, -1); lua_pop(L, 1);
    lua_getfield(L, idx, "oz"); view->plane.origin[2] = luaL_checknumber(L, -1); lua_pop(L, 1);

    lua_getfield(L, idx, "nx"); view->plane.normal[0] = luaL_checknumber(L, -1); lua_pop(L, 1);
    lua_getfield(L, idx, "ny"); view->plane.normal[1] = luaL_checknumber(L, -1); lua_pop(L, 1);
    lua_getfield(L, idx, "nz"); view->plane.normal[2] = luaL_checknumber(L, -1); lua_pop(L, 1);

    lua_getfield(L, idx, "ux"); view->plane.u_axis[0] = luaL_checknumber(L, -1); lua_pop(L, 1);
    lua_getfield(L, idx, "uy"); view->plane.u_axis[1] = luaL_checknumber(L, -1); lua_pop(L, 1);
    lua_getfield(L, idx, "uz"); view->plane.u_axis[2] = luaL_checknumber(L, -1); lua_pop(L, 1);

    lua_getfield(L, idx, "vx"); view->plane.v_axis[0] = luaL_checknumber(L, -1); lua_pop(L, 1);
    lua_getfield(L, idx, "vy"); view->plane.v_axis[1] = luaL_checknumber(L, -1); lua_pop(L, 1);
    lua_getfield(L, idx, "vz"); view->plane.v_axis[2] = luaL_checknumber(L, -1); lua_pop(L, 1);

    lua_getfield(L, idx, "u_min"); view->u_min = luaL_checknumber(L, -1); lua_pop(L, 1);
    lua_getfield(L, idx, "u_max"); view->u_max = luaL_checknumber(L, -1); lua_pop(L, 1);
    lua_getfield(L, idx, "v_min"); view->v_min = luaL_checknumber(L, -1); lua_pop(L, 1);
    lua_getfield(L, idx, "v_max"); view->v_max = luaL_checknumber(L, -1); lua_pop(L, 1);
}

/* ============================================================================
 * Module-level functions
 * ============================================================================ */

/* alea.slice_view_axis(axis, val, u0, u1, v0, v1) -> table */
static int l_slice_view_axis(lua_State* L) {
    int axis    = (int)luaL_checkinteger(L, 1);
    double val  = luaL_checknumber(L, 2);
    double u0   = luaL_checknumber(L, 3);
    double u1   = luaL_checknumber(L, 4);
    double v0   = luaL_checknumber(L, 5);
    double v1   = luaL_checknumber(L, 6);

    alea_slice_view_t view;
    alea_slice_view_axis(&view, axis, val, u0, u1, v0, v1);
    push_slice_view(L, &view);
    return 1;
}

/* alea.slice_view(ox,oy,oz,nx,ny,nz,ux,uy,uz,u0,u1,v0,v1) -> table */
static int l_slice_view_init(lua_State* L) {
    double ox = luaL_checknumber(L, 1);
    double oy = luaL_checknumber(L, 2);
    double oz = luaL_checknumber(L, 3);
    double nx = luaL_checknumber(L, 4);
    double ny = luaL_checknumber(L, 5);
    double nz = luaL_checknumber(L, 6);
    double ux = luaL_checknumber(L, 7);
    double uy = luaL_checknumber(L, 8);
    double uz = luaL_checknumber(L, 9);
    double u0 = luaL_checknumber(L, 10);
    double u1 = luaL_checknumber(L, 11);
    double v0 = luaL_checknumber(L, 12);
    double v1 = luaL_checknumber(L, 13);

    alea_slice_view_t view;
    alea_slice_view_init(&view, ox, oy, oz, nx, ny, nz, ux, uy, uz, u0, u1, v0, v1);
    push_slice_view(L, &view);
    return 1;
}

/* alea.slice_curve_set_debug(bool) */
static int l_slice_curve_set_debug(lua_State* L) {
    alea_slice_curve_set_debug(lua_toboolean(L, 1));
    return 0;
}

/* alea.slice_point_trace_set_debug(bool) */
static int l_slice_point_trace_set_debug(lua_State* L) {
    alea_slice_point_trace_set_debug(lua_toboolean(L, 1));
    return 0;
}

/* ============================================================================
 * System methods
 * ============================================================================ */

/* sys:find_cells_grid(view, nu, nv, depth) -> {cell_ids, material_ids, errors} */
static int l_find_cells_grid(lua_State* L) {
    alea_system_t* sys = alea_get_sys(L, 1);
    alea_slice_view_t view;
    lua_to_slice_view(L, 2, &view);
    int nu = (int)luaL_checkinteger(L, 3);
    int nv = (int)luaL_checkinteger(L, 4);
    int depth = (int)luaL_optinteger(L, 5, -1);

    size_t total = (size_t)nu * (size_t)nv;
    int* cell_ids = (int*)calloc(total, sizeof(int));
    int* mat_ids  = (int*)calloc(total, sizeof(int));
    uint8_t* errors = (uint8_t*)calloc(total, sizeof(uint8_t));
    if (!cell_ids || !mat_ids || !errors) {
        free(cell_ids); free(mat_ids); free(errors);
        return luaL_error(L, "find_cells_grid: out of memory");
    }

    int rc = alea_find_cells_grid(sys, &view, nu, nv, depth,
                                   cell_ids, mat_ids, errors);
    if (rc != 0) {
        free(cell_ids); free(mat_ids); free(errors);
        return luaL_error(L, "find_cells_grid failed: %s", alea_error());
    }

    lua_createtable(L, 0, 3);

    /* cell_ids flat table */
    lua_createtable(L, (int)total, 0);
    for (size_t i = 0; i < total; i++) {
        lua_pushinteger(L, cell_ids[i]);
        lua_rawseti(L, -2, (lua_Integer)(i + 1));
    }
    lua_setfield(L, -2, "cell_ids");

    /* material_ids flat table */
    lua_createtable(L, (int)total, 0);
    for (size_t i = 0; i < total; i++) {
        lua_pushinteger(L, mat_ids[i]);
        lua_rawseti(L, -2, (lua_Integer)(i + 1));
    }
    lua_setfield(L, -2, "material_ids");

    /* errors flat table */
    lua_createtable(L, (int)total, 0);
    for (size_t i = 0; i < total; i++) {
        lua_pushinteger(L, errors[i]);
        lua_rawseti(L, -2, (lua_Integer)(i + 1));
    }
    lua_setfield(L, -2, "errors");

    free(cell_ids);
    free(mat_ids);
    free(errors);
    return 1;
}

/* sys:check_grid_overlaps(view, nu, nv, depth, cell_ids, errors) -> errors table */
static int l_check_grid_overlaps(lua_State* L) {
    alea_system_t* sys = alea_get_sys(L, 1);
    alea_slice_view_t view;
    lua_to_slice_view(L, 2, &view);
    int nu = (int)luaL_checkinteger(L, 3);
    int nv = (int)luaL_checkinteger(L, 4);
    int depth = (int)luaL_checkinteger(L, 5);
    luaL_checktype(L, 6, LUA_TTABLE);
    luaL_checktype(L, 7, LUA_TTABLE);

    size_t total = (size_t)nu * (size_t)nv;
    int* cell_ids = (int*)calloc(total, sizeof(int));
    uint8_t* errors = (uint8_t*)calloc(total, sizeof(uint8_t));
    if (!cell_ids || !errors) {
        free(cell_ids); free(errors);
        return luaL_error(L, "check_grid_overlaps: out of memory");
    }

    for (size_t i = 0; i < total; i++) {
        lua_rawgeti(L, 6, (lua_Integer)(i + 1));
        cell_ids[i] = (int)lua_tointeger(L, -1);
        lua_pop(L, 1);
        lua_rawgeti(L, 7, (lua_Integer)(i + 1));
        errors[i] = (uint8_t)lua_tointeger(L, -1);
        lua_pop(L, 1);
    }

    int rc = alea_check_grid_overlaps(sys, &view, nu, nv, depth, cell_ids, errors);
    free(cell_ids);

    if (rc != 0) {
        free(errors);
        return luaL_error(L, "check_grid_overlaps failed: %s", alea_error());
    }

    lua_createtable(L, (int)total, 0);
    for (size_t i = 0; i < total; i++) {
        lua_pushinteger(L, errors[i]);
        lua_rawseti(L, -2, (lua_Integer)(i + 1));
    }
    free(errors);
    return 1;
}

/* sys:get_slice_curves(view) -> Curves */
static int l_get_slice_curves(lua_State* L) {
    alea_system_t* sys = alea_get_sys(L, 1);
    alea_slice_view_t view;
    lua_to_slice_view(L, 2, &view);

    alea_lua_slice_curves_t* ud = (alea_lua_slice_curves_t*)lua_newuserdata(
        L, sizeof(alea_lua_slice_curves_t));
    ud->ptr = NULL;
    luaL_setmetatable(L, ALEA_CURVES_MT);

    ud->ptr = alea_get_slice_curves(sys, &view);
    if (!ud->ptr)
        return luaL_error(L, "get_slice_curves failed: %s", alea_error());
    return 1;
}

/* ============================================================================
 * Curves methods
 * ============================================================================ */

/* curves:count() -> n */
static int l_curves_count(lua_State* L) {
    alea_lua_slice_curves_t* ud = check_curves(L, 1);
    if (!ud->ptr) return luaL_error(L, "curves freed");
    lua_pushinteger(L, (lua_Integer)alea_slice_curves_count(ud->ptr));
    return 1;
}

static const char* curve_type_name(alea_curve_type_t type) {
    switch (type) {
        case ALEA_CURVE_NONE:           return "none";
        case ALEA_CURVE_POINT:          return "point";
        case ALEA_CURVE_LINE:           return "line";
        case ALEA_CURVE_LINE_SEGMENT:   return "line_segment";
        case ALEA_CURVE_RAY:            return "ray";
        case ALEA_CURVE_CIRCLE:         return "circle";
        case ALEA_CURVE_ARC:            return "arc";
        case ALEA_CURVE_ELLIPSE:        return "ellipse";
        case ALEA_CURVE_ELLIPSE_ARC:    return "ellipse_arc";
        case ALEA_CURVE_PARABOLA:       return "parabola";
        case ALEA_CURVE_HYPERBOLA:      return "hyperbola";
        case ALEA_CURVE_POLYGON:        return "polygon";
        case ALEA_CURVE_QUARTIC:        return "quartic";
        case ALEA_CURVE_PARALLEL_LINES: return "parallel_lines";
        default:                        return "unknown";
    }
}

/* curves:get(i) -> table (1-based) */
static int l_curves_get(lua_State* L) {
    alea_lua_slice_curves_t* ud = check_curves(L, 1);
    if (!ud->ptr) return luaL_error(L, "curves freed");
    int i = (int)luaL_checkinteger(L, 2);
    size_t idx = (size_t)(i - 1);

    alea_curve_t curve;
    if (alea_slice_curves_get(ud->ptr, idx, &curve) != 0)
        return luaL_error(L, "curves:get: invalid index %d", i);

    lua_createtable(L, 0, 5);
    lua_pushstring(L, curve_type_name(curve.type)); lua_setfield(L, -2, "type");
    lua_pushinteger(L, curve.surface_id);           lua_setfield(L, -2, "surface_id");
    lua_pushnumber(L, curve.t_min);                 lua_setfield(L, -2, "t_min");
    lua_pushnumber(L, curve.t_max);                 lua_setfield(L, -2, "t_max");

    /* Type-specific data */
    switch (curve.type) {
        case ALEA_CURVE_LINE:
        case ALEA_CURVE_LINE_SEGMENT:
        case ALEA_CURVE_RAY: {
            lua_createtable(L, 0, 4);
            lua_pushnumber(L, curve.data.line.point[0]);     lua_setfield(L, -2, "px");
            lua_pushnumber(L, curve.data.line.point[1]);     lua_setfield(L, -2, "py");
            lua_pushnumber(L, curve.data.line.direction[0]); lua_setfield(L, -2, "dx");
            lua_pushnumber(L, curve.data.line.direction[1]); lua_setfield(L, -2, "dy");
            lua_setfield(L, -2, "data");
            break;
        }
        case ALEA_CURVE_CIRCLE:
        case ALEA_CURVE_ARC: {
            lua_createtable(L, 0, 3);
            lua_pushnumber(L, curve.data.circle.center[0]); lua_setfield(L, -2, "cx");
            lua_pushnumber(L, curve.data.circle.center[1]); lua_setfield(L, -2, "cy");
            lua_pushnumber(L, curve.data.circle.radius);    lua_setfield(L, -2, "radius");
            lua_setfield(L, -2, "data");
            break;
        }
        case ALEA_CURVE_ELLIPSE:
        case ALEA_CURVE_ELLIPSE_ARC: {
            lua_createtable(L, 0, 5);
            lua_pushnumber(L, curve.data.ellipse.center[0]); lua_setfield(L, -2, "cx");
            lua_pushnumber(L, curve.data.ellipse.center[1]); lua_setfield(L, -2, "cy");
            lua_pushnumber(L, curve.data.ellipse.semi_a);    lua_setfield(L, -2, "semi_a");
            lua_pushnumber(L, curve.data.ellipse.semi_b);    lua_setfield(L, -2, "semi_b");
            lua_pushnumber(L, curve.data.ellipse.angle);     lua_setfield(L, -2, "angle");
            lua_setfield(L, -2, "data");
            break;
        }
        case ALEA_CURVE_POLYGON: {
            lua_createtable(L, 0, 3);
            lua_pushinteger(L, curve.data.polygon.count);  lua_setfield(L, -2, "count");
            lua_pushboolean(L, curve.data.polygon.closed); lua_setfield(L, -2, "closed");
            lua_createtable(L, curve.data.polygon.count, 0);
            for (int j = 0; j < curve.data.polygon.count; j++) {
                lua_createtable(L, 2, 0);
                lua_pushnumber(L, curve.data.polygon.vertices[j][0]);
                lua_rawseti(L, -2, 1);
                lua_pushnumber(L, curve.data.polygon.vertices[j][1]);
                lua_rawseti(L, -2, 2);
                lua_rawseti(L, -2, j + 1);
            }
            lua_setfield(L, -2, "vertices");
            lua_setfield(L, -2, "data");
            break;
        }
        default:
            break;
    }
    return 1;
}

/* curves:bounds() -> u0, u1, v0, v1 */
static int l_curves_bounds(lua_State* L) {
    alea_lua_slice_curves_t* ud = check_curves(L, 1);
    if (!ud->ptr) return luaL_error(L, "curves freed");
    double u0, u1, v0, v1;
    alea_slice_curves_bounds(ud->ptr, &u0, &u1, &v0, &v1);
    lua_pushnumber(L, u0);
    lua_pushnumber(L, u1);
    lua_pushnumber(L, v0);
    lua_pushnumber(L, v1);
    return 4;
}

/* __gc for curves */
static int l_curves_gc(lua_State* L) {
    alea_lua_slice_curves_t* ud = check_curves(L, 1);
    if (ud->ptr) {
        alea_slice_curves_free(ud->ptr);
        ud->ptr = NULL;
    }
    return 0;
}

/* ============================================================================
 * Label functions (module-level)
 * ============================================================================ */

/* alea.find_label_positions(ids, w, h, min_px) -> labels table */
static int l_find_label_positions(lua_State* L) {
    luaL_checktype(L, 1, LUA_TTABLE);
    int w = (int)luaL_checkinteger(L, 2);
    int h = (int)luaL_checkinteger(L, 3);
    int min_px = (int)luaL_optinteger(L, 4, 100);

    size_t total = (size_t)w * (size_t)h;
    int* ids = (int*)calloc(total, sizeof(int));
    if (!ids) return luaL_error(L, "out of memory");

    for (size_t i = 0; i < total; i++) {
        lua_rawgeti(L, 1, (lua_Integer)(i + 1));
        ids[i] = (int)lua_tointeger(L, -1);
        lua_pop(L, 1);
    }

    alea_label_position_t* labels = NULL;
    int count = 0;
    int rc = alea_find_label_positions(ids, w, h, min_px, &labels, &count);
    free(ids);
    if (rc != 0)
        return luaL_error(L, "find_label_positions failed");

    lua_createtable(L, count, 0);
    for (int i = 0; i < count; i++) {
        lua_createtable(L, 0, 4);
        lua_pushinteger(L, labels[i].id);          lua_setfield(L, -2, "id");
        lua_pushinteger(L, labels[i].px);          lua_setfield(L, -2, "px");
        lua_pushinteger(L, labels[i].py);          lua_setfield(L, -2, "py");
        lua_pushinteger(L, labels[i].pixel_count); lua_setfield(L, -2, "pixel_count");
        lua_rawseti(L, -2, i + 1);
    }
    free(labels);
    return 1;
}

/* alea.find_surface_label_positions(curves, x_min, x_max, y_min, y_max, w, h, margin[, boundary_ids]) -> labels */
static int l_find_surface_label_positions(lua_State* L) {
    alea_lua_slice_curves_t* ud = check_curves(L, 1);
    if (!ud->ptr) return luaL_error(L, "curves freed");
    double x_min = luaL_checknumber(L, 2);
    double x_max = luaL_checknumber(L, 3);
    double y_min = luaL_checknumber(L, 4);
    double y_max = luaL_checknumber(L, 5);
    int w = (int)luaL_checkinteger(L, 6);
    int h = (int)luaL_checkinteger(L, 7);
    int margin = (int)luaL_optinteger(L, 8, 10);

    int* boundary_ids = NULL;
    if (!lua_isnoneornil(L, 9)) {
        luaL_checktype(L, 9, LUA_TTABLE);
        size_t n = (size_t)w * (size_t)h;
        boundary_ids = malloc(n * sizeof(int));
        if (!boundary_ids) return luaL_error(L, "find_surface_label_positions: out of memory");
        for (size_t i = 0; i < n; i++) {
            lua_rawgeti(L, 9, (lua_Integer)i + 1);
            boundary_ids[i] = (int)luaL_checkinteger(L, -1);
            lua_pop(L, 1);
        }
    }

    alea_label_position_t* labels = NULL;
    int count = 0;
    int rc = boundary_ids
        ? alea_find_surface_label_positions_on_boundaries(
              ud->ptr, boundary_ids, x_min, x_max, y_min, y_max,
              w, h, margin, &labels, &count)
        : alea_find_surface_label_positions(
              ud->ptr, x_min, x_max, y_min, y_max,
              w, h, margin, &labels, &count);
    free(boundary_ids);
    if (rc != 0)
        return luaL_error(L, "find_surface_label_positions failed");

    lua_createtable(L, count, 0);
    for (int i = 0; i < count; i++) {
        lua_createtable(L, 0, 4);
        lua_pushinteger(L, labels[i].id);          lua_setfield(L, -2, "id");
        lua_pushinteger(L, labels[i].px);          lua_setfield(L, -2, "px");
        lua_pushinteger(L, labels[i].py);          lua_setfield(L, -2, "py");
        lua_pushinteger(L, labels[i].pixel_count); lua_setfield(L, -2, "pixel_count");
        lua_rawseti(L, -2, i + 1);
    }
    free(labels);
    return 1;
}

/* ============================================================================
 * Registration
 * ============================================================================ */

static const luaL_Reg slice_system_methods[] = {
    {"find_cells_grid",    l_find_cells_grid},
    {"check_grid_overlaps",l_check_grid_overlaps},
    {"get_slice_curves",   l_get_slice_curves},
    {NULL, NULL}
};

static const luaL_Reg curves_meta[] = {
    {"__gc", l_curves_gc},
    {NULL, NULL}
};

static const luaL_Reg curves_methods[] = {
    {"count",  l_curves_count},
    {"get",    l_curves_get},
    {"bounds", l_curves_bounds},
    {NULL, NULL}
};

int luaopen_alea_slice(lua_State* L) {
    /* Add system methods */
    luaL_getmetatable(L, ALEA_SYSTEM_MT);
    lua_getfield(L, -1, "__index");
    luaL_setfuncs(L, slice_system_methods, 0);
    lua_pop(L, 2);

    /* Create SliceCurves metatable */
    luaL_newmetatable(L, ALEA_CURVES_MT);
    luaL_setfuncs(L, curves_meta, 0);
    lua_newtable(L);
    luaL_setfuncs(L, curves_methods, 0);
    lua_setfield(L, -2, "__index");
    lua_pop(L, 1);

    /* Add module-level functions to alea table (top of stack) */
    lua_pushcfunction(L, l_slice_view_axis);
    lua_setfield(L, -2, "slice_view_axis");
    lua_pushcfunction(L, l_slice_view_init);
    lua_setfield(L, -2, "slice_view");
    lua_pushcfunction(L, l_find_label_positions);
    lua_setfield(L, -2, "find_label_positions");
    lua_pushcfunction(L, l_find_surface_label_positions);
    lua_setfield(L, -2, "find_surface_label_positions");
    lua_pushcfunction(L, l_slice_curve_set_debug);
    lua_setfield(L, -2, "slice_curve_set_debug");
    lua_pushcfunction(L, l_slice_point_trace_set_debug);
    lua_setfield(L, -2, "slice_point_trace_set_debug");

    return 0;
}
