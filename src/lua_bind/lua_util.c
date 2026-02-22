// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

#include "alea_lua.h"
#include <stdlib.h>

/* ============================================================================
 * Universe operations
 * ============================================================================ */

/* sys:flatten(universe_id) */
static int l_flatten(lua_State* L) {
    alea_system_t* sys = alea_get_sys(L, 1);
    int uid = (int)luaL_optinteger(L, 2, 0);
    if (alea_flatten(sys, uid) < 0)
        return luaL_error(L, "flatten failed: %s", alea_error());
    return 0;
}

/* sys:split_union_cells() -> number of new cells */
static int l_split_union_cells(lua_State* L) {
    alea_system_t* sys = alea_get_sys(L, 1);
    int n = alea_split_union_cells(sys);
    if (n < 0)
        return luaL_error(L, "split_union_cells failed: %s", alea_error());
    lua_pushinteger(L, n);
    return 1;
}

/* sys:extract_universe(id) -> new System */
static int l_extract_universe(lua_State* L) {
    alea_system_t* sys = alea_get_sys(L, 1);
    int uid = (int)luaL_checkinteger(L, 2);
    alea_lua_system_t* ud = (alea_lua_system_t*)lua_newuserdata(L, sizeof(alea_lua_system_t));
    ud->sys = NULL;
    ud->owned = 1;
    luaL_setmetatable(L, ALEA_SYSTEM_MT);
    ud->sys = alea_extract_universe(sys, uid);
    if (!ud->sys)
        return luaL_error(L, "extract_universe failed: %s", alea_error());
    return 1;
}

/* sys:merge(other, offset) */
static int l_merge(lua_State* L) {
    alea_system_t* sys = alea_get_sys(L, 1);
    alea_system_t* other = alea_get_sys(L, 2);
    int offset = (int)luaL_optinteger(L, 3, 0);
    if (alea_merge(sys, other, offset) < 0)
        return luaL_error(L, "merge failed: %s", alea_error());
    return 0;
}

/* sys:clone() -> new System */
static int l_clone(lua_State* L) {
    alea_system_t* sys = alea_get_sys(L, 1);
    alea_lua_system_t* ud = (alea_lua_system_t*)lua_newuserdata(L, sizeof(alea_lua_system_t));
    ud->sys = NULL;
    ud->owned = 1;
    luaL_setmetatable(L, ALEA_SYSTEM_MT);
    ud->sys = alea_clone(sys);
    if (!ud->sys)
        return luaL_error(L, "clone failed: %s", alea_error());
    return 1;
}

/* ============================================================================
 * Volume estimation
 * ============================================================================ */

/* sys:bounding_sphere(tol) -> cx, cy, cz, radius */
static int l_bounding_sphere(lua_State* L) {
    alea_system_t* sys = alea_get_sys(L, 1);
    double tol = luaL_optnumber(L, 2, 1.0);
    double cx, cy, cz, r;
    if (alea_compute_bounding_sphere(sys, tol, &cx, &cy, &cz, &r) != 0)
        return luaL_error(L, "bounding_sphere failed: %s", alea_error());
    lua_pushnumber(L, cx);
    lua_pushnumber(L, cy);
    lua_pushnumber(L, cz);
    lua_pushnumber(L, r);
    return 4;
}

/* sys:estimate_volumes(n_rays, cx, cy, cz, radius) -> table of volumes */
static int l_estimate_volumes(lua_State* L) {
    alea_system_t* sys = alea_get_sys(L, 1);
    int n_rays = (int)luaL_checkinteger(L, 2);
    double cx  = luaL_checknumber(L, 3);
    double cy  = luaL_checknumber(L, 4);
    double cz  = luaL_checknumber(L, 5);
    double r   = luaL_checknumber(L, 6);

    size_t nc = alea_cell_count(sys);
    double* volumes = (double*)calloc(nc, sizeof(double));
    double* errors  = (double*)calloc(nc, sizeof(double));
    if (!volumes || !errors) {
        free(volumes);
        free(errors);
        return luaL_error(L, "out of memory");
    }

    int rc = alea_estimate_cell_volumes(sys, cx, cy, cz, r, n_rays, volumes, errors);
    if (rc != 0) {
        free(volumes);
        free(errors);
        return luaL_error(L, "estimate_volumes failed: %s", alea_error());
    }

    lua_createtable(L, (int)nc, 0);
    for (size_t i = 0; i < nc; i++) {
        lua_createtable(L, 0, 2);
        lua_pushnumber(L, volumes[i]); lua_setfield(L, -2, "volume");
        lua_pushnumber(L, errors[i]);  lua_setfield(L, -2, "rel_error");
        lua_rawseti(L, -2, (lua_Integer)(i + 1));
    }

    free(volumes);
    free(errors);
    return 1;
}

/* sys:remove_cells_by_volume(volumes_table, threshold) -> count removed */
static int l_remove_cells_by_volume(lua_State* L) {
    alea_system_t* sys = alea_get_sys(L, 1);
    luaL_checktype(L, 2, LUA_TTABLE);
    double threshold = luaL_checknumber(L, 3);

    size_t nc = alea_cell_count(sys);
    double* volumes = (double*)calloc(nc, sizeof(double));
    if (!volumes)
        return luaL_error(L, "out of memory");

    for (size_t i = 0; i < nc; i++) {
        lua_rawgeti(L, 2, (lua_Integer)(i + 1));
        if (lua_istable(L, -1)) {
            lua_getfield(L, -1, "volume");
            volumes[i] = lua_tonumber(L, -1);
            lua_pop(L, 1);
        } else {
            volumes[i] = lua_tonumber(L, -1);
        }
        lua_pop(L, 1);
    }

    int removed = alea_remove_cells_by_volume(sys, volumes, threshold);
    free(volumes);
    if (removed < 0)
        return luaL_error(L, "remove_cells_by_volume failed: %s", alea_error());
    lua_pushinteger(L, removed);
    return 1;
}

/* ============================================================================
 * Bbox tightening
 * ============================================================================ */

/* sys:tighten_all_bboxes(tol) -> count */
static int l_tighten_all_bboxes(lua_State* L) {
    alea_system_t* sys = alea_get_sys(L, 1);
    double tol = luaL_optnumber(L, 2, 1.0);
    int n = alea_tighten_all_bboxes(sys, tol);
    if (n < 0)
        return luaL_error(L, "tighten_all_bboxes failed: %s", alea_error());
    lua_pushinteger(L, n);
    return 1;
}

/* sys:tighten_bbox_numerical(cell_index) */
static int l_tighten_bbox_numerical(lua_State* L) {
    alea_system_t* sys = alea_get_sys(L, 1);
    int idx = (int)luaL_checkinteger(L, 2);
    if (alea_tighten_cell_bbox_numerical(sys, idx) != 0)
        return luaL_error(L, "tighten_bbox_numerical failed: %s", alea_error());
    return 0;
}

/* ============================================================================
 * Renumbering
 * ============================================================================ */

static int l_renumber_cells(lua_State* L) {
    alea_system_t* sys = alea_get_sys(L, 1);
    int start = (int)luaL_checkinteger(L, 2);
    if (alea_renumber_cells(sys, start) != 0)
        return luaL_error(L, "renumber_cells failed: %s", alea_error());
    return 0;
}

static int l_renumber_surfaces(lua_State* L) {
    alea_system_t* sys = alea_get_sys(L, 1);
    int start = (int)luaL_checkinteger(L, 2);
    if (alea_renumber_surfaces(sys, start) != 0)
        return luaL_error(L, "renumber_surfaces failed: %s", alea_error());
    return 0;
}

static int l_offset_cell_ids(lua_State* L) {
    alea_system_t* sys = alea_get_sys(L, 1);
    int offset = (int)luaL_checkinteger(L, 2);
    if (alea_offset_cell_ids(sys, offset) != 0)
        return luaL_error(L, "offset_cell_ids failed: %s", alea_error());
    return 0;
}

static int l_offset_surface_ids(lua_State* L) {
    alea_system_t* sys = alea_get_sys(L, 1);
    int offset = (int)luaL_checkinteger(L, 2);
    if (alea_offset_surface_ids(sys, offset) != 0)
        return luaL_error(L, "offset_surface_ids failed: %s", alea_error());
    return 0;
}

static int l_offset_material_ids(lua_State* L) {
    alea_system_t* sys = alea_get_sys(L, 1);
    int offset = (int)luaL_checkinteger(L, 2);
    if (alea_offset_material_ids(sys, offset) != 0)
        return luaL_error(L, "offset_material_ids failed: %s", alea_error());
    return 0;
}

/* ============================================================================
 * Macrobody expansion
 * ============================================================================ */

static int l_expand_macrobodies(lua_State* L) {
    alea_system_t* sys = alea_get_sys(L, 1);
    int n = alea_expand_macrobodies_in_system(sys);
    if (n < 0)
        return luaL_error(L, "expand_macrobodies failed: %s", alea_error());
    lua_pushinteger(L, n);
    return 1;
}

/* ============================================================================
 * Void generation
 * ============================================================================ */

typedef struct {
    void_result_t* ptr;
    alea_system_t* sys;
} alea_lua_void_result_t;

static int l_void_generate(lua_State* L) {
    alea_system_t* sys = alea_get_sys(L, 1);
    luaL_checktype(L, 2, LUA_TTABLE);

    alea_bbox_t bbox;
    lua_getfield(L, 2, "min_x"); bbox.min_x = luaL_checknumber(L, -1); lua_pop(L, 1);
    lua_getfield(L, 2, "max_x"); bbox.max_x = luaL_checknumber(L, -1); lua_pop(L, 1);
    lua_getfield(L, 2, "min_y"); bbox.min_y = luaL_checknumber(L, -1); lua_pop(L, 1);
    lua_getfield(L, 2, "max_y"); bbox.max_y = luaL_checknumber(L, -1); lua_pop(L, 1);
    lua_getfield(L, 2, "min_z"); bbox.min_z = luaL_checknumber(L, -1); lua_pop(L, 1);
    lua_getfield(L, 2, "max_z"); bbox.max_z = luaL_checknumber(L, -1); lua_pop(L, 1);

    alea_lua_void_result_t* ud = (alea_lua_void_result_t*)lua_newuserdata(
        L, sizeof(alea_lua_void_result_t));
    ud->ptr = NULL;
    ud->sys = sys;
    luaL_setmetatable(L, ALEA_VOIDRESULT_MT);

    ud->ptr = alea_void_generate(sys, &bbox);
    if (!ud->ptr)
        return luaL_error(L, "void_generate failed: %s", alea_error());
    return 1;
}

static int l_void_gc(lua_State* L) {
    alea_lua_void_result_t* ud = (alea_lua_void_result_t*)luaL_checkudata(L, 1, ALEA_VOIDRESULT_MT);
    if (ud->ptr) {
        alea_void_free(ud->ptr);
        ud->ptr = NULL;
    }
    return 0;
}

static int l_void_count(lua_State* L) {
    alea_lua_void_result_t* ud = (alea_lua_void_result_t*)luaL_checkudata(L, 1, ALEA_VOIDRESULT_MT);
    if (!ud->ptr) return luaL_error(L, "void result freed");
    lua_pushinteger(L, (lua_Integer)alea_void_count(ud->ptr));
    return 1;
}

static int l_void_add_cells(lua_State* L) {
    alea_lua_void_result_t* ud = (alea_lua_void_result_t*)luaL_checkudata(L, 1, ALEA_VOIDRESULT_MT);
    if (!ud->ptr) return luaL_error(L, "void result freed");
    int n = alea_void_add_cells(ud->sys, ud->ptr);
    if (n < 0)
        return luaL_error(L, "void_add_cells failed: %s", alea_error());
    lua_pushinteger(L, n);
    return 1;
}

static int l_void_merge(lua_State* L) {
    alea_lua_void_result_t* ud = (alea_lua_void_result_t*)luaL_checkudata(L, 1, ALEA_VOIDRESULT_MT);
    if (!ud->ptr) return luaL_error(L, "void result freed");
    if (alea_void_merge(ud->sys, ud->ptr) < 0)
        return luaL_error(L, "void_merge failed: %s", alea_error());
    return 0;
}

/* ============================================================================
 * Registration
 * ============================================================================ */

static const luaL_Reg util_methods[] = {
    {"flatten",              l_flatten},
    {"split_union_cells",    l_split_union_cells},
    {"extract_universe",     l_extract_universe},
    {"merge",                l_merge},
    {"clone",                l_clone},
    {"bounding_sphere",      l_bounding_sphere},
    {"estimate_volumes",     l_estimate_volumes},
    {"remove_cells_by_volume", l_remove_cells_by_volume},
    {"tighten_all_bboxes",   l_tighten_all_bboxes},
    {"tighten_bbox_numerical", l_tighten_bbox_numerical},
    {"renumber_cells",       l_renumber_cells},
    {"renumber_surfaces",    l_renumber_surfaces},
    {"offset_cell_ids",      l_offset_cell_ids},
    {"offset_surface_ids",   l_offset_surface_ids},
    {"offset_material_ids",  l_offset_material_ids},
    {"expand_macrobodies",   l_expand_macrobodies},
    {"void_generate",        l_void_generate},
    {NULL, NULL}
};

static const luaL_Reg void_meta[] = {
    {"__gc", l_void_gc},
    {NULL, NULL}
};

static const luaL_Reg void_methods[] = {
    {"count",     l_void_count},
    {"add_cells", l_void_add_cells},
    {"merge",     l_void_merge},
    {NULL, NULL}
};

int luaopen_alea_util(lua_State* L) {
    /* Add util methods to System metatable's __index */
    luaL_getmetatable(L, ALEA_SYSTEM_MT);
    lua_getfield(L, -1, "__index");
    luaL_setfuncs(L, util_methods, 0);
    lua_pop(L, 2);

    /* Create VoidResult metatable */
    luaL_newmetatable(L, ALEA_VOIDRESULT_MT);
    luaL_setfuncs(L, void_meta, 0);
    lua_newtable(L);
    luaL_setfuncs(L, void_methods, 0);
    lua_setfield(L, -2, "__index");
    lua_pop(L, 1);

    return 0;
}
