// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

#include "alea_lua.h"
#include "alea_raycast.h"
#include <stdlib.h>

/* ============================================================================
 * Userdata
 * ============================================================================ */

typedef struct {
    alea_raycast_result_t* ptr;
} alea_lua_raycast_result_t;

static alea_lua_raycast_result_t* check_rayresult(lua_State* L, int idx) {
    return (alea_lua_raycast_result_t*)luaL_checkudata(L, idx, ALEA_RAYRESULT_MT);
}

/* ============================================================================
 * System methods
 * ============================================================================ */

/* sys:raycast(ox, oy, oz, dx, dy, dz, t_max) -> Result */
static int l_raycast(lua_State* L) {
    alea_system_t* sys = alea_get_sys(L, 1);
    double ox = luaL_checknumber(L, 2);
    double oy = luaL_checknumber(L, 3);
    double oz = luaL_checknumber(L, 4);
    double dx = luaL_checknumber(L, 5);
    double dy = luaL_checknumber(L, 6);
    double dz = luaL_checknumber(L, 7);
    double t_max = luaL_optnumber(L, 8, 0.0);

    alea_lua_raycast_result_t* ud = (alea_lua_raycast_result_t*)lua_newuserdata(
        L, sizeof(alea_lua_raycast_result_t));
    ud->ptr = NULL;
    luaL_setmetatable(L, ALEA_RAYRESULT_MT);

    ud->ptr = alea_raycast_result_create();
    if (!ud->ptr)
        return luaL_error(L, "raycast: out of memory");

    if (alea_raycast(sys, ox, oy, oz, dx, dy, dz, t_max, ud->ptr) != 0) {
        return luaL_error(L, "raycast failed: %s", alea_error());
    }
    return 1;
}

/* sys:raycast_cell_aware(ox, oy, oz, dx, dy, dz, t_max) -> Result
 * Semantic equivalent of raycast(); lattice models use the DDA-aware path. */
static int l_raycast_cell_aware(lua_State* L) {
    alea_system_t* sys = alea_get_sys(L, 1);
    double ox = luaL_checknumber(L, 2);
    double oy = luaL_checknumber(L, 3);
    double oz = luaL_checknumber(L, 4);
    double dx = luaL_checknumber(L, 5);
    double dy = luaL_checknumber(L, 6);
    double dz = luaL_checknumber(L, 7);
    double t_max = luaL_optnumber(L, 8, 0.0);

    alea_lua_raycast_result_t* ud = (alea_lua_raycast_result_t*)lua_newuserdata(
        L, sizeof(alea_lua_raycast_result_t));
    ud->ptr = NULL;
    luaL_setmetatable(L, ALEA_RAYRESULT_MT);

    ud->ptr = alea_raycast_result_create();
    if (!ud->ptr)
        return luaL_error(L, "raycast_cell_aware: out of memory");

    if (alea_raycast_cell_aware(sys, ox, oy, oz, dx, dy, dz, t_max, ud->ptr) != 0) {
        return luaL_error(L, "raycast_cell_aware failed: %s", alea_error());
    }
    return 1;
}

/* sys:ray_first_cell(ox, oy, oz, dx, dy, dz, t_max) -> cell_id, t or nil */
static int l_ray_first_cell(lua_State* L) {
    alea_system_t* sys = alea_get_sys(L, 1);
    double ox = luaL_checknumber(L, 2);
    double oy = luaL_checknumber(L, 3);
    double oz = luaL_checknumber(L, 4);
    double dx = luaL_checknumber(L, 5);
    double dy = luaL_checknumber(L, 6);
    double dz = luaL_checknumber(L, 7);
    double t_max = luaL_optnumber(L, 8, 0.0);
    double t;

    int cell_id = alea_ray_first_cell(sys, ox, oy, oz, dx, dy, dz, t_max, &t);
    if (cell_id < 0) {
        lua_pushnil(L);
        return 1;
    }
    lua_pushinteger(L, cell_id);
    lua_pushnumber(L, t);
    return 2;
}

/* ============================================================================
 * Result methods
 * ============================================================================ */

/* result:segment_count() -> n */
static int l_rayresult_segment_count(lua_State* L) {
    alea_lua_raycast_result_t* ud = check_rayresult(L, 1);
    if (!ud->ptr) return luaL_error(L, "raycast result freed");
    lua_pushinteger(L, (lua_Integer)alea_raycast_segment_count(ud->ptr));
    return 1;
}

/* result:segment(i) -> table {t_enter, t_exit, cell_id, material_id, density, enter_surface_id, exit_surface_id} */
static int l_rayresult_segment(lua_State* L) {
    alea_lua_raycast_result_t* ud = check_rayresult(L, 1);
    if (!ud->ptr) return luaL_error(L, "raycast result freed");
    int i = (int)luaL_checkinteger(L, 2);
    /* Lua 1-based to C 0-based */
    size_t idx = (size_t)(i - 1);

    double t_enter, t_exit, density;
    int cell_id, material_id, enter_surface_id, exit_surface_id;
    if (alea_raycast_segment_get(ud->ptr, idx, &t_enter, &t_exit,
                                  &cell_id, &material_id, &density,
                                  &enter_surface_id, &exit_surface_id) != 0)
        return luaL_error(L, "segment: invalid index %d", i);

    lua_createtable(L, 0, 7);
    lua_pushnumber(L, t_enter);     lua_setfield(L, -2, "t_enter");
    lua_pushnumber(L, t_exit);      lua_setfield(L, -2, "t_exit");
    lua_pushinteger(L, cell_id);    lua_setfield(L, -2, "cell_id");
    lua_pushinteger(L, material_id);lua_setfield(L, -2, "material_id");
    lua_pushnumber(L, density);     lua_setfield(L, -2, "density");
    lua_pushinteger(L, enter_surface_id); lua_setfield(L, -2, "enter_surface_id");
    lua_pushinteger(L, exit_surface_id);  lua_setfield(L, -2, "exit_surface_id");
    return 1;
}

/* result:segments() -> table of all segments */
static int l_rayresult_segments(lua_State* L) {
    alea_lua_raycast_result_t* ud = check_rayresult(L, 1);
    if (!ud->ptr) return luaL_error(L, "raycast result freed");
    size_t n = alea_raycast_segment_count(ud->ptr);

    lua_createtable(L, (int)n, 0);
    for (size_t i = 0; i < n; i++) {
        double t_enter, t_exit, density;
        int cell_id, material_id, enter_surface_id, exit_surface_id;
        alea_raycast_segment_get(ud->ptr, i, &t_enter, &t_exit,
                                  &cell_id, &material_id, &density,
                                  &enter_surface_id, &exit_surface_id);

        lua_createtable(L, 0, 7);
        lua_pushnumber(L, t_enter);     lua_setfield(L, -2, "t_enter");
        lua_pushnumber(L, t_exit);      lua_setfield(L, -2, "t_exit");
        lua_pushinteger(L, cell_id);    lua_setfield(L, -2, "cell_id");
        lua_pushinteger(L, material_id);lua_setfield(L, -2, "material_id");
        lua_pushnumber(L, density);     lua_setfield(L, -2, "density");
        lua_pushinteger(L, enter_surface_id); lua_setfield(L, -2, "enter_surface_id");
        lua_pushinteger(L, exit_surface_id);  lua_setfield(L, -2, "exit_surface_id");
        lua_rawseti(L, -2, (lua_Integer)(i + 1));
    }
    return 1;
}

/* result:path_length(mat_id) -> double */
static int l_rayresult_path_length(lua_State* L) {
    alea_lua_raycast_result_t* ud = check_rayresult(L, 1);
    if (!ud->ptr) return luaL_error(L, "raycast result freed");
    int mat_id = (int)luaL_optinteger(L, 2, -1);
    lua_pushnumber(L, alea_raycast_path_length(ud->ptr, mat_id));
    return 1;
}

/* __gc */
static int l_rayresult_gc(lua_State* L) {
    alea_lua_raycast_result_t* ud = check_rayresult(L, 1);
    if (ud->ptr) {
        alea_raycast_result_destroy(ud->ptr);
        ud->ptr = NULL;
    }
    return 0;
}

/* ============================================================================
 * Registration
 * ============================================================================ */

static const luaL_Reg raycast_system_methods[] = {
    {"raycast",            l_raycast},
    {"raycast_cell_aware", l_raycast_cell_aware},
    {"ray_first_cell",     l_ray_first_cell},
    {NULL, NULL}
};

static const luaL_Reg rayresult_meta[] = {
    {"__gc", l_rayresult_gc},
    {NULL, NULL}
};

static const luaL_Reg rayresult_methods[] = {
    {"segment_count", l_rayresult_segment_count},
    {"segment",       l_rayresult_segment},
    {"segments",      l_rayresult_segments},
    {"path_length",   l_rayresult_path_length},
    {NULL, NULL}
};

int luaopen_alea_raycast(lua_State* L) {
    /* Add system methods */
    luaL_getmetatable(L, ALEA_SYSTEM_MT);
    lua_getfield(L, -1, "__index");
    luaL_setfuncs(L, raycast_system_methods, 0);
    lua_pop(L, 2);

    /* Create RaycastResult metatable */
    luaL_newmetatable(L, ALEA_RAYRESULT_MT);
    luaL_setfuncs(L, rayresult_meta, 0);
    lua_newtable(L);
    luaL_setfuncs(L, rayresult_methods, 0);
    lua_setfield(L, -2, "__index");
    lua_pop(L, 1);

    return 0;
}
