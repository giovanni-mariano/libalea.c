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

/* sys:first_visible(ox, oy, oz, dx, dy, dz [, t_max]) -> table or nil */
static int l_first_visible(lua_State* L) {
    alea_system_t* sys = alea_get_sys(L, 1);
    alea_ray_first_visible_options_t options;
    alea_ray_first_visible_options_init(&options);
    options.fields = ALEA_RAY_FIRST_VISIBLE_SURFACE_ID |
                     ALEA_RAY_FIRST_VISIBLE_SURFACE_NORMAL;
    if (lua_istable(L, 8)) {
        lua_getfield(L, 8, "t_min"); if (!lua_isnil(L, -1)) options.t_min = lua_tonumber(L, -1); lua_pop(L, 1);
        lua_getfield(L, 8, "t_max"); if (!lua_isnil(L, -1)) options.t_max = lua_tonumber(L, -1); lua_pop(L, 1);
        lua_getfield(L, 8, "material_filter"); if (!lua_isnil(L, -1)) options.material_filter = (int)lua_tointeger(L, -1); lua_pop(L, 1);
        lua_getfield(L, 8, "surface_id"); if (lua_isboolean(L, -1) && !lua_toboolean(L, -1)) options.fields &= ~ALEA_RAY_FIRST_VISIBLE_SURFACE_ID; lua_pop(L, 1);
        lua_getfield(L, 8, "normal"); if (lua_isboolean(L, -1) && !lua_toboolean(L, -1)) options.fields &= ~ALEA_RAY_FIRST_VISIBLE_SURFACE_NORMAL; lua_pop(L, 1);
    } else options.t_max = luaL_optnumber(L, 8, 0.0);
    alea_ray_first_visible_query_result_t* result =
        alea_ray_first_visible_query_result_create();
    if (!result) return luaL_error(L, "first_visible: out of memory");
    int rc = alea_ray_first_visible_query(sys,
        luaL_checknumber(L, 2), luaL_checknumber(L, 3), luaL_checknumber(L, 4),
        luaL_checknumber(L, 5), luaL_checknumber(L, 6), luaL_checknumber(L, 7),
        &options, result);
    if (rc != 0) { alea_ray_first_visible_query_result_destroy(result); return luaL_error(L, "first_visible failed: %s", alea_error()); }
    if (!alea_ray_first_visible_found(result)) { alea_ray_first_visible_query_result_destroy(result); lua_pushnil(L); return 1; }
    double nx, ny, nz;
    lua_createtable(L, 0, 7);
    lua_pushnumber(L, alea_ray_first_visible_t(result)); lua_setfield(L, -2, "t");
    lua_pushinteger(L, alea_ray_first_visible_cell_id(result)); lua_setfield(L, -2, "cell_id");
    lua_pushinteger(L, alea_ray_first_visible_material_id(result)); lua_setfield(L, -2, "material_id");
    lua_pushnumber(L, alea_ray_first_visible_density(result)); lua_setfield(L, -2, "density");
    lua_pushinteger(L, alea_ray_first_visible_surface_id(result)); lua_setfield(L, -2, "surface_id");
    if (alea_ray_first_visible_normal(result, &nx, &ny, &nz) == 0) { lua_createtable(L, 3, 0); lua_pushnumber(L,nx); lua_rawseti(L,-2,1); lua_pushnumber(L,ny); lua_rawseti(L,-2,2); lua_pushnumber(L,nz); lua_rawseti(L,-2,3); lua_setfield(L,-2,"normal"); }
    alea_ray_first_visible_query_result_destroy(result);
    return 1;
}

/* sys:boundary_events(ox, oy, oz, dx, dy, dz [, t_max]) -> {events...} */
static int l_boundary_events(lua_State* L) {
    alea_system_t* sys = alea_get_sys(L, 1);
    alea_ray_boundary_event_options_t options;
    alea_ray_boundary_event_options_init(&options);
    options.fields = ALEA_RAY_BOUNDARY_EVENT_PRIMITIVE_ID |
                     ALEA_RAY_BOUNDARY_EVENT_NORMAL;
    if (lua_istable(L, 8)) {
        lua_getfield(L, 8, "t_min"); if (!lua_isnil(L, -1)) options.t_min = lua_tonumber(L, -1); lua_pop(L, 1);
        lua_getfield(L, 8, "t_max"); if (!lua_isnil(L, -1)) options.t_max = lua_tonumber(L, -1); lua_pop(L, 1);
        lua_getfield(L, 8, "include_all_coincident_physical"); options.include_all_coincident_physical = lua_toboolean(L, -1); lua_pop(L, 1);
        lua_getfield(L, 8, "max_events"); if (!lua_isnil(L, -1)) options.max_events = (uint64_t)lua_tointeger(L, -1); lua_pop(L, 1);
        lua_getfield(L, 8, "max_output_bytes"); if (!lua_isnil(L, -1)) options.max_output_bytes = (uint64_t)lua_tointeger(L, -1); lua_pop(L, 1);
        lua_getfield(L, 8, "primitive_id"); if (lua_isboolean(L, -1) && !lua_toboolean(L, -1)) options.fields &= ~ALEA_RAY_BOUNDARY_EVENT_PRIMITIVE_ID; lua_pop(L, 1);
        lua_getfield(L, 8, "normal"); if (lua_isboolean(L, -1) && !lua_toboolean(L, -1)) options.fields &= ~ALEA_RAY_BOUNDARY_EVENT_NORMAL; lua_pop(L, 1);
    } else options.t_max = luaL_optnumber(L, 8, 0.0);
    alea_ray_boundary_event_query_result_t* result =
        alea_ray_boundary_event_query_result_create();
    if (!result) return luaL_error(L, "boundary_events: out of memory");
    int rc = alea_ray_boundary_event_query(sys,
        luaL_checknumber(L, 2), luaL_checknumber(L, 3), luaL_checknumber(L, 4),
        luaL_checknumber(L, 5), luaL_checknumber(L, 6), luaL_checknumber(L, 7),
        &options, result);
    if (rc != 0) { alea_ray_boundary_event_query_result_destroy(result); return luaL_error(L, "boundary_events failed: %s", alea_error()); }
    size_t count = alea_ray_boundary_event_count(result);
    lua_createtable(L, (int)count, 0);
    for (size_t i = 0; i < count; i++) {
        double t, nx, ny, nz; int kind, surface_id, before, after, mat_before, mat_after;
        uint32_t flags, primitive_id;
        alea_ray_boundary_event_get(result, i, &t, &kind, &surface_id, &before, &after,
                                    &mat_before, &mat_after, &flags, &primitive_id,
                                    &nx, &ny, &nz);
        lua_createtable(L, 0, 11);
        lua_pushnumber(L, t); lua_setfield(L, -2, "t");
        lua_pushinteger(L, kind); lua_setfield(L, -2, "kind");
        lua_pushinteger(L, surface_id); lua_setfield(L, -2, "surface_id");
        lua_pushinteger(L, before); lua_setfield(L, -2, "cell_before");
        lua_pushinteger(L, after); lua_setfield(L, -2, "cell_after");
        lua_pushinteger(L, mat_before); lua_setfield(L, -2, "material_before");
        lua_pushinteger(L, mat_after); lua_setfield(L, -2, "material_after");
        lua_pushinteger(L, primitive_id); lua_setfield(L, -2, "primitive_id");
        lua_pushinteger(L, flags); lua_setfield(L, -2, "resolution_flags");
        lua_createtable(L, 3, 0); lua_pushnumber(L,nx); lua_rawseti(L,-2,1); lua_pushnumber(L,ny); lua_rawseti(L,-2,2); lua_pushnumber(L,nz); lua_rawseti(L,-2,3); lua_setfield(L,-2,"normal");
        lua_rawseti(L, -2, (lua_Integer)(i + 1));
    }
    alea_ray_boundary_event_query_result_destroy(result);
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
    {"first_visible",      l_first_visible},
    {"boundary_events",    l_boundary_events},
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
