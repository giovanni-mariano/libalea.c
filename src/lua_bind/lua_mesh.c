// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

#include "alea_lua.h"
#include "alea_mesh.h"
#include <stdlib.h>

/* ============================================================================
 * Userdata
 * ============================================================================ */

typedef struct {
    alea_mesh_result_t* ptr;
} alea_lua_mesh_result_t;

static alea_lua_mesh_result_t* check_meshresult(lua_State* L, int idx) {
    return (alea_lua_mesh_result_t*)luaL_checkudata(L, idx, ALEA_MESHRESULT_MT);
}

/* ============================================================================
 * Config helper: Lua table -> alea_mesh_config_t
 * ============================================================================ */

static void lua_to_mesh_config(lua_State* L, int idx, alea_mesh_config_t* cfg) {
    luaL_checktype(L, idx, LUA_TTABLE);

    lua_getfield(L, idx, "x_min");
    if (!lua_isnil(L, -1)) cfg->x_min = lua_tonumber(L, -1);
    lua_pop(L, 1);
    lua_getfield(L, idx, "x_max");
    if (!lua_isnil(L, -1)) cfg->x_max = lua_tonumber(L, -1);
    lua_pop(L, 1);
    lua_getfield(L, idx, "y_min");
    if (!lua_isnil(L, -1)) cfg->y_min = lua_tonumber(L, -1);
    lua_pop(L, 1);
    lua_getfield(L, idx, "y_max");
    if (!lua_isnil(L, -1)) cfg->y_max = lua_tonumber(L, -1);
    lua_pop(L, 1);
    lua_getfield(L, idx, "z_min");
    if (!lua_isnil(L, -1)) cfg->z_min = lua_tonumber(L, -1);
    lua_pop(L, 1);
    lua_getfield(L, idx, "z_max");
    if (!lua_isnil(L, -1)) cfg->z_max = lua_tonumber(L, -1);
    lua_pop(L, 1);

    lua_getfield(L, idx, "nx");
    if (!lua_isnil(L, -1)) cfg->nx = (int)lua_tointeger(L, -1);
    lua_pop(L, 1);
    lua_getfield(L, idx, "ny");
    if (!lua_isnil(L, -1)) cfg->ny = (int)lua_tointeger(L, -1);
    lua_pop(L, 1);
    lua_getfield(L, idx, "nz");
    if (!lua_isnil(L, -1)) cfg->nz = (int)lua_tointeger(L, -1);
    lua_pop(L, 1);

    lua_getfield(L, idx, "format");
    if (!lua_isnil(L, -1)) cfg->format = (alea_mesh_format_t)(int)lua_tointeger(L, -1);
    lua_pop(L, 1);

    lua_getfield(L, idx, "void_material_id");
    if (!lua_isnil(L, -1)) cfg->void_material_id = (int)lua_tointeger(L, -1);
    lua_pop(L, 1);

    lua_getfield(L, idx, "auto_pad");
    if (!lua_isnil(L, -1)) cfg->auto_pad = lua_tonumber(L, -1);
    lua_pop(L, 1);
}

/* ============================================================================
 * System methods
 * ============================================================================ */

/* sys:mesh_sample(config) -> MeshResult */
static int l_mesh_sample(lua_State* L) {
    alea_system_t* sys = alea_get_sys(L, 1);
    alea_mesh_config_t cfg;
    alea_mesh_config_init(&cfg);
    if (lua_istable(L, 2))
        lua_to_mesh_config(L, 2, &cfg);

    alea_lua_mesh_result_t* ud = (alea_lua_mesh_result_t*)lua_newuserdata(
        L, sizeof(alea_lua_mesh_result_t));
    ud->ptr = NULL;
    luaL_setmetatable(L, ALEA_MESHRESULT_MT);

    ud->ptr = alea_mesh_sample(sys, &cfg);
    if (!ud->ptr)
        return luaL_error(L, "mesh_sample failed: %s", alea_error());
    return 1;
}

/* sys:mesh_export(config, filename) */
static int l_mesh_export_system(lua_State* L) {
    alea_system_t* sys = alea_get_sys(L, 1);
    alea_mesh_config_t cfg;
    alea_mesh_config_init(&cfg);
    if (lua_istable(L, 2))
        lua_to_mesh_config(L, 2, &cfg);
    const char* filename = luaL_checkstring(L, 3);

    if (alea_mesh_export_system(sys, &cfg, filename) != 0)
        return luaL_error(L, "mesh_export failed: %s", alea_error());
    return 0;
}

/* ============================================================================
 * MeshResult methods
 * ============================================================================ */

/* mesh:export(format, filename) */
static int l_meshresult_export(lua_State* L) {
    alea_lua_mesh_result_t* ud = check_meshresult(L, 1);
    if (!ud->ptr) return luaL_error(L, "mesh result freed");
    int fmt = (int)luaL_checkinteger(L, 2);
    const char* filename = luaL_checkstring(L, 3);
    if (alea_mesh_export(ud->ptr, (alea_mesh_format_t)fmt, filename) != 0)
        return luaL_error(L, "mesh:export failed: %s", alea_error());
    return 0;
}

/* mesh:info() -> table */
static int l_meshresult_info(lua_State* L) {
    alea_lua_mesh_result_t* ud = check_meshresult(L, 1);
    if (!ud->ptr) return luaL_error(L, "mesh result freed");
    alea_mesh_result_t* m = ud->ptr;

    lua_createtable(L, 0, 5);
    lua_pushinteger(L, m->nx);            lua_setfield(L, -2, "nx");
    lua_pushinteger(L, m->ny);            lua_setfield(L, -2, "ny");
    lua_pushinteger(L, m->nz);            lua_setfield(L, -2, "nz");
    lua_pushinteger(L, m->num_materials); lua_setfield(L, -2, "num_materials");

    /* unique_materials */
    lua_createtable(L, m->num_materials, 0);
    for (int i = 0; i < m->num_materials; i++) {
        lua_pushinteger(L, m->unique_materials[i]);
        lua_rawseti(L, -2, i + 1);
    }
    lua_setfield(L, -2, "unique_materials");
    return 1;
}

/* mesh:material_ids() -> flat table */
static int l_meshresult_material_ids(lua_State* L) {
    alea_lua_mesh_result_t* ud = check_meshresult(L, 1);
    if (!ud->ptr) return luaL_error(L, "mesh result freed");
    alea_mesh_result_t* m = ud->ptr;
    size_t total = (size_t)m->nx * (size_t)m->ny * (size_t)m->nz;
    lua_createtable(L, (int)total, 0);
    for (size_t i = 0; i < total; i++) {
        lua_pushinteger(L, m->material_ids[i]);
        lua_rawseti(L, -2, (lua_Integer)(i + 1));
    }
    return 1;
}

/* mesh:cell_ids() -> flat table */
static int l_meshresult_cell_ids(lua_State* L) {
    alea_lua_mesh_result_t* ud = check_meshresult(L, 1);
    if (!ud->ptr) return luaL_error(L, "mesh result freed");
    alea_mesh_result_t* m = ud->ptr;
    size_t total = (size_t)m->nx * (size_t)m->ny * (size_t)m->nz;
    lua_createtable(L, (int)total, 0);
    for (size_t i = 0; i < total; i++) {
        lua_pushinteger(L, m->cell_ids[i]);
        lua_rawseti(L, -2, (lua_Integer)(i + 1));
    }
    return 1;
}

/* __gc */
static int l_meshresult_gc(lua_State* L) {
    alea_lua_mesh_result_t* ud = check_meshresult(L, 1);
    if (ud->ptr) {
        alea_mesh_result_free(ud->ptr);
        ud->ptr = NULL;
    }
    return 0;
}

/* ============================================================================
 * Registration
 * ============================================================================ */

static const luaL_Reg mesh_system_methods[] = {
    {"mesh_sample",  l_mesh_sample},
    {"mesh_export",  l_mesh_export_system},
    {NULL, NULL}
};

static const luaL_Reg meshresult_meta[] = {
    {"__gc", l_meshresult_gc},
    {NULL, NULL}
};

static const luaL_Reg meshresult_methods[] = {
    {"export",       l_meshresult_export},
    {"info",         l_meshresult_info},
    {"material_ids", l_meshresult_material_ids},
    {"cell_ids",     l_meshresult_cell_ids},
    {NULL, NULL}
};

int luaopen_alea_mesh(lua_State* L) {
    /* Add system methods */
    luaL_getmetatable(L, ALEA_SYSTEM_MT);
    lua_getfield(L, -1, "__index");
    luaL_setfuncs(L, mesh_system_methods, 0);
    lua_pop(L, 2);

    /* Create MeshResult metatable */
    luaL_newmetatable(L, ALEA_MESHRESULT_MT);
    luaL_setfuncs(L, meshresult_meta, 0);
    lua_newtable(L);
    luaL_setfuncs(L, meshresult_methods, 0);
    lua_setfield(L, -2, "__index");
    lua_pop(L, 1);

    return 0;
}
