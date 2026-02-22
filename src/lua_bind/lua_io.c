// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

#include "alea_lua.h"

/* ============================================================================
 * Load functions (return new owned System userdata)
 * ============================================================================ */

/* alea.load_mcnp(filename) -> System */
static int l_load_mcnp(lua_State* L) {
    const char* filename = luaL_checkstring(L, 1);
    alea_lua_system_t* ud = (alea_lua_system_t*)lua_newuserdata(L, sizeof(alea_lua_system_t));
    ud->sys = NULL;
    ud->owned = 1;
    luaL_setmetatable(L, ALEA_SYSTEM_MT);

    ud->sys = alea_load_mcnp(filename);
    if (!ud->sys)
        return luaL_error(L, "load_mcnp failed: %s", alea_error());
    return 1;
}

/* alea.load_mcnp_string(str) -> System */
static int l_load_mcnp_string(lua_State* L) {
    size_t len;
    const char* str = luaL_checklstring(L, 1, &len);
    alea_lua_system_t* ud = (alea_lua_system_t*)lua_newuserdata(L, sizeof(alea_lua_system_t));
    ud->sys = NULL;
    ud->owned = 1;
    luaL_setmetatable(L, ALEA_SYSTEM_MT);

    ud->sys = alea_load_mcnp_string(str, len);
    if (!ud->sys)
        return luaL_error(L, "load_mcnp_string failed: %s", alea_error());
    return 1;
}

/* alea.load_openmc(filename) -> System */
static int l_load_openmc(lua_State* L) {
    const char* filename = luaL_checkstring(L, 1);
    alea_lua_system_t* ud = (alea_lua_system_t*)lua_newuserdata(L, sizeof(alea_lua_system_t));
    ud->sys = NULL;
    ud->owned = 1;
    luaL_setmetatable(L, ALEA_SYSTEM_MT);

    ud->sys = alea_load_openmc(filename);
    if (!ud->sys)
        return luaL_error(L, "load_openmc failed: %s", alea_error());
    return 1;
}

/* ============================================================================
 * Export functions (System methods)
 * ============================================================================ */

/* sys:export_mcnp(filename) */
static int l_export_mcnp(lua_State* L) {
    alea_system_t* sys = alea_get_sys(L, 1);
    const char* filename = luaL_checkstring(L, 2);
    if (alea_export_mcnp(sys, filename) != 0)
        return luaL_error(L, "export_mcnp failed: %s", alea_error());
    return 0;
}

/* sys:export_openmc(filename) */
static int l_export_openmc(lua_State* L) {
    alea_system_t* sys = alea_get_sys(L, 1);
    const char* filename = luaL_checkstring(L, 2);
    if (alea_export_openmc(sys, filename) != 0)
        return luaL_error(L, "export_openmc failed: %s", alea_error());
    return 0;
}

/* ============================================================================
 * Registration
 * ============================================================================ */

static const luaL_Reg system_io_methods[] = {
    {"export_mcnp",  l_export_mcnp},
    {"export_openmc", l_export_openmc},
    {NULL, NULL}
};

int luaopen_alea_io(lua_State* L) {
    /* Add export methods to System metatable's __index */
    luaL_getmetatable(L, ALEA_SYSTEM_MT);
    lua_getfield(L, -1, "__index");
    luaL_setfuncs(L, system_io_methods, 0);
    lua_pop(L, 2);

    /* Add load functions to alea table (top of stack) */
    lua_pushcfunction(L, l_load_mcnp);
    lua_setfield(L, -2, "load_mcnp");
    lua_pushcfunction(L, l_load_mcnp_string);
    lua_setfield(L, -2, "load_mcnp_string");
    lua_pushcfunction(L, l_load_openmc);
    lua_setfield(L, -2, "load_openmc");

    return 0;
}
