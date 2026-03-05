// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

#include "alea_lua.h"
#include "alea_mcnp.h"
#include "alea_openmc.h"
#include <string.h>

/* ============================================================================
 * System lifecycle
 * ============================================================================ */

static int l_system_create(lua_State* L) {
    alea_lua_system_t* ud = (alea_lua_system_t*)lua_newuserdata(L, sizeof(alea_lua_system_t));
    ud->sys = alea_create();
    ud->mcnp_model = NULL;
    ud->openmc_model = NULL;
    ud->owned = 1;
    if (!ud->sys)
        return luaL_error(L, "alea_create failed: %s", alea_error());
    luaL_setmetatable(L, ALEA_SYSTEM_MT);
    return 1;
}

static int l_system_destroy(lua_State* L) {
    alea_lua_system_t* ud = alea_check_system(L, 1);
    if (ud->owned) {
        if (ud->mcnp_model) {
            mcnp_model_destroy((mcnp_model_t*)ud->mcnp_model);
            ud->mcnp_model = NULL;
        } else if (ud->openmc_model) {
            openmc_model_destroy((openmc_model_t*)ud->openmc_model);
            ud->openmc_model = NULL;
        } else if (ud->sys) {
            alea_destroy(ud->sys);
        }
    }
    ud->sys = NULL;
    return 0;
}

static int l_system_gc(lua_State* L) {
    return l_system_destroy(L);
}

static int l_system_tostring(lua_State* L) {
    alea_lua_system_t* ud = alea_check_system(L, 1);
    if (!ud->sys) {
        lua_pushliteral(L, "System(destroyed)");
    } else {
        lua_pushfstring(L, "System(%d cells, %d surfaces)",
                        (int)alea_cell_count(ud->sys),
                        (int)alea_surface_count(ud->sys));
    }
    return 1;
}

static int l_system_reset(lua_State* L) {
    alea_system_t* sys = alea_get_sys(L, 1);
    alea_reset(sys);
    return 0;
}

static int l_version(lua_State* L) {
    lua_pushstring(L, alea_version());
    return 1;
}

/* ============================================================================
 * Logging and error (module-level functions)
 * ============================================================================ */

static int l_log_set_level(lua_State* L) {
    int level = (int)luaL_checkinteger(L, 1);
    alea_log_set_level((alea_log_level_t)level);
    return 0;
}

static int l_log_get_level(lua_State* L) {
    lua_pushinteger(L, (lua_Integer)alea_log_get_level());
    return 1;
}

static int l_error(lua_State* L) {
    lua_pushstring(L, alea_error());
    return 1;
}

static int l_error_code(lua_State* L) {
    lua_pushinteger(L, (lua_Integer)alea_error_code());
    return 1;
}

static int l_error_clear(lua_State* L) {
    (void)L;
    alea_error_clear();
    return 0;
}

static int l_set_debug_trace(lua_State* L) {
    int enable = lua_toboolean(L, 1);
    alea_set_debug_trace(enable);
    return 0;
}

/* ============================================================================
 * Configuration
 * ============================================================================ */

static int l_system_get_config(lua_State* L) {
    alea_system_t* sys = alea_get_sys(L, 1);
    alea_config_t cfg = alea_get_config(sys);

    lua_newtable(L);

    lua_pushnumber(L, cfg.abs_tol);        lua_setfield(L, -2, "abs_tol");
    lua_pushnumber(L, cfg.rel_tol);        lua_setfield(L, -2, "rel_tol");
    lua_pushnumber(L, cfg.zero_threshold); lua_setfield(L, -2, "zero_threshold");
    lua_pushboolean(L, cfg.dedup);         lua_setfield(L, -2, "dedup");
    lua_pushinteger(L, cfg.log_level);     lua_setfield(L, -2, "log_level");
    lua_pushboolean(L, cfg.export_materials); lua_setfield(L, -2, "export_materials");
    lua_pushboolean(L, cfg.export_transforms);lua_setfield(L, -2, "export_transforms");
    lua_pushinteger(L, cfg.universe_depth);lua_setfield(L, -2, "universe_depth");
    lua_pushinteger(L, cfg.fill_depth);    lua_setfield(L, -2, "fill_depth");

    return 1;
}

/* Helper: read optional field from table at stack top */
#define READ_OPT_NUM(field) do { \
    lua_getfield(L, 2, #field); \
    if (!lua_isnil(L, -1)) cfg.field = luaL_checknumber(L, -1); \
    lua_pop(L, 1); \
} while(0)

#define READ_OPT_INT(field) do { \
    lua_getfield(L, 2, #field); \
    if (!lua_isnil(L, -1)) cfg.field = (int)luaL_checkinteger(L, -1); \
    lua_pop(L, 1); \
} while(0)

#define READ_OPT_BOOL(field) do { \
    lua_getfield(L, 2, #field); \
    if (!lua_isnil(L, -1)) cfg.field = lua_toboolean(L, -1); \
    lua_pop(L, 1); \
} while(0)

static int l_system_set_config(lua_State* L) {
    alea_system_t* sys = alea_get_sys(L, 1);
    luaL_checktype(L, 2, LUA_TTABLE);

    alea_config_t cfg = alea_get_config(sys);

    READ_OPT_NUM(abs_tol);
    READ_OPT_NUM(rel_tol);
    READ_OPT_NUM(zero_threshold);
    READ_OPT_BOOL(dedup);
    READ_OPT_INT(log_level);
    READ_OPT_BOOL(export_materials);
    READ_OPT_BOOL(export_transforms);
    READ_OPT_INT(universe_depth);
    READ_OPT_INT(fill_depth);

    alea_set_config(sys, &cfg);
    return 0;
}

#undef READ_OPT_NUM
#undef READ_OPT_INT
#undef READ_OPT_BOOL

/* ============================================================================
 * System methods table (populated by submodules)
 * ============================================================================ */

/* Forward: we register a shared method table in the registry */
static const luaL_Reg system_meta[] = {
    {"__gc",       l_system_gc},
    {"__tostring", l_system_tostring},
    {NULL, NULL}
};

static const luaL_Reg system_methods[] = {
    {"destroy",    l_system_destroy},
    {"reset",      l_system_reset},
    {"get_config", l_system_get_config},
    {"set_config", l_system_set_config},
    {NULL, NULL}
};

static const luaL_Reg alea_funcs[] = {
    {"create",          l_system_create},
    {"version",         l_version},
    {"log_set_level",   l_log_set_level},
    {"log_get_level",   l_log_get_level},
    {"error",           l_error},
    {"error_code",      l_error_code},
    {"error_clear",     l_error_clear},
    {"set_debug_trace", l_set_debug_trace},
    {NULL, NULL}
};

/* ============================================================================
 * Module opener
 * ============================================================================ */

int luaopen_alea(lua_State* L) {
    /* Create System metatable */
    luaL_newmetatable(L, ALEA_SYSTEM_MT);
    luaL_setfuncs(L, system_meta, 0);

    /* __index = methods table */
    lua_newtable(L);
    luaL_setfuncs(L, system_methods, 0);
    lua_setfield(L, -2, "__index");
    lua_pop(L, 1); /* pop metatable */

    /* Create the 'alea' global table */
    luaL_newlib(L, alea_funcs);

    /* Register submodules: each adds methods to System metatable
       and optionally functions to the alea table (top of stack) */
    luaopen_alea_ops(L);
    luaopen_alea_surfaces(L);
    luaopen_alea_io(L);
    luaopen_alea_query(L);
    luaopen_alea_util(L);
    luaopen_alea_raycast(L);
    luaopen_alea_slice(L);
    luaopen_alea_render(L);
    luaopen_alea_mesh(L);
    luaopen_alea_materials(L);

    /* Set as global */
    lua_pushvalue(L, -1);
    lua_setglobal(L, "alea");

    return 1;
}
