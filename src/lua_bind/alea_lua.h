// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

#ifndef ALEA_LUA_H
#define ALEA_LUA_H

#include "lua.h"
#include "lauxlib.h"
#include "lualib.h"
#include "alea.h"
#include "alea_slice.h"

/* ============================================================================
 * Metatable names
 * ============================================================================ */

#define ALEA_SYSTEM_MT      "alea.System"
#define ALEA_NODE_MT        "alea.Node"
#define ALEA_RAYRESULT_MT   "alea.RaycastResult"
#define ALEA_CURVES_MT      "alea.SliceCurves"
#define ALEA_MESHRESULT_MT  "alea.MeshResult"
#define ALEA_VOIDRESULT_MT  "alea.VoidResult"
#define ALEA_FRAMEBUFFER_MT "alea.Framebuffer"
#define ALEA_GEOMRESULT_MT  "alea.GeometryValidationResult"

/* ============================================================================
 * Userdata types
 * ============================================================================ */

typedef struct {
    alea_system_t* sys;
    void* mcnp_model;    /* mcnp_model_t* if loaded from MCNP, NULL otherwise */
    void* openmc_model;  /* openmc_model_t* if loaded from OpenMC, NULL otherwise */
    int owned;            /* 1 if we should destroy on __gc */
    int destroy_pending;  /* destroy requested while dependent userdata exists */
    int active_void_results;
} alea_lua_system_t;

typedef struct {
    alea_system_t* sys; /* borrowed pointer */
    alea_node_id_t id;
} alea_lua_node_t;

/* ============================================================================
 * Helpers
 * ============================================================================ */

static inline alea_lua_system_t* alea_check_system(lua_State* L, int idx) {
    return (alea_lua_system_t*)luaL_checkudata(L, idx, ALEA_SYSTEM_MT);
}

static inline void alea_lua_system_init(alea_lua_system_t* ud) {
    ud->sys = NULL;
    ud->mcnp_model = NULL;
    ud->openmc_model = NULL;
    ud->owned = 1;
    ud->destroy_pending = 0;
    ud->active_void_results = 0;
}

static inline alea_system_t* alea_get_sys(lua_State* L, int idx) {
    alea_lua_system_t* ud = alea_check_system(L, idx);
    if (!ud->sys || ud->destroy_pending)
        luaL_error(L, "system has been destroyed");
    return ud->sys;
}

void alea_lua_system_release_if_pending(alea_lua_system_t* ud);

static inline alea_lua_node_t* alea_check_node(lua_State* L, int idx) {
    return (alea_lua_node_t*)luaL_checkudata(L, idx, ALEA_NODE_MT);
}

static inline void alea_push_node(lua_State* L, alea_system_t* sys, alea_node_id_t id) {
    if (id == ALEA_NODE_ID_INVALID) {
        lua_pushnil(L);
        return;
    }
    alea_lua_node_t* ud = (alea_lua_node_t*)lua_newuserdata(L, sizeof(alea_lua_node_t));
    ud->sys = sys;
    ud->id = id;
    luaL_setmetatable(L, ALEA_NODE_MT);
}

/* ============================================================================
 * Submodule registration functions
 * ============================================================================ */

int luaopen_alea_ops(lua_State* L);
int luaopen_alea_surfaces(lua_State* L);
int luaopen_alea_io(lua_State* L);
int luaopen_alea_query(lua_State* L);
int luaopen_alea_util(lua_State* L);
int luaopen_alea_raycast(lua_State* L);
int luaopen_alea_geo_validator(lua_State* L);
int luaopen_alea_slice(lua_State* L);

/* Read a Lua slice-view table into a view struct (defined in lua_slice.c). */
void alea_lua_to_slice_view(lua_State* L, int idx, alea_slice_view_t* view);
int luaopen_alea_render(lua_State* L);
int luaopen_alea_mesh(lua_State* L);
int luaopen_alea_materials(lua_State* L);
int luaopen_alea_nucdata(lua_State* L);
int luaopen_alea_svg(lua_State* L);

/* The main module opener (called from lua_main.c) */
int luaopen_alea(lua_State* L);

#endif /* ALEA_LUA_H */
