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
#define ALEA_DIRECTIONAL_TRACE_CACHE_MT "alea.DirectionalTraceCache"
#define ALEA_FREE_GUARD_MT   "alea.FreeGuard"

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
    int active_directional_trace_caches;
    uint64_t node_generation;
} alea_lua_system_t;

typedef struct {
    alea_lua_system_t* owner; /* retained by uservalue */
    alea_node_id_t id;
    uint64_t generation;
} alea_lua_node_t;

typedef struct { void* ptr; } alea_lua_free_guard_t;

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
    ud->active_directional_trace_caches = 0;
    ud->node_generation = 1;
}

static inline alea_system_t* alea_get_sys(lua_State* L, int idx) {
    alea_lua_system_t* ud = alea_check_system(L, idx);
    if (!ud->sys || ud->destroy_pending)
        luaL_error(L, "system has been destroyed");
    return ud->sys;
}

static inline void alea_lua_require_no_dependents(lua_State* L,
                                                   alea_lua_system_t* ud,
                                                   const char* operation) {
    if (ud->active_void_results || ud->active_directional_trace_caches)
        luaL_error(L, "cannot %s while dependent results exist", operation);
}

void alea_lua_system_release_if_pending(alea_lua_system_t* ud);

static inline alea_lua_free_guard_t* alea_lua_push_free_guard(lua_State* L) {
    alea_lua_free_guard_t* guard = lua_newuserdatauv(L, sizeof(*guard), 0);
    guard->ptr = NULL;
    luaL_setmetatable(L, ALEA_FREE_GUARD_MT);
    return guard;
}

static inline alea_lua_node_t* alea_check_node(lua_State* L, int idx) {
    alea_lua_node_t* node = (alea_lua_node_t*)luaL_checkudata(L, idx, ALEA_NODE_MT);
    if (!node->owner || !node->owner->sys || node->owner->destroy_pending)
        luaL_error(L, "node belongs to a destroyed system");
    if (node->generation != node->owner->node_generation)
        luaL_error(L, "node is stale after system mutation");
    return node;
}

static inline alea_lua_node_t* alea_check_node_for_system(
    lua_State* L, int node_idx, int system_idx)
{
    alea_lua_node_t* node = alea_check_node(L, node_idx);
    alea_lua_system_t* owner = alea_check_system(L, system_idx);
    if (node->owner != owner)
        luaL_error(L, "node belongs to a different system");
    return node;
}

static inline alea_lua_node_t* alea_check_node_same_owner(
    lua_State* L, int node_idx, int other_node_idx)
{
    alea_lua_node_t* node = alea_check_node(L, node_idx);
    alea_lua_node_t* other = alea_check_node(L, other_node_idx);
    if (node->owner != other->owner)
        luaL_error(L, "nodes belong to different systems");
    return node;
}

static inline void alea_push_node(lua_State* L, int owner_idx, alea_node_id_t id) {
    if (id == ALEA_NODE_ID_INVALID) {
        lua_pushnil(L);
        return;
    }
    owner_idx = lua_absindex(L, owner_idx);
    alea_lua_system_t* owner = alea_check_system(L, owner_idx);
    alea_lua_node_t* ud = (alea_lua_node_t*)lua_newuserdatauv(
        L, sizeof(alea_lua_node_t), 1);
    ud->owner = owner;
    ud->id = id;
    ud->generation = owner->node_generation;
    luaL_setmetatable(L, ALEA_NODE_MT);
    lua_pushvalue(L, owner_idx);
    lua_setiuservalue(L, -2, 1);
}

static inline void alea_push_node_from_node(lua_State* L, int node_idx,
                                            alea_node_id_t id) {
    node_idx = lua_absindex(L, node_idx);
    (void)alea_check_node(L, node_idx);
    lua_getiuservalue(L, node_idx, 1);
    alea_push_node(L, -1, id);
    lua_remove(L, -2);
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
