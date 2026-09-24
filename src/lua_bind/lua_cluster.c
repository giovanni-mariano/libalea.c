// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
// SPDX-License-Identifier: MPL-2.0

#include "alea_lua.h"
#include "alea_cluster.h"

#include <stdlib.h>

#define ALEA_CLUSTER_MT "alea.ClusterContext"

typedef struct {
    alea_cluster_t* ptr;
    int busy;
} alea_lua_cluster_t;

static int runtime_active;
static int open_contexts;

static alea_lua_cluster_t* check_cluster(lua_State* L, int index) {
    alea_lua_cluster_t* ctx = luaL_checkudata(L, index, ALEA_CLUSTER_MT);
    if (!runtime_active) luaL_error(L, "cluster runtime is not initialized");
    if (!ctx->ptr) luaL_error(L, "cluster context is closed");
    if (ctx->busy) luaL_error(L, "cluster operation already active");
    return ctx;
}

static int cluster_error(lua_State* L, const char* operation,
                         alea_cluster_status_t status) {
    return luaL_error(L, "%s failed: %s", operation,
                      alea_cluster_status_string(status));
}

static int l_cluster_initialize(lua_State* L) {
    (void)L;
    if (runtime_active) return 0;
    alea_cluster_status_t status = alea_cluster_initialize(NULL, NULL);
    if (status != ALEA_CLUSTER_OK)
        return cluster_error(L, "cluster_initialize", status);
    runtime_active = 1;
    return 0;
}

static int l_cluster_finalize(lua_State* L) {
    if (!runtime_active) return luaL_error(L, "cluster runtime is not initialized");
    if (open_contexts)
        return luaL_error(L, "close every cluster context before finalizing");
    alea_cluster_status_t status = alea_cluster_finalize();
    if (status != ALEA_CLUSTER_OK)
        return cluster_error(L, "cluster_finalize", status);
    runtime_active = 0;
    return 0;
}

static int l_cluster_create(lua_State* L) {
    if (!runtime_active) return luaL_error(L, "call cluster_initialize first");
    if (open_contexts) return luaL_error(L, "only one cluster context may be open");
    alea_cluster_t* native = alea_cluster_create();
    if (!native) return luaL_error(L, "collective cluster creation failed");
    alea_lua_cluster_t* ctx = lua_newuserdatauv(L, sizeof(*ctx), 0);
    ctx->ptr = native;
    ctx->busy = 0;
    luaL_setmetatable(L, ALEA_CLUSTER_MT);
    open_contexts++;
    return 1;
}

static int l_cluster_close(lua_State* L) {
    alea_lua_cluster_t* ctx = check_cluster(L, 1);
    alea_cluster_destroy(ctx->ptr);
    ctx->ptr = NULL;
    open_contexts--;
    return 0;
}

static int l_cluster_gc(lua_State* L) {
    alea_lua_cluster_t* ctx = luaL_checkudata(L, 1, ALEA_CLUSTER_MT);
    if (ctx->ptr)
        lua_warning(L, "unclosed cluster context leaked; call close() collectively", 0);
    return 0;
}

static int l_cluster_rank(lua_State* L) {
    lua_pushinteger(L, alea_cluster_rank(check_cluster(L, 1)->ptr));
    return 1;
}

static int l_cluster_size(lua_State* L) {
    lua_pushinteger(L, alea_cluster_size(check_cluster(L, 1)->ptr));
    return 1;
}

static int l_cluster_is_root(lua_State* L) {
    lua_pushboolean(L, alea_cluster_is_root(check_cluster(L, 1)->ptr));
    return 1;
}

static int l_cluster_backend(lua_State* L) {
    lua_pushstring(L, alea_cluster_backend(check_cluster(L, 1)->ptr));
    return 1;
}

static int l_cluster_agree(lua_State* L) {
    alea_lua_cluster_t* ctx = check_cluster(L, 1);
    alea_cluster_status_t local = (alea_cluster_status_t)luaL_optinteger(
        L, 2, ALEA_CLUSTER_OK);
    ctx->busy = 1;
    alea_cluster_status_t status = alea_cluster_agree(ctx->ptr, local);
    ctx->busy = 0;
    lua_pushinteger(L, status);
    lua_pushstring(L, alea_cluster_status_string(status));
    return 2;
}

static int cluster_read(lua_State* L, int mcnp) {
    alea_lua_cluster_t* ctx = check_cluster(L, 1);
    const char* path = NULL;
    alea_cluster_status_t local = ALEA_CLUSTER_OK;
    if (alea_cluster_is_root(ctx->ptr)) {
        if (lua_type(L, 2) == LUA_TSTRING) path = lua_tostring(L, 2);
        else local = ALEA_CLUSTER_INVALID_ARGUMENT;
    }
    char* data = NULL;
    size_t length = 0;
    alea_cluster_status_t status = alea_cluster_agree(ctx->ptr, local);
    if (status != ALEA_CLUSTER_OK)
        return cluster_error(L, mcnp ? "cluster read_mcnp_input preparation"
                                     : "cluster read_file preparation", status);
    ctx->busy = 1;
    status = mcnp
        ? alea_cluster_read_mcnp_input(ctx->ptr, path, &data, &length)
        : alea_cluster_read_file(ctx->ptr, path, &data, &length);
    ctx->busy = 0;
    if (status != ALEA_CLUSTER_OK) {
        free(data);
        return cluster_error(L, mcnp ? "cluster read_mcnp_input" : "cluster read_file", status);
    }
    alea_lua_free_guard_t* guard = alea_lua_push_free_guard(L);
    guard->ptr = data;
    lua_pushlstring(L, data, length);
    free(data);
    guard->ptr = NULL;
    lua_remove(L, -2);
    return 1;
}

static int l_cluster_read_file(lua_State* L) { return cluster_read(L, 0); }
static int l_cluster_read_mcnp(lua_State* L) { return cluster_read(L, 1); }

static size_t option_size(lua_State* L, int table, const char* name, size_t fallback) {
    lua_getfield(L, table, name);
    lua_Integer value = luaL_optinteger(L, -1, (lua_Integer)fallback);
    lua_pop(L, 1);
    if (value < 0) luaL_error(L, "%s must be non-negative", name);
    return (size_t)value;
}

static double option_number(lua_State* L, int table, const char* name, double fallback) {
    lua_getfield(L, table, name);
    double value = luaL_optnumber(L, -1, fallback);
    lua_pop(L, 1);
    return value;
}

static int l_cluster_estimate_volumes(lua_State* L) {
    alea_lua_cluster_t* ctx = check_cluster(L, 1);
    alea_system_t* sys = alea_get_sys(L, 2);
    luaL_checktype(L, 3, LUA_TTABLE);
    alea_volume_estimate_options_t options;
    alea_volume_estimate_options_init(&options);
    options.max_rays = option_size(L, 3, "max_rays", 0);
    options.seed = (uint64_t)option_size(L, 3, "seed", 1);
    options.requested_workers = option_size(L, 3, "workers", 0);
    options.batch_size = option_size(L, 3, "batch_size", 0);
    options.max_parallel_scratch_bytes =
        option_size(L, 3, "max_parallel_scratch_bytes",
                    options.max_parallel_scratch_bytes);
    options.target_rel_error = option_number(L, 3, "target_rel_error", 0.0);
    if (!options.max_rays) return luaL_error(L, "max_rays must be positive");

    size_t count = alea_volume_path_count(sys);
    double* volumes = count ? calloc(count, sizeof(*volumes)) : NULL;
    double* errors = count ? calloc(count, sizeof(*errors)) : NULL;
    alea_cluster_status_t local = count && volumes && errors
        ? ALEA_CLUSTER_OK : (count ? ALEA_CLUSTER_OUT_OF_MEMORY
                                   : ALEA_CLUSTER_INVALID_ARGUMENT);
    alea_cluster_status_t status = alea_cluster_agree(ctx->ptr, local);
    if (status != ALEA_CLUSTER_OK) {
        free(volumes); free(errors);
        return cluster_error(L, "cluster estimate_volumes preparation", status);
    }

    alea_cluster_volume_stats_t stats;
    ctx->busy = 1;
    status = alea_cluster_estimate_volumes(ctx->ptr, sys, &options,
                                           volumes, errors, &stats);
    ctx->busy = 0;
    if (status != ALEA_CLUSTER_OK) {
        free(volumes); free(errors);
        return cluster_error(L, "cluster estimate_volumes", status);
    }

    lua_createtable(L, 0, 10);
    lua_createtable(L, (int)count, 0);
    for (size_t i = 0; i < count; ++i) {
        lua_pushnumber(L, volumes[i]);
        lua_rawseti(L, -2, (lua_Integer)i + 1);
    }
    lua_setfield(L, -2, "volumes");
    lua_createtable(L, (int)count, 0);
    for (size_t i = 0; i < count; ++i) {
        lua_pushnumber(L, errors[i]);
        lua_rawseti(L, -2, (lua_Integer)i + 1);
    }
    lua_setfield(L, -2, "rel_errors");
    free(volumes); free(errors);
#define SET_INTEGER(name, value) do { lua_pushinteger(L, (lua_Integer)(value)); lua_setfield(L, -2, name); } while (0)
    SET_INTEGER("rays_completed", stats.volume.rays_completed);
    SET_INTEGER("actual_workers", stats.volume.actual_workers);
    SET_INTEGER("rank_count", stats.rank_count);
    SET_INTEGER("local_rays_completed", stats.local_rays_completed);
    SET_INTEGER("local_workers", stats.local_workers);
    SET_INTEGER("seed", stats.volume.seed);
    SET_INTEGER("parallel_scratch_limit_bytes",
                stats.volume.parallel_scratch_limit_bytes);
    SET_INTEGER("parallel_scratch_bytes",
                stats.volume.parallel_scratch_bytes);
    SET_INTEGER("worker_scratch_bytes", stats.volume.worker_scratch_bytes);
#undef SET_INTEGER
    lua_pushnumber(L, stats.volume.maximum_relative_error);
    lua_setfield(L, -2, "maximum_relative_error");
    lua_pushboolean(L, stats.volume.converged); lua_setfield(L, -2, "converged");
    lua_pushboolean(L, stats.volume.cancelled); lua_setfield(L, -2, "cancelled");
    return 1;
}

static int try_read_vec3(lua_State* L, int index, double out[3]) {
    index = lua_absindex(L, index);
    if (!lua_istable(L, index)) return 0;
    for (int i = 0; i < 3; ++i) {
        lua_geti(L, index, i + 1);
        if (!lua_isnumber(L, -1)) {
            lua_pop(L, 1);
            return 0;
        }
        out[i] = lua_tonumber(L, -1);
        lua_pop(L, 1);
    }
    return 1;
}

/* ctx:raycast_first_segments(system, rays [, t_max]) -> results on root, nil elsewhere
 * Each ray is {origin={x,y,z}, direction={x,y,z}}. */
static int l_cluster_raycast_first_segments(lua_State* L) {
    alea_lua_cluster_t* ctx = check_cluster(L, 1);
    alea_system_t* sys = alea_get_sys(L, 2);
    size_t count = 0;
    double t_max = 0.0;
    double *origins = NULL, *directions = NULL, *t_enter = NULL, *t_exit = NULL;
    unsigned char* hit = NULL;
    int32_t* cell_ids = NULL;
    alea_cluster_status_t local = ALEA_CLUSTER_OK;

    if (alea_cluster_is_root(ctx->ptr)) {
        if (!lua_istable(L, 3) ||
            (!lua_isnoneornil(L, 4) && !lua_isnumber(L, 4)))
            local = ALEA_CLUSTER_INVALID_ARGUMENT;
        if (local == ALEA_CLUSTER_OK) {
            count = lua_rawlen(L, 3);
            t_max = luaL_optnumber(L, 4, 0.0);
        }
        if (count > SIZE_MAX / (3 * sizeof(double)))
            local = ALEA_CLUSTER_INVALID_ARGUMENT;
        if (local == ALEA_CLUSTER_OK && count) {
            origins = malloc(count * 3 * sizeof(*origins));
            directions = malloc(count * 3 * sizeof(*directions));
            hit = malloc(count * sizeof(*hit));
            cell_ids = malloc(count * sizeof(*cell_ids));
            t_enter = malloc(count * sizeof(*t_enter));
            t_exit = malloc(count * sizeof(*t_exit));
            if (!origins || !directions || !hit || !cell_ids || !t_enter || !t_exit)
                local = ALEA_CLUSTER_OUT_OF_MEMORY;
        }
        for (size_t i = 0; local == ALEA_CLUSTER_OK && i < count; ++i) {
            lua_geti(L, 3, (lua_Integer)i + 1);
            if (!lua_istable(L, -1)) local = ALEA_CLUSTER_INVALID_ARGUMENT;
            else {
                lua_getfield(L, -1, "origin");
                if (!try_read_vec3(L, -1, origins + 3 * i))
                    local = ALEA_CLUSTER_INVALID_ARGUMENT;
                lua_pop(L, 1);
                lua_getfield(L, -1, "direction");
                if (!try_read_vec3(L, -1, directions + 3 * i))
                    local = ALEA_CLUSTER_INVALID_ARGUMENT;
                lua_pop(L, 1);
            }
            lua_pop(L, 1);
        }
    }

    alea_cluster_status_t status = alea_cluster_agree(ctx->ptr, local);
    if (status == ALEA_CLUSTER_OK) {
        ctx->busy = 1;
        status = alea_cluster_raycast_first_segments(ctx->ptr, sys,
            origins, directions, count, t_max, hit, cell_ids, t_enter, t_exit);
        ctx->busy = 0;
    }
    if (status != ALEA_CLUSTER_OK) {
        free(origins); free(directions); free(hit); free(cell_ids);
        free(t_enter); free(t_exit);
        return cluster_error(L, "cluster raycast_first_segments", status);
    }
    if (!alea_cluster_is_root(ctx->ptr)) {
        lua_pushnil(L);
        return 1;
    }

    lua_createtable(L, (int)count, 0);
    for (size_t i = 0; i < count; ++i) {
        if (!hit[i]) lua_pushboolean(L, 0);
        else {
            lua_createtable(L, 0, 3);
            lua_pushinteger(L, cell_ids[i]); lua_setfield(L, -2, "cell_id");
            lua_pushnumber(L, t_enter[i]); lua_setfield(L, -2, "t_enter");
            lua_pushnumber(L, t_exit[i]); lua_setfield(L, -2, "t_exit");
        }
        lua_rawseti(L, -2, (lua_Integer)i + 1);
    }
    free(origins); free(directions); free(hit); free(cell_ids);
    free(t_enter); free(t_exit);
    return 1;
}

static const luaL_Reg cluster_methods[] = {
    {"close", l_cluster_close},
    {"rank", l_cluster_rank},
    {"size", l_cluster_size},
    {"is_root", l_cluster_is_root},
    {"backend", l_cluster_backend},
    {"agree", l_cluster_agree},
    {"read_file", l_cluster_read_file},
    {"read_mcnp_input", l_cluster_read_mcnp},
    {"estimate_volumes", l_cluster_estimate_volumes},
    {"raycast_first_segments", l_cluster_raycast_first_segments},
    {NULL, NULL}
};

static const luaL_Reg cluster_module[] = {
    {"initialize", l_cluster_initialize},
    {"finalize", l_cluster_finalize},
    {"create", l_cluster_create},
    {NULL, NULL}
};

int luaopen_alea_cluster(lua_State* L) {
    luaL_requiref(L, "alea", luaopen_alea, 1);
    lua_pop(L, 1);
    luaL_newmetatable(L, ALEA_CLUSTER_MT);
    lua_pushcfunction(L, l_cluster_gc); lua_setfield(L, -2, "__gc");
    lua_newtable(L);
    luaL_setfuncs(L, cluster_methods, 0);
    lua_setfield(L, -2, "__index");
    lua_pop(L, 1);
    luaL_newlib(L, cluster_module);
    return 1;
}
