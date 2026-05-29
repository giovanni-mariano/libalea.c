// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

#include "alea_lua.h"

#include "alea_geo_validator.h"

#include <stdlib.h>

typedef struct {
    alea_geom_validator_result_t* result;
} alea_lua_geom_result_t;

static alea_lua_geom_result_t* check_geom_result(lua_State* L, int idx) {
    return (alea_lua_geom_result_t*)luaL_checkudata(L, idx, ALEA_GEOMRESULT_MT);
}

static int get_bool_field(lua_State* L, int table_idx, const char* name, int* out) {
    lua_getfield(L, table_idx, name);
    if (lua_isnil(L, -1)) {
        lua_pop(L, 1);
        return 0;
    }
    *out = lua_toboolean(L, -1);
    lua_pop(L, 1);
    return 1;
}

static void parse_validator_options(lua_State* L,
                                    int idx,
                                    alea_geom_validator_options_t* options) {
    alea_geom_validator_options_init(options);
    if (idx > lua_gettop(L) || lua_isnoneornil(L, idx))
        return;

    idx = lua_absindex(L, idx);
    luaL_checktype(L, idx, LUA_TTABLE);

    lua_getfield(L, idx, "ray_count");
    if (!lua_isnil(L, -1))
        options->ray_count = (int)luaL_checkinteger(L, -1);
    lua_pop(L, 1);

    lua_getfield(L, idx, "seed");
    if (!lua_isnil(L, -1))
        options->seed = (uint64_t)luaL_checkinteger(L, -1);
    lua_pop(L, 1);

    lua_getfield(L, idx, "max_errors");
    if (!lua_isnil(L, -1))
        options->max_errors = (size_t)luaL_checkinteger(L, -1);
    lua_pop(L, 1);

    lua_getfield(L, idx, "max_crossings");
    if (!lua_isnil(L, -1))
        options->max_crossings = (size_t)luaL_checkinteger(L, -1);
    lua_pop(L, 1);

    lua_getfield(L, idx, "sample_offset");
    if (!lua_isnil(L, -1))
        options->sample_offset = luaL_checknumber(L, -1);
    lua_pop(L, 1);

    lua_getfield(L, idx, "t_max");
    if (!lua_isnil(L, -1))
        options->t_max = luaL_checknumber(L, -1);
    lua_pop(L, 1);

    lua_getfield(L, idx, "universe_depth");
    if (!lua_isnil(L, -1))
        options->universe_depth = (int)luaL_checkinteger(L, -1);
    lua_pop(L, 1);

    int value = 0;
    if (get_bool_field(L, idx, "allow_exterior_void", &value) ||
        get_bool_field(L, idx, "allow_exterior", &value)) {
        if (value)
            options->flags |= ALEA_GEOM_VALIDATE_ALLOW_EXTERIOR_VOID;
        else
            options->flags &= ~ALEA_GEOM_VALIDATE_ALLOW_EXTERIOR_VOID;
    }

    if (get_bool_field(L, idx, "hierarchical", &value) ||
        get_bool_field(L, idx, "hier", &value)) {
        if (value)
            options->flags |= ALEA_GEOM_VALIDATE_HIERARCHICAL;
        else
            options->flags &= ~ALEA_GEOM_VALIDATE_HIERARCHICAL;
    }

    if (get_bool_field(L, idx, "strict", &value)) {
        if (value)
            options->flags |= ALEA_GEOM_VALIDATE_STRICT_ADJACENCY;
        else
            options->flags &= ~ALEA_GEOM_VALIDATE_STRICT_ADJACENCY;
    }
}

static void push_vec3(lua_State* L, const double v[3]) {
    lua_createtable(L, 3, 0);
    for (int i = 0; i < 3; i++) {
        lua_pushnumber(L, v[i]);
        lua_rawseti(L, -2, i + 1);
    }
}

static void push_geom_error(lua_State* L, const alea_geom_error_t* error) {
    lua_createtable(L, 0, 15);
    lua_pushstring(L, alea_geom_error_type_name(error->type));
    lua_setfield(L, -2, "type");
    lua_pushinteger(L, error->type);
    lua_setfield(L, -2, "type_id");
    lua_pushinteger(L, error->previous_cell_id);
    lua_setfield(L, -2, "previous_cell_id");
    lua_pushinteger(L, error->found_cell_id);
    lua_setfield(L, -2, "found_cell_id");
    lua_pushinteger(L, error->expected_neighbor_cell_id);
    lua_setfield(L, -2, "expected_neighbor_cell_id");
    lua_pushinteger(L, error->secondary_cell_id);
    lua_setfield(L, -2, "secondary_cell_id");
    lua_pushinteger(L, error->found_cell_count);
    lua_setfield(L, -2, "found_cell_count");
    lua_pushinteger(L, error->surface_id);
    lua_setfield(L, -2, "surface_id");
    lua_pushinteger(L, error->universe_id);
    lua_setfield(L, -2, "universe_id");
    lua_pushinteger(L, error->universe_depth);
    lua_setfield(L, -2, "universe_depth");
    lua_pushnumber(L, error->t);
    lua_setfield(L, -2, "t");
    lua_pushnumber(L, error->offset);
    lua_setfield(L, -2, "offset");
    lua_pushinteger(L, error->flags);
    lua_setfield(L, -2, "flags");
    push_vec3(L, error->crossing_point);
    lua_setfield(L, -2, "crossing_point");
    push_vec3(L, error->sample_point);
    lua_setfield(L, -2, "sample_point");
    push_vec3(L, error->direction);
    lua_setfield(L, -2, "direction");
}

static alea_lua_geom_result_t* push_new_result(lua_State* L) {
    alea_lua_geom_result_t* ud =
        (alea_lua_geom_result_t*)lua_newuserdata(L, sizeof(*ud));
    ud->result = (alea_geom_validator_result_t*)calloc(1, sizeof(*ud->result));
    if (!ud->result)
        luaL_error(L, "validate_geometry: out of memory");
    alea_geom_validator_result_init(ud->result);
    luaL_setmetatable(L, ALEA_GEOMRESULT_MT);
    return ud;
}

/* sys:validate_geometry(options?) -> GeometryValidationResult */
static int l_validate_geometry(lua_State* L) {
    alea_system_t* sys = alea_get_sys(L, 1);
    alea_geom_validator_options_t options;
    parse_validator_options(L, 2, &options);

    alea_lua_geom_result_t* ud = push_new_result(L);
    if (alea_validate_geometry(sys, &options, ud->result) != 0) {
        alea_geom_validator_result_free(ud->result);
        free(ud->result);
        ud->result = NULL;
        return luaL_error(L, "validate_geometry: %s", alea_error());
    }
    return 1;
}

/* sys:validate_geometry_ray(ox,oy,oz, dx,dy,dz, t_max, options?) */
static int l_validate_geometry_ray(lua_State* L) {
    alea_system_t* sys = alea_get_sys(L, 1);
    double ox = luaL_checknumber(L, 2);
    double oy = luaL_checknumber(L, 3);
    double oz = luaL_checknumber(L, 4);
    double dx = luaL_checknumber(L, 5);
    double dy = luaL_checknumber(L, 6);
    double dz = luaL_checknumber(L, 7);
    double t_max = luaL_checknumber(L, 8);

    alea_geom_validator_options_t options;
    parse_validator_options(L, 9, &options);

    alea_lua_geom_result_t* ud = push_new_result(L);
    if (alea_validate_geometry_ray(sys, &options, ox, oy, oz, dx, dy, dz,
                                   t_max, ud->result) != 0) {
        alea_geom_validator_result_free(ud->result);
        free(ud->result);
        ud->result = NULL;
        return luaL_error(L, "validate_geometry_ray: %s", alea_error());
    }
    return 1;
}

static int l_geom_error_count(lua_State* L) {
    alea_lua_geom_result_t* ud = check_geom_result(L, 1);
    lua_pushinteger(L, (lua_Integer)alea_geom_validator_error_count(ud->result));
    return 1;
}

static int l_geom_error(lua_State* L) {
    alea_lua_geom_result_t* ud = check_geom_result(L, 1);
    size_t index = (size_t)luaL_checkinteger(L, 2);
    if (index == 0)
        return luaL_error(L, "error index is 1-based");

    alea_geom_error_t error;
    if (alea_geom_validator_error_get(ud->result, index - 1, &error) != 0) {
        lua_pushnil(L);
        return 1;
    }
    push_geom_error(L, &error);
    return 1;
}

static int l_geom_errors(lua_State* L) {
    alea_lua_geom_result_t* ud = check_geom_result(L, 1);
    size_t n = alea_geom_validator_error_count(ud->result);
    lua_createtable(L, (int)n, 0);
    for (size_t i = 0; i < n; i++) {
        alea_geom_error_t error;
        if (alea_geom_validator_error_get(ud->result, i, &error) == 0) {
            push_geom_error(L, &error);
            lua_rawseti(L, -2, (lua_Integer)i + 1);
        }
    }
    return 1;
}

static int l_geom_summary(lua_State* L) {
    alea_lua_geom_result_t* ud = check_geom_result(L, 1);
    lua_createtable(L, 0, 6);
    for (size_t i = 0; i < alea_geom_validator_error_count(ud->result); i++) {
        alea_geom_error_t error;
        if (alea_geom_validator_error_get(ud->result, i, &error) != 0)
            continue;
        const char* name = alea_geom_error_type_name(error.type);
        lua_getfield(L, -1, name);
        lua_Integer count = lua_isnil(L, -1) ? 0 : lua_tointeger(L, -1);
        lua_pop(L, 1);
        lua_pushinteger(L, count + 1);
        lua_setfield(L, -2, name);
    }
    lua_pushinteger(L, (lua_Integer)alea_geom_validator_error_count(ud->result));
    lua_setfield(L, -2, "total");
    return 1;
}

static int l_geom_stats(lua_State* L) {
    alea_lua_geom_result_t* ud = check_geom_result(L, 1);
    const alea_geom_validator_result_t* result = ud->result;

    lua_createtable(L, 0, 7);
    lua_pushinteger(L, (lua_Integer)result->error_count);
    lua_setfield(L, -2, "error_count");
    lua_pushinteger(L, (lua_Integer)result->crossings_checked);
    lua_setfield(L, -2, "crossings_checked");
    lua_pushinteger(L, (lua_Integer)result->adjacency_hits);
    lua_setfield(L, -2, "adjacency_hits");
    lua_pushinteger(L, (lua_Integer)result->exact_queries);
    lua_setfield(L, -2, "exact_queries");
    lua_pushinteger(L, (lua_Integer)result->ambiguous_crossings);
    lua_setfield(L, -2, "ambiguous_crossings");
    lua_pushboolean(L, result->truncated);
    lua_setfield(L, -2, "truncated");
    return 1;
}

static int l_geom_gc(lua_State* L) {
    alea_lua_geom_result_t* ud = check_geom_result(L, 1);
    if (ud->result) {
        alea_geom_validator_result_free(ud->result);
        free(ud->result);
        ud->result = NULL;
    }
    return 0;
}

static const luaL_Reg geom_system_methods[] = {
    {"validate_geometry",     l_validate_geometry},
    {"validate_geometry_ray", l_validate_geometry_ray},
    {NULL, NULL}
};

static const luaL_Reg geom_result_methods[] = {
    {"error_count", l_geom_error_count},
    {"error",       l_geom_error},
    {"errors",      l_geom_errors},
    {"summary",     l_geom_summary},
    {"stats",       l_geom_stats},
    {NULL, NULL}
};

static const luaL_Reg geom_result_meta[] = {
    {"__gc", l_geom_gc},
    {NULL, NULL}
};

int luaopen_alea_geo_validator(lua_State* L) {
    luaL_getmetatable(L, ALEA_SYSTEM_MT);
    lua_getfield(L, -1, "__index");
    luaL_setfuncs(L, geom_system_methods, 0);
    lua_pop(L, 2);

    luaL_newmetatable(L, ALEA_GEOMRESULT_MT);
    luaL_setfuncs(L, geom_result_meta, 0);
    lua_newtable(L);
    luaL_setfuncs(L, geom_result_methods, 0);
    lua_setfield(L, -2, "__index");
    lua_pop(L, 1);

    return 0;
}
