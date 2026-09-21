// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

#include "alea_lua.h"
#include <math.h>

static int read_transform_values(lua_State* L, int index, double values[13]) {
    luaL_checktype(L, index, LUA_TTABLE);
    lua_Integer count = luaL_len(L, index);
    if (count < 3 || count > 13)
        luaL_argerror(L, index,
            "values must contain between 3 and 13 MCNP transform values");
    for (lua_Integer i = 1; i <= count; ++i) {
        lua_geti(L, index, i);
        values[i - 1] = luaL_checknumber(L, -1);
        lua_pop(L, 1);
        if (!isfinite(values[i - 1]))
            luaL_argerror(L, index, "transform values must be finite");
    }
    return (int)count;
}

static int option_boolean(lua_State* L, int table, const char* field,
                          int default_value) {
    if (lua_isnoneornil(L, table)) return default_value;
    luaL_checktype(L, table, LUA_TTABLE);
    lua_getfield(L, table, field);
    int result = lua_isnil(L, -1) ? default_value : lua_toboolean(L, -1);
    lua_pop(L, 1);
    return result;
}

/* sys:add_transform(id, values[, options]) -> id */
static int l_system_add_transform(lua_State* L) {
    alea_system_t* sys = alea_get_sys(L, 1);
    int transform_id = (int)luaL_checkinteger(L, 2);
    double values[13];
    int count = read_transform_values(L, 3, values);
    int degrees = option_boolean(L, 4, "degrees", 0);
    if (alea_add_transform(sys, transform_id, values, count, degrees) != 0)
        return luaL_error(L, "add_transform: %s", alea_error());
    lua_pushinteger(L, transform_id);
    return 1;
}

/* sys:add_inline_transform(values[, options]) -> assigned id */
static int l_system_add_inline_transform(lua_State* L) {
    alea_system_t* sys = alea_get_sys(L, 1);
    double values[13];
    int count = read_transform_values(L, 2, values);
    int degrees = option_boolean(L, 3, "degrees", 0);
    int cell_id = 0;
    const char* role = "fill";
    if (!lua_isnoneornil(L, 3)) {
        lua_getfield(L, 3, "cell_id");
        if (!lua_isnil(L, -1)) cell_id = (int)luaL_checkinteger(L, -1);
        lua_pop(L, 1);
        lua_getfield(L, 3, "role");
        if (!lua_isnil(L, -1)) role = luaL_checkstring(L, -1);
        lua_pop(L, 1);
    }
    int transform_id = alea_add_inline_transform(
        sys, values, count, degrees, cell_id, role);
    if (transform_id < 0)
        return luaL_error(L, "add_inline_transform: %s", alea_error());
    lua_pushinteger(L, transform_id);
    return 1;
}

static int l_system_set_comment(lua_State* L) {
    alea_system_t* sys = alea_get_sys(L, 1);
    int cell = (int)luaL_checkinteger(L, 2);
    const char* value = luaL_optstring(L, 3, NULL);
    if (alea_cell_set_comment(sys, cell, value) != 0)
        return luaL_error(L, "set_comment: invalid cell index %d", cell);
    return 0;
}

static int l_system_set_inline_comment(lua_State* L) {
    alea_system_t* sys = alea_get_sys(L, 1);
    int cell = (int)luaL_checkinteger(L, 2);
    const char* value = luaL_optstring(L, 3, NULL);
    if (alea_cell_set_inline_comment(sys, cell, value) != 0)
        return luaL_error(L, "set_inline_comment: invalid cell index %d", cell);
    return 0;
}

static int l_system_cell_set_material(lua_State* L) {
    alea_system_t* sys = alea_get_sys(L, 1);
    int cell = (int)luaL_checkinteger(L, 2);
    int material = (int)luaL_checkinteger(L, 3);
    if (alea_cell_set_material(sys, cell, material) != 0)
        return luaL_error(L, "cell_set_material: invalid cell or material index");
    return 0;
}

static int l_system_cell_set_density(lua_State* L) {
    alea_system_t* sys = alea_get_sys(L, 1);
    int cell = (int)luaL_checkinteger(L, 2);
    double density = luaL_checknumber(L, 3);
    if (!isfinite(density)) return luaL_argerror(L, 3, "density must be finite");
    if (alea_cell_set_density(sys, cell, density) != 0)
        return luaL_error(L, "cell_set_density: invalid cell index %d", cell);
    return 0;
}

static int l_system_cell_set_temperature(lua_State* L) {
    alea_system_t* sys = alea_get_sys(L, 1);
    int cell = (int)luaL_checkinteger(L, 2);
    double temperature = luaL_checknumber(L, 3);
    if (alea_cell_set_temperature(sys, cell, temperature) != 0)
        return luaL_error(L,
            "cell_set_temperature: temperature must be finite and positive and the cell must exist");
    return 0;
}

static int l_system_cell_clear_temperature(lua_State* L) {
    alea_system_t* sys = alea_get_sys(L, 1);
    int cell = (int)luaL_checkinteger(L, 2);
    if (alea_cell_clear_temperature(sys, cell) != 0)
        return luaL_error(L, "cell_clear_temperature: invalid cell index %d", cell);
    return 0;
}

static int l_system_cell_set_universe(lua_State* L) {
    alea_system_t* sys = alea_get_sys(L, 1);
    int cell = (int)luaL_checkinteger(L, 2);
    int universe = (int)luaL_checkinteger(L, 3);
    if (alea_cell_set_universe(sys, cell, universe) != 0)
        return luaL_error(L, "cell_set_universe: invalid cell index %d", cell);
    return 0;
}

static int l_system_cell_remove(lua_State* L) {
    alea_system_t* sys = alea_get_sys(L, 1);
    int cell = (int)luaL_checkinteger(L, 2);
    if (alea_cell_remove(sys, cell) != 0)
        return luaL_error(L, "cell_remove: invalid cell index %d", cell);
    return 0;
}

/* ============================================================================
 * Node operators
 * ============================================================================ */

static int l_node_mul(lua_State* L) {
    alea_lua_node_t* a = alea_check_node(L, 1);
    alea_lua_node_t* b = alea_check_node_same_owner(L, 2, 1);
    alea_node_id_t result = alea_intersection(a->owner->sys, a->id, b->id);
    if (result == ALEA_NODE_ID_INVALID)
        return luaL_error(L, "intersection failed: %s", alea_error());
    alea_push_node_from_node(L, 1, result);
    return 1;
}

static int l_node_add(lua_State* L) {
    alea_lua_node_t* a = alea_check_node(L, 1);
    alea_lua_node_t* b = alea_check_node_same_owner(L, 2, 1);
    alea_node_id_t result = alea_union(a->owner->sys, a->id, b->id);
    if (result == ALEA_NODE_ID_INVALID)
        return luaL_error(L, "union failed: %s", alea_error());
    alea_push_node_from_node(L, 1, result);
    return 1;
}

static int l_node_sub(lua_State* L) {
    alea_lua_node_t* a = alea_check_node(L, 1);
    alea_lua_node_t* b = alea_check_node_same_owner(L, 2, 1);
    alea_node_id_t result = alea_difference(a->owner->sys, a->id, b->id);
    if (result == ALEA_NODE_ID_INVALID)
        return luaL_error(L, "difference failed: %s", alea_error());
    alea_push_node_from_node(L, 1, result);
    return 1;
}

static int l_node_bnot(lua_State* L) {
    alea_lua_node_t* a = alea_check_node(L, 1);
    alea_node_id_t result = alea_complement(a->owner->sys, a->id);
    if (result == ALEA_NODE_ID_INVALID)
        return luaL_error(L, "complement failed: %s", alea_error());
    alea_push_node_from_node(L, 1, result);
    return 1;
}

static int l_node_tostring(lua_State* L) {
    alea_lua_node_t* nd = (alea_lua_node_t*)luaL_checkudata(L, 1, ALEA_NODE_MT);
    if (!nd->owner || !nd->owner->sys || nd->owner->destroy_pending) {
        lua_pushliteral(L, "Node(destroyed)");
        return 1;
    }
    if (nd->generation != nd->owner->node_generation) {
        lua_pushfstring(L, "Node(stale, id=%d)", (int)nd->id);
        return 1;
    }
    alea_operation_t op = alea_node_operation(nd->owner->sys, nd->id);
    const char* name;
    switch (op) {
        case ALEA_OP_PRIMITIVE:    name = "PRIMITIVE";    break;
        case ALEA_OP_UNION:        name = "UNION";        break;
        case ALEA_OP_INTERSECTION: name = "INTERSECTION"; break;
        case ALEA_OP_DIFFERENCE:   name = "DIFFERENCE";   break;
        case ALEA_OP_COMPLEMENT:   name = "COMPLEMENT";   break;
        default:                   name = "UNKNOWN";      break;
    }
    lua_pushfstring(L, "Node(%s, id=%d)", name, (int)nd->id);
    return 1;
}

/* ============================================================================
 * System methods for node creation
 * ============================================================================ */

/* sys:half(surface_index, sense) -> Node */
static int l_system_half(lua_State* L) {
    alea_system_t* sys = alea_get_sys(L, 1);
    int idx = (int)luaL_checkinteger(L, 2);
    int sense = (int)luaL_checkinteger(L, 3);
    alea_node_id_t id = alea_halfspace(sys, idx, sense);
    if (id == ALEA_NODE_ID_INVALID)
        return luaL_error(L, "halfspace failed: invalid surface index %d", idx);
    alea_push_node(L, 1, id);
    return 1;
}

/* sys:inside(surface_index) -> Node (negative sense) */
static int l_system_inside(lua_State* L) {
    alea_system_t* sys = alea_get_sys(L, 1);
    int idx = (int)luaL_checkinteger(L, 2);
    alea_node_id_t id = alea_halfspace(sys, idx, -1);
    if (id == ALEA_NODE_ID_INVALID)
        return luaL_error(L, "inside failed: invalid surface index %d", idx);
    alea_push_node(L, 1, id);
    return 1;
}

/* sys:outside(surface_index) -> Node (positive sense) */
static int l_system_outside(lua_State* L) {
    alea_system_t* sys = alea_get_sys(L, 1);
    int idx = (int)luaL_checkinteger(L, 2);
    alea_node_id_t id = alea_halfspace(sys, idx, +1);
    if (id == ALEA_NODE_ID_INVALID)
        return luaL_error(L, "outside failed: invalid surface index %d", idx);
    alea_push_node(L, 1, id);
    return 1;
}

/* sys:material(id) -> material_index */
static int l_system_material(lua_State* L) {
    alea_system_t* sys = alea_get_sys(L, 1);
    int mat_id = (int)luaL_checkinteger(L, 2);
    int idx = alea_add_material(sys, mat_id);
    if (idx < 0)
        return luaL_error(L, "add_material failed: %s", alea_error());
    lua_pushinteger(L, idx);
    return 1;
}

/* sys:cell{id=N, region=Node, material=M, density=D, universe=U} */
static int l_system_cell(lua_State* L) {
    alea_system_t* sys = alea_get_sys(L, 1);
    luaL_checktype(L, 2, LUA_TTABLE);

    /* Required: region */
    lua_getfield(L, 2, "region");
    alea_lua_node_t* node = alea_check_node_for_system(L, -1, 1);
    alea_node_id_t root = node->id;
    lua_pop(L, 1);

    /* Optional fields */
    lua_getfield(L, 2, "id");
    int cell_id = lua_isnil(L, -1) ? 0 : (int)lua_tointeger(L, -1);
    lua_pop(L, 1);

    lua_getfield(L, 2, "material");
    int material = lua_isnil(L, -1) ? ALEA_MATERIAL_VOID : (int)lua_tointeger(L, -1);
    lua_pop(L, 1);

    lua_getfield(L, 2, "density");
    double density = lua_isnil(L, -1) ? 0.0 : lua_tonumber(L, -1);
    lua_pop(L, 1);

    lua_getfield(L, 2, "universe");
    int universe = lua_isnil(L, -1) ? 0 : (int)lua_tointeger(L, -1);
    lua_pop(L, 1);

    int idx = alea_add_cell(sys, cell_id, root, material, density, universe);
    if (idx < 0)
        return luaL_error(L, "add_cell failed: %s", alea_error());

    /* Optional: fill */
    lua_getfield(L, 2, "fill");
    if (!lua_isnil(L, -1)) {
        int fill_u = (int)lua_tointeger(L, -1);
        lua_pop(L, 1);
        lua_getfield(L, 2, "fill_transform");
        int fill_tr = lua_isnil(L, -1) ? 0 : (int)lua_tointeger(L, -1);
        lua_pop(L, 1);
        if (alea_set_fill(sys, idx, fill_u, fill_tr) < 0)
            return luaL_error(L, "set_fill failed: %s", alea_error());
    } else {
        lua_pop(L, 1);
    }

    /* Optional: comments */
    lua_getfield(L, 2, "comment");
    if (lua_isstring(L, -1)) {
        alea_cell_set_comment(sys, idx, lua_tostring(L, -1));
    }
    lua_pop(L, 1);

    lua_getfield(L, 2, "inline_comment");
    if (lua_isstring(L, -1)) {
        alea_cell_set_inline_comment(sys, idx, lua_tostring(L, -1));
    }
    lua_pop(L, 1);

    lua_pushinteger(L, idx);
    return 1;
}

/* sys:set_fill(cell_index, fill_universe, transform) */
static int l_system_set_fill(lua_State* L) {
    alea_system_t* sys = alea_get_sys(L, 1);
    int cell_idx = (int)luaL_checkinteger(L, 2);
    int fill_u = (int)luaL_checkinteger(L, 3);
    int transform = (int)luaL_optinteger(L, 4, 0);
    if (alea_set_fill(sys, cell_idx, fill_u, transform) < 0)
        return luaL_error(L, "set_fill failed: %s", alea_error());
    return 0;
}

/* alea.union_n(sys, {node1, node2, ...}) -> Node */
static int l_union_n(lua_State* L) {
    alea_system_t* sys = alea_get_sys(L, 1);
    luaL_checktype(L, 2, LUA_TTABLE);
    lua_Integer n = luaL_len(L, 2);
    if (n < 1)
        return luaL_error(L, "union_n: empty table");

    alea_node_id_t* ids = (alea_node_id_t*)lua_newuserdata(L, sizeof(alea_node_id_t) * (size_t)n);
    for (lua_Integer i = 1; i <= n; i++) {
        lua_geti(L, 2, i);
        alea_lua_node_t* nd = alea_check_node_for_system(L, -1, 1);
        ids[i - 1] = nd->id;
        lua_pop(L, 1);
    }

    alea_node_id_t result = alea_union_n(sys, ids, (size_t)n);
    if (result == ALEA_NODE_ID_INVALID)
        return luaL_error(L, "union_n failed: %s", alea_error());
    alea_push_node(L, 1, result);
    return 1;
}

/* alea.intersection_n(sys, {node1, node2, ...}) -> Node */
static int l_intersection_n(lua_State* L) {
    alea_system_t* sys = alea_get_sys(L, 1);
    luaL_checktype(L, 2, LUA_TTABLE);
    lua_Integer n = luaL_len(L, 2);
    if (n < 1)
        return luaL_error(L, "intersection_n: empty table");

    alea_node_id_t* ids = (alea_node_id_t*)lua_newuserdata(L, sizeof(alea_node_id_t) * (size_t)n);
    for (lua_Integer i = 1; i <= n; i++) {
        lua_geti(L, 2, i);
        alea_lua_node_t* nd = alea_check_node_for_system(L, -1, 1);
        ids[i - 1] = nd->id;
        lua_pop(L, 1);
    }

    alea_node_id_t result = alea_intersection_n(sys, ids, (size_t)n);
    if (result == ALEA_NODE_ID_INVALID)
        return luaL_error(L, "intersection_n failed: %s", alea_error());
    alea_push_node(L, 1, result);
    return 1;
}

/* ============================================================================
 * Registration
 * ============================================================================ */

static const luaL_Reg node_meta[] = {
    {"__mul",      l_node_mul},
    {"__add",      l_node_add},
    {"__sub",      l_node_sub},
    {"__bnot",     l_node_bnot},
    {"__tostring", l_node_tostring},
    {NULL, NULL}
};

static const luaL_Reg system_ops_methods[] = {
    {"half",     l_system_half},
    {"inside",   l_system_inside},
    {"outside",  l_system_outside},
    {"material", l_system_material},
    {"cell",     l_system_cell},
    {"set_fill", l_system_set_fill},
    {"add_transform", l_system_add_transform},
    {"add_inline_transform", l_system_add_inline_transform},
    {"set_comment", l_system_set_comment},
    {"set_inline_comment", l_system_set_inline_comment},
    {"cell_set_material", l_system_cell_set_material},
    {"cell_set_density", l_system_cell_set_density},
    {"cell_set_temperature", l_system_cell_set_temperature},
    {"cell_clear_temperature", l_system_cell_clear_temperature},
    {"cell_set_universe", l_system_cell_set_universe},
    {"cell_remove", l_system_cell_remove},
    {NULL, NULL}
};

int luaopen_alea_ops(lua_State* L) {
    /* Create Node metatable */
    luaL_newmetatable(L, ALEA_NODE_MT);
    luaL_setfuncs(L, node_meta, 0);
    lua_pop(L, 1);

    /* Add ops methods to System metatable's __index */
    luaL_getmetatable(L, ALEA_SYSTEM_MT);
    lua_getfield(L, -1, "__index");
    luaL_setfuncs(L, system_ops_methods, 0);
    lua_pop(L, 2);

    /* Add union_n, intersection_n to alea table (top of stack) */
    lua_pushcfunction(L, l_union_n);
    lua_setfield(L, -2, "union_n");
    lua_pushcfunction(L, l_intersection_n);
    lua_setfield(L, -2, "intersection_n");

    return 0;
}
