// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

#include "alea_lua.h"

/* ============================================================================
 * Node operators
 * ============================================================================ */

static int l_node_mul(lua_State* L) {
    alea_lua_node_t* a = alea_check_node(L, 1);
    alea_lua_node_t* b = alea_check_node(L, 2);
    alea_node_id_t result = alea_intersection(a->sys, a->id, b->id);
    if (result == ALEA_NODE_ID_INVALID)
        return luaL_error(L, "intersection failed: %s", alea_error());
    alea_push_node(L, a->sys, result);
    return 1;
}

static int l_node_add(lua_State* L) {
    alea_lua_node_t* a = alea_check_node(L, 1);
    alea_lua_node_t* b = alea_check_node(L, 2);
    alea_node_id_t result = alea_union(a->sys, a->id, b->id);
    if (result == ALEA_NODE_ID_INVALID)
        return luaL_error(L, "union failed: %s", alea_error());
    alea_push_node(L, a->sys, result);
    return 1;
}

static int l_node_sub(lua_State* L) {
    alea_lua_node_t* a = alea_check_node(L, 1);
    alea_lua_node_t* b = alea_check_node(L, 2);
    alea_node_id_t result = alea_difference(a->sys, a->id, b->id);
    if (result == ALEA_NODE_ID_INVALID)
        return luaL_error(L, "difference failed: %s", alea_error());
    alea_push_node(L, a->sys, result);
    return 1;
}

static int l_node_bnot(lua_State* L) {
    alea_lua_node_t* a = alea_check_node(L, 1);
    alea_node_id_t result = alea_complement(a->sys, a->id);
    if (result == ALEA_NODE_ID_INVALID)
        return luaL_error(L, "complement failed: %s", alea_error());
    alea_push_node(L, a->sys, result);
    return 1;
}

static int l_node_tostring(lua_State* L) {
    alea_lua_node_t* nd = alea_check_node(L, 1);
    alea_operation_t op = alea_node_operation(nd->sys, nd->id);
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
    alea_push_node(L, sys, id);
    return 1;
}

/* sys:inside(surface_index) -> Node (negative sense) */
static int l_system_inside(lua_State* L) {
    alea_system_t* sys = alea_get_sys(L, 1);
    int idx = (int)luaL_checkinteger(L, 2);
    alea_node_id_t id = alea_halfspace(sys, idx, -1);
    if (id == ALEA_NODE_ID_INVALID)
        return luaL_error(L, "inside failed: invalid surface index %d", idx);
    alea_push_node(L, sys, id);
    return 1;
}

/* sys:outside(surface_index) -> Node (positive sense) */
static int l_system_outside(lua_State* L) {
    alea_system_t* sys = alea_get_sys(L, 1);
    int idx = (int)luaL_checkinteger(L, 2);
    alea_node_id_t id = alea_halfspace(sys, idx, +1);
    if (id == ALEA_NODE_ID_INVALID)
        return luaL_error(L, "outside failed: invalid surface index %d", idx);
    alea_push_node(L, sys, id);
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
    alea_lua_node_t* node = (alea_lua_node_t*)luaL_checkudata(L, -1, ALEA_NODE_MT);
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
        alea_lua_node_t* nd = (alea_lua_node_t*)luaL_checkudata(L, -1, ALEA_NODE_MT);
        ids[i - 1] = nd->id;
        lua_pop(L, 1);
    }

    alea_node_id_t result = alea_union_n(sys, ids, (size_t)n);
    if (result == ALEA_NODE_ID_INVALID)
        return luaL_error(L, "union_n failed: %s", alea_error());
    alea_push_node(L, sys, result);
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
        alea_lua_node_t* nd = (alea_lua_node_t*)luaL_checkudata(L, -1, ALEA_NODE_MT);
        ids[i - 1] = nd->id;
        lua_pop(L, 1);
    }

    alea_node_id_t result = alea_intersection_n(sys, ids, (size_t)n);
    if (result == ALEA_NODE_ID_INVALID)
        return luaL_error(L, "intersection_n failed: %s", alea_error());
    alea_push_node(L, sys, result);
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
