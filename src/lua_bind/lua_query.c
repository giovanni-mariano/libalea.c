// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

#include "alea_lua.h"

/* ============================================================================
 * Indexing
 * ============================================================================ */

static int l_build_universe_index(lua_State* L) {
    alea_system_t* sys = alea_get_sys(L, 1);
    if (alea_build_universe_index(sys) != 0)
        return luaL_error(L, "build_universe_index failed: %s", alea_error());
    return 0;
}

static int l_build_spatial_index(lua_State* L) {
    alea_system_t* sys = alea_get_sys(L, 1);
    if (alea_build_spatial_index(sys) != 0)
        return luaL_error(L, "build_spatial_index failed: %s", alea_error());
    return 0;
}

/* ============================================================================
 * Point queries
 * ============================================================================ */

/* sys:find_cell(x, y, z) -> cell_index or nil */
static int l_find_cell(lua_State* L) {
    alea_system_t* sys = alea_get_sys(L, 1);
    double x = luaL_checknumber(L, 2);
    double y = luaL_checknumber(L, 3);
    double z = luaL_checknumber(L, 4);
    int idx = alea_find_cell(sys, x, y, z);
    if (idx < 0) {
        lua_pushnil(L);
    } else {
        lua_pushinteger(L, idx);
    }
    return 1;
}

/* sys:material_at(x, y, z) -> material_id or nil */
static int l_material_at(lua_State* L) {
    alea_system_t* sys = alea_get_sys(L, 1);
    double x = luaL_checknumber(L, 2);
    double y = luaL_checknumber(L, 3);
    double z = luaL_checknumber(L, 4);
    int mat = alea_material_at(sys, x, y, z);
    if (mat < 0) {
        lua_pushnil(L);
    } else {
        lua_pushinteger(L, mat);
    }
    return 1;
}

/* sys:find_all_cells(x, y, z) -> table of {cell_id, material_id, universe_id, ...} */
static int l_find_all_cells(lua_State* L) {
    alea_system_t* sys = alea_get_sys(L, 1);
    double x = luaL_checknumber(L, 2);
    double y = luaL_checknumber(L, 3);
    double z = luaL_checknumber(L, 4);

    alea_cell_hit_t hits[256];
    int n = alea_find_all_cells(sys, x, y, z, hits, 256);
    if (n < 0) n = 0;

    lua_createtable(L, n, 0);
    for (int i = 0; i < n; i++) {
        lua_createtable(L, 0, 6);
        lua_pushinteger(L, hits[i].cell_id);       lua_setfield(L, -2, "cell_id");
        lua_pushinteger(L, hits[i].cell_index);     lua_setfield(L, -2, "cell_index");
        lua_pushinteger(L, hits[i].material_id);    lua_setfield(L, -2, "material_id");
        lua_pushinteger(L, hits[i].universe_id);    lua_setfield(L, -2, "universe_id");
        lua_pushinteger(L, hits[i].fill_universe);  lua_setfield(L, -2, "fill_universe");
        lua_pushinteger(L, hits[i].depth);          lua_setfield(L, -2, "depth");
        lua_rawseti(L, -2, i + 1);
    }
    return 1;
}

/* sys:point_inside(node, x, y, z) -> boolean */
static int l_point_inside(lua_State* L) {
    alea_system_t* sys = alea_get_sys(L, 1);
    alea_lua_node_t* nd = alea_check_node(L, 2);
    double x = luaL_checknumber(L, 3);
    double y = luaL_checknumber(L, 4);
    double z = luaL_checknumber(L, 5);
    bool inside = alea_point_inside(sys, nd->id, x, y, z);
    lua_pushboolean(L, inside);
    return 1;
}

/* sys:find_overlaps() -> table of {cell_a, cell_b} pairs */
static int l_find_overlaps(lua_State* L) {
    alea_system_t* sys = alea_get_sys(L, 1);
    int pairs[2048];
    int n = alea_find_overlaps(sys, pairs, 1024);
    if (n < 0) n = 0;

    lua_createtable(L, n, 0);
    for (int i = 0; i < n; i++) {
        lua_createtable(L, 2, 0);
        lua_pushinteger(L, pairs[2 * i]);
        lua_rawseti(L, -2, 1);
        lua_pushinteger(L, pairs[2 * i + 1]);
        lua_rawseti(L, -2, 2);
        lua_rawseti(L, -2, i + 1);
    }
    return 1;
}

/* ============================================================================
 * Information
 * ============================================================================ */

static int l_cell_count(lua_State* L) {
    alea_system_t* sys = alea_get_sys(L, 1);
    lua_pushinteger(L, (lua_Integer)alea_cell_count(sys));
    return 1;
}

static int l_surface_count(lua_State* L) {
    alea_system_t* sys = alea_get_sys(L, 1);
    lua_pushinteger(L, (lua_Integer)alea_surface_count(sys));
    return 1;
}

static int l_universe_count(lua_State* L) {
    alea_system_t* sys = alea_get_sys(L, 1);
    lua_pushinteger(L, (lua_Integer)alea_universe_count(sys));
    return 1;
}

static int l_print_summary(lua_State* L) {
    alea_system_t* sys = alea_get_sys(L, 1);
    alea_print_summary(sys);
    return 0;
}

static int l_validate(lua_State* L) {
    alea_system_t* sys = alea_get_sys(L, 1);
    int issues = alea_validate(sys);
    lua_pushinteger(L, issues);
    return 1;
}

/* sys:cell_info(index) -> table */
static int l_cell_info(lua_State* L) {
    alea_system_t* sys = alea_get_sys(L, 1);
    int index = (int)luaL_checkinteger(L, 2);
    alea_cell_info_t info;
    if (alea_cell_get_info(sys, (size_t)index, &info) != 0)
        return luaL_error(L, "cell_info: invalid index %d", index);

    lua_createtable(L, 0, 10);
    lua_pushinteger(L, info.cell_id);        lua_setfield(L, -2, "cell_id");
    lua_pushinteger(L, info.material_id);    lua_setfield(L, -2, "material_id");
    lua_pushnumber(L, info.density);         lua_setfield(L, -2, "density");
    lua_pushinteger(L, info.universe_id);    lua_setfield(L, -2, "universe_id");
    lua_pushinteger(L, info.fill_universe);  lua_setfield(L, -2, "fill_universe");
    lua_pushinteger(L, info.fill_transform); lua_setfield(L, -2, "fill_transform");
    lua_pushinteger(L, info.lat_type);       lua_setfield(L, -2, "lat_type");

    /* bbox subtable */
    lua_createtable(L, 0, 6);
    lua_pushnumber(L, info.bbox.min_x); lua_setfield(L, -2, "min_x");
    lua_pushnumber(L, info.bbox.max_x); lua_setfield(L, -2, "max_x");
    lua_pushnumber(L, info.bbox.min_y); lua_setfield(L, -2, "min_y");
    lua_pushnumber(L, info.bbox.max_y); lua_setfield(L, -2, "max_y");
    lua_pushnumber(L, info.bbox.min_z); lua_setfield(L, -2, "min_z");
    lua_pushnumber(L, info.bbox.max_z); lua_setfield(L, -2, "max_z");
    lua_setfield(L, -2, "bbox");

    /* root node */
    alea_push_node(L, sys, info.root);
    lua_setfield(L, -2, "root");

    return 1;
}

/* sys:cell_find(cell_id) -> index or nil */
static int l_cell_find(lua_State* L) {
    alea_system_t* sys = alea_get_sys(L, 1);
    int cell_id = (int)luaL_checkinteger(L, 2);
    int idx = alea_cell_find(sys, cell_id);
    if (idx < 0)
        lua_pushnil(L);
    else
        lua_pushinteger(L, idx);
    return 1;
}

/* sys:get_cell_id(cell_index) -> cell_id */
static int l_get_cell_id(lua_State* L) {
    alea_system_t* sys = alea_get_sys(L, 1);
    int index = (int)luaL_checkinteger(L, 2);
    int id = alea_get_cell_id(sys, index);
    lua_pushinteger(L, id);
    return 1;
}

/* sys:cells_by_material(mat_id) -> table of indices */
static int l_cells_by_material(lua_State* L) {
    alea_system_t* sys = alea_get_sys(L, 1);
    int mat_id = (int)luaL_checkinteger(L, 2);
    int buf[4096];
    size_t n = alea_get_cells_by_material(sys, mat_id, buf, 4096);

    lua_createtable(L, (int)n, 0);
    for (size_t i = 0; i < n; i++) {
        lua_pushinteger(L, buf[i]);
        lua_rawseti(L, -2, (lua_Integer)(i + 1));
    }
    return 1;
}

/* sys:cells_by_universe(universe_id) -> table of indices */
static int l_cells_by_universe(lua_State* L) {
    alea_system_t* sys = alea_get_sys(L, 1);
    int uid = (int)luaL_checkinteger(L, 2);
    int buf[4096];
    size_t n = alea_get_cells_by_universe(sys, uid, buf, 4096);

    lua_createtable(L, (int)n, 0);
    for (size_t i = 0; i < n; i++) {
        lua_pushinteger(L, buf[i]);
        lua_rawseti(L, -2, (lua_Integer)(i + 1));
    }
    return 1;
}

/* ============================================================================
 * Node inspection
 * ============================================================================ */

/* sys:node_operation(node) -> string */
static int l_node_operation(lua_State* L) {
    alea_system_t* sys = alea_get_sys(L, 1);
    alea_lua_node_t* nd = alea_check_node(L, 2);
    alea_operation_t op = alea_node_operation(sys, nd->id);
    const char* name;
    switch (op) {
        case ALEA_OP_PRIMITIVE:    name = "primitive";    break;
        case ALEA_OP_UNION:        name = "union";        break;
        case ALEA_OP_INTERSECTION: name = "intersection"; break;
        case ALEA_OP_DIFFERENCE:   name = "difference";   break;
        case ALEA_OP_COMPLEMENT:   name = "complement";   break;
        default:                   name = "unknown";      break;
    }
    lua_pushstring(L, name);
    return 1;
}

/* sys:node_left(node) -> Node or nil */
static int l_node_left(lua_State* L) {
    alea_system_t* sys = alea_get_sys(L, 1);
    alea_lua_node_t* nd = alea_check_node(L, 2);
    alea_push_node(L, sys, alea_node_left(sys, nd->id));
    return 1;
}

/* sys:node_right(node) -> Node or nil */
static int l_node_right(lua_State* L) {
    alea_system_t* sys = alea_get_sys(L, 1);
    alea_lua_node_t* nd = alea_check_node(L, 2);
    alea_push_node(L, sys, alea_node_right(sys, nd->id));
    return 1;
}

/* sys:node_sense(node) -> +1 or -1 or 0 */
static int l_node_sense(lua_State* L) {
    alea_system_t* sys = alea_get_sys(L, 1);
    alea_lua_node_t* nd = alea_check_node(L, 2);
    lua_pushinteger(L, alea_node_sense(sys, nd->id));
    return 1;
}

/* sys:node_surface_id(node) -> int */
static int l_node_surface_id(lua_State* L) {
    alea_system_t* sys = alea_get_sys(L, 1);
    alea_lua_node_t* nd = alea_check_node(L, 2);
    lua_pushinteger(L, alea_node_surface_id(sys, nd->id));
    return 1;
}

/* ============================================================================
 * Registration
 * ============================================================================ */

static const luaL_Reg query_methods[] = {
    {"build_universe_index", l_build_universe_index},
    {"build_spatial_index",  l_build_spatial_index},
    {"find_cell",            l_find_cell},
    {"material_at",          l_material_at},
    {"find_all_cells",       l_find_all_cells},
    {"point_inside",         l_point_inside},
    {"find_overlaps",        l_find_overlaps},
    {"cell_count",           l_cell_count},
    {"surface_count",        l_surface_count},
    {"universe_count",       l_universe_count},
    {"print_summary",        l_print_summary},
    {"validate",             l_validate},
    {"cell_info",            l_cell_info},
    {"cell_find",            l_cell_find},
    {"get_cell_id",          l_get_cell_id},
    {"cells_by_material",    l_cells_by_material},
    {"cells_by_universe",    l_cells_by_universe},
    {"node_operation",       l_node_operation},
    {"node_left",            l_node_left},
    {"node_right",           l_node_right},
    {"node_sense",           l_node_sense},
    {"node_surface_id",      l_node_surface_id},
    {NULL, NULL}
};

int luaopen_alea_query(lua_State* L) {
    luaL_getmetatable(L, ALEA_SYSTEM_MT);
    lua_getfield(L, -1, "__index");
    luaL_setfuncs(L, query_methods, 0);
    lua_pop(L, 2);
    return 0;
}
