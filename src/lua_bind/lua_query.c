// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

#include "alea_lua.h"
#include "core/alea_system.h"
#include <stdlib.h>

#if defined(__GNUC__) || defined(__clang__)
#define ALEA_LUA_DEPRECATED_CALL_BEGIN \
    _Pragma("GCC diagnostic push") \
    _Pragma("GCC diagnostic ignored \"-Wdeprecated-declarations\"")
#define ALEA_LUA_DEPRECATED_CALL_END \
    _Pragma("GCC diagnostic pop")
#else
#define ALEA_LUA_DEPRECATED_CALL_BEGIN
#define ALEA_LUA_DEPRECATED_CALL_END
#endif

/* ============================================================================
 * Indexing
 * ============================================================================ */

static int l_build_universe_index(lua_State* L) {
    alea_system_t* sys = alea_get_sys(L, 1);
    if (alea_build_universe_index(sys) != 0)
        return luaL_error(L, "build_universe_index failed: %s", alea_error());
    return 0;
}

static int l_prepare_query_acceleration(lua_State* L) {
    alea_system_t* sys = alea_get_sys(L, 1);
    luaL_checktype(L, 2, LUA_TNONE);
    if (alea_prepare_query_acceleration(sys) != 0)
        return luaL_error(L, "prepare_query_acceleration failed: %s", alea_error());
    return 0;
}

/* sys:query_acceleration_stats() -> table */
static int l_query_acceleration_stats(lua_State* L) {
    alea_system_t* sys = alea_get_sys(L, 1);
    alea_query_acceleration_stats_t stats;
    if (alea_query_acceleration_stats(sys, &stats) != 0)
        return luaL_error(L, "query_acceleration_stats failed: %s", alea_error());

    lua_createtable(L, 0, 21);
    lua_pushboolean(L, stats.built);                             lua_setfield(L, -2, "built");
    lua_pushinteger(L, (lua_Integer)stats.hier_universe_count);   lua_setfield(L, -2, "hier_universe_count");
    lua_pushinteger(L, (lua_Integer)stats.hier_blas_count);       lua_setfield(L, -2, "hier_blas_count");
    lua_pushinteger(L, (lua_Integer)stats.hier_linear_universe_count); lua_setfield(L, -2, "hier_linear_universe_count");
    lua_pushinteger(L, (lua_Integer)stats.hier_blas_cell_count);  lua_setfield(L, -2, "hier_blas_cell_count");
    lua_pushinteger(L, (lua_Integer)stats.hier_blas_node_count);  lua_setfield(L, -2, "hier_blas_node_count");
    lua_pushinteger(L, (lua_Integer)stats.hier_fill_cell_count);  lua_setfield(L, -2, "hier_fill_cell_count");
    lua_pushinteger(L, (lua_Integer)stats.hier_lattice_cell_count); lua_setfield(L, -2, "hier_lattice_cell_count");
    lua_pushinteger(L, (lua_Integer)stats.hier_transform_count);  lua_setfield(L, -2, "hier_transform_count");
    lua_pushinteger(L, (lua_Integer)stats.hier_placement_count);  lua_setfield(L, -2, "hier_placement_count");
    lua_pushinteger(L, (lua_Integer)stats.hier_root_placement_count); lua_setfield(L, -2, "hier_root_placement_count");
    lua_pushinteger(L, (lua_Integer)stats.hier_fill_placement_count); lua_setfield(L, -2, "hier_fill_placement_count");
    lua_pushinteger(L, (lua_Integer)stats.hier_lattice_placement_count); lua_setfield(L, -2, "hier_lattice_placement_count");
    lua_pushinteger(L, (lua_Integer)stats.hier_max_placement_depth); lua_setfield(L, -2, "hier_max_placement_depth");
    lua_pushinteger(L, (lua_Integer)stats.hier_max_universe_cells); lua_setfield(L, -2, "hier_max_universe_cells");
    lua_pushinteger(L, (lua_Integer)stats.hier_largest_universe_id); lua_setfield(L, -2, "hier_largest_universe_id");
    lua_pushinteger(L, (lua_Integer)stats.memory_bytes);          lua_setfield(L, -2, "memory_bytes");
    lua_pushinteger(L, (lua_Integer)stats.point_queries);         lua_setfield(L, -2, "point_queries");
    return 1;
}

/* ============================================================================
 * Point queries
 * ============================================================================ */

/* sys:find_cell(x, y, z) -> cell_index or nil
 * Legacy convenience wrapper. Prefer sys:find_cell_at(x, y, z). */
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

/* sys:material_at(x, y, z) -> material_id or nil
 * Legacy convenience wrapper. Prefer sys:find_cell_at(x, y, z). */
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
    if (info.has_temperature) {
        lua_pushnumber(L, info.temperature); lua_setfield(L, -2, "temperature");
    }

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

    /* comments */
    if (info.comments) {
        lua_pushstring(L, info.comments);  lua_setfield(L, -2, "comments");
    }
    if (info.inline_comment) {
        lua_pushstring(L, info.inline_comment); lua_setfield(L, -2, "inline_comment");
    }

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
 * Extended queries
 * ============================================================================ */

/* sys:find_cell_at(x,y,z) -> cell_id, material or nil */
static int l_find_cell_at(lua_State* L) {
    alea_system_t* sys = alea_get_sys(L, 1);
    double x = luaL_checknumber(L, 2);
    double y = luaL_checknumber(L, 3);
    double z = luaL_checknumber(L, 4);
    int cell_id, material;
    int rc = alea_find_cell_at(sys, x, y, z, &cell_id, &material);
    if (rc < 0) {
        lua_pushnil(L);
        return 1;
    }
    lua_pushinteger(L, cell_id);
    lua_pushinteger(L, material);
    return 2;
}

/* sys:cell_find_info(cell_id) -> info table or nil */
static int l_cell_find_info(lua_State* L) {
    alea_system_t* sys = alea_get_sys(L, 1);
    int cell_id = (int)luaL_checkinteger(L, 2);
    alea_cell_info_t info;
    if (alea_cell_find_info(sys, cell_id, &info) != 0) {
        lua_pushnil(L);
        return 1;
    }

    lua_createtable(L, 0, 10);
    lua_pushinteger(L, info.cell_id);        lua_setfield(L, -2, "cell_id");
    lua_pushinteger(L, info.material_id);    lua_setfield(L, -2, "material_id");
    lua_pushnumber(L, info.density);         lua_setfield(L, -2, "density");
    lua_pushinteger(L, info.universe_id);    lua_setfield(L, -2, "universe_id");
    lua_pushinteger(L, info.fill_universe);  lua_setfield(L, -2, "fill_universe");
    lua_pushinteger(L, info.fill_transform); lua_setfield(L, -2, "fill_transform");
    lua_pushinteger(L, info.lat_type);       lua_setfield(L, -2, "lat_type");
    if (info.has_temperature) {
        lua_pushnumber(L, info.temperature); lua_setfield(L, -2, "temperature");
    }

    lua_createtable(L, 0, 6);
    lua_pushnumber(L, info.bbox.min_x); lua_setfield(L, -2, "min_x");
    lua_pushnumber(L, info.bbox.max_x); lua_setfield(L, -2, "max_x");
    lua_pushnumber(L, info.bbox.min_y); lua_setfield(L, -2, "min_y");
    lua_pushnumber(L, info.bbox.max_y); lua_setfield(L, -2, "max_y");
    lua_pushnumber(L, info.bbox.min_z); lua_setfield(L, -2, "min_z");
    lua_pushnumber(L, info.bbox.max_z); lua_setfield(L, -2, "max_z");
    lua_setfield(L, -2, "bbox");

    alea_push_node(L, sys, info.root);
    lua_setfield(L, -2, "root");

    if (info.comments) {
        lua_pushstring(L, info.comments);  lua_setfield(L, -2, "comments");
    }
    if (info.inline_comment) {
        lua_pushstring(L, info.inline_comment); lua_setfield(L, -2, "inline_comment");
    }

    return 1;
}

/* sys:cells_in_universe(uid) -> table of indices */
static int l_cells_in_universe(lua_State* L) {
    alea_system_t* sys = alea_get_sys(L, 1);
    int uid = (int)luaL_checkinteger(L, 2);
    int buf[4096];
    int n = alea_cells_in_universe(sys, uid, buf, 4096);
    if (n < 0) n = 0;

    lua_createtable(L, n, 0);
    for (int i = 0; i < n; i++) {
        lua_pushinteger(L, buf[i]);
        lua_rawseti(L, -2, i + 1);
    }
    return 1;
}

/* sys:surface_info(idx) -> table */
static int l_surface_info(lua_State* L) {
    alea_system_t* sys = alea_get_sys(L, 1);
    int index = (int)luaL_checkinteger(L, 2);
    int surface_id;
    alea_primitive_type_t type;
    alea_node_id_t pos_node, neg_node;
    alea_boundary_type_t btype;
    if (alea_surface_get(sys, (size_t)index, &surface_id, &type,
                          &pos_node, &neg_node, &btype) != 0)
        return luaL_error(L, "surface_info: invalid index %d", index);

    lua_createtable(L, 0, 5);
    lua_pushinteger(L, surface_id);  lua_setfield(L, -2, "surface_id");
    lua_pushinteger(L, (int)type);   lua_setfield(L, -2, "type");
    lua_pushinteger(L, (int)btype);  lua_setfield(L, -2, "boundary_type");
    alea_push_node(L, sys, pos_node);
    lua_setfield(L, -2, "pos_node");
    alea_push_node(L, sys, neg_node);
    lua_setfield(L, -2, "neg_node");
    return 1;
}

/* sys:surface_find(sid) -> index or nil */
static int l_surface_find(lua_State* L) {
    alea_system_t* sys = alea_get_sys(L, 1);
    int sid = (int)luaL_checkinteger(L, 2);
    int idx = alea_surface_find(sys, sid);
    if (idx < 0)
        lua_pushnil(L);
    else
        lua_pushinteger(L, idx);
    return 1;
}

/* sys:universe_info(idx) -> table */
static int l_universe_info(lua_State* L) {
    alea_system_t* sys = alea_get_sys(L, 1);
    int index = (int)luaL_checkinteger(L, 2);
    int uid;
    size_t cc;
    alea_bbox_t bbox;
    if (alea_universe_get(sys, (size_t)index, &uid, &cc, &bbox) != 0)
        return luaL_error(L, "universe_info: invalid index %d", index);

    lua_createtable(L, 0, 3);
    lua_pushinteger(L, uid);          lua_setfield(L, -2, "universe_id");
    lua_pushinteger(L, (lua_Integer)cc); lua_setfield(L, -2, "cell_count");

    lua_createtable(L, 0, 6);
    lua_pushnumber(L, bbox.min_x); lua_setfield(L, -2, "min_x");
    lua_pushnumber(L, bbox.max_x); lua_setfield(L, -2, "max_x");
    lua_pushnumber(L, bbox.min_y); lua_setfield(L, -2, "min_y");
    lua_pushnumber(L, bbox.max_y); lua_setfield(L, -2, "max_y");
    lua_pushnumber(L, bbox.min_z); lua_setfield(L, -2, "min_z");
    lua_pushnumber(L, bbox.max_z); lua_setfield(L, -2, "max_z");
    lua_setfield(L, -2, "bbox");
    return 1;
}

/* sys:universe_find(uid) -> index or nil */
static int l_universe_find(lua_State* L) {
    alea_system_t* sys = alea_get_sys(L, 1);
    int uid = (int)luaL_checkinteger(L, 2);
    int idx = alea_universe_find(sys, uid);
    if (idx < 0)
        lua_pushnil(L);
    else
        lua_pushinteger(L, idx);
    return 1;
}

/* sys:node_primitive_type(node) -> string */
static int l_node_primitive_type(lua_State* L) {
    alea_system_t* sys = alea_get_sys(L, 1);
    alea_lua_node_t* nd = alea_check_node(L, 2);
    alea_primitive_type_t type = alea_node_primitive_type(sys, nd->id);
    const char* name;
    switch (type) {
        case ALEA_PRIMITIVE_PLANE:       name = "plane";       break;
        case ALEA_PRIMITIVE_SPHERE:      name = "sphere";      break;
        case ALEA_PRIMITIVE_CYLINDER_X:  name = "cylinder_x";  break;
        case ALEA_PRIMITIVE_CYLINDER_Y:  name = "cylinder_y";  break;
        case ALEA_PRIMITIVE_CYLINDER_Z:  name = "cylinder_z";  break;
        case ALEA_PRIMITIVE_CONE_X:      name = "cone_x";      break;
        case ALEA_PRIMITIVE_CONE_Y:      name = "cone_y";      break;
        case ALEA_PRIMITIVE_CONE_Z:      name = "cone_z";      break;
        case ALEA_PRIMITIVE_RPP:         name = "rpp";         break;
        case ALEA_PRIMITIVE_QUADRIC:     name = "quadric";     break;
        case ALEA_PRIMITIVE_TORUS_X:     name = "torus_x";     break;
        case ALEA_PRIMITIVE_TORUS_Y:     name = "torus_y";     break;
        case ALEA_PRIMITIVE_TORUS_Z:     name = "torus_z";     break;
        case ALEA_PRIMITIVE_RCC:         name = "rcc";         break;
        case ALEA_PRIMITIVE_BOX:         name = "box";         break;
        case ALEA_PRIMITIVE_SPH:         name = "sph";         break;
        case ALEA_PRIMITIVE_TRC:         name = "trc";         break;
        case ALEA_PRIMITIVE_ELL:         name = "ell";         break;
        case ALEA_PRIMITIVE_REC:         name = "rec";         break;
        case ALEA_PRIMITIVE_WED:         name = "wed";         break;
        case ALEA_PRIMITIVE_RHP:         name = "rhp";         break;
        case ALEA_PRIMITIVE_ARB:         name = "arb";         break;
        default:                         name = "unknown";     break;
    }
    lua_pushstring(L, name);
    return 1;
}

/* sys:node_primitive_id(node) -> int */
static int l_node_primitive_id(lua_State* L) {
    alea_system_t* sys = alea_get_sys(L, 1);
    alea_lua_node_t* nd = alea_check_node(L, 2);
    alea_primitive_id_t pid = alea_node_primitive_id(sys, nd->id);
    if (pid == ALEA_PRIMITIVE_ID_INVALID)
        lua_pushnil(L);
    else
        lua_pushinteger(L, (lua_Integer)pid);
    return 1;
}

/* sys:node_primitive_data(node) -> table */
static int l_node_primitive_data(lua_State* L) {
    alea_system_t* sys = alea_get_sys(L, 1);
    alea_lua_node_t* nd = alea_check_node(L, 2);
    alea_primitive_data_t data;
    if (alea_node_primitive_data(sys, nd->id, &data) != 0) {
        lua_pushnil(L);
        return 1;
    }

    alea_primitive_type_t type = alea_node_primitive_type(sys, nd->id);
    lua_newtable(L);

    switch (type) {
        case ALEA_PRIMITIVE_PLANE:
            lua_pushnumber(L, data.plane.a); lua_setfield(L, -2, "a");
            lua_pushnumber(L, data.plane.b); lua_setfield(L, -2, "b");
            lua_pushnumber(L, data.plane.c); lua_setfield(L, -2, "c");
            lua_pushnumber(L, data.plane.d); lua_setfield(L, -2, "d");
            break;
        case ALEA_PRIMITIVE_SPHERE:
            lua_pushnumber(L, data.sphere.center_x); lua_setfield(L, -2, "cx");
            lua_pushnumber(L, data.sphere.center_y); lua_setfield(L, -2, "cy");
            lua_pushnumber(L, data.sphere.center_z); lua_setfield(L, -2, "cz");
            lua_pushnumber(L, data.sphere.radius);   lua_setfield(L, -2, "r");
            break;
        case ALEA_PRIMITIVE_CYLINDER_Z:
            lua_pushnumber(L, data.cyl_z.center_x); lua_setfield(L, -2, "cx");
            lua_pushnumber(L, data.cyl_z.center_y); lua_setfield(L, -2, "cy");
            lua_pushnumber(L, data.cyl_z.radius);   lua_setfield(L, -2, "r");
            break;
        case ALEA_PRIMITIVE_CYLINDER_X:
            lua_pushnumber(L, data.cyl_x.center_y); lua_setfield(L, -2, "cy");
            lua_pushnumber(L, data.cyl_x.center_z); lua_setfield(L, -2, "cz");
            lua_pushnumber(L, data.cyl_x.radius);   lua_setfield(L, -2, "r");
            break;
        case ALEA_PRIMITIVE_CYLINDER_Y:
            lua_pushnumber(L, data.cyl_y.center_x); lua_setfield(L, -2, "cx");
            lua_pushnumber(L, data.cyl_y.center_z); lua_setfield(L, -2, "cz");
            lua_pushnumber(L, data.cyl_y.radius);   lua_setfield(L, -2, "r");
            break;
        case ALEA_PRIMITIVE_RPP:
            lua_pushnumber(L, data.box.min_x); lua_setfield(L, -2, "min_x");
            lua_pushnumber(L, data.box.max_x); lua_setfield(L, -2, "max_x");
            lua_pushnumber(L, data.box.min_y); lua_setfield(L, -2, "min_y");
            lua_pushnumber(L, data.box.max_y); lua_setfield(L, -2, "max_y");
            lua_pushnumber(L, data.box.min_z); lua_setfield(L, -2, "min_z");
            lua_pushnumber(L, data.box.max_z); lua_setfield(L, -2, "max_z");
            break;
        default:
            /* For other types, return type name only */
            break;
    }
    return 1;
}

/* sys:cells_filling_universe(uid) -> table */
static int l_cells_filling_universe(lua_State* L) {
    alea_system_t* sys = alea_get_sys(L, 1);
    int uid = (int)luaL_checkinteger(L, 2);
    int buf[4096];
    size_t n = alea_get_cells_filling_universe(sys, uid, buf, 4096);

    lua_createtable(L, (int)n, 0);
    for (size_t i = 0; i < n; i++) {
        lua_pushinteger(L, buf[i]);
        lua_rawseti(L, -2, (lua_Integer)(i + 1));
    }
    return 1;
}

/* sys:stats() -> table */
static int l_stats(lua_State* L) {
    alea_system_t* sys = alea_get_sys(L, 1);
    alea_stats_t stats;
    alea_stats(sys, &stats);

    lua_createtable(L, 0, 9);
    lua_pushinteger(L, (lua_Integer)stats.unique_primitives);    lua_setfield(L, -2, "unique_primitives");
    lua_pushinteger(L, (lua_Integer)stats.dedup_hits);           lua_setfield(L, -2, "dedup_hits");
    lua_pushinteger(L, (lua_Integer)stats.dedup_saved_bytes);    lua_setfield(L, -2, "dedup_saved_bytes");
    lua_pushinteger(L, (lua_Integer)stats.current_memory);       lua_setfield(L, -2, "current_memory");
    lua_pushinteger(L, (lua_Integer)stats.peak_memory);          lua_setfield(L, -2, "peak_memory");
    lua_pushinteger(L, (lua_Integer)stats.reallocations);        lua_setfield(L, -2, "reallocations");
    lua_pushinteger(L, (lua_Integer)stats.surfaces_converted);   lua_setfield(L, -2, "surfaces_converted");
    lua_pushinteger(L, (lua_Integer)stats.failed_surfaces);      lua_setfield(L, -2, "failed_surfaces");
    lua_pushinteger(L, (lua_Integer)stats.cells_converted);      lua_setfield(L, -2, "cells_converted");
    return 1;
}

/* sys:tree_print(node) */
static int l_tree_print(lua_State* L) {
    alea_system_t* sys = alea_get_sys(L, 1);
    alea_lua_node_t* nd = alea_check_node(L, 2);
    alea_tree_print(sys, nd->id);
    return 0;
}

/* ============================================================================
 * Registration
 * ============================================================================ */

static const luaL_Reg query_methods[] = {
    {"build_universe_index", l_build_universe_index},
    {"prepare_query_acceleration", l_prepare_query_acceleration},
    {"query_acceleration_stats", l_query_acceleration_stats},
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
    {"node_left",               l_node_left},
    {"node_right",              l_node_right},
    {"node_sense",              l_node_sense},
    {"node_surface_id",         l_node_surface_id},
    {"find_cell_at",            l_find_cell_at},
    {"cell_find_info",          l_cell_find_info},
    {"cells_in_universe",       l_cells_in_universe},
    {"surface_info",            l_surface_info},
    {"surface_find",            l_surface_find},
    {"universe_info",           l_universe_info},
    {"universe_find",           l_universe_find},
    {"node_primitive_type",     l_node_primitive_type},
    {"node_primitive_id",       l_node_primitive_id},
    {"node_primitive_data",     l_node_primitive_data},
    {"cells_filling_universe",  l_cells_filling_universe},
    {"stats",                   l_stats},
    {"tree_print",              l_tree_print},
    {NULL, NULL}
};

int luaopen_alea_query(lua_State* L) {
    luaL_getmetatable(L, ALEA_SYSTEM_MT);
    lua_getfield(L, -1, "__index");
    luaL_setfuncs(L, query_methods, 0);
    lua_pop(L, 2);
    return 0;
}
