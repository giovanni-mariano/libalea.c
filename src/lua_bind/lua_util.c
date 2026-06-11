// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

#include "alea_lua.h"
#include "core/alea_system.h"
#include <stdlib.h>
#include <string.h>

#if defined(__GNUC__) || defined(__clang__)
#define ALEA_LUA_UTIL_DEPRECATED_CALL_BEGIN \
    _Pragma("GCC diagnostic push") \
    _Pragma("GCC diagnostic ignored \"-Wdeprecated-declarations\"")
#define ALEA_LUA_UTIL_DEPRECATED_CALL_END \
    _Pragma("GCC diagnostic pop")
#else
#define ALEA_LUA_UTIL_DEPRECATED_CALL_BEGIN
#define ALEA_LUA_UTIL_DEPRECATED_CALL_END
#endif

/* ============================================================================
 * Universe operations
 * ============================================================================ */

/* sys:flatten(universe_id) */
static int l_flatten(lua_State* L) {
    alea_system_t* sys = alea_get_sys(L, 1);
    int uid = (int)luaL_optinteger(L, 2, 0);
    if (alea_flatten(sys, uid) < 0)
        return luaL_error(L, "flatten failed: %s", alea_error());
    return 0;
}

/* sys:split_union_cells() -> number of new cells */
static int l_split_union_cells(lua_State* L) {
    alea_system_t* sys = alea_get_sys(L, 1);
    int n = alea_split_union_cells(sys);
    if (n < 0)
        return luaL_error(L, "split_union_cells failed: %s", alea_error());
    lua_pushinteger(L, n);
    return 1;
}

/* sys:extract_universe(id) -> new System */
static int l_extract_universe(lua_State* L) {
    alea_system_t* sys = alea_get_sys(L, 1);
    int uid = (int)luaL_checkinteger(L, 2);
    alea_lua_system_t* ud = (alea_lua_system_t*)lua_newuserdata(L, sizeof(alea_lua_system_t));
    alea_lua_system_init(ud);
    luaL_setmetatable(L, ALEA_SYSTEM_MT);
    ud->sys = alea_extract_universe(sys, uid);
    if (!ud->sys)
        return luaL_error(L, "extract_universe failed: %s", alea_error());
    return 1;
}

/* sys:merge(other, offset) */
static int l_merge(lua_State* L) {
    alea_system_t* sys = alea_get_sys(L, 1);
    alea_system_t* other = alea_get_sys(L, 2);
    int offset = (int)luaL_optinteger(L, 3, 0);
    if (alea_merge(sys, other, offset) < 0)
        return luaL_error(L, "merge failed: %s", alea_error());
    return 0;
}

/* sys:clone() -> new System */
static int l_clone(lua_State* L) {
    alea_system_t* sys = alea_get_sys(L, 1);
    alea_lua_system_t* ud = (alea_lua_system_t*)lua_newuserdata(L, sizeof(alea_lua_system_t));
    alea_lua_system_init(ud);
    luaL_setmetatable(L, ALEA_SYSTEM_MT);
    ud->sys = alea_clone(sys);
    if (!ud->sys)
        return luaL_error(L, "clone failed: %s", alea_error());
    return 1;
}

/* ============================================================================
 * Volume estimation
 * ============================================================================ */

/* sys:bounding_sphere(tol) -> cx, cy, cz, radius */
static int l_bounding_sphere(lua_State* L) {
    alea_system_t* sys = alea_get_sys(L, 1);
    double tol = luaL_optnumber(L, 2, 1.0);
    double cx, cy, cz, r;
    if (alea_compute_bounding_sphere(sys, tol, &cx, &cy, &cz, &r) != 0)
        return luaL_error(L, "bounding_sphere failed: %s", alea_error());
    lua_pushnumber(L, cx);
    lua_pushnumber(L, cy);
    lua_pushnumber(L, cz);
    lua_pushnumber(L, r);
    return 4;
}

/* sys:estimate_volumes(n_rays, cx, cy, cz, radius) -> table of volumes */
static int l_estimate_volumes(lua_State* L) {
    alea_system_t* sys = alea_get_sys(L, 1);
    int n_rays = (int)luaL_checkinteger(L, 2);
    double cx  = luaL_checknumber(L, 3);
    double cy  = luaL_checknumber(L, 4);
    double cz  = luaL_checknumber(L, 5);
    double r   = luaL_checknumber(L, 6);

    size_t nc = alea_cell_count(sys);
    double* volumes = (double*)calloc(nc, sizeof(double));
    double* errors  = (double*)calloc(nc, sizeof(double));
    if (!volumes || !errors) {
        free(volumes);
        free(errors);
        return luaL_error(L, "out of memory");
    }

    int rc = alea_estimate_cell_volumes(sys, cx, cy, cz, r, n_rays, volumes, errors);
    if (rc != 0) {
        free(volumes);
        free(errors);
        return luaL_error(L, "estimate_volumes failed: %s", alea_error());
    }

    lua_createtable(L, (int)nc, 0);
    for (size_t i = 0; i < nc; i++) {
        lua_createtable(L, 0, 2);
        lua_pushnumber(L, volumes[i]); lua_setfield(L, -2, "volume");
        lua_pushnumber(L, errors[i]);  lua_setfield(L, -2, "rel_error");
        lua_rawseti(L, -2, (lua_Integer)(i + 1));
    }

    free(volumes);
    free(errors);
    return 1;
}

/* sys:remove_cells_by_volume(volumes_table, threshold) -> count removed */
static int l_remove_cells_by_volume(lua_State* L) {
    alea_system_t* sys = alea_get_sys(L, 1);
    luaL_checktype(L, 2, LUA_TTABLE);
    double threshold = luaL_checknumber(L, 3);

    size_t nc = alea_cell_count(sys);
    double* volumes = (double*)calloc(nc, sizeof(double));
    if (!volumes)
        return luaL_error(L, "out of memory");

    for (size_t i = 0; i < nc; i++) {
        lua_rawgeti(L, 2, (lua_Integer)(i + 1));
        if (lua_istable(L, -1)) {
            lua_getfield(L, -1, "volume");
            volumes[i] = lua_tonumber(L, -1);
            lua_pop(L, 1);
        } else {
            volumes[i] = lua_tonumber(L, -1);
        }
        lua_pop(L, 1);
    }

    int removed = alea_remove_cells_by_volume(sys, volumes, threshold);
    free(volumes);
    if (removed < 0)
        return luaL_error(L, "remove_cells_by_volume failed: %s", alea_error());
    lua_pushinteger(L, removed);
    return 1;
}

/* ============================================================================
 * Bbox tightening
 * ============================================================================ */

/* sys:tighten_all_bboxes(tol) -> count */
static int l_tighten_all_bboxes(lua_State* L) {
    alea_system_t* sys = alea_get_sys(L, 1);
    double tol = luaL_optnumber(L, 2, 1.0);
    int n = alea_tighten_all_bboxes(sys, tol);
    if (n < 0)
        return luaL_error(L, "tighten_all_bboxes failed: %s", alea_error());
    lua_pushinteger(L, n);
    return 1;
}

/* sys:tighten_bbox_numerical(cell_index) */
static int l_tighten_bbox_numerical(lua_State* L) {
    alea_system_t* sys = alea_get_sys(L, 1);
    int idx = (int)luaL_checkinteger(L, 2);
    if (alea_tighten_cell_bbox_numerical(sys, idx) != 0)
        return luaL_error(L, "tighten_bbox_numerical failed: %s", alea_error());
    return 0;
}

/* ============================================================================
 * Renumbering
 * ============================================================================ */

static int l_renumber_cells(lua_State* L) {
    alea_system_t* sys = alea_get_sys(L, 1);
    int start = (int)luaL_checkinteger(L, 2);
    if (alea_renumber_cells(sys, start) != 0)
        return luaL_error(L, "renumber_cells failed: %s", alea_error());
    return 0;
}

static int l_renumber_surfaces(lua_State* L) {
    alea_system_t* sys = alea_get_sys(L, 1);
    int start = (int)luaL_checkinteger(L, 2);
    if (alea_renumber_surfaces(sys, start) != 0)
        return luaL_error(L, "renumber_surfaces failed: %s", alea_error());
    return 0;
}

static int l_offset_cell_ids(lua_State* L) {
    alea_system_t* sys = alea_get_sys(L, 1);
    int offset = (int)luaL_checkinteger(L, 2);
    if (alea_offset_cell_ids(sys, offset) != 0)
        return luaL_error(L, "offset_cell_ids failed: %s", alea_error());
    return 0;
}

static int l_offset_surface_ids(lua_State* L) {
    alea_system_t* sys = alea_get_sys(L, 1);
    int offset = (int)luaL_checkinteger(L, 2);
    if (alea_offset_surface_ids(sys, offset) != 0)
        return luaL_error(L, "offset_surface_ids failed: %s", alea_error());
    return 0;
}

static int l_offset_material_ids(lua_State* L) {
    alea_system_t* sys = alea_get_sys(L, 1);
    int offset = (int)luaL_checkinteger(L, 2);
    if (alea_offset_material_ids(sys, offset) != 0)
        return luaL_error(L, "offset_material_ids failed: %s", alea_error());
    return 0;
}

/* ============================================================================
 * Flatten all cells
 * ============================================================================ */

/* sys:flatten_all() -> stats table */
static int l_flatten_all(lua_State* L) {
    alea_system_t* sys = alea_get_sys(L, 1);
    alea_simplify_stats_t stats;
    memset(&stats, 0, sizeof(stats));
    alea_flatten_all_cells(sys, &stats);

    lua_createtable(L, 0, 14);
    lua_pushinteger(L, (lua_Integer)stats.nodes_before);            lua_setfield(L, -2, "nodes_before");
    lua_pushinteger(L, (lua_Integer)stats.nodes_after);             lua_setfield(L, -2, "nodes_after");
    lua_pushinteger(L, (lua_Integer)stats.complements_eliminated);  lua_setfield(L, -2, "complements_eliminated");
    lua_pushinteger(L, (lua_Integer)stats.double_negations);        lua_setfield(L, -2, "double_negations");
    lua_pushinteger(L, (lua_Integer)stats.idempotent_reductions);   lua_setfield(L, -2, "idempotent_reductions");
    lua_pushinteger(L, (lua_Integer)stats.absorption_reductions);   lua_setfield(L, -2, "absorption_reductions");
    lua_pushinteger(L, (lua_Integer)stats.subtrees_deduplicated);   lua_setfield(L, -2, "subtrees_deduplicated");
    lua_pushinteger(L, (lua_Integer)stats.cell_complements_expanded);lua_setfield(L, -2, "cell_complements_expanded");
    lua_pushinteger(L, (lua_Integer)stats.contradictions_found);    lua_setfield(L, -2, "contradictions_found");
    lua_pushinteger(L, (lua_Integer)stats.tautologies_found);       lua_setfield(L, -2, "tautologies_found");
    lua_pushinteger(L, (lua_Integer)stats.empty_cells_removed);     lua_setfield(L, -2, "empty_cells_removed");
    lua_pushinteger(L, (lua_Integer)stats.union_branches_absorbed); lua_setfield(L, -2, "union_branches_absorbed");
    lua_pushinteger(L, (lua_Integer)stats.union_common_factors);    lua_setfield(L, -2, "union_common_factors");
    lua_pushinteger(L, (lua_Integer)stats.union_branches_subsumed); lua_setfield(L, -2, "union_branches_subsumed");
    return 1;
}

/* ============================================================================
 * Volume estimation (instance-based)
 * ============================================================================ */

static void lua_push_volume_path(lua_State* L,
                                 const alea_volume_path_t* path,
                                 const double* volume,
                                 const double* rel_error) {
    lua_createtable(L, 0, 12);
    lua_pushinteger(L, (lua_Integer)(path->path_id + 1)); lua_setfield(L, -2, "path_id");
    lua_pushinteger(L, (lua_Integer)path->terminal_cell_index); lua_setfield(L, -2, "cell_index");
    lua_pushinteger(L, (lua_Integer)path->terminal_cell_id); lua_setfield(L, -2, "cell_id");
    lua_pushinteger(L, (lua_Integer)path->material_id); lua_setfield(L, -2, "material_id");
    lua_pushinteger(L, (lua_Integer)path->universe_id); lua_setfield(L, -2, "universe_id");
    lua_pushinteger(L, (lua_Integer)path->depth); lua_setfield(L, -2, "depth");
    lua_pushboolean(L, path->truncated != 0); lua_setfield(L, -2, "truncated");

    if (volume) {
        lua_pushnumber(L, *volume); lua_setfield(L, -2, "volume");
    }
    if (rel_error) {
        lua_pushnumber(L, *rel_error); lua_setfield(L, -2, "rel_error");
    }

    lua_createtable(L, path->ancestor_count, 0);
    for (uint8_t i = 0; i < path->ancestor_count; i++) {
        lua_createtable(L, 0, 2);
        lua_pushinteger(L, (lua_Integer)path->ancestor_cell_indices[i]);
        lua_setfield(L, -2, "cell_index");
        lua_pushinteger(L, (lua_Integer)path->ancestor_universe_ids[i]);
        lua_setfield(L, -2, "universe_id");
        lua_rawseti(L, -2, (lua_Integer)(i + 1));
    }
    lua_setfield(L, -2, "ancestors");

    lua_createtable(L, path->lattice_step_count, 0);
    for (uint8_t i = 0; i < path->lattice_step_count; i++) {
        const alea_volume_lattice_step_t* step = &path->lattice_steps[i];
        lua_createtable(L, 0, 6);
        lua_pushinteger(L, (lua_Integer)step->lattice_cell_index); lua_setfield(L, -2, "cell_index");
        lua_pushinteger(L, (lua_Integer)step->fill_universe); lua_setfield(L, -2, "fill_universe");
        lua_pushinteger(L, (lua_Integer)step->i); lua_setfield(L, -2, "i");
        lua_pushinteger(L, (lua_Integer)step->j); lua_setfield(L, -2, "j");
        lua_pushinteger(L, (lua_Integer)step->k); lua_setfield(L, -2, "k");
        lua_pushinteger(L, (lua_Integer)step->linear_index); lua_setfield(L, -2, "linear_index");
        lua_rawseti(L, -2, (lua_Integer)(i + 1));
    }
    lua_setfield(L, -2, "lattice");
}

/* sys:volume_paths() -> table */
static int l_volume_paths(lua_State* L) {
    alea_system_t* sys = alea_get_sys(L, 1);
    size_t n = alea_volume_path_count(sys);

    alea_volume_path_t* paths = n ? (alea_volume_path_t*)calloc(n, sizeof(*paths)) : NULL;
    if (n && !paths) return luaL_error(L, "out of memory");

    size_t got = alea_volume_paths_get(sys, paths, n);
    if (got > n) got = n;

    lua_createtable(L, (int)got, 0);
    for (size_t i = 0; i < got; i++) {
        lua_push_volume_path(L, &paths[i], NULL, NULL);
        lua_rawseti(L, -2, (lua_Integer)(i + 1));
    }

    free(paths);
    return 1;
}

/* sys:estimate_path_volumes(n_rays) -> table of path records */
static int l_estimate_path_volumes(lua_State* L) {
    alea_system_t* sys = alea_get_sys(L, 1);
    int n_rays = (int)luaL_checkinteger(L, 2);
    size_t n = alea_volume_path_count(sys);

    alea_volume_path_t* paths = n ? (alea_volume_path_t*)calloc(n, sizeof(*paths)) : NULL;
    double* volumes = n ? (double*)calloc(n, sizeof(double)) : NULL;
    double* errors = n ? (double*)calloc(n, sizeof(double)) : NULL;
    if (n && (!paths || !volumes || !errors)) {
        free(paths); free(volumes); free(errors);
        return luaL_error(L, "out of memory");
    }

    size_t got = alea_volume_paths_get(sys, paths, n);
    if (got > n) got = n;

    if (alea_estimate_path_volumes(sys, n_rays, volumes, errors) != 0) {
        free(paths); free(volumes); free(errors);
        return luaL_error(L, "estimate_path_volumes failed: %s", alea_error());
    }

    lua_createtable(L, (int)got, 0);
    for (size_t i = 0; i < got; i++) {
        lua_push_volume_path(L, &paths[i], &volumes[i], &errors[i]);
        lua_rawseti(L, -2, (lua_Integer)(i + 1));
    }

    free(paths);
    free(volumes);
    free(errors);
    return 1;
}

/* ============================================================================
 * Tighten single cell bbox
 * ============================================================================ */

/* sys:tighten_cell_bbox(idx, tol) -> bbox table */
static int l_tighten_cell_bbox(lua_State* L) {
    alea_system_t* sys = alea_get_sys(L, 1);
    int idx = (int)luaL_checkinteger(L, 2);
    double tol = luaL_optnumber(L, 3, 1.0);
    alea_bbox_t out;
    if (alea_tighten_cell_bbox(sys, (size_t)idx, tol, &out) != 0)
        return luaL_error(L, "tighten_cell_bbox failed: %s", alea_error());

    lua_createtable(L, 0, 6);
    lua_pushnumber(L, out.min_x); lua_setfield(L, -2, "min_x");
    lua_pushnumber(L, out.max_x); lua_setfield(L, -2, "max_x");
    lua_pushnumber(L, out.min_y); lua_setfield(L, -2, "min_y");
    lua_pushnumber(L, out.max_y); lua_setfield(L, -2, "max_y");
    lua_pushnumber(L, out.min_z); lua_setfield(L, -2, "min_z");
    lua_pushnumber(L, out.max_z); lua_setfield(L, -2, "max_z");
    return 1;
}

/* ============================================================================
 * Mixture creation
 * ============================================================================ */

/* sys:create_mixture(ids, fracs, new_id) -> id */
static int l_create_mixture(lua_State* L) {
    alea_system_t* sys = alea_get_sys(L, 1);
    luaL_checktype(L, 2, LUA_TTABLE);
    luaL_checktype(L, 3, LUA_TTABLE);
    int new_id = (int)luaL_optinteger(L, 4, 0);

    lua_Integer n = luaL_len(L, 2);
    int* ids = (int*)calloc((size_t)n, sizeof(int));
    double* fracs = (double*)calloc((size_t)n, sizeof(double));
    if (!ids || !fracs) {
        free(ids); free(fracs);
        return luaL_error(L, "out of memory");
    }

    for (lua_Integer i = 1; i <= n; i++) {
        lua_rawgeti(L, 2, i);
        ids[i - 1] = (int)lua_tointeger(L, -1);
        lua_pop(L, 1);
        lua_rawgeti(L, 3, i);
        fracs[i - 1] = lua_tonumber(L, -1);
        lua_pop(L, 1);
    }

    int result = alea_create_mixture(sys, ids, fracs, (size_t)n, new_id);
    free(ids);
    free(fracs);
    if (result < 0)
        return luaL_error(L, "create_mixture failed: %s", alea_error());
    lua_pushinteger(L, result);
    return 1;
}

/* ============================================================================
 * Extract / filter
 * ============================================================================ */

/* sys:extract_region(bbox) -> new System */
static int l_extract_region(lua_State* L) {
    alea_system_t* sys = alea_get_sys(L, 1);
    luaL_checktype(L, 2, LUA_TTABLE);

    alea_bbox_t bbox;
    lua_getfield(L, 2, "min_x"); bbox.min_x = luaL_checknumber(L, -1); lua_pop(L, 1);
    lua_getfield(L, 2, "max_x"); bbox.max_x = luaL_checknumber(L, -1); lua_pop(L, 1);
    lua_getfield(L, 2, "min_y"); bbox.min_y = luaL_checknumber(L, -1); lua_pop(L, 1);
    lua_getfield(L, 2, "max_y"); bbox.max_y = luaL_checknumber(L, -1); lua_pop(L, 1);
    lua_getfield(L, 2, "min_z"); bbox.min_z = luaL_checknumber(L, -1); lua_pop(L, 1);
    lua_getfield(L, 2, "max_z"); bbox.max_z = luaL_checknumber(L, -1); lua_pop(L, 1);

    alea_lua_system_t* ud = (alea_lua_system_t*)lua_newuserdata(L, sizeof(alea_lua_system_t));
    alea_lua_system_init(ud);
    luaL_setmetatable(L, ALEA_SYSTEM_MT);

    ud->sys = alea_extract_region(sys, &bbox);
    if (!ud->sys)
        return luaL_error(L, "extract_region failed: %s", alea_error());
    return 1;
}

/* sys:cells_in_bbox(bbox) -> table of indices */
static int l_cells_in_bbox(lua_State* L) {
    alea_system_t* sys = alea_get_sys(L, 1);
    luaL_checktype(L, 2, LUA_TTABLE);

    alea_bbox_t bbox;
    lua_getfield(L, 2, "min_x"); bbox.min_x = luaL_checknumber(L, -1); lua_pop(L, 1);
    lua_getfield(L, 2, "max_x"); bbox.max_x = luaL_checknumber(L, -1); lua_pop(L, 1);
    lua_getfield(L, 2, "min_y"); bbox.min_y = luaL_checknumber(L, -1); lua_pop(L, 1);
    lua_getfield(L, 2, "max_y"); bbox.max_y = luaL_checknumber(L, -1); lua_pop(L, 1);
    lua_getfield(L, 2, "min_z"); bbox.min_z = luaL_checknumber(L, -1); lua_pop(L, 1);
    lua_getfield(L, 2, "max_z"); bbox.max_z = luaL_checknumber(L, -1); lua_pop(L, 1);

    int buf[4096];
    size_t n = alea_get_cells_in_bbox(sys, &bbox, buf, 4096);

    lua_createtable(L, (int)n, 0);
    for (size_t i = 0; i < n; i++) {
        lua_pushinteger(L, buf[i]);
        lua_rawseti(L, -2, (lua_Integer)(i + 1));
    }
    return 1;
}

/* ============================================================================
 * Node-level macrobody expansion
 * ============================================================================ */

/* sys:expand_macrobody_node(node) -> Node */
static int l_expand_macrobody_node(lua_State* L) {
    alea_system_t* sys = alea_get_sys(L, 1);
    alea_lua_node_t* nd = alea_check_node(L, 2);
    alea_node_id_t result = alea_expand_macrobody(sys, nd->id);
    alea_push_node(L, sys, result);
    return 1;
}

/* sys:expand_all_macrobodies_node(node) -> Node */
static int l_expand_all_macrobodies_node(lua_State* L) {
    alea_system_t* sys = alea_get_sys(L, 1);
    alea_lua_node_t* nd = alea_check_node(L, 2);
    alea_node_id_t result = alea_expand_all_macrobodies(sys, nd->id);
    alea_push_node(L, sys, result);
    return 1;
}

/* ============================================================================
 * Macrobody expansion (system-level, already present)
 * ============================================================================ */

static int l_expand_macrobodies(lua_State* L) {
    alea_system_t* sys = alea_get_sys(L, 1);
    int n = alea_expand_macrobodies_in_system(sys);
    if (n < 0)
        return luaL_error(L, "expand_macrobodies failed: %s", alea_error());
    lua_pushinteger(L, n);
    return 1;
}

/* ============================================================================
 * Void generation
 * ============================================================================ */

typedef struct {
    void_result_t* ptr;
    alea_system_t* sys;
    alea_lua_system_t* owner;
} alea_lua_void_result_t;

static int l_void_generate(lua_State* L) {
    alea_lua_system_t* owner = alea_check_system(L, 1);
    alea_system_t* sys = alea_get_sys(L, 1);
    luaL_checktype(L, 2, LUA_TTABLE);

    alea_bbox_t bbox;
    lua_getfield(L, 2, "min_x"); bbox.min_x = luaL_checknumber(L, -1); lua_pop(L, 1);
    lua_getfield(L, 2, "max_x"); bbox.max_x = luaL_checknumber(L, -1); lua_pop(L, 1);
    lua_getfield(L, 2, "min_y"); bbox.min_y = luaL_checknumber(L, -1); lua_pop(L, 1);
    lua_getfield(L, 2, "max_y"); bbox.max_y = luaL_checknumber(L, -1); lua_pop(L, 1);
    lua_getfield(L, 2, "min_z"); bbox.min_z = luaL_checknumber(L, -1); lua_pop(L, 1);
    lua_getfield(L, 2, "max_z"); bbox.max_z = luaL_checknumber(L, -1); lua_pop(L, 1);

    alea_lua_void_result_t* ud = (alea_lua_void_result_t*)lua_newuserdata(
        L, sizeof(alea_lua_void_result_t));
    ud->ptr = NULL;
    ud->sys = sys;
    ud->owner = owner;
    luaL_setmetatable(L, ALEA_VOIDRESULT_MT);
    lua_pushvalue(L, 1);
    lua_setiuservalue(L, -2, 1);

    ud->ptr = alea_void_generate_in_bbox(sys, &bbox);
    if (!ud->ptr)
        return luaL_error(L, "void_generate failed: %s", alea_error());
    owner->active_void_results++;
    return 1;
}

static int l_void_gc(lua_State* L) {
    alea_lua_void_result_t* ud = (alea_lua_void_result_t*)luaL_checkudata(L, 1, ALEA_VOIDRESULT_MT);
    if (ud->ptr) {
        alea_void_free(ud->ptr);
        ud->ptr = NULL;
    }
    if (ud->owner && ud->owner->active_void_results > 0) {
        ud->owner->active_void_results--;
        alea_lua_system_release_if_pending(ud->owner);
        ud->owner = NULL;
    }
    return 0;
}

static int l_void_count(lua_State* L) {
    alea_lua_void_result_t* ud = (alea_lua_void_result_t*)luaL_checkudata(L, 1, ALEA_VOIDRESULT_MT);
    if (!ud->ptr) return luaL_error(L, "void result freed");
    lua_pushinteger(L, (lua_Integer)alea_void_count(ud->ptr));
    return 1;
}

static int l_void_add_cells(lua_State* L) {
    alea_lua_void_result_t* ud = (alea_lua_void_result_t*)luaL_checkudata(L, 1, ALEA_VOIDRESULT_MT);
    if (!ud->ptr) return luaL_error(L, "void result freed");
    int n = alea_void_add_cells(ud->sys, ud->ptr);
    if (n < 0)
        return luaL_error(L, "void_add_cells failed: %s", alea_error());
    lua_pushinteger(L, n);
    return 1;
}

static int l_void_merge(lua_State* L) {
    alea_lua_void_result_t* ud = (alea_lua_void_result_t*)luaL_checkudata(L, 1, ALEA_VOIDRESULT_MT);
    if (!ud->ptr) return luaL_error(L, "void result freed");
    if (alea_void_merge(ud->sys, ud->ptr) < 0)
        return luaL_error(L, "void_merge failed: %s", alea_error());
    return 0;
}

/* void_result:get(idx) -> bbox table (1-based) */
static int l_void_get(lua_State* L) {
    alea_lua_void_result_t* ud = (alea_lua_void_result_t*)luaL_checkudata(L, 1, ALEA_VOIDRESULT_MT);
    if (!ud->ptr) return luaL_error(L, "void result freed");
    int i = (int)luaL_checkinteger(L, 2);
    alea_bbox_t box;
    if (alea_void_get(ud->ptr, (size_t)(i - 1), &box) != 0)
        return luaL_error(L, "void:get: invalid index %d", i);

    lua_createtable(L, 0, 6);
    lua_pushnumber(L, box.min_x); lua_setfield(L, -2, "min_x");
    lua_pushnumber(L, box.max_x); lua_setfield(L, -2, "max_x");
    lua_pushnumber(L, box.min_y); lua_setfield(L, -2, "min_y");
    lua_pushnumber(L, box.max_y); lua_setfield(L, -2, "max_y");
    lua_pushnumber(L, box.min_z); lua_setfield(L, -2, "min_z");
    lua_pushnumber(L, box.max_z); lua_setfield(L, -2, "max_z");
    return 1;
}

/* void_result:to_node() -> Node */
static int l_void_to_node(lua_State* L) {
    alea_lua_void_result_t* ud = (alea_lua_void_result_t*)luaL_checkudata(L, 1, ALEA_VOIDRESULT_MT);
    if (!ud->ptr) return luaL_error(L, "void result freed");
    alea_node_id_t node = alea_void_to_node(ud->sys, ud->ptr);
    alea_push_node(L, ud->sys, node);
    return 1;
}

/* alea.is_macrobody(type_int) -> bool */
static int l_is_macrobody(lua_State* L) {
    int type = (int)luaL_checkinteger(L, 1);
    lua_pushboolean(L, alea_is_macrobody((alea_primitive_type_t)type));
    return 1;
}

/* ============================================================================
 * Registration
 * ============================================================================ */

static const luaL_Reg util_methods[] = {
    {"flatten",                    l_flatten},
    {"flatten_all",                l_flatten_all},
    {"split_union_cells",          l_split_union_cells},
    {"extract_universe",           l_extract_universe},
    {"merge",                      l_merge},
    {"clone",                      l_clone},
    {"bounding_sphere",            l_bounding_sphere},
    {"estimate_volumes",           l_estimate_volumes},
    {"volume_paths",               l_volume_paths},
    {"estimate_path_volumes",      l_estimate_path_volumes},
    {"remove_cells_by_volume",     l_remove_cells_by_volume},
    {"tighten_all_bboxes",         l_tighten_all_bboxes},
    {"tighten_cell_bbox",          l_tighten_cell_bbox},
    {"tighten_bbox_numerical",     l_tighten_bbox_numerical},
    {"renumber_cells",             l_renumber_cells},
    {"renumber_surfaces",          l_renumber_surfaces},
    {"offset_cell_ids",            l_offset_cell_ids},
    {"offset_surface_ids",         l_offset_surface_ids},
    {"offset_material_ids",        l_offset_material_ids},
    {"expand_macrobodies",         l_expand_macrobodies},
    {"expand_macrobody_node",      l_expand_macrobody_node},
    {"expand_all_macrobodies_node",l_expand_all_macrobodies_node},
    {"create_mixture",             l_create_mixture},
    {"extract_region",             l_extract_region},
    {"cells_in_bbox",              l_cells_in_bbox},
    {"void_generate",              l_void_generate},
    {NULL, NULL}
};

static const luaL_Reg void_meta[] = {
    {"__gc", l_void_gc},
    {NULL, NULL}
};

static const luaL_Reg void_methods[] = {
    {"count",     l_void_count},
    {"add_cells", l_void_add_cells},
    {"merge",     l_void_merge},
    {"get",       l_void_get},
    {"to_node",   l_void_to_node},
    {NULL, NULL}
};

int luaopen_alea_util(lua_State* L) {
    /* Add util methods to System metatable's __index */
    luaL_getmetatable(L, ALEA_SYSTEM_MT);
    lua_getfield(L, -1, "__index");
    luaL_setfuncs(L, util_methods, 0);
    lua_pop(L, 2);

    /* Create VoidResult metatable */
    luaL_newmetatable(L, ALEA_VOIDRESULT_MT);
    luaL_setfuncs(L, void_meta, 0);
    lua_newtable(L);
    luaL_setfuncs(L, void_methods, 0);
    lua_setfield(L, -2, "__index");
    lua_pop(L, 1);

    /* Module-level functions on alea table (top of stack) */
    lua_pushcfunction(L, l_is_macrobody);
    lua_setfield(L, -2, "is_macrobody");

    return 0;
}
