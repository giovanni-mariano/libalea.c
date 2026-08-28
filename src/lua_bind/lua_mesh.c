// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

#include "alea_lua.h"
#include "alea_mesh.h"
#include <stdlib.h>

/* ============================================================================
 * Userdata
 * ============================================================================ */

typedef struct {
    alea_mesh_result_t* ptr;
} alea_lua_mesh_result_t;

#define ALEA_ADAPTIVE_GRID_MT "alea.AdaptiveGridResult"
typedef struct {
    alea_adaptive_grid_result_t* ptr;
} alea_lua_adaptive_grid_t;

static alea_lua_mesh_result_t* check_meshresult(lua_State* L, int idx) {
    return (alea_lua_mesh_result_t*)luaL_checkudata(L, idx, ALEA_MESHRESULT_MT);
}

/* ============================================================================
 * Config helper: Lua table -> alea_mesh_config_t
 * ============================================================================ */

static void lua_to_mesh_config(lua_State* L, int idx, alea_mesh_config_t* cfg) {
    luaL_checktype(L, idx, LUA_TTABLE);

    lua_getfield(L, idx, "x_min");
    if (!lua_isnil(L, -1)) cfg->x_min = lua_tonumber(L, -1);
    lua_pop(L, 1);
    lua_getfield(L, idx, "x_max");
    if (!lua_isnil(L, -1)) cfg->x_max = lua_tonumber(L, -1);
    lua_pop(L, 1);
    lua_getfield(L, idx, "y_min");
    if (!lua_isnil(L, -1)) cfg->y_min = lua_tonumber(L, -1);
    lua_pop(L, 1);
    lua_getfield(L, idx, "y_max");
    if (!lua_isnil(L, -1)) cfg->y_max = lua_tonumber(L, -1);
    lua_pop(L, 1);
    lua_getfield(L, idx, "z_min");
    if (!lua_isnil(L, -1)) cfg->z_min = lua_tonumber(L, -1);
    lua_pop(L, 1);
    lua_getfield(L, idx, "z_max");
    if (!lua_isnil(L, -1)) cfg->z_max = lua_tonumber(L, -1);
    lua_pop(L, 1);

    lua_getfield(L, idx, "nx");
    if (!lua_isnil(L, -1)) cfg->nx = (int)lua_tointeger(L, -1);
    lua_pop(L, 1);
    lua_getfield(L, idx, "ny");
    if (!lua_isnil(L, -1)) cfg->ny = (int)lua_tointeger(L, -1);
    lua_pop(L, 1);
    lua_getfield(L, idx, "nz");
    if (!lua_isnil(L, -1)) cfg->nz = (int)lua_tointeger(L, -1);
    lua_pop(L, 1);

    lua_getfield(L, idx, "format");
    if (!lua_isnil(L, -1)) cfg->format = (alea_mesh_format_t)(int)lua_tointeger(L, -1);
    lua_pop(L, 1);

    lua_getfield(L, idx, "void_material_id");
    if (!lua_isnil(L, -1)) cfg->void_material_id = (int)lua_tointeger(L, -1);
    lua_pop(L, 1);

    lua_getfield(L, idx, "auto_pad");
    if (!lua_isnil(L, -1)) cfg->auto_pad = lua_tonumber(L, -1);
    lua_pop(L, 1);

    lua_getfield(L, idx, "bounds_mode");
    if (!lua_isnil(L, -1))
        cfg->bounds_mode = (alea_mesh_bounds_mode_t)(int)luaL_checkinteger(L, -1);
    lua_pop(L, 1);

    lua_getfield(L, idx, "sampling_mode");
    if (!lua_isnil(L, -1))
        cfg->sampling_mode = (alea_mesh_sampling_mode_t)(int)lua_tointeger(L, -1);
    lua_pop(L, 1);

    lua_getfield(L, idx, "subsamples_per_axis");
    if (!lua_isnil(L, -1)) cfg->subsamples_per_axis = (int)lua_tointeger(L, -1);
    lua_pop(L, 1);

    lua_getfield(L, idx, "mixed_threshold");
    if (!lua_isnil(L, -1)) cfg->mixed_threshold = lua_tonumber(L, -1);
    lua_pop(L, 1);

    lua_getfield(L, idx, "target_error");
    if (!lua_isnil(L, -1)) cfg->target_error = luaL_checknumber(L, -1);
    lua_pop(L, 1);
    lua_getfield(L, idx, "max_refine_depth");
    if (!lua_isnil(L, -1)) cfg->max_refine_depth = (int)luaL_checkinteger(L, -1);
    lua_pop(L, 1);
    lua_getfield(L, idx, "max_samples_per_voxel");
    if (!lua_isnil(L, -1)) cfg->max_samples_per_voxel = (uint32_t)luaL_checkinteger(L, -1);
    lua_pop(L, 1);
    lua_getfield(L, idx, "max_total_samples");
    if (!lua_isnil(L, -1)) cfg->max_total_samples = (uint64_t)luaL_checkinteger(L, -1);
    lua_pop(L, 1);
    lua_getfield(L, idx, "sampling_seed");
    if (!lua_isnil(L, -1)) cfg->sampling_seed = (uint64_t)luaL_checkinteger(L, -1);
    lua_pop(L, 1);
    lua_getfield(L, idx, "workers");
    if (!lua_isnil(L, -1)) cfg->workers = (int)luaL_checkinteger(L, -1);
    lua_pop(L, 1);

    lua_getfield(L, idx, "fields");
    if (!lua_isnil(L, -1)) cfg->fields = (uint32_t)luaL_checkinteger(L, -1);
    lua_pop(L, 1);
}

static alea_lua_adaptive_grid_t* check_adaptive_grid(lua_State* L, int idx) {
    return (alea_lua_adaptive_grid_t*)luaL_checkudata(L, idx,
                                                       ALEA_ADAPTIVE_GRID_MT);
}

static void lua_to_mesh_export_options(lua_State* L, int idx,
                                       alea_mesh_export_options_t* options) {
    if (!lua_istable(L, idx)) return;
    lua_getfield(L, idx, "export_fields");
    if (!lua_isnil(L, -1))
        options->fields = (uint32_t)luaL_checkinteger(L, -1);
    lua_pop(L, 1);
    lua_getfield(L, idx, "max_fraction_materials");
    if (!lua_isnil(L, -1))
        options->max_fraction_materials = (int)luaL_checkinteger(L, -1);
    lua_pop(L, 1);
}

/* ============================================================================
 * System methods
 * ============================================================================ */

/* sys:mesh_sample(config) -> MeshResult */
static int l_mesh_sample(lua_State* L) {
    alea_system_t* sys = alea_get_sys(L, 1);
    alea_mesh_config_t cfg;
    alea_mesh_config_init(&cfg);
    if (lua_istable(L, 2))
        lua_to_mesh_config(L, 2, &cfg);

    alea_lua_mesh_result_t* ud = (alea_lua_mesh_result_t*)lua_newuserdata(
        L, sizeof(alea_lua_mesh_result_t));
    ud->ptr = NULL;
    luaL_setmetatable(L, ALEA_MESHRESULT_MT);

    ud->ptr = alea_mesh_sample(sys, &cfg);
    if (!ud->ptr)
        return luaL_error(L, "mesh_sample failed: %s", alea_error());
    return 1;
}

/* sys:mesh_export(config, filename) */
static int l_mesh_export_system(lua_State* L) {
    alea_system_t* sys = alea_get_sys(L, 1);
    alea_mesh_config_t cfg;
    alea_mesh_config_init(&cfg);
    if (lua_istable(L, 2))
        lua_to_mesh_config(L, 2, &cfg);
    const char* filename = luaL_checkstring(L, 3);
    alea_mesh_export_options_t options;
    alea_mesh_export_options_init(&options);
    lua_to_mesh_export_options(L, 2, &options);
    alea_mesh_result_t* mesh = alea_mesh_sample(sys, &cfg);
    if (!mesh)
        return luaL_error(L, "mesh_export sampling failed: %s", alea_error());
    int rc = alea_mesh_export_ex(mesh, cfg.format, filename, &options);
    alea_mesh_result_free(mesh);
    if (rc != 0)
        return luaL_error(L, "mesh_export failed: %s", alea_error());
    return 0;
}

/* ============================================================================
 * MeshResult methods
 * ============================================================================ */

/* mesh:export(format, filename) */
static int l_meshresult_export(lua_State* L) {
    alea_lua_mesh_result_t* ud = check_meshresult(L, 1);
    if (!ud->ptr) return luaL_error(L, "mesh result freed");
    int fmt = (int)luaL_checkinteger(L, 2);
    const char* filename = luaL_checkstring(L, 3);
    alea_mesh_export_options_t options;
    alea_mesh_export_options_init(&options);
    if (lua_istable(L, 4)) lua_to_mesh_export_options(L, 4, &options);
    if (alea_mesh_export_ex(ud->ptr, (alea_mesh_format_t)fmt, filename,
                            &options) != 0)
        return luaL_error(L, "mesh:export failed: %s", alea_error());
    return 0;
}

/* mesh:info() -> table */
static int l_meshresult_info(lua_State* L) {
    alea_lua_mesh_result_t* ud = check_meshresult(L, 1);
    if (!ud->ptr) return luaL_error(L, "mesh result freed");
    alea_mesh_result_t* m = ud->ptr;

    lua_createtable(L, 0, 13);
    lua_pushinteger(L, m->nx);            lua_setfield(L, -2, "nx");
    lua_pushinteger(L, m->ny);            lua_setfield(L, -2, "ny");
    lua_pushinteger(L, m->nz);            lua_setfield(L, -2, "nz");
    lua_pushinteger(L, m->num_materials); lua_setfield(L, -2, "num_materials");
    lua_pushinteger(L, m->mixed_count);   lua_setfield(L, -2, "mixed_count");
    lua_pushinteger(L, (lua_Integer)m->fraction_count);
    lua_setfield(L, -2, "fraction_count");
    lua_pushinteger(L, m->sampling_mode); lua_setfield(L, -2, "sampling_mode");
    lua_pushinteger(L, (lua_Integer)m->sampling_seed);
    lua_setfield(L, -2, "sampling_seed");
    lua_pushnumber(L, m->target_error); lua_setfield(L, -2, "target_error");
    lua_pushinteger(L, m->bounds_source); lua_setfield(L, -2, "bounds_source");
    lua_pushnumber(L, m->bounds_padding); lua_setfield(L, -2, "bounds_padding");

    /* unique_materials */
    lua_createtable(L, m->num_materials, 0);
    for (int i = 0; i < m->num_materials; i++) {
        lua_pushinteger(L, m->unique_materials[i]);
        lua_rawseti(L, -2, i + 1);
    }
    lua_setfield(L, -2, "unique_materials");
    return 1;
}

static void push_fraction_entries(lua_State* L, const alea_mesh_result_t* m,
                                  size_t voxel_index) {
    if (!m->fraction_spans || !m->fractions) {
        lua_newtable(L);
        return;
    }

    alea_mesh_fraction_span_t span = m->fraction_spans[voxel_index];
    lua_createtable(L, (int)span.count, 0);
    for (uint32_t i = 0; i < span.count; i++) {
        const alea_mesh_material_fraction_t* f =
            &m->fractions[(size_t)span.offset + (size_t)i];
        lua_createtable(L, 0, 2);
        lua_pushinteger(L, f->material_id);
        lua_setfield(L, -2, "material_id");
        lua_pushnumber(L, f->fraction);
        lua_setfield(L, -2, "fraction");
        lua_pushnumber(L, f->fraction);
        lua_setfield(L, -2, "sampled_fraction");
        lua_rawseti(L, -2, (lua_Integer)i + 1);
    }
}

/* mesh:material_fractions([voxel_index]) -> fraction table */
static int l_meshresult_material_fractions(lua_State* L) {
    alea_lua_mesh_result_t* ud = check_meshresult(L, 1);
    if (!ud->ptr) return luaL_error(L, "mesh result freed");
    alea_mesh_result_t* m = ud->ptr;
    size_t total = (size_t)m->nx * (size_t)m->ny * (size_t)m->nz;

    if (!lua_isnoneornil(L, 2)) {
        lua_Integer idx = luaL_checkinteger(L, 2);
        if (idx < 1 || (size_t)idx > total)
            return luaL_error(L, "voxel index out of range");
        push_fraction_entries(L, m, (size_t)idx - 1);
        return 1;
    }

    lua_createtable(L, (int)total, 0);
    for (size_t i = 0; i < total; i++) {
        push_fraction_entries(L, m, i);
        lua_rawseti(L, -2, (lua_Integer)i + 1);
    }
    return 1;
}

/* mesh:material_ids() -> flat table */
static int l_meshresult_material_ids(lua_State* L) {
    alea_lua_mesh_result_t* ud = check_meshresult(L, 1);
    if (!ud->ptr) return luaL_error(L, "mesh result freed");
    alea_mesh_result_t* m = ud->ptr;
    size_t total = (size_t)m->nx * (size_t)m->ny * (size_t)m->nz;
    lua_createtable(L, (int)total, 0);
    for (size_t i = 0; i < total; i++) {
        lua_pushinteger(L, m->material_ids[i]);
        lua_rawseti(L, -2, (lua_Integer)(i + 1));
    }
    return 1;
}

/* mesh:cell_ids() -> flat table */
static int l_meshresult_cell_ids(lua_State* L) {
    alea_lua_mesh_result_t* ud = check_meshresult(L, 1);
    if (!ud->ptr) return luaL_error(L, "mesh result freed");
    alea_mesh_result_t* m = ud->ptr;
    size_t total = (size_t)m->nx * (size_t)m->ny * (size_t)m->nz;
    lua_createtable(L, (int)total, 0);
    for (size_t i = 0; i < total; i++) {
        lua_pushinteger(L, m->cell_ids[i]);
        lua_rawseti(L, -2, (lua_Integer)(i + 1));
    }
    return 1;
}

static int l_adaptive_grid_sample(lua_State* L) {
    alea_system_t* sys = alea_get_sys(L, 1);
    alea_adaptive_grid_config_t cfg;
    alea_adaptive_grid_config_init(&cfg);
    if (lua_istable(L, 2)) {
        lua_to_mesh_config(L, 2, &cfg.sampling);
        lua_getfield(L, 2, "max_grid_depth");
        if (!lua_isnil(L, -1)) cfg.max_grid_depth = (uint32_t)luaL_checkinteger(L, -1);
        lua_pop(L, 1);
        lua_getfield(L, 2, "max_cells");
        if (!lua_isnil(L, -1)) cfg.max_cells = (size_t)luaL_checkinteger(L, -1);
        lua_pop(L, 1);
        lua_getfield(L, 2, "refine_mixed");
        if (!lua_isnil(L, -1)) cfg.refine_mixed = lua_toboolean(L, -1);
        lua_pop(L, 1);
        lua_getfield(L, 2, "refine_high_error");
        if (!lua_isnil(L, -1)) cfg.refine_high_error = lua_toboolean(L, -1);
        lua_pop(L, 1);
    }
    alea_lua_adaptive_grid_t* ud = (alea_lua_adaptive_grid_t*)lua_newuserdata(
        L, sizeof(*ud));
    ud->ptr = alea_adaptive_grid_sample(sys, &cfg);
    if (!ud->ptr) return luaL_error(L, "adaptive_grid_sample failed: %s", alea_error());
    luaL_setmetatable(L, ALEA_ADAPTIVE_GRID_MT);
    return 1;
}

/* mesh:sample_counts() -> flat table */
static int l_meshresult_sample_counts(lua_State* L) {
    alea_lua_mesh_result_t* ud = check_meshresult(L, 1);
    if (!ud->ptr) return luaL_error(L, "mesh result freed");
    alea_mesh_result_t* m = ud->ptr;
    if (!m->sample_counts)
        return luaL_error(L, "sample counts were not retained by the field mask");
    size_t total = (size_t)m->nx * (size_t)m->ny * (size_t)m->nz;
    lua_createtable(L, (int)total, 0);
    for (size_t i = 0; i < total; i++) {
        lua_pushinteger(L, m->sample_counts[i]);
        lua_rawseti(L, -2, (lua_Integer)(i + 1));
    }
    return 1;
}

static int l_meshresult_estimated_errors(lua_State* L) {
    alea_lua_mesh_result_t* ud = check_meshresult(L, 1);
    if (!ud->ptr) return luaL_error(L, "mesh result freed");
    alea_mesh_result_t* m = ud->ptr;
    if (!m->estimated_errors)
        return luaL_error(L, "estimated errors were not retained by the field mask");
    size_t total = (size_t)m->nx * (size_t)m->ny * (size_t)m->nz;
    lua_createtable(L, (int)total, 0);
    for (size_t i = 0; i < total; i++) {
        lua_pushnumber(L, m->estimated_errors[i]);
        lua_rawseti(L, -2, (lua_Integer)(i + 1));
    }
    return 1;
}

static int l_meshresult_refinement_flags(lua_State* L) {
    alea_lua_mesh_result_t* ud = check_meshresult(L, 1);
    if (!ud->ptr) return luaL_error(L, "mesh result freed");
    alea_mesh_result_t* m = ud->ptr;
    if (!m->refinement_flags)
        return luaL_error(L, "refinement flags were not retained by the field mask");
    size_t total = (size_t)m->nx * (size_t)m->ny * (size_t)m->nz;
    lua_createtable(L, (int)total, 0);
    for (size_t i = 0; i < total; i++) {
        lua_pushinteger(L, m->refinement_flags[i]);
        lua_rawseti(L, -2, (lua_Integer)(i + 1));
    }
    return 1;
}

/* mesh:tie_flags() -> flat table */
static int l_meshresult_tie_flags(lua_State* L) {
    alea_lua_mesh_result_t* ud = check_meshresult(L, 1);
    if (!ud->ptr) return luaL_error(L, "mesh result freed");
    alea_mesh_result_t* m = ud->ptr;
    if (!m->tie_flags)
        return luaL_error(L, "tie flags were not retained by the field mask");
    size_t total = (size_t)m->nx * (size_t)m->ny * (size_t)m->nz;
    lua_createtable(L, (int)total, 0);
    for (size_t i = 0; i < total; i++) {
        lua_pushinteger(L, m->tie_flags[i]);
        lua_rawseti(L, -2, (lua_Integer)(i + 1));
    }
    return 1;
}

/* __gc */
static int l_meshresult_gc(lua_State* L) {
    alea_lua_mesh_result_t* ud = check_meshresult(L, 1);
    if (ud->ptr) {
        alea_mesh_result_free(ud->ptr);
        ud->ptr = NULL;
    }
    return 0;
}

static int l_adaptive_grid_info(lua_State* L) {
    alea_lua_adaptive_grid_t* ud = check_adaptive_grid(L, 1);
    if (!ud->ptr) return luaL_error(L, "adaptive grid result freed");
    lua_createtable(L, 0, 5);
    lua_pushinteger(L, (lua_Integer)ud->ptr->cell_count); lua_setfield(L, -2, "cell_count");
    lua_pushinteger(L, (lua_Integer)ud->ptr->leaf_count); lua_setfield(L, -2, "leaf_count");
    lua_pushinteger(L, (lua_Integer)ud->ptr->root_count); lua_setfield(L, -2, "root_count");
    lua_pushinteger(L, ud->ptr->max_level); lua_setfield(L, -2, "max_level");
    lua_pushboolean(L, ud->ptr->balanced); lua_setfield(L, -2, "balanced");
    return 1;
}

static int l_adaptive_grid_cells(lua_State* L) {
    alea_lua_adaptive_grid_t* ud = check_adaptive_grid(L, 1);
    if (!ud->ptr) return luaL_error(L, "adaptive grid result freed");
    lua_createtable(L, (int)ud->ptr->cell_count, 0);
    for (size_t i = 0; i < ud->ptr->cell_count; i++) {
        const alea_adaptive_grid_cell_t* c = &ud->ptr->cells[i];
        lua_createtable(L, 8, 16);
        lua_pushinteger(L, (lua_Integer)c->id); lua_setfield(L, -2, "id");
        lua_pushinteger(L, (lua_Integer)c->parent_id); lua_setfield(L, -2, "parent_id");
        lua_pushinteger(L, c->level); lua_setfield(L, -2, "level");
        lua_pushboolean(L, c->is_leaf); lua_setfield(L, -2, "is_leaf");
        lua_pushinteger(L, c->flags); lua_setfield(L, -2, "flags");
        lua_pushinteger(L, c->material_id); lua_setfield(L, -2, "material_id");
        lua_pushinteger(L, c->cell_id); lua_setfield(L, -2, "cell_id");
        lua_pushboolean(L, c->mixed); lua_setfield(L, -2, "mixed");
        lua_pushnumber(L, c->dominant_fraction); lua_setfield(L, -2, "dominant_fraction");
        lua_pushnumber(L, c->estimated_error); lua_setfield(L, -2, "estimated_error");
        lua_pushinteger(L, c->sample_count); lua_setfield(L, -2, "sample_count");
        lua_pushnumber(L, c->x_min); lua_setfield(L, -2, "x_min");
        lua_pushnumber(L, c->x_max); lua_setfield(L, -2, "x_max");
        lua_pushnumber(L, c->y_min); lua_setfield(L, -2, "y_min");
        lua_pushnumber(L, c->y_max); lua_setfield(L, -2, "y_max");
        lua_pushnumber(L, c->z_min); lua_setfield(L, -2, "z_min");
        lua_pushnumber(L, c->z_max); lua_setfield(L, -2, "z_max");
        lua_createtable(L, 8, 0);
        for (int child = 0; child < 8; child++) {
            lua_pushinteger(L, (lua_Integer)c->child_ids[child]);
            lua_rawseti(L, -2, child + 1);
        }
        lua_setfield(L, -2, "child_ids");
        lua_rawseti(L, -2, (lua_Integer)i + 1);
    }
    return 1;
}

static int l_adaptive_grid_export(lua_State* L) {
    alea_lua_adaptive_grid_t* ud = check_adaptive_grid(L, 1);
    if (!ud->ptr) return luaL_error(L, "adaptive grid result freed");
    int fmt = (int)luaL_checkinteger(L, 2);
    const char* filename = luaL_checkstring(L, 3);
    if (alea_adaptive_grid_export(ud->ptr, (alea_mesh_format_t)fmt, filename) != 0)
        return luaL_error(L, "adaptive grid export failed: %s", alea_error());
    return 0;
}

static int l_adaptive_grid_gc(lua_State* L) {
    alea_lua_adaptive_grid_t* ud = check_adaptive_grid(L, 1);
    if (ud->ptr) alea_adaptive_grid_result_free(ud->ptr);
    ud->ptr = NULL;
    return 0;
}

/* ============================================================================
 * Registration
 * ============================================================================ */

static const luaL_Reg mesh_system_methods[] = {
    {"mesh_sample",  l_mesh_sample},
    {"mesh_export",  l_mesh_export_system},
    {"adaptive_grid_sample", l_adaptive_grid_sample},
    {NULL, NULL}
};

static const luaL_Reg meshresult_meta[] = {
    {"__gc", l_meshresult_gc},
    {NULL, NULL}
};

static const luaL_Reg meshresult_methods[] = {
    {"export",             l_meshresult_export},
    {"info",               l_meshresult_info},
    {"material_ids",       l_meshresult_material_ids},
    {"cell_ids",           l_meshresult_cell_ids},
    {"sample_counts",      l_meshresult_sample_counts},
    {"tie_flags",          l_meshresult_tie_flags},
    {"estimated_errors",   l_meshresult_estimated_errors},
    {"refinement_flags",   l_meshresult_refinement_flags},
    {"material_fractions", l_meshresult_material_fractions},
    {NULL, NULL}
};

int luaopen_alea_mesh(lua_State* L) {
    /* Add system methods */
    luaL_getmetatable(L, ALEA_SYSTEM_MT);
    lua_getfield(L, -1, "__index");
    luaL_setfuncs(L, mesh_system_methods, 0);
    lua_pop(L, 2);

    /* Create MeshResult metatable */
    luaL_newmetatable(L, ALEA_MESHRESULT_MT);
    luaL_setfuncs(L, meshresult_meta, 0);
    lua_newtable(L);
    luaL_setfuncs(L, meshresult_methods, 0);
    lua_setfield(L, -2, "__index");
    lua_pop(L, 1);

    luaL_newmetatable(L, ALEA_ADAPTIVE_GRID_MT);
    lua_pushcfunction(L, l_adaptive_grid_gc);
    lua_setfield(L, -2, "__gc");
    lua_newtable(L);
    lua_pushcfunction(L, l_adaptive_grid_info); lua_setfield(L, -2, "info");
    lua_pushcfunction(L, l_adaptive_grid_cells); lua_setfield(L, -2, "cells");
    lua_pushcfunction(L, l_adaptive_grid_export); lua_setfield(L, -2, "export");
    lua_setfield(L, -2, "__index");
    lua_pop(L, 1);

    return 0;
}
