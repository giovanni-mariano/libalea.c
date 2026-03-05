// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

#include "alea_lua.h"

/* ============================================================================
 * Material queries
 * ============================================================================ */

/* sys:material_count() -> int */
static int l_material_count(lua_State* L) {
    alea_system_t* sys = alea_get_sys(L, 1);
    lua_pushinteger(L, (lua_Integer)alea_material_count(sys));
    return 1;
}

/* sys:material_get_id(mat_index) -> int */
static int l_material_get_id(lua_State* L) {
    alea_system_t* sys = alea_get_sys(L, 1);
    int idx = (int)luaL_checkinteger(L, 2);
    int id = alea_material_get_id(sys, idx);
    if (id < 0)
        return luaL_error(L, "material_get_id: invalid index %d", idx);
    lua_pushinteger(L, id);
    return 1;
}

/* sys:find_material(material_id) -> index or nil */
static int l_find_material(lua_State* L) {
    alea_system_t* sys = alea_get_sys(L, 1);
    int mat_id = (int)luaL_checkinteger(L, 2);
    int idx = alea_find_material_by_id(sys, mat_id);
    if (idx < 0)
        lua_pushnil(L);
    else
        lua_pushinteger(L, idx);
    return 1;
}

/* ============================================================================
 * Material modification
 * ============================================================================ */

/* sys:material_add_nuclide(mat_index, zaid, library, fraction) -> 0 */
static int l_material_add_nuclide(lua_State* L) {
    alea_system_t* sys = alea_get_sys(L, 1);
    int mat_idx = (int)luaL_checkinteger(L, 2);
    int zaid = (int)luaL_checkinteger(L, 3);
    const char* lib = lua_isnil(L, 4) ? NULL : luaL_checkstring(L, 4);
    double frac = luaL_checknumber(L, 5);
    if (alea_material_add_nuclide(sys, mat_idx, zaid, lib, frac) != 0)
        return luaL_error(L, "material_add_nuclide failed: %s", alea_error());
    return 0;
}

/* sys:material_add_element(mat_index, Z, library, fraction) -> 0 */
static int l_material_add_element(lua_State* L) {
    alea_system_t* sys = alea_get_sys(L, 1);
    int mat_idx = (int)luaL_checkinteger(L, 2);
    int Z = (int)luaL_checkinteger(L, 3);
    const char* lib = lua_isnil(L, 4) ? NULL : luaL_checkstring(L, 4);
    double frac = luaL_checknumber(L, 5);
    if (alea_material_add_element(sys, mat_idx, Z, lib, frac) != 0)
        return luaL_error(L, "material_add_element failed: %s", alea_error());
    return 0;
}

/* sys:material_set_density(mat_index, density) */
static int l_material_set_density(lua_State* L) {
    alea_system_t* sys = alea_get_sys(L, 1);
    int mat_idx = (int)luaL_checkinteger(L, 2);
    double density = luaL_checknumber(L, 3);
    if (alea_material_set_density(sys, mat_idx, density) != 0)
        return luaL_error(L, "material_set_density failed: %s", alea_error());
    return 0;
}

/* sys:material_set_weight_fraction(mat_index, is_weight) */
static int l_material_set_weight_fraction(lua_State* L) {
    alea_system_t* sys = alea_get_sys(L, 1);
    int mat_idx = (int)luaL_checkinteger(L, 2);
    luaL_checktype(L, 3, LUA_TBOOLEAN);
    bool is_weight = lua_toboolean(L, 3);
    if (alea_material_set_weight_fraction(sys, mat_idx, is_weight) != 0)
        return luaL_error(L, "material_set_weight_fraction failed: %s", alea_error());
    return 0;
}

/* sys:material_expand_elements(mat_index) */
static int l_material_expand_elements(lua_State* L) {
    alea_system_t* sys = alea_get_sys(L, 1);
    int mat_idx = (int)luaL_checkinteger(L, 2);
    if (alea_material_expand_elements(sys, mat_idx) != 0)
        return luaL_error(L, "material_expand_elements failed: %s", alea_error());
    return 0;
}

/* ============================================================================
 * Material nuclide/element inspection
 * ============================================================================ */

/* sys:material_nuclide_count(mat_index) -> int */
static int l_material_nuclide_count(lua_State* L) {
    alea_system_t* sys = alea_get_sys(L, 1);
    int mat_idx = (int)luaL_checkinteger(L, 2);
    lua_pushinteger(L, (lua_Integer)alea_material_nuclide_count(sys, mat_idx));
    return 1;
}

/* sys:material_nuclide_get(mat_index, nuc_index) -> {zaid, library, fraction} */
static int l_material_nuclide_get(lua_State* L) {
    alea_system_t* sys = alea_get_sys(L, 1);
    int mat_idx = (int)luaL_checkinteger(L, 2);
    int nuc_idx = (int)luaL_checkinteger(L, 3) - 1; /* 1-based Lua -> 0-based C */
    int zaid;
    const char* library;
    double fraction;
    if (alea_material_nuclide_get(sys, mat_idx, (size_t)nuc_idx,
                                   &zaid, &library, &fraction) != 0)
        return luaL_error(L, "material_nuclide_get: invalid indices mat=%d nuc=%d",
                          mat_idx, nuc_idx);
    lua_createtable(L, 0, 3);
    lua_pushinteger(L, zaid);              lua_setfield(L, -2, "zaid");
    if (library)
        lua_pushstring(L, library);
    else
        lua_pushnil(L);
    lua_setfield(L, -2, "library");
    lua_pushnumber(L, fraction);           lua_setfield(L, -2, "fraction");
    return 1;
}

/* sys:material_element_count(mat_index) -> int */
static int l_material_element_count(lua_State* L) {
    alea_system_t* sys = alea_get_sys(L, 1);
    int mat_idx = (int)luaL_checkinteger(L, 2);
    lua_pushinteger(L, (lua_Integer)alea_material_element_count(sys, mat_idx));
    return 1;
}

/* sys:material_element_get(mat_index, elem_index) -> {Z, library, fraction} */
static int l_material_element_get(lua_State* L) {
    alea_system_t* sys = alea_get_sys(L, 1);
    int mat_idx = (int)luaL_checkinteger(L, 2);
    int elem_idx = (int)luaL_checkinteger(L, 3) - 1; /* 1-based */
    int Z;
    const char* library;
    double fraction;
    if (alea_material_element_get(sys, mat_idx, (size_t)elem_idx,
                                   &Z, &library, &fraction) != 0)
        return luaL_error(L, "material_element_get: invalid indices mat=%d elem=%d",
                          mat_idx, elem_idx);
    lua_createtable(L, 0, 3);
    lua_pushinteger(L, Z);                 lua_setfield(L, -2, "Z");
    if (library)
        lua_pushstring(L, library);
    else
        lua_pushnil(L);
    lua_setfield(L, -2, "library");
    lua_pushnumber(L, fraction);           lua_setfield(L, -2, "fraction");
    return 1;
}

/* sys:material_get_density(mat_index) -> density, has_density */
static int l_material_get_density(lua_State* L) {
    alea_system_t* sys = alea_get_sys(L, 1);
    int mat_idx = (int)luaL_checkinteger(L, 2);
    double density;
    bool has_density;
    if (alea_material_get_density(sys, mat_idx, &density, &has_density) != 0)
        return luaL_error(L, "material_get_density: invalid index %d", mat_idx);
    lua_pushnumber(L, density);
    lua_pushboolean(L, has_density);
    return 2;
}

/* sys:material_is_weight_fraction(mat_index) -> bool */
static int l_material_is_weight_fraction(lua_State* L) {
    alea_system_t* sys = alea_get_sys(L, 1);
    int mat_idx = (int)luaL_checkinteger(L, 2);
    lua_pushboolean(L, alea_material_is_weight_fraction(sys, mat_idx));
    return 1;
}

/* ============================================================================
 * Mixture queries
 * ============================================================================ */

/* sys:mixture_count() -> int */
static int l_mixture_count(lua_State* L) {
    alea_system_t* sys = alea_get_sys(L, 1);
    lua_pushinteger(L, (lua_Integer)alea_mixture_count(sys));
    return 1;
}

/* sys:mixture_get_id(mix_index) -> int */
static int l_mixture_get_id(lua_State* L) {
    alea_system_t* sys = alea_get_sys(L, 1);
    int idx = (int)luaL_checkinteger(L, 2);
    int id = alea_mixture_get_id(sys, idx);
    if (id < 0)
        return luaL_error(L, "mixture_get_id: invalid index %d", idx);
    lua_pushinteger(L, id);
    return 1;
}

/* sys:find_mixture(mixture_id) -> index or nil */
static int l_find_mixture(lua_State* L) {
    alea_system_t* sys = alea_get_sys(L, 1);
    int mix_id = (int)luaL_checkinteger(L, 2);
    int idx = alea_find_mixture_by_id(sys, mix_id);
    if (idx < 0)
        lua_pushnil(L);
    else
        lua_pushinteger(L, idx);
    return 1;
}

/* sys:mixture_component_count(mix_index) -> int */
static int l_mixture_component_count(lua_State* L) {
    alea_system_t* sys = alea_get_sys(L, 1);
    int mix_idx = (int)luaL_checkinteger(L, 2);
    lua_pushinteger(L, (lua_Integer)alea_mixture_component_count(sys, mix_idx));
    return 1;
}

/* sys:mixture_component_get(mix_index, comp_index) -> {material_id, fraction} */
static int l_mixture_component_get(lua_State* L) {
    alea_system_t* sys = alea_get_sys(L, 1);
    int mix_idx = (int)luaL_checkinteger(L, 2);
    int comp_idx = (int)luaL_checkinteger(L, 3) - 1; /* 1-based */
    int material_id;
    double fraction;
    if (alea_mixture_component_get(sys, mix_idx, (size_t)comp_idx,
                                    &material_id, &fraction) != 0)
        return luaL_error(L, "mixture_component_get: invalid indices mix=%d comp=%d",
                          mix_idx, comp_idx);
    lua_createtable(L, 0, 2);
    lua_pushinteger(L, material_id);       lua_setfield(L, -2, "material_id");
    lua_pushnumber(L, fraction);           lua_setfield(L, -2, "fraction");
    return 1;
}

/* ============================================================================
 * Cell-mixture assignment
 * ============================================================================ */

/* sys:cell_set_mixture(cell_index, mixture_id) */
static int l_cell_set_mixture(lua_State* L) {
    alea_system_t* sys = alea_get_sys(L, 1);
    int cell_idx = (int)luaL_checkinteger(L, 2);
    int mix_id = (int)luaL_checkinteger(L, 3);
    if (alea_cell_set_mixture(sys, cell_idx, mix_id) != 0)
        return luaL_error(L, "cell_set_mixture failed: %s", alea_error());
    return 0;
}

/* ============================================================================
 * Registration
 * ============================================================================ */

static const luaL_Reg material_methods[] = {
    {"material_count",              l_material_count},
    {"material_get_id",             l_material_get_id},
    {"find_material",               l_find_material},
    {"material_add_nuclide",        l_material_add_nuclide},
    {"material_add_element",        l_material_add_element},
    {"material_set_density",        l_material_set_density},
    {"material_set_weight_fraction",l_material_set_weight_fraction},
    {"material_expand_elements",    l_material_expand_elements},
    {"material_nuclide_count",      l_material_nuclide_count},
    {"material_nuclide_get",        l_material_nuclide_get},
    {"material_element_count",      l_material_element_count},
    {"material_element_get",        l_material_element_get},
    {"material_get_density",        l_material_get_density},
    {"material_is_weight_fraction", l_material_is_weight_fraction},
    {"mixture_count",               l_mixture_count},
    {"mixture_get_id",              l_mixture_get_id},
    {"find_mixture",                l_find_mixture},
    {"mixture_component_count",     l_mixture_component_count},
    {"mixture_component_get",       l_mixture_component_get},
    {"cell_set_mixture",            l_cell_set_mixture},
    {NULL, NULL}
};

int luaopen_alea_materials(lua_State* L) {
    luaL_getmetatable(L, ALEA_SYSTEM_MT);
    lua_getfield(L, -1, "__index");
    luaL_setfuncs(L, material_methods, 0);
    lua_pop(L, 2);
    return 0;
}
