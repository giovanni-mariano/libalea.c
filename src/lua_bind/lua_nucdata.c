// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/**
 * @file lua_nucdata.c
 * @brief Lua bindings for the nuclear data module
 *
 * Provides three userdata types:
 *   alea.Xsdir     — xsdir directory (load nuclides)
 *   alea.Nuclide   — decoded ACE nuclide (cross-section queries)
 *   alea.NucMaterial — material composition (macroscopic XS)
 *
 * Example:
 *   local xsdir = alea.nuc_load_xsdir("/path/to/xsdir")
 *   local u235 = xsdir:load_nuclide("92235.80c")
 *   print(u235:xs_total(1.0))
 *
 *   local mat = alea.nuc_material()
 *   mat:add(u235, 0.048)
 *   print(mat:xs_total(1.0))
 *   print(mat:mean_free_path(1.0))
 */

#include "alea_lua.h"
#include "alea_nucdata.h"
#include "core/alea_system.h"
#include <stdlib.h>
#include <string.h>

#define ALEA_XSDIR_MT      "alea.Xsdir"
#define ALEA_NUCLIDE_MT    "alea.Nuclide"
#define ALEA_NUCMAT_MT     "alea.NucMaterial"
#define ALEA_MULTIGROUP_MT "alea.Multigroup"

/* Nuclide userdata — declared early because xsdir methods use it */
typedef struct {
    alea_nuc_nuclide_t* nuc;
    int owned;   /* 1 = caller must free, 0 = cached in xsdir */
} lua_nuc_ud_t;

/* ============================================================================
 * Xsdir userdata
 * ============================================================================ */

static int l_nuc_load_xsdir(lua_State* L) {
    const char* path = luaL_checkstring(L, 1);
    alea_nuc_xsdir_t* xsdir = alea_nuc_xsdir_load(path);
    if (!xsdir)
        return luaL_error(L, "failed to load xsdir: %s", path);

    alea_nuc_xsdir_t** ud = lua_newuserdata(L, sizeof(alea_nuc_xsdir_t*));
    *ud = xsdir;
    luaL_setmetatable(L, ALEA_XSDIR_MT);
    return 1;
}

static int l_nuc_load_xsdir_dir(lua_State* L) {
    const char* path = luaL_checkstring(L, 1);
    alea_nuc_xsdir_t* xsdir = alea_nuc_xsdir_load_dir(path);
    if (!xsdir)
        return luaL_error(L, "failed to load xsdir directory: %s", path);

    alea_nuc_xsdir_t** ud = lua_newuserdata(L, sizeof(alea_nuc_xsdir_t*));
    *ud = xsdir;
    luaL_setmetatable(L, ALEA_XSDIR_MT);
    return 1;
}

static alea_nuc_xsdir_t* check_xsdir(lua_State* L, int idx) {
    alea_nuc_xsdir_t** ud = luaL_checkudata(L, idx, ALEA_XSDIR_MT);
    if (!*ud) luaL_error(L, "xsdir has been destroyed");
    return *ud;
}

/* xsdir:count() -> int */
static int l_xsdir_count(lua_State* L) {
    alea_nuc_xsdir_t* xsdir = check_xsdir(L, 1);
    lua_pushinteger(L, (lua_Integer)alea_nuc_xsdir_count(xsdir));
    return 1;
}

/* xsdir:find(zaid) -> {zaid, awr, temperature, filename} or nil */
static int l_xsdir_find(lua_State* L) {
    alea_nuc_xsdir_t* xsdir = check_xsdir(L, 1);
    const char* zaid = luaL_checkstring(L, 2);
    const alea_nuc_xsdir_entry_t* e = alea_nuc_xsdir_find(xsdir, zaid);
    if (!e) {
        lua_pushnil(L);
        return 1;
    }
    lua_createtable(L, 0, 4);
    lua_pushstring(L, e->zaid);      lua_setfield(L, -2, "zaid");
    lua_pushnumber(L, e->awr);       lua_setfield(L, -2, "awr");
    lua_pushnumber(L, e->temperature); lua_setfield(L, -2, "temperature");
    lua_pushstring(L, e->filename);  lua_setfield(L, -2, "filename");
    return 1;
}

/* xsdir:load_nuclide(zaid) -> Nuclide (cached in xsdir)
 *
 * Returns a shared, cached pointer. Multiple calls with the same ZAID
 * return the same nuclide. Do NOT modify in place — use nuc:broadened(kT)
 * to get an independent Doppler-broadened copy. */
static int l_xsdir_load_nuclide(lua_State* L) {
    alea_nuc_xsdir_t* xsdir = check_xsdir(L, 1);
    const char* zaid = luaL_checkstring(L, 2);
    alea_nuc_nuclide_t* nuc = alea_nuc_xsdir_get_nuclide(xsdir, zaid);
    if (!nuc)
        return luaL_error(L, "failed to load nuclide: %s", zaid);

    lua_nuc_ud_t* ud = lua_newuserdata(L, sizeof(lua_nuc_ud_t));
    ud->nuc = nuc;
    ud->owned = 0;  /* cached in xsdir */
    luaL_setmetatable(L, ALEA_NUCLIDE_MT);

    /* Keep xsdir alive as long as this nuclide exists */
    lua_pushvalue(L, 1);
    lua_setuservalue(L, -2);

    return 1;
}

static int l_xsdir_gc(lua_State* L) {
    alea_nuc_xsdir_t** ud = luaL_checkudata(L, 1, ALEA_XSDIR_MT);
    if (*ud) { alea_nuc_xsdir_free(*ud); *ud = NULL; }
    return 0;
}

/* ============================================================================
 * Nuclide userdata
 * ============================================================================ */

static alea_nuc_nuclide_t* check_nuclide(lua_State* L, int idx) {
    lua_nuc_ud_t* ud = luaL_checkudata(L, idx, ALEA_NUCLIDE_MT);
    if (!ud->nuc) luaL_error(L, "nuclide has been destroyed");
    return ud->nuc;
}

/* nuc:zaid() -> string */
static int l_nuclide_zaid(lua_State* L) {
    alea_nuc_nuclide_t* nuc = check_nuclide(L, 1);
    lua_pushstring(L, nuc->zaid);
    return 1;
}

/* nuc:Z() -> int */
static int l_nuclide_Z(lua_State* L) {
    alea_nuc_nuclide_t* nuc = check_nuclide(L, 1);
    lua_pushinteger(L, nuc->Z);
    return 1;
}

/* nuc:A() -> int */
static int l_nuclide_A(lua_State* L) {
    alea_nuc_nuclide_t* nuc = check_nuclide(L, 1);
    lua_pushinteger(L, nuc->A);
    return 1;
}

/* nuc:awr() -> number */
static int l_nuclide_awr(lua_State* L) {
    alea_nuc_nuclide_t* nuc = check_nuclide(L, 1);
    lua_pushnumber(L, nuc->awr);
    return 1;
}

/* nuc:temperature() -> number (kT in MeV) */
static int l_nuclide_temperature(lua_State* L) {
    alea_nuc_nuclide_t* nuc = check_nuclide(L, 1);
    lua_pushnumber(L, nuc->temperature);
    return 1;
}

/* nuc:n_energies() -> int */
static int l_nuclide_n_energies(lua_State* L) {
    alea_nuc_nuclide_t* nuc = check_nuclide(L, 1);
    lua_pushinteger(L, nuc->n_energies);
    return 1;
}

/* nuc:n_reactions() -> int */
static int l_nuclide_n_reactions(lua_State* L) {
    alea_nuc_nuclide_t* nuc = check_nuclide(L, 1);
    lua_pushinteger(L, nuc->n_reactions);
    return 1;
}

/* nuc:xs_total(energy) -> number (barns) */
static int l_nuclide_xs_total(lua_State* L) {
    alea_nuc_nuclide_t* nuc = check_nuclide(L, 1);
    double E = luaL_checknumber(L, 2);
    lua_pushnumber(L, alea_nuc_xs_total(nuc, E));
    return 1;
}

/* nuc:xs_absorption(energy) -> number */
static int l_nuclide_xs_absorption(lua_State* L) {
    alea_nuc_nuclide_t* nuc = check_nuclide(L, 1);
    double E = luaL_checknumber(L, 2);
    lua_pushnumber(L, alea_nuc_xs_absorption(nuc, E));
    return 1;
}

/* nuc:xs_elastic(energy) -> number */
static int l_nuclide_xs_elastic(lua_State* L) {
    alea_nuc_nuclide_t* nuc = check_nuclide(L, 1);
    double E = luaL_checknumber(L, 2);
    lua_pushnumber(L, alea_nuc_xs_elastic(nuc, E));
    return 1;
}

/* nuc:xs_reaction(mt, energy) -> number */
static int l_nuclide_xs_reaction(lua_State* L) {
    alea_nuc_nuclide_t* nuc = check_nuclide(L, 1);
    int mt = (int)luaL_checkinteger(L, 2);
    double E = luaL_checknumber(L, 3);
    lua_pushnumber(L, alea_nuc_xs_reaction(nuc, mt, E));
    return 1;
}

/* nuc:xs_heating(energy) -> number (MeV·barn) */
static int l_nuclide_xs_heating(lua_State* L) {
    alea_nuc_nuclide_t* nuc = check_nuclide(L, 1);
    double E = luaL_checknumber(L, 2);
    lua_pushnumber(L, alea_nuc_xs_heating(nuc, E));
    return 1;
}

/* nuc:nu_bar(energy) -> number */
static int l_nuclide_nu_bar(lua_State* L) {
    alea_nuc_nuclide_t* nuc = check_nuclide(L, 1);
    double E = luaL_checknumber(L, 2);
    lua_pushnumber(L, alea_nuc_nu_bar(nuc, E));
    return 1;
}

/* nuc:broadened(kT_target) -> new Nuclide at higher temperature
 *
 * Returns a new, independent nuclide Doppler-broadened to kT_target.
 * The original nuclide is not modified. Requires the xsdir uservalue
 * to reload a fresh copy from disk before broadening. */
static int l_nuclide_broadened(lua_State* L) {
    lua_nuc_ud_t* src_ud = luaL_checkudata(L, 1, ALEA_NUCLIDE_MT);
    if (!src_ud->nuc) return luaL_error(L, "nuclide has been destroyed");
    double kT = luaL_checknumber(L, 2);

    /* Get the xsdir from uservalue to reload a fresh copy */
    lua_getuservalue(L, 1);
    alea_nuc_xsdir_t** xsdir_ud = luaL_testudata(L, -1, ALEA_XSDIR_MT);
    if (!xsdir_ud || !*xsdir_ud)
        return luaL_error(L, "broadened: nuclide has no associated xsdir");
    lua_pop(L, 1);

    /* Load a fresh (non-cached) copy */
    alea_nuc_nuclide_t* fresh = alea_nuc_load_nuclide(*xsdir_ud, src_ud->nuc->zaid);
    if (!fresh)
        return luaL_error(L, "broadened: failed to reload nuclide %s", src_ud->nuc->zaid);

    /* Broaden the fresh copy */
    alea_error_t err = alea_nuc_doppler_broaden(fresh, kT);
    if (err != ALEA_OK) {
        alea_nuc_nuclide_free(fresh);
        return luaL_error(L, "broadened: %s", alea_error_string(err));
    }

    /* Wrap as owned userdata */
    lua_nuc_ud_t* ud = lua_newuserdata(L, sizeof(lua_nuc_ud_t));
    ud->nuc = fresh;
    ud->owned = 1;
    luaL_setmetatable(L, ALEA_NUCLIDE_MT);

    return 1;
}

/* nuc:energy_range() -> E_min, E_max */
static int l_nuclide_energy_range(lua_State* L) {
    alea_nuc_nuclide_t* nuc = check_nuclide(L, 1);
    if (nuc->n_energies <= 0) return luaL_error(L, "no energy grid");
    lua_pushnumber(L, nuc->energy[0]);
    lua_pushnumber(L, nuc->energy[nuc->n_energies - 1]);
    return 2;
}

/* nuc:reactions() -> table of {mt, q_value} */
static int l_nuclide_reactions(lua_State* L) {
    alea_nuc_nuclide_t* nuc = check_nuclide(L, 1);
    lua_createtable(L, nuc->n_reactions, 0);
    for (int i = 0; i < nuc->n_reactions; i++) {
        lua_createtable(L, 0, 2);
        lua_pushinteger(L, nuc->reactions[i].mt);
        lua_setfield(L, -2, "mt");
        lua_pushnumber(L, nuc->reactions[i].q_value);
        lua_setfield(L, -2, "q_value");
        lua_rawseti(L, -2, i + 1);
    }
    return 1;
}

/* nuc:photon_xs(energy) -> {incoherent, coherent, photoelectric, pair} or nil */
static int l_nuclide_photon_xs(lua_State* L) {
    alea_nuc_nuclide_t* nuc = check_nuclide(L, 1);
    if (!nuc->photon) { lua_pushnil(L); return 1; }
    double E = luaL_checknumber(L, 2);
    lua_createtable(L, 0, 4);
    lua_pushnumber(L, alea_nuc_photon_xs_incoherent(nuc, E));
    lua_setfield(L, -2, "incoherent");
    lua_pushnumber(L, alea_nuc_photon_xs_coherent(nuc, E));
    lua_setfield(L, -2, "coherent");
    lua_pushnumber(L, alea_nuc_photon_xs_photoelectric(nuc, E));
    lua_setfield(L, -2, "photoelectric");
    lua_pushnumber(L, alea_nuc_photon_xs_pair(nuc, E));
    lua_setfield(L, -2, "pair");
    return 1;
}

static int l_nuclide_gc(lua_State* L) {
    lua_nuc_ud_t* ud = luaL_checkudata(L, 1, ALEA_NUCLIDE_MT);
    if (ud->nuc && ud->owned) {
        alea_nuc_nuclide_free(ud->nuc);
    }
    ud->nuc = NULL;
    return 0;
}

static int l_nuclide_tostring(lua_State* L) {
    alea_nuc_nuclide_t* nuc = check_nuclide(L, 1);
    lua_pushfstring(L, "Nuclide(%s, Z=%d, A=%d, %d energies, %d reactions)",
                    nuc->zaid, nuc->Z, nuc->A,
                    nuc->n_energies, nuc->n_reactions);
    return 1;
}

/* ============================================================================
 * NucMaterial userdata
 * ============================================================================ */

static int l_nuc_material_create(lua_State* L) {
    alea_nuc_material_t* mat = alea_nuc_material_create();
    if (!mat) return luaL_error(L, "failed to create material");

    alea_nuc_material_t** ud = lua_newuserdata(L, sizeof(alea_nuc_material_t*));
    *ud = mat;
    luaL_setmetatable(L, ALEA_NUCMAT_MT);
    return 1;
}

static alea_nuc_material_t* check_nucmat(lua_State* L, int idx) {
    alea_nuc_material_t** ud = luaL_checkudata(L, idx, ALEA_NUCMAT_MT);
    if (!*ud) luaL_error(L, "material has been destroyed");
    return *ud;
}

/* alea.nuc_material_from_cell(system, cell_index, xsdir) -> NucMaterial
 *
 * Builds a NucMaterial from a core cell's material definition.
 * Loads ACE nuclides via the xsdir and computes number densities
 * from the cell's composition and density.
 */
static int l_nuc_material_from_cell(lua_State* L) {
    alea_system_t* sys = alea_get_sys(L, 1);
    int cell_index = (int)luaL_checkinteger(L, 2) - 1; /* 1-based Lua */
    alea_nuc_xsdir_t* xsdir = check_xsdir(L, 3);

    /* Validate cell index */
    if (cell_index < 0 || (size_t)cell_index >= alea_cell_count(sys))
        return luaL_error(L, "cell index %d out of range [1,%d]",
                          cell_index + 1, (int)alea_cell_count(sys));

    const alea_cell_entry_t* cell = &sys->cells.data[cell_index];
    if (cell->material_index < 0)
        return luaL_error(L, "cell %d is void (no material)", cell_index + 1);

    alea_material_t* mat = &sys->materials.data[cell->material_index];

    /* Expand elements to nuclides if needed */
    if (mat->elements.count > 0 && mat->nuclides.count == 0)
        alea_mat_expand_elements(mat);

    if (mat->nuclides.count == 0)
        return luaL_error(L, "material %d has no nuclides", mat->material_id);

    /* Determine density: use cell density, fall back to material standard */
    double density = cell->density;
    bool is_mass = cell->is_mass_density;
    if (density == 0.0 && mat->has_standard_density) {
        density = mat->standard_density;
        is_mass = (density > 0.0);
        if (density < 0.0) density = -density;
    }
    if (density == 0.0)
        return luaL_error(L, "no density set for cell %d / material %d",
                          cell_index + 1, mat->material_id);

    alea_nuc_material_t* nmat = alea_nuc_material_from_core(mat, density, is_mass, xsdir);
    if (!nmat)
        return luaL_error(L, "failed to build nuclear material from cell %d", cell_index + 1);

    alea_nuc_material_t** ud = lua_newuserdata(L, sizeof(alea_nuc_material_t*));
    *ud = nmat;
    luaL_setmetatable(L, ALEA_NUCMAT_MT);
    return 1;
}

/* mat:add(nuclide, number_density) */
static int l_nucmat_add(lua_State* L) {
    alea_nuc_material_t* mat = check_nucmat(L, 1);
    alea_nuc_nuclide_t* nuc = check_nuclide(L, 2);
    double nd = luaL_checknumber(L, 3);
    alea_error_t err = alea_nuc_material_add(mat, nuc, nd);
    if (err != ALEA_OK)
        return luaL_error(L, "material_add failed: %s", alea_error_string(err));
    return 0;
}

/* mat:xs_total(energy) -> number (cm⁻¹) */
static int l_nucmat_xs_total(lua_State* L) {
    alea_nuc_material_t* mat = check_nucmat(L, 1);
    double E = luaL_checknumber(L, 2);
    lua_pushnumber(L, alea_nuc_mat_xs_total(mat, E));
    return 1;
}

/* mat:xs_absorption(energy) -> number */
static int l_nucmat_xs_absorption(lua_State* L) {
    alea_nuc_material_t* mat = check_nucmat(L, 1);
    double E = luaL_checknumber(L, 2);
    lua_pushnumber(L, alea_nuc_mat_xs_absorption(mat, E));
    return 1;
}

/* mat:xs_elastic(energy) -> number */
static int l_nucmat_xs_elastic(lua_State* L) {
    alea_nuc_material_t* mat = check_nucmat(L, 1);
    double E = luaL_checknumber(L, 2);
    lua_pushnumber(L, alea_nuc_mat_xs_elastic(mat, E));
    return 1;
}

/* mat:mean_free_path(energy) -> number (cm) */
static int l_nucmat_mean_free_path(lua_State* L) {
    alea_nuc_material_t* mat = check_nucmat(L, 1);
    double E = luaL_checknumber(L, 2);
    lua_pushnumber(L, alea_nuc_mean_free_path(mat, E));
    return 1;
}

/* mat:sample_distance(energy, xi) -> number (cm) */
static int l_nucmat_sample_distance(lua_State* L) {
    alea_nuc_material_t* mat = check_nucmat(L, 1);
    double E = luaL_checknumber(L, 2);
    double xi = luaL_checknumber(L, 3);
    lua_pushnumber(L, alea_nuc_sample_distance(mat, E, xi));
    return 1;
}

/* mat:n_components() -> int */
static int l_nucmat_n_components(lua_State* L) {
    alea_nuc_material_t* mat = check_nucmat(L, 1);
    lua_pushinteger(L, mat->n_components);
    return 1;
}

static int l_nucmat_gc(lua_State* L) {
    alea_nuc_material_t** ud = luaL_checkudata(L, 1, ALEA_NUCMAT_MT);
    if (*ud) { alea_nuc_material_destroy(*ud); *ud = NULL; }
    return 0;
}

/* ============================================================================
 * Multigroup userdata
 * ============================================================================ */

/* alea.nuc_mg_create(n_groups, bounds_table) -> Multigroup */
static int l_nuc_mg_create(lua_State* L) {
    int n_groups = (int)luaL_checkinteger(L, 1);
    luaL_checktype(L, 2, LUA_TTABLE);

    int n_bounds = (int)lua_rawlen(L, 2);
    if (n_bounds != n_groups + 1)
        return luaL_error(L, "expected %d group boundaries, got %d",
                          n_groups + 1, n_bounds);

    double* bounds = malloc((size_t)(n_groups + 1) * sizeof(double));
    if (!bounds) return luaL_error(L, "out of memory");

    for (int i = 0; i < n_groups + 1; i++) {
        lua_rawgeti(L, 2, i + 1);
        bounds[i] = luaL_checknumber(L, -1);
        lua_pop(L, 1);
    }

    alea_nuc_multigroup_t* mg = alea_nuc_mg_create(n_groups, bounds);
    free(bounds);
    if (!mg) return luaL_error(L, "failed to create multigroup structure");

    alea_nuc_multigroup_t** ud = lua_newuserdata(L, sizeof(alea_nuc_multigroup_t*));
    *ud = mg;
    luaL_setmetatable(L, ALEA_MULTIGROUP_MT);
    return 1;
}

static alea_nuc_multigroup_t* check_mg(lua_State* L, int idx) {
    alea_nuc_multigroup_t** ud = luaL_checkudata(L, idx, ALEA_MULTIGROUP_MT);
    if (!*ud) luaL_error(L, "multigroup has been destroyed");
    return *ud;
}

/* mg:collapse(nuclide) */
static int l_mg_collapse(lua_State* L) {
    alea_nuc_multigroup_t* mg = check_mg(L, 1);
    alea_nuc_nuclide_t* nuc = check_nuclide(L, 2);
    alea_error_t err = alea_nuc_mg_collapse(mg, nuc);
    if (err != ALEA_OK)
        return luaL_error(L, "mg_collapse failed: %s", alea_error_string(err));
    return 0;
}

/* mg:sigma_t(group) -> number */
static int l_mg_sigma_t(lua_State* L) {
    alea_nuc_multigroup_t* mg = check_mg(L, 1);
    int g = (int)luaL_checkinteger(L, 2) - 1; /* 1-based Lua */
    if (g < 0 || g >= mg->n_groups)
        return luaL_error(L, "group index %d out of range [1,%d]", g + 1, mg->n_groups);
    lua_pushnumber(L, mg->sigma_t[g]);
    return 1;
}

/* mg:sigma_a(group) -> number */
static int l_mg_sigma_a(lua_State* L) {
    alea_nuc_multigroup_t* mg = check_mg(L, 1);
    int g = (int)luaL_checkinteger(L, 2) - 1;
    if (g < 0 || g >= mg->n_groups)
        return luaL_error(L, "group index out of range");
    lua_pushnumber(L, mg->sigma_a[g]);
    return 1;
}

/* mg:sigma_s(group) -> number */
static int l_mg_sigma_s(lua_State* L) {
    alea_nuc_multigroup_t* mg = check_mg(L, 1);
    int g = (int)luaL_checkinteger(L, 2) - 1;
    if (g < 0 || g >= mg->n_groups)
        return luaL_error(L, "group index out of range");
    lua_pushnumber(L, mg->sigma_s[g]);
    return 1;
}

/* mg:sigma_f(group) -> number */
static int l_mg_sigma_f(lua_State* L) {
    alea_nuc_multigroup_t* mg = check_mg(L, 1);
    int g = (int)luaL_checkinteger(L, 2) - 1;
    if (g < 0 || g >= mg->n_groups)
        return luaL_error(L, "group index out of range");
    lua_pushnumber(L, mg->sigma_f[g]);
    return 1;
}

/* mg:nu_sigma_f(group) -> number */
static int l_mg_nu_sigma_f(lua_State* L) {
    alea_nuc_multigroup_t* mg = check_mg(L, 1);
    int g = (int)luaL_checkinteger(L, 2) - 1;
    if (g < 0 || g >= mg->n_groups)
        return luaL_error(L, "group index out of range");
    lua_pushnumber(L, mg->nu_sigma_f[g]);
    return 1;
}

/* mg:chi(group) -> number */
static int l_mg_chi(lua_State* L) {
    alea_nuc_multigroup_t* mg = check_mg(L, 1);
    int g = (int)luaL_checkinteger(L, 2) - 1;
    if (g < 0 || g >= mg->n_groups)
        return luaL_error(L, "group index out of range");
    lua_pushnumber(L, mg->chi[g]);
    return 1;
}

/* mg:scatter(g_from, g_to) -> number */
static int l_mg_scatter(lua_State* L) {
    alea_nuc_multigroup_t* mg = check_mg(L, 1);
    int gf = (int)luaL_checkinteger(L, 2) - 1;
    int gt = (int)luaL_checkinteger(L, 3) - 1;
    lua_pushnumber(L, alea_nuc_mg_scatter(mg, gf, gt));
    return 1;
}

/* mg:scatter_adjoint(g_from, g_to) -> number */
static int l_mg_scatter_adjoint(lua_State* L) {
    alea_nuc_multigroup_t* mg = check_mg(L, 1);
    int gf = (int)luaL_checkinteger(L, 2) - 1;
    int gt = (int)luaL_checkinteger(L, 3) - 1;
    lua_pushnumber(L, alea_nuc_mg_scatter_adjoint(mg, gf, gt));
    return 1;
}

/* mg:n_groups() -> int */
static int l_mg_n_groups(lua_State* L) {
    alea_nuc_multigroup_t* mg = check_mg(L, 1);
    lua_pushinteger(L, mg->n_groups);
    return 1;
}

static int l_mg_gc(lua_State* L) {
    alea_nuc_multigroup_t** ud = luaL_checkudata(L, 1, ALEA_MULTIGROUP_MT);
    if (*ud) { alea_nuc_mg_destroy(*ud); *ud = NULL; }
    return 0;
}

/* ============================================================================
 * Registration
 * ============================================================================ */

static const luaL_Reg xsdir_methods[] = {
    {"count",         l_xsdir_count},
    {"find",          l_xsdir_find},
    {"load_nuclide",  l_xsdir_load_nuclide},
    {NULL, NULL}
};

static const luaL_Reg xsdir_meta[] = {
    {"__gc", l_xsdir_gc},
    {NULL, NULL}
};

static const luaL_Reg nuclide_methods[] = {
    {"zaid",            l_nuclide_zaid},
    {"Z",               l_nuclide_Z},
    {"A",               l_nuclide_A},
    {"awr",             l_nuclide_awr},
    {"temperature",     l_nuclide_temperature},
    {"n_energies",      l_nuclide_n_energies},
    {"n_reactions",     l_nuclide_n_reactions},
    {"xs_total",        l_nuclide_xs_total},
    {"xs_absorption",   l_nuclide_xs_absorption},
    {"xs_elastic",      l_nuclide_xs_elastic},
    {"xs_reaction",     l_nuclide_xs_reaction},
    {"xs_heating",      l_nuclide_xs_heating},
    {"nu_bar",          l_nuclide_nu_bar},
    {"broadened",       l_nuclide_broadened},
    {"energy_range",    l_nuclide_energy_range},
    {"reactions",       l_nuclide_reactions},
    {"photon_xs",       l_nuclide_photon_xs},
    {NULL, NULL}
};

static const luaL_Reg nuclide_meta[] = {
    {"__gc",       l_nuclide_gc},
    {"__tostring", l_nuclide_tostring},
    {NULL, NULL}
};

static const luaL_Reg nucmat_methods[] = {
    {"add",              l_nucmat_add},
    {"xs_total",         l_nucmat_xs_total},
    {"xs_absorption",    l_nucmat_xs_absorption},
    {"xs_elastic",       l_nucmat_xs_elastic},
    {"mean_free_path",   l_nucmat_mean_free_path},
    {"sample_distance",  l_nucmat_sample_distance},
    {"n_components",     l_nucmat_n_components},
    {NULL, NULL}
};

static const luaL_Reg nucmat_meta[] = {
    {"__gc", l_nucmat_gc},
    {NULL, NULL}
};

static const luaL_Reg mg_methods[] = {
    {"collapse",         l_mg_collapse},
    {"sigma_t",          l_mg_sigma_t},
    {"sigma_a",          l_mg_sigma_a},
    {"sigma_s",          l_mg_sigma_s},
    {"sigma_f",          l_mg_sigma_f},
    {"nu_sigma_f",       l_mg_nu_sigma_f},
    {"chi",              l_mg_chi},
    {"scatter",          l_mg_scatter},
    {"scatter_adjoint",  l_mg_scatter_adjoint},
    {"n_groups",         l_mg_n_groups},
    {NULL, NULL}
};

static const luaL_Reg mg_meta[] = {
    {"__gc", l_mg_gc},
    {NULL, NULL}
};

static void register_type(lua_State* L, const char* mt_name,
                           const luaL_Reg* methods, const luaL_Reg* meta) {
    luaL_newmetatable(L, mt_name);
    if (meta) luaL_setfuncs(L, meta, 0);
    lua_newtable(L);
    luaL_setfuncs(L, methods, 0);
    lua_setfield(L, -2, "__index");
    lua_pop(L, 1);
}

int luaopen_alea_nucdata(lua_State* L) {
    /* Register userdata types */
    register_type(L, ALEA_XSDIR_MT, xsdir_methods, xsdir_meta);
    register_type(L, ALEA_NUCLIDE_MT, nuclide_methods, nuclide_meta);
    register_type(L, ALEA_NUCMAT_MT, nucmat_methods, nucmat_meta);
    register_type(L, ALEA_MULTIGROUP_MT, mg_methods, mg_meta);

    /* Add factory functions to the alea table (on top of stack) */
    lua_pushcfunction(L, l_nuc_load_xsdir);
    lua_setfield(L, -2, "nuc_load_xsdir");

    lua_pushcfunction(L, l_nuc_load_xsdir_dir);
    lua_setfield(L, -2, "nuc_load_xsdir_dir");

    lua_pushcfunction(L, l_nuc_material_create);
    lua_setfield(L, -2, "nuc_material");

    lua_pushcfunction(L, l_nuc_mg_create);
    lua_setfield(L, -2, "nuc_mg_create");

    lua_pushcfunction(L, l_nuc_material_from_cell);
    lua_setfield(L, -2, "nuc_material_from_cell");

    return 0;
}
