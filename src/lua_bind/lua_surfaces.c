// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

#include "alea_lua.h"

/* ============================================================================
 * Surface creation — all are System methods returning surface index (int)
 * ============================================================================ */

/* sys:plane(id, a, b, c, d) */
static int l_plane(lua_State* L) {
    alea_system_t* sys = alea_get_sys(L, 1);
    int id    = (int)luaL_checkinteger(L, 2);
    double a  = luaL_checknumber(L, 3);
    double b  = luaL_checknumber(L, 4);
    double c  = luaL_checknumber(L, 5);
    double d  = luaL_checknumber(L, 6);
    int idx = alea_plane_surface(sys, id, a, b, c, d);
    if (idx < 0) return luaL_error(L, "plane: %s", alea_error());
    lua_pushinteger(L, idx);
    return 1;
}

/* sys:sphere(id, cx, cy, cz, r) */
static int l_sphere(lua_State* L) {
    alea_system_t* sys = alea_get_sys(L, 1);
    int id    = (int)luaL_checkinteger(L, 2);
    double cx = luaL_checknumber(L, 3);
    double cy = luaL_checknumber(L, 4);
    double cz = luaL_checknumber(L, 5);
    double r  = luaL_checknumber(L, 6);
    int idx = alea_sphere_surface(sys, id, cx, cy, cz, r);
    if (idx < 0) return luaL_error(L, "sphere: %s", alea_error());
    lua_pushinteger(L, idx);
    return 1;
}

/* sys:cylinder_z(id, cx, cy, r) */
static int l_cylinder_z(lua_State* L) {
    alea_system_t* sys = alea_get_sys(L, 1);
    int id    = (int)luaL_checkinteger(L, 2);
    double cx = luaL_checknumber(L, 3);
    double cy = luaL_checknumber(L, 4);
    double r  = luaL_checknumber(L, 5);
    int idx = alea_cylinder_z_surface(sys, id, cx, cy, r);
    if (idx < 0) return luaL_error(L, "cylinder_z: %s", alea_error());
    lua_pushinteger(L, idx);
    return 1;
}

/* sys:cylinder_x(id, cy, cz, r) */
static int l_cylinder_x(lua_State* L) {
    alea_system_t* sys = alea_get_sys(L, 1);
    int id    = (int)luaL_checkinteger(L, 2);
    double cy = luaL_checknumber(L, 3);
    double cz = luaL_checknumber(L, 4);
    double r  = luaL_checknumber(L, 5);
    int idx = alea_cylinder_x_surface(sys, id, cy, cz, r);
    if (idx < 0) return luaL_error(L, "cylinder_x: %s", alea_error());
    lua_pushinteger(L, idx);
    return 1;
}

/* sys:cylinder_y(id, cx, cz, r) */
static int l_cylinder_y(lua_State* L) {
    alea_system_t* sys = alea_get_sys(L, 1);
    int id    = (int)luaL_checkinteger(L, 2);
    double cx = luaL_checknumber(L, 3);
    double cz = luaL_checknumber(L, 4);
    double r  = luaL_checknumber(L, 5);
    int idx = alea_cylinder_y_surface(sys, id, cx, cz, r);
    if (idx < 0) return luaL_error(L, "cylinder_y: %s", alea_error());
    lua_pushinteger(L, idx);
    return 1;
}

/* sys:cone_z(id, cx, cy, cz, t2) */
static int l_cone_z(lua_State* L) {
    alea_system_t* sys = alea_get_sys(L, 1);
    int id    = (int)luaL_checkinteger(L, 2);
    double cx = luaL_checknumber(L, 3);
    double cy = luaL_checknumber(L, 4);
    double cz = luaL_checknumber(L, 5);
    double t2 = luaL_checknumber(L, 6);
    int idx = alea_cone_z_surface(sys, id, cx, cy, cz, t2);
    if (idx < 0) return luaL_error(L, "cone_z: %s", alea_error());
    lua_pushinteger(L, idx);
    return 1;
}

/* sys:cone_x(id, cx, cy, cz, t2) */
static int l_cone_x(lua_State* L) {
    alea_system_t* sys = alea_get_sys(L, 1);
    int id    = (int)luaL_checkinteger(L, 2);
    double cx = luaL_checknumber(L, 3);
    double cy = luaL_checknumber(L, 4);
    double cz = luaL_checknumber(L, 5);
    double t2 = luaL_checknumber(L, 6);
    int idx = alea_cone_x_surface(sys, id, cx, cy, cz, t2);
    if (idx < 0) return luaL_error(L, "cone_x: %s", alea_error());
    lua_pushinteger(L, idx);
    return 1;
}

/* sys:cone_y(id, cx, cy, cz, t2) */
static int l_cone_y(lua_State* L) {
    alea_system_t* sys = alea_get_sys(L, 1);
    int id    = (int)luaL_checkinteger(L, 2);
    double cx = luaL_checknumber(L, 3);
    double cy = luaL_checknumber(L, 4);
    double cz = luaL_checknumber(L, 5);
    double t2 = luaL_checknumber(L, 6);
    int idx = alea_cone_y_surface(sys, id, cx, cy, cz, t2);
    if (idx < 0) return luaL_error(L, "cone_y: %s", alea_error());
    lua_pushinteger(L, idx);
    return 1;
}

/* sys:torus_z(id, cx, cy, cz, R, r) */
static int l_torus_z(lua_State* L) {
    alea_system_t* sys = alea_get_sys(L, 1);
    int id     = (int)luaL_checkinteger(L, 2);
    double cx  = luaL_checknumber(L, 3);
    double cy  = luaL_checknumber(L, 4);
    double cz  = luaL_checknumber(L, 5);
    double R   = luaL_checknumber(L, 6);
    double r   = luaL_checknumber(L, 7);
    int idx = alea_torus_z_surface(sys, id, cx, cy, cz, R, r);
    if (idx < 0) return luaL_error(L, "torus_z: %s", alea_error());
    lua_pushinteger(L, idx);
    return 1;
}

/* sys:torus_x(id, cx, cy, cz, R, r) */
static int l_torus_x(lua_State* L) {
    alea_system_t* sys = alea_get_sys(L, 1);
    int id     = (int)luaL_checkinteger(L, 2);
    double cx  = luaL_checknumber(L, 3);
    double cy  = luaL_checknumber(L, 4);
    double cz  = luaL_checknumber(L, 5);
    double R   = luaL_checknumber(L, 6);
    double r   = luaL_checknumber(L, 7);
    int idx = alea_torus_x_surface(sys, id, cx, cy, cz, R, r);
    if (idx < 0) return luaL_error(L, "torus_x: %s", alea_error());
    lua_pushinteger(L, idx);
    return 1;
}

/* sys:torus_y(id, cx, cy, cz, R, r) */
static int l_torus_y(lua_State* L) {
    alea_system_t* sys = alea_get_sys(L, 1);
    int id     = (int)luaL_checkinteger(L, 2);
    double cx  = luaL_checknumber(L, 3);
    double cy  = luaL_checknumber(L, 4);
    double cz  = luaL_checknumber(L, 5);
    double R   = luaL_checknumber(L, 6);
    double r   = luaL_checknumber(L, 7);
    int idx = alea_torus_y_surface(sys, id, cx, cy, cz, R, r);
    if (idx < 0) return luaL_error(L, "torus_y: %s", alea_error());
    lua_pushinteger(L, idx);
    return 1;
}

/* sys:box(id, xmin, xmax, ymin, ymax, zmin, zmax) */
static int l_box(lua_State* L) {
    alea_system_t* sys = alea_get_sys(L, 1);
    int id        = (int)luaL_checkinteger(L, 2);
    double xmin   = luaL_checknumber(L, 3);
    double xmax   = luaL_checknumber(L, 4);
    double ymin   = luaL_checknumber(L, 5);
    double ymax   = luaL_checknumber(L, 6);
    double zmin   = luaL_checknumber(L, 7);
    double zmax   = luaL_checknumber(L, 8);
    int idx = alea_box_surface(sys, id, xmin, xmax, ymin, ymax, zmin, zmax);
    if (idx < 0) return luaL_error(L, "box: %s", alea_error());
    lua_pushinteger(L, idx);
    return 1;
}

/* sys:quadric(id, A, B, C, D, E, F, G, H, I, J) */
static int l_quadric(lua_State* L) {
    alea_system_t* sys = alea_get_sys(L, 1);
    int id   = (int)luaL_checkinteger(L, 2);
    double A = luaL_checknumber(L, 3);
    double B = luaL_checknumber(L, 4);
    double C = luaL_checknumber(L, 5);
    double D = luaL_checknumber(L, 6);
    double E = luaL_checknumber(L, 7);
    double F = luaL_checknumber(L, 8);
    double G = luaL_checknumber(L, 9);
    double H = luaL_checknumber(L, 10);
    double I = luaL_checknumber(L, 11);
    double J = luaL_checknumber(L, 12);
    int idx = alea_quadric_surface(sys, id, A, B, C, D, E, F, G, H, I, J);
    if (idx < 0) return luaL_error(L, "quadric: %s", alea_error());
    lua_pushinteger(L, idx);
    return 1;
}

/* ============================================================================
 * Macrobodies
 * ============================================================================ */

/* sys:rcc(id, bx, by, bz, hx, hy, hz, r) */
static int l_rcc(lua_State* L) {
    alea_system_t* sys = alea_get_sys(L, 1);
    int id    = (int)luaL_checkinteger(L, 2);
    double bx = luaL_checknumber(L, 3);
    double by = luaL_checknumber(L, 4);
    double bz = luaL_checknumber(L, 5);
    double hx = luaL_checknumber(L, 6);
    double hy = luaL_checknumber(L, 7);
    double hz = luaL_checknumber(L, 8);
    double r  = luaL_checknumber(L, 9);
    int idx = alea_rcc_surface(sys, id, bx, by, bz, hx, hy, hz, r);
    if (idx < 0) return luaL_error(L, "rcc: %s", alea_error());
    lua_pushinteger(L, idx);
    return 1;
}

/* sys:sph(id, cx, cy, cz, r) */
static int l_sph(lua_State* L) {
    alea_system_t* sys = alea_get_sys(L, 1);
    int id    = (int)luaL_checkinteger(L, 2);
    double cx = luaL_checknumber(L, 3);
    double cy = luaL_checknumber(L, 4);
    double cz = luaL_checknumber(L, 5);
    double r  = luaL_checknumber(L, 6);
    int idx = alea_sph_surface(sys, id, cx, cy, cz, r);
    if (idx < 0) return luaL_error(L, "sph: %s", alea_error());
    lua_pushinteger(L, idx);
    return 1;
}

/* sys:trc(id, bx, by, bz, hx, hy, hz, r1, r2) */
static int l_trc(lua_State* L) {
    alea_system_t* sys = alea_get_sys(L, 1);
    int id    = (int)luaL_checkinteger(L, 2);
    double bx = luaL_checknumber(L, 3);
    double by = luaL_checknumber(L, 4);
    double bz = luaL_checknumber(L, 5);
    double hx = luaL_checknumber(L, 6);
    double hy = luaL_checknumber(L, 7);
    double hz = luaL_checknumber(L, 8);
    double r1 = luaL_checknumber(L, 9);
    double r2 = luaL_checknumber(L, 10);
    int idx = alea_trc_surface(sys, id, bx, by, bz, hx, hy, hz, r1, r2);
    if (idx < 0) return luaL_error(L, "trc: %s", alea_error());
    lua_pushinteger(L, idx);
    return 1;
}

/* sys:ell(id, v1x, v1y, v1z, v2x, v2y, v2z, major_len) */
static int l_ell(lua_State* L) {
    alea_system_t* sys = alea_get_sys(L, 1);
    int id     = (int)luaL_checkinteger(L, 2);
    double v1x = luaL_checknumber(L, 3);
    double v1y = luaL_checknumber(L, 4);
    double v1z = luaL_checknumber(L, 5);
    double v2x = luaL_checknumber(L, 6);
    double v2y = luaL_checknumber(L, 7);
    double v2z = luaL_checknumber(L, 8);
    double ml  = luaL_checknumber(L, 9);
    int idx = alea_ell_surface(sys, id, v1x, v1y, v1z, v2x, v2y, v2z, ml);
    if (idx < 0) return luaL_error(L, "ell: %s", alea_error());
    lua_pushinteger(L, idx);
    return 1;
}

/* sys:rec(id, bx, by, bz, hx, hy, hz, a1x, a1y, a1z, a2x, a2y, a2z) */
static int l_rec(lua_State* L) {
    alea_system_t* sys = alea_get_sys(L, 1);
    int id     = (int)luaL_checkinteger(L, 2);
    double bx  = luaL_checknumber(L, 3);
    double by  = luaL_checknumber(L, 4);
    double bz  = luaL_checknumber(L, 5);
    double hx  = luaL_checknumber(L, 6);
    double hy  = luaL_checknumber(L, 7);
    double hz  = luaL_checknumber(L, 8);
    double a1x = luaL_checknumber(L, 9);
    double a1y = luaL_checknumber(L, 10);
    double a1z = luaL_checknumber(L, 11);
    double a2x = luaL_checknumber(L, 12);
    double a2y = luaL_checknumber(L, 13);
    double a2z = luaL_checknumber(L, 14);
    int idx = alea_rec_surface(sys, id, bx, by, bz, hx, hy, hz,
                               a1x, a1y, a1z, a2x, a2y, a2z);
    if (idx < 0) return luaL_error(L, "rec: %s", alea_error());
    lua_pushinteger(L, idx);
    return 1;
}

/* sys:wed(id, vx, vy, vz, v1x, v1y, v1z, v2x, v2y, v2z, v3x, v3y, v3z) */
static int l_wed(lua_State* L) {
    alea_system_t* sys = alea_get_sys(L, 1);
    int id     = (int)luaL_checkinteger(L, 2);
    double vx  = luaL_checknumber(L, 3);
    double vy  = luaL_checknumber(L, 4);
    double vz  = luaL_checknumber(L, 5);
    double v1x = luaL_checknumber(L, 6);
    double v1y = luaL_checknumber(L, 7);
    double v1z = luaL_checknumber(L, 8);
    double v2x = luaL_checknumber(L, 9);
    double v2y = luaL_checknumber(L, 10);
    double v2z = luaL_checknumber(L, 11);
    double v3x = luaL_checknumber(L, 12);
    double v3y = luaL_checknumber(L, 13);
    double v3z = luaL_checknumber(L, 14);
    int idx = alea_wed_surface(sys, id, vx, vy, vz,
                               v1x, v1y, v1z, v2x, v2y, v2z, v3x, v3y, v3z);
    if (idx < 0) return luaL_error(L, "wed: %s", alea_error());
    lua_pushinteger(L, idx);
    return 1;
}

/* sys:box_general(id, cx,cy,cz, v1x,v1y,v1z, v2x,v2y,v2z, v3x,v3y,v3z) */
static int l_box_general(lua_State* L) {
    alea_system_t* sys = alea_get_sys(L, 1);
    int id     = (int)luaL_checkinteger(L, 2);
    double cx  = luaL_checknumber(L, 3);
    double cy  = luaL_checknumber(L, 4);
    double cz  = luaL_checknumber(L, 5);
    double v1x = luaL_checknumber(L, 6);
    double v1y = luaL_checknumber(L, 7);
    double v1z = luaL_checknumber(L, 8);
    double v2x = luaL_checknumber(L, 9);
    double v2y = luaL_checknumber(L, 10);
    double v2z = luaL_checknumber(L, 11);
    double v3x = luaL_checknumber(L, 12);
    double v3y = luaL_checknumber(L, 13);
    double v3z = luaL_checknumber(L, 14);
    int idx = alea_box_general_surface(sys, id, cx, cy, cz,
                                       v1x, v1y, v1z, v2x, v2y, v2z, v3x, v3y, v3z);
    if (idx < 0) return luaL_error(L, "box_general: %s", alea_error());
    lua_pushinteger(L, idx);
    return 1;
}

/* sys:rhp(id, bx,by,bz, hx,hy,hz, r1x,r1y,r1z, r2x,r2y,r2z, r3x,r3y,r3z) */
static int l_rhp(lua_State* L) {
    alea_system_t* sys = alea_get_sys(L, 1);
    int id     = (int)luaL_checkinteger(L, 2);
    double bx  = luaL_checknumber(L, 3);
    double by  = luaL_checknumber(L, 4);
    double bz  = luaL_checknumber(L, 5);
    double hx  = luaL_checknumber(L, 6);
    double hy  = luaL_checknumber(L, 7);
    double hz  = luaL_checknumber(L, 8);
    double r1x = luaL_checknumber(L, 9);
    double r1y = luaL_checknumber(L, 10);
    double r1z = luaL_checknumber(L, 11);
    double r2x = luaL_checknumber(L, 12);
    double r2y = luaL_checknumber(L, 13);
    double r2z = luaL_checknumber(L, 14);
    double r3x = luaL_checknumber(L, 15);
    double r3y = luaL_checknumber(L, 16);
    double r3z = luaL_checknumber(L, 17);
    int idx = alea_rhp_surface(sys, id, bx, by, bz, hx, hy, hz,
                               r1x, r1y, r1z, r2x, r2y, r2z, r3x, r3y, r3z);
    if (idx < 0) return luaL_error(L, "rhp: %s", alea_error());
    lua_pushinteger(L, idx);
    return 1;
}

/* ============================================================================
 * Registration
 * ============================================================================ */

static const luaL_Reg surface_methods[] = {
    {"plane",       l_plane},
    {"sphere",      l_sphere},
    {"cylinder_z",  l_cylinder_z},
    {"cylinder_x",  l_cylinder_x},
    {"cylinder_y",  l_cylinder_y},
    {"cone_z",      l_cone_z},
    {"cone_x",      l_cone_x},
    {"cone_y",      l_cone_y},
    {"torus_z",     l_torus_z},
    {"torus_x",     l_torus_x},
    {"torus_y",     l_torus_y},
    {"box",         l_box},
    {"quadric",     l_quadric},
    {"rcc",         l_rcc},
    {"sph",         l_sph},
    {"trc",         l_trc},
    {"ell",         l_ell},
    {"rec",         l_rec},
    {"wed",         l_wed},
    {"box_general", l_box_general},
    {"rhp",         l_rhp},
    {NULL, NULL}
};

int luaopen_alea_surfaces(lua_State* L) {
    /* Add surface methods to System metatable's __index */
    luaL_getmetatable(L, ALEA_SYSTEM_MT);
    lua_getfield(L, -1, "__index");
    luaL_setfuncs(L, surface_methods, 0);
    lua_pop(L, 2);
    return 0;
}
