// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

#include "alea_lua.h"
#include <math.h>
#include <string.h>

static int check_finite_sequence(lua_State* L, int index, const char* name,
                                 double* values, int expected) {
    luaL_checktype(L, index, LUA_TTABLE);
    lua_Integer count = luaL_len(L, index);
    if (count != expected)
        return luaL_error(L, "%s must contain exactly %d values", name, expected);
    for (int i = 0; i < expected; ++i) {
        lua_geti(L, index, i + 1);
        values[i] = luaL_checknumber(L, -1);
        lua_pop(L, 1);
        if (!isfinite(values[i]))
            return luaL_error(L, "%s values must be finite", name);
    }
    return 0;
}

static int primitive_parameter_count(int type) {
    switch ((alea_primitive_type_t)type) {
        case ALEA_PRIMITIVE_PLANE:
        case ALEA_PRIMITIVE_SPHERE:
        case ALEA_PRIMITIVE_SPH: return 4;
        case ALEA_PRIMITIVE_CYLINDER_X:
        case ALEA_PRIMITIVE_CYLINDER_Y:
        case ALEA_PRIMITIVE_CYLINDER_Z: return 3;
        case ALEA_PRIMITIVE_CONE_X:
        case ALEA_PRIMITIVE_CONE_Y:
        case ALEA_PRIMITIVE_CONE_Z:
        case ALEA_PRIMITIVE_TORUS_X:
        case ALEA_PRIMITIVE_TORUS_Y:
        case ALEA_PRIMITIVE_TORUS_Z: return 5;
        case ALEA_PRIMITIVE_RPP: return 6;
        case ALEA_PRIMITIVE_QUADRIC: return 10;
        case ALEA_PRIMITIVE_RCC:
        case ALEA_PRIMITIVE_ELL: return 7;
        case ALEA_PRIMITIVE_TRC: return 8;
        case ALEA_PRIMITIVE_BOX:
        case ALEA_PRIMITIVE_REC:
        case ALEA_PRIMITIVE_WED: return 12;
        case ALEA_PRIMITIVE_RHP: return 15;
        default: return -1;
    }
}

/* alea.primitive_evaluate(type, parameters, {x, y, z}) -> number */
static int l_primitive_evaluate(lua_State* L) {
    int type = (int)luaL_checkinteger(L, 1);
    int count = primitive_parameter_count(type);
    if (count < 0) return luaL_error(L, "unsupported primitive type");

    double p[15];
    double xyz[3];
    check_finite_sequence(L, 2, "parameters", p, count);
    check_finite_sequence(L, 3, "point", xyz, 3);

    alea_primitive_data_t data = {0};
#define COPY_FIELDS(member, field_count) \
    memcpy(&(member), p, (field_count) * sizeof(double))
    switch ((alea_primitive_type_t)type) {
        case ALEA_PRIMITIVE_PLANE: COPY_FIELDS(data.plane, 4); break;
        case ALEA_PRIMITIVE_SPHERE: COPY_FIELDS(data.sphere, 4); break;
        case ALEA_PRIMITIVE_SPH: COPY_FIELDS(data.sph, 4); break;
        case ALEA_PRIMITIVE_CYLINDER_X: COPY_FIELDS(data.cyl_x, 3); break;
        case ALEA_PRIMITIVE_CYLINDER_Y: COPY_FIELDS(data.cyl_y, 3); break;
        case ALEA_PRIMITIVE_CYLINDER_Z: COPY_FIELDS(data.cyl_z, 3); break;
        case ALEA_PRIMITIVE_CONE_X:
            COPY_FIELDS(data.cone_x, 4);
            if (p[4] != -1.0 && p[4] != 0.0 && p[4] != 1.0)
                return luaL_error(L, "cone sheet must be -1, 0, or 1");
            data.cone_x.sheet_selection = (int)p[4];
            break;
        case ALEA_PRIMITIVE_CONE_Y:
            COPY_FIELDS(data.cone_y, 4);
            if (p[4] != -1.0 && p[4] != 0.0 && p[4] != 1.0)
                return luaL_error(L, "cone sheet must be -1, 0, or 1");
            data.cone_y.sheet_selection = (int)p[4];
            break;
        case ALEA_PRIMITIVE_CONE_Z:
            COPY_FIELDS(data.cone_z, 4);
            if (p[4] != -1.0 && p[4] != 0.0 && p[4] != 1.0)
                return luaL_error(L, "cone sheet must be -1, 0, or 1");
            data.cone_z.sheet_selection = (int)p[4];
            break;
        case ALEA_PRIMITIVE_RPP: COPY_FIELDS(data.box, 6); break;
        case ALEA_PRIMITIVE_QUADRIC:
            memcpy(data.quadric.coeffs, p, 10 * sizeof(double));
            break;
        case ALEA_PRIMITIVE_TORUS_X:
        case ALEA_PRIMITIVE_TORUS_Y:
        case ALEA_PRIMITIVE_TORUS_Z:
            data.torus.axis = type == ALEA_PRIMITIVE_TORUS_X ? ALEA_AXIS_X :
                              type == ALEA_PRIMITIVE_TORUS_Y ? ALEA_AXIS_Y : ALEA_AXIS_Z;
            data.torus.center_x = p[0];
            data.torus.center_y = p[1];
            data.torus.center_z = p[2];
            data.torus.major_radius = p[3];
            data.torus.minor_radius = p[4];
            data.torus.axial_semiwidth_B = p[4];
            break;
        case ALEA_PRIMITIVE_RCC: COPY_FIELDS(data.rcc, 7); break;
        case ALEA_PRIMITIVE_BOX: COPY_FIELDS(data.box_general, 12); break;
        case ALEA_PRIMITIVE_TRC: COPY_FIELDS(data.trc, 8); break;
        case ALEA_PRIMITIVE_ELL: COPY_FIELDS(data.ell, 7); break;
        case ALEA_PRIMITIVE_REC: COPY_FIELDS(data.rec, 12); break;
        case ALEA_PRIMITIVE_WED: COPY_FIELDS(data.wed, 12); break;
        case ALEA_PRIMITIVE_RHP: COPY_FIELDS(data.rhp, 15); break;
        default: return luaL_error(L, "unsupported primitive type");
    }
#undef COPY_FIELDS

    double value;
    alea_error_t error = alea_primitive_evaluate_checked(
        (alea_primitive_type_t)type, &data, xyz[0], xyz[1], xyz[2], &value);
    if (error != ALEA_OK)
        return luaL_error(L, "primitive_evaluate: %s", alea_error());
    lua_pushnumber(L, value);
    return 1;
}

static int check_cone_sheet(lua_State* L, int index) {
    int sheet = (int)luaL_optinteger(L, index, 0);
    if (sheet < -1 || sheet > 1)
        luaL_argerror(L, index, "sheet must be -1, 0, or 1");
    return sheet;
}

static int boundary_from_name(lua_State* L, int index) {
    const char* name = luaL_checkstring(L, index);
    if (strcmp(name, "transmissive") == 0) return ALEA_BOUNDARY_TRANSMISSIVE;
    if (strcmp(name, "reflective") == 0) return ALEA_BOUNDARY_REFLECTIVE;
    if (strcmp(name, "white") == 0) return ALEA_BOUNDARY_WHITE;
    if (strcmp(name, "periodic") == 0) return ALEA_BOUNDARY_PERIODIC;
    if (strcmp(name, "vacuum") == 0) return ALEA_BOUNDARY_VACUUM;
    return luaL_argerror(L, index,
        "boundary must be transmissive, reflective, white, periodic, or vacuum");
}

static const char* boundary_name(alea_boundary_type_t boundary) {
    switch (boundary) {
        case ALEA_BOUNDARY_REFLECTIVE: return "reflective";
        case ALEA_BOUNDARY_WHITE: return "white";
        case ALEA_BOUNDARY_PERIODIC: return "periodic";
        case ALEA_BOUNDARY_VACUUM: return "vacuum";
        default: return "transmissive";
    }
}

/* sys:surface_set_boundary(surface_id, name) */
static int l_surface_set_boundary(lua_State* L) {
    alea_system_t* sys = alea_get_sys(L, 1);
    int surface_id = (int)luaL_checkinteger(L, 2);
    alea_boundary_type_t boundary =
        (alea_boundary_type_t)boundary_from_name(L, 3);
    if (alea_surface_set_boundary(sys, surface_id, boundary) != 0)
        return luaL_error(L, "surface_set_boundary: surface ID %d not found", surface_id);
    return 0;
}

/* sys:surface_set_periodic_pair(first_surface_id, second_surface_id) */
static int l_surface_set_periodic_pair(lua_State* L) {
    alea_system_t* sys = alea_get_sys(L, 1);
    int first = (int)luaL_checkinteger(L, 2);
    int second = (int)luaL_checkinteger(L, 3);
    if (alea_surface_set_periodic_pair(sys, first, second) != 0)
        return luaL_error(L,
            "surface_set_periodic_pair: expected two unpaired plane surfaces");
    return 0;
}

/* sys:surface_get_boundary(surface_id) -> name */
static int l_surface_get_boundary(lua_State* L) {
    alea_system_t* sys = alea_get_sys(L, 1);
    int surface_id = (int)luaL_checkinteger(L, 2);
    int index = alea_surface_find(sys, surface_id);
    alea_boundary_type_t boundary;
    if (index < 0 || alea_surface_get(sys, (size_t)index, NULL, NULL, NULL, NULL,
                                     &boundary) != 0)
        return luaL_error(L, "surface_get_boundary: surface ID %d not found", surface_id);
    lua_pushstring(L, boundary_name(boundary));
    return 1;
}

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

/* sys:cone_z(id, cx, cy, cz, t2[, sheet]) */
static int l_cone_z(lua_State* L) {
    alea_system_t* sys = alea_get_sys(L, 1);
    int id    = (int)luaL_checkinteger(L, 2);
    double cx = luaL_checknumber(L, 3);
    double cy = luaL_checknumber(L, 4);
    double cz = luaL_checknumber(L, 5);
    double t2 = luaL_checknumber(L, 6);
    int sheet = check_cone_sheet(L, 7);
    int idx = sheet == 0 ? alea_cone_z_surface(sys, id, cx, cy, cz, t2) :
                          alea_cone_z_surface_sheet(sys, id, cx, cy, cz, t2, sheet);
    if (idx < 0) return luaL_error(L, "cone_z: %s", alea_error());
    lua_pushinteger(L, idx);
    return 1;
}

/* sys:cone_x(id, cx, cy, cz, t2[, sheet]) */
static int l_cone_x(lua_State* L) {
    alea_system_t* sys = alea_get_sys(L, 1);
    int id    = (int)luaL_checkinteger(L, 2);
    double cx = luaL_checknumber(L, 3);
    double cy = luaL_checknumber(L, 4);
    double cz = luaL_checknumber(L, 5);
    double t2 = luaL_checknumber(L, 6);
    int sheet = check_cone_sheet(L, 7);
    int idx = sheet == 0 ? alea_cone_x_surface(sys, id, cx, cy, cz, t2) :
                          alea_cone_x_surface_sheet(sys, id, cx, cy, cz, t2, sheet);
    if (idx < 0) return luaL_error(L, "cone_x: %s", alea_error());
    lua_pushinteger(L, idx);
    return 1;
}

/* sys:cone_y(id, cx, cy, cz, t2[, sheet]) */
static int l_cone_y(lua_State* L) {
    alea_system_t* sys = alea_get_sys(L, 1);
    int id    = (int)luaL_checkinteger(L, 2);
    double cx = luaL_checknumber(L, 3);
    double cy = luaL_checknumber(L, 4);
    double cz = luaL_checknumber(L, 5);
    double t2 = luaL_checknumber(L, 6);
    int sheet = check_cone_sheet(L, 7);
    int idx = sheet == 0 ? alea_cone_y_surface(sys, id, cx, cy, cz, t2) :
                          alea_cone_y_surface_sheet(sys, id, cx, cy, cz, t2, sheet);
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
    {"surface_set_boundary", l_surface_set_boundary},
    {"surface_set_periodic_pair", l_surface_set_periodic_pair},
    {"surface_get_boundary", l_surface_get_boundary},
    {NULL, NULL}
};

int luaopen_alea_surfaces(lua_State* L) {
    /* Add surface methods to System metatable's __index */
    luaL_getmetatable(L, ALEA_SYSTEM_MT);
    lua_getfield(L, -1, "__index");
    luaL_setfuncs(L, surface_methods, 0);
    lua_pop(L, 2);

    lua_pushcfunction(L, l_primitive_evaluate);
    lua_setfield(L, -2, "primitive_evaluate");

#define SET_PRIMITIVE(short_name, value) do { \
    lua_pushinteger(L, (lua_Integer)(value)); \
    lua_setfield(L, -2, "PRIMITIVE_" #short_name); \
} while (0)
    SET_PRIMITIVE(PLANE, ALEA_PRIMITIVE_PLANE);
    SET_PRIMITIVE(SPHERE, ALEA_PRIMITIVE_SPHERE);
    SET_PRIMITIVE(CYLINDER_X, ALEA_PRIMITIVE_CYLINDER_X);
    SET_PRIMITIVE(CYLINDER_Y, ALEA_PRIMITIVE_CYLINDER_Y);
    SET_PRIMITIVE(CYLINDER_Z, ALEA_PRIMITIVE_CYLINDER_Z);
    SET_PRIMITIVE(CONE_X, ALEA_PRIMITIVE_CONE_X);
    SET_PRIMITIVE(CONE_Y, ALEA_PRIMITIVE_CONE_Y);
    SET_PRIMITIVE(CONE_Z, ALEA_PRIMITIVE_CONE_Z);
    SET_PRIMITIVE(RPP, ALEA_PRIMITIVE_RPP);
    SET_PRIMITIVE(QUADRIC, ALEA_PRIMITIVE_QUADRIC);
    SET_PRIMITIVE(TORUS_X, ALEA_PRIMITIVE_TORUS_X);
    SET_PRIMITIVE(TORUS_Y, ALEA_PRIMITIVE_TORUS_Y);
    SET_PRIMITIVE(TORUS_Z, ALEA_PRIMITIVE_TORUS_Z);
    SET_PRIMITIVE(RCC, ALEA_PRIMITIVE_RCC);
    SET_PRIMITIVE(BOX, ALEA_PRIMITIVE_BOX);
    SET_PRIMITIVE(SPH, ALEA_PRIMITIVE_SPH);
    SET_PRIMITIVE(TRC, ALEA_PRIMITIVE_TRC);
    SET_PRIMITIVE(ELL, ALEA_PRIMITIVE_ELL);
    SET_PRIMITIVE(REC, ALEA_PRIMITIVE_REC);
    SET_PRIMITIVE(WED, ALEA_PRIMITIVE_WED);
    SET_PRIMITIVE(RHP, ALEA_PRIMITIVE_RHP);
#undef SET_PRIMITIVE
    return 0;
}
