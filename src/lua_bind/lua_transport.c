// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

#include "alea_lua.h"
#include "alea_transport.h"
#include "alea_source.h"
#include <math.h>
#include <limits.h>
#include <stdlib.h>
#include <string.h>

#define TRANSPORT_GUARD_MT "alea.TransportGuard"
#define SOURCE_MT "alea.Source"
#define SOURCE_BUFFERS_MT "alea.SourceBuffers"

typedef struct {
    double *angle_mu, *angle_pdf, *energy_values, *energy_weights;
    double *r_edges, *z_edges, *emissivity;
} source_buffers_t;

static void source_buffers_release(source_buffers_t* buffers) {
    free(buffers->angle_mu); free(buffers->angle_pdf);
    free(buffers->energy_values); free(buffers->energy_weights);
    free(buffers->r_edges); free(buffers->z_edges);
    free(buffers->emissivity);
    memset(buffers, 0, sizeof(*buffers));
}

static int source_buffers_gc(lua_State* L) {
    source_buffers_t* buffers = luaL_checkudata(L, 1, SOURCE_BUFFERS_MT);
    source_buffers_release(buffers);
    return 0;
}

typedef struct {
    alea_tally_plan_t* plan;
    alea_nuc_cell_bindings_t* bindings;
    alea_source_t* owned_source;
    alea_transport_result_t result;
} transport_guard_t;

static int guard_gc(lua_State* L) {
    transport_guard_t* guard = luaL_checkudata(L, 1, TRANSPORT_GUARD_MT);
    alea_transport_result_free(&guard->result);
    alea_nuc_cell_bindings_free(guard->bindings);
    alea_source_free(guard->owned_source);
    alea_tally_plan_free(guard->plan);
    return 0;
}

static double number_field(lua_State* L, int table, const char* key, double fallback) {
    lua_getfield(L, table, key);
    double value = lua_isnil(L, -1) ? fallback : luaL_checknumber(L, -1);
    lua_pop(L, 1);
    return value;
}

static lua_Integer integer_field(lua_State* L, int table, const char* key,
                                 lua_Integer fallback) {
    lua_getfield(L, table, key);
    lua_Integer value = lua_isnil(L, -1) ? fallback : luaL_checkinteger(L, -1);
    lua_pop(L, 1);
    return value;
}

static int enum_field(lua_State* L, int table, const char* key,
                      const char* const* names, int count, int fallback) {
    lua_getfield(L, table, key);
    if (lua_isnil(L, -1)) { lua_pop(L, 1); return fallback; }
    const char* value = luaL_checkstring(L, -1);
    for (int i = 0; i < count; ++i) {
        if (strcmp(value, names[i]) == 0) { lua_pop(L, 1); return i; }
    }
    return luaL_error(L, "unknown %s: %s", key, value);
}

static void vec3_field(lua_State* L, int table, const char* key,
                       double out[3], int required) {
    lua_getfield(L, table, key);
    if (lua_isnil(L, -1) && !required) { lua_pop(L, 1); return; }
    luaL_checktype(L, -1, LUA_TTABLE);
    if (lua_rawlen(L, -1) != 3) luaL_error(L, "%s must have three elements", key);
    for (int i = 0; i < 3; ++i) {
        lua_rawgeti(L, -1, i + 1);
        out[i] = luaL_checknumber(L, -1);
        lua_pop(L, 1);
    }
    lua_pop(L, 1);
}

static double* source_number_array(lua_State* L, int table, const char* key,
                                   size_t* count) {
    lua_getfield(L, table, key);
    luaL_checktype(L, -1, LUA_TTABLE);
    size_t n = lua_rawlen(L, -1);
    if (n < 2 || n > UINT32_MAX || n > SIZE_MAX / sizeof(double))
        luaL_error(L, "%s needs at least two values", key);
    double* values = malloc(n * sizeof(double));
    if (!values) luaL_error(L, "out of memory");
    for (size_t i = 0; i < n; ++i) {
        lua_rawgeti(L, -1, (lua_Integer)i + 1);
        if (!lua_isnumber(L, -1)) {
            free(values); luaL_error(L, "%s values must be numeric", key);
        }
        values[i] = lua_tonumber(L, -1);
        lua_pop(L, 1);
    }
    lua_pop(L, 1);
    *count = n;
    return values;
}

static void set_integer(lua_State* L, const char* key, lua_Integer value) {
    lua_pushinteger(L, value); lua_setfield(L, -2, key);
}

static void set_number(lua_State* L, const char* key, double value) {
    lua_pushnumber(L, value); lua_setfield(L, -2, key);
}

static void push_floats(lua_State* L, const double* values, size_t count) {
    lua_createtable(L, (int)count, 0);
    for (size_t i = 0; i < count; ++i) {
        lua_pushnumber(L, values[i]); lua_rawseti(L, -2, (lua_Integer)i + 1);
    }
}

static void push_vec3(lua_State* L, const double values[3]) {
    push_floats(L, values, 3);
}

static const char* const scores[] = {"track_length", "collision", "reaction_event",
    "reaction_rate", "heating", "local_deposition"};
static const char* const domains[] = {"cell", "universe", "mesh"};
static const char* const particles[] = {"neutron", "photon", "all"};

static int has_field(lua_State* L, int table, const char* key) {
    lua_getfield(L, table, key);
    int present = !lua_isnil(L, -1);
    lua_pop(L, 1);
    return present;
}

static alea_source_t* source_from_table_depth(lua_State* L, int idx,
                                               unsigned depth);

static alea_source_t* source_mixture_from_table(lua_State* L, int idx,
                                                 unsigned depth) {
    lua_getfield(L, idx, "components");
    luaL_checktype(L, -1, LUA_TTABLE);
    int list = lua_absindex(L, -1);
    size_t count = lua_rawlen(L, list);
    if (!count || count > UINT32_MAX ||
        count > SIZE_MAX / sizeof(alea_source_t*) ||
        count > SIZE_MAX / sizeof(double))
        luaL_error(L, "mixture needs at least one component");
    lua_newtable(L);
    int guards = lua_absindex(L, -1);
    alea_source_t** parts = lua_newuserdatauv(L,
        count*sizeof(alea_source_t*), 0);
    double* strengths = lua_newuserdatauv(L, count*sizeof(double), 0);
    for (size_t i = 0; i < count; ++i) {
        lua_rawgeti(L, list, (lua_Integer)i + 1);
        luaL_checktype(L, -1, LUA_TTABLE);
        int component = lua_absindex(L, -1);
        lua_getfield(L, component, "strength");
        strengths[i] = luaL_checknumber(L, -1);
        lua_pop(L, 1);
        lua_getfield(L, component, "source");
        luaL_checktype(L, -1, LUA_TTABLE);
        alea_source_t** guard = lua_newuserdatauv(L, sizeof(*guard), 0);
        *guard = NULL;
        luaL_setmetatable(L, SOURCE_MT);
        *guard = source_from_table_depth(L, lua_absindex(L, -2), depth + 1);
        parts[i] = *guard;
        lua_rawseti(L, guards, (lua_Integer)i + 1);
        lua_pop(L, 2);
    }
    alea_source_t* mixture = NULL;
    alea_error_t err = alea_source_mixture_prepare(parts, strengths, count,
                                                   &mixture);
    if (err != ALEA_OK)
        luaL_error(L, "invalid source mixture: %s", alea_error_string(err));
    for (size_t i = 0; i < count; ++i) {
        lua_rawgeti(L, guards, (lua_Integer)i + 1);
        *(alea_source_t**)luaL_checkudata(L, -1, SOURCE_MT) = NULL;
        lua_pop(L, 1);
    }
    lua_pop(L, 4);
    return mixture;
}

static alea_source_t* source_from_table_depth(lua_State* L, int idx,
                                               unsigned depth) {
    idx = lua_absindex(L, idx);
    luaL_checktype(L, idx, LUA_TTABLE);
    if (depth > 8) luaL_error(L, "source mixtures are nested too deeply");
    lua_getfield(L, idx, "type");
    if (!lua_isnil(L, -1)) {
        const char* type = luaL_checkstring(L, -1);
        if (strcmp(type, "mixture") != 0)
            luaL_error(L, "unknown source type: %s", type);
        lua_pop(L, 1);
        return source_mixture_from_table(L, idx, depth);
    }
    lua_pop(L, 1);
    source_buffers_t* buffers = lua_newuserdatauv(L, sizeof(*buffers), 0);
    memset(buffers, 0, sizeof(*buffers));
    luaL_setmetatable(L, SOURCE_BUFFERS_MT);
    alea_source_spec_t spec = {0};
    double *angle_mu = NULL, *angle_pdf = NULL;
    double *r_edges = NULL, *z_edges = NULL, *emissivity = NULL;
    spec.weight = number_field(L, idx, "weight", 1);
    spec.particle = (alea_nuc_particle_t)enum_field(L, idx, "particle", particles, 2, 0);
    if (has_field(L, idx, "kind") || has_field(L, idx, "position") ||
        has_field(L, idx, "direction") || has_field(L, idx, "lower") ||
        has_field(L, idx, "upper"))
        luaL_error(L, "put spatial and angular fields inside space and angle");
    {
        lua_getfield(L, idx, "space");
        luaL_checktype(L, -1, LUA_TTABLE);
        static const char* const spaces[] = {"point", "box", "line", "sphere", "cylinder", "tokamak_rz"};
        int spatial = enum_field(L, lua_absindex(L, -1), "type", spaces, 6, -1);
        if (spatial < 0) luaL_error(L, "space requires type");
        spec.space = (alea_source_space_t)spatial;
        if (spatial == 0) vec3_field(L, lua_absindex(L, -1), "position", spec.position, 1);
        else if (spatial == 1) {
            vec3_field(L, lua_absindex(L, -1), "lower", spec.lower, 1);
            vec3_field(L, lua_absindex(L, -1), "upper", spec.upper, 1);
        } else if (spatial == 2) {
            vec3_field(L, lua_absindex(L, -1), "start", spec.start, 1);
            vec3_field(L, lua_absindex(L, -1), "end", spec.end, 1);
        } else if (spatial == 5) {
            spec.phi_min = number_field(L, lua_absindex(L, -1), "phi_min", 0);
            spec.phi_max = number_field(L, lua_absindex(L, -1), "phi_max", 0);
        } else {
            if (spatial == 3)
                vec3_field(L, lua_absindex(L, -1), "center", spec.center, 1);
            else {
                vec3_field(L, lua_absindex(L, -1), "base", spec.base, 1);
                vec3_field(L, lua_absindex(L, -1), "axis", spec.axis, 1);
            }
            lua_getfield(L, -1, "outer_radius");
            spec.outer_radius = luaL_checknumber(L, -1);
            lua_pop(L, 1);
            spec.inner_radius = number_field(L, lua_absindex(L, -1), "inner_radius", 0);
        }
        lua_pop(L, 1);
        lua_getfield(L, idx, "angle");
        luaL_checktype(L, -1, LUA_TTABLE);
        static const char* const angles[] = {
            "monodirectional", "isotropic", "cone", "cosine", "tabulated_mu", "radial"
        };
        int angular = enum_field(L, lua_absindex(L, -1), "type", angles, 6, -1);
        if (angular < 0) luaL_error(L, "angle requires type");
        spec.angle = (alea_source_angle_t)angular;
        if (angular == 0 || angular == 2 || angular == 3 || angular == 4)
            vec3_field(L, lua_absindex(L, -1), "direction", spec.direction, 1);
        if (angular == 2) {
            lua_getfield(L, -1, "half_angle");
            spec.cone_half_angle = luaL_checknumber(L, -1);
            lua_pop(L, 1);
        } else if (angular == 4) {
            static const char* const modes[] = {"histogram", "linear"};
            int mode = enum_field(L, lua_absindex(L, -1), "interpolation", modes, 2, -1);
            if (mode < 0) luaL_error(L, "tabulated_mu requires interpolation");
            spec.angle_interpolation = (alea_source_pdf_t)mode;
            lua_getfield(L, -1, "mu");
            luaL_checktype(L, -1, LUA_TTABLE);
            lua_getfield(L, -2, "pdf");
            luaL_checktype(L, -1, LUA_TTABLE);
            size_t count = lua_rawlen(L, -2);
            if (count < 2 || count != lua_rawlen(L, -1) ||
                count > UINT32_MAX || count > SIZE_MAX / sizeof(double))
                luaL_error(L, "tabulated_mu needs equal arrays of at least two entries");
            angle_mu = malloc(count * sizeof(double));
            angle_pdf = malloc(count * sizeof(double));
            buffers->angle_mu = angle_mu;
            buffers->angle_pdf = angle_pdf;
            if (!angle_mu || !angle_pdf) luaL_error(L, "out of memory");
            for (size_t i = 0; i < count; ++i) {
                lua_rawgeti(L, -2, (lua_Integer)i + 1);
                if (!lua_isnumber(L, -1)) {
                    luaL_error(L, "mu values must be numeric");
                }
                angle_mu[i] = lua_tonumber(L, -1);
                lua_pop(L, 1);
                lua_rawgeti(L, -1, (lua_Integer)i + 1);
                if (!lua_isnumber(L, -1)) {
                    luaL_error(L, "pdf values must be numeric");
                }
                angle_pdf[i] = lua_tonumber(L, -1);
                lua_pop(L, 1);
            }
            spec.angle_mu = angle_mu;
            spec.angle_pdf = angle_pdf;
            spec.angle_count = count;
            lua_pop(L, 2);
        } else if (angular == 5) {
            vec3_field(L, lua_absindex(L, -1), "origin", spec.angle_origin, 1);
            lua_getfield(L, -1, "inward");
            if (!lua_isnil(L, -1)) {
                luaL_checktype(L, -1, LUA_TBOOLEAN);
                spec.radial_inward = lua_toboolean(L, -1);
            }
            lua_pop(L, 1);
        }
        lua_pop(L, 1);
    }
    lua_getfield(L, idx, "time");
    if (lua_istable(L, -1)) {
        static const char* const time_types[] = {"constant"};
        if (enum_field(L, lua_absindex(L, -1), "type", time_types, 1, -1) < 0)
            luaL_error(L, "time requires type");
        lua_getfield(L, -1, "value");
        spec.time = luaL_checknumber(L, -1);
        lua_pop(L, 1);
    } else if (!lua_isnil(L, -1)) spec.time = luaL_checknumber(L, -1);
    lua_pop(L, 1);
    double *values = NULL, *weights = NULL;
    lua_getfield(L, idx, "energy");
    if (lua_istable(L, -1)) {
        static const char* const energy_types[] = {"mono", "lines", "tabulated"};
        int type = enum_field(L, lua_absindex(L, -1), "type", energy_types, 3, -1);
        if (type < 0) luaL_error(L, "energy requires type");
        spec.energy_type = (alea_source_energy_t)type;
        if (type == 0) {
            lua_getfield(L, -1, "value");
            spec.energy = luaL_checknumber(L, -1);
            lua_pop(L, 1);
        } else {
            if (type == 2) {
                static const char* const modes[] = {"histogram", "linear"};
                int mode = enum_field(L, lua_absindex(L, -1),
                                      "interpolation", modes, 2, -1);
                if (mode < 0) luaL_error(L, "tabulated energy requires interpolation");
                spec.energy_interpolation = (alea_source_pdf_t)mode;
            }
            lua_getfield(L, -1, "values");
            luaL_checktype(L, -1, LUA_TTABLE);
            lua_getfield(L, -2, type == 1 ? "weights" : "pdf");
            luaL_checktype(L, -1, LUA_TTABLE);
            size_t count = lua_rawlen(L, -2);
            if (count < (type == 1 ? 1u : 2u) || count != lua_rawlen(L, -1) ||
                count > UINT32_MAX || count > SIZE_MAX / sizeof(double))
                luaL_error(L, "energy arrays have invalid lengths");
            values = malloc(count * sizeof(double));
            weights = malloc(count * sizeof(double));
            buffers->energy_values = values;
            buffers->energy_weights = weights;
            if (!values || !weights) luaL_error(L, "out of memory");
            for (size_t i = 0; i < count; ++i) {
                lua_rawgeti(L, -2, (lua_Integer)i + 1);
                if (!lua_isnumber(L, -1)) {
                    luaL_error(L, "energy values must be numeric");
                }
                values[i] = lua_tonumber(L, -1);
                lua_pop(L, 1);
                lua_rawgeti(L, -1, (lua_Integer)i + 1);
                if (!lua_isnumber(L, -1)) {
                    luaL_error(L, "energy weights must be numeric");
                }
                weights[i] = lua_tonumber(L, -1);
                lua_pop(L, 1);
            }
            spec.energy_values = values;
            spec.energy_weights = weights;
            spec.energy_count = count;
            lua_pop(L, 2);
        }
    } else spec.energy = luaL_checknumber(L, -1);
    lua_pop(L, 1);
    if (spec.space == ALEA_SOURCE_TOKAMAK_RZ) {
        lua_getfield(L, idx, "space");
        int space = lua_absindex(L, -1);
        r_edges = source_number_array(L, space, "r_edges", &spec.r_edge_count);
        buffers->r_edges = r_edges;
        z_edges = source_number_array(L, space, "z_edges", &spec.z_edge_count);
        buffers->z_edges = z_edges;
        size_t nr = spec.r_edge_count - 1, nz = spec.z_edge_count - 1;
        if (nr > SIZE_MAX / nz || nr*nz > UINT32_MAX) {
            luaL_error(L, "tokamak_rz grid is too large");
        }
        lua_getfield(L, space, "emissivity");
        luaL_checktype(L, -1, LUA_TTABLE);
        if (lua_rawlen(L, -1) != nr) {
            luaL_error(L, "emissivity needs one row per R bin");
        }
        emissivity = malloc(nr*nz*sizeof(double));
        buffers->emissivity = emissivity;
        if (!emissivity) luaL_error(L, "out of memory");
        for (size_t ir = 0; ir < nr; ++ir) {
            lua_rawgeti(L, -1, (lua_Integer)ir + 1);
            if (!lua_istable(L, -1) || lua_rawlen(L, -1) != nz) {
                luaL_error(L, "emissivity row needs one value per Z bin");
            }
            for (size_t iz = 0; iz < nz; ++iz) {
                lua_rawgeti(L, -1, (lua_Integer)iz + 1);
                if (!lua_isnumber(L, -1)) {
                    luaL_error(L, "emissivity values must be numeric");
                }
                emissivity[ir*nz + iz] = lua_tonumber(L, -1);
                lua_pop(L, 1);
            }
            lua_pop(L, 1);
        }
        lua_pop(L, 2);
        spec.r_edges = r_edges;
        spec.z_edges = z_edges;
        spec.rz_emissivity = emissivity;
    }
    alea_source_t* source = NULL;
    alea_error_t err = alea_source_prepare(&spec, &source);
    source_buffers_release(buffers);
    lua_pop(L, 1);
    if (err != ALEA_OK) luaL_error(L, "invalid source: %s", alea_error_string(err));
    return source;
}

static alea_source_t* source_from_table(lua_State* L, int idx) {
    return source_from_table_depth(L, idx, 0);
}

static int source_gc(lua_State* L) {
    alea_source_t** source = luaL_checkudata(L, 1, SOURCE_MT);
    alea_source_free(*source);
    *source = NULL;
    return 0;
}

static int source_prepare(lua_State* L) {
    luaL_checktype(L, 1, LUA_TTABLE);
    alea_source_t** source = lua_newuserdatauv(L, sizeof(*source), 0);
    *source = NULL;
    luaL_setmetatable(L, SOURCE_MT);
    *source = source_from_table(L, 1);
    return 1;
}

static int source_integrated_emissivity(lua_State* L) {
    alea_source_t** source = luaL_checkudata(L, 1, SOURCE_MT);
    double value;
    if (alea_source_integrated_emissivity(*source, &value) != ALEA_OK)
        lua_pushnil(L);
    else lua_pushnumber(L, value);
    return 1;
}

static int sample_source(lua_State* L) {
    alea_source_t* owned = NULL;
    alea_source_t* source;
    if (luaL_testudata(L, 1, SOURCE_MT)) {
        source = *(alea_source_t**)luaL_checkudata(L, 1, SOURCE_MT);
        if (!source) return luaL_error(L, "source has been destroyed");
    } else {
        /* The userdata guard survives Lua errors while building output tables. */
        alea_source_t** guard = lua_newuserdatauv(L, sizeof(*guard), 0);
        *guard = NULL;
        luaL_setmetatable(L, SOURCE_MT);
        owned = source = source_from_table(L, 1);
        *guard = owned;
    }
    lua_Integer count = luaL_optinteger(L, 2, 1);
    lua_Integer seed = luaL_optinteger(L, 3, 1);
    lua_Integer offset = luaL_optinteger(L, 4, 0);
    if (count < 0 || count > 1000000 || seed < 0 || offset < 0 ||
        (uint64_t)offset > UINT32_MAX ||
        (count && (uint64_t)(count - 1) > UINT32_MAX - (uint64_t)offset))
        return luaL_error(L, "invalid source preview range");
    lua_createtable(L, 0, 7);
    const char* names[] = {"position", "direction", "energy", "time",
                           "weight", "particle", "history_id"};
    for (int k = 0; k < 7; ++k) {
        lua_createtable(L, (int)count, 0);
        lua_setfield(L, -2, names[k]);
    }
    int result = lua_absindex(L, -1);
    for (lua_Integer i = 0; i < count; ++i) {
        alea_transport_source_t sample;
        alea_error_t err = alea_source_sample(source, (uint64_t)seed,
                                              (uint32_t)(offset + i), &sample);
        if (err != ALEA_OK)
            return luaL_error(L, "source sampling failed: %s", alea_error_string(err));
        lua_getfield(L, result, "position");
        push_vec3(L, sample.position); lua_rawseti(L, -2, i + 1); lua_pop(L, 1);
        lua_getfield(L, result, "direction");
        push_vec3(L, sample.particle.direction); lua_rawseti(L, -2, i + 1); lua_pop(L, 1);
        const double values[] = {sample.particle.energy, sample.particle.time,
                                 sample.particle.weight};
        for (int k = 0; k < 3; ++k) {
            lua_getfield(L, result, names[k + 2]);
            lua_pushnumber(L, values[k]); lua_rawseti(L, -2, i + 1); lua_pop(L, 1);
        }
        lua_getfield(L, result, "particle");
        lua_pushinteger(L, sample.particle.type); lua_rawseti(L, -2, i + 1); lua_pop(L, 1);
        lua_getfield(L, result, "history_id");
        lua_pushinteger(L, offset + i); lua_rawseti(L, -2, i + 1); lua_pop(L, 1);
    }
    return 1;
}

static void add_tally(lua_State* L, int idx, transport_guard_t* guard) {
    luaL_checktype(L, idx, LUA_TTABLE);
    alea_tally_spec_t spec = {0};
    spec.score = (alea_tally_score_t)enum_field(L, idx, "score", scores, 6, 0);
    spec.domain = (alea_tally_domain_t)enum_field(L, idx, "domain", domains, 3, 0);
    int particle = enum_field(L, idx, "particle", particles, 3, 2);
    spec.particle_mask = particle == 2 ? 0 : (1u << particle);
    vec3_field(L, idx, "lower", spec.lower, spec.domain == ALEA_TALLY_CARTESIAN_MESH);
    vec3_field(L, idx, "upper", spec.upper, spec.domain == ALEA_TALLY_CARTESIAN_MESH);
    if (spec.domain == ALEA_TALLY_CARTESIAN_MESH) {
        lua_getfield(L, idx, "dimensions");
        luaL_checktype(L, -1, LUA_TTABLE);
        if (lua_rawlen(L, -1) != 3) luaL_error(L, "dimensions must have three elements");
        for (int i = 0; i < 3; ++i) {
            lua_rawgeti(L, -1, i + 1);
            lua_Integer n = luaL_checkinteger(L, -1);
            if (n < 0 || (uint64_t)n > UINT32_MAX) luaL_error(L, "invalid mesh dimension");
            spec.dimensions[i] = (uint32_t)n;
            lua_pop(L, 1);
        }
        lua_pop(L, 1);
    }
    spec.energy_min = number_field(L, idx, "energy_min", 0);
    spec.energy_max = number_field(L, idx, "energy_max", 0);
    spec.time_min = number_field(L, idx, "time_min", 0);
    spec.time_max = number_field(L, idx, "time_max", 0);
    lua_Integer n = integer_field(L, idx, "material_id", 0);
    if (n < -1 || n > INT_MAX) luaL_error(L, "invalid material_id");
    spec.material_id = (int)n;
    n = integer_field(L, idx, "reaction_mt", 0);
    if (n < 0 || n > INT_MAX) luaL_error(L, "invalid reaction_mt");
    spec.reaction_mt = (int)n;
    n = integer_field(L, idx, "nuclide_zaid", 0);
    if (n < 0 || n > INT_MAX) luaL_error(L, "invalid nuclide_zaid");
    spec.nuclide_zaid = (int)n;
    lua_getfield(L, idx, "energy_edges");
    double* edges = NULL;
    if (!lua_isnil(L, -1)) {
        luaL_checktype(L, -1, LUA_TTABLE);
        size_t len = lua_rawlen(L, -1);
        if (len < 2) luaL_error(L, "energy_edges needs at least two entries");
        edges = malloc(len * sizeof(double));
        if (!edges) luaL_error(L, "out of memory");
        for (size_t i = 0; i < len; ++i) {
            lua_rawgeti(L, -1, (lua_Integer)i + 1);
            if (!lua_isnumber(L, -1)) {
                free(edges); luaL_error(L, "energy_edges must contain numbers");
            }
            edges[i] = lua_tonumber(L, -1);
            lua_pop(L, 1);
        }
        spec.energy_edges = edges;
        spec.energy_group_count = len - 1;
    }
    lua_pop(L, 1);
    alea_error_t err = alea_tally_plan_add(guard->plan, &spec, NULL);
    free(edges);
    if (err != ALEA_OK) luaL_error(L, "invalid tally: %s", alea_error_string(err));
}

static void push_result(lua_State* L, const alea_transport_result_t* result) {
    lua_createtable(L, 0, 17);
    set_integer(L, "histories", result->histories);
    set_integer(L, "absorbed", result->absorbed);
    set_integer(L, "replaced", result->replaced);
    set_integer(L, "leaked", result->leaked);
    set_integer(L, "collisions", result->collisions);
    set_integer(L, "boundary_crossings", result->boundary_crossings);
    set_integer(L, "reflections", result->reflections);
    set_integer(L, "emitted_neutrons", result->emitted_neutrons);
    set_integer(L, "emitted_photons", result->emitted_photons);
    set_integer(L, "photon_collisions", result->photon_collisions);
    set_integer(L, "photon_absorbed", result->photon_absorbed);
    set_integer(L, "photon_replaced", result->photon_replaced);
    set_integer(L, "photon_leaked", result->photon_leaked);
    push_floats(L, result->track_length, result->cell_count);
    lua_setfield(L, -2, "track_length");
    push_floats(L, result->track_length_squared, result->cell_count);
    lua_setfield(L, -2, "track_length_squared");
    size_t count = alea_tally_results_count(result->tallies);
    lua_createtable(L, (int)count, 0);
    for (size_t i = 0; i < count; ++i) {
        alea_tally_view_t v;
        if (alea_tally_results_view(result->tallies, i, &v) != ALEA_OK)
            luaL_error(L, "cannot read tally result");
        lua_createtable(L, 0, 19);
        lua_pushstring(L, scores[v.score]); lua_setfield(L, -2, "score");
        lua_pushstring(L, domains[v.domain]); lua_setfield(L, -2, "domain");
        lua_pushstring(L, v.particle_mask == ALEA_TALLY_NEUTRON ? "neutron" :
            v.particle_mask == ALEA_TALLY_PHOTON ? "photon" : "all");
        lua_setfield(L, -2, "particle");
        set_integer(L, "histories", v.histories);
        set_integer(L, "bin_count", (lua_Integer)v.bin_count);
        set_integer(L, "spatial_bin_count", (lua_Integer)v.spatial_bin_count);
        set_integer(L, "energy_group_count", (lua_Integer)v.energy_group_count);
        set_integer(L, "material_id", v.material_id);
        set_integer(L, "reaction_mt", v.reaction_mt);
        set_integer(L, "nuclide_zaid", v.nuclide_zaid);
        set_number(L, "energy_min", v.energy_min);
        set_number(L, "energy_max", v.energy_max);
        set_number(L, "time_min", v.time_min);
        set_number(L, "time_max", v.time_max);
        push_floats(L, v.sum, v.bin_count); lua_setfield(L, -2, "sum");
        push_floats(L, v.sum_squared, v.bin_count); lua_setfield(L, -2, "sum_squared");
        lua_createtable(L, (int)v.bin_count, 0);
        for (size_t k = 0; k < v.bin_count; ++k) {
            lua_pushnumber(L, v.sum[k] / v.histories);
            lua_rawseti(L, -2, (lua_Integer)k + 1);
        }
        lua_setfield(L, -2, "mean");
        lua_createtable(L, (int)v.bin_count, 0);
        for (size_t k = 0; k < v.bin_count; ++k) {
            double h = (double)v.histories;
            double variance = h > 1 ? (v.sum_squared[k] - v.sum[k] * v.sum[k] / h) / (h - 1) : 0;
            lua_pushnumber(L, sqrt(fmax(0.0, variance / h)));
            lua_rawseti(L, -2, (lua_Integer)k + 1);
        }
        lua_setfield(L, -2, "standard_error");
        if (v.bin_ids) {
            lua_createtable(L, (int)v.spatial_bin_count, 0);
            for (size_t k = 0; k < v.spatial_bin_count; ++k) {
                lua_pushinteger(L, v.bin_ids[k]); lua_rawseti(L, -2, (lua_Integer)k + 1);
            }
        } else lua_pushnil(L);
        lua_setfield(L, -2, "bin_ids");
        if (v.energy_edges) push_floats(L, v.energy_edges, v.energy_group_count + 1);
        else lua_pushnil(L);
        lua_setfield(L, -2, "energy_edges");
        push_vec3(L, v.lower); lua_setfield(L, -2, "lower");
        push_vec3(L, v.upper); lua_setfield(L, -2, "upper");
        lua_createtable(L, 3, 0);
        for (int k = 0; k < 3; ++k) { lua_pushinteger(L, v.dimensions[k]); lua_rawseti(L, -2, k + 1); }
        lua_setfield(L, -2, "dimensions");
        lua_rawseti(L, -2, (lua_Integer)i + 1);
    }
    lua_setfield(L, -2, "tallies");
}

static int transport_run(lua_State* L) {
    alea_system_t* sys = alea_get_sys(L, 1);
    alea_nuc_xsdir_t** xsdir = luaL_checkudata(L, 2, "alea.Xsdir");
    if (!*xsdir) return luaL_error(L, "xsdir has been destroyed");
    luaL_checktype(L, 3, LUA_TTABLE);
    transport_guard_t* guard = lua_newuserdatauv(L, sizeof(*guard), 0);
    memset(guard, 0, sizeof(*guard));
    luaL_setmetatable(L, TRANSPORT_GUARD_MT);
    alea_transport_options_t options = {.histories = 1, .seed = 1,
        .max_events_per_history = 100000, .max_segment_distance = 100.0};
    lua_Integer n = integer_field(L, 3, "histories", 1);
    if (n < 1 || (uint64_t)n > UINT32_MAX) return luaL_error(L, "histories must be 1..UINT32_MAX");
    options.histories = (uint32_t)n;
    n = integer_field(L, 3, "seed", 1);
    if (n < 0) return luaL_error(L, "seed must be nonnegative");
    options.seed = (uint64_t)n;
    n = integer_field(L, 3, "history_offset", 0);
    if (n < 0 || (uint64_t)n > UINT32_MAX) return luaL_error(L, "invalid history_offset");
    options.history_offset = (uint32_t)n;
    n = integer_field(L, 3, "max_events_per_history", 100000);
    if (n < 0 || (uint64_t)n > UINT32_MAX) return luaL_error(L, "invalid max_events_per_history");
    options.max_events_per_history = (uint32_t)n;
    n = integer_field(L, 3, "max_pending_particles", 0);
    if (n < 0) return luaL_error(L, "invalid max_pending_particles");
    options.max_pending_particles = (size_t)n;
    options.max_segment_distance = number_field(L, 3, "max_segment_distance", 100);
    lua_getfield(L, 3, "coupled");
    int coupled = lua_toboolean(L, -1); lua_pop(L, 1);
    lua_getfield(L, 3, "source");
    int source_idx = lua_absindex(L, -1);
    alea_source_t* source;
    if (luaL_testudata(L, source_idx, SOURCE_MT)) {
        source = *(alea_source_t**)luaL_checkudata(L, source_idx, SOURCE_MT);
        if (!source) return luaL_error(L, "source has been destroyed");
    } else source = guard->owned_source = source_from_table(L, source_idx);
    lua_pop(L, 1);
    lua_getfield(L, 3, "tallies");
    if (!lua_isnil(L, -1)) {
        luaL_checktype(L, -1, LUA_TTABLE);
        guard->plan = alea_tally_plan_create(sys);
        if (!guard->plan) return luaL_error(L, "cannot create tally plan");
        size_t count = lua_rawlen(L, -1);
        for (size_t i = 0; i < count; ++i) {
            lua_rawgeti(L, -1, (lua_Integer)i + 1);
            add_tally(L, lua_absindex(L, -1), guard);
            lua_pop(L, 1);
        }
    }
    lua_pop(L, 1);
    options.tally_plan = guard->plan;
    alea_nuc_prepare_requirements_t req = {.required_capabilities =
        ALEA_NUC_CAP_CONTINUOUS_NEUTRON |
        (coupled ? ALEA_NUC_CAP_PHOTON_PRODUCTION : 0)};
    uint32_t mask = alea_source_particle_mask(source) |
        (coupled ? ALEA_NUC_BIND_NEUTRON | ALEA_NUC_BIND_PHOTON : 0);
    alea_error_t err = alea_nuc_cell_bindings_prepare(sys, *xsdir, mask,
        mask & ALEA_NUC_BIND_NEUTRON ? &req : NULL, NULL, &guard->bindings);
    if (err != ALEA_OK)
        return luaL_error(L, "nuclear-data binding failed: %s", alea_error_string(err));
    alea_transport_failure_t failure = {0};
    err = alea_transport_run_sampled_source(sys, guard->bindings,
        alea_source_sample, source, &options, &guard->result, &failure);
    if (err != ALEA_OK)
        return luaL_error(L, "transport failed: %s (history %u, cell %d, position %.6g %.6g %.6g)",
            alea_error_string(err), failure.history_id, failure.cell_id,
            failure.position[0], failure.position[1], failure.position[2]);
    push_result(L, &guard->result);
    return 1;
}

int luaopen_alea_transport(lua_State* L) {
    luaL_newmetatable(L, SOURCE_BUFFERS_MT);
    lua_pushcfunction(L, source_buffers_gc); lua_setfield(L, -2, "__gc");
    lua_pop(L, 1);
    luaL_newmetatable(L, SOURCE_MT);
    lua_pushcfunction(L, source_gc); lua_setfield(L, -2, "__gc");
    lua_pop(L, 1);
    luaL_newmetatable(L, TRANSPORT_GUARD_MT);
    lua_pushcfunction(L, guard_gc); lua_setfield(L, -2, "__gc");
    lua_pop(L, 1);
    lua_pushcfunction(L, transport_run);
    lua_setfield(L, -2, "transport_run");
    lua_pushcfunction(L, source_prepare);
    lua_setfield(L, -2, "source_prepare");
    lua_pushcfunction(L, source_integrated_emissivity);
    lua_setfield(L, -2, "source_integrated_emissivity");
    lua_pushcfunction(L, sample_source);
    lua_setfield(L, -2, "sample_source");
    return 0;
}
