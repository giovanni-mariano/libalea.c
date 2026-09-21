// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

#include "alea_lua.h"
#include "alea_mcnp.h"
#include "alea_openmc.h"
#include "alea_serpent.h"
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static FILE* new_export_stream(lua_State* L) {
    FILE* stream = tmpfile();
    if (!stream) luaL_error(L, "could not create temporary export stream: %s",
                            strerror(errno));
    return stream;
}

/* A FreeGuard must already be on the stack so allocation remains protected if
 * lua_pushlstring raises an out-of-memory error. */
static int push_exported_string(lua_State* L, FILE* stream,
                                alea_lua_free_guard_t* guard,
                                const char* operation) {
    if (fflush(stream) != 0 || fseek(stream, 0, SEEK_END) != 0) {
        int error_number = errno;
        fclose(stream);
        return luaL_error(L, "%s: %s", operation, strerror(error_number));
    }
    long length = ftell(stream);
    if (length < 0 || fseek(stream, 0, SEEK_SET) != 0) {
        int error_number = errno;
        fclose(stream);
        return luaL_error(L, "%s: %s", operation, strerror(error_number));
    }
    if ((unsigned long)length > (unsigned long)SIZE_MAX - 1) {
        fclose(stream);
        return luaL_error(L, "%s: output is too large", operation);
    }
    guard->ptr = malloc((size_t)length + 1);
    if (!guard->ptr) {
        fclose(stream);
        return luaL_error(L, "%s: out of memory", operation);
    }
    size_t got = fread(guard->ptr, 1, (size_t)length, stream);
    int read_error = ferror(stream);
    fclose(stream);
    if (got != (size_t)length || read_error)
        return luaL_error(L, "%s: could not read temporary export", operation);
    ((char*)guard->ptr)[length] = '\0';
    lua_pushlstring(L, (const char*)guard->ptr, (size_t)length);
    free(guard->ptr);
    guard->ptr = NULL;
    lua_remove(L, -2); /* remove guard and leave the string */
    return 1;
}

/* ============================================================================
 * Load functions (return new owned System userdata)
 * ============================================================================ */

/* alea.load_mcnp(filename) -> System */
static int l_load_mcnp(lua_State* L) {
    const char* filename = luaL_checkstring(L, 1);
    alea_lua_system_t* ud = (alea_lua_system_t*)lua_newuserdata(L, sizeof(alea_lua_system_t));
    alea_lua_system_init(ud);
    luaL_setmetatable(L, ALEA_SYSTEM_MT);

    mcnp_model_t* model = mcnp_load(filename);
    if (!model)
        return luaL_error(L, "load_mcnp failed: %s", alea_error());
    ud->mcnp_model = model;
    ud->sys = model->sys;
    return 1;
}

/* alea.load_mcnp_string(str) -> System */
static int l_load_mcnp_string(lua_State* L) {
    size_t len;
    const char* str = luaL_checklstring(L, 1, &len);
    alea_lua_system_t* ud = (alea_lua_system_t*)lua_newuserdata(L, sizeof(alea_lua_system_t));
    alea_lua_system_init(ud);
    luaL_setmetatable(L, ALEA_SYSTEM_MT);

    mcnp_model_t* model = mcnp_load_string(str, len);
    if (!model)
        return luaL_error(L, "load_mcnp_string failed: %s", alea_error());
    ud->mcnp_model = model;
    ud->sys = model->sys;
    return 1;
}

/* alea.load_openmc(filename) -> System */
static int l_load_openmc(lua_State* L) {
    const char* filename = luaL_checkstring(L, 1);
    alea_lua_system_t* ud = (alea_lua_system_t*)lua_newuserdata(L, sizeof(alea_lua_system_t));
    alea_lua_system_init(ud);
    luaL_setmetatable(L, ALEA_SYSTEM_MT);

    openmc_model_t* model = openmc_load(filename);
    if (!model)
        return luaL_error(L, "load_openmc failed: %s", alea_error());
    ud->openmc_model = model;
    ud->sys = model->sys;
    return 1;
}

/* alea.load_openmc_string(str) -> System */
static int l_load_openmc_string(lua_State* L) {
    size_t len;
    const char* str = luaL_checklstring(L, 1, &len);
    alea_lua_system_t* ud =
        (alea_lua_system_t*)lua_newuserdata(L, sizeof(alea_lua_system_t));
    alea_lua_system_init(ud);
    luaL_setmetatable(L, ALEA_SYSTEM_MT);

    openmc_model_t* model = openmc_load_string(str, len);
    if (!model)
        return luaL_error(L, "load_openmc_string failed: %s", alea_error());
    ud->openmc_model = model;
    ud->sys = openmc_model_system(model);
    if (!ud->sys)
        return luaL_error(L, "load_openmc_string failed: model has no geometry system");
    return 1;
}

/* ============================================================================
 * Export functions (System methods)
 * ============================================================================ */

/* sys:export_mcnp(filename) */
static int l_export_mcnp(lua_State* L) {
    alea_lua_system_t* ud = alea_check_system(L, 1);
    alea_system_t* sys = alea_get_sys(L, 1);
    const char* filename = luaL_checkstring(L, 2);

    int rc;
    if (ud->mcnp_model) {
        rc = mcnp_export((const mcnp_model_t*)ud->mcnp_model, filename);
    } else {
        rc = mcnp_export_system(sys, filename);
    }
    if (rc != 0)
        return luaL_error(L, "export_mcnp failed: %s", alea_error());
    return 0;
}

/* sys:export_openmc(filename) */
static int l_export_openmc(lua_State* L) {
    alea_lua_system_t* ud = alea_check_system(L, 1);
    alea_system_t* sys = alea_get_sys(L, 1);
    const char* filename = luaL_checkstring(L, 2);

    int rc;
    if (ud->openmc_model) {
        rc = openmc_export((const openmc_model_t*)ud->openmc_model, filename);
    } else {
        rc = openmc_export_system(sys, filename);
    }
    if (rc != 0)
        return luaL_error(L, "export_openmc failed: %s", alea_error());
    return 0;
}

/* sys:export_serpent(filename) */
static int l_export_serpent(lua_State* L) {
    alea_system_t* sys = alea_get_sys(L, 1);
    const char* filename = luaL_checkstring(L, 2);

    int rc = serpent_export_system(sys, filename);
    if (rc != 0)
        return luaL_error(L, "export_serpent failed: %s", alea_error());
    return 0;
}

static int l_export_mcnp_string(lua_State* L) {
    alea_lua_system_t* ud = alea_check_system(L, 1);
    alea_system_t* sys = alea_get_sys(L, 1);
    alea_lua_free_guard_t* guard = alea_lua_push_free_guard(L);
    FILE* stream = new_export_stream(L);
    int rc = ud->mcnp_model
        ? mcnp_export_stream((const mcnp_model_t*)ud->mcnp_model, stream)
        : mcnp_export_system_stream(sys, stream);
    if (rc != 0) {
        fclose(stream);
        return luaL_error(L, "export_mcnp_string failed: %s", alea_error());
    }
    return push_exported_string(L, stream, guard, "export_mcnp_string");
}

static int l_export_openmc_string(lua_State* L) {
    alea_system_t* sys = alea_get_sys(L, 1);
    alea_lua_free_guard_t* guard = alea_lua_push_free_guard(L);
    FILE* stream = new_export_stream(L);
    if (openmc_export_system_stream(sys, stream) != 0) {
        fclose(stream);
        return luaL_error(L, "export_openmc_string failed: %s", alea_error());
    }
    return push_exported_string(L, stream, guard, "export_openmc_string");
}

static int l_export_serpent_string(lua_State* L) {
    alea_system_t* sys = alea_get_sys(L, 1);
    alea_lua_free_guard_t* guard = alea_lua_push_free_guard(L);
    FILE* stream = new_export_stream(L);
    if (serpent_export_system_stream(sys, stream) != 0) {
        fclose(stream);
        return luaL_error(L, "export_serpent_string failed: %s", alea_error());
    }
    return push_exported_string(L, stream, guard, "export_serpent_string");
}

/* ============================================================================
 * Registration
 * ============================================================================ */

static const luaL_Reg system_io_methods[] = {
    {"export_mcnp",  l_export_mcnp},
    {"export_openmc", l_export_openmc},
    {"export_serpent", l_export_serpent},
    {"export_mcnp_string", l_export_mcnp_string},
    {"export_openmc_string", l_export_openmc_string},
    {"export_serpent_string", l_export_serpent_string},
    {NULL, NULL}
};

int luaopen_alea_io(lua_State* L) {
    /* Add export methods to System metatable's __index */
    luaL_getmetatable(L, ALEA_SYSTEM_MT);
    lua_getfield(L, -1, "__index");
    luaL_setfuncs(L, system_io_methods, 0);
    lua_pop(L, 2);

    /* Add load functions to alea table (top of stack) */
    lua_pushcfunction(L, l_load_mcnp);
    lua_setfield(L, -2, "load_mcnp");
    lua_pushcfunction(L, l_load_mcnp_string);
    lua_setfield(L, -2, "load_mcnp_string");
    lua_pushcfunction(L, l_load_openmc);
    lua_setfield(L, -2, "load_openmc");
    lua_pushcfunction(L, l_load_openmc_string);
    lua_setfield(L, -2, "load_openmc_string");

    return 0;
}
