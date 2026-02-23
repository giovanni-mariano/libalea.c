// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

#include "alea_lua.h"
#include "alea_render.h"
#include <stdlib.h>
#include <string.h>

/* ============================================================================
 * Userdata
 * ============================================================================ */

typedef struct {
    render_framebuffer_t* fb;
    uint8_t* pixels;  /* tonemapped RGB pixels, allocated on first write */
} alea_lua_framebuffer_t;

static alea_lua_framebuffer_t* check_framebuffer(lua_State* L, int idx) {
    return (alea_lua_framebuffer_t*)luaL_checkudata(L, idx, ALEA_FRAMEBUFFER_MT);
}

/* ============================================================================
 * Config helper: Lua table -> render_config_t
 * ============================================================================ */

static void lua_to_render_config(lua_State* L, int idx, render_config_t* cfg) {
    luaL_checktype(L, idx, LUA_TTABLE);

#define GET_OPT_INT(field) do { \
    lua_getfield(L, idx, #field); \
    if (!lua_isnil(L, -1)) cfg->field = (int)lua_tointeger(L, -1); \
    lua_pop(L, 1); \
} while(0)

#define GET_OPT_NUM(field) do { \
    lua_getfield(L, idx, #field); \
    if (!lua_isnil(L, -1)) cfg->field = lua_tonumber(L, -1); \
    lua_pop(L, 1); \
} while(0)

#define GET_OPT_FLOAT(field) do { \
    lua_getfield(L, idx, #field); \
    if (!lua_isnil(L, -1)) cfg->field = (float)lua_tonumber(L, -1); \
    lua_pop(L, 1); \
} while(0)

    GET_OPT_INT(width);
    GET_OPT_INT(height);
    GET_OPT_NUM(fov);
    GET_OPT_NUM(ortho_height);
    GET_OPT_INT(shadows);
    GET_OPT_INT(edges);
    GET_OPT_FLOAT(ambient);
    GET_OPT_FLOAT(diffuse);
    GET_OPT_FLOAT(specular);
    GET_OPT_FLOAT(shininess);
    GET_OPT_FLOAT(cross_section_tint);
    GET_OPT_INT(aa_samples);
    GET_OPT_INT(aa_adaptive);
    GET_OPT_INT(universe_depth);
    GET_OPT_INT(threads);
    GET_OPT_INT(tile_size);
    GET_OPT_INT(preview);
    GET_OPT_INT(aux_output);
    GET_OPT_FLOAT(xray_density_scale);
    GET_OPT_INT(dedup);
    GET_OPT_INT(log_level);

    /* color_mode / render_mode as integers */
    lua_getfield(L, idx, "color_mode");
    if (!lua_isnil(L, -1)) cfg->color_mode = (render_color_mode_t)(int)lua_tointeger(L, -1);
    lua_pop(L, 1);

    lua_getfield(L, idx, "render_mode");
    if (!lua_isnil(L, -1)) cfg->render_mode = (render_mode_t)(int)lua_tointeger(L, -1);
    lua_pop(L, 1);

    /* eye = {x, y, z} */
    lua_getfield(L, idx, "eye");
    if (lua_istable(L, -1)) {
        lua_rawgeti(L, -1, 1); cfg->eye[0] = lua_tonumber(L, -1); lua_pop(L, 1);
        lua_rawgeti(L, -1, 2); cfg->eye[1] = lua_tonumber(L, -1); lua_pop(L, 1);
        lua_rawgeti(L, -1, 3); cfg->eye[2] = lua_tonumber(L, -1); lua_pop(L, 1);
        cfg->eye_set = 1;
    }
    lua_pop(L, 1);

    /* target = {x, y, z} */
    lua_getfield(L, idx, "target");
    if (lua_istable(L, -1)) {
        lua_rawgeti(L, -1, 1); cfg->target[0] = lua_tonumber(L, -1); lua_pop(L, 1);
        lua_rawgeti(L, -1, 2); cfg->target[1] = lua_tonumber(L, -1); lua_pop(L, 1);
        lua_rawgeti(L, -1, 3); cfg->target[2] = lua_tonumber(L, -1); lua_pop(L, 1);
        cfg->target_set = 1;
    }
    lua_pop(L, 1);

    /* up = {x, y, z} */
    lua_getfield(L, idx, "up");
    if (lua_istable(L, -1)) {
        lua_rawgeti(L, -1, 1); cfg->up[0] = lua_tonumber(L, -1); lua_pop(L, 1);
        lua_rawgeti(L, -1, 2); cfg->up[1] = lua_tonumber(L, -1); lua_pop(L, 1);
        lua_rawgeti(L, -1, 3); cfg->up[2] = lua_tonumber(L, -1); lua_pop(L, 1);
    }
    lua_pop(L, 1);

    /* light_dir = {x, y, z} */
    lua_getfield(L, idx, "light_dir");
    if (lua_istable(L, -1)) {
        lua_rawgeti(L, -1, 1); cfg->light_dir[0] = lua_tonumber(L, -1); lua_pop(L, 1);
        lua_rawgeti(L, -1, 2); cfg->light_dir[1] = lua_tonumber(L, -1); lua_pop(L, 1);
        lua_rawgeti(L, -1, 3); cfg->light_dir[2] = lua_tonumber(L, -1); lua_pop(L, 1);
        cfg->light_set = 1;
    }
    lua_pop(L, 1);

    /* background = {r, g, b} */
    lua_getfield(L, idx, "background");
    if (lua_istable(L, -1)) {
        lua_rawgeti(L, -1, 1); cfg->background[0] = (float)lua_tonumber(L, -1); lua_pop(L, 1);
        lua_rawgeti(L, -1, 2); cfg->background[1] = (float)lua_tonumber(L, -1); lua_pop(L, 1);
        lua_rawgeti(L, -1, 3); cfg->background[2] = (float)lua_tonumber(L, -1); lua_pop(L, 1);
    }
    lua_pop(L, 1);

    /* clips = {{nx, ny, nz, d}, ...} */
    lua_getfield(L, idx, "clips");
    if (lua_istable(L, -1)) {
        int n = (int)luaL_len(L, -1);
        if (n > RENDER_MAX_CLIPS) n = RENDER_MAX_CLIPS;
        cfg->num_clips = n;
        for (int i = 0; i < n; i++) {
            lua_rawgeti(L, -1, i + 1);
            if (lua_istable(L, -1)) {
                lua_rawgeti(L, -1, 1); cfg->clips[i].normal[0] = lua_tonumber(L, -1); lua_pop(L, 1);
                lua_rawgeti(L, -1, 2); cfg->clips[i].normal[1] = lua_tonumber(L, -1); lua_pop(L, 1);
                lua_rawgeti(L, -1, 3); cfg->clips[i].normal[2] = lua_tonumber(L, -1); lua_pop(L, 1);
                lua_rawgeti(L, -1, 4); cfg->clips[i].d = lua_tonumber(L, -1); lua_pop(L, 1);
            }
            lua_pop(L, 1);
        }
    }
    lua_pop(L, 1);

#undef GET_OPT_INT
#undef GET_OPT_NUM
#undef GET_OPT_FLOAT
}

/* Push render_config_t as Lua table */
static void push_render_config(lua_State* L, const render_config_t* cfg) {
    lua_createtable(L, 0, 30);
    lua_pushinteger(L, cfg->width);          lua_setfield(L, -2, "width");
    lua_pushinteger(L, cfg->height);         lua_setfield(L, -2, "height");
    lua_pushnumber(L, cfg->fov);             lua_setfield(L, -2, "fov");
    lua_pushnumber(L, cfg->ortho_height);    lua_setfield(L, -2, "ortho_height");
    lua_pushinteger(L, cfg->shadows);        lua_setfield(L, -2, "shadows");
    lua_pushinteger(L, cfg->edges);          lua_setfield(L, -2, "edges");
    lua_pushnumber(L, cfg->ambient);         lua_setfield(L, -2, "ambient");
    lua_pushnumber(L, cfg->diffuse);         lua_setfield(L, -2, "diffuse");
    lua_pushnumber(L, cfg->specular);        lua_setfield(L, -2, "specular");
    lua_pushnumber(L, cfg->shininess);       lua_setfield(L, -2, "shininess");
    lua_pushnumber(L, cfg->cross_section_tint); lua_setfield(L, -2, "cross_section_tint");
    lua_pushinteger(L, cfg->aa_samples);     lua_setfield(L, -2, "aa_samples");
    lua_pushinteger(L, cfg->aa_adaptive);    lua_setfield(L, -2, "aa_adaptive");
    lua_pushinteger(L, cfg->universe_depth); lua_setfield(L, -2, "universe_depth");
    lua_pushinteger(L, cfg->threads);        lua_setfield(L, -2, "threads");
    lua_pushinteger(L, cfg->tile_size);      lua_setfield(L, -2, "tile_size");
    lua_pushinteger(L, cfg->preview);        lua_setfield(L, -2, "preview");
    lua_pushinteger(L, cfg->aux_output);     lua_setfield(L, -2, "aux_output");
    lua_pushinteger(L, cfg->color_mode);     lua_setfield(L, -2, "color_mode");
    lua_pushinteger(L, cfg->render_mode);    lua_setfield(L, -2, "render_mode");
    lua_pushnumber(L, cfg->xray_density_scale); lua_setfield(L, -2, "xray_density_scale");
    lua_pushinteger(L, cfg->dedup);          lua_setfield(L, -2, "dedup");
    lua_pushinteger(L, cfg->log_level);      lua_setfield(L, -2, "log_level");

    /* eye */
    lua_createtable(L, 3, 0);
    lua_pushnumber(L, cfg->eye[0]); lua_rawseti(L, -2, 1);
    lua_pushnumber(L, cfg->eye[1]); lua_rawseti(L, -2, 2);
    lua_pushnumber(L, cfg->eye[2]); lua_rawseti(L, -2, 3);
    lua_setfield(L, -2, "eye");

    /* target */
    lua_createtable(L, 3, 0);
    lua_pushnumber(L, cfg->target[0]); lua_rawseti(L, -2, 1);
    lua_pushnumber(L, cfg->target[1]); lua_rawseti(L, -2, 2);
    lua_pushnumber(L, cfg->target[2]); lua_rawseti(L, -2, 3);
    lua_setfield(L, -2, "target");

    /* up */
    lua_createtable(L, 3, 0);
    lua_pushnumber(L, cfg->up[0]); lua_rawseti(L, -2, 1);
    lua_pushnumber(L, cfg->up[1]); lua_rawseti(L, -2, 2);
    lua_pushnumber(L, cfg->up[2]); lua_rawseti(L, -2, 3);
    lua_setfield(L, -2, "up");

    /* background */
    lua_createtable(L, 3, 0);
    lua_pushnumber(L, cfg->background[0]); lua_rawseti(L, -2, 1);
    lua_pushnumber(L, cfg->background[1]); lua_rawseti(L, -2, 2);
    lua_pushnumber(L, cfg->background[2]); lua_rawseti(L, -2, 3);
    lua_setfield(L, -2, "background");
}

/* ============================================================================
 * Module-level functions
 * ============================================================================ */

/* alea.render_config() -> default config table */
static int l_render_config(lua_State* L) {
    render_config_t cfg;
    render_config_init(&cfg);
    push_render_config(L, &cfg);
    return 1;
}

/* alea.render_load_colors(file, cfg_table) -> updated cfg_table */
static int l_render_load_colors(lua_State* L) {
    const char* filename = luaL_checkstring(L, 1);
    render_config_t cfg;
    render_config_init(&cfg);
    if (lua_istable(L, 2))
        lua_to_render_config(L, 2, &cfg);

    if (render_load_colors(filename, &cfg) != 0) {
        render_config_free(&cfg);
        return luaL_error(L, "render_load_colors failed: %s", alea_error());
    }
    push_render_config(L, &cfg);
    render_config_free(&cfg);
    return 1;
}

/* ============================================================================
 * System methods
 * ============================================================================ */

/* Ensure pixels buffer is allocated and tonemapped */
static void ensure_pixels(alea_lua_framebuffer_t* ud) {
    if (!ud->pixels && ud->fb) {
        size_t sz = (size_t)ud->fb->width * (size_t)ud->fb->height * 3;
        ud->pixels = (uint8_t*)malloc(sz);
        if (ud->pixels)
            render_tonemap(ud->fb, ud->pixels);
    }
}

/* sys:render(config) -> Framebuffer */
static int l_render(lua_State* L) {
    alea_system_t* sys = alea_get_sys(L, 1);
    render_config_t cfg;
    render_config_init(&cfg);
    if (lua_istable(L, 2))
        lua_to_render_config(L, 2, &cfg);

    render_camera_t cam;
    if (render_camera_setup(&cam, &cfg, sys) != 0) {
        render_config_free(&cfg);
        return luaL_error(L, "render_camera_setup failed: %s", alea_error());
    }

    alea_lua_framebuffer_t* ud = (alea_lua_framebuffer_t*)lua_newuserdata(
        L, sizeof(alea_lua_framebuffer_t));
    ud->fb = NULL;
    ud->pixels = NULL;
    luaL_setmetatable(L, ALEA_FRAMEBUFFER_MT);

    ud->fb = render_framebuffer_create(cfg.width, cfg.height, cfg.aux_output);
    if (!ud->fb) {
        render_config_free(&cfg);
        return luaL_error(L, "render_framebuffer_create failed");
    }

    int rc = render_scene(sys, &cfg, &cam, ud->fb);
    render_config_free(&cfg);
    if (rc != 0)
        return luaL_error(L, "render_scene failed: %s", alea_error());

    return 1;
}

/* sys:render_camera_setup(config) -> camera table */
static int l_render_camera_setup(lua_State* L) {
    alea_system_t* sys = alea_get_sys(L, 1);
    render_config_t cfg;
    render_config_init(&cfg);
    if (lua_istable(L, 2))
        lua_to_render_config(L, 2, &cfg);

    render_camera_t cam;
    int rc = render_camera_setup(&cam, &cfg, sys);
    render_config_free(&cfg);
    if (rc != 0)
        return luaL_error(L, "render_camera_setup failed: %s", alea_error());

    lua_createtable(L, 0, 10);

    lua_createtable(L, 3, 0);
    lua_pushnumber(L, cam.eye[0]); lua_rawseti(L, -2, 1);
    lua_pushnumber(L, cam.eye[1]); lua_rawseti(L, -2, 2);
    lua_pushnumber(L, cam.eye[2]); lua_rawseti(L, -2, 3);
    lua_setfield(L, -2, "eye");

    lua_createtable(L, 3, 0);
    lua_pushnumber(L, cam.target[0]); lua_rawseti(L, -2, 1);
    lua_pushnumber(L, cam.target[1]); lua_rawseti(L, -2, 2);
    lua_pushnumber(L, cam.target[2]); lua_rawseti(L, -2, 3);
    lua_setfield(L, -2, "target");

    lua_createtable(L, 3, 0);
    lua_pushnumber(L, cam.forward[0]); lua_rawseti(L, -2, 1);
    lua_pushnumber(L, cam.forward[1]); lua_rawseti(L, -2, 2);
    lua_pushnumber(L, cam.forward[2]); lua_rawseti(L, -2, 3);
    lua_setfield(L, -2, "forward");

    lua_pushnumber(L, cam.fov_degrees);   lua_setfield(L, -2, "fov_degrees");
    lua_pushnumber(L, cam.ortho_height);  lua_setfield(L, -2, "ortho_height");
    lua_pushboolean(L, cam.auto_fit);     lua_setfield(L, -2, "auto_fit");

    return 1;
}

/* ============================================================================
 * Framebuffer methods
 * ============================================================================ */

/* fb:width() */
static int l_fb_width(lua_State* L) {
    alea_lua_framebuffer_t* ud = check_framebuffer(L, 1);
    if (!ud->fb) return luaL_error(L, "framebuffer freed");
    lua_pushinteger(L, ud->fb->width);
    return 1;
}

/* fb:height() */
static int l_fb_height(lua_State* L) {
    alea_lua_framebuffer_t* ud = check_framebuffer(L, 1);
    if (!ud->fb) return luaL_error(L, "framebuffer freed");
    lua_pushinteger(L, ud->fb->height);
    return 1;
}

/* fb:edge_darken() */
static int l_fb_edge_darken(lua_State* L) {
    alea_lua_framebuffer_t* ud = check_framebuffer(L, 1);
    if (!ud->fb) return luaL_error(L, "framebuffer freed");
    render_edge_darken(ud->fb);
    /* Invalidate cached pixels */
    if (ud->pixels) { free(ud->pixels); ud->pixels = NULL; }
    return 0;
}

/* fb:write(filename) */
static int l_fb_write(lua_State* L) {
    alea_lua_framebuffer_t* ud = check_framebuffer(L, 1);
    if (!ud->fb) return luaL_error(L, "framebuffer freed");
    const char* filename = luaL_checkstring(L, 2);
    ensure_pixels(ud);
    if (!ud->pixels) return luaL_error(L, "tonemap failed: out of memory");
    if (render_write_image(filename, ud->pixels, ud->fb->width, ud->fb->height) != 0)
        return luaL_error(L, "write failed: %s", alea_error());
    return 0;
}

/* fb:write_png(filename) */
static int l_fb_write_png(lua_State* L) {
    alea_lua_framebuffer_t* ud = check_framebuffer(L, 1);
    if (!ud->fb) return luaL_error(L, "framebuffer freed");
    const char* filename = luaL_checkstring(L, 2);
    ensure_pixels(ud);
    if (!ud->pixels) return luaL_error(L, "tonemap failed: out of memory");
    if (render_write_png(filename, ud->pixels, ud->fb->width, ud->fb->height) != 0)
        return luaL_error(L, "write_png failed: %s", alea_error());
    return 0;
}

/* fb:write_bmp(filename) */
static int l_fb_write_bmp(lua_State* L) {
    alea_lua_framebuffer_t* ud = check_framebuffer(L, 1);
    if (!ud->fb) return luaL_error(L, "framebuffer freed");
    const char* filename = luaL_checkstring(L, 2);
    ensure_pixels(ud);
    if (!ud->pixels) return luaL_error(L, "tonemap failed: out of memory");
    if (render_write_bmp(filename, ud->pixels, ud->fb->width, ud->fb->height) != 0)
        return luaL_error(L, "write_bmp failed: %s", alea_error());
    return 0;
}

/* fb:write_ppm(filename) */
static int l_fb_write_ppm(lua_State* L) {
    alea_lua_framebuffer_t* ud = check_framebuffer(L, 1);
    if (!ud->fb) return luaL_error(L, "framebuffer freed");
    const char* filename = luaL_checkstring(L, 2);
    ensure_pixels(ud);
    if (!ud->pixels) return luaL_error(L, "tonemap failed: out of memory");
    if (render_write_ppm(filename, ud->pixels, ud->fb->width, ud->fb->height) != 0)
        return luaL_error(L, "write_ppm failed: %s", alea_error());
    return 0;
}

/* fb:write_aux(base) */
static int l_fb_write_aux(lua_State* L) {
    alea_lua_framebuffer_t* ud = check_framebuffer(L, 1);
    if (!ud->fb) return luaL_error(L, "framebuffer freed");
    const char* base = luaL_checkstring(L, 2);
    if (render_write_aux(base, ud->fb) != 0)
        return luaL_error(L, "write_aux failed: %s", alea_error());
    return 0;
}

/* __gc */
static int l_fb_gc(lua_State* L) {
    alea_lua_framebuffer_t* ud = check_framebuffer(L, 1);
    if (ud->pixels) { free(ud->pixels); ud->pixels = NULL; }
    if (ud->fb) { render_framebuffer_free(ud->fb); ud->fb = NULL; }
    return 0;
}

/* ============================================================================
 * Registration
 * ============================================================================ */

static const luaL_Reg render_system_methods[] = {
    {"render",              l_render},
    {"render_camera_setup", l_render_camera_setup},
    {NULL, NULL}
};

static const luaL_Reg fb_meta[] = {
    {"__gc", l_fb_gc},
    {NULL, NULL}
};

static const luaL_Reg fb_methods[] = {
    {"width",       l_fb_width},
    {"height",      l_fb_height},
    {"edge_darken", l_fb_edge_darken},
    {"write",       l_fb_write},
    {"write_png",   l_fb_write_png},
    {"write_bmp",   l_fb_write_bmp},
    {"write_ppm",   l_fb_write_ppm},
    {"write_aux",   l_fb_write_aux},
    {NULL, NULL}
};

int luaopen_alea_render(lua_State* L) {
    /* Add system methods */
    luaL_getmetatable(L, ALEA_SYSTEM_MT);
    lua_getfield(L, -1, "__index");
    luaL_setfuncs(L, render_system_methods, 0);
    lua_pop(L, 2);

    /* Create Framebuffer metatable */
    luaL_newmetatable(L, ALEA_FRAMEBUFFER_MT);
    luaL_setfuncs(L, fb_meta, 0);
    lua_newtable(L);
    luaL_setfuncs(L, fb_methods, 0);
    lua_setfield(L, -2, "__index");
    lua_pop(L, 1);

    /* Module-level functions */
    lua_pushcfunction(L, l_render_config);
    lua_setfield(L, -2, "render_config");
    lua_pushcfunction(L, l_render_load_colors);
    lua_setfield(L, -2, "render_load_colors");

    return 0;
}
