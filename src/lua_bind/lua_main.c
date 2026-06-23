// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <signal.h>

#include "util/compat.h"
#define IS_TTY() alea_file_is_tty(stdin)

#include "lua.h"
#include "lauxlib.h"
#include "lualib.h"
#include "alea.h"

#ifndef _WIN32
#include "linenoise.h"
#endif

/* Forward declaration */
int luaopen_alea(lua_State* L);

/* ============================================================================
 * SIGINT handling
 * ============================================================================ */

static lua_State* g_lua_state = NULL;

static void sigint_hook(lua_State* L, lua_Debug* ar) {
    (void)ar;
    lua_sethook(L, NULL, 0, 0);
    luaL_error(L, "interrupted!");
}

static void sigint_handler(int sig) {
    (void)sig;
    alea_interrupt();
    if (g_lua_state)
        lua_sethook(g_lua_state, sigint_hook, LUA_MASKCALL | LUA_MASKRET | LUA_MASKCOUNT, 1);
}

/* ============================================================================
 * Usage
 * ============================================================================ */

static void print_usage(const char* prog) {
    fprintf(stderr,
        "Usage: %s [options] [script.lua [args...]]\n"
        "       %s                    (interactive REPL)\n"
        "\n"
        "Options:\n"
        "  -e 'chunk'     Execute string\n"
        "  --set key=val  Set alea.params.key = val\n"
        "  -l level       Set log level (0-4)\n"
        "  -v             Print version\n"
        "  -h             Print this help\n",
        prog, prog);
}

/* ============================================================================
 * Param injection: --set key=val -> alea.params.key
 * ============================================================================ */

static void inject_param(lua_State* L, const char* kv) {
    lua_getglobal(L, "alea");
    lua_getfield(L, -1, "params");
    if (lua_isnil(L, -1)) {
        lua_pop(L, 1);
        lua_newtable(L);
        lua_pushvalue(L, -1);
        lua_setfield(L, -3, "params");
    }

    const char* eq = strchr(kv, '=');
    if (!eq) {
        fprintf(stderr, "Warning: --set requires key=value, got '%s'\n", kv);
        lua_pop(L, 2);
        return;
    }

    lua_pushlstring(L, kv, (size_t)(eq - kv));

    const char* val = eq + 1;
    char* end;
    double num = strtod(val, &end);
    if (end != val && *end == '\0') {
        lua_pushnumber(L, num);
    } else {
        lua_pushstring(L, val);
    }
    lua_settable(L, -3);
    lua_pop(L, 2);
}

/* ============================================================================
 * Arg table: alea.arg
 * ============================================================================ */

static void create_arg_table(lua_State* L, int argc, char** argv, int script_idx) {
    lua_newtable(L);
    if (script_idx < argc) {
        lua_pushstring(L, argv[script_idx]);
        lua_rawseti(L, -2, 0);
    }
    for (int i = script_idx + 1; i < argc; i++) {
        lua_pushstring(L, argv[i]);
        lua_rawseti(L, -2, i - script_idx);
    }

    /* Set as global 'arg' (standard Lua convention) */
    lua_pushvalue(L, -1);
    lua_setglobal(L, "arg");

    /* Also set as alea.arg */
    lua_getglobal(L, "alea");
    lua_pushvalue(L, -2);
    lua_setfield(L, -2, "arg");
    lua_pop(L, 2); /* pop alea table and original table */
}

/* ============================================================================
 * Error reporting
 * ============================================================================ */

static int msghandler(lua_State* L) {
    const char* msg = lua_tostring(L, 1);
    if (!msg) {
        if (luaL_callmeta(L, 1, "__tostring") && lua_type(L, -1) == LUA_TSTRING)
            return 1;
        msg = lua_pushfstring(L, "(error object is a %s value)", luaL_typename(L, 1));
    }
    luaL_traceback(L, L, msg, 1);
    return 1;
}

static int docall(lua_State* L, int narg, int nres) {
    int base = lua_gettop(L) - narg;
    lua_pushcfunction(L, msghandler);
    lua_insert(L, base);
    int status = lua_pcall(L, narg, nres, base);
    lua_remove(L, base);
    return status;
}

/* ============================================================================
 * REPL helpers
 * ============================================================================ */

/* Check if a load error is an incomplete-input (<eof>) error */
static int incomplete(lua_State* L, int status) {
    if (status != LUA_ERRSYNTAX)
        return 0;
    const char* msg = lua_tostring(L, -1);
    if (!msg)
        return 0;
    size_t len = strlen(msg);
    if (len >= 5 && strcmp(msg + len - 5, "<eof>") == 0)
        return 1;
    return 0;
}

/* Append src to *buf (realloc as needed). Caller frees *buf. */
static void append(char** buf, size_t* len, size_t* cap, const char* src) {
    size_t slen = strlen(src);
    while (*len + slen + 1 > *cap) {
        *cap = *cap ? *cap * 2 : 256;
        *buf = (char*)realloc(*buf, *cap);
    }
    memcpy(*buf + *len, src, slen);
    *len += slen;
    (*buf)[*len] = '\0';
}

/* Try to compile and execute a complete input line/block */
static void repl_execute(lua_State* L, const char* source, int first_line) {
    int top = lua_gettop(L);

    int status;
    if (first_line) {
        /* Try as expression "return <line>" */
        lua_pushfstring(L, "return %s", source);
        const char* expr = lua_tostring(L, -1);
        status = luaL_loadbuffer(L, expr, strlen(expr), "=stdin");
        lua_remove(L, -2);
        if (status != LUA_OK) {
            lua_pop(L, 1);
            status = luaL_loadbuffer(L, source, strlen(source), "=stdin");
        }
    } else {
        status = luaL_loadbuffer(L, source, strlen(source), "=stdin");
    }

    if (status == LUA_OK) {
        int base = lua_gettop(L) - 1;
        status = docall(L, 0, LUA_MULTRET);
        if (status == LUA_OK) {
            int nresults = lua_gettop(L) - base;
            if (nresults > 0) {
                luaL_checkstack(L, LUA_MINSTACK, "too many results");
                lua_getglobal(L, "print");
                lua_insert(L, base + 1);
                if (lua_pcall(L, nresults, 0, 0) != LUA_OK) {
                    fprintf(stderr, "%s\n", lua_tostring(L, -1));
                    lua_pop(L, 1);
                }
            }
        }
    }

    if (status != LUA_OK) {
        const char* msg = lua_tostring(L, -1);
        if (msg)
            fprintf(stderr, "%s\n", msg);
        lua_pop(L, 1);
    }

    lua_settop(L, top);
    alea_clear_interrupt();
}

/* ============================================================================
 * REPL — linenoise (Unix) / fgets fallback (Windows)
 *
 * - Multi-line: on <eof> syntax error, continue with ">> " prompt
 * - Expression mode: try "return <line>" first, then raw <line>
 * - Ctrl-C breaks running Lua; Ctrl-D exits
 * - Print non-nil results
 * ============================================================================ */

#ifndef _WIN32
/* ---- Unix REPL: linenoise with history ---- */

#define HISTORY_FILE ".alea_history"

static const char* get_history_path(void) {
    static char path[1024];
    const char* home = getenv("HOME");
    if (home) {
        snprintf(path, sizeof(path), "%s/%s", home, HISTORY_FILE);
        return path;
    }
    return HISTORY_FILE;
}

static void repl(lua_State* L) {
    const char* hist = get_history_path();
    linenoiseHistorySetMaxLen(1000);
    linenoiseHistoryLoad(hist);

    printf("alea %s (Lua 5.5) -- type Ctrl-D to exit\n", alea_version());

    char* accum = NULL;
    size_t accum_len = 0;
    size_t accum_cap = 0;

    for (;;) {
        const char* prompt = accum ? ">> " : "alea> ";
        char* line = linenoise(prompt);

        if (!line) {
            if (accum) {
                free(accum);
                accum = NULL;
                accum_len = accum_cap = 0;
                printf("\n");
                continue;
            }
            printf("\n");
            break;
        }

        if (line[0] == '\0' && !accum) {
            free(line);
            continue;
        }

        /* Build the full source string */
        const char* source;
        if (!accum) {
            source = line;
        } else {
            append(&accum, &accum_len, &accum_cap, "\n");
            append(&accum, &accum_len, &accum_cap, line);
            source = accum;
        }

        /* Test-compile to check for incomplete input */
        int status = luaL_loadbuffer(L, source, strlen(source), "=stdin");
        if (incomplete(L, status)) {
            lua_pop(L, 1); /* pop error */
            if (!accum) {
                accum_len = 0;
                accum_cap = 0;
                accum = NULL;
                append(&accum, &accum_len, &accum_cap, line);
            }
            linenoiseHistoryAdd(line);
            free(line);
            continue;
        }
        lua_pop(L, 1); /* pop test-compile result (chunk or error) */

        /* Input is complete — execute */
        linenoiseHistoryAdd(accum ? accum : line);
        repl_execute(L, accum ? accum : line, accum == NULL);

        free(line);
        if (accum) {
            free(accum);
            accum = NULL;
            accum_len = accum_cap = 0;
        }
    }

    linenoiseHistorySave(hist);
}

#else
/* ---- Windows REPL: plain fgets, no line editing ---- */

static char* read_line(const char* prompt) {
    static char buf[4096];
    fputs(prompt, stdout);
    fflush(stdout);
    if (!fgets(buf, sizeof(buf), stdin))
        return NULL;
    /* Strip trailing newline */
    size_t len = strlen(buf);
    if (len > 0 && buf[len - 1] == '\n')
        buf[--len] = '\0';
    if (len > 0 && buf[len - 1] == '\r')
        buf[--len] = '\0';
    /* Return heap copy (caller frees) */
    char* copy = (char*)malloc(len + 1);
    if (copy)
        memcpy(copy, buf, len + 1);
    return copy;
}

static void repl(lua_State* L) {
    printf("alea %s (Lua 5.5) -- type Ctrl-Z then Enter to exit\n", alea_version());

    char* accum = NULL;
    size_t accum_len = 0;
    size_t accum_cap = 0;

    for (;;) {
        const char* prompt = accum ? ">> " : "alea> ";
        char* line = read_line(prompt);

        if (!line) {
            if (accum) {
                free(accum);
                accum = NULL;
                accum_len = accum_cap = 0;
                printf("\n");
                continue;
            }
            printf("\n");
            break;
        }

        if (line[0] == '\0' && !accum) {
            free(line);
            continue;
        }

        const char* source;
        if (!accum) {
            source = line;
        } else {
            append(&accum, &accum_len, &accum_cap, "\n");
            append(&accum, &accum_len, &accum_cap, line);
            source = accum;
        }

        int status = luaL_loadbuffer(L, source, strlen(source), "=stdin");
        if (incomplete(L, status)) {
            lua_pop(L, 1);
            if (!accum) {
                accum_len = 0;
                accum_cap = 0;
                accum = NULL;
                append(&accum, &accum_len, &accum_cap, line);
            }
            free(line);
            continue;
        }
        lua_pop(L, 1);

        int was_first = (accum == NULL);
        const char* final_source = accum ? accum : line;

        repl_execute(L, final_source, was_first);

        free(line);
        if (accum) {
            free(accum);
            accum = NULL;
            accum_len = accum_cap = 0;
        }
    }
}

#endif /* _WIN32 */

/* ============================================================================
 * Main
 * ============================================================================ */

int main(int argc, char** argv) {
    const char* exec_chunk = NULL;
    const char* script = NULL;
    int log_level = -1;
    int script_idx = -1;

    const char* params[64];
    int n_params = 0;

    int i;
    for (i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-v") == 0 || strcmp(argv[i], "--version") == 0) {
            printf("alea %s (Lua 5.5)\n", alea_version());
            return 0;
        }
        if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0) {
            print_usage(argv[0]);
            return 0;
        }
        if (strcmp(argv[i], "-e") == 0) {
            if (++i >= argc) {
                fprintf(stderr, "Error: -e requires an argument\n");
                return 1;
            }
            exec_chunk = argv[i];
            continue;
        }
        if (strcmp(argv[i], "--set") == 0) {
            if (++i >= argc) {
                fprintf(stderr, "Error: --set requires key=value\n");
                return 1;
            }
            if (n_params < 64)
                params[n_params++] = argv[i];
            continue;
        }
        if (strcmp(argv[i], "-l") == 0) {
            if (++i >= argc) {
                fprintf(stderr, "Error: -l requires a level\n");
                return 1;
            }
            log_level = atoi(argv[i]);
            continue;
        }
        if (argv[i][0] == '-') {
            fprintf(stderr, "Unknown option: %s\n", argv[i]);
            print_usage(argv[0]);
            return 1;
        }
        script = argv[i];
        script_idx = i;
        break;
    }

    (void)log_level;

    lua_State* L = luaL_newstate();
    if (!L) {
        fprintf(stderr, "Error: cannot create Lua state\n");
        return 1;
    }

    luaL_openlibs(L);
    luaL_requiref(L, "alea", luaopen_alea, 1);
    lua_pop(L, 1);

    for (int j = 0; j < n_params; j++)
        inject_param(L, params[j]);

    g_lua_state = L;
    signal(SIGINT, sigint_handler);

    int status;

    if (exec_chunk) {
        status = luaL_loadstring(L, exec_chunk);
        if (status == LUA_OK)
            status = docall(L, 0, 0);
    } else if (script) {
        create_arg_table(L, argc, argv, script_idx);
        status = luaL_loadfile(L, script);
        if (status == LUA_OK)
            status = docall(L, 0, LUA_MULTRET);
    } else if (IS_TTY()) {
        /* Interactive REPL */
        repl(L);
        status = LUA_OK;
    } else {
        /* Pipe/stdin mode */
        status = luaL_loadfile(L, NULL);
        if (status == LUA_OK)
            status = docall(L, 0, LUA_MULTRET);
    }

    if (status != LUA_OK) {
        const char* msg = lua_tostring(L, -1);
        if (msg)
            fprintf(stderr, "%s\n", msg);
        lua_pop(L, 1);
    }

    lua_close(L);
    alea_clear_interrupt();

    return status != LUA_OK ? 1 : 0;
}
