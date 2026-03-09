// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/**
 * @file xsdir.c
 * @brief xsdir/xsdata directory file parser
 *
 * Parses MCNP-style xsdir files that map ZAIDs to ACE data files.
 *
 * Format:
 *   - Optional "datapath=" line specifying base directory
 *   - "directory" keyword marks start of entries
 *   - Each entry: ZAID AWR filename access file_type address length recl entries temp
 *   - Lines ending with '+' are continuation lines
 */

#include "alea_nucdata.h"
#include "util/alea_log.h"
#include "util/compat.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#define XSDIR_INITIAL_CAPACITY 1024
#define XSDIR_MAX_LINE 4096

/** Trim leading and trailing whitespace in-place, return pointer to start */
static char* trim(char* s) {
    while (isspace((unsigned char)*s)) s++;
    if (*s == '\0') return s;
    char* end = s + strlen(s) - 1;
    while (end > s && isspace((unsigned char)*end)) *end-- = '\0';
    return s;
}

/** Determine table type from ZAID suffix character */
static alea_nuc_table_type_t table_type_from_zaid(const char* zaid) {
    const char* dot = strchr(zaid, '.');
    if (!dot) return ALEA_NUC_TABLE_CONTINUOUS_NEUTRON;
    const char* p = dot + 1;
    while (*p >= '0' && *p <= '9') p++;
    switch (*p) {
    case 'c': return ALEA_NUC_TABLE_CONTINUOUS_NEUTRON;
    case 'p': return ALEA_NUC_TABLE_PHOTOATOMIC;
    case 'u': return ALEA_NUC_TABLE_PHOTONUCLEAR;
    case 't': return ALEA_NUC_TABLE_THERMAL_SAB;
    case 'e': return ALEA_NUC_TABLE_ELECTRON;
    }
    return ALEA_NUC_TABLE_CONTINUOUS_NEUTRON;
}

static alea_error_t xsdir_add_entry(alea_nuc_xsdir_t* xsdir, const char* line) {
    /* Grow if needed */
    if (xsdir->count >= xsdir->capacity) {
        size_t new_cap = xsdir->capacity ? xsdir->capacity * 2 : XSDIR_INITIAL_CAPACITY;
        alea_nuc_xsdir_entry_t* p = realloc(xsdir->entries, new_cap * sizeof(*p));
        if (!p) return ALEA_ERR_OUT_OF_MEMORY;
        xsdir->entries = p;
        xsdir->capacity = new_cap;
    }

    alea_nuc_xsdir_entry_t* e = &xsdir->entries[xsdir->count];
    memset(e, 0, sizeof(*e));

    /*
     * Fields: ZAID  AWR  filename  access  file_type  address  length  recl  entries  temp
     * Some xsdir files have fewer fields — temperature may be missing.
     */
    int n = sscanf(line, "%23s %lf %511s %d %d %d %d %d %d %lf",
                   e->zaid, &e->awr, e->filename,
                   &e->access_route, &e->file_type, &e->address,
                   &e->table_length, &e->record_length, &e->num_entries,
                   &e->temperature);

    if (n < 7) return ALEA_ERR_PARSE_ERROR;

    e->type = table_type_from_zaid(e->zaid);
    xsdir->count++;
    return ALEA_OK;
}

/** Internal: parse xsdir file contents into an existing xsdir struct */
static alea_error_t xsdir_parse_file(alea_nuc_xsdir_t* xsdir, const char* path) {
    FILE* fp = fopen(path, "r");
    if (!fp) return ALEA_ERR_FILE_NOT_FOUND;

    /* Extract directory from path for relative filenames */
    const char* last_sep = strrchr(path, '/');
    if (last_sep) {
        size_t len = (size_t)(last_sep - path);
        if (len >= sizeof(xsdir->datapath)) len = sizeof(xsdir->datapath) - 1;
        memcpy(xsdir->datapath, path, len);
        xsdir->datapath[len] = '\0';
    }

    char line[XSDIR_MAX_LINE];
    char accumulated[XSDIR_MAX_LINE * 4]; /* for continuation lines */
    bool in_directory = false;
    bool has_continuation = false;
    accumulated[0] = '\0';

    while (fgets(line, sizeof(line), fp)) {
        char* s = trim(line);

        /* Skip empty lines and comments */
        if (*s == '\0' || *s == '#') continue;

        /* Look for datapath */
        if (alea_strncasecmp(s, "datapath", 8) == 0) {
            char* eq = strchr(s, '=');
            if (eq) {
                char* dp = trim(eq + 1);
                strncpy(xsdir->datapath, dp, sizeof(xsdir->datapath) - 1);
            }
            continue;
        }

        /* Look for "directory" keyword */
        if (alea_strncasecmp(s, "directory", 9) == 0) {
            in_directory = true;
            continue;
        }

        /*
         * If no "directory" keyword seen yet, try to parse lines anyway.
         * FENDL .xsd files and some xsdir variants have bare entries
         * without a "directory" header. Skip known non-entry keywords.
         */
        if (!in_directory) {
            /* Skip known header keywords */
            if (alea_strncasecmp(s, "atomic", 6) == 0) continue;
            /* Try parsing — if it works, treat as bare entry */
        }

        /* Handle continuation lines (ending with '+') */
        size_t slen = strlen(s);
        if (slen > 0 && s[slen - 1] == '+') {
            s[slen - 1] = '\0'; /* remove '+' */
            if (has_continuation) {
                strncat(accumulated, " ", sizeof(accumulated) - strlen(accumulated) - 1);
                strncat(accumulated, s, sizeof(accumulated) - strlen(accumulated) - 1);
            } else {
                strncpy(accumulated, s, sizeof(accumulated) - 1);
                accumulated[sizeof(accumulated) - 1] = '\0';
            }
            has_continuation = true;
            continue;
        }

        if (has_continuation) {
            strncat(accumulated, " ", sizeof(accumulated) - strlen(accumulated) - 1);
            strncat(accumulated, s, sizeof(accumulated) - strlen(accumulated) - 1);
            s = accumulated;
            has_continuation = false;
        }

        xsdir_add_entry(xsdir, s);

        if (!has_continuation)
            accumulated[0] = '\0';
    }

    fclose(fp);
    return ALEA_OK;
}

alea_nuc_xsdir_t* alea_nuc_xsdir_load(const char* path) {
    if (!path) return NULL;

    alea_nuc_xsdir_t* xsdir = calloc(1, sizeof(*xsdir));
    if (!xsdir) return NULL;

    alea_error_t err = xsdir_parse_file(xsdir, path);
    if (err != ALEA_OK || xsdir->count == 0) {
        free(xsdir->entries);
        free(xsdir);
        return NULL;
    }

    ALEA_LOG_INFO("loaded xsdir with %zu entries from %s", xsdir->count, path);
    return xsdir;
}

alea_nuc_xsdir_t* alea_nuc_xsdir_load_dir(const char* dirpath) {
    if (!dirpath) return NULL;

    alea_dir_t* dir = alea_dir_open(dirpath);
    if (!dir) return NULL;

    alea_nuc_xsdir_t* xsdir = calloc(1, sizeof(*xsdir));
    if (!xsdir) { alea_dir_close(dir); return NULL; }

    /* Set datapath to the directory so relative filenames resolve */
    strncpy(xsdir->datapath, dirpath, sizeof(xsdir->datapath) - 1);
    xsdir->datapath[sizeof(xsdir->datapath) - 1] = '\0';

    size_t dirpath_len = strlen(dirpath);
    char line[XSDIR_MAX_LINE];
    size_t loaded = 0;

    const char* name;
    while ((name = alea_dir_next(dir)) != NULL) {
        /* Only process .xsd files */
        const char* ext = strrchr(name, '.');
        if (!ext || strcmp(ext, ".xsd") != 0) continue;

        size_t name_len = strlen(name);
        size_t path_len = dirpath_len + 1 + name_len + 1;
        char filepath_buf[1024];
        char* filepath = filepath_buf;
        char* filepath_alloc = NULL;
        if (path_len > sizeof(filepath_buf)) {
            filepath_alloc = malloc(path_len);
            if (!filepath_alloc) continue;
            filepath = filepath_alloc;
        }
        snprintf(filepath, path_len, "%s/%s", dirpath, name);

        FILE* fp = fopen(filepath, "r");
        if (!fp) { free(filepath_alloc); continue; }

        while (fgets(line, sizeof(line), fp)) {
            char* s = trim(line);
            if (*s == '\0' || *s == '#') continue;

            if (xsdir_add_entry(xsdir, s) == ALEA_OK)
                loaded++;
        }
        fclose(fp);
        free(filepath_alloc);
    }
    alea_dir_close(dir);

    if (loaded == 0) {
        free(xsdir->entries);
        free(xsdir);
        return NULL;
    }

    ALEA_LOG_INFO("loaded %zu entries from .xsd files in %s", loaded, dirpath);
    return xsdir;
}

void alea_nuc_xsdir_free(alea_nuc_xsdir_t* xsdir) {
    if (!xsdir) return;
    /* Free cached nuclides */
    if (xsdir->cache) {
        for (size_t i = 0; i < xsdir->cache_capacity; i++) {
            if (xsdir->cache[i].zaid[0])
                alea_nuc_nuclide_free(xsdir->cache[i].nuclide);
        }
        free(xsdir->cache);
    }
    free(xsdir->entries);
    free(xsdir);
}

const alea_nuc_xsdir_entry_t* alea_nuc_xsdir_find(const alea_nuc_xsdir_t* xsdir,
                                                   const char* zaid) {
    if (!xsdir || !zaid) return NULL;
    for (size_t i = 0; i < xsdir->count; i++) {
        if (strcmp(xsdir->entries[i].zaid, zaid) == 0)
            return &xsdir->entries[i];
    }
    return NULL;
}

size_t alea_nuc_xsdir_count(const alea_nuc_xsdir_t* xsdir) {
    return xsdir ? xsdir->count : 0;
}

/* ============================================================================
 * Nuclide cache — open-addressing hash table
 * ============================================================================ */

#define CACHE_INITIAL_CAPACITY 64   /* must be power of 2 */
#define CACHE_MAX_LOAD_FACTOR  0.7

static uint32_t cache_hash(const char* s) {
    /* FNV-1a */
    uint32_t h = 2166136261u;
    for (; *s; s++) {
        h ^= (uint8_t)*s;
        h *= 16777619u;
    }
    return h;
}

static int cache_grow(alea_nuc_xsdir_t* xsdir) {
    size_t new_cap = xsdir->cache_capacity * 2;
    alea_nuc_cache_entry_t* new_cache = calloc(new_cap, sizeof(*new_cache));
    if (!new_cache) return -1;

    /* Re-insert all existing entries */
    uint32_t mask = (uint32_t)(new_cap - 1);
    for (size_t i = 0; i < xsdir->cache_capacity; i++) {
        if (!xsdir->cache[i].zaid[0]) continue;
        uint32_t idx = cache_hash(xsdir->cache[i].zaid) & mask;
        while (new_cache[idx].zaid[0])
            idx = (idx + 1) & mask;
        new_cache[idx] = xsdir->cache[i];
    }

    free(xsdir->cache);
    xsdir->cache = new_cache;
    xsdir->cache_capacity = new_cap;
    return 0;
}

alea_nuc_nuclide_t* alea_nuc_xsdir_get_nuclide(alea_nuc_xsdir_t* xsdir, const char* zaid) {
    if (!xsdir || !zaid) return NULL;

    /* Initialize cache on first use */
    if (!xsdir->cache) {
        xsdir->cache = calloc(CACHE_INITIAL_CAPACITY, sizeof(*xsdir->cache));
        if (!xsdir->cache) return NULL;
        xsdir->cache_capacity = CACHE_INITIAL_CAPACITY;
        xsdir->cache_count = 0;
    }

    /* Lookup */
    uint32_t mask = (uint32_t)(xsdir->cache_capacity - 1);
    uint32_t idx = cache_hash(zaid) & mask;
    while (xsdir->cache[idx].zaid[0]) {
        if (strcmp(xsdir->cache[idx].zaid, zaid) == 0)
            return xsdir->cache[idx].nuclide;
        idx = (idx + 1) & mask;
    }

    /* Cache miss — load the nuclide */
    alea_nuc_nuclide_t* nuc = alea_nuc_load_nuclide(xsdir, zaid);
    if (!nuc) return NULL;

    /* Grow if needed before inserting */
    if ((double)(xsdir->cache_count + 1) >
        (double)xsdir->cache_capacity * CACHE_MAX_LOAD_FACTOR) {
        if (cache_grow(xsdir) != 0) {
            alea_nuc_nuclide_free(nuc);
            return NULL;
        }
        /* Recompute insertion slot after grow */
        mask = (uint32_t)(xsdir->cache_capacity - 1);
        idx = cache_hash(zaid) & mask;
        while (xsdir->cache[idx].zaid[0])
            idx = (idx + 1) & mask;
    }

    /* Insert */
    strncpy(xsdir->cache[idx].zaid, zaid, sizeof(xsdir->cache[idx].zaid) - 1);
    xsdir->cache[idx].zaid[sizeof(xsdir->cache[idx].zaid) - 1] = '\0';
    xsdir->cache[idx].nuclide = nuc;
    xsdir->cache_count++;

    return nuc;
}
