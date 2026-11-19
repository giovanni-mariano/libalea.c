// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
// SPDX-License-Identifier: MPL-2.0

/* Build-host tool. Hash each relative path, NUL separator and raw file contents
 * in component-sorted order, as required by the snapshot compatibility stamp. */
#ifndef _WIN32
#define _POSIX_C_SOURCE 200809L
#endif
#include "util/alea_digest.h"
#include "util/compat.h"
#include <errno.h>
#include <stdlib.h>
#include <string.h>
#ifdef _WIN32
#include <windows.h>
#else
#include <sys/stat.h>
#endif

static void fail(const char* path) {
    fprintf(stderr, "Cannot generate snapshot compatibility hash at %s\n", path);
    exit(EXIT_FAILURE);
}

static int compare_names(const void* lhs, const void* rhs) {
    const char* a = *(char* const*)lhs;
    const char* b = *(char* const*)rhs;
    return strcmp(a, b);
}

static void hash_file(alea_sha256_t* hash, const char* path) {
    FILE* input = fopen(path, "rb");
    if (!input) fail(path);
    alea_sha256_update(hash, path, strlen(path) + 1);
    unsigned char chunk[16384];
    size_t n;
    while ((n = fread(chunk, 1, sizeof(chunk), input)) != 0)
        alea_sha256_update(hash, chunk, n);
    if (ferror(input) || fclose(input)) fail(path);
}

static void hash_directory(alea_sha256_t* hash, const char* path) {
    alea_dir_t* dir = alea_dir_open(path);
    if (!dir) fail(path);
    char** names = NULL;
    size_t count = 0, capacity = 0;
    const char* name;
    for (;;) {
        errno = 0;
        name = alea_dir_next(dir);
        if (!name) {
#ifdef _WIN32
            if (GetLastError() != ERROR_NO_MORE_FILES) fail(path);
#else
            if (errno) fail(path);
#endif
            break;
        }
        if (!strcmp(name, ".") || !strcmp(name, "..")) continue;
        if (count == capacity) {
            capacity = capacity ? 2 * capacity : 32;
            names = realloc(names, capacity * sizeof(*names));
            if (!names) fail(path);
        }
        names[count] = malloc(strlen(path) + strlen(name) + 2);
        if (!names[count]) fail(path);
        sprintf(names[count++], "%s/%s", path, name);
    }
    alea_dir_close(dir);
    if (count) qsort(names, count, sizeof(*names), compare_names);
    for (size_t i = 0; i < count; ++i) {
        char* child = names[i];
#ifdef _WIN32
        DWORD attr = GetFileAttributesA(child);
        if (attr == INVALID_FILE_ATTRIBUTES) fail(child);
        bool directory = (attr & FILE_ATTRIBUTE_DIRECTORY) != 0;
        bool link = (attr & FILE_ATTRIBUTE_REPARSE_POINT) != 0;
#else
        struct stat st;
        if (lstat(child, &st)) fail(child);
        bool link = S_ISLNK(st.st_mode);
        if (link && stat(child, &st)) fail(child);
        bool directory = S_ISDIR(st.st_mode);
#endif
        if (directory) {
            if (!link) hash_directory(hash, child);
        } else {
            const char* ext = strrchr(child, '.');
            if (ext && (!strcmp(ext, ".c") || !strcmp(ext, ".h"))) hash_file(hash, child);
        }
        free(child);
    }
    free(names);
}

static void write_changed(const char* path, const char* line) {
    char old[128];
    FILE* file = fopen(path, "rb");
    if (file) {
        size_t n = fread(old, 1, sizeof(old), file);
        if (ferror(file) || fclose(file)) fail(path);
        if (n == strlen(line) && !memcmp(old, line, n)) return;
    }
    size_t capacity = strlen(path) + 256;
    char* temporary = malloc(capacity);
    if (!temporary) fail(path);
    file = alea_sibling_tmpfile(path, temporary, capacity);
    if (!file) fail(path);
    bool written = fputs(line, file) >= 0;
    if (fclose(file) || !written || alea_replace_file(temporary, path)) {
        remove(temporary);
        fail(path);
    }
    free(temporary);
}

int main(int argc, char** argv) {
    if (argc != 2) return EXIT_FAILURE;
    alea_sha256_t hash;
    unsigned char digest[32];
    alea_sha256_init(&hash);
    hash_directory(&hash, "include");
    hash_directory(&hash, "src");
    alea_sha256_final(&hash, digest);
    char hex[65], header[128];
    for (size_t i = 0; i < sizeof(digest); ++i) sprintf(hex + 2*i, "%02x", (unsigned)digest[i]);
    snprintf(header, sizeof(header), "#define ALEA_ENGINE_HASH \"%s\"\n", hex);
    write_changed(argv[1], header);
    return EXIT_SUCCESS;
}
