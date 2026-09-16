// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

#include "cluster_internal.h"

#include <ctype.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#define MCNP_INCLUDE_DEPTH_LIMIT 32u

typedef struct {
    char* data;
    size_t length;
    size_t capacity;
} input_builder_t;

static alea_cluster_status_t input_append(input_builder_t* out,
                                          const char* bytes, size_t count) {
    if (count > SIZE_MAX - out->length - 1)
        return ALEA_CLUSTER_OUT_OF_MEMORY;
    size_t needed = out->length + count + 1;
    if (needed > out->capacity) {
        size_t capacity = out->capacity ? out->capacity : 4096;
        while (capacity < needed) {
            if (capacity > SIZE_MAX / 2) { capacity = needed; break; }
            capacity *= 2;
        }
        char* grown = realloc(out->data, capacity);
        if (!grown) return ALEA_CLUSTER_OUT_OF_MEMORY;
        out->data = grown;
        out->capacity = capacity;
    }
    if (count) memcpy(out->data + out->length, bytes, count);
    out->length += count;
    out->data[out->length] = '\0';
    return ALEA_CLUSTER_OK;
}

static int keyword(const char* text, const char* end, const char* word) {
    while (*word) {
        if (text == end || tolower((unsigned char)*text++) != *word++)
            return 0;
    }
    return text == end || isspace((unsigned char)*text) || *text == '=';
}

/* Returns 1 for READ FILE=, 0 for an ordinary line, -1 for a malformed
 * READ card, and -2 for allocation failure. The caller owns *name. */
static int mcnp_include_name(const char* begin, const char* end,
                             char** name) {
    *name = NULL;
    while (begin < end && isspace((unsigned char)*begin)) ++begin;
    if (begin == end || (tolower((unsigned char)*begin) == 'c' &&
                         begin + 1 < end && isspace((unsigned char)begin[1])))
        return 0;
    if (!keyword(begin, end, "read")) return 0;
    begin += 4;
    if (begin == end || !isspace((unsigned char)*begin)) return 0;
    while (begin < end && isspace((unsigned char)*begin)) ++begin;
    if (!keyword(begin, end, "file")) return -1;
    begin += 4;
    while (begin < end && isspace((unsigned char)*begin)) ++begin;
    if (begin == end || *begin++ != '=') return -1;
    while (begin < end && isspace((unsigned char)*begin)) ++begin;
    if (begin == end) return -1;
    const char* start = begin;
    size_t count;
    if (*begin == '\'' || *begin == '"') {
        char quote = *begin++;
        start = begin;
        while (begin < end && *begin != quote) ++begin;
        if (begin == end) return -1;
        count = (size_t)(begin - start);
        ++begin;
    } else {
        while (begin < end && !isspace((unsigned char)*begin) &&
               *begin != '$') ++begin;
        count = (size_t)(begin - start);
    }
    while (begin < end && isspace((unsigned char)*begin)) ++begin;
    if (count == 0 || (begin < end && *begin != '$')) return -1;
    char* copy = malloc(count + 1);
    if (!copy) return -2;
    memcpy(copy, start, count);
    copy[count] = '\0';
    *name = copy;
    return 1;
}

static char* resolve_path(const char* parent, const char* child) {
    if (*child == '/') {
        size_t length = strlen(child);
        char* copy = malloc(length + 1);
        if (copy) memcpy(copy, child, length + 1);
        return copy;
    }
    const char* slash = strrchr(parent, '/');
    size_t prefix = slash ? (size_t)(slash - parent + 1) : 0;
    size_t suffix = strlen(child);
    if (prefix > SIZE_MAX - suffix - 1) return NULL;
    char* path = malloc(prefix + suffix + 1);
    if (!path) return NULL;
    if (prefix) memcpy(path, parent, prefix);
    memcpy(path + prefix, child, suffix + 1);
    return path;
}

static alea_cluster_status_t expand_file(const char* path,
                                          input_builder_t* out,
                                          size_t depth) {
    if (depth >= MCNP_INCLUDE_DEPTH_LIMIT)
        return ALEA_CLUSTER_COMPUTE_ERROR;
    char* input = NULL;
    size_t length = 0;
    alea_cluster_status_t status = alea_cluster_read_path_root(
        path, &input, &length);
    if (status != ALEA_CLUSTER_OK) return status;
    if (memchr(input, '\0', length)) {
        free(input);
        return ALEA_CLUSTER_COMPUTE_ERROR;
    }
    const char* cursor = input;
    const char* end = input + length;
    int title_seen = depth != 0;
    while (cursor < end && status == ALEA_CLUSTER_OK) {
        int append_newline = 1;
        const char* line_end = memchr(cursor, '\n', (size_t)(end - cursor));
        line_end = line_end ? line_end : end;
        const char* logical_end = line_end;
        if (logical_end > cursor && logical_end[-1] == '\r') --logical_end;
        const char* trimmed = cursor;
        while (trimmed < logical_end && isspace((unsigned char)*trimmed))
            ++trimmed;
        if (!title_seen && trimmed < logical_end) {
            title_seen = 1;
            status = input_append(out, cursor,
                                  (size_t)(line_end - cursor));
        } else {
            char* name = NULL;
            int included = mcnp_include_name(cursor, logical_end, &name);
            if (included == -2) status = ALEA_CLUSTER_OUT_OF_MEMORY;
            else if (included < 0) status = ALEA_CLUSTER_COMPUTE_ERROR;
            else if (included == 1) {
                char* child_path = resolve_path(path, name);
                free(name);
                if (!child_path) status = ALEA_CLUSTER_OUT_OF_MEMORY;
                else {
                    status = expand_file(child_path, out, depth + 1);
                    free(child_path);
                    if (status == ALEA_CLUSTER_OK && out->length &&
                        out->data[out->length - 1] == '\n')
                        append_newline = 0;
                }
            } else {
                status = input_append(out, cursor,
                                      (size_t)(line_end - cursor));
            }
        }
        if (status == ALEA_CLUSTER_OK && append_newline && line_end < end)
            status = input_append(out, "\n", 1);
        cursor = line_end < end ? line_end + 1 : end;
    }
    free(input);
    return status;
}

alea_cluster_status_t alea_cluster_read_mcnp_input(
        alea_cluster_t* cluster, const char* path,
        char** data, size_t* length) {
    if (!cluster) return ALEA_CLUSTER_INVALID_ARGUMENT;
    input_builder_t expanded = {0};
    alea_cluster_status_t local = ALEA_CLUSTER_OK;
    if (cluster->rank == 0) {
        local = path ? expand_file(path, &expanded, 0)
                     : ALEA_CLUSTER_INVALID_ARGUMENT;
        if (local == ALEA_CLUSTER_OK && !expanded.data)
            local = input_append(&expanded, "", 0);
    }
    return alea_cluster_broadcast_owned_bytes(cluster,
        expanded.data, expanded.length, local, data, length);
}
