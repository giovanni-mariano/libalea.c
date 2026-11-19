// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
// SPDX-License-Identifier: MPL-2.0
#ifndef ALEA_SNAPSHOT_H
#define ALEA_SNAPSHOT_H
#include "alea.h"
#ifdef __cplusplus
extern "C" {
#endif
/** Separate file and decoded-allocation limits. Zero selects the defaults:
 * 2 GiB file size, 8 GiB cumulative allocation requests (capped at SIZE_MAX).
 * Set struct_size to sizeof(alea_snapshot_read_options_t).
 * Decoded bytes include model/wrapper storage, rebuilt indexes, scratch arrays
 * and full vector reallocation capacities. Freed bytes are not refunded.
 * This is not a process RSS limit: stack, allocator overhead, stdio buffers,
 * mapped file pages and acceleration built after the read are excluded. */
typedef struct {
    size_t struct_size;
    size_t max_file_bytes;
    size_t max_decoded_bytes;
} alea_snapshot_read_options_t;
#define ALEA_SNAPSHOT_READ_OPTIONS_INIT {sizeof(alea_snapshot_read_options_t),0,0}

/** Explicit binary geometry snapshots. Restore preserves IDs and model order.
 * Model-wrapper metadata requires the format-specific snapshot API. The caller
 * must exclude concurrent mutation during writing. max_bytes=0 uses 2 GiB.
 * Reads use buffered I/O by default. ALEA_SNAPSHOT_IO=mapped opts into mapped
 * decoding, with buffered fallback if mapping is unavailable. The input file
 * must not be overwritten or truncated while mapped. Restore still allocates
 * an owned editable model and releases the mapping before returning. */
int alea_snapshot_write(const alea_system_t* sys,const char* path);
alea_system_t* alea_snapshot_read(const char* path,size_t max_bytes);
/** NULL options selects defaults. Returns NULL on invalid input, allocation
 * failure or budget exhaustion; no partial model is returned. */
alea_system_t* alea_snapshot_read_with_options(
    const char* path,const alea_snapshot_read_options_t* options);
#ifdef __cplusplus
}
#endif
#endif
