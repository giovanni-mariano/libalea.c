// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
// SPDX-License-Identifier: MPL-2.0
#ifndef ALEA_SNAPSHOT_INTERNAL_H
#define ALEA_SNAPSHOT_INTERNAL_H
#include "alea_system.h"
#include "alea_snapshot.h"
#include "util/alea_digest.h"
#include "util/compat.h"
#include <stdio.h>
typedef struct { FILE* file; int read,failed; uint64_t remaining; alea_sha256_t digest; size_t position, buffered; unsigned char buffer[65536]; alea_mapped_file_t mapping; int mapped,profile; double decode_hash_seconds; } alea_snapshot_io;
void alea_snapshot_bytes(alea_snapshot_io* io,void* data,size_t bytes);
uint64_t alea_snapshot_uint(alea_snapshot_io* io,uint64_t value);
float alea_snapshot_float(alea_snapshot_io* io,float value);
double alea_snapshot_double(alea_snapshot_io* io,double value);
void alea_snapshot_string(alea_snapshot_io* io,char** text);
size_t alea_snapshot_count(alea_snapshot_io* io,size_t count,size_t element_size);
const char* alea_snapshot_engine_id(void);
int alea_snapshot_begin(alea_snapshot_io* io,FILE* file,int read,size_t max_bytes,unsigned kind);
/* Opt-in ALEA_SNAPSHOT_IO=mapped; otherwise buffered, with mapping-failure fallback.
 * Mapping is released before returning a restored, fully owned editable model. */
int alea_snapshot_open_read(alea_snapshot_io* io,const char* path,size_t max_bytes,unsigned kind);
void alea_snapshot_close_read(alea_snapshot_io* io);
int alea_snapshot_end(alea_snapshot_io* io);
int alea_snapshot_hierarchy(alea_snapshot_io* io,alea_system_t* sys);
int alea_snapshot_system(alea_snapshot_io* io,alea_system_t* sys);
int alea_snapshot_read_limits(const alea_snapshot_read_options_t* options,
                              size_t* file_bytes,size_t* decoded_bytes);
void alea_snapshot_profile_allocations(const alea_alloc_budget_t* budget);
/* Convert to the destination width explicitly, then check that the field
 * represents the complete decoded value (including signed or bit fields). */
#define SU(v) do { \
    uint64_t sv = alea_snapshot_uint(io, (uint64_t)(v)); \
    if (io->read) { \
        (v) = _Generic((v), bool: (bool)sv, int8_t: (int8_t)sv, uint8_t: (uint8_t)sv, \
                       int16_t: (int16_t)sv, uint16_t: (uint16_t)sv, \
                       int32_t: (int32_t)sv, uint32_t: (uint32_t)sv, \
                       int64_t: (int64_t)sv, uint64_t: sv, default: (unsigned int)sv); \
        if ((uint64_t)(v) != sv) io->failed = 1; \
    } \
} while (0)
#define SD(v) do { double sv=alea_snapshot_double(io,(v));if(io->read){(v)=sv;if((double)(v)!=sv)io->failed=1;} } while(0)
#define SF(v) do { float sv=alea_snapshot_float(io,(v));if(io->read)(v)=sv; } while(0)
#define SS(v) alea_snapshot_string(io,&(v))
#define SC(count,ptr) do { size_t sn=alea_snapshot_count(io,(count),sizeof(*(ptr))); if(io->read){(ptr)=alea_calloc(sn?sn:1,sizeof(*(ptr)));if(!(ptr)){io->failed=1;sn=0;}(count)=sn;} } while(0)
#endif
