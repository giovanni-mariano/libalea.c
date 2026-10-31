// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
// SPDX-License-Identifier: MPL-2.0
#ifndef ALEA_COMPUTE_ROOT_REUSE_H
#define ALEA_COMPUTE_ROOT_REUSE_H
#include "compute_layout.h"

/* Diagnostic only: simulate direct-mapped caches, never substitute roots.
 * Snapshot ownership is implicit. Reset keys for every ray/AA sample and
 * counters for every frame. Exact local-ray bits distinguish placements. */
#define ROOT_REUSE_FAMILIES 14
#define ROOT_REUSE_CAPACITIES 4
#define ROOT_REUSE_SLOTS (8+32+128+512)
#define ROOT_REUSE_MAX_PIXELS 1024
#define ROOT_REUSE_COUNTERS (3+2*ROOT_REUSE_CAPACITIES)
enum { ROOT_REQUESTS, ROOT_ERRORS, ROOT_NONFINITE, ROOT_HITS, ROOT_EVICTIONS=ROOT_HITS+ROOT_REUSE_CAPACITIES };
typedef struct { compute_uint valid, key[7]; } compute_root_key;
typedef struct {
    compute_root_key slots[ROOT_REUSE_SLOTS];
    compute_uint counts[ROOT_REUSE_FAMILIES][ROOT_REUSE_COUNTERS];
    compute_uint saturated;
} compute_root_reuse;
#ifdef __METAL_VERSION__
#define REUSE_DEVICE device
#define REUSE_THREAD thread
#else
#define REUSE_DEVICE
#define REUSE_THREAD
#endif
static inline compute_uint root_reuse_capacity(compute_uint c) { return 8u<<(2*c); }
static inline void root_reuse_reset(REUSE_DEVICE compute_root_reuse* s) {
    for(compute_uint i=0;i<ROOT_REUSE_SLOTS;++i)s->slots[i].valid=0;
}
static inline void root_reuse_add(REUSE_DEVICE compute_root_reuse* s,compute_uint family,compute_uint counter) {
    if(s->counts[family][counter]==COMPUTE_NONE)s->saturated=1;
    else ++s->counts[family][counter];
}
static inline void root_reuse_record(REUSE_DEVICE compute_root_reuse* s,
        REUSE_THREAD const compute_uint* key,compute_uint family,compute_uint rejected) {
    root_reuse_add(s,family,ROOT_REQUESTS);
    if(rejected){root_reuse_add(s,family,rejected==1?ROOT_ERRORS:ROOT_NONFINITE);return;}
    compute_uint hash=2166136261u;
    for(compute_uint j=0;j<7;++j)hash=(hash^key[j])*16777619u;
    /* Mix high float bits into the small table index. */
    hash^=hash>>16;hash*=0x7feb352du;hash^=hash>>15;
    compute_uint offset=0;
    for(compute_uint c=0;c<ROOT_REUSE_CAPACITIES;++c) {
        compute_uint capacity=root_reuse_capacity(c);
        REUSE_DEVICE compute_root_key* slot=&s->slots[offset+(hash&(capacity-1))];
        compute_uint same=slot->valid;
        for(compute_uint j=0;j<7;++j)same&=slot->key[j]==key[j];
        if(same)root_reuse_add(s,family,ROOT_HITS+c);
        else {
            if(slot->valid)root_reuse_add(s,family,ROOT_EVICTIONS+c);
            for(compute_uint j=0;j<7;++j)slot->key[j]=key[j];
            slot->valid=1;
        }
        offset+=capacity;
    }
}
#undef REUSE_DEVICE
#undef REUSE_THREAD
#endif
