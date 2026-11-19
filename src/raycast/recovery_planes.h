// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
// SPDX-License-Identifier: MPL-2.0
#ifndef ALEA_RECOVERY_PLANES_H
#define ALEA_RECOVERY_PLANES_H

#include "core/alea_system.h"
#include "alea.h"
#include "ray_intersect.h"
#include "recovery_plane_blocks.h"
#include "util/ray_profile.h"
#include "util/alea_alloc.h"
#include <float.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* Optional index for nearest recovery, owned with the universe surface list.
 * Parallel planes share an exact offset search. Root construction can also
 * reject blocks of other planes. The mask preserves original encounter order,
 * including aliases, so equal-distance hits keep their original provenance. */
#define RECOVERY_PLANE_MAX_CARDS 65536u
#define RECOVERY_PLANE_MASK_WORDS (RECOVERY_PLANE_MAX_CARDS / 64)
typedef struct { double a, b, c; uint32_t first, count; } recovery_plane_group_t;
typedef struct { double d; uint32_t ordinal, pad; } recovery_plane_member_t;
typedef struct {
    ray_quadric_packet_t packet;
    uint32_t ordinals[RAY_QUADRIC_LANES];
} recovery_quadric_group_t;

typedef struct recovery_plane_index {
    size_t bytes, quadric_offset;
    uint32_t cards, groups, members, words;
    uint32_t torus_cards, blocks, block_members, quadric_count;
    /* One allocation: groups, members, base mask, block ordinals, blocks, then optional quadrics. */
    recovery_plane_group_t data[];
} recovery_plane_index_t;

static inline const recovery_plane_member_t* recovery_plane_members(const recovery_plane_index_t* p) {
    return (const recovery_plane_member_t*)(p->data + p->groups);
}
static inline const uint64_t* recovery_plane_base(const recovery_plane_index_t* p) {
    return (const uint64_t*)(recovery_plane_members(p) + p->members);
}
/* The ordinal array is padded to an even count to align the blocks. */
static inline const uint32_t* recovery_plane_block_ordinals(const recovery_plane_index_t* p) {
    return (const uint32_t*)(recovery_plane_base(p) + p->words);
}
static inline const recovery_plane_block_t* recovery_plane_blocks(const recovery_plane_index_t* p) {
    return (const recovery_plane_block_t*)(recovery_plane_block_ordinals(p) + ((p->block_members + 1u) & ~1u));
}

static inline int recovery_plane_compare(const void* a, const void* b) {
    const recovery_plane_sort_t *x=a, *y=b;
    int cmp = memcmp(&x->plane, &y->plane, 3 * sizeof(double));
    if (cmp) return cmp;
    if (x->plane.d != y->plane.d) return x->plane.d < y->plane.d ? -1 : 1;
    return (x->ordinal > y->ordinal) - (x->ordinal < y->ordinal);
}

static inline recovery_plane_index_t* recovery_plane_build(
    const alea_system_t* sys, const uint32_t* cards, size_t count, size_t budget) {
    if (count < 128 || count > RECOVERY_PLANE_MAX_CARDS) return NULL;
    recovery_plane_sort_t* sorted = alea_malloc(count * sizeof(*sorted));
    if (!sorted) return NULL;
    size_t n=0;
    uint32_t torus_cards=0;
    for (size_t i=0;i<count;++i) {
        if (cards[i] >= sys->surfaces.count) { free(sorted); return NULL; }
        uint32_t id = sys->surfaces.data[cards[i]].primitive_id;
        if (id >= sys->primitives.count) { free(sorted); return NULL; }
        alea_primitive_type_t type = sys->primitives.data[id].type;
        if (type >= ALEA_PRIMITIVE_TORUS_X && type <= ALEA_PRIMITIVE_TORUS_Z) ++torus_cards;
        if (type != ALEA_PRIMITIVE_PLANE) continue;
        const alea_plane_data_t* p = alea_primitive_payload_const(sys, id);
        if (!p || !isfinite(p->a) || !isfinite(p->b) || !isfinite(p->c) || !isfinite(p->d))
            continue;
        sorted[n++] = (recovery_plane_sort_t){.plane = *p, .ordinal = (uint32_t)i};
    }
    qsort(sorted, n, sizeof(*sorted), recovery_plane_compare);
    size_t groups=0, members=0;
    for (size_t i=0;i<n;) {
        size_t j=i+1;
        while (j<n && !memcmp(&sorted[i].plane, &sorted[j].plane, 3*sizeof(double))) ++j;
        if (j-i >= 8) { ++groups; members += j-i; }
        i=j;
    }
    size_t words = (count+63)/64;
    size_t bytes = sizeof(recovery_plane_index_t) + groups*sizeof(recovery_plane_group_t) +
        members*sizeof(recovery_plane_member_t) + words*sizeof(uint64_t);
    if (!groups || bytes > budget || alea_interrupted()) { free(sorted); return NULL; }
    /* A non-root leaf of the median partition has at least 12 members.
     * Keep the parallel-plane index when extra storage exceeds the budget. */
    size_t other = n - members;
    size_t ordinal_bytes = ((other + 1) & ~(size_t)1) * sizeof(uint32_t);
    size_t block_capacity = other / (RECOVERY_PLANE_BLOCK_SIZE / 2);
    size_t extra = ordinal_bytes + block_capacity * sizeof(recovery_plane_block_t);
    bool add_blocks = other >= RECOVERY_PLANE_BLOCK_SIZE && extra <= budget - bytes;
    size_t allocation = bytes + (add_blocks ? extra : 0);
    recovery_plane_index_t* out = alea_malloc(allocation);
    if (!out && add_blocks) {
        add_blocks = false;
        out = alea_malloc(bytes);
    }
    if (!out) { free(sorted); return NULL; }
    *out = (recovery_plane_index_t){.bytes=bytes, .cards=(uint32_t)count, .groups=(uint32_t)groups,
        .members=(uint32_t)members, .words=(uint32_t)words, .torus_cards=torus_cards};
    recovery_plane_member_t* values = (recovery_plane_member_t*)(out->data + groups);
    uint64_t* base = (uint64_t*)(values + members);
    memset(base, 0xff, words*sizeof(*base));
    size_t g=0, k=0;
    for (size_t i=0;i<n;) {
        size_t j=i+1;
        while (j<n && !memcmp(&sorted[i].plane, &sorted[j].plane, 3*sizeof(double))) ++j;
        if (j-i >= 8) {
            out->data[g++] = (recovery_plane_group_t){sorted[i].plane.a,
                sorted[i].plane.b, sorted[i].plane.c, (uint32_t)k, (uint32_t)(j-i)};
            for (size_t m=i;m<j;++m) {
                uint32_t ordinal = sorted[m].ordinal;
                values[k++] = (recovery_plane_member_t){sorted[m].plane.d, ordinal, 0};
                base[ordinal/64] &= ~(UINT64_C(1) << (ordinal%64));
            }
        }
        i=j;
    }
    if (add_blocks) {
        size_t kept = 0;
        double offset_scale = 1;
        for (size_t i = 0; i < n; i++) {
            uint32_t card = sorted[i].ordinal;
            if (!(base[card / 64] & (UINT64_C(1) << (card % 64)))) continue;
            offset_scale = fmax(offset_scale, fabs(sorted[i].plane.d));
            sorted[kept++] = sorted[i];
        }
        uint32_t* ordinals = (uint32_t*)(base + words);
        recovery_plane_block_t* blocks = (recovery_plane_block_t*)((char*)ordinals + ordinal_bytes);
        recovery_plane_partition(sorted, 0, (uint32_t)kept, offset_scale,
                                 blocks, &out->blocks, ordinals);
        out->bytes = allocation;
        out->block_members = (uint32_t)kept;
    }
    free(sorted);
    return out;
}

static inline void recovery_plane_activate(uint64_t* mask,
    const recovery_plane_member_t* values, size_t begin, size_t end) {
    for (size_t i=begin;i<end;++i)
        mask[values[i].ordinal/64] |= UINT64_C(1) << (values[i].ordinal%64);
}

static inline void recovery_plane_candidates(const recovery_plane_index_t* index,
    const alea_ray_t* ray, double lo, double hi, uint64_t* mask) {
    memcpy(mask, recovery_plane_base(index), index->words*sizeof(*mask));
    const recovery_plane_member_t* values = recovery_plane_members(index);
    for (uint32_t i=0;i<index->groups;++i) {
        const recovery_plane_group_t* g = &index->data[i];
        const size_t end = g->first + g->count;
        /* These bounds rule out overflow/NaN in every kernel evaluation.
         * Uncertain groups stay active; parallel rays need no special proof. */
        double scale = fabs(g->a*ray->ox) + fabs(g->b*ray->oy) + fabs(g->c*ray->oz) +
            fmax(fabs(values[g->first].d), fabs(values[end-1].d));
        double direction = fabs(g->a*ray->dx) + fabs(g->b*ray->dy) + fabs(g->c*ray->dz);
        if (!(scale < DBL_MAX/4) || !(direction < DBL_MAX/4) || !isfinite(lo) || !isfinite(hi) || hi<lo) {
            recovery_plane_activate(mask, values, g->first, end);
            continue;
        }
        /* With a fixed normal and local ray, the established rounded kernel
         * is monotone in d. Search with that very kernel, including endpoint
         * equality, rather than approximate roots or geometric tolerances. */
        alea_plane_data_t p = {g->a, g->b, g->c, values[g->first].d};
        double first;
        RAY_PROFILE_COUNT(plane_index_groups);
        if (!ray_intersect_plane(ray, &p, &first, NULL, NULL, NULL)) {
#ifdef ALEA_RAY_PROFILE
            if (alea_ray_profile.enabled) alea_ray_profile.plane_index_pruned += g->count;
#endif
            continue;
        }
        p.d=values[end-1].d;
        double last;
        ray_intersect_plane(ray, &p, &last, NULL, NULL, NULL);
        bool reverse = first > last;
        const double smallest = reverse ? last : first;
        const double largest = reverse ? first : last;
        if (largest < lo || smallest > hi) {
#ifdef ALEA_RAY_PROFILE
            if (alea_ray_profile.enabled) alea_ray_profile.plane_index_pruned += g->count;
#endif
            continue;
        }
        if (smallest >= lo && largest <= hi) {
            recovery_plane_activate(mask, values, g->first, end);
            continue;
        }
        size_t bounds[2];
        for (int upper=0;upper<2;++upper) {
            size_t a=0, b=g->count;
            while (a<b) {
                size_t middle=a+(b-a)/2;
                size_t at=g->first+(reverse ? g->count-1-middle : middle);
                p.d=values[at].d;
                double t;
                ray_intersect_plane(ray, &p, &t, NULL, NULL, NULL);
                if (upper ? t<=hi : t<lo) a=middle+1;
                else b=middle;
            }
            bounds[upper]=a;
        }
        size_t begin = g->first + (reverse ? g->count-bounds[1] : bounds[0]);
        size_t finish = g->first + (reverse ? g->count-bounds[0] : bounds[1]);
        recovery_plane_activate(mask, values, begin, finish);
#ifdef ALEA_RAY_PROFILE
        if (alea_ray_profile.enabled) alea_ray_profile.plane_index_pruned += g->count-(finish-begin);
#endif
    }
}
/* The wider construction window amortizes the additional block tests. */
static inline void recovery_plane_root_candidates(const recovery_plane_index_t* p,
        const alea_ray_t* ray, double lo, double hi, uint64_t* mask) {
    recovery_plane_candidates(p, ray, lo, hi, mask);
    if (p->blocks)
        recovery_plane_prune_blocks(recovery_plane_blocks(p), p->blocks,
                                    recovery_plane_block_ordinals(p), ray, lo, hi, mask);
}
#endif
