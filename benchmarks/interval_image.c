// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
// SPDX-License-Identifier: MPL-2.0
#include "raycast/regional/interval_internal.h"
#include "interval_image.h"
#include <stdio.h>

/* Native, same-executable scratch image. The benchmark launcher verifies its
 * digest and binds it to the exact snapshot, options and executable. The
 * fixed arrays and all vectors below are immutable after preparation. */
#define IMAGE_VECTORS(X) X(nodes) X(patches) X(packets) X(extra_axes) X(tiles) \
    X(point_items) X(point_nodes) X(decisions) X(wide) X(boxes) X(surface_cards) X(patch_cards) X(guards.pairs)
#define IMAGE_ARRAYS(X) X(universe_roots,universes) X(domains,universes) \
    X(point_grids,universes) X(cell_cards,cells) X(planes,planes) X(quadrics,quadrics)

typedef struct {
    uint64_t magic, planes, quadrics, selective, cells, universes, primitives;
    alea_bbox_t domain, envelope;
    struct exterior_region exterior[6];
    double margin;
    uint32_t exterior_ready;
} interval_header_t;

typedef struct {
    FILE* out;
    const alea_mapped_file_t* in;
    size_t pos;
} interval_io_t;

/* Every section starts with its item count and is padded to eight bytes. */
static void* image_section(interval_io_t* io, const char* name, void* data,
                           size_t* count, size_t width, bool fixed) {
    uint64_t n = *count;
    if (io->in) {
        if (io->in->size - io->pos < sizeof(n)) return NULL;
        memcpy(&n,io->in->data + io->pos,sizeof(n));
        if ((fixed && n != *count) || n > UINT32_MAX) return NULL;
    }
    if (n > (SIZE_MAX - 7) / width) return NULL;
    size_t bytes = n * width, pad = (8 - bytes % 8) % 8;
    uint64_t zero = 0;
    if (io->out) {
        if (fwrite(&n,sizeof(n),1,io->out) != 1 ||
            (bytes && fwrite(data,1,bytes,io->out) != bytes) ||
            fwrite(&zero,1,pad,io->out) != pad) return NULL;
        if (getenv("ALEA_INTERVAL_PROFILE")) fprintf(stderr,"[interval-table] %s %zu bytes\n",name,bytes);
        return io;
    }
    io->pos += sizeof(n);
    if (bytes + pad > io->in->size - io->pos) return NULL;
    const char* base = io->in->data + io->pos;
    io->pos += bytes + pad;
    *count = n;
    return (void*)base;
}

static int image_tables(interval_index_t* idx, interval_io_t* io) {
    size_t universes = idx->sys->universes.count + 1, cells = idx->sys->cells.count + 1;
    size_t planes = idx->plane_count + 1, quadrics = idx->quadric_count + 1;
#define VECTOR(field) { \
    void* p = image_section(io,#field,idx->field.data,&idx->field.count,sizeof(*idx->field.data),false); \
    if (!p) return -1; \
    if (io->in) { idx->field.data = p; idx->field.capacity = idx->field.count; } }
    IMAGE_VECTORS(VECTOR)
#undef VECTOR
#define ARRAY(field,count) { \
    void* p = image_section(io,#field,idx->field,&count,sizeof(*idx->field),true); \
    if (!p) return -1; \
    if (io->in) idx->field = p; }
    IMAGE_ARRAYS(ARRAY)
#undef ARRAY
    return 0;
}

static void detach_image(interval_index_t* idx) {
#define CLEAR_VEC(field) idx->field.data = NULL;
    IMAGE_VECTORS(CLEAR_VEC)
#undef CLEAR_VEC
#define CLEAR_PTR(field,count) idx->field = NULL;
    IMAGE_ARRAYS(CLEAR_PTR)
#undef CLEAR_PTR
    alea_file_unmap(&idx->mapping);
}

#define INTERVAL_IMAGE_MAGIC UINT64_C(0x414c454149580006)
int interval_index_save(interval_index_t* idx, const char* path) {
    if (!idx->domains) return -1;
    interval_header_t h = {.magic=INTERVAL_IMAGE_MAGIC,.planes=idx->plane_count,
        .quadrics=idx->quadric_count,.selective=idx->selective,.domain=idx->domain,
        .cells=idx->sys->cells.count,.universes=idx->sys->universes.count,
        .primitives=idx->sys->primitives.count,.envelope=idx->envelope,
        .margin=idx->proof_margin,.exterior_ready=idx->exterior_ready};
    memcpy(h.exterior,idx->exterior,sizeof(h.exterior));
    FILE* f = fopen(path,"wbx");
    if (!f) return -1;
    interval_io_t io = {.out=f};
    int rc = fwrite(&h,sizeof(h),1,f) != 1 || image_tables(idx,&io);
    if (fclose(f)) rc = 1;
    return rc ? -1 : 0;
}

interval_index_t* interval_index_restore(alea_system_t* sys, const alea_bbox_t* domain,
                                          const char* path) {
    interval_index_t* idx = calloc(1,sizeof(*idx));
    if (!idx) return NULL;
    idx->sys = sys;
    idx->geometry_generation = alea_system_geometry_generation(sys);
    idx->material_generation = alea_system_material_generation(sys);
    if (!alea_file_map(path,&idx->mapping)) {
        free(idx);
        return NULL;
    }
    idx->detach_mapping = detach_image;
    interval_header_t h;
    if (idx->mapping.size < sizeof(h)) goto fail;
    memcpy(&h,idx->mapping.data,sizeof(h));
    if (!domain || h.magic != INTERVAL_IMAGE_MAGIC ||
        h.cells != sys->cells.count || h.universes != sys->universes.count ||
        h.primitives != sys->primitives.count || h.planes + h.quadrics > h.primitives ||
        h.planes > h.primitives || h.quadrics > h.primitives ||
        h.selective != (getenv("ALEA_INTERVAL_SELECTIVE") != NULL) ||
        memcmp(domain,&h.domain,sizeof(*domain))) goto fail;
    idx->domain=h.domain; idx->envelope=h.envelope; idx->proof_margin=h.margin;
    idx->exterior_ready=h.exterior_ready; idx->selective=h.selective;
    idx->plane_count=h.planes; idx->quadric_count=h.quadrics;
    memcpy(idx->exterior,h.exterior,sizeof(h.exterior));
    interval_io_t io = {.in=&idx->mapping,.pos=sizeof(h)};
    if (image_tables(idx,&io) || io.pos != idx->mapping.size) goto fail;
    if (alea_regional_reserve_workers(idx, 1)) goto fail;
    return idx;
fail:
    alea_prepared_queries_destroy(idx);
    return NULL;
}
