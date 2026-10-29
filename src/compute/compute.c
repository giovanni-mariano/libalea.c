// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
// SPDX-License-Identifier: MPL-2.0

#include "compute_internal.h"
#include "alea.h"
#include "core/alea_system.h"
#include "core/alea_transform.h"
#include "util/alea_parallel.h"
#include <float.h>
#include <limits.h>
#include <stdlib.h>

struct alea_compute {
    alea_compute_backend_t backend;
    alea_system_t* cpu;
    void* device;
    int flat, lattice;
};

/* Validate before compiling CSG or constructing render spatial programs. */
int alea_compute_validate_rect_lattices(const alea_system_t* sys) {
    for (size_t i = 0; i < sys->cells.count; ++i) {
        const alea_cell_entry_t* c = &sys->cells.data[i];
        if (!c->lat_type) continue;
        if (c->lat_type != 1 || c->fill_transform > 0) {
            alea_set_error_detail(ALEA_ERR_UNSUPPORTED,
                "Compute requires rectangular lattices without lattice-local fill transforms (cell %d)", c->mc_cell_id);
            return -1;
        }
        uint64_t count = 1;
        for (int a = 0; a < 3; ++a) {
            int64_t length = (int64_t)c->lat_fill_dims[2*a+1] - c->lat_fill_dims[2*a] + 1;
            if (length < 1 || length > INT32_MAX || count > UINT32_MAX / (uint64_t)length ||
                !isfinite(c->lat_pitch[a]) || c->lat_pitch[a] < 0 ||
                !isfinite(c->lat_lower_left[a]) ||
                ((a < 2 || length > 1) && c->lat_pitch[a] <= 0)) {
                alea_set_error_detail(ALEA_ERR_INVALID_STATE, "Invalid compute lattice bounds or pitch (cell %d)", c->mc_cell_id);
                return -1;
            }
            count *= (uint64_t)length;
        }
        if (!c->lat_fill || c->lat_fill_count != (c->lat_fill_repeating ? 1 : count)) {
            alea_set_error_detail(ALEA_ERR_INVALID_STATE, "Invalid compute lattice fill table (cell %d)", c->mc_cell_id);
            return -1;
        }
    }
    return 0;
}

int alea_compute_available(alea_compute_backend_t backend) {
    if (backend == ALEA_COMPUTE_CPU) return 1;
#ifdef ALEA_USE_METAL
    if (backend == ALEA_COMPUTE_METAL) return alea_compute_metal_available();
#endif
    return 0;
}

alea_compute_t* alea_compute_create(const alea_system_t* sys,
                                   alea_compute_backend_t backend,
                                   alea_compute_precision_t precision) {
    if (!sys) {
        alea_set_error_detail(ALEA_ERR_NULL_ARG, "Compute requires a geometry");
        return NULL;
    }
    if ((backend != ALEA_COMPUTE_CPU && backend != ALEA_COMPUTE_METAL) ||
        (precision != ALEA_COMPUTE_FP64 && precision != ALEA_COMPUTE_FP32)) {
        alea_set_error_detail(ALEA_ERR_INVALID_ARG, "Invalid compute backend or precision");
        return NULL;
    }
    if ((backend == ALEA_COMPUTE_CPU && precision != ALEA_COMPUTE_FP64) ||
        (backend == ALEA_COMPUTE_METAL && precision != ALEA_COMPUTE_FP32)) {
        alea_set_error_detail(ALEA_ERR_UNSUPPORTED,
            "CPU requires FP64; Metal requires explicit approximate FP32");
        return NULL;
    }
    if (!alea_compute_available(backend)) {
        alea_set_error_detail(ALEA_ERR_UNSUPPORTED, "Requested compute backend is unavailable");
        return NULL;
    }
    if (sys->cells.count > UINT32_MAX) {
        alea_set_error_detail(ALEA_ERR_OVERFLOW, "Too many compute cells");
        return NULL;
    }
    int flat = 1, lattice = 0;
    for (size_t i = 0; i < sys->cells.count; ++i) {
        const alea_cell_entry_t* c = &sys->cells.data[i];
        if (c->lat_type) lattice = 1;
        if (c->universe_id != 0 || c->fill_universe > 0 || c->lat_type != 0) flat = 0;
        if (c->root_node_id >= sys->nodes.count) {
            alea_set_error_detail(ALEA_ERR_INVALID_ID, "Invalid root for compute cell %d", c->mc_cell_id);
            return NULL;
        }
    }
    if (backend == ALEA_COMPUTE_METAL && lattice && alea_compute_validate_rect_lattices(sys)) return NULL;
    alea_compute_t* ctx = calloc(1, sizeof(*ctx));
    if (!ctx) {
        alea_set_error_detail(ALEA_ERR_OUT_OF_MEMORY, "Allocating compute context");
        return NULL;
    }
    ctx->backend = backend;
    ctx->flat = flat;
    ctx->lattice = backend == ALEA_COMPUTE_CPU ? lattice : 0;
    if (backend == ALEA_COMPUTE_CPU) {
        ctx->cpu = alea_clone(sys);
        if (!ctx->cpu) {
            alea_set_error_detail(ALEA_ERR_OUT_OF_MEMORY, "Copying compute geometry");
            free(ctx);
            return NULL;
        }
    }
#ifdef ALEA_USE_METAL
    else {
        ctx->device = alea_compute_metal_create(sys);
        if (!ctx->device) { free(ctx); return NULL; }
    }
#endif
    return ctx;
}

void alea_compute_destroy(alea_compute_t* ctx) {
    if (!ctx) return;
    alea_destroy(ctx->cpu);
#ifdef ALEA_USE_METAL
    if (ctx->device) alea_compute_metal_destroy(ctx->device);
#endif
    free(ctx);
}

typedef struct {
    const alea_system_t* sys;
    const double* xyz;
    alea_compute_result_t* results;
    int hierarchical;
} cpu_batch_t;

/* Independent FP64 element lookup. Empty/nonpositive fills are gaps. Active
 * axes use floor (the upper element owns a shared face). As in core queries,
 * finite singleton Y/Z windows include their upper face; untiled Z is unshifted
 * in center coordinates. The enclosing fill CSG supplies ancestor clipping. */
static int cpu_lattice_point(const alea_cell_entry_t* c, double p[3], int* child) {
    int64_t relative[3];
    int outside = 0;
    for (int a = 0; a < 3; ++a) {
        int64_t length = (int64_t)c->lat_fill_dims[2*a+1] - c->lat_fill_dims[2*a] + 1;
        double pitch = c->lat_pitch[a], lower = c->lat_lower_left[a];
        int tiled = a == 0 || length > 1 || (c->lat_fill_repeating && (a != 2 || pitch > 0));
        if (!c->lat_fill_repeating && c->lat_outer_universe <= 0 && length == 1 && pitch > 0 &&
            (p[a] < lower || p[a] > lower + pitch)) return 0;
        double index = tiled ? floor((p[a] - lower) / pitch) : 0;
        if (!isfinite(index) || index < INT32_MIN || index > INT32_MAX) return -1;
        relative[a] = (int64_t)index;
        int64_t logical = relative[a] + c->lat_fill_dims[2*a];
        if (logical < INT32_MIN || logical > INT32_MAX) return -1;
        outside |= relative[a] < 0 || relative[a] >= length;
        double shift = c->lat_fill_zero_element_coords ? logical * pitch :
            a == 2 && !tiled ? 0 : lower + (index + .5) * pitch;
        p[a] -= shift;
        if (!isfinite(p[a])) return -1;
    }
    if (c->lat_fill_repeating) *child = c->lat_fill[0];
    else if (outside) *child = c->lat_outer_universe;
    else {
        size_t ny = (size_t)((int64_t)c->lat_fill_dims[3] - c->lat_fill_dims[2] + 1);
        size_t nz = (size_t)((int64_t)c->lat_fill_dims[5] - c->lat_fill_dims[4] + 1);
        size_t at = ((size_t)relative[0] * ny + (size_t)relative[1]) * nz + (size_t)relative[2];
        *child = c->lat_fill[at];
    }
    return *child > 0;
}

/* Reference enumeration: count terminal occurrences, not fill containers.
 * Keep this independent of Metal packing and of selected-owner traversal. */
static int cpu_classify_universe(const alea_system_t* sys, int universe,
                                double x, double y, double z, unsigned depth,
                                alea_compute_result_t* result) {
    if (depth >= 32 || alea_interrupted()) return -1;
    const alea_universe_t* u = alea_get_universe(sys, universe);
    if (!u) return -1;
    for (size_t i = 0; i < u->cell_indices.count; ++i) {
        if ((i & 1023) == 0 && alea_interrupted()) return -1;
        const alea_cell_entry_t* c = &sys->cells.data[u->cell_indices.data[i]];
        if (c->lat_type) {
            double p[3] = {x, y, z}; int child;
            int found = cpu_lattice_point(c, p, &child);
            if (found < 0) return -1;
            if (found && cpu_classify_universe(sys, child, p[0], p[1], p[2], depth + 1, result)) return -1;
            continue;
        }
        if (!alea_point_inside(sys, c->root_node_id, x, y, z)) continue;
        if (c->fill_universe > 0) {
            double lx = x, ly = y, lz = z;
            if (c->fill_transform > 0) {
                const alea_transform_t* tr = alea_get_transform(sys, c->fill_transform);
                alea_matrix_t matrix;
                if (!tr || !alea_matrix_from_mcnp(&matrix, tr->cosines,
                        tr->value_count, false) || !alea_matrix_invert(&matrix)) return -1;
                alea_matrix_transform_point_inverse(&matrix, &lx, &ly, &lz);
                if (!isfinite(lx) || !isfinite(ly) || !isfinite(lz)) return -1;
            }
            if (cpu_classify_universe(sys, c->fill_universe, lx, ly, lz,
                                      depth + 1, result)) return -1;
        } else {
            if (result->owner_count == UINT32_MAX - 1) return -1;
            ++result->owner_count;
            result->cell_id = c->mc_cell_id;
            result->material_id = c->material_id;
        }
    }
    return 0;
}

static int cpu_range(void* arg, size_t worker, size_t begin, size_t end) {
    (void)worker;
    cpu_batch_t* batch = arg;
    for (size_t i = begin; i < end; ++i) {
        if (alea_interrupted()) return -1;
        const double* p = batch->xyz + 3 * i;
        alea_compute_result_t r = {-1, -1, 0};
        if (batch->hierarchical) {
            if (cpu_classify_universe(batch->sys, 0, p[0], p[1], p[2], 0, &r)) return -1;
        } else for (size_t j = 0; j < batch->sys->cells.count; ++j) {
            const alea_cell_entry_t* c = &batch->sys->cells.data[j];
            if (alea_point_inside(batch->sys, c->root_node_id, p[0], p[1], p[2])) {
                ++r.owner_count;
                r.cell_id = c->mc_cell_id;
                r.material_id = c->material_id;
            }
        }
        if (r.owner_count != 1) r.cell_id = r.material_id = -1;
        batch->results[i] = r;
    }
    return 0;
}

int alea_compute_classify(alea_compute_t* ctx, const double* xyz,
                          size_t count, alea_compute_result_t* results) {
    if (!ctx || (count && (!xyz || !results))) {
        alea_set_error_detail(ALEA_ERR_NULL_ARG, "Invalid compute batch storage");
        return -1;
    }
    if (count > SIZE_MAX / (3 * sizeof(double)) ||
        count > SIZE_MAX / sizeof(*results) ||
        (ctx->backend == ALEA_COMPUTE_METAL && count > UINT32_MAX)) {
        alea_set_error_detail(ALEA_ERR_OVERFLOW, "Compute batch is too large");
        return -1;
    }
    if (count == 0) return 0;
    if (ctx->lattice) {
        if (alea_compute_validate_rect_lattices(ctx->cpu)) return -1;
        ctx->lattice = 0;
    }
    ALEA_CHECK_INTERRUPTED(-1);
    for (size_t i = 0; i < count * 3; ++i) {
        if (!isfinite(xyz[i]) ||
            (ctx->backend == ALEA_COMPUTE_METAL && fabs(xyz[i]) > FLT_MAX)) {
            alea_set_error_detail(ALEA_ERR_INVALID_ARG, "Compute coordinate %zu is out of range", i);
            return -1;
        }
    }
#ifdef ALEA_USE_METAL
    if (ctx->backend == ALEA_COMPUTE_METAL)
        return alea_compute_metal_classify(ctx->device, xyz, count, results);
#endif
    if (!ctx->flat && !ctx->cpu->universe_index_built && alea_build_universe_index(ctx->cpu))
        return -1;
    cpu_batch_t batch = {ctx->cpu, xyz, results, !ctx->flat};
    alea_parallel_status_t status = alea_parallel_for(count, 256, 0,
        ALEA_PARALLEL_STATIC_BLOCK, cpu_range, &batch, NULL);
    if (status != ALEA_PARALLEL_OK) {
        alea_set_error_detail(alea_interrupted() ? ALEA_ERR_INTERRUPTED : ALEA_ERR_INVALID_STATE,
                             "CPU compute: %s", alea_parallel_status_string(status));
        return -1;
    }
    return 0;
}

/* Ordered CPU reference, independent of Metal packing/candidate traversal. */
static int cpu_select_universe(const alea_system_t* sys,int universe,const double p[3],
    unsigned depth,alea_compute_selection_t* result) {
    if(depth>=32 || alea_interrupted())return -1;
    const alea_universe_t* u=alea_get_universe(sys,universe);
    if(!u)return 0;
    for(size_t i=0;i<u->cell_indices.count;++i) {
        if(alea_interrupted())return -1;
        const alea_cell_entry_t* c=&sys->cells.data[u->cell_indices.data[i]];
        double local[3]={p[0],p[1],p[2]};int child=c->fill_universe;
        if(c->lat_type) {
            double checked[3]={p[0],p[1],p[2]};int checked_child;
            if(cpu_lattice_point(c,checked,&checked_child)<0)return -1;
            alea_lattice_location_t loc;
            int found=alea_lattice_locate_point(sys,c,p[0],p[1],p[2],&loc);
            if(found<0)return -1;
            if(!found)continue;
            local[0]-=loc.ox;local[1]-=loc.oy;local[2]-=loc.oz;child=loc.fill_universe;
        } else {
            if(!alea_point_inside(sys,c->root_node_id,p[0],p[1],p[2]))continue;
            if(child>0 && c->fill_transform>0) {
                const alea_transform_t* tr=alea_get_transform(sys,c->fill_transform);alea_matrix_t matrix;
                if(!tr || !alea_matrix_from_mcnp(&matrix,tr->cosines,tr->value_count,false) ||
                   !alea_matrix_invert(&matrix))return -1;
                alea_matrix_transform_point_inverse(&matrix,&local[0],&local[1],&local[2]);
            }
        }
        if(!isfinite(local[0]) || !isfinite(local[1]) || !isfinite(local[2]))return -1;
        *result=(alea_compute_selection_t){c->mc_cell_id,c->material_id,(int32_t)depth,
            c->lat_type || child>0?ALEA_COMPUTE_UNDEFINED_FILL:ALEA_COMPUTE_FOUND};
        if(child>0) {
            alea_compute_selection_t descendant={-1,-1,-1,ALEA_COMPUTE_NOT_FOUND};
            if(cpu_select_universe(sys,child,local,depth+1,&descendant))return -1;
            if(descendant.status!=ALEA_COMPUTE_NOT_FOUND)*result=descendant;
        }
        return 0;
    }
    return 0;
}
int alea_compute_select(alea_compute_t* ctx,const double* xyz,size_t count,alea_compute_selection_t* results) {
    if(!ctx || (count && (!xyz || !results))) {
        alea_set_error_detail(ALEA_ERR_NULL_ARG,"Invalid selection storage");return -1;
    }
    if(count>SIZE_MAX/(3*sizeof(double)) || count>SIZE_MAX/sizeof(*results) ||
       (ctx->backend==ALEA_COMPUTE_METAL && count>UINT32_MAX)) {
        alea_set_error_detail(ALEA_ERR_OVERFLOW,"Selection batch is too large");return -1;
    }
    if(!count)return 0;
    ALEA_CHECK_INTERRUPTED(-1);
    for(size_t i=0;i<count*3;++i)if(!isfinite(xyz[i]) ||
        (ctx->backend==ALEA_COMPUTE_METAL && fabs(xyz[i])>FLT_MAX)) {
        alea_set_error_detail(ALEA_ERR_INVALID_ARG,"Selection coordinate out of range");return -1;
    }
#ifdef ALEA_USE_METAL
    if(ctx->backend==ALEA_COMPUTE_METAL)return alea_compute_metal_select(ctx->device,xyz,count,results);
#endif
    if(ctx->lattice) {
        if(alea_compute_validate_rect_lattices(ctx->cpu))return -1;
        ctx->lattice=0;
    }
    if(!ctx->cpu->universe_index_built && alea_build_universe_index(ctx->cpu))return -1;
    for(size_t i=0;i<count;++i) {
        results[i]=(alea_compute_selection_t){-1,-1,-1,ALEA_COMPUTE_NOT_FOUND};
        if(cpu_select_universe(ctx->cpu,0,xyz+3*i,0,&results[i])) {
            alea_set_error_detail(alea_interrupted()?ALEA_ERR_INTERRUPTED:ALEA_ERR_OVERFLOW,
                                  "CPU selection failed at point %zu",i);return -1;
        }
    }
    return 0;
}

static int compute_render(alea_compute_t* ctx, const render_config_t* cfg,
                        const render_camera_t* cam, render_framebuffer_t* fb,
                        alea_compute_render_progress_t progress, void* user) {
    if (!ctx || !cfg || !cam || !fb || !fb->color || !fb->cell_id) {
        alea_set_error_detail(ALEA_ERR_NULL_ARG, "Invalid compute render arguments");
        return -1;
    }
    if (fb->width < 1 || fb->height < 1 || fb->width > 16384 || fb->height > 16384 ||
        cfg->width != fb->width || cfg->height != fb->height) {
        alea_set_error_detail(ALEA_ERR_INVALID_ARG, "Invalid or mismatched render dimensions");
        return -1;
    }
    ALEA_CHECK_INTERRUPTED(-1);
#ifdef ALEA_USE_METAL
    if (ctx->backend == ALEA_COMPUTE_METAL)
        return alea_compute_metal_render(ctx->device, cfg, cam, fb, progress, user);
#endif
    (void)user;
    if(progress){alea_set_error_detail(ALEA_ERR_UNSUPPORTED,"Progressive rendering requires a Metal ray context");return -1;}
    if (alea_prepare_query_acceleration(ctx->cpu)) return -1;
    return render_scene(ctx->cpu, cfg, cam, fb);
}

int alea_compute_render(alea_compute_t* ctx, const render_config_t* cfg,
        const render_camera_t* cam, render_framebuffer_t* fb) {
    return compute_render(ctx,cfg,cam,fb,NULL,NULL);
}
int alea_compute_render_progressive(alea_compute_t* ctx, const render_config_t* cfg,
        const render_camera_t* cam, render_framebuffer_t* fb,
        alea_compute_render_progress_t progress, void* user) {
    if(!progress){alea_set_error_detail(ALEA_ERR_NULL_ARG,"Progress callback is required");return -1;}
    return compute_render(ctx,cfg,cam,fb,progress,user);
}
