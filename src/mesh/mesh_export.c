// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/**
 * @file mesh_export.c
 * @brief Structured mesh sampling and export (Gmsh / VTK)
 *
 * Samples CSG geometry onto a hexahedral grid and writes mesh files
 * suitable for OpenSn (FromFileMeshGenerator) and other codes.
 */

#include "alea_mesh.h"
#include "alea.h"
#include "core/alea_system.h"
#include "core/alea_universe.h"
#include "primitives/bbox.h"
#include "raycast/raycast.h"
#include "util/compat.h"

#include <errno.h>
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include <stdint.h>
#ifdef _OPENMP
#include <omp.h>
#endif

/* ============================================================================
 * Config
 * ============================================================================ */

void alea_mesh_config_init(alea_mesh_config_t *cfg) {
    if (!cfg) return;
    memset(cfg, 0, sizeof(*cfg));
    cfg->nx = cfg->ny = cfg->nz = 10;
    cfg->format = ALEA_MESH_GMSH;
    cfg->void_material_id = 0;
    cfg->auto_pad = 0.01;
    cfg->sampling_mode = ALEA_MESH_SAMPLE_SUBCELL;
    cfg->subsamples_per_axis = 2;
    cfg->mixed_threshold = 0.0;
    cfg->target_error = 0.05;
    cfg->max_refine_depth = 3;
    cfg->max_samples_per_voxel = 32768;
    cfg->sampling_seed = UINT64_C(0x6a09e667f3bcc909);
    cfg->workers = 1;
    cfg->ray_grid_u = 4;
    cfg->ray_grid_v = 4;
    cfg->ray_directions = ALEA_MESH_RAY_XYZ;
    cfg->fields = ALEA_MESH_FIELD_MATERIAL_ID |
                  ALEA_MESH_FIELD_CELL_ID |
                  ALEA_MESH_FIELD_MIXED_FLAG |
                  ALEA_MESH_FIELD_DOMINANT_FRACTION |
                  ALEA_MESH_FIELD_SAMPLED_FRACTIONS |
                  ALEA_MESH_FIELD_SAMPLE_COUNT |
                  ALEA_MESH_FIELD_TIE_FLAG |
                  ALEA_MESH_FIELD_ESTIMATED_ERROR |
                  ALEA_MESH_FIELD_REFINEMENT_FLAG |
                  ALEA_MESH_FIELD_CELL_FRACTIONS;
}

/* ============================================================================
 * Helpers
 * ============================================================================ */

/** Build uniform node array: n+1 values from lo to hi */
static double *build_uniform_nodes(int n, double lo, double hi) {
    double *nodes = malloc((size_t)(n + 1) * sizeof(double));
    if (!nodes) return NULL;
    double step = (hi - lo) / n;
    for (int i = 0; i <= n; i++)
        nodes[i] = lo + i * step;
    return nodes;
}

/** Copy n+1 doubles from user array */
static double *copy_nodes(const double *src, int n) {
    size_t sz = (size_t)(n + 1) * sizeof(double);
    double *dst = malloc(sz);
    if (dst) memcpy(dst, src, sz);
    return dst;
}

static int checked_mul_size(size_t a, size_t b, size_t *out) {
    if (!out) return -1;
    if (a != 0 && b > SIZE_MAX / a) {
        alea_set_error_detail(ALEA_ERR_OVERFLOW,
                              "mesh dimension product overflows size_t");
        return -1;
    }
    *out = a * b;
    return 0;
}

static int checked_add_size(size_t a, size_t b, size_t *out) {
    if (!out) return -1;
    if (b > SIZE_MAX - a) {
        alea_set_error_detail(ALEA_ERR_OVERFLOW,
                              "mesh dimension sum overflows size_t");
        return -1;
    }
    *out = a + b;
    return 0;
}

static int checked_mesh_cell_count(int nx, int ny, int nz, size_t *out) {
    size_t xy = 0;
    if (checked_mul_size((size_t)nx, (size_t)ny, &xy) != 0) return -1;
    return checked_mul_size(xy, (size_t)nz, out);
}

static int checked_mesh_node_count(int nx, int ny, int nz, size_t *out) {
    size_t nx1 = 0, ny1 = 0, nz1 = 0, xy = 0;
    if (checked_add_size((size_t)nx, 1, &nx1) != 0) return -1;
    if (checked_add_size((size_t)ny, 1, &ny1) != 0) return -1;
    if (checked_add_size((size_t)nz, 1, &nz1) != 0) return -1;
    if (checked_mul_size(nx1, ny1, &xy) != 0) return -1;
    return checked_mul_size(xy, nz1, out);
}

static int validate_nodes(const double *nodes, int n, const char *axis) {
    if (!nodes) return 0;
    for (int i = 0; i <= n; i++) {
        if (!isfinite(nodes[i])) {
            alea_set_error_detail(ALEA_ERR_INVALID_ARG,
                                  "mesh %s node %d is not finite", axis, i);
            return -1;
        }
        if (i > 0 && !(nodes[i] > nodes[i - 1])) {
            alea_set_error_detail(ALEA_ERR_INVALID_ARG,
                                  "mesh %s nodes must be strictly increasing", axis);
            return -1;
        }
    }
    return 0;
}

static int validate_uniform_axis(double lo, double hi, const char *axis) {
    if (!isfinite(lo) || !isfinite(hi)) {
        alea_set_error_detail(ALEA_ERR_INVALID_ARG,
                              "mesh %s bounds are not finite", axis);
        return -1;
    }
    if (!(lo < hi)) {
        alea_set_error_detail(ALEA_ERR_INVALID_ARG,
                              "mesh %s bounds must satisfy min < max", axis);
        return -1;
    }
    return 0;
}

static int validate_config_basic(const alea_mesh_config_t *cfg) {
    const uint32_t known_fields = ALEA_MESH_FIELD_MATERIAL_ID |
                                  ALEA_MESH_FIELD_CELL_ID |
                                  ALEA_MESH_FIELD_MIXED_FLAG |
                                  ALEA_MESH_FIELD_DOMINANT_FRACTION |
                                  ALEA_MESH_FIELD_SAMPLED_FRACTIONS |
                                  ALEA_MESH_FIELD_SAMPLE_COUNT |
                                  ALEA_MESH_FIELD_TIE_FLAG |
                                  ALEA_MESH_FIELD_ESTIMATED_ERROR |
                                  ALEA_MESH_FIELD_REFINEMENT_FLAG |
                                  ALEA_MESH_FIELD_CELL_FRACTIONS;
    if (!cfg) {
        alea_set_error_detail(ALEA_ERR_NULL_ARG, "mesh config is NULL");
        return -1;
    }
    if (cfg->nx <= 0 || cfg->ny <= 0 || cfg->nz <= 0) {
        alea_set_error_detail(ALEA_ERR_INVALID_ARG,
                              "mesh dimensions must be positive");
        return -1;
    }
    if (cfg->sampling_mode != ALEA_MESH_SAMPLE_CENTER &&
        cfg->sampling_mode != ALEA_MESH_SAMPLE_CORNERS &&
        cfg->sampling_mode != ALEA_MESH_SAMPLE_SUBCELL &&
        cfg->sampling_mode != ALEA_MESH_SAMPLE_STRATIFIED &&
        cfg->sampling_mode != ALEA_MESH_SAMPLE_ADAPTIVE &&
        cfg->sampling_mode != ALEA_MESH_SAMPLE_RAY) {
        alea_set_error_detail(ALEA_ERR_INVALID_ARG,
                              "mesh sampling_mode is invalid");
        return -1;
    }
    if (cfg->subsamples_per_axis <= 0) {
        alea_set_error_detail(ALEA_ERR_INVALID_ARG,
                              "mesh subsamples_per_axis must be positive");
        return -1;
    }
    if (!isfinite(cfg->target_error) || cfg->target_error <= 0.0 ||
        cfg->target_error > 1.0) {
        alea_set_error_detail(ALEA_ERR_INVALID_ARG,
                              "mesh target_error must be finite and in (0,1]");
        return -1;
    }
    if (cfg->max_refine_depth < 0 || cfg->max_refine_depth > 20) {
        alea_set_error_detail(ALEA_ERR_INVALID_ARG,
                              "mesh max_refine_depth must be in [0,20]");
        return -1;
    }
    if (cfg->max_samples_per_voxel == 0) {
        alea_set_error_detail(ALEA_ERR_INVALID_ARG,
                              "mesh max_samples_per_voxel must be positive");
        return -1;
    }
    if (cfg->workers < 0) {
        alea_set_error_detail(ALEA_ERR_INVALID_ARG,
                              "mesh workers must be non-negative");
        return -1;
    }
    if (cfg->ray_grid_u <= 0 || cfg->ray_grid_v <= 0) {
        alea_set_error_detail(ALEA_ERR_INVALID_ARG,
                              "mesh ray grid dimensions must be positive");
        return -1;
    }
    if (cfg->ray_directions == 0 ||
        (cfg->ray_directions & ~ALEA_MESH_RAY_XYZ) != 0) {
        alea_set_error_detail(ALEA_ERR_INVALID_ARG,
                              "mesh ray direction mask is invalid");
        return -1;
    }
    if (!isfinite(cfg->mixed_threshold) || cfg->mixed_threshold < 0.0 ||
        cfg->mixed_threshold >= 1.0) {
        alea_set_error_detail(ALEA_ERR_INVALID_ARG,
                              "mesh mixed_threshold must be finite and in [0,1)");
        return -1;
    }
    if (!isfinite(cfg->auto_pad) || cfg->auto_pad < 0.0) {
        alea_set_error_detail(ALEA_ERR_INVALID_ARG,
                              "mesh auto_pad must be finite and non-negative");
        return -1;
    }
    if (cfg->bounds_mode != ALEA_MESH_BOUNDS_LEGACY &&
        cfg->bounds_mode != ALEA_MESH_BOUNDS_AUTO &&
        cfg->bounds_mode != ALEA_MESH_BOUNDS_EXPLICIT) {
        alea_set_error_detail(ALEA_ERR_INVALID_ARG,
                              "mesh bounds mode is invalid");
        return -1;
    }
    if (cfg->fields & ~known_fields) {
        alea_set_error_detail(ALEA_ERR_INVALID_ARG,
                              "mesh result field mask is invalid");
        return -1;
    }
    if (validate_nodes(cfg->x_nodes, cfg->nx, "X") != 0) return -1;
    if (validate_nodes(cfg->y_nodes, cfg->ny, "Y") != 0) return -1;
    if (validate_nodes(cfg->z_nodes, cfg->nz, "Z") != 0) return -1;
    return 0;
}

static int compare_ints(const void *lhs, const void *rhs) {
    const int a = *(const int *)lhs;
    const int b = *(const int *)rhs;
    return (a > b) - (a < b);
}

static void *mesh_alloc_array(size_t count, size_t element_size, int clear) {
    size_t bytes = 0;
    if (checked_mul_size(count, element_size, &bytes) != 0) return NULL;
    void *result = clear ? calloc(count, element_size) : malloc(bytes);
    if (!result) {
        alea_set_error_detail(ALEA_ERR_OUT_OF_MEMORY,
                              "mesh array allocation failed");
    }
    return result;
}

void alea_mesh_export_options_init(alea_mesh_export_options_t *options) {
    if (!options) return;
    options->fields = ALEA_MESH_EXPORT_MIXED_FLAG |
                      ALEA_MESH_EXPORT_DOMINANT_FRACTION |
                      ALEA_MESH_EXPORT_TIE_FLAG |
                      ALEA_MESH_EXPORT_SAMPLE_COUNT |
                      ALEA_MESH_EXPORT_ESTIMATED_ERROR |
                      ALEA_MESH_EXPORT_REFINEMENT_FLAG;
    options->max_fraction_materials = 64;
}

typedef struct {
    int *keys;
    unsigned char *used;
    size_t capacity;
    size_t count;
} mesh_int_set_t;

static size_t mesh_int_hash(int value) {
    uint32_t x = (uint32_t)value;
    x ^= x >> 16;
    x *= UINT32_C(0x7feb352d);
    x ^= x >> 15;
    x *= UINT32_C(0x846ca68b);
    x ^= x >> 16;
    return (size_t)x;
}

static void mesh_int_set_free(mesh_int_set_t *set) {
    free(set->keys);
    free(set->used);
    memset(set, 0, sizeof(*set));
}

static int mesh_int_set_rehash(mesh_int_set_t *set, size_t capacity) {
    int *keys = mesh_alloc_array(capacity, sizeof(*keys), 0);
    unsigned char *used = mesh_alloc_array(capacity, sizeof(*used), 1);
    if (!keys || !used) {
        free(keys); free(used);
        alea_set_error_detail(ALEA_ERR_OUT_OF_MEMORY,
                              "mesh failed to grow material set");
        return -1;
    }
    for (size_t i = 0; i < set->capacity; i++) {
        if (!set->used[i]) continue;
        size_t slot = mesh_int_hash(set->keys[i]) & (capacity - 1);
        while (used[slot]) slot = (slot + 1) & (capacity - 1);
        keys[slot] = set->keys[i];
        used[slot] = 1;
    }
    free(set->keys); free(set->used);
    set->keys = keys;
    set->used = used;
    set->capacity = capacity;
    return 0;
}

static int mesh_int_set_insert(mesh_int_set_t *set, int value) {
    if (set->capacity == 0 && mesh_int_set_rehash(set, 16) != 0) return -1;
    if ((set->count + 1) * 10 >= set->capacity * 7) {
        if (set->capacity > SIZE_MAX / 2 ||
            mesh_int_set_rehash(set, set->capacity * 2) != 0) return -1;
    }
    size_t slot = mesh_int_hash(value) & (set->capacity - 1);
    while (set->used[slot]) {
        if (set->keys[slot] == value) return 0;
        slot = (slot + 1) & (set->capacity - 1);
    }
    set->keys[slot] = value;
    set->used[slot] = 1;
    set->count++;
    return 0;
}

static int *mesh_int_set_sorted_array(const mesh_int_set_t *set, int *out_num) {
    if (set->count == 0 || set->count > (size_t)INT_MAX) return NULL;
    int *result = mesh_alloc_array(set->count, sizeof(*result), 0);
    if (!result) {
        alea_set_error_detail(ALEA_ERR_OUT_OF_MEMORY,
                              "mesh failed to materialize material table");
        return NULL;
    }
    size_t n = 0;
    for (size_t i = 0; i < set->capacity; i++)
        if (set->used[i]) result[n++] = set->keys[i];
    qsort(result, n, sizeof(*result), compare_ints);
    *out_num = (int)n;
    return result;
}

static void sample_owner_at_point(alea_system_t *sys,
                                  double x, double y, double z,
                                  int void_mat,
                                  int *out_material_id,
                                  int *out_cell_id) {
    int cell_idx = alea_identify_cell_at_point(sys, x, y, z);
    if (cell_idx < 0) {
        *out_material_id = void_mat;
        *out_cell_id = -1;
        return;
    }
    *out_material_id = sys->cells.data[cell_idx].material_id;
    *out_cell_id = sys->cells.data[cell_idx].mc_cell_id;
}

static void tally_material(int material_id,
                           int *materials,
                           int *counts,
                           int *num_materials) {
    for (int i = 0; i < *num_materials; i++) {
        if (materials[i] == material_id) {
            counts[i]++;
            return;
        }
    }
    materials[*num_materials] = material_id;
    counts[*num_materials] = 1;
    (*num_materials)++;
}

static void tally_owner(int material_id, int cell_id,
                        int *materials, int *cells, int *counts,
                        int *num_owners) {
    for (int i = 0; i < *num_owners; i++) {
        if (materials[i] == material_id && cells[i] == cell_id) {
            counts[i]++;
            return;
        }
    }
    materials[*num_owners] = material_id;
    cells[*num_owners] = cell_id;
    counts[*num_owners] = 1;
    (*num_owners)++;
}

static uint64_t mesh_mix64(uint64_t x) {
    x ^= x >> 30;
    x *= UINT64_C(0xbf58476d1ce4e5b9);
    x ^= x >> 27;
    x *= UINT64_C(0x94d049bb133111eb);
    return x ^ (x >> 31);
}

static double mesh_jitter(const alea_mesh_config_t *cfg, int i, int j, int k,
                          int ii, int jj, int kk, int axis) {
    uint64_t x = cfg->sampling_seed;
    x = mesh_mix64(x ^ (uint32_t)i);
    x = mesh_mix64(x ^ ((uint64_t)(uint32_t)j << 21));
    x = mesh_mix64(x ^ ((uint64_t)(uint32_t)k << 42));
    x = mesh_mix64(x ^ (uint32_t)ii ^ ((uint64_t)(uint32_t)jj << 20) ^
                   ((uint64_t)(uint32_t)kk << 40) ^ ((uint64_t)(uint32_t)axis << 52));
    return (double)(x >> 11) * 0x1.0p-53;
}

static int mesh_sample_voxel_materials(alea_system_t *sys,
                                       const alea_mesh_config_t *cfg,
                                       const double *xn,
                                       const double *yn,
                                       const double *zn,
                                       int i, int j, int k, int n,
                                       int *materials,
                                       int *counts,
                                       int *out_num_materials,
                                       int *owner_materials,
                                       int *owner_cells,
                                       int *owner_counts,
                                       int *out_num_owners,
                                       int *out_num_samples) {
    int nsamples = n * n * n;
    int n_materials = 0;
    int n_owners = 0;

    double x0 = xn[i], x1 = xn[i + 1];
    double y0 = yn[j], y1 = yn[j + 1];
    double z0 = zn[k], z1 = zn[k + 1];
    double dx = x1 - x0, dy = y1 - y0, dz = z1 - z0;
    double eps = 1e-9;

    memset(counts, 0, (size_t)nsamples * sizeof(int));
    memset(owner_counts, 0, (size_t)nsamples * sizeof(int));

    for (int kk = 0; kk < n; kk++) {
        for (int jj = 0; jj < n; jj++) {
            for (int ii = 0; ii < n; ii++) {
                double fx, fy, fz;
                if (cfg->sampling_mode == ALEA_MESH_SAMPLE_CENTER) {
                    fx = fy = fz = 0.5;
                } else if (cfg->sampling_mode == ALEA_MESH_SAMPLE_CORNERS) {
                    fx = (ii == 0) ? eps : 1.0 - eps;
                    fy = (jj == 0) ? eps : 1.0 - eps;
                    fz = (kk == 0) ? eps : 1.0 - eps;
                } else if (cfg->sampling_mode == ALEA_MESH_SAMPLE_STRATIFIED ||
                           cfg->sampling_mode == ALEA_MESH_SAMPLE_ADAPTIVE) {
                    fx = ((double)ii + mesh_jitter(cfg, i, j, k, ii, jj, kk, n)) /
                         (double)n;
                    fy = ((double)jj + mesh_jitter(cfg, i, j, k, jj, kk, ii,
                                                   n ^ 0x155)) / (double)n;
                    fz = ((double)kk + mesh_jitter(cfg, i, j, k, kk, ii, jj,
                                                   n ^ 0x2aa)) / (double)n;
                } else {
                    fx = ((double)ii + 0.5) / (double)n;
                    fy = ((double)jj + 0.5) / (double)n;
                    fz = ((double)kk + 0.5) / (double)n;
                }

                double x = x0 + fx * dx;
                double y = y0 + fy * dy;
                double z = z0 + fz * dz;
                int mat = cfg->void_material_id;
                int cell = -1;
                sample_owner_at_point(sys, x, y, z, cfg->void_material_id,
                                      &mat, &cell);
                tally_material(mat, materials, counts, &n_materials);
                tally_owner(mat, cell, owner_materials, owner_cells,
                            owner_counts, &n_owners);
            }
        }
    }

    *out_num_materials = n_materials;
    *out_num_owners = n_owners;
    *out_num_samples = nsamples;
    return 0;
}

static void mesh_fraction_stats(const int *materials,
                                const int *counts,
                                int num_materials,
                                int num_samples,
                                const alea_mesh_config_t *cfg,
                                const int *owner_materials,
                                const int *owner_cells,
                                const int *owner_counts,
                                int num_owners,
                                int *out_dominant_material,
                                int *out_dominant_cell,
                                unsigned char *out_mixed,
                                double *out_dominant_fraction,
                                uint8_t *out_tie_flags) {
    int max_count = 0;
    int dominant_material = cfg->void_material_id;
    for (int m = 0; m < num_materials; m++) {
        if (counts[m] > max_count ||
            (counts[m] == max_count && materials[m] < dominant_material)) {
            max_count = counts[m];
            dominant_material = materials[m];
        }
    }

    uint8_t ties = 0;
    int max_material_count_matches = 0;
    for (int m = 0; m < num_materials; m++)
        if (counts[m] == max_count) max_material_count_matches++;
    if (max_material_count_matches > 1) ties |= ALEA_MESH_TIE_MATERIAL;

    int dominant_cell = -1;
    int max_cell_count = 0;
    for (int o = 0; o < num_owners; o++) {
        if (owner_materials[o] != dominant_material) continue;
        if (owner_counts[o] > max_cell_count ||
            (owner_counts[o] == max_cell_count && owner_cells[o] < dominant_cell)) {
            max_cell_count = owner_counts[o];
            dominant_cell = owner_cells[o];
        }
    }
    int max_cell_count_matches = 0;
    for (int o = 0; o < num_owners; o++)
        if (owner_materials[o] == dominant_material &&
            owner_counts[o] == max_cell_count) max_cell_count_matches++;
    if (max_cell_count_matches > 1) ties |= ALEA_MESH_TIE_CELL;

    double dominant = num_samples > 0 ? (double)max_count / (double)num_samples : 1.0;
    *out_dominant_material = dominant_material;
    *out_dominant_cell = dominant_cell;
    *out_dominant_fraction = dominant;
    *out_mixed = (dominant < 1.0 - cfg->mixed_threshold) ? 1 : 0;
    *out_tie_flags = ties;
}

static double mesh_distribution_error(const int *old_materials,
                                      const int *old_counts, int old_n,
                                      int old_samples,
                                      const int *new_materials,
                                      const int *new_counts, int new_n,
                                      int new_samples) {
    double l1 = 0.0;
    for (int m = 0; m < old_n; m++) {
        int count = 0;
        for (int n = 0; n < new_n; n++)
            if (new_materials[n] == old_materials[m]) { count = new_counts[n]; break; }
        l1 += fabs((double)old_counts[m] / old_samples - (double)count / new_samples);
    }
    for (int n = 0; n < new_n; n++) {
        int found = 0;
        for (int m = 0; m < old_n; m++)
            if (old_materials[m] == new_materials[n]) { found = 1; break; }
        if (!found) l1 += (double)new_counts[n] / new_samples;
    }
    return 0.5 * l1;
}

static int mesh_cube_samples(int axis, uint32_t *out) {
    if (axis <= 0) return -1;
    uint64_t a = (uint64_t)axis;
    if (a > UINT64_MAX / a) return -1;
    uint64_t square = a * a;
    if (a > UINT64_MAX / square) return -1;
    uint64_t cube = square * a;
    if (cube > UINT32_MAX) return -1;
    *out = (uint32_t)cube;
    return 0;
}

static int ensure_fraction_capacity(alea_mesh_material_fraction_t **fractions,
                                    size_t *capacity,
                                    size_t needed) {
    if (needed <= *capacity) return 0;

    size_t new_cap = (*capacity == 0) ? 1024 : *capacity;
    while (new_cap < needed) {
        if (new_cap > SIZE_MAX / 2) {
            alea_set_error_detail(ALEA_ERR_OVERFLOW,
                                  "mesh fraction storage overflows size_t");
            return -1;
        }
        new_cap *= 2;
    }

    alea_mesh_material_fraction_t *tmp = NULL;
    size_t bytes = 0;
    if (checked_mul_size(new_cap, sizeof(alea_mesh_material_fraction_t),
                         &bytes) != 0) return -1;
    tmp = realloc(*fractions, bytes);
    if (!tmp) {
        alea_set_error_detail(ALEA_ERR_OUT_OF_MEMORY,
                              "mesh failed to grow material fraction storage");
        return -1;
    }

    *fractions = tmp;
    *capacity = new_cap;
    return 0;
}

static int ensure_cell_fraction_capacity(alea_mesh_cell_fraction_t **fractions,
                                         size_t *capacity, size_t needed) {
    if (needed <= *capacity) return 0;
    size_t new_cap = (*capacity == 0) ? 1024 : *capacity;
    while (new_cap < needed) {
        if (new_cap > SIZE_MAX / 2) {
            alea_set_error_detail(ALEA_ERR_OVERFLOW,
                                  "mesh cell-fraction storage overflows size_t");
            return -1;
        }
        new_cap *= 2;
    }
    size_t bytes = 0;
    if (checked_mul_size(new_cap, sizeof(**fractions), &bytes) != 0) return -1;
    alea_mesh_cell_fraction_t *tmp = realloc(*fractions, bytes);
    if (!tmp) {
        alea_set_error_detail(ALEA_ERR_OUT_OF_MEMORY,
                              "mesh failed to grow cell-fraction storage");
        return -1;
    }
    *fractions = tmp;
    *capacity = new_cap;
    return 0;
}

typedef struct {
    int cell_id;
    int material_id;
    double value;
} mesh_ray_component_t;

typedef struct {
    mesh_ray_component_t *data;
    uint32_t count;
    uint32_t capacity;
} mesh_ray_voxel_accum_t;

static int mesh_ray_accum_add(mesh_ray_voxel_accum_t *accum, int cell_id,
                              int material_id, double value) {
    if (!(value > 0.0)) return 0;
    for (uint32_t i = 0; i < accum->count; i++) {
        if (accum->data[i].cell_id == cell_id &&
            accum->data[i].material_id == material_id) {
            accum->data[i].value += value;
            return 0;
        }
    }
    if (accum->count == accum->capacity) {
        uint32_t capacity = accum->capacity ? accum->capacity * 2u : 4u;
        if (capacity < accum->capacity) return -1;
        mesh_ray_component_t *data = realloc(
            accum->data, (size_t)capacity * sizeof(*data));
        if (!data) return -1;
        accum->data = data;
        accum->capacity = capacity;
    }
    accum->data[accum->count++] = (mesh_ray_component_t){
        .cell_id = cell_id, .material_id = material_id, .value = value};
    return 0;
}

static int mesh_ray_component_compare(const void *lhs, const void *rhs) {
    const mesh_ray_component_t *a = lhs, *b = rhs;
    if (a->material_id != b->material_id)
        return (a->material_id > b->material_id) -
               (a->material_id < b->material_id);
    return (a->cell_id > b->cell_id) - (a->cell_id < b->cell_id);
}

static int mesh_ray_find_bin(const double *nodes, int count, double x) {
    if (x <= nodes[0]) return 0;
    if (x >= nodes[count]) return count - 1;
    int lo = 0, hi = count;
    while (lo + 1 < hi) {
        int mid = lo + (hi - lo) / 2;
        if (x < nodes[mid]) hi = mid;
        else lo = mid;
    }
    return lo;
}

static size_t mesh_ray_voxel_index(int axis, int along, int a, int b,
                                   int nx, int ny) {
    if (axis == 0) return (size_t)b * (size_t)nx * (size_t)ny +
                           (size_t)a * (size_t)nx + (size_t)along;
    if (axis == 1) return (size_t)b * (size_t)nx * (size_t)ny +
                           (size_t)along * (size_t)nx + (size_t)a;
    return (size_t)along * (size_t)nx * (size_t)ny +
           (size_t)b * (size_t)nx + (size_t)a;
}

static int mesh_ray_accumulate_segment(mesh_ray_voxel_accum_t *accums,
                                       const double *nodes, int count,
                                       int axis, int a, int b, int nx, int ny,
                                       double begin, double end,
                                       int cell_id, int material_id) {
    if (begin < nodes[0]) begin = nodes[0];
    if (end > nodes[count]) end = nodes[count];
    if (!(end > begin)) return 0;
    int bin = mesh_ray_find_bin(nodes, count, begin);
    while (bin < count && nodes[bin] < end) {
        double lo = begin > nodes[bin] ? begin : nodes[bin];
        double hi = end < nodes[bin + 1] ? end : nodes[bin + 1];
        if (hi > lo) {
            size_t voxel = mesh_ray_voxel_index(axis, bin, a, b, nx, ny);
            double value = (hi - lo) / (nodes[bin + 1] - nodes[bin]);
            if (mesh_ray_accum_add(&accums[voxel], cell_id, material_id,
                                   value) != 0) return -1;
        }
        bin++;
    }
    return 0;
}

static int mesh_ray_trace_direction(alea_system_t *sys,
                                    const alea_mesh_config_t *cfg,
                                    const double *xn, const double *yn,
                                    const double *zn, int axis,
                                    mesh_ray_voxel_accum_t *accums) {
    const double *along_nodes = axis == 0 ? xn : axis == 1 ? yn : zn;
    const double *u_nodes = axis == 0 ? yn : xn;
    const double *v_nodes = axis == 2 ? yn : zn;
    const int along_count = axis == 0 ? cfg->nx : axis == 1 ? cfg->ny : cfg->nz;
    const int u_count = axis == 0 ? cfg->ny : cfg->nx;
    const int v_count = axis == 2 ? cfg->ny : cfg->nz;
    const size_t column_count = (size_t)u_count * (size_t)v_count;
#ifdef _OPENMP
    int workers = cfg->workers > 0 ? cfg->workers : omp_get_max_threads();
    if ((size_t)workers > column_count) workers = (int)column_count;
    if (workers < 1 || omp_in_parallel()) workers = 1;
#endif
    int failed = 0;
#ifdef _OPENMP
    #pragma omp parallel for schedule(static) num_threads(workers)
#endif
    for (size_t column = 0; column < column_count; column++) {
        int a = (int)(column % (size_t)u_count);
        int b = (int)(column / (size_t)u_count);
        alea_raycast_result_t trace;
        alea_raycast_result_init(&trace);
        for (int gv = 0; gv < cfg->ray_grid_v; gv++) {
            for (int gu = 0; gu < cfg->ray_grid_u; gu++) {
                const double u = u_nodes[a] +
                    ((double)gu + 0.5) / (double)cfg->ray_grid_u *
                    (u_nodes[a + 1] - u_nodes[a]);
                const double v = v_nodes[b] +
                    ((double)gv + 0.5) / (double)cfg->ray_grid_v *
                    (v_nodes[b + 1] - v_nodes[b]);
                double origin[3] = {0.0, 0.0, 0.0};
                double direction[3] = {0.0, 0.0, 0.0};
                origin[axis] = nextafter(along_nodes[0], -INFINITY);
                if (axis == 0) { origin[1] = u; origin[2] = v; }
                else if (axis == 1) { origin[0] = u; origin[2] = v; }
                else { origin[0] = u; origin[1] = v; }
                direction[axis] = 1.0;
                alea_ray_t ray;
                alea_raycast_result_clear(&trace);
                if (alea_ray_init(&ray, origin[0], origin[1], origin[2],
                                  direction[0], direction[1], direction[2]) != 0 ||
                    alea_raycast_hier_segments_nocache(
                        sys, &ray, along_nodes[along_count] - origin[axis],
                        &trace) != 0) {
#ifdef _OPENMP
                    #pragma omp atomic write
#endif
                    failed = 1;
                    continue;
                }
                for (size_t s = 0; s < trace.segments.count; s++) {
                    const alea_ray_segment_t *segment = &trace.segments.data[s];
                    if (mesh_ray_accumulate_segment(
                            accums, along_nodes, along_count, axis, a, b,
                            cfg->nx, cfg->ny,
                            origin[axis] + segment->t_enter,
                            origin[axis] + segment->t_exit,
                            segment->cell_id, segment->material_id) != 0) {
#ifdef _OPENMP
                        #pragma omp atomic write
#endif
                        failed = 1;
                        break;
                    }
                }
                if (alea_interrupted()) {
#ifdef _OPENMP
                    #pragma omp atomic write
#endif
                    failed = 1;
                }
            }
        }
        alea_raycast_result_free(&trace);
    }
    if (failed) {
        if (alea_interrupted())
            alea_set_error_detail(ALEA_ERR_INTERRUPTED,
                                  "mesh ray sampling interrupted");
        else
            alea_set_error_detail(ALEA_ERR_INVALID_STATE,
                                  "mesh ray sampling failed");
        return -1;
    }
    return 0;
}

static alea_mesh_result_t *mesh_sample_rays(
    alea_system_t *sys, const alea_mesh_config_t *cfg,
    double *xn, double *yn, double *zn, size_t ncells,
    alea_mesh_bounds_source_t bounds_source, double bounds_padding) {
    alea_mesh_result_t *res = NULL;
    mesh_ray_voxel_accum_t *accums = calloc(ncells, sizeof(*accums));
    int *mat_ids = NULL, *cell_ids = NULL, *unique_mats = NULL;
    unsigned char *mixed_flags = NULL;
    double *dominant_fractions = NULL, *estimated_errors = NULL;
    uint32_t *sample_counts = NULL;
    uint8_t *tie_flags = NULL, *refinement_flags = NULL;
    alea_mesh_fraction_span_t *fraction_spans = NULL, *cell_fraction_spans = NULL;
    alea_mesh_material_fraction_t *fractions = NULL;
    alea_mesh_cell_fraction_t *cell_fractions = NULL;
    size_t fraction_count = 0, fraction_capacity = 0;
    size_t cell_fraction_count = 0, cell_fraction_capacity = 0;
    int *scratch_mats = NULL;
    double *scratch_values = NULL;
    mesh_int_set_t material_set = {0};
    int mixed_count = 0, num_mats = 0;
    const uint32_t directions = cfg->ray_directions;
    const uint32_t direction_count =
        ((directions & ALEA_MESH_RAY_X) ? 1u : 0u) +
        ((directions & ALEA_MESH_RAY_Y) ? 1u : 0u) +
        ((directions & ALEA_MESH_RAY_Z) ? 1u : 0u);
    const uint64_t rays_per_direction =
        (uint64_t)cfg->ray_grid_u * (uint64_t)cfg->ray_grid_v;
    const uint64_t denominator64 = rays_per_direction * direction_count;
    if (!accums || denominator64 == 0 || denominator64 > UINT32_MAX) {
        alea_set_error_detail(ALEA_ERR_OVERFLOW,
                              "mesh ray samples per voxel overflow uint32_t");
        goto fail;
    }
    if (cfg->visit) {
        alea_set_error_detail(ALEA_ERR_UNSUPPORTED,
                              "streaming visitors are not yet supported for ray meshes");
        goto fail;
    }
    if (alea_system_prepare_query_caches(
            sys, ALEA_CACHE_HIER_SPATIAL | ALEA_CACHE_CELL_SURFACES) != 0)
        goto fail;
    if ((directions & ALEA_MESH_RAY_X) &&
        mesh_ray_trace_direction(sys, cfg, xn, yn, zn, 0, accums) != 0) goto fail;
    if ((directions & ALEA_MESH_RAY_Y) &&
        mesh_ray_trace_direction(sys, cfg, xn, yn, zn, 1, accums) != 0) goto fail;
    if ((directions & ALEA_MESH_RAY_Z) &&
        mesh_ray_trace_direction(sys, cfg, xn, yn, zn, 2, accums) != 0) goto fail;

    size_t max_components = 1;
    for (size_t v = 0; v < ncells; v++) {
        double sum = 0.0;
        for (uint32_t i = 0; i < accums[v].count; i++)
            sum += accums[v].data[i].value;
        double residual = (double)denominator64 - sum;
        if (residual > 1e-10 && mesh_ray_accum_add(
                &accums[v], -1, cfg->void_material_id, residual) != 0) {
            alea_set_error_detail(ALEA_ERR_OUT_OF_MEMORY,
                                  "mesh ray composition allocation failed");
            goto fail;
        }
        if (sum > (double)denominator64 + 1e-8) {
            alea_set_error_detail(ALEA_ERR_INVALID_STATE,
                                  "mesh ray ownership %.17g exceeds voxel measure %llu at voxel %zu",
                                  sum, (unsigned long long)denominator64, v);
            goto fail;
        }
        if (accums[v].count > max_components) max_components = accums[v].count;
    }
    scratch_mats = malloc(max_components * sizeof(*scratch_mats));
    scratch_values = malloc(max_components * sizeof(*scratch_values));
    if (!scratch_mats || !scratch_values) goto fail;

#define RAY_ALLOC(FIELD, TARGET, TYPE, CLEAR) do { \
    if (cfg->fields & (FIELD)) { \
        TARGET = mesh_alloc_array(ncells, sizeof(TYPE), (CLEAR)); \
        if (!(TARGET)) goto fail; \
    } \
} while (0)
    RAY_ALLOC(ALEA_MESH_FIELD_MATERIAL_ID, mat_ids, int, 0);
    RAY_ALLOC(ALEA_MESH_FIELD_CELL_ID, cell_ids, int, 0);
    RAY_ALLOC(ALEA_MESH_FIELD_MIXED_FLAG, mixed_flags, unsigned char, 1);
    RAY_ALLOC(ALEA_MESH_FIELD_DOMINANT_FRACTION, dominant_fractions, double, 0);
    RAY_ALLOC(ALEA_MESH_FIELD_ESTIMATED_ERROR, estimated_errors, double, 1);
    RAY_ALLOC(ALEA_MESH_FIELD_SAMPLE_COUNT, sample_counts, uint32_t, 0);
    RAY_ALLOC(ALEA_MESH_FIELD_TIE_FLAG, tie_flags, uint8_t, 1);
    RAY_ALLOC(ALEA_MESH_FIELD_REFINEMENT_FLAG, refinement_flags, uint8_t, 1);
    RAY_ALLOC(ALEA_MESH_FIELD_SAMPLED_FRACTIONS, fraction_spans,
              alea_mesh_fraction_span_t, 1);
    RAY_ALLOC(ALEA_MESH_FIELD_CELL_FRACTIONS, cell_fraction_spans,
              alea_mesh_fraction_span_t, 1);
#undef RAY_ALLOC

    for (size_t v = 0; v < ncells; v++) {
        mesh_ray_voxel_accum_t *accum = &accums[v];
        qsort(accum->data, accum->count, sizeof(*accum->data),
              mesh_ray_component_compare);
        int material_count = 0;
        for (uint32_t i = 0; i < accum->count; i++) {
            int m = 0;
            while (m < material_count && scratch_mats[m] != accum->data[i].material_id)
                m++;
            if (m == material_count) {
                scratch_mats[m] = accum->data[i].material_id;
                scratch_values[m] = 0.0;
                material_count++;
            }
            scratch_values[m] += accum->data[i].value;
        }
        for (int i = 1; i < material_count; i++) {
            int id = scratch_mats[i]; double value = scratch_values[i]; int j = i;
            while (j > 0 && scratch_mats[j - 1] > id) {
                scratch_mats[j] = scratch_mats[j - 1];
                scratch_values[j] = scratch_values[j - 1]; j--;
            }
            scratch_mats[j] = id; scratch_values[j] = value;
        }
        int dominant_material = cfg->void_material_id;
        double dominant_value = -1.0;
        int material_ties = 0;
        for (int m = 0; m < material_count; m++) {
            if (scratch_values[m] > dominant_value + 1e-12) {
                dominant_value = scratch_values[m];
                dominant_material = scratch_mats[m];
                material_ties = 1;
            } else if (fabs(scratch_values[m] - dominant_value) <= 1e-12) {
                material_ties++;
                if (scratch_mats[m] < dominant_material)
                    dominant_material = scratch_mats[m];
            }
            if (mesh_int_set_insert(&material_set, scratch_mats[m]) != 0) goto fail;
        }
        int dominant_cell = -1, cell_ties = 0;
        double dominant_cell_value = -1.0;
        for (uint32_t i = 0; i < accum->count; i++) {
            if (accum->data[i].material_id != dominant_material) continue;
            if (accum->data[i].value > dominant_cell_value + 1e-12) {
                dominant_cell_value = accum->data[i].value;
                dominant_cell = accum->data[i].cell_id;
                cell_ties = 1;
            } else if (fabs(accum->data[i].value - dominant_cell_value) <= 1e-12) {
                cell_ties++;
                if (accum->data[i].cell_id < dominant_cell)
                    dominant_cell = accum->data[i].cell_id;
            }
        }
        double dominant_fraction = dominant_value / (double)denominator64;
        unsigned char mixed = dominant_fraction < 1.0 - cfg->mixed_threshold;
        uint8_t ties = (material_ties > 1 ? ALEA_MESH_TIE_MATERIAL : 0) |
                       (cell_ties > 1 ? ALEA_MESH_TIE_CELL : 0);
        if (mat_ids) mat_ids[v] = dominant_material;
        if (cell_ids) cell_ids[v] = dominant_cell;
        if (mixed_flags) mixed_flags[v] = mixed;
        if (dominant_fractions) dominant_fractions[v] = dominant_fraction;
        if (sample_counts) sample_counts[v] = (uint32_t)denominator64;
        if (tie_flags) tie_flags[v] = ties;
        mixed_count += mixed ? 1 : 0;

        if (fraction_spans) {
            if (fraction_count > UINT32_MAX ||
                (size_t)material_count > UINT32_MAX - fraction_count ||
                ensure_fraction_capacity(&fractions, &fraction_capacity,
                                         fraction_count + (size_t)material_count) != 0)
                goto fail;
            fraction_spans[v] = (alea_mesh_fraction_span_t){
                (uint32_t)fraction_count, (uint32_t)material_count};
            for (int m = 0; m < material_count; m++) {
                fractions[fraction_count++] = (alea_mesh_material_fraction_t){
                    scratch_mats[m], scratch_values[m] / (double)denominator64};
            }
        }
        if (cell_fraction_spans) {
            if (cell_fraction_count > UINT32_MAX ||
                (size_t)accum->count > UINT32_MAX - cell_fraction_count ||
                ensure_cell_fraction_capacity(
                    &cell_fractions, &cell_fraction_capacity,
                    cell_fraction_count + accum->count) != 0) goto fail;
            cell_fraction_spans[v] = (alea_mesh_fraction_span_t){
                (uint32_t)cell_fraction_count, accum->count};
            for (uint32_t i = 0; i < accum->count; i++) {
                cell_fractions[cell_fraction_count++] = (alea_mesh_cell_fraction_t){
                    accum->data[i].cell_id, accum->data[i].material_id,
                    accum->data[i].value / (double)denominator64};
            }
        }
    }
    unique_mats = mesh_int_set_sorted_array(&material_set, &num_mats);
    if (!unique_mats) goto fail;
    res = calloc(1, sizeof(*res));
    if (!res) goto fail;
    res->nx = cfg->nx; res->ny = cfg->ny; res->nz = cfg->nz;
    res->fields = cfg->fields; res->bounds_source = bounds_source;
    res->bounds_padding = bounds_padding; res->sampling_mode = cfg->sampling_mode;
    res->sampling_seed = cfg->sampling_seed; res->target_error = cfg->target_error;
    res->x_nodes = xn; res->y_nodes = yn; res->z_nodes = zn;
    res->material_ids = mat_ids; res->cell_ids = cell_ids;
    res->num_materials = num_mats; res->unique_materials = unique_mats;
    res->mixed_flags = mixed_flags; res->dominant_fractions = dominant_fractions;
    res->estimated_errors = estimated_errors; res->sample_counts = sample_counts;
    res->tie_flags = tie_flags; res->refinement_flags = refinement_flags;
    res->mixed_count = mixed_count; res->fraction_spans = fraction_spans;
    res->fractions = fractions; res->fraction_count = fraction_count;
    res->cell_fraction_spans = cell_fraction_spans;
    res->cell_fractions = cell_fractions;
    res->cell_fraction_count = cell_fraction_count;
    for (size_t v = 0; v < ncells; v++) free(accums[v].data);
    free(accums); free(scratch_mats); free(scratch_values);
    mesh_int_set_free(&material_set);
    return res;

fail:
    if (accums) for (size_t v = 0; v < ncells; v++) free(accums[v].data);
    free(accums); free(scratch_mats); free(scratch_values);
    mesh_int_set_free(&material_set);
    free(xn); free(yn); free(zn); free(mat_ids); free(cell_ids); free(unique_mats);
    free(mixed_flags); free(dominant_fractions); free(estimated_errors);
    free(sample_counts); free(tie_flags); free(refinement_flags);
    free(fraction_spans); free(fractions);
    free(cell_fraction_spans); free(cell_fractions); free(res);
    return NULL;
}

/* ============================================================================
 * Sampling
 * ============================================================================ */

static int mesh_root_world_aabb(const alea_system_t *sys, alea_bbox_t *out);

alea_mesh_result_t *alea_mesh_sample(alea_system_t *sys,
                                             const alea_mesh_config_t *cfg) {
    if (!sys) {
        alea_set_error_detail(ALEA_ERR_NULL_ARG, "mesh system is NULL");
        return NULL;
    }
    if (validate_config_basic(cfg) != 0) return NULL;

    int nx = cfg->nx, ny = cfg->ny, nz = cfg->nz;

    size_t ncells = 0;
    if (checked_mesh_cell_count(nx, ny, nz, &ncells) != 0) return NULL;
    if (ncells > (size_t)INT_MAX) {
        alea_set_error_detail(ALEA_ERR_OVERFLOW,
                              "mesh has too many cells for legacy exporters");
        return NULL;
    }

    size_t nnodes = 0;
    if (checked_mesh_node_count(nx, ny, nz, &nnodes) != 0) return NULL;
    if (nnodes > (size_t)INT_MAX) {
        alea_set_error_detail(ALEA_ERR_OVERFLOW,
                              "mesh has too many nodes for legacy exporters");
        return NULL;
    }

    /* Determine bounds */
    double xlo = cfg->x_min, xhi = cfg->x_max;
    double ylo = cfg->y_min, yhi = cfg->y_max;
    double zlo = cfg->z_min, zhi = cfg->z_max;

    int all_bounds_zero =
        xlo == 0.0 && xhi == 0.0 && ylo == 0.0 && yhi == 0.0 &&
        zlo == 0.0 && zhi == 0.0;
    int needs_uniform_axis = !cfg->x_nodes || !cfg->y_nodes || !cfg->z_nodes;
    int auto_bounds = cfg->bounds_mode == ALEA_MESH_BOUNDS_AUTO ||
        (cfg->bounds_mode == ALEA_MESH_BOUNDS_LEGACY && all_bounds_zero);
    alea_mesh_bounds_source_t bounds_source =
        (!needs_uniform_axis) ? ALEA_MESH_BOUNDS_SOURCE_CUSTOM_NODES :
        (auto_bounds ? ALEA_MESH_BOUNDS_SOURCE_INFERRED_ROOT_AABB :
                       ALEA_MESH_BOUNDS_SOURCE_EXPLICIT);

    if (auto_bounds && needs_uniform_axis) {
        alea_bbox_t root_bounds;
        if (mesh_root_world_aabb(sys, &root_bounds) != 0) return NULL;
        const double xpad = (root_bounds.max_x - root_bounds.min_x) * cfg->auto_pad;
        const double ypad = (root_bounds.max_y - root_bounds.min_y) * cfg->auto_pad;
        const double zpad = (root_bounds.max_z - root_bounds.min_z) * cfg->auto_pad;
        if (!cfg->x_nodes) { xlo = root_bounds.min_x - xpad; xhi = root_bounds.max_x + xpad; }
        if (!cfg->y_nodes) { ylo = root_bounds.min_y - ypad; yhi = root_bounds.max_y + ypad; }
        if (!cfg->z_nodes) { zlo = root_bounds.min_z - zpad; zhi = root_bounds.max_z + zpad; }
    }

    if (!cfg->x_nodes && validate_uniform_axis(xlo, xhi, "X") != 0) return NULL;
    if (!cfg->y_nodes && validate_uniform_axis(ylo, yhi, "Y") != 0) return NULL;
    if (!cfg->z_nodes && validate_uniform_axis(zlo, zhi, "Z") != 0) return NULL;

    /* Build node arrays */
    double *xn = cfg->x_nodes ? copy_nodes(cfg->x_nodes, nx) : build_uniform_nodes(nx, xlo, xhi);
    double *yn = cfg->y_nodes ? copy_nodes(cfg->y_nodes, ny) : build_uniform_nodes(ny, ylo, yhi);
    double *zn = cfg->z_nodes ? copy_nodes(cfg->z_nodes, nz) : build_uniform_nodes(nz, zlo, zhi);
    if (!xn || !yn || !zn) {
        alea_set_error_detail(ALEA_ERR_OUT_OF_MEMORY,
                              "mesh failed to allocate node arrays");
        free(xn); free(yn); free(zn);
        return NULL;
    }

    if (cfg->sampling_mode == ALEA_MESH_SAMPLE_RAY)
        return mesh_sample_rays(sys, cfg, xn, yn, zn, ncells, bounds_source,
                                auto_bounds ? cfg->auto_pad : 0.0);

    int *mat_ids = (cfg->fields & ALEA_MESH_FIELD_MATERIAL_ID) ?
        mesh_alloc_array(ncells, sizeof(int), 0) : NULL;
    int *cell_ids = (cfg->fields & ALEA_MESH_FIELD_CELL_ID) ?
        mesh_alloc_array(ncells, sizeof(int), 0) : NULL;
    if (((cfg->fields & ALEA_MESH_FIELD_MATERIAL_ID) && !mat_ids) ||
        ((cfg->fields & ALEA_MESH_FIELD_CELL_ID) && !cell_ids)) {
        alea_set_error_detail(ALEA_ERR_OUT_OF_MEMORY,
                              "mesh failed to allocate cell data arrays");
        free(xn); free(yn); free(zn); free(mat_ids); free(cell_ids);
        return NULL;
    }

    unsigned char *mixed_flags = NULL;
    double *dominant_fractions = NULL;
    double *estimated_errors = NULL;
    uint32_t *sample_counts = NULL;
    uint8_t *tie_flags = NULL;
    uint8_t *refinement_flags = NULL;
    alea_mesh_fraction_span_t *fraction_spans = NULL;
    alea_mesh_material_fraction_t *fractions = NULL;
    size_t fraction_count = 0;
    size_t fraction_capacity = 0;
    alea_mesh_fraction_span_t *cell_fraction_spans = NULL;
    alea_mesh_cell_fraction_t *cell_fractions = NULL;
    size_t cell_fraction_count = 0;
    size_t cell_fraction_capacity = 0;
    if (cfg->fields & ALEA_MESH_FIELD_MIXED_FLAG)
        mixed_flags = mesh_alloc_array(ncells, sizeof(unsigned char), 1);
    if (cfg->fields & ALEA_MESH_FIELD_DOMINANT_FRACTION)
        dominant_fractions = mesh_alloc_array(ncells, sizeof(double), 0);
    if (cfg->fields & ALEA_MESH_FIELD_ESTIMATED_ERROR)
        estimated_errors = mesh_alloc_array(ncells, sizeof(double), 1);
    if (cfg->fields & ALEA_MESH_FIELD_SAMPLE_COUNT)
        sample_counts = mesh_alloc_array(ncells, sizeof(uint32_t), 0);
    if (cfg->fields & ALEA_MESH_FIELD_TIE_FLAG)
        tie_flags = mesh_alloc_array(ncells, sizeof(uint8_t), 1);
    if (cfg->fields & ALEA_MESH_FIELD_REFINEMENT_FLAG)
        refinement_flags = mesh_alloc_array(ncells, sizeof(uint8_t), 1);
    if (cfg->fields & ALEA_MESH_FIELD_SAMPLED_FRACTIONS)
        fraction_spans = mesh_alloc_array(ncells,
                                          sizeof(alea_mesh_fraction_span_t), 1);
    if (cfg->fields & ALEA_MESH_FIELD_CELL_FRACTIONS)
        cell_fraction_spans = mesh_alloc_array(
            ncells, sizeof(alea_mesh_fraction_span_t), 1);
    if (((cfg->fields & ALEA_MESH_FIELD_MIXED_FLAG) && !mixed_flags) ||
        ((cfg->fields & ALEA_MESH_FIELD_DOMINANT_FRACTION) && !dominant_fractions) ||
        ((cfg->fields & ALEA_MESH_FIELD_ESTIMATED_ERROR) && !estimated_errors) ||
        ((cfg->fields & ALEA_MESH_FIELD_SAMPLE_COUNT) && !sample_counts) ||
        ((cfg->fields & ALEA_MESH_FIELD_TIE_FLAG) && !tie_flags) ||
        ((cfg->fields & ALEA_MESH_FIELD_REFINEMENT_FLAG) && !refinement_flags) ||
        ((cfg->fields & ALEA_MESH_FIELD_SAMPLED_FRACTIONS) && !fraction_spans) ||
        ((cfg->fields & ALEA_MESH_FIELD_CELL_FRACTIONS) && !cell_fraction_spans)) {
        alea_set_error_detail(ALEA_ERR_OUT_OF_MEMORY,
                              "mesh failed to allocate composition arrays");
        free(xn); free(yn); free(zn);
        free(mat_ids); free(cell_ids);
        free(mixed_flags); free(dominant_fractions); free(estimated_errors);
        free(sample_counts); free(tie_flags); free(refinement_flags);
        free(fraction_spans); free(cell_fraction_spans);
        return NULL;
    }

    /* Ensure universe lookup is ready before point sampling. Mesh sampling
     * uses recursive point lookup, not the expanded flat spatial index. */
    if (alea_system_prepare_query_caches(sys, ALEA_CACHE_UNIVERSE) != 0) {
        free(xn); free(yn); free(zn);
        free(mat_ids); free(cell_ids);
        free(mixed_flags); free(dominant_fractions); free(estimated_errors);
        free(sample_counts); free(tie_flags); free(refinement_flags);
        free(fraction_spans); free(cell_fraction_spans);
        return NULL;
    }

    int sample_axis = cfg->sampling_mode == ALEA_MESH_SAMPLE_CENTER ? 1 :
                      cfg->sampling_mode == ALEA_MESH_SAMPLE_CORNERS ? 2 :
                      cfg->subsamples_per_axis;
    uint32_t base_samples = 0;
    if (mesh_cube_samples(sample_axis, &base_samples) != 0 ||
        base_samples > cfg->max_samples_per_voxel) {
        alea_set_error_detail(ALEA_ERR_INVALID_ARG,
                              "mesh base sampling exceeds max_samples_per_voxel");
        free(xn); free(yn); free(zn); free(mat_ids); free(cell_ids);
        free(mixed_flags); free(dominant_fractions); free(estimated_errors);
        free(sample_counts); free(tie_flags); free(refinement_flags);
        free(fraction_spans); free(cell_fraction_spans);
        return NULL;
    }
    if (cfg->max_total_samples != 0 &&
        ((uint64_t)ncells > UINT64_MAX / base_samples ||
         (uint64_t)ncells * base_samples > cfg->max_total_samples)) {
        alea_set_error_detail(ALEA_ERR_INVALID_ARG,
                              "mesh max_total_samples cannot cover the base grid");
        free(xn); free(yn); free(zn); free(mat_ids); free(cell_ids);
        free(mixed_flags); free(dominant_fractions); free(estimated_errors);
        free(sample_counts); free(tie_flags); free(refinement_flags);
        free(fraction_spans); free(cell_fraction_spans);
        return NULL;
    }
    int max_sample_axis = sample_axis;
    if (cfg->sampling_mode == ALEA_MESH_SAMPLE_ADAPTIVE) {
        uint64_t cumulative = base_samples;
        for (int depth = 0; depth < cfg->max_refine_depth; depth++) {
            if (max_sample_axis > INT_MAX / 2) break;
            uint32_t next_samples = 0;
            int next_axis = max_sample_axis * 2;
            if (mesh_cube_samples(next_axis, &next_samples) != 0 ||
                cumulative + next_samples > cfg->max_samples_per_voxel) break;
            cumulative += next_samples;
            max_sample_axis = next_axis;
        }
    }
    size_t scratch_count = (size_t)max_sample_axis * (size_t)max_sample_axis *
                           (size_t)max_sample_axis;
    int worker_count = 1;
    int parallel_sampling = 0;
#ifdef _OPENMP
    if (cfg->workers != 1 && !fraction_spans && !cell_fraction_spans &&
        !cfg->visit && !cfg->progress &&
        cfg->sampling_mode != ALEA_MESH_SAMPLE_ADAPTIVE &&
        cfg->max_total_samples == 0) {
        worker_count = cfg->workers > 0 ? cfg->workers : omp_get_max_threads();
        if ((size_t)worker_count > ncells) worker_count = (int)ncells;
        if (worker_count > 1 && !omp_in_parallel()) parallel_sampling = 1;
        else worker_count = 1;
    }
#endif
    size_t scratch_bytes = 0;
    size_t worker_scratch_count = 0;
    if (checked_mul_size(scratch_count, (size_t)worker_count,
                         &worker_scratch_count) != 0 ||
        checked_mul_size(worker_scratch_count, 7 * sizeof(int),
                         &scratch_bytes) != 0) {
        free(xn); free(yn); free(zn); free(mat_ids); free(cell_ids);
        free(mixed_flags); free(dominant_fractions); free(estimated_errors);
        free(sample_counts); free(tie_flags); free(refinement_flags);
        free(fraction_spans); free(cell_fraction_spans);
        return NULL;
    }
    int *scratch = malloc(scratch_bytes);
    if (!scratch) {
        alea_set_error_detail(ALEA_ERR_OUT_OF_MEMORY,
                              "mesh failed to allocate sampling scratch");
        free(xn); free(yn); free(zn); free(mat_ids); free(cell_ids);
        free(mixed_flags); free(dominant_fractions); free(estimated_errors);
        free(sample_counts); free(tie_flags); free(refinement_flags);
        free(fraction_spans); free(cell_fraction_spans);
        return NULL;
    }
    int *scratch_materials = scratch;
    int *scratch_counts = scratch + worker_scratch_count;
    int *scratch_owner_materials = scratch + 2 * worker_scratch_count;
    int *scratch_owner_cells = scratch + 3 * worker_scratch_count;
    int *scratch_owner_counts = scratch + 4 * worker_scratch_count;
    int *scratch_previous_materials = scratch + 5 * worker_scratch_count;
    int *scratch_previous_counts = scratch + 6 * worker_scratch_count;
    alea_mesh_material_fraction_t *voxel_fractions = cfg->visit ?
        mesh_alloc_array(scratch_count, sizeof(*voxel_fractions), 0) : NULL;
    alea_mesh_cell_fraction_t *voxel_cell_fractions = cfg->visit ?
        mesh_alloc_array(scratch_count, sizeof(*voxel_cell_fractions), 0) : NULL;
    if (cfg->visit && (!voxel_fractions || !voxel_cell_fractions)) {
        free(scratch);
        free(xn); free(yn); free(zn); free(mat_ids); free(cell_ids);
        free(mixed_flags); free(dominant_fractions); free(estimated_errors);
        free(sample_counts); free(tie_flags); free(refinement_flags);
        free(fraction_spans); free(cell_fraction_spans);
        free(voxel_fractions); free(voxel_cell_fractions);
        return NULL;
    }

    /* Sample voxel composition and coherent dominant cell/material IDs. */
    size_t nxy = (size_t)nx * (size_t)ny;
    int mixed_count = 0;
    uint64_t total_sample_work = 0;
    mesh_int_set_t material_set = {0};

    if (parallel_sampling) {
#ifdef _OPENMP
        int parallel_error = 0;
        #pragma omp parallel for schedule(static) num_threads(worker_count) reduction(+:mixed_count)
        for (size_t idx = 0; idx < ncells; idx++) {
            int tid = omp_get_thread_num();
            size_t offset = (size_t)tid * scratch_count;
            int *materials = scratch_materials + offset;
            int *counts = scratch_counts + offset;
            int *owner_materials = scratch_owner_materials + offset;
            int *owner_cells = scratch_owner_cells + offset;
            int *owner_counts = scratch_owner_counts + offset;
            int i = (int)(idx % (size_t)nx);
            int j = (int)((idx / (size_t)nx) % (size_t)ny);
            int k = (int)(idx / nxy);
            int n_materials = 0, n_owners = 0, nsamples = 0;
            int dominant_material = cfg->void_material_id, dominant_cell = -1;
            unsigned char mixed = 0;
            uint8_t ties = 0;
            double dominant = 1.0;

            mesh_sample_voxel_materials(sys, cfg, xn, yn, zn, i, j, k,
                                        sample_axis, materials, counts,
                                        &n_materials, owner_materials,
                                        owner_cells, owner_counts,
                                        &n_owners, &nsamples);
            mesh_fraction_stats(materials, counts, n_materials, nsamples, cfg,
                                owner_materials, owner_cells, owner_counts,
                                n_owners, &dominant_material, &dominant_cell,
                                &mixed, &dominant, &ties);
            if (mat_ids) mat_ids[idx] = dominant_material;
            if (cell_ids) cell_ids[idx] = dominant_cell;
            if (mixed_flags) mixed_flags[idx] = mixed;
            if (dominant_fractions) dominant_fractions[idx] = dominant;
            if (sample_counts) sample_counts[idx] = (uint32_t)nsamples;
            if (tie_flags) tie_flags[idx] = ties;
            mixed_count += mixed ? 1 : 0;
            if (alea_interrupted()) {
                #pragma omp atomic write
                parallel_error = 1;
            }
            for (int m = 0; m < n_materials; m++) {
                int rc;
                #pragma omp critical(alea_mesh_material_set)
                rc = mesh_int_set_insert(&material_set, materials[m]);
                if (rc != 0) {
                    #pragma omp atomic write
                    parallel_error = 1;
                }
            }
        }
        if (parallel_error) {
            if (alea_interrupted())
                alea_set_error_detail(ALEA_ERR_INTERRUPTED,
                                      "parallel mesh sampling interrupted");
            free(xn); free(yn); free(zn); free(mat_ids); free(cell_ids);
            free(mixed_flags); free(dominant_fractions); free(estimated_errors);
            free(sample_counts); free(tie_flags); free(refinement_flags);
            free(fraction_spans); free(fractions);
            free(cell_fraction_spans); free(cell_fractions); free(scratch);
            free(voxel_fractions); free(voxel_cell_fractions);
            mesh_int_set_free(&material_set);
            return NULL;
        }
#endif
    } else for (int k = 0; k < nz; k++) {
        for (int j = 0; j < ny; j++) {
            for (int i = 0; i < nx; i++) {
                size_t idx = (size_t)k * nxy + (size_t)j * (size_t)nx + (size_t)i;
                int *materials = scratch_materials;
                int *counts = scratch_counts;
                int *owner_materials = scratch_owner_materials;
                int *owner_cells = scratch_owner_cells;
                int *owner_counts = scratch_owner_counts;
                int n_materials = 0;
                int n_owners = 0;
                int nsamples = 0;
                int dominant_material = cfg->void_material_id;
                int dominant_cell = -1;
                unsigned char mixed = 0;
                uint8_t ties = 0;
                uint8_t refinement = 0;
                double dominant = 1.0;
                double estimated_error = 0.0;
                uint32_t voxel_sample_work = base_samples;

                mesh_sample_voxel_materials(sys, cfg, xn, yn, zn, i, j, k,
                                            sample_axis,
                                            materials, counts,
                                            &n_materials, owner_materials,
                                            owner_cells, owner_counts,
                                            &n_owners, &nsamples);
                total_sample_work += (uint32_t)nsamples;

                if (cfg->sampling_mode == ALEA_MESH_SAMPLE_ADAPTIVE) {
                    estimated_error = 1.0;
                    int axis = sample_axis;
                    for (int depth = 0; depth < cfg->max_refine_depth; depth++) {
                        memcpy(scratch_previous_materials, materials,
                               (size_t)n_materials * sizeof(int));
                        memcpy(scratch_previous_counts, counts,
                               (size_t)n_materials * sizeof(int));
                        int previous_n_materials = n_materials;
                        int previous_samples = nsamples;
                        uint32_t next_samples = 0;
                        if (axis > INT_MAX / 2 ||
                            mesh_cube_samples(axis * 2, &next_samples) != 0 ||
                            voxel_sample_work + (uint64_t)next_samples >
                                cfg->max_samples_per_voxel ||
                            (cfg->max_total_samples != 0 &&
                             total_sample_work + next_samples >
                                cfg->max_total_samples)) {
                            refinement |= ALEA_MESH_REFINEMENT_LIMIT_REACHED;
                            break;
                        }
                        axis *= 2;
                        mesh_sample_voxel_materials(sys, cfg, xn, yn, zn,
                                                    i, j, k, axis,
                                                    materials, counts,
                                                    &n_materials,
                                                    owner_materials,
                                                    owner_cells, owner_counts,
                                                    &n_owners, &nsamples);
                        voxel_sample_work += next_samples;
                        total_sample_work += next_samples;
                        estimated_error = mesh_distribution_error(
                            scratch_previous_materials, scratch_previous_counts,
                            previous_n_materials, previous_samples,
                            materials, counts, n_materials, nsamples);
                        if (estimated_error <= cfg->target_error) break;
                        if (depth + 1 == cfg->max_refine_depth)
                            refinement |= ALEA_MESH_REFINEMENT_LIMIT_REACHED;
                    }
                    if (cfg->max_refine_depth == 0)
                        refinement |= ALEA_MESH_REFINEMENT_LIMIT_REACHED;
                }
                mesh_fraction_stats(materials, counts, n_materials, nsamples, cfg,
                                    owner_materials, owner_cells, owner_counts,
                                    n_owners, &dominant_material, &dominant_cell,
                                    &mixed, &dominant, &ties);

                if (cfg->visit) {
                    for (int m = 0; m < n_materials; m++) {
                        voxel_fractions[m].material_id = materials[m];
                        voxel_fractions[m].fraction =
                            (double)counts[m] / (double)nsamples;
                    }
                    for (int o = 0; o < n_owners; o++) {
                        voxel_cell_fractions[o].cell_id = owner_cells[o];
                        voxel_cell_fractions[o].material_id = owner_materials[o];
                        voxel_cell_fractions[o].fraction =
                            (double)owner_counts[o] / (double)nsamples;
                    }
                    alea_mesh_voxel_sample_t sample = {
                        .i = i, .j = j, .k = k,
                        .x_min = xn[i], .x_max = xn[i + 1],
                        .y_min = yn[j], .y_max = yn[j + 1],
                        .z_min = zn[k], .z_max = zn[k + 1],
                        .material_id = dominant_material,
                        .cell_id = dominant_cell,
                        .mixed = mixed,
                        .tie_flags = ties,
                        .dominant_fraction = dominant,
                        .estimated_error = estimated_error,
                        .sample_count = voxel_sample_work,
                        .refinement_flags = refinement,
                        .fractions = voxel_fractions,
                        .fraction_count = (uint32_t)n_materials,
                        .cell_fractions = voxel_cell_fractions,
                        .cell_fraction_count = (uint32_t)n_owners
                    };
                    if (cfg->visit(&sample, cfg->visit_user_data) != 0) {
                        alea_set_error_detail(ALEA_ERR_INTERRUPTED,
                                              "mesh voxel visitor cancelled sampling");
                        free(xn); free(yn); free(zn); free(mat_ids); free(cell_ids);
                        free(mixed_flags); free(dominant_fractions); free(estimated_errors);
                        free(sample_counts); free(tie_flags); free(refinement_flags);
                        free(fraction_spans); free(fractions);
                        free(cell_fraction_spans); free(cell_fractions);
                        free(scratch); free(voxel_fractions);
                        free(voxel_cell_fractions);
                        mesh_int_set_free(&material_set);
                        return NULL;
                    }
                }
                if (mat_ids) mat_ids[idx] = dominant_material;
                if (cell_ids) cell_ids[idx] = dominant_cell;
                if (sample_counts) sample_counts[idx] = voxel_sample_work;
                if (tie_flags) tie_flags[idx] = ties;
                if (estimated_errors) estimated_errors[idx] = estimated_error;
                if (refinement_flags) refinement_flags[idx] = refinement;

                for (int m = 0; m < n_materials; m++) {
                    if (mesh_int_set_insert(&material_set, materials[m]) != 0) {
                        free(xn); free(yn); free(zn);
                        free(mat_ids); free(cell_ids);
                        free(mixed_flags); free(dominant_fractions); free(estimated_errors);
                        free(sample_counts); free(tie_flags); free(refinement_flags);
                        free(fraction_spans); free(fractions);
                        free(cell_fraction_spans); free(cell_fractions);
                        free(scratch); free(voxel_fractions);
                        free(voxel_cell_fractions);
                        mesh_int_set_free(&material_set);
                        return NULL;
                    }
                }

                if (fraction_spans &&
                    (fraction_count > (size_t)UINT32_MAX ||
                     (size_t)n_materials > (size_t)UINT32_MAX - fraction_count)) {
                    alea_set_error_detail(ALEA_ERR_OVERFLOW,
                                          "mesh fraction span offset exceeds uint32_t");
                    free(xn); free(yn); free(zn);
                    free(mat_ids); free(cell_ids);
                    free(mixed_flags); free(dominant_fractions); free(estimated_errors);
                    free(sample_counts); free(tie_flags); free(refinement_flags);
                    free(fraction_spans); free(fractions);
                    free(cell_fraction_spans); free(cell_fractions);
                    free(scratch); free(voxel_fractions); free(voxel_cell_fractions);
                    mesh_int_set_free(&material_set);
                    return NULL;
                }
                if (fraction_spans &&
                    ensure_fraction_capacity(&fractions, &fraction_capacity,
                                             fraction_count + (size_t)n_materials) != 0) {
                    free(xn); free(yn); free(zn);
                    free(mat_ids); free(cell_ids);
                    free(mixed_flags); free(dominant_fractions); free(estimated_errors);
                    free(sample_counts); free(tie_flags); free(refinement_flags);
                    free(fraction_spans); free(fractions);
                    free(cell_fraction_spans); free(cell_fractions);
                    free(scratch); free(voxel_fractions); free(voxel_cell_fractions);
                    mesh_int_set_free(&material_set);
                    return NULL;
                }

                if (fraction_spans) {
                    fraction_spans[idx].offset = (uint32_t)fraction_count;
                    fraction_spans[idx].count = (uint32_t)n_materials;
                    for (int m = 0; m < n_materials; m++) {
                        fractions[fraction_count].material_id = materials[m];
                        fractions[fraction_count].fraction =
                            nsamples > 0 ? (double)counts[m] / (double)nsamples : 0.0;
                        fraction_count++;
                    }
                }

                if (cell_fraction_spans &&
                    (cell_fraction_count > (size_t)UINT32_MAX ||
                     (size_t)n_owners > (size_t)UINT32_MAX - cell_fraction_count)) {
                    alea_set_error_detail(ALEA_ERR_OVERFLOW,
                                          "mesh cell-fraction span offset exceeds uint32_t");
                    free(xn); free(yn); free(zn); free(mat_ids); free(cell_ids);
                    free(mixed_flags); free(dominant_fractions); free(estimated_errors);
                    free(sample_counts); free(tie_flags); free(refinement_flags);
                    free(fraction_spans); free(fractions);
                    free(cell_fraction_spans); free(cell_fractions);
                    free(scratch); free(voxel_fractions); free(voxel_cell_fractions);
                    mesh_int_set_free(&material_set);
                    return NULL;
                }
                if (cell_fraction_spans &&
                    ensure_cell_fraction_capacity(
                        &cell_fractions, &cell_fraction_capacity,
                        cell_fraction_count + (size_t)n_owners) != 0) {
                    free(xn); free(yn); free(zn); free(mat_ids); free(cell_ids);
                    free(mixed_flags); free(dominant_fractions); free(estimated_errors);
                    free(sample_counts); free(tie_flags); free(refinement_flags);
                    free(fraction_spans); free(fractions);
                    free(cell_fraction_spans); free(cell_fractions);
                    free(scratch); free(voxel_fractions); free(voxel_cell_fractions);
                    mesh_int_set_free(&material_set);
                    return NULL;
                }
                if (cell_fraction_spans) {
                    cell_fraction_spans[idx].offset = (uint32_t)cell_fraction_count;
                    cell_fraction_spans[idx].count = (uint32_t)n_owners;
                    for (int o = 0; o < n_owners; o++) {
                        cell_fractions[cell_fraction_count].cell_id = owner_cells[o];
                        cell_fractions[cell_fraction_count].material_id = owner_materials[o];
                        cell_fractions[cell_fraction_count].fraction =
                            nsamples > 0 ? (double)owner_counts[o] / (double)nsamples : 0.0;
                        cell_fraction_count++;
                    }
                }

                if (mixed_flags) mixed_flags[idx] = mixed;
                if (dominant_fractions) dominant_fractions[idx] = dominant;
                mixed_count += mixed ? 1 : 0;
            }
        }
        int cancelled = alea_interrupted();
        if (!cancelled && cfg->progress) {
            const size_t completed = (size_t)(k + 1) * nxy;
            cancelled = cfg->progress(completed, ncells,
                                      cfg->progress_user_data) != 0;
        }
        if (cancelled) {
                alea_set_error_detail(ALEA_ERR_INTERRUPTED,
                                      "mesh sampling interrupted");
                free(xn); free(yn); free(zn);
                free(mat_ids); free(cell_ids);
                free(mixed_flags); free(dominant_fractions); free(estimated_errors);
                free(sample_counts); free(tie_flags); free(refinement_flags);
                free(fraction_spans); free(fractions);
                free(cell_fraction_spans); free(cell_fractions);
                free(scratch); free(voxel_fractions); free(voxel_cell_fractions);
                mesh_int_set_free(&material_set);
                return NULL;
        }
    }

    free(scratch); free(voxel_fractions); free(voxel_cell_fractions);

    /* Collect unique materials */
    int num_mats = 0;
    int *unique_mats = mesh_int_set_sorted_array(&material_set, &num_mats);
    mesh_int_set_free(&material_set);
    if (!unique_mats) {
        free(xn); free(yn); free(zn);
        free(mat_ids); free(cell_ids);
        free(mixed_flags); free(dominant_fractions); free(estimated_errors);
        free(sample_counts); free(tie_flags); free(refinement_flags);
        free(fraction_spans); free(fractions);
        free(cell_fraction_spans); free(cell_fractions);
        return NULL;
    }

    /* Build result */
    alea_mesh_result_t *res = calloc(1, sizeof(*res));
    if (!res) {
        free(xn); free(yn); free(zn); free(mat_ids); free(cell_ids); free(unique_mats);
        free(mixed_flags); free(dominant_fractions); free(estimated_errors);
        free(sample_counts); free(tie_flags); free(refinement_flags);
        free(fraction_spans); free(fractions);
        free(cell_fraction_spans); free(cell_fractions);
        return NULL;
    }
    res->nx = nx;  res->ny = ny;  res->nz = nz;
    res->fields = cfg->fields;
    res->bounds_source = bounds_source;
    res->bounds_padding = auto_bounds ? cfg->auto_pad : 0.0;
    res->sampling_mode = cfg->sampling_mode;
    res->sampling_seed = cfg->sampling_seed;
    res->target_error = cfg->target_error;
    res->x_nodes = xn;  res->y_nodes = yn;  res->z_nodes = zn;
    res->material_ids = mat_ids;
    res->cell_ids = cell_ids;
    res->num_materials = num_mats;
    res->unique_materials = unique_mats;
    res->mixed_flags = mixed_flags;
    res->dominant_fractions = dominant_fractions;
    res->estimated_errors = estimated_errors;
    res->sample_counts = sample_counts;
    res->tie_flags = tie_flags;
    res->refinement_flags = refinement_flags;
    res->mixed_count = mixed_count;
    res->fraction_spans = fraction_spans;
    res->fractions = fractions;
    res->fraction_count = fraction_count;
    res->cell_fraction_spans = cell_fraction_spans;
    res->cell_fractions = cell_fractions;
    res->cell_fraction_count = cell_fraction_count;
    return res;
}

int alea_mesh_visit(alea_system_t *sys, const alea_mesh_config_t *cfg,
                    alea_mesh_voxel_visit_fn visit, void *user_data) {
    if (!cfg || !visit) {
        alea_set_error_detail(ALEA_ERR_NULL_ARG,
                              "mesh visit requires config and callback");
        return -1;
    }
    alea_mesh_config_t streaming = *cfg;
    streaming.fields = 0;
    streaming.visit = visit;
    streaming.visit_user_data = user_data;
    alea_mesh_result_t *result = alea_mesh_sample(sys, &streaming);
    if (!result) return -1;
    alea_mesh_result_free(result);
    return 0;
}

void alea_mesh_result_free(alea_mesh_result_t *mesh) {
    if (!mesh) return;
    free(mesh->x_nodes);
    free(mesh->y_nodes);
    free(mesh->z_nodes);
    free(mesh->material_ids);
    free(mesh->cell_ids);
    free(mesh->unique_materials);
    free(mesh->mixed_flags);
    free(mesh->dominant_fractions);
    free(mesh->estimated_errors);
    free(mesh->sample_counts);
    free(mesh->tie_flags);
    free(mesh->refinement_flags);
    free(mesh->fraction_spans);
    free(mesh->fractions);
    free(mesh->cell_fraction_spans);
    free(mesh->cell_fractions);
    free(mesh);
}

/* ============================================================================
 * Gmsh .msh v2.2 ASCII Writer
 * ============================================================================ */

static int mesh_writef(FILE *out, const char *format, ...) {
    va_list args;
    va_start(args, format);
    int rc = vfprintf(out, format, args);
    va_end(args);
    if (rc < 0) {
        alea_set_error_detail(ALEA_ERR_FILE_WRITE,
                              "mesh stream write failed: %s", strerror(errno));
        return -1;
    }
    return 0;
}

static int mesh_root_world_aabb(const alea_system_t *sys, alea_bbox_t *out) {
    const double sentinel = 9.0e9; /* Internal unbounded-box representation. */
    int found = 0;
    *out = (alea_bbox_t){INFINITY, -INFINITY, INFINITY, -INFINITY,
                         INFINITY, -INFINITY};
    for (size_t i = 0; i < alea_vec_count(&sys->cells); i++) {
        const alea_cell_entry_t *cell = &sys->cells.data[i];
        if (cell->universe_id != 0 || cell->root_node_id == ALEA_NODE_ID_INVALID ||
            cell->root_node_id >= alea_vec_count(&sys->nodes)) continue;
        const alea_bbox_t box =
            alea_node_bbox_get(&sys->nodes.data[cell->root_node_id].bbox);
        if (!isfinite(box.min_x) || !isfinite(box.max_x) ||
            !isfinite(box.min_y) || !isfinite(box.max_y) ||
            !isfinite(box.min_z) || !isfinite(box.max_z) ||
            fabs(box.min_x) >= sentinel || fabs(box.max_x) >= sentinel ||
            fabs(box.min_y) >= sentinel || fabs(box.max_y) >= sentinel ||
            fabs(box.min_z) >= sentinel || fabs(box.max_z) >= sentinel ||
            !(box.min_x < box.max_x) || !(box.min_y < box.max_y) ||
            !(box.min_z < box.max_z)) continue;
        if (box.min_x < out->min_x) out->min_x = box.min_x;
        if (box.max_x > out->max_x) out->max_x = box.max_x;
        if (box.min_y < out->min_y) out->min_y = box.min_y;
        if (box.max_y > out->max_y) out->max_y = box.max_y;
        if (box.min_z < out->min_z) out->min_z = box.min_z;
        if (box.max_z > out->max_z) out->max_z = box.max_z;
        found = 1;
    }
    if (!found) {
        alea_set_error_detail(ALEA_ERR_INVALID_STATE,
                              "mesh auto-bounds found no bounded root-universe cell");
        return -1;
    }
    return 0;
}

#define MESH_WRITE(...) do { if (mesh_writef(out, __VA_ARGS__) != 0) return -1; } while (0)

/** Node ID: 1-based, (i,j,k) -> k*(ny+1)*(nx+1) + j*(nx+1) + i + 1 */
#define GMSH_NODE(i, j, k, nx1, ny1) \
    ((k) * (ny1) * (nx1) + (j) * (nx1) + (i) + 1)

static int gmsh_material_tag(const alea_mesh_result_t *mesh, int material_id) {
    for (int i = 0; i < mesh->num_materials; i++) {
        if (mesh->unique_materials[i] == material_id)
            return i + 1;
    }
    return 0;
}

static double mesh_sampled_fraction_at(const alea_mesh_result_t *mesh,
                                       size_t voxel, int material_id) {
    if (!mesh->fraction_spans || !mesh->fractions)
        return mesh->material_ids[voxel] == material_id ? 1.0 : 0.0;
    const alea_mesh_fraction_span_t span = mesh->fraction_spans[voxel];
    for (uint32_t i = 0; i < span.count; i++) {
        const alea_mesh_material_fraction_t *fraction =
            &mesh->fractions[(size_t)span.offset + i];
        if (fraction->material_id == material_id) return fraction->fraction;
    }
    return 0.0;
}

static int write_gmsh_element_data_begin(FILE *out, const char *name,
                                         int element_count) {
    MESH_WRITE("$ElementData\n");
    MESH_WRITE("1\n\"%s\"\n", name);
    MESH_WRITE("1\n0.0\n");
    MESH_WRITE("3\n0\n1\n%d\n", element_count);
    return 0;
}

static int write_gmsh(const alea_mesh_result_t *mesh, FILE *out,
                      const alea_mesh_export_options_t *options) {
    int nx = mesh->nx, ny = mesh->ny, nz = mesh->nz;
    int nx1 = nx + 1, ny1 = ny + 1, nz1 = nz + 1;
    size_t nnodes_sz = 0, nelems_sz = 0;
    if (checked_mesh_node_count(nx, ny, nz, &nnodes_sz) != 0) return -1;
    if (checked_mesh_cell_count(nx, ny, nz, &nelems_sz) != 0) return -1;
    if (nnodes_sz > (size_t)INT_MAX || nelems_sz > (size_t)INT_MAX) {
        alea_set_error_detail(ALEA_ERR_OVERFLOW,
                              "mesh is too large for Gmsh v2.2 export");
        return -1;
    }
    int nnodes = (int)nnodes_sz;
    int nelems = (int)nelems_sz;
    size_t nxy = (size_t)nx * (size_t)ny;

    /* Header */
    MESH_WRITE("$MeshFormat\n2.2 0 8\n$EndMeshFormat\n");

    /* Physical names use positive sequential tags; names keep original IDs. */
    MESH_WRITE("$PhysicalNames\n%d\n", mesh->num_materials);
    for (int m = 0; m < mesh->num_materials; m++) {
        MESH_WRITE("3 %d \"material_%d\"\n",
                m + 1, mesh->unique_materials[m]);
    }
    MESH_WRITE("$EndPhysicalNames\n");

    /* Nodes */
    MESH_WRITE("$Nodes\n%d\n", nnodes);
    for (int k = 0; k < nz1; k++) {
        for (int j = 0; j < ny1; j++) {
            for (int i = 0; i < nx1; i++) {
                int nid = GMSH_NODE(i, j, k, nx1, ny1);
                MESH_WRITE("%d %.15g %.15g %.15g\n",
                        nid, mesh->x_nodes[i], mesh->y_nodes[j], mesh->z_nodes[k]);
            }
        }
    }
    MESH_WRITE("$EndNodes\n");

    /* Elements — type 5 = 8-node hex */
    /* Gmsh hex ordering: bottom face CCW (SW,SE,NE,NW) then top face CCW */
    MESH_WRITE("$Elements\n%d\n", nelems);
    int eid = 1;
    for (int k = 0; k < nz; k++) {
        for (int j = 0; j < ny; j++) {
            for (int i = 0; i < nx; i++) {
                int mat = mesh->material_ids[(size_t)k * nxy + (size_t)j * (size_t)nx + (size_t)i];
                int tag = gmsh_material_tag(mesh, mat);
                if (tag <= 0) {
                    alea_set_error_detail(ALEA_ERR_INVALID_STATE,
                                          "mesh material %d is missing from unique material table", mat);
                    return -1;
                }

                int n0 = GMSH_NODE(i,     j,     k,     nx1, ny1);
                int n1 = GMSH_NODE(i + 1, j,     k,     nx1, ny1);
                int n2 = GMSH_NODE(i + 1, j + 1, k,     nx1, ny1);
                int n3 = GMSH_NODE(i,     j + 1, k,     nx1, ny1);
                int n4 = GMSH_NODE(i,     j,     k + 1, nx1, ny1);
                int n5 = GMSH_NODE(i + 1, j,     k + 1, nx1, ny1);
                int n6 = GMSH_NODE(i + 1, j + 1, k + 1, nx1, ny1);
                int n7 = GMSH_NODE(i,     j + 1, k + 1, nx1, ny1);

                /* 5 = hex8, 2 tags: physical_group, elementary_entity */
                MESH_WRITE("%d 5 2 %d %d %d %d %d %d %d %d %d %d\n",
                        eid++, tag, tag,
                        n0, n1, n2, n3, n4, n5, n6, n7);
            }
        }
    }
    MESH_WRITE("$EndElements\n");

    if ((options->fields & ALEA_MESH_EXPORT_MIXED_FLAG) && mesh->mixed_flags) {
        if (write_gmsh_element_data_begin(out, "mixed_flag", nelems) != 0) return -1;
        for (int v = 0; v < nelems; v++) MESH_WRITE("%d %d\n", v + 1, mesh->mixed_flags[v]);
        MESH_WRITE("$EndElementData\n");
    }
    if ((options->fields & ALEA_MESH_EXPORT_DOMINANT_FRACTION) &&
        mesh->dominant_fractions) {
        if (write_gmsh_element_data_begin(out, "dominant_sampled_fraction", nelems) != 0) return -1;
        for (int v = 0; v < nelems; v++)
            MESH_WRITE("%d %.15g\n", v + 1, mesh->dominant_fractions[v]);
        MESH_WRITE("$EndElementData\n");
    }
    if ((options->fields & ALEA_MESH_EXPORT_TIE_FLAG) && mesh->tie_flags) {
        if (write_gmsh_element_data_begin(out, "tie_flag", nelems) != 0) return -1;
        for (int v = 0; v < nelems; v++)
            MESH_WRITE("%d %u\n", v + 1, (unsigned)mesh->tie_flags[v]);
        MESH_WRITE("$EndElementData\n");
    }
    if ((options->fields & ALEA_MESH_EXPORT_SAMPLE_COUNT) && mesh->sample_counts) {
        if (write_gmsh_element_data_begin(out, "sample_count", nelems) != 0) return -1;
        for (int v = 0; v < nelems; v++)
            MESH_WRITE("%d %lu\n", v + 1,
                       (unsigned long)mesh->sample_counts[v]);
        MESH_WRITE("$EndElementData\n");
    }
    if ((options->fields & ALEA_MESH_EXPORT_ESTIMATED_ERROR) &&
        mesh->estimated_errors) {
        if (write_gmsh_element_data_begin(out, "estimated_error", nelems) != 0) return -1;
        for (int v = 0; v < nelems; v++)
            MESH_WRITE("%d %.15g\n", v + 1, mesh->estimated_errors[v]);
        MESH_WRITE("$EndElementData\n");
    }
    if ((options->fields & ALEA_MESH_EXPORT_REFINEMENT_FLAG) &&
        mesh->refinement_flags) {
        if (write_gmsh_element_data_begin(out, "refinement_flag", nelems) != 0) return -1;
        for (int v = 0; v < nelems; v++)
            MESH_WRITE("%d %u\n", v + 1, (unsigned)mesh->refinement_flags[v]);
        MESH_WRITE("$EndElementData\n");
    }
    if (options->fields & ALEA_MESH_EXPORT_MATERIAL_FRACTIONS) {
        char name[96];
        for (int m = 0; m < mesh->num_materials; m++) {
            snprintf(name, sizeof(name), "sampled_fraction_material_%d",
                     mesh->unique_materials[m]);
            if (write_gmsh_element_data_begin(out, name, nelems) != 0) return -1;
            for (int v = 0; v < nelems; v++)
                MESH_WRITE("%d %.15g\n", v + 1,
                           mesh_sampled_fraction_at(mesh, (size_t)v,
                                                    mesh->unique_materials[m]));
            MESH_WRITE("$EndElementData\n");
        }
    }

    return 0;
}

/* ============================================================================
 * VTK Legacy Writer
 * ============================================================================ */

/** Check if nodes are uniformly spaced */
static int is_uniform(const double *nodes, int n, double *out_spacing) {
    if (n < 1) return 0;
    double d = nodes[1] - nodes[0];
    for (int i = 2; i <= n; i++) {
        double di = nodes[i] - nodes[i - 1];
        if (fabs(di - d) > 1e-12 * fabs(d))
            return 0;
    }
    *out_spacing = d;
    return 1;
}

static int write_vtk(const alea_mesh_result_t *mesh, FILE *out,
                     const alea_mesh_export_options_t *options) {
    int nx = mesh->nx, ny = mesh->ny, nz = mesh->nz;
    int nx1 = nx + 1, ny1 = ny + 1, nz1 = nz + 1;
    size_t nelems_sz = 0;
    if (checked_mesh_cell_count(nx, ny, nz, &nelems_sz) != 0) return -1;
    if (nelems_sz > (size_t)INT_MAX) {
        alea_set_error_detail(ALEA_ERR_OVERFLOW,
                              "mesh is too large for VTK legacy export");
        return -1;
    }
    int nelems = (int)nelems_sz;
    size_t nxy = (size_t)nx * (size_t)ny;

    MESH_WRITE("# vtk DataFile Version 3.0\n");
    MESH_WRITE("Alea mesh export\n");
    MESH_WRITE("ASCII\n");

    double dx = 0, dy = 0, dz = 0;
    int ux = is_uniform(mesh->x_nodes, nx, &dx);
    int uy = is_uniform(mesh->y_nodes, ny, &dy);
    int uz = is_uniform(mesh->z_nodes, nz, &dz);

    if (ux && uy && uz) {
        /* STRUCTURED_POINTS — most compact */
        MESH_WRITE("DATASET STRUCTURED_POINTS\n");
        MESH_WRITE("DIMENSIONS %d %d %d\n", nx1, ny1, nz1);
        MESH_WRITE("ORIGIN %.15g %.15g %.15g\n",
                mesh->x_nodes[0], mesh->y_nodes[0], mesh->z_nodes[0]);
        MESH_WRITE("SPACING %.15g %.15g %.15g\n", dx, dy, dz);
    } else {
        /* RECTILINEAR_GRID — non-uniform spacing */
        MESH_WRITE("DATASET RECTILINEAR_GRID\n");
        MESH_WRITE("DIMENSIONS %d %d %d\n", nx1, ny1, nz1);

        MESH_WRITE("X_COORDINATES %d double\n", nx1);
        for (int i = 0; i < nx1; i++) MESH_WRITE("%.15g\n", mesh->x_nodes[i]);

        MESH_WRITE("Y_COORDINATES %d double\n", ny1);
        for (int j = 0; j < ny1; j++) MESH_WRITE("%.15g\n", mesh->y_nodes[j]);

        MESH_WRITE("Z_COORDINATES %d double\n", nz1);
        for (int k = 0; k < nz1; k++) MESH_WRITE("%.15g\n", mesh->z_nodes[k]);
    }

    /* Cell data */
    MESH_WRITE("CELL_DATA %d\n", nelems);

    /* Material IDs */
    MESH_WRITE("SCALARS material_id int 1\n");
    MESH_WRITE("LOOKUP_TABLE default\n");
    /* VTK structured grid cell ordering: x fastest, then y, then z (same as ours) */
    for (int k = 0; k < nz; k++)
        for (int j = 0; j < ny; j++)
            for (int i = 0; i < nx; i++)
                MESH_WRITE("%d\n", mesh->material_ids[(size_t)k * nxy + (size_t)j * (size_t)nx + (size_t)i]);

    /* Cell IDs */
    MESH_WRITE("SCALARS cell_id int 1\n");
    MESH_WRITE("LOOKUP_TABLE default\n");
    for (int k = 0; k < nz; k++)
        for (int j = 0; j < ny; j++)
            for (int i = 0; i < nx; i++)
                MESH_WRITE("%d\n", mesh->cell_ids[(size_t)k * nxy + (size_t)j * (size_t)nx + (size_t)i]);

    if ((options->fields & ALEA_MESH_EXPORT_MIXED_FLAG) && mesh->mixed_flags) {
        MESH_WRITE("SCALARS mixed_flag int 1\n");
        MESH_WRITE("LOOKUP_TABLE default\n");
        for (int k = 0; k < nz; k++)
            for (int j = 0; j < ny; j++)
                for (int i = 0; i < nx; i++)
                    MESH_WRITE("%d\n", (int)mesh->mixed_flags[(size_t)k * nxy + (size_t)j * (size_t)nx + (size_t)i]);
    }

    if ((options->fields & ALEA_MESH_EXPORT_DOMINANT_FRACTION) &&
        mesh->dominant_fractions) {
        MESH_WRITE("SCALARS dominant_sampled_fraction double 1\n");
        MESH_WRITE("LOOKUP_TABLE default\n");
        for (int k = 0; k < nz; k++)
            for (int j = 0; j < ny; j++)
                for (int i = 0; i < nx; i++)
                    MESH_WRITE("%.15g\n", mesh->dominant_fractions[(size_t)k * nxy + (size_t)j * (size_t)nx + (size_t)i]);
    }

    if ((options->fields & ALEA_MESH_EXPORT_TIE_FLAG) && mesh->tie_flags) {
        MESH_WRITE("SCALARS tie_flag int 1\n");
        MESH_WRITE("LOOKUP_TABLE default\n");
        for (int v = 0; v < nelems; v++)
            MESH_WRITE("%u\n", (unsigned)mesh->tie_flags[v]);
    }

    if ((options->fields & ALEA_MESH_EXPORT_SAMPLE_COUNT) && mesh->sample_counts) {
        MESH_WRITE("SCALARS sample_count int 1\n");
        MESH_WRITE("LOOKUP_TABLE default\n");
        for (int v = 0; v < nelems; v++)
            MESH_WRITE("%lu\n", (unsigned long)mesh->sample_counts[v]);
    }

    if ((options->fields & ALEA_MESH_EXPORT_ESTIMATED_ERROR) &&
        mesh->estimated_errors) {
        MESH_WRITE("SCALARS estimated_error double 1\n");
        MESH_WRITE("LOOKUP_TABLE default\n");
        for (int v = 0; v < nelems; v++)
            MESH_WRITE("%.15g\n", mesh->estimated_errors[v]);
    }

    if ((options->fields & ALEA_MESH_EXPORT_REFINEMENT_FLAG) &&
        mesh->refinement_flags) {
        MESH_WRITE("SCALARS refinement_flag int 1\n");
        MESH_WRITE("LOOKUP_TABLE default\n");
        for (int v = 0; v < nelems; v++)
            MESH_WRITE("%u\n", (unsigned)mesh->refinement_flags[v]);
    }

    if (options->fields & ALEA_MESH_EXPORT_MATERIAL_FRACTIONS) {
        for (int m = 0; m < mesh->num_materials; m++) {
            MESH_WRITE("SCALARS sampled_fraction_material_%d double 1\n",
                       mesh->unique_materials[m]);
            MESH_WRITE("LOOKUP_TABLE default\n");
            for (int v = 0; v < nelems; v++)
                MESH_WRITE("%.15g\n", mesh_sampled_fraction_at(
                    mesh, (size_t)v, mesh->unique_materials[m]));
        }
    }

    return 0;
}

/* ============================================================================
 * Export API
 * ============================================================================ */

static int mesh_result_material_known(const alea_mesh_result_t *mesh,
                                      int material_id) {
    for (int i = 0; i < mesh->num_materials; i++)
        if (mesh->unique_materials[i] == material_id) return 1;
    return 0;
}

static int validate_mesh_result(const alea_mesh_result_t *mesh) {
    if (!mesh) {
        alea_set_error_detail(ALEA_ERR_NULL_ARG, "mesh result is NULL");
        return -1;
    }

    if (mesh->nx <= 0 || mesh->ny <= 0 || mesh->nz <= 0) {
        alea_set_error_detail(ALEA_ERR_INVALID_ARG,
                              "mesh result dimensions must be positive");
        return -1;
    }
    size_t ncells = 0, nnodes = 0;
    if (checked_mesh_cell_count(mesh->nx, mesh->ny, mesh->nz, &ncells) != 0 ||
        checked_mesh_node_count(mesh->nx, mesh->ny, mesh->nz, &nnodes) != 0)
        return -1;
    if (ncells > (size_t)INT_MAX || nnodes > (size_t)INT_MAX) {
        alea_set_error_detail(ALEA_ERR_OVERFLOW,
                              "mesh result exceeds legacy exporter limits");
        return -1;
    }
    if (!mesh->x_nodes || !mesh->y_nodes || !mesh->z_nodes ||
        !mesh->material_ids || !mesh->cell_ids) {
        alea_set_error_detail(ALEA_ERR_INVALID_ARG,
                              "mesh result is missing required arrays");
        return -1;
    }
    if (validate_nodes(mesh->x_nodes, mesh->nx, "X") != 0 ||
        validate_nodes(mesh->y_nodes, mesh->ny, "Y") != 0 ||
        validate_nodes(mesh->z_nodes, mesh->nz, "Z") != 0)
        return -1;
    if (mesh->num_materials <= 0 || !mesh->unique_materials) {
        alea_set_error_detail(ALEA_ERR_INVALID_ARG,
                              "mesh result has no material table");
        return -1;
    }
    for (int i = 0; i < mesh->num_materials; i++) {
        if (i > 0 && mesh->unique_materials[i] <= mesh->unique_materials[i - 1]) {
            alea_set_error_detail(ALEA_ERR_INVALID_ARG,
                                  "mesh material table is not strictly increasing");
            return -1;
        }
    }
    for (size_t v = 0; v < ncells; v++) {
        if (!mesh_result_material_known(mesh, mesh->material_ids[v])) {
            alea_set_error_detail(ALEA_ERR_INVALID_ARG,
                                  "mesh voxel material is absent from material table");
            return -1;
        }
        if (mesh->dominant_fractions &&
            (!isfinite(mesh->dominant_fractions[v]) ||
             mesh->dominant_fractions[v] < 0.0 ||
             mesh->dominant_fractions[v] > 1.0)) {
            alea_set_error_detail(ALEA_ERR_INVALID_ARG,
                                  "mesh dominant sampled fraction is invalid");
            return -1;
        }
        if (mesh->sample_counts && mesh->sample_counts[v] == 0) {
            alea_set_error_detail(ALEA_ERR_INVALID_ARG,
                                  "mesh voxel has zero point samples");
            return -1;
        }
        if (mesh->estimated_errors &&
            (!isfinite(mesh->estimated_errors[v]) ||
             mesh->estimated_errors[v] < 0.0 || mesh->estimated_errors[v] > 1.0)) {
            alea_set_error_detail(ALEA_ERR_INVALID_ARG,
                                  "mesh voxel has invalid estimated error");
            return -1;
        }
        if (mesh->refinement_flags &&
            (mesh->refinement_flags[v] & ~ALEA_MESH_REFINEMENT_LIMIT_REACHED)) {
            alea_set_error_detail(ALEA_ERR_INVALID_ARG,
                                  "mesh voxel has unknown refinement flags");
            return -1;
        }
        if (mesh->tie_flags &&
            (mesh->tie_flags[v] & ~(ALEA_MESH_TIE_MATERIAL | ALEA_MESH_TIE_CELL))) {
            alea_set_error_detail(ALEA_ERR_INVALID_ARG,
                                  "mesh voxel has unknown tie flags");
            return -1;
        }
    }

    if (!!mesh->fraction_spans != !!mesh->fractions ||
        ((mesh->fraction_spans || mesh->fractions) && mesh->fraction_count == 0)) {
        alea_set_error_detail(ALEA_ERR_INVALID_ARG,
                              "mesh sampled-fraction storage is inconsistent");
        return -1;
    }
    if (mesh->fraction_spans) {
        for (size_t v = 0; v < ncells; v++) {
            const alea_mesh_fraction_span_t span = mesh->fraction_spans[v];
            const size_t begin = span.offset;
            const size_t count = span.count;
            if (count == 0 || begin > mesh->fraction_count ||
                count > mesh->fraction_count - begin) {
                alea_set_error_detail(ALEA_ERR_INVALID_ARG,
                                      "mesh sampled-fraction span is invalid");
                return -1;
            }
            double sum = 0.0;
            for (size_t i = 0; i < count; i++) {
                const alea_mesh_material_fraction_t *f = &mesh->fractions[begin + i];
                if (!mesh_result_material_known(mesh, f->material_id) ||
                    !isfinite(f->fraction) || f->fraction < 0.0 || f->fraction > 1.0) {
                    alea_set_error_detail(ALEA_ERR_INVALID_ARG,
                                          "mesh sampled-fraction entry is invalid");
                    return -1;
                }
                sum += f->fraction;
            }
            if (fabs(sum - 1.0) > 1e-12) {
                alea_set_error_detail(ALEA_ERR_INVALID_ARG,
                                      "mesh sampled fractions do not sum to one");
                return -1;
            }
        }
    }
    if (!!mesh->cell_fraction_spans != !!mesh->cell_fractions ||
        ((mesh->cell_fraction_spans || mesh->cell_fractions) &&
         mesh->cell_fraction_count == 0)) {
        alea_set_error_detail(ALEA_ERR_INVALID_ARG,
                              "mesh cell-fraction storage is inconsistent");
        return -1;
    }
    if (mesh->cell_fraction_spans) {
        for (size_t v = 0; v < ncells; v++) {
            const alea_mesh_fraction_span_t span = mesh->cell_fraction_spans[v];
            const size_t begin = span.offset;
            const size_t count = span.count;
            if (count == 0 || begin > mesh->cell_fraction_count ||
                count > mesh->cell_fraction_count - begin) {
                alea_set_error_detail(ALEA_ERR_INVALID_ARG,
                                      "mesh cell-fraction span is invalid");
                return -1;
            }
            double sum = 0.0;
            for (size_t i = 0; i < count; i++) {
                const alea_mesh_cell_fraction_t *f =
                    &mesh->cell_fractions[begin + i];
                if (!mesh_result_material_known(mesh, f->material_id) ||
                    !isfinite(f->fraction) || f->fraction < 0.0 ||
                    f->fraction > 1.0) {
                    alea_set_error_detail(ALEA_ERR_INVALID_ARG,
                                          "mesh cell-fraction entry is invalid");
                    return -1;
                }
                sum += f->fraction;
            }
            if (fabs(sum - 1.0) > 1e-12) {
                alea_set_error_detail(ALEA_ERR_INVALID_ARG,
                                      "mesh cell fractions do not sum to one");
                return -1;
            }
        }
    }
    return 0;
}

static int validate_export_options(const alea_mesh_result_t *mesh,
                                   const alea_mesh_export_options_t *options) {
    const uint32_t known = ALEA_MESH_EXPORT_MIXED_FLAG |
                           ALEA_MESH_EXPORT_DOMINANT_FRACTION |
                           ALEA_MESH_EXPORT_TIE_FLAG |
                           ALEA_MESH_EXPORT_SAMPLE_COUNT |
                           ALEA_MESH_EXPORT_MATERIAL_FRACTIONS |
                           ALEA_MESH_EXPORT_ESTIMATED_ERROR |
                           ALEA_MESH_EXPORT_REFINEMENT_FLAG;
    if (!options) {
        alea_set_error_detail(ALEA_ERR_NULL_ARG,
                              "mesh export options are NULL");
        return -1;
    }
    if (options->fields & ~known) {
        alea_set_error_detail(ALEA_ERR_INVALID_ARG,
                              "mesh export options contain unknown fields");
        return -1;
    }
    if (options->max_fraction_materials < 0) {
        alea_set_error_detail(ALEA_ERR_INVALID_ARG,
                              "mesh export material limit is negative");
        return -1;
    }
    if ((options->fields & ALEA_MESH_EXPORT_MATERIAL_FRACTIONS) &&
        options->max_fraction_materials > 0 &&
        mesh->num_materials > options->max_fraction_materials) {
        alea_set_error_detail(ALEA_ERR_INVALID_ARG,
                              "mesh fraction export exceeds material-array limit");
        return -1;
    }
    return 0;
}

int alea_mesh_export_stream_ex(const alea_mesh_result_t *mesh,
                               alea_mesh_format_t fmt, FILE *out,
                               const alea_mesh_export_options_t *options) {
    if (!out) {
        alea_set_error_detail(ALEA_ERR_NULL_ARG, "mesh export stream is NULL");
        return -1;
    }
    if (validate_mesh_result(mesh) != 0) return -1;
    if (validate_export_options(mesh, options) != 0) return -1;
    int rc = -1;
    switch (fmt) {
        case ALEA_MESH_GMSH: rc = write_gmsh(mesh, out, options); break;
        case ALEA_MESH_VTK:  rc = write_vtk(mesh, out, options); break;
        default:
            alea_set_error_detail(ALEA_ERR_INVALID_ARG,
                                  "mesh export format is invalid");
            return -1;
    }
    if (rc != 0) return -1;
    if (fflush(out) != 0 || ferror(out)) {
        alea_set_error_detail(ALEA_ERR_FILE_WRITE,
                              "mesh stream flush failed: %s", strerror(errno));
        return -1;
    }
    return 0;
}

int alea_mesh_export_stream(const alea_mesh_result_t *mesh,
                            alea_mesh_format_t fmt, FILE *out) {
    alea_mesh_export_options_t options;
    alea_mesh_export_options_init(&options);
    return alea_mesh_export_stream_ex(mesh, fmt, out, &options);
}

int alea_mesh_export_ex(const alea_mesh_result_t *mesh,
                        alea_mesh_format_t fmt, const char *filename,
                        const alea_mesh_export_options_t *options) {
    if (!filename) {
        alea_set_error_detail(ALEA_ERR_NULL_ARG, "mesh export filename is NULL");
        return -1;
    }
    if (validate_mesh_result(mesh) != 0) return -1;
    if (validate_export_options(mesh, options) != 0) return -1;
    size_t temp_path_size = 0;
    if (checked_add_size(strlen(filename), 32, &temp_path_size) != 0) return -1;
    char *temp_path = malloc(temp_path_size);
    if (!temp_path) {
        alea_set_error_detail(ALEA_ERR_OUT_OF_MEMORY,
                              "failed to allocate mesh temporary path");
        return -1;
    }
    FILE *f = alea_sibling_tmpfile(filename, temp_path, temp_path_size);
    if (!f) {
        alea_set_error_detail(ALEA_ERR_FILE_WRITE,
                              "failed to create mesh export temporary file for '%s': %s",
                              filename, strerror(errno));
        free(temp_path);
        return -1;
    }
    int rc = alea_mesh_export_stream_ex(mesh, fmt, f, options);
    if (fclose(f) != 0) {
        alea_set_error_detail(ALEA_ERR_FILE_WRITE,
                              "failed to close mesh export file '%s': %s",
                              filename, strerror(errno));
        rc = -1;
    }
    if (rc == 0 && alea_replace_file(temp_path, filename) != 0) {
        alea_set_error_detail(ALEA_ERR_FILE_WRITE,
                              "failed to replace mesh export file '%s': %s",
                              filename, strerror(errno));
        rc = -1;
    }
    if (rc != 0) remove(temp_path);
    free(temp_path);
    return rc;
}

int alea_mesh_export(const alea_mesh_result_t *mesh,
                     alea_mesh_format_t fmt, const char *filename) {
    alea_mesh_export_options_t options;
    alea_mesh_export_options_init(&options);
    return alea_mesh_export_ex(mesh, fmt, filename, &options);
}

int alea_mesh_export_system(alea_system_t *sys,
                             const alea_mesh_config_t *cfg,
                             const char *filename) {
    if (!sys || !cfg || !filename) return -1;

    alea_mesh_result_t *mesh = alea_mesh_sample(sys, cfg);
    if (!mesh) return -1;

    int rc = alea_mesh_export(mesh, cfg->format, filename);
    alea_mesh_result_free(mesh);
    return rc;
}
