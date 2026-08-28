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
#include "util/compat.h"

#include <errno.h>
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include <stdint.h>

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
        cfg->sampling_mode != ALEA_MESH_SAMPLE_SUBCELL) {
        alea_set_error_detail(ALEA_ERR_INVALID_ARG,
                              "mesh sampling_mode is invalid");
        return -1;
    }
    if (cfg->subsamples_per_axis <= 0 || cfg->subsamples_per_axis > 8) {
        alea_set_error_detail(ALEA_ERR_INVALID_ARG,
                              "mesh subsamples_per_axis must be in [1,8]");
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

void alea_mesh_export_options_init(alea_mesh_export_options_t *options) {
    if (!options) return;
    options->fields = ALEA_MESH_EXPORT_MIXED_FLAG |
                      ALEA_MESH_EXPORT_DOMINANT_FRACTION |
                      ALEA_MESH_EXPORT_TIE_FLAG |
                      ALEA_MESH_EXPORT_SAMPLE_COUNT;
    options->max_fraction_materials = 64;
}

/** Collect every material present in the sampled-fraction entries. */
static int *collect_unique_fractions(
    const alea_mesh_material_fraction_t *fractions,
    size_t count, int *out_num) {
    if (count == 0) { *out_num = 0; return NULL; }
    if (count > (size_t)INT_MAX) {
        alea_set_error_detail(ALEA_ERR_OVERFLOW,
                              "mesh has too many cells for material collection");
        *out_num = 0;
        return NULL;
    }

    size_t bytes = 0;
    if (checked_mul_size(count, sizeof(int), &bytes) != 0) {
        *out_num = 0;
        return NULL;
    }
    int *tmp = malloc(bytes);
    if (!tmp) {
        alea_set_error_detail(ALEA_ERR_OUT_OF_MEMORY,
                              "mesh failed to collect unique materials");
        *out_num = 0;
        return NULL;
    }
    for (size_t i = 0; i < count; i++)
        tmp[i] = fractions[i].material_id;
    qsort(tmp, count, sizeof(*tmp), compare_ints);

    /* dedup */
    size_t n = 1;
    for (size_t i = 1; i < count; i++) {
        if (tmp[i] != tmp[i - 1])
            tmp[n++] = tmp[i];
    }

    int *result = realloc(tmp, (size_t)n * sizeof(int));
    if (!result) result = tmp;  /* realloc to smaller should not fail, but just in case */
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

static int mesh_sample_voxel_materials(alea_system_t *sys,
                                       const alea_mesh_config_t *cfg,
                                       const double *xn,
                                       const double *yn,
                                       const double *zn,
                                       int i, int j, int k,
                                       int *materials,
                                       int *counts,
                                       int *out_num_materials,
                                       int *owner_materials,
                                       int *owner_cells,
                                       int *owner_counts,
                                       int *out_num_owners,
                                       int *out_num_samples) {
    int n = cfg->subsamples_per_axis;
    if (cfg->sampling_mode == ALEA_MESH_SAMPLE_CENTER)
        n = 1;
    else if (cfg->sampling_mode == ALEA_MESH_SAMPLE_CORNERS)
        n = 2;
    int nsamples = n * n * n;
    int n_materials = 0;
    int n_owners = 0;

    double x0 = xn[i], x1 = xn[i + 1];
    double y0 = yn[j], y1 = yn[j + 1];
    double z0 = zn[k], z1 = zn[k + 1];
    double dx = x1 - x0, dy = y1 - y0, dz = z1 - z0;
    double eps = 1e-9;

    memset(counts, 0, 512 * sizeof(int));
    memset(owner_counts, 0, 512 * sizeof(int));

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

    alea_mesh_material_fraction_t *tmp =
        realloc(*fractions, new_cap * sizeof(alea_mesh_material_fraction_t));
    if (!tmp) {
        alea_set_error_detail(ALEA_ERR_OUT_OF_MEMORY,
                              "mesh failed to grow material fraction storage");
        return -1;
    }

    *fractions = tmp;
    *capacity = new_cap;
    return 0;
}

/* ============================================================================
 * Sampling
 * ============================================================================ */

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

    if (all_bounds_zero && needs_uniform_axis) {
        /* Auto-detect from bounding sphere */
        double cx, cy, cz, r;
        if (alea_compute_bounding_sphere(sys, 1.0, &cx, &cy, &cz, &r) != 0) {
            alea_set_error_detail(ALEA_ERR_INVALID_STATE,
                                  "mesh auto-bounds failed; provide explicit bounds");
            return NULL;
        }
        double pad = r * cfg->auto_pad;
        xlo = cx - r - pad;  xhi = cx + r + pad;
        ylo = cy - r - pad;  yhi = cy + r + pad;
        zlo = cz - r - pad;  zhi = cz + r + pad;
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

    int *mat_ids = malloc(ncells * sizeof(int));
    int *cell_ids = malloc(ncells * sizeof(int));
    if (!mat_ids || !cell_ids) {
        alea_set_error_detail(ALEA_ERR_OUT_OF_MEMORY,
                              "mesh failed to allocate cell data arrays");
        free(xn); free(yn); free(zn); free(mat_ids); free(cell_ids);
        return NULL;
    }

    unsigned char *mixed_flags = NULL;
    double *dominant_fractions = NULL;
    uint32_t *sample_counts = NULL;
    uint8_t *tie_flags = NULL;
    alea_mesh_fraction_span_t *fraction_spans = NULL;
    alea_mesh_material_fraction_t *fractions = NULL;
    size_t fraction_count = 0;
    size_t fraction_capacity = 0;
    mixed_flags = calloc(ncells, sizeof(unsigned char));
    dominant_fractions = malloc(ncells * sizeof(double));
    sample_counts = malloc(ncells * sizeof(uint32_t));
    tie_flags = calloc(ncells, sizeof(uint8_t));
    fraction_spans = calloc(ncells, sizeof(alea_mesh_fraction_span_t));
    if (!mixed_flags || !dominant_fractions || !sample_counts || !tie_flags ||
        !fraction_spans) {
        alea_set_error_detail(ALEA_ERR_OUT_OF_MEMORY,
                              "mesh failed to allocate composition arrays");
        free(xn); free(yn); free(zn);
        free(mat_ids); free(cell_ids);
        free(mixed_flags); free(dominant_fractions); free(sample_counts);
        free(tie_flags); free(fraction_spans);
        return NULL;
    }

    /* Ensure universe lookup is ready before point sampling. Mesh sampling
     * uses recursive point lookup, not the expanded flat spatial index. */
    if (alea_system_prepare_query_caches(sys, ALEA_CACHE_UNIVERSE) != 0) {
        free(xn); free(yn); free(zn);
        free(mat_ids); free(cell_ids);
        free(mixed_flags); free(dominant_fractions); free(sample_counts);
        free(tie_flags); free(fraction_spans);
        return NULL;
    }

    /* Sample voxel composition and coherent dominant cell/material IDs. */
    size_t nxy = (size_t)nx * (size_t)ny;
    int mixed_count = 0;

    for (int k = 0; k < nz; k++) {
        for (int j = 0; j < ny; j++) {
            for (int i = 0; i < nx; i++) {
                size_t idx = (size_t)k * nxy + (size_t)j * (size_t)nx + (size_t)i;
                int materials[512];
                int counts[512];
                int owner_materials[512];
                int owner_cells[512];
                int owner_counts[512];
                int n_materials = 0;
                int n_owners = 0;
                int nsamples = 0;
                int dominant_material = cfg->void_material_id;
                int dominant_cell = -1;
                unsigned char mixed = 0;
                double dominant = 1.0;

                mesh_sample_voxel_materials(sys, cfg, xn, yn, zn, i, j, k,
                                            materials, counts,
                                            &n_materials, owner_materials,
                                            owner_cells, owner_counts,
                                            &n_owners, &nsamples);
                mesh_fraction_stats(materials, counts, n_materials, nsamples, cfg,
                                    owner_materials, owner_cells, owner_counts,
                                    n_owners, &dominant_material, &dominant_cell,
                                    &mixed, &dominant, &tie_flags[idx]);
                mat_ids[idx] = dominant_material;
                cell_ids[idx] = dominant_cell;
                sample_counts[idx] = (uint32_t)nsamples;

                if (fraction_count > (size_t)UINT32_MAX ||
                    (size_t)n_materials > (size_t)UINT32_MAX - fraction_count) {
                    alea_set_error_detail(ALEA_ERR_OVERFLOW,
                                          "mesh fraction span offset exceeds uint32_t");
                    free(xn); free(yn); free(zn);
                    free(mat_ids); free(cell_ids);
                    free(mixed_flags); free(dominant_fractions); free(sample_counts);
                    free(tie_flags);
                    free(fraction_spans); free(fractions);
                    return NULL;
                }
                if (ensure_fraction_capacity(&fractions, &fraction_capacity,
                                             fraction_count + (size_t)n_materials) != 0) {
                    free(xn); free(yn); free(zn);
                    free(mat_ids); free(cell_ids);
                    free(mixed_flags); free(dominant_fractions); free(sample_counts);
                    free(tie_flags);
                    free(fraction_spans); free(fractions);
                    return NULL;
                }

                fraction_spans[idx].offset = (uint32_t)fraction_count;
                fraction_spans[idx].count = (uint32_t)n_materials;
                for (int m = 0; m < n_materials; m++) {
                    fractions[fraction_count].material_id = materials[m];
                    fractions[fraction_count].fraction =
                        nsamples > 0 ? (double)counts[m] / (double)nsamples : 0.0;
                    fraction_count++;
                }

                mixed_flags[idx] = mixed;
                dominant_fractions[idx] = dominant;
                mixed_count += mixed ? 1 : 0;
            }
        }
    }

    /* Collect unique materials */
    int num_mats = 0;
    int *unique_mats = collect_unique_fractions(fractions, fraction_count,
                                               &num_mats);
    if (!unique_mats) {
        free(xn); free(yn); free(zn);
        free(mat_ids); free(cell_ids);
        free(mixed_flags); free(dominant_fractions); free(sample_counts);
        free(tie_flags);
        free(fraction_spans); free(fractions);
        return NULL;
    }

    /* Build result */
    alea_mesh_result_t *res = calloc(1, sizeof(*res));
    if (!res) {
        free(xn); free(yn); free(zn); free(mat_ids); free(cell_ids); free(unique_mats);
        free(mixed_flags); free(dominant_fractions); free(sample_counts);
        free(tie_flags);
        free(fraction_spans); free(fractions);
        return NULL;
    }
    res->nx = nx;  res->ny = ny;  res->nz = nz;
    res->x_nodes = xn;  res->y_nodes = yn;  res->z_nodes = zn;
    res->material_ids = mat_ids;
    res->cell_ids = cell_ids;
    res->num_materials = num_mats;
    res->unique_materials = unique_mats;
    res->mixed_flags = mixed_flags;
    res->dominant_fractions = dominant_fractions;
    res->sample_counts = sample_counts;
    res->tie_flags = tie_flags;
    res->mixed_count = mixed_count;
    res->fraction_spans = fraction_spans;
    res->fractions = fractions;
    res->fraction_count = fraction_count;
    return res;
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
    free(mesh->sample_counts);
    free(mesh->tie_flags);
    free(mesh->fraction_spans);
    free(mesh->fractions);
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
    return 0;
}

static int validate_export_options(const alea_mesh_result_t *mesh,
                                   const alea_mesh_export_options_t *options) {
    const uint32_t known = ALEA_MESH_EXPORT_MIXED_FLAG |
                           ALEA_MESH_EXPORT_DOMINANT_FRACTION |
                           ALEA_MESH_EXPORT_TIE_FLAG |
                           ALEA_MESH_EXPORT_SAMPLE_COUNT |
                           ALEA_MESH_EXPORT_MATERIAL_FRACTIONS;
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
