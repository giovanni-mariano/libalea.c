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
    if (cfg->format != ALEA_MESH_GMSH && cfg->format != ALEA_MESH_VTK) {
        alea_set_error_detail(ALEA_ERR_INVALID_ARG,
                              "mesh format is invalid");
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

/** Collect unique values from array, return sorted array and count */
static int *collect_unique(const int *arr, size_t count, int *out_num) {
    if (count == 0) { *out_num = 0; return NULL; }
    if (count > (size_t)INT_MAX) {
        alea_set_error_detail(ALEA_ERR_OVERFLOW,
                              "mesh has too many cells for material collection");
        *out_num = 0;
        return NULL;
    }

    int *tmp = malloc(count * sizeof(int));
    if (!tmp) {
        alea_set_error_detail(ALEA_ERR_OUT_OF_MEMORY,
                              "mesh failed to collect unique materials");
        *out_num = 0;
        return NULL;
    }
    memcpy(tmp, arr, count * sizeof(int));

    /* sort */
    for (size_t i = 1; i < count; i++) {
        int key = tmp[i];
        size_t j = i;
        while (j > 0 && tmp[j - 1] > key) {
            tmp[j] = tmp[j - 1];
            j--;
        }
        tmp[j] = key;
    }

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

static int sample_material_at_point(alea_system_t *sys,
                                    double x, double y, double z,
                                    int void_mat) {
    int cell_idx = alea_identify_cell_at_point(sys, x, y, z);
    if (cell_idx < 0) return void_mat;
    return sys->cells.data[cell_idx].material_id;
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

static int mesh_sample_voxel_materials(alea_system_t *sys,
                                       const alea_mesh_config_t *cfg,
                                       const double *xn,
                                       const double *yn,
                                       const double *zn,
                                       int i, int j, int k,
                                       int *materials,
                                       int *counts,
                                       int *out_num_materials,
                                       int *out_num_samples) {
    int n = cfg->subsamples_per_axis;
    if (cfg->sampling_mode == ALEA_MESH_SAMPLE_CENTER)
        n = 1;
    else if (cfg->sampling_mode == ALEA_MESH_SAMPLE_CORNERS)
        n = 2;
    int nsamples = n * n * n;
    int n_materials = 0;

    double x0 = xn[i], x1 = xn[i + 1];
    double y0 = yn[j], y1 = yn[j + 1];
    double z0 = zn[k], z1 = zn[k + 1];
    double dx = x1 - x0, dy = y1 - y0, dz = z1 - z0;
    double eps = 1e-9;

    memset(counts, 0, 512 * sizeof(int));

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
                int mat = sample_material_at_point(sys, x, y, z,
                                                   cfg->void_material_id);
                tally_material(mat, materials, counts, &n_materials);
            }
        }
    }

    *out_num_materials = n_materials;
    *out_num_samples = nsamples;
    return 0;
}

static void mesh_fraction_stats(const int *materials,
                                const int *counts,
                                int num_materials,
                                int num_samples,
                                const alea_mesh_config_t *cfg,
                                int *out_dominant_material,
                                unsigned char *out_mixed,
                                double *out_dominant_fraction) {
    int max_count = 0;
    int dominant_material = cfg->void_material_id;
    for (int m = 0; m < num_materials; m++) {
        if (counts[m] > max_count) {
            max_count = counts[m];
            dominant_material = materials[m];
        }
    }

    double dominant = num_samples > 0 ? (double)max_count / (double)num_samples : 1.0;
    *out_dominant_material = dominant_material;
    *out_dominant_fraction = dominant;
    *out_mixed = (dominant < 1.0 - cfg->mixed_threshold) ? 1 : 0;
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
    alea_mesh_fraction_span_t *fraction_spans = NULL;
    alea_mesh_material_fraction_t *fractions = NULL;
    size_t fraction_count = 0;
    size_t fraction_capacity = 0;
    mixed_flags = calloc(ncells, sizeof(unsigned char));
    dominant_fractions = malloc(ncells * sizeof(double));
    fraction_spans = calloc(ncells, sizeof(alea_mesh_fraction_span_t));
    if (!mixed_flags || !dominant_fractions || !fraction_spans) {
        alea_set_error_detail(ALEA_ERR_OUT_OF_MEMORY,
                              "mesh failed to allocate composition arrays");
        free(xn); free(yn); free(zn);
        free(mat_ids); free(cell_ids);
        free(mixed_flags); free(dominant_fractions); free(fraction_spans);
        return NULL;
    }

    /* Ensure universe lookup is ready before point sampling. Mesh sampling
     * uses recursive point lookup, not the expanded flat spatial index. */
    if (alea_system_prepare_query_caches(sys, ALEA_CACHE_UNIVERSE) != 0) {
        free(xn); free(yn); free(zn);
        free(mat_ids); free(cell_ids);
        free(mixed_flags); free(dominant_fractions); free(fraction_spans);
        return NULL;
    }

    /* Sample voxel composition and center cell IDs. */
    size_t nxy = (size_t)nx * (size_t)ny;
    int mixed_count = 0;

    for (int k = 0; k < nz; k++) {
        for (int j = 0; j < ny; j++) {
            for (int i = 0; i < nx; i++) {
                double xc = (xn[i] + xn[i + 1]) * 0.5;
                double yc = (yn[j] + yn[j + 1]) * 0.5;
                double zc = (zn[k] + zn[k + 1]) * 0.5;
                size_t idx = (size_t)k * nxy + (size_t)j * (size_t)nx + (size_t)i;
                int materials[512];
                int counts[512];
                int n_materials = 0;
                int nsamples = 0;
                int dominant_material = cfg->void_material_id;
                unsigned char mixed = 0;
                double dominant = 1.0;

                int cell_idx = alea_identify_cell_at_point(sys, xc, yc, zc);
                if (cell_idx < 0) {
                    cell_ids[idx] = -1;
                } else {
                    cell_ids[idx] = sys->cells.data[cell_idx].mc_cell_id;
                }

                mesh_sample_voxel_materials(sys, cfg, xn, yn, zn, i, j, k,
                                            materials, counts,
                                            &n_materials, &nsamples);
                mesh_fraction_stats(materials, counts, n_materials, nsamples, cfg,
                                    &dominant_material, &mixed, &dominant);
                mat_ids[idx] = dominant_material;

                if (fraction_count > (size_t)UINT32_MAX ||
                    (size_t)n_materials > (size_t)UINT32_MAX - fraction_count) {
                    alea_set_error_detail(ALEA_ERR_OVERFLOW,
                                          "mesh fraction span offset exceeds uint32_t");
                    free(xn); free(yn); free(zn);
                    free(mat_ids); free(cell_ids);
                    free(mixed_flags); free(dominant_fractions);
                    free(fraction_spans); free(fractions);
                    return NULL;
                }
                if (ensure_fraction_capacity(&fractions, &fraction_capacity,
                                             fraction_count + (size_t)n_materials) != 0) {
                    free(xn); free(yn); free(zn);
                    free(mat_ids); free(cell_ids);
                    free(mixed_flags); free(dominant_fractions);
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
    int *unique_mats = collect_unique(mat_ids, ncells, &num_mats);
    if (!unique_mats) {
        free(xn); free(yn); free(zn);
        free(mat_ids); free(cell_ids);
        free(mixed_flags); free(dominant_fractions);
        free(fraction_spans); free(fractions);
        return NULL;
    }

    /* Build result */
    alea_mesh_result_t *res = calloc(1, sizeof(*res));
    if (!res) {
        free(xn); free(yn); free(zn); free(mat_ids); free(cell_ids); free(unique_mats);
        free(mixed_flags); free(dominant_fractions);
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
    free(mesh->fraction_spans);
    free(mesh->fractions);
    free(mesh);
}

/* ============================================================================
 * Gmsh .msh v2.2 ASCII Writer
 * ============================================================================ */

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

static int write_gmsh(const alea_mesh_result_t *mesh, FILE *out) {
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
    fprintf(out, "$MeshFormat\n2.2 0 8\n$EndMeshFormat\n");

    /* Physical names use positive sequential tags; names keep original IDs. */
    fprintf(out, "$PhysicalNames\n%d\n", mesh->num_materials);
    for (int m = 0; m < mesh->num_materials; m++) {
        fprintf(out, "3 %d \"material_%d\"\n",
                m + 1, mesh->unique_materials[m]);
    }
    fprintf(out, "$EndPhysicalNames\n");

    /* Nodes */
    fprintf(out, "$Nodes\n%d\n", nnodes);
    for (int k = 0; k < nz1; k++) {
        for (int j = 0; j < ny1; j++) {
            for (int i = 0; i < nx1; i++) {
                int nid = GMSH_NODE(i, j, k, nx1, ny1);
                fprintf(out, "%d %.15g %.15g %.15g\n",
                        nid, mesh->x_nodes[i], mesh->y_nodes[j], mesh->z_nodes[k]);
            }
        }
    }
    fprintf(out, "$EndNodes\n");

    /* Elements — type 5 = 8-node hex */
    /* Gmsh hex ordering: bottom face CCW (SW,SE,NE,NW) then top face CCW */
    fprintf(out, "$Elements\n%d\n", nelems);
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
                fprintf(out, "%d 5 2 %d %d %d %d %d %d %d %d %d %d\n",
                        eid++, tag, tag,
                        n0, n1, n2, n3, n4, n5, n6, n7);
            }
        }
    }
    fprintf(out, "$EndElements\n");

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

static int write_vtk(const alea_mesh_result_t *mesh, FILE *out) {
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

    fprintf(out, "# vtk DataFile Version 3.0\n");
    fprintf(out, "Alea mesh export\n");
    fprintf(out, "ASCII\n");

    double dx = 0, dy = 0, dz = 0;
    int ux = is_uniform(mesh->x_nodes, nx, &dx);
    int uy = is_uniform(mesh->y_nodes, ny, &dy);
    int uz = is_uniform(mesh->z_nodes, nz, &dz);

    if (ux && uy && uz) {
        /* STRUCTURED_POINTS — most compact */
        fprintf(out, "DATASET STRUCTURED_POINTS\n");
        fprintf(out, "DIMENSIONS %d %d %d\n", nx1, ny1, nz1);
        fprintf(out, "ORIGIN %.15g %.15g %.15g\n",
                mesh->x_nodes[0], mesh->y_nodes[0], mesh->z_nodes[0]);
        fprintf(out, "SPACING %.15g %.15g %.15g\n", dx, dy, dz);
    } else {
        /* RECTILINEAR_GRID — non-uniform spacing */
        fprintf(out, "DATASET RECTILINEAR_GRID\n");
        fprintf(out, "DIMENSIONS %d %d %d\n", nx1, ny1, nz1);

        fprintf(out, "X_COORDINATES %d double\n", nx1);
        for (int i = 0; i < nx1; i++) fprintf(out, "%.15g\n", mesh->x_nodes[i]);

        fprintf(out, "Y_COORDINATES %d double\n", ny1);
        for (int j = 0; j < ny1; j++) fprintf(out, "%.15g\n", mesh->y_nodes[j]);

        fprintf(out, "Z_COORDINATES %d double\n", nz1);
        for (int k = 0; k < nz1; k++) fprintf(out, "%.15g\n", mesh->z_nodes[k]);
    }

    /* Cell data */
    fprintf(out, "CELL_DATA %d\n", nelems);

    /* Material IDs */
    fprintf(out, "SCALARS material_id int 1\n");
    fprintf(out, "LOOKUP_TABLE default\n");
    /* VTK structured grid cell ordering: x fastest, then y, then z (same as ours) */
    for (int k = 0; k < nz; k++)
        for (int j = 0; j < ny; j++)
            for (int i = 0; i < nx; i++)
                fprintf(out, "%d\n", mesh->material_ids[(size_t)k * nxy + (size_t)j * (size_t)nx + (size_t)i]);

    /* Cell IDs */
    fprintf(out, "SCALARS cell_id int 1\n");
    fprintf(out, "LOOKUP_TABLE default\n");
    for (int k = 0; k < nz; k++)
        for (int j = 0; j < ny; j++)
            for (int i = 0; i < nx; i++)
                fprintf(out, "%d\n", mesh->cell_ids[(size_t)k * nxy + (size_t)j * (size_t)nx + (size_t)i]);

    if (mesh->mixed_flags) {
        fprintf(out, "SCALARS mixed_flag int 1\n");
        fprintf(out, "LOOKUP_TABLE default\n");
        for (int k = 0; k < nz; k++)
            for (int j = 0; j < ny; j++)
                for (int i = 0; i < nx; i++)
                    fprintf(out, "%d\n", (int)mesh->mixed_flags[(size_t)k * nxy + (size_t)j * (size_t)nx + (size_t)i]);
    }

    if (mesh->dominant_fractions) {
        fprintf(out, "SCALARS dominant_fraction double 1\n");
        fprintf(out, "LOOKUP_TABLE default\n");
        for (int k = 0; k < nz; k++)
            for (int j = 0; j < ny; j++)
                for (int i = 0; i < nx; i++)
                    fprintf(out, "%.15g\n", mesh->dominant_fractions[(size_t)k * nxy + (size_t)j * (size_t)nx + (size_t)i]);
    }

    return 0;
}

/* ============================================================================
 * Export API
 * ============================================================================ */

int alea_mesh_export_stream(const alea_mesh_result_t *mesh,
                                alea_mesh_format_t fmt, FILE *out) {
    if (!mesh || !out) return -1;
    switch (fmt) {
        case ALEA_MESH_GMSH: return write_gmsh(mesh, out);
        case ALEA_MESH_VTK:  return write_vtk(mesh, out);
        default: return -1;
    }
}

int alea_mesh_export(const alea_mesh_result_t *mesh,
                         alea_mesh_format_t fmt, const char *filename) {
    if (!mesh || !filename) return -1;
    FILE *f = fopen(filename, "w");
    if (!f) return -1;
    int rc = alea_mesh_export_stream(mesh, fmt, f);
    fclose(f);
    return rc;
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
