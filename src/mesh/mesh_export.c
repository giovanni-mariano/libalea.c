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
#include "core/alea_spatial.h"
#include "core/alea_universe.h"
#include "primitives/bbox.h"

#include <stdlib.h>
#include <string.h>
#include <math.h>

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

/** Collect unique values from array, return sorted array and count */
static int *collect_unique(const int *arr, int count, int *out_num) {
    if (count == 0) { *out_num = 0; return NULL; }

    int *tmp = malloc((size_t)count * sizeof(int));
    if (!tmp) { *out_num = 0; return NULL; }
    memcpy(tmp, arr, (size_t)count * sizeof(int));

    /* sort */
    for (int i = 1; i < count; i++) {
        int key = tmp[i];
        int j = i - 1;
        while (j >= 0 && tmp[j] > key) { tmp[j + 1] = tmp[j]; j--; }
        tmp[j + 1] = key;
    }

    /* dedup */
    int n = 1;
    for (int i = 1; i < count; i++) {
        if (tmp[i] != tmp[i - 1])
            tmp[n++] = tmp[i];
    }

    int *result = realloc(tmp, (size_t)n * sizeof(int));
    if (!result) result = tmp;  /* realloc to smaller should not fail, but just in case */
    *out_num = n;
    return result;
}

/* ============================================================================
 * Sampling
 * ============================================================================ */

alea_mesh_result_t *alea_mesh_sample(const alea_system_t *sys,
                                             const alea_mesh_config_t *cfg) {
    if (!sys || !cfg) return NULL;
    if (cfg->nx <= 0 || cfg->ny <= 0 || cfg->nz <= 0) return NULL;

    int nx = cfg->nx, ny = cfg->ny, nz = cfg->nz;

    /* Determine bounds */
    double xlo = cfg->x_min, xhi = cfg->x_max;
    double ylo = cfg->y_min, yhi = cfg->y_max;
    double zlo = cfg->z_min, zhi = cfg->z_max;

    if (xlo == 0.0 && xhi == 0.0 && ylo == 0.0 && yhi == 0.0 &&
        zlo == 0.0 && zhi == 0.0) {
        /* Auto-detect from bounding sphere */
        double cx, cy, cz, r;
        if (alea_compute_bounding_sphere(sys, 1.0, &cx, &cy, &cz, &r) != 0)
            return NULL;
        double pad = r * cfg->auto_pad;
        xlo = cx - r - pad;  xhi = cx + r + pad;
        ylo = cy - r - pad;  yhi = cy + r + pad;
        zlo = cz - r - pad;  zhi = cz + r + pad;
    }

    /* Build node arrays */
    double *xn = cfg->x_nodes ? copy_nodes(cfg->x_nodes, nx) : build_uniform_nodes(nx, xlo, xhi);
    double *yn = cfg->y_nodes ? copy_nodes(cfg->y_nodes, ny) : build_uniform_nodes(ny, ylo, yhi);
    double *zn = cfg->z_nodes ? copy_nodes(cfg->z_nodes, nz) : build_uniform_nodes(nz, zlo, zhi);
    if (!xn || !yn || !zn) { free(xn); free(yn); free(zn); return NULL; }

    size_t ncells = (size_t)nx * (size_t)ny * (size_t)nz;
    int *mat_ids = malloc(ncells * sizeof(int));
    int *cell_ids = malloc(ncells * sizeof(int));
    if (!mat_ids || !cell_ids) {
        free(xn); free(yn); free(zn); free(mat_ids); free(cell_ids);
        return NULL;
    }

    /* Ensure indices are built (cast away const — same pattern as slice_api.c) */
    alea_system_t *sys_mut = (alea_system_t *)sys;
    if (!sys_mut->universe_index_built)
        alea_build_universe_index(sys_mut);
    if (!sys_mut->cell_adjacency_built)
        alea_build_cell_adjacency(sys_mut);
    alea_spatial_index_build(sys_mut);

    /* Sample cell centers */
    int void_mat = cfg->void_material_id;

#ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic, 4)
#endif
    for (int k = 0; k < nz; k++) {
        for (int j = 0; j < ny; j++) {
            for (int i = 0; i < nx; i++) {
                double xc = (xn[i] + xn[i + 1]) * 0.5;
                double yc = (yn[j] + yn[j + 1]) * 0.5;
                double zc = (zn[k] + zn[k + 1]) * 0.5;
                size_t idx = (size_t)k * (size_t)(ny * nx) + (size_t)j * (size_t)nx + (size_t)i;

                int cell_idx = alea_identify_cell_at_point(sys, xc, yc, zc);
                if (cell_idx < 0) {
                    mat_ids[idx] = void_mat;
                    cell_ids[idx] = -1;
                } else {
                    mat_ids[idx] = sys->cells.data[cell_idx].material_id;
                    cell_ids[idx] = sys->cells.data[cell_idx].mcnp_cell_id;
                }
            }
        }
    }

    /* Collect unique materials */
    int num_mats = 0;
    int *unique_mats = collect_unique(mat_ids, (int)ncells, &num_mats);

    /* Build result */
    alea_mesh_result_t *res = calloc(1, sizeof(*res));
    if (!res) {
        free(xn); free(yn); free(zn); free(mat_ids); free(cell_ids); free(unique_mats);
        return NULL;
    }
    res->nx = nx;  res->ny = ny;  res->nz = nz;
    res->x_nodes = xn;  res->y_nodes = yn;  res->z_nodes = zn;
    res->material_ids = mat_ids;
    res->cell_ids = cell_ids;
    res->num_materials = num_mats;
    res->unique_materials = unique_mats;
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
    free(mesh);
}

/* ============================================================================
 * Gmsh .msh v2.2 ASCII Writer
 * ============================================================================ */

/** Node ID: 1-based, (i,j,k) -> k*(ny+1)*(nx+1) + j*(nx+1) + i + 1 */
#define GMSH_NODE(i, j, k, nx1, ny1) \
    ((k) * (ny1) * (nx1) + (j) * (nx1) + (i) + 1)

static int write_gmsh(const alea_mesh_result_t *mesh, FILE *out) {
    int nx = mesh->nx, ny = mesh->ny, nz = mesh->nz;
    int nx1 = nx + 1, ny1 = ny + 1, nz1 = nz + 1;
    int nnodes = nx1 * ny1 * nz1;
    int nelems = nx * ny * nz;

    /* Header */
    fprintf(out, "$MeshFormat\n2.2 0 8\n$EndMeshFormat\n");

    /* Physical names (one per unique material) */
    fprintf(out, "$PhysicalNames\n%d\n", mesh->num_materials);
    for (int m = 0; m < mesh->num_materials; m++) {
        fprintf(out, "3 %d \"material_%d\"\n",
                mesh->unique_materials[m], mesh->unique_materials[m]);
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
                int mat = mesh->material_ids[(size_t)k * (ny * nx) + (size_t)j * nx + i];

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
                        eid++, mat, mat,
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
    int nelems = nx * ny * nz;

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
                fprintf(out, "%d\n", mesh->material_ids[(size_t)k * (ny * nx) + (size_t)j * nx + i]);

    /* Cell IDs */
    fprintf(out, "SCALARS cell_id int 1\n");
    fprintf(out, "LOOKUP_TABLE default\n");
    for (int k = 0; k < nz; k++)
        for (int j = 0; j < ny; j++)
            for (int i = 0; i < nx; i++)
                fprintf(out, "%d\n", mesh->cell_ids[(size_t)k * (ny * nx) + (size_t)j * nx + i]);

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

int alea_mesh_export_system(const alea_system_t *sys,
                             const alea_mesh_config_t *cfg,
                             const char *filename) {
    if (!sys || !cfg || !filename) return -1;

    alea_mesh_result_t *mesh = alea_mesh_sample(sys, cfg);
    if (!mesh) return -1;

    int rc = alea_mesh_export(mesh, cfg->format, filename);
    alea_mesh_result_free(mesh);
    return rc;
}
