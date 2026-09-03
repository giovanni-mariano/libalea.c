// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

#include "alea_mesh.h"
#include "util/compat.h"

#include <errno.h>
#include <math.h>
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>

typedef struct {
    alea_adaptive_grid_result_t result;
    size_t capacity;
    size_t fraction_capacity;
    size_t cell_fraction_capacity;
    alea_system_t *sys;
    const alea_adaptive_grid_config_t *config;
} adaptive_builder_t;

void alea_adaptive_grid_config_init(alea_adaptive_grid_config_t *cfg) {
    if (!cfg) return;
    memset(cfg, 0, sizeof(*cfg));
    alea_mesh_config_init(&cfg->sampling);
    cfg->sampling.sampling_mode = ALEA_MESH_SAMPLE_ADAPTIVE;
    cfg->max_grid_depth = 4;
    cfg->max_cells = 1000000;
    cfg->refine_mixed = 1;
    cfg->refine_high_error = 1;
}

static int adaptive_reserve(adaptive_builder_t *builder, size_t needed) {
    if (needed <= builder->capacity) return 0;
    size_t capacity = builder->capacity ? builder->capacity : 64;
    while (capacity < needed) {
        if (capacity > SIZE_MAX / 2) {
            alea_set_error_detail(ALEA_ERR_OVERFLOW,
                                  "adaptive grid cell storage overflows");
            return -1;
        }
        capacity *= 2;
    }
    if (capacity > SIZE_MAX / sizeof(*builder->result.cells)) {
        alea_set_error_detail(ALEA_ERR_OVERFLOW,
                              "adaptive grid allocation overflows");
        return -1;
    }
    void *cells = realloc(builder->result.cells,
                          capacity * sizeof(*builder->result.cells));
    if (!cells) {
        alea_set_error_detail(ALEA_ERR_OUT_OF_MEMORY,
                              "adaptive grid allocation failed");
        return -1;
    }
    builder->result.cells = cells;
    builder->capacity = capacity;
    return 0;
}

static int adaptive_reserve_fractions(adaptive_builder_t *builder,
                                      size_t material_needed,
                                      size_t cell_needed) {
    if (material_needed > UINT32_MAX || cell_needed > UINT32_MAX) {
        alea_set_error_detail(ALEA_ERR_OVERFLOW,
                              "adaptive composition offsets exceed uint32_t");
        return -1;
    }
#define ADAPTIVE_GROW(TARGET, CAPACITY, NEEDED, TYPE, LABEL) do { \
    if ((NEEDED) > (CAPACITY)) { \
        size_t next = (CAPACITY) ? (CAPACITY) : 64; \
        while (next < (NEEDED)) { \
            if (next > SIZE_MAX / 2) { \
                alea_set_error_detail(ALEA_ERR_OVERFLOW, \
                                      "adaptive " LABEL " storage overflows"); \
                return -1; \
            } \
            next *= 2; \
        } \
        if (next > SIZE_MAX / sizeof(TYPE)) { \
            alea_set_error_detail(ALEA_ERR_OVERFLOW, \
                                  "adaptive " LABEL " allocation overflows"); \
            return -1; \
        } \
        void *grown = realloc((TARGET), next * sizeof(TYPE)); \
        if (!grown) { \
            alea_set_error_detail(ALEA_ERR_OUT_OF_MEMORY, \
                                  "adaptive " LABEL " allocation failed"); \
            return -1; \
        } \
        (TARGET) = grown; \
        (CAPACITY) = next; \
    } \
} while (0)
    ADAPTIVE_GROW(builder->result.fractions, builder->fraction_capacity,
                  material_needed, alea_mesh_material_fraction_t,
                  "material-fraction");
    ADAPTIVE_GROW(builder->result.cell_fractions,
                  builder->cell_fraction_capacity, cell_needed,
                  alea_mesh_cell_fraction_t, "cell-fraction");
#undef ADAPTIVE_GROW
    return 0;
}

static int adaptive_copy_sample(adaptive_builder_t *builder,
                                alea_adaptive_grid_cell_t *cell,
                                const alea_mesh_result_t *sample,
                                size_t voxel) {
    cell->material_id = sample->material_ids[voxel];
    cell->cell_id = sample->cell_ids[voxel];
    cell->mixed = sample->mixed_flags[voxel];
    cell->tie_flags = sample->tie_flags[voxel];
    cell->refinement_flags = sample->refinement_flags[voxel];
    cell->dominant_fraction = sample->dominant_fractions[voxel];
    cell->estimated_error = sample->estimated_errors[voxel];
    cell->sample_count = sample->sample_counts[voxel];
    const alea_mesh_fraction_span_t materials = sample->fraction_spans[voxel];
    const alea_mesh_fraction_span_t owners = sample->cell_fraction_spans[voxel];
    size_t material_needed = builder->result.fraction_count + materials.count;
    size_t cell_needed = builder->result.cell_fraction_count + owners.count;
    if (material_needed < builder->result.fraction_count ||
        cell_needed < builder->result.cell_fraction_count ||
        adaptive_reserve_fractions(builder, material_needed, cell_needed) != 0)
        return -1;
    cell->fraction_span.offset = (uint32_t)builder->result.fraction_count;
    cell->fraction_span.count = materials.count;
    memcpy(&builder->result.fractions[builder->result.fraction_count],
           &sample->fractions[materials.offset],
           (size_t)materials.count * sizeof(*sample->fractions));
    builder->result.fraction_count = material_needed;
    cell->cell_fraction_span.offset =
        (uint32_t)builder->result.cell_fraction_count;
    cell->cell_fraction_span.count = owners.count;
    memcpy(&builder->result.cell_fractions[builder->result.cell_fraction_count],
           &sample->cell_fractions[owners.offset],
           (size_t)owners.count * sizeof(*sample->cell_fractions));
    builder->result.cell_fraction_count = cell_needed;
    return 0;
}

static alea_mesh_result_t *adaptive_sample_box(adaptive_builder_t *builder,
                                                double x0, double x1,
                                                double y0, double y1,
                                                double z0, double z1) {
    alea_mesh_config_t cfg = builder->config->sampling;
    cfg.nx = cfg.ny = cfg.nz = 1;
    cfg.x_nodes = cfg.y_nodes = cfg.z_nodes = NULL;
    cfg.x_min = x0; cfg.x_max = x1;
    cfg.y_min = y0; cfg.y_max = y1;
    cfg.z_min = z0; cfg.z_max = z1;
    cfg.bounds_mode = ALEA_MESH_BOUNDS_EXPLICIT;
    cfg.fields = ALEA_MESH_FIELD_MATERIAL_ID |
                 ALEA_MESH_FIELD_CELL_ID |
                 ALEA_MESH_FIELD_MIXED_FLAG |
                 ALEA_MESH_FIELD_DOMINANT_FRACTION |
                 ALEA_MESH_FIELD_SAMPLE_COUNT |
                 ALEA_MESH_FIELD_TIE_FLAG |
                 ALEA_MESH_FIELD_ESTIMATED_ERROR |
                 ALEA_MESH_FIELD_REFINEMENT_FLAG |
                 ALEA_MESH_FIELD_SAMPLED_FRACTIONS |
                 ALEA_MESH_FIELD_CELL_FRACTIONS;
    cfg.max_total_samples = 0;
    cfg.progress = NULL;
    cfg.progress_user_data = NULL;
    cfg.visit = NULL;
    cfg.visit_user_data = NULL;
    return alea_mesh_sample(builder->sys, &cfg);
}

static int adaptive_wants_refinement(const adaptive_builder_t *builder,
                                     const alea_adaptive_grid_cell_t *cell) {
    if (builder->config->refine_mixed && cell->mixed) return 1;
    if (builder->config->refine_high_error &&
        (cell->estimated_error > builder->config->sampling.target_error ||
         (cell->refinement_flags & ALEA_MESH_REFINEMENT_LIMIT_REACHED))) return 1;
    return 0;
}

static int adaptive_refine(adaptive_builder_t *builder, size_t parent_index) {
    alea_adaptive_grid_cell_t parent = builder->result.cells[parent_index];
    if (!adaptive_wants_refinement(builder, &parent)) return 0;
    if (parent.level >= builder->config->max_grid_depth) {
        builder->result.cells[parent_index].flags |=
            ALEA_ADAPTIVE_GRID_DEPTH_LIMIT_REACHED;
        return 0;
    }
    if (builder->result.cell_count > builder->config->max_cells ||
        builder->config->max_cells - builder->result.cell_count < 8) {
        builder->result.cells[parent_index].flags |=
            ALEA_ADAPTIVE_GRID_CELL_LIMIT_REACHED;
        return 0;
    }
    if (adaptive_reserve(builder, builder->result.cell_count + 8) != 0) return -1;

    const double xm = 0.5 * (parent.x_min + parent.x_max);
    const double ym = 0.5 * (parent.y_min + parent.y_max);
    const double zm = 0.5 * (parent.z_min + parent.z_max);
    size_t children[8];
    for (int child = 0; child < 8; child++) {
        const double x0 = (child & 1) ? xm : parent.x_min;
        const double x1 = (child & 1) ? parent.x_max : xm;
        const double y0 = (child & 2) ? ym : parent.y_min;
        const double y1 = (child & 2) ? parent.y_max : ym;
        const double z0 = (child & 4) ? zm : parent.z_min;
        const double z1 = (child & 4) ? parent.z_max : zm;
        alea_mesh_result_t *sample = adaptive_sample_box(builder, x0, x1,
                                                          y0, y1, z0, z1);
        if (!sample) return -1;
        size_t index = builder->result.cell_count++;
        children[child] = index;
        alea_adaptive_grid_cell_t *cell = &builder->result.cells[index];
        memset(cell, 0, sizeof(*cell));
        cell->id = index + 1;
        cell->parent_id = parent_index + 1;
        cell->level = parent.level + 1;
        cell->is_leaf = 1;
        cell->x_min = x0; cell->x_max = x1;
        cell->y_min = y0; cell->y_max = y1;
        cell->z_min = z0; cell->z_max = z1;
        if (adaptive_copy_sample(builder, cell, sample, 0) != 0) {
            alea_mesh_result_free(sample);
            return -1;
        }
        alea_mesh_result_free(sample);
        builder->result.leaf_count++;
        if (cell->level > builder->result.max_level)
            builder->result.max_level = cell->level;
        builder->result.cells[parent_index].child_ids[child] = cell->id;
    }
    builder->result.cells[parent_index].is_leaf = 0;
    builder->result.leaf_count--;

    for (int child = 0; child < 8; child++)
        if (adaptive_refine(builder, children[child]) != 0) return -1;
    return 0;
}

alea_adaptive_grid_result_t *alea_adaptive_grid_sample(
    alea_system_t *sys, const alea_adaptive_grid_config_t *cfg) {
    if (!sys || !cfg) {
        alea_set_error_detail(ALEA_ERR_NULL_ARG,
                              "adaptive grid requires system and config");
        return NULL;
    }
    if (cfg->max_cells == 0 || cfg->max_grid_depth > 30) {
        alea_set_error_detail(ALEA_ERR_INVALID_ARG,
                              "adaptive grid limits are invalid");
        return NULL;
    }

    alea_mesh_config_t initial_cfg = cfg->sampling;
    initial_cfg.fields = ALEA_MESH_FIELD_MATERIAL_ID |
                         ALEA_MESH_FIELD_CELL_ID |
                         ALEA_MESH_FIELD_MIXED_FLAG |
                         ALEA_MESH_FIELD_DOMINANT_FRACTION |
                         ALEA_MESH_FIELD_SAMPLE_COUNT |
                         ALEA_MESH_FIELD_TIE_FLAG |
                         ALEA_MESH_FIELD_ESTIMATED_ERROR |
                         ALEA_MESH_FIELD_REFINEMENT_FLAG |
                         ALEA_MESH_FIELD_SAMPLED_FRACTIONS |
                         ALEA_MESH_FIELD_CELL_FRACTIONS;
    initial_cfg.visit = NULL;
    initial_cfg.visit_user_data = NULL;
    alea_mesh_result_t *initial = alea_mesh_sample(sys, &initial_cfg);
    if (!initial) return NULL;
    size_t root_count = (size_t)initial->nx * (size_t)initial->ny *
                        (size_t)initial->nz;
    if (root_count > cfg->max_cells) {
        alea_set_error_detail(ALEA_ERR_INVALID_ARG,
                              "adaptive grid root count exceeds max_cells");
        alea_mesh_result_free(initial);
        return NULL;
    }

    adaptive_builder_t builder = {0};
    builder.sys = sys;
    builder.config = cfg;
    builder.result.root_count = root_count;
    builder.result.leaf_count = root_count;
    if (adaptive_reserve(&builder, root_count) != 0) {
        alea_mesh_result_free(initial);
        return NULL;
    }
    builder.result.cell_count = root_count;
    size_t nxy = (size_t)initial->nx * (size_t)initial->ny;
    for (int k = 0; k < initial->nz; k++) {
        for (int j = 0; j < initial->ny; j++) {
            for (int i = 0; i < initial->nx; i++) {
                size_t index = (size_t)k * nxy + (size_t)j * initial->nx + i;
                alea_adaptive_grid_cell_t *cell = &builder.result.cells[index];
                memset(cell, 0, sizeof(*cell));
                cell->id = index + 1;
                cell->is_leaf = 1;
                cell->x_min = initial->x_nodes[i];
                cell->x_max = initial->x_nodes[i + 1];
                cell->y_min = initial->y_nodes[j];
                cell->y_max = initial->y_nodes[j + 1];
                cell->z_min = initial->z_nodes[k];
                cell->z_max = initial->z_nodes[k + 1];
                if (adaptive_copy_sample(&builder, cell, initial, index) != 0) {
                    alea_mesh_result_free(initial);
                    free(builder.result.cells);
                    free(builder.result.fractions);
                    free(builder.result.cell_fractions);
                    return NULL;
                }
            }
        }
    }
    alea_mesh_result_free(initial);

    for (size_t root = 0; root < root_count; root++) {
        if (adaptive_refine(&builder, root) != 0) {
            free(builder.result.cells);
            free(builder.result.fractions);
            free(builder.result.cell_fractions);
            return NULL;
        }
    }

    alea_adaptive_grid_result_t *result = malloc(sizeof(*result));
    if (!result) {
        alea_set_error_detail(ALEA_ERR_OUT_OF_MEMORY,
                              "adaptive grid result allocation failed");
        free(builder.result.cells);
        free(builder.result.fractions);
        free(builder.result.cell_fractions);
        return NULL;
    }
    *result = builder.result;
    return result;
}

void alea_adaptive_grid_result_free(alea_adaptive_grid_result_t *grid) {
    if (!grid) return;
    free(grid->cells);
    free(grid->fractions);
    free(grid->cell_fractions);
    free(grid);
}

static int adaptive_writef(FILE *out, const char *format, ...) {
    va_list args;
    va_start(args, format);
    int rc = vfprintf(out, format, args);
    va_end(args);
    if (rc < 0) {
        alea_set_error_detail(ALEA_ERR_FILE_WRITE,
                              "adaptive grid stream write failed: %s",
                              strerror(errno));
        return -1;
    }
    return 0;
}

#define ADAPTIVE_WRITE(...) do { if (adaptive_writef(out, __VA_ARGS__) != 0) return -1; } while (0)

static int adaptive_validate(const alea_adaptive_grid_result_t *grid) {
    if (!grid || !grid->cells || grid->cell_count == 0 || grid->leaf_count == 0) {
        alea_set_error_detail(ALEA_ERR_INVALID_ARG,
                              "adaptive grid result is empty or malformed");
        return -1;
    }
    size_t leaves = 0;
    for (size_t i = 0; i < grid->cell_count; i++) {
        const alea_adaptive_grid_cell_t *cell = &grid->cells[i];
        if (cell->id != i + 1 || !(cell->x_min < cell->x_max) ||
            !(cell->y_min < cell->y_max) || !(cell->z_min < cell->z_max)) {
            alea_set_error_detail(ALEA_ERR_INVALID_ARG,
                                  "adaptive grid cell topology is malformed");
            return -1;
        }
        if (cell->is_leaf) leaves++;
        if ((size_t)cell->fraction_span.offset + cell->fraction_span.count >
                grid->fraction_count ||
            (size_t)cell->cell_fraction_span.offset +
                cell->cell_fraction_span.count > grid->cell_fraction_count) {
            alea_set_error_detail(ALEA_ERR_INVALID_ARG,
                                  "adaptive composition spans are malformed");
            return -1;
        }
    }
    if (leaves != grid->leaf_count) {
        alea_set_error_detail(ALEA_ERR_INVALID_ARG,
                              "adaptive grid leaf count is inconsistent");
        return -1;
    }
    return 0;
}

static void adaptive_corners(const alea_adaptive_grid_cell_t *c,
                             double xyz[8][3]) {
    const double x[2] = {c->x_min, c->x_max};
    const double y[2] = {c->y_min, c->y_max};
    const double z[2] = {c->z_min, c->z_max};
    const int bits[8][3] = {{0,0,0},{1,0,0},{1,1,0},{0,1,0},
                            {0,0,1},{1,0,1},{1,1,1},{0,1,1}};
    for (int n = 0; n < 8; n++) {
        xyz[n][0] = x[bits[n][0]];
        xyz[n][1] = y[bits[n][1]];
        xyz[n][2] = z[bits[n][2]];
    }
}

static int adaptive_write_gmsh(const alea_adaptive_grid_result_t *grid,
                               FILE *out) {
    size_t nodes = grid->leaf_count * 8;
    ADAPTIVE_WRITE("$MeshFormat\n2.2 0 8\n$EndMeshFormat\n");
    ADAPTIVE_WRITE("$Nodes\n%lu\n", (unsigned long)nodes);
    size_t node = 1;
    for (size_t i = 0; i < grid->cell_count; i++) if (grid->cells[i].is_leaf) {
        double xyz[8][3]; adaptive_corners(&grid->cells[i], xyz);
        for (int n = 0; n < 8; n++, node++)
            ADAPTIVE_WRITE("%lu %.15g %.15g %.15g\n", (unsigned long)node,
                           xyz[n][0], xyz[n][1], xyz[n][2]);
    }
    ADAPTIVE_WRITE("$EndNodes\n$Elements\n%lu\n",
                   (unsigned long)grid->leaf_count);
    size_t element = 1; node = 1;
    for (size_t i = 0; i < grid->cell_count; i++) if (grid->cells[i].is_leaf) {
        ADAPTIVE_WRITE("%lu 5 0 %lu %lu %lu %lu %lu %lu %lu %lu\n",
                       (unsigned long)element++, (unsigned long)node,
                       (unsigned long)(node + 1), (unsigned long)(node + 2),
                       (unsigned long)(node + 3), (unsigned long)(node + 4),
                       (unsigned long)(node + 5), (unsigned long)(node + 6),
                       (unsigned long)(node + 7));
        node += 8;
    }
    ADAPTIVE_WRITE("$EndElements\n");
    const char *names[] = {"material_id", "cell_id", "level", "mixed_flag"};
    for (int field = 0; field < 4; field++) {
        ADAPTIVE_WRITE("$ElementData\n1\n\"%s\"\n1\n0.0\n3\n0\n1\n%lu\n",
                       names[field], (unsigned long)grid->leaf_count);
        element = 1;
        for (size_t i = 0; i < grid->cell_count; i++) if (grid->cells[i].is_leaf) {
            const alea_adaptive_grid_cell_t *c = &grid->cells[i];
            int value = field == 0 ? c->material_id : field == 1 ? c->cell_id :
                        field == 2 ? (int)c->level : (int)c->mixed;
            ADAPTIVE_WRITE("%lu %d\n", (unsigned long)element++, value);
        }
        ADAPTIVE_WRITE("$EndElementData\n");
    }
    return 0;
}

static int adaptive_write_vtk(const alea_adaptive_grid_result_t *grid,
                              FILE *out) {
    size_t points = grid->leaf_count * 8;
    ADAPTIVE_WRITE("# vtk DataFile Version 3.0\nAlea adaptive grid\nASCII\n");
    ADAPTIVE_WRITE("DATASET UNSTRUCTURED_GRID\nPOINTS %lu double\n",
                   (unsigned long)points);
    for (size_t i = 0; i < grid->cell_count; i++) if (grid->cells[i].is_leaf) {
        double xyz[8][3]; adaptive_corners(&grid->cells[i], xyz);
        for (int n = 0; n < 8; n++)
            ADAPTIVE_WRITE("%.15g %.15g %.15g\n", xyz[n][0], xyz[n][1], xyz[n][2]);
    }
    ADAPTIVE_WRITE("CELLS %lu %lu\n", (unsigned long)grid->leaf_count,
                   (unsigned long)(grid->leaf_count * 9));
    size_t point = 0;
    for (size_t i = 0; i < grid->cell_count; i++) if (grid->cells[i].is_leaf) {
        ADAPTIVE_WRITE("8 %lu %lu %lu %lu %lu %lu %lu %lu\n",
                       (unsigned long)point, (unsigned long)(point + 1),
                       (unsigned long)(point + 2), (unsigned long)(point + 3),
                       (unsigned long)(point + 4), (unsigned long)(point + 5),
                       (unsigned long)(point + 6), (unsigned long)(point + 7));
        point += 8;
    }
    ADAPTIVE_WRITE("CELL_TYPES %lu\n", (unsigned long)grid->leaf_count);
    for (size_t i = 0; i < grid->leaf_count; i++) ADAPTIVE_WRITE("12\n");
    ADAPTIVE_WRITE("CELL_DATA %lu\n", (unsigned long)grid->leaf_count);
    const char *names[] = {"material_id", "cell_id", "level", "mixed_flag"};
    for (int field = 0; field < 4; field++) {
        ADAPTIVE_WRITE("SCALARS %s int 1\nLOOKUP_TABLE default\n", names[field]);
        for (size_t i = 0; i < grid->cell_count; i++) if (grid->cells[i].is_leaf) {
            const alea_adaptive_grid_cell_t *c = &grid->cells[i];
            int value = field == 0 ? c->material_id : field == 1 ? c->cell_id :
                        field == 2 ? (int)c->level : (int)c->mixed;
            ADAPTIVE_WRITE("%d\n", value);
        }
    }
    return 0;
}

int alea_adaptive_grid_export_stream(const alea_adaptive_grid_result_t *grid,
                                     alea_mesh_format_t fmt, FILE *out) {
    if (!out) {
        alea_set_error_detail(ALEA_ERR_NULL_ARG,
                              "adaptive grid output stream is NULL");
        return -1;
    }
    if (adaptive_validate(grid) != 0) return -1;
    int rc = fmt == ALEA_MESH_GMSH ? adaptive_write_gmsh(grid, out) :
             fmt == ALEA_MESH_VTK ? adaptive_write_vtk(grid, out) : -1;
    if (fmt != ALEA_MESH_GMSH && fmt != ALEA_MESH_VTK)
        alea_set_error_detail(ALEA_ERR_INVALID_ARG,
                              "adaptive grid export format is invalid");
    if (rc == 0 && (fflush(out) != 0 || ferror(out))) {
        alea_set_error_detail(ALEA_ERR_FILE_WRITE,
                              "adaptive grid stream flush failed");
        rc = -1;
    }
    return rc;
}

int alea_adaptive_grid_export(const alea_adaptive_grid_result_t *grid,
                              alea_mesh_format_t fmt, const char *filename) {
    if (!filename) {
        alea_set_error_detail(ALEA_ERR_NULL_ARG,
                              "adaptive grid filename is NULL");
        return -1;
    }
    char temporary[4096];
    FILE *out = alea_sibling_tmpfile(filename, temporary, sizeof(temporary));
    if (!out) {
        alea_set_error_detail(ALEA_ERR_FILE_WRITE,
                              "cannot create adaptive grid temporary file");
        return -1;
    }
    int rc = alea_adaptive_grid_export_stream(grid, fmt, out);
    if (fclose(out) != 0 && rc == 0) {
        alea_set_error_detail(ALEA_ERR_FILE_WRITE,
                              "cannot close adaptive grid temporary file");
        rc = -1;
    }
    if (rc == 0 && alea_replace_file(temporary, filename) != 0) {
        alea_set_error_detail(ALEA_ERR_FILE_WRITE,
                              "cannot replace adaptive grid output file");
        rc = -1;
    }
    if (rc != 0) remove(temporary);
    return rc;
}
