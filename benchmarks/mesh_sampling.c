// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

#include "alea.h"
#include "alea_mesh.h"

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

static alea_system_t *make_scene(int kind) {
    alea_system_t *sys = alea_create();
    if (!sys) return NULL;
    if (kind == 0) {
        int box = alea_box_surface(sys, 1, -10, 10, -10, 10, -10, 10);
        int material = alea_add_material(sys, 1);
        alea_add_cell(sys, 1, alea_halfspace(sys, box, -1), material, -1, 0);
    } else if (kind == 1) {
        int plane = alea_plane_surface(sys, 1, 1, 0, 0, 0);
        int left = alea_add_material(sys, 1);
        int right = alea_add_material(sys, 2);
        alea_add_cell(sys, 1, alea_halfspace(sys, plane, -1), left, -1, 0);
        alea_add_cell(sys, 2, alea_halfspace(sys, plane, 1), right, -1, 0);
    } else {
        for (int i = 0; i < 16; i++) {
            double x0 = -10.0 + 1.25 * i;
            double x1 = x0 + 1.25;
            int box = alea_box_surface(sys, i + 1, x0, x1,
                                       -10, 10, -10, 10);
            int material = alea_add_material(sys, i + 1);
            alea_add_cell(sys, i + 1, alea_halfspace(sys, box, -1),
                          material, -1, 0);
        }
    }
    if (alea_build_universe_index(sys) != 0) {
        alea_destroy(sys);
        return NULL;
    }
    return sys;
}

static size_t retained_bytes(const alea_mesh_result_t *mesh) {
    size_t voxels = (size_t)mesh->nx * mesh->ny * mesh->nz;
    size_t bytes = sizeof(*mesh) +
        ((size_t)mesh->nx + mesh->ny + mesh->nz + 3) * sizeof(double) +
        (size_t)mesh->num_materials * sizeof(int) +
        mesh->fraction_count * sizeof(*mesh->fractions);
    if (mesh->material_ids) bytes += voxels * sizeof(int);
    if (mesh->cell_ids) bytes += voxels * sizeof(int);
    if (mesh->mixed_flags) bytes += voxels;
    if (mesh->dominant_fractions) bytes += voxels * sizeof(double);
    if (mesh->estimated_errors) bytes += voxels * sizeof(double);
    if (mesh->sample_counts) bytes += voxels * sizeof(uint32_t);
    if (mesh->tie_flags) bytes += voxels;
    if (mesh->refinement_flags) bytes += voxels;
    if (mesh->fraction_spans) bytes += voxels * sizeof(*mesh->fraction_spans);
    return bytes;
}

int main(int argc, char **argv) {
    int n = argc > 1 ? atoi(argv[1]) : 64;
    if (n < 1) n = 1;
    const char *names[] = {"uniform", "interface", "many-material"};
    const alea_mesh_sampling_mode_t modes[] = {
        ALEA_MESH_SAMPLE_CENTER, ALEA_MESH_SAMPLE_SUBCELL,
        ALEA_MESH_SAMPLE_STRATIFIED
    };
    for (int scene = 0; scene < 3; scene++) {
        alea_system_t *sys = make_scene(scene);
        if (!sys) { fprintf(stderr, "%s\n", alea_error()); return 1; }
        for (int mode = 0; mode < 3; mode++) {
            alea_mesh_config_t cfg;
            alea_mesh_config_init(&cfg);
            cfg.nx = cfg.ny = cfg.nz = n;
            cfg.x_min = cfg.y_min = cfg.z_min = -10;
            cfg.x_max = cfg.y_max = cfg.z_max = 10;
            cfg.bounds_mode = ALEA_MESH_BOUNDS_EXPLICIT;
            cfg.sampling_mode = modes[mode];
            cfg.fields &= ~ALEA_MESH_FIELD_SAMPLED_FRACTIONS;
            clock_t begin = clock();
            alea_mesh_result_t *mesh = alea_mesh_sample(sys, &cfg);
            double seconds = (double)(clock() - begin) / CLOCKS_PER_SEC;
            if (!mesh) { fprintf(stderr, "%s\n", alea_error()); return 1; }
            size_t voxels = (size_t)n * n * n;
            printf("%-13s mode=%d voxels=%lu seconds=%.6f ns/voxel=%.1f bytes/voxel=%.1f\n",
                   names[scene], (int)modes[mode], (unsigned long)voxels,
                   seconds, seconds * 1e9 / voxels,
                   (double)retained_bytes(mesh) / voxels);
            alea_mesh_result_free(mesh);
        }
        alea_destroy(sys);
    }

    alea_system_t *sys = make_scene(1);
    if (!sys) return 1;
    const uint32_t masks[] = {
        ALEA_MESH_FIELD_MATERIAL_ID,
        ALEA_MESH_FIELD_MATERIAL_ID | ALEA_MESH_FIELD_CELL_ID |
            ALEA_MESH_FIELD_MIXED_FLAG | ALEA_MESH_FIELD_DOMINANT_FRACTION,
        UINT32_C(0x1ff)
    };
    const char *mask_names[] = {"material-only", "diagnostics", "complete"};
    for (int m = 0; m < 3; m++) {
        alea_mesh_config_t cfg;
        alea_mesh_config_init(&cfg);
        cfg.nx = cfg.ny = cfg.nz = n;
        cfg.x_min = cfg.y_min = cfg.z_min = -10;
        cfg.x_max = cfg.y_max = cfg.z_max = 10;
        cfg.bounds_mode = ALEA_MESH_BOUNDS_EXPLICIT;
        cfg.fields = masks[m];
        alea_mesh_result_t *mesh = alea_mesh_sample(sys, &cfg);
        if (!mesh) { fprintf(stderr, "%s\n", alea_error()); return 1; }
        size_t voxels = (size_t)n * n * n;
        printf("memory %-13s retained=%lu bytes bytes/voxel=%.1f\n",
               mask_names[m], (unsigned long)retained_bytes(mesh),
               (double)retained_bytes(mesh) / voxels);
        alea_mesh_result_free(mesh);
    }
    alea_destroy(sys);
    return 0;
}
