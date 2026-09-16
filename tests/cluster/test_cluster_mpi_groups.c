// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

#include "alea_cluster_mpi.h"
#include "alea_cluster_raycast.h"
#include "alea_cluster_volume.h"

#include <mpi.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static void require(int condition, const char* message) {
    if (condition) return;
    int world_rank = -1;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    fprintf(stderr, "world rank %d: %s\n", world_rank, message);
    MPI_Abort(MPI_COMM_WORLD, 1);
}

int main(int argc, char** argv) {
    require(alea_cluster_initialize(&argc, &argv) == ALEA_CLUSTER_OK,
        "cluster initialization failed");
    int world_rank = -1, world_size = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    require(world_size == 4, "test needs four MPI ranks");
    const int color = world_rank % 2;
    MPI_Comm caller_group = MPI_COMM_NULL;
    require(MPI_Comm_split(MPI_COMM_WORLD, color, world_rank,
                          &caller_group) == MPI_SUCCESS,
        "communicator split failed");
    int group_rank = -1, group_size = 0;
    MPI_Comm_rank(caller_group, &group_rank);
    MPI_Comm_size(caller_group, &group_size);
    require(group_size == 2, "wrong split-group size");

    alea_cluster_t* single = group_rank == 0
        ? alea_cluster_create_mpi(MPI_COMM_SELF) : NULL;
    require(group_rank != 0 || single != NULL,
        "single-rank context creation failed");
    alea_cluster_t* rejected = alea_cluster_create_mpi(caller_group);
    require(rejected == NULL,
        "asymmetric live-context state was not agreed by group");
    if (group_rank == 0) alea_cluster_destroy(single);

    alea_cluster_t* cluster = alea_cluster_create_mpi(caller_group);
    require(cluster && alea_cluster_rank(cluster) == group_rank &&
            alea_cluster_size(cluster) == group_size,
        "cluster context did not use caller's group");

    char* input = NULL;
    size_t length = 0;
    alea_cluster_status_t status = color == 0
        ? alea_cluster_read_file(cluster,
            group_rank == 0 ? "README.md" : NULL, &input, &length)
        : alea_cluster_read_mcnp_input(cluster,
            group_rank == 0 ? "tests/cluster/mcnp_deps_main.inp" : NULL,
            &input, &length);
    require(status == ALEA_CLUSTER_OK && input && length > 0 &&
            (color == 0 ? strstr(input, "# libalea.c") != NULL
                        : strstr(input, "1 so 1\n") != NULL),
        "group-specific collective input failed");
    free(input);

    alea_system_t* sys = alea_create();
    require(sys != NULL, "model creation failed");
    int material = alea_add_material(sys, 1);
    int surface = alea_sphere_surface(sys, 1, 0, 0, 0,
                                     color == 0 ? 1.0 : 2.0);
    require(material >= 0 && surface >= 0 &&
            alea_add_cell(sys, 7, alea_halfspace(sys, surface, -1),
                          material, 1.0, 0) >= 0,
        "group model setup failed");
    if (color == 0) {
        alea_volume_estimate_options_t options;
        alea_volume_estimate_options_init(&options);
        options.max_rays = 1003;
        options.batch_size = 127;
        options.use_sampling_sphere = true;
        options.sampling_radius = 1.1;
        double volume = 0, error = 0;
        status = alea_cluster_estimate_volumes(cluster, sys, &options,
                                               &volume, &error, NULL);
        require(status == ALEA_CLUSTER_OK && isfinite(volume) &&
                volume > 0 && isfinite(error),
            "group volume estimate failed");
    } else {
        const double origin[3] = {0, 0, -3};
        const double direction[3] = {0, 0, 1};
        unsigned char hit = 0;
        int32_t cell_id = -1;
        double enter = 0, exit = 0;
        status = alea_cluster_raycast_first_segments(cluster, sys,
            group_rank == 0 ? origin : NULL,
            group_rank == 0 ? direction : NULL,
            group_rank == 0 ? 1 : 0, group_rank == 0 ? 6.0 : 0.0,
            group_rank == 0 ? &hit : NULL,
            group_rank == 0 ? &cell_id : NULL,
            group_rank == 0 ? &enter : NULL,
            group_rank == 0 ? &exit : NULL);
        require(status == ALEA_CLUSTER_OK &&
                (group_rank != 0 || (hit && cell_id == 7 &&
                                     fabs(enter - 1.0) < 1e-9 &&
                                     fabs(exit - 5.0) < 1e-9)),
            "group raycast failed");
    }
    alea_destroy(sys);
    alea_cluster_destroy(cluster);

    int owned_rank = -1, rank_sum = -1;
    require(MPI_Comm_rank(caller_group, &owned_rank) == MPI_SUCCESS &&
            MPI_Allreduce(&owned_rank, &rank_sum, 1, MPI_INT, MPI_SUM,
                          caller_group) == MPI_SUCCESS &&
            owned_rank == group_rank && rank_sum == 1,
        "caller communicator was invalidated by cluster destruction");
    require(alea_cluster_create_mpi(MPI_COMM_NULL) == NULL,
        "MPI_COMM_NULL was accepted");
    require(MPI_Comm_free(&caller_group) == MPI_SUCCESS,
        "caller communicator cleanup failed");
    require(alea_cluster_finalize() == ALEA_CLUSTER_OK,
        "cluster finalization failed");
    if (world_rank == 0) puts("independent MPI groups passed");
    return 0;
}
