// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

#include "alea.h"
#include "alea_cluster.h"

#include <math.h>
#include <stdio.h>

static int fail(int rank, const char* message) {
    fprintf(stderr, "rank %d: %s\n", rank, message);
    return 1;
}

int main(int argc, char** argv) {
    alea_cluster_status_t status = alea_cluster_initialize(&argc, &argv);
    if (status != ALEA_CLUSTER_OK) return fail(-1, "cluster initialization failed");
    alea_cluster_t* cluster = alea_cluster_create();
    if (!cluster) return fail(-1, "cluster context creation failed");
    const int rank = alea_cluster_rank(cluster);

    alea_system_t* sys = alea_create();
    if (!sys) return fail(rank, "system creation failed");
    int material = alea_add_material(sys, 1);
    int surface = alea_sphere_surface(sys, 1, 0.0, 0.0, 0.0, 1.0);
    alea_node_id_t interior = alea_halfspace(sys, surface, -1);
    if (material < 0 || surface < 0 ||
        alea_add_cell(sys, 1, interior, material, 1.0, 0) < 0)
        return fail(rank, "sphere model creation failed");

    alea_volume_estimate_options_t options;
    alea_volume_estimate_options_init(&options);
    options.max_rays = 1003;
    options.batch_size = 211;
    options.seed = UINT64_C(987654321);
    options.requested_workers = 1;
    options.use_sampling_sphere = true;
    options.sampling_radius = 1.1;

    double volume = 0.0, error = 0.0;
    alea_cluster_volume_stats_t stats;
    status = alea_cluster_estimate_volumes(
        cluster, sys, &options, &volume, &error, &stats);
    if (status != ALEA_CLUSTER_OK)
        return fail(rank, alea_cluster_status_string(status));
    const double exact = 4.0 * 3.14159265358979323846 / 3.0;
    if (!(volume > exact * 0.7 && volume < exact * 1.3))
        return fail(rank, "volume estimate outside tolerance");
    if (!(error >= 0.0) || stats.volume.rays_completed != options.max_rays ||
        stats.rank_count != alea_cluster_size(cluster))
        return fail(rank, "invalid execution statistics");

    if (alea_cluster_is_root(cluster))
        printf("cluster=%s ranks=%d volume=%.8g rel_error=%.5g\n",
               alea_cluster_backend(cluster), alea_cluster_size(cluster),
               volume, error);

    if (alea_cluster_size(cluster) > 1) {
        alea_system_t* mismatch = alea_create();
        if (!mismatch) return fail(rank, "mismatch system creation failed");
        int mismatch_material = alea_add_material(mismatch, 1);
        int mismatch_surface = alea_sphere_surface(
            mismatch, 1, 0.0, 0.0, 0.0, 1.0);
        alea_node_id_t mismatch_interior =
            alea_halfspace(mismatch, mismatch_surface, -1);
        if (mismatch_material < 0 || mismatch_surface < 0 ||
            alea_add_cell(mismatch, rank == 0 ? 1 : 2, mismatch_interior,
                          mismatch_material,
                          1.0, 0) < 0)
            return fail(rank, "mismatch model creation failed");
        status = alea_cluster_estimate_volumes(
            cluster, mismatch, &options, &volume, &error, NULL);
        alea_destroy(mismatch);
        if (status != ALEA_CLUSTER_MODEL_MISMATCH)
            return fail(rank, "model mismatch was not detected");
    }

    alea_destroy(sys);
    alea_cluster_destroy(cluster);
    status = alea_cluster_finalize();
    return status == ALEA_CLUSTER_OK ? 0 : fail(rank, "cluster finalize failed");
}
