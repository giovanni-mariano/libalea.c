// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

#include "cluster_internal.h"

#include <limits.h>
#include <mpi.h>
#include <stdlib.h>

typedef struct {
    MPI_Comm communicator;
} alea_cluster_mpi_state_t;

static int g_owns_mpi;

int alea_cluster_backend_initialize(int* argc, char*** argv) {
    int initialized = 0;
    int finalized = 0;
    if (MPI_Initialized(&initialized) != MPI_SUCCESS ||
        MPI_Finalized(&finalized) != MPI_SUCCESS || finalized) return -1;
    if (!initialized) {
        int provided = MPI_THREAD_SINGLE;
        if (MPI_Init_thread(argc, argv, MPI_THREAD_FUNNELED, &provided) !=
                MPI_SUCCESS) return -1;
        g_owns_mpi = 1;
        if (provided < MPI_THREAD_FUNNELED) {
            MPI_Finalize();
            g_owns_mpi = 0;
            return -1;
        }
    } else {
        int provided = MPI_THREAD_SINGLE;
        if (MPI_Query_thread(&provided) != MPI_SUCCESS ||
            provided < MPI_THREAD_FUNNELED) return -1;
    }
    return 0;
}

int alea_cluster_backend_finalize(void) {
    if (!g_owns_mpi) return 0;
    g_owns_mpi = 0;
    return MPI_Finalize() == MPI_SUCCESS ? 0 : -1;
}

int alea_cluster_backend_create(void** out, int* rank, int* size) {
    if (!out || !rank || !size) return -1;
    alea_cluster_mpi_state_t* state = calloc(1, sizeof(*state));
    if (!state) return -1;
    state->communicator = MPI_COMM_NULL;
    if (MPI_Comm_dup(MPI_COMM_WORLD, &state->communicator) != MPI_SUCCESS ||
        MPI_Comm_set_errhandler(state->communicator, MPI_ERRORS_RETURN) !=
            MPI_SUCCESS ||
        MPI_Comm_rank(state->communicator, rank) != MPI_SUCCESS ||
        MPI_Comm_size(state->communicator, size) != MPI_SUCCESS) {
        if (state->communicator != MPI_COMM_NULL)
            MPI_Comm_free(&state->communicator);
        free(state);
        return -1;
    }
    *out = state;
    return 0;
}

void alea_cluster_backend_destroy(void* opaque) {
    alea_cluster_mpi_state_t* state = opaque;
    if (!state) return;
    if (state->communicator != MPI_COMM_NULL)
        MPI_Comm_free(&state->communicator);
    free(state);
}

const char* alea_cluster_backend_name(void) { return "mpi"; }

int alea_cluster_backend_agree_status(void* opaque, int local, int* global) {
    alea_cluster_mpi_state_t* state = opaque;
    return state && global && MPI_Allreduce(
        &local, global, 1, MPI_INT, MPI_MAX, state->communicator) == MPI_SUCCESS
        ? 0 : -1;
}

int alea_cluster_backend_u64_minmax(void* opaque, uint64_t local,
                                    uint64_t* minimum, uint64_t* maximum) {
    alea_cluster_mpi_state_t* state = opaque;
    if (!state || !minimum || !maximum) return -1;
    if (MPI_Allreduce(&local, minimum, 1, MPI_UINT64_T, MPI_MIN,
                      state->communicator) != MPI_SUCCESS) return -1;
    return MPI_Allreduce(&local, maximum, 1, MPI_UINT64_T, MPI_MAX,
                         state->communicator) == MPI_SUCCESS ? 0 : -1;
}

int alea_cluster_backend_sum_doubles(void* opaque, double* values,
                                     size_t count) {
    alea_cluster_mpi_state_t* state = opaque;
    if (!state || (!values && count)) return -1;
    while (count) {
        int chunk = count > (size_t)INT_MAX ? INT_MAX : (int)count;
        if (MPI_Allreduce(MPI_IN_PLACE, values, chunk, MPI_DOUBLE, MPI_SUM,
                          state->communicator) != MPI_SUCCESS) return -1;
        values += chunk;
        count -= (size_t)chunk;
    }
    return 0;
}

int alea_cluster_backend_broadcast_int(void* opaque, int* value, int root) {
    alea_cluster_mpi_state_t* state = opaque;
    return state && value && MPI_Bcast(
        value, 1, MPI_INT, root, state->communicator) == MPI_SUCCESS ? 0 : -1;
}
