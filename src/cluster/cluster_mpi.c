// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

#include "cluster_internal.h"
#include "alea_cluster_mpi.h"

#include <limits.h>
#include <mpi.h>
#include <stdlib.h>
#include <string.h>

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

static int create_mpi_backend(void** out, int* rank, int* size,
                              MPI_Comm communicator, int ready) {
    if (!out || !rank || !size || communicator == MPI_COMM_NULL) return -1;
    int intercommunicator = 0;
    if (MPI_Comm_test_inter(communicator, &intercommunicator) != MPI_SUCCESS ||
        intercommunicator) return -1;
    *out = NULL;
    int local_ready = ready != 0, all_ready = 0;
    if (MPI_Allreduce(&local_ready, &all_ready, 1, MPI_INT, MPI_MIN,
                      communicator) != MPI_SUCCESS || !all_ready)
        return -1;
    alea_cluster_mpi_state_t* state = calloc(1, sizeof(*state));
    local_ready = state != NULL;
    if (MPI_Allreduce(&local_ready, &all_ready, 1, MPI_INT, MPI_MIN,
                      communicator) != MPI_SUCCESS || !all_ready) {
        free(state);
        return -1;
    }
    state->communicator = MPI_COMM_NULL;
    if (MPI_Comm_dup(communicator, &state->communicator) != MPI_SUCCESS ||
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

int alea_cluster_backend_create(void** out, int* rank, int* size, int ready) {
    return create_mpi_backend(out, rank, size, MPI_COMM_WORLD, ready);
}

static int create_supplied_backend(void** out, int* rank, int* size,
                                   void* user_data, int ready) {
    return create_mpi_backend(out, rank, size,
        *(const MPI_Comm*)user_data, ready);
}

alea_cluster_t* alea_cluster_create_mpi(MPI_Comm communicator) {
    if (communicator == MPI_COMM_NULL) return NULL;
    return alea_cluster_create_with_backend(create_supplied_backend,
        &communicator);
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

int alea_cluster_backend_broadcast_u64(void* opaque, uint64_t* value, int root) {
    alea_cluster_mpi_state_t* state = opaque;
    return state && value && MPI_Bcast(
        value, 1, MPI_UINT64_T, root, state->communicator) == MPI_SUCCESS
        ? 0 : -1;
}

int alea_cluster_backend_broadcast_bytes(void* opaque, void* bytes,
                                         size_t count, int root) {
    alea_cluster_mpi_state_t* state = opaque;
    if (!state || (!bytes && count)) return -1;
    unsigned char* cursor = bytes;
    while (count) {
        int chunk = count > (size_t)INT_MAX ? INT_MAX : (int)count;
        if (MPI_Bcast(cursor, chunk, MPI_BYTE, root, state->communicator) !=
            MPI_SUCCESS) return -1;
        cursor += chunk;
        count -= (size_t)chunk;
    }
    return 0;
}

static size_t rank_begin(size_t count, size_t rank, size_t ranks) {
    const size_t base = count / ranks;
    const size_t extra = count % ranks;
    return rank * base + (rank < extra ? rank : extra);
}

static size_t rank_count(size_t count, size_t rank, size_t ranks) {
    return count / ranks + (rank < count % ranks ? 1u : 0u);
}

int alea_cluster_backend_scatter_blocks(void* opaque, const void* root_data,
                                        void* local_data, size_t item_size,
                                        size_t item_count, int root) {
    alea_cluster_mpi_state_t* state = opaque;
    int rank, ranks;
    if (!state || !item_size ||
        MPI_Comm_rank(state->communicator, &rank) != MPI_SUCCESS ||
        MPI_Comm_size(state->communicator, &ranks) != MPI_SUCCESS ||
        ranks <= 0) return -1;
    const size_t local_count = rank_count(item_count, (size_t)rank,
                                          (size_t)ranks);
    if (local_count && !local_data) return -1;
    if (rank == root) {
        if (item_count && !root_data) return -1;
        const unsigned char* source = root_data;
        for (int peer = 0; peer < ranks; ++peer) {
            size_t bytes = rank_count(item_count, (size_t)peer,
                                      (size_t)ranks) * item_size;
            const unsigned char* cursor = source +
                rank_begin(item_count, (size_t)peer, (size_t)ranks) * item_size;
            if (peer == root) {
                if (bytes) memcpy(local_data, cursor, bytes);
                continue;
            }
            while (bytes) {
                int chunk = bytes > (size_t)INT_MAX ? INT_MAX : (int)bytes;
                if (MPI_Send(cursor, chunk, MPI_BYTE, peer, 31,
                             state->communicator) != MPI_SUCCESS) return -1;
                cursor += chunk;
                bytes -= (size_t)chunk;
            }
        }
    } else {
        size_t bytes = local_count * item_size;
        unsigned char* cursor = local_data;
        while (bytes) {
            int chunk = bytes > (size_t)INT_MAX ? INT_MAX : (int)bytes;
            if (MPI_Recv(cursor, chunk, MPI_BYTE, root, 31,
                         state->communicator, MPI_STATUS_IGNORE) != MPI_SUCCESS)
                return -1;
            cursor += chunk;
            bytes -= (size_t)chunk;
        }
    }
    return 0;
}

int alea_cluster_backend_gather_u64(void* opaque, uint64_t local,
                                    uint64_t* root_values, int root) {
    alea_cluster_mpi_state_t* state = opaque;
    return state && MPI_Gather(&local, 1, MPI_UINT64_T, root_values, 1,
                               MPI_UINT64_T, root, state->communicator) ==
                               MPI_SUCCESS ? 0 : -1;
}

int alea_cluster_backend_gather_bytes(void* opaque, const void* local_data,
                                      size_t local_bytes, void* root_data,
                                      const size_t* root_bytes_by_rank,
                                      int root) {
    alea_cluster_mpi_state_t* state = opaque;
    int rank, ranks;
    if (!state || (local_bytes && !local_data) ||
        MPI_Comm_rank(state->communicator, &rank) != MPI_SUCCESS ||
        MPI_Comm_size(state->communicator, &ranks) != MPI_SUCCESS) return -1;
    if (rank == root) {
        if (!root_bytes_by_rank || !root_data ||
            root_bytes_by_rank[root] != local_bytes) return -1;
        unsigned char* cursor = root_data;
        for (int peer = 0; peer < ranks; ++peer) {
            size_t bytes = root_bytes_by_rank[peer];
            if (peer == root) {
                if (bytes) memcpy(cursor, local_data, bytes);
                cursor += bytes;
                continue;
            }
            while (bytes) {
                int chunk = bytes > (size_t)INT_MAX ? INT_MAX : (int)bytes;
                if (MPI_Recv(cursor, chunk, MPI_BYTE, peer, 32,
                             state->communicator, MPI_STATUS_IGNORE) !=
                    MPI_SUCCESS) return -1;
                cursor += chunk;
                bytes -= (size_t)chunk;
            }
        }
    } else {
        const unsigned char* cursor = local_data;
        while (local_bytes) {
            int chunk = local_bytes > (size_t)INT_MAX
                ? INT_MAX : (int)local_bytes;
            if (MPI_Send(cursor, chunk, MPI_BYTE, root, 32,
                         state->communicator) != MPI_SUCCESS) return -1;
            cursor += chunk;
            local_bytes -= (size_t)chunk;
        }
    }
    return 0;
}
