// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

#ifndef ALEA_CLUSTER_MPI_H
#define ALEA_CLUSTER_MPI_H

#include <mpi.h>
#include "alea_cluster_base.h"

#ifdef __cplusplus
extern "C" {
#endif

/** Collectively create a context on the supplied intracommunicator. Libalea
 * duplicates it and frees only the duplicate in alea_cluster_destroy().
 * The caller retains and must eventually free its communicator. Every
 * member must create, call operations, and destroy in matching order within
 * this group; ranks outside the group do not participate. A process may
 * have one live cluster context. MPI must first be initialized through
 * alea_cluster_initialize() or by the application with MPI_THREAD_FUNNELED.
 * Returns NULL on failure or for MPI_COMM_NULL. */
alea_cluster_t* alea_cluster_create_mpi(MPI_Comm communicator);

#ifdef __cplusplus
}
#endif

#endif /* ALEA_CLUSTER_MPI_H */
