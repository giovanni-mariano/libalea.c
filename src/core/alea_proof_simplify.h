// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

#ifndef ALEA_PROOF_SIMPLIFY_H
#define ALEA_PROOF_SIMPLIFY_H

#include "alea.h"

/* Internal bridge used by the symbolic simplifier. The branch and every
 * remaining term are existing nodes; no temporary CSG nodes are published. */
alea_proof_status_t alea_prove_union_branch_redundant(
    const alea_system_t* sys,
    alea_node_id_t branch,
    const alea_node_id_t* remaining,
    size_t remaining_count,
    const alea_bbox_t* complete_branch_bounds,
    int max_depth,
    size_t max_nodes);

#endif
