// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

#ifndef ALEA_SIMPLIFY_H
#define ALEA_SIMPLIFY_H

/**
 * @file alea_simplify.h
 * @brief CSG tree simplification and normalization
 * 
 * Location: src/core/alea_simplify.h
 * 
 * Simplifies CSG trees by:
 * 1. Converting to Negation Normal Form (pushing complements to leaves)
 * 2. Eliminating double negations
 * 3. Applying idempotent rules (A ∩ A → A, A ∪ A → A)
 * 4. Applying absorption rules (A ∩ (A ∪ B) → A, A ∪ (A ∩ B) → A)
 * 5. Hash-based subtree deduplication
 * 
 * Each cell is simplified independently - no cross-cell node sharing.
 */

#include "alea_types.h"
#include <stddef.h>
#include <stdbool.h>


typedef struct alea_system alea_system_t;

/**
 * @brief Statistics from simplification pass
 * (also defined in include/alea.h for public API)
 */
#ifndef ALEA_SIMPLIFY_STATS_DEFINED
#define ALEA_SIMPLIFY_STATS_DEFINED
typedef struct {
    size_t nodes_before;
    size_t nodes_after;
    size_t complements_eliminated;
    size_t double_negations;
    size_t idempotent_reductions;
    size_t absorption_reductions;
    size_t subtrees_deduplicated;
    size_t cell_complements_expanded;
    size_t contradictions_found;      // A ∩ ¬A = ∅
    size_t tautologies_found;         // A ∪ ¬A = U
    size_t empty_cells_removed;
    size_t union_branches_absorbed;   /* Phase 1: algebraic absorption */
    size_t union_common_factors;      /* Phase 2: common factor extraction */
    size_t union_branches_subsumed;   /* Phase 3: geometric subsumption */
} alea_simplify_stats_t;
#endif

/**
 * @brief Simplify a CSG tree
 * 
 * Applies all simplification rules to the tree rooted at root_id.
 * Creates new nodes as needed; does not modify existing nodes.
 * 
 * @param sys CSG system
 * @param root_id Root node of tree to simplify
 * @param stats Optional stats output (can be NULL)
 * @return Root node ID of simplified tree
 */
alea_node_id_t alea_tree_simplify(
    alea_system_t* sys,
    alea_node_id_t root_id,
    alea_simplify_stats_t* stats
);

/**
 * @brief Simplify all cells in the system
 * 
 * Iterates through all cells and simplifies each tree.
 * Updates cell root pointers to point to simplified trees.
 * 
 * @param sys CSG system
 * @param stats Optional stats output (can be NULL) - accumulated
 */
void alea_simplify_all_cells(
    alea_system_t* sys,
    alea_simplify_stats_t* stats
);

/**
 * @brief Convert tree to Negation Normal Form only
 * 
 * Pushes all complements to leaves without other simplifications.
 * Useful if you only want NNF without full simplification.
 * 
 * @param sys CSG system
 * @param root_id Root node of tree
 * @return Root node ID of NNF tree
 */
alea_node_id_t alea_tree_to_nnf(
    alea_system_t* sys,
    alea_node_id_t root_id
);

/**
 * @brief Check if tree is in Negation Normal Form
 * 
 * A tree is in NNF if it contains no COMPLEMENT nodes -
 * all negations are pushed to primitive sense values.
 * 
 * @param sys CSG system
 * @param root_id Root node to check
 * @return true if tree is in NNF
 */
bool alea_tree_is_nnf(
    const alea_system_t* sys,
    alea_node_id_t root_id
);

/**
 * @brief Print simplification statistics
 *
 * @param stats Statistics to print
 */
void alea_simplify_stats_print(const alea_simplify_stats_t* stats);

/**
 * @brief Split cells with top-level unions into multiple simpler cells.
 *
 * For each cell whose root is a union (T1 ∪ T2 ∪ ... ∪ Tk),
 * creates k new cells (one per branch), inheriting material, density,
 * and universe from the original. The original cell is removed.
 *
 * @param sys CSG system
 * @return Number of new cells created, or -1 on error.
 */
int alea_split_union_cells(alea_system_t* sys);

/**
 * @brief Apply De Morgan's law to create complement without COMPLEMENT nodes
 *
 * Instead of wrapping in a COMPLEMENT node, this recursively transforms:
 *   #(A ∩ B) → #A ∪ #B
 *   #(A ∪ B) → #A ∩ #B
 *   #(primitive) → primitive with flipped sense
 *
 * Useful for macrobody expansion where we want the exterior representation
 * as a union of positive-sense half-spaces rather than #(intersection).
 *
 * @param sys CSG system
 * @param node_id Node to complement
 * @return New node ID representing the complement
 */
alea_node_id_t alea_apply_demorgan(alea_system_t* sys, alea_node_id_t node_id);

/**
 * @brief Simplify all cells and remove cells proven empty
 * 
 * Applies full optimization to all cells, including:
 * - NNF conversion
 * - Associative Boolean normalization and balancing
 * - Contradiction detection
 * - Empty cell removal
 * 
 * Cells that simplify to empty (geometric impossibility) are marked
 * with root_node_id = ALEA_NODE_ID_INVALID.
 * 
 * @param sys CSG system
 * @param stats Optional accumulated stats output (can be NULL)
 */
void alea_simplify_and_prune_cells(
    alea_system_t* sys,
    alea_simplify_stats_t* stats
);


#endif // ALEA_SIMPLIFY_H
