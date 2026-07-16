// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

#ifndef ALEA_CELL_H
#define ALEA_CELL_H

/**
 * @file alea_cell.h
 * @brief Cell types and operations
 *
 * Cell entry struct, cell neighbor struct, and cell API functions.
 */

#include "alea_types.h"
#include "util/alea_vec.h"


// ============================================================================
// CELL ADJACENCY
// ============================================================================

/**
 * @brief Neighbor relationship between two cells
 *
 * When crossing the specified surface from this cell, you enter neighbor_cell.
 */
typedef struct alea_cell_neighbor {
    int surface_id;           // MCNP surface ID that separates the cells
    uint32_t surface_index;   // Index into sys->surfaces[]
    int neighbor_cell_id;     // MCNP cell ID of the neighbor
    uint32_t neighbor_index;  // Index into sys->cells.data[]
    int8_t our_sense;         // Our sense of the surface (-1 or +1)
} alea_cell_neighbor_t;

// ============================================================================
// CELL ENTRY
// ============================================================================

/**
 * @brief Cell information
 */
typedef struct alea_cell_entry {
    int mc_cell_id;
    uint32_t root_node_id;
    uint32_t original_root_node_id;  // Pre-TRCL root (ALEA_NODE_ID_INVALID if none)
    int material_id;
    int material_index;     /* Index into sys->materials (-1 for void) */
    double density;
    unsigned int is_mass_density : 1;  // 1 = g/cm³, 0 = atoms/b-cm

    // Temperature (Kelvin, 0.0 = not set)
    double temperature;
    unsigned int has_temperature : 1;

    int universe_id;        // Which universe this cell belongs to (0 = base)
    int fill_universe;      // If >0, this cell is filled with this universe
    int fill_transform;     // Transform ID to apply to filled universe

    // Lattice support
    int lat_type;               // 0=none, 1=rect, 2=hex
    int lat_fill_dims[6];       // imin, imax, jmin, jmax, kmin, kmax
    int* lat_fill;              // Array of universe IDs for lattice
    size_t lat_fill_count;      // Total count of universes in array
    double lat_pitch[3];        // Element pitch in each dimension
    double lat_lower_left[3];   // Lower-left corner of the lattice

    // Cached surface references (built by alea_build_cell_surface_index)
    uint32_t* surface_indices;      // Array of indices into sys->surfaces[]
    size_t surface_index_count;     // Number of surface indices

    // Cell adjacency (built by alea_build_cell_adjacency)
    struct alea_cell_neighbor* neighbors;  // Array of neighboring cells
    size_t neighbor_count;                // Number of neighbors

    // Comments (preserved from MCNP input)
    char* comments;        // "C" comment lines before the cell (NULL if none)
    char* inline_comment;  // Inline "$" comment from last line (NULL if none)
} alea_cell_entry_t;

ALEA_VEC_DEFINE(alea_cell_vec, alea_cell_entry_t);

/** True for cells that are never a terminal answer: fill containers and
 *  lattice cells. When a point resolver returns one, the fill chain below
 *  it is unresolved (ALEA_RESOLVE_UNDEFINED_FILL). */
static inline int alea_cell_entry_is_container(const alea_cell_entry_t* cell) {
    return cell->fill_universe > 0 ||
           (cell->lat_type != 0 && cell->lat_fill != NULL);
}

// ============================================================================
// API - CELL OPERATIONS
// ============================================================================

/**
 * @brief Find cell index by cell ID
 *
 * @param sys CSG system
 * @param cell_id Cell ID to find
 * @return Cell index (0-based) or -1 if not found
 */
int alea_find_cell_by_id(const alea_system_t* sys, int cell_id);

/**
 * @brief Add a cell to the system (programmatic use)
 *
 * For programmatic cell creation. If cell_id <= 0, auto-assigns a new ID.
 * If cell_id > 0 but already exists, logs a warning and auto-assigns a new ID.
 *
 * Use alea_add_cell_with_id() for file loading where duplicate IDs are errors.
 *
 * @param sys CSG system
 * @param cell_id Requested cell ID (0 or negative for auto-assign)
 * @param root_node Root node of CSG tree
 * @param material_index Material index from alea_add_material(), or ALEA_MATERIAL_VOID
 * @param density Density (g/cm3, negative for atom density)
 * @param universe_id Universe this cell belongs to
 * @return Result with cell index, or error code
 */
int alea_add_cell(alea_system_t* sys, int cell_id, alea_node_id_t root_node,
                 int material_index, double density, int universe_id);

/**
 * @brief Add a cell with explicit ID (for file loading)
 *
 * Used by MCNP/OpenMC parsers. Fast path without per-insertion duplicate check.
 * Call alea_validate_cell_ids() after loading all cells.
 *
 * @param sys CSG system
 * @param cell_id Cell ID
 * @param root_node Root node of CSG tree
 * @param material_index Material index from alea_add_material(), or ALEA_MATERIAL_VOID
 * @param density Density (g/cm3, negative for atom density)
 * @param universe_id Universe this cell belongs to
 * @return Cell index, or -1 on error
 */
int alea_add_cell_with_id(alea_system_t* sys, int cell_id, alea_node_id_t root_node,
                         int material_index, double density, int universe_id);

/**
 * @brief Validate that all cell IDs are unique
 *
 * Call after loading all cells from a file. Reports all duplicate IDs found.
 *
 * @param sys CSG system
 * @return 0 if all IDs unique, -1 if duplicates found
 */
int alea_validate_cell_ids(alea_system_t* sys);

/**
 * @brief Get the maximum cell ID in the system
 * @param sys CSG system
 * @return Maximum cell ID, or 0 if no cells exist
 */
int alea_max_cell_id(const alea_system_t* sys);

/* alea_cell_get is the public API (individual pointer args) — see alea.h */
/* alea_cell_get_info is the struct-based variant — see alea_public_api.c */

/**
 * @brief Set cell FILL reference
 * @param sys CSG system
 * @param cell_index Cell index
 * @param fill_universe Universe ID to fill with
 * @param fill_transform Transform ID (0 for none)
 * @return 0 on success, -1 on error
 */
int alea_set_cell_fill(alea_system_t* sys, int cell_index, int fill_universe, int fill_transform);

/**
 * @brief Get number of cells in the system
 */
size_t alea_cell_count(const alea_system_t* sys);


#endif /* ALEA_CELL_H */
