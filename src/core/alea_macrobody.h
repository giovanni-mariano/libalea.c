// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

#ifndef ALEA_MACROBODY_H
#define ALEA_MACROBODY_H

#include "core/alea_system.h"


/**
 * @file alea_macrobody.h
 * @brief Macrobody expansion - convert macrobodies to primitive surface CSG trees
 *
 * MCNP macrobodies (RCC, BOX, TRC, etc.) are single surfaces that represent
 * bounded volumes. This module provides functions to expand them into equivalent
 * CSG trees built from primitive surfaces (planes, cylinders, cones).
 *
 * Expansion is useful for:
 * - Export to codes that don't support macrobodies
 * - CSG simplification that works at the primitive level
 * - Visualization that needs to render individual faces
 */

/**
 * @brief Check if a primitive type is a macrobody
 *
 * @param type Primitive type
 * @return true if the primitive is a macrobody (RCC, BOX, TRC, etc.)
 */
bool alea_is_macrobody(alea_primitive_type_t type);

/**
 * @brief Check if a node is a 1-sheet cone (needs expansion)
 *
 * 1-sheet cones have sheet_selection != 0. They need to be expanded
 * to double-sheet cone + plane for some export formats.
 *
 * @param sys CSG system
 * @param node_id Node ID to check
 * @return true if node is a 1-sheet cone
 */
bool alea_is_1sheet_cone(const alea_system_t* sys, alea_node_id_t node_id);

/**
 * @brief Expand a macrobody node into primitive surfaces
 *
 * Creates a CSG intersection tree equivalent to the macrobody:
 * - RCC -> cylinder ∩ plane(base) ∩ plane(top)
 * - BOX -> 6 planes
 * - TRC -> cone ∩ plane(base) ∩ plane(top)
 * - SPH -> sphere (no expansion, returns same)
 * - WED -> 5 planes
 * - RHP -> 8 planes
 * - ELL -> quadric (no expansion, returns same)
 * - REC -> quadric ∩ 2 planes
 * - ARB -> N planes
 *
 * The original node is NOT modified. The caller can replace references to
 * the original node with the returned node.
 *
 * @param sys CSG system
 * @param node_id Node ID of macrobody primitive
 * @return Root node ID of expanded CSG tree, or ALEA_NODE_ID_INVALID on error
 */
alea_node_id_t alea_expand_macrobody(alea_system_t* sys, alea_node_id_t node_id);

/**
 * @brief Expand all macrobodies in a CSG tree
 *
 * Recursively walks the tree and expands any macrobody primitives.
 * Returns a new tree with all macrobodies replaced by primitive surfaces.
 *
 * @param sys CSG system
 * @param root_id Root node of CSG tree
 * @return New root node ID with expanded macrobodies, or ALEA_NODE_ID_INVALID on error
 */
alea_node_id_t alea_expand_all_macrobodies(alea_system_t* sys, alea_node_id_t root_id);

/**
 * @brief Expand a macrobody at creation time (immediate expansion)
 *
 * Called from surface conversion when a macrobody is parsed. Creates the
 * expanded CSG tree (component primitives) immediately and returns both
 * sense nodes.
 *
 * The component primitives do NOT get MCNP surface IDs assigned here -
 * that happens at export time.
 *
 * The pos_node is created using De Morgan transformation:
 *   neg_node = -cyl ∩ -plane_base ∩ -plane_top  (interior)
 *   pos_node = +cyl ∪ +plane_base ∪ +plane_top  (exterior, no COMPLEMENT)
 *
 * @param sys CSG system
 * @param type Macrobody primitive type (must be a macrobody)
 * @param data Macrobody primitive data
 * @param out_neg_node Output: negative sense node (interior)
 * @param out_pos_node Output: positive sense node (exterior)
 * @return 0 on success, -1 on error (not a macrobody or expansion failed)
 */
int alea_expand_macrobody_immediate(alea_system_t* sys,
                                    alea_primitive_type_t type,
                                    const alea_primitive_data_t* data,
                                    alea_node_id_t* out_neg_node,
                                    alea_node_id_t* out_pos_node);

/**
 * @brief Expand a 1-sheet cone at creation time (immediate expansion)
 *
 * Called from surface conversion when a 1-sheet cone is parsed.
 * A 1-sheet cone has sheet_selection != 0 and is expanded to:
 * - A double-sheet cone (sheet_selection = 0)
 * - A plane at the apex that selects the appropriate sheet
 *
 * @param sys CSG system
 * @param type Cone primitive type (CONE_X, CONE_Y, or CONE_Z)
 * @param data Cone primitive data with sheet_selection != 0
 * @param out_neg_node Output: negative sense node (interior)
 * @param out_pos_node Output: positive sense node (exterior)
 * @return 0 on success, -1 on error
 */
int alea_expand_1sheet_cone_immediate(alea_system_t* sys,
                                      alea_primitive_type_t type,
                                      const alea_primitive_data_t* data,
                                      alea_node_id_t* out_neg_node,
                                      alea_node_id_t* out_pos_node);


#endif // ALEA_MACROBODY_H
