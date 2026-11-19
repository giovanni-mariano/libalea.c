// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
// SPDX-License-Identifier: MPL-2.0
#ifndef ALEA_REGIONAL_H
#define ALEA_REGIONAL_H

#include "core/alea_system.h"

/* Strictly positive throughout a closed box under the scalar evaluator. */
bool alea_regional_proves_outside(const alea_system_t* sys, alea_node_id_t root,
                                  const alea_bbox_t* box, size_t* remaining);
int alea_regional_build(alea_system_t* sys);
void alea_regional_free(struct alea_regional_cache* cache);
/* Statistics for the benchmark adapter, including retained allocation capacity. */
void alea_regional_stats(const alea_system_t* sys, size_t* roots, size_t* bytes);
/* Returns false when the caller must evaluate the full expression. */
bool alea_regional_contains(const alea_system_t* sys, alea_node_id_t root,
                            double x, double y, double z, bool* inside);

#endif
