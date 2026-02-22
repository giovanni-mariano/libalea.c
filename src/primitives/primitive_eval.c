// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

#include "primitive_eval.h"
#include "primitive_desc.h"
#include "util/alea_log.h"

/**
 * @file primitive_eval.c
 * @brief Implementation of signed distance evaluation for primitives
 *
 * Uses the primitive descriptor table for dispatch.
 */

double alea_primitive_eval(
    alea_primitive_type_t type,
    const alea_primitive_data_t* data,
    double x, double y, double z
) {
    if (!data) return 1.0;

    const alea_primitive_desc_t* desc = alea_primitive_get_desc(type);
    if (!desc || !desc->eval) {
        ALEA_LOG_WARN("Unknown primitive type %d in evaluation", type);
        return 1.0;
    }

    return desc->eval(data, x, y, z);
}

double alea_primitive_eval_with_sense(
    alea_primitive_type_t type,
    const alea_primitive_data_t* data,
    int8_t sense,
    double x, double y, double z
) {
    double result = alea_primitive_eval(type, data, x, y, z);
    return sense < 0 ? result : -result;
}

