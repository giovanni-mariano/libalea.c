// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
// SPDX-License-Identifier: MPL-2.0
#ifndef ALEA_PRIMITIVE_CODEC_H
#define ALEA_PRIMITIVE_CODEC_H
#include "alea_types.h"
#include <stddef.h>
int alea_primitive_from_values(alea_primitive_type_t type,const double* values,alea_primitive_data_t* data);
size_t alea_primitive_to_values(alea_primitive_type_t type,const alea_primitive_data_t* data,double* values);
#endif
