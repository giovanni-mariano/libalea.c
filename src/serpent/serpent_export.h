// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

#ifndef SERPENT_EXPORT_H
#define SERPENT_EXPORT_H

#include "alea_types.h"

typedef struct export_context export_context_t;

int export_serpent(alea_system_t* sys, export_context_t* ctx);

#endif /* SERPENT_EXPORT_H */
