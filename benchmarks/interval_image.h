// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
// SPDX-License-Identifier: MPL-2.0
#ifndef ALEA_BENCH_INTERVAL_IMAGE_H
#define ALEA_BENCH_INTERVAL_IMAGE_H
#include "raycast/regional/interval_index.h"

int interval_index_save(interval_index_t* index,const char* path);
interval_index_t* interval_index_restore(alea_system_t* sys,const alea_bbox_t* domain,const char* path);

#endif
