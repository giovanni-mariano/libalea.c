// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
// SPDX-License-Identifier: MPL-2.0
#define ALEA_RAY_EXECUTION_EXPERIMENT 1
#define RAY_EXECUTE_BIND(e) do { \
    [e setBuffer:s->group_maps offset:0 atIndex:19]; \
    [e setBuffer:s->group_refs offset:0 atIndex:20]; } while(0)
#include "metal_ray_groups.m"
