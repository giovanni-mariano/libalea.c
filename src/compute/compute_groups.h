// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
// SPDX-License-Identifier: MPL-2.0
#ifndef ALEA_COMPUTE_GROUPS_H
#define ALEA_COMPUTE_GROUPS_H
#include "compute_internal.h"
typedef struct {uint32_t begin,end;} compute_group_span;
typedef int (*compute_group_emit)(void*,compute_group_span);
typedef struct {uint32_t root,count,offset,groups;} compute_group_program;
typedef struct {
    compute_group_program* programs;
    compute_group_span* spans;
    size_t program_count,program_capacity,span_count,span_capacity;
} compute_group_plan;
/* Parse unthreaded bytecode without reassociation. Callback failure returns -1
 * without replacing its error. Structural failures and interruption set errors. */
int alea_compute_group_program(const compute_scene*,compute_cell,compute_group_emit,void*);
/* 1: owned plan; 0: no groups or optional allocation budget unavailable;
 * -1: invalid bytecode or interruption. Output is empty on failure/fallback. */
int alea_compute_groups_plan(const compute_scene*,size_t budget,compute_group_plan*);
void alea_compute_groups_free(compute_group_plan*);
#endif
