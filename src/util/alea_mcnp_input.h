// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
// SPDX-License-Identifier: MPL-2.0
#ifndef ALEA_MCNP_INPUT_H
#define ALEA_MCNP_INPUT_H
#include <stddef.h>
/* Expand READ FILE cards. 0=success, -1=I/O, -2=allocation, -3=syntax/depth.
 * cache_conflict flags any source path aliasing the proposed cache file. */
int alea_read_mcnp_input(const char* path,const char* cache_path,char** data,size_t* length,int* cache_conflict);
#endif
