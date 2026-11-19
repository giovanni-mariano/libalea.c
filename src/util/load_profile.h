// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
// SPDX-License-Identifier: MPL-2.0
#ifndef ALEA_LOAD_PROFILE_H
#define ALEA_LOAD_PROFILE_H
#include "compat.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
/* Stage timers only; no per-query instrumentation in production hot paths. */
static inline double alea_load_profile_start(void) {
    const char* value=getenv("ALEA_PROFILE_LOAD");
    return value && *value && strcmp(value,"0") ? alea_monotonic_seconds() : -1;
}
static inline void alea_load_profile_end(const char* stage,double start) {
    if(start>=0)fprintf(stderr,"[alea-load-profile] %-28s %.9f s\n",
                        stage,alea_monotonic_seconds()-start);
}
#endif
