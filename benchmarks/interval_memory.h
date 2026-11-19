// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
// SPDX-License-Identifier: MPL-2.0
#ifndef ALEA_INTERVAL_MEMORY_H
#define ALEA_INTERVAL_MEMORY_H
#include "util/compat.h"
#include <stdio.h>
#include <stdlib.h>
static void interval_memory(const char* stage) {
    if (!getenv("ALEA_INTERVAL_PROFILE")) return;
    FILE* f = fopen("/proc/self/status","r");
    if (!f) return;
    unsigned long rss = 0, peak = 0;
    char line[256];
    while (fgets(line,sizeof(line),f)) {
        if (sscanf(line,"VmRSS: %lu",&rss) == 1) break;
        (void)sscanf(line,"VmHWM: %lu",&peak);
    }
    fclose(f);
    fprintf(stderr,"[interval] %.6f %s rss=%lu KiB peak=%lu KiB\n",
            alea_monotonic_seconds(),stage,rss,peak);
}
#endif
