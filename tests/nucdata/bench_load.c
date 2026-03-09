#define _POSIX_C_SOURCE 200809L

// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/**
 * @file bench_load.c
 * @brief Benchmark: load all nuclides from FENDL xsdir directory
 *
 * Measures wall-clock time for:
 *   1. xsdir parsing
 *   2. Loading all nuclides (ACE read + decode)
 *   3. Total
 *
 * Also reports memory usage and per-nuclide statistics.
 */

#include "alea_nucdata.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

static double wtime(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + ts.tv_nsec * 1e-9;
}

/* Read peak RSS from /proc/self/status (Linux) */
static long peak_rss_kb(void) {
    FILE* f = fopen("/proc/self/status", "r");
    if (!f) return -1;
    char line[256];
    long rss = -1;
    while (fgets(line, sizeof(line), f)) {
        if (strncmp(line, "VmHWM:", 6) == 0) {
            sscanf(line + 6, "%ld", &rss);
            break;
        }
    }
    fclose(f);
    return rss;
}

int main(int argc, char** argv) {
    const char* datadir = NULL;

    if (argc > 1) {
        datadir = argv[1];
    } else {
        /* Default FENDL path */
        datadir = "fendl-FENDL-3.2c-neutron-ace/neutron/ace";
    }

    printf("=== nucdata loading benchmark ===\n\n");
    printf("Data directory: %s\n\n", datadir);

    double t0, t1;

    /* Phase 1: xsdir */
    t0 = wtime();
    alea_nuc_xsdir_t* xsdir = alea_nuc_xsdir_load_dir(datadir);
    t1 = wtime();
    if (!xsdir) {
        fprintf(stderr, "alea_nuc_xsdir_load_dir() failed\n");
        return 1;
    }

    size_t n_entries = alea_nuc_xsdir_count(xsdir);
    double t_xsdir = t1 - t0;
    printf("Phase 1 — xsdir parsing\n");
    printf("  Entries:  %zu\n", n_entries);
    printf("  Time:     %.3f ms\n\n", t_xsdir * 1000.0);

    /* Collect ZAIDs by scanning .xsd files */
    char cmd[512];
    snprintf(cmd, sizeof(cmd), "grep -h '' %s/*.xsd 2>/dev/null | awk '{print $1}'", datadir);
    FILE* pipe = popen(cmd, "r");
    if (!pipe) {
        fprintf(stderr, "Cannot list ZAIDs\n");
        alea_nuc_xsdir_free(xsdir);
        return 1;
    }

    char zaids[500][32];
    int n_zaids = 0;
    char buf[256];
    while (fgets(buf, sizeof(buf), pipe) && n_zaids < 500) {
        buf[strcspn(buf, "\n")] = '\0';
        if (strlen(buf) > 0) {
            strncpy(zaids[n_zaids], buf, 31);
            zaids[n_zaids][31] = '\0';
            n_zaids++;
        }
    }
    pclose(pipe);

    printf("Phase 2 — loading %d nuclides (ACE read + decode)\n", n_zaids);

    /* Phase 2: load all nuclides */
    long total_energies = 0;
    long total_reactions = 0;
    int loaded = 0;
    int failed = 0;
    double t_min = 1e30, t_max = 0.0;
    const char* slowest_zaid = "";

    t0 = wtime();
    for (int i = 0; i < n_zaids; i++) {
        double ti = wtime();
        alea_nuc_nuclide_t* nuc = alea_nuc_load_nuclide(xsdir, zaids[i]);
        double dt = wtime() - ti;

        if (nuc) {
            loaded++;
            total_energies += nuc->n_energies;
            total_reactions += nuc->n_reactions;
            if (dt < t_min) t_min = dt;
            if (dt > t_max) { t_max = dt; slowest_zaid = zaids[i]; }
            alea_nuc_nuclide_free(nuc);
        } else {
            failed++;
        }
    }
    t1 = wtime();
    double t_load = t1 - t0;

    printf("  Loaded:   %d / %d\n", loaded, n_zaids);
    if (failed > 0)
        printf("  Failed:   %d\n", failed);
    printf("  Time:     %.3f ms  (%.3f ms/nuclide avg)\n",
           t_load * 1000.0, (t_load / loaded) * 1000.0);
    printf("  Fastest:  %.3f ms\n", t_min * 1000.0);
    printf("  Slowest:  %.3f ms  (%s)\n", t_max * 1000.0, slowest_zaid);
    printf("\n");

    /* Statistics */
    printf("Data statistics:\n");
    printf("  Total energy grid points: %ld  (avg %ld/nuclide)\n",
           total_energies, total_energies / loaded);
    printf("  Total reactions:          %ld  (avg %ld/nuclide)\n",
           total_reactions, total_reactions / loaded);
    printf("\n");

    /* Memory */
    long rss = peak_rss_kb();
    if (rss > 0) {
        printf("Memory:\n");
        printf("  Peak RSS: %.1f MB\n\n", rss / 1024.0);
    }

    /* Total */
    double t_total = t_xsdir + t_load;
    printf("TOTAL: %.3f ms  (xsdir %.1f%% + load %.1f%%)\n",
           t_total * 1000.0,
           (t_xsdir / t_total) * 100.0,
           (t_load / t_total) * 100.0);

    alea_nuc_xsdir_free(xsdir);
    return 0;
}
