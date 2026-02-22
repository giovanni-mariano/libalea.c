// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/*
 * mcnp_roundtrip.c - Read an MCNP input and re-export it in MCNP format
 *
 * Usage:
 *   ./mcnp_roundtrip <input.inp> [output.inp]
 *
 * If output is not specified, writes to stdout.
 *
 * Build:
 *   gcc -o mcnp_roundtrip mcnp_roundtrip.c -I../include ../bin/libalea.a -lm
 */

#include <stdio.h>
#include <string.h>
#include <alea.h>

int main(int argc, char** argv) {
    if (argc < 2 || strcmp(argv[1], "-h") == 0 || strcmp(argv[1], "--help") == 0) {
        fprintf(stderr, "Usage: %s <input.inp> [output.inp]\n", argv[0]);
        fprintf(stderr, "Reads an MCNP input and re-exports it in MCNP format.\n");
        fprintf(stderr, "If output is not specified, writes to stdout.\n");
        return argc < 2 ? 1 : 0;
    }

    const char* input_file = argv[1];
    const char* output_file = argc >= 3 ? argv[2] : NULL;

    printf("Alea %s - MCNP roundtrip\n\n", alea_version());

    /* Load MCNP input */
    printf("Loading: %s\n", input_file);
    alea_system_t* sys = alea_load_mcnp(input_file);
    if (!sys) {
        fprintf(stderr, "Error: %s\n", alea_error());
        return 1;
    }

    alea_print_summary(sys);

    /* Export back to MCNP */
    int rc;
    if (output_file) {
        printf("\nExporting to: %s\n", output_file);
        rc = alea_export_mcnp(sys, output_file);
    } else {
        printf("\n--- MCNP output ---\n");
        rc = alea_export_mcnp_stream(sys, stdout);
    }

    if (rc != 0) {
        fprintf(stderr, "Export failed: %s\n", alea_error());
        alea_destroy(sys);
        return 1;
    }

    printf("Done.\n");
    alea_destroy(sys);
    return 0;
}
