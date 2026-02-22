// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/*
 * mc_convert.c - Convert between MCNP and OpenMC geometry formats
 *
 * Usage:
 *   ./mc_convert [options] <input> [output]
 *
 * Options:
 *   -if, --input-format  FORMAT   Input format:  mcnp | openmc (auto-detected)
 *   -of, --output-format FORMAT   Output format: mcnp | openmc (auto-detected)
 *   --no-dedup                     Disable surface deduplication
 *   -h,  --help                   Show this help
 *
 * When formats are not specified, auto-detection uses file extension:
 *   .xml  -> OpenMC input,  MCNP output  (default output: output.inp)
 *   other -> MCNP input,    OpenMC output (default output: model.xml)
 *
 * Build:
 *   make mc_convert
 */

#include <stdio.h>
#include <string.h>
#include <alea.h>

enum format { FMT_AUTO, FMT_MCNP, FMT_OPENMC };

static void print_usage(const char* prog) {
    fprintf(stderr, "Usage: %s [options] <input> [output]\n", prog);
    fprintf(stderr, "\nConvert between MCNP and OpenMC geometry formats.\n");
    fprintf(stderr, "\nOptions:\n");
    fprintf(stderr, "  -if, --input-format  FORMAT   mcnp | openmc (auto-detected)\n");
    fprintf(stderr, "  -of, --output-format FORMAT   mcnp | openmc (auto-detected)\n");
    fprintf(stderr, "  --no-dedup                     Disable surface deduplication\n");
    fprintf(stderr, "  -h,  --help                   Show this help\n");
    fprintf(stderr, "\nAuto-detection uses file extension (.xml = OpenMC, otherwise MCNP).\n");
}

static int ends_with(const char* str, const char* suffix) {
    size_t slen = strlen(str);
    size_t xlen = strlen(suffix);
    if (slen < xlen) return 0;
    return strcmp(str + slen - xlen, suffix) == 0;
}

static enum format parse_format(const char* s) {
    if (strcmp(s, "mcnp") == 0)   return FMT_MCNP;
    if (strcmp(s, "openmc") == 0) return FMT_OPENMC;
    return FMT_AUTO;
}

int main(int argc, char** argv) {
    enum format in_fmt = FMT_AUTO;
    enum format out_fmt = FMT_AUTO;
    const char* input_file = NULL;
    const char* output_file = NULL;
    int no_dedup = 0;

    /* Parse arguments */
    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0) {
            print_usage(argv[0]);
            return 0;
        } else if (strcmp(argv[i], "--no-dedup") == 0) {
            no_dedup = 1;
        } else if (strcmp(argv[i], "-if") == 0 || strcmp(argv[i], "--input-format") == 0) {
            if (++i >= argc) { fprintf(stderr, "Error: %s requires an argument\n", argv[i-1]); return 1; }
            in_fmt = parse_format(argv[i]);
            if (in_fmt == FMT_AUTO) { fprintf(stderr, "Error: unknown format '%s' (use mcnp or openmc)\n", argv[i]); return 1; }
        } else if (strcmp(argv[i], "-of") == 0 || strcmp(argv[i], "--output-format") == 0) {
            if (++i >= argc) { fprintf(stderr, "Error: %s requires an argument\n", argv[i-1]); return 1; }
            out_fmt = parse_format(argv[i]);
            if (out_fmt == FMT_AUTO) { fprintf(stderr, "Error: unknown format '%s' (use mcnp or openmc)\n", argv[i]); return 1; }
        } else if (!input_file) {
            input_file = argv[i];
        } else if (!output_file) {
            output_file = argv[i];
        } else {
            fprintf(stderr, "Error: unexpected argument '%s'\n", argv[i]);
            return 1;
        }
    }

    if (!input_file) {
        print_usage(argv[0]);
        return 1;
    }

    /* Auto-detect input format from extension */
    if (in_fmt == FMT_AUTO)
        in_fmt = ends_with(input_file, ".xml") ? FMT_OPENMC : FMT_MCNP;

    /* Auto-detect output format as the opposite of input */
    if (out_fmt == FMT_AUTO)
        out_fmt = (in_fmt == FMT_MCNP) ? FMT_OPENMC : FMT_MCNP;

    /* Default output filename */
    if (!output_file)
        output_file = (out_fmt == FMT_OPENMC) ? "model.xml" : "output.inp";

    printf("Alea %s - MC format converter\n\n", alea_version());

    /* Load */
    alea_system_t* sys;
    if (in_fmt == FMT_OPENMC) {
        printf("Loading OpenMC XML: %s\n", input_file);
        sys = alea_load_openmc(input_file);
    } else {
        printf("Loading MCNP input: %s\n", input_file);
        sys = alea_load_mcnp(input_file);
    }
    if (!sys) {
        fprintf(stderr, "Error: %s\n", alea_error());
        return 1;
    }

    if (no_dedup) {
        alea_config_t cfg = alea_get_config(sys);
        cfg.dedup = false;
        alea_set_config(sys, &cfg);
    }

    alea_print_summary(sys);

    /* Export */
    int rc;
    if (out_fmt == FMT_OPENMC) {
        printf("\nExporting to OpenMC XML: %s\n", output_file);
        rc = alea_export_openmc(sys, output_file);
    } else {
        printf("\nExporting to MCNP: %s\n", output_file);
        rc = alea_export_mcnp(sys, output_file);
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
