// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/*
 * mc_convert.c - Convert between MCNP, OpenMC, and native ALEA formats
 *
 * Usage:
 *   ./mc_convert [options] <input> [output]
 *
 * Options:
 *   -if, --input-format  FORMAT   Input format:  mcnp | openmc (auto-detected)
 *   -of, --output-format FORMAT   Output format: mcnp | openmc (auto-detected)
 *   --no-dedup                     Disable surface deduplication
 *   -v,  --verbose                Set log level to INFO
 *   -vv                           Set log level to DEBUG
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
#include <alea_mcnp.h>
#include <alea_openmc.h>
#include <alea_serpent.h>
#include <alea_xml.h>

enum format { FMT_AUTO, FMT_MCNP, FMT_OPENMC, FMT_ALEA, FMT_SERPENT };

static void print_usage(const char* prog) {
    fprintf(stderr, "Usage: %s [options] <input> [output]\n", prog);
    fprintf(stderr, "\nConvert between MCNP, OpenMC, ALEA XML, and Serpent geometry formats.\n");
    fprintf(stderr, "\nOptions:\n");
    fprintf(stderr, "  -if, --input-format  FORMAT   mcnp | openmc | alea (auto-detected)\n");
    fprintf(stderr, "  -of, --output-format FORMAT   mcnp | openmc | alea | serpent\n");
    fprintf(stderr, "  --no-dedup                     Disable surface deduplication\n");
    fprintf(stderr, "  -v,  --verbose                Set log level to INFO\n");
    fprintf(stderr, "  -vv                           Set log level to DEBUG\n");
    fprintf(stderr, "  -h,  --help                   Show this help\n");
    fprintf(stderr, "\nXML input is detected from its root element; .alea.xml selects ALEA output.\n");
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
    if (strcmp(s, "alea") == 0) return FMT_ALEA;
    if (strcmp(s, "serpent") == 0) return FMT_SERPENT;
    return FMT_AUTO;
}

static enum format detect_input_format(const char* filename) {
    if (!ends_with(filename, ".xml")) return FMT_MCNP;
    FILE* f = fopen(filename, "rb");
    if (!f) return ends_with(filename, ".alea.xml") ? FMT_ALEA : FMT_OPENMC;
    char buf[4097];
    size_t n = fread(buf, 1, sizeof(buf) - 1, f);
    fclose(f);
    buf[n] = '\0';
    return strstr(buf, "<alea") ? FMT_ALEA : FMT_OPENMC;
}

int main(int argc, char** argv) {
    enum format in_fmt = FMT_AUTO;
    enum format out_fmt = FMT_AUTO;
    const char* input_file = NULL;
    const char* output_file = NULL;
    int no_dedup = 0;
    int verbosity = 0;

    /* Parse arguments */
    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0) {
            print_usage(argv[0]);
            return 0;
        } else if (strcmp(argv[i], "--no-dedup") == 0) {
            no_dedup = 1;
        } else if (strcmp(argv[i], "-vv") == 0) {
            verbosity = 2;
        } else if (strcmp(argv[i], "-v") == 0 || strcmp(argv[i], "--verbose") == 0) {
            if (verbosity < 1) verbosity = 1;
        } else if (strcmp(argv[i], "-if") == 0 || strcmp(argv[i], "--input-format") == 0) {
            if (++i >= argc) { fprintf(stderr, "Error: %s requires an argument\n", argv[i-1]); return 1; }
            in_fmt = parse_format(argv[i]);
            if (in_fmt == FMT_AUTO || in_fmt == FMT_SERPENT) { fprintf(stderr, "Error: unknown or unsupported input format '%s' (use mcnp, openmc, or alea)\n", argv[i]); return 1; }
        } else if (strcmp(argv[i], "-of") == 0 || strcmp(argv[i], "--output-format") == 0) {
            if (++i >= argc) { fprintf(stderr, "Error: %s requires an argument\n", argv[i-1]); return 1; }
            out_fmt = parse_format(argv[i]);
            if (out_fmt == FMT_AUTO) { fprintf(stderr, "Error: unknown format '%s' (use mcnp, openmc, alea, or serpent)\n", argv[i]); return 1; }
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
        in_fmt = detect_input_format(input_file);

    /* Prefer an explicit output suffix, otherwise choose a useful opposite. */
    if (out_fmt == FMT_AUTO) {
        if (output_file && ends_with(output_file, ".alea.xml")) out_fmt = FMT_ALEA;
        else if (output_file && ends_with(output_file, ".xml")) out_fmt = FMT_OPENMC;
        else if (output_file && ends_with(output_file, ".serp")) out_fmt = FMT_SERPENT;
        else out_fmt = (in_fmt == FMT_MCNP) ? FMT_OPENMC : FMT_MCNP;
    }

    /* Default output filename */
    if (!output_file)
        output_file = (out_fmt == FMT_OPENMC) ? "model.xml" :
                      (out_fmt == FMT_ALEA) ? "model.alea.xml" :
                      (out_fmt == FMT_SERPENT) ? "model.serp" : "output.inp";

    printf("Alea %s - MC format converter\n\n", alea_version());

    if (verbosity == 1)
        alea_log_set_level(3);  /* INFO */
    else if (verbosity >= 2)
        alea_log_set_level(4);  /* DEBUG */

    /* Load */
    mcnp_model_t* mcnp_model = NULL;
    openmc_model_t* omc_model = NULL;
    alea_model_t* alea_model = NULL;
    alea_system_t* sys = NULL;

    if (in_fmt == FMT_ALEA) {
        printf("Loading ALEA XML: %s\n", input_file);
        alea_model = alea_xml_load(input_file);
        if (!alea_model) {
            fprintf(stderr, "Error: %s\n", alea_error());
            return 1;
        }
        sys = alea_model_system(alea_model);
    } else if (in_fmt == FMT_OPENMC) {
        printf("Loading OpenMC XML: %s\n", input_file);
        omc_model = openmc_load(input_file);
        if (!omc_model) {
            fprintf(stderr, "Error: %s\n", alea_error());
            return 1;
        }
        sys = omc_model->sys;
    } else {
        printf("Loading MCNP input: %s\n", input_file);
        mcnp_model = mcnp_load(input_file);
        if (!mcnp_model) {
            fprintf(stderr, "Error: %s\n", alea_error());
            return 1;
        }
        sys = mcnp_model->sys;
    }

    if (no_dedup) {
        alea_config_t cfg = alea_get_config(sys);
        cfg.dedup = false;
        alea_set_config(sys, &cfg);
    }

    alea_print_summary(sys);

    /* Export */
    int rc;
    if (out_fmt == FMT_ALEA) {
        printf("\nExporting to ALEA XML: %s\n", output_file);
        alea_model_t* export_model = alea_model;
        if (!export_model && mcnp_model) export_model = mcnp_model_to_alea_model(mcnp_model);
        if (!export_model) export_model = alea_model_adopt(alea_clone(sys));
        rc = export_model ? alea_xml_export(export_model, output_file) : -1;
        if (export_model != alea_model) alea_model_destroy(export_model);
    } else if (out_fmt == FMT_OPENMC) {
        printf("\nExporting to OpenMC XML: %s\n", output_file);
        rc = openmc_export_system(sys, output_file);
    } else if (out_fmt == FMT_SERPENT) {
        printf("\nExporting to Serpent: %s\n", output_file);
        rc = serpent_export_system(sys, output_file);
    } else if (mcnp_model) {
        printf("\nExporting to MCNP: %s\n", output_file);
        rc = mcnp_export(mcnp_model, output_file);
    } else if (alea_model) {
        printf("\nExporting to MCNP: %s\n", output_file);
        mcnp_model_t* export_model = mcnp_model_from_alea_model(alea_model);
        rc = export_model ? mcnp_export(export_model, output_file) : -1;
        mcnp_model_destroy(export_model);
    } else {
        printf("\nExporting to MCNP: %s\n", output_file);
        rc = mcnp_export_system(sys, output_file);
    }
    if (rc != 0) {
        fprintf(stderr, "Export failed: %s\n", alea_error());
    }

    if (mcnp_model)
        mcnp_model_destroy(mcnp_model);
    if (omc_model)
        openmc_model_destroy(omc_model);
    if (alea_model)
        alea_model_destroy(alea_model);

    if (rc != 0) return 1;
    printf("Done.\n");
    return 0;
}
