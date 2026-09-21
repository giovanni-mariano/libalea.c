// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/* Fixed-source fusion-neutron example. The source box must lie inside the
 * modeled geometry. Output is per source neutron, without volume or source-
 * rate normalization. Use --coupled only with neutron photon-production and
 * photoatomic ACE tables available for the geometry materials. */

#include "alea_mcnp.h"
#include "alea_transport.h"

#include <errno.h>
#include <limits.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int parse_double(const char* text, double* value) {
    char* end = NULL;
    errno = 0;
    *value = strtod(text, &end);
    return errno == 0 && end != text && *end == '\0' && isfinite(*value);
}

static int parse_histories(const char* text, uint32_t* value) {
    char* end = NULL;
    errno = 0;
    unsigned long long parsed = strtoull(text, &end, 10);
    if (errno || end == text || *end || text[0] == '-' ||
        parsed == 0 || parsed > UINT32_MAX) return 0;
    *value = (uint32_t)parsed;
    return 1;
}

static int parse_surface_id(const char* text, int* value) {
    char* end = NULL;
    errno = 0;
    long parsed = strtol(text, &end, 10);
    if (errno || end == text || *end || parsed <= 0 || parsed > INT_MAX)
        return 0;
    *value = (int)parsed;
    return 1;
}

static double standard_error(const alea_tally_view_t* view, size_t bin) {
    if (view->histories < 2) return 0.0;
    double n = (double)view->histories;
    double variance = (view->sum_squared[bin] -
        view->sum[bin] * view->sum[bin] / n) / (n - 1.0);
    return sqrt(fmax(0.0, variance / n));
}

int main(int argc, char** argv) {
    if (argc < 12) {
        fprintf(stderr, "Usage: %s geometry.i xsdir-or-xsd-dir histories energy_MeV "
            "xmin xmax ymin ymax zmin zmax output.csv "
            "[--coupled] [--vacuum-surface ID]\n", argv[0]);
        return 2;
    }
    int coupled = 0;
    int vacuum_surface = 0;
    for (int arg = 12; arg < argc; ++arg) {
        if (strcmp(argv[arg], "--coupled") == 0 && !coupled) {
            coupled = 1;
        } else if (strcmp(argv[arg], "--vacuum-surface") == 0 &&
                   vacuum_surface == 0 && arg + 1 < argc &&
                   parse_surface_id(argv[arg + 1], &vacuum_surface)) {
            ++arg;
        } else {
            fprintf(stderr, "Invalid source option: %s\n", argv[arg]);
            return 2;
        }
    }
    uint32_t histories;
    double energy;
    alea_transport_box_source_t source = {
        .particle_type = ALEA_NUC_PARTICLE_NEUTRON, .weight = 1.0
    };
    if (!parse_histories(argv[3], &histories) ||
        !parse_double(argv[4], &energy) || energy <= 0.0 ||
        !parse_double(argv[5], &source.lower[0]) ||
        !parse_double(argv[6], &source.upper[0]) ||
        !parse_double(argv[7], &source.lower[1]) ||
        !parse_double(argv[8], &source.upper[1]) ||
        !parse_double(argv[9], &source.lower[2]) ||
        !parse_double(argv[10], &source.upper[2])) {
        fprintf(stderr, "Invalid history count, energy, or source box\n");
        return 2;
    }
    for (int axis = 0; axis < 3; ++axis) {
        if (source.upper[axis] < source.lower[axis]) {
            fprintf(stderr, "Source box upper bound is below lower bound\n");
            return 2;
        }
    }
    source.energy = energy;
    int status = 1;
    mcnp_model_t* model = NULL;
    alea_nuc_xsdir_t* xsdir = NULL;
    alea_nuc_cell_bindings_t* bindings = NULL;
    alea_tally_plan_t* plan = NULL;
    alea_transport_result_t result = {0};
    FILE* csv = NULL;

    model = mcnp_load(argv[1]);
    if (!model) {
        fprintf(stderr, "Geometry load failed: %s\n", argv[1]);
        goto done;
    }
    if (vacuum_surface && alea_surface_set_boundary(model->sys,
            vacuum_surface, ALEA_BOUNDARY_VACUUM) != 0) {
        fprintf(stderr, "Vacuum surface %d was not found\n", vacuum_surface);
        goto done;
    }
    xsdir = alea_nuc_xsdir_load(argv[2]);
    if (!xsdir) xsdir = alea_nuc_xsdir_load_dir(argv[2]);
    if (!xsdir) {
        fprintf(stderr, "xsdir or .xsd directory load failed: %s\n", argv[2]);
        goto done;
    }
    alea_nuc_prepare_requirements_t requirements = {
        .required_capabilities = ALEA_NUC_CAP_CONTINUOUS_NEUTRON |
            (coupled ? ALEA_NUC_CAP_PHOTON_PRODUCTION : 0)
    };
    uint32_t particles = ALEA_NUC_BIND_NEUTRON |
        (coupled ? ALEA_NUC_BIND_PHOTON : 0);
    alea_error_t err = alea_nuc_cell_bindings_prepare(model->sys, xsdir,
        particles, &requirements, NULL, &bindings);
    if (err != ALEA_OK) {
        fprintf(stderr, "Nuclear-data binding failed: %s\n",
                alea_error_string(err));
        goto done;
    }
    plan = alea_tally_plan_create(model->sys);
    if (!plan) {
        fprintf(stderr, "Tally plan allocation failed\n");
        goto done;
    }
    alea_tally_spec_t score = {
        .domain = ALEA_TALLY_CELL, .score = ALEA_TALLY_TRACK_LENGTH,
        .particle_mask = ALEA_TALLY_NEUTRON
    };
    err = alea_tally_plan_add(plan, &score, NULL);
    if (err == ALEA_OK && coupled) {
        score.score = ALEA_TALLY_LOCAL_DEPOSITION;
        score.particle_mask = ALEA_TALLY_PHOTON;
        err = alea_tally_plan_add(plan, &score, NULL);
    }
    if (err != ALEA_OK) {
        fprintf(stderr, "Tally definition failed: %s\n",
                alea_error_string(err));
        goto done;
    }
    alea_transport_options_t options = {
        .histories = histories, .seed = 1,
        .max_events_per_history = 100000,
        .max_segment_distance = 100.0,
        .tally_plan = plan
    };
    alea_transport_failure_t failure = {0};
    err = alea_transport_run_sampled_source(model->sys, bindings,
        alea_transport_sample_box_isotropic, &source, &options,
        &result, &failure);
    if (err != ALEA_OK) {
        fprintf(stderr, "Transport failed: %s, history %u, cell %d, "
            "position (%.6g, %.6g, %.6g)\n",
            alea_error_string(err), failure.history_id, failure.cell_id,
            failure.position[0], failure.position[1], failure.position[2]);
        goto done;
    }
    csv = fopen(argv[11], "w");
    if (!csv) {
        fprintf(stderr, "Cannot open output: %s\n", argv[11]);
        goto done;
    }
    fprintf(csv, "score,cell_id,mean_per_source,standard_error,unit\n");
    for (size_t tally = 0; tally < alea_tally_results_count(result.tallies);
         ++tally) {
        alea_tally_view_t view;
        err = alea_tally_results_view(result.tallies, tally, &view);
        if (err != ALEA_OK) goto done;
        const char* name = tally == 0 ? "neutron_track_length" :
                                         "photon_local_deposition";
        const char* unit = tally == 0 ? "cm" : "MeV";
        for (size_t bin = 0; bin < view.bin_count; ++bin)
            fprintf(csv, "%s,%d,%.17g,%.17g,%s\n", name,
                view.bin_ids[bin], view.sum[bin] / view.histories,
                standard_error(&view, bin), unit);
    }
    if (ferror(csv)) {
        fprintf(stderr, "Failed writing output: %s\n", argv[11]);
        goto done;
    }
    printf("%u source histories, %llu collisions, %llu leaked; wrote %s\n",
        histories, (unsigned long long)result.collisions,
        (unsigned long long)result.leaked, argv[11]);
    status = 0;
done:
    if (csv && fclose(csv) != 0) status = 1;
    alea_transport_result_free(&result);
    alea_tally_plan_free(plan);
    alea_nuc_cell_bindings_free(bindings);
    alea_nuc_xsdir_free(xsdir);
    mcnp_model_destroy(model);
    return status;
}
