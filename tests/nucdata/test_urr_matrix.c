// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/** @file test_urr_matrix.c
 *  @brief Full-library unresolved-resonance representation matrix
 */

#include "alea_nucdata.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct {
    double value;
    int calls;
} fixed_rng_t;

typedef struct {
    size_t tables;
    size_t evaluations;
    size_t expected_rejections;
    int interp2;
    int interp5;
    int absolute;
    int factors;
    int ilf_negative;
    int ilf_zero;
    int ilf_positive;
    int ioa_negative;
    int ioa_zero;
    int ioa_positive;
} matrix_t;

static double fixed_uniform(void* context) {
    fixed_rng_t* rng = context;
    rng->calls++;
    return rng->value;
}

static int ends_with(const char* text, const char* suffix) {
    size_t n = strlen(text), m = strlen(suffix);
    return n >= m && strcmp(text + n - m, suffix) == 0;
}

static alea_nuc_prepare_issue_t expected_endfb71_rejection(
    const char* zaid) {
    if (strcmp(zaid, "32070.71c") == 0)
        return ALEA_NUC_PREP_UNSUPPORTED_URR;
    if (strcmp(zaid, "63153.71c") == 0)
        return ALEA_NUC_PREP_INVALID_ENERGY_DISTRIBUTION;
    return ALEA_NUC_PREP_OK;
}

static void record_representation(matrix_t* matrix, const alea_nuc_urr_t* urr) {
    matrix->interp2 |= urr->interp == 2;
    matrix->interp5 |= urr->interp == 5;
    matrix->absolute |= !urr->multiply_smooth;
    matrix->factors |= urr->multiply_smooth;
    matrix->ilf_negative |= urr->inelastic_flag < 0;
    matrix->ilf_zero |= urr->inelastic_flag == 0;
    matrix->ilf_positive |= urr->inelastic_flag > 0;
    matrix->ioa_negative |= urr->absorption_flag < 0;
    matrix->ioa_zero |= urr->absorption_flag == 0;
    matrix->ioa_positive |= urr->absorption_flag > 0;
}

static int valid_evaluation(const alea_nuc_evaluation_t* evaluation) {
    const double sum = evaluation->macro_elastic + evaluation->macro_thermal +
        evaluation->macro_absorption + evaluation->macro_neutron_emission;
    const double scale = fmax(1.0, evaluation->macro_total);
    return isfinite(evaluation->macro_total) && evaluation->macro_total >= 0.0 &&
           isfinite(evaluation->macro_elastic) &&
           evaluation->macro_elastic >= 0.0 &&
           isfinite(evaluation->macro_thermal) &&
           evaluation->macro_thermal >= 0.0 &&
           isfinite(evaluation->macro_absorption) &&
           evaluation->macro_absorption >= 0.0 &&
           isfinite(evaluation->macro_neutron_emission) &&
           evaluation->macro_neutron_emission >= 0.0 &&
           fabs(evaluation->macro_total - sum) <= 1e-12 * scale;
}

static int evaluate_table(const char* zaid, alea_nuc_nuclide_t* nuclide,
                          matrix_t* matrix) {
    static const double quantiles[] = {
        0.0, 0.125, 0.25, 0.5, 0.75, 0.875, 0.999999999999
    };
    alea_nuc_mat_component_t component = {nuclide, 1.0};
    alea_nuc_material_t material = {&component, 1, 1, NULL};
    alea_nuc_prepare_requirements_t requirements = {
        .required_capabilities = ALEA_NUC_CAP_CONTINUOUS_NEUTRON |
                                 ALEA_NUC_CAP_URR
    };
    alea_nuc_capability_report_t report;
    alea_nuc_prepared_material_t* prepared = NULL;
    alea_error_t error = alea_nuc_prepare_material(
        &material, &requirements, &report, &prepared);
    if (error != ALEA_OK) {
        fprintf(stderr, "%s: preparation failed: %s (MT=%d)\n",
                zaid, report.detail, report.mt);
        return 0;
    }

    alea_nuc_urr_sample_t sample;
    alea_nuc_evaluation_workspace_t workspace = {&sample, 1};
    for (int i = 0; i < nuclide->urr->n_energies; i++) {
        double energies[2] = {nuclide->urr->energy[i], 0.0};
        int count = 1;
        if (i + 1 < nuclide->urr->n_energies) {
            energies[1] = sqrt(nuclide->urr->energy[i] *
                               nuclide->urr->energy[i + 1]);
            count = 2;
        }
        for (int j = 0; j < count; j++) {
            for (size_t q = 0; q < sizeof(quantiles) / sizeof(quantiles[0]); q++) {
                alea_nuc_particle_state_t incident = {
                    ALEA_NUC_PARTICLE_NEUTRON, energies[j],
                    {0.0, 0.0, 1.0}, 1.0, 0.0
                };
                fixed_rng_t rng = {quantiles[q], 0};
                alea_nuc_evaluation_t evaluation;
                error = alea_nuc_evaluate_urr(
                    prepared, &incident, fixed_uniform, &rng,
                    &workspace, &evaluation);
                matrix->evaluations++;
                if (error != ALEA_OK || rng.calls != 1 || !sample.active ||
                    !valid_evaluation(&evaluation)) {
                    fprintf(stderr,
                            "%s: evaluation failed at %.17g, q=%.12g "
                            "(status=%d calls=%d active=%d)\n",
                            zaid, energies[j], quantiles[q], error, rng.calls,
                            sample.active ? 1 : 0);
                    alea_nuc_prepared_material_free(prepared);
                    return 0;
                }
            }
        }
    }
    alea_nuc_prepared_material_free(prepared);
    return 1;
}

static int scan_library(const char* label, const char* path, const char* suffix,
                        size_t expected_tables, int allow_endfb71_rejections,
                        matrix_t* matrix) {
    alea_nuc_xsdir_t* xsdir = alea_nuc_xsdir_load(path);
    if (!xsdir) {
        fprintf(stderr, "%s: cannot load %s\n", label, path);
        return 0;
    }

    int ok = 1;
    for (size_t i = 0; i < xsdir->count; i++) {
        const alea_nuc_xsdir_entry_t* entry = &xsdir->entries[i];
        if (entry->type != ALEA_NUC_TABLE_CONTINUOUS_NEUTRON ||
            !ends_with(entry->zaid, suffix))
            continue;
        alea_nuc_nuclide_t* nuclide = alea_nuc_load_nuclide(xsdir, entry->zaid);
        if (!nuclide) {
            fprintf(stderr, "%s: cannot decode %s\n", label, entry->zaid);
            ok = 0;
            continue;
        }
        if (!nuclide->urr) {
            alea_nuc_nuclide_free(nuclide);
            continue;
        }
        matrix->tables++;
        record_representation(matrix, nuclide->urr);

        alea_nuc_prepare_issue_t expected_issue =
            allow_endfb71_rejections
                ? expected_endfb71_rejection(entry->zaid)
                : ALEA_NUC_PREP_OK;
        if (expected_issue != ALEA_NUC_PREP_OK) {
            alea_nuc_mat_component_t component = {nuclide, 1.0};
            alea_nuc_material_t material = {&component, 1, 1, NULL};
            alea_nuc_prepare_requirements_t requirements = {
                .required_capabilities = ALEA_NUC_CAP_CONTINUOUS_NEUTRON |
                                         ALEA_NUC_CAP_URR
            };
            alea_nuc_capability_report_t report;
            alea_nuc_prepared_material_t* prepared = NULL;
            if (alea_nuc_prepare_material(&material, &requirements, &report,
                                          &prepared) != ALEA_ERR_UNSUPPORTED ||
                prepared != NULL || report.issue != expected_issue) {
                fprintf(stderr, "%s: expected %s to fail closed\n",
                        label, entry->zaid);
                alea_nuc_prepared_material_free(prepared);
                ok = 0;
            } else {
                matrix->expected_rejections++;
            }
        } else if (!evaluate_table(entry->zaid, nuclide, matrix)) {
            ok = 0;
        }
        alea_nuc_nuclide_free(nuclide);
    }
    alea_nuc_xsdir_free(xsdir);

    if (matrix->tables != expected_tables) {
        fprintf(stderr, "%s: expected %zu URR tables, found %zu\n",
                label, expected_tables, matrix->tables);
        ok = 0;
    }
    printf("%s: %zu URR tables, %zu coordinated evaluations, "
           "%zu expected rejections\n",
           label, matrix->tables, matrix->evaluations,
           matrix->expected_rejections);
    return ok;
}

int main(void) {
    const char* endfb80 = getenv("ALEA_ENDFB80_XSDIR");
    const char* endfb71 = getenv("ALEA_ENDFB71_XSDIR");
    if (!endfb80 || !endfb80[0] || !endfb71 || !endfb71[0]) {
        puts("SKIP: ALEA_ENDFB80_XSDIR and ALEA_ENDFB71_XSDIR are not set");
        return 0;
    }

    matrix_t lib80 = {0}, endf71 = {0};
    int ok = scan_library("ENDF/B-VIII.0 Lib80", endfb80, ".00c", 343, 0,
                          &lib80);
    ok &= lib80.interp2 && lib80.absolute && lib80.factors &&
          lib80.ilf_negative && lib80.ilf_positive && lib80.ioa_negative &&
          lib80.ioa_zero && lib80.ioa_positive;
    ok &= lib80.expected_rejections == 0;
    ok &= scan_library("ENDF/B-VII.1", endfb71, ".71c", 268, 1, &endf71);
    ok &= endf71.interp2 && endf71.interp5 && endf71.absolute &&
          endf71.factors && endf71.ilf_negative && endf71.ilf_positive &&
          endf71.ioa_negative;
    ok &= endf71.expected_rejections == 2;
    if (!ok) fputs("URR representation matrix is incomplete\n", stderr);
    return ok ? 0 : 1;
}
