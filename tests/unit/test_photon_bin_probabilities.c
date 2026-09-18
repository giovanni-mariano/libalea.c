// SPDX-FileCopyrightText: 2026 F4E
// SPDX-License-Identifier: EUPL-1.2

#include "alea_nucdata.h"
#include <assert.h>
#include <math.h>
#include <stdio.h>

static alea_nuc_prepared_photon_production_t prepared;
static double incident_energy;

static void close_to(double actual, double expected) {
    if (fabs(actual - expected) > 1e-12) {
        fprintf(stderr, "expected %.17g, got %.17g\n", expected, actual);
        assert(0);
    }
}

static void check_bin(double low, double high, int final, double expected) {
    static int check_number = 0;
    check_number++;
    double edges[] = {low, high, high + 1.0};
    double probabilities[] = {-1.0, -1.0};
    assert(alea_nuc_prepared_photon_bin_probabilities(
        &prepared, incident_energy, edges, final ? 2 : 3,
        probabilities) == ALEA_OK);
    if (fabs(probabilities[0] - expected) > 1e-12)
        fprintf(stderr, "bin check %d: ", check_number);
    close_to(probabilities[0], expected);
}

int main(void) {
    alea_nuc_nuclide_t nuclide = {0};
    alea_nuc_photon_production_t production = {0};
    prepared.nuclide = &nuclide;
    prepared.production = &production;
    double incident_one[] = {1.0};
    int histogram[] = {1}, linear[] = {2}, no_lines[] = {0};
    int two_out[] = {2};
    double energy[] = {0.0, 2.0};
    double flat[] = {0.5, 0.5}, rising[] = {0.0, 1.0};
    double cdf[] = {0.0, 1.0};
    double* out[] = {energy}, *pdf[] = {flat}, *cumulative[] = {cdf};
    alea_nuc_energy_dist_t law = {0};
    production.spectrum = &law;
    law.law = ALEA_NUC_ELAW_CONT_TABULAR;
    law.tab.n_ein = 1;
    law.tab.ein = incident_one;
    law.tab.interpolation = histogram;
    law.tab.n_discrete = no_lines;
    law.tab.n_eout = two_out;
    law.tab.eout = out;
    law.tab.pdf = pdf;
    law.tab.cdf = cumulative;
    incident_energy = 1.0;
    check_bin(0.0, 0.5, 0, 0.25);
    check_bin(0.5, 2.0, 1, 0.75);
    check_bin(2.0, 3.0, 1, 0.0);

    law.tab.interpolation = linear;
    pdf[0] = rising;
    incident_energy = 1.0;
    check_bin(0.0, 1.0, 0, 0.25);
    check_bin(1.0, 2.0, 1, 0.75);

    double mixed_out[] = {0.5, 1.0, 3.0};
    double mixed_pdf[] = {0.0, 0.25, 0.25};
    double mixed_cdf[] = {0.5, 0.5, 1.0};
    int one_line[] = {1}, three_out[] = {3};
    out[0] = mixed_out;
    pdf[0] = mixed_pdf;
    cumulative[0] = mixed_cdf;
    law.tab.interpolation = histogram;
    law.tab.n_discrete = one_line;
    law.tab.n_eout = three_out;
    incident_energy = 1.0;
    check_bin(0.0, 0.5, 0, 0.0);
    check_bin(0.5, 1.0, 0, 0.5);
    check_bin(1.0, 2.0, 0, 0.25);
    check_bin(2.0, 3.0, 1, 0.25);

    double incident_two[] = {1.0, 3.0};
    double upper_out[] = {2.0, 4.0};
    double upper_pdf[] = {0.5, 0.5};
    double upper_cdf[] = {0.0, 1.0};
    double* pair_out[] = {energy, upper_out};
    double* pair_pdf[] = {flat, upper_pdf};
    double* pair_cdf[] = {cdf, upper_cdf};
    int pair_interp[] = {1, 1}, pair_no_lines[] = {0, 0};
    int pair_nout[] = {2, 2};
    law.tab.n_ein = 2;
    law.tab.ein = incident_two;
    law.tab.interpolation = pair_interp;
    law.tab.n_discrete = pair_no_lines;
    law.tab.n_eout = pair_nout;
    law.tab.eout = pair_out;
    law.tab.pdf = pair_pdf;
    law.tab.cdf = pair_cdf;
    incident_energy = 2.0;
    check_bin(0.0, 1.0, 0, 0.0);
    check_bin(1.0, 2.0, 0, 0.5);
    check_bin(2.0, 3.0, 1, 0.5);

    double shifted_mixed_out[] = {1.5, 3.0, 5.0};
    double* mixed_pair_out[] = {mixed_out, shifted_mixed_out};
    double* mixed_pair_pdf[] = {mixed_pdf, mixed_pdf};
    double* mixed_pair_cdf[] = {mixed_cdf, mixed_cdf};
    int mixed_pair_lines[] = {1, 1}, mixed_pair_count[] = {3, 3};
    law.tab.n_discrete = mixed_pair_lines;
    law.tab.n_eout = mixed_pair_count;
    law.tab.eout = mixed_pair_out;
    law.tab.pdf = mixed_pair_pdf;
    law.tab.cdf = mixed_pair_cdf;
    incident_energy = 2.0;
    check_bin(0.0, 1.0, 0, 0.0);
    check_bin(1.0, 2.0, 0, 0.5);
    check_bin(2.0, 3.0, 0, 0.25);
    check_bin(3.0, 4.0, 1, 0.25);

    double lower_line[] = {1.0}, upper_line[] = {3.0};
    double line_cdf[] = {1.0}, line_pdf[] = {0.0};
    double* lines_out[] = {lower_line, upper_line};
    double* lines_pdf[] = {line_pdf, line_pdf};
    double* lines_cdf[] = {line_cdf, line_cdf};
    int lines_count[] = {1, 1}, line_interpolation[] = {0, 0};
    law.tab.interpolation = line_interpolation;
    law.tab.n_discrete = lines_count;
    law.tab.n_eout = lines_count;
    law.tab.eout = lines_out;
    law.tab.pdf = lines_pdf;
    law.tab.cdf = lines_cdf;
    incident_energy = 2.0;
    check_bin(1.0, 2.0, 0, 0.0);
    check_bin(2.0, 3.0, 0, 1.0);

    int nbt[] = {2}, incident_histogram[] = {1};
    law.tab.n_regions = 1;
    law.tab.nbt = nbt;
    law.tab.interp = incident_histogram;
    incident_energy = 2.0;
    check_bin(1.0, 2.0, 0, 1.0);
    check_bin(2.0, 3.0, 0, 0.0);

    alea_nuc_energy_dist_t next = {0};
    law.next = &next;
    double chain_edges[] = {0.0, 1.0};
    double chain_probability = -1.0;
    assert(alea_nuc_prepared_photon_bin_probabilities(
        &prepared, 2.0, chain_edges, 2, &chain_probability) ==
        ALEA_ERR_UNSUPPORTED);
    law.next = NULL;
    law.law = ALEA_NUC_ELAW_DISCRETE_PHOTON;
    law.discrete_photon_energy = 1.0;
    law.discrete_photon_primary = 2;
    law.discrete_photon_awr = 1.0;
    check_bin(1.0, 2.0, 0, 0.0);
    check_bin(1.0, 2.0, 1, 1.0);
    law.law = ALEA_NUC_ELAW_LEVEL;
    law.level_A = 0.5;
    law.level_Q = 2.0;
    check_bin(2.0, 3.0, 1, 1.0);
    double bad_edges[] = {1.0, 1.0};
    assert(alea_nuc_prepared_photon_bin_probabilities(
        &prepared, 2.0, bad_edges, 2, &chain_probability) ==
        ALEA_ERR_INVALID_ARG);
    puts("analytic photon bin probability checks passed");
    return 0;
}
