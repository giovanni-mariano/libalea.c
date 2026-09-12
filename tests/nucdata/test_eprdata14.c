// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/**
 * @file test_eprdata14.c
 * @brief Tests for photoatomic data loading from eprdata14
 */

#include "alea_nucdata.h"
#include "alea_test.h"

#include <stdint.h>

#ifndef EPRDATA_XSDIR
#define EPRDATA_XSDIR "eprdata14/xsdir"
#endif
#ifndef EPRDATA_SUFFIX
#define EPRDATA_SUFFIX ".14p"
#endif
#ifndef EPRDATA_LABEL
#define EPRDATA_LABEL "eprdata14"
#endif
#ifndef EPRDATA_FORMAT
#define EPRDATA_FORMAT 3
#endif

static alea_nuc_xsdir_t* xsdir;
static alea_nuc_nuclide_t* h;
static alea_nuc_nuclide_t* pb;

typedef struct { uint64_t state; } test_rng_t;

static double test_rng(void* context) {
    test_rng_t* rng = context;
    rng->state = rng->state * UINT64_C(6364136223846793005) +
                 UINT64_C(1442695040888963407);
    return (double)(rng->state >> 11) * 0x1.0p-53;
}

static void setup(void) {
    xsdir = alea_nuc_xsdir_load(EPRDATA_XSDIR);
    if (!xsdir) return;
    h  = alea_nuc_xsdir_get_nuclide(xsdir, "1000" EPRDATA_SUFFIX);
    pb = alea_nuc_xsdir_get_nuclide(xsdir, "82000" EPRDATA_SUFFIX);
}

/* --- xsdir --- */

TEST(load_xsdir) {
    if (!xsdir) SKIP("eprdata14 not downloaded (make data-eprdata14)");
}

TEST(xsdir_has_100_entries) {
    if (!xsdir) SKIP("no data");
    ASSERT_EQ(alea_nuc_xsdir_count(xsdir), 100);
}

/* --- Hydrogen (Z=1) --- */

TEST(load_hydrogen) {
    if (!xsdir) SKIP("no data");
    ASSERT_NOT_NULL(h);
}

TEST(hydrogen_is_photon) {
    if (!h) SKIP("no data");
    ASSERT_EQ(h->particle, ALEA_NUC_PARTICLE_PHOTON);
}

TEST(hydrogen_energy_grid_larger_than_mcplib84) {
    if (!h) SKIP("no data");
    ASSERT(h->n_energies > 278);
}

TEST(hydrogen_compton_1MeV) {
    if (!h) SKIP("no data");
    double sig_c = alea_nuc_photon_xs_incoherent(h, 1.0);
    ASSERT(sig_c > 0.10 && sig_c < 0.30);
}

TEST(hydrogen_total_1MeV) {
    if (!h) SKIP("no data");
    double sig_t = alea_nuc_xs_total(h, 1.0);
    ASSERT(sig_t > 0.10 && sig_t < 0.35);
}

TEST(hydrogen_coherent_1MeV_small) {
    if (!h) SKIP("no data");
    double sig_coh = alea_nuc_photon_xs_coherent(h, 1.0);
    ASSERT(sig_coh >= 0.0 && sig_coh < 0.01);
}

TEST(hydrogen_pe_10keV) {
    if (!h) SKIP("no data");
    double sig_pe = alea_nuc_photon_xs_photoelectric(h, 0.01);
    ASSERT(sig_pe >= 0.0);
}

TEST(hydrogen_pair_below_threshold) {
    if (!h) SKIP("no data");
    double sig_pp = alea_nuc_photon_xs_pair(h, 0.5);
    ASSERT(sig_pp < 0.001);
}

TEST(hydrogen_heating_1MeV) {
    if (!h) SKIP("no data");
    double heat = alea_nuc_xs_heating(h, 1.0);
    ASSERT(heat > 0.0);
}

/* --- Lead (Z=82) --- */

TEST(load_lead) {
    if (!xsdir) SKIP("no data");
    ASSERT_NOT_NULL(pb);
}

TEST(lead_energy_grid) {
    if (!pb) SKIP("no data");
    ASSERT(pb->n_energies > 0);
}

TEST(lead_total_1MeV) {
    if (!pb) SKIP("no data");
    double sig_t = alea_nuc_xs_total(pb, 1.0);
    ASSERT(sig_t > 5.0 && sig_t < 50.0);
}

TEST(lead_pe_dominates_50keV) {
    if (!pb) SKIP("no data");
    double sig_pe = alea_nuc_photon_xs_photoelectric(pb, 0.05);
    double sig_tot = alea_nuc_xs_total(pb, 0.05);
    ASSERT(sig_tot > 0.0);
    ASSERT(sig_pe / sig_tot > 0.5);
}

TEST(lead_pair_5MeV) {
    if (!pb) SKIP("no data");
    double sig_pp = alea_nuc_photon_xs_pair(pb, 5.0);
    ASSERT(sig_pp > 0.0);
}

TEST(lead_incoherent_ff) {
    if (!pb) SKIP("no data");
    ASSERT_NOT_NULL(pb->photon);
    ASSERT(pb->photon->n_incoherent_ff > 0);
}

TEST(lead_coherent_ff) {
    if (!pb) SKIP("no data");
    ASSERT_NOT_NULL(pb->photon);
    ASSERT(pb->photon->n_coherent_ff > 0);
}

TEST(lead_decodes_detailed_subshell_relaxation_data) {
    if (!pb) SKIP("no data");
    ASSERT_EQ(pb->photon->epr_format, EPRDATA_FORMAT);
    ASSERT(pb->photon->n_subshells > 10);
    int transitions = 0;
    size_t max_photons = 0;
    for (int i = 0; i < pb->photon->n_subshells; i++) {
        const alea_nuc_atomic_subshell_t* shell = &pb->photon->subshells[i];
        ASSERT(shell->designator > 0);
        ASSERT(shell->binding_energy > 0.0);
        transitions += shell->n_transitions;
        if (shell->max_relaxation_photons > max_photons)
            max_photons = shell->max_relaxation_photons;
    }
    ASSERT(transitions > 100);
    ASSERT(max_photons > 0 && max_photons < 256);
}

TEST(lead_decodes_bound_compton_profiles) {
    if (!pb) SKIP("no data");
    ASSERT(pb->photon->n_compton_shells > 10);
    ASSERT(pb->photon->n_compton_profiles > 10);
    ASSERT(pb->photon->n_compton_shells >=
           pb->photon->n_compton_profiles);
    ASSERT_NEAR(pb->photon->compton_shells[
                    pb->photon->n_compton_shells - 1]
                    .cumulative_probability, 1.0, 1e-12);
    for (int i = 0; i < pb->photon->n_compton_profiles; i++) {
        const alea_nuc_compton_profile_t* profile =
            &pb->photon->compton_profiles[i];
        ASSERT_EQ(profile->interpolation, 2);
        ASSERT(profile->n_momenta >= 2);
        ASSERT_NEAR(profile->cdf[0], 0.0, 1e-12);
        ASSERT_NEAR(profile->cdf[profile->n_momenta - 1], 1.0, 1e-12);
    }
}

TEST(lead_subshell_cross_sections_sum_to_photoelectric_total) {
    if (!pb) SKIP("no data");
    const double energies[] = {0.05, 0.1, 1.0};
    for (size_t j = 0; j < sizeof(energies) / sizeof(energies[0]); j++) {
        double sum = 0.0;
        for (int i = 0; i < pb->photon->n_subshells; i++)
            sum += alea_nuc_photon_xs_photoelectric_subshell(
                pb, pb->photon->subshells[i].designator, energies[j]);
        ASSERT_NEAR(sum, alea_nuc_photon_xs_photoelectric(pb, energies[j]),
                    5e-6 * fmax(1.0, sum));
    }
}

TEST(lead_photoelectric_cascades_conserve_energy) {
    if (!pb) SKIP("no data");
    size_t capacity = 0;
    for (int i = 0; i < pb->photon->n_subshells; i++)
        if (pb->photon->subshells[i].max_relaxation_photons > capacity)
            capacity = pb->photon->subshells[i].max_relaxation_photons;
    if (capacity < 2) SKIP("no multi-photon relaxation cascade");
    alea_nuc_particle_state_t* particles = calloc(capacity, sizeof(*particles));
    ASSERT_NOT_NULL(particles);
    alea_nuc_secondary_buffer_t buffer = {particles, capacity, 0};
    alea_nuc_particle_state_t incident = {
        ALEA_NUC_PARTICLE_PHOTON, 0.05, {0.0, 0.0, 1.0}, 0.75, 2.0
    };
    test_rng_t rng = {UINT64_C(0x123456789abcdef0)};
    int photoelectric_events = 0, relaxation_photons = 0;
    for (int history = 0; history < 1000; history++) {
        buffer.count = 0;
        alea_nuc_collision_result_t result;
        ASSERT_EQ(alea_nuc_sample_photon_collision_with_secondaries(
                      pb, &incident, test_rng, &rng, &buffer, &result), ALEA_OK);
        if (result.mt != 522) continue;
        photoelectric_events++;
        double photon_energy = 0.0;
        ASSERT_EQ(buffer.count, result.n_emitted);
        for (size_t i = 0; i < buffer.count; i++) {
            photon_energy += buffer.particles[i].energy;
            ASSERT_EQ(buffer.particles[i].type, ALEA_NUC_PARTICLE_PHOTON);
            ASSERT_EQ(buffer.particles[i].weight, incident.weight);
            ASSERT_EQ(buffer.particles[i].time, incident.time);
        }
        relaxation_photons += (int)buffer.count;
        ASSERT_NEAR(photon_energy + result.local_energy_deposition,
                    incident.energy, 2e-13);
    }
    free(particles);
    ASSERT(photoelectric_events > 500);
    ASSERT(relaxation_photons > 0);
}

TEST(lead_bound_compton_events_conserve_energy) {
    if (!pb) SKIP("no data");
    alea_nuc_particle_state_t incident = {
        ALEA_NUC_PARTICLE_PHOTON, 1.0, {0.0, 0.0, 1.0}, 0.75, 2.0
    };
    size_t capacity = alea_nuc_photon_secondary_capacity(pb, incident.energy);
    ASSERT(capacity > 0 && capacity < 256);
    alea_nuc_particle_state_t* particles = calloc(capacity, sizeof(*particles));
    ASSERT_NOT_NULL(particles);
    alea_nuc_secondary_buffer_t buffer = {particles, capacity, 0};
    test_rng_t rng = {UINT64_C(0x6a09e667f3bcc909)};
    int compton_events = 0, broadened_events = 0, relaxation_photons = 0;
    for (int history = 0; history < 2000; history++) {
        buffer.count = 0;
        alea_nuc_collision_result_t result;
        ASSERT_EQ(alea_nuc_sample_photon_collision_with_secondaries(
                      pb, &incident, test_rng, &rng, &buffer, &result), ALEA_OK);
        if (result.mt != 504) continue;
        compton_events++;
        double banked = 0.0;
        for (size_t i = 0; i < buffer.count; i++)
            banked += buffer.particles[i].energy;
        relaxation_photons += (int)buffer.count;
        double free_energy = incident.energy /
            (1.0 + incident.energy / 0.51099895069 * (1.0 - result.mu_lab));
        if (fabs(result.outgoing.energy - free_energy) > 1e-10)
            broadened_events++;
        ASSERT_NEAR(result.outgoing.energy + banked +
                    result.local_energy_deposition, incident.energy, 2e-10);
        ASSERT_EQ(buffer.count, result.n_emitted);
    }
    free(particles);
    ASSERT(compton_events > 100);
    ASSERT(broadened_events > 100);
    ASSERT(relaxation_photons > 0);
}

TEST(lead_heating_per_collision_1MeV) {
    if (!pb) SKIP("no data");
    double hpc = alea_nuc_heating_per_collision(pb, 1.0);
    ASSERT(hpc > 0.0 && hpc < 1.0);
}

/* --- Bulk load all 100 elements --- */

TEST(load_all_100_elements) {
    if (!xsdir) SKIP("no data");
    int loaded = 0;
    for (int z = 1; z <= 100; z++) {
        char zaid[24];
        snprintf(zaid, sizeof(zaid), "%d000%s", z, EPRDATA_SUFFIX);
        alea_nuc_nuclide_t* nuc = alea_nuc_xsdir_get_nuclide(xsdir, zaid);
        if (!nuc || !nuc->photon || nuc->n_energies <= 0) continue;
        double sig = alea_nuc_xs_total(nuc, 1.0);
        if (sig > 0.0) loaded++;
    }
    ASSERT_EQ(loaded, 100);
}

/* --- Runner --- */

alea_test_entry_t *alea_test_list = NULL;
alea_test_entry_t **alea_test_tail = &alea_test_list;
int alea_test_passed = 0;
int alea_test_failed = 0;
int alea_test_skipped = 0;
int alea_test_current_failed = 0;
int alea_test_current_skipped = 0;
const char *alea_test_current_name = NULL;

int main(int argc, char **argv) {
    const char *filter = argc > 1 ? argv[1] : NULL;

    printf("nucdata %s photon tests\n", EPRDATA_LABEL);
    printf("===============================\n");

    setup();

    for (alea_test_entry_t *t = alea_test_list; t; t = t->next) {
        if (filter && strstr(t->name, filter) == NULL) continue;
        alea_test_current_failed = 0;
        alea_test_current_skipped = 0;
        alea_test_current_name = t->name;
        printf("  %-50s ", t->name);
        fflush(stdout);
        t->fn();
        if (alea_test_current_skipped) {
            alea_test_skipped++;
        } else if (alea_test_current_failed) {
            alea_test_failed++;
        } else {
            printf("OK\n");
            alea_test_passed++;
        }
    }

    printf("\n----------------------------------------\n");
    printf("Results: %d passed, %d failed, %d skipped\n",
           alea_test_passed, alea_test_failed, alea_test_skipped);
    printf("----------------------------------------\n\n");

    alea_nuc_xsdir_free(xsdir);
    return alea_test_failed > 0 ? 1 : 0;
}
