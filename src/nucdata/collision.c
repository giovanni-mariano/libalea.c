// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/**
 * @file collision.c
 * @brief Prepared stationary-target neutron collision physics
 */

#include "nuclear_internal.h"
#include "core/alea_materials.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static alea_error_t draw_uniform(alea_nuc_random_fn random, void* context,
                                 double* value);

typedef struct {
    const alea_nuc_mat_component_t* source;
    const alea_nuc_thermal_t* thermal;
    int* event_reactions;
    int n_event_reactions;
} prepared_component_t;

struct alea_nuc_prepared_material {
    const alea_nuc_material_t* source;
    prepared_component_t* components;
    int n_components;
    uint32_t capabilities;
    uint32_t requested_capabilities;
    alea_nuc_particle_t particle;
};

static void report_clear(alea_nuc_capability_report_t* report) {
    if (!report) return;
    memset(report, 0, sizeof(*report));
    report->issue = ALEA_NUC_PREP_OK;
    report->component_index = -1;
}

static alea_error_t report_failure(alea_nuc_capability_report_t* report,
                                   alea_nuc_prepare_issue_t issue,
                                   uint32_t missing, int component, int mt,
                                   const char* detail) {
    if (report) {
        report->issue = issue;
        report->missing_capabilities = missing;
        report->component_index = component;
        report->mt = mt;
        snprintf(report->detail, sizeof(report->detail), "%s", detail);
    }
    return ALEA_ERR_UNSUPPORTED;
}

static int reaction_has_positive_xs(const alea_nuc_reaction_t* reaction) {
    if (!reaction || !reaction->xs || reaction->n_energies <= 0) return 0;
    for (int i = 0; i < reaction->n_energies; i++)
        if (isfinite(reaction->xs[i]) && reaction->xs[i] > 0.0) return 1;
    return 0;
}

static int is_naturally_occurring_isotope(int z, int a) {
    const alea_element_t* element = alea_get_element(z);
    if (!element) return 0;
    for (size_t i = 0; i < element->isotope_count; i++)
        if (element->isotopes[i].mass_number == a &&
            element->isotopes[i].abundance > 0.0) return 1;
    return 0;
}

static int reaction_is_supported_absorption_mt(int mt) {
    if (mt == 102 || (mt >= 103 && mt <= 109) ||
        (mt >= 111 && mt <= 117)) return 1;
    return mt >= 600 && mt <= 849;
}

static int aggregate_mt_for(int mt) {
    if (mt >= 600 && mt <= 649) return 103;
    if (mt >= 650 && mt <= 699) return 104;
    if (mt >= 700 && mt <= 749) return 105;
    if (mt >= 750 && mt <= 799) return 106;
    if (mt >= 800 && mt <= 849) return 107;
    return 0;
}

static int nuclide_has_active_mt(const alea_nuc_nuclide_t* nuc, int mt) {
    for (int i = 0; i < nuc->n_reactions; i++)
        if (nuc->reactions[i].mt == mt &&
            reaction_has_positive_xs(&nuc->reactions[i])) return 1;
    return 0;
}

static int reaction_is_redundant(const alea_nuc_nuclide_t* nuc,
                                 const alea_nuc_reaction_t* reaction) {
    /* These are derived production, heating, and damage responses rather
     * than mutually exclusive collision channels. */
    if (reaction->mt == 3 ||
        (reaction->mt >= 201 && reaction->mt <= 207) ||
        reaction->mt == 301 || reaction->mt == 444) return 1;
    if (reaction->mt == 4) {
        for (int mt = 51; mt <= 91; mt++)
            if (nuclide_has_active_mt(nuc, mt)) return 1;
    }
    if ((reaction->mt == 19 || reaction->mt == 20 || reaction->mt == 21 ||
         reaction->mt == 38) && nuclide_has_active_mt(nuc, 18)) return 1;
    int aggregate = aggregate_mt_for(reaction->mt);
    return aggregate != 0 && nuclide_has_active_mt(nuc, aggregate);
}

static int composite_neutron_charged_particle_mt(int mt) {
    return (mt >= 22 && mt <= 25) || (mt >= 28 && mt <= 30) ||
           (mt >= 32 && mt <= 36);
}

static int photon_parent_maps_to_event(
    const alea_nuc_nuclide_t* nuc,
    const alea_nuc_photon_production_t* production, int reaction_mt) {
    int parent_mt = production->parent_mt;
    if (aggregate_mt_for(parent_mt) == reaction_mt) return 1;
    return production->production_xs && reaction_mt == 5 &&
           composite_neutron_charged_particle_mt(parent_mt) &&
           !nuclide_has_active_mt(nuc, parent_mt);
}

static int photon_production_applies(
    const alea_nuc_nuclide_t* nuc,
    const alea_nuc_photon_production_t* production, int reaction_mt) {
    int parent = production->parent_mt;
    if (parent == reaction_mt) return 1;
    if (photon_parent_maps_to_event(nuc, production, reaction_mt)) return 1;
    if (parent == 3) return reaction_mt != 2;
    if (parent == 4) return reaction_mt >= 51 && reaction_mt <= 91;
    if (parent == 18)
        return reaction_mt == 19 || reaction_mt == 20 ||
               reaction_mt == 21 || reaction_mt == 38;
    return 0;
}

static double photon_yield_for_event(
    const alea_nuc_nuclide_t* nuc,
    const alea_nuc_photon_production_t* production,
    int reaction_mt, double energy) {
    double yield = alea_nuc_photon_production_yield(
        nuc, production, energy);
    if (!photon_parent_maps_to_event(
            nuc, production, reaction_mt))
        return yield;
    return alea_nuc_photon_production_event_yield(
        nuc, production, reaction_mt, energy);
}

static double reaction_xs(const alea_nuc_nuclide_t* nuc,
                          const alea_nuc_reaction_t* reaction,
                          double energy) {
    double f;
    int ie = alea_nuc_energy_lookup_trusted(nuc->energy, nuc->n_energies,
                                            energy, &f);
    if (ie < 0) return 0.0;
    int start = reaction->threshold_index - 1;
    int ri = ie - start;
    if (ri < 0 || !reaction->xs || reaction->n_energies <= 0) return 0.0;
    if (ri >= reaction->n_energies - 1)
        return ri == reaction->n_energies - 1
            ? reaction->xs[reaction->n_energies - 1] : 0.0;
    return reaction->xs[ri] +
           f * (reaction->xs[ri + 1] - reaction->xs[ri]);
}

static int validate_angular_point(const alea_nuc_angular_point_t* point) {
    if (point->type == ALEA_NUC_ANG_ISOTROPIC) return 1;
    if (point->type == ALEA_NUC_ANG_EQUIPROBABLE) {
        if (point->n_cosines != 33 || !point->cosine) return 0;
        for (int i = 0; i < 33; i++) {
            if (!isfinite(point->cosine[i]) || point->cosine[i] < -1.0 ||
                point->cosine[i] > 1.0 ||
                (i > 0 && point->cosine[i] < point->cosine[i - 1])) return 0;
        }
        return 1;
    }
    if (point->type != ALEA_NUC_ANG_TABULAR || point->n_cosines < 2 ||
        !point->cosine || !point->pdf || !point->cdf ||
        (point->interpolation != 1 && point->interpolation != 2)) return 0;

    int n = point->n_cosines;
    for (int i = 0; i < n; i++) {
        if (point->cosine[i] < -1.0 || point->cosine[i] > 1.0) return 0;
    }
    return alea_nuc_tabular_pdf_valid(point->cosine, point->pdf, point->cdf,
                                      n, point->interpolation, 0);
}

int alea_nuc_validate_angular_internal(
    const alea_nuc_angular_dist_t* angular) {
    if (!angular) return 1; /* absent LAND means isotropic */
    if (angular->n_energies <= 0 || !angular->energy || !angular->data)
        return 0;
    for (int i = 0; i < angular->n_energies; i++) {
        if (!isfinite(angular->energy[i]) ||
            (i > 0 && angular->energy[i] < angular->energy[i - 1]) ||
            !validate_angular_point(&angular->data[i])) return 0;
    }
    return 1;
}

static int validate_urr(const alea_nuc_urr_t* urr) {
    if (!urr) return 1;
    if (urr->n_energies < 2 || urr->n_bands <= 0 || !urr->energy ||
        !urr->table || (urr->interp != 2 && urr->interp != 5)) return 0;
    int n = urr->n_energies, m = urr->n_bands;
    for (int i = 0; i < n; i++) {
        if (!isfinite(urr->energy[i]) || urr->energy[i] <= 0.0 ||
            (i > 0 && urr->energy[i] <= urr->energy[i - 1])) return 0;
        const double* row = &urr->table[i * 6 * m];
        for (int j = 0; j < m; j++) {
            if (!isfinite(row[j]) || row[j] < 0.0 || row[j] > 1.0 ||
                (j > 0 && row[j] < row[j - 1])) return 0;
            for (int q = 1; q < 6; q++)
                if (!isfinite(row[q * m + j]))
                    return 0;
        }
        if (fabs(row[m - 1] - 1.0) > 1e-6) return 0;
    }
    return 1;
}

static int validate_delayed(const alea_nuc_nuclide_t* nuc) {
    if (!nuc->fission) return 1;
    int n = nuc->fission->n_delayed_groups;
    if (n == 0)
        return nuc->fission->delayed == NULL &&
               nuc->fission->delayed_groups == NULL;
    if (!nuc->fission->delayed || !nuc->fission->delayed_groups) return 0;
    for (int g = 0; g < n; g++) {
        const alea_nuc_delayed_group_t* group =
            &nuc->fission->delayed_groups[g];
        if (!isfinite(group->decay_rate) || group->decay_rate <= 0.0 ||
            !alea_nuc_interp_regions_valid(group->nbt, group->interp,
                                           group->n_regions,
                                           group->n_energies) ||
            !group->energy || !group->probability ||
            alea_nuc_energy_dist_validate(group->spectrum, NULL) != ALEA_OK)
            return 0;
        for (int i = 0; i < group->n_energies; i++)
            if (!isfinite(group->energy[i]) ||
                !isfinite(group->probability[i]) ||
                group->probability[i] < 0.0 ||
                (i > 0 && group->energy[i] < group->energy[i - 1]))
                return 0;
    }
    for (int i = 0; i < nuc->n_energies; i++) {
        double sum = 0.0;
        for (int g = 0; g < n; g++) {
            const alea_nuc_delayed_group_t* group =
                &nuc->fission->delayed_groups[g];
            double probability;
            if (alea_nuc_interp_eval(
                    group->energy, group->probability, group->n_energies,
                    group->nbt, group->interp, group->n_regions,
                    nuc->energy[i], &probability) != ALEA_OK)
                return 0;
            sum += probability;
        }
        if (!(sum > 0.0) || !isfinite(sum)) return 0;
    }
    return 1;
}

static alea_error_t inspect_photon(const alea_nuc_nuclide_t* nuc,
                                   int component,
                                   alea_nuc_capability_report_t* report) {
    const alea_nuc_photon_data_t* ph = nuc ? nuc->photon : NULL;
    if (!ph || ph->n_energies < 2 || !ph->energy || !ph->ln_energy ||
        !ph->ln_sigma_coherent || !ph->ln_sigma_incoherent ||
        !ph->ln_sigma_photoelectric || !ph->ln_sigma_pair)
        return report_failure(report, ALEA_NUC_PREP_INVALID_CROSS_SECTIONS,
            ALEA_NUC_CAP_PHOTON, component, 0,
            "photoatomic cross-section tables are incomplete");
    int coherent = 0, incoherent = 0;
    for (int i = 0; i < ph->n_energies; i++) {
        if (!isfinite(ph->energy[i]) || ph->energy[i] <= 0.0 ||
            (i > 0 && ph->energy[i] < ph->energy[i - 1]))
            return report_failure(report,
                ALEA_NUC_PREP_INVALID_CROSS_SECTIONS, ALEA_NUC_CAP_PHOTON,
                component, 0, "photoatomic energy grid is invalid");
        coherent |= ph->ln_sigma_coherent[i] > -1e30;
        incoherent |= ph->ln_sigma_incoherent[i] > -1e30;
    }
    if ((coherent && (ph->n_coherent_ff < 2 ||
                      !ph->coherent_momentum ||
                      !ph->coherent_ff_cumulative)) ||
        (incoherent && (ph->n_incoherent_ff < 2 ||
                        !ph->incoherent_momentum || !ph->incoherent_ff)))
        return report_failure(report, ALEA_NUC_PREP_INVALID_ANGULAR,
            ALEA_NUC_CAP_PHOTON, component, coherent ? 502 : 504,
            "photoatomic scattering factors are incomplete");
    if (coherent)
        for (int i = 0; i < ph->n_coherent_ff; i++)
            if (!isfinite(ph->coherent_momentum[i]) ||
                !isfinite(ph->coherent_ff_cumulative[i]) ||
                ph->coherent_momentum[i] < 0.0 ||
                ph->coherent_ff_cumulative[i] < 0.0 ||
                (i > 0 && (ph->coherent_momentum[i] <=
                               ph->coherent_momentum[i - 1] ||
                           ph->coherent_ff_cumulative[i] <
                               ph->coherent_ff_cumulative[i - 1])))
                return report_failure(report, ALEA_NUC_PREP_INVALID_ANGULAR,
                    ALEA_NUC_CAP_PHOTON, component, 502,
                    "coherent form-factor integral is invalid");
    if (incoherent)
        for (int i = 0; i < ph->n_incoherent_ff; i++)
            if (!isfinite(ph->incoherent_momentum[i]) ||
                !isfinite(ph->incoherent_ff[i]) ||
                ph->incoherent_momentum[i] < 0.0 ||
                ph->incoherent_ff[i] < 0.0 ||
                (i > 0 && ph->incoherent_momentum[i] <=
                              ph->incoherent_momentum[i - 1]))
                return report_failure(report, ALEA_NUC_PREP_INVALID_ANGULAR,
                    ALEA_NUC_CAP_PHOTON, component, 504,
                    "incoherent scattering function is invalid");
    if (report) report->available_capabilities |= ALEA_NUC_CAP_PHOTON;
    return ALEA_OK;
}

static alea_error_t inspect_nuclide(const alea_nuc_nuclide_t* nuc,
                                    int component,
                                    uint32_t requested,
                                    alea_nuc_capability_report_t* report) {
    if (nuc && nuc->particle == ALEA_NUC_PARTICLE_PHOTON)
        return inspect_photon(nuc, component, report);
    if (!nuc || nuc->particle != ALEA_NUC_PARTICLE_NEUTRON)
        return report_failure(report, ALEA_NUC_PREP_UNSUPPORTED_PARTICLE,
            ALEA_NUC_CAP_RESTRICTED_NEUTRON, component, 0,
            "continuous collision preparation accepts neutron tables only");
    if (!validate_urr(nuc->urr))
        return report_failure(report, ALEA_NUC_PREP_UNSUPPORTED_URR,
            ALEA_NUC_CAP_URR, component, 0,
            "unresolved-resonance probability table is invalid");
    if (!validate_delayed(nuc))
        return report_failure(report,
            ALEA_NUC_PREP_INVALID_ENERGY_DISTRIBUTION,
            ALEA_NUC_CAP_DELAYED_NEUTRON, component, 18,
            "delayed-neutron group data is invalid");
    if (nuc->n_energies < 2 || !nuc->energy || !nuc->sigma_total ||
        !nuc->sigma_elastic || !nuc->sigma_abs || !isfinite(nuc->awr) ||
        nuc->awr <= 0.0)
        return report_failure(report, ALEA_NUC_PREP_INVALID_CROSS_SECTIONS,
            ALEA_NUC_CAP_RESTRICTED_NEUTRON, component, 0,
            "principal neutron cross sections or AWR are incomplete");
    if (nuc->n_reactions < 0 || (nuc->n_reactions > 0 && !nuc->reactions))
        return report_failure(report, ALEA_NUC_PREP_INVALID_CROSS_SECTIONS,
            ALEA_NUC_CAP_RESTRICTED_NEUTRON, component, 0,
            "reaction table is incomplete");
    if (!alea_nuc_validate_angular_internal(nuc->elastic_angular))
        return report_failure(report, ALEA_NUC_PREP_INVALID_ANGULAR,
            ALEA_NUC_CAP_STATIONARY_ELASTIC, component, 2,
            "elastic angular distribution is not sampleable");

    if ((requested & ALEA_NUC_CAP_PHOTON_PRODUCTION) &&
        (nuc->n_photon_productions <= 0 || !nuc->photon_productions))
        return report_failure(report, ALEA_NUC_PREP_UNSUPPORTED_CAPABILITY,
            ALEA_NUC_CAP_PHOTON_PRODUCTION, component, 0,
            "neutron table has no photon-production channels");
    for (int p = 0; p < nuc->n_photon_productions &&
                    (requested & ALEA_NUC_CAP_PHOTON_PRODUCTION); p++) {
        const alea_nuc_photon_production_t* production =
            &nuc->photon_productions[p];
        int yield_valid = production->n_energies > 0 && production->values;
        if (production->production_xs) {
            yield_valid = yield_valid && production->mf == 13 &&
                production->threshold_index >= 1 &&
                production->threshold_index <= nuc->n_energies &&
                production->n_energies <=
                    nuc->n_energies - production->threshold_index + 1;
        } else {
            yield_valid = yield_valid &&
                (production->mf == 12 || production->mf == 16) &&
                production->n_energies >= 2 && production->energy &&
                alea_nuc_interp_regions_valid(
                    production->nbt, production->interp,
                    production->n_regions, production->n_energies);
        }
        for (int j = 0; yield_valid && j < production->n_energies; j++) {
            if (!isfinite(production->values[j]) ||
                production->values[j] < 0.0 ||
                (!production->production_xs &&
                 (!isfinite(production->energy[j]) ||
                  (j > 0 && production->energy[j] <
                                production->energy[j - 1]))))
                yield_valid = 0;
        }
        if (!yield_valid)
            return report_failure(report,
                ALEA_NUC_PREP_INVALID_CROSS_SECTIONS,
                ALEA_NUC_CAP_PHOTON_PRODUCTION, component,
                production->parent_mt,
                "photon-production yield or cross section is invalid");
        int selectable = 0;
        for (int r = 0; r < nuc->n_reactions; r++)
            if (photon_production_applies(
                    nuc, production, nuc->reactions[r].mt) &&
                reaction_has_positive_xs(&nuc->reactions[r]) &&
                !reaction_is_redundant(nuc, &nuc->reactions[r]))
                selectable = 1;
        if (!selectable)
            return report_failure(report,
                ALEA_NUC_PREP_UNSUPPORTED_REACTION,
                ALEA_NUC_CAP_PHOTON_PRODUCTION, component,
                production->parent_mt,
                "photon production is attached to a nonselectable reaction");
        int unsupported_law = 0;
        if (!alea_nuc_validate_angular_internal(production->angular))
            return report_failure(report, ALEA_NUC_PREP_INVALID_ANGULAR,
                ALEA_NUC_CAP_PHOTON_PRODUCTION, component,
                production->parent_mt,
                "photon-production angular distribution is not sampleable");
        alea_error_t status = alea_nuc_energy_dist_validate(
            production->spectrum, &unsupported_law);
        if (status != ALEA_OK) {
            if (report) report->law = unsupported_law;
            return report_failure(report,
                ALEA_NUC_PREP_INVALID_ENERGY_DISTRIBUTION,
                ALEA_NUC_CAP_PHOTON_PRODUCTION, component,
                production->parent_mt,
                status == ALEA_ERR_UNSUPPORTED
                    ? "photon-production energy law is unsupported"
                    : "photon-production energy distribution is invalid");
        }
    }
    if (nuc->n_photon_productions > 0 &&
        (requested & ALEA_NUC_CAP_PHOTON_PRODUCTION))
        report->available_capabilities |= ALEA_NUC_CAP_PHOTON_PRODUCTION;

    for (int r = 0; r < nuc->n_reactions; r++) {
        const alea_nuc_reaction_t* reaction = &nuc->reactions[r];
        int maximum = nuc->n_energies - reaction->threshold_index + 1;
        if (reaction->threshold_index < 1 ||
            reaction->threshold_index > nuc->n_energies ||
            reaction->n_energies < 0 || reaction->n_energies > maximum ||
            (reaction->n_energies > 0 && !reaction->xs))
            return report_failure(report, ALEA_NUC_PREP_INVALID_CROSS_SECTIONS,
                ALEA_NUC_CAP_RESTRICTED_NEUTRON, component, reaction->mt,
                "reaction threshold or cross-section extent is invalid");
        for (int j = 0; j < reaction->n_energies; j++) {
            if (!isfinite(reaction->xs[j]) || reaction->xs[j] < 0.0)
                return report_failure(report, ALEA_NUC_PREP_INVALID_CROSS_SECTIONS,
                    ALEA_NUC_CAP_RESTRICTED_NEUTRON, component, reaction->mt,
                    "reaction cross sections must be finite and nonnegative");
        }
        if (!reaction_has_positive_xs(reaction)) continue;
        if (reaction_is_redundant(nuc, reaction)) continue;
        if (reaction->ty == 0) {
            if (reaction_is_supported_absorption_mt(reaction->mt)) continue;
            return report_failure(report, ALEA_NUC_PREP_UNSUPPORTED_REACTION,
                ALEA_NUC_CAP_ABSORPTION, component, reaction->mt,
                "zero-yield reaction is not a supported absorption channel");
        }
        if (!alea_nuc_validate_angular_internal(reaction->angular))
            return report_failure(report, ALEA_NUC_PREP_INVALID_ANGULAR,
                ALEA_NUC_CAP_NEUTRON_EMISSION, component, reaction->mt,
                "neutron-emission angular distribution is not sampleable");
        int unsupported_law = 0;
        alea_error_t distribution_status = alea_nuc_energy_dist_validate(
            reaction->energy, &unsupported_law);
        if (distribution_status != ALEA_OK) {
            if (report) report->law = unsupported_law;
            return report_failure(report,
                ALEA_NUC_PREP_INVALID_ENERGY_DISTRIBUTION,
                ALEA_NUC_CAP_NEUTRON_EMISSION, component, reaction->mt,
                distribution_status == ALEA_ERR_UNSUPPORTED
                    ? "neutron-emission energy law is not implemented"
                    : "neutron-emission energy distribution is invalid");
        }
        for (int i = 0; i < nuc->n_energies; i++) {
            double yield = alea_nuc_reaction_yield(nuc, reaction->mt,
                                                   nuc->energy[i]);
            if (!isfinite(yield) || yield < 0.0)
                return report_failure(report,
                    ALEA_NUC_PREP_UNSUPPORTED_REACTION,
                    ALEA_NUC_CAP_NEUTRON_EMISSION, component, reaction->mt,
                    "neutron multiplicity is invalid");
        }
        if (report) {
            report->available_capabilities |= ALEA_NUC_CAP_NEUTRON_EMISSION;
            if (reaction->mt == 18 || reaction->mt == 19 ||
                reaction->mt == 20 || reaction->mt == 21 || reaction->mt == 38)
                report->available_capabilities |= ALEA_NUC_CAP_FISSION;
        }
    }

    for (int i = 0; i < nuc->n_energies; i++) {
        double total = nuc->sigma_total[i];
        double elastic = nuc->sigma_elastic[i];
        double absorption = nuc->sigma_abs[i];
        if (!isfinite(nuc->energy[i]) || !isfinite(total) ||
            !isfinite(elastic) || !isfinite(absorption) || total < 0.0 ||
            elastic < 0.0 || absorption < 0.0 ||
            (i > 0 && nuc->energy[i] < nuc->energy[i - 1]))
            return report_failure(report, ALEA_NUC_PREP_INVALID_CROSS_SECTIONS,
                ALEA_NUC_CAP_RESTRICTED_NEUTRON, component, 0,
                "cross-section grids must be finite, nonnegative and nondecreasing");

        double scale = fmax(1.0, total);
        double reaction_total = 0.0;
        for (int r = 0; r < nuc->n_reactions; r++) {
            const alea_nuc_reaction_t* reaction = &nuc->reactions[r];
            if (!reaction_is_redundant(nuc, reaction))
                reaction_total += reaction_xs(nuc, reaction, nuc->energy[i]);
        }
        if (fabs(total - elastic - reaction_total) > 1e-6 * scale)
            return report_failure(report, ALEA_NUC_PREP_INVALID_CROSS_SECTIONS,
                ALEA_NUC_CAP_RESTRICTED_NEUTRON, component, 0,
                "total cross section is inconsistent with selectable event channels");

        double reaction_absorption = 0.0;
        for (int r = 0; r < nuc->n_reactions; r++) {
            const alea_nuc_reaction_t* reaction = &nuc->reactions[r];
            if (reaction_is_supported_absorption_mt(reaction->mt) &&
                reaction->ty == 0 && !reaction_is_redundant(nuc, reaction))
                reaction_absorption += reaction_xs(nuc, reaction, nuc->energy[i]);
        }
        if (fabs(absorption - reaction_absorption) > 1e-6 * scale)
            return report_failure(report, ALEA_NUC_PREP_INVALID_CROSS_SECTIONS,
                ALEA_NUC_CAP_ABSORPTION, component, 0,
                "absorption cross section is inconsistent with event channels");
    }

    if (report) {
        report->available_capabilities |= ALEA_NUC_CAP_RESTRICTED_NEUTRON;
        if (nuc->temperature > 0.0 && isfinite(nuc->temperature))
            report->available_capabilities |= ALEA_NUC_CAP_FREE_GAS;
        if (nuc->urr) report->available_capabilities |= ALEA_NUC_CAP_URR;
        if (nuc->fission && nuc->fission->n_delayed_groups > 0)
            report->available_capabilities |= ALEA_NUC_CAP_DELAYED_NEUTRON;
    }
    return ALEA_OK;
}

alea_error_t alea_nuc_capabilities(const alea_nuc_nuclide_t* nuc,
                                   alea_nuc_capability_report_t* report) {
    if (!nuc || !report) return ALEA_ERR_NULL_ARG;
    report_clear(report);
    uint32_t requested = UINT32_MAX;
    if (nuc->n_photon_productions <= 0)
        requested &= ~ALEA_NUC_CAP_PHOTON_PRODUCTION;
    return inspect_nuclide(nuc, -1, requested, report);
}

static void prepared_free_components(alea_nuc_prepared_material_t* prepared) {
    if (!prepared) return;
    for (int i = 0; i < prepared->n_components; i++)
        free(prepared->components[i].event_reactions);
    free(prepared->components);
}

void alea_nuc_prepared_material_free(alea_nuc_prepared_material_t* prepared) {
    if (!prepared) return;
    prepared_free_components(prepared);
    free(prepared);
}

alea_error_t alea_nuc_prepare_material(
    const alea_nuc_material_t* material,
    const alea_nuc_prepare_requirements_t* requirements,
    alea_nuc_capability_report_t* report,
    alea_nuc_prepared_material_t** output) {
    if (!material || !requirements || !report || !output)
        return ALEA_ERR_NULL_ARG;
    *output = NULL;
    report_clear(report);

    uint32_t supported = ALEA_NUC_CAP_CONTINUOUS_NEUTRON |
                         ALEA_NUC_CAP_FISSION | ALEA_NUC_CAP_URR |
                         ALEA_NUC_CAP_DELAYED_NEUTRON |
                         ALEA_NUC_CAP_FREE_GAS | ALEA_NUC_CAP_THERMAL_SAB |
                         ALEA_NUC_CAP_PHOTON |
                         ALEA_NUC_CAP_PHOTON_PRODUCTION;
    uint32_t unknown = requirements->required_capabilities & ~supported;
    if (unknown)
        return report_failure(report, ALEA_NUC_PREP_UNSUPPORTED_CAPABILITY,
            unknown, -1, 0, "requested collision capability is not implemented");
    if (material->n_components <= 0 || !material->components)
        return report_failure(report, ALEA_NUC_PREP_EMPTY_MATERIAL,
            requirements->required_capabilities, -1, 0,
            "material has no nuclide components");
    if (!isfinite(requirements->thermal_temperature_tolerance) ||
        requirements->thermal_temperature_tolerance < 0.0 ||
        (requirements->n_thermal_associations > 0 &&
         !requirements->thermal_associations))
        return report_failure(report,
            ALEA_NUC_PREP_INVALID_THERMAL_ASSOCIATION,
            ALEA_NUC_CAP_THERMAL_SAB, -1, 0,
            "thermal associations or temperature tolerance are invalid");

    alea_nuc_prepared_material_t* prepared = alea_nuc_calloc(1, sizeof(*prepared));
    if (!prepared) return ALEA_ERR_OUT_OF_MEMORY;
    prepared->components = alea_nuc_calloc((size_t)material->n_components,
                                  sizeof(*prepared->components));
    if (!prepared->components) {
        free(prepared);
        return ALEA_ERR_OUT_OF_MEMORY;
    }
    prepared->source = material;
    prepared->n_components = material->n_components;
    prepared->capabilities = 0;
    if (!material->components[0].nuclide) {
        alea_nuc_prepared_material_free(prepared);
        return report_failure(report, ALEA_NUC_PREP_UNSUPPORTED_PARTICLE,
            requirements->required_capabilities, 0, 0,
            "material component has no nuclear-data table");
    }
    prepared->particle = material->components[0].nuclide->particle;

    int populated_components = 0;
    for (int i = 0; i < material->n_components; i++) {
        const alea_nuc_mat_component_t* source = &material->components[i];
        if (!source->nuclide ||
            source->nuclide->particle != prepared->particle) {
            alea_nuc_prepared_material_free(prepared);
            return report_failure(report, ALEA_NUC_PREP_UNSUPPORTED_PARTICLE,
                requirements->required_capabilities, i, 0,
                "prepared materials cannot mix particle table types");
        }
        if (!isfinite(source->number_density) || source->number_density < 0.0) {
            prepared_free_components(prepared);
            free(prepared);
            return report_failure(report, ALEA_NUC_PREP_INVALID_CROSS_SECTIONS,
                requirements->required_capabilities, i, 0,
                "component number density is invalid");
        }
        if (source->number_density > 0.0) populated_components++;
        alea_error_t err = inspect_nuclide(
            source->nuclide, i, requirements->required_capabilities, report);
        if (err != ALEA_OK) {
            prepared_free_components(prepared);
            free(prepared);
            return err;
        }
        if ((requirements->required_capabilities & ALEA_NUC_CAP_FREE_GAS) &&
            (!(source->nuclide->temperature > 0.0) ||
             !isfinite(source->nuclide->temperature))) {
            prepared_free_components(prepared);
            free(prepared);
            return report_failure(report,
                ALEA_NUC_PREP_UNSUPPORTED_CAPABILITY,
                ALEA_NUC_CAP_FREE_GAS, i, 2,
                "free-gas scattering requires a positive table temperature");
        }
        prepared_component_t* component = &prepared->components[i];
        component->source = source;
        int maximum = source->nuclide->n_reactions;
        if (maximum > 0) {
            component->event_reactions = alea_nuc_malloc((size_t)maximum * sizeof(int));
            if (!component->event_reactions) {
                alea_nuc_prepared_material_free(prepared);
                return ALEA_ERR_OUT_OF_MEMORY;
            }
        }
        for (int r = 0; r < maximum; r++) {
            const alea_nuc_reaction_t* reaction = &source->nuclide->reactions[r];
            if (reaction_has_positive_xs(reaction) &&
                !reaction_is_redundant(source->nuclide, reaction))
                component->event_reactions[component->n_event_reactions++] = r;
        }
    }

    if (populated_components == 0) {
        alea_nuc_prepared_material_free(prepared);
        return report_failure(report, ALEA_NUC_PREP_EMPTY_MATERIAL,
            requirements->required_capabilities, -1, 0,
            "material has no positive number density");
    }

    if (requirements->required_capabilities & ALEA_NUC_CAP_THERMAL_SAB) {
        int populated_associations = 0;
        for (size_t i = 0; i < requirements->n_thermal_associations; i++) {
            const alea_nuc_thermal_association_t* association =
                &requirements->thermal_associations[i];
            int index = association->component_index;
            if (index < 0 || index >= material->n_components ||
                !association->thermal ||
                prepared->components[index].thermal) {
                alea_nuc_prepared_material_free(prepared);
                return report_failure(report,
                    ALEA_NUC_PREP_INVALID_THERMAL_ASSOCIATION,
                    ALEA_NUC_CAP_THERMAL_SAB, index, 0,
                    "thermal association has an invalid or duplicate component");
            }
            const alea_nuc_thermal_t* thermal = association->thermal;
            const alea_nuc_nuclide_t* nuc =
                material->components[index].nuclide;
            if (prepared->particle != ALEA_NUC_PARTICLE_NEUTRON ||
                !alea_nuc_thermal_validate_internal(thermal)) {
                alea_nuc_prepared_material_free(prepared);
                return report_failure(report,
                    ALEA_NUC_PREP_INVALID_THERMAL_ASSOCIATION,
                    ALEA_NUC_CAP_THERMAL_SAB, index, 0,
                    "thermal table is invalid or is not associated with neutron data");
            }
            int applicable = 0;
            int64_t nuclide_zaid = (int64_t)1000 * nuc->Z + nuc->A;
            for (int j = 0; j < thermal->n_applicable_zaids; j++) {
                int table_zaid = thermal->applicable_zaids[j];
                int table_z = table_zaid / 1000;
                /* ACE readers conventionally expand non-H/Fe entries to the
                 * naturally occurring isotopes of the listed element. */
                if (table_zaid == nuclide_zaid ||
                    (table_z == nuc->Z && table_z != 1 && table_z != 26 &&
                     is_naturally_occurring_isotope(nuc->Z, nuc->A))) {
                    applicable = 1;
                    break;
                }
            }
            if (!applicable) {
                alea_nuc_prepared_material_free(prepared);
                return report_failure(report,
                    ALEA_NUC_PREP_INVALID_THERMAL_ASSOCIATION,
                    ALEA_NUC_CAP_THERMAL_SAB, index, 0,
                    "thermal table does not list the component nuclide as applicable");
            }
            if (!isfinite(nuc->temperature) || nuc->temperature <= 0.0 ||
                fabs(thermal->temperature - nuc->temperature) >
                    requirements->thermal_temperature_tolerance) {
                alea_nuc_prepared_material_free(prepared);
                return report_failure(report,
                    ALEA_NUC_PREP_INVALID_THERMAL_ASSOCIATION,
                    ALEA_NUC_CAP_THERMAL_SAB, index, 0,
                    "thermal and neutron table temperatures do not match");
            }
            prepared->components[index].thermal = thermal;
            if (material->components[index].number_density > 0.0)
                populated_associations++;
        }
        if (populated_associations == 0) {
            alea_nuc_prepared_material_free(prepared);
            return report_failure(report,
                ALEA_NUC_PREP_UNSUPPORTED_CAPABILITY,
                ALEA_NUC_CAP_THERMAL_SAB, -1, 0,
                "thermal scattering requires an associated populated component");
        }
        report->available_capabilities |= ALEA_NUC_CAP_THERMAL_SAB;
    }

    uint32_t available = report->available_capabilities;
    report->missing_capabilities = requirements->required_capabilities & ~available;
    if (report->missing_capabilities) {
        alea_nuc_prepared_material_free(prepared);
        return report_failure(report, ALEA_NUC_PREP_UNSUPPORTED_CAPABILITY,
            report->missing_capabilities, -1, 0,
            "material does not provide a requested collision capability");
    }
    prepared->capabilities = available;
    prepared->requested_capabilities = requirements->required_capabilities;
    *output = prepared;
    return ALEA_OK;
}

static int particle_valid(const alea_nuc_particle_state_t* particle,
                          alea_nuc_particle_t expected) {
    if (!particle || particle->type != expected ||
        !isfinite(particle->energy) || particle->energy <= 0.0 ||
        !isfinite(particle->weight) || particle->weight <= 0.0 ||
        !isfinite(particle->time)) return 0;
    double norm2 = 0.0;
    for (int i = 0; i < 3; i++) {
        if (!isfinite(particle->direction[i])) return 0;
        norm2 += particle->direction[i] * particle->direction[i];
    }
    return fabs(norm2 - 1.0) <= 1e-8;
}

static int fission_mt(int mt) {
    return mt == 18 || mt == 19 || mt == 20 || mt == 21 || mt == 38;
}

static double sampled_reaction_xs(const alea_nuc_nuclide_t* nuc,
                                  const alea_nuc_reaction_t* reaction,
                                  double energy,
                                  const alea_nuc_urr_sample_t* urr) {
    double xs = reaction_xs(nuc, reaction, energy);
    if (!urr || !urr->active) return xs;
    if (fission_mt(reaction->mt)) return xs * urr->factors[2];
    if (reaction->mt == 102) return xs * urr->factors[3];
    return xs;
}

static const alea_nuc_urr_sample_t* component_urr(
    const alea_nuc_evaluation_workspace_t* workspace, int component) {
    if (!workspace || !workspace->components || component < 0 ||
        (size_t)component >= workspace->capacity) return NULL;
    return &workspace->components[component];
}

static int component_thermal_active(const prepared_component_t* component,
                                    double energy) {
    if (!component || !component->thermal || !isfinite(energy) ||
        energy <= 0.0 ||
        energy > component->thermal->inelastic_energy[
            component->thermal->n_inelastic_energies - 1])
        return 0;
    return 1;
}

static alea_error_t evaluate_checked(
    const alea_nuc_prepared_material_t* prepared,
    const alea_nuc_particle_state_t* incident,
    const alea_nuc_evaluation_workspace_t* workspace,
    alea_nuc_evaluation_t* candidate) {
    if (!prepared || !particle_valid(incident, prepared->particle))
        return ALEA_ERR_INVALID_ARG;
    memset(candidate, 0, sizeof(*candidate));
    candidate->prepared = prepared;
    candidate->incident = *incident;
    candidate->workspace = workspace;

    for (int i = 0; i < prepared->n_components; i++) {
        const alea_nuc_mat_component_t* component = prepared->components[i].source;
        const alea_nuc_nuclide_t* nuc = component->nuclide;
        if (component->number_density == 0.0) continue;
        if (incident->energy < nuc->energy[0] ||
            incident->energy > nuc->energy[nuc->n_energies - 1])
            return ALEA_ERR_INVALID_ARG;
        if (prepared->particle == ALEA_NUC_PARTICLE_PHOTON) {
            double coherent = alea_nuc_photon_xs_coherent(nuc, incident->energy);
            double incoherent = alea_nuc_photon_xs_incoherent(nuc, incident->energy);
            double absorption = alea_nuc_photon_xs_photoelectric(
                nuc, incident->energy) +
                alea_nuc_photon_xs_pair(nuc, incident->energy);
            double total = coherent + incoherent + absorption;
            candidate->macro_total += component->number_density * total;
            candidate->macro_elastic += component->number_density *
                                         (coherent + incoherent);
            candidate->macro_absorption += component->number_density * absorption;
            continue;
        }
        const prepared_component_t* prepared_component =
            &prepared->components[i];
        int thermal_active = component_thermal_active(
            prepared_component, incident->energy);
        double microscopic_thermal = thermal_active ?
            alea_nuc_thermal_xs_total(prepared_component->thermal,
                                      incident->energy) : 0.0;
        const alea_nuc_urr_sample_t* urr = component_urr(workspace, i);
        int urr_applies = nuc->urr && incident->energy >= nuc->urr->energy[0] &&
            incident->energy <= nuc->urr->energy[nuc->urr->n_energies - 1];
        if (thermal_active && urr_applies)
            return ALEA_ERR_UNSUPPORTED;
        if (urr_applies && (!urr || !urr->active)) return ALEA_ERR_UNSUPPORTED;
        double microscopic_elastic = thermal_active ? 0.0 :
            alea_nuc_xs_elastic(nuc, incident->energy);
        if (urr && urr->active) microscopic_elastic *= urr->factors[1];
        double microscopic_absorption = 0.0;
        double microscopic_emission = 0.0;
        for (int j = 0; j < prepared_component->n_event_reactions; j++) {
            const alea_nuc_reaction_t* reaction = &nuc->reactions[
                prepared_component->event_reactions[j]];
            double xs = sampled_reaction_xs(nuc, reaction, incident->energy,
                                            urr);
            if (reaction->ty == 0) microscopic_absorption += xs;
            else microscopic_emission += xs;
        }
        double microscopic_total = microscopic_thermal + microscopic_elastic +
            microscopic_absorption + microscopic_emission;
        /* Use the sum of the sampled partial cross sections for transport.
         * Processed probability-table totals can differ slightly from that
         * sum because each column is independently rounded.  The partial sum
         * keeps flight and reaction selection exactly normalized. */
        if (urr && urr->active) {
            double table_total = alea_nuc_xs_total(
                nuc, incident->energy) * urr->factors[0];
            if (fabs(microscopic_total - table_total) >
                1e-3 * fmax(1.0, table_total))
                return ALEA_ERR_UNSUPPORTED;
        }
        candidate->macro_total += component->number_density * microscopic_total;
        candidate->macro_elastic += component->number_density * microscopic_elastic;
        candidate->macro_thermal += component->number_density * microscopic_thermal;
        candidate->macro_absorption += component->number_density * microscopic_absorption;
        candidate->macro_neutron_emission += component->number_density * microscopic_emission;
    }
    if (!isfinite(candidate->macro_total) || candidate->macro_total <= 0.0 ||
        !isfinite(candidate->macro_elastic) ||
        !isfinite(candidate->macro_thermal) ||
        !isfinite(candidate->macro_absorption) ||
        !isfinite(candidate->macro_neutron_emission) ||
        candidate->macro_neutron_emission < -1e-12 ||
        fabs(candidate->macro_total - candidate->macro_elastic -
             candidate->macro_thermal -
             candidate->macro_absorption -
             candidate->macro_neutron_emission) >
            1e-8 * fmax(1.0, candidate->macro_total))
        return ALEA_ERR_INVALID_STATE;
    return ALEA_OK;
}

static int evaluation_matches(const alea_nuc_evaluation_t* evaluation) {
    if (!evaluation || !evaluation->prepared) return 0;
    alea_nuc_evaluation_t current;
    if (evaluate_checked(evaluation->prepared, &evaluation->incident,
                         evaluation->workspace,
                         &current) != ALEA_OK) return 0;
    double scale = fmax(1.0, current.macro_total);
    return fabs(evaluation->macro_total - current.macro_total) <= 1e-12 * scale &&
           fabs(evaluation->macro_elastic - current.macro_elastic) <= 1e-12 * scale &&
           fabs(evaluation->macro_thermal - current.macro_thermal) <= 1e-12 * scale &&
           fabs(evaluation->macro_absorption - current.macro_absorption) <=
               1e-12 * scale &&
           fabs(evaluation->macro_neutron_emission -
                current.macro_neutron_emission) <= 1e-12 * scale;
}

alea_error_t alea_nuc_evaluate(
    const alea_nuc_prepared_material_t* prepared,
    const alea_nuc_particle_state_t* incident,
    alea_nuc_evaluation_t* evaluation) {
    if (!prepared || !incident || !evaluation) return ALEA_ERR_NULL_ARG;
    alea_nuc_evaluation_t candidate;
    alea_error_t err = evaluate_checked(prepared, incident, NULL, &candidate);
    if (err != ALEA_OK) return err;
    *evaluation = candidate;
    return ALEA_OK;
}


alea_error_t alea_nuc_evaluate_urr(
    const alea_nuc_prepared_material_t* prepared,
    const alea_nuc_particle_state_t* incident, alea_nuc_random_fn random,
    void* random_context, alea_nuc_evaluation_workspace_t* workspace,
    alea_nuc_evaluation_t* evaluation) {
    if (!prepared || !incident || !random || !workspace || !evaluation)
        return ALEA_ERR_NULL_ARG;
    if (prepared->particle != ALEA_NUC_PARTICLE_NEUTRON)
        return ALEA_ERR_UNSUPPORTED;
    if (workspace->capacity < (size_t)prepared->n_components ||
        (prepared->n_components > 0 && !workspace->components))
        return ALEA_ERR_INVALID_ARG;
    for (int i = 0; i < prepared->n_components; i++) {
        alea_nuc_urr_sample_t sample;
        memset(&sample, 0, sizeof(sample));
        const alea_nuc_nuclide_t* nuc =
            prepared->components[i].source->nuclide;
        if (nuc->urr && incident->energy >= nuc->urr->energy[0] &&
            incident->energy <= nuc->urr->energy[nuc->urr->n_energies - 1]) {
            double u;
            alea_error_t err = draw_uniform(random, random_context, &u);
            if (err != ALEA_OK) return err;
            if (!alea_nuc_urr_factors(nuc, incident->energy, u,
                                      sample.factors))
                return ALEA_ERR_INVALID_STATE;
            sample.active = true;
        }
        workspace->components[i] = sample;
    }
    alea_nuc_evaluation_t candidate;
    alea_error_t err = evaluate_checked(prepared, incident, workspace,
                                        &candidate);
    if (err != ALEA_OK) return err;
    *evaluation = candidate;
    return ALEA_OK;
}

static alea_error_t draw_uniform(alea_nuc_random_fn random, void* context,
                                 double* value) {
    if (!random || !value) return ALEA_ERR_NULL_ARG;
    double candidate = random(context);
    if (!isfinite(candidate) || candidate < 0.0 || candidate >= 1.0)
        return ALEA_ERR_INVALID_ARG;
    *value = candidate;
    return ALEA_OK;
}

alea_error_t alea_nuc_sample_flight(const alea_nuc_evaluation_t* evaluation,
                                    alea_nuc_random_fn random,
                                    void* random_context,
                                    double* distance) {
    if (!evaluation || !distance || !random) return ALEA_ERR_NULL_ARG;
    if (!evaluation_matches(evaluation)) return ALEA_ERR_INVALID_STATE;
    double xi;
    alea_error_t err = draw_uniform(random, random_context, &xi);
    if (err != ALEA_OK) return err;
    *distance = -log1p(-xi) / evaluation->macro_total;
    return ALEA_OK;
}

static double sample_tabular_mu(const alea_nuc_angular_point_t* point,
                                double xi) {
    double value = point->cosine[0];
    (void)alea_nuc_tabular_pdf_sample(point->cosine, point->pdf, point->cdf,
                                      point->n_cosines, point->interpolation,
                                      0, xi, &value, NULL);
    return value;
}

static double sample_point_mu(const alea_nuc_angular_point_t* point,
                              double xi) {
    if (!point || point->type == ALEA_NUC_ANG_ISOTROPIC)
        return 2.0 * xi - 1.0;
    if (point->type == ALEA_NUC_ANG_EQUIPROBABLE) {
        double scaled = 32.0 * xi;
        int bin = (int)scaled;
        if (bin > 31) bin = 31;
        double f = scaled - bin;
        return point->cosine[bin] +
               f * (point->cosine[bin + 1] - point->cosine[bin]);
    }
    return sample_tabular_mu(point, xi);
}

double alea_nuc_sample_angular_mu_internal(
    const alea_nuc_angular_dist_t* angular, double energy,
    double select_xi, double sample_xi) {
    if (!angular) return 2.0 * sample_xi - 1.0;
    int point = 0;
    if (angular->n_energies > 1) {
        double f;
        int lo = alea_nuc_energy_lookup_trusted(angular->energy,
            angular->n_energies, energy, &f);
        point = lo + (select_xi < f ? 1 : 0);
        if (point >= angular->n_energies) point = angular->n_energies - 1;
    }
    return sample_point_mu(&angular->data[point], sample_xi);
}

static double sample_elastic_mu_cm(const alea_nuc_nuclide_t* nuc,
                                   double energy, double select_xi,
                                   double sample_xi) {
    return alea_nuc_sample_angular_mu_internal(
        nuc->elastic_angular, energy, select_xi, sample_xi);
}

alea_error_t alea_nuc_sample_collision(const alea_nuc_nuclide_t* nuc, int mt,
                                       double energy, const double xi[3],
                                       alea_nuc_interaction_t* result) {
    if (!nuc || !xi || !result) return ALEA_ERR_NULL_ARG;
    if (!isfinite(energy) || energy <= 0.0 || !isfinite(nuc->awr) ||
        nuc->awr <= 0.0) return ALEA_ERR_INVALID_ARG;
    if (nuc->n_energies < 2 || !nuc->energy || energy < nuc->energy[0] ||
        energy > nuc->energy[nuc->n_energies - 1]) return ALEA_ERR_INVALID_ARG;
    for (int i = 0; i < 3; i++)
        if (!isfinite(xi[i]) || xi[i] < 0.0 || xi[i] >= 1.0)
            return ALEA_ERR_INVALID_ARG;

    alea_nuc_interaction_t candidate;
    memset(&candidate, 0, sizeof(candidate));
    candidate.mt = mt;
    candidate.weight_factor = 1.0;
    if (mt == 2) {
        if (!alea_nuc_validate_angular_internal(nuc->elastic_angular))
            return ALEA_ERR_INVALID_STATE;
        candidate.mu = sample_elastic_mu_cm(nuc, energy, xi[0], xi[1]);
        double A = nuc->awr;
        candidate.energy_out = energy *
            (A * A + 2.0 * A * candidate.mu + 1.0) /
            ((A + 1.0) * (A + 1.0));
        candidate.n_secondary = 1;
    } else if (reaction_is_supported_absorption_mt(mt)) {
        int found = 0;
        for (int i = 0; i < nuc->n_reactions; i++) {
            if (nuc->reactions[i].mt == mt && nuc->reactions[i].ty == 0 &&
                reaction_xs(nuc, &nuc->reactions[i], energy) > 0.0) {
                found = 1;
                break;
            }
        }
        if (!found) return ALEA_ERR_INVALID_ARG;
        candidate.energy_out = 0.0;
        candidate.n_secondary = 0;
    } else {
        return ALEA_ERR_UNSUPPORTED;
    }
    *result = candidate;
    return ALEA_OK;
}

void alea_nuc_rotate_direction_internal(
    const double incident[3], double mu, double phi, double output[3]) {
    double sin_theta = sqrt(fmax(0.0, 1.0 - mu * mu));
    double cos_phi = cos(phi);
    double sin_phi = sin(phi);
    double ux = incident[0], uy = incident[1], uz = incident[2];
    if (fabs(uz) < 1.0 - 1e-12) {
        double transverse = sqrt(1.0 - uz * uz);
        output[0] = mu * ux + sin_theta *
            (ux * uz * cos_phi - uy * sin_phi) / transverse;
        output[1] = mu * uy + sin_theta *
            (uy * uz * cos_phi + ux * sin_phi) / transverse;
        output[2] = mu * uz - sin_theta * transverse * cos_phi;
    } else {
        output[0] = sin_theta * cos_phi;
        output[1] = sin_theta * sin_phi;
        output[2] = mu * copysign(1.0, uz);
    }
    double norm = sqrt(output[0] * output[0] + output[1] * output[1] +
                       output[2] * output[2]);
    if (norm > 0.0) {
        output[0] /= norm;
        output[1] /= norm;
        output[2] /= norm;
    }
}

static void evaluated_component_xs(const alea_nuc_evaluation_t* evaluation,
                                   int component_index, double* total,
                                   double* elastic, double* thermal) {
    const prepared_component_t* component =
        &evaluation->prepared->components[component_index];
    const alea_nuc_nuclide_t* nuc = component->source->nuclide;
    const alea_nuc_urr_sample_t* urr = component_urr(
        evaluation->workspace, component_index);
    int thermal_active = component_thermal_active(
        component, evaluation->incident.energy);
    *thermal = thermal_active ? alea_nuc_thermal_xs_total(
        component->thermal, evaluation->incident.energy) : 0.0;
    *elastic = thermal_active ? 0.0 :
        alea_nuc_xs_elastic(nuc, evaluation->incident.energy);
    if (urr && urr->active) *elastic *= urr->factors[1];
    *total = *thermal + *elastic;
    for (int i = 0; i < component->n_event_reactions; i++) {
        const alea_nuc_reaction_t* reaction =
            &nuc->reactions[component->event_reactions[i]];
        *total += sampled_reaction_xs(nuc, reaction,
                                      evaluation->incident.energy, urr);
    }
}

static alea_error_t sample_prepared_elastic(
    const alea_nuc_evaluation_t* evaluation,
    const alea_nuc_nuclide_t* nuc, alea_nuc_random_fn random,
    void* random_context, alea_nuc_collision_result_t* candidate) {
    candidate->outcome = ALEA_NUC_OUTCOME_SCATTERED;
    candidate->mt = 2;

    if ((evaluation->prepared->requested_capabilities &
         ALEA_NUC_CAP_FREE_GAS) &&
        evaluation->incident.energy <= 400.0 * nuc->temperature) {
        alea_nuc_free_gas_result_t free_gas;
        alea_error_t err = alea_nuc_sample_free_gas_elastic(
            nuc, &evaluation->incident, nuc->temperature,
            random, random_context, &free_gas);
        if (err != ALEA_OK) return err;
        candidate->outgoing = free_gas.outgoing;
        candidate->mu_cm = free_gas.mu_cm;
        candidate->mu_lab = free_gas.mu_lab;
    } else {
        double angular_select, angular_sample, azimuth;
        alea_error_t err;
        if ((err = draw_uniform(random, random_context,
                                &angular_select)) != ALEA_OK ||
            (err = draw_uniform(random, random_context,
                                &angular_sample)) != ALEA_OK ||
            (err = draw_uniform(random, random_context, &azimuth)) != ALEA_OK)
            return err;
        double xi[3] = {angular_select, angular_sample, 0.0};
        alea_nuc_interaction_t interaction;
        err = alea_nuc_sample_collision(nuc, 2, evaluation->incident.energy,
                                        xi, &interaction);
        if (err != ALEA_OK) return err;
        candidate->mu_cm = interaction.mu;
        double ratio = interaction.energy_out / evaluation->incident.energy;
        candidate->mu_lab = ratio > 0.0
            ? (1.0 + nuc->awr * candidate->mu_cm) /
              ((nuc->awr + 1.0) * sqrt(ratio))
            : 1.0;
        candidate->mu_lab = fmin(1.0, fmax(-1.0, candidate->mu_lab));
        candidate->outgoing = evaluation->incident;
        candidate->outgoing.energy = interaction.energy_out;
        alea_nuc_rotate_direction_internal(
            evaluation->incident.direction, candidate->mu_lab,
            2.0 * M_PI * azimuth, candidate->outgoing.direction);
    }
    candidate->deposition_available = true;
    candidate->local_energy_deposition =
        evaluation->incident.energy - candidate->outgoing.energy;
    return ALEA_OK;
}

static alea_error_t collide_photon(const alea_nuc_evaluation_t* evaluation,
                                   alea_nuc_random_fn random,
                                   void* random_context,
                                   alea_nuc_secondary_buffer_t* secondaries,
                                   alea_nuc_collision_result_t* result) {
    if (secondaries) {
        for (int i = 0; i < evaluation->prepared->n_components; i++) {
            const alea_nuc_nuclide_t* element =
                evaluation->prepared->components[i].source->nuclide;
            size_t required = alea_nuc_photon_secondary_capacity(
                element, evaluation->incident.energy);
            if (secondaries->capacity - secondaries->count < required)
                return ALEA_ERR_OUT_OF_MEMORY;
        }
    }
    double u;
    alea_error_t err = draw_uniform(random, random_context, &u);
    if (err != ALEA_OK) return err;
    double threshold = u * evaluation->macro_total;
    double sum = 0.0;
    int selected = -1;
    for (int i = 0; i < evaluation->prepared->n_components; i++) {
        const alea_nuc_mat_component_t* component =
            evaluation->prepared->components[i].source;
        double contribution = component->number_density *
            alea_nuc_xs_total(component->nuclide,
                              evaluation->incident.energy);
        if (contribution > 0.0) selected = i;
        sum += contribution;
        if (threshold < sum) { selected = i; break; }
    }
    if (selected < 0) return ALEA_ERR_INVALID_STATE;
    alea_nuc_collision_result_t candidate;
    if (secondaries)
        err = alea_nuc_sample_photon_collision_with_secondaries(
            evaluation->prepared->components[selected].source->nuclide,
            &evaluation->incident, random, random_context, secondaries,
            &candidate);
    else
        err = alea_nuc_sample_photon_collision(
            evaluation->prepared->components[selected].source->nuclide,
            &evaluation->incident, random, random_context, &candidate);
    if (err != ALEA_OK) return err;
    candidate.component_index = selected;
    *result = candidate;
    return ALEA_OK;
}

alea_error_t alea_nuc_collide(const alea_nuc_evaluation_t* evaluation,
                              alea_nuc_random_fn random,
                              void* random_context,
                              alea_nuc_collision_result_t* result) {
    if (!evaluation || !random || !result) return ALEA_ERR_NULL_ARG;
    const alea_nuc_prepared_material_t* prepared = evaluation->prepared;
    if (!prepared || !evaluation_matches(evaluation))
        return ALEA_ERR_INVALID_STATE;
    if (prepared->particle == ALEA_NUC_PARTICLE_PHOTON)
        return collide_photon(evaluation, random, random_context, NULL, result);

    double target_xi, reaction_xi;
    alea_error_t err = draw_uniform(random, random_context, &target_xi);
    if (err != ALEA_OK) return err;
    err = draw_uniform(random, random_context, &reaction_xi);
    if (err != ALEA_OK) return err;

    double target_threshold = target_xi * evaluation->macro_total;
    double cumulative = 0.0;
    int selected = -1;
    double microscopic_total = 0.0;
    for (int i = 0; i < prepared->n_components; i++) {
        const alea_nuc_mat_component_t* component = prepared->components[i].source;
        double ignored_elastic, ignored_thermal;
        evaluated_component_xs(evaluation, i, &microscopic_total,
                               &ignored_elastic, &ignored_thermal);
        double contribution = component->number_density * microscopic_total;
        if (contribution > 0.0) selected = i;
        cumulative += contribution;
        if (target_threshold < cumulative) {
            selected = i;
            break;
        }
    }
    if (selected < 0) return ALEA_ERR_INVALID_STATE;

    const prepared_component_t* component = &prepared->components[selected];
    const alea_nuc_nuclide_t* nuc = component->source->nuclide;
    double elastic, thermal;
    evaluated_component_xs(evaluation, selected, &microscopic_total, &elastic,
                           &thermal);
    double reaction_threshold = reaction_xi * microscopic_total;

    alea_nuc_collision_result_t candidate;
    memset(&candidate, 0, sizeof(candidate));
    candidate.component_index = selected;
    candidate.mt = 2;

    if (reaction_threshold < thermal) {
        err = alea_nuc_sample_thermal_collision(
            component->thermal, &evaluation->incident, random,
            random_context, &candidate);
        if (err != ALEA_OK) return err;
        candidate.component_index = selected;
    } else if (reaction_threshold < thermal + elastic) {
        err = sample_prepared_elastic(
            evaluation, nuc, random, random_context, &candidate);
        if (err != ALEA_OK) return err;
    } else {
        double nonelastic = microscopic_total - thermal - elastic;
        if (nonelastic <= 0.0 || component->n_event_reactions <= 0)
            return ALEA_ERR_INVALID_STATE;
        double residual = reaction_threshold - thermal - elastic;
        double sum = 0.0;
        int reaction_index = component->event_reactions[
            component->n_event_reactions - 1];
        for (int i = 0; i < component->n_event_reactions; i++) {
            int index = component->event_reactions[i];
            sum += sampled_reaction_xs(
                nuc, &nuc->reactions[index], evaluation->incident.energy,
                component_urr(evaluation->workspace, selected));
            if (sum >= residual) {
                reaction_index = index;
                break;
            }
        }
        if (nuc->reactions[reaction_index].ty != 0)
            return ALEA_ERR_UNSUPPORTED;
        candidate.outcome = ALEA_NUC_OUTCOME_ABSORBED;
        candidate.mt = nuc->reactions[reaction_index].mt;
        candidate.deposition_available = false;
        candidate.local_energy_deposition = NAN;
    }

    *result = candidate;
    return ALEA_OK;
}

static alea_error_t maximum_emission_count(
    const alea_nuc_evaluation_t* evaluation, size_t* maximum) {
    *maximum = 0;
    const alea_nuc_prepared_material_t* prepared = evaluation->prepared;
    for (int i = 0; i < prepared->n_components; i++) {
        const prepared_component_t* component = &prepared->components[i];
        const alea_nuc_nuclide_t* nuc = component->source->nuclide;
        if (component->source->number_density <= 0.0) continue;
        for (int j = 0; j < component->n_event_reactions; j++) {
            const alea_nuc_reaction_t* reaction =
                &nuc->reactions[component->event_reactions[j]];
            if (reaction_xs(nuc, reaction, evaluation->incident.energy) <= 0.0)
                continue;
            size_t bound = 0;
            if (reaction->ty != 0) {
                double yield = alea_nuc_reaction_yield(
                    nuc, reaction->mt, evaluation->incident.energy);
                if (!isfinite(yield) || yield < 0.0 ||
                    yield > (double)SIZE_MAX)
                    return ALEA_ERR_INVALID_STATE;
                bound = (size_t)ceil(yield);
            }
            if (reaction->ty != 0 && fission_mt(reaction->mt) &&
                (prepared->requested_capabilities &
                 ALEA_NUC_CAP_DELAYED_NEUTRON)) {
                double delayed = alea_nuc_delayed_nu_bar(
                    nuc, evaluation->incident.energy);
                if (!isfinite(delayed) || delayed < 0.0 ||
                    delayed > (double)SIZE_MAX - (double)bound)
                    return ALEA_ERR_INVALID_STATE;
                bound += (size_t)ceil(delayed);
            }
            if (prepared->requested_capabilities &
                ALEA_NUC_CAP_PHOTON_PRODUCTION) {
                for (int p = 0; p < nuc->n_photon_productions; p++) {
                    const alea_nuc_photon_production_t* production =
                        &nuc->photon_productions[p];
                    if (!photon_production_applies(
                            nuc, production, reaction->mt))
                        continue;
                    double photons = photon_yield_for_event(
                        nuc, production, reaction->mt,
                        evaluation->incident.energy);
                    if (!isfinite(photons) || photons < 0.0 ||
                        photons > (double)(SIZE_MAX - bound))
                        return ALEA_ERR_INVALID_STATE;
                    bound += (size_t)ceil(photons);
                }
            }
            if (bound > *maximum) *maximum = bound;
        }
    }
    return ALEA_OK;
}

static void cm_neutron_to_lab(double incident_energy, double awr,
                              double energy_cm, double mu_cm,
                              double* energy_lab, double* mu_lab) {
    double translation = sqrt(incident_energy) / (awr + 1.0);
    double outgoing = sqrt(fmax(0.0, energy_cm));
    *energy_lab = energy_cm + translation * translation +
                  2.0 * translation * outgoing * mu_cm;
    if (*energy_lab > 0.0)
        *mu_lab = (outgoing * mu_cm + translation) / sqrt(*energy_lab);
    else
        *mu_lab = 1.0;
    *mu_lab = fmin(1.0, fmax(-1.0, *mu_lab));
}

static alea_error_t sample_delayed_group(
    const alea_nuc_nuclide_t* nuc, double incident_energy,
    alea_nuc_random_fn random, void* random_context,
    const alea_nuc_delayed_group_t** selected) {
    int n = nuc->fission->n_delayed_groups;
    double total = 0.0;
    for (int g = 0; g < n; g++) {
        const alea_nuc_delayed_group_t* group =
            &nuc->fission->delayed_groups[g];
        double probability;
        alea_error_t err = alea_nuc_interp_eval(
            group->energy, group->probability, group->n_energies,
            group->nbt, group->interp, group->n_regions, incident_energy,
            &probability);
        if (err != ALEA_OK || probability < 0.0) return ALEA_ERR_INVALID_STATE;
        total += probability;
    }
    if (!(total > 0.0) || !isfinite(total)) return ALEA_ERR_INVALID_STATE;
    double u;
    alea_error_t err = draw_uniform(random, random_context, &u);
    if (err != ALEA_OK) return err;
    double threshold = u * total, cumulative = 0.0;
    *selected = &nuc->fission->delayed_groups[n - 1];
    for (int g = 0; g < n; g++) {
        const alea_nuc_delayed_group_t* group =
            &nuc->fission->delayed_groups[g];
        double probability;
        (void)alea_nuc_interp_eval(
            group->energy, group->probability, group->n_energies,
            group->nbt, group->interp, group->n_regions, incident_energy,
            &probability);
        cumulative += probability;
        if (threshold < cumulative) { *selected = group; break; }
    }
    return ALEA_OK;
}

static alea_error_t sample_delayed_particle(
    const alea_nuc_nuclide_t* nuc,
    const alea_nuc_particle_state_t* incident,
    alea_nuc_random_fn random, void* random_context,
    alea_nuc_particle_state_t* particle) {
    const alea_nuc_delayed_group_t* group;
    alea_error_t err = sample_delayed_group(
        nuc, incident->energy, random, random_context, &group);
    if (err != ALEA_OK) return err;
    double energy, mu = 0.0;
    bool correlated = false;
    err = alea_nuc_sample_energy_angle_distribution(
        group->spectrum, incident->energy, random, random_context,
        &energy, &mu, &correlated);
    if (err != ALEA_OK) return err;
    if (!correlated) {
        double u;
        if ((err = draw_uniform(random, random_context, &u)) != ALEA_OK)
            return err;
        mu = 2.0 * u - 1.0;
    }
    double azimuth, lifetime_draw;
    if ((err = draw_uniform(random, random_context, &azimuth)) != ALEA_OK ||
        (err = draw_uniform(random, random_context, &lifetime_draw)) != ALEA_OK)
        return err;
    *particle = *incident;
    particle->energy = energy;
    particle->time += -log1p(-lifetime_draw) / group->decay_rate;
    alea_nuc_rotate_direction_internal(
        incident->direction, mu, 2.0 * M_PI * azimuth, particle->direction);
    return ALEA_OK;
}

alea_error_t alea_nuc_collide_with_secondaries(
    const alea_nuc_evaluation_t* evaluation, alea_nuc_random_fn random,
    void* random_context, alea_nuc_secondary_buffer_t* secondaries,
    alea_nuc_collision_result_t* result) {
    if (!evaluation || !random || !secondaries || !result)
        return ALEA_ERR_NULL_ARG;
    if (!evaluation->prepared || !evaluation_matches(evaluation) ||
        secondaries->count > secondaries->capacity ||
        (secondaries->capacity > 0 && !secondaries->particles))
        return ALEA_ERR_INVALID_STATE;
    if (evaluation->prepared->particle == ALEA_NUC_PARTICLE_PHOTON)
        return collide_photon(evaluation, random, random_context, secondaries,
                              result);

    size_t maximum;
    alea_error_t err = maximum_emission_count(evaluation, &maximum);
    if (err != ALEA_OK) return err;
    if (maximum > secondaries->capacity - secondaries->count)
        return ALEA_ERR_OUT_OF_MEMORY;

    const alea_nuc_prepared_material_t* prepared = evaluation->prepared;
    double target_xi, reaction_xi;
    if ((err = draw_uniform(random, random_context, &target_xi)) != ALEA_OK ||
        (err = draw_uniform(random, random_context, &reaction_xi)) != ALEA_OK)
        return err;

    double target_threshold = target_xi * evaluation->macro_total;
    double cumulative = 0.0;
    int selected = -1;
    for (int i = 0; i < prepared->n_components; i++) {
        const alea_nuc_mat_component_t* source = prepared->components[i].source;
        double component_total, ignored_elastic, ignored_thermal;
        evaluated_component_xs(evaluation, i, &component_total,
                               &ignored_elastic, &ignored_thermal);
        double contribution = source->number_density * component_total;
        if (contribution > 0.0) selected = i;
        cumulative += contribution;
        if (target_threshold < cumulative) { selected = i; break; }
    }
    if (selected < 0) return ALEA_ERR_INVALID_STATE;

    const prepared_component_t* component = &prepared->components[selected];
    const alea_nuc_nuclide_t* nuc = component->source->nuclide;
    double microscopic_total, elastic, thermal;
    evaluated_component_xs(evaluation, selected, &microscopic_total, &elastic,
                           &thermal);
    double threshold = reaction_xi * microscopic_total;
    alea_nuc_collision_result_t candidate;
    memset(&candidate, 0, sizeof(candidate));
    candidate.component_index = selected;

    if (threshold < thermal) {
        err = alea_nuc_sample_thermal_collision(
            component->thermal, &evaluation->incident, random,
            random_context, &candidate);
        if (err != ALEA_OK) return err;
        candidate.component_index = selected;
        *result = candidate;
        return ALEA_OK;
    }

    if (threshold < thermal + elastic) {
        err = sample_prepared_elastic(
            evaluation, nuc, random, random_context, &candidate);
        if (err != ALEA_OK) return err;
        *result = candidate;
        return ALEA_OK;
    }

    double residual = threshold - thermal - elastic;
    double sum = 0.0;
    int reaction_index = -1;
    for (int i = 0; i < component->n_event_reactions; i++) {
        int index = component->event_reactions[i];
        double xs = sampled_reaction_xs(
            nuc, &nuc->reactions[index], evaluation->incident.energy,
            component_urr(evaluation->workspace, selected));
        if (xs > 0.0) reaction_index = index;
        sum += xs;
        if (residual < sum) { reaction_index = index; break; }
    }
    if (reaction_index < 0) return ALEA_ERR_INVALID_STATE;
    const alea_nuc_reaction_t* reaction = &nuc->reactions[reaction_index];
    candidate.mt = reaction->mt;
    candidate.deposition_available = false;
    candidate.local_energy_deposition = NAN;
    size_t count = 0;
    if (reaction->ty != 0) {
        double yield = alea_nuc_reaction_yield(
            nuc, reaction->mt, evaluation->incident.energy);
        count = (size_t)floor(yield);
        double fraction = yield - floor(yield);
        if (fraction > 0.0) {
            double u;
            if ((err = draw_uniform(random, random_context, &u)) != ALEA_OK)
                return err;
            if (u < fraction) count++;
        }
    }
    size_t delayed_count = 0;
    if (reaction->ty != 0 && fission_mt(reaction->mt) &&
        (prepared->requested_capabilities & ALEA_NUC_CAP_DELAYED_NEUTRON)) {
        double delayed_yield = alea_nuc_delayed_nu_bar(
            nuc, evaluation->incident.energy);
        delayed_count = (size_t)floor(delayed_yield);
        double delayed_fraction = delayed_yield - floor(delayed_yield);
        if (delayed_fraction > 0.0) {
            double u;
            if ((err = draw_uniform(random, random_context, &u)) != ALEA_OK)
                return err;
            if (u < delayed_fraction) delayed_count++;
        }
    }

    size_t start = secondaries->count;
    for (size_t i = 0; i < count; i++) {
        double energy, mu = 0.0, azimuth;
        bool correlated = false;
        err = alea_nuc_sample_energy_angle_distribution(
            reaction->energy, evaluation->incident.energy, random,
            random_context, &energy, &mu, &correlated);
        if (err != ALEA_OK) return err;
        if (!correlated) {
            double angular_select, angular_sample;
            if ((err = draw_uniform(random, random_context,
                                    &angular_select)) != ALEA_OK ||
                (err = draw_uniform(random, random_context,
                                    &angular_sample)) != ALEA_OK)
                return err;
            mu = alea_nuc_sample_angular_mu_internal(
                reaction->angular, evaluation->incident.energy,
                angular_select, angular_sample);
        }
        if ((err = draw_uniform(random, random_context, &azimuth)) != ALEA_OK)
            return err;
        double energy_lab = energy;
        double mu_lab = mu;
        if (reaction->center_of_mass)
            cm_neutron_to_lab(evaluation->incident.energy, nuc->awr,
                              energy, mu, &energy_lab, &mu_lab);
        alea_nuc_particle_state_t particle = evaluation->incident;
        particle.energy = energy_lab;
        alea_nuc_rotate_direction_internal(
            evaluation->incident.direction, mu_lab,
            2.0 * M_PI * azimuth, particle.direction);
        secondaries->particles[start + i] = particle;
    }
    for (size_t i = 0; i < delayed_count; i++) {
        err = sample_delayed_particle(
            nuc, &evaluation->incident, random, random_context,
            &secondaries->particles[start + count + i]);
        if (err != ALEA_OK) return err;
    }
    size_t photon_count = 0;
    if (prepared->requested_capabilities &
        ALEA_NUC_CAP_PHOTON_PRODUCTION) {
        for (int p = 0; p < nuc->n_photon_productions; p++) {
            const alea_nuc_photon_production_t* production =
                &nuc->photon_productions[p];
            if (!photon_production_applies(
                    nuc, production, reaction->mt)) continue;
            double yield = photon_yield_for_event(
                nuc, production, reaction->mt,
                evaluation->incident.energy);
            size_t photons = (size_t)floor(yield);
            double fraction = yield - floor(yield);
            if (fraction > 0.0) {
                double u;
                if ((err = draw_uniform(random, random_context, &u)) != ALEA_OK)
                    return err;
                if (u < fraction) photons++;
            }
            for (size_t i = 0; i < photons; i++) {
                size_t destination = start + count + delayed_count +
                                     photon_count;
                err = alea_nuc_sample_photon_production(
                    nuc, production, &evaluation->incident, random,
                    random_context, &secondaries->particles[destination]);
                if (err != ALEA_OK) return err;
                photon_count++;
            }
        }
    }
    size_t emitted = count + delayed_count + photon_count;
    secondaries->count = start + emitted;
    candidate.outcome = emitted > 0 ? ALEA_NUC_OUTCOME_REPLACED
                                  : ALEA_NUC_OUTCOME_ABSORBED;
    candidate.n_emitted = emitted;
    *result = candidate;
    return ALEA_OK;
}
