// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/**
 * @file collision.c
 * @brief Prepared stationary-target neutron collision physics
 */

#include "nuclear_internal.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

typedef struct {
    const alea_nuc_mat_component_t* source;
    int* absorption_reactions;
    int n_absorption_reactions;
} prepared_component_t;

struct alea_nuc_prepared_material {
    const alea_nuc_material_t* source;
    prepared_component_t* components;
    int n_components;
    uint32_t capabilities;
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
    int aggregate = aggregate_mt_for(reaction->mt);
    return aggregate != 0 && nuclide_has_active_mt(nuc, aggregate);
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
        if (!isfinite(point->cosine[i]) || !isfinite(point->pdf[i]) ||
            !isfinite(point->cdf[i]) || point->cosine[i] < -1.0 ||
            point->cosine[i] > 1.0 || point->pdf[i] < 0.0 ||
            point->cdf[i] < 0.0 || point->cdf[i] > 1.0 ||
            (i > 0 && (point->cosine[i] <= point->cosine[i - 1] ||
                       point->cdf[i] < point->cdf[i - 1]))) return 0;
    }
    return fabs(point->cdf[0]) <= 1e-8 &&
           fabs(point->cdf[n - 1] - 1.0) <= 1e-6;
}

static int validate_angular(const alea_nuc_angular_dist_t* angular) {
    if (!angular) return 1; /* absent LAND means isotropic */
    if (angular->n_energies <= 0 || !angular->energy || !angular->data)
        return 0;
    for (int i = 0; i < angular->n_energies; i++) {
        if (!isfinite(angular->energy[i]) ||
            (i > 0 && angular->energy[i] <= angular->energy[i - 1]) ||
            !validate_angular_point(&angular->data[i])) return 0;
    }
    return 1;
}

static alea_error_t inspect_nuclide(const alea_nuc_nuclide_t* nuc,
                                    int component,
                                    alea_nuc_capability_report_t* report) {
    if (!nuc || nuc->particle != ALEA_NUC_PARTICLE_NEUTRON)
        return report_failure(report, ALEA_NUC_PREP_UNSUPPORTED_PARTICLE,
            ALEA_NUC_CAP_RESTRICTED_NEUTRON, component, 0,
            "the restricted collision model accepts neutron tables only");
    if (nuc->urr)
        return report_failure(report, ALEA_NUC_PREP_UNSUPPORTED_URR,
            ALEA_NUC_CAP_URR, component, 0,
            "URR data requires coordinated probability-table evaluation");
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
    if (!validate_angular(nuc->elastic_angular))
        return report_failure(report, ALEA_NUC_PREP_INVALID_ANGULAR,
            ALEA_NUC_CAP_STATIONARY_ELASTIC, component, 2,
            "elastic angular distribution is not sampleable");

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
        alea_nuc_reaction_class_t cls = alea_nuc_reaction_classify(reaction->mt);
        if (!reaction_is_supported_absorption_mt(reaction->mt) ||
            reaction->ty != 0) {
            char detail[192];
            if (reaction->ty != 0)
                snprintf(detail, sizeof(detail),
                         "MT %d produces neutrons and is outside the restricted model",
                         reaction->mt);
            else
                snprintf(detail, sizeof(detail),
                         "MT %d is not a supported absorption event channel",
                         reaction->mt);
            return report_failure(report, ALEA_NUC_PREP_UNSUPPORTED_REACTION,
                cls == ALEA_NUC_RXN_MULTIPLY ? ALEA_NUC_CAP_FISSION
                                             : ALEA_NUC_CAP_RESTRICTED_NEUTRON,
                component, reaction->mt, detail);
        }
    }

    for (int i = 0; i < nuc->n_energies; i++) {
        double total = nuc->sigma_total[i];
        double elastic = nuc->sigma_elastic[i];
        double absorption = nuc->sigma_abs[i];
        if (!isfinite(nuc->energy[i]) || !isfinite(total) ||
            !isfinite(elastic) || !isfinite(absorption) || total < 0.0 ||
            elastic < 0.0 || absorption < 0.0 ||
            (i > 0 && nuc->energy[i] <= nuc->energy[i - 1]))
            return report_failure(report, ALEA_NUC_PREP_INVALID_CROSS_SECTIONS,
                ALEA_NUC_CAP_RESTRICTED_NEUTRON, component, 0,
                "cross-section grids must be finite, nonnegative and ascending");

        double scale = fmax(1.0, total);
        if (fabs(total - elastic - absorption) > 1e-6 * scale)
            return report_failure(report, ALEA_NUC_PREP_INVALID_CROSS_SECTIONS,
                ALEA_NUC_CAP_RESTRICTED_NEUTRON, component, 0,
                "total cross section is inconsistent with elastic plus absorption");

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

    if (report) report->available_capabilities |= ALEA_NUC_CAP_RESTRICTED_NEUTRON;
    return ALEA_OK;
}

alea_error_t alea_nuc_capabilities(const alea_nuc_nuclide_t* nuc,
                                   alea_nuc_capability_report_t* report) {
    if (!nuc || !report) return ALEA_ERR_NULL_ARG;
    report_clear(report);
    return inspect_nuclide(nuc, -1, report);
}

static void prepared_free_components(alea_nuc_prepared_material_t* prepared) {
    if (!prepared) return;
    for (int i = 0; i < prepared->n_components; i++)
        free(prepared->components[i].absorption_reactions);
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

    uint32_t supported = ALEA_NUC_CAP_RESTRICTED_NEUTRON;
    uint32_t unknown = requirements->required_capabilities & ~supported;
    if (unknown)
        return report_failure(report, ALEA_NUC_PREP_UNSUPPORTED_CAPABILITY,
            unknown, -1, 0, "requested collision capability is not implemented");
    if (material->n_components <= 0 || !material->components)
        return report_failure(report, ALEA_NUC_PREP_EMPTY_MATERIAL,
            requirements->required_capabilities, -1, 0,
            "material has no nuclide components");

    alea_nuc_prepared_material_t* prepared = calloc(1, sizeof(*prepared));
    if (!prepared) return ALEA_ERR_OUT_OF_MEMORY;
    prepared->components = calloc((size_t)material->n_components,
                                  sizeof(*prepared->components));
    if (!prepared->components) {
        free(prepared);
        return ALEA_ERR_OUT_OF_MEMORY;
    }
    prepared->source = material;
    prepared->n_components = material->n_components;
    prepared->capabilities = supported;

    int populated_components = 0;
    for (int i = 0; i < material->n_components; i++) {
        const alea_nuc_mat_component_t* source = &material->components[i];
        if (!isfinite(source->number_density) || source->number_density < 0.0) {
            prepared_free_components(prepared);
            free(prepared);
            return report_failure(report, ALEA_NUC_PREP_INVALID_CROSS_SECTIONS,
                requirements->required_capabilities, i, 0,
                "component number density is invalid");
        }
        if (source->number_density > 0.0) populated_components++;
        alea_error_t err = inspect_nuclide(source->nuclide, i, report);
        if (err != ALEA_OK) {
            prepared_free_components(prepared);
            free(prepared);
            return err;
        }
        prepared_component_t* component = &prepared->components[i];
        component->source = source;
        int maximum = source->nuclide->n_reactions;
        if (maximum > 0) {
            component->absorption_reactions = malloc((size_t)maximum * sizeof(int));
            if (!component->absorption_reactions) {
                alea_nuc_prepared_material_free(prepared);
                return ALEA_ERR_OUT_OF_MEMORY;
            }
        }
        for (int r = 0; r < maximum; r++) {
            const alea_nuc_reaction_t* reaction = &source->nuclide->reactions[r];
            if (reaction_is_supported_absorption_mt(reaction->mt) &&
                reaction->ty == 0 && reaction_has_positive_xs(reaction) &&
                !reaction_is_redundant(source->nuclide, reaction))
                component->absorption_reactions[component->n_absorption_reactions++] = r;
        }
    }

    if (populated_components == 0) {
        alea_nuc_prepared_material_free(prepared);
        return report_failure(report, ALEA_NUC_PREP_EMPTY_MATERIAL,
            requirements->required_capabilities, -1, 0,
            "material has no positive number density");
    }

    report->available_capabilities = supported;
    report->missing_capabilities = requirements->required_capabilities & ~supported;
    *output = prepared;
    return ALEA_OK;
}

static int particle_valid(const alea_nuc_particle_state_t* particle) {
    if (!particle || particle->type != ALEA_NUC_PARTICLE_NEUTRON ||
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

static alea_error_t evaluate_checked(
    const alea_nuc_prepared_material_t* prepared,
    const alea_nuc_particle_state_t* incident,
    alea_nuc_evaluation_t* candidate) {
    if (!prepared || !particle_valid(incident)) return ALEA_ERR_INVALID_ARG;
    memset(candidate, 0, sizeof(*candidate));
    candidate->prepared = prepared;
    candidate->incident = *incident;

    for (int i = 0; i < prepared->n_components; i++) {
        const alea_nuc_mat_component_t* component = prepared->components[i].source;
        const alea_nuc_nuclide_t* nuc = component->nuclide;
        if (component->number_density == 0.0) continue;
        if (incident->energy < nuc->energy[0] ||
            incident->energy > nuc->energy[nuc->n_energies - 1])
            return ALEA_ERR_INVALID_ARG;
        candidate->macro_total += component->number_density *
            alea_nuc_xs_total(nuc, incident->energy);
        candidate->macro_elastic += component->number_density *
            alea_nuc_xs_elastic(nuc, incident->energy);
        candidate->macro_absorption += component->number_density *
            alea_nuc_xs_absorption(nuc, incident->energy);
    }
    if (!isfinite(candidate->macro_total) || candidate->macro_total <= 0.0 ||
        !isfinite(candidate->macro_elastic) ||
        !isfinite(candidate->macro_absorption) ||
        fabs(candidate->macro_total - candidate->macro_elastic -
             candidate->macro_absorption) >
            1e-8 * fmax(1.0, candidate->macro_total))
        return ALEA_ERR_INVALID_STATE;
    return ALEA_OK;
}

static int evaluation_matches(const alea_nuc_evaluation_t* evaluation) {
    if (!evaluation || !evaluation->prepared) return 0;
    alea_nuc_evaluation_t current;
    if (evaluate_checked(evaluation->prepared, &evaluation->incident,
                         &current) != ALEA_OK) return 0;
    double scale = fmax(1.0, current.macro_total);
    return fabs(evaluation->macro_total - current.macro_total) <= 1e-12 * scale &&
           fabs(evaluation->macro_elastic - current.macro_elastic) <= 1e-12 * scale &&
           fabs(evaluation->macro_absorption - current.macro_absorption) <=
               1e-12 * scale;
}

alea_error_t alea_nuc_evaluate(
    const alea_nuc_prepared_material_t* prepared,
    const alea_nuc_particle_state_t* incident,
    alea_nuc_evaluation_t* evaluation) {
    if (!prepared || !incident || !evaluation) return ALEA_ERR_NULL_ARG;
    alea_nuc_evaluation_t candidate;
    alea_error_t err = evaluate_checked(prepared, incident, &candidate);
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
    int hi = 1;
    while (hi < point->n_cosines && point->cdf[hi] < xi) hi++;
    if (hi >= point->n_cosines) return point->cosine[point->n_cosines - 1];
    int lo = hi - 1;
    double probability = xi - point->cdf[lo];
    double width = point->cosine[hi] - point->cosine[lo];
    if (probability <= 0.0 || width <= 0.0) return point->cosine[lo];
    if (point->interpolation == 1 ||
        fabs(point->pdf[hi] - point->pdf[lo]) <= 1e-14) {
        if (point->pdf[lo] <= 0.0) return point->cosine[lo];
        return fmin(point->cosine[hi],
                    point->cosine[lo] + probability / point->pdf[lo]);
    }
    double slope = (point->pdf[hi] - point->pdf[lo]) / width;
    double discriminant = point->pdf[lo] * point->pdf[lo] +
                          2.0 * slope * probability;
    if (discriminant < 0.0) discriminant = 0.0;
    double denominator = point->pdf[lo] + sqrt(discriminant);
    double offset = denominator > 0.0 ? 2.0 * probability / denominator : 0.0;
    return fmin(point->cosine[hi],
                fmax(point->cosine[lo], point->cosine[lo] + offset));
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

static double sample_elastic_mu_cm(const alea_nuc_nuclide_t* nuc,
                                   double energy, double select_xi,
                                   double sample_xi) {
    const alea_nuc_angular_dist_t* angular = nuc->elastic_angular;
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
        if (!validate_angular(nuc->elastic_angular)) return ALEA_ERR_INVALID_STATE;
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

static void rotate_direction(const double incident[3], double mu,
                             double phi, double output[3]) {
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
        output[2] = copysign(mu, uz);
    }
    double norm = sqrt(output[0] * output[0] + output[1] * output[1] +
                       output[2] * output[2]);
    if (norm > 0.0) {
        output[0] /= norm;
        output[1] /= norm;
        output[2] /= norm;
    }
}

alea_error_t alea_nuc_collide(const alea_nuc_evaluation_t* evaluation,
                              alea_nuc_random_fn random,
                              void* random_context,
                              alea_nuc_collision_result_t* result) {
    if (!evaluation || !random || !result) return ALEA_ERR_NULL_ARG;
    const alea_nuc_prepared_material_t* prepared = evaluation->prepared;
    if (!prepared || !evaluation_matches(evaluation))
        return ALEA_ERR_INVALID_STATE;

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
        microscopic_total = alea_nuc_xs_total(component->nuclide,
                                               evaluation->incident.energy);
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
    microscopic_total = alea_nuc_xs_total(nuc, evaluation->incident.energy);
    double elastic = alea_nuc_xs_elastic(nuc, evaluation->incident.energy);
    double reaction_threshold = reaction_xi * microscopic_total;

    alea_nuc_collision_result_t candidate;
    memset(&candidate, 0, sizeof(candidate));
    candidate.component_index = selected;
    candidate.mt = 2;

    if (reaction_threshold < elastic) {
        double angular_select, angular_sample, azimuth;
        if ((err = draw_uniform(random, random_context, &angular_select)) != ALEA_OK ||
            (err = draw_uniform(random, random_context, &angular_sample)) != ALEA_OK ||
            (err = draw_uniform(random, random_context, &azimuth)) != ALEA_OK)
            return err;
        double xi[3] = {angular_select, angular_sample, 0.0};
        alea_nuc_interaction_t interaction;
        err = alea_nuc_sample_collision(nuc, 2, evaluation->incident.energy,
                                        xi, &interaction);
        if (err != ALEA_OK) return err;
        candidate.outcome = ALEA_NUC_OUTCOME_SCATTERED;
        candidate.mu_cm = interaction.mu;
        double ratio = interaction.energy_out / evaluation->incident.energy;
        candidate.mu_lab = ratio > 0.0
            ? (1.0 + nuc->awr * candidate.mu_cm) /
              ((nuc->awr + 1.0) * sqrt(ratio))
            : 1.0;
        if (candidate.mu_lab > 1.0) candidate.mu_lab = 1.0;
        if (candidate.mu_lab < -1.0) candidate.mu_lab = -1.0;
        candidate.outgoing = evaluation->incident;
        candidate.outgoing.energy = interaction.energy_out;
        rotate_direction(evaluation->incident.direction, candidate.mu_lab,
                         2.0 * M_PI * azimuth, candidate.outgoing.direction);
        candidate.deposition_available = true;
        candidate.local_energy_deposition =
            evaluation->incident.energy - interaction.energy_out;
    } else {
        double absorption = microscopic_total - elastic;
        if (absorption <= 0.0 || component->n_absorption_reactions <= 0)
            return ALEA_ERR_INVALID_STATE;
        double residual = reaction_threshold - elastic;
        double sum = 0.0;
        int reaction_index = component->absorption_reactions[
            component->n_absorption_reactions - 1];
        for (int i = 0; i < component->n_absorption_reactions; i++) {
            int index = component->absorption_reactions[i];
            sum += reaction_xs(nuc, &nuc->reactions[index],
                               evaluation->incident.energy);
            if (sum >= residual) {
                reaction_index = index;
                break;
            }
        }
        candidate.outcome = ALEA_NUC_OUTCOME_ABSORBED;
        candidate.mt = nuc->reactions[reaction_index].mt;
        candidate.deposition_available = false;
        candidate.local_energy_deposition = NAN;
    }

    *result = candidate;
    return ALEA_OK;
}
