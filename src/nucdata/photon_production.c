// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/**
 * @file photon_production.c
 * @brief Neutron-induced photon-production ACE data
 */

#include "nuclear_internal.h"

#include <math.h>
#include <stdlib.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define ALEA_NUC_PHOTON_YIELD_TOL 1e-12

static alea_error_t draw_uniform(alea_nuc_random_fn random, void* context,
                                 double* value) {
    double sampled = random(context);
    if (!isfinite(sampled) || sampled < 0.0 || sampled >= 1.0)
        return ALEA_ERR_INVALID_STATE;
    *value = sampled;
    return ALEA_OK;
}

static alea_error_t decode_tabulated_yield(
    const alea_nuc_ace_table_t* table, int start,
    alea_nuc_photon_production_t* production) {
    int nr = xss_int(table, start);
    if (nr < 0 || nr > 100000 ||
        !xss_range_valid(table, start + 1, 2 * nr + 1))
        return ALEA_ERR_PARSE_ERROR;
    int ne_pos = start + 1 + 2 * nr;
    int ne = xss_int(table, ne_pos);
    if (ne < 2 || ne > 100000 ||
        !xss_range_valid(table, ne_pos + 1, 2 * ne))
        return ALEA_ERR_PARSE_ERROR;

    production->n_regions = nr;
    production->n_energies = ne;
    if (nr > 0) {
        production->nbt = alea_nuc_malloc((size_t)nr * sizeof(int));
        production->interp = alea_nuc_malloc((size_t)nr * sizeof(int));
        if (!production->nbt || !production->interp)
            return ALEA_ERR_OUT_OF_MEMORY;
        for (int i = 0; i < nr; i++) {
            production->nbt[i] = xss_int(table, start + 1 + i);
            production->interp[i] = xss_int(table, start + 1 + nr + i);
        }
    }
    production->energy = xss_copy(table, ne_pos + 1, ne);
    production->values = xss_copy(table, ne_pos + 1 + ne, ne);
    if (!production->energy || !production->values)
        return table->allocation_error ? ALEA_ERR_OUT_OF_MEMORY
                                       : ALEA_ERR_PARSE_ERROR;
    if (!alea_nuc_interp_regions_valid(production->nbt, production->interp,
                                       nr, ne))
        return ALEA_ERR_PARSE_ERROR;
    for (int i = 0; i < ne; i++) {
        if (production->values[i] < 0.0 &&
            production->values[i] >= -ALEA_NUC_PHOTON_YIELD_TOL)
            production->values[i] = 0.0;
        if (!isfinite(production->energy[i]) ||
            !isfinite(production->values[i]) || production->values[i] < 0.0 ||
            (i > 0 && production->energy[i] < production->energy[i - 1]))
            return ALEA_ERR_PARSE_ERROR;
    }
    return ALEA_OK;
}

alea_error_t alea_nuc_decode_photon_production(
    alea_nuc_nuclide_t* nuc, const alea_nuc_ace_table_t* table) {
    int count = table->nxs[5]; /* NXS(6): photon-production reactions */
    int gpd = table->jxs[11];
    if (gpd > 0) {
        nuc->total_photon_production_xs =
            xss_copy(table, gpd, nuc->n_energies);
        if (!nuc->total_photon_production_xs)
            return table->allocation_error ? ALEA_ERR_OUT_OF_MEMORY
                                           : ALEA_ERR_PARSE_ERROR;
        for (int i = 0; i < nuc->n_energies; i++)
            if (!isfinite(nuc->total_photon_production_xs[i]) ||
                nuc->total_photon_production_xs[i] < 0.0)
                return ALEA_ERR_PARSE_ERROR;
    }
    if (count == 0) return ALEA_OK;
    int mtrp = table->jxs[12], lsigp = table->jxs[13];
    int sigp = table->jxs[14], landp = table->jxs[15];
    int andp = table->jxs[16], ldlwp = table->jxs[17];
    int dlwp = table->jxs[18], yp = table->jxs[19];
    if (count < 0 || count > 100000 || mtrp <= 0 || lsigp <= 0 ||
        sigp <= 0 || landp <= 0 || ldlwp <= 0 || dlwp <= 0 || yp <= 0 ||
        !xss_range_valid(table, mtrp, count) ||
        !xss_range_valid(table, lsigp, count) ||
        !xss_range_valid(table, landp, count) ||
        !xss_range_valid(table, ldlwp, count))
        return ALEA_ERR_PARSE_ERROR;

    if (!xss_range_valid(table, yp, 1)) return ALEA_ERR_PARSE_ERROR;
    int nyp = xss_int(table, yp);
    if (nyp < 0 || nyp > count || !xss_range_valid(table, yp + 1, nyp))
        return ALEA_ERR_PARSE_ERROR;
    if (nyp > 0) {
        nuc->photon_yield_multipliers = alea_nuc_malloc((size_t)nyp * sizeof(int));
        if (!nuc->photon_yield_multipliers) return ALEA_ERR_OUT_OF_MEMORY;
        nuc->n_photon_yield_multipliers = nyp;
        for (int i = 0; i < nyp; i++) {
            int mt = xss_int(table, yp + 1 + i);
            if (mt <= 0) return ALEA_ERR_PARSE_ERROR;
            nuc->photon_yield_multipliers[i] = mt;
        }
    }

    nuc->photon_productions = alea_nuc_calloc(
        (size_t)count, sizeof(*nuc->photon_productions));
    if (!nuc->photon_productions) return ALEA_ERR_OUT_OF_MEMORY;
    nuc->n_photon_productions = count;

    for (int i = 0; i < count; i++) {
        alea_nuc_photon_production_t* production =
            &nuc->photon_productions[i];
        production->mt = xss_int(table, mtrp + i);
        production->parent_mt = production->mt / 1000;
        if (production->mt <= 0 || production->parent_mt <= 0)
            return ALEA_ERR_PARSE_ERROR;

        int locator = xss_int(table, lsigp + i);
        int position = xss_relative_loc(table, sigp, locator);
        if (position == 0 || !xss_range_valid(table, position, 3))
            return ALEA_ERR_PARSE_ERROR;
        production->mf = xss_int(table, position++);
        if (production->mf == 12 || production->mf == 16) {
            int multiplier = xss_int(table, position++);
            if (multiplier != production->parent_mt)
                return ALEA_ERR_PARSE_ERROR;
            int listed = 0;
            for (int j = 0; j < nyp; j++)
                if (nuc->photon_yield_multipliers[j] == multiplier) listed = 1;
            if (!listed) return ALEA_ERR_PARSE_ERROR;
            alea_error_t err = decode_tabulated_yield(
                table, position, production);
            if (err != ALEA_OK) return err;
        } else if (production->mf == 13) {
            production->production_xs = true;
            production->threshold_index = xss_int(table, position++);
            production->n_energies = xss_int(table, position++);
            if (production->threshold_index < 1 ||
                production->threshold_index > nuc->n_energies ||
                production->n_energies <= 0 ||
                production->n_energies >
                    nuc->n_energies - production->threshold_index + 1)
                return ALEA_ERR_PARSE_ERROR;
            production->values = xss_copy(
                table, position, production->n_energies);
            if (!production->values)
                return table->allocation_error ? ALEA_ERR_OUT_OF_MEMORY
                                               : ALEA_ERR_PARSE_ERROR;
            for (int j = 0; j < production->n_energies; j++) {
                if (production->values[j] < 0.0 &&
                    production->values[j] >= -ALEA_NUC_PHOTON_YIELD_TOL)
                    production->values[j] = 0.0;
                if (!isfinite(production->values[j]) ||
                    production->values[j] < 0.0)
                    return ALEA_ERR_PARSE_ERROR;
            }
        } else {
            return ALEA_ERR_UNSUPPORTED;
        }

        int spectrum_locator = xss_int(table, ldlwp + i);
        production->spectrum = alea_nuc_decode_energy_dist_base(
            table, spectrum_locator, dlwp);
        if (!production->spectrum)
            return table->allocation_error ? ALEA_ERR_OUT_OF_MEMORY
                                           : ALEA_ERR_PARSE_ERROR;
        int angular_locator = xss_int(table, landp + i);
        if (angular_locator != 0) {
            if (andp <= 0) return ALEA_ERR_PARSE_ERROR;
            production->angular = alea_nuc_decode_angular_base(
                table, angular_locator, andp);
            if (!production->angular)
                return table->allocation_error ? ALEA_ERR_OUT_OF_MEMORY
                                               : ALEA_ERR_PARSE_ERROR;
        }
    }
    return ALEA_OK;
}

double alea_nuc_xs_photon_production_total(
    const alea_nuc_nuclide_t* nuc, double energy) {
    if (!nuc || nuc->particle != ALEA_NUC_PARTICLE_NEUTRON ||
        !isfinite(energy)) return 0.0;
    if (nuc->total_photon_production_xs) {
        double fraction;
        int index = alea_nuc_energy_lookup_trusted(
            nuc->energy, nuc->n_energies, energy, &fraction);
        if (index < 0) return 0.0;
        return nuc->total_photon_production_xs[index] + fraction *
            (nuc->total_photon_production_xs[index + 1] -
             nuc->total_photon_production_xs[index]);
    }
    double total = 0.0;
    for (int i = 0; i < nuc->n_photon_productions; i++) {
        const alea_nuc_photon_production_t* production =
            &nuc->photon_productions[i];
        total += alea_nuc_photon_production_yield(nuc, production, energy) *
                 alea_nuc_xs_reaction(nuc, production->parent_mt, energy);
    }
    return total;
}

static double photon_production_xs_value(
    const alea_nuc_nuclide_t* nuc,
    const alea_nuc_photon_production_t* production, double energy) {
    int start = production->threshold_index - 1;
    double fraction;
    int index = alea_nuc_energy_lookup_trusted(
        nuc->energy, nuc->n_energies, energy, &fraction);
    int local = index - start;
    if (index < 0 || local < 0 || local >= production->n_energies)
        return 0.0;
    double production_xs = production->values[local];
    if (local + 1 < production->n_energies)
        production_xs += fraction *
            (production->values[local + 1] - production->values[local]);
    return production_xs;
}

double alea_nuc_photon_production_yield(
    const alea_nuc_nuclide_t* nuc,
    const alea_nuc_photon_production_t* production,
    double energy) {
    if (!nuc || !production || !isfinite(energy)) return 0.0;
    if (!production->production_xs) {
        double value;
        if (alea_nuc_interp_eval(
                production->energy, production->values,
                production->n_energies, production->nbt,
                production->interp, production->n_regions, energy,
                &value) != ALEA_OK)
            return 0.0;
        return value;
    }
    double parent_xs = alea_nuc_xs_reaction(
        nuc, production->parent_mt, energy);
    return parent_xs > 0.0
        ? photon_production_xs_value(nuc, production, energy) / parent_xs
        : 0.0;
}

double alea_nuc_photon_production_event_yield(
    const alea_nuc_nuclide_t* nuc,
    const alea_nuc_photon_production_t* production,
    int event_mt, double energy) {
    if (!nuc || !production || !isfinite(energy)) return 0.0;
    double event_xs = alea_nuc_xs_reaction(nuc, event_mt, energy);
    if (!(event_xs > 0.0)) return 0.0;
    if (production->production_xs)
        return photon_production_xs_value(nuc, production, energy) / event_xs;
    double parent_xs = alea_nuc_xs_reaction(
        nuc, production->parent_mt, energy);
    return alea_nuc_photon_production_yield(nuc, production, energy) *
           parent_xs / event_xs;
}

alea_error_t alea_nuc_sample_photon_production(
    const alea_nuc_nuclide_t* nuc,
    const alea_nuc_photon_production_t* production,
    const alea_nuc_particle_state_t* incident,
    alea_nuc_random_fn random, void* random_context,
    alea_nuc_particle_state_t* photon) {
    if (!nuc || !production || !incident || !random || !photon)
        return ALEA_ERR_NULL_ARG;
    if (nuc->particle != ALEA_NUC_PARTICLE_NEUTRON ||
        incident->type != ALEA_NUC_PARTICLE_NEUTRON ||
        !production->spectrum)
        return ALEA_ERR_INVALID_ARG;
    int owned = 0;
    for (int i = 0; i < nuc->n_photon_productions; i++)
        if (&nuc->photon_productions[i] == production) owned = 1;
    if (!owned) return ALEA_ERR_INVALID_ARG;
    if (!alea_nuc_validate_angular_internal(production->angular) ||
        alea_nuc_energy_dist_validate(production->spectrum, NULL) != ALEA_OK)
        return ALEA_ERR_INVALID_ARG;

    double energy, mu = 0.0;
    bool correlated = false;
    alea_error_t err = alea_nuc_sample_energy_angle_distribution(
        production->spectrum, incident->energy, random, random_context,
        &energy, &mu, &correlated);
    if (err != ALEA_OK) return err;
    if (!correlated) {
        double select, sample;
        if ((err = draw_uniform(random, random_context, &select)) != ALEA_OK ||
            (err = draw_uniform(random, random_context, &sample)) != ALEA_OK)
            return err;
        mu = alea_nuc_sample_angular_mu_internal(
            production->angular, incident->energy, select, sample);
    }
    double azimuth;
    if ((err = draw_uniform(random, random_context, &azimuth)) != ALEA_OK)
        return err;
    alea_nuc_particle_state_t candidate = *incident;
    candidate.type = ALEA_NUC_PARTICLE_PHOTON;
    candidate.energy = energy;
    alea_nuc_rotate_direction_internal(
        incident->direction, mu, 2.0 * M_PI * azimuth,
        candidate.direction);
    *photon = candidate;
    return ALEA_OK;
}
