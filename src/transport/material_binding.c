// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

#include "alea_transport.h"
#include "core/alea_system.h"
#include "core/alea_cell.h"
#include "core/alea_materials.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define AVOGADRO_BARN_CM 0.602214076
#define KELVIN_TO_MEV 8.617333262145e-11

typedef struct {
    alea_nuc_material_t* material;
    alea_nuc_prepared_material_t* prepared;
} bound_particle_t;

typedef struct {
    bound_particle_t neutron;
    bound_particle_t photon;
} bound_cell_t;

typedef struct {
    char zaid[24];
    alea_nuc_thermal_t* table;
} cached_thermal_t;

struct alea_nuc_cell_bindings {
    const alea_system_t* system;
    uint64_t geometry_generation;
    size_t cell_count;
    bound_cell_t* cells;
    alea_nuc_binding_notice_t* notices;
    size_t notice_count;
    size_t notice_capacity;
    alea_nuc_binding_selection_t* selections;
    size_t selection_count;
    size_t selection_capacity;
    cached_thermal_t* thermals;
    size_t thermal_count;
    size_t thermal_capacity;
};

typedef struct {
    int zaid;
    const char* library;
    double fraction;
    double mass;
    int source_element; /* -1 for explicit nuclides */
    double atom_abundance;
    double element_fraction;
    bool omitted;
} isotope_t;

static alea_error_t add_notice(alea_nuc_cell_bindings_t* bindings,
                               const alea_nuc_binding_notice_t* notice) {
    if (bindings->notice_count == bindings->notice_capacity) {
        size_t capacity = bindings->notice_capacity ? bindings->notice_capacity * 2 : 8;
        if (capacity < bindings->notice_capacity ||
            capacity > SIZE_MAX / sizeof(*bindings->notices))
            return ALEA_ERR_OVERFLOW;
        void* next = realloc(bindings->notices, capacity * sizeof(*bindings->notices));
        if (!next) return ALEA_ERR_OUT_OF_MEMORY;
        bindings->notices = next;
        bindings->notice_capacity = capacity;
    }
    bindings->notices[bindings->notice_count++] = *notice;
    return ALEA_OK;
}

static alea_error_t add_selection(alea_nuc_cell_bindings_t* bindings,
                                  const alea_nuc_binding_selection_t* selection) {
    if (bindings->selection_count == bindings->selection_capacity) {
        size_t capacity = bindings->selection_capacity ?
            bindings->selection_capacity * 2 : 8;
        if (capacity < bindings->selection_capacity ||
            capacity > SIZE_MAX / sizeof(*bindings->selections))
            return ALEA_ERR_OVERFLOW;
        void* next = realloc(bindings->selections,
                             capacity * sizeof(*bindings->selections));
        if (!next) return ALEA_ERR_OUT_OF_MEMORY;
        bindings->selections = next;
        bindings->selection_capacity = capacity;
    }
    bindings->selections[bindings->selection_count++] = *selection;
    return ALEA_OK;
}

size_t alea_nuc_cell_bindings_selection_count(
    const alea_nuc_cell_bindings_t* bindings) {
    return bindings ? bindings->selection_count : 0;
}

const alea_nuc_binding_selection_t* alea_nuc_cell_bindings_selection(
    const alea_nuc_cell_bindings_t* bindings, size_t index) {
    return bindings && index < bindings->selection_count ?
        &bindings->selections[index] : NULL;
}

size_t alea_nuc_cell_bindings_notice_count(const alea_nuc_cell_bindings_t* bindings) {
    return bindings ? bindings->notice_count : 0;
}

const alea_nuc_binding_notice_t* alea_nuc_cell_bindings_notice(
    const alea_nuc_cell_bindings_t* bindings, size_t index) {
    return bindings && index < bindings->notice_count ?
        &bindings->notices[index] : NULL;
}

static void free_particle(bound_particle_t* part) {
    alea_nuc_prepared_material_free(part->prepared);
    alea_nuc_material_destroy(part->material);
}

void alea_nuc_cell_bindings_free(alea_nuc_cell_bindings_t* bindings) {
    if (!bindings) return;
    for (size_t i = 0; i < bindings->cell_count; ++i) {
        free_particle(&bindings->cells[i].neutron);
        free_particle(&bindings->cells[i].photon);
    }
    free(bindings->cells);
    free(bindings->notices);
    free(bindings->selections);
    for (size_t i = 0; i < bindings->thermal_count; ++i)
        alea_nuc_thermal_free(bindings->thermals[i].table);
    free(bindings->thermals);
    free(bindings);
}

const alea_nuc_prepared_material_t* alea_nuc_cell_bindings_get(
    const alea_nuc_cell_bindings_t* bindings, size_t cell_index,
    alea_nuc_particle_t particle) {
    if (!bindings || cell_index >= bindings->cell_count ||
        alea_system_geometry_generation(bindings->system) !=
            bindings->geometry_generation)
        return NULL;
    const bound_cell_t* cell = &bindings->cells[cell_index];
    if (particle == ALEA_NUC_PARTICLE_NEUTRON) return cell->neutron.prepared;
    if (particle == ALEA_NUC_PARTICLE_PHOTON) return cell->photon.prepared;
    return NULL;
}

bool alea_nuc_cell_bindings_matches_geometry(
    const alea_nuc_cell_bindings_t* bindings, const alea_system_t* sys) {
    return bindings && sys && bindings->system == sys &&
        bindings->cell_count == sys->cells.count &&
        bindings->geometry_generation == alea_system_geometry_generation(sys);
}

static double isotope_mass(int zaid) {
    int z = alea_zaid_to_Z(zaid);
    int a = alea_zaid_to_A(zaid);
    const alea_element_t* element = alea_get_element(z);
    if (element) {
        for (size_t i = 0; i < element->isotope_count; ++i)
            if (element->isotopes[i].mass_number == a)
                return element->isotopes[i].atomic_mass;
        if (a == 0 && element->standard_atomic_weight > 0.0)
            return element->standard_atomic_weight;
    }
    return a > 0 ? (double)a : 0.0;
}

static alea_error_t collect_isotopes(const alea_material_t* mat,
                                     isotope_t** output, size_t* count) {
    size_t capacity = mat->nuclides.count;
    for (size_t i = 0; i < mat->elements.count; ++i) {
        const alea_element_t* e = alea_get_element(mat->elements.data[i].atomic_number);
        size_t extra = e && e->isotope_count ? e->isotope_count : 1;
        if (capacity > SIZE_MAX - extra) return ALEA_ERR_OVERFLOW;
        capacity += extra;
    }
    if (capacity == 0 || capacity > SIZE_MAX / sizeof(isotope_t))
        return ALEA_ERR_EMPTY;
    isotope_t* list = calloc(capacity, sizeof(*list));
    if (!list) return ALEA_ERR_OUT_OF_MEMORY;
    size_t n = 0;
    for (size_t i = 0; i < mat->nuclides.count; ++i) {
        const alea_nuclide_t* nu = &mat->nuclides.data[i];
        list[n++] = (isotope_t){.zaid = nu->zaid, .library = nu->library,
            .fraction = nu->fraction, .mass = isotope_mass(nu->zaid),
            .source_element = -1};
    }
    for (size_t i = 0; i < mat->elements.count; ++i) {
        const alea_element_comp_t* ec = &mat->elements.data[i];
        const alea_element_t* e = alea_get_element(ec->atomic_number);
        if (!e || e->isotope_count == 0) {
            int zaid = ec->atomic_number * 1000;
            list[n++] = (isotope_t){.zaid = zaid, .library = ec->library,
                .fraction = ec->fraction, .mass = isotope_mass(zaid),
                .source_element = (int)i, .atom_abundance = 1.0,
                .element_fraction = ec->fraction};
            continue;
        }
        double mean_mass = 0.0;
        for (size_t j = 0; j < e->isotope_count; ++j)
            if (e->isotopes[j].abundance > 0.0)
                mean_mass += e->isotopes[j].abundance *
                             e->isotopes[j].atomic_mass;
        if (mat->is_weight_fraction &&
            (!isfinite(mean_mass) || mean_mass <= 0.0)) {
            free(list);
            return ALEA_ERR_INVALID_ARG;
        }
        for (size_t j = 0; j < e->isotope_count; ++j) {
            const alea_isotope_t* iso = &e->isotopes[j];
            if (iso->abundance <= 0.0) continue;
            double fraction = ec->fraction * iso->abundance;
            if (mat->is_weight_fraction)
                fraction *= iso->atomic_mass / mean_mass;
            list[n++] = (isotope_t){
                .zaid = ec->atomic_number * 1000 + iso->mass_number,
                .library = ec->library, .fraction = fraction,
                .mass = iso->atomic_mass, .source_element = (int)i,
                .atom_abundance = iso->abundance,
                .element_fraction = ec->fraction};
        }
    }
    *output = list;
    *count = n;
    return n ? ALEA_OK : ALEA_ERR_EMPTY;
}

static const alea_nuc_xsdir_entry_t* find_unique(
    const alea_nuc_xsdir_t* xsdir, int numeric_zaid,
    alea_nuc_table_type_t type) {
    const alea_nuc_xsdir_entry_t* found = NULL;
    for (size_t i = 0; i < xsdir->count; ++i) {
        const alea_nuc_xsdir_entry_t* e = &xsdir->entries[i];
        if (e->type != type || atoi(e->zaid) != numeric_zaid) continue;
        if (found) return NULL;
        found = e;
    }
    return found;
}

typedef struct {
    const alea_nuc_xsdir_entry_t* lower;
    const alea_nuc_xsdir_entry_t* upper;
    double upper_fraction;
    alea_nuc_binding_notice_kind_t notice;
    bool has_notice;
    bool missing_table;
} selected_neutron_t;

static alea_error_t select_neutron(const alea_nuc_xsdir_t* xsdir,
    const isotope_t* iso, const alea_cell_entry_t* cell,
    const alea_nuc_binding_policy_t* policy, selected_neutron_t* selected) {
    memset(selected, 0, sizeof(*selected));
    const bool explicit_library = iso->library && iso->library[0];
    char anchor[32];
    if (explicit_library) {
        int len = snprintf(anchor, sizeof(anchor), "%d%s%s", iso->zaid,
            iso->library[0] == '.' ? "" : ".", iso->library);
        if (len < 0 || (size_t)len >= sizeof(anchor)) return ALEA_ERR_INVALID_ARG;
        selected->lower = alea_nuc_xsdir_find(xsdir, anchor);
        if (!selected->lower ||
            selected->lower->type != ALEA_NUC_TABLE_CONTINUOUS_NEUTRON) {
            selected->missing_table = true;
            return ALEA_ERR_NOT_FOUND;
        }
    } else {
        const alea_nuc_xsdir_entry_t* first = NULL;
        int matches = 0;
        for (size_t i = 0; i < xsdir->count; ++i) {
            const alea_nuc_xsdir_entry_t* e = &xsdir->entries[i];
            if (e->type != ALEA_NUC_TABLE_CONTINUOUS_NEUTRON ||
                atoi(e->zaid) != iso->zaid) continue;
            if (!first) first = e;
            matches++;
        }
        if (!first) {
            selected->missing_table = true;
            return ALEA_ERR_NOT_FOUND;
        }
        if (!cell->has_temperature && matches != 1) return ALEA_ERR_INVALID_STATE;
        selected->lower = first;
        snprintf(anchor, sizeof(anchor), "%s", first->zaid);
    }
    selected->upper = selected->lower;
    if (!cell->has_temperature) return ALEA_OK;
    if (!isfinite(cell->temperature) || cell->temperature <= 0.0)
        return ALEA_ERR_INVALID_ARG;
    const double requested = cell->temperature * KELVIN_TO_MEV;
    if (explicit_library) {
        double delta_kelvin = fabs(selected->lower->temperature - requested) /
                              KELVIN_TO_MEV;
        if (delta_kelvin <= 1e-6) return ALEA_OK;
        if (policy->nearest_temperature_tolerance > 0.0 &&
            delta_kelvin <= policy->nearest_temperature_tolerance) {
            selected->has_notice = true;
            selected->notice = ALEA_NUC_BINDING_TEMPERATURE_NEAREST;
            return ALEA_OK;
        }
        return ALEA_ERR_NOT_FOUND;
    }
    alea_error_t err = alea_nuc_xsdir_find_temperature_bracket(xsdir,
        anchor, requested, &selected->lower, &selected->upper,
        &selected->upper_fraction);
    if (err == ALEA_OK) {
        if (selected->lower != selected->upper) {
            selected->has_notice = true;
            selected->notice = ALEA_NUC_BINDING_TEMPERATURE_INTERPOLATED;
        }
        return ALEA_OK;
    }
    if (err != ALEA_ERR_NOT_FOUND ||
        policy->nearest_temperature_tolerance == 0.0) return err;
    const alea_nuc_xsdir_entry_t* nearest = NULL;
    err = alea_nuc_xsdir_find_temperature(xsdir, anchor, requested,
        policy->nearest_temperature_tolerance * KELVIN_TO_MEV, &nearest);
    if (err != ALEA_OK) return err;
    selected->lower = selected->upper = nearest;
    selected->upper_fraction = 0.0;
    selected->has_notice = true;
    selected->notice = ALEA_NUC_BINDING_TEMPERATURE_NEAREST;
    return ALEA_OK;
}

static alea_error_t select_thermal(const alea_nuc_xsdir_t* xsdir,
    const char* identifier, double neutron_kT,
    const alea_nuc_binding_policy_t* policy,
    const alea_nuc_xsdir_entry_t** output, bool* nearest_notice) {
    *output = NULL;
    *nearest_notice = false;
    if (!identifier || !identifier[0] || !isfinite(neutron_kT) ||
        neutron_kT <= 0.0) return ALEA_ERR_INVALID_ARG;
    const bool fixed = strchr(identifier, '.') != NULL;
    const alea_nuc_xsdir_entry_t* best = NULL;
    double best_delta = HUGE_VAL;
    bool ambiguous = false;
    size_t prefix_length = strlen(identifier);
    for (size_t i = 0; i < xsdir->count; ++i) {
        const alea_nuc_xsdir_entry_t* entry = &xsdir->entries[i];
        if (entry->type != ALEA_NUC_TABLE_THERMAL_SAB ||
            !isfinite(entry->temperature) || entry->temperature <= 0.0)
            continue;
        if (fixed ? strcmp(identifier, entry->zaid) != 0 :
            (strncmp(identifier, entry->zaid, prefix_length) != 0 ||
             entry->zaid[prefix_length] != '.')) continue;
        double delta = fabs(entry->temperature - neutron_kT);
        if (!best || delta < best_delta - 1e-18) {
            best = entry;
            best_delta = delta;
            ambiguous = false;
        } else if (fabs(delta - best_delta) <= 1e-18) {
            ambiguous = true;
        }
    }
    if (!best) return ALEA_ERR_NOT_FOUND;
    if (ambiguous) return ALEA_ERR_INVALID_STATE;
    double delta_kelvin = best_delta / KELVIN_TO_MEV;
    if (delta_kelvin > 1e-6) {
        if (policy->nearest_temperature_tolerance == 0.0 ||
            delta_kelvin > policy->nearest_temperature_tolerance)
            return ALEA_ERR_NOT_FOUND;
        *nearest_notice = true;
    }
    *output = best;
    return ALEA_OK;
}

static alea_nuc_thermal_t* get_thermal(alea_nuc_cell_bindings_t* bindings,
    const alea_nuc_xsdir_t* xsdir, const char* zaid) {
    for (size_t i = 0; i < bindings->thermal_count; ++i)
        if (strcmp(bindings->thermals[i].zaid, zaid) == 0)
            return bindings->thermals[i].table;
    alea_nuc_thermal_t* table = alea_nuc_load_thermal(xsdir, zaid);
    if (!table) return NULL;
    if (bindings->thermal_count == bindings->thermal_capacity) {
        size_t capacity = bindings->thermal_capacity ?
            bindings->thermal_capacity * 2 : 4;
        if (capacity < bindings->thermal_capacity ||
            capacity > SIZE_MAX / sizeof(*bindings->thermals)) {
            alea_nuc_thermal_free(table);
            return NULL;
        }
        void* next = realloc(bindings->thermals,
                             capacity * sizeof(*bindings->thermals));
        if (!next) {
            alea_nuc_thermal_free(table);
            return NULL;
        }
        bindings->thermals = next;
        bindings->thermal_capacity = capacity;
    }
    cached_thermal_t* cached = &bindings->thermals[bindings->thermal_count++];
    snprintf(cached->zaid, sizeof(cached->zaid), "%s", zaid);
    cached->table = table;
    return table;
}

static bool thermal_applies(const alea_nuc_thermal_t* thermal,
                            const alea_nuc_nuclide_t* nuc) {
    int zaid = 1000 * nuc->Z + nuc->A;
    for (int i = 0; i < thermal->n_applicable_zaids; ++i) {
        int listed = thermal->applicable_zaids[i];
        if (listed == zaid) return true;
        if (listed / 1000 != nuc->Z || nuc->Z == 1 || nuc->Z == 26)
            continue;
        const alea_element_t* element = alea_get_element(nuc->Z);
        if (!element) continue;
        for (size_t j = 0; j < element->isotope_count; ++j)
            if (element->isotopes[j].mass_number == nuc->A &&
                element->isotopes[j].abundance > 0.0)
                return true;
    }
    return false;
}

static const char* resolve_thermal_name(
    const alea_nuc_binding_policy_t* policy, const char* name) {
    if (!name) return NULL;
    for (size_t i = 0; i < policy->thermal_name_map_count; ++i)
        if (strcmp(policy->thermal_name_map[i].material_law, name) == 0)
            return policy->thermal_name_map[i].xsdir_table;
    return name;
}

static alea_error_t bind_cell(alea_nuc_cell_bindings_t* bindings,
    alea_nuc_xsdir_t* xsdir, size_t index, uint32_t particles,
    const alea_nuc_prepare_requirements_t* neutron_requirements,
    const alea_nuc_binding_policy_t* policy) {
    const alea_cell_entry_t* cell = &bindings->system->cells.data[index];
    if (cell->material_index < 0 || alea_cell_entry_is_container(cell))
        return ALEA_OK;
    if ((size_t)cell->material_index >= bindings->system->materials.count) {
        alea_set_error_detail(ALEA_ERR_INVALID_STATE,
            "cell %d references invalid material index", cell->mc_cell_id);
        return ALEA_ERR_INVALID_STATE;
    }
    const alea_material_t* mat = &bindings->system->materials.data[cell->material_index];
    isotope_t* isotopes = NULL;
    selected_neutron_t* selected = NULL;
    alea_nuc_thermal_association_t* associations = NULL;
    size_t association_count = 0;
    size_t count = 0;
    char context[96] = "composition or density";
    alea_error_t err = collect_isotopes(mat, &isotopes, &count);
    if (err != ALEA_OK) goto fail;
    if (particles & ALEA_NUC_BIND_NEUTRON) {
        selected = calloc(count, sizeof(*selected));
        if (!selected) { err = ALEA_ERR_OUT_OF_MEMORY; goto fail; }
        for (size_t i = 0; i < count; ++i) {
            isotope_t* iso = &isotopes[i];
            if (iso->fraction == 0.0) continue;
            snprintf(context, sizeof(context), "neutron ZAID %d", iso->zaid);
            err = select_neutron(xsdir, iso, cell, policy, &selected[i]);
            if (err == ALEA_OK) continue;
            if (err == ALEA_ERR_NOT_FOUND && selected[i].missing_table &&
                iso->source_element >= 0 &&
                policy->max_omitted_natural_abundance > 0.0) {
                iso->omitted = true;
                continue;
            }
            goto fail;
        }
        for (size_t e = 0; e < mat->elements.count; ++e) {
            double original = 0.0, retained = 0.0, retained_mass = 0.0;
            for (size_t i = 0; i < count; ++i) {
                if (isotopes[i].source_element != (int)e) continue;
                original += isotopes[i].atom_abundance;
                if (!isotopes[i].omitted) {
                    retained += isotopes[i].atom_abundance;
                    retained_mass += isotopes[i].atom_abundance * isotopes[i].mass;
                }
            }
            if (!(original > 0.0)) continue;
            double omitted_fraction = (original - retained) / original;
            if (omitted_fraction > policy->max_omitted_natural_abundance + 1e-15 ||
                !(retained > 0.0) || !(retained_mass > 0.0)) {
                snprintf(context, sizeof(context),
                    "natural element Z %d omits %.6g atom fraction (limit %.6g)",
                    mat->elements.data[e].atomic_number, omitted_fraction,
                    policy->max_omitted_natural_abundance);
                err = ALEA_ERR_NOT_FOUND;
                goto fail;
            }
            for (size_t i = 0; i < count; ++i) {
                isotope_t* iso = &isotopes[i];
                if (iso->source_element != (int)e) continue;
                if (iso->omitted) {
                    iso->fraction = 0.0;
                    alea_nuc_binding_notice_t notice = {
                        .kind = ALEA_NUC_BINDING_NATURAL_ISOTOPE_OMITTED,
                        .table_type = ALEA_NUC_TABLE_CONTINUOUS_NEUTRON,
                        .cell_index = index, .cell_id = cell->mc_cell_id,
                        .material_id = mat->material_id, .zaid = iso->zaid,
                        .omitted_abundance = iso->atom_abundance / original
                    };
                    err = add_notice(bindings, &notice);
                    if (err != ALEA_OK) goto fail;
                } else {
                    double share = iso->atom_abundance / retained;
                    if (mat->is_weight_fraction)
                        share = iso->atom_abundance * iso->mass / retained_mass;
                    iso->fraction = iso->element_fraction * share;
                }
            }
        }
    }

    double density = fabs(cell->density);
    bool mass_density = cell->is_mass_density;
    if (density == 0.0 && mat->has_standard_density) {
        density = fabs(mat->standard_density);
        mass_density = mat->standard_density > 0.0;
    }
    if (!isfinite(density) || density <= 0.0) {
        err = ALEA_ERR_INVALID_STATE;
        goto fail;
    }
    double fraction_sum = 0.0, weighted_mass = 0.0, weight_over_mass = 0.0;
    for (size_t i = 0; i < count; ++i) {
        isotope_t* iso = &isotopes[i];
        if (!isfinite(iso->fraction) || iso->fraction < 0.0 ||
            !isfinite(iso->mass) || iso->mass <= 0.0) {
            err = ALEA_ERR_INVALID_ARG;
            goto fail;
        }
        fraction_sum += iso->fraction;
        weighted_mass += iso->fraction * iso->mass;
        weight_over_mass += iso->fraction / iso->mass;
    }
    if (!isfinite(fraction_sum) || fraction_sum <= 0.0 ||
        !isfinite(weighted_mass) || !isfinite(weight_over_mass)) {
        err = ALEA_ERR_INVALID_ARG;
        goto fail;
    }

    bound_cell_t* bound = &bindings->cells[index];
    if (particles & ALEA_NUC_BIND_NEUTRON) {
        bound->neutron.material = alea_nuc_material_create();
        if (!bound->neutron.material) { err = ALEA_ERR_OUT_OF_MEMORY; goto fail; }
    }
    if (particles & ALEA_NUC_BIND_PHOTON) {
        bound->photon.material = alea_nuc_material_create();
        if (!bound->photon.material) { err = ALEA_ERR_OUT_OF_MEMORY; goto fail; }
    }
    for (size_t i = 0; i < count; ++i) {
        const isotope_t* iso = &isotopes[i];
        if (iso->fraction == 0.0) continue;
        double ni;
        if (mat->is_weight_fraction) {
            double share = (iso->fraction / iso->mass) / weight_over_mass;
            ni = mass_density ? density * AVOGADRO_BARN_CM *
                    iso->fraction / fraction_sum / iso->mass : density * share;
        } else {
            ni = mass_density ? density * AVOGADRO_BARN_CM *
                    iso->fraction / weighted_mass : density *
                    iso->fraction / fraction_sum;
        }
        if (!isfinite(ni) || ni <= 0.0) { err = ALEA_ERR_INVALID_ARG; goto fail; }

        if (particles & ALEA_NUC_BIND_NEUTRON) {
            snprintf(context, sizeof(context), "neutron ZAID %d", iso->zaid);
            const selected_neutron_t* choice = &selected[i];
            const alea_nuc_xsdir_entry_t* low = choice->lower;
            const alea_nuc_xsdir_entry_t* high = choice->upper;
            alea_nuc_nuclide_t* lower =
                alea_nuc_xsdir_get_nuclide(xsdir, low->zaid);
            alea_nuc_nuclide_t* upper = low == high ? lower :
                alea_nuc_xsdir_get_nuclide(xsdir, high->zaid);
            if (!lower || !upper) { err = ALEA_ERR_FILE_READ; goto fail; }
            err = alea_nuc_material_add_temperature_mix(
                bound->neutron.material, lower, upper, choice->upper_fraction, ni);
            if (err != ALEA_OK) goto fail;
            alea_nuc_binding_selection_t selection = {
                .cell_index = index, .cell_id = cell->mc_cell_id,
                .material_id = mat->material_id,
                .particle = ALEA_NUC_PARTICLE_NEUTRON,
                .table_type = ALEA_NUC_TABLE_CONTINUOUS_NEUTRON,
                .zaid = iso->zaid,
                .requested_kelvin = cell->has_temperature ? cell->temperature : 0.0,
                .selected_kelvin = low->temperature / KELVIN_TO_MEV,
                .upper_kelvin = high->temperature / KELVIN_TO_MEV,
                .upper_fraction = choice->upper_fraction,
                .number_density = ni
            };
            snprintf(selection.table, sizeof(selection.table), "%s", low->zaid);
            if (low != high)
                snprintf(selection.upper_table, sizeof(selection.upper_table),
                         "%s", high->zaid);
            err = add_selection(bindings, &selection);
            if (err != ALEA_OK) goto fail;
            if (choice->has_notice) {
                alea_nuc_binding_notice_t notice = {
                    .kind = choice->notice,
                    .table_type = ALEA_NUC_TABLE_CONTINUOUS_NEUTRON,
                    .cell_index = index,
                    .cell_id = cell->mc_cell_id,
                    .material_id = mat->material_id, .zaid = iso->zaid,
                    .requested_kelvin = cell->temperature,
                    .selected_kelvin = low->temperature / KELVIN_TO_MEV,
                    .upper_kelvin = high->temperature / KELVIN_TO_MEV,
                    .upper_fraction = choice->upper_fraction
                };
                snprintf(notice.table, sizeof(notice.table), "%s", low->zaid);
                if (low != high)
                    snprintf(notice.upper_table, sizeof(notice.upper_table),
                             "%s", high->zaid);
                err = add_notice(bindings, &notice);
                if (err != ALEA_OK) goto fail;
            }
        }
        if (particles & ALEA_NUC_BIND_PHOTON) {
            int z = alea_zaid_to_Z(iso->zaid);
            snprintf(context, sizeof(context), "photoatomic Z %d", z);
            const alea_nuc_xsdir_entry_t* photo = find_unique(xsdir,
                z * 1000, ALEA_NUC_TABLE_PHOTOATOMIC);
            if (!photo) { err = ALEA_ERR_NOT_FOUND; goto fail; }
            alea_nuc_nuclide_t* table =
                alea_nuc_xsdir_get_nuclide(xsdir, photo->zaid);
            if (!table) { err = ALEA_ERR_FILE_READ; goto fail; }
            err = alea_nuc_material_add(bound->photon.material, table, ni);
            if (err != ALEA_OK) goto fail;
            alea_nuc_binding_selection_t selection = {
                .cell_index = index, .cell_id = cell->mc_cell_id,
                .material_id = mat->material_id,
                .particle = ALEA_NUC_PARTICLE_PHOTON,
                .table_type = ALEA_NUC_TABLE_PHOTOATOMIC,
                .zaid = iso->zaid,
                .number_density = ni
            };
            snprintf(selection.table, sizeof(selection.table), "%s", photo->zaid);
            err = add_selection(bindings, &selection);
            if (err != ALEA_OK) goto fail;
        }
    }
    if (bound->neutron.material) {
        snprintf(context, sizeof(context), "neutron collision capability");
        alea_nuc_prepare_requirements_t default_req = {
            .required_capabilities = ALEA_NUC_CAP_CONTINUOUS_NEUTRON
        };
        alea_nuc_prepare_requirements_t req = neutron_requirements ?
            *neutron_requirements : default_req;
        if (mat->thermal_laws.count) {
            size_t external_count = req.n_thermal_associations;
            size_t components = (size_t)bound->neutron.material->n_components;
            if (mat->thermal_laws.count &&
                components > (SIZE_MAX - external_count) / mat->thermal_laws.count) {
                err = ALEA_ERR_OVERFLOW;
                goto fail;
            }
            size_t capacity = external_count +
                              components * mat->thermal_laws.count;
            if (capacity > SIZE_MAX / sizeof(*associations)) {
                err = ALEA_ERR_OVERFLOW;
                goto fail;
            }
            associations = calloc(capacity ? capacity : 1,
                                  sizeof(*associations));
            if (!associations) { err = ALEA_ERR_OUT_OF_MEMORY; goto fail; }
            if (external_count) {
                if (!req.thermal_associations) {
                    err = ALEA_ERR_INVALID_ARG;
                    goto fail;
                }
                memcpy(associations, req.thermal_associations,
                       external_count * sizeof(*associations));
                association_count = external_count;
            }
            for (size_t l = 0; l < mat->thermal_laws.count; ++l) {
                const alea_thermal_law_t* law = &mat->thermal_laws.data[l];
                const char* thermal_name =
                    resolve_thermal_name(policy, law->identifier);
                snprintf(context, sizeof(context), "thermal law %.60s",
                         law->identifier ? law->identifier : "(null)");
                const alea_nuc_xsdir_entry_t* anchor = NULL;
                for (size_t j = 0; j < xsdir->count; ++j) {
                    const alea_nuc_xsdir_entry_t* candidate = &xsdir->entries[j];
                    if (candidate->type != ALEA_NUC_TABLE_THERMAL_SAB ||
                        !thermal_name) continue;
                    size_t len = strlen(thermal_name);
                    bool fixed = strchr(thermal_name, '.') != NULL;
                    if (fixed ? strcmp(thermal_name, candidate->zaid) == 0 :
                        (strncmp(thermal_name, candidate->zaid, len) == 0 &&
                         candidate->zaid[len] == '.')) {
                        anchor = candidate;
                        break;
                    }
                }
                if (!anchor) { err = ALEA_ERR_NOT_FOUND; goto fail; }
                alea_nuc_thermal_t* anchor_table = get_thermal(bindings, xsdir,
                                                                 anchor->zaid);
                if (!anchor_table) { err = ALEA_ERR_FILE_READ; goto fail; }
                size_t matches = 0;
                for (size_t c = 0; c < components; ++c) {
                    const alea_nuc_mat_component_t* component =
                        &bound->neutron.material->components[c];
                    alea_nuc_nuclide_t* nuc = component->nuclide;
                    int component_zaid = 1000 * nuc->Z + nuc->A;
                    if (law->zaid_match > 0 &&
                        law->zaid_match != component_zaid) continue;
                    if (!thermal_applies(anchor_table, nuc)) continue;
                    const alea_nuc_xsdir_entry_t* entry = NULL;
                    bool nearest = false;
                    err = select_thermal(xsdir, thermal_name,
                        nuc->temperature, policy, &entry, &nearest);
                    if (err != ALEA_OK) goto fail;
                    alea_nuc_thermal_t* thermal = get_thermal(bindings, xsdir,
                                                                entry->zaid);
                    if (!thermal) { err = ALEA_ERR_FILE_READ; goto fail; }
                    if (!thermal_applies(thermal, nuc)) {
                        err = ALEA_ERR_INVALID_STATE;
                        goto fail;
                    }
                    if (association_count >= capacity) {
                        err = ALEA_ERR_OVERFLOW;
                        goto fail;
                    }
                    associations[association_count++] =
                        (alea_nuc_thermal_association_t){(int)c, thermal};
                    matches++;
                    alea_nuc_binding_selection_t selection = {
                        .cell_index = index, .cell_id = cell->mc_cell_id,
                        .material_id = mat->material_id,
                        .particle = ALEA_NUC_PARTICLE_NEUTRON,
                        .table_type = ALEA_NUC_TABLE_THERMAL_SAB,
                        .zaid = component_zaid,
                        .requested_kelvin = nuc->temperature / KELVIN_TO_MEV,
                        .selected_kelvin = thermal->temperature / KELVIN_TO_MEV,
                        .number_density = component->number_density
                    };
                    snprintf(selection.table, sizeof(selection.table), "%s",
                             entry->zaid);
                    err = add_selection(bindings, &selection);
                    if (err != ALEA_OK) goto fail;
                    if (nearest) {
                        alea_nuc_binding_notice_t notice = {
                            .kind = ALEA_NUC_BINDING_TEMPERATURE_NEAREST,
                            .table_type = ALEA_NUC_TABLE_THERMAL_SAB,
                            .cell_index = index, .cell_id = cell->mc_cell_id,
                            .material_id = mat->material_id,
                            .zaid = component_zaid,
                            .requested_kelvin = nuc->temperature / KELVIN_TO_MEV,
                            .selected_kelvin = thermal->temperature / KELVIN_TO_MEV
                        };
                        snprintf(notice.table, sizeof(notice.table), "%s",
                                 entry->zaid);
                        err = add_notice(bindings, &notice);
                        if (err != ALEA_OK) goto fail;
                    }
                }
                if (matches == 0) {
                    err = ALEA_ERR_NOT_FOUND;
                    goto fail;
                }
            }
            req.required_capabilities |= ALEA_NUC_CAP_THERMAL_SAB;
            req.thermal_associations = associations;
            req.n_thermal_associations = association_count;
            double tolerance = policy->nearest_temperature_tolerance *
                               KELVIN_TO_MEV;
            if (req.thermal_temperature_tolerance < tolerance)
                req.thermal_temperature_tolerance = tolerance;
        }
        alea_nuc_capability_report_t report;
        err = alea_nuc_prepare_material(bound->neutron.material,
            &req,
            &report, &bound->neutron.prepared);
        if (err != ALEA_OK) {
            if (report.detail[0])
                snprintf(context, sizeof(context), "neutron: %.80s", report.detail);
            goto fail;
        }
    }
    if (bound->photon.material) {
        snprintf(context, sizeof(context), "photon collision capability");
        alea_nuc_prepare_requirements_t photon_req = {
            .required_capabilities = ALEA_NUC_CAP_PHOTON
        };
        alea_nuc_capability_report_t report;
        err = alea_nuc_prepare_material(bound->photon.material,
            &photon_req, &report, &bound->photon.prepared);
        if (err != ALEA_OK) {
            if (report.detail[0])
                snprintf(context, sizeof(context), "photon: %.81s", report.detail);
            goto fail;
        }
    }
    free(associations);
    free(selected);
    free(isotopes);
    return ALEA_OK;
fail:
    free(associations);
    free(selected);
    free(isotopes);
    alea_set_error_detail(err,
        "transport binding failed for cell %d / material %d: %s",
        cell->mc_cell_id, mat->material_id, context);
    return err;
}

alea_error_t alea_nuc_cell_bindings_prepare(
    alea_system_t* sys, alea_nuc_xsdir_t* xsdir, uint32_t particles,
    const alea_nuc_prepare_requirements_t* neutron_requirements,
    const alea_nuc_binding_policy_t* requested_policy,
    alea_nuc_cell_bindings_t** output) {
    if (!output) return ALEA_ERR_NULL_ARG;
    *output = NULL;
    if (!sys || !xsdir) return ALEA_ERR_NULL_ARG;
    alea_nuc_binding_policy_t strict = {0};
    const alea_nuc_binding_policy_t* policy = requested_policy ?
        requested_policy : &strict;
    if (!isfinite(policy->nearest_temperature_tolerance) ||
        policy->nearest_temperature_tolerance < 0.0 ||
        !isfinite(policy->max_omitted_natural_abundance) ||
        policy->max_omitted_natural_abundance < 0.0 ||
        policy->max_omitted_natural_abundance >= 1.0 ||
        (policy->thermal_name_map_count && !policy->thermal_name_map))
        return ALEA_ERR_INVALID_ARG;
    for (size_t i = 0; i < policy->thermal_name_map_count; ++i) {
        const alea_nuc_thermal_name_map_t* entry = &policy->thermal_name_map[i];
        if (!entry->material_law || !entry->material_law[0] ||
            !entry->xsdir_table || !entry->xsdir_table[0])
            return ALEA_ERR_INVALID_ARG;
        for (size_t j = 0; j < i; ++j)
            if (strcmp(entry->material_law,
                       policy->thermal_name_map[j].material_law) == 0)
                return ALEA_ERR_INVALID_ARG;
    }
    if (!(particles & (ALEA_NUC_BIND_NEUTRON | ALEA_NUC_BIND_PHOTON)) ||
        (particles & ~(ALEA_NUC_BIND_NEUTRON | ALEA_NUC_BIND_PHOTON)))
        return ALEA_ERR_INVALID_ARG;
    alea_nuc_cell_bindings_t* bindings = calloc(1, sizeof(*bindings));
    if (!bindings) return ALEA_ERR_OUT_OF_MEMORY;
    bindings->system = sys;
    bindings->geometry_generation = alea_system_geometry_generation(sys);
    bindings->cell_count = sys->cells.count;
    bindings->cells = calloc(bindings->cell_count ? bindings->cell_count : 1,
                             sizeof(*bindings->cells));
    if (!bindings->cells) {
        free(bindings);
        return ALEA_ERR_OUT_OF_MEMORY;
    }
    for (size_t i = 0; i < bindings->cell_count; ++i) {
        alea_error_t err = bind_cell(bindings, xsdir, i, particles,
                                    neutron_requirements, policy);
        if (err != ALEA_OK) {
            alea_nuc_cell_bindings_free(bindings);
            return err;
        }
    }
    if (alea_system_geometry_generation(sys) != bindings->geometry_generation) {
        alea_nuc_cell_bindings_free(bindings);
        return ALEA_ERR_INVALID_STATE;
    }
    *output = bindings;
    return ALEA_OK;
}
