// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

#include "tally_internal.h"
#include "core/alea_system.h"
#include "core/alea_cell.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

typedef struct {
    alea_tally_spec_t spec;
    size_t bin_count;
    size_t spatial_bin_count;
    size_t cell_count;
    int* ids;
    size_t* cell_to_bin;
    double width[3];
} tally_definition_t;

struct alea_tally_plan {
    const alea_system_t* system;
    uint64_t geometry_generation;
    size_t count, capacity;
    tally_definition_t* definitions;
    int has_track_scores;
    int has_collision_scores;
    int has_deposition_scores;
    int needs_target;
};

typedef struct {
    alea_tally_view_t view;
    double* sum;
    double* sum_squared;
    double* values;
    uint8_t* used;
    size_t* touched;
    size_t touched_count;
} tally_accumulator_t;

struct alea_tally_results {
    size_t count;
    tally_accumulator_t* tallies;
};

static int compare_int(const void* left, const void* right) {
    int a = *(const int*)left, b = *(const int*)right;
    return (a > b) - (a < b);
}

static size_t find_id(const int* ids, size_t count, int id) {
    size_t lo = 0, hi = count;
    while (lo < hi) {
        size_t mid = lo + (hi - lo) / 2;
        if (ids[mid] < id) lo = mid + 1;
        else hi = mid;
    }
    return lo < count && ids[lo] == id ? lo : SIZE_MAX;
}

static void definition_free(tally_definition_t* definition) {
    free(definition->ids);
    free(definition->cell_to_bin);
    free((void*)definition->spec.energy_edges);
    memset(definition, 0, sizeof(*definition));
}

alea_tally_plan_t* alea_tally_plan_create(const alea_system_t* sys) {
    if (!sys) return NULL;
    alea_tally_plan_t* plan = calloc(1, sizeof(*plan));
    if (!plan) return NULL;
    plan->system = sys;
    plan->geometry_generation = alea_system_geometry_generation(sys);
    return plan;
}

void alea_tally_plan_free(alea_tally_plan_t* plan) {
    if (!plan) return;
    for (size_t i = 0; i < plan->count; ++i)
        definition_free(&plan->definitions[i]);
    free(plan->definitions);
    free(plan);
}

int alea_tally_plan_matches_geometry(const alea_tally_plan_t* plan,
    const alea_system_t* sys) {
    return plan && sys && plan->system == sys &&
        plan->geometry_generation == alea_system_geometry_generation(sys);
}

int alea_tally_plan_has_track_scores(const alea_tally_plan_t* plan) {
    return plan && plan->has_track_scores;
}

int alea_tally_plan_has_collision_scores(const alea_tally_plan_t* plan) {
    return plan && plan->has_collision_scores;
}

int alea_tally_plan_has_deposition_scores(const alea_tally_plan_t* plan) {
    return plan && plan->has_deposition_scores;
}

int alea_tally_plan_needs_target(const alea_tally_plan_t* plan) {
    return plan && plan->needs_target;
}

size_t alea_tally_plan_count(const alea_tally_plan_t* plan) {
    return plan ? plan->count : 0;
}

alea_error_t alea_tally_plan_add(alea_tally_plan_t* plan,
    const alea_tally_spec_t* spec, size_t* index) {
    if (!plan || !spec) return ALEA_ERR_NULL_ARG;
    if (!alea_tally_plan_matches_geometry(plan, plan->system))
        return ALEA_ERR_INVALID_STATE;
    if (spec->score < ALEA_TALLY_TRACK_LENGTH ||
        spec->score > ALEA_TALLY_LOCAL_DEPOSITION ||
        spec->domain < ALEA_TALLY_CELL ||
        spec->domain > ALEA_TALLY_CARTESIAN_MESH ||
        (spec->particle_mask & ~(ALEA_TALLY_NEUTRON | ALEA_TALLY_PHOTON)) ||
        !isfinite(spec->energy_min) || !isfinite(spec->energy_max) ||
        spec->energy_min < 0.0 || spec->energy_max < 0.0 ||
        ((spec->energy_min != 0.0 || spec->energy_max != 0.0) &&
         spec->energy_max <= spec->energy_min) ||
        !isfinite(spec->time_min) || !isfinite(spec->time_max) ||
        ((spec->time_min != 0.0 || spec->time_max != 0.0) &&
         spec->time_max <= spec->time_min) ||
        spec->material_id < -1 || spec->reaction_mt < 0 ||
        spec->nuclide_zaid < 0 ||
        (spec->energy_group_count == 0 && spec->energy_edges) ||
        (spec->energy_group_count > 0 && !spec->energy_edges) ||
        spec->energy_group_count == SIZE_MAX ||
        (spec->score != ALEA_TALLY_REACTION_EVENT &&
         spec->score != ALEA_TALLY_REACTION_RATE &&
         spec->score != ALEA_TALLY_LOCAL_DEPOSITION &&
         spec->reaction_mt != 0) ||
        (spec->score == ALEA_TALLY_TRACK_LENGTH && spec->nuclide_zaid != 0))
        return ALEA_ERR_INVALID_ARG;
    for (size_t group = 0; spec->energy_group_count > 0 &&
                           group <= spec->energy_group_count; ++group) {
        double edge = spec->energy_edges[group];
        if (!isfinite(edge) || edge < 0.0 ||
            (group && edge <= spec->energy_edges[group - 1]))
            return ALEA_ERR_INVALID_ARG;
    }

    tally_definition_t candidate = {.spec = *spec};
    candidate.spec.energy_edges = NULL;
    const size_t cell_count = plan->system->cells.count;
    candidate.cell_count = cell_count;
    if (spec->domain == ALEA_TALLY_CARTESIAN_MESH) {
        candidate.bin_count = 1;
        for (int axis = 0; axis < 3; ++axis) {
            uint32_t n = spec->dimensions[axis];
            double low = spec->lower[axis], high = spec->upper[axis];
            if (!n || !isfinite(low) || !isfinite(high) || high <= low ||
                candidate.bin_count > SIZE_MAX / n) return ALEA_ERR_INVALID_ARG;
            candidate.bin_count *= n;
            candidate.width[axis] = (high - low) / n;
            if (!(candidate.width[axis] > 0.0) ||
                !isfinite(candidate.width[axis])) return ALEA_ERR_INVALID_ARG;
        }
    } else {
        if (cell_count > SIZE_MAX / sizeof(int) ||
            cell_count > SIZE_MAX / sizeof(size_t)) return ALEA_ERR_OVERFLOW;
        candidate.ids = malloc((cell_count ? cell_count : 1) * sizeof(int));
        candidate.cell_to_bin = malloc((cell_count ? cell_count : 1) *
                                       sizeof(size_t));
        if (!candidate.ids || !candidate.cell_to_bin) {
            definition_free(&candidate);
            return ALEA_ERR_OUT_OF_MEMORY;
        }
        if (spec->domain == ALEA_TALLY_CELL) {
            candidate.bin_count = cell_count;
            for (size_t i = 0; i < cell_count; ++i) {
                candidate.ids[i] = plan->system->cells.data[i].mc_cell_id;
                candidate.cell_to_bin[i] = i;
            }
        } else {
            for (size_t i = 0; i < cell_count; ++i)
                candidate.ids[i] = plan->system->cells.data[i].universe_id;
            qsort(candidate.ids, cell_count, sizeof(int), compare_int);
            for (size_t i = 0; i < cell_count; ++i)
                if (i == 0 || candidate.ids[i] != candidate.ids[i - 1])
                    candidate.ids[candidate.bin_count++] = candidate.ids[i];
            for (size_t i = 0; i < cell_count; ++i)
                candidate.cell_to_bin[i] = find_id(candidate.ids,
                    candidate.bin_count,
                    plan->system->cells.data[i].universe_id);
        }
    }
    candidate.spatial_bin_count = candidate.bin_count;
    if (spec->energy_group_count > 0) {
        if (spec->energy_group_count + 1 > SIZE_MAX / sizeof(double) ||
            candidate.bin_count > SIZE_MAX / spec->energy_group_count) {
            definition_free(&candidate);
            return ALEA_ERR_OVERFLOW;
        }
        double* edges = malloc((spec->energy_group_count + 1) *
                               sizeof(double));
        if (!edges) {
            definition_free(&candidate);
            return ALEA_ERR_OUT_OF_MEMORY;
        }
        memcpy(edges, spec->energy_edges,
               (spec->energy_group_count + 1) * sizeof(double));
        candidate.spec.energy_edges = edges;
        candidate.bin_count *= spec->energy_group_count;
    }
    if (candidate.bin_count > SIZE_MAX / sizeof(double) ||
        candidate.bin_count > SIZE_MAX / sizeof(size_t)) {
        definition_free(&candidate);
        return ALEA_ERR_OVERFLOW;
    }
    if (plan->count == plan->capacity) {
        size_t next = plan->capacity ? plan->capacity * 2 : 4;
        if (next < plan->capacity || next > SIZE_MAX / sizeof(candidate)) {
            definition_free(&candidate);
            return ALEA_ERR_OVERFLOW;
        }
        void* resized = realloc(plan->definitions, next * sizeof(candidate));
        if (!resized) {
            definition_free(&candidate);
            return ALEA_ERR_OUT_OF_MEMORY;
        }
        plan->definitions = resized;
        plan->capacity = next;
    }
    if (index) *index = plan->count;
    plan->definitions[plan->count++] = candidate;
    if (spec->score == ALEA_TALLY_TRACK_LENGTH ||
        spec->score == ALEA_TALLY_REACTION_RATE ||
        spec->score == ALEA_TALLY_HEATING) plan->has_track_scores = 1;
    else if (spec->score == ALEA_TALLY_LOCAL_DEPOSITION)
        plan->has_deposition_scores = 1;
    else plan->has_collision_scores = 1;
    if (spec->nuclide_zaid && spec->score != ALEA_TALLY_REACTION_RATE &&
        spec->score != ALEA_TALLY_HEATING)
        plan->needs_target = 1;
    return ALEA_OK;
}

void alea_tally_results_free(alea_tally_results_t* results) {
    if (!results) return;
    for (size_t i = 0; i < results->count; ++i) {
        tally_accumulator_t* tally = &results->tallies[i];
        free((void*)tally->view.bin_ids);
        free((void*)tally->view.energy_edges);
        free(tally->sum);
        free(tally->sum_squared);
        free(tally->values);
        free(tally->used);
        free(tally->touched);
    }
    free(results->tallies);
    free(results);
}

alea_tally_results_t* alea_tally_results_create(
    const alea_tally_plan_t* plan) {
    if (!plan) return NULL;
    alea_tally_results_t* results = calloc(1, sizeof(*results));
    if (!results) return NULL;
    results->count = plan->count;
    results->tallies = calloc(plan->count ? plan->count : 1,
                              sizeof(*results->tallies));
    if (!results->tallies) { free(results); return NULL; }
    for (size_t i = 0; i < plan->count; ++i) {
        const tally_definition_t* definition = &plan->definitions[i];
        tally_accumulator_t* tally = &results->tallies[i];
        size_t bins = definition->bin_count;
        size_t allocation = bins ? bins : 1;
        tally->view.score = definition->spec.score;
        tally->view.domain = definition->spec.domain;
        tally->view.bin_count = bins;
        tally->view.spatial_bin_count = definition->spatial_bin_count;
        tally->view.energy_group_count = definition->spec.energy_group_count;
        tally->view.particle_mask = definition->spec.particle_mask;
        tally->view.energy_min = definition->spec.energy_min;
        tally->view.energy_max = definition->spec.energy_max;
        tally->view.time_min = definition->spec.time_min;
        tally->view.time_max = definition->spec.time_max;
        tally->view.material_id = definition->spec.material_id;
        tally->view.reaction_mt = definition->spec.reaction_mt;
        tally->view.nuclide_zaid = definition->spec.nuclide_zaid;
        memcpy(tally->view.lower, definition->spec.lower, sizeof(tally->view.lower));
        memcpy(tally->view.upper, definition->spec.upper, sizeof(tally->view.upper));
        memcpy(tally->view.dimensions, definition->spec.dimensions,
               sizeof(tally->view.dimensions));
        tally->sum = calloc(allocation, sizeof(double));
        tally->sum_squared = calloc(allocation, sizeof(double));
        tally->view.sum = tally->sum;
        tally->view.sum_squared = tally->sum_squared;
        tally->values = calloc(allocation, sizeof(double));
        tally->used = calloc(allocation, sizeof(uint8_t));
        tally->touched = malloc(allocation * sizeof(size_t));
        if (definition->ids) {
            size_t id_allocation = definition->spatial_bin_count ?
                definition->spatial_bin_count : 1;
            int* ids = malloc(id_allocation * sizeof(int));
            if (ids) memcpy(ids, definition->ids,
                            definition->spatial_bin_count * sizeof(int));
            tally->view.bin_ids = ids;
        }
        if (definition->spec.energy_group_count) {
            size_t edge_count = definition->spec.energy_group_count + 1;
            double* edges = malloc(edge_count * sizeof(double));
            if (edges) memcpy(edges, definition->spec.energy_edges,
                              edge_count * sizeof(double));
            tally->view.energy_edges = edges;
        }
        if (!tally->sum || !tally->sum_squared || !tally->values ||
            !tally->used || !tally->touched ||
            (definition->ids && !tally->view.bin_ids) ||
            (definition->spec.energy_group_count &&
             !tally->view.energy_edges)) {
            alea_tally_results_free(results);
            return NULL;
        }
    }
    return results;
}

size_t alea_tally_results_count(const alea_tally_results_t* results) {
    return results ? results->count : 0;
}

alea_error_t alea_tally_results_view(const alea_tally_results_t* results,
    size_t index, alea_tally_view_t* view) {
    if (!results || !view) return ALEA_ERR_NULL_ARG;
    if (index >= results->count) return ALEA_ERR_INVALID_ARG;
    *view = results->tallies[index].view;
    return ALEA_OK;
}

static int filter_accepts(const alea_tally_spec_t* spec,
    const alea_nav_location_t* location, double energy, double time,
    uint32_t particle_mask, int reaction_mt, int nuclide_zaid) {
    if (spec->particle_mask && !(spec->particle_mask & particle_mask)) return 0;
    if ((spec->energy_min != 0.0 || spec->energy_max != 0.0) &&
        (energy < spec->energy_min || energy >= spec->energy_max)) return 0;
    if ((spec->time_min != 0.0 || spec->time_max != 0.0) &&
        (time < spec->time_min || time >= spec->time_max)) return 0;
    if (spec->material_id > 0 && location->material_id != spec->material_id)
        return 0;
    if (spec->material_id == -1 && location->material_id != 0) return 0;
    if ((spec->score == ALEA_TALLY_REACTION_EVENT ||
         spec->score == ALEA_TALLY_LOCAL_DEPOSITION) &&
        spec->reaction_mt &&
        spec->reaction_mt != reaction_mt) return 0;
    if (spec->score != ALEA_TALLY_REACTION_RATE &&
        spec->score != ALEA_TALLY_HEATING && spec->nuclide_zaid &&
        spec->nuclide_zaid != nuclide_zaid) return 0;
    return 1;
}

static alea_error_t add_score(tally_accumulator_t* tally,
    size_t bin, double score) {
    if (!isfinite(score) || bin >= tally->view.bin_count)
        return ALEA_ERR_INVALID_STATE;
    if (score == 0.0) return ALEA_OK;
    if (!tally->used[bin]) {
        tally->used[bin] = 1;
        tally->touched[tally->touched_count++] = bin;
    }
    tally->values[bin] += score;
    return isfinite(tally->values[bin]) ? ALEA_OK : ALEA_ERR_OVERFLOW;
}

static size_t energy_bin(const tally_definition_t* definition,
                         double energy) {
    size_t groups = definition->spec.energy_group_count;
    if (!groups) return 0;
    const double* edges = definition->spec.energy_edges;
    if (energy < edges[0] || energy >= edges[groups]) return SIZE_MAX;
    size_t lo = 0, hi = groups;
    while (lo < hi) {
        size_t mid = lo + (hi - lo) / 2;
        if (energy < edges[mid + 1]) hi = mid;
        else lo = mid + 1;
    }
    return lo;
}

static alea_error_t score_mesh(tally_accumulator_t* tally,
    const tally_definition_t* definition, const double start[3],
    const double direction[3], double length, double weight,
    size_t energy_offset) {
    const alea_tally_spec_t* spec = &definition->spec;
    double begin = 0.0, end = length;
    for (int axis = 0; axis < 3; ++axis) {
        double d = direction[axis], x = start[axis];
        if (d == 0.0) {
            if (x < spec->lower[axis] || x >= spec->upper[axis])
                return ALEA_OK;
        } else {
            double a = (spec->lower[axis] - x) / d;
            double b = (spec->upper[axis] - x) / d;
            if (a > b) { double tmp = a; a = b; b = tmp; }
            begin = fmax(begin, a);
            end = fmin(end, b);
        }
    }
    if (!(end > begin)) return ALEA_OK;
    double t = begin;
    while (t < end) {
        size_t bin[3];
        double next = end;
        for (int axis = 0; axis < 3; ++axis) {
            double x = fma(direction[axis], t, start[axis]);
            double probe = direction[axis] == 0.0 ? x :
                nextafter(x, direction[axis] > 0.0 ? INFINITY : -INFINITY);
            double scaled = (probe - spec->lower[axis]) / definition->width[axis];
            double index = floor(scaled);
            if (!isfinite(index)) return ALEA_ERR_INVALID_STATE;
            if (index < 0.0) index = 0.0;
            if (index >= spec->dimensions[axis])
                index = spec->dimensions[axis] - 1;
            bin[axis] = (size_t)index;
            if (direction[axis] != 0.0) {
                double edge = spec->lower[axis] +
                    (bin[axis] + (direction[axis] > 0.0 ? 1 : 0)) *
                    definition->width[axis];
                double crossing = (edge - start[axis]) / direction[axis];
                if (crossing > t && crossing < next) next = crossing;
            }
        }
        if (!(next > t)) return ALEA_ERR_INVALID_STATE;
        size_t flat = (bin[2] * spec->dimensions[1] + bin[1]) *
            spec->dimensions[0] + bin[0];
        alea_error_t err = add_score(tally, energy_offset + flat,
                                     weight * (next - t));
        if (err != ALEA_OK) return err;
        t = next;
    }
    return ALEA_OK;
}

static alea_error_t score_location(tally_accumulator_t* tally,
    const tally_definition_t* definition, const alea_nav_location_t* location,
    const double position[3], const double direction[3], double distance,
    double weight, size_t energy_offset) {
    if (definition->spec.domain == ALEA_TALLY_CARTESIAN_MESH) {
        if (distance > 0.0)
            return score_mesh(tally, definition, position, direction,
                              distance, weight, energy_offset);
        size_t bin[3];
        for (int axis = 0; axis < 3; ++axis) {
            double x = position[axis];
            if (x < definition->spec.lower[axis] ||
                x >= definition->spec.upper[axis]) return ALEA_OK;
            bin[axis] = (size_t)floor((x - definition->spec.lower[axis]) /
                                      definition->width[axis]);
            if (bin[axis] >= definition->spec.dimensions[axis])
                bin[axis] = definition->spec.dimensions[axis] - 1;
        }
        size_t flat = (bin[2] * definition->spec.dimensions[1] + bin[1]) *
            definition->spec.dimensions[0] + bin[0];
        return add_score(tally, energy_offset + flat, weight);
    }
    if (location->cell_index < 0) return ALEA_OK;
    size_t cell = (size_t)location->cell_index;
    if (cell >= definition->cell_count) return ALEA_ERR_INVALID_STATE;
    size_t bin = definition->cell_to_bin[cell];
    return add_score(tally, energy_offset + bin,
                     weight * (distance > 0.0 ? distance : 1.0));
}

alea_error_t alea_tally_record_track(alea_tally_results_t* results,
    const alea_tally_plan_t* plan, const alea_nav_location_t* location,
    const double start[3], const double direction[3], double distance,
    double weight, double energy, double time, double speed,
    uint32_t particle_mask, const alea_nuc_evaluation_t* evaluation) {
    if (!results || !plan || !location || !start || !direction)
        return ALEA_ERR_NULL_ARG;
    if (!isfinite(distance) || distance < 0.0 || !isfinite(weight) ||
        !isfinite(energy) || !isfinite(time) || !isfinite(speed) ||
        speed <= 0.0) return ALEA_ERR_INVALID_ARG;
    if (distance == 0.0) return ALEA_OK;
    for (size_t i = 0; i < plan->count; ++i) {
        const tally_definition_t* definition = &plan->definitions[i];
        if (definition->spec.score != ALEA_TALLY_TRACK_LENGTH &&
            definition->spec.score != ALEA_TALLY_REACTION_RATE &&
            definition->spec.score != ALEA_TALLY_HEATING) continue;
        double begin = 0.0, end = distance;
        if (definition->spec.time_min != 0.0 ||
            definition->spec.time_max != 0.0) {
            begin = fmax(begin, (definition->spec.time_min - time) * speed);
            end = fmin(end, (definition->spec.time_max - time) * speed);
        }
        double filter_time =
            definition->spec.time_min != 0.0 ||
            definition->spec.time_max != 0.0 ?
                definition->spec.time_min : time;
        if (!(end > begin) ||
            !filter_accepts(&definition->spec, location, energy,
                            filter_time, particle_mask, 0, 0))
            continue;
        size_t group = energy_bin(definition, energy);
        if (group == SIZE_MAX) continue;
        double response = 1.0;
        if (definition->spec.score == ALEA_TALLY_REACTION_RATE ||
            definition->spec.score == ALEA_TALLY_HEATING) {
            if (particle_mask != ALEA_TALLY_NEUTRON || !evaluation)
                continue; /* neutron-only response; void has no rate */
            alea_error_t response_err =
                definition->spec.score == ALEA_TALLY_HEATING ?
                alea_nuc_evaluation_macro_heating(evaluation,
                    definition->spec.nuclide_zaid, &response) :
                alea_nuc_evaluation_macro_reaction_rate(evaluation,
                    definition->spec.reaction_mt,
                    definition->spec.nuclide_zaid, &response);
            if (response_err != ALEA_OK) return response_err;
        }
        double clipped_start[3] = {
            fma(direction[0], begin, start[0]),
            fma(direction[1], begin, start[1]),
            fma(direction[2], begin, start[2])
        };
        alea_error_t err = score_location(&results->tallies[i], definition,
            location, clipped_start, direction, end - begin,
            weight * response, group * definition->spatial_bin_count);
        if (err != ALEA_OK) return err;
    }
    return ALEA_OK;
}

alea_error_t alea_tally_record_collision(alea_tally_results_t* results,
    const alea_tally_plan_t* plan, const alea_nav_location_t* location,
    const double position[3], double weight, double energy, double time,
    uint32_t particle_mask, int reaction_mt, int nuclide_zaid) {
    if (!results || !plan || !location || !position) return ALEA_ERR_NULL_ARG;
    if (!isfinite(weight) || !isfinite(energy) || !isfinite(time))
        return ALEA_ERR_INVALID_ARG;
    for (size_t i = 0; i < plan->count; ++i) {
        const tally_definition_t* definition = &plan->definitions[i];
        if (definition->spec.score == ALEA_TALLY_TRACK_LENGTH ||
            definition->spec.score == ALEA_TALLY_REACTION_RATE ||
            definition->spec.score == ALEA_TALLY_HEATING ||
            definition->spec.score == ALEA_TALLY_LOCAL_DEPOSITION ||
            !filter_accepts(&definition->spec, location, energy, time,
                            particle_mask, reaction_mt,
                            nuclide_zaid)) continue;
        size_t group = energy_bin(definition, energy);
        if (group == SIZE_MAX) continue;
        alea_error_t err = score_location(&results->tallies[i], definition,
            location, position, NULL, 0.0, weight,
            group * definition->spatial_bin_count);
        if (err != ALEA_OK) return err;
    }
    return ALEA_OK;
}

alea_error_t alea_tally_record_deposition(alea_tally_results_t* results,
    const alea_tally_plan_t* plan, const alea_nav_location_t* location,
    const double position[3], double weight, double energy, double time,
    uint32_t particle_mask, int reaction_mt, int nuclide_zaid,
    double deposited_energy) {
    if (!results || !plan || !location || !position) return ALEA_ERR_NULL_ARG;
    if (!isfinite(weight) || !isfinite(energy) || !isfinite(time) ||
        !isfinite(deposited_energy)) return ALEA_ERR_INVALID_ARG;
    for (size_t i = 0; i < plan->count; ++i) {
        const tally_definition_t* definition = &plan->definitions[i];
        if (definition->spec.score != ALEA_TALLY_LOCAL_DEPOSITION ||
            !filter_accepts(&definition->spec, location, energy, time,
                            particle_mask, reaction_mt, nuclide_zaid))
            continue;
        size_t group = energy_bin(definition, energy);
        if (group == SIZE_MAX) continue;
        alea_error_t err = score_location(&results->tallies[i], definition,
            location, position, NULL, 0.0, weight * deposited_energy,
            group * definition->spatial_bin_count);
        if (err != ALEA_OK) return err;
    }
    return ALEA_OK;
}

alea_error_t alea_tally_commit_history(alea_tally_results_t* results) {
    if (!results) return ALEA_ERR_NULL_ARG;
    for (size_t i = 0; i < results->count; ++i) {
        tally_accumulator_t* tally = &results->tallies[i];
        for (size_t j = 0; j < tally->touched_count; ++j) {
            size_t bin = tally->touched[j];
            double value = tally->values[bin];
            double sum = tally->sum[bin] + value;
            double squared = tally->sum_squared[bin] + value * value;
            if (!isfinite(sum) || !isfinite(squared)) return ALEA_ERR_OVERFLOW;
            tally->sum[bin] = sum;
            tally->sum_squared[bin] = squared;
            tally->values[bin] = 0.0;
            tally->used[bin] = 0;
        }
        tally->touched_count = 0;
        tally->view.histories++;
    }
    return ALEA_OK;
}

void alea_tally_results_finalize(alea_tally_results_t* results) {
    if (!results) return;
    for (size_t i = 0; i < results->count; ++i) {
        tally_accumulator_t* tally = &results->tallies[i];
        free(tally->values); tally->values = NULL;
        free(tally->used); tally->used = NULL;
        free(tally->touched); tally->touched = NULL;
    }
}
