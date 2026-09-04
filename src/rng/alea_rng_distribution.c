// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

#include "rng/alea_rng_distribution.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

#define ALEA_RNG_EVENT_DRAW_CAPACITY (UINT64_C(1) << 32)
#define ALEA_RNG_DEFAULT_REJECTION_LIMIT UINT32_C(128)
#define ALEA_TWO_PI 6.283185307179586476925286766559

struct alea_rng_discrete_table {
    size_t count;
    double* cdf;
    double* alias_probability;
    uint32_t* alias;
};

struct alea_rng_tabular_table {
    size_t point_count;
    alea_rng_tabular_interpolation_t interpolation;
    double total_area;
    double* x;
    double* pdf;
    double* cumulative_area;
};

static int rng_algorithm_valid(alea_rng_algorithm_t algorithm) {
    /* The compatibility LCG intentionally receives no new production APIs. */
    return algorithm == ALEA_RNG_PHILOX4X32_10;
}

static void rng_invalid(const char* message) {
    alea_set_error_detail(ALEA_ERR_INVALID_ARG, message);
}

static void rng_overflow(const char* message) {
    alea_set_error_detail(ALEA_ERR_OVERFLOW, message);
}

static int rng_allocate_array(size_t count, size_t element_size, void** output) {
    if (!output || count == 0 || element_size > SIZE_MAX / count) {
        rng_overflow("RNG distribution table size overflows");
        return -1;
    }
    *output = calloc(count, element_size);
    if (!*output) {
        alea_set_error_detail(ALEA_ERR_OUT_OF_MEMORY,
                              "failed to allocate RNG distribution table");
        return -1;
    }
    return 0;
}

uint64_t alea_rng_transport_entity_id(uint32_t history_id,
                                      uint32_t particle_ordinal) {
    return ((uint64_t)history_id << 32) | (uint64_t)particle_ordinal;
}

int alea_rng_event_init(alea_rng_event_t* event,
                        alea_rng_algorithm_t algorithm,
                        uint64_t seed,
                        uint32_t domain,
                        uint64_t entity_id,
                        uint32_t event_index) {
    if (!event) {
        alea_set_error_detail(ALEA_ERR_NULL_ARG,
                              "RNG event descriptor is required");
        return -1;
    }
    if (!rng_algorithm_valid(algorithm) || domain == 0) {
        rng_invalid("RNG event algorithm or domain is invalid");
        return -1;
    }
    memset(event, 0, sizeof(*event));
    event->algorithm = algorithm;
    event->seed = seed;
    event->domain = domain;
    event->entity_id = entity_id;
    event->event_index = event_index;
    return 0;
}

uint64_t alea_rng_event_position(const alea_rng_event_t* event) {
    return event ? event->local_draw : 0;
}

int alea_rng_event_seek(alea_rng_event_t* event, uint64_t local_draw) {
    if (!event || local_draw > ALEA_RNG_EVENT_DRAW_CAPACITY) {
        rng_invalid("RNG event draw position is outside its 32-bit range");
        return -1;
    }
    event->local_draw = local_draw;
    event->cached_block_valid = 0;
    return 0;
}

int alea_rng_event_next_u32(alea_rng_event_t* event, uint32_t* output) {
    if (!event || !output) {
        alea_set_error_detail(ALEA_ERR_NULL_ARG,
                              "RNG event and output are required");
        return -1;
    }
    if (event->local_draw >= ALEA_RNG_EVENT_DRAW_CAPACITY) {
        rng_overflow("RNG event exhausted its 32-bit local draw space");
        return -1;
    }
    const uint64_t address =
        ((uint64_t)event->event_index << 32) | event->local_draw;
    const uint64_t block_index = address >> 2;
    if (!event->cached_block_valid ||
        event->cached_block_index != block_index) {
        if (alea_rng_block_at(event->algorithm, event->seed, event->domain,
                              event->entity_id, block_index,
                              event->cached_block) != 0) {
            rng_invalid("RNG event contains an unsupported algorithm");
            return -1;
        }
        event->cached_block_index = block_index;
        event->cached_block_valid = 1;
    }
    *output = event->cached_block[address & UINT64_C(3)];
    event->local_draw++;
    return 0;
}

int alea_rng_event_next_uniform(alea_rng_event_t* event, double* output) {
    if (!event || !output) {
        alea_set_error_detail(ALEA_ERR_NULL_ARG,
                              "RNG event and uniform output are required");
        return -1;
    }
    alea_rng_event_t saved = *event;
    uint32_t high, low;
    if (alea_rng_event_next_u32(event, &high) != 0 ||
        alea_rng_event_next_u32(event, &low) != 0) {
        *event = saved;
        return -1;
    }
    *output = alea_rng_uniform53_from_u32(high, low);
    return 0;
}

int alea_rng_event_next_uniform_open(alea_rng_event_t* event, double* output) {
    if (!event || !output) {
        alea_set_error_detail(ALEA_ERR_NULL_ARG,
                              "RNG event and open-uniform output are required");
        return -1;
    }
    alea_rng_event_t saved = *event;
    uint32_t high, low;
    if (alea_rng_event_next_u32(event, &high) != 0 ||
        alea_rng_event_next_u32(event, &low) != 0) {
        *event = saved;
        return -1;
    }
    *output = alea_rng_uniform53_open_from_u32(high, low);
    return 0;
}

int alea_rng_sample_u32_bounded(alea_rng_event_t* event,
                                uint32_t bound,
                                uint32_t maximum_attempts,
                                uint32_t* output) {
    if (!event || !output || bound == 0) {
        rng_invalid("bounded integer sampling requires an event, output, and positive bound");
        return -1;
    }
    if (maximum_attempts == 0)
        maximum_attempts = ALEA_RNG_DEFAULT_REJECTION_LIMIT;
    alea_rng_event_t saved = *event;
    const uint32_t threshold = (uint32_t)(-bound) % bound;
    for (uint32_t attempt = 0; attempt < maximum_attempts; attempt++) {
        uint32_t bits;
        if (alea_rng_event_next_u32(event, &bits) != 0) {
            *event = saved;
            return -1;
        }
        const uint64_t product = (uint64_t)bits * (uint64_t)bound;
        if ((uint32_t)product >= threshold) {
            *output = (uint32_t)(product >> 32);
            return 0;
        }
    }
    *event = saved;
    alea_set_error_detail(ALEA_ERR_INVALID_STATE,
                          "bounded integer rejection limit exhausted");
    return -1;
}

int alea_rng_sample_uniform_range(alea_rng_event_t* event,
                                  double lower,
                                  double upper,
                                  double* output) {
    if (!event || !output || !isfinite(lower) || !isfinite(upper) ||
        !(lower < upper) || !isfinite(upper - lower)) {
        rng_invalid("uniform interval requires finite increasing bounds");
        return -1;
    }
    double u;
    if (alea_rng_event_next_uniform(event, &u) != 0) return -1;
    *output = lower + (upper - lower) * u;
    return 0;
}

int alea_rng_sample_exponential(alea_rng_event_t* event,
                                double rate,
                                double* output) {
    if (!event || !output || !isfinite(rate) || !(rate > 0.0)) {
        rng_invalid("exponential sampling requires a finite positive rate");
        return -1;
    }
    double u;
    if (alea_rng_event_next_uniform_open(event, &u) != 0) return -1;
    *output = -log(u) / rate;
    return 0;
}

int alea_rng_sample_isotropic_direction(alea_rng_event_t* event,
                                        double direction[3]) {
    if (!event || !direction) {
        alea_set_error_detail(ALEA_ERR_NULL_ARG,
                              "isotropic sampling requires an event and output");
        return -1;
    }
    alea_rng_event_t saved = *event;
    double azimuth_u, cosine_u;
    if (alea_rng_event_next_uniform(event, &azimuth_u) != 0 ||
        alea_rng_event_next_uniform(event, &cosine_u) != 0) {
        *event = saved;
        return -1;
    }
    const double phi = ALEA_TWO_PI * azimuth_u;
    const double z = 2.0 * cosine_u - 1.0;
    const double radial = sqrt(fmax(0.0, 1.0 - z * z));
    direction[0] = radial * cos(phi);
    direction[1] = radial * sin(phi);
    direction[2] = z;
    return 0;
}

int alea_rng_sample_cosine_direction(alea_rng_event_t* event,
                                     double direction[3]) {
    if (!event || !direction) {
        alea_set_error_detail(ALEA_ERR_NULL_ARG,
                              "cosine sampling requires an event and output");
        return -1;
    }
    alea_rng_event_t saved = *event;
    double radial_u, azimuth_u;
    if (alea_rng_event_next_uniform(event, &radial_u) != 0 ||
        alea_rng_event_next_uniform(event, &azimuth_u) != 0) {
        *event = saved;
        return -1;
    }
    const double radial = sqrt(radial_u);
    const double phi = ALEA_TWO_PI * azimuth_u;
    direction[0] = radial * cos(phi);
    direction[1] = radial * sin(phi);
    direction[2] = sqrt(fmax(0.0, 1.0 - radial_u));
    return 0;
}

alea_rng_discrete_table_t* alea_rng_discrete_table_create(
        const double* weights, size_t count) {
    if (!weights || count == 0 || count > UINT32_MAX) {
        rng_invalid("discrete table requires 1..2^32-1 weights");
        return NULL;
    }
    alea_rng_discrete_table_t* table = calloc(1, sizeof(*table));
    double* scaled = NULL;
    uint32_t* small = NULL;
    uint32_t* large = NULL;
    if (!table) {
        alea_set_error_detail(ALEA_ERR_OUT_OF_MEMORY,
                              "failed to allocate discrete distribution table");
        return NULL;
    }
    if (
        rng_allocate_array(count, sizeof(*table->cdf),
                           (void**)&table->cdf) != 0 ||
        rng_allocate_array(count, sizeof(*table->alias_probability),
                           (void**)&table->alias_probability) != 0 ||
        rng_allocate_array(count, sizeof(*table->alias),
                           (void**)&table->alias) != 0 ||
        rng_allocate_array(count, sizeof(*scaled), (void**)&scaled) != 0 ||
        rng_allocate_array(count, sizeof(*small), (void**)&small) != 0 ||
        rng_allocate_array(count, sizeof(*large), (void**)&large) != 0) goto fail;

    double sum = 0.0, correction = 0.0;
    for (size_t i = 0; i < count; i++) {
        if (!isfinite(weights[i]) || weights[i] < 0.0) {
            rng_invalid("discrete weights must be finite and nonnegative");
            goto fail;
        }
        const double adjusted = weights[i] - correction;
        const double next = sum + adjusted;
        correction = (next - sum) - adjusted;
        sum = next;
    }
    if (!(sum > 0.0) || !isfinite(sum)) {
        rng_invalid("discrete weights must have a finite positive sum");
        goto fail;
    }

    double cumulative = 0.0;
    size_t small_count = 0, large_count = 0;
    for (size_t i = 0; i < count; i++) {
        cumulative += weights[i] / sum;
        table->cdf[i] = cumulative;
        scaled[i] = (weights[i] / sum) * (double)count;
        if (scaled[i] < 1.0) small[small_count++] = (uint32_t)i;
        else large[large_count++] = (uint32_t)i;
    }
    table->cdf[count - 1] = 1.0;

    while (small_count && large_count) {
        const uint32_t s = small[--small_count];
        const uint32_t l = large[--large_count];
        table->alias_probability[s] = scaled[s];
        if (table->alias_probability[s] < 0.0)
            table->alias_probability[s] = 0.0;
        if (table->alias_probability[s] > 1.0)
            table->alias_probability[s] = 1.0;
        table->alias[s] = l;
        scaled[l] = (scaled[l] + scaled[s]) - 1.0;
        if (scaled[l] < 1.0) small[small_count++] = l;
        else large[large_count++] = l;
    }
    while (large_count) {
        const uint32_t i = large[--large_count];
        table->alias_probability[i] = 1.0;
        table->alias[i] = i;
    }
    while (small_count) {
        const uint32_t i = small[--small_count];
        table->alias_probability[i] = 1.0;
        table->alias[i] = i;
    }
    table->count = count;
    free(scaled); free(small); free(large);
    return table;

fail:
    free(scaled); free(small); free(large);
    alea_rng_discrete_table_destroy(table);
    return NULL;
}

void alea_rng_discrete_table_destroy(alea_rng_discrete_table_t* table) {
    if (!table) return;
    free(table->cdf);
    free(table->alias_probability);
    free(table->alias);
    free(table);
}

size_t alea_rng_discrete_table_count(const alea_rng_discrete_table_t* table) {
    return table ? table->count : 0;
}

int alea_rng_sample_discrete_cdf(alea_rng_event_t* event,
                                 const alea_rng_discrete_table_t* table,
                                 size_t* output_index) {
    if (!event || !table || !output_index || table->count == 0) {
        rng_invalid("CDF sampling requires an event, table, and output");
        return -1;
    }
    double u;
    if (alea_rng_event_next_uniform(event, &u) != 0) return -1;
    size_t low = 0, high = table->count;
    while (low < high) {
        const size_t middle = low + (high - low) / 2;
        if (u < table->cdf[middle]) high = middle;
        else low = middle + 1;
    }
    *output_index = low < table->count ? low : table->count - 1;
    return 0;
}

int alea_rng_sample_discrete_alias(alea_rng_event_t* event,
                                   const alea_rng_discrete_table_t* table,
                                   size_t* output_index) {
    if (!event || !table || !output_index || table->count == 0) {
        rng_invalid("alias sampling requires an event, table, and output");
        return -1;
    }
    alea_rng_event_t saved = *event;
    uint32_t column, bits;
    if (alea_rng_sample_u32_bounded(event, (uint32_t)table->count, 0,
                                    &column) != 0 ||
        alea_rng_event_next_u32(event, &bits) != 0) {
        *event = saved;
        return -1;
    }
    const double u = alea_rng_uniform32_from_u32(bits);
    *output_index = u < table->alias_probability[column]
        ? (size_t)column : (size_t)table->alias[column];
    return 0;
}

alea_rng_tabular_table_t* alea_rng_tabular_table_create(
        const double* x, const double* pdf, size_t point_count,
        alea_rng_tabular_interpolation_t interpolation) {
    if (!x || !pdf || point_count < 2 ||
        (interpolation != ALEA_RNG_TABULAR_HISTOGRAM &&
         interpolation != ALEA_RNG_TABULAR_LIN_LIN)) {
        rng_invalid("tabular distribution requires valid points and interpolation");
        return NULL;
    }
    alea_rng_tabular_table_t* table = calloc(1, sizeof(*table));
    if (!table) {
        alea_set_error_detail(ALEA_ERR_OUT_OF_MEMORY,
                              "failed to allocate tabular distribution table");
        return NULL;
    }
    if (
        rng_allocate_array(point_count, sizeof(*table->x),
                           (void**)&table->x) != 0 ||
        rng_allocate_array(point_count, sizeof(*table->pdf),
                           (void**)&table->pdf) != 0 ||
        rng_allocate_array(point_count, sizeof(*table->cumulative_area),
                           (void**)&table->cumulative_area) != 0) goto fail;
    memcpy(table->x, x, point_count * sizeof(*x));
    memcpy(table->pdf, pdf, point_count * sizeof(*pdf));
    table->point_count = point_count;
    table->interpolation = interpolation;

    double sum = 0.0, correction = 0.0;
    for (size_t i = 0; i < point_count; i++) {
        if (!isfinite(x[i]) || !isfinite(pdf[i]) || pdf[i] < 0.0 ||
            (i > 0 && !(x[i] > x[i - 1]))) {
            rng_invalid("tabular coordinates must increase and PDF values must be nonnegative");
            goto fail;
        }
        if (i + 1 == point_count) continue;
        const double width = x[i + 1] - x[i];
        const double height = interpolation == ALEA_RNG_TABULAR_HISTOGRAM
            ? pdf[i] : 0.5 * (pdf[i] + pdf[i + 1]);
        const double area = width * height;
        if (!isfinite(area) || area < 0.0) {
            rng_invalid("tabular distribution area is invalid");
            goto fail;
        }
        const double adjusted = area - correction;
        const double next = sum + adjusted;
        correction = (next - sum) - adjusted;
        sum = next;
        table->cumulative_area[i + 1] = sum;
    }
    if (!(sum > 0.0) || !isfinite(sum)) {
        rng_invalid("tabular distribution must have finite positive area");
        goto fail;
    }
    table->total_area = sum;
    return table;

fail:
    alea_rng_tabular_table_destroy(table);
    return NULL;
}

void alea_rng_tabular_table_destroy(alea_rng_tabular_table_t* table) {
    if (!table) return;
    free(table->x);
    free(table->pdf);
    free(table->cumulative_area);
    free(table);
}

int alea_rng_sample_tabular(alea_rng_event_t* event,
                            const alea_rng_tabular_table_t* table,
                            double* output) {
    if (!event || !table || !output || table->point_count < 2) {
        rng_invalid("tabular sampling requires an event, table, and output");
        return -1;
    }
    double u;
    if (alea_rng_event_next_uniform(event, &u) != 0) return -1;
    double target = u * table->total_area;
    if (!(target < table->total_area))
        target = nextafter(table->total_area, 0.0);
    size_t low = 0, high = table->point_count - 1;
    while (low < high) {
        const size_t middle = low + (high - low) / 2;
        if (target < table->cumulative_area[middle + 1]) high = middle;
        else low = middle + 1;
    }
    size_t segment = low;
    if (segment + 1 >= table->point_count)
        segment = table->point_count - 2;
    const double residual = target - table->cumulative_area[segment];
    const double width = table->x[segment + 1] - table->x[segment];
    double offset;
    if (residual <= 0.0) {
        offset = 0.0;
    } else if (table->interpolation == ALEA_RNG_TABULAR_HISTOGRAM ||
        table->pdf[segment] == table->pdf[segment + 1]) {
        offset = residual / table->pdf[segment];
    } else {
        const double slope =
            (table->pdf[segment + 1] - table->pdf[segment]) / width;
        if (table->pdf[segment] == 0.0 && slope > 0.0) {
            offset = sqrt(2.0 * residual / slope);
        } else {
            const double discriminant = fmax(
                0.0, table->pdf[segment] * table->pdf[segment] +
                     2.0 * slope * residual);
            const double denominator =
                table->pdf[segment] + sqrt(discriminant);
            offset = denominator > 0.0
                ? 2.0 * residual / denominator : 0.0;
        }
    }
    if (offset < 0.0) offset = 0.0;
    if (offset > width) offset = width;
    *output = table->x[segment] + offset;
    return 0;
}
