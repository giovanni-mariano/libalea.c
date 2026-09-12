// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

#include "nuclear_internal.h"

#include <stdio.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static void thermal_clear(alea_nuc_thermal_t* thermal) {
    if (!thermal) return;
    free(thermal->inelastic_energy);
    free(thermal->inelastic_xs);
    free(thermal->inelastic_outgoing_offset);
    free(thermal->inelastic_energy_out);
    free(thermal->inelastic_pdf);
    free(thermal->inelastic_cdf);
    free(thermal->inelastic_mu);
    free(thermal->coherent_edge);
    free(thermal->coherent_factor);
    free(thermal->incoherent_energy);
    free(thermal->incoherent_xs);
    free(thermal->incoherent_mu);
}

void alea_nuc_thermal_free(alea_nuc_thermal_t* thermal) {
    if (!thermal) return;
    thermal_clear(thermal);
    free(thermal);
}

static int finite_ascending(const double* values, int n, int positive) {
    if (!values || n < 1) return 0;
    for (int i = 0; i < n; i++)
        if (!isfinite(values[i]) || (positive && values[i] <= 0.0) ||
            (i > 0 && values[i] <= values[i - 1]))
            return 0;
    return 1;
}

static int finite_nonnegative(const double* values, int n) {
    if (!values || n < 1) return 0;
    for (int i = 0; i < n; i++)
        if (!isfinite(values[i]) || values[i] < 0.0) return 0;
    return 1;
}

int alea_nuc_thermal_validate_internal(const alea_nuc_thermal_t* thermal) {
    if (!thermal || !isfinite(thermal->awr) || thermal->awr <= 0.0 ||
        !isfinite(thermal->temperature) || thermal->temperature <= 0.0 ||
        thermal->n_applicable_zaids < 1 ||
        thermal->n_applicable_zaids > 16 ||
        thermal->n_inelastic_energies < 2 ||
        thermal->n_inelastic_cosines < 1 ||
        !finite_ascending(thermal->inelastic_energy,
                          thermal->n_inelastic_energies, 1) ||
        !finite_nonnegative(thermal->inelastic_xs,
                            thermal->n_inelastic_energies) ||
        !thermal->inelastic_energy_out || !thermal->inelastic_mu)
        return 0;
    size_t records;
    if (thermal->inelastic_continuous) {
        if (thermal->n_inelastic_outgoing != 0 ||
            thermal->n_inelastic_outgoing_total < 2 ||
            !thermal->inelastic_outgoing_offset ||
            !thermal->inelastic_pdf || !thermal->inelastic_cdf ||
            thermal->inelastic_outgoing_offset[0] != 0 ||
            thermal->inelastic_outgoing_offset[
                thermal->n_inelastic_energies] !=
                    thermal->n_inelastic_outgoing_total)
            return 0;
        records = (size_t)thermal->n_inelastic_outgoing_total;
        for (int i = 0; i < thermal->n_inelastic_energies; i++) {
            int begin = thermal->inelastic_outgoing_offset[i];
            int end = thermal->inelastic_outgoing_offset[i + 1];
            if (begin < 0 || end - begin < 2 ||
                end > thermal->n_inelastic_outgoing_total)
                return 0;
            for (int j = begin; j < end; j++) {
                if (!isfinite(thermal->inelastic_energy_out[j]) ||
                    thermal->inelastic_energy_out[j] < 0.0 ||
                    !isfinite(thermal->inelastic_pdf[j]) ||
                    thermal->inelastic_pdf[j] < 0.0 ||
                    !isfinite(thermal->inelastic_cdf[j]) ||
                    thermal->inelastic_cdf[j] < 0.0 ||
                    (j > begin &&
                     (thermal->inelastic_energy_out[j] <=
                          thermal->inelastic_energy_out[j - 1] ||
                      thermal->inelastic_cdf[j] <
                          thermal->inelastic_cdf[j - 1])))
                    return 0;
            }
            if (thermal->inelastic_cdf[begin] != 0.0 ||
                fabs(thermal->inelastic_cdf[end - 1] - 1.0) > 1e-6)
                return 0;
        }
    } else {
        if (thermal->n_inelastic_outgoing < 1 ||
            (size_t)thermal->n_inelastic_energies >
                SIZE_MAX / (size_t)thermal->n_inelastic_outgoing)
            return 0;
        records = (size_t)thermal->n_inelastic_energies *
                  (size_t)thermal->n_inelastic_outgoing;
    }
    if (records > SIZE_MAX / (size_t)thermal->n_inelastic_cosines) return 0;
    for (int i = 0; i < thermal->n_applicable_zaids; i++)
        if (thermal->applicable_zaids[i] <= 0) return 0;
    for (size_t i = 0; i < records; i++) {
        if (!isfinite(thermal->inelastic_energy_out[i]) ||
            thermal->inelastic_energy_out[i] < 0.0)
            return 0;
        if (!thermal->inelastic_continuous &&
            i % (size_t)thermal->n_inelastic_outgoing != 0 &&
            thermal->inelastic_energy_out[i] <
                thermal->inelastic_energy_out[i - 1])
            return 0;
    }
    size_t nmu = records * (size_t)thermal->n_inelastic_cosines;
    for (size_t i = 0; i < nmu; i++) {
        if (!isfinite(thermal->inelastic_mu[i]) ||
            thermal->inelastic_mu[i] < -1.0 || thermal->inelastic_mu[i] > 1.0)
            return 0;
        if (i % (size_t)thermal->n_inelastic_cosines != 0 &&
            thermal->inelastic_mu[i] < thermal->inelastic_mu[i - 1])
            return 0;
    }

    int coherent = thermal->elastic_mode == ALEA_NUC_THERMAL_ELASTIC_COHERENT ||
                   thermal->elastic_mode == ALEA_NUC_THERMAL_ELASTIC_MIXED;
    int incoherent =
        thermal->elastic_mode == ALEA_NUC_THERMAL_ELASTIC_INCOHERENT ||
        thermal->elastic_mode == ALEA_NUC_THERMAL_ELASTIC_MIXED;
    if (thermal->elastic_mode != ALEA_NUC_THERMAL_ELASTIC_NONE &&
        !coherent && !incoherent)
        return 0;
    if (coherent) {
        if (thermal->n_coherent_edges < 1 ||
            !finite_ascending(thermal->coherent_edge,
                              thermal->n_coherent_edges, 1) ||
            !finite_nonnegative(thermal->coherent_factor,
                                thermal->n_coherent_edges))
            return 0;
        for (int i = 1; i < thermal->n_coherent_edges; i++)
            if (thermal->coherent_factor[i] < thermal->coherent_factor[i - 1])
                return 0;
    }
    if (incoherent) {
        if (thermal->n_incoherent_energies < 2 ||
            thermal->n_incoherent_cosines < 1 ||
            !finite_ascending(thermal->incoherent_energy,
                              thermal->n_incoherent_energies, 1) ||
            !finite_nonnegative(thermal->incoherent_xs,
                                thermal->n_incoherent_energies) ||
            !thermal->incoherent_mu ||
            (size_t)thermal->n_incoherent_energies >
                SIZE_MAX / (size_t)thermal->n_incoherent_cosines)
            return 0;
        size_t count = (size_t)thermal->n_incoherent_energies *
                       (size_t)thermal->n_incoherent_cosines;
        for (size_t i = 0; i < count; i++) {
            if (!isfinite(thermal->incoherent_mu[i]) ||
                thermal->incoherent_mu[i] < -1.0 ||
                thermal->incoherent_mu[i] > 1.0)
                return 0;
            if (i % (size_t)thermal->n_incoherent_cosines != 0 &&
                thermal->incoherent_mu[i] < thermal->incoherent_mu[i - 1])
                return 0;
        }
    }
    return 1;
}

static int compare_double(const void* left, const void* right) {
    double a = *(const double*)left;
    double b = *(const double*)right;
    return (a > b) - (a < b);
}

static alea_error_t decode_continuous_inelastic(
    const alea_nuc_ace_table_t* table, alea_nuc_thermal_t* thermal,
    int itxe, int ne) {
    int dimension = table->nxs[2]; /* NXS(3) */
    if (dimension < 2 || dimension > INT_MAX - 2)
        return ALEA_ERR_PARSE_ERROR;
    int nmu = dimension - 1;
    if (!xss_range_valid(table, itxe, 2 * ne))
        return ALEA_ERR_PARSE_ERROR;
    thermal->inelastic_continuous = true;
    thermal->n_inelastic_cosines = nmu;
    thermal->inelastic_outgoing_offset =
        alea_nuc_calloc((size_t)ne + 1, sizeof(int));
    if (!thermal->inelastic_outgoing_offset) return ALEA_ERR_OUT_OF_MEMORY;

    int total = 0;
    int stride = nmu + 3;
    for (int i = 0; i < ne; i++) {
        int locator = xss_int(table, itxe + i);
        int count = xss_int(table, itxe + ne + i);
        if (locator < 1 || locator >= table->xss_length || count < 2 ||
            count > INT_MAX / stride ||
            !xss_range_valid(table, locator + 1, count * stride))
            return ALEA_ERR_PARSE_ERROR;
        double first_energy = xss(table, locator + 1);
        double first_cdf = xss(table, locator + 3);
        int insert_zero = first_cdf > 0.0;
        if (insert_zero && (!(first_energy > 0.0) || !isfinite(first_energy)))
            return ALEA_ERR_PARSE_ERROR;
        if (count > INT_MAX - total - insert_zero) return ALEA_ERR_OVERFLOW;
        total += count + insert_zero;
        thermal->inelastic_outgoing_offset[i + 1] = total;
    }
    thermal->n_inelastic_outgoing_total = total;
    thermal->inelastic_energy_out = alea_nuc_malloc((size_t)total * sizeof(double));
    thermal->inelastic_pdf = alea_nuc_malloc((size_t)total * sizeof(double));
    thermal->inelastic_cdf = alea_nuc_malloc((size_t)total * sizeof(double));
    if ((size_t)total > SIZE_MAX / (size_t)nmu)
        return ALEA_ERR_OVERFLOW;
    thermal->inelastic_mu =
        alea_nuc_malloc((size_t)total * (size_t)nmu * sizeof(double));
    if (!thermal->inelastic_energy_out || !thermal->inelastic_pdf ||
        !thermal->inelastic_cdf || !thermal->inelastic_mu)
        return ALEA_ERR_OUT_OF_MEMORY;

    for (int i = 0; i < ne; i++) {
        int locator = xss_int(table, itxe + i);
        int raw_count = xss_int(table, itxe + ne + i);
        int output = thermal->inelastic_outgoing_offset[i];
        if (xss(table, locator + 3) > 0.0) {
            thermal->inelastic_energy_out[output] = 0.0;
            thermal->inelastic_pdf[output] = 0.0;
            thermal->inelastic_cdf[output] = 0.0;
            double* mu = &thermal->inelastic_mu[(size_t)output * (size_t)nmu];
            for (int k = 0; k < nmu; k++)
                mu[k] = -1.0 + (2.0 * k + 1.0) / nmu;
            output++;
        }
        for (int j = 0; j < raw_count; j++, output++) {
            int base = locator + 1 + j * stride;
            thermal->inelastic_energy_out[output] = xss(table, base);
            double pdf = xss(table, base + 1);
            if (pdf < 0.0 && pdf >= -1e-12) pdf = 0.0;
            thermal->inelastic_pdf[output] = pdf;
            thermal->inelastic_cdf[output] = xss(table, base + 2);
            double* mu = &thermal->inelastic_mu[
                (size_t)output * (size_t)nmu];
            for (int k = 0; k < nmu; k++)
                mu[k] = xss(table, base + 3 + k);
            qsort(mu, (size_t)nmu, sizeof(double), compare_double);
        }
    }
    return ALEA_OK;
}

static alea_error_t decode_inelastic(const alea_nuc_ace_table_t* table,
                                     alea_nuc_thermal_t* thermal) {
    int itie = table->jxs[0]; /* JXS(1) */
    int itxe = table->jxs[2]; /* JXS(3) */
    int ne = itie > 0 ? xss_int(table, itie) : 0;
    int ifeng = table->nxs[6]; /* NXS(7) */
    if (itie <= 0 || itxe <= 0 || ne < 2 || ne > INT_MAX / 2 ||
        (ifeng != 0 && ifeng != 1 && ifeng != 2) ||
        !xss_range_valid(table, itie + 1, 2 * ne))
        return ALEA_ERR_PARSE_ERROR;
    thermal->n_inelastic_energies = ne;
    thermal->inelastic_energy = xss_copy(table, itie + 1, ne);
    thermal->inelastic_xs = xss_copy(table, itie + 1 + ne, ne);
    if (!thermal->inelastic_energy || !thermal->inelastic_xs)
        return ALEA_ERR_OUT_OF_MEMORY;
    if (!finite_ascending(thermal->inelastic_energy, ne, 1) ||
        !finite_nonnegative(thermal->inelastic_xs, ne))
        return ALEA_ERR_PARSE_ERROR;
    if (ifeng == 2)
        return decode_continuous_inelastic(table, thermal, itxe, ne);

    int nout = table->nxs[3]; /* NXS(4) */
    int dimension = table->nxs[2]; /* NXS(3) */
    if (nout < 1 || dimension < 0 || dimension == INT_MAX ||
        (ifeng == 1 && nout < 4))
        return ALEA_ERR_PARSE_ERROR;
    int nmu = dimension + 1; /* NXS(4) stores number minus one */
    if (ne > INT_MAX / nout || ne * nout > INT_MAX / (nmu + 1))
        return ALEA_ERR_OVERFLOW;
    int records = ne * nout;
    int stride = nmu + 1;
    if (!xss_range_valid(table, itxe, records * stride))
        return ALEA_ERR_PARSE_ERROR;

    thermal->n_inelastic_outgoing = nout;
    thermal->n_inelastic_cosines = nmu;
    thermal->n_inelastic_outgoing_total = records;
    thermal->inelastic_skewed = ifeng == 1;
    thermal->inelastic_energy_out =
        alea_nuc_malloc((size_t)records * sizeof(double));
    thermal->inelastic_mu =
        alea_nuc_malloc((size_t)records * (size_t)nmu * sizeof(double));
    if (!thermal->inelastic_energy_out || !thermal->inelastic_mu)
        return ALEA_ERR_OUT_OF_MEMORY;

    for (int record = 0; record < records; record++) {
        int base = itxe + record * stride;
        double eout = xss(table, base);
        if (!isfinite(eout) || eout < 0.0) return ALEA_ERR_PARSE_ERROR;
        thermal->inelastic_energy_out[record] = eout;
        double* mu = &thermal->inelastic_mu[(size_t)record * (size_t)nmu];
        for (int k = 0; k < nmu; k++) {
            mu[k] = xss(table, base + 1 + k);
            if (!isfinite(mu[k]) || mu[k] < -1.0 || mu[k] > 1.0)
                return ALEA_ERR_PARSE_ERROR;
        }
        qsort(mu, (size_t)nmu, sizeof(double), compare_double);
    }
    return ALEA_OK;
}

static alea_error_t decode_coherent(const alea_nuc_ace_table_t* table,
                                    alea_nuc_thermal_t* thermal) {
    int itce = table->jxs[3]; /* JXS(4) */
    int n = itce > 0 ? xss_int(table, itce) : 0;
    if (n < 1 || !xss_range_valid(table, itce + 1, 2 * n))
        return ALEA_ERR_PARSE_ERROR;
    thermal->n_coherent_edges = n;
    thermal->coherent_edge = xss_copy(table, itce + 1, n);
    thermal->coherent_factor = xss_copy(table, itce + 1 + n, n);
    if (!thermal->coherent_edge || !thermal->coherent_factor)
        return ALEA_ERR_OUT_OF_MEMORY;
    if (!finite_ascending(thermal->coherent_edge, n, 1) ||
        !finite_nonnegative(thermal->coherent_factor, n))
        return ALEA_ERR_PARSE_ERROR;
    for (int i = 1; i < n; i++)
        if (thermal->coherent_factor[i] < thermal->coherent_factor[i - 1])
            return ALEA_ERR_PARSE_ERROR;
    return ALEA_OK;
}

static alea_error_t decode_incoherent(const alea_nuc_ace_table_t* table,
                                      alea_nuc_thermal_t* thermal,
                                      int mixed) {
    int itce = mixed ? table->jxs[6] : table->jxs[3];
    int itca = mixed ? table->jxs[8] : table->jxs[5];
    int nmu = (mixed ? table->nxs[7] : table->nxs[5]) + 1;
    int ne = itce > 0 ? xss_int(table, itce) : 0;
    if (itce <= 0 || itca <= 0 || ne < 2 || nmu < 1 ||
        ne > INT_MAX / nmu ||
        !xss_range_valid(table, itce + 1, 2 * ne) ||
        !xss_range_valid(table, itca, ne * nmu))
        return ALEA_ERR_PARSE_ERROR;
    thermal->n_incoherent_energies = ne;
    thermal->n_incoherent_cosines = nmu;
    thermal->incoherent_energy = xss_copy(table, itce + 1, ne);
    thermal->incoherent_xs = xss_copy(table, itce + 1 + ne, ne);
    thermal->incoherent_mu = xss_copy(table, itca, ne * nmu);
    if (!thermal->incoherent_energy || !thermal->incoherent_xs ||
        !thermal->incoherent_mu)
        return ALEA_ERR_OUT_OF_MEMORY;
    if (!finite_ascending(thermal->incoherent_energy, ne, 1) ||
        !finite_nonnegative(thermal->incoherent_xs, ne))
        return ALEA_ERR_PARSE_ERROR;
    for (int i = 0; i < ne; i++) {
        double* mu = &thermal->incoherent_mu[(size_t)i * (size_t)nmu];
        for (int k = 0; k < nmu; k++)
            if (!isfinite(mu[k]) || mu[k] < -1.0 || mu[k] > 1.0)
                return ALEA_ERR_PARSE_ERROR;
        qsort(mu, (size_t)nmu, sizeof(double), compare_double);
    }
    return ALEA_OK;
}

alea_error_t alea_nuc_decode_thermal_internal(
    const alea_nuc_ace_table_t* table, alea_nuc_thermal_t** output) {
    if (!table || !output) return ALEA_ERR_NULL_ARG;
    *output = NULL;
    if (table->type != ALEA_NUC_TABLE_THERMAL_SAB)
        return ALEA_ERR_INVALID_ARG;
    alea_nuc_thermal_t* thermal = alea_nuc_calloc(1, sizeof(*thermal));
    if (!thermal) return ALEA_ERR_OUT_OF_MEMORY;
    snprintf(thermal->zaid, sizeof(thermal->zaid), "%s", table->zaid);
    thermal->awr = table->awr;
    thermal->temperature = table->temperature;
    for (int i = 0; i < 16; i++)
        if (table->iz[i] > 0)
            thermal->applicable_zaids[thermal->n_applicable_zaids++] =
                table->iz[i];

    alea_error_t err = decode_inelastic(table, thermal);
    int mode = table->nxs[4]; /* NXS(5) */
    if (err == ALEA_OK && mode != 0 && mode != 3 && mode != 4 && mode != 5)
        err = ALEA_ERR_UNSUPPORTED;
    thermal->elastic_mode = (alea_nuc_thermal_elastic_mode_t)mode;
    if (err == ALEA_OK && (mode == 4 || mode == 5))
        err = decode_coherent(table, thermal);
    if (err == ALEA_OK && (mode == 3 || mode == 5))
        err = decode_incoherent(table, thermal, mode == 5);
    if (err == ALEA_OK && !alea_nuc_thermal_validate_internal(thermal))
        err = ALEA_ERR_PARSE_ERROR;
    if (err != ALEA_OK || table->decode_error) {
        alea_nuc_thermal_free(thermal);
        return table->decode_error ? ALEA_ERR_PARSE_ERROR : err;
    }
    *output = thermal;
    return ALEA_OK;
}

static alea_error_t read_thermal_raw(const alea_nuc_xsdir_t* xsdir,
                                     const alea_nuc_xsdir_entry_t* entry,
                                     alea_nuc_ace_table_t* raw) {
    char path[1024];
    int written;
    if (entry->filename[0] == '/')
        written = snprintf(path, sizeof(path), "%s", entry->filename);
    else if (xsdir->datapath[0])
        written = snprintf(path, sizeof(path), "%s/%s", xsdir->datapath,
                           entry->filename);
    else
        written = snprintf(path, sizeof(path), "%s", entry->filename);
    if (written < 0 || (size_t)written >= sizeof(path)) return ALEA_ERR_OVERFLOW;
    return alea_nuc_ace_read(path, entry->address, entry->file_type, raw);
}

alea_nuc_thermal_t* alea_nuc_load_thermal(const alea_nuc_xsdir_t* xsdir,
                                          const char* zaid) {
    if (!xsdir || !zaid) return NULL;
    const alea_nuc_xsdir_entry_t* entry = alea_nuc_xsdir_find(xsdir, zaid);
    if (!entry || entry->type != ALEA_NUC_TABLE_THERMAL_SAB) return NULL;
    alea_nuc_ace_table_t raw;
    alea_error_t err = read_thermal_raw(xsdir, entry, &raw);
    if (err != ALEA_OK) return NULL;
    alea_nuc_thermal_t* thermal = NULL;
    err = alea_nuc_decode_thermal_internal(&raw, &thermal);
    alea_nuc_ace_free(&raw);
    return err == ALEA_OK ? thermal : NULL;
}

static int energy_interval(const double* grid, int n, double energy,
                           double* fraction) {
    if (energy <= grid[0]) { *fraction = 0.0; return 0; }
    if (energy >= grid[n - 1]) { *fraction = 1.0; return n - 2; }
    int lo = 0, hi = n - 1;
    while (hi - lo > 1) {
        int mid = lo + (hi - lo) / 2;
        if (grid[mid] <= energy) lo = mid;
        else hi = mid;
    }
    *fraction = (energy - grid[lo]) / (grid[lo + 1] - grid[lo]);
    return lo;
}

static double tabulated_xs(const double* grid, const double* xs, int n,
                           double energy) {
    if (!grid || !xs || n < 2 || !isfinite(energy) || energy <= 0.0 ||
        energy > grid[n - 1])
        return 0.0;
    double f;
    int i = energy_interval(grid, n, energy, &f);
    return xs[i] + f * (xs[i + 1] - xs[i]);
}

double alea_nuc_thermal_xs_inelastic(const alea_nuc_thermal_t* thermal,
                                     double energy) {
    if (!thermal) return 0.0;
    return tabulated_xs(thermal->inelastic_energy, thermal->inelastic_xs,
                        thermal->n_inelastic_energies, energy);
}

static double coherent_xs(const alea_nuc_thermal_t* thermal, double energy) {
    if (!thermal || !thermal->coherent_edge || !thermal->coherent_factor ||
        thermal->n_coherent_edges < 1 || !isfinite(energy) || energy <= 0.0)
        return 0.0;
    int edge = -1;
    for (int i = 0; i < thermal->n_coherent_edges; i++) {
        if (thermal->coherent_edge[i] > energy) break;
        edge = i;
    }
    return edge >= 0 ? thermal->coherent_factor[edge] / energy : 0.0;
}

double alea_nuc_thermal_xs_elastic(const alea_nuc_thermal_t* thermal,
                                   double energy) {
    if (!thermal) return 0.0;
    double value = 0.0;
    if (thermal->elastic_mode == ALEA_NUC_THERMAL_ELASTIC_COHERENT ||
        thermal->elastic_mode == ALEA_NUC_THERMAL_ELASTIC_MIXED)
        value += coherent_xs(thermal, energy);
    if (thermal->elastic_mode == ALEA_NUC_THERMAL_ELASTIC_INCOHERENT ||
        thermal->elastic_mode == ALEA_NUC_THERMAL_ELASTIC_MIXED)
        value += tabulated_xs(thermal->incoherent_energy,
                              thermal->incoherent_xs,
                              thermal->n_incoherent_energies, energy);
    return value;
}

double alea_nuc_thermal_xs_total(const alea_nuc_thermal_t* thermal,
                                 double energy) {
    return alea_nuc_thermal_xs_inelastic(thermal, energy) +
           alea_nuc_thermal_xs_elastic(thermal, energy);
}

static alea_error_t draw(alea_nuc_random_fn random, void* context,
                         double* value) {
    if (!random || !value) return ALEA_ERR_NULL_ARG;
    double sampled = random(context);
    if (!isfinite(sampled) || sampled < 0.0 || sampled >= 1.0)
        return ALEA_ERR_INVALID_STATE;
    *value = sampled;
    return ALEA_OK;
}

static int sample_outgoing_index(int n, int skewed, double xi) {
    if (!skewed) return (int)(xi * n);
    double r = xi * (n - 3);
    if (r > 1.0) return (int)r + 1;
    if (r > 0.6) return n - 2;
    if (r > 0.5) return n - 1;
    if (r > 0.1) return 1;
    return 0;
}

static alea_error_t sample_discrete_inelastic(
    const alea_nuc_thermal_t* thermal, double energy,
    alea_nuc_random_fn random, void* context, double* energy_out,
    double* mu) {
    double xi_e, xi_mu;
    alea_error_t err = draw(random, context, &xi_e);
    if (err == ALEA_OK) err = draw(random, context, &xi_mu);
    if (err != ALEA_OK) return err;
    double f;
    int i = energy_interval(thermal->inelastic_energy,
                            thermal->n_inelastic_energies, energy, &f);
    int j = sample_outgoing_index(thermal->n_inelastic_outgoing,
                                  thermal->inelastic_skewed, xi_e);
    int k = (int)(xi_mu * thermal->n_inelastic_cosines);
    size_t r0 = (size_t)i * thermal->n_inelastic_outgoing + (size_t)j;
    size_t r1 = r0 + (size_t)thermal->n_inelastic_outgoing;
    *energy_out = thermal->inelastic_energy_out[r0] + f *
        (thermal->inelastic_energy_out[r1] - thermal->inelastic_energy_out[r0]);
    size_t m0 = r0 * (size_t)thermal->n_inelastic_cosines + (size_t)k;
    size_t m1 = r1 * (size_t)thermal->n_inelastic_cosines + (size_t)k;
    *mu = thermal->inelastic_mu[m0] +
          f * (thermal->inelastic_mu[m1] - thermal->inelastic_mu[m0]);
    return ALEA_OK;
}

static alea_error_t sample_continuous_inelastic(
    const alea_nuc_thermal_t* thermal, double energy,
    alea_nuc_random_fn random, void* context, double* energy_out,
    double* mu) {
    double incident_fraction;
    int incident_index = energy_interval(
        thermal->inelastic_energy, thermal->n_inelastic_energies,
        energy, &incident_fraction);
    int selected_incident = incident_fraction > 0.5 ?
        incident_index + 1 : incident_index;
    int begin = thermal->inelastic_outgoing_offset[selected_incident];
    int end = thermal->inelastic_outgoing_offset[selected_incident + 1];

    double xi_energy;
    alea_error_t err = draw(random, context, &xi_energy);
    if (err != ALEA_OK) return err;
    double target = xi_energy * thermal->inelastic_cdf[end - 1];
    int outgoing = begin;
    while (outgoing < end - 2 &&
           target >= thermal->inelastic_cdf[outgoing + 1])
        outgoing++;
    double e0 = thermal->inelastic_energy_out[outgoing];
    double e1 = thermal->inelastic_energy_out[outgoing + 1];
    double p0 = thermal->inelastic_pdf[outgoing];
    double p1 = thermal->inelastic_pdf[outgoing + 1];
    double c0 = thermal->inelastic_cdf[outgoing];
    double c1 = thermal->inelastic_cdf[outgoing + 1];
    double delta_cdf = c1 - c0;
    double slope = (p1 - p0) / (e1 - e0);
    if (!(delta_cdf > 0.0)) return ALEA_ERR_INVALID_STATE;
    if (slope == 0.0) {
        if (!(p0 > 0.0)) return ALEA_ERR_INVALID_STATE;
        *energy_out = e0 + (target - c0) / p0;
    } else {
        double discriminant = p0 * p0 + 2.0 * slope * (target - c0);
        if (discriminant < 0.0 && discriminant > -1e-14 * p0 * p0)
            discriminant = 0.0;
        if (discriminant < 0.0) return ALEA_ERR_INVALID_STATE;
        *energy_out = e0 + (sqrt(discriminant) - p0) / slope;
    }
    *energy_out = fmin(e1, fmax(e0, *energy_out));

    double selected_energy = thermal->inelastic_energy[selected_incident];
    if (*energy_out < 0.5 * selected_energy)
        *energy_out *= 2.0 * energy / selected_energy - 1.0;
    else
        *energy_out += energy - selected_energy;

    double outgoing_fraction = (target - c0) / delta_cdf;
    outgoing_fraction = fmin(1.0, fmax(0.0, outgoing_fraction));
    double xi_mu, xi_smear;
    if ((err = draw(random, context, &xi_mu)) != ALEA_OK ||
        (err = draw(random, context, &xi_smear)) != ALEA_OK)
        return err;
    int nmu = thermal->n_inelastic_cosines;
    int k = (int)(xi_mu * nmu);
    const double* mu0 = &thermal->inelastic_mu[
        (size_t)outgoing * (size_t)nmu];
    const double* mu1 = &thermal->inelastic_mu[
        (size_t)(outgoing + 1) * (size_t)nmu];
    *mu = mu0[k] + outgoing_fraction * (mu1[k] - mu0[k]);
    double left = k == 0 ? -2.0 - *mu :
        mu0[k - 1] + outgoing_fraction * (mu1[k - 1] - mu0[k - 1]);
    double right = k == nmu - 1 ? 2.0 - *mu :
        mu0[k + 1] + outgoing_fraction * (mu1[k + 1] - mu0[k + 1]);
    double width = fmin(*mu - left, right - *mu);
    if (width < 0.0) return ALEA_ERR_INVALID_STATE;
    *mu += width * (xi_smear - 0.5);
    return ALEA_OK;
}

static alea_error_t sample_inelastic(const alea_nuc_thermal_t* thermal,
                                     double energy, alea_nuc_random_fn random,
                                     void* context, double* energy_out,
                                     double* mu) {
    if (thermal->inelastic_continuous)
        return sample_continuous_inelastic(
            thermal, energy, random, context, energy_out, mu);
    return sample_discrete_inelastic(
        thermal, energy, random, context, energy_out, mu);
}

static alea_error_t sample_coherent(const alea_nuc_thermal_t* thermal,
                                    double energy, alea_nuc_random_fn random,
                                    void* context, double* mu) {
    int last = -1;
    for (int i = 0; i < thermal->n_coherent_edges; i++) {
        if (thermal->coherent_edge[i] > energy) break;
        last = i;
    }
    if (last < 0 || thermal->coherent_factor[last] <= 0.0)
        return ALEA_ERR_NOT_FOUND;
    double xi;
    alea_error_t err = draw(random, context, &xi);
    if (err != ALEA_OK) return err;
    double target = xi * thermal->coherent_factor[last];
    int edge = 0;
    while (edge < last && thermal->coherent_factor[edge] < target) edge++;
    *mu = 1.0 - 2.0 * thermal->coherent_edge[edge] / energy;
    return ALEA_OK;
}

static alea_error_t sample_incoherent(const alea_nuc_thermal_t* thermal,
                                      double energy, alea_nuc_random_fn random,
                                      void* context, double* mu) {
    double xi;
    alea_error_t err = draw(random, context, &xi);
    if (err != ALEA_OK) return err;
    double f;
    int i = energy_interval(thermal->incoherent_energy,
                            thermal->n_incoherent_energies, energy, &f);
    int k = (int)(xi * thermal->n_incoherent_cosines);
    size_t m0 = (size_t)i * thermal->n_incoherent_cosines + (size_t)k;
    size_t m1 = m0 + (size_t)thermal->n_incoherent_cosines;
    *mu = thermal->incoherent_mu[m0] +
          f * (thermal->incoherent_mu[m1] - thermal->incoherent_mu[m0]);
    return ALEA_OK;
}

alea_error_t alea_nuc_sample_thermal_collision(
    const alea_nuc_thermal_t* thermal,
    const alea_nuc_particle_state_t* incident, alea_nuc_random_fn random,
    void* random_context, alea_nuc_collision_result_t* output) {
    if (!thermal || !incident || !random || !output) return ALEA_ERR_NULL_ARG;
    if (incident->type != ALEA_NUC_PARTICLE_NEUTRON ||
        !isfinite(incident->energy) || incident->energy <= 0.0 ||
        !isfinite(incident->weight) || incident->weight <= 0.0 ||
        !isfinite(incident->time))
        return ALEA_ERR_INVALID_ARG;
    double norm2 = 0.0;
    for (int i = 0; i < 3; i++) {
        if (!isfinite(incident->direction[i])) return ALEA_ERR_INVALID_ARG;
        norm2 += incident->direction[i] * incident->direction[i];
    }
    if (fabs(norm2 - 1.0) > 1e-8) return ALEA_ERR_INVALID_ARG;

    double inelastic = alea_nuc_thermal_xs_inelastic(thermal, incident->energy);
    double coherent = coherent_xs(thermal, incident->energy);
    double incoherent = 0.0;
    if (thermal->elastic_mode == ALEA_NUC_THERMAL_ELASTIC_INCOHERENT ||
        thermal->elastic_mode == ALEA_NUC_THERMAL_ELASTIC_MIXED)
        incoherent = tabulated_xs(thermal->incoherent_energy,
                                  thermal->incoherent_xs,
                                  thermal->n_incoherent_energies,
                                  incident->energy);
    if (thermal->elastic_mode == ALEA_NUC_THERMAL_ELASTIC_NONE ||
        thermal->elastic_mode == ALEA_NUC_THERMAL_ELASTIC_INCOHERENT)
        coherent = 0.0;
    if (inelastic > 0.0 &&
        (!thermal->inelastic_energy_out || !thermal->inelastic_mu ||
         (thermal->inelastic_continuous &&
          (!thermal->inelastic_outgoing_offset || !thermal->inelastic_pdf ||
           !thermal->inelastic_cdf ||
           thermal->n_inelastic_outgoing_total < 2)) ||
         (!thermal->inelastic_continuous &&
          thermal->n_inelastic_outgoing < 1) ||
         thermal->n_inelastic_cosines < 1))
        return ALEA_ERR_INVALID_STATE;
    if (incoherent > 0.0 &&
        (!thermal->incoherent_mu || thermal->n_incoherent_cosines < 1))
        return ALEA_ERR_INVALID_STATE;
    double total = inelastic + coherent + incoherent;
    if (!(total > 0.0) || !isfinite(total)) return ALEA_ERR_NOT_FOUND;

    double xi_channel;
    alea_error_t err = draw(random, random_context, &xi_channel);
    if (err != ALEA_OK) return err;
    double energy_out = incident->energy;
    double mu;
    int mt;
    double selected = xi_channel * total;
    if (selected < inelastic) {
        mt = 4;
        err = sample_inelastic(thermal, incident->energy, random,
                               random_context, &energy_out, &mu);
    } else if (selected < inelastic + coherent) {
        mt = 2;
        err = sample_coherent(thermal, incident->energy, random,
                              random_context, &mu);
    } else {
        mt = 2;
        err = sample_incoherent(thermal, incident->energy, random,
                                random_context, &mu);
    }
    if (err != ALEA_OK) return err;
    if (!isfinite(energy_out) || energy_out <= 0.0 ||
        !isfinite(mu) || mu < -1.0 || mu > 1.0)
        return ALEA_ERR_INVALID_STATE;
    double xi_phi;
    err = draw(random, random_context, &xi_phi);
    if (err != ALEA_OK) return err;

    alea_nuc_collision_result_t result;
    memset(&result, 0, sizeof(result));
    result.outcome = ALEA_NUC_OUTCOME_SCATTERED;
    result.component_index = -1;
    result.mt = mt;
    result.mu_cm = NAN;
    result.mu_lab = mu;
    result.deposition_available = true;
    result.local_energy_deposition = incident->energy - energy_out;
    result.outgoing = *incident;
    result.outgoing.energy = energy_out;
    alea_nuc_rotate_direction_internal(incident->direction, mu,
                                       2.0 * M_PI * xi_phi,
                                       result.outgoing.direction);
    *output = result;
    return ALEA_OK;
}
