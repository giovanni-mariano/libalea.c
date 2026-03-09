// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0


/**
 * @file nuc_plot.c
 * @brief Nuclear data visualizer — SVG plot generator
 *
 * Usage:
 *   nuc_plot --xsdir <path> --zaid <zaid> [options]
 *
 * Plot types (--plot):
 *   xs         Cross sections vs energy (default)
 *   reactions  Per-reaction cross sections
 *   angular    Angular distribution at selected energy
 *   fission    Fission spectrum χ(E)
 *   nubar      ν̄(E) vs energy
 *   doppler    Cross sections at multiple temperatures
 *   heating    Heating numbers vs energy
 *   photon-xs  Photon cross section components vs energy
 *   photon-ff  Incoherent/coherent form factors vs momentum transfer
 */

#include "alea_nucdata.h"
#include "util/alea_svg.h"
#include "util/compat.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* Simple LCG PRNG */
static inline double local_rng(uint64_t* state) {
    *state = *state * 6364136223846793005ULL + 1442695040888963407ULL;
    return (double)(*state >> 11) / (double)(1ULL << 53);
}

/* Maximum MT numbers the user can request */
#define MAX_MT_FILTER 32

/* Parsed --mt filter */
typedef struct {
    int mt[MAX_MT_FILTER];
    int count;              /* 0 = no filter (show all) */
} mt_filter_t;

/** Parse comma-separated MT list, e.g. "18,102,16" */
static void mt_filter_parse(mt_filter_t* f, const char* str) {
    f->count = 0;
    if (!str) return;
    char buf[256];
    strncpy(buf, str, sizeof(buf) - 1);
    buf[sizeof(buf) - 1] = '\0';
    char* tok = strtok(buf, ",");
    while (tok && f->count < MAX_MT_FILTER) {
        f->mt[f->count++] = atoi(tok);
        tok = strtok(NULL, ",");
    }
}

/** Check if MT is in the filter (returns 1 if no filter or if MT matches) */
static int mt_filter_match(const mt_filter_t* f, int mt) {
    if (f->count == 0) return 1; /* no filter = show all */
    for (int i = 0; i < f->count; i++)
        if (f->mt[i] == mt) return 1;
    return 0;
}

/* Forward declaration */
static const char* mt_name(int mt);

/** Print available reactions for a nuclide */
static void print_reactions(const alea_nuc_nuclide_t* nuc) {
    printf("Available reactions for %s:\n", nuc->zaid);
    printf("  MT=2    elastic\n");
    for (int i = 0; i < nuc->n_reactions; i++) {
        const alea_nuc_reaction_t* r = &nuc->reactions[i];
        if (r->n_energies <= 0) continue;
        printf("  MT=%-4d %s", r->mt, mt_name(r->mt));
        if (r->q_value != 0.0) printf("  Q=%.4f MeV", r->q_value);
        printf("\n");
    }
}

/* ============================================================================
 * PLOT: Cross Sections vs Energy
 * ============================================================================ */

static int plot_xs(alea_nuc_nuclide_t* nuc, const mt_filter_t* filter,
                    const char* output) {
    int N = nuc->n_energies;

    alea_svg_plot_t plot;
    alea_svg_plot_init(&plot);
    plot.x_label = "Energy (MeV)";
    plot.y_label = "Cross section (b)";

    char title[128];
    snprintf(title, sizeof(title), "%s — Cross Sections", nuc->zaid);
    plot.title = title;

    if (filter->count == 0) {
        /* No filter: show total, elastic, absorption + fission if present */
        alea_svg_plot_add(&plot, nuc->energy, nuc->sigma_total, N,
                     "Total", "#2266cc");
        alea_svg_plot_add(&plot, nuc->energy, nuc->sigma_elastic, N,
                     "Elastic", "#22aa44");
        alea_svg_plot_add(&plot, nuc->energy, nuc->sigma_abs, N,
                     "Absorption", "#cc2222");

        for (int i = 0; i < nuc->n_reactions; i++) {
            if (nuc->reactions[i].mt == 18 && nuc->reactions[i].n_energies > 0) {
                double* sigma_f = calloc((size_t)N, sizeof(double));
                if (!sigma_f) break;
                int ie_start = nuc->reactions[i].threshold_index - 1;
                for (int j = 0; j < nuc->reactions[i].n_energies && ie_start + j < N; j++)
                    sigma_f[ie_start + j] = nuc->reactions[i].xs[j];
                alea_svg_plot_add(&plot, nuc->energy, sigma_f, N,
                             "Fission (MT=18)", "#cc8800");
                break;
            }
        }
    } else {
        /* User selected specific MTs */
        for (int f = 0; f < filter->count && plot.n_curves < ALEA_SVG_MAX_CURVES; f++) {
            int mt = filter->mt[f];

            if (mt == 0) {
                /* MT=0 is a shorthand for total */
                alea_svg_plot_add(&plot, nuc->energy, nuc->sigma_total, N,
                             "Total", NULL);
            } else if (mt == 2) {
                alea_svg_plot_add(&plot, nuc->energy, nuc->sigma_elastic, N,
                             "Elastic (MT=2)", NULL);
            } else {
                /* Find reaction */
                for (int i = 0; i < nuc->n_reactions; i++) {
                    alea_nuc_reaction_t* r = &nuc->reactions[i];
                    if (r->mt != mt || r->n_energies <= 0) continue;

                    double* sigma = calloc((size_t)N, sizeof(double));
                    if (!sigma) break;
                    int ie_start = r->threshold_index - 1;
                    for (int j = 0; j < r->n_energies && ie_start + j < N; j++)
                        sigma[ie_start + j] = r->xs[j];

                    char label[32];
                    snprintf(label, sizeof(label), "%s (MT=%d)", mt_name(mt), mt);
                    alea_svg_plot_add(&plot, nuc->energy, sigma, N,
                                 alea_strdup(label), NULL);
                    break;
                }
            }
        }
    }

    int rc = alea_svg_plot_write(&plot, output);
    if (rc == 0)
        printf("Written: %s\n", output);
    return rc;
}

/* ============================================================================
 * PLOT: Per-Reaction Cross Sections
 * ============================================================================ */

static const char* mt_name(int mt) {
    switch (mt) {
    case 2: return "elastic";
    case 4: return "(n,n')";
    case 5: return "(n,misc)";
    case 16: return "(n,2n)";
    case 17: return "(n,3n)";
    case 18: return "fission";
    case 22: return "(n,na)";
    case 28: return "(n,np)";
    case 37: return "(n,4n)";
    case 51: return "(n,n'1)";
    case 52: return "(n,n'2)";
    case 53: return "(n,n'3)";
    case 91: return "(n,n'c)";
    case 102: return "(n,g)";
    case 103: return "(n,p)";
    case 104: return "(n,d)";
    case 105: return "(n,t)";
    case 107: return "(n,a)";
    default: break;
    }
    static char buf[16];
    snprintf(buf, sizeof(buf), "MT=%d", mt);
    return buf;
}

static int plot_reactions(alea_nuc_nuclide_t* nuc, const mt_filter_t* filter,
                           const char* output) {
    int N = nuc->n_energies;

    alea_svg_plot_t plot;
    alea_svg_plot_init(&plot);
    plot.x_label = "Energy (MeV)";
    plot.y_label = "Cross section (b)";

    char title[128];
    snprintf(title, sizeof(title), "%s — Reaction Cross Sections", nuc->zaid);
    plot.title = title;

    /* Add total as reference */
    alea_svg_plot_add(&plot, nuc->energy, nuc->sigma_total, N,
                 "Total", "#cccccc");
    plot.curves[0].width = 2.5;

    /* Add individual reactions (filtered if --mt given) */
    int added = 0;
    for (int i = 0; i < nuc->n_reactions && added < ALEA_SVG_MAX_CURVES - 2; i++) {
        alea_nuc_reaction_t* r = &nuc->reactions[i];
        if (r->n_energies <= 0 || !r->xs) continue;
        if (!mt_filter_match(filter, r->mt)) continue;

        /* Build on full grid */
        double* sigma = calloc((size_t)N, sizeof(double));
        if (!sigma) continue;
        int ie_start = r->threshold_index - 1;
        for (int j = 0; j < r->n_energies && ie_start + j < N; j++)
            sigma[ie_start + j] = r->xs[j];

        char label[32];
        snprintf(label, sizeof(label), "%s", mt_name(r->mt));

        alea_svg_plot_add(&plot, nuc->energy, sigma, N, alea_strdup(label), NULL);
        added++;
    }

    int rc = alea_svg_plot_write(&plot, output);
    if (rc == 0)
        printf("Written: %s (%d reactions)\n", output, added);
    return rc;
}

/* ============================================================================
 * PLOT: Angular Distribution
 * ============================================================================ */

static int plot_angular(alea_nuc_nuclide_t* nuc, double energy, const char* output) {
    if (!nuc->elastic_angular) {
        fprintf(stderr, "No elastic angular distribution data\n");
        return -1;
    }

    /* Sample angular distribution to build histogram */
    int N_bins = 200;
    double* cosines = malloc((size_t)N_bins * sizeof(double));
    double* pdf = calloc((size_t)N_bins, sizeof(double));
    if (!cosines || !pdf) return -1;

    for (int i = 0; i < N_bins; i++)
        cosines[i] = -1.0 + 2.0 * (i + 0.5) / N_bins;

    /* Monte Carlo sampling */
    int N_samples = 500000;
    uint64_t rng_state = 0xDEADBEEF12345678ULL;
    double bin_width = 2.0 / N_bins;

    for (int s = 0; s < N_samples; s++) {
        rng_state = rng_state * 6364136223846793005ULL + 1442695040888963407ULL;
        double xi = (double)(rng_state >> 11) / (double)(1ULL << 53);

        /* Inline angle sampling — find energy bracket */
        alea_nuc_angular_dist_t* ang = nuc->elastic_angular;
        double f;
        int ie = alea_nuc_energy_lookup(ang->energy, ang->n_energies, energy, &f);
        if (ie < 0) ie = 0;

        /* Use lower energy point for simplicity */
        alea_nuc_angular_point_t* ap = &ang->data[ie];
        double mu;

        if (ap->type == ALEA_NUC_ANG_ISOTROPIC || !ap->cosine) {
            mu = 2.0 * xi - 1.0;
        } else if (ap->cdf && ap->n_cosines > 0) {
            /* Sample from CDF */
            int lo = 0, hi = ap->n_cosines - 2;
            while (lo < hi) {
                int mid = (lo + hi + 1) / 2;
                if (ap->cdf[mid] <= xi) lo = mid;
                else hi = mid - 1;
            }
            double dcdf = ap->cdf[lo + 1] - ap->cdf[lo];
            double frac = (dcdf > 0.0) ? (xi - ap->cdf[lo]) / dcdf : 0.0;
            mu = ap->cosine[lo] + frac * (ap->cosine[lo + 1] - ap->cosine[lo]);
        } else {
            /* Equiprobable bins */
            int bin = (int)(xi * 32);
            if (bin >= 32) bin = 31;
            double frac = xi * 32.0 - bin;
            mu = ap->cosine[bin] + frac * (ap->cosine[bin + 1] - ap->cosine[bin]);
        }

        if (mu < -1.0) mu = -1.0;
        if (mu > 1.0) mu = 1.0;

        int b = (int)((mu + 1.0) / bin_width);
        if (b >= N_bins) b = N_bins - 1;
        if (b < 0) b = 0;
        pdf[b] += 1.0;
    }

    /* Normalize to PDF */
    for (int i = 0; i < N_bins; i++)
        pdf[i] /= (N_samples * bin_width);

    alea_svg_plot_t plot;
    alea_svg_plot_init(&plot);
    plot.x_scale = ALEA_SVG_SCALE_LINEAR;
    plot.y_scale = ALEA_SVG_SCALE_LINEAR;
    plot.x_label = "Scattering cosine (mu)";
    plot.y_label = "P(mu)";

    char title[128];
    snprintf(title, sizeof(title), "%s — Elastic Angular Distribution at %.3g MeV",
             nuc->zaid, energy);
    plot.title = title;

    alea_svg_plot_add(&plot, cosines, pdf, N_bins, NULL, "#2266cc");

    /* Add isotropic reference */
    double* iso = malloc((size_t)N_bins * sizeof(double));
    if (iso) {
        for (int i = 0; i < N_bins; i++) iso[i] = 0.5;
        alea_svg_plot_add(&plot, cosines, iso, N_bins, "Isotropic", "#999999");
        plot.curves[1].dashed = 1;
    }

    int rc = alea_svg_plot_write(&plot, output);
    if (rc == 0)
        printf("Written: %s\n", output);

    free(cosines);
    free(pdf);
    free(iso);
    return rc;
}

/* ============================================================================
 * PLOT: Fission Spectrum
 * ============================================================================ */

/**
 * Sample Watt fission spectrum: f(E') = C * exp(-E'/a) * sinh(√(b·E'))
 * Pure math — no nuclear data needed.
 */
static double sample_watt_inline(double a, double b, double E_restrict,
                                  uint64_t* rng_state) {
    double K = 1.0 + b * a / 8.0;
    double Kk = K * K - 1.0;
    if (Kk < 0.0) Kk = 0.0;
    double L = (K + sqrt(Kk)) / a;
    double M_param = a * L - 1.0;

    for (int iter = 0; iter < 1000; iter++) {
        double r1 = local_rng(rng_state);
        double r2 = local_rng(rng_state);

        double x = -log(r1) / L;
        double eta = -log(r2);
        double diff = sqrt(b * x) - M_param;

        if (eta >= diff * diff / 4.0) {
            if (x > 0.0 && x <= E_restrict)
                return x;
        }
    }
    return a;
}

static int plot_fission(alea_nuc_nuclide_t* nuc, const char* output) {
    if (!nuc->fission) {
        fprintf(stderr, "Nuclide %s is not fissile\n", nuc->zaid);
        return -1;
    }

    int N_bins = 200;
    double E_max = 15.0; /* MeV */
    double bin_width = E_max / N_bins;
    double* energies = malloc((size_t)N_bins * sizeof(double));
    double* spectrum = calloc((size_t)N_bins, sizeof(double));
    if (!energies || !spectrum) return -1;

    for (int i = 0; i < N_bins; i++)
        energies[i] = (i + 0.5) * bin_width;

    /* Sample Watt fission spectrum (U-235 thermal: a≈0.988, b≈2.249) */
    double a = 0.988, b = 2.249;
    int N_samples = 500000;
    uint64_t rng_state = 0xCAFEBABE42424242ULL;

    for (int s = 0; s < N_samples; s++) {
        double E_out = sample_watt_inline(a, b, E_max, &rng_state);

        if (E_out > 0.0 && E_out < E_max) {
            int bin = (int)(E_out / bin_width);
            if (bin >= N_bins) bin = N_bins - 1;
            spectrum[bin] += 1.0;
        }
    }

    /* Normalize */
    double total = 0.0;
    for (int i = 0; i < N_bins; i++) total += spectrum[i] * bin_width;
    if (total > 0.0)
        for (int i = 0; i < N_bins; i++) spectrum[i] /= total;

    /* Watt analytic reference */
    double* watt_ref = malloc((size_t)N_bins * sizeof(double));
    if (watt_ref) {
        double sum = 0.0;
        for (int i = 0; i < N_bins; i++) {
            double E = energies[i];
            watt_ref[i] = exp(-E / a) * sinh(sqrt(b * E));
            sum += watt_ref[i] * bin_width;
        }
        if (sum > 0.0)
            for (int i = 0; i < N_bins; i++) watt_ref[i] /= sum;
    }

    alea_svg_plot_t plot;
    alea_svg_plot_init(&plot);
    plot.x_scale = ALEA_SVG_SCALE_LINEAR;
    plot.y_scale = ALEA_SVG_SCALE_LINEAR;
    plot.x_label = "Energy (MeV)";
    plot.y_label = "chi(E) (1/MeV)";
    plot.x_min = 0.0;
    plot.x_max = E_max;

    char title[128];
    snprintf(title, sizeof(title), "%s — Fission Spectrum (thermal)", nuc->zaid);
    plot.title = title;

    alea_svg_plot_add(&plot, energies, spectrum, N_bins, "Sampled (Watt)", "#2266cc");
    if (watt_ref) {
        alea_svg_plot_add(&plot, energies, watt_ref, N_bins, "Watt (a=0.988, b=2.249)", "#cc2222");
        plot.curves[1].dashed = 1;
    }

    int rc = alea_svg_plot_write(&plot, output);
    if (rc == 0)
        printf("Written: %s (%d samples)\n", output, N_samples);

    free(energies);
    free(spectrum);
    free(watt_ref);
    return rc;
}

/* ============================================================================
 * PLOT: ν̄(E) vs Energy
 * ============================================================================ */

static int plot_nubar(alea_nuc_nuclide_t* nuc, const char* output) {
    if (!nuc->fission) {
        fprintf(stderr, "Nuclide %s is not fissile\n", nuc->zaid);
        return -1;
    }

    int N = 500;
    double* energies = malloc((size_t)N * sizeof(double));
    double* nubar = malloc((size_t)N * sizeof(double));
    if (!energies || !nubar) return -1;

    /* Log-spaced energy grid from 1e-9 to 20 MeV */
    for (int i = 0; i < N; i++) {
        double t = (double)i / (N - 1);
        energies[i] = 1e-9 * pow(20.0 / 1e-9, t);
        nubar[i] = alea_nuc_nu_bar(nuc, energies[i]);
    }

    alea_svg_plot_t plot;
    alea_svg_plot_init(&plot);
    plot.y_scale = ALEA_SVG_SCALE_LINEAR;
    plot.x_label = "Energy (MeV)";
    plot.y_label = "nu-bar";

    char title[128];
    snprintf(title, sizeof(title), "%s — Average Neutrons per Fission", nuc->zaid);
    plot.title = title;

    alea_svg_plot_add(&plot, energies, nubar, N, NULL, "#2266cc");

    int rc = alea_svg_plot_write(&plot, output);
    if (rc == 0)
        printf("Written: %s\n", output);

    free(energies);
    free(nubar);
    return rc;
}


/* ============================================================================
 * PLOT: Doppler Broadening Comparison
 * ============================================================================ */

static int plot_doppler(const alea_nuc_xsdir_t* xsdir, const char* zaid,
                         double temp_K, const char* output) {
    /* Load nuclide for cold XS */
    alea_nuc_nuclide_t* nuc_cold = alea_nuc_load_nuclide(xsdir, zaid);
    if (!nuc_cold) {
        fprintf(stderr, "Failed to load %s\n", zaid);
        return -1;
    }

    /* Save cold XS on a sub-grid around resonances */
    int N = 2000;
    double* energies = malloc((size_t)N * sizeof(double));
    double* sigma_cold = malloc((size_t)N * sizeof(double));
    double* sigma_hot = malloc((size_t)N * sizeof(double));
    if (!energies || !sigma_cold || !sigma_hot) return -1;

    /* Focus on resolved resonance region: 1e-5 to 1e-2 MeV (10 eV to 10 keV) */
    for (int i = 0; i < N; i++) {
        double t = (double)i / (N - 1);
        energies[i] = 1e-6 * pow(1e-2 / 1e-6, t);
        sigma_cold[i] = alea_nuc_xs_total(nuc_cold, energies[i]);
    }

    /* Broaden to target temperature */
    double kT = temp_K * 8.617333e-11; /* K to MeV */
    alea_error_t err = alea_nuc_doppler_broaden(nuc_cold, kT);
    if (err != ALEA_OK) {
        fprintf(stderr, "Doppler broadening failed: %s\n", alea_error_string(err));
        free(energies); free(sigma_cold); free(sigma_hot);
        return -1;
    }

    for (int i = 0; i < N; i++)
        sigma_hot[i] = alea_nuc_xs_total(nuc_cold, energies[i]);

    alea_svg_plot_t plot;
    alea_svg_plot_init(&plot);
    plot.x_label = "Energy (MeV)";
    plot.y_label = "Cross section (b)";

    char title[128];
    char label_cold[64], label_hot[64];
    snprintf(title, sizeof(title), "%s — Doppler Broadening", zaid);
    snprintf(label_cold, sizeof(label_cold), "Original (%.0f K)",
             nuc_cold->temperature / 8.617333e-11); /* approx */
    snprintf(label_hot, sizeof(label_hot), "Broadened (%.0f K)", temp_K);
    plot.title = title;

    /* Cold will show the data BEFORE broadening only if we had saved it.
     * Since we broadened in-place, sigma_cold was captured before broadening. */
    alea_svg_plot_add(&plot, energies, sigma_cold, N, "Original", "#2266cc");
    alea_svg_plot_add(&plot, energies, sigma_hot, N, alea_strdup(label_hot), "#cc2222");

    int rc = alea_svg_plot_write(&plot, output);
    if (rc == 0)
        printf("Written: %s\n", output);

    free(energies);
    free(sigma_cold);
    free(sigma_hot);
    return rc;
}

/* ============================================================================
 * PLOT: Heating Numbers vs Energy
 * ============================================================================ */

static int plot_heating(alea_nuc_nuclide_t* nuc, const char* output) {
    if (!nuc->heating && !(nuc->photon && nuc->photon->heating)) {
        fprintf(stderr, "Nuclide %s has no heating data\n", nuc->zaid);
        return -1;
    }

    /* Use photon energy grid for photon tables, main grid otherwise */
    const double* egrid = (nuc->particle == ALEA_NUC_PARTICLE_PHOTON && nuc->photon)
                          ? nuc->photon->energy : nuc->energy;
    int N = (nuc->particle == ALEA_NUC_PARTICLE_PHOTON && nuc->photon)
            ? nuc->photon->n_energies : nuc->n_energies;

    double* h_vals = malloc((size_t)N * sizeof(double));
    double* h_per_coll = malloc((size_t)N * sizeof(double));
    if (!h_vals || !h_per_coll) { free(h_vals); free(h_per_coll); return -1; }

    for (int i = 0; i < N; i++) {
        h_vals[i] = alea_nuc_xs_heating(nuc, egrid[i]);
        h_per_coll[i] = alea_nuc_heating_per_collision(nuc, egrid[i]);
    }

    alea_svg_plot_t plot;
    alea_svg_plot_init(&plot);
    plot.x_label = "Energy (MeV)";
    plot.y_label = "Heating number (MeV/collision)";

    char title[128];
    snprintf(title, sizeof(title), "%s — Heating Numbers", nuc->zaid);
    plot.title = title;

    alea_svg_plot_add(&plot, egrid, h_vals, N, "Heating (MeV·b)", "#2266cc");
    alea_svg_plot_add(&plot, egrid, h_per_coll, N, "Heating/collision (MeV)", "#cc2222");

    int rc = alea_svg_plot_write(&plot, output);
    if (rc == 0)
        printf("Written: %s\n", output);

    free(h_vals);
    free(h_per_coll);
    return rc;
}

/* ============================================================================
 * PLOT: Photon Cross Section Components
 * ============================================================================ */

static int plot_photon_xs(alea_nuc_nuclide_t* nuc, const char* output) {
    if (!nuc->photon) {
        fprintf(stderr, "Nuclide %s has no photon data\n", nuc->zaid);
        return -1;
    }
    alea_nuc_photon_data_t* ph = nuc->photon;
    int N = ph->n_energies;

    alea_svg_plot_t plot;
    alea_svg_plot_init(&plot);
    plot.x_label = "Energy (MeV)";
    plot.y_label = "Cross section (b)";

    char title[128];
    snprintf(title, sizeof(title), "%s — Photon Cross Sections", nuc->zaid);
    plot.title = title;

    /* Total (from main grid) */
    if (nuc->sigma_total)
        alea_svg_plot_add(&plot, nuc->energy, nuc->sigma_total, nuc->n_energies,
                     "Total", "#2266cc");

    alea_svg_plot_add(&plot, ph->energy, ph->sigma_incoherent, N,
                 "Compton (incoherent)", "#22aa44");
    alea_svg_plot_add(&plot, ph->energy, ph->sigma_coherent, N,
                 "Rayleigh (coherent)", "#cc8800");
    alea_svg_plot_add(&plot, ph->energy, ph->sigma_photoelectric, N,
                 "Photoelectric", "#cc2222");
    alea_svg_plot_add(&plot, ph->energy, ph->sigma_pair, N,
                 "Pair production", "#8822cc");

    int rc = alea_svg_plot_write(&plot, output);
    if (rc == 0)
        printf("Written: %s\n", output);
    return rc;
}

/* ============================================================================
 * PLOT: Photon Form Factors
 * ============================================================================ */

static int plot_photon_ff(alea_nuc_nuclide_t* nuc, const char* output) {
    if (!nuc->photon) {
        fprintf(stderr, "Nuclide %s has no photon data\n", nuc->zaid);
        return -1;
    }
    alea_nuc_photon_data_t* ph = nuc->photon;

    if (ph->n_incoherent_ff <= 0 && ph->n_coherent_ff <= 0) {
        fprintf(stderr, "No form factor data available\n");
        return -1;
    }

    alea_svg_plot_t plot;
    alea_svg_plot_init(&plot);
    plot.x_label = "Momentum transfer q (1/A)";
    plot.y_label = "Form factor";

    char title[128];
    snprintf(title, sizeof(title), "%s — Photon Form Factors", nuc->zaid);
    plot.title = title;

    if (ph->n_incoherent_ff > 0 && ph->incoherent_momentum && ph->incoherent_ff)
        alea_svg_plot_add(&plot, ph->incoherent_momentum, ph->incoherent_ff,
                     ph->n_incoherent_ff, "S(q,Z) incoherent", "#2266cc");

    if (ph->n_coherent_ff > 0 && ph->coherent_momentum && ph->coherent_ff)
        alea_svg_plot_add(&plot, ph->coherent_momentum, ph->coherent_ff,
                     ph->n_coherent_ff, "F(q,Z) coherent", "#cc2222");

    if (ph->n_coherent_ff > 0 && ph->coherent_momentum && ph->coherent_ff_cumulative)
        alea_svg_plot_add(&plot, ph->coherent_momentum, ph->coherent_ff_cumulative,
                     ph->n_coherent_ff, "Integrated F(q,Z)", "#22aa44");

    int rc = alea_svg_plot_write(&plot, output);
    if (rc == 0)
        printf("Written: %s\n", output);
    return rc;
}

/* ============================================================================
 * PLOT: Compton Scattered Photon Spectrum
 * ============================================================================ */

/**
 * Klein-Nishina Compton sampling — returns energy ratio x = E'/E.
 * Pure physics (no nuclear data needed).
 */
static double sample_compton_x(double alpha, uint64_t* rng_state) {
    double x_min = 1.0 / (1.0 + 2.0 * alpha);
    double a1 = log(1.0 + 2.0 * alpha);
    double a2 = 0.5 * (1.0 - x_min * x_min);

    for (int iter = 0; iter < 1000; iter++) {
        double r1 = local_rng(rng_state);
        double r2 = local_rng(rng_state);
        double r3 = local_rng(rng_state);

        double x;
        if (r1 * (a1 + a2) < a1)
            x = pow(x_min, 1.0 - r2);
        else
            x = sqrt(x_min * x_min + (1.0 - x_min * x_min) * r2);

        double t = (1.0 - x) / (alpha * x);
        double sin2 = t * (2.0 - t);
        double g = 1.0 - sin2 / (x + 1.0 / x);
        if (r3 <= g) return x;
    }
    return 0.5;
}

static int plot_photon_compton(alea_nuc_nuclide_t* nuc, double energy,
                                const char* output) {
    if (!nuc->photon) {
        fprintf(stderr, "Nuclide %s has no photon data\n", nuc->zaid);
        return -1;
    }

    double m_e = 0.511; /* MeV */
    double E_min = energy / (1.0 + 2.0 * energy / m_e);
    double E_max = energy;

    int N_bins = 200;
    double bin_width = (E_max - E_min) / N_bins;
    if (bin_width <= 0.0) {
        fprintf(stderr, "Energy too low for Compton spectrum\n");
        return -1;
    }

    double* energies = malloc((size_t)N_bins * sizeof(double));
    double* spectrum = calloc((size_t)N_bins, sizeof(double));
    if (!energies || !spectrum) {
        free(energies); free(spectrum);
        return -1;
    }

    for (int i = 0; i < N_bins; i++)
        energies[i] = E_min + (i + 0.5) * bin_width;

    /* Monte Carlo sampling via inline Klein-Nishina */
    int N_samples = 500000;
    uint64_t rng_state = 0xC0DE70BEEF123456ULL;
    double alpha = energy / m_e;

    for (int s = 0; s < N_samples; s++) {
        double x = sample_compton_x(alpha, &rng_state);
        double E_out = energy * x;

        if (E_out >= E_min && E_out <= E_max) {
            int b = (int)((E_out - E_min) / bin_width);
            if (b >= N_bins) b = N_bins - 1;
            if (b < 0) b = 0;
            spectrum[b] += 1.0;
        }
    }

    /* Normalize to PDF */
    double total = 0.0;
    for (int i = 0; i < N_bins; i++) total += spectrum[i] * bin_width;
    if (total > 0.0)
        for (int i = 0; i < N_bins; i++) spectrum[i] /= total;

    /* Klein-Nishina analytic reference */
    double* kn_ref = malloc((size_t)N_bins * sizeof(double));
    if (kn_ref) {
        double kn_sum = 0.0;
        double k = energy / m_e;
        for (int i = 0; i < N_bins; i++) {
            double r = energies[i] / energy;
            double cos_theta = 1.0 + 1.0 / k - 1.0 / (k * r);
            double sin2 = 1.0 - cos_theta * cos_theta;
            if (sin2 < 0.0) sin2 = 0.0;
            kn_ref[i] = (r + 1.0 / r - sin2) * r * r / (energy * k * k);
            kn_sum += kn_ref[i] * bin_width;
        }
        if (kn_sum > 0.0)
            for (int i = 0; i < N_bins; i++) kn_ref[i] /= kn_sum;
    }

    alea_svg_plot_t plot;
    alea_svg_plot_init(&plot);
    plot.x_scale = ALEA_SVG_SCALE_LINEAR;
    plot.y_scale = ALEA_SVG_SCALE_LINEAR;
    plot.x_label = "Scattered photon energy (MeV)";
    plot.y_label = "P(E') (1/MeV)";
    plot.x_min = E_min;
    plot.x_max = E_max;

    char title[128];
    snprintf(title, sizeof(title),
             "%s — Compton Spectrum at %.3g MeV", nuc->zaid, energy);
    plot.title = title;

    alea_svg_plot_add(&plot, energies, spectrum, N_bins, "Sampled", "#2266cc");
    if (kn_ref) {
        alea_svg_plot_add(&plot, energies, kn_ref, N_bins,
                     "Klein-Nishina (free e)", "#cc2222");
        plot.curves[1].dashed = 1;
    }

    int rc = alea_svg_plot_write(&plot, output);
    if (rc == 0)
        printf("Written: %s (%d samples)\n", output, N_samples);

    free(energies);
    free(spectrum);
    free(kn_ref);
    return rc;
}

/* ============================================================================
 * MAIN
 * ============================================================================ */

static void usage(const char* prog) {
    fprintf(stderr,
        "Usage: %s --xsdir <path> --zaid <zaid> [options]\n"
        "\n"
        "Options:\n"
        "  --xsdir <path>     Path to xsdir file or FENDL ace/ directory\n"
        "  --zaid <zaid>      Nuclide ZAID (e.g., 92235.32c)\n"
        "  --plot <type>      Plot type (default: xs)\n"
        "  --mt <list>        Comma-separated MT numbers (e.g., 18,102,16)\n"
        "  --list             List available reactions and exit\n"
        "  --output <file>    Output SVG file (default: <zaid>_<type>.svg)\n"
        "  --energy <MeV>     Energy for angular plot (default: 1.0)\n"
        "  --temp <K>         Temperature for Doppler plot (default: 1200)\n"
        "\n"
        "Plot types:\n"
        "  xs          Total/elastic/absorption cross sections\n"
        "              With --mt: plots only the specified reactions\n"
        "  reactions   Individual reaction cross sections\n"
        "              With --mt: filters to selected reactions\n"
        "  angular     Elastic angular distribution at given energy\n"
        "  fission     Fission spectrum chi(E)\n"
        "  nubar       Average neutrons per fission vs energy\n"
        "  doppler     Cross section at original vs broadened temperature\n"
        "  heating     Heating numbers vs energy (neutron and photon)\n"
        "\n"
        "Photon plot types (for .p tables):\n"
        "  photon-xs       Component cross sections vs energy\n"
        "  photon-ff       Form factors S(q,Z) and F(q,Z) vs momentum\n"
        "  photon-compton  Compton scattered energy spectrum at --energy\n"
        "\n"
        "Examples:\n"
        "  %s --xsdir ace/ --zaid 92235.32c --plot xs\n"
        "  %s --xsdir ace/ --zaid 92235.32c --plot reactions --mt 18,102,16\n"
        "  %s --xsdir ace/ --zaid 92235.32c --list\n"
        "  %s --xsdir mcplib84/ --zaid 82000.84p --plot photon-xs\n"
        "  %s --xsdir mcplib84/ --zaid 82000.84p --plot photon-compton -E 1.0\n",
        prog, prog, prog, prog, prog, prog);
}

int main(int argc, char** argv) {
    const char* xsdir_path = NULL;
    const char* zaid = NULL;
    const char* plot_type = "xs";
    const char* output = NULL;
    const char* mt_str = NULL;
    double energy = 1.0;
    double temp_K = 1200.0;
    int list_reactions = 0;

    /* Parse arguments */
    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "--xsdir") == 0 && i + 1 < argc)
            xsdir_path = argv[++i];
        else if (strcmp(argv[i], "--zaid") == 0 && i + 1 < argc)
            zaid = argv[++i];
        else if (strcmp(argv[i], "--plot") == 0 && i + 1 < argc)
            plot_type = argv[++i];
        else if (strcmp(argv[i], "--mt") == 0 && i + 1 < argc)
            mt_str = argv[++i];
        else if (strcmp(argv[i], "--list") == 0)
            list_reactions = 1;
        else if (strcmp(argv[i], "--output") == 0 && i + 1 < argc)
            output = argv[++i];
        else if ((strcmp(argv[i], "--energy") == 0 || strcmp(argv[i], "-E") == 0)
                 && i + 1 < argc)
            energy = atof(argv[++i]);
        else if ((strcmp(argv[i], "--temp") == 0 || strcmp(argv[i], "-T") == 0)
                 && i + 1 < argc)
            temp_K = atof(argv[++i]);
        else if (strcmp(argv[i], "--help") == 0 || strcmp(argv[i], "-h") == 0) {
            usage(argv[0]);
            return 0;
        } else {
            fprintf(stderr, "Unknown option: %s\n", argv[i]);
            usage(argv[0]);
            return 1;
        }
    }

    mt_filter_t mt_filter = {.count = 0};
    mt_filter_parse(&mt_filter, mt_str);

    if (!xsdir_path || !zaid) {
        usage(argv[0]);
        return 1;
    }

    /* Default output filename */
    char default_output[256];
    if (!output) {
        /* Sanitize ZAID for filename (replace . with _) */
        char safe_zaid[32];
        strncpy(safe_zaid, zaid, sizeof(safe_zaid) - 1);
        safe_zaid[sizeof(safe_zaid) - 1] = '\0';
        for (char* p = safe_zaid; *p; p++)
            if (*p == '.') *p = '_';
        snprintf(default_output, sizeof(default_output),
                 "%s_%s.svg", safe_zaid, plot_type);
        output = default_output;
    }

    /* Load nuclear data — try as directory first (FENDL), then as file (xsdir) */
    alea_nuc_xsdir_t* xsdir = alea_nuc_xsdir_load_dir(xsdir_path);
    if (!xsdir)
        xsdir = alea_nuc_xsdir_load(xsdir_path);
    if (!xsdir) {
        fprintf(stderr, "Failed to load xsdir from '%s'\n", xsdir_path);
        return 1;
    }

    printf("Loaded %zu xsdir entries\n", alea_nuc_xsdir_count(xsdir));

    alea_nuc_nuclide_t* nuc = alea_nuc_load_nuclide(xsdir, zaid);
    if (!nuc) {
        fprintf(stderr, "Failed to load nuclide '%s'\n", zaid);
        alea_nuc_xsdir_free(xsdir);
        return 1;
    }

    printf("Nuclide: %s (Z=%d A=%d, %d energies, %d reactions)\n",
           nuc->zaid, nuc->Z, nuc->A, nuc->n_energies, nuc->n_reactions);

    /* --list: print available reactions and exit */
    if (list_reactions) {
        if (nuc->particle == ALEA_NUC_PARTICLE_PHOTON && nuc->photon) {
            printf("Photon table %s (Z=%d):\n", nuc->zaid, nuc->Z);
            printf("  %d energy points\n", nuc->photon->n_energies);
            printf("  Interactions:\n");
            printf("    MT=504  Compton (incoherent) scattering\n");
            printf("    MT=502  Rayleigh (coherent) scattering\n");
            printf("    MT=522  Photoelectric absorption\n");
            printf("    MT=516  Pair production\n");
            printf("  Form factors:\n");
            printf("    S(q,Z) incoherent: %d points\n", nuc->photon->n_incoherent_ff);
            printf("    F(q,Z) coherent:   %d points\n", nuc->photon->n_coherent_ff);
        } else {
            print_reactions(nuc);
        }
        alea_nuc_nuclide_free(nuc);
        alea_nuc_xsdir_free(xsdir);
        return 0;
    }

    /* Generate plot */
    int rc;
    if (strcmp(plot_type, "xs") == 0)
        rc = plot_xs(nuc, &mt_filter, output);
    else if (strcmp(plot_type, "reactions") == 0)
        rc = plot_reactions(nuc, &mt_filter, output);
    else if (strcmp(plot_type, "angular") == 0)
        rc = plot_angular(nuc, energy, output);
    else if (strcmp(plot_type, "fission") == 0)
        rc = plot_fission(nuc, output);
    else if (strcmp(plot_type, "nubar") == 0)
        rc = plot_nubar(nuc, output);
    else if (strcmp(plot_type, "doppler") == 0)
        rc = plot_doppler(xsdir, zaid, temp_K, output);
    else if (strcmp(plot_type, "heating") == 0)
        rc = plot_heating(nuc, output);
    else if (strcmp(plot_type, "photon-xs") == 0)
        rc = plot_photon_xs(nuc, output);
    else if (strcmp(plot_type, "photon-ff") == 0)
        rc = plot_photon_ff(nuc, output);
    else if (strcmp(plot_type, "photon-compton") == 0)
        rc = plot_photon_compton(nuc, energy, output);
    else {
        fprintf(stderr, "Unknown plot type: %s\n", plot_type);
        rc = 1;
    }

    alea_nuc_nuclide_free(nuc);
    alea_nuc_xsdir_free(xsdir);
    return rc;
}
