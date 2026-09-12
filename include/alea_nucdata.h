// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/**
 * @file alea_nucdata.h
 * @brief Nuclear data library — public API
 *
 * Reads ACE-format nuclear data files for use in Monte Carlo transport codes.
 * Supports continuous-energy neutron (.c) and photoatomic (.p) tables.
 *
 * All objects are user-owned. No hidden state or caches.
 *
 * Basic usage:
 *   alea_nuc_xsdir_t* xsdir = alea_nuc_xsdir_load("/path/to/xsdir");
 *   alea_nuc_nuclide_t* u235 = alea_nuc_load_nuclide(xsdir, "92235.80c");
 *   double sigma = alea_nuc_xs_total(u235, 1.0);  // 1 MeV
 *   alea_nuc_nuclide_free(u235);
 *   alea_nuc_xsdir_free(xsdir);
 */

#ifndef ALEA_NUCDATA_H
#define ALEA_NUCDATA_H

#include "alea_nucdata_types.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ============================================================================
 * XSDIR
 * ============================================================================ */

/**
 * @brief Load an xsdir/xsdata directory file
 *
 * Caller owns the returned xsdir and must free it with alea_nuc_xsdir_free().
 *
 * @return Allocated xsdir, or NULL on error
 */
alea_nuc_xsdir_t* alea_nuc_xsdir_load(const char* path);

/**
 * @brief Load all .xsd files from a directory (FENDL-style per-nuclide xsdir)
 *
 * Each .xsd file contains one xsdir-format entry line. The directory path
 * is used as the datapath for resolving relative ACE filenames.
 *
 * Caller owns the returned xsdir and must free it with alea_nuc_xsdir_free().
 *
 * @return Allocated xsdir, or NULL on error
 */
alea_nuc_xsdir_t* alea_nuc_xsdir_load_dir(const char* dirpath);

/**
 * @brief Free an xsdir and all its entries
 */
void alea_nuc_xsdir_free(alea_nuc_xsdir_t* xsdir);

/**
 * @brief Find an xsdir entry by ZAID
 * @return Pointer to entry (owned by xsdir), or NULL if not found
 */
const alea_nuc_xsdir_entry_t* alea_nuc_xsdir_find(const alea_nuc_xsdir_t* xsdir,
                                                   const char* zaid);

/**
 * @brief Find a table in the same ZAID family at an exact temperature
 *
 * The ZAID is used as a family anchor: the part before the dot and the table
 * type suffix must match, while the numeric evaluation suffix may differ.
 * For example, "1001.80c" can select "1001.81c", and "lwtr.20t" can
 * select another light-water thermal table. The closest entry whose kT is
 * within abs_tolerance is returned. Temperatures and tolerance are in MeV.
 *
 * On failure, *entry is set to NULL. Equidistant matching entries are
 * rejected because choosing between distinct evaluated tables is ambiguous.
 */
alea_error_t alea_nuc_xsdir_find_temperature(
    const alea_nuc_xsdir_t* xsdir,
    const char* zaid,
    double kT,
    double abs_tolerance,
    const alea_nuc_xsdir_entry_t** entry);

/**
 * @brief Find evaluated tables bracketing a requested temperature
 *
 * Matching uses the same ZAID-family and table-type rules as
 * alea_nuc_xsdir_find_temperature(). Extrapolation is rejected. For an exact
 * match, lower and upper identify the same entry and upper_fraction is zero.
 * Otherwise upper_fraction is (kT - lower->temperature) divided by the
 * bracketing temperature interval. Duplicate entries at either selected
 * temperature are rejected as ambiguous.
 */
alea_error_t alea_nuc_xsdir_find_temperature_bracket(
    const alea_nuc_xsdir_t* xsdir,
    const char* zaid,
    double kT,
    const alea_nuc_xsdir_entry_t** lower,
    const alea_nuc_xsdir_entry_t** upper,
    double* upper_fraction);

/**
 * @brief Get a nuclide from the xsdir cache, loading it if needed
 *
 * Unlike alea_nuc_load_nuclide() which always creates a fresh copy owned
 * by the caller, this function caches loaded nuclides in the xsdir.
 * Subsequent calls with the same ZAID return the cached pointer.
 *
 * The returned nuclide is owned by the xsdir and must NOT be freed or
 * modified by the caller. It remains valid until the xsdir is freed. Load a
 * fresh caller-owned nuclide with alea_nuc_load_nuclide() before applying an
 * in-place operation such as Doppler broadening.
 *
 * Uses an open-addressing hash table for O(1) amortized lookup,
 * suitable for large nuclide counts (burnup compositions, etc.).
 *
 * @return Borrowed pointer to nuclide, or NULL on error
 */
alea_nuc_nuclide_t* alea_nuc_xsdir_get_nuclide(alea_nuc_xsdir_t* xsdir, const char* zaid);

/**
 * @brief Get number of entries in loaded xsdir
 */
size_t alea_nuc_xsdir_count(const alea_nuc_xsdir_t* xsdir);

/* ============================================================================
 * ACE TABLE LOADING
 * ============================================================================ */

/**
 * @brief Load raw ACE table from file
 *
 * Reads header (NXS, JXS) and XSS data array. Does not decode
 * physics blocks — use alea_nuc_load_nuclide() for that.
 *
 * @param path      Path to ACE file
 * @param address   Start line (1-based, Type 1) or record number (Type 2)
 * @param file_type 1 = ASCII, 2 = binary
 * @param table     Output table (caller manages lifetime)
 * @return ALEA_OK on success
 */
alea_error_t alea_nuc_ace_read(const char* path, int address, int file_type,
                          alea_nuc_ace_table_t* table);

/**
 * @brief Free internals of an ACE table (does not free the struct itself)
 */
void alea_nuc_ace_free(alea_nuc_ace_table_t* table);

/* ============================================================================
 * NUCLIDE LOADING
 * ============================================================================ */

/**
 * @brief Load a nuclide by ZAID, decoding cross sections from ACE data
 *
 * Looks up the ZAID in the xsdir, reads the ACE file, and decodes all
 * physics blocks. Caller owns the returned nuclide and must free it
 * with alea_nuc_nuclide_free().
 *
 * @return Allocated nuclide, or NULL on error
 */
alea_nuc_nuclide_t* alea_nuc_load_nuclide(const alea_nuc_xsdir_t* xsdir, const char* zaid);

/**
 * @brief Free a nuclide and all its data
 */
void alea_nuc_nuclide_free(alea_nuc_nuclide_t* nuc);

/**
 * Load and decode a thermal scattering ACE table. Discrete IFENG=0/1 and
 * continuous correlated IFENG=2 inelastic representations are supported.
 * The caller owns the result.
 */
alea_nuc_thermal_t* alea_nuc_load_thermal(const alea_nuc_xsdir_t* xsdir,
                                          const char* zaid);

void alea_nuc_thermal_free(alea_nuc_thermal_t* thermal);

/** Thermal incoherent-inelastic, elastic, and total cross sections (barns). */
double alea_nuc_thermal_xs_inelastic(const alea_nuc_thermal_t* thermal,
                                     double energy);
double alea_nuc_thermal_xs_elastic(const alea_nuc_thermal_t* thermal,
                                   double energy);
double alea_nuc_thermal_xs_total(const alea_nuc_thermal_t* thermal,
                                 double energy);

/**
 * Sample one collision from a decoded thermal table. The outgoing
 * neutron retains the incident weight and time. MT is 2 for elastic and 4 for
 * incoherent inelastic scattering.
 */
alea_error_t alea_nuc_sample_thermal_collision(
    const alea_nuc_thermal_t* thermal,
    const alea_nuc_particle_state_t* incident,
    alea_nuc_random_fn random,
    void* random_context,
    alea_nuc_collision_result_t* result);

/* ============================================================================
 * CROSS-SECTION LOOKUP (microscopic)
 *
 * All energies in MeV. Returns cross section in barns.
 * Uses binary search + interpolation on the ACE energy grid.
 * ============================================================================ */

double alea_nuc_xs_total(const alea_nuc_nuclide_t* nuc, double energy);
double alea_nuc_xs_absorption(const alea_nuc_nuclide_t* nuc, double energy);
double alea_nuc_xs_elastic(const alea_nuc_nuclide_t* nuc, double energy);
double alea_nuc_xs_reaction(const alea_nuc_nuclide_t* nuc, int mt, double energy);

/** Heating cross section at given energy (MeV·barn). Works for neutron and photon. */
double alea_nuc_xs_heating(const alea_nuc_nuclide_t* nuc, double energy);

/** Average energy deposited per collision (MeV) = heating / σ_total */
double alea_nuc_heating_per_collision(const alea_nuc_nuclide_t* nuc, double energy);

/**
 * @brief Find energy grid index for interpolation
 *
 * Returns i such that energy[i] <= E < energy[i+1].
 * Also computes interpolation fraction f = (E - E[i]) / (E[i+1] - E[i]).
 * Values outside the grid are clamped to the nearest endpoint interval.
 *
 * @param energy    Energy grid (ascending)
 * @param n         Grid size
 * @param E         Query energy (MeV)
 * @param frac      Output interpolation fraction [0,1]
 * @return Grid index, or -1 if the grid/input is invalid
 */
int alea_nuc_energy_lookup(const double* energy, int n, double E, double* frac);

/** Log-log interpolation on a grid */
double alea_nuc_interp_loglog(const double* grid, const double* values, int n,
                          double x);

/**
 * @brief Get URR-modified cross section factors at given energy
 *
 * Samples a probability band using the inverse CDF and returns cross section
 * factors for total, elastic, fission, capture, and heating. Absolute ACE
 * probability tables are normalized to the corresponding smooth cross
 * sections before these factors are returned. The ACE total is diagnostic;
 * transport reconstructs its total from sampled partials and the smooth
 * competition channels selected by the table's ILF and IOA flags.
 *
 * @param nuc     Nuclide with URR data
 * @param energy  Incident energy (MeV)
 * @param xi      Random number [0,1) for band selection
 * @param factors Output array of 5 factors [total, elastic, fission, capture, heating]
 * @return 1 if URR applies at this energy, 0 otherwise
 */
int alea_nuc_urr_factors(const alea_nuc_nuclide_t* nuc, double energy, double xi,
                     double factors[5]);

/* ============================================================================
 * MATERIAL
 * ============================================================================ */

alea_nuc_material_t* alea_nuc_material_create(void);
void            alea_nuc_material_destroy(alea_nuc_material_t* mat);

alea_error_t alea_nuc_material_add(alea_nuc_material_t* mat, alea_nuc_nuclide_t* nuclide,
                              double number_density);

/**
 * Add one nuclide represented by two bracketing temperature tables. The
 * material borrows both tables and adds densities (1-upper_fraction)*N and
 * upper_fraction*N atomically. This gives linear expected cross sections;
 * collision selection retains the distributions and URR data belonging to
 * the selected table. Coordinated URR evaluation uses the same probability
 * quantile for both bounding tables. Endpoint fractions add only one component.
 */
alea_error_t alea_nuc_material_add_temperature_mix(
    alea_nuc_material_t* mat,
    alea_nuc_nuclide_t* lower,
    alea_nuc_nuclide_t* upper,
    double upper_fraction,
    double number_density);

/**
 * Build a nuclear material from a zero-based geometry cell index.
 *
 * The caller owns the returned material. Nuclides referenced by it remain
 * owned by xsdir, so xsdir must outlive the returned material.
 */
struct alea_system;
alea_nuc_material_t* alea_nuc_material_from_cell(
    struct alea_system* sys, int cell_index, alea_nuc_xsdir_t* xsdir);

/** Macroscopic cross sections (cm⁻¹) */
double alea_nuc_mat_xs_total(const alea_nuc_material_t* mat, double energy);
double alea_nuc_mat_xs_absorption(const alea_nuc_material_t* mat, double energy);
double alea_nuc_mat_xs_elastic(const alea_nuc_material_t* mat, double energy);

/** Mean free path (cm) */
double alea_nuc_mean_free_path(const alea_nuc_material_t* mat, double energy);

/** Sample distance to next collision: -ln(1-ξ)/Σ_t */
double alea_nuc_sample_distance(const alea_nuc_material_t* mat, double energy, double xi);

/* ============================================================================
 * COLLISION SAMPLING
 * ============================================================================ */

/**
 * @brief Sample which nuclide in a material is hit
 * @param xi  Random number [0,1)
 * @param out_nuclide  Output: selected nuclide
 * @return Component index, or -1 on error
 */
int alea_nuc_sample_nuclide(const alea_nuc_material_t* mat, double energy, double xi,
                        alea_nuc_nuclide_t** out_nuclide);

/**
 * @brief Sample which reaction (MT) occurs on a given nuclide
 * @param xi  Random number [0,1)
 * @param out_mt  Output: selected MT number
 * @return Reaction index (0=elastic), or -1 on error
 */
int alea_nuc_sample_reaction(const alea_nuc_nuclide_t* nuc, double energy, double xi,
                         int* out_mt);

/**
 * Sample a single-nuclide collision using caller-supplied random values.
 * The initial implementation supports stationary-target elastic scattering
 * (MT=2) and absorption reactions. Energy is in MeV and xi values must lie in
 * [0,1). For elastic scattering, result->mu is in the center-of-mass frame.
 * The result is unchanged on failure.
 */
alea_error_t alea_nuc_sample_collision(const alea_nuc_nuclide_t* nuc, int mt,
                                       double energy, const double xi[3],
                                       alea_nuc_interaction_t* result);

/* ============================================================================
 * PREPARED CONTINUOUS-ENERGY COLLISION PHYSICS
 * ============================================================================ */

/** Initialize a Philox event stream from src/rng. Particle ordinals must be
 * unique within a history; event indices must be unique within a particle.
 * Separate domains isolate flight, collision, and URR draws. No allocation.
 */
alea_error_t alea_nuc_rng_init(alea_nuc_rng_t* rng, uint64_t seed,
    uint32_t history_id, uint32_t particle_ordinal, uint32_t event_index,
    alea_nuc_rng_domain_t domain);

/** Callback for nuclear-data samplers: 53-bit uniform in [0,1).
 * Returns NAN on invalid state or exhaustion of the event's 2^32 words.
 */
double alea_nuc_rng_uniform(void* context);

/**
 * Sample an outgoing energy from an ACE energy-distribution chain.
 * Supported laws are level scattering (3), continuous tabular (4), general
 * evaporation (5), Maxwell (7), evaporation (9), Watt (11), and the energy
 * marginal of Kalbach-Mann (44), correlated law 61, N-body phase space law 66,
 * and laboratory angle-energy law 67. The RNG must return finite values in
 * [0,1). The output is unchanged on failure.
 */
alea_error_t alea_nuc_sample_energy_distribution(
    const alea_nuc_energy_dist_t* distribution,
    double incident_energy,
    alea_nuc_random_fn random,
    void* random_context,
    double* energy_out);

/**
 * Sample an outgoing energy and any angle correlated with it. For laws 44,
 * 61, 66, and 67, mu_out is the sampled cosine and angle_is_correlated is
 * true. Other supported laws leave mu_out unchanged and report false so the
 * caller can sample the reaction's separate angular distribution.
 */
alea_error_t alea_nuc_sample_energy_angle_distribution(
    const alea_nuc_energy_dist_t* distribution,
    double incident_energy,
    alea_nuc_random_fn random,
    void* random_context,
    double* energy_out,
    double* mu_out,
    bool* angle_is_correlated);

/**
 * Sample elastic scattering from a free target at kT (MeV), using the
 * collision-conditioned constant-cross-section target-velocity model.
 */
alea_error_t alea_nuc_sample_free_gas_elastic(
    const alea_nuc_nuclide_t* nuc,
    const alea_nuc_particle_state_t* incident,
    double kT,
    alea_nuc_random_fn random,
    void* random_context,
    alea_nuc_free_gas_result_t* result);

/** Inspect whether a nuclide supports the restricted collision model. */
alea_error_t alea_nuc_capabilities(const alea_nuc_nuclide_t* nuc,
                                   alea_nuc_capability_report_t* report);

/**
 * Prepare an immutable material for collision sampling. The returned object
 * borrows the material, nuclides, and any thermal tables named by the
 * requirements; all must outlive it. A thermal association replaces the
 * component's free-atom elastic cross section through the highest incident
 * energy in the thermal table. Preparation fails when active reaction physics
 * cannot be represented faithfully.
 */
alea_error_t alea_nuc_prepare_material(
    const alea_nuc_material_t* material,
    const alea_nuc_prepare_requirements_t* requirements,
    alea_nuc_capability_report_t* report,
    alea_nuc_prepared_material_t** prepared);

void alea_nuc_prepared_material_free(alea_nuc_prepared_material_t* prepared);

/** Evaluate one incident neutron in a prepared material without sampling. */
alea_error_t alea_nuc_evaluate(
    const alea_nuc_prepared_material_t* prepared,
    const alea_nuc_particle_state_t* incident,
    alea_nuc_evaluation_t* evaluation);

/**
 * Evaluate a neutron while sampling one coordinated unresolved-resonance
 * realization per material component. The caller-owned workspace must remain
 * unchanged and alive through flight and collision sampling. Components made
 * by one temperature-mix call share a probability quantile.
 */
alea_error_t alea_nuc_evaluate_urr(
    const alea_nuc_prepared_material_t* prepared,
    const alea_nuc_particle_state_t* incident,
    alea_nuc_random_fn random,
    void* random_context,
    alea_nuc_evaluation_workspace_t* workspace,
    alea_nuc_evaluation_t* evaluation);

/**
 * Sample a flight distance in cm from an existing evaluation. The RNG must
 * return finite values in [0,1). The distance is unchanged on failure. A call
 * may consume RNG values before reporting an error.
 */
alea_error_t alea_nuc_sample_flight(const alea_nuc_evaluation_t* evaluation,
                                    alea_nuc_random_fn random,
                                    void* random_context,
                                    double* distance);

/**
 * Sample target, reaction and collision outcome from one evaluation. The RNG
 * must return finite values in [0,1). No allocation occurs in this call. The
 * result is unchanged on failure; caller RNG values may already be consumed.
 */
alea_error_t alea_nuc_collide(const alea_nuc_evaluation_t* evaluation,
                              alea_nuc_random_fn random,
                              void* random_context,
                              alea_nuc_collision_result_t* result);

/**
 * Sample a collision that may emit neutrons and, when requested during
 * preparation, reaction-conditioned photons. Emissions are appended to
 * caller-owned storage. Capacity is
 * checked before RNG is consumed, and buffer count/result are unchanged on
 * failure. Ordinary calls allocate no memory.
 */
alea_error_t alea_nuc_collide_with_secondaries(
    const alea_nuc_evaluation_t* evaluation,
    alea_nuc_random_fn random,
    void* random_context,
    alea_nuc_secondary_buffer_t* secondaries,
    alea_nuc_collision_result_t* result);

/* ============================================================================
 * FISSION DATA
 * ============================================================================ */

/** Evaluate ν̄(E) — average neutrons per fission */
double alea_nuc_nu_bar(const alea_nuc_nuclide_t* nuc, double energy);

/** Evaluate prompt ν̄(E), falling back to the sole ν̄ representation. */
double alea_nuc_prompt_nu_bar(const alea_nuc_nuclide_t* nuc, double energy);

/** Evaluate delayed ν̄(E), or zero when no delayed-yield data are present. */
double alea_nuc_delayed_nu_bar(const alea_nuc_nuclide_t* nuc, double energy);

/* ============================================================================
 * REACTION CLASSIFICATION
 * ============================================================================ */

/** Classify an MT reaction number into scatter/multiply/absorption */
alea_nuc_reaction_class_t alea_nuc_reaction_classify(int mt);

/** Get neutron yield for a reaction at given energy (from TYR field) */
double alea_nuc_reaction_yield(const alea_nuc_nuclide_t* nuc, int mt, double energy);

/** Evaluate the mean photon yield for one decoded production channel. */
double alea_nuc_photon_production_yield(
    const alea_nuc_nuclide_t* nuc,
    const alea_nuc_photon_production_t* production,
    double energy);

/** Evaluate the aggregate neutron-induced photon-production cross section. */
double alea_nuc_xs_photon_production_total(
    const alea_nuc_nuclide_t* nuc, double energy);

/** Sample one photon from a decoded production channel without allocation. */
alea_error_t alea_nuc_sample_photon_production(
    const alea_nuc_nuclide_t* nuc,
    const alea_nuc_photon_production_t* production,
    const alea_nuc_particle_state_t* incident,
    alea_nuc_random_fn random,
    void* random_context,
    alea_nuc_particle_state_t* photon);

/* ============================================================================
 * PHOTON CROSS SECTIONS
 * ============================================================================ */

/** Photon cross section by component (barns). Log-log interpolation. */
double alea_nuc_photon_xs_incoherent(const alea_nuc_nuclide_t* nuc, double energy);
double alea_nuc_photon_xs_coherent(const alea_nuc_nuclide_t* nuc, double energy);
double alea_nuc_photon_xs_photoelectric(const alea_nuc_nuclide_t* nuc, double energy);
/** Photoelectric cross section for one EPR subshell designator (barns). */
double alea_nuc_photon_xs_photoelectric_subshell(
    const alea_nuc_nuclide_t* nuc, int designator, double energy);
double alea_nuc_photon_xs_pair(const alea_nuc_nuclide_t* nuc, double energy);

/** Maximum secondary slots required by any photoatomic event at this energy. */
size_t alea_nuc_photon_secondary_capacity(
    const alea_nuc_nuclide_t* element, double energy);

/**
 * Sample one photoatomic collision. MT 502 is coherent scattering, MT 504
 * incoherent scattering, MT 517 pair production, and MT 522 photoelectric
 * absorption. Charged-particle and annihilation energy is deposited locally.
 */
alea_error_t alea_nuc_sample_photon_collision(
    const alea_nuc_nuclide_t* element,
    const alea_nuc_particle_state_t* incident,
    alea_nuc_random_fn random,
    void* random_context,
    alea_nuc_collision_result_t* result);

/**
 * Sample one photoatomic collision and append explicit secondaries. Pair
 * production emits two back-to-back 0.51099895069 MeV annihilation photons.
 * Photoelectric absorption emits the radiative part of a detailed EPR atomic
 * relaxation cascade when those data are present. Charged-particle energy is
 * deposited locally. Capacity for every possible secondary is required before
 * RNG consumption.
 */
alea_error_t alea_nuc_sample_photon_collision_with_secondaries(
    const alea_nuc_nuclide_t* element,
    const alea_nuc_particle_state_t* incident,
    alea_nuc_random_fn random,
    void* random_context,
    alea_nuc_secondary_buffer_t* secondaries,
    alea_nuc_collision_result_t* result);

/* ============================================================================
 * MULTIGROUP
 * ============================================================================ */

/**
 * @brief Create a multigroup data structure
 *
 * Group convention: bounds[0] = highest energy, bounds[G] = lowest.
 * Boundaries are in descending order.
 *
 * @param n_groups  Number of energy groups
 * @param bounds    Group boundaries [n_groups+1], descending (MeV)
 * @return Allocated structure, or NULL on error
 */
alea_nuc_multigroup_t* alea_nuc_mg_create(int n_groups, const double* bounds);
void alea_nuc_mg_destroy(alea_nuc_multigroup_t* mg);

/**
 * @brief Set user-defined weighting spectrum for multigroup collapse
 *
 * If not set, the default spectrum is used: Maxwellian (E < 0.625 eV),
 * 1/E (0.625 eV to 100 keV), and fission spectrum (> 100 keV).
 *
 * @param mg       Multigroup structure
 * @param fn       Spectrum function φ(E, ctx), or NULL to reset to default
 * @param ctx      User context pointer passed to fn
 */
void alea_nuc_mg_set_spectrum(alea_nuc_multigroup_t* mg, alea_nuc_spectrum_fn fn,
                              void* ctx);

/**
 * @brief Collapse continuous-energy cross sections into multigroup constants
 *
 * Computes group-averaged total, absorption, elastic, fission cross sections,
 * fission spectrum χ, and scattering transfer matrix (elastic + inelastic).
 * Uses the weighting spectrum set via alea_nuc_mg_set_spectrum(), or the
 * default Maxwellian+1/E+fission spectrum if none was set.
 *
 * @param mg   Multigroup structure (must have bounds set)
 * @param nuc  Nuclide with pointwise cross sections
 * @return ALEA_OK on success
 */
alea_error_t alea_nuc_mg_collapse(alea_nuc_multigroup_t* mg, const alea_nuc_nuclide_t* nuc);

/** Forward scattering matrix element: σ_s(g_from → g_to) */
double alea_nuc_mg_scatter(const alea_nuc_multigroup_t* mg, int g_from, int g_to);

/** Adjoint scattering matrix element: transpose of forward */
double alea_nuc_mg_scatter_adjoint(const alea_nuc_multigroup_t* mg, int g_from, int g_to);

/**
 * @brief Sample outgoing group from scattering
 *
 * @param mg       Multigroup data
 * @param g_from   Incoming group
 * @param xi       Random number [0,1)
 * @param adjoint  0 = forward, 1 = adjoint (transposed matrix)
 * @return Outgoing group index
 */
int alea_nuc_mg_sample_scatter(const alea_nuc_multigroup_t* mg, int g_from,
                           double xi, int adjoint);

/* ============================================================================
 * DOPPLER BROADENING
 * ============================================================================ */

/**
 * @brief Doppler-broaden cross sections to a higher temperature
 *
 * Convolves pointwise cross sections with the exact Doppler kernel
 * to account for thermal motion of target nuclei. Modifies the nuclide
 * in-place (total, absorption, elastic, heating, and per-reaction XS).
 *
 * Can only broaden to a higher temperature than the current one. Nuclides
 * with unresolved-resonance probability tables are rejected because those
 * tables retain the temperature at which they were processed.
 * The nuclide must be caller-owned; cached pointers returned by
 * alea_nuc_xsdir_get_nuclide() are shared and immutable. The decoded URR
 * probability table and raw ACE data retain their source-table temperature.
 *
 * @param nuc        Nuclide to broaden
 * @param kT_target  Target temperature in MeV (e.g., 2.53e-8 for 293.6 K)
 * @return ALEA_OK on success, ALEA_ERR_UNSUPPORTED if kT_target <= current
 *         or a URR probability table is attached
 */
alea_error_t alea_nuc_doppler_broaden(alea_nuc_nuclide_t* nuc, double kT_target);

/* ============================================================================
 * UTILITY
 * ============================================================================ */

/**
 * @brief Parse ZAID string into Z, A, metastable, and table type
 * @param zaid      e.g. "92235.80c"
 * @param Z         Output atomic number
 * @param A         Output mass number
 * @param meta      Output metastable state (0 = ground)
 * @param type      Output table type
 * @return ALEA_OK on success
 */
alea_error_t alea_nuc_parse_zaid(const char* zaid, int* Z, int* A, int* meta,
                            alea_nuc_table_type_t* type);

/* Error strings: use alea_error_string() from alea_types.h */

#ifdef __cplusplus
}
#endif

#endif /* ALEA_NUCDATA_H */
