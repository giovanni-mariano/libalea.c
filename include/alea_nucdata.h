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
 * @brief Get a nuclide from the xsdir cache, loading it if needed
 *
 * Unlike alea_nuc_load_nuclide() which always creates a fresh copy owned
 * by the caller, this function caches loaded nuclides in the xsdir.
 * Subsequent calls with the same ZAID return the cached pointer.
 *
 * The returned nuclide is owned by the xsdir and must NOT be freed
 * by the caller. It remains valid until the xsdir is freed.
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
 * @param address   Start line (1-based, Type 1) or byte offset (Type 2)
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

/** Heating number at given energy (MeV·barn). Works for neutron and photon. */
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
 * Samples a probability band using xi, returns cross section factors
 * for total, elastic, fission, capture, heating.
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
 * @brief Create a nucdata material from a core material definition
 *
 * Loads ACE nuclides via the xsdir and computes number densities from
 * the material's composition (atom or weight fractions) and density.
 * Elements must be expanded to nuclides before calling this function
 * (see alea_mat_expand_elements).
 *
 * The returned material owns the loaded nuclides — they are freed
 * when the material is destroyed.
 *
 * @param mat           Core material (with nuclides and fraction type)
 * @param cell_density  Cell density (positive = g/cm³, negative = atoms/b-cm)
 * @param is_mass_density  true if density is in g/cm³
 * @param xsdir         xsdir for resolving ZAIDs and loading ACE data
 * @return Allocated nucdata material, or NULL on error
 */
struct alea_material;
alea_nuc_material_t* alea_nuc_material_from_core(
    const struct alea_material* mat,
    double cell_density,
    bool is_mass_density,
    alea_nuc_xsdir_t* xsdir);

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

/* ============================================================================
 * FISSION DATA
 * ============================================================================ */

/** Evaluate ν̄(E) — average neutrons per fission */
double alea_nuc_nu_bar(const alea_nuc_nuclide_t* nuc, double energy);

/* ============================================================================
 * REACTION CLASSIFICATION
 * ============================================================================ */

/** Classify an MT reaction number into scatter/multiply/absorption */
alea_nuc_reaction_class_t alea_nuc_reaction_classify(int mt);

/** Get neutron yield for a reaction at given energy (from TYR field) */
double alea_nuc_reaction_yield(const alea_nuc_nuclide_t* nuc, int mt, double energy);

/* ============================================================================
 * PHOTON CROSS SECTIONS
 * ============================================================================ */

/** Photon cross section by component (barns). Log-log interpolation. */
double alea_nuc_photon_xs_incoherent(const alea_nuc_nuclide_t* nuc, double energy);
double alea_nuc_photon_xs_coherent(const alea_nuc_nuclide_t* nuc, double energy);
double alea_nuc_photon_xs_photoelectric(const alea_nuc_nuclide_t* nuc, double energy);
double alea_nuc_photon_xs_pair(const alea_nuc_nuclide_t* nuc, double energy);

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
 * Can only broaden to a higher temperature than the current one.
 *
 * @param nuc        Nuclide to broaden
 * @param kT_target  Target temperature in MeV (e.g., 2.53e-8 for 293.6 K)
 * @return ALEA_OK on success, ALEA_ERR_INVALID_ARG if kT_target <= current
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
