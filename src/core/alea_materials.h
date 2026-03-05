// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

#ifndef ALEA_MATERIALS_H
#define ALEA_MATERIALS_H

/**
 * @file alea_materials.h
 * @brief Material, element, and mixture definitions for CSG system
 *
 * This module provides:
 * - Nuclide: Individual isotope with library suffix
 * - Element: Chemical element with natural isotopic composition
 * - Material: Composition of nuclides or elements with density
 * - Mixture: Combination of materials preserving component structure
 */

#include <stddef.h>
#include <stdbool.h>
#include "util/alea_vec.h"


/* ============================================================================
 * NUCLIDE DEFINITION
 * ============================================================================ */

/**
 * @brief Individual nuclide specification
 *
 * ZAID format: ZZAAA or ZZAAAI
 * - ZZ: atomic number (1-99)
 * - AAA: mass number (1-999)
 * - I: isomeric state (0=ground, 1,2,3=excited)
 *
 * Examples: 92235 (U-235), 1001 (H-1), 922351 (U-235m1)
 */
typedef struct {
    int zaid;               /* ZZAAA or ZZAAAI format */
    char* library;          /* Cross-section library suffix (e.g., ".80c") */
    double fraction;        /* Atom or weight fraction (always positive) */
} alea_nuclide_t;

/* ============================================================================
 * ELEMENT DEFINITION
 * ============================================================================ */

/**
 * @brief Natural isotope abundance entry
 */
typedef struct {
    int mass_number;        /* A (e.g., 235 for U-235) */
    double abundance;       /* Natural abundance (0-1, sum to 1) */
    double atomic_mass;     /* Atomic mass in amu */
} alea_isotope_t;

/**
 * @brief Chemical element with natural isotopic composition
 *
 * Elements store natural abundance data used to expand element-based
 * material definitions into explicit nuclide lists.
 */
typedef struct {
    int atomic_number;              /* Z (e.g., 92 for Uranium) */
    char symbol[4];                 /* Element symbol (e.g., "U", "Fe") */
    char name[32];                  /* Full name (e.g., "Uranium") */
    double standard_atomic_weight;  /* Standard atomic weight in amu */

    alea_isotope_t* isotopes;        /* Natural isotopes */
    size_t isotope_count;
} alea_element_t;

/**
 * @brief Element component in a material
 *
 * Used for element-based material composition. On export, elements
 * are expanded to their natural isotopic composition.
 */
typedef struct {
    int atomic_number;      /* Z of element */
    char* library;          /* Library suffix for generated nuclides */
    double fraction;        /* Weight or atom fraction */
} alea_element_comp_t;

/* ============================================================================
 * THERMAL SCATTERING
 * ============================================================================ */

/**
 * @brief Thermal scattering law reference (S(alpha,beta))
 *
 * References an S(a,b) table for thermal neutron scattering.
 */
typedef struct {
    char* identifier;       /* e.g., "lwtr.20t", "grph.20t" */
    int zaid_match;         /* Which ZAID this applies to (0 = match by Z) */
} alea_thermal_law_t;

/* Dynamic array types for material sub-arrays */
ALEA_VEC_DEFINE(alea_nuclide_vec, alea_nuclide_t);
ALEA_VEC_DEFINE(alea_element_comp_vec, alea_element_comp_t);
ALEA_VEC_DEFINE(alea_thermal_law_vec, alea_thermal_law_t);

/* ============================================================================
 * MATERIAL DEFINITION
 * ============================================================================ */

/**
 * @brief Material definition
 *
 * Materials can be defined either by:
 * 1. Explicit nuclides (nuclides array)
 * 2. Elements with natural composition (elements array)
 *
 * The standard_density provides a reference density used when cells
 * don't specify their own density.
 */
typedef struct alea_material {
    int material_id;                /* MCNP Mn number */

    /* Explicit nuclide composition */
    alea_nuclide_vec_t nuclides;

    /* Element-based composition (alternative) */
    alea_element_comp_vec_t elements;

    /* Fraction interpretation */
    bool is_weight_fraction;        /* true = weight, false = atom fractions */

    /* Standard density */
    double standard_density;        /* g/cm³ (>0) or atoms/b-cm (<0) */
    bool has_standard_density;

    /* Thermal scattering (MT card) */
    alea_thermal_law_vec_t thermal_laws;

    /* Metadata */
    char* name;                     /* User-friendly name */
    char* comments;                 /* Preserved comments */

    /* Computed properties (cached) */
    double avg_atomic_mass;
    bool properties_valid;
} alea_material_t;

/* ============================================================================
 * MIXTURE DEFINITION
 * ============================================================================ */

/**
 * @brief Component of a mixture
 */
typedef struct {
    int material_id;        /* Reference to a base material */
    double fraction;        /* Fraction in mixture (weight or atom) */
} alea_mixture_comp_t;

ALEA_VEC_DEFINE(alea_mixture_comp_vec, alea_mixture_comp_t);

/**
 * @brief Mixture of materials
 *
 * A mixture combines multiple materials. When exported to MCNP:
 * - Components are listed separately
 * - Identical nuclides are NOT combined
 * - Fractions are NOT renormalized
 *
 * This preserves the physical structure (e.g., a coating on a substrate).
 */
typedef struct {
    int mixture_id;                 /* Unique ID for this mixture */

    alea_mixture_comp_vec_t components;

    bool is_weight_fraction;        /* Interpretation of fractions */

    /* When exported, gets assigned an MCNP material number */
    int mcnp_material_id;

    /* Metadata */
    char* name;
    char* comments;
} alea_mixture_t;

/* ============================================================================
 * ELEMENT DATA
 * ============================================================================ */

/**
 * @brief Global element table (defined in alea_elements_data.c)
 * Index by atomic number: g_elements[Z] for Z=1..118
 * g_elements[0] is a placeholder
 */
extern const alea_element_t g_elements[119];
extern const size_t g_elements_count;

/**
 * @brief Get element data by atomic number
 * @param Z Atomic number (1-118)
 * @return Pointer to element data, or NULL if not found
 */
const alea_element_t* alea_get_element(int Z);

/**
 * @brief Get element data by symbol
 * @param symbol Element symbol (e.g., "H", "U", "Fe")
 * @return Pointer to element data, or NULL if not found
 */
const alea_element_t* alea_get_element_by_symbol(const char* symbol);

/**
 * @brief Get atomic number from ZAID
 */
static inline int alea_zaid_to_Z(int zaid) {
    if (zaid > 100000) zaid /= 10;  /* Remove isomer digit */
    return zaid / 1000;
}

/**
 * @brief Get mass number from ZAID
 */
static inline int alea_zaid_to_A(int zaid) {
    if (zaid > 100000) zaid /= 10;
    return zaid % 1000;
}

/**
 * @brief Create ZAID from Z and A
 */
static inline int alea_make_zaid(int Z, int A) {
    return Z * 1000 + A;
}

/* ============================================================================
 * MATERIAL OPERATIONS
 * ============================================================================ */

/**
 * @brief Create an empty material
 * @param material_id MCNP material number
 * @return New material (caller must free with alea_material_destroy)
 */
alea_material_t* alea_material_create(int material_id);

/**
 * @brief Destroy a material and free all memory
 */
void alea_material_destroy(alea_material_t* mat);

/**
 * @brief Add a nuclide to a material (internal, operates on material pointer)
 */
int alea_mat_add_nuclide(alea_material_t* mat, int zaid,
                         const char* library, double fraction);

/**
 * @brief Add an element to a material (internal, operates on material pointer)
 */
int alea_mat_add_element(alea_material_t* mat, int Z,
                         const char* library, double fraction);

/**
 * @brief Set material density (internal, operates on material pointer)
 */
void alea_mat_set_density(alea_material_t* mat, double density);

/**
 * @brief Expand element components to explicit nuclides (internal)
 */
int alea_mat_expand_elements(alea_material_t* mat);

/* ============================================================================
 * MIXTURE OPERATIONS
 * ============================================================================ */

/**
 * @brief Create an empty mixture
 * @param mixture_id Unique mixture ID
 * @return New mixture (caller must free with alea_mixture_destroy)
 */
alea_mixture_t* alea_mixture_create(int mixture_id);

/**
 * @brief Destroy a mixture
 */
void alea_mixture_destroy(alea_mixture_t* mix);

/**
 * @brief Add a material component to a mixture
 * @param mix Mixture to modify
 * @param material_id ID of material to add
 * @param fraction Fraction of this material in mixture
 * @return 0 on success, -1 on error
 */
int alea_mixture_add_component(alea_mixture_t* mix, int material_id,
                              double fraction);


#endif /* ALEA_MATERIALS_H */
