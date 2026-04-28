// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/**
 * @file alea_materials.c
 * @brief Material, element, and mixture implementation
 *
 * Element data is in alea_elements_data.c
 */

#include "alea_materials.h"
#include "util/compat.h"
#include <stdlib.h>
#include <string.h>

/* ============================================================================
 * ELEMENT LOOKUP
 * ============================================================================ */

const alea_element_t* alea_get_element(int Z) {
    if (Z < 1 || Z > 118) return NULL;
    return &g_elements[Z];
}

const alea_element_t* alea_get_element_by_symbol(const char* symbol) {
    if (!symbol) return NULL;

    for (size_t i = 1; i < g_elements_count; i++) {
        if (alea_strcasecmp(g_elements[i].symbol, symbol) == 0) {
            return &g_elements[i];
        }
    }
    return NULL;
}

/* ============================================================================
 * MATERIAL OPERATIONS
 * ============================================================================ */

alea_material_t* alea_material_create(int material_id) {
    alea_material_t* mat = calloc(1, sizeof(alea_material_t));
    if (!mat) return NULL;

    mat->material_id = material_id;
    int r = alea_vec_reserve(&mat->nuclides, 8, alea_nuclide_t);
    if (r != 0) {
        free(mat);
        return NULL;
    }

    r = alea_vec_reserve(&mat->elements, 4, alea_element_comp_t);
    if (r != 0) {
        alea_vec_free(&mat->nuclides);
        free(mat);
        return NULL;
    }

    return mat;
}

void alea_material_destroy(alea_material_t* mat) {
    if (!mat) return;

    for (size_t i = 0; i < alea_vec_count(&mat->nuclides); i++) {
        free(mat->nuclides.data[i].library);
    }
    alea_vec_free(&mat->nuclides);

    for (size_t i = 0; i < alea_vec_count(&mat->elements); i++) {
        free(mat->elements.data[i].library);
    }
    alea_vec_free(&mat->elements);

    for (size_t i = 0; i < alea_vec_count(&mat->thermal_laws); i++) {
        free(mat->thermal_laws.data[i].identifier);
    }
    alea_vec_free(&mat->thermal_laws);

    free(mat->name);
    free(mat->comments);
    free(mat);
}

int alea_mat_add_nuclide(alea_material_t* mat, int zaid,
                             const char* library, double fraction) {
    if (!mat) return -1;

    alea_nuclide_t* nuc = alea_vec_push_uninit(&mat->nuclides, alea_nuclide_t);
    if (!nuc) return -1;

    nuc->zaid = zaid;
    nuc->library = library ? alea_strdup(library) : NULL;
    nuc->fraction = fraction;

    mat->properties_valid = false;
    return 0;
}

int alea_mat_add_element(alea_material_t* mat, int Z,
                             const char* library, double fraction) {
    if (!mat) return -1;

    /* Verify element exists */
    if (!alea_get_element(Z)) return -1;

    alea_element_comp_t* elem = alea_vec_push_uninit(&mat->elements, alea_element_comp_t);
    if (!elem) return -1;

    elem->atomic_number = Z;
    elem->library = library ? alea_strdup(library) : NULL;
    elem->fraction = fraction;

    mat->properties_valid = false;
    return 0;
}

void alea_mat_set_density(alea_material_t* mat, double density) {
    if (!mat) return;
    mat->standard_density = density;
    mat->has_standard_density = true;
}

int alea_mat_expand_elements(alea_material_t* mat) {
    if (!mat || alea_vec_count(&mat->elements) == 0) return 0;

    /* For each element, add its isotopes as nuclides */
    for (size_t i = 0; i < alea_vec_count(&mat->elements); i++) {
        const alea_element_comp_t* ec = &mat->elements.data[i];
        const alea_element_t* elem = alea_get_element(ec->atomic_number);

        if (!elem || !elem->isotopes || elem->isotope_count == 0) {
            /* No isotope data - add as single nuclide with A=0 (natural) */
            int zaid = ec->atomic_number * 1000;  /* ZZ000 format for natural */
            if (alea_mat_add_nuclide(mat, zaid, ec->library,
                                         ec->fraction) < 0) {
                return -1;
            }
            continue;
        }

        /* Add each isotope with scaled fraction */
        for (size_t j = 0; j < elem->isotope_count; j++) {
            const alea_isotope_t* iso = &elem->isotopes[j];

            /* Skip isotopes with zero abundance (synthetic) */
            if (iso->abundance <= 0.0) continue;

            int zaid = alea_make_zaid(ec->atomic_number, iso->mass_number);
            double frac = ec->fraction * iso->abundance;

            if (alea_mat_add_nuclide(mat, zaid, ec->library, frac) < 0) {
                return -1;
            }
        }
    }

    /* Clear element list (already expanded) */
    for (size_t i = 0; i < alea_vec_count(&mat->elements); i++) {
        free(mat->elements.data[i].library);
    }
    alea_vec_clear(&mat->elements);

    return 0;
}

/* ============================================================================
 * MIXTURE OPERATIONS
 * ============================================================================ */

alea_mixture_t* alea_mixture_create(int mixture_id) {
    alea_mixture_t* mix = calloc(1, sizeof(alea_mixture_t));
    if (!mix) return NULL;

    mix->mixture_id = mixture_id;
    int r = alea_vec_reserve(&mix->components, 4, alea_mixture_comp_t);
    if (r != 0) {
        free(mix);
        return NULL;
    }

    return mix;
}

void alea_mixture_destroy(alea_mixture_t* mix) {
    if (!mix) return;
    alea_vec_free(&mix->components);
    free(mix->name);
    free(mix->comments);
    free(mix);
}

int alea_mixture_add_component(alea_mixture_t* mix, int material_id,
                              double fraction) {
    if (!mix) return -1;

    alea_mixture_comp_t* comp = alea_vec_push_uninit(&mix->components, alea_mixture_comp_t);
    if (!comp) return -1;

    comp->material_id = material_id;
    comp->fraction = fraction;

    return 0;
}
