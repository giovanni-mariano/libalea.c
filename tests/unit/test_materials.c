// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/*
 * test_materials.c - Tests for public material and mixture API
 */

#include "alea_test.h"
#include "alea.h"

/* ========================================================================= */
/* Material lifecycle                                                        */
/* ========================================================================= */

TEST(material_add_and_count) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);

    ASSERT_EQ(alea_material_count(sys), 0);

    int m0 = alea_add_material(sys, 1);
    ASSERT(m0 >= 0);
    ASSERT_EQ(alea_material_count(sys), 1);

    int m1 = alea_add_material(sys, 5);
    ASSERT(m1 >= 0);
    ASSERT_EQ(alea_material_count(sys), 2);

    ASSERT_EQ(alea_material_get_id(sys, m0), 1);
    ASSERT_EQ(alea_material_get_id(sys, m1), 5);

    alea_destroy(sys);
}

TEST(material_find_by_id) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);

    alea_add_material(sys, 10);
    alea_add_material(sys, 20);

    ASSERT_EQ(alea_find_material_by_id(sys, 10), 0);
    ASSERT_EQ(alea_find_material_by_id(sys, 20), 1);
    ASSERT_EQ(alea_find_material_by_id(sys, 99), -1);

    alea_destroy(sys);
}

TEST(material_auto_id) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);

    int m0 = alea_add_material(sys, 0);  /* auto-assign */
    int m1 = alea_add_material(sys, 0);  /* auto-assign */
    ASSERT(m0 >= 0);
    ASSERT(m1 >= 0);

    int id0 = alea_material_get_id(sys, m0);
    int id1 = alea_material_get_id(sys, m1);
    ASSERT(id0 > 0);
    ASSERT(id1 > 0);
    ASSERT(id0 != id1);

    alea_destroy(sys);
}

/* ========================================================================= */
/* Nuclide operations                                                        */
/* ========================================================================= */

TEST(material_add_nuclide) {
    alea_system_t* sys = alea_create();
    int m = alea_add_material(sys, 1);

    /* Add U-235 and U-238 */
    ASSERT_EQ(alea_material_add_nuclide(sys, m, 92235, ".80c", 0.04), 0);
    ASSERT_EQ(alea_material_add_nuclide(sys, m, 92238, ".80c", 0.96), 0);

    ASSERT_EQ(alea_material_nuclide_count(sys, m), 2);

    int zaid;
    const char* lib;
    double frac;

    ASSERT_EQ(alea_material_nuclide_get(sys, m, 0, &zaid, &lib, &frac), 0);
    ASSERT_EQ(zaid, 92235);
    ASSERT_NOT_NULL(lib);
    ASSERT_STR_EQ(lib, ".80c");
    ASSERT_NEAR(frac, 0.04, 1e-10);

    ASSERT_EQ(alea_material_nuclide_get(sys, m, 1, &zaid, &lib, &frac), 0);
    ASSERT_EQ(zaid, 92238);
    ASSERT_NEAR(frac, 0.96, 1e-10);

    /* Out of range */
    ASSERT_EQ(alea_material_nuclide_get(sys, m, 2, &zaid, NULL, NULL), -1);

    alea_destroy(sys);
}

TEST(material_add_nuclide_null_library) {
    alea_system_t* sys = alea_create();
    int m = alea_add_material(sys, 1);

    ASSERT_EQ(alea_material_add_nuclide(sys, m, 1001, NULL, 0.5), 0);

    const char* lib = (const char*)0x1;  /* sentinel */
    ASSERT_EQ(alea_material_nuclide_get(sys, m, 0, NULL, &lib, NULL), 0);
    ASSERT_NULL(lib);

    alea_destroy(sys);
}

TEST(material_add_nuclide_invalid) {
    alea_system_t* sys = alea_create();

    /* Invalid material index */
    ASSERT_EQ(alea_material_add_nuclide(sys, -1, 92235, NULL, 0.5), -1);
    ASSERT_EQ(alea_material_add_nuclide(sys, 0, 92235, NULL, 0.5), -1);  /* no materials yet */
    ASSERT_EQ(alea_material_add_nuclide(NULL, 0, 92235, NULL, 0.5), -1);

    alea_destroy(sys);
}

/* ========================================================================= */
/* Element operations                                                        */
/* ========================================================================= */

TEST(material_add_element) {
    alea_system_t* sys = alea_create();
    int m = alea_add_material(sys, 1);

    /* Add iron (Z=26) */
    ASSERT_EQ(alea_material_add_element(sys, m, 26, ".80c", 1.0), 0);
    ASSERT_EQ(alea_material_element_count(sys, m), 1);

    int Z;
    double frac;
    ASSERT_EQ(alea_material_element_get(sys, m, 0, &Z, NULL, &frac), 0);
    ASSERT_EQ(Z, 26);
    ASSERT_NEAR(frac, 1.0, 1e-10);

    alea_destroy(sys);
}

TEST(material_expand_elements) {
    alea_system_t* sys = alea_create();
    int m = alea_add_material(sys, 1);

    /* Add hydrogen (Z=1) — has 2 natural isotopes (H-1, H-2) */
    ASSERT_EQ(alea_material_add_element(sys, m, 1, NULL, 1.0), 0);
    ASSERT_EQ(alea_material_element_count(sys, m), 1);
    ASSERT_EQ(alea_material_nuclide_count(sys, m), 0);

    /* Expand to nuclides */
    ASSERT_EQ(alea_material_expand_elements(sys, m), 0);
    ASSERT_EQ(alea_material_element_count(sys, m), 0);  /* elements cleared */
    ASSERT(alea_material_nuclide_count(sys, m) >= 2);  /* at least H-1, H-2 */

    /* First nuclide should be H-1 (zaid=1001) */
    int zaid;
    ASSERT_EQ(alea_material_nuclide_get(sys, m, 0, &zaid, NULL, NULL), 0);
    ASSERT_EQ(zaid, 1001);

    alea_destroy(sys);
}

TEST(material_add_element_invalid_Z) {
    alea_system_t* sys = alea_create();
    int m = alea_add_material(sys, 1);

    ASSERT_EQ(alea_material_add_element(sys, m, 0, NULL, 1.0), -1);
    ASSERT_EQ(alea_material_add_element(sys, m, 999, NULL, 1.0), -1);

    alea_destroy(sys);
}

/* ========================================================================= */
/* Density and fraction type                                                 */
/* ========================================================================= */

TEST(material_density) {
    alea_system_t* sys = alea_create();
    int m = alea_add_material(sys, 1);

    /* Initially no density set */
    double dens;
    bool has;
    ASSERT_EQ(alea_material_get_density(sys, m, &dens, &has), 0);
    ASSERT_FALSE(has);

    /* Set density */
    ASSERT_EQ(alea_material_set_density(sys, m, 7.87), 0);
    ASSERT_EQ(alea_material_get_density(sys, m, &dens, &has), 0);
    ASSERT_TRUE(has);
    ASSERT_NEAR(dens, 7.87, 1e-10);

    alea_destroy(sys);
}

TEST(material_weight_fraction) {
    alea_system_t* sys = alea_create();
    int m = alea_add_material(sys, 1);

    /* Default: atom fractions */
    ASSERT_FALSE(alea_material_is_weight_fraction(sys, m));

    ASSERT_EQ(alea_material_set_weight_fraction(sys, m, true), 0);
    ASSERT_TRUE(alea_material_is_weight_fraction(sys, m));

    ASSERT_EQ(alea_material_set_weight_fraction(sys, m, false), 0);
    ASSERT_FALSE(alea_material_is_weight_fraction(sys, m));

    alea_destroy(sys);
}

TEST(material_density_invalid) {
    alea_system_t* sys = alea_create();

    ASSERT_EQ(alea_material_set_density(sys, 0, 1.0), -1);
    ASSERT_EQ(alea_material_set_density(NULL, 0, 1.0), -1);

    alea_destroy(sys);
}

/* ========================================================================= */
/* Mixture operations                                                        */
/* ========================================================================= */

TEST(mixture_create_and_query) {
    alea_system_t* sys = alea_create();

    /* Create two base materials */
    alea_add_material(sys, 1);
    alea_add_material(sys, 2);

    ASSERT_EQ(alea_mixture_count(sys), 0);

    /* Create mixture of mat 1 (60%) and mat 2 (40%) */
    int mat_ids[] = {1, 2};
    double fracs[] = {0.6, 0.4};
    int mix_id = alea_create_mixture(sys, mat_ids, fracs, 2, 100);
    ASSERT(mix_id > 0);

    ASSERT_EQ(alea_mixture_count(sys), 1);
    ASSERT_EQ(alea_mixture_get_id(sys, 0), mix_id);
    ASSERT_EQ(alea_mixture_component_count(sys, 0), 2);

    int comp_mat;
    double comp_frac;
    ASSERT_EQ(alea_mixture_component_get(sys, 0, 0, &comp_mat, &comp_frac), 0);
    ASSERT_EQ(comp_mat, 1);
    ASSERT_NEAR(comp_frac, 0.6, 1e-10);

    ASSERT_EQ(alea_mixture_component_get(sys, 0, 1, &comp_mat, &comp_frac), 0);
    ASSERT_EQ(comp_mat, 2);
    ASSERT_NEAR(comp_frac, 0.4, 1e-10);

    /* Out of range */
    ASSERT_EQ(alea_mixture_component_get(sys, 0, 2, NULL, NULL), -1);

    alea_destroy(sys);
}

TEST(mixture_auto_id) {
    alea_system_t* sys = alea_create();
    alea_add_material(sys, 1);

    int ids[] = {1};
    double fracs[] = {1.0};
    int id1 = alea_create_mixture(sys, ids, fracs, 1, 0);
    int id2 = alea_create_mixture(sys, ids, fracs, 1, 0);
    ASSERT(id1 > 0);
    ASSERT(id2 > 0);
    ASSERT(id1 != id2);

    ASSERT_EQ(alea_mixture_count(sys), 2);

    alea_destroy(sys);
}

TEST(mixture_invalid_queries) {
    alea_system_t* sys = alea_create();

    ASSERT_EQ(alea_mixture_count(sys), 0);
    ASSERT_EQ(alea_mixture_get_id(sys, 0), -1);
    ASSERT_EQ(alea_mixture_component_count(sys, 0), 0);
    ASSERT_EQ(alea_mixture_component_get(sys, 0, 0, NULL, NULL), -1);
    ASSERT_EQ(alea_mixture_count(NULL), 0);

    alea_destroy(sys);
}

/* ========================================================================= */
/* Full roundtrip: create material, populate, query back                     */
/* ========================================================================= */

TEST(material_full_roundtrip) {
    alea_system_t* sys = alea_create();

    /* Create UO2 material */
    int m = alea_add_material(sys, 10);
    alea_material_set_weight_fraction(sys, m, true);
    alea_material_set_density(sys, m, 10.97);
    alea_material_add_nuclide(sys, m, 92235, ".80c", 0.04);
    alea_material_add_nuclide(sys, m, 92238, ".80c", 0.84);
    alea_material_add_nuclide(sys, m, 8016, ".80c", 0.12);

    /* Verify everything */
    ASSERT_EQ(alea_material_get_id(sys, m), 10);
    ASSERT_TRUE(alea_material_is_weight_fraction(sys, m));

    double dens;
    bool has;
    alea_material_get_density(sys, m, &dens, &has);
    ASSERT_TRUE(has);
    ASSERT_NEAR(dens, 10.97, 1e-10);

    ASSERT_EQ(alea_material_nuclide_count(sys, m), 3);

    int zaid;
    double frac;
    alea_material_nuclide_get(sys, m, 2, &zaid, NULL, &frac);
    ASSERT_EQ(zaid, 8016);
    ASSERT_NEAR(frac, 0.12, 1e-10);

    alea_destroy(sys);
}

TEST_MAIN()
