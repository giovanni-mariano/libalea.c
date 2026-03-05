// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/*
 * test_materials.c - Tests for public material and mixture API
 */

#include "alea_test.h"
#include "alea.h"
#include "alea_mcnp.h"
#include "alea_openmc.h"
#include <string.h>

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

TEST(mixture_find_by_id) {
    alea_system_t* sys = alea_create();
    alea_add_material(sys, 1);
    alea_add_material(sys, 2);

    int mat_ids[] = {1, 2};
    double fracs[] = {0.5, 0.5};
    int mix_id = alea_create_mixture(sys, mat_ids, fracs, 2, 100);
    ASSERT_EQ(mix_id, 100);

    ASSERT_EQ(alea_find_mixture_by_id(sys, 100), 0);
    ASSERT_EQ(alea_find_mixture_by_id(sys, 999), -1);
    ASSERT_EQ(alea_find_mixture_by_id(NULL, 100), -1);

    alea_destroy(sys);
}

TEST(mixture_with_populated_materials) {
    alea_system_t* sys = alea_create();

    /* Material 1: U-235 enriched uranium */
    int m0 = alea_add_material(sys, 1);
    alea_material_add_nuclide(sys, m0, 92235, ".80c", 0.04);
    alea_material_add_nuclide(sys, m0, 92238, ".80c", 0.96);

    /* Material 2: oxygen */
    int m1 = alea_add_material(sys, 2);
    alea_material_add_nuclide(sys, m1, 8016, ".80c", 1.0);

    /* Mixture: 70% mat1 + 30% mat2 */
    int mat_ids[] = {1, 2};
    double fracs[] = {0.7, 0.3};
    int mix_id = alea_create_mixture(sys, mat_ids, fracs, 2, 100);
    ASSERT_EQ(mix_id, 100);

    /* Query back mixture components */
    int idx = alea_find_mixture_by_id(sys, 100);
    ASSERT(idx >= 0);
    ASSERT_EQ(alea_mixture_component_count(sys, idx), 2);

    int comp_mat;
    double comp_frac;
    alea_mixture_component_get(sys, idx, 0, &comp_mat, &comp_frac);
    ASSERT_EQ(comp_mat, 1);
    ASSERT_NEAR(comp_frac, 0.7, 1e-10);

    alea_mixture_component_get(sys, idx, 1, &comp_mat, &comp_frac);
    ASSERT_EQ(comp_mat, 2);
    ASSERT_NEAR(comp_frac, 0.3, 1e-10);

    alea_destroy(sys);
}

/* Helper: read FILE* into malloc'd string */
static char* read_file_to_string(FILE* f) {
    rewind(f);
    fseek(f, 0, SEEK_END);
    long len = ftell(f);
    rewind(f);
    char* buf = malloc(len + 1);
    if (!buf) return NULL;
    fread(buf, 1, len, f);
    buf[len] = '\0';
    return buf;
}

/* Helper: create a minimal exportable system with a mixture */
static alea_system_t* create_mixture_system(void) {
    alea_system_t* sys = alea_create();

    /* Materials */
    int m0 = alea_add_material(sys, 1);
    alea_material_add_nuclide(sys, m0, 92235, "80c", 0.04);
    alea_material_add_nuclide(sys, m0, 92238, "80c", 0.96);
    alea_material_set_weight_fraction(sys, m0, true);

    int m1 = alea_add_material(sys, 2);
    alea_material_add_nuclide(sys, m1, 8016, "80c", 1.0);
    alea_material_set_weight_fraction(sys, m1, true);

    /* Mixture: 60% mat1 + 40% mat2, ID=100 */
    int mat_ids[] = {1, 2};
    double mix_fracs[] = {0.6, 0.4};
    alea_create_mixture(sys, mat_ids, mix_fracs, 2, 100);

    /* Geometry: one sphere, two cells */
    int si = alea_sphere_surface(sys, 0, 0, 0, 0, 5.0);
    alea_node_id_t inside = alea_halfspace(sys, si, -1);
    alea_node_id_t outside = alea_halfspace(sys, si, +1);

    /* Cell 1: inside sphere, use mixture with density */
    int c0 = alea_add_cell(sys, 1, inside, -1, 0, 0);
    alea_cell_set_mixture(sys, c0, 100);
    alea_cell_set_density(sys, c0, -7.0);

    /* Cell 2: outside sphere, void */
    alea_add_cell(sys, 2, outside, -1, 0, 0);

    return sys;
}

TEST(mixture_mcnp_export) {
    alea_system_t* sys = create_mixture_system();

    FILE* f = tmpfile();
    ASSERT_NOT_NULL(f);
    int rc = mcnp_export_system_stream(sys, f);
    ASSERT_EQ(rc, 0);

    char* output = read_file_to_string(f);
    fclose(f);
    ASSERT_NOT_NULL(output);

    /* Verify mixture material card M100 appears */
    ASSERT(strstr(output, "M100") != NULL);

    /* Verify base materials are also exported */
    ASSERT(strstr(output, "M1") != NULL);
    ASSERT(strstr(output, "M2") != NULL);

    /* Verify mixture header comment */
    ASSERT(strstr(output, "Mixture M100") != NULL);

    /* Verify nuclides from component materials appear in mixture */
    /* U-235 from M1 should be in the mixture card */
    ASSERT(strstr(output, "92235") != NULL);
    /* O-16 from M2 should be in the mixture card */
    ASSERT(strstr(output, "8016") != NULL);

    free(output);
    alea_destroy(sys);
}

TEST(mixture_openmc_export) {
    alea_system_t* sys = create_mixture_system();

    FILE* f = tmpfile();
    ASSERT_NOT_NULL(f);
    int rc = openmc_export_system_stream(sys, f);
    ASSERT_EQ(rc, 0);

    char* output = read_file_to_string(f);
    fclose(f);
    ASSERT_NOT_NULL(output);

    /* Verify the output is valid XML with materials */
    ASSERT(strstr(output, "<materials>") != NULL);

    /* Verify nuclides appear */
    ASSERT(strstr(output, "U235") != NULL || strstr(output, "92235") != NULL);
    ASSERT(strstr(output, "O16") != NULL || strstr(output, "8016") != NULL);

    free(output);
    alea_destroy(sys);
}

TEST(mixture_multiple_components) {
    alea_system_t* sys = alea_create();

    /* Create 3 materials */
    int m0 = alea_add_material(sys, 1);
    alea_material_add_nuclide(sys, m0, 1001, NULL, 0.5);
    alea_material_add_nuclide(sys, m0, 1002, NULL, 0.5);

    int m1 = alea_add_material(sys, 2);
    alea_material_add_nuclide(sys, m1, 8016, NULL, 1.0);

    int m2 = alea_add_material(sys, 3);
    alea_material_add_nuclide(sys, m2, 26056, NULL, 1.0);

    /* Mixture of all three */
    int mat_ids[] = {1, 2, 3};
    double fracs[] = {0.5, 0.3, 0.2};
    int mix_id = alea_create_mixture(sys, mat_ids, fracs, 3, 50);
    ASSERT_EQ(mix_id, 50);

    int idx = alea_find_mixture_by_id(sys, 50);
    ASSERT_EQ(alea_mixture_component_count(sys, idx), 3);

    int comp_mat;
    double comp_frac;
    alea_mixture_component_get(sys, idx, 2, &comp_mat, &comp_frac);
    ASSERT_EQ(comp_mat, 3);
    ASSERT_NEAR(comp_frac, 0.2, 1e-10);

    alea_destroy(sys);
}

TEST(mixture_duplicate_id) {
    alea_system_t* sys = alea_create();
    alea_add_material(sys, 1);

    int ids[] = {1};
    double fracs[] = {1.0};

    /* First mixture with explicit ID */
    int id1 = alea_create_mixture(sys, ids, fracs, 1, 50);
    ASSERT_EQ(id1, 50);

    /* Second mixture with same ID — should still succeed (no dedup check) */
    int id2 = alea_create_mixture(sys, ids, fracs, 1, 50);
    ASSERT_EQ(id2, 50);
    ASSERT_EQ(alea_mixture_count(sys), 2);

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
