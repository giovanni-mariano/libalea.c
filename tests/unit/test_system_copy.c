// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
// SPDX-License-Identifier: MPL-2.0
#include "alea_test.h"
#include "alea.h"
#include "alea_xml.h"
#include "core/alea_system.h"
#include "util/alea_alloc.h"
#include <limits.h>

/* Multiple owners at each level expose failures before and after a prior
 * entry has been copied. Prepared caches exercise borrowed index pointers. */
static alea_system_t* copy_fixture(void) {
    alea_system_t* sys = alea_create();
    if (!sys) return NULL;
    for (int i = 0; i < 2; i++) {
        int mat = alea_add_material(sys, 20 + 10*i);
        if (mat < 0 || alea_material_add_nuclide(sys, mat, 92235 + 3*i, "80c", 0.75) ||
            alea_material_add_nuclide(sys, mat, 1001, "81c", 0.25) ||
            alea_material_add_element(sys, mat, 8, "80c", 0.5)) goto fail;
        alea_material_t* material = &sys->materials.data[mat];
        if (alea_vec_reserve(&material->thermal_laws, 2, alea_thermal_law_t)) goto fail;
        for (int j = 0; j < 2; j++) {
            alea_thermal_law_t law = {alea_alloc_strdup("lwtr.20t"), 1001};
            if (!law.identifier) goto fail;
            material->thermal_laws.data[material->thermal_laws.count++] = law;
        }
        material->name = alea_alloc_strdup("fuel & moderator");
        material->comments = alea_alloc_strdup("material comment");
        if (!material->name || !material->comments) goto fail;
        material->has_standard_density = true;
        material->standard_density = 2.5;
        material->is_weight_fraction = true;
        int surface = alea_sphere_surface(sys, i+1, 10.0*(i+1), 0, 0, 1);
        if (surface < 0 || alea_add_cell(sys, i+1, alea_halfspace(sys, surface, -1),
                                       mat, -1.5, 0) < 0) goto fail;
        alea_cell_entry_t* cell = &sys->cells.data[i];
        cell->lat_fill = alea_calloc(3, sizeof(int));
        cell->lat_fill_count = 3;
        cell->comments = alea_alloc_strdup("cell comment");
        cell->inline_comment = alea_alloc_strdup("inline comment");
        if (!cell->lat_fill || !cell->comments || !cell->inline_comment) goto fail;
        cell->lat_fill[0] = 7;
    }
    int ids[] = {20, 30};
    double fractions[] = {0.25, 0.75};
    if (alea_create_mixture(sys, ids, fractions, 2, 40) != 40) goto fail;
    sys->mixtures.data[0].name = alea_alloc_strdup("mixture name");
    sys->mixtures.data[0].comments = alea_alloc_strdup("mixture comment");
    sys->mixtures.data[0].mc_material_id = 50;
    if (!sys->mixtures.data[0].name || !sys->mixtures.data[0].comments) goto fail;
    int surface = alea_sphere_surface(sys, 3, 30, 0, 0, 1);
    if (surface < 0 || alea_add_cell(sys, 3, alea_halfspace(sys, surface, -1),
                                   ALEA_MATERIAL_VOID, -1, 0) < 0) goto fail;
    sys->cells.data[2].material_id = 40;
    if (alea_prepare_query_acceleration(sys)) goto fail;
    return sys;
fail:
    alea_destroy(sys);
    return NULL;
}

static void check_fixture(const alea_system_t* sys) {
    ASSERT_EQ(sys->cells.count, 3);
    for (size_t i = 0; i < 2; i++) {
        ASSERT_STR_EQ(sys->cells.data[i].comments, "cell comment");
        ASSERT_STR_EQ(sys->cells.data[i].inline_comment, "inline comment");
        ASSERT_EQ(sys->cells.data[i].lat_fill[0], 7);
        ASSERT_STR_EQ(sys->materials.data[i].name, "fuel & moderator");
        ASSERT_STR_EQ(sys->materials.data[i].comments, "material comment");
        ASSERT_STR_EQ(sys->materials.data[i].nuclides.data[1].library, "81c");
        ASSERT_STR_EQ(sys->materials.data[i].elements.data[0].library, "80c");
        ASSERT_STR_EQ(sys->materials.data[i].thermal_laws.data[1].identifier, "lwtr.20t");
    }
    ASSERT_STR_EQ(sys->mixtures.data[0].name, "mixture name");
    ASSERT_STR_EQ(sys->mixtures.data[0].comments, "mixture comment");
    ASSERT_EQ(sys->mixtures.data[0].components.data[1].material_id, 30);
}

TEST(clone_failure_never_frees_source_owners) {
    alea_system_t* source = copy_fixture();
    ASSERT_NOT_NULL(source);
    long failures = 0;
    alea_system_t* clone = NULL;
    for (; failures < 256; failures++) {
        alea_alloc_fail_after(failures);
        clone = alea_clone(source);
        alea_alloc_fail_after(-1);
        check_fixture(source);
        ASSERT_EQ(alea_find_cell(source, 10, 0, 0), 0);
        if (clone) break;
        ASSERT_EQ(alea_error_code(), ALEA_ERR_OUT_OF_MEMORY);
    }
    ASSERT_NOT_NULL(clone);
    ASSERT(failures > 30);
    for (size_t i = 0; i < 2; i++) {
        ASSERT(clone->cells.data[i].lat_fill != source->cells.data[i].lat_fill);
        ASSERT(clone->cells.data[i].comments != source->cells.data[i].comments);
        ASSERT(clone->materials.data[i].nuclides.data != source->materials.data[i].nuclides.data);
        ASSERT(clone->materials.data[i].thermal_laws.data[1].identifier !=
               source->materials.data[i].thermal_laws.data[1].identifier);
    }
    ASSERT(clone->mixtures.data[0].components.data != source->mixtures.data[0].components.data);
    alea_destroy(source);
    check_fixture(clone);
    ASSERT_EQ(alea_prepare_query_acceleration(clone), 0);
    ASSERT_EQ(alea_find_cell(clone, 20, 0, 0), 1);
    alea_destroy(clone);
}

TEST(clone_preserves_automatic_ids) {
    alea_system_t* source = copy_fixture();
    ASSERT_NOT_NULL(source);
    alea_system_t* clone = alea_clone(source);
    ASSERT_NOT_NULL(clone);
    int index = alea_add_cell(clone, 0, clone->cells.data[0].root_node_id,
                              ALEA_MATERIAL_VOID, 0, 0);
    ASSERT_EQ(index, 3);
    ASSERT_EQ(clone->cells.data[index].mc_cell_id, source->next_auto_cell_id);
    ASSERT_EQ(alea_find_cell_by_id(clone, 1), 0);
    ASSERT_EQ(alea_validate_cell_ids(clone), 0);
    int material = alea_add_material(clone, 0);
    ASSERT(material >= 0);
    ASSERT_EQ(alea_material_get_id(clone, material), source->next_auto_material_id);
    alea_destroy(source);
    alea_destroy(clone);
}

TEST(merge_preserves_composition_and_xml_roundtrip) {
    alea_system_t* source = copy_fixture();
    alea_system_t* target = alea_create();
    ASSERT_NOT_NULL(source); ASSERT_NOT_NULL(target);
    int hydrogen = alea_add_material(target, 10);
    ASSERT(hydrogen >= 0);
    ASSERT_EQ(alea_material_add_nuclide(target, hydrogen, 1001, "80c", 1), 0);
    ASSERT_EQ(alea_merge(target, source, 100), 3);
    ASSERT_EQ(target->materials.count, 3);
    ASSERT_EQ(target->mixtures.count, 1);
    for (int i = 0; i < 2; i++) {
        const alea_cell_entry_t* cell = &target->cells.data[i];
        ASSERT_EQ(cell->material_id, 120 + 10*i);
        ASSERT_EQ(cell->material_index, i+1);
        ASSERT_EQ(target->materials.data[i+1].material_id, cell->material_id);
        ASSERT_EQ(target->materials.data[i+1].nuclides.data[0].zaid, 92235 + 3*i);
    }
    ASSERT_EQ(target->cells.data[2].material_id, 140);
    ASSERT_EQ(target->cells.data[2].material_index, -1);
    ASSERT_EQ(target->mixtures.data[0].mc_material_id, 150);
    ASSERT_EQ(target->mixtures.data[0].components.data[0].material_id, 120);
    ASSERT_EQ(target->mixtures.data[0].components.data[1].material_id, 130);
    ASSERT_EQ(target->next_auto_material_id, 131);
    alea_destroy(source);
    const char* path = "test_system_copy.xml";
    ASSERT_EQ(alea_xml_export_system(target, path), 0);
    alea_model_t* model = alea_xml_load(path);
    remove(path);
    ASSERT_NOT_NULL(model);
    alea_system_t* loaded = alea_model_system(model);
    ASSERT_EQ(alea_material_count(loaded), 3);
    int index = alea_find_material_by_id(loaded, 120);
    ASSERT(index >= 0);
    ASSERT_EQ(loaded->materials.data[index].nuclides.data[0].zaid, 92235);
    ASSERT_STR_EQ(loaded->materials.data[index].thermal_laws.data[1].identifier, "lwtr.20t");
    ASSERT_STR_EQ(loaded->materials.data[index].comments, "material comment");
    ASSERT_TRUE(loaded->materials.data[index].is_weight_fraction);
    ASSERT_NEAR(loaded->materials.data[index].standard_density, 2.5, 1e-12);
    ASSERT_EQ(loaded->mixtures.data[0].components.data[1].material_id, 130);
    ASSERT_EQ(loaded->cells.data[0].material_index, index);
    ASSERT_EQ(loaded->cells.data[2].material_index, -1);
    alea_model_destroy(model);
    alea_destroy(target);
}

TEST(merge_rejects_material_conflicts_before_mutation) {
    alea_system_t* source = copy_fixture();
    alea_system_t* target = alea_create();
    ASSERT_NOT_NULL(source); ASSERT_NOT_NULL(target);
    ASSERT(alea_add_material(target, 120) >= 0);
    uint64_t generation = alea_system_geometry_generation(target);
    ASSERT_EQ(alea_merge(target, source, 100), -1);
    ASSERT_EQ(alea_error_code(), ALEA_ERR_INVALID_ARG);
    ASSERT_EQ(target->materials.count, 1);
    ASSERT_EQ(target->cells.count, 0);
    ASSERT_EQ(target->nodes.count, 0);
    ASSERT_EQ(alea_system_geometry_generation(target), generation);
    ASSERT_EQ(alea_merge(target, source, INT_MAX), -1);
    ASSERT_EQ(alea_merge(target, source, -30), -1);
    alea_destroy(target);
    target = alea_create();
    ASSERT_NOT_NULL(target);
    source->cells.data[0].material_index = 99;
    ASSERT_EQ(alea_merge(target, source, 100), -1);
    ASSERT_EQ(target->materials.count, 0);
    alea_destroy(target);
    alea_destroy(source);
}

TEST(merge_rejects_invalid_mixture_references_and_conflicts) {
    alea_system_t* source = copy_fixture();
    alea_system_t* target = alea_create();
    ASSERT_NOT_NULL(source); ASSERT_NOT_NULL(target);
    source->mixtures.data[0].components.data[1].material_id = 999;
    ASSERT_EQ(alea_merge(target, source, 100), -1);
    ASSERT_EQ(target->materials.count, 0);
    source->mixtures.data[0].components.data[1].material_id = 30;
    source->mixtures.data[0].mc_material_id = INT_MAX;
    ASSERT_EQ(alea_merge(target, source, 100), -1);
    source->mixtures.data[0].mc_material_id = 0;
    alea_mixture_t* conflict = alea_mixture_create(140);
    ASSERT_NOT_NULL(conflict);
    ASSERT_EQ(alea_add_mixture(target, conflict), 0);
    alea_mixture_destroy(conflict);
    ASSERT_EQ(alea_merge(target, source, 100), -1);
    ASSERT_EQ(target->materials.count, 0);
    ASSERT_EQ(target->mixtures.count, 1);
    ASSERT_EQ(target->nodes.count, 0);
    alea_destroy(source);
    alea_destroy(target);
}

TEST(merge_failure_rolls_back_owned_materials_and_allows_retry) {
    alea_system_t* source = copy_fixture();
    ASSERT_NOT_NULL(source);
    int succeeded = 0;
    long failures = 0;
    for (; failures < 256; failures++) {
        alea_system_t* target = alea_create();
        ASSERT_NOT_NULL(target);
        ASSERT(alea_add_material(target, 10) >= 0);
        int surface = alea_sphere_surface(target, 1, 0, 0, 0, 1);
        ASSERT(surface >= 0);
        ASSERT_EQ(alea_add_cell(target, 1, alea_halfspace(target, surface, -1), 0, -1, 0), 0);
        ASSERT_EQ(alea_prepare_query_acceleration(target), 0);
        size_t nodes = target->nodes.count, surfaces = target->surfaces.count;
        int material_id = target->next_auto_material_id;
        alea_error_clear();
        alea_alloc_fail_after(failures);
        int result = alea_merge(target, source, 100);
        alea_alloc_fail_after(-1);
        check_fixture(source);
        if (result == 3) {
            succeeded = 1;
            alea_destroy(target);
            break;
        }
        ASSERT_EQ(result, -1);
        ASSERT_EQ(alea_error_code(), ALEA_ERR_OUT_OF_MEMORY);
        ASSERT_EQ(target->cells.count, 1);
        ASSERT_EQ(target->surfaces.count, surfaces);
        ASSERT_EQ(target->nodes.count, nodes);
        ASSERT_EQ(target->materials.count, 1);
        ASSERT_EQ(target->mixtures.count, 0);
        ASSERT_EQ(target->next_auto_material_id, material_id);
        ASSERT_EQ(alea_prepare_query_acceleration(target), 0);
        ASSERT_EQ(alea_find_cell(target, 0, 0, 0), 0);
        ASSERT_EQ(alea_find_cell_by_id(target, 101), -1);
        ASSERT_EQ(alea_merge(target, source, 100), 3);
        alea_destroy(target);
    }
    ASSERT(succeeded);
    ASSERT(failures > 25);
    alea_destroy(source);
}

TEST_MAIN()
