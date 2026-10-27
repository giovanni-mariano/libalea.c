// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
// SPDX-License-Identifier: MPL-2.0
#include "alea_test.h"
#include "alea_model.h"
#include "alea_mcnp.h"
#include "core/alea_system.h"
#include "util/alea_alloc.h"
#include "util/cell_parameter.h"

typedef struct {
    alea_system_t* sys;
    alea_model_t* alea;
    mcnp_model_t* mcnp;
} wrapper_t;

static wrapper_t wrap(alea_system_t* sys, int kind) {
    wrapper_t result = {.sys = sys};
    if (kind) result.mcnp = mcnp_model_wrap(sys);
    else result.alea = alea_model_wrap(sys);
    return result;
}
static void destroy(wrapper_t* wrapper) {
    alea_model_destroy(wrapper->alea);
    mcnp_model_destroy(wrapper->mcnp);
    alea_destroy(wrapper->sys);
}
static size_t count(const wrapper_t* w) {
    return w->alea ? alea_model_cell_metadata_count(w->alea) : w->mcnp->cell_params_count;
}
static int set_name(wrapper_t* w, size_t i, const char* name) {
    return w->alea ? alea_model_cell_set_name(w->alea, i, name)
                   : mcnp_model_cell_set_name(w->mcnp, i, name);
}
static const char* name_at(const wrapper_t* w, size_t i) {
    return w->alea ? alea_model_cell_metadata(w->alea, i)->name
                   : mcnp_model_cell_name(w->mcnp, i);
}
static double importance(const wrapper_t* w, size_t i) {
    return w->alea ? alea_model_cell_metadata(w->alea, i)->importance_neutron
                   : mcnp_cell_params_const(w->mcnp, i)->imp_n;
}
static void set_importance(wrapper_t* w, size_t i, double value) {
    if (w->alea) alea_model_cell_metadata_mut(w->alea, i)->importance_neutron = value;
    else mcnp_cell_params(w->mcnp, i)->imp_n = value;
}
static alea_node_id_t sphere(alea_system_t* sys, int id, double x) {
    int surface = alea_sphere_surface(sys, id, x, 0, 0, 1);
    return surface < 0 ? ALEA_NODE_ID_INVALID : alea_halfspace(sys, surface, -1);
}

TEST(model_merge_reset_and_later_additions_keep_metadata_aligned) {
    for (int kind = 0; kind < 2; kind++) {
        alea_system_t* sys = alea_create();
        alea_system_t* source = alea_create();
        ASSERT_NOT_NULL(sys); ASSERT_NOT_NULL(source);
        ASSERT_EQ(alea_add_cell(sys, 1, sphere(sys, 1, 0), ALEA_MATERIAL_VOID, 0, 0), 0);
        ASSERT_EQ(alea_add_cell(source, 1, sphere(source, 1, 10), ALEA_MATERIAL_VOID, 0, 0), 0);
        wrapper_t w = wrap(sys, kind);
        ASSERT(w.alea || w.mcnp);
        ASSERT_EQ(set_name(&w, 0, "original"), 0);
        set_importance(&w, 0, 0);
        ASSERT_EQ(alea_merge(sys, source, 100), 1);
        ASSERT_EQ(count(&w), 2);
        ASSERT_STR_EQ(name_at(&w, 0), "original");
        ASSERT_NULL(name_at(&w, 1));
        ASSERT_NEAR(importance(&w, 1), 1, 0);
        ASSERT_EQ(set_name(&w, 1, "merged"), 0);
        ASSERT_EQ(alea_add_cell(sys, 200, sphere(sys, 200, 20), ALEA_MATERIAL_VOID, 0, 0), 2);
        ASSERT_EQ(count(&w), 3);
        alea_reset(sys);
        ASSERT_EQ(count(&w), 0);
        ASSERT_EQ(alea_add_cell(sys, 1, sphere(sys, 1, 30), ALEA_MATERIAL_VOID, 0, 0), 0);
        ASSERT_EQ(count(&w), 1);
        ASSERT_NULL(name_at(&w, 0));
        ASSERT_NEAR(importance(&w, 0), 1, 0);
        alea_destroy(source);
        destroy(&w);
    }
}

TEST(model_reservation_failures_prevent_cell_publication_and_allow_retry) {
    for (int kind = 0; kind < 2; kind++) {
        alea_system_t* sys = alea_create();
        alea_system_t* source = alea_create();
        ASSERT_NOT_NULL(sys); ASSERT_NOT_NULL(source);
        alea_node_id_t root = sphere(sys, 1, 0);
        int original_count = kind ? 64 : 8;
        for (int i = 0; i < original_count; i++)
            ASSERT_EQ(alea_add_cell(sys, i+1, root, ALEA_MATERIAL_VOID, 0, 0), i);
        wrapper_t w = wrap(sys, kind);
        ASSERT(w.alea || w.mcnp);
        ASSERT_EQ(set_name(&w, 0, "preserved"), 0);
        ASSERT_EQ(alea_add_cell(source, 1, sphere(source, 1, 10), ALEA_MATERIAL_VOID, 0, 0), 0);
        for (long failure = 0; failure < (kind ? 2 : 1); failure++) {
            alea_alloc_fail_after(failure);
            int result = alea_add_cell(sys, 0, root, ALEA_MATERIAL_VOID, 0, 0);
            alea_alloc_fail_after(-1);
            ASSERT_EQ(result, -1);
            ASSERT_EQ(sys->cells.count, original_count);
            ASSERT_EQ(count(&w), original_count);
            ASSERT_STR_EQ(name_at(&w, 0), "preserved");
            alea_alloc_fail_after(failure);
            result = alea_merge(sys, source, 1000);
            alea_alloc_fail_after(-1);
            ASSERT_EQ(result, -1);
            ASSERT_EQ(sys->cells.count, original_count);
            ASSERT_EQ(count(&w), original_count);
        }
        ASSERT_EQ(alea_merge(sys, source, 1000), 1);
        ASSERT_EQ(count(&w), (size_t)original_count + 1);
        ASSERT_EQ(alea_add_cell(sys, 0, root, ALEA_MATERIAL_VOID, 0, 0), original_count+1);
        ASSERT_EQ(count(&w), (size_t)original_count + 2);
        alea_destroy(source);
        destroy(&w);
    }
}

TEST(model_merge_rolls_back_without_publishing_sidecar_entries) {
    for (int kind = 0; kind < 2; kind++) {
        int succeeded = 0;
        for (long failure = 0; failure < 32; failure++) {
            alea_system_t* sys = alea_create();
            alea_system_t* source = alea_create();
            ASSERT_NOT_NULL(sys); ASSERT_NOT_NULL(source);
            ASSERT_EQ(alea_add_cell(sys, 1, sphere(sys, 1, 0), ALEA_MATERIAL_VOID, 0, 0), 0);
            ASSERT_EQ(alea_add_cell(source, 1, sphere(source, 1, 10), ALEA_MATERIAL_VOID, 0, 0), 0);
            ASSERT_EQ(alea_cell_set_comment(source, 0, "owned comment"), 0);
            wrapper_t w = wrap(sys, kind);
            ASSERT(w.alea || w.mcnp);
            ASSERT_EQ(set_name(&w, 0, "preserved"), 0);
            alea_alloc_fail_after(failure);
            int result = alea_merge(sys, source, 100);
            alea_alloc_fail_after(-1);
            if (result == 1) {
                succeeded = 1;
                ASSERT_EQ(count(&w), 2);
            } else {
                ASSERT_EQ(result, -1);
                ASSERT_EQ(count(&w), 1);
                ASSERT_EQ(sys->cells.count, 1);
                ASSERT_STR_EQ(name_at(&w, 0), "preserved");
                ASSERT_EQ(alea_merge(sys, source, 100), 1);
                ASSERT_EQ(count(&w), 2);
            }
            alea_destroy(source);
            destroy(&w);
            if (succeeded) break;
        }
        ASSERT(succeeded);
    }
}

TEST(model_split_copy_failures_restore_original_geometry_and_metadata) {
    for (int kind = 0; kind < 2; kind++) {
        int succeeded = 0;
        for (long failure = 0; failure < 32; failure++) {
            alea_system_t* sys = alea_create();
            ASSERT_NOT_NULL(sys);
            alea_node_id_t left = sphere(sys, 1, -3), right = sphere(sys, 2, 3);
            alea_node_id_t root = alea_union(sys, left, right);
            ASSERT_EQ(alea_add_cell(sys, 1, root, ALEA_MATERIAL_VOID, 0, 0), 0);
            wrapper_t w = wrap(sys, kind);
            ASSERT(w.alea || w.mcnp);
            ASSERT_EQ(set_name(&w, 0, "copied name"), 0);
            set_importance(&w, 0, 0.25);
            if (w.alea) ASSERT_EQ(alea_model_cell_set_detector_probability(w.alea, 0, 4, 0.5), 0);
            else {
                mcnp_scoped_cell_params_t* params = &mcnp_cell_params(w.mcnp, 0)->scoped;
                ASSERT_EQ(alea_detector_probability_store(&params->detector_probabilities,
                    &params->detector_probability_count, 4, 0.5), 0);
            }
            alea_alloc_fail_after(failure);
            int result = alea_split_union_cells(sys);
            alea_alloc_fail_after(-1);
            if (result < 0) {
                ASSERT_EQ(sys->cells.count, 1);
                ASSERT_EQ(count(&w), 1);
                ASSERT_EQ(sys->cells.data[0].root_node_id, root);
                ASSERT_EQ(alea_find_cell_by_id(sys, 1), 0);
                ASSERT_STR_EQ(name_at(&w, 0), "copied name");
                ASSERT_EQ(alea_split_union_cells(sys), 2);
            } else { ASSERT_EQ(result, 2); succeeded = 1; }
            ASSERT_EQ(sys->cells.count, 2);
            ASSERT_EQ(count(&w), 2);
            ASSERT_EQ(alea_find_cell_by_id(sys, 1), -1);
            for (size_t i = 0; i < 2; i++) {
                ASSERT_STR_EQ(name_at(&w, i), "copied name");
                ASSERT_NEAR(importance(&w, i), 0.25, 0);
            }
            const alea_detector_probability_t* entries[2];
            for (size_t i = 0; i < 2; i++) {
                size_t n;
                if (w.alea) {
                    const alea_model_cell_metadata_t* meta = alea_model_cell_metadata(w.alea, i);
                    entries[i] = meta->detector_probabilities; n = meta->detector_probability_count;
                } else {
                    const mcnp_scoped_cell_params_t* params = &mcnp_cell_params_const(w.mcnp, i)->scoped;
                    entries[i] = params->detector_probabilities; n = params->detector_probability_count;
                }
                ASSERT_EQ(n, 1);
                ASSERT_EQ(entries[i][0].tally, 4);
                ASSERT_NEAR(entries[i][0].probability, 0.5, 0);
            }
            ASSERT(entries[0] != entries[1]);
            ASSERT_EQ(set_name(&w, 0, "independent"), 0);
            ASSERT_STR_EQ(name_at(&w, 1), "copied name");
            destroy(&w);
            if (succeeded) break;
        }
        ASSERT(succeeded);
    }
}

TEST(model_pruning_removes_matching_metadata_and_updates_id_lookup) {
    for (int kind = 0; kind < 2; kind++) {
        alea_system_t* sys = alea_create();
        ASSERT_NOT_NULL(sys);
        alea_node_id_t root = sphere(sys, 1, 0);
        ASSERT_EQ(alea_add_cell(sys, 1, root, ALEA_MATERIAL_VOID, 0, 0), 0);
        ASSERT_EQ(alea_add_cell(sys, 2, root, ALEA_MATERIAL_VOID, 0, 0), 1);
        wrapper_t w = wrap(sys, kind);
        ASSERT(w.alea || w.mcnp);
        ASSERT_EQ(set_name(&w, 0, "removed"), 0);
        ASSERT_EQ(set_name(&w, 1, "survivor"), 0);
        ASSERT_EQ(alea_prepare_query_acceleration(sys), 0);
        sys->cells.data[0].root_node_id = ALEA_NODE_ID_INVALID;
        alea_simplify_and_prune_cells(sys, NULL);
        ASSERT_EQ(sys->cells.count, 1);
        ASSERT_EQ(count(&w), 1);
        ASSERT_STR_EQ(name_at(&w, 0), "survivor");
        ASSERT_EQ(alea_find_cell_by_id(sys, 1), -1);
        ASSERT_EQ(alea_find_cell_by_id(sys, 2), 0);
        destroy(&w);
    }
}

TEST(model_batch_pruning_preserves_owned_entries_for_alternating_cells) {
    for (int kind = 0; kind < 2; kind++) {
        alea_system_t* sys = alea_create();
        ASSERT_NOT_NULL(sys);
        alea_node_id_t root = sphere(sys, 1, 0);
        for (int i = 0; i < 32; i++)
            ASSERT_EQ(alea_add_cell(sys, i+1, root, ALEA_MATERIAL_VOID, 0, 0), i);
        wrapper_t w = wrap(sys, kind);
        ASSERT(w.alea || w.mcnp);
        for (int i = 0; i < 32; i++) {
            char name[32]; snprintf(name, sizeof(name), "cell %d", i+1);
            ASSERT_EQ(set_name(&w, i, name), 0);
            set_importance(&w, i, i+1);
        }
        for (int i = 0; i < 32; i += 2) sys->cells.data[i].root_node_id = ALEA_NODE_ID_INVALID;
        alea_simplify_and_prune_cells(sys, NULL);
        ASSERT_EQ(count(&w), 16);
        for (int i = 0; i < 16; i++) {
            char name[32]; snprintf(name, sizeof(name), "cell %d", 2*i+2);
            ASSERT_STR_EQ(name_at(&w, i), name);
            ASSERT_NEAR(importance(&w, i), 2*i+2, 0);
            ASSERT_EQ(alea_find_cell_by_id(sys, 2*i+2), i);
        }
        destroy(&w);
    }
}

TEST_MAIN()
