// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
// SPDX-License-Identifier: MPL-2.0
#include "alea_test.h"
#include "alea_openmc.h"
#include "core/alea_system.h"
#include "util/alea_alloc.h"

static openmc_model_t* rect(const char* dims, const char* pitch, const char* universes) {
    char text[2048];
    snprintf(text, sizeof(text), "<geometry><lattice id='10'><dimension>%s</dimension>"
        "<pitch>%s</pitch><lower_left>0 0</lower_left><universes>%s</universes>"
        "</lattice></geometry>", dims, pitch, universes);
    return openmc_load_string(text, 0);
}

TEST(rectangular_lattice_rejects_invalid_sizes_pitch_and_counts) {
    const char* dimensions[] = {"65536 65536", "2147483647 2147483647", "-2147483648 1",
        "0 1", "-1 1", "1001 1000", "1 1 1 1", "1 garbage", "4294967297 1"};
    for (size_t i = 0; i < sizeof(dimensions)/sizeof(dimensions[0]); i++)
        ASSERT_NULL(rect(dimensions[i], "1 1", "1"));
    const char* pitches[] = {"0 0", "-1 1", "nan 1", "inf 1", "1", "1 1 1", "1 1x"};
    for (size_t i = 0; i < sizeof(pitches)/sizeof(pitches[0]); i++)
        ASSERT_NULL(rect("1 1", pitches[i], "1"));
    const char* universes[] = {"", "1 1", "-1", "4294967297", "1garbage"};
    for (size_t i = 0; i < sizeof(universes)/sizeof(universes[0]); i++)
        ASSERT_NULL(rect("1 1", "1 1", universes[i]));
    ASSERT_NULL(rect("2 1", "1 1", "1"));
    openmc_model_t* good = rect("1 1", "1 1", "1");
    ASSERT_NOT_NULL(good);
    ASSERT_EQ(good->sys->cells.data[0].lat_fill_count, 1);
    openmc_model_destroy(good);
}

TEST(hex_lattice_rejects_overflow_and_incomplete_visual_input) {
    const char* cases[] = {
        "<hex_lattice id='10' n_rings='2147483647'><center>0 0</center><pitch>1</pitch><universes>1</universes></hex_lattice>",
        "<hex_lattice id='10' n_rings='2' n_axial='2147483647'><center>0 0 0</center><pitch>1 1</pitch><universes>1</universes></hex_lattice>",
        "<hex_lattice id='10' n_rings='1' n_axial='-2147483648'><center>0 0</center><pitch>1</pitch><universes>1</universes></hex_lattice>",
        "<hex_lattice id='10' n_rings='2'><center>0 0</center><pitch>1</pitch><universes>1</universes></hex_lattice>",
        "<hex_lattice id='10' n_rings='1'><center>0 0</center><pitch>0</pitch><universes>1</universes></hex_lattice>"
    };
    for (size_t i = 0; i < sizeof(cases)/sizeof(cases[0]); i++) {
        char text[2048];
        snprintf(text, sizeof(text), "<geometry>%s</geometry>", cases[i]);
        ASSERT_NULL(openmc_load_string(text, 0));
    }
    const char* valid = "<geometry><hex_lattice id='10' n_rings='2' n_axial='1'>"
        "<center>0 0 0</center><pitch>1 1</pitch><universes>1 1 1 1 1 1 1</universes>"
        "</hex_lattice></geometry>";
    openmc_model_t* model = openmc_load_string(valid, 0);
    ASSERT_NOT_NULL(model);
    ASSERT_EQ(model->sys->cells.data[0].lat_fill_count, 9);
    openmc_model_destroy(model);
}

TEST(lattice_failure_discards_prior_lattices_without_double_free) {
    const char* text = "<geometry><lattice id='10'><dimension>1 1</dimension>"
        "<pitch>1 1</pitch><lower_left>0 0</lower_left><universes>1</universes></lattice>"
        "<lattice id='11'><dimension>1 1</dimension><pitch>0 0</pitch>"
        "<lower_left>0 0</lower_left><universes>1</universes></lattice></geometry>";
    ASSERT_NULL(openmc_load_string(text, 0));
}

TEST(lattice_conversion_allocation_failures_are_atomic_and_retryable) {
    const char* text = "<geometry><surface id='1' type='sphere' coeffs='0 0 0 1'/>"
        "<cell id='1' material='void' region='-1' universe='1'/>"
        "<cell id='2' fill='10' region='-1'/>"
        "<lattice id='10'><dimension>1 1</dimension><pitch>1 1</pitch>"
        "<lower_left>0 0</lower_left><universes>1</universes></lattice>"
        "<hex_lattice id='20' n_rings='1'><center>0 0</center><pitch>1</pitch>"
        "<universes>1</universes></hex_lattice></geometry>";
    long failure;
    for (failure = 0; failure < 256; failure++) {
        alea_alloc_fail_after(failure);
        openmc_model_t* model = openmc_load_string(text, 0);
        alea_alloc_fail_after(-1);
        if (model) {
            ASSERT_EQ(model->sys->cells.count, 4);
            openmc_model_destroy(model);
            break;
        }
    }
    ASSERT(failure > 10 && failure < 256);
    openmc_model_t* retry = openmc_load_string(text, 0);
    ASSERT_NOT_NULL(retry);
    openmc_model_destroy(retry);
}

TEST_MAIN()
