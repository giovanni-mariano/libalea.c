// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/*
 * test_openmc.c - OpenMC XML parse, conversion, and export tests
 */

#include "alea_test.h"
#include "alea.h"
#include "alea_mcnp.h"
#include "alea_openmc.h"
#include "core/alea_system.h"
#include "core/alea_export.h"
#include "core/alea_universe.h"
#include <string.h>
#include <stdio.h>

/* ========================================================================= */
/* Parse tests                                                               */
/* ========================================================================= */

TEST(openmc_parse_simple) {
    openmc_model_t* omc = openmc_load("tests/data/openmc_simple.xml");
    if (!omc) SKIP("Test data file not found");

    ASSERT(alea_cell_count(omc->sys) > 0);
    openmc_model_destroy(omc);
}

TEST(openmc_parse_lattice_rect) {
    openmc_model_t* omc = openmc_load("tests/data/openmc_lattice.xml");
    if (!omc) SKIP("Test data file not found");

    ASSERT(alea_cell_count(omc->sys) > 0);
    openmc_model_destroy(omc);
}

TEST(openmc_parse_lattice_hex) {
    openmc_model_t* omc = openmc_load("tests/data/openmc_hex_lattice.xml");
    if (!omc) SKIP("Test data file not found");

    ASSERT(alea_cell_count(omc->sys) > 0);
    openmc_model_destroy(omc);
}

/* ========================================================================= */
/* Inline OpenMC XML tests                                                   */
/* ========================================================================= */

TEST(openmc_sphere) {
    const char* xml =
        "<?xml version='1.0'?>\n"
        "<model>\n"
        " <geometry>\n"
        "  <surface id=\"1\" type=\"sphere\" coeffs=\"0.0 0.0 0.0 5.0\" />\n"
        "  <cell id=\"1\" material=\"1\" region=\"-1\" />\n"
        "  <cell id=\"2\" material=\"void\" region=\"1\" />\n"
        " </geometry>\n"
        " <materials>\n"
        "  <material id=\"1\"><density value=\"1.0\" />"
        "  <nuclide name=\"H1\" ao=\"1.0\" /></material>\n"
        " </materials>\n"
        "</model>\n";
    openmc_model_t* omc = openmc_load_string(xml, strlen(xml));
    if (!omc) SKIP("OpenMC string parse not supported");
    alea_system_t* sys = omc->sys;

    alea_build_universe_index(sys);
    ASSERT_EQ(alea_material_at(sys, 0, 0, 0), 1);
    ASSERT_EQ(alea_material_at(sys, 10, 0, 0), 0);

    openmc_model_destroy(omc);
}

TEST(openmc_cylinder_z) {
    const char* xml =
        "<?xml version='1.0'?>\n"
        "<model>\n"
        " <geometry>\n"
        "  <surface id=\"1\" type=\"z-cylinder\" coeffs=\"0.0 0.0 3.0\" />\n"
        "  <cell id=\"1\" material=\"1\" region=\"-1\" />\n"
        "  <cell id=\"2\" material=\"void\" region=\"1\" />\n"
        " </geometry>\n"
        " <materials>\n"
        "  <material id=\"1\"><density value=\"1.0\" />"
        "  <nuclide name=\"H1\" ao=\"1.0\" /></material>\n"
        " </materials>\n"
        "</model>\n";
    openmc_model_t* omc = openmc_load_string(xml, strlen(xml));
    if (!omc) SKIP("OpenMC string parse not supported");
    alea_system_t* sys = omc->sys;

    alea_build_universe_index(sys);
    ASSERT_EQ(alea_material_at(sys, 0, 0, 0), 1);
    ASSERT_EQ(alea_material_at(sys, 5, 0, 0), 0);

    openmc_model_destroy(omc);
}

TEST(openmc_plane) {
    const char* xml =
        "<?xml version='1.0'?>\n"
        "<model>\n"
        " <geometry>\n"
        "  <surface id=\"1\" type=\"z-plane\" coeffs=\"5.0\" />\n"
        "  <cell id=\"1\" material=\"1\" region=\"-1\" />\n"
        "  <cell id=\"2\" material=\"void\" region=\"1\" />\n"
        " </geometry>\n"
        " <materials>\n"
        "  <material id=\"1\"><density value=\"1.0\" />"
        "  <nuclide name=\"H1\" ao=\"1.0\" /></material>\n"
        " </materials>\n"
        "</model>\n";
    openmc_model_t* omc = openmc_load_string(xml, strlen(xml));
    if (!omc) SKIP("OpenMC string parse not supported");
    alea_system_t* sys = omc->sys;

    alea_build_universe_index(sys);
    ASSERT_EQ(alea_material_at(sys, 0, 0, 3), 1);
    ASSERT_EQ(alea_material_at(sys, 0, 0, 7), 0);

    openmc_model_destroy(omc);
}

TEST(openmc_region_intersection) {
    /* Region with intersection: -1 -2 (inside both surfaces) */
    const char* xml =
        "<?xml version='1.0'?>\n"
        "<model>\n"
        " <geometry>\n"
        "  <surface id=\"1\" type=\"sphere\" coeffs=\"0.0 0.0 0.0 5.0\" />\n"
        "  <surface id=\"2\" type=\"z-plane\" coeffs=\"0.0\" />\n"
        "  <cell id=\"1\" material=\"1\" region=\"-1 -2\" />\n"
        "  <cell id=\"2\" material=\"void\" region=\"1 | 2\" />\n"
        " </geometry>\n"
        " <materials>\n"
        "  <material id=\"1\"><density value=\"1.0\" />"
        "  <nuclide name=\"H1\" ao=\"1.0\" /></material>\n"
        " </materials>\n"
        "</model>\n";
    openmc_model_t* omc = openmc_load_string(xml, strlen(xml));
    if (!omc) SKIP("OpenMC string parse not supported");
    alea_system_t* sys = omc->sys;

    alea_build_universe_index(sys);
    /* Below z=0 and inside sphere → mat 1 */
    ASSERT_EQ(alea_material_at(sys, 0, 0, -3), 1);
    /* Above z=0 → void */
    ASSERT_EQ(alea_material_at(sys, 0, 0, 3), 0);

    openmc_model_destroy(omc);
}

TEST(openmc_region_union) {
    /* Region with union: (-1 | -2) */
    const char* xml =
        "<?xml version='1.0'?>\n"
        "<model>\n"
        " <geometry>\n"
        "  <surface id=\"1\" type=\"sphere\" coeffs=\"-5.0 0.0 0.0 3.0\" />\n"
        "  <surface id=\"2\" type=\"sphere\" coeffs=\"5.0 0.0 0.0 3.0\" />\n"
        "  <cell id=\"1\" material=\"1\" region=\"(-1 | -2)\" />\n"
        "  <cell id=\"2\" material=\"void\" region=\"(1 2)\" />\n"
        " </geometry>\n"
        " <materials>\n"
        "  <material id=\"1\"><density value=\"1.0\" />"
        "  <nuclide name=\"H1\" ao=\"1.0\" /></material>\n"
        " </materials>\n"
        "</model>\n";
    openmc_model_t* omc = openmc_load_string(xml, strlen(xml));
    if (!omc) SKIP("OpenMC string parse not supported");
    alea_system_t* sys = omc->sys;

    alea_build_universe_index(sys);
    ASSERT_EQ(alea_material_at(sys, -5, 0, 0), 1);
    ASSERT_EQ(alea_material_at(sys, 5, 0, 0), 1);
    ASSERT_EQ(alea_material_at(sys, 0, 0, 0), 0);

    openmc_model_destroy(omc);
}

TEST(openmc_complement) {
    /* Region with complement: ~(-1) = NOT(inside sphere) = outside sphere
     * Note: ~1 would mean NOT(outside surface 1) = inside, which equals -1.
     * To get a true complement, use ~(-1) or just 1 (positive half-space). */
    const char* xml =
        "<?xml version='1.0'?>\n"
        "<model>\n"
        " <geometry>\n"
        "  <surface id=\"1\" type=\"sphere\" coeffs=\"0.0 0.0 0.0 5.0\" />\n"
        "  <cell id=\"1\" material=\"1\" region=\"-1\" />\n"
        "  <cell id=\"2\" material=\"2\" region=\"~(-1)\" />\n"
        " </geometry>\n"
        " <materials>\n"
        "  <material id=\"1\"><density value=\"1.0\" />"
        "  <nuclide name=\"H1\" ao=\"1.0\" /></material>\n"
        "  <material id=\"2\"><density value=\"1.0\" />"
        "  <nuclide name=\"H1\" ao=\"1.0\" /></material>\n"
        " </materials>\n"
        "</model>\n";
    openmc_model_t* omc = openmc_load_string(xml, strlen(xml));
    if (!omc) SKIP("OpenMC string parse not supported");
    alea_system_t* sys = omc->sys;

    alea_build_universe_index(sys);
    ASSERT_EQ(alea_material_at(sys, 0, 0, 0), 1);
    /* Complement of inside sphere → outside sphere → mat 2 */
    ASSERT_EQ(alea_material_at(sys, 10, 0, 0), 2);

    openmc_model_destroy(omc);
}

/* ========================================================================= */
/* Fill/universe tests                                                       */
/* ========================================================================= */

TEST(openmc_fill_universe) {
    const char* xml =
        "<?xml version='1.0'?>\n"
        "<model>\n"
        " <geometry>\n"
        "  <surface id=\"1\" type=\"sphere\" coeffs=\"0.0 0.0 0.0 10.0\" />\n"
        "  <surface id=\"2\" type=\"sphere\" coeffs=\"0.0 0.0 0.0 3.0\" />\n"
        "  <cell id=\"1\" fill=\"2\" region=\"-1\" />\n"
        "  <cell id=\"10\" universe=\"2\" material=\"1\" region=\"-2\" />\n"
        "  <cell id=\"11\" universe=\"2\" material=\"void\" region=\"2\" />\n"
        " </geometry>\n"
        " <materials>\n"
        "  <material id=\"1\"><density value=\"1.0\" />"
        "  <nuclide name=\"H1\" ao=\"1.0\" /></material>\n"
        " </materials>\n"
        "</model>\n";
    openmc_model_t* omc = openmc_load_string(xml, strlen(xml));
    if (!omc) SKIP("OpenMC string parse not supported");
    alea_system_t* sys = omc->sys;

    alea_build_universe_index(sys);
    /* Verify cell 1 has fill_universe=2 */
    alea_cell_info_t info;
    ASSERT_EQ(alea_cell_find_info(sys, 1, &info), 0);
    ASSERT_EQ(info.fill_universe, 2);
    /* Cell 10 in universe 2 has material 1 */
    ASSERT_EQ(alea_cell_find_info(sys, 10, &info), 0);
    ASSERT_EQ(info.material_id, 1);
    /* Lazy resolution at origin → material 1 (inside inner sphere) */
    int cell_id, material;
    ASSERT_EQ(alea_find_cell_lazy(sys, 0, 0, 0, &cell_id, &material, NULL), 0);
    ASSERT_EQ(material, 1);

    openmc_model_destroy(omc);
}

/* ========================================================================= */
/* Export roundtrip                                                          */
/* ========================================================================= */

/* ========================================================================= */
/* 1-sheet cone expansion (MCNP → OpenMC roundtrip)                          */
/* ========================================================================= */

/*
 * Helper: parse MCNP string, export to OpenMC XML, re-parse, verify that
 * containment at test points matches the original.  Returns 1 on success.
 */
static int mcnp_to_openmc_roundtrip(const char* mcnp_input, const char* label,
                                     double x_in, double y_in, double z_in,
                                     double x_out, double y_out, double z_out) {
    /* Parse MCNP */
    mcnp_model_t* model1 = mcnp_load_string(mcnp_input, strlen(mcnp_input));
    if (!model1) {
        printf("    FAIL: %s: MCNP parse failed: %s\n", label, alea_error());
        return 0;
    }
    alea_build_universe_index(model1->sys);

    int mat_in  = alea_material_at(model1->sys, x_in, y_in, z_in);
    int mat_out = alea_material_at(model1->sys, x_out, y_out, z_out);

    /* Export to OpenMC */
    const char* tmpfile = "test_cone_sheet_rt_tmp.xml";
    int rc = openmc_export_system(model1->sys, tmpfile);
    if (rc != 0) {
        printf("    FAIL: %s: OpenMC export failed: %s\n", label, alea_error());
        mcnp_model_destroy(model1);
        return 0;
    }
    mcnp_model_destroy(model1);

    /* Re-parse OpenMC XML */
    openmc_model_t* omc2 = openmc_load(tmpfile);
    if (!omc2) {
        printf("    FAIL: %s: OpenMC re-parse failed: %s\n", label, alea_error());
        remove(tmpfile);
        return 0;
    }
    alea_build_universe_index(omc2->sys);

    int mat_in2  = alea_material_at(omc2->sys, x_in, y_in, z_in);
    int mat_out2 = alea_material_at(omc2->sys, x_out, y_out, z_out);

    openmc_model_destroy(omc2);
    remove(tmpfile);

    if (mat_in != mat_in2) {
        printf("    FAIL: %s: inside point mat mismatch: %d vs %d\n",
               label, mat_in, mat_in2);
        return 0;
    }
    if (mat_out != mat_out2) {
        printf("    FAIL: %s: outside point mat mismatch: %d vs %d\n",
               label, mat_out, mat_out2);
        return 0;
    }
    return 1;
}

TEST(openmc_1sheet_cone_kz_pos) {
    /* KZ apex=0 t²=1 sheet=+1  →  only z>0 nappe
       Cell 1: inside cone AND z>0 (material 1)
       Cell 2: everything else (void)
       Inside point:  (0.5, 0, 2)  → on positive nappe, dist from axis = 0.5 < cone radius 2
       Outside point: (0.5, 0, -2) → on negative nappe, should NOT be inside */
    const char* input =
        "Test KZ sheet +1\n"
        "1 1 -1.0 -1 2\n"
        "2 0 1 : -2\n"
        "\n"
        "1 KZ 0.0 1.0 1\n"
        "2 PZ 0.0\n"
        "\n"
        "M1 92235.80c 1.0\n";
    ASSERT(mcnp_to_openmc_roundtrip(input, "KZ+1", 0.5, 0, 2.0, 0.5, 0, -2.0));
}

TEST(openmc_1sheet_cone_kz_neg) {
    /* KZ apex=0 t²=1 sheet=-1  →  only z<0 nappe
       Inside point:  (0.5, 0, -2) → on negative nappe
       Outside point: (0.5, 0, 2)  → on positive nappe */
    const char* input =
        "Test KZ sheet -1\n"
        "1 1 -1.0 -1 2\n"
        "2 0 1 : -2\n"
        "\n"
        "1 KZ 0.0 1.0 -1\n"
        "2 PZ 0.0\n"
        "\n"
        "M1 92235.80c 1.0\n";
    ASSERT(mcnp_to_openmc_roundtrip(input, "KZ-1", 0.5, 0, -2.0, 0.5, 0, 2.0));
}

TEST(openmc_1sheet_cone_kx_pos) {
    /* KX apex=0 t²=1 sheet=+1  →  only x>0 nappe
       Inside point:  (2, 0.5, 0) → positive nappe
       Outside point: (-2, 0.5, 0) → negative nappe */
    const char* input =
        "Test KX sheet +1\n"
        "1 1 -1.0 -1 2\n"
        "2 0 1 : -2\n"
        "\n"
        "1 KX 0.0 1.0 1\n"
        "2 PX 0.0\n"
        "\n"
        "M1 92235.80c 1.0\n";
    ASSERT(mcnp_to_openmc_roundtrip(input, "KX+1", 2.0, 0.5, 0, -2.0, 0.5, 0));
}

TEST(openmc_1sheet_cone_ky_neg) {
    /* KY apex=5 t²=1 sheet=-1  →  only y<5 nappe
       Inside point:  (0.5, 3, 0) → y=3 < 5, dist from axis at y=3 is 0.5, cone radius = 2
       Outside point: (0.5, 7, 0) → y=7 > 5, wrong nappe */
    const char* input =
        "Test KY sheet -1\n"
        "1 1 -1.0 -1 2\n"
        "2 0 1 : -2\n"
        "\n"
        "1 K/Y 0.0 5.0 0.0 1.0 -1\n"
        "2 PY 5.0\n"
        "\n"
        "M1 92235.80c 1.0\n";
    ASSERT(mcnp_to_openmc_roundtrip(input, "KY-1", 0.5, 3.0, 0, 0.5, 7.0, 0));
}

TEST(openmc_union_intersection_parens) {
    /* Test that (A∩B) | (C∩D) intersection branches are parenthesized in OpenMC.
       Without parens, OpenMC misparsed "A B | C D" inside a larger intersection.

       Geometry: two boxes defined via plane intersections, unioned, then
       intersected with a z-slab.
         Cell 1: ((-1 2) | (-3 4)) -5 6   = (left_box | right_box) ∩ z_slab
         Cell 2: everything else

       Surfaces:
         1: PX -2     2: PX -1     → left box:  -2 < x < -1
         3: PX  1     4: PX  2     → right box:  1 < x <  2
         5: PZ  5     6: PZ -5     → z slab:   -5 < z <  5
         7: SO 100                 → outer sphere

       Inside left box:  (-1.5, 0, 0) → mat 1
       Inside right box: ( 1.5, 0, 0) → mat 1
       Between boxes:    ( 0,   0, 0) → void (outside both boxes)
       Above z-slab:     (-1.5, 0, 8) → void (outside z-slab)
    */
    const char* input =
        "Union parens test\n"
        "1 1 -1.0 (-1 2 : -3 4) -5 6\n"
        "2 0 (1:-2) (3:-4) : 5 : -6\n"
        "\n"
        "1 PX -2.0\n"
        "2 PX -1.0\n"
        "3 PX  1.0\n"
        "4 PX  2.0\n"
        "5 PZ  5.0\n"
        "6 PZ -5.0\n"
        "\n"
        "M1 92235.80c 1.0\n";

    /* Roundtrip through OpenMC: parse MCNP → export OpenMC → re-parse */
    ASSERT(mcnp_to_openmc_roundtrip(input, "union_parens(left_box)",
           -1.5, 0, 0,   /* inside left box → mat 1 */
            0,   0, 0));  /* between boxes → void */

    /* Also verify right box and z-slab boundary */
    ASSERT(mcnp_to_openmc_roundtrip(input, "union_parens(right_box)",
            1.5, 0, 0,    /* inside right box → mat 1 */
           -1.5, 0, 8));  /* above z-slab → void */
}

TEST(openmc_export_roundtrip) {
    openmc_model_t* omc = openmc_load("tests/data/openmc_simple.xml");
    if (!omc) SKIP("Test data file not found");
    alea_system_t* sys = omc->sys;

    alea_build_universe_index(sys);
    int mat_origin = alea_material_at(sys, 0, 0, 0);

    /* Export to MCNP */
    const char* tmpfile = "test_openmc_rt_tmp.mcnp";
    int rc = mcnp_export_system(sys, tmpfile);
    ASSERT_EQ(rc, 0);
    openmc_model_destroy(omc);

    /* Re-parse as MCNP */
    mcnp_model_t* model2 = mcnp_load(tmpfile);
    ASSERT_NOT_NULL(model2);
    alea_build_universe_index(model2->sys);

    int mat_origin2 = alea_material_at(model2->sys, 0, 0, 0);
    ASSERT_EQ(mat_origin, mat_origin2);

    mcnp_model_destroy(model2);
    remove(tmpfile);
}

TEST(openmc_hex_lattice_roundtrip) {
    openmc_model_t* omc = openmc_load("tests/data/openmc_hex_lattice.xml");
    if (!omc) SKIP("Test data file not found");
    alea_system_t* sys = omc->sys;
    alea_build_universe_index(sys);

    /* Verify parse: origin point query */
    int cell_id, mat;
    ASSERT_EQ(alea_find_cell_lazy(sys, 0, 0, 0, &cell_id, &mat, NULL), 0);
    ASSERT_EQ(mat, 1);

    /* Export to OpenMC XML */
    const char* tmpfile = "test_hex_lat_rt_tmp.xml";
    int rc = openmc_export_system(sys, tmpfile);
    ASSERT_EQ(rc, 0);
    openmc_model_destroy(omc);

    /* Re-parse the exported file */
    openmc_model_t* omc2 = openmc_load(tmpfile);
    ASSERT_NOT_NULL(omc2);
    alea_build_universe_index(omc2->sys);

    /* Same point queries should produce the same materials */
    ASSERT_EQ(alea_find_cell_lazy(omc2->sys, 0, 0, 0, &cell_id, &mat, NULL), 0);
    ASSERT_EQ(mat, 1);

    /* Element (1,0): center at (2,0) → univ 3 → mat 3 */
    ASSERT_EQ(alea_find_cell_lazy(omc2->sys, 2, 0, 0, &cell_id, &mat, NULL), 0);
    ASSERT_EQ(mat, 3);

    openmc_model_destroy(omc2);
    remove(tmpfile);

    /* --- MCNP roundtrip: parse OpenMC → export MCNP → re-parse → verify --- */
    openmc_model_t* omc3 = openmc_load("tests/data/openmc_hex_lattice.xml");
    ASSERT_NOT_NULL(omc3);
    alea_build_universe_index(omc3->sys);

    const char* mcnp_file = "test_hex_lat_rt_tmp.mcnp";
    rc = mcnp_export_system(omc3->sys, mcnp_file);
    ASSERT_EQ(rc, 0);
    openmc_model_destroy(omc3);

    mcnp_model_t* model4 = mcnp_load(mcnp_file);
    ASSERT_NOT_NULL(model4);
    alea_build_universe_index(model4->sys);

    /* Same point queries should work after MCNP roundtrip */
    ASSERT_EQ(alea_find_cell_lazy(model4->sys, 0, 0, 0, &cell_id, &mat, NULL), 0);
    ASSERT_EQ(mat, 1);

    ASSERT_EQ(alea_find_cell_lazy(model4->sys, 2, 0, 0, &cell_id, &mat, NULL), 0);
    ASSERT_EQ(mat, 3);

    mcnp_model_destroy(model4);
    remove(mcnp_file);
}

TEST_MAIN()
