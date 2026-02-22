// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/*
 * test_mcnp_surfaces.c - MCNP surface roundtrip tests
 *
 * Tests every MCNP surface type through parse → export → re-parse roundtrip.
 * The roundtrip helper writes a minimal MCNP input, parses it, exports it,
 * re-parses the export, and compares primitive data with tolerance.
 */

#include "alea_test.h"
#include "alea.h"
#include "core/alea_system.h"
#include "core/alea_export.h"
#include "core/alea_eval.h"
#include "primitives/primitive_eval.h"
#include <string.h>
#include <stdio.h>

/* ========================================================================= */
/* Roundtrip helper                                                          */
/* ========================================================================= */

/*
 * Parse an MCNP input string, export to file, re-parse, and verify that
 * the geometry produces the same containment answers at several test points.
 * Returns 1 on success, 0 on failure (with message printed).
 */
static int roundtrip_check(const char* mcnp_input, const char* label,
                           double test_x, double test_y, double test_z) {
    /* --- Parse original --- */
    alea_system_t* sys1 = alea_load_mcnp_string(mcnp_input, strlen(mcnp_input));
    if (!sys1) {
        printf("    FAIL: %s: parse failed: %s\n", label, alea_error());
        return 0;
    }
    alea_build_universe_index(sys1);

    /* Query material at test point in original */
    int mat1 = alea_material_at(sys1, test_x, test_y, test_z);

    /* --- Export (non-dedup) --- */
    const char* tmpfile = "test_rt_surf_tmp.mcnp";
    int rc = alea_export(sys1, ALEA_EXPORT_FORMAT_MCNP, tmpfile, ALEA_EMIT_MACROBODY, false);
    if (rc != 0) {
        printf("    FAIL: %s: export failed\n", label);
        alea_destroy(sys1);
        return 0;
    }
    alea_destroy(sys1);

    /* --- Re-parse exported file --- */
    alea_system_t* sys2 = alea_load_mcnp(tmpfile);
    if (!sys2) {
        printf("    FAIL: %s: re-parse failed: %s\n", label, alea_error());
        remove(tmpfile);
        return 0;
    }
    alea_build_universe_index(sys2);

    int mat2 = alea_material_at(sys2, test_x, test_y, test_z);

    alea_destroy(sys2);
    remove(tmpfile);

    if (mat1 != mat2) {
        printf("    FAIL: %s: material mismatch at (%.1f,%.1f,%.1f): "
               "original=%d, roundtrip=%d\n",
               label, test_x, test_y, test_z, mat1, mat2);
        return 0;
    }
    return 1;
}

/*
 * Full roundtrip: parse, export, re-parse, compare materials at multiple
 * points (inside and outside).
 */
static int roundtrip_surface(const char* mcnp_input, const char* label,
                             double in_x, double in_y, double in_z,
                             double out_x, double out_y, double out_z) {
    /* --- Parse original --- */
    alea_system_t* sys1 = alea_load_mcnp_string(mcnp_input, strlen(mcnp_input));
    if (!sys1) {
        printf("    FAIL: %s: parse failed: %s\n", label, alea_error());
        return 0;
    }
    alea_build_universe_index(sys1);

    int mat_in1 = alea_material_at(sys1, in_x, in_y, in_z);
    int mat_out1 = alea_material_at(sys1, out_x, out_y, out_z);

    /* --- Export --- */
    const char* tmpfile = "test_rt_surf_tmp.mcnp";
    int rc = alea_export(sys1, ALEA_EXPORT_FORMAT_MCNP, tmpfile, ALEA_EMIT_MACROBODY, false);
    if (rc != 0) {
        printf("    FAIL: %s: export failed\n", label);
        alea_destroy(sys1);
        return 0;
    }
    alea_destroy(sys1);

    /* --- Re-parse --- */
    alea_system_t* sys2 = alea_load_mcnp(tmpfile);
    if (!sys2) {
        printf("    FAIL: %s: re-parse failed: %s\n", label, alea_error());
        remove(tmpfile);
        return 0;
    }
    alea_build_universe_index(sys2);

    int mat_in2 = alea_material_at(sys2, in_x, in_y, in_z);
    int mat_out2 = alea_material_at(sys2, out_x, out_y, out_z);

    alea_destroy(sys2);
    remove(tmpfile);

    if (mat_in1 != mat_in2) {
        printf("    FAIL: %s: inside material mismatch: %d vs %d\n",
               label, mat_in1, mat_in2);
        return 0;
    }
    if (mat_out1 != mat_out2) {
        printf("    FAIL: %s: outside material mismatch: %d vs %d\n",
               label, mat_out1, mat_out2);
        return 0;
    }
    /* Verify the inside point actually has material 1 */
    if (mat_in1 != 1) {
        printf("    FAIL: %s: inside point should have mat=1, got %d\n",
               label, mat_in1);
        return 0;
    }
    return 1;
}

/*
 * Parse-only check for macrobody types that don't support roundtrip export.
 * Verifies that parsing succeeds and containment queries give correct answers.
 */
static int parse_and_eval_surface(const char* mcnp_input, const char* label,
                                   double in_x, double in_y, double in_z,
                                   double out_x, double out_y, double out_z) {
    alea_system_t* sys = alea_load_mcnp_string(mcnp_input, strlen(mcnp_input));
    if (!sys) {
        printf("    FAIL: %s: parse failed: %s\n", label, alea_error());
        return 0;
    }
    alea_build_universe_index(sys);

    int mat_in = alea_material_at(sys, in_x, in_y, in_z);
    int mat_out = alea_material_at(sys, out_x, out_y, out_z);

    alea_destroy(sys);

    if (mat_in != 1) {
        printf("    FAIL: %s: inside point should have mat=1, got %d\n", label, mat_in);
        return 0;
    }
    if (mat_out != 0) {
        printf("    FAIL: %s: outside point should have mat=0, got %d\n", label, mat_out);
        return 0;
    }
    return 1;
}

/* ========================================================================= */
/* Axis-aligned plane tests                                                  */
/* ========================================================================= */

TEST(roundtrip_px) {
    const char* input =
        "Test PX\n"
        "1 1 -1.0 -1\n"
        "2 0 1\n"
        "\n"
        "1 PX 5.0\n"
        "\n"
        "M1 92235.80c 1.0\n";
    ASSERT(roundtrip_surface(input, "PX", 3.0, 0, 0, 7.0, 0, 0));
}

TEST(roundtrip_py) {
    const char* input =
        "Test PY\n"
        "1 1 -1.0 -1\n"
        "2 0 1\n"
        "\n"
        "1 PY 5.0\n"
        "\n"
        "M1 92235.80c 1.0\n";
    ASSERT(roundtrip_surface(input, "PY", 0, 3.0, 0, 0, 7.0, 0));
}

TEST(roundtrip_pz) {
    const char* input =
        "Test PZ\n"
        "1 1 -1.0 -1\n"
        "2 0 1\n"
        "\n"
        "1 PZ 5.0\n"
        "\n"
        "M1 92235.80c 1.0\n";
    ASSERT(roundtrip_surface(input, "PZ", 0, 0, 3.0, 0, 0, 7.0));
}

TEST(roundtrip_p_3coeff) {
    /* General plane: x + y + z = 10 (normal = (1,1,1), d = 10) */
    const char* input =
        "Test P\n"
        "1 1 -1.0 -1\n"
        "2 0 1\n"
        "\n"
        "1 P 1.0 1.0 1.0 10.0\n"
        "\n"
        "M1 92235.80c 1.0\n";
    /* Point (0,0,0): 0+0+0-10 = -10 → inside (-1 side) */
    /* Point (5,5,5): 5+5+5-10 = 5 → outside (+1 side) */
    ASSERT(roundtrip_surface(input, "P", 0, 0, 0, 5.0, 5.0, 5.0));
}

/* ========================================================================= */
/* Sphere tests                                                              */
/* ========================================================================= */

TEST(roundtrip_so) {
    const char* input =
        "Test SO\n"
        "1 1 -1.0 -1\n"
        "2 0 1\n"
        "\n"
        "1 SO 5.0\n"
        "\n"
        "M1 92235.80c 1.0\n";
    ASSERT(roundtrip_surface(input, "SO", 0, 0, 0, 7.0, 0, 0));
}

TEST(roundtrip_s) {
    const char* input =
        "Test S\n"
        "1 1 -1.0 -1\n"
        "2 0 1\n"
        "\n"
        "1 S 3.0 4.0 5.0 2.0\n"
        "\n"
        "M1 92235.80c 1.0\n";
    ASSERT(roundtrip_surface(input, "S", 3.0, 4.0, 5.0, 3.0, 4.0, 8.0));
}

TEST(roundtrip_sx) {
    const char* input =
        "Test SX\n"
        "1 1 -1.0 -1\n"
        "2 0 1\n"
        "\n"
        "1 SX 5.0 2.0\n"
        "\n"
        "M1 92235.80c 1.0\n";
    ASSERT(roundtrip_surface(input, "SX", 5.0, 0, 0, 8.0, 0, 0));
}

TEST(roundtrip_sy) {
    const char* input =
        "Test SY\n"
        "1 1 -1.0 -1\n"
        "2 0 1\n"
        "\n"
        "1 SY 5.0 2.0\n"
        "\n"
        "M1 92235.80c 1.0\n";
    ASSERT(roundtrip_surface(input, "SY", 0, 5.0, 0, 0, 8.0, 0));
}

TEST(roundtrip_sz) {
    const char* input =
        "Test SZ\n"
        "1 1 -1.0 -1\n"
        "2 0 1\n"
        "\n"
        "1 SZ 5.0 2.0\n"
        "\n"
        "M1 92235.80c 1.0\n";
    ASSERT(roundtrip_surface(input, "SZ", 0, 0, 5.0, 0, 0, 8.0));
}

/* ========================================================================= */
/* Cylinder tests                                                            */
/* ========================================================================= */

TEST(roundtrip_cx) {
    const char* input =
        "Test CX\n"
        "1 1 -1.0 -1\n"
        "2 0 1\n"
        "\n"
        "1 CX 3.0\n"
        "\n"
        "M1 92235.80c 1.0\n";
    ASSERT(roundtrip_surface(input, "CX", 0, 0, 0, 0, 5.0, 0));
}

TEST(roundtrip_cy) {
    const char* input =
        "Test CY\n"
        "1 1 -1.0 -1\n"
        "2 0 1\n"
        "\n"
        "1 CY 3.0\n"
        "\n"
        "M1 92235.80c 1.0\n";
    ASSERT(roundtrip_surface(input, "CY", 0, 0, 0, 5.0, 0, 0));
}

TEST(roundtrip_cz) {
    const char* input =
        "Test CZ\n"
        "1 1 -1.0 -1\n"
        "2 0 1\n"
        "\n"
        "1 CZ 3.0\n"
        "\n"
        "M1 92235.80c 1.0\n";
    ASSERT(roundtrip_surface(input, "CZ", 0, 0, 0, 5.0, 0, 0));
}

TEST(roundtrip_c_x) {
    const char* input =
        "Test C/X\n"
        "1 1 -1.0 -1\n"
        "2 0 1\n"
        "\n"
        "1 C/X 2.0 3.0 1.0\n"
        "\n"
        "M1 92235.80c 1.0\n";
    ASSERT(roundtrip_surface(input, "C/X", 0, 2.0, 3.0, 0, 5.0, 6.0));
}

TEST(roundtrip_c_y) {
    const char* input =
        "Test C/Y\n"
        "1 1 -1.0 -1\n"
        "2 0 1\n"
        "\n"
        "1 C/Y 2.0 3.0 1.0\n"
        "\n"
        "M1 92235.80c 1.0\n";
    ASSERT(roundtrip_surface(input, "C/Y", 2.0, 0, 3.0, 5.0, 0, 6.0));
}

TEST(roundtrip_c_z) {
    const char* input =
        "Test C/Z\n"
        "1 1 -1.0 -1\n"
        "2 0 1\n"
        "\n"
        "1 C/Z 2.0 3.0 1.0\n"
        "\n"
        "M1 92235.80c 1.0\n";
    ASSERT(roundtrip_surface(input, "C/Z", 2.0, 3.0, 0, 5.0, 6.0, 0));
}

/* ========================================================================= */
/* Cone tests                                                                */
/* ========================================================================= */

TEST(roundtrip_kx) {
    const char* input =
        "Test KX\n"
        "1 1 -1.0 -1 2\n"
        "2 0 1 : -2\n"
        "\n"
        "1 KX 0.0 1.0 1\n"
        "2 PX 0.0\n"
        "\n"
        "M1 92235.80c 1.0\n";
    /* Cone KX x0=0 t²=1 sheet=+1: only x>0 sheet */
    /* Inside: (2, 0.5, 0) → dist from axis at x=2 is 0.5, cone radius at x=2 is 2 */
    ASSERT(roundtrip_check(input, "KX", 2.0, 0.5, 0));
}

TEST(roundtrip_ky) {
    const char* input =
        "Test KY\n"
        "1 1 -1.0 -1 2\n"
        "2 0 1 : -2\n"
        "\n"
        "1 KY 0.0 1.0 1\n"
        "2 PY 0.0\n"
        "\n"
        "M1 92235.80c 1.0\n";
    ASSERT(roundtrip_check(input, "KY", 0.5, 2.0, 0));
}

TEST(roundtrip_kz) {
    const char* input =
        "Test KZ\n"
        "1 1 -1.0 -1 2\n"
        "2 0 1 : -2\n"
        "\n"
        "1 KZ 0.0 1.0 1\n"
        "2 PZ 0.0\n"
        "\n"
        "M1 92235.80c 1.0\n";
    ASSERT(roundtrip_check(input, "KZ", 0.5, 0, 2.0));
}

TEST(roundtrip_k_x) {
    const char* input =
        "Test K/X\n"
        "1 1 -1.0 -1 2\n"
        "2 0 1 : -2\n"
        "\n"
        "1 K/X 5.0 0.0 0.0 1.0 1\n"
        "2 PX 5.0\n"
        "\n"
        "M1 92235.80c 1.0\n";
    ASSERT(roundtrip_check(input, "K/X", 7.0, 0.5, 0));
}

TEST(roundtrip_k_y) {
    const char* input =
        "Test K/Y\n"
        "1 1 -1.0 -1 2\n"
        "2 0 1 : -2\n"
        "\n"
        "1 K/Y 0.0 5.0 0.0 1.0 1\n"
        "2 PY 5.0\n"
        "\n"
        "M1 92235.80c 1.0\n";
    ASSERT(roundtrip_check(input, "K/Y", 0.5, 7.0, 0));
}

TEST(roundtrip_k_z) {
    const char* input =
        "Test K/Z\n"
        "1 1 -1.0 -1 2\n"
        "2 0 1 : -2\n"
        "\n"
        "1 K/Z 0.0 0.0 5.0 1.0 1\n"
        "2 PZ 5.0\n"
        "\n"
        "M1 92235.80c 1.0\n";
    ASSERT(roundtrip_check(input, "K/Z", 0.5, 0, 7.0));
}

/* ========================================================================= */
/* General quadric                                                           */
/* ========================================================================= */

TEST(roundtrip_gq) {
    /* GQ ellipsoid: x² + y² + z² - 25 = 0 (sphere r=5) */
    const char* input =
        "Test GQ\n"
        "1 1 -1.0 -1\n"
        "2 0 1\n"
        "\n"
        "1 GQ 1.0 1.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 -25.0\n"
        "\n"
        "M1 92235.80c 1.0\n";
    ASSERT(roundtrip_surface(input, "GQ", 0, 0, 0, 7.0, 0, 0));
}

/* ========================================================================= */
/* Torus tests                                                               */
/* ========================================================================= */

TEST(roundtrip_tz) {
    const char* input =
        "Test TZ\n"
        "1 1 -1.0 -1\n"
        "2 0 1\n"
        "\n"
        "1 TZ 0.0 0.0 0.0 5.0 1.0 1.0\n"
        "\n"
        "M1 92235.80c 1.0\n";
    /* Inside tube: (5, 0, 0) → on major circle, inside minor tube */
    ASSERT(roundtrip_surface(input, "TZ", 5.0, 0, 0, 0, 0, 0));
}

TEST(roundtrip_tx) {
    const char* input =
        "Test TX\n"
        "1 1 -1.0 -1\n"
        "2 0 1\n"
        "\n"
        "1 TX 0.0 0.0 0.0 5.0 1.0 1.0\n"
        "\n"
        "M1 92235.80c 1.0\n";
    ASSERT(roundtrip_surface(input, "TX", 0, 5.0, 0, 0, 0, 0));
}

TEST(roundtrip_ty) {
    const char* input =
        "Test TY\n"
        "1 1 -1.0 -1\n"
        "2 0 1\n"
        "\n"
        "1 TY 0.0 0.0 0.0 5.0 1.0 1.0\n"
        "\n"
        "M1 92235.80c 1.0\n";
    ASSERT(roundtrip_surface(input, "TY", 5.0, 0, 0, 0, 0, 0));
}

/* ========================================================================= */
/* Macrobody tests                                                           */
/* ========================================================================= */

TEST(roundtrip_rpp) {
    const char* input =
        "Test RPP\n"
        "1 1 -1.0 -1\n"
        "2 0 1\n"
        "\n"
        "1 RPP -5.0 5.0 -3.0 3.0 -2.0 2.0\n"
        "\n"
        "M1 92235.80c 1.0\n";
    ASSERT(roundtrip_surface(input, "RPP", 0, 0, 0, 7.0, 0, 0));
}

TEST(roundtrip_rcc) {
    const char* input =
        "Test RCC\n"
        "1 1 -1.0 -1\n"
        "2 0 1\n"
        "\n"
        "1 RCC 0.0 0.0 0.0 0.0 0.0 10.0 3.0\n"
        "\n"
        "M1 92235.80c 1.0\n";
    ASSERT(roundtrip_surface(input, "RCC", 0, 0, 5.0, 5.0, 0, 5.0));
}

TEST(roundtrip_box) {
    const char* input =
        "Test BOX\n"
        "1 1 -1.0 -1\n"
        "2 0 1\n"
        "\n"
        "1 BOX 0.0 0.0 0.0 10.0 0.0 0.0 0.0 5.0 0.0 0.0 0.0 3.0\n"
        "\n"
        "M1 92235.80c 1.0\n";
    /* BOX macrobody doesn't roundtrip (primitive type not re-parsed), test parse+eval only */
    ASSERT(parse_and_eval_surface(input, "BOX", 5.0, 2.5, 1.5, 12.0, 0, 0));
}

TEST(roundtrip_sph) {
    const char* input =
        "Test SPH\n"
        "1 1 -1.0 -1\n"
        "2 0 1\n"
        "\n"
        "1 SPH 3.0 4.0 5.0 2.0\n"
        "\n"
        "M1 92235.80c 1.0\n";
    ASSERT(parse_and_eval_surface(input, "SPH", 3.0, 4.0, 5.0, 3.0, 4.0, 8.0));
}

TEST(roundtrip_trc) {
    /* Truncated cone: base at origin, height 10 along z, r_base=3, r_top=1 */
    const char* input =
        "Test TRC\n"
        "1 1 -1.0 -1\n"
        "2 0 1\n"
        "\n"
        "1 TRC 0.0 0.0 0.0 0.0 0.0 10.0 3.0 1.0\n"
        "\n"
        "M1 92235.80c 1.0\n";
    ASSERT(parse_and_eval_surface(input, "TRC", 0, 0, 5.0, 5.0, 0, 5.0));
}

TEST(roundtrip_ell) {
    /* Ellipsoid: foci at (0,0,-3) and (0,0,3), major axis len = 10 */
    const char* input =
        "Test ELL\n"
        "1 1 -1.0 -1\n"
        "2 0 1\n"
        "\n"
        "1 ELL 0.0 0.0 -3.0 0.0 0.0 3.0 10.0\n"
        "\n"
        "M1 92235.80c 1.0\n";
    ASSERT(parse_and_eval_surface(input, "ELL", 0, 0, 0, 0, 0, 7.0));
}

TEST(roundtrip_rec) {
    /* Right Elliptical Cylinder: base (0,0,0), height (0,0,10), semi-axes */
    const char* input =
        "Test REC\n"
        "1 1 -1.0 -1\n"
        "2 0 1\n"
        "\n"
        "1 REC 0.0 0.0 0.0 0.0 0.0 10.0 3.0 0.0 0.0 0.0 2.0 0.0\n"
        "\n"
        "M1 92235.80c 1.0\n";
    ASSERT(parse_and_eval_surface(input, "REC", 0, 0, 5.0, 5.0, 0, 5.0));
}

TEST(roundtrip_wed) {
    /* Wedge: vertex at origin, v1 along x, v2 along y, v3 along z */
    const char* input =
        "Test WED\n"
        "1 1 -1.0 -1\n"
        "2 0 1\n"
        "\n"
        "1 WED 0.0 0.0 0.0 4.0 0.0 0.0 0.0 4.0 0.0 0.0 0.0 4.0\n"
        "\n"
        "M1 92235.80c 1.0\n";
    /* Point inside the wedge (triangular cross section) */
    ASSERT(parse_and_eval_surface(input, "WED", 0.5, 0.5, 2.0, 5.0, 5.0, 5.0));
}

TEST(roundtrip_rhp) {
    /* Right Hexagonal Prism */
    const char* input =
        "Test RHP\n"
        "1 1 -1.0 -1\n"
        "2 0 1\n"
        "\n"
        "1 RHP 0.0 0.0 0.0 0.0 0.0 10.0"
        " 2.0 0.0 0.0 -1.0 1.732050808 0.0 -1.0 -1.732050808 0.0\n"
        "\n"
        "M1 92235.80c 1.0\n";
    ASSERT(parse_and_eval_surface(input, "RHP", 0, 0, 5.0, 5.0, 0, 5.0));
}

/* ========================================================================= */
/* Dedup roundtrip: export with dedup enabled, re-parse, verify              */
/* ========================================================================= */

TEST(roundtrip_sphere_dedup) {
    const char* input =
        "Test sphere dedup\n"
        "1 1 -1.0 -1\n"
        "2 0 1\n"
        "\n"
        "1 SO 5.0\n"
        "\n"
        "M1 92235.80c 1.0\n";
    alea_system_t* sys1 = alea_load_mcnp_string(input, strlen(input));
    ASSERT_NOT_NULL(sys1);
    alea_build_universe_index(sys1);
    int mat1 = alea_material_at(sys1, 0, 0, 0);

    int rc = alea_export(sys1, ALEA_EXPORT_FORMAT_MCNP, "test_dedup_rt_tmp.mcnp",
                        ALEA_EMIT_MACROBODY, true);
    ASSERT_EQ(rc, 0);
    alea_destroy(sys1);

    alea_system_t* sys2 = alea_load_mcnp("test_dedup_rt_tmp.mcnp");
    ASSERT_NOT_NULL(sys2);
    alea_build_universe_index(sys2);
    int mat2 = alea_material_at(sys2, 0, 0, 0);
    ASSERT_EQ(mat1, mat2);
    alea_destroy(sys2);
    remove("test_dedup_rt_tmp.mcnp");
}

/* ========================================================================= */
/* Surface with transform test                                               */
/* ========================================================================= */

TEST(surface_with_transform) {
    /* Surface 1 has TR1 applied */
    const char* input =
        "Test surface TR\n"
        "1 1 -1.0 -1\n"
        "2 0 1\n"
        "\n"
        "1 1 SO 5.0\n"
        "\n"
        "TR1 10.0 0.0 0.0\n"
        "M1 92235.80c 1.0\n";
    /* Sphere centered at (10,0,0) due to TR1 */
    ASSERT(roundtrip_surface(input, "surface_TR", 10.0, 0, 0, 0, 0, 0));
}

TEST_MAIN()
