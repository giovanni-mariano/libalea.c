// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

#include "alea_xml.h"
#include "alea_mcnp.h"
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static const char input[] =
    "<?xml version=\"1.0\"?>\n"
    "<alea version=\"1\" length_units=\"cm\" name=\"mixed library model\" title=\"round trip\" comments=\"document comment\">\n"
    "  <materials>\n"
    "    <material id=\"7\" name=\"fuel\" comments=\"material comment\" fraction_basis=\"atom\" density=\"10.4\" density_units=\"g/cm3\">\n"
    "      <nuclide zaid=\"92235\" fraction=\"0.04\" library=\"80c\"/>\n"
    "      <nuclide zaid=\"92238\" fraction=\"0.96\" library=\"32c\"/>\n"
    "    </material>\n"
    "  </materials>\n"
    "  <surfaces>\n"
    "    <surface id=\"100000001\" type=\"sphere\" coeffs=\"0 0 0 5\" boundary=\"vacuum\"/>\n"
    "  </surfaces>\n"
    "  <cells>\n"
    "    <cell id=\"10\" name=\"fuel cell\" comments=\"cell comment\" inline_comment=\"inline\" universe=\"0\" material=\"7\" density=\"10.2\" density_units=\"g/cm3\" region=\"-100000001\" bbox=\"-5 5 -5 5 -5 5\">\n"
    "      <importance particle=\"neutron\" value=\"1\"/>\n"
    "      <importance particle=\"photon\" value=\"0.5\"/>\n"
    "      <parameter name=\"volume\" value=\"523.598775598\"/>\n"
    "      <parameter name=\"nonu\" value=\"1\"/>\n"
    "    </cell>\n"
    "  </cells>\n"
    "</alea>\n";

static alea_model_t* round_trip(const alea_model_t* model) {
    FILE* f = tmpfile();
    assert(f);
    assert(alea_xml_export_stream(model, f) == 0);
    assert(fseek(f, 0, SEEK_END) == 0);
    long size = ftell(f);
    assert(size > 0);
    rewind(f);
    char* text = malloc((size_t)size + 1);
    assert(text);
    assert(fread(text, 1, (size_t)size, f) == (size_t)size);
    text[size] = '\0';
    assert(strstr(text, " bbox=\"") != NULL);
    fclose(f);
    alea_model_t* result = alea_xml_load_string(text, (size_t)size);
    free(text);
    return result;
}

static void check_model(const alea_model_t* model) {
    assert(model);
    assert(strcmp(alea_model_name(model), "mixed library model") == 0);
    assert(strcmp(alea_model_title(model), "round trip") == 0);
    assert(strcmp(alea_model_comments(model), "document comment") == 0);
    const alea_system_t* sys = alea_model_system_const(model);
    assert(alea_material_count(sys) == 1);
    assert(alea_material_nuclide_count(sys, 0) == 2);
    int zaid = 0;
    const char* library = NULL;
    double fraction = 0.0;
    assert(alea_material_nuclide_get(sys, 0, 0, &zaid, &library, &fraction) == 0);
    assert(zaid == 92235 && strcmp(library, ".80c") == 0 && fabs(fraction - 0.04) < 1e-14);
    assert(alea_material_nuclide_get(sys, 0, 1, &zaid, &library, &fraction) == 0);
    assert(zaid == 92238 && strcmp(library, ".32c") == 0 && fabs(fraction - 0.96) < 1e-14);
    assert(alea_cell_count(sys) == 1);
    alea_cell_info_t info;
    assert(alea_cell_get_info(sys, 0, &info) == 0);
    assert(info.cell_id == 10 && info.material_id == 7 && info.is_mass_density);
    assert(fabs(info.density - 10.2) < 1e-14);
    assert(strcmp(info.comments, "cell comment") == 0);
    assert(strcmp(info.inline_comment, "inline") == 0);
    const alea_model_cell_metadata_t* meta = alea_model_cell_metadata(model, 0);
    assert(meta && strcmp(meta->name, "fuel cell") == 0);
    assert(meta->has_importance_neutron && meta->has_importance_photon);
    assert(meta->importance_neutron == 1.0 && meta->importance_photon == 0.5);
    assert((meta->parameter_flags & ALEA_CELL_PARAM_VOLUME) != 0);
    assert((meta->parameter_flags & ALEA_CELL_PARAM_NONU) != 0);
}

int main(void) {
    alea_model_t* first = alea_xml_load_string(input, sizeof(input) - 1);
    if (!first) fprintf(stderr, "%s\n", alea_get_error_detail());
    check_model(first);
    alea_model_t* second = round_trip(first);
    if (!second) fprintf(stderr, "%s\n", alea_get_error_detail());
    check_model(second);
    mcnp_model_t* mcnp = mcnp_model_from_alea_model(second);
    assert(mcnp);
    const mcnp_cell_params_t* params = mcnp_cell_params_const(mcnp, 0);
    assert(params && params->has_imp_n && params->has_imp_p);
    assert(params->has_vol && params->has_nonu);
    alea_model_t* third = mcnp_model_to_alea_model(mcnp);
    assert(third);
    const alea_model_cell_metadata_t* converted = alea_model_cell_metadata(third, 0);
    assert(converted && converted->has_importance_neutron);
    assert((converted->parameter_flags & ALEA_CELL_PARAM_VOLUME) != 0);
    alea_model_destroy(third);
    mcnp_model_destroy(mcnp);

    const char with_dtd[] = "<!DOCTYPE alea><alea version=\"1\" length_units=\"cm\"><cells/></alea>";
    assert(alea_xml_load_string(with_dtd, sizeof(with_dtd) - 1) == NULL);
    alea_model_destroy(second);
    alea_model_destroy(first);
    puts("test_alea_xml: OK");
    return 0;
}
