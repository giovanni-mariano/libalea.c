// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

#include "alea_xml.h"
#include "alea_mcnp.h"
#include "core/alea_system.h"
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static alea_model_t* load_parameters(const char* parameters) {
    char xml[8192];
    int size = snprintf(xml, sizeof(xml),
        "<alea version=\"1\" length_units=\"cm\">"
        "<surfaces><surface id=\"1\" type=\"sphere\" coeffs=\"0 0 0 2\"/></surfaces>"
        "<cells><cell id=\"1\" material=\"0\" region=\"-1\">%s</cell></cells></alea>", parameters);
    assert(size > 0 && (size_t)size < sizeof(xml));
    return alea_xml_load_string(xml, (size_t)size);
}

static char* read_stream(FILE* stream) {
    assert(fseek(stream, 0, SEEK_END) == 0);
    long size = ftell(stream);
    assert(size >= 0);
    char* text = malloc((size_t)size + 1);
    assert(text);
    rewind(stream);
    assert(fread(text, 1, (size_t)size, stream) == (size_t)size);
    text[size] = '\0';
    return text;
}

static char* export_xml(const alea_model_t* model) {
    FILE* stream = tmpfile();
    assert(stream && alea_xml_export_stream(model, stream) == 0);
    char* text = read_stream(stream);
    fclose(stream);
    return text;
}

static char* export_mcnp(const mcnp_model_t* model) {
    FILE* stream = tmpfile();
    assert(stream && mcnp_export_stream(model, stream) == 0);
    char* text = read_stream(stream);
    fclose(stream);
    return text;
}

static void check_defaults(void) {
    const char* inputs[] = {"", "<parameters/>",
        "<parameters><importance particle=\"neutron\" value=\"0\"/></parameters>"};
    for (size_t i = 0; i < sizeof(inputs) / sizeof(inputs[0]); i++) {
        alea_model_t* model = load_parameters(inputs[i]);
        assert(model);
        const alea_model_cell_metadata_t* meta = alea_model_cell_metadata(model, 0);
        assert(meta->importance_neutron == (i == 2 ? 0 : 1));
        assert(meta->importance_photon == 1 && meta->importance_electron == 1);
        assert(meta->has_importance_neutron == (i == 2));
        assert(!meta->has_importance_photon && !meta->has_importance_electron);
        assert(meta->fission_mode == ALEA_FISSION_NORMAL);
        assert(meta->secondary_collision_state == ALEA_SECONDARY_UNCOLLIDED);
        assert(meta->detector_contribution == 1);
        char* xml = export_xml(model);
        assert((strstr(xml, "<parameters>") != NULL) == (i == 2));
        alea_model_t* copy = alea_xml_load_string(xml, 0);
        assert(copy && alea_model_cell_metadata(copy, 0)->importance_electron == 1);
        mcnp_model_t* mcnp = mcnp_model_from_alea_model(copy);
        assert(mcnp);
        char* mcnp_text = export_mcnp(mcnp);
        assert(strstr(mcnp_text, i == 2 ? "IMP:N=0" : "IMP:N=1"));
        assert(strstr(mcnp_text, "IMP:P=1"));
        free(mcnp_text);
        mcnp_model_destroy(mcnp);
        alea_model_destroy(copy);
        free(xml);
        alea_model_destroy(model);
    }
}

static void check_scoped_metadata(const alea_model_t* model) {
    assert(model);
    const alea_model_cell_metadata_t* m = alea_model_cell_metadata(model, 0);
    assert(m->parameter_flags & ALEA_CELL_PARAM_VOLUME);
    assert(m->user_volume == 12.5);
    assert(m->photon_production.mode == ALEA_PHOTON_PRODUCTION_THRESHOLD);
    assert(m->photon_production.weight_basis == ALEA_PHOTON_WEIGHT_SOURCE_RELATIVE);
    assert(m->photon_production.weight_threshold == 0.5);
    assert(m->fission_mode == ALEA_FISSION_CAPTURE_WITHOUT_PHOTONS);
    assert(m->importance_photon == 0 && m->importance_neutron == 1);
    assert(m->energy_cutoff_particles == ((1u << ALEA_PARTICLE_NEUTRON) | (1u << ALEA_PARTICLE_ELECTRON)));
    assert(m->particle_energy_cutoff[ALEA_PARTICLE_NEUTRON] == 0.01);
    assert(m->particle_energy_cutoff[ALEA_PARTICLE_ELECTRON] == 0.002);
    assert(m->secondary_state_particles == ((1u << ALEA_PARTICLE_PHOTON) | (1u << ALEA_PARTICLE_ELECTRON)));
    assert(m->particle_secondary_state[ALEA_PARTICLE_PHOTON] == ALEA_SECONDARY_COLLIDED);
    assert(m->particle_secondary_state[ALEA_PARTICLE_ELECTRON] == ALEA_SECONDARY_UNCOLLIDED);
    assert(m->detector_contribution == 0.75);
    assert(m->detector_probability_count == 2);
    assert(m->detector_probabilities[0].tally == 5 && m->detector_probabilities[0].probability == 0.25);
    assert(m->detector_probabilities[1].tally == 15 && m->detector_probabilities[1].probability == 0);
    assert(m->magnetic_field == 2);
    alea_cell_info_t info;
    assert(alea_cell_get_info(alea_model_system_const(model), 0, &info) == 0);
    assert(info.has_temperature && fabs(info.temperature - 600) < 1e-8);
}

static void check_scoped_round_trip(void) {
    alea_model_t* model = load_parameters(
        "<parameters>"
        "<temperature value=\"600\" units=\"K\"/>"
        "<volume value=\"12.5\"/>"
        "<importance particle=\"photon\" value=\"0\"/>"
        "<photon_production mode=\"threshold\" weight_threshold=\"0.5\" weight_basis=\"source_relative\"/>"
        "<fission_mode value=\"capture_without_photons\"/>"
        "<detector_contribution_probability value=\"0.75\"/>"
        "<detector_contribution_probability tally=\"15\" value=\"0\"/>"
        "<detector_contribution_probability tally=\"5\" value=\"0.25\"/>"
        "<energy_cutoff particle=\"neutron\" value=\"0.01\" units=\"MeV\"/>"
        "<energy_cutoff particle=\"electron\" value=\"0.002\"/>"
        "<secondary_collision_state particle=\"photon\" value=\"collided\"/>"
        "<secondary_collision_state particle=\"electron\" value=\"uncollided\"/>"
        "<magnetic_field ref=\"2\"/>"
        "</parameters>");
    check_scoped_metadata(model);
    char* xml = export_xml(model);
    assert(strstr(xml, "<parameters>") && strstr(xml, "<photon_production"));
    assert(!strstr(xml, "<parameter ") && !strstr(xml, "-1000000"));
    alea_model_t* copy = alea_xml_load_string(xml, 0);
    check_scoped_metadata(copy);
    free(xml);
    mcnp_model_t* mcnp = mcnp_model_from_alea_model(copy);
    assert(mcnp);
    const mcnp_cell_params_t* params = mcnp_cell_params_const(mcnp, 0);
    assert(params->has_pwt && params->pwt == -0.5);
    assert(params->has_nonu && params->nonu == 2);
    char* text = export_mcnp(mcnp);
    assert(strstr(text, "PD0=0.75") && strstr(text, "PD5=0.25") && strstr(text, "PD15=0"));
    assert(strstr(text, "ELPT:N=0.01") && strstr(text, "ELPT:E=0.002"));
    assert(strstr(text, "UNC:P=0") && strstr(text, "UNC:E=1"));
    assert(strstr(text, "TMP="));
    mcnp_model_t* reloaded = mcnp_load_string(text, 0);
    if (!reloaded) fprintf(stderr, "%s\n", alea_get_error_detail());
    assert(reloaded);
    alea_model_t* returned = mcnp_model_to_alea_model(reloaded);
    check_scoped_metadata(returned);
    free(text);
    alea_model_destroy(returned);
    mcnp_model_destroy(reloaded);
    mcnp_model_destroy(mcnp);
    alea_model_destroy(copy);
    alea_model_destroy(model);
}

static void check_photon_and_fission_modes(void) {
    const double weights[] = {1, -1, 0, -1e6, 0.5, -2.5, 0.12345678912345678};
    const char* modes[] = {"threshold", "threshold", "one_per_collision", "off", "threshold", "threshold", "threshold"};
    const char* fission_modes[] = {"capture_with_photons", "normal", "capture_without_photons"};
    for (size_t i = 0; i < sizeof(weights) / sizeof(weights[0]); i++) {
        for (int fission = 0; fission < 3; fission++) {
            char input[512];
            snprintf(input, sizeof(input),
                "Photon and fission modes\n1 0 -1 PWT=%.17g NONU=%d\n\n1 SO 2\n\n",
                weights[i], fission);
            mcnp_model_t* imported = mcnp_load_string(input, 0);
            assert(imported);
            alea_model_t* converted = mcnp_model_to_alea_model(imported);
            assert(converted);
            char* xml = export_xml(converted);
            assert(strstr(xml, modes[i]) && strstr(xml, fission_modes[fission]));
            assert(!strstr(xml, "<parameter "));
            alea_model_t* native = alea_xml_load_string(xml, 0);
            assert(native);
            mcnp_model_t* mcnp = mcnp_model_from_alea_model(native);
            assert(mcnp);
            const mcnp_cell_params_t* params = mcnp_cell_params_const(mcnp, 0);
            assert(params->has_pwt && params->pwt == weights[i]);
            assert(params->has_nonu && params->nonu == fission);
            char* mcnp_text = export_mcnp(mcnp);
            mcnp_model_t* reloaded = mcnp_load_string(mcnp_text, 0);
            assert(reloaded && mcnp_cell_params_const(reloaded, 0)->pwt == weights[i]);
            assert(mcnp_cell_params_const(reloaded, 0)->nonu == fission);
            free(mcnp_text);
            mcnp_model_destroy(reloaded);
            mcnp_model_destroy(mcnp);
            alea_model_destroy(native);
            free(xml);
            alea_model_destroy(converted);
            mcnp_model_destroy(imported);
        }
    }
    /* A valid native threshold that clashes with MCNP's special sentinel fails conversion. */
    alea_model_t* model = load_parameters(
        "<parameters><photon_production mode=\"threshold\" weight_threshold=\"1000000\" weight_basis=\"source_relative\"/></parameters>");
    assert(model);
    char* xml = export_xml(model);
    free(xml);
    assert(!mcnp_model_from_alea_model(model));
    assert(strstr(alea_get_error_detail(), "cannot be represented"));
    alea_model_destroy(model);
}

static void check_all_particle_defaults(void) {
    alea_model_t* model = load_parameters(
        "<parameters>"
        "<importance particle=\"photon\" value=\"0\"/>"
        "<volume value=\"0\"/>"
        "<detector_contribution_probability value=\"0.5\"/>"
        "<energy_cutoff value=\"0.001\"/>"
        "<secondary_collision_state value=\"collided\"/>"
        "<magnetic_field ref=\"0\"/>"
        "</parameters>");
    assert(model);
    char* xml = export_xml(model);
    assert(strstr(xml, "<volume value=\"0\"") && strstr(xml, "<magnetic_field ref=\"0\""));
    assert(strstr(xml, "<detector_contribution_probability") && strstr(xml, "<secondary_collision_state"));
    mcnp_model_t* mcnp = mcnp_model_from_alea_model(model);
    assert(mcnp);
    char* text = export_mcnp(mcnp);
    assert(strstr(text, "ELPT:N=") && strstr(text, "ELPT:P=") && strstr(text, "ELPT:E="));
    assert(strstr(text, "UNC:N=0") && strstr(text, "UNC:P=0") && strstr(text, "UNC:E=0"));
    assert(strstr(text, "PD0=0.5"));
    mcnp_model_t* copy = mcnp_load_string(text, 0);
    assert(copy);
    mcnp_model_destroy(copy);
    free(text);
    mcnp_model_destroy(mcnp);
    free(xml);
    alea_model_destroy(model);
}

static void check_invalid_parameters(void) {
    const char* invalid[] = {
        "<parameters/><parameters/>",
        "<parameters value=\"1\"/>",
        "<parameters>text</parameters>",
        "<volume value=\"1\"/>",
        "<parameters><parameter name=\"volume\" value=\"1\"/></parameters>",
        "<parameter name=\"volume\" value=\"1\"/>",
        "<parameter name=\"pwt\" value=\"-1\"/>",
        "<parameter name=\"nonu\" value=\"1\"/>",
        "<parameter name=\"pd\" value=\"0.5\"/>",
        "<parameter name=\"elpt\" value=\"0.001\"/>",
        "<parameter name=\"unc\" value=\"1\"/>",
        "<parameter name=\"bflcl\" value=\"0\"/>",
        "<importance particle=\"neutron\" value=\"1\"/>",
        "<parameters><unknown/></parameters>",
        "<parameters><volume value=\"1\"><volume/></volume></parameters>",
        "<parameters><volume value=\"-1\"/></parameters>",
        "<parameters><volume value=\"nan\"/></parameters>",
        "<parameters><volume/></parameters>",
        "<parameters><volume value=\" \"/></parameters>",
        "<parameters><magnetic_field ref=\" \"/></parameters>",
        "<parameters><volume value=\"1\"/><volume value=\"2\"/></parameters>",
        "<parameters><importance particle=\"neutron\" value=\"1\"/><importance particle=\"neutron\" value=\"0\"/></parameters>",
        "<parameters><importance particle=\"unknown\" value=\"1\"/></parameters>",
        "<parameters><importance particle=\"neutron\" value=\"-1\"/></parameters>",
        "<parameters><photon_production mode=\"threshold\" weight_threshold=\"0\" weight_basis=\"absolute\"/></parameters>",
        "<parameters><photon_production mode=\"threshold\" weight_threshold=\"1\"/></parameters>",
        "<parameters><photon_production mode=\"off\" weight_threshold=\"1\"/></parameters>",
        "<parameters><photon_production mode=\"unknown\"/></parameters>",
        "<parameters><fission_mode value=\"1\"/></parameters>",
        "<parameters><fission_mode value=\"unknown\"/></parameters>",
        "<parameters><detector_contribution_probability value=\"1.1\"/></parameters>",
        "<parameters><detector_contribution_probability tally=\"-5\" value=\"0.5\"/></parameters>",
        "<parameters><detector_contribution_probability tally=\"5\" value=\"0.5\"/><detector_contribution_probability tally=\"5\" value=\"0.25\"/></parameters>",
        "<parameters><energy_cutoff particle=\"unknown\" value=\"1\"/></parameters>",
        "<parameters><energy_cutoff value=\"1\" units=\"eV\"/></parameters>",
        "<parameters><energy_cutoff particle=\"neutron\" value=\"1\"/><energy_cutoff particle=\"neutron\" value=\"2\"/></parameters>",
        "<parameters><secondary_collision_state value=\"1\"/></parameters>",
        "<parameters><secondary_collision_state particle=\"photon\" value=\"collided\"/><secondary_collision_state particle=\"photon\" value=\"uncollided\"/></parameters>",
        "<parameters><magnetic_field ref=\"-1\"/></parameters>",
        "<parameters><temperature value=\"600\" units=\"MeV\"/></parameters>",
        "<parameters><temperature value=\"600\"/><temperature value=\"300\"/></parameters>"
    };
    for (size_t i = 0; i < sizeof(invalid) / sizeof(invalid[0]); i++) {
        alea_model_t* model = load_parameters(invalid[i]);
        if (model) fprintf(stderr, "unexpectedly accepted: %s\n", invalid[i]);
        assert(!model);
    }
    const char old_temperature[] =
        "<alea version=\"1\" length_units=\"cm\">"
        "<surfaces><surface id=\"1\" type=\"sphere\" coeffs=\"0 0 0 2\"/></surfaces>"
        "<cells><cell id=\"1\" region=\"-1\" temperature=\"600\"/></cells></alea>";
    assert(!alea_xml_load_string(old_temperature, 0));
    alea_model_t* model = load_parameters("<parameters><volume value=\"1\"/></parameters>");
    assert(model);
    alea_model_cell_metadata_mut(model, 0)->user_volume = NAN;
    FILE* stream = tmpfile();
    assert(stream && alea_xml_export_stream(model, stream) == -1);
    assert(ftell(stream) == 0); /* Validation precedes output. */
    fclose(stream);
    alea_model_destroy(model);
}

static void check_copy_and_removal(void) {
    alea_model_t* model = load_parameters(
        "<parameters><detector_contribution_probability tally=\"5\" value=\"0.5\"/></parameters>");
    assert(model);
    alea_system_t* sys = alea_model_system(model);
    int copy = alea_add_cell(sys, 2, sys->cells.data[0].root_node_id, ALEA_MATERIAL_VOID, 0, 0);
    assert(copy == 1);
    sys->on_cell_copied(sys->cell_hook_userdata, 1, 0);
    assert(alea_model_cell_metadata(model, 1)->detector_probabilities !=
           alea_model_cell_metadata(model, 0)->detector_probabilities);
    assert(alea_model_cell_set_detector_probability(model, 1, 5, 0.25) == 0);
    assert(alea_model_cell_metadata(model, 0)->detector_probabilities[0].probability == 0.5);
    assert(alea_cell_remove(sys, 0) == 0);
    assert(alea_model_cell_metadata_count(model) == 1);
    assert(alea_model_cell_metadata(model, 0)->detector_probabilities[0].probability == 0.25);
    mcnp_model_t* mcnp = mcnp_model_from_alea_model(model);
    assert(mcnp);
    copy = alea_add_cell(mcnp->sys, 3, mcnp->sys->cells.data[0].root_node_id, ALEA_MATERIAL_VOID, 0, 0);
    assert(copy == 1);
    mcnp->sys->on_cell_copied(mcnp->sys->cell_hook_userdata, 1, 0);
    mcnp_cell_params(mcnp, 1)->scoped.detector_probabilities[0].probability = 0.125;
    assert(mcnp_cell_params_const(mcnp, 0)->scoped.detector_probabilities[0].probability == 0.25);
    assert(alea_cell_remove(mcnp->sys, 0) == 0);
    assert(mcnp_cell_params_const(mcnp, 0)->scoped.detector_probabilities[0].probability == 0.125);
    mcnp_model_destroy(mcnp);
    alea_model_destroy(model);
}

static void check_mcnp_particle_lists_and_like(void) {
    const char* input = "Scoped inheritance\n"
        "1 0 -1 PWT=-1 NONU=0 VOL=2 TMP=5.1703999572e-8 ELPT:N,P=0.125 UNC:P,E=0 PD5=0.5 PD15=0.25\n"
        "2 LIKE 1 BUT ELPT:P=0.25 UNC:E=1 PD5=0.125\n"
        "\n1 SO 2\n\n";
    mcnp_model_t* model = mcnp_load_string(input, 0);
    if (!model) fprintf(stderr, "%s\n", alea_get_error_detail());
    assert(model);
    const mcnp_cell_params_t* first = mcnp_cell_params_const(model, 0);
    const mcnp_cell_params_t* second = mcnp_cell_params_const(model, 1);
    assert(first->scoped.energy_cutoff[ALEA_PARTICLE_PHOTON] == 0.125);
    assert(second->scoped.energy_cutoff[ALEA_PARTICLE_NEUTRON] == 0.125);
    assert(second->scoped.energy_cutoff[ALEA_PARTICLE_PHOTON] == 0.25);
    assert(second->scoped.secondary_state[ALEA_PARTICLE_PHOTON] == 0);
    assert(second->scoped.secondary_state[ALEA_PARTICLE_ELECTRON] == 1);
    assert(second->scoped.detector_probability_count == 2);
    assert(second->scoped.detector_probabilities[0].probability == 0.125);
    assert(second->scoped.detector_probabilities[1].probability == 0.25);
    assert(second->has_nonu && second->nonu == 0 && second->has_vol && second->vol == 2);
    alea_cell_info_t info;
    assert(alea_cell_get_info(model->sys, 1, &info) == 0);
    assert(info.has_temperature && fabs(info.temperature - 600) < 1e-8);
    char* text = export_mcnp(model);
    mcnp_model_t* copy = mcnp_load_string(text, 0);
    assert(copy);
    assert(mcnp_cell_params_const(copy, 1)->scoped.detector_probability_count == 2);
    free(text);
    mcnp_model_destroy(copy);
    mcnp_model_destroy(model);
}

int main(void) {
    check_defaults();
    check_scoped_round_trip();
    check_photon_and_fission_modes();
    check_all_particle_defaults();
    check_invalid_parameters();
    check_copy_and_removal();
    check_mcnp_particle_lists_and_like();
    puts("test_alea_cell_parameters: OK");
    return 0;
}
