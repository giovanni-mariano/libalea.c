// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/**
 * @file nuc_inventory.c
 * @brief Report decoded and sampleable capabilities for every xsdir entry.
 */

#include "alea_nucdata.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static const char* table_type_name(alea_nuc_table_type_t type) {
    switch (type) {
    case ALEA_NUC_TABLE_CONTINUOUS_NEUTRON: return "neutron";
    case ALEA_NUC_TABLE_PHOTOATOMIC: return "photoatomic";
    case ALEA_NUC_TABLE_PHOTONUCLEAR: return "photonuclear";
    case ALEA_NUC_TABLE_THERMAL_SAB: return "thermal";
    case ALEA_NUC_TABLE_ELECTRON: return "electron";
    }
    return "unknown";
}

static const char* issue_name(alea_nuc_prepare_issue_t issue) {
    switch (issue) {
    case ALEA_NUC_PREP_OK: return "ok";
    case ALEA_NUC_PREP_EMPTY_MATERIAL: return "empty-material";
    case ALEA_NUC_PREP_UNSUPPORTED_CAPABILITY: return "unsupported-capability";
    case ALEA_NUC_PREP_UNSUPPORTED_PARTICLE: return "unsupported-particle";
    case ALEA_NUC_PREP_UNSUPPORTED_URR: return "unsupported-urr";
    case ALEA_NUC_PREP_UNSUPPORTED_REACTION: return "unsupported-reaction";
    case ALEA_NUC_PREP_INVALID_ENERGY_DISTRIBUTION:
        return "invalid-energy-distribution";
    case ALEA_NUC_PREP_INVALID_ANGULAR: return "invalid-angular";
    case ALEA_NUC_PREP_INVALID_CROSS_SECTIONS: return "invalid-cross-sections";
    case ALEA_NUC_PREP_INVALID_THERMAL_ASSOCIATION:
        return "invalid-thermal-association";
    }
    return "unknown";
}

static void append_capability(char* text, size_t capacity, const char* name) {
    size_t used = strlen(text);
    if (used && used + 1 < capacity) text[used++] = ',';
    if (used < capacity)
        snprintf(text + used, capacity - used, "%s", name);
}

static void format_capabilities(uint32_t capabilities, char text[192]) {
    static const struct {
        uint32_t bit;
        const char* name;
    } names[] = {
        {ALEA_NUC_CAP_STATIONARY_ELASTIC, "stationary-elastic"},
        {ALEA_NUC_CAP_ABSORPTION, "absorption"},
        {ALEA_NUC_CAP_FREE_GAS, "free-gas"},
        {ALEA_NUC_CAP_THERMAL_SAB, "thermal-sab"},
        {ALEA_NUC_CAP_FISSION, "fission"},
        {ALEA_NUC_CAP_PHOTON, "photon"},
        {ALEA_NUC_CAP_URR, "urr"},
        {ALEA_NUC_CAP_NEUTRON_EMISSION, "neutron-emission"},
        {ALEA_NUC_CAP_DELAYED_NEUTRON, "delayed-neutron"},
        {ALEA_NUC_CAP_PHOTON_PRODUCTION, "photon-production"},
    };
    text[0] = '\0';
    for (size_t i = 0; i < sizeof(names) / sizeof(names[0]); i++)
        if (capabilities & names[i].bit)
            append_capability(text, 192, names[i].name);
    if (!text[0]) snprintf(text, 192, "none");
}

static int selected(int argc, char** argv, const char* zaid) {
    if (argc == 2) return 1;
    for (int i = 2; i < argc; i++)
        if (strcmp(argv[i], zaid) == 0) return 1;
    return 0;
}

static void print_row(const alea_nuc_xsdir_entry_t* entry,
                      const char* decoded, const char* cross_sections,
                      const char* collision, uint32_t capabilities,
                      const char* issue, int mt, int law, const char* detail) {
    char names[192];
    format_capabilities(capabilities, names);
    printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%d\t%d\t%s\n",
           entry->zaid, table_type_name(entry->type), decoded, cross_sections,
           collision, names, issue, mt, law, detail ? detail : "");
}

int main(int argc, char** argv) {
    if (argc < 2) {
        fprintf(stderr, "Usage: %s XSDIR [ZAID ...]\n", argv[0]);
        return 2;
    }

    alea_nuc_xsdir_t* xsdir = alea_nuc_xsdir_load(argv[1]);
    if (!xsdir) {
        fprintf(stderr, "Cannot load xsdir: %s\n", argv[1]);
        return 1;
    }
    for (int i = 2; i < argc; i++) {
        if (!alea_nuc_xsdir_find(xsdir, argv[i])) {
            fprintf(stderr, "ZAID is not present in xsdir: %s\n", argv[i]);
            alea_nuc_xsdir_free(xsdir);
            return 3;
        }
    }

    printf("zaid\ttype\tdecoded\tcross_sections\tcollision\tcapabilities"
           "\tissue\tmt\tlaw\tdetail\n");
    for (size_t i = 0; i < xsdir->count; i++) {
        const alea_nuc_xsdir_entry_t* entry = &xsdir->entries[i];
        if (!selected(argc, argv, entry->zaid)) continue;

        if (entry->type == ALEA_NUC_TABLE_THERMAL_SAB) {
            alea_nuc_thermal_t* thermal =
                alea_nuc_load_thermal(xsdir, entry->zaid);
            if (thermal) {
                print_row(entry, "yes", "yes", "yes",
                          ALEA_NUC_CAP_THERMAL_SAB, "ok", 0, 0, "");
                alea_nuc_thermal_free(thermal);
            } else {
                print_row(entry, "no", "no", "no", 0, "decode-failed",
                          0, 0, "thermal table could not be decoded");
            }
            continue;
        }

        if (entry->type != ALEA_NUC_TABLE_CONTINUOUS_NEUTRON &&
            entry->type != ALEA_NUC_TABLE_PHOTOATOMIC) {
            print_row(entry, "no", "no", "no", 0, "unsupported-table-type",
                      0, 0, "no decoder is implemented for this table type");
            continue;
        }

        alea_nuc_nuclide_t* nuc =
            alea_nuc_load_nuclide(xsdir, entry->zaid);
        if (!nuc) {
            print_row(entry, "no", "no", "no", 0, "decode-failed", 0, 0,
                      "table could not be decoded");
            continue;
        }

        alea_nuc_capability_report_t report;
        alea_error_t status = alea_nuc_capabilities(nuc, &report);
        int collision_ready = status == ALEA_OK;
        uint32_t capabilities = report.available_capabilities;
        if (entry->type == ALEA_NUC_TABLE_CONTINUOUS_NEUTRON) {
            /* Photon production is opt-in. Report its defect when present,
             * but do not let it hide otherwise valid neutron collisions. */
            alea_nuc_mat_component_t component = {nuc, 1.0};
            alea_nuc_material_t material = {&component, 1, 1};
            alea_nuc_prepare_requirements_t requirements = {
                .required_capabilities = ALEA_NUC_CAP_RESTRICTED_NEUTRON
            };
            alea_nuc_capability_report_t base_report;
            alea_nuc_prepared_material_t* prepared = NULL;
            alea_error_t base_status = alea_nuc_prepare_material(
                &material, &requirements, &base_report, &prepared);
            collision_ready = base_status == ALEA_OK;
            if (base_status == ALEA_OK) {
                capabilities |= base_report.available_capabilities;
                alea_nuc_prepared_material_free(prepared);
            } else {
                report = base_report;
                capabilities = base_report.available_capabilities;
            }
        }
        print_row(entry, "yes", "yes", collision_ready ? "yes" : "no",
                  capabilities, issue_name(report.issue), report.mt,
                  report.law, report.detail);
        alea_nuc_nuclide_free(nuc);
    }

    alea_nuc_xsdir_free(xsdir);
    return 0;
}
