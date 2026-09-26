// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
// SPDX-License-Identifier: MPL-2.0

#include "alea_xml.h"
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>

int LLVMFuzzerTestOneInput(const uint8_t* data, size_t size) {
    if (!size || size > 1024 * 1024) return 0;
    alea_model_t* model = alea_xml_load_string((const char*)data, size);
    if (!model) return 0;
    const alea_system_t* sys = alea_model_system_const(model);
    (void)alea_cell_count(sys);
    (void)alea_surface_count(sys);
    FILE* stream = tmpfile();
    if (stream) {
        (void)alea_xml_export_stream(model, stream);
        fclose(stream);
    }
    alea_model_destroy(model);
    return 0;
}
