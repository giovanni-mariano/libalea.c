// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
// SPDX-License-Identifier: MPL-2.0

#include "alea_wasm.h"

#include <stdint.h>
#include <stdio.h>
#include <string.h>

static const char input[] =
    "WASM MCNP test sphere\n"
    "1 1 -1.0 -1\n"
    "99 0 1\n\n"
    "1 so 5\n\n"
    "m1 1001.80c 1.0\n";

int main(void) {
    if (!alea_wasm_version() || alea_wasm_version()[0] == '\0') {
        fprintf(stderr, "binding returned no version\n");
        return 1;
    }
    if (alea_wasm_init(160, 90) != 0) {
        fprintf(stderr, "init failed: %s\n", alea_wasm_last_error());
        return 1;
    }
    if (alea_wasm_pick(0, 0) >= 0 ||
        alea_wasm_set_picking_enabled(1) != 0 ||
        alea_wasm_load_mcnp(input, (int)strlen(input)) != 0) {
        fprintf(stderr, "load failed: %s\n", alea_wasm_last_error());
        return 1;
    }
    if (alea_wasm_set_target(1.0, 2.0, 3.0) != 0 ||
        alea_wasm_resize(80, 45) != 0 ||
        alea_wasm_render(0.7, 0.35, 2.7, 0.5, 0) != 0 ||
        alea_wasm_resize(160, 90) != 0 ||
        alea_wasm_render(0.7, 0.35, 2.7, 0.5, 1) != 0) {
        fprintf(stderr, "render failed: %s\n", alea_wasm_last_error());
        alea_wasm_destroy();
        return 1;
    }
    const uint8_t* pixels = alea_wasm_pixels();
    const int bytes = alea_wasm_pixel_bytes();
    uint64_t checksum = 1469598103934665603ULL;
    for (int i = 0; i < bytes; ++i) {
        checksum ^= pixels[i];
        checksum *= 1099511628211ULL;
    }
    if (!pixels || bytes != 160 * 90 * 4 || alea_wasm_cell_count() != 2 ||
        alea_wasm_surface_count() != 1 || alea_wasm_radius() < 4.9) {
        fprintf(stderr, "binding returned invalid model data\n");
        alea_wasm_destroy();
        return 1;
    }
    int picked = 0;
    for (int y = 0; y < alea_wasm_height() && !picked; ++y) {
        for (int x = 0; x < alea_wasm_width(); ++x) {
            int hit = alea_wasm_pick(x, y);
            if (hit < 0) {
                fprintf(stderr, "pick failed: %s\n", alea_wasm_last_error());
                alea_wasm_destroy();
                return 1;
            }
            if (hit) {
                picked = alea_wasm_pick_cell_id() > 0 &&
                         alea_wasm_pick_material_id() > 0 &&
                         alea_wasm_pick_depth() > 0.0;
                break;
            }
        }
    }
    if (!picked) {
        fprintf(stderr, "rendered model did not expose a pickable pixel\n");
        alea_wasm_destroy();
        return 1;
    }
    printf("wasm MCNP binding: %d bytes, checksum %llu\n", bytes,
           (unsigned long long)checksum);
    alea_wasm_destroy();
    return 0;
}
