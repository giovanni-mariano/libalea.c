// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
// SPDX-License-Identifier: MPL-2.0

import {readFile} from "node:fs/promises";

const moduleName = process.argv[2] ?? "alea.js";
const expectOpenmp = process.argv[3] === "openmp";
const largeInputMb = Number(process.argv[4] ?? 0);
const imported = await import(new URL(`../dist/${moduleName}`, import.meta.url));
const module = await imported.default({
  locateFile: (name) => new URL(`../dist/${name}`, import.meta.url).pathname,
});

const input = await readFile(new URL("../web/models/pin-cluster.mcnp", import.meta.url));

function errorMessage() {
  const pointer = module._alea_wasm_last_error();
  let end = pointer;
  while (module.HEAPU8[end]) end += 1;
  return new TextDecoder().decode(module.HEAPU8.subarray(pointer, end));
}

function loadMcnp(bytes) {
  const pointer = module._malloc(bytes.length);
  module.HEAPU8.set(bytes, pointer);
  const result = module._alea_wasm_load_mcnp(pointer, bytes.length);
  module._free(pointer);
  return result;
}

const openmp = module._alea_wasm_openmp_enabled() !== 0;
if (openmp !== expectOpenmp) {
  throw new Error(`unexpected OpenMP state: expected ${expectOpenmp}, got ${openmp}`);
}
let deck = input;
if (largeInputMb > 0) {
  const targetBytes = Math.floor(largeInputMb * 1024 * 1024);
  const padding = Buffer.alloc(Math.max(0, targetBytes - input.length), "c wasm heap growth test padding\n");
  deck = Buffer.concat([input, padding]);
}

if (module._alea_wasm_init(160, 90) !== 0 || loadMcnp(deck) !== 0) {
  throw new Error(`WASM MCNP load failed: ${errorMessage()}`);
}
if (module._alea_wasm_set_target(1.0, 2.0, 3.0) !== 0 ||
    module._alea_wasm_render(0.7, 0.35, 2.7, 0.5) !== 0) {
  throw new Error(`WASM render failed: ${errorMessage()}`);
}
const pointer = module._alea_wasm_pixels();
const length = module._alea_wasm_pixel_bytes();
const pixels = module.HEAPU8.subarray(pointer, pointer + length);
const colors = new Set();
for (let i = 0; i < pixels.length; i += 4) {
  colors.add(`${pixels[i]},${pixels[i + 1]},${pixels[i + 2]}`);
}
if (length !== 160 * 90 * 4 || colors.size < 4 ||
    module._alea_wasm_cell_count() !== 9 || module._alea_wasm_surface_count() !== 13) {
  throw new Error(`invalid MCNP render: ${length} bytes, ${colors.size} colors`);
}
const mode = openmp ? `${module._alea_wasm_openmp_max_threads()} OpenMP threads` : "single thread";
const inputDescription = largeInputMb > 0 ? `, ${largeInputMb} MiB input` : "";
console.log(`WASM MCNP smoke test: ${length} bytes, ${colors.size} colors, ${mode}${inputDescription}`);
module._alea_wasm_destroy();
