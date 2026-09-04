// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
// SPDX-License-Identifier: MPL-2.0

const WIDTH = 320;
const HEIGHT = 180;
const PREVIEW_WIDTH = 160;
const PREVIEW_HEIGHT = 90;

let module;
let threaded = false;

function readCString(pointer) {
  if (!pointer) return "unknown libalea error";
  let end = pointer;
  while (module.HEAPU8[end]) end += 1;
  return new TextDecoder().decode(module.HEAPU8.subarray(pointer, end));
}

function loadMcnp(input) {
  const bytes = typeof input === "string"
    ? new TextEncoder().encode(input)
    : new Uint8Array(input);
  const pointer = module._malloc(bytes.length);
  if (!pointer) {
    throw new Error(`not enough WASM memory for the ${(bytes.length / 1048576).toFixed(1)} MB MCNP input`);
  }
  module.HEAPU8.set(bytes, pointer);
  let result;
  try {
    result = module._alea_wasm_load_mcnp(pointer, bytes.length);
  } finally {
    module._free(pointer);
  }
  if (result !== 0) throw new Error(readCString(module._alea_wasm_last_error()));
  return {
    cells: module._alea_wasm_cell_count(),
    surfaces: module._alea_wasm_surface_count(),
    radius: module._alea_wasm_radius(),
    center: [
      module._alea_wasm_center_x(),
      module._alea_wasm_center_y(),
      module._alea_wasm_center_z(),
    ],
  };
}

async function createModule(wantThreads) {
  let moduleName = wantThreads ? "alea-openmp.js" : "alea.js";
  try {
    const imported = await import(new URL(`../dist/${moduleName}`, import.meta.url));
    module = await imported.default({
      locateFile: (name) => new URL(`../dist/${name}`, import.meta.url).href,
    });
    threaded = module._alea_wasm_parallel_enabled() !== 0;
  } catch (error) {
    if (!wantThreads) throw error;
    moduleName = "alea.js";
    const imported = await import(new URL("../dist/alea.js", import.meta.url));
    module = await imported.default({
      locateFile: (name) => new URL(`../dist/${name}`, import.meta.url).href,
    });
    threaded = false;
  }
  if (module._alea_wasm_init(WIDTH, HEIGHT) !== 0) {
    throw new Error(readCString(module._alea_wasm_last_error()));
  }
}

self.onmessage = async ({data}) => {
  try {
    if (data.type === "init") {
      await createModule(data.threaded);
      const response = await fetch("models/pin-cluster.mcnp");
      if (!response.ok) throw new Error(`could not load default MCNP model (${response.status})`);
      const model = loadMcnp(await response.arrayBuffer());
      postMessage({type: "ready", threaded, threads: module._alea_wasm_parallel_max_threads(), model});
      return;
    }
    if (data.type === "load" && module) {
      const model = loadMcnp(data.buffer);
      postMessage({type: "loaded", name: data.name, model});
      return;
    }
    if (data.type === "target" && module) {
      if (module._alea_wasm_set_target(data.x, data.y, data.z) !== 0) {
        throw new Error(readCString(module._alea_wasm_last_error()));
      }
      postMessage({type: "targeted", target: [data.x, data.y, data.z]});
      return;
    }
    if (data.type === "render" && module) {
      const startedAt = performance.now();
      const width = data.preview ? PREVIEW_WIDTH : WIDTH;
      const height = data.preview ? PREVIEW_HEIGHT : HEIGHT;
      if (module._alea_wasm_resize(width, height) !== 0) {
        throw new Error(readCString(module._alea_wasm_last_error()));
      }
      const clip = data.clip ? data.clipFraction : -1.0;
      if (module._alea_wasm_render(
        data.azimuth, data.elevation, data.distance, clip, data.preview ? 0 : 1) !== 0) {
        throw new Error(readCString(module._alea_wasm_last_error()));
      }
      const pointer = module._alea_wasm_pixels();
      const length = module._alea_wasm_pixel_bytes();
      const pixels = module.HEAPU8.slice(pointer, pointer + length);
      postMessage({
        type: "frame",
        pixels: pixels.buffer,
        renderMs: performance.now() - startedAt,
        width,
        height,
        preview: data.preview,
      }, [pixels.buffer]);
    }
  } catch (error) {
    postMessage({type: "error", message: error?.message ?? String(error)});
  }
};
