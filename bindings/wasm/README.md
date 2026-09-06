<!--
SPDX-FileCopyrightText: 2026 Giovanni MARIANO
SPDX-License-Identifier: MPL-2.0
-->

# libalea MCNP WebAssembly binding

This directory provides a browser-facing MCNP geometry loader and 3D renderer.
MCNP input text is parsed by `libalea_mcnp` inside WebAssembly, accelerated by
libalea, and rendered from a perspective orbit camera. The included viewer has
an MCNP file picker, coordinate targeting, interactive orbit/zoom, a cutaway
plane, and WebM capture.

For responsive interaction, camera motion uses a 160x90 preview without edge
post-processing and is capped at 15 frames per second. Pausing the orbit or
leaving an interaction idle for 180 ms produces an edge-enhanced 320x180
frame. Recording keeps the full 320x180 resolution. Both WASM variants are
built with SIMD enabled and use the renderer's native 32x32 tile size.

`web/models/pin-cluster.mcnp` is a small default model, not geometry hard-coded
into the binding. Loading another `.mcnp`, `.inp`, `.i`, or text input replaces
it at runtime without rebuilding the module.

## Prerequisites

- Emscripten 6.0.3 or newer for the threaded build (tested with 6.0.9)
- GNU Make
- Python 3 for the local development server

Activate Emscripten, then build and serve the portable module:

```sh
source /path/to/emsdk/emsdk_env.sh
make single
python3 serve.py
```

Open <http://localhost:8000/web/>.

For browser threads, also build the TinyPar variant:

```sh
make threaded PTHREAD_POOL_SIZE=8
python3 serve.py
```

The viewer uses `alea-threaded.js` when the page is cross-origin isolated and
falls back to `alea.js`. `serve.py` sends the COOP/COEP headers required for
`SharedArrayBuffer`; production hosting must do the same for all assets.

Both builds use a growable WASM heap with a 2 GiB ceiling. The threaded build
starts at 256 MiB so large MCNP inputs can be parsed without an immediate heap
resize; override `WASM_INITIAL_MEMORY` or `WASM_MAXIMUM_MEMORY` at build time
when a deployment needs different limits. Uploaded files are transferred to
the render worker as raw bytes to avoid an additional UTF-16 text copy.

## Tests

```sh
make test
make test-threaded
make test-large-threaded LARGE_INPUT_MB=64
```

The tests parse MCNP text, prepare geometry acceleration, render a cutaway 3D
frame, and validate the RGBA output. The Node test uses the same reactor-pin
input as the browser.

## JavaScript-facing API

The facade owns one parsed `mcnp_model_t` and framebuffer. Its main operations
are `alea_wasm_init`, `alea_wasm_load_mcnp`, and `alea_wasm_render`. It also
exports pixels, dimensions, model counts/bounds, parallel capabilities, error
reporting, and destruction. JavaScript copies UTF-8 MCNP text into WASM memory
with the exported `malloc`/`free` functions before calling the loader.
