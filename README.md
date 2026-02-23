<!--
SPDX-FileCopyrightText: 2026 Giovanni MARIANO

SPDX-License-Identifier: MPL-2.0
-->

# libalea.c

A C library for debugging and analyzing large Constructive Solid Geometry (CSG) models used in neutron and gamma transport simulations. 

**The library is under active development. The API may change.**

## Hello World

```c
#include <alea.h>
#include <stdio.h>

int main(void) {
    alea_system_t* sys = alea_load_mcnp("geometry.inp");
    alea_build_universe_index(sys);

    int cell = alea_find_cell(sys, 100.0, 0.0, 0.0);
    int mat  = alea_material_at(sys, 100.0, 0.0, 0.0);
    printf("Cell %d, material %d\n", cell, mat);

    alea_destroy(sys);
}
```

```bash
gcc -o hello hello.c -Iinclude bin/libalea_full.a -lm
```

## What It Does

- **Load** MCNP and OpenMC geometry files
- **Query** which cell and material exists at any point
- **Detect** overlapping cells and undefined regions
- **Trace** rays through the model and report every cell crossing
- **Visualize** 2D cross-sections with exact analytical surface boundaries
- **Generate** void regions to fill gaps in the geometry
- **Convert** between MCNP and OpenMC formats
- **Build** geometry programmatically with boolean operations



## Installation

### Pre-built Binaries

Download pre-built binaries from [GitHub Releases](https://github.com/giovanni-mariano/libalea.c/releases):

| Platform | Archive |
|----------|---------|
| Linux x64 | `alea-linux-x64.tar.gz` |
| Linux ARM64 | `alea-linux-arm64.tar.gz` |
| macOS Intel | `alea-macos-x64.tar.gz` |
| macOS Apple Silicon | `alea-macos-arm64.tar.gz` |
| Windows x64 | `alea-windows-x64.zip` |

Each release includes the `alea` CLI, `mc_convert` and `mc_plotter` tools, static libraries, and headers.

### Building from Source

```bash
git clone https://github.com/giovanni-mariano/libalea.c.git
cd libalea.c
```

Build the library, CLI, and tools:

```bash
make lib          # Build all libraries
make cli          # Build the alea CLI tool
make tools        # Build mc_convert and mc_plotter
make test         # Build and run tests
```

Optional flags:

```bash
make USE_OPENMP=1 lib cli    # Enable OpenMP parallelization
make RELEASE=1 lib cli       # Optimized build
```

### Dependencies

- C11 compiler (gcc or clang)
- Standard math library (`-lm`)
- Optional: OpenMP for parallel slice rendering and ray tracing

### Libraries Produced

| Library | Contents |
|---------|----------|
| `libalea.a` | Core engine (CSG evaluation, export, dedup) |
| `libalea_mcnp.a` | MCNP parser and exporter |
| `libalea_openmc.a` | OpenMC XML parser and exporter |
| `libalea_raycast.a` | Ray tracing module |
| `libalea_slice.a` | 2D slice and visualization module |
| `libalea_full.a` | Everything in one archive |

Link against `libalea_full.a` unless you need to minimize binary size.

## CLI Tool

The `alea` CLI provides an interactive Lua environment for geometry analysis:

```bash
bin/alea                           # Interactive REPL
bin/alea script.lua                # Run a Lua script
bin/alea script.lua arg1 arg2      # Pass arguments to script
```

Example session:

```lua
> sys = alea.load_mcnp("model.inp")
> sys:build_universe_index()
> sys:find_cell(0, 0, 0)
1
> sys:material_at(0, 0, 0)
1
```

See `examples/lua/` for complete Lua scripts demonstrating all features.

## Tools

Pre-built tools ship with every release and are also built via `make tools`:

| Tool | Description |
|------|-------------|
| `mc_convert` | Convert between MCNP and OpenMC geometry formats |
| `mc_plotter` | Render 2D cross-section slices of CSG geometry to PNG/BMP |

```bash
bin/mc_convert model.inp model.xml          # MCNP -> OpenMC
bin/mc_plotter model.inp Z 0 -100 100 -100 100 800 output.png
```

## Examples

The `examples/` directory contains complete working programs:

| Example | What it shows |
|---------|---------------|
| `basic.c` | Build geometry from scratch, point queries, void generation, MCNP export |
| `mcnp_roundtrip.c` | Parse and re-export an MCNP input file |
| `mcnp_volume.c` | Estimate cell volumes via Monte Carlo ray tracing |
| `render3d.c` | Render a 3D image of the geometry |

```bash
cd examples/c && make
./basic
```

## Documentation

| Document | Audience | Purpose |
|----------|----------|---------|
| [Tutorial](docs/TUTORIAL.md) | New users | C API walk-through from loading a model to exporting results |
| [Lua Tutorial](docs/LUA_TUTORIAL.md) | New users | Lua API for scripting and interactive use |
| [Concepts](docs/CONCEPTS.md) | All users | Surfaces, sense, cells, universes, lattices, and other domain concepts |
| [Architecture](docs/ARCHITECTURE.md) | Contributors | Internal data model, algorithms, and design decisions |
| [API Reference](docs/API.md) | All users | Every public function, grouped by task |

Start with the **Tutorial** (C) or **Lua Tutorial** if you're new. Refer to **Concepts** when something doesn't behave as you expect. The **API Reference** is for when you know what you want but forgot the function name.

## Project Structure

```
include/               Public headers
  alea.h               Main API
  alea_types.h         Type definitions
  alea_raycast.h       Ray tracing
  alea_slice.h         2D visualization
src/
  core/                CSG engine, evaluation, export, dedup
  primitives/          Geometric primitives (plane, sphere, cylinder, ...)
  mcnp/                MCNP parser, converter, exporter
  openmc/              OpenMC XML parser, converter, exporter
  raycast/             Ray-geometry intersection, BVH
  slice/               2D slice curves, grid queries
  render/              Image output (PPM, BMP)
  lua_bind/            Lua bindings for CLI
  util/                Arena allocator, logging, vectors, math
tools/               mc_convert, mc_plotter
examples/
  lua/                 Lua example scripts
tests/
  unit/                Unit tests
  integration/         Integration tests
  lua/                 Lua tests
  fuzz/                Fuzz testing (MCNP and OpenMC parsers)
  data/                Test geometry files (MCNP, OpenMC XML)
vendor/                Third-party dependencies (Lua, linenoise)
```

## Disclaimer

This package was developed with support of AI tools.

## License

Mozilla Public License 2.0 (MPL-2.0). See [LICENSES/MPL-2.0.txt](LICENSE) for details.
