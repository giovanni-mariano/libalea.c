<!--
SPDX-FileCopyrightText: 2026 Giovanni MARIANO

SPDX-License-Identifier: MPL-2.0
-->

# libalea.c

A C library for building, debugging, and analyzing Constructive Solid Geometry (CSG) models used in neutron and gamma transport simulations.

**The library is under active development. The API may change.**

## Hello World

```c
#include <alea.h>
#include <alea_mcnp.h>
#include <stdio.h>

int main(void) {
    mcnp_model_t* model = mcnp_load("geometry.inp");
    if (!model) {
        fprintf(stderr, "load failed: %s\n", alea_error());
        return 1;
    }

    alea_system_t* sys = model->sys;
    alea_build_universe_index(sys);

    int cell = 0;
    int mat = 0;
    if (alea_find_cell_at(sys, 100.0, 0.0, 0.0, &cell, &mat) == 0)
        printf("Cell %d, material %d\n", cell, mat);

    mcnp_model_destroy(model);
    return 0;
}
```

```bash
gcc -o hello hello.c -Iinclude bin/libalea_mcnp.a bin/libalea.a -lm
```

## What It Does

- **Load** MCNP and OpenMC geometry files
- **Query** which cell and material exists at any point
- **Detect** overlapping cells and undefined regions
- **Trace** rays through the model and report every cell crossing
- **Visualize** 2D cross-sections with exact analytical surface boundaries
- **Render** 3D images with Phong shading, cutaway views, and shadow rays
- **Export** structured hex meshes to Gmsh (.msh) and VTK (.vtk) formats (exp.)
- **Generate** void regions to fill gaps in the geometry
- **Convert** between MCNP and OpenMC formats
- **Build** geometry programmatically with boolean operations
- **Materials** definition with nuclide/element composition and mixture support
- **Nuclear data** read ACE-format cross sections (neutron, photon), sample free paths, target nuclides, reaction MTs, and multigroup scattering, Doppler broaden, and collapse to multigroup constants

## Installation

### Pre-built Binaries

Download pre-built binaries from [GitHub Releases](https://github.com/giovanni-mariano/libalea/releases):

| Platform | Archive |
|----------|---------|
| Linux x64 | `alea-linux-x64.tar.gz` |
| Linux ARM64 | `alea-linux-arm64.tar.gz` |
| macOS Intel | `alea-macos-x64.tar.gz` |
| macOS Apple Silicon | `alea-macos-arm64.tar.gz` |
| Windows x64 | `alea-windows-x64.zip` |

The release workflow packages the `alea` CLI, `mc_convert`, `mc_plotter`, static libraries, and headers. The `nuc_plot` tool is built from source with `make tools`.

### Building from Source

```bash
git clone --recursive https://github.com/giovanni-mariano/libalea.git
cd libalea
```

If you already cloned without submodules, initialize the vendored Lua and linenoise sources before building the CLI:

```bash
git submodule update --init --recursive
```

Build the library, CLI, and tools:

```bash
make              # Build core library (bin/libalea.a)
make modules      # Build format modules (libalea_mcnp.a, libalea_openmc.a, libalea_nucdata.a)
make full         # Build everything into libalea_full.a
make cli          # Build the alea CLI tool
make tools        # Build mc_convert, mc_plotter, and nuc_plot
make test         # Build and run tests
make test-lua     # Build the CLI and run Lua tests
```

Optional flags:

```bash
make USE_OPENMP=1 full cli    # Enable OpenMP parallelization
make RELEASE=1 full cli       # Optimized build
```

### Dependencies

- C11 compiler (gcc or clang)
- Standard math library (`-lm`)
- Optional: OpenMP for parallel rendering and ray tracing

### Libraries Produced

| Library | Contents |
|---------|----------|
| `libalea.a` | Core engine: CSG evaluation, primitives, raycast, slice, 3D render, mesh export |
| `libalea_mcnp.a` | MCNP parser, converter, and exporter |
| `libalea_openmc.a` | OpenMC XML parser, converter, and exporter |
| `libalea_nucdata.a` | Nuclear data: ACE reader, cross-section lookup, free-path/nuclide/reaction sampling, Doppler broadening, multigroup collapse |
| `libalea_full.a` | Everything in one archive |

Link against the core library plus the format modules you need:

```bash
# MCNP support
gcc -o myapp myapp.c -Iinclude bin/libalea_mcnp.a bin/libalea.a -lm

# Nuclear data
gcc -o myapp myapp.c -Iinclude bin/libalea_nucdata.a bin/libalea.a -lm

# Full library (core + all formats + nucdata)
gcc -o myapp myapp.c -Iinclude bin/libalea_full.a -lm
```

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

Tools are built via `make tools`:

| Tool | Description |
|------|-------------|
| `mc_convert` | Convert between MCNP and OpenMC geometry formats |
| `mc_plotter` | Render 2D cross-section slices of CSG geometry to PNG/BMP |
| `nuc_plot` | Generate SVG plots of nuclear cross sections, angular distributions, fission spectra, and more |

```bash
bin/mc_convert model.inp model.xml
bin/mc_plotter model.inp Z 0 -100 100 -100 100 800x800 output.png
bin/nuc_plot --xsdir /path/to/xsdir --zaid 92235.80c --plot xs --output u235.svg
```

## Sampling

The tracked public sampling API is in `include/alea_nucdata.h` and `include/alea_mesh.h`:

| API | What it samples |
|-----|-----------------|
| `alea_nuc_sample_distance` | Distance to the next collision from macroscopic total cross section |
| `alea_nuc_sample_nuclide` | Target nuclide in a material |
| `alea_nuc_sample_reaction` | Reaction MT on a selected nuclide |
| `alea_nuc_urr_factors` | Unresolved-resonance probability-table factors |
| `alea_nuc_mg_sample_scatter` | Outgoing multigroup scatter group |
| `alea_mesh_sample` | Structured mesh cells from CSG geometry |

Full outgoing collision kinematics are not exposed by the tracked public headers or Makefile.

## Examples

### C Examples

The `examples/c/` directory contains complete working programs:

| Example | What it shows |
|---------|---------------|
| `basic.c` | Build geometry from scratch, point queries, void generation, MCNP export |
| `mcnp_roundtrip.c` | Parse and re-export an MCNP input file |
| `mcnp_volume.c` | Estimate cell volumes via Monte Carlo ray tracing |
| `mcnp_mesh.c` | Export MCNP geometry as a structured hex mesh (Gmsh/VTK) |
| `render3d.c` | Render a 3D image of the geometry |
| `void_demo.c` | Void generation with explicit bounds and multiple cells |

```bash
cd examples/c && make
./basic
```

The examples Makefile currently builds `basic`, `mcnp_roundtrip`, `mcnp_volume`, `render3d`, and `void_demo`. `mcnp_mesh.c` can be compiled directly against the tracked libraries.

### Lua Examples

The `examples/lua/` directory contains scripts for use with the `alea` CLI:

| Script | What it shows |
|--------|---------------|
| `01_hello.lua` | Load a model and query a point |
| `02_build_geometry.lua` | Programmatic geometry construction |
| `03_point_queries.lua` | Cell and material queries |
| `04_model_inspection.lua` | Inspect cells, surfaces, universes |
| `05_format_conversion.lua` | MCNP/OpenMC conversion |
| `06_volume_estimation.lua` | Monte Carlo volume estimation |
| `07_overlap_check.lua` | Detect overlapping cells |
| `08_flatten_and_simplify.lua` | Flatten universes and simplify CSG |
| `09_universe_extract_merge.lua` | Extract and merge universes |
| `10_parametric_geometry.lua` | Parametric model generation |
| `11_raycast.lua` | Ray tracing through geometry |
| `12_slice.lua` | 2D cross-section slicing |
| `13_render.lua` | 3D rendering |
| `14_mesh.lua` | Mesh export |
| `15_analysis_pipeline.lua` | Full analysis workflow |
| `16_cell_comments.lua` | Cell comment handling |
| `17_build_with_comments.lua` | Building geometry with comments |
| `18_materials_and_mixtures.lua` | Material and mixture definitions |
| `19_nucdata.lua` | Nuclear data: load ACE cross sections, query, build materials |
| `20_nucdata_plots.lua` | Nuclear data SVG plotting |

## Documentation

| Document | Audience | Purpose |
|----------|----------|---------|
| [Tutorial](docs/TUTORIAL.md) | New users | C API walk-through from loading a model to exporting results |
| [Lua Tutorial](docs/LUA_TUTORIAL.md) | New users | Lua API for scripting and interactive use |
| [Concepts](docs/CONCEPTS.md) | All users | Surfaces, sense, cells, universes, lattices, and other domain concepts |
| [Architecture](docs/ARCHITECTURE.md) | Contributors | Internal data model, algorithms, and design decisions |
| [API Reference](docs/API.md) | All users | Every public function, grouped by task |
| [Mesh Sampling Plan](docs/PLAN_MESH_SAMPLING.md) | Contributors | Plan for hardening mesh sampling and adding mixed-material support |

Start with the **Tutorial** (C) or **Lua Tutorial** if you're new. Refer to **Concepts** when something doesn't behave as you expect. The **API Reference** is for when you know what you want but forgot the function name.

## Project Structure

```
include/               Public headers
  alea.h               Main API (CSG engine, queries, void, materials)
  alea_types.h         Type definitions
  alea_raycast.h       Ray tracing API
  alea_slice.h         2D slice/visualization API
  alea_render.h        3D rendering API
  alea_mesh.h          Mesh export API
  alea_nucdata.h       Nuclear data API
  alea_nucdata_types.h Nuclear data type definitions
  alea_mcnp.h          MCNP module API
  alea_openmc.h        OpenMC module API
src/
  core/                CSG engine, evaluation, export, dedup, void, materials
  primitives/          Geometric primitives (plane, sphere, cylinder, cone, torus, ...)
  mcnp/                MCNP parser, converter, exporter
    parser/            Lexer and parser
    conversion/        Surface and cell conversion
    exporter/          MCNP output formatting
  nucdata/             Nuclear data: ACE reader, XS lookup, public sampling APIs, Doppler, multigroup
  openmc/              OpenMC XML parser, converter, exporter
  raycast/             Ray-geometry intersection, BVH
  slice/               2D slice curves, analytical intersection
  render/              3D batch renderer (Phong, shadows, cutaway)
  mesh/                Structured hex mesh export (Gmsh, VTK)
  lua_bind/            Lua bindings for CLI
  util/                Arena allocator, logging, vectors, math
tools/               mc_convert, mc_plotter, nuc_plot
examples/
  c/                   C example programs
  lua/                 Lua example scripts
tests/
  unit/                Unit tests
  integration/         Integration tests
  nucdata/             Nuclear data tests
  lua/                 Lua tests
  fuzz/                Fuzz testing and corpora for MCNP and OpenMC parsers
  data/                Test geometry files (MCNP, OpenMC XML)
vendor/                Third-party dependencies (Lua, linenoise)
```

## Disclaimer

This package was developed with support of AI tools.

## License

Mozilla Public License 2.0 (MPL-2.0). See [LICENSE](LICENSES/MPL-2.0.txt) for details.
