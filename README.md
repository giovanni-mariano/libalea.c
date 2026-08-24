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
- **Convert** between MCNP, OpenMC, and Serpent geometry formats
- **Build** geometry programmatically with boolean operations
- **Materials** definition with nuclide/element composition and mixture support
- **Nuclear data** read ACE-format cross sections (neutron, photon), sample free paths, target nuclides, reaction MTs, and multigroup scattering, Doppler broaden, and collapse to multigroup constants

## Installation

### Pre-built Binaries

Download pre-built binaries from [GitHub Releases](https://github.com/giovanni-mariano/libalea/releases):

| Platform | Archive |
|----------|---------|
| Linux x64 | `alea-linux-x64.tar.gz` |
| Linux x64 (conda + OpenMP) | `alea-linux-conda-openmp-x64.tar.gz` |
| Linux ARM64 | `alea-linux-arm64.tar.gz` |
| macOS Intel | `alea-macos-x64.tar.gz` |
| macOS Apple Silicon | `alea-macos-arm64.tar.gz` |
| Windows x64 (MinGW/UCRT) | `alea-windows-x64.zip` |
| Windows x64 (MSVC) | `alea-windows-msvc-x64.zip` |
| Windows x64 (MSVC + OpenMP) | `alea-windows-msvc-openmp-x64.zip` |

The Linux, macOS, and MinGW/UCRT Windows archives package the `alea` CLI,
`mc_convert`, `mc_plotter`, `nuc_plot`, `large_model_probe`, static libraries,
and headers. The MSVC archives package the `.lib` static libraries and headers.
Applications linked against the conda OpenMP Linux archive or the OpenMP-enabled
MSVC archive require the LLVM OpenMP runtime at run time:

```powershell
conda install -c conda-forge llvm-openmp
```

This runtime requirement does not require Visual Studio, Windows SDK headers, or
compiler tools; those are needed only when compiling or linking applications.
On Linux, run from an activated conda environment or make sure the conda
environment's `lib` directory is on the runtime library search path.

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
make modules      # Build format modules (libalea_mcnp.a, libalea_openmc.a, libalea_serpent.a, libalea_nucdata.a)
make full         # Build everything into libalea_full.a
make cli          # Build the alea CLI tool
make tools        # Build mc_convert, mc_plotter, nuc_plot, and large_model_probe
make test         # Build and run tests
make test-lua     # Build the CLI and run Lua tests
make install      # Install libraries, headers, CLI, tools, and docs
```

### Build Options by Platform

The default build does not enable OpenMP. Add `USE_OPENMP=1` to compile the
same targets with OpenMP parallelization enabled. Add `RELEASE=1` for an
optimized build.

Common Makefile variables:

```bash
make CC=clang full cli tools                  # Select compiler
make PREFIX=/opt/libalea install             # Install prefix
make DESTDIR=/tmp/pkg PREFIX=/usr install    # Package/stage install
make LIBOMP_PREFIX=/path/to/libomp USE_OPENMP=1 full cli tools
```

`make install` builds and installs the static libraries, public headers, `alea`
CLI, tools, README, and license files. Use `install-libs`, `install-cli`, or
`install-tools` to install only one part.

#### Linux

Non-OpenMP build:

```bash
make full cli tools
make test
```

OpenMP build with gcc:

```bash
make clean
make USE_OPENMP=1 RELEASE=1 full cli tools
make USE_OPENMP=1 test
```

OpenMP build with clang requires an OpenMP runtime such as `libomp`:

```bash
make clean
make CC=clang USE_OPENMP=1 RELEASE=1 full cli tools
make CC=clang USE_OPENMP=1 test
```

If `libomp` is installed in a non-standard location, pass its prefix:

```bash
make CC=clang USE_OPENMP=1 LIBOMP_PREFIX=/opt/libomp full cli tools
```

#### macOS

Non-OpenMP build:

```bash
make full cli tools
make test
```

OpenMP builds with Apple Clang require Homebrew `libomp`:

```bash
brew install libomp
make clean
make USE_OPENMP=1 RELEASE=1 full cli tools
make USE_OPENMP=1 test
```

The Makefile detects `libomp` with `brew --prefix libomp` and falls back to
`/usr/local/opt/libomp` if Homebrew is not on `PATH`.

Override the detected path when needed:

```bash
make USE_OPENMP=1 LIBOMP_PREFIX=/opt/homebrew/opt/libomp full cli tools
```

#### Windows with MinGW/UCRT

Use the MSYS2 UCRT64 environment with the MinGW-w64 GCC toolchain. Install the
build tools from an MSYS2 shell:

```bash
pacman -S --needed mingw-w64-ucrt-x86_64-gcc mingw-w64-ucrt-x86_64-make make
```

Open a UCRT64 shell, clone the repository with submodules, and build with the
GNU Makefile:

```bash
make full cli tools
make test-unit test-integration
```

OpenMP builds require the MinGW OpenMP runtime:

```bash
pacman -S --needed mingw-w64-ucrt-x86_64-omp
make clean
make USE_OPENMP=1 RELEASE=1 full cli tools
make USE_OPENMP=1 test-unit test-integration
```

#### Windows with conda clang-cl

Use this path when users can install conda packages without admin rights and the
machine already has Windows SDK/MSVC headers and import libraries available.
The compiler, archive tool, OpenMP runtime, and `jom` are supplied by conda, but
`clang-cl` still targets the MSVC ABI.

```powershell
conda create -n libalea-clang -c conda-forge clang_win-64 llvm-openmp jom
conda activate libalea-clang

jom /J 1 /f Makefile.msvc CONDA_CLANG=1 USE_OPENMP=1 full
jom /J 1 /f Makefile.msvc CONDA_CLANG=1 USE_OPENMP=1 test
```

If `clang-cl` reports missing headers such as `vcruntime.h` or `windows.h`, the
Windows SDK/MSVC C++ components are not visible to the shell. Open a VS x64
developer prompt before activating conda, or install the required SDK/toolset
components. The `/J 1` option keeps `jom` serial because its dependency handling
differs from `nmake` for this makefile.

#### Windows with conda MinGW/UCRT

For machines without Windows SDK/MSVC headers, use the conda-forge MinGW-w64
UCRT toolchain. It provides its own compiler, headers, and runtime inside the
conda environment.

```powershell
conda create -n libalea-ucrt `
  -c conda-forge/label/m2w64-experimental -c conda-forge `
  mingw-w64-ucrt-x86_64-gcc mingw-w64-ucrt-x86_64-binutils make
conda activate libalea-ucrt

$cc = (Get-Command x86_64-w64-mingw32-gcc).Source
$ar = & $cc -print-prog-name=ar
make WINDOWS_GNU=1 CC="$cc" AR="$ar" USE_OPENMP=1 full cli tools
make WINDOWS_GNU=1 CC="$cc" AR="$ar" USE_OPENMP=1 test-unit test-integration test-lua
```

#### Windows with MSVC

Use Visual Studio 2019 or newer with the C++ x64 toolset installed. The wrapper
scripts locate Visual Studio, enter the x64 developer environment, and run
`nmake /f Makefile.msvc`. The MSVC build currently covers the static libraries
and tests; the Lua CLI target is built by the GNU Makefile.

Non-OpenMP build from PowerShell:

```powershell
.\build-msvc.ps1 full
.\build-msvc.ps1 test
```

OpenMP build from PowerShell:

```powershell
.\build-msvc.ps1 USE_OPENMP=1 RELEASE=1 full
.\build-msvc.ps1 USE_OPENMP=1 test
```

The same commands are available from `cmd.exe`:

```bat
build-msvc.bat full
build-msvc.bat USE_OPENMP=1 RELEASE=1 full
```

MSVC OpenMP uses `/openmp:llvm`, so install the Visual Studio LLVM/OpenMP
runtime component if OpenMP executables cannot find `libomp140.x86_64.dll`.

### Dependencies

- C11 compiler (gcc, clang, or MSVC)
- `make` on Linux/macOS/MinGW, or Visual Studio `nmake` on Windows
- Standard math library (`-lm`)
- Optional: OpenMP for parallel rendering and ray tracing

### Libraries Produced

| Library | Contents |
|---------|----------|
| `libalea.a` | Core engine: CSG evaluation, primitives, raycast, slice, 3D render, mesh export |
| `libalea_mcnp.a` | MCNP parser, converter, and exporter |
| `libalea_openmc.a` | OpenMC XML parser, converter, and exporter |
| `libalea_serpent.a` | Serpent exporter |
| `libalea_nucdata.a` | Nuclear data: ACE reader, cross-section lookup, free-path/nuclide/reaction sampling, Doppler broadening, multigroup collapse |
| `libalea_full.a` | Core, MCNP, OpenMC, Serpent, and nuclear data in one archive |

Link against the core library plus the format modules you need:

```bash
# MCNP support
gcc -o myapp myapp.c -Iinclude bin/libalea_mcnp.a bin/libalea.a -lm

# OpenMC support
gcc -o myapp myapp.c -Iinclude bin/libalea_openmc.a bin/libalea.a -lm

# Nuclear data
gcc -o myapp myapp.c -Iinclude bin/libalea_nucdata.a bin/libalea.a -lm

# Serpent export
gcc -o myapp myapp.c -Iinclude bin/libalea_serpent.a bin/libalea.a -lm

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
| `mc_convert` | Convert between MCNP, OpenMC, and Serpent geometry formats |
| `mc_plotter` | Render 2D cross-section slices of CSG geometry to PNG/BMP |
| `nuc_plot` | Generate SVG plots of nuclear cross sections, angular distributions, fission spectra, and more |
| `large_model_probe` | Inspect large MCNP models and benchmark hierarchy query/raycast behavior |

```bash
bin/mc_convert model.inp model.xml
bin/mc_convert model.inp model.serp --output-format serpent
bin/mc_plotter model.inp Z 0 -100 100 -100 100 800x800 output.png
bin/nuc_plot --xsdir /path/to/xsdir --zaid 92235.80c --plot xs --output u235.svg
bin/large_model_probe model.inp --queries 10000 --hier-build
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

### First-visible and boundary-event queries

Use `sys:first_visible()` when only the frontmost material interval is needed;
it avoids building a full Lua raycast result. `sys:boundary_events()` returns
ordered ownership transitions for diagnostics and surface provenance.

```lua
local hit = sys:first_visible(-10, 0, 0, 1, 0, 0,
    {t_max = 100, normal = true})
if hit then print(hit.cell_id, hit.t, hit.surface_id) end

for _, event in ipairs(sys:boundary_events(-10, 0, 0, 1, 0, 0,
                                           {t_max = 100})) do
    print(event.t, event.surface_id, event.cell_before, event.cell_after)
end
```

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
  alea_serpent.h       Serpent exporter API
src/
  core/                CSG engine, evaluation, export, dedup, void, materials
  primitives/          Geometric primitives (plane, sphere, cylinder, cone, torus, ...)
  mcnp/                MCNP parser, converter, exporter
    parser/            Lexer and parser
    conversion/        Surface and cell conversion
    exporter/          MCNP output formatting
  nucdata/             Nuclear data: ACE reader, XS lookup, public sampling APIs, Doppler, multigroup
  openmc/              OpenMC XML parser, converter, exporter
  serpent/             Serpent exporter
  raycast/             Ray-geometry intersection, BVH
  slice/               2D slice curves, analytical intersection
  render/              3D batch renderer (Phong, shadows, cutaway)
  mesh/                Structured hex mesh export (Gmsh, VTK)
  lua_bind/            Lua bindings for CLI
  util/                Arena allocator, logging, vectors, math
tools/               mc_convert, mc_plotter, nuc_plot
                     large_model_probe
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
