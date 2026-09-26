<!--
SPDX-FileCopyrightText: 2026 Giovanni MARIANO

SPDX-License-Identifier: MPL-2.0
-->

# libalea.c

A C library for building, debugging, and analyzing Constructive Solid Geometry (CSG) models used in neutron and gamma transport simulations.

**The library is under active development. The API may change.**

**[Read the documentation](https://giovanni-mariano.github.io/libalea.c/)**

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
- **Detect** overlapping cells and undefined regions with sampled diagnostics
  or bounded, verified slice scans
- **Trace** rays through the model and report every cell crossing
- **Visualize** 2D cross-sections with exact analytical surface boundaries
- **Render** 3D images with Phong shading, cutaway views, and shadow rays
- **Sample and export** CSG onto structured rectilinear grids in Gmsh (.msh) and VTK (.vtk) formats (exp.)
- **Generate** void regions to fill gaps in the geometry
- **Convert** between MCNP, OpenMC, and Serpent geometry formats
- **Build** geometry programmatically with boolean operations
- **Materials** definition with nuclide/element composition and mixture support
- **Nuclear data** read ACE-format cross sections (neutron, photon), sample free paths, target nuclides, reaction MTs, and multigroup scattering, Doppler broaden, and collapse to multigroup constants

## Installation

### Pre-built Binaries

Download pre-built binaries from [GitHub Releases](https://github.com/giovanni-mariano/libalea.c/releases):

| Platform | Archive |
|----------|---------|
| Linux x64 | `alea-linux-x64.tar.gz` |
| Linux ARM64 | `alea-linux-arm64.tar.gz` |
| macOS Intel | `alea-macos-x64.tar.gz` |
| macOS Apple Silicon | `alea-macos-arm64.tar.gz` |
| Windows x64 (MinGW/UCRT) | `alea-windows-x64.zip` |
| Windows x64 (MSVC) | `alea-windows-msvc-x64.zip` |

The Linux, macOS, and MinGW/UCRT Windows archives package the `alea` CLI,
`mc_convert`, `mc_plotter`, `nuc_plot`, `nuc_inventory`, `large_model_probe`, static libraries,
and headers. The MSVC archives package the `.lib` static libraries and headers.
All builds use the vendored TinyPar backend by default and have no OpenMP
runtime dependency.

Release pages also contain the native Python bindings as `pyalea` archives for
CPython 3.10–3.14 on Linux x86_64/aarch64, macOS x86_64/arm64, and Windows
x86_64. Each archive is specific to its Python ABI and platform and includes a
matching `.sha256` file. Applications can bundle the contained `pyalea/`
directory and use `import pyalea`; see
[`bindings/python/README.md`](bindings/python/README.md) for the artifact layout
and local build commands.

### Building from Source

```bash
git clone --recursive https://github.com/giovanni-mariano/libalea.c.git
cd libalea.c
```

If you already cloned without submodules, initialize the vendored Lua and linenoise sources before building the CLI:

```bash
git submodule update --init --recursive
```

Build the library, CLI, and tools:

```bash
make              # Build core library (bin/libalea.a)
make modules      # Build optional format, nuclear-data, and transport modules
make transport    # Build transport with its geometry and nuclear-data dependencies
make full         # Build everything into libalea_full.a
make cli          # Build the alea CLI tool
make lua-module   # Build the Lua 5.5 require("alea") module
make tools        # Build command-line conversion, plotting, and inspection tools
make test         # Build and run tests
make test-lua     # Build the CLI and run Lua tests
make test-lua-module   # Test the module in a plain Lua host
make test-lua-cluster  # Test the optional one-process cluster binding
make install      # Install libraries, headers, CLI, tools, and docs
make -C bindings/python test PYTHON=python3  # Build and test pyAlea
```

### Build Options by Platform

The default build uses the TinyPar `native` backend (POSIX threads on
Linux/macOS, Win32 threads on Windows). Set `TINYPAR_BACKEND=serial` for an
explicitly serial build. Add `RELEASE=1` for an optimized build.

Threaded builds reuse a process-wide TinyPar executor across parallel library
operations. Set `ALEA_NUM_THREADS=n`, or call
`alea_parallel_set_threads(n)` before the first parallel operation, to select
its worker count. An explicit API call overrides the environment; zero restores
the hardware default. Calls submitted concurrently share the bounded executor
rather than creating independent native thread teams.
Nested calls and operations that need only one worker take a direct serial
path without creating or waking the process executor.
Per-operation worker arguments limit active participants and scratch usage;
they do not shrink the persistent native team. Set the process limit before
first use when the number of created native threads must also be bounded.

On POSIX, a child created with `fork()` abandons the inherited worker state and
lazily creates its own executor on the first parallel operation.
Inherited heap storage is retained by the fork handler and freed before that
replacement is created, or during normal process-exit cleanup. Copied locks
and condition variables are never destroyed or used during this reclamation.

Common Makefile variables:

For GNU Make release builds, `PORTABLE=1` keeps the compiler's default CPU
target; `PORTABLE=0` (or leaving it unset) enables `-march=native`. Both use
`-O3`. Environment variables and command-line assignments have the same
meaning. Use `PORTABLE=1` for binaries distributed to other CPU types.

```bash
make CC=clang full cli tools                  # Select compiler
make PREFIX=/opt/libalea install             # Install prefix
make DESTDIR=/tmp/pkg PREFIX=/usr install    # Package/stage install
make TINYPAR_BACKEND=native RELEASE=1 full cli tools
make cluster USE_MPI=1                        # Optional MPI cluster module
make test-cluster USE_MPI=1                   # Two-rank cluster test
```

`make install` builds and installs the static libraries, public headers, `alea`
CLI, tools, README, and license files. Use `install-libs`, `install-cli`, or
`install-tools` to install only one part.

### Optional cluster module

The cluster module splits work across multiple processes, on one computer or
several, using MPI. It supports volume estimation, ray tracing, rendering,
slices, mesh sampling, and geometry checks. It is optional: the rest of libalea
works without MPI.

To try it, install MPI with `mpicc` and `mpiexec` available, then build and run
the volume example with four processes:

```bash
make cluster modules USE_MPI=1
make -C examples/c cluster_volumes
mpiexec -n 4 examples/c/cluster_volumes \
    --rays 1000000 --radius 250 --center 0 0 0 model.inp
```

Replace `model.inp` with your MCNP input or OpenMC geometry XML file. `--rays`
sets the total number of sampling rays shared by all processes. Set `--center`
and `--radius` so the sampling sphere encloses the cells you want to measure.
Add `--csv -o volumes.csv` to save the results as CSV.

The volume command writes startup stages for each rank and progress after each
global batch to stderr, leaving stdout/CSV for results. Progress includes elapsed
time, rays per second, estimated time to the ray limit, unsampled instance count,
and the largest sampled relative error. Updates wait for every rank; a slow ray
or rank can delay them. For an initial diagnostic run, use `--rays 1000 --batch 100`
and set `--workers` to the CPUs allocated per rank. Keep the stage messages when
reporting a stall: they distinguish input reading, parsing, instance enumeration,
and cache preparation/sampling. These messages are not a timed heartbeat.

`--target-rel-error` requires every instance to have a measurable uncertainty.
Unsampled instances (including unreachable or zero-volume instances) prevent
early convergence, so the run may consume its entire ray budget. A warning is
printed if the requested error was not reached. One ray cannot establish an
uncertainty. The sampling sphere must enclose the intended geometry; increasing
it unnecessarily makes small instances harder to sample.

For your own C program, include [`alea_cluster.h`](include/alea_cluster.h) and
link `bin/libalea_cluster.a` with the core library. Each process loads the same
model and calls the cluster operations in the same order. See
[`cluster_volumes.c`](examples/c/cluster_volumes.c) for a complete example and
the included headers for API details.

Set `ALEA_NUM_THREADS` to the number of CPU threads available per process.
Volume estimation limits dense worker scratch to 256 MiB per rank by default;
use `--worker-memory-mib` to set a different per-rank limit. If one worker does
not fit, the command fails before sampling. The effective worker count can be
lower than `--workers` when the memory limit requires it. The command reports
the cached path table, fixed dense rank arrays, and worker scratch before
sampling.
For local development without MPI, build with `make cluster USE_MPI=0`.

#### Linux

Portable serial build:

```bash
make full cli tools
make test
```

Threaded build with GCC or Clang:

```bash
make clean
make TINYPAR_BACKEND=native RELEASE=1 full cli tools
make TINYPAR_BACKEND=native test
```

#### macOS

Portable serial build:

```bash
make full cli tools
make test
```

Threaded build:

```bash
make clean
make TINYPAR_BACKEND=native RELEASE=1 full cli tools
make TINYPAR_BACKEND=native test
```

#### Windows with MinGW/UCRT

Use the MSYS2 UCRT64 environment with the MinGW-w64 GCC toolchain. Install the
build tools from an MSYS2 shell:

```bash
pacman -S --needed mingw-w64-ucrt-x86_64-gcc mingw-w64-ucrt-x86_64-make make
```

Open a UCRT64 shell, clone the repository with submodules, and build with the
GNU Makefile:

MinGW/MSYS2 and `OS=Windows_NT` are detected automatically: TinyPar uses
Win32 threads and executables use the `.exe` suffix. For cross-compilation
from a Unix host, pass `WINDOWS_GNU=1` explicitly.

```bash
make full cli tools
make test-unit test-integration
```

Threaded build:

```bash
make clean
make TINYPAR_BACKEND=native RELEASE=1 full cli tools
make TINYPAR_BACKEND=native test-unit test-integration
```

#### Windows with conda clang-cl

Use this path when users can install conda packages and the machine has
Visual Studio 2019+ or Build Tools with the C++ x64 toolset and a Windows SDK.
Open an x64 Native Tools Command Prompt for VS before activating conda; it
sets the include and library paths for the MSVC/CRT and SDK files. Conda
provides the compiler, archive tool, and `jom`, but `clang-cl` still targets
the MSVC ABI.

```powershell
conda create -n libalea-clang -c conda-forge clang_win-64 jom
conda activate libalea-clang
where.exe clang-cl
where.exe llvm-lib
where.exe jom

jom /J 1 /f Makefile.msvc WINDOWS_CLANG_CL=1 TINYPAR_BACKEND=native full
jom /J 1 /f Makefile.msvc WINDOWS_CLANG_CL=1 TINYPAR_BACKEND=native test
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
make WINDOWS_GNU=1 CC="$cc" AR="$ar" TINYPAR_BACKEND=native full cli tools
make WINDOWS_GNU=1 CC="$cc" AR="$ar" TINYPAR_BACKEND=native test-unit test-integration test-lua
```

#### Windows with MSVC

Use Visual Studio 2019 or newer with the C++ x64 toolset installed. The wrapper
scripts locate Visual Studio, enter the x64 developer environment, and run
`nmake /f Makefile.msvc`. The MSVC build currently covers the static libraries
and tests; the Lua CLI target is built by the GNU Makefile.

Portable serial build from PowerShell:

```powershell
.\build-msvc.ps1 full
.\build-msvc.ps1 test
```

Threaded build from PowerShell:

```powershell
.\build-msvc.ps1 TINYPAR_BACKEND=native RELEASE=1 full
.\build-msvc.ps1 TINYPAR_BACKEND=native test
```

The same commands are available from `cmd.exe`:

```bat
build-msvc.bat full
build-msvc.bat TINYPAR_BACKEND=native RELEASE=1 full
```

No separate threading runtime is required.

### Dependencies

- C11 compiler (gcc, clang, or MSVC)
- `make` on Linux/macOS/MinGW, or Visual Studio `nmake` on Windows
- Standard math library (`-lm`)
- Vendored tinypar backend; native threaded builds use POSIX or Win32 threads

### Libraries Produced

| Library | Contents |
|---------|----------|
| `libalea.a` | Core engine: CSG evaluation, primitives, raycast, slice, 3D render, mesh export |
| `libalea_mcnp.a` | MCNP parser, converter, and exporter |
| `libalea_openmc.a` | OpenMC XML parser, converter, and exporter |
| `libalea_serpent.a` | Serpent exporter |
| `libalea_nucdata.a` | Nuclear data: ACE reader, cross-section lookup, free-path/nuclide/reaction sampling, Doppler broadening, multigroup collapse |
| `libalea_transport.a` | Geometry/material bindings and fixed-source neutron/photon histories; depends on nuclear data and core geometry |
| `libalea_full.a` | Core, MCNP, OpenMC, Serpent, nuclear data, and transport in one archive |

Transport applications include `alea_transport.h` and `alea_source.h` and link in dependency order:

```sh
gcc -o myapp myapp.c -Iinclude bin/libalea_transport.a bin/libalea_nucdata.a bin/libalea.a -lm -pthread
```

The transport module supports sampled fixed-source neutron and coupled
neutron-photon histories, reusable spatial sources, and optional cell,
terminal-universe, and Cartesian-mesh tallies. Geometry/material bindings are
declared in `alea_transport.h`; nuclear-data preparation and collisions remain
in `alea_nucdata.h`. See the [fixed-source workflow](docs/FIXED_SOURCE_WORKFLOW.md)
for setup and normalization, and [transport tallies](docs/TRANSPORT_TALLIES.md)
for scoring and filter semantics.

For detector-driven multigroup calculations, use `alea_adjoint.h` and
`alea_adjoint_run()`. The current adjoint kernel is steady-state and isotropic,
supports separate neutron or photon runs, and accepts vacuum and specular
boundaries. See the one-group
[`adjoint_neutron.c`](examples/c/adjoint_neutron.c) example. Continuous-energy
adjoint sampling and coupled neutron-photon adjoint transport are not yet
supported.

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

For a standard Lua 5.5 interpreter, build `make lua-module` and use
`local alea = require("alea")`. The optional `alea_cluster` module has both a
one-process development backend and an MPI backend; build and test it with
`make test-lua-cluster USE_MPI=0` or `USE_MPI=1`.

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

Lua geometry indices are zero-based and can be passed directly between APIs;
this includes cell, surface, and material indices returned by construction and
inspection functions. Domain identifiers such as MCNP cell IDs are preserved
without adjustment. Positions within Lua sequences, including ray segments,
slice curves, material components, and multigroup entries, are one-based.
CSG nodes are opaque objects tied to their originating System and cannot be
combined with nodes from another System.

## Tools

Tools are built via `make tools`:

| Tool | Description |
|------|-------------|
| `mc_convert` | Convert between MCNP, OpenMC, and Serpent geometry formats |
| `mc_plotter` | Render 2D cross-section slices of CSG geometry to PNG/BMP |
| `nuc_plot` | Generate SVG plots of nuclear cross sections, angular distributions, fission spectra, and more |
| `nuc_inventory` | Report which tables in an xsdir can be decoded, evaluated, and sampled |
| `large_model_probe` | Inspect large MCNP models and benchmark hierarchy query/raycast behavior |

```bash
bin/mc_convert model.inp model.xml
bin/mc_convert model.inp model.serp --output-format serpent
bin/mc_plotter model.inp Z 0 -100 100 -100 100 800x800 output.png
bin/nuc_plot --xsdir /path/to/xsdir --zaid 92235.80c --plot xs --output u235.svg
bin/nuc_inventory /path/to/xsdir 92235.80c
bin/large_model_probe model.inp --queries 10000 --hier-build
```

## Sampling

The tracked public sampling API is in `include/alea_nucdata.h` and `include/alea_mesh.h`:

| API | What it samples |
|-----|-----------------|
| `alea_nuc_sample_distance` | Distance to the next collision from macroscopic total cross section |
| `alea_nuc_sample_nuclide` | Target nuclide in a material |
| `alea_nuc_sample_reaction` | Reaction MT on a selected nuclide |
| `alea_nuc_sample_energy_angle_distribution` | Correlated outgoing energy and angle for supported neutron laws |
| `alea_nuc_prepare_material` / `alea_nuc_evaluate` / `alea_nuc_collide` | Capability-checked neutron or photoatomic mixture collisions |
| `alea_nuc_collide_with_secondaries` | Neutron, delayed-neutron, and photon emission using a caller-owned particle buffer |
| `alea_nuc_evaluate_urr` / `alea_nuc_urr_factors` | Coordinated evaluation or standalone unresolved-resonance factors |
| `alea_nuc_sample_photon_collision` | Coherent, bound-electron Compton, photoelectric, or pair interaction on one element |
| `alea_nuc_xs_photon_production_total` | Aggregate neutron-induced photon-production cross section |
| `alea_nuc_photon_production_audit` | Native-grid GPD versus decoded-channel consistency report |
| `alea_nuc_sample_thermal_collision` | Discrete ACE bound thermal elastic or inelastic collision |
| `alea_nuc_mg_sample_scatter` | Outgoing multigroup scatter group |
| `alea_mesh_sample` / `alea_mesh_visit` | Fixed structured-grid composition estimates |
| `alea_adaptive_grid_sample` | Nonconforming adaptive octree voxels |

Build `make mesh-benchmark` (or `make TINYPAR_BACKEND=native mesh-benchmark`) to measure
uniform, interface, and many-material grids across center, regular subcell, and
stratified modes. The benchmark also reports retained bytes per voxel for
material-only, diagnostic, and complete result masks.

See [nuclear-data transport capabilities](docs/NUCDATA_CAPABILITIES.md) for the
supported physics, ownership contract, units, and restricted slab example.

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

### Verified slice error scans

For bounded gap and overlap analysis, `alea_slice_error_query_create()` splits
a required slice rectangle into independently runnable pages. Each page reports
verified defective regions and boundary intervals or circles, confirmed point
witnesses, unresolved areas, and a receipt describing completeness and resource
use. Pages can be run individually or in a memory-bounded parallel batch with
`alea_slice_error_query_run_pages()`.

Certification currently requires a coordinate-aligned slice and a supported
combination of planes, spheres, rectangular-lattice seams, and hierarchy
occurrences. Unsupported or numerically ambiguous areas are returned as
unresolved rather than silently treated as clean. Check both `scope_classified`
and `output_complete` in every page receipt before interpreting the absence of
reported defects as a clean result. The Python binding exposes the same workflow
as `System.slice_error_query()` and `System.slice_error_page()`; see the
[`pyalea` README](bindings/python/README.md) and
[`alea_geo_validator.h`](include/alea_geo_validator.h) for the complete API and
evidence semantics.

## Documentation

| Document | Audience | Purpose |
|----------|----------|---------|
| [Tutorial](docs/TUTORIAL.md) | New users | C API walk-through from loading a model to exporting results |
| [Lua Tutorial](docs/LUA_TUTORIAL.md) | New users | Lua API for scripting and interactive use |
| [Concepts](docs/CONCEPTS.md) | All users | Surfaces, sense, cells, universes, lattices, and other domain concepts |
| [Architecture](docs/ARCHITECTURE.md) | Contributors | Internal data model, algorithms, and design decisions |
| [API Reference](docs/API.md) | All users | Every public function, grouped by task |

Start with the **Tutorial** (C) or **Lua Tutorial** if you're new. Refer to **Concepts** when something doesn't behave as you expect. The **API Reference** is for when you know what you want but forgot the function name.

Build and preview the documentation site locally with:

```bash
python3 -m venv .venv-docs
.venv-docs/bin/pip install -r requirements-docs.txt
.venv-docs/bin/mkdocs serve
```

Use `.venv-docs/bin/mkdocs build --strict` for the same strict build intended
for publication. Only the public documents listed in `mkdocs.yml` are included;
other local notes under `docs/` are excluded.

The documentation website is published at
<https://giovanni-mariano.github.io/libalea.c/> by the **Documentation** GitHub
Actions workflow. Documentation changes in pull requests targeting `main` are
built with `--strict`; changes pushed to `main` are built and deployed.

To enable publication, open the repository's **Settings → Pages → Build and
deployment** and select **GitHub Actions** as the source. Then push the workflow
to `main`, or select **Actions → Documentation → Run workflow** on `main`.
The workflow uses GitHub's built-in token; no personal access token or
`gh-pages` branch is needed.

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
tools/               mc_convert, mc_plotter, nuc_plot, nuc_inventory
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
