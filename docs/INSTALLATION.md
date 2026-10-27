<!--
SPDX-FileCopyrightText: 2026 Giovanni MARIANO

SPDX-License-Identifier: MPL-2.0
-->

# Installation

## Release archives

Pre-built archives are available from
[GitHub Releases](https://github.com/giovanni-mariano/libalea.c/releases) for
Linux, macOS, and Windows on the supported x86-64 and Arm platforms. The Unix
and MinGW archives contain the command-line tools, static libraries, and public
headers. MSVC archives contain static libraries and headers.

Native `pyalea` archives are published separately for each supported CPython
ABI and platform. See [Python bindings](PYTHON.md) for their layout and runtime
requirements.

## Build from source

A C11 compiler and GNU Make are required on Linux, macOS, and MinGW. Clone the
repository with its vendored dependencies:

```bash
git clone --recursive https://github.com/giovanni-mariano/libalea.c.git
cd libalea.c
```

Build only what the application needs:

```bash
make                 # Core library
make modules         # Format, nuclear-data, and transport modules
make full            # Combined static library
make cli tools       # Interactive CLI and command-line tools
make test            # Build and run the test suite
```

The default TinyPar backend uses native threads. Select an explicitly serial
build with `TINYPAR_BACKEND=serial`, or add `RELEASE=1` for optimization:

```bash
make TINYPAR_BACKEND=native RELEASE=1 full cli tools
```

Use `ALEA_NUM_THREADS` or `alea_parallel_set_threads()` to limit the process-wide
worker pool. Applications distributing binaries to other CPU types should add
`PORTABLE=1` to release builds.

Install under `/usr/local` by default, or select another prefix:

```bash
make install
make PREFIX=/opt/libalea install
```

Windows MSVC builds use `build-msvc.ps1` or `build-msvc.bat`.
See the platform-specific instructions below for toolchain setup.

The experimental [GPU rendering and batch compute module](GPU_COMPUTE.md) is built separately:
`make compute` for CPU, or `make compute USE_METAL=1` for CPU and Metal on macOS.
It is not included in the default or full archives, or the pre-built releases.
For the interactive Metal viewer, use `make tools USE_METAL=1 USE_SDL=1` and
`bin/render3d model.i --backend metal --interactive`. See the module guide for
supported geometry and rendering features.

## Link an application

Libraries must precede their dependencies. For example, an MCNP application
links the parser before the core library:

```bash
cc -o myapp myapp.c -Iinclude bin/libalea_mcnp.a bin/libalea.a -lm
```

The combined archive is convenient when several modules are needed:

```bash
cc -o myapp myapp.c -Iinclude bin/libalea_full.a -lm
```

Continue with the [C tutorial](TUTORIAL.md) or [Lua tutorial](LUA_TUTORIAL.md).

## Build options

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

## Platform-specific builds

### Linux

Default threaded build:

```bash
make full cli tools
make test
```

Optimized build with GCC or Clang:

```bash
make clean
make TINYPAR_BACKEND=native RELEASE=1 full cli tools
make TINYPAR_BACKEND=native test
```

### macOS

Default threaded build:

```bash
make full cli tools
make test
```

Optimized build:

```bash
make clean
make TINYPAR_BACKEND=native RELEASE=1 full cli tools
make TINYPAR_BACKEND=native test
```

### Windows with MinGW/UCRT

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

Optimized build:

```bash
make clean
make TINYPAR_BACKEND=native RELEASE=1 full cli tools
make TINYPAR_BACKEND=native test-unit test-integration
```

### Windows with conda clang-cl

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

### Windows with conda MinGW/UCRT

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

### Windows with MSVC

Use Visual Studio 2019 or newer with the C++ x64 toolset installed. The wrapper
scripts locate Visual Studio, enter the x64 developer environment, and run
`nmake /f Makefile.msvc`. The MSVC build currently covers the static libraries
and tests; the Lua CLI target is built by the GNU Makefile.

Default threaded build from PowerShell:

```powershell
.\build-msvc.ps1 full
.\build-msvc.ps1 test
```

Optimized build from PowerShell:

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

## Optional cluster module

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

For your own C program, include [`alea_cluster.h`](https://github.com/giovanni-mariano/libalea.c/blob/main/include/alea_cluster.h) and
link `bin/libalea_cluster.a` with the core library. Each process loads the same
model and calls the cluster operations in the same order. See
[`cluster_volumes.c`](https://github.com/giovanni-mariano/libalea.c/blob/main/examples/c/cluster_volumes.c) for a complete example and
the included headers for API details.

Set `ALEA_NUM_THREADS` to the number of CPU threads available per process.
Volume estimation limits dense worker scratch to 256 MiB per rank by default;
use `--worker-memory-mib` to set a different per-rank limit. If one worker does
not fit, the command fails before sampling. The effective worker count can be
lower than `--workers` when the memory limit requires it. The command reports
the cached path table, fixed dense rank arrays, and worker scratch before
sampling.
For local development without MPI, build with `make cluster USE_MPI=0`.

## Libraries produced

| Library | Contents |
|---------|----------|
| `libalea.a` | Core engine: CSG evaluation, primitives, raycast, slice, 3D render, mesh export |
| `libalea_mcnp.a` | MCNP parser, converter, and exporter |
| `libalea_openmc.a` | OpenMC XML parser, converter, and exporter |
| `libalea_xml.a` | Native ALEA model XML importer and exporter |
| `libalea_serpent.a` | Serpent exporter |
| `libalea_nucdata.a` | Nuclear data: ACE reader, cross-section lookup, free-path/nuclide/reaction sampling, Doppler broadening, multigroup collapse |
| `libalea_transport.a` | Geometry/material bindings and fixed-source neutron/photon histories; depends on nuclear data and core geometry |
| `libalea_full.a` | Core, MCNP, OpenMC, ALEA XML, Serpent, nuclear data, and transport in one archive |

Transport applications include `alea_transport.h` and `alea_source.h` and link in dependency order:

```sh
gcc -o myapp myapp.c -Iinclude bin/libalea_transport.a bin/libalea_nucdata.a bin/libalea.a -lm -pthread
```

The transport module supports sampled fixed-source neutron and coupled
neutron-photon histories, reusable spatial sources, and optional cell,
terminal-universe, and Cartesian-mesh tallies. Geometry/material bindings are
declared in `alea_transport.h`; nuclear-data preparation and collisions remain
in `alea_nucdata.h`. See the [fixed-source workflow](FIXED_SOURCE_WORKFLOW.md)
for setup and normalization, and [transport tallies](TRANSPORT_TALLIES.md)
for scoring and filter semantics.

For detector-driven multigroup calculations, use `alea_adjoint.h` and
`alea_adjoint_run()`. The current adjoint kernel is steady-state and isotropic,
supports separate neutron or photon runs, and accepts vacuum and specular
boundaries. See the one-group
[`adjoint_neutron.c`](https://github.com/giovanni-mariano/libalea.c/blob/main/examples/c/adjoint_neutron.c) example. Continuous-energy
adjoint sampling and coupled neutron-photon adjoint transport are not yet
supported.

Link against the core library plus the format modules you need:

```bash
# MCNP support
gcc -o myapp myapp.c -Iinclude bin/libalea_mcnp.a bin/libalea.a -lm

# OpenMC support
gcc -o myapp myapp.c -Iinclude bin/libalea_openmc.a bin/libalea.a -lm

# Native ALEA XML model support
gcc -o myapp myapp.c -Iinclude bin/libalea_xml.a bin/libalea.a -lm

# Nuclear data
gcc -o myapp myapp.c -Iinclude bin/libalea_nucdata.a bin/libalea.a -lm

# Serpent export
gcc -o myapp myapp.c -Iinclude bin/libalea_serpent.a bin/libalea.a -lm

# Full library (core + all formats + nucdata)
gcc -o myapp myapp.c -Iinclude bin/libalea_full.a -lm
```
