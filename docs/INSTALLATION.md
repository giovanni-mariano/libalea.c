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

Windows MSVC builds use `build-msvc.ps1` or `build-msvc.bat` from a Visual
Studio x64 developer environment. The repository
[README](https://github.com/giovanni-mariano/libalea.c#building-from-source)
contains the complete platform-specific commands and packaging options.

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
