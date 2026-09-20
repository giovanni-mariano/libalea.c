<!--
SPDX-FileCopyrightText: 2026 Giovanni MARIANO

SPDX-License-Identifier: MPL-2.0
-->

# pyAlea native bindings

This directory builds the native Python interface to libalea. The package is
named `pyalea`; its implementation extension is the private module
`pyalea._alea`. It uses only headers from `include/` and statically links the
complete libalea archive.

The release archives are intended for applications to bundle inside their
platform-specific distributions. They are tied to a CPython ABI,
operating system, architecture, and the platform baseline used by the release
builder. Installing NumPy is still required at runtime.

One archive is released per CPython ABI and platform, covering CPython
3.10–3.14. Separate archives are not needed for each NumPy minor release. The
extension is built against the oldest practical NumPy 2.x release for each
Python version: 2.0.2 for Python 3.10–3.12, 2.1.3 for Python 3.13, and 2.3.3
for Python 3.14. The Python 3.10–3.12 artifacts are also tested at runtime with
NumPy 1.26.4 on every release platform. Every artifact is tested against the
latest compatible NumPy 2.x release without rebuilding the extension.
Consumers should currently declare `numpy>=1.26.4,<3`; Python's package
resolver will select a newer NumPy where an older release does not support that
Python version.

## Unix development build

Install NumPy and pytest for the selected interpreter, then run:

```sh
make -C bindings/python PYTHON=python3
make -C bindings/python test PYTHON=python3
make -C bindings/python package PYTHON=python3
```

`pytest` is needed only to run the binding tests. Released applications require
NumPy but do not require pytest.

Build products remain below `bindings/python/build/`; the ordinary repository
`build/` and `bin/` directories are not modified.

## Windows MSVC build

From an x64 Native Tools Command Prompt:

```bat
cd bindings\python
nmake /f Makefile.msvc PYTHON=C:\path\to\python.exe
nmake /f Makefile.msvc PYTHON=C:\path\to\python.exe test
nmake /f Makefile.msvc PYTHON=C:\path\to\python.exe package
```

The Python installation and NumPy must match the target architecture.

## Artifact layout

`make package` writes a `.tar.gz` archive on Unix or `.zip` on Windows, plus a
SHA-256 checksum. GitHub Releases also include `pyalea-manifest.json`, which
lists every binding archive with its size and digest. Extracting an archive
produces:

```text
pyalea-<version>-<python>-<abi>-<platform>/
    pyalea/
        __init__.py
        _alea<extension suffix>
        _build_info.json
    LICENSES/
    README.md
```

Add the extracted top-level directory to `PYTHONPATH`, or copy its `pyalea/`
directory into an application's wheel staging area. Applications should import
the public package:

```python
import pyalea

system = pyalea.System("example")
```
