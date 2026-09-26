<!--
SPDX-FileCopyrightText: 2026 Giovanni MARIANO

SPDX-License-Identifier: MPL-2.0
-->

# Python bindings

`pyalea` is the native Python interface to ALEA. It statically links the
complete library and provides typed access to geometry construction, format
I/O, inspection, rendering, validation, nuclear data, and transport.

## Install a release archive

Download the archive matching the CPython ABI, operating system, and
architecture from [GitHub Releases](https://github.com/giovanni-mariano/libalea.c/releases).
Extract it and either add its top-level directory to `PYTHONPATH` or bundle the
contained `pyalea/` directory with the application. NumPy is required at
runtime; supported releases currently use `numpy>=1.26.4,<3`.

```python
import pyalea

system = pyalea.load_mcnp("geometry.inp")
system.build_universe_index()
```

Release archives are platform-specific application artifacts rather than
universal wheels. Their accompanying checksum and `pyalea-manifest.json` can
be used to verify a download.

## Build locally

Install NumPy and pytest for the selected interpreter, then run:

```bash
make -C bindings/python PYTHON=python3
make -C bindings/python test PYTHON=python3
make -C bindings/python package PYTHON=python3
```

Build products remain under `bindings/python/build/`. On Windows, use
`Makefile.msvc` from an x64 Visual Studio developer prompt.

## Paged geometry validation

One native query can analyze every page of a required slice rectangle:

```python
with system.slice_error_query(
    origin, normal, up, view_bounds, required_bounds,
    tile_columns=8, tile_rows=8,
) as query:
    for index in range(query.page_count):
        page = query.run_page(index)
        receipt = page["receipt"]
        if not receipt["scope_classified"] or not receipt["output_complete"]:
            handle_incomplete_page(page)
        consume(page["regions"], page["intervals"], page["circles"])
```

Pages can be requested in any order and retain Python-owned snapshots after
the query closes. `System.slice_error_page()` is a convenience method for a
single page. See [Geometry validation](GEOMETRY_VALIDATION.md) for the evidence
and completeness rules.

The repository's
[binding README](https://github.com/giovanni-mariano/libalea.c/tree/main/bindings/python)
documents artifact contents, primitive parameter layouts, and additional
development details.
