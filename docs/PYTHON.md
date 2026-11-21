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
runtime. Release archives target CPython 3.12–3.14 and use `numpy>=2.0,<3`.

```python
import pyalea

system = pyalea.load_mcnp("geometry.inp")
system.build_universe_index()
```

Native model XML preserves cell metadata and per-nuclide library extensions:

```python
system = pyalea.load_alea("model.alea.xml")
xml_text = system.export_alea_string()
system.export_alea("copy.alea.xml")
```

`load_alea_string()` loads the same format from memory. See the
[ALEA XML format guide](ALEA_XML_FORMAT.md) for the schema and round-trip
contract.

Release archives are platform-specific application artifacts rather than
universal wheels. Their accompanying checksum and `pyalea-manifest.json` can
be used to verify a download.

For a suspect ray crossing inside a fill, call
`system.diagnose_substep_fill(ox, oy, oz, dx, dy, dz, t)`. A finding gives
the parent cell, fill universe, entry and exit surface IDs in flight order,
rounded distances, and independently computed width. Equal rounded distances
can enclose a physical gap. The method returns `None` for an inconclusive query.

## Build locally

Install NumPy and pytest for the selected interpreter, then run:

```bash
make -C bindings/python PYTHON=python3
make -C bindings/python test PYTHON=python3
make -C bindings/python package PYTHON=python3
```

Build products remain under `bindings/python/build/`. On Windows, use
`Makefile.msvc` from an x64 Visual Studio developer prompt.

## Repeated point lookups

```python
system = pyalea.load_mcnp("geometry.inp")
owner = system.find_cell_at(x, y, z)  # (cell ID, material ID), or None
```

`find_cell_at` and `material_at` share the model's point index. They start
with the spatial hierarchy. After repeated lookups have spent an estimated
125 ms there, one caller builds the regional index. Its lookup pays the
build cost and other readers can keep using the hierarchy. The model chooses
the world box and adds child bounds on demand. A world without a usable box
continues to use the hierarchy.

Call `system.prepare_point_acceleration()` to pay the build cost up front.
Geometry and material edits discard the index and restart observation.
Editing a model requires exclusive access.

Python's `find_cell` uses the hierarchy to preserve its overlap ordering
and resolution flags.

## Repeated finite ray batches

```python
system.prepare_ray_acceleration()
segments = system.raycast_hier_fast_batch(origins, directions, t_max=1000, include_surfaces=False)
```

`raycast_hier_fast_batch(..., t_max=distance, include_surfaces=False)` also
works on a freshly loaded model. During a costly batch, it can prepare
hierarchy recovery lists while retaining the rows already traced. If you
expect sustained finite-ray work, call `prepare_ray_acceleration()` once.
Subsequent eligible calls on the same model reuse the index. The build can
delay the first result substantially.

The method derives a box from bounded world cells. To specify the working
volume, pass `bounds=(xmin, xmax, ymin, ymax, zmin, zmax)` during preparation.
Each segment reports entry and exit distance, cell, material, density and
resolution flags. The default `include_surfaces=True`, path requests and
unbounded rays use hierarchy traversal. Finish active queries before preparing
or editing. Call preparation again after edits to accelerate the updated geometry.

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
