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
builder. Installing NumPy is still required at runtime. The package includes a
`py.typed` marker and a public type stub; native names remain available through
`pyalea._alea`, while `pyalea.__all__` is explicitly maintained.

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

## Programmatic fill transforms

Named and inline transforms are separate because named transforms replace a
caller-selected ID, while inline transforms are deduplicated and receive an
automatically assigned ID. Both return the ID accepted by `System.set_fill()`:

```python
transform_id = system.add_transform(10, (5.0, 0.0, 0.0))
system.set_fill(cell_index, universe_id, transform=transform_id)

inline_id = system.add_inline_transform(
    (5.0, 0.0, 0.0), cell_id=cell_id, role="fill"
)
```

Set `degrees=True` when rotation entries are MCNP angles instead of direction
cosines. Transform sequences use the MCNP displacement, optional rotation, and
optional trailing direction flag accepted by libalea.

Models can also be exported without a temporary path:

```python
mcnp_text = system.export_mcnp_string()
openmc_xml = system.export_openmc_string()
serpent_text = system.export_serpent_string()
```

## Paged slice error analysis

Use one native query for all pages of a required validation rectangle:

```python
with system.slice_error_query(
    origin, normal, up, view_bounds, required_bounds,
    tile_columns=8, tile_rows=8,
) as query:
    for index in range(query.page_count):
        page = query.run_page(index)
        consume(page["receipt"], page["intervals"],
                page["regions"], page["unresolved"])
```

Pages may be requested in any order. Their receipts share a process-local
`query_id`; each result remains a Python-owned snapshot after the query closes.
The query keeps the system alive and rejects pages after a geometry change.
Keep the system unchanged while a page is running. The query holds Python's
GIL during each bounded page call to prevent concurrent Python mutation.
`System.slice_error_page()` remains a convenience call for one page and creates
a fresh native query each time.

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

Standalone primitives can be evaluated without creating a `System`:

```python
value = pyalea.primitive_evaluate(
    pyalea.PRIMITIVE_SPHERE,
    (0.0, 0.0, 0.0, 2.0),
    (1.0, 0.0, 0.0),
)
```

The result is a signed implicit value, with negative values on the primitive's
interior side and zero on its boundary. Its magnitude is not generally a
Euclidean distance. Parameter sequence lengths follow the corresponding
public `alea_primitive_data_t` geometry fields. Incorrect sequence lengths,
degenerate geometry, non-finite parameters or point coordinates, and unsupported
primitive types raise `ValueError`. Non-sequence inputs or elements that cannot
be converted to numbers raise `TypeError`. If finite inputs produce a non-finite
native result, evaluation raises `ArithmeticError`.

The accepted parameter layouts are:

| Primitive constants | Parameters |
|---|---|
| `PRIMITIVE_PLANE` | `a, b, c, d` for `ax + by + cz + d = 0` |
| `PRIMITIVE_SPHERE`, `PRIMITIVE_SPH` | `cx, cy, cz, radius` |
| `PRIMITIVE_CYLINDER_X` | `center_y, center_z, radius` |
| `PRIMITIVE_CYLINDER_Y` | `center_x, center_z, radius` |
| `PRIMITIVE_CYLINDER_Z` | `center_x, center_y, radius` |
| `PRIMITIVE_CONE_X/Y/Z` | `apex_x, apex_y, apex_z, tan_angle_sq, sheet` |
| `PRIMITIVE_RPP` | `xmin, xmax, ymin, ymax, zmin, zmax` |
| `PRIMITIVE_QUADRIC` | `A, B, C, D, E, F, G, H, I, J` |
| `PRIMITIVE_TORUS_X/Y/Z` | `cx, cy, cz, major_radius, minor_radius` |
| `PRIMITIVE_RCC` | base (3), height (3), radius |
| `PRIMITIVE_BOX` | corner (3), then edge vectors v1, v2, v3 (3 each) |
| `PRIMITIVE_TRC` | base (3), height (3), base radius, top radius |
| `PRIMITIVE_ELL` | first focus (3), second focus (3), full major-axis length |
| `PRIMITIVE_REC` | base (3), height (3), ellipse axes a1 and a2 (3 each) |
| `PRIMITIVE_WED` | vertex (3), then edge vectors v1, v2, v3 (3 each) |
| `PRIMITIVE_RHP` | base (3), height (3), radial vectors r1, r2, r3 (3 each) |

Cone `sheet` is `-1`, `0`, or `1`. RHP permits the regular-hex shorthand in
which r2 and r3 are zero vectors. `PRIMITIVE_ARB` is not accepted by this
sequence interface.

The registered cone constructors also accept an optional final positional
argument `sheet`, defaulting to `0`:

```python
surface_index, positive_node, negative_node = system.cone_z_surface(
    10, 0.0, 0.0, 0.0, 1.0, 1,
)
```

The arguments are `surface_id, apex_x, apex_y, apex_z, tan_angle_sq, sheet`;
`cone_x_surface()` and `cone_y_surface()` use the same layout. A sheet of `1`
selects the positive axis direction from the apex, `-1` the negative direction,
and `0` both sheets. Omitting the final argument preserves the two-sheet
behavior. These constructors return the surface index and its positive and
negative halfspace node IDs, and raise `RuntimeError` if native surface creation
fails, including an invalid sheet selection.
