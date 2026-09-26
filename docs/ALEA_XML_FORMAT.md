<!--
SPDX-FileCopyrightText: 2026 Giovanni MARIANO

SPDX-License-Identifier: MPL-2.0
-->

# ALEA XML model format

ALEA XML is the native interchange format for a complete ALEA model. It keeps
the geometry in an `alea_system_t` and the document and cell transport metadata
in an `alea_model_t`. The conventional suffix is `.alea.xml`.

```xml
<?xml version="1.0" encoding="utf-8"?>
<alea version="1" length_units="cm" name="example">
  <materials>
    <material id="7" name="fuel" fraction_basis="atom">
      <nuclide zaid="92235" fraction="0.04" library="80c" />
      <nuclide zaid="92238" fraction="0.96" library="32c" />
    </material>
  </materials>
  <surfaces>
    <surface id="1" type="sphere" coeffs="0 0 0 10" boundary="vacuum" />
  </surfaces>
  <cells>
    <cell id="1" name="fuel cell" universe="0" material="7"
          density="10.2" density_units="g/cm3" region="-1">
      <importance particle="neutron" value="1" />
    </cell>
  </cells>
</alea>
```

## C API and ownership

Include `alea_xml.h` and link `libalea_xml.a` before `libalea.a`:

```c
alea_model_t *model = alea_xml_load("model.alea.xml");
if (!model) {
    fprintf(stderr, "%s\n", alea_error());
    return 1;
}

alea_system_t *sys = alea_model_system(model); /* borrowed */
/* geometry queries use sys */

alea_xml_export(model, "copy.alea.xml");
alea_model_destroy(model); /* also destroys the loaded system */
```

`alea_xml_load_string()` accepts a byte count; zero means a NUL-terminated
string. Stream exporters do not close the supplied stream. The convenience
`alea_xml_export_system*()` functions write geometry without model metadata.

`alea_model_wrap()` is non-owning and `alea_model_adopt()` takes ownership of
the system. A system can have only one model wrapper at a time because the
wrapper tracks cell insertions, copies, and removals.

## Document structure

The root must be exactly `alea` with `version="1"` and
`length_units="cm"`. `name`, `title`, and `comments` are optional. The sections
are `materials`, `transforms`, `surfaces`, and the required `cells` section.
References use positive integer external IDs.

The parser accepts ordinary elements, attributes, text, XML comments, the five
predefined XML escapes, and numeric character references. DTD and entity
declarations are rejected and no network or filesystem resource is resolved.

## Materials and density

A material has an `id`, `fraction_basis="atom|weight"`, and optional `name`,
`comments`, `density`, and `density_units`. Density units are `g/cm3` and
`atom/b-cm`; values in XML are always nonnegative.

Every `nuclide` has `zaid`, `fraction`, and a required string `library`
extension. The extension is serialized without a leading dot, but ALEA's
in-memory material representation uses one. Library selection is per nuclide,
so a material can deliberately mix `31c`, `32c`, `80c`, or other extensions.
An `element` similarly uses `z`, `fraction`, and `library`. A `thermal` entry
has an `identifier` and optional `zaid_match`.

A mixture has an `id`, fraction basis, optional `mc_material_id`, `name`, and
`comments`, and contains `component` entries with `material` and `fraction`.
A cell referring to one uses `material_kind="mixture"`.

Cell density is independent of a material's optional standard density. Every
non-void cell explicitly records both its effective `density` and
`density_units`. A void cell uses `material="0"` and omits them. Loading XML
does not open ACE files or prepare nuclear data.

## Surfaces and regions

Surface types and coefficient order are:

| Type | Coefficients |
|---|---|
| `plane` | `a b c d` for `ax + by + cz + d = 0` |
| `sphere`, `sph` | `cx cy cz radius` |
| `cylinder-x/y/z` | the two transverse center coordinates, then radius |
| `cone-x/y/z` | `apex-x apex-y apex-z tan-angle-squared sheet` |
| `rpp` | `xmin xmax ymin ymax zmin zmax` |
| `quadric` | `A B C D E F G H I J` |
| `torus-x/y/z` | `cx cy cz major-radius minor-radius axial-semiwidth` |
| `rcc` | base vector, height vector, radius |
| `box` | corner and three edge vectors |
| `trc` | base vector, height vector, base radius, top radius |
| `ell` | two three-vectors and major-axis length |
| `rec`, `wed` | four three-vectors |
| `rhp` | five three-vectors |
| `arb` | corner count, face count, 24 corner values, 24 face indices |

`boundary` is `transmissive` (the default), `reflective`, `white`, `periodic`,
or `vacuum`. Periodic surfaces use `periodic_surface`. Optional `transform` and
`transform_applied` retain transform provenance; coefficients represent the
stored ALEA primitive and are not transformed again during import.

The `region` grammar uses signed surface IDs, whitespace intersection, `|`
union, `~` complement, and parentheses. Precedence is complement,
intersection, then union. For example, `-1 2 | ~(-3 4)`.

An optional `bbox="xmin xmax ymin ymax zmin zmax"` is a finite, ordered cache
hint only. The importer validates it but computes geometry bounds itself. A
bounding box never clips a cell and cannot change containment.

## Cell metadata and hierarchy

Cells preserve `name`, leading `comments`, `inline_comment`, temperature in
kelvin, universe and fill IDs, and fill transforms. Importances are explicit:

```xml
<importance particle="photon" value="0" />
```

The supported typed parameters are `volume`, `pwt`, `nonu`, `pd`, `elpt`,
`unc`, and `bflcl`. Presence is distinct from a zero value:

```xml
<parameter name="volume" value="12.5" />
```

Rectangular and hexagonal lattice cells use `lattice_type`, six integer
`lattice_dims`, flattened `lattice_fill`, three-value `lattice_pitch` and
`lattice_lower_left`, and optional `lattice_outer`, `lattice_repeating`, and
`lattice_zero_coords`. Fill values use ALEA's existing x-fastest internal
ordering.

MCNP source spelling such as `LIKE/BUT` is not native model state. Conversion
preserves its resolved geometry and the supported normalized metadata instead.

## Command-line conversion

`mc_convert` accepts `alea` as an input and output format:

```bash
mc_convert -if mcnp -of alea model.inp model.alea.xml
mc_convert -if alea -of mcnp model.alea.xml model.inp
```

For input, XML files are distinguished by their root element. For output, the
`.alea.xml` suffix selects the native format when the format is not explicit.
