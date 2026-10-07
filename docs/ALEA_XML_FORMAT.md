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
      <parameters>
        <volume value="4188.790204786391" />
        <importance particle="photon" value="0" />
      </parameters>
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

Cells preserve `name`, leading `comments`, `inline_comment`, universe and fill
IDs, and fill transforms. Supplementary physical and transport settings belong
in one optional `<parameters>` container:

```xml
<cell id="1" material="0" region="-1">
  <parameters>
    <volume value="12.5" />
    <temperature value="600" units="K" />
    <importance particle="neutron" value="0" />
    <photon_production mode="threshold" weight_threshold="0.5"
                       weight_basis="source_relative" />
    <fission_mode value="capture_without_photons" />
    <detector_contribution_probability value="1" />
    <detector_contribution_probability tally="5" value="0.1" />
    <energy_cutoff particle="electron" value="0.001" units="MeV" />
    <secondary_collision_state particle="photon" value="collided" />
    <magnetic_field ref="2" />
  </parameters>
</cell>
```

The container may be empty or omitted. Parameter order has no meaning.
Singleton parameters cannot be repeated; particle-specific entries must be
unique per particle, and detector probabilities must be unique per tally.
Unknown elements, attributes, enum values, nonfinite numbers, and duplicate
entries are rejected. Supported particle names are `neutron`, `photon`, and
`electron`.

| Element | Attributes and semantics | Omitted setting |
|---|---|---|
| `volume` | Nonnegative `value`, in cm³; supplied volume for tally normalization, independent of computed geometry | No supplied volume |
| `temperature` | Nonnegative `value`; optional `units="K"` | No supplied temperature |
| `importance` | Required `particle` and nonnegative `value` | **1 for each particle** |
| `photon_production` | Required `mode`; threshold attributes described below | Threshold 1, relative to the source neutron weight |
| `fission_mode` | `value="normal|capture_with_photons|capture_without_photons"` | `normal` |
| `detector_contribution_probability` | `value` in [0,1]; optional `tally` | 1, unless overridden by the all-tally default |
| `energy_cutoff` | Nonnegative `value`; optional `particle` and `units="MeV"` | Use the global transport cutoff |
| `secondary_collision_state` | `value="collided|uncollided"`; optional `particle` | `uncollided`, unless overridden by the all-particle default |
| `magnetic_field` | Nonnegative integer `ref`; zero explicitly means no field | No field |

An omitted importance has an effective value of 1. Explicit `value="0"`
is retained and overrides that default. Explicit settings equal to defaults
also retain their presence through ALEA XML round trips.

Without `particle`, an energy cutoff or secondary state defines a default for
all supported particles in the cell. A particle-specific entry overrides that
cell default. A cell energy cutoff supplies a lower bound; the higher of the
cell cutoff and the global transport cutoff applies. Omitting `tally` (or
using `tally="0"`) defines the detector probability default for all tallies;
a positive tally ID selects a specific override. References preserve external
IDs; this format currently does not define magnetic fields or detector tallies.

### Photon production

`photon_production` controls neutron-induced photon production:

- `mode="threshold"` requires a positive `weight_threshold` and a
  `weight_basis="absolute|source_relative"`. Both bases include neutron
  importance scaling between the source and collision cells; `source_relative`
  additionally scales by the source neutron's starting weight. Photons below
  the resulting threshold undergo Russian roulette.
- `mode="one_per_collision"` produces one photon per neutron collision when
  photon production is possible.
- `mode="off"` disables neutron-induced photon production in the cell.

`weight_threshold` and `weight_basis` are accepted only in threshold mode.
The native representation uses a positive threshold and a descriptive mode,
without MCNP's sign conventions or sentinel values.

### MCNP conversion

MCNP conversion maps the normalized settings as follows:

| ALEA setting | MCNP |
|---|---|
| Photon threshold, `absolute` | Positive `PWT` |
| Photon threshold, `source_relative` | Negative `PWT` |
| Photon mode `one_per_collision` | `PWT=0` |
| Photon mode `off` | `PWT=-1000000` |
| Fission mode `normal` | `NONU=1` |
| Fission mode `capture_with_photons` | `NONU=0` |
| Fission mode `capture_without_photons` | `NONU=2` |
| Detector probability | `PD0` default and `PDn` tally overrides |
| Energy cutoff | `ELPT:N`, `ELPT:P`, or `ELPT:E` |
| Secondary state | `UNC:N`, `UNC:P`, or `UNC:E`; collided=0, uncollided=1 |
| Magnetic field reference | `BFLCL` |

A source-relative threshold of exactly 1000000 is valid in native ALEA XML
but cannot be converted to MCNP because it collides with the production-off
sentinel. Conversion reports an error rather than changing its meaning.
All-particle defaults are expanded into particle-specific MCNP entries, with
explicit particle overrides taking precedence. MCNP cell-card imports support
comma-separated particle designators for `ELPT` and `UNC`, and per-tally `PDn`
entries; supported parameters also survive `LIKE/BUT` resolution.

Cell parameters must use the named elements inside `<parameters>`. Direct
`<importance>` children, generic `<parameter name="...">` entries, and the
cell `temperature` attribute are rejected.

The C model exposes `alea_photon_production_t`, `alea_fission_mode_t`, per-particle
cutoff/state arrays with presence masks, and owned per-tally detector entries.
Use `alea_model_cell_set_detector_probability()` to set probabilities. Cell
copy and removal hooks maintain those owned entries. The original parameter
flag names remain aliases for the descriptive flag names; clients accessing
model metadata structures should rebuild against the updated header.

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
