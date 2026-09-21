<!--
SPDX-FileCopyrightText: 2026 Giovanni MARIANO

SPDX-License-Identifier: MPL-2.0
-->

# Transport source implementation plan

Status: phase 1 implemented for point/box space, isotropic/monodirectional
angle, constant energy/time/weight, prepared source reuse, and preview in
C/Python/Lua. Line, spherical-volume and cylindrical-volume sampling,
cone/cosine/tabulated-polar/radial angles, and discrete energy lines from
phase 2 are also implemented. Other phase 2 and later
capabilities remain proposed. The phase descriptions below distinguish
implemented features from the remaining work.

## Objective and terminology

Fixed-source transport is the simulation mode: an externally prescribed source,
as opposed to a criticality/eigenvalue calculation. It is not a spatial source
type. Sources describe a joint distribution of particle species, position,
direction, energy, and emission time, together with statistical weight.

The priority is fusion neutronics: simple sources for verification, followed by
cylindrical and toroidal volumes and spatially varying tokamak emissivity.
Provide the same built-in capabilities through C, Python, and Lua. Support
sampling and plotting sources without loading any xsdir or ACE data.

## Existing foundation

- `include/alea_transport.h` defines one primary particle state and a sampler
  callback receiving seed and global history ID.
- `src/transport/fixed_source.c` supports a repeated particle state and a
  callback for independently sampled primary histories. Prepared spatial and
  angular distributions live in `src/transport/source.c`.
- Python and Lua expose prepared descriptions through `space`, `angle`, and
  `energy` fields. The former combined `kind` selectors have been removed.
- `src/rng/alea_rng_distribution.*` already implements discrete CDF/alias tables,
  histogram and linear PDF sampling, and isotropic/cosine angular sampling.
- Position, direction, and energy already have separate transport RNG domains.
  Existing tests exercise replay and splitting histories into batches.

Retain the collision and tally engines. A prepared source will adapt to the
existing sampler callback; there is no need to redesign the transport loop.

## Source description and prepared representation

The canonical source description has these fields:

| Component | Initial choices | Later choices |
| --- | --- | --- |
| Particle | neutron or photon | weighted mixtures of source components |
| Space | point, line, box, sphere, cylinder; rectangle/disk/sphere/cylinder surfaces | torus, weighted mesh, tokamak emissivity, sampled CSG region |
| Angle | monodirectional, isotropic, uniform-solid-angle cone, cosine hemisphere, tabulated polar PDF, radial from sampled position | general correlated models |
| Energy | monoenergetic, discrete lines, tabulated PDF | position-dependent and validated fusion emission models |
| Time | constant, uniform interval, tabulated PDF | pulse trains and correlations |
| Weight | positive constant, default one | explicit importance-sampling correction |

Use `space`, `angle`, `energy`, and `time` sub-descriptions with a `type` tag.
For example, a point source may use isotropic or monodirectional emission;
a cylindrical volume may use either a line spectrum or tabulated energy PDF.
Require energy explicitly. Make angular defaults explicit in the documented
schema; use isotropic for new independent source descriptions.

Independent components are the default, not an architectural restriction.
Reserve a separate correlated-source variant that generates a complete state
and can use spatial context such as a plasma-region index or surface normal.

Proposed C interface, with names finalized during phase 1:

- `alea_source_spec_t`: tagged host-side configuration.
- `alea_source_t`: opaque, immutable prepared source.
- Prepare/free operations that validate and copy configuration and build tables.
- Sample one history, and sample a batch into caller-owned arrays.
- Inspect particle mask, source normalization, and supported execution backends.
- An adapter to `alea_transport_run_sampled_source()`.

Preparation owns its tables; sampling performs no allocation or configuration
parsing. A sample contains world position and direction, energy in MeV, time in
seconds, weight, and neutron/photon identity. A shared local frame handles
translations and rotations for cylinders, disks, cones, and tori. Validate
frames once and reject singular transforms. Initially support rigid transforms;
scaling belongs in each distribution's explicit dimensions.

Keep host pointers behind the opaque handle. Store built-in parameters and
sampling tables so they can later be packed into tagged records with offsets.
An arbitrary C callback remains a CPU extension; it is not automatically GPU
or WASM portable. No GPU kernels are required in this plan.

## Sampling, normalization, and failure contracts

1. Default sampling follows the normalized physical source distribution with
   weight one. Specify a physical source rate separately from statistical
   weight; existing tally means remain per sampled primary history. Scaling
   to particles/s is explicit. For a pulse, distinguish total yield from rate.
2. Source mixtures select a component according to its normalized strength.
   Do not multiply its strength into the weight a second time. If biased
   component sampling is added later, apply the physical/sampling probability
   ratio explicitly and require nonzero sampling support everywhere needed.
3. State whether a table contains probability masses or probability densities.
   Include interval widths, area, or volume when converting densities to masses.
4. Sample spatial distributions with the appropriate measure: line length,
   surface area, or volume. Uniform radius is not uniform cylindrical volume;
   uniform coordinates in a parameterized surface are not generally uniform area.
5. A surface sample supplies a local normal. Define inward/outward sides and
   distinguish uniform hemispherical emission from cosine-law emission.
   Geometry-bound surfaces must initialize the navigator on the selected side;
   avoid a universal displacement epsilon that can cross thin regions.
6. Seed plus global history ID determine a sample. Allocate distinct domains
   for component selection and time, keeping existing domain values stable.
   Rejection draws must not perturb direction, energy, other histories, or
   collision streams. Record the source schema and sampling algorithm versions.
   Do not promise bitwise identity across all future CPU/GPU math backends.
7. Bound rejection attempts. On exhaustion, return an error with history ID,
   source component, and reason. Never silently discard a source history or
   count it as a zero contribution.
8. Geometry restrictions are explicit conditioning of the source. A point
   outside the permitted region is a failure unless the caller requested a
   conditional distribution. Report acceptance statistics; do not silently
   preserve an unconditioned physical rate after restricting the source.
9. Retain one primary particle per source history initially. Correlated
   multiparticle emission events require a later history/bank contract because
   they change how tally second moments must be accumulated.

## Phase 1 — source API and migration

Implement `include/alea_source.h` and a source implementation under
`src/transport/`, with preparation, sampling, diagnostics, and a driver adapter.
Keep source sampling independent of nuclear-data loading. Reuse the existing
RNG utilities rather than adding language-specific samplers.

Initially implement point and box spatial distributions, isotropic and
monodirectional angles, and constant energy/time/weight. Use the prepared
source sampler throughout the scripting bindings and migrate examples and tests
to the component schema. Combined source selectors are removed.

Add matching prepared-source objects in Python and Lua, accepted by
`transport_run`; a configuration dictionary/table remains a convenience that
prepares a source internally. Determine required primary particle bindings from
source metadata, extending it for requested neutron/photon coupling. Source
mixtures must not depend on the old single `source.particle` lookup.

Add `sample_source` access without a geometry or xsdir. Python returns NumPy
arrays: position/direction `(N,3)`, energy/time/weight `(N,)`, and integer species
and history IDs. Lua provides equivalent fields in tables, with bounded batches
for large previews. Share sample generation with transport so previews reproduce
the same history IDs exactly. Preparation and preview do not prepare materials.

Acceptance: migrated examples work; both new point-angle combinations work;
preview and transport use the same source histories; all three languages agree.

## Phase 2 — analytic spatial and independent distributions

Add line segments, spheres/spherical shells, and cylinders/annular cylinders as
volume sources. Add finite rectangles, disks/annuli, spherical surfaces, and
cylinder side/cap surfaces as distinct spatial types. If a closed cylinder is
requested, select its patches by area. Support arbitrary orientations through
the shared frame. Do not accept an infinite emitting surface without bounds.

Cone, cosine-hemisphere, radial, and tabulated-polar angular laws and discrete
energy lines are implemented. Add tabulated
energy PDFs, and constant/uniform/tabulated emission times. Reuse prepared CDF
and alias tables. Clearly specify histogram versus linear PDF interpolation.
The tabulated polar law uses density per unit `mu = cos(theta)` and uniform
azimuth; it is independent of position. Radial direction is an explicit
position-dependent law, while general angle-position correlations remain later work.

Acceptance: analytic CDF/moment checks validate length, area, volume and solid
angle sampling; invalid dimensions, spectra, and frames fail at preparation.
Add matching binding examples in the same changes as the native features.

## Phase 3 — mixtures and tabulated fusion emissivity

Add weighted source components and a weighted Cartesian spatial mesh. Specify
whether mesh entries are emission density or integrated voxel strengths; use
the corresponding voxel volumes when preparing selection probabilities.

Prioritize axisymmetric tokamak sources represented by emissivity on an `(R,Z)`
grid and a toroidal angle interval. Start with piecewise-constant emissivity in
each bin. Prepare cell masses using cylindrical volume: for a bin, volume is
`0.5 * (R_hi^2 - R_lo^2) * (Z_hi - Z_lo) * delta_phi`. Sample radius within a
selected cell according to that volume measure, not uniformly in `R`.
Require nonnegative radii, valid extents and positive total source strength.

Add an analytic circular torus and toroidal sectors for simple fusion studies.
For the initial ring-torus model require major radius greater than minor radius.
Uniform volume must include the toroidal-coordinate Jacobian. Treat a toroidal
surface as a separate future area sampler rather than reusing the volume law.

These models accept user-supplied emissivity. They do not silently compute fusion
reaction rates from temperatures and densities or require ACE data.

Acceptance: zero-strength bins are never selected; frequency ratios follow
integrated emissivity; constant `(R,Z)` emissivity reproduces volume weighting;
sampled sources can be visualized and transported through existing geometry.

## Phase 4 — geometry restrictions and shaped plasma models

Add optional conditioning to cells/materials/regions with finite proposal bounds,
explicit attempt limits, and acceptance diagnostics. Geometry-dependent prepared
sources retain/check geometry lifetime and generation, like tally plans.
Repeated universes need an explicit definition of occurrence versus cell ID.

Add shaped tokamak profiles only after phase 3 is validated: specify the mapping
from flux-like coordinates to physical `(R,Z)`, its validity domain, emissivity
measure, and Jacobian. A user-supplied `(R,Z)` emissivity grid remains the general
fallback for imported plasma calculations. Add non-axisymmetric toroidal
dependence through an additional angle axis or source-component mixture.

Keep generic geometry-conforming surface sampling separate from the analytic
surface samplers. CSG intersections, trimmed patches, transformations, and
normal orientation require a dedicated area/preparation strategy; they should
not be presented as an automatically supported consequence of a surface ID.

Acceptance: conditional distributions have documented normalization; empty or
low-acceptance regions give actionable failures; sampling on geometry boundaries
does not introduce lost particles or cross thin cells by arbitrary offsets.

## Phase 5 — correlations and external phase-space sources

Implement conditional tables so a sampled spatial bin can choose an energy
spectrum or angular law. Add joint distributions or a dedicated sampler where
independent components cannot express the required correlations.

Introduce DT/DD source presets first as clearly labeled monoenergetic models.
Thermal broadening, plasma flow, and energy-angle correlations are separate
physics implementations requiring an explicit model, validity range, and
reference-based verification. Do not label an arbitrary Gaussian a complete
fusion source model.

Add finite particle banks with specified replay versus probabilistic sampling,
weights, exhaustion behavior, and original history grouping. Grouped emission
events must preserve the definition of one statistical history. Keep the C
callback for specialized models; Python callbacks inside the per-history loop
are not required for built-in parity.

## WASM and GPU integration

After the native schema and initial samplers are stable, expose source preparation
and preview in WASM with the same field names and validation. Return TypedArrays
and run large sampling/transport batches in a Web Worker. Initially copy results
into JavaScript-owned buffers, matching Python's ownership guarantees.

For future GPU work, preserve immutable prepared tables, explicit RNG addresses,
bounded sampling work, and batched arrays. Declare backend support per source
type: geometry-conditioned rejection and arbitrary callbacks can remain CPU-only
until an actual device implementation exists. Source sampling is once per primary
history, so measure it separately from the collision/navigation cost.

## Verification and completion criteria

- Deterministic replay, nonoverlapping batch equivalence, and reordered-history
  sampling for every built-in source. Keep the current uint32 history-ID limit
  until a separate coordinated widening of transport and RNG identities.
- Distribution checks against known CDFs/moments using fixed seeds and justified
  statistical tolerances; never rely only on samples lying inside their bounds.
- Independent quadrature/reference checks for toroidal and shaped-profile weights.
- Units, normalization, zero-probability bins, invalid input, and rejection errors.
- Rotated/transformed sources and surface starts aimed inward/outward.
- C/Python/Lua parity, preview/result lifetime, NumPy dtype/shape, and source reuse.
- End-to-end void leakage and material transport checks, including mixed neutron
  and photon primary components and unchanged tally history normalization.
- Preparation time, samples/s, rejection acceptance, and memory measured on point,
  cylinder, weighted mesh, and tokamak sources before optimizing.

Update GNU/MSVC builds and installed public headers, add focused source unit tests,
and extend `docs/FIXED_SOURCE_WORKFLOW.md` at each phase. Existing parser failures
or missing test dependencies must be reported separately from source validation.

## Recommended first implementation

Complete phase 1, then cylinder/line/sphere spatial sampling and discrete/tabulated
energy distributions from phase 2. Next deliver the `(R,Z)` emissivity source from
phase 3. This provides useful fusion studies early while keeping one consistent
source API. Surface variants and the remaining independent distributions can be
added through that same API without delaying the first tokamak workflow.
