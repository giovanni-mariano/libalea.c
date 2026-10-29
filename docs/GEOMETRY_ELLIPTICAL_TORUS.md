<!--
SPDX-FileCopyrightText: 2026 Giovanni MARIANO
SPDX-License-Identifier: MPL-2.0
-->

# Elliptical torus intersections, normals, and slice curves

For the consolidated decision and current limitations, see the
[geometry query status](GEOMETRY_QUERY_STATUS.md). This report retains the
baseline and conclusions of its own experiment.

Correctness follow-up to `ba5a3ba` and the [torus recovery experiment](GEOMETRY_TORUS_RECOVERY.md).
All new fixtures are independently specified, with no private model geometry.

## Defect and exact transformation

The point evaluator represents a Z-axis torus by

    ((sqrt(x*x + y*y) - R)/C)^2 + (z/B)^2 = 1.

CPU ray intersections and normals previously used C in both directions,
ignoring the axial semiwidth B. Analytical slice curves had the same assumption.
For R=5, C=1, B=2, the axial line through radial coordinate 5 crosses at z=+2
and z=-2. The old ray solver instead describes crossings at z=+1 and z=-1.
Point ownership and crossing distances could therefore disagree.

Set z'=(C/B)z. In the transformed coordinates the equation is exactly

    (sqrt(x*x + y*y) - R)^2 + z'*z' = C*C.

Scale both the ray's axial origin and axial direction. Normalize the resulting
direction for the quartic calculation, then divide the shifted roots by that
direction's speed to restore the original ray parameter. Normalizing without
this final division would give incorrect distances. This is an exact coordinate
transformation for finite positive widths, not a geometric approximation.

## Implementation

- CPU non-circular intersections use the transformed ray and existing quartic
  solver. Circular intersections retain their previous arithmetic path.
- Non-circular normals come from the gradient of the original evaluator
  equation, including the axial width. Robust normalization handles very wide
  tubes whose gradient magnitude is below the old absolute cutoff.
- The constant polynomial coefficient uses its factored form on the non-circular
  path, avoiding subtraction of two terms of order R^4 in thin rings. A test
  with R=1e6, C=1, B=2 reproduced an axial-distance error before this refinement.
- Non-circular zero-major-radius cases reduce to a quadratic. Spindle cases
  reject the extraneous sheet introduced by squaring the equation. The ray
  path also respects the evaluator's non-positive axial-width convention by
  removing axial dependence; an axial ray then has no discrete crossing.
- Finite elliptical slice curves use the correct perpendicular-section radii.
  Central parallel sections use the general quartic representation. General
  scanlines scale the axial components and restore the original scanline
  parameter. Curve boxes retain the unscaled plane basis and true axial width.
- Nearest-recovery rejection uses the true axial semiwidth. The earlier envelope
  covering both the incorrect circular surface and the ellipse is no longer
  needed. Existing conservative guards and fallbacks remain.

The existing Metal intersection and normal routines already apply axial scaling;
they are unchanged. This work does not claim new GPU validation. Public APIs
are unchanged; the curve's extra axial-width field is internal.

## Independent validation

The first four regression tests fail against the previous implementation and
pass after correction. The complete new suite has twelve tests covering:

- All three torus axes, axial widths from 1e-6 to 1e9, and known axial roots.
- Four-root radial sections, tangent roots, and a nearby miss. Tangent-root
  positions use an absolute tolerance of 1e-7 because repeated roots are
  ill-conditioned; ordinary roots use tighter checks.
- Inner and outer normals against explicit gradients and finite differences,
  plus unit normals for axial width 1e15.
- 1,200 deterministic oblique rays at translated coordinates, checking every
  returned intersection against the point evaluator and rejection envelope.
- A new MCNP fixture with unique inside/outside materials, forward/reverse
  tracks, and analytical material lengths of 1, 4, and 1 along a six-unit ray.
- Strict navigation through a rotated, translated occurrence, with boundary
  distances 1 and 4 and the expected material transitions.
- Zero-major-radius, spindle, non-positive axial-width ray cases, thin rings,
  analytical slice intersections and boxes on all axes, and strongly scaled
  axial scanlines.
- Both split and fused 32-by-32 raster paths, checking every cell, material,
  density, and resolution flag against an independent pixel-center equation.

A separate diagnostic checks every pixel at 512 and 1,024 resolution in both
raster modes against that same independent equation at axial offset 1.5. All
fields match, including the region the previous circular interpretation missed.
A deterministic comparison of 12,288 circular rays gives bitwise-identical root
and normal digests before and after the change.

The existing thirteen focused geometry tests and characteristic-track example
also pass. Public repeating, translated-fill, and transformed-lattice complete
slices match the previous implementation in all fields at both resolutions and
in both modes. These controls contain no targeted elliptical-torus correction;
old/new equality is a compatibility check. Elliptical correctness is established
by analytical references, not by requiring equality to the old wrong results.

## Private checks and limits

All private inputs, coordinates, raw digests, diagnostic sources, logs, and
binaries remain under ignored `build/elliptical-torus/` or outside the repository.
Only aggregate observations are recorded here.

The bounded 32-row sets at each resolution complete without errors and retain
their earlier output digests. Those views were already observed to miss torus
surfaces in the recovery windows, so their equality does not establish the
elliptical correction. A separate local primitive audit constructs axial rays
from the imported non-circular ring-torus parameters and checks the known roots
in both directions. It is an intersection check, not a validation of full-model
ownership, overlaps, or the complete occurrence hierarchy.

Numerical limits of the shared polynomial solver, extremely ill-conditioned
parameters, and singular surface points remain broader library concerns. The
finite elliptical tests here do not establish exhaustive support for every
degenerate analytical slice configuration. This change makes the supported
finite elliptical geometry consistent across the tested CPU consumers.

The private primitive audit contains four axial-ray checks: the previous solver
fails two and the corrected solver passes all four. This is a small, targeted
check of the imported parameters, not an exhaustive private-model validation.
The final native unit/integration/Lua suite and serial unit/integration suite
pass. The corrected circular controls show no geometry-output changes.

The final-build private row recheck also completes all 64 rows with unchanged
output digests and no query errors. A final cooperative-cancellation check
returns at approximately 100–105 ms for 100 ms requests, with subsequent queries
remaining usable. These checks support integration compatibility in the sampled
views; the independent ellipse equations remain the correctness reference.
