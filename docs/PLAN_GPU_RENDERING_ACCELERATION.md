<!--
SPDX-FileCopyrightText: 2026 Giovanni MARIANO

SPDX-License-Identifier: MPL-2.0
-->

# GPU rendering acceleration: baseline and next gates

Status: 17 October 2026. The current renderer and supported controls are described
in [GPU compute](GPU_COMPUTE.md). Experimental results and retirement decisions
are recorded in [GPU boundary analysis](GPU_BOUNDARY_ANALYSIS.md).

## Retained implementation

- Shared universe BVHs and procedural rectangular lattices.
- Automatic conservative bounds and adaptive CSG program simplification.
- Selective deferred parent clipping, preserving interval intersections and
  parent-created boundary normals. Expensive child universes retain eager clipping.
- Slot refilling and bounded progressive rendering with cancellation.
- Optional SAH/weighted BVH rebuilds, pending broader performance evidence.

The manual root-universe grid was removed because automatic spatial
specialization covers its purpose without user-selected bounds. Shared candidate
masks were removed after isolated, alternating five-frame measurements made both
ITER cutaway views slower. Surface patches and uniform boundary compilation were
also retired. Their evidence remains available; these are not pending features.

## Accuracy gate before further acceleration

The pre-fix 32x24 ITER comparison had 14 CPU/GPU cell-ID disagreements and 11
matching-ID depth differences greater than 0.1 model units. The CPU viewer is not
an independent truth oracle: bounded full interval queries confirm that it
misses earlier material in several inspected rays. One near-clipping-plane case
crosses a sub-millimetric same-material interval; a different case assigned
water instead of steel because near-singular quadrics were translated to distant
origins before FP32 packing. The packing fix passes an isolated regression that
fails with the old implementation. The complete frame now has 10 ID differences
and 11 matching-ID depth differences above 0.1. Bounded interval checks now
classify all 21 disputed pixels: six ID disagreements and all eleven large depth
differences involve later CPU viewer entries; four ID disagreements concern
same-material clipped intervals only 0.000386–0.000466 units wide. The largest
local Metal entry error in the 17 earlier-material cases is 0.021598 units.
These are sampled local checks, not a whole-ray or whole-model proof.

Before introducing a root cache:

1. Preserve the completed disputed-ray classification and the four unresolved
   FP32 ownership limitations. Extend independent checks to new views and
   representative ordinary image sizes; do not repeat CPU-viewer comparisons
   as if they were ground truth.
2. Reduce confirmed renderer defects to small maintained fixtures with analytic
   or independent flattened references. Distinguish FP32 input loss, arithmetic
   error, traversal defects, overlaps and CPU viewer defects.
3. Fix correctness defects before comparing acceleration. Report unresolved FP32
   limitations explicitly; do not use agreement between Metal modes as proof.

The local diagnostic programs and full ITER captures are ignored build artifacts,
not another maintained tracing subsystem.

## Corrected timing baseline

At `a3e66ec`, five warm default progressive frames per 16x12 ITER view measured
19.511 s median for the cutaway and 21.571 s for the opposite cutaway, with
bitwise-stable repeated output. Scene preparation took about ten seconds after
about eight seconds of loading. This establishes a current small-image baseline;
ordinary-resolution utilization and broader-view independent accuracy checks
remain separate work. See the findings report for ranges and methodology.

## Conditional root-cache experiment

Historical four-ray simulations found reuse, but mostly for inexpensive planes.
Only after the accuracy gate, compare a small cache of expensive primitive roots
against no cache. Include complete local-ray identity, transformed instances and
progressive continuation state in the design. Use actual frame time, state memory
and preparation cost; cache hit rate alone is insufficient.

Measure multiple views and ordinary image sizes, with a warm-up and alternating
repeated samples. Retain the change only with repeatable improvement and unchanged
reference results. Keep CUDA/HIP, GPU transport and new meshing architectures
outside this baseline.
