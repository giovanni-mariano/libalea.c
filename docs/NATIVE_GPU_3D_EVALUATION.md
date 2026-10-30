<!--
SPDX-FileCopyrightText: 2026 Giovanni MARIANO
SPDX-License-Identifier: MPL-2.0
-->

# Bounded native Metal 3D evaluation

Measured 2026-10-30 on Apple M1 after `92d8a77`. This completes the fresh bounded
three-view observation requested by M0; it does not close M0, establish CPU parity
or validate physics navigation. The private model, cameras, images, probes and
raw logs remain ignored/outside Git. Only aggregates are recorded here.

## Method and limits

Use the established local exterior, cutaway and component-close-up cameras.
Normal timing runs are 16x12, AA1, solid rendering, four CPU threads, compact Metal
instructions enabled, default spatial/selective-lazy/refill settings and profiling
disabled. Each mode/view runs in its own process with one loaded model and one
Metal context. Loading is about 9.5 seconds; snapshot preparation is 9.2-9.5 seconds.
Lazy render preparation/upload is included in first-frame latency.

A local copy of the native host bridge records command-buffer start/end timestamps
and allocated buffer lengths, using the unchanged production metallib. Valid
camera/configuration calls use the same native renderer entry point. No production
source or kernel is modified for these observations. GPU timings describe command
durations, not individual shader operations or hardware stalls.

Each process has a 90-second wall budget. Request up to three frames, but stop
before another frame when elapsed time plus the preceding frame estimate exceeds
85 seconds. Leave time for a CPU comparison and cleanup. All six normal runs
complete without hitting the supervisor timeout; process walls are 60.8-85.7 s.
An additional process measures five warm CPU renders per view on one reused model.
The CPU hierarchy is already warm after the first exterior render. These samples
are engineering observations, not release p95/worst-case bounds.

## Complete-frame timings and submission behavior

| View | Ordinary first s | Ordinary warm samples s | Progressive first s | Progressive warm samples s | CPU warm median s, five calls |
| --- | ---: | --- | ---: | --- | ---: |
| Exterior | 35.796 | None within budget | 27.379 | 27.123 | 0.441 |
| Cutaway | 23.224 | 23.219 | 20.219 | 20.018, 20.334 | 0.698 |
| Close-up | 13.643 | 13.634, 25.070 | 15.756 | 15.346, 15.726 | 0.0312 |

CPU warm ranges are 0.433-0.453 s exterior, 0.599-0.967 s cutaway and
0.0309-0.0325 s close-up. The first exterior CPU render takes 5.473 s, including
first-use preparation. Comparing that cold call with warm Metal would hide much
of the practical performance gap. CPU and Metal accuracy differences below also
mean these measurements do not establish interchangeable outputs.

Ordinary frames submit one long complete-ray command, except the last close-up
frame, which submits two. Maximum command durations are 35.752 s exterior,
23.210 s cutaway and 13.661 s close-up. In the close-up run, the existing adaptive
policy reduces its initial 512-pixel batch limit to 256 and then 128 after slow
frames. The third 192-pixel image therefore needs two submissions and takes
25.070 s. This illustrates a scheduling tradeoff; reducing batch size does not
ensure lower complete-frame latency. No controlled causal timing study of that
policy was performed here.

| Progressive view | Submissions per frame | Largest observed command ms | First completed pixel s | Largest observed callback gap ms |
| --- | --- | ---: | --- | ---: |
| Exterior | 1,620-1,676 | 26.8 | 3.386-3.420 | 43.7 |
| Cutaway | 1,223-1,297 | 30.8 | 0.589-0.629 | 31.2 |
| Close-up | 925-1,011 | 31.1 | 3.312-3.552 | 59.8 |

First completed pixel means any pixel, not a useful preview image. Small
submissions permit frequent callbacks while final images remain slow. Progressive
summed GPU durations account for roughly 94-97% of complete render-call wall time;
submission/readback overhead alone does not explain the result. This is not a
measurement of memory stalls, divergence, occupancy or continuation traffic.

Observed distinct context buffer lengths are 430,845,932 bytes (410.89 MiB) after
ordinary rendering and 460,828,652 bytes (439.48 MiB) after progressive rendering.
The six process high-water marks range from 2.616 to 3.284 GiB. Both modes fit the
provisional 1 GiB context, 4 GiB process and 15 s snapshot-preparation guards.
Resident and buffer figures refer to shared system memory and must not be added.
The renderer does not meet the declared preview/full-frame performance targets.
No larger private resolution was attempted.

## Correctness observations

All three ordinary/progressive Metal final buffers are **byte-identical** across
cell/material IDs, depth, normals and RGB. This establishes consistency between
those two paths on these samples, not CPU or analytic accuracy.

| View, 192 pixels | CPU/Metal ID or material differences | Maximum depth difference with matching hit IDs | Raw normal differences above 0.01 | Remaining after orienting CPU normals toward camera |
| --- | ---: | ---: | ---: | ---: |
| Exterior | 2 | 0.0007324 | 2 | 0 |
| Cutaway | 1 | 5.89746 | 36 | 2 |
| Close-up | 0 | 0.0005341 | 161 | 0 |

Source inspection confirms that Metal flips normals toward the camera before both
shading and auxiliary output. CPU shading performs the same orientation locally,
but its auxiliary buffer retains the original normal. Applying the CPU shading
orientation when comparing buffers removes all opposite-normal cases in exterior
and close-up. RGB differences on matching hit IDs are below 1e-6 in those views.
This is an auxiliary-output convention mismatch, not evidence that either output
is a validated outward physical-boundary normal. The contract needs an explicit
compatibility fix and public tests.

The cutaway differences remain unresolved. Besides the one ID/material mismatch,
two matching-ID pixels have depth differences above 0.01 (0.14697 and 5.89746).
The differing-ID pixel has a depth difference of 31.2417. Two matching-ID pixels
also retain normal/RGB differences; the largest RGB component difference is
0.51094. These cannot all be dismissed as normal orientation or near-equal depths.
No decision has yet established whether the CPU path, Metal arithmetic, clipping,
spatial specialization or traversal precedence is responsible. Preserve the local
pixel reproductions and diagnose them before claiming parity or optimizing them
away. The two exterior ID mismatches have near-equal depths (0 and 0.0001221
absolute differences); their cause is also not certified by this observation.

## Separate software attribution at 4x3

Run one instrumented ordinary frame per view at 4x3 to keep attribution bounded.
These are different rays from the 16x12 timing grid. Do not extrapolate their
counts or container fractions to every full-size image. Instrumented calls take
11.508, 8.402 and 4.628 seconds respectively; those are profiling costs, not normal
performance comparisons.

| Work, 12 rays | Exterior | Cutaway | Close-up |
| --- | ---: | ---: | ---: |
| Visited packed cells | 25,649 | 18,909 | 21,353 |
| Cell interval calls | 98,633 | 70,565 | 93,419 |
| CSG instructions | 14,540,303 | 10,898,462 | 5,811,815 |
| Primitive intersection evaluations | 2,672,665 | 1,981,525 | 1,113,849 |
| Interval steps | 338,332 | 200,871 | 148,827 |
| Container CSG instructions | 12,663,840 | 9,559,715 | 2,800,090 |
| Container share of instructions | 87.09% | 87.72% | 48.18% |
| Largest single-cell instruction share | 0.82% | 0.71% | 0.63% |
| Top 100 cells' instruction share | 17.07% | 18.67% | 25.53% |

Per-ray and per-cell totals agree for calls, instructions, primitive evaluations
and interval steps. No reported counter is saturated. Work remains distributed;
there is no single dominant cell. Container predicates deserve attention for
exterior/cutaway, while close-up also has substantial terminal work. Container
work includes all hierarchy levels; it is not the root-depth statistic from the
point-query report. These are work fractions rather than elapsed-time fractions.
The tested M1 public counter API still exposes timestamps only.

## External cancellation and reuse

A separate 16x12 exterior progressive run first warms the context (26.558 s).
Another host thread requests interruption 100 ms into the next render. The call
returns `ALEA_ERR_INTERRUPTED` about 5.5 ms after the request, before completing
any pixels. Clearing interruption and reusing that context completes all 192
pixels in 26.505 s, with a largest command of 28.1 ms. The reused output exactly
matches the completed ordinary image. The entire process takes 72.7 s.

This is one timed cancellation observation, not a p95 guarantee. It supports the
existing continuation/lifetime design while leaving useful-preview latency and
ordinary-render cancellation limitations open.

## Decision and next work

Keep Metal experimental and the optimized CPU renderer as the practical baseline.
Prioritize the reproducible correctness differences: (1) align auxiliary normal
semantics with the CPU contract while retaining shading behavior; (2) isolate
cutaway hit/depth differences and near-coincident exterior ownership. Use public
independently constructed regressions, never extracted private geometry.

Then evaluate ray-specific work reduction using the measured container/terminal
cost distribution. Retain complete-frame, first-pixel, submission, cancellation
and memory checks. The point-packet prototype's pruning counts do not justify a
ray index or a new caching policy. Final M0 accuracy classification and latency
acceptance remain open; CUDA execution still awaits NVIDIA hardware.

## Auxiliary-normal compatibility fix and reference investigation

Follow-up to `8fc7868`: both ordinary and resumable Metal shading now preserve
an independent copy of the source-oriented auxiliary normal. Only the shading
copy is oriented toward the camera. This changes auxiliary normal output while
retaining the existing RGB calculation, hit IDs and depth. The public API comment
now states this visualization convention explicitly; it is not an outward physical
boundary orientation guarantee.

An independently constructed public regression uses equivalent plane halfspaces
with opposite coefficients, a translated occurrence and a clipping face. It
checks analytic depths and CPU reference IDs/materials/normals/RGB in both render
paths and both instruction formats. Before the fix, the first case returns -Y
instead of the CPU reference's +Y. The full double-precision reference helper now
preserves source orientation separately from shading; an older parent-hole test
also checks the exposed sphere gradient (+Y) rather than its facing shading copy.

The translated test uses an independently expanded world-space CPU reference.
Rendering the filled source directly with the optimized CPU renderer misses its
material child after an undefined prefix. This is a public counterexample to using
that optimized path as the only geometry oracle, consistent with limitations
already noted by the existing full-segment test helper. The corresponding CPU
walk needs a dedicated regression/fix, including transformed boundary provenance,
before using it to decide all private hit discrepancies.

Private point checks, made recursively in FP64 on both original and quantized
rays, confirm material at samples near two earlier Metal hits that the optimized
CPU first-visible query skips. A third case is sensitive to ray quantization at a
sample just past the Metal entry. These samples support investigating the CPU
reference; they do not prove the exact first boundary or complete intervening
interval. An exhaustive full-segment reference exceeded the 90-second process
budget and provides no completed accuracy result.

The production cutaway image retains its saved ID/depth discrepancies after the
normal fix. A local isolated-ray diagnostic gives different results from the full
entry point even after checking matching geometry-buffer and parameter
fingerprints; it is not a validated substitute for the full-frame reproduction.
A diagnostic shader returning ray components also changes results, so no conclusion
about a particular compiler defect is drawn. A local variant disabling implicit
floating-point contraction does not resolve the full-frame discrepancy and is not
retained. Production compiler settings and traversal remain unchanged.

A bounded full cutaway run also completes with [Apple's Metal API/shader validation](https://developer.apple.com/documentation/xcode/validating-your-apps-metal-shader-usage)
enabled, with no invalid access reported; its hit discrepancies remain. This is a
validation observation, not a proof of correct numerical or interval behavior.
All raw captures, fingerprints, parameter dumps and probe variants remain ignored.

The next correctness dependency is a reliable CPU first-visible/interval oracle,
starting with the independently constructed undefined-fill prefix. Continue the
private cutaway classification afterward; do not waive the discrepancy or alter
Metal to imitate an unverified reference.


Acceptance for the normal correction: all three private 16x12 views were rerun
in ordinary and progressive modes with compact storage. Each mode pair is
byte-identical after the change, and IDs, depth and RGB are byte-identical to
its saved pre-fix image. Matching-hit normal differences above 0.01 fall from
2 to 0 exterior, 161 to 0 close-up and 36 to 2 cutaway. The two remaining cutaway
normal differences accompany the unresolved hit/depth investigation; they are
not simple orientation flips. Public checks pass: 74 Metal compute tests in each
instruction format, 5 builder tests, 6 resource tests, 16 CPU compute tests and
the independent analytic 2D/three-view 3D benchmark. Viewer event-loop,
progressive publication, cancellation and reuse checks pass. Both viewers are
rebuilt. No speculative floating-point or ray-traversal change is retained.
