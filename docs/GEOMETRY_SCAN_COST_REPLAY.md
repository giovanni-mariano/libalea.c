<!--
SPDX-FileCopyrightText: 2026 Giovanni MARIANO
SPDX-License-Identifier: MPL-2.0
-->

# Indexed scan cost: captured-context replay

Measured on 2026-10-29 at documentation checkpoint `3b27f4a`, with production
code still at `9375ed9`. This investigation follows the
[plane-group experiment](GEOMETRY_PLANE_GROUP_EXPERIMENT.md). Production source
is unchanged. Diagnostic source, executables, and raw results remain under
ignored `build/scan-cost/`.

The next bounded implementation experiment should call the existing plane
kernel directly from the indexed scan, using the original payload and copying
it only when a normal is needed. This reduces replay time by 16.0–16.4% across
three samples without a prepared table. It has not yet been measured in an
optimized complete frame and is not a retained production change.

## Method

A diagnostic copy captures complete indexed-scan contexts during actual nearest
recovery in a broad private 512x512 slice. Capture occurs after the existing
context-reuse cache, before the scan clips its upper bound to an existing hit.
Each context preserves the local ray, original ordered card list, input window,
breakpoint threshold, occurrence transform, and initial nearest-hit state.

Each of three sequential processes uses a different deterministic seed for a
uniform reservoir of 256 scan calls. Every process observes 256,517 calls and
815,382,047 references. No calls are excluded by the capture eligibility guards.
The three samples contain 820,083, 764,192, and 803,898 references respectively,
covering 84, 85, and 80 distinct lists. Sampling is uniform over calls, not over
references; each selected call is replayed with its entire list.

Planes account for 61.0–61.6% of sampled references and quadrics for 27.2–29.7%,
compared with 60.67% and 27.76% in the earlier full-frame attribution. This is
useful workload coverage, not a guarantee that every expensive context appears.

All variants process identical contexts, in original card order, with the same
existing scalar kernels, entry-window clipping, torus rejection, normal work,
and nearest-hit insertion. They differ as follows:

| Variant | Work changed |
| --- | --- |
| Indexed | Original indexed scan copied into a non-inlined replay helper |
| Resolved, copied | Pre-resolved records bypass surface/primitive lookup and payload validation; payloads still copy into the union before generic dispatch |
| Resolved, typed | Same records, with direct typed kernel calls on original payload pointers; copies remain for normal calculation and fallback types |
| Packed, typed | Same typed path with aligned payload copies packed by list |
| Indexed, direct plane | Original index and payload lookup path; only planes bypass the union copy and generic kernel dispatch, copying on an in-window root for normal calculation |

The three resolved variants move validation and preparation outside the timed
loop. Their differences therefore do not isolate memory latency, copy cost,
or dispatch cost individually. The direct-plane variant retains the original
lookup checks and does not use the prepared records.

The benchmark-only records and packed copies together occupy 21.1–22.6 MiB,
under a 64 MiB cap. Preparing them takes 5.5–6.3 ms for these sampled lists;
that time is excluded from replay timing. This is not a whole-model preparation
or cache-lifecycle measurement. Captured geometry and resolved payloads remain
in process memory and are not exported as fixtures.

## Results

Apple Clang 21, macOS ARM64, `-O3 -march=native -DNDEBUG`, one configured worker.
Each process performs five passes with rotating variant order, approximately
16 million references per variant per pass. Values below are the within-process
median nanoseconds per reference; reductions compare each variant with that
process's indexed median.

| Variant | Sample 1 | Sample 2 | Sample 3 | Reduction range |
| --- | ---: | ---: | ---: | ---: |
| Indexed | 22.929 | 22.206 | 20.859 | — |
| Resolved, copied | 13.111 | 12.921 | 12.490 | 40.1–42.8% |
| Resolved, typed | 10.783 | 10.485 | 10.197 | 51.1–53.0% |
| Packed, typed | 9.237 | 9.154 | 9.217 | 55.8–59.7% |
| Indexed, direct plane | 19.163 | 18.653 | 17.495 | 16.0–16.4% |

These results support investigating work surrounding the primitive kernels.
They do not show that the remaining time is pure intersection arithmetic or
that any single lookup accounts for a particular percentage. Code generation,
validation, copying, dispatch, and locality change together in the resolved
variants. Their percentage reductions cannot be added.

Repeated replay also warms data and omits the intervening work of a complete
frame. All variants run in a process containing the diagnostic prepared data.
The original replay helper's layout and inlining differ from its production
caller. Consequently these are hypothesis-selection results, not frame-time
predictions, hardware-counter attribution, or evidence for a GPU/SIMD backend.

## Validation and limits

Every sampled context compares all four alternatives with the original scan:
failure state, nearest-hit count, tested-surface count, distance, surface and
primitive identifiers, and all normal components agree. Independently, every
sampled primitive is solved through generic, typed original-payload, and typed
packed-payload paths, including roots outside the active window. Root counts
and returned root arrays agree bit-for-bit: 1,390,945 returned roots across the
three private samples. Timed work-count checksums also agree within each run;
those checksums are not substitutes for the explicit result comparisons.

The three instrumented private frames and four public control frames match
the established cell/material/density/flag digests. Public repeating and
transformed fixtures capture 62 and 3 contexts; the deep-fill and ellipse views
capture none. Those latter controls validate capture-build frame behavior,
not optimized replay coverage. The complete frames still execute the original
scan; only subsequent replay executes the alternatives.

These checks establish compatibility for the samples, not arbitrary numerical
or malformed-input behavior. Native/serial suites were not rerun because no
production source changed. A production candidate still needs focused boundary,
failure, transform, and alias tests, ordinary all-hit compatibility, and full
integration checks before acceptance.

## Next experiment

Start with the direct plane call, preserving the existing payload accessor,
scalar root formula, surface order, inclusive window, tie handling, and failure
semantics. Keep the change local to nearest indexed scans; add no persistent
metadata or new ownership/invalidation obligations. Do not tighten the upper
bound within the loop or combine this with the rejected plane-group filter.

Measure complete first/warm frames against retained production in three
alternating process pairs at one and four workers, along with public controls,
sampled rows, detail views, cancellation, and independent numerical regressions.
Retain only a clear net complete-frame benefit with matching outputs and no
repeatable public-control regression above 5%. Keep all slower runs in the
analysis. If the replay saving does not survive, remove the candidate before
considering more elaborate preparation or a change in recovery breadth.

The local implementation plan is `plans/geometry-direct-plane-scans.md`.
The one-second objective remains a warm complete 512x512 broad 2D slice after
loading/preparation, including tracing and rasterization but excluding display
and hashing. It remains unmet. Private inputs, paths, coordinates, raw hashes,
and diagnostics stay ignored; only aggregate observations are tracked here.
This report is not added to the public documentation-site allowlist.
