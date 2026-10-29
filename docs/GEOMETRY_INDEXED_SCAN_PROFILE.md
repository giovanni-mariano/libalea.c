<!--
SPDX-FileCopyrightText: 2026 Giovanni MARIANO
SPDX-License-Identifier: MPL-2.0
-->

# Indexed surface scan attribution after child-cell lists

Measurement-only follow-up at `9375ed9`, following the
[ordered child-cell implementation](GEOMETRY_RECOVERY_CHILD_CELLS.md).
Production code is unchanged. The objective is to identify a bounded next
experiment, rather than choose a hardware backend from a total timing share.

## What is scanned

Nearest recovery uses ordered surface-card lists prepared from all cells in a
universe. These lists remove repeated cards but are not ray-specific spatial
candidate sets. Reuse of identical ray/list/transform contexts already avoids
many complete scans. The remaining calls still visit every card in their list.

An ignored diagnostic copy attributes scans, visits, time, roots, and exact
plane-normal repetition to each list. It uses the existing release libraries,
Apple Clang 21, `-O3 -march=native -DNDEBUG`, native TinyPar, and one worker.
Four synthetic 512 frames and the complete broad private 512 frame match their
release cell/material/density/flag digests. Existing work counters are unchanged
in this attribution run. Its 28.476 s frame time is instrumented and must not
be treated as a new release throughput measurement.

Of 187 prepared lists, 111 are scanned during actual nearest recovery. The
256,516 scans visit 815,376,379 card references. Cost is distributed: the largest
contributor accounts for 7.49% of references and 8.04% of indexed diagnostic
time; the top five account for 20.14%/19.60%, and the top ten for 32.21%/31.25%.
Optimizing a single exceptional universe would leave most of the work intact.

| Primitive class | Reference visits | Share |
| --- | ---: | ---: |
| Planes | 494.68 million | 60.67% |
| General quadrics | 226.35 million | 27.76% |
| Axis-aligned cylinders | 76.56 million | 9.39% |
| Remaining types | 17.79 million | 2.18% |

Only 5,457 returned roots lie inside the scan-entry window, including 4,932
plane roots and 339 quadric roots. These are diagnostic candidate roots,
not distinct selected physical crossings. General quadrics return no roots
at all in about 97% of their solver calls. This low yield motivates candidate
reduction, but does not prove that an inexpensive safe rejection test exists.

Two plausible small changes have limited reach in this workload. Distinct
cards sharing a primitive account for only 0.81% of weighted references;
removing them also requires care with aliases. After a nearer hit is found
inside a scan, only 0.60% of references are visited with an upper bound larger
than that new hit. Tightening that bound alone cannot remove most scan work.

## Why terminal-cell bounds are insufficient

The current recovery contract enumerates referenced surfaces, including their
crossings outside the cell's Boolean region. A synthetic probe demonstrates
this without private geometry: a child region is the intersection of x >= 0
with a box spanning x = 10 to 20. Translate its occurrence by +3. From world
x = 1, a positive-x ray over t = 0 to 5 misses the child cell's box, yet recovery
reports the referenced plane at world x = 3, t = 2. The unplaced plane is behind
the ray. Applying the terminal-cell box as a filter would remove that result.

Therefore the existing cell spatial index cannot simply replace these scans
while preserving their contract. Primitive-surface bounds require their own
conservativeness audit; current plane bbox placeholders are not a suitable
finite geometric certificate. This does not rule out spatial acceleration
under a deliberately different, independently validated recovery contract.

## Exact-normal plane groups: observational rejection test

Many planes have exactly equal normal coefficients and different offsets.
Weighted over the measured scans, planes in repeated-normal families account
for 40.92% of all references. Large enough groups may amortize a shared test.
A second ignored diagnostic groups planes by exact numerical equality of their
normal coefficients, without approximate merging. It tests the extreme offsets
using the existing plane kernel but still executes every original surface scan.
Thus it observes hypothetical exclusions without changing results.

For groups with at least sixteen members, this shadow run has:

| Observation | Result |
| --- | ---: |
| Prepared eligible groups across the model | 2,481 |
| Group tests over the frame | 3,013,027 |
| Added extreme-offset kernel calls | 6,026,054 |
| Plane visits that would be excluded | 145,753,655 |
| Share of all observed reference visits | 17.87% |
| Share of observed plane visits | 29.46% |
| Excluded roots lying inside the active window | 0 |

A group is hypothetically rejected only when both endpoint calls report
parallelism, or both finite roots lie strictly on the same outside side of the
window. Equality remains eligible. Other cases retain the original scan.
The production proposal must explicitly guard finite inputs and uncertain
arithmetic and establish monotonicity of the actual kernel, not merely of the
ideal plane equation.

The shadow run matches all four synthetic frames and the complete private
frame. It visits 815,415,528 references through 256,527 scans: eleven fewer
context-cache hits than the attribution run produce eleven additional scans.
This small change is consistent with the bounded context cache hashing list
addresses, which can make reuse sensitive to process address/layout variation; no recovery, fill-recursion, ownership-attempt, or
in-window-root counts change. Use the shadow run's own denominator for its
17.87% figure. Its instrumented frame time, 30.176 s, is not a candidate speedup.

A separate deterministic synthetic probe checks 42,955 group/window cases,
including 22,955 windows on or next to existing roots, scales from 1e-200 to
1e200, near-parallel directions, and nonfinite/overflow fallbacks. It observes
457,970 in-window roots and zero false exclusions. This supports the hypothesis
but is not an exhaustive floating-point correctness proof.

Contiguous runs of sixteen equal-normal planes account for only about 1.36%
of visits. Restricting the experiment to such runs would miss most of the
potential. Conversely, groups as small as two would require roughly 45.5 million
group tests in this frame. Sixteen is an initial amortization threshold,
not a proven optimum. Group lookup and map traffic may erase a theoretical gain.

## Next experiment and limits

The saved next plan is optional exact-normal plane groups for nearest recovery,
initially with at least sixteen members, while retaining the original surface
order. Cap additional maps/descriptors at 8 MiB per model and group flags at
4 KiB per active scan helper. Missing metadata or uncertain arithmetic uses the
existing scan. Cache publication, invalidation, direct rebuilds, clones,
interruption, aliases, transformed occurrences, and equal-distance ties all
need explicit coverage before retaining an implementation.

The decision must come from paired release measurements after paying for group
tests and memory access. A 17.87% potential reference reduction is not a 17.87%
frame-time prediction. The current one-second goal concerns a prepared 512x512
broad 2D slice, excluding loading, preparation, and display; this investigation
does not establish that goal, 3D performance, transport coverage, or GPU benefit.

Raw diagnostics, private paths, coordinates, input data, and hashes remain
ignored under `build/`. Only aggregate private observations are tracked here.
No production implementation or public documentation-site publication occurs
in this measurement step.
