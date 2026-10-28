<!--
SPDX-FileCopyrightText: 2026 Giovanni MARIANO
SPDX-License-Identifier: MPL-2.0
-->

# Prepared universe surfaces for nearest recovery

Follow-up to [slice and row profiling](GEOMETRY_SLICE_VALIDATION.md), using
`b859b99` as the baseline on `experiment/geometry-recovery-search`.

## Change and contracts

The occurrence-scoped bitset already avoided repeated primitive solves, but
nearest recovery still scanned each cell's surface references on every query.
Preparation now optionally stores each universe's distinct surface cards in
first-occurrence order. Recovery traverses that list with the current local ray
and occurrence transform. Distinct cards sharing a primitive remain distinct;
list order preserves coincident-hit tie behavior. Transformed occurrences and
lattice elements share card identities, never solved intersections.

Only nearest-boundary recovery consumes the list. Other enumeration retains
its existing behavior. Missing cell indexes, invalid references, insufficient
memory, and list allocation failure retain the previous per-cell scan, including
its bounded occurrence bitset and tree fallback. Single-cell universes use the
existing cell index.

Retained list storage is capped at 64 MiB per model. Preparation uses at most
1 MiB for its temporary seen-card bitset, plus the list under construction;
allocator overhead and transient reallocations are additional. Failed list
shrinking discards the optional list. A separate cache-ready bit publishes the
lists after construction under the existing preparation lock. Universe and
cell-index invalidation discard them; clones prepare their own storage.

Recovery checks interruption every 1,024 indexed references, at universe-cell
boundaries, and during tree traversal. Cancellation remains cooperative and
whole-application latency can include other operations.

## Validation method

Independent tests compare prepared recovery with tree fallback for close planes,
coincident aliases, translated fills, and repeated lattices. The fallback tests
explicitly disable both preparation levels. A lifecycle test checks cloning,
direct cell-index rebuilding, geometry edits, and destruction order.

Complete public slices compare all pixel fields at 512 and 1,024 resolution in
split and fused modes, over three alternating baseline/candidate pairs with a
cold and warm frame per process. Private validation uses the same 32 stratified
rows at each resolution as the preceding report, with exact completed-output
and work-counter comparisons. Repeated timing uses three alternating process
pairs at 1,024 rows, one worker, and the same cancellation budgets. These are
systematic samples, not a random estimate of a full frame.

All private paths, view coordinates, hashes, raw logs, diagnostic sources, and
model inputs remain ignored under `build/` or outside the repository. Only
aggregate measurements belong in this report.

## Measurements and decision

All 64 initial private row samples complete and match the preceding checkpoint
in geometry digests, segment/interval counts, snapshot counts, and ownership
query counters. All 192 completed rows in the three alternating timing pairs
also match exactly between builds.

For the 32 sampled rows at resolution 1,024, median total tracing time across
three processes falls from 17.148 s to 10.281 s: approximately 40% less time,
or 1.67 times the throughput for this workload. Individual baseline totals are
17.123–17.183 s and candidate totals are 10.218–11.002 s. Pooled row medians
are approximately 627 ms and 369 ms. This is a bounded recovery result, not a
whole-model or full-frame speed claim.

Preparation has a cost. Median reported preparation is 4.946 s at baseline and
5.678 s with the lists, approximately 0.732 s additional in these samples.
Candidate preparation ranges from 4.906 to 6.351 s, so the estimate is noisy.
Model loading is separate. Process peak RSS varies substantially across runs;
these measurements do not establish a reduction in whole-process memory.

A 100 ms cooperative cancellation request now returns at approximately
100–105 ms on the measured expensive rows. The one-second profiling budget
finishes in 1.004 s with seven completed and nine cancelled rows, no query
errors, and no digest for cancelled rows. Completed rows after cancelled ones
exercise cleanup and reuse. This is observed responsiveness, not a universal
hard deadline.

The local preparation diagnostic reports 187 prepared universes, containing
9,823,721 cell surface references and 604,854 retained ordered cards in aggregate.
The lists occupy 2,419,416 bytes (approximately 2.31 MiB), below the 64 MiB cap.
This counts each universe separately; it is not a count of globally unique
surfaces. It explains the benefit: subsequent recovery avoids those repeated
references while preserving the primitive solves required for that occurrence.

The final public slice comparisons match all fields in both modes at both
resolutions for repeating, deep translated, and transformed lattice models.
Warm total-time changes range from approximately -7% to +3%; no tested control
shows a regression above the 5% screening threshold. Small timings remain
sensitive to noise. The final native unit/integration/Lua suite and serial
unit/integration suite pass, including nine focused geometry tests and the
characteristic-track example.

The private material-containing complete slices also match all fields at both
resolutions: approximately 91.7% material pixels and no flagged pixels. These
single cold frames are completion checks, not evidence of repeatable frame
speedup. The final build includes the publication check for queries that
prepare only a subset of caches. Private timing comparisons of sampled rows
were taken before this additional read-side check; the final full-slice and
public comparisons exercise it.

The broader focused private view at resolution 512 still exceeds the 60-second
whole-process screen, including load and preparation, without completing the
trace. No full-field comparison or final peak-memory claim is made for that
cancelled process. The optimization is worth retaining for its repeatable row
benefit and small retained storage, but this view is still not interactive.
The next investigation should measure the remaining recovery intersection and
ownership work with the prepared list active before choosing another index or
parallel scheduling change.
