<!--
SPDX-FileCopyrightText: 2026 Giovanni MARIANO
SPDX-License-Identifier: MPL-2.0
-->

# Shared geometry query experiment

Initial implementation and navigation follow-up, 28 October 2026. Branch: `experiment/geometry-queries`.
Baseline: `8c585eb48e1cdba6f26e02345ab4b5df1c21a798`.

## Finding and implementation

The existing selected-interval walker already serves persistent particle
navigation, slice rows, first-visible rendering, and ordered segments. It uses
prepared caches on the source system and separate traversal scratch. The source
must remain alive and unmodified while navigation is active. An additional
prepared-geometry object or evaluator was unnecessary for this experiment.

Compact batch queries previously enriched each crossing with coincident local
boundary provenance even when no surface output was requested. The implementation
now selects that work from `ALEA_RAY_BATCH_SURFACES`, in both the common worker
arena and full-path/projected-owner batch paths. It reuses the walker's existing
enrichment control, also used by fast navigation. Ownership resolution, midpoint
verification, path capture, and resource budgets retain their existing policies.

Surface-requesting batches continue to enrich crossings. Other scalar consumers,
mesh sampling, particle boundary actions, and first-visible rendering retain
their previous enrichment behavior. Public interfaces are unchanged. The change
adds no geometry buffers or per-ray caches.

The experiment deliberately targets unnecessary query work. It does not add
hardware backends or establish transport accuracy for the approximate Metal
renderer. Selected-owner tracks also retain their existing semantics: complete
coverage queries are needed to detect hidden overlaps.

## Consumers and regression checks

[The characteristic-track example](../examples/c/characteristic_tracks.c)
uses bounded public compact batches with full occurrence paths. Three adjacent
unit-width regions deliberately share a material but remain distinct cells.
Forward and reverse tracks start a quarter unit inside the model and traverse
2.5 units. Independently expected lengths are 0.75, 1, and 0.75.

The example checks complete ownership coverage, compares each segment's cell,
length, and occurrence against persistent strict navigation, and evaluates
uncollided attenuation with prescribed removal coefficients 0.2, 0.4, and 0.8
per model length. Both directions give `exp(-1.15)`, approximately
0.316636769379. It requires no external nuclear data and is not a scattering or
MOC solver. Region-specific coefficients illustrate an external physics mapping.

[New regression tests](../tests/unit/test_geometry_queries.c) preserve an
analytically specified 0.00001-wide region between same-material neighbors and
verify occurrence identities under track reversal. Repeating-lattice tests
compare ownership, lengths, densities, resolution flags, and path/projection
identities with and without requested surface output, exercising both batch
allocation paths. Existing tests continue to cover boundary actions, transformed
fills, gaps, overlaps, cancellation, budgets, first-visible output, and raster
pixel ownership.

Run the focused regression tests and analytic example with:

```sh
make -j4 test-geometry-queries
```

## Measurements

Release baseline and experimental libraries were built separately from original
and modified sources with `RELEASE=1`: `-O3 -march=native -DNDEBUG`, using the
native TinyPar backend on macOS ARM64. Benchmark clients used the same release
optimization flags. `ALEA_NUM_THREADS=1` fixed the worker count.

[The layered benchmark](../benchmarks/geometry_queries.c) traces 512 rays through
32 adjacent boxes, checking each of 16,384 returned segments against analytic
cell IDs and endpoints after every query. Each process performs one warm-up and
seven measured queries. Four process pairs alternate baseline/experimental order,
giving 28 measured queries per configuration. Times include synchronous output
publication and exclude geometry construction and preparation.

| Output | Baseline median | Experimental median | Reduction |
| --- | --- | --- | --- |
| Material segments | 7.723 ms | 7.271 ms | 5.9% |
| Material segments with full paths | 8.275 ms | 7.931 ms | 4.2% |

Material-only ranges were 7.593–11.563 ms and 7.119–11.853 ms respectively;
full-path ranges were 8.130–8.559 ms and 7.660–8.279 ms. Scheduling noise remains
visible. These measurements support a modest improvement on this workload.

The existing sphere slice benchmark was also run at 1024 by 1024, with 20
iterations per process and four alternating process pairs. Medians of the four
reported process means were 1.031 versus 1.003 ms for compact tracing, and 1.806
versus 1.757 ms for fused tracing and rasterization. Its initial allocations and
warm-up are included, so this is an end-to-end compatibility screen rather than
a precise speedup attribution. Each configuration returned 2,664 segments and
used 17,825,792 bytes of raster field storage.

Build and run the layered benchmark with:

```sh
make RELEASE=1 BUILD_DIR=build/query-release BIN_DIR=build/query-release-bin \
    build/query-release-bin/geometry_queries
ALEA_NUM_THREADS=1 build/query-release-bin/geometry_queries 0
ALEA_NUM_THREADS=1 build/query-release-bin/geometry_queries 1
```

For a baseline comparison, build the recorded base commit in a separate checkout
with the same flags and link this benchmark source against its archive. Use
separate build/output directories when changing release or threading settings;
Make does not infer object invalidation from changed compiler flags alone.

## Validation completed

- Native unit and integration suites, including numerical geometry, navigation,
  transport, slices, CPU rendering, and the new geometry-query regression cases.
- Lua suite and public-header checks through the standard test targets.
- Unit and integration suites in an isolated `TINYPAR_BACKEND=serial` build.
- Analytic characteristic example against both native and serial libraries.
- Release benchmark endpoint checks and the existing sphere raster benchmark.

Commands used for the full checks:

```sh
make -j4 test-unit test-integration test-lua
make -j4 BUILD_DIR=build/geometry-query-serial \
    BIN_DIR=build/geometry-query-serial-bin TINYPAR_BACKEND=serial \
    test-unit test-integration
```

GPU kernels and approximate GPU precision were not changed or validated by this
experiment. The representative-model pass below extends the CPU evidence;
ordinary-resolution GPU utilization and accelerated particle navigation remain
separate work.

## Representative model validation

[The model benchmark](../benchmarks/geometry_model_queries.c) loads MCNP inputs,
separates loading and preparation from queries, and records process peak RSS.
It exercises ownership-only and surface-requesting batches, with and without
full hierarchy paths, changing-direction strict navigation, and ray-row
rasterization. Output digests include endpoints, IDs, density, resolution flags,
full paths, and requested boundary IDs. These are comparisons against the
baseline, not independent proofs of imported geometry.

[The comparison driver](../tools/check_geometry_models.py) alternates baseline
and candidate processes. For each public model it uses three process pairs,
1,024 alternating forward/reverse rays, one warm-up and five measured queries
per output configuration, and a 1,024 by 1,024 slice. Navigation uses up to 16
material sources and 32 direction changes per source. The tests use isolated
release builds, the same compiler settings, and one worker, as above.

| Public model | Ownership segments baseline / candidate | Full-path segments baseline / candidate | Raster baseline / candidate |
| --- | --- | --- | --- |
| Repeating lattice | 13.435 / 12.603 ms | 15.470 / 14.661 ms | 13.431 / 12.565 ms |
| Eight translated fills | 3.506 / 3.280 ms | 4.362 / 4.062 ms | 3.797 / 3.586 ms |
| Transformed nested lattice | 16.701 / 15.071 ms | 18.724 / 16.947 ms | 16.550 / 14.748 ms |

Medians are taken over 15 measured queries per configuration. Ownership-only
queries improve by approximately 5–10%, and raster timings by 6–11%. Requested
surface output remains on the original enrichment path. Its timing variation
does not establish a performance change. Navigation is also unchanged; its
sub-millisecond measurements are a compatibility screen, not a speedup claim.

All compared public output digests match, including the surface-requesting
configurations. Process peak RSS is approximately 17–26 MiB across the public
workloads, with no evidence of added query storage. Allocator and process-level
RSS variation should not be attributed to the optimization.

The eight-level model is a new analytic fixture, independently authored from
nested spheres, translations, and a plane. Eight translations of 0.25 place the
terminal interface at world x=2. A regression checks forward/reverse path
lengths, nine-entry occurrence paths, surface identity, strict navigation to
the interface, and reversal into the incident region. It passes in native and
serial builds.

### Navigation boundary repair

The targeted 16-row transformed-lattice case originally stopped two source
histories after direction changes: baseline and the selective-enrichment change
both produced 478 events and two unresolved-ownership failures. The three main
public workloads completed 512 navigation events each.

The cause was a forward location sample crossing a nearby surface before the
walker selected the interval owner. The walker then assigned the short interval
before that surface to the owner after it. Strict point verification correctly
rejected the mismatch. A bounded interior resampling path already existed for
first-visible queries; navigation now uses that path when its location sample
reaches or exceeds the next crossing. Entry and exit distances are preserved.
No displacement epsilon or interval suppression is added. Failure to establish
an interior sample still fails closed.

Two additional regressions cover the public cylinder reproducer and an analytic
translated plane. The latter preserves a flight of exactly `2^-49` model units
before world x=2, in both directions, after restart and direction changes, in
all three validation modes. It checks incident/outgoing cells, surface ID,
positive crossing distance, and subsequent collision position. Both new tests
fail against the original baseline and pass with the repair.

Native unit, integration, and Lua suites and serial unit/integration suites pass.
Three alternating release process pairs now complete all 512 targeted navigation
events without ambiguity or failure. The main workloads retain matching event
digests; all compared track, requested-surface, and raster digests also match.
The targeted navigation digest intentionally changes as previously interrupted
histories finish. This is evidence for the reproduced cases, not a general
numerical guarantee.

Reproduce the targeted workload using the model benchmark with:

```sh
bin/geometry_model_queries tests/data/mcnp_nested_lattice_transformed.mcnp \
    16 1 3 9 -0.8 0.8 0
```

The comparison driver requires explicit `--expect-navigation-repair` to accept
navigation-only differences from a failing baseline. Every pair must preserve
non-navigation output and complete the same sources with at least as many events
and no ambiguity or failures. Other differences and incomplete runs fail the
gate. Exact output agreement, baseline navigation cleanliness, candidate
cleanliness, and accepted repair are recorded separately. A zero process exit
status also requires all expected warm-up and measured records with output
digests; missing or duplicate measurement records fail the gate.

### Private large model

ITER inputs remain outside the repository. Raw logs, query digests, and any
model-derived diagnostic output stay under Git-ignored `build/`. No model,
extracted geometry, rendered image, or model-derived fixture is added to the
publishable files. The optional driver argument takes a user-supplied path;
only aggregate observations belong in this report.

A one-row center probe completed in both builds with identical recorded
ownership, full-path, navigation, and raster digests. Its measured
queries were about 24–25 ms; one process pair and four navigation events do
not establish a speedup. Loading took about seven seconds and preparation about
five seconds. Peak process RSS was approximately 2.6–2.9 GiB and is not evidence
of a memory improvement.

A wider eight-row pilot took about 33 seconds per ownership-only batch. The
large variation between sampled rays prevents extrapolating from the center
probe. A denser 32-row focused pass hit the 180-second process timeout in both
baseline and candidate before the first warm-up batch completed. The timeout
includes loading and preparation. Those runs contain no completed track,
navigation, or raster comparison, so they establish neither agreement nor
disagreement and provide no speedup evidence. Large-model acceptance requires
completed bounded workloads, including expensive rays, before any
ordinary-resolution performance claim.

Build the benchmark in an isolated release directory and compare two binaries:

```sh
make RELEASE=1 BUILD_DIR=build/model-release BIN_DIR=build/model-release-bin \
    build/model-release-bin/geometry_model_queries
python3 tools/check_geometry_models.py \
    --baseline /path/to/baseline/geometry_model_queries \
    --candidate build/model-release-bin/geometry_model_queries \
    --expect-navigation-repair
```

To include a privately supplied ITER model, add `--private-model /local/model`
and `--private-view XMIN XMAX YMIN YMAX Z` using local coordinates. The optional
`--private-rows N` defaults to 32. Row count and timeout are explicit resource
limits; private slice coordinates are not embedded in the driver.
The driver rejects a nonignored output directory when a private model is
supplied. Keep output under ignored local storage. The public comparison
requires `--expect-navigation-repair` when comparing the repaired navigator
against the recorded failing baseline.

## Bounded individual-ray profiling

A follow-up on the repaired branch separates ray preparation from query timing
in the existing `large_model_probe` tool and prints its existing ownership,
interval, enrichment, snapshot, and allocation-growth counters. Three public
fixtures retain identical cold/warm segment and work counts compared with the
previous probe linked against the same library. The probe's scalar
`--hier-fast-segments` path retains boundary enrichment; these are diagnostic
cost measurements, not ownership-only batch timings.

Four selected private rays each completed two evaluations within a 75-second
process limit, including loading and preparation. Processes ran sequentially
with the same optimized libraries and one worker. Across this small sample,
second-evaluation times range from 0.023 to 22.606 seconds. Two selected rays
with no recovery take about 23–32 ms; two with three recovery snapshots take
about 20.6–22.6 seconds. This is a targeted sample, not a latency distribution
for the model. Preparation is excluded from these query times.

A separate diagnostic build times only `raycast_recovery_sample` on one
expensive ray, with a 60-second process limit. It completes the query in about
21.733 seconds, of which three recovery calls consume about 21.704 seconds
(over 99%). Those calls perform approximately 932 million primitive tests in
total. The ordinary walker records fewer than 20,000 tests for that query:
recovery uses a separate scratch result, so ordinary counters omit that work.
The diagnostic timers and progress logging remain local ignored artifacts.
They are used for attribution; instrumented time is not a speedup measurement.

Code inspection identifies the exceptional nearest-crossing search as the
next acceleration candidate. `raycast_recovery_sample` enumerates validation
breakpoints through filled universes; `raycast_universe_surfaces` visits every
primitive leaf of every cell's Boolean tree, including repeated references.
An existing `raycast_cell_indexed_surface_hits` helper can use prepared cell
surface lists with a tree fallback when no index is available. This finding
selected the nearest-only substitution measured in the experiment below.

Before implementation, the experiment required a public, independently authored
fixture with repeated Boolean references and analytic crossings. Preserve the
nearest crossing distance, surface identity, occurrence frame, and coincident
boundary policy, including transformed fills and lattices. Do not discard
recovery or replace it with a larger displacement. Compare the indexed and
original paths in alternating release runs: require at least a twofold reduction
on completed recovery-heavy workloads, no added persistent geometry storage,
and investigate any repeatable regression above 5% on recovery-free controls.
Public correctness checks must pass regardless of timing. Keep the general
complete-coverage enumerator's semantics outside this narrowly scoped change.

All private inputs, per-ray parameters, raw counters, and diagnostic source
copies remain under ignored local storage or outside the repository. Only
aggregate observations are recorded here.

## Indexed nearest-crossing experiment

Nearest-only recovery now calls the existing cell-indexed surface enumerator
when it visits cells in a filled universe. Prepared lists contain each surface
card once in first tree-occurrence order. Distinct cards that share a primitive
remain distinct, preserving the original first-hit tie policy. Cells without
an index use the existing tree fallback. The general full-breakpoint enumerator
retains its original tree path. No new geometry storage or public API is added.

The independent regression extends the public translated-plane fixture with
coincident cards, a distinct plane 0.00001 units away, and 4,096 repeated
references to a shared Boolean subexpression. It checks analytic nearest
crossings from both directions and between the nearby planes, exact card and
primitive identity, and world normals against the tree fallback. It also checks
that indexed enumeration performs at least ten times fewer primitive tests on
this fixture. The old implementation fails that work-reduction assertion.
A second regression compares indexed and fallback nearest results in both
directions for repeating lattices, transformed nested lattices, and deep fills.

Native unit/integration and Lua suites and serial unit/integration suites pass;
the two added regressions pass in both builds. Three alternating release process
pairs against the immediate pre-change baseline retain all public workload
track, surface, navigation, and raster digests, including the repaired 512-event
navigation workload. These comparisons isolate this optimization from the
previous boundary repair.

A separate bounded private comparison uses full hierarchy paths and requested
surface output. The baseline is the immediately preceding repaired branch;
only the nearest-only universe enumeration differs. Reproducing that baseline
requires retaining the original `raycast_tree_primitives` call in
`raycast_universe_surfaces`, then rebuilding in a separate release directory.
Both clients use identical `-O3 -march=native -DNDEBUG` settings and one worker.
Three process pairs alternate baseline/candidate order, with one warm-up per
process. Each heavy process has one measured query; each control process has
30. Loading and preparation are outside query timing; processes run sequentially
under a 90-second limit.

| Targeted private workload | Baseline median (range) | Indexed median (range) |
| --- | --- | --- |
| Recovery-heavy query, 3 measured samples each | 20.922 s (20.358–20.987) | 7.267 s (7.247–7.637) |
| Recovery-free control, 90 measured samples each | 23.618 ms (23.006–25.526) | 23.809 ms (23.058–26.672) |

Every warm-up and measured result has matching endpoint, cell/material,
density, resolution, full-path, occurrence, and requested surface-ID digests
across baseline and candidate. The heavy median improves by 2.88 times (65.3%
less time), exceeding the predeclared twofold gate. The control median increases
by about 0.8%, below the 5% regression threshold. No new geometry buffers or
persistent indexes are introduced. These targeted results support retaining the
change; a query still takes seconds, and ordinary-resolution large-model plotting
or rendering has not been established. Private inputs, clients, parameters,
and raw comparison output remain local and ignored.

## Integration decision and follow-up

Retain the small query-controlled enrichment change for branch review: it
reduces unnecessary work, adds no per-query storage, passes the compatibility
checks, and shows modest release improvements in the measured clients. No
minimum numerical speedup threshold was fixed before this exploratory timing;
these results must not be presented as satisfying a predeclared percentage gate.

The public representative-model pass supports retaining the selective-enrichment
optimization for review. The follow-up repairs the reproduced navigation
boundary failure with an independent short-flight regression. Completed bounded
ray probes identify exceptional ownership recovery as the dominant measured
cost. The indexed nearest-crossing experiment meets its scoped acceptance
criteria and is retained for review. Further acceleration should begin with
profiling the remaining recovery cost; broader large-model acceptance and
ordinary-resolution rendering remain unresolved.

The common-mechanism and deterministic-consumer milestones use existing
infrastructure. Benchmark results remain local evidence, not proof of a
universal acceleration architecture or general solver readiness.

Future experiments should be independently scoped from measured costs. A second
hardware backend, root cache, Boolean dependency engine, or mesh representation
requires its own hypothesis, resource budget, and acceptance criteria. Preserve
the independent GPU accuracy gate in the existing GPU plan.

The [nearest-boundary recovery follow-up](GEOMETRY_RECOVERY_SEARCH.md) profiles
the repaired, indexed baseline and evaluates sharing surface solves across cells
within each universe occurrence. Its measurements and acceptance decision are
recorded separately from this experiment.
