<!--
SPDX-FileCopyrightText: 2026 Giovanni MARIANO
SPDX-License-Identifier: MPL-2.0
-->

# Ordered point queries through prepared spatial candidates

Follow-up to `77eabe0`. The public single-point resolver previously scanned
all cells to decide whether any fill/lattice existed, then scanned the root
universe again for flat models. Preparation could not remove that per-query
linear cost, even when the first cell contained the point.

## Selection and acceleration

The resolver now uses the existing recursive selection rule for all models:
select the first containing cell in the current universe's definition order,
then descend into its fill or lattice. An unresolved fill returns its container
as before. No model-wide hierarchy flag or new invalidation state is needed.
The existing hit metadata and transform/undefined-fill handling remain in the
same resolver. Public headers now describe the selection rule.

At each universe, a prepared spatial index supplies an ordered, unhinted
lookup. It uses the existing exact containment predicate. No previous-hit cache,
bounding-volume preference, or crossed-surface hint can change the winner.
Unavailable acceleration or failed candidate collection retains the ordered
linear scan. A universe index alone remains sufficient for the fallback;
`alea_prepare_query_acceleration()` supplies the full spatial index.

An early exact test of the first universe cell preserves cheap early winners.
Dense queries with more than 64 candidates and more than half the universe's
cells fall back to the ordered scan before sorting. Other large candidate sets
use comparison sorting above 64 entries; small sets retain insertion sorting.
This avoids quadratic sorting as candidate counts grow. Dense queries still
pay for candidate collection before falling back.

## Bounding-box scope

The evaluator rejects a point outside its stored root-node box before evaluating
CSG. For ordinary flat cells, the prepared spatial box widens those same stored
bounds, so it can conservatively filter candidates. Invalid root boxes cannot
pass the exact predicate.

The compact spatial representation clamps bounds to the finite float range.
A reproduced box around coordinate 1e39 can therefore pass the flat evaluator
while being absent from compact spatial candidates. Ordered point lookup now
uses its linear fallback for non-finite coordinates or coordinates beyond
FLT_MAX. The shared box representation remains unchanged: changing it to
infinities broke a pre-existing lattice region-transform regression. This fix
preserves point-query results through fallback rather than expanding every
spatial consumer's numerical contract.

## Validation

Five new public synthetic regression tests cover:

- Model order differing from cell-ID order and a previous-hit cache primed in
  the later overlapping cell, plus points on and one ULP beside boundaries.
- Positive/negative large coordinates, an unbounded halfspace, and shared faces.
- A late hit among 1,024 disjoint spheres, with fewer than sixteen spatial
  candidates and primitive evaluations, and universe-index-only fallback.
- Adding/removing a fill, undefined-fill selection, independent clones, region
  edits, and moving a cell out of and back into universe zero.
- A late winner among 128 overlapping shell boxes, exercising dense fallback.

Flat answers are compared to an independent ordered loop using the original
containment predicate. Hierarchical/lattice semantics also retain the existing
integration and spatial tests. The final native unit/integration/Lua suite and
serial unit/integration suites pass, including the region-transform regression
that rejected the initial shared-bound conversion experiment.

All test fixtures and measurements are independently synthetic. No private
models were needed. Diagnostic sources, raw timing output, and executables
remain ignored under build/flat-point/. This change supplies no new rendering,
transport-throughput, or GPU performance claim.

## Public API timing and limits

Release comparison against the library at 77eabe0 uses Apple Clang on macOS
ARM64 with -O3 -march=native, sequential processes, and three alternating
baseline/candidate process pairs. Preparation is outside query timing. Each
case repeats one fixed point, with 1,000 warm-up calls and 100,000 measured
calls (1,000 measured calls for the dense shell case). Every answer is checked.
Medians below are nanoseconds per query across the three processes.

| Synthetic model | Query | Baseline | Candidate |
| --- | --- | ---: | ---: |
| 256 separated spheres | First cell | 149.9 | 49.3 |
| 256 separated spheres | Last cell | 1,242.5 | 117.6 |
| 256 separated spheres | Miss | 1,229.2 | 33.4 |
| 4,096 separated spheres | First cell | 4,520.6 | 28.9 |
| 4,096 separated spheres | Last cell | 22,424.1 | 87.5 |
| 4,096 separated spheres | Miss | 22,491.5 | 20.1 |
| 4,096 dense shell/terminal cells | First cell | 4,479.0 | 54.0 |
| 4,096 dense shell/terminal cells | Last cell | 98,064.0 | 107,929.0 |
| 4,096 dense shell/terminal cells | Miss outside all boxes | 23,858.0 | 21.0 |

Separated spheres have radius one and centers spaced four units apart. Dense
cells share a shell between radii one and ten, followed by a terminal inner
sphere. The center query rejects all shells before finding the last cell.
That dense late-winner case is about 10% slower: candidate collection adds
work before the linear fallback. Worst-case point queries remain linear or
worse for candidate sorting; the fix makes no universal logarithmic-time claim.
The initial dense sorting path measured about 247 microseconds before fallback
was added; the final result is about 108 microseconds.

Repeated fixed points favor cache locality. Sub-microsecond differences across
model sizes are sensitive to timing noise and code/data placement. These results
establish removal of the global scan and a strong benefit on spatially selective
queries; they are not random-point, whole-application, or imported-model results.
No new persistent model buffers or cached hierarchy flags are introduced.
