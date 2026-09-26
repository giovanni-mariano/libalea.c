<!--
SPDX-FileCopyrightText: 2026 Giovanni MARIANO

SPDX-License-Identifier: MPL-2.0
-->

# Geometry validation

ALEA provides complementary validation paths. Choose one according to the
claim the result must support.

| Method | Best use | Interpretation |
|---|---|---|
| `alea_find_overlaps()` | Fast initial screen | Bounded heuristic; a clean result is not proof |
| `alea_validate_geometry()` | Broad transport-style diagnostics | Sampled rays with structured crossing findings |
| `alea_validate_geometry_slice()` | Diagnose analytical slice boundaries | Sampled curves with slice provenance |
| Slice error query | Certify a bounded supported 2D domain | Verified regions and boundaries, with explicit unresolved output |

## Sampled transport-style validation

The occurrence-aware validator reports gaps, overlaps, non-adjacent
transitions, missing neighbors, ambiguous boundaries, and incomplete samples:

```c
#include <alea_geo_validator.h>

alea_geom_validator_options_t options;
alea_geom_validator_options_init(&options);
options.ray_count = 20000;
options.seed = 12345;
options.max_errors = 1000;

alea_geom_validator_result_t result;
alea_geom_validator_result_init(&result);

if (alea_validate_geometry(sys, &options, &result) != 0) {
    fprintf(stderr, "validation failed: %s\n", alea_error());
}

/* Read findings with alea_geom_validator_error_get(). */
alea_geom_validator_result_free(&result);
```

Set `ALEA_GEOM_VALIDATE_DOMAIN_BOUNDS` and `validation_bounds` when unowned
space inside a closed world box must count as an interior gap. Before treating
an empty finding list as meaningful, inspect `truncated`, `incomplete_rays`,
and `incomplete_slice_samples`.

## Verified slice error queries

A slice error query partitions a required slice rectangle into pages. Each
page can contain:

- verified gap or overlap regions;
- verified straight boundary intervals or circular arcs;
- locally confirmed point witnesses;
- contextual diagnostic findings; and
- rectangles that remain unresolved.

The receipt separates successful execution from a complete geometrical claim.
`scope_classified` is true only when every face and adjacent physical boundary
in the page core was classified. `output_complete` says that every finding
discovered by the executed scan was retained. Both must be true on every page,
and no defective region or boundary may be present, before the covered domain
can be described as clean under the supported proof policy.

Pages are independent and may be run individually with
`alea_slice_error_query_run_page()` or as a scratch-bounded batch with
`alea_slice_error_query_run_pages()`. A failed page can be retried without
discarding earlier completed pages. Do not mutate the borrowed system during a
query; changing its geometry generation invalidates remaining pages.

## Current certification scope

Certification currently requires a coordinate-aligned slice. Supported proofs
include combinations of world-axis and projected plane boundaries, selected
sphere sections, rectangular-lattice seams, transformed fills, and hierarchy
occurrences. More complex curved arrangements, hexagonal lattices, ambiguous
coincidences, and numerically inseparable boundaries may be returned as
unresolved.

Unresolved output is intentional: it prevents an unsupported area from being
silently interpreted as valid geometry. Consult
[`alea_geo_validator.h`](https://github.com/giovanni-mariano/libalea.c/blob/main/include/alea_geo_validator.h)
for the exact structures, limits, and evidence semantics.
