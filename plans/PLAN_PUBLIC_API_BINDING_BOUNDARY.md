# Plan: Close the Public API Boundary for Native Bindings

**Status:** Proposed
**Primary consumer:** AleaTHOR CPython binding
**Related plans:** `plans/BINDINGS.md`, `docs/PLAN_PUBLIC_API_OPAQUE_SYSTEM.md`

## Summary

Make the installed libalea headers the complete and exclusive C contract used
by native language bindings. AleaTHOR must compile with `-Iinclude` and without
any path below `src/` while retaining its current raycast, logging, model-load,
and geometry functionality.

The migration is additive first: introduce the missing public operations,
migrate AleaTHOR, enforce the boundary in CI, and only then remove leaked
storage types and unusable public macros.

## Problem

The public API is mostly sufficient for AleaTHOR, but the binding still relies
on selected implementation details:

1. It includes `raycast/raycast.h` to read raycast result vectors directly.
2. It forward-declares the private `alea_log_set_callback()` contract.
3. It transfers loaded systems by modifying `mcnp_model_t` and
   `openmc_model_t` fields directly.
4. libalea publicly declares storage details such as `alea_surface_entry_t`
   and the unusable `alea_surface_at` macro even though `alea_system_t` is
   intended to be opaque.

This coupling makes bindings sensitive to internal layout changes and prevents
the installed headers from being the authoritative API definition.

## Current AleaTHOR Usage

The binding in `/home/giovanni/projects/aleathor_1208/aleathor` currently:

- uses `alea_halfspace()` four times for public geometry construction;
- uses `alea_surface_get()` 23 times for positive/negative node inspection;
- does not use `alea_surface_at`;
- does not name `alea_surface_entry_t`;
- includes `raycast/raycast.h` once;
- reads `result->segments`, `result->hits`, and their vector counts directly;
- forward-declares and calls `alea_log_set_callback()`;
- reads `model->sys` and writes `model->owns_sys` when taking ownership from
  MCNP and OpenMC model wrappers.

Removing `alea_surface_at` therefore does not require an AleaTHOR source
change. The raycast, logging, and model-ownership dependencies must be resolved
before the binding can compile against public headers only.

## Goals

1. Build AleaTHOR using only headers installed from `include/`.
2. Preserve all current Python-visible behavior and result fields.
3. Keep `alea_system_t` opaque to external consumers.
4. Provide efficient low-level access without exposing internal vectors.
5. Make ownership transfers explicit and safe.
6. Add checks that prevent future internal-header dependencies.
7. Preserve source compatibility during the additive migration phases.

## Non-Goals

1. Do not redesign raycasting algorithms or result storage.
2. Do not move the AleaTHOR binding into libalea as part of this work.
3. Do not redesign the Python API.
4. Do not promise ABI stability for headers below `src/`.
5. Do not introduce a shared-library target solely for this migration.
6. Do not make all public structs opaque in one patch.

## Design Principles

- Installed headers define the supported external contract.
- Opaque handles may expose efficient accessors and borrowed array views.
- Public accessors validate bounds and document pointer lifetimes.
- Internal modules and white-box tests may retain direct storage access through
  explicitly internal headers.
- Add replacement APIs before removing leaked declarations.
- Binding requirements should be demonstrated by real consumers, not inferred
  from internal structs.

## Phase 1: Complete the Public Raycast Result API

### Required API additions

Add to `include/alea_raycast.h`:

```c
int alea_raycast_segment_resolution_flags(
    const alea_raycast_result_t* result,
    size_t index,
    uint8_t* out_flags);

size_t alea_raycast_hit_count(
    const alea_raycast_result_t* result);

int alea_raycast_hit_get(
    const alea_raycast_result_t* result,
    size_t index,
    double* out_t,
    int* out_surface_id);
```

### Semantics

- Accessors return `-1` for a null result, invalid index, or required null
  output pointer.
- `alea_raycast_hit_count(NULL)` returns `0`, matching other count accessors.
- Accessors return value copies and never expose internal vector pointers.
- Existing `alea_raycast_segment_get()` remains unchanged.
- A raycast mode that does not promise a full hit list may return zero hits;
  its existing contract remains authoritative.

### Implementation

Implement thin accessors next to the existing public raycast result accessors.
Do not duplicate tracing logic or allocate memory.

### Tests

Add tests for:

- segment flags agreeing with the internal segment produced by a known trace;
- valid hit count and ordered hit values;
- null result handling;
- out-of-range segment and hit indexes;
- null output pointers;
- result reuse after a second raycast.

### Acceptance criteria

- AleaTHOR can produce its existing `segments` and `hits` Python objects
  without dereferencing `alea_raycast_result_t`.
- No binding source includes `raycast/raycast.h`.

## Phase 2: Publish the Logging Callback Contract

### Required API additions

Publish the callback type and setter in an installed header, either directly in
`alea.h` or in a new `include/alea_log.h` included by `alea.h`:

```c
typedef void (*alea_log_callback_t)(
    alea_log_level_t level,
    const char* file,
    int line,
    const char* message,
    void* user_data);

void alea_log_set_callback(
    alea_log_callback_t callback,
    void* user_data);
```

Prefer a dedicated public `alea_log.h` if additional logging configuration is
expected to become public. Keep internal logging macros in
`src/util/alea_log.h`.

### Required documentation

Document that:

- the callback and `user_data` are borrowed and retained until replaced;
- `alea_log_set_callback(NULL, NULL)` unregisters the callback;
- callbacks may be invoked from worker threads;
- the consumer owns synchronization and callback-state lifetime;
- language bindings must acquire their runtime lock, such as Python's GIL,
  before calling runtime APIs;
- bindings must unregister callbacks before interpreter shutdown or module
  teardown;
- recursive logging from inside the callback is unsupported unless explicitly
  made safe by the implementation.

### Tests

Add tests for registration, message delivery, user data, replacement,
unregistration, and log-level filtering.

### Acceptance criteria

- AleaTHOR deletes its local callback typedef and `extern` declaration.
- The logging bridge compiles using installed headers only.

## Phase 3: Add Explicit Model System Access and Ownership Transfer

### Required API additions

Add equivalent borrowed-access and ownership-transfer functions for MCNP and
OpenMC models:

```c
alea_system_t* mcnp_model_system(mcnp_model_t* model);
alea_system_t* mcnp_model_take_system(mcnp_model_t* model);

alea_system_t* openmc_model_system(openmc_model_t* model);
alea_system_t* openmc_model_take_system(openmc_model_t* model);
```

`*_model_system()` returns a borrowed pointer. `*_model_take_system()` returns
the pointer and leaves the model detached so destroying the model cannot
destroy the transferred system.

### Transfer semantics

Conceptually:

```c
alea_system_t* mcnp_model_take_system(mcnp_model_t* model) {
    if (!model) return NULL;
    alea_system_t* sys = model->sys;
    model->sys = NULL;
    model->owns_sys = 0;
    return sys;
}
```

The exact internal representation may differ after models become opaque.

### Tests

Add tests proving:

- borrowed access does not alter ownership;
- taking returns the original system;
- destroying a detached model does not destroy the system;
- destroying the transferred system remains valid;
- taking twice returns `NULL` on the second call;
- null model arguments return `NULL`;
- behavior is identical for loaded and wrapped models.

### Acceptance criteria

- AleaTHOR no longer reads `model->sys` or writes `model->owns_sys`.
- The binding uses `*_model_take_system()` before destroying the model shell.

## Phase 4: Migrate and Isolate the AleaTHOR Binding

Update the AleaTHOR binding to:

1. Replace direct segment flag reads with
   `alea_raycast_segment_resolution_flags()`.
2. Replace direct hit vector reads with `alea_raycast_hit_count()` and
   `alea_raycast_hit_get()`.
3. Remove `#include "raycast/raycast.h"`.
4. Use the installed logging callback declaration.
5. Replace model field mutation with `mcnp_model_take_system()` and
   `openmc_model_take_system()`.
6. Remove libalea `src/` directories from the extension include path.

Build the extension against the vendored static archives if desired, but give
the compiler only public include roots for binding source compilation.

### AleaTHOR compatibility tests

Run or add tests covering:

- primitive and surface construction;
- rich raycast segments including `resolution_flags`;
- global ray hit output;
- MCNP and OpenMC load from file and string;
- ownership and repeated destruction paths;
- Python logging delivery and interpreter teardown;
- SIGINT interruption;
- debug and release extension builds.

### Acceptance criteria

- No binding source includes any path below libalea `src/`.
- No binding source forward-declares an undeclared libalea symbol.
- No binding source dereferences an opaque libalea result.
- AleaTHOR's existing test suite passes.

## Phase 5: Finish the Opaque Surface Boundary

This phase follows `docs/PLAN_PUBLIC_API_OPAQUE_SYSTEM.md`.

### Public construction

- Remove `alea_surface_at` from `include/alea.h`.
- Use `alea_halfspace(sys, surface, -1)` for the negative halfspace.
- Use `alea_halfspace(sys, surface, +1)` for the positive halfspace.
- Optionally add `alea_surface_inside()` and `alea_surface_outside()` only if
  they materially improve examples and bindings.

### Public inspection

Expand `alea_surface_info_t` into a read-only value view and add:

```c
int alea_surface_info(
    const alea_system_t* sys,
    size_t index,
    alea_surface_info_t* out);
```

Keep `alea_surface_get()` as a compatibility wrapper until its removal is
separately decided.

### Internal storage

- Move `alea_surface_entry_t` out of `include/alea_types.h`.
- Define it in an internal surface header.
- If useful, provide an explicitly internal helper:

  ```c
  #define alea_surface_entry_at(sys, idx) (&(sys)->surfaces.data[(idx)])
  ```

- Convert public-style tests to `alea_halfspace()` or public inspection.
- Keep direct entry access only in tests that intentionally validate storage or
  acceleration internals.

### Acceptance criteria

- `rg "alea_surface_at" include examples` returns no matches.
- `rg "alea_surface_entry_t" include` returns no matches.
- Public construction and inspection examples compile with only `-Iinclude`.
- AleaTHOR requires no change because it already uses `alea_halfspace()` and
  `alea_surface_get()`.

## Phase 6: Consider Opaque Format Models

After all consumers use model accessors, decide whether `mcnp_model_t` and
`openmc_model_t` should become opaque.

If made opaque, first provide public operations for every intentional use:

- system borrow/take;
- destruction;
- export configuration get/set;
- MCNP cell parameter inspection and controlled mutation;
- inline transform inspection;
- non-owning model wrapping.

Do not make the models opaque until all supported consumers have replacements.
This phase is optional and should be a separate compatibility change.

## Boundary Enforcement

### Public consumer compile test

Retain the header smoke test, and add a real consumer program that:

- creates and destroys a system;
- creates a surface and selects both halfspaces;
- inspects surface information;
- performs a raycast and reads segments, flags, and hits;
- registers and unregisters a logging callback;
- loads or wraps a format model and transfers its system.

Compile it with only:

```sh
-Iinclude
```

The test must not include `src/` or rely on the repository root include path.

### Static source checks

Reject internal includes in binding directories:

```text
#include "core/..."
#include "util/..."
#include "raycast/..."
#include "slice/..."
#include "mcnp/..."
#include "openmc/..."
```

Also audit for:

- manual `extern` declarations of libalea functions;
- direct fields on opaque result types;
- copied definitions of private structs;
- build scripts adding libalea `src/` to binding include paths.

Static checks complement compilation; they do not replace it.

## Compatibility Strategy

Phases 1 through 4 are additive and should not break existing C consumers.
Phase 5 intentionally removes declarations that were public in location but
not usable through the documented opaque contract.

Compatibility rules:

- do not change `alea_raycast_segment_get()` incompatibly;
- retain `alea_surface_get()` while introducing structured inspection;
- add ownership-transfer operations before hiding model fields;
- do not renumber public enums or change existing result ownership;
- document borrowed pointer lifetimes explicitly;
- keep static archive internals outside the public contract even though their
  symbols may be visible to link tools.

## Suggested Patch Sequence

1. Add raycast segment-flag and hit accessors with tests.
2. Publish logging callbacks with lifecycle documentation and tests.
3. Add MCNP/OpenMC system borrow/take APIs with ownership tests.
4. Migrate AleaTHOR and remove its internal include path.
5. Add the complete public-consumer compile test and CI source audit.
6. Remove `alea_surface_at` and introduce an internal surface-entry helper.
7. Add/expand `alea_surface_info()` and move `alea_surface_entry_t` internal.
8. Evaluate opaque MCNP/OpenMC models as a separate change.

Keep these patches separate enough that API additions, binding migration, and
intentional compatibility removals can be reviewed independently.

## Completion Criteria

The boundary work is complete when:

1. AleaTHOR builds and tests with only libalea installed headers.
2. No supported binding includes a libalea internal header.
3. No supported binding forward-declares private libalea functions.
4. No supported binding reads opaque result storage directly.
5. Model ownership transfer uses named public operations.
6. `alea_system_t` remains opaque and no installed macro dereferences it.
7. Surface storage entries are internal-only.
8. Public consumer compilation and source-boundary checks run in CI.
9. Full libalea unit, integration, Lua, and relevant AleaTHOR tests pass.

## Expected Final Binding Pattern

```c
#include <alea.h>
#include <alea_mcnp.h>
#include <alea_openmc.h>
#include <alea_raycast.h>

mcnp_model_t* model = mcnp_load("input.i");
alea_system_t* sys = mcnp_model_take_system(model);
mcnp_model_destroy(model);

alea_raycast_result_t* result = alea_raycast_result_create();
alea_raycast(sys, ox, oy, oz, dx, dy, dz, t_max, result);

for (size_t i = 0; i < alea_raycast_segment_count(result); ++i) {
    uint8_t flags;
    alea_raycast_segment_resolution_flags(result, i, &flags);
}

for (size_t i = 0; i < alea_raycast_hit_count(result); ++i) {
    double t;
    int surface_id;
    alea_raycast_hit_get(result, i, &t, &surface_id);
}

alea_raycast_result_destroy(result);
alea_destroy(sys);
```

No internal header, vector layout, ownership flag, or private symbol is needed
by the binding.
