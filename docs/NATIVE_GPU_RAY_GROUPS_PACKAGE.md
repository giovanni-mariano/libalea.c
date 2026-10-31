# Packaged Metal ray-group execution

Native ordered intersection-group evaluation is a supported option in
`libalea_compute` and the normal `render3d` tool. It applies to ordinary rendering
through `alea_compute_render()`.

```sh
make compute USE_METAL=1 RELEASE=1
make -C tools ../bin/render3d USE_METAL=1 USE_SDL=1 RELEASE=1
ALEA_METAL_RAY_GROUPS=1 bin/render3d model.i --backend metal \
    --width 640 --height 480 -o render.png
```

Set the option before creating a Metal context or starting the tool. Its default
is `0`; only `0` and `1` are accepted. Both wide/compact code and jump threading
are supported. No benchmark shader or driver is required.

| Setting | Operations |
| --- | --- |
| `ALEA_METAL_CSG_GROUPS=1` | Point classification and ordered selection |
| `ALEA_METAL_RAY_GROUPS=1` | Ordinary ray rendering |

The settings are independent. When both are enabled, one context-owned set of
maps and signed references serves both operations. The existing optional 64-MiB
preparation/storage budget is shared, not doubled. Unavailable optional storage,
allocation/device limits or no eligible groups retain generic Metal evaluation.
Malformed programs and interrupted preparation remain errors. Destruction and
later construction failures release all owned buffers.

The source model may be freed or edited after snapshot creation. Existing contexts
retain their geometry and captured settings; create a new context to incorporate
geometry edits. Independent contexts have independent data and pipelines.

## Rendering scope and diagnostics

The native loop preserves emitted leaf order, duplicate occurrences, half-space
sense, finite errors, short-circuiting, interval validity, boundary ties and
at-cursor identity. The existing primitive solvers, numerical precision, spatial
certificates, hierarchy/lattice traversal and normal construction remain in use.
Generic bytecode remains available for unrecognized programs and fallback.

`alea_compute_render_progressive()` continues to use the generic resumable
evaluator. It retains bounded per-submission work and its existing cancellation
behavior. Enabling ordinary ray groups therefore does not accelerate interactive
progressive frames. Independent tests check final ordinary/progressive equality.
Native bounded group resumptions remain a separate implementation milestone.

Ray and cell CSV profiles append `csg_groups` and `instruction_unit` columns.
Grouped ordinary rendering reports `1,bytecode_and_group_ops`; generic and
progressive rendering report `0,bytecode`. Existing columns retain their order.
Primitive and interval counts keep their meanings; grouped instruction counts
include bytecode, group dispatches and visited leaves and use a different unit.
The bounded root-reuse diagnostic also records grouped leaf evaluations.

The option is off by default because small mixed expressions can be slower, and
snapshot preparation/storage has a cost. The preceding
[controlled evaluation](NATIVE_GPU_RAY_EXECUTION.md) measured 18–20% improvements
on private 4x3 samples, a public short-case regression, and an incomplete exterior
throughput comparison. Packaging follows explicit user direction; it does not
assert universal speedup, interactive/full-resolution targets or physics accuracy.

## Integration validation

The packaged evaluator passes 406 analytic cursor/expression cases in both
formats with groups enabled/disabled, plus visited/skipped nonfinite errors and
multiple material-interval cases. Tests compile the actual production shader;
there is no second implementation in the benchmark.

Public API checks cover flat and mixed CSG, nested/transformed fills,
finite/repeating/outer lattices, unbounded halfspaces, clipping and AA. Contexts
are queried after source destruction and environment changes. Ordinary grouped,
ordinary generic and generic progressive results match. Existing analytic
render/primitive tests also pass with the option enabled.

The Metal compute suite passes 76 tests in default, ray-only wide, and combined
compact/threaded configurations. Sixteen resource tests pass with Metal API
validation, including enabled ordinary rays, generic progressive profiler units,
independent point/ray settings, no-group fallback and shared allocation cleanup.
The CPU-only suite passes sixteen tests and public headers compile.

All three bounded private views match fresh generic output buffers through public
APIs after source destruction; their primitive counts also match. The rebuilt
normal renderer produces byte-identical public-fixture images with ray groups
off/on. Raw models, cameras, profiles and images remain ignored.
