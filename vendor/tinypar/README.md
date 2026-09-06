<!-- SPDX-FileCopyrightText: 2026 Giovanni MARIANO -->
<!-- SPDX-License-Identifier: MPL-2.0 -->

# tinypar

tinypar is a small C11 library for running an independent index range in
parallel. It provides dynamic chunk scheduling, deterministic static-block
scheduling, caller participation, cooperative callback failure, an optional
reusable executor, and a compile-time serial backend. It is not a task graph,
event loop, or general OpenMP implementation.

On POSIX platforms tinypar uses pthreads. On Windows it uses `_beginthreadex`
and native Windows thread handles. The public API does not expose either
backend's types. Windows workers are distributed across active processor groups
when the machine has more than one group, so the automatic worker default can
use all active logical processors.

## Build and test

```sh
make
make test
make test-serial
make benchmark
```

Runnable examples live in `examples/`, each demonstrating one usage pattern:

| File | Shows |
| --- | --- |
| `array_map.c` | Minimal parallel-for over an array (independent writes). |
| `pi_integration.c` | Reduction via per-worker partial sums. |
| `parallel_find.c` | Cooperative cancellation using a non-zero callback result. |
| `mandelbrot.c` | Dynamic scheduling balancing an uneven per-item workload. |
| `neutron_slab.c` | A reproducible Monte Carlo simulation, one history per index. |

Build them all (binaries land in `build/`):

```sh
make examples
```

`make benchmark` runs a serial/static/dynamic/executor comparison, including
repeated short jobs. Optional arguments select the item count, worker limit,
chunk size, repetition count, and short-job item count:

```sh
build/benchmark_parallel_for 10000000 8 4096 1000 64
```

With an MSVC developer prompt:

```text
nmake /f Makefile.msvc
nmake /f Makefile.msvc test
nmake /f Makefile.msvc clean
nmake /f Makefile.msvc TINYPAR_THREADS=0 test
nmake /f Makefile.msvc examples
```

POSIX consumers must compile and link with `-pthread`:

```sh
cc -std=c11 -Iinclude examples/array_map.c lib/libtinypar.a -pthread
```

For a library with no operating-system thread dependency, build with
`TINYPAR_THREADS=0`. The API remains the same and every invocation uses worker
zero:

```sh
make clean
make TINYPAR_THREADS=0
```

## API model

`tinypar_parallel_for()` partitions `[0, item_count)` into half-open ranges.
With `TINYPAR_SCHEDULE_DYNAMIC`, workers claim `chunk_size` items at a time.
With `TINYPAR_SCHEDULE_STATIC_BLOCK`, every worker receives one deterministic,
contiguous, balanced range; `chunk_size` is then a grain hint used to cap the
worker count. Static callers that want up to one worker per item can use a
`chunk_size` of one.

The callback receives a stable worker index. A dynamic callback may run more
than once for the same worker. A non-zero callback result cancels work that has
not started and returns `TINYPAR_CALLBACK_FAILED` after all workers join.

Each one-shot top-level invocation owns its queue and cancellation state, so
independent calls may run concurrently. Native workers report ready at a start
gate before the complete team is released. Dynamic scheduling reserves one
initial chunk for each effective worker, then distributes all remaining chunks
on demand. Every worker therefore participates in a successful job that has at
least one chunk per worker.

Nested calls made from a tinypar callback execute serially automatically. The
public `tinypar_in_parallel()` query reports whether the current callback
belongs to a multi-worker operation; it remains false in a one-worker callback.
`tinypar_in_callback()` reports lexical callback execution in both cases and is
useful for consumers that must size nested scratch for the actual serial path.

```c
static int process(void* context, size_t worker, size_t begin, size_t end) {
    (void)worker;
    for (size_t i = begin; i < end; i++) {
        /* Process item i. */
    }
    return 0;
}

tinypar_config_t config = {
    .item_count = count,
    .chunk_size = 16,
    .max_workers = 0,  /* process worker default */
    .schedule = TINYPAR_SCHEDULE_DYNAMIC
};
tinypar_status_t status = tinypar_parallel_for(&config, process, context);
```

Static-block scheduling is useful for per-worker floating-point reductions
whose partition and final reduction order must remain reproducible:

```c
tinypar_config_t config = {
    .item_count = count,
    .chunk_size = 1024,
    .max_workers = 0,
    .schedule = TINYPAR_SCHEDULE_STATIC_BLOCK
};
```

Nested invocations do not need to override `max_workers`; TinyPar applies the
serial fallback before choosing a team.

The calling thread participates as worker zero. Successful invocations use the
number of workers returned by `tinypar_effective_workers()`. A serial build
always returns one effective worker for non-empty valid work.

## Reusable executor

Repeated short operations should use an explicit executor so native thread
creation is paid once:

```c
tinypar_executor_config_t executor_config;
tinypar_executor_config_init(&executor_config);
executor_config.max_workers = 8;

tinypar_executor_t* executor = NULL;
tinypar_status_t status = tinypar_executor_create(
    &executor_config, &executor);
if (status == TINYPAR_OK) {
    status = tinypar_executor_parallel_for(
        executor, &config, process, context);
    tinypar_status_t destroy_status = tinypar_executor_destroy(&executor);
    if (status == TINYPAR_OK) status = destroy_status;
}
```

The executor has a stable participant count returned by
`tinypar_executor_workers()`. It can be smaller than requested if only part of
the native worker team could be created. Worker indices remain stable and are
always smaller than that count. Calls on one executor serialize; separate
executors may run concurrently. A job with a smaller per-call worker limit
wakes only that job's active workers.

Executor destruction must not race with new submissions. If destruction
returns an error, the pointer remains non-null for a later destruction attempt,
except when all workers terminated and only native handle cleanup failed; in
that case the executor is safely released and the pointer is null.

On POSIX, an executor must not be used or destroyed in a child created by
`fork()` after the executor's workers exist. The child should call `exec()` or
call `tinypar_executor_abandon_after_fork()` before creating a fresh executor.
Abandoning only clears the child copy of the pointer; it intentionally does not
touch synchronization state inherited from vanished worker threads.

## Worker defaults and failures

`tinypar_set_default_workers(n)` sets the process default used when a config or
executor requests zero workers. Passing zero restores the hardware-derived
default. An explicit nonzero `max_workers` always wins. This controls team size,
not the aggregate size of unrelated executors or concurrent one-shot calls.

All parallel-for calls are synchronous: after they return, no worker can access
the callback or its context. Callback failure cancels work not yet started and
waits for every participant. A platform failure that makes worker termination
unknowable cannot safely satisfy this contract and is treated as an
unrecoverable runtime invariant failure: TinyPar prints a diagnostic and
aborts rather than returning while callback state may still be live.

## License

tinypar is licensed under the Mozilla Public License, version 2.0. See
`LICENSES/MPL-2.0.txt`.
