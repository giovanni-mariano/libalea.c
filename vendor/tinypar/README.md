<!-- SPDX-FileCopyrightText: 2026 Giovanni MARIANO -->
<!-- SPDX-License-Identifier: MPL-2.0 -->

# tinypar

tinypar is a small C11 library for running an independent index range in
parallel. It provides dynamic chunk scheduling, deterministic static-block
scheduling, caller participation, cooperative callback failure, and a
compile-time serial backend. It is not a task graph, event loop, or persistent
worker pool.

On POSIX platforms tinypar uses pthreads. On Windows it uses `_beginthreadex`
and native Windows thread handles. The public API does not expose either
backend's types.

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

`make benchmark` runs a small serial/static/dynamic comparison. Optional
arguments select the item count, worker limit, and chunk size:

```sh
build/benchmark_parallel_for 10000000 8 4096
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

Each invocation owns its queue and cancellation state, so calls may run
concurrently or recursively. tinypar does not use a global worker pool.
`tinypar_in_parallel()` lets a consumer choose a one-worker configuration for
nested calls when oversubscription is undesirable.

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
    .max_workers = 0,  /* available logical processor count */
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
    .max_workers = tinypar_in_parallel() ? 1 : 0,
    .schedule = TINYPAR_SCHEDULE_STATIC_BLOCK
};
```

The calling thread participates as worker zero. Successful invocations use the
number of workers returned by `tinypar_effective_workers()`. A serial build
always returns one effective worker for non-empty valid work.

## License

tinypar is licensed under the Mozilla Public License, version 2.0. See
`LICENSES/MPL-2.0.txt`.
