// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/* Estimate every concrete cell-instance volume in an MCNP or OpenMC model. */

#include "alea.h"
#include "alea_cluster.h"
#include "alea_mcnp.h"
#include "alea_openmc.h"

#include <errno.h>
#include <limits.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

typedef enum { FORMAT_AUTO, FORMAT_MCNP, FORMAT_OPENMC } input_format_t;

typedef struct {
    const char* input_path;
    const char* output_path;
    input_format_t format;
    size_t rays;
    size_t batch_size;
    size_t workers;
    size_t worker_memory_mib;
    uint64_t seed;
    double target_rel_error;
    double center[3];
    double radius;
    int csv;
} arguments_t;

static void usage(FILE* stream, const char* program) {
    fprintf(stream,
        "Usage: %s [options] INPUT\n"
        "\nEstimate every concrete cell-instance volume in MCNP or OpenMC input.\n"
        "Run directly for one rank or with mpiexec for distributed execution.\n"
        "\nOptions:\n"
        "  --format auto|mcnp|openmc  Input format (default: extension)\n"
        "  --rays N                   Maximum global rays (required)\n"
        "  --radius R                 Sampling-sphere radius (required)\n"
        "  --center X Y Z             Sampling-sphere center (default: 0 0 0)\n"
        "  --batch N                  Global rays per reduction (default: 10000)\n"
        "  --workers N                TinyPar workers per rank (default: automatic)\n"
        "  --worker-memory-mib N      Total worker scratch per rank (default: 256)\n"
        "  --seed N                   Sampling seed (default: 42)\n"
        "  --target-rel-error X       Stop when every path reaches X (default: off)\n"
        "  --csv                      Write machine-readable CSV\n"
        "  -o, --output FILE          Write rank-zero output to FILE\n"
        "  -h, --help                 Show this help\n", program);
}

static int ends_with(const char* text, const char* suffix) {
    size_t text_length = strlen(text), suffix_length = strlen(suffix);
    return text_length >= suffix_length &&
        strcmp(text + text_length - suffix_length, suffix) == 0;
}

static int parse_size(const char* text, int allow_zero, size_t* output) {
    if (!*text || strspn(text, "0123456789") != strlen(text)) return -1;
    char* end = NULL;
    errno = 0;
    unsigned long long value = strtoull(text, &end, 10);
    if (errno || end == text || *end || (!allow_zero && value == 0) ||
        value > (unsigned long long)SIZE_MAX) return -1;
    *output = (size_t)value;
    return 0;
}

static int parse_u64(const char* text, uint64_t* output) {
    if (!*text || strspn(text, "0123456789") != strlen(text)) return -1;
    char* end = NULL;
    errno = 0;
    unsigned long long value = strtoull(text, &end, 10);
    if (errno || end == text || *end) return -1;
    *output = (uint64_t)value;
    return 0;
}

static int parse_double(const char* text, double* output) {
    char* end = NULL;
    errno = 0;
    double value = strtod(text, &end);
    if (errno || end == text || *end || !isfinite(value)) return -1;
    *output = value;
    return 0;
}

static int parse_arguments(int argc, char** argv, arguments_t* arguments) {
    *arguments = (arguments_t){
        .format = FORMAT_AUTO, .batch_size = 10000, .seed = 42
    };
    for (int i = 1; i < argc; ++i) {
        const char* option = argv[i];
        if (!strcmp(option, "-h") || !strcmp(option, "--help")) return 1;
        if (!strcmp(option, "--csv")) { arguments->csv = 1; continue; }
        if (!strcmp(option, "--center")) {
            if (i + 3 >= argc || parse_double(argv[i + 1], &arguments->center[0]) ||
                parse_double(argv[i + 2], &arguments->center[1]) ||
                parse_double(argv[i + 3], &arguments->center[2])) return -1;
            i += 3;
            continue;
        }
        if (!strcmp(option, "--format")) {
            if (++i >= argc) return -1;
            if (!strcmp(argv[i], "auto")) arguments->format = FORMAT_AUTO;
            else if (!strcmp(argv[i], "mcnp")) arguments->format = FORMAT_MCNP;
            else if (!strcmp(argv[i], "openmc")) arguments->format = FORMAT_OPENMC;
            else return -1;
            continue;
        }
        if (!strcmp(option, "--rays") || !strcmp(option, "--radius") ||
            !strcmp(option, "--batch") ||
            !strcmp(option, "--workers") ||
            !strcmp(option, "--worker-memory-mib") || !strcmp(option, "--seed") ||
            !strcmp(option, "--target-rel-error") || !strcmp(option, "-o") ||
            !strcmp(option, "--output")) {
            if (++i >= argc) return -1;
            if (!strcmp(option, "--rays") &&
                parse_size(argv[i], 0, &arguments->rays)) return -1;
            if (!strcmp(option, "--radius") &&
                (parse_double(argv[i], &arguments->radius) ||
                 arguments->radius <= 0.0)) return -1;
            if (!strcmp(option, "--batch") &&
                parse_size(argv[i], 0, &arguments->batch_size)) return -1;
            if (!strcmp(option, "--workers") &&
                (parse_size(argv[i], 1, &arguments->workers) ||
                 arguments->workers > INT_MAX)) return -1;
            if (!strcmp(option, "--worker-memory-mib") &&
                (parse_size(argv[i], 0, &arguments->worker_memory_mib) ||
                 arguments->worker_memory_mib > SIZE_MAX / (1024u * 1024u)))
                return -1;
            if (!strcmp(option, "--seed") &&
                parse_u64(argv[i], &arguments->seed)) return -1;
            if (!strcmp(option, "--target-rel-error") &&
                (parse_double(argv[i], &arguments->target_rel_error) ||
                 arguments->target_rel_error < 0.0 ||
                 arguments->target_rel_error > 1.0)) return -1;
            if (!strcmp(option, "-o") || !strcmp(option, "--output"))
                arguments->output_path = argv[i];
            continue;
        }
        if (option[0] == '-' || arguments->input_path) return -1;
        arguments->input_path = option;
    }
    if (!arguments->input_path || arguments->rays == 0 ||
        arguments->radius <= 0.0) return -1;
    if (arguments->batch_size > arguments->rays)
        arguments->batch_size = arguments->rays;
    if (arguments->format == FORMAT_AUTO)
        arguments->format = ends_with(arguments->input_path, ".xml")
            ? FORMAT_OPENMC : FORMAT_MCNP;
    return 0;
}

typedef struct {
    time_t started;
    alea_system_t* sys;
    const double* errors;
    size_t count;
} progress_t;

static int report_progress(size_t completed, size_t maximum,
                           double maximum_error, void* opaque) {
    const progress_t* progress = opaque;
    double elapsed = difftime(time(NULL), progress->started);
    size_t unseen = 0, worst = 0;
    double worst_error = -1.0;
    for (size_t i = 0; i < progress->count; ++i) {
        if (progress->errors[i] < 0.0) ++unseen;
        else if (progress->errors[i] > worst_error) {
            worst_error = progress->errors[i];
            worst = i;
        }
    }
    fprintf(stderr, "[volume] rays=%zu/%zu (%.1f%%) elapsed=%.0fs",
            completed, maximum, 100.0 * (double)completed / (double)maximum,
            elapsed);
    if (elapsed > 0.0 && completed)
        fprintf(stderr, " rate=%.1f rays/s ETA-to-ray-limit=%.0fs",
                (double)completed / elapsed,
                elapsed * (double)(maximum - completed) / (double)completed);
    if (isfinite(maximum_error))
        fprintf(stderr, " max-rel-error=%.5g", maximum_error);
    else
        fputs(" max-rel-error=unavailable", stderr);
    fprintf(stderr, " unsampled=%zu/%zu", unseen, progress->count);
    if (worst_error >= 0.0) {
        alea_volume_path_t path;
        if (alea_volume_paths_get_range(progress->sys, worst, &path, 1) == 1)
            fprintf(stderr, " worst-sampled-path=%llu cell=%d rel-error=%.5g",
                    (unsigned long long)path.path_id,
                    path.terminal_cell_id, worst_error);
        else
            fprintf(stderr, " worst-sampled-path=%zu rel-error=%.5g",
                    worst, worst_error);
    }
    fputc('\n', stderr);
    fflush(stderr);
    return 0;
}

static int cell_id_from_index(const alea_system_t* sys, int cell_index) {
    alea_cell_info_t info;
    if (cell_index < 0 ||
        alea_cell_get_info(sys, (size_t)cell_index, &info) != 0) return -1;
    return info.cell_id;
}

static void print_instance(FILE* output, const alea_system_t* sys,
                           const alea_volume_path_t* path) {
    fputs("root", output);
    for (size_t i = 0; i < path->ancestor_count; ++i) {
        int cell_id = cell_id_from_index(sys, path->ancestor_cell_indices[i]);
        if (cell_id >= 0)
            fprintf(output, "/u%d:c%d", path->ancestor_universe_ids[i], cell_id);
        else
            fprintf(output, "/u%d:cidx%d", path->ancestor_universe_ids[i],
                    path->ancestor_cell_indices[i]);
    }
    for (size_t i = 0; i < path->lattice_step_count; ++i) {
        const alea_volume_lattice_step_t* step = &path->lattice_steps[i];
        int cell_id = cell_id_from_index(sys, step->lattice_cell_index);
        fprintf(output, "/lat-c%d[%d:%d:%d]=>u%d",
                cell_id >= 0 ? cell_id : step->lattice_cell_index,
                step->i, step->j, step->k, step->fill_universe);
    }
    fprintf(output, "/u%d:c%d", path->universe_id, path->terminal_cell_id);
}

static int write_report(FILE* output, int csv,
                        alea_system_t* sys,
                        const double* volumes, const double* errors,
                        size_t count) {
    enum { PATH_CHUNK = 256 };
    alea_volume_path_t* paths = calloc(PATH_CHUNK, sizeof(*paths));
    if (!paths) return -1;
    if (csv) {
        fputs("path_id,cell_id,material_id,universe_id,depth,instance,"
              "world_to_local_tx,world_to_local_ty,world_to_local_tz,"
              "volume,relative_error\n", output);
    } else {
        fprintf(output, "%6s %8s %8s %8s %5s %16s %10s  %s\n",
                "Path", "Cell", "Material", "Universe", "Depth", "Volume",
                "Rel.err", "Instance");
    }
    for (size_t base = 0; base < count; base += PATH_CHUNK) {
        size_t wanted = count - base;
        if (wanted > PATH_CHUNK) wanted = PATH_CHUNK;
        if (alea_volume_paths_get_range(sys, base, paths, wanted) != wanted) {
            free(paths);
            return -1;
        }
        for (size_t j = 0; j < wanted; ++j) {
            const size_t i = base + j;
            const alea_volume_path_t* path = &paths[j];
            if (csv) {
                fprintf(output, "%llu,%d,%d,%d,%d,\"",
                        (unsigned long long)path->path_id,
                        path->terminal_cell_id, path->material_id,
                        path->universe_id, path->depth);
                print_instance(output, sys, path);
                fprintf(output, "\",%.17g,%.17g,%.17g,%.17g,",
                        path->world_to_local[3], path->world_to_local[7],
                        path->world_to_local[11], volumes[i]);
                if (errors[i] >= 0.0) fprintf(output, "%.17g", errors[i]);
                fputc('\n', output);
            } else {
                fprintf(output, "%6llu %8d %8d %8d %5d ",
                        (unsigned long long)path->path_id,
                        path->terminal_cell_id, path->material_id,
                        path->universe_id, path->depth);
                fprintf(output, "%16.8e ", volumes[i]);
                if (errors[i] >= 0.0) fprintf(output, "%10.4g  ", errors[i]);
                else fprintf(output, "%10s  ", "-");
                print_instance(output, sys, path);
                fputc('\n', output);
            }
        }
    }
    free(paths);
    return 0;
}

int main(int argc, char** argv) {
    /* Batch schedulers redirect stderr; keep diagnostics immediately visible. */
    setvbuf(stderr, NULL, _IONBF, 0);
    int result = 1;
    mcnp_model_t* mcnp_model = NULL;
    openmc_model_t* openmc_model = NULL;
    alea_cluster_t* cluster = NULL;
    double* volumes = NULL;
    double* errors = NULL;
    char* input_data = NULL;
    size_t input_length = 0;
    FILE* output = NULL;

    alea_cluster_status_t status = alea_cluster_initialize(&argc, &argv);
    if (status != ALEA_CLUSTER_OK) {
        fprintf(stderr, "cluster initialization: %s\n",
                alea_cluster_status_string(status));
        return 1;
    }
    cluster = alea_cluster_create();
    if (!cluster) { fprintf(stderr, "cluster context creation failed\n"); goto done; }
    int rank = alea_cluster_rank(cluster);

    arguments_t arguments;
    int parsed = parse_arguments(argc, argv, &arguments);
    if (parsed != 0) {
        if (rank == 0) usage(parsed > 0 ? stdout : stderr, argv[0]);
        result = parsed > 0 ? 0 : 2;
        goto done;
    }

    if (rank == 0)
        fprintf(stderr, "[volume] backend=%s ranks=%d rays=%zu batch=%zu "
                "requested-workers/rank=%zu (0=automatic)\n"
                "[volume] reading %s\n",
                alea_cluster_backend(cluster), alea_cluster_size(cluster),
                arguments.rays, arguments.batch_size, arguments.workers,
                arguments.input_path);
    status = arguments.format == FORMAT_OPENMC
        ? alea_cluster_read_file(cluster, arguments.input_path,
                                 &input_data, &input_length)
        : alea_cluster_read_mcnp_input(cluster, arguments.input_path,
                                       &input_data, &input_length);
    if (status != ALEA_CLUSTER_OK) {
        if (rank == 0)
            fprintf(stderr, "cannot read %s: %s\n", arguments.input_path,
                    alea_cluster_status_string(status));
        goto done;
    }

    fprintf(stderr, "[volume rank %d] parsing %zu input bytes\n", rank, input_length);
    alea_system_t* sys = NULL;
    if (arguments.format == FORMAT_OPENMC) {
        openmc_model = openmc_load_string(input_data, input_length);
        if (openmc_model) sys = openmc_model_system(openmc_model);
    } else {
        mcnp_model = mcnp_load_string(input_data, input_length);
        if (mcnp_model) sys = mcnp_model_system(mcnp_model);
    }
    free(input_data);
    input_data = NULL;
    status = alea_cluster_agree(cluster, sys ? ALEA_CLUSTER_OK
                                           : ALEA_CLUSTER_COMPUTE_ERROR);
    if (!sys) {
        fprintf(stderr, "rank %d: cannot load %s: %s\n", rank,
                arguments.input_path, alea_error());
    }
    if (status != ALEA_CLUSTER_OK) {
        if (rank == 0 && sys)
            fprintf(stderr, "model parsing failed on another rank\n");
        goto done;
    }

    fprintf(stderr, "[volume rank %d] enumerating physical cell instances\n", rank);
    size_t path_count = alea_volume_path_count(sys);
    if (path_count) {
        volumes = calloc(path_count, sizeof(*volumes));
        errors = calloc(path_count, sizeof(*errors));
    }
    status = alea_cluster_agree(cluster,
        path_count == 0 ? ALEA_CLUSTER_COMPUTE_ERROR :
        (!volumes || !errors) ? ALEA_CLUSTER_OUT_OF_MEMORY : ALEA_CLUSTER_OK);
    if (status != ALEA_CLUSTER_OK) {
        if (rank == 0)
            fprintf(stderr, "cannot enumerate volume paths: %s\n",
                    alea_cluster_status_string(status));
        goto done;
    }

    alea_volume_estimate_options_t options;
    alea_volume_estimate_options_init(&options);
    options.max_rays = arguments.rays;
    options.batch_size = arguments.batch_size;
    options.requested_workers = arguments.workers;
    if (arguments.worker_memory_mib)
        options.max_parallel_scratch_bytes =
            arguments.worker_memory_mib * 1024u * 1024u;
    options.seed = arguments.seed;
    options.target_rel_error = arguments.target_rel_error;
    options.use_sampling_sphere = true;
    options.sampling_center[0] = arguments.center[0];
    options.sampling_center[1] = arguments.center[1];
    options.sampling_center[2] = arguments.center[2];
    options.sampling_radius = arguments.radius;
    progress_t progress = {time(NULL), sys, errors, path_count};
    options.progress = report_progress;
    options.progress_user_data = &progress;
    fprintf(stderr, "[volume rank %d] %zu paths; preparing caches and sampling\n",
            rank, path_count);
    if (rank == 0) {
        const size_t worker_bytes_per_path =
            3u * sizeof(double) + sizeof(size_t);
        const size_t base_bytes_per_path = 4u * sizeof(double);
        const double path_table_mib =
            path_count <= SIZE_MAX / sizeof(alea_volume_path_t)
                ? (double)(path_count * sizeof(alea_volume_path_t)) /
                    (1024.0 * 1024.0) : INFINITY;
        const double worker_mib =
            path_count <= SIZE_MAX / worker_bytes_per_path
                ? (double)(path_count * worker_bytes_per_path) /
                    (1024.0 * 1024.0) : INFINITY;
        const double base_mib = path_count <= SIZE_MAX / base_bytes_per_path
            ? (double)(path_count * base_bytes_per_path) /
                (1024.0 * 1024.0) : INFINITY;
        fprintf(stderr, "[volume] path metadata table=%.1f MiB/rank; "
                "dense rank arrays=%.1f MiB; "
                "scratch=%.1f MiB/worker; "
                "worker-scratch-limit=%.1f MiB/rank\n",
                path_table_mib, base_mib, worker_mib,
                (double)options.max_parallel_scratch_bytes /
                    (1024.0 * 1024.0));
        fputs("[volume] progress follows each global batch, after all ranks finish; "
              "use a smaller --batch for more frequent updates\n", stderr);
        if (arguments.target_rel_error > 0.0)
            fprintf(stderr, "[volume] target-rel-error=%.5g requires every path; "
                    "unsampled paths prevent early convergence\n",
                    arguments.target_rel_error);
    }
    alea_cluster_volume_stats_t stats;
    status = alea_cluster_estimate_volumes(
        cluster, sys, &options, volumes, errors, &stats);
    if (status != ALEA_CLUSTER_OK) {
        if (rank == 0) fprintf(stderr, "volume estimation: %s\n",
                               alea_cluster_status_string(status));
        goto done;
    }

    if (rank == 0) {
        fprintf(stderr, "[volume] finished: %s; rays=%zu root-workers=%zu; "
                "writing report to %s\n",
                stats.volume.converged ? "target reached" : "ray limit reached",
                stats.volume.rays_completed, stats.local_workers,
                arguments.output_path ? arguments.output_path : "stdout");
        fprintf(stderr, "[volume] root worker scratch=%.1f MiB (%zu workers, "
                "limit %.1f MiB)\n",
                (double)stats.volume.parallel_scratch_bytes /
                    (1024.0 * 1024.0), stats.local_workers,
                (double)stats.volume.parallel_scratch_limit_bytes /
                    (1024.0 * 1024.0));
        if (arguments.target_rel_error > 0.0 && !stats.volume.converged)
            fputs("[volume] warning: requested relative error was not reached\n", stderr);
        output = stdout;
        if (arguments.output_path) {
            output = fopen(arguments.output_path, "w");
            if (!output) {
                fprintf(stderr, "cannot open %s: %s\n", arguments.output_path,
                        strerror(errno));
                goto done;
            }
        }
        if (!arguments.csv)
            fprintf(output, "# backend=%s ranks=%d rays=%zu paths=%zu seed=%llu "
                    "sphere=(%.9g,%.9g,%.9g; %.9g)%s\n",
                    alea_cluster_backend(cluster), stats.rank_count,
                    stats.volume.rays_completed, path_count,
                    (unsigned long long)stats.volume.seed,
                    arguments.center[0], arguments.center[1], arguments.center[2],
                    arguments.radius,
                    stats.volume.converged ? " converged" : "");
        if (write_report(output, arguments.csv, sys, volumes, errors,
                         path_count) != 0) {
            fprintf(stderr, "cannot stream volume path metadata\n");
            goto done;
        }
        if (ferror(output)) { fprintf(stderr, "error writing report\n"); goto done; }
        if (output != stdout && fclose(output) != 0) {
            output = NULL;
            fprintf(stderr, "error closing volume report\n");
            goto done;
        }
        output = NULL;
    }
    result = 0;

done:
    free(input_data);
    if (output && output != stdout) fclose(output);
    free(volumes);
    free(errors);
    openmc_model_destroy(openmc_model);
    mcnp_model_destroy(mcnp_model);
    alea_cluster_destroy(cluster);
    if (alea_cluster_finalize() != ALEA_CLUSTER_OK) result = 1;
    return result;
}
