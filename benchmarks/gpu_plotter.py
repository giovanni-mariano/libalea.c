#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2026 Giovanni MARIANO
# SPDX-License-Identifier: MPL-2.0
"""Compare complete CPU and native Metal plotter slices on public geometry.

Build bin/mc_plotter with USE_METAL=1 first. Each trial starts a fresh process;
query timings include candidate preparation and all frame batches. Process wall
also includes model load, acceleration preparation, drawing and image output.
This is not the independent CPU compute enumeration benchmark.
"""
import argparse
import os
from pathlib import Path
import re
import statistics
import subprocess
import tempfile
import time


def model():
    cells, surfaces = [], ['1 so .75']
    for y in range(8):
        for x in range(8):
            i, cx, cy = 8*y+x, 2*x-7, 2*y-7
            cells.append(f'{1000+i} 0 -{100+i} fill=1 ({cx} {cy} 0)')
            surfaces.append(f'{100+i} rpp {cx-1} {cx+1} {cy-1} {cy+1} -1 1')
    cells.extend(['42 7 -1 -1 u=1', '43 0 1 u=1'])
    return 'Public translated spheres and void\n' + '\n'.join(cells) + '\n\n' + '\n'.join(surfaces) + '\n\n'


def main():
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plotter', type=Path, default=root / 'bin/mc_plotter')
    parser.add_argument('--trials', type=int, default=4)
    parser.add_argument('--size', type=int, default=512)
    parser.add_argument('--threads', type=int, default=4)
    args = parser.parse_args()
    if args.trials < 1 or args.threads < 1 or not 16 <= args.size <= 16384 or args.size % 16:
        parser.error('positive trials/threads and size a multiple of 16 in 16..16384 required')
    timings = {m: [] for m in ('cpu-grid', 'cpu-ray', 'metal-scan', 'metal-candidates', 'metal-selected', 'metal-selected-candidates')}
    environment = {k: v for k, v in os.environ.items()
                   if not k.startswith(('ALEA_METAL_', 'ALEA_PLOT_', 'ALEA_GRID_'))}
    environment['ALEA_NUM_THREADS'] = str(args.threads)
    reference = None
    with tempfile.TemporaryDirectory(prefix='alea-plotter-bench-') as directory:
        tmp = Path(directory)
        source = tmp / 'public.mcnp'
        source.write_text(model())
        for trial in range(args.trials):
            modes = list(timings)
            if trial % 2:
                modes.reverse()
            for mode in modes:
                env = dict(environment, ALEA_METAL_POINT_CANDIDATES=str(int(mode.endswith('candidates'))))
                options = ['--backend=metal'] if mode.startswith('metal') else ['--backend=cpu']
                if mode.startswith('metal-selected'):
                    options.append('--ownership=selected')
                if mode == 'cpu-ray':
                    options.append('--method=ray')
                output = tmp / 'slice.bmp'
                command = [str(args.plotter.resolve()), str(source), 'Z', '0', '-8', '8', '-8', '8',
                           str(args.size), str(output), '--color=materials', *options]
                start = time.perf_counter()
                result = subprocess.run(command, env=env, capture_output=True, text=True, timeout=120)
                elapsed = time.perf_counter() - start
                if result.returncode:
                    raise RuntimeError(f'{mode}: {result.stdout}\n{result.stderr}')
                match = re.search(r'query: ([0-9.]+) ms', result.stdout)
                if not match:
                    raise RuntimeError(f'Missing query timing: {mode}')
                timings[mode].append((elapsed, float(match.group(1))))
                image = output.read_bytes()
                if reference is None:
                    reference = image
                if image != reference:
                    raise RuntimeError(f'Image mismatch: {mode}, trial {trial}')
    print(f'size={args.size} trials={args.trials} cpu_threads={args.threads} profiling=off')
    for mode, values in timings.items():
        query = [v[1] for v in values]
        print(f'mode={mode} process_median_s={statistics.median(v[0] for v in values):.4f} '
              f'query_median_ms={statistics.median(query):.3f} '
              f'query_min_ms={min(query):.3f} query_max_ms={max(query):.3f}')
    print('Every complete image matches across modes and trials.')


if __name__ == '__main__':
    main()
