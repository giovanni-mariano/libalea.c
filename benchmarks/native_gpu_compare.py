#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2026 Giovanni MARIANO
# SPDX-License-Identifier: MPL-2.0
"""Run serialized, alternating native GPU comparisons with bounded processes.

A config contains variants: [{name, argv: [...], env: {...}}]. Each executable
owns its warmup, within-process repetitions and correctness checks. Raw logs,
commands and metadata must remain under an ignored directory inside build/.
This runner records observations; it does not infer a statistically valid speedup.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import platform
import signal
import shutil
import subprocess
import time

ROOT = Path(__file__).resolve().parents[1]


def capture(argv):
    return subprocess.check_output(argv, cwd=ROOT, text=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--blocks', type=int, default=2)
    parser.add_argument('--timeout', type=float, default=90)
    args = parser.parse_args()
    output = args.output.resolve()
    if not output.is_relative_to((ROOT / 'build').resolve()):
        parser.error('output must be inside the ignored repository build directory')
    if args.blocks < 1 or not 0 < args.timeout <= 90:
        parser.error('blocks must be positive; timeout must be in (0, 90] seconds')
    variants = json.loads(args.config.read_text())['variants']
    if not variants or any(not isinstance(v.get('argv'), list) or not v['argv'] or
                           any(not isinstance(a, str) for a in v['argv']) or
                           not isinstance(v.get('name'), str) or
                           not isinstance(v.get('env', {}), dict) or
                           any(not isinstance(k, str) or not isinstance(x, str)
                               for k, x in v.get('env', {}).items()) for v in variants):
        parser.error('each variant requires a name, argv string array and string environment')
    if len({v['name'] for v in variants}) != len(variants):
        parser.error('variant names must be unique')
    if output.exists():
        parser.error('output must be a new directory to preserve earlier evidence')
    ignored = subprocess.run(['git', 'check-ignore', '-q', str(output / 'report.json')], cwd=ROOT)
    if ignored.returncode:
        parser.error('output is not ignored by Git')
    output.mkdir(parents=True)
    # Private paths may occur in a local config or diff: store only in ignored output.
    identities = []
    for variant in variants:
        executable = shutil.which(variant['argv'][0])
        if executable:
            digest = hashlib.sha256()
            with open(executable, 'rb') as source:
                for chunk in iter(lambda: source.read(1024 * 1024), b''):
                    digest.update(chunk)
            identities.append({'path': str(Path(executable).resolve()), 'sha256': digest.hexdigest()})
        else:
            identities.append(None)
    report = {'schema_version': 1, 'platform': platform.platform(),
              'machine': platform.machine(), 'revision': capture(['git', 'rev-parse', 'HEAD']).strip(),
              'variants': variants, 'executables': identities, 'cpu_threads_default': 4,
              'blocks': args.blocks, 'timeout_seconds': args.timeout,
              'runs': [], 'status': 'running'}
    (output / 'working-tree.patch').write_text(capture(['git', 'diff', 'HEAD']))
    (output / 'status.txt').write_text(capture(['git', 'status', '--short']))
    def save():
        (output / 'report.json').write_text(json.dumps(report, indent=2) + '\n')
    save()
    for block in range(args.blocks):
        order = range(len(variants)) if block % 2 == 0 else reversed(range(len(variants)))
        for index in order:
            variant = variants[index]
            env = {k: v for k, v in os.environ.items()
                   if not k.startswith(('ALEA_METAL_', 'ALEA_PLOT_', 'ALEA_GRID_', 'PROBE_', 'MTL_'))}
            env['ALEA_NUM_THREADS'] = '4'
            env.update(variant.get('env', {}))
            start = time.monotonic()
            status = 'complete'
            with (output / f'block-{block}-variant-{index}.log').open('w') as log:
                try:
                    process = subprocess.Popen(variant['argv'], cwd=ROOT, env=env,
                                               stdout=log, stderr=subprocess.STDOUT,
                                               start_new_session=True)
                except OSError as error:
                    log.write(str(error)); status = 'launch_error'; rc = None
                else:
                    try:
                        rc = process.wait(timeout=args.timeout)
                    except subprocess.TimeoutExpired:
                        status = 'timeout'
                        os.killpg(process.pid, signal.SIGKILL)
                        rc = process.wait()
                    except KeyboardInterrupt:
                        os.killpg(process.pid, signal.SIGKILL)
                        process.wait()
                        report['status'] = 'interrupted'; save()
                        raise
            report['runs'].append({'block': block, 'variant': index, 'status': status,
                                   'returncode': rc, 'wall_seconds': time.monotonic() - start})
            save()
            print(f'block={block} variant={index} status={status} returncode={rc}', flush=True)
            if rc != 0 or status != 'complete':
                report['status'] = 'failed'; save(); return 1
    report['status'] = 'complete'; save()
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
