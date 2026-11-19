#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2026 Giovanni MARIANO
# SPDX-License-Identifier: MPL-2.0
"""CLI profile and JSON contract checks; no third-party dependencies."""
import json
from pathlib import Path
import subprocess
import sys
import tempfile

exe = str(Path(sys.argv[1] if len(sys.argv) > 1 else 'bin/mc_compare').resolve())
with tempfile.TemporaryDirectory() as directory:
    paths = [Path(directory) / name for name in ('before.i', 'after.i')]
    for path, fractions in zip(paths, ('1001 1 8016 2', '1001 2 8016 4')):
        path.write_text(f'profiles\n1 1 -1 -1\n\n1 so 5\n\nm1 {fractions}\n')
    for profile, code in [('legacy', 1), ('semantic', 0)]:
        result = subprocess.run([exe, '--json', '--profile', profile, *map(str, paths)],
                                capture_output=True, text=True)
        assert result.returncode == code, result.stderr
        report = json.loads(result.stdout)
        assert report['profile'] == profile
        assert report['version'] == (2 if profile == 'legacy' else 3)
        assert {c['change'] for c in report['changes']} == ({'content'} if code else {'unchanged'})
    for args in [('--profile',), ('--profile', 'unknown'), ('--unknown',)]:
        result = subprocess.run([exe, *args], capture_output=True, text=True)
        assert result.returncode == 2
print('Comparison CLI profile checks passed')
