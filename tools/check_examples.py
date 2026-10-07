#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2026 Giovanni MARIANO
# SPDX-License-Identifier: MPL-2.0
"""Smoke-test examples without external nuclear data in a temporary directory.

python3 tools/check_examples.py --lua bin/alea
PYTHONPATH=<binding package directory> python3 tools/check_examples.py --python
"""
import argparse
import os
from pathlib import Path
import subprocess
import sys
import tempfile


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--lua", type=Path, help="Path to the alea CLI")
    parser.add_argument("--python", action="store_true", help="Test the imported pyalea package")
    args = parser.parse_args()
    if not args.lua and not args.python:
        parser.error("select --lua and/or --python")
    root = Path(__file__).resolve().parents[1]
    failures = []
    with tempfile.TemporaryDirectory(prefix="alea-examples-") as directory:
        env = dict(os.environ, ALEA_EXAMPLE_OUTPUT=directory)
        commands = []
        if args.lua:
            cli = str(args.lua.resolve())
            for path in sorted((root / "examples/lua").glob("*.lua")):
                if int(path.name[:2]) >= 19:  # Require user-supplied ACE data.
                    continue
                command = [cli, str(path)]
                if "Standalone: yes" not in path.read_text():
                    command.append(str(Path(directory) / "alea_shielded_sphere.i"))
                commands.append((path.name, command))
        if args.python:
            for path in sorted((root / "examples/python").glob("[0-9][0-9]_*.py")):
                commands.append((path.name, [sys.executable, str(path), directory]))
        for name, command in commands:
            try:
                result = subprocess.run(command, env=env, capture_output=True,
                                        text=True, timeout=60)
                if result.returncode:
                    failures.append(name)
                    print(result.stdout + result.stderr)
                print(f"{name}: {'PASS' if not result.returncode else 'FAIL'}", flush=True)
            except subprocess.TimeoutExpired:
                failures.append(name)
                print(f"{name}: TIMEOUT", flush=True)
    return bool(failures)


if __name__ == "__main__":
    sys.exit(main())
