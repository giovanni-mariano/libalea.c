#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2026 Giovanni MARIANO
# SPDX-License-Identifier: MPL-2.0
"""Run the Linux benchmark adapter with a private, same-binary regional image."""
import argparse
import errno
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import time


def sha(stream):
    stream.seek(0)
    digest = hashlib.sha256()
    buffer = memoryview(bytearray(1024 * 1024))
    while size := stream.readinto(buffer):
        digest.update(buffer[:size])
    stream.seek(0)
    return digest.hexdigest()


def run():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--cache", type=Path, required=True)
    p.add_argument("--selective", action="store_true")
    p.add_argument("adapter", type=Path)
    p.add_argument("mode", choices=["snapshot"])
    p.add_argument("model")
    p.add_argument("workload")
    p.add_argument("snapshot", type=Path)
    p.add_argument("--point-grid", nargs=6, required=True)
    args = p.parse_args()
    start = time.monotonic()
    args.cache.mkdir(parents=True, exist_ok=True, mode=0o700)
    permissions = args.cache.stat()
    if permissions.st_uid != os.getuid() or permissions.st_mode & 0o022:
        p.error("cache directory must be owned by this user and writable only by its owner")
    env = os.environ.copy()
    for key in ("ALEA_INTERVAL_IMAGE_READ", "ALEA_INTERVAL_IMAGE_WRITE", "ALEA_INTERVAL_SELECTIVE"):
        env.pop(key, None)
    if args.selective:
        env["ALEA_INTERVAL_SELECTIVE"] = "1"
    with args.snapshot.open("rb") as snapshot, args.adapter.open("rb") as adapter:
        key = [sha(adapter), sha(snapshot), [format(float(v), ".17g") for v in args.point_grid], args.selective]
        entry = args.cache / hashlib.sha256(json.dumps(key).encode()).hexdigest()
        fds = [adapter.fileno(), snapshot.fileno()]
        cmd = [f"/proc/self/fd/{fds[0]}", args.mode, args.model, args.workload,
               f"/proc/self/fd/{fds[1]}", "--point-grid", *args.point_grid]
        if entry.exists():
            with (entry / "index").open("rb") as image:
                if sha(image) != (entry / "sha256").read_text().strip():
                    raise ValueError(f"regional cache checksum mismatch in {entry}")
                fd = image.fileno()
                env["ALEA_INTERVAL_IMAGE_READ"] = f"/proc/self/fd/{fd}"
                verify_s = time.monotonic() - start
                print(f"verified regional image in {verify_s:.6f} s", file=sys.stderr, flush=True)
                # Keep the verified files open across exec.
                for inherited in (*fds, fd):
                    os.set_inheritable(inherited, True)
                os.execve(cmd[0], cmd, env)
        with tempfile.TemporaryDirectory(dir=args.cache) as tmp:
            staged = Path(tmp) / "entry"
            staged.mkdir(mode=0o700)
            env["ALEA_INTERVAL_IMAGE_WRITE"] = str(staged / "index")
            proc = subprocess.run(cmd, env=env, pass_fds=fds)
            if proc.returncode:
                return proc.returncode
            with (staged / "index").open("rb") as image:
                (staged / "sha256").write_text(sha(image) + "\n")
            try:
                staged.rename(entry)
            except OSError as e:
                if e.errno not in (errno.EEXIST, errno.ENOTEMPTY):
                    raise
    return 0


if __name__ == "__main__":
    sys.exit(run())
