# SPDX-FileCopyrightText: 2026 Giovanni MARIANO
#
# SPDX-License-Identifier: MPL-2.0

"""Create the machine-readable index for pyAlea GitHub Release assets."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("artifact_root", type=Path)
    parser.add_argument("output", type=Path)
    args = vars(parser.parse_args())

    artifact_root = args["artifact_root"]
    assets = []
    for path in sorted(artifact_root.glob("pyalea-*/*")):
        if not path.is_file() or path.name.endswith(".sha256"):
            continue
        assets.append({
            "filename": path.name,
            "sha256": sha256(path),
            "size": path.stat().st_size,
        })
    if not assets:
        raise RuntimeError(f"no pyAlea archives found below {artifact_root}")

    args["output"].write_text(
        json.dumps({"schema_version": 1, "assets": assets}, indent=2) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
