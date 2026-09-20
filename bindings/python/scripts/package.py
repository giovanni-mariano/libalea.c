# SPDX-FileCopyrightText: 2026 Giovanni MARIANO
#
# SPDX-License-Identifier: MPL-2.0

"""Create a self-contained pyAlea release archive."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import platform
import re
import subprocess
import sys
import sysconfig
import tarfile
import zipfile

import numpy


def libalea_version(root: Path) -> str:
    header = (root / "include" / "alea.h").read_text(encoding="utf-8")
    parts = []
    for name in ("MAJOR", "MINOR", "PATCH"):
        match = re.search(rf"^#define ALEA_VERSION_{name} (\d+)$", header, re.MULTILINE)
        if not match:
            raise RuntimeError(f"ALEA_VERSION_{name} is missing from include/alea.h")
        parts.append(match.group(1))
    return ".".join(parts)


def git_commit(root: Path) -> str:
    result = subprocess.run(
        ["git", "rev-parse", "HEAD"], cwd=root, text=True,
        capture_output=True, check=True,
    )
    return result.stdout.strip()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def include_archive_member(info: tarfile.TarInfo):
    parts = Path(info.name).parts
    if "__pycache__" in parts or info.name.endswith((".pyc", ".pyo")):
        return None
    return info


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--package-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--platform-tag")
    args = parser.parse_args()

    root = args.root.resolve()
    package_dir = args.package_dir.resolve()
    output_dir = args.output_dir.resolve()
    version = libalea_version(root)
    python_tag = f"cp{sys.version_info.major}{sys.version_info.minor}"
    platform_tag = args.platform_tag or sysconfig.get_platform()
    platform_tag = platform_tag.replace("-", "_").replace(".", "_")
    archive_root = f"pyalea-{version}-{python_tag}-{python_tag}-{platform_tag}"

    metadata = {
        "abi_tag": python_tag,
        "compiler": platform.python_compiler(),
        "features": {"mpi": False, "parallel_backend": "tinypar"},
        "libalea_commit": git_commit(root),
        "libalea_version": version,
        "numpy_build_version": numpy.__version__,
        "platform": sysconfig.get_platform(),
        "platform_tag": platform_tag,
        "python_implementation": platform.python_implementation(),
        "python_tag": python_tag,
        "python_version": platform.python_version(),
    }
    (package_dir / "_build_info.json").write_text(
        json.dumps(metadata, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )

    output_dir.mkdir(parents=True, exist_ok=True)
    if sys.platform == "win32":
        archive = output_dir / f"{archive_root}.zip"
        with zipfile.ZipFile(archive, "w", compression=zipfile.ZIP_DEFLATED) as bundle:
            for path in sorted(package_dir.rglob("*")):
                if path.is_file() and "__pycache__" not in path.parts and path.suffix not in (".pyc", ".pyo"):
                    bundle.write(path, Path(archive_root) / "pyalea" / path.relative_to(package_dir))
            for path in sorted((root / "LICENSES").glob("*")):
                bundle.write(path, Path(archive_root) / "LICENSES" / path.name)
            bundle.write(root / "bindings" / "python" / "README.md", Path(archive_root) / "README.md")
    else:
        archive = output_dir / f"{archive_root}.tar.gz"
        with tarfile.open(archive, "w:gz") as bundle:
            bundle.add(
                package_dir, arcname=str(Path(archive_root) / "pyalea"),
                filter=include_archive_member,
            )
            bundle.add(root / "LICENSES", arcname=str(Path(archive_root) / "LICENSES"))
            bundle.add(
                root / "bindings" / "python" / "README.md",
                arcname=str(Path(archive_root) / "README.md"),
            )

    checksum = archive.with_name(archive.name + ".sha256")
    checksum.write_text(f"{sha256(archive)}  {archive.name}\n", encoding="ascii")
    print(archive)


if __name__ == "__main__":
    main()
