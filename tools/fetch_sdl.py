#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2026 Giovanni MARIANO
# SPDX-License-Identifier: MPL-2.0
"""Fetch and build the native SDL2 fallback for optional interactive viewers."""

import argparse
import hashlib
import json
import os
from pathlib import Path
import shlex
import shutil
import subprocess
import sys
import tarfile
import urllib.request

VERSION = "2.32.10"
SHA256 = "5f5993c530f084535c65a6879e9b26ad441169b3e25d789d83287040a9ca5165"
URL = f"https://www.libsdl.org/release/SDL2-{VERSION}.tar.gz"


def digest(path):
    checksum = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            checksum.update(chunk)
    return checksum.hexdigest()


def extract(archive, destination):
    # Validate paths and link targets before extracting (also on Python < 3.12).
    with tarfile.open(archive, "r:gz") as source:
        for member in source.getmembers():
            path = (destination / member.name).resolve()
            path.relative_to(destination)
            if not member.name.startswith(f"SDL2-{VERSION}/"):
                raise RuntimeError("unexpected SDL archive layout")
            if member.issym() or member.islnk():
                base = path.parent if member.issym() else destination
                (base / member.linkname).resolve().relative_to(destination)
            elif not (member.isfile() or member.isdir()):
                raise RuntimeError("unexpected SDL archive member")
        if hasattr(tarfile, "data_filter"):
            source.extractall(destination, filter="data")
        else:
            source.extractall(destination)


def build(args):
    directory = args.directory.resolve()
    prefix = directory / "install"
    ready = directory / "ready.json"
    configuration = {"version": VERSION, "compiler": args.compiler}
    if (ready.is_file() and (prefix / "bin/sdl2-config").is_file()
            and (prefix / "lib/libSDL2.a").is_file()
            and json.loads(ready.read_text()) == configuration):
        return
    if not shutil.which("cmake"):
        raise RuntimeError("SDL2 fallback requires CMake; install CMake or SDL2 development files")
    directory.mkdir(parents=True, exist_ok=True)
    archive = directory / f"SDL2-{VERSION}.tar.gz"
    if not archive.exists():
        print(f"Downloading {URL}", flush=True)
        temporary = archive.with_suffix(".part")
        try:
            with urllib.request.urlopen(URL, timeout=60) as response, temporary.open("wb") as output:
                shutil.copyfileobj(response, output)
            if digest(temporary) != SHA256:
                raise RuntimeError("SDL2 download failed SHA-256 verification")
            temporary.replace(archive)
        finally:
            temporary.unlink(missing_ok=True)
    if digest(archive) != SHA256:
        raise RuntimeError(f"SDL2 archive failed SHA-256 verification: {archive}; remove it and retry")
    source = directory / f"SDL2-{VERSION}"
    extracted = directory / "extracted"
    if not extracted.exists():
        extract(archive, directory)
        extracted.touch()
    cmake_build = directory / "cmake-build"
    compiler = shlex.split(args.compiler)
    if not compiler:
        raise RuntimeError("empty C compiler")
    print(f"Building SDL2 {VERSION} locally in {directory}", flush=True)
    subprocess.run([
        "cmake", "-S", str(source), "-B", str(cmake_build),
        "-DCMAKE_BUILD_TYPE=Release", f"-DCMAKE_INSTALL_PREFIX={prefix}",
        "-DCMAKE_INSTALL_LIBDIR=lib", "-DCMAKE_C_COMPILER=" + ";".join(compiler),
        "-DSDL_SHARED=OFF", "-DSDL_STATIC=ON", "-DSDL_TEST=OFF", "-DSDL_TESTS=OFF",
    ], check=True)
    # The SDL build has its own parallelism, independent of the parent Make jobserver.
    environment = dict(os.environ)
    environment.pop("MAKEFLAGS", None)
    environment.pop("MFLAGS", None)
    subprocess.run(["cmake", "--build", str(cmake_build), "--config", "Release",
                    "--parallel", str(args.jobs)], check=True, env=environment)
    subprocess.run(["cmake", "--install", str(cmake_build), "--config", "Release"], check=True)
    ready.write_text(json.dumps(configuration) + "\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--directory", required=True, type=Path)
    parser.add_argument("--compiler", default="cc")
    parser.add_argument("--jobs", type=int, default=4)
    args = parser.parse_args()
    if args.jobs < 1:
        parser.error("--jobs must be positive")
    try:
        build(args)
    except (OSError, ValueError, RuntimeError, subprocess.CalledProcessError, tarfile.TarError) as error:
        print(f"SDL2 fallback: {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
