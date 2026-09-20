# SPDX-FileCopyrightText: 2026 Giovanni MARIANO
#
# SPDX-License-Identifier: MPL-2.0

"""Report interpreter build settings used by the binding Makefiles."""

from __future__ import annotations

import argparse
from pathlib import Path
import shlex
import sys
import sysconfig

import numpy


def value(name: str) -> str:
    values = {
        "abi-tag": f"cp{sys.version_info.major}{sys.version_info.minor}",
        "cc": shlex.split(sysconfig.get_config_var("CC") or "cc")[0],
        "ext-suffix": sysconfig.get_config_var("EXT_SUFFIX") or ".so",
        "numpy-include": numpy.get_include(),
        "platform": sysconfig.get_platform().replace("-", "_").replace(".", "_"),
        "python-include": sysconfig.get_path("include"),
        "python-tag": f"cp{sys.version_info.major}{sys.version_info.minor}",
    }
    return values[name]


def write_msvc_responses(compile_path: Path, link_path: Path, root: Path) -> None:
    compile_path.parent.mkdir(parents=True, exist_ok=True)
    python_include = Path(sysconfig.get_path("include")).resolve()
    numpy_include = Path(numpy.get_include()).resolve()
    compile_path.write_text(
        f'/I"{root.resolve() / "include"}"\n'
        f'/I"{python_include}"\n'
        f'/I"{numpy_include}"\n',
        encoding="utf-8",
    )

    library = sysconfig.get_config_var("LDLIBRARY")
    library_dir = sysconfig.get_config_var("LIBDIR")
    if not library_dir:
        library_dir = str(Path(sys.base_prefix) / "libs")
    if not library:
        library = f"python{sys.version_info.major}{sys.version_info.minor}.lib"
    link_path.write_text(
        f'/LIBPATH:"{Path(library_dir).resolve()}"\n{library}\n',
        encoding="utf-8",
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="command", required=True)
    value_parser = subparsers.add_parser("value")
    value_parser.add_argument("name", choices=(
        "abi-tag", "cc", "ext-suffix", "numpy-include", "platform",
        "python-include", "python-tag",
    ))
    msvc_parser = subparsers.add_parser("write-msvc-responses")
    msvc_parser.add_argument("--compile", type=Path, required=True)
    msvc_parser.add_argument("--link", type=Path, required=True)
    msvc_parser.add_argument("--root", type=Path, required=True)
    args = parser.parse_args()

    if args.command == "value":
        print(value(args.name))
    else:
        write_msvc_responses(args.compile, args.link, args.root)


if __name__ == "__main__":
    main()
