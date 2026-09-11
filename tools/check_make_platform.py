#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2026 Giovanni MARIANO
# SPDX-License-Identifier: MPL-2.0

"""Check GNU Make platform selection without a target compiler or build."""

import os
from pathlib import Path
import subprocess


def main():
    root = Path(__file__).resolve().parents[1]
    names = (
        "WINDOWS_GNU", "TINYPAR_PLATFORM_SRC", "EXEEXT", "PICFLAGS",
        "CFLAGS", "LDFLAGS", "LUA_PLAT_FLAGS", "LINENOISE_OBJ",
    )
    probe = "\n".join("$(info {}=$({}))".format(key, key) for key in names)
    probe += "\n.PHONY: platform_probe\nplatform_probe:\n\t@:\n"
    cases = [
        ("Linux", "", None, 1, False),
        ("Darwin", "", None, 1, False),
        ("MINGW64_NT-10.0", "", None, 1, True),
        ("MINGW32_NT-10.0", "", None, 1, True),
        ("MSYS_NT-10.0", "", None, 1, True),
        ("", "Windows_NT", None, 1, True),
        ("Linux", "", 1, 1, True),
        ("MINGW64_NT-10.0", "Windows_NT", 0, 1, False),
        ("Linux", "", None, 0, False),
        ("MINGW64_NT-10.0", "Windows_NT", None, 0, True),
    ]
    env = os.environ.copy()
    for key in (*names, "MAKEFLAGS", "MFLAGS", "MAKEOVERRIDES", "PORTABLE", "RELEASE"):
        env.pop(key, None)
    for uname, system, override, threads, windows in cases:
        args = ["make", "--no-print-directory", "-s", "-f", "Makefile",
                "-f", "-", "platform_probe", "UNAME_S=" + uname,
                "OS=" + system, "USE_TINYPAR=" + str(threads)]
        if override is not None:
            args.append("WINDOWS_GNU=" + str(override))
        result = subprocess.run(args, input=probe, text=True, cwd=root,
                                env=env, capture_output=True, check=True)
        values = dict(line.split("=", 1) for line in result.stdout.splitlines()
                      if "=" in line)
        backend = "serial" if not threads else "win32" if windows else "posix"
        expected = {
            "WINDOWS_GNU": str(int(windows)),
            "TINYPAR_PLATFORM_SRC": "vendor/tinypar/src/tinypar_" + backend + ".c",
            "EXEEXT": ".exe" if windows else "",
            "PICFLAGS": "" if windows else "-fPIC",
            "LUA_PLAT_FLAGS": "-DLUA_USE_WINDOWS" if windows else "-DLUA_USE_POSIX",
            "LINENOISE_OBJ": "" if windows else "build/linenoise/linenoise.o",
        }
        for key, value in expected.items():
            assert values[key] == value, (args, key, values[key], value)
        for key in ("CFLAGS", "LDFLAGS"):
            assert ("-pthread" in values[key].split()) == bool(threads and not windows), (args, key)
        assert ("-DTINYPAR_NO_THREADS" in values["CFLAGS"].split()) == (not threads), args
        tool_names = ("EXEEXT", "THREAD_FLAGS", "PROBE_LIBS")
        tool_probe = "\n".join("$(info {}=$({}))".format(key, key) for key in tool_names)
        tool_probe += "\n.PHONY: platform_probe\nplatform_probe:\n\t@:\n"
        tool_result = subprocess.run(args, input=tool_probe, text=True,
                                     cwd=root / "tools", env=env,
                                     capture_output=True, check=True)
        tool_values = dict(line.split("=", 1) for line in tool_result.stdout.splitlines()
                           if "=" in line)
        assert tool_values["EXEEXT"] == expected["EXEEXT"], args
        assert tool_values["PROBE_LIBS"] == ("-lpsapi" if windows else ""), args
        assert tool_values["THREAD_FLAGS"] == ("-pthread" if threads and not windows else ""), args
    print("GNU Make platform checks: {} cases passed".format(len(cases)))

    # Test actual Makefile flags for command-line and inherited environment
    # values. AleaTHOR passes PORTABLE=1 by default, but inherits PORTABLE=0
    # from the environment when the user requests a native build.
    flag_probe = "$(info FLAGS=$(CFLAGS))\n.PHONY: flags_probe\nflags_probe:\n\t@:\n"
    count = 0
    for directory in (root, root / "tools"):
        for release in (False, True):
            for portable, source in ((None, "unset"), ("0", "command"),
                                     ("1", "command"), ("0", "environment"),
                                     ("1", "environment")):
                args = ["make", "--no-print-directory", "-s", "-f", "Makefile",
                        "-f", "-", "flags_probe", "UNAME_S=Linux", "OS="]
                case_env = env.copy()
                if release:
                    args.append("RELEASE=1")
                if source == "command":
                    args.append("PORTABLE=" + portable)
                elif source == "environment":
                    case_env["PORTABLE"] = portable
                result = subprocess.run(args, input=flag_probe, text=True,
                                        cwd=directory, env=case_env,
                                        capture_output=True, check=True)
                flags = next(line.removeprefix("FLAGS=").split()
                             for line in result.stdout.splitlines()
                             if line.startswith("FLAGS="))
                context = (str(directory), release, portable, source, flags)
                assert ("-march=native" in flags) == (release and portable != "1"), context
                assert ("-O3" in flags) == release, context
                assert ("-DNDEBUG" in flags) == release, context
                count += 1
    print("GNU Make portability checks: {} cases passed".format(count))


if __name__ == "__main__":
    main()
