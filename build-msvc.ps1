# SPDX-FileCopyrightText: 2026 Giovanni MARIANO
#
# SPDX-License-Identifier: MPL-2.0

<#
.SYNOPSIS
    Bootstrap an MSVC x64 build of libalea through Makefile.msvc.

.DESCRIPTION
    Single entry point for building with MSVC — for CI and for contributors who
    don't already have a VS dev prompt open. It:

      1. Locates Visual Studio with vswhere, called by its ABSOLUTE path (so it
         works even when the VS Installer directory isn't on PATH).
      2. Requires the x64 C++ toolset (VC.Tools.x86.x64).
      3. Enters the x64 dev environment via vcvarsall.bat (no toolset pin — uses
         the install's default/latest toolset), then runs
         `nmake /f Makefile.msvc <args>`.

    When USE_OPENMP=1 is among the arguments, it also locates
    libomp140.x86_64.dll under the VS install and prepends its directory to PATH
    so the OpenMP-enabled test executables can load it at run time (vcvars does
    not put the LLVM OpenMP redist on PATH).

.EXAMPLE
    ./build-msvc.ps1 full
    ./build-msvc.ps1 test
    ./build-msvc.ps1 USE_OPENMP=1 RELEASE=1 full
#>
[CmdletBinding()]
param(
    [Parameter(ValueFromRemainingArguments = $true)]
    [string[]] $NmakeArgs
)

$ErrorActionPreference = 'Stop'

# --- locate vswhere by absolute path (avoids "Installer not on PATH") ----------
$vswhere = Join-Path ${env:ProgramFiles(x86)} 'Microsoft Visual Studio\Installer\vswhere.exe'
if (-not (Test-Path $vswhere)) {
    throw "vswhere not found at '$vswhere'. Install Visual Studio 2019+ with the C++ x64 toolset."
}

# --- find a VS install that has the x64 C++ toolset ---------------------------
$vsPath = & $vswhere -latest -products * `
    -requires Microsoft.VisualStudio.Component.VC.Tools.x86.x64 `
    -property installationPath
if (-not $vsPath) {
    throw "No VS install with the C++ x64 toolset (VC.Tools.x86.x64). Add it via the VS Installer."
}

$vcvars = Join-Path $vsPath 'VC\Auxiliary\Build\vcvarsall.bat'
if (-not (Test-Path $vcvars)) {
    throw "vcvarsall.bat not found at '$vcvars'."
}

$repoRoot = Split-Path -Parent $MyInvocation.MyCommand.Path
$argStr   = ($NmakeArgs -join ' ')

# --- OpenMP: make libomp140.x86_64.dll loadable by the test exes --------------
# /openmp:llvm exes load this DLL at startup, but vcvars doesn't put the LLVM
# OpenMP redist on PATH. Find it under the VS tree and prepend its directory.
$ompPrep = ''
if ($NmakeArgs -contains 'USE_OPENMP=1') {
    $roots = @(
        (Join-Path $vsPath 'VC\Redist'),
        (Join-Path $vsPath 'VC\Tools\Llvm')
    ) | Where-Object { Test-Path $_ }

    $dll = $null
    foreach ($root in $roots) {
        $dll = Get-ChildItem -Path $root -Recurse -Filter 'libomp140.x86_64.dll' -ErrorAction SilentlyContinue |
               Select-Object -First 1
        if ($dll) { break }
    }

    if ($dll) {
        Write-Host "OpenMP runtime: $($dll.FullName)"
        $ompPrep = "set `"PATH=$($dll.DirectoryName);%PATH%`" && "
    } else {
        Write-Warning "libomp140.x86_64.dll not found under '$vsPath'. OpenMP test exes may fail to start; install the VS 'C++ Clang/OpenMP' redist component."
    }
}

# --- vcvars + (optional PATH prep) + nmake, all in one cmd session -------------
$cmd = "call `"$vcvars`" x64 && cd /d `"$repoRoot`" && $ompPrep nmake /nologo /f Makefile.msvc $argStr"
Write-Host ">> $cmd"
& cmd.exe /c $cmd
exit $LASTEXITCODE
