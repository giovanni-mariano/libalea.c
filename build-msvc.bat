@echo off
rem SPDX-FileCopyrightText: 2026 Giovanni MARIANO
rem
rem SPDX-License-Identifier: MPL-2.0
rem
rem Thin shim: delegate to build-msvc.ps1 (the real logic) for cmd users / CI.
rem Usage:  build-msvc.bat full   |   build-msvc.bat test   |   build-msvc.bat USE_OPENMP=1 full
powershell -NoProfile -ExecutionPolicy Bypass -File "%~dp0build-msvc.ps1" %*
