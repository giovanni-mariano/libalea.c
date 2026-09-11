# SPDX-FileCopyrightText: 2026 Giovanni MARIANO
#
# SPDX-License-Identifier: MPL-2.0

"""Limit MkDocs builds to the explicitly public documentation sources."""

from pathlib import Path


PUBLIC_DOCUMENTS = frozenset(
    {
        "index.md",
        "TUTORIAL.md",
        "LUA_TUTORIAL.md",
        "CONCEPTS.md",
        "API.md",
        "ARCHITECTURE.md",
        "SURFACE_DEDUP.md",
    }
)


def on_files(files, *, config):
    """Remove non-public files originating in docs_dir, preserving theme assets."""

    docs_dir = Path(config.docs_dir).resolve()
    for file in list(files):
        if file.src_dir is None:
            continue
        if Path(file.src_dir).resolve() == docs_dir and file.src_uri not in PUBLIC_DOCUMENTS:
            files.remove(file)
    return files
