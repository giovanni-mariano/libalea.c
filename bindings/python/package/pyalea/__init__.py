# SPDX-FileCopyrightText: 2026 Giovanni MARIANO
#
# SPDX-License-Identifier: MPL-2.0

"""Native Python interface to libalea."""

from . import _alea

__version__ = _alea.version()
__all__ = sorted(name for name in dir(_alea) if not name.startswith("_"))


def __getattr__(name: str):
    """Forward public names to the native module.

    Forwarding on access keeps mutable runtime constants, such as
    ``PARALLEL_MAX_THREADS``, synchronized with the extension module.
    """
    return getattr(_alea, name)


def __dir__():
    return sorted(set(globals()) | set(dir(_alea)))
