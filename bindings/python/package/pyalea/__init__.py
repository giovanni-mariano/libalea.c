# SPDX-FileCopyrightText: 2026 Giovanni MARIANO
#
# SPDX-License-Identifier: MPL-2.0

"""Native Python interface to libalea."""

from . import _alea

__version__ = _alea.version()

# Keep the supported top-level surface explicit. __getattr__ still forwards
# extension attributes for compatibility, but adding an internal native name
# no longer publishes it accidentally through ``from pyalea import *``.
_PUBLIC_NAMES = {
    "System", "VoidResult", "XsDir", "Nuclide", "ThermalScattering",
    "NucMaterial", "Multigroup",
    "load_mcnp", "load_mcnp_string", "load_openmc", "load_openmc_string",
    "generate_void", "version", "parallel_max_threads", "set_parallel_threads",
    "get_error", "clear_error", "set_log_level", "get_log_level",
    "enable_logging", "disable_logging", "parse_zaid", "reaction_classify",
    "NODE_INVALID", "MATERIAL_NONE", "PARALLEL_BACKEND", "PARALLEL_MAX_THREADS",
    "VERSION_MAJOR", "VERSION_MINOR", "VERSION_PATCH",
    "LOG_NONE", "LOG_ERROR", "LOG_WARN", "LOG_INFO", "LOG_DEBUG", "LOG_TRACE",
    "ALEA_OP_PRIMITIVE", "ALEA_OP_UNION", "ALEA_OP_INTERSECTION",
    "ALEA_OP_DIFFERENCE", "ALEA_OP_COMPLEMENT",
    "PRIMITIVE_PLANE", "PRIMITIVE_SPHERE", "PRIMITIVE_CYLINDER_X",
    "PRIMITIVE_CYLINDER_Y", "PRIMITIVE_CYLINDER_Z", "PRIMITIVE_CONE_X",
    "PRIMITIVE_CONE_Y", "PRIMITIVE_CONE_Z", "PRIMITIVE_RPP",
    "PRIMITIVE_QUADRIC", "PRIMITIVE_TORUS_X", "PRIMITIVE_TORUS_Y",
    "PRIMITIVE_TORUS_Z", "PRIMITIVE_RCC", "PRIMITIVE_BOX", "PRIMITIVE_SPH",
    "PRIMITIVE_TRC", "PRIMITIVE_ELL", "PRIMITIVE_REC", "PRIMITIVE_WED",
    "PRIMITIVE_RHP", "PRIMITIVE_ARB",
    "ClusterContext", "cluster_initialize", "cluster_finalize", "cluster_create",
}
__all__ = sorted(name for name in _PUBLIC_NAMES if hasattr(_alea, name))


def __getattr__(name: str):
    """Forward public names to the native module.

    Forwarding on access keeps mutable runtime constants, such as
    ``PARALLEL_MAX_THREADS``, synchronized with the extension module.
    """
    return getattr(_alea, name)


def __dir__():
    return sorted(set(globals()) | set(dir(_alea)))
