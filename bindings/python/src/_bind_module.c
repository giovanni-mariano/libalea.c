// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/* This file is part of pyalea_binding.c — do NOT compile separately.
 * It is #included from pyalea_binding.c to keep a single translation unit.
 *
 * Contents: Module-level functions (load_mcnp, load_openmc, load_mcnp_string,
 *           version, error functions), logging integration, module definition,
 *           PyInit__alea.
 */

/* ============================================================================
 * Module-level Functions
 * ============================================================================ */

static PyObject* mod_load_mcnp(PyObject* self, PyObject* args) {
    (void)self;
    const char* filename;
    if (!PyArg_ParseTuple(args, "s", &filename)) return NULL;

    mcnp_model_t* model;
    sighandler_func old_sigint = install_sigint();
    Py_BEGIN_ALLOW_THREADS
    model = mcnp_load(filename);
    Py_END_ALLOW_THREADS
    if (restore_sigint(old_sigint)) {
        if (model) mcnp_model_destroy(model);
        return NULL;
    }

    if (!model) {
        PyErr_Format(PyExc_IOError, "Failed to load %s: %s", filename, alea_error());
        return NULL;
    }

    /* Keep the MCNP sidecar alive: it owns IMP:N/P and other MCNP-only data. */
    alea_system_t* sys = mcnp_model_system(model);

    if (!sys) {
        mcnp_model_destroy(model);
        PyErr_SetString(PyExc_RuntimeError, "MCNP model did not contain a system");
        return NULL;
    }

    PyAleaSystemObject* obj = (PyAleaSystemObject*)PyAleaSystemType.tp_alloc(&PyAleaSystemType, 0);
    if (!obj) {
        mcnp_model_destroy(model);
        return NULL;
    }

    obj->sys = sys;
    obj->owns_sys = 0;
    obj->mcnp_model = model;
    return (PyObject*)obj;
}

static PyObject* mod_load_openmc(PyObject* self, PyObject* args) {
    (void)self;
    const char* filename;
    if (!PyArg_ParseTuple(args, "s", &filename)) return NULL;

    openmc_model_t* model;
    sighandler_func old_sigint = install_sigint();
    Py_BEGIN_ALLOW_THREADS
    model = openmc_load(filename);
    Py_END_ALLOW_THREADS
    if (restore_sigint(old_sigint)) {
        if (model) openmc_model_destroy(model);
        return NULL;
    }

    if (!model) {
        PyErr_Format(PyExc_IOError, "Failed to load %s: %s", filename, alea_error());
        return NULL;
    }

    /* Take ownership of the system, then free the detached model shell. */
    alea_system_t* sys = openmc_model_take_system(model);
    openmc_model_destroy(model);

    if (!sys) {
        PyErr_SetString(PyExc_RuntimeError, "OpenMC model did not contain a system");
        return NULL;
    }

    PyAleaSystemObject* obj = (PyAleaSystemObject*)PyAleaSystemType.tp_alloc(&PyAleaSystemType, 0);
    if (!obj) {
        alea_destroy(sys);
        return NULL;
    }

    obj->sys = sys;
    obj->owns_sys = 1;
    obj->mcnp_model = NULL;
    return (PyObject*)obj;
}

static PyObject* mod_load_mcnp_string(PyObject* self, PyObject* args) {
    (void)self;
    const char* input;
    Py_ssize_t length;

    if (!PyArg_ParseTuple(args, "s#", &input, &length)) return NULL;

    mcnp_model_t* model;
    sighandler_func old_sigint = install_sigint();
    Py_BEGIN_ALLOW_THREADS
    model = mcnp_load_string(input, (size_t)length);
    Py_END_ALLOW_THREADS
    if (restore_sigint(old_sigint)) {
        if (model) mcnp_model_destroy(model);
        return NULL;
    }

    if (!model) {
        PyErr_Format(PyExc_ValueError, "Failed to parse MCNP input: %s", alea_error());
        return NULL;
    }

    /* Keep the MCNP sidecar alive: it owns IMP:N/P and other MCNP-only data. */
    alea_system_t* sys = mcnp_model_system(model);

    if (!sys) {
        mcnp_model_destroy(model);
        PyErr_SetString(PyExc_RuntimeError, "MCNP model did not contain a system");
        return NULL;
    }

    PyAleaSystemObject* obj = (PyAleaSystemObject*)PyAleaSystemType.tp_alloc(&PyAleaSystemType, 0);
    if (!obj) {
        mcnp_model_destroy(model);
        return NULL;
    }

    obj->sys = sys;
    obj->owns_sys = 0;
    obj->mcnp_model = model;
    return (PyObject*)obj;
}

static PyObject* mod_load_openmc_string(PyObject* self, PyObject* args) {
    (void)self;
    const char* input;
    Py_ssize_t length;

    if (!PyArg_ParseTuple(args, "s#", &input, &length)) return NULL;

    openmc_model_t* model;
    sighandler_func old_sigint = install_sigint();
    Py_BEGIN_ALLOW_THREADS
    model = openmc_load_string(input, (size_t)length);
    Py_END_ALLOW_THREADS
    if (restore_sigint(old_sigint)) {
        if (model) openmc_model_destroy(model);
        return NULL;
    }

    if (!model) {
        PyErr_Format(PyExc_ValueError, "Failed to parse OpenMC input: %s", alea_error());
        return NULL;
    }

    /* Take ownership of the system, then free the detached model shell. */
    alea_system_t* sys = openmc_model_take_system(model);
    openmc_model_destroy(model);

    if (!sys) {
        PyErr_SetString(PyExc_RuntimeError, "OpenMC model did not contain a system");
        return NULL;
    }

    PyAleaSystemObject* obj = (PyAleaSystemObject*)PyAleaSystemType.tp_alloc(&PyAleaSystemType, 0);
    if (!obj) {
        alea_destroy(sys);
        return NULL;
    }

    obj->sys = sys;
    obj->owns_sys = 1;
    obj->mcnp_model = NULL;
    return (PyObject*)obj;
}

static PyObject* mod_version(PyObject* self, PyObject* Py_UNUSED(ignored)) {
    (void)self;
    return PyUnicode_FromString(alea_version());
}

static PyObject* mod_parallel_max_threads(
        PyObject* self, PyObject* Py_UNUSED(ignored)) {
    (void)self;
    return PyLong_FromLong(alea_parallel_max_threads());
}

static PyObject* mod_set_parallel_threads(PyObject* self, PyObject* args) {
    int threads;
    if (!PyArg_ParseTuple(args, "i", &threads)) return NULL;
    if (threads < 0) {
        PyErr_SetString(PyExc_ValueError, "threads must be non-negative");
        return NULL;
    }
    if (alea_parallel_set_threads(threads) != 0) {
        PyErr_SetString(
            PyExc_RuntimeError,
            "parallel worker count must be set before the first parallel operation");
        return NULL;
    }
    const int configured = alea_parallel_max_threads();
    PyObject* value = PyLong_FromLong(configured);
    if (!value) return NULL;
    if (PyObject_SetAttrString(self, "PARALLEL_MAX_THREADS", value) < 0) {
        Py_DECREF(value);
        return NULL;
    }
    return value;
}

static PyObject* mod_get_error(PyObject* self, PyObject* Py_UNUSED(ignored)) {
    (void)self;
    const char* err = alea_error();
    if (err && *err) return PyUnicode_FromString(err);
    Py_RETURN_NONE;
}

static PyObject* mod_clear_error(PyObject* self, PyObject* Py_UNUSED(ignored)) {
    (void)self;
    alea_error_clear();
    Py_RETURN_NONE;
}

/* ============================================================================
 * Logging Integration with Python's logging module
 * ============================================================================ */

/* Cached Python logging module and logger */
static PyObject* g_logging_module = NULL;
static PyObject* g_logger = NULL;

/* Map C log levels to Python logging method names */
static const char* level_to_method_name(alea_log_level_t level) {
    switch (level) {
        case ALEA_LOG_LEVEL_ERROR: return "error";
        case ALEA_LOG_LEVEL_WARN:  return "warning";
        case ALEA_LOG_LEVEL_INFO:  return "info";
        case ALEA_LOG_LEVEL_DEBUG: return "debug";
        case ALEA_LOG_LEVEL_TRACE: return "debug";  /* Python has no TRACE, use DEBUG */
        default:            return "debug";
    }
}

/* Callback that routes C library logs to Python logging module */
static void python_log_callback(alea_log_level_t level, const char* file,
                                 int line, const char* message, void* user_data) {
    (void)user_data;

    /* Must hold GIL to call Python */
    PyGILState_STATE gstate = PyGILState_Ensure();

    if (g_logger) {
        /* Build log message with file:line for debug levels */
        PyObject* log_msg;
        if (file && (level >= ALEA_LOG_LEVEL_DEBUG)) {
            log_msg = PyUnicode_FromFormat("[%s:%d] %s", file, line, message);
        } else {
            log_msg = PyUnicode_FromString(message);
        }

        if (log_msg) {
            const char* method = level_to_method_name(level);
            PyObject* result = PyObject_CallMethod(g_logger, method, "O", log_msg);
            Py_XDECREF(result);
            Py_DECREF(log_msg);
        }

        /* Clear any Python exception (don't let logging errors propagate) */
        if (PyErr_Occurred()) {
            PyErr_Clear();
        }
    }

    PyGILState_Release(gstate);
}

/* Initialize logging integration */
static int init_logging(void) {
    /* Import logging module */
    g_logging_module = PyImport_ImportModule("logging");
    if (!g_logging_module) {
        PyErr_Clear();
        return -1;
    }

    /* Use a logger owned by the binding package. */
    g_logger = PyObject_CallMethod(g_logging_module, "getLogger", "s", "pyalea");
    if (!g_logger) {
        PyErr_Clear();
        Py_CLEAR(g_logging_module);
        return -1;
    }

    /* Gate at the C level so sub-threshold messages never pay the cost of the
     * Python callback (GIL re-acquire + string alloc + logging call). Default
     * to WARN to match Python's default logger; mod_set_log_level keeps the C
     * level in sync when callers change it. */
    alea_log_set_level(ALEA_LOG_LEVEL_WARN);

    /* Set the C library callback to route to Python */
    alea_log_set_callback(python_log_callback, NULL);

    return 0;
}

/* Cleanup logging on module unload */
static void cleanup_logging(void) {
    alea_log_set_callback(NULL, NULL);
    Py_CLEAR(g_logger);
    Py_CLEAR(g_logging_module);
}

/* Map our log levels to Python logging levels */
static int to_python_log_level(int level) {
    switch (level) {
        case ALEA_LOG_LEVEL_NONE:  return 100;  /* Higher than CRITICAL, effectively disabled */
        case ALEA_LOG_LEVEL_ERROR: return 40;   /* logging.ERROR */
        case ALEA_LOG_LEVEL_WARN:  return 30;   /* logging.WARNING */
        case ALEA_LOG_LEVEL_INFO:  return 20;   /* logging.INFO */
        case ALEA_LOG_LEVEL_DEBUG: return 10;   /* logging.DEBUG */
        case ALEA_LOG_LEVEL_TRACE: return 5;    /* Below DEBUG */
        default:            return 30;   /* Default to WARNING */
    }
}

/* Map Python logging level back to our levels */
static int from_python_log_level(int py_level) {
    if (py_level >= 100) return ALEA_LOG_LEVEL_NONE;
    if (py_level >= 40) return ALEA_LOG_LEVEL_ERROR;
    if (py_level >= 30) return ALEA_LOG_LEVEL_WARN;
    if (py_level >= 20) return ALEA_LOG_LEVEL_INFO;
    if (py_level >= 10) return ALEA_LOG_LEVEL_DEBUG;
    return ALEA_LOG_LEVEL_TRACE;
}

/* Python API: set_log_level(level) - sets the Python logger level */
static PyObject* mod_set_log_level(PyObject* self, PyObject* args) {
    (void)self;
    int level;
    if (!PyArg_ParseTuple(args, "i", &level)) return NULL;

    if (level < ALEA_LOG_LEVEL_NONE || level > ALEA_LOG_LEVEL_TRACE) {
        PyErr_SetString(PyExc_ValueError, "Invalid log level (must be 0-5)");
        return NULL;
    }

    /* Gate at the C level too: messages below this never invoke the Python
     * callback, which is essential for performance on large models (the
     * callback re-acquires the GIL and calls into the logging module per
     * message). Keep the C level and the Python logger level in lockstep. */
    alea_log_set_level((alea_log_level_t)level);

    if (g_logger) {
        int py_level = to_python_log_level(level);
        PyObject* result = PyObject_CallMethod(g_logger, "setLevel", "i", py_level);
        Py_XDECREF(result);
        if (PyErr_Occurred()) {
            return NULL;
        }
    }

    Py_RETURN_NONE;
}

/* Python API: get_log_level() -> int */
static PyObject* mod_get_log_level(PyObject* self, PyObject* Py_UNUSED(ignored)) {
    (void)self;

    if (g_logger) {
        PyObject* effective_level = PyObject_CallMethod(g_logger, "getEffectiveLevel", NULL);
        if (effective_level && PyLong_Check(effective_level)) {
            int py_level = (int)PyLong_AsLong(effective_level);
            Py_DECREF(effective_level);
            return PyLong_FromLong(from_python_log_level(py_level));
        }
        Py_XDECREF(effective_level);
        PyErr_Clear();
    }

    return PyLong_FromLong(ALEA_LOG_LEVEL_WARN);  /* Default */
}

/* Python API: disable_logging() - restore default stderr output */
static PyObject* mod_disable_logging(PyObject* self, PyObject* Py_UNUSED(ignored)) {
    (void)self;
    alea_log_set_callback(NULL, NULL);
    Py_RETURN_NONE;
}

/* Python API: enable_logging() - route to Python logging module */
static PyObject* mod_enable_logging(PyObject* self, PyObject* Py_UNUSED(ignored)) {
    (void)self;
    if (!g_logger) {
        if (init_logging() < 0) {
            PyErr_SetString(PyExc_RuntimeError, "Failed to initialize Python logging integration");
            return NULL;
        }
    }
    alea_log_set_callback(python_log_callback, NULL);
    Py_RETURN_NONE;
}

/* ============================================================================
 * Module Definition
 * ============================================================================ */

static PyMethodDef mod_methods[] = {
#ifdef PYALEA_USE_MPI
    {"cluster_initialize", mod_cluster_initialize, METH_NOARGS, NULL},
    {"cluster_finalize", mod_cluster_finalize, METH_NOARGS, NULL},
    {"cluster_create", mod_cluster_create, METH_O, NULL},
#endif
    {"load_mcnp", mod_load_mcnp, METH_VARARGS,
     "load_mcnp(filename) -> System\n\nLoad MCNP input file."},
    {"load_mcnp_string", mod_load_mcnp_string, METH_VARARGS,
     "load_mcnp_string(input) -> System\n\nLoad MCNP from string."},
    {"load_openmc", mod_load_openmc, METH_VARARGS,
     "load_openmc(filename) -> System\n\nLoad OpenMC XML geometry file."},
    {"load_openmc_string", mod_load_openmc_string, METH_VARARGS,
     "load_openmc_string(input) -> System\n\nLoad OpenMC from XML string."},
    {"generate_void", (PyCFunction)mod_generate_void, METH_VARARGS | METH_KEYWORDS,
     "generate_void(system, bounds=None, max_depth=8, min_size=0.1, probes_per_axis=3, bounds_region=None, workers=0, max_parallel_scratch_bytes=67108864) -> VoidResult\n\n"
     "Generate void regions using octree algorithm."},
    {"version", mod_version, METH_NOARGS,
     "version() -> str\n\nGet library version string."},
    {"parallel_max_threads", mod_parallel_max_threads, METH_NOARGS,
     "parallel_max_threads() -> int\n\nGet the process parallel worker limit."},
    {"set_parallel_threads", mod_set_parallel_threads, METH_VARARGS,
     "set_parallel_threads(threads) -> int\n\n"
     "Set the persistent process worker limit before the first parallel operation; "
     "zero uses hardware detection."},
    {"get_error", mod_get_error, METH_NOARGS,
     "get_error() -> str or None\n\nGet last error message."},
    {"clear_error", mod_clear_error, METH_NOARGS,
     "clear_error()\n\nClear error state."},
    /* Logging */
    {"set_log_level", mod_set_log_level, METH_VARARGS,
     "set_log_level(level)\n\n"
     "Set log level for the 'pyalea' logger.\n\n"
     "Levels: LOG_NONE=0, LOG_ERROR=1, LOG_WARN=2, LOG_INFO=3, LOG_DEBUG=4, LOG_TRACE=5\n\n"
     "Example:\n"
     "    import pyalea\n"
     "    pyalea.set_log_level(pyalea.LOG_DEBUG)  # See debug messages"},
    {"get_log_level", mod_get_log_level, METH_NOARGS,
     "get_log_level() -> int\n\nGet current log level of the 'pyalea' logger."},
    {"enable_logging", mod_enable_logging, METH_NOARGS,
     "enable_logging()\n\n"
     "Route C library logs to Python's logging module (logger: 'pyalea').\n"
     "This is enabled by default on module import."},
    {"disable_logging", mod_disable_logging, METH_NOARGS,
     "disable_logging()\n\n"
     "Disable Python logging integration, restore default stderr output."},
    /* Nuclear data */
    {"parse_zaid", mod_parse_zaid, METH_VARARGS,
     "parse_zaid(zaid) -> dict\n\n"
     "Parse ZAID string (e.g. '92235.80c') into Z, A, metastable, type."},
    {"reaction_classify", mod_reaction_classify, METH_VARARGS,
     "reaction_classify(mt) -> str\n\n"
     "Classify reaction MT into 'absorption', 'scatter', or 'multiply'."},
    {NULL, NULL, 0, NULL}
};

/* Module cleanup function */
static void mod_module_free(void* module) {
    (void)module;
    cleanup_logging();
}

static struct PyModuleDef mod_module = {
    PyModuleDef_HEAD_INIT,
    .m_name = "_alea",
    .m_doc = PyDoc_STR("pyAlea: Python bindings for libalea.\n\n"
                       "This module provides native bindings for libalea,\n"
                       "which parses MCNP input files and provides CSG tree operations.\n\n"
                       "Logging is automatically routed to Python's logging module (logger: 'pyalea')."),
    .m_size = 0,  /* Use 0 instead of -1 to enable m_free */
    .m_methods = mod_methods,
    .m_free = mod_module_free,
};

PyMODINIT_FUNC PyInit__alea(void) {
    import_array();
#ifdef PYALEA_USE_MPI
    if (PyType_Ready(&PyAleaClusterType) < 0) return NULL;
#endif
    /* Initialize types */
    if (PyType_Ready(&PyAleaSystemType) < 0) return NULL;
    if (PyType_Ready(&PyAleaVoidResultType) < 0) return NULL;
    if (PyType_Ready(&PyAleaXsDirType) < 0) return NULL;
    if (PyType_Ready(&PyAleaNuclideType) < 0) return NULL;
    if (PyType_Ready(&PyAleaThermalType) < 0) return NULL;
    if (PyType_Ready(&PyAleaNucMaterialType) < 0) return NULL;
    if (PyType_Ready(&PyAleaMultigroupType) < 0) return NULL;

    /* Create module */
    PyObject* m = PyModule_Create(&mod_module);
    if (!m) return NULL;
#ifdef PYALEA_USE_MPI
    Py_INCREF(&PyAleaClusterType);
    if (PyModule_AddObject(m, "ClusterContext", (PyObject*)&PyAleaClusterType) < 0) {
        Py_DECREF(&PyAleaClusterType);
        Py_DECREF(m);
        return NULL;
    }
#endif

    /* Add types */
    Py_INCREF(&PyAleaSystemType);
    if (PyModule_AddObject(m, "System", (PyObject*)&PyAleaSystemType) < 0) {
        Py_DECREF(&PyAleaSystemType);
        Py_DECREF(m);
        return NULL;
    }

    Py_INCREF(&PyAleaVoidResultType);
    if (PyModule_AddObject(m, "VoidResult", (PyObject*)&PyAleaVoidResultType) < 0) {
        Py_DECREF(&PyAleaVoidResultType);
        Py_DECREF(m);
        return NULL;
    }

    /* Nuclear data types */
    Py_INCREF(&PyAleaXsDirType);
    if (PyModule_AddObject(m, "XsDir", (PyObject*)&PyAleaXsDirType) < 0) {
        Py_DECREF(&PyAleaXsDirType);
        Py_DECREF(m);
        return NULL;
    }

    Py_INCREF(&PyAleaNuclideType);
    if (PyModule_AddObject(m, "Nuclide", (PyObject*)&PyAleaNuclideType) < 0) {
        Py_DECREF(&PyAleaNuclideType);
        Py_DECREF(m);
        return NULL;
    }

    Py_INCREF(&PyAleaThermalType);
    if (PyModule_AddObject(m, "ThermalScattering", (PyObject*)&PyAleaThermalType) < 0) {
        Py_DECREF(&PyAleaThermalType);
        Py_DECREF(m);
        return NULL;
    }

    Py_INCREF(&PyAleaNucMaterialType);
    if (PyModule_AddObject(m, "NucMaterial", (PyObject*)&PyAleaNucMaterialType) < 0) {
        Py_DECREF(&PyAleaNucMaterialType);
        Py_DECREF(m);
        return NULL;
    }

    Py_INCREF(&PyAleaMultigroupType);
    if (PyModule_AddObject(m, "Multigroup", (PyObject*)&PyAleaMultigroupType) < 0) {
        Py_DECREF(&PyAleaMultigroupType);
        Py_DECREF(m);
        return NULL;
    }

    /* Add constants */
    PyModule_AddIntConstant(m, "NODE_INVALID", ALEA_NODE_ID_INVALID);
    PyModule_AddIntConstant(m, "MATERIAL_NONE", ALEA_MATERIAL_NONE);

    /* Add log level constants */
    PyModule_AddIntConstant(m, "LOG_NONE", ALEA_LOG_LEVEL_NONE);
    PyModule_AddIntConstant(m, "LOG_ERROR", ALEA_LOG_LEVEL_ERROR);
    PyModule_AddIntConstant(m, "LOG_WARN", ALEA_LOG_LEVEL_WARN);
    PyModule_AddIntConstant(m, "LOG_INFO", ALEA_LOG_LEVEL_INFO);
    PyModule_AddIntConstant(m, "LOG_DEBUG", ALEA_LOG_LEVEL_DEBUG);
    PyModule_AddIntConstant(m, "LOG_TRACE", ALEA_LOG_LEVEL_TRACE);

    /* Add CSG operation types */
    PyModule_AddIntConstant(m, "ALEA_OP_PRIMITIVE", ALEA_OP_PRIMITIVE);
    PyModule_AddIntConstant(m, "ALEA_OP_UNION", ALEA_OP_UNION);
    PyModule_AddIntConstant(m, "ALEA_OP_INTERSECTION", ALEA_OP_INTERSECTION);
    PyModule_AddIntConstant(m, "ALEA_OP_DIFFERENCE", ALEA_OP_DIFFERENCE);
    PyModule_AddIntConstant(m, "ALEA_OP_COMPLEMENT", ALEA_OP_COMPLEMENT);

    /* Add primitive types */
    PyModule_AddIntConstant(m, "PRIMITIVE_PLANE", ALEA_PRIMITIVE_PLANE);
    PyModule_AddIntConstant(m, "PRIMITIVE_SPHERE", ALEA_PRIMITIVE_SPHERE);
    PyModule_AddIntConstant(m, "PRIMITIVE_CYLINDER_X", ALEA_PRIMITIVE_CYLINDER_X);
    PyModule_AddIntConstant(m, "PRIMITIVE_CYLINDER_Y", ALEA_PRIMITIVE_CYLINDER_Y);
    PyModule_AddIntConstant(m, "PRIMITIVE_CYLINDER_Z", ALEA_PRIMITIVE_CYLINDER_Z);
    PyModule_AddIntConstant(m, "PRIMITIVE_CONE_X", ALEA_PRIMITIVE_CONE_X);
    PyModule_AddIntConstant(m, "PRIMITIVE_CONE_Y", ALEA_PRIMITIVE_CONE_Y);
    PyModule_AddIntConstant(m, "PRIMITIVE_CONE_Z", ALEA_PRIMITIVE_CONE_Z);
    PyModule_AddIntConstant(m, "PRIMITIVE_RPP", ALEA_PRIMITIVE_RPP);
    PyModule_AddIntConstant(m, "PRIMITIVE_QUADRIC", ALEA_PRIMITIVE_QUADRIC);
    PyModule_AddIntConstant(m, "PRIMITIVE_TORUS_X", ALEA_PRIMITIVE_TORUS_X);
    PyModule_AddIntConstant(m, "PRIMITIVE_TORUS_Y", ALEA_PRIMITIVE_TORUS_Y);
    PyModule_AddIntConstant(m, "PRIMITIVE_TORUS_Z", ALEA_PRIMITIVE_TORUS_Z);
    PyModule_AddIntConstant(m, "PRIMITIVE_RCC", ALEA_PRIMITIVE_RCC);
    PyModule_AddIntConstant(m, "PRIMITIVE_BOX", ALEA_PRIMITIVE_BOX);
    PyModule_AddIntConstant(m, "PRIMITIVE_SPH", ALEA_PRIMITIVE_SPH);
    PyModule_AddIntConstant(m, "PRIMITIVE_TRC", ALEA_PRIMITIVE_TRC);
    PyModule_AddIntConstant(m, "PRIMITIVE_ELL", ALEA_PRIMITIVE_ELL);
    PyModule_AddIntConstant(m, "PRIMITIVE_REC", ALEA_PRIMITIVE_REC);
    PyModule_AddIntConstant(m, "PRIMITIVE_WED", ALEA_PRIMITIVE_WED);
    PyModule_AddIntConstant(m, "PRIMITIVE_RHP", ALEA_PRIMITIVE_RHP);
    PyModule_AddIntConstant(m, "PRIMITIVE_ARB", ALEA_PRIMITIVE_ARB);

    /* Add version info */
    PyModule_AddIntConstant(m, "VERSION_MAJOR", ALEA_VERSION_MAJOR);
    PyModule_AddIntConstant(m, "VERSION_MINOR", ALEA_VERSION_MINOR);
    PyModule_AddIntConstant(m, "VERSION_PATCH", ALEA_VERSION_PATCH);
    PyModule_AddStringConstant(m, "PARALLEL_BACKEND", "tinypar");
    int parallel_max_threads = alea_parallel_max_threads();
    PyModule_AddIntConstant(m, "PARALLEL_MAX_THREADS", parallel_max_threads);

    /* Initialize logging integration (route C logs to Python logging module) */
    init_logging();

    return m;
}
