// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/* This file is part of pyalea_binding.c — do NOT compile separately.
 * It is #included from pyalea_binding.c to keep a single translation unit.
 *
 * Contents: Utilities (set_verbose, validate, print_summary, set_tolerance,
 *           clone, reset), config (get_config, set_config, set_log_level).
 */

/* ============================================================================
 * Utilities
 * ============================================================================ */

static PyObject* PyAleaSystem_set_verbose(PyAleaSystemObject* self, PyObject* args) {
    int verbose;

    if (!PyArg_ParseTuple(args, "p", &verbose)) {
        return NULL;
    }

    if (!self->sys) {
        PyErr_SetString(PyExc_RuntimeError, "System not initialized");
        return NULL;
    }

    alea_config_t cfg = alea_get_config(self->sys);
    cfg.log_level = verbose ? ALEA_LOG_LEVEL_INFO : ALEA_LOG_LEVEL_WARN;
    alea_set_config(self->sys, &cfg);
    Py_RETURN_NONE;
}

static PyObject* PyAleaSystem_validate(PyAleaSystemObject* self, PyObject* Py_UNUSED(ignored)) {
    if (!self->sys) {
        PyErr_SetString(PyExc_RuntimeError, "System not initialized");
        return NULL;
    }

    int issues = alea_validate(self->sys);
    return PyLong_FromLong(issues);
}

static PyObject* PyAleaSystem_request_interrupt(
        PyAleaSystemObject* self, PyObject* Py_UNUSED(ignored)) {
    if (!self->sys) {
        PyErr_SetString(PyExc_RuntimeError, "System not initialized");
        return NULL;
    }
    /* The native flag is process-wide and intentionally broadcasts to every
     * active cooperative operation. */
    alea_interrupt();
    Py_RETURN_NONE;
}

static PyObject* PyAleaSystem_print_summary(PyAleaSystemObject* self, PyObject* Py_UNUSED(ignored)) {
    if (!self->sys) {
        PyErr_SetString(PyExc_RuntimeError, "System not initialized");
        return NULL;
    }

    alea_print_summary(self->sys);
    Py_RETURN_NONE;
}

static PyObject* PyAleaSystem_set_tolerance(PyAleaSystemObject* self, PyObject* args, PyObject* kwds) {
    double abs_tol = 1e-6;
    double rel_tol = 1e-9;
    double zero_thresh = 1e-10;
    static char* kwlist[] = {"abs_tol", "rel_tol", "zero_thresh", NULL};

    if (!PyArg_ParseTupleAndKeywords(args, kwds, "|ddd", kwlist, &abs_tol, &rel_tol, &zero_thresh)) {
        return NULL;
    }

    if (!self->sys) {
        PyErr_SetString(PyExc_RuntimeError, "System not initialized");
        return NULL;
    }

    alea_config_t cfg = alea_get_config(self->sys);
    cfg.abs_tol = abs_tol;
    cfg.rel_tol = rel_tol;
    cfg.zero_threshold = zero_thresh;
    alea_set_config(self->sys, &cfg);
    Py_RETURN_NONE;
}

static PyObject* PyAleaSystem_clone(PyAleaSystemObject* self, PyObject* Py_UNUSED(ignored)) {
    if (!self->sys) {
        PyErr_SetString(PyExc_RuntimeError, "System not initialized");
        return NULL;
    }

    alea_system_t* cloned = alea_clone(self->sys);
    if (!cloned) {
        PyErr_SetString(PyExc_RuntimeError, "Failed to clone system");
        return NULL;
    }

    PyAleaSystemObject* obj = (PyAleaSystemObject*)PyAleaSystemType.tp_alloc(&PyAleaSystemType, 0);
    if (!obj) {
        alea_destroy(cloned);
        return NULL;
    }

    obj->sys = cloned;
    obj->owns_sys = 1;
    obj->mcnp_model = NULL;
    if (copy_mcnp_importance_sidecar(self, obj) < 0) {
        Py_DECREF(obj);
        return NULL;
    }
    return (PyObject*)obj;
}

static PyObject* PyAleaSystem_reset(PyAleaSystemObject* self, PyObject* Py_UNUSED(ignored)) {
    if (!self->sys) {
        PyErr_SetString(PyExc_RuntimeError, "System not initialized");
        return NULL;
    }

    alea_reset(self->sys);
    Py_RETURN_NONE;
}

/* ============================================================================
 * Config
 * ============================================================================ */

static PyObject* PyAleaSystem_set_log_level(PyAleaSystemObject* self, PyObject* args) {
    int level;
    if (!PyArg_ParseTuple(args, "i", &level)) return NULL;
    if (!self->sys) { PyErr_SetString(PyExc_RuntimeError, "System not initialized"); return NULL; }

    alea_config_t cfg = alea_get_config(self->sys);
    cfg.log_level = level;
    alea_set_config(self->sys, &cfg);
    Py_RETURN_NONE;
}

static PyObject* PyAleaSystem_get_config(PyAleaSystemObject* self, PyObject* Py_UNUSED(ignored)) {
    if (!self->sys) { PyErr_SetString(PyExc_RuntimeError, "System not initialized"); return NULL; }

    alea_config_t cfg = alea_get_config(self->sys);

    PyObject* dict = PyDict_New();
    dict_set_new(dict, "abs_tol", PyFloat_FromDouble(cfg.abs_tol));
    dict_set_new(dict, "rel_tol", PyFloat_FromDouble(cfg.rel_tol));
    dict_set_new(dict, "zero_threshold", PyFloat_FromDouble(cfg.zero_threshold));
    dict_set_new(dict, "dedup", PyBool_FromLong(cfg.dedup));
    dict_set_new(dict, "log_level", PyLong_FromLong(cfg.log_level));
    dict_set_new(dict, "export_materials", PyBool_FromLong(cfg.export_materials));
    dict_set_new(dict, "export_transforms", PyBool_FromLong(cfg.export_transforms));
    dict_set_new(dict, "universe_depth", PyLong_FromLong(cfg.universe_depth));
    dict_set_new(dict, "fill_depth", PyLong_FromLong(cfg.fill_depth));
    dict_set_new(dict, "void_max_depth", PyLong_FromLong(cfg.void_max_depth));
    dict_set_new(dict, "void_min_size", PyFloat_FromDouble(cfg.void_min_size));
    dict_set_new(dict, "void_probes_per_axis", PyLong_FromLong(cfg.void_probes_per_axis));
    dict_set_new(dict, "merge_cell_weight", PyFloat_FromDouble(cfg.merge_cell_weight));
    dict_set_new(dict, "merge_surface_weight", PyFloat_FromDouble(cfg.merge_surface_weight));
    dict_set_new(dict, "merge_max_surfaces", PyLong_FromLong(cfg.merge_max_surfaces));
    dict_set_new(dict, "merge_min_cells", PyLong_FromLong(cfg.merge_min_cells));
    dict_set_new(dict, "merge_use_greedy", PyBool_FromLong(cfg.merge_use_greedy));
    dict_set_new(dict, "void_consolidate", PyLong_FromLong(cfg.void_consolidate));
    dict_set_new(dict, "flatten_max_depth", PyLong_FromLong(cfg.flatten_max_depth));
    return dict;
}

static PyObject* PyAleaSystem_set_config(PyAleaSystemObject* self, PyObject* args) {
    PyObject* dict;
    if (!PyArg_ParseTuple(args, "O!", &PyDict_Type, &dict)) return NULL;
    if (!self->sys) { PyErr_SetString(PyExc_RuntimeError, "System not initialized"); return NULL; }

    alea_config_t cfg = alea_get_config(self->sys);
    PyObject* val;

#define SET_INT(field)    if ((val = PyDict_GetItemString(dict, #field))) cfg.field = (int)PyLong_AsLong(val)
#define SET_DOUBLE(field) if ((val = PyDict_GetItemString(dict, #field))) cfg.field = PyFloat_AsDouble(val)
#define SET_BOOL(field)   if ((val = PyDict_GetItemString(dict, #field))) cfg.field = PyObject_IsTrue(val)

    SET_DOUBLE(abs_tol);
    SET_DOUBLE(rel_tol);
    SET_DOUBLE(zero_threshold);
    SET_BOOL(dedup);
    SET_INT(log_level);
    SET_BOOL(export_materials);
    SET_BOOL(export_transforms);
    SET_INT(universe_depth);
    SET_INT(fill_depth);
    SET_INT(void_max_depth);
    SET_DOUBLE(void_min_size);
    SET_INT(void_probes_per_axis);
    SET_DOUBLE(merge_cell_weight);
    SET_DOUBLE(merge_surface_weight);
    SET_INT(merge_max_surfaces);
    SET_INT(merge_min_cells);
    SET_BOOL(merge_use_greedy);
    SET_INT(void_consolidate);
    SET_INT(flatten_max_depth);

#undef SET_INT
#undef SET_DOUBLE
#undef SET_BOOL

    if (PyErr_Occurred()) return NULL;

    alea_set_config(self->sys, &cfg);
    Py_RETURN_NONE;
}
