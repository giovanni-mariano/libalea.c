// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/* This file is part of pyalea_binding.c — do NOT compile separately.
 * It is #included from pyalea_binding.c to keep a single translation unit.
 *
 * Contents: VoidResult Python type (struct, dealloc, methods, type spec),
 *           generate_void module function.
 */

/* ============================================================================
 * VoidResult Python Type
 * ============================================================================ */

typedef struct {
    PyObject_HEAD
    void_result_t* result;
    PyAleaSystemObject* sys_ref;  /* Keep system alive */
} PyAleaVoidResultObject;

static void PyAleaVoidResult_dealloc(PyAleaVoidResultObject* self) {
    if (self->result) {
        alea_void_free(self->result);
    }
    Py_XDECREF(self->sys_ref);
    Py_TYPE(self)->tp_free((PyObject*)self);
}

static PyObject* PyAleaVoidResult_get_box_count(PyAleaVoidResultObject* self, void* closure) {
    (void)closure;
    if (!self->result) {
        PyErr_SetString(PyExc_RuntimeError, "Result not initialized");
        return NULL;
    }
    return PyLong_FromSize_t(alea_void_count(self->result));
}

static PyObject* PyAleaVoidResult_get_execution_stats(
        PyAleaVoidResultObject* self, void* closure) {
    (void)closure;
    if (!self->result) {
        PyErr_SetString(PyExc_RuntimeError, "Result not initialized");
        return NULL;
    }
    alea_void_execution_stats_t stats;
    if (alea_void_get_execution_stats(self->result, &stats) != 0) {
        PyErr_SetString(PyExc_RuntimeError, alea_error());
        return NULL;
    }
    return Py_BuildValue(
        "{s:K,s:K,s:K,s:K,s:K,s:K}",
        "requested_workers", (unsigned long long)stats.requested_workers,
        "actual_workers", (unsigned long long)stats.actual_workers,
        "frontier_task_count", (unsigned long long)stats.frontier_task_count,
        "parallel_batch_count", (unsigned long long)stats.parallel_batch_count,
        "scratch_bytes_per_worker",
            (unsigned long long)stats.scratch_bytes_per_worker,
        "reserved_parallel_scratch_bytes",
            (unsigned long long)stats.reserved_parallel_scratch_bytes);
}

static PyObject* PyAleaVoidResult_get_box(PyAleaVoidResultObject* self, PyObject* args) {
    size_t index;
    if (!PyArg_ParseTuple(args, "n", &index)) return NULL;

    if (!self->result) {
        PyErr_SetString(PyExc_RuntimeError, "Result not initialized");
        return NULL;
    }

    alea_bbox_t box;
    if (alea_void_get(self->result, index, &box) < 0) {
        PyErr_Format(PyExc_IndexError, "Box index %zu out of range", index);
        return NULL;
    }

    return Py_BuildValue("(dddddd)",
        box.min_x, box.max_x,
        box.min_y, box.max_y,
        box.min_z, box.max_z);
}

static PyObject* PyAleaVoidResult_get_boxes(PyAleaVoidResultObject* self, PyObject* Py_UNUSED(ignored)) {
    if (!self->result) {
        PyErr_SetString(PyExc_RuntimeError, "Result not initialized");
        return NULL;
    }

    size_t count = alea_void_count(self->result);
    PyObject* list = PyList_New(count);
    if (!list) return NULL;

    for (size_t i = 0; i < count; i++) {
        alea_bbox_t box;
        if (alea_void_get(self->result, i, &box) < 0) {
            Py_DECREF(list);
            PyErr_SetString(PyExc_RuntimeError, "Failed to get box");
            return NULL;
        }

        PyObject* bbox = Py_BuildValue("(dddddd)",
            box.min_x, box.max_x,
            box.min_y, box.max_y,
            box.min_z, box.max_z);
        if (!bbox) {
            Py_DECREF(list);
            return NULL;
        }
        PyList_SET_ITEM(list, i, bbox);
    }

    return list;
}

static PyObject* PyAleaVoidResult_to_node(PyAleaVoidResultObject* self, PyObject* Py_UNUSED(ignored)) {
    if (!self->result || !self->sys_ref || !self->sys_ref->sys) {
        PyErr_SetString(PyExc_RuntimeError, "Result or system not initialized");
        return NULL;
    }

    alea_node_id_t node = alea_void_to_node(self->sys_ref->sys, self->result);
    if (node == ALEA_NODE_ID_INVALID) {
        PyErr_SetString(PyExc_RuntimeError, "Failed to create void node");
        return NULL;
    }

    return PyLong_FromUnsignedLong(node);
}

static PyObject* PyAleaVoidResult_merge(PyAleaVoidResultObject* self, PyObject* Py_UNUSED(ignored)) {
    if (!self->result || !self->sys_ref || !self->sys_ref->sys) {
        PyErr_SetString(PyExc_RuntimeError, "Result or system not initialized");
        return NULL;
    }

    int result = alea_void_merge(self->sys_ref->sys, self->result);
    if (result < 0) {
        PyErr_SetString(PyExc_RuntimeError, alea_error());
        return NULL;
    }
    return PyLong_FromLong(result);
}

static PyObject* PyAleaVoidResult_add_cells(PyAleaVoidResultObject* self, PyObject* Py_UNUSED(ignored)) {
    if (!self->result || !self->sys_ref || !self->sys_ref->sys) {
        PyErr_SetString(PyExc_RuntimeError, "Result or system not initialized");
        return NULL;
    }

    int added = alea_void_add_cells(self->sys_ref->sys, self->result);
    if (added < 0) {
        PyErr_SetString(PyExc_RuntimeError, alea_error());
        return NULL;
    }
    return PyLong_FromLong(added);
}

static PyObject* PyAleaVoidResult_add_graveyard(PyAleaVoidResultObject* self, PyObject* Py_UNUSED(ignored)) {
    if (!self->result || !self->sys_ref || !self->sys_ref->sys) {
        PyErr_SetString(PyExc_RuntimeError, "Result or system not initialized");
        return NULL;
    }

    int result = alea_void_add_graveyard(self->sys_ref->sys, self->result);
    if (result < 0) {
        PyErr_SetString(PyExc_RuntimeError, alea_error());
        return NULL;
    }
    return PyLong_FromLong(result);
}

static PyGetSetDef PyAleaVoidResult_getsetters[] = {
    {"box_count", (getter)PyAleaVoidResult_get_box_count, NULL, "Number of void boxes", NULL},
    {"execution_stats", (getter)PyAleaVoidResult_get_execution_stats, NULL,
     "Parallel execution statistics", NULL},
    {NULL}
};

static PyMethodDef PyAleaVoidResult_methods[] = {
    {"get_box", (PyCFunction)PyAleaVoidResult_get_box, METH_VARARGS,
     "get_box(index) -> (xmin, xmax, ymin, ymax, zmin, zmax)"},
    {"get_boxes", (PyCFunction)PyAleaVoidResult_get_boxes, METH_NOARGS,
     "get_boxes() -> list of bbox tuples"},
    {"to_node", (PyCFunction)PyAleaVoidResult_to_node, METH_NOARGS,
     "to_node() -> node_id\n\nConvert void result to CSG node (union of boxes)."},
    {"merge", (PyCFunction)PyAleaVoidResult_merge, METH_NOARGS,
     "merge() -> int\n\nMerge void cells to reduce count while balancing complexity."},
    {"add_cells", (PyCFunction)PyAleaVoidResult_add_cells, METH_NOARGS,
     "add_cells() -> int\n\nAdd void boxes as cells to the system. Returns number added."},
    {"add_graveyard", (PyCFunction)PyAleaVoidResult_add_graveyard, METH_NOARGS,
     "add_graveyard() -> int\n\nAdd a graveyard cell (sphere + outside) enclosing the void bounds."},
    {NULL}
};

static PyTypeObject PyAleaVoidResultType = {
    PyVarObject_HEAD_INIT(NULL, 0)
    .tp_name = "pyalea._alea.VoidResult",
    .tp_doc = PyDoc_STR("Result of void generation."),
    .tp_basicsize = sizeof(PyAleaVoidResultObject),
    .tp_itemsize = 0,
    .tp_flags = Py_TPFLAGS_DEFAULT,
    .tp_dealloc = (destructor)PyAleaVoidResult_dealloc,
    .tp_methods = PyAleaVoidResult_methods,
    .tp_getset = PyAleaVoidResult_getsetters,
};

static PyObject* mod_generate_void(PyObject* self, PyObject* args, PyObject* kwds) {
    (void)self;
    PyAleaSystemObject* sys_obj;
    PyObject* bounds_obj = Py_None;
    PyObject* bounds_region_obj = Py_None;
    int max_depth = 8;
    double min_size = 0.1;
    int probes_per_axis = 3;
    Py_ssize_t workers = 0;
    unsigned long long max_parallel_scratch_bytes = 64ULL * 1024ULL * 1024ULL;
    static char* kwlist[] = {
        "system", "bounds", "max_depth", "min_size", "probes_per_axis",
        "bounds_region", "workers", "max_parallel_scratch_bytes", NULL
    };

    if (!PyArg_ParseTupleAndKeywords(args, kwds, "O!|OidiOnK", kwlist,
                                     &PyAleaSystemType, &sys_obj,
                                     &bounds_obj, &max_depth, &min_size, &probes_per_axis,
                                     &bounds_region_obj, &workers,
                                     &max_parallel_scratch_bytes)) {
        return NULL;
    }

    if (!sys_obj->sys) {
        PyErr_SetString(PyExc_RuntimeError, "System not initialized");
        return NULL;
    }

    if (bounds_obj != Py_None && bounds_region_obj != Py_None) {
        PyErr_SetString(PyExc_ValueError, "bounds and bounds_region are mutually exclusive");
        return NULL;
    }
    if (workers < 0) {
        PyErr_SetString(PyExc_ValueError, "workers must be non-negative");
        return NULL;
    }

    /* Parse bounds if provided */
    alea_bbox_t bounds;
    alea_bbox_t* bounds_ptr = NULL;
    alea_node_id_t bounds_region = ALEA_NODE_ID_INVALID;

    if (bounds_obj != Py_None) {
        if (!PyTuple_Check(bounds_obj) || PyTuple_Size(bounds_obj) != 6) {
            PyErr_SetString(PyExc_TypeError, "bounds must be (xmin, xmax, ymin, ymax, zmin, zmax)");
            return NULL;
        }
        bounds.min_x = PyFloat_AsDouble(PyTuple_GetItem(bounds_obj, 0));
        bounds.max_x = PyFloat_AsDouble(PyTuple_GetItem(bounds_obj, 1));
        bounds.min_y = PyFloat_AsDouble(PyTuple_GetItem(bounds_obj, 2));
        bounds.max_y = PyFloat_AsDouble(PyTuple_GetItem(bounds_obj, 3));
        bounds.min_z = PyFloat_AsDouble(PyTuple_GetItem(bounds_obj, 4));
        bounds.max_z = PyFloat_AsDouble(PyTuple_GetItem(bounds_obj, 5));
        if (PyErr_Occurred()) return NULL;
        bounds_ptr = &bounds;
    }

    if (bounds_region_obj != Py_None) {
        unsigned long node_id = PyLong_AsUnsignedLong(bounds_region_obj);
        if (PyErr_Occurred()) return NULL;
        if (node_id > (unsigned long)UINT32_MAX) {
            PyErr_SetString(PyExc_OverflowError, "bounds_region node id is out of range");
            return NULL;
        }
        bounds_region = (alea_node_id_t)node_id;
    }

    alea_void_options_t options;
    alea_void_options_init(&options);
    options.max_depth = max_depth;
    options.min_size = min_size;
    options.probes_per_axis = probes_per_axis;
    options.requested_workers = (size_t)workers;
    options.max_parallel_scratch_bytes =
        (uint64_t)max_parallel_scratch_bytes;

    void_result_t* result;
    sighandler_func old_sigint = install_sigint();
    Py_BEGIN_ALLOW_THREADS
    if (bounds_region != ALEA_NODE_ID_INVALID) {
        result = alea_void_generate_in_region_ex(
            sys_obj->sys, bounds_region, &options);
    } else {
        result = alea_void_generate_in_bbox_ex(
            sys_obj->sys, bounds_ptr, &options);
    }
    Py_END_ALLOW_THREADS

    if (restore_sigint(old_sigint)) {
        if (result) alea_void_free(result);
        return NULL;
    }

    if (!result) {
        const char* detail = alea_error();
        PyErr_SetString(PyExc_RuntimeError,
                        detail && detail[0] ? detail : "Failed to generate void");
        return NULL;
    }

    /* Create Python object */
    PyAleaVoidResultObject* obj = PyObject_New(PyAleaVoidResultObject, &PyAleaVoidResultType);
    if (!obj) {
        alea_void_free(result);
        return NULL;
    }

    obj->result = result;
    obj->sys_ref = sys_obj;
    Py_INCREF(sys_obj);

    return (PyObject*)obj;
}
