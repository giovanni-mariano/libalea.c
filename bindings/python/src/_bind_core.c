// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/* This file is part of pyalea_binding.c — do NOT compile separately.
 * It is #included from pyalea_binding.c to keep a single translation unit.
 *
 * Contents: Type lifecycle (new/init/dealloc), property getters,
 *           point queries, cell operations, universe operations.
 */

/* ============================================================================
 * Type Lifecycle
 * ============================================================================ */

static void PyAleaSystem_dealloc(PyAleaSystemObject* self) {
    int model_owns_sys = self->mcnp_model && self->mcnp_model->owns_sys;
    if (self->mcnp_model)
        mcnp_model_destroy(self->mcnp_model);
    if (self->sys && self->owns_sys && !model_owns_sys) {
        alea_destroy(self->sys);
    }
    Py_TYPE(self)->tp_free((PyObject*)self);
}

static PyObject* PyAleaSystem_new(PyTypeObject* type, PyObject* args, PyObject* kwds) {
    (void)args; (void)kwds;
    PyAleaSystemObject* self = (PyAleaSystemObject*)type->tp_alloc(type, 0);
    if (self) {
        self->sys = NULL;
        self->owns_sys = 1;
        self->mcnp_model = NULL;
    }
    return (PyObject*)self;
}

static int PyAleaSystem_init(PyAleaSystemObject* self, PyObject* args, PyObject* kwds) {
    static char* kwlist[] = {NULL};
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "", kwlist)) {
        return -1;
    }

    self->sys = alea_create();
    if (!self->sys) {
        PyErr_SetString(PyExc_RuntimeError, "Failed to create CSG system");
        return -1;
    }
    self->owns_sys = 1;
    return 0;
}

/* ============================================================================
 * Property Getters
 * ============================================================================ */

static PyObject* PyAleaSystem_get_cell_count(PyAleaSystemObject* self, void* closure) {
    (void)closure;
    if (!self->sys) {
        PyErr_SetString(PyExc_RuntimeError, "System not initialized");
        return NULL;
    }
    return PyLong_FromSize_t(alea_cell_count(self->sys));
}

static PyObject* PyAleaSystem_get_surface_count(PyAleaSystemObject* self, void* closure) {
    (void)closure;
    if (!self->sys) {
        PyErr_SetString(PyExc_RuntimeError, "System not initialized");
        return NULL;
    }
    return PyLong_FromSize_t(alea_surface_count(self->sys));
}

static PyObject* PyAleaSystem_get_universe_count(PyAleaSystemObject* self, void* closure) {
    (void)closure;
    if (!self->sys) {
        PyErr_SetString(PyExc_RuntimeError, "System not initialized");
        return NULL;
    }
    return PyLong_FromSize_t(alea_universe_count(self->sys));
}

static PyObject* PyAleaSystem_get_volume_path_count(PyAleaSystemObject* self, void* closure) {
    (void)closure;
    if (!self->sys) {
        return PyLong_FromLong(0);
    }
    if (ensure_query_acceleration(self) < 0) return NULL;
    return PyLong_FromSize_t(alea_volume_path_count(self->sys));
}

/* ============================================================================
 * Point Queries
 * ============================================================================ */

static PyObject* PyAleaSystem_find_cell(PyAleaSystemObject* self, PyObject* args) {
    double x, y, z;

    if (!PyArg_ParseTuple(args, "ddd", &x, &y, &z)) {
        return NULL;
    }

    if (!self->sys) {
        PyErr_SetString(PyExc_RuntimeError, "System not initialized");
        return NULL;
    }
    if (ensure_query_acceleration(self) < 0) return NULL;

    /* Deepest hit along the first containing (DFS) chain: identical answer
     * to alea_find_cell, but the hit carries resolution_flags (undefined
     * fill regions resolve to their container, flagged). */
    alea_cell_hit_t hits[32];
    int n = alea_find_all_cells(self->sys, x, y, z, hits, 32);
    if (n <= 0) {
        Py_RETURN_NONE;  /* Point in void or outside */
    }
    int ti = 0;
    while (ti + 1 < n && hits[ti + 1].depth == hits[ti].depth + 1) ti++;
    return Py_BuildValue("(iii)", hits[ti].cell_index, hits[ti].material_id,
                         (int)hits[ti].resolution_flags);
}

static PyObject* PyAleaSystem_point_inside(PyAleaSystemObject* self, PyObject* args) {
    unsigned long node_id;
    double x, y, z;

    if (!PyArg_ParseTuple(args, "kddd", &node_id, &x, &y, &z)) {
        return NULL;
    }

    if (!self->sys) {
        PyErr_SetString(PyExc_RuntimeError, "System not initialized");
        return NULL;
    }

    bool inside = alea_point_inside(self->sys, (alea_node_id_t)node_id, x, y, z);
    return PyBool_FromLong(inside);
}

static PyObject* PyAleaSystem_material_at(PyAleaSystemObject* self, PyObject* args) {
    double x, y, z;

    if (!PyArg_ParseTuple(args, "ddd", &x, &y, &z)) {
        return NULL;
    }

    if (!self->sys) {
        PyErr_SetString(PyExc_RuntimeError, "System not initialized");
        return NULL;
    }
    if (ensure_query_acceleration(self) < 0) return NULL;

    alea_material_id_t mat = alea_material_at(self->sys, x, y, z);
    if (mat == ALEA_MATERIAL_NONE) {
        Py_RETURN_NONE;
    }
    return PyLong_FromUnsignedLong(mat);
}

/* ============================================================================
 * Cell Operations
 * ============================================================================ */

static PyObject* PyAleaSystem_get_cell(PyAleaSystemObject* self, PyObject* args) {
    int cell_id;

    if (!PyArg_ParseTuple(args, "i", &cell_id)) {
        return NULL;
    }

    if (!self->sys) {
        PyErr_SetString(PyExc_RuntimeError, "System not initialized");
        return NULL;
    }

    alea_cell_info_t info;
    if (alea_cell_find_info(self->sys, cell_id, &info) < 0) {
        PyErr_Format(PyExc_KeyError, "Cell %d not found", cell_id);
        return NULL;
    }

    PyObject* dict = Py_BuildValue("{s:i, s:i, s:d, s:N, s:i, s:i, s:i, s:k, s:(dddddd)}",
        "cell_id", info.cell_id,
        "material_id", info.material_id,
        "density", info.density,
        "is_mass_density", PyBool_FromLong(info.is_mass_density),
        "universe_id", info.universe_id,
        "fill_universe", info.fill_universe,
        "fill_transform", info.fill_transform,
        "root_node", (unsigned long)info.root,
        "bbox", info.bbox.min_x, info.bbox.max_x,
                info.bbox.min_y, info.bbox.max_y,
                info.bbox.min_z, info.bbox.max_z);
    if (!dict) return NULL;

    int cell_index = alea_cell_find(self->sys, cell_id);
    if (add_importance_fields(self, dict, (size_t)cell_index) < 0) {
        Py_DECREF(dict);
        return NULL;
    }

    if (add_lattice_fields(dict, &info) < 0) {
        Py_DECREF(dict);
        return NULL;
    }
    if (add_comment_fields(dict, &info) < 0) {
        Py_DECREF(dict);
        return NULL;
    }
    if (add_temperature_fields(dict, &info) < 0) {
        Py_DECREF(dict);
        return NULL;
    }
    return dict;
}

static PyObject* PyAleaSystem_get_cell_by_index(PyAleaSystemObject* self, PyObject* args) {
    size_t index;

    if (!PyArg_ParseTuple(args, "n", &index)) {
        return NULL;
    }

    if (!self->sys) {
        PyErr_SetString(PyExc_RuntimeError, "System not initialized");
        return NULL;
    }

    alea_cell_info_t info;
    if (alea_cell_get_info(self->sys, index, &info) < 0) {
        PyErr_Format(PyExc_IndexError, "Cell index %zu out of range", index);
        return NULL;
    }

    PyObject* dict = Py_BuildValue("{s:i, s:i, s:d, s:N, s:i, s:i, s:i, s:k, s:(dddddd)}",
        "cell_id", info.cell_id,
        "material_id", info.material_id,
        "density", info.density,
        "is_mass_density", PyBool_FromLong(info.is_mass_density),
        "universe_id", info.universe_id,
        "fill_universe", info.fill_universe,
        "fill_transform", info.fill_transform,
        "root_node", (unsigned long)info.root,
        "bbox", info.bbox.min_x, info.bbox.max_x,
                info.bbox.min_y, info.bbox.max_y,
                info.bbox.min_z, info.bbox.max_z);
    if (!dict) return NULL;

    if (add_importance_fields(self, dict, index) < 0) {
        Py_DECREF(dict);
        return NULL;
    }

    if (add_lattice_fields(dict, &info) < 0) {
        Py_DECREF(dict);
        return NULL;
    }
    if (add_comment_fields(dict, &info) < 0) {
        Py_DECREF(dict);
        return NULL;
    }
    if (add_temperature_fields(dict, &info) < 0) {
        Py_DECREF(dict);
        return NULL;
    }
    return dict;
}

static PyObject* PyAleaSystem_get_cells(PyAleaSystemObject* self, PyObject* Py_UNUSED(ignored)) {
    if (!self->sys) {
        PyErr_SetString(PyExc_RuntimeError, "System not initialized");
        return NULL;
    }

    size_t count = alea_cell_count(self->sys);
    PyObject* list = PyList_New(count);
    if (!list) return NULL;

    for (size_t i = 0; i < count; i++) {
        alea_cell_info_t info;
        if (alea_cell_get_info(self->sys, i, &info) < 0) {
            Py_DECREF(list);
            PyErr_SetString(PyExc_RuntimeError, "Failed to get cell info");
            return NULL;
        }

        PyObject* cell_dict = Py_BuildValue("{s:i, s:i, s:d, s:N, s:i, s:i, s:k}",
            "cell_id", info.cell_id,
            "material_id", info.material_id,
            "density", info.density,
            "is_mass_density", PyBool_FromLong(info.is_mass_density),
            "universe_id", info.universe_id,
            "fill_universe", info.fill_universe,
            "root_node", (unsigned long)info.root);

        if (!cell_dict) {
            Py_DECREF(list);
            return NULL;
        }

        if (add_importance_fields(self, cell_dict, i) < 0) {
            Py_DECREF(cell_dict);
            Py_DECREF(list);
            return NULL;
        }

        if (add_temperature_fields(cell_dict, &info) < 0) {
            Py_DECREF(cell_dict);
            Py_DECREF(list);
            return NULL;
        }

        if (add_lattice_fields(cell_dict, &info) < 0) {
            Py_DECREF(cell_dict);
            Py_DECREF(list);
            return NULL;
        }
        if (add_comment_fields(cell_dict, &info) < 0) {
            Py_DECREF(cell_dict);
            Py_DECREF(list);
            return NULL;
        }

        PyList_SET_ITEM(list, i, cell_dict);
    }

    return list;
}

/* Gather only categorical colour IDs. This deliberately avoids constructing
 * one rich Python cell dictionary per entity during model loading. */
static int color_id_add(PyObject* values, int value) {
    PyObject* identifier = PyLong_FromLong(value);
    if (!identifier) return -1;
    int result = PySet_Add(values, identifier);
    Py_DECREF(identifier);
    return result;
}

static PyObject* PyAleaSystem_color_identifiers(
    PyAleaSystemObject* self, PyObject* Py_UNUSED(ignored)) {
    if (!self->sys) {
        PyErr_SetString(PyExc_RuntimeError, "System not initialized");
        return NULL;
    }

    PyObject* cells = PySet_New(NULL);
    PyObject* materials = PySet_New(NULL);
    PyObject* universes = PySet_New(NULL);
    PyObject* fills = PySet_New(NULL);
    PyObject* result = NULL;
    if (!cells || !materials || !universes || !fills) goto error;

    size_t count = alea_cell_count(self->sys);
    for (size_t i = 0; i < count; i++) {
        alea_cell_info_t info;
        if (alea_cell_get_info(self->sys, i, &info) < 0) {
            PyErr_SetString(PyExc_RuntimeError, "Failed to scan cell identifiers");
            goto error;
        }
        if (color_id_add(cells, info.cell_id) < 0 ||
            color_id_add(materials, info.material_id) < 0 ||
            color_id_add(universes, info.universe_id) < 0) goto error;
        if (info.fill_universe > 0 && color_id_add(fills, info.fill_universe) < 0)
            goto error;
        for (size_t j = 0; j < info.lat_fill_count; j++) {
            if (info.lat_fill[j] > 0 && color_id_add(fills, info.lat_fill[j]) < 0)
                goto error;
        }
    }
    count = alea_material_count(self->sys);
    for (size_t i = 0; i < count; i++) {
        int material_id = alea_material_get_id(self->sys, (int)i);
        if (material_id < 0) {
            PyErr_SetString(PyExc_RuntimeError, "Failed to scan material identifiers");
            goto error;
        }
        if (color_id_add(materials, material_id) < 0) goto error;
    }

    result = PyDict_New();
    if (!result ||
        PyDict_SetItemString(result, "cell", cells) < 0 ||
        PyDict_SetItemString(result, "material", materials) < 0 ||
        PyDict_SetItemString(result, "universe", universes) < 0 ||
        PyDict_SetItemString(result, "fill", fills) < 0) goto error;
    Py_DECREF(cells);
    Py_DECREF(materials);
    Py_DECREF(universes);
    Py_DECREF(fills);
    return result;

error:
    Py_XDECREF(result);
    Py_XDECREF(cells);
    Py_XDECREF(materials);
    Py_XDECREF(universes);
    Py_XDECREF(fills);
    return NULL;
}

/* ============================================================================
 * Universe Operations
 * ============================================================================ */

static PyObject* PyAleaSystem_build_universe_index(PyAleaSystemObject* self, PyObject* Py_UNUSED(ignored)) {
    if (!self->sys) {
        PyErr_SetString(PyExc_RuntimeError, "System not initialized");
        return NULL;
    }

    if (alea_build_universe_index(self->sys) < 0) {
        PyErr_SetString(PyExc_RuntimeError, alea_error());
        return NULL;
    }

    Py_RETURN_NONE;
}

static PyObject* PyAleaSystem_flatten_universe(PyAleaSystemObject* self, PyObject* args) {
    int universe_id = 0;

    if (!PyArg_ParseTuple(args, "|i", &universe_id)) {
        return NULL;
    }

    if (!self->sys) {
        PyErr_SetString(PyExc_RuntimeError, "System not initialized");
        return NULL;
    }

    int result;
    sighandler_func old_sigint = install_sigint();
    Py_BEGIN_ALLOW_THREADS
    result = alea_flatten(self->sys, universe_id);
    Py_END_ALLOW_THREADS
    if (restore_sigint(old_sigint)) return NULL;

    if (result < 0) {
        PyErr_SetString(PyExc_RuntimeError, alea_error());
        return NULL;
    }

    return PyLong_FromLong(result);
}

static PyObject* PyAleaSystem_build_spatial_index(PyAleaSystemObject* self, PyObject* Py_UNUSED(ignored)) {
    if (!self->sys) {
        PyErr_SetString(PyExc_RuntimeError, "System not initialized");
        return NULL;
    }

    int result;
    sighandler_func old_sigint = install_sigint();
    Py_BEGIN_ALLOW_THREADS
    result = alea_prepare_query_acceleration(self->sys);
    Py_END_ALLOW_THREADS
    if (restore_sigint(old_sigint)) return NULL;

    if (result < 0) {
        PyErr_SetString(PyExc_RuntimeError, alea_error());
        return NULL;
    }

    Py_RETURN_NONE;
}

static PyObject* PyAleaSystem_prepare_query_acceleration(PyAleaSystemObject* self, PyObject* Py_UNUSED(ignored)) {
    if (ensure_query_acceleration(self) < 0) return NULL;
    Py_RETURN_NONE;
}

static PyObject* PyAleaSystem_query_acceleration_stats(PyAleaSystemObject* self, PyObject* Py_UNUSED(ignored)) {
    if (ensure_query_acceleration(self) < 0) return NULL;

    alea_query_acceleration_stats_t s;
    if (alea_query_acceleration_stats(self->sys, &s) < 0) {
        PyErr_SetString(PyExc_RuntimeError, alea_error());
        return NULL;
    }

    return Py_BuildValue(
        "{s:O,"
        "s:n,s:n,s:n,s:n,s:n,s:n,s:n,s:n,s:n,s:n,s:n,s:n,"
        "s:i,s:i,s:i,s:n}",
        "built", s.built ? Py_True : Py_False,
        "hier_universe_count", (Py_ssize_t)s.hier_universe_count,
        "hier_blas_count", (Py_ssize_t)s.hier_blas_count,
        "hier_linear_universe_count", (Py_ssize_t)s.hier_linear_universe_count,
        "hier_blas_cell_count", (Py_ssize_t)s.hier_blas_cell_count,
        "hier_blas_node_count", (Py_ssize_t)s.hier_blas_node_count,
        "hier_fill_cell_count", (Py_ssize_t)s.hier_fill_cell_count,
        "hier_lattice_cell_count", (Py_ssize_t)s.hier_lattice_cell_count,
        "hier_transform_count", (Py_ssize_t)s.hier_transform_count,
        "hier_placement_count", (Py_ssize_t)s.hier_placement_count,
        "hier_root_placement_count", (Py_ssize_t)s.hier_root_placement_count,
        "hier_fill_placement_count", (Py_ssize_t)s.hier_fill_placement_count,
        "hier_lattice_placement_count", (Py_ssize_t)s.hier_lattice_placement_count,
        "hier_max_placement_depth", s.hier_max_placement_depth,
        "hier_max_universe_cells", s.hier_max_universe_cells,
        "hier_largest_universe_id", s.hier_largest_universe_id,
        "memory_bytes", (Py_ssize_t)s.memory_bytes);
}

static PyObject* PyAleaSystem_surface_reference_stats(PyAleaSystemObject* self, PyObject* Py_UNUSED(ignored)) {
    if (ensure_query_acceleration(self) < 0) return NULL;

    alea_surface_reference_stats_t s;
    if (alea_surface_reference_stats(self->sys, &s) < 0) {
        PyErr_SetString(PyExc_RuntimeError, alea_error());
        return NULL;
    }

    return Py_BuildValue(
        "{s:O,s:n,s:n,s:n,s:n}",
        "built", s.built ? Py_True : Py_False,
        "surface_count", (Py_ssize_t)s.surface_count,
        "reference_count", (Py_ssize_t)s.reference_count,
        "max_references_per_surface", (Py_ssize_t)s.max_references_per_surface,
        "memory_bytes", (Py_ssize_t)s.memory_bytes);
}

static PyObject* PyAleaSystem_surface_cell_references(
        PyAleaSystemObject* self, PyObject* args) {
    int surface_id;
    if (!PyArg_ParseTuple(args, "i", &surface_id)) return NULL;
    if (ensure_query_acceleration(self) < 0) return NULL;

    size_t count = 0;
    if (alea_surface_cell_references(
            self->sys, surface_id, NULL, 0, &count) < 0) {
        PyErr_SetString(PyExc_RuntimeError, alea_error());
        return NULL;
    }
    alea_surface_cell_reference_t* refs = NULL;
    if (count > 0) {
        refs = PyMem_Malloc(count * sizeof(*refs));
        if (!refs) return PyErr_NoMemory();
        size_t copied_count = 0;
        if (alea_surface_cell_references(
                self->sys, surface_id, refs, count, &copied_count) < 0) {
            PyMem_Free(refs);
            PyErr_SetString(PyExc_RuntimeError, alea_error());
            return NULL;
        }
        count = copied_count;
    }

    PyObject* result = PyList_New((Py_ssize_t)count);
    if (!result) {
        PyMem_Free(refs);
        return NULL;
    }
    for (size_t i = 0; i < count; i++) {
        PyObject* item = Py_BuildValue(
            "{s:I,s:i,s:i,s:i}",
            "cell_index", refs[i].cell_index,
            "cell_id", refs[i].cell_id,
            "universe_id", refs[i].universe_id,
            "sense", (int)refs[i].sense);
        if (!item) {
            Py_DECREF(result);
            PyMem_Free(refs);
            return NULL;
        }
        PyList_SET_ITEM(result, (Py_ssize_t)i, item);
    }
    PyMem_Free(refs);
    return result;
}

static PyObject* PyAleaSystem_get_universe(PyAleaSystemObject* self, PyObject* args) {
    int universe_id;

    if (!PyArg_ParseTuple(args, "i", &universe_id)) {
        return NULL;
    }

    if (!self->sys) {
        PyErr_SetString(PyExc_RuntimeError, "System not initialized");
        return NULL;
    }

    int idx = alea_universe_find(self->sys, universe_id);
    if (idx < 0) {
        PyErr_Format(PyExc_KeyError, "Universe %d not found", universe_id);
        return NULL;
    }

    int uid;
    size_t cell_count;
    alea_bbox_t bbox;
    if (alea_universe_get(self->sys, (size_t)idx, &uid, &cell_count, &bbox) < 0) {
        PyErr_Format(PyExc_RuntimeError, "Failed to get universe %d info", universe_id);
        return NULL;
    }

    return Py_BuildValue("{s:i, s:n, s:(dddddd)}",
        "universe_id", uid,
        "cell_count", (Py_ssize_t)cell_count,
        "bbox", bbox.min_x, bbox.max_x,
                bbox.min_y, bbox.max_y,
                bbox.min_z, bbox.max_z);
}
