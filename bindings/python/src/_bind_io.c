// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/* This file is part of pyalea_binding.c — do NOT compile separately.
 * It is #included from pyalea_binding.c to keep a single translation unit.
 *
 * Contents: Export (MCNP, OpenMC), merge, extract/filter operations,
 *           material operations, renumbering, split/expand.
 */

/* ============================================================================
 * PyAleaSystem Methods - Export
 * ============================================================================ */

static PyObject* export_stream_to_string(FILE* stream) {
    if (fflush(stream) != 0 || fseek(stream, 0, SEEK_END) != 0) {
        PyErr_SetFromErrno(PyExc_OSError);
        return NULL;
    }
    long length = ftell(stream);
    if (length < 0 || fseek(stream, 0, SEEK_SET) != 0) {
        PyErr_SetFromErrno(PyExc_OSError);
        return NULL;
    }
    if ((unsigned long)length > (unsigned long)PY_SSIZE_T_MAX) {
        PyErr_SetString(PyExc_OverflowError, "exported text is too large for Python");
        return NULL;
    }
    char* text = PyMem_Malloc((size_t)length + 1);
    if (!text) return PyErr_NoMemory();
    size_t read = fread(text, 1, (size_t)length, stream);
    if (read != (size_t)length) {
        PyMem_Free(text);
        PyErr_SetFromErrno(PyExc_OSError);
        return NULL;
    }
    text[length] = '\0';
    PyObject* result = PyUnicode_DecodeUTF8(text, (Py_ssize_t)length, "strict");
    PyMem_Free(text);
    return result;
}

static FILE* open_export_stream(void) {
    FILE* stream = tmpfile();
    if (!stream) PyErr_SetFromErrno(PyExc_OSError);
    return stream;
}

static PyObject* PyAleaSystem_export_mcnp(PyAleaSystemObject* self, PyObject* args, PyObject* kwds) {
    const char* filename;
    int deduplicate = 1;
    int universe_depth = -1;  /* -1 = all universes */
    int fill_depth = 0;       /* 0 = no expansion */
    static char* kwlist[] = {"filename", "deduplicate", "universe_depth", "fill_depth", NULL};

    if (!PyArg_ParseTupleAndKeywords(args, kwds, "s|pii", kwlist,
            &filename, &deduplicate, &universe_depth, &fill_depth)) {
        return NULL;
    }

    if (!self->sys) {
        PyErr_SetString(PyExc_RuntimeError, "System not initialized");
        return NULL;
    }

    /* Save original config and apply export settings */
    alea_config_t orig = alea_get_config(self->sys);
    alea_config_t cfg = orig;
    cfg.dedup = deduplicate;
    cfg.universe_depth = universe_depth;
    cfg.fill_depth = fill_depth;
    alea_set_config(self->sys, &cfg);

    int result;
    sighandler_func old_sigint = install_sigint();
    Py_BEGIN_ALLOW_THREADS
    if (self->mcnp_model) result = mcnp_export(self->mcnp_model, filename);
    else if (self->alea_model) {
        mcnp_model_t* model = mcnp_model_from_alea_model(self->alea_model);
        result = model ? mcnp_export(model, filename) : -1;
        mcnp_model_destroy(model);
    } else result = mcnp_export_system(self->sys, filename);
    Py_END_ALLOW_THREADS

    /* Restore original config */
    alea_set_config(self->sys, &orig);

    if (restore_sigint(old_sigint)) return NULL;

    if (result < 0) {
        PyErr_Format(PyExc_IOError, "Failed to export to %s: %s", filename, alea_error());
        return NULL;
    }

    Py_RETURN_NONE;
}

static PyObject* PyAleaSystem_export_mcnp_string(
        PyAleaSystemObject* self, PyObject* args, PyObject* kwds) {
    int deduplicate = 1;
    int universe_depth = -1;
    int fill_depth = 0;
    static char* kwlist[] = {
        "deduplicate", "universe_depth", "fill_depth", NULL
    };
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "|pii", kwlist,
                                     &deduplicate, &universe_depth, &fill_depth))
        return NULL;
    if (!self->sys) {
        PyErr_SetString(PyExc_RuntimeError, "System not initialized");
        return NULL;
    }
    FILE* stream = open_export_stream();
    if (!stream) return NULL;

    alea_config_t original = alea_get_config(self->sys);
    alea_config_t config = original;
    config.dedup = deduplicate;
    config.universe_depth = universe_depth;
    config.fill_depth = fill_depth;
    alea_set_config(self->sys, &config);

    int rc;
    sighandler_func old_sigint = install_sigint();
    Py_BEGIN_ALLOW_THREADS
    if (self->mcnp_model) rc = mcnp_export_stream(self->mcnp_model, stream);
    else if (self->alea_model) {
        mcnp_model_t* model = mcnp_model_from_alea_model(self->alea_model);
        rc = model ? mcnp_export_stream(model, stream) : -1;
        mcnp_model_destroy(model);
    } else rc = mcnp_export_system_stream(self->sys, stream);
    Py_END_ALLOW_THREADS
    alea_set_config(self->sys, &original);

    if (restore_sigint(old_sigint)) {
        fclose(stream);
        return NULL;
    }
    if (rc != 0) {
        fclose(stream);
        PyErr_Format(PyExc_RuntimeError, "MCNP export failed: %s", alea_error());
        return NULL;
    }
    PyObject* result = export_stream_to_string(stream);
    fclose(stream);
    return result;
}

static PyObject* PyAleaSystem_export_openmc(PyAleaSystemObject* self, PyObject* args, PyObject* kwds) {
    const char* filename;
    static char* kwlist[] = {"filename", NULL};

    if (!PyArg_ParseTupleAndKeywords(args, kwds, "s", kwlist, &filename)) {
        return NULL;
    }

    if (!self->sys) {
        PyErr_SetString(PyExc_RuntimeError, "System not initialized");
        return NULL;
    }

    int result;
    sighandler_func old_sigint = install_sigint();
    Py_BEGIN_ALLOW_THREADS
    result = openmc_export_system(self->sys, filename);
    Py_END_ALLOW_THREADS
    if (restore_sigint(old_sigint)) return NULL;

    if (result < 0) {
        PyErr_Format(PyExc_IOError, "Failed to export to %s: %s", filename, alea_error());
        return NULL;
    }

    Py_RETURN_NONE;
}

static PyObject* PyAleaSystem_export_openmc_string(
        PyAleaSystemObject* self, PyObject* Py_UNUSED(ignored)) {
    if (!self->sys) {
        PyErr_SetString(PyExc_RuntimeError, "System not initialized");
        return NULL;
    }
    FILE* stream = open_export_stream();
    if (!stream) return NULL;
    int rc;
    sighandler_func old_sigint = install_sigint();
    Py_BEGIN_ALLOW_THREADS
    rc = openmc_export_system_stream(self->sys, stream);
    Py_END_ALLOW_THREADS
    if (restore_sigint(old_sigint)) {
        fclose(stream);
        return NULL;
    }
    if (rc != 0) {
        fclose(stream);
        PyErr_Format(PyExc_RuntimeError, "OpenMC export failed: %s", alea_error());
        return NULL;
    }
    PyObject* result = export_stream_to_string(stream);
    fclose(stream);
    return result;
}

static PyObject* PyAleaSystem_export_serpent(PyAleaSystemObject* self, PyObject* args, PyObject* kwds) {
    const char* filename;
    static char* kwlist[] = {"filename", NULL};

    if (!PyArg_ParseTupleAndKeywords(args, kwds, "s", kwlist, &filename)) {
        return NULL;
    }

    if (!self->sys) {
        PyErr_SetString(PyExc_RuntimeError, "System not initialized");
        return NULL;
    }

    int result;
    sighandler_func old_sigint = install_sigint();
    Py_BEGIN_ALLOW_THREADS
    result = serpent_export_system(self->sys, filename);
    Py_END_ALLOW_THREADS
    if (restore_sigint(old_sigint)) return NULL;

    if (result < 0) {
        PyErr_Format(PyExc_IOError, "Failed to export to %s: %s", filename, alea_error());
        return NULL;
    }

    Py_RETURN_NONE;
}

static PyObject* PyAleaSystem_export_serpent_string(
        PyAleaSystemObject* self, PyObject* Py_UNUSED(ignored)) {
    if (!self->sys) {
        PyErr_SetString(PyExc_RuntimeError, "System not initialized");
        return NULL;
    }
    FILE* stream = open_export_stream();
    if (!stream) return NULL;
    int rc;
    sighandler_func old_sigint = install_sigint();
    Py_BEGIN_ALLOW_THREADS
    rc = serpent_export_system_stream(self->sys, stream);
    Py_END_ALLOW_THREADS
    if (restore_sigint(old_sigint)) {
        fclose(stream);
        return NULL;
    }
    if (rc != 0) {
        fclose(stream);
        PyErr_Format(PyExc_RuntimeError, "Serpent export failed: %s", alea_error());
        return NULL;
    }
    PyObject* result = export_stream_to_string(stream);
    fclose(stream);
    return result;
}

static PyObject* PyAleaSystem_export_alea(PyAleaSystemObject* self,
                                          PyObject* args, PyObject* kwds) {
    const char* filename;
    static char* kwlist[] = {"filename", NULL};
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "s", kwlist, &filename)) return NULL;
    if (!self->sys) { PyErr_SetString(PyExc_RuntimeError, "System not initialized"); return NULL; }
    int rc;
    if (self->alea_model) rc = alea_xml_export(self->alea_model, filename);
    else if (self->mcnp_model) {
        alea_model_t* model = mcnp_model_to_alea_model(self->mcnp_model);
        rc = model ? alea_xml_export(model, filename) : -1;
        alea_model_destroy(model);
    } else rc = alea_xml_export_system(self->sys, filename);
    if (rc) { PyErr_Format(PyExc_IOError, "ALEA XML export failed: %s", alea_error()); return NULL; }
    Py_RETURN_NONE;
}

static PyObject* PyAleaSystem_export_alea_string(
        PyAleaSystemObject* self, PyObject* Py_UNUSED(ignored)) {
    if (!self->sys) { PyErr_SetString(PyExc_RuntimeError, "System not initialized"); return NULL; }
    FILE* stream = open_export_stream();
    if (!stream) return NULL;
    alea_model_t* temporary = NULL;
    int rc;
    if (self->alea_model) rc = alea_xml_export_stream(self->alea_model, stream);
    else if (self->mcnp_model) {
        temporary = mcnp_model_to_alea_model(self->mcnp_model);
        rc = temporary ? alea_xml_export_stream(temporary, stream) : -1;
    } else rc = alea_xml_export_system_stream(self->sys, stream);
    alea_model_destroy(temporary);
    if (rc) { fclose(stream); PyErr_Format(PyExc_RuntimeError, "ALEA XML export failed: %s", alea_error()); return NULL; }
    PyObject* result = export_stream_to_string(stream);
    fclose(stream);
    return result;
}

/* ============================================================================
 * PyAleaSystem Methods - Merge
 * ============================================================================ */

static PyObject* PyAleaSystem_merge(PyAleaSystemObject* self, PyObject* args, PyObject* kwds) {
    PyAleaSystemObject* other;
    int id_offset = 0;
    static char* kwlist[] = {"other", "id_offset", NULL};

    if (!PyArg_ParseTupleAndKeywords(args, kwds, "O!|i", kwlist, &PyAleaSystemType, &other, &id_offset)) {
        return NULL;
    }

    if (!self->sys || !other->sys) {
        PyErr_SetString(PyExc_RuntimeError, "System not initialized");
        return NULL;
    }

    if (alea_merge(self->sys, other->sys, id_offset) < 0) {
        PyErr_SetString(PyExc_RuntimeError, alea_error());
        return NULL;
    }

    Py_RETURN_NONE;
}

/* ============================================================================
 * PyAleaSystem Methods - Extract / Filter
 * ============================================================================ */

static PyObject* PyAleaSystem_get_cells_by_material(PyAleaSystemObject* self, PyObject* args) {
    int material_id;
    if (!PyArg_ParseTuple(args, "i", &material_id)) return NULL;

    if (!self->sys) {
        PyErr_SetString(PyExc_RuntimeError, "System not initialized");
        return NULL;
    }

    /* First call to get count */
    size_t count = alea_get_cells_by_material(self->sys, material_id, NULL, 0);

    if (count == 0) {
        return PyList_New(0);
    }

    int* indices = malloc(count * sizeof(int));
    if (!indices) return PyErr_NoMemory();

    alea_get_cells_by_material(self->sys, material_id, indices, count);

    PyObject* result = PyList_New(count);
    for (size_t i = 0; i < count; i++) {
        PyList_SET_ITEM(result, i, PyLong_FromLong(indices[i]));
    }

    free(indices);
    return result;
}

static PyObject* PyAleaSystem_get_cells_by_universe(PyAleaSystemObject* self, PyObject* args) {
    int universe_id;
    if (!PyArg_ParseTuple(args, "i", &universe_id)) return NULL;

    if (!self->sys) {
        PyErr_SetString(PyExc_RuntimeError, "System not initialized");
        return NULL;
    }

    size_t count = alea_get_cells_by_universe(self->sys, universe_id, NULL, 0);

    if (count == 0) {
        return PyList_New(0);
    }

    int* indices = malloc(count * sizeof(int));
    if (!indices) return PyErr_NoMemory();

    alea_get_cells_by_universe(self->sys, universe_id, indices, count);

    PyObject* result = PyList_New(count);
    for (size_t i = 0; i < count; i++) {
        PyList_SET_ITEM(result, i, PyLong_FromLong(indices[i]));
    }

    free(indices);
    return result;
}

static PyObject* PyAleaSystem_get_cells_filling_universe(PyAleaSystemObject* self, PyObject* args) {
    int universe_id;
    if (!PyArg_ParseTuple(args, "i", &universe_id)) return NULL;

    if (!self->sys) {
        PyErr_SetString(PyExc_RuntimeError, "System not initialized");
        return NULL;
    }

    size_t count = alea_get_cells_filling_universe(self->sys, universe_id, NULL, 0);

    if (count == 0) {
        return PyList_New(0);
    }

    int* indices = malloc(count * sizeof(int));
    if (!indices) return PyErr_NoMemory();

    alea_get_cells_filling_universe(self->sys, universe_id, indices, count);

    PyObject* result = PyList_New(count);
    for (size_t i = 0; i < count; i++) {
        PyList_SET_ITEM(result, i, PyLong_FromLong(indices[i]));
    }

    free(indices);
    return result;
}

static PyObject* PyAleaSystem_get_cells_in_bbox(PyAleaSystemObject* self, PyObject* args) {
    double x_min, x_max, y_min, y_max, z_min, z_max;
    if (!PyArg_ParseTuple(args, "dddddd", &x_min, &x_max, &y_min, &y_max, &z_min, &z_max)) {
        return NULL;
    }

    if (!self->sys) {
        PyErr_SetString(PyExc_RuntimeError, "System not initialized");
        return NULL;
    }

    alea_bbox_t bbox = { x_min, x_max, y_min, y_max, z_min, z_max };
    size_t count = alea_get_cells_in_bbox(self->sys, &bbox, NULL, 0);

    if (count == 0) {
        return PyList_New(0);
    }

    int* indices = malloc(count * sizeof(int));
    if (!indices) return PyErr_NoMemory();

    alea_get_cells_in_bbox(self->sys, &bbox, indices, count);

    PyObject* result = PyList_New(count);
    for (size_t i = 0; i < count; i++) {
        PyList_SET_ITEM(result, i, PyLong_FromLong(indices[i]));
    }

    free(indices);
    return result;
}

static PyObject* PyAleaSystem_extract_universe(PyAleaSystemObject* self, PyObject* args) {
    int universe_id;
    if (!PyArg_ParseTuple(args, "i", &universe_id)) return NULL;

    if (!self->sys) {
        PyErr_SetString(PyExc_RuntimeError, "System not initialized");
        return NULL;
    }

    alea_system_t* new_sys = alea_extract_universe(self->sys, universe_id);
    if (!new_sys) {
        PyErr_Format(PyExc_RuntimeError, "Failed to extract universe: %s", alea_error());
        return NULL;
    }

    /* Create new Python object wrapping the extracted system */
    PyAleaSystemObject* new_obj = (PyAleaSystemObject*)PyAleaSystem_new(&PyAleaSystemType, NULL, NULL);
    if (!new_obj) {
        alea_destroy(new_sys);
        return NULL;
    }

    new_obj->sys = new_sys;
    new_obj->owns_sys = 1;
    if (copy_mcnp_importance_sidecar(self, new_obj) < 0) {
        Py_DECREF(new_obj);
        return NULL;
    }

    return (PyObject*)new_obj;
}

static PyObject* PyAleaSystem_extract_region(PyAleaSystemObject* self, PyObject* args) {
    double x_min, x_max, y_min, y_max, z_min, z_max;
    if (!PyArg_ParseTuple(args, "dddddd", &x_min, &x_max, &y_min, &y_max, &z_min, &z_max)) {
        return NULL;
    }

    if (!self->sys) {
        PyErr_SetString(PyExc_RuntimeError, "System not initialized");
        return NULL;
    }

    alea_bbox_t bbox = { x_min, x_max, y_min, y_max, z_min, z_max };
    alea_system_t* new_sys = alea_extract_region(self->sys, &bbox);
    if (!new_sys) {
        PyErr_Format(PyExc_RuntimeError, "Failed to extract region: %s", alea_error());
        return NULL;
    }

    PyAleaSystemObject* new_obj = (PyAleaSystemObject*)PyAleaSystem_new(&PyAleaSystemType, NULL, NULL);
    if (!new_obj) {
        alea_destroy(new_sys);
        return NULL;
    }

    new_obj->sys = new_sys;
    new_obj->owns_sys = 1;
    if (copy_mcnp_importance_sidecar(self, new_obj) < 0) {
        Py_DECREF(new_obj);
        return NULL;
    }

    return (PyObject*)new_obj;
}

/* ============================================================================
 * PyAleaSystem Methods - Material Operations
 * ============================================================================ */

static PyObject* PyAleaSystem_create_mixture(PyAleaSystemObject* self, PyObject* args, PyObject* kwds) {
    PyObject* mat_list;
    PyObject* frac_list;
    int new_mat_id = 0;
    static char* kwlist[] = {"material_ids", "fractions", "new_id", NULL};

    if (!PyArg_ParseTupleAndKeywords(args, kwds, "OO|i", kwlist, &mat_list, &frac_list, &new_mat_id)) {
        return NULL;
    }

    if (!self->sys) {
        PyErr_SetString(PyExc_RuntimeError, "System not initialized");
        return NULL;
    }

    if (!PyList_Check(mat_list) || !PyList_Check(frac_list)) {
        PyErr_SetString(PyExc_TypeError, "material_ids and fractions must be lists");
        return NULL;
    }

    Py_ssize_t count = PyList_Size(mat_list);
    if (count != PyList_Size(frac_list)) {
        PyErr_SetString(PyExc_ValueError, "material_ids and fractions must have same length");
        return NULL;
    }

    if (count == 0) {
        PyErr_SetString(PyExc_ValueError, "At least one material required");
        return NULL;
    }

    int* mat_ids = malloc(count * sizeof(int));
    double* fractions = malloc(count * sizeof(double));
    if (!mat_ids || !fractions) {
        free(mat_ids);
        free(fractions);
        return PyErr_NoMemory();
    }

    for (Py_ssize_t i = 0; i < count; i++) {
        mat_ids[i] = (int)PyLong_AsLong(PyList_GetItem(mat_list, i));
        fractions[i] = PyFloat_AsDouble(PyList_GetItem(frac_list, i));
    }

    if (PyErr_Occurred()) {
        free(mat_ids);
        free(fractions);
        return NULL;
    }

    int result = alea_create_mixture(self->sys, mat_ids, fractions, (size_t)count, new_mat_id);

    free(mat_ids);
    free(fractions);

    if (result < 0) {
        PyErr_Format(PyExc_RuntimeError, "Failed to create mixture: %s", alea_error());
        return NULL;
    }

    return PyLong_FromLong(result);
}

/* ============================================================================
 * PyAleaSystem Methods - Renumbering
 * ============================================================================ */

static PyObject* PyAleaSystem_renumber_cells(PyAleaSystemObject* self, PyObject* args) {
    int start_id;
    if (!PyArg_ParseTuple(args, "i", &start_id)) return NULL;
    if (!self->sys) { PyErr_SetString(PyExc_RuntimeError, "System not initialized"); return NULL; }

    int result = alea_renumber_cells(self->sys, start_id);
    if (result < 0) { PyErr_SetString(PyExc_RuntimeError, alea_error()); return NULL; }
    return PyLong_FromLong(result);
}

static PyObject* PyAleaSystem_renumber_surfaces(PyAleaSystemObject* self, PyObject* args) {
    int start_id;
    if (!PyArg_ParseTuple(args, "i", &start_id)) return NULL;
    if (!self->sys) { PyErr_SetString(PyExc_RuntimeError, "System not initialized"); return NULL; }

    int result = alea_renumber_surfaces(self->sys, start_id);
    if (result < 0) { PyErr_SetString(PyExc_RuntimeError, alea_error()); return NULL; }
    return PyLong_FromLong(result);
}

static PyObject* PyAleaSystem_offset_cell_ids(PyAleaSystemObject* self, PyObject* args) {
    int offset;
    if (!PyArg_ParseTuple(args, "i", &offset)) return NULL;
    if (!self->sys) { PyErr_SetString(PyExc_RuntimeError, "System not initialized"); return NULL; }

    if (alea_offset_cell_ids(self->sys, offset) < 0) {
        PyErr_SetString(PyExc_RuntimeError, alea_error()); return NULL;
    }
    Py_RETURN_NONE;
}

static PyObject* PyAleaSystem_offset_surface_ids(PyAleaSystemObject* self, PyObject* args) {
    int offset;
    if (!PyArg_ParseTuple(args, "i", &offset)) return NULL;
    if (!self->sys) { PyErr_SetString(PyExc_RuntimeError, "System not initialized"); return NULL; }

    if (alea_offset_surface_ids(self->sys, offset) < 0) {
        PyErr_SetString(PyExc_RuntimeError, alea_error()); return NULL;
    }
    Py_RETURN_NONE;
}

static PyObject* PyAleaSystem_offset_material_ids(PyAleaSystemObject* self, PyObject* args) {
    int offset;
    if (!PyArg_ParseTuple(args, "i", &offset)) return NULL;
    if (!self->sys) { PyErr_SetString(PyExc_RuntimeError, "System not initialized"); return NULL; }

    if (alea_offset_material_ids(self->sys, offset) < 0) {
        PyErr_SetString(PyExc_RuntimeError, alea_error()); return NULL;
    }
    Py_RETURN_NONE;
}

/* ============================================================================
 * PyAleaSystem Methods - Split / Expand
 * ============================================================================ */

static PyObject* PyAleaSystem_split_union_cells(PyAleaSystemObject* self, PyObject* Py_UNUSED(ignored)) {
    if (!self->sys) { PyErr_SetString(PyExc_RuntimeError, "System not initialized"); return NULL; }

    int result = alea_split_union_cells(self->sys);
    if (result < 0) {
        PyErr_SetString(PyExc_RuntimeError, alea_error());
        return NULL;
    }
    return PyLong_FromLong(result);
}

static PyObject* PyAleaSystem_expand_macrobodies(PyAleaSystemObject* self, PyObject* Py_UNUSED(ignored)) {
    if (!self->sys) { PyErr_SetString(PyExc_RuntimeError, "System not initialized"); return NULL; }

    int result = alea_expand_macrobodies_in_system(self->sys);
    if (result < 0) {
        PyErr_SetString(PyExc_RuntimeError, alea_error());
        return NULL;
    }
    return PyLong_FromLong(result);
}
