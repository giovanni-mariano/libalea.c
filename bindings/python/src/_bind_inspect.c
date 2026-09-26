// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/* This file is part of pyalea_binding.c — do NOT compile separately.
 * It is #included from pyalea_binding.c to keep a single translation unit.
 *
 * Contents: Cell lookup, find all cells at point, set fill, cell/surface/node info,
 *           node inspection, bbox tightening, CSG simplification, numerical bbox,
 *           primitive data, CSG node tree.
 */

/* ============================================================================
 * PyAleaSystem Methods - Cell Lookup by ID
 * ============================================================================ */

static PyObject* PyAleaSystem_cell_find(PyAleaSystemObject* self, PyObject* args) {
    int cell_id;

    if (!PyArg_ParseTuple(args, "i", &cell_id)) {
        return NULL;
    }

    if (!self->sys) {
        PyErr_SetString(PyExc_RuntimeError, "System not initialized");
        return NULL;
    }

    int idx = alea_cell_find(self->sys, cell_id);
    if (idx < 0) {
        Py_RETURN_NONE;
    }
    return PyLong_FromLong(idx);
}

/* ============================================================================
 * PyAleaSystem Methods - Find All Cells at Point (Hierarchy)
 * ============================================================================ */

static PyObject* PyAleaSystem_find_all_cells(PyAleaSystemObject* self, PyObject* args) {
    double x, y, z;

    if (!PyArg_ParseTuple(args, "ddd", &x, &y, &z)) {
        return NULL;
    }

    if (!self->sys) {
        PyErr_SetString(PyExc_RuntimeError, "System not initialized");
        return NULL;
    }
    if (ensure_query_acceleration(self) < 0) return NULL;

    /* Allocate hits on stack - 64 levels should be more than enough */
    alea_cell_hit_t hits[64];
    int nhits = alea_find_all_cells(self->sys, x, y, z, hits, 64);
    if (nhits < 0) {
        PyErr_SetString(PyExc_RuntimeError, "find_all_cells failed");
        return NULL;
    }

    PyObject* list = PyList_New(nhits);
    if (!list) return NULL;

    for (int i = 0; i < nhits; i++) {
        PyObject* hit = Py_BuildValue("{s:i, s:i, s:i, s:i, s:i, s:i, s:d, s:d, s:d}",
            "cell_id", hits[i].cell_id,
            "cell_index", hits[i].cell_index,
            "material_id", hits[i].material_id,
            "universe_id", hits[i].universe_id,
            "fill_universe", hits[i].fill_universe,
            "depth", hits[i].depth,
            "local_x", hits[i].local_x,
            "local_y", hits[i].local_y,
            "local_z", hits[i].local_z);
        if (!hit) {
            Py_DECREF(list);
            return NULL;
        }
        PyList_SET_ITEM(list, i, hit);
    }

    return list;
}

/* Complete diagnostic ownership query. Unlike find_all_cells(), this always
 * uses the recursive coverage resolver and retains occurrence-parent links. */
static PyObject* coverage_hit_dict(const alea_cell_hit_t* hit,
                                   uint64_t occurrence_key,
                                   uint64_t parent_key) {
    return Py_BuildValue(
        "{s:i,s:i,s:i,s:i,s:i,s:i,s:d,s:d,s:d,s:i,s:K,s:K}",
        "cell_id", hit->cell_id, "cell_index", hit->cell_index,
        "material_id", hit->material_id, "universe_id", hit->universe_id,
        "fill_universe", hit->fill_universe, "depth", hit->depth,
        "local_x", hit->local_x, "local_y", hit->local_y,
        "local_z", hit->local_z,
        "resolution_flags", (int)hit->resolution_flags,
        "occurrence_key", (unsigned long long)occurrence_key,
        "parent_occurrence_key", (unsigned long long)parent_key);
}

static const char* point_coverage_kind_name(alea_point_coverage_kind_t kind) {
    switch (kind) {
        case ALEA_POINT_COVERAGE_UNIQUE: return "unique";
        case ALEA_POINT_COVERAGE_GAP: return "gap";
        case ALEA_POINT_COVERAGE_OVERLAP: return "overlap";
        case ALEA_POINT_COVERAGE_UNDEFINED_FILL: return "undefined_fill";
        case ALEA_POINT_COVERAGE_UNRESOLVED: return "unresolved";
    }
    return "unresolved";
}

static PyObject* PyAleaSystem_find_all_cells_coverage(
    PyAleaSystemObject* self, PyObject* args, PyObject* kwds) {
    double x, y, z;
    Py_ssize_t max_hits = 256;
    int universe_id = 0;
    int universe_depth = -1;
    static char* kwlist[] = {
        "x", "y", "z", "max_hits", "universe_id", "universe_depth", NULL
    };
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "ddd|nii", kwlist,
                                     &x, &y, &z, &max_hits,
                                     &universe_id, &universe_depth))
        return NULL;
    if (!self->sys) {
        PyErr_SetString(PyExc_RuntimeError, "System not initialized");
        return NULL;
    }
    if (max_hits <= 0 || max_hits > 16384) {
        PyErr_SetString(PyExc_ValueError, "max_hits must be between 1 and 16384");
        return NULL;
    }
    if (ensure_query_acceleration(self) < 0) return NULL;

    const size_t capacity = (size_t)max_hits;
    alea_cell_hit_t* hits = PyMem_Malloc(capacity * sizeof(*hits));
    uint64_t* occurrence_keys = PyMem_Malloc(capacity * sizeof(*occurrence_keys));
    uint64_t* parent_keys = PyMem_Malloc(capacity * sizeof(*parent_keys));
    uint8_t* owner_mask = PyMem_Calloc(capacity, sizeof(*owner_mask));
    if (!hits || !occurrence_keys || !parent_keys || !owner_mask) {
        PyMem_Free(hits); PyMem_Free(occurrence_keys); PyMem_Free(parent_keys);
        PyMem_Free(owner_mask);
        return PyErr_NoMemory();
    }
    int nhits = alea_find_all_cells_in_universe_coverage_chain(
        self->sys, universe_id, x, y, z,
        hits, occurrence_keys, parent_keys, capacity);
    if (nhits < 0) {
        PyMem_Free(hits); PyMem_Free(occurrence_keys); PyMem_Free(parent_keys);
        PyMem_Free(owner_mask);
        PyErr_SetString(PyExc_RuntimeError, alea_error());
        return NULL;
    }
    PyObject* list = PyList_New(nhits);
    PyObject* owners = PyList_New(0);
    if (!list || !owners) {
        Py_XDECREF(list); Py_XDECREF(owners);
        PyMem_Free(hits); PyMem_Free(occurrence_keys); PyMem_Free(parent_keys);
        PyMem_Free(owner_mask);
        return NULL;
    }
    for (int i = 0; i < nhits; i++) {
        PyObject* hit = coverage_hit_dict(&hits[i], occurrence_keys[i], parent_keys[i]);
        if (!hit) {
            Py_DECREF(list); Py_DECREF(owners);
            PyMem_Free(hits); PyMem_Free(occurrence_keys); PyMem_Free(parent_keys);
            PyMem_Free(owner_mask);
            return NULL;
        }
        PyList_SET_ITEM(list, i, hit);
    }

    const int truncated = nhits >= max_hits;
    alea_point_coverage_classification_t classification = {
        .kind = ALEA_POINT_COVERAGE_UNRESOLVED,
        .target_depth = universe_depth,
        .owner_count = 0
    };
    if (!truncated && alea_classify_point_coverage_chain(
            hits, occurrence_keys, parent_keys, (size_t)nhits,
            universe_depth, owner_mask, &classification) != 0) {
        Py_DECREF(list); Py_DECREF(owners);
        PyMem_Free(hits); PyMem_Free(occurrence_keys); PyMem_Free(parent_keys);
        PyMem_Free(owner_mask);
        PyErr_SetString(PyExc_RuntimeError, "point coverage classification failed");
        return NULL;
    }
    if (!truncated) {
        for (int i = 0; i < nhits; i++) {
            if (!owner_mask[i]) continue;
            PyObject* owner = coverage_hit_dict(
                &hits[i], occurrence_keys[i], parent_keys[i]);
            if (!owner || PyList_Append(owners, owner) != 0) {
                Py_XDECREF(owner); Py_DECREF(list); Py_DECREF(owners);
                PyMem_Free(hits); PyMem_Free(occurrence_keys);
                PyMem_Free(parent_keys); PyMem_Free(owner_mask);
                return NULL;
            }
            Py_DECREF(owner);
        }
    }
    const char* kind = truncated ? "truncated"
        : point_coverage_kind_name(classification.kind);
    PyMem_Free(hits); PyMem_Free(occurrence_keys); PyMem_Free(parent_keys);
    PyMem_Free(owner_mask);
    return Py_BuildValue(
        "{s:N,s:N,s:s,s:i,s:i,s:O,s:i,s:n,s:n}",
        "hits", list, "owners", owners, "kind", kind,
        "universe_id", universe_id,
        "target_depth", classification.target_depth,
        "truncated", truncated ? Py_True : Py_False,
        "hit_count", nhits, "max_hits", max_hits,
        "owner_count", (Py_ssize_t)classification.owner_count);
}

/* ============================================================================
 * PyAleaSystem Methods - Transforms
 * ============================================================================ */

static int parse_transform_values(PyObject* object, double values[13],
                                  int* value_count) {
    PyObject* sequence = PySequence_Fast(
        object, "values must be a sequence of MCNP transform values");
    if (!sequence) return -1;

    Py_ssize_t count = PySequence_Fast_GET_SIZE(sequence);
    if (count < 3 || count > 13) {
        Py_DECREF(sequence);
        PyErr_SetString(
            PyExc_ValueError,
            "values must contain between 3 and 13 MCNP transform values");
        return -1;
    }

    for (Py_ssize_t i = 0; i < count; ++i) {
        values[i] = PyFloat_AsDouble(PySequence_Fast_GET_ITEM(sequence, i));
        if (PyErr_Occurred()) {
            Py_DECREF(sequence);
            return -1;
        }
    }
    Py_DECREF(sequence);
    *value_count = (int)count;
    return 0;
}

static PyObject* PyAleaSystem_add_transform(
        PyAleaSystemObject* self, PyObject* args, PyObject* kwds) {
    int transform_id;
    PyObject* values_object;
    int degrees = 0;
    static char* kwlist[] = {
        "transform_id", "values", "degrees", NULL
    };
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "iO|p", kwlist,
                                     &transform_id, &values_object, &degrees))
        return NULL;
    if (!self->sys) {
        PyErr_SetString(PyExc_RuntimeError, "System not initialized");
        return NULL;
    }

    double values[13];
    int value_count;
    if (parse_transform_values(values_object, values, &value_count) < 0)
        return NULL;
    if (alea_add_transform(self->sys, transform_id, values,
                           value_count, degrees) != 0) {
        PyErr_SetString(PyExc_ValueError, alea_error());
        return NULL;
    }
    return PyLong_FromLong(transform_id);
}

static PyObject* PyAleaSystem_add_inline_transform(
        PyAleaSystemObject* self, PyObject* args, PyObject* kwds) {
    PyObject* values_object;
    int degrees = 0;
    int cell_id = 0;
    const char* role = "fill";
    static char* kwlist[] = {
        "values", "degrees", "cell_id", "role", NULL
    };
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "O|piz", kwlist,
                                     &values_object, &degrees, &cell_id, &role))
        return NULL;
    if (!self->sys) {
        PyErr_SetString(PyExc_RuntimeError, "System not initialized");
        return NULL;
    }

    double values[13];
    int value_count;
    if (parse_transform_values(values_object, values, &value_count) < 0)
        return NULL;
    int transform_id = alea_add_inline_transform(
        self->sys, values, value_count, degrees, cell_id, role);
    if (transform_id < 0) {
        PyErr_SetString(PyExc_ValueError, alea_error());
        return NULL;
    }
    return PyLong_FromLong(transform_id);
}

/* ============================================================================
 * PyAleaSystem Methods - Set Fill
 * ============================================================================ */

static PyObject* PyAleaSystem_set_fill(PyAleaSystemObject* self, PyObject* args) {
    int cell_index, fill_universe, transform = 0;
    if (!PyArg_ParseTuple(args, "ii|i", &cell_index, &fill_universe, &transform)) return NULL;
    if (!self->sys) { PyErr_SetString(PyExc_RuntimeError, "System not initialized"); return NULL; }

    if (alea_set_fill(self->sys, cell_index, fill_universe, transform) < 0) {
        PyErr_Format(PyExc_RuntimeError, "Failed to set fill: %s", alea_error());
        return NULL;
    }
    Py_RETURN_NONE;
}

/* ============================================================================
 * PyAleaSystem Methods - Cell Comments
 * ============================================================================ */

static PyObject* PyAleaSystem_set_comment(PyAleaSystemObject* self, PyObject* args) {
    int cell_index;
    const char* comment;
    if (!PyArg_ParseTuple(args, "iz", &cell_index, &comment)) return NULL;
    if (!self->sys) { PyErr_SetString(PyExc_RuntimeError, "System not initialized"); return NULL; }

    if (alea_cell_set_comment(self->sys, cell_index, comment) < 0) {
        PyErr_Format(PyExc_IndexError, "Cell index %d out of range", cell_index);
        return NULL;
    }
    Py_RETURN_NONE;
}

static PyObject* PyAleaSystem_set_inline_comment(PyAleaSystemObject* self, PyObject* args) {
    int cell_index;
    const char* comment;
    if (!PyArg_ParseTuple(args, "iz", &cell_index, &comment)) return NULL;
    if (!self->sys) { PyErr_SetString(PyExc_RuntimeError, "System not initialized"); return NULL; }

    if (alea_cell_set_inline_comment(self->sys, cell_index, comment) < 0) {
        PyErr_Format(PyExc_IndexError, "Cell index %d out of range", cell_index);
        return NULL;
    }
    Py_RETURN_NONE;
}

/* ============================================================================
 * PyAleaSystem Methods - Cell Property Setters
 * ============================================================================ */

static PyObject* PyAleaSystem_cell_set_material(PyAleaSystemObject* self, PyObject* args) {
    int cell_index, material_index;
    if (!PyArg_ParseTuple(args, "ii", &cell_index, &material_index)) return NULL;
    if (!self->sys) { PyErr_SetString(PyExc_RuntimeError, "System not initialized"); return NULL; }

    if (alea_cell_set_material(self->sys, cell_index, material_index) < 0) {
        PyErr_Format(PyExc_IndexError, "Cell index %d out of range", cell_index);
        return NULL;
    }
    Py_RETURN_NONE;
}

static PyObject* PyAleaSystem_cell_set_density(PyAleaSystemObject* self, PyObject* args) {
    int cell_index;
    double density;
    if (!PyArg_ParseTuple(args, "id", &cell_index, &density)) return NULL;
    if (!self->sys) { PyErr_SetString(PyExc_RuntimeError, "System not initialized"); return NULL; }

    if (alea_cell_set_density(self->sys, cell_index, density) < 0) {
        PyErr_Format(PyExc_IndexError, "Cell index %d out of range", cell_index);
        return NULL;
    }
    Py_RETURN_NONE;
}

static PyObject* PyAleaSystem_cell_set_temperature(
        PyAleaSystemObject* self, PyObject* args) {
    int cell_index;
    double temperature_K;
    if (!PyArg_ParseTuple(args, "id", &cell_index, &temperature_K)) return NULL;
    if (!self->sys) {
        PyErr_SetString(PyExc_RuntimeError, "System not initialized");
        return NULL;
    }
    if (alea_cell_set_temperature(self->sys, cell_index, temperature_K) < 0) {
        PyErr_SetString(PyExc_ValueError,
                        "temperature must be finite and positive and the cell must exist");
        return NULL;
    }
    Py_RETURN_NONE;
}

static PyObject* PyAleaSystem_cell_clear_temperature(
        PyAleaSystemObject* self, PyObject* args) {
    int cell_index;
    if (!PyArg_ParseTuple(args, "i", &cell_index)) return NULL;
    if (!self->sys) {
        PyErr_SetString(PyExc_RuntimeError, "System not initialized");
        return NULL;
    }
    if (alea_cell_clear_temperature(self->sys, cell_index) < 0) {
        PyErr_Format(PyExc_IndexError, "Cell index %d out of range", cell_index);
        return NULL;
    }
    Py_RETURN_NONE;
}

static PyObject* PyAleaSystem_cell_set_importance(
    PyAleaSystemObject* self, PyObject* args) {
    int cell_index;
    const char* particle;
    double importance;
    if (!PyArg_ParseTuple(args, "isd", &cell_index, &particle, &importance)) return NULL;
    if (self->alea_model) {
        alea_model_cell_metadata_t* meta =
            alea_model_cell_metadata_mut(self->alea_model, (size_t)cell_index);
        if (!meta || !isfinite(importance)) {
            PyErr_SetString(PyExc_ValueError, "invalid cell index or particle importance");
            return NULL;
        }
        if (strcmp(particle, "neutron") == 0) {
            meta->importance_neutron = importance;
            meta->has_importance_neutron = 1;
        } else if (strcmp(particle, "photon") == 0) {
            meta->importance_photon = importance;
            meta->has_importance_photon = 1;
        } else if (strcmp(particle, "electron") == 0) {
            meta->importance_electron = importance;
            meta->has_importance_electron = 1;
        } else {
            PyErr_SetString(PyExc_ValueError,
                            "particle must be 'neutron', 'photon', or 'electron'");
            return NULL;
        }
        Py_RETURN_NONE;
    }
    mcnp_model_t* model = ensure_mcnp_sidecar(self);
    mcnp_cell_params_t* params = model
        ? mcnp_cell_params(model, (size_t)cell_index) : NULL;
    if (!params || !isfinite(importance)) {
        PyErr_SetString(PyExc_ValueError, "invalid cell index or particle importance");
        return NULL;
    }
    if (strcmp(particle, "neutron") == 0) {
        params->imp_n = importance;
        params->has_imp_n = 1;
    } else if (strcmp(particle, "photon") == 0) {
        params->imp_p = importance;
        params->has_imp_p = 1;
    } else if (strcmp(particle, "electron") == 0) {
        params->imp_e = importance;
        params->has_imp_e = 1;
    } else {
        PyErr_SetString(PyExc_ValueError,
                        "particle must be 'neutron', 'photon', or 'electron'");
        return NULL;
    }
    Py_RETURN_NONE;
}

static PyObject* PyAleaSystem_cell_set_universe(PyAleaSystemObject* self, PyObject* args) {
    int cell_index, universe_id;
    if (!PyArg_ParseTuple(args, "ii", &cell_index, &universe_id)) return NULL;
    if (!self->sys) { PyErr_SetString(PyExc_RuntimeError, "System not initialized"); return NULL; }

    if (alea_cell_set_universe(self->sys, cell_index, universe_id) < 0) {
        PyErr_Format(PyExc_IndexError, "Cell index %d out of range", cell_index);
        return NULL;
    }
    Py_RETURN_NONE;
}

static PyObject* PyAleaSystem_cell_set_region(PyAleaSystemObject* self,
                                              PyObject* args) {
    int cell_index;
    unsigned int root;
    if (!PyArg_ParseTuple(args, "iI", &cell_index, &root)) return NULL;
    if (!self->sys) {
        PyErr_SetString(PyExc_RuntimeError, "System not initialized");
        return NULL;
    }
    if (alea_cell_set_region(self->sys, cell_index, (alea_node_id_t)root) < 0) {
        PyErr_SetString(PyExc_ValueError, "invalid cell index or CSG root");
        return NULL;
    }
    Py_RETURN_NONE;
}

static PyObject* PyAleaSystem_cell_remove(PyAleaSystemObject* self, PyObject* args) {
    int cell_index;
    if (!PyArg_ParseTuple(args, "i", &cell_index)) return NULL;
    if (!self->sys) { PyErr_SetString(PyExc_RuntimeError, "System not initialized"); return NULL; }

    if (alea_cell_remove(self->sys, cell_index) < 0) {
        PyErr_Format(PyExc_IndexError, "Cell index %d out of range", cell_index);
        return NULL;
    }
    Py_RETURN_NONE;
}

/* ============================================================================
 * PyAleaSystem Methods - Cell/Surface/Node Info
 * ============================================================================ */

static PyObject* PyAleaSystem_get_cell_id(PyAleaSystemObject* self, PyObject* args) {
    int cell_index;
    if (!PyArg_ParseTuple(args, "i", &cell_index)) return NULL;
    if (!self->sys) { PyErr_SetString(PyExc_RuntimeError, "System not initialized"); return NULL; }

    int cell_id = alea_get_cell_id(self->sys, cell_index);
    if (cell_id < 0) {
        PyErr_Format(PyExc_IndexError, "Cell index %d out of range", cell_index);
        return NULL;
    }
    return PyLong_FromLong(cell_id);
}

static PyObject* PyAleaSystem_surface_find(PyAleaSystemObject* self, PyObject* args) {
    int surface_id;
    if (!PyArg_ParseTuple(args, "i", &surface_id)) return NULL;
    if (!self->sys) { PyErr_SetString(PyExc_RuntimeError, "System not initialized"); return NULL; }

    int idx = alea_surface_find(self->sys, surface_id);
    if (idx < 0) Py_RETURN_NONE;
    return PyLong_FromLong(idx);
}

static PyObject* PyAleaSystem_surface_node(PyAleaSystemObject* self, PyObject* args) {
    int surface_id, sense;
    if (!PyArg_ParseTuple(args, "ii", &surface_id, &sense)) return NULL;
    if (!self->sys) { PyErr_SetString(PyExc_RuntimeError, "System not initialized"); return NULL; }

    int idx = alea_surface_find(self->sys, surface_id);
    if (idx < 0) Py_RETURN_NONE;
    alea_node_id_t pos_node, neg_node;
    alea_surface_get(self->sys, idx, NULL, NULL, &pos_node, &neg_node, NULL);
    return PyLong_FromUnsignedLong(sense > 0 ? pos_node : neg_node);
}

static int inspect_parse_vec3(PyObject* object, double out[3],
                              const char* name) {
    PyObject* sequence = PySequence_Fast(object, name);
    if (!sequence) return -1;
    if (PySequence_Fast_GET_SIZE(sequence) != 3) {
        Py_DECREF(sequence);
        PyErr_Format(PyExc_ValueError, "%s must contain exactly 3 values", name);
        return -1;
    }
    for (Py_ssize_t index = 0; index < 3; ++index) {
        out[index] = PyFloat_AsDouble(PySequence_Fast_GET_ITEM(sequence, index));
        if (PyErr_Occurred()) {
            Py_DECREF(sequence);
            return -1;
        }
    }
    Py_DECREF(sequence);
    return 0;
}

static PyObject* PyAleaSystem_surface_project_along(
        PyAleaSystemObject* self, PyObject* args) {
    int surface_id;
    PyObject *point_obj, *direction_obj;
    if (!PyArg_ParseTuple(args, "iOO", &surface_id, &point_obj, &direction_obj))
        return NULL;
    if (!self->sys) {
        PyErr_SetString(PyExc_RuntimeError, "System not initialized");
        return NULL;
    }
    double point[3], direction[3];
    if (inspect_parse_vec3(point_obj, point, "point") < 0 ||
        inspect_parse_vec3(direction_obj, direction, "direction") < 0)
        return NULL;
    double parameter, projected[3];
    alea_primitive_type_t primitive_type = 0;
    int rc = alea_surface_project_along(
        self->sys, surface_id, point, direction, &parameter, projected,
        &primitive_type);
    if (rc < 0) {
        PyErr_SetString(PyExc_ValueError,
                        "invalid direction or unknown surface ID");
        return NULL;
    }
    if (rc == 0) Py_RETURN_NONE;
    return Py_BuildValue("{s:d,s:(ddd),s:i}",
                         "trajectory_parameter", parameter,
                         "point", projected[0], projected[1], projected[2],
                         "primitive_type", (int)primitive_type);
}

static PyObject* PyAleaSystem_surface_id_at(PyAleaSystemObject* self, PyObject* args) {
    Py_ssize_t idx;
    if (!PyArg_ParseTuple(args, "n", &idx)) return NULL;
    if (!self->sys) { PyErr_SetString(PyExc_RuntimeError, "System not initialized"); return NULL; }

    int sid = alea_surface_id_at(self->sys, (size_t)idx);
    if (sid < 0) {
        PyErr_Format(PyExc_IndexError, "Surface index %zd out of range", idx);
        return NULL;
    }
    return PyLong_FromLong(sid);
}

static PyObject* PyAleaSystem_get_surface_ids(PyAleaSystemObject* self, PyObject* Py_UNUSED(ignored)) {
    if (!self->sys) { PyErr_SetString(PyExc_RuntimeError, "System not initialized"); return NULL; }

    size_t count = alea_surface_count(self->sys);
    if (count == 0) return PyList_New(0);

    int* buf = (int*)PyMem_Malloc(count * sizeof(int));
    if (!buf) return PyErr_NoMemory();

    size_t n = alea_get_surface_ids(self->sys, buf);
    PyObject* list = PyList_New((Py_ssize_t)n);
    if (!list) { PyMem_Free(buf); return NULL; }
    for (size_t i = 0; i < n; i++) {
        PyList_SET_ITEM(list, (Py_ssize_t)i, PyLong_FromLong(buf[i]));
    }
    PyMem_Free(buf);
    return list;
}

static PyObject* PyAleaSystem_cells_in_universe(PyAleaSystemObject* self, PyObject* args) {
    int universe_id;
    if (!PyArg_ParseTuple(args, "i", &universe_id)) return NULL;
    if (!self->sys) { PyErr_SetString(PyExc_RuntimeError, "System not initialized"); return NULL; }

    /* First call to get count */
    int count = alea_cells_in_universe(self->sys, universe_id, NULL, 0);
    if (count <= 0) return PyList_New(0);

    int* indices = malloc(count * sizeof(int));
    if (!indices) return PyErr_NoMemory();
    alea_cells_in_universe(self->sys, universe_id, indices, count);

    PyObject* result = PyList_New(count);
    for (int i = 0; i < count; i++) {
        PyList_SET_ITEM(result, i, PyLong_FromLong(indices[i]));
    }
    free(indices);
    return result;
}

static PyObject* PyAleaSystem_find_cell_at(PyAleaSystemObject* self, PyObject* args) {
    double x, y, z;
    if (!PyArg_ParseTuple(args, "ddd", &x, &y, &z)) return NULL;
    if (!self->sys) { PyErr_SetString(PyExc_RuntimeError, "System not initialized"); return NULL; }
    if (ensure_query_acceleration(self) < 0) return NULL;

    int cell_id, material;
    int result = alea_find_cell_at(self->sys, x, y, z, &cell_id, &material);
    if (result < 0) Py_RETURN_NONE;
    return Py_BuildValue("(ii)", cell_id, material);
}

/* ============================================================================
 * PyAleaSystem Methods - Node Inspection
 * ============================================================================ */

static PyObject* PyAleaSystem_node_primitive_type(PyAleaSystemObject* self, PyObject* args) {
    unsigned long node_id;
    if (!PyArg_ParseTuple(args, "k", &node_id)) return NULL;
    if (!self->sys) { PyErr_SetString(PyExc_RuntimeError, "System not initialized"); return NULL; }

    alea_primitive_type_t ptype = alea_node_primitive_type(self->sys, (alea_node_id_t)node_id);
    return PyLong_FromLong((int)ptype);
}

static PyObject* PyAleaSystem_node_primitive_id(PyAleaSystemObject* self, PyObject* args) {
    unsigned long node_id;
    if (!PyArg_ParseTuple(args, "k", &node_id)) return NULL;
    if (!self->sys) { PyErr_SetString(PyExc_RuntimeError, "System not initialized"); return NULL; }

    alea_primitive_id_t pid = alea_node_primitive_id(self->sys, (alea_node_id_t)node_id);
    if (pid == ALEA_PRIMITIVE_ID_INVALID) Py_RETURN_NONE;
    return PyLong_FromUnsignedLong(pid);
}

static PyObject* PyAleaSystem_node_sense(PyAleaSystemObject* self, PyObject* args) {
    unsigned long node_id;
    if (!PyArg_ParseTuple(args, "k", &node_id)) return NULL;
    if (!self->sys) { PyErr_SetString(PyExc_RuntimeError, "System not initialized"); return NULL; }

    int sense = alea_node_sense(self->sys, (alea_node_id_t)node_id);
    return PyLong_FromLong(sense);
}

static PyObject* PyAleaSystem_node_surface_id(PyAleaSystemObject* self, PyObject* args) {
    unsigned long node_id;
    if (!PyArg_ParseTuple(args, "k", &node_id)) return NULL;
    if (!self->sys) { PyErr_SetString(PyExc_RuntimeError, "System not initialized"); return NULL; }

    int sid = alea_node_surface_id(self->sys, (alea_node_id_t)node_id);
    return PyLong_FromLong(sid);
}

static PyObject* PyAleaSystem_node_operation(PyAleaSystemObject* self, PyObject* args) {
    unsigned long node_id;
    if (!PyArg_ParseTuple(args, "k", &node_id)) return NULL;
    if (!self->sys) { PyErr_SetString(PyExc_RuntimeError, "System not initialized"); return NULL; }

    alea_operation_t op = alea_node_operation(self->sys, (alea_node_id_t)node_id);
    return PyLong_FromLong((int)op);
}

static PyObject* PyAleaSystem_node_left(PyAleaSystemObject* self, PyObject* args) {
    unsigned long node_id;
    if (!PyArg_ParseTuple(args, "k", &node_id)) return NULL;
    if (!self->sys) { PyErr_SetString(PyExc_RuntimeError, "System not initialized"); return NULL; }

    alea_node_id_t left = alea_node_left(self->sys, (alea_node_id_t)node_id);
    if (left == ALEA_NODE_ID_INVALID) Py_RETURN_NONE;
    return PyLong_FromUnsignedLong(left);
}

static PyObject* PyAleaSystem_node_right(PyAleaSystemObject* self, PyObject* args) {
    unsigned long node_id;
    if (!PyArg_ParseTuple(args, "k", &node_id)) return NULL;
    if (!self->sys) { PyErr_SetString(PyExc_RuntimeError, "System not initialized"); return NULL; }

    alea_node_id_t right = alea_node_right(self->sys, (alea_node_id_t)node_id);
    if (right == ALEA_NODE_ID_INVALID) Py_RETURN_NONE;
    return PyLong_FromUnsignedLong(right);
}

/* ============================================================================
 * PyAleaSystem Methods - BBox Tightening
 * ============================================================================ */

static PyObject* PyAleaSystem_tighten_cell_bbox(PyAleaSystemObject* self, PyObject* args) {
    size_t cell_index;
    double tol = 1.0;
    if (!PyArg_ParseTuple(args, "n|d", &cell_index, &tol)) return NULL;
    if (!self->sys) { PyErr_SetString(PyExc_RuntimeError, "System not initialized"); return NULL; }

    alea_bbox_t out;
    if (alea_tighten_cell_bbox(self->sys, cell_index, tol, &out) < 0) {
        PyErr_SetString(PyExc_RuntimeError, alea_error());
        return NULL;
    }
    return Py_BuildValue("(dddddd)", out.min_x, out.max_x, out.min_y, out.max_y, out.min_z, out.max_z);
}

static PyObject* PyAleaSystem_tighten_all_bboxes(PyAleaSystemObject* self, PyObject* args) {
    double tol = 1.0;
    if (!PyArg_ParseTuple(args, "|d", &tol)) return NULL;
    if (!self->sys) { PyErr_SetString(PyExc_RuntimeError, "System not initialized"); return NULL; }

    int result;
    sighandler_func old_sigint = install_sigint();
    Py_BEGIN_ALLOW_THREADS
    result = alea_tighten_all_bboxes(self->sys, tol);
    Py_END_ALLOW_THREADS
    if (restore_sigint(old_sigint)) return NULL;

    if (result < 0) {
        PyErr_SetString(PyExc_RuntimeError, alea_error());
        return NULL;
    }
    return PyLong_FromLong(result);
}

/* ============================================================================
 * PyAleaSystem Methods - CSG Simplification
 * ============================================================================ */

static PyObject* PyAleaSystem_simplify_all_cells(PyAleaSystemObject* self, PyObject* Py_UNUSED(ignored)) {
    if (!self->sys) {
        PyErr_SetString(PyExc_RuntimeError, "System not initialized");
        return NULL;
    }

    alea_simplify_stats_t stats;
    memset(&stats, 0, sizeof(stats));

    sighandler_func old_sigint = install_sigint();
    Py_BEGIN_ALLOW_THREADS
    alea_simplify_and_prune_cells(self->sys, &stats);
    Py_END_ALLOW_THREADS
    if (restore_sigint(old_sigint)) return NULL;

    return Py_BuildValue(
        "{s:n, s:n, s:n, s:n, s:n, s:n, s:n, s:n, s:n, s:n, s:n, s:n, s:n, s:n}",
        "nodes_before",             (Py_ssize_t)stats.nodes_before,
        "nodes_after",              (Py_ssize_t)stats.nodes_after,
        "complements_eliminated",   (Py_ssize_t)stats.complements_eliminated,
        "double_negations",         (Py_ssize_t)stats.double_negations,
        "idempotent_reductions",    (Py_ssize_t)stats.idempotent_reductions,
        "absorption_reductions",    (Py_ssize_t)stats.absorption_reductions,
        "subtrees_deduplicated",    (Py_ssize_t)stats.subtrees_deduplicated,
        "cell_complements_expanded",(Py_ssize_t)stats.cell_complements_expanded,
        "contradictions_found",     (Py_ssize_t)stats.contradictions_found,
        "tautologies_found",        (Py_ssize_t)stats.tautologies_found,
        "empty_cells_removed",      (Py_ssize_t)stats.empty_cells_removed,
        "union_branches_absorbed",  (Py_ssize_t)stats.union_branches_absorbed,
        "union_common_factors",     (Py_ssize_t)stats.union_common_factors,
        "union_branches_subsumed",  (Py_ssize_t)stats.union_branches_subsumed
    );
}

static const char* proof_bounds_source_name(alea_proof_bounds_source_t source) {
    switch (source) {
        case ALEA_PROOF_BOUNDS_STORED: return "stored";
        case ALEA_PROOF_BOUNDS_PLANE_CONSTRAINTS: return "plane_constraints";
        case ALEA_PROOF_BOUNDS_EXPLICIT: return "explicit";
        default: return "unknown";
    }
}

static const char* proof_limit_name(alea_proof_limit_t limit) {
    switch (limit) {
        case ALEA_PROOF_LIMIT_NONE: return "none";
        case ALEA_PROOF_LIMIT_DEPTH: return "depth";
        case ALEA_PROOF_LIMIT_NODES: return "nodes";
        case ALEA_PROOF_LIMIT_MEMORY: return "memory";
        case ALEA_PROOF_LIMIT_DOMAIN: return "domain";
        case ALEA_PROOF_LIMIT_UNSUPPORTED: return "unsupported";
        default: return "unknown";
    }
}

static PyObject* proof_simplification_result_dict(
        const alea_cell_simplify_proof_result_t* result,
        bool has_explicit_bounds) {
    PyObject* witness = result->has_witness
        ? Py_BuildValue("(ddd)", result->witness[0], result->witness[1],
                        result->witness[2])
        : Py_NewRef(Py_None);
    PyObject* stats = Py_BuildValue(
        "{s:n,s:n,s:n,s:n,s:n,s:n,s:n,s:n,s:n,s:n,s:n,s:n,s:n,s:n,s:K}",
        "nodes_before", (Py_ssize_t)result->nodes_before,
        "nodes_after", (Py_ssize_t)result->nodes_after,
        "surfaces_before", (Py_ssize_t)result->surfaces_before,
        "surfaces_after", (Py_ssize_t)result->surfaces_after,
        "depth_before", (Py_ssize_t)result->depth_before,
        "depth_after", (Py_ssize_t)result->depth_after,
        "patterns_collected", (Py_ssize_t)result->patterns_collected,
        "candidates_proposed", (Py_ssize_t)result->candidates_proposed,
        "candidates_proven", (Py_ssize_t)result->candidates_proven,
        "candidates_disproven", (Py_ssize_t)result->candidates_disproven,
        "candidates_inconclusive", (Py_ssize_t)result->candidates_inconclusive,
        "proof_nodes", (Py_ssize_t)result->proof_nodes,
        "mixed_leaf_nodes", (Py_ssize_t)result->mixed_leaf_nodes,
        "parallel_batch_count", (Py_ssize_t)result->parallel_batch_count,
        "reserved_parallel_scratch_bytes",
            (unsigned long long)result->reserved_parallel_scratch_bytes);
    if (!witness || !stats) {
        Py_XDECREF(witness);
        Py_XDECREF(stats);
        return NULL;
    }
    return Py_BuildValue(
        "{s:O,s:O,s:O,s:O,s:I,s:(dddddd),s:s,s:O,s:s,s:s,s:N,s:n,s:n,s:n,s:N}",
        "changed", result->changed ? Py_True : Py_False,
        "applied", result->applied ? Py_True : Py_False,
        "proven_empty", result->proven_empty ? Py_True : Py_False,
        "complete", result->complete ? Py_True : Py_False,
        "root_node", (unsigned int)result->root_node_id,
        "bounds", result->bounds.min_x, result->bounds.max_x,
            result->bounds.min_y, result->bounds.max_y,
            result->bounds.min_z, result->bounds.max_z,
        "bounds_source", proof_bounds_source_name(result->bounds_source),
        "bounds_verified", result->bounds_verified ? Py_True : Py_False,
        "bounds_assumption", has_explicit_bounds
            ? "caller_asserted_complete_domain"
            : (result->bounds_verified
               ? "verified_complete_support" : "unresolved_complete_support"),
        "last_limit", proof_limit_name(result->last_limit),
        "witness", witness,
        "requested_workers", (Py_ssize_t)result->requested_workers,
        "actual_workers", (Py_ssize_t)result->actual_workers,
        "frontier_task_count", (Py_ssize_t)result->frontier_task_count,
        "stats", stats);
}

static PyObject* PyAleaSystem_simplify_cell_proven(
        PyAleaSystemObject* self, PyObject* args, PyObject* kwds) {
    Py_ssize_t cell_index, max_nodes = 250000, max_patterns = 64;
    Py_ssize_t max_candidates = 256, workers = 0;
    PyObject* bounds_obj = Py_None;
    int apply = 0, max_depth = 12;
    unsigned long long max_parallel_scratch_bytes = 64u * 1024u * 1024u;
    static char* kwlist[] = {
        "cell_index", "bounds", "apply", "max_depth", "max_nodes",
        "max_patterns", "max_candidates", "workers",
        "max_parallel_scratch_bytes", NULL
    };
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "n|OpinnnnK", kwlist,
            &cell_index, &bounds_obj, &apply, &max_depth, &max_nodes,
            &max_patterns, &max_candidates, &workers,
            &max_parallel_scratch_bytes)) return NULL;
    if (!self->sys) {
        PyErr_SetString(PyExc_RuntimeError, "System not initialized");
        return NULL;
    }
    if (cell_index < 0 || max_nodes <= 0 || max_patterns < 0 ||
        max_candidates <= 0 || workers < 0) {
        PyErr_SetString(PyExc_ValueError,
            "indices and limits must be non-negative, with positive node and candidate limits");
        return NULL;
    }

    alea_cell_simplify_proof_options_t options;
    alea_cell_simplify_proof_options_init(&options);
    options.apply = apply != 0;
    options.max_depth = max_depth;
    options.max_nodes = (size_t)max_nodes;
    options.max_patterns = (size_t)max_patterns;
    options.max_candidates = (size_t)max_candidates;
    options.requested_workers = (size_t)workers;
    options.max_parallel_scratch_bytes = (uint64_t)max_parallel_scratch_bytes;
    if (bounds_obj != Py_None) {
        if (parse_cell_volume_bounds(bounds_obj, &options.bounds) != 0) return NULL;
        options.has_bounds = true;
    }

    alea_cell_simplify_proof_result_t result;
    int rc;
    sighandler_func old_sigint = install_sigint();
    Py_BEGIN_ALLOW_THREADS
    rc = alea_cell_simplify_proven(self->sys, (size_t)cell_index,
                                   &options, &result);
    Py_END_ALLOW_THREADS
    if (restore_sigint(old_sigint)) return NULL;
    if (rc != 0) {
        PyErr_SetString(PyExc_RuntimeError, alea_error());
        return NULL;
    }

    return proof_simplification_result_dict(&result, options.has_bounds);
}

static PyObject* PyAleaSystem_simplify_cells_proven(
        PyAleaSystemObject* self, PyObject* args, PyObject* kwds) {
    PyObject* indices_obj;
    PyObject* bounds_obj = Py_None;
    Py_ssize_t max_nodes = 250000, max_patterns = 64;
    Py_ssize_t max_candidates = 256, workers = 0;
    int apply = 0, max_depth = 12;
    unsigned long long max_parallel_scratch_bytes = 64u * 1024u * 1024u;
    static char* kwlist[] = {
        "cell_indices", "bounds", "apply", "max_depth", "max_nodes",
        "max_patterns", "max_candidates", "workers",
        "max_parallel_scratch_bytes", NULL
    };
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "O|OpinnnnK", kwlist,
            &indices_obj, &bounds_obj, &apply, &max_depth, &max_nodes,
            &max_patterns, &max_candidates, &workers,
            &max_parallel_scratch_bytes)) return NULL;
    if (!self->sys) {
        PyErr_SetString(PyExc_RuntimeError, "System not initialized");
        return NULL;
    }
    if (max_depth < 0 || max_nodes <= 0 || max_patterns < 0 ||
        max_candidates <= 0 || workers < 0) {
        PyErr_SetString(PyExc_ValueError,
            "limits must be non-negative, with positive node and candidate limits");
        return NULL;
    }

    PyObject* indices = PySequence_Fast(
        indices_obj, "cell_indices must be a sequence");
    if (!indices) return NULL;
    Py_ssize_t count = PySequence_Fast_GET_SIZE(indices);
    PyObject* bounds = NULL;
    if (bounds_obj != Py_None) {
        bounds = PySequence_Fast(bounds_obj, "bounds must be a sequence");
        if (!bounds) { Py_DECREF(indices); return NULL; }
        if (PySequence_Fast_GET_SIZE(bounds) != count) {
            PyErr_SetString(PyExc_ValueError,
                            "bounds must have one entry per cell index");
            Py_DECREF(bounds);
            Py_DECREF(indices);
            return NULL;
        }
    }

    alea_cell_simplify_request_t* requests = count
        ? PyMem_Calloc((size_t)count, sizeof(*requests)) : NULL;
    alea_cell_simplify_proof_result_t* results = count
        ? PyMem_Calloc((size_t)count, sizeof(*results)) : NULL;
    if (count && (!requests || !results)) {
        PyErr_NoMemory();
        PyMem_Free(requests); PyMem_Free(results);
        Py_XDECREF(bounds); Py_DECREF(indices);
        return NULL;
    }
    for (Py_ssize_t i = 0; i < count; i++) {
        Py_ssize_t index = PyLong_AsSsize_t(PySequence_Fast_GET_ITEM(indices, i));
        if (index < 0 || PyErr_Occurred()) {
            if (!PyErr_Occurred())
                PyErr_SetString(PyExc_ValueError, "cell indices must be non-negative");
            PyMem_Free(requests); PyMem_Free(results);
            Py_XDECREF(bounds); Py_DECREF(indices);
            return NULL;
        }
        requests[i].cell_index = (size_t)index;
        if (bounds) {
            PyObject* item = PySequence_Fast_GET_ITEM(bounds, i);
            if (item != Py_None) {
                if (parse_cell_volume_bounds(item, &requests[i].bounds) != 0) {
                    PyMem_Free(requests); PyMem_Free(results);
                    Py_DECREF(bounds); Py_DECREF(indices);
                    return NULL;
                }
                requests[i].has_bounds = true;
            }
        }
    }
    Py_XDECREF(bounds);
    Py_DECREF(indices);

    alea_cells_simplify_proof_options_t options;
    alea_cells_simplify_proof_options_init(&options);
    options.apply = apply != 0;
    options.max_depth = max_depth;
    options.max_nodes_per_cell = (size_t)max_nodes;
    options.max_patterns_per_cell = (size_t)max_patterns;
    options.max_candidates_per_cell = (size_t)max_candidates;
    options.requested_workers = (size_t)workers;
    options.max_parallel_scratch_bytes = (uint64_t)max_parallel_scratch_bytes;
    alea_cells_simplify_proof_summary_t summary;
    int rc;
    sighandler_func old_sigint = install_sigint();
    Py_BEGIN_ALLOW_THREADS
    rc = alea_cells_simplify_proven(self->sys, requests, (size_t)count,
                                    &options, results, &summary);
    Py_END_ALLOW_THREADS
    if (restore_sigint(old_sigint)) {
        PyMem_Free(requests); PyMem_Free(results);
        return NULL;
    }
    if (rc != 0) {
        PyErr_SetString(PyExc_RuntimeError, alea_error());
        PyMem_Free(requests); PyMem_Free(results);
        return NULL;
    }

    PyObject* result_list = PyList_New(count);
    if (!result_list) {
        PyMem_Free(requests); PyMem_Free(results);
        return NULL;
    }
    for (Py_ssize_t i = 0; i < count; i++) {
        PyObject* item = proof_simplification_result_dict(
            &results[i], requests[i].has_bounds);
        if (!item) {
            Py_DECREF(result_list);
            PyMem_Free(requests); PyMem_Free(results);
            return NULL;
        }
        PyList_SET_ITEM(result_list, i, item);
    }
    PyObject* summary_dict = Py_BuildValue(
        "{s:n,s:n,s:n,s:n,s:n,s:n,s:n,s:n,s:n,s:K}",
        "selected_cells", (Py_ssize_t)summary.selected_cells,
        "changed_cells", (Py_ssize_t)summary.changed_cells,
        "applied_cells", (Py_ssize_t)summary.applied_cells,
        "proven_empty_cells", (Py_ssize_t)summary.proven_empty_cells,
        "complete_cells", (Py_ssize_t)summary.complete_cells,
        "inconclusive_cells", (Py_ssize_t)summary.inconclusive_cells,
        "requested_workers", (Py_ssize_t)summary.requested_workers,
        "actual_workers", (Py_ssize_t)summary.actual_workers,
        "parallel_batch_count", (Py_ssize_t)summary.parallel_batch_count,
        "reserved_parallel_scratch_bytes",
            (unsigned long long)summary.reserved_parallel_scratch_bytes);
    PyMem_Free(requests);
    PyMem_Free(results);
    if (!summary_dict) { Py_DECREF(result_list); return NULL; }
    return Py_BuildValue("{s:N,s:N}", "summary", summary_dict,
                         "results", result_list);
}

static PyObject* PyAleaSystem_carve_universe(PyAleaSystemObject* self,
                                                 PyObject* args,
                                                 PyObject* kwds) {
    int universe_id;
    unsigned int carve_root;
    int simplify = 1;
    int cell_limit = -1;
    static char* kwlist[] = {"universe_id", "carve_root", "simplify", "cell_limit", NULL};
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "iI|pi", kwlist,
                                     &universe_id, &carve_root, &simplify,
                                     &cell_limit)) {
        return NULL;
    }
    if (!self->sys) {
        PyErr_SetString(PyExc_RuntimeError, "System not initialized");
        return NULL;
    }

    int modified = 0, removed = 0;
    if (alea_carve_universe(self->sys, universe_id,
                            (alea_node_id_t)carve_root, simplify, cell_limit,
                            &modified, &removed) < 0) {
        PyErr_SetString(PyExc_RuntimeError, alea_error());
        return NULL;
    }
    return Py_BuildValue("{s:i,s:i}", "modified", modified,
                         "removed", removed);
}

/* ============================================================================
 * PyAleaSystem Methods - Numerical BBox Tightening
 * ============================================================================ */

static PyObject* PyAleaSystem_tighten_cell_bbox_numerical(PyAleaSystemObject* self, PyObject* args) {
    int cell_index;
    if (!PyArg_ParseTuple(args, "i", &cell_index)) return NULL;
    if (!self->sys) {
        PyErr_SetString(PyExc_RuntimeError, "System not initialized");
        return NULL;
    }

    int result;
    sighandler_func old_sigint = install_sigint();
    Py_BEGIN_ALLOW_THREADS
    result = alea_tighten_cell_bbox_numerical(self->sys, cell_index);
    Py_END_ALLOW_THREADS
    if (restore_sigint(old_sigint)) return NULL;

    if (result < 0) {
        PyErr_SetString(PyExc_RuntimeError, alea_error());
        return NULL;
    }

    Py_RETURN_NONE;
}

/* ============================================================================
 * PyAleaSystem Methods - Primitive Data
 * ============================================================================ */

static PyObject* PyAleaSystem_node_primitive_data(PyAleaSystemObject* self, PyObject* args) {
    unsigned long node_id;
    if (!PyArg_ParseTuple(args, "k", &node_id)) return NULL;
    if (!self->sys) {
        PyErr_SetString(PyExc_RuntimeError, "System not initialized");
        return NULL;
    }

    alea_primitive_data_t data;
    int rc = alea_node_primitive_data(self->sys, (alea_node_id_t)node_id, &data);
    if (rc < 0) {
        PyErr_SetString(PyExc_ValueError, "Node is not a primitive or invalid");
        return NULL;
    }

    alea_primitive_type_t ptype = alea_node_primitive_type(self->sys, (alea_node_id_t)node_id);

    switch (ptype) {
    case ALEA_PRIMITIVE_PLANE:
        return Py_BuildValue("{s:i, s:d, s:d, s:d, s:d}",
            "type", (int)ptype,
            "a", data.plane.a, "b", data.plane.b,
            "c", data.plane.c, "d", data.plane.d);

    case ALEA_PRIMITIVE_SPHERE:
        return Py_BuildValue("{s:i, s:d, s:d, s:d, s:d}",
            "type", (int)ptype,
            "center_x", data.sphere.center_x,
            "center_y", data.sphere.center_y,
            "center_z", data.sphere.center_z,
            "radius", data.sphere.radius);

    case ALEA_PRIMITIVE_CYLINDER_X:
        return Py_BuildValue("{s:i, s:d, s:d, s:d}",
            "type", (int)ptype,
            "center_y", data.cyl_x.center_y,
            "center_z", data.cyl_x.center_z,
            "radius", data.cyl_x.radius);

    case ALEA_PRIMITIVE_CYLINDER_Y:
        return Py_BuildValue("{s:i, s:d, s:d, s:d}",
            "type", (int)ptype,
            "center_x", data.cyl_y.center_x,
            "center_z", data.cyl_y.center_z,
            "radius", data.cyl_y.radius);

    case ALEA_PRIMITIVE_CYLINDER_Z:
        return Py_BuildValue("{s:i, s:d, s:d, s:d}",
            "type", (int)ptype,
            "center_x", data.cyl_z.center_x,
            "center_y", data.cyl_z.center_y,
            "radius", data.cyl_z.radius);

    case ALEA_PRIMITIVE_CONE_X:
        return Py_BuildValue("{s:i, s:d, s:d, s:d, s:d, s:i}",
            "type", (int)ptype,
            "apex_x", data.cone_x.apex_x,
            "apex_y", data.cone_x.apex_y,
            "apex_z", data.cone_x.apex_z,
            "tan_angle_sq", data.cone_x.tan_angle_sq,
            "sheet_selection", data.cone_x.sheet_selection);

    case ALEA_PRIMITIVE_CONE_Y:
        return Py_BuildValue("{s:i, s:d, s:d, s:d, s:d, s:i}",
            "type", (int)ptype,
            "apex_x", data.cone_y.apex_x,
            "apex_y", data.cone_y.apex_y,
            "apex_z", data.cone_y.apex_z,
            "tan_angle_sq", data.cone_y.tan_angle_sq,
            "sheet_selection", data.cone_y.sheet_selection);

    case ALEA_PRIMITIVE_CONE_Z:
        return Py_BuildValue("{s:i, s:d, s:d, s:d, s:d, s:i}",
            "type", (int)ptype,
            "apex_x", data.cone_z.apex_x,
            "apex_y", data.cone_z.apex_y,
            "apex_z", data.cone_z.apex_z,
            "tan_angle_sq", data.cone_z.tan_angle_sq,
            "sheet_selection", data.cone_z.sheet_selection);

    case ALEA_PRIMITIVE_RPP:
        return Py_BuildValue("{s:i, s:d, s:d, s:d, s:d, s:d, s:d}",
            "type", (int)ptype,
            "min_x", data.box.min_x, "max_x", data.box.max_x,
            "min_y", data.box.min_y, "max_y", data.box.max_y,
            "min_z", data.box.min_z, "max_z", data.box.max_z);

    case ALEA_PRIMITIVE_QUADRIC: {
        PyObject* coeffs = PyList_New(10);
        if (!coeffs) return NULL;
        for (int i = 0; i < 10; i++) {
            PyList_SET_ITEM(coeffs, i, PyFloat_FromDouble(data.quadric.coeffs[i]));
        }
        PyObject* result = Py_BuildValue("{s:i, s:N}",
            "type", (int)ptype, "coeffs", coeffs);
        return result;
    }

    case ALEA_PRIMITIVE_TORUS_X:
    case ALEA_PRIMITIVE_TORUS_Y:
    case ALEA_PRIMITIVE_TORUS_Z:
        return Py_BuildValue("{s:i, s:i, s:d, s:d, s:d, s:d, s:d, s:d}",
            "type", (int)ptype,
            "axis", (int)data.torus.axis,
            "center_x", data.torus.center_x,
            "center_y", data.torus.center_y,
            "center_z", data.torus.center_z,
            "major_radius", data.torus.major_radius,
            "minor_radius", data.torus.minor_radius,
            "axial_semiwidth_B", data.torus.axial_semiwidth_B);

    case ALEA_PRIMITIVE_RCC:
        return Py_BuildValue("{s:i, s:d, s:d, s:d, s:d, s:d, s:d, s:d}",
            "type", (int)ptype,
            "base_x", data.rcc.base_x, "base_y", data.rcc.base_y, "base_z", data.rcc.base_z,
            "height_x", data.rcc.height_x, "height_y", data.rcc.height_y, "height_z", data.rcc.height_z,
            "radius", data.rcc.radius);

    case ALEA_PRIMITIVE_BOX:
        return Py_BuildValue("{s:i, s:d, s:d, s:d, s:d, s:d, s:d, s:d, s:d, s:d, s:d, s:d, s:d}",
            "type", (int)ptype,
            "corner_x", data.box_general.corner_x,
            "corner_y", data.box_general.corner_y,
            "corner_z", data.box_general.corner_z,
            "v1_x", data.box_general.v1_x, "v1_y", data.box_general.v1_y, "v1_z", data.box_general.v1_z,
            "v2_x", data.box_general.v2_x, "v2_y", data.box_general.v2_y, "v2_z", data.box_general.v2_z,
            "v3_x", data.box_general.v3_x, "v3_y", data.box_general.v3_y, "v3_z", data.box_general.v3_z);

    case ALEA_PRIMITIVE_SPH:
        return Py_BuildValue("{s:i, s:d, s:d, s:d, s:d}",
            "type", (int)ptype,
            "center_x", data.sph.center_x,
            "center_y", data.sph.center_y,
            "center_z", data.sph.center_z,
            "radius", data.sph.radius);

    case ALEA_PRIMITIVE_TRC:
        return Py_BuildValue("{s:i, s:d, s:d, s:d, s:d, s:d, s:d, s:d, s:d}",
            "type", (int)ptype,
            "base_x", data.trc.base_x, "base_y", data.trc.base_y, "base_z", data.trc.base_z,
            "height_x", data.trc.height_x, "height_y", data.trc.height_y, "height_z", data.trc.height_z,
            "base_radius", data.trc.base_radius, "top_radius", data.trc.top_radius);

    case ALEA_PRIMITIVE_ELL:
        return Py_BuildValue("{s:i, s:d, s:d, s:d, s:d, s:d, s:d, s:d}",
            "type", (int)ptype,
            "v1_x", data.ell.v1_x, "v1_y", data.ell.v1_y, "v1_z", data.ell.v1_z,
            "v2_x", data.ell.v2_x, "v2_y", data.ell.v2_y, "v2_z", data.ell.v2_z,
            "major_axis_len", data.ell.major_axis_len);

    case ALEA_PRIMITIVE_REC:
        return Py_BuildValue("{s:i, s:d, s:d, s:d, s:d, s:d, s:d, s:d, s:d, s:d, s:d, s:d, s:d}",
            "type", (int)ptype,
            "base_x", data.rec.base_x, "base_y", data.rec.base_y, "base_z", data.rec.base_z,
            "height_x", data.rec.height_x, "height_y", data.rec.height_y, "height_z", data.rec.height_z,
            "axis1_x", data.rec.axis1_x, "axis1_y", data.rec.axis1_y, "axis1_z", data.rec.axis1_z,
            "axis2_x", data.rec.axis2_x, "axis2_y", data.rec.axis2_y, "axis2_z", data.rec.axis2_z);

    case ALEA_PRIMITIVE_WED:
        return Py_BuildValue("{s:i, s:d, s:d, s:d, s:d, s:d, s:d, s:d, s:d, s:d, s:d, s:d, s:d}",
            "type", (int)ptype,
            "vertex_x", data.wed.vertex_x, "vertex_y", data.wed.vertex_y, "vertex_z", data.wed.vertex_z,
            "v1_x", data.wed.v1_x, "v1_y", data.wed.v1_y, "v1_z", data.wed.v1_z,
            "v2_x", data.wed.v2_x, "v2_y", data.wed.v2_y, "v2_z", data.wed.v2_z,
            "v3_x", data.wed.v3_x, "v3_y", data.wed.v3_y, "v3_z", data.wed.v3_z);

    case ALEA_PRIMITIVE_RHP:
        return Py_BuildValue("{s:i, s:d, s:d, s:d, s:d, s:d, s:d, s:d, s:d, s:d, s:d, s:d, s:d, s:d, s:d, s:d}",
            "type", (int)ptype,
            "base_x", data.rhp.base_x, "base_y", data.rhp.base_y, "base_z", data.rhp.base_z,
            "height_x", data.rhp.height_x, "height_y", data.rhp.height_y, "height_z", data.rhp.height_z,
            "r1_x", data.rhp.r1_x, "r1_y", data.rhp.r1_y, "r1_z", data.rhp.r1_z,
            "r2_x", data.rhp.r2_x, "r2_y", data.rhp.r2_y, "r2_z", data.rhp.r2_z,
            "r3_x", data.rhp.r3_x, "r3_y", data.rhp.r3_y, "r3_z", data.rhp.r3_z);

    case ALEA_PRIMITIVE_ARB: {
        PyObject* corners = PyList_New(data.arb.num_corners);
        if (!corners) return NULL;
        for (int i = 0; i < data.arb.num_corners; i++) {
            PyObject* pt = Py_BuildValue("(ddd)",
                data.arb.corners[i][0], data.arb.corners[i][1], data.arb.corners[i][2]);
            if (!pt) { Py_DECREF(corners); return NULL; }
            PyList_SET_ITEM(corners, i, pt);
        }
        PyObject* faces = PyList_New(data.arb.num_faces);
        if (!faces) { Py_DECREF(corners); return NULL; }
        for (int i = 0; i < data.arb.num_faces; i++) {
            PyObject* face = Py_BuildValue("(iiii)",
                data.arb.faces[i][0], data.arb.faces[i][1],
                data.arb.faces[i][2], data.arb.faces[i][3]);
            if (!face) { Py_DECREF(corners); Py_DECREF(faces); return NULL; }
            PyList_SET_ITEM(faces, i, face);
        }
        PyObject* result = Py_BuildValue("{s:i, s:N, s:N, s:i, s:i}",
            "type", (int)ptype,
            "corners", corners, "faces", faces,
            "num_corners", data.arb.num_corners,
            "num_faces", data.arb.num_faces);
        return result;
    }

    default:
        return Py_BuildValue("{s:i}", "type", (int)ptype);
    }
}

/* ============================================================================
 * PyAleaSystem Methods - Cell Expression
 * ============================================================================ */

static PyObject* PyAleaSystem_cell_expr(PyAleaSystemObject* self, PyObject* args, PyObject* kwds) {
    size_t cell_index;
    const char* union_op = ":";
    const char* inter_op = " ";
    const char* compl_op = "#";
    static char* kwlist[] = {"cell_index", "union_op", "inter_op", "compl_op", NULL};

    if (!PyArg_ParseTupleAndKeywords(args, kwds, "n|sss", kwlist,
            &cell_index, &union_op, &inter_op, &compl_op)) {
        return NULL;
    }

    if (!self->sys) {
        PyErr_SetString(PyExc_RuntimeError, "System not initialized");
        return NULL;
    }

    char* expr = alea_cell_expr(self->sys, cell_index, union_op, inter_op, compl_op);
    if (!expr) {
        PyErr_Format(PyExc_IndexError, "Cell index %zu out of range", cell_index);
        return NULL;
    }

    PyObject* result = PyUnicode_FromString(expr);
    free(expr);
    return result;
}

static PyObject* PyAleaSystem_node_tree(PyAleaSystemObject* self, PyObject* args) {
    unsigned long node_id;
    if (!PyArg_ParseTuple(args, "k", &node_id)) return NULL;
    if (!self->sys) {
        PyErr_SetString(PyExc_RuntimeError, "System not initialized");
        return NULL;
    }
    return build_node_tree(self->sys, (alea_node_id_t)node_id);
}

/* A bounded CSG-leaf receipt for interactive explanations. */
static PyObject* PyAleaSystem_cell_surface_ids(
    PyAleaSystemObject* self, PyObject* args, PyObject* kwds) {
    Py_ssize_t cell_index, max_surfaces = 32, max_nodes = 256;
    static char* kwlist[] = {"cell_index", "max_surfaces", "max_nodes", NULL};
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "n|nn", kwlist,
                                     &cell_index, &max_surfaces, &max_nodes)) return NULL;
    if (!self->sys) { PyErr_SetString(PyExc_RuntimeError, "System not initialized"); return NULL; }
    if (cell_index < 0 || max_surfaces <= 0 || max_nodes <= 0 ||
        max_surfaces > 4096 || max_nodes > 65536) {
        PyErr_SetString(PyExc_ValueError, "invalid bounded cell-surface limits");
        return NULL;
    }
    alea_cell_info_t info;
    if (alea_cell_get_info(self->sys, (size_t)cell_index, &info) < 0) {
        PyErr_Format(PyExc_IndexError, "Cell index %zd out of range", cell_index);
        return NULL;
    }
    PyObject* result = PyDict_New();
    PyObject* ids = PyList_New(0);
    if (!result || !ids) { Py_XDECREF(result); Py_XDECREF(ids); return NULL; }
    alea_node_id_t* stack = PyMem_Malloc((size_t)max_nodes * sizeof(*stack));
    int* seen = PyMem_Malloc((size_t)max_surfaces * sizeof(*seen));
    if (!stack || !seen) {
        PyMem_Free(stack); PyMem_Free(seen); Py_DECREF(ids); Py_DECREF(result);
        return PyErr_NoMemory();
    }
    Py_ssize_t top = 0, count = 0, visited = 0;
    int truncated = 0;
    if (info.root != ALEA_NODE_ID_INVALID) stack[top++] = info.root;
    while (top > 0) {
        if (visited >= max_nodes) { truncated = 1; break; }
        alea_node_id_t node = stack[--top];
        visited++;
        alea_operation_t operation = alea_node_operation(self->sys, node);
        if (operation == ALEA_OP_PRIMITIVE) {
            int surface_id = alea_node_surface_id(self->sys, node);
            int duplicate = 0;
            for (Py_ssize_t i = 0; i < count; i++)
                if (seen[i] == surface_id) { duplicate = 1; break; }
            if (!duplicate) {
                if (count >= max_surfaces) { truncated = 1; break; }
                seen[count++] = surface_id;
                PyObject* value = PyLong_FromLong(surface_id);
                if (!value || PyList_Append(ids, value) < 0) {
                    Py_XDECREF(value); PyMem_Free(stack); PyMem_Free(seen);
                    Py_DECREF(ids); Py_DECREF(result); return NULL;
                }
                Py_DECREF(value);
            }
            continue;
        }
        int children = operation == ALEA_OP_COMPLEMENT ? 1 : 2;
        if (top + children > max_nodes) { truncated = 1; break; }
        stack[top++] = alea_node_left(self->sys, node);
        if (children == 2) stack[top++] = alea_node_right(self->sys, node);
    }
    PyMem_Free(stack); PyMem_Free(seen);
    PyDict_SetItemString(result, "surface_ids", ids);
    dict_set_new(result, "visited_nodes", PyLong_FromSsize_t(visited));
    dict_set_new(result, "truncated", PyBool_FromLong(truncated));
    Py_DECREF(ids);
    return result;
}
