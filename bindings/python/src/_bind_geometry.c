// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/* This file is part of pyalea_binding.c — do NOT compile separately.
 * It is #included from pyalea_binding.c to keep a single translation unit.
 *
 * Contents: Primitive creation, boolean operations, surface creation,
 *           cell registration (add_cell).
 */

/* ============================================================================
 * PyAleaSystem Methods - Primitive Creation
 * ============================================================================ */

static PyObject* PyAleaSystem_create_plane(PyAleaSystemObject* self, PyObject* args) {
    double a, b, c, d;
    int sense;
    if (!PyArg_ParseTuple(args, "ddddi", &a, &b, &c, &d, &sense)) return NULL;
    if (!self->sys) { PyErr_SetString(PyExc_RuntimeError, "System not initialized"); return NULL; }

    int idx = alea_plane_surface(self->sys, 0, a, b, c, d);
    if (idx < 0) {
        PyErr_SetString(PyExc_RuntimeError, "Failed to create plane");
        return NULL;
    }
    alea_node_id_t node = alea_halfspace(self->sys, idx, sense);
    if (node == ALEA_NODE_ID_INVALID) {
        PyErr_SetString(PyExc_RuntimeError, "Failed to create plane half-space");
        return NULL;
    }
    return PyLong_FromUnsignedLong(node);
}

static PyObject* PyAleaSystem_create_sphere(PyAleaSystemObject* self, PyObject* args) {
    double cx, cy, cz, radius;
    int sense;
    if (!PyArg_ParseTuple(args, "ddddi", &cx, &cy, &cz, &radius, &sense)) return NULL;
    if (!self->sys) { PyErr_SetString(PyExc_RuntimeError, "System not initialized"); return NULL; }

    int idx = alea_sphere_surface(self->sys, 0, cx, cy, cz, radius);
    if (idx < 0) {
        PyErr_SetString(PyExc_RuntimeError, "Failed to create sphere");
        return NULL;
    }
    alea_node_id_t node = alea_halfspace(self->sys, idx, sense);
    if (node == ALEA_NODE_ID_INVALID) {
        PyErr_SetString(PyExc_RuntimeError, "Failed to create sphere half-space");
        return NULL;
    }
    return PyLong_FromUnsignedLong(node);
}

static PyObject* PyAleaSystem_create_box(PyAleaSystemObject* self, PyObject* args) {
    double xmin, xmax, ymin, ymax, zmin, zmax;
    int sense;
    if (!PyArg_ParseTuple(args, "ddddddi", &xmin, &xmax, &ymin, &ymax, &zmin, &zmax, &sense)) return NULL;
    if (!self->sys) { PyErr_SetString(PyExc_RuntimeError, "System not initialized"); return NULL; }

    int idx = alea_box_surface(self->sys, 0, xmin, xmax, ymin, ymax, zmin, zmax);
    if (idx < 0) {
        PyErr_SetString(PyExc_RuntimeError, "Failed to create box");
        return NULL;
    }
    alea_node_id_t node = alea_halfspace(self->sys, idx, sense);
    if (node == ALEA_NODE_ID_INVALID) {
        PyErr_SetString(PyExc_RuntimeError, "Failed to create box half-space");
        return NULL;
    }
    return PyLong_FromUnsignedLong(node);
}

static PyObject* PyAleaSystem_create_cylinder_z(PyAleaSystemObject* self, PyObject* args) {
    double cx, cy, radius;
    int sense;
    if (!PyArg_ParseTuple(args, "dddi", &cx, &cy, &radius, &sense)) return NULL;
    if (!self->sys) { PyErr_SetString(PyExc_RuntimeError, "System not initialized"); return NULL; }

    int idx = alea_cylinder_z_surface(self->sys, 0, cx, cy, radius);
    if (idx < 0) {
        PyErr_SetString(PyExc_RuntimeError, "Failed to create cylinder");
        return NULL;
    }
    alea_node_id_t node = alea_halfspace(self->sys, idx, sense);
    if (node == ALEA_NODE_ID_INVALID) {
        PyErr_SetString(PyExc_RuntimeError, "Failed to create cylinder half-space");
        return NULL;
    }
    return PyLong_FromUnsignedLong(node);
}

/* ============================================================================
 * PyAleaSystem Methods - Boolean Operations
 * ============================================================================ */

static PyObject* PyAleaSystem_create_union(PyAleaSystemObject* self, PyObject* args) {
    unsigned long a, b;
    if (!PyArg_ParseTuple(args, "kk", &a, &b)) return NULL;
    if (!self->sys) { PyErr_SetString(PyExc_RuntimeError, "System not initialized"); return NULL; }

    alea_node_id_t node = alea_union(self->sys, (alea_node_id_t)a, (alea_node_id_t)b);
    if (node == ALEA_NODE_ID_INVALID) {
        PyErr_SetString(PyExc_RuntimeError, "Failed to create union");
        return NULL;
    }
    return PyLong_FromUnsignedLong(node);
}

static PyObject* PyAleaSystem_create_intersection(PyAleaSystemObject* self, PyObject* args) {
    unsigned long a, b;
    if (!PyArg_ParseTuple(args, "kk", &a, &b)) return NULL;
    if (!self->sys) { PyErr_SetString(PyExc_RuntimeError, "System not initialized"); return NULL; }

    alea_node_id_t node = alea_intersection(self->sys, (alea_node_id_t)a, (alea_node_id_t)b);
    if (node == ALEA_NODE_ID_INVALID) {
        PyErr_SetString(PyExc_RuntimeError, "Failed to create intersection");
        return NULL;
    }
    return PyLong_FromUnsignedLong(node);
}

static PyObject* PyAleaSystem_create_difference(PyAleaSystemObject* self, PyObject* args) {
    unsigned long a, b;
    if (!PyArg_ParseTuple(args, "kk", &a, &b)) return NULL;
    if (!self->sys) { PyErr_SetString(PyExc_RuntimeError, "System not initialized"); return NULL; }

    alea_node_id_t node = alea_difference(self->sys, (alea_node_id_t)a, (alea_node_id_t)b);
    if (node == ALEA_NODE_ID_INVALID) {
        PyErr_SetString(PyExc_RuntimeError, "Failed to create difference");
        return NULL;
    }
    return PyLong_FromUnsignedLong(node);
}

static PyObject* PyAleaSystem_create_complement(PyAleaSystemObject* self, PyObject* args) {
    unsigned long a;
    if (!PyArg_ParseTuple(args, "k", &a)) return NULL;
    if (!self->sys) { PyErr_SetString(PyExc_RuntimeError, "System not initialized"); return NULL; }

    alea_node_id_t node = alea_complement(self->sys, (alea_node_id_t)a);
    if (node == ALEA_NODE_ID_INVALID) {
        PyErr_SetString(PyExc_RuntimeError, "Failed to create complement");
        return NULL;
    }
    return PyLong_FromUnsignedLong(node);
}

static PyObject* PyAleaSystem_create_union_many(PyAleaSystemObject* self, PyObject* args) {
    PyObject* nodes_list;
    if (!PyArg_ParseTuple(args, "O!", &PyList_Type, &nodes_list)) return NULL;
    if (!self->sys) { PyErr_SetString(PyExc_RuntimeError, "System not initialized"); return NULL; }

    Py_ssize_t count = PyList_Size(nodes_list);
    if (count < 1) {
        PyErr_SetString(PyExc_ValueError, "Need at least one node");
        return NULL;
    }

    alea_node_id_t* nodes = malloc(count * sizeof(alea_node_id_t));
    if (!nodes) {
        PyErr_NoMemory();
        return NULL;
    }

    for (Py_ssize_t i = 0; i < count; i++) {
        PyObject* item = PyList_GetItem(nodes_list, i);
        if (!PyLong_Check(item)) {
            free(nodes);
            PyErr_SetString(PyExc_TypeError, "All items must be integers (node IDs)");
            return NULL;
        }
        nodes[i] = (alea_node_id_t)PyLong_AsUnsignedLong(item);
    }

    alea_node_id_t result = alea_union_n(self->sys, nodes, count);
    free(nodes);

    if (result == ALEA_NODE_ID_INVALID) {
        PyErr_SetString(PyExc_RuntimeError, "Failed to create union");
        return NULL;
    }
    return PyLong_FromUnsignedLong(result);
}

static PyObject* PyAleaSystem_create_intersection_many(PyAleaSystemObject* self, PyObject* args) {
    PyObject* nodes_list;
    if (!PyArg_ParseTuple(args, "O!", &PyList_Type, &nodes_list)) return NULL;
    if (!self->sys) { PyErr_SetString(PyExc_RuntimeError, "System not initialized"); return NULL; }

    Py_ssize_t count = PyList_Size(nodes_list);
    if (count < 1) {
        PyErr_SetString(PyExc_ValueError, "Need at least one node");
        return NULL;
    }

    alea_node_id_t* nodes = malloc(count * sizeof(alea_node_id_t));
    if (!nodes) {
        PyErr_NoMemory();
        return NULL;
    }

    for (Py_ssize_t i = 0; i < count; i++) {
        PyObject* item = PyList_GetItem(nodes_list, i);
        if (!PyLong_Check(item)) {
            free(nodes);
            PyErr_SetString(PyExc_TypeError, "All items must be integers (node IDs)");
            return NULL;
        }
        nodes[i] = (alea_node_id_t)PyLong_AsUnsignedLong(item);
    }

    alea_node_id_t result = alea_intersection_n(self->sys, nodes, count);
    free(nodes);

    if (result == ALEA_NODE_ID_INVALID) {
        PyErr_SetString(PyExc_RuntimeError, "Failed to create intersection");
        return NULL;
    }
    return PyLong_FromUnsignedLong(result);
}

/* ============================================================================
 * PyAleaSystem Methods - Surface Creation (with automatic registration)
 *
 * These functions use alea_*_surface() which creates surfaces with both
 * positive and negative halfspace nodes, properly registered for raycast.
 * ============================================================================ */

static PyObject* PyAleaSystem_sphere_surface(PyAleaSystemObject* self, PyObject* args) {
    int surface_id;
    double cx, cy, cz, radius;
    if (!PyArg_ParseTuple(args, "idddd", &surface_id, &cx, &cy, &cz, &radius)) return NULL;
    if (!self->sys) { PyErr_SetString(PyExc_RuntimeError, "System not initialized"); return NULL; }

    int idx = alea_sphere_surface(self->sys, surface_id, cx, cy, cz, radius);
    if (idx < 0) {
        PyErr_SetString(PyExc_RuntimeError, "Failed to create sphere surface");
        return NULL;
    }

    alea_node_id_t pos_node, neg_node;
    alea_surface_get(self->sys, idx, NULL, NULL, &pos_node, &neg_node, NULL);
    return Py_BuildValue("(ikk)", idx, (unsigned long)pos_node, (unsigned long)neg_node);
}

static PyObject* PyAleaSystem_cylinder_z_surface(PyAleaSystemObject* self, PyObject* args) {
    int surface_id;
    double cx, cy, radius;
    if (!PyArg_ParseTuple(args, "iddd", &surface_id, &cx, &cy, &radius)) return NULL;
    if (!self->sys) { PyErr_SetString(PyExc_RuntimeError, "System not initialized"); return NULL; }

    int idx = alea_cylinder_z_surface(self->sys, surface_id, cx, cy, radius);
    if (idx < 0) {
        PyErr_SetString(PyExc_RuntimeError, "Failed to create cylinder surface");
        return NULL;
    }

    alea_node_id_t pos_node, neg_node;
    alea_surface_get(self->sys, idx, NULL, NULL, &pos_node, &neg_node, NULL);
    return Py_BuildValue("(ikk)", idx, (unsigned long)pos_node, (unsigned long)neg_node);
}

static PyObject* PyAleaSystem_box_surface(PyAleaSystemObject* self, PyObject* args) {
    int surface_id;
    double xmin, xmax, ymin, ymax, zmin, zmax;
    if (!PyArg_ParseTuple(args, "idddddd", &surface_id, &xmin, &xmax, &ymin, &ymax, &zmin, &zmax)) return NULL;
    if (!self->sys) { PyErr_SetString(PyExc_RuntimeError, "System not initialized"); return NULL; }

    int idx = alea_box_surface(self->sys, surface_id, xmin, xmax, ymin, ymax, zmin, zmax);
    if (idx < 0) {
        PyErr_SetString(PyExc_RuntimeError, "Failed to create box surface");
        return NULL;
    }

    alea_node_id_t pos_node, neg_node;
    alea_surface_get(self->sys, idx, NULL, NULL, &pos_node, &neg_node, NULL);
    return Py_BuildValue("(ikk)", idx, (unsigned long)pos_node, (unsigned long)neg_node);
}

static PyObject* PyAleaSystem_plane_surface(PyAleaSystemObject* self, PyObject* args) {
    int surface_id;
    double a, b, c, d;
    if (!PyArg_ParseTuple(args, "idddd", &surface_id, &a, &b, &c, &d)) return NULL;
    if (!self->sys) { PyErr_SetString(PyExc_RuntimeError, "System not initialized"); return NULL; }

    int idx = alea_plane_surface(self->sys, surface_id, a, b, c, d);
    if (idx < 0) {
        PyErr_SetString(PyExc_RuntimeError, "Failed to create plane surface");
        return NULL;
    }

    alea_node_id_t pos_node, neg_node;
    alea_surface_get(self->sys, idx, NULL, NULL, &pos_node, &neg_node, NULL);
    return Py_BuildValue("(ikk)", idx, (unsigned long)pos_node, (unsigned long)neg_node);
}

static PyObject* PyAleaSystem_cylinder_x_surface(PyAleaSystemObject* self, PyObject* args) {
    int surface_id;
    double cy, cz, radius;
    if (!PyArg_ParseTuple(args, "iddd", &surface_id, &cy, &cz, &radius)) return NULL;
    if (!self->sys) { PyErr_SetString(PyExc_RuntimeError, "System not initialized"); return NULL; }

    int idx = alea_cylinder_x_surface(self->sys, surface_id, cy, cz, radius);
    if (idx < 0) {
        PyErr_SetString(PyExc_RuntimeError, "Failed to create cylinder surface");
        return NULL;
    }

    alea_node_id_t pos_node, neg_node;
    alea_surface_get(self->sys, idx, NULL, NULL, &pos_node, &neg_node, NULL);
    return Py_BuildValue("(ikk)", idx, (unsigned long)pos_node, (unsigned long)neg_node);
}

static PyObject* PyAleaSystem_cylinder_y_surface(PyAleaSystemObject* self, PyObject* args) {
    int surface_id;
    double cx, cz, radius;
    if (!PyArg_ParseTuple(args, "iddd", &surface_id, &cx, &cz, &radius)) return NULL;
    if (!self->sys) { PyErr_SetString(PyExc_RuntimeError, "System not initialized"); return NULL; }

    int idx = alea_cylinder_y_surface(self->sys, surface_id, cx, cz, radius);
    if (idx < 0) {
        PyErr_SetString(PyExc_RuntimeError, "Failed to create cylinder surface");
        return NULL;
    }

    alea_node_id_t pos_node, neg_node;
    alea_surface_get(self->sys, idx, NULL, NULL, &pos_node, &neg_node, NULL);
    return Py_BuildValue("(ikk)", idx, (unsigned long)pos_node, (unsigned long)neg_node);
}

static PyObject* PyAleaSystem_cone_z_surface(PyAleaSystemObject* self, PyObject* args) {
    int surface_id;
    double cx, cy, cz, t_squared;
    if (!PyArg_ParseTuple(args, "idddd", &surface_id, &cx, &cy, &cz, &t_squared)) return NULL;
    if (!self->sys) { PyErr_SetString(PyExc_RuntimeError, "System not initialized"); return NULL; }

    int idx = alea_cone_z_surface(self->sys, surface_id, cx, cy, cz, t_squared);
    if (idx < 0) {
        PyErr_SetString(PyExc_RuntimeError, "Failed to create cone surface");
        return NULL;
    }

    alea_node_id_t pos_node, neg_node;
    alea_surface_get(self->sys, idx, NULL, NULL, &pos_node, &neg_node, NULL);
    return Py_BuildValue("(ikk)", idx, (unsigned long)pos_node, (unsigned long)neg_node);
}

static PyObject* PyAleaSystem_cone_x_surface(PyAleaSystemObject* self, PyObject* args) {
    int surface_id;
    double cx, cy, cz, t_squared;
    if (!PyArg_ParseTuple(args, "idddd", &surface_id, &cx, &cy, &cz, &t_squared)) return NULL;
    if (!self->sys) { PyErr_SetString(PyExc_RuntimeError, "System not initialized"); return NULL; }

    int idx = alea_cone_x_surface(self->sys, surface_id, cx, cy, cz, t_squared);
    if (idx < 0) {
        PyErr_SetString(PyExc_RuntimeError, "Failed to create cone surface");
        return NULL;
    }

    alea_node_id_t pos_node, neg_node;
    alea_surface_get(self->sys, idx, NULL, NULL, &pos_node, &neg_node, NULL);
    return Py_BuildValue("(ikk)", idx, (unsigned long)pos_node, (unsigned long)neg_node);
}

static PyObject* PyAleaSystem_cone_y_surface(PyAleaSystemObject* self, PyObject* args) {
    int surface_id;
    double cx, cy, cz, t_squared;
    if (!PyArg_ParseTuple(args, "idddd", &surface_id, &cx, &cy, &cz, &t_squared)) return NULL;
    if (!self->sys) { PyErr_SetString(PyExc_RuntimeError, "System not initialized"); return NULL; }

    int idx = alea_cone_y_surface(self->sys, surface_id, cx, cy, cz, t_squared);
    if (idx < 0) {
        PyErr_SetString(PyExc_RuntimeError, "Failed to create cone surface");
        return NULL;
    }

    alea_node_id_t pos_node, neg_node;
    alea_surface_get(self->sys, idx, NULL, NULL, &pos_node, &neg_node, NULL);
    return Py_BuildValue("(ikk)", idx, (unsigned long)pos_node, (unsigned long)neg_node);
}

static PyObject* PyAleaSystem_torus_z_surface(PyAleaSystemObject* self, PyObject* args) {
    int surface_id;
    double cx, cy, cz, major_radius, minor_radius;
    if (!PyArg_ParseTuple(args, "iddddd", &surface_id, &cx, &cy, &cz, &major_radius, &minor_radius)) return NULL;
    if (!self->sys) { PyErr_SetString(PyExc_RuntimeError, "System not initialized"); return NULL; }

    int idx = alea_torus_z_surface(self->sys, surface_id, cx, cy, cz, major_radius, minor_radius);
    if (idx < 0) {
        PyErr_SetString(PyExc_RuntimeError, "Failed to create torus surface");
        return NULL;
    }

    alea_node_id_t pos_node, neg_node;
    alea_surface_get(self->sys, idx, NULL, NULL, &pos_node, &neg_node, NULL);
    return Py_BuildValue("(ikk)", idx, (unsigned long)pos_node, (unsigned long)neg_node);
}

static PyObject* PyAleaSystem_torus_x_surface(PyAleaSystemObject* self, PyObject* args) {
    int surface_id;
    double cx, cy, cz, major_radius, minor_radius;
    if (!PyArg_ParseTuple(args, "iddddd", &surface_id, &cx, &cy, &cz, &major_radius, &minor_radius)) return NULL;
    if (!self->sys) { PyErr_SetString(PyExc_RuntimeError, "System not initialized"); return NULL; }

    int idx = alea_torus_x_surface(self->sys, surface_id, cx, cy, cz, major_radius, minor_radius);
    if (idx < 0) {
        PyErr_SetString(PyExc_RuntimeError, "Failed to create torus surface");
        return NULL;
    }

    alea_node_id_t pos_node, neg_node;
    alea_surface_get(self->sys, idx, NULL, NULL, &pos_node, &neg_node, NULL);
    return Py_BuildValue("(ikk)", idx, (unsigned long)pos_node, (unsigned long)neg_node);
}

static PyObject* PyAleaSystem_torus_y_surface(PyAleaSystemObject* self, PyObject* args) {
    int surface_id;
    double cx, cy, cz, major_radius, minor_radius;
    if (!PyArg_ParseTuple(args, "iddddd", &surface_id, &cx, &cy, &cz, &major_radius, &minor_radius)) return NULL;
    if (!self->sys) { PyErr_SetString(PyExc_RuntimeError, "System not initialized"); return NULL; }

    int idx = alea_torus_y_surface(self->sys, surface_id, cx, cy, cz, major_radius, minor_radius);
    if (idx < 0) {
        PyErr_SetString(PyExc_RuntimeError, "Failed to create torus surface");
        return NULL;
    }

    alea_node_id_t pos_node, neg_node;
    alea_surface_get(self->sys, idx, NULL, NULL, &pos_node, &neg_node, NULL);
    return Py_BuildValue("(ikk)", idx, (unsigned long)pos_node, (unsigned long)neg_node);
}

static PyObject* PyAleaSystem_quadric_surface(PyAleaSystemObject* self, PyObject* args) {
    int surface_id;
    double A, B, C, D, E, F, G, H, I, J;
    if (!PyArg_ParseTuple(args, "idddddddddd", &surface_id, &A, &B, &C, &D, &E, &F, &G, &H, &I, &J)) return NULL;
    if (!self->sys) { PyErr_SetString(PyExc_RuntimeError, "System not initialized"); return NULL; }

    int idx = alea_quadric_surface(self->sys, surface_id, A, B, C, D, E, F, G, H, I, J);
    if (idx < 0) {
        PyErr_SetString(PyExc_RuntimeError, "Failed to create quadric surface");
        return NULL;
    }

    alea_node_id_t pos_node, neg_node;
    alea_surface_get(self->sys, idx, NULL, NULL, &pos_node, &neg_node, NULL);
    return Py_BuildValue("(ikk)", idx, (unsigned long)pos_node, (unsigned long)neg_node);
}

static PyObject* PyAleaSystem_rcc_surface(PyAleaSystemObject* self, PyObject* args) {
    int surface_id;
    double base_x, base_y, base_z, height_x, height_y, height_z, radius;
    if (!PyArg_ParseTuple(args, "iddddddd", &surface_id, &base_x, &base_y, &base_z,
                          &height_x, &height_y, &height_z, &radius)) return NULL;
    if (!self->sys) { PyErr_SetString(PyExc_RuntimeError, "System not initialized"); return NULL; }

    int idx = alea_rcc_surface(self->sys, surface_id, base_x, base_y, base_z,
                                   height_x, height_y, height_z, radius);
    if (idx < 0) {
        PyErr_SetString(PyExc_RuntimeError, "Failed to create RCC surface");
        return NULL;
    }

    alea_node_id_t pos_node, neg_node;
    alea_surface_get(self->sys, idx, NULL, NULL, &pos_node, &neg_node, NULL);
    return Py_BuildValue("(ikk)", idx, (unsigned long)pos_node, (unsigned long)neg_node);
}

static PyObject* PyAleaSystem_box_general_surface(PyAleaSystemObject* self, PyObject* args) {
    int surface_id;
    double corner_x, corner_y, corner_z;
    double v1_x, v1_y, v1_z, v2_x, v2_y, v2_z, v3_x, v3_y, v3_z;
    if (!PyArg_ParseTuple(args, "idddddddddddd", &surface_id,
                          &corner_x, &corner_y, &corner_z,
                          &v1_x, &v1_y, &v1_z,
                          &v2_x, &v2_y, &v2_z,
                          &v3_x, &v3_y, &v3_z)) return NULL;
    if (!self->sys) { PyErr_SetString(PyExc_RuntimeError, "System not initialized"); return NULL; }

    int idx = alea_box_general_surface(self->sys, surface_id, corner_x, corner_y, corner_z,
                                           v1_x, v1_y, v1_z, v2_x, v2_y, v2_z, v3_x, v3_y, v3_z);
    if (idx < 0) {
        PyErr_SetString(PyExc_RuntimeError, "Failed to create general box surface");
        return NULL;
    }

    alea_node_id_t pos_node, neg_node;
    alea_surface_get(self->sys, idx, NULL, NULL, &pos_node, &neg_node, NULL);
    return Py_BuildValue("(ikk)", idx, (unsigned long)pos_node, (unsigned long)neg_node);
}

static PyObject* PyAleaSystem_sph_surface(PyAleaSystemObject* self, PyObject* args) {
    int surface_id;
    double cx, cy, cz, r;
    if (!PyArg_ParseTuple(args, "idddd", &surface_id, &cx, &cy, &cz, &r)) return NULL;
    if (!self->sys) { PyErr_SetString(PyExc_RuntimeError, "System not initialized"); return NULL; }

    int idx = alea_sph_surface(self->sys, surface_id, cx, cy, cz, r);
    if (idx < 0) {
        PyErr_SetString(PyExc_RuntimeError, "Failed to create SPH surface");
        return NULL;
    }

    alea_node_id_t pos_node, neg_node;
    alea_surface_get(self->sys, idx, NULL, NULL, &pos_node, &neg_node, NULL);
    return Py_BuildValue("(ikk)", idx, (unsigned long)pos_node, (unsigned long)neg_node);
}

static PyObject* PyAleaSystem_trc_surface(PyAleaSystemObject* self, PyObject* args) {
    int surface_id;
    double base_x, base_y, base_z, height_x, height_y, height_z, base_radius, top_radius;
    if (!PyArg_ParseTuple(args, "idddddddd", &surface_id, &base_x, &base_y, &base_z,
                          &height_x, &height_y, &height_z, &base_radius, &top_radius)) return NULL;
    if (!self->sys) { PyErr_SetString(PyExc_RuntimeError, "System not initialized"); return NULL; }

    int idx = alea_trc_surface(self->sys, surface_id, base_x, base_y, base_z,
                                   height_x, height_y, height_z, base_radius, top_radius);
    if (idx < 0) {
        PyErr_SetString(PyExc_RuntimeError, "Failed to create TRC surface");
        return NULL;
    }

    alea_node_id_t pos_node, neg_node;
    alea_surface_get(self->sys, idx, NULL, NULL, &pos_node, &neg_node, NULL);
    return Py_BuildValue("(ikk)", idx, (unsigned long)pos_node, (unsigned long)neg_node);
}

static PyObject* PyAleaSystem_ell_surface(PyAleaSystemObject* self, PyObject* args) {
    int surface_id;
    double v1_x, v1_y, v1_z, v2_x, v2_y, v2_z, major_axis_len;
    if (!PyArg_ParseTuple(args, "iddddddd", &surface_id, &v1_x, &v1_y, &v1_z,
                          &v2_x, &v2_y, &v2_z, &major_axis_len)) return NULL;
    if (!self->sys) { PyErr_SetString(PyExc_RuntimeError, "System not initialized"); return NULL; }

    int idx = alea_ell_surface(self->sys, surface_id, v1_x, v1_y, v1_z,
                                   v2_x, v2_y, v2_z, major_axis_len);
    if (idx < 0) {
        PyErr_SetString(PyExc_RuntimeError, "Failed to create ELL surface");
        return NULL;
    }

    alea_node_id_t pos_node, neg_node;
    alea_surface_get(self->sys, idx, NULL, NULL, &pos_node, &neg_node, NULL);
    return Py_BuildValue("(ikk)", idx, (unsigned long)pos_node, (unsigned long)neg_node);
}

static PyObject* PyAleaSystem_rec_surface(PyAleaSystemObject* self, PyObject* args) {
    int surface_id;
    double base_x, base_y, base_z, height_x, height_y, height_z;
    double axis1_x, axis1_y, axis1_z, axis2_x, axis2_y, axis2_z;
    if (!PyArg_ParseTuple(args, "idddddddddddd", &surface_id,
                          &base_x, &base_y, &base_z,
                          &height_x, &height_y, &height_z,
                          &axis1_x, &axis1_y, &axis1_z,
                          &axis2_x, &axis2_y, &axis2_z)) return NULL;
    if (!self->sys) { PyErr_SetString(PyExc_RuntimeError, "System not initialized"); return NULL; }

    int idx = alea_rec_surface(self->sys, surface_id, base_x, base_y, base_z,
                                   height_x, height_y, height_z,
                                   axis1_x, axis1_y, axis1_z,
                                   axis2_x, axis2_y, axis2_z);
    if (idx < 0) {
        PyErr_SetString(PyExc_RuntimeError, "Failed to create REC surface");
        return NULL;
    }

    alea_node_id_t pos_node, neg_node;
    alea_surface_get(self->sys, idx, NULL, NULL, &pos_node, &neg_node, NULL);
    return Py_BuildValue("(ikk)", idx, (unsigned long)pos_node, (unsigned long)neg_node);
}

static PyObject* PyAleaSystem_wed_surface(PyAleaSystemObject* self, PyObject* args) {
    int surface_id;
    double vertex_x, vertex_y, vertex_z;
    double v1_x, v1_y, v1_z, v2_x, v2_y, v2_z, v3_x, v3_y, v3_z;
    if (!PyArg_ParseTuple(args, "idddddddddddd", &surface_id,
                          &vertex_x, &vertex_y, &vertex_z,
                          &v1_x, &v1_y, &v1_z,
                          &v2_x, &v2_y, &v2_z,
                          &v3_x, &v3_y, &v3_z)) return NULL;
    if (!self->sys) { PyErr_SetString(PyExc_RuntimeError, "System not initialized"); return NULL; }

    int idx = alea_wed_surface(self->sys, surface_id, vertex_x, vertex_y, vertex_z,
                                   v1_x, v1_y, v1_z, v2_x, v2_y, v2_z, v3_x, v3_y, v3_z);
    if (idx < 0) {
        PyErr_SetString(PyExc_RuntimeError, "Failed to create WED surface");
        return NULL;
    }

    alea_node_id_t pos_node, neg_node;
    alea_surface_get(self->sys, idx, NULL, NULL, &pos_node, &neg_node, NULL);
    return Py_BuildValue("(ikk)", idx, (unsigned long)pos_node, (unsigned long)neg_node);
}

static PyObject* PyAleaSystem_rhp_surface(PyAleaSystemObject* self, PyObject* args) {
    int surface_id;
    double base_x, base_y, base_z, height_x, height_y, height_z;
    double r1_x, r1_y, r1_z, r2_x, r2_y, r2_z, r3_x, r3_y, r3_z;
    if (!PyArg_ParseTuple(args, "idddddddddddddddd", &surface_id,
                          &base_x, &base_y, &base_z,
                          &height_x, &height_y, &height_z,
                          &r1_x, &r1_y, &r1_z,
                          &r2_x, &r2_y, &r2_z,
                          &r3_x, &r3_y, &r3_z)) return NULL;
    if (!self->sys) { PyErr_SetString(PyExc_RuntimeError, "System not initialized"); return NULL; }

    int idx = alea_rhp_surface(self->sys, surface_id, base_x, base_y, base_z,
                                   height_x, height_y, height_z,
                                   r1_x, r1_y, r1_z,
                                   r2_x, r2_y, r2_z,
                                   r3_x, r3_y, r3_z);
    if (idx < 0) {
        PyErr_SetString(PyExc_RuntimeError, "Failed to create RHP surface");
        return NULL;
    }

    alea_node_id_t pos_node, neg_node;
    alea_surface_get(self->sys, idx, NULL, NULL, &pos_node, &neg_node, NULL);
    return Py_BuildValue("(ikk)", idx, (unsigned long)pos_node, (unsigned long)neg_node);
}

static PyObject* PyAleaSystem_get_surface_nodes(PyAleaSystemObject* self, PyObject* args) {
    int surface_idx;
    if (!PyArg_ParseTuple(args, "i", &surface_idx)) return NULL;
    if (!self->sys) { PyErr_SetString(PyExc_RuntimeError, "System not initialized"); return NULL; }

    if (surface_idx < 0 || (size_t)surface_idx >= alea_surface_count(self->sys)) {
        PyErr_SetString(PyExc_IndexError, "Surface index out of range");
        return NULL;
    }

    alea_node_id_t pos_node, neg_node;
    alea_surface_get(self->sys, surface_idx, NULL, NULL, &pos_node, &neg_node, NULL);
    return Py_BuildValue("(kk)", (unsigned long)pos_node, (unsigned long)neg_node);
}

/* ============================================================================
 * PyAleaSystem Methods - Cell Registration
 * ============================================================================ */

static PyObject* PyAleaSystem_add_cell(PyAleaSystemObject* self, PyObject* args, PyObject* kwds) {
    int cell_id;
    unsigned long root_node;
    int material_index = -1;  /* ALEA_MATERIAL_VOID */
    double density = 0.0;
    int universe_id = 0;
    static char* kwlist[] = {"cell_id", "root_node", "material_index", "density", "universe_id", NULL};

    if (!PyArg_ParseTupleAndKeywords(args, kwds, "ik|idi", kwlist,
                                     &cell_id, &root_node, &material_index, &density, &universe_id)) {
        return NULL;
    }
    if (!self->sys) { PyErr_SetString(PyExc_RuntimeError, "System not initialized"); return NULL; }

    int idx = alea_add_cell(self->sys, cell_id, (alea_node_id_t)root_node, material_index, density, universe_id);
    if (idx < 0) {
        PyErr_SetString(PyExc_RuntimeError, "Failed to add cell");
        return NULL;
    }
    return PyLong_FromLong(idx);
}

/* ============================================================================
 * PyAleaSystem Methods - Material Registration
 * ============================================================================ */

static PyObject* PyAleaSystem_add_material(PyAleaSystemObject* self, PyObject* args) {
    int material_id;
    if (!PyArg_ParseTuple(args, "i", &material_id)) return NULL;
    if (!self->sys) { PyErr_SetString(PyExc_RuntimeError, "System not initialized"); return NULL; }

    int idx = alea_add_material(self->sys, material_id);
    if (idx < 0) {
        PyErr_SetString(PyExc_RuntimeError, "Failed to add material");
        return NULL;
    }
    return PyLong_FromLong(idx);
}

static PyObject* PyAleaSystem_find_material_by_id(PyAleaSystemObject* self, PyObject* args) {
    int material_id;
    if (!PyArg_ParseTuple(args, "i", &material_id)) return NULL;
    if (!self->sys) { PyErr_SetString(PyExc_RuntimeError, "System not initialized"); return NULL; }

    int idx = alea_find_material_by_id(self->sys, material_id);
    if (idx < 0) {
        Py_RETURN_NONE;
    }
    return PyLong_FromLong(idx);
}

static PyObject* PyAleaSystem_cell_set_mixture(PyAleaSystemObject* self, PyObject* args) {
    int cell_index, mixture_id;
    if (!PyArg_ParseTuple(args, "ii", &cell_index, &mixture_id)) return NULL;
    if (!self->sys) { PyErr_SetString(PyExc_RuntimeError, "System not initialized"); return NULL; }

    if (alea_cell_set_mixture(self->sys, cell_index, mixture_id) < 0) {
        PyErr_Format(PyExc_RuntimeError, "Failed to set mixture %d on cell index %d", mixture_id, cell_index);
        return NULL;
    }
    Py_RETURN_NONE;
}

static PyObject* PyAleaSystem_find_mixture_by_id(PyAleaSystemObject* self, PyObject* args) {
    int mixture_id;
    if (!PyArg_ParseTuple(args, "i", &mixture_id)) return NULL;
    if (!self->sys) { PyErr_SetString(PyExc_RuntimeError, "System not initialized"); return NULL; }

    int idx = alea_find_mixture_by_id(self->sys, mixture_id);
    if (idx < 0) Py_RETURN_NONE;
    return PyLong_FromLong(idx);
}

/* ============================================================================
 * PyAleaSystem Methods - Material Composition API
 * ============================================================================ */

static PyObject* PyAleaSystem_material_count(PyAleaSystemObject* self, PyObject* Py_UNUSED(args)) {
    if (!self->sys) { PyErr_SetString(PyExc_RuntimeError, "System not initialized"); return NULL; }
    return PyLong_FromSize_t(alea_material_count(self->sys));
}

static PyObject* PyAleaSystem_material_get_id(PyAleaSystemObject* self, PyObject* args) {
    int mat_index;
    if (!PyArg_ParseTuple(args, "i", &mat_index)) return NULL;
    if (!self->sys) { PyErr_SetString(PyExc_RuntimeError, "System not initialized"); return NULL; }

    int mid = alea_material_get_id(self->sys, mat_index);
    if (mid < 0) {
        PyErr_Format(PyExc_IndexError, "Material index %d out of range", mat_index);
        return NULL;
    }
    return PyLong_FromLong(mid);
}

static PyObject* PyAleaSystem_material_add_nuclide(PyAleaSystemObject* self, PyObject* args) {
    int mat_index, zaid;
    const char* library = NULL;
    double fraction;
    if (!PyArg_ParseTuple(args, "iizd", &mat_index, &zaid, &library, &fraction)) return NULL;
    if (!self->sys) { PyErr_SetString(PyExc_RuntimeError, "System not initialized"); return NULL; }

    if (alea_material_add_nuclide(self->sys, mat_index, zaid, library, fraction) < 0) {
        PyErr_Format(PyExc_RuntimeError, "Failed to add nuclide %d to material index %d", zaid, mat_index);
        return NULL;
    }
    Py_RETURN_NONE;
}

static PyObject* PyAleaSystem_material_add_element(PyAleaSystemObject* self, PyObject* args) {
    int mat_index, Z;
    const char* library = NULL;
    double fraction;
    if (!PyArg_ParseTuple(args, "iizd", &mat_index, &Z, &library, &fraction)) return NULL;
    if (!self->sys) { PyErr_SetString(PyExc_RuntimeError, "System not initialized"); return NULL; }

    if (alea_material_add_element(self->sys, mat_index, Z, library, fraction) < 0) {
        PyErr_Format(PyExc_RuntimeError, "Failed to add element Z=%d to material index %d", Z, mat_index);
        return NULL;
    }
    Py_RETURN_NONE;
}

static PyObject* PyAleaSystem_material_set_density(PyAleaSystemObject* self, PyObject* args) {
    int mat_index;
    double density;
    if (!PyArg_ParseTuple(args, "id", &mat_index, &density)) return NULL;
    if (!self->sys) { PyErr_SetString(PyExc_RuntimeError, "System not initialized"); return NULL; }

    if (alea_material_set_density(self->sys, mat_index, density) < 0) {
        PyErr_Format(PyExc_IndexError, "Material index %d out of range", mat_index);
        return NULL;
    }
    Py_RETURN_NONE;
}

static PyObject* PyAleaSystem_material_set_weight_fraction(PyAleaSystemObject* self, PyObject* args) {
    int mat_index, is_weight;
    if (!PyArg_ParseTuple(args, "ip", &mat_index, &is_weight)) return NULL;
    if (!self->sys) { PyErr_SetString(PyExc_RuntimeError, "System not initialized"); return NULL; }

    if (alea_material_set_weight_fraction(self->sys, mat_index, is_weight) < 0) {
        PyErr_Format(PyExc_IndexError, "Material index %d out of range", mat_index);
        return NULL;
    }
    Py_RETURN_NONE;
}

static PyObject* PyAleaSystem_material_expand_elements(PyAleaSystemObject* self, PyObject* args) {
    int mat_index;
    if (!PyArg_ParseTuple(args, "i", &mat_index)) return NULL;
    if (!self->sys) { PyErr_SetString(PyExc_RuntimeError, "System not initialized"); return NULL; }

    if (alea_material_expand_elements(self->sys, mat_index) < 0) {
        PyErr_Format(PyExc_RuntimeError, "Failed to expand elements for material index %d", mat_index);
        return NULL;
    }
    Py_RETURN_NONE;
}

static PyObject* PyAleaSystem_material_nuclide_count(PyAleaSystemObject* self, PyObject* args) {
    int mat_index;
    if (!PyArg_ParseTuple(args, "i", &mat_index)) return NULL;
    if (!self->sys) { PyErr_SetString(PyExc_RuntimeError, "System not initialized"); return NULL; }
    return PyLong_FromSize_t(alea_material_nuclide_count(self->sys, mat_index));
}

static PyObject* PyAleaSystem_material_get_nuclides(PyAleaSystemObject* self, PyObject* args) {
    int mat_index;
    if (!PyArg_ParseTuple(args, "i", &mat_index)) return NULL;
    if (!self->sys) { PyErr_SetString(PyExc_RuntimeError, "System not initialized"); return NULL; }

    size_t count = alea_material_nuclide_count(self->sys, mat_index);
    PyObject* list = PyList_New(count);
    if (!list) return NULL;

    for (size_t i = 0; i < count; i++) {
        int zaid;
        const char* library;
        double fraction;
        if (alea_material_nuclide_get(self->sys, mat_index, i, &zaid, &library, &fraction) < 0) {
            Py_DECREF(list);
            PyErr_Format(PyExc_RuntimeError, "Failed to get nuclide %zu from material %d", i, mat_index);
            return NULL;
        }
        PyObject* entry = Py_BuildValue("{s:i,s:s,s:d}",
            "zaid", zaid,
            "library", library ? library : "",
            "fraction", fraction);
        if (!entry) { Py_DECREF(list); return NULL; }
        PyList_SET_ITEM(list, i, entry);
    }
    return list;
}

static PyObject* PyAleaSystem_material_element_count(PyAleaSystemObject* self, PyObject* args) {
    int mat_index;
    if (!PyArg_ParseTuple(args, "i", &mat_index)) return NULL;
    if (!self->sys) { PyErr_SetString(PyExc_RuntimeError, "System not initialized"); return NULL; }
    return PyLong_FromSize_t(alea_material_element_count(self->sys, mat_index));
}

static PyObject* PyAleaSystem_material_get_elements(PyAleaSystemObject* self, PyObject* args) {
    int mat_index;
    if (!PyArg_ParseTuple(args, "i", &mat_index)) return NULL;
    if (!self->sys) { PyErr_SetString(PyExc_RuntimeError, "System not initialized"); return NULL; }

    size_t count = alea_material_element_count(self->sys, mat_index);
    PyObject* list = PyList_New(count);
    if (!list) return NULL;

    for (size_t i = 0; i < count; i++) {
        int Z;
        const char* library;
        double fraction;
        if (alea_material_element_get(self->sys, mat_index, i, &Z, &library, &fraction) < 0) {
            Py_DECREF(list);
            PyErr_Format(PyExc_RuntimeError, "Failed to get element %zu from material %d", i, mat_index);
            return NULL;
        }
        PyObject* entry = Py_BuildValue("{s:i,s:s,s:d}",
            "Z", Z,
            "library", library ? library : "",
            "fraction", fraction);
        if (!entry) { Py_DECREF(list); return NULL; }
        PyList_SET_ITEM(list, i, entry);
    }
    return list;
}

static PyObject* PyAleaSystem_material_get_density(PyAleaSystemObject* self, PyObject* args) {
    int mat_index;
    if (!PyArg_ParseTuple(args, "i", &mat_index)) return NULL;
    if (!self->sys) { PyErr_SetString(PyExc_RuntimeError, "System not initialized"); return NULL; }

    double density;
    bool has_density;
    if (alea_material_get_density(self->sys, mat_index, &density, &has_density) < 0) {
        PyErr_Format(PyExc_IndexError, "Material index %d out of range", mat_index);
        return NULL;
    }
    if (!has_density) Py_RETURN_NONE;
    return PyFloat_FromDouble(density);
}

static PyObject* PyAleaSystem_material_is_weight_fraction(PyAleaSystemObject* self, PyObject* args) {
    int mat_index;
    if (!PyArg_ParseTuple(args, "i", &mat_index)) return NULL;
    if (!self->sys) { PyErr_SetString(PyExc_RuntimeError, "System not initialized"); return NULL; }

    bool is_weight = alea_material_is_weight_fraction(self->sys, mat_index);
    return PyBool_FromLong(is_weight);
}

/* ============================================================================
 * PyAleaSystem Methods - Mixture Query API
 * ============================================================================ */

static PyObject* PyAleaSystem_mixture_count(PyAleaSystemObject* self, PyObject* Py_UNUSED(args)) {
    if (!self->sys) { PyErr_SetString(PyExc_RuntimeError, "System not initialized"); return NULL; }
    return PyLong_FromSize_t(alea_mixture_count(self->sys));
}

static PyObject* PyAleaSystem_mixture_get_id(PyAleaSystemObject* self, PyObject* args) {
    int mix_index;
    if (!PyArg_ParseTuple(args, "i", &mix_index)) return NULL;
    if (!self->sys) { PyErr_SetString(PyExc_RuntimeError, "System not initialized"); return NULL; }

    int mid = alea_mixture_get_id(self->sys, mix_index);
    if (mid < 0) {
        PyErr_Format(PyExc_IndexError, "Mixture index %d out of range", mix_index);
        return NULL;
    }
    return PyLong_FromLong(mid);
}

static PyObject* PyAleaSystem_mixture_get_components(PyAleaSystemObject* self, PyObject* args) {
    int mix_index;
    if (!PyArg_ParseTuple(args, "i", &mix_index)) return NULL;
    if (!self->sys) { PyErr_SetString(PyExc_RuntimeError, "System not initialized"); return NULL; }

    size_t count = alea_mixture_component_count(self->sys, mix_index);
    PyObject* list = PyList_New(count);
    if (!list) return NULL;

    for (size_t i = 0; i < count; i++) {
        int material_id;
        double fraction;
        if (alea_mixture_component_get(self->sys, mix_index, i, &material_id, &fraction) < 0) {
            Py_DECREF(list);
            PyErr_Format(PyExc_RuntimeError, "Failed to get component %zu from mixture %d", i, mix_index);
            return NULL;
        }
        PyObject* entry = Py_BuildValue("{s:i,s:d}",
            "material_id", material_id,
            "fraction", fraction);
        if (!entry) { Py_DECREF(list); return NULL; }
        PyList_SET_ITEM(list, i, entry);
    }
    return list;
}
