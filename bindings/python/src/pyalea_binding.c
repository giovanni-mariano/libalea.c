// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/**
 * @file pyalea_binding.c
 * @brief Python C extension for libalea
 *
 * Uses ONLY the public API from alea.h (alea_* functions).
 *
 * Build: make -C bindings/python
 *
 * Implementation is split across multiple files that are #included below
 * to keep a single compilation unit.
 */

#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <limits.h>
#include <math.h>
#include <structmember.h>
#include <stdint.h>
#include <string.h>
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/arrayobject.h>
/* numpy/complex.h defines I as the imaginary-unit macro; libalea's public
 * quadric API has a parameter named I. */
#ifdef I
#undef I
#endif
/* Public API headers */
#include "alea.h"
#include "alea_mcnp.h"
#include "alea_openmc.h"
#include "alea_serpent.h"
#include "alea_slice.h"
#include "alea_raycast.h"
#include "alea_render.h"
#include "alea_geo_validator.h"
#include "alea_mesh.h"
#include "alea_nucdata.h"
#include "alea_transport.h"
#include "alea_source.h"

#include <signal.h>

/* SIGINT cooperative-interruption helpers.
 * While the GIL is released Python cannot run its own SIGINT handler.
 * We temporarily install a handler that sets the C library's interrupt flag
 * so that long-running C operations return early with ALEA_ERR_INTERRUPTED. */
static void sigint_handler(int sig) { (void)sig; alea_interrupt(); }

typedef void (*sighandler_func)(int);

/* These helpers are entered with the GIL held, so the scope count and saved
 * handler are serialized even while earlier native calls are running without
 * the GIL. SIGINT is intentionally a broadcast to every active cooperative
 * operation; only the final scope restores the Python handler and clears it. */
static size_t active_sigint_scopes;
static sighandler_func previous_sigint_handler;

static inline sighandler_func install_sigint(void) {
    if (active_sigint_scopes++ == 0)
        previous_sigint_handler = signal(SIGINT, sigint_handler);
    return previous_sigint_handler;
}

/* Restore the previous handler, check whether the interrupt flag was set,
 * and if so raise KeyboardInterrupt.  Returns 1 if interrupted, 0 otherwise. */
static inline int restore_sigint(sighandler_func old) {
    (void)old;
    int was_interrupted = alea_interrupted();
    if (active_sigint_scopes > 0 && --active_sigint_scopes == 0) {
        signal(SIGINT, previous_sigint_handler);
        alea_clear_interrupt();
    }
    if (was_interrupted) {
        PyErr_SetNone(PyExc_KeyboardInterrupt);
    }
    return was_interrupted;
}

/* ============================================================================
 * Forward Declarations
 * ============================================================================ */

static PyTypeObject PyAleaSystemType;
static PyTypeObject PyAleaVoidResultType;

/* ============================================================================
 * PyAleaSystem Python Type
 * ============================================================================ */

typedef struct {
    PyObject_HEAD
    alea_system_t* sys;
    int owns_sys;
    mcnp_model_t* mcnp_model;  /* Optional MCNP-only metadata sidecar. */
} PyAleaSystemObject;

/* ============================================================================
 * Shared Helpers
 * ============================================================================ */

/* Set a dict item from a new reference and always release the caller's
 * reference. PyDict_SetItemString() retains its own reference on success. */
static int dict_set_new(PyObject* dict, const char* key, PyObject* value) {
    if (!value) return -1;
    if (!dict) {
        Py_DECREF(value);
        return -1;
    }
    int rc = PyDict_SetItemString(dict, key, value);
    Py_DECREF(value);
    return rc;
}

/* Helper: add lattice fields to an existing cell dict.
 * Only adds keys when lat_type != 0 to avoid clutter for non-lattice cells.
 * Returns 0 on success, -1 on failure (with Python exception set). */
static int add_lattice_fields(PyObject* dict, const alea_cell_info_t* info) {
    if (info->lat_type == 0)
        return 0;

    if (dict_set_new(dict, "lat_type",
            PyLong_FromLong(info->lat_type)) < 0)
        return -1;

    PyObject* dims = Py_BuildValue("(iiiiii)",
        info->lat_fill_dims[0], info->lat_fill_dims[1],
        info->lat_fill_dims[2], info->lat_fill_dims[3],
        info->lat_fill_dims[4], info->lat_fill_dims[5]);
    if (!dims) return -1;
    if (PyDict_SetItemString(dict, "lat_fill_dims", dims) < 0) {
        Py_DECREF(dims); return -1;
    }
    Py_DECREF(dims);

    PyObject* pitch = Py_BuildValue("(ddd)",
        info->lat_pitch[0], info->lat_pitch[1], info->lat_pitch[2]);
    if (!pitch) return -1;
    if (PyDict_SetItemString(dict, "lat_pitch", pitch) < 0) {
        Py_DECREF(pitch); return -1;
    }
    Py_DECREF(pitch);

    PyObject* ll = Py_BuildValue("(ddd)",
        info->lat_lower_left[0], info->lat_lower_left[1], info->lat_lower_left[2]);
    if (!ll) return -1;
    if (PyDict_SetItemString(dict, "lat_lower_left", ll) < 0) {
        Py_DECREF(ll); return -1;
    }
    Py_DECREF(ll);

    PyObject* fill_count = PyLong_FromSize_t(info->lat_fill_count);
    if (!fill_count) return -1;
    if (PyDict_SetItemString(dict, "lat_fill_count", fill_count) < 0) {
        Py_DECREF(fill_count); return -1;
    }
    Py_DECREF(fill_count);

    PyObject* outer_universe = info->lat_outer_universe >= 0
        ? PyLong_FromLong(info->lat_outer_universe)
        : Py_NewRef(Py_None);
    if (!outer_universe) return -1;
    if (PyDict_SetItemString(dict, "lat_outer_universe", outer_universe) < 0) {
        Py_DECREF(outer_universe); return -1;
    }
    Py_DECREF(outer_universe);

    if (info->lat_fill && info->lat_fill_count > 0) {
        PyObject* fill = PyList_New((Py_ssize_t)info->lat_fill_count);
        if (!fill) return -1;
        for (size_t i = 0; i < info->lat_fill_count; i++) {
            PyList_SET_ITEM(fill, (Py_ssize_t)i, PyLong_FromLong(info->lat_fill[i]));
        }
        if (PyDict_SetItemString(dict, "lat_fill", fill) < 0) {
            Py_DECREF(fill); return -1;
        }
        Py_DECREF(fill);
    } else {
        if (PyDict_SetItemString(dict, "lat_fill", Py_None) < 0)
            return -1;
    }

    return 0;
}

/* Helper: add comment fields to an existing cell dict.
 * Only adds keys when a comment is present to avoid clutter.
 * Returns 0 on success, -1 on failure (with Python exception set). */
static int add_comment_fields(PyObject* dict, const alea_cell_info_t* info) {
    if (info->comments) {
        PyObject* val = PyUnicode_FromString(info->comments);
        if (!val) return -1;
        if (PyDict_SetItemString(dict, "comments", val) < 0) {
            Py_DECREF(val); return -1;
        }
        Py_DECREF(val);
    }
    if (info->inline_comment) {
        PyObject* val = PyUnicode_FromString(info->inline_comment);
        if (!val) return -1;
        if (PyDict_SetItemString(dict, "inline_comment", val) < 0) {
            Py_DECREF(val); return -1;
        }
        Py_DECREF(val);
    }
    return 0;
}

/* Add generic cell temperature metadata. Temperature is stored in Kelvin;
 * None distinguishes an unset value from a physical temperature. */
static int add_temperature_fields(PyObject* dict,
                                  const alea_cell_info_t* info) {
    PyObject* value = info->has_temperature
        ? PyFloat_FromDouble(info->temperature)
        : Py_NewRef(Py_None);
    if (!value) return -1;
    int rc = PyDict_SetItemString(dict, "temperature", value);
    Py_DECREF(value);
    if (rc < 0) return -1;
    return PyDict_SetItemString(
        dict, "has_temperature",
        info->has_temperature ? Py_True : Py_False);
}

static mcnp_model_t* ensure_mcnp_sidecar(PyAleaSystemObject* self) {
    if (!self->mcnp_model && self->sys)
        self->mcnp_model = mcnp_model_wrap(self->sys);
    return self->mcnp_model;
}

static int add_importance_fields(PyAleaSystemObject* self, PyObject* dict,
                                 size_t cell_index) {
    if (!self->mcnp_model)
        return 0;
    const mcnp_cell_params_t* params = mcnp_cell_params_const(
        self->mcnp_model, cell_index);
    if (!params)
        return 0;
    PyObject* importance = Py_BuildValue("{s:d,s:d,s:d}",
                                         "neutron", params->imp_n,
                                         "photon", params->imp_p,
                                         "electron", params->imp_e);
    if (!importance)
        return -1;
    int rc = PyDict_SetItemString(dict, "importance", importance);
    Py_DECREF(importance);
    return rc;
}

/* Copy only MCNP particle importance into a derived geometry object.  The
 * generic alea_system_t deliberately remains unaware of MCNP-only metadata.
 * Cell IDs are used instead of array indices because extraction may compact
 * or reorder cells. */
static int copy_mcnp_importance_sidecar(PyAleaSystemObject* source,
                                        PyAleaSystemObject* target) {
    if (!source->mcnp_model || !target->sys)
        return 0;

    mcnp_model_t* copied = mcnp_model_wrap(target->sys);
    if (!copied) {
        PyErr_NoMemory();
        return -1;
    }

    const size_t target_count = alea_cell_count(target->sys);
    for (size_t target_index = 0; target_index < target_count; target_index++) {
        alea_cell_info_t info;
        if (alea_cell_get_info(target->sys, target_index, &info) < 0)
            continue;
        const int source_index = alea_cell_find(source->sys, info.cell_id);
        if (source_index < 0)
            continue;
        const mcnp_cell_params_t* from = mcnp_cell_params_const(
            source->mcnp_model, (size_t)source_index);
        mcnp_cell_params_t* to = mcnp_cell_params(copied, target_index);
        if (!from || !to)
            continue;
        to->imp_n = from->imp_n;
        to->imp_p = from->imp_p;
        to->imp_e = from->imp_e;
        to->has_imp_n = from->has_imp_n;
        to->has_imp_p = from->has_imp_p;
        to->has_imp_e = from->has_imp_e;
    }

    target->mcnp_model = copied;
    return 0;
}

static const char* curve_type_to_string(alea_curve_type_t type) {
    switch (type) {
        case ALEA_CURVE_LINE: return "line";
        case ALEA_CURVE_LINE_SEGMENT: return "line_segment";
        case ALEA_CURVE_CIRCLE: return "circle";
        case ALEA_CURVE_ARC: return "arc";
        case ALEA_CURVE_ELLIPSE: return "ellipse";
        case ALEA_CURVE_ELLIPSE_ARC: return "ellipse_arc";
        case ALEA_CURVE_PARABOLA: return "parabola";
        case ALEA_CURVE_HYPERBOLA: return "hyperbola";
        case ALEA_CURVE_POLYGON: return "polygon";
        case ALEA_CURVE_QUARTIC: return "quartic";
        case ALEA_CURVE_PARALLEL_LINES: return "parallel_lines";
        default: return "none";
    }
}

static PyObject* build_node_tree(const alea_system_t* sys, alea_node_id_t node) {
    alea_operation_t op = alea_node_operation(sys, node);

    if (op == ALEA_OP_PRIMITIVE) {
        int sid = alea_node_surface_id(sys, node);
        int sense = alea_node_sense(sys, node);
        return Py_BuildValue("(iii)", (int)op, sid, sense);
    }

    if (op == ALEA_OP_COMPLEMENT) {
        alea_node_id_t child = alea_node_left(sys, node);
        PyObject* child_tree = build_node_tree(sys, child);
        if (!child_tree) return NULL;
        PyObject* result = PyTuple_New(2);
        if (!result) { Py_DECREF(child_tree); return NULL; }
        PyTuple_SET_ITEM(result, 0, PyLong_FromLong((int)op));
        PyTuple_SET_ITEM(result, 1, child_tree);
        return result;
    }

    /* Binary: UNION, INTERSECTION, DIFFERENCE */
    alea_node_id_t left = alea_node_left(sys, node);
    alea_node_id_t right = alea_node_right(sys, node);
    PyObject* left_tree = build_node_tree(sys, left);
    if (!left_tree) return NULL;
    PyObject* right_tree = build_node_tree(sys, right);
    if (!right_tree) { Py_DECREF(left_tree); return NULL; }
    PyObject* result = PyTuple_New(3);
    if (!result) { Py_DECREF(left_tree); Py_DECREF(right_tree); return NULL; }
    PyTuple_SET_ITEM(result, 0, PyLong_FromLong((int)op));
    PyTuple_SET_ITEM(result, 1, left_tree);
    PyTuple_SET_ITEM(result, 2, right_tree);
    return result;
}

static int ensure_query_acceleration(PyAleaSystemObject* self) {
    if (!self->sys) {
        PyErr_SetString(PyExc_RuntimeError, "System not initialized");
        return -1;
    }
    int rc;
    sighandler_func old_sigint = install_sigint();
    Py_BEGIN_ALLOW_THREADS
    rc = alea_prepare_query_acceleration(self->sys);
    Py_END_ALLOW_THREADS
    if (restore_sigint(old_sigint)) return -1;
    if (rc != 0) {
        PyErr_SetString(PyExc_RuntimeError, alea_error());
        return -1;
    }
    return 0;
}

static alea_slice_curves_t* get_slice_curves_allow_threads(
        alea_system_t* sys, const alea_slice_view_t* view) {
    alea_slice_curves_t* curves = NULL;
    sighandler_func old_sigint = install_sigint();
    Py_BEGIN_ALLOW_THREADS
    curves = alea_get_slice_curves(sys, view);
    Py_END_ALLOW_THREADS
    if (restore_sigint(old_sigint)) {
        alea_slice_curves_free(curves);
        return NULL;
    }
    return curves;
}

/* ============================================================================
 * Implementation Files
 * ============================================================================ */

#include "_bind_core.c"
#include "_bind_geometry.c"
#include "_bind_io.c"
#include "_bind_util.c"
#include "_bind_raycast.c"
#include "_bind_render.c"
#include "_bind_slice.c"
#include "_bind_inspect.c"
#include "_bind_mesh.c"
#include "_bind_geo_validator.c"
#include "_bind_cluster.c"

/* ============================================================================
 * PyAleaSystem Getters/Setters Table
 * ============================================================================ */

static PyGetSetDef PyAleaSystem_getsetters[] = {
    {"cell_count", (getter)PyAleaSystem_get_cell_count, NULL, "Number of cells", NULL},
    {"surface_count", (getter)PyAleaSystem_get_surface_count, NULL, "Number of surfaces", NULL},
    {"universe_count", (getter)PyAleaSystem_get_universe_count, NULL, "Number of universes", NULL},
    {"volume_path_count", (getter)PyAleaSystem_get_volume_path_count, NULL, "Number of concrete hierarchical volume paths", NULL},
    {NULL}
};

/* ============================================================================
 * PyAleaSystem Method Table
 * ============================================================================ */

static PyMethodDef PyAleaSystem_methods[] = {
    /* Queries */
    {"find_cell", (PyCFunction)PyAleaSystem_find_cell, METH_VARARGS,
     "find_cell(x, y, z) -> (cell_index, material_id, resolution_flags) or None\n\n"
     "Find cell containing point. In undefined fill regions (a fill/lattice\n"
     "container whose filling universe has no cell at the point) the container\n"
     "is returned with resolution_flags bit 0x1 set."},
    {"point_inside", (PyCFunction)PyAleaSystem_point_inside, METH_VARARGS,
     "point_inside(node_id, x, y, z) -> bool\n\nTest if point is inside CSG tree."},
    {"material_at", (PyCFunction)PyAleaSystem_material_at, METH_VARARGS,
     "material_at(x, y, z) -> int or None\n\nGet material ID at point."},
    {"find_overlaps", (PyCFunction)PyAleaSystem_find_overlaps, METH_VARARGS,
     "find_overlaps(max_pairs=100) -> list of (cell_idx, cell_idx)\n\nFind overlapping cells."},
    {"request_interrupt", (PyCFunction)PyAleaSystem_request_interrupt, METH_NOARGS,
     "request_interrupt()\n\nBroadcast cancellation to all active cooperative native operations. "
     "The interrupt flag is cleared after the final active operation returns."},

    /* Cell operations */
    {"get_cell", (PyCFunction)PyAleaSystem_get_cell, METH_VARARGS,
     "get_cell(cell_id) -> dict\n\nGet cell by MCNP ID."},
    {"get_cell_by_index", (PyCFunction)PyAleaSystem_get_cell_by_index, METH_VARARGS,
     "get_cell_by_index(index) -> dict\n\nGet cell by array index."},
    {"get_cells", (PyCFunction)PyAleaSystem_get_cells, METH_NOARGS,
     "get_cells() -> list of dict\n\nGet all cells."},
    {"_color_identifiers", (PyCFunction)PyAleaSystem_color_identifiers, METH_NOARGS,
     "Return unique categorical identifiers for deterministic color initialization."},
    {"cell_set_importance", (PyCFunction)PyAleaSystem_cell_set_importance,
     METH_VARARGS, "cell_set_importance(index, particle, value)"},
    {"cell_set_temperature", (PyCFunction)PyAleaSystem_cell_set_temperature,
     METH_VARARGS, "cell_set_temperature(index, temperature_K)"},
    {"cell_clear_temperature", (PyCFunction)PyAleaSystem_cell_clear_temperature,
     METH_VARARGS, "cell_clear_temperature(index)"},
    {"cell_find", (PyCFunction)PyAleaSystem_cell_find, METH_VARARGS,
     "cell_find(cell_id) -> int or None\n\nFind cell index by MCNP cell ID. O(1) lookup."},
    {"find_all_cells", (PyCFunction)PyAleaSystem_find_all_cells, METH_VARARGS,
     "find_all_cells(x, y, z) -> list of dict\n\n"
     "Find all cells at a point across hierarchy depths.\n"
     "Returns list of dicts with cell_id, cell_index, material_id, universe_id, fill_universe, depth, local_x/y/z."},
    {"find_all_cells_coverage", (PyCFunction)PyAleaSystem_find_all_cells_coverage,
     METH_VARARGS | METH_KEYWORDS,
     "find_all_cells_coverage(x, y, z, max_hits=256, universe_id=0, universe_depth=-1) -> dict\n\n"
     "Complete occurrence-aware diagnostic ownership query and native classification. "
     "never interpret hits as a complete owner set when truncated is true."},

    /* CSG Node Inspection */
    {"node_tree", (PyCFunction)PyAleaSystem_node_tree, METH_VARARGS,
     "node_tree(node_id) -> tuple\n\n"
     "Walk CSG tree and return nested tuple structure.\n"
     "Primitive: (0, surface_id, sense)\n"
     "Union: (1, left, right)\n"
     "Intersection: (2, left, right)\n"
     "Difference: (3, left, right)\n"
     "Complement: (4, child)\n"},

    /* Universe operations */
    {"build_universe_index", (PyCFunction)PyAleaSystem_build_universe_index, METH_NOARGS,
     "build_universe_index()\n\nBuild universe lookup tables."},
    {"flatten_universe", (PyCFunction)PyAleaSystem_flatten_universe, METH_VARARGS,
     "flatten_universe(universe_id=0) -> int\n\nFlatten universe hierarchy."},
    {"build_spatial_index", (PyCFunction)PyAleaSystem_build_spatial_index, METH_NOARGS,
     "build_spatial_index()\n\nBuild the hierarchical spatial index / query caches for fast queries."},
    {"prepare_query_acceleration", (PyCFunction)PyAleaSystem_prepare_query_acceleration, METH_NOARGS,
     "prepare_query_acceleration()\n\nBuild all query caches used by point queries, grids, and raycasts."},
    {"get_universe", (PyCFunction)PyAleaSystem_get_universe, METH_VARARGS,
     "get_universe(universe_id) -> dict\n\nGet universe info."},

    /* Primitive creation */
    {"create_plane", (PyCFunction)PyAleaSystem_create_plane, METH_VARARGS,
     "create_plane(a, b, c, d, sense) -> node_id\n\nCreate plane halfspace."},
    {"create_sphere", (PyCFunction)PyAleaSystem_create_sphere, METH_VARARGS,
     "create_sphere(cx, cy, cz, radius, sense) -> node_id\n\nCreate sphere halfspace."},
    {"create_box", (PyCFunction)PyAleaSystem_create_box, METH_VARARGS,
     "create_box(xmin, xmax, ymin, ymax, zmin, zmax, sense) -> node_id\n\nCreate box halfspace."},
    {"create_cylinder_z", (PyCFunction)PyAleaSystem_create_cylinder_z, METH_VARARGS,
     "create_cylinder_z(cx, cy, radius, sense) -> node_id\n\nCreate Z-cylinder halfspace."},

    /* Boolean operations */
    {"create_union", (PyCFunction)PyAleaSystem_create_union, METH_VARARGS,
     "create_union(a, b) -> node_id\n\nCreate union of two nodes."},
    {"create_intersection", (PyCFunction)PyAleaSystem_create_intersection, METH_VARARGS,
     "create_intersection(a, b) -> node_id\n\nCreate intersection of two nodes."},
    {"create_difference", (PyCFunction)PyAleaSystem_create_difference, METH_VARARGS,
     "create_difference(a, b) -> node_id\n\nCreate difference (a - b)."},
    {"create_complement", (PyCFunction)PyAleaSystem_create_complement, METH_VARARGS,
     "create_complement(a) -> node_id\n\nCreate complement (not a)."},
    {"create_union_many", (PyCFunction)PyAleaSystem_create_union_many, METH_VARARGS,
     "create_union_many([nodes]) -> node_id\n\nCreate union of multiple nodes."},
    {"create_intersection_many", (PyCFunction)PyAleaSystem_create_intersection_many, METH_VARARGS,
     "create_intersection_many([nodes]) -> node_id\n\nCreate intersection of multiple nodes."},

    /* Surface creation (with automatic registration for raycast) */
    {"sphere_surface", (PyCFunction)PyAleaSystem_sphere_surface, METH_VARARGS,
     "sphere_surface(surface_id, cx, cy, cz, radius) -> (index, pos_node, neg_node)\n\n"
     "Create sphere surface with both halfspace nodes registered for raycast."},
    {"cylinder_z_surface", (PyCFunction)PyAleaSystem_cylinder_z_surface, METH_VARARGS,
     "cylinder_z_surface(surface_id, cx, cy, radius) -> (index, pos_node, neg_node)\n\n"
     "Create Z-cylinder surface with both halfspace nodes registered for raycast."},
    {"box_surface", (PyCFunction)PyAleaSystem_box_surface, METH_VARARGS,
     "box_surface(surface_id, xmin, xmax, ymin, ymax, zmin, zmax) -> (index, pos_node, neg_node)\n\n"
     "Create box surface with both halfspace nodes registered for raycast."},
    {"plane_surface", (PyCFunction)PyAleaSystem_plane_surface, METH_VARARGS,
     "plane_surface(surface_id, a, b, c, d) -> (index, pos_node, neg_node)\n\n"
     "Create plane surface with both halfspace nodes registered for raycast."},
    {"cylinder_x_surface", (PyCFunction)PyAleaSystem_cylinder_x_surface, METH_VARARGS,
     "cylinder_x_surface(surface_id, cy, cz, radius) -> (index, pos_node, neg_node)\n\n"
     "Create X-cylinder surface with both halfspace nodes registered for raycast."},
    {"cylinder_y_surface", (PyCFunction)PyAleaSystem_cylinder_y_surface, METH_VARARGS,
     "cylinder_y_surface(surface_id, cx, cz, radius) -> (index, pos_node, neg_node)\n\n"
     "Create Y-cylinder surface with both halfspace nodes registered for raycast."},
    {"cone_z_surface", (PyCFunction)PyAleaSystem_cone_z_surface, METH_VARARGS,
     "cone_z_surface(surface_id, cx, cy, cz, t_squared, sheet=0) -> (index, pos_node, neg_node)\n\n"
     "Create Z-cone surface with both halfspace nodes registered for raycast."},
    {"cone_x_surface", (PyCFunction)PyAleaSystem_cone_x_surface, METH_VARARGS,
     "cone_x_surface(surface_id, cx, cy, cz, t_squared, sheet=0) -> (index, pos_node, neg_node)\n\n"
     "Create X-cone surface with both halfspace nodes registered for raycast."},
    {"cone_y_surface", (PyCFunction)PyAleaSystem_cone_y_surface, METH_VARARGS,
     "cone_y_surface(surface_id, cx, cy, cz, t_squared, sheet=0) -> (index, pos_node, neg_node)\n\n"
     "Create Y-cone surface with both halfspace nodes registered for raycast."},
    {"torus_z_surface", (PyCFunction)PyAleaSystem_torus_z_surface, METH_VARARGS,
     "torus_z_surface(surface_id, cx, cy, cz, major_radius, minor_radius) -> (index, pos_node, neg_node)\n\n"
     "Create Z-torus surface with both halfspace nodes registered for raycast."},
    {"torus_x_surface", (PyCFunction)PyAleaSystem_torus_x_surface, METH_VARARGS,
     "torus_x_surface(surface_id, cx, cy, cz, major_radius, minor_radius) -> (index, pos_node, neg_node)\n\n"
     "Create X-torus surface with both halfspace nodes registered for raycast."},
    {"torus_y_surface", (PyCFunction)PyAleaSystem_torus_y_surface, METH_VARARGS,
     "torus_y_surface(surface_id, cx, cy, cz, major_radius, minor_radius) -> (index, pos_node, neg_node)\n\n"
     "Create Y-torus surface with both halfspace nodes registered for raycast."},
    {"quadric_surface", (PyCFunction)PyAleaSystem_quadric_surface, METH_VARARGS,
     "quadric_surface(surface_id, A, B, C, D, E, F, G, H, I, J) -> (index, pos_node, neg_node)\n\n"
     "Create general quadric surface (Ax² + By² + Cz² + Dxy + Eyz + Fzx + Gx + Hy + Iz + J = 0)."},
    {"rcc_surface", (PyCFunction)PyAleaSystem_rcc_surface, METH_VARARGS,
     "rcc_surface(surface_id, base_x, base_y, base_z, height_x, height_y, height_z, radius) -> (index, pos_node, neg_node)\n\n"
     "Create RCC (Right Circular Cylinder) macrobody surface."},
    {"box_general_surface", (PyCFunction)PyAleaSystem_box_general_surface, METH_VARARGS,
     "box_general_surface(surface_id, corner_x, corner_y, corner_z, v1_x, v1_y, v1_z, v2_x, v2_y, v2_z, v3_x, v3_y, v3_z) -> (index, pos_node, neg_node)\n\n"
     "Create general box macrobody from corner and three edge vectors."},
    {"sph_surface", (PyCFunction)PyAleaSystem_sph_surface, METH_VARARGS,
     "sph_surface(surface_id, cx, cy, cz, r) -> (index, pos_node, neg_node)\n\n"
     "Create SPH (sphere) macrobody surface."},
    {"trc_surface", (PyCFunction)PyAleaSystem_trc_surface, METH_VARARGS,
     "trc_surface(surface_id, base_x, base_y, base_z, height_x, height_y, height_z, base_radius, top_radius) -> (index, pos_node, neg_node)\n\n"
     "Create TRC (Truncated Right Cone) macrobody surface."},
    {"ell_surface", (PyCFunction)PyAleaSystem_ell_surface, METH_VARARGS,
     "ell_surface(surface_id, v1_x, v1_y, v1_z, v2_x, v2_y, v2_z, major_axis_len) -> (index, pos_node, neg_node)\n\n"
     "Create ELL (ellipsoid) macrobody from two foci and major axis length."},
    {"rec_surface", (PyCFunction)PyAleaSystem_rec_surface, METH_VARARGS,
     "rec_surface(surface_id, base_x, base_y, base_z, height_x, height_y, height_z, axis1_x, axis1_y, axis1_z, axis2_x, axis2_y, axis2_z) -> (index, pos_node, neg_node)\n\n"
     "Create REC (Right Elliptical Cylinder) macrobody surface."},
    {"wed_surface", (PyCFunction)PyAleaSystem_wed_surface, METH_VARARGS,
     "wed_surface(surface_id, vertex_x, vertex_y, vertex_z, v1_x, v1_y, v1_z, v2_x, v2_y, v2_z, v3_x, v3_y, v3_z) -> (index, pos_node, neg_node)\n\n"
     "Create WED (wedge) macrobody from vertex and three edge vectors."},
    {"rhp_surface", (PyCFunction)PyAleaSystem_rhp_surface, METH_VARARGS,
     "rhp_surface(surface_id, base_x, base_y, base_z, height_x, height_y, height_z, r1_x, r1_y, r1_z, r2_x, r2_y, r2_z, r3_x, r3_y, r3_z) -> (index, pos_node, neg_node)\n\n"
     "Create RHP (Right Hexagonal Prism) macrobody surface."},
    {"surface_set_boundary", (PyCFunction)PyAleaSystem_surface_set_boundary,
     METH_VARARGS,
     "surface_set_boundary(surface_id, boundary) -> None\n\n"
     "Set transmissive, reflective, or vacuum boundary metadata."},
    {"surface_get_boundary", (PyCFunction)PyAleaSystem_surface_get_boundary,
     METH_VARARGS,
     "surface_get_boundary(surface_id) -> str\n\nReturn boundary metadata."},
    {"get_surface_nodes", (PyCFunction)PyAleaSystem_get_surface_nodes, METH_VARARGS,
     "get_surface_nodes(surface_index) -> (pos_node, neg_node)\n\n"
     "Get halfspace nodes for a registered surface."},

    /* Cell registration */
    {"add_cell", (PyCFunction)PyAleaSystem_add_cell,
     METH_VARARGS | METH_KEYWORDS,
     "add_cell(cell_id, root_node, material_index=-1, density=0.0, universe_id=0) -> index\n\n"
     "Register a cell. material_index is a material index from add_material(), or -1 for void."},

    /* Material registration */
    {"add_material", (PyCFunction)PyAleaSystem_add_material, METH_VARARGS,
     "add_material(material_id) -> index\n\nRegister a material by MCNP ID. Returns material index."},
    {"find_material_by_id", (PyCFunction)PyAleaSystem_find_material_by_id, METH_VARARGS,
     "find_material_by_id(material_id) -> index or None\n\nFind material index by MCNP ID."},

    /* Material composition */
    {"material_count", (PyCFunction)PyAleaSystem_material_count, METH_NOARGS,
     "material_count() -> int\n\nGet number of materials."},
    {"material_get_id", (PyCFunction)PyAleaSystem_material_get_id, METH_VARARGS,
     "material_get_id(mat_index) -> int\n\nGet material MCNP ID by index."},
    {"material_add_nuclide", (PyCFunction)PyAleaSystem_material_add_nuclide, METH_VARARGS,
     "material_add_nuclide(mat_index, zaid, library, fraction)\n\n"
     "Add a nuclide to a material. library can be None."},
    {"material_add_element", (PyCFunction)PyAleaSystem_material_add_element, METH_VARARGS,
     "material_add_element(mat_index, Z, library, fraction)\n\n"
     "Add an element (natural composition) to a material. library can be None."},
    {"material_set_density", (PyCFunction)PyAleaSystem_material_set_density, METH_VARARGS,
     "material_set_density(mat_index, density)\n\nSet material standard density."},
    {"material_set_weight_fraction", (PyCFunction)PyAleaSystem_material_set_weight_fraction, METH_VARARGS,
     "material_set_weight_fraction(mat_index, is_weight)\n\nSet whether fractions are weight (True) or atom (False)."},
    {"material_expand_elements", (PyCFunction)PyAleaSystem_material_expand_elements, METH_VARARGS,
     "material_expand_elements(mat_index)\n\nExpand element entries to explicit nuclides."},
    {"material_nuclide_count", (PyCFunction)PyAleaSystem_material_nuclide_count, METH_VARARGS,
     "material_nuclide_count(mat_index) -> int\n\nGet number of nuclides in a material."},
    {"material_get_nuclides", (PyCFunction)PyAleaSystem_material_get_nuclides, METH_VARARGS,
     "material_get_nuclides(mat_index) -> list of {zaid, library, fraction}\n\n"
     "Get all nuclides in a material."},
    {"material_element_count", (PyCFunction)PyAleaSystem_material_element_count, METH_VARARGS,
     "material_element_count(mat_index) -> int\n\nGet number of elements in a material."},
    {"material_get_elements", (PyCFunction)PyAleaSystem_material_get_elements, METH_VARARGS,
     "material_get_elements(mat_index) -> list of {Z, library, fraction}\n\n"
     "Get all elements in a material."},
    {"material_get_density", (PyCFunction)PyAleaSystem_material_get_density, METH_VARARGS,
     "material_get_density(mat_index) -> float or None\n\nGet material density, or None if not set."},
    {"material_is_weight_fraction", (PyCFunction)PyAleaSystem_material_is_weight_fraction, METH_VARARGS,
     "material_is_weight_fraction(mat_index) -> bool\n\nTrue if material uses weight fractions."},

    /* Cell-mixture assignment */
    {"cell_set_mixture", (PyCFunction)PyAleaSystem_cell_set_mixture, METH_VARARGS,
     "cell_set_mixture(cell_index, mixture_id)\n\nSet cell material to a mixture."},

    /* Mixture queries */
    {"find_mixture_by_id", (PyCFunction)PyAleaSystem_find_mixture_by_id, METH_VARARGS,
     "find_mixture_by_id(mixture_id) -> index or None\n\nFind mixture index by ID."},
    {"mixture_count", (PyCFunction)PyAleaSystem_mixture_count, METH_NOARGS,
     "mixture_count() -> int\n\nGet number of mixtures."},
    {"mixture_get_id", (PyCFunction)PyAleaSystem_mixture_get_id, METH_VARARGS,
     "mixture_get_id(mix_index) -> int\n\nGet mixture ID by index."},
    {"mixture_get_components", (PyCFunction)PyAleaSystem_mixture_get_components, METH_VARARGS,
     "mixture_get_components(mix_index) -> list of {material_id, fraction}\n\n"
     "Get components of a mixture."},

    /* Export */
    {"export_mcnp", (PyCFunction)PyAleaSystem_export_mcnp,
     METH_VARARGS | METH_KEYWORDS,
     "export_mcnp(filename, deduplicate=True, universe_depth=-1, fill_depth=0)\n\n"
     "Low-level MCNP export. The public Python API is Model.export_mcnp(filename).\n\n"
     "Args:\n"
     "  filename: Output file path\n"
     "  deduplicate: Deduplicate surfaces (default True)\n"
     "  universe_depth: Universe filter (-1=all, 0=base only, N=N levels deep)\n"
     "  fill_depth: FILL expansion (0=none, N=N levels, -1=full flatten)"},
    {"export_mcnp_string", (PyCFunction)PyAleaSystem_export_mcnp_string,
     METH_VARARGS | METH_KEYWORDS,
     "export_mcnp_string(deduplicate=True, universe_depth=-1, fill_depth=0) -> str\n\n"
     "Export MCNP text in memory."},
    {"export_openmc", (PyCFunction)PyAleaSystem_export_openmc,
     METH_VARARGS | METH_KEYWORDS,
     "export_openmc(filename)\n\nExport to OpenMC XML format."},
    {"export_openmc_string", (PyCFunction)PyAleaSystem_export_openmc_string,
     METH_NOARGS, "export_openmc_string() -> str\n\nExport OpenMC XML in memory."},
    {"export_serpent", (PyCFunction)PyAleaSystem_export_serpent,
     METH_VARARGS | METH_KEYWORDS,
     "export_serpent(filename)\n\nExport to Serpent input format."},
    {"export_serpent_string", (PyCFunction)PyAleaSystem_export_serpent_string,
     METH_NOARGS, "export_serpent_string() -> str\n\nExport Serpent text in memory."},

    /* Merge */
    {"merge", (PyCFunction)PyAleaSystem_merge, METH_VARARGS | METH_KEYWORDS,
     "merge(other, id_offset=0)\n\nMerge another system into this one."},

    /* Utilities */
    {"set_verbose", (PyCFunction)PyAleaSystem_set_verbose, METH_VARARGS,
     "set_verbose(enabled)\n\nEnable/disable verbose output."},
    {"validate", (PyCFunction)PyAleaSystem_validate, METH_NOARGS,
     "validate() -> int\n\nValidate system integrity."},
    {"print_summary", (PyCFunction)PyAleaSystem_print_summary, METH_NOARGS,
     "print_summary()\n\nPrint system summary."},
    {"set_tolerance", (PyCFunction)PyAleaSystem_set_tolerance,
     METH_VARARGS | METH_KEYWORDS,
     "set_tolerance(abs_tol=1e-6, rel_tol=1e-9, zero_thresh=1e-10)\n\nSet tolerances."},
    {"clone", (PyCFunction)PyAleaSystem_clone, METH_NOARGS,
     "clone() -> System\n\nCreate deep copy of system."},
    {"reset", (PyCFunction)PyAleaSystem_reset, METH_NOARGS,
     "reset()\n\nReset system to empty state."},

    /* Raycast */
    {"raycast", (PyCFunction)PyAleaSystem_raycast, METH_VARARGS | METH_KEYWORDS,
     "raycast(ox, oy, oz, dx, dy, dz, t_max=0) -> {'segments': [...], 'hits': [...]}\n\n"
     "Cast ray through geometry (global surface-BVH pipeline). 'segments' is\n"
     "the precedence-resolved cell list; 'hits' is the raw global crossing\n"
     "list [{t, surface_id}, ...] underneath it."},
    {"classify_ray_intervals", (PyCFunction)PyAleaSystem_classify_ray_intervals, METH_VARARGS,
     "classify_ray_intervals(ox, oy, oz, dx, dy, dz, t_max) -> list of findings\n\n"
     "Classify elementary ray intervals by complete cell ownership:\n"
     "kind 0=ok, 1=gap, 2=overlap (cell_id/overlap_cell_id = claimant pair),\n"
     "3=undefined fill. Catches overlaps in isolation that segments hide."},
    {"ray_first_cell", (PyCFunction)PyAleaSystem_ray_first_cell, METH_VARARGS,
     "ray_first_cell(ox, oy, oz, dx, dy, dz, t_max=0) -> (cell_id, t) or None"},
    {"ray_first_visible", (PyCFunction)PyAleaSystem_ray_first_visible,
     METH_VARARGS | METH_KEYWORDS,
     "ray_first_visible(ox, oy, oz, dx, dy, dz, t_min=0, t_max=0, material_filter=-1,\n"
     "                  include_surface=False, include_normal=False) -> dict or None\n\n"
     "Return the first matching non-void interval using libalea's unified ray query.\n"
     "material_filter < 0 accepts every material; otherwise it selects one material."},
    {"ray_boundary_events", (PyCFunction)PyAleaSystem_ray_boundary_events,
     METH_VARARGS | METH_KEYWORDS,
     "ray_boundary_events(ox, oy, oz, dx, dy, dz, t_min=0, t_max=0, max_events=0,\n"
     "                    max_output_bytes=0, include_primitive_id=False,\n"
     "                    include_normal=False, include_all_coincident_physical=False,\n"
     "                    include_occurrence_provenance=False)\n"
     "-> list[dict]\n\n"
     "Return physical, synthetic-lattice, or unresolved boundary events. The\n"
     "occurrence-provenance option uses hierarchical selected traversal and adds\n"
     "the active universe-local transition frame."},
    {"render_3d", (PyCFunction)PyAleaSystem_render_3d,
     METH_VARARGS | METH_KEYWORDS,
     "render_3d(width=800, height=600, eye=None, target=None, up=None, fov=45,\n"
     "          ortho_height=0, color_by='material', mode='solid', background=None,\n"
     "          shadows=True, edges=False, aa_samples=1,\n"
     "          auxiliary=False, clips=None, clip_mode='and',\n"
     "          material_filter_mode='all', material_ids=None,\n"
     "          cell_filter_mode='all', cell_ids=None) -> {'rgb': ndarray, ...}\n\n"
     "Render the CSG model with libalea's native CPU ray renderer. Set auxiliary\n"
     "to also return depth, cell_ids, material_ids, and world-space normals.\n"
     "clips contains (nx, ny, nz, d) retained half-spaces combined by clip_mode.\n"
     "Cell and material filters also remove hidden geometry from shadows."},

    /* Slice curves API (for matplotlib) */
    {"get_slice_curves_z", (PyCFunction)PyAleaSystem_get_slice_curves_z, METH_VARARGS | METH_KEYWORDS,
     "get_slice_curves_z(z, x_min, x_max, y_min, y_max) -> dict\n\n"
     "Get analytical curves from surface-plane intersections (XY plane at z).\n\n"
     "Returns dict with:\n"
     "  'curves': list of curve dicts (type, surface_id, center/radius/etc)\n"
     "  'u_min', 'u_max', 'v_min', 'v_max': bounding box of curves"},
    {"get_slice_curves_y", (PyCFunction)PyAleaSystem_get_slice_curves_y, METH_VARARGS | METH_KEYWORDS,
     "get_slice_curves_y(y, x_min, x_max, z_min, z_max) -> dict\n\n"
     "Get analytical curves (XZ plane at y)."},
    {"get_slice_curves_x", (PyCFunction)PyAleaSystem_get_slice_curves_x, METH_VARARGS | METH_KEYWORDS,
     "get_slice_curves_x(x, y_min, y_max, z_min, z_max) -> dict\n\n"
     "Get analytical curves (YZ plane at x)."},

    /* Grid cell queries (for matplotlib fills) */
    {"find_cells_grid_z", (PyCFunction)PyAleaSystem_find_cells_grid_z, METH_VARARGS | METH_KEYWORDS,
     "find_cells_grid_z(z, x_min, x_max, y_min, y_max, nx, ny, universe_depth=-1, error_mode='none') -> dict\n\n"
     "Find cells on a 2D grid (XY plane at z).\n\n"
     "Args:\n"
     "  universe_depth: Universe depth (-1=innermost, 0=root, N=depth N)\n"
     "  error_mode: 'none', 'fast', or 'full'\n\n"
     "Returns dict with:\n"
     "  'cell_ids': list of cell IDs (nx*ny, -1 for void)\n"
     "  'material_ids': list of material IDs\n"
     "  'errors': list of error codes (if error_mode != 'none'): 0=ok, 1=overlap, 2=undefined, 3=undefined fill\n"
     "  'coverage', 'secondary_cell_ids', 'error_components': plot diagnostics (if error_mode != 'none')\n"
     "  'error_lines': analytical error segments (if error_mode='full' and available)\n"
     "  'nx', 'ny': grid dimensions"},
    {"find_cells_grid_y", (PyCFunction)PyAleaSystem_find_cells_grid_y, METH_VARARGS | METH_KEYWORDS,
     "find_cells_grid_y(y, x_min, x_max, z_min, z_max, nx, nz, universe_depth=-1, error_mode='none') -> dict\n\n"
     "Find cells on XZ grid at y. Same optional args as find_cells_grid_z."},
    {"find_cells_grid_x", (PyCFunction)PyAleaSystem_find_cells_grid_x, METH_VARARGS | METH_KEYWORDS,
     "find_cells_grid_x(x, y_min, y_max, z_min, z_max, ny, nz, universe_depth=-1, error_mode='none') -> dict\n\n"
     "Find cells on YZ grid at x. Same optional args as find_cells_grid_z."},

    /* Arbitrary plane slice curves and grid */
    {"get_slice_curves", (PyCFunction)PyAleaSystem_get_slice_curves, METH_VARARGS | METH_KEYWORDS,
     "get_slice_curves(origin, normal, up, u_min, u_max, v_min, v_max) -> dict\n\n"
     "Get analytical curves on arbitrary slice plane.\n\n"
     "Args:\n"
     "  origin: Point on the plane as (x, y, z) tuple\n"
     "  normal: Plane normal as (nx, ny, nz) tuple\n"
     "  up: Up vector hint as (ux, uy, uz) tuple\n"
     "  u_min, u_max: Horizontal bounds in plane coordinates\n"
     "  v_min, v_max: Vertical bounds in plane coordinates\n\n"
     "Returns dict with:\n"
     "  'curves': list of curve dicts\n"
     "  'u_min', 'u_max', 'v_min', 'v_max': bounds"},
    {"find_cells_grid", (PyCFunction)PyAleaSystem_find_cells_grid, METH_VARARGS | METH_KEYWORDS,
     "find_cells_grid(origin, normal, up, u_min, u_max, v_min, v_max, nu, nv, universe_depth=-1, error_mode='none') -> dict\n\n"
     "Find cells on arbitrary slice plane grid.\n\n"
     "Args:\n"
     "  origin: Point on the plane as (x, y, z) tuple\n"
     "  normal: Plane normal as (nx, ny, nz) tuple\n"
     "  up: Up vector hint as (ux, uy, uz) tuple\n"
     "  u_min, u_max: Horizontal bounds in plane coordinates\n"
     "  v_min, v_max: Vertical bounds in plane coordinates\n"
     "  nu, nv: Grid resolution\n"
     "  universe_depth: Universe depth (-1=innermost, 0=root, N=depth N)\n"
     "  error_mode: 'none', 'fast', or 'full'\n\n"
     "Returns dict with:\n"
     "  'cell_ids': list of cell IDs (nu*nv, -1 for void)\n"
     "  'material_ids': list of material IDs\n"
     "  'errors': list of error codes (if error_mode != 'none')\n"
     "  'coverage', 'secondary_cell_ids', 'error_components': plot diagnostics (if error_mode != 'none')\n"
     "  'error_lines': analytical error segments (if error_mode='full' and available)\n"
     "  'nu', 'nv': grid dimensions"},
    {"find_local_coverage_components", (PyCFunction)PyAleaSystem_find_local_coverage_components,
     METH_VARARGS | METH_KEYWORDS,
     "find_local_coverage_components(origin, normal, up, u_min, u_max, v_min, v_max, nu, nv, universe_depth=-1, max_pixels=0, max_scratch_bytes=0) -> dict\n\n"
     "Run exact local coverage and return compact connected components only. "
     "No pixel grid or analytical error-line payload is retained."},
    {"refine_grid_coverage_tiles", (PyCFunction)PyAleaSystem_refine_grid_coverage_tiles,
     METH_VARARGS | METH_KEYWORDS,
     "refine_grid_coverage_tiles(origin, normal, up, u_min, u_max, v_min, v_max,\n"
     "                           nu, nv, coverage, errors, secondary_cell_ids, *,\n"
     "                           universe_depth=-1, tile_w=16, tile_h=16,\n"
     "                           max_candidates=4096, max_evaluated_candidates=4096,\n"
     "                           max_exact_fallback_pixels=0,\n"
     "                           chain_candidates=False, tile_mask=None) -> dict\n\n"
     "Run bounded tile coverage refinement on an existing fast grid. Returns\n"
     "updated coverage/error arrays plus stats; when stats['incomplete'] is true,\n"
     "skipped pixels remain provisional and must not be treated as clean."},
    {"refine_grid_coverage_paths", (PyCFunction)PyAleaSystem_refine_grid_coverage_paths,
     METH_VARARGS | METH_KEYWORDS,
     "refine_grid_coverage_paths(origin, normal, up, u_min, u_max, v_min, v_max,\n"
     "                           nu, nv, cell_ids, coverage, errors,\n"
     "                           secondary_cell_ids, *, universe_depth=-1,\n"
     "                           tile_w=16, tile_h=16, max_candidates=4096,\n"
     "                           max_exact_fallback_pixels=0,\n"
     "                           use_path_2d_index=False, tile_mask=None) -> dict\n\n"
     "Rebuild compact concrete path IDs and refine selected tiles within each\n"
     "path's local universe. Returns provisional pixels and bounded-work stats."},

    /* Label positioning */
    {"find_label_positions", (PyCFunction)PyAleaSystem_find_label_positions, METH_VARARGS | METH_KEYWORDS,
     "find_label_positions(ids, width, height, min_pixels=100) -> list\n\n"
     "Find optimal label positions for regions in a cell/material grid.\n\n"
     "Args:\n"
     "  ids: List of cell or material IDs from find_cells_grid_*()\n"
     "  width, height: Grid dimensions\n"
     "  min_pixels: Minimum region size to include (default 100)\n\n"
     "Returns list of dicts with:\n"
     "  'id': Cell/material ID\n"
     "  'px', 'py': Pixel coordinates for label\n"
     "  'pixel_count': Region size in pixels"},

    /* Cell filtering */
    {"get_cells_by_material", (PyCFunction)PyAleaSystem_get_cells_by_material, METH_VARARGS,
     "get_cells_by_material(material_id) -> list of indices\n\n"
     "Get indices of all cells with the given material ID."},
    {"get_cells_by_universe", (PyCFunction)PyAleaSystem_get_cells_by_universe, METH_VARARGS,
     "get_cells_by_universe(universe_id) -> list of indices\n\n"
     "Get indices of all cells in the given universe."},
    {"get_cells_filling_universe", (PyCFunction)PyAleaSystem_get_cells_filling_universe, METH_VARARGS,
     "get_cells_filling_universe(universe_id) -> list of indices\n\n"
     "Get indices of all cells that FILL the given universe."},
    {"get_cells_in_bbox", (PyCFunction)PyAleaSystem_get_cells_in_bbox, METH_VARARGS,
     "get_cells_in_bbox(x_min, x_max, y_min, y_max, z_min, z_max) -> list of indices\n\n"
     "Get indices of cells whose bounding boxes intersect the given region."},

    /* Extract operations */
    {"extract_universe", (PyCFunction)PyAleaSystem_extract_universe, METH_VARARGS,
     "extract_universe(universe_id) -> System\n\n"
     "Extract a universe and all universes it references into a new system."},
    {"extract_region", (PyCFunction)PyAleaSystem_extract_region, METH_VARARGS,
     "extract_region(x_min, x_max, y_min, y_max, z_min, z_max) -> System\n\n"
     "Extract cells in a bounding box region into a new system."},

    /* Material operations */
    {"create_mixture", (PyCFunction)PyAleaSystem_create_mixture, METH_VARARGS | METH_KEYWORDS,
     "create_mixture(material_ids, fractions, new_id=0) -> int\n\n"
     "Create a mixture of materials.\n\n"
     "Args:\n"
     "  material_ids: List of material IDs to mix\n"
     "  fractions: List of fractions (normalized automatically)\n"
     "  new_id: ID for new mixture (0 for auto-assign)\n\n"
     "Returns: Assigned material ID"},

    /* Cell fill / ID operations */
    {"add_transform", (PyCFunction)PyAleaSystem_add_transform,
     METH_VARARGS | METH_KEYWORDS,
     "add_transform(transform_id, values, degrees=False) -> int\n\n"
     "Add or replace a named MCNP transform and return transform_id. "
     "values accepts a translation or full/partial MCNP transform sequence."},
    {"add_inline_transform", (PyCFunction)PyAleaSystem_add_inline_transform,
     METH_VARARGS | METH_KEYWORDS,
     "add_inline_transform(values, degrees=False, cell_id=0, role='fill') -> int\n\n"
     "Add a deduplicated inline MCNP transform and return its assigned ID. "
     "cell_id and role provide diagnostic source context."},
    {"set_fill", (PyCFunction)PyAleaSystem_set_fill, METH_VARARGS,
     "set_fill(cell_index, fill_universe, transform=0)\n\nSet fill universe for a cell."},
    {"set_comment", (PyCFunction)PyAleaSystem_set_comment, METH_VARARGS,
     "set_comment(cell_index, comment)\n\nSet cell comment (None to clear)."},
    {"set_inline_comment", (PyCFunction)PyAleaSystem_set_inline_comment, METH_VARARGS,
     "set_inline_comment(cell_index, comment)\n\nSet cell inline comment (None to clear)."},
    {"cell_set_material", (PyCFunction)PyAleaSystem_cell_set_material, METH_VARARGS,
     "cell_set_material(cell_index, material_index)\n\nSet cell material index (-1 for void)."},
    {"cell_set_density", (PyCFunction)PyAleaSystem_cell_set_density, METH_VARARGS,
     "cell_set_density(cell_index, density)\n\nSet cell density (signed: negative=g/cm3, positive=atoms/b-cm)."},
    {"cell_set_universe", (PyCFunction)PyAleaSystem_cell_set_universe, METH_VARARGS,
     "cell_set_universe(cell_index, universe_id)\n\nSet cell universe membership."},
    {"cell_remove", (PyCFunction)PyAleaSystem_cell_remove, METH_VARARGS,
     "cell_remove(cell_index)\n\nRemove cell by index."},
    {"get_cell_id", (PyCFunction)PyAleaSystem_get_cell_id, METH_VARARGS,
     "get_cell_id(cell_index) -> int\n\nGet MCNP cell ID from cell index."},
    {"cells_in_universe", (PyCFunction)PyAleaSystem_cells_in_universe, METH_VARARGS,
     "cells_in_universe(universe_id) -> list of indices\n\nGet cell indices in a universe."},
    {"find_cell_at", (PyCFunction)PyAleaSystem_find_cell_at, METH_VARARGS,
     "find_cell_at(x, y, z) -> (cell_id, material_id) or None\n\nFind cell and get both cell ID and material."},

    /* Surface operations */
    {"surface_find", (PyCFunction)PyAleaSystem_surface_find, METH_VARARGS,
     "surface_find(surface_id) -> int or None\n\nFind surface index by MCNP surface ID."},
    {"surface_node", (PyCFunction)PyAleaSystem_surface_node, METH_VARARGS,
     "surface_node(surface_id, sense) -> node_id or None\n\nGet node by surface ID and sense (-1 or +1)."},
    {"surface_project_along", (PyCFunction)PyAleaSystem_surface_project_along,
     METH_VARARGS,
     "surface_project_along(surface_id, point, direction) -> dict or None\n\n"
     "Project onto the exact surface along the nearest signed line intersection."},
    {"surface_id_at", (PyCFunction)PyAleaSystem_surface_id_at, METH_VARARGS,
     "surface_id_at(idx) -> int\n\nGet the surface ID at storage index `idx` in [0, surface_count)."},
    {"get_surface_ids", (PyCFunction)PyAleaSystem_get_surface_ids, METH_NOARGS,
     "get_surface_ids() -> list[int]\n\nReturn all surface IDs in storage order — single C call, no tree walk."},

    /* CSG Node Inspection (individual) */
    {"node_operation", (PyCFunction)PyAleaSystem_node_operation, METH_VARARGS,
     "node_operation(node_id) -> int\n\nGet operation type: 0=primitive, 1=union, 2=intersection, 3=difference, 4=complement."},
    {"node_left", (PyCFunction)PyAleaSystem_node_left, METH_VARARGS,
     "node_left(node_id) -> node_id or None\n\nGet left child of a boolean node."},
    {"node_right", (PyCFunction)PyAleaSystem_node_right, METH_VARARGS,
     "node_right(node_id) -> node_id or None\n\nGet right child of a boolean node."},
    {"node_primitive_type", (PyCFunction)PyAleaSystem_node_primitive_type, METH_VARARGS,
     "node_primitive_type(node_id) -> int\n\nGet primitive type for a leaf node."},
    {"node_primitive_id", (PyCFunction)PyAleaSystem_node_primitive_id, METH_VARARGS,
     "node_primitive_id(node_id) -> int or None\n\nGet primitive ID for a leaf node."},
    {"node_sense", (PyCFunction)PyAleaSystem_node_sense, METH_VARARGS,
     "node_sense(node_id) -> int\n\nGet sense for a primitive node (+1 or -1)."},
    {"node_surface_id", (PyCFunction)PyAleaSystem_node_surface_id, METH_VARARGS,
     "node_surface_id(node_id) -> int\n\nGet MCNP surface ID associated with a primitive node."},

    /* Cell expression */
    {"cell_expr", (PyCFunction)PyAleaSystem_cell_expr, METH_VARARGS | METH_KEYWORDS,
     "cell_expr(cell_index, union_op=':', inter_op=' ', compl_op='#') -> str\n\n"
     "Get CSG expression string for a cell.\n"
     "Default operators produce MCNP-style output. Use union_op=' | ', compl_op='~' for OpenMC-style."},
    {"cell_surface_ids", (PyCFunction)PyAleaSystem_cell_surface_ids,
     METH_VARARGS | METH_KEYWORDS,
     "cell_surface_ids(cell_index, max_surfaces=32, max_nodes=256) -> dict\n\n"
     "Return bounded unique primitive surface IDs for one cell CSG root."},

    /* Renumbering */
    {"renumber_cells", (PyCFunction)PyAleaSystem_renumber_cells, METH_VARARGS,
     "renumber_cells(start_id) -> int\n\nRenumber all cells starting from start_id."},
    {"renumber_surfaces", (PyCFunction)PyAleaSystem_renumber_surfaces, METH_VARARGS,
     "renumber_surfaces(start_id) -> int\n\nRenumber all surfaces starting from start_id."},
    {"offset_cell_ids", (PyCFunction)PyAleaSystem_offset_cell_ids, METH_VARARGS,
     "offset_cell_ids(offset)\n\nAdd offset to all cell IDs."},
    {"offset_surface_ids", (PyCFunction)PyAleaSystem_offset_surface_ids, METH_VARARGS,
     "offset_surface_ids(offset)\n\nAdd offset to all surface IDs."},
    {"offset_material_ids", (PyCFunction)PyAleaSystem_offset_material_ids, METH_VARARGS,
     "offset_material_ids(offset)\n\nAdd offset to all material IDs."},

    /* Split / Expand */
    {"split_union_cells", (PyCFunction)PyAleaSystem_split_union_cells, METH_NOARGS,
     "split_union_cells() -> int\n\nSplit cells with top-level unions into multiple simpler cells.\n"
     "Returns number of new cells created."},
    {"expand_macrobodies", (PyCFunction)PyAleaSystem_expand_macrobodies, METH_NOARGS,
     "expand_macrobodies() -> int\n\nExpand all macrobodies to primitive surfaces."},
    {"carve_universe", (PyCFunction)PyAleaSystem_carve_universe,
     METH_VARARGS | METH_KEYWORDS,
     "carve_universe(universe_id, carve_root, simplify=True, cell_limit=-1) -> dict\n\n"
     "Subtract a CSG region from ordinary cells in one universe."},

    /* Volume estimation */
    {"compute_bounding_sphere", (PyCFunction)PyAleaSystem_compute_bounding_sphere, METH_VARARGS,
     "compute_bounding_sphere(tol=1.0) -> (cx, cy, cz, radius)\n\n"
     "Compute a tight bounding sphere for the entire model."},
    {"estimate_volumes", (PyCFunction)PyAleaSystem_estimate_volumes,
     METH_VARARGS | METH_KEYWORDS,
     "estimate_volumes(n_rays=100000, seed=42, workers=0, "
     "target_rel_error=None, max_rays=None, batch_size=10000, "
     "progress=None, rng='philox4x32-10') -> dict\n\n"
     "Estimate physical volumes per concrete hierarchical placement.\n"
     "Returns dict with 'volumes', 'rel_errors'\n"
     "and 'paths' (one path-identity dict per entry)."},
    {"estimate_cell_volume", (PyCFunction)PyAleaSystem_estimate_cell_volume,
     METH_VARARGS | METH_KEYWORDS,
     "estimate_cell_volume(cell_index, bounds=None, relative_tolerance=1e-3,\n"
     "                     absolute_tolerance=0, max_depth=10, min_size=0,\n"
     "                     samples_per_axis=2, workers=0,\n"
     "                     max_parallel_scratch_bytes=67108864) -> dict\n\n"
     "Estimate one universe-local cell CSG volume with a deterministic octree."},
    {"volume_path_at_point", (PyCFunction)PyAleaSystem_volume_path_at_point, METH_VARARGS,
     "volume_path_at_point(x, y, z) -> dict or None\n\n"
     "Resolve a world-space point to its concrete hierarchical volume path."},
    {"volume_path_resolve_at_point", (PyCFunction)PyAleaSystem_volume_path_resolve_at_point, METH_VARARGS,
     "volume_path_resolve_at_point(x, y, z) -> dict or None\n\n"
     "Resolve one structural hierarchy path without global path enumeration."},
    {"volume_path_resolve_cell_at_point",
     (PyCFunction)PyAleaSystem_volume_path_resolve_cell_at_point, METH_VARARGS,
     "volume_path_resolve_cell_at_point(x, y, z, cell_id, universe_id) -> dict\n\n"
     "Resolve zero, one, or multiple occurrences of a specified hierarchy cell."},
    {"volume_path_resolve_cell_point_sets",
     (PyCFunction)PyAleaSystem_volume_path_resolve_cell_point_sets,
     METH_VARARGS | METH_KEYWORDS,
     "volume_path_resolve_cell_point_sets(point_sets, requested_workers=0, "
     "max_parallel_scratch_bytes=67108864) -> dict\n\n"
     "Resolve and aggregate independent point sets for specified hierarchy cells."},
    {"volume_path_resolve_cell_from_transform_evidence",
     (PyCFunction)PyAleaSystem_volume_path_resolve_cell_from_transform_evidence,
     METH_VARARGS,
     "volume_path_resolve_cell_from_transform_evidence(cell_id, universe_id, "
     "world_point, world_direction, local_point, local_direction, "
     "world_point_precision, world_direction_precision, local_point_precision, "
     "local_direction_precision) -> dict\n\n"
     "Resolve a hierarchy occurrence from corresponding coordinate frames "
     "without requiring cell containment."},
    {"query_acceleration_stats", (PyCFunction)PyAleaSystem_query_acceleration_stats, METH_NOARGS,
     "query_acceleration_stats() -> dict\n\n"
     "Inspect the built query-acceleration structure.\n"
     "Includes hierarchical placement/BLAS/transform/memory counts."},
    {"surface_reference_stats", (PyCFunction)PyAleaSystem_surface_reference_stats, METH_NOARGS,
     "surface_reference_stats() -> dict\n\n"
     "Inspect the exact MCNP surface-card to cell-reference CSR."},
    {"surface_cell_references", (PyCFunction)PyAleaSystem_surface_cell_references, METH_VARARGS,
     "surface_cell_references(surface_id) -> list[dict]\n\n"
     "Return exact MCNP surface-card references from the reverse index."},
    /* BBox tightening */
    {"tighten_cell_bbox", (PyCFunction)PyAleaSystem_tighten_cell_bbox, METH_VARARGS,
     "tighten_cell_bbox(cell_index, tol=1.0) -> (xmin, xmax, ymin, ymax, zmin, zmax)\n\n"
     "Tighten a single cell's bounding box via interval arithmetic."},
    {"tighten_all_bboxes", (PyCFunction)PyAleaSystem_tighten_all_bboxes, METH_VARARGS,
     "tighten_all_bboxes(tol=1.0) -> int\n\n"
     "Tighten all cell bounding boxes. Returns number of cells tightened."},

    /* Cell-aware raycast */
    {"raycast_cell_aware", (PyCFunction)PyAleaSystem_raycast_cell_aware,
     METH_VARARGS | METH_KEYWORDS,
     "raycast_cell_aware(ox, oy, oz, dx, dy, dz, t_max=0) -> list of segments\n\n"
     "Cell-aware raycast using per-cell surface index. More efficient than global raycast."},
    {"raycast_hier_fast", (PyCFunction)PyAleaSystem_raycast_hier_fast,
     METH_VARARGS | METH_KEYWORDS,
     "raycast_hier_fast(ox, oy, oz, dx, dy, dz, t_max=0, include_paths=False) -> list of segments\n\n"
     "Fast hierarchical segment tracer. Walks the hierarchical spatial index and\n"
     "tests only surfaces of cells along the ray (not the global surface BVH).\n"
     "Dramatically faster on large models; the default path used by Model.trace().\n"
     "Set include_paths=True to attach resolved hierarchy paths to segments."},
    {"raycast_hier_fast_batch", (PyCFunction)PyAleaSystem_raycast_hier_fast_batch,
     METH_VARARGS | METH_KEYWORDS,
     "raycast_hier_fast_batch(origins, directions, t_max=0, include_paths=False) -> list[list[dict]]\n\n"
     "Trace independent rays in parallel using the hierarchical segment tracer."},
    {"trace_ray_slice_compact", (PyCFunction)PyAleaSystem_trace_ray_slice_compact,
     METH_VARARGS | METH_KEYWORDS,
     "trace_ray_slice_compact(origin, normal, up, u_min, u_max, v_min, v_max, row_count,\n"
     "                        projected_depth=-1, include_paths=False, max_segments=0,\n"
     "                        max_path_entries=0, max_output_bytes=0) -> dict[str, numpy.ndarray]\n\n"
     "Trace an arbitrary-plane ray slice into compact CSR arrays. Distances are\n"
     "clipped view-U coordinates; row_offsets indexes increasing-V rows."},
    {"trace_ray_slice_grid", (PyCFunction)PyAleaSystem_trace_ray_slice_grid,
     METH_VARARGS | METH_KEYWORDS,
     "trace_ray_slice_grid(origin, normal, up, u_min, u_max, v_min, v_max, nu, nv,\n"
     "                     projected_depth=-1, max_segments=0, max_output_bytes=0) -> dict\n\n"
     "Trace and rasterize a ray slice into fresh caller-owned NumPy buffers.\n"
     "The result follows the viewer grid-field contract."},
    {"validate_ray_slice_compact", (PyCFunction)PyAleaSystem_validate_ray_slice_compact,
     METH_VARARGS | METH_KEYWORDS,
     "validate_ray_slice_compact(origin, normal, up, u_min, u_max, v_min, v_max, row_count,\n"
     "                           projected_depth=-1, include_agreements=False,\n"
     "                           max_trace_intervals=0, max_path_entries=0,\n"
     "                           max_output_intervals=0, max_output_bytes=0,\n"
     "                           include_provenance=False, cache_width=0,\n"
     "                           bidirectional=True, coverage=False, coverage_domain=None,\n"
     "                           coverage_unowned_is_exterior=False,\n"
     "                           coverage_report_exterior=False, max_coverage_owners=0,\n"
     "                           max_coverage_rows=0, refinement_depth=0, refine_signals=0,\n"
     "                           refine_displacement=0.0, refine_crossing_density=0,\n"
     "                           min_transverse_spacing=0.0) -> dict\n\n"
     "Run native ray-slice validation and return diagnostics plus the forward compact\n"
     "render result. Set include_provenance=True to use a directional cache and include\n"
     "canonical endpoint surface provenance.\n\n"
     "The bidirectional check is trace-consistency evidence only. Only coverage=True\n"
     "classifies gaps and overlaps, and it requires exactly one unowned-space policy:\n"
     "coverage_domain=(u_min, u_max) or coverage_unowned_is_exterior=True.\n"
     "refinement_depth>0 publishes refined rows between the requested ones, so the\n"
     "returned row_count exceeds the requested row_count; read row_base_indices\n"
     "(-1 marks a refined row) instead of assuming a one-to-one row mapping."},
    {"ray_coverage_slice", (PyCFunction)PyAleaSystem_ray_coverage_slice,
     METH_VARARGS | METH_KEYWORDS,
     "ray_coverage_slice(origins, directions, t_max, direction_tags=None,\n"
     "                   transverse_coordinates=None, domain=None, report_exterior=False,\n"
     "                   max_rows=0, max_intervals=0, max_owners=0, max_output_bytes=0,\n"
     "                   max_refinement_depth=0, refine_signals=0, endpoint_displacement=0.0,\n"
     "                   crossing_density=0, min_transverse_spacing=0.0) -> dict\n\n"
     "Complete ownership coverage for packed (n_rays, 3) origins and directions,\n"
     "published as input-order CSR: row_offsets indexes intervals, owner_offsets\n"
     "indexes every concrete owner occurrence claiming an interval. Unlike the tracing\n"
     "methods, which report one selected owner, this is what reveals gaps and overlaps.\n"
     "kinds are 0=unique, 1=gap, 2=allowed_exterior, 3=overlap, 4=undefined_fill,\n"
     "5=unresolved, 6=truncated. domain=(t_min, t_max) marks where ownership is\n"
     "required; report_exterior requires domain. TRUNCATED owners are a retained\n"
     "prefix: use owner_count_lower_bounds, not owner_offsets, for the known count."},

    /* Geometry validator */
    {"validate_geometry", (PyCFunction)PyAleaSystem_validate_geometry,
     METH_VARARGS | METH_KEYWORDS,
     "validate_geometry(options=None) -> dict\n\n"
     "Transport-style geometry validation with auto-generated rays.\n"
     "options is a dict: ray_count, seed, universe_depth, max_errors,\n"
     "max_samples_per_signature, max_samples_per_curve, max_crossings,\n"
     "sample_offset, t_max, strict, allow_exterior_void, domain_bounds.\n"
     "domain_bounds is (min_x, max_x, min_y, max_y, min_z, max_z); unowned space inside\n"
     "it is reported as an interior_gap error.\n"
     "Returns {errors, crossings_checked, adjacency_hits, exact_queries, ambiguous_crossings, truncated}."},
    {"validate_geometry_ray", (PyCFunction)PyAleaSystem_validate_geometry_ray,
     METH_VARARGS | METH_KEYWORDS,
     "validate_geometry_ray(ox, oy, oz, dx, dy, dz, t_max=0, options=None) -> dict\n\n"
     "Validate geometry along a single ray."},
    {"validate_geometry_slice", (PyCFunction)PyAleaSystem_validate_geometry_slice,
     METH_VARARGS | METH_KEYWORDS,
     "validate_geometry_slice(origin, normal, up, u_min, u_max, v_min, v_max, options=None) -> dict\n\n"
     "Validate geometry by sampling the analytical boundary curves on a slice plane."},
    {"check_transition", (PyCFunction)PyAleaSystem_check_transition,
     METH_VARARGS | METH_KEYWORDS,
     "check_transition(universe, current_cell, surface, point, direction, tied_surfaces=None, "
     "probe_distance=0, max_probe_distance=0, max_coverage_hits=256) -> dict\n\n"
     "Classify a supplied universe-local boundary transition witness."},
    {"transition_slice_screen",
     (PyCFunction)PyAleaSystem_transition_slice_screen,
     METH_VARARGS | METH_KEYWORDS,
     "transition_slice_screen(origin, normal, up, u_min, u_max, v_min, v_max, "
     "options=None) -> dict\n\n"
     "Stream bounded horizontal/vertical slice rays, discard valid transitions, "
     "and retain only transition or requested all-owner coverage findings with "
     "bounded component links, resource, and completion statistics."},
    {"transition_slice_screen_batch",
     (PyCFunction)PyAleaSystem_transition_slice_screen_batch,
     METH_VARARGS | METH_KEYWORDS,
     "transition_slice_screen_batch(origin, normal, up, bounds, options=None, "
     "requested_workers=0, max_parallel_scratch_bytes=0) -> dict\n\n"
     "Screen bounded page views with deterministic parallel ordinal output."},
    {"slice_error_page", (PyCFunction)PyAleaSystem_slice_error_page,
     METH_VARARGS | METH_KEYWORDS,
     "slice_error_page(origin, normal, up, view_bounds, required_bounds, "
     "tile_columns=1, tile_rows=1, page_index=0, options=None) -> dict\n\n"
     "Return verified intervals and defect regions for a supported slice tile, "
     "or an unresolved tile with its receipt."},

    /* Grid overlap check */
    {"check_grid_overlaps", (PyCFunction)PyAleaSystem_check_grid_overlaps,
     METH_VARARGS | METH_KEYWORDS,
     "check_grid_overlaps(origin, normal, up, u_min, u_max, v_min, v_max, nu, nv, cell_ids, errors, universe_depth=-1) -> list\n\n"
     "Check grid for overlapping cells (comprehensive). Returns updated error list."},

    /* Surface label positions */
    {"find_surface_label_positions", (PyCFunction)PyAleaSystem_find_surface_label_positions,
     METH_VARARGS | METH_KEYWORDS,
     "find_surface_label_positions(origin, normal, up, u_min, u_max, v_min, v_max, width, height, margin=20, boundary_ids=None) -> list\n\n"
     "Find label positions for surfaces on a slice plane, optionally filtered "
     "against a rendered cell/material boundary grid."},
    {"find_surface_labels_on_boundary_map", (PyCFunction)PyAleaSystem_find_surface_labels_on_boundary_map,
     METH_VARARGS | METH_KEYWORDS,
     "find_surface_labels_on_boundary_map(origin, normal, up, u_min, u_max, v_min, v_max, width, height, grid_ids, margin=20, boundary_by='cell', universe_depth=-1) -> list\n\n"
     "Find label positions from traced surface provenance for the transitions "
     "in a rendered grid. One label per connected arc, with its edge_count."},
    {"find_surface_labels_sparse_on_grid", (PyCFunction)PyAleaSystem_find_surface_labels_sparse_on_grid,
     METH_VARARGS | METH_KEYWORDS,
     "find_surface_labels_sparse_on_grid(origin, normal, up, u_min, u_max, v_min, v_max, width, height, grid_ids, margin=20, boundary_by='cell', universe_depth=-1, max_queries=128, max_labels=256) -> list\n\n"
     "Find ranked surface labels using a bounded sample of exact changed-edge "
     "provenance queries, without analytical curves or a full boundary map."},
    {"sparse_surface_label_stats", (PyCFunction)PyAleaSystem_sparse_surface_label_stats,
     METH_NOARGS,
     "sparse_surface_label_stats() -> dict\n\n"
     "Counters from the most recent sparse surface-label query."},

    /* Config */
    {"get_config", (PyCFunction)PyAleaSystem_get_config, METH_NOARGS,
     "get_config() -> dict\n\nGet current system configuration."},
    {"set_config", (PyCFunction)PyAleaSystem_set_config, METH_VARARGS,
     "set_config(dict)\n\nUpdate system configuration from dict. Only provided keys are changed."},
    {"set_log_level", (PyCFunction)PyAleaSystem_set_log_level, METH_VARARGS,
     "set_log_level(level)\n\nSet log level: 0=none, 1=error, 2=warn, 3=info, 4=debug, 5=trace."},

    /* CSG Simplification */
    {"simplify_all_cells", (PyCFunction)PyAleaSystem_simplify_all_cells, METH_NOARGS,
     "simplify_all_cells() -> dict\n\nSimplify all CSG cells and remove proven-empty cells.\n"
     "Returns dict with simplification statistics."},
    {"simplify_cell_proven", (PyCFunction)PyAleaSystem_simplify_cell_proven,
     METH_VARARGS | METH_KEYWORDS,
     "simplify_cell_proven(cell_index, bounds=None, apply=False, max_depth=12,\n"
     "                     max_nodes=250000, max_patterns=64,\n"
     "                     max_candidates=256, workers=0,\n"
     "                     max_parallel_scratch_bytes=67108864) -> dict\n\n"
     "Discover a simpler cell expression and apply it only after an exact proof."},
    {"simplify_cells_proven", (PyCFunction)PyAleaSystem_simplify_cells_proven,
     METH_VARARGS | METH_KEYWORDS,
     "simplify_cells_proven(cell_indices, bounds=None, apply=False, max_depth=12,\n"
     "                      max_nodes=250000, max_patterns=64,\n"
     "                      max_candidates=256, workers=0,\n"
     "                      max_parallel_scratch_bytes=67108864) -> dict\n\n"
     "Analyze cells in parallel and transactionally apply exact simplifications."},

    /* Numerical BBox tightening */
    {"tighten_cell_bbox_numerical", (PyCFunction)PyAleaSystem_tighten_cell_bbox_numerical, METH_VARARGS,
     "tighten_cell_bbox_numerical(cell_index) -> None\n\n"
     "Tighten a cell's bounding box using numerical sampling (fallback for complex cells)."},

    /* Primitive data */
    {"node_primitive_data", (PyCFunction)PyAleaSystem_node_primitive_data, METH_VARARGS,
     "node_primitive_data(node_id) -> dict\n\n"
     "Get full primitive geometry data for a leaf node.\n"
     "Returns dict with 'type' and type-specific fields (center, radius, coefficients, etc.)."},

    /* Mesh module */
    {"mesh_export", (PyCFunction)PyAleaSystem_mesh_export, METH_VARARGS | METH_KEYWORDS,
     "mesh_export(filename, nx=10, ny=10, nz=10, ...) -> None\n\n"
     "Export geometry as structured mesh (Gmsh or VTK).\n"
     "Optional kwargs: x_min, x_max, y_min, y_max, z_min, z_max (auto if 0),\n"
     "format ('gmsh'/'vtk'), void_material_id, auto_pad."},
    {"mesh_sample", (PyCFunction)PyAleaSystem_mesh_sample, METH_VARARGS | METH_KEYWORDS,
     "mesh_sample(nx=10, ny=10, nz=10, ...) -> dict\n\n"
     "Sample geometry on structured mesh. Returns dict with:\n"
     "  'material_ids', 'cell_ids': flat lists (nx*ny*nz, Z-major)\n"
     "  'x_nodes', 'y_nodes', 'z_nodes': node positions (n+1 each)\n"
     "  'nx', 'ny', 'nz': grid dimensions."},
    {"mesh_sample_configured", (PyCFunction)PyAleaSystem_mesh_sample_configured, METH_O,
     "mesh_sample_configured(options) -> dict\n\n"
     "Sample a structured mesh with configurable sampling and retained fields."},
    {"mesh_export_configured", (PyCFunction)PyAleaSystem_mesh_export_configured, METH_VARARGS,
     "mesh_export_configured(filename, options) -> None\n\n"
     "Sample and export a structured mesh with configurable diagnostic fields."},
    {"mesh_sample_adaptive", (PyCFunction)PyAleaSystem_mesh_sample_adaptive, METH_O,
     "mesh_sample_adaptive(options) -> dict\n\n"
     "Build a nonconforming adaptive octree mesh."},
    {"mesh_export_adaptive", (PyCFunction)PyAleaSystem_mesh_export_adaptive, METH_VARARGS,
     "mesh_export_adaptive(filename, options) -> None\n\n"
     "Build and export a nonconforming adaptive octree mesh."},

    {NULL}
};

/* ============================================================================
 * PyAleaSystem Type Definition
 * ============================================================================ */

static PyTypeObject PyAleaSystemType = {
    PyVarObject_HEAD_INIT(NULL, 0)
    .tp_name = "pyalea._alea.System",
    .tp_doc = PyDoc_STR("pyAlea geometry system.\n\n"
                        "Create an empty system with System(), or use load_mcnp() to load a file."),
    .tp_basicsize = sizeof(PyAleaSystemObject),
    .tp_itemsize = 0,
    .tp_flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE,
    .tp_new = PyAleaSystem_new,
    .tp_init = (initproc)PyAleaSystem_init,
    .tp_dealloc = (destructor)PyAleaSystem_dealloc,
    .tp_methods = PyAleaSystem_methods,
    .tp_getset = PyAleaSystem_getsetters,
};

/* ============================================================================
 * VoidResult Type + Module Init (included files)
 * ============================================================================ */

#include "_bind_void.c"
#include "_bind_nucdata.c"
#include "_bind_transport.c"
#include "_bind_module.c"
