// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/* This file is part of pyalea_binding.c — do NOT compile separately.
 * It is #included from pyalea_binding.c to keep a single translation unit.
 *
 * Contents: Nuclear data Python types (XsDir, Nuclide, NucMaterial, Multigroup)
 *           and module-level functions for the nucdata API.
 */

/* ============================================================================
 * Forward Declarations
 * ============================================================================ */

static PyTypeObject PyAleaXsDirType;
static PyTypeObject PyAleaNuclideType;
static PyTypeObject PyAleaThermalType;
static PyTypeObject PyAleaNucMaterialType;
static PyTypeObject PyAleaMultigroupType;

/* ============================================================================
 * XsDir Python Type
 * ============================================================================ */

typedef struct PyAleaPhotonAuditCacheEntry {
    const alea_nuc_nuclide_t* nuc;
    alea_error_t error;
    alea_nuc_photon_production_audit_t report;
    struct PyAleaPhotonAuditCacheEntry* next;
} PyAleaPhotonAuditCacheEntry;

typedef struct {
    PyObject_HEAD
    alea_nuc_xsdir_t* xsdir;
    PyAleaPhotonAuditCacheEntry* photon_audits;
} PyAleaXsDirObject;

static PyObject* PyAleaXsDir_entry_dict(
        const alea_nuc_xsdir_entry_t* entry) {
    if (!entry) Py_RETURN_NONE;
    return Py_BuildValue("{s:s, s:d, s:s, s:i, s:i, s:d}",
        "zaid", entry->zaid,
        "awr", entry->awr,
        "filename", entry->filename,
        "file_type", entry->file_type,
        "address", entry->address,
        "temperature", entry->temperature);
}

static void PyAleaXsDir_clear(PyAleaXsDirObject* self) {
    while (self->photon_audits) {
        PyAleaPhotonAuditCacheEntry* next = self->photon_audits->next;
        free(self->photon_audits);
        self->photon_audits = next;
    }
    if (self->xsdir) {
        alea_nuc_xsdir_free(self->xsdir);
        self->xsdir = NULL;
    }
}

static void PyAleaXsDir_dealloc(PyAleaXsDirObject* self) {
    PyAleaXsDir_clear(self);
    Py_TYPE(self)->tp_free((PyObject*)self);
}

static PyObject* PyAleaXsDir_new(PyTypeObject* type, PyObject* args, PyObject* kwds) {
    (void)args; (void)kwds;
    PyAleaXsDirObject* self = (PyAleaXsDirObject*)type->tp_alloc(type, 0);
    if (self) {
        self->xsdir = NULL;
        self->photon_audits = NULL;
    }
    return (PyObject*)self;
}

static int PyAleaXsDir_init(PyAleaXsDirObject* self, PyObject* args, PyObject* kwds) {
    static char* kwlist[] = {"path", "directory", NULL};
    const char* path = NULL;
    int directory = 0;

    if (!PyArg_ParseTupleAndKeywords(args, kwds, "s|p", kwlist, &path, &directory))
        return -1;

    /* Cached Nuclide objects borrow entries owned by this directory. Replacing
     * it in place would invalidate those objects, so reject explicit re-init. */
    if (self->xsdir) {
        PyErr_SetString(PyExc_RuntimeError, "XsDir is already initialized");
        return -1;
    }

    alea_nuc_xsdir_t* xsdir = directory
        ? alea_nuc_xsdir_load_dir(path)
        : alea_nuc_xsdir_load(path);

    if (!xsdir) {
        PyErr_Format(PyExc_IOError, "Failed to load xsdir from '%s'", path);
        return -1;
    }

    self->xsdir = xsdir;
    return 0;
}

static PyObject* PyAleaXsDir_find(PyAleaXsDirObject* self, PyObject* args) {
    const char* zaid;
    if (!PyArg_ParseTuple(args, "s", &zaid)) return NULL;
    if (!self->xsdir) {
        PyErr_SetString(PyExc_RuntimeError, "XsDir not initialized");
        return NULL;
    }

    const alea_nuc_xsdir_entry_t* entry = alea_nuc_xsdir_find(self->xsdir, zaid);
    return PyAleaXsDir_entry_dict(entry);
}

static PyObject* PyAleaXsDir_find_temperature(
        PyAleaXsDirObject* self, PyObject* args, PyObject* kwds) {
    static char* kwlist[] = {"zaid", "kT", "abs_tolerance", NULL};
    const char* zaid;
    double kT, abs_tolerance = 0.0;
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "sd|d", kwlist,
                                     &zaid, &kT, &abs_tolerance))
        return NULL;
    if (!self->xsdir) {
        PyErr_SetString(PyExc_RuntimeError, "XsDir not initialized");
        return NULL;
    }
    const alea_nuc_xsdir_entry_t* entry = NULL;
    alea_error_t err = alea_nuc_xsdir_find_temperature(
        self->xsdir, zaid, kT, abs_tolerance, &entry);
    if (err == ALEA_ERR_NOT_FOUND) Py_RETURN_NONE;
    if (err != ALEA_OK) {
        PyErr_Format(PyExc_ValueError,
                     "temperature lookup failed for '%s': %s",
                     zaid, alea_error_string(err));
        return NULL;
    }
    return PyAleaXsDir_entry_dict(entry);
}

static PyObject* PyAleaXsDir_find_temperature_bracket(
        PyAleaXsDirObject* self, PyObject* args) {
    const char* zaid;
    double kT;
    if (!PyArg_ParseTuple(args, "sd", &zaid, &kT)) return NULL;
    if (!self->xsdir) {
        PyErr_SetString(PyExc_RuntimeError, "XsDir not initialized");
        return NULL;
    }
    const alea_nuc_xsdir_entry_t *lower = NULL, *upper = NULL;
    double upper_fraction = 0.0;
    alea_error_t err = alea_nuc_xsdir_find_temperature_bracket(
        self->xsdir, zaid, kT, &lower, &upper, &upper_fraction);
    if (err == ALEA_ERR_NOT_FOUND) Py_RETURN_NONE;
    if (err != ALEA_OK) {
        PyErr_Format(PyExc_ValueError,
                     "temperature bracket lookup failed for '%s': %s",
                     zaid, alea_error_string(err));
        return NULL;
    }
    PyObject* lower_dict = PyAleaXsDir_entry_dict(lower);
    PyObject* upper_dict = PyAleaXsDir_entry_dict(upper);
    if (!lower_dict || !upper_dict) {
        Py_XDECREF(lower_dict);
        Py_XDECREF(upper_dict);
        return NULL;
    }
    return Py_BuildValue("{s:N,s:N,s:d}",
        "lower", lower_dict, "upper", upper_dict,
        "upper_fraction", upper_fraction);
}

static PyObject* PyAleaXsDir_get_count(PyAleaXsDirObject* self, void* closure) {
    (void)closure;
    if (!self->xsdir) {
        PyErr_SetString(PyExc_RuntimeError, "XsDir not initialized");
        return NULL;
    }
    return PyLong_FromSize_t(alea_nuc_xsdir_count(self->xsdir));
}

static PyGetSetDef PyAleaXsDir_getsetters[] = {
    {"count", (getter)PyAleaXsDir_get_count, NULL, "Number of entries", NULL},
    {NULL}
};

static PyMethodDef PyAleaXsDir_methods[] = {
    {"find", (PyCFunction)PyAleaXsDir_find, METH_VARARGS,
     "find(zaid) -> dict or None\n\nFind xsdir entry by ZAID string (e.g. '92235.80c')."},
    {"find_temperature", (PyCFunction)PyAleaXsDir_find_temperature,
     METH_VARARGS | METH_KEYWORDS,
     "find_temperature(zaid, kT, abs_tolerance=0.0) -> dict or None\n\n"
     "Find an unambiguous table in the same ZAID family at kT (MeV)."},
    {"find_temperature_bracket",
     (PyCFunction)PyAleaXsDir_find_temperature_bracket, METH_VARARGS,
     "find_temperature_bracket(zaid, kT) -> dict or None\n\n"
     "Find bounded temperature tables and the upper interpolation fraction."},
    {NULL}
};

static PyTypeObject PyAleaXsDirType = {
    PyVarObject_HEAD_INIT(NULL, 0)
    .tp_name = "pyalea._alea.XsDir",
    .tp_doc = PyDoc_STR("Nuclear data cross-section directory (xsdir/xsdata)."),
    .tp_basicsize = sizeof(PyAleaXsDirObject),
    .tp_itemsize = 0,
    .tp_flags = Py_TPFLAGS_DEFAULT,
    .tp_new = PyAleaXsDir_new,
    .tp_init = (initproc)PyAleaXsDir_init,
    .tp_dealloc = (destructor)PyAleaXsDir_dealloc,
    .tp_methods = PyAleaXsDir_methods,
    .tp_getset = PyAleaXsDir_getsetters,
};

/* ============================================================================
 * Nuclide Python Type
 * ============================================================================ */

typedef struct {
    PyObject_HEAD
    alea_nuc_nuclide_t* nuc;
    int owned;                    /* 1 = we own it (must free), 0 = borrowed */
    PyObject* xsdir_ref;         /* Keep xsdir alive if borrowed */
    int default_photon_audit_cached;
    alea_error_t default_photon_audit_error;
    alea_nuc_photon_production_audit_t default_photon_audit;
} PyAleaNuclideObject;

static void PyAleaNuclide_dealloc(PyAleaNuclideObject* self) {
    if (self->nuc && self->owned) {
        alea_nuc_nuclide_free(self->nuc);
    }
    Py_XDECREF(self->xsdir_ref);
    Py_TYPE(self)->tp_free((PyObject*)self);
}

static PyObject* PyAleaNuclide_new(PyTypeObject* type, PyObject* args, PyObject* kwds) {
    (void)args; (void)kwds;
    PyAleaNuclideObject* self = (PyAleaNuclideObject*)type->tp_alloc(type, 0);
    if (self) {
        self->nuc = NULL;
        self->owned = 0;
        self->xsdir_ref = NULL;
        self->default_photon_audit_cached = 0;
        self->default_photon_audit_error = ALEA_OK;
        self->default_photon_audit =
            (alea_nuc_photon_production_audit_t){0};
    }
    return (PyObject*)self;
}

static int PyAleaNuclide_init(PyAleaNuclideObject* self, PyObject* args, PyObject* kwds) {
    static char* kwlist[] = {"xsdir", "zaid", "cached", NULL};
    PyAleaXsDirObject* xsdir_obj;
    const char* zaid;
    int cached = 1;

    if (!PyArg_ParseTupleAndKeywords(args, kwds, "O!s|p", kwlist,
                                     &PyAleaXsDirType, &xsdir_obj, &zaid, &cached))
        return -1;

    if (!xsdir_obj->xsdir) {
        PyErr_SetString(PyExc_RuntimeError, "XsDir not initialized");
        return -1;
    }

    alea_nuc_nuclide_t* nuc;
    int owned;
    PyObject* xsdir_ref = NULL;
    if (cached) {
        /* Borrowed from xsdir cache */
        nuc = alea_nuc_xsdir_get_nuclide(xsdir_obj->xsdir, zaid);
        if (!nuc) {
            PyErr_Format(PyExc_ValueError, "Failed to load nuclide '%s'", zaid);
            return -1;
        }
        owned = 0;
        xsdir_ref = (PyObject*)xsdir_obj;
        Py_INCREF(xsdir_ref);
    } else {
        /* Owned copy */
        nuc = alea_nuc_load_nuclide(xsdir_obj->xsdir, zaid);
        if (!nuc) {
            PyErr_Format(PyExc_ValueError, "Failed to load nuclide '%s'", zaid);
            return -1;
        }
        owned = 1;
    }

    if (self->nuc && self->owned) alea_nuc_nuclide_free(self->nuc);
    Py_XDECREF(self->xsdir_ref);
    self->nuc = nuc;
    self->owned = owned;
    self->xsdir_ref = xsdir_ref;
    self->default_photon_audit_cached = 0;
    self->default_photon_audit_error = ALEA_OK;
    self->default_photon_audit = (alea_nuc_photon_production_audit_t){0};
    return 0;
}

static alea_error_t PyAleaNuclide_default_photon_audit(
        PyAleaNuclideObject* self,
        alea_nuc_photon_production_audit_t* report) {
    if (self->xsdir_ref) {
        PyAleaXsDirObject* xsdir =
            (PyAleaXsDirObject*)self->xsdir_ref;
        PyAleaPhotonAuditCacheEntry* entry = xsdir->photon_audits;
        while (entry && entry->nuc != self->nuc) entry = entry->next;
        if (!entry) {
            entry = (PyAleaPhotonAuditCacheEntry*)calloc(1, sizeof(*entry));
            if (!entry) return ALEA_ERR_OUT_OF_MEMORY;
            entry->nuc = self->nuc;
            entry->error = alea_nuc_photon_production_audit(
                self->nuc, 5e-4, 1e-14, &entry->report);
            entry->next = xsdir->photon_audits;
            xsdir->photon_audits = entry;
        }
        if (report) *report = entry->report;
        return entry->error;
    }
    if (!self->default_photon_audit_cached) {
        self->default_photon_audit_error =
            alea_nuc_photon_production_audit(
                self->nuc, 5e-4, 1e-14, &self->default_photon_audit);
        self->default_photon_audit_cached = 1;
    }
    if (report) *report = self->default_photon_audit;
    return self->default_photon_audit_error;
}

/* Cross-section lookups */
static PyObject* PyAleaNuclide_xs_total(PyAleaNuclideObject* self, PyObject* args) {
    double energy;
    if (!PyArg_ParseTuple(args, "d", &energy)) return NULL;
    if (!self->nuc) { PyErr_SetString(PyExc_RuntimeError, "Nuclide not initialized"); return NULL; }
    return PyFloat_FromDouble(alea_nuc_xs_total(self->nuc, energy));
}

static PyObject* PyAleaNuclide_xs_absorption(PyAleaNuclideObject* self, PyObject* args) {
    double energy;
    if (!PyArg_ParseTuple(args, "d", &energy)) return NULL;
    if (!self->nuc) { PyErr_SetString(PyExc_RuntimeError, "Nuclide not initialized"); return NULL; }
    return PyFloat_FromDouble(alea_nuc_xs_absorption(self->nuc, energy));
}

static PyObject* PyAleaNuclide_xs_elastic(PyAleaNuclideObject* self, PyObject* args) {
    double energy;
    if (!PyArg_ParseTuple(args, "d", &energy)) return NULL;
    if (!self->nuc) { PyErr_SetString(PyExc_RuntimeError, "Nuclide not initialized"); return NULL; }
    return PyFloat_FromDouble(alea_nuc_xs_elastic(self->nuc, energy));
}

static PyObject* PyAleaNuclide_xs_reaction(PyAleaNuclideObject* self, PyObject* args) {
    int mt;
    double energy;
    if (!PyArg_ParseTuple(args, "id", &mt, &energy)) return NULL;
    if (!self->nuc) { PyErr_SetString(PyExc_RuntimeError, "Nuclide not initialized"); return NULL; }
    return PyFloat_FromDouble(alea_nuc_xs_reaction(self->nuc, mt, energy));
}

static PyObject* PyAleaNuclide_xs_heating(PyAleaNuclideObject* self, PyObject* args) {
    double energy;
    if (!PyArg_ParseTuple(args, "d", &energy)) return NULL;
    if (!self->nuc) { PyErr_SetString(PyExc_RuntimeError, "Nuclide not initialized"); return NULL; }
    return PyFloat_FromDouble(alea_nuc_xs_heating(self->nuc, energy));
}

static PyObject* PyAleaNuclide_xs_values(PyAleaNuclideObject* self, PyObject* args) {
    int kind, mt;
    PyObject* energies_obj;
    if (!PyArg_ParseTuple(args, "iiO", &kind, &mt, &energies_obj)) return NULL;
    if (!self->nuc) { PyErr_SetString(PyExc_RuntimeError, "Nuclide not initialized"); return NULL; }
    if (kind < 0 || kind > 5) { PyErr_SetString(PyExc_ValueError, "invalid cross-section kind"); return NULL; }
    PyArrayObject* energies = (PyArrayObject*)PyArray_FROM_OTF(energies_obj, NPY_DOUBLE, NPY_ARRAY_CARRAY_RO);
    if (!energies) return NULL;
    npy_intp size = PyArray_SIZE(energies), dims[1] = {size};
    PyArrayObject* values = (PyArrayObject*)PyArray_SimpleNew(1, dims, NPY_DOUBLE);
    if (!values) { Py_DECREF(energies); return NULL; }
    const double* input = PyArray_DATA(energies);
    double* output = PyArray_DATA(values);
    Py_BEGIN_ALLOW_THREADS
    for (npy_intp i = 0; i < size; i++) {
        switch (kind) {
            case 0: output[i] = alea_nuc_xs_total(self->nuc, input[i]); break;
            case 1: output[i] = alea_nuc_xs_absorption(self->nuc, input[i]); break;
            case 2: output[i] = alea_nuc_xs_elastic(self->nuc, input[i]); break;
            case 3: output[i] = alea_nuc_xs_reaction(self->nuc, mt, input[i]); break;
            case 4: output[i] = alea_nuc_xs_heating(self->nuc, input[i]); break;
            default: output[i] = alea_nuc_xs_photon_production_total(self->nuc, input[i]); break;
        }
    }
    Py_END_ALLOW_THREADS
    Py_DECREF(energies);
    return (PyObject*)values;
}

static PyObject* PyAleaNuclide_heating_per_collision(PyAleaNuclideObject* self, PyObject* args) {
    double energy;
    if (!PyArg_ParseTuple(args, "d", &energy)) return NULL;
    if (!self->nuc) { PyErr_SetString(PyExc_RuntimeError, "Nuclide not initialized"); return NULL; }
    return PyFloat_FromDouble(alea_nuc_heating_per_collision(self->nuc, energy));
}

static PyObject* PyAleaNuclide_nu_bar(PyAleaNuclideObject* self, PyObject* args) {
    double energy;
    if (!PyArg_ParseTuple(args, "d", &energy)) return NULL;
    if (!self->nuc) { PyErr_SetString(PyExc_RuntimeError, "Nuclide not initialized"); return NULL; }
    if (!self->nuc->fission) {
        PyErr_SetString(PyExc_ValueError, "Nuclide is not fissile");
        return NULL;
    }
    return PyFloat_FromDouble(alea_nuc_nu_bar(self->nuc, energy));
}

static PyObject* PyAleaNuclide_prompt_nu_bar(PyAleaNuclideObject* self, PyObject* args) {
    double energy;
    if (!PyArg_ParseTuple(args, "d", &energy)) return NULL;
    if (!self->nuc) { PyErr_SetString(PyExc_RuntimeError, "Nuclide not initialized"); return NULL; }
    if (!self->nuc->fission) {
        PyErr_SetString(PyExc_ValueError, "Nuclide is not fissile");
        return NULL;
    }
    return PyFloat_FromDouble(alea_nuc_prompt_nu_bar(self->nuc, energy));
}

static PyObject* PyAleaNuclide_delayed_nu_bar(PyAleaNuclideObject* self, PyObject* args) {
    double energy;
    if (!PyArg_ParseTuple(args, "d", &energy)) return NULL;
    if (!self->nuc) { PyErr_SetString(PyExc_RuntimeError, "Nuclide not initialized"); return NULL; }
    if (!self->nuc->fission) {
        PyErr_SetString(PyExc_ValueError, "Nuclide is not fissile");
        return NULL;
    }
    return PyFloat_FromDouble(alea_nuc_delayed_nu_bar(self->nuc, energy));
}

static PyObject* PyAleaNuclide_capabilities(
        PyAleaNuclideObject* self, PyObject* Py_UNUSED(ignored)) {
    if (!self->nuc) {
        PyErr_SetString(PyExc_RuntimeError, "Nuclide not initialized");
        return NULL;
    }
    alea_nuc_capability_report_t report;
    alea_error_t err = alea_nuc_capabilities(self->nuc, &report);
    if (err != ALEA_OK && err != ALEA_ERR_UNSUPPORTED) {
        PyErr_Format(PyExc_RuntimeError, "Capability inspection failed: %s",
                     alea_error_string(err));
        return NULL;
    }
    return Py_BuildValue("{s:k,s:k,s:i,s:i,s:i,s:i,s:s}",
        "available_mask", (unsigned long)report.available_capabilities,
        "missing_mask", (unsigned long)report.missing_capabilities,
        "issue", (int)report.issue,
        "component_index", report.component_index,
        "mt", report.mt,
        "law", report.law,
        "detail", report.detail);
}

static PyObject* PyAleaNuclide_doppler_broaden(PyAleaNuclideObject* self, PyObject* args) {
    double kT_target;
    if (!PyArg_ParseTuple(args, "d", &kT_target)) return NULL;
    if (!self->nuc) { PyErr_SetString(PyExc_RuntimeError, "Nuclide not initialized"); return NULL; }
    if (!self->owned) {
        PyErr_SetString(PyExc_RuntimeError, "Cannot Doppler-broaden a cached (borrowed) nuclide. Load with cached=False.");
        return NULL;
    }
    alea_error_t err = alea_nuc_doppler_broaden(self->nuc, kT_target);
    if (err != ALEA_OK) {
        PyErr_Format(PyExc_ValueError, "Doppler broadening failed: %s", alea_error_string(err));
        return NULL;
    }
    self->default_photon_audit_cached = 0;
    Py_RETURN_NONE;
}

/* Photon cross sections */
static PyObject* PyAleaNuclide_photon_xs_incoherent(PyAleaNuclideObject* self, PyObject* args) {
    double energy;
    if (!PyArg_ParseTuple(args, "d", &energy)) return NULL;
    if (!self->nuc) { PyErr_SetString(PyExc_RuntimeError, "Nuclide not initialized"); return NULL; }
    return PyFloat_FromDouble(alea_nuc_photon_xs_incoherent(self->nuc, energy));
}

static PyObject* PyAleaNuclide_photon_xs_coherent(PyAleaNuclideObject* self, PyObject* args) {
    double energy;
    if (!PyArg_ParseTuple(args, "d", &energy)) return NULL;
    if (!self->nuc) { PyErr_SetString(PyExc_RuntimeError, "Nuclide not initialized"); return NULL; }
    return PyFloat_FromDouble(alea_nuc_photon_xs_coherent(self->nuc, energy));
}

static PyObject* PyAleaNuclide_photon_xs_photoelectric(PyAleaNuclideObject* self, PyObject* args) {
    double energy;
    if (!PyArg_ParseTuple(args, "d", &energy)) return NULL;
    if (!self->nuc) { PyErr_SetString(PyExc_RuntimeError, "Nuclide not initialized"); return NULL; }
    return PyFloat_FromDouble(alea_nuc_photon_xs_photoelectric(self->nuc, energy));
}

static PyObject* PyAleaNuclide_photon_xs_pair(PyAleaNuclideObject* self, PyObject* args) {
    double energy;
    if (!PyArg_ParseTuple(args, "d", &energy)) return NULL;
    if (!self->nuc) { PyErr_SetString(PyExc_RuntimeError, "Nuclide not initialized"); return NULL; }
    return PyFloat_FromDouble(alea_nuc_photon_xs_pair(self->nuc, energy));
}

static PyObject* PyAleaNuclide_photon_production_xs(
        PyAleaNuclideObject* self, PyObject* args) {
    double energy;
    if (!PyArg_ParseTuple(args, "d", &energy)) return NULL;
    if (!self->nuc) {
        PyErr_SetString(PyExc_RuntimeError, "Nuclide not initialized");
        return NULL;
    }
    return PyFloat_FromDouble(
        alea_nuc_xs_photon_production_total(self->nuc, energy));
}

static PyObject* PyAleaNuclide_photon_productions(
        PyAleaNuclideObject* self, PyObject* Py_UNUSED(ignored)) {
    if (!self->nuc) {
        PyErr_SetString(PyExc_RuntimeError, "Nuclide not initialized");
        return NULL;
    }
    PyObject* result = PyList_New(self->nuc->n_photon_productions);
    if (!result) return NULL;
    for (int i = 0; i < self->nuc->n_photon_productions; i++) {
        const alea_nuc_photon_production_t* production =
            &self->nuc->photon_productions[i];
        PyObject* item = Py_BuildValue("{s:i,s:i,s:i,s:i,s:O}",
            "index", i,
            "mt", production->mt,
            "parent_mt", production->parent_mt,
            "mf", production->mf,
            "production_xs", production->production_xs ? Py_True : Py_False);
        if (!item) {
            Py_DECREF(result);
            return NULL;
        }
        if (production->spectrum) {
            PyObject* law = PyLong_FromLong(production->spectrum->law);
            if (!law || PyDict_SetItemString(item, "energy_law", law) < 0) {
                Py_XDECREF(law);
                Py_DECREF(item);
                Py_DECREF(result);
                return NULL;
            }
            Py_DECREF(law);
        }
        PyList_SET_ITEM(result, i, item);
    }
    return result;
}

static PyObject* PyAleaNuclide_photon_production_audit(
        PyAleaNuclideObject* self, PyObject* args, PyObject* kwds) {
    static char* kwlist[] = {
        "relative_tolerance", "absolute_tolerance", NULL
    };
    double relative_tolerance = 5e-4;
    double absolute_tolerance = 1e-14;
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "|dd", kwlist,
            &relative_tolerance, &absolute_tolerance))
        return NULL;
    if (!self->nuc) {
        PyErr_SetString(PyExc_RuntimeError, "Nuclide not initialized");
        return NULL;
    }

    alea_nuc_photon_production_audit_t report;
    alea_error_t err =
        relative_tolerance == 5e-4 && absolute_tolerance == 1e-14
        ? PyAleaNuclide_default_photon_audit(self, &report)
        : alea_nuc_photon_production_audit(
            self->nuc, relative_tolerance, absolute_tolerance, &report);
    if (err != ALEA_OK) {
        PyErr_Format(PyExc_ValueError,
                     "Photon-production audit failed: %s",
                     alea_error_string(err));
        return NULL;
    }
    return Py_BuildValue("{s:O,s:O,s:d,s:d,s:d,s:i}",
        "aggregate_available",
        report.aggregate_available ? Py_True : Py_False,
        "native_grid_consistent",
        report.native_grid_consistent ? Py_True : Py_False,
        "maximum_absolute_difference",
        report.maximum_absolute_difference,
        "maximum_relative_difference",
        report.maximum_relative_difference,
        "worst_energy", report.worst_energy,
        "worst_energy_index", report.worst_energy_index);
}

/* URR factors */
static PyObject* PyAleaNuclide_urr_factors(PyAleaNuclideObject* self, PyObject* args) {
    double energy, xi;
    if (!PyArg_ParseTuple(args, "dd", &energy, &xi)) return NULL;
    if (!self->nuc) { PyErr_SetString(PyExc_RuntimeError, "Nuclide not initialized"); return NULL; }

    double factors[5];
    int applies = alea_nuc_urr_factors(self->nuc, energy, xi, factors);
    if (!applies) Py_RETURN_NONE;

    return Py_BuildValue("{s:d, s:d, s:d, s:d, s:d}",
        "total", factors[0], "elastic", factors[1],
        "fission", factors[2], "capture", factors[3],
        "heating", factors[4]);
}

/* Energy grid access */
static PyObject* PyAleaNuclide_get_energy_grid(PyAleaNuclideObject* self, PyObject* args) {
    (void)args;
    if (!self->nuc) { PyErr_SetString(PyExc_RuntimeError, "Nuclide not initialized"); return NULL; }

    PyObject* list = PyList_New(self->nuc->n_energies);
    if (!list) return NULL;
    for (int i = 0; i < self->nuc->n_energies; i++) {
        PyList_SET_ITEM(list, i, PyFloat_FromDouble(self->nuc->energy[i]));
    }
    return list;
}

/* Reactions list */
static PyObject* PyAleaNuclide_get_reactions(PyAleaNuclideObject* self, PyObject* args) {
    (void)args;
    if (!self->nuc) { PyErr_SetString(PyExc_RuntimeError, "Nuclide not initialized"); return NULL; }

    PyObject* list = PyList_New(self->nuc->n_reactions);
    if (!list) return NULL;
    for (int i = 0; i < self->nuc->n_reactions; i++) {
        const alea_nuc_reaction_t* rxn = &self->nuc->reactions[i];
        PyObject* d = Py_BuildValue("{s:i, s:d, s:i}",
            "mt", rxn->mt,
            "q_value", rxn->q_value,
            "ty", rxn->ty);
        if (!d) { Py_DECREF(list); return NULL; }
        PyList_SET_ITEM(list, i, d);
    }
    return list;
}

/* Reaction yield and classification */
static PyObject* PyAleaNuclide_reaction_yield(PyAleaNuclideObject* self, PyObject* args) {
    int mt;
    double energy;
    if (!PyArg_ParseTuple(args, "id", &mt, &energy)) return NULL;
    if (!self->nuc) { PyErr_SetString(PyExc_RuntimeError, "Nuclide not initialized"); return NULL; }
    return PyFloat_FromDouble(alea_nuc_reaction_yield(self->nuc, mt, energy));
}

/* Properties */
static PyObject* PyAleaNuclide_get_zaid(PyAleaNuclideObject* self, void* closure) {
    (void)closure;
    if (!self->nuc) { PyErr_SetString(PyExc_RuntimeError, "Nuclide not initialized"); return NULL; }
    return PyUnicode_FromString(self->nuc->zaid);
}

static PyObject* PyAleaNuclide_get_Z(PyAleaNuclideObject* self, void* closure) {
    (void)closure;
    if (!self->nuc) { PyErr_SetString(PyExc_RuntimeError, "Nuclide not initialized"); return NULL; }
    return PyLong_FromLong(self->nuc->Z);
}

static PyObject* PyAleaNuclide_get_A(PyAleaNuclideObject* self, void* closure) {
    (void)closure;
    if (!self->nuc) { PyErr_SetString(PyExc_RuntimeError, "Nuclide not initialized"); return NULL; }
    return PyLong_FromLong(self->nuc->A);
}

static PyObject* PyAleaNuclide_get_awr(PyAleaNuclideObject* self, void* closure) {
    (void)closure;
    if (!self->nuc) { PyErr_SetString(PyExc_RuntimeError, "Nuclide not initialized"); return NULL; }
    return PyFloat_FromDouble(self->nuc->awr);
}

static PyObject* PyAleaNuclide_get_temperature(PyAleaNuclideObject* self, void* closure) {
    (void)closure;
    if (!self->nuc) { PyErr_SetString(PyExc_RuntimeError, "Nuclide not initialized"); return NULL; }
    return PyFloat_FromDouble(self->nuc->temperature);
}

static PyObject* PyAleaNuclide_get_is_fissile(PyAleaNuclideObject* self, void* closure) {
    (void)closure;
    if (!self->nuc) { PyErr_SetString(PyExc_RuntimeError, "Nuclide not initialized"); return NULL; }
    return PyBool_FromLong(self->nuc->fission != NULL);
}

static PyObject* PyAleaNuclide_get_has_urr(PyAleaNuclideObject* self, void* closure) {
    (void)closure;
    if (!self->nuc) { PyErr_SetString(PyExc_RuntimeError, "Nuclide not initialized"); return NULL; }
    return PyBool_FromLong(self->nuc->urr != NULL);
}

static PyObject* PyAleaNuclide_get_has_photon(PyAleaNuclideObject* self, void* closure) {
    (void)closure;
    if (!self->nuc) { PyErr_SetString(PyExc_RuntimeError, "Nuclide not initialized"); return NULL; }
    return PyBool_FromLong(self->nuc->photon != NULL);
}

static PyObject* PyAleaNuclide_get_n_reactions(PyAleaNuclideObject* self, void* closure) {
    (void)closure;
    if (!self->nuc) { PyErr_SetString(PyExc_RuntimeError, "Nuclide not initialized"); return NULL; }
    return PyLong_FromLong(self->nuc->n_reactions);
}

static PyObject* PyAleaNuclide_get_n_energies(PyAleaNuclideObject* self, void* closure) {
    (void)closure;
    if (!self->nuc) { PyErr_SetString(PyExc_RuntimeError, "Nuclide not initialized"); return NULL; }
    return PyLong_FromLong(self->nuc->n_energies);
}

static PyGetSetDef PyAleaNuclide_getsetters[] = {
    {"zaid", (getter)PyAleaNuclide_get_zaid, NULL, "ZAID string (e.g. '92235.80c')", NULL},
    {"Z", (getter)PyAleaNuclide_get_Z, NULL, "Atomic number", NULL},
    {"A", (getter)PyAleaNuclide_get_A, NULL, "Mass number", NULL},
    {"awr", (getter)PyAleaNuclide_get_awr, NULL, "Atomic weight ratio", NULL},
    {"temperature", (getter)PyAleaNuclide_get_temperature, NULL, "Temperature kT in MeV", NULL},
    {"is_fissile", (getter)PyAleaNuclide_get_is_fissile, NULL, "True if nuclide has fission data", NULL},
    {"has_urr", (getter)PyAleaNuclide_get_has_urr, NULL, "True if nuclide has URR probability tables", NULL},
    {"has_photon", (getter)PyAleaNuclide_get_has_photon, NULL, "True if nuclide has photon data", NULL},
    {"n_reactions", (getter)PyAleaNuclide_get_n_reactions, NULL, "Number of non-elastic reactions", NULL},
    {"n_energies", (getter)PyAleaNuclide_get_n_energies, NULL, "Number of energy grid points", NULL},
    {NULL}
};

static PyMethodDef PyAleaNuclide_methods[] = {
    {"xs_total", (PyCFunction)PyAleaNuclide_xs_total, METH_VARARGS,
     "xs_total(energy) -> float\n\nTotal cross section in barns at energy (MeV)."},
    {"xs_absorption", (PyCFunction)PyAleaNuclide_xs_absorption, METH_VARARGS,
     "xs_absorption(energy) -> float\n\nAbsorption cross section in barns."},
    {"xs_elastic", (PyCFunction)PyAleaNuclide_xs_elastic, METH_VARARGS,
     "xs_elastic(energy) -> float\n\nElastic scattering cross section in barns."},
    {"xs_reaction", (PyCFunction)PyAleaNuclide_xs_reaction, METH_VARARGS,
     "xs_reaction(mt, energy) -> float\n\nCross section for reaction MT at energy (MeV)."},
    {"xs_heating", (PyCFunction)PyAleaNuclide_xs_heating, METH_VARARGS,
     "xs_heating(energy) -> float\n\nHeating number (MeV·barn) at energy."},
    {"xs_values", (PyCFunction)PyAleaNuclide_xs_values, METH_VARARGS,
     "xs_values(kind, mt, energies) -> numpy.ndarray\n\nBulk cross-section evaluation."},
    {"heating_per_collision", (PyCFunction)PyAleaNuclide_heating_per_collision, METH_VARARGS,
     "heating_per_collision(energy) -> float\n\nAverage energy deposited per collision (MeV)."},
    {"nu_bar", (PyCFunction)PyAleaNuclide_nu_bar, METH_VARARGS,
     "nu_bar(energy) -> float\n\nAverage neutrons per fission at energy (MeV)."},
    {"prompt_nu_bar", (PyCFunction)PyAleaNuclide_prompt_nu_bar, METH_VARARGS,
     "prompt_nu_bar(energy) -> float\n\nAverage prompt neutrons per fission."},
    {"delayed_nu_bar", (PyCFunction)PyAleaNuclide_delayed_nu_bar, METH_VARARGS,
     "delayed_nu_bar(energy) -> float\n\nAverage delayed neutrons per fission."},
    {"capabilities", (PyCFunction)PyAleaNuclide_capabilities, METH_NOARGS,
     "capabilities() -> dict\n\nInspect decoded collision-physics capabilities."},
    {"doppler_broaden", (PyCFunction)PyAleaNuclide_doppler_broaden, METH_VARARGS,
     "doppler_broaden(kT_target) -> None\n\n"
     "Doppler-broaden cross sections in-place to temperature kT_target (MeV).\n"
     "Only works on owned (non-cached) nuclides."},
    {"urr_factors", (PyCFunction)PyAleaNuclide_urr_factors, METH_VARARGS,
     "urr_factors(energy, xi) -> dict or None\n\n"
     "Get URR probability table factors. Returns None if URR doesn't apply.\n"
     "Returns dict with keys: total, elastic, fission, capture, heating."},
    {"photon_xs_incoherent", (PyCFunction)PyAleaNuclide_photon_xs_incoherent, METH_VARARGS,
     "photon_xs_incoherent(energy) -> float\n\nCompton (incoherent) scattering cross section (barns)."},
    {"photon_xs_coherent", (PyCFunction)PyAleaNuclide_photon_xs_coherent, METH_VARARGS,
     "photon_xs_coherent(energy) -> float\n\nRayleigh (coherent) scattering cross section (barns)."},
    {"photon_xs_photoelectric", (PyCFunction)PyAleaNuclide_photon_xs_photoelectric, METH_VARARGS,
     "photon_xs_photoelectric(energy) -> float\n\nPhotoelectric cross section (barns)."},
    {"photon_xs_pair", (PyCFunction)PyAleaNuclide_photon_xs_pair, METH_VARARGS,
     "photon_xs_pair(energy) -> float\n\nPair production cross section (barns)."},
    {"photon_production_xs", (PyCFunction)PyAleaNuclide_photon_production_xs, METH_VARARGS,
     "photon_production_xs(energy) -> float\n\nAggregate neutron-induced photon-production cross section (barns)."},
    {"photon_productions", (PyCFunction)PyAleaNuclide_photon_productions, METH_NOARGS,
     "photon_productions() -> list of dict\n\nDecoded neutron-induced photon-production channel metadata."},
    {"photon_production_audit",
     (PyCFunction)PyAleaNuclide_photon_production_audit,
     METH_VARARGS | METH_KEYWORDS,
     "photon_production_audit(relative_tolerance=5e-4, absolute_tolerance=1e-14) -> dict\n\n"
     "Compare the GPD aggregate with decoded channels on the native energy grid."},
    {"energy_grid", (PyCFunction)PyAleaNuclide_get_energy_grid, METH_NOARGS,
     "energy_grid() -> list of float\n\nGet the energy grid (MeV)."},
    {"reactions", (PyCFunction)PyAleaNuclide_get_reactions, METH_NOARGS,
     "reactions() -> list of dict\n\nGet list of reactions with mt, q_value, ty."},
    {"reaction_yield", (PyCFunction)PyAleaNuclide_reaction_yield, METH_VARARGS,
     "reaction_yield(mt, energy) -> float\n\nGet neutron yield for reaction MT at energy."},
    {NULL}
};

static PyTypeObject PyAleaNuclideType = {
    PyVarObject_HEAD_INIT(NULL, 0)
    .tp_name = "pyalea._alea.Nuclide",
    .tp_doc = PyDoc_STR("Decoded ACE nuclide with cross-section lookup."),
    .tp_basicsize = sizeof(PyAleaNuclideObject),
    .tp_itemsize = 0,
    .tp_flags = Py_TPFLAGS_DEFAULT,
    .tp_new = PyAleaNuclide_new,
    .tp_init = (initproc)PyAleaNuclide_init,
    .tp_dealloc = (destructor)PyAleaNuclide_dealloc,
    .tp_methods = PyAleaNuclide_methods,
    .tp_getset = PyAleaNuclide_getsetters,
};

/* ============================================================================
 * Thermal Scattering Python Type
 * ============================================================================ */

typedef struct {
    PyObject_HEAD
    alea_nuc_thermal_t* thermal;
    PyObject* xsdir_ref;
} PyAleaThermalObject;

static void PyAleaThermal_dealloc(PyAleaThermalObject* self) {
    if (self->thermal) alea_nuc_thermal_free(self->thermal);
    Py_XDECREF(self->xsdir_ref);
    Py_TYPE(self)->tp_free((PyObject*)self);
}

static PyObject* PyAleaThermal_new(PyTypeObject* type, PyObject* args, PyObject* kwds) {
    (void)args; (void)kwds;
    PyAleaThermalObject* self = (PyAleaThermalObject*)type->tp_alloc(type, 0);
    if (self) { self->thermal = NULL; self->xsdir_ref = NULL; }
    return (PyObject*)self;
}

static int PyAleaThermal_init(PyAleaThermalObject* self, PyObject* args, PyObject* kwds) {
    static char* kwlist[] = {"xsdir", "zaid", NULL};
    PyAleaXsDirObject* xsdir_obj;
    const char* zaid;
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "O!s", kwlist,
                                     &PyAleaXsDirType, &xsdir_obj, &zaid))
        return -1;
    if (!xsdir_obj->xsdir) {
        PyErr_SetString(PyExc_RuntimeError, "XsDir not initialized");
        return -1;
    }
    alea_nuc_thermal_t* thermal = alea_nuc_load_thermal(xsdir_obj->xsdir, zaid);
    if (!thermal) {
        PyErr_Format(PyExc_ValueError, "Failed to load thermal scattering table '%s'", zaid);
        return -1;
    }
    PyObject* xsdir_ref = (PyObject*)xsdir_obj;
    Py_INCREF(xsdir_ref);
    if (self->thermal) alea_nuc_thermal_free(self->thermal);
    Py_XDECREF(self->xsdir_ref);
    self->thermal = thermal;
    self->xsdir_ref = xsdir_ref;
    return 0;
}

static PyObject* PyAleaThermal_xs_values(PyAleaThermalObject* self, PyObject* args) {
    int kind;
    PyObject* energies_obj;
    if (!PyArg_ParseTuple(args, "iO", &kind, &energies_obj)) return NULL;
    if (!self->thermal) { PyErr_SetString(PyExc_RuntimeError, "Thermal table not initialized"); return NULL; }
    if (kind < 0 || kind > 2) { PyErr_SetString(PyExc_ValueError, "invalid thermal cross-section kind"); return NULL; }
    PyArrayObject* energies = (PyArrayObject*)PyArray_FROM_OTF(energies_obj, NPY_DOUBLE, NPY_ARRAY_CARRAY_RO);
    if (!energies) return NULL;
    npy_intp size = PyArray_SIZE(energies), dims[1] = {size};
    PyArrayObject* values = (PyArrayObject*)PyArray_SimpleNew(1, dims, NPY_DOUBLE);
    if (!values) { Py_DECREF(energies); return NULL; }
    const double* input = PyArray_DATA(energies);
    double* output = PyArray_DATA(values);
    Py_BEGIN_ALLOW_THREADS
    for (npy_intp i = 0; i < size; i++) {
        if (kind == 0) output[i] = alea_nuc_thermal_xs_total(self->thermal, input[i]);
        else if (kind == 1) output[i] = alea_nuc_thermal_xs_elastic(self->thermal, input[i]);
        else output[i] = alea_nuc_thermal_xs_inelastic(self->thermal, input[i]);
    }
    Py_END_ALLOW_THREADS
    Py_DECREF(energies);
    return (PyObject*)values;
}

#define THERMAL_SCALAR_METHOD(name, function) \
static PyObject* PyAleaThermal_##name(PyAleaThermalObject* self, PyObject* args) { \
    double energy; \
    if (!PyArg_ParseTuple(args, "d", &energy)) return NULL; \
    if (!self->thermal) { PyErr_SetString(PyExc_RuntimeError, "Thermal table not initialized"); return NULL; } \
    return PyFloat_FromDouble(function(self->thermal, energy)); \
}

THERMAL_SCALAR_METHOD(xs_total, alea_nuc_thermal_xs_total)
THERMAL_SCALAR_METHOD(xs_elastic, alea_nuc_thermal_xs_elastic)
THERMAL_SCALAR_METHOD(xs_inelastic, alea_nuc_thermal_xs_inelastic)

static PyObject* PyAleaThermal_get_zaid(PyAleaThermalObject* self, void* closure) {
    (void)closure;
    return PyUnicode_FromString(self->thermal->zaid);
}
static PyObject* PyAleaThermal_get_awr(PyAleaThermalObject* self, void* closure) {
    (void)closure;
    return PyFloat_FromDouble(self->thermal->awr);
}
static PyObject* PyAleaThermal_get_temperature(PyAleaThermalObject* self, void* closure) {
    (void)closure;
    return PyFloat_FromDouble(self->thermal->temperature);
}

static PyMethodDef PyAleaThermal_methods[] = {
    {"xs_total", (PyCFunction)PyAleaThermal_xs_total, METH_VARARGS, "Total thermal scattering cross section (barns)."},
    {"xs_elastic", (PyCFunction)PyAleaThermal_xs_elastic, METH_VARARGS, "Elastic thermal scattering cross section (barns)."},
    {"xs_inelastic", (PyCFunction)PyAleaThermal_xs_inelastic, METH_VARARGS, "Inelastic thermal scattering cross section (barns)."},
    {"xs_values", (PyCFunction)PyAleaThermal_xs_values, METH_VARARGS, "Bulk thermal cross-section evaluation."},
    {NULL}
};

static PyGetSetDef PyAleaThermal_getsetters[] = {
    {"zaid", (getter)PyAleaThermal_get_zaid, NULL, "Thermal table ZAID", NULL},
    {"awr", (getter)PyAleaThermal_get_awr, NULL, "Atomic weight ratio", NULL},
    {"temperature", (getter)PyAleaThermal_get_temperature, NULL, "Temperature kT in MeV", NULL},
    {NULL}
};

static PyTypeObject PyAleaThermalType = {
    PyVarObject_HEAD_INIT(NULL, 0)
    .tp_name = "pyalea._alea.ThermalScattering",
    .tp_doc = PyDoc_STR("Decoded ACE thermal scattering table."),
    .tp_basicsize = sizeof(PyAleaThermalObject),
    .tp_flags = Py_TPFLAGS_DEFAULT,
    .tp_new = PyAleaThermal_new,
    .tp_init = (initproc)PyAleaThermal_init,
    .tp_dealloc = (destructor)PyAleaThermal_dealloc,
    .tp_methods = PyAleaThermal_methods,
    .tp_getset = PyAleaThermal_getsetters,
};

/* ============================================================================
 * NucMaterial Python Type
 * ============================================================================ */

typedef struct {
    PyObject_HEAD
    alea_nuc_material_t* mat;
    PyObject* refs;   /* List to keep nuclide refs alive */
} PyAleaNucMaterialObject;

static alea_error_t PyAleaNucMaterial_default_photon_audit(
        PyAleaNucMaterialObject* self,
        const alea_nuc_nuclide_t* nuc,
        alea_nuc_photon_production_audit_t* report) {
    Py_ssize_t count = self->refs ? PyList_GET_SIZE(self->refs) : 0;
    for (Py_ssize_t i = 0; i < count; i++) {
        PyAleaNuclideObject* item = (PyAleaNuclideObject*)
            PyList_GET_ITEM(self->refs, i);
        if (item->nuc == nuc)
            return PyAleaNuclide_default_photon_audit(item, report);
    }
    return alea_nuc_photon_production_audit(
        nuc, 5e-4, 1e-14, report);
}

static void PyAleaNucMaterial_dealloc(PyAleaNucMaterialObject* self) {
    if (self->mat) {
        alea_nuc_material_destroy(self->mat);
    }
    Py_XDECREF(self->refs);
    Py_TYPE(self)->tp_free((PyObject*)self);
}

static PyObject* PyAleaNucMaterial_new(PyTypeObject* type, PyObject* args, PyObject* kwds) {
    (void)args; (void)kwds;
    PyAleaNucMaterialObject* self = (PyAleaNucMaterialObject*)type->tp_alloc(type, 0);
    if (self) {
        self->mat = NULL;
        self->refs = NULL;
    }
    return (PyObject*)self;
}

static int PyAleaNucMaterial_init(PyAleaNucMaterialObject* self, PyObject* args, PyObject* kwds) {
    (void)args; (void)kwds;
    alea_nuc_material_t* mat = alea_nuc_material_create();
    if (!mat) {
        PyErr_SetString(PyExc_MemoryError, "Failed to create nuclear material");
        return -1;
    }
    PyObject* refs = PyList_New(0);
    if (!refs) {
        alea_nuc_material_destroy(mat);
        return -1;
    }
    if (self->mat) alea_nuc_material_destroy(self->mat);
    Py_XDECREF(self->refs);
    self->mat = mat;
    self->refs = refs;
    return 0;
}

static PyObject* PyAleaNucMaterial_add(PyAleaNucMaterialObject* self, PyObject* args) {
    PyAleaNuclideObject* nuc_obj;
    double number_density;
    if (!PyArg_ParseTuple(args, "O!d", &PyAleaNuclideType, &nuc_obj, &number_density))
        return NULL;

    if (!self->mat) {
        PyErr_SetString(PyExc_RuntimeError, "Material not initialized");
        return NULL;
    }
    if (!nuc_obj->nuc) {
        PyErr_SetString(PyExc_RuntimeError, "Nuclide not initialized");
        return NULL;
    }

    Py_ssize_t refs_before = PyList_GET_SIZE(self->refs);
    if (PyList_Append(self->refs, (PyObject*)nuc_obj) < 0)
        return NULL;

    alea_error_t err = alea_nuc_material_add(self->mat, nuc_obj->nuc, number_density);
    if (err != ALEA_OK) {
        PyList_SetSlice(self->refs, refs_before,
                        PyList_GET_SIZE(self->refs), NULL);
        PyErr_Format(PyExc_RuntimeError, "Failed to add nuclide: %s", alea_error_string(err));
        return NULL;
    }
    Py_RETURN_NONE;
}

static PyObject* PyAleaNucMaterial_add_temperature_mix(
        PyAleaNucMaterialObject* self, PyObject* args) {
    PyAleaNuclideObject *lower_obj, *upper_obj;
    double upper_fraction, number_density;
    if (!PyArg_ParseTuple(args, "O!O!dd",
            &PyAleaNuclideType, &lower_obj,
            &PyAleaNuclideType, &upper_obj,
            &upper_fraction, &number_density))
        return NULL;
    if (!self->mat || !lower_obj->nuc || !upper_obj->nuc) {
        PyErr_SetString(PyExc_RuntimeError,
                        "Material and nuclides must be initialized");
        return NULL;
    }

    /* Retain borrowed tables before mutating the native material. */
    Py_ssize_t refs_before = PyList_GET_SIZE(self->refs);
    if (PyList_Append(self->refs, (PyObject*)lower_obj) < 0)
        return NULL;
    if (PyList_Append(self->refs, (PyObject*)upper_obj) < 0) {
        PyList_SetSlice(self->refs, refs_before,
                        PyList_GET_SIZE(self->refs), NULL);
        return NULL;
    }

    alea_error_t err = alea_nuc_material_add_temperature_mix(
        self->mat, lower_obj->nuc, upper_obj->nuc,
        upper_fraction, number_density);
    if (err != ALEA_OK) {
        PyList_SetSlice(self->refs, refs_before,
                        PyList_GET_SIZE(self->refs), NULL);
        PyErr_Format(PyExc_ValueError,
                     "Failed to add temperature mixture: %s",
                     alea_error_string(err));
        return NULL;
    }
    Py_RETURN_NONE;
}

/* Macroscopic cross sections */
static PyObject* PyAleaNucMaterial_xs_total(PyAleaNucMaterialObject* self, PyObject* args) {
    double energy;
    if (!PyArg_ParseTuple(args, "d", &energy)) return NULL;
    if (!self->mat) { PyErr_SetString(PyExc_RuntimeError, "Material not initialized"); return NULL; }
    return PyFloat_FromDouble(alea_nuc_mat_xs_total(self->mat, energy));
}

static PyObject* PyAleaNucMaterial_xs_absorption(PyAleaNucMaterialObject* self, PyObject* args) {
    double energy;
    if (!PyArg_ParseTuple(args, "d", &energy)) return NULL;
    if (!self->mat) { PyErr_SetString(PyExc_RuntimeError, "Material not initialized"); return NULL; }
    return PyFloat_FromDouble(alea_nuc_mat_xs_absorption(self->mat, energy));
}

static PyObject* PyAleaNucMaterial_xs_elastic(PyAleaNucMaterialObject* self, PyObject* args) {
    double energy;
    if (!PyArg_ParseTuple(args, "d", &energy)) return NULL;
    if (!self->mat) { PyErr_SetString(PyExc_RuntimeError, "Material not initialized"); return NULL; }
    return PyFloat_FromDouble(alea_nuc_mat_xs_elastic(self->mat, energy));
}

static PyObject* PyAleaNucMaterial_mean_free_path(PyAleaNucMaterialObject* self, PyObject* args) {
    double energy;
    if (!PyArg_ParseTuple(args, "d", &energy)) return NULL;
    if (!self->mat) { PyErr_SetString(PyExc_RuntimeError, "Material not initialized"); return NULL; }
    return PyFloat_FromDouble(alea_nuc_mean_free_path(self->mat, energy));
}

static PyObject* PyAleaNucMaterial_sample_distance(PyAleaNucMaterialObject* self, PyObject* args) {
    double energy, xi;
    if (!PyArg_ParseTuple(args, "dd", &energy, &xi)) return NULL;
    if (!self->mat) { PyErr_SetString(PyExc_RuntimeError, "Material not initialized"); return NULL; }
    return PyFloat_FromDouble(alea_nuc_sample_distance(self->mat, energy, xi));
}

static PyObject* PyAleaNucMaterial_sample_nuclide(PyAleaNucMaterialObject* self, PyObject* args) {
    double energy, xi;
    if (!PyArg_ParseTuple(args, "dd", &energy, &xi)) return NULL;
    if (!self->mat) { PyErr_SetString(PyExc_RuntimeError, "Material not initialized"); return NULL; }

    alea_nuc_nuclide_t* out_nuc = NULL;
    int idx = alea_nuc_sample_nuclide(self->mat, energy, xi, &out_nuc);
    if (idx < 0 || !out_nuc) {
        PyErr_SetString(PyExc_RuntimeError, "Failed to sample nuclide");
        return NULL;
    }

    return Py_BuildValue("(is)", idx, out_nuc->zaid);
}

static uint64_t photon_response_mix64(uint64_t value) {
    value += UINT64_C(0x9e3779b97f4a7c15);
    value = (value ^ (value >> 30)) * UINT64_C(0xbf58476d1ce4e5b9);
    value = (value ^ (value >> 27)) * UINT64_C(0x94d049bb133111eb);
    return value ^ (value >> 31);
}

static int photon_response_bin(const double* edges, npy_intp n_bins,
                               double energy) {
    if (!isfinite(energy) || energy < edges[0] || energy > edges[n_bins])
        return -1;
    if (energy == edges[n_bins]) return (int)n_bins - 1;
    npy_intp low = 0, high = n_bins;
    while (low + 1 < high) {
        npy_intp middle = low + (high - low) / 2;
        if (energy < edges[middle]) high = middle;
        else low = middle;
    }
    return (int)low;
}


/* Evaluate each production channel with deterministic bin probabilities when
 * supported, or stratified energy sampling. Exact production responses weight
 * both paths, so rare reaction selection does not waste samples. */
static PyObject* PyAleaNucMaterial_photon_response(
        PyAleaNucMaterialObject* self, PyObject* args, PyObject* kwds) {
    static char* kwlist[] = {"neutron_energies", "photon_edges",
                             "samples_per_channel", "seed", "by_channel",
                             "deterministic_lines", "energy_offset",
                             "component_filter", "parent_mt_filter",
                             "deterministic_only", NULL};
    PyObject* energies_obj;
    PyObject* edges_obj;
    unsigned int samples = 2048;
    unsigned long long seed = 1;
    int by_channel = 0;
    int deterministic_lines = 0;
    unsigned long long energy_offset = 0;
    int component_filter = -1;
    int parent_mt_filter = -1;
    int deterministic_only = 0;
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "OO|IKppKiip", kwlist,
            &energies_obj, &edges_obj, &samples, &seed, &by_channel,
            &deterministic_lines, &energy_offset,
            &component_filter, &parent_mt_filter, &deterministic_only))
        return NULL;
    if (!self->mat) {
        PyErr_SetString(PyExc_RuntimeError, "Material not initialized");
        return NULL;
    }
    if (samples == 0) {
        PyErr_SetString(PyExc_ValueError,
                        "samples_per_channel must be positive");
        return NULL;
    }

    PyArrayObject* energies = (PyArrayObject*)PyArray_FROM_OTF(
        energies_obj, NPY_DOUBLE, NPY_ARRAY_CARRAY_RO);
    PyArrayObject* edges = (PyArrayObject*)PyArray_FROM_OTF(
        edges_obj, NPY_DOUBLE, NPY_ARRAY_CARRAY_RO);
    if (!energies || !edges) {
        Py_XDECREF(energies);
        Py_XDECREF(edges);
        return NULL;
    }
    if (PyArray_NDIM(energies) != 1 || PyArray_NDIM(edges) != 1) {
        Py_DECREF(energies);
        Py_DECREF(edges);
        PyErr_SetString(PyExc_ValueError,
                        "neutron_energies and photon_edges must be one-dimensional");
        return NULL;
    }
    npy_intp n_energies = PyArray_SIZE(energies);
    npy_intp n_edges = PyArray_SIZE(edges);
    if (n_energies == 0 || n_edges < 2 || n_edges - 1 > INT_MAX) {
        Py_DECREF(energies);
        Py_DECREF(edges);
        PyErr_SetString(PyExc_ValueError,
                        "need at least one neutron energy and two photon edges");
        return NULL;
    }
    const double* energy_data = (const double*)PyArray_DATA(energies);
    const double* edge_data = (const double*)PyArray_DATA(edges);
    for (npy_intp i = 0; i < n_energies; i++) {
        if (!isfinite(energy_data[i]) || energy_data[i] < 0.0) {
            Py_DECREF(energies);
            Py_DECREF(edges);
            PyErr_SetString(PyExc_ValueError,
                            "neutron energies must be finite and non-negative");
            return NULL;
        }
    }
    for (npy_intp i = 0; i < n_edges; i++) {
        if (!isfinite(edge_data[i]) || edge_data[i] < 0.0 ||
            (i > 0 && !(edge_data[i] > edge_data[i - 1]))) {
            Py_DECREF(energies);
            Py_DECREF(edges);
            PyErr_SetString(PyExc_ValueError,
                            "photon edges must be finite, non-negative, and strictly increasing");
            return NULL;
        }
    }

    int aggregate_available = 1;
    int native_grid_consistent = 1;
    for (int c = 0; c < self->mat->n_components; c++) {
        const alea_nuc_mat_component_t* component =
            &self->mat->components[c];
        const alea_nuc_nuclide_t* nuc = component->nuclide;
        if (!nuc || component->number_density <= 0.0 ||
            (nuc->n_photon_productions <= 0 &&
             !nuc->total_photon_production_xs))
            continue;
        alea_nuc_photon_production_audit_t audit;
        alea_error_t audit_error =
            PyAleaNucMaterial_default_photon_audit(self, nuc, &audit);
        if (audit_error != ALEA_OK) {
            Py_DECREF(energies);
            Py_DECREF(edges);
            PyErr_Format(PyExc_RuntimeError,
                         "Photon-production audit failed for component %d: %s",
                         c, alea_error_string(audit_error));
            return NULL;
        }
        if (!audit.aggregate_available)
            aggregate_available = 0;
        else if (!audit.native_grid_consistent)
            native_grid_consistent = 0;
    }

    npy_intp n_bins = n_edges - 1;
    npy_intp n_channels = 0;
    npy_intp* channel_offsets = (npy_intp*)calloc(
        (size_t)self->mat->n_components + 1, sizeof(npy_intp));
    if (!channel_offsets) {
        Py_DECREF(energies);
        Py_DECREF(edges);
        return PyErr_NoMemory();
    }
    for (int c = 0; c < self->mat->n_components; c++) {
        const alea_nuc_nuclide_t* nuc = self->mat->components[c].nuclide;
        npy_intp count = nuc && nuc->n_photon_productions > 0
            ? (npy_intp)nuc->n_photon_productions : 0;
        if (n_channels > NPY_MAX_INTP - count) {
            free(channel_offsets);
            Py_DECREF(energies);
            Py_DECREF(edges);
            PyErr_SetString(PyExc_OverflowError,
                            "too many photon-production channels");
            return NULL;
        }
        n_channels += count;
        channel_offsets[c + 1] = n_channels;
    }
    alea_nuc_prepared_photon_production_t* prepared = n_channels > 0
        ? (alea_nuc_prepared_photon_production_t*)calloc(
            (size_t)n_channels, sizeof(*prepared))
        : NULL;
    if (n_channels > 0 && !prepared) {
        free(channel_offsets);
        Py_DECREF(energies);
        Py_DECREF(edges);
        return PyErr_NoMemory();
    }
    for (int c = 0; c < self->mat->n_components; c++) {
        const alea_nuc_mat_component_t* component = &self->mat->components[c];
        const alea_nuc_nuclide_t* nuc = component->nuclide;
        if (!nuc || component->number_density <= 0.0) continue;
        for (int p = 0; p < nuc->n_photon_productions; p++) {
            alea_error_t prepare_error = alea_nuc_prepare_photon_production(
                nuc, &nuc->photon_productions[p],
                &prepared[channel_offsets[c] + p]);
            if (prepare_error != ALEA_OK) {
                free(prepared);
                free(channel_offsets);
                Py_DECREF(energies);
                Py_DECREF(edges);
                PyErr_Format(PyExc_RuntimeError,
                    "photon-production preparation failed for component %d, production %d: %s",
                    c, p, alea_error_string(prepare_error));
                return NULL;
            }
        }
    }
    npy_intp matrix_dims[2] = {n_energies, n_bins};
    npy_intp curve_dims[1] = {n_energies};
    npy_intp channel_dims[3] = {n_channels, n_energies, n_bins};
    PyArrayObject* values = (PyArrayObject*)PyArray_ZEROS(
        2, matrix_dims, NPY_DOUBLE, 0);
    PyArrayObject* variances = (PyArrayObject*)PyArray_ZEROS(
        2, matrix_dims, NPY_DOUBLE, 0);
    PyArrayObject* production_total = (PyArrayObject*)PyArray_ZEROS(
        1, curve_dims, NPY_DOUBLE, 0);
    PyArrayObject* accounted_total = (PyArrayObject*)PyArray_ZEROS(
        1, curve_dims, NPY_DOUBLE, 0);
    PyArrayObject* collision_total = (PyArrayObject*)PyArray_ZEROS(
        1, curve_dims, NPY_DOUBLE, 0);
    PyArrayObject* channel_values = by_channel
        ? (PyArrayObject*)PyArray_ZEROS(3, channel_dims, NPY_DOUBLE, 0)
        : NULL;
    PyArrayObject* channel_variances = by_channel
        ? (PyArrayObject*)PyArray_ZEROS(3, channel_dims, NPY_DOUBLE, 0)
        : NULL;
    size_t* counts = (size_t*)calloc((size_t)n_bins, sizeof(size_t));
    double* probabilities = (double*)malloc((size_t)n_bins * sizeof(double));
    if (!values || !variances || !production_total || !accounted_total ||
        !collision_total || !counts || !probabilities ||
        (by_channel && (!channel_values || !channel_variances))) {
        Py_DECREF(energies);
        Py_DECREF(edges);
        Py_XDECREF(values);
        Py_XDECREF(variances);
        Py_XDECREF(production_total);
        Py_XDECREF(accounted_total);
        Py_XDECREF(collision_total);
        Py_XDECREF(channel_values);
        Py_XDECREF(channel_variances);
        free(channel_offsets);
        free(prepared);
        free(counts);
        free(probabilities);
        return PyErr_NoMemory();
    }

    double* value_data = (double*)PyArray_DATA(values);
    double* variance_data = (double*)PyArray_DATA(variances);
    double* published_data = (double*)PyArray_DATA(production_total);
    double* accounted_data = (double*)PyArray_DATA(accounted_total);
    double* collision_data = (double*)PyArray_DATA(collision_total);
    double* channel_value_data = channel_values
        ? (double*)PyArray_DATA(channel_values) : NULL;
    double* channel_variance_data = channel_variances
        ? (double*)PyArray_DATA(channel_variances) : NULL;
    alea_error_t native_error = ALEA_OK;
    int error_component = -1, error_production = -1;
    npy_intp error_energy = -1;
    sighandler_func old_handler = install_sigint();

    Py_BEGIN_ALLOW_THREADS
    for (npy_intp e = 0; e < n_energies && native_error == ALEA_OK; e++) {
        double incident_energy = energy_data[e];
        collision_data[e] = alea_nuc_mat_xs_total(self->mat, incident_energy);
        for (int c = 0; c < self->mat->n_components; c++) {
            if (component_filter >= 0 && c != component_filter) continue;
            const alea_nuc_mat_component_t* component =
                &self->mat->components[c];
            const alea_nuc_nuclide_t* nuc = component->nuclide;
            if (!nuc || component->number_density <= 0.0) continue;
            if (nuc->total_photon_production_xs)
                published_data[e] += component->number_density *
                    alea_nuc_xs_photon_production_total(
                        nuc, incident_energy);

            for (int p = 0; p < nuc->n_photon_productions; p++) {
                const alea_nuc_photon_production_t* production =
                    &nuc->photon_productions[p];
                if (parent_mt_filter >= 0 &&
                    production->parent_mt != parent_mt_filter) continue;
                const alea_nuc_prepared_photon_production_t* sampler =
                    &prepared[channel_offsets[c] + p];
                double channel_xs =
                    alea_nuc_prepared_photon_production_response(
                        sampler, incident_energy);
                double channel_response =
                    component->number_density * channel_xs;
                if (!(channel_response > 0.0)) continue;
                accounted_data[e] += channel_response;
                if (!nuc->total_photon_production_xs)
                    published_data[e] += channel_response;
                if (!production->spectrum) {
                    native_error = ALEA_ERR_INVALID_STATE;
                    error_component = c;
                    error_production = p;
                    error_energy = e;
                    break;
                }

                if (deterministic_lines) {
                    alea_error_t probability_error =
                        alea_nuc_prepared_photon_bin_probabilities(
                            sampler, incident_energy, edge_data,
                            (size_t)n_bins + 1, probabilities);
                    if (probability_error == ALEA_OK) {
                        for (npy_intp b = 0; b < n_bins; b++) {
                            double fraction = probabilities[b];
                            npy_intp index = e * n_bins + b;
                            value_data[index] += channel_response * fraction;
                            if (by_channel) {
                                npy_intp detailed_index =
                                    ((channel_offsets[c] + p) * n_energies + e) *
                                    n_bins + b;
                                channel_value_data[detailed_index] =
                                    channel_response * fraction;
                            }
                        }
                        continue;
                    }
                    if (probability_error != ALEA_ERR_UNSUPPORTED ||
                        deterministic_only) {
                        native_error = probability_error;
                        error_component = c;
                        error_production = p;
                        error_energy = e;
                        break;
                    }
                }

                memset(counts, 0, (size_t)n_bins * sizeof(size_t));
                const alea_nuc_particle_state_t incident = {
                    ALEA_NUC_PARTICLE_NEUTRON, incident_energy,
                    {0.0, 0.0, 1.0}, 1.0, 0.0
                };
                uint64_t stream_seed = (uint64_t)seed ^
                    photon_response_mix64((uint64_t)e +
                                          (uint64_t)energy_offset + 1) ^
                    photon_response_mix64(((uint64_t)(uint32_t)c << 32) |
                                           (uint32_t)p);
                for (unsigned int s = 0; s < samples; s++) {
                    alea_nuc_rng_t rng;
                    native_error = alea_nuc_rng_init(
                        &rng, stream_seed, s, (uint32_t)c, (uint32_t)p,
                        ALEA_NUC_RNG_COLLISION);
                    if (native_error != ALEA_OK) break;
                    double emitted_energy;
                    if (deterministic_lines) {
                        double mu;
                        bool correlated;
                        native_error = alea_nuc_sample_energy_angle_distribution(
                            production->spectrum, incident_energy,
                            alea_nuc_rng_uniform, &rng, &emitted_energy,
                            &mu, &correlated);
                    } else {
                        alea_nuc_particle_state_t photon;
                        native_error = alea_nuc_sample_prepared_photon_production(
                            sampler, &incident,
                            alea_nuc_rng_uniform, &rng, &photon);
                        emitted_energy = photon.energy;
                    }
                    if (native_error != ALEA_OK) break;
                    int bin = photon_response_bin(edge_data, n_bins,
                                                  emitted_energy);
                    if (bin >= 0) counts[bin]++;
                }
                if (native_error != ALEA_OK) {
                    error_component = c;
                    error_production = p;
                    error_energy = e;
                    break;
                }
                for (npy_intp b = 0; b < n_bins; b++) {
                    double fraction = (double)counts[b] / (double)samples;
                    npy_intp index = e * n_bins + b;
                    npy_intp channel_index = by_channel
                        ? channel_offsets[c] + p : 0;
                    npy_intp detailed_index =
                        (channel_index * n_energies + e) * n_bins + b;
                    value_data[index] += channel_response * fraction;
                    if (by_channel)
                        channel_value_data[detailed_index] =
                            channel_response * fraction;
                    if (samples > 1) {
                        double variance = channel_response * channel_response *
                            fraction * (1.0 - fraction) /
                            (double)(samples - 1);
                        variance_data[index] += variance;
                        if (by_channel)
                            channel_variance_data[detailed_index] = variance;
                    }
                }
            }
            if (native_error != ALEA_OK) break;
        }
        if (alea_interrupted() && native_error == ALEA_OK)
            native_error = ALEA_ERR_INTERRUPTED;
    }
    Py_END_ALLOW_THREADS

    free(counts);
    free(probabilities);
    free(prepared);
    if (restore_sigint(old_handler)) {
        Py_DECREF(energies);
        Py_DECREF(edges);
        Py_DECREF(values);
        Py_DECREF(variances);
        Py_DECREF(production_total);
        Py_DECREF(accounted_total);
        Py_DECREF(collision_total);
        Py_XDECREF(channel_values);
        Py_XDECREF(channel_variances);
        free(channel_offsets);
        return NULL;
    }
    if (native_error != ALEA_OK) {
        Py_DECREF(energies);
        Py_DECREF(edges);
        Py_DECREF(values);
        Py_DECREF(variances);
        Py_DECREF(production_total);
        Py_DECREF(accounted_total);
        Py_DECREF(collision_total);
        Py_XDECREF(channel_values);
        Py_XDECREF(channel_variances);
        free(channel_offsets);
        PyErr_Format(PyExc_RuntimeError,
            "photon response evaluation failed at neutron energy index %zd, component %d, production %d: %s",
            error_energy, error_component, error_production,
            alea_error_string(native_error));
        return NULL;
    }

    for (npy_intp i = 0; i < n_energies * n_bins; i++)
        variance_data[i] = sqrt(variance_data[i]);
    if (by_channel) {
        npy_intp size = PyArray_SIZE(channel_variances);
        for (npy_intp i = 0; i < size; i++)
            channel_variance_data[i] = sqrt(channel_variance_data[i]);
    }

    PyObject* result = NULL;
    PyObject* channel_metadata = NULL;
    if (by_channel) {
        channel_metadata = PyList_New(n_channels);
        if (!channel_metadata) goto result_error;
        for (int c = 0; c < self->mat->n_components; c++) {
            const alea_nuc_mat_component_t* component =
                &self->mat->components[c];
            const alea_nuc_nuclide_t* nuc = component->nuclide;
            if (!nuc) continue;
            for (int p = 0; p < nuc->n_photon_productions; p++) {
                const alea_nuc_photon_production_t* production =
                    &nuc->photon_productions[p];
                PyObject* item = Py_BuildValue(
                    "{s:i,s:s,s:i,s:i,s:i,s:i,s:d}",
                    "component_index", c,
                    "zaid", nuc->zaid,
                    "channel_index", p,
                    "mt", production->mt,
                    "parent_mt", production->parent_mt,
                    "mf", production->mf,
                    "number_density", component->number_density);
                if (!item) {
                    Py_DECREF(channel_metadata);
                    channel_metadata = NULL;
                    goto result_error;
                }
                PyList_SET_ITEM(channel_metadata, channel_offsets[c] + p, item);
            }
        }
    }

    result = PyDict_New();
    if (!result ||
        PyDict_SetItemString(result, "values", (PyObject*)values) < 0 ||
        PyDict_SetItemString(result, "standard_error", (PyObject*)variances) < 0 ||
        PyDict_SetItemString(result, "production_total", (PyObject*)production_total) < 0 ||
        PyDict_SetItemString(result, "accounted_total", (PyObject*)accounted_total) < 0 ||
        PyDict_SetItemString(result, "collision_total", (PyObject*)collision_total) < 0 ||
        PyDict_SetItemString(result, "aggregate_available",
                            aggregate_available ? Py_True : Py_False) < 0 ||
        PyDict_SetItemString(result, "native_grid_consistent",
                            native_grid_consistent ? Py_True : Py_False) < 0 ||
        (by_channel &&
         (PyDict_SetItemString(result, "channel_values", (PyObject*)channel_values) < 0 ||
          PyDict_SetItemString(result, "channel_standard_error",
                               (PyObject*)channel_variances) < 0 ||
          PyDict_SetItemString(result, "channels", channel_metadata) < 0))) {
        Py_XDECREF(result);
        result = NULL;
    }
result_error:
    if (by_channel && !channel_metadata) {
        Py_XDECREF(result);
        result = NULL;
    }
    Py_DECREF(energies);
    Py_DECREF(edges);
    Py_DECREF(values);
    Py_DECREF(variances);
    Py_DECREF(production_total);
    Py_DECREF(accounted_total);
    Py_DECREF(collision_total);
    Py_XDECREF(channel_values);
    Py_XDECREF(channel_variances);
    Py_XDECREF(channel_metadata);
    free(channel_offsets);
    return result;
}

static PyMethodDef PyAleaNucMaterial_methods[] = {
    {"add", (PyCFunction)PyAleaNucMaterial_add, METH_VARARGS,
     "add(nuclide, number_density)\n\nAdd a nuclide component with number density (atoms/barn-cm)."},
    {"add_temperature_mix",
     (PyCFunction)PyAleaNucMaterial_add_temperature_mix, METH_VARARGS,
     "add_temperature_mix(lower, upper, upper_fraction, number_density)\n\n"
     "Add bracketing neutron tables as a correlated temperature mixture."},
    {"xs_total", (PyCFunction)PyAleaNucMaterial_xs_total, METH_VARARGS,
     "xs_total(energy) -> float\n\nMacroscopic total cross section (cm^-1)."},
    {"xs_absorption", (PyCFunction)PyAleaNucMaterial_xs_absorption, METH_VARARGS,
     "xs_absorption(energy) -> float\n\nMacroscopic absorption cross section (cm^-1)."},
    {"xs_elastic", (PyCFunction)PyAleaNucMaterial_xs_elastic, METH_VARARGS,
     "xs_elastic(energy) -> float\n\nMacroscopic elastic scattering cross section (cm^-1)."},
    {"mean_free_path", (PyCFunction)PyAleaNucMaterial_mean_free_path, METH_VARARGS,
     "mean_free_path(energy) -> float\n\nMean free path (cm) at energy (MeV)."},
    {"sample_distance", (PyCFunction)PyAleaNucMaterial_sample_distance, METH_VARARGS,
     "sample_distance(energy, xi) -> float\n\nSample distance to next collision (cm). xi in [0,1)."},
    {"sample_nuclide", (PyCFunction)PyAleaNucMaterial_sample_nuclide, METH_VARARGS,
     "sample_nuclide(energy, xi) -> (index, zaid)\n\nSample which nuclide is hit. xi in [0,1)."},
    {"photon_response", (PyCFunction)PyAleaNucMaterial_photon_response,
     METH_VARARGS | METH_KEYWORDS,
     "photon_response(neutron_energies, photon_edges, samples_per_channel=2048, seed=1, by_channel=False, deterministic_lines=False) -> dict\n\n"
     "Estimate the macroscopic neutron-to-photon response matrix (cm^-1), "
     "optionally preserving individual channel contributions."},
    {NULL}
};

static PyTypeObject PyAleaNucMaterialType = {
    PyVarObject_HEAD_INIT(NULL, 0)
    .tp_name = "pyalea._alea.NucMaterial",
    .tp_doc = PyDoc_STR("Nuclear material composition for transport calculations."),
    .tp_basicsize = sizeof(PyAleaNucMaterialObject),
    .tp_itemsize = 0,
    .tp_flags = Py_TPFLAGS_DEFAULT,
    .tp_new = PyAleaNucMaterial_new,
    .tp_init = (initproc)PyAleaNucMaterial_init,
    .tp_dealloc = (destructor)PyAleaNucMaterial_dealloc,
    .tp_methods = PyAleaNucMaterial_methods,
};

/* ============================================================================
 * Multigroup Python Type
 * ============================================================================ */

typedef struct {
    PyObject_HEAD
    alea_nuc_multigroup_t* mg;
    PyObject* spectrum_callable;   /* Python callable for custom spectrum, or NULL */
} PyAleaMultigroupObject;

/* Trampoline: called from C, dispatches to the Python callable stored in ctx */
static double multigroup_spectrum_trampoline(double E, void* ctx) {
    PyObject* callable = (PyObject*)ctx;
    PyObject* result = PyObject_CallFunction(callable, "d", E);
    if (!result) {
        /* Can't propagate Python exceptions through C; return 0 and print */
        PyErr_Print();
        return 0.0;
    }
    double val = PyFloat_AsDouble(result);
    Py_DECREF(result);
    if (PyErr_Occurred()) {
        PyErr_Print();
        return 0.0;
    }
    return val;
}

static void PyAleaMultigroup_dealloc(PyAleaMultigroupObject* self) {
    if (self->mg) {
        alea_nuc_mg_destroy(self->mg);
    }
    Py_XDECREF(self->spectrum_callable);
    Py_TYPE(self)->tp_free((PyObject*)self);
}

static PyObject* PyAleaMultigroup_new(PyTypeObject* type, PyObject* args, PyObject* kwds) {
    (void)args; (void)kwds;
    PyAleaMultigroupObject* self = (PyAleaMultigroupObject*)type->tp_alloc(type, 0);
    if (self) {
        self->mg = NULL;
        self->spectrum_callable = NULL;
    }
    return (PyObject*)self;
}

static int PyAleaMultigroup_init(PyAleaMultigroupObject* self, PyObject* args, PyObject* kwds) {
    static char* kwlist[] = {"bounds", NULL};
    PyObject* bounds_obj;

    if (!PyArg_ParseTupleAndKeywords(args, kwds, "O", kwlist, &bounds_obj))
        return -1;

    /* Convert bounds list to double array */
    PyObject* bounds_seq = PySequence_Fast(bounds_obj, "bounds must be a sequence");
    if (!bounds_seq) return -1;

    Py_ssize_t n = PySequence_Fast_GET_SIZE(bounds_seq);
    if (n < 2) {
        Py_DECREF(bounds_seq);
        PyErr_SetString(PyExc_ValueError, "Need at least 2 group boundaries");
        return -1;
    }

    double* bounds = (double*)malloc(sizeof(double) * n);
    if (!bounds) {
        Py_DECREF(bounds_seq);
        PyErr_NoMemory();
        return -1;
    }

    for (Py_ssize_t i = 0; i < n; i++) {
        bounds[i] = PyFloat_AsDouble(PySequence_Fast_GET_ITEM(bounds_seq, i));
        if (PyErr_Occurred()) {
            free(bounds);
            Py_DECREF(bounds_seq);
            return -1;
        }
    }
    Py_DECREF(bounds_seq);

    int n_groups = (int)(n - 1);
    alea_nuc_multigroup_t* mg = alea_nuc_mg_create(n_groups, bounds);
    free(bounds);

    if (!mg) {
        PyErr_SetString(PyExc_RuntimeError, "Failed to create multigroup structure");
        return -1;
    }

    if (self->mg) alea_nuc_mg_destroy(self->mg);
    Py_CLEAR(self->spectrum_callable);
    self->mg = mg;
    return 0;
}

static PyObject* PyAleaMultigroup_collapse(PyAleaMultigroupObject* self, PyObject* args) {
    PyAleaNuclideObject* nuc_obj;
    if (!PyArg_ParseTuple(args, "O!", &PyAleaNuclideType, &nuc_obj)) return NULL;
    if (!self->mg) { PyErr_SetString(PyExc_RuntimeError, "Multigroup not initialized"); return NULL; }
    if (!nuc_obj->nuc) { PyErr_SetString(PyExc_RuntimeError, "Nuclide not initialized"); return NULL; }

    alea_error_t err = alea_nuc_mg_collapse(self->mg, nuc_obj->nuc);
    if (err != ALEA_OK) {
        PyErr_Format(PyExc_RuntimeError, "Multigroup collapse failed: %s", alea_error_string(err));
        return NULL;
    }
    Py_RETURN_NONE;
}

static PyObject* PyAleaMultigroup_scatter(PyAleaMultigroupObject* self, PyObject* args) {
    int g_from, g_to;
    if (!PyArg_ParseTuple(args, "ii", &g_from, &g_to)) return NULL;
    if (!self->mg) { PyErr_SetString(PyExc_RuntimeError, "Multigroup not initialized"); return NULL; }
    return PyFloat_FromDouble(alea_nuc_mg_scatter(self->mg, g_from, g_to));
}

static PyObject* PyAleaMultigroup_scatter_adjoint(PyAleaMultigroupObject* self, PyObject* args) {
    int g_from, g_to;
    if (!PyArg_ParseTuple(args, "ii", &g_from, &g_to)) return NULL;
    if (!self->mg) { PyErr_SetString(PyExc_RuntimeError, "Multigroup not initialized"); return NULL; }
    return PyFloat_FromDouble(alea_nuc_mg_scatter_adjoint(self->mg, g_from, g_to));
}

static PyObject* PyAleaMultigroup_sample_scatter(PyAleaMultigroupObject* self, PyObject* args) {
    int g_from, adjoint = 0;
    double xi;
    if (!PyArg_ParseTuple(args, "id|p", &g_from, &xi, &adjoint)) return NULL;
    if (!self->mg) { PyErr_SetString(PyExc_RuntimeError, "Multigroup not initialized"); return NULL; }
    return PyLong_FromLong(alea_nuc_mg_sample_scatter(self->mg, g_from, xi, adjoint));
}

static PyObject* PyAleaMultigroup_set_spectrum(PyAleaMultigroupObject* self, PyObject* args) {
    PyObject* callable;
    if (!PyArg_ParseTuple(args, "O", &callable)) return NULL;
    if (!self->mg) { PyErr_SetString(PyExc_RuntimeError, "Multigroup not initialized"); return NULL; }

    if (callable == Py_None) {
        /* Reset to default spectrum */
        alea_nuc_mg_set_spectrum(self->mg, NULL, NULL);
        Py_CLEAR(self->spectrum_callable);
    } else {
        if (!PyCallable_Check(callable)) {
            PyErr_SetString(PyExc_TypeError, "spectrum must be callable or None");
            return NULL;
        }
        Py_INCREF(callable);
        Py_XDECREF(self->spectrum_callable);
        self->spectrum_callable = callable;
        alea_nuc_mg_set_spectrum(self->mg, multigroup_spectrum_trampoline, callable);
    }
    Py_RETURN_NONE;
}

/* Get group constants as dict */
static PyObject* PyAleaMultigroup_get_data(PyAleaMultigroupObject* self, PyObject* args) {
    (void)args;
    if (!self->mg) { PyErr_SetString(PyExc_RuntimeError, "Multigroup not initialized"); return NULL; }

    int ng = self->mg->n_groups;
    PyObject* result = PyDict_New();
    if (!result) return NULL;

    /* Helper macro to create list from double array */
    #define MG_LIST(name, arr, size) do { \
        PyObject* list = PyList_New(size); \
        if (!list) { Py_DECREF(result); return NULL; } \
        for (int i = 0; i < (size); i++) \
            PyList_SET_ITEM(list, i, PyFloat_FromDouble((arr)[i])); \
        PyDict_SetItemString(result, name, list); \
        Py_DECREF(list); \
    } while(0)

    MG_LIST("bounds", self->mg->bounds, ng + 1);
    MG_LIST("sigma_t", self->mg->sigma_t, ng);
    MG_LIST("sigma_a", self->mg->sigma_a, ng);
    MG_LIST("sigma_s", self->mg->sigma_s, ng);
    MG_LIST("sigma_f", self->mg->sigma_f, ng);
    MG_LIST("nu_sigma_f", self->mg->nu_sigma_f, ng);
    MG_LIST("chi", self->mg->chi, ng);

    /* Scattering matrix as list of lists */
    PyObject* scatter_mat = PyList_New(ng);
    if (!scatter_mat) { Py_DECREF(result); return NULL; }
    for (int g = 0; g < ng; g++) {
        PyObject* row = PyList_New(ng);
        if (!row) { Py_DECREF(scatter_mat); Py_DECREF(result); return NULL; }
        for (int gp = 0; gp < ng; gp++) {
            PyList_SET_ITEM(row, gp, PyFloat_FromDouble(self->mg->scatter[g * ng + gp]));
        }
        PyList_SET_ITEM(scatter_mat, g, row);
    }
    PyDict_SetItemString(result, "scatter_matrix", scatter_mat);
    Py_DECREF(scatter_mat);

    #undef MG_LIST

    dict_set_new(result, "n_groups", PyLong_FromLong(ng));

    return result;
}

static PyObject* PyAleaMultigroup_get_n_groups(PyAleaMultigroupObject* self, void* closure) {
    (void)closure;
    if (!self->mg) { PyErr_SetString(PyExc_RuntimeError, "Multigroup not initialized"); return NULL; }
    return PyLong_FromLong(self->mg->n_groups);
}

static PyGetSetDef PyAleaMultigroup_getsetters[] = {
    {"n_groups", (getter)PyAleaMultigroup_get_n_groups, NULL, "Number of energy groups", NULL},
    {NULL}
};

static PyMethodDef PyAleaMultigroup_methods[] = {
    {"collapse", (PyCFunction)PyAleaMultigroup_collapse, METH_VARARGS,
     "collapse(nuclide)\n\nCollapse continuous-energy cross sections into multigroup constants."},
    {"set_spectrum", (PyCFunction)PyAleaMultigroup_set_spectrum, METH_VARARGS,
     "set_spectrum(fn)\n\n"
     "Set a custom weighting spectrum for group collapse.\n"
     "fn must be a callable taking energy (MeV) and returning the spectrum value,\n"
     "or None to reset to the default (Maxwellian + 1/E + fission)."},
    {"scatter", (PyCFunction)PyAleaMultigroup_scatter, METH_VARARGS,
     "scatter(g_from, g_to) -> float\n\nForward scattering matrix element."},
    {"scatter_adjoint", (PyCFunction)PyAleaMultigroup_scatter_adjoint, METH_VARARGS,
     "scatter_adjoint(g_from, g_to) -> float\n\nAdjoint scattering matrix element."},
    {"sample_scatter", (PyCFunction)PyAleaMultigroup_sample_scatter, METH_VARARGS,
     "sample_scatter(g_from, xi, adjoint=False) -> int\n\nSample outgoing group from scattering."},
    {"get_data", (PyCFunction)PyAleaMultigroup_get_data, METH_NOARGS,
     "get_data() -> dict\n\nGet all multigroup constants as a dictionary."},
    {NULL}
};

static PyTypeObject PyAleaMultigroupType = {
    PyVarObject_HEAD_INIT(NULL, 0)
    .tp_name = "pyalea._alea.Multigroup",
    .tp_doc = PyDoc_STR("Multigroup cross sections and scattering matrix."),
    .tp_basicsize = sizeof(PyAleaMultigroupObject),
    .tp_itemsize = 0,
    .tp_flags = Py_TPFLAGS_DEFAULT,
    .tp_new = PyAleaMultigroup_new,
    .tp_init = (initproc)PyAleaMultigroup_init,
    .tp_dealloc = (destructor)PyAleaMultigroup_dealloc,
    .tp_methods = PyAleaMultigroup_methods,
    .tp_getset = PyAleaMultigroup_getsetters,
};

/* ============================================================================
 * Module-level nucdata functions
 * ============================================================================ */

static PyObject* mod_parse_zaid(PyObject* self, PyObject* args) {
    (void)self;
    const char* zaid;
    if (!PyArg_ParseTuple(args, "s", &zaid)) return NULL;

    int Z, A, meta;
    alea_nuc_table_type_t type;
    alea_error_t err = alea_nuc_parse_zaid(zaid, &Z, &A, &meta, &type);
    if (err != ALEA_OK) {
        PyErr_Format(PyExc_ValueError, "Invalid ZAID '%s': %s", zaid, alea_error_string(err));
        return NULL;
    }

    const char* type_str;
    switch (type) {
        case ALEA_NUC_TABLE_CONTINUOUS_NEUTRON: type_str = "continuous_neutron"; break;
        case ALEA_NUC_TABLE_PHOTOATOMIC:        type_str = "photoatomic"; break;
        case ALEA_NUC_TABLE_PHOTONUCLEAR:       type_str = "photonuclear"; break;
        case ALEA_NUC_TABLE_THERMAL_SAB:        type_str = "thermal_sab"; break;
        case ALEA_NUC_TABLE_ELECTRON:           type_str = "electron"; break;
        default:                                 type_str = "unknown"; break;
    }

    return Py_BuildValue("{s:i, s:i, s:i, s:s}",
        "Z", Z, "A", A, "metastable", meta, "type", type_str);
}

static PyObject* mod_reaction_classify(PyObject* self, PyObject* args) {
    (void)self;
    int mt;
    if (!PyArg_ParseTuple(args, "i", &mt)) return NULL;

    alea_nuc_reaction_class_t cls = alea_nuc_reaction_classify(mt);
    const char* cls_str;
    switch (cls) {
        case ALEA_NUC_RXN_ABSORPTION: cls_str = "absorption"; break;
        case ALEA_NUC_RXN_SCATTER:    cls_str = "scatter"; break;
        case ALEA_NUC_RXN_MULTIPLY:   cls_str = "multiply"; break;
        default:                       cls_str = "unknown"; break;
    }

    return PyUnicode_FromString(cls_str);
}
