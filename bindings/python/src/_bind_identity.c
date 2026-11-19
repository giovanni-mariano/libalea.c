// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
// SPDX-License-Identifier: MPL-2.0
/* Included by pyalea_binding.c. */
static int identity_options(const char* profile,alea_fingerprint_options_t* options) {
    *options=(alea_fingerprint_options_t){sizeof(*options),ALEA_FINGERPRINT_LEGACY};
    if(!strcmp(profile,"semantic"))options->profile=ALEA_FINGERPRINT_SEMANTIC;
    else if(strcmp(profile,"legacy")){PyErr_SetString(PyExc_ValueError,"profile must be legacy or semantic");return -1;}
    return 0;
}
static PyObject* PyAleaSystem_fingerprints(PyAleaSystemObject* self,PyObject* args,PyObject* kwargs) {
    const char* profile="legacy";
    static char* keywords[]={"profile",NULL};
    if(!PyArg_ParseTupleAndKeywords(args,kwargs,"|$s",keywords,&profile))return NULL;
    alea_fingerprint_options_t options;
    if(identity_options(profile,&options))return NULL;
    if(!self->sys) {
        PyErr_SetString(PyExc_RuntimeError,"System not initialized");
        return NULL;
    }
    alea_fingerprint_set_t* set;
    Py_BEGIN_ALLOW_THREADS
    set=alea_fingerprint_set_create(self->sys,&options);
    Py_END_ALLOW_THREADS
    const alea_fingerprints_t* f=alea_fingerprint_set_values(set);
    if(!f) {
        PyErr_SetString(PyExc_RuntimeError,alea_error());
        return NULL;
    }
    PyObject* cells=PyDict_New();
    PyObject* universes=PyDict_New();
    PyObject* result=NULL;
    if(!cells || !universes)goto done;
    for(int u=0;u<2;u++) {
        size_t count=u?f->universe_count:f->cell_count;
        alea_fingerprint_t* entries=u?f->universes:f->cells;
        for(size_t i=0;i<count;i++) {
            char geometry[65],content[65];
            alea_digest_hex(&entries[i].geometry,geometry);
            alea_digest_hex(&entries[i].content,content);
            PyObject* key=PyLong_FromLong(entries[i].id);
            PyObject* value=Py_BuildValue("{s:s,s:s}","geometry",geometry,"content",content);
            int failed=!key || !value || PyDict_SetItem(u?universes:cells,key,value);
            Py_XDECREF(key);
            Py_XDECREF(value);
            if(failed)goto done;
        }
    }
    if(options.profile==ALEA_FINGERPRINT_LEGACY)
        result=Py_BuildValue("{s:I,s:O,s:O}","version",f->version,"cells",cells,"universes",universes);
    else
        result=Py_BuildValue("{s:I,s:s,s:O,s:O}","version",f->version,"profile",profile,"cells",cells,"universes",universes);
    done:
    Py_XDECREF(cells);
    Py_XDECREF(universes);
    alea_fingerprint_set_destroy(set);
    return result;
}
static PyObject* PyAleaSystem_compare(PyAleaSystemObject* self,PyObject* args,PyObject* kwargs) {
    PyObject* other;
    const char* profile="legacy";
    static char* keywords[]={"other","profile",NULL};
    if(!PyArg_ParseTupleAndKeywords(args,kwargs,"O|$s",keywords,&other,&profile))return NULL;
    alea_fingerprint_options_t options;
    if(identity_options(profile,&options))return NULL;
    if(!PyObject_TypeCheck(other,&PyAleaSystemType)) {
        PyErr_SetString(PyExc_TypeError,"other must be a System");
        return NULL;
    }
    if(!self->sys || !((PyAleaSystemObject*)other)->sys) {
        PyErr_SetString(PyExc_RuntimeError,"System not initialized");
        return NULL;
    }
    alea_fingerprint_set_t *a,*b;
    alea_model_diff_t* diff;
    Py_BEGIN_ALLOW_THREADS
    a=alea_fingerprint_set_create(self->sys,&options);
    b=alea_fingerprint_set_create(((PyAleaSystemObject*)other)->sys,&options);
    diff=alea_fingerprint_set_diff(a,b);
    Py_END_ALLOW_THREADS
    alea_fingerprint_set_destroy(a);
    alea_fingerprint_set_destroy(b);
    if(!diff) {
        PyErr_SetString(PyExc_RuntimeError,"Cannot compare models");
        return NULL;
    }
    PyObject* result=PyList_New(0);
    static const char* kinds[]= {
        "unchanged","geometry","content","renumbered","removed","added","ambiguous"
    };
    for(size_t i=0;result && i<diff->count;i++) {
        alea_model_change_t* c=&diff->changes[i];
        PyObject* value=Py_BuildValue("{s:s,s:i,s:i,s:s}","kind",c->universe?"universe":"cell","old_id",c->old_id,"new_id",c->new_id,"change",kinds[c->kind]);
        if(!value || PyList_Append(result,value)) {
            Py_CLEAR(result);
        }
        Py_XDECREF(value);
    }
    alea_model_diff_destroy(diff);
    return result;
}
