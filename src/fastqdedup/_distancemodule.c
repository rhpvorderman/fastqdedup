// Copyright (C) 2022 Leiden University Medical Center
// This file is part of fastqdedup
//
// fastqdedup is free software: you can redistribute it and/or modify
// it under the terms of the GNU Affero General Public License as
// published by the Free Software Foundation, either version 3 of the
// License, or (at your option) any later version.
//
// fastqdedup is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Affero General Public License for more details.
//
// You should have received a copy of the GNU Affero General Public License
// along with fastqdedup.  If not, see <https://www.gnu.org/licenses/

#define PY_SSIZE_T_CLEAN
#include <Python.h>

#include "distances.h"


PyDoc_STRVAR(within_distance__doc__,
"within_distance($module, string1, string2, /, max_distance, use_edit_distance=False)\n"
"--\n"
"\n"
"Calculates whether two strings are within the specified distance.\n"
"\n"
"  string1\n"
"    An ASCII string.\n"
"  string2\n"
"    Another ASCII string\n"
"  max_distance\n"
"     The maximum distance\n"
"  use_edit_distance\n"
"    Use edit (Levenshtein) distance instead of Hamming distance"
"\n"
"Returns an integer representing the hamming distance.\n"
"Raises a ValueError when strings are not of the same length.\n"
"\n");

#define WITHIN_DISTANCE_METHODDEF    \
    {"within_distance", (PyCFunction)(void(*)(void))within_distance, \
    METH_VARARGS | METH_KEYWORDS, within_distance__doc__}

PyObject *
within_distance(PyObject *module,
                PyObject *args,
                PyObject *kwargs)
{
    PyObject *string1 = NULL;
    PyObject *string2 = NULL;
    int max_distance = 0;
    int use_edit_distance = 0;
    char *keywords[] = {"", "", "max_distance", "use_edit_distance", NULL};
    char *format = "O!O!i|p:within_distance";
    if (!PyArg_ParseTupleAndKeywords(
            args, kwargs, format, keywords,
            &PyUnicode_Type, &string1, &PyUnicode_Type, &string2,
            &max_distance, &use_edit_distance)) {
                return NULL;
    }

    if (!(PyUnicode_KIND(string1) == PyUnicode_1BYTE_KIND)) {
        PyErr_SetString(PyExc_ValueError, 
                        "string1 must be ASCII or latin-1 encoded.");
        return NULL;
    }
    if (!(PyUnicode_KIND(string2) == PyUnicode_1BYTE_KIND)) {
        PyErr_SetString(PyExc_ValueError, 
                        "string2 must be ASCII or latin-1 encoded.");
        return NULL;
    }

    uint8_t *string1chars = PyUnicode_1BYTE_DATA(string1);
    uint8_t *string2chars = PyUnicode_1BYTE_DATA(string2);
    Py_ssize_t string1len = PyUnicode_GET_LENGTH(string1);
    Py_ssize_t string2len = PyUnicode_GET_LENGTH(string2);
    ssize_t ret;

    if (use_edit_distance) {
        ret = within_edit_distance(string1chars, string1len,
                                   string2chars, string2len, 
                                   max_distance);

    }
    else {
        ret = within_hamming_distance(string1chars, string1len,
                                      string2chars, string2len, 
                                      max_distance);
    }
    return PyBool_FromLong(ret);
}


static PyMethodDef _distance_functions[] = {
    WITHIN_DISTANCE_METHODDEF,
    {NULL}
};

static struct PyModuleDef _distance_module = {
    PyModuleDef_HEAD_INIT,
    "_distance",   /* name of module */
    NULL, /* module documentation, may be NULL */
    -1,
    _distance_functions  /* module methods */
};

PyMODINIT_FUNC
PyInit__distance(void)
{
    PyObject *m;

    m = PyModule_Create(&_distance_module);
    if (m == NULL) {
        return NULL;
    }
    return m;
}