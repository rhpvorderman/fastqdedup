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

PyDoc_STRVAR(hamming_distance__doc__,
"hamming_distance($module, string1, string2, /)\n"
"--\n"
"\n"
"Calculates the hamming distance between two strings.\n"
"\n"
"  string1\n"
"    An ASCII string.\n"
"  string2\n"
"    Another ASCII string\n"
"\n"
"Returns an integer representing the hamming distance.\n"
"Raises a ValueError when strings are not of the same length.\n"
"\n");

#define HAMMING_DISTANCE_METHODDEF    \
    {"hamming_distance", (PyCFunction)(void(*)(void))hamming_distance, METH_FASTCALL, \
     hamming_distance__doc__}

PyObject *
hamming_distance(PyObject *module,
                 PyObject *const *args,
                 Py_ssize_t nargs)
{
    if (nargs != 2) {
        PyErr_Format(
            PyExc_TypeError,
            "hamming distance expects exactly two arguments, got %zd", 
            nargs);
        return NULL;
    }
    Py_ssize_t i;
    PyObject *arg; 
    for (i=0; i<2; i+=1) {
        arg = args[i];
        if (!PyUnicode_CheckExact(arg)) {
            PyErr_Format(PyExc_TypeError, 
                         "string%zd must be of type str, got %s.",
                         i + 1, Py_TYPE(arg)->tp_name);
            return NULL;
        }
        if (!(PyUnicode_KIND(arg) == PyUnicode_1BYTE_KIND)) {
            PyErr_Format(PyExc_ValueError, 
                        "string%zd must be ASCII or latin-1 encoded.", 
                        i + 1);
            return NULL;
        }
    }
    PyObject *string1 = args[0];
    PyObject *string2 = args[1];
    uint8_t *string1chars = PyUnicode_1BYTE_DATA(string1);
    uint8_t *string2chars = PyUnicode_1BYTE_DATA(string2);
    Py_ssize_t string1len = PyUnicode_GET_LENGTH(string1);
    Py_ssize_t string2len = PyUnicode_GET_LENGTH(string2);

    if (string1len != string2len) {
        PyErr_SetString(PyExc_ValueError, 
                        "string1 and string2 must be of the same length");
        return NULL;
    }
    Py_ssize_t hamming_distance = 0;
    for (i=0; i<string1len; i+=1) {
        if (string1chars[i] != string2chars[i]) {
            hamming_distance += 1;
        }
    }
    return PyLong_FromSsize_t(hamming_distance);
}

static PyMethodDef _distance_functions[] = {
    HAMMING_DISTANCE_METHODDEF,
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