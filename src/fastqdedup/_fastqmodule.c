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

#include "score_to_error_rate.h"
#define MAXIMUM_PHRED_SCORE 126

PyDoc_STRVAR(average_error_rate__doc__, 
"average_error_rate($self, phred_scores, /, phred_offset=33)\n"
"--\n"
"\n"
"Returns the average error rate as a float. \n"
"\n"
"  phred_scores\n"
"    ASCII string with the phred scores.\n"
);

#define AVERAGE_ERROR_RATE_METHODDEF    \
    {"average_error_rate", (PyCFunction)(void(*)(void))average_error_rate, \
     METH_VARARGS | METH_KEYWORDS, average_error_rate__doc__}

static PyObject *
average_error_rate(PyObject *module, PyObject *args, PyObject *kwargs) 
{
    PyObject *phred_scores = NULL;
    uint8_t phred_offset = 33;
    const char *kwarg_names[] = {"", "phred_offset", NULL};
    const char *format = "O!|$b=:QualityFilter.__new__";
    if (!PyArg_ParseTupleAndKeywords(
        args, kwargs, format, kwarg_names, 
        &PyUnicode_Type, 
        &phred_scores, 
        &phred_offset)) {
            return NULL;
    }

    if (!PyUnicode_IS_COMPACT_ASCII(phred_scores)) {
        PyErr_SetString(PyExc_ValueError, 
                        "phred_scores must be ASCII encoded.");
        return NULL;
    }
    double total_error_rate = 0.0;
    uint8_t *scores = PyUnicode_DATA(phred_scores);
    uint8_t score;
    uint8_t max_score = MAXIMUM_PHRED_SCORE - phred_offset;
    Py_ssize_t length = PyUnicode_GET_LENGTH(phred_scores);
    for (Py_ssize_t i; i<length; i+=1) {
        score = scores[i] - phred_offset;
        if (score > max_score) {
            PyErr_Format(
                PyExc_ValueError, 
                "Character %c outside of valid phred range %c-%c", 
                scores[i], phred_offset, MAXIMUM_PHRED_SCORE);
            return NULL;
        }
        total_error_rate += SCORE_TO_ERROR_RATE[score];
    }
    double average_error = total_error_rate / (double)length;
    return PyFloat_FromDouble(average_error);
}

static PyMethodDef _fastq_functions[] = {
    AVERAGE_ERROR_RATE_METHODDEF,
    {NULL}
};

static struct PyModuleDef _fastq_module = {
    PyModuleDef_HEAD_INIT,
    "_fastq",   /* name of module */
    NULL, /* module documentation, may be NULL */
    -1,
    _fastq_functions  /* module methods */
};

PyMODINIT_FUNC
PyInit__fastq(void)
{
    PyObject *m;

    m = PyModule_Create(&_fastq_module);
    if (m == NULL) {
        return NULL;
    }
    return m;
}