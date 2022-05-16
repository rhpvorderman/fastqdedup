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

#define MAXIMUM_PHRED_SCORE 126
typedef struct {
    PyObject_HEAD 
    size_t total;
    size_t pass;
    double threshold;
    uint8_t phred_offset;
} QualityFilter;

static PyObject *
QualityFilter_new(PyTypeObject *type, PyObject *args, PyObject *kwargs) 
{
    double threshold = 0.0L;
    uint8_t phred_offset = 33;
    const char *kwarg_names[] = {"threshold", "phred_offset", NULL};
    const char *format = "d|$b=:QualityFilter.__new__";
    if (!PyArg_ParseTupleAndKeywords(
        args, kwargs, format, kwarg_names, threshold, phred_offset)) {
            return NULL;
    }
    QualityFilter *self = PyObject_New(QualityFilter, type);
    self->phred_offset = phred_offset;
    self->threshold = threshold;
    self->total = 0;
    self->pass = 0;
    return self;
}

PyDoc_STRVAR(QualityFilter_passes_filter__doc__, 
"passes_filter($self, phred_scores, /)\n"
"--\n"
"\n"
"Check if the \n"
"\n"
"  phred_scores\n"
"    ASCII string with the phred scores.\n"
"\n"
"Returns True if the phred_scores pass, false otherwise.\n"
"\n");

#define HAMMING_DISTANCE_METHODDEF    \
    {"hamming_distance", (PyCFunction)(void(*)(void))hamming_distance, \
     METH_O, QualityFilter_passes_filter__doc__}

static PyObject *
QualityFilter_passes_filter(QualityFilter *self, PyObject *phred_scores) 
{
    if (!PyUnicode_CheckExact(phred_scores)) {
        PyErr_Format(PyExc_TypeError, 
                        "phred_scores must be of type str, got %s.",
                        Py_TYPE(phred_scores)->tp_name);
        return NULL;
    }
    if (!PyUnicode_IS_COMPACT_ASCII(phred_scores)) {
        PyErr_SetString(PyExc_ValueError, 
                        "phred_scores must be ASCII encoded.");
        return NULL;
    }
    double total_score = 0.0;
    uint8_t *scores = PyUnicode_DATA(phred_scores);
    uint8_t score;
    uint8_t phred_offset = self->phred_offset;
    Py_ssize_t length = PyUnicode_GET_LENGTH(phred_scores);
    for (Py_ssize_t i; i<length; i+=1) {
        score = scores[i];
        if (score < phred_offset || score > MAXIMUM_PHRED_SCORE ) {
            return NULL;
        }
    }
}

static struct PyModuleDef _fastq_module = {
    PyModuleDef_HEAD_INIT,
    "_fastq",   /* name of module */
    NULL, /* module documentation, may be NULL */
    -1,
    NULL  /* module methods */
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