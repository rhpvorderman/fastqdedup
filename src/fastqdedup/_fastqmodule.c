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