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

#include <stdint.h>

typedef struct {
    uint32_t alphabet_size;
    uint32_t count;
    void * children[0];
} TrieNode;

#define TRIE_NODE_TERMINAL_FLAG     0x80000000
#define TRIE_NODE_SUFFIX_SIZE_MASK  0x7FFFFFFF
#define TRIE_NODE_SUFFIX_MAX_SIZE   0x7FFFFFFF
#define TRIE_NODE_ALPHABET_MAX_SIZE 254
#define TrieNode_IS_TERMINAL(n) (n->alphabet_size & TRIE_NODE_TERMINAL_FLAG)
#define TrieNode_GET_SUFFIX_SIZE(n) (assert (TrieNode_IS_TERMINAL(n)), \
    n->alphabet_size & TRIE_NODE_SUFFIX_SIZE_MASK)
#define _TrieNode_SET_SUFFIX_SIZE(n, s) (n->alphabet_size = TRIE_NODE_TERMINAL_FLAG | s)
#define TrieNode_GET_SUFFIX(n) (assert (TrieNode_IS_TERMINAL(n)), \
    (uint8_t *)n->children)
#define TrieNode_GET_CHILD(n, i) ((TrieNode *)(n->children[i]))

static inline TrieNode * 
TrieNode_GetChild(TrieNode * parent, size_t index) {
    if TrieNode_IS_TERMINAL(parent){
        return NULL;
    }
    if (index >= parent->alphabet_size) {
        return NULL;
    }
    return (TrieNode *)(parent->children[index]);
}


static TrieNode * 
TrieNode_Resize(TrieNode * trie_node, uint32_t alphabet_size) {
    if (alphabet_size > TRIE_NODE_ALPHABET_MAX_SIZE) {
        PyErr_SetString(PyExc_SystemError, "TrieNode resized with excessive alphabet size");
        return NULL;
    }
    if (alphabet_size == trie_node->alphabet_size) {
        return trie_node;
    }
    size_t old_alphabet_size = TrieNode_IS_TERMINAL(trie_node) ? 0 : trie_node->alphabet_size;
    size_t new_size = sizeof(TrieNode) + sizeof(TrieNode *) * alphabet_size;
    TrieNode * new_node = PyMem_Realloc(trie_node, new_size);
    if (new_node == NULL) {
        PyErr_NoMemory();
        return NULL;
    }
    new_node->alphabet_size = alphabet_size;
    if (alphabet_size > old_alphabet_size) {
        size_t empty_size = (alphabet_size - old_alphabet_size) * sizeof(TrieNode *);
        memset(new_node->children + old_alphabet_size, 0, empty_size);
    }
    return new_node;
}

static void
TrieNode_Destroy(TrieNode * trie_node) {
    if (trie_node == NULL) {
        return;
    }
    if (!TrieNode_IS_TERMINAL(trie_node)){
        // Destroy children first
        uint32_t alphabet_size = trie_node->alphabet_size;
        uint32_t i;
        for (i=0; i < alphabet_size; i += 1){
            TrieNode_Destroy(TrieNode_GET_CHILD(trie_node, i));
        }
    }
    PyMem_Free(trie_node);
    return;
}

static TrieNode *
TrieNode_NewLeaf(uint8_t * suffix, uint32_t suffix_size, uint32_t sequence_count) {
    if (suffix_size > TRIE_NODE_SUFFIX_MAX_SIZE) {
        PyErr_SetString(PyExc_SystemError, "TrieNode initialized with excessive suffix size");
        return NULL;
    }
    size_t leaf_node_size = sizeof(TrieNode) + suffix_size;
    TrieNode *  new = PyMem_Malloc(leaf_node_size);
    if (new == NULL) {
        PyErr_NoMemory();
        return NULL;
    }
    _TrieNode_SET_SUFFIX_SIZE(new, suffix_size);
    new->count = sequence_count;
    memcpy(TrieNode_GET_SUFFIX(new), suffix, suffix_size);
    return new;
}

static int
TrieNode_AddSequence(TrieNode ** trie_node,
                     uint8_t * sequence, 
                     uint32_t sequence_size, 
                     uint8_t *alphabet_size, 
                     uint8_t * charmap,
                     uint8_t * alphabet,
                     uint32_t sequence_count) {
    TrieNode * this_node = trie_node[0];
    if (this_node == NULL) {
        trie_node[0] = TrieNode_NewLeaf(sequence, sequence_size, sequence_count);
        return 0;
    }
    if (TrieNode_IS_TERMINAL(this_node)) {
        uint32_t suffix_size = TrieNode_GET_SUFFIX_SIZE(this_node);
        if (sequence_size == suffix_size) {
            if (memcmp(TrieNode_GET_SUFFIX(this_node), sequence, sequence_size) == 0){
                this_node->count += sequence_count;
                return 0;
            }
        }
        uint8_t * suffix = TrieNode_GET_SUFFIX(this_node);
        uint8_t * tmp = PyMem_Malloc(suffix_size);
        if (tmp == NULL) {
            PyErr_NoMemory();
            return -1;
        }
        memcpy(tmp, suffix, suffix_size);
        this_node->alphabet_size = 0;
        uint32_t count = this_node->count;
        this_node->count = 0;
        int ret = TrieNode_AddSequence(
            trie_node, tmp, suffix_size, alphabet_size, charmap, alphabet, count);
        this_node = trie_node[0];
        PyMem_Free(tmp);
        if (ret != 0){
            return ret;
        }
    }
    if (sequence_size == 0) {
        this_node->count += sequence_count;
        return 0;
    }

    uint8_t character = sequence[0];
    uint8_t node_index = charmap[character];
    if (node_index == 255) {
        node_index = *alphabet_size;
        charmap[character] = node_index;
        alphabet[node_index] = character;
        *alphabet_size = node_index + 1;
    }
    if (node_index >= this_node->alphabet_size) {
        TrieNode * new_node = TrieNode_Resize(this_node, node_index + 1);
        if (new_node == NULL) {
            return -1;
        }
        this_node = new_node;
        trie_node[0] = this_node;
    }
    return TrieNode_AddSequence((TrieNode **)&(this_node->children[node_index]), 
                                 sequence + 1, 
                                 sequence_size - 1, 
                                 alphabet_size, 
                                 charmap,
                                 alphabet,
                                 sequence_count);
}

static ssize_t
TrieNode_DeleteSequence(
    TrieNode ** trie_node,
    uint8_t * sequence,
    uint32_t sequence_size,
    const uint8_t * charmap)
{
    TrieNode * this_node = trie_node[0];
    uint32_t count;
    if (TrieNode_IS_TERMINAL(this_node)) {
        uint32_t suffix_size = TrieNode_GET_SUFFIX_SIZE(this_node);
        if (!(sequence_size == suffix_size)) {
            return -1;
        }
        if (!memcmp(TrieNode_GET_SUFFIX(this_node), sequence, sequence_size) == 0){
            return -1;
        }
        count = this_node->count;
        PyMem_Free(this_node);
        *trie_node = NULL;
        return count;
    }

    if (sequence_size == 0) {
        uint32_t count = this_node->count;
        if (count == 0) {
            return -1;
        }
        this_node->count = 0;
        return count;
    }

    uint8_t character = sequence[0];
    uint8_t node_index = charmap[character];
    if (node_index == 255) {
        return -1;
    }

    if (node_index >= this_node->alphabet_size) {
       return -1;
    }
    TrieNode * next_node = TrieNode_GetChild(this_node, node_index);
    if (next_node == NULL) {
        return -1;
    }
    ssize_t ret = TrieNode_DeleteSequence(&(this_node->children[node_index]), sequence + 1, sequence_size - 1, charmap);
    if (ret > -1) {
        for (size_t i=0; i < this_node->alphabet_size; i+=1) {
            if (TrieNode_GET_CHILD(this_node, i) != NULL) {
                return ret;
            }
        }
        // All children are null
        if (this_node->count) {
            trie_node[0] = TrieNode_NewLeaf(NULL, 0, this_node->count);
        }
        else {
            trie_node[0] = NULL;
        }
        PyMem_Free(this_node);
    }
    return ret;
}

static ssize_t
TrieNode_FindNearest(
    TrieNode * trie_node, 
    const uint8_t * sequence,
    uint32_t sequence_length,
    int max_distance, 
    const uint8_t * charmap, 
    uint8_t *buffer,
    const uint8_t * alphabet) 
{
    if (max_distance < 0) {
        return 0;
    }
    if TrieNode_IS_TERMINAL(trie_node) {
        uint32_t suffix_length = TrieNode_GET_SUFFIX_SIZE(trie_node); 
        if (sequence_length != suffix_length) {
            // Hamming is technically only valid for sequences with the same
            // length. 
            return 0;
        }
        uint8_t * suffix = TrieNode_GET_SUFFIX(trie_node);
        for (uint32_t i=0; i < suffix_length; i++) {
            if (sequence[i] != suffix[i]) {
                max_distance -= 1;
                if (max_distance < 0) {
                    return 0;
                }
            }
        }
        if (buffer) {
            memcpy(buffer, suffix, suffix_length);
        }
        return trie_node->count;
    }
    if (sequence_length == 0) {
        return trie_node->count;
    }

    uint8_t character = sequence[0];
    uint8_t node_index = charmap[character];
    uint8_t *new_buffer = NULL; 
    if (buffer) {
        new_buffer = buffer + 1;
    }
    TrieNode * child = TrieNode_GetChild(trie_node, node_index);
    if (child == NULL) {
        // Mismatch, try all children and deduct a point from the max distance.
        max_distance -= 1;
        for (uint32_t i=0; i < trie_node->alphabet_size; i++) {
            child = TrieNode_GET_CHILD(trie_node, i);
            if (child == NULL) {
                continue;
            }
            if (buffer) {
                buffer[0] = alphabet[i];
            }
            int ret = TrieNode_FindNearest(
                child, sequence + 1, sequence_length -1, max_distance, charmap, 
                new_buffer, alphabet);
            if (ret) {
                return ret; 
            }
        }
        return 0;
    }
    // Found a match, continue the computation with the child node.
    if (buffer) {
        buffer[0] = character;
    }
    return TrieNode_FindNearest(
        child, sequence + 1, sequence_length -1, max_distance, charmap, 
        new_buffer, alphabet);
}

/**
 * @brief Get a sequence from the trie node. 
 * 
 * This function takes the first sequence based on the order of the 
 * provided alphabet.
 * 
 * @param trie_node 
 * @param alphabet 
 * @param buffer A buffer where 
 * @param buffer_size 
 * @return ssize_t The size of the sequence. 0 if none found. -1 if buffer was
 *                 too small.
 */
static ssize_t
TrieNode_GetSequence(
    TrieNode *trie_node, 
    const uint8_t *alphabet, 
    uint8_t *buffer, 
    uint32_t buffer_size) 
{
    if (TrieNode_IS_TERMINAL(trie_node)) {
        uint32_t suffix_size = TrieNode_GET_SUFFIX_SIZE(trie_node);
        if (suffix_size > buffer_size) {
            return -1;
        }
        uint8_t *suffix = TrieNode_GET_SUFFIX(trie_node);
        memcpy(buffer, suffix, suffix_size);
        return suffix_size;
    }
    // Node has children but we cannot add these to the buffer anymore.
    if (buffer_size < 1) {
        return -1;
    }
    uint32_t alphabet_size = trie_node->alphabet_size;
    uint32_t i;
    TrieNode * child;
    for (i=0; i<alphabet_size; i+=1) {
        child = TrieNode_GET_CHILD(trie_node, i);
        if (child == NULL) {
            continue;
        }
        buffer[0] = alphabet[i];
        ssize_t ret = TrieNode_GetSequence(child, alphabet, buffer + 1, buffer_size - 1);
        if (ret == -1) {
            return ret;
        }
        return 1 + ret;
    }
    // No children found.
    if (trie_node->count > 0) {
        return 0;
    }
    // Fail when the node count is 0 so this does not store a sequence.
    return -1;
}


typedef struct {
    PyObject_HEAD
    uint8_t charmap[256];
    uint8_t alphabet[256];
    uint8_t alphabet_size;
    Py_ssize_t number_of_sequences;
    uint32_t max_sequence_size;
    TrieNode * root;
} Trie;

static void 
Trie_Dealloc(Trie * self) {
    TrieNode_Destroy(self->root);
}

static PyObject *
Trie__new__(PyTypeObject *type, PyObject *args, PyObject *kwargs) {
    char * keywords[] = {NULL};
    const char *format = "|:Trie.__new__";
    if (!PyArg_ParseTupleAndKeywords(
            args, kwargs, format, keywords)) {
        return NULL;
    }
    Trie * new_trie = PyObject_New(Trie, type);
    new_trie->alphabet_size = 0;
    memset(new_trie->charmap, 255, 256);
    memset(new_trie->alphabet, 0, 256);
    new_trie->root = NULL;
    new_trie->number_of_sequences = 0;
    new_trie->max_sequence_size = 0; 
    return (PyObject *)new_trie;
}

PyDoc_STRVAR(Trie_add_sequence__doc__,
"add_sequence($self, sequence, /)\n"
"--\n"
"\n"
"Adds a sequence to the trie.\n"
"\n"
"  sequence\n"
"    An ASCII string.\n"
"\n");

#define TRIE_ADD_SEQUENCE_METHODDEF    \
    {"add_sequence", (PyCFunction)(void(*)(void))Trie_add_sequence, METH_O, \
     Trie_add_sequence__doc__}

static PyObject * 
Trie_add_sequence(Trie *self, PyObject * sequence) {
    if (!PyUnicode_CheckExact(sequence)) {
        PyErr_Format(PyExc_TypeError, "Sequence must be a str, got %s", 
            Py_TYPE(sequence)->tp_name);
        return NULL;
    }
    if (!PyUnicode_IS_COMPACT_ASCII(sequence)) {
        PyErr_SetString(PyExc_ValueError, 
                        "Sequence must consist only of ASCII characters");
        return NULL;
    }
    uint8_t * seq = PyUnicode_DATA(sequence);
    Py_ssize_t seq_size = PyUnicode_GET_LENGTH(sequence);
    if (seq_size > TRIE_NODE_SUFFIX_MAX_SIZE) {
        PyErr_Format(
            PyExc_ValueError, 
            "Sequences larger than %d can not be stored in the Trie",
            TRIE_NODE_SUFFIX_MAX_SIZE);
        return NULL;
    }
    if (TrieNode_AddSequence(&(self->root), seq, seq_size, &(self->alphabet_size), self->charmap, self->alphabet, 1) == 0) {
        self->number_of_sequences += 1;
        if (seq_size > self->max_sequence_size) {
            self->max_sequence_size = seq_size;
        }
        Py_RETURN_NONE;
    }
    return NULL;
}

PyDoc_STRVAR(Trie_contains_sequence__doc__,
"contains_sequence($self, sequence, /, max_hamming_distance=0)\n"
"--\n"
"\n"
"Check if a sequence is present in the trie.\n"
"\n"
"Optionally check if a similar sequence is present at the specified\n"
"maximum hamming distance.\n"
"Sequences with unequal size are considered unequal.\n"
"\n"
"  sequence\n"
"    An ASCII string.\n"
"  max_hamming_distance\n"
"    The maximal Hamming distance\n"
"\n");

#define TRIE_CONTAINS_SEQUENCE_METHODDEF    \
    {"contains_sequence", (PyCFunction)(void(*)(void))Trie_contains_sequence, \
    METH_VARARGS | METH_KEYWORDS, Trie_contains_sequence__doc__}

static PyObject *
Trie_contains_sequence(Trie *self, PyObject *args, PyObject* kwargs) {
    PyObject * sequence = NULL;
    int max_distance = 0;
    char * keywords[] = {"", "max_hamming_distance", NULL};
    const char *format = "O|i:Trie.contains_sequence";
    if (!PyArg_ParseTupleAndKeywords(
            args, kwargs, format, keywords,
            &sequence, &max_distance)) {
        return NULL;
    }
    if (!PyUnicode_CheckExact(sequence)) {
        PyErr_Format(PyExc_TypeError, "Sequence must be a str, got %s",
            Py_TYPE(sequence)->tp_name);
        return NULL;
    }
    if (!PyUnicode_IS_COMPACT_ASCII(sequence)) {
        PyErr_SetString(PyExc_ValueError, "sequence must contain only ASCII characters");
        return NULL;
    }
    uint8_t * seq = (uint8_t *)PyUnicode_DATA(sequence);
    Py_ssize_t seq_size = PyUnicode_GET_LENGTH(sequence);
    if (seq_size > TRIE_NODE_SUFFIX_MAX_SIZE) {
        PyErr_Format(
            PyExc_ValueError,
            "Sequences larger than %d can not be stored in the Trie",
            TRIE_NODE_SUFFIX_MAX_SIZE);
        return NULL;
    }
    int ret = TrieNode_FindNearest(self->root, seq, seq_size, max_distance, 
              self->charmap, NULL, NULL);
    return PyBool_FromLong(ret);
}

PyDoc_STRVAR(Trie_pop_cluster__doc__,
"pop_cluster($self, max_hamming_distance, /)\n"
"--\n"
"\n"
"Find a cluster of sequences within the same hamming distance and remove them from the trie.\n"
"\n"
"Optionally check if a similar sequence is present at the specified\n"
"maximum hamming distance.\n"
"Sequences with unequal size are considered unequal.\n"
"\n"
"  max_hamming_distance\n"
"    The maximal Hamming distance\n"
"\n");

#define TRIE_POP_CLUSTER_METHODDEF    \
    {"pop_cluster", (PyCFunction)(void(*)(void))Trie_pop_cluster, \
    METH_O, Trie_pop_cluster__doc__}

static PyObject *
Trie_pop_cluster(Trie *self, PyObject *max_hamming_distance) {
    if (!PyLong_CheckExact(max_hamming_distance)) {
        PyErr_Format(PyExc_TypeError, 
                     "max_hamming_distance expected an int not %s", 
                     Py_TYPE(max_hamming_distance)->tp_name);
        return NULL;
    }
    int max_distance = PyLong_AsLong(max_hamming_distance);
    if (max_distance < 0) {
        PyErr_SetString(PyExc_ValueError, 
                        "max_hamming distance should be larger than 0");
        return NULL;
    }
    if (self->root == NULL) {
        PyErr_SetString(PyExc_LookupError, "No sequences left in Trie.");
        return NULL;
    }
    uint32_t buffer_size = self->max_sequence_size;
    uint8_t *buffer = PyMem_Malloc(buffer_size);
    if (buffer == NULL) {
        return PyErr_NoMemory();
    }
    ssize_t sequence_size = TrieNode_GetSequence(self->root, self->alphabet, 
                                                buffer, buffer_size);
    if (sequence_size == - 1){
        PyErr_SetString(PyExc_RuntimeError, "Incorrect buffer size used.");
        PyMem_Free(buffer);
        return NULL;
    }
    PyObject * first_sequence_obj = PyUnicode_New(sequence_size, 127);
    if (first_sequence_obj == NULL) {
        PyMem_Free(buffer);
        return PyErr_NoMemory();
    }
    memcpy(PyUnicode_DATA(first_sequence_obj), buffer, sequence_size);
    buffer = PyMem_Realloc(buffer, sequence_size);
    if (buffer == NULL){
        Py_DECREF(first_sequence_obj);
        return PyErr_NoMemory();
    }
    uint8_t * template_sequence = PyUnicode_DATA(first_sequence_obj);
    uint32_t template_size = sequence_size;
    ssize_t template_count = TrieNode_DeleteSequence(&(self->root), 
                              template_sequence, template_size, self->charmap);
    if (template_count == -1) {
        PyErr_SetString(PyExc_RuntimeError, "Retrieved undeletable sequence.");
        PyMem_Free(buffer);
        Py_DECREF(first_sequence_obj);
        return NULL;
    }
    PyObject * cluster = PyList_New(1);
    PyObject * tup = PyTuple_New(2);
    PyTuple_SET_ITEM(tup, 0, PyLong_FromSsize_t(template_count));
    PyTuple_SET_ITEM(tup, 1, first_sequence_obj);
    PyList_SET_ITEM(cluster, 0, tup);
    if (max_distance == 0) {
        return cluster;
    }
    Py_ssize_t cluster_index = 0;
    Py_ssize_t cluster_size = 1;
    PyObject * template_tup;
    PyObject * template;
    PyObject * sequence;
    ssize_t sequence_count;
    ssize_t ret;
    while ((cluster_index != cluster_size) && (self->root != NULL)) {
        template_tup = PyList_GET_ITEM(cluster, cluster_index);
        template = PyTuple_GET_ITEM(template_tup, 1);
        template_sequence = PyUnicode_DATA(template);
        template_size = PyUnicode_GET_LENGTH(template);
        sequence_count = TrieNode_FindNearest(self->root, template_sequence, template_size,
            max_distance, self->charmap, buffer, self->alphabet);
        if (sequence_count) {
            sequence = PyUnicode_New(sequence_size, 127);
            memcpy(PyUnicode_DATA(sequence), buffer, template_size);
            ret = TrieNode_DeleteSequence(&(self->root), buffer,
                                          template_size, self->charmap);
            if (ret == -1) {
                PyErr_SetString(PyExc_RuntimeError, "Retrieved undeletable sequence.");
                PyMem_Free(buffer);
                Py_DECREF(sequence);
                Py_DECREF(cluster);
                return NULL;
            }
            tup = PyTuple_New(2);
            PyTuple_SET_ITEM(tup, 0, PyLong_FromSsize_t(sequence_count));
            PyTuple_SET_ITEM(tup, 1, sequence);
            PyList_Append(cluster, tup);
            cluster_size += 1;
        }
        else {
            // Use the next sequence in the appended list as template. This way
            // we traverse all sequences to create a cluster.
            cluster_index +=1;
        }
    }
    return cluster;
}


static PyMethodDef Trie_methods[] = {
    TRIE_ADD_SEQUENCE_METHODDEF,
    TRIE_CONTAINS_SEQUENCE_METHODDEF,
    TRIE_POP_CLUSTER_METHODDEF,
    {NULL}
};

static PyTypeObject Trie_Type = {
    PyVarObject_HEAD_INIT(NULL, 0)
    .tp_name = "_trie.Trie",
    .tp_basicsize = sizeof(Trie),
    .tp_dealloc = (destructor)Trie_Dealloc,
    .tp_new = Trie__new__,
    .tp_methods = Trie_methods,
};


static struct PyModuleDef _trie_module = {
    PyModuleDef_HEAD_INIT,
    "_trie",   /* name of module */
    NULL, /* module documentation, may be NULL */
    -1,
    NULL  /* module methods */
};

PyMODINIT_FUNC
PyInit__trie(void)
{
    PyObject *m;

    m = PyModule_Create(&_trie_module);
    if (m == NULL)
        return NULL;

    if (PyType_Ready(&Trie_Type) < 0)
        return NULL;
    PyObject * TrieType = (PyObject *)&Trie_Type;
    Py_INCREF(TrieType);
    if (PyModule_AddObject(m, "Trie", TrieType) < 0)
        return NULL;
    return m;
}