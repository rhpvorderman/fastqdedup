#define PY_SSIZE_T_CLEAN
#include <Python.h>

#include <stdint.h>

typedef struct {
    uint32_t alphabet_size;
    uint32_t count;
    void * children[0];
} TrieNode;

#define TRIE_NODE_TERMINAL_FLAG     0x80000000
#define TRIE_NODE_SUFFIX_SIZE_MASK  0x0FFFFFFF
#define TRIE_NODE_SUFFIX_MAX_SIZE   0x0FFFFFFF
#define TRIE_NODE_ALPHABET_MAX_SIZE 254
#define TrieNode_IS_TERMINAL(n) (n->alphabet_size & TRIE_NODE_TERMINAL_FLAG)
#define TrieNode_GET_SUFFIX_SIZE(n) (assert (TrieNode_IS_TERMINAL(n)), \
    n->alphabet_size & TRIE_NODE_SUFFIX_SIZE_MASK)
#define _TrieNode_SET_SUFFIX_SIZE(n, s) (n->alphabet_size = TRIE_NODE_TERMINAL_FLAG | s)
#define TrieNode_GET_SUFFIX(n) (assert (TrieNode_IS_TERMINAL(n)), \
    (uint8_t *)n->children)

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

#define TrieNode_GET_CHILD(n, i) ((TrieNode *)(n->children[i]))
static TrieNode *
TrieNode_New(uint32_t alphabet_size) {
    if (alphabet_size > TRIE_NODE_ALPHABET_MAX_SIZE) {
        PyErr_SetString(PyExc_SystemError, "TrieNode initialized with excessive alphabet size");
        return NULL;
    }
    size_t trie_node_size = sizeof(TrieNode) + sizeof(TrieNode *) * alphabet_size;
    TrieNode * new = PyMem_Malloc(trie_node_size);
    if (new == NULL) {
        PyErr_NoMemory();
        return NULL;
    }
    memset(new, 0, trie_node_size);
    new->alphabet_size = alphabet_size;
    return new; 
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
TrieNode_NewLeaf(uint8_t * suffix, uint32_t suffix_size) {
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
    new->count = 1;
    memcpy(TrieNode_GET_SUFFIX(new), suffix, suffix_size);
    return new;
}

static int
TrieNode_AddSequence(TrieNode ** trie_node,
                     uint8_t * sequence, 
                     uint32_t sequence_size, 
                     uint8_t *alphabet_size, 
                     uint8_t * charmap) {
    TrieNode * this_node = trie_node[0];
    if (sequence_size == 0) {
        return 0;
    }
    if (TrieNode_IS_TERMINAL(this_node)) {
        uint32_t suffix_size = TrieNode_GET_SUFFIX_SIZE(this_node);
        if (sequence_size == suffix_size) {
            if (memcmp(TrieNode_GET_SUFFIX(this_node), sequence, sequence_size) == 0){
                this_node->count += 1;
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
        int ret = TrieNode_AddSequence(trie_node, tmp, suffix_size, alphabet_size, charmap);
        this_node = trie_node[0];
        PyMem_Free(tmp);
        if (ret != 0){
            return ret;
        }
    }

    uint8_t character = sequence[0];
    uint8_t node_index = charmap[character];
    if (node_index == 255) {
        node_index = *alphabet_size;
        charmap[character] = node_index;
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
    TrieNode * next_node = TrieNode_GetChild(this_node, node_index);
    if (next_node == NULL) {
        next_node = TrieNode_NewLeaf(sequence + 1, sequence_size -1);
        if (next_node == NULL) {
            return -1;
        }
        this_node->children[node_index] = next_node;
        return 0;
    }
    return TrieNode_AddSequence(&(this_node->children[node_index]), sequence + 1, sequence_size - 1, alphabet_size, charmap);
}

typedef struct {
    PyObject_HEAD
    uint8_t charmap[256];
    uint8_t alphabet_size;
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
    new_trie->root = TrieNode_New(0);
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
    if (TrieNode_AddSequence(&(self->root), seq, seq_size, &(self->alphabet_size), self->charmap) == 0) {
        Py_RETURN_NONE;
    }
    return NULL;
}

static PyMethodDef Trie_methods[] = {
    TRIE_ADD_SEQUENCE_METHODDEF,
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