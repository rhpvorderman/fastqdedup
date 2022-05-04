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

#define TRIE_NODE_ALPHABET_MAX_SIZE 254

/**
 * @brief Simple struct to store an alphabet and conversion operations. 
 * 
 * from_index stores the actual alphabet and should be initialized with zeroes.
 * to_index stores the indexes at each character's position, with 255 for 
 * unknown characters. It should be initialized with 255s (0xff) 
 * size stores the alphabet size.
 * */
typedef struct {
    uint8_t from_index[256];
    uint8_t to_index[256];
    uint8_t size;
} Alphabet;

static int 
Alphabet_InitializeFromString(Alphabet *alphabet, uint8_t *string) {
    memset(alphabet->from_index, 0, 256);
    memset(alphabet->to_index, 255, 256);
    alphabet->size = 0;
    uint8_t c;
    uint8_t current_to_index; 
    while (1) {
        if (alphabet->size > TRIE_NODE_ALPHABET_MAX_SIZE) {
            PyErr_SetString(PyExc_ValueError, "Maximum alphabet length exceeded");
            return -1;
        }
        c = string[alphabet->size];
        if (c == 0) {
            return 0;
        }
        current_to_index = alphabet->to_index[c];
        if (current_to_index != 255) {
            PyErr_Format(PyExc_ValueError,
                "Alphabet should consist of unique characters." 
                "Character %c was repeated. ",
                c);
            return -1;
        }
        alphabet->to_index[c] = alphabet->size;
        alphabet->from_index[alphabet->size] = c;
        alphabet->size += 1;
    }

}

/**
 * @brief A node in a trie.
 * 
 * Each node stores its children in children. These are of type 'TrieNode *'. 
 * 
 * The alphabet_size signifies how many children are stored. This does not have
 * to reflect the alphabet size in the application. Say the alphabet is ACGT 
 * (DNA letters) with indexes 0, 1, 2 and 3. If we have a node that has only 
 * one child 'C', we only need to have children[1] populated. children[0] can
 * be set to NULL. We only need to reserve memory for 2 adresses (0 and 1) and 
 * can set the alphabet_size to 2. Implicitly assuming that any index of 2 and
 * higher is NULL. The TrieNode_GetChild function has codes this behaviour. 
 * When nodes are sparsely populated this saves a lot of memory.
 * 
 * The highest bit of alphabet_size is used to store whether the node is 
 * terminal, in other words a leaf node. When a node is a leaf node, it has a
 * suffix, which is stored at the address of children. This suffix is of type
 * uint8_t. It allows for storing sequences efficiently as in a radix tree, at
 * the terminal ends of the trie.
 * When this is the case, the remaining 31 bits of the alphabet_size store the
 * suffix size.
 * 
 * A count higher than 0 signifies that there are sequences that have this node 
 * as last node. The amount of sequences is stored in this count. Nodes with
 * a count are not necessarily terminal (a leaf node) as sequences stored in 
 * the trie may be of unequal size.
 * 
 * The advantage of storing both normal and leaf nodes in the same struct is 
 * that it makes it easy to expand or shrink the tree as desired.
 */
typedef struct {
    uint32_t alphabet_size;
    uint32_t count;
    void *children[0];
} TrieNode;

#define TRIE_NODE_TERMINAL_FLAG     0x80000000
#define TRIE_NODE_SUFFIX_SIZE_MASK  0x7FFFFFFF
#define TRIE_NODE_SUFFIX_MAX_SIZE   0x7FFFFFFF
#define TrieNode_IS_TERMINAL(n) (n->alphabet_size & TRIE_NODE_TERMINAL_FLAG)
#define TrieNode_GET_SUFFIX_SIZE(n) (assert (TrieNode_IS_TERMINAL(n)), \
    n->alphabet_size & TRIE_NODE_SUFFIX_SIZE_MASK)
#define _TrieNode_SET_SUFFIX_SIZE(n, s) (n->alphabet_size = TRIE_NODE_TERMINAL_FLAG | s)
#define TrieNode_GET_SUFFIX(n) (assert (TrieNode_IS_TERMINAL(n)), \
    (uint8_t *)n->children)
#define TrieNode_GET_CHILD(n, i) ((TrieNode *)(n->children[i]))

/**
 * @brief Get the child of parent for a given index. Always return NULL when
 *        child is node present or the node is terminal.
 */
static inline TrieNode *
TrieNode_GetChild(TrieNode *parent, size_t index) {
    if TrieNode_IS_TERMINAL(parent){
        return NULL;
    }
    if (index >= parent->alphabet_size) {
        return NULL;
    }
    return (TrieNode *)(parent->children[index]);
}


/**
 * @brief Resize a trie_node so it can contain a given alphabet size. Returns
 *       a new pointer.
 */
static TrieNode *
TrieNode_Resize(TrieNode *trie_node, uint32_t alphabet_size) {
    if (alphabet_size > TRIE_NODE_ALPHABET_MAX_SIZE) {
        PyErr_SetString(PyExc_SystemError, 
                        "TrieNode resized with excessive alphabet size");
        return NULL;
    }
    if (alphabet_size == trie_node->alphabet_size) {
        return trie_node;
    }
    size_t old_alphabet_size = 
        TrieNode_IS_TERMINAL(trie_node) ? 0 : trie_node->alphabet_size;
    size_t new_size = sizeof(TrieNode) + sizeof(TrieNode *) * alphabet_size;
    TrieNode *new_node = PyMem_Realloc(trie_node, new_size);
    if (new_node == NULL) {
        PyErr_NoMemory();
        return NULL;
    }
    new_node->alphabet_size = alphabet_size;
    if (alphabet_size > old_alphabet_size) {
        // Set the memory block containing the pointers to the new children to NULL
        size_t empty_size = (alphabet_size - old_alphabet_size) * sizeof(TrieNode *);
        memset(new_node->children + old_alphabet_size, 0, empty_size);
    }
    return new_node;
}

/**
 * @brief Free a TrieNode and all its children.
 */
static void
TrieNode_Destroy(TrieNode *trie_node) {
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

/**
 * @brief Create a new TrieNode. The TrieNode is a leaf node (terminal node).
 * 
 * @param suffix The sequence to store in the node, may be NULL.
 * @param suffix_size The suffix size, in case of a NULL suffix, this should be 0.
 * @param sequence_count The amount of sequences stored in the node.
 * @return TrieNode* Pointer to the new leaf node.
 */
static TrieNode *
TrieNode_NewLeaf(uint8_t *suffix, uint32_t suffix_size, uint32_t sequence_count) {
    if (suffix_size > TRIE_NODE_SUFFIX_MAX_SIZE) {
        PyErr_SetString(PyExc_SystemError, "TrieNode initialized with excessive suffix size");
        return NULL;
    }
    size_t leaf_node_size = sizeof(TrieNode) + suffix_size;
    TrieNode *new = PyMem_Malloc(leaf_node_size);
    if (new == NULL) {
        PyErr_NoMemory();
        return NULL;
    }
    _TrieNode_SET_SUFFIX_SIZE(new, suffix_size);
    new->count = sequence_count;
    memcpy(TrieNode_GET_SUFFIX(new), suffix, suffix_size);
    return new;
}

/**
 * @brief Add a sequence to the TrieNode at the trie_node_address. May resize
 *        the TrieNode accordingly.
 * 
 * @param trie_node_address Pointer to a TrieNode pointer. 
 * @param sequence The sequence to add.
 * @param sequence_size The size of the sequence.
 * @param sequence_count How many sequence should be added.
 * @param alphabet The alphabet. The alphabet may be updated if new characters
 *                 are in the sequence.
 * @return 0 on success, -1 on failure.
 */
static int
TrieNode_AddSequence(TrieNode **trie_node_address,
                     uint8_t *sequence, 
                     uint32_t sequence_size, 
                     uint32_t sequence_count,
                     Alphabet *alphabet) {
    TrieNode *this_node = trie_node_address[0];
    if (this_node == NULL) {
        trie_node_address[0] = TrieNode_NewLeaf(sequence, sequence_size, sequence_count);
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
        // Store the suffix in a temporary space and add  it to this node.
        uint8_t *suffix = TrieNode_GET_SUFFIX(this_node);
        uint8_t *tmp = PyMem_Malloc(suffix_size);
        if (tmp == NULL) {
            PyErr_NoMemory();
            return -1;
        }
        memcpy(tmp, suffix, suffix_size);
        this_node->alphabet_size = 0;
        uint32_t count = this_node->count;
        this_node->count = 0;
        int ret = TrieNode_AddSequence(
            trie_node_address, tmp, suffix_size, count, alphabet);
        this_node = trie_node_address[0];
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
    uint8_t node_index = alphabet->to_index[character];
    if (node_index == 255) { // Letter not present in the alphabet. Add it.
        node_index = alphabet->size;
        alphabet->to_index[character] = node_index;
        alphabet->from_index[node_index] = character;
        alphabet->size = node_index + 1;
    }

    if (node_index >= this_node->alphabet_size) {
        TrieNode *new_node = TrieNode_Resize(this_node, node_index + 1);
        if (new_node == NULL) {
            return -1;
        }
        this_node = new_node;
        trie_node_address[0] = this_node;
    }
    return TrieNode_AddSequence((TrieNode **)&(this_node->children[node_index]), 
                                 sequence + 1, 
                                 sequence_size - 1, 
                                 sequence_count,
                                 alphabet);
}

/**
 * @brief Delete the sequence from the TrieNode at this address. If this leaves
 *        the TrieNode with only empty (NULL) children, the address is set
 *        to NULL and the trie node is freed.
 * 
 * @param trie_node_address a pointer to a TrieNode pointer.
 * @param sequence the sequence to remove
 * @param sequence_size the size of the sequence
 * @param alphabet the alphabet.
 * @return ssize_t the number of removed sequences or -1 for error.
 */
static ssize_t
TrieNode_DeleteSequence(
    TrieNode **trie_node_address,
    uint8_t *sequence,
    uint32_t sequence_size,
    Alphabet *alphabet)
{
    TrieNode *this_node = trie_node_address[0];
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
        *trie_node_address = NULL;
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
    uint8_t node_index = alphabet->to_index[character];
    if (node_index == 255) {
        return -1;
    }

    if (node_index >= this_node->alphabet_size) {
       return -1;
    }
    TrieNode *next_node = TrieNode_GetChild(this_node, node_index);
    if (next_node == NULL) {
        return -1;
    }
    ssize_t ret = TrieNode_DeleteSequence(
        (TrieNode **)&(this_node->children[node_index]), 
        sequence + 1, sequence_size - 1, alphabet);
    
    // Make sure the tree is properly pruned so there are no dead-end nodes
    // that will mess up the search algorithms. 
    if (ret > -1) {
        for (size_t i=0; i < this_node->alphabet_size; i+=1) {
            if (TrieNode_GET_CHILD(this_node, i) != NULL) {
                return ret;
            }
        }
        // All children are NULL. If the node has a count it can be converted
        // into a leave node. If not, it can be freed to.
        if (this_node->count) {
            trie_node_address[0] = TrieNode_NewLeaf(NULL, 0, this_node->count);
        }
        else {
            trie_node_address[0] = NULL;
        }
        PyMem_Free(this_node);
    }
    return ret;
}

/**
 * @brief Find the nearest sequence in the trie at the specified hamming 
 *        distance. Optionally store this sequence in a buffer.
 * 
 * @param trie_node the trie to search.
 * @param sequence the sequence
 * @param sequence_length the length of the sequence
 * @param max_distance the maximum hamming distance. Use 0 for exact match.
 * @param alphabet the used alphabet in the trie.
 * @param buffer a buffer to store the found sequence in. Can be set to NULL 
 *               for no storage. The buffer should have a memory size of 
 *               sequence_length.
 * @return ssize_t the count of the found sequences. If 0 this 
 *                indicates the sequence was not found.
 */
static ssize_t
TrieNode_FindNearest(
    TrieNode *trie_node, 
    const uint8_t *sequence,
    uint32_t sequence_length,
    int max_distance, 
    Alphabet *alphabet, 
    uint8_t *buffer) 
{
    if TrieNode_IS_TERMINAL(trie_node) {
        uint32_t suffix_length = TrieNode_GET_SUFFIX_SIZE(trie_node); 
        if (sequence_length != suffix_length) {
            // Hamming is technically only valid for sequences with the same
            // length. 
            return 0;
        }
        uint8_t *suffix = TrieNode_GET_SUFFIX(trie_node);
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
    uint8_t node_index = alphabet->to_index[character];
    uint8_t *new_buffer = NULL; 
    if (buffer) {
        new_buffer = buffer + 1;
    }
    TrieNode *child = TrieNode_GetChild(trie_node, node_index);
    ssize_t result;
    if (child != NULL) {
        // Found a match, continue the computation with the child node.
        if (buffer) {
            buffer[0] = character;
        }
        result = TrieNode_FindNearest(
            child, sequence + 1, sequence_length -1, max_distance, alphabet, 
            new_buffer);
        if (result) {
            return result;
        }
    }
    // Mismatch, try all children and deduct a point from the max distance.
    max_distance -= 1;
    if (max_distance < 0) {
        return 0;
    }
    for (uint32_t i=0; i < trie_node->alphabet_size; i++) {
        if (i == node_index){
            continue;  // Already tried this route.
        }
        child = TrieNode_GET_CHILD(trie_node, i);
        if (child == NULL) {
            continue;
        }
        if (buffer) {
            buffer[0] = alphabet->from_index[i];
        }
        result = TrieNode_FindNearest(
            child, sequence + 1, sequence_length -1, max_distance, alphabet, 
            new_buffer);
        if (result) {
            return result; 
        }
    }
    // All children searched but could not find anything.
    return 0;
}

/**
 * @brief Get a sequence from the trie node. 
 * 
 * This function takes the first sequence based on the order of the 
 * provided alphabet.
 * 
 * @param trie_node The trie to traverse.
 * @param alphabet The used alphabet.
 * @param buffer A buffer where the sequence is stored.
 * @param buffer_size The size of the buffer.
 * @return ssize_t The size of the sequence. 0 if none found. -1 if buffer was
 *                 too small.
 */
static ssize_t
TrieNode_GetSequence(
    TrieNode *trie_node, 
    const Alphabet *alphabet, 
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
    TrieNode *child;
    for (i=0; i<alphabet_size; i+=1) {
        child = TrieNode_GET_CHILD(trie_node, i);
        if (child == NULL) {
            continue;
        }
        buffer[0] = alphabet->from_index[i];
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

static size_t 
TrieNode_GetMemorySize(TrieNode *trie_node) {
    if (trie_node == NULL) {
        return 0;
    }
    size_t size = sizeof(TrieNode);  // Basic size
    // Calculate size allocated to children member
    if (TrieNode_IS_TERMINAL(trie_node)) {
        return size + TrieNode_GET_SUFFIX_SIZE(trie_node);
    }
    size += (sizeof(TrieNode *) * trie_node->alphabet_size);

    // Calculate size of children
    TrieNode *child;
    for (uint32_t i=0; i<trie_node->alphabet_size; i+=1) {
        child = TrieNode_GET_CHILD(trie_node, i);
        size += TrieNode_GetMemorySize(child);
    }
    return size;
}

static void
TrieNode_GetStats(TrieNode *trie_node,
                  size_t layer,
                  size_t alphabet_size,
                  size_t *stats) {
    if (trie_node == NULL) {
        return;
    }
    size_t *layer_stats = stats + (alphabet_size + 1) * layer;
    if (TrieNode_IS_TERMINAL(trie_node)) {
        // Store terminal nodes at index 0
        layer_stats[0] += 1;
        return;
    }
    layer_stats[trie_node->alphabet_size] += 1;
    for (size_t i=0; i<trie_node->alphabet_size; i+=1) {
        TrieNode_GetStats(TrieNode_GET_CHILD(trie_node, i),
                          layer + 1,
                          alphabet_size,
                          stats);
    }
    return;
}

typedef struct {
    PyObject_HEAD
    Alphabet alphabet;
    Py_ssize_t number_of_sequences;
    uint32_t max_sequence_size;
    TrieNode *root;
    uint8_t *sequence_buffer;
    Py_ssize_t sequence_buffer_size;
} Trie;

static void 
Trie_Dealloc(Trie *self) {
    TrieNode_Destroy(self->root);
    PyMem_Free(self->sequence_buffer);
    Py_TYPE(self)->tp_free((PyObject *)self);
}

static PyObject *
Trie__new__(PyTypeObject *type, PyObject *args, PyObject *kwargs) {
    PyObject *alphabet = NULL;
    char *keywords[] = {"alphabet", NULL};
    const char *format = "|O!:Trie.__new__";
    if (!PyArg_ParseTupleAndKeywords(
            args, kwargs, format, keywords,
            &PyUnicode_Type, &alphabet)) {
        return NULL;
    }
    uint8_t *alphabet_string = (uint8_t *)"";
    if (alphabet != NULL) {
        if (!PyUnicode_IS_COMPACT_ASCII(alphabet)) {
            PyErr_SetString(PyExc_ValueError, "Alphabet should be an ASCII string.");
            return NULL;
        }
        alphabet_string = PyUnicode_1BYTE_DATA(alphabet);
    }
    Trie *self = PyObject_New(Trie, type);
    self->root = NULL;
    self->number_of_sequences = 0;
    self->max_sequence_size = 0;
    self->sequence_buffer = NULL;
    self->sequence_buffer_size = 0;
    if (Alphabet_InitializeFromString(&self->alphabet, alphabet_string) != 0) {
        Py_DECREF(self);
        return NULL;
    }
    return (PyObject *)self;
}

static PyObject *
Trie_get_alphabet(Trie *self, void *closure) {
    return PyUnicode_DecodeLatin1((char *)(self->alphabet.from_index), 
                                   self->alphabet.size, NULL);
}

static PyObject *
Trie_get_number_of_sequences(Trie *self, void *closure) {
    return PyLong_FromSsize_t(self->number_of_sequences);
}

static PyGetSetDef Trie_properties[] = {
    {"alphabet", (getter)Trie_get_alphabet, NULL, NULL, 
    "The alphabet this trie uses."},
    {"number_of_sequences", (getter)Trie_get_number_of_sequences, NULL, NULL, 
     "The number of sequences stored in the trie."},
    {NULL}
};

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
Trie_add_sequence(Trie *self, PyObject *sequence) {
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
    uint8_t *seq = PyUnicode_DATA(sequence);
    Py_ssize_t seq_size = PyUnicode_GET_LENGTH(sequence);
    if (seq_size > TRIE_NODE_SUFFIX_MAX_SIZE) {
        PyErr_Format(
            PyExc_ValueError, 
            "Sequences larger than %d can not be stored in the Trie",
            TRIE_NODE_SUFFIX_MAX_SIZE);
        return NULL;
    }
    if (TrieNode_AddSequence(&(self->root), seq, seq_size, 1, &(self->alphabet)) == 0) {
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
    PyObject *sequence = NULL;
    int max_distance = 0;
    char *keywords[] = {"", "max_hamming_distance", NULL};
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
    uint8_t *seq = (uint8_t *)PyUnicode_DATA(sequence);
    Py_ssize_t seq_size = PyUnicode_GET_LENGTH(sequence);
    if (seq_size > TRIE_NODE_SUFFIX_MAX_SIZE) {
        PyErr_Format(
            PyExc_ValueError,
            "Sequences larger than %d can not be stored in the Trie",
            TRIE_NODE_SUFFIX_MAX_SIZE);
        return NULL;
    }
    ssize_t ret = TrieNode_FindNearest(self->root, seq, seq_size, max_distance, 
                                       &(self->alphabet), NULL);
    return PyBool_FromLong(ret);
}

PyDoc_STRVAR(Trie_pop_cluster__doc__,
"pop_cluster($self, max_hamming_distance, /)\n"
"--\n"
"\n"
"Find a cluster of sequences within the same hamming distance and remove them\n"
"from the trie.\n"
"A list of tuples of (count, sequence) is returned.\n"
"\n"
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

    // By making the buffer max_sequence_size we should not run into buffer 
    // overflows.
    if (self->sequence_buffer_size != self->max_sequence_size) {
        uint8_t *tmp = PyMem_Realloc(self->sequence_buffer, self->max_sequence_size);
        if (tmp == NULL) {
            return PyErr_NoMemory();
        }
        self->sequence_buffer = tmp;
        self->sequence_buffer_size = self->max_sequence_size;
    }
    uint8_t *buffer = self->sequence_buffer;
    uint32_t buffer_size = self->sequence_buffer_size;

    // Get an initial sequence to build the cluster around.
    ssize_t sequence_size = TrieNode_GetSequence(self->root, &(self->alphabet), 
                                                 buffer, buffer_size);
    if (sequence_size == -1){
        PyErr_SetString(PyExc_RuntimeError, "Incorrect buffer size used.");
        return NULL;
    }
    // PyUnicode_New + memcpy is faster than PyUnicode_DecodeXXXX family.
    PyObject *first_sequence_obj = PyUnicode_New(sequence_size, 127);
    if (first_sequence_obj == NULL) {
        return PyErr_NoMemory();
    }
    memcpy(PyUnicode_DATA(first_sequence_obj), buffer, sequence_size);

    // The sequence is saved in a Python object and can be removed from the 
    // trie.
    uint8_t *template_sequence = PyUnicode_DATA(first_sequence_obj);
    uint32_t template_size = sequence_size;
    ssize_t template_count = TrieNode_DeleteSequence(&(self->root), 
                              template_sequence, template_size, &(self->alphabet));
    if (template_count == -1) {
        PyErr_SetString(PyExc_RuntimeError, "Retrieved undeletable sequence.");
        Py_DECREF(first_sequence_obj);
        return NULL;
    }
    self->number_of_sequences -= template_count;
    // Initiate a cluster from the obtained sequence.
    PyObject *cluster = PyList_New(1);
    PyObject *tup = PyTuple_New(2);
    PyTuple_SET_ITEM(tup, 0, PyLong_FromSsize_t(template_count));
    PyTuple_SET_ITEM(tup, 1, first_sequence_obj);
    PyList_SET_ITEM(cluster, 0, tup);
    if (max_distance == 0) {
        return cluster;
    }

    // For each sequence in the cluster.
    // - Find a neighbour at the specified Hamming distance. 
    // - If found:
    //      - Store it in a PyObject *. Store the count. 
    //      - Add it to the cluster
    //      - Delete the neighbour from the trie.
    //      - Look for a new neighbour and repeat the above steps.
    // - If not found, go to the next sequence in the cluster and repeat the
    //   above steps.
    // This way we find the neighbours for all sequences and keep expanding the
    // cluster until no other sequences can be found. 
    Py_ssize_t cluster_index = 0;
    Py_ssize_t cluster_size = 1;
    PyObject *template_tup;
    PyObject *template;
    PyObject *sequence;
    ssize_t sequence_count;
    ssize_t deleted_count;
    while ((cluster_index != cluster_size) && (self->root != NULL)) {
        template_tup = PyList_GET_ITEM(cluster, cluster_index);
        template = PyTuple_GET_ITEM(template_tup, 1);
        template_sequence = PyUnicode_DATA(template);
        template_size = PyUnicode_GET_LENGTH(template);
        sequence_count = TrieNode_FindNearest(self->root, template_sequence, template_size,
            max_distance, &(self->alphabet), buffer);
        if (sequence_count) {
            sequence = PyUnicode_New(sequence_size, 127);
            memcpy(PyUnicode_DATA(sequence), buffer, template_size);
            deleted_count = TrieNode_DeleteSequence(&(self->root), buffer,
                                          template_size, &(self->alphabet));
            if (deleted_count == -1) {
                PyErr_SetString(PyExc_RuntimeError, "Retrieved undeletable sequence.");
                Py_DECREF(sequence);
                Py_DECREF(cluster);
                return NULL;
            }
            self->number_of_sequences -= deleted_count;
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

PyDoc_STRVAR(Trie_memory_size__doc__,
"memory_size($self)\n"
"--\n"
"\n"
"Traverse the trie to get the memory size.\n");

#define TRIE_MEMORY_SIZE_METHODDEF    \
    {"memory_size", (PyCFunction)(void(*)(void))Trie_memory_size, METH_NOARGS, \
     Trie_memory_size__doc__}

static PyObject *
Trie_memory_size(Trie *self, PyObject *Py_UNUSED(ignore)) {
    size_t memory_size = TrieNode_GetMemorySize(self->root);
    return PyLong_FromSize_t(memory_size);
}


PyDoc_STRVAR(Trie_raw_stats__doc__,
"raw_stats($self)\n"
"--\n"
"\n"
"Traverse the trie and for each layer create a list with \n"
"at index 0 the amount of terminal nodes. At the other indexes the amount of\n"
"nodes for index alphabet size\n."
"Returns a list of above lists.");

#define TRIE_RAW_STATS_METHODDEF    \
    {"raw_stats", (PyCFunction)(void(*)(void))Trie_raw_stats, METH_NOARGS, \
     Trie_raw_stats__doc__}

static PyObject *
Trie_raw_stats(Trie *self, PyObject *Py_UNUSED(ignore)) {
    size_t layer_size = self->alphabet.size + 1;
    size_t number_of_layers = self->max_sequence_size + 1;
    size_t stats_size = number_of_layers * layer_size;
    size_t *stats = PyMem_Calloc(stats_size, sizeof(size_t));
    if (stats == NULL) {
        return PyErr_NoMemory();
    }
    
    TrieNode_GetStats(self->root, 0, self->alphabet.size, stats);

    PyObject *return_val = PyList_New(number_of_layers);
    if (return_val == NULL) {
        PyMem_Free(stats);
        return PyErr_NoMemory();
    }

    PyObject *layer_list;
    size_t *layer_stats;
    for (size_t i=0; i<(number_of_layers); i+=1) {
        layer_list = PyList_New(layer_size);
        layer_stats = stats + layer_size * i;
        if (layer_list == NULL) {
            Py_DECREF(return_val);
            PyMem_Free(stats);
            return PyErr_NoMemory();
        }
        for (size_t j=0; j<layer_size; j+=1) {
            PyList_SET_ITEM(layer_list, j, PyLong_FromSize_t(layer_stats[j]));
        }
        PyList_SET_ITEM(return_val, i, layer_list);
    }
    PyMem_Free(stats);
    return return_val;
}

static PyMethodDef Trie_methods[] = {
    TRIE_ADD_SEQUENCE_METHODDEF,
    TRIE_CONTAINS_SEQUENCE_METHODDEF,
    TRIE_POP_CLUSTER_METHODDEF,
    TRIE_MEMORY_SIZE_METHODDEF,
    TRIE_RAW_STATS_METHODDEF,
    {NULL}
};

static PyTypeObject Trie_Type = {
    PyVarObject_HEAD_INIT(NULL, 0)
    .tp_name = "_trie.Trie",
    .tp_basicsize = sizeof(Trie),
    .tp_dealloc = (destructor)Trie_Dealloc,
    .tp_new = Trie__new__,
    .tp_methods = Trie_methods,
    .tp_getset = Trie_properties,
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
    PyObject *TrieType = (PyObject *)&Trie_Type;
    Py_INCREF(TrieType);
    if (PyModule_AddObject(m, "Trie", TrieType) < 0)
        return NULL;
    return m;
}