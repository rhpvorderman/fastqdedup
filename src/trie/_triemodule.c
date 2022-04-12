#define PY_SSIZE_T_CLEAN
#include <Python.h>

#include <stdint.h>

typedef struct {
    uint32_t alphabet_size;
    uint32_t count;
    TrieNode * children[0];
} TrieNode;

#define TRIE_NODE_TERMINAL_FLAG     0x80000000
#define TRIE_NODE_SUFFIX_SIZE_MASK  0x0FFFFFFF
#define TRIE_NODE_SUFFIX_MAX_SIZE   0x0FFFFFFF
#define TRIE_NODE_ALPHABET_MAX_SIZE 254
#define TrieNode_IS_TERMINAL(n) (n->alphabet_size & TRIE_NODE_TERMINAL_FLAG)
#define TrieNode_GET_SUFFIX_SIZE(n) (assert (TrieNode_IS_TERMINAL(n)), \
    n->alphabet_size & TRIE_NODE_SUFFIX_SIZE_MASK)
#define _TrieNode_SET_SUFFIX_SIZE(n, s) (n->alphabet_size = TRIE_NODE_TERMINAL_FLAG & s)
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
    return parent->children[index];
}

static TrieNode *
TrieNode_New(uint32_t alphabet_size) {
    if (alphabet_size > TRIE_NODE_ALPHABET_MAX_SIZE) {
        PyErr_SetString(PyExc_SystemError, "TrieNode initialized with excessive alphabet size");
        return NULL;
    }
    size_t trie_node_size = sizeof(TrieNode) + sizeof(TrieNode *) * alphabet_size;
    TrieNode * new = PyMem_Malloc(trie_node_size);
    if (new == NULL) {
        return PyErr_NoMemory();
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
        return PyErr_NoMemory();
    }
    new_node->alphabet_size = alphabet_size;
    if (alphabet_size > old_alphabet_size) {
        size_t empty_size = (alphabet_size - old_alphabet_size) * sizeof(TrieNode *);
        memset(new_node->children + old_alphabet_size, 0, empty_size);
    }
    return new_node;
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
        return PyErr_NoMemory();
    }
    _TrieNode_SET_SUFFIX_SIZE(new, suffix_size);
    new->count = 1;
    memcpy(TrieNode_GET_SUFFIX(new), suffix, suffix_size);
    return new;
}

static int
TrieNode_AddSequence(TrieNode * trie_node, 
                     uint8_t * sequence, 
                     uint32_t sequence_size, 
                     uint8_t *alphabet_size, 
                     uint8_t * charmap) {
    if (sequence_size == 0) {
        return 0;
    }

    if (TrieNode_IS_TERMINAL(trie_node)) {
        uint32_t suffix_size = TrieNode_GET_SUFFIX_SIZE(trie_node);
        if (sequence_size == suffix_size) {
            if (memcmp(TrieNode_GET_SUFFIX(trie_node), sequence, sequence_size) == 0){
                trie_node->count += 1;
                return 0;
            }
        }
        uint8_t * suffix = TrieNode_GET_SUFFIX(trie_node);
        uint8_t * tmp = PyMem_Malloc(suffix_size);
        if (tmp == NULL) {
            PyErr_NoMemory();
            return -1;
        }
        memcpy(tmp, suffix, suffix_size);
        TrieNode_Resize(trie_node, alphabet_size);
        int ret = TrieNode_AddSequence(trie_node, tmp, suffix_size, alphabet_size, charmap);
        PyMem_Free(tmp);
        if (ret != 0){
            return ret;
        }
    }

    uint8_t character = sequence[0];
    uint8_t node_character = charmap[character];
    if (node_character == 255) {
        charmap[character] = *alphabet_size;
        *alphabet_size = *alphabet_size + 1;
    }
    if (node_character >= trie_node->alphabet_size) {
        TrieNode_Resize(trie_node, node_character + 1);
    }
    TrieNode * next_node = trie_node->children[node_character];
    if (next_node == NULL) {
        next_node = TrieNode_NewLeaf(sequence + 1, sequence_size -1);
        if (next_node == NULL) {
            return -1;
        }
        trie_node->children[node_character] = next_node;
        return 0;
    }
    return TrieNode_AddSequence(next_node, sequence + 1, sequence_size - 1, alphabet_size, charmap);
}