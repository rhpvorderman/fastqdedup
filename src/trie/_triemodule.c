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
    (char *)n->children)

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

static uint8_t
TrieNode_AddSequence(TrieNode * trie_node, 
                     uint8_t * sequence, 
                     size_t sequence_size, 
                     uint8_t *alphabet_size, 
                     uint8_t * charmap) {
    if (sequence_size == 0) {
        return 0;
    }
    uint8_t character = sequence[0];
    uint8_t node_character = charmap[character];
    if (node_character == 255) {
        charmap[character] = *alphabet_size;
        *alphabet_size = *alphabet_size + 1;
    }


    if (TrieNode_IS_TERMINAL(trie_node)) {
        uint32_t suffix_size = TrieNode_GET_SUFFIX_SIZE(trie_node);

    }
}