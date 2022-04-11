from typing import List


class FastqTrieNode:
    children: List[...]

    __slots__ = ["children"]

    def __init__(self, size: int):
        children = [None] * size

    def is_leaf(self):
        return False

class FastqTrieLeaf:
    count: int
    suffix: str

    __slots__ = ["count", "suffix"]

    def __init__(self, suffix):
        self.count = 1
        self.suffix = suffix

    def is_leaf(self):
        return True

class FastqTrie(FastqTrieNode):
    __slots__ = ["charmap", "alphabet_size"]

    def __init__(self):
        self.charmap = {}
        self.alphabet_size = 0
        super().__init__(self.alphabet_size)

    def add_sequence(self, sequence: str):
        node = self
        for i, char in enumerate(sequence):
            try:
                index = self.charmap[char]
            except KeyError:
                self.charmap[char] = self.alphabet_size
                index = self.alphabet_size
                self.alphabet_size += 1
            try:
                next_node = node.children[index]
            except IndexError:
                node.children.extend(
                    [None] * (self.alphabet_size - len(node.children)))
                next_node = None
            if next_node is None:
                node.children[index] = FastqTrieLeaf(sequence[i + 1:])
            if isinstance(next_node, FastqTrieLeaf):
                if next_node.suffix == sequence[i + 1:]:
                    next_node.count += 1
                else:


            elif isinstance(next_node, FastqTrieNode):