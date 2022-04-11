from typing import Dict


class FastqTrieNode:
    children: Dict
    count: int
    suffix: str

    __slots__ = ["count", "suffix", "children"]

    def is_terminal(self):
        return bool(self.children)

    def __init__(self):
        self.children: Dict[FastqTrieNode] = {}
        self.count = 0
        self.suffix = ""

    def add_sequence(self, sequence: str):
        next_node = self
        for i, char in enumerate(sequence):
            node = next_node
            node.count += 1
            try:
                next_node = node.children[char]
            except KeyError:
                next_node = FastqTrieNode()
                next_node.count = 1
                next_node.suffix = sequence[i + 1:]
                next_node.suffix = ""
                node.children[char] = next_node
                break
            if next_node.is_terminal():
                end_sequence = sequence[i + 1:]
                if next_node.suffix == end_sequence:
                    next_node.count += 1
                    break
                next_node.add_sequence(next_node.suffix)
