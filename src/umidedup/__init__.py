from ._trie import Trie

import dnaio

import sys


def main():
    trie = Trie()
    with dnaio.open(sys.argv[1], mode="r") as fastq_reader:
        for record in fastq_reader:
            trie.add_sequence(record.sequence)
    input("Finished?")
    del(trie)
    input("Checked mem usage?")

if __name__ == "__main__":
    main()
