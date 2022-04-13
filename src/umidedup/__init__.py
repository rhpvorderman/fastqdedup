import argparse
import sys
from typing import List

import dnaio

from ._trie import Trie


def argument_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "fastq", metavar="FASTQ", nargs="+",
        help="Forward FASTQ and optional reverse and UMI FASTQ files."
    )
    parser.add_argument(
        "--check-lengths",
        help="Comma-separated string with the maximum string check length of "
             "each file. For example "
             "``umidedup --check-lengths 16,8 R1.fastq R2.fastq`` only checks "
             "the first 16 bases of R1 and the first 8 bases of fastq for "
             "duplication.")
    return parser


def main():
    trie = Trie()
    with dnaio.open(sys.argv[1], mode="r") as fastq_reader:
        for record in fastq_reader:
            trie.add_sequence(record.sequence)


if __name__ == "__main__":
    main()
