import argparse
import contextlib
import functools
from typing import List, Iterable, Iterator, Optional, Tuple

import dnaio

import xopen

from ._trie import Trie

DEFAULT_PREFIX="umidedup_R"


def file_to_fastq_reader(filename: str) -> Iterator[dnaio.SequenceRecord]:
    opener = functools.partial(xopen.xopen, threads=0)
    with dnaio.open(filename, mode="r", opener=opener) as fastqreader:
        yield from fastqreader


def _key_from_records(records: Iterable[dnaio.SequenceRecord],
                      check_lengths: Optional[Iterable[int]]):
    """
    Creates a key from several records. This is an internal function.
    This function does not check if records and check_lengths have the same
    length.
    """
    if check_lengths:
        return "".join(record.sequence[:length] for record, length in zip(records, check_lengths))
    return "".join(record.sequence for record in records)


def deduplicate(input_files: List[str],
                output_files: List[str],
                check_lengths: Optional[List[int]],
                max_distance: int = 0):
    if len(input_files) != len(output_files):
        raise ValueError(f"Amount of output files ({len(output_files)}) "
                         f"must be equal to the amount of input files "
                         f"({len(input_files)}). ")
    if check_lengths and len(input_files) != len(check_lengths):
        raise ValueError(f"Amount of check lengths ({len(check_lengths)}) "
                         f"must be equal to the amount of input files "
                         f"({len(input_files)}). ")
    input_readers = [file_to_fastq_reader(f) for f in input_files]
    output_stack = contextlib.ExitStack()
    output_opener = functools.partial(xopen.xopen, mode="wb",
                                      compresslevel=1, threads=0)
    output_files: List[xopen.PipedCompressionWriter] = [
        output_stack.enter_context(output_opener(x)) for x in output_files]
    trie = Trie()
    for records in zip(*input_readers):  # type: Tuple[dnaio.SequenceRecord, ...]
        key = _key_from_records(records, check_lengths)
        if trie.sequence_present_hamming(key, max_distance):
            continue
        trie.add_sequence(key)
        for record, output in zip(records, output_files):
            output.write(record.fastq_bytes())


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
    parser.add_argument(
        "-o", "--output", action="append", required=False,
        help="Output file, must be specified multiple times for multiple "
             "input files. For example ``umidedup -o dedupR1.fastq "
             "-o dedupR2.fastq R1.fastq R2.fastq``.")
    parser.add_argument("-p", "--prefix", default=DEFAULT_PREFIX,
                        help=f"Prefix for the output files. "
                             f"Default: '{DEFAULT_PREFIX}'")
    parser.add_argument("-d", "--max-distance", type=int, default=0,
                        help="The Hamming distance at which inputs are "
                             "considered different.")
    return parser


def main():
    args = argument_parser().parse_args()
    input_files: List[str] = args.fastq
    if args.check_lengths:
        check_lengths = [int(x) for x in args.check_lengths.split(",")]
    else:
        check_lengths = None
    if args.output:
        output_files = args.output
    else:
        output_files = [args.prefix + str(x) + ".fastq.gz"
                        for x in range(1, len(input_files) + 1)]
    max_distance = args.max_distance
    deduplicate(input_files, output_files, check_lengths, max_distance)


if __name__ == "__main__":
    main()
