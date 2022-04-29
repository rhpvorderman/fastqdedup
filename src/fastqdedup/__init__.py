# Copyright (C) 2022 Leiden University Medical Center
# This file is part of fastqdedup
#
# fastqdedup is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.
#
# fastqdedup is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with fastqdedup.  If not, see <https://www.gnu.org/licenses/

import argparse
import contextlib
import functools
from typing import List, Iterable, Iterator, Optional, Tuple

import dnaio

import xopen

from ._trie import Trie

DEFAULT_PREFIX = "fastqdedup_R"
DEFAULT_MAX_DISTANCE = 1


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


def deduplicate_naive(input_files: List[str],
                      output_files: List[str],
                      check_lengths: Optional[List[int]],
                      max_distance: int = DEFAULT_MAX_DISTANCE):
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
        if not trie.contains_sequence(key, max_distance):
            for record, output in zip(records, output_files):
                output.write(record.fastq_bytes())
        # Always add sequence to the trie. This way clusters use only the first
        # read that showed up.
        trie.add_sequence(key)


def deduplicate_cluster(input_files: List[str],
                        output_files: List[str],
                        check_lengths: Optional[List[int]],
                        max_distance: int = DEFAULT_MAX_DISTANCE):
    if len(input_files) != len(output_files):
        raise ValueError(f"Amount of output files ({len(output_files)}) "
                         f"must be equal to the amount of input files "
                         f"({len(input_files)}). ")
    if check_lengths and len(input_files) != len(check_lengths):
        raise ValueError(f"Amount of check lengths ({len(check_lengths)}) "
                         f"must be equal to the amount of input files "
                         f"({len(input_files)}). ")
    input_readers = [file_to_fastq_reader(f) for f in input_files]
    trie = Trie()
    for records in zip(*input_readers):  # type: Tuple[dnaio.SequenceRecord, ...]
        key = _key_from_records(records, check_lengths)
        trie.add_sequence(key)
    deduplicated_set = set()
    while True:
        try:
            cluster = trie.pop_cluster()
        except LookupError:
            break
        if len(cluster) > 1:
            cluster.sort()
        # Hash the key first before storing in the set to save memory.
        deduplicated_set.add(hash(cluster[0][0]))

    # Read the fastq files again and filter against the deduplicated set
    input_readers = [file_to_fastq_reader(f) for f in input_files]
    output_stack = contextlib.ExitStack()
    output_opener = functools.partial(xopen.xopen, mode="wb",
                                      compresslevel=1, threads=0)
    output_files: List[xopen.PipedCompressionWriter] = [
        output_stack.enter_context(output_opener(x)) for x in output_files]
    for records in zip(*input_readers):  # type: Tuple[dnaio.SequenceRecord, ...]
        key = _key_from_records(records, check_lengths)
        if hash(key) in deduplicated_set:
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
             "``fastqdedup --check-lengths 16,8 R1.fastq R2.fastq`` only "
             "checks the first 16 bases of R1 and the first 8 bases of R2 for "
             "duplication.")
    parser.add_argument(
        "-o", "--output", action="append", required=False,
        help="Output file (optional), must be specified multiple times for "
             "multiple input files. For example ``fastqdedup -o dedupR1.fastq "
             "-o dedupR2.fastq R1.fastq R2.fastq``.")
    parser.add_argument("-p", "--prefix", default=DEFAULT_PREFIX,
                        help=f"Prefix for the output files. "
                             f"Default: '{DEFAULT_PREFIX}'")
    parser.add_argument("-d", "--max-distance", type=int,
                        default=DEFAULT_MAX_DISTANCE,
                        help="The Hamming distance at which inputs are "
                             f"considered different. "
                             f"Default: {DEFAULT_MAX_DISTANCE}.")
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
    deduplicate_cluster(input_files, output_files, check_lengths, max_distance)


if __name__ == "__main__":
    main()
