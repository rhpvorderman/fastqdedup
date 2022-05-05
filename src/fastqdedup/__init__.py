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
import io
from typing import Any, Callable, IO, Iterable, Iterator, List, Optional, Set, Tuple

import dnaio

import xopen

from ._trie import Trie

DEFAULT_PREFIX = "fastqdedup_R"
DEFAULT_MAX_DISTANCE = 1


def file_to_fastq_reader(filename: str) -> Iterator[dnaio.SequenceRecord]:
    opener = functools.partial(xopen.xopen, threads=0)
    with dnaio.open(filename, mode="r", opener=opener) as fastqreader:  # type: ignore
        yield from fastqreader


def trie_stats(trie: Trie) -> str:
    outbuffer = io.StringIO()
    raw_stats = trie.raw_stats()
    layer_size = len(trie.alphabet) + 1
    all_totals = [0 for _ in range(layer_size + 1)]
    outbuffer.write("layer     terminal  " +
                    "".join(f"{i:10}" for i in range(1, layer_size)) +
                    "     total\n")
    for i, layer_stats in enumerate(raw_stats):
        total = sum(layer_stats)
        for j in range(layer_size):
            all_totals[j] += layer_stats[j]
        all_totals[layer_size] += total
        line = [str(i)] + layer_stats + [total]  # type: ignore
        outbuffer.write("".join(f"{i:10}" for i in line) + "\n")
    last_line = ["total"] + all_totals  # type: ignore
    outbuffer.write("".join(f"{i:10}" for i in last_line) + "\n")
    node_memory_usage = sum((8 + 8 * i) * all_totals[i] for i in range(layer_size))
    total_memory_usage = trie.memory_size()
    suffix_memory_usage = total_memory_usage - node_memory_usage
    gb = 1024 ** 3
    outbuffer.write(f"Node memory usage: {node_memory_usage / gb:.2} GiB\n"
                    f"Suffix memory usage: {suffix_memory_usage / gb:.2} GiB\n"
                    f"Total memory usage: {total_memory_usage / gb:.2} GiB\n")
    return outbuffer.getvalue()


def _key_from_records(records: Iterable[dnaio.SequenceRecord],
                      check_lengths: Optional[Iterable[int]]):
    """
    Creates a key from several records. This is an internal function.
    This function does not check if records and check_lengths have the same
    length.
    """
    if check_lengths:
        return "".join(record.sequence[:length]
                       for record, length
                       in zip(records, check_lengths))
    return "".join(record.sequence for record in records)


def cluster_keys_hamming(keys: Iterable[str],
                         max_distance: int = DEFAULT_MAX_DISTANCE,
                         initial_alphabet: str = "ACGTN",
                         ) -> Iterator[str]:
    """
    Cluster the keys based on hamming distance. Return the most common key
    per cluster.
    :param keys: An iterable of keys that need to be clustered.
    :param max_distance: The maximal hamming distance
    :param initial_alphabet: The starting alphabet for the trie.
    :return: An iterator of most common key per cluster
    """
    trie = Trie(alphabet=initial_alphabet)
    for key in keys:
        trie.add_sequence(key)
    while trie.number_of_sequences:
        cluster = trie.pop_cluster(max_distance)
        cluster.sort(reverse=True)
        count, key = cluster[0]
        yield key


def fastq_files_to_keys(
        input_files: List[str],
        keyfunc: Callable[[Iterable[dnaio.SequenceRecord]], str]
) -> Iterator[str]:
    """

    :param input_files:
    :param keyfunc:
    :return:
    """
    input_readers = [file_to_fastq_reader(f) for f in input_files]
    for records in zip(*input_readers):  # type: Tuple[dnaio.SequenceRecord, ...]
        first_record = records[0]
        for record in records[1:]:
            if not first_record.is_mate(record):
                raise dnaio.FastqFormatError(
                    f"FASTQ files not in sync: {first_record.name} is not a "
                    f"mate of {record.name}", None)
        yield keyfunc(records)


def filter_fastq_files_on_set(
        input_files: List[str],
        output_files: List[str],
        filter_set: Set[Any],
        keyfunc: Callable[[Iterable[dnaio.SequenceRecord]], Any]
):
    input_readers = [file_to_fastq_reader(f) for f in input_files]
    output_stack = contextlib.ExitStack()
    output_opener = functools.partial(xopen.xopen, mode="wb",
                                      compresslevel=1, threads=0)
    output_writers: List[IO[Any]] = [
        output_stack.enter_context(output_opener(x)) for x in output_files]
    for records in zip(*input_readers):
        key = keyfunc(records)
        if key in filter_set:
            filter_set.remove(key)
            for output, record in zip(output_writers, records):
                output.write(record.fastq_bytes())


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
    trie = Trie(alphabet="ACGTN")
    for records in zip(*input_readers):  # type: Tuple[dnaio.SequenceRecord, ...]
        key = _key_from_records(records, check_lengths)
        trie.add_sequence(key)
    deduplicated_set = set()
    while trie.number_of_sequences:
        cluster = trie.pop_cluster(max_distance)
        if len(cluster) > 1:
            # Reverse sort so read with highest count is first.
            cluster.sort(reverse=True)
        count, key = cluster[0]
        # Hash the key first before storing in the set to save memory.
        deduplicated_set.add(hash(key))
    del(trie)
    # Read the fastq files again and filter against the deduplicated set
    input_readers = [file_to_fastq_reader(f) for f in input_files]
    output_stack = contextlib.ExitStack()
    output_opener = functools.partial(xopen.xopen, mode="wb",
                                      compresslevel=1, threads=0)
    outputs: List[IO[Any]] = [
        output_stack.enter_context(output_opener(x)) for x in output_files]
    for records in zip(*input_readers):
        key = _key_from_records(records, check_lengths)
        h = hash(key)
        if h in deduplicated_set:
            deduplicated_set.remove(h)
            for record, output in zip(records, outputs):
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
