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
import datetime
import functools
import io
import logging
import resource
import time
from typing import (Any, Callable, Dict, IO, Iterable, Iterator, List,
                    Optional, Set, Tuple)

import dnaio

import xopen

from ._distance import hamming_distance
from ._fastq import average_error_rate as fastq_average_error_rate
from ._trie import Trie

DEFAULT_PREFIX = "fastqdedup_R"
DEFAULT_MAX_DISTANCE = 1
DEFAULT_CLUSTER_DISSECTION = "directional"
DEFAULT_MAX_AVERAGE_ERROR_RATE = 0.001


class Timer:
    """Simple timer object to reduce timing boilerplate"""
    def __init__(self):
        self.start_time = time.time()

    def get_difference(self) -> datetime.timedelta:
        current_time = time.time()
        delta = datetime.timedelta(seconds=round(current_time - self.start_time))
        self.start_time = current_time
        return delta


def file_to_fastq_reader(filename: str) -> Iterator[dnaio.SequenceRecord]:
    opener = functools.partial(xopen.xopen, threads=0)
    with dnaio.open(filename, mode="r", opener=opener) as fastqreader:  # type: ignore
        yield from fastqreader


def cluster_dissection_directional(cluster: List[Tuple[int, str]],
                                   max_distance: int = DEFAULT_MAX_DISTANCE
                                   ) -> Iterator[str]:
    """Take the read with the highest count. Test all other reads, if they are
    within hamming distance and have a count for which 2n-1 is lower than the
    template count, assume they are derived trough PCR artifact. In that case
    add them to the template chain for testing."""
    cluster = sorted(cluster)
    while cluster:
        # The last read has the highest count since we sorted (ascending order
        # is default).
        original_item = cluster.pop()
        _, original_string = original_item
        template_list = [original_item]

        # Appending to the list while iterating over it is safe. Each __next__
        # call retrieves an updated length value from the list.
        for template_count, template_string in template_list:
            if not cluster:
                break  # fast exit when nothing to compare
            distinct_list = []
            for item in cluster:
                compare_count, compare_string = item
                if (2 * compare_count - 1) <= template_count:
                    if hamming_distance(template_string, compare_string
                                        ) <= max_distance:
                        template_list.append(item)
                        continue
                distinct_list.append(item)
            cluster = distinct_list
        yield original_string


def cluster_dissection_highest_count(cluster: List[Tuple[int, str]],
                                     max_distance: int = DEFAULT_MAX_DISTANCE
                                     ) -> Iterator[str]:
    """Select the read with the highest count. Only yields 1 read."""
    cluster = sorted(cluster, reverse=True)
    # We sorted in descending order, so the first read has the highest count.
    _, string = cluster[0]
    yield string


def cluster_dissection_adjacency(cluster: List[Tuple[int, str]],
                                 max_distance: int = DEFAULT_MAX_DISTANCE
                                 ) -> Iterator[str]:
    """Take the read with the highest count, find all the reads are not
    directly adjacent within max distance and repeat."""
    cluster = sorted(cluster, reverse=True)
    while cluster:
        # We sorted in descending order, so the first read has the highest count.
        _, template_string = cluster[0]
        distinct_list = []
        for item in cluster[1:]:
            _, compare_string = item
            if hamming_distance(template_string, compare_string) > max_distance:
                distinct_list.append(item)
        yield template_string
        cluster = distinct_list[:]


ClusterDissectionFunc = Callable[[List[Tuple[int, str]], int], Iterator[str]]
CLUSTER_DISSECTION_METHODS: Dict[str, ClusterDissectionFunc] = {
    "highest_count": cluster_dissection_highest_count,
    "adjacency": cluster_dissection_adjacency,
    "directional": cluster_dissection_directional,
}


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


def joinfunc_from_check_slices(
        check_slices: Iterable[slice]
) -> Callable[[Iterable[str]], str]:
    def joinfunc(strings: Iterable[str]):
        return "".join(string[slc]
                       for string, slc
                       in zip(strings, check_slices))
    return joinfunc


def fastq_files_to_records(
        input_files: List[str]
) -> Iterator[Tuple[dnaio.SequenceRecord, ...]]:
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
        yield records


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


def deduplicate_cluster(
    input_files: List[str],
    output_files: List[str],
    check_slices: Optional[List[slice]],
    max_distance: int = DEFAULT_MAX_DISTANCE,
    max_average_error_rate: float = DEFAULT_MAX_AVERAGE_ERROR_RATE,
    cluster_dissection_func: ClusterDissectionFunc = cluster_dissection_directional,
):
    if len(input_files) != len(output_files):
        raise ValueError(f"Amount of output files ({len(output_files)}) "
                         f"must be equal to the amount of input files "
                         f"({len(input_files)}). ")
    if check_slices and len(input_files) != len(check_slices):
        raise ValueError(f"Amount of check lengths ({len(check_slices)}) "
                         f"must be equal to the amount of input files "
                         f"({len(input_files)}). ")

    # Create a keyfunc in order to collapse multiple FASTQ records into
    # one key that can be used to determine the clusters.
    if check_slices:
        joinfunc = joinfunc_from_check_slices(check_slices)
    else:
        joinfunc = "".join

    record_tuples = fastq_files_to_records(input_files)
    filter_on_quality = max_average_error_rate < 1.0
    total_records = 0
    discarded_records = 0
    timer = Timer()
    logger = logging.getLogger("fastqdedup")
    trie = Trie(alphabet="ACGTN")

    for record_tuple in record_tuples:
        qualities = joinfunc(record.qualities
                             for record in record_tuple
                             if record.qualities is not None)
        total_records += 1
        if (filter_on_quality and
                fastq_average_error_rate(qualities) > max_average_error_rate):
            discarded_records += 1
            continue
        key = joinfunc(record.sequence for record in record_tuple)
        trie.add_sequence(key)
    if filter_on_quality:
        logger.info(
            f"{discarded_records} records out of {total_records} "
            f"records had an error rate higher than {max_average_error_rate} "
            f"and were discarded.")
    logger.info(f"Processed {trie.number_of_sequences} sequences. "
                f"({timer.get_difference()})")
    if logger.level <= logging.DEBUG:
        # Do not perform expensive stats calc when not requested.
        stats = trie_stats(trie)
        logger.debug(f"Calculated stats. ({timer.get_difference()})")
        logger.debug("\n" + stats)

    # Create a deduplicated set by popping of clusters from the trie and
    # selecting the most prevalent read per cluster.
    # Not the keys, but the hash values of the keys are stored in the set.
    # This saves a lot of memory.
    deduplicated_set: Set[int] = set()
    number_of_clusters = 0
    while trie.number_of_sequences:
        cluster = trie.pop_cluster(max_distance)
        number_of_clusters += 1
        for key in cluster_dissection_func(cluster, max_distance):
            deduplicated_set.add(hash(key))

    del(trie)
    logger.info(f"Found {len(deduplicated_set)} distinct reads "
                f"in {number_of_clusters} clusters."
                f"({timer.get_difference()})")

    def hashfunc(records):
        return hash(joinfunc(record.sequence for record in records))

    filter_fastq_files_on_set(input_files, output_files, deduplicated_set, hashfunc)
    logger.info(f"Filtered FASTQ files based on distinct reads from each cluster. "
                f"({timer.get_difference()}) ")


def initiate_logger(verbose: int = 0, quiet: int = 0):
    log_level = logging.INFO - 10 * (verbose - quiet)
    logger = logging.getLogger("fastqdedup")
    logger.setLevel(log_level)
    console_handler = logging.StreamHandler()
    console_handler.setLevel(log_level)
    formatter = logging.Formatter(
        "{asctime}:{levelname}:{name}: {message}",
        datefmt="%m/%d/%Y %I:%M:%S",
        style="{")
    console_handler.setFormatter(formatter)
    logger.addHandler(console_handler)


def argument_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "fastq", metavar="FASTQ", nargs="+",
        help="Forward FASTQ and optional reverse and UMI FASTQ files."
    )
    parser.add_argument(
        "-l", "--check-lengths",
        help="Comma-separated string with the maximum string check length of "
             "each file. For example "
             "'fastqdedup --check-lengths 16,8 R1.fastq R2.fastq' only "
             "checks the first 16 bases of R1 and the first 8 bases of R2 for "
             "duplication. Supports slice notation such as '4:8' or '::8'.")
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
    parser.add_argument("-e", "--max-average-error-rate", type=float,
                        default=DEFAULT_MAX_AVERAGE_ERROR_RATE,
                        help=f"The maximum average per base error rate for "
                             f"each FASTQ record. Average is evaluated over "
                             f"bases taken into account by --check-lengths."
                             f"Default: {DEFAULT_MAX_AVERAGE_ERROR_RATE}")
    parser.add_argument("-E", "--no-average-error-rate-filter",
                        action="store_const", dest="max_average_error_rate",
                        const=1.0,
                        help="Do not filter on average per base error rate.")
    parser.add_argument(
        "-c", "--cluster-dissection-method",
        choices=CLUSTER_DISSECTION_METHODS.keys(),
        default=DEFAULT_CLUSTER_DISSECTION,
        help="How to approach clusters with multiple reads. "
             "'highest_count' selects only one read, the one with the highest "
             "count. "
             "'adjacency' starts from the read with the highest count and "
             "selects all reads that are within the specified distance. "
             "The process is repeated for the remaining reads. "
             "'directional' is similar to adjacency but uses counts to "
             "determine if an error is a PCR/sequencing artifact or derived "
             "from a difference in the molecule (default).")
    parser.add_argument("-v", "--verbose", action="count", default=0,
                        help="Increase log verbosity.")
    parser.add_argument("-q", "--quiet", action="count", default=0,
                        help="Reduce log verbosity.")
    return parser


def length_string_to_slices(length_string: str) -> List[slice]:
    """
    Converts a comma-separated string of lengths or slices such as 8,8,8 or
    8:16,8,24:8:-1 to a list of slice objects.
    """
    parts = length_string.split(",")
    slices = []
    for part in parts:
        values = [None if (x == "None" or x == "") else int(x)
                  for x in part.split(":")]
        slices.append(slice(*values))
    return slices


def main():
    args = argument_parser().parse_args()
    initiate_logger(args.verbose, args.quiet)
    logger = logging.getLogger("fastqdedup")

    input_files: List[str] = args.fastq
    if args.check_lengths:
        check_slices = length_string_to_slices(args.check_lengths)
    else:
        check_slices = None
    if args.output:
        output_files = args.output
    else:
        output_files = [args.prefix + str(x) + ".fastq.gz"
                        for x in range(1, len(input_files) + 1)]
    max_distance = args.max_distance
    max_average_error_rate = args.max_average_error_rate
    cluster_dissection_func = CLUSTER_DISSECTION_METHODS[
        args.cluster_dissection_method]
    timer = Timer()
    logger.info(f"Input files: {', '.join(input_files)}")
    logger.info(f"Output files: {', '.join(output_files)}")
    logger.info(f"Check lengths: {args.check_lengths}")
    logger.info(f"Maximum hamming distance: {max_distance}")
    logger.info(f"Maximum average error rate: {max_average_error_rate}")
    logger.info(f"Cluster dissection method: {args.cluster_dissection_method}")
    deduplicate_cluster(input_files, output_files, check_slices, max_distance,
                        max_average_error_rate,
                        cluster_dissection_func)
    resources = resource.getrusage(resource.RUSAGE_SELF)
    logger.info(f"Finished. Total time: {timer.get_difference()}. "
                f"Memory usage: {resources.ru_maxrss / (1024 ** 2):.2} GiB")


if __name__ == "__main__":
    main()
