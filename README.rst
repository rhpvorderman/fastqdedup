fastqdedup
==========

Fastqdedup deduplicates FASTQ files with UMI data while taking sequencing
errors into account.

Fastqdedup can operate on one or multiple FASTQ files. So single reads with
prepended UMIs, paired reads with prepended UMIs or paired reads with separate
UMI file are all possible.

It executes the following steps:

+ Read all sequences
+ Determine clusters based on hamming distance
+ Take the most frequently read sequence from each cluster. Put these in a set.
+ Read all sequence records again from the FASTQ files. Check against the set.
  If a sequence is present the sequence record is output to the output FASTQ
  files and the sequence is removed from the set so that it is not included
  again.

Installation
------------
To install the latest development version:

``pip install git+https://github.com/rhpvorderman/fastqdedup.git``

Usage
-----

Since all the sequences are stored in memory the memory usage is between
approximately 45-110% of the size of the original uncompressed FASTQ files
depending on read length and similarities between reads.
The memory usage can be substantially reduced by setting ``--check-lengths``.

.. code-block::

    usage: fastqdedup [-h] [-l CHECK_LENGTHS] [-o OUTPUT] [-p PREFIX]
                      [-d MAX_DISTANCE] [-e MAX_AVERAGE_ERROR_RATE] [-E] [--edit]
                      [-c {highest_count,adjacency,directional}] [-v] [-q]
                      FASTQ [FASTQ ...]

    positional arguments:
      FASTQ                 Forward FASTQ and optional reverse and UMI FASTQ
                            files.

    optional arguments:
      -h, --help            show this help message and exit
      -l CHECK_LENGTHS, --check-lengths CHECK_LENGTHS
                            Comma-separated string with the maximum string check
                            length of each file. For example 'fastqdedup --check-
                            lengths 16,8 R1.fastq R2.fastq' only checks the first
                            16 bases of R1 and the first 8 bases of R2 for
                            duplication. Supports slice notation such as '4:8' or
                            '::8'.
      -o OUTPUT, --output OUTPUT
                            Output file (optional), must be specified multiple
                            times for multiple input files. For example
                            ``fastqdedup -o dedupR1.fastq -o dedupR2.fastq
                            R1.fastq R2.fastq``.
      -p PREFIX, --prefix PREFIX
                            Prefix for the output files. Default: 'fastqdedup_R'
      -d MAX_DISTANCE, --max-distance MAX_DISTANCE
                            The Hamming distance at which inputs are considered
                            different. Default: 1.
      -e MAX_AVERAGE_ERROR_RATE, --max-average-error-rate MAX_AVERAGE_ERROR_RATE
                            The maximum average per base error rate for each FASTQ
                            record. Average is evaluated over bases taken into
                            account by --check-lengths.Default: 0.001
      -E, --no-average-error-rate-filter
                            Do not filter on average per base error rate.
      --edit                Use edit (Levenshtein) distance instead of Hamming
                            distance.
      -c {highest_count,adjacency,directional}, --cluster-dissection-method {highest_count,adjacency,directional}
                            How to approach clusters with multiple reads.
                            'highest_count' selects only one read, the one with
                            the highest count. 'adjacency' starts from the read
                            with the highest count and selects all reads that are
                            within the specified distance. The process is repeated
                            for the remaining reads. 'directional' is similar to
                            adjacency but uses counts to determine if an error is
                            a PCR/sequencing artifact or derived from a difference
                            in the molecule (default).
      -v, --verbose         Increase log verbosity.
      -q, --quiet           Reduce log verbosity.

Methodology
-----------
Fastqdedup was heavily inspired by `@jfjlaros's trie implementation
<https://github.com/jfjlaros/trie>`_. A `trie
<https://en.wikipedia.org/wiki/Trie>`_ is a way to store sequences in a tree
format. A disadvantage of this format is that memory usage is quite high since
every character is stored as a node in the tree, and a node uses a lot more
space than a single character.

Fastqdedup utilizes the following optimizations.

1. `Radix tree <https://en.wikipedia.org/wiki/Radix_tree>`_ suffixes.
   While a lot of similarities between FASTQ files are
   basically guaranteed across the first 8 characters (65536 possible sequences)
   the amount of possible sequences explodes after that.
   Therefore a trie of FASTQ sequences contains a lot of branches that consists
   of nodes that only link one parent to one child for quite a long length.
   Fastqdedup stores terminating branches as strings rather than nodes, saving
   a lot of memory.
2. Variable size nodes. The trie is initialized with an ``ACGTN`` alphabet, where
   ``A`` has index ``0``, ``C`` has index ``1`` etc. Nodes are sized
   such that they can contain the character with the highest index that is in
   the node.
3. Fast fail Hamming distances. Since the most likely result is that the
   sequence is not seen before, the algorithm gives up as quickly as it finds
   the Hamming distance will be exceeded.

Background
----------
Illumina reads have an error rate of approximately 1 in 1000. This means that
the chance of a sequencing error in a 150bp illumina read is
``1 - ((1 - (1/1000)) ^ 150) ~= 14%``. This means that most (86%) illumina
reads are correct and that it is hard to distinguish between reads with one
SNP and reads with one sequencing error. It is also hard to distinguish between
sequencing replicates and actual biological replicates.

Unique Molecular Identifiers (UMIs) solve this problem by prepending each read
with a short sequence of bases. With a UMI length of 8, there are 65536
(``4^8``) possible UMIs. The chance of two biological replicates having the
same UMI is 1 in 65536. Therefore two sequences with the same UMI and only one
base pair different are probably derived from the same biological replicate
since there is a 14% chance of a sequencing error.

