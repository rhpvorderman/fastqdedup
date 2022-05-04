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

from fastqdedup import Trie

import pytest


def test_trie_one_seq():
    trie = Trie()
    trie.add_sequence("GATTACA")
    assert trie.contains_sequence("GATTACA", 0)
    assert trie.contains_sequence("AATTACA", 1)
    assert trie.contains_sequence("GATTACC", 1)
    assert trie.contains_sequence("GACCACA", 2)
    assert not trie.contains_sequence("GACCACA", 1)
    assert not trie.contains_sequence("GATTACC", 0)


def test_trie_subseq():
    trie = Trie()
    trie.add_sequence("GATTACA")
    trie.add_sequence("GATTA")
    assert trie.contains_sequence("GATTA")
    assert trie.contains_sequence("GATTACA")
    assert not trie.contains_sequence("GATTAC")


def test_trie_pop_cluster():
    trie = Trie()
    trie.add_sequence("AAAA")
    trie.add_sequence("AAAA")
    trie.add_sequence("AAAC")
    trie.add_sequence("AAGC")
    trie.add_sequence("AGGC")
    trie.add_sequence("CCCG")
    trie.add_sequence("CCCG")
    trie.add_sequence("TTCA")
    trie.add_sequence("TTCC")
    trie.add_sequence("TTTA")
    trie.add_sequence("TTT")
    trie.add_sequence("TTC")
    cluster_list = []
    while True:
        try:
            cluster = trie.pop_cluster(1)
        except LookupError:
            break
        cluster_list.append(cluster)
    cluster_set = [set(cluster) for cluster in cluster_list]
    expected_clusters = [
        {(2, "AAAA"), (1, "AAGC"), (1, "AAAC"), (1, "AGGC")},
        {(2, "CCCG")},
        {(1, "TTCA"), (1, "TTCC"), (1, "TTTA")},
        {(1, "TTT"), (1, "TTC")},  # Hamming distance only for equal size.
    ]
    for expected_cluster in expected_clusters:
        assert expected_cluster in cluster_set
        cluster_set.remove(expected_cluster)
    assert not cluster_set


def test_trie_new_with_alphabet():
    trie = Trie(alphabet="acd")
    assert trie.alphabet == "acd"


def test_trie_alphabet_repeated():
    with pytest.raises(ValueError) as error:
        Trie(alphabet="abcc")
    error.match("c was repeated")


def test_trie_alphabet_during_adding():
    trie = Trie()
    trie.add_sequence("abc")
    # No alphabet yet. The above simply creates a terminal node
    trie.add_sequence("badabccdaafacb")
    # The trie now has two branches one with a and one with b.
    assert trie.alphabet == "ab"
    trie.add_sequence("bcadac")
    assert trie.alphabet == "abc"


def test_trie_number_of_sequences():
    trie = Trie()
    trie.add_sequence("abc")
    trie.add_sequence("ab")
    trie.add_sequence("abcd")
    assert trie.number_of_sequences == 3
    while True:
        try:
            trie.pop_cluster(0)
        except LookupError:
            break
    assert trie.number_of_sequences == 0
