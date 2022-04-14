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


def test_trie_one_seq():
    trie = Trie()
    trie.add_sequence("GATTACA")
    assert trie.sequence_present_hamming("GATTACA", 0)
    assert trie.sequence_present_hamming("AATTACA", 1)
    assert trie.sequence_present_hamming("GATTACC", 1)
    assert trie.sequence_present_hamming("GACCACA", 2)
    assert not trie.sequence_present_hamming("GACCACA", 1)
    assert not trie.sequence_present_hamming("GATTACC", 0)


def test_trie_subseq():
    trie = Trie()
    trie.add_sequence("GATTACA")
    trie.add_sequence("GATTA")
    assert trie.sequence_present_hamming("GATTA")
    assert trie.sequence_present_hamming("GATTACA")
    assert not trie.sequence_present_hamming("GATTAC")
