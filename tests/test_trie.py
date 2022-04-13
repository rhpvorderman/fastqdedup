from umidedup._trie import Trie


def test_trie_one_seq():
    trie = Trie()
    trie.add_sequence("GATTACA")
    assert trie.sequence_present_hamming("GATTACA", 0)
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
