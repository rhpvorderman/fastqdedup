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

from fastqdedup import (
    cluster_dissection_adjacency,
    cluster_dissection_directional,
    cluster_dissection_most_reads,
    length_string_to_slices,
)

import pytest


@pytest.mark.parametrize(["string", "result"], [
    ("5,6,7", [slice(5), slice(6), slice(7)]),
    ("5:8,3,-5:3:-1", [slice(5, 8), slice(3), slice(-5, 3, -1)]),
    ("None:None:16", [slice(None, None, 16)]),
    ("::16", [slice(None, None, 16)])
])
def test_length_string_to_slices(string, result):
    assert length_string_to_slices(string) == result


class TestClusterDissection:
    TEST_CLUSTER = [
        (3, "AAAGT"),   # Derived
        (10, "AAAAT"),  # Derived
        (50, "AACAA"),  # Origin read
        (60, "AAAAA"),  # Origin read
        (10, "CAAAA"),  # Derived
        (30, "CTAAA")   # Origin read
    ]

    def test_most_reads(self):
        dissected = list(cluster_dissection_most_reads(self.TEST_CLUSTER))
        assert len(dissected) == 1
        assert dissected[0] == "AAAAA"

    def test_adjacency(self):
        dissected = set(cluster_dissection_adjacency(self.TEST_CLUSTER))
        assert len(dissected) == 3
        assert dissected == {"AAAAA", "CTAAA", "AAAGT"}

    def test_directional(self):
        dissected = set(cluster_dissection_directional(self.TEST_CLUSTER))
        assert len(dissected) == 3
        assert dissected == {"AACAA", "AAAAA", "CTAAA"}

    @pytest.mark.parametrize("function", [cluster_dissection_directional,
                                          cluster_dissection_adjacency,
                                          cluster_dissection_most_reads])
    def test_no_list_aliasing(self, function):
        cluster = self.TEST_CLUSTER[:]
        old_cluster = cluster[:]
        list(function(cluster))
        assert old_cluster == cluster
