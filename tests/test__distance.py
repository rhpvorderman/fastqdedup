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

from fastqdedup._distance import within_distance

import pytest


@pytest.mark.parametrize(
    ["string1", "string2", "max_distance", "result"],
    [
        ("AAAA", "AAAA", 0, True),
        ("AAAA", "AAA", 3, False),
        ("AAAA", "AAAC", 1, True),
        ("AAAA", "AAAC", 0, False),
        ("AACA", "AAAC", 2, True),
        ("AACC", "CCAA", 3, False)
    ]
)
def test_within_distance_hamming(string1, string2, max_distance, result):
    assert within_distance(string1, string2, max_distance) is result


@pytest.mark.parametrize(
    ["string1", "string2", "max_distance", "result"],
    [
        ("AAAA", "AAAA", 0, True),
        ("AAAA", "AAA", 1, True),
        ("AAAA", "A", 3, True),
        ("AAA", "C", 2, False),  # One substitution 2 deletions
        ("AAA", "C", 3, True),
        ("AAAA", "AAAC", 1, True),
        ("AAAA", "AAAC", 0, False),
        ("AACA", "AAAC", 2, True),
        ("AACC", "CCAA", 3, False)
    ]
)
def test_within_distance_levenstein(string1, string2, max_distance, result):
    assert within_distance(string1, string2, max_distance,
                           use_edit_distance=True) is result
