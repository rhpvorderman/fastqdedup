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

from fastqdedup import length_string_to_slices

import pytest


@pytest.mark.parametrize(["string", "result"], [
    ("5,6,7", [slice(5), slice(6), slice(7)]),
    ("5:8,3,-5:3:-1", [slice(5, 8), slice(3), slice(-5, 3, -1)]),
    ("None:None:16", [slice(None, None, 16)]),
    ("::16", [slice(None, None, 16)])
])
def test_length_string_to_slices(string, result):
    assert length_string_to_slices(string) == result
