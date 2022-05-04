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

from typing import List, Tuple


class Trie:
    def __init__(self, alphabet: str = ""): ...

    def add_sequence(self, sequence: str): ...

    def contains_sequence(self, 
                          sequence: str,
                          max_hamming_distance: int = 0
                          ) -> bool: ...

    def pop_cluster(self, max_hamming_distance: int) -> List[Tuple[int, str]]: ...

    def memory_size(self) -> int: ...

    def raw_stats(self) -> List[List[int]]: ...

    @property
    def alphabet(self) -> str: ...

