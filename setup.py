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

from pathlib import Path

from setuptools import Extension, find_packages, setup

LONG_DESCRIPTION = Path("README.rst").read_text()

setup(
    name="fastqdedup",
    version="0.1.0-dev",
    description="Alignment-free UMI deduplication",
    author="Leiden University Medical Center",
    author_email="r.h.p.vorderman@lumc.nl",
    long_description=LONG_DESCRIPTION,
    long_description_content_type="text/x-rst",
    license="MIT",
    keywords="FASTQ UMI",
    zip_safe=False,
    packages=find_packages('src'),
    package_dir={'': 'src'},
    package_data={'fastqdedup': ['py.typed', '*.pyi', '*.c']},
    license_file="LICENSE",
    url="https://github.com/rhpvorderman/fastqdedup",
    classifiers=[
        "Programming Language :: Python :: 3 :: Only",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: C",
        "License :: OSI Approved :: MIT License",
    ],
    python_requires=">=3.6",
    install_requires=["dnaio >=0.8.0"],
    ext_modules=[
        Extension("fastqdedup._trie", ["src/fastqdedup/_triemodule.c"])
    ],
    entry_points={"console_scripts": [
        "fastqdedup = fastqdedup:main"]}
)
