from pathlib import Path

from setuptools import Extension, find_packages, setup

LONG_DESCRIPTION = Path("README.rst").read_text()

setup(
    name="umidedup",
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
    package_data={'umidedup': ['py.typed', '*.pyi']},
    license_file="LICENSE",
    url="https://github.com/rhpvorderman/umidedup",
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
        Extension("umidedup._trie", ["src/umidedup/_triemodule.c"])
    ],
    entry_points={"console_scripts": [
        "umidedup = umidedup:main"]}
)
