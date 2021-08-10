""" MetaProFi setup module
"""

from setuptools import setup, find_packages
from Cython.Build import cythonize
import numpy as np
from metaprofi.version import __version__

with open("README.md", "r") as readme:
    long_description = readme.read()

setup(
    name="MetaProFi",
    version=__version__,
    license="GPLv2+",
    author="Sanjay Kumar Srikakulam",
    maintainer="Sanjay Kumar Srikakulam",
    description="MetaProFi is a protein-based Bloom filter tool for storing and querying sequence data for accurate identification of functionally relevant genetic variants",
    long_description=long_description,
    long_description_content_type="text/markdown",
    keywords="metagenome protein finder bloom filters index sequence search",
    classifiers=[
        "Development Status :: 3 - Alpha",
        "License :: OSI Approved :: GNU General Public License v2 or later (GPLv2+)",
        "Programming Language :: Python :: 3.7",
        "Operating System :: POSIX :: Linux",
        "Natural Language :: English",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    setup_requires=["setuptools_scm", "Cython==0.29.23", "numpy==1.20.3"],
    include_package_data=True,
    python_requires=">=3.7",
    packages=find_packages(),
    install_requires=[
        "numpy==1.20.3",
        "Cython==0.29.23",
        "zarr==2.8.3",
        "pyfastx==0.8.4",
        "bitarray==2.2.1",
        "humanfriendly==9.1",
        "pyyaml==5.4.1",
        "zstd==1.5.0.2",
        "psutil==5.8.0",
        "tqdm==4.61.2",
        "SharedArray==3.2.1",
        "indexed-gzip==1.6.1",
        "lmdb==1.2.1",
        "msgpack==1.0.2",
    ],
    ext_modules=cythonize(
        [
            "metaprofi/lib/bloomfilter_cython.pyx",
            "metaprofi/lib/utilities_cython.pyx",
        ],
        compiler_directives={
            "language_level": "3",
            "boundscheck": False,
            "wraparound": False,
        },
    ),
    include_dirs=[np.get_include()],
    entry_points={
        "console_scripts": ["metaprofi = metaprofi.metaprofi_main:main"],
    },
)
