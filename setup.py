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
        "Programming Language :: Python :: 3.8",
        "Operating System :: POSIX :: Linux",
        "Natural Language :: English",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    include_package_data=True,
    python_requires=">=3.8",
    packages=find_packages(),
    install_requires=[
        "Cython==0.29.28",
        "numpy==1.22.3",
        "zarr==2.11.1",
        "pyfastx==0.8.4",
        "bitarray==2.4.0",
        "humanfriendly==9.1",
        "pyyaml==5.4.1",
        "zstd==1.5.1.0",
        "psutil==5.9.0",
        "tqdm==4.61.2",
        "SharedArray==3.2.1",
        "indexed-gzip==1.6.4",
        "lmdb==1.3.0",
        "msgpack==1.0.3",
        "msgspec==0.5.0",
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
    zip_safe=False,
)
