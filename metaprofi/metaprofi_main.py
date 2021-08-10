"""MetaProFi input handling module
"""

import argparse
import sys
import time
from metaprofi.lib.zarrstore import ZarrStore
from metaprofi.lib.zarrstore_seq import ZarrStoreSeq
from metaprofi.lib.search_index import search_index
import SharedArray as sarray
from metaprofi.lib.build_index_zarr import build_index, update_index
from metaprofi.lib.utilities import (
    get_config,
    get_summary_from_index_store,
)
from metaprofi.version import __version__


def main():
    """Handles all inputs and related"""
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(
        dest="subcommand",
        title="Available subcommands:\nNOTE: Use only 1 of these subcommands at a time and to know further information about a subcommand, please use the subcommand followed by -h",
        required=True,
    )

    # For building both zarr store and index database
    build_parser = subparsers.add_parser(
        "build",
        help="To build bloom filter matrix and to create index store. Note: All in one command (builds bloom filter matrix and index store)",
    )
    build_parser.add_argument(
        dest="input_file",
        help="Input file containing lines each of which specifies a sample identifier and one or more compressed or uncompressed fasta or fastq file path. Refer README for format",
    )
    build_parser.add_argument(dest="config_file", help="Configuration YAML file")

    # For building both zarr store and index at sequence level
    build_parser = subparsers.add_parser(
        "build-seq",
        help="To build bloom filter matrix (every sequence in the input file will be considered as a sample) and to create index store. Note: All in one command (builds bloom filter matrix and index store)",
    )
    build_parser.add_argument(
        dest="input_file",
        help="FASTA/FASTQ (.GZ) file containing one or more sequences",
    )
    build_parser.add_argument(dest="config_file", help="Configuration YAML file")

    # For appending new data to the zarr store and the index
    update_parser = subparsers.add_parser(
        "update",
        help="To build bloom filter matrix for the new samples and to append/update the index with the new data",
    )
    update_parser.add_argument(
        dest="input_file",
        help="Input file containing lines each of which specifies a sample identifier and one or more compressed or uncompressed fasta or fastq file path. Refer README for format",
    )
    update_parser.add_argument(dest="config_file", help="Configuration YAML file")

    # For updating both zarr store and index at sequence level
    update_parser = subparsers.add_parser(
        "update-seq",
        help="To build bloom filter matrix (every sequence in the input file will be considered as a sample) for the new samples to append/update the index with the new data.",
    )
    update_parser.add_argument(
        dest="input_file",
        help="FASTA/FASTQ (.GZ) file containing one or more sequences",
    )
    update_parser.add_argument(dest="config_file", help="Configuration YAML file")

    # For building BF matrix store
    build_zarr_parser = subparsers.add_parser(
        "build_matrix",
        help="To build bloom filter matrix",
    )
    build_zarr_parser.add_argument(
        dest="input_file",
        help="Input file containing lines each of which specifies a sample identifier and one or more compressed or uncompressed fasta or fastq file path. Refer README for format",
    )
    build_zarr_parser.add_argument(dest="config_file", help="Configuration YAML file")

    # For building the index store using the pre-built zarr store
    build_db_parser = subparsers.add_parser(
        "build_index",
        help="To build the index",
    )
    build_db_parser.add_argument(dest="config_file", help="Configuration YAML file")

    # For appending new data to the the index store
    update_db_parser = subparsers.add_parser(
        "update_index",
        help="To update index with the new data",
    )
    update_db_parser.add_argument(dest="config_file", help="Configuration YAML file")

    # For searching/querying the index store
    search_db_parser = subparsers.add_parser(
        "search_index", help="To search/query the index"
    )
    group = search_db_parser.add_mutually_exclusive_group()
    group.add_argument(
        "-s",
        dest="sequence",
        help="Input sequence to search against the index store",
    )
    group.add_argument(
        "-f",
        dest="input_file",
        help="An uncompressed input fasta/fastq file containing sequence(s) to search against the index database",
    )
    search_db_parser.add_argument(
        "-i",
        dest="seq_type",
        required=True,
        choices=["nucleotide", "aminoacid"],
        help="Type of the query sequence (nucleotide or aminoacid)",
    )
    search_db_parser.add_argument(
        "-t",
        dest="threshold",
        type=int,
        help="Threshold value between 0 and 100 to invoke approximate search",
        default=100,
    )
    search_db_parser.add_argument(dest="config_file", help="Configuration YAML file")

    # For extracting summary about the data from the index store
    get_zarr_summary = subparsers.add_parser(
        "summary",
        help="For summary",
    )
    get_zarr_summary.add_argument(dest="config_file", help="Configuration YAML file")

    # For MetaProFi version
    parser.add_argument(
        "--version",
        action="version",
        help="Show MetaProFi version number",
        version=f"MetaProFi version: {__version__}",
    )

    # ----------------------------------------------------------------------- #
    # process arguments and function call
    args = parser.parse_args(args=None if sys.argv[1:] else ["--help"])

    # For Building both the database and the zarr store
    if args.subcommand == "build":
        config = get_config(args.config_file, config_check=True)
        posix_shared_memory_array = config["POSIX_SHARED_MEMORY_ARRAY"]
        try:
            start_time = time.time()
            zarr_store = ZarrStore(args.input_file, config, update=False)
            print("Building Bloom Filter matrix")
            zarr_store.create_zarr()
            print(f"Matrix construction took: {time.time() - start_time}s")
            print("Building the index store")
            start_time = time.time()
            build_index(config)
            print(f"Index construction took: {time.time() - start_time}s")
        except KeyboardInterrupt:
            sarray.delete(posix_shared_memory_array)
            sys.exit()

    # For building both zarr store and index at sequence level
    if args.subcommand == "build-seq":
        config = get_config(args.config_file, config_check=True)
        posix_shared_memory_array = config["POSIX_SHARED_MEMORY_ARRAY"]
        try:
            start_time = time.time()
            zarr_store = ZarrStoreSeq(args.input_file, config, update=False)
            print("Building Bloom Filter matrix")
            zarr_store.create_zarr()
            print(f"Matrix construction took: {time.time() - start_time}s")
            print("Building the index store")
            start_time = time.time()
            build_index(config)
            print(f"Index construction took: {time.time() - start_time}s")
        except KeyboardInterrupt:
            sarray.delete(posix_shared_memory_array)
            sys.exit()

    # For updating the BF matrix (we build new matrix for every new dataset) and the database
    if args.subcommand == "update":
        config = get_config(args.config_file, config_check=True)
        posix_shared_memory_array = config["POSIX_SHARED_MEMORY_ARRAY"]
        try:
            zarr_store = ZarrStore(args.input_file, config, update=True)
            print("Building Bloom Filter matrix for the new samples")
            zarr_store.create_zarr()
            print("Updating the index store with the new data")
            update_index(config)
        except KeyboardInterrupt:
            sarray.delete(posix_shared_memory_array)
            sys.exit()

    # For updating both zarr store and index at sequence level
    if args.subcommand == "update-seq":
        config = get_config(args.config_file, config_check=True)
        posix_shared_memory_array = config["POSIX_SHARED_MEMORY_ARRAY"]
        try:
            start_time = time.time()
            zarr_store = ZarrStoreSeq(args.input_file, config, update=True)
            print("Building Bloom Filter matrix for the new samples")
            zarr_store.create_zarr()
            print(f"Matrix construction took: {time.time() - start_time}s")
            print("Updating the index store with the new data")
            start_time = time.time()
            update_index(config)
            print(f"Index update took: {time.time() - start_time}s")
        except KeyboardInterrupt:
            sarray.delete(posix_shared_memory_array)
            sys.exit()

    # For Building BF matrix
    if args.subcommand == "build_matrix":
        config = get_config(args.config_file, config_check=True)
        posix_shared_memory_array = config["POSIX_SHARED_MEMORY_ARRAY"]
        try:
            zarr_store = ZarrStore(args.input_file, config, update=False)
            print("Building Bloom Filter matrix")
            zarr_store.create_zarr()
        except KeyboardInterrupt:
            sarray.delete(posix_shared_memory_array)
            sys.exit()

    # For Building the index using the pre-built BF matrix
    if args.subcommand == "build_index":
        config = get_config(args.config_file, config_check=True)
        posix_shared_memory_array = config["POSIX_SHARED_MEMORY_ARRAY"]
        try:
            print("Building the index store")
            build_index(config)
        except KeyboardInterrupt:
            sarray.delete(posix_shared_memory_array)
            sys.exit()

    # For updating the the index using the newly built BF matrix
    if args.subcommand == "update_index":
        config = get_config(args.config_file, config_check=True)
        posix_shared_memory_array = config["POSIX_SHARED_MEMORY_ARRAY"]
        try:
            print("Updating the index store with the new data")
            update_index(config)
        except KeyboardInterrupt:
            sarray.delete(posix_shared_memory_array)
            sys.exit()

    # For searching the index
    if args.subcommand == "search_index":
        assert (
            0 <= args.threshold <= 100
        ), f"Threshold should be a value between 0 and 100, you have provided {args.threshold}"

        if args.sequence:
            query = args.sequence
            inp_seq = True
        elif args.input_file:
            inp_seq = False
            query = args.input_file
        else:
            raise ValueError(
                "search_db requires either a sequence or a FASTA/FASTQ (.GZ) file containing sequence(s)"
            )

        config = get_config(args.config_file, config_check=True)
        search_index(
            config,
            query,
            inp_seq,
            args.seq_type,
            args.threshold,
            write_to_file=True,
        )

    # For getting summary about the data using zarr store
    if args.subcommand == "summary":
        print(get_summary_from_index_store(args.config_file))


if __name__ == "__main__":
    main()
