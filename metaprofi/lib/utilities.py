"""MetaProFi utilities module
"""

import gzip
import math
import os
import shutil
import uuid
from functools import reduce
import numpy as np
import yaml
import zarr
import zstd
from bitarray import bitarray
from humanfriendly import format_size, parse_size
from metaprofi.lib.constants import (
    BLOOMFILTER_DATASET_NAME,
    BLOOMFILTER_INDEX_DATASET_NAME,
    INDEX_STORE,
    MATRIX_STORE,
    METADATA_DATASET_NAME,
)
from numcodecs import PackBits
from psutil import cpu_count, virtual_memory
from pyfastx import Fasta, Fastq


# ---------- Config related ---------- #
def get_config(config_file, config_check=True):
    """Parse YAML config file

    Args:
      config_file: Configuration file
      config_check: To validate the user provided config file (Default value = True)

    Returns: Configuration dictionary

    """
    with open(config_file, "r") as infile:
        loaded_config = yaml.load(infile, Loader=yaml.FullLoader)

    if config_check:
        config = check_config(loaded_config)
    else:
        config = loaded_config

    # Determine dtype
    dtype = check_dtype_to_use(config["m"])
    config["DTYPE"] = dtype
    return config


# ---------- Kmers/Fasta/Fastq related ---------- #
def get_fasta_object(fasta_file, build_index=True):
    """Generate an index for the given fasta file and return a fasta object

    Args:
      fasta_file: Fasta sequence file
      build_index: Whether to build fasta index (Default value = True)

    Returns: Fasta object

    """
    fasta_object = Fasta(fasta_file, build_index=build_index)
    return fasta_object


def extract_kmers_from_kmers_txt(kmers_txt):
    """Extract kmers from a kmer text file

    Args:
      kmers_txt: K-mers text file

    Returns: Generator object of K-mers

    """
    for kmer in open(kmers_txt):
        yield kmer.strip()


# ---------- (De)Encoding/(De)Compression ---------- #
def zstd_compress(np_array):
    """Compress the given boolean numpy array using Zstandard

    Args:
      np_array: 1D numpy boolean array

    Returns: Z-standard compressed array in bytes

    """
    return zstd.compress(PackBits().encode(np_array), 9, 2)


def zstd_compress_cat(old_array, new_array):
    """Concatenate 2 1D NumPy array and compress them

    Args:
      old_array: Compressed (using zstd_compress) 1D NumPy boolean array
      new_array: 1D numpy boolean array

    Returns: Z-standard compressed array in bytes

    """
    return zstd_compress(
        np.append(zstd_decompress(old_array, np_array=True), new_array)
    )


def zstd_decompress(bytes_string, np_array=False):
    """Decompress the Zstandard compressed numpy array and return a bitarray

    Args:
      bytes_string: Z-standard compressed array in bytes
      np_array: Boolean (True: to return as a NumPy boolean array, False: to return as a bitarray)

    Returns: NumPy boolean array or Bitarray

    """
    decompressed = PackBits().decode(zstd.decompress(bytes_string))

    if np_array:
        return decompressed

    bit_array = bitarray()
    bit_array.pack(decompressed.tobytes())
    return bit_array


# ---------- Calculations reated ---------- #
def calculate_chunksize(config, number_of_samples):
    """Calculates the chunk size based on the user provided configuration in the config file

    Args:
      config: Config dictionary containing all the necessary configuration
      number_of_samples: Number of samples to Index

    Returns: Chunk size for Zarr store

    """
    # We want to fit the number of samples processed per iteration into the
    # user allocated memory, so we apply chunking (this is also used for
    # parallelization)
    # 1) Since Python does not have bit dtype, we can use bool but that still
    # consumes 1 byte, so the best would be to use uint8 dtype as it is also 1
    # byte, we can pack the bloom filter and reduce the memory requirement to
    # 1/8th and process lot of samples per iteration and save memory. May
    # affect the performance while indexing, that's a compromise we have to
    # make!
    packed_bytes_per_bloomfilter = config["packed_bytes_per_bloomfilter"]
    total_required_bytes = packed_bytes_per_bloomfilter * number_of_samples
    num_chunks = math.ceil(
        total_required_bytes
        / (parse_size(config["max_memory"], binary=True) * (85 / 100))
    )
    chunk_size_samples = round(number_of_samples / num_chunks)

    # 2) Calculate chunk size rows such that when we read the chunk from zarr
    # and unpack them into bits (chunk_size_rows * 8) we still are within the
    # memory allocated by the user
    chunk_size_rows = round((packed_bytes_per_bloomfilter / num_chunks) / 8)

    # 3) Make sure to utilize all available cores while writing row chunks to
    # the disk (downside: may reduce compression ratio because of smaller chunk
    # size, but gains performance)
    # if math.ceil(packed_bytes_per_bloom_filter / chunk_size_rows) < config["nproc"]:
    #     chunk_size_rows = math.ceil(packed_bytes_per_bloom_filter / config["nproc"])

    # 4) Zarr numcodecs Blosc limitation: Codec does not support buffer size > 2147483631 bytes (chunk size should be less than 2 or 2.1 GiB)
    # Currently commented out because of switching to Zstd compressor for Zarr from Blosc
    # chunk_size = chunk_size_rows * chunk_size_samples
    # if chunk_size > 2147481647:
    #     rows = math.floor(2147481647 / chunk_size_samples)
    #     columns = chunk_size_samples
    # else:
    rows = chunk_size_rows
    columns = chunk_size_samples

    return rows, columns, chunk_size_rows, chunk_size_samples


def calculate_index_chunksize(config):
    """Calculates index chunk size based on the user provided configuration in the config file

    Args:
      config: Config dictionary containing all the necessary configuration

    Returns: Index chunk size for Zarr store

    """
    # TODO: Need a better way to calculate the row chunk size for index
    # store. Number of chunks is equal to the number of files, so here we
    # try to limit it to max 500_000. We do not want to produce more than
    # 500_000 files due to various disk limitations (inodes).
    # This chunk size also affects the query performance during the retrieval
    # of the rows from the index store.
    if config["m"] <= 10_000_000:
        index_chunk_rows = math.ceil(config["m"] / 5_000)
    if 10_000_000 <= config["m"] < 100_000_000:
        index_chunk_rows = math.ceil(config["m"] / 10_000)
    if 100_000_000 <= config["m"] < 1_000_000_000:
        index_chunk_rows = math.ceil(config["m"] / 100_000)
    if 1_000_000_000 <= config["m"] < 10_000_000_000:
        index_chunk_rows = math.ceil(config["m"] / 400_000)
    if config["m"] >= 10_000_000_000:
        index_chunk_rows = math.ceil(config["m"] / 500_000)

    return index_chunk_rows


# ---------- Checks related ---------- #
def check_config(config):
    """Check/Validate the user provided config in the config file

    Args:
      config: Config dictionary

    Returns: Validated config dictionary

    """
    num_procs = cpu_count() - 1
    system_memory = virtual_memory().total
    session_id = uuid.uuid4().hex

    # When the user did not set max_memory to use in the config file, we try to use 50% of the available system memory as default
    max_memory = round(((system_memory * (9.31 * 10 ** -10)) * (50 / 100)), 1)

    # Mandatory config to be present in the config file
    # 1) Number of hash functions to apply on each k-mer
    if not "h" in config or not config["h"]:
        raise KeyError("Number of hash functions 'h' is not set in the config file")

    if config["h"] == 0:
        raise ValueError(
            "Number of hash functions 'h' in the config file cannot be zero"
        )

    # 2) K-mer size
    if not "k" in config or not config["k"]:
        raise KeyError("Size of kmer 'k' is not set in the config file")

    if config["k"] == 0:
        raise ValueError("Size of kmer 'k' in the config file  cannot be zero")

    # 3) Size of the bloom filter and related checks
    if not "m" in config or not config["m"]:
        raise KeyError("Size of bloom filter 'm' is not set in the config file")

    if config["m"] == 0:
        raise ValueError(
            "Size of the bloom filter 'm' in the config file cannot be zero"
        )

    # Find the number of zeros to be padded to make bloom filter a multiple of 8
    if config["m"] % 8 == 0:
        number_of_zeros_to_pad = 0
    else:
        number_of_zeros_to_pad = 8 - (config["m"] % 8)

    packed_bytes_per_bloomfilter = int((config["m"] + number_of_zeros_to_pad) / 8)
    config["number_of_zeros_to_pad"] = number_of_zeros_to_pad
    config["packed_bytes_per_bloomfilter"] = packed_bytes_per_bloomfilter

    # 4) Type of the sequences being indexed (aminoacid or nucleotide)
    if not "sequence_type" in config or not config["sequence_type"]:
        raise KeyError(
            "Type of sequence in the input file is not set in the config file"
        )

    # 5) Output directory to store all the outputs and also to look for other input files and metaprofi files
    if not "output_directory" in config or not config["output_directory"]:
        raise KeyError("Output directory is not set in the config file")

    config["output_directory"] = os.path.abspath(config["output_directory"])

    if not os.path.exists(config["output_directory"]):
        os.mkdir(config["output_directory"])

    # Optional config expected to be present in the config file
    # 6) Matrix store file name
    if not "matrix_store_name" in config or not config["matrix_store_name"]:
        matrix_store = f"{config['output_directory']}/{MATRIX_STORE}"
        config["matrix_store_name"] = matrix_store
        print(
            f"INFO: Matrix store filename is not available in the config file, using default value: {matrix_store}"
        )
    else:
        matrix_store_basename = (
            os.path.splitext(os.path.basename(config["matrix_store_name"]))[0] + ".zarr"
        )
        config[
            "matrix_store_name"
        ] = f"{config['output_directory']}/{matrix_store_basename}"

    # 7) Index store file name
    if not "index_store_name" in config or not config["index_store_name"]:
        config["index_store_name"] = f"{config['output_directory']}/{INDEX_STORE}"
        print(
            f"INFO: Index store filename is not available in the config file, using default value: {config['index_store_name']}"
        )
    else:
        index_store_basename = (
            os.path.splitext(os.path.basename(config["index_store_name"]))[0] + ".zarr"
        )
        config[
            "index_store_name"
        ] = f"{config['output_directory']}/{index_store_basename}"

    # 8) Number of cores to use for parallelization
    if not "nproc" in config or not config["nproc"]:
        config["nproc"] = num_procs
        print(
            f"INFO: Number of cores (nproc) to use is not set in the config file, using default value: {num_procs}"
        )

    # 9) Max memory to be used by MetaProFi
    if not "max_memory" in config or not config["max_memory"]:
        config["max_memory"] = f"{max_memory}GiB"
        print(
            f"INFO: Maximum system memory (max_memory) to use is not set in the config file, using default value: {max_memory} GiB"
        )

    # # Make sure that the user does not set more than 85% of the system memory
    # if not system_memory > (
    #     parse_size(config["max_memory"], binary=True)
    #     + (parse_size(config["max_memory"], binary=True) * (15 / 100))
    # ):
    #     config["max_memory"] = f"{max_memory}GiB"
    #     print(
    #         f"INFO: Don't set more than 85% of the system memory in the config file. Using default value: {max_memory} GiB"
    #     )

    # 10) Check if the max_memory is enough to process at least 1 bloom filter and if not raise error
    if (
        not (
            parse_size(config["max_memory"], binary=True)
            // packed_bytes_per_bloomfilter
        )
        >= 1
    ):
        raise ValueError(
            f"Not enough memory to process even a single sample. Requires at least {int(packed_bytes_per_bloomfilter * (110/100))} bytes to build a single bloom filter but given only {parse_size(config['max_memory'], binary=True)} bytes."
        )

    # 11) Calculate the availability of the shared memory "/dev/shm" which will be used as the shared memory for storing the NumPy array
    config[
        "POSIX_SHARED_MEMORY_ARRAY"
    ] = f"shm://metaprofi_shared_memory_array-{session_id}"
    stat_vfs = os.statvfs("/dev/shm")
    available_shared_mem = (stat_vfs.f_bsize * stat_vfs.f_bavail) * (90 / 100)
    if parse_size(config["max_memory"], binary=True) > available_shared_mem:
        print(
            f"WARNING: Given max_memory ({parse_size(config['max_memory'], binary=True)} bytes) in config file exceeds memory available in '/dev/shm' ({available_shared_mem} bytes). So, using '/tmp' as the shared memory store and this may slow down MetaProFi's performance! You maybe able to free up space by deleting files in '/dev/shm' or cancel the current MetaProFi session and then re-start it with less max_memory set in the config file"
        )
        config[
            "POSIX_SHARED_MEMORY_ARRAY"
        ] = f"file:///tmp/metaprofi_shared_memory_array-{session_id}"

    # Check if required database config exists in the config file
    # if not config["storage-config"]:
    #     raise KeyError(
    #         "Database required configuration is not available in the input config file"
    #     )

    # if not all(
    #     key in config["storage-config"]
    #     for key in [
    #         "db_create",
    #         "db_hostname",
    #         "db_username",
    #         "db_password",
    #         "db_databasename",
    #     ]
    # ):
    #     raise KeyError(
    #         "Not all configuration required for the database connection is available in the config file"
    #     )

    return config


def check_dtype_to_use(number):
    """Determine which unsigned integer data type to be used

    Args:
      number: Number (bloom filter size for example)

    Returns: Numpy unsigned integer data type to use

    """
    if number < 255:
        dtype = "uint8"
    elif number < 65535:
        dtype = "uint16"
    elif number < 4294967295:
        dtype = "uint32"
    else:
        dtype = "uint64"
    return dtype


def check_file_type(in_file):
    """Checks and returns the input sequence file type

    Args:
      in_file: User provided sequence input file path

    Returns: File type (fasta or fastq)

    """
    compressed = is_gzip_file(in_file)
    line = ""

    if compressed:
        with gzip.open(in_file) as in_f:
            line = in_f.readline()
    else:
        with open(in_file, "rb") as in_f:
            line = in_f.readline()

    if line.startswith(b">"):
        return "fasta"
    elif line.startswith(b"@"):
        return "fastq"
    else:
        raise ValueError("Input file should be either a Fasta or a Fastq file format")


def check_zarr_store(zarr_store_dir):
    """Check the Zarr store and its metadata

    Args:
      zarr_store_dir: Zarr store directory path

    Returns: None

    """
    if not os.path.exists(zarr_store_dir):
        raise FileNotFoundError(
            f"Zarr store direcotry '{zarr_store_dir}' does not exist"
        )

    with zarr.open(zarr_store_dir, mode="r") as zarr_store:
        if not zarr_store.get(BLOOMFILTER_DATASET_NAME) or not zarr_store.get(
            METADATA_DATASET_NAME
        ):
            raise ValueError(
                f"The zarr store ({MATRIX_STORE}) does not contain the required {BLOOMFILTER_DATASET_NAME} and/or {METADATA_DATASET_NAME} datasets."
            )


def check_zarr_index_store(zarr_index_store_dir):
    """Check the Zarr store and its metadata

    Args:
      zarr_store_dir: Zarr store directory path

    Returns: None

    """
    if not os.path.exists(zarr_index_store_dir):
        raise FileNotFoundError(
            f"Zarr store direcotry '{zarr_index_store_dir}' does not exist"
        )

    with zarr.open(zarr_index_store_dir, mode="r") as zarr_store:
        if not zarr_store.get(BLOOMFILTER_INDEX_DATASET_NAME) or not zarr_store.get(
            METADATA_DATASET_NAME
        ):
            raise ValueError(
                f"The zarr store ({INDEX_STORE}) does not contain the required {BLOOMFILTER_INDEX_DATASET_NAME} and/or {METADATA_DATASET_NAME} datasets."
            )


def check_seq(sequence):
    """Checks if a given sequence is DNA or protein

    Args:
      sequence: Sequence

    Returns: Type of the sequence (nucleotide or aminoacid)

    """
    dna_nucleotides = "ATGCN"

    if all(i in dna_nucleotides for i in sequence):
        return "nucleotide"

    return "aminoacid"


def check_dir(dir_name):
    """Checks if given directory exists

    Args:
      dir_name: Name of the directory to check

    Returns: Boolean

    """
    return bool(os.path.isdir(dir_name))


# ---------- Others ---------- #
def get_file_reader(input_file):
    """Return a file reader (Fasta/Fastq) based on the input file

    Args:
      input_file: Path to the input file (Fasta/Fastq)

    Returns: Fasta/Fastq file reader

    """
    # Based on the file type use the appropriate file reader
    file_type = check_file_type(input_file)
    if file_type == "fastq":
        file_reader = Fastq
    elif file_type == "fasta":
        file_reader = Fasta
    return file_reader


def is_gzip_file(input_file):
    """Check if a file is gzip compressed

    Args:
      in_file: Input file to check

    Returns: Boolean

    """
    gzip_magic_number = b"\x1f\x8b"

    with open(input_file, "rb") as in_f:
        return bool(in_f.read(2) == gzip_magic_number)


def splitjobs(items_list, num_chunks):
    """Given a list of items and number of chunks, will return the list in N chunks

    Args:
      items_list: List object
      num_chunks: Number of chunks to split the items list into

    Returns: N small chunks of the input items list

    """
    # If len(items_list) is < num_chunks, np.array_split will provide some empty lists and we don't want that
    num_items = len(items_list)
    if num_items < num_chunks:
        num_chunks = num_items
    return [
        items_list[i * num_items // num_chunks : (i + 1) * num_items // num_chunks]
        for i in range(num_chunks)
    ]


def chunks(samples_list, chunk_size):
    """Given a sample list and chunk size, returns a chunk from the list

    Args:
      samples_list: List object
      chunk_size: Number of items to have in each chunk

    Returns: Generator object of chunks of the input list

    """
    for i in range(0, len(samples_list), chunk_size):
        yield samples_list[i : i + chunk_size]


def bitwise_and(bloomfilters):
    """Performs a bitwise and operation given a list of bloomfilters/bitarrays

    Args:
      bloomfilters: A list of bloom filters (bitarrays)

    Returns: Bit array

    """
    if len(bloomfilters) == 1:
        return bloomfilters[0]
    return reduce(lambda x, y: x & y, bloomfilters)


def get_summary_from_index_store(config_file, as_dict=False):
    """Extract the configuration and metadata info from the index store

    Args:
      config_file: Config file
      as_dict:  Format of the return results [False]

    Returns: Yaml or Dictionary formatted summary of the index and samples

    """
    config = get_config(config_file, config_check=True)
    index_store_dir = config["index_store_name"]

    # Check if zarr store contains required datasets
    check_zarr_index_store(index_store_dir)

    with zarr.open(index_store_dir, mode="r") as index_store:
        metadata = index_store.get(METADATA_DATASET_NAME)
        index_dataset = index_store.get(BLOOMFILTER_INDEX_DATASET_NAME)
        info = {
            "Bloom filter configuration": {
                "Kmer size": metadata.attrs["kmer_size"],
                "Bloom filter size": metadata.attrs["bloom_filter_size"],
                "No. of hashes": metadata.attrs["number_of_hash_functions"],
            },
            "Metadata": {
                "Zarr store created on": metadata.attrs["created_on"],
                "No. of updates": metadata.attrs["update_count"],
                "Last updated on": metadata.attrs["last_updated_on"],
                "No. of datasets": len(list(index_store.array_keys())),
                "Dataset names": ", ".join(list(index_store.array_keys())),
                "No. of samples": metadata.attrs["number_of_samples"],
                "Compressor configuration": metadata.attrs["compressor"],
                "No. of bytes required to store without compression": f"{index_dataset.nbytes} ({format_size(index_dataset.nbytes, binary=True)})",
                "No. of bytes used to store with compression": f"{index_dataset.nbytes_stored} ({format_size(index_dataset.nbytes_stored, binary=True)})",
            },
        }

    if as_dict:
        return info
    return yaml.dump(info, sort_keys=False)


def reverse_complement(seq):
    """Reverse complement of the given sequence

    Args:
      seq: Nucleotide sequence

    Returns: Reverse complement of the given sequence
    """
    # https://arep.med.harvard.edu/labgc/adnan/projects/Utilities/revcomp.html
    # DNA IUPAC ambiguity codes
    table = str.maketrans("ACGTMRWSYKVHDBXN", "TGCAKYWSRMBDHVXN")
    return seq.translate(table)[::-1]


# ---------- Cleanup related ---------- #
def cleanup_dir(directory):
    """Removes the given directory

    Args:
      directory: Directory to delete

    Returns: None

    """
    try:
        if os.path.exists(directory):
            shutil.rmtree(directory)
    except OSError as excep:
        print("Error: %s : %s" % (directory, excep.strerror))


def cleanup_file(filename):
    """Remove the given file

    Args:
      filename: File to delete

    Returns: None

    """
    try:
        if os.path.exists(filename):
            os.remove(filename)
    except OSError as excep:
        print("Error: %s : %s" % (filename, excep.strerror))


# ---------- Utility Classes ---------- #
class BColors:
    """ANSI escape sequences for colored print statements
    Code taken from https://stackoverflow.com/questions/287871/how-to-print-colored-text-in-python

    Args:

    Returns:

    """

    HEADER = "\033[95m"
    OKBLUE = "\033[94m"
    OKGREEN = "\033[92m"
    WARNING = "\033[93m"
    FAIL = "\033[91m"
    ENDC = "\033[0m"
    BOLD = "\033[1m"
    UNDERLINE = "\033[4m"
