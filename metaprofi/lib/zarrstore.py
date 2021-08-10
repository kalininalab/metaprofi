"""MetaProFi zarr store Bloom Filter (BF) matrix file-level module
"""

import os
from datetime import datetime
import re
import sys
import numpy as np
import gzip
import zarr
import SharedArray as sarray
from multiprocessing import Process
from tqdm import tqdm
from numcodecs import blosc, Zstd
from metaprofi.lib.bloomfilter_cython import bloomfilter_cython, check_kmer_exists
from metaprofi.lib.constants import (
    BLOOMFILTER_DATASET_NAME,
    COMPRESSION_ALGO,
    COMPRESSION_LEVEL,
    METADATA_DATASET_NAME,
)
from metaprofi.lib.utilities import (
    calculate_chunksize,
    calculate_index_chunksize,
    check_dir,
    check_zarr_index_store,
    check_zarr_store,
    chunks,
    cleanup_dir,
    splitjobs,
    is_gzip_file,
)

# Global
POSIX_SHARED_MEMORY_ARRAY = ""


class ZarrStore:
    """Zarr BF matrix store

    Args:
      input_file: Input file (.GZ) contianing a sample identifier and fasta/fastq file path. One line per sample.
      config: Config dictionary containing all the necessary configuration
      update: Boolean indicating whether this is run for an update index session or build session.

    NOTE:
      This is used at file-level meaning every FASTA/FASTQ file will be considered as a sample (or multiple files per sample) so 1 Bloom Filter (BF) per file or 1 BF per sample will be created in the BF matrix.

    """

    def __init__(self, input_file, config, update=False):
        self.config = config
        self.input_file = input_file
        self.zarr_dir = self.config["matrix_store_name"]

        global POSIX_SHARED_MEMORY_ARRAY
        POSIX_SHARED_MEMORY_ARRAY = self.config["POSIX_SHARED_MEMORY_ARRAY"]

        # Check for endianness during update sessions
        if update:
            index_store_dir = self.config["index_store_name"]
            check_zarr_index_store(index_store_dir)
            index_store = zarr.open(index_store_dir, mode="r")
            index_metadata_dataset = index_store.get(METADATA_DATASET_NAME)
            current_endianness = sys.byteorder
            index_endianness = index_metadata_dataset.attrs["endianness"]

            if index_endianness != current_endianness:
                raise ValueError(
                    f"Can't proceed with update request. Index endianness ({index_endianness} endian) does not match the current hardware endianness ({current_endianness} endian)."
                )

        # Get data from input file
        self.samples_path_list, self.sample_identifiers = self.parse_input_file()

        # Calculate chunksize
        (
            self.rows,
            self.columns,
            self.chunk_size_rows,
            self.chunk_size_samples,
        ) = calculate_chunksize(config, len(self.sample_identifiers))

        # During index update sessions, we need to keep the memory under the
        # limit in scenarios where a super large dataset (millions) was indexed
        # first and then in the update session, a very small dataset (few
        # thousands) is being indexed or in the reverse order. In such
        # scenarios when reading data from the index store for merging/
        # appending with the new data MetaProFi may require/consume extra RAM
        # than what the user set, and this is not ideal. To avoid this it is
        # best to simply use the min chunk size rows this way we can keep
        # memory consumption still under the user set limit but may lose some
        # performance and compression ratio (tradeoff!).
        if update:
            self.rows = min(index_metadata_dataset.attrs["chunk_size_rows"], self.rows)

        # Packed bloom filter size and number of zeros to pad
        self.number_of_zeros_to_pad = self.config["number_of_zeros_to_pad"]
        self.packed_bf_size = self.config["packed_bytes_per_bloomfilter"]

    def parse_input_file(self):
        """Process and parse the input file containing sample ids and sample path"""

        print("MetaProFi: checking and reading input file ...")

        if os.stat(self.input_file).st_size == 0:
            raise ValueError(f"Empty input file: {self.input_file}")

        # fmt: off
        sample_id_not_allowed = ['/', '\\', '*', '.', '?', ',', '[', ']', '{', '}', '!', '@', '#', '$', '%', '^', '&', '(', ')', '+', '=', '<', '>']
        # fmt: on

        nonexistent_files = []
        samples_path_list = []

        compressed = is_gzip_file(self.input_file)

        if compressed:
            file_handle = gzip.open(self.input_file)
        else:
            file_handle = open(self.input_file, "r")

        for line in file_handle:
            if line.strip().startswith("#"):
                continue

            if line.count(":") > 1:
                raise ValueError(
                    "Only one sample identifier per line is allowed in the input file. Refer README for the format."
                )

            if line.count(":") < 1:
                raise ValueError(
                    "Every line in the input file should contain one sample identifier followed by the path. Refer README for the format."
                )

            row = re.split(": |; ", line.strip())
            if any(i in row[0] for i in sample_id_not_allowed):
                raise ValueError(
                    "Sample identifiers in the input file contain special characters that are not allowed. Refer README for the format"
                )

            if row[1:]:
                for file_path in row[1:]:
                    if not os.path.exists(file_path):
                        nonexistent_files.append(file_path)

            # Make sure there is at least 1 sequence containing a k-mer
            # If no k-mer found in the input file then we will just skip that input
            if row[1:]:
                if check_kmer_exists(row[1:], self.config["k"]):
                    size = sum([os.stat(f).st_size for f in row[1:]])
                    row.append(size)
                    samples_path_list.append(row)

        # Sort input files by their size
        samples_path_list.sort(key=lambda x: x[-1])

        # Delete the size value (NR after sorting)
        _ = [sample.pop() for sample in samples_path_list]

        # Extract sample identifiers
        sample_identifiers = [sample.pop(0) for sample in samples_path_list]

        if len(sample_identifiers) != len(samples_path_list):
            raise ValueError(
                f"Number of samples '{len(samples_path_list)}' is not equal to the number of sample identifiers '{len(sample_identifiers)}'"
            )

        if nonexistent_files:
            out_file = f"{self.config['output_directory']}/nonexistent_file_paths.txt"
            with open(out_file, "w") as non_file:
                for non_f in nonexistent_files:
                    non_file.write(f"{non_f}\n")
                raise FileNotFoundError(
                    f"Found non existant input file paths, check: {out_file}, please fix the input file and re-run MetaProFi."
                )

        file_handle.close()
        return samples_path_list, sample_identifiers

    def create_zarr(self):
        """Create and builds the Zarr store"""
        # Check if zarr store already exists
        if check_dir(self.zarr_dir):
            raise FileExistsError(f"Matrix store '{self.zarr_dir}' already exists")

        # Create Zarr store and necessary arrays/groups
        store = zarr.DirectoryStore(self.zarr_dir)
        sync_dir = os.path.splitext(self.zarr_dir)[0] + ".sync"
        synchronizer = zarr.ProcessSynchronizer(sync_dir)
        blosc.use_threads = False
        # compressor = blosc.Blosc(
        #     cname=COMPRESSION_ALGO,
        #     clevel=COMPRESSION_LEVEL,
        #     shuffle=blosc.Blosc.BITSHUFFLE,
        # )
        compressor = Zstd(level=COMPRESSION_LEVEL)
        bf_group = zarr.group(store=store)
        bf_dataset = bf_group.create_dataset(
            name=BLOOMFILTER_DATASET_NAME,
            shape=(self.packed_bf_size, len(self.samples_path_list)),
            chunks=(self.rows, self.columns),
            compressor=compressor,
            synchronizer=synchronizer,
            dtype="uint8",
        )
        metadata = bf_group.create_dataset(
            name=METADATA_DATASET_NAME,
            shape=(len(self.sample_identifiers)),
            chunks=(50_000),
            compressor=compressor,
            dtype="str",
        )

        # Insert sample identifiers into the zarr store
        metadata[:] = np.array(self.sample_identifiers, dtype="object")

        # Index store chunk size
        index_chunk_rows = calculate_index_chunksize(self.config)

        # Add some metadata
        metadata.attrs["created_on"] = datetime.today().strftime("%d-%m-%Y")
        metadata.attrs["last_updated_on"] = datetime.today().strftime("%d-%m-%Y")
        metadata.attrs["update_count"] = 0
        metadata.attrs["kmer_size"] = self.config["k"]
        metadata.attrs["bloom_filter_size"] = self.config["m"]
        metadata.attrs["packed_bf_size"] = self.packed_bf_size
        metadata.attrs["number_of_zeros_padded"] = self.number_of_zeros_to_pad
        metadata.attrs["number_of_hash_functions"] = self.config["h"]
        metadata.attrs["number_of_samples"] = len(self.sample_identifiers)
        metadata.attrs["compressor"] = str(bf_dataset.compressor)
        metadata.attrs["chunk_size_rows"] = self.rows
        metadata.attrs["chunk_size_samples"] = self.columns
        metadata.attrs["sequence_type"] = self.config["sequence_type"]
        metadata.attrs["endianness"] = sys.byteorder
        metadata.attrs["index_chunk_rows"] = index_chunk_rows

        # Compute and insert Bloom filter matrix into the Zarr store
        samples = list(
            list(sample)
            for sample in zip(
                range(len(self.samples_path_list)), self.samples_path_list
            )
        )

        # Delete and recover some memory (sample identifier is already written to disk, so not needed here anymore)
        del self.sample_identifiers

        self.insert_bf_into_zarr(samples, bf_dataset)

        # Clean up
        cleanup_dir(sync_dir)

        # Log output
        print("Bloom filter matrix store is now created with all the samples")

    def insert_bf_into_zarr(self, samples, bf_dataset):
        """Create and insert the BloomFilters as columns into the Zarr store

        Args:
          samples: Samples path list
          bf_dataset: Bloom filter matrix dataset handle from Zarr store
          resize: Whether to resize the Zarr store dataset (Default value = False)

        Returns: None

        """
        # Progress bar object
        pbar = tqdm(
            total=len(samples),
            bar_format="{l_bar}{bar}{n_fmt}/{total_fmt}",
            desc="Number of samples processed",
            leave=True,
            ncols=80,
            colour="green",
        )

        # Create a POSIX shared memory NumPy boolean array
        shared_array = sarray.create(
            POSIX_SHARED_MEMORY_ARRAY,
            (self.packed_bf_size, self.chunk_size_samples),
            "uint8",
        )

        # Build bloom filters in parallel
        for chunked_samples in chunks(samples, self.chunk_size_samples):
            # shared_array reset and re-initialization
            if len(chunked_samples) != self.chunk_size_samples:
                # Delete the existing shared array
                del shared_array
                sarray.delete(POSIX_SHARED_MEMORY_ARRAY)

                # Re-create shared array
                shared_array = sarray.create(
                    POSIX_SHARED_MEMORY_ARRAY,
                    (self.packed_bf_size, len(chunked_samples)),
                    "uint8",
                )

                # Initialize shared array
                shared_array[:] = 0
            else:
                # Reset shared array
                shared_array[:] = 0

            # Column index of the shared array
            for idx, sample_list in enumerate(chunked_samples):
                sample_list.append(idx)

            # Build bloom filters in parallel
            bf_processes = []
            for samp_list in splitjobs(chunked_samples, self.config["nproc"]):
                proc = Process(
                    target=mp_bloomfilter,
                    args=(samp_list, self.config),
                )
                proc.start()
                bf_processes.append(proc)

            for proc in bf_processes:
                proc.join()

            # Write the bloom filter matrix into the Zarr store
            # Chunk rows so chunks can be inserted parallely
            row_list = list(chunks(range(self.packed_bf_size), self.rows))

            # Column index start and stop for correct insert into the Zarr store
            col_idx = (chunked_samples[0][0], chunked_samples[-1][0] + 1)

            # Write chunks in parallel
            zarr_processes = []
            for rows in splitjobs(row_list, self.config["nproc"]):
                zarr_proc = Process(
                    target=write_to_zarr_from_shared_array,
                    args=(rows, bf_dataset, col_idx),
                )
                zarr_proc.start()
                zarr_processes.append(zarr_proc)

            for zarr_proc in zarr_processes:
                zarr_proc.join()

            # Update progress bar
            pbar.update(len(chunked_samples))

        # Clean up
        del shared_array
        sarray.delete(POSIX_SHARED_MEMORY_ARRAY)

        # Close progress bar
        pbar.close()


def mp_bloomfilter(samples_list, config):
    """Helper function for creating bloom filters in parallel

    Args:
      samples_list: Samples list
      config: Config dictionary

    """
    # Attach to the POSIX shared memory NumPy boolean array
    shared_array = sarray.attach(POSIX_SHARED_MEMORY_ARRAY)

    # Bloom filters
    for sample_info in samples_list:
        file_paths = sample_info[1]
        sample_column_idx = sample_info[2]
        bloomfilter_cython(file_paths, config, shared_array, sample_column_idx)


def write_to_zarr_from_shared_array(rows, bf_dataset, columns_idx):
    """Helper function for inserting bloomfilters chunkwise in parallel

    Args:
      rows: Bloom filter rows to insert into the Zarr store
      bf_dataset: Bloom filter dataset handle from Zarr store
      columns_idx: Column indices

    """
    # Attach to the POSIX shared memory NumPy boolean array
    shared_array = sarray.attach(POSIX_SHARED_MEMORY_ARRAY)

    # Write row chunks into the Zarr store
    for row in rows:
        bf_dataset[
            row.start : row.stop, columns_idx[0] : columns_idx[1]
        ] = shared_array[row.start : row.stop, :]


def zarr_store_get_metadata(zarr_store_dir):
    """Extracts the Metadata dataset from the Zarr store

    Args:
      zarr_store_dir: Zarr store directory

    Returns: Extracted metadata dataset from the Zarr store

    """
    with zarr.open(zarr_store_dir, mode="r") as zarr_store:
        metadata = zarr_store.get(METADATA_DATASET_NAME)[0]
    return list((idx, name) for idx, name in enumerate(metadata))


def zarr_store_get_config_attributes(zarr_store_dir):
    """Extract the config attributes stored in Metadata dataset from the Zarr store

    Args:
      zarr_store_dir: Zarr store directory

    Returns: Extracted config from the metadata dataset from the Zarr store directory

    """
    config_data = list()
    with zarr.open(zarr_store_dir, mode="r") as zarr_store:
        metadata = zarr_store.get(METADATA_DATASET_NAME)
        config_data.append(("k", metadata.attrs["kmer_size"]))
        config_data.append(("m", metadata.attrs["bloom_filter_size"]))
        config_data.append(("h", metadata.attrs["number_of_hash_functions"]))
        config_data.append(("zarr_store_created_on", metadata.attrs["created_on"]))
        config_data.append(("db_built_on", datetime.today().strftime("%d-%m-%Y")))
        config_data.append(("num_updates", metadata.attrs["update_count"]))
        config_data.append(("last_updated_on", metadata.attrs["last_updated_on"]))
        config_data.append(("num_samples", metadata.attrs["number_of_samples"]))
        config_data.append(("sequence_type", metadata.attrs["sequence_type"]))

    return config_data


def zarr_store_get_chunksize(zarr_store_dir):
    """Extracts the chunk size metadata from the Zarr store

    Args:
      zarr_store_dir: Zarr store directory

    Returns: Row chunk size used in the zarr store

    """
    with zarr.open(zarr_store_dir, mode="r") as zarr_store:
        metadata = zarr_store.get(METADATA_DATASET_NAME)
        chunk_size_rows = metadata.attrs["chunk_size_rows"]
    return chunk_size_rows
