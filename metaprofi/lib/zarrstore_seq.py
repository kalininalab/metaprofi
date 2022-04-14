"""MetaProFi zarr store Bloom Filter (BF) matrix sequence-level module
"""

import os
import sys
from datetime import datetime
from multiprocessing import Process
import numpy as np
import SharedArray as sarray
import zarr
from metaprofi.lib.bloomfilter_cython import seq_bloomfilter_cython
from metaprofi.lib.constants import (
    BLOOMFILTER_DATASET_NAME,
    COMPRESSION_LEVEL,
    METADATA_DATASET_NAME,
)
from metaprofi.lib.lmdb_faq_index import LMDBFAQStore
from metaprofi.lib.utilities import (
    calculate_chunksize,
    calculate_index_chunksize,
    check_dir,
    check_zarr_index_store,
    chunks,
    cleanup_dir,
    splitjobs,
)
from numcodecs import Zstd, blosc
from tqdm import tqdm

# Global
POSIX_SHARED_MEMORY_ARRAY = ""


class ZarrStoreSeq:
    """Zarr BF matrix store

    Args:
      input_file: FASTA/FASTQ (.GZ) input file contianing one or more sequences
      config: Config dictionary containing all the necessary configuration
      update: Boolean indicating whether this is run for an update index session or build session.

    NOTE:
      1. This is used at sequence-level meaning every sequence in the input FASTA/FASTQ (.GZ) file will be considered as a sample, so 1 Bloom Filter (BF) per sequence will be created in the BF matrix.
      2. A intermediate FASTA/FASTQ index file will be created (in-terms of storage may cost at least half the size of the original input FASTA/FASTQ (.GZ) file) and later will be cleaned automatically

    """

    def __init__(self, input_file, config, update=False):
        self.config = config
        self.input_file = os.path.abspath(input_file)
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

        # Remove any previously built input file index if present
        lmdb_faq_folder = f"{self.config['output_directory']}/{os.path.splitext(os.path.basename(self.input_file))[0]}.lmdb"
        cleanup_dir(lmdb_faq_folder)

        # Building the input file index
        print("INFO: Building index for the input file ...")
        self.lmdb_faq_handle = LMDBFAQStore(lmdb_faq_folder, self.config["k"])
        self.lmdb_faq_handle.build_index(self.input_file)
        self.lmdb_faq_handle.flush()
        self.lmdb_faq_handle.close()

        # Re-open LMDB index in readonly mode
        self.lmdb_faq_handle = LMDBFAQStore(
            lmdb_faq_folder,
            self.config["k"],
            readonly=True,
            lock=False,
            max_readers=self.config["nproc"],
        )

        # Get the total number of sequences that has at least 1 k-mer
        self.number_of_seqs = len(self.lmdb_faq_handle)

        # Calculate chunksize
        (
            self.rows,
            self.columns,
            self.chunk_size_rows,
            self.chunk_size_samples,
        ) = calculate_chunksize(config, self.number_of_seqs)

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
            shape=(self.packed_bf_size, self.number_of_seqs),
            chunks=(self.rows, self.columns),
            compressor=compressor,
            synchronizer=synchronizer,
            dtype="uint8",
        )
        metadata = bf_group.create_dataset(
            name=METADATA_DATASET_NAME,
            shape=(self.number_of_seqs),
            chunks=(50_000),
            compressor=compressor,
            dtype="str",
        )
        metadata.attrs["created_on"] = datetime.today().strftime("%d-%m-%Y")
        metadata.attrs["last_updated_on"] = datetime.today().strftime("%d-%m-%Y")
        metadata.attrs["update_count"] = 0
        metadata.attrs["kmer_size"] = self.config["k"]
        metadata.attrs["bloom_filter_size"] = self.config["m"]
        metadata.attrs["packed_bf_size"] = self.packed_bf_size
        metadata.attrs["number_of_zeros_padded"] = self.number_of_zeros_to_pad
        metadata.attrs["number_of_hash_functions"] = self.config["h"]
        metadata.attrs["number_of_samples"] = self.number_of_seqs
        metadata.attrs["compressor"] = str(bf_dataset.compressor)
        metadata.attrs["chunk_size_rows"] = self.rows
        metadata.attrs["chunk_size_samples"] = self.columns
        metadata.attrs["sequence_type"] = self.config["sequence_type"]
        metadata.attrs["endianness"] = sys.byteorder

        # Index store chunk size
        index_chunk_rows = calculate_index_chunksize(self.config)
        metadata.attrs["index_chunk_rows"] = index_chunk_rows

        # Insert sample identifiers into the zarr store
        write_sample_identifiers_to_zarr(self.lmdb_faq_handle, metadata)

        # Insert the bloomfilter data into the zarr store
        self.insert_bf_into_zarr(bf_dataset)

        # Clean up
        cleanup_dir(sync_dir)
        self.lmdb_faq_handle.close()
        self.lmdb_faq_handle.cleanup()

        # Log output
        print("Bloom filter matrix store is now created with all the samples")

    def insert_bf_into_zarr(self, bf_dataset):
        """Create and insert the BloomFilters as columns into the Zarr store

        Args:
          samples: Samples path list
          bf_dataset: Bloom filter matrix dataset handle from Zarr store
          resize: Whether to resize the Zarr store dataset (Default value = False)

        Returns: None

        """
        # Progress bar object
        pbar = tqdm(
            total=self.number_of_seqs,
            bar_format="{l_bar}{bar}{n_fmt}/{total_fmt}",
            desc="Number of sequences processed",
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
        for chunked_samples in chunks(
            range(self.number_of_seqs), self.chunk_size_samples
        ):
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

            # Build bloom filters in parallel
            bf_processes = []
            for samp_list in splitjobs(
                list(zip(chunked_samples, range(len(chunked_samples)))),
                self.config["nproc"],
            ):
                proc = Process(
                    target=mp_bloomfilter,
                    args=(
                        samp_list,
                        self.config,
                        self.lmdb_faq_handle,
                    ),
                )
                proc.start()
                bf_processes.append(proc)

            for proc in bf_processes:
                proc.join()

            # Write the bloom filter matrix into the Zarr store
            # Chunk rows so chunks can be inserted parallely
            row_list = list(chunks(range(self.packed_bf_size), self.rows))

            # Column index start and stop for correct insert into the Zarr store
            col_idx = (chunked_samples.start, chunked_samples.stop)

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


def mp_bloomfilter(samples_list, config, lmdb_faq_handle):
    """Helper function for creating bloom filters in parallel

    Args:
      samples_list: Samples list
      config: Config dictionary
      lmdb_faq_handle: LMDB handle for the faq database

    """
    # Attach to the POSIX shared memory NumPy boolean array
    shared_array = sarray.attach(POSIX_SHARED_MEMORY_ARRAY)

    # Bloom filters
    for sample_info in samples_list:
        _, seq, _ = lmdb_faq_handle[sample_info[0]]
        seq_bloomfilter_cython(seq, config, shared_array, sample_info[1])
        del seq


def write_sample_identifiers_to_zarr(lmdb_faq_handle, metadata_dataset):
    """Function for inserting sample identifiers into the zarr store

    Args:
      lmdb_faq_handle: LMDB handle for the faq database
      metadata_dataset: BF matrix store's metadata dataset handle

    """
    keys = []
    bpoint = 0

    num_keys = len(lmdb_faq_handle)
    if num_keys > 5_000_000:
        bpoint = 5_000_000
    else:
        bpoint = num_keys

    start_idx = 0
    for idx, key in enumerate(lmdb_faq_handle.get_seq_names(), start=1):
        keys.append(key)
        if len(keys) == bpoint:
            metadata_dataset[start_idx:idx] = np.array(keys, dtype="object")
            keys = []
            start_idx = idx
        elif num_keys == idx:
            metadata_dataset[start_idx:idx] = np.array(keys, dtype="object")
            keys = []


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
