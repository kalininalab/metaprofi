"""MetaProFi zarr store based index module
"""

import os
import sys
import warnings
from datetime import datetime
from multiprocessing import Process
import numpy as np
import SharedArray as sarray
import zarr
from metaprofi.lib.constants import (
    BLOOMFILTER_DATASET_NAME,
    BLOOMFILTER_INDEX_DATASET_NAME,
    COMPRESSION_ALGO,
    COMPRESSION_LEVEL,
    METADATA_DATASET_NAME,
)
from metaprofi.lib.utilities import (
    check_dir,
    check_zarr_index_store,
    check_zarr_store,
    chunks,
    cleanup_dir,
    splitjobs,
    zstd_compress,
    zstd_compress_cat,
)
from numcodecs import Blosc, Zstd
from tqdm import tqdm

# Global
POSIX_SHARED_MEMORY_ARRAY = ""


def build_index(config):
    """Read bloom filter matrix from zarr store and create index store

    Args:
      config: Config dictionary containing all the necessary configuration

    """
    matrix_store_dir = config["matrix_store_name"]
    index_store_dir = config["index_store_name"]

    global POSIX_SHARED_MEMORY_ARRAY
    POSIX_SHARED_MEMORY_ARRAY = config["POSIX_SHARED_MEMORY_ARRAY"]

    # Check if zarr store contains required datasets
    check_zarr_store(matrix_store_dir)

    # Check if zarr store already exists
    if check_dir(index_store_dir):
        raise FileExistsError(
            f"Index zarr store '{index_store_dir}' already exists and build related command is used. Please use appropriate MetaProFi command or delete the existing index store."
        )

    # Zarr store related
    Blosc.use_threads = False
    # compressor = Blosc(cname=COMPRESSION_ALGO, clevel=COMPRESSION_LEVEL)
    compressor = Zstd(level=COMPRESSION_LEVEL)

    # Matrix store
    matrix_store = zarr.open(matrix_store_dir, mode="r")

    # Dataset
    bf_dataset = matrix_store.get(BLOOMFILTER_DATASET_NAME)
    metadata = matrix_store.get(METADATA_DATASET_NAME)
    number_of_samples = metadata.attrs["number_of_samples"]
    chunk_size_rows = int(bf_dataset.chunks[0] // 1.5)
    chunk_size_samples = bf_dataset.chunks[1]
    index_chunk_rows = metadata.attrs["index_chunk_rows"]
    packed_bf_size = metadata.attrs["packed_bf_size"]
    number_of_zeros_padded = metadata.attrs["number_of_zeros_padded"]

    # Index store
    index_sync_dir = os.path.splitext(index_store_dir)[0] + ".sync"
    index_synchronizer = zarr.ProcessSynchronizer(index_sync_dir)
    index_store = zarr.NestedDirectoryStore(index_store_dir)
    index_group = zarr.group(store=index_store)

    # Zarr currently creates a FutureWarning when copy or copy_all functions
    # are used and this is safe to ignore, for further information check,
    # https://github.com/zarr-developers/zarr-python/issues/257
    with warnings.catch_warnings():
        warnings.simplefilter(action="ignore", category=FutureWarning)
        zarr.copy(metadata, index_group, without_attrs=False)

    # Index store datasets
    index_dataset = index_group.create_dataset(
        name=BLOOMFILTER_INDEX_DATASET_NAME,
        shape=(config["m"],),
        chunks=(index_chunk_rows,),
        compressor=compressor,
        synchronizer=index_synchronizer,
        dtype=bytes,
    )

    # Progress bar object
    pbar = tqdm(
        total=config["m"],
        bar_format="{l_bar}{bar}{n_fmt}/{total_fmt}",
        desc="Indexing (Row/Bit-slices)",
        leave=True,
        ncols=80,
        colour="green",
    )

    # Create shared array
    shared_array = sarray.create(
        POSIX_SHARED_MEMORY_ARRAY, (chunk_size_rows * 8, number_of_samples), "uint8"
    )
    shared_array[:] = False

    for row_chunk, mrow_chunk in zip(
        chunks(range(packed_bf_size), chunk_size_rows),
        chunks(range(config["m"] + number_of_zeros_padded), chunk_size_rows * 8),
    ):
        # shared_array reset and re-initialization
        if len(row_chunk) != chunk_size_rows:
            # Delete the existing shared array
            del shared_array
            sarray.delete(POSIX_SHARED_MEMORY_ARRAY)

            # Re-create shared array
            shared_array = sarray.create(
                POSIX_SHARED_MEMORY_ARRAY,
                (len(row_chunk) * 8, number_of_samples),
                "uint8",
            )

            # Initialize shared array
            shared_array[:] = False

        # Create a columns index list to extract based on the sample chunksize
        columns_list = list(
            chunks(range(number_of_samples), chunk_size_samples),
        )

        # Chunk rows so chunks can be read in parallel
        # read_row_list = list(
        #     zip(
        #         chunks(range(len(row_chunk)), chunk_size_rows),
        #         chunks(row_chunk, chunk_size_rows),
        #     )
        # )

        # Extract rows from Zarr store
        processes = []
        for sub_col in splitjobs(columns_list, config["nproc"]):
            proc = Process(
                target=extract_rows_from_zarr,
                args=(bf_dataset, row_chunk, sub_col),
            )
            proc.start()
            processes.append(proc)

        for proc in processes:
            proc.join()

        # Write the rows (encode and compress) from matrix into the Zarr index store
        # Chunk rows so chunks can be inserted parallely
        if mrow_chunk.stop == config["m"] + number_of_zeros_padded:
            row_list = list(
                zip(
                    chunks(
                        range(len(mrow_chunk) - number_of_zeros_padded),
                        index_chunk_rows,
                    ),
                    chunks(
                        range(
                            mrow_chunk.start, mrow_chunk.stop - number_of_zeros_padded
                        ),
                        index_chunk_rows,
                    ),
                )
            )
        else:
            row_list = list(
                zip(
                    chunks(range(len(mrow_chunk)), index_chunk_rows),
                    chunks(mrow_chunk, index_chunk_rows),
                )
            )

        # Encode and compress every row and insert into the index store
        zarr_processes = []
        for sub_row in splitjobs(row_list, config["nproc"]):
            zarr_proc = Process(
                target=write_to_zarr_index_store,
                args=(index_dataset, sub_row, False),
            )
            zarr_proc.start()
            zarr_processes.append(zarr_proc)

        for zarr_proc in zarr_processes:
            zarr_proc.join()

        # Update progress bar
        if mrow_chunk.stop == config["m"] + number_of_zeros_padded:
            pbar.update(len(mrow_chunk) - number_of_zeros_padded)
        else:
            pbar.update(len(mrow_chunk))

    # Close progress bar
    pbar.close()

    # Clean up
    del shared_array
    sarray.delete(POSIX_SHARED_MEMORY_ARRAY)
    cleanup_dir(index_sync_dir)

    # Log output
    print("Index store is ready for querying/searching")


def update_index(config):
    """Update index store with new samples bloom filter data

    Args:
      config: Config dictionary containing all the necessary configuration

    """
    matrix_store_dir = config["matrix_store_name"]
    index_store_dir = config["index_store_name"]

    global POSIX_SHARED_MEMORY_ARRAY
    POSIX_SHARED_MEMORY_ARRAY = config["POSIX_SHARED_MEMORY_ARRAY"]

    # Check if zarr store(s) contains required datasets
    check_zarr_store(matrix_store_dir)
    check_zarr_index_store(index_store_dir)

    # Zarr store related
    Blosc.use_threads = False

    # Matrix store
    matrix_store = zarr.open(matrix_store_dir, mode="r")

    # Matrix store datasets
    bf_dataset = matrix_store.get(BLOOMFILTER_DATASET_NAME)
    metadata = matrix_store.get(METADATA_DATASET_NAME)
    matrix_store_number_of_samples = metadata.attrs["number_of_samples"]
    packed_bf_size = metadata.attrs["packed_bf_size"]
    number_of_zeros_padded = metadata.attrs["number_of_zeros_padded"]
    chunk_size_rows = int(bf_dataset.chunks[0] // 1.5)
    chunk_size_samples = bf_dataset.chunks[1]
    bf_endianness = metadata.attrs["endianness"]

    # Index store
    index_sync_dir = os.path.splitext(index_store_dir)[0] + ".sync"
    index_synchronizer = zarr.ProcessSynchronizer(index_sync_dir)
    index_store = zarr.open(index_store_dir, synchronizer=index_synchronizer)

    # Index store datasets
    index_dataset = index_store.get(BLOOMFILTER_INDEX_DATASET_NAME)
    index_metadata_dataset = index_store.get(METADATA_DATASET_NAME)
    index_chunk_rows = index_metadata_dataset.attrs["index_chunk_rows"]
    index_endianness = index_metadata_dataset.attrs["endianness"]

    # Check the metadata in the matrix store matches the metadata in the index store
    matrix_metadata_attrs = [
        metadata.attrs["kmer_size"],
        metadata.attrs["bloom_filter_size"],
        metadata.attrs["number_of_hash_functions"],
        metadata.attrs["sequence_type"],
    ]
    index_metadata_attrs = [
        index_metadata_dataset.attrs["kmer_size"],
        index_metadata_dataset.attrs["bloom_filter_size"],
        index_metadata_dataset.attrs["number_of_hash_functions"],
        index_metadata_dataset.attrs["sequence_type"],
    ]

    if matrix_metadata_attrs != index_metadata_attrs:
        raise ValueError(
            f"Can't perform the index store update. The configuration used for creating the index store does not match with the configuration used for the matrix store \nConfig from matrix store:\n\tk:{matrix_metadata_attrs[0]}, m:{matrix_metadata_attrs[1]}, h:{matrix_metadata_attrs[2]}, sequence_type:{matrix_metadata_attrs[3]}\nConfig from index store:\n\tk:{index_metadata_attrs[0]}, m:{index_metadata_attrs[1]}, h:{index_metadata_attrs[2]}, sequence_type:{index_metadata_attrs[3]}"
        )

    if bf_endianness != index_endianness:
        raise ValueError(
            f"Index cannot be updated. Bloom filter matrix endianness ({bf_endianness} endian) and index endianness ({index_endianness} endian) does not match."
        )

    # Update/append metadata in Metadata dataset and update certain attributes
    row_idx = index_metadata_dataset.shape[0]
    index_metadata_dataset.resize(row_idx + metadata.shape[0])
    index_metadata_dataset[row_idx : row_idx + metadata.shape[0]] = metadata[:]

    # Update metadata attributes
    index_metadata_dataset.attrs["number_of_samples"] = (
        index_metadata_dataset.attrs["number_of_samples"]
        + metadata.attrs["number_of_samples"]
    )
    index_metadata_dataset.attrs["last_updated_on"] = datetime.today().strftime(
        "%d-%m-%Y"
    )
    index_metadata_dataset.attrs["update_count"] = (
        index_metadata_dataset.attrs["update_count"] + 1
    )
    index_metadata_dataset.attrs["chunk_size_rows"] = metadata.attrs["chunk_size_rows"]

    # New number of samples to be updated/added to the index store
    number_of_new_samples = matrix_store_number_of_samples

    # Progress bar object
    pbar = tqdm(
        total=config["m"],
        bar_format="{l_bar}{bar}{n_fmt}/{total_fmt}",
        desc="Updating index (Row/Bit-slices)",
        leave=True,
        ncols=80,
        colour="green",
    )

    # # Rather than wasting the user allocated memory we try to read in as many rows as possible
    # bytes_per_bloom_filter = config["m"]
    # total_required_bytes = bytes_per_bloom_filter * number_of_new_samples
    # num_chunks = math.ceil(
    #     total_required_bytes
    #     / (parse_size(config["max_memory"], binary=True) * (85 / 100))
    # )
    # possible_chunk_size_rows = math.floor(config["m"] / num_chunks)

    # Create shared array
    shared_array = sarray.create(
        POSIX_SHARED_MEMORY_ARRAY, (chunk_size_rows * 8, number_of_new_samples), "uint8"
    )
    shared_array[:] = False

    for row_chunk, mrow_chunk in zip(
        chunks(range(packed_bf_size), chunk_size_rows),
        chunks(range(config["m"] + number_of_zeros_padded), chunk_size_rows * 8),
    ):
        # shared_array reset and re-initialization
        if len(row_chunk) != chunk_size_rows:
            # Delete the existing shared array
            del shared_array
            sarray.delete(POSIX_SHARED_MEMORY_ARRAY)

            # Re-create shared array
            shared_array = sarray.create(
                POSIX_SHARED_MEMORY_ARRAY,
                (len(row_chunk) * 8, number_of_new_samples),
                "uint8",
            )

            # Initialize shared array
            shared_array[:] = False

        # Create a columns index list to extract based on the sample chunksize
        columns_list = list(
            chunks(range(number_of_new_samples), chunk_size_samples),
        )

        # Chunk rows so chunks can be read in parallel
        # read_row_list = list(
        #     zip(
        #         chunks(range(len(row_chunk)), chunk_size_rows),
        #         chunks(row_chunk, chunk_size_rows),
        #     )
        # )

        # Extract rows from Zarr store
        processes = []
        for sub_col in splitjobs(columns_list, config["nproc"]):
            proc = Process(
                target=extract_rows_from_zarr,
                args=(bf_dataset, row_chunk, sub_col),
            )
            proc.start()
            processes.append(proc)

        for proc in processes:
            proc.join()

        # Write the rows (encode and compress) from matrix into the Zarr index store
        # Chunk rows so chunks can be inserted parallely
        if mrow_chunk.stop == config["m"] + number_of_zeros_padded:
            row_list = list(
                zip(
                    chunks(
                        range(len(mrow_chunk) - number_of_zeros_padded),
                        index_chunk_rows,
                    ),
                    chunks(
                        range(
                            mrow_chunk.start, mrow_chunk.stop - number_of_zeros_padded
                        ),
                        index_chunk_rows,
                    ),
                )
            )
        else:
            row_list = list(
                zip(
                    chunks(range(len(mrow_chunk)), index_chunk_rows),
                    chunks(mrow_chunk, index_chunk_rows),
                )
            )

        # Encode and compress every row and insert into the index store
        zarr_processes = []
        for sub_row in splitjobs(row_list, config["nproc"]):
            zarr_proc = Process(
                target=write_to_zarr_index_store,
                args=(index_dataset, sub_row, True),
            )
            zarr_proc.start()
            zarr_processes.append(zarr_proc)

        for zarr_proc in zarr_processes:
            zarr_proc.join()

        # Update progress bar
        if mrow_chunk.stop == config["m"] + number_of_zeros_padded:
            pbar.update(len(mrow_chunk) - number_of_zeros_padded)
        else:
            pbar.update(len(mrow_chunk))

    # Close progress bar
    pbar.close()

    # Clean up
    del shared_array
    sarray.delete(POSIX_SHARED_MEMORY_ARRAY)
    cleanup_dir(index_sync_dir)

    # Log output
    print("Index store is updated with new samples and ready for querying/searching")


def extract_rows_from_zarr(bf_dataset, row_chunk, columns_list):
    """Extract rows chunkwise from Zarr store and store it in the shared array

    Args:
      bf_dataset: Bloom filter dataset object handle from Zarr store
      read_row_list: Row indices to extract data from Zarr store
      columns_list: Column indices to extract data from Zarr store

    """
    # Attach to the POSIX shared memory NumPy boolean array
    shared_array = sarray.attach(POSIX_SHARED_MEMORY_ARRAY)

    # for shared_array_idx, zarr_idx in read_row_list:
    #     for sub_col in columns_list:
    #         # shared_array[
    #         #     shared_array_idx.start : shared_array_idx.stop,
    #         #     sub_col.start : sub_col.stop,
    #         # ] = bf_dataset[zarr_idx.start : zarr_idx.stop, sub_col.start : sub_col.stop]
    #         shared_array[:, sub_col.start : sub_col.stop] = np.unpackbits(
    #             bf_dataset[
    #                 zarr_idx.start : zarr_idx.stop, sub_col.start : sub_col.stop
    #             ],
    #             axis=0,
    #         )

    for sub_col in columns_list:
        shared_array[:, sub_col.start : sub_col.stop] = np.unpackbits(
            bf_dataset[row_chunk.start : row_chunk.stop, sub_col.start : sub_col.stop],
            axis=0,
            bitorder="little",
        )


def write_to_zarr_index_store(index_dataset, row_chunks, update=False):
    """Set multiple key, value pairs at once

    Args:
        index_dataset: Index dataset object handle from Zarr store
        row_chunks: List of row indices to write to the index store chunkwise

    """
    # Attach to the POSIX shared memory NumPy boolean array
    shared_array = sarray.attach(POSIX_SHARED_MEMORY_ARRAY)

    if not update:
        for sub_chunk, zarr_idx in row_chunks:
            np_obj_array = np.zeros(len(sub_chunk), dtype=object)
            for idx, obj in enumerate(
                map(zstd_compress, shared_array[sub_chunk.start : sub_chunk.stop])
            ):
                np_obj_array[idx] = obj
            index_dataset[zarr_idx.start : zarr_idx.stop] = np_obj_array
    else:
        for sub_chunk, zarr_idx in row_chunks:
            np_obj_array = np.zeros(len(sub_chunk), dtype=object)
            for idx, obj in enumerate(
                map(
                    zstd_compress_cat,
                    index_dataset[zarr_idx.start : zarr_idx.stop],
                    shared_array[sub_chunk.start : sub_chunk.stop],
                )
            ):
                np_obj_array[idx] = obj
            index_dataset[zarr_idx.start : zarr_idx.stop] = np_obj_array
