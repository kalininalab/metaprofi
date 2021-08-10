"""MetaProFi index search module
"""

import itertools
import math
import os
import sys
from datetime import datetime
from multiprocessing import Manager, Process
from operator import itemgetter
import numpy as np
import zarr
from metaprofi.lib.bloomfilter_cython import sequence_to_hash_cython
from metaprofi.lib.constants import (
    BLOOMFILTER_INDEX_DATASET_NAME,
    METADATA_DATASET_NAME,
)
from metaprofi.lib.translate_fasta import six_frame_translation
from metaprofi.lib.utilities import (
    bitwise_and,
    check_seq,
    check_zarr_index_store,
    cleanup_dir,
    cleanup_file,
    splitjobs,
    zstd_decompress,
)
from metaprofi.lib.lmdb_index import LMDBStore


def search_index(config, query, inp_seq, seq_type, threshold=100, write_to_file=True):
    """Given a query and a threshold, will search for the samples
    containing that sequence/K-mers and returns the sample names when found

    Args:
      config: Config dictionary containing all the necessary configuration
      query: Query sequence or Fasta/Fastq file containing one or more sequences (could be nucleotides or aminoacids)
      inp_seq: If input query is sequence, set it to True
      seq_type: Type of the sequence. Either nucleotides or aminoacids
      threshold: Threshold value to perform exact or approximate search (Default value = 100)
      write_to_file: Boolean whether to write results to a file or return (default: True)

    Returns: A file containing results of the sample names/ids in which the query sequence was found or returns the results

    """
    index_store_dir = config["index_store_name"]

    # Not valid
    if seq_type == "aminoacid" and config["sequence_type"] == "nucleotide":
        raise ValueError(
            "Given input sequence is of type aminoacid but the index store is of type nucleotide"
        )

    # Check if index store contains required datasets
    check_zarr_index_store(index_store_dir)

    # Check if the output_directory is set in the config file
    output_dir = config["output_directory"]

    # Index store
    index_store = zarr.open(index_store_dir, mode="r")

    # Index store datasets
    index_dataset = index_store.get(BLOOMFILTER_INDEX_DATASET_NAME)
    index_metadata_dataset = index_store.get(METADATA_DATASET_NAME)
    index_endianness = index_metadata_dataset.attrs["endianness"]
    current_endianness = sys.byteorder

    # Validate config with config stored in index store
    index_store_config = [
        index_metadata_dataset.attrs["kmer_size"],
        index_metadata_dataset.attrs["number_of_hash_functions"],
        index_metadata_dataset.attrs["bloom_filter_size"],
        index_metadata_dataset.attrs["sequence_type"],
    ]

    if (
        config["k"] != index_store_config[0]
        or config["h"] != index_store_config[1]
        or config["m"] != index_store_config[2]
        or config["sequence_type"] != index_store_config[3]
    ):
        raise ValueError(
            f"K-mer size (k) and/or number of hash functions (h) and/or bloom filter size (m) and/or sequence type (sequence_type) in the input config file does not match the config used for building the index store. Please check.\nConfig used for building index store\n\tk: {index_store_config[0]}, h: {index_store_config[1]}, m: {index_store_config[2]}, sequence_type: {index_store_config[3]}"
        )

    if index_endianness != current_endianness:
        raise ValueError(
            f"Index cannot be queried. Index endianness ({index_endianness} endian) does not match the current hardware endianness ({current_endianness} endian)."
        )

    # Final results container
    all_results = {}
    orfs_filename = ""

    # Scenario 1: One sequence only via command line input
    if inp_seq:
        # When input sequence is of type nucleotide, and the index store is of type aminoacid, we do a six frame translation and use them for searching/querying against index store
        query = query.upper()
        if seq_type == "nucleotide" and config["sequence_type"] == "aminoacid":
            seq_type = "aminoacid"
            # Validate the sequences are actually nucleotides
            if check_seq(query) != "nucleotide":
                raise ValueError(
                    f"Given input sequence does not match the sequence type: '{seq_type}'"
                )

            # Six frame translation
            open_reading_frames = six_frame_translation(
                query, "user_input_seq", config["k"]
            )

            if not open_reading_frames:
                print(
                    f"Given nucleotide sequence does not contain any amino acid sequence with size >= {config['k']} (K-mer size)"
                )
                sys.exit()

            # Do index search/query
            for key, val in open_reading_frames.items():
                # hash_kmers = seq_to_kmers_to_hash_array(val, config, seq_type)
                hash_kmers = sequence_to_hash_cython(val, seq_type, config)
                all_results[key] = sequence_query(
                    hash_kmers, index_dataset, index_metadata_dataset, threshold
                )

            # Output processing
            query_results = process_results(all_results)

            if write_to_file:
                write_query_results(query_results, threshold, output_dir)
                return
            else:
                return query_results, threshold
        else:
            # Validate the sequence
            if seq_type != check_seq(query):
                raise ValueError(
                    f"Given input sequence does not match the given sequence type: '{seq_type}'"
                )

            # hash_kmers = seq_to_kmers_to_hash_array(query, config, seq_type)
            hash_kmers = sequence_to_hash_cython(query, seq_type, config)
            results = sequence_query(
                hash_kmers, index_dataset, index_metadata_dataset, threshold
            )

            if write_to_file:
                write_query_results(results, threshold, output_dir)
                return
            else:
                return results, threshold

    # Scenario 2: One Fasta/Fastq file with 1 or more sequence(s)
    # Create a LMDB index for the query file
    query = os.path.abspath(query)
    lmdb_folder = f"{config['output_directory']}/{os.path.splitext(os.path.basename(query))[0]}.lmdb"

    # Remove the index if it already exists
    cleanup_dir(lmdb_folder)

    # Create the index
    lmdb_handle = LMDBStore(lmdb_folder, config["k"])
    lmdb_handle.build_index(query)
    lmdb_handle.flush()
    lmdb_handle.close()

    # Re-open LMDB index in readonly mode (for parallel reads)
    lmdb_handle = LMDBStore(
        lmdb_folder,
        config["k"],
        readonly=True,
        lock=False,
        max_readers=config["nproc"],
    )

    # If input fasta/fastq file is of type nucleotide, and the index store is of type aminoacid we do a 6 frame translation to aminoacids and use those sequences to do the index search/query
    if seq_type == "nucleotide" and config["sequence_type"] == "aminoacid":
        seq_type = "aminoacid"
        no_orfs_seq_names = []
        orfs_filename = f"{os.path.splitext(query)[0]}_orfs_translated.fasta"

        # If the tranlated seqs file already exists (e.g, from previous runs) we just remove it
        if os.path.isfile(orfs_filename):
            cleanup_file(orfs_filename)

        # Validate the sequences are actually nucleotides
        if lmdb_handle.type != "nucleotide":
            lmdb_handle.close()
            lmdb_handle.cleanup()
            raise ValueError(
                f"Given input fasta/fastq file sequences do not match the given sequence type: '{seq_type}'"
            )

        # Six frame translation and writing to an uncompressed fasta file
        num_seqs = lmdb_handle.__len__()
        for i in range(0, num_seqs):
            name, seq, _ = lmdb_handle[i]
            translated = six_frame_translation(
                seq,
                name,
                config["k"],
                orfs_filename,
                out_file=True,
            )

            # Keep the name of the seqs which did not have any K-mers/ORFs
            if not translated:
                no_orfs_seq_names.append(name)

        # Cleanup fasta/fastq index and create an index for the ORF translated seq file
        lmdb_handle.close()
        lmdb_handle.cleanup()

        # If none of the sequence in the input file has a K-mer
        if len(no_orfs_seq_names) == num_seqs:
            print(
                "None of the sequence in the input file contains amino acid sequence (after six-frame translation) with size >= {config['k']} (K-mer size)"
            )
            sys.exit()

        query = os.path.abspath(orfs_filename)
        lmdb_folder = f"{config['output_directory']}/{os.path.splitext(os.path.basename(query))[0]}.lmdb"

        # Remove the index if it already exists
        cleanup_dir(lmdb_folder)

        # Create the index
        lmdb_handle = LMDBStore(lmdb_folder, config["k"])
        lmdb_handle.build_index(query)
        lmdb_handle.flush()
        lmdb_handle.close()

        # Re-open LMDB index in readonly mode
        lmdb_handle = LMDBStore(
            lmdb_folder,
            config["k"],
            readonly=True,
            lock=False,
            max_readers=config["nproc"],
        )

    # Validate the sequences
    if (
        seq_type == "nucleotide"
        and lmdb_handle.type != "nucleotide"
        or seq_type == "aminoacid"
        and lmdb_handle.type == "nucleotide"
    ):
        lmdb_handle.close()
        lmdb_handle.cleanup()
        raise ValueError(
            f"Given input fasta/fastq file sequences does not match the given sequence type: '{seq_type}'"
        )

    num_sequences = lmdb_handle.__len__()

    if num_sequences == 0:
        lmdb_handle.close()
        lmdb_handle.cleanup()
        raise ValueError("Input query file is empty")

    # Setup lock and queue for collecting results of every sequence
    lock = Manager().Lock()
    out_queue = Manager().Queue()
    procs = []

    # Calculate hash for every kmer in every sequence
    for seq_indices in splitjobs(range(num_sequences), config["nproc"]):
        proc = Process(
            target=file_query,
            args=(
                lmdb_handle,
                seq_indices,
                config,
                seq_type,
                threshold,
                index_dataset,
                index_metadata_dataset,
                lock,
                out_queue,
            ),
        )
        procs.append(proc)
        proc.start()

    for proc in procs:
        proc.join()

    # Dequeue the results and make a dictionary
    while not out_queue.empty():
        for seq_id, sample_ids in out_queue.get():
            all_results[seq_id] = sample_ids

    # Process results
    query_results = process_results(all_results)

    # Cleanup
    if orfs_filename:
        cleanup_file(orfs_filename)
    lmdb_handle.close()
    lmdb_handle.cleanup()

    # Write results to a file
    if write_to_file:
        write_query_results(query_results, threshold, output_dir)
        return
    else:
        return query_results, threshold


def sequence_query(kmer_hash_array, index_dataset, index_metadata_dataset, threshold):
    """Get sample identifiers for sequence queries

    Args:
      kmer_hash_array: Hash array of the query sequence
      index_dataset: Index dataset object handle from Zarr store
      index_metadata_dataset: Index store metadata dataset handle
      threshold: Threshold to perfrom exact or approxiamte search
      seq_id: Sequence identifier

    Returns: Sample names/ids in which the query sequence was found

    """
    kmer_rows = {}
    num_kmers = len(kmer_hash_array)

    # Retrieve the rows from the index store and decode the
    # encoded bitarray and perform a bitwise and operation for
    # each kmer
    for key, val in enumerate(kmer_hash_array):
        kmer_rows[key] = bitwise_and(
            [zstd_decompress(row) for row in index_dataset.vindex[val]]
        )

    # Search
    if threshold != 100:
        results = approx_search(num_kmers, kmer_rows, threshold, index_metadata_dataset)
    else:
        results = exact_search(kmer_rows, index_metadata_dataset)

    return results


def file_query(
    lmdb_handle,
    seq_indices,
    config,
    seq_type,
    threshold,
    index_dataset,
    index_metadata_dataset,
    lock=None,
    queue=None,
):
    """Get sample identifiers for given FASTA/FASTQ file

    Args:
      lmdb_handle: LMDB index DB handle
      seq_indices: Subsequence index list
      config: Config dictionary
      seq_type: Type of the sequence in the input file
      threshold: Threshold to perfrom exact or approxiamte search
      index_dataset: Index dataset object handle from Zarr store
      index_metadata_dataset: Index store metadata dataset handle
      lock: Lock object to lock the Queue
      queue: Queue object to store the results

    Note:
      Writes sample identifier(s) found for every sequence in the given sequence list to the queue

    """
    results = []
    hash_arrays = []

    # Convert sequence(s) to kmers to hash array(s)
    for seq_index in seq_indices:
        name, seq, _ = lmdb_handle[seq_index]
        hash_array = sequence_to_hash_cython(seq, seq_type, config)
        hash_arrays.append([name, hash_array])

    # Deduplicate hashes and retrieve the corresponding bit-slices from the index store
    index_keys = list(
        set(itertools.chain(*itertools.chain(*map(itemgetter(1), hash_arrays))))
    )
    index_keys.sort()
    index_dict = dict(zip(index_keys, index_dataset.vindex[index_keys]))

    # For every sequence (hash array of the sequence) apply exact or approximate query search
    for seq_id, hash_array in hash_arrays:
        kmer_rows = {}
        num_kmers = len(hash_array)

        for key, val in enumerate(hash_array):
            kmer_rows[key] = bitwise_and(
                [zstd_decompress(index_dict[row]) for row in val]
            )

        # Based on the threshold apply approximate or exact search
        if threshold != 100:
            results.append(
                [
                    seq_id,
                    approx_search(
                        num_kmers, kmer_rows, threshold, index_metadata_dataset
                    ),
                ]
            )
        else:
            results.append([seq_id, exact_search(kmer_rows, index_metadata_dataset)])
        del kmer_rows

    # Write the sample identifiers found to the queue
    with lock:
        queue.put(results)
    return


def approx_search(num_kmers, kmer_rows, threshold, index_metadata_dataset):
    """Approximate search

    Args:
      num_kmers: Number of K-mers in the given query sequence
      kmer_rows: Dictionary object (keys=kmer index and values=bitarrays)
      threshold: Threshold to perfrom approxiamte search
      index_metadata_dataset: Index store metadata dataset handle

    Returns: Sample names/ids in which the query sequence was found

    """
    # Find if minimum number of expected kmers (threshold_to_min_kmers) are
    # present in any sample and then retrieve the sample names
    threshold_to_min_kmers = math.ceil((num_kmers / 100) * threshold)

    # Sum the bitarrays of all kmers, element wise
    kmers_summed = map(sum, zip(*kmer_rows.values()))
    del kmer_rows

    # Filter out the kmers that are less than the minimum number of
    # expected kmers and then get the sample names from the index store
    filtered_samples = {
        key: f"{val} ({round((val / num_kmers)*100, 2)}%)"
        for key, val in enumerate(kmers_summed, start=0)
        if val >= threshold_to_min_kmers
    }

    if not filtered_samples:
        return None

    sample_identifiers = list(
        index_metadata_dataset.vindex[list(filtered_samples.keys())]
    )

    return dict(zip(sample_identifiers, filtered_samples.values()))


def exact_search(kmer_rows, index_metadata_dataset):
    """Exact search

    Args:
      kmer_rows: Dictionary object (keys=kmer index and values=bitarrays)
      index_metadata_dataset: Index store metadata dataset handle

    Returns: Sample names/ids in which the query sequence was found

    """
    # Find if all kmers of the query are present and then retrieve
    # sample names (when threshold == 100)
    sample_idxs = bitwise_and(kmer_rows.values()).search(1)
    del kmer_rows

    if not sample_idxs:
        return None

    sample_identifiers = list(index_metadata_dataset.vindex[sample_idxs])
    return sample_identifiers


def process_results(query_results):
    """Removes empty results

    Args:
      query_results: Disctionary object containing the name(s) of the sequences and the sample idnetifiers

    Returns: A dictionary with results after clean up

    """
    results = {}
    # no_results = []
    for key, val in query_results.items():
        if val:
            results[key] = val
        # else:
        #     no_results.append(key)
    return results


def write_query_results(results, threshold, output_dir):
    """Format the output results and write it to an output file

    Args:
      results: Sample names/ids in which the query sequence was found as a dictionary or list object
      threshold: Threshold value to perform exact or approximate query search
      output_dir: Path to the output directory to store the query results

    Returns: A file containing results of the sample names/ids in which the query sequence(s) were found

    """
    outfile_name = f"{output_dir}/metaprofi_query_results-{datetime.now().strftime('%d_%m_%Y-%I_%M_%S')}_t{threshold}.txt"
    if results:
        with open(outfile_name, "w") as output_writer:
            output_writer.write(
                "The query sequence(s) were found in the following sample(s):\n"
            )
            if threshold == 100:
                if isinstance(results, dict):
                    for seq_id, names in results.items():
                        output_writer.write(f"Query sequence id: {seq_id}\n")
                        output_writer.write("\tSample identifier(s):\n")
                        for name in names:
                            output_writer.write(f"\t\t{name}\n")
                else:
                    for name in results:
                        output_writer.write(f"{name}\n")
            else:
                if any(isinstance(i, dict) for i in results.values()):
                    for seq_id, val in results.items():
                        output_writer.write(f"Query sequence id/name: {seq_id}\n")
                        output_writer.write("\tSample identifier(s):\n")
                        for names, num_kmers in val.items():
                            output_writer.write(
                                f"\t\t{names}, Number of kmers found: {num_kmers}\n"
                            )
                else:
                    output_writer.write("Sample identifier(s):\n")
                    for key, val in results.items():
                        output_writer.write(f"\t{key}, Number of K-mers found: {val}\n")

        print(f"\nMetaProFi search results can be found in the file: {outfile_name}\n")
    else:
        print(
            f"MetaProFi: Query sequence(s) were not found in any sample in the index store. Maybe reduce the threshold (current threshold = {threshold}) and search again."
        )
