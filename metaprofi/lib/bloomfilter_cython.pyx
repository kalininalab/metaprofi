#cython: language_level=3, boundscheck=False, wraparound=False
"""NumPy boolean bloom filter in cython
"""

cimport cython
cimport numpy as cnp
from libc.stdlib cimport free, malloc
from libc.stdint cimport uint8_t, uint16_t, uint32_t, uint64_t
from metaprofi.lib.utilities import get_file_reader


cdef extern from "Python.h":
    const char* PyUnicode_AsUTF8AndSize(object, Py_ssize_t *size)


def check_kmer_exists(list input_files, uint32_t kmer_size):
    """Check if a k-mer exists in any of the given input files
    
    Used for filtering out the input sample files

    Args:
      input_files: List of Fasta/Fastq nucleotide/aminoacid sequence(s) file path corresponding to the same sample id
      kmer_size: Size of the k-mer

    Returns: Boolean (k-mer exists or not)

    """
    # Declaration and initialization
    cdef tuple read_obj
    cdef str input_file
    cdef Py_ssize_t seq_len
    cdef const char *seq_cstr
    cdef bint exists = False

    # Check if at least 1 of the sequence in any of the given input files contains a k-mer
    for input_file in input_files:
        file_reader = get_file_reader(input_file)
        for read_obj in file_reader(input_file, build_index=False):
            seq_cstr = PyUnicode_AsUTF8AndSize(read_obj[1], &seq_len)
            if seq_len >= kmer_size:
                exists = True
                break
        if exists:
            break
    return exists


cdef uint8_t setbit(uint8_t np_num, int bit_idx):
    """Flips a specific bit on a given NumPy uint8 integer 

    Args:
      np_num: NumPy uint8 integer
      bit_idx: Index of the bit to be flipped after unpacking of the given uint8 [0...7]

    Note:
      Little endian bitorder

    Returns: Bit flipped integer
    """
    return np_num | (1 << bit_idx) & 255


@cython.cdivision(True)
def bloomfilter_cython(list input_files, dict config, cnp.npy_uint8[:, :] shared_array, Py_ssize_t column_idx):
    """Computes hash for every k-mer in ever sequence in the input file(s) and creates the indexed bloom filter

    Args:
    input_files: List of Fasta/Fastq nucleotide/aminoacid sequence(s) file path corresponding to the same sample id
    config: Config dictionary
    shared_array: SharedArray (Boolean NumPy array)
    column_idx: Column index in the shared array to use as bloom filter for the given sample

    Returns: None

    """
    # Declaration and initialization
    cdef uint64_t kmer_size = config["k"]
    cdef uint64_t num_hash = config["h"]
    cdef Py_ssize_t bf_size = config["m"]
    cdef str sequence_type = config["sequence_type"]
    cdef tuple read_obj
    cdef Py_ssize_t i, row_idx, seq_len
    cdef int h, s
    cdef str input_file
    cdef const char *seq_cstr
    cdef char *seq_pointer
    cdef char *kmer
    cdef char *canonical
    cdef uint64_t hash_value
    cdef uint8_t bit_idx, value
    cdef char *reverse_kmer = <char *>malloc(kmer_size)
    cdef uint64_t *hashes = <uint64_t *>malloc(num_hash * sizeof(uint64_t))
    cdef uint64_t *seeds = <uint64_t *>malloc(num_hash * sizeof(uint64_t))

    for s in range(0, num_hash):
        seeds[s] = s

    # Bloom filter
    if sequence_type == "nucleotide":
        for input_file in input_files:
            file_reader = get_file_reader(input_file)
            for read_obj in file_reader(input_file, build_index=False):
                seq_cstr = PyUnicode_AsUTF8AndSize(read_obj[1], &seq_len)
                seq_pointer = &seq_cstr[0]
                if seq_len >= kmer_size:
                    for i in range(seq_len - kmer_size + 1):
                        kmer = seq_pointer + i
                        canonical = canonical_kmer(kmer, kmer_size, reverse_kmer)
                        mmh2_hash_64(canonical, kmer_size, seeds, num_hash, hashes)
                        for h in range(0, num_hash):
                            hash_value = hashes[h] % bf_size
                            bit_idx = hash_value % 8
                            row_idx = hash_value // 8
                            value = setbit(shared_array[row_idx, column_idx], bit_idx)
                            shared_array[row_idx, column_idx] = value
    elif sequence_type == "aminoacid":
        for input_file in input_files:
            file_reader = get_file_reader(input_file)
            for read_obj in file_reader(input_file, build_index=False):
                seq_cstr = PyUnicode_AsUTF8AndSize(read_obj[1], &seq_len)
                seq_pointer = &seq_cstr[0]
                if seq_len >= kmer_size:
                    for i in range(seq_len - kmer_size + 1):
                        kmer = seq_pointer + i
                        mmh2_hash_64(kmer, kmer_size, seeds, num_hash, hashes)
                        for h in range(0, num_hash):
                            hash_value = hashes[h] % bf_size
                            bit_idx = hash_value % 8
                            row_idx = hash_value // 8
                            value = setbit(shared_array[row_idx, column_idx], bit_idx)
                            shared_array[row_idx, column_idx] = value
    
    # Deallocate memory
    free(reverse_kmer)
    free(hashes)
    free(seeds)


@cython.cdivision(True)
def seq_bloomfilter_cython(str seq, dict config, cnp.npy_uint8[:, :] shared_array, int column_idx):
    """Computes hash for every k-mer in ever sequence in the input file(s) and creates the indexed bloom filter

    Args:
    seq_obj: Pyfastx's sequence object
    config: Config dictionary
    shared_array: SharedArray (Boolean NumPy array)
    column_idx: Column index in the shared array to use as bloom filter for the given sample

    Returns: None

    """
    # Declaration and initialization
    cdef uint64_t kmer_size = config["k"]
    cdef uint64_t num_hash = config["h"]
    cdef Py_ssize_t bf_size = config["m"]
    cdef str sequence_type = config["sequence_type"]
    cdef Py_ssize_t i, row_idx, seq_len
    cdef int h, s
    cdef const char *seq_cstr
    cdef char *seq_pointer
    cdef char *kmer
    cdef char *canonical
    cdef uint64_t hash_value
    cdef uint8_t bit_idx, value
    cdef char *reverse_kmer = <char *>malloc(kmer_size)
    cdef uint64_t *hashes = <uint64_t *>malloc(num_hash * sizeof(uint64_t))
    cdef uint64_t *seeds = <uint64_t *>malloc(num_hash * sizeof(uint64_t))

    for s in range(0, num_hash):
        seeds[s] = s

    # Bloom filter
    if sequence_type == "nucleotide":
        seq_cstr = PyUnicode_AsUTF8AndSize(seq, &seq_len)
        seq_pointer = &seq_cstr[0]
        for i in range(seq_len - kmer_size + 1):
            kmer = seq_pointer + i
            canonical = canonical_kmer(kmer, kmer_size, reverse_kmer)
            mmh2_hash_64(canonical, kmer_size, seeds, num_hash, hashes)
            for h in range(0, num_hash):
                hash_value = hashes[h] % bf_size
                bit_idx = hash_value % 8
                row_idx = hash_value // 8
                value = setbit(shared_array[row_idx, column_idx], bit_idx)
                shared_array[row_idx, column_idx] = value
    elif sequence_type == "aminoacid":
        seq_cstr = PyUnicode_AsUTF8AndSize(seq, &seq_len)
        seq_pointer = &seq_cstr[0]
        for i in range(seq_len - kmer_size + 1):
            kmer = seq_pointer + i
            mmh2_hash_64(kmer, kmer_size, seeds, num_hash, hashes)
            for h in range(0, num_hash):
                hash_value = hashes[h] % bf_size
                bit_idx = hash_value % 8
                row_idx = hash_value // 8
                value = setbit(shared_array[row_idx, column_idx], bit_idx)
                shared_array[row_idx, column_idx] = value
    
    # Deallocate memory
    free(reverse_kmer)
    free(hashes)
    free(seeds)


@cython.cdivision(True)
def sequence_to_hash_cython(str seq, str seq_type, dict config):
    """Computes hash for every k-mer in the given input sequence and return the hash value for every k-mer

    Args:
      sequence: Amino acid or nucleotide sequence
      seq_type: aminoacid or nucleotide

    Returns: List of hash values

    """
    # Declaration and initialization
    cdef uint64_t kmer_size = config["k"]
    cdef uint64_t num_hash = config["h"]
    cdef Py_ssize_t bf_size = config["m"]
    cdef Py_ssize_t i, seq_len
    cdef int h, s
    cdef const char *seq_cstr
    cdef char *seq_pointer
    cdef char *kmer
    cdef char *canonical
    cdef uint64_t hash_value    
    cdef list kmer_array = []
    cdef list hash_array = []
    cdef char *reverse_kmer = <char *>malloc(kmer_size)
    cdef uint64_t *hashes = <uint64_t *>malloc(num_hash * sizeof(uint64_t))
    cdef uint64_t *seeds = <uint64_t *>malloc(num_hash * sizeof(uint64_t))

    for s in range(0, num_hash):
        seeds[s] = s

    if seq_type == "nucleotide":
        seq_cstr = PyUnicode_AsUTF8AndSize(seq, &seq_len)
        seq_pointer = &seq_cstr[0]
        for i in range(seq_len - kmer_size + 1):
            kmer = seq_pointer + i
            canonical = canonical_kmer(kmer, kmer_size, reverse_kmer)
            kmer_array = []
            mmh2_hash_64(canonical, kmer_size, seeds, num_hash, hashes)
            for h in range(0, num_hash):
                hash_value = hashes[h] % bf_size
                kmer_array.append(hash_value)
            hash_array.append(kmer_array)

        # Deallocate memory
        free(reverse_kmer)
        free(hashes)
        free(seeds)
        return hash_array
    elif seq_type == "aminoacid":
        seq_cstr = PyUnicode_AsUTF8AndSize(seq, &seq_len)
        seq_pointer = &seq_cstr[0]
        for i in range(seq_len - kmer_size + 1):
            kmer = seq_pointer + i
            kmer_array = []
            mmh2_hash_64(kmer, kmer_size, seeds, num_hash, hashes)
            for h in range(0, num_hash):
                hash_value = hashes[h] % bf_size
                kmer_array.append(hash_value)
            hash_array.append(kmer_array)

        # Deallocate memory
        free(reverse_kmer)
        free(hashes)
        free(seeds)
        return hash_array


cdef char *canonical_kmer(char *kmer, int kmer_size, char *rev_comp):
    """To get canonical (lexicographically small) k-mer
    Args:
      kmer: K-mer
      kmer_size: Size of the k-mer
      rev_comp: memory to hold the reverse complement of k-mer
    Returns: Canonical k-mer
    NOTE: 
      1) Creates reverse complement to the given k-mer
      2) Returns the canonical k-mer
    """
    cdef char* basemap = [b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0',  b'T', b'\0', b'G', b'\0', b'\0', b'\0', b'C', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'N', b'\0', b'\0', b'\0', b'\0', b'\0', b'A', b'A', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b't', b'\0', b'g', b'\0', b'\0', b'\0', b'c', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'a',  b'a']
    cdef bint use_rev_comp = False
    cdef size_t end = kmer_size - 1
    cdef Py_ssize_t i

    # start filling reverse complement and stop once we know what is smaller
    for i in range(kmer_size):
        rev_comp[i] = basemap[<int>kmer[end - i]]
        if kmer[i] < rev_comp[i]:
            break
        elif rev_comp[i] < kmer[i]:
            use_rev_comp = True
            break

    if not use_rev_comp:
        return kmer

    # fill remainder of rev_comp
    for i in range(i + 1, kmer_size):
        rev_comp[i] = basemap[<int>kmer[end - i]]
    return rev_comp


cdef uint32_t mmh3_hash(char *data, uint32_t kmer_size, uint32_t *seeds, Py_ssize_t num_seeds, uint32_t *out_hashes):
    """Cython implementation of MurmurHash3_x86_32 unsigned
    Adapted from the python version https://github.com/wc-duck/pymmh3/blob/master/pymmh3.py

    Args:
      data: C byte string
      kmer_size: Size of the k-mer
      seeds: Memory pointer to the uint32_t seed(s)
      num_seeds: Number of seed(s) (equal to the number of hash functions to be applied on k-mer)
      out_hashes: Memory pointer to the uint32_t hash value(s) (to store the resulting hash value(s))

    Returns: Unsigned int 32 bit hash value(s)

    """
    cdef Py_ssize_t nblocks, block_start, tail_size
    cdef const uint32_t *data32
    cdef uint32_t length, c1, c2, k1, h, tail_index, x2, x3, x4, x5, result
    cdef signed long int x1
    data32 = <const uint32_t *> data

    length = kmer_size
    nblocks = int(length / 4)

    c1 = 0xCC9E2D51
    c2 = 0x1B873593
    x1 = 0xFFFFFFFF
    x2 = 0xE6546B64
    x3 = 0x85EBCA6B
    x4 = 0xC2B2AE35
    x5 = 0x80000000

    for s in range(0, num_seeds):
        out_hashes[s] = seeds[s]

    # body
    for block_start in range(0, nblocks * 4, 4):
        k1 = data32[block_start]
        k1 = (c1 * k1) & x1
        k1 = (k1 << 15 | k1 >> 17) & x1  # inlined ROTL32
        k1 = (c2 * k1) & x1

        for s in range(0, num_seeds):
            out_hashes[s] ^= k1
            out_hashes[s] = (out_hashes[s] << 13 | out_hashes[s] >> 19) & x1  # inlined ROTL32
            out_hashes[s] = (out_hashes[s] * 5 + x2) & x1

    # tail
    tail_index = nblocks * 4
    k1 = 0
    tail_size = length & 3

    if tail_size >= 3:
        k1 ^= data[tail_index + 2] << 16
    if tail_size >= 2:
        k1 ^= data[tail_index + 1] << 8
    if tail_size >= 1:
        k1 ^= data[tail_index + 0]

    if tail_size > 0:
        k1 = (k1 * c1) & x1
        k1 = (k1 << 15 | k1 >> 17) & x1  # inlined ROTL32
        k1 = (k1 * c2) & x1

        for s in range(0, num_seeds):
            out_hashes[s] ^= k1

    # finalization
    for s in range(0, num_seeds):
        h = out_hashes[s] ^ length
        h ^= h >> 16
        h = (h * x3) & x1
        h ^= h >> 13
        h = (h * x4) & x1
        h ^= h >> 16
        if h & x5 == 0:
            out_hashes[s] = h
        else:
            out_hashes[s] = -((h ^ x1) + 1)


cdef uint64_t mmh3_hash_64(const char *data, const uint64_t kmer_size, const uint64_t *seeds, const Py_ssize_t num_seeds, uint64_t *out_hashes):
    """Cython implementation of MurmurHash3 unsigned 64 bit

    128 bit implementation was used as reference

    Args:
      data: C byte string
      kmer_size: Size of the k-mer
      seeds: Memory pointer to the uint64_t seed(s)
      num_seeds: Number of seed(s) (equal to the number of hash functions to be applied on k-mer)
      out_hashes: Memory pointer to the uint64_t hash values (to store the resulting hash values)
    """
    cdef Py_ssize_t length, nblocks, block_start, tail_size, tail_index
    cdef const uint64_t *data64
    cdef uint64_t c1, c2, x2, x4, x5, x6, k1
    cdef int s
    data64 = <const uint64_t *> data

    for s in range(0, num_seeds):
        out_hashes[s] = seeds[s]

    length = kmer_size
    nblocks = int(length / 8)

    c1 = 0x87c37b91114253d5
    c2 = 0x4cf5ad432745937f
    x2 = 0x52dce729
    x4 = 0xff51afd7ed558ccd
    x5 = 0xc4ceb9fe1a85ec53

    # body
    for block_start in range(0, nblocks):
        k1 = data64[block_start]
        k1 = c1 * k1
        k1 = k1 << 31 | k1 >> 33 # inlined ROTL64
        k1 = c2 * k1
        for s in range(0, num_seeds):
            out_hashes[s] ^= k1
            out_hashes[s] = out_hashes[s] << 27 | out_hashes[s] >> 37 # inlined ROTL64
            out_hashes[s] = out_hashes[s] * 5 + x2
    
    # tail
    tail_index = nblocks * 8
    k1 = 0
    tail_size = length & 7

    if tail_size == 7:
        k1 ^= <uint64_t> data[tail_index + 6] << 48
        k1 ^= <uint64_t> data[tail_index + 5] << 40
        k1 ^= <uint64_t> data[tail_index + 4] << 32
        k1 ^= data[tail_index + 3] << 24
        k1 ^= data[tail_index + 2] << 16
        k1 ^= data[tail_index + 1] << 8
        k1 ^= data[tail_index + 0]
    elif tail_size == 6:
        k1 ^= <uint64_t> data[tail_index + 5] << 40
        k1 ^= <uint64_t> data[tail_index + 4] << 32
        k1 ^= data[tail_index + 3] << 24
        k1 ^= data[tail_index + 2] << 16
        k1 ^= data[tail_index + 1] << 8
        k1 ^= data[tail_index + 0]
    elif tail_size == 5:
        k1 ^= <uint64_t> data[tail_index + 4] << 32
        k1 ^= data[tail_index + 3] << 24
        k1 ^= data[tail_index + 2] << 16
        k1 ^= data[tail_index + 1] << 8
        k1 ^= data[tail_index + 0]
    elif tail_size == 4:
        k1 ^= data[tail_index + 3] << 24
        k1 ^= data[tail_index + 2] << 16
        k1 ^= data[tail_index + 1] << 8
        k1 ^= data[tail_index + 0]
    elif tail_size == 3:
        k1 ^= data[tail_index + 2] << 16
        k1 ^= data[tail_index + 1] << 8
        k1 ^= data[tail_index + 0]
    elif tail_size == 2:
        k1 ^= data[tail_index + 1] << 8
        k1 ^= data[tail_index + 0]
    elif tail_size == 1:
        k1 ^= data[tail_index + 0]
    if tail_size > 0:
        k1 = k1 * c1
        k1 = k1 << 31 | k1 >> 33  # inlined ROTL64
        k1 = k1 * c2
        for s in range(0, num_seeds):
            out_hashes[s] ^= k1

    #finalization
    for s in range(0, num_seeds):
        k1 = out_hashes[s] ^ length
        k1 ^= k1 >> 33
        k1 = k1 * x4
        k1 ^= k1 >> 33
        k1 = k1 * x5
        k1 ^= k1 >> 33
        out_hashes[s] = k1


cdef void mmh2_hash_64(const char *data, const uint64_t kmer_size, const uint64_t *seeds, const Py_ssize_t num_seeds, uint64_t *out_hashes):
    """Cython implementation of MurmurHash2_x64_64

    Adapted from https://github.com/aappleby/smhasher/blob/master/src/MurmurHash2.cpp

    Args:
      data: C byte string
      kmer_size: Size of the k-mer
      seeds: Memory pointer to the uint64_t seed(s)
      num_seeds: Number of seed(s) (equal to the number of hash functions to be applied on k-mer)
      out_hashes: Memory pointer to the uint64_t hash values (to store the resulting hash values)

    """
    cdef Py_ssize_t length, nblocks, block_start, tail_size, tail_index, s
    cdef const uint64_t *data64 = <const uint64_t *> data
    cdef const uint32_t *tail32
    cdef const uint16_t *tail16
    cdef const uint8_t *tail8
    cdef uint64_t m = 0xc6a4a7935bd1e995
    cdef uint64_t k
    cdef uint8_t *k8
    cdef uint16_t *k16

    # compile time detection of endianness
    cdef uint16_t endian16 = 1
    cdef uint8_t *endian8 = <uint8_t *> &endian16
    cdef bint little_endian = endian8[0]

    length = kmer_size
    nblocks = length // 8

    for s in range(0, num_seeds):
        out_hashes[s] = seeds[s] ^ (m * length)

    # body
    for block_start in range(0, nblocks):
        k = data64[block_start]
        k *= m
        k ^= k >> 47
        k *= m

        for s in range(0, num_seeds):
            out_hashes[s] ^= k
            out_hashes[s] *= m

    # tail
    tail_size = length & 7

    if tail_size > 0:
        # Read the remaining bytes into the least
        # significant bytes of a 64 bit integer.
        # The order within these least significant
        # bytes does not matter, because other
        # operations in the hash function already
        # are endianness dependent.

        tail8 = <uint8_t *> &data[nblocks * 8]
        tail16 = <uint16_t *> tail8
        tail32 = <uint32_t *> tail8

        k8 = <uint8_t *> &k
        k16 = <uint16_t *> &k

        if tail_size == 7:
            if little_endian:
                k = tail32[0]
                k16[2] = tail16[2]
                k8[6] = tail8[6]
            else:
                k = tail32[0]
                k16[1] = tail16[2]
                k8[1] = tail8[6]
        elif tail_size == 6:
            if little_endian:
                k = tail32[0]
                k16[2] = tail16[2]
            else:
                k = tail32[0]
                k16[1] = tail16[2]
        elif tail_size == 5:
            if little_endian:
                k = tail32[0]
                k8[4] = tail8[4]
            else:
                k = tail32[0]
                k8[3] = tail8[3]
        elif tail_size == 4:
            k = tail32[0]
        elif tail_size == 3:
            if little_endian:
                k = tail16[0]
                k8[2] = tail8[2]
            else:
                k = tail16[0]
                k8[5] = tail8[2]
        elif tail_size == 2:
            k = tail16[0]
        elif tail_size == 1:
            k = tail8[0]

        for s in range(0, num_seeds):
            out_hashes[s] ^= k
            out_hashes[s] *= m

    #finalization
    for s in range(0, num_seeds):
        out_hashes[s] ^= out_hashes[s] >> 47
        out_hashes[s] *= m
        out_hashes[s] ^= out_hashes[s] >> 47
