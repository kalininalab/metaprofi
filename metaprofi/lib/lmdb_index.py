"""LMDB backed key-value store for FASTA/FASTQ file index
"""

import gzip
import itertools
import os
import subprocess
import sys
from multiprocessing import Process
import indexed_gzip as igzip
import lmdb
import msgpack
import psutil
import zstd
from metaprofi.lib.utilities import cleanup_dir, is_gzip_file, check_seq


class LMDBStore:
    """LMDB backed key-value store

    Args:
      path: Location of LMDB directory to use as the root folder to store fasta and gzip index
      map_size: Maximum size database may grow to, in bytes [2199023255552]
      min_seq_len: Minimum length of the sequence to index (for example: k-mer size)
      readonly: Whether to open the LMDB as read only [False]
      lock: Whether to lock LMDB during concurrent writing [True]
      **kwargs: Keyword arguments passed through to the 'lmdb.open' function

    Notes:
      1) Map size is set to 2 TiB assuming it to be used on 64-bit Linux OS, change it to <2GB for 32-bit
      2) MessagePack is used for serializing keys and values
      3) Zstandard is used for compressing the serialized values
      4) Flush and close methods must be called exclusively to make sure that the writes are synced to the disk

    Index:
      FASTA:
        1) Creates a Gzip index for the Gzip compressed input file for fast seeking, otherwise uses the input file directly for seeking
        2) For every sequence larger than the min_seq_len we create a seek index which contains the sequence id, name of the sequence, sequence seek location, byte length of the sequence and the number of bases in the sequence
      FASTQ:
        1) Creates a Gzip index for the Gzip compressed input file for fast seeking, otherwise uses the input file directly for seeking
        2) For every read  larger than the min_seq_len we create a seek index which contains the sequence id, name of the read, read seek location, byte length of the read and the number of bases in the read

    """

    def __init__(
        self,
        path,
        min_seq_len,
        map_size=2199023255552,
        readonly=False,
        lock=True,
        **kwargs,
    ):

        # Check if command line tool pigz is installed
        self.check_program()

        # Set LMDB defaults and user provided config to 'lmdb.open'
        # Default map size is set to 2TiB (only on 64 bit systems, should be < 2GB for 32 bit systems)
        kwargs.setdefault("map_size", map_size)
        kwargs.setdefault("readonly", readonly)
        kwargs.setdefault("metasync", False)
        kwargs.setdefault("sync", False)

        # Enable writemap only on Linux systems
        # LMDB doc: This option may cause filesystems that don’t support
        # sparse files, such as OSX, to immediately preallocate map_size= bytes
        # of underlying storage when the environment is opened or closed for
        # the first time.
        if sys.platform.startswith("linux"):
            kwargs.setdefault("writemap", True)

        kwargs.setdefault("map_async", False)
        kwargs.setdefault("max_spare_txns", psutil.cpu_count())
        kwargs.setdefault("lock", lock)

        self.readonly = readonly
        self.path = os.path.abspath(path)
        self.min_seq_len = min_seq_len
        self.env = lmdb.open(self.path, **kwargs)
        self.file_type = None
        self.compressed = None
        self.file = None

        if self.readonly:
            with self.env.begin() as txn:
                self.compressed = self.deserialize(
                    txn.get(self.serialize("compressed"))
                )
                self.input_file = self.deserialize(
                    txn.get(self.serialize("input_file"))
                )
                if self.compressed:
                    self.gzip_index_file = self.deserialize(
                        txn.get(self.serialize("gzip_index_file"))
                    )
                    self.file = igzip.IndexedGzipFile(
                        self.input_file, index_file=self.gzip_index_file
                    )
                else:
                    self.file = open(self.input_file, "rb")
                self.type = check_seq(self.__getitem__(0)[1])

    def __getitem__(self, key):
        """Get item from the DB

        Args:
          key: key

        Returns: Value if key is present else an empty list is returned as a default value

        """
        with self.env.begin() as txn:
            value = txn.get(self.serialize(key))

        if not value:
            return []
        result = self.decompress(value)
        name, seq_len = result[1], result[4]
        self.file.seek(result[2])
        seq = self.file.read(result[3]).decode().replace("\n", "")
        return name, seq.upper(), seq_len

    def __len__(self):
        """Total number of sequences present in the DB"""
        # """Total number of index items present in the DB"""
        # Subtract 2 keys as they correspond to the Gzip index file names
        # return self.env.stat()["entries"] - 2
        with self.env.begin() as txn:
            length = self.deserialize(txn.get(self.serialize("SeqCount")))
        return int(length)

    def serialize(self, item):
        """Serialization using MessagePack

        Args:
          item: Any item (int, str, list, dict) to serialize

        Returns: Serialized object

        """
        return msgpack.dumps(item)

    def deserialize(self, item):
        """Deserialization using MessagePack

        Args:
          item: Any MessagePack serialized object

        Returns: Deserialized object

        """
        return msgpack.loads(item)

    def compress(self, value):
        """Serialize and Zstandard compress the value

        Args:
          value: Any value, could be list, dict, string, int

        Returns: Compressed serialized bytes

        """
        return zstd.compress(self.serialize(value), 9, 2)

    def decompress(self, compressed_value):
        """Zstandard decompress and deserialize the compressed value

        Args:
          compressed_value: Any bytes that was compressed using 'compress' function

        Returns: Decompressed and deserialized value

        """
        return self.deserialize(zstd.decompress(compressed_value))

    def get_seq_names(self):
        """To get the names of the sequences

        Returns: Generator object with the sequence names

        """
        with self.env.begin() as txn:
            with txn.cursor() as cursor:
                for i in range(0, self.__len__()):
                    yield self.decompress(cursor.get(self.serialize(i)))[1]

    def build_index(self, input_file):
        """Builds index for a input file (FASTA/FASTQ (.GZ))

        Args:
        input_file: FASTA/FASTQ(.GZ) file

        """
        self.input_file = os.path.abspath(input_file)
        self.file_type, self.compressed = self.check_file_type(self.input_file)
        processes = []

        if self.file_type == "FASTA":
            lmdb_index_process = Process(target=self.__build_fasta_index)
            lmdb_index_process.start()
            processes.append(lmdb_index_process)
        elif self.file_type == "FASTQ":
            lmdb_index_process = Process(target=self.__build_fastq_index)
            lmdb_index_process.start()
            processes.append(lmdb_index_process)

        if self.compressed:
            gzip_index = Process(target=self.__build_gzip_index)
            gzip_index.start()
            processes.append(gzip_index)
        else:
            with self.env.begin(write=True) as txn:
                with txn.cursor() as cursor:
                    cursor.put(self.serialize("input_file"), self.serialize(input_file))
                    cursor.put(self.serialize("compressed"), self.serialize(0))

        for proc in processes:
            proc.join()

    def __build_fasta_index(self):
        """Create a LMDB index of the sequence file (FASTA)"""
        input_file = self.input_file
        file_stdout = self.get_file_handle(input_file, self.compressed)
        results = []
        seq_idx = 1
        previous_offset = 0
        header, seq = None, []
        key = 0

        with self.env.begin(write=True) as txn:
            with txn.cursor() as cursor:
                for line in file_stdout.stdout:
                    line = line.decode()
                    if line.startswith(">"):
                        if header:
                            seq_bytes = len("".join(seq))
                            seq_len = len("".join(s.strip() for s in seq))
                            previous_offset = current_line_offset + seq_bytes
                            if seq_len >= self.min_seq_len:
                                results[0].extend([seq_bytes, seq_len])
                                cursor.put(
                                    self.serialize(key), self.compress(results[0])
                                )
                                key += 1
                            results = []
                        header, seq = line, []
                        current_line_offset = previous_offset + len(header)
                        results.append(
                            [seq_idx, header.split(" ")[0][1:], current_line_offset]
                        )
                        seq_idx += 1
                    else:
                        seq.append(line)
                if header:
                    seq_bytes = len("".join(seq))
                    seq_len = len("".join(s.strip() for s in seq))
                    previous_offset = current_line_offset + seq_bytes
                    if seq_len >= self.min_seq_len:
                        results[0].extend([seq_bytes, seq_len])
                        cursor.put(self.serialize(key), self.compress(results[0]))
                        key += 1
                cursor.put(self.serialize("SeqCount"), self.serialize(key))

    def __build_fastq_index(self):
        """Create a LMDB index of the sequence file (FASTQ)"""
        input_file = self.input_file
        file_stdout = self.get_file_handle(input_file, self.compressed)
        seq_idx = 1
        previous_read_offset = 0
        current_read_offset = 0
        key = 0

        with self.env.begin(write=True) as txn:
            with txn.cursor() as cursor:
                for header, seq, sep, qual in itertools.zip_longest(
                    *[file_stdout.stdout] * 4
                ):
                    header = header.decode()
                    seq = seq.decode()
                    seq_bytes = len(seq)
                    seq_len = seq_bytes - seq.count("\n")
                    current_read_offset = previous_read_offset + len(header)
                    if seq_len >= self.min_seq_len:
                        cursor.put(
                            self.serialize(key),
                            self.compress(
                                [
                                    seq_idx,
                                    header.split(" ")[0][1:],
                                    current_read_offset,
                                    seq_bytes,
                                    seq_len,
                                ]
                            ),
                        )
                        key += 1
                    seq_idx += 1
                    previous_read_offset = (
                        current_read_offset + len(seq) + len(sep) + len(qual)
                    )
                cursor.put(self.serialize("SeqCount"), self.serialize(key))

    def __build_gzip_index(self):
        """Builds Gzip index for the Gzip compressed input file using indexed_gzip library"""
        input_file = self.input_file
        gzip_index_file = (
            f"{self.path}/{os.path.splitext(os.path.basename(input_file))[0]}.gzidx"
        )
        with self.env.begin(write=True) as txn:
            with txn.cursor() as cursor:
                cursor.put(self.serialize("input_file"), self.serialize(input_file))
                cursor.put(
                    self.serialize("gzip_index_file"),
                    self.serialize(gzip_index_file),
                )
                cursor.put(self.serialize("compressed"), self.serialize(1))

        gzip_obj = igzip.IndexedGzipFile(input_file)
        gzip_obj.build_full_index()
        gzip_obj.export_index(gzip_index_file)

    def check_file_type(self, input_file):
        """Checks and returns the input sequence file type

        Args:
        input_file: User provided sequence input file

        Returns: File type (FASTA or FASTQ) and a boolean for compressed file or not

        """
        if not os.path.exists(input_file):
            self.close()
            self.cleanup()
            raise FileNotFoundError("Input file does not exist")
        if not os.path.isfile(input_file):
            self.close()
            self.cleanup()
            raise IsADirectoryError("Input should be a file (FASTA/FASTQ (.GZ))")
        if not os.stat(input_file).st_size > 0:
            self.close()
            self.cleanup()
            raise ValueError("Input file is empty")

        compressed = is_gzip_file(input_file)
        line = ""

        if compressed:
            with gzip.open(input_file) as in_f:
                line = in_f.readline()
        else:
            with open(input_file, "rb") as in_f:
                line = in_f.readline()

        if line.startswith(b">"):
            return "FASTA", compressed
        if line.startswith(b"@"):
            return "FASTQ", compressed

        # If the input file is not a FASTA/FASTQ format
        self.close()
        cleanup_dir(self.path)
        raise ValueError("Input file should be a FASTA/FASTQ (.GZ) file")

    @staticmethod
    def get_file_handle(input_file, compressed):
        """Opens the input file using systems command line tools (pigz/cat) and returns the stdout

        Args:
        input_file: User provided sequence input file
        compressed: Boolean

        Returns: Subprocess handle

        """
        if compressed:
            file_out = subprocess.Popen(
                ["pigz", "-cd", input_file], stdout=subprocess.PIPE, bufsize=-1
            )
        else:
            file_out = subprocess.Popen(
                ["cat", input_file], stdout=subprocess.PIPE, bufsize=-1
            )
        return file_out

    @staticmethod
    def check_program():
        """Check if command line tool (pigz) is installed"""
        status, _ = subprocess.getstatusoutput("which pigz")
        if status != 0:
            raise ValueError(
                "The command line tool 'pigz' is not available. Please install!"
            )

    def flush(self):
        """Sync/Write data to the file system"""
        self.env.sync()

    def close(self):
        """Close the DB and opened file handle"""
        self.env.close()

        if self.readonly:
            self.file.close()

    def cleanup(self):
        """Deletes the LMDB index folder"""
        cleanup_dir(self.path)
