"""Key-value store using LMDB as backend"""

import os
import shutil
import sys
import lmdb
import msgspec
import psutil


class LMDBKVStore:
    """LMDB backed key-value store

    Args:
      path: Location of LMDB directory to use as the root folder (will be created if it does not exist)
      map_size: Maximum size database may grow to, in bytes [2199023255552]
      readonly: Whether to open the LMDB as read only [False]
      lock: Whether to lock LMDB during concurrent writing or reading [True]
      **kwargs: Keyword arguments passed through to the 'lmdb.open' function

    Notes:
      1) Map size is set to 2 TiB assuming it to be used on 64-bit Linux OS, change it to <2GB for 32-bit
      2) msgpec.MessagePack is used for serializing keys and values (Numpy arrays and data types are not supported as keys and values)
      3) Flush method (flush()) must be called exclusively to make sure that the writes are synced to the disk
      4) Value should be of type bytes as we do not do any serialization of the values

    """

    def __init__(
        self, path, map_size=2199023255552, readonly=False, lock=True, **kwargs
    ):
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
        self.env = lmdb.open(self.path, **kwargs)

        # Serializer and deserializer initialization
        self.encoder = msgspec.msgpack.Encoder()
        self.decoder = msgspec.msgpack.Decoder()

    def __getitem__(self, key):
        """Get item from the DB

        Args:
          key: key

        Returns: Value if key is present else None

        """
        with self.env.begin() as txn:
            value = txn.get(self.serialize(key))

        if value:
            return value
        else:
            return None

    def __setitem__(self, key, value):
        """Set item in the DB

        Args:
          key: key
          value: value (Only bytes type is allowed as we do not do any serialization)

        """
        with self.env.begin(write=True) as txn:
            txn.put(self.serialize(key), value)

    def __delitem__(self, key):
        """Delete item from the DB

        Args:
          key: key

        """
        with self.env.begin(write=True) as txn:
            txn.delete(self.serialize(key))

    def __contains__(self, key):
        """Check if key is present in the DB

        Args:
          key: key

        Returns: True if key is present else False

        """
        with self.env.begin() as txn:
            return txn.get(self.serialize(key)) is not None

    def __iter__(self):
        """Iterate over the DB keys and values"""
        with self.env.begin() as txn:
            with txn.cursor() as cursor:
                for key, value in cursor.iternext(keys=True, values=True):
                    yield self.deserialize(key), value

    def __len__(self):
        """Total number of items present in the DB"""
        with self.env.begin() as txn:
            return txn.stat()["entries"]

    def __enter__(self):
        """Context manager entry point"""
        return self

    def __exit__(self, *args):
        """Context manager exit point"""
        self.close()

    def keys(self):
        """Get all keys from the DB"""
        keys_list = []
        with self.env.begin() as txn:
            with txn.cursor() as cursor:
                for key in cursor.iternext(keys=True, values=False):
                    keys_list.append(self.deserialize(key))

        return keys_list

    def values(self):
        """Get all values from the DB"""
        values_list = []
        with self.env.begin() as txn:
            with txn.cursor() as cursor:
                for value in cursor.iternext(keys=False, values=True):
                    values_list.append(value)

        return values_list

    def pop(self, key):
        """Pop an item from the DB

        Args:
          key: key

        Returns: Value if key is present else None

        """
        with self.env.begin(write=True) as txn:
            value = txn.get(self.serialize(key))
            if value:
                txn.delete(self.serialize(key))
                return value
            else:
                return None

    def get_multi(self, keys):
        """Get multiple items from the DB

        Args:
          keys: List of keys

        Returns: Dictionary of key-value pairs. If key is not present in the DB, it is not present in the returned dictionary

        """
        local_dict = {}

        with self.env.begin() as txn:
            with txn.cursor() as cursor:
                for key in keys:
                    value = cursor.get(self.serialize(key))
                    if value:
                        local_dict[key] = value

        return local_dict

    def set_multi(self, keys, values):
        """Set multiple items in the DB

        Values should be of type bytes as we do not do any serialization

        Args:
          keys: List of keys to insert in the DB. Duplicate keys will be overwritten
          values: List of values to insert in the DB

        """
        if len(keys) != len(values):
            raise ValueError("Keys and values should be of same length")

        with self.env.begin(write=True) as txn:
            with txn.cursor() as cursor:
                for key, value in zip(keys, values):
                    cursor.put(self.serialize(key), value)

    def serialize(self, item):
        """Serialization using msgpsec.msgpack

        Args:
          item: Any item (int, str, list, dict) to serialize

        Returns: Serialized object

        """
        return self.encoder.encode(item)

    def deserialize(self, item):
        """Deserialization using MessagePack

        Args:
          item: Any MessagePack serialized object

        Returns: Deserialized object

        """
        return self.decoder.decode(item)

    def flush(self):
        """Sync/Write data to the file system"""
        self.env.sync()

    def close(self):
        """Close the DB"""
        self.env.close()

    def cleanup(self):
        """Deletes the LMDB index folder"""
        # Check if directory exists and delete it
        if os.path.exists(self.path):
            shutil.rmtree(self.path)
