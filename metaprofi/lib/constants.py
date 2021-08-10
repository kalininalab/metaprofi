"""MetaProFi constant variables module
"""

# Constants for Zarr functions
BLOOMFILTER_DATASET_NAME = "BloomFilters"
METADATA_DATASET_NAME = "Metadata"
BLOOMFILTER_INDEX_DATASET_NAME = "MetaProFiIndex"
COMPRESSION_LEVEL = 3
COMPRESSION_ALGO = "zstd"

# File names
MATRIX_STORE = "metaprofi_bfmatrix_store.zarr"
INDEX_STORE = "metaprofi_index_store.zarr"
LMDB_INDEX_STORE = "metaprofi_index_store.lmdb"
