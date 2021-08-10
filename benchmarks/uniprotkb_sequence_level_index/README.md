# UniProtKB sequence-level bacterial dataset index construction

This directory contains all the information on how we constructed **sequence-level** MetaProFi index for the UniProtKB bacterial dataset.

### Data download

* UniPrtoKB bacterial dataset downloaded for the organism-level indexing was reused for the sequence-level indexing

### Dataset preparation

* A single compressed FASTA file containing both Swiss-Prot and TrEMBL sequences were constructed using the same [custom script](https://github.com/kalininalab/metaprofi/blob/master/benchmarks/uniprotkb_organism_level_index/uniprotkb_parser.py) used in the data preparation of the organism-level indexing

### MetaProFi commands

- Config [file](https://github.com/kalininalab/metaprofi/blob/master/benchmarks/uniprotkb_sequence_level_index/config.yml)

``` bash
# Activate conda environment
conda activate metaprofi

# Run MetaProFi
metaprofi build-seq /data/uniprotkb_datasets/uniprot_bacteria.fasta.gz config.yml
```

### Results

#### Bloom filter matrix construction

|Time (min) | RAM (GiB) | CPU cores | Disk (GiB) |
| --- | --- | --- | --- |
| 66.33 | <60 | 64 | 232 |


#### Index construction

| Time (min) | RAM (GiB) | CPU cores | Disk (GiB) |
| --- | --- | --- | --- |
| 862.57 | <60 | 64 | 210 |