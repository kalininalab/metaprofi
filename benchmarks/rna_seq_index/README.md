# Human RNA-seq dataset index construction

This directory contains all the information on how we constructed MetaProFi index for RNA-seq dataset.

### Data download

``` bash
# Activate conda environment
conda activate metaprofi

# Install dependency
conda install parallel-fastq-dump
```

- Run data download custom python [script](https://github.com/kalininalab/metaprofi/blob/master/benchmarks/rna_seq_index/rna_seq_data_download.py)
- Size of the compressed data is 2.7 TiB

### MetaProFi commands

- [config](https://github.com/kalininalab/metaprofi/blob/master/benchmarks/rna_seq_index/config.yml) and [input](https://github.com/kalininalab/metaprofi/blob/master/benchmarks/rna_seq_index/rna_seq_input.txt) files were created

``` bash
# Run MetaProFi
metaprofi build input_data.txt config.yml
```

### Results

#### Bloom filter matrix construction

| Time (min) | RAM (GiB) | CPU cores | Disk (GiB) |
| --- | --- | --- | --- |
| 1399.15 | 59 | 64 | 295 |

#### Index construction

| Time (min) | RAM (GiB) | CPU cores | Disk (GiB) |
| --- | --- | --- | --- |
| 240.168 | 56 | 64 | 333 |

#### Querying

- Update `nproc` to 20 in config.yml

``` bash
# Activate conda environment
conda activate metaprofi

# Data download (https://github.com/kamimrcht/REINDEER/blob/master/reproduce_manuscript_results/queries.sh)
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/refMrna.fa.gz

# Randomly extract 1000 sequences
pyfastx sample /data/rna_seq/query/refMrna.fa.gz -n 1000 -o /data/rna_seq/query/query_1000.fa
```

- Query reads can be found [here](https://github.com/kalininalab/metaprofi/blob/master/benchmarks/rna_seq_index/query_1000.fa)

``` bash
# Exact sequence search
metaprofi search_index config.yml -f /data/rna_seq/query/query_1000.fa -i nucleotide

# Approximate sequence search
metaprofi search_index config.yml -f /data/rna_seq/query/query_1000.fa -i nucleotide -t 75
```

| Search type | Time (seconds) | RAM (GiB) | CPU cores | Results |
| --- | --- | --- | --- | --- |
| Exact search (T=100%) | 349 | 1.9 | 20 | [results](https://github.com/kalininalab/metaprofi/blob/master/benchmarks/rna_seq_index/metaprofi_query_results-05_08_2021-06_56_31_t100.txt) |
| Approximate search (T=75%) | 359 | 1.9 | 20 | [results](https://github.com/kalininalab/metaprofi/blob/master/benchmarks/rna_seq_index/metaprofi_query_results-05_08_2021-07_17_31_t75.txt) |