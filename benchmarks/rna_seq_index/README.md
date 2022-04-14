# Human RNA-seq dataset index construction

This directory contains all the information on how we constructed MetaProFi and COBS index for RNA-seq dataset.

### Data download

``` bash
# Activate conda environment
conda activate metaprofi

# Install dependency
conda install parallel-fastq-dump
```

- Run data download custom python [script](https://github.com/kalininalab/metaprofi/blob/master/benchmarks/rna_seq_index/rna_seq_data_download.py)
- Size of the compressed data is 2.7 TiB

### Results

#### Index construction
| Tool | RAM (GiB) | CPU cores | Disk BF (GiB) | Disk index (GiB) | Disk total (GiB) | Time BF (min) | Time index (min) | Time total (min) |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| MetaProFi | 59 | 64 | 295 | 333 | 628 | 1108 | 127 | 1235 |
| COBS | 69.4 | 64 | - | 935 | 996 | - | 1000 | 1000 |

_- = N/A_

#### Querying

* Data download (https://github.com/kamimrcht/REINDEER/blob/master/reproduce_manuscript_results/queries.sh)


    ``` bash
    wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/refMrna.fa.gz
    ```

* Randomly extract 1000 sequences

    ``` bash
    pyfastx sample /data/rna_seq/query/refMrna.fa.gz -n 1000 -o /data/rna_seq/query/query_1000.fa
    ```
* Query reads can be found [here](https://github.com/kalininalab/metaprofi/blob/master/benchmarks/rna_seq_index/query_1000.fa)

| Tool | RAM (GiB) | CPU cores | Time (s) (T = 100) | Time (s) (T = 75) |
| --- | --- | --- | --- | --- |
| MetaProFi | 3.4 | 64 | 43 | 48 |
| COBS | 92.5 | 64 | 290 | 290 |

### Commands

#### MetaProFi
- [config](https://github.com/kalininalab/metaprofi/blob/master/benchmarks/rna_seq_index/config.yml) and [input](https://github.com/kalininalab/metaprofi/blob/master/benchmarks/rna_seq_index/rna_seq_input.txt) files were created

``` bash
# Activate conda environment
conda activate metaprofi

# Index construction
metaprofi build input_data.txt config.yml

# Querying
# Exact sequence search
metaprofi search_index config.yml -f /data/rna_seq/query/query_1000.fa -i nucleotide

# Approximate sequence search
metaprofi search_index config.yml -f /data/rna_seq/query/query_1000.fa -i nucleotide -t 75
```

#### COBS

``` bash
# Installation
## Create conda environment
conda create --name cobs -c conda-forge gxx compilers cmake libboost

## Activate conda environment
conda activate cobs

## Install COBS
git clone --recursive https://github.com/bingmann/cobs.git
mkdir cobs/build
cd cobs/build
cmake ..
make -j8

# Input file creation
ls /data/rna_seq/*.fasta.gz > ./benchmarks/cobs/cobs_rna_seq_input.list

# Index construction
./cobs/build/src/cobs compact-construct --file-type list -T 64 -k 21 -m 68720000000 ./benchmarks/cobs/cobs_rna_seq_input.list ./benchmarks/cobs/rna_seq.cobs_compact

# Querying
## Exact sequence search
./cobs/build/src/cobs query -i ./benchmarks/cobs/rna_seq.cobs_compact -T 64 -f /data/rna_seq/query/query_1000.fa -t 1.0 > ./benchmarks/cobs/rna_seq_query_t100.txt

## Approximate sequence search
./cobs/build/src/cobs query -i ./benchmarks/cobs/rna_seq.cobs_compact -T 64 -f /data/rna_seq/query/query_1000.fa -t 0.75 > ./benchmarks/cobs/rna_seq_query_t75.txt
```
