# Tara Oceans dataset index construction

This directory contains all the information on how we constructed MetaProFi index for Tara Oceans dataset.

### Data download

- Download the Tara Oceans dataset [PRJEB1787](https://www.ebi.ac.uk/ena/browser/view/PRJEB1787) (Run accession ids are provided in [accessions.txt](https://github.com/kalininalab/metaprofi/blob/master/benchmarks/tara_oceans_index/accessions.txt) file)
- Contains 249 samples and a total of 495 FASTQ files
- Size of the compressed data is 4 TiB

### MetaProFi commands

- [config](https://github.com/kalininalab/metaprofi/blob/master/benchmarks/tara_oceans_index/config.yml) and [input](https://github.com/kalininalab/metaprofi/blob/master/benchmarks/tara_oceans_index/input_data.txt) files were created
    - MetaProFi input file was created using the [custom script](https://github.com/kalininalab/metaprofi/blob/master/benchmarks/tara_oceans_index/tara_oceans_input_prep.py)

``` bash
# Activate conda environment
conda activate metaprofi

# Run MetaProFi
metaprofi build input_data.txt config.yml
```

### Results

#### Index construction

| Tool | RAM (GiB) | CPU cores | Disk BF (GiB) | Disk index (GiB) | Disk total (GiB) | Time BF (min) | Time index (min) | Time total (min) |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| MetaProFi | 68 | 64 | 643 | 750 | 1393 | 2642 | 279 | 2921 |
| kmtricks + HowDeSBT | 47 | 64 | 1228.8 | >390* | 2344 + 390* | 865 | >7217.42* | >8082.42* |

_BF: Bloom filter, Disk total is the total storage used for BF, index and intermediate files, Time total: total time used to construct BF and index, *: terminated after 120 hrs, data reported as it is at the time of termination._

#### Querying

- 1000 reads were extracted randomly from 5 different accession ids (200 reads per accession id)

- Run accessions list: ERR598955, ERR599059, ERR598964, ERR599064, ERR599165

    ``` bash
    zcat /data/tara_oceans/ERR598955/ERR598955_1.fastq.gz | head -n 800 > query_reads.fastq
    zcat /data/tara_oceans/ERR599059/ERR599059_1.fastq.gz | head -n 800 >> query_reads.fastq
    zcat /data/tara_oceans/ERR598964/ERR598964_1.fastq.gz | head -n 800 >> query_reads.fastq
    zcat /data/tara_oceans/ERR599064/ERR599064_1.fastq.gz | head -n 800 >> query_reads.fastq
    zcat /data/tara_oceans/ERR599165/ERR599165_1.fastq.gz | head -n 800 >> query_reads.fastq
    ```

- Query reads can be found [here](https://github.com/kalininalab/metaprofi/blob/master/benchmarks/tara_oceans_index/query_reads.fastq)

``` bash
# Exact sequence search
metaprofi search_index config.yml -f query_reads.fastq -i nucleotide

# Approximate sequence search
metaprofi search_index config.yml -f query_reads.fastq -i nucleotide -t 75
```

* MetaProFi query results

    | Search type | Time (seconds) | RAM (GiB) | CPU cores |
    | --- | --- | --- | --- |
    | Exact search (T=100%) | 164 | 14.3 | 64 |
    | Approximate search (T=75%) | 166 | 14.3 | 64 |

* kmtricks + HowDeSBT query results: We were unable to compare the query results with the kmtricks + HowDeSBT setup as the index construction had to be terminated after 120 hrs.

### kmtricks

#### Installation

``` bash
# 1: Create conda environment
conda create --name kmtricks kmtricks -c tlemane

# 2: Activate conda environment
conda activate kmtricks

# 3: Version check
kmtricks --version

# Output
kmtricks v1.1.1
```

#### Bloom filter and index construction as outlined [here](https://github.com/pierrepeterlongo/kmtricks_benchmarks/tree/master/tara-metag-bacterial/expe_kmtricks) and [here](https://github.com/tlemane/kmtricks/wiki)

- kmtricks input file can be found [here](https://github.com/kalininalab/metaprofi/blob/master/benchmarks/tara_oceans_index/kmtricks_input.txt)

``` bash
# Activate conda environment
conda activate kmtricks

# Bloom filter construction
kmtricks pipeline --run-dir benchmarks/kmtricks_tara_oceans --file kmtricks_input.txt --hard-min 1 --kmer-size 31 --threads 64 --mode hash:bft:bin --cpr --bloom-size 40000000000 --bf-format howdesbt

# Index construction
kmtricks index --run-dir benchmarks/kmtricks_tara_oceans --bits 400000000 --howde -t 64
```
