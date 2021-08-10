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

#### Bloom filter matrix construction

| Time (min) | RAM (GiB) | CPU cores | Disk (GiB) |
| --- | --- | --- | --- |
| 2998.34 | <58 | 64 | 640 |

#### Index construction

| Time (min) | RAM (GiB) | CPU cores | Disk (GiB) |
| --- | --- | --- | --- |
| 252.24 | 70 | 64 | 751 |

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
metaprofi search_index config.yml -f query_reads.fastq -i nucleotide -t 40
```

| Search type | Time (seconds) | RAM (GiB) | CPU cores | Results |
| --- | --- | --- | --- | --- |
| Exact search (T=100%) | 12.877 | 1.3 | 64 | [results](https://github.com/kalininalab/metaprofi/blob/master/benchmarks/tara_oceans_index/metaprofi_query_results-16_07_2021-11_54_08_t100.txt) |
| Approximate search (T=40%) | 13.048 | 1.3 | 64 | [results](https://github.com/kalininalab/metaprofi/blob/master/benchmarks/tara_oceans_index/metaprofi_query_results-16_07_2021-11_54_40_t40.txt) |

### kmtricks

#### Installation

``` bash
# 1: Create conda environment
conda create --name kmtricks python==3.7.7 kmtricks -c tlemane

# 2: Activate conda environment
conda activate kmtricks

# 3: Version check
kmtricks.py --version

# Output
kmtricks v0.0.6, git_sha1 : 8539f16
```

#### Bloom filter construction as outlined [here](https://github.com/pierrepeterlongo/kmtricks_benchmarks/tree/master/tara-metag-bacterial/expe_kmtricks) and [here](https://github.com/tlemane/kmtricks)

- kmtricks input file can be found [here](https://github.com/kalininalab/metaprofi/blob/master/benchmarks/tara_oceans_index/kmtricks_input.txt)
- Two sets of parameters were used
    1.  Most of the options were left to the default values including `--max-memory`

    ``` bash
    kmtricks.py run --file kmtricks_input.txt --kmer-size 31 --run-dir /data/kmtricks/tara_oceans/ --count-abundance-min 1 --mode bf_trp --nb-cores 64 --lz4 --log-files merge --max-hash 40000000000 --hasher sabuhash --split howde
    ```

    2. Using parameters reported in [kmtricks](https://github.com/pierrepeterlongo/kmtricks_benchmarks/tree/master/tara-metag-bacterial/expe_kmtricks)

    ``` bash
    kmtricks.py run --file kmtricks_input.txt --kmer-size 31 --run-dir /data/kmtricks/tara_oceans/ --count-abundance-min 1 --max-count 256 --max-memory 8000 --mode bf_trp --nb-cores 64 --lz4 --merge-abundance-min 3 --recurrence-min 1 --save-if 1 --log-files merge --max-hash 40000000000 --hasher sabuhash --split howde
    ```

- Three run attempts were made to construct Bloom filters
    1. Run1: Using parameter set 1 (defaults)
        1. Got terminated after 9 hrs due to the consumption of the entire disk space of 2.9 TiB
        2. Raised the following error

        ``` bash
        OSError: [Errno 28] No space left on device
        ```

    2. Run2: Using parameter set 1 (defaults)
        1. Just changed the output directory to a large non-RAID NVME NFS file system
        2. Observed a peak usage of 4.8 TiB of storage before manually terminating the run after 100 hrs
    3. Run3: Using parameter set 2 (same as applied in [kmtricks paper](https://github.com/pierrepeterlongo/kmtricks_benchmarks/tree/master/tara-metag-bacterial/expe_kmtricks))
        1. Got terminated after 11 hrs
        2. Used the large non-RAID NVME NFS file system for storing the output
        3. Raised the following error

        ``` bash
        Repartition: 1/1, Superkmer: 143/249, Count: 62749/655119, Merge: 0/2631, Output: 0/1

        Signal SIGSEGV received from km_reads_to_superk with the following arguments:
        {'id': 174, 'f': '/data/tara_oceans/ERR599054/ERR599054_1.fastq.gz,/data/tara_oceans/ERR599054/ERR599054_2.fastq.gz', 'fof': <__main__.Fof obje
        ct at 0x7f3aaf162e90>, 'log': None, 'exp_id': 'ERR599054', 'verbose': False, 'debug': False, 'cmd': 'run', 'file': 'input.txt', 'run_dir': '/data/kmtricks/tara_oceans/', 'kmer_size': 31, 'count_abundance_min': 1, 'abundance_max': 3000000000, 'max_count': 256, 'max_memory': 8000, 'mode': 'bf_trp', 'kff_output': 0, 'log_files': 'merg
        e', 'nb_cores': 1, 'merge_abundance_min': '3', 'recurrence_min': 1, 'save_if': 1, 'skip_merge': 0, 'until': 'all', 'only': 'all', 'minimizer_type': 0, 'minimizer_size': 10, 'repartition_typ
        e': 0, 'nb_partitions': 2631, 'hasher': 'sabuhash', 'max_hash': 40000000000, 'split': 'howde', 'keep_tmp': 0, 'lz4': 1, 'abundance_min': 1}.
        All children are killed. Check your inputs. If the problem persists, please contact us with a description of your run and the following files: ./km_backtrace/backtrace.log and ~/anaconda3/envs/kmtricks/bin/build/build_infos.txt.
        ```

        4. Observed a peak usage of 2.1 TiB of storage

#### Results

- Using the second run's partial results

| Time (min) | RAM (GiB) | CPU cores | Disk (TiB) |
| --- | --- | --- | --- |
| - | <50 | 64| 4.8 |