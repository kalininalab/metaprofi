# Human RNA-seq-mini dataset index construction

This directory contains all the information on how we constructed indexes for RNA-seq-mini dataset using several tools.

### Data download

* We randomly selected 650 sequence files from the previously downloaded data (from RNA-seq index construction).
* Their accession id's can be found [here](https://github.com/kalininalab/metaprofi/blob/main/benchmarks/rna_seq_mini_index/rna_seq_mini_sra_accessions.txt).

### Input preparation

* We first constructed compacted DBG's using BCALM2
* The compacted DBG's were then used as input for all the tools (except for Mantis)

```bash
# BCALM 2 installation
conda create --name bcalm bcalm -c bioconda

# Activate conda environment
conda activate bcalm

# Run BCALM2 (provided: 170 cores, 512 GiB RAM)
# The file _rna_seq_650_list.txt_ can be found in this folder
cat rna_seq_650_list.txt | while read line; do name=$(basename $line .fasta.gz); bcalm -in $line -kmer-size 21 -abundance-min 2 -nb-cores 170 -out ./benchmarks/bcalm_cdbgs/$name -out-dir ./benchmarks/bcalm_cdbgs -minimizer-type 0 -repartition-type 0 -max-memory 524288; done
```

### Results

#### Index construction
| Tool | RAM (GiB) | CPU cores | Disk BF (GiB) | Disk index (GiB) | Disk total (GiB) | Time BF (min) | Time index (min) | Time total (min) |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| HowDeSBT | 2.4 | 64 | 152 | 4.4 |	168.4 | 51.94 |	93.1 | 145.1 |
| kmtricks + HowDeSBT |	286.39 + 4 | 64	| 156 | 4.4 | 308 + 12 | 54.49 | 383.9 | 438.39 |
| MetaProFi | 12 | 64 |	8.2 | 9.4 |	17.6 | 4.37	| 20.42 | 24.79 |
| COBS | 12 | 64 | - | 51 | 54 | - | 8.59 | 8.59 |
| Squeakr + MANTIS | 7.3 + 49 |	64 | 28 | 14 | 42 | 145.22 | 33.27 | 178.49 |

_BF: Bloom filter, Disk total: total storage used for BF, index and intermediate files, Time total: total time for constructing BF and index, - = N/A_

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
| HowDeSBT | 0.61 | - | 22 | 558 |
| kmtricks + HowDeSBT |	0.64 | 20 |	2952 | 2957 |
| MetaProFi | 1.9 |	20 | 29 | 33 |
| COBS | 25.1 | 20 | 234 | 228 |
| Squeakr + MANTIS | 14 | - | 37 | - |

_- = N/A_


### Commands

#### HowDeSBT

```bash
# Installation
conda create --name howdesbt howdesbt -c bioconda

# Activate conda environment
conda activate howdesbt

# Index construction
## BF
for file in $(ls ./benchmarks/bcalm_cdbgs/unitigs_files/*); do
	base=$(basename $file .unitigs.fa)
	howdesbt makebf $file --k=21 --min=1 --bits=2000000000 --threads=64 --hashes=1 --out=./benchmarks/howdesbt/bloomfilters/${base}.bf
done

## Cluster
ls ./benchmarks/howdesbt/bloomfilters*.bf > ./benchmarks/howdesbt/bf_list.txt

howdesbt cluster --list=./benchmarks/howdesbt/bf_list.txt --bits=20000000 --tree=./benchmarks/howdesbt/uncompressed.culled.sbt --nodename=node{number}.bf

## Index
howdesbt build --determined,brief --rrr --tree=./benchmarks/howdesbt/uncompressed.culled.sbt --outtree=./benchmarks/howdesbt/howde.culled.rrr.sbt

# Querying
## Extract sequence search
howdesbt query --threshold=1.0 --tree=./benchmarks/howdesbt/howde.culled.rrr.sbt /data/rna_seq/query/query_1000.fa > ./benchmarks/howdesbt/rna_seq_mini_query_t100.txt

## Approximate sequence search
howdesbt query --threshold=0.75 --tree=./benchmarks/howdesbt/howde.culled.rrr.sbt /data/rna_seq/query/query_1000.fa > ./benchmarks/howdesbt/rna_seq_mini_query_t75.txt

```

#### kmtricks + HowDeSBT

```bash
# Installation
conda create --name kmtricks kmtricks -c tlemane

# Activate conda environment
conda activate kmtricks

# Index construction
## BF
kmtricks pipeline --run-dir ./benchmarks/kmtricks_rna_seq_mini --file .benchmarks/kmtricks_rna_seq_mini/kmtricks_input.txt --hard-min 1 --kmer-size 21 --threads 64 --mode hash:bft:bin --cpr --bloom-size 2000000000 --bf-format howdesbt

## Index
kmtricks index --run-dir ./benchmarks/kmtricks_rna_seq_mini/ --bits 20000000 --howde -t 64

# Querying
## Extract sequence search
kmtricks query --run-dir ./benchmarks/kmtricks_rna_seq_mini/ --query /data/rna_seq/query/query_1000.fa --threshold 1.0 --sort -t 20 > ./benchmarks/kmtricks_rna_seq_mini/rna_seq_mini_query_t100.txt

## Approximate sequence search
kmtricks query --run-dir ./benchmarks/kmtricks_rna_seq_mini/ --query /data/rna_seq/query/query_1000.fa --threshold 1.0 --sort -t 20 > ./benchmarks/kmtricks_rna_seq_mini/rna_seq_mini_query_t100.txt
```


#### MetaProFi
- [config](https://github.com/kalininalab/metaprofi/blob/master/benchmarks/rna_seq_index/config.yml) file was created

``` bash
# Activate conda environment
conda activate metaprofi

# Index construction
metaprofi build ./benchmarks/metaprofi/input.txt config.yml

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
ls ./benchmarks/bcalm_cdbgs/unitigs_files/* > ./benchmarks/cobs_mini/cobs_rna_seq_mini_input.list

# Index construction
./cobs/build/src/cobs compact-construct --file-type list -T 64 -k 21 -m 68720000000 ./benchmarks/cobs_mini/cobs_rna_seq_mini_input.list ./benchmarks/cobs_mini/rna_seq_mini.cobs_compact

# Querying
## Exact sequence search
./cobs/build/src/cobs query -i ./benchmarks/cobs_mini/rna_seq_mini.cobs_compact -T 20 -f /data/rna_seq/query/query_1000.fa -t 1.0 > ./benchmarks/cobs_mini/rna_seq_mini_query_t100.txt

## Approximate sequence search
./cobs/build/src/cobs query -i ./benchmarks/cobs_mini/rna_seq_mini.cobs_compact -T 20 -f /data/rna_seq/query/query_1000.fa -t 0.75 > ./benchmarks/cobs_mini/rna_seq_mini_query_t75.txt
```

#### Squeakr + MANTIS

```bash
# Installation
conda create --name mantis -c bioconda squeakr mantis

# Activate conda environment
conda activate mantis

# Index construction
## Run squeakr
for file in $(ls ./benchmarks/datasets/rna_seq_mini_files/*fastq.gz); do
	base=$(basename $file .fastq.gz)
	squeakr count -e -k 21 -c 2 -n -s 32 -t 64 -o ./benchmarks/mantis/squeaker_files/${base}.squeakr ${file}
done

## Index
ls ./benchmarks/mantis/squeaker_files/* > ./benchmarks/mantis/mantis_input.txt

mantis build -s 32 -i ./benchmarks/mantis/mantis_input.txt -o ./benchmarks/mantis/mantis_index

# Querying
## Exact sequence search
mantis query -1 -p ./benchmarks/mantis/mantis_index -o ./benchmarks/mantis/rna_seq_mini_query_t100.res /data/rna_seq/query/query_1000.fa
```